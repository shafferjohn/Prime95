/*----------------------------------------------------------------------
| gwnum.c
|
| This file contains the C routines and global variables that are used
| in the multi-precision arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
| 
|  Copyright 2002-2019 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <time.h>
#include "cpuid.h"
#include "gwnum.h"
#include "gwtables.h"
#include "gwini.h"
#include "gwutil.h"
#include "gwdbldbl.h"
#include "gwbench.h"

//#define GDEBUG_MEM	1			// Print out memory used

/* Option to share sin/cos data where possible amongst several gwnum callers */
/* There is little downside to enabling this option (more memory allocated) and */
/* little upside (very slight decrease in working set size and needed memory bandwidth). */

#ifndef GDEBUG_MEM
#define SHARE_SINCOS_DATA
#endif

/* Include a random number generator.  For reasons we'll discuss later */
/* we do not want to use the C runtime library's random number generator */
/* to initialize GW_RANDOM */

#include "mt19937ar.c"

/* Handy macros to improve readability */

#define log2(n)				(log((double)(n)) / log (2.0))
#define logb(n)				(log((double)(n)) / log ((double)(b)))
#define divide_rounding_up(a,b)		((a + (b) - 1) / (b))
#define divide_rounding_down(a,b)	((a) / (b))
#define round_up_to_multiple_of(a,b)	(divide_rounding_up (a, b) * (b))
#define round_down_to_multiple_of(a,b)	(divide_rounding_down (a, b) * (b))
#define fltmod(a,b)			((double)(a) - floor ((double)(a) / (double)(b)) * (double)(b))
#define round_to_cache_line(p)		(void *) (((intptr_t)(p) + 63) & ~63)

/* MSVC6 has trouble with the pow function using integer arguments. */
/* For example, "(unsigned long) pow (5.0, 7.0)" returns 78124 instead */
/* of the correct 78125.  This macro, works around this trouble. */

#define intpow(b,n)	((long) floor (pow ((double)(b), (double)(n)) + 0.1))

/* global variables */

/* When debugging gwnum and giants, I sometimes write code that "cheats" */
/* by calling a routine that is part of prime95 rather than the gwnum and */
/* giants library.  Prime95 will set this routine pointer so that gwnum */
/* code can cheat while keeping the gwnum library interface clean. */

void (*OutputBothRoutine)(int, const char *) = NULL;
#define OutputBoth(t,x)	(*OutputBothRoutine)(t,x)

/* Assembly helper routines */

struct gwinfo1_data {
	const struct gwasm_jmptab *sse2_cyclic_fft_info; /* SSE2 mersenne mod FFT info */
	const struct gwasm_jmptab *sse2_negacyclic_fft_info; /* SSE2 2^N+1 mod FFT */
	const struct gwasm_jmptab *x86_cyclic_fft_info; /* x86 mersenne mod FFT */
	const struct gwasm_jmptab *x86_negacyclic_fft_info; /* x86 2^N+1 mod FFT */
	const struct gwasm_jmptab *avx_cyclic_fft_info; /* AVX mersenne mod FFT */
	const struct gwasm_jmptab *avx_negacyclic_fft_info; /* AVX 2^N+1 mod FFT */
	const struct gwasm_jmptab *avx512_cyclic_fft_info; /* AVX mersenne mod FFT */
	const struct gwasm_jmptab *avx512_negacyclic_fft_info; /* AVX 2^N+1 mod FFT */
	uint32_t version;			/* Gwnum lib version number */
};

void gwinfo1 (struct gwinfo1_data *);
void prefetchL2 (void *addr, int count);
void pause_for_count (int count);

/* gwnum assembly routine pointers */

#define gw_fft(h,a)	(*(h)->GWPROCPTRS[0])(a)
#define gw_add(h,a)	(*(h)->GWPROCPTRS[1])(a)
#define gw_addq(h,a)	(*(h)->GWPROCPTRS[2])(a)
#define gw_addf(h,a)	(*(h)->GWPROCPTRS[2])(a)
#define gw_sub(h,a)	(*(h)->GWPROCPTRS[3])(a)
#define gw_subq(h,a)	(*(h)->GWPROCPTRS[4])(a)
#define gw_subf(h,a)	(*(h)->GWPROCPTRS[4])(a)
#define gw_addsub(h,a)	(*(h)->GWPROCPTRS[5])(a)
#define gw_addsubq(h,a)	(*(h)->GWPROCPTRS[6])(a)
#define gw_addsubf(h,a)	(*(h)->GWPROCPTRS[6])(a)
#define gw_copyzero(h,a) (*(h)->GWPROCPTRS[7])(a)
#define gw_adds(h,a)	(*(h)->GWPROCPTRS[8])(a)
#define gw_muls(h,a)	(*(h)->GWPROCPTRS[9])(a)
#define norm_routines	10
#define zerohigh_routines 14

#define addr_gw_fft(h)		((h)->GWPROCPTRS[0])
#define addr_gw_add(h)		((h)->GWPROCPTRS[1])
#define addr_gw_addq(h)		((h)->GWPROCPTRS[2])
#define addr_gw_addf(h)		((h)->GWPROCPTRS[2])
#define addr_gw_sub(h)		((h)->GWPROCPTRS[3])
#define addr_gw_subq(h)		((h)->GWPROCPTRS[4])
#define addr_gw_subf(h)		((h)->GWPROCPTRS[4])
#define addr_gw_addsub(h)	((h)->GWPROCPTRS[5])
#define addr_gw_addsubq(h)	((h)->GWPROCPTRS[6])
#define addr_gw_addsubf(h)	((h)->GWPROCPTRS[6])
#define addr_gw_copyzero(h)	((h)->GWPROCPTRS[7])
#define addr_gw_adds(h)		((h)->GWPROCPTRS[8])
#define addr_gw_muls(h)		((h)->GWPROCPTRS[9])
#define call_op(p,a)		(*p)(a)

/* Macro to aid in build prctabs (a prctab is an array of pointers to assembly routines) */

#define extern_decl(name)		extern void (*name)(void);
#define array_entry(name)		&name,

/* Build an array of the assembly auxiliary routines (addition, subtraction, etc.) */

#define aux_decl(name)			extern_decl(name##1)	extern_decl(name##2)
#define aux_decl3(name)			extern_decl(name##1)	extern_decl(name##2)	extern_decl(name##3)
#define aux_decl3y(name)		aux_decl3(name)	aux_decl3(name##r) aux_decl3(name##zp) aux_decl3(name##rzp) aux_decl3(name##n) aux_decl3(name##nr) aux_decl3(name##nzp) aux_decl3(name##nrzp)
#define aux_decl3z(name)		aux_decl3(name)	aux_decl3(name##r) aux_decl3(name##zp) aux_decl3(name##rzp)

#ifndef X86_64
aux_decl(gwadd)		aux_decl(gwaddq)
aux_decl(gwsub)		aux_decl(gwsubq)
aux_decl(gwaddsub)	aux_decl(gwaddsubq)
aux_decl(gwcopyzero)	aux_decl(gwmuls)

void *x87_aux_prctab[] = {
	&gwadd1, &gwaddq1, &gwsub1, &gwsubq1, &gwaddsub1, &gwaddsubq1, &gwcopyzero1, NULL, &gwmuls1,
	&gwadd2, &gwaddq2, &gwsub2, &gwsubq2, &gwaddsub2, &gwaddsubq2, &gwcopyzero2, NULL, &gwmuls2};
#endif

aux_decl3(gwxadd)	aux_decl(gwxaddq)
aux_decl3(gwxsub)	aux_decl(gwxsubq)
aux_decl3(gwxaddsub)	aux_decl(gwxaddsubq)
aux_decl(gwxcopyzero)	aux_decl3(gwxadds)	aux_decl3(gwxmuls)

void *sse2_aux_prctab[] = {
	&gwxadd1, &gwxaddq1, &gwxsub1, &gwxsubq1, &gwxaddsub1, &gwxaddsubq1, &gwxcopyzero1, &gwxadds1, &gwxmuls1,
	&gwxadd2, &gwxaddq2, &gwxsub2, &gwxsubq2, &gwxaddsub2, &gwxaddsubq2, &gwxcopyzero2, &gwxadds2, &gwxmuls2,
	&gwxadd3, &gwxaddq2, &gwxsub3, &gwxsubq2, &gwxaddsub3, &gwxaddsubq2, &gwxcopyzero2, &gwxadds3, &gwxmuls3};

aux_decl3y(gwyadd)	aux_decl3(gwyaddq)
aux_decl3y(gwysub)	aux_decl3(gwysubq)
aux_decl3y(gwyaddsub)	aux_decl3(gwyaddsubq)
aux_decl3(gwycopyzero)	aux_decl3y(gwyadds)	aux_decl3y(gwymuls)

void *avx_aux_prctab[] = {
	&gwyaddq1, &gwysubq1, &gwyaddsubq1, &gwycopyzero1,
	&gwyadd1, &gwysub1, &gwyaddsub1, &gwyadds1, &gwymuls1,
	&gwyaddr1, &gwysubr1, &gwyaddsubr1, &gwyadds1, &gwymulsr1,
	&gwyaddzp1, &gwysubzp1, &gwyaddsubzp1, &gwyadds1, &gwymulszp1,
	&gwyaddrzp1, &gwysubrzp1, &gwyaddsubrzp1, &gwyadds1, &gwymulsrzp1,
	&gwyaddn1, &gwysubn1, &gwyaddsubn1, &gwyaddsn1, &gwymulsn1,
	&gwyaddnr1, &gwysubnr1, &gwyaddsubnr1, &gwyaddsn1, &gwymulsnr1,
	&gwyaddnzp1, &gwysubnzp1, &gwyaddsubnzp1, &gwyaddsn1, &gwymulsnzp1,
	&gwyaddnrzp1, &gwysubnrzp1, &gwyaddsubnrzp1, &gwyaddsn1, &gwymulsnrzp1,

	&gwyaddq3, &gwysubq3, &gwyaddsubq3, &gwycopyzero3,
	&gwyadd3, &gwysub3, &gwyaddsub3, &gwyadds3, &gwymuls3,
	&gwyaddr3, &gwysubr3, &gwyaddsubr3, &gwyaddsr3, &gwymulsr3,
	&gwyaddzp3, &gwysubzp3, &gwyaddsubzp3, &gwyadds3, &gwymulszp3,
	&gwyaddrzp3, &gwysubrzp3, &gwyaddsubrzp3, &gwyaddsr3, &gwymulsrzp3,
	&gwyaddn3, &gwysubn3, &gwyaddsubn3, &gwyaddsn3, &gwymulsn3,
	&gwyaddnr3, &gwysubnr3, &gwyaddsubnr3, &gwyaddsnr3, &gwymulsnr3,
	&gwyaddnzp3, &gwysubnzp3, &gwyaddsubnzp3, &gwyaddsn3, &gwymulsnzp3,
	&gwyaddnrzp3, &gwysubnrzp3, &gwyaddsubnrzp3, &gwyaddsnr3, &gwymulsnrzp3};

aux_decl3y(ygw_carries_wpn)

void *avx_carries_prctab[] = {
	&ygw_carries_wpn3, &ygw_carries_wpnr3, &ygw_carries_wpnzp3, &ygw_carries_wpnrzp3,
	&ygw_carries_wpnn3, &ygw_carries_wpnnr3, &ygw_carries_wpnnzp3, &ygw_carries_wpnnrzp3};

#ifdef X86_64
aux_decl3z(gwzadd)	aux_decl3(gwzaddq)
aux_decl3z(gwzsub)	aux_decl3(gwzsubq)
aux_decl3z(gwzaddsub)	aux_decl3(gwzaddsubq)
aux_decl3(gwzcopyzero)	aux_decl3y(gwzadds)	aux_decl3y(gwzmuls)

void *avx512_aux_prctab[] = {
	&gwzaddq3, &gwzsubq3, &gwzaddsubq3, &gwzcopyzero3,
	&gwzadd3, &gwzsub3, &gwzaddsub3, &gwzadds3, &gwzmuls3,
	&gwzaddr3, &gwzsubr3, &gwzaddsubr3, &gwzaddsr3, &gwzmulsr3,
	&gwzaddzp3, &gwzsubzp3, &gwzaddsubzp3, &gwzadds3, &gwzmulszp3,
	&gwzaddrzp3, &gwzsubrzp3, &gwzaddsubrzp3, &gwzaddsr3, &gwzmulsrzp3};

aux_decl3z(zgw_carries_wpn)
extern_decl(zgw_carries_op_wpnzp3)
extern_decl(zgw_carries_op_wpnrzp3)
void gwz3_apply_carries (void *);

void *avx512_carries_prctab[] = {
	&zgw_carries_wpn3, &zgw_carries_wpnr3, &zgw_carries_wpnzp3, &zgw_carries_wpnrzp3, &zgw_carries_op_wpnzp3, &zgw_carries_op_wpnrzp3};
#else
#define gwz3_apply_carries(a)
#endif

/* Now we put the normalization routines in an array so we can easily pick */
/* the normalization routines to use.  There is one table for AVX-512, AVX, SSE2 and */
/* x87.  The SSE2 table has these 820 combinations: */
/*	r or i		(rational or irrational) */
/*	1 or 2 or 2AMD or 3 or 3AMD (1 pass or 2 pass or 2 pass with partial normalization - with optional AMD prefetching) */
/*	z or zp or blank (zero upper half of result or zero-padded FFT or normal FFT */
/*	e or blank	(roundoff error checking or not) */
/*	b or blank	(b > 2 or not, not used in AVX-512) */
/*	s4 or blank	(SSE4 or not) */
/*	k or blank	(k for XMM_K_HI is zero or not) */
/*	c1 or cm1 or blank (c=1, c=-1, abs(c)!=1, not used in AVX-512) */
/* We also define a macro that will pick the correct entry from the array. */

/* First is the AVX-512 table */

#ifdef X86_64
#define avx512_explode(macro)			avx512_explode1(macro,zr)		avx512_explode1(macro,zi)
#define avx512_explode1(macro,name)		avx512_explode2(macro,name##3,SKX)
#define avx512_explode2(macro,name,suff)	avx512_zero_explode(macro,name##z,suff)	avx512_explode3(macro,name,suff,notzp)	avx512_explode3(macro,name##zp,suff,zp)
#define avx512_zero_explode(macro,name,suff)	avx512_explode9##suff(macro,name)	avx512_explode9##suff(macro,name##e)
#define avx512_explode3(macro,name,suff,zp)	avx512_explode4(macro,name,suff,zp)	avx512_explode4(macro,name##e,suff,zp)
#define avx512_explode4(macro,name,suff,zp)	avx512_explode5(macro,name,suff,zp,notc) avx512_explode5(macro,name##c,suff,zp,c)
#define avx512_explode5(macro,name,suff,zp,c)	avx512_explode6(macro,name,suff,zp,c)
#define avx512_explode6(macro,name,suff,zp,c)	avx512_explode7##zp(macro,name,suff,c)
#define avx512_explode7notzp(macro,name,suff,c)	avx512_explode9##suff(macro,name)
#define avx512_explode7zp(macro,name,suff,c)	avx512_explode8##c(macro,name,suff)	avx512_explode8##c(macro,name##k,suff)
#define avx512_explode8c(macro,name,suff)	avx512_explode9##suff(macro,name)
#define avx512_explode8notc(macro,name,suff)	avx512_explode9##suff(macro,name)
#define avx512_explode9SKX(macro,name)		macro(name##SKX)

avx512_explode(extern_decl)
void *avx512_prctab[] = { avx512_explode(array_entry) NULL };

int avx512_prctab_index (gwhandle *gwdata, int z, int e, int c)
{
	int	index = 0;

	if (! gwdata->RATIONAL_FFT) index += 14;		/* Irrational FFTs are after the rational FFTs */

	if (z) {
		if (e) index += 1;
		return (index);
	}
	index += 2;
	if (! gwdata->ZERO_PADDED_FFT) {
		if (e) index += 2;
		if (c) index += 1;
	} else {
		struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
		index += 4;
		if (e) index += 4;
		if (!c) {
			if (asm_data->u.zmm.ZMM_K_HI_OVER_SMALL_BASE == 0.0) index += 1;
		} else {
			index += 2;
			if (asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_SMALL_BASE == 0.0) index += 1;
		}
	}
	return (index);
}
#endif

/* Now the AVX table */

#define avx_explode(macro)			avx_explode1(macro,yr)			avx_explode1(macro,yi)
#define avx_explode1(macro,name)		avx_explode2(macro,name##1,CORE)	avx_explode2(macro,name##1,FMA3)	avx_explode2(macro,name##3,CORE)	avx_explode2(macro,name##3,FMA3)
#define avx_explode2(macro,name,suff)		avx_zero_explode(macro,name##z,suff)	avx_explode3(macro,name,suff,notzp)	avx_explode3(macro,name##zp,suff,zp)
#define avx_zero_explode(macro,name,suff)	avx_explode9##suff(macro,name)		avx_explode9##suff(macro,name##e)
#define avx_explode3(macro,name,suff,zp)	avx_explode4(macro,name,suff,zp)	avx_explode4(macro,name##e,suff,zp)
#define avx_explode4(macro,name,suff,zp)	avx_explode5(macro,name,suff,zp,notc)	avx_explode5(macro,name##c,suff,zp,c)
#define avx_explode5(macro,name,suff,zp,c)	avx_explode6(macro,name,suff,zp,c)	avx_explode6(macro,name##b,suff,zp,c)
#define avx_explode6(macro,name,suff,zp,c)	avx_explode7##zp(macro,name,suff,c)
#define avx_explode7notzp(macro,name,suff,c)	avx_explode9##suff(macro,name)
#define avx_explode7zp(macro,name,suff,c)	avx_explode8##c(macro,name,suff)	avx_explode8##c(macro,name##k,suff)
#define avx_explode8c(macro,name,suff)		avx_explode9##suff(macro,name)
#define avx_explode8notc(macro,name,suff)	avx_explode9##suff(macro,name)		avx_explode9##suff(macro,name##c1)	avx_explode9##suff(macro,name##cm1)
#define avx_explode9BLEND(macro,name)		macro(name##BLEND)
#define avx_explode9CORE(macro,name)		macro(name##CORE)
#ifdef X86_64
#define avx_explode9FMA3(macro,name)		macro(name##FMA3)
#else
#define avx_explode9FMA3(macro,name)		macro(name##CORE)			/* We don't support FMA FFTs on 32-bit builds */
#endif

avx_explode(extern_decl)
void *avx_prctab[] = { avx_explode(array_entry) NULL };

int avx_prctab_index (gwhandle *gwdata, int z, int e, int c)
{
	int	index = 0;

	if (! gwdata->RATIONAL_FFT) index += 168;		/* Irrational FFTs are after the 4 sets of rational FFTs */
	if (gwdata->PASS2_SIZE) index += 84;			/* Two-pass FFTs are after the one-pass FFTs */
	if (gwdata->cpu_flags & CPU_FMA3) index += 42;		/* FMA3 FFTs are after the CORE FFTs */

	if (z) {
		if (e) index += 1;
		return (index);
	}
	index += 2;
	if (! gwdata->ZERO_PADDED_FFT) {
		if (e) index += 4;
		if (c) index += 2;
		if (gwdata->b > 2) index += 1;
	} else {
		struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
		index += 8;
		if (e) index += 16;
		if (!c) {
			if (gwdata->b > 2) index += 6;
			if (asm_data->u.ymm.YMM_K_HI[0] == 0.0) index += 3;
			if (gwdata->c == 1) index += 1;  if (gwdata->c == -1) index += 2;
		} else {
			index += 12;
			if (gwdata->b > 2) index += 2;
			if (asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[0] == 0.0) index += 1;
		}
	}
	return (index);
}				    

/* Now the SSE2 table */

#define sse2_explode(macro)			sse2_explode1(macro,xr)			sse2_explode1(macro,xi)
#define sse2_explode1(macro,name)		sse2_explode2(macro,name##1,BLEND)	sse2_explode2(macro,name##2,CORE)	sse2_explode2(macro,name##2,K8)		sse2_explode2(macro,name##3,CORE)	sse2_explode2(macro,name##3,K8)
#define sse2_explode2(macro,name,suff)		sse2_zero_explode(macro,name##z,suff)	sse2_explode3(macro,name,suff,notzp)	sse2_explode3(macro,name##zp,suff,zp)
#define sse2_zero_explode(macro,name,suff)	sse2_explode9##suff(macro,name)		sse2_explode9##suff(macro,name##e)
#define sse2_explode3(macro,name,suff,zp)	sse2_explode4(macro,name,suff,zp)	sse2_explode4(macro,name##e,suff,zp)
#define sse2_explode4(macro,name,suff,zp)	sse2_explode5(macro,name,suff,zp,notc)	sse2_explode5(macro,name##c,suff,zp,c)
#define sse2_explode5(macro,name,suff,zp,c)	sse2_explode6(macro,name,suff,zp,c)	sse2_explode6(macro,name##b,suff,zp,c)
#define sse2_explode6(macro,name,suff,zp,c)	sse2_explode7##zp(macro,name,suff,c)	sse2_explode7##zp(macro,name##s4,suff,c)
#define sse2_explode7notzp(macro,name,suff,c)	sse2_explode9##suff(macro,name)
#define sse2_explode7zp(macro,name,suff,c)	sse2_explode8##c(macro,name,suff)	sse2_explode8##c(macro,name##k,suff)
#define sse2_explode8c(macro,name,suff)		sse2_explode9##suff(macro,name)
#define sse2_explode8notc(macro,name,suff)	sse2_explode9##suff(macro,name)		sse2_explode9##suff(macro,name##c1)	sse2_explode9##suff(macro,name##cm1)
#define sse2_explode9BLEND(macro,name)		macro(name##BLEND)
#define sse2_explode9CORE(macro,name)		macro(name##CORE)
#ifndef __APPLE__
#define sse2_explode9K8(macro,name)		macro(name##K8)
#else
#define sse2_explode9K8(macro,name)		macro(name##CORE)		/* Macs do not have AMD CPUs */
#endif

sse2_explode(extern_decl)
void *sse2_prctab[] = { sse2_explode(array_entry) NULL };

int sse2_prctab_index (gwhandle *gwdata, int z, int e, int c)
{
	int	index = 0;

	if (! gwdata->RATIONAL_FFT) index += 410;
	if (gwdata->PASS2_SIZE) {
		index += 82;
		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) index += 164;
		if (gwdata->cpu_flags & CPU_3DNOW_PREFETCH) index += 82;
	}
	if (z) {
		if (e) index += 1;
		return (index);
	}
	index += 2;
	if (! gwdata->ZERO_PADDED_FFT) {
		if (e) index += 8;
		if (c) index += 4;
		if (gwdata->b > 2) index += 2;
		if (gwdata->cpu_flags & CPU_SSE41) index += 1;
	} else {
		struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
		index += 16;
		if (e) index += 32;
		if (!c) {
			if (gwdata->b > 2) index += 12;
			if (gwdata->cpu_flags & CPU_SSE41) index += 6;
			if (asm_data->u.xmm.XMM_K_HI[0] == 0.0) index += 3;
			if (gwdata->c == 1) index += 1;  if (gwdata->c == -1) index += 2;
		} else {
			index += 24;
			if (gwdata->b > 2) index += 4;
			if (gwdata->cpu_flags & CPU_SSE41) index += 2;
			if (asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[0] == 0.0) index += 1;
		}
	}
	return (index);
}				    

/* The x87 normalization routines array has 40 combinations: */
/*	r or i		(rational or irrational) */
/*	1 or 2		(1 or 2 pass FFTs) */
/*	z or zp or blank (zero upper half of result or zero-padded FFT or normal FFT */
/*	e or blank	(roundoff error checking or not) */
/*	c or blank	(mul by small const or not)  */
/* We also define a macro that will pick the correct entry from the array. */

#ifndef X86_64
#define x87_explode(macro)		x87_explode1(macro,r)		x87_explode1(macro,i)
#define x87_explode1(macro,name)	x87_explode2(macro,name##1)	x87_explode2(macro,name##2)
#define x87_explode2(macro,name)	x87_zero_explode(macro,name##z)	x87_explode3(macro,name)	x87_explode3(macro,name##zp)
#define x87_zero_explode(macro,name)	macro(name)			macro(name##e)
#define x87_explode3(macro,name)	x87_explode4(macro,name)	x87_explode4(macro,name##e)
#define x87_explode4(macro,name)	macro(name)			macro(name##c)

x87_explode(extern_decl)
void *x87_prctab[] = { x87_explode(array_entry) NULL };

#define	x87_prctab_index(gwdata, z, e, c)  \
	    ((gwdata)->RATIONAL_FFT ? 0 : 20) + \
	    ((gwdata)->PASS2_SIZE ? 10 : 0) + \
	    (z ? (e ? 1 : 0) : 2 + \
	     ((gwdata)->ZERO_PADDED_FFT ? 4 : 0) + \
	     (e ? 2 : 0) + \
	     (c ? 1 : 0))
#endif

/* Helper macros */

#ifndef isinf
#define isinf(x)		((x) != 0.0 && (x) == 2.0*(x))
#endif
#ifndef isnan
#define isnan(x)		((x) != (x))
#endif
#define is_valid_double(x)	(! isnan (x) && ! isinf (x))

/* More #defines */

#define	GWINIT_WAS_CALLED_VALUE	0xAA12BB34	/* Special value checked by gwsetup to ensure gwinit was called */
#define GWFREED_TEMPORARILY	0x80000000	/* Special flag value set in the gwnum freeable field */

/* Forward declarations */

int convert_giant_to_k2ncd (
	giant	g,		/* Giant to examine */
	double	*k,		/* K in (K*2^N+C)/D. */
	unsigned long *n,	/* N in (K*2^N+C)/D. */
	signed long *c,		/* C in (K*2^N+C)/D. */
	unsigned long *d);	/* D in (K*2^N+C)/D. */
int internal_gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c);		/* C in K*B^N+C. Must be rel. prime to K. */
long nonbase2_gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	a,
	gwnum	g,
	unsigned long limit,	/* How many FFT words to set */
	unsigned long offset,	/* Offset into FFT array of words to set */
	long	carry);		/* Carry to add into this section */
void raw_gwsetaddin (gwhandle *gwdata, unsigned long word, double val);
int multithread_init (gwhandle *gwdata);
void multithread_term (gwhandle *gwdata);
void do_multithread_op_work (gwhandle *gwdata, struct gwasm_data *asm_data);
void pass1_aux_entry_point (void*);
void pass2_aux_entry_point (void*);
void create_auxiliary_hyperthread (struct gwasm_data *);
void auxiliary_hyperthread (void *);

/* Routine to split a r4dwpn FFT word into column and group multiplier indexes */
/* We remove bit(s) associated with the upper SSE2/AVX/AVX-512 words because those are */
/* handled by the group multipliers. */

unsigned long dwpn_col (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long word)
{
	unsigned long high, low;

	ASSERTG (! (gwdata->cpu_flags & CPU_AVX512F));	// Except for ADDIN_VALUE, AVX-512 FFT data is normalized by the time most C code sees it

	if (gwdata->cpu_flags & CPU_AVX) {
		low = word % gwdata->PASS2_SIZE;
		high = (((word / gwdata->PASS2_SIZE) & ~3) % gwdata->wpn_count) * gwdata->PASS2_SIZE;
	} else {
		unsigned long upper_sse2_word = gwdata->PASS2_SIZE >> 1;
		high = word / upper_sse2_word;
		low = word - high * upper_sse2_word;

		high = high >> 1;
		high = high % gwdata->wpn_count;
		high = high * gwdata->PASS2_SIZE;
	}
	return (high + low);
}

/* Return true if half of the words would have the same pattern of big */
/* and little words.  */

int match_pathological_pattern (
	unsigned long num_big_words,
	unsigned long total_length,
	double	pathological_fraction)
{
	double	half_length;
	unsigned long pathological_base, actual_base;

/* Compute the gwfft_base you would get of the word half way into the FFT if */
/* you had the pathological fraction of big and little words */

	half_length = (double) total_length * 0.5;
	pathological_base = (unsigned long) ceil (half_length * pathological_fraction);

/* Compute the base you would get given the actual fraction of big words */

	actual_base = (unsigned long) ceil (half_length * (double) num_big_words / (double) total_length);

/* Return pathological (true) if the actual_base is close to the pathological_base */

	return (actual_base >= pathological_base && actual_base <= pathological_base + 1);
}

/* Here is a particularly nasty routine.  It tries to detect whether the distribution */
/* of big and little words is "pathological".  We want the distribution to be random. */
/* If, for example, there are an equal number of big words and little words then the */
/* every other FFT word consists of big word * big word products, while the other half */
/* contains big word * small word products.  This greatly increases the round off error */
/* especially when b is large (big words are much larger than small words).  This */
/* ugliness was added to handle these cases that where the wrong FFT length was selected: */
/* 211*210^2047-1, 211*210^2687-1, 211*210^7679-1.  There are undoubtedly many others. */

int is_pathological_distribution (
	unsigned long num_big_words,
	unsigned long num_small_words)
{
	unsigned long total_length;

/* Handle cases that we really should never see (rational FFTs) */

	if (num_big_words == 0 || num_small_words == 0) return (FALSE);

/* While the remaining number of big words and small words is even, this */
/* represents a case of a big repeating pattern (the pattern in the upper half */
/* of the remaining words is the same as the pattern in the lower half). */

	total_length = num_big_words + num_small_words;
	while ((num_big_words & 1) == 0 && (total_length & 1) == 0) {
		num_big_words >>= 1;
		total_length >>= 1;
	}

/* The bad patterns occur when the number of big words divided by the FFT length */
/* is close to a small rational number like 1/2, 2/5, 3/4, etc.	 We'll define a */
/* pathological bit pattern as one where more than half of the FFT repeats the */
/* same cycle of big words and small words.  This definition may require some */
/* tweaking over time. */

	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 2.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 4.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 3.0 / 4.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 1.0 / 8.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 3.0 / 8.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 5.0 / 8.0)) return (TRUE);
	if (match_pathological_pattern (num_big_words, total_length, 7.0 / 8.0)) return (TRUE);
	if (total_length % 3 == 0) {
		if (match_pathological_pattern (num_big_words, total_length, 1.0 / 3.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 2.0 / 3.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 1.0 / 6.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 5.0 / 6.0)) return (TRUE);
	}
	if (total_length % 5 == 0) {
		if (match_pathological_pattern (num_big_words, total_length, 1.0 / 5.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 2.0 / 5.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 3.0 / 5.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 4.0 / 5.0)) return (TRUE);
	}
	if (total_length % 7 == 0) {
		if (match_pathological_pattern (num_big_words, total_length, 1.0 / 7.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 2.0 / 7.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 3.0 / 7.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 4.0 / 7.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 5.0 / 7.0)) return (TRUE);
		if (match_pathological_pattern (num_big_words, total_length, 6.0 / 7.0)) return (TRUE);
	}

/* That's all the cases we test for now */

	return (FALSE);
}

/* Determine the "bif" value we will look for.  Often this is a straight-forward mapping from */
/* the CPU_ARCHITECTURE.  However, for some CPU architectures, like Pentium M and Core Solo, we */
/* don't have jmptable entries detailing the fastest FFT implementations for those architectures. */

/* BIF values copied from mult.asm */
#define BIF_CORE2		0	// Core 2 CPUs, big L2 caches
#define BIF_CORE2_512		1	// Core 2 Celerons, 512K L2 cache, no L3 cache
#define BIF_I7			2	// Core i3/5/7 CPUs, 256K L2, big L3 caches
#define BIF_FMA3		3	// Core i3/5/7 CPUs - Haswell architecture with FMA3 support
#define BIF_P4_1024		4	// Pentium 4, 1MB cache
#define BIF_P4TP_512		5	// Pentium 4, 512K cache (did not support EMT64)
#define BIF_P4TP_256		6	// Pentium 4, 256K cache (did not support EMT64)
#define BIF_SKX			7	// Intel CPUs that support AVX-512 (Skylake-X)
#define BIF_K8			8	// AMD K8 CPUs
#define BIF_K10			9	// AMD K10 CPUs
#define BIF_RYZEN		11	// AMD Ryzen CPUs

/* Architecture values copied from mult.asm */
#define ARCH_P4TP		1	/* Ancient CPU - architecture value should be retired */
#define ARCH_P4			2
#define ARCH_CORE		3
#define ARCH_FMA3		4
#define ARCH_K8			5
#define ARCH_K10		6
#define ARCH_SKX		8
#define ARCH_BLEND		0

int calculate_bif (
	gwhandle *gwdata,	/* Gwnum global data */
	unsigned long fftlen)
{
	int	cpu_arch, retval;

/* If the architecture and CPU flags are inconsistent, correct the architecture.  This shouldn't */
/* ever happen unless the results from CPUID are overridden. */

	cpu_arch = CPU_ARCHITECTURE;
	if (cpu_arch == CPU_ARCHITECTURE_AMD_K8 && ! (gwdata->cpu_flags & CPU_3DNOW_PREFETCH))
		cpu_arch = CPU_ARCHITECTURE_CORE_2;
	if (cpu_arch == CPU_ARCHITECTURE_AMD_K10 && ! (gwdata->cpu_flags & CPU_3DNOW_PREFETCH))
		cpu_arch = CPU_ARCHITECTURE_CORE_2;
	if (cpu_arch == CPU_ARCHITECTURE_AMD_BULLDOZER && ! (gwdata->cpu_flags & CPU_3DNOW_PREFETCH))
		cpu_arch = CPU_ARCHITECTURE_CORE_2;

/* Map the CPU architecture as determined by CPUID to one of the CPU architectures */
/* that the FFT assembly code is optimized for. */

	switch (cpu_arch) {
	case CPU_ARCHITECTURE_PENTIUM_M:	/* Not sure what is best for these three architectures */
	case CPU_ARCHITECTURE_CORE:
	case CPU_ARCHITECTURE_ATOM:
		retval = BIF_P4_1024;		/* Look for FFTs optimized for large cache P4s */
		break;
	case CPU_ARCHITECTURE_PENTIUM_4:
		if (CPU_L2_CACHE_SIZE <= 128)	/* We haven't optimized for these yet */
			retval = BIF_P4TP_256;	/* Look for FFTs optimized for 256K cache P4s */
		else if (CPU_L2_CACHE_SIZE <= 256)
			retval = BIF_P4TP_256;	/* Look for FFTs optimized for 256K cache P4s */
		else if (CPU_L2_CACHE_SIZE <= 512)
			retval = BIF_P4TP_512;	/* Look for FFTs optimized for 512K cache P4s */
		else if (gwdata->cpu_flags & CPU_TLB_PRIMING)
			retval = BIF_P4TP_512;	/* Look for FFTs optimized for 512K cache P4s */
		else
			retval = BIF_P4_1024;	/* Look for FFTs optimized for large cache P4s */
		break;
	case CPU_ARCHITECTURE_CORE_2:
		if (CPU_L2_CACHE_SIZE <= 1024)
			retval = BIF_CORE2_512;	/* Look for FFTs optimized for Core 2 Celerons */
		else
			retval = BIF_CORE2;	/* Look for FFTs optimized for Core 2 */
		break;
	case CPU_ARCHITECTURE_CORE_I7:
		if (gwdata->cpu_flags & CPU_AVX512F)
			retval = BIF_SKX;	/* Look for Intel-optimized AVX-512 FFT. */
		else if (gwdata->cpu_flags & CPU_FMA3)
			retval = BIF_FMA3;	/* Look for FFTs optimized for Haswell CPUs with FMA3 support */
		else
			retval = BIF_I7;	/* Look for FFTs optimized for Core i3/i5/i7 */
		break;
	case CPU_ARCHITECTURE_PHI:		/* Intel's Xeon Phi CPUs */
		if (gwdata->cpu_flags & CPU_AVX512F)
			retval = BIF_SKX;	/* Look for Intel-optimized AVX-512 FFT. */
		else if (gwdata->cpu_flags & CPU_FMA3)
			retval = BIF_FMA3;	/* Look for FFTs optimized for Haswell CPUs with FMA3 support */
		else
			retval = BIF_I7;	/* Look for FFTs optimized for Core i3/i5/i7 */
		break;
	case CPU_ARCHITECTURE_INTEL_OTHER:	/* This is probably one of Intel's next generation CPUs */
		if (gwdata->cpu_flags & CPU_AVX512F)
			retval = BIF_SKX;	/* We don't know which AVX-512 FFT is fastest.  Try this one. */
		else if (gwdata->cpu_flags & CPU_FMA3)
			retval = BIF_FMA3;	/* Look for FFTs optimized for Haswell CPUs with FMA3 support */
		else
			retval = BIF_I7;	/* Look for FFTs optimized for Core i3/i5/i7 */
		break;
	case CPU_ARCHITECTURE_AMD_K8:
		retval = BIF_K8;		/* Look for FFTs optimized for K8 */
		break;
	case CPU_ARCHITECTURE_AMD_K10:
		retval = BIF_K10;		/* Look for FFTs optimized for K10 */
		break;
	case CPU_ARCHITECTURE_AMD_BULLDOZER:	/* Bulldozer is horrible at AVX.  Gwinit turns off AVX & FMA3 to get K10 optimized. */
		if (gwdata->cpu_flags & CPU_FMA3)  /* Should only happen during torture test */
			retval = BIF_FMA3;	/* FMA3 FFTs */
		else if (gwdata->cpu_flags & CPU_AVX)  /* Should only happen during torture test */
			retval = BIF_I7;	/* AVX without FMA3 FFTs */
		else
			retval = BIF_K10;	/* Look for FFTs optimized for K10 */
		break;
	case CPU_ARCHITECTURE_AMD_ZEN:		/* Look for FFTs optimized for Ryzen */
		if (! (gwdata->cpu_flags & CPU_FMA3))
			retval = BIF_I7;	/* Shouldn't happen */
		else
			retval = BIF_RYZEN;
		break;
	case CPU_ARCHITECTURE_AMD_OTHER:	/* For no particularly good reason, assume future AMD processors do well with Intel FFTs */
		if (! (gwdata->cpu_flags & CPU_FMA3))
			retval = BIF_I7;	/* Look for FFTs optimized for Core i3/i5/i7 */
		else
			retval = BIF_FMA3;	/* Look for FFTs optimized for Intel FMA3 CPUs */
		break;
	case CPU_ARCHITECTURE_OTHER:		/* Probably a VIA processor */
		if (gwdata->cpu_flags & CPU_AVX512F)
			retval = BIF_SKX;	/* We don't know which AVX-512 FFT is fastest.  Try this one. */
		else if (gwdata->cpu_flags & CPU_FMA3)
			retval = BIF_FMA3;	/* Look for FFTs optimized for Haswell CPUs with FMA3 support */
		else if (gwdata->cpu_flags & CPU_AVX)
			retval = BIF_I7;	/* Look for FFTs optimized for Core i7 (AVX without FMA3) */
		else
			retval = BIF_CORE2;	/* We don't know which FFT is fastest.  For no particularly good reason, use Core 2 FFTs */
		break;
	case CPU_ARCHITECTURE_PRE_SSE2:		/* Cannot happen, gwinfo should have selected x87 FFTs */
	default:				/* For no particularly good reason, look for FFTs optimized for Core 2 */
		retval = BIF_CORE2;
		break;
	}

/* We never implemented FMA3 FFTs for 32-bit builds (lack of YMM registers). */

#ifndef X86_64
	if (retval == BIF_FMA3 || retval == BIF_RYZEN) retval = BIF_I7;
#endif

/* For slower CPU architectures we didn't bother to find the best FFT implementation */
/* for the larger FFTs.  This was done to reduce the size of the executable.  If we */
/* are asked to run one of these large FFTs, select an FFT optimized for a different */
/* CPU architecture. */

	if (fftlen > 1572864 && retval == BIF_P4TP_256)
		retval = BIF_P4_1024;	/* Tiny cache P4s have best FFT implementations up to 1536K */
	if (fftlen > 4194304 && retval == BIF_P4TP_512)
		retval = BIF_P4_1024;	/* Small cache P4s have best FFT implementations up to 4M */
	if (fftlen > 6291456 && retval == BIF_P4_1024)
		retval = BIF_CORE2;	/* P4s have best FFT implementations up to 6M */
	if (fftlen > 6291456 && retval == BIF_K8)
		retval = BIF_K10;	/* K8s have best FFT implementations up to 6M */
	if (fftlen > 4194304 && retval == BIF_CORE2_512)
		retval = BIF_CORE2;	/* Small cache Core 2 Celerons have best FFT implementations up to 4M */

/* Return the result */

	return (retval);
}

/* Ugly little macros to bump jmptable pointer to next procedure entry or to next count */
static __inline const struct gwasm_jmptab *INC_JMPTAB (const struct gwasm_jmptab *x)
{
	if (x->flags & 0x40000000) return ((const struct gwasm_jmptab *) ((const char *)(x) + sizeof(uint32_t) + 2*sizeof(void*) + sizeof(uint32_t)));
	return ((const struct gwasm_jmptab *) ((const char *)(x) + sizeof(uint32_t) + sizeof(void*) + sizeof(uint32_t)));
}
static __inline const struct gwasm_jmptab *LAST_JMPTAB (const struct gwasm_jmptab *x)
{
	const struct gwasm_jmptab *next;
	for ( ; (next = INC_JMPTAB (x))->flags & 0x80000000; x = next);
	return (x);
}
static __inline const struct gwasm_jmptab *NEXT_SET_OF_JMPTABS (const struct gwasm_jmptab *x)
{
	// Handle case where there are no FFT implementations
	if (x->flags == 0) return ((const struct gwasm_jmptab *) &x->proc_ptr);
	// Get pointer to the last FFT implementation
	x = LAST_JMPTAB (x);
	// Adjust for jmptab containing two proc ptrs
	if (x->flags & 0x40000000) x = (const struct gwasm_jmptab *) ((const char *)(x) + sizeof(void*));
	// Skip non-zero counts
	while (x->counts[0]) x = (const struct gwasm_jmptab *) ((const char *)(x) + sizeof(int32_t));
	// Next set of jmptabs is after the zero count
	return ((const struct gwasm_jmptab *) &x->counts[1]);
}

/* This routine checks to see if there is an FFT implementation for this FFT length and */
/* CPU architecture.  For example, when the FFT length is just less than a power of two, on */
/* some CPUs it may be better to use the larger power-of-two FFT length and thus there */
/* will not be an FFT implementation for this slightly smaller FFT length. */ 

int is_fft_implemented (
	gwhandle *gwdata,		/* Gwnum global data */
	int	all_complex,		/* TRUE if this jmptab entry if from the all-complex FFT table */
	const struct gwasm_jmptab *jmptab, /* Jmptable entry from mult.asm to examine */
	int	*best_impl_id)		/* Returned ID of the best FFT implementation as determined by the benchmark database */
{
	int	desired_bif;		/* The "best implementation for" value we will look for. */
					/* See mult.asm for defined BIF_ values. */

/* If there are no implementations, return false.  This can happen for SSE2 and earlier implementations where there are no best-impl-for */
/* settings or only for an architecture (K8) we can't run.  Examples include SSE2 FFT lengths 2000K and 3456K in 32-bit mode. */

	if (! (gwdata->cpu_flags & CPU_3DNOW_PREFETCH)) {
		for ( ; jmptab->flags & 0x80000000; jmptab = INC_JMPTAB (jmptab)) {
			int	arch = (jmptab->flags >> 17) & 0xF;
			if (arch != ARCH_K8 && arch != ARCH_K10) break;
		}
	}
	if (!(jmptab->flags & 0x80000000)) return (FALSE);

/* Determine the "bif" value we will look for.  Often this is a straight-forward mapping from */
/* the CPU_ARCHITECTURE.  However, for some CPU architectures, like Pentium M and Core Solo, we */
/* don't have jmptable entries detailing the fastest FFT implementations for those architectures. */

	desired_bif = calculate_bif (gwdata, jmptab->fftlen);

/* Assume we will use the old school method (no benchmark database) where the */
/* best FFT implementation is hardwired into the mult.asm jmptable. */

	*best_impl_id = -1;

/* For small SSE2 FFTs as well as all x87 FFTs there is only one implementation -- use it. */

	if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))) {
		if (gwdata->cpu_flags & CPU_SSE2) {
			if (jmptab->fftlen < 7168) return (TRUE);
		} else
			return (TRUE);
	}

/* If we are benchmarking all FFT implementations or we are doing QA, then we want to test */
/* this FFT length even if it isn't optimal */

	if (gwdata->bench_pick_nth_fft || gwdata->qa_pick_nth_fft) return (TRUE);

/* If we are using benchmark data to select fastest FFT implementation check the benchmark database to see */
/* if there is benchmark data available and that a larger FFT size will not be faster. */

	if (gwdata->use_benchmarks) {
		int	arch, num_cores, num_workers, num_hyperthreads, i;

/* Calculate the architecture, number of cores, workers, and hyperthreads for looking up the right benchmark data */

		switch (desired_bif) {			/* Only get benchmark data for relevant CPU architecture */
		case BIF_SKX:
			arch = ARCH_SKX;
			break;
		case BIF_FMA3:
		case BIF_RYZEN:
			arch = ARCH_FMA3;
			break;
		case BIF_I7:
		case BIF_CORE2:
		case BIF_CORE2_512:
			arch = ARCH_CORE;
			break;
		case BIF_K10:
			arch = ARCH_K10;
			break;
		case BIF_K8:
			arch = ARCH_K8;
			break;
		default:
			arch = ARCH_P4;
		}
		if (gwdata->will_hyperthread) {
			num_hyperthreads = gwdata->will_hyperthread;		/* User can tell us more than 2 hyperthreads will be used */
			if (num_hyperthreads <= 1) num_hyperthreads = 2;	/* Minimum hyperthread count is two */
		} else
			num_hyperthreads = 1;
		num_cores = gwdata->bench_num_cores;				/* Use suggested value from caller */
		if (num_cores == 0) num_cores = CPU_CORES;			/* Else default to all cores will be busy */
		num_workers = gwdata->bench_num_workers;			/* Use suggested value from caller */
		if (num_workers == 0) num_workers = num_cores * num_hyperthreads / gwdata->num_threads; /* Else default worker count */

/* See if the benchmark database has bench data either with or without error checking. */
/* Once we have throughput data from the benchmark database, make sure a slightly larger */
/* FFT length will not offer even more throughput. */

		for (i = 0; i <= 1; i++) {
			const struct gwasm_jmptab *next_jmptab;
			int	error_check, no_r4dwpn, impl, next_impl;
			double	throughput, next_throughput;

			if (gwdata->will_error_check == 0) error_check = i;	/* Look for no-error-checking benchmarks first */
			if (gwdata->will_error_check == 1) error_check = !i;	/* Look for error-checking benchmarks first */
			if (gwdata->will_error_check == 2) error_check = i;	/* Need a more sophisticated approach in this case */
			no_r4dwpn = (gwdata->sum_inputs_checking && ! gwdata->ALL_COMPLEX_FFT && ! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX)));
			gwbench_get_max_throughput (jmptab->fftlen, arch, num_cores, num_workers, num_hyperthreads,
						    all_complex, error_check, no_r4dwpn, &impl, &throughput);
			if (throughput <= 0.0) continue;

			for (next_jmptab = NEXT_SET_OF_JMPTABS(jmptab); ; next_jmptab = NEXT_SET_OF_JMPTABS(next_jmptab)) {
				if (next_jmptab->fftlen == 0 ||				/* There is no next FFT length */
				    next_jmptab->fftlen > 1.03 * jmptab->fftlen) {	/* Next FFT length is much bigger (and therefore slower) */
					*best_impl_id = impl;
					return (TRUE);
				}
				gwbench_get_max_throughput (next_jmptab->fftlen, arch, num_cores, num_workers, num_hyperthreads,
							    all_complex, error_check, no_r4dwpn, &next_impl, &next_throughput);
				if (next_throughput <= 0.0) break;		/* We don't have enough bench data to know if slightly larger FFT length will be faster */
				if (next_throughput > throughput) return (FALSE); /* Larger FFT length is faster */
			}
		}
	}

/* Loop through the FFT implementations to see if we find an implementation that matches our desired "bif" value. */

	while (jmptab->flags & 0x80000000) {
		if (((jmptab->flags >> 13) & 0xF) == desired_bif) return (TRUE);
		jmptab = INC_JMPTAB (jmptab);
	}

/* FFT implementation not found.  A larger FFT length should be faster. */

	return (FALSE);
}

/* Some gwinfo calculations depend on whether this is a one-pass or two-pass FFT. */
/* However, some FFTs have both a one-pass and two-pass implementation.  In such cases, */
/* we must make sure that jmptab points to the implementation that we will actually use. */

const struct gwasm_jmptab *choose_one_pass_or_two_pass_impl (
	gwhandle *gwdata,	/* Gwnum global data */
	const struct gwasm_jmptab *jmptab)
{
	const struct gwasm_jmptab *orig_jmptab;	
	int	desired_bif;		/* The "best implementation for" value we will look for. */
					/* See mult.asm for defined BIF_ values. */

/* For small SSE2 FFTs as well as all x87 FFTs, there is one implementation and it is always available */

	if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))) {
		if (gwdata->cpu_flags & CPU_SSE2) {
			if (jmptab->fftlen < 7168) return (jmptab);
		} else
			return (jmptab);
	}

/* If we are benchmarking all FFT implementations, then we always want to start with the first implementation */

	if (gwdata->bench_pick_nth_fft) return (jmptab);

/* If the first entry is a two-pass implementation, then all the implementations are two-pass */
/* Use the first entry, as any of the entries will do for our purposes */

	if ((jmptab->flags & 0x1FF) != 0) return (jmptab);

/* Determine the "bif" value we will look for.  Often this is a straight-forward mapping from */
/* the CPU_ARCHITECTURE.  However, for some CPU architectures, like Pentium M and Core Solo, we */
/* don't have jmptable entries detailing the fastest FFT implementations for those architectures. */

	desired_bif = calculate_bif (gwdata, jmptab->fftlen);

/* Loop through the FFT implementations to see if we find a one-pass implementation */
/* that matches our desired "bif" value.  Otherwise, return first two-pass implementation. */

	orig_jmptab = jmptab;
	while (jmptab->flags & 0x80000000) {
		if ((jmptab->flags & 0x1FF) != 0) return (jmptab);
		if (((jmptab->flags >> 13) & 0xF) == desired_bif) return (orig_jmptab);
		jmptab = INC_JMPTAB (jmptab);
	}

/* FFT implementation not found.  Can't happen as is_fft_implemented should have been called. */

	return (orig_jmptab);
}

/* FMA3 FFTs compared to AVX FFTs have on average 7% less round off error.  This lets us have a higher max exponent. */
/* Increase max_exp by (1-SQRT(.93)) * FFTlen.  We'll be a little more conservative for smaller FFT lengths. */

unsigned long adjusted_max_exponent (
	const gwhandle *gwdata,			/* Gwnum global data */
	const struct gwasm_jmptab *jmptab)	/* Jmptable entry */
{
	unsigned long max_exp;

	max_exp = jmptab->max_exp;
	if (max_exp == 0) return (0);
	if (!(gwdata->cpu_flags & CPU_AVX512F) && gwdata->cpu_flags & CPU_FMA3) {
		if (jmptab->fftlen >= 8192) max_exp += (int) (0.035635 * (double) jmptab->fftlen);
		else if (jmptab->fftlen >= 2048) max_exp += (int) (0.017656 * (double) jmptab->fftlen);
	}
	return (max_exp);
}

/* This routine used to be in assembly language.  It scans the assembly */
/* code arrays looking for the best FFT size to implement our k*b^n+c FFT. */
/* Returns 0 for IBDWT FFTs, 1 for zero padded FFTs, or a gwsetup error */
/* code. */

int gwinfo (			/* Return zero-padded fft flag or error code */
	gwhandle *gwdata,	/* Gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* N in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	struct gwinfo1_data asm_info;
	const struct gwasm_jmptab *jmptab, *zpad_jmptab, *generic_jmptab, *impl_jmptab;
	int	best_impl_id, zpad_best_impl_id;
	double	log2k, logbk, log2b, log2c, log2maxmulbyconst;
	double	max_bits_per_input_word, max_bits_per_output_word;
	double	max_weighted_bits_per_output_word;
	int	num_b_in_big_word, num_small_words, num_big_words;
	double	b_per_input_word, bits_per_output_word;
	double	weighted_bits_per_output_word;
	unsigned long max_exp;
	char	buf[20];
	int	qa_nth_fft, desired_bif;
	void	*prev_proc_ptrs[5];
	uint32_t flags;
	struct gwasm_data *asm_data;

/* Get pointer to 6 assembly jmptables and the version number */

	gwinfo1 (&asm_info);

/* Make sure that the assembly code version number matches the C version */
/* number.  If they do not match, then the user linked in the wrong gwnum */
/* object files! */

	sprintf (buf, "%d.%d", asm_info.version / 100, asm_info.version % 100);
	if (strcmp (buf, GWNUM_VERSION)) return (GWERROR_VERSION);

/* Precalculate some needed values */

	log2k = log2 (k);
	logbk = logb (k);
	log2b = log2 (b);
	log2c = log2 (labs (c));
	log2maxmulbyconst = log2 (gwdata->maxmulbyconst);

/* The smallest AVX-512F FFT is length 1K.  Limited testing has indicated that less than 5 bits per FFT word */
/* can result in carry propagation errors.  For small k*b^n+c values switch to not using AVX-512 instructions. */
/* Don't catch the n==0 case from gwmap_with_cpu_flags_fftlen_to_max_exponent. */

	if (n && log2b * (double) n < 5.0 * 1024.0)
		gwdata->cpu_flags &= ~CPU_AVX512F;

/* First, see what FFT length we would get if we emulate the k*b^n+c modulo */
/* with a zero padded FFT.  If k is 1 and abs (c) is 1 then we can skip this */
/* loop as we're sure to find an IBDWT that will do the job. Also skip if called from */
/* gwmap_fftlen_to_max_exponent (n = 0) or we are QAing IBDWT FFTs (qa_pick_nth_fft >= 1000) */

again:	zpad_jmptab = NULL;
	generic_jmptab = NULL;
	if (! gwdata->force_general_mod && (k > 1.0 || (n > 0 && n < 500) || labs (c) > 1) && gwdata->qa_pick_nth_fft < 1000) {

/* Use the proper 2^N-1 jmptable */

		if (gwdata->cpu_flags & CPU_AVX512F)
			zpad_jmptab = asm_info.avx512_cyclic_fft_info;
		else if (gwdata->cpu_flags & CPU_AVX)
			zpad_jmptab = asm_info.avx_cyclic_fft_info;
		else if (gwdata->cpu_flags & CPU_SSE2)
			zpad_jmptab = asm_info.sse2_cyclic_fft_info;
		else
			zpad_jmptab = asm_info.x86_cyclic_fft_info;

/* If no jmptable was found (should only happen if caller has disabled all but the x87 CPU_FLAGS in 64-bit mode */

		if (zpad_jmptab == NULL) return (GWERROR_TOO_LARGE);

/* Find the table entry for the FFT that can do a mod 2^2n FFT, handling */
/* k and c in the normalization routines.  We will compare this to the */
/* non-zero-padded FFT length later.  The zeroes in the upper half of FFT */
/* input data let us get about another 0.3 bits per input word. */

		while ((max_exp = adjusted_max_exponent (gwdata, zpad_jmptab)) != 0) {

/* Do a quick check on the suitability of this FFT */

			if ((double) (n + n) * log2b / (double) zpad_jmptab->fftlen > 27.0 - 0.25 * log2 (zpad_jmptab->fftlen)) goto next1;
			if (zpad_jmptab->fftlen < gwdata->minimum_fftlen) goto next1;

/* Don't bother looking at this FFT length if the generic reduction would be faster */

			if (generic_jmptab != NULL && zpad_jmptab->timing > 3.0 * generic_jmptab->timing) goto next1;

/* Make sure this FFT length is implemented and benchmarking does not show that a larger FFT will be faster */

			if (! is_fft_implemented (gwdata, FALSE, zpad_jmptab, &zpad_best_impl_id)) goto next1;

/* See if this is the FFT length that would be used for a generic modulo reduction */

			if (generic_jmptab == NULL &&
			    2.0 * (log2k + n * log2b) + 128.0 < max_exp + 0.3 * zpad_jmptab->fftlen)
				generic_jmptab = zpad_jmptab;

/* Compare the maximum number of bits allowed in the FFT input word */
/* with the number of bits we would use.  Break when we find an acceptable */
/* FFT length. */
//  This is the old code which only supported b == 2
//			max_bits_per_word = (double) max_exp / zpad_jmptab->fftlen;
//			max_bits_per_word -= gwdata->safety_margin;
//			bits_per_word = (double) (n + n) * log2b / zpad_jmptab->fftlen;
//			if (bits_per_word < max_bits_per_word + 0.3) {
//				break;
//			}

/* In version 25.11, we need to handle b != 2.  See comments later on in this routine */
/* for a description of the concepts involved. */

/* Compute the maximum number of bits allowed in the FFT input word */

			max_bits_per_input_word = (double) max_exp / zpad_jmptab->fftlen;
			max_bits_per_input_word -= gwdata->safety_margin;

/* Apply our new formula (described later) to the maximum Mersenne exponent for this FFT length. */

			num_b_in_big_word = (int) ceil (max_bits_per_input_word);
			num_small_words = (int) ((num_b_in_big_word - max_bits_per_input_word) * zpad_jmptab->fftlen);
			num_big_words = zpad_jmptab->fftlen - num_small_words;
			max_bits_per_output_word =
				2 * (num_b_in_big_word - 1) +
				0.6 * log2 (num_big_words + num_small_words / 3.174802103936252);

/* Apply our new formula (described later) to the number we are planning to test.  */
/* This is different for the zero-pad case because only 4 words in the upper half */
/* of the FFT contain any data.  We can't use the FFT length if the k value will */
/* not fit in 4 words. */

			b_per_input_word = (double) (n + n) / zpad_jmptab->fftlen;
			if (logbk > 4.0 * b_per_input_word) goto next1;
			num_b_in_big_word = (int) ceil (b_per_input_word);
			num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * (zpad_jmptab->fftlen / 2 + 4));
			num_big_words = (zpad_jmptab->fftlen / 2 + 4) - num_small_words;
			bits_per_output_word =
				2.0 * (num_b_in_big_word * log2b - 1.0) +
				0.6 * log2 (num_big_words + num_small_words / pow (2.0, log2b / 0.6));

/* And compute the weighted values as per the formulas described later.  In 29.7 we added the log2b < 12.5 case (to match */
/* the non-zero-padded code) as 27904^53415-7 was choosing a zero-padded 200K AVX-512 FFT and getting roundoffs frequently */
/* above 0.4 and rarely above 0.5.  We need to do a careful analysis using average round error to fine tune this adjustment further. */

			max_weighted_bits_per_output_word =
				2.0 * max_bits_per_input_word + 0.6 * log2 (zpad_jmptab->fftlen / 2 + 4);
			weighted_bits_per_output_word =
				2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
				0.6 * log2 (zpad_jmptab->fftlen / 2 + 4);
			if ((n + n) % zpad_jmptab->fftlen == 0)
				weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);
			else if (! is_pathological_distribution (num_big_words, num_small_words))
				weighted_bits_per_output_word -=
					((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
					 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
					 (log2b <= 12.5) ? 2.0 + (log2b - 6.0) / 6.5 : 3.0);

/* See if this FFT length might work */

			if ((weighted_bits_per_output_word <= max_weighted_bits_per_output_word || gwdata->minimum_fftlen) &&

/* Result words are multiplied by k and the mul-by-const and any carry spread over 6 words. */
/* Thus, the multiplied FFT result word cannot be more than 7 times bits-per-input-word */
/* (bits-per-input-word are stored in the current word and the 6 words we propagate carries to). */

			    bits_per_output_word + log2k + log2maxmulbyconst <= floor (7.0 * b_per_input_word) * log2b &&

/* The high part of upper result words are multiplied by c and the mul-by-const.  This must not exceed 51 bits. */

			    bits_per_output_word - floor (b_per_input_word) * log2b + log2c + log2maxmulbyconst <= 51.0) {

/* We can use this FFT.  Look for a non-zero-padded FFT that might be even faster. */

				break;
			}

/* Move past procedure entries and counts to next jmptable entry */

next1:			zpad_jmptab = NEXT_SET_OF_JMPTABS (zpad_jmptab);
		}

/* Some k/b/n/c values can't be handled by AVX-512 FFTs because there are relatively few small FFT sizes available. */
/* In these cases, we'll retry using FMA3 FFTs.  One example is 15312853462553*2^5257+1.  NOTE:  There may be non-base2 cases */
/* where we could use a zero-padded FMA3 FFT, but instead we'll end up selecting a generic reduction AVX-512 FFT. */

		if (zpad_jmptab->max_exp == 0 && gwdata->cpu_flags & CPU_AVX512F && b == 2 && n < 100000) {
			gwdata->cpu_flags &= ~CPU_AVX512F;
			goto again;
		}
	}

/* Now see what FFT length we would use if a DWT does the k*b^n+c modulo. */

/* Use the proper 2^N+1 or 2^N-1 jmptable */

	if (c < 0) {
		if (gwdata->cpu_flags & CPU_AVX512F)
			jmptab = asm_info.avx512_cyclic_fft_info;
		else if (gwdata->cpu_flags & CPU_AVX)
			jmptab = asm_info.avx_cyclic_fft_info;
		else if (gwdata->cpu_flags & CPU_SSE2)
			jmptab = asm_info.sse2_cyclic_fft_info;
		else
			jmptab = asm_info.x86_cyclic_fft_info;
	} else {
		if (gwdata->cpu_flags & CPU_AVX512F)
			jmptab = asm_info.avx512_negacyclic_fft_info;
		else if (gwdata->cpu_flags & CPU_AVX)
			jmptab = asm_info.avx_negacyclic_fft_info;
		else if (gwdata->cpu_flags & CPU_SSE2)
			jmptab = asm_info.sse2_negacyclic_fft_info;
		else
			jmptab = asm_info.x86_negacyclic_fft_info;
	}

/* If no jmptable was found (should only happen if caller has disabled all but the x87 CPU_FLAGS in 64-bit mode */

	if (jmptab == NULL) return (GWERROR_TOO_LARGE);

/* Find the table entry using either the specified fft length or */
/* the smallest FFT length that can handle the k,b,n,c being tested. */

	while ((max_exp = adjusted_max_exponent (gwdata, jmptab)) != 0) {

/* Do a quick check on the suitability of this FFT */

		if ((double) n * log2b / (double) jmptab->fftlen > 26.0 - 0.25 * log2 (jmptab->fftlen)) goto next2;
		if (jmptab->fftlen < gwdata->minimum_fftlen) goto next2;

/* Top carry adjust can only handle k values of 34 bits or less */

		if (log2k >= 34.0) goto next2;

/* Make sure this FFT length is implemented and benchmarking does not show that a larger FFT will be faster */

		if (! is_fft_implemented (gwdata, c > 0, jmptab, &best_impl_id)) goto next2;

/* Always use the minimum_fftlen if n is zero (a special case call from gwmap_fftlen_to_max_exponent) */

		if (gwdata->minimum_fftlen && n == 0) break;

/* Some calculations below depend on whether this is a one-pass or two-pass FFT. */
/* However, some FFTs have both a one-pass and two-pass implementation.  In such cases, */
/* we must make sure that jmptab points to the implementation that we will actually use. */

		impl_jmptab = choose_one_pass_or_two_pass_impl (gwdata, jmptab);

/* Check if this FFT length will work with this k,n,c combo */

//  This is the old code which only supported b == 2
//		double max_bits_per_word;
//		double bits_per_word;
//
/* Compute the maximum number of bits allowed in the FFT input word */
//
//		max_bits_per_word = (double) max_exp / jmptab->fftlen;
//		max_bits_per_word -= gwdata->safety_margin;
//
/* For historical reasons, the jmptable computes maximum exponent based on */
/* a Mersenne-mod FFT (i.e k=1.0, c=-1).  Handle more complex cases here. */
/* A Mersenne-mod FFT produces 2 * bits_per_word in each FFT result word. */
/* The more general case yields 2 * bits_per_word + log2(k) + 1.7 * log2(c) */
/* in each FFT result word.  NOTE: From the data I've observed, doubling c */
/* about triples the roundoff error (that is log2(3) = 1.585 output bits). */
/* However, when I used 1.585 in the formula it was not hard to find cases */
/* where the roundoff error was too high, so we use 1.7 here for extra */
/* safety. */
//
//		bits_per_word = (log2k + n * log2b) / jmptab->fftlen;
//		if (2.0 * bits_per_word + log2k + 1.7 * log2c <=
//					2.0 * max_bits_per_word) {
/* Because carries are spread over 4 words, there is a minimum limit on */
/* the bits per word.  An FFT result word cannot be more than 5 times */
/* bits-per-word (bits-per-word are stored in the current word and the */
/* 4 words we propagate carries to).  How many bits are in an FFT result */
/* word?  Well, because of balanced representation the abs(input word) is */
/* (bits_per_word-1) bits long. An FFT result word contains multiplied data */
/* words, that's (bits_per_word-1)*2 bits.  Adding up many multiplied data */
/* words adds some bits proportional to the size of the FFT.  Experience */
/* has shown this to be 0.6 * log (FFTLEN).  This entire result is */
/* multiplied by k in the normalization code, so add another log2(k) bits. */
/* Finally, the mulbyconst does not affect our chance of getting a round off */
/* error, but does add to the size of the carry. */
//
//		loglen = log2 (jmptab->fftlen);
//		total_bits = (bits_per_word - 1.0) * 2.0 +
//				     1.7 * log2c + loglen * 0.6 +
//				     log2k + log2maxmulbyconst;
//		if (total_bits > 5.0 * bits_per_word) {

/* In version 25.11, we now need to handle b != 2.  Consider the case */
/* where b is ~4000.  If small words contain one b (~12 bits) and large words */
/* contain two b (~24 bits), then the number of bits in a result word is */
/* dominated by big words * big word products (~48 bits).  The old code above */
/* tested average bits per word (~18 bits) and underestimates a result word as */
/* containing ~36 bits.  So here's our upgraded model.  We calculate the number */
/* of big and little words.  A result word adds up many big word times big word */
/* products and big word times small word products.  Let base = b ^ num_b_in_big_word. */
/* Because of balanced representation, a big word times big word */
/* product is 2 * (log2(base) - 1) bits.  Summing them up adds about */
/* 0.6 * log2 (num_big_words) more bits.  Now for the big words times */
/* small words products that are also added in, the more bits in a small word the more */
/* it impacts the result word.  A big word times small word product has log2(b) fewer */
/* bits in it.  If we add two big word times small word products, the sum is */
/* about 0.6 bits bigger, add four products to get 1.2 bits bigger, etc. -- do this until you */
/* overcome the log2(b) bit difference.  That is, 2^(log2(b)/0.6) small */
/* products equals one big product.  Putting it all together, a result word contains */
/* 2 * (log2(base) - 1) + 0.6 * log2 (num_big_words + num_small_words / 2^(log2(b)/0.6)) */
/* bits plus the k and c adjustments noted above. */

/* Compute the maximum number of bits allowed in the FFT input word */

		max_bits_per_input_word = (double) max_exp / jmptab->fftlen;
		max_bits_per_input_word -= gwdata->safety_margin;

/* Apply our new formula above to the maximum Mersenne exponent for this FFT length. */

		num_b_in_big_word = (int) ceil (max_bits_per_input_word);
		num_small_words = (int) ((num_b_in_big_word - max_bits_per_input_word) * jmptab->fftlen);
		num_big_words = jmptab->fftlen - num_small_words;
		max_bits_per_output_word =
				2 * (num_b_in_big_word - 1) +
				0.6 * log2 (num_big_words + num_small_words / 3.174802103936252);

/* Apply our new formula to the number we are planning to test.  In version 26.3 we changed */
/* "2.0 * (num_b_in_big_word * log2b - 1.0)" to "floor (2.0 * b_per_input_word) * log2b - 2.0" */
/* because of examples like 10024*603^153-1 which has num_b_in_big_word = 1 and b_per_input_word = 0.8. */
/* This means weighted values are as large as b^1.8.  Squaring that large value yields b^3.6.  Apply */
/* the inverse weight of b^-.6 and we have output values as large as b^3.  The improved formula */
/* reflects these larger output values. */

		b_per_input_word = (logbk + n) / jmptab->fftlen;
		num_b_in_big_word = (int) ceil (b_per_input_word);
		num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * jmptab->fftlen);
		num_big_words = jmptab->fftlen - num_small_words;
		if (k == 1.0 && n % jmptab->fftlen == 0)
			bits_per_output_word = 
				2.0 * (num_b_in_big_word * log2b - 1.0) +
				0.6 * log2 (num_big_words + num_small_words / pow (2.0, log2b / 0.6)) +
				log2k + 1.7 * log2c;
		else
			bits_per_output_word = 
				floor (2.0 * (b_per_input_word + 1.0)) * log2b - 2.0 +
				0.6 * log2 (num_big_words + num_small_words / pow (2.0, log2b / 0.6)) +
				log2k + 1.7 * log2c;

/* Unfortunately, the story does not end there.  The weights applied to each FFT word */
/* range from 1 to b.  These extra bits impact the round off error.  Thus, we calculate */
/* the weighted_bits_per_output_word for irrational FFTs as using another log2b bits. */

		max_weighted_bits_per_output_word = 2.0 * max_bits_per_input_word + 0.6 * log2 (jmptab->fftlen);
		if (k == 1.0 && n % jmptab->fftlen == 0) {
			/* New in 29.7: 1889^20480+1 was generating too large a roundoff using a 10K AVX-512 FFT. */
			/* Testing shows that an irrational FFT has about twice the average roundoff error as a */
			/* rational FFT (equals one bit in the output word).  But comparing 2^204800+1 and 2^204801+1 */
			/* in a 10K FFT, our code calculates a two bit difference between weighted_bits_per_output_word */
			/* and max_weighted_bits_per_output_word.  The proper solution is to correct the calculation of */
			/* max_weighted_bits_per_output_word and irrational weighted_bits_per_output_word.  However, */
			/* I'm a little leery of changing that working code so I'm just adding one to the rational case here. */
			weighted_bits_per_output_word =	bits_per_output_word + 1.0;
		} else {
			weighted_bits_per_output_word =
				2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
				0.6 * log2 (jmptab->fftlen) + log2k + 1.7 * log2c;

/* A pathological case occurs when num_big_words is one and k is greater than one. */
/* The FFT weights for the small words will not range from 1 to b.  Depending on the */
/* fractional part of logb(k).  In the worst case scenario, the small word weights */
/* range from b - epsilon to b.  The example that raised this issue is 28*3^12285-1. */

			if (num_big_words == 1 && k > 1.0)
				weighted_bits_per_output_word += log2b;

/* Furthermore, testing shows us that larger b values don't quite need the full log2b */
/* bits added (except for some pathological cases), probably because there are fewer */
/* extra bits generated by adding products because the smallest weighted words have */
/* fewer bits.  The correction is if log2b is 3 you can get 1 more output bit than */
/* expected, if log2b is 6 you get about 2 extra bits, if log2b is 12 you can get */
/* 3 extra bits.  This formula was found to be a bit too aggressive, at least for */
/* large b.  Two examples:  2*22563^22563-1 and  2*22576^22576-1 fail in a 40K FFT. */
/* Thus, we're changing the formula in v27.9 to log2b of 12.5 gets the max of 3 extra bits. */
/* Also, some examples such as 19464*19^31895+1 and 245*830^492-1 (worst case we */
/* know of) still raise round off errors.  For added safety we assume an extra */
/* 0.3 bits of output are needed when base is not 2. */

			else if (! is_pathological_distribution (num_big_words, num_small_words)) {
				weighted_bits_per_output_word -=
					((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
					 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
					 (log2b <= 12.5) ? 2.0 + (log2b - 6.0) / 6.5 : 3.0);
				if (b != 2) weighted_bits_per_output_word += 0.3;
			}
		}

/* If the bits in an output word is less than the maximum allowed (or the user is trying to force us */
/* to use this FFT), then we can probably use this FFT length -- though we need to do a few more tests. */

		if (weighted_bits_per_output_word <= max_weighted_bits_per_output_word || gwdata->minimum_fftlen) {
			double	carries_spread_over;

/* Originally, carries were spread over 4 FFT words.  Some FFT code has been */
/* upgraded to spread the carry over 6 FFT words.  Handle that here.  Note that */
/* FFT lengths 80 and 112 were not upgraded. */

			if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2)))
				carries_spread_over = 4.0;
			else if (jmptab->fftlen == 80 || jmptab->fftlen == 112)
				carries_spread_over = 4.0;
			else if ((impl_jmptab->flags & 0x1FF) == 0)	// One pass AVX-512/AVX/SSE2 FFTs
				carries_spread_over = 6.0;
			else						// Two pass AVX-512/AVX/SSE2 FFTs
				carries_spread_over = 6.0;
				
/* Because carries are spread over 4 words, there is a minimum value for the bits */
/* per FFT word.  An FFT result word must fit in the floor(bits-per-input-word) stored */
/* in the current word plus ceil (4 * bits-per-input-word) for the carries to */
/* propagate into.  The mul-by-const during the normalization process adds to */
/* the size of the result word. */

			if (bits_per_output_word + log2maxmulbyconst >
					(floor (b_per_input_word) + ceil (carries_spread_over * b_per_input_word)) * log2b) {
// This assert was designed to find any cases where using more carry words
// would use a shorter FFT than using a zero-padded FFT.  There are many such
// cases, especially with larger bases.  One example is 5001*500^100000-1.
// It would use 132K FFT length if I implemented 7 carry words
// and requires a 144K FFT length for a zero-padded FFT.
//				ASSERTG (zpad_jmptab == NULL || jmptab->fftlen >= zpad_jmptab->fftlen);
				goto next2;
			}

/* Because of limitations in the top_carry_adjust code, there is a limit */
/* on the size of k that can be handled.  This isn't a big deal since the */
/* zero-padded implementation should use the same FFT length.  Check to see */
/* if this k can be handled.  K must fit in the top three words for */
/* one-pass FFTs and within the top two words of two-pass FFTs. */

			if ((impl_jmptab->flags & 0x1FF) == 0 &&
			    logbk > ceil (jmptab->fftlen * b_per_input_word) - ceil ((jmptab->fftlen-3) * b_per_input_word)) {
// This assert is designed to find any cases where using 4 or more carry adjust words
// would use a shorter FFT than using a zero-padded FFT.  We found an example:
// 102233*299^239-1.  It would use length 256 FTT versus a length 320 zero-padded FFT.
//				ASSERTG (zpad_jmptab == NULL || jmptab->fftlen >= zpad_jmptab->fftlen);
				goto next2;
			}

			if ((impl_jmptab->flags & 0x1FF) != 0 &&
			    logbk > ceil (jmptab->fftlen * b_per_input_word) - ceil ((jmptab->fftlen-2) * b_per_input_word)) {
// This assert is designed to find any cases where using 3 or more carry adjust words
// would use a shorter FFT than using a zero-padded FFT.  One example: 501*500^100000-1
// It would use a 112K FFT instead of a 144K zero-padded FFT.
//				ASSERTG (zpad_jmptab == NULL || jmptab->fftlen >= zpad_jmptab->fftlen);
				goto next2;
			}

/* We've found an FFT length to use */

			break;
		}

/* Move past procedure entries and counts to next jmptable entry */

next2:		jmptab = NEXT_SET_OF_JMPTABS (jmptab);
	}

/* If the zero pad FFT length is less than the DWT FFT length OR we */
/* are QA'ing every FFT implementation, then use the zero pad FFT length. */

	if (zpad_jmptab != NULL && zpad_jmptab->max_exp &&
	    (jmptab->max_exp == 0 || zpad_jmptab->fftlen < jmptab->fftlen || gwdata->qa_pick_nth_fft)) {
		gwdata->ZERO_PADDED_FFT = TRUE;
		gwdata->ALL_COMPLEX_FFT = FALSE;
		jmptab = zpad_jmptab;
		best_impl_id = zpad_best_impl_id;
	}

/* If we found a DWT table entry then use it. */

	else if (jmptab->max_exp) {
		gwdata->ZERO_PADDED_FFT = FALSE;
		gwdata->ALL_COMPLEX_FFT = (c > 0);
	}

/* Error - neither method could handle this huge number */

	else
		return (GWERROR_TOO_LARGE);

/* See if the user requested a larger than normal FFT size */

	if (gwdata->larger_fftlen_count) {
		gwdata->larger_fftlen_count--;
		gwdata->minimum_fftlen = jmptab->fftlen + 1;
		goto again;
	}

/* We've found the right "jump" table entry, save the pointer and FFT length */

	gwdata->jmptab = jmptab;
	gwdata->FFTLEN = jmptab->fftlen;

/************************************************************************/
/* Decide which implementation of this FFT length is best for this CPU. */
/************************************************************************/

/* Loop through all the implementations for this FFT length until we find */
/* the one best suited to this CPU. */

	qa_nth_fft = gwdata->ZERO_PADDED_FFT ? 100 : 1000;
	desired_bif = calculate_bif (gwdata, gwdata->FFTLEN);
	prev_proc_ptrs[0] = NULL;
	prev_proc_ptrs[1] = NULL;
	prev_proc_ptrs[2] = NULL;
	prev_proc_ptrs[3] = NULL;
	prev_proc_ptrs[4] = NULL;
	for ( ; ; ) {
		int	arch;		/* (blend=0,p3=1,p4=2,core=3,fma3=4,k8=5,etc.) */
		int	best_impl_for;	/* (CORE2=0,I7=1,etc.  See BIF_ definitions in mult.asm) */
		int	fft_type;	/* (home-grown=0, radix-4=1, r4delay=2, r4dwpn=3) */

/* Handle an FFT implementation not found condition.  Should only happen */
/* if we're benchmarking or QA'ing and we've tested every implementation. */

		if (! (jmptab->flags & 0x80000000)) {

/* If we are QA'ing every FFT implementation and we did not find another */
/* zero-padded FFT implementation, then go find a non-zero-padded one. */

			if (gwdata->qa_pick_nth_fft && gwdata->qa_pick_nth_fft < 1000) {
				gwdata->qa_pick_nth_fft = 1000;
				goto again;
			}

/* Else return an error */

			return (GWERROR_TOO_LARGE);
		}

/* If this CPU will crash running this FFT then skip this entry. */
/* Our K8 and K10 optimized FFTs requires prefetchw (3DNow!) capability. */
/* If this is an Intel CPU, skip over these implementations. */

		arch = (jmptab->flags >> 17) & 0xF;
		if ((arch == ARCH_K8 || arch == ARCH_K10) && ! (gwdata->cpu_flags & CPU_3DNOW_PREFETCH))
			goto next3;

/* If this CPU will crash running this FFT then skip this entry. */
/* Our FMA3 FFTs require support of the Intel FMA3 instruction. */

		if (arch == ARCH_FMA3 && ! (gwdata->cpu_flags & CPU_FMA3))
			goto next3;

/* Handle benchmarking case that selects the nth FFT implementation */
/* without regard to any other consideration.  NOTE: Due to the extreme */
/* penalty a K8 pays for using the movaps instruction that the Core and P4 */
/* implementations use, we will not benchmark these on a K8.  Also, since */
/* FMA3 FFTs always outperform their non-FMA3 alternatives, we'll skip the */
/* non-FMA3 FFTs. */

		if (gwdata->bench_pick_nth_fft) {
			if (jmptab->proc_ptr == prev_proc_ptrs[0]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[1]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[2]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[3]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[4]) goto next3;
			if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_K8 &&
			    (jmptab->flags & 0x1FF) != 0 &&
			    (arch == ARCH_P4 || arch == ARCH_P4TP || arch == ARCH_CORE))
				goto next3;
			if (gwdata->cpu_flags & CPU_FMA3 && (arch == ARCH_P4 || arch == ARCH_P4TP || arch == ARCH_CORE))
				goto next3;
			gwdata->bench_pick_nth_fft--;
			if (gwdata->bench_pick_nth_fft) goto next3;
			break;
		}

/* Handle the QA case that tests every possible FFT implementation */
/* Remember the FFT returned so that we can return a different FFT to */
/* the QA code next time. */

		if (gwdata->qa_pick_nth_fft) {
			if (jmptab->proc_ptr == prev_proc_ptrs[0]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[1]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[2]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[3]) goto next3;
			if (jmptab->proc_ptr == prev_proc_ptrs[4]) goto next3;
			if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_K8 &&
			    (arch == ARCH_P4 || arch == ARCH_P4TP || arch == ARCH_CORE))
				goto next3;
			qa_nth_fft++;
			if (qa_nth_fft <= gwdata->qa_pick_nth_fft) goto next3;
			gwdata->qa_picked_nth_fft = qa_nth_fft;
			break;
		}

/* If this is a small SSE2 FFT or an x87 FFT, then there is only one implementation */

		if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))) {
			if (gwdata->cpu_flags & CPU_SSE2) {
				if (gwdata->FFTLEN < 7168) break;
			} else
				break;
		}

/* The Radix-4/8 DJB FFT with partial normalization saves a few multiplies by doing part */
/* of the normalization during the forward and inverse FFT.  Unfortunately, this optimization */
/* makes the SUM(INPUTS) != SUM(OUTPUTS) error check impossible.  If the user prefers */
/* SUM(INPUTS) != SUM(OUTPUTS) error checking, then skip r4dwpn FFT type (except all-complex */
/* FFTs which never could support SUM(INPUTS) != SUM(OUTPUTS) error checking. */
/* Ignore this preference for AVX-512/AVX FFTs as we do not support two-pass FFTs with this error check. */

		fft_type = (jmptab->flags >> 21) & 0xF;
		if (gwdata->sum_inputs_checking &&
		    ! gwdata->ALL_COMPLEX_FFT &&
		    ! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX)) &&
		    fft_type == FFT_TYPE_RADIX_4_DWPN)
			goto next3;

/* If we got a best_impl_id using the benchmark_database, then see if this jmptable matches the best_impl_id */
/* Unfortunately, this duplicates much of the flags parsing found later on in this routine. */

		if (best_impl_id != -1) {
			int	flags, fft_type, arch, clm, p2size, no_prefetch, in_place;
			flags = jmptab->flags;
			no_prefetch = (flags >> 27) & 0x00000001;
			in_place = (flags >> 26) & 0x00000001;
			fft_type = (flags >> 21) & 0x0000000F;
			arch = (flags >> 17) & 0x0000000F;
			clm = (flags >> 9) & 0x0000000F;
			if ((flags & 0x0000001FF) == 511)
				p2size = 48;
			else if ((flags & 0x0000001FF) == 510)
				p2size = 80;
			else
				p2size = (flags & 0x0000001FF) << 6;
			if (internal_implementation_ids_match (best_impl_id, gwdata->FFTLEN, fft_type, no_prefetch, in_place, p2size, arch, clm)) break;
		}

/* Otherwise use old school method where best FFT implementation is hardwired into mult.asm's jmptable. */
/* See if this is the best implementation for this CPU architecture */

		else {
			best_impl_for = (jmptab->flags >> 13) & 0xF;
			if (best_impl_for == desired_bif) break;
		}

/* Move onto the next FFT implementation */

next3:		prev_proc_ptrs[4] = prev_proc_ptrs[3];
		prev_proc_ptrs[3] = prev_proc_ptrs[2];
		prev_proc_ptrs[2] = prev_proc_ptrs[1];
		prev_proc_ptrs[1] = prev_proc_ptrs[0];
		prev_proc_ptrs[0] = jmptab->proc_ptr;
		jmptab = INC_JMPTAB (jmptab);
	}

/* Remember the information from the chosen FFT implementation */

	flags = jmptab->flags;
	gwdata->GWPROCPTRS[0] = jmptab->proc_ptr;
	if (jmptab->flags & 0x40000000) {
		struct gwasm_alt_jmptab *altjmptab = (struct gwasm_alt_jmptab *) jmptab;
		gwdata->mem_needed = altjmptab->mem_needed;
	} else
		gwdata->mem_needed = jmptab->mem_needed;

/* Break the flags word from the jmptable entry into its constituent parts */
/* The 32-bit flags word is as follows (copied from mult.asm): */
/*	80000000h		always on */
/*	40000000h		on if there are 2 proc ptrs (main FFT entry point, address of routine to do second FFT pass) */
/*	2 SHL 26		(no prefetching - not used by gwnum) */
/*	1 SHL 26		(in_place) */
/*	fft_type SHL 21		(hg=0, r4=1, r4delay=2) */
/*	arch SHL 17		(blend=0,p3=1,p4=2,core=3,fma3=4,k8=5,etc.) */
/*	best_impl_for SHL 13	(CORE2=0,I7=1,P4_1024=2,etc.) */
/*	clm SHL 9		(1,2,4,8) */
/*	pass2size_over_64	(many valid values) */

	gwdata->NO_PREFETCH_FFT = (flags >> 27) & 0x00000001;
	gwdata->IN_PLACE_FFT = (flags >> 26) & 0x00000001;
	gwdata->FFT_TYPE = (flags >> 21) & 0x0000000F;
	gwdata->ARCH = (flags >> 17) & 0x0000000F;
	if (gwdata->cpu_flags & CPU_AVX512F) gwdata->PASS1_CACHE_LINES = ((flags >> 9) & 0x0000000F) * 8;
	else if (gwdata->cpu_flags & CPU_AVX) gwdata->PASS1_CACHE_LINES = ((flags >> 9) & 0x0000000F) * 4;
	else gwdata->PASS1_CACHE_LINES = ((flags >> 9) & 0x0000000F) * 2;
	if ((flags & 0x0000001FF) == 511) gwdata->PASS2_SIZE = 48;
	else if ((flags & 0x0000001FF) == 510) gwdata->PASS2_SIZE = 80;
	else if ((flags & 0x0000001FF) == 509) gwdata->PASS2_SIZE = 32768;
	else gwdata->PASS2_SIZE = (flags & 0x0000001FF) << 6;
	if (gwdata->PASS2_SIZE) gwdata->PASS1_SIZE = gwdata->FFTLEN / gwdata->PASS2_SIZE; /* Real values in a pass1 section */
	if (gwdata->PASS1_SIZE == 2) gwdata->PASS1_SIZE = 0;	/* Don't treat AVX-512 one-pass wrapper as a true pass 1 */

/* Set more info so that addr_offset called from gwmap_to_estimated_size can work without a call to gwsetup. */

	gwdata->GW_ALIGNMENT = 4096;	/* Guess an alignment so gwsize can return a reasonable value for */
					/* large page support in gwsetup */

/* Calculate the scratch area size -- needed by gwmemused without calling gwsetup */

	if (gwdata->PASS1_SIZE && !gwdata->IN_PLACE_FFT) {
		if (gwdata->cpu_flags & CPU_AVX512F) {		// AVX-512 scratch area size
			int	pass1_size, num_clmblks, clmblkdst, gaps;
			// Mimic the clmblkdst calculations in zmult.mac.  
			pass1_size = gwdata->PASS1_SIZE;
			num_clmblks = pass1_size >> 4;
			clmblkdst = gwdata->PASS1_CACHE_LINES * 128;
			if (pass1_size == 192 || pass1_size == 640 || pass1_size == 768 || pass1_size == 896 || pass1_size == 960 ||
			    pass1_size == 1152 || pass1_size == 1280 || pass1_size == 1344 || pass1_size == 1536 || pass1_size == 1920 ||
			    pass1_size == 2304) { // Oddball pass 1 sizes
				gwdata->SCRATCH_SIZE = num_clmblks * (clmblkdst + 64) - 64;
			} else if (gwdata->PASS1_CACHE_LINES == 8) {	// clm=1
				gaps = num_clmblks / 8 - 1;
				gwdata->SCRATCH_SIZE = num_clmblks * clmblkdst + gaps * 128;
			} else if (gwdata->PASS1_CACHE_LINES == 16) {	// clm=2
				gaps = num_clmblks / 8 - 1;
				gwdata->SCRATCH_SIZE = num_clmblks * clmblkdst + gaps * 192;
			} else {					// clm=4 or more
				gaps = num_clmblks / 8;
				gwdata->SCRATCH_SIZE = num_clmblks * (clmblkdst + 64) + gaps * -64;
			}
		}
		else if (gwdata->cpu_flags & CPU_AVX) {		// AVX scratch area size
			int	pass1_size, pass1_chunks, gaps;
			// For small clms, AVX pads 64 bytes every 8 chunks (where pass1_chunks is the number of
			// cache lines that a pass 1 data set occupies).  For large clms, AVX pads 64 bytes every
			// chunk, and -64 bytes every 8 chunks.
			pass1_size = gwdata->PASS1_SIZE;
			pass1_chunks = pass1_size >> 3;
			if ((pass1_size == 384 && gwdata->ALL_COMPLEX_FFT) ||
			    (pass1_size == 640 && gwdata->ALL_COMPLEX_FFT) ||
			    (pass1_size == 1536 && gwdata->ALL_COMPLEX_FFT)) { // Oddball pass 1 sizes
				gwdata->SCRATCH_SIZE = pass1_chunks * (gwdata->PASS1_CACHE_LINES * 64 + 64) - 64;
			} else if (gwdata->PASS1_CACHE_LINES == 4) {	// clm=1
				gaps = pass1_chunks / 8 - 1;
				gwdata->SCRATCH_SIZE = pass1_chunks * gwdata->PASS1_CACHE_LINES * 64 + gaps * 64;
			} else if (gwdata->PASS1_CACHE_LINES == 8) {	// clm=2
				gaps = pass1_chunks / 8 - 1;
				gwdata->SCRATCH_SIZE = pass1_chunks * gwdata->PASS1_CACHE_LINES * 64 + gaps * 192;
			} else {					// clm=4 or more
				gaps = pass1_chunks / 8;
				gwdata->SCRATCH_SIZE = pass1_chunks * (gwdata->PASS1_CACHE_LINES * 64 + 64) + gaps * -64;
			}
		} else if (gwdata->cpu_flags & CPU_SSE2) {	// SSE2 scratch area size
			int	pass1_chunks, gaps;
			// SSE2 pads 128 bytes every 8 chunks
			pass1_chunks = gwdata->PASS1_SIZE >> 2;
			gaps = pass1_chunks / 8 - 1;
			gwdata->SCRATCH_SIZE = pass1_chunks * gwdata->PASS1_CACHE_LINES * 64 + gaps * 128;
		} else {					// x87 scratch area size
			int	pass1_chunks, gaps;
			// x87 pads 64 bytes every 32 chunks
			pass1_chunks = gwdata->PASS1_SIZE >> 1;
			gaps = pass1_chunks / 32 - 1;
			if (gaps < 3) gaps = 0;
			gwdata->SCRATCH_SIZE = pass1_chunks * gwdata->PASS1_CACHE_LINES * 32 + gaps * 64;
		}
	} else
		gwdata->SCRATCH_SIZE = 0;

/* Calculate the gaps between data blocks.  This is used by addr_offset called from */
/* gwmap_to_estimated_size without a call to gwsetup.  This must match the calculations */
/* done in the set_FFT_constants macro in zmult.mac, ymult.mac, and xmult.mac */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		gwdata->FOURKBGAPSIZE = 0;		/* BUG  No padding for now  */
		gwdata->PASS2GAPSIZE = 64;		/* BUG  64 bytes padding for now, might want zero for small pass 2 sizes */
		if (gwdata->PASS2_SIZE == 3072 || gwdata->PASS2_SIZE == 3584 || gwdata->PASS2_SIZE == 3840 || gwdata->PASS2_SIZE == 4096 ||
		    gwdata->PASS2_SIZE == 5120 || gwdata->PASS2_SIZE == 5376 || gwdata->PASS2_SIZE == 6144 || gwdata->PASS2_SIZE >= 7168) {
			gwdata->FOURKBGAPSIZE = 64;	/* Proven good choice for 7680*//* BUG  Try 128, others? */
			gwdata->PASS2GAPSIZE = -64;	/* Proven good choice for 7680,clm=2*/
							/* +64 was bad for 7680,clm=2 since 1951 cache lines mod 64 = 31  causing */
							/* 4KB overlaps when clm=2 reads 4 cache lines whereas 1949 mod 61 = 29 did not? */
							/* BUG  try 4kbgap = 128 to work better with hardware prefetcher's adjacent line try 64, others? */
		}
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		// Pass 2 sizes avoid the 4KB distance problem with 64, 128, or 192 pad bytes
		if (gwdata->PASS2_SIZE == 48 || gwdata->PASS2_SIZE == 64 || gwdata->PASS2_SIZE == 80 || gwdata->PASS2_SIZE == 192 || gwdata->PASS2_SIZE == 320)
			gwdata->FOURKBGAPSIZE = 0;
		if (gwdata->PASS2_SIZE == 256 || gwdata->PASS2_SIZE == 768 || gwdata->PASS2_SIZE == 2048 || gwdata->PASS2_SIZE == 2304 ||
		    gwdata->PASS2_SIZE == 3072 || gwdata->PASS2_SIZE == 3840 || gwdata->PASS2_SIZE == 4096 || gwdata->PASS2_SIZE == 5120 ||
		    gwdata->PASS2_SIZE == 6144 || gwdata->PASS2_SIZE == 6400 || gwdata->PASS2_SIZE == 8192 || gwdata->PASS2_SIZE == 9216 ||
		    gwdata->PASS2_SIZE == 10240 || gwdata->PASS2_SIZE == 12288 || gwdata->PASS2_SIZE == 15360 || gwdata->PASS2_SIZE == 16384 ||
		    gwdata->PASS2_SIZE == 20480 || gwdata->PASS2_SIZE == 25600)
			gwdata->FOURKBGAPSIZE = 64;
		if (gwdata->PASS2_SIZE == 1024 || gwdata->PASS2_SIZE == 1280 || gwdata->PASS2_SIZE == 1536 || gwdata->PASS2_SIZE == 2560 ||
		    gwdata->PASS2_SIZE == 4608 || gwdata->PASS2_SIZE == 7680 || gwdata->PASS2_SIZE == 12800)
			gwdata->FOURKBGAPSIZE = 128;
		// Make sure pass 2 block gapsize matches the computation of blkdst in ymult.mac
		if (gwdata->FOURKBGAPSIZE == 0)
			gwdata->PASS2GAPSIZE = 0;
		else if (gwdata->PASS2_SIZE == 2304 || gwdata->PASS2_SIZE == 9216)
			gwdata->PASS2GAPSIZE = 4*64;
		else if (gwdata->PASS2_SIZE * 2 * 8 / 4096 * gwdata->FOURKBGAPSIZE % 128 == 0)
			gwdata->PASS2GAPSIZE = -64;
		else
			gwdata->PASS2GAPSIZE = 0;
		// For one-pass FFTs, select the optimal padding frequency
		if (gwdata->PASS2_SIZE == 0) {
			if (gwdata->FFTLEN == 1536 || gwdata->FFTLEN == 2560)
				gwdata->FOURKBGAPSIZE = 16;		/* Pad every 1KB */
			else if (gwdata->FFTLEN == 1024 || gwdata->FFTLEN == 3072 || gwdata->FFTLEN == 5120)
				gwdata->FOURKBGAPSIZE = 32;		/* Pad every 2KB */
			else if (gwdata->FFTLEN == 2048 || gwdata->FFTLEN == 4096 || gwdata->FFTLEN == 6144 || gwdata->FFTLEN == 8192 ||
				 gwdata->FFTLEN == 10240 || gwdata->FFTLEN == 12288 || gwdata->FFTLEN == 16384 || gwdata->FFTLEN == 18432 ||
				 gwdata->FFTLEN == 20480 || gwdata->FFTLEN == 24576 || gwdata->FFTLEN == 32768)
				gwdata->FOURKBGAPSIZE = 64;		/* Pad every 4KB */
			else
				gwdata->FOURKBGAPSIZE = 0;		/* No padding */
		}
	} else {
		if (gwdata->PASS2_SIZE * 2 * 2 * 8 / 8192 * 128 % 256 == 0)
			gwdata->PASS2GAPSIZE = -128;
		else
			gwdata->PASS2GAPSIZE = 0;
	}

/* Copy or calculate various counts */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	if (asm_data != NULL) {
		if (gwdata->cpu_flags & CPU_AVX512F) {
			/* Copy the second pass routine ptr */
			struct gwasm_alt_jmptab *altjmptab = (struct gwasm_alt_jmptab *) jmptab;
			asm_data->u.zmm.ZMM_PASS2_ROUTINE = altjmptab->pass2_proc_ptr;

			if (gwdata->PASS1_SIZE == 0) {
				/* Count of sections for add, sub, addsub, zgw_carries_wpn in a one pass FFT */
				asm_data->addcount1 = gwdata->FFTLEN / 128;
			}
			else {
				/* Count of sections for add, sub, addsub, zgw_carries_wpn.  Value is pass 1 size over 16 doubles */
				asm_data->addcount1 = gwdata->PASS1_SIZE / 16;
			}
		}
		else if (gwdata->cpu_flags & CPU_AVX) {
			if (gwdata->PASS2_SIZE == 0) {
				if (gwdata->FOURKBGAPSIZE) {				/* Pad every 1KB, 2KB, or 4KB */
					asm_data->count1 = gwdata->FOURKBGAPSIZE;	/* Cache lines before a padding occurs */
					asm_data->normcount1 = gwdata->FFTLEN / 32 / gwdata->FOURKBGAPSIZE; /* Number of padding groups in each of the 4 sections */
					asm_data->count2 = asm_data->count1 / 2;	/* Counter for add/sub quick functions */
					asm_data->addcount1 = asm_data->normcount1 * 4;	/* Number of add/sub quick padding groups */
				} else {						/* No padding */
					asm_data->count1 = gwdata->FFTLEN / 32;		/* Cache lines in a section */
					asm_data->normcount1 = 0;			/* Number of padding groups in each of the 4 sections */
					asm_data->count2 = gwdata->FFTLEN / 16;		/* Counter for add/sub quick functions (4 AVX words processed at a time) */
					asm_data->addcount1 = 1;			/* Number of add/sub quick padding groups */
				}
			}
			else {
				/* Count of sections for add, sub, addsub, ygw_carries_wpn */
				asm_data->addcount1 = (gwdata->FFTLEN / 2) / (gwdata->PASS2_SIZE * 4);
				/* NOTE: more counts (count2, count3) for 2-pass FFTs are set in yr4dwpn_build_pass1_table */

				/* Copy the second pass routine ptr.  This is only present when first pass code is shared. */
				if (jmptab->flags & 0x40000000) {
					struct gwasm_alt_jmptab *altjmptab = (struct gwasm_alt_jmptab *) jmptab;
					asm_data->u.ymm.YMM_PASS2_ROUTINE = altjmptab->pass2_proc_ptr;
				}
			}
		} else if (gwdata->cpu_flags & CPU_SSE2) {
			if (gwdata->PASS2_SIZE == 0) {
				const struct gwasm_jmptab *last_jmptab;
				/* Copy 7 counts for one pass SSE2 FFTs.  The counts are after the last jmptab entry */
				last_jmptab = LAST_JMPTAB (jmptab);
				asm_data->addcount1 = last_jmptab->counts[0];
				asm_data->normcount1 = last_jmptab->counts[1];
				asm_data->count1 = last_jmptab->counts[2];
				asm_data->count2 = last_jmptab->counts[3];
				asm_data->count3 = last_jmptab->counts[4];
				asm_data->count4 = last_jmptab->counts[5];
				asm_data->count5 = last_jmptab->counts[6];
			} else if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
				int	pfa, pfa_shift;
				/* Count of pass 2 sections */
				asm_data->addcount1 = gwdata->PASS1_SIZE / 4;
				/* Count of all-complex pass 2 sections */
				asm_data->count1 = asm_data->addcount1;
				if (gwdata->ZERO_PADDED_FFT || c < 0) asm_data->count1--;
				/* Up to three section counts for gw_carries.  Examples: */
				/* 320 pass 2 sections: (256 << 11) + 64 */
				/* 384 pass 2 sections: 384 */
				/* 448 pass 2 sections: (256 << 22) + (128 << 11) + 64 */
				/* 512 pass 2 sections: 512 */
				for (pfa = asm_data->addcount1, pfa_shift = 0;
				     pfa > 8;
				     pfa >>= 1, pfa_shift++);
				if (pfa == 5)
					asm_data->count3 = ((4 << 11) + 1) << pfa_shift;
				else if (pfa == 7)
					asm_data->count3 = ((4 << 22) + (2 << 11) + 1) << pfa_shift;
				else
					asm_data->count3 = asm_data->addcount1;
			} else {
				/* Count of pass 2 sections */
				asm_data->addcount1 = gwdata->PASS1_SIZE / 4;
				/* Count of all-complex pass 2 sections */
				asm_data->count1 = asm_data->addcount1;
				if (gwdata->ZERO_PADDED_FFT || c < 0) asm_data->count1--;
				/* Only one section counts for gw_carries */
				asm_data->count3 = asm_data->addcount1;
			}
		} else {
			if (gwdata->PASS2_SIZE == 0) {
				/* 5 counts for one pass x87 FFTs */
				asm_data->count1 = jmptab->counts[0];
				asm_data->count2 = jmptab->counts[1];
				asm_data->count3 = jmptab->counts[2];
				asm_data->count4 = jmptab->counts[3];
				asm_data->count5 = jmptab->counts[4];
			} else {
				/* Count of pass 2 sections */
				asm_data->addcount1 = gwdata->PASS1_SIZE / 2;
				/* Count of all-complex pass 2 sections */
				asm_data->count1 = asm_data->addcount1;
				if (gwdata->ZERO_PADDED_FFT || c < 0) asm_data->count1--;
			}
		}
	}

/* All done */

	return (0);
}

/* Code to manage sharing sin/cos data where possible amongst several gwnum callers */

#define FIXED_PASS1_SINCOS_DATA			1
#define PASS2_REAL_SINCOS_DATA			2
#define PASS2_COMPLEX_SINCOS_DATA		3

struct shareable_sincos_data {
	struct shareable_sincos_data *next;	/* Next in linked list of shareable data blocks */
	double	*data;				/* Shareable data */
	size_t	data_size;			/* Size of the shareable data */
	int	use_count;			/* Count of times data is shared */
};
struct shareable_sincos_data *shareable_data = NULL;  /* Linked list of shareable data blocks */
gwmutex	shareable_lock;				/* This mutex limits one caller into sharing routines */
int	shareable_lock_initialized = FALSE;	/* Whether shareable mutex is initialized */

/* Share sin/cos data where possible amongst several gwnum callers */

double *share_sincos_data (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	int	table_type,	/* Type of the sin/cos table defined above */
	double *table,		/* Sin/cos data to share */
	size_t	table_size)	/* Size of the sin/cos data */
{
#ifdef SHARE_SINCOS_DATA
	struct shareable_sincos_data *p;	/* Ptr to a shareable data block */

/* No need to share zero-sized tables */

	if (table_size == 0) return (table);

/* Initialize the mutex if necessary, then grab the lock */

	if (!shareable_lock_initialized) {
		gwmutex_init (&shareable_lock);
		shareable_lock_initialized = TRUE;
	}
	gwmutex_lock (&shareable_lock);

/* Look through the list of shareable blocks looking for a match.  If we find a match, use it! */

	for (p = shareable_data; p != NULL; p = p->next) {
		if (p->data_size == table_size && !memcmp (p->data, table, table_size)) {
			p->use_count++;
			gwmutex_unlock (&shareable_lock);
			return (p->data);
		}
	}

/* Unfortunately, no match.  Copy the data so that a future gwnum caller can share the data */

	p = (struct shareable_sincos_data *) malloc (sizeof (struct shareable_sincos_data));
	if (p == NULL) {
		gwmutex_unlock (&shareable_lock);
		return (table);
	}
	p->data = (double *) aligned_malloc (table_size, 4096);
	if (p->data == NULL) {
		free (p);
		gwmutex_unlock (&shareable_lock);
		return (table);
	}
	memcpy (p->data, table, table_size);
	p->data_size = table_size;
	p->use_count = 1;
	p->next = shareable_data;
	shareable_data = p;
	gwmutex_unlock (&shareable_lock);
	return (p->data);
#else
	return (table);
#endif
}

/* Free shared sin/cos data */

void unshare_sincos_data (
	double *table)		/* Possibly shared sin/cos data */
{
#ifdef SHARE_SINCOS_DATA
	struct shareable_sincos_data *p;		/* Ptr to a shareable data block */
	struct shareable_sincos_data **ptr_to_p;	/* Linked list pointer to patch */

/* Ignore NULL table.  Should only happen when there are errors during gwsetup. */

	if (table == NULL) return;

/* Grab the lock */

	if (!shareable_lock_initialized) return;
	gwmutex_lock (&shareable_lock);

/* Look through the list of shareable blocks to find this shareable block */
/* Decrement the use count and when it reaches zero, free the memory */

	for (ptr_to_p = &shareable_data; (p = *ptr_to_p) != NULL; ptr_to_p = &p->next) {
		if (p->data != table) continue;
		if (--p->use_count == 0) {
			aligned_free (p->data);
			*ptr_to_p = p->next;
			free (p);
		}
		break;
	}

/* Free the lock and return */

	gwmutex_unlock (&shareable_lock);
#endif
}

/* Initialize gwhandle for a future gwsetup call. */
/* The gwinit function has been superceeded by gwinit2.  By passing in the */
/* version number we can verify the caller used the same gwnum.h file as the */
/* one he eventually links with.  The sizeof (gwhandle) structure is used */
/* to verify he compiles with the same structure alignment options that */
/* were used when compiling gwnum.c.  For compatibility with existing code */
/* we delay reporting any compatibility problems until gwsetup is called. */

void gwinit2 (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	int	struct_size,	/* Size of the gwdata structure */
	const char *version_string)
{

/* See if caller is using the same gwnum.h file that was used when */
/* this file was compiled.  Checking structure size also verifies he */
/* used the same compiler switches - especially regarding alignment. */
/* As a hack, we use GWERROR to record delayed error messages. */

	if (strcmp (version_string, GWNUM_VERSION)) {
		gwdata->GWERROR = GWERROR_VERSION_MISMATCH;
		return;
	}
	if (struct_size != sizeof (gwhandle)) {
		gwdata->GWERROR = GWERROR_STRUCT_SIZE_MISMATCH;
		return;
	}

/* Read the gwnum.txt INI file */

	gwbench_read_data ();

/* Initialize gwhandle structure with the default values */

	memset (gwdata, 0, sizeof (gwhandle));
	gwdata->safety_margin = 0.0;
	gwdata->maxmulbyconst = 3;
	gwdata->minimum_fftlen = 0;
	gwdata->larger_fftlen_count = 0;
	gwdata->num_threads = 1;
	gwdata->force_general_mod = 0;
	gwdata->use_irrational_general_mod = 0;
	gwdata->use_large_pages = 0;
	gwdata->use_benchmarks = 1;
	gwdata->mem_needed = GWINIT_WAS_CALLED_VALUE;	/* Special code checked by gwsetup to ensure gwinit was called */

/* Init structure that allows giants and gwnum code to share allocated memory */

	init_ghandle (&gwdata->gdata);
	gwdata->gdata.allocate = &gwgiantalloc;
	gwdata->gdata.free = &gwgiantfree;
	gwdata->gdata.deallocate = &gwgiantdealloc;
	gwdata->gdata.handle = (void *) gwdata;

/* If CPU type and speed have not been initialized by the caller, do so now. */

	if (CPU_FLAGS == 0 && CPU_SPEED == 0.0) {
		guessCpuType ();
		guessCpuSpeed ();
	}
	gwdata->cpu_flags = CPU_FLAGS;

/* We have not and will not write FMA3 or AVX-512 FFTs for 32-bit OSes */

#ifndef X86_64
	gwdata->cpu_flags &= ~(CPU_AVX512F | CPU_FMA3);
#endif

/* FMA3 FFTs require both AVX and FMA3 instructions.  This will always be the case when CPUID */
/* is queried.  However, prime95 has an option to turn off just the AVX bit with CpuSupportsAVX=0. */
/* Handle, this oddball case by also turning off FMA3. */

	if (! (gwdata->cpu_flags & CPU_AVX)) gwdata->cpu_flags &= ~CPU_FMA3;

/* AMD Bulldozer is faster using SSE2 rather than AVX. */
/* Why do we do this here when calculate_bif selects K10 FFTs???  Is it so that gwnum_map_to_timing and other */
/* informational routines return more accurate information?   Since the code below seems to work, leave it as is. */
/* 2019: A user reports that 3rd & 4th generation Bulldozer (Steamroller & Excavator) is better at AVX. */
/* I'm confident that 1st and 2nd are dismal (Bulldozer and Piledriver).  Consequently the code below was changed to */
/* check the extended model number using the chart at https://en.wikipedia.org/wiki/List_of_AMD_CPU_microarchitectures */

	if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_BULLDOZER && ((CPU_SIGNATURE >> 16) & 0xF) < 3)
		gwdata->cpu_flags &= ~(CPU_AVX512F | CPU_AVX | CPU_FMA3);
}

/* Allocate memory and initialize assembly code for arithmetic */
/* modulo k*b^n+c */

int gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. */
{
	int	gcd, error_code, setup_completed;
	double	orig_k;
	unsigned long orig_n;

/* Make sure gwinit was called */

	if (gwdata->mem_needed != GWINIT_WAS_CALLED_VALUE) return (GWERROR_NO_INIT);

/* Return delayed errors from gwinit2 */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Sanity check the k,b,n,c values */

	if (k < 1.0) return (GWERROR_K_TOO_SMALL);
	if (k > 9007199254740991.0) return (GWERROR_K_TOO_LARGE);
	if (gwdata->minimum_fftlen == 0) {
		if (gwdata->cpu_flags & CPU_AVX512F) {
			if (log2(b) * (double) n > MAX_PRIME_AVX512) return (GWERROR_TOO_LARGE);
		} else if (gwdata->cpu_flags & CPU_FMA3) {
			if (log2(b) * (double) n > MAX_PRIME_FMA3) return (GWERROR_TOO_LARGE);
		} else if (gwdata->cpu_flags & CPU_AVX) {
			if (log2(b) * (double) n > MAX_PRIME_AVX) return (GWERROR_TOO_LARGE);
		} else if (gwdata->cpu_flags & CPU_SSE2) {
			if (log2(b) * (double) n > MAX_PRIME_SSE2) return (GWERROR_TOO_LARGE);
		} else {
			if (log2(b) * (double) n > MAX_PRIME) return (GWERROR_TOO_LARGE);
		}
	}
	if ((k == 1.0 && n == 0 && c == 0) ||
	    (c < 0 && n * log ((double) b) + log (k) <= log ((double) 1-c)))
		return (GWERROR_TOO_SMALL);

/* Init */

	setup_completed = FALSE;
	orig_k = k;
	orig_n = n;

/* Our code fails if k is a power of b.  For example, 3481*59^805-1 which */
/* equals 59^807-1.  I think this is because gwfft_base(FFTLEN) is off by one */
/* because even quad-precision floats won't calculate FFTLEN * num_b_per_word */
/* correctly.  There is an easy fix, if k is divisible by b we divide k by b */
/* and add one to n. */

	while (k > 1.0 && b > 1 && fmod (k, (double) b) == 0.0) {
		k = k / (double) b;
		n = n + 1;
	}

/* Our code fast code fails if k and c are not relatively prime.  This */
/* is because we cannot calculate 1/k.  Although the user shouldn't call */
/* us with this case, we handle it anyway by reverting to the slow general */
/* purpose multiply routines. */

	if (c == 0)
		gcd = 0;
	else if (k == 1.0 || labs (c) == 1)
		gcd = 1;
	else {
		stackgiant(kg,2);
		stackgiant(cg,2);
		dbltog (k, kg);
		itog (labs (c), cg);
		gcdg (kg, cg);
		gcd = cg->n[0];
	}

/* Call the internal setup routine when we can.  Gcd (k, c) must be 1, */
/* k * mulbyconst and c * mulbyconst cannot be too large.  Also, the FFT */
/* code has bugs when there are too few bits per FFT.  Rather than make */
/* difficult fixes we simply force these small numbers to use the generic */
/* reduction.  In truth, the caller should use a different math package for */
/* these small numbers. */

	if (gcd == 1 &&
	    k * gwdata->maxmulbyconst <= MAX_ZEROPAD_K &&
	    labs (c) * gwdata->maxmulbyconst <= MAX_ZEROPAD_C &&
	    log2(b) * (double) n >= 350.0 &&
	    (b == 2 || (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2))) &&
	    !gwdata->force_general_mod) {
		error_code = internal_gwsetup (gwdata, k, b, n, c);
		if (error_code == 0) setup_completed = TRUE;
		else if (b == 2) return (error_code);
		gwdata->GENERAL_MOD = FALSE;
	}

/* Emulate k not relatively prime to c, small n values, and */
/* large k or c values with a call to the general purpose modulo setup code. */

	if (!setup_completed) {
		/* If we've already copied the modulus, use it.  For example, */
		/* gwsetup_general_mod_giant on (2^313+1)/3 will call this routine */
		/* to try an IBDWT on 2^313+1.  This number is too small and */
		/* we need to revert back to a general mod on (2^313+1)/3. */
		if (gwdata->GW_MODULUS != NULL) {
			gwdata->force_general_mod = TRUE;
			error_code = gwsetup_general_mod_giant (gwdata, gwdata->GW_MODULUS);
			if (error_code) return (error_code);
		} else {
			double	bits;
			giant	g;
			bits = (double) n * log2 (b);
			g = allocgiant (((unsigned long) bits >> 5) + 4);
			if (g == NULL) return (GWERROR_MALLOC);
			ultog (b, g);
			power (g, n);
			dblmulg (k, g);
			iaddg (c, g);
			gwdata->force_general_mod = TRUE;
			error_code = gwsetup_general_mod_giant (gwdata, g);
			free (g);
			if (error_code) return (error_code);
		}
	}

/* For future messages, format the input number as a string */

	gw_as_string (gwdata->GWSTRING_REP, orig_k, b, orig_n, c);

/* Return success */

	return (0);
}

/* These setup routines are for operations modulo an arbitrary binary number. */
/* This is three times slower than the special forms above. */

int gwsetup_general_mod (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	const uint32_t *array,	/* The modulus as an array of 32-bit values */
	uint32_t arraylen)	/* Number of values in the array */
{
	giantstruct tmp;
	tmp.sign = arraylen;
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	return (gwsetup_general_mod_giant (gwdata, &tmp));
}

int gwsetup_general_mod_64 (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	const uint64_t *array,	/* The modulus as an array of 64-bit values */
	uint64_t arraylen)	/* Number of values in the array */
{
	giantstruct tmp;
	tmp.sign = (int) arraylen * 2;
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	return (gwsetup_general_mod_giant (gwdata, &tmp));
}

/* Setup the FFT code for generic reduction */

int gwsetup_general_mod_giant (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	giant	g)		/* The modulus */
{
	unsigned long bits;	/* Bit length of modulus */
	int	convertible;	/* Value can be converted to (k*2^n+c)/d */
	double	k;
	unsigned long n;
	signed long c;
	unsigned long d;
	unsigned long safety_bits;
	const struct gwasm_jmptab *info;
	int	error_code;
	unsigned long fftlen, max_exponent, desired_n;
	giant	modified_modulus, tmp;

/* Make sure gwinit was called */

	if (gwdata->mem_needed != GWINIT_WAS_CALLED_VALUE) return (GWERROR_NO_INIT);

/* Return delayed errors from gwinit2 */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Init */

	bits = bitlen (g);

/* Examine the giant to see if it a (k*2^n+c)/d value that we can better optimize. */
/* Also detect nasty bit patterns, like Phi (82730,2), where multiplying by a small d */
/* results in a less nasty bit pattern for the modulus. */

	d = 1;
	convertible = (!gwdata->force_general_mod &&
		       bits > 300 &&
		       convert_giant_to_k2ncd (g, &k, &n, &c, &d));

/* Copy the modulus except for convertible k*2^n+c values */

	if (!convertible || d != 1) {
		/* Yuk, we may already have saved the modulus.  For example, (2^313+1)/3 */
		/* will come through here and save the modulus.  But gwsetup of 2^313+1 */
		/* is too small for IBDWT, so this routine is recursively called. */
		/* We cannot reallocate because that will cause a memory leak. */
		if (gwdata->GW_MODULUS == NULL) {
			gwdata->GW_MODULUS = allocgiant ((bits >> 5) + 1);
			if (gwdata->GW_MODULUS == NULL) {
				gwdone (gwdata);
				return (GWERROR_MALLOC);
			}
			gtog (g, gwdata->GW_MODULUS);
		}
	}

/* Setup for values we are converting to use faster k*2^n+c FFTs. */

	if (convertible) {
		error_code = gwsetup (gwdata, k, 2, n, c);
		if (error_code) return (error_code);
		if (d != 1) {
			char	buf[60];
			strcpy (buf, gwdata->GWSTRING_REP);
			if (isdigit (buf[0]))
				sprintf (gwdata->GWSTRING_REP, "(%s)/%lu", buf, d);
			else
				sprintf (gwdata->GWSTRING_REP, "%s/%lu", buf, d);
		}
		return (0);
	}

/* If we need to multiply the modulus by a small d value, do so here */

	if (d != 1) {
		modified_modulus = allocgiant ((bits >> 5) + 2);
		if (modified_modulus == NULL) {
			gwdone (gwdata);
			return (GWERROR_MALLOC);
		}
		gtog (g, modified_modulus);
		ulmulg (d, modified_modulus);
		g = modified_modulus;
		bits = bitlen (g);
	} else
		modified_modulus = NULL;

/* We will need twice the number of input bits plus some padding */

	n = bits + bits + 128;

/* Setup the FFT code in much the same way that gwsetup_without_mod does. */
/* Unless the user insists, we try for an integral number of bits per word. */
/* There are pathological bit patterns that generate huge roundoff errors. */
/* For example, if we test (10^828809-1)/9 and put exactly 18 bits into */
/* each FFT word, then every FFT word in GW_MODULUS_FFT will contain the */
/* same value!  Not exactly, the random data our FFTs require for small */
/* roundoff errors.  Thus, the caller may need to insist we use an */
/* irrational FFT on occasion. */

/* Call gwinfo and have it figure out the FFT length we will use. */
/* Since we zero the upper half of FFT input data, the FFT */
/* outputs will be smaller.  This lets us get about another 0.3 bits */
/* per input word. */

	gwdata->safety_margin -= 0.3;
	error_code = gwinfo (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	if (error_code) return (error_code);
	info = gwdata->jmptab;
	fftlen = info->fftlen;
	max_exponent = adjusted_max_exponent (gwdata, info);

/* Our FFTs don't handle cases where there are few bits per word because */
/* carries must be propagated over too many words.  Arbitrarily insist */
/* that n is at least 12 * fftlen.  */

	if (n < 12 * fftlen) n = 12 * fftlen;

/* Let the user request rational FFTs as they are a few percent faster */

	if (!gwdata->use_irrational_general_mod) {

/* If possible, increase n to the next multiple of FFT length.  This is */
/* because rational FFTs are faster than irrational FFTs (no FFT weights). */

		desired_n = round_up_to_multiple_of (n, fftlen);
		if (desired_n < max_exponent) n = desired_n;
	}

/* If the user requested irrational FFTs, then make sure the bits */
/* per FFT word will distribute the big and little words of the modulus */
/* semi-randomly.  For example, in the (10^828809-1)/9 case above, if */
/* bits-per-word is 18.5 or 18.25 you will still get non-random patterns */
/* in the FFT words. */

	else {
		double	prime_number, bits_per_word;

/* Round bits_per_word up to the next half-multiple of 1/prime_number */

		prime_number = 53.0;
		bits_per_word = (double) n / (double) fftlen;
		bits_per_word = (ceil (bits_per_word * prime_number) + 0.5)/ prime_number;

/* If possible, use the n associated with the just-computed bits-per-word */

		desired_n = (unsigned long) ceil (bits_per_word * (double) fftlen);
		if (desired_n < max_exponent) n = desired_n;
	}

/* If possible, increase n to the next multiple of FFT length. */
/* The extra bits allow gwsmallmul to avoid emulate_mod calls more often. */
/* We hope the 0.3 safety_limit increase above will avoid getting too */
/* close to the FFT limit as many users of this library turn on error */
/* checking (slower) when near the FFT limit. */
/* If that doesn't work, try adding a half FFT length instead. */

	if (n + fftlen < max_exponent)
		n = n + fftlen;
	else if (gwdata->use_irrational_general_mod && n + fftlen / 2 < max_exponent)
		n = n + fftlen / 2;

/* Now setup the assembly code */

	gwdata->safety_margin -= 0.3;
	error_code = internal_gwsetup (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	if (error_code) return (error_code);

// BUG - setting the bit_length to the modulus size will break gwtogiant.
// we need a better/more-consistent way of dealing with the various 
// needed bit_lengths.  Also, PFGW should not be reading the bit_length
// value in integer.cpp.
//	gwdata->bit_length = bits;
	
/* Allocate memory for an FFTed copy of the modulus. */

	gwdata->GW_MODULUS_FFT = gwalloc (gwdata);
	if (gwdata->GW_MODULUS_FFT == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	gianttogw (gwdata, g, gwdata->GW_MODULUS_FFT);
	gwfft (gwdata, gwdata->GW_MODULUS_FFT, gwdata->GW_MODULUS_FFT);

/* Calculate number of words to zero during the copy prior to calculating the quotient. */

	gwdata->GW_ZEROWORDSLOW = (unsigned long) floor ((double) bits / gwdata->avg_num_b_per_word);

/* A quick emulate_mod refresher: */
/* 1) The maximum size of a quotient is gwdata->n/2 bits due to the */
/*    zeroing of high words in the normalization routine.  Obviously the */
/*    reciprocal needs to be accurate to at least gwdata->n/2 bits. */
/* 2) So that the quotient doesn't overflow, the maximum size of a value */
/*    entering emulate_mod is gwdata->n/2+bits bits */
/* 3) So that gwsquare and gwmul don't send emulate_mod a value that is */
/*    too large, the maximum input to these routines should be (allowing */
/*    for an 8 bit mulbyconstant) is (gwdata->n/2+bits-8)/2 bits.  This is */
/*    used by gwsmallmul to know when an emulate_mod is required */
/* 4) We cannot quite push to the limits calculated above because we have to */
/*    make sure the quotient calculation does not produce more than gwdata->n */
/*    bits of result -- otherwise the *high* order bits of the quotient will be */
/*    corrupted.  We allow ourselves a small safety margin by decreasing the */
/*    number of reciprocal bits calculated.  The safety margin must be larger */
/*    than the number of unzeroed bits caused by using the floor function in */
/*    calculating GW_ZEROWORDSLOW. */

	safety_bits = bits - (unsigned long) ((double) gwdata->GW_ZEROWORDSLOW * gwdata->avg_num_b_per_word) + 3;

/* Precompute the reciprocal of the modified modulus. */

	gwdata->GW_RECIP_FFT = gwalloc (gwdata);
	if (gwdata->GW_RECIP_FFT == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	tmp = allocgiant ((gwdata->n >> 5) + 1);
	if (tmp == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	itog (1, tmp);
	gshiftleft (bits + gwdata->n / 2 - safety_bits, tmp);
	divg (g, tmp);		/* computes gwdata->n/2-safety_margin+1 bits of reciprocal */
	gshiftleft (gwdata->n - (bits + gwdata->n / 2 - safety_bits), tmp);
				/* shift so gwmul routines wrap */
				/* quotient to lower end of fft */
	gianttogw (gwdata, tmp, gwdata->GW_RECIP_FFT);
	gwfft (gwdata, gwdata->GW_RECIP_FFT, gwdata->GW_RECIP_FFT);
	free (tmp);
	free (modified_modulus);

/* Calculate the maximum allowable size of a number used as input */
/* to gwmul.  We will make sure gwsmallmul does not generate any */
/* results bigger than this. */

	gwdata->GW_GEN_MOD_MAX = (unsigned long)
		 floor ((double)((gwdata->n/2-safety_bits+bits-8)/2) / gwdata->avg_num_b_per_word);
	gwdata->GW_GEN_MOD_MAX_OFFSET = addr_offset (gwdata, gwdata->GW_GEN_MOD_MAX-1);

/* Set flag indicating general-purpose modulo operations are in force */

	gwdata->GENERAL_MOD = TRUE;

/* Reciprocals in generic modular reduction may well have a nasty bit pattern. */
/* Two test cases that brought this to light are 3*2^77574+3, and 3*8^86103+1 */
/* (when forced to use generic reduction).  These nasty patterns can trigger */
/* spurious SUM(INPUTS) != SUM(OUTPUTS) errors.  To counter this we increase */
/* the MAXDIFF setting. */
	
	gwdata->MAXDIFF *= 1000.0;

/* Create dummy string representation. Calling gtoc to get the first */
/* several digits would be better, but it is too slow. */

	sprintf (gwdata->GWSTRING_REP, "A %ld-bit number", bits);

/* Return success */

	return (0);
}

/* This setup routine is for operations without a modulo. In essence, */
/* you are using gwnums as a general-purpose FFT multiply library. */

int gwsetup_without_mod (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	unsigned long n)	/* Maximum number of bits in OUTPUT numbers. */
{
	const struct gwasm_jmptab *info;
	unsigned long fftlen, max_exponent, desired_n;
	int	error_code;

/* Make sure gwinit was called */

	if (gwdata->mem_needed != GWINIT_WAS_CALLED_VALUE) return (GWERROR_NO_INIT);

/* Return delayed errors from gwinit2 */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Call gwinfo and have it figure out the FFT length we will use. */
/* Since the user must zero the upper half of FFT input data, the FFT */
/* outputs will be smaller.  This lets us get about another 0.3 bits */
/* per input word. */

	gwdata->safety_margin -= 0.3;
	error_code = gwinfo (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	info = gwdata->jmptab;
	if (error_code) return (error_code);

	max_exponent = adjusted_max_exponent (gwdata, info);
	fftlen = info->fftlen;

/* If possible, increase n to the next multiple of FFT length.  This is */
/* because rational FFTs are faster than irrational FFTs (no FFT weights). */

	desired_n = round_up_to_multiple_of (n, fftlen);
	if (desired_n < max_exponent) n = desired_n;

/* Our FFTs don't handle cases where there are few bits per word because */
/* carries must be propagated over too many words.  Arbitrarily insist */
/* that n is at least 12 * fftlen.  */

	if (n < 12 * fftlen) n = 12 * fftlen;

/* Now setup the assembly code */

	gwdata->safety_margin -= 0.3;
	error_code = internal_gwsetup (gwdata, 1.0, 2, n, -1);
	gwdata->safety_margin += 0.3;
	if (error_code) return (error_code);

/* Set flag indicating general-purpose modulo operations are not in force */

	gwdata->GENERAL_MOD = FALSE;

/* Create dummy string representation. */

	strcpy (gwdata->GWSTRING_REP, "No modulus");

/* Return success */

	return (0);
}


/* Examine a giant to see if it a (k*2^n+c)/d value. */
/* Returns TRUE if conversion was successful. */

int convert_giant_to_k2ncd (
	giant	g,		/* Giant to examine */
	double	*k,		/* K in (K*2^N+C)/D. */
	unsigned long *n,	/* N in (K*2^N+C)/D. */
	signed long *c,		/* C in (K*2^N+C)/D. */
	unsigned long *d)	/* D in (K*2^N+C)/D. */
{
	unsigned long less_nasty_d;
	int	i;
	uint32_t quick_test;
	giant	test_g, alloc_g;
	stackgiant(tmp,3);

/* Loop through a lot of small d values in hopes of finding a multiplier that */
/* will convert the input value into a k*2^n+c value.  We do this because */
/* small d values are the ones that generate repeating bit patterns in */
/* the modulus and unexpectedly large round off errors during operations */

	alloc_g = NULL;
	less_nasty_d = 1;
	for (*d = 1; *d <= 999; *d += 2) {

/* Do a quick test to see if this is a viable candidate */

		quick_test = (g->n[1] * *d) & 0xFFFFFC00;
		if (quick_test != 0 && quick_test != 0xFFFFFC00) continue;

/* Compute g * d to see if it has the proper k*2^n+c bit pattern */

		if (*d == 1) {
			test_g = g;
		} else {
			if (alloc_g == NULL) {
				alloc_g = allocgiant (((bitlen (g) + 10) >> 5) + 1);
				if (alloc_g == NULL) return (FALSE);
			}
			ultog (*d, alloc_g);
			mulg (g, alloc_g);
			test_g = alloc_g;
		}

/* See if this d value might result in a less nasty bit pattern for */
/* emulate_mod.  For example, Phi(82730,2) behaves much better if you */
/* multiply the modulus by 11.  We'll assume that d produces a better */
/* modulus candidate if more than a quarter of the words are zero or */
/* minus one -- at least until someone improves on this scheme. */

		if (*d >= 3 && less_nasty_d == 1) {
			int	count = 0;
			for (i = 0; i < test_g->sign; i++)
				if (test_g->n[i] == 0 || test_g->n[i] == -1) count++;
			if (count >= (test_g->sign >> 2)) less_nasty_d = *d;
		}

/* See if low order 2 words are viable for a k*2^n+c candidate */

		*c = (int32_t) test_g->n[0];

		if (test_g->n[1] == 0 && test_g->n[0] <= MAX_ZEROPAD_C);
		else if (test_g->n[1] == 0xFFFFFFFF && *c < 0 && *c >= -MAX_ZEROPAD_C);
		else continue;

/* Examine the middle words */
	
		for (i = 2; i < test_g->sign - 1; i++)
			if (test_g->n[i] != test_g->n[1]) break;

/* Now see if the high bits can fit in a 51-bit k */

		tmp->n[0] = test_g->n[i]; tmp->sign = 1;
		if (test_g->sign - i >= 2) { tmp->n[1] = test_g->n[i+1]; tmp->sign = 2; }
		if (test_g->sign - i >= 3) { tmp->n[2] = test_g->n[i+2]; tmp->sign = 3; }
		if (test_g->sign - i >= 4) continue;
		if (test_g->n[1] == 0xFFFFFFFF) iaddg (1, tmp);

		*n = i * 32;
		while ((tmp->n[0] & 0x1) == 0) {
			gshiftright(1, tmp);
			(*n)++;
		}
		if (bitlen (tmp) > 51) continue;

/* Set k and return success */

		*k = tmp->n[0];
		if (tmp->sign == 2) *k += (double) tmp->n[1] * 4294967296.0;
		free (alloc_g);
		return (TRUE);
	}

/* No luck in finding a (k*2^n+c)/d equivalent for the input value */

	*d = less_nasty_d;
	free (alloc_g);
	return (FALSE);
}

/* Common setup routine for the three different user-visible setup routines */
/* Allocate memory and initialize assembly code for arithmetic */
/* modulo k*b^n+c */

int internal_gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	const struct gwasm_jmptab *info;
	void	*asm_data_alloc;
	struct gwasm_data *asm_data;
	int	error_code;
	unsigned long mem_needed;
	double	*tables;		/* Pointer tables we are building */
	unsigned long pass1_size;
	double	small_word, big_word, temp, asm_values[50];

/* Remember the arguments */

	gwdata->k = k;
	gwdata->b = b;
	gwdata->n = n;
	gwdata->c = c;

/* Init the FPU to assure we are in 64-bit precision mode */

	fpu_init ();

/* Allocate space for the assembly code global data.  This area is preceded */
/* by a temporary stack.  This allows the assembly code to access the global */
/* data using offsets from the stack pointer. */

	asm_data_alloc = aligned_malloc (sizeof (struct gwasm_data) + NEW_STACK_SIZE, 4096);
	if (asm_data_alloc == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	gwdata->asm_data = (char *) asm_data_alloc + NEW_STACK_SIZE;
	asm_data = (struct gwasm_data *) gwdata->asm_data;
	memset (asm_data, 0, sizeof (struct gwasm_data));

/* Select the proper FFT size for this k,b,n,c combination */

	error_code = gwinfo (gwdata, k, b, n, c);
	if (error_code) {
		gwdone (gwdata);
		return (error_code);
	}
	info = gwdata->jmptab;

/* Get pointer to fft info and allocate needed memory.  If we are */
/* trying to allocate large pages, then also allocate space for the */
/* one gwnum value that will be stored in large pages.  We allocate */
/* a little extra space to align the gwnum on a cache line boundary */
/* and because gwnum_size may not produce an accurate value prior */
/* to gwsetup completing. */

	tables = NULL;
	gwdata->large_pages_ptr = NULL;
	gwdata->large_pages_gwnum = NULL;
	mem_needed = gwdata->mem_needed + gwdata->SCRATCH_SIZE;
	if (gwdata->use_large_pages) {
		tables = (double *) large_pages_malloc (mem_needed + gwnum_size (gwdata) + 4096);
		if (tables != NULL) {
			/* Save large pages pointer for later freeing */
			gwdata->large_pages_ptr = tables;
			tables = align_ptr (tables, 128);
			/* Save pointer to the gwnum we also allocated, so */
			/* that first gwalloc call can return it. */
			gwdata->large_pages_gwnum = align_ptr ((char *) tables + mem_needed, 128);
		}
	}
	if (tables == NULL) {
		tables = (double *) aligned_malloc (mem_needed, 4096);
		if (tables == NULL) return (GWERROR_MALLOC);
	}
	gwdata->gwnum_memory = tables;

/* Do a seemingly pointless memset! */
/* The memset will walk through the allocated memory sequentially, which */
/* increases the likelihood that contiguous virtual memory will map to */
/* contiguous physical memory. */

	memset (tables, 0, mem_needed);

/* Copy values for asm code to use */

	asm_data->FFTLEN = gwdata->FFTLEN;
	asm_data->ZERO_PADDED_FFT = gwdata->ZERO_PADDED_FFT;
	asm_data->ALL_COMPLEX_FFT = gwdata->ALL_COMPLEX_FFT;
	asm_data->B_IS_2 = (b == 2);

/* Initialize the extended precision code that computes the FFT weights */

	gwdata->dd_data = gwdbldbl_data_alloc ();
	if (gwdata->dd_data == NULL) {
		gwdone (gwdata);
		return (GWERROR_MALLOC);
	}
	gwfft_weight_setup (gwdata->dd_data, gwdata->ZERO_PADDED_FFT, k, b, n, c, gwdata->FFTLEN);

/* Calculate the number of bits in k*b^n.  This will be helpful in */
/* determining how much memory to allocate for giants. */

	gwdata->bit_length = log2 (k) + n * log2 (b);

/* Calculate the average number of base b's stored in each FFT word.  The total */
/* number of base b's the underlying FFT works with (i.e. the point at which data */
/* wraps around to the low FFT word) is 2*n for a zero pad FFT and logb(k) + n */
/* otherwise. */	

	gwdata->avg_num_b_per_word = (gwdata->ZERO_PADDED_FFT ? n * 2.0 : (logb (k) + n)) / gwdata->FFTLEN;

/* Calculate the number of base b's stored in each small FFT word. */

	gwdata->NUM_B_PER_SMALL_WORD = (unsigned long) gwdata->avg_num_b_per_word;
	small_word = pow ((double) b, gwdata->NUM_B_PER_SMALL_WORD);
	big_word = (double) b * small_word;

/* Set a flag if this is a rational FFT.  That is, an FFT where all the weighting factors are 1.0. */
/* This happens when zero-padded or abs(c) is 1, and every FFT word has the same number of b's. */
/* The assembly code can make some obvious optimizations when all the FFT weights are one. */

	gwdata->RATIONAL_FFT = asm_data->RATIONAL_FFT =
		((double) gwdata->NUM_B_PER_SMALL_WORD == gwdata->avg_num_b_per_word) && (gwdata->ZERO_PADDED_FFT || labs (c) == 1);

/* Remember the maximum number of bits per word that this FFT length supports.  We use this in gwnear_fft_limit. */
/* Note that zero padded FFTs can support an extra 0.3 bits per word because of the all the zeroes. */

	gwdata->fft_max_bits_per_word = (double) adjusted_max_exponent (gwdata, info) / (double) gwdata->FFTLEN;
	if (gwdata->ZERO_PADDED_FFT) gwdata->fft_max_bits_per_word += 0.3;

/* Compute extra bits - the number of adds we can tolerate without */
/* a normalization operation. Under normal circumstances, max_bits */
/* will be greater than virtual bits, but playing with the safety margin */
/* or forcing use of a specific FFT length could change that. */

	gwdata->EXTRA_BITS = (unsigned long) pow (2.0, (gwdata->fft_max_bits_per_word - virtual_bits_per_word (gwdata)) / 2.0);

/* Determine the pass 1 size.  This affects how we build many of the sin/cos tables. */

	if (gwdata->PASS2_SIZE)
		pass1_size = gwdata->PASS1_SIZE;
	else
		pass1_size = gwdata->FFTLEN;

/* Initialize tables for r4dwpn AVX-512 FFT code */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		double	rndval_base;
		int	cnt;

/* Compute the normalization constants.  The rounding constant, a value larger than 3*2^51, is chosen such that */
/* it is a multiple of big_word and RNDVAL * big_word - RNDVAL fits in 53 bits. */
/* Let's figure out how to compute RNDVAL, by imagining a CPU with 20 bit mantissas instead of 53-bits */
/* First example is bigbase = 2^8, RNDVAL is 3*2^18: */
/*	RNDVAL = 1100 0000 0000 0000 0000 */
/*	RNDVAL * bigbase - RNDVAL = 1100 0000 0000 0000 0000 0000 0000 */
/*					    - 1100 0000 0000 0000 0000, which easily fits in 20 bits. */
/* Next example is a with 6-bit odd bigbase.  We choose RNDVAL such that it is a multiple of bigbase and a multiple of 2^6. */
/*	RNDVAL = 1100 0000 xxxx xx00 0000 */
/*	RNDVAL * bigbase - RNDVAL = yy yyyy yyyy yyyy yyyy yy00 0000 */
/*					  - 1100 0000 xxxx xx00 0000, which fits in 20 bits. */
/* Next example is a with 6-bit even bigbase.  We choose RNDVAL such that it is a multiple of bigbase and a multiple of 2^5. */
/*	RNDVAL = 1100 0000 0xxx xxx0 0000 */
/*	RNDVAL * bigbase - RNDVAL = yy yyyy yyyy yyyy yyyy yy00 0000 */
/*					  - 1100 0000 0xxx xxx0 0000, */
/* which does NOT fit in 20 bits, we must choose RNDVAL as a multiple of 2^6. */
/* This leaves us with the following simple algorithm: */
/*	cnt = count of bits in bigbase */
/*	RNDVAL = 3*2^51 / 2^cnt */
/*	RNDVAL += (bigbase - RNDVAL % bigbase) % bigbase */
/*	RNDVAL *= 2^cnt */

		asm_data->u.zmm.ZMM_LARGE_BASE = big_word;				/* Upper limit */
		asm_data->u.zmm.ZMM_LARGE_BASE_INVERSE = 1.0 / big_word;		/* Upper limit inverse */
		asm_data->u.zmm.ZMM_SMALL_BASE = small_word;				/* Lower limit */
		asm_data->u.zmm.ZMM_SMALL_BASE_INVERSE = 1.0 / small_word;		/* Lower limit inverse */

		// Use small_word for FFTs where there are no big words (and big_word might exceed fatal 2^26)
		// This is more than just RATIONAL_FFTs (e.g. 117^3072-5).
		rndval_base = ((double) gwdata->NUM_B_PER_SMALL_WORD == gwdata->avg_num_b_per_word ? small_word : big_word);
		ASSERTG (rndval_base < 67108864.0);
		cnt = (int) log2 (rndval_base) + 1;
		asm_data->u.zmm.ZMM_RNDVAL = 3.0 * 131072.0 * 131072.0 * 131072.0;	/* Rounding value */
		asm_data->u.zmm.ZMM_RNDVAL /= pow (2.0, cnt);
		asm_data->u.zmm.ZMM_RNDVAL += fltmod (rndval_base - fltmod (asm_data->u.zmm.ZMM_RNDVAL, rndval_base), rndval_base);
		asm_data->u.zmm.ZMM_RNDVAL *= pow (2.0, cnt);

		asm_data->u.zmm.ZMM_RNDVAL_TIMES_LARGE_BASE = asm_data->u.zmm.ZMM_RNDVAL * big_word - asm_data->u.zmm.ZMM_RNDVAL;
		asm_data->u.zmm.ZMM_RNDVAL_TIMES_SMALL_BASE = asm_data->u.zmm.ZMM_RNDVAL * small_word - asm_data->u.zmm.ZMM_RNDVAL;
		asm_data->u.zmm.ZMM_RNDVAL_OVER_LARGE_BASE = asm_data->u.zmm.ZMM_RNDVAL / big_word - asm_data->u.zmm.ZMM_RNDVAL;
		asm_data->u.zmm.ZMM_RNDVAL_OVER_SMALL_BASE = asm_data->u.zmm.ZMM_RNDVAL / small_word - asm_data->u.zmm.ZMM_RNDVAL;

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 64-byte boundaries */

		ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
		gwdata->pass1_var_data = tables;
		tables = zr4dwpn_build_pass1_table (gwdata, tables);
		tables = round_to_cache_line (tables);
		/* The wrapper for "one-pass" FFTs does not use a fixed sin/cos table, but it does access the variable data using sincos2 */
		if (gwdata->PASS1_SIZE == 0) {
			asm_data->sincos2 = gwdata->pass1_var_data;
		} else {
			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->sincos2 = tables;
			tables = zr4dwpn_build_fixed_pass1_table (gwdata, tables);
			asm_data->sincos2 = share_sincos_data (gwdata, FIXED_PASS1_SINCOS_DATA, asm_data->sincos2, (char *) tables - (char *) asm_data->sincos2);
			tables = round_to_cache_line (tables);
		}

/* Build the sin/cos table used in complex pass 2 blocks */
/* The pass 2 tables are the same as for a traditional radix-4 FFT */		

		ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
		asm_data->xsincos_complex = tables;
		tables = zr4_build_pass2_complex_table (gwdata, tables);
		asm_data->xsincos_complex = share_sincos_data (gwdata, PASS2_COMPLEX_SINCOS_DATA, asm_data->xsincos_complex, (char *) tables - (char *) asm_data->xsincos_complex);
		tables = round_to_cache_line (tables);
		ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
		asm_data->sincos3 = tables;
		tables = zr4_build_pass2_real_table (gwdata, tables);
		asm_data->sincos3 = share_sincos_data (gwdata, PASS2_REAL_SINCOS_DATA, asm_data->sincos3, (char *) tables - (char *) asm_data->sincos3);
		tables = round_to_cache_line (tables);

/* Allocate a table for carries.  For better distribution of data in the caches, make this table contiguous with all */
/* other data used in the first pass (scratch area, normalization tables, etc.)  Note that we put the tables that */
/* are only partly loaded (column multipliers and big/lit table) after the tables that are loaded throughout the first pass. */

		{
			int	i, carry_table_size;
			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->carries = tables;
			if (gwdata->PASS1_SIZE == 0) carry_table_size = gwdata->FFTLEN / 8; // Every 8th one-pass FFT element needs a carry
			else carry_table_size = gwdata->PASS1_SIZE;	// Pass 1 size - each needs a carry
			for (i = 0; i < carry_table_size; i++) *tables++ = asm_data->u.zmm.ZMM_RNDVAL;
		}

/* Build the group muliplier normalization table.  Keep this table contiguous with other data used in pass 1. */

		ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
		asm_data->norm_grp_mults = tables;
		tables = zr4dwpn_build_norm_table (gwdata, tables);
		tables = round_to_cache_line (tables);

/* Reserve room for the pass 1 scratch area. */

		if (gwdata->SCRATCH_SIZE) {
			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->scratch_area = tables;
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
			tables = round_to_cache_line (tables);
		}

/* Build the table of big vs. little flags.  Build the table of fudge factor flags */

		ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
		tables = zr4dwpn_build_biglit_table (gwdata, tables);
		tables = round_to_cache_line (tables);
		tables = zr4dwpn_build_fudge_table (gwdata, tables);
		tables = round_to_cache_line (tables);

#ifdef GDEBUG_MEM
		{
		char buf[80];
		sprintf (buf, "FFTlen: %d, clm: %d\n", (int) gwdata->FFTLEN, (int) gwdata->PASS1_CACHE_LINES); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, scratch area: %d\n", (int) gwdata->FFTLEN, (int) gwdata->SCRATCH_SIZE); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass1 var s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->sincos2 - (intptr_t) gwdata->pass1_var_data)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass1 fixed s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->xsincos_complex - (intptr_t) asm_data->sincos2)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass2 complex s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->sincos3 - (intptr_t) asm_data->xsincos_complex)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass2 real s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->carries - (intptr_t) asm_data->sincos3)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, carries: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->norm_grp_mults - (intptr_t) asm_data->carries)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, norm grp mults: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->scratch_area - (intptr_t) asm_data->norm_grp_mults)); OutputBoth (0, buf);
		}
#endif

/* Create offsets for carry propagation code to step through eight FFT data elements */

		asm_data->u.zmm.ZMM_SRC_INCR = (intptr_t) addr_offset (gwdata, 1) - (intptr_t) addr_offset (gwdata, 0);
		ASSERTG (asm_data->u.zmm.ZMM_SRC_INCR == (intptr_t) addr_offset (gwdata, 2) - (intptr_t) addr_offset (gwdata, 1));
		ASSERTG (asm_data->u.zmm.ZMM_SRC_INCR == (intptr_t) addr_offset (gwdata, 3) - (intptr_t) addr_offset (gwdata, 2));
		ASSERTG (asm_data->u.zmm.ZMM_SRC_INCR == (intptr_t) addr_offset (gwdata, 4) - (intptr_t) addr_offset (gwdata, 3));
		ASSERTG (asm_data->u.zmm.ZMM_SRC_INCR == (intptr_t) addr_offset (gwdata, 5) - (intptr_t) addr_offset (gwdata, 4));
		ASSERTG (asm_data->u.zmm.ZMM_SRC_INCR == (intptr_t) addr_offset (gwdata, 6) - (intptr_t) addr_offset (gwdata, 5));
		ASSERTG (asm_data->u.zmm.ZMM_SRC_INCR == (intptr_t) addr_offset (gwdata, 7) - (intptr_t) addr_offset (gwdata, 6));
	}

/* Initialize tables for AVX FFT code */

	else if (gwdata->cpu_flags & CPU_AVX) {

/* Initialize tables for the one pass radix-4 FFT assembly code. */

		if (gwdata->PASS2_SIZE == 0) {

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->sincos1 = tables;
			tables = yr4_build_onepass_sincos_table (gwdata, tables);

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->norm_col_mults = tables;
			tables = yr4_build_onepass_norm_table (gwdata, tables);

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->norm_biglit_array = tables;
			tables = yr4_build_onepass_biglit_table (gwdata, tables);
		}

/* Initialize tables for an r4delay FFT with partial normalization (r4dwpn). */

		else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 64-byte boundaries */

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			tables = yr4dwpn_build_pass1_table (gwdata, tables);
			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->sincos2 = tables;
			tables = yr4dwpn_build_fixed_pass1_table (gwdata, tables);
			asm_data->sincos2 = share_sincos_data (gwdata, FIXED_PASS1_SINCOS_DATA, asm_data->sincos2, (char *) tables - (char *) asm_data->sincos2);

/* Build the sin/cos table used in complex pass 2 blocks */
/* The pass 2 tables are the same as for a traditional radix-4 FFT */		

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->xsincos_complex = tables;
			tables = yr4_build_pass2_complex_table (gwdata, tables);
			asm_data->xsincos_complex = share_sincos_data (gwdata, PASS2_COMPLEX_SINCOS_DATA, asm_data->xsincos_complex, (char *) tables - (char *) asm_data->xsincos_complex);
			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->sincos3 = tables;
			tables = yr4_build_pass2_real_table (gwdata, tables);
			asm_data->sincos3 = share_sincos_data (gwdata, PASS2_REAL_SINCOS_DATA, asm_data->sincos3, (char *) tables - (char *) asm_data->sincos3);

/* Allocate a table for carries.  Init with YMM_BIGVAL.  For better distribution of data */
/* in the L2 cache, make this table contiguous with all other data used in the first pass */
/* (scratch area, normalization tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table) after the tables that are loaded throughout the first pass. */

			if (gwdata->PASS2_SIZE) {
				int	i, carry_table_size;
				double	carryval;
				ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
				asm_data->carries = tables;
				carry_table_size = gwdata->PASS1_SIZE;
				if (b == 2) carryval = 3.0 * 131072.0 * 131072.0 * 131072.0;
				else carryval = 0.0;
				for (i = 0; i < carry_table_size; i++)
					if (gwdata->ZERO_PADDED_FFT && (i & 4)) *tables++ = 0.0;
					else *tables++ = carryval;
				tables += (8 - (tables - gwdata->gwnum_memory)) & 7;
			}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->norm_grp_mults = tables;
			tables = yr4dwpn_build_norm_table (gwdata, tables);

/* Reserve room for the pass 1 scratch area. */

			if (gwdata->SCRATCH_SIZE) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
				asm_data->scratch_area = tables;
				tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
			}

/* Build the table of big vs. little flags. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 7) == 0);
			asm_data->norm_biglit_array = tables;
			tables = yr4dwpn_build_biglit_table (gwdata, tables);
			tables += (8 - (tables - gwdata->gwnum_memory)) & 7;

/* Build the column normalization multiplier table. */

			asm_data->norm_col_mults = NULL;
		}

#ifdef GDEBUG_MEM
	if (gwdata->PASS2_SIZE) {
		char buf[80];
		sprintf (buf, "FFTlen: %d, clm: %d, count3: %d, count2: %d\n", (int) gwdata->FFTLEN, (int) gwdata->PASS1_CACHE_LINES, (int) asm_data->count3, (int) asm_data->count2); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, scratch area: %d\n", (int) gwdata->FFTLEN, (int) gwdata->SCRATCH_SIZE); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass1 var s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->sincos2 - (intptr_t) gwdata->pass1_var_data)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass1 fixed s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->xsincos_complex - (intptr_t) asm_data->sincos2)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass2 complex s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->sincos3 - (intptr_t) asm_data->xsincos_complex)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, pass2 real s/c data: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->carries - (intptr_t) asm_data->sincos3)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, carries: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->norm_grp_mults - (intptr_t) asm_data->carries)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, norm grp mults: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) asm_data->scratch_area - (intptr_t) asm_data->norm_grp_mults)); OutputBoth (0, buf);
		sprintf (buf, "FFTlen: %d, norm biglit: %d\n", (int) gwdata->FFTLEN, (int) ((intptr_t) tables - (intptr_t) asm_data->norm_biglit_array)); OutputBoth (0, buf);
	}
#endif

/* Create offsets for carry propagation code to step through norm array */

		asm_data->u.ymm.YMM_SRC_INCR7 = (intptr_t) addr_offset (gwdata, 7) - (intptr_t) addr_offset (gwdata, 6);
		asm_data->u.ymm.YMM_SRC_INCR6 = (intptr_t) addr_offset (gwdata, 6) - (intptr_t) addr_offset (gwdata, 5);
		asm_data->u.ymm.YMM_SRC_INCR5 = (intptr_t) addr_offset (gwdata, 5) - (intptr_t) addr_offset (gwdata, 4);
		asm_data->u.ymm.YMM_SRC_INCR4 = (intptr_t) addr_offset (gwdata, 4) - (intptr_t) addr_offset (gwdata, 3);
		asm_data->u.ymm.YMM_SRC_INCR3 = (intptr_t) addr_offset (gwdata, 3) - (intptr_t) addr_offset (gwdata, 2);
		asm_data->u.ymm.YMM_SRC_INCR2 = (intptr_t) addr_offset (gwdata, 2) - (intptr_t) addr_offset (gwdata, 1);
		asm_data->u.ymm.YMM_SRC_INCR1 = (intptr_t) addr_offset (gwdata, 1) - (intptr_t) addr_offset (gwdata, 0);
	}

/* Initialize tables for SSE2 FFTs */

	else if (gwdata->cpu_flags & CPU_SSE2) {

/* Initialize tables for the radix-4 FFT assembly code.  These are guaranteed to be */
/* two-pass FFTs as we use the older FFT code for one pass FFTs. */

		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 128-byte boundaries */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			tables = r4_build_pass1_table (gwdata, tables);

/* Build the sin/cos table used in complex pass 2 blocks */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->xsincos_complex = tables;
			tables = r4_build_pass2_complex_table (gwdata, tables);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos3 = tables;
			tables = r4_build_pass2_real_table (gwdata, tables);

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table) after the tables that are */
/* loaded throughout the first pass. */

			if (gwdata->PASS2_SIZE) {
				int	i, carry_table_size;
				double	xmm_bigval;
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->carries = tables;
				carry_table_size = gwdata->PASS1_SIZE * 2;
				xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
				for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
			}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_grp_mults = tables;
			tables = r4_build_norm_table (gwdata, tables, 0);

/* Reserve room for the pass 1 scratch area. */

			if (gwdata->SCRATCH_SIZE) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->scratch_area = tables;
				tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
			}

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_biglit_array = tables;
			tables = r4_build_biglit_table (gwdata, tables);
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_col_mults = tables;
			tables = r4_build_norm_table (gwdata, tables, 1);
		}

/* Initialize tables for a modified radix-4 FFT.  This particular FFT */
/* uses a radix-8 building block when there are an odd number of FFT */
/* levels in a pass. */

		else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DELAYED) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 128-byte boundaries */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			tables = r4delay_build_pass1_table (gwdata, tables);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos1 = tables;
			tables = r4delay_build_fixed_premult_table (gwdata, tables);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos2 = tables;
			tables = r4delay_build_fixed_pass1_table (gwdata, tables);

/* Build the sin/cos table used in complex pass 2 blocks */
/* The pass 2 tables are the same as for a traditional radix-4 FFT */		

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->xsincos_complex = tables;
			tables = r4_build_pass2_complex_table (gwdata, tables);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos3 = tables;
			tables = r4_build_pass2_real_table (gwdata, tables);

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

			if (gwdata->PASS2_SIZE) {
				int	i, carry_table_size;
				double	xmm_bigval;
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->carries = tables;
				carry_table_size = gwdata->PASS1_SIZE * 2;
				xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
				for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
			}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_grp_mults = tables;
			tables = r4_build_norm_table (gwdata, tables, 0);

/* Reserve room for the pass 1 scratch area. */

			if (gwdata->SCRATCH_SIZE) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->scratch_area = tables;
				tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
			}

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_biglit_array = tables;
			tables = r4_build_biglit_table (gwdata, tables);
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_col_mults = tables;
			tables = r4_build_norm_table (gwdata, tables, 1);
		}

/* Initialize tables for an r4delay FFT with partial normalization (r4dwpn). */

		else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {

/* Build sin/cos and premultiplier tables used in pass 1 of two pass FFTs. */
/* For best prefetching, make sure tables remain on 128-byte boundaries */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			tables = r4dwpn_build_pass1_table (gwdata, tables);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos1 = tables;
			tables = r4delay_build_fixed_premult_table (gwdata, tables);
			asm_data->sincos1 = share_sincos_data (gwdata, FIXED_PASS1_SINCOS_DATA, asm_data->sincos1, (char *) tables - (char *) asm_data->sincos1);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos2 = tables;
			tables = r4delay_build_fixed_pass1_table (gwdata, tables);
			asm_data->sincos2 = share_sincos_data (gwdata, FIXED_PASS1_SINCOS_DATA, asm_data->sincos2, (char *) tables - (char *) asm_data->sincos2);

/* Build the sin/cos table used in complex pass 2 blocks */
/* The pass 2 tables are the same as for a traditional radix-4 FFT */		

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->xsincos_complex = tables;
			tables = r4_build_pass2_complex_table (gwdata, tables);
			asm_data->xsincos_complex = share_sincos_data (gwdata, PASS2_COMPLEX_SINCOS_DATA, asm_data->xsincos_complex, (char *) tables - (char *) asm_data->xsincos_complex);
			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos3 = tables;
			tables = r4_build_pass2_real_table (gwdata, tables);
			asm_data->sincos3 = share_sincos_data (gwdata, PASS2_REAL_SINCOS_DATA, asm_data->sincos3, (char *) tables - (char *) asm_data->sincos3);

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

			if (gwdata->PASS2_SIZE) {
				int	i, carry_table_size;
				double	xmm_bigval;
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->carries = tables;
				carry_table_size = gwdata->PASS1_SIZE * 2;
				xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
				for (i = 0; i < carry_table_size; i++) *tables++ = xmm_bigval;
				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
			}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_grp_mults = tables;
			tables = r4dwpn_build_norm_table (gwdata, tables);

/* Reserve room for the pass 1 scratch area. */

			if (gwdata->SCRATCH_SIZE) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->scratch_area = tables;
				tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
			}

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_biglit_array = tables;
			tables = r4dwpn_build_biglit_table (gwdata, tables);
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

			asm_data->norm_col_mults = NULL;
		}

/* Initialize tables for the home-grown SSE2 FFT code. */

		else {

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but build_sin_cos_table wants the number of reals in a section. */
/* However, we build a 1/4-sized table by mixing some of the sin/cos */
/* data into the premultiplier table.  So, divide pass2_size by 2 instead of */
/* multiplying pass2_size by 2. */

/* For best prefetching, make sure tables remain on 128-byte boundaries */

			if (gwdata->PASS2_SIZE) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->u.xmm.pass2_premults = tables;
				tables = hg_build_premult_table (gwdata, tables);

/* Build the rest of the tables */

				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->xsincos_complex = tables;
				tables = hg_build_sin_cos_table (tables, gwdata->PASS2_SIZE/2, 0, 1);

				if (!gwdata->ALL_COMPLEX_FFT) {
					ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
					asm_data->u.xmm.sincos6 = tables;
					tables = hg_build_sin_cos_table (tables, gwdata->PASS2_SIZE * 4, 1, 2);
					asm_data->u.xmm.sincos7 = tables;
					tables = hg_build_sin_cos_table (tables, gwdata->PASS2_SIZE, 1, 1);
					tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
				}

				asm_data->u.xmm.sincos8 = asm_data->u.xmm.sincos7;
				asm_data->u.xmm.sincos9 = asm_data->u.xmm.sincos7;
				asm_data->u.xmm.sincos10 = asm_data->u.xmm.sincos7;
				asm_data->u.xmm.sincos11 = asm_data->u.xmm.sincos7;
			}

/* Allocate a table for carries.  Init with XMM_BIGVAL.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with all other data used in the first pass (scratch area, normalization */
/* tables, etc.)  Note that we put the tables that are only partly loaded */
/* (column multipliers and big/lit table after the tables that are */
/* loaded throughout the first pass. */

			if (gwdata->PASS2_SIZE) {
				int	i, carry_table_size;
				double	xmm_bigval;
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->carries = tables;
				carry_table_size = gwdata->PASS1_SIZE * 2;
				xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
				for (i = 0; i < carry_table_size; i++)
					*tables++ = xmm_bigval;
				tables += (16 - (tables - gwdata->gwnum_memory)) & 15;
			}

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_grp_mults = tables;
			tables = hg_build_norm_table (gwdata, tables, 0);

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

			if (gwdata->ALL_COMPLEX_FFT) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->plus1_premults = tables;
				tables = hg_build_plus1_table (gwdata, tables);
			}

/* Reserve room for the pass 1 scratch area. */

			if (gwdata->SCRATCH_SIZE) {
				ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
				asm_data->scratch_area = tables;
				tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);
			}

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->sincos1 = tables;
			tables = hg_build_sin_cos_table (tables, pass1_size, !gwdata->ALL_COMPLEX_FFT, gwdata->PASS2_SIZE == 0 ? 2 : 1);

			if (gwdata->PASS2_SIZE && pass1_size == pow_two_above_or_equal (pass1_size))
				asm_data->sincos2 = asm_data->sincos1;
			else {
				asm_data->sincos2 = tables;
				tables = hg_build_sin_cos_table (tables, pass1_size/4, !gwdata->ALL_COMPLEX_FFT, 1);
			}

			if (pass1_size == pow_two_above_or_equal (pass1_size)) {
				asm_data->sincos3 = asm_data->sincos2;
				asm_data->sincos4 = asm_data->sincos2;
				asm_data->sincos5 = asm_data->sincos2;
			} else {
				asm_data->sincos3 = tables;
				tables = hg_build_sin_cos_table (tables, pass1_size/16, !gwdata->ALL_COMPLEX_FFT, 1);
				asm_data->sincos4 = tables;
				tables = hg_build_sin_cos_table (tables, pass1_size/64, !gwdata->ALL_COMPLEX_FFT, 1);
				asm_data->sincos5 = tables;
				tables = hg_build_sin_cos_table (tables, pass1_size/256, !gwdata->ALL_COMPLEX_FFT, 1);
			}
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the table of big vs. little flags.  This cannot be last table */
/* built as xnorm_2d macro reads 2 bytes past the end of the array. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_biglit_array = tables;
			tables = hg_build_biglit_table (gwdata, tables);
			tables += (16 - (tables - gwdata->gwnum_memory)) & 15;

/* Build the column normalization multiplier table. */

			ASSERTG (((tables - gwdata->gwnum_memory) & 15) == 0);
			asm_data->norm_col_mults = tables;
			tables = hg_build_norm_table (gwdata, tables, 1);
		}
	}

/* Initialze table for the x87 assembly code. */

#ifndef X86_64

	else {

/* Allocate a table for carries.  Init with zero.  For best */
/* distribution of data in the L2 cache, make this table contiguous */
/* with the scratch area which is also used in the first pass. */

		if (gwdata->PASS2_SIZE) {
			int	i, carry_table_size;
			asm_data->carries = tables;
			carry_table_size = gwdata->PASS1_SIZE;
			for (i = 0; i < carry_table_size; i++) *tables++ = 0.0;
		}

/* Reserve room for the pass 1 scratch area. */

		asm_data->scratch_area = tables;
		if (gwdata->SCRATCH_SIZE)
			tables = (double *) ((char *) tables + gwdata->SCRATCH_SIZE);

/* Build the group muliplier normalization table.  Keep this table */
/* contiguous with other data used in pass 1. */

		asm_data->norm_grp_mults = tables;
		tables = x87_build_norm_table (gwdata, tables, 0);

/* Build sin/cos tables used in pass 1.  If FFTLEN is a power of two, */
/* many of the sin/cos tables can be shared. */

		asm_data->sincos1 = tables;
		tables = x87_build_sin_cos_table (tables, pass1_size, !gwdata->ALL_COMPLEX_FFT);

		if (pass1_size == pow_two_above_or_equal (pass1_size)) {
			asm_data->sincos2 = asm_data->sincos1;
			asm_data->sincos3 = asm_data->sincos1;
			asm_data->sincos4 = asm_data->sincos1;
			asm_data->sincos5 = asm_data->sincos1;
		} else {
			asm_data->sincos2 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/4, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos3 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/16, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos4 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/64, !gwdata->ALL_COMPLEX_FFT);
			asm_data->sincos5 = tables;
			tables = x87_build_sin_cos_table (tables, pass1_size/256, !gwdata->ALL_COMPLEX_FFT);
		}

/* Build sin/cos and premultiplier tables used in pass 2 of two pass FFTs */
/* Remember that pass2_size is the number of complex values in a pass 2 */
/* section, but x87_build_sin_cos_table wants the number of reals in a section. */

		if (gwdata->PASS2_SIZE) {
			asm_data->u.x87.pass2_premults = tables;
			tables = x87_build_premult_table (gwdata, tables);
			asm_data->xsincos_complex = tables;
			tables = x87_build_sin_cos_table (tables, gwdata->PASS2_SIZE*2, 0);

			if (!gwdata->ALL_COMPLEX_FFT) {
				asm_data->u.x87.sincos6 = tables;
				tables = x87_build_sin_cos_table (tables, gwdata->PASS2_SIZE*2, 1);
				asm_data->u.x87.sincos7 = asm_data->u.x87.sincos6;
				asm_data->u.x87.sincos8 = asm_data->u.x87.sincos6;
				asm_data->u.x87.sincos9 = asm_data->u.x87.sincos6;
				asm_data->u.x87.sincos10 = asm_data->u.x87.sincos6;
			}
		}

/* Build the plus1-pre-multiplier table (complex weights applied when c > 0 */
/* and we are doing a all-complex FFT rather than emulating it with a */
/* zero-padded FFT. */

		if (gwdata->ALL_COMPLEX_FFT) {
			asm_data->plus1_premults = tables;
			tables = x87_build_plus1_table (gwdata, tables);
		}

/* Build the column normalization multiplier table. */

		asm_data->norm_col_mults = tables;
		tables = x87_build_norm_table (gwdata, tables, 1);

/* Build the table of big vs. little flags. */

		asm_data->norm_biglit_array = tables;
		tables = x87_build_biglit_table (gwdata, tables);
	}

#endif

/* Return "impossible" internal errors from building tables */

	if (gwdata->GWERROR) return (gwdata->GWERROR);

/* Finish verifying table size */

#ifdef GDEBUG
	if (!gwdata->ZERO_PADDED_FFT && !gwdata->RATIONAL_FFT) {
		char buf[80];
		intptr_t mem_used = (intptr_t) tables - (intptr_t) gwdata->gwnum_memory;
		if (mem_used != mem_needed) {
			// AVX-512 FFTs (due to rounding compressed fudges/biglits within each variable block to a cache line boundary)
			// may allocate a little too much memory when clm != 1.  Allow a small discrepancy without complaining.
			if (gwdata->cpu_flags & CPU_AVX512F && gwdata->PASS1_CACHE_LINES > 8) {
				if ((int) mem_used > (int) mem_needed || (int) mem_used < (int) (mem_needed * .98)) {
					sprintf (buf, "FFTlen: %d, mem needed should be approximately: %d\n",
						 (int) gwdata->FFTLEN, (int) (mem_used - gwdata->SCRATCH_SIZE));
					OutputBoth (0, buf);
				}
			} else {
				sprintf (buf, "FFTlen: %d, mem needed should be: %d\n",
					 (int) gwdata->FFTLEN, (int) (mem_used - gwdata->SCRATCH_SIZE));
				OutputBoth (0, buf);
			}
		}
	}
#endif

/* Copy the count of cache lines grouped in pass 1.  This affects how we build the normalization tables. */
/* Note that cache line sizes are different in the x87 (16 bytes) and AVX-512/AVX/SSE2 code (64 bytes). */
/* In the assembly code this is called clm or cache_line_multiplier. */

	asm_data->cache_line_multiplier = gwdata->PASS1_CACHE_LINES;

/* The sumout value is FFTLEN/2 times larger than it should be.  Create an */
/* inverse to properly calculate the sumout when a multiplication ends. */

	asm_data->ttmp_ff_inv = 2.0 / (double) gwdata->FFTLEN;

/* Now compute a number of constants the assembly code needs.  These will be copied to properly aligned (SSE2 constants */
/* must be on 16-byte boundaries, AVX constants must be on 32-byte boundaries, AVX-512 are unaligned due to cheap vbroadcastsd) */
/* and grouped (for better cache line locality) assembly global variables. */

	gwasm_constants ((double *) &asm_values);

/* Init constants.  Some of these could be true global variables but I elected not */
/* to so that constants could be close to other variables used at the */
/* same time.  This might free up a cache line or two. */

/* Init constants for AVX-512 FFT routines */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		asm_data->u.zmm.ZMM_HALF = 0.5;
		asm_data->u.zmm.ZMM_TWO = 2.0;
		asm_data->u.zmm.ZMM_ONE = 1.0;
		asm_data->u.zmm.ZMM_SQRTHALF = sqrt (0.5);
		asm_data->u.zmm.ZMM_SQRT2 = sqrt (2.0);
		asm_data->u.zmm.ZMM_ABSVAL[0] = 0xFFFFFFFF;
		asm_data->u.zmm.ZMM_ABSVAL[1] = 0x7FFFFFFF;
		asm_data->u.zmm.ZMM_B = (double) b;
		asm_data->u.zmm.ZMM_ONE_OVER_B = 1.0 / asm_data->u.zmm.ZMM_B;

/* Negate c and store as a double */

		asm_data->u.zmm.ZMM_MINUS_C = (double) -c;

/* Split k for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k if k * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		if (k * big_word < 562949953421312.0) {
			asm_data->u.zmm.ZMM_K_LO = k;
			asm_data->u.zmm.ZMM_K_HI_OVER_SMALL_BASE = 0.0;
			asm_data->u.zmm.ZMM_K_HI_OVER_LARGE_BASE = 0.0;
		} else {
			asm_data->u.zmm.ZMM_K_HI_OVER_LARGE_BASE = floor (k / big_word);
			asm_data->u.zmm.ZMM_K_HI_OVER_SMALL_BASE = floor (k / big_word) * (double) b;
			asm_data->u.zmm.ZMM_K_LO = k - asm_data->u.zmm.ZMM_K_HI_OVER_LARGE_BASE * big_word;
		}
		gwsetmulbyconst (gwdata, gwdata->maxmulbyconst);

/* Copy computed assembly constants */

		asm_data->u.zmm.ZMM_P951 = asm_values[0];
		asm_data->u.zmm.ZMM_P588_P951 = asm_values[1];
		asm_data->u.zmm.ZMM_P309 = asm_values[2];
		asm_data->u.zmm.ZMM_P809 = -asm_values[6];
		asm_data->u.zmm.ZMM_P866 = asm_values[8];
		asm_data->u.zmm.ZMM_P623 = asm_values[14];
		asm_data->u.zmm.ZMM_P223 = -asm_values[17];
		asm_data->u.zmm.ZMM_P901 = -asm_values[18];
		asm_data->u.zmm.ZMM_P383 = asm_values[21];
		asm_data->u.zmm.ZMM_P383_1 = asm_values[21];
		asm_data->u.zmm.ZMM_P975 = asm_values[11];
		asm_data->u.zmm.ZMM_P434_P975 = asm_values[12];
		asm_data->u.zmm.ZMM_P782_P975 = asm_values[42];
		asm_data->u.zmm.ZMM_P901_P975 = asm_values[38];
		asm_data->u.zmm.ZMM_P623_P975 = asm_values[39];
		asm_data->u.zmm.ZMM_P223_P975 = asm_values[40];
		asm_data->u.zmm.ZMM_P1_P975 = asm_values[41];
		asm_data->u.zmm.ZMM_P259_P707 = asm_values[28];
		asm_data->u.zmm.ZMM_P966_P707 = asm_values[29];
		asm_data->u.zmm.ZMM_P924_P383 = asm_values[30];
		asm_data->u.zmm.ZMM_P981_P195 = asm_values[31];
		asm_data->u.zmm.ZMM_P195 = asm_values[33];
		asm_data->u.zmm.ZMM_P831_P556 = asm_values[34];
		asm_data->u.zmm.ZMM_P556_P195 = asm_values[37];

/* Swizzling indices used in vpermt2pd */

		asm_data->u.zmm.ZMM_PERMUTE1 = 0x0c040e0608000a02ULL;
		asm_data->u.zmm.ZMM_PERMUTE2 = 0x0d050f0709010b03ULL;
	}

/* Init constants for AVX FFT routines */

	else if (gwdata->cpu_flags & CPU_AVX) {
		int	i;

		asm_data->u.ymm.YMM_HALF[0] = asm_data->u.ymm.YMM_HALF[1] =
		asm_data->u.ymm.YMM_HALF[2] = asm_data->u.ymm.YMM_HALF[3] = 0.5;
		asm_data->u.ymm.YMM_ONE[0] = asm_data->u.ymm.YMM_ONE[1] =
		asm_data->u.ymm.YMM_ONE[2] = asm_data->u.ymm.YMM_ONE[3] = 1.0;
		asm_data->u.ymm.YMM_TWO[0] = asm_data->u.ymm.YMM_TWO[1] =
		asm_data->u.ymm.YMM_TWO[2] = asm_data->u.ymm.YMM_TWO[3] = 2.0;
		asm_data->u.ymm.YMM_SQRTHALF[0] = asm_data->u.ymm.YMM_SQRTHALF[1] =
		asm_data->u.ymm.YMM_SQRTHALF[2] = asm_data->u.ymm.YMM_SQRTHALF[3] = sqrt (0.5);
		asm_data->u.ymm.YMM_ABSVAL[0] = asm_data->u.ymm.YMM_ABSVAL[2] =
		asm_data->u.ymm.YMM_ABSVAL[4] = asm_data->u.ymm.YMM_ABSVAL[6] = 0xFFFFFFFF;
		asm_data->u.ymm.YMM_ABSVAL[1] = asm_data->u.ymm.YMM_ABSVAL[3] =
		asm_data->u.ymm.YMM_ABSVAL[5] = asm_data->u.ymm.YMM_ABSVAL[7] = 0x7FFFFFFF;

/* Compute the AVX (53-bit) rounding constants */

		asm_data->u.ymm.YMM_BIGVAL[0] = asm_data->u.ymm.YMM_BIGVAL[1] =
		asm_data->u.ymm.YMM_BIGVAL[2] = asm_data->u.ymm.YMM_BIGVAL[3] = 3.0 * pow (2.0, 51.0);
		asm_data->u.ymm.YMM_BIGBIGVAL[0] = asm_data->u.ymm.YMM_BIGBIGVAL[1] =
		asm_data->u.ymm.YMM_BIGBIGVAL[2] = asm_data->u.ymm.YMM_BIGBIGVAL[3] = big_word * asm_data->u.ymm.YMM_BIGVAL[0];

/* Negate c and store as a double */

		asm_data->u.ymm.YMM_MINUS_C[0] = asm_data->u.ymm.YMM_MINUS_C[1] =
		asm_data->u.ymm.YMM_MINUS_C[2] = asm_data->u.ymm.YMM_MINUS_C[3] = (double) -c;

/* Compute constant that converts fft_weight_over_fftlen found in the */
/* two-to-minus-phi tables into the true fft_weight value.  This is usually */
/* FFTLEN / 2, but when doing a non-zero-padded FFT this is FFTLEN / 2k. */

		asm_data->u.ymm.YMM_NORM012_FF[0] = asm_data->u.ymm.YMM_NORM012_FF[1] =
		asm_data->u.ymm.YMM_NORM012_FF[2] = asm_data->u.ymm.YMM_NORM012_FF[3] =
			(gwdata->ZERO_PADDED_FFT) ?
				(double) (gwdata->FFTLEN / 2) :
				(double) (gwdata->FFTLEN / 2) / k;

/* Split k for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k if k * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		if (k * big_word < 562949953421312.0) {
			asm_data->u.ymm.YMM_K_HI[0] = asm_data->u.ymm.YMM_K_HI[1] =
			asm_data->u.ymm.YMM_K_HI[2] = asm_data->u.ymm.YMM_K_HI[3] = 0.0;
			asm_data->u.ymm.YMM_K_LO[0] = asm_data->u.ymm.YMM_K_LO[1] =
			asm_data->u.ymm.YMM_K_LO[2] = asm_data->u.ymm.YMM_K_LO[3] = k;
		} else {
			asm_data->u.ymm.YMM_K_HI[0] = asm_data->u.ymm.YMM_K_HI[1] =
			asm_data->u.ymm.YMM_K_HI[2] = asm_data->u.ymm.YMM_K_HI[3] = floor (k / big_word) * big_word;
			asm_data->u.ymm.YMM_K_LO[0] = asm_data->u.ymm.YMM_K_LO[1] =
			asm_data->u.ymm.YMM_K_LO[2] = asm_data->u.ymm.YMM_K_LO[3] = k - asm_data->u.ymm.YMM_K_HI[0];
		}
		gwsetmulbyconst (gwdata, gwdata->maxmulbyconst);

/* Compute the normalization constants indexed by biglit array entries.  This simple */
/* method is used in one pass FFTs. */

		for (i = 0; i <= 63; i++) {
			int	ymmword, bit;
			ymmword = i >> 2;
			bit = i & 3;
			if (! (ymmword & (1 << bit))) {		/* A small word */
				asm_data->u.ymm.YMM_LIMIT_INVERSE[i] = 1.0 / small_word;	/* Lower limit inverse */
				if (b > 2) temp = small_word;					/* Compute lower limit bigmax */
				else temp = small_word * asm_data->u.ymm.YMM_BIGVAL[0] - asm_data->u.ymm.YMM_BIGVAL[0];
				asm_data->u.ymm.YMM_LIMIT_BIGMAX[i] = temp;
			} else {				/* A big word */
				asm_data->u.ymm.YMM_LIMIT_INVERSE[i] = 1.0 / big_word;		/* Upper limit inverse */
				if (b > 2) temp = big_word;					/* Compute upper limit bigmax */
				else temp = big_word * asm_data->u.ymm.YMM_BIGVAL[0] - asm_data->u.ymm.YMM_BIGVAL[0];
				asm_data->u.ymm.YMM_LIMIT_BIGMAX[i] = temp;
			}
		}

/* Two-pass FFTs use a clever mechanism to reduce big/lit flags data. */
/* Output the valid combinations that were precomputed in r4dwpn_build_biglit_table. */

		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
			double	inv1, inv2, bmax1, bmax2;

			inv1 = asm_data->u.ymm.YMM_LIMIT_INVERSE[0];
			inv2 = asm_data->u.ymm.YMM_LIMIT_INVERSE[63];
			bmax1 = asm_data->u.ymm.YMM_LIMIT_BIGMAX[0];
			bmax2 = asm_data->u.ymm.YMM_LIMIT_BIGMAX[63];
			for (i = 0; i < 48; i++) {
				asm_data->u.ymm.YMM_LIMIT_INVERSE[i*4] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? inv2 : inv1;
				asm_data->u.ymm.YMM_LIMIT_INVERSE[i*4+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? inv2 : inv1;
				asm_data->u.ymm.YMM_LIMIT_INVERSE[i*4+2] = (((char *)gwdata->ASM_TIMERS)[i] & 4) ? inv2 : inv1;
				asm_data->u.ymm.YMM_LIMIT_INVERSE[i*4+3] = (((char *)gwdata->ASM_TIMERS)[i] & 8) ? inv2 : inv1;
				asm_data->u.ymm.YMM_LIMIT_BIGMAX[i*4] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? bmax2 : bmax1;
				asm_data->u.ymm.YMM_LIMIT_BIGMAX[i*4+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? bmax2 : bmax1;
				asm_data->u.ymm.YMM_LIMIT_BIGMAX[i*4+2] = (((char *)gwdata->ASM_TIMERS)[i] & 4) ? bmax2 : bmax1;
				asm_data->u.ymm.YMM_LIMIT_BIGMAX[i*4+3] = (((char *)gwdata->ASM_TIMERS)[i] & 8) ? bmax2 : bmax1;
			}
		}

/* Copy computed assembly constants */

		asm_data->u.ymm.YMM_P951[0] = asm_data->u.ymm.YMM_P951[1] =
		asm_data->u.ymm.YMM_P951[2] = asm_data->u.ymm.YMM_P951[3] = asm_values[0];
		asm_data->u.ymm.YMM_P618[0] = asm_data->u.ymm.YMM_P618[1] =
		asm_data->u.ymm.YMM_P618[2] = asm_data->u.ymm.YMM_P618[3] = asm_values[1];
		asm_data->u.ymm.YMM_P309[0] = asm_data->u.ymm.YMM_P309[1] =
		asm_data->u.ymm.YMM_P309[2] = asm_data->u.ymm.YMM_P309[3] = asm_values[2];
		asm_data->u.ymm.YMM_P588[0] = asm_data->u.ymm.YMM_P588[1] =
		asm_data->u.ymm.YMM_P588[2] = asm_data->u.ymm.YMM_P588[3] = asm_values[4];
		asm_data->u.ymm.YMM_P809[0] = asm_data->u.ymm.YMM_P809[1] =
		asm_data->u.ymm.YMM_P809[2] = asm_data->u.ymm.YMM_P809[3] = -asm_values[6];
		asm_data->u.ymm.YMM_P866[0] = asm_data->u.ymm.YMM_P866[1] =
		asm_data->u.ymm.YMM_P866[2] = asm_data->u.ymm.YMM_P866[3] = asm_values[8];
		asm_data->u.ymm.YMM_P975[0] = asm_data->u.ymm.YMM_P975[1] =
		asm_data->u.ymm.YMM_P975[2] = asm_data->u.ymm.YMM_P975[3] = asm_values[11];
		asm_data->u.ymm.YMM_P623[0] = asm_data->u.ymm.YMM_P623[1] =
		asm_data->u.ymm.YMM_P623[2] = asm_data->u.ymm.YMM_P623[3] = asm_values[14];
		asm_data->u.ymm.YMM_P223[0] = asm_data->u.ymm.YMM_P223[1] =
		asm_data->u.ymm.YMM_P223[2] = asm_data->u.ymm.YMM_P223[3] = -asm_values[17];
		asm_data->u.ymm.YMM_P901[0] = asm_data->u.ymm.YMM_P901[1] =
		asm_data->u.ymm.YMM_P901[2] = asm_data->u.ymm.YMM_P901[3] = -asm_values[18];
		asm_data->u.ymm.YMM_P924[0] = asm_data->u.ymm.YMM_P924[1] =
		asm_data->u.ymm.YMM_P924[2] = asm_data->u.ymm.YMM_P924[3] = asm_values[20];
		asm_data->u.ymm.YMM_P383[0] = asm_data->u.ymm.YMM_P383[1] =
		asm_data->u.ymm.YMM_P383[2] = asm_data->u.ymm.YMM_P383[3] = asm_values[21];
		asm_data->u.ymm.YMM_P782[0] = asm_data->u.ymm.YMM_P782[1] =
		asm_data->u.ymm.YMM_P782[2] = asm_data->u.ymm.YMM_P782[3] = asm_values[22];
		asm_data->u.ymm.YMM_P434[0] = asm_data->u.ymm.YMM_P434[1] =
		asm_data->u.ymm.YMM_P434[2] = asm_data->u.ymm.YMM_P434[3] = asm_values[23];
		asm_data->u.ymm.YMM_P975_P434[0] = asm_data->u.ymm.YMM_P975_P434[1] =
		asm_data->u.ymm.YMM_P975_P434[2] = asm_data->u.ymm.YMM_P975_P434[3] = asm_values[24];
		asm_data->u.ymm.YMM_P782_P434[0] = asm_data->u.ymm.YMM_P782_P434[1] =
		asm_data->u.ymm.YMM_P782_P434[2] = asm_data->u.ymm.YMM_P782_P434[3] = asm_values[25];
		asm_data->u.ymm.YMM_P259[0] = asm_data->u.ymm.YMM_P259[1] =
		asm_data->u.ymm.YMM_P259[2] = asm_data->u.ymm.YMM_P259[3] = asm_values[26];
		asm_data->u.ymm.YMM_P966[0] = asm_data->u.ymm.YMM_P966[1] =
		asm_data->u.ymm.YMM_P966[2] = asm_data->u.ymm.YMM_P966[3] = asm_values[27];
		asm_data->u.ymm.YMM_P259_P707[0] = asm_data->u.ymm.YMM_P259_P707[1] =
		asm_data->u.ymm.YMM_P259_P707[2] = asm_data->u.ymm.YMM_P259_P707[3] = asm_values[28];
		asm_data->u.ymm.YMM_P966_P707[0] = asm_data->u.ymm.YMM_P966_P707[1] =
		asm_data->u.ymm.YMM_P966_P707[2] = asm_data->u.ymm.YMM_P966_P707[3] = asm_values[29];
		asm_data->u.ymm.YMM_P924_P383[0] = asm_data->u.ymm.YMM_P924_P383[1] =
		asm_data->u.ymm.YMM_P924_P383[2] = asm_data->u.ymm.YMM_P924_P383[3] = asm_values[30];
	}

/* Init constants for SSE2 code */

	else if (gwdata->cpu_flags & CPU_SSE2) {
		asm_data->u.xmm.XMM_TWO[0] = asm_data->u.xmm.XMM_TWO[1] = 2.0;
		asm_data->u.xmm.XMM_HALF[0] = asm_data->u.xmm.XMM_HALF[1] = 0.5;
		asm_data->u.xmm.XMM_SQRTHALF[0] = asm_data->u.xmm.XMM_SQRTHALF[1] = sqrt (0.5);
		asm_data->u.xmm.XMM_ABSVAL[0] = asm_data->u.xmm.XMM_ABSVAL[2] = 0xFFFFFFFF;
		asm_data->u.xmm.XMM_ABSVAL[1] = asm_data->u.xmm.XMM_ABSVAL[3] = 0x7FFFFFFF;

		asm_data->u.xmm.XMM_TTP_FUDGE[0] = asm_data->u.xmm.XMM_TTP_FUDGE[1] =
		asm_data->u.xmm.XMM_TTP_FUDGE[2] = asm_data->u.xmm.XMM_TTP_FUDGE[3] =
		asm_data->u.xmm.XMM_TTP_FUDGE[4] = asm_data->u.xmm.XMM_TTP_FUDGE[5] =
		asm_data->u.xmm.XMM_TTP_FUDGE[6] = asm_data->u.xmm.XMM_TTP_FUDGE[7] =
		asm_data->u.xmm.XMM_TTP_FUDGE[9] = asm_data->u.xmm.XMM_TTP_FUDGE[11] =
		asm_data->u.xmm.XMM_TTP_FUDGE[13] = asm_data->u.xmm.XMM_TTP_FUDGE[15] =
		asm_data->u.xmm.XMM_TTP_FUDGE[16] = asm_data->u.xmm.XMM_TTP_FUDGE[18] =
		asm_data->u.xmm.XMM_TTP_FUDGE[20] = asm_data->u.xmm.XMM_TTP_FUDGE[22] = 1.0;
		asm_data->u.xmm.XMM_TTP_FUDGE[8] = asm_data->u.xmm.XMM_TTP_FUDGE[10] =
		asm_data->u.xmm.XMM_TTP_FUDGE[12] = asm_data->u.xmm.XMM_TTP_FUDGE[14] =
		asm_data->u.xmm.XMM_TTP_FUDGE[17] = asm_data->u.xmm.XMM_TTP_FUDGE[19] =
		asm_data->u.xmm.XMM_TTP_FUDGE[21] = asm_data->u.xmm.XMM_TTP_FUDGE[23] =
		asm_data->u.xmm.XMM_TTP_FUDGE[24] = asm_data->u.xmm.XMM_TTP_FUDGE[25] =
		asm_data->u.xmm.XMM_TTP_FUDGE[26] = asm_data->u.xmm.XMM_TTP_FUDGE[27] =
		asm_data->u.xmm.XMM_TTP_FUDGE[28] = asm_data->u.xmm.XMM_TTP_FUDGE[29] =
		asm_data->u.xmm.XMM_TTP_FUDGE[30] = asm_data->u.xmm.XMM_TTP_FUDGE[31] = 1.0 / (double) b;

		asm_data->u.xmm.XMM_TTMP_FUDGE[0] = asm_data->u.xmm.XMM_TTMP_FUDGE[1] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[2] = asm_data->u.xmm.XMM_TTMP_FUDGE[3] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[4] = asm_data->u.xmm.XMM_TTMP_FUDGE[5] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[6] = asm_data->u.xmm.XMM_TTMP_FUDGE[7] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[9] = asm_data->u.xmm.XMM_TTMP_FUDGE[11] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[13] = asm_data->u.xmm.XMM_TTMP_FUDGE[15] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[16] = asm_data->u.xmm.XMM_TTMP_FUDGE[18] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[20] = asm_data->u.xmm.XMM_TTMP_FUDGE[22] = 1.0;
		asm_data->u.xmm.XMM_TTMP_FUDGE[8] = asm_data->u.xmm.XMM_TTMP_FUDGE[10] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[12] = asm_data->u.xmm.XMM_TTMP_FUDGE[14] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[17] = asm_data->u.xmm.XMM_TTMP_FUDGE[19] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[21] = asm_data->u.xmm.XMM_TTMP_FUDGE[23] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[24] = asm_data->u.xmm.XMM_TTMP_FUDGE[25] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[26] = asm_data->u.xmm.XMM_TTMP_FUDGE[27] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[28] = asm_data->u.xmm.XMM_TTMP_FUDGE[29] =
		asm_data->u.xmm.XMM_TTMP_FUDGE[30] = asm_data->u.xmm.XMM_TTMP_FUDGE[31] = (double) b;

/* Compute the SSE2 (53-bit) rounding constants */

		asm_data->u.xmm.XMM_BIGVAL[0] = asm_data->u.xmm.XMM_BIGVAL[1] = 3.0 * pow (2.0, 51.0);
		asm_data->u.xmm.XMM_BIGVAL_NEG[0] = asm_data->u.xmm.XMM_BIGVAL_NEG[1] = -asm_data->u.xmm.XMM_BIGVAL[0];
		asm_data->u.xmm.XMM_BIGBIGVAL[0] = asm_data->u.xmm.XMM_BIGBIGVAL[1] = big_word * asm_data->u.xmm.XMM_BIGVAL[0];

/* Negate c and store as a double */

		asm_data->u.xmm.XMM_MINUS_C[0] = asm_data->u.xmm.XMM_MINUS_C[1] = (double) -c;

/* Compute constant that converts fft_weight_over_fftlen found in the */
/* two-to-minus-phi tables into the true fft_weight value.  This is usually */
/* FFTLEN / 2, but when doing a non-zero-padded FFT this is FFTLEN / 2k. */

		asm_data->u.xmm.XMM_NORM012_FF[0] = asm_data->u.xmm.XMM_NORM012_FF[1] =
			(gwdata->ZERO_PADDED_FFT) ?
				(double) (gwdata->FFTLEN / 2) :
				(double) (gwdata->FFTLEN / 2) / k;

/* Split k for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k if k * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		if (k * big_word < 562949953421312.0) {
			asm_data->u.xmm.XMM_K_HI[0] = asm_data->u.xmm.XMM_K_HI[1] = 0.0;
			asm_data->u.xmm.XMM_K_LO[0] = asm_data->u.xmm.XMM_K_LO[1] = k;
		} else {
			asm_data->u.xmm.XMM_K_HI[0] = asm_data->u.xmm.XMM_K_HI[1] = floor (k / big_word) * big_word;
			asm_data->u.xmm.XMM_K_LO[0] = asm_data->u.xmm.XMM_K_LO[1] = k - asm_data->u.xmm.XMM_K_HI[0];
		}
		gwsetmulbyconst (gwdata, gwdata->maxmulbyconst);

/* Compute the normalization constants indexed by biglit array entries */

		asm_data->u.xmm.XMM_LIMIT_INVERSE[0] =
		asm_data->u.xmm.XMM_LIMIT_INVERSE[1] =
		asm_data->u.xmm.XMM_LIMIT_INVERSE[3] =
		asm_data->u.xmm.XMM_LIMIT_INVERSE[4] = 1.0 / small_word;	/* Lower limit inverse */
		asm_data->u.xmm.XMM_LIMIT_INVERSE[2] =
		asm_data->u.xmm.XMM_LIMIT_INVERSE[5] =
		asm_data->u.xmm.XMM_LIMIT_INVERSE[6] =
		asm_data->u.xmm.XMM_LIMIT_INVERSE[7] = 1.0 / big_word;		/* Upper limit inverse */

		/* Compute lower limit bigmax */
		if (b > 2)
			temp = small_word;
		else
			temp = small_word * asm_data->u.xmm.XMM_BIGVAL[0] - asm_data->u.xmm.XMM_BIGVAL[0];
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[0] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[1] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[3] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[4] = temp;

		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[0] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[1] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[3] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[4] = -temp;	/* Negative lower limit bigmax */

		/* Compute upper limit bigmax */
		if (b > 2)
			temp = big_word;
		else
			temp = big_word * asm_data->u.xmm.XMM_BIGVAL[0] - asm_data->u.xmm.XMM_BIGVAL[0];
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[2] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[5] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[6] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX[7] = temp;

		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[2] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[5] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[6] =
		asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[7] = -temp;	/* Negative upper limit bigmax */

		memcpy (asm_data->u.xmm.XMM_LIMIT_INVERSE+8, asm_data->u.xmm.XMM_LIMIT_INVERSE, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_INVERSE+16, asm_data->u.xmm.XMM_LIMIT_INVERSE, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_INVERSE+24, asm_data->u.xmm.XMM_LIMIT_INVERSE, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_BIGMAX+8, asm_data->u.xmm.XMM_LIMIT_BIGMAX, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_BIGMAX+16, asm_data->u.xmm.XMM_LIMIT_BIGMAX, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_BIGMAX+24, asm_data->u.xmm.XMM_LIMIT_BIGMAX, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG+8, asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG+16, asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG, 8 * sizeof (double));
		memcpy (asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG+24, asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG, 8 * sizeof (double));

/* Newer FFTs use a clever mechanism to reduce big/lit flags data. */
/* Output the valid combinations that were precomputed in r4dwpn_build_biglit_table. */

		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
			int	i;
			double	inv1, inv2, bmax1, bmax2;

			inv1 = asm_data->u.xmm.XMM_LIMIT_INVERSE[0];
			inv2 = asm_data->u.xmm.XMM_LIMIT_INVERSE[7];
			bmax1 = asm_data->u.xmm.XMM_LIMIT_BIGMAX[0];
			bmax2 = asm_data->u.xmm.XMM_LIMIT_BIGMAX[7];
			for (i = 0; i < 48; i++) {
				asm_data->u.xmm.XMM_LIMIT_INVERSE[i+i] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? inv2 : inv1;
				asm_data->u.xmm.XMM_LIMIT_INVERSE[i+i+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? inv2 : inv1;
				asm_data->u.xmm.XMM_LIMIT_BIGMAX[i+i] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? bmax2 : bmax1;
				asm_data->u.xmm.XMM_LIMIT_BIGMAX[i+i+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? bmax2 : bmax1;
				asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[i+i] = (((char *)gwdata->ASM_TIMERS)[i] & 1) ? -bmax2 : -bmax1;
				asm_data->u.xmm.XMM_LIMIT_BIGMAX_NEG[i+i+1] = (((char *)gwdata->ASM_TIMERS)[i] & 2) ? -bmax2 : -bmax1;
			}
		}

/* Copy computed assembly constants */

		asm_data->u.xmm.XMM_P951[0] = asm_data->u.xmm.XMM_P951[1] = asm_values[0];
		asm_data->u.xmm.XMM_P618[0] = asm_data->u.xmm.XMM_P618[1] = asm_values[1];
		asm_data->u.xmm.XMM_P309[0] = asm_data->u.xmm.XMM_P309[1] = asm_values[2];
		asm_data->u.xmm.XMM_M262[0] = asm_data->u.xmm.XMM_M262[1] = asm_values[3];
		asm_data->u.xmm.XMM_P588[0] = asm_data->u.xmm.XMM_P588[1] = asm_values[4];
		asm_data->u.xmm.XMM_M162[0] = asm_data->u.xmm.XMM_M162[1] = asm_values[5];
		asm_data->u.xmm.XMM_P809[0] = asm_data->u.xmm.XMM_P809[1] = -asm_values[6];
		asm_data->u.xmm.XMM_M809[0] = asm_data->u.xmm.XMM_M809[1] = asm_values[6];
		asm_data->u.xmm.XMM_M382[0] = asm_data->u.xmm.XMM_M382[1] = asm_values[7];
		asm_data->u.xmm.XMM_P866[0] = asm_data->u.xmm.XMM_P866[1] = asm_values[8];
		asm_data->u.xmm.XMM_P975[0] = asm_data->u.xmm.XMM_P975[1] = asm_values[11];
		asm_data->u.xmm.XMM_P445[0] = asm_data->u.xmm.XMM_P445[1] = asm_values[12];
		asm_data->u.xmm.XMM_P180[0] = asm_data->u.xmm.XMM_P180[1] = asm_values[13];
		asm_data->u.xmm.XMM_P623[0] = asm_data->u.xmm.XMM_P623[1] = asm_values[14];
		asm_data->u.xmm.XMM_M358[0] = asm_data->u.xmm.XMM_M358[1] = asm_values[15];
		asm_data->u.xmm.XMM_P404[0] = asm_data->u.xmm.XMM_P404[1] = asm_values[16];
		asm_data->u.xmm.XMM_P223[0] = asm_data->u.xmm.XMM_P223[1] = -asm_values[17];
		asm_data->u.xmm.XMM_P901[0] = asm_data->u.xmm.XMM_P901[1] = -asm_values[18];
		asm_data->u.xmm.XMM_P924[0] = asm_data->u.xmm.XMM_P924[1] = asm_values[20];
		asm_data->u.xmm.XMM_P383[0] = asm_data->u.xmm.XMM_P383[1] = asm_values[21];
		asm_data->u.xmm.XMM_P782[0] = asm_data->u.xmm.XMM_P782[1] = asm_values[22];
		asm_data->u.xmm.XMM_P434[0] = asm_data->u.xmm.XMM_P434[1] = asm_values[23];
	}

/* Init constants for x87 code */


	else {
		asm_data->u.x87.HALF = 0.5;
		asm_data->u.x87.SQRTHALF = sqrt (0.5);
		asm_data->u.x87.P25 = 0.25;
		asm_data->u.x87.P75 = 0.75;
		asm_data->u.x87.P3 = 3.0;

		asm_data->u.x87.TTP_FUDGE[0] = asm_data->u.x87.TTP_FUDGE[1] =
		asm_data->u.x87.TTP_FUDGE[2] = asm_data->u.x87.TTP_FUDGE[3] =
		asm_data->u.x87.TTP_FUDGE[4] = asm_data->u.x87.TTP_FUDGE[5] =
		asm_data->u.x87.TTP_FUDGE[6] = asm_data->u.x87.TTP_FUDGE[7] =
		asm_data->u.x87.TTP_FUDGE[9] = asm_data->u.x87.TTP_FUDGE[11] =
		asm_data->u.x87.TTP_FUDGE[13] = asm_data->u.x87.TTP_FUDGE[15] =
		asm_data->u.x87.TTP_FUDGE[16] = asm_data->u.x87.TTP_FUDGE[18] =
		asm_data->u.x87.TTP_FUDGE[20] = asm_data->u.x87.TTP_FUDGE[22] = 1.0;
		asm_data->u.x87.TTP_FUDGE[8] = asm_data->u.x87.TTP_FUDGE[10] =
		asm_data->u.x87.TTP_FUDGE[12] = asm_data->u.x87.TTP_FUDGE[14] =
		asm_data->u.x87.TTP_FUDGE[17] = asm_data->u.x87.TTP_FUDGE[19] =
		asm_data->u.x87.TTP_FUDGE[21] = asm_data->u.x87.TTP_FUDGE[23] =
		asm_data->u.x87.TTP_FUDGE[24] = asm_data->u.x87.TTP_FUDGE[25] =
		asm_data->u.x87.TTP_FUDGE[26] = asm_data->u.x87.TTP_FUDGE[27] =
		asm_data->u.x87.TTP_FUDGE[28] = asm_data->u.x87.TTP_FUDGE[29] =
		asm_data->u.x87.TTP_FUDGE[30] = asm_data->u.x87.TTP_FUDGE[31] = 1.0 / (double) b;

		asm_data->u.x87.TTMP_FUDGE[0] = asm_data->u.x87.TTMP_FUDGE[1] =
		asm_data->u.x87.TTMP_FUDGE[2] = asm_data->u.x87.TTMP_FUDGE[3] =
		asm_data->u.x87.TTMP_FUDGE[4] = asm_data->u.x87.TTMP_FUDGE[5] =
		asm_data->u.x87.TTMP_FUDGE[6] = asm_data->u.x87.TTMP_FUDGE[7] =
		asm_data->u.x87.TTMP_FUDGE[9] = asm_data->u.x87.TTMP_FUDGE[11] =
		asm_data->u.x87.TTMP_FUDGE[13] = asm_data->u.x87.TTMP_FUDGE[15] =
		asm_data->u.x87.TTMP_FUDGE[16] = asm_data->u.x87.TTMP_FUDGE[18] =
		asm_data->u.x87.TTMP_FUDGE[20] = asm_data->u.x87.TTMP_FUDGE[22] = 1.0;
		asm_data->u.x87.TTMP_FUDGE[8] = asm_data->u.x87.TTMP_FUDGE[10] =
		asm_data->u.x87.TTMP_FUDGE[12] = asm_data->u.x87.TTMP_FUDGE[14] =
		asm_data->u.x87.TTMP_FUDGE[17] = asm_data->u.x87.TTMP_FUDGE[19] =
		asm_data->u.x87.TTMP_FUDGE[21] = asm_data->u.x87.TTMP_FUDGE[23] =
		asm_data->u.x87.TTMP_FUDGE[24] = asm_data->u.x87.TTMP_FUDGE[25] =
		asm_data->u.x87.TTMP_FUDGE[26] = asm_data->u.x87.TTMP_FUDGE[27] =
		asm_data->u.x87.TTMP_FUDGE[28] = asm_data->u.x87.TTMP_FUDGE[29] =
		asm_data->u.x87.TTMP_FUDGE[30] = asm_data->u.x87.TTMP_FUDGE[31] = (double) b;

/* Compute the x87 (64-bit) rounding constants */

		asm_data->u.x87.BIGVAL = (float) (3.0 * pow (2.0, 62.0));
		asm_data->u.x87.BIGBIGVAL = (float) (big_word * asm_data->u.x87.BIGVAL);

/* Negate c and store as a double */

		asm_data->u.x87.MINUS_C = (double) -c;

/* Compute constant that converts fft_weight_over_fftlen found in the */
/* two-to-minus-phi tables into the true fft_weight value.  This is usually */
/* FFTLEN / 2, but when doing a non-zero-padded FFT this is FFTLEN / 2k. */

		asm_data->u.x87.NORM012_FF =
			(gwdata->ZERO_PADDED_FFT) ?
				(double) (gwdata->FFTLEN / 2) :
				(double) (gwdata->FFTLEN / 2) / k;

/* Split k for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k if k * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		asm_data->u.x87.K_HI = floor (k / big_word) * big_word;
		asm_data->u.x87.K_LO = k - asm_data->u.x87.K_HI;
		asm_data->u.x87.K_HI_2 = floor (k / big_word / big_word) * big_word * big_word;
		asm_data->u.x87.K_HI_1 = asm_data->u.x87.K_HI - asm_data->u.x87.K_HI_2;
		gwsetmulbyconst (gwdata, gwdata->maxmulbyconst);

/* Compute the normalization constants indexed by biglit array entries */

		asm_data->u.x87.LIMIT_INVERSE[0] =
		asm_data->u.x87.LIMIT_INVERSE[1] =
		asm_data->u.x87.LIMIT_INVERSE[3] =
		asm_data->u.x87.LIMIT_INVERSE[4] = 1.0 / small_word;	/* Lower limit inverse */
		asm_data->u.x87.LIMIT_INVERSE[2] =
		asm_data->u.x87.LIMIT_INVERSE[5] =
		asm_data->u.x87.LIMIT_INVERSE[6] =
		asm_data->u.x87.LIMIT_INVERSE[7] = 1.0 / big_word;	/* Upper limit inverse */

		/* Compute lower limit bigmax */
		temp = small_word * asm_data->u.x87.BIGVAL - asm_data->u.x87.BIGVAL;
		asm_data->u.x87.LIMIT_BIGMAX[0] =
		asm_data->u.x87.LIMIT_BIGMAX[1] =
		asm_data->u.x87.LIMIT_BIGMAX[3] =
		asm_data->u.x87.LIMIT_BIGMAX[4] = temp;

		asm_data->u.x87.LIMIT_BIGMAX_NEG[0] =
		asm_data->u.x87.LIMIT_BIGMAX_NEG[1] =
		asm_data->u.x87.LIMIT_BIGMAX_NEG[3] =
		asm_data->u.x87.LIMIT_BIGMAX_NEG[4] = -temp;	/* Negative lower limit bigmax */

		/* Compute upper limit bigmax */
		temp = big_word * asm_data->u.x87.BIGVAL - asm_data->u.x87.BIGVAL;
		asm_data->u.x87.LIMIT_BIGMAX[2] =
		asm_data->u.x87.LIMIT_BIGMAX[5] =
		asm_data->u.x87.LIMIT_BIGMAX[6] =
		asm_data->u.x87.LIMIT_BIGMAX[7] = temp;

		asm_data->u.x87.LIMIT_BIGMAX_NEG[2] =
		asm_data->u.x87.LIMIT_BIGMAX_NEG[5] =
		asm_data->u.x87.LIMIT_BIGMAX_NEG[6] =
		asm_data->u.x87.LIMIT_BIGMAX_NEG[7] = -temp;	/* Negative upper limit bigmax */

		memcpy (asm_data->u.x87.LIMIT_INVERSE+8, asm_data->u.x87.LIMIT_INVERSE, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_INVERSE+16, asm_data->u.x87.LIMIT_INVERSE, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_INVERSE+24, asm_data->u.x87.LIMIT_INVERSE, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_BIGMAX+8, asm_data->u.x87.LIMIT_BIGMAX, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_BIGMAX+16, asm_data->u.x87.LIMIT_BIGMAX, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_BIGMAX+24, asm_data->u.x87.LIMIT_BIGMAX, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_BIGMAX_NEG+8, asm_data->u.x87.LIMIT_BIGMAX_NEG, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_BIGMAX_NEG+16, asm_data->u.x87.LIMIT_BIGMAX_NEG, 8 * sizeof (double));
		memcpy (asm_data->u.x87.LIMIT_BIGMAX_NEG+24, asm_data->u.x87.LIMIT_BIGMAX_NEG, 8 * sizeof (double));

/* Copy computed assembly constants */

		asm_data->u.x87.P951 = asm_values[0];
		asm_data->u.x87.P618 = asm_values[1];
		asm_data->u.x87.P309 = asm_values[2];
		asm_data->u.x87.M262 = asm_values[3];
		asm_data->u.x87.P588 = asm_values[4];
		asm_data->u.x87.M162 = asm_values[5];
		asm_data->u.x87.M809 = asm_values[6];
		asm_data->u.x87.M382 = asm_values[7];
		asm_data->u.x87.P866 = asm_values[8];
		asm_data->u.x87.P433 = asm_values[9];
		asm_data->u.x87.P577 = asm_values[10];
		asm_data->u.x87.P975 = asm_values[11];
		asm_data->u.x87.P445 = asm_values[12];
		asm_data->u.x87.P180 = asm_values[13];
		asm_data->u.x87.P623 = asm_values[14];
		asm_data->u.x87.M358 = asm_values[15];
		asm_data->u.x87.P404 = asm_values[16];
		asm_data->u.x87.M223 = asm_values[17];
		asm_data->u.x87.M901 = asm_values[18];
		asm_data->u.x87.M691 = asm_values[19];
	}

/* More AVX-512 initialization code */

	if (gwdata->cpu_flags & CPU_AVX512F) {

		/* One pass AVX-512 FFTs use the same 4KB padding scheme as two-pass FFTs */
		if (gwdata->PASS1_SIZE == 0) {
			asm_data->fourKBgapsize = gwdata->FOURKBGAPSIZE;
			asm_data->normblkdst = 0;		/* Small pass 1's with PFA have no padding every 8 clmblkdsts */
			asm_data->normblkdst4 = gwdata->FOURKBGAPSIZE;
			asm_data->normblkdst8 = 0;
			asm_data->pass1blkdst = 8*128;		/* Allows one-pass FFTs to use two-pass carry propagate code */
			asm_data->pass2blkdst = 128;		/* Allows one-pass FFTs to use two-pass carry propagate code */
			asm_data->normval4 = 1;			/* Allows one-pass FFTs to use two-pass add/sub/addsub code */
			asm_data->normval2 = 1;			/* Allows one-pass FFTs to use two-pass add/sub/addsub code */
			asm_data->normval1 = 1;			/* Allows one-pass FFTs to use two-pass add/sub/addsub quick code */
			asm_data->pass2gapsize = 0;		/* Allows one-pass FFTs to use two-pass add/sub/addsub quick code */
		}

		/* Calculate padding constants for two-pass FFTs */
		else {
			int	pass1_size;

			/* There are 64 pad bytes every 4KB and optional pad bytes between each pass 1 block. */
			asm_data->fourKBgapsize = gwdata->FOURKBGAPSIZE;
			asm_data->pass2gapsize = gwdata->PASS2GAPSIZE;
			asm_data->pass2blkdst = gwdata->PASS2_SIZE * 16;
			asm_data->pass2blkdst += (asm_data->pass2blkdst >> 12) * asm_data->fourKBgapsize + asm_data->pass2gapsize;
			asm_data->pass1blkdst = asm_data->pass2blkdst * 8;

			/* Calculate normblk distances used in normalizing two-pass FFTs */
			pass1_size = gwdata->PASS1_SIZE;
			if (pass1_size == 192 || pass1_size == 640 || pass1_size == 768 || pass1_size == 896 || pass1_size == 960 ||
			    pass1_size == 1152 || pass1_size == 1280 || pass1_size == 1344 || pass1_size == 1536 || pass1_size == 1920 ||
			    pass1_size == 2304) {
				asm_data->normblkdst = 64;		/* Many pass 1's do not support padding every 8 clmblkdsts */
				asm_data->normblkdst8 = 0;
			} else if (gwdata->PASS1_CACHE_LINES == 8) {	/* If clm=1 */
				asm_data->normblkdst = 0;		/* Pad for clmblkdst is 64 bytes */
				asm_data->normblkdst8 = 128;		/* Pad for clmblkdst8 is 64 bytes */
			} else if (gwdata->PASS1_CACHE_LINES == 16) {	/* If clm=2 */
				asm_data->normblkdst = 0;		/* Pad for clmblkdst is 0 bytes */
				asm_data->normblkdst8 = 192;		/* Pad for clmblkdst8 is 192 bytes */
			} else {					/* If clm=4 or more */
				asm_data->normblkdst = 64;		/* Pad for clmblkdst is 64 bytes */
				asm_data->normblkdst8 = -64;		/* Pad for clmblkdst8 is -64 bytes */
			}
			asm_data->normblkdst4 = 0;

			/* If we are padding every 4KB, then calculate number of 4KB pages in a pass 2 block */
			/* and "FFT loword/hiword pairs in a 4KB page. */
			/* Otherwise calculate "FFT loword/hiword pairs in a pass 2 block" */
			if (gwdata->FOURKBGAPSIZE) {
				asm_data->normval4 = (gwdata->PASS2_SIZE * 16) / 4096;
				asm_data->normval1 = 4096 / 128;
			} else {
				asm_data->normval4 = 1;
				asm_data->normval1 = gwdata->PASS2_SIZE * 16 / 128;
			}
			/* When adjusting big/lit ptr in add/sub/addsub/smallmul routines, we instead need to know */
			/* "number of "clms in a 4KB page".  We also need to calculate the amount to bump the big/lit pointer. */
			asm_data->normval2 = asm_data->normval1 / (asm_data->cache_line_multiplier / 8);
			if (gwdata->ZERO_PADDED_FFT) asm_data->normval3 = gwdata->pass1_var_data_size -  asm_data->cache_line_multiplier / 2;
			else asm_data->normval3 = gwdata->pass1_var_data_size - asm_data->cache_line_multiplier;
		}
	}

/* More AVX initialization code */

	else if (gwdata->cpu_flags & CPU_AVX) {

		if (gwdata->PASS2_SIZE) {
			int	pass1_size;

			/* There are 64 pad bytes every 4KB and optional */
			/* pad bytes between each pass 1 block. */
			asm_data->fourKBgapsize = gwdata->FOURKBGAPSIZE;
			asm_data->pass2gapsize = gwdata->PASS2GAPSIZE;
			asm_data->pass2blkdst = gwdata->PASS2_SIZE * 16;
			asm_data->pass2blkdst += (asm_data->pass2blkdst >> 12) * asm_data->fourKBgapsize + asm_data->pass2gapsize;
			asm_data->pass1blkdst = asm_data->pass2blkdst * 4;

			/* Calculate normblk distances used in normalizing two-pass FFTs */
			pass1_size = gwdata->PASS1_SIZE;
			if ((pass1_size == 384 && gwdata->ALL_COMPLEX_FFT) ||
			    (pass1_size == 640 && gwdata->ALL_COMPLEX_FFT) ||
			    (pass1_size == 1536 && gwdata->ALL_COMPLEX_FFT)) {
				asm_data->normblkdst = 64;		/* Small pass 1's with PFA have no padding every 8 clmblkdsts */
				asm_data->normblkdst8 = 0;
			} else if (gwdata->PASS1_CACHE_LINES == 4) {	/* If clm=1 */
				asm_data->normblkdst = 0;		/* Pad for clmblkdst is 64 bytes */
				asm_data->normblkdst8 = 64;		/* Pad for clmblkdst8 is 64 bytes */
			} else if (gwdata->PASS1_CACHE_LINES == 8) {	/* If clm=2 */
				asm_data->normblkdst = 0;		/* Pad for clmblkdst is 0 bytes */
				asm_data->normblkdst8 = 192;		/* Pad for clmblkdst8 is 192 bytes */
			} else {					/* If clm=4 or more */
				asm_data->normblkdst = 64;		/* Pad for clmblkdst is 64 bytes */
				asm_data->normblkdst8 = -64;		/* Pad for clmblkdst8 is -64 bytes */
			}

			/* If we are padding 4KB pages, then calculate the number of 4KB pages in a pass 2 block */
			/* and "cache lines in a 4KB page" / clm. */
			/* Otherwise calculate "cache lines in a pass 2 block" / clm. */
			if (gwdata->FOURKBGAPSIZE) {
				asm_data->normval4 = (gwdata->PASS2_SIZE * 16) >> 12;
				asm_data->normval1 = (4096 / 64) / (asm_data->cache_line_multiplier / 4);
			} else {
				asm_data->normval4 = 1;
				asm_data->normval1 = (gwdata->PASS2_SIZE * 16 / 64) / (asm_data->cache_line_multiplier / 4);
			}

			/* Big/lit flags ptr fudge factor in add/sub */
			asm_data->normval2 = asm_data->cache_line_multiplier * 2 * (asm_data->addcount1 - 1);

			/* Big/lit flags ptr fudge factor #2 in add/sub */
			asm_data->normval3 = asm_data->cache_line_multiplier * 2;
		}
	}

/* SSE2 initialization code formerly done in assembly language */

	else if (gwdata->cpu_flags & CPU_SSE2) {

		/* Data size = pass2_size complex numbers * 2 (for SSE2) */
		/* Calculate number of 8KB pages */
		asm_data->normval4 = (gwdata->PASS2_SIZE * 16 * 2) >> 13;

		/* Pass 2 block distance includes 128 pad bytes every 8KB */
		asm_data->pass2gapsize = gwdata->PASS2GAPSIZE;
		asm_data->pass2blkdst = gwdata->PASS2_SIZE * 2 * 8 * 2;
		asm_data->pass2blkdst += (asm_data->pass2blkdst >> 13) * 128 + asm_data->pass2gapsize;
		asm_data->pass1blkdst = asm_data->pass2blkdst;

		/* Calculate normblk distances */
		if (gwdata->SCRATCH_SIZE == 0) {
			/* Normblkdst = pass1blkdst - clm*64 */
			asm_data->normblkdst = asm_data->pass1blkdst - asm_data->cache_line_multiplier * 64;
			asm_data->normblkdst8 = 0;
		} else if ((gwdata->FFT_TYPE == FFT_TYPE_RADIX_4 ||
			    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DELAYED ||
			    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) &&
			   (gwdata->PASS1_SIZE == 80 ||
			    gwdata->PASS1_SIZE == 96 ||
			    gwdata->PASS1_SIZE == 112 ||
			    gwdata->PASS1_SIZE == 224)) {
			/* Small pass 1's with PFA have no padding every 8 clmblkdsts */
			asm_data->normblkdst = 0;
			asm_data->normblkdst8 = 0;
		} else {
			/* Pad in clmblkdst is zero, clmblkdst8 is 128 */
			asm_data->normblkdst = 0;
			asm_data->normblkdst8 = 128;
		}

		if (gwdata->PASS2_SIZE) {
			/* Calculate "cache lines in a chunk" / clm */
			asm_data->normval1 = (8192/64) / asm_data->cache_line_multiplier;

			/* Big/lit flags ptr fudge factor in add/sub */
			asm_data->normval2 = asm_data->cache_line_multiplier * 4 * (asm_data->addcount1 - 1);

			/* Big/lit flags ptr fudge factor #2 in add/sub */
			asm_data->normval3 = asm_data->cache_line_multiplier * 4 -
					     (asm_data->cache_line_multiplier * 4 +
					      asm_data->normval2) *
					     asm_data->normval1 * asm_data->normval4;

			/* This FFT type uses only half the big/lit data */
			if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
				asm_data->normval2 /= 2;
				asm_data->normval3 /= 2;
			}
		}
	}

/* x87 initialization. Calculate constants used in two pass FFTs. */
/* Foremost is the pass 1 blkdst and normalize blkdst for auxiliary mult */
/* routines.  The two values are the same except for larger FFTs which */
/* use a scratch area. */

	else {
		if (gwdata->PASS2_SIZE) {
			/* Pass 2 data size: pass2_size complex numbers */
			/* Compute number of 4KB pages, used in normalized */
			/* add/sub */
			asm_data->normval4 = (gwdata->PASS2_SIZE * 16) >> 12;

			/* pad 64 bytes every 4KB + 64 bytes between blocks */
			asm_data->pass2blkdst =
		        asm_data->pass1blkdst =
			asm_data->normblkdst = asm_data->normval4 * (4096+64) + 64;
			asm_data->normblkdst8 = 0;

			if (gwdata->SCRATCH_SIZE) {
				/* Compute scratch area normblkdst */
				asm_data->normblkdst = asm_data->cache_line_multiplier * 32;

				/* Only larger pass1's have padding */
				if (asm_data->addcount1 >= 128) asm_data->normblkdst8 = 64;
			}

			/* Compute "cache lines in a page" / clm */
			asm_data->normval1 = (4096/32) / asm_data->cache_line_multiplier;

			/* Big/lit flags ptr fudge factor in add/sub */
			asm_data->normval2 = asm_data->cache_line_multiplier * 2 * (asm_data->addcount1 - 1);

			/* Big/lit flags ptr fudge factor #2 in add/sub */
			asm_data->normval3 = asm_data->cache_line_multiplier * 2 -
					     ((asm_data->cache_line_multiplier * 2 +
					       asm_data->normval2) *
					      asm_data->normval1 * asm_data->normval4);
		}
	}

/* For x87 and SSE2 FFTs:  if the carry must be spread over more than 2 words, then set */
/* variable so that assembly code knows this.  We check to see if three small FFT words */
/* can absorb the expected number of bits in a result word.  We are not */
/* aggressive in pushing this limit (we assume no big words will absorb */
/* any of the carry) as it is not a major performance penalty to do 4 or 6 word */
/* carry propagations.  In fact, we might should do 4 or 6 words all the time. */

	if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))) {
		if (gwdata->ZERO_PADDED_FFT ||
		    3.0 * gwdata->NUM_B_PER_SMALL_WORD * log2 (b) >
				2.0 * ((gwdata->NUM_B_PER_SMALL_WORD + 1) * log2 (b) - 1) +
				0.6 * log2 (gwdata->FFTLEN) + log2 (k) + 1.7 * log2 (labs (c)))
			asm_data->SPREAD_CARRY_OVER_EXTRA_WORDS = FALSE;
		else
			asm_data->SPREAD_CARRY_OVER_EXTRA_WORDS = TRUE;
	}

/* Set some global variables that make life easier in the assembly code */
/* that wraps carry out of top FFT word into the bottom FFT word. */
/* This is needed when k > 1 and we are not doing a zero padded FFT. */

	asm_data->TOP_CARRY_NEEDS_ADJUSTING = (k > 1.0 && !gwdata->ZERO_PADDED_FFT);
	if (asm_data->TOP_CARRY_NEEDS_ADJUSTING) {
		unsigned long num_b_in_top_word, num_b_in_second_top_word, num_b_in_third_top_word, num_b_in_k;
		double	carry_adjust_1_mod_k;

/* Copy k and inverted k */

		asm_data->K = k;
		asm_data->INVERSE_K = 1.0 / k;

/* Calculate top carry adjusting constants */

		num_b_in_top_word = gwdata->NUM_B_PER_SMALL_WORD;
		if (is_big_word (gwdata, gwdata->FFTLEN-1)) num_b_in_top_word++;
		num_b_in_second_top_word = gwdata->NUM_B_PER_SMALL_WORD;
		if (is_big_word (gwdata, gwdata->FFTLEN-2)) num_b_in_second_top_word++;
		num_b_in_third_top_word = gwdata->NUM_B_PER_SMALL_WORD;
		if (is_big_word (gwdata, gwdata->FFTLEN-3)) num_b_in_third_top_word++;

		num_b_in_k = (unsigned long) ceil (logb (k));
		asm_data->CARRY_ADJUST1 = pow ((double) b, num_b_in_k);
		carry_adjust_1_mod_k = fltmod (asm_data->CARRY_ADJUST1, k);
		asm_data->CARRY_ADJUST1_HI = floor (carry_adjust_1_mod_k / 131072.0);
		asm_data->CARRY_ADJUST1_LO = carry_adjust_1_mod_k - asm_data->CARRY_ADJUST1_HI * 131072.0;
		asm_data->TWO_TO_17 = 131072.0;
		asm_data->CARRY_ADJUST2 = pow ((double) b, num_b_in_top_word) / asm_data->CARRY_ADJUST1;
		asm_data->CARRY_ADJUST4 = pow ((double) b, num_b_in_second_top_word);
		asm_data->CARRY_ADJUST6 = pow ((double) b, num_b_in_third_top_word);
		if (! (gwdata->cpu_flags & CPU_AVX512F)) {
			if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
				asm_data->CARRY_ADJUST3 = gwfft_partial_weight (gwdata->dd_data, gwdata->FFTLEN-1, dwpn_col (gwdata, gwdata->FFTLEN-1));
				asm_data->CARRY_ADJUST5 = gwfft_partial_weight (gwdata->dd_data, gwdata->FFTLEN-2, dwpn_col (gwdata, gwdata->FFTLEN-2));
				asm_data->CARRY_ADJUST7 = gwfft_partial_weight (gwdata->dd_data, gwdata->FFTLEN-3, dwpn_col (gwdata, gwdata->FFTLEN-3));
			} else {
				asm_data->CARRY_ADJUST3 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-1);
				asm_data->CARRY_ADJUST5 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-2);
				asm_data->CARRY_ADJUST7 = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN-3);
			}
		}

		// It's not worth our time to upgrade the old x87 code to match the AVX-512/AVX/SSE2 code.
		// So, for x87 generate the same constants used prior to version 25.11
		if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2))) {
			unsigned long kbits, kbits_lo;
			kbits = (unsigned long) ceil (log2 (k));
			kbits_lo = kbits / 2;
			asm_data->u.x87.ALT_K_HI = ((unsigned long) k) & ~((1 << kbits_lo) - 1);
			asm_data->u.x87.ALT_K_LO = ((unsigned long) k) &  ((1 << kbits_lo) - 1);
			asm_data->CARRY_ADJUST1_HI = ((unsigned long) asm_data->CARRY_ADJUST1) & ~((1 << kbits_lo) - 1);
			asm_data->CARRY_ADJUST1_LO = ((unsigned long) asm_data->CARRY_ADJUST1) &  ((1 << kbits_lo) - 1);
			if (gwdata->PASS2_SIZE)
				asm_data->CARRY_ADJUST4 *= asm_data->CARRY_ADJUST5;
			else
				asm_data->CARRY_ADJUST6 *= asm_data->CARRY_ADJUST7;
		}

/* In two-pass FFTs, we only support tweaking the top two words. In one-pass FFTs, */
/* we adjust the top three words.  Make sure this works.  A test case that fails: */
/* 489539*3^72778+1.  We should consider supporting tweaking the top 3 words. */

		ASSERTG ((gwdata->PASS2_SIZE &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word) ||
			 (!gwdata->PASS2_SIZE &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word + num_b_in_third_top_word));
		if (!	((gwdata->PASS2_SIZE &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word) ||
			 (!gwdata->PASS2_SIZE &&
			  num_b_in_k <= num_b_in_top_word + num_b_in_second_top_word + num_b_in_third_top_word))) return (GWERROR_INTERNAL + 8);

/* Get the addr of the top three words.  This is funky because in two-pass */
/* FFTs we want the scratch area offset when normalizing after a multiply, */
/* but the FFT data when normalizing after an add/sub.  For one-pass FFTs, */
/* we always want the FFT data offset. */

		asm_data->HIGH_WORD1_OFFSET = addr_offset (gwdata, gwdata->FFTLEN-1);
		asm_data->HIGH_WORD2_OFFSET = addr_offset (gwdata, gwdata->FFTLEN-2);
		asm_data->HIGH_WORD3_OFFSET = addr_offset (gwdata, gwdata->FFTLEN-3);

		raw_gwsetaddin (gwdata, gwdata->FFTLEN-1, 0.0);
		asm_data->HIGH_SCRATCH1_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN-2, 0.0);
		asm_data->HIGH_SCRATCH2_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN-3, 0.0);
		asm_data->HIGH_SCRATCH3_OFFSET = asm_data->ADDIN_OFFSET;
	}

/* Set some global variables that make life easier in the assembly code */
/* that handles zero padded FFTs. */

	if (gwdata->ZERO_PADDED_FFT) {
		unsigned long num_b_in_k, num_b_in_word_0, num_b_in_word_1;
		unsigned long num_b_in_word_2, num_b_in_word_3, num_b_in_word_4, num_b_in_word_5;

		if (gwdata->cpu_flags & CPU_SSE2 && ! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))) {
			asm_data->u.xmm.ZPAD_WORD5_OFFSET = addr_offset (gwdata, 4);
			if (asm_data->u.xmm.ZPAD_WORD5_OFFSET == 8) {  /* FFTLEN = 80 and 112 */
				asm_data->u.xmm.ZPAD_WORD5_RBP_OFFSET = 8;
			} else {
				asm_data->u.xmm.ZPAD_WORD5_RBP_OFFSET = 256;
			}
		}

		asm_data->HIGH_WORD1_OFFSET = addr_offset (gwdata, gwdata->FFTLEN/2-1);
		asm_data->HIGH_WORD2_OFFSET = addr_offset (gwdata, gwdata->FFTLEN/2-2);
		asm_data->HIGH_WORD3_OFFSET = addr_offset (gwdata, gwdata->FFTLEN/2-3);

		raw_gwsetaddin (gwdata, gwdata->FFTLEN/2-1, 0.0);
		asm_data->HIGH_SCRATCH1_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN/2-2, 0.0);
		asm_data->HIGH_SCRATCH2_OFFSET = asm_data->ADDIN_OFFSET;
		raw_gwsetaddin (gwdata, gwdata->FFTLEN/2-3, 0.0);
		asm_data->HIGH_SCRATCH3_OFFSET = asm_data->ADDIN_OFFSET;

		num_b_in_k = (unsigned long) ceil (logb (k));
		num_b_in_word_0 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 0)) num_b_in_word_0++;
		num_b_in_word_1 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 1)) num_b_in_word_1++;
		num_b_in_word_2 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 2)) num_b_in_word_2++;
		num_b_in_word_3 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 3)) num_b_in_word_3++;
		num_b_in_word_4 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 4)) num_b_in_word_4++;
		num_b_in_word_5 = gwdata->NUM_B_PER_SMALL_WORD; if (is_big_word (gwdata, 5)) num_b_in_word_5++;

		asm_data->ZPAD_SHIFT1 = pow ((double) b, (int) num_b_in_word_0);
		asm_data->ZPAD_SHIFT2 = pow ((double) b, (int) num_b_in_word_1);
		asm_data->ZPAD_SHIFT3 = pow ((double) b, (int) num_b_in_word_2);
		asm_data->ZPAD_SHIFT4 = pow ((double) b, (int) num_b_in_word_3);
		asm_data->ZPAD_SHIFT5 = pow ((double) b, (int) num_b_in_word_4);
		asm_data->ZPAD_SHIFT6 = pow ((double) b, (int) num_b_in_word_5);

		if (num_b_in_k <= num_b_in_word_0) asm_data->ZPAD_TYPE = 1;
		else if (num_b_in_k <= num_b_in_word_0 + num_b_in_word_1) asm_data->ZPAD_TYPE = 2;
		else asm_data->ZPAD_TYPE = 3;

		if (asm_data->ZPAD_TYPE == 1) {
			asm_data->ZPAD_K1_LO = k;
			asm_data->ZPAD_INVERSE_K1 = 1.0 / k;
		}

		if (asm_data->ZPAD_TYPE == 2) {
			asm_data->ZPAD_K1_HI = floor (k / asm_data->ZPAD_SHIFT1);
			asm_data->ZPAD_K1_LO = k - asm_data->ZPAD_K1_HI * asm_data->ZPAD_SHIFT1;
			asm_data->ZPAD_INVERSE_K1 = asm_data->ZPAD_SHIFT1 / k;
			asm_data->ZPAD_K2_HI = floor (k / asm_data->ZPAD_SHIFT2);
			asm_data->ZPAD_K2_LO = k - asm_data->ZPAD_K2_HI * asm_data->ZPAD_SHIFT2;
			asm_data->ZPAD_INVERSE_K2 = asm_data->ZPAD_SHIFT2 / k;
			asm_data->ZPAD_K3_HI = floor (k / asm_data->ZPAD_SHIFT3);
			asm_data->ZPAD_K3_LO = k - asm_data->ZPAD_K3_HI * asm_data->ZPAD_SHIFT3;
			asm_data->ZPAD_INVERSE_K3 = asm_data->ZPAD_SHIFT3 / k;
			asm_data->ZPAD_K4_HI = floor (k / asm_data->ZPAD_SHIFT4);
			asm_data->ZPAD_K4_LO = k - asm_data->ZPAD_K4_HI * asm_data->ZPAD_SHIFT4;
			asm_data->ZPAD_INVERSE_K4 = asm_data->ZPAD_SHIFT4 / k;
			asm_data->ZPAD_K5_HI = floor (k / asm_data->ZPAD_SHIFT5);
			asm_data->ZPAD_K5_LO = k - asm_data->ZPAD_K5_HI * asm_data->ZPAD_SHIFT5;
			asm_data->ZPAD_INVERSE_K5 = asm_data->ZPAD_SHIFT5 / k;
			asm_data->ZPAD_K6_HI = floor (k / asm_data->ZPAD_SHIFT6);
			asm_data->ZPAD_K6_LO = k - asm_data->ZPAD_K6_HI * asm_data->ZPAD_SHIFT6;
			asm_data->ZPAD_INVERSE_K6 = asm_data->ZPAD_SHIFT6 / k;
		}

		if (asm_data->ZPAD_TYPE == 3) {
			double	powb, bigpowb;
			powb = pow ((double) b, (int) num_b_in_word_0);
			bigpowb = pow ((double) b, (int) (num_b_in_word_0 + num_b_in_word_1));
			asm_data->ZPAD_K2_HI = floor (k / bigpowb);
			asm_data->ZPAD_K2_MID = floor ((k - asm_data->ZPAD_K2_HI*bigpowb) / powb);
			asm_data->ZPAD_K2_LO = k - asm_data->ZPAD_K2_HI*bigpowb - asm_data->ZPAD_K2_MID*powb;
			asm_data->ZPAD_INVERSE_K2 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_1);
			bigpowb = pow ((double) b, (int) (num_b_in_word_1 + num_b_in_word_2));
			asm_data->ZPAD_K3_HI = floor (k / bigpowb);
			asm_data->ZPAD_K3_MID = floor ((k - asm_data->ZPAD_K3_HI*bigpowb) / powb);
			asm_data->ZPAD_K3_LO = k - asm_data->ZPAD_K3_HI*bigpowb - asm_data->ZPAD_K3_MID*powb;
			asm_data->ZPAD_INVERSE_K3 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_2);
			bigpowb = pow ((double) b, (int) (num_b_in_word_2 + num_b_in_word_3));
			asm_data->ZPAD_K4_HI = floor (k / bigpowb);
			asm_data->ZPAD_K4_MID = floor ((k - asm_data->ZPAD_K4_HI*bigpowb) / powb);
			asm_data->ZPAD_K4_LO = k - asm_data->ZPAD_K4_HI*bigpowb - asm_data->ZPAD_K4_MID*powb;
			asm_data->ZPAD_INVERSE_K4 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_3);
			bigpowb = pow ((double) b, (int) (num_b_in_word_3 + num_b_in_word_4));
			asm_data->ZPAD_K5_HI = floor (k / bigpowb);
			asm_data->ZPAD_K5_MID = floor ((k - asm_data->ZPAD_K5_HI*bigpowb) / powb);
			asm_data->ZPAD_K5_LO = k - asm_data->ZPAD_K5_HI*bigpowb - asm_data->ZPAD_K5_MID*powb;
			asm_data->ZPAD_INVERSE_K5 = powb / k;
			powb = pow ((double) b, (int) num_b_in_word_4);
			bigpowb = pow ((double) b, (int) (num_b_in_word_4 + num_b_in_word_5));
			asm_data->ZPAD_K6_HI = floor (k / bigpowb);
			asm_data->ZPAD_K6_MID = floor ((k - asm_data->ZPAD_K6_HI*bigpowb) / powb);
			asm_data->ZPAD_K6_LO = k - asm_data->ZPAD_K6_HI*bigpowb - asm_data->ZPAD_K6_MID*powb;
			asm_data->ZPAD_INVERSE_K6 = bigpowb / k;
		}

/* Pre-compute the adjustments to copying the 7 words around the halfway point */
/* In a radix-4 delay with full or partial normalization, we must apply an adjustment */
/* so that the copied words are fully weighted. */

		if (gwdata->cpu_flags & CPU_AVX512F) {
			int	i;
			for (i = 0; i < 7; i++) {
				/* Weights for words we are copying before an FFT is performed */
				gwdata->ZPAD_COPY7_ADJUST[i] = gwfft_weight (gwdata->dd_data, gwdata->FFTLEN / 2 - 3 + i);
				/* Inverse weights for ZPAD0 - ZPAD6 applied after an inverse FFT is performed */
				gwdata->ZPAD_0_6_ADJUST[i] = gwfft_weight_inverse (gwdata->dd_data, i);
			}
		} else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
			int	i;
			for (i = 0; i < 7; i++) {
				/* Partial weights for words we are copying before an FFT is performed */
				gwdata->ZPAD_COPY7_ADJUST[i] = gwfft_weight (gwdata->dd_data, dwpn_col (gwdata, gwdata->FFTLEN / 2 - 3 + i));
				/* Partial inverse weights for ZPAD0 - ZPAD6 applied after an inverse FFT is performed */
				gwdata->ZPAD_0_6_ADJUST[i] = gwfft_weight_inverse (gwdata->dd_data, dwpn_col (gwdata, i));
			}
		} else {
			int	i;
			for (i = 0; i < 7; i++) gwdata->ZPAD_COPY7_ADJUST[i] = 1.0;
		}
	}

/* Set the procedure pointers from the proc tables */

#ifdef X86_64
	if (gwdata->cpu_flags & CPU_AVX512F) {
		int	index;
		index = 0;
		if (gwdata->ZERO_PADDED_FFT) index += 2;
		if (gwdata->RATIONAL_FFT) index += 1;
		asm_data->u.zmm.ZMM_CARRIES_ROUTINE = avx512_carries_prctab[index]; // Two-pass square/multiply carry propagation routine
		asm_data->u.zmm.ZMM_OP_CARRIES_ROUTINE = avx512_carries_prctab[gwdata->ZERO_PADDED_FFT ? index + 2 : index]; // Two-pass add/sub/smallmul carry propagation routine
		index = 4 + 5 * index;
		gwdata->GWPROCPTRS[1] = avx512_aux_prctab[index];	// Add
		gwdata->GWPROCPTRS[2] = avx512_aux_prctab[0];		// Add quick
		gwdata->GWPROCPTRS[3] = avx512_aux_prctab[index+1];	// Subtract
		gwdata->GWPROCPTRS[4] = avx512_aux_prctab[1];		// Sub quick
		gwdata->GWPROCPTRS[5] = avx512_aux_prctab[index+2];	// AddSub
		gwdata->GWPROCPTRS[6] = avx512_aux_prctab[2];		// AddSub quick
		gwdata->GWPROCPTRS[7] = avx512_aux_prctab[3];		// Copyzero
		gwdata->GWPROCPTRS[8] = avx512_aux_prctab[index+3];	// Add small
		gwdata->GWPROCPTRS[9] = avx512_aux_prctab[index+4];	// Mul small
		gwdata->GWPROCPTRS[norm_routines] = avx512_prctab[avx512_prctab_index (gwdata, 0, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = avx512_prctab[avx512_prctab_index (gwdata, 0, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = avx512_prctab[avx512_prctab_index (gwdata, 0, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = avx512_prctab[avx512_prctab_index (gwdata, 0, 1, 1)];  // Error, mulbyconst
		gwdata->GWPROCPTRS[zerohigh_routines] = avx512_prctab[avx512_prctab_index (gwdata, 1, 0, 0)];  // Zero, no error
		gwdata->GWPROCPTRS[zerohigh_routines+1] = avx512_prctab[avx512_prctab_index (gwdata, 1, 1, 0)];  // Zero, error
	} else
#endif
	if (gwdata->cpu_flags & CPU_AVX) {
		int	index, aux_base;
		index = 0;
		if (b != 2) index += 4;
		if (gwdata->ZERO_PADDED_FFT) index += 2;
		if (gwdata->RATIONAL_FFT) index += 1;
		asm_data->u.ymm.YMM_CARRIES_ROUTINE = avx_carries_prctab[index]; // Two-pass carry propagation routine
		aux_base = (gwdata->PASS2_SIZE == 0) ? 0 : 44;
		index = aux_base + 4 + 5 * index;
		gwdata->GWPROCPTRS[1] = avx_aux_prctab[index];		// Add
		gwdata->GWPROCPTRS[2] = avx_aux_prctab[aux_base];	// Add quick
		gwdata->GWPROCPTRS[3] = avx_aux_prctab[index+1];	// Subtract
		gwdata->GWPROCPTRS[4] = avx_aux_prctab[aux_base+1];	// Sub quick
		gwdata->GWPROCPTRS[5] = avx_aux_prctab[index+2];	// AddSub
		gwdata->GWPROCPTRS[6] = avx_aux_prctab[aux_base+2];	// AddSub quick
		gwdata->GWPROCPTRS[7] = avx_aux_prctab[aux_base+3];	// Copyzero
		gwdata->GWPROCPTRS[8] = avx_aux_prctab[index+3];	// Add small
		gwdata->GWPROCPTRS[9] = avx_aux_prctab[index+4];	// Mul small
		gwdata->GWPROCPTRS[norm_routines] = avx_prctab[avx_prctab_index (gwdata, 0, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = avx_prctab[avx_prctab_index (gwdata, 0, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = avx_prctab[avx_prctab_index (gwdata, 0, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = avx_prctab[avx_prctab_index (gwdata, 0, 1, 1)];  // Error, mulbyconst
		gwdata->GWPROCPTRS[zerohigh_routines] = avx_prctab[avx_prctab_index (gwdata, 1, 0, 0)];  // Zero, no error
		gwdata->GWPROCPTRS[zerohigh_routines+1] = avx_prctab[avx_prctab_index (gwdata, 1, 1, 0)];  // Zero, error
	}
	else if (gwdata->cpu_flags & CPU_SSE2) {
		memcpy (gwdata->GWPROCPTRS+1, &sse2_aux_prctab[gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN ? 18 : gwdata->PASS2_SIZE ? 9 : 0], 9 * sizeof (void *));
		gwdata->GWPROCPTRS[norm_routines] = sse2_prctab[sse2_prctab_index (gwdata, 0, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = sse2_prctab[sse2_prctab_index (gwdata, 0, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = sse2_prctab[sse2_prctab_index (gwdata, 0, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = sse2_prctab[sse2_prctab_index (gwdata, 0, 1, 1)];  // Error, mulbyconst
		gwdata->GWPROCPTRS[zerohigh_routines] = sse2_prctab[sse2_prctab_index (gwdata, 1, 0, 0)];  // Zero, no error
		gwdata->GWPROCPTRS[zerohigh_routines+1] = sse2_prctab[sse2_prctab_index (gwdata, 1, 1, 0)];  // Zero, error
	}
#ifndef X86_64
	else {
		memcpy (gwdata->GWPROCPTRS+1, &x87_aux_prctab[gwdata->PASS2_SIZE ? 9 : 0], 9 * sizeof (void *));
		gwdata->GWPROCPTRS[norm_routines] = x87_prctab[x87_prctab_index (gwdata, 0, 0, 0)];  // No error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+1] = x87_prctab[x87_prctab_index (gwdata, 0, 1, 0)];  // Error, no mulbyconst
		gwdata->GWPROCPTRS[norm_routines+2] = x87_prctab[x87_prctab_index (gwdata, 0, 0, 1)];  // No error, mulbyconst
		gwdata->GWPROCPTRS[norm_routines+3] = x87_prctab[x87_prctab_index (gwdata, 0, 1, 1)];  // Error, mulbyconst
		gwdata->GWPROCPTRS[zerohigh_routines] = x87_prctab[x87_prctab_index (gwdata, 1, 0, 0)];  // Zero, no error
		gwdata->GWPROCPTRS[zerohigh_routines+1] = x87_prctab[x87_prctab_index (gwdata, 1, 1, 0)];  // Zero, error
	}
#endif

/* Default normalization routines and behaviors */

	gwsetnormroutine (gwdata, 0, 0, 0);
	gwstartnextfft (gwdata, 0);
	raw_gwsetaddin (gwdata, 0, 0.0);
	if (gwdata->square_carefully_count == -1) gwset_square_carefully_count (gwdata, -1);

/* Clear globals */

	asm_data->MAXERR = 0.0;
	gwdata->saved_copyz_n = 0;
	gwdata->GWERROR = 0;
	gwdata->GW_RANDOM = NULL;

/* Compute maximum allowable difference for error checking */
/* This error check is disabled for mod B^N+1 arithmetic */
/* and for radix-4 delay with partial normalization FFTs. */

	if (gwdata->ALL_COMPLEX_FFT || gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
		gwdata->MAXDIFF = 1.0E80;

/* We have observed that the difference seems to vary based on the size */
/* the FFT result word.  This is two times the number of bits per double. */
/* Subtract 1 from bits per double because one bit is the sign bit. */
/* Add log2(b) for the FFT weight that range from 1 to b. */
/* Add in a percentage of the log(FFTLEN) to account for carries. */
/* We use a different threshold for AVX-512/AVX/SSE2 which uses 64-bit instead of */
/* 80-bit doubles during the FFT */

	else {
		double bits_per_double, total_bits, loglen;
		bits_per_double = gwdata->avg_num_b_per_word * log2 (b) - 1.0;
		if (!gwdata->RATIONAL_FFT)
			bits_per_double += log2 (b);
		if (!gwdata->ZERO_PADDED_FFT)
			bits_per_double += log2 (-c);
		loglen = log2 (gwdata->FFTLEN);
		loglen *= 0.69;
		total_bits = bits_per_double * 2.0 + loglen * 2.0;
		gwdata->MAXDIFF = pow ((double) 2.0, total_bits -
				((gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2)) ? 49.08 : 49.65));
	}

/* Clear counters, init internal timers */

	gwdata->fft_count = 0.0;
	asm_data->ASM_TIMERS = (uint32_t *) &gwdata->ASM_TIMERS;

/* Default size of gwnum_alloc array is 50 */

	gwdata->gwnum_alloc = NULL;
	gwdata->gwnum_alloc_count = 0;
	gwdata->gwnum_alloc_array_size = 50;
	gwdata->gwnum_free = NULL;
	gwdata->gwnum_free_count = 0;

/* Activate giants / gwnum shared cached memory allocates */

	gwdata->gdata.blksize = gwnum_datasize (gwdata);

/* Compute alignment for allocated data.  Strangely enough assembly */
/* prefetching works best in pass 1 on a P4 if the data is allocated */
/* on an odd cache line.  An optimal 31 of the 32 cache lines on a 4KB */
/* page will be prefetchable.  Page aligned data would only prefetch */
/* 28 of the 32 cache lines. */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		if (gwdata->PASS2_SIZE == 0) {		/* One pass */
			gwdata->GW_ALIGNMENT = 128;	/* Not tested yet */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else {				/* Two passes */
			gwdata->GW_ALIGNMENT = 1024;	/* Clmblkdst (up to 8) */
			gwdata->GW_ALIGNMENT_MOD = 0;	/* Not tested yet */
		}
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		if (gwdata->PASS2_SIZE == 0) {		/* One pass */
			gwdata->GW_ALIGNMENT = 64;	/* Sandy Bridge cache line alignment */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else if (gwdata->SCRATCH_SIZE == 0) {	/* Small two passes */
			gwdata->GW_ALIGNMENT = 4096;	/* Page alignment */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else {				/* Large two passes */
			gwdata->GW_ALIGNMENT = 1024;	/* Clmblkdst (up to 8) */
			gwdata->GW_ALIGNMENT_MOD = 64;	/* + 1 cache line */
		}
	}
	else if (gwdata->cpu_flags & CPU_SSE2) {
		if (gwdata->PASS2_SIZE == 0) {		/* One pass */
			gwdata->GW_ALIGNMENT = 128;	/* P4 cache line alignment */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else if (gwdata->SCRATCH_SIZE == 0) {	/* Small two passes */
			gwdata->GW_ALIGNMENT = 4096;	/* Page alignment */
			gwdata->GW_ALIGNMENT_MOD = 0;
		} else {				/* Large two passes */
			gwdata->GW_ALIGNMENT = 1024;	/* Clmblkdst (up to 8) */
			gwdata->GW_ALIGNMENT_MOD = 128; /* + 1 cache line */
		}
	} else {
		if (gwdata->PASS2_SIZE == 0)		/* One pass */
			gwdata->GW_ALIGNMENT = 128;	/* P4 cache line alignment */
		else				/* Two passes */
			gwdata->GW_ALIGNMENT = 4096;	/* Page alignment */
		gwdata->GW_ALIGNMENT_MOD = 0;
	}

/* If we are going to use multiple threads for multiplications, then do */
/* the required multi-thread initializations.  Someday, we might allow */
/* setting num_threads after gwsetup so we put all the multi-thread */
/* initialization in its own routine. */

	return (multithread_init (gwdata));
}


/* Utility routines to deal with the 7 words near the half-way point */
/* in a zero-padded AVX-512/AVX/SSE2 FFT.  This used to be all done in assembly code, */
/* but I moved it to C code when the multithread code was written. */

/* When doing zero-padded FFTs, the 7 words around the halfway point must be copied (after applying weights) */
/* for later processing.  This macro does that before a typical forward FFT. */

void xcopy_7_words (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	double	*srcarg, *destarg;

	gwdata = asm_data->gwdata;
	destarg = (double *) asm_data->DESTARG;
	srcarg = (double *) ((char *) destarg + asm_data->DIST_TO_FFTSRCARG);
	if (gwdata->cpu_flags & CPU_AVX512F) {
		char	*srcp;
		srcp = (char *) srcarg + 64;
		destarg[-8] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[3];	/* Copy 1st word above halfway point */
		srcp += asm_data->u.zmm.ZMM_SRC_INCR;
		destarg[-7] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
		srcp += asm_data->u.zmm.ZMM_SRC_INCR;
		destarg[-6] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
		srcp += asm_data->u.zmm.ZMM_SRC_INCR;
		destarg[-5] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		char	*srcp;
		srcp = (char *) srcarg + 32;
		destarg[-5] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[3];	/* Copy 1st word above halfway point */
		srcp += asm_data->u.ymm.YMM_SRC_INCR1;
		destarg[-6] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
		srcp += asm_data->u.ymm.YMM_SRC_INCR2;
		destarg[-7] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
		srcp += asm_data->u.ymm.YMM_SRC_INCR3;
		destarg[-8] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	} else {
		destarg[-5] = srcarg[4] * gwdata->ZPAD_COPY7_ADJUST[3];	/* Copy 1st word above halfway point */
		destarg[-6] = srcarg[12] * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
		destarg[-7] = srcarg[20] * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
		destarg[-8] = srcarg[28] * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	}
	destarg[-9] = * (double *)	/* Copy 1st word below halfway point */
		((char *) srcarg + asm_data->HIGH_WORD1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
	destarg[-10] = * (double *)	/* Copy 2nd word below */
		((char *) srcarg + asm_data->HIGH_WORD2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
	destarg[-11] = * (double *)	/* Copy 3rd word below */
		((char *) srcarg + asm_data->HIGH_WORD3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
}

/* When POSTFFT is set, we must copy the 7 words at two different spots. */
/* These two routines copy the four values above the half-way point after */
/* carries have been propagated and copy the three words just below the */
/* half-way point right after the last NORMRTN has been called. */

void xcopy_4_words (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	double	*srcarg;

	gwdata = asm_data->gwdata;
	srcarg = (double *) asm_data->DESTARG;
	if (gwdata->cpu_flags & CPU_AVX512F) {
		char	*srcp;
		srcp = (char *) srcarg + 64;
		srcarg[-8] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[3];	/* Copy 1st word above halfway point */
		srcp += asm_data->u.zmm.ZMM_SRC_INCR;
		srcarg[-7] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
		srcp += asm_data->u.zmm.ZMM_SRC_INCR;
		srcarg[-6] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
		srcp += asm_data->u.zmm.ZMM_SRC_INCR;
		srcarg[-5] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		char	*srcp;
		srcp = (char *) srcarg + 32;
		srcarg[-5] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[3];	/* Copy 1st word above halfway point */
		srcp += asm_data->u.ymm.YMM_SRC_INCR1;
		srcarg[-6] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
		srcp += asm_data->u.ymm.YMM_SRC_INCR2;
		srcarg[-7] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
		srcp += asm_data->u.ymm.YMM_SRC_INCR3;
		srcarg[-8] = * (double *) srcp * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	} else {
		srcarg[-5] = srcarg[4] * gwdata->ZPAD_COPY7_ADJUST[3]; /* Copy 1st word above halfway point */
		srcarg[-6] = srcarg[12] * gwdata->ZPAD_COPY7_ADJUST[4]; /* Copy 2nd word */
		srcarg[-7] = srcarg[20] * gwdata->ZPAD_COPY7_ADJUST[5]; /* Copy 3rd word */
		srcarg[-8] = srcarg[28] * gwdata->ZPAD_COPY7_ADJUST[6]; /* Copy 4th word */
	}
}

void xcopy_3_words (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	double	*srcarg;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* Copy 1st and 2nd words below the halfway point from the last block */

	srcarg = (double *) asm_data->DESTARG;
	if (asm_data->this_block == asm_data->last_pass1_block) {
		if (gwdata->SCRATCH_SIZE) {
			srcarg[-9] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
			srcarg[-10] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
		} else {
			srcarg[-9] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
			srcarg[-10] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
		}
	}

/* Copy 3rd word below the halfway point */

	if (asm_data->this_block <= gwdata->num_pass1_blocks - 3 &&
	    asm_data->this_block + asm_data->cache_line_multiplier > gwdata->num_pass1_blocks - 3) {
		if (gwdata->SCRATCH_SIZE) {
			srcarg[-11] = * (double *)
					((char *) asm_data->scratch_area +
					 asm_data->HIGH_SCRATCH3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
		} else {
			srcarg[-11] = * (double *)
					((char *) srcarg +
					 asm_data->HIGH_WORD3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
		}
	}
}

void xcopy_3_words_after_gwcarries (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;
	double	*srcarg;

/* Get pointer to the gwdata structure */

	gwdata = asm_data->gwdata;

/* Copy three words below the halfway point from the FFT data */

	srcarg = (double *) asm_data->DESTARG;
	srcarg[-9] = * (double *)
			((char *) srcarg + asm_data->HIGH_WORD1_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[2];
	srcarg[-10] = * (double *)
			((char *) srcarg + asm_data->HIGH_WORD2_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[1];
	srcarg[-11] = * (double *)
			((char *) srcarg + asm_data->HIGH_WORD3_OFFSET) * gwdata->ZPAD_COPY7_ADJUST[0];
}

/* Possible states for the (mis-named) gwdata->pass1_state variable */
/* 0 = pass 1 forward fft */
/* 1 = pass 1 inverse fft */
#define	PASS1_STATE_PASS2		999		/* Auxiliary thread is doing pass 2 work */
#define	PASS1_STATE_MULTITHREAD_OP	4000		/* Auxiliary thread is doing add/sub/addsub/smallmul work */

/* Inline routines for the processing blocks in multi-thread code below */

/* Calculate pass 1 block address.  Note that SSE2 block addresses are */
/* computed based on 128 pad bytes after every 128 cache lines.  That is: */
/* DESTARG + (block * 64) + (block >> 7 * 128).  AVX-512/AVX addresses are based */
/* on 64 to 192 pad bytes after every 64 cache lines. */

static __inline void *pass1_data_addr (
	gwhandle *gwdata,
	struct gwasm_data *asm_data,
	unsigned long block)
{
	if (gwdata->cpu_flags & CPU_AVX512F) {
		block >>= 3;
		return ((char *) asm_data->DESTARG + (block << 7) + ((block >> 5) * gwdata->FOURKBGAPSIZE));
	} else if (gwdata->cpu_flags & CPU_AVX) {
		block >>= 2;
		return ((char *) asm_data->DESTARG + (block << 6) + ((block >> 6) * gwdata->FOURKBGAPSIZE));
	} else
		return ((char *) asm_data->DESTARG + (block << 6) + ((block >> 7) << 7));
}

/* Calculate pass 1 sin/cos/premult address (for those FFTs that do not use */
/* the same sin/cos table for every pass 1 group). */

static __inline void *pass1_premult_addr (
	gwhandle *gwdata,
	unsigned long block)
{
	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4 ||
	    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DELAYED ||
	    gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
		return ((char *) gwdata->pass1_var_data + block / gwdata->PASS1_CACHE_LINES * gwdata->pass1_var_data_size);
	return (NULL);
}

/* Calculate pass 1 state 1 normalization ptr addresses. */

static __inline void pass1_state1_norm_addrs (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	if (gwdata->cpu_flags & CPU_AVX512F) {
		// AVX-512 puts biglit flags in the variable sin/cos data
	} else if (gwdata->cpu_flags & CPU_AVX) {
		asm_data->norm_ptr1 = (char *) asm_data->norm_biglit_array + (asm_data->addcount1 * asm_data->this_block * 2);
	} else {
		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
			asm_data->norm_ptr1 = (char *) asm_data->norm_biglit_array + (asm_data->addcount1 * asm_data->this_block * 2);
		else {
			asm_data->norm_ptr1 = (char *) asm_data->norm_biglit_array + (asm_data->addcount1 * asm_data->this_block * 4);
			asm_data->norm_ptr2 = (char *) asm_data->norm_col_mults + (asm_data->this_block * 32);
		}
	}
}

static __inline void pass1_state1_carry_addrs (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	if (gwdata->cpu_flags & CPU_AVX512F) {
		// Calculate address of biglit data in this pass 1 block's variable data.
		asm_data->norm_ptr1 = (char *) gwdata->pass1_var_data +
				      gwdata->pass1_var_data_size * asm_data->this_block / gwdata->PASS1_CACHE_LINES + gwdata->biglit_data_offset;
	} else if (gwdata->cpu_flags & CPU_AVX) {
		asm_data->norm_ptr1 = (char *) asm_data->norm_biglit_array + (asm_data->addcount1 * asm_data->this_block * 2);
	} else {
		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
			asm_data->norm_ptr1 = (char *) asm_data->norm_biglit_array + (asm_data->addcount1 * asm_data->this_block * 2);
		else {
			asm_data->norm_ptr1 = (char *) asm_data->norm_biglit_array + (asm_data->addcount1 * asm_data->this_block * 4);
			asm_data->norm_ptr2 = (char *) asm_data->norm_col_mults + (asm_data->this_block * 32);
		}
	}
}

/* Calculate pass 2 block address */

static __inline void *pass2_data_addr (
	gwhandle *gwdata,
	struct gwasm_data *asm_data,
	unsigned long block)
{
	return ((char *) asm_data->DESTARG + block * (unsigned long) asm_data->pass2blkdst);
}

/* Calculate pass 2 premultiplier address */

static __inline void * pass2_premult_addr (
	gwhandle *gwdata,
	unsigned long block)
{
	return ((char *) gwdata->adjusted_pass2_premults + block * gwdata->pass2_premult_block_size);
}

/* Assign a thread's first block to process in pass 1 state 0.  These are assigned in */
/* sequential order. */

static __inline void pass1_state0_assign_first_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	asm_data->this_block = gwdata->next_block;
	asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
	asm_data->premult_addr = pass1_premult_addr (gwdata, asm_data->this_block);
	gwdata->next_block += asm_data->cache_line_multiplier;
}

/* Assign next available block in pass 1 state 0.  These are assigned in */
/* sequential order. */

static __inline void pass1_state0_assign_next_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	if (gwdata->next_block < gwdata->num_pass1_blocks) {
		asm_data->next_block = gwdata->next_block;
		gwdata->next_block += asm_data->cache_line_multiplier;
		/* Init prefetching for the next block */
		if (gwdata->hyperthread_prefetching) {
			asm_data->data_prefetch = asm_data->data_addr;
			asm_data->premult_prefetch = asm_data->premult_addr;
			gwevent_signal (&asm_data->hyperthread_work_to_do);
		} else {
			asm_data->data_prefetch = pass1_data_addr (gwdata, asm_data, asm_data->next_block);
			asm_data->premult_prefetch = pass1_premult_addr (gwdata, asm_data->next_block);
		}
	} else {
		asm_data->next_block = asm_data->this_block;
		asm_data->data_prefetch = asm_data->data_addr;
		asm_data->premult_prefetch = asm_data->premult_addr;
	}
}

/* Assign a thread's first block in pass 1 state 1. */
/* Returns TRUE if successful, FALSE if no free sections were found */

int pass1_state1_assign_first_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	int	i;

	/* Get thread number / section number */
	i = asm_data->thread_num;

	/* If there are zero blocks in the section, see if we can find (part of) another section we can work on */
	if (gwdata->pass1_carry_sections[i].next_block == gwdata->pass1_carry_sections[i].last_block) {
		unsigned int j, largest_unfinished, largest_unfinished_size, min_split_size;

		/* Calculate minimum split size.  It must be at least num_postfft_blocks. */
		/* It must also be a multiple of 8 if AVX zero-padded (because of using YMM_SRC_INCR[0-7] in */
		/* ynorm012 macros) or if SSE2 because of using BIGLIT_INCR4 in xnorm012 macros. */
		/* AVX-512 will already be a multiple of 8. */
		min_split_size = gwdata->num_postfft_blocks;
		if (((gwdata->cpu_flags & CPU_AVX) && gwdata->ZERO_PADDED_FFT) ||
		    !(gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX)))
			min_split_size = round_up_to_multiple_of (min_split_size, 8);

		/* Find largest unfinished section for splitting */
		largest_unfinished = i;
		largest_unfinished_size = 0;
		for (j = 0; j < gwdata->num_threads; j++) {
			unsigned int size;

			if (gwdata->pass1_carry_sections[j].section_state >= 2) continue;
			/* Since new section must be at least num_postfft_blocks in size and the existing */
			/* section is prefetching the next block, make sure we are splitting at least */
			/* num_postfft_blocks + asm_data->cache_line_multiplier in size. */
			/* Note: deadlocks can occur if we don't take over sections that haven't started. */
			size = gwdata->pass1_carry_sections[j].last_block - gwdata->pass1_carry_sections[j].next_block;
			if (gwdata->pass1_carry_sections[j].section_state != 0 &&
			    size < gwdata->num_postfft_blocks + asm_data->cache_line_multiplier) continue;
			/* See if this section is big enough and the largest one thusfar. */
			if (size > largest_unfinished_size && size >= min_split_size) {
				largest_unfinished = j;
				largest_unfinished_size = size;
			}
		}

		/* If we found no sections that we can split, then return FALSE */
		if (largest_unfinished_size == 0) return (FALSE);

		/* If the target section hasn't even started, then take the entire section under */
		/* the theory that the thread must be busy servicing some other process. */
		if (gwdata->pass1_carry_sections[largest_unfinished].section_state == 0) {
			memcpy (&gwdata->pass1_carry_sections[i],
				&gwdata->pass1_carry_sections[largest_unfinished],
				sizeof (struct pass1_carry_sections));
			gwdata->pass1_carry_sections[largest_unfinished].start_block =
			gwdata->pass1_carry_sections[largest_unfinished].next_block =
				gwdata->pass1_carry_sections[largest_unfinished].last_block;
			gwdata->pass1_carry_sections[largest_unfinished].dependent_section = -1;
		}

		/* Split the target section */
		else {
			unsigned int new_size;

			new_size = largest_unfinished_size / 2;
			new_size = round_down_to_multiple_of (new_size, asm_data->cache_line_multiplier);
			if (((gwdata->cpu_flags & CPU_AVX512F) && gwdata->ZERO_PADDED_FFT) ||
			    (!(gwdata->cpu_flags & CPU_AVX512F) && (gwdata->cpu_flags & CPU_AVX) && gwdata->ZERO_PADDED_FFT) ||
			    (!(gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))))
				new_size = round_up_to_multiple_of (new_size, 8);
			if (new_size < gwdata->num_postfft_blocks) new_size = gwdata->num_postfft_blocks;

			// Make sure the new size is a multiple of num_postfft_blocks.  This is necessary because
			// a multihreaded PRP test of 66666*5^1560000-1 will fail if we don't.  It fails because the
			// we need more than 4 carry words and the FFT has a clm of 1 (which is 4 words).  If the
			// start_block in a section is not a multiple of num_postfft_blocks, then normval2
			// used in ynorm012_wpn will not properly correct the big/lit flags pointer.
			new_size = round_up_to_multiple_of (new_size, gwdata->num_postfft_blocks);

			// Peel off the ending blocks of largest_unfinished section
			gwdata->pass1_carry_sections[i].start_block =
			gwdata->pass1_carry_sections[i].next_block =
				gwdata->pass1_carry_sections[largest_unfinished].last_block - new_size;
			gwdata->pass1_carry_sections[i].last_block =
				gwdata->pass1_carry_sections[largest_unfinished].last_block;
			gwdata->pass1_carry_sections[i].can_carry_into_next =
				gwdata->pass1_carry_sections[largest_unfinished].can_carry_into_next;

			// Shrink largest_unfinished section
			gwdata->pass1_carry_sections[largest_unfinished].last_block -= new_size;
			gwdata->pass1_carry_sections[largest_unfinished].can_carry_into_next = FALSE;
			gwdata->pass1_carry_sections[i].dependent_section = largest_unfinished;
		}

		/* Finally, patch dependent_section if necessary */
		if (! gwdata->pass1_carry_sections[i].can_carry_into_next) {
			for (j = 0; j < gwdata->num_threads; j++) {
				if (gwdata->pass1_carry_sections[j].start_block == gwdata->pass1_carry_sections[j].last_block) continue;
				if (gwdata->pass1_carry_sections[j].start_block == gwdata->pass1_carry_sections[i].last_block ||
				    (gwdata->pass1_carry_sections[j].start_block == 0 &&
				     gwdata->pass1_carry_sections[i].last_block == gwdata->num_pass1_blocks)) {
					gwdata->pass1_carry_sections[j].dependent_section = i;
					break;
				}
			}
			ASSERTG (j != gwdata->num_threads);
		}
	}

	/* Return the first block in this section */
	asm_data->this_block = gwdata->pass1_carry_sections[i].next_block;
	asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
	asm_data->premult_addr = pass1_premult_addr (gwdata, asm_data->this_block);
	pass1_state1_norm_addrs (gwdata, asm_data);

	/* Bump section's next available block pointer, update section state */
	gwdata->pass1_carry_sections[i].next_block += asm_data->cache_line_multiplier;
	gwdata->pass1_carry_sections[i].section_state = 1;

	/* In SSE2 zero-padded FFT, set the upper half of the carries to zero instead */
	/* of XMM_BIGVAL.  This saves us a couple of clocks per FFT data element in xnorm_2d_zpad. */
	/* Ideally, we'd eliminate this code by fixing the add/sub and carry propagation code to also */
	/* expect zero in the upper half of each carry cache line  (like we did in the AVX code). */
	if (gwdata->ZERO_PADDED_FFT && ! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))) {
		int	i, carry_table_size;
		double	*table;
		carry_table_size = (gwdata->PASS1_SIZE) << 1;
		table = (double *) asm_data->carries;
		for (i = 0; i < carry_table_size; i += 8, table += 8) {
			table[4] = 0.0;
			table[5] = 0.0;
			table[6] = 0.0;
			table[7] = 0.0;
		}
	}

	/* Return success */
	return (TRUE);
}

/* Assign a thread's next block in pass 1 state 1. */

void pass1_state1_assign_next_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	int	i;

	/* Get thread number / section number */
	i = asm_data->thread_num;

	/* Return the next block in this section (or first block in next section */
	/* for carry propagation and possible postfft processing) */
	asm_data->next_block = gwdata->pass1_carry_sections[i].next_block;
	if (asm_data->next_block == gwdata->pass1_carry_sections[i].last_block && gwdata->pass1_carry_sections[i].section_state == 3)
		asm_data->next_block = asm_data->this_block;
	else if (asm_data->next_block == gwdata->num_pass1_blocks)
		asm_data->next_block = 0;

	/* Init prefetching for the next block */
	if (gwdata->hyperthread_prefetching) {
		asm_data->data_prefetch = asm_data->data_addr;
		asm_data->premult_prefetch = asm_data->premult_addr;
		gwevent_signal (&asm_data->hyperthread_work_to_do);
	} else {
		asm_data->data_prefetch = pass1_data_addr (gwdata, asm_data, asm_data->next_block);
		asm_data->premult_prefetch = pass1_premult_addr (gwdata, asm_data->next_block);
	}
}

/* Assign a thread's first pass 2 block.  These are assigned in */
/* sequential order. */

static __inline void pass2_assign_first_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	asm_data->this_block = gwdata->next_block++;
	asm_data->data_addr = pass2_data_addr (gwdata, asm_data, asm_data->this_block);
	asm_data->premult_addr = pass2_premult_addr (gwdata, asm_data->this_block);
}

/* Assign next available block in pass 2.  These are assigned in sequential order. */

static __inline void pass2_assign_next_block (
	gwhandle *gwdata,
	struct gwasm_data *asm_data)
{
	if (gwdata->next_block < gwdata->num_pass2_blocks) {
		asm_data->next_block = gwdata->next_block++;
		/* Init prefetching for the next block */
		if (gwdata->hyperthread_prefetching) {
			asm_data->data_prefetch = asm_data->data_addr;
			asm_data->premult_prefetch = asm_data->premult_addr;
			gwevent_signal (&asm_data->hyperthread_work_to_do);
		} else {
			asm_data->data_prefetch = pass2_data_addr (gwdata, asm_data, asm_data->next_block);
			asm_data->premult_prefetch = pass2_premult_addr (gwdata, asm_data->next_block);
		}
	} else {
		asm_data->next_block = asm_data->this_block;
		asm_data->data_prefetch = asm_data->data_addr;
		asm_data->premult_prefetch = asm_data->premult_addr;
	}
}


/* Routine for auxiliary threads */

struct thread_data {
	gwhandle *gwdata;
	int	thread_num;
	void	*asm_data_alloc;
	void	*scratch_area;
	void	*carries;
};

void auxiliary_thread (void *arg)
{
	gwhandle *gwdata;
	struct gwasm_data *asm_data, *main_thread_asm_data;
	struct thread_data *info;

/* Get pointers to various structures */

	info = (struct thread_data *) arg;
	gwdata = info->gwdata;
	main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Each thread needs its own copy of the asm_data.  Each thread needs */
/* its own stack too. */

	asm_data = (struct gwasm_data *) ((char *) info->asm_data_alloc + NEW_STACK_SIZE);
	memcpy (asm_data, main_thread_asm_data, sizeof (struct gwasm_data));

/* Set the thread number so that the assembly code can differentiate between */
/* the main thread and an auxiliary thread. */

	asm_data->thread_num = info->thread_num;

/* Each auxiliary thread needs its own pass 1 scratch area */

	asm_data->scratch_area = info->scratch_area;

/* Init each auxiliary thread's carries area */
	
	asm_data->carries = info->carries;
	if (gwdata->cpu_flags & CPU_AVX512F) {
		int	i, carry_table_size;
		double	*table;
		carry_table_size = gwdata->PASS1_SIZE;
		table = asm_data->carries;
		for (i = 0; i < carry_table_size; i++) *table++ = asm_data->u.zmm.ZMM_RNDVAL;
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		int	i, carry_table_size;
		double	carryval, *table;
		carry_table_size = gwdata->PASS1_SIZE;
		if (gwdata->b == 2) carryval = 3.0 * 131072.0 * 131072.0 * 131072.0;
		else carryval = 0.0;
		table = asm_data->carries;
		for (i = 0; i < carry_table_size; i++)
			if (gwdata->ZERO_PADDED_FFT && (i & 4)) *table++ = 0.0;
			else *table++ = carryval;
	} else {
		int	i, carry_table_size;
		double	xmm_bigval, *table;
		carry_table_size = gwdata->PASS1_SIZE * 2;
		xmm_bigval = 3.0 * 131072.0 * 131072.0 * 131072.0;
		table = asm_data->carries;
		for (i = 0; i < carry_table_size; i++) *table++ = xmm_bigval;
	}

/* Call optional user provided callback routine so that the caller can set */
/* the thread's priority and affinity */

	if (gwdata->thread_callback != NULL)
		(*gwdata->thread_callback) (info->thread_num, 0, gwdata->thread_callback_data);

/* Each auxiliary compute thread optionally creates its own hyperthread for prefetching */

	create_auxiliary_hyperthread (asm_data);

/* This describes the various valid states for synchronizing the main and auxiliary threads.  We use four variables/mutexes. */
/* work_to_do	all_threads_done    catch_stragglers	num_active_threads */
/* FALSE	FALSE			FALSE		0			inital state */
/* TRUE		FALSE			FALSE		0			start aux workers */
/* TRUE		FALSE			FALSE		> 0			one or more aux workers active */
/* FALSE	FALSE			FALSE		> 0			some aux workers active, some finished */
/* FALSE	TRUE			TRUE		0			aux workers finished (but some could still wake up) */

/* Loop waiting for work to do.  The main thread signals the thread_work_to_do event whenever there is work */
/* for the auxiliary thread(s) to do. */

	for ( ; ; ) {
		gwevent_wait (&gwdata->thread_work_to_do, 0);

/* If threads are to exit, break out of this work loop */

		if (gwdata->threads_must_exit) break;

/* WARNING: We can get here with no work to do (and with the main thread moved on thinking all auxiliary threads are done). */
/* For example, consider an 8-thread example.  The main thread signals the 7 auxiliary threads to start.  Six of the threads */
/* can wake up, complete the work, and signal gwdata->all_threads_done when gwdata->num_active_threads returns to zero. */
/* The main thread catches the signal and continues on.  The seventh auxiliary thread can now get scheduled by the OS */
/* and end up here.  The catch_straggler_threads flag handles was created to handle this situation. */

		gwmutex_lock (&gwdata->thread_lock);
		if (gwdata->catch_straggler_threads) {
			gwmutex_unlock (&gwdata->thread_lock);
			continue;
		}

/* Increment the number of active threads so that we can tell when all auxiliary routines have finished. */

		gwdata->num_active_threads++;

/* If we are to do add/sub/addsub/smallmul work, copy a little bit of state and go do the work */

		if (gwdata->pass1_state == PASS1_STATE_MULTITHREAD_OP) {
			asm_data->DBLARG = main_thread_asm_data->DBLARG;
			gwmutex_unlock (&gwdata->thread_lock);
			do_multithread_op_work (gwdata, asm_data);
			goto aux_out_of_work_unlocked;
		}

/* Copy the main thread's asm_data's DESTARG for proper next_block */
/* address calculations.  We'll copy more asm_data later. */

		asm_data->DESTARG = main_thread_asm_data->DESTARG;

/* Get an available block for this thread to process (store it in the */
/* next_block field).  Note we set this_block to a dummy value so that */
/* get_next_block knows there is more work to do.  NOTE: There is no */
/* guarantee that there is an available block to process. */

		asm_data->this_block = 0;
		if (gwdata->pass1_state == 0) {
			if (gwdata->next_block >= gwdata->num_pass1_blocks) goto aux_out_of_work_locked;
			pass1_state0_assign_first_block (gwdata, asm_data);
			pass1_state0_assign_next_block (gwdata, asm_data);
		} else if (gwdata->pass1_state == 1) {
			if (! pass1_state1_assign_first_block (gwdata, asm_data)) goto aux_out_of_work_locked;
			pass1_state1_assign_next_block (gwdata, asm_data);
		} else {
			if (gwdata->next_block >= gwdata->num_pass2_blocks) goto aux_out_of_work_locked;
			pass2_assign_first_block (gwdata, asm_data);
			pass2_assign_next_block (gwdata, asm_data);
		}

/* Copy some data from the main thread's asm_data to this thread.  We only */
/* need to copy data that can change for each multiply call */

		asm_data->DIST_TO_FFTSRCARG = main_thread_asm_data->DIST_TO_FFTSRCARG;
		asm_data->DIST_TO_MULSRCARG = main_thread_asm_data->DIST_TO_MULSRCARG;
		asm_data->ffttype = main_thread_asm_data->ffttype;
		asm_data->thread_work_routine = main_thread_asm_data->thread_work_routine;
		if (gwdata->pass1_state == 1) {
			asm_data->NORMRTN = main_thread_asm_data->NORMRTN;
			asm_data->zero_fft = main_thread_asm_data->zero_fft;
			asm_data->const_fft = main_thread_asm_data->const_fft;
			asm_data->add_sub_smallmul_op = main_thread_asm_data->add_sub_smallmul_op;
			asm_data->ADDIN_ROW = main_thread_asm_data->ADDIN_ROW;
			asm_data->ADDIN_OFFSET = main_thread_asm_data->ADDIN_OFFSET;
			asm_data->ADDIN_VALUE = main_thread_asm_data->ADDIN_VALUE;
			if (gwdata->cpu_flags & CPU_AVX512F) {
				asm_data->u.zmm.ZMM_MULCONST = main_thread_asm_data->u.zmm.ZMM_MULCONST;
				asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE = main_thread_asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE;
				asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_SMALL_BASE = main_thread_asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_SMALL_BASE;
				asm_data->u.zmm.ZMM_K_TIMES_MULCONST_LO = main_thread_asm_data->u.zmm.ZMM_K_TIMES_MULCONST_LO;
				asm_data->u.zmm.ZMM_MINUS_C_TIMES_MULCONST = main_thread_asm_data->u.zmm.ZMM_MINUS_C_TIMES_MULCONST;
			} else if (gwdata->cpu_flags & CPU_AVX) {
				memcpy (asm_data->u.ymm.YMM_MULCONST, main_thread_asm_data->u.ymm.YMM_MULCONST, 4 * sizeof (double));
				memcpy (asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI, main_thread_asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI, 4 * sizeof (double));
				memcpy (asm_data->u.ymm.YMM_K_TIMES_MULCONST_LO, main_thread_asm_data->u.ymm.YMM_K_TIMES_MULCONST_LO, 4 * sizeof (double));
				memcpy (asm_data->u.ymm.YMM_MINUS_C_TIMES_MULCONST, main_thread_asm_data->u.ymm.YMM_MINUS_C_TIMES_MULCONST, 4 * sizeof (double));
			} else if (gwdata->cpu_flags & CPU_SSE2) {
				memcpy (asm_data->u.xmm.XMM_MULCONST, main_thread_asm_data->u.xmm.XMM_MULCONST, 2 * sizeof (double));
				memcpy (asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI, main_thread_asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI, 2 * sizeof (double));
				memcpy (asm_data->u.xmm.XMM_K_TIMES_MULCONST_LO, main_thread_asm_data->u.xmm.XMM_K_TIMES_MULCONST_LO, 2 * sizeof (double));
				memcpy (asm_data->u.xmm.XMM_MINUS_C_TIMES_MULCONST, main_thread_asm_data->u.xmm.XMM_MINUS_C_TIMES_MULCONST, 2 * sizeof (double));
			}
			// Copy maxerr and sumout
			if (gwdata->cpu_flags & CPU_AVX512F) {
				asm_data->MAXERR = main_thread_asm_data->MAXERR;
			} else if (gwdata->cpu_flags & CPU_AVX) {
				memcpy (asm_data->u.ymm.YMM_MAXERR, &main_thread_asm_data->u.ymm.YMM_MAXERR, 4 * sizeof (double));
			} else {
				asm_data->u.xmm.XMM_SUMOUT[0] = 0.0;
				asm_data->u.xmm.XMM_SUMOUT[1] = 0.0;
				memcpy (asm_data->u.xmm.XMM_MAXERR, &main_thread_asm_data->u.xmm.XMM_MAXERR, 2 * sizeof (double));
			}
		}

/* Now call the assembly code to do some work! */

		gwmutex_unlock (&gwdata->thread_lock);
		if (gwdata->pass1_state < PASS1_STATE_PASS2)
			pass1_aux_entry_point (asm_data);
		else
			pass2_aux_entry_point (asm_data);

/* The auxiliary thread has run out of work.  Decrement the count of number of active auxiliary threads. */
/* Signal all threads done when last auxiliary thread is done. */

aux_out_of_work_unlocked:
		gwmutex_lock (&gwdata->thread_lock);
aux_out_of_work_locked:
		gwdata->num_active_threads--;
		if (gwdata->num_active_threads == 0) {
			gwdata->catch_straggler_threads = TRUE;
			gwevent_signal (&gwdata->all_threads_done);
		}

/* Reset thread_work_to_do event before looping to wait for more work. */

		gwevent_reset (&gwdata->thread_work_to_do);
		gwmutex_unlock (&gwdata->thread_lock);
	}

/* Call optional user provided callback routine so that the caller can */
/* do any necessary cleanup. */

	if (gwdata->thread_callback != NULL)
		(*gwdata->thread_callback) (info->thread_num, 1, gwdata->thread_callback_data);

/* Wait for the optional prefetching hyperthread to finish before we delete the shared asm_data structure */

	if (asm_data->hyperthread_id) {
		gwevent_signal (&asm_data->hyperthread_work_to_do);
		gwthread_wait_for_exit (&asm_data->hyperthread_id);
		gwevent_destroy (&asm_data->hyperthread_work_to_do);
	}

/* Free the allocated memory and exit the auxiliary thread */

	aligned_free (info->scratch_area);
	aligned_free (info->carries);
	aligned_free (info->asm_data_alloc);
	free (arg);
}


/* This routine is called by the main thread assembly code to */
/* fire up all the auxiliary worker threads in pass 1. */

void pass1_wake_up_threads (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;

/* Obtain the lock if multithreading. */
/* Why?  The other threads should be off right now, but we've seen in the debugger that it */
/* is possible to get here with the lock held by an auxiliary thread that is finishing up. */

	if (gwdata->num_threads > 1) {
		gwmutex_lock (&gwdata->thread_lock);
	}

/* Init pass1_state (kludge: passed in next_block parameter) */
/* State is either 0 (forward FFT only) or 1 (inverse FFT and */
/* optional next forward FFT) */

	gwdata->pass1_state = asm_data->next_block;

/* Initialization for the forward FFT case */

	if (gwdata->pass1_state == 0) {

/* Prior to doing a forward zero-padded FFT copy 7 FFT elements */

		if (gwdata->ZERO_PADDED_FFT) xcopy_7_words (asm_data);

/* Set up the this_block and next_block values for the main thread */

		gwdata->next_block = 0;
		pass1_state0_assign_first_block (gwdata, asm_data);
		pass1_state0_assign_next_block (gwdata, asm_data);
	}

/* Initialization for the inverse FFT (and optional next forward FFT) case */

	if (gwdata->pass1_state == 1) {
		int	i, num_blocks_allocated, num_blocks_remaining, num_threads_remaining;

/* Initialize values that the normalization / carry propagation code uses. */

		if (gwdata->cpu_flags & CPU_AVX512F) {
		}
		else if (gwdata->cpu_flags & CPU_AVX) {
			asm_data->u.ymm.YMM_MAXERR[0] =
			asm_data->u.ymm.YMM_MAXERR[1] =
			asm_data->u.ymm.YMM_MAXERR[2] =
			asm_data->u.ymm.YMM_MAXERR[3] = asm_data->MAXERR;
		} else {
			((double *) asm_data->DESTARG)[-3] = 0.0; // For accumulating SUMOUT
			asm_data->u.xmm.XMM_SUMOUT[0] =
			asm_data->u.xmm.XMM_SUMOUT[1] = 0.0;
			asm_data->u.xmm.XMM_MAXERR[0] =
			asm_data->u.xmm.XMM_MAXERR[1] = asm_data->MAXERR;
		}

/* Split the available pass 1 blocks into contiguous sections for each thread to process. */

		num_blocks_allocated = 0;
		num_blocks_remaining = gwdata->num_pass1_blocks;
		num_threads_remaining = gwdata->num_threads;
		for (i = 0; i < (int) gwdata->num_threads; i++) {
			int	num_blocks;

			/* Calculate number of blocks in this section.  It must be a multiple of cache_line_multiplier. */
			/* It must also be a multiple of 8 if AVX-512 zero-padded (because of using ZMM_SRC_INCR[0-7] in */
			/* znorm012 macros) or if AVX zero-padded (because of using YMM_SRC_INCR[0-7] in ynorm012 macros) */
			/* or if SSE2 because of using BIGLIT_INCR4 in xnorm012 macros. */
			num_blocks = divide_rounding_up (num_blocks_remaining, num_threads_remaining);
			num_blocks = round_up_to_multiple_of (num_blocks, asm_data->cache_line_multiplier);
			if (((gwdata->cpu_flags & CPU_AVX512F) && gwdata->ZERO_PADDED_FFT) ||
			    (!(gwdata->cpu_flags & CPU_AVX512F) && (gwdata->cpu_flags & CPU_AVX) && gwdata->ZERO_PADDED_FFT) ||
			    (!(gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX))))
				num_blocks = round_up_to_multiple_of (num_blocks, 8);
			if (num_blocks < (int) gwdata->num_postfft_blocks) num_blocks = gwdata->num_postfft_blocks;

			// Make sure the new size is a multiple of num_postfft_blocks.  This is necessary because
			// a 5 or 7 threaded PRP test of 66666*5^1560000-1 will fail if we don't.  It fails because the
			// we need more than 4 carry words and the FFT has a clm of 1 (which is 4 words).  If the
			// start_block in a section is not a multiple of num_postfft_blocks, then normval2
			// used in ynorm012_wpn will not properly correct the big/lit flags pointer.
			num_blocks = round_up_to_multiple_of (num_blocks, gwdata->num_postfft_blocks);

			/* Init section info */
			gwdata->pass1_carry_sections[i].start_block = num_blocks_allocated;
			gwdata->pass1_carry_sections[i].last_block = num_blocks_allocated + num_blocks;
			gwdata->pass1_carry_sections[i].next_block = num_blocks_allocated;
			gwdata->pass1_carry_sections[i].section_state = 0;
			gwdata->pass1_carry_sections[i].can_carry_into_next = FALSE;
			gwdata->pass1_carry_sections[i].dependent_section = i - 1;

			/* Adjust counters */
			num_blocks_allocated += num_blocks;
			num_blocks_remaining -= num_blocks;
			num_threads_remaining--;

			/* The last section that allocates blocks propagates its carries into the first section */
			if (num_blocks && num_blocks_remaining == 0)
				gwdata->pass1_carry_sections[0].dependent_section = i;
		}

/* Set up the this_block and next_block values for the main thread */

		pass1_state1_assign_first_block (gwdata, asm_data);
		pass1_state1_assign_next_block (gwdata, asm_data);
	}

/* Signal the auxiliary threads to resume working */

	if (gwdata->num_threads > 1) {
		gwdata->catch_straggler_threads = FALSE;
		gwevent_reset (&gwdata->all_threads_done);
		gwevent_signal (&gwdata->thread_work_to_do);
		gwmutex_unlock (&gwdata->thread_lock);
	}
}

/* This routine is called before normalizing a block of data. */
/* This gives us an opportunity to make some adjustments for zero-padded FFTs. */

void pass1_pre_carries (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;

/* In zero-padded FFTs, we must subtract ZPAD0-6 from the first seven result words. */
/* In the radix-4 delay with partial normalization also apply the partial */
/* normalization multipliers to ZPAD0-6. */

	if (gwdata->ZERO_PADDED_FFT && asm_data->this_block <= 6) {
		struct gwasm_data *main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
		unsigned int i;
		double	*fft_data;

		if (gwdata->SCRATCH_SIZE) fft_data = (double *) asm_data->scratch_area;
		else fft_data = (double *) asm_data->DESTARG;
		for (i = asm_data->this_block; i <= 6 && i < asm_data->this_block + asm_data->cache_line_multiplier; i++) {
			if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
				main_thread_asm_data->ZPAD0_6[i] *= gwdata->ZPAD_0_6_ADJUST[i];
				*fft_data -= main_thread_asm_data->ZPAD0_6[i];
			} else {
				*fft_data -= main_thread_asm_data->ZPAD0_6[i] * asm_data->u.xmm.XMM_NORM012_FF[0];
				asm_data->u.xmm.XMM_SUMOUT[0] += main_thread_asm_data->ZPAD0_6[i] * asm_data->u.xmm.XMM_NORM012_FF[0];
			}
			if (gwdata->cpu_flags & CPU_AVX512F)
				fft_data += 16;		// Next low-word/high-word pair (two cache lines)
			else
				fft_data += 8;		// Next cache line
		}
	}
}

/* This routine is called by the pass 1 threads when they are done */
/* with the normalization code. */

#define PASS1_CARRIES_NO_FORWARD_FFT	0
#define PASS1_CARRIES_FORWARD_FFT	1

int pass1_post_carries (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;
	int	i;

/* If postfft is not set, then do not perform forward FFT on the result */

	if (!gwdata->POSTFFT) return (PASS1_CARRIES_NO_FORWARD_FFT);

/* We must delay the forward FFT for the first few data blocks in each */
/* section (until the carries are added back in). */

	i = asm_data->thread_num;
	if (asm_data->this_block < gwdata->pass1_carry_sections[i].start_block + gwdata->num_postfft_blocks)
		return (PASS1_CARRIES_NO_FORWARD_FFT);

/* If this is a zero padded FFT and we are doing the last 3 blocks, then */
/* we need to copy a few of the data values before they are FFTed. */

	if (gwdata->ZERO_PADDED_FFT && asm_data->this_block + asm_data->cache_line_multiplier >= gwdata->num_pass1_blocks - 3)
		xcopy_3_words (asm_data);

/* Perform the forward FFT on these data blocks while they are still in the cache */

	return (PASS1_CARRIES_FORWARD_FFT);
}

/* This routine is called by assembly code threads to get the next */
/* pass 1 block to process.  It returns a code indicating what part */
/* of the pass 1 to do next. */

/* Return codes: */

#define PASS1_DO_MORE_INVERSE_FFT	0
#define PASS1_DO_MORE_FORWARD_FFT	1
#define PASS1_COMPLETE			2
#define PASS1_EXIT_THREAD		3
#define PASS1_START_PASS2		4
#define PASS1_DO_GWCARRIES		5

/* Pass 1 states: */
/* 0 = forward fft */
/* 1 = inverse fft */
/* 999 = not in pass 1, we're doing pass 2 */

int pass1_get_next_block (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;

/* If state is zero, then we are doing the forward FFT.  In this case we */
/* can process the blocks in any order (much like we do in pass 2). */

	if (gwdata->pass1_state == 0) {

/* There are no more blocks to process when the next block is the same */
/* as the block just processed. */

		if (asm_data->this_block == asm_data->next_block) {
			asm_data->DIST_TO_FFTSRCARG = 0;
			return (PASS1_START_PASS2);
		}

/* There is more pass 1 work to do.  Get next available block (if any) for prefetching. */

		asm_data->this_block = asm_data->next_block;
		asm_data->data_addr = asm_data->data_prefetch;
		asm_data->premult_addr = asm_data->premult_prefetch;
		pass1_state0_assign_next_block (gwdata, asm_data);
		return (PASS1_DO_MORE_FORWARD_FFT);
	}

/* Otherwise, pass1_state is one and we are doing the inverse FFT (and */
/* if POSTFFT is set pass 1 of the forward FFT on the result). */

/* Pass 1 state 1 section states: */
/* 0 = section not yet started */
/* 1 = section being processed */
/* 2 = section is adding carries into the next section */
/* 3 = section is doing postfft processing of next section */
/* 4 = section is complete */

/* Handle the common case, we're in the middle of processing a section */

	ASSERTG (gwdata->pass1_carry_sections[0].section_state >= 1 && gwdata->pass1_carry_sections[0].section_state <= 3);
	if (gwdata->pass1_carry_sections[0].section_state == 1) {

/* If there is another block to process in the section, let's process it */

		if (gwdata->pass1_carry_sections[0].next_block != gwdata->pass1_carry_sections[0].last_block) {
			asm_data->this_block = gwdata->pass1_carry_sections[0].next_block;
			gwdata->pass1_carry_sections[0].next_block += asm_data->cache_line_multiplier;
			asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
			asm_data->premult_addr = pass1_premult_addr (gwdata, asm_data->this_block);
			pass1_state1_norm_addrs (gwdata, asm_data);
			pass1_state1_assign_next_block (gwdata, asm_data);
			return (PASS1_DO_MORE_INVERSE_FFT);
		}

/* Handle the case where we just did the last block.  We now need to apply the */
/* carries back to the first blocks. */

		gwdata->pass1_carry_sections[0].section_state = 2;
		asm_data->this_block = 0;
		asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
		pass1_state1_carry_addrs (gwdata, asm_data);
		return (PASS1_DO_GWCARRIES);
	}

/* Handle case where we just propagated carries into the next block. */

	if (gwdata->pass1_carry_sections[0].section_state == 2 && gwdata->POSTFFT) {
		gwdata->pass1_carry_sections[0].section_state = 3;

/* After the carries are added back in to the first block, save the lowest four words of a */
/* zero padded FFTs (before POSTFFT processing destroys the data). */

		if (gwdata->ZERO_PADDED_FFT) xcopy_4_words (asm_data);

/* Set blocks needing postfft processing */

		gwdata->pass1_carry_sections[0].start_block = 0;
		gwdata->pass1_carry_sections[0].last_block = gwdata->num_postfft_blocks;
		gwdata->pass1_carry_sections[0].next_block = 0;
	}

/* Do the forward FFT on the result blocks affected by our carries */

	if (gwdata->pass1_carry_sections[0].section_state == 3 &&
	    gwdata->pass1_carry_sections[0].next_block != gwdata->pass1_carry_sections[0].last_block) {
		asm_data->this_block = gwdata->pass1_carry_sections[0].next_block;
		gwdata->pass1_carry_sections[0].next_block += asm_data->cache_line_multiplier;
		asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
		asm_data->premult_addr = pass1_premult_addr (gwdata, asm_data->this_block);
		pass1_state1_assign_next_block (gwdata, asm_data);
		return (PASS1_DO_MORE_FORWARD_FFT);
	}

/* We're finished, do final cleanup */

	if (gwdata->cpu_flags & CPU_AVX512F) {
	} else if (gwdata->cpu_flags & CPU_AVX) {
		if (asm_data->u.ymm.YMM_MAXERR[0] > asm_data->MAXERR) asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[0];
		if (asm_data->u.ymm.YMM_MAXERR[1] > asm_data->MAXERR) asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[1];
		if (asm_data->u.ymm.YMM_MAXERR[2] > asm_data->MAXERR) asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[2];
		if (asm_data->u.ymm.YMM_MAXERR[3] > asm_data->MAXERR) asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[3];
	} else if (gwdata->cpu_flags & CPU_SSE2) {
		if (asm_data->u.xmm.XMM_MAXERR[0] > asm_data->MAXERR) asm_data->MAXERR = asm_data->u.xmm.XMM_MAXERR[0];
		if (asm_data->u.xmm.XMM_MAXERR[1] > asm_data->MAXERR) asm_data->MAXERR = asm_data->u.xmm.XMM_MAXERR[1];
		/* Accumulate SUMOUT */
		((double *) asm_data->DESTARG)[-3] = asm_data->u.xmm.XMM_SUMOUT[0] + asm_data->u.xmm.XMM_SUMOUT[1];
		/* Normalize SUMOUT value by multiplying by 1 / (fftlen/2). */
		((double *) asm_data->DESTARG)[-3] *= asm_data->ttmp_ff_inv;
	}

	return (PASS1_COMPLETE);
}

int pass1_get_next_block_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;
	int	i;

/* Grab lock before reading or writing any gwdata values. */

	gwmutex_lock (&gwdata->thread_lock);

/* If state is zero, then we are doing the forward FFT.  In this case we */
/* can process the blocks in any order (much like we do in pass 2). */

	if (gwdata->pass1_state == 0) {

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  In that case, if this is the main thread */
/* wait for auxiliary threads to complete.  If this is an auxiliary thread */
/* then return code telling the assembly code to exit. */

		if (asm_data->this_block == asm_data->next_block) {
			if (asm_data->thread_num) {
				gwmutex_unlock (&gwdata->thread_lock);
				return (PASS1_EXIT_THREAD);
			}
			gwmutex_unlock (&gwdata->thread_lock);
			gwevent_wait (&gwdata->all_threads_done, 0);
			asm_data->DIST_TO_FFTSRCARG = 0;
			return (PASS1_START_PASS2);
		}

/* There is more pass 1 work to do.  Get next available block (if any) for prefetching. */

		asm_data->this_block = asm_data->next_block;
		asm_data->data_addr = asm_data->data_prefetch;
		asm_data->premult_addr = asm_data->premult_prefetch;
		pass1_state0_assign_next_block (gwdata, asm_data);
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_DO_MORE_FORWARD_FFT);
	}

/* Otherwise, pass1_state is one and we are doing the inverse FFT (and */
/* if POSTFFT is set pass 1 of the forward FFT on the result). */

/* Pass 1 state 1 section states: */
/* 0 = section not yet started */
/* 1 = section being processed */
/* 2 = section is adding carries into the next section */
/* 3 = section is doing postfft processing of next section */
/* 4 = section complete */

/* Handle the common case, we're in the middle of processing a section */

	i = asm_data->thread_num;
	ASSERTG (gwdata->pass1_carry_sections[i].section_state >= 1 && gwdata->pass1_carry_sections[i].section_state <= 3);
	if (gwdata->pass1_carry_sections[i].section_state == 1) {

/* If we just finished off the first blocks in the section, then set dependent section's OK to carry into flag */
/* Signal an event in case the dependent section is waiting on us to process our first blocks. */

		if (asm_data->this_block + asm_data->cache_line_multiplier ==
				gwdata->pass1_carry_sections[i].start_block + gwdata->num_postfft_blocks) {
			int	dep = gwdata->pass1_carry_sections[i].dependent_section;
			ASSERTG (gwdata->pass1_carry_sections[dep].last_block == gwdata->pass1_carry_sections[i].start_block ||
				 (gwdata->pass1_carry_sections[dep].last_block == gwdata->num_pass1_blocks &&
				  gwdata->pass1_carry_sections[i].start_block == 0));
			gwdata->pass1_carry_sections[dep].can_carry_into_next = TRUE;
			gwdata->pass1_carry_sections[i].dependent_section = -1;
			gwevent_signal (&gwdata->can_carry_into);
		}

/* If there is another block to process in the section, let's process it */

		if (gwdata->pass1_carry_sections[i].next_block != gwdata->pass1_carry_sections[i].last_block) {
			asm_data->this_block = gwdata->pass1_carry_sections[i].next_block;
			gwdata->pass1_carry_sections[i].next_block += asm_data->cache_line_multiplier;
			asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
			asm_data->premult_addr = pass1_premult_addr (gwdata, asm_data->this_block);
			pass1_state1_norm_addrs (gwdata, asm_data);
			pass1_state1_assign_next_block (gwdata, asm_data);
			gwmutex_unlock (&gwdata->thread_lock);
			return (PASS1_DO_MORE_INVERSE_FFT);
		}

/* Handle the case where we just did the last block in the section.  We now need to apply this */
/* section's carries to the next section.  Unfortunately, the next section may not have processed */
/* it's first blocks yet forcing us to wait.  Technically, each section should have its own event */
/* as there is no guarantee that the signal will wake up the thread waiting here.  However, */
/* waiting here should be rare as it only happens when one thread completes an entire section */
/* before the dependent thread has processed its first blocks. */

		while (! gwdata->pass1_carry_sections[i].can_carry_into_next) {
			gwevent_reset (&gwdata->can_carry_into);
			gwmutex_unlock (&gwdata->thread_lock);
			gwevent_wait (&gwdata->can_carry_into, 0);
			gwmutex_lock (&gwdata->thread_lock);
		}

/* Now apply this section's carries to the next section. */

		gwdata->pass1_carry_sections[i].section_state = 2;
		asm_data->this_block = gwdata->pass1_carry_sections[i].last_block;
		if (asm_data->this_block == gwdata->num_pass1_blocks) asm_data->this_block = 0;
		/* Copy the 7 ZPAD values that were computed in the main thread */
		if (gwdata->ZERO_PADDED_FFT && asm_data->this_block == 0) {
			struct gwasm_data *main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
			memcpy (asm_data->ZPAD0_6, main_thread_asm_data->ZPAD0_6, 7 * sizeof (double));
		}
		asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
		pass1_state1_carry_addrs (gwdata, asm_data);
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_DO_GWCARRIES);
	}

/* Handle case where we just propagated carries into the next block and we need to start the */
/* forward FFT on the affected blocks. */

	if (gwdata->pass1_carry_sections[i].section_state == 2 && gwdata->POSTFFT) {
		gwdata->pass1_carry_sections[i].section_state = 3;

/* After the carries are added back in to the first blocks, save the lowest four words of a */
/* zero padded FFTs (before POSTFFT processing destroys the data). */

		if (gwdata->ZERO_PADDED_FFT) {
			if (asm_data->this_block == 0) xcopy_4_words (asm_data);

/* After the carries are added back in to the last blocks, save the three words below the halfway */
/* point in a zero padded FFTs (before POSTFFT processing destroys the data).  This can only happen */
/* due to section splitting or a huge number of threads. */

			if (asm_data->this_block + gwdata->num_postfft_blocks >= gwdata->num_pass1_blocks - 3)
				xcopy_3_words_after_gwcarries (asm_data);
		}

/* Figure out which blocks need postfft processing */

		gwdata->pass1_carry_sections[i].start_block = gwdata->pass1_carry_sections[i].next_block;
		if (gwdata->pass1_carry_sections[i].start_block == gwdata->num_pass1_blocks) gwdata->pass1_carry_sections[i].start_block = 0;
		gwdata->pass1_carry_sections[i].last_block = gwdata->pass1_carry_sections[i].start_block + gwdata->num_postfft_blocks;
		gwdata->pass1_carry_sections[i].next_block = gwdata->pass1_carry_sections[i].start_block;
	}

/* Do the forward FFT on the postfft blocks */

	if (gwdata->pass1_carry_sections[i].section_state == 3 &&
	    gwdata->pass1_carry_sections[i].next_block != gwdata->pass1_carry_sections[i].last_block) {
		asm_data->this_block = gwdata->pass1_carry_sections[i].next_block;
		gwdata->pass1_carry_sections[i].next_block += asm_data->cache_line_multiplier;
		asm_data->data_addr = pass1_data_addr (gwdata, asm_data, asm_data->this_block);
		asm_data->premult_addr = pass1_premult_addr (gwdata, asm_data->this_block);
		pass1_state1_assign_next_block (gwdata, asm_data);
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_DO_MORE_FORWARD_FFT);
	}

/* We've finished this section */

	gwdata->pass1_carry_sections[i].section_state = 4;

/* See if there is some other section we can help finish off which will set the section_state back to 1 */

	if (pass1_state1_assign_first_block (gwdata, asm_data)) {
		pass1_state1_assign_next_block (gwdata, asm_data);
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_DO_MORE_INVERSE_FFT);
	}
   
/* We've finished this thread, merge this thread's maxerr and sumout with the main thread */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		struct gwasm_data *main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
		if (asm_data->MAXERR > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->MAXERR;
	} else if (gwdata->cpu_flags & CPU_AVX) {
		struct gwasm_data *main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
		if (asm_data->u.ymm.YMM_MAXERR[0] > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[0];
		if (asm_data->u.ymm.YMM_MAXERR[1] > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[1];
		if (asm_data->u.ymm.YMM_MAXERR[2] > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[2];
		if (asm_data->u.ymm.YMM_MAXERR[3] > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->u.ymm.YMM_MAXERR[3];
	} else if (gwdata->cpu_flags & CPU_SSE2) {
		struct gwasm_data *main_thread_asm_data = (struct gwasm_data *) gwdata->asm_data;
		if (asm_data->u.xmm.XMM_MAXERR[0] > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->u.xmm.XMM_MAXERR[0];
		if (asm_data->u.xmm.XMM_MAXERR[1] > main_thread_asm_data->MAXERR) main_thread_asm_data->MAXERR = asm_data->u.xmm.XMM_MAXERR[1];
		/* Accumulate SUMOUT */
		((double *) asm_data->DESTARG)[-3] += asm_data->u.xmm.XMM_SUMOUT[0] + asm_data->u.xmm.XMM_SUMOUT[1];
	}

/* There are no more blocks to process.  If this is an auxiliary thread */
/* then return code telling the assembly code to exit. */

	if (asm_data->thread_num) {
		gwmutex_unlock (&gwdata->thread_lock);
		return (PASS1_EXIT_THREAD);
	}

/* There are no more blocks to process.  This is the main thread, */
/* wait for auxiliary threads to complete. */

	gwmutex_unlock (&gwdata->thread_lock);
	gwevent_wait (&gwdata->all_threads_done, 0);

/* We're finished, do final cleanup */

	/* Normalize SUMOUT value by multiplying by 1 / (fftlen/2). */
	((double *) asm_data->DESTARG)[-3] *= asm_data->ttmp_ff_inv;

	return (PASS1_COMPLETE);
}

/* This callback routine is called by the main thread assembly code to */
/* fire up all the auxiliary worker threads in pass 2. */

void pass2_wake_up_threads (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;

/* Obtain the lock if multithreading. */
/* Why?  The other threads should be off right now, but we've seen in the debugger that it */
/* is possible to get here with the lock held by an auxiliary thread that is finishing up. */

	if (gwdata->num_threads > 1 && gwdata->PASS1_SIZE) {
		gwmutex_lock (&gwdata->thread_lock);
	}

/* Call assign_block twice to set up the main thread's this_block and next_block values. */

	gwdata->pass1_state = PASS1_STATE_PASS2;
	gwdata->next_block = 0;
	pass2_assign_first_block (gwdata, asm_data);
	pass2_assign_next_block (gwdata, asm_data);

/* Signal the auxiliary threads to resume working */

	if (gwdata->num_threads > 1 && gwdata->PASS1_SIZE) {
		gwdata->catch_straggler_threads = FALSE;
		gwevent_reset (&gwdata->all_threads_done);
		gwevent_signal (&gwdata->thread_work_to_do);
		gwmutex_unlock (&gwdata->thread_lock);
	}
}

/* This routine is called by assembly code threads to get the next */
/* pass 2 block to process. */

int pass2_get_next_block (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata;

	if (asm_data->this_block == asm_data->next_block) {
		asm_data->DIST_TO_FFTSRCARG = 0;
		return (TRUE);
	}

	gwdata = asm_data->gwdata;

	asm_data->this_block = asm_data->next_block;
	asm_data->data_addr = asm_data->data_prefetch;
	asm_data->premult_addr = asm_data->premult_prefetch;
	pass2_assign_next_block (gwdata, asm_data);

	return (FALSE);
}

/* This is the multi-thread version of the routine above */

int pass2_get_next_block_mt (
	struct gwasm_data *asm_data)
{
	gwhandle *gwdata = asm_data->gwdata;

/* There are no more blocks to process when the next block is the same */
/* as the block just processed.  In that case, if this is the main thread */
/* wait for auxiliary threads to complete.  If this is an auxiliary thread */
/* then return code telling the assembly code to exit. */

	if (asm_data->this_block == asm_data->next_block) {
		if (asm_data->thread_num) return (TRUE);
		gwevent_wait (&gwdata->all_threads_done, 0);
		asm_data->DIST_TO_FFTSRCARG = 0;
		return (TRUE);
	}

/* Grab lock before reading or writing any gwdata values. */

	gwmutex_lock (&gwdata->thread_lock);

/* Copy prefetched block and addresses to this block.  Get next available */
/* block to prefetch. */

	asm_data->this_block = asm_data->next_block;
	asm_data->data_addr = asm_data->data_prefetch;
	asm_data->premult_addr = asm_data->premult_prefetch;
	pass2_assign_next_block (gwdata, asm_data);
	gwmutex_unlock (&gwdata->thread_lock);

/* Return code indicating more work to do */

	return (FALSE);
}


/* Perform initializations required for multi-threaded operation */

int multithread_init (
	gwhandle *gwdata)
{
	struct gwasm_data *asm_data;
	unsigned long i;

/* Only two pass AVX-512/AVX/SSE2 FFTs support multi-threaded execution */

	if (gwdata->PASS2_SIZE == 0 || !(gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2))) {
		gwdata->num_threads = 1;
		return (0);
	}

/* Get pointer to assembly structure */

	asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Save gwdata pointer in asm_data so that C callback routines can */
/* access gwdata.  Set flag indicating this is the main thread. */

	asm_data->gwdata = gwdata;
	asm_data->thread_num = 0;

/* Init other variables */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		if (gwdata->PASS1_SIZE == 0) {
			gwdata->num_pass1_blocks = 1;
			gwdata->num_pass2_blocks = 1;
			asm_data->last_pass1_block = 0;
		} else {
			gwdata->num_pass1_blocks = gwdata->PASS2_SIZE;
			gwdata->num_pass2_blocks = gwdata->PASS1_SIZE >> 1;
			asm_data->last_pass1_block = gwdata->num_pass1_blocks - asm_data->cache_line_multiplier;
		}
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		gwdata->num_pass1_blocks = gwdata->PASS2_SIZE;
		gwdata->num_pass2_blocks = gwdata->PASS1_SIZE >> 1;
		asm_data->last_pass1_block = gwdata->num_pass1_blocks - asm_data->cache_line_multiplier;
	} else {
		gwdata->num_pass1_blocks = gwdata->PASS2_SIZE >> 1;
		gwdata->num_pass2_blocks = gwdata->PASS1_SIZE >> 2;
		asm_data->last_pass1_block = gwdata->num_pass1_blocks - asm_data->cache_line_multiplier;
	}

/* Place a limit on the number of threads */

	if (gwdata->num_threads > gwdata->num_pass1_blocks / asm_data->cache_line_multiplier)
		gwdata->num_threads = gwdata->num_pass1_blocks / asm_data->cache_line_multiplier;
	if (gwdata->num_threads > gwdata->num_pass1_blocks / 8)
		gwdata->num_threads = gwdata->num_pass1_blocks / 8;

/* Determine how many data blocks are affected by carries out of pass 1 section.  Zero-padded FFTs require 8 words */
/* to propagate carries into.  For AVX-512 FFTs, znorm012_wpn can spread the carry over a maximum of 8 words. */
/* For AVX FFTs, ynorm012_wpn can spread the carry over a maximum of either 4 or 8 words. */
/* For SSE2 FFTs, xnorm012_2d and xnorm012_2d_wpn spreads carries over either 2 or 6 words. */

	if (gwdata->ZERO_PADDED_FFT)
		gwdata->num_postfft_blocks = 8;
	else if (gwdata->cpu_flags & CPU_AVX512F) {
		gwdata->num_postfft_blocks = (int) floor (53.0 / (gwdata->avg_num_b_per_word * log2 (gwdata->b)));
		ASSERTG (gwdata->num_postfft_blocks <= 8);
	} else if (gwdata->cpu_flags & CPU_AVX) {
		gwdata->num_postfft_blocks = (int) floor (53.0 / (gwdata->avg_num_b_per_word * log2 (gwdata->b)));
		ASSERTG (gwdata->num_postfft_blocks <= 8);
		asm_data->SPREAD_CARRY_OVER_EXTRA_WORDS = (gwdata->num_postfft_blocks > 4);
	} else {
		if (asm_data->SPREAD_CARRY_OVER_EXTRA_WORDS) gwdata->num_postfft_blocks = 6;
		else gwdata->num_postfft_blocks = 2;
	}
	gwdata->num_postfft_blocks = round_up_to_multiple_of (gwdata->num_postfft_blocks, asm_data->cache_line_multiplier);

/* Calculate the values used to compute pass 2 premultier pointers. */
/* This only happens for our home-grown FFTs.  Non-power-of-2 pass 2 */
/* sizes are not supported. */
/* We calculate an adjusted starting address of the premultiplier data */
/* so that both real and all-complex FFTs can use the same formula to */
/* calculate the proper address given a block number. */

	if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
		gwdata->pass2_premult_block_size =
			(gwdata->PASS2_SIZE == 256) ? 32 * 128 :
			(gwdata->PASS2_SIZE == 1024) ? 64 * 128 :
			(gwdata->PASS2_SIZE == 2048) ? 96 * 128 :
			(gwdata->PASS2_SIZE == 4096) ? 128 * 128 :
			(gwdata->PASS2_SIZE == 8192) ? 192 * 128 : 0;
		if (gwdata->ALL_COMPLEX_FFT)
			gwdata->adjusted_pass2_premults = asm_data->u.xmm.pass2_premults;
		else
			gwdata->adjusted_pass2_premults =
				(char *) asm_data->u.xmm.pass2_premults - gwdata->pass2_premult_block_size;
	}

/* Allocate and init carry section array */

	gwdata->pass1_carry_sections = (struct pass1_carry_sections *) malloc (gwdata->num_threads * sizeof (struct pass1_carry_sections));
	if (gwdata->pass1_carry_sections == NULL) return (GWERROR_MALLOC);

/* Create the prefetching hyperthread for the main compute thread */

	create_auxiliary_hyperthread (asm_data);

/* If we aren't multithreading, use the simpler version of routines */

	if (gwdata->num_threads <= 1) {
		asm_data->pass1_wake_up_threads = pass1_wake_up_threads;
		asm_data->pass1_pre_carries = pass1_pre_carries;
		asm_data->pass1_post_carries = pass1_post_carries;
		asm_data->pass1_get_next_block = pass1_get_next_block;
		asm_data->pass2_wake_up_threads = pass2_wake_up_threads;
		asm_data->pass2_get_next_block = pass2_get_next_block;
		return (0);
	}

/* Init thread arrays */

	gwdata->thread_ids = (gwthread *) malloc (gwdata->num_threads * sizeof (gwthread));
	if (gwdata->thread_ids == NULL) return (GWERROR_MALLOC);
	memset (gwdata->thread_ids, 0, gwdata->num_threads * sizeof (gwthread));

/* Init mutexes and events used to control auxiliary threads */

	gwmutex_init (&gwdata->thread_lock);
	gwevent_init (&gwdata->thread_work_to_do);
	gwevent_init (&gwdata->all_threads_done);
	gwevent_init (&gwdata->can_carry_into);
	gwdata->num_active_threads = 0;
	gwevent_signal (&gwdata->all_threads_done);

/* Set ptrs to call back routines in structure used by assembly code */

	asm_data->pass1_wake_up_threads = pass1_wake_up_threads;
	asm_data->pass1_pre_carries = pass1_pre_carries;
	asm_data->pass1_post_carries = pass1_post_carries;
	asm_data->pass1_get_next_block = pass1_get_next_block_mt;
	asm_data->pass2_wake_up_threads = pass2_wake_up_threads;
	asm_data->pass2_get_next_block = pass2_get_next_block_mt;

/* Pre-create each auxiliary thread used in multiplication code. */
/* We allocate the memory here so that error recovery is easier. */

	gwdata->threads_must_exit = FALSE;
	for (i = 0; i < gwdata->num_threads - 1; i++) {
		struct thread_data *info;

		info = (struct thread_data *) malloc (sizeof (struct thread_data));
		if (info == NULL) return (GWERROR_MALLOC);
		info->gwdata = gwdata;
		info->thread_num = i+1;
		/* Allocate the asm_data area and thread stack */
		info->asm_data_alloc = aligned_malloc (sizeof (struct gwasm_data) + NEW_STACK_SIZE, 4096);
		if (info->asm_data_alloc == NULL) {
			free (info);
			return (GWERROR_MALLOC);
		}
		/* Allocate the scratch area */
		info->scratch_area = aligned_malloc (gwdata->SCRATCH_SIZE, 128);
		if (info->scratch_area == NULL) {
			aligned_free (info->asm_data_alloc);
			free (info);
			return (GWERROR_MALLOC);
		}
		/* Allocate the carries area */
		if (gwdata->cpu_flags & CPU_AVX512F) {
			int	carry_table_size = gwdata->PASS1_SIZE;
			info->carries = aligned_malloc (carry_table_size * sizeof (double), 128);
		} else if (gwdata->cpu_flags & CPU_AVX) {
			int	carry_table_size = gwdata->PASS1_SIZE;
			info->carries = aligned_malloc (carry_table_size * sizeof (double), 128);
		} else {
			int	carry_table_size = gwdata->PASS1_SIZE * 2;
			info->carries = aligned_malloc (carry_table_size * sizeof (double), 128);
		}
		if (info->carries == NULL) {
			aligned_free (info->asm_data_alloc);
			aligned_free (info->scratch_area);
			free (info);
			return (GWERROR_MALLOC);
		}
		/* Launch the auxiliary thread */
		gwthread_create_waitable (&gwdata->thread_ids[i], &auxiliary_thread, info);
	}

/* Return success */

	return (0);
}

/* Routines for auxiliary hyperthreads that do prefetching */

void create_auxiliary_hyperthread (			/* Create the prefetching hyperthread for one of the compute threads */
	struct gwasm_data *asm_data)
{
#ifdef AUX_HYPER
	gwhandle *gwdata;

	gwdata = asm_data->gwdata;
	if (gwdata->hyperthread_prefetching) {
		gwevent_init (&asm_data->hyperthread_work_to_do);
		gwthread_create_waitable (&asm_data->hyperthread_id, &auxiliary_hyperthread, asm_data);
	}
}

void auxiliary_hyperthread (void *arg)
{
	struct gwasm_data *asm_data;
	gwhandle *gwdata;
	char	*prefetch_ptr;
static int wakeups = 0;
static int prefetches1 = 0;
static int prefetches2 = 0;

/* Get pointers to various structures.  The compute thread and corresponding prefetching hyperthread share the same asm_data. */

	asm_data = (struct gwasm_data *) arg;
	gwdata = asm_data->gwdata;

/* Call optional user provided callback routine so that the caller can set */
/* the thread's priority and affinity */

//BUG - if thread_callback is null, disable hyperthread_prefetching?? or maybe just prefetch to the L3 cache -- might be useful on a Mac?
	if (gwdata->thread_callback != NULL)
		(*gwdata->thread_callback) (asm_data->thread_num, 10, gwdata->thread_callback_data);

/* Loop waiting for work to do.  Each compute thread will signal there is prefething work to do. */

	for ( ; ; ) {
		gwevent_wait (&asm_data->hyperthread_work_to_do, 0);

/* If threads are to exit, break out of this work loop */

		if (gwdata->threads_must_exit) break;

/* Reset hyperthread_work_to_do event so we can wait for more work when we loop */

		gwevent_reset (&asm_data->hyperthread_work_to_do);
wakeups++;

//BUG -  do I need locks in case where compute thread completes before hyperthread even starts up?
// do we need to detect early completion of the compute thread??  Or is the worst case scenario the
// hyperthread prefetches the wrong data?


/* If this is the last block in pass 1 or pass 2, then there is no need to do any prefetching. */
/* Actually, this should not happen.  We don't signal in this case. */
//BUG - in theory, the last pass 1 block could prefetch the first pass 2 block (and vice versa).  This would require
// knowing each thread_number's first block...

		if (asm_data->this_block == asm_data->next_block) continue;

/* Determine which cache lines need to be prefetched.  Obviously, this depends on whether we are in pass 1 or 2. */

//BUG - test different prefetching options:  1) to L2, L3, L2 with write intent,  2) prefetch current block
// data (plus variable & fixed s/c data and group multipliers) before next block data for optimal LRU settings
// (assumes prefetch of data already in cache is quick), 3) use NOPs or delays between prefetch instructions

		if (gwdata->pass1_state < PASS1_STATE_PASS2) {			// In pass 1
			int	pass1_size, i;
			int	current_block;

			// Remember the current block locally so we can detect the compute hyperthread outpacing this prefetch hyperthread
			current_block = asm_data->this_block;

			// Prefetch the FFT data
			prefetch_ptr = pass1_data_addr (gwdata, asm_data, asm_data->next_block);
			pass1_size = gwdata->PASS1_SIZE;
			for (i = 0; i < pass1_size / 2; i++) {
				// Check for change in compute thread's current block
				if (gwdata->pass1_state >= PASS1_STATE_PASS2 || current_block != asm_data->this_block) break;
				// Prefetch some cache lines
				prefetchL2 (prefetch_ptr, asm_data->cache_line_multiplier / 4);
prefetches1 += asm_data->cache_line_multiplier / 4;
				prefetch_ptr += asm_data->pass2blkdst;
			}

			// Check for change in compute thread's current block
			if (gwdata->pass1_state >= PASS1_STATE_PASS2 || current_block != asm_data->this_block) continue;

			// Prefetch the sin/cos data
			prefetch_ptr = pass1_premult_addr (gwdata, asm_data->next_block);
// BUG - break this into smaller chunks?
			prefetchL2 (prefetch_ptr, gwdata->pass1_var_data_size / 64);
prefetches1 += gwdata->pass1_var_data_size / 64;
		}

		// Pass 2 prefetching
		else {
			int	current_block;

			// Remember the current block locally so we can detect the compute hyperthread outpacing this prefetch hyperthread
			current_block = asm_data->this_block;

			// Check for change in compute thread's current block
			if (gwdata->pass1_state < PASS1_STATE_PASS2 || current_block != asm_data->this_block) continue;

			// Prefetch the FFT data
// BUG - break this into smaller chunks?
			prefetch_ptr = pass2_data_addr (gwdata, asm_data, asm_data->next_block);
			prefetchL2 (prefetch_ptr, (int) (asm_data->pass2blkdst / 64));
			//BUG - skip over the FOURKBGAPS
prefetches2 += (int) asm_data->pass2blkdst / 64;
		}
	}

/* Call optional user provided callback routine so that the caller can do any necessary cleanup. */

	if (gwdata->thread_callback != NULL)
		(*gwdata->thread_callback) (asm_data->thread_num, 11, gwdata->thread_callback_data);
#endif
}


/* Perform cleanup required by multi-threaded operation */

void multithread_term (
	gwhandle *gwdata)
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	int	i;

/* If we are multithreading AND multithreading was initialized properly, then wait for compute threads to exit */
/* before we free memory. */

	if (gwdata->num_threads > 1 && gwdata->thread_lock != NULL) {

/* Set variable that tells all auxiliary threads and hyperthreads to terminate. */
/* Fire up auxiliary compute threads so they can exit. */

		gwmutex_lock (&gwdata->thread_lock);
		gwdata->threads_must_exit = TRUE;
		gwevent_signal (&gwdata->thread_work_to_do);
		gwmutex_unlock (&gwdata->thread_lock);

/* Wait for all compute threads to exit.  We must do this so */
/* that this thread can safely delete the gwdata structure */

		for (i = 0; i < (int) gwdata->num_threads - 1; i++)
			if (gwdata->thread_ids[i]) gwthread_wait_for_exit (&gwdata->thread_ids[i]);

/* Free up multithreading memory */

		free (gwdata->thread_ids);
		gwdata->thread_ids = NULL;

/* Free up the multithread resources */

		gwmutex_destroy (&gwdata->thread_lock);
		gwevent_destroy (&gwdata->thread_work_to_do);
		gwevent_destroy (&gwdata->all_threads_done);
		gwevent_destroy (&gwdata->can_carry_into);
		gwdata->thread_lock = NULL;
	}

/* Wait for the optional prefetching hyperthread to finish */

	if (asm_data != NULL && asm_data->hyperthread_id) {
		gwevent_signal (&asm_data->hyperthread_work_to_do);
		gwthread_wait_for_exit (&asm_data->hyperthread_id);
		gwevent_destroy (&asm_data->hyperthread_work_to_do);
	}

/* Free up memory allocated by multithread_init */

	free (gwdata->pass1_carry_sections);
	gwdata->pass1_carry_sections = NULL;
}

/* Cleanup any memory allocated for multi-precision math */

void gwdone (
	gwhandle *gwdata)	/* Handle returned by gwsetup */
{
	unsigned int i;

	multithread_term (gwdata);

	term_ghandle (&gwdata->gdata);
	if (gwdata->asm_data != NULL) {
		struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
		unshare_sincos_data (asm_data->sincos1);			// SSE2
		unshare_sincos_data (asm_data->sincos2);			// SSE2 & AVX & AVX512
		unshare_sincos_data (asm_data->xsincos_complex);		// SSE2 & AVX & AVX512
		unshare_sincos_data (asm_data->sincos3);			// SSE2 & AVX & AVX512
		aligned_free ((char *) gwdata->asm_data - NEW_STACK_SIZE);
		gwdata->asm_data = NULL;
	}
	free (gwdata->dd_data);
	gwdata->dd_data = NULL;
	free (gwdata->gwnum_free);
	gwdata->gwnum_free = NULL;
	if (gwdata->gwnum_alloc != NULL) {
		for (i = 0; i < gwdata->gwnum_alloc_count; i++) {
			char	*p;
			int32_t	freeable;
			p = (char *) gwdata->gwnum_alloc[i];
			freeable = * (int32_t *) (p - 32) & ~GWFREED_TEMPORARILY;
			if (freeable) aligned_free ((char *) p - GW_HEADER_SIZE);
		}
		free (gwdata->gwnum_alloc);
		gwdata->gwnum_alloc = NULL;
	}
	free (gwdata->GW_MODULUS);
	gwdata->GW_MODULUS = NULL;
	if (gwdata->large_pages_ptr != NULL) {
		large_pages_free (gwdata->large_pages_ptr);
		gwdata->large_pages_ptr = NULL;
	} else {
		aligned_free (gwdata->gwnum_memory);
		gwdata->gwnum_memory = NULL;
	}
}

/* Routine to allocate aligned memory for our big numbers */
/* Memory is allocated on 128-byte boundaries, with an additional */
/* 32 bytes prior to the data for storing useful stuff */

gwnum gwalloc (
	gwhandle *gwdata)
{
	unsigned long size, aligned_size;
	char	*p, *q;
	int32_t	freeable;

/* Return cached gwnum if possible */

	if (gwdata->gwnum_free_count)
		return (gwdata->gwnum_free[--gwdata->gwnum_free_count]);

/* Allocate arrays if necessary */

	if (gwdata->gwnum_alloc == NULL) {
		gwdata->gwnum_free = (gwnum *)
			malloc (gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_free == NULL) return (NULL);
		gwdata->gwnum_alloc = (gwnum *)
			malloc (gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_alloc == NULL) return (NULL);
	} else if (gwdata->gwnum_alloc_count == gwdata->gwnum_alloc_array_size) {
		gwdata->gwnum_alloc_array_size += gwdata->gwnum_alloc_array_size >> 1;
		gwdata->gwnum_free = (gwnum *)
			realloc (gwdata->gwnum_free,
				 gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_free == NULL) return (NULL);
		gwdata->gwnum_alloc = (gwnum *)
			realloc (gwdata->gwnum_alloc,
				 gwdata->gwnum_alloc_array_size * sizeof (gwnum));
		if (gwdata->gwnum_alloc == NULL) return (NULL);
	}

/* Compute the amount of memory to allocate. */
/* Allocate 96 extra bytes for header information and align the data */
/* appropriately.  When allocating memory out of the big buffer for */
/* the torture test, then only allocate data on 128 byte boundaries */
/* to maximize the number of gwnums allocated. */

	size = gwnum_datasize (gwdata);
	aligned_size = (size + GW_HEADER_SIZE + 127) & ~127;
	if (gwdata->large_pages_gwnum != NULL) {
		p = (char *) gwdata->large_pages_gwnum;
		gwdata->large_pages_gwnum = NULL;
		p += -GW_HEADER_SIZE & 127;
		freeable = 0;
	}
	else if (gwdata->GW_BIGBUF_SIZE >= size + aligned_size) {
		p = gwdata->GW_BIGBUF;
		gwdata->GW_BIGBUF += aligned_size;
		gwdata->GW_BIGBUF_SIZE -= aligned_size;
		p += -GW_HEADER_SIZE & 127;
		freeable = 0;
	}
/* FreeBSD supports supports large pages (superpages) automatically.  Allocate the */
/* first gwnum on a 4MB to minimize fragmentation across superpage boundaries.  The */
/* theory is this will maximize the number of superpage promotions. */
#ifdef __FreeBSD__
	else if (size >= 0x400000 && gwdata->gwnum_alloc_count == 0) {
		p = (char *) aligned_offset_malloc (
				size + GW_HEADER_SIZE, 0x400000,
				(GW_HEADER_SIZE - gwdata->GW_ALIGNMENT_MOD) & 0x3FFFFF);
		if (p == NULL) return (NULL);
		freeable = 1;
	}
#endif
	else {
		p = (char *) aligned_offset_malloc (
				size + GW_HEADER_SIZE, gwdata->GW_ALIGNMENT,
				(GW_HEADER_SIZE - gwdata->GW_ALIGNMENT_MOD) &
					(gwdata->GW_ALIGNMENT - 1));
		if (p == NULL) return (NULL);
		freeable = 1;
	}

/* Do a seemingly pointless memset!  This actually is very important. */
/* The memset will walk through the allocated memory sequentially, which */
/* increases the likelihood that contiguous virtual memory will map to */
/* contiguous physical memory.  The FFTs, especially the larger ones, */
/* optimizes L2 cache line collisions on the assumption that the FFT data */
/* is in contiguous physical memory.  Failure to do this results in as */
/* much as a 30% performance hit in an SSE2 2M FFT. */

	q = p + GW_HEADER_SIZE;
	memset (q, 0, size);

/* Initialize the header */

	* (uint32_t *) (q - 8) = size;	/* Size in bytes */
	* (uint32_t *) (q - 4) = 1;	/* Unnormalized adds count */
	* (uint32_t *) (q - 28) = 0;	/* Has-been-pre-ffted flag */
	* (int32_t *) (q - 32) = freeable; /* Mem should be freed flag */
	* (double *) (q - 16) = 0.0;
	* (double *) (q - 24) = 0.0;

/* Save pointer for easier cleanup */

	gwdata->gwnum_alloc[gwdata->gwnum_alloc_count++] = (gwnum) q;

/* Return the gwnum */

	return ((gwnum) q);
}

/* Free one of our special numbers */

void gwfree (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	q)		/* Number to free */
{
	if (gwdata->gwnum_free != NULL && q != NULL)
		gwdata->gwnum_free[gwdata->gwnum_free_count++] = q;
}

/* Typeless gwalloc and gwfree routines for giants code to call */

void *gwgiantalloc (
	void	*gwdata)
{
	return ((void *) gwalloc ((gwhandle *) gwdata));
}
void gwgiantfree (
	void	*gwdata,
	void	*q)
{
	gwfree ((gwhandle *) gwdata, (gwnum) q);
}

/* Specialized routine that allows giants code to deallocate one */
/* cached gwnum to free up space for allocating FFT giant's sincos data. */

void gwgiantdealloc (
	void	*gwdata_arg)
{
	gwhandle *gwdata;
	gwnum	p;
	int32_t	freeable;
	unsigned long i, j;

	gwdata = (gwhandle *) gwdata_arg;
	for (i = 0; i < gwdata->gwnum_free_count; i++) {
		p = gwdata->gwnum_free[i];
		freeable = * (int32_t *) ((char *) p - 32);
		if (freeable & GWFREED_TEMPORARILY) continue;
		if (!freeable) continue;

		for (j = 0; j < gwdata->gwnum_alloc_count; j++) {
			if (gwdata->gwnum_alloc[j] != p) continue;

			aligned_free ((char *) p - GW_HEADER_SIZE);
			gwdata->gwnum_free[i] = gwdata->gwnum_free[--gwdata->gwnum_free_count];
			gwdata->gwnum_alloc[j] = gwdata->gwnum_alloc[--gwdata->gwnum_alloc_count];
			return;
		}
	}
}


/* Specialized routines that let the giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	q)
{
	* (int32_t *) ((char *) q - 32) |= GWFREED_TEMPORARILY;
	gwfree (gwdata, q);
}
void gwrealloc_temporarily (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	q)
{
	unsigned long i, j;

	ASSERTG (* (int32_t *) ((char *) q - 32) & GWFREED_TEMPORARILY);

	* (int32_t *) ((char *) q - 32) &= ~GWFREED_TEMPORARILY;

	for (i = j = 0; i < gwdata->gwnum_free_count; i++)
		if (gwdata->gwnum_free[i] != q)
			gwdata->gwnum_free[j++] = gwdata->gwnum_free[i];
	gwdata->gwnum_free_count = j;
}

/* Free all user allocated numbers */

void gwfreeall (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	unsigned int i;
	char	*p;
	int32_t	freeable;

/* Go through the allocated list and free any user allocated gwnums that */
/* are freeable.  In other words, unless the user is using the BIGBUF */
/* kludge, free all possible memory. */

	gwdata->gwnum_free_count = 0;
	for (i = 0; i < gwdata->gwnum_alloc_count; i++) {
		if (gwdata->gwnum_alloc[i] == gwdata->GW_MODULUS_FFT) continue;
		if (gwdata->gwnum_alloc[i] == gwdata->GW_RECIP_FFT) continue;
		if (gwdata->gwnum_alloc[i] == gwdata->GW_RANDOM) continue;
		p = (char *) gwdata->gwnum_alloc[i];
		freeable = * (int32_t *) (p - 32) & ~GWFREED_TEMPORARILY;
		if (freeable) {
			aligned_free ((char *) p - GW_HEADER_SIZE);
			gwdata->gwnum_alloc[i--] = gwdata->gwnum_alloc[--gwdata->gwnum_alloc_count];
		}
		else
			gwdata->gwnum_free[gwdata->gwnum_free_count++] = gwdata->gwnum_alloc[i];
	}
}

/* To optimize use of the L1 cache we scramble the FFT data. */
/* Consult the assembly language code for better descriptions of this */
/* shuffling process.  This C code must accurately reflect the shuffling */
/* the assembly language code is expecting. */

unsigned long addr_offset (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long fftlen;
	unsigned long addr, i1, i2, i3, i6;

	fftlen = gwdata->FFTLEN;

/* Memory layouts for AVX-512 machines */

	if (gwdata->cpu_flags & CPU_AVX512F) {

/* Small FFTs use one pass, not very convoluted.  This example shows low-word/high-word pairs for a length 1280 FFT */
/*	0	8	16	24	32	40	48	56	640	+8	...	*/
/*	1	...										*/
/*	...											*/
/*	7											*/
/*	64											*/
/*	...											*/

		if (gwdata->PASS1_SIZE == 0) {
			unsigned long top_bit, over64;
			top_bit = i / (fftlen >> 1); i -= top_bit * (fftlen >> 1); // 640 in example above
			over64 = i / 64; i -= over64 * 64;
			addr = (over64 << 7) + ((i & 7) << 4) + (top_bit << 3) + (i >> 3);
			addr = addr * sizeof (double);
			/* Now optionally add 64 pad bytes every 4KB */
			if (gwdata->FOURKBGAPSIZE) addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}

/* Larger FFTs use two passes.  This example shows low-word/high-word pairs for a length 128K FFT (pass2_size = 1K complex or 2K reals) */
/*	0	+1K	+1K	+1K	+1K	+1K	+1K	+1K	64K	+1K	...	*/
/*	8	...										*/
/*	...											*/
/*	1016	...										*/
/*	1	...										*/
/*	...											*/
/*	1017	...										*/
/*	...											*/
/*	7	...										*/
/*	...											*/
/*	1023	...										*/
/*	8K	...										*/
/*	...											*/

		else {
			unsigned long top_bit, grp, row;
			top_bit = i / (fftlen >> 1); i -= top_bit * (fftlen >> 1); // 64K in example above
			i1 = i % gwdata->PASS2_SIZE;	// Element within pass 2, 0-1K in example above
			i = i / gwdata->PASS2_SIZE;	// Pass 1 block number, which 1K block in example above
			grp = (i >> 3) * 8 + (i1 & 7);	// Which 8K block plus which 1-7 block in example above
			row = grp * (gwdata->PASS2_SIZE >> 3) + (i1 >> 3); // Row number, where row is 128-bytes
			addr = (row << 4) + (top_bit << 3) + (i & 7);
			addr = addr * sizeof (double);
			/* Now optionally add bytes for every 4KB page and every pass 2 block */
			addr = addr + (addr >> 12) * gwdata->FOURKBGAPSIZE + grp * gwdata->PASS2GAPSIZE;
		}
	}

/* Memory layouts for AVX machines */

	else if (gwdata->cpu_flags & CPU_AVX) {

/* Small FFTs use one pass, not very convoluted.  This example is for	*/
/* a length 2048 FFT:							*/
/*	0	64	128	192	1024	+64	+64	+64	*/
/*	1	...							*/
/*	...								*/
/*	63								*/
/*	256								*/
/*	...								*/

		if (gwdata->PASS2_SIZE == 0) {
			unsigned long top5bits;
			top5bits = i / (fftlen >> 5); i -= top5bits * (fftlen >> 5);
			addr = ((top5bits >> 2) & 3) * (fftlen >> 2) + (i << 3) + ((top5bits >> 4) << 2) + (top5bits & 3);
			addr = addr * sizeof (double);
			/* Now optionally add 64 pad bytes every 1KB, 2KB or 4KB */
			if (gwdata->FOURKBGAPSIZE)
				addr = addr + (addr / (gwdata->FOURKBGAPSIZE << 6)) * 64;
		}

/* Larger FFTs use two passes.  This example for a length 64K FFT (pass2_size = 1K): */
/*	0	+1K	+1K	+1K	32K	+1K	+1K	+1K	*/
/*	4	...							*/
/*	...								*/
/*	1020	...							*/
/*	1	...							*/
/*	...								*/
/*	1021	...							*/
/*	2	...							*/
/*	...								*/
/*	1022	...							*/
/*	3	...							*/
/*	...								*/
/*	1023	...							*/
/*	4K	...							*/
/*	...								*/

		else {
			unsigned long top_bit, grps, row;
			top_bit = (i << 1) / fftlen; i -= (top_bit * fftlen) >> 1;
			i1 = i % gwdata->PASS2_SIZE; i = i / gwdata->PASS2_SIZE;
			grps = (i >> 2) * 4 + (i1 & 3);
			row = grps * (gwdata->PASS2_SIZE >> 2) + (i1 >> 2);
			addr = (row << 3) + (top_bit << 2) + (i & 3);
			addr = addr * sizeof (double);
			/* Now add 64 bytes every 4KB and one pass2gapsize */
			/* for every pass 2 block. */
			addr = addr + (addr >> 12) * gwdata->FOURKBGAPSIZE + grps * gwdata->PASS2GAPSIZE;
		}
	}

/* Memory layouts for SSE2 machines */

	else if (gwdata->cpu_flags & CPU_SSE2) {
		unsigned long sets, pfa, temp;

/* Small FFTs use one pass, not very convoluted.  This example is for	*/
/* a length 2048 FFT:							*/
/*	0	512	1	513	1024	1536	1025	1537	*/
/*	2	...							*/
/*	...								*/
/*	510								*/
/* PFA-style FFTs are a little tricker.  See assembly code for example.	*/

		if (gwdata->PASS2_SIZE == 0) {
			sets = fftlen >> 3;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			addr = (((((i3 << 1) + i6) << 1) + i1) << 1) + i;
			addr = addr * sizeof (double);
		}

/* Larger FFTs use two passes.  This example is for a length 64K FFT (pass2_size = 2K): */
/*	0	1K	16K	17K	32K	33K	48K	49K	*/
/*	1	...							*/
/*	...								*/
/*	1023	...							*/
/*	2K	...							*/
/*	...								*/
/* and PFA layouts are even funkier.					*/

		else if (gwdata->FFT_TYPE == FFT_TYPE_HOME_GROWN) {
			sets = (fftlen / gwdata->PASS2_SIZE) >> 2;
			if (i >= (fftlen >> 1)) {
				i6 = 1;
				i -= (fftlen >> 1);
			} else
				i6 = 0;
			i1 = i % (gwdata->PASS2_SIZE >> 1);
			i = i / (gwdata->PASS2_SIZE >> 1);
			i2 = i & 1; i >>= 1;
			i3 = 0;
			for (pfa = sets; pfa > 8; pfa >>= 1);
			if (pfa == 5) {
				temp = sets / 5;
				if (i < temp * 2) {
					sets = temp;
				} else {
					i3 = temp; i -= temp * 2;
					sets = temp * 4;
				}
			} else if (pfa == 7) {
				temp = sets / 7;
				if (i < temp * 2) {
					sets = temp;
				} else if (i < temp * 6) {
					i3 = temp; i -= temp * 2;
					sets = temp * 2;
				} else {
					i3 = temp * 3; i -= temp * 6;
					sets = temp * 4;
				}
			}
			i3 += i % sets; i /= sets;
			addr = i3 * (gwdata->PASS2_SIZE >> 1);
			addr = ((((((addr + i1) << 1) + i6) << 1) + i) << 1) + i2;
			addr = addr * sizeof (double);
			/* Now add 128 bytes every 8KB and one pass2gapsize */
			/* for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + i3 * gwdata->PASS2GAPSIZE;
		}

/* Newer traditional radix-4 large FFTs use don't have a special layout for PFA. */

		else {
			unsigned long top2, row;
			top2 = (i << 2) / fftlen; i -= (top2 * fftlen) >> 2;
			i1 = i % (gwdata->PASS2_SIZE >> 1); i = i / (gwdata->PASS2_SIZE >> 1);
			i2 = i & 1; i >>= 1;
			row = i * (gwdata->PASS2_SIZE >> 1) + i1;
			addr = (((row << 2) + top2) << 1) + i2;
			addr = addr * sizeof (double);
			/* Now add 128 bytes every 8KB and one pass2gapsize */
			/* for every pass 2 block. */
			addr = addr + (addr >> 13) * 128 + i * gwdata->PASS2GAPSIZE;
		}
	}

/* One pass x87 FFTs use a near flat memory model. */

	else if (gwdata->PASS2_SIZE == 0) {
		if (i >= (fftlen >> 1)) {
			i2 = 1;
			i -= (fftlen >> 1);
		} else
			i2 = 0;
		addr = i * 16 + i2 * 8;
	}

/* Two pass x87 FFTs use a near flat memory model.  Waste 64 bytes */
/* between 4KB.  Waste 64 bytes between every block (4KB, 16KB, or 64KB). */

	else {
		if (i >= (fftlen >> 1)) {
			i2 = 1;
			i -= (fftlen >> 1);
		} else
			i2 = 0;
		addr = i * 16 + i2 * 8 + (i >> 8) * 64 + (i / gwdata->PASS2_SIZE) * 64;
	}

/* Return the offset */

	return (addr);
}

/* Return the address of ith element in the FFT array */

double *addr (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	g,
	unsigned long i)
{
	return ((double *) ((char *) g + addr_offset (gwdata, i)));
}

/* Return the amount of data allocated by gwsetup */

unsigned long gwmemused (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	return (gwdata->mem_needed + gwdata->SCRATCH_SIZE);
}

/* Get the amount of memory required for the gwnum's raw FFT data.  This */
/* does not include the GW_HEADER_SIZE bytes for the header or any pad */
/* bytes that might be allocated for alignment.  I see little need for */
/* programs to use this routine. */

unsigned long gwnum_datasize (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	return (addr_offset (gwdata, gwdata->FFTLEN - 1) + sizeof (double));
}

/* Get the amount of memory likely to be allocated a gwnum.  This includes */
/* FFT data, headers, and pad bytes for alignment. */

unsigned long gwnum_size (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	unsigned long mem;
	mem = gwnum_datasize (gwdata) + GW_HEADER_SIZE;
	mem += gwdata->GW_ALIGNMENT-1;
	return (mem - mem % (gwdata->GW_ALIGNMENT));
}

/* Each FFT word is multiplied by a two-to-phi value.  These */
/* routines set and get the FFT value without the two-to-phi */
/* multiplier. */

int get_fft_value (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	g,
	unsigned long i,
	long	*retval)
{
	double	val;

	if (((uint32_t *) g)[-7] == 3) return (GWERROR_FFT); /* Test the FFTed flag */
	if (((uint32_t *) g)[-7] == 1) return (GWERROR_PARTIAL_FFT); /* Test the FFT-started flag */

/* Get the FFT data and validate it */

	val = * addr (gwdata, g, i);
	if (! is_valid_double (val)) return (GWERROR_BAD_FFT_DATA);

/* Rational and AVX-512 FFTs are not weighted */

	if (!gwdata->RATIONAL_FFT && !(gwdata->cpu_flags & CPU_AVX512F)) {

/* Handle r4dwpn FFTs which are only partially normalized */

		if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN)
			val = val * gwfft_partial_weight_inverse_sloppy (gwdata->dd_data, i, dwpn_col (gwdata, i));

/* Multiply by two-to-minus-phi to generate an integer. */

		else
			val = val * gwfft_weight_inverse_sloppy (gwdata->dd_data, i);
	}

/* Round the value to the nearest integer */

	if (val < -0.5)
		*retval = (long) (val - 0.5);
	else
		*retval = (long) (val + 0.5);

/* Return success */

	return (0);
}

void set_fft_value (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	g,
	unsigned long i,
	long	val)
{

/* Handle the rational and AVX-512 FFT case quickly (not weighted) */

	if (gwdata->RATIONAL_FFT || (gwdata->cpu_flags & CPU_AVX512F) || val == 0) {
		* addr (gwdata, g, i) = val;
		return;
	}

/* Handle r4dwpn FFTs which are only partially normalized */

	if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		* addr (gwdata, g, i) = val * gwfft_partial_weight_sloppy (gwdata->dd_data, i, dwpn_col (gwdata, i));
		return;
	}

/* Multiply by two-to-phi to generate the proper double. */

	* addr (gwdata, g, i) = val * gwfft_weight_sloppy (gwdata->dd_data, i);
}

/* Some words in the FFT data contain floor(p/N), some words contain */
/* floor(p/N)+1 bits.  This function returns TRUE in the latter case. */

int is_big_word (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long i)
{
	unsigned long base, next_base;

/* Compute the number of b in this word.  It is a big word if */
/* the number of b is more than NUM_B_PER_SMALL_WORD. */

	base = gwfft_base (gwdata->dd_data, i);
	next_base = gwfft_base (gwdata->dd_data, i+1);
	return ((next_base - base) > gwdata->NUM_B_PER_SMALL_WORD);
}

/* Routine map a "bit" number into an FFT word and "bit" within that word */
/* If b != 2, this routine locates the nth b amongst the FFT words. */

void bitaddr (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long bit,
	unsigned long *word,
	unsigned long *bit_in_word)
{

/* What word is the bit in? */

	*word = (unsigned long) ((double) bit / gwdata->avg_num_b_per_word);
	if (*word >= gwdata->FFTLEN) *word = gwdata->FFTLEN - 1;

/* Compute the bit within the word. */

	*bit_in_word = bit - gwfft_base (gwdata->dd_data, *word);
}

/* Map a gwerror code into human readable text */

void gwerror_text (
	gwhandle *gwdata,	/* Handle used in all gwnum calls */
	int	error_code,	/* Error code to turn ino text */
	char	*buf,		/* Buffer to write text to */
	int	buflen)		/* Sizeof the text buffer */
{
	char	localbuf[512];

/* Map error code to text */

	switch (error_code) {
	case GWERROR_VERSION:
		strcpy (localbuf, "Improperly compiled and linked.  Gwnum.h and FFT assembly code version numbers do not match.");
		break;
	case GWERROR_TOO_SMALL:
		strcpy (localbuf, "Number sent to gwsetup is less than or equal to one.");
		break;
	case GWERROR_TOO_LARGE:
		strcpy (localbuf, "Number sent to gwsetup is too large for the FFTs to handle.");
		break;
	case GWERROR_K_TOO_SMALL:
		strcpy (localbuf, "Value of k in k*b^n+c is too small.  Values less than one are not supported.");
		break;
	case GWERROR_K_TOO_LARGE:
		strcpy (localbuf, "Value of k in k*b^n+c is too large.  Values greater than 2251799813685247 are not supported.");
		break;
	case GWERROR_MALLOC:
		strcpy (localbuf, "Unable to allocate memory.  One possible cause is the operating system's swap area is too small.");
		break;
	case GWERROR_VERSION_MISMATCH:
		strcpy (localbuf, "GWNUM_VERSION from gwinit call doesn't match GWNUM_VERSION when gwnum.c was compiled.  Recompile and relink.");
		break;
	case GWERROR_STRUCT_SIZE_MISMATCH:
		strcpy (localbuf, "Gwhandle structure size from gwinit call doesn't match size when gwnum.c was compiled.  Check compiler alignment switches, recompile and relink.");
		break;
	default:
		if (error_code >= GWERROR_INTERNAL && error_code <= GWERROR_INTERNAL+100)
			sprintf (localbuf, "Internal error #%d.  Please contact the program's author.", error_code - GWERROR_INTERNAL);
		else
			sprintf (localbuf, "Unknown gwnum error code: %d", error_code);
		break;
	}

/* Copy error message to caller's buffer */

	if ((int) strlen (localbuf) >= buflen) localbuf[buflen-1] = 0;
	strcpy (buf, localbuf);
}

/* Return a description of the FFT type chosen */

void gwfft_description (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	char	*buf)		/* Buffer to return string in */
{
	char	*arch, *ffttype;

	arch = "";
	if (gwdata->cpu_flags & CPU_AVX512F) {
		// No output means Intel Skylake-X or Blend optimized
	}
	else if (gwdata->cpu_flags & CPU_AVX) {
		// No output means Intel Core2 or FMA3 or Blend optimized
	}
	else if (gwdata->cpu_flags & CPU_SSE2) {
		// No output means Intel Core2 or Blend optimized
		if (gwdata->ARCH == ARCH_P4 || gwdata->ARCH == ARCH_P4TP) arch = "Pentium4 ";
		if (gwdata->ARCH == ARCH_K8) arch = "AMD K8 ";
		if (gwdata->ARCH == ARCH_K10) arch = "AMD K10 ";
	}

	ffttype = "";
	if (gwdata->PASS2_SIZE) {
		if (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX)) {
			// No output means r4dwpn (type-3) FFT
		}
		else if (gwdata->cpu_flags & CPU_SSE2) {
			// No output means r4dwpn (type-3) FFT
			if (gwdata->FFT_TYPE == 0) ffttype = "type-0 ";
			if (gwdata->FFT_TYPE == 1) ffttype = "type-1 ";
			if (gwdata->FFT_TYPE == 2) ffttype = "type-2 ";
		}
	}

	sprintf (buf, "%s%s%s%sFFT length %lu%s",
		 gwdata->ALL_COMPLEX_FFT ? "all-complex " :
		 gwdata->ZERO_PADDED_FFT ? "zero-padded " :
		 gwdata->GENERAL_MOD ? "generic reduction " : "",
		 (gwdata->ARCH == ARCH_FMA3) ? "FMA3 " :
		 (gwdata->cpu_flags & CPU_AVX512F) ? "AVX-512 " :
		 (gwdata->cpu_flags & CPU_AVX) ? "AVX " :
		 (gwdata->cpu_flags & CPU_SSE2) ? "" : "x87 ",
		 arch, ffttype,
		 gwdata->FFTLEN >= 1048576 && (gwdata->FFTLEN & 0xFFFFF) == 0 ? gwdata->FFTLEN / 1048576 :
		 gwdata->FFTLEN >= 1024 && (gwdata->FFTLEN & 0x3FF) == 0 ? gwdata->FFTLEN / 1024 : gwdata->FFTLEN,
		 gwdata->FFTLEN >= 1048576 && (gwdata->FFTLEN & 0xFFFFF) == 0 ? "M" :
		 gwdata->FFTLEN >= 1024 && (gwdata->FFTLEN & 0x3FF) == 0 ? "K" : "");

	if (gwdata->PASS1_SIZE) {
		int	clm;
		char	p1buf[20], p2buf[20];

		if (gwdata->cpu_flags & CPU_AVX512F) clm = gwdata->PASS1_CACHE_LINES / 8;
		else if (gwdata->cpu_flags & CPU_AVX) clm = gwdata->PASS1_CACHE_LINES / 4;
		else clm = gwdata->PASS1_CACHE_LINES / 2;
		if (gwdata->PASS1_SIZE % 1024) sprintf (p1buf, "%d", (int) gwdata->PASS1_SIZE);
		else sprintf (p1buf, "%dK", (int) (gwdata->PASS1_SIZE / 1024));
		if (gwdata->PASS2_SIZE % 1024) sprintf (p2buf, "%d", (int) gwdata->PASS2_SIZE);
		else sprintf (p2buf, "%dK", (int) (gwdata->PASS2_SIZE / 1024));
		sprintf (buf + strlen (buf), ", Pass1=%s, Pass2=%s, clm=%d", p1buf, p2buf, clm);
	}

	if (gwdata->num_threads > 1)
		sprintf (buf + strlen (buf), ", %d threads", (int) gwdata->num_threads);

	if (gw_using_large_pages (gwdata))
		strcat (buf, " using large pages");
}

/* Return a string representation of a k/b/n/c combination */

void gw_as_string (
	char	*buf,		/* Buffer to return string in */
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c)		/* C in K*B^N+C */
{
	if (n == 0)
		sprintf (buf, "%.0f", k + c);
	else if (k != 1.0)
		sprintf (buf, "%.0f*%lu^%lu%c%lu", k, b, n, c < 0 ? '-' : '+', (unsigned long) labs (c));
	else if (b == 2 && c == -1)
		sprintf (buf, "M%lu", n);
	else {
		unsigned long cnt, temp_n;
		for (cnt = 0, temp_n = n; !(temp_n & 1); temp_n >>= 1, cnt++);
		if (b == 2 && temp_n == 1 && c == 1)
			sprintf (buf, "F%lu", cnt);
		else
			sprintf (buf, "%lu^%lu%c%lu", b, n, c < 0 ? '-' : '+', (unsigned long) labs (c));
	}
}

/* Get or clear the roundoff error.  Remember that if the roundoff error */
/* exceeds 0.5 then the FFT results will be wrong.  It is prudent to watch */
/* the roundoff error to make sure the roundoff error does not get close */
/* to 0.5. */

double gw_get_maxerr (
	gwhandle *gwdata)
{
	return (((struct gwasm_data *) gwdata->asm_data)->MAXERR);
}
void gw_clear_maxerr (
	gwhandle *gwdata)
{
	((struct gwasm_data *) gwdata->asm_data)->MAXERR = 0.0;
}

/* Return TRUE if we are operating near the limit of this FFT length */
/* Input argument is the percentage to consider as near the limit. */
/* For example, if percent is 0.1 and the FFT can handle 20 bits per word, */
/* then if there are more than 19.98 bits per word this function will */
/* return TRUE. */

int gwnear_fft_limit (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	pct)
{

/* Return TRUE if the virtual bits per word is near the maximum bits */
/* per word. */

	return (virtual_bits_per_word (gwdata) >
			(100.0 - pct) / 100.0 * gwdata->fft_max_bits_per_word);
}

/* Compute the virtual bits per word.  That is, the mersenne-mod-equivalent */
/* bits that this k,b,c combination uses.  This code must carefully invert */
/* the calculations gwinfo uses in determining whether a k,b,n,c combination */
/* will work for a given FFT size. */

double virtual_bits_per_word (
	gwhandle *gwdata)	/* Handle initialized by gwsetup */
{
	double	log2b, b_per_input_word, weighted_bits_per_output_word;
	double	max_weighted_bits_per_output_word;
	int	num_b_in_big_word, num_small_words, num_big_words;

	log2b = log2 (gwdata->b);

/* Compute our bits per output word exactly like gwinfo does for a zero padded FFT. */

	if (gwdata->ZERO_PADDED_FFT) {
		b_per_input_word = (double) (gwdata->n + gwdata->n) / gwdata->FFTLEN;
		num_b_in_big_word = (int) ceil (b_per_input_word);
		num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * (gwdata->FFTLEN / 2 + 4));
		num_big_words = (gwdata->FFTLEN / 2 + 4) - num_small_words;
		max_weighted_bits_per_output_word =
			2.0 * gwdata->fft_max_bits_per_word + 0.6 * log2 (gwdata->FFTLEN / 2 + 4);
		weighted_bits_per_output_word =
		       2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
		       0.6 * log2 (gwdata->FFTLEN / 2 + 4);
		if ((gwdata->n + gwdata->n) % gwdata->FFTLEN == 0)
			weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);
		else if (! is_pathological_distribution (num_big_words, num_small_words))
			weighted_bits_per_output_word -=
				((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
				 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
						  2.0 + (log2b - 6.0) / 6.0);
	}

/* Compute our bits per output word exactly like gwinfo does for a non-zero-padded FFT. */

	else {
		b_per_input_word = gwdata->avg_num_b_per_word;
		num_b_in_big_word = (int) ceil (b_per_input_word);
		num_small_words = (int) ((num_b_in_big_word - b_per_input_word) * gwdata->FFTLEN);
		num_big_words = gwdata->FFTLEN - num_small_words;
		max_weighted_bits_per_output_word =
			2.0 * gwdata->fft_max_bits_per_word + 0.6 * log2 (gwdata->FFTLEN);
		weighted_bits_per_output_word =
			2.0 * ((b_per_input_word + 1.0) * log2b - 1.0) +
			0.6 * log2 (gwdata->FFTLEN) +
			log2 (gwdata->k) + 1.7 * log2 (labs (gwdata->c));
		if (gwdata->k == 1.0 && gwdata->n % gwdata->FFTLEN == 0)
			weighted_bits_per_output_word -= ((log2b <= 4.0) ? log2b : 1.4 * log2b);
		else if (num_big_words == 1 && gwdata->k > 1.0)
			weighted_bits_per_output_word += log2b;
		else if (! is_pathological_distribution (num_big_words, num_small_words))
			weighted_bits_per_output_word -=
				((log2b <= 3.0) ? (log2b - 1.0) / 2.0 :
				 (log2b <= 6.0) ? 1.0 + (log2b - 3.0) / 3.0 :
						  2.0 + (log2b - 6.0) / 6.0);
	}

/* Now generate a value that can compared to gwdata->fft_max_bits_per_word */

	return (weighted_bits_per_output_word / max_weighted_bits_per_output_word * gwdata->fft_max_bits_per_word);
}

/* Given k,b,n,c determine the fft length.  If k,b,n,c is not supported */
/* then return zero.  Does not use benchmarking data. */

unsigned long gwmap_to_fftlen (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	return (gwmap_with_cpu_flags_to_fftlen (0, k, b, n, c));
}

unsigned long gwmap_with_cpu_flags_to_fftlen (
	int	cpu_flags,
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the FFT length */

	gwinit (&gwdata);
	if (cpu_flags) gwdata.cpu_flags = cpu_flags;
	gwclear_use_benchmarks (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (0);
	return (gwdata.jmptab->fftlen);
}

/* Given an fft length, determine the maximum allowable exponent.  If fftlen */
/* is not supported then return zero.  Does not use benchmarking data. */

unsigned long gwmap_fftlen_to_max_exponent (
	unsigned long fftlen)
{
	return (gwmap_with_cpu_flags_fftlen_to_max_exponent (0, fftlen));
}

unsigned long gwmap_with_cpu_flags_fftlen_to_max_exponent (
	int	cpu_flags,
	unsigned long fftlen)
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the maximum exponent for the FFT length */

	gwinit (&gwdata);
	if (cpu_flags) gwdata.cpu_flags = cpu_flags;
	gwclear_use_benchmarks (&gwdata);
	gwset_minimum_fftlen (&gwdata, fftlen);
	if (gwinfo (&gwdata, 1.0, 2, 0, -1)) return (0);
	return (adjusted_max_exponent (&gwdata, gwdata.jmptab));
}

/* Given an fft length, determine how much memory is used for normalization */
/* and sin/cos tables.  If k,b,n,c is not supported, then kludgily return */
/* 100 million bytes used.  Does not use benchmarking data. */

unsigned long gwmap_to_memused (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the memory used */

	gwinit (&gwdata);
	gwclear_use_benchmarks (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (100000000L);
	return (gwdata.mem_needed + gwdata.SCRATCH_SIZE);
}

/* Return the estimated size of a gwnum.  Does not use benchmarking data. */

unsigned long gwmap_to_estimated_size (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */

/* Get pointer to fft info and return the memory used */

	gwinit (&gwdata);
	gwclear_use_benchmarks (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (100000000L);
	return (addr_offset (&gwdata, gwdata.FFTLEN - 1) + sizeof (double));
}

/* Speed of other x87 processors compared to a Pentium II */

#define REL_486_SPEED	8.4	/* 486 is over 8 times slower than PII */
#define REL_K6_SPEED	3.0	/* K6 is 3 times slower than PII */
#define REL_P3_SPEED	0.8	/* Pentium III is 20% faster than PII */
#define REL_K7_SPEED	0.6	/* Athlons are much faster than a PII */

/* Speed of other SSE2 processors compared to a Pentium 4 */

#define REL_AMD64_SPEED	1.1	/* AMD64 is slightly slower than a P4 */
#define REL_PM_SPEED	1.4	/* Pentium M, Core are much slower than a P4 */
#define REL_ATOM_SPEED	5.0	/* Atoms are much, much slower than a P4 */
#define REL_CORE2_SPEED	0.625	/* Core 2 is much faster than a P4 */
#define REL_I7_SPEED	0.59	/* Core i7 is even faster than a Core 2 */
#define REL_PHENOM_SPEED 0.67	/* AMD Phenom is faster that a P4 */
#define REL_FMA3_SPEED	0.95	/* Haswell FMA3 is faster than a Ivy/Sandy Bridge AVX CPU */

/* Speed of other AVX processors compared to a Sandy Bridge */

#define REL_BULLDOZER_SPEED	1.9	/* Bulldozer is slower than Sandy Bridge */
#define REL_ZEN_SPEED		1.4	/* Zen is likely slower than Sandy Bridge which has true 256-bit AVX support whereas Zen has FMA3 */

/* Make a guess as to how long a squaring will take.  If the number cannot */
/* be handled, then kludgily return 100.0.  Does not use benchmarking data. */

double gwmap_to_timing (
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
{
	gwhandle gwdata;	/* Temporary gwnum handle */
	double	timing;

/* Get pointer to fft info */

	gwinit (&gwdata);
	gwclear_use_benchmarks (&gwdata);
	if (gwinfo (&gwdata, k, b, n, c)) return (100.0);

/* Use my PII-400 or P4-1400 timings as a guide. */

	timing = gwdata.jmptab->timing;

/* Since the program is about 10% memory bound, the program will not */
/* speed up linearly with increase in chip speed.  Note, no attempt is */
/* made to differentiate between memory bus speed - we're */
/* just returning an educated guess here. */

/* Adjust timing for various CPU architectures. */
/* For Intel, 486s were very slow.  Pentium, Pentium Pro, Pentium II, */
/* and old celerons were slow because they did not support prefetch. */
/* AMD64s and Pentium Ms are slower than P4s. */

	if (gwdata.cpu_flags & CPU_AVX512F) {
		timing = 0.10 * timing + 0.90 * timing * 4100.0 / CPU_SPEED;	/* Calibrated for Skylake */	//BUG - not timed yet
		if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_OTHER) timing *= REL_ZEN_SPEED;  /* Complete guess for future AMD CPUs */
	} else if (gwdata.cpu_flags & CPU_AVX) {
		timing = 0.10 * timing + 0.90 * timing * 4100.0 / CPU_SPEED;	/* Calibrated for Sandy Bridge */
		if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_OTHER) timing *= REL_ZEN_SPEED;  /* Complete guess for future AMD CPUs */
		else if (CPU_ARCHITECTURE == CPU_ARCHITECTURE_AMD_ZEN) timing *= REL_ZEN_SPEED;
		else if (strstr (CPU_BRAND, "AMD")) timing *= REL_BULLDOZER_SPEED;
		if (gwdata.cpu_flags & CPU_FMA3) timing *= REL_FMA3_SPEED;
	} else if (gwdata.cpu_flags & CPU_SSE2) {
		timing = 0.10 * timing + 0.90 * timing * 1400.0 / CPU_SPEED;
		if (strstr (CPU_BRAND, "Phenom")) timing *= REL_PHENOM_SPEED;
		else if (strstr (CPU_BRAND, "AMD")) timing *= REL_AMD64_SPEED;
		else if (strstr (CPU_BRAND, "Atom")) timing *= REL_ATOM_SPEED;
		else if (strstr (CPU_BRAND, "Core 2")) timing *= REL_CORE2_SPEED;
		else if (strstr (CPU_BRAND, "Core(TM)2")) timing *= REL_CORE2_SPEED;
		else if (strstr (CPU_BRAND, "Core(TM) i7")) timing *= REL_I7_SPEED;
		else if (strstr (CPU_BRAND, "Pentium(R) M")) timing *= REL_PM_SPEED;
		else if (strstr (CPU_BRAND, "Core")) timing *= REL_PM_SPEED;
	} else {
		timing = 0.10 * timing + 0.90 * timing * 400.0 / CPU_SPEED;
		if (strstr (CPU_BRAND, "486")) timing *= REL_486_SPEED;
		else if (strstr (CPU_BRAND, "Intel")) {
			if (gwdata.cpu_flags & CPU_PREFETCH) timing *= REL_P3_SPEED;
		} else if (strstr (CPU_BRAND, "AMD")) {
			if (strstr (CPU_BRAND, "Unknown")) timing *= REL_486_SPEED;
			else if (strstr (CPU_BRAND, "K5")) timing *= REL_486_SPEED;
			else if (strstr (CPU_BRAND, "K6")) timing *= REL_K6_SPEED;
			else timing *= REL_K7_SPEED;
		} else
			timing *= REL_486_SPEED;
	}
	return (timing);
}

/* Given k,b,n,c determine the fft length and zero-padding state to be */
/* used.  Caller peers into the gwdata structure to get this info. */
// DEPRECATED
//int gwmap_to_fft_info (
//	gwhandle *gwdata,	/* Uninitialized gwnum global data */
//	double	k,		/* K in K*B^N+C. Must be a positive integer. */
//	unsigned long b,	/* B in K*B^N+C. */
//	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
//	signed long c)		/* C in K*B^N+C. Must be rel. prime to K. */
//{
//	gwinit (gwdata);
//	return (gwinfo (gwdata, k, b, n, c));
//}


/* Internal routine to help gwcopyzero */

void calc8ptrs (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long n,	/* Number of words to zero */
	uint32_t *ptrs)
{
	unsigned long i, j, k;

/* AVX-512F processes one double cache line at a time.  The one-pass and two-pass mameory layouts are different. */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		int	block_size, upper_avx512_word;

/* First the one-pass case.  The copyzero array contains masks and counts */
/* for when to load new masks.  Here is an example for a 2K FFT (layout below) */
/*	0	+8	+8	+8	+8	+8	+8	+8	1K	+8	...	*/
/*	1	...										*/
/*	...											*/
/*	7	...										*/
/*	64	...										*/
/*	...											*/
/* The assembly code processes 1 block at a time.  Block 0 processes 0 to 63 and 1K to 1K+63. */
/* Block 1 processes 64 to 127 and 1K+64 to 1K+127, etc.  If n = 64 + 18, then the masks in block 0 are */
/* 00000000, 11111111, count = inf (zero low words, copy high words).  The masks in block 1 are */
/* 11111000, 11111111, count = 2, 11111100, 11111111, count = inf.  The masks in blocks 2+ are */
/* 11111111, 11111111, count = inf (copy low and high words). */

/* The two-pass case is a little more complex. Here is an example for a 128K FFT (layout below) */
/*	0	+1K	+1K	+1K	+1K	+1K	+1K	+1K	64K	+1K	...	*/
/*	8	...										*/
/*	...											*/
/*	1016	...										*/
/*	1	...										*/
/*	...											*/
/*	1017	...										*/
/*	...											*/
/*	7	...										*/
/*	...											*/
/*	1023	...										*/
/*	8K	...										*/
/* The assembly code processes 1 block at a time.  Block 0 processes 0 to 8K-1 and 64K to 72K-1. */
/* Block 1 processes 8K to 16K-1 and 72K to 80K-1, etc.  If n = 9K + 20, then the masks in block 0 are */
/* 00000000, 11111111, count = inf (zero low words, copy high words).  The masks in block 1 are */
/* 11111100, 11111111, count = 20, 11111110, 11111111, count = inf.  The masks in blocks 2+ are */
/* 11111111, 11111111, count = inf (copy low and high words). */

		if (gwdata->PASS1_SIZE == 0) {
			upper_avx512_word = 8;
			block_size = 64;
		} else {
			upper_avx512_word = gwdata->PASS2_SIZE;
			block_size = upper_avx512_word * 8;
		}

		// Create masks assuming n is below the half way point, we correct for this later if not true
		// Masks and count for the "processing block before n" case.
		ptrs[0] = 0xFF00;
		ptrs[1] = 0xFFFFFFFF;
		// Masks and count for the "processing block containing n" case.
		i = n % block_size;				// Offset within block (0 to 63 or 8K-1 in examples above)
		j = i / upper_avx512_word;			// Column number (which multiple of 8 or 1K in examples above)
		ptrs[2] = 0xFFFF << (j+1);			// Zero j+1 columns mask
		ptrs[3] = i - j * upper_avx512_word;		// Count
		ptrs[4] = 0xFFFF << j;				// Zero j columns mask
		ptrs[5] = 0xFFFFFFFF;				// Count
		// Masks and count for the "processing block after n" case.
		ptrs[6] = 0xFFFF;
		ptrs[7] = 0xFFFFFFFF;
		/* If n > half way point, change masks to always zero low words */
		if (n >= gwdata->FFTLEN/2) ptrs[0] <<= 8, ptrs[2] <<= 8, ptrs[4] <<= 8, ptrs[6] <<= 8;
		/* Assembly code does not like first count of zero, compensate */
		if (ptrs[3] == 0) ptrs[2] = ptrs[4], ptrs[3] = ptrs[5];
	}

/* This is a grossly inefficient way to do this.  However, it should be called rarely. */

/* AVX offsets do not always increase, use a new method.  Remember the */
/* offset of the first word in each column that doesn't need zeroing. */

	else if (gwdata->cpu_flags & CPU_AVX) {
		for (i = 0; i < 8; i++) ptrs[i] = 0xFFFFFFFF;	/* Start by clearing the eight ptrs. */
		for (i = n; i < gwdata->FFTLEN; i++) {
			j = addr_offset (gwdata, i);
			k = (j & 63) >> 3;
			if (ptrs[k] == 0xFFFFFFFF) ptrs[k] = j - (k << 3);
		}
	}

/* The original version that works for offsets that are always increasing */

	else {
		for (i = 0; i < 8; i++) ptrs[i] = 0;	/* Start by clearing the eight ptrs. */
		for (i = 0; i < n; i++) {
			j = addr_offset (gwdata, i);
			k = (j & 63) >> 3;
			if (j >= ptrs[k]) ptrs[k] = j - (k << 3) + 64;
		}
	}
}

/* Routine that sets up and calls assembly code to copy a gwnum from */
/* source to dest while zeroing some lower FFT words */

void gwcopyzero (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,
	gwnum	d,
	unsigned long n)
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s)[-7] != 3);	// Number not FFTed
	ASSERTG (((uint32_t *) s)[-7] != 1);	// Number not partially FFTed

/* Handle case where no words are zeroed.  Some of the assembly routines */
/* do not like a word count of zero. */

	if (n == 0) {
		gwcopy (gwdata, s, d);
		return;
	}

/* Calculate pointers for ASM routines to do there zero-copy */

	if (n != gwdata->saved_copyz_n) {
		if (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2))
			calc8ptrs (gwdata, n, asm_data->COPYZERO);
		gwdata->saved_copyz_n = n;
	}

/* Do an AVX-512 FFT copyzero one block at a time */
/* This is in preparation for multi-thread support */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		//BUG - multithread this?? support AVX too?? change asm to process more than 16 doubles at a time?
		int	i, num_blks, blk_containing_n;

		if (gwdata->PASS1_SIZE == 0)
			blk_containing_n = (n % (gwdata->FFTLEN / 2)) / 64;
		else
			blk_containing_n = (n % (gwdata->FFTLEN / 2)) / (gwdata->PASS2_SIZE * 8);
		num_blks = asm_data->addcount1;
		for (i = 0; i < num_blks; i++) {
			asm_data->SRCARG = (char *) s + i * asm_data->pass1blkdst + (i/4) * asm_data->normblkdst4;
			asm_data->DESTARG = (char *) d + i * asm_data->pass1blkdst + (i/4) * asm_data->normblkdst4;
			if (i < blk_containing_n) asm_data->SRC2ARG = (char *) asm_data->COPYZERO;		// This block is before n
			else if (i == blk_containing_n) asm_data->SRC2ARG = (char *) asm_data->COPYZERO + 8;	// This block contains n
			else asm_data->SRC2ARG = (char *) asm_data->COPYZERO + 24;				// This block is past n
			gw_copyzero (gwdata, asm_data);
		}
	}

/* Do the copyzero the old-fashioned way -- in one call */
/* Call assembly language copy-and-zeroing routine */

	else {
		asm_data->SRCARG = s;
		asm_data->DESTARG = d;
		asm_data->NUMARG = n;
		gw_copyzero (gwdata, asm_data);
	}

/* Copy the unnormalized add counter and clear the has been completely/partially FFTed flag. */

	((uint32_t *) d)[-1] = ((uint32_t *) s)[-1];
	((uint32_t *) d)[-7] = 0;
}

/* Set the constant which the results of a multiplication should be */
/* multiplied by.  Use this macro in conjunction with the c argument of */
/* gwsetnormroutine. */

void gwsetmulbyconst (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	long	val)
{
	struct gwasm_data *asm_data;
	double	ktimesval, big_word;

/* Perform common computations */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	ktimesval = gwdata->k * (double) val;
	big_word = pow ((double) gwdata->b, gwdata->NUM_B_PER_SMALL_WORD + 1);

/* Adjust the AVX-512 assembly language constants affected by mulbyconst */

	if (gwdata->cpu_flags & CPU_AVX512F) {

/* Save mulbyconst and -c * mulbyconst as a double */

		asm_data->u.zmm.ZMM_MULCONST = (double) val;
		asm_data->u.zmm.ZMM_MINUS_C_TIMES_MULCONST = (double) -gwdata->c * (double) val;

/* Split k*mulconst for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k*mulconst if k*mulconst * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		if (fabs (ktimesval) * big_word < 562949953421312.0)
			asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE = 0.0;
		else if (ktimesval > 0.0)
			asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE = floor (ktimesval / big_word);
		else
			asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE = -floor (-ktimesval / big_word);
		asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_SMALL_BASE = asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE * (double) gwdata->b;
		asm_data->u.zmm.ZMM_K_TIMES_MULCONST_LO = ktimesval - asm_data->u.zmm.ZMM_K_TIMES_MULCONST_HI_OVER_LARGE_BASE * big_word;
	}

/* Adjust the AVX assembly language constants affected by mulbyconst */

	else if (gwdata->cpu_flags & CPU_AVX) {

/* Save mulbyconst and -c * mulbyconst as a double */

		asm_data->u.ymm.YMM_MULCONST[0] = asm_data->u.ymm.YMM_MULCONST[1] =
		asm_data->u.ymm.YMM_MULCONST[2] = asm_data->u.ymm.YMM_MULCONST[3] = (double) val;
		asm_data->u.ymm.YMM_MINUS_C_TIMES_MULCONST[0] = asm_data->u.ymm.YMM_MINUS_C_TIMES_MULCONST[1] =
		asm_data->u.ymm.YMM_MINUS_C_TIMES_MULCONST[2] = asm_data->u.ymm.YMM_MINUS_C_TIMES_MULCONST[3] = (double) -gwdata->c * (double) val;

/* Split k*mulconst for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k*mulconst if k*mulconst * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		if (fabs (ktimesval) * big_word < 562949953421312.0)
			asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[0] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[1] =
			asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[2] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[3] = 0.0;
		else if (ktimesval > 0.0)
			asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[0] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[1] =
			asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[2] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[3] =
				floor (ktimesval / big_word) * big_word;
		else
			asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[0] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[1] =
			asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[2] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[3] =
				-floor (-ktimesval / big_word) * big_word;
		asm_data->u.ymm.YMM_K_TIMES_MULCONST_LO[0] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_LO[1] =
		asm_data->u.ymm.YMM_K_TIMES_MULCONST_LO[2] = asm_data->u.ymm.YMM_K_TIMES_MULCONST_LO[3] =
			ktimesval - asm_data->u.ymm.YMM_K_TIMES_MULCONST_HI[0];
	}

/* Adjust the SSE2 assembly language constants affected by mulbyconst */

	else if (gwdata->cpu_flags & CPU_SSE2) {

/* Save mulbyconst and -c * mulbyconst as a double */

		asm_data->u.xmm.XMM_MULCONST[0] = asm_data->u.xmm.XMM_MULCONST[1] = (double) val;
		asm_data->u.xmm.XMM_MINUS_C_TIMES_MULCONST[0] = asm_data->u.xmm.XMM_MINUS_C_TIMES_MULCONST[1] = (double) -gwdata->c * (double) val;

/* Split k*mulconst for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k*mulconst if k*mulconst * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		if (fabs (ktimesval) * big_word < 562949953421312.0)
			asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[0] = asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[1] = 0.0;
		else if (ktimesval > 0.0)
			asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[0] = asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[1] =
				floor (ktimesval / big_word) * big_word;
		else
			asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[0] = asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[1] =
				-floor (-ktimesval / big_word) * big_word;
		asm_data->u.xmm.XMM_K_TIMES_MULCONST_LO[0] =
		asm_data->u.xmm.XMM_K_TIMES_MULCONST_LO[1] = ktimesval - asm_data->u.xmm.XMM_K_TIMES_MULCONST_HI[0];
	}

/* Adjust the x87 assembly language constants affected by mulbyconst */

	else {

/* Save mulbyconst and -c * mulbyconst as a double */

		asm_data->u.x87.MULCONST = (double) val;
		asm_data->u.x87.MINUS_C_TIMES_MULCONST = (double) -gwdata->c * (double) val;

/* Split k*mulconst for zero-padded FFTs emulating modulo k*b^n+c.  Note we only need */
/* to split k*mulconst if k*mulconst * big_word exceeds what a floating point register can hold. */
/* 2^49 should give us enough room to handle a carry and still round properly. */	

		asm_data->u.x87.K_TIMES_MULCONST_HI = floor (ktimesval / big_word) * big_word;
		asm_data->u.x87.K_TIMES_MULCONST_LO = ktimesval - asm_data->u.x87.K_TIMES_MULCONST_HI;
		big_word = big_word * big_word;
		asm_data->u.x87.K_TIMES_MULCONST_HI_2 = floor (asm_data->u.x87.K_TIMES_MULCONST_HI / big_word) * big_word;
		asm_data->u.x87.K_TIMES_MULCONST_HI_1 = asm_data->u.x87.K_TIMES_MULCONST_HI - asm_data->u.x87.K_TIMES_MULCONST_HI_2;
	}
}

/* Set the flag which controls whether the multiply code should begin the */
/* forward FFT of the results of a multiply. This is a little faster than */
/* than doing a full forward FFT later on.  The downside is the caller */
/* cannot convert the results of the multiply to an integer. */

void gwstartnextfft (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	int	state)		/* New value for POSTFFT */
{
	/* Option is not supported for generic modular reduction and one-pass FFTs */
	if (!gwdata->GENERAL_MOD && gwdata->PASS1_SIZE) {
		gwdata->POSTFFT = state;
		if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2))) {	// Copy state to asm_data for x87 FFTs
			struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
			asm_data->u.x87.POSTFFT = state;
		}
	}
}

/* Add a small constant at the specified power of b after the */
/* next multiplication.  That is, value*b^power_of_b is added to */
/* the next multiplication result.  This only works if k=1. */

void gwsetaddinatpowerofb (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	long	value,
	unsigned long power_of_b)
{
	unsigned long word, b_in_word;

	ASSERTG (gwdata->k == 1.0);

/* Tell assembly code to add the shifted value to the multiplication result. */

	bitaddr (gwdata, power_of_b, &word, &b_in_word);
	raw_gwsetaddin (gwdata, word, value * pow ((double) gwdata->b, b_in_word));
}

/* Routine that tells the assembly code to add a small value to the */
/* results of each multiply. */

void gwsetaddin (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	long	value)
{
	unsigned long word, b_in_word;

	ASSERTG (gwdata->k == 1.0 || labs (gwdata->c) == 1);

/* In a zero-padded FFT, the value is added into ZPAD0 */

	if (gwdata->ZERO_PADDED_FFT) {
		((struct gwasm_data *) gwdata->asm_data)->ADDIN_VALUE = (double) value;
		return;
	}

/* If k is 1, add the value into the first FFT word */

	if (gwdata->k == 1.0) {
		raw_gwsetaddin (gwdata, 0, value);
		return;
	}

/* If value is a multiple of b, "shift" it right and increment b count.  This */
/* will ensure that we modify the proper FFT word. */

	for (b_in_word = 0; value && value % (long) gwdata->b == 0; value = value / (long) gwdata->b)
		b_in_word++;

/* Convert the input value to 1/k format.  Case 1 (k*b^n-1): Inverse of k */
/* is b^n.  Case 3 (k*b^n+1): Inverse of k is -b^n.  No other cases can */
/* be handled. */

	if (gwdata->c == -1) {
		bitaddr (gwdata, gwdata->n + b_in_word, &word, &b_in_word);
	}
	else if (gwdata->c == 1) {
		bitaddr (gwdata, gwdata->n + b_in_word, &word, &b_in_word);
		value = -value;
	}

/* Tell assembly code to add the shifted value to the multiplication result. */

	raw_gwsetaddin (gwdata, word, value * pow ((double) gwdata->b, b_in_word));
}

/* Routine that tells the assembly code to add a small value to the results of each multiply */

void raw_gwsetaddin (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	unsigned long word,
	double	val)
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	unsigned long row;

/* Handle calculations for AVX-512 FFTs.  All two-pass FFTs use a scratch area. */

	if (gwdata->cpu_flags & CPU_AVX512F) {

/* For one pass FFTs, use the offset to the FFT data value */

		if (gwdata->PASS1_SIZE == 0) {
			asm_data->ADDIN_ROW = 0;
			asm_data->ADDIN_OFFSET = addr_offset (gwdata, word);
		}

/* For two-pass FFTs, compute the first FFT value in the pass 1 block (ADDIN_ROW) */
/* and the offset when the target data is in the scratch area. */

		else {
			row = word % gwdata->PASS2_SIZE;
			asm_data->ADDIN_ROW = row & ~(gwdata->PASS1_CACHE_LINES - 1);
			asm_data->ADDIN_OFFSET = (row - asm_data->ADDIN_ROW) * 128;
			if (word >= gwdata->FFTLEN / 2) asm_data->ADDIN_OFFSET += 64;
			asm_data->ADDIN_OFFSET += ((word / gwdata->PASS2_SIZE) & 7) * 8;
			row = (word % (gwdata->FFTLEN / 2)) / (8 * gwdata->PASS2_SIZE);
			asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 128 + (uint32_t) asm_data->normblkdst);
			asm_data->ADDIN_OFFSET += (row >> 3) * (uint32_t) asm_data->normblkdst8;
		}
	}

/* Handle calculations for AVX FFTs.  All two-pass FFTs use a scratch area. */

	else if (gwdata->cpu_flags & CPU_AVX) {

/* For two-pass FFTs, compute the first FFT value in the pass 1 block (ADDIN_ROW) */
/* and the offset when the target data is in the scratch area. */

		if (gwdata->PASS2_SIZE) {
			row = word % gwdata->PASS2_SIZE;
			asm_data->ADDIN_ROW = row & ~(gwdata->PASS1_CACHE_LINES - 1);
			asm_data->ADDIN_OFFSET = (row - asm_data->ADDIN_ROW) * 64;
			if (word >= gwdata->FFTLEN / 2) asm_data->ADDIN_OFFSET += 32;
			asm_data->ADDIN_OFFSET += ((word / gwdata->PASS2_SIZE) & 3) * 8;
			row = (word % (gwdata->FFTLEN / 2)) / (4 * gwdata->PASS2_SIZE);
			asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 64 + (uint32_t) asm_data->normblkdst);
			asm_data->ADDIN_OFFSET += (row >> 3) * (uint32_t) asm_data->normblkdst8;
		}

/* For one pass FFTs, use the offset to the FFT data value */

		else
			asm_data->ADDIN_OFFSET = addr_offset (gwdata, word);

	}

/* Handle calculations for SSE2 FFTs */

	else if (gwdata->cpu_flags & CPU_SSE2) {

/* Compute the offset to the FFT data value */

		asm_data->ADDIN_OFFSET = addr_offset (gwdata, word);

/* If this is a one-pass SSE2 FFT, then we need to tell the assembly code */
/* the affected "row", that is which set of data the add-in will */
/* take place */

		if (gwdata->PASS2_SIZE == 0) {
			row = asm_data->ADDIN_OFFSET & 31;
			if (row == 8) asm_data->ADDIN_OFFSET += 8;
			if (row == 16) asm_data->ADDIN_OFFSET -= 8;
		}

/* Factor in the blkdst value in xfft3.mac to compute the two pass */
/* SSE2 addin_offset. */

		else {
			unsigned long num_rows;

			num_rows = gwdata->PASS2_SIZE >> 1;
			row = word % num_rows;
			asm_data->ADDIN_ROW =
					row & ~(gwdata->PASS1_CACHE_LINES - 1);
			asm_data->ADDIN_OFFSET -= (row >> 7) * 128 +
					(row / gwdata->PASS1_CACHE_LINES) *
					gwdata->PASS1_CACHE_LINES * 64;

/* This case is particularly nasty as we have to convert the FFT data offset */
/* into a scratch area offset.  In assembly language terms, this means */
/* subtracting out multiples of blkdst and adding in multiples of clmblkdst */
/* and clmblkdst8. */

			if (gwdata->SCRATCH_SIZE) {
				unsigned long blkdst;

				blkdst = addr_offset (gwdata, gwdata->PASS2_SIZE);
				row = asm_data->ADDIN_OFFSET / blkdst;
				asm_data->ADDIN_OFFSET -= row * blkdst;
				asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 64);
				asm_data->ADDIN_OFFSET += (row >> 3) * (uint32_t) asm_data->normblkdst8;
			}
		}
	}

/* Handle calculations for x87 FFTs */

	else {

/* Compute the offset to the FFT data value */

		asm_data->ADDIN_OFFSET = addr_offset (gwdata, word);

/* x87 FFTs use a scratch area.  Like the SSE2 code */
/* we have to convert the FFT data offsets for two-pass FFTs. */

		if (gwdata->PASS2_SIZE) {
			unsigned long num_cache_lines, cache_line;

			num_cache_lines = gwdata->PASS2_SIZE >> 1;
			cache_line = ((word >> 1) & (num_cache_lines - 1));

			asm_data->ADDIN_ROW = ((num_cache_lines>>7) - (cache_line>>7)) * 65536 +
					(128 / gwdata->PASS1_CACHE_LINES -
					(cache_line & 127) / gwdata->PASS1_CACHE_LINES) * 256;
			asm_data->ADDIN_OFFSET -= (cache_line >> 7) * 64 +
					(cache_line / gwdata->PASS1_CACHE_LINES) *
					gwdata->PASS1_CACHE_LINES * 32;

/* This case is particularly nasty as we have to convert the FFT data offset */
/* into a scratch area offset.  In assembly language terms, this means */
/* subtracting out multiples of blkdst and adding in multiples of clmblkdst */
/* and clmblkdst32. */

			if (gwdata->SCRATCH_SIZE) {
				unsigned long blkdst;

				blkdst = addr_offset (gwdata, gwdata->PASS2_SIZE);
				row = asm_data->ADDIN_OFFSET / blkdst;
				asm_data->ADDIN_OFFSET -= row * blkdst;
				asm_data->ADDIN_OFFSET += row * (gwdata->PASS1_CACHE_LINES * 32);

/* Handle the FFTs where clmblkdst32 is used */

				if (((gwdata->FFTLEN / gwdata->PASS2_SIZE) >> 1) >= 128)
					asm_data->ADDIN_OFFSET += (row >> 5) * 64;
			}
		}
	}

/* Handle AVX-512 FFTs.  They are normalized differently than FMA3 and earlier FFTs */

	if (gwdata->cpu_flags & CPU_AVX512F) {

/* Rational FFTs are not normalized, don't need to apply a weight to the addin value. */
/* Exceptions are two-pass all-complex FFTs which have a sine component as part of the inverse group multiplier. */

		if (gwdata->RATIONAL_FFT && !(gwdata->PASS1_SIZE && gwdata->ALL_COMPLEX_FFT)) {
			asm_data->ADDIN_VALUE = val;
		}

/* AVX-512 FFT's ADDIN_VALUE is partially normalized (since we moved the weighting from the last_unfft macro to the normalize */
/* code to take advantage of an FMA opportunity.  Apply a partial weight to the addin value.  Also, if it's an all-complex */
/* FFT then we need to adjust for the roots-of-minus-one sine optimization too. */

		else {
			unsigned long upper_avx512_word, grp;
			double	ttmp;
			// Mimic code from building inverse multipliers in zr4dwpn_build_norm_table
			upper_avx512_word = (gwdata->PASS1_SIZE == 0) ? 8 : gwdata->PASS2_SIZE;
			grp = word - word % upper_avx512_word;
			if (gwdata->PASS1_SIZE && gwdata->ALL_COMPLEX_FFT) {
				unsigned long sine_word;
				sine_word = grp % (gwdata->FFTLEN / 2);
				gwfft_weights_times_sine (gwdata->dd_data, grp, sine_word, gwdata->FFTLEN * 2, NULL, &ttmp);
			} else
				ttmp = gwfft_weight_inverse (gwdata->dd_data, grp);
			asm_data->ADDIN_VALUE = val / ttmp;
		}
	}

/* Handle r4dwpn FFTs which are only partially normalized */

	else if (gwdata->FFT_TYPE == FFT_TYPE_RADIX_4_DWPN) {
		asm_data->ADDIN_VALUE = val * gwfft_partial_weight_sloppy (gwdata->dd_data, word, dwpn_col (gwdata, word));
	}

/* Set the addin value - multiply it by two-to-phi and FFTLEN/2/k. */

	else
		asm_data->ADDIN_VALUE = val * gwfft_weight_sloppy (gwdata->dd_data, word) * gwdata->FFTLEN * 0.5 / gwdata->k;
}


/********************************************************/
/* Routines to convert between gwnums and other formats */
/********************************************************/

/* Convert a double to a gwnum */

void dbltogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	d,		/* Input number */
	gwnum	g)		/* Gwnum value to set */
{
	stackgiant(tmp, 2);

	dbltog (d, tmp);
	gianttogw (gwdata, tmp, g);
}

/* Convert a binary value to a gwnum */

void binarytogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const uint32_t *array,	/* Array containing the binary value */
	uint32_t arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = arraylen;
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (gwdata, &tmp, n);
}

/* Convert a binary value to a gwnum */

void binary64togw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const uint64_t *array,	/* Array containing the binary value */
	uint64_t arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = (int) arraylen * sizeof (uint64_t) / sizeof (uint32_t);
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (gwdata, &tmp, n);
}

/* Convert a binary value to a gwnum */

void binarylongstogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const unsigned long *array, /* Array containing the binary value */
	unsigned long arraylen,	/* Length of the array */
	gwnum	n)		/* Destination gwnum */
{
	giantstruct tmp;
	tmp.sign = arraylen * sizeof (unsigned long) / sizeof (uint32_t);
	tmp.n = (uint32_t *) array;
	while (tmp.sign && tmp.n[tmp.sign-1] == 0) tmp.sign--;
	gianttogw (gwdata, &tmp, n);
}

/* Convert a giant to gwnum FFT format.  Giant must be a positive number. */

void gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	a,
	gwnum	g)
{
	ASSERTG (a->sign >= 0);		/* We only handle positive numbers */

/* Jean Penne requested that we optimize the small number cases. */
/* Setting the gwnum to zero is real easy. */

	if (a->sign == 0) {
		memset (g, 0, gwnum_datasize (gwdata));
	}

/* Small numbers can also be optimized for many moduli by zeroing all the */
/* FFT data using memset and then setting only the affected FFT elements. */

	else if (a->sign == 1 && (gwdata->k == 1.0 || labs (gwdata->c) == 1)) {
		uint32_t low_addin;
		int	i;

/* Zero the FFT data */

		memset (g, 0, gwnum_datasize (gwdata));

/* To make the mod k*b^n+c step faster, gwnum's are pre-multiplied by 1/k. */
/* Case 1 (k*b^n-1): Inverse of k is b^n.  Case 2 (k*b^n+1): Inverse of k is -b^n. */
/* Since the value we are adding (or subtracting) from the top FFT words can be */
/* greater than k, we must divide value by k and wrap-around to add it in to the */
/* lowest FFT words.  The remainder gets added to (or subtracted from) the upper FFT words. */

		if (gwdata->k > 1.0) {
			double	quotient, remainder, high_addin;
			unsigned long word, bit_in_word;

			quotient = floor ((double) a->n[0] / gwdata->k);
			remainder = (double) a->n[0] - quotient * gwdata->k;
			if (gwdata->c == -1) {
				low_addin = (uint32_t) quotient;
				high_addin = remainder;
			} else {
				low_addin = (uint32_t) quotient + 1;
				high_addin = gwdata->k - remainder;
			}

/* Add the high_addin value across the top FFT words */

			bitaddr (gwdata, gwdata->n, &word, &bit_in_word);
			for (i = word; high_addin > 0.0; i++, bit_in_word = 0) {
				long	value, maxval;

				maxval = intpow (gwdata->b, gwdata->NUM_B_PER_SMALL_WORD);
				if (is_big_word (gwdata, i)) maxval *= gwdata->b;

				if (bit_in_word == 0) {
					quotient = floor (high_addin / (double) maxval);
					remainder = high_addin - quotient * (double) maxval;
					value = (long) remainder;
					high_addin = quotient;
				} else {
					long	pow_fudge, fudged_maxval;
					pow_fudge = intpow (gwdata->b, bit_in_word);
					fudged_maxval = maxval / pow_fudge;
					quotient = floor (high_addin / (double) fudged_maxval);
					remainder = high_addin - quotient * (double) fudged_maxval;
					value = (long) remainder * pow_fudge;
					high_addin = quotient;
				}

				if (value > (maxval >> 1)) {
					value = value - maxval;
					high_addin += 1.0;
				}
				if (i == gwdata->FFTLEN - 1) {
					value = value + (long) high_addin * maxval;
					high_addin = 0.0;
				}
				set_fft_value (gwdata, g, i, value);
			}
		}

/* If k is 1, simply copy the giant value to the low FFT words */

		else
			low_addin = a->n[0];

/* Spread the low_addin value across the lowest FFT words as necessary */

		for (i = 0; low_addin; i++) {
			long	value, maxval;

			maxval = intpow (gwdata->b, gwdata->NUM_B_PER_SMALL_WORD);
			if (is_big_word (gwdata, i)) maxval *= gwdata->b;

			value = low_addin % maxval;
			low_addin = low_addin / maxval;
			if (value > (maxval >> 1)) {
				value = value - maxval;
				low_addin++;
			}
			set_fft_value (gwdata, g, i, value);
		}
	}

/* To make the mod k*b^n+c step faster, gwnum's are pre-multiplied by 1/k */
/* If k is greater than 1, then we calculate the inverse of k, multiply */
/* the giant by the inverse of k, and do a mod k*b^n+c. */

	else {
		unsigned long i, limit, carry;
		giant	newg = NULL;

		if (gwdata->k > 1.0) {
			newg = popg (&gwdata->gdata, (((unsigned long) gwdata->bit_length >> 5) + 1) * 2);

			/* Easy case 1 (k*b^n-1): Inverse of k is b^n */

			if (gwdata->c == -1) {
				if (gwdata->b == 2) {
					gtog (a, newg);
					gshiftleft (gwdata->n, newg);
				} else {
					itog (gwdata->b, newg);
					power (newg, gwdata->n);
					mulgi (&gwdata->gdata, a, newg);
				}
			}

			/* Easy case 2 (k*b^n+1): Inverse of k is -b^n */

			else if (gwdata->c == 1) {
				if (gwdata->b == 2) {
					gtog (a, newg);
					gshiftleft (gwdata->n, newg);
					negg (newg);
				} else {
					itog (gwdata->b, newg);
					power (newg, gwdata->n);
					negg (newg);
					mulgi (&gwdata->gdata, a, newg);
				}
			}

			else {				/* General inverse case */
				giant	n;
				n = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 1);
				ultog (gwdata->b, n);		/* Compute k*b^n+c */
				power (n, gwdata->n);
				dblmulg (gwdata->k, n);
				iaddg (gwdata->c, n);
				dbltog (gwdata->k, newg);	/* Compute 1/k */
				invg (n, newg);
				ASSERTG (newg->sign > 0);	/* Assert inverse found */
				mulgi (&gwdata->gdata, a, newg); /* Multiply input num by 1/k */
				pushg (&gwdata->gdata, 1);
			}

			specialmodg (gwdata, newg);
			a = newg;
		}

/* Figure out how many FFT words we will need to set */

		limit = (unsigned long) ceil ((double) bitlen (a) / (gwdata->avg_num_b_per_word * log2 (gwdata->b)));
		if (limit > gwdata->FFTLEN) limit = gwdata->FFTLEN;
		if (gwdata->ZERO_PADDED_FFT && limit > gwdata->FFTLEN / 2 +4) limit = gwdata->FFTLEN / 2 + 4;

/* Now convert the giant to FFT format.  For base 2 we simply copy bits.  */

		if (gwdata->b == 2) {
			unsigned long mask1, mask2, e1len;
			int	bits1, bits2, bits_in_next_binval;
			unsigned long binval;
			uint32_t *e1;

			e1len = a->sign;
			e1 = a->n;

			bits1 = gwdata->NUM_B_PER_SMALL_WORD;
			bits2 = bits1 + 1;
			mask1 = (1L << bits1) - 1;
			mask2 = (1L << bits2) - 1;
			if (e1len) {binval = *e1++; e1len--; bits_in_next_binval = 32;}
			else binval = 0;
			carry = 0;
			for (i = 0; i < limit; i++) {
				int	big_word, bits;
				long	value, mask;
				big_word = is_big_word (gwdata, i);
				bits = big_word ? bits2 : bits1;
				mask = big_word ? mask2 : mask1;
				if (i == limit - 1) value = binval;
				else value = binval & mask;
				value = value + carry;
				if (value > (mask >> 1) && bits > 1 && i != limit - 1) {
					value = value - (mask + 1);
					carry = 1;
				} else {
					carry = 0;
				}
				set_fft_value (gwdata, g, i, value);

				binval >>= bits;
				if (e1len == 0) continue;
				if (bits_in_next_binval < bits) {
					if (bits_in_next_binval)
						binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
					bits -= bits_in_next_binval;
					e1++; e1len--; bits_in_next_binval = 32;
					if (e1len == 0) continue;
				}
				if (bits) {
					binval |= (*e1 >> (32 - bits_in_next_binval)) << (32 - bits);
					bits_in_next_binval -= bits;
				}
			}
			for ( ; i < gwdata->FFTLEN; i++) {
				set_fft_value (gwdata, g, i, 0);
			}
		}

/* Otherwise (non-base 2), we do a recursive divide and conquer radix conversion. */
/* The recursive routine writes on a, so make a copy before calling */

		else {
			if (a != newg) {
				newg = popg (&gwdata->gdata, a->sign * 2);
				gtog (a, newg);
				a = newg;
			}
			carry = nonbase2_gianttogw (gwdata, a, g, limit, 0, 0);
		}

/* Write carry, if any, to FFT data */

		if (carry) set_fft_value (gwdata, g, limit++, carry);

/* Clear the upper words */

		for (i = limit; i < gwdata->FFTLEN; i++)
			set_fft_value (gwdata, g, i, 0);

/* Free allocated memory */

		if (a == newg) pushg (&gwdata->gdata, 1);
	}

/* Clear various flags */

	((uint32_t *) g)[-1] = 1; /* Set unnormalized add counter */
	((uint32_t *) g)[-7] = 0; /* Clear has been completely/partially FFTed flag */
}

/* Internal recursive routine to convert a giant to gwnum FFT format. */

long nonbase2_gianttogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	a,
	gwnum	g,
	unsigned long limit,	/* How many FFT words to set */
	unsigned long offset,	/* Offset into FFT array of words to set */
	long	carry)		/* Carry to add into this section */		 
{
	ASSERTG (a->sign >= 0);		/* We only handle positive numbers */

/* If we are converting a lot of words, divide and conquer. */

	if (limit >= 50) {
		giant	upper, tmp;
		int	num_b;
		unsigned long half_limit = limit >> 1;

		tmp = popg (&gwdata->gdata, limit * 2);
		upper = popg (&gwdata->gdata, a->sign * 2);
		num_b = gwfft_base (gwdata->dd_data, offset + half_limit) - gwfft_base (gwdata->dd_data, offset);
		itog (gwdata->b, tmp);
		power (tmp, num_b);
		gtog (a, upper);
		divg (tmp, upper);
		mulgi (&gwdata->gdata, upper, tmp);
		subg (tmp, a);
		carry = nonbase2_gianttogw (gwdata, a, g, half_limit, offset, carry);
		carry = nonbase2_gianttogw (gwdata, upper, g, limit - half_limit, offset + half_limit, carry);
		pushg (&gwdata->gdata, 2);
	}

/* Convert the giant to FFT format */

	else {
		giant	newg, tmp;
		unsigned long i, mask1, mask2;
		long	value;

		newg = popg (&gwdata->gdata, a->sign * 2);
		tmp = popg (&gwdata->gdata, a->sign * 2);

		mask1 = intpow (gwdata->b, gwdata->NUM_B_PER_SMALL_WORD);
		mask2 = gwdata->b * mask1;
		for (i = offset; i < offset + limit; i++) {
			unsigned long mask;

			mask = is_big_word (gwdata, i) ? mask2 : mask1;

			gtog (a, newg);
			if (i != gwdata->FFTLEN - 1) {
				ultog (mask, tmp);
				divg (tmp, a);
				mulgi (&gwdata->gdata, a, tmp);
				subg (tmp, newg);
			}
			value = (newg->sign) ? newg->n[0] : 0;
			value += carry;

			if (value > (long) (mask >> 1) && i != gwdata->FFTLEN - 1) {
				value = value - mask;
				carry = 1;
			} else {
				carry = 0;
			}
			set_fft_value (gwdata, g, i, value);
		}
		pushg (&gwdata->gdata, 2);
	}

/* Return carry for next section */

	return (carry);
}


/* Convert a gwnum to a binary value.  Returns the number of 32-bit values */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinary (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	uint32_t *array,	/* Array to contain the binary value */
	uint32_t arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	if (((uint32_t *) n)[-7] == 3) return (GWERROR_FFT); /* Test the FFTed flag */
	if (((uint32_t *) n)[-7] == 1) return (GWERROR_PARTIAL_FFT); /* Test the FFT-started flag */
	tmp = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
	err_code = gwtogiant (gwdata, n, tmp);
	if (err_code < 0) return (err_code);
	ASSERTG ((unsigned long) tmp->sign <= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (uint32_t));
	pushg (&gwdata->gdata, 1);
	return (tmp->sign);
}

/* Convert a gwnum to a binary value.  Returns the number of 64-bit values */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinary64 (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	uint64_t *array,	/* Array to contain the binary value */
	uint64_t arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	if (((uint32_t *) n)[-7] == 3) return (GWERROR_FFT); /* Test the FFTed flag */
	if (((uint32_t *) n)[-7] == 1) return (GWERROR_PARTIAL_FFT); /* Test the FFT-started flag */
	tmp = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
	err_code = gwtogiant (gwdata, n, tmp);
	if (err_code < 0) return (err_code);
	tmp->n[tmp->sign] = 0;
	tmp->sign = (tmp->sign + 1) / 2;
	ASSERTG ((unsigned long) tmp->sign <= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (uint64_t));
	pushg (&gwdata->gdata, 1);
	return (tmp->sign);
}

/* Convert a gwnum to a binary value.  Returns the number of longs */
/* written to the array.  The array is NOT zero-padded.  Returns a */
/* negative number if an error occurs during the conversion.  An error */
/* can happen if the FFT data contains a NaN or infinity value. */

long gwtobinarylongs (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	unsigned long *array,	/* Array to contain the binary value */
	unsigned long arraylen)	/* Maximum size of the array */
{
	giant	tmp;
	int	err_code;

	if (((uint32_t *) n)[-7] == 3) return (GWERROR_FFT); /* Test the FFTed flag */
	if (((uint32_t *) n)[-7] == 1) return (GWERROR_PARTIAL_FFT); /* Test the FFT-started flag */
	tmp = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
	err_code = gwtogiant (gwdata, n, tmp);
	if (err_code < 0) return (err_code);
	if (sizeof (unsigned long) > sizeof (uint32_t)) {
		tmp->n[tmp->sign] = 0;
		tmp->sign = (tmp->sign + 1) / 2;
	}
	ASSERTG ((unsigned long) tmp->sign <= arraylen);
	memcpy (array, tmp->n, tmp->sign * sizeof (unsigned long));
	pushg (&gwdata->gdata, 1);
	return (tmp->sign);
}

/* Convert a gwnum to a giant.  WARNING: Caller must allocate an array that */
/* is several words larger than the maximum result that can be returned. */
/* This is a gross kludge that lets gwtogiant use the giant for intermediate */
/* calculations. */

int gwtogiant (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	gg,
	giant	v)
{
	int	err_code;
	unsigned long limit;

/* Set result to zero in case of error.  If caller does not check the returned error code */
/* a result value of zero is less likely to cause problems/crashes. */

	v->sign = 0;

/* Sanity check the input number */

	if (((uint32_t *) gg)[-7] == 3) return (GWERROR_FFT); /* Test the FFTed flag */
	if (((uint32_t *) gg)[-7] == 1) return (GWERROR_PARTIAL_FFT); /* Test the FFT-started flag */

/* If this is a general-purpose mod, then only convert the needed words */
/* which will be less than half the FFT length.  If this is a zero padded */
/* FFT, then only convert a little more than half of the FFT data words. */
/* For a DWT, convert all the FFT data. */

	if (gwdata->GENERAL_MOD) limit = gwdata->GW_GEN_MOD_MAX + 3;
	else if (gwdata->ZERO_PADDED_FFT) limit = gwdata->FFTLEN / 2 + 4;
	else limit = gwdata->FFTLEN;

/* GENERAL_MOD has some strange cases we must handle.  In particular the */
/* last fft word translated can be 2^bits and the next word could be -1, */
/* this must be translated into zero, zero. */

	if (gwdata->GENERAL_MOD) {
		long	val, prev_val;
		while (limit < gwdata->FFTLEN) {
			err_code = get_fft_value (gwdata, gg, limit, &val);
			if (err_code) return (err_code);
			if (val == -1 || val == 0) break;
			limit++;
			ASSERTG (limit <= gwdata->FFTLEN / 2 + 2);
			if (limit > gwdata->FFTLEN / 2 + 2) return (GWERROR_INTERNAL + 9);
		}
		while (limit > 1) {		/* Find top word */
			err_code = get_fft_value (gwdata, gg, limit-1, &prev_val);
			if (err_code) return (err_code);
			if (val != prev_val || val < -1 || val > 0) break;
			limit--;
		}
		limit++;
	}

/* If base is 2 we can simply copy the bits out of each FFT word */

	if (gwdata->b == 2) {
		long	val;
		int	j, bits, bitsout, carry;
		unsigned long i;
		uint32_t *outptr;

/* Collect bits until we have all of them */

		carry = 0;
		bitsout = 0;
		outptr = v->n;
		*outptr = 0;
		for (i = 0; i < limit; i++) {
			err_code = get_fft_value (gwdata, gg, i, &val);
			if (err_code) return (err_code);
			bits = gwdata->NUM_B_PER_SMALL_WORD;
			if (is_big_word (gwdata, i)) bits++;
			val += carry;

			carry = (val >> bits);
			val -= (carry << bits);
			*outptr += (val << bitsout);
			bitsout += bits;
			if (bitsout >= 32) {
				bitsout -= 32;
				*++outptr = (val >> (bits - bitsout));
			}
		}

/* Finish outputting the last word and any carry data */

		if (bitsout) {
			*outptr++ += (carry << bitsout);
			carry >>= (32 - bitsout);
		}
		if (carry != -1 && carry != 0) {
			*outptr++ = carry;
			carry >>= 31;
		}

/* Set the length */

		v->sign = (long) (outptr - v->n);
		while (v->sign && v->n[v->sign-1] == 0) v->sign--;

/* If carry is -1, the gwnum is negative.  Ugh.  Flip the bits and sign. */

		if (carry == -1) {
			for (j = 0; j < v->sign; j++) v->n[j] = ~v->n[j];
			while (v->sign && v->n[v->sign-1] == 0) v->sign--;
			iaddg (1, v);
			v->sign = -v->sign;
		}
	}

/* Otherwise (base is not 2) we must do a radix conversion */

	else {
		giantstruct *array = NULL;
		uint32_t *buf = NULL;
		giant	small_base = NULL;
		giant	large_base = NULL;
		unsigned long i, gap, small_size, last_small_size;

		array = (giantstruct *) malloc (limit * sizeof (giantstruct));
		if (array == NULL) { memerr: err_code = GWERROR_MALLOC; goto err; }
		buf = (uint32_t *) malloc (limit * sizeof (uint32_t));
		if (buf == NULL) goto memerr;
		small_base = popg (&gwdata->gdata, limit);
		if (small_base == NULL) goto memerr;
		large_base = popg (&gwdata->gdata, limit);
		if (large_base == NULL) goto memerr;
		for (i = 0; i < limit; i++) {
			long	val;
			err_code = get_fft_value (gwdata, gg, i, &val);
			if (err_code) {
err:				free (array);
				free (buf);
				if (small_base != NULL) pushg (&gwdata->gdata, 1);
				if (large_base != NULL) pushg (&gwdata->gdata, 1);
				return (err_code);
			}
			array[i].n = &buf[i];
			setmaxsize(&array[i], limit);
			itog (val, &array[i]);
		}

/* Loop combining pairs into ever larger and larger numbers.  Do all but last combining pass. */

		gap = 1;
		while (gap + gap < limit) {
			small_size = gwfft_base (gwdata->dd_data, gap) - 1;
			if (gap == 1)
				itog (intpow (gwdata->b, small_size), small_base);
			else if (small_size == last_small_size * 2)
				squaregi (&gwdata->gdata, small_base);
			else
				mulgi (&gwdata->gdata, large_base, small_base);
			itog (gwdata->b, large_base);
			mulgi (&gwdata->gdata, small_base, large_base);
			for (i = 0; i + gap < limit; i += gap + gap) {
				gtog (&array[i+gap], v);
				if (gwfft_base (gwdata->dd_data, i+gap) - gwfft_base (gwdata->dd_data, i) == small_size)
					mulgi (&gwdata->gdata, small_base, v);
				else
					mulgi (&gwdata->gdata, large_base, v);
				addg (v, &array[i]);
			}
			gap = gap << 1;
			last_small_size = small_size;
		}

/* Do the last combining pass, outputting result directly to v. */

		if (gwfft_base (gwdata->dd_data, gap) == small_size * 2 + 1)
			mulgi (&gwdata->gdata, small_base, large_base);
		else
			squaregi (&gwdata->gdata, large_base);
		gtog (&array[gap], v);
		mulgi (&gwdata->gdata, large_base, v);
		addg (&array[0], v);

/* Clean up */

		free (array);
		free (buf);
		pushg (&gwdata->gdata, 2);
	}

/* Since all gwnums are premultiplied by the inverse of k, we must now */
/* multiply by k to get the true result. */

	if (gwdata->k > 1.0) {
		stackgiant(k,2);
		dbltog (gwdata->k, k);
		mulgi (&gwdata->gdata, k, v);
	}

/* The gwnum is not guaranteed to be smaller than k*b^n+c.  Handle this */
/* possibility.  This also converts negative values to positive. */

	specialmodg (gwdata, v);

/* Return success */

	return (0);
}

/* Special modg.  This is a fast implementation of mod k*2^n+c using just */
/* shifts, adds, and divide and mul by small numbers.  All others moduli */
/* call the slow giants code. */

void specialmodg (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	giant	g)
{
	int	neg, count;
	giant	n;

/* If the modulus is a general-purpose number, then let the giants code */
/* do the work.  This is done for both GENERAL_MOD and (k*b^n+c)/d cases. */

	if (gwdata->GW_MODULUS != NULL) {
		modgi (&gwdata->gdata, gwdata->GW_MODULUS, g);
		return;
	}

/* Calculate the modulo number - k*b^n+c */

	n = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 1);
	ultog (gwdata->b, n);
	power (n, gwdata->n);
	dblmulg (gwdata->k, n);
	iaddg (gwdata->c, n);

/* If b is not 2 let the giants code do the work. */

	if (gwdata->b != 2) {
		modgi (&gwdata->gdata, n, g);
		pushg (&gwdata->gdata, 1);
		return;
	}

/* Do the quick modulus code twice because in the case where */
/* abs(c) > k once won't get us close enough. */

	neg = FALSE;
	for (count = 0; count < 2; count++) {

/* Handle negative input values */

	    neg ^= (g->sign < 0);
	    g->sign = abs (g->sign);

/* If number is bigger than the modulus, do a mod using shifts and adds */
/* This will get us close to the right answer. */

	    if (gcompg (g, n) > 0) {
		giant	tmp1;

/* Allocate temporary */

		tmp1 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);

/* Calculate the modulo by dividing the upper bits of k, multiplying by */
/* c and subtracting that from the bottom bits. */

		gtogshiftright (gwdata->n, g, tmp1);	// Upper bits
		gmaskbits (gwdata->n, g);		// Lower bits

		if (gwdata->k == 1.0) {
			imulg (gwdata->c, tmp1);	// Upper bits times C
			subg (tmp1, g);
		} else {
			giant	tmp2, tmp3;

			tmp2 = popg (&gwdata->gdata, (((unsigned long) gwdata->bit_length >> 5) + 5) * 2);
			tmp3 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);

			gtog (tmp1, tmp2);
			dbltog (gwdata->k, tmp3);
			divgi (&gwdata->gdata, tmp3, tmp1);	// Upper bits over K
			mulgi (&gwdata->gdata, tmp1, tmp3);
			subg (tmp3, tmp2);	// Upper bits mod K

			gshiftleft (gwdata->n, tmp2);
			addg (tmp2, g);		// Upper bits mod K+lower bits

			imulg (gwdata->c, tmp1);// Upper bits over K times C
			subg (tmp1, g);
			pushg (&gwdata->gdata, 2);
		}

		pushg (&gwdata->gdata, 1);
	    }
	}

/* Add or subtract n until the g is between 0 and n-1 */

	while (g->sign < 0) addg (n, g);
	while (gcompg (g, n) >= 0) subg (n, g);

/* If input was negative, return k*b^n+c - g */

	if (neg && g->sign) {
		g->sign = -g->sign;
		addg (n, g);
	}

/* Free memory */

	pushg (&gwdata->gdata, 1);
}

/* Test if a gwnum is zero.  This routine was originally written by Jean Penne. */
/* It has not been adequately tested and MAY NOT BE BUG-FREE.  Use at your own risk! */

int gwiszero (
	gwhandle *gwdata,
	gwnum 	gg)
{
	unsigned long j;
	int 	result, count;

	ASSERTG (((uint32_t *) gg)[-1] >= 1);
	ASSERTG (((uint32_t *) gg)[-7] == 0);

/* If the input number is the result of an unormalized addition or subtraction, then */
/* we had better normalize the number! */

	if (((uint32_t *) gg)[-1] > 1) {
		struct gwasm_data *asm_data;
		gwnum	gwnorm;

		gwnorm = gwalloc (gwdata);
		if (gwnorm == NULL) return (-GWERROR_MALLOC);
		dbltogw (gwdata, 0.0, gwnorm);
		asm_data = (struct gwasm_data *) gwdata->asm_data;
		asm_data->SRCARG = gg;
		asm_data->SRC2ARG = gwnorm;
		asm_data->DESTARG = gwnorm;
		gw_add (gwdata, asm_data);
		result = gwiszero (gwdata, gwnorm);
		gwfree (gwdata, gwnorm);
		return (result);
	}

/* CONCERN!!!  Could the result of a normalized multiply be greater than k*b^n+c? */
/* If so, we should test the top FFT word and if it is bigger than the maximum valid */
/* value, do a normalizing add identical to the code above. */	

/* Look through all the FFT data.  If each FFT word is zero, then the gwnum is zero. */
/* If we run into just a few non-zero FFT elements, then the gwnum might still be zero */
/* because of the way carries are propagated in the assembly code.  If we run into */
/* a large number of non-zero FFT words then the gwnum is not zero. */

#define	MAX_NZ_COUNT 16
	count = 0;
	for (j = 0; j < gwdata->FFTLEN; j++) {
		double	val;
		val = * addr (gwdata, gg, j);
		if (! is_valid_double (val)) return (GWERROR_BAD_FFT_DATA);
		if (val == 0.0) continue;
		if (++count > MAX_NZ_COUNT) return (FALSE);	// Too many non zero words, the gwnum is not zero.
	}
	if (count) {			// The gwnum may be zero but needs a more accurate test...
		giant	gtest;
		gtest = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 5) + 5);
		if (gtest == NULL) return (-GWERROR_MALLOC);
		gwtogiant (gwdata, gg, gtest);
		result = isZero (gtest);
		pushg (&gwdata->gdata, 1);
		return (result);
	}
	else
		return (TRUE);			// The gwnum is zero
}

/* Test two gwnums for equality.  Written by Jean Penne.  Uses the gwiszero routine */
/* which MAY NOT BE BUG-FREE.  Use this routine at your own risk. */

int gwequal (
	gwhandle *gwdata,
	gwnum	gw1,
	gwnum	gw2)
{
	struct gwasm_data *asm_data;
	gwnum	gwdiff;
	int	result;

	ASSERTG (((uint32_t *) gw1)[-1] >= 1);
	ASSERTG (((uint32_t *) gw2)[-1] >= 1);
	ASSERTG (((uint32_t *) gw1)[-7] == 0);
	ASSERTG (((uint32_t *) gw2)[-7] == 0);

/* Allocate memory for the difference */

	gwdiff = gwalloc (gwdata);

/* Do a normalized subtract */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->SRCARG = gw1;
	asm_data->SRC2ARG = gw2;
	asm_data->DESTARG = gwdiff;
	gw_sub (gwdata, asm_data);
	((uint32_t *) gwdiff)[-1] = 1;
	((uint32_t *) gwdiff)[-7] = 0;

/* The two input numbers are equal if the difference is zero */

	result = gwiszero (gwdata, gwdiff);

/* Cleanup and return result */

	gwfree (gwdata, gwdiff);
	return (result);
}

/******************************************************************/
/* Wrapper routines for the multiplication assembly code routines */
/******************************************************************/

/* Internal wrapper routine to call fftmul assembly code. */
/* Caller must set NORMRTN prior to calling this routine.*/

void raw_gwfftmul (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,
	gwnum	d)
{
	struct gwasm_data *asm_data;
	uint32_t norm_count1, norm_count2;
	double	sumdiff;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) d)[-1] >= 1);

/* Get the unnormalized add count for later use */

	norm_count1 = ((uint32_t *) s)[-1];
	norm_count2 = ((uint32_t *) d)[-1];

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->DESTARG = d;
	asm_data->DIST_TO_FFTSRCARG = 0;
	asm_data->DIST_TO_MULSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->ffttype = 3;
	asm_data->zero_fft = (asm_data->NORMRTN == gwdata->GWPROCPTRS[zerohigh_routines + (gwdata->NORMNUM & 1)]);
	asm_data->const_fft = gwdata->NORMNUM & 2;
	asm_data->add_sub_smallmul_op = 0;
//if (rand() % 100 < 1) *s += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *s += 1.0E200 * 1.0E200;
	gw_fft (gwdata, asm_data);
//if (rand() % 100 < 1) *d += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *d += 1.0E200 * 1.0E200;
	if (! is_valid_double (gwsumout (gwdata, d))) gwdata->GWERROR |= 1;
	gwdata->fft_count += 2.0;
	((uint32_t *) d)[-7] = gwdata->POSTFFT; // Set has-been-partially-FFTed flag

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to */
/* compensate. */

	if (norm_count1 != 1 || norm_count2 != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double)norm_count1 * (double)norm_count2);
		gwsuminp (gwdata, d) *= adjustment;
		gwsumout (gwdata, d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (gwdata, d) - gwsumout (gwdata, d);
	if (fabs (sumdiff) > gwdata->MAXDIFF) gwdata->GWERROR |= 2; 

/* Reset the unnormalized add count */

	((uint32_t *) d)[-1] = 1;
}

/* Common code to emulate the modulo with two multiplies in the general purpose case */

void emulate_mod (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s)		/* Source and destination */
{
	struct gwasm_data *asm_data;
	gwnum	tmp;
	double	saved_addin_value;

	ASSERTG (* addr (gwdata, s, gwdata->FFTLEN-1) > -2.0 && * addr (gwdata, s, gwdata->FFTLEN-1) <= 0.0);
	if (* addr (gwdata, s, gwdata->FFTLEN-1) <= -2.0 || * addr (gwdata, s, gwdata->FFTLEN-1) > 0.0)
		gwdata->GWERROR = GWERROR_INTERNAL + 10;

/* Save and clear the addin value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	saved_addin_value = asm_data->ADDIN_VALUE;
	asm_data->ADDIN_VALUE = 0.0;

/* Copy the number and zero out the low words. */

	tmp = gwalloc (gwdata);
	gwcopyzero (gwdata, s, tmp, gwdata->GW_ZEROWORDSLOW);

/* Multiply by the reciprocal that has been carefully shifted so that the */
/* integer part of the result wraps to the lower FFT words.  Adjust the */
/* normalization routine so that the FFT code zeroes the high FFT words */
/* and we are left with just the quotient! */

	asm_data->NORMRTN = gwdata->GWPROCPTRS[zerohigh_routines + (gwdata->NORMNUM & 1)];
	raw_gwfftmul (gwdata, gwdata->GW_RECIP_FFT, tmp);
	ASSERTG (* addr (gwdata, tmp, gwdata->FFTLEN/2-1) > -2.0 && * addr (gwdata, tmp, gwdata->FFTLEN/2-1) <= 0.0);
	if (* addr (gwdata, tmp, gwdata->FFTLEN/2-1) <= -2.0 || * addr (gwdata, tmp, gwdata->FFTLEN/2-1) > 0.0)
		gwdata->GWERROR = GWERROR_INTERNAL + 11;

/* Muliply quotient and modulus.  Select normalization routine that does */
/* not zero the high FFT words. */

	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + (gwdata->NORMNUM & 1)];
	raw_gwfftmul (gwdata, gwdata->GW_MODULUS_FFT, tmp);

/* Subtract from the original number to get the remainder */

	gwsub (gwdata, tmp, s);
	ASSERTG (* addr (gwdata, s, gwdata->FFTLEN-1) == 0.0);
	if (* addr (gwdata, s, gwdata->FFTLEN-1) != 0.0) gwdata->GWERROR = GWERROR_INTERNAL + 12;
	gwfree (gwdata, tmp);

/* Restore the addin value */

	asm_data->ADDIN_VALUE = saved_addin_value;
}

/* User-visible routines */

void gwfft (			/* Forward FFT */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number */
	gwnum	d)		/* Destination (can overlap source) */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) s)[-7] != 3);	// Make sure input has not already been completely FFTed (we could turn this into a nop or copy)

/* Copy the unnormalized add count */

	((uint32_t *) d)[-1] = ((uint32_t *) s)[-1];

/* If this is a zero-padded FFT and the source has been partially FFTed */
/* (by a prior POSTFFT setting) and the destination is different than the */
/* source, then we must also copy the 7 words around the halfway point from */
/* the source to the destination. */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	if (asm_data->ZERO_PADDED_FFT &&
	    ((uint32_t *) s)[-7] == 1 && s != d) {
		d[-5] = s[-5];
		d[-6] = s[-6];
		d[-7] = s[-7];
		d[-8] = s[-8];
		d[-9] = s[-9];
		d[-10] = s[-10];
		d[-11] = s[-11];
	}

/* Call the assembly code */

	asm_data->DESTARG = d;
	asm_data->DIST_TO_FFTSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->DIST_TO_MULSRCARG = 0;
	asm_data->ffttype = 1;
	gw_fft (gwdata, asm_data);
	gwdata->fft_count += 1.0;
	((uint32_t *) d)[-7] = 3; // Set has-been-FFTed flag
}

void gwsquare2 (		/* Square a number */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;
	uint32_t norm_count;
	double	sumdiff;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) s)[-7] != 3);	// Make sure input has not already been completely FFTed (we could turn this into a gwfftfftmul)

/* If we are converting gwsquare calls into gwsquare_carefully calls */
/* do so now.  Turn off option to do a partial forward FFT on the result. */
/* NOTE: We must clear count since gwsquare_carefully calls back to this */
/* gwsquare routine (well, it used to -- we'll still clear the count to be safe) */

	if (gwdata->square_carefully_count) {
		int	n = gwdata->square_carefully_count;
		if (n > 1) gwstartnextfft (gwdata, 0);
		gwdata->square_carefully_count = 0;
		gwsquare2_carefully (gwdata, s, d);
		gwdata->square_carefully_count = n - 1;
		return;
	}

/* Get the unnormalized add count for later use */

	norm_count = ((uint32_t *) s)[-1];

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + gwdata->NORMNUM];
	asm_data->DESTARG = d;
	asm_data->DIST_TO_FFTSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->DIST_TO_MULSRCARG = 0;
	asm_data->ffttype = 2;
	asm_data->zero_fft = 0;
	asm_data->const_fft = gwdata->NORMNUM & 2;
	asm_data->add_sub_smallmul_op = 0;
//if (rand() % 100 < 1) *s += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *s += 1.0E200 * 1.0E200;
	gw_fft (gwdata, asm_data);
//if (rand() % 100 < 1) *d += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *d += 1.0E200 * 1.0E200;
	if (! is_valid_double (gwsumout (gwdata, d))) gwdata->GWERROR |= 1;
	gwdata->fft_count += 2.0;
	((uint32_t *) d)[-7] = gwdata->POSTFFT; // Set has-been-partially-FFTed flag

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to compensate. */

	if (norm_count != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double) norm_count * (double) norm_count);
		gwsuminp (gwdata, d) *= adjustment;
		gwsumout (gwdata, d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (gwdata, d) - gwsumout (gwdata, d);
	if (fabs (sumdiff) > gwdata->MAXDIFF) gwdata->GWERROR |= 2; 

/* Reset the unnormalized add count */

	((uint32_t *) d)[-1] = 1;

/* Emulate mod with 2 multiplies case */

	if (gwdata->GENERAL_MOD) emulate_mod (gwdata, d);
}

void gwfftmul (			/* Multiply already FFTed source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Already FFTed source number */
	gwnum	d)		/* Non-FFTed source. Also destination */
{
	struct gwasm_data *asm_data;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) d)[-1] >= 1);
	ASSERTG (((uint32_t *) s)[-7] == 3);	// Make sure input has been completely FFTed
	ASSERTG (((uint32_t *) d)[-7] != 3);	// Make sure input has not already been completely FFTed (we could turn this into a gwfftfftmul)

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + gwdata->NORMNUM];
	raw_gwfftmul (gwdata, s, d);

/* Emulate mod with 2 multiplies case */

	if (gwdata->GENERAL_MOD) emulate_mod (gwdata, d);
}

void gwfftfftmul (		/* Multiply two already FFTed sources */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Already FFTed source number */
	gwnum	s2,		/* Already FFTed source number */
	gwnum	d)		/* Destination (can overlap sources) */
{
	struct gwasm_data *asm_data;
	uint32_t norm_count1, norm_count2;
	double	sumdiff;

	ASSERTG (((uint32_t *) s)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s)[-7] == 3);	// Make sure input has been completely FFTed
	ASSERTG (((uint32_t *) s2)[-7] == 3);	// Make sure input has been completely FFTed

/* Get the unnormalized add count for later use */

	norm_count1 = ((uint32_t *) s)[-1];
	norm_count2 = ((uint32_t *) s2)[-1];

/* Call the assembly code */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + gwdata->NORMNUM];
	asm_data->DESTARG = d;
	asm_data->DIST_TO_MULSRCARG = (intptr_t) s - (intptr_t) d;
	asm_data->DIST_TO_FFTSRCARG = (intptr_t) s2 - (intptr_t) d;
	asm_data->ffttype = 4;
	asm_data->zero_fft = 0;
	asm_data->const_fft = gwdata->NORMNUM & 2;
	asm_data->add_sub_smallmul_op = 0;
//if (rand() % 100 < 1) *s += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *s += 1.0E200 * 1.0E200;
//if (rand() % 100 < 1) *s2 += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *s2 += 1.0E200 * 1.0E200;
	gw_fft (gwdata, asm_data);
//if (rand() % 100 < 1) *d += 1.0;			// Generate random errors (for caller to test error recovery)
//if (rand() % 1000 < 2) *d += 1.0E200 * 1.0E200;
	if (! is_valid_double (gwsumout (gwdata, d))) gwdata->GWERROR |= 1;
	gwdata->fft_count += 1.0;
	((uint32_t *) d)[-7] = gwdata->POSTFFT; // Set has-been-partially-FFTed flag

/* Adjust if necessary the SUM(INPUTS) vs. SUM(OUTPUTS).  If norm_count */
/* is more than one, then the sums will be larger than normal.  This */
/* could trigger a spurious MAXDIFF warning.  Shrink the two SUMS to compensate. */

	if (norm_count1 != 1 || norm_count2 != 1) {
		double	adjustment;
		adjustment = 1.0 / ((double)norm_count1 * (double)norm_count2);
		gwsuminp (gwdata, d) *= adjustment;
		gwsumout (gwdata, d) *= adjustment;
	}

/* Test SUM(INPUTS) vs. SUM(OUTPUTS) */

	sumdiff = gwsuminp (gwdata, d) - gwsumout (gwdata, d);
	if (fabs (sumdiff) > gwdata->MAXDIFF) gwdata->GWERROR |= 2; 

/* Reset the unnormalized add count */

	((uint32_t *) d)[-1] = 1;

/* Emulate mod with 2 multiplies case */

	if (gwdata->GENERAL_MOD) emulate_mod (gwdata, d);
}

void gwmul (			/* Multiply source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number (changed to FFTed source!) */
	gwnum	d)		/* Source and destination */
{
	gwfft (gwdata, s, s);
	gwfftmul (gwdata, s, d);
}

void gwsafemul (		/* Multiply source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number (not changed) */
	gwnum	d)		/* Source and destination */
{
	gwnum	qqq;

	qqq = gwalloc (gwdata);
	gwfft (gwdata, s, qqq);
	gwfftmul (gwdata, qqq, d);
	gwfree (gwdata, qqq);
}

/* Generate random FFT data.  We used to use the C runtime library. */
/* However, when a caller discovered a bug in gwsquare_carefully it */
/* very difficult to track down because the bug was not reproducible. */
/* We could make bugs reproducible by calling srand with a fixed value, */
/* but it is bad form for a library to do this.  Thus, we found a */
/* public domain random number generator to use. */

void gw_random_number (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	x)
{
	struct mt_state rand_info;
	giant	g;
	unsigned long i, len;

/* Init the random generator to a reproducible state for gwsquare_carefully. */
/* For other users of this routine, initialize to a more random state. */

	if (x == gwdata->GW_RANDOM)
		init_genrand (&rand_info, 5489);
	else
		init_genrand (&rand_info, (unsigned long) time (NULL));

/* Generate the random number */

	len = (((unsigned long) gwdata->bit_length) >> 5) + 1;
	g = popg (&gwdata->gdata, len);
	for (i = 0; i < len; i++) {
		g->n[i] = genrand_int32 (&rand_info);
	}
	g->sign = len;
	specialmodg (gwdata, g);
	gianttogw (gwdata, g, x);
	pushg (&gwdata->gdata, 1);
}

/* The FFT selection code assumes FFT data will essentially be random data */
/* yielding pretty well understood maximum round off errors.  When working */
/* with some numbers, especially at the start of a PRP exponentiation, the */
/* FFT data is decidedly not random, leading to much larger than expected */
/* roundoff errors.  In my own PRP code, I call gwsquare_carefully for the */
/* first 30 iterations.  To make this easier (and code more readable) you */
/* can call this routine and the next n gwsquare calls will be replaced by */
/* gwsquare_carefully calls.  If you pass an n of -1, the gwnum code will */
/* use a default value for n that should be suitable for getting a PRP */
/* exponentiation into a "random data state".  This routine can be called */
/* before gwsetup is called. */

void gwset_square_carefully_count (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	int	n)		/* Number of gwsquare calls to do carefully. */
				/* If n is -1, a default value is used */
{

/* Set the default count to enough iterations for a PRP exponentiation to */
/* generate a number larger than the modulus.  Then do an extra dozen */
/* iterations to hopefully scramble the data into a nice random looking */
/* pattern. */

	if (n == -1 && gwdata->FFTLEN) n = (int) ceil (log2 (gwdata->bit_length)) + 12;

/* Now remember the count */

	gwdata->square_carefully_count = n;
}

/* Square a number using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input number has not been partially or fully FFTed. */

void gwsquare2_carefully (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data;
	gwnum	tmp1;
	double	saved_addin_value;
	unsigned long saved_extra_bits;

/* Generate a random number, if we have't already done so */

	if (gwdata->GW_RANDOM == NULL) {
		gwdata->GW_RANDOM = gwalloc (gwdata);
		gw_random_number (gwdata, gwdata->GW_RANDOM);
	}

/* Save the addin value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	saved_addin_value = asm_data->ADDIN_VALUE;

/* Make sure we do not do addquick when computing s+random. */
/* If we do not do this, then the non-randomness of s can swamp */
/* the randomness of tmp1.  An example is the first PRP iterations */
/* of 2*3^599983-1 -- s is all positive values and gwadd3 thinks */
/* there are enough extra bits to do an add quick.  This generates */
/* temps with nearly all positive values -- very bad. */

	saved_extra_bits = gwdata->EXTRA_BITS;
	gwdata->EXTRA_BITS = 0;

/* Now do the squaring using two multiplies and several adds. */

	tmp1 = gwalloc (gwdata);
	gwstartnextfft (gwdata, 0);			/* Disable POSTFFT */
	gwaddsub4 (gwdata, s, gwdata->GW_RANDOM, tmp1, d); /* Compute s+random and s-random */
	gwmul (gwdata, tmp1, d);			/* Compute (s+random)(s-random) */
	asm_data->ADDIN_VALUE = 0.0;			/* Clear the addin value */
	gwfft (gwdata, gwdata->GW_RANDOM, tmp1);
	gwfftfftmul (gwdata, tmp1, tmp1, tmp1);		/* Compute random^2 (we could compute this once and save it) */
	gwadd3 (gwdata, d, tmp1, d);			/* Calc s^2 from 2 results */

/* Restore state, free memory and return */

	asm_data->ADDIN_VALUE = saved_addin_value;	/* Restore the addin value */
	gwdata->EXTRA_BITS = saved_extra_bits;
	gwfree (gwdata, tmp1);
}

/* Multiply numbers using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input numbers have not been partially or fully FFTed. */

void gwmul_carefully (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source -- preserved */
	gwnum	t)		/* Source and destination */
{
	struct gwasm_data *asm_data;
	gwnum	tmp1, tmp3, tmp4;
	double	saved_addin_value;
	unsigned long saved_extra_bits;

/* Generate a random number, if we have't already done so */

	if (gwdata->GW_RANDOM == NULL) {
		gwdata->GW_RANDOM = gwalloc (gwdata);
		gw_random_number (gwdata, gwdata->GW_RANDOM);
	}

/* Save and clear the addin value */

	asm_data = (struct gwasm_data *) gwdata->asm_data;
	saved_addin_value = asm_data->ADDIN_VALUE;

/* Make sure we do not do addquick when computing s+random and t+random. */

	saved_extra_bits = gwdata->EXTRA_BITS;
	gwdata->EXTRA_BITS = 0;

/* Now do the multiply using three multiplies and several adds */

	tmp1 = gwalloc (gwdata);
	tmp3 = gwalloc (gwdata);
	tmp4 = gwalloc (gwdata);

	gwadd3 (gwdata, s, gwdata->GW_RANDOM, tmp1);	/* Compute s+random */
	gwadd3 (gwdata, t, gwdata->GW_RANDOM, t);	/* Compute t+random */
	gwadd3 (gwdata, tmp1, gwdata->GW_RANDOM, tmp3);	/* Compute s+2*random */
	gwadd3 (gwdata, t, gwdata->GW_RANDOM, tmp4);	/* Compute t+2*random */

	gwstartnextfft (gwdata, 0);			/* Disable POSTFFT */
	gwmul (gwdata, tmp1, t);			/* Compute (s+r)*(t+r) = st + rs + rt + rr + addin */
	gwmul (gwdata, tmp3, tmp4);			/* Compute (s+2r)*(t+2r) = st + 2rs + 2rt + 4rr + addin */

	asm_data->ADDIN_VALUE = 0.0;
	gwfft (gwdata, gwdata->GW_RANDOM, tmp1);
	gwfftfftmul (gwdata, tmp1, tmp1, tmp1);		/* Compute random^2 (we could compute this once and save it) */

	gwaddquick (gwdata, t, t);			/* Compute 2st + 2rs + 2rt + 2rr + 2addin */
	gwsubquick (gwdata, tmp4, t);			/* Compute st - 2rr + addin */
	gwaddquick (gwdata, tmp1, t);			/* Compute st - rr + addin */
	gwadd (gwdata, tmp1, t);			/* Compute st + addin */

/* Restore state, free memory and return */

	gwdata->EXTRA_BITS = saved_extra_bits;
	asm_data->ADDIN_VALUE = saved_addin_value;
	gwfree (gwdata, tmp1);
	gwfree (gwdata, tmp3);
	gwfree (gwdata, tmp4);
}

/********************************************************************/
/* Helper routines for the multithreaded add/sub/smallmul/etc. code */
/********************************************************************/

/* Return a pointer to the big/lit flags within the first pass 1 block's variable data */

static __inline void *biglit_data_ptr (
	gwhandle *gwdata,
	int	i)
{
	char	*p;

	// This code only works on AVX-512 FFTs.  Earlier FFT implementations did not store big/lit flags in the variable data
	ASSERTG (gwdata->cpu_flags & CPU_AVX512F && gwdata->PASS2_SIZE);

	// Calculate address of biglit data in the first pass 1 block's variable data.
	p = (char *) gwdata->pass1_var_data + gwdata->biglit_data_offset;

	// Point to the correct set of big/lit flags
	if (gwdata->ZERO_PADDED_FFT) return (p + (i * gwdata->PASS1_CACHE_LINES * 2 / 8 / 2));
	else return (p + (i * gwdata->PASS1_CACHE_LINES * 2 / 8));
}

/* Structure to share add/sub/addsub/smallmul data amongst all the compute threads */

struct multithread_op_data {
	gwnum	s1;			/* Source #1 */
	gwnum	s2;			/* Source #2 */
	gwnum	d1;			/* Destination #1 */
	gwnum	d2;			/* Destination #2 (if operation is addsub) */
	void	(*asm_proc)(void*);	/* Pointer to assembly routine that will perform the operation */
	int	is_quick;		/* TRUE if this is a "quick" operation that does not require carry propagation */
	int	num_blks;		/* Number of "blocks" to process */
	void	*d1_carries;		/* Carries area for destination #1 calculations */
	void	*d2_carries;		/* Carries area for destination #2 calculations */
};

/* Perform a multithreaded add/sub/addsub/smallmul operation */

void multithread_op (
	gwhandle *gwdata,		/* Handle initialized by gwsetup */
	gwnum	s1,			/* Source #1 */
	gwnum	s2,			/* Source #2 */
	gwnum	d1,			/* Destination #1 */
	gwnum	d2,			/* Destination #2 (if operation is addsub) */
	void	(*asm_proc)(void*),	/* Pointer to assembly routine that will perform the operation */
	int	is_quick)		/* TRUE if this is a "quick" operation that does not require carry propagation */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	struct multithread_op_data data;

// Common data structure initialization for both AVX and AVX-512

	data.s1 = s1;
	data.s2 = s2;
	data.d1 = d1;
	data.d2 = d2;
	data.asm_proc = asm_proc;
	data.is_quick = is_quick;

/* Handle gwcopy for all architectures */
//BUG - is there a better alternative to memcpy esp. built with our MSVC 2005 environment?  Like using SSE/AVX/AVX-512 loads and stores.
//BUG - is there any particular advantage to aligning moves on 128, 256, or 4096 byte boundary?

	if (asm_proc == NULL) {
		int	size, movesize;

		size = ((uint32_t *) s1)[-2] + GW_HEADER_SIZE;
		data.s1 = (gwnum) ((char *) data.s1 - GW_HEADER_SIZE);
		data.d1 = (gwnum) ((char *) data.d1 - GW_HEADER_SIZE);

		// Move any bytes before a 128-byte boundary
		if ((intptr_t) data.s1 & 127) {
			movesize = 128 - (int) ((intptr_t) data.s1 & 127);
			memcpy (data.d1, data.s1, movesize);
			size -= movesize;
			data.s1 = (gwnum) ((char *) data.s1 + movesize);
			data.d1 = (gwnum) ((char *) data.d1 + movesize);
		}

		// Calculate number of 4KB blocks
		data.num_blks = size >> 12;

		// Move any bytes after last 4KB block
		if (size & 4095) {
			memcpy ((char *) data.d1 + data.num_blks * 4096, (char *) data.s1 + data.num_blks * 4096, size & 4095);
		}
	}

/* Initialize the structure to share necessary info amongst the AVX-512 compute threads */

	else if (gwdata->cpu_flags & CPU_AVX512F) {
		data.num_blks = asm_data->addcount1;
		data.d1_carries = asm_data->carries;
		// BUG/OPT - use scratch area rather than aligned_malloc?  one-pass FFTs are an issue.  We'd need to make
		// sure scratch area is big enough and pre-allocate area for one-pass FFTs
		if (d2 != NULL && !is_quick) data.d2_carries = aligned_malloc (data.num_blks * 16 * sizeof (double), 128);
	}

/* Initialize the structure to share necessary info amongst the AVX compute threads */

	else {
		data.num_blks = asm_data->addcount1;
//		data.d1_carries = asm_data->carries;
		// BUG/OPT - use scratch area rather than aligned_malloc?  one-pass FFTs are an issue.  We'd need to make
		// sure scratch area is big enough and pre-allocate area for one-pass FFTs
//		if (d2 != NULL && !is_quick) data.d2_carries = aligned_malloc (data.num_blks * 16 * sizeof (double), 128);
	}

/* Wake up the auxiliary threads */

	gwdata->pass1_state = PASS1_STATE_MULTITHREAD_OP;
	gwdata->multithread_op_data = &data;
	gwdata->next_block = 0;
	if (gwdata->num_threads > 1) {
		gwmutex_lock (&gwdata->thread_lock);
		gwdata->catch_straggler_threads = FALSE;
		gwevent_reset (&gwdata->all_threads_done);
		gwevent_signal (&gwdata->thread_work_to_do);
		gwmutex_unlock (&gwdata->thread_lock);
	}

/* Call subroutine to work on blocks just like the auxiliary threads */

	do_multithread_op_work (gwdata, asm_data);

/* Wait for auxiliary threads to finish */

	if (gwdata->num_threads > 1)
		gwevent_wait (&gwdata->all_threads_done, 0);

/* If no carry propagation is required then we're done */

	if (is_quick) return;

/* Do the AVX-512 carry propagation */

	if (gwdata->cpu_flags & CPU_AVX512F) {

/* Set flags used in processing carries */

		asm_data->zero_fft = 0;
		asm_data->const_fft = 0;
		asm_data->add_sub_smallmul_op = 1;

/* Do the final carry propagations for the dest #1 */

		asm_data->DESTARG = d1;		// Start of the FFT data
		asm_data->data_addr = d1;	// FFT addr to apply carries
		asm_data->norm_ptr1 = biglit_data_ptr (gwdata, 0);
		asm_data->this_block = 0;
		gwz3_apply_carries (asm_data);

/* Do the final carry propagations for the dest #2 */

		if (d2 != NULL) {
			asm_data->carries = data.d2_carries;
			asm_data->DESTARG = d2;		// Start of the FFT data
			asm_data->data_addr = d2;	// FFT addr to apply carries
			asm_data->norm_ptr1 = biglit_data_ptr (gwdata, 0);
			asm_data->this_block = 0;
			gwz3_apply_carries (asm_data);
			asm_data->carries = data.d1_carries;
			aligned_free (data.d2_carries);
		}
	}

/* Do the AVX carry propagation */

#ifdef AVX_ONLY_DOES_QUICK_RIGHT_NOW
	else {

/* Set flags used in processing carries */

		asm_data->zero_fft = 0;
		asm_data->const_fft = 0;
		asm_data->add_sub_smallmul_op = 1;

/* Do the final carry propagations for the dest #1 */

		asm_data->DESTARG = d1;		// Start of the FFT data
		asm_data->data_addr = d1;	// FFT addr to apply carries
		asm_data->norm_ptr1 = biglit_data_ptr (gwdata, 0);
		asm_data->this_block = 0;
		gwy3_apply_carries (asm_data);

/* Do the final carry propagations for the dest #2 */

		if (d2 != NULL) {
			asm_data->carries = data.d2_carries;
			asm_data->DESTARG = d2;		// Start of the FFT data
			asm_data->data_addr = d2;	// FFT addr to apply carries
			asm_data->norm_ptr1 = biglit_data_ptr (gwdata, 0);
			asm_data->this_block = 0;
			gwy3_apply_carries (asm_data);
			asm_data->carries = data.d1_carries;
			aligned_free (data.d2_carries);
		}
	}
#endif
}

/* Routine for the main thread and auxiliary threads to do add/sub/addsub/smallmul work */

void do_multithread_op_work (
	gwhandle *gwdata,		/* Handle initialized by gwsetup */
	struct gwasm_data *asm_data)	/* This thread's asm_data */
{
	struct multithread_op_data *data = (struct multithread_op_data *) gwdata->multithread_op_data;

/* Loop processing gwcopy blocks (4KB) */

	if (data->asm_proc == NULL) {
		for ( ; ; ) {
			int	i;

/* Get next block to process */

			if (gwdata->num_threads > 1) gwmutex_lock (&gwdata->thread_lock);
			i = gwdata->next_block;
			gwdata->next_block++;
			if (gwdata->num_threads > 1) gwmutex_unlock (&gwdata->thread_lock);

/* Break out of loop when there are no more blocks to process */

			if (i >= data->num_blks) break;

/* Move a 4KB block */

			memcpy ((char *) data->d1 + i * 4096, (char *) data->s1 + i * 4096, 4096);
		}
	}

/* Loop processing AVX-512 blocks */

	else if (gwdata->cpu_flags & CPU_AVX512F) {
		for ( ; ; ) {
			int	i;

/* Get next block to process */

			if (gwdata->num_threads > 1) gwmutex_lock (&gwdata->thread_lock);
			i = gwdata->next_block;
			gwdata->next_block += (data->is_quick ? 1 : 4);
			if (gwdata->num_threads > 1) gwmutex_unlock (&gwdata->thread_lock);

/* Break out of loop when there are no more blocks to process */

			if (i >= data->num_blks) break;

/* Process the block */

			asm_data->SRCARG = (char *) data->s1 + i * asm_data->pass1blkdst + (i/4) * asm_data->normblkdst4;
			asm_data->SRC2ARG = (char *) data->s2 + i * asm_data->pass1blkdst + (i/4) * asm_data->normblkdst4;
			asm_data->DESTARG = (char *) data->d1 + i * asm_data->pass1blkdst + (i/4) * asm_data->normblkdst4;
			asm_data->DEST2ARG = (char *) data->d2 + i * asm_data->pass1blkdst + (i/4) * asm_data->normblkdst4;
			asm_data->data_addr = (char *) data->d1_carries + i * 128;	// Addr to store dest #1 carries
			asm_data->premult_addr = (char *) data->d2_carries + i * 128;	// Addr to store dest #2 carries
			asm_data->norm_ptr1 = biglit_data_ptr (gwdata, i);
			call_op (data->asm_proc, asm_data);
		}
	}


/* Loop processing AVX blocks */

	else {
		for ( ; ; ) {
			int	i;

/* Get next block to process */

			if (gwdata->num_threads > 1) gwmutex_lock (&gwdata->thread_lock);
			i = gwdata->next_block;
			gwdata->next_block += (data->is_quick ? 1 : 2);
			if (gwdata->num_threads > 1) gwmutex_unlock (&gwdata->thread_lock);

/* Break out of loop when there are no more blocks to process */

			if (i >= data->num_blks) break;

/* Process the block */

			asm_data->SRCARG = (char *) data->s1 + i * asm_data->pass1blkdst;
			asm_data->SRC2ARG = (char *) data->s2 + i * asm_data->pass1blkdst;
			asm_data->DESTARG = (char *) data->d1 + i * asm_data->pass1blkdst;
			asm_data->DEST2ARG = (char *) data->d2 + i * asm_data->pass1blkdst;
//			asm_data->data_addr = (char *) data->d1_carries + i * 128;	// Addr to store dest #1 carries
//			asm_data->premult_addr = (char *) data->d2_carries + i * 128;	// Addr to store dest #2 carries
//			asm_data->norm_ptr1 = biglit_data_ptr (gwdata, i);
			call_op (data->asm_proc, asm_data);
		}
	}
}

/* Copy one gwnum to another gwnum */

void gwcopy (			/* Copy a gwnum */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d)		/* Dest */
{
	uint32_t free_offset;

/* Load the one piece of information that should not be copied over */

	free_offset = ((uint32_t *) d)[-8];

/* Copy the data and 96-byte header */

//if (rand() % 100 < 3) *s += 1.0;			// Generate random errors (for caller to test error recovery)
	if (gwdata->num_threads > 1)
		multithread_op (gwdata, s, NULL, d, NULL, NULL, TRUE);
	else
		memmove ((char *) d - GW_HEADER_SIZE, (char *) s - GW_HEADER_SIZE, ((uint32_t *) s)[-2] + GW_HEADER_SIZE);
//if (rand() % 100 < 5) *d += 1.0;			// Generate random errors (for caller to test error recovery)

/* Restore the one piece of information that should not be copied over */

	((uint32_t *) d)[-8] = free_offset;
}

/*********************************************************/
/* Wrapper routines for the add and sub assembly code    */
/*********************************************************/

void gwadd3quick (		/* Add two numbers without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* We could support completely/partially FFTed inputs if we updated the 7 zero pad */
/* values.  Until then, assert inputs are not completely/partially FFTed. */

	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-completely/partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* Do an AVX-512 or two-pass AVX addquick */

	if ((gwdata->cpu_flags & CPU_AVX512F) ||
	    (gwdata->cpu_flags & CPU_AVX && gwdata->PASS2_SIZE)) {
		multithread_op (gwdata, s1, s2, d, NULL, addr_gw_addq (gwdata), TRUE);
	}

/* Do the add the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s1;
		asm_data->SRC2ARG = s2;
		asm_data->DESTARG = d;
		gw_addq (gwdata, asm_data);
	}
}

void gwsub3quick (		/* Compute s1 - s2 without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* We could support completely/partially FFTed inputs if we updated the 7 zero pad */
/* values.  Until then, assert inputs are not partially FFTed. */

	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-completely/partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* Do an AVX-512 or two-pass AVX subquick */

	if ((gwdata->cpu_flags & CPU_AVX512F) ||
	    (gwdata->cpu_flags & CPU_AVX && gwdata->PASS2_SIZE)) {
		multithread_op (gwdata, s2, s1, d, NULL, addr_gw_subq (gwdata), TRUE);
	}

/* Do the subtract the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s2;
		asm_data->SRC2ARG = s1;
		asm_data->DESTARG = d;
		gw_subq (gwdata, asm_data);
	}
}

void gwaddsub4quick (		/* Add & sub two numbers without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);

/* We could support completely/partially FFTed inputs if we updated the 7 zero pad */
/* values.  Until then, assert inputs are not partially FFTed. */

	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Update the counts of unnormalized adds and subtracts */

	((uint32_t *) d1)[-1] =
	((uint32_t *) d2)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-completely/partially-FFTed flag */

	((uint32_t *) d1)[-7] =
	((uint32_t *) d2)[-7] = ((uint32_t *) s1)[-7];

/* Do an AVX-512 or two-pass AVX addsubquick */

	if ((gwdata->cpu_flags & CPU_AVX512F) ||
	    (gwdata->cpu_flags & CPU_AVX && gwdata->PASS2_SIZE)) {
		multithread_op (gwdata, s1, s2, d1, d2, addr_gw_addsubq (gwdata), TRUE);
	}

/* Do the add & subtract the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s1;
		asm_data->SRC2ARG = s2;
		asm_data->DESTARG = d1;
		asm_data->DEST2ARG = d2;
		gw_addsubq (gwdata, asm_data);
	}
}


void gwadd3 (			/* Add two numbers normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	uint32_t normcnt1, normcnt2;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */
/* If we can do an addquick, do that - it should be faster */

	normcnt1 = ((uint32_t *) s1)[-1];
	normcnt2 = ((uint32_t *) s2)[-1];
	if (normcnt1 + normcnt2 <= gwdata->EXTRA_BITS) {
		gwadd3quick (gwdata, s1, s2, d);
		return;
	}

/* Do an AVX-512 add */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		multithread_op (gwdata, s1, s2, d, NULL, addr_gw_add (gwdata), FALSE);
	}

/* Do the add the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s1;
		asm_data->SRC2ARG = s2;
		asm_data->DESTARG = d;
		gw_add (gwdata, asm_data);
	}

/* Clear the has-been-completely/partially-FFTed flag, reset the unnormalized adds and subtracts count */

	((uint32_t *) d)[-7] = 0;
	((uint32_t *) d)[-1] = 1;
}

void gwsub3 (			/* Compute s1 - s2 normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	uint32_t normcnt1, normcnt2;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */
/* If we can do a subquick, do that - it should be faster */

	normcnt1 = ((uint32_t *) s1)[-1];
	normcnt2 = ((uint32_t *) s2)[-1];
	if (normcnt1 + normcnt2 <= gwdata->EXTRA_BITS) {
		gwsub3quick (gwdata, s1, s2, d);
		return;
	}

/* Do an AVX-512 subtract */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		multithread_op (gwdata, s2, s1, d, NULL, addr_gw_sub (gwdata), FALSE);
	}

/* Do the add the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s2;
		asm_data->SRC2ARG = s1;
		asm_data->DESTARG = d;
		gw_sub (gwdata, asm_data);
	}

/* Clear the has-been-completely/partially-FFTed flag, reset the unnormalized adds and subtracts count */

	((uint32_t *) d)[-7] = 0;
	((uint32_t *) d)[-1] = 1;
}

void gwaddsub4 (		/* Add & sub two nums normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;
	uint32_t normcnt1, normcnt2;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == 0);
	ASSERTG (((uint32_t *) s2)[-7] == 0);

/* Get counts of unnormalized adds and subtracts */
/* If we can do an addsubquick, do that - it should be faster */

	normcnt1 = ((uint32_t *) s1)[-1];
	normcnt2 = ((uint32_t *) s2)[-1];
	if (normcnt1 + normcnt2 <= gwdata->EXTRA_BITS) {
		gwaddsub4quick (gwdata, s1, s2, d1, d2);
		return;
	}

/* Do an AVX-512 add/sub */

	if (gwdata->cpu_flags & CPU_AVX512F) {
		multithread_op (gwdata, s1, s2, d1, d2, addr_gw_addsub (gwdata), FALSE);
	}

/* Do the add & subtract the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s1;
		asm_data->SRC2ARG = s2;
		asm_data->DESTARG = d1;
		asm_data->DEST2ARG = d2;
		gw_addsub (gwdata, asm_data);
	}

/* Clear the has-been-completely/partially-FFTed flag, reset the unnormalized adds and subtracts count */

	((uint32_t *) d1)[-7] = 0;
	((uint32_t *) d2)[-7] = 0;
	((uint32_t *) d1)[-1] = 1;
	((uint32_t *) d2)[-1] = 1;
}


void gwfftadd3 (		/* Add two FFTed numbers */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == ((uint32_t *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-completely/partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* If this is a zero-padded FFT, then also add the 7 copied doubles in the gwnum header */

	if (gwdata->ZERO_PADDED_FFT) {
		d[-5] = s1[-5] + s2[-5];
		d[-6] = s1[-6] + s2[-6];
		d[-7] = s1[-7] + s2[-7];
		d[-8] = s1[-8] + s2[-8];
		d[-9] = s1[-9] + s2[-9];
		d[-10] = s1[-10] + s2[-10];
		d[-11] = s1[-11] + s2[-11];
	}

/* Do an AVX-512 or two-pass AVX addquick */

	if ((gwdata->cpu_flags & CPU_AVX512F) ||
	    (gwdata->cpu_flags & CPU_AVX && gwdata->PASS2_SIZE)) {
		multithread_op (gwdata, s1, s2, d, NULL, addr_gw_addf (gwdata), TRUE);
	}

/* Do the add the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s1;
		asm_data->SRC2ARG = s2;
		asm_data->DESTARG = d;
		gw_addf (gwdata, asm_data);
	}
}

void gwfftsub3 (		/* Compute FFTed s1 - FFTed s2 */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d)		/* Destination */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == ((uint32_t *) s2)[-7]);

/* Update the count of unnormalized adds and subtracts */

	((uint32_t *) d)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-completely/partially-FFTed flag */

	((uint32_t *) d)[-7] = ((uint32_t *) s1)[-7];

/* If this is a zero-padded FFT, then also subtract the 7 copied doubles in the gwnum header */

	if (gwdata->ZERO_PADDED_FFT) {
		d[-5] = s1[-5] - s2[-5];
		d[-6] = s1[-6] - s2[-6];
		d[-7] = s1[-7] - s2[-7];
		d[-8] = s1[-8] - s2[-8];
		d[-9] = s1[-9] - s2[-9];
		d[-10] = s1[-10] - s2[-10];
		d[-11] = s1[-11] - s2[-11];
	}

/* Do an AVX-512 or two-pass AVX subquick */

	if ((gwdata->cpu_flags & CPU_AVX512F) ||
	    (gwdata->cpu_flags & CPU_AVX && gwdata->PASS2_SIZE)) {
		multithread_op (gwdata, s2, s1, d, NULL, addr_gw_subf (gwdata), TRUE);
	}

/* Do the subtract the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s2;
		asm_data->SRC2ARG = s1;
		asm_data->DESTARG = d;
		gw_subf (gwdata, asm_data);
	}
}

void gwfftaddsub4 (		/* Add & sub two FFTed numbers */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2)		/* Destination #2 */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

	ASSERTG (((uint32_t *) s1)[-1] >= 1);
	ASSERTG (((uint32_t *) s2)[-1] >= 1);
	ASSERTG (((uint32_t *) s1)[-7] == ((uint32_t *) s2)[-7]);

/* Update the counts of unnormalized adds and subtracts */

	((uint32_t *) d1)[-1] =
	((uint32_t *) d2)[-1] = ((uint32_t *) s1)[-1] + ((uint32_t *) s2)[-1];

/* Copy the has-been-completely/partially-FFTed flag */

	((uint32_t *) d1)[-7] =
	((uint32_t *) d2)[-7] = ((uint32_t *) s1)[-7];

/* If this is a zero-padded FFT, then also add & sub the 7 copied doubles in */
/* the gwnum header.  Copy data to temporaries first in case s1, s2 pointers */
/* are equal to the d1, d2 pointers! */

	if (gwdata->ZERO_PADDED_FFT) {
		double	v1, v2;
		v1 = s1[-5]; v2 = s2[-5]; d1[-5] = v1 + v2; d2[-5] = v1 - v2;
		v1 = s1[-6]; v2 = s2[-6]; d1[-6] = v1 + v2; d2[-6] = v1 - v2;
		v1 = s1[-7]; v2 = s2[-7]; d1[-7] = v1 + v2; d2[-7] = v1 - v2;
		v1 = s1[-8]; v2 = s2[-8]; d1[-8] = v1 + v2; d2[-8] = v1 - v2;
		v1 = s1[-9]; v2 = s2[-9]; d1[-9] = v1 + v2; d2[-9] = v1 - v2;
		v1 = s1[-10]; v2 = s2[-10]; d1[-10] = v1+v2; d2[-10] = v1-v2;
		v1 = s1[-11]; v2 = s2[-11]; d1[-11] = v1+v2; d2[-11] = v1-v2;
	}

/* Do an AVX-512 or two-pass AVX addsubquick */

	if ((gwdata->cpu_flags & CPU_AVX512F) ||
	    (gwdata->cpu_flags & CPU_AVX && gwdata->PASS2_SIZE)) {
		multithread_op (gwdata, s1, s2, d1, d2, addr_gw_addsubf (gwdata), TRUE);
	}

/* Do the add & subtract the old way -- single-threaded in assembly code */

	else {
		asm_data->SRCARG = s1;
		asm_data->SRC2ARG = s2;
		asm_data->DESTARG = d1;
		asm_data->DEST2ARG = d2;
		gw_addsubf (gwdata, asm_data);
	}
}

/* Routine to add a small number to a gwnum.  Some day, */
/* I might optimize this routine for the cases where just one or two */
/* doubles need to be modified in the gwnum */

void gwsmalladd (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	addin,		/* Small value to add to g */
	gwnum	g)		/* Gwnum to add a value into */
{
	int	cant_handle;

/* Assert unnormalized add count valid, input not completely/partially FFTed. */

	ASSERTG (((uint32_t *) g)[-1] >= 1);
	ASSERTG (((uint32_t *) g)[-7] == 0);

/* AVX-512 one-pass implementation spreads carry over four words.  If there aren't many */
/* bits-per-word (e.g. M388 in a 128 word FFT) then you can end up with a giant carry */
/* into the last FFT word.  We may want to update this code for other FFT implementations. */

	cant_handle = (addin > 32768.0 || addin < -32768.0) &&
		      log2 (fabs (addin)) > 5.0 * log2 (gwdata->b) * gwdata->avg_num_b_per_word;

/* A simple brute-force implementation.  We put k > 1 numbers through */
/* here because multiplying the addin value by 1/k complicates matters. */
/* If this routine is ever used much, we can try to optimize this. */
/* We also make sure there is at least one b value per FFT word, so */
/* that any carry can successfully spread over 3 FFT words. */

	if (gwdata->GWPROCPTRS[8] == NULL || gwdata->k > 1.0 || gwdata->NUM_B_PER_SMALL_WORD < 1 || cant_handle) {
		gwnum	tmp;
		tmp = gwalloc (gwdata);
		if (addin >= 0.0) {
			dbltogw (gwdata, addin, tmp);
			gwaddquick (gwdata, tmp, g);
		} else {
			dbltogw (gwdata, -addin, tmp);
			gwsubquick (gwdata, tmp, g);
		}
		gwfree (gwdata, tmp);
	}

/* The assembler optimized version */

	else {
		struct gwasm_data *asm_data;

		asm_data = (struct gwasm_data *) gwdata->asm_data;
		asm_data->DESTARG = g;
		asm_data->DBLARG = addin;
		gw_adds (gwdata, asm_data);
	}

}

/* This routine multiplies a gwnum by a small positive value.  This lets us apply some */
/* optimizations that cannot be performed by a full FFT multiplication. */

void gwsmallmul (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	double	mult,		/* Small value to multiply g by */
	gwnum	g)		/* Gwnum to multiply */
{
	struct gwasm_data *asm_data = (struct gwasm_data *) gwdata->asm_data;

/* Assert unnormalized add count valid, input not completely/partially FFTed. */

	ASSERTG (((uint32_t *) g)[-1] >= 1);
	ASSERTG (((uint32_t *) g)[-7] == 0);

/* The x87 assembly version won't spread carries over multiple words. */
/* Also, the x87 assembly version won't guard against carries out of the */
/* critical top words in a zero-padded number. */
/* A simple brute-force implementation.  Note we cannot call gwmul as this */
/* will use the caller's last gwsetnormroutine value which could incorrectly */
/* multiply by a constant. */
	
	if (! (gwdata->cpu_flags & (CPU_AVX512F | CPU_AVX | CPU_SSE2)) && (mult > 1024.0 || gwdata->ZERO_PADDED_FFT)) {
		gwnum	tmp;
		tmp = gwalloc (gwdata);
		if (mult == 1.0);
		else if (mult == 2.0) {
			gwadd (gwdata, g, g);
		}
		else if (mult == 3.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwadd (gwdata, tmp, g);
		}
		else if (mult == 4.0) {
			gwaddquick (gwdata, g, g);
			gwadd (gwdata, g, g);
		}
		else if (mult == 5.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwaddquick (gwdata, tmp, tmp);
			gwadd (gwdata, tmp, g);
		}
		else if (mult == 6.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwaddquick (gwdata, tmp, g);
			gwadd (gwdata, g, g);
		}
		else if (mult == 8.0) {
			gwaddquick (gwdata, g, g);
			gwaddquick (gwdata, g, g);
			gwadd (gwdata, g, g);
		}
		else if (mult == 9.0) {
			gwadd3quick (gwdata, g, g, tmp);
			gwaddquick (gwdata, tmp, tmp);
			gwaddquick (gwdata, tmp, tmp);
			gwadd (gwdata, tmp, g);
		}
		else {
			dbltogw (gwdata, mult, tmp);
			gwfft (gwdata, tmp, tmp);
			asm_data->NORMRTN = gwdata->GWPROCPTRS[norm_routines + (gwdata->NORMNUM & 1)];
			raw_gwfftmul (gwdata, tmp, g);
		}
		gwfree (gwdata, tmp);
	}

/* Do an AVX-512 small mul */

	else if (gwdata->cpu_flags & CPU_AVX512F) {
		asm_data->DBLARG = mult;
		multithread_op (gwdata, g, NULL, g, NULL, addr_gw_muls (gwdata), FALSE);
	}

/* Do the small mul the old way -- single-threaded in assembly code */

	else {
		asm_data->DESTARG = g;
		asm_data->DBLARG = mult;
		gw_muls (gwdata, asm_data);
		((uint32_t *) g)[-1] = 1;
	}

/* If the number has gotten too large (high words should all be */
/* weighted -1 or 0) then emulate general mod with 2 multiplies */

	if (gwdata->GENERAL_MOD &&
	    (* (double *) ((char *) g + gwdata->GW_GEN_MOD_MAX_OFFSET) <= -2.0 ||
	     * (double *) ((char *) g + gwdata->GW_GEN_MOD_MAX_OFFSET) > 0.0))
		emulate_mod (gwdata, g);
}
