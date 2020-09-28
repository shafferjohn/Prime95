/*----------------------------------------------------------------------
| common.h
|
| This file contains handy #defines that I use in all my projects
| 
|  Copyright 2005-2020 Mersenne Research, Inc.
|  All Rights Reserved.
+---------------------------------------------------------------------*/

#ifndef _COMMON_H
#define _COMMON_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#define TRUE	1
#define FALSE	0

/* Define the world/group/owner read/write attributes for creating files */
/* I've always used 0666 in Unix (everyone gets R/W access), but MSVC 8 */
/* now refuses to work with that setting -- insisting on 0600 instead. */

#ifdef _WIN32
#define	CREATE_FILE_ACCESS	0600
#else
#define	CREATE_FILE_ACCESS	0666
#endif

/* Define the ASSERT macro I use while debugging */

#ifdef GDEBUG
#include <assert.h>
#define ASSERTG		assert
#else
#define ASSERTG(a)
#endif

/* Bit manipulation macros */

#define bitset(a,i)	{ ((char *)a)[(i) >> 3] |= (1 << ((i) & 7)); }
#define bitclr(a,i)	{ ((char *)a)[(i) >> 3] &= ~(1 << ((i) & 7)); }
#define bittst(a,i)	(((char *)a)[(i) >> 3] & (1 << ((i) & 7)))

/* Handy macros to improve readability (copied from gwnum.c) */

#define _log2(n)			(log ((double)(n)) / log ((double)2.0))
#define _log10(n)			(log ((double)(n)) / log ((double)10.0))
#define divide_rounding_up(a,b)		((a + (b) - 1) / (b))
#define divide_rounding_down(a,b)	((a) / (b))
#define round_up_to_multiple_of(a,b)	(divide_rounding_up (a, b) * (b))
#define round_down_to_multiple_of(a,b)	(divide_rounding_down (a, b) * (b))
#define is_multiple_of(a,b)		((a) % (b) == 0)
#define one_based_modulo(a,b)		((((a) - 1) % (b)) + 1)					// Return value between 1 and b
#define _intmin(a,b)			(((int)(a) < (int)(b)) ? (int)(a) : (int)(b))
#define _intmax(a,b)			(((int)(a) > (int)(b)) ? (int)(a) : (int)(b))

/* Define a "safe" strcpy.  The official C runtime library says that overlapping */
/* buffers produce undefined results.  This safe strcpy allows overlapping */
/* buffers by using memmove instead. */

#define safe_strcpy(d,s)	memmove (d, s, strlen (s) + 1)
#ifdef GDEBUG
#undef strcpy
//#define strcpy(d,s)	assert((d) >= ((s)+strlen(s)+1) || (s) >= (d)+strlen(s)+1), safe_strcpy(d,s)
#define strcpy(d,s)	debug_strcpy(d,s)
__inline char *debug_strcpy(char *d, const char *s) {assert((d) >= ((s)+strlen(s)+1) || (s) >= (d)+strlen(s)+1); return (char *) safe_strcpy(d,s); }
#endif

/* Routines missing from GMP */

#define mpz_add_si(d,s,addin)	if (addin >= 0) mpz_add_ui(d,s,(unsigned int)addin); else mpz_sub_ui(d,s,(unsigned int)-addin);
#define mpz_sub_si(d,s,addin)	if (addin >= 0) mpz_sub_ui(d,s,(unsigned int)addin); else mpz_add_ui(d,s,(unsigned int)-addin);
#define mpz_mul_d(d,s,flt)	{ mpz_t t; mpz_init_set_d(t,flt); mpz_mul(d,s,t); mpz_clear(t); }
#define mpz_eq(a,b)		(mpz_cmp(a,b) == 0)
#define mpz_eq_ui(a,b)		(mpz_cmp_ui(a,b) == 0)

#endif
