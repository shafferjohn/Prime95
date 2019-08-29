/*----------------------------------------------------------------------
| gwnum.h
|
| This file contains the headers and definitions that are used in the
| multi-precision IBDWT arithmetic routines.  That is, all routines
| that deal with the gwnum data type.
|
| Gwnums are great for applications that do a lot of multiplies modulo
| a number.  Only Intel x86-platforms are supported.  Add and subtract
| are also pretty fast.
|
| Gwnums are not suited to applications that need to convert to and from
| binary frequently or need to change the modulus frequently.
|
| MULTI-THREAD WARNING: You CAN perform gwnum operations in different
| threads IF AND ONLY IF each uses a different gwhandle structure
| initialized by gwinit.
| 
|  Copyright 2002-2019 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWNUM_H
#define _GWNUM_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "gwcommon.h"
#include "giants.h"
#include "gwthread.h"

/* To support multithreading, callers of the gwnum routines must allocate */
/* a gwhandle (on the heap or stack) and pass it to all gwnum routines. */
/* gwinit and gwsetup fill this structure up with lots of data that used to */
/* be stored in global variables. */

typedef struct gwhandle_struct gwhandle;

/* The gwnum data type.  A gwnum points to an array of doubles - the */
/* FFT data.  In practice, there is data stored before the doubles. */
/* See the internals section below if you really must know. */

typedef double *gwnum;

/*---------------------------------------------------------------------+
|                     SETUP AND TERMINATION ROUTINES                   |
+---------------------------------------------------------------------*/

/* This is the version number for the gwnum libraries. It changes whenever */
/* there is a change to the gwnum code.  Since Prime95 also uses the same */
/* version numbering scheme, you will see some strange jumps in gwnum version */
/* numbers when there are new prime95 versions without any changes in the gwnum code. */
/* This version number is also embedded in the assembly code and */
/* gwsetup verifies that the version numbers match.  This prevents bugs */
/* from accidentally linking in the wrong gwnum library. */

#define GWNUM_VERSION		"29.8"
#define GWNUM_MAJOR_VERSION	29
#define GWNUM_MINOR_VERSION	8

/* Error codes returned by the three gwsetup routines */

#define GWERROR_VERSION		1001	/* GWNUM.H and FFT assembly code version numbers do not match. */
#define GWERROR_TOO_LARGE	1002	/* Number too large for the FFTs */
#define GWERROR_K_TOO_SMALL	1003	/* k < 1 is not supported */
#define GWERROR_K_TOO_LARGE	1004	/* k > 53 bits is not supported */
#define GWERROR_MALLOC		1005	/* Insufficient memory available */
#define GWERROR_VERSION_MISMATCH 1006	/* GWNUM_VERSION from gwinit call doesn't match GWNUM_VERSION when */
					/* gwnum.c was compiled. */
#define GWERROR_STRUCT_SIZE_MISMATCH 1007 /* Gwhandle structure size from gwinit call doesn't match size */
					/* when gwnum.c was compiled.  Check compiler alignment switches. */
#define GWERROR_TOO_SMALL	1008	/* Gwsetup called on a number <= 1 */
#define GWERROR_NO_INIT		1009	/* gwinit was not called prior to gwsetup */
#define GWERROR_INTERNAL	2000	/* 2000 and up are "impossible" internal errors. */

/* Error codes returned by gwtobinary, gwtogiant, and get_fft_value */

#define GWERROR_BAD_FFT_DATA	-1	/* Nan or inf data encountered */
#define GWERROR_PARTIAL_FFT	-1009	/* Attempt to convert a partially FFTed number to binary */
#define GWERROR_FFT		-1010	/* Attempt to convert an FFTed number to binary */

/* Prior to calling gwsetup, you MUST CALL gwinit. This initializes the */
/* gwhandle structure. It gives us a place to set rarely used gwsetup */
/* options prior to calling gwsetup. */
#define gwinit(h)		gwinit2 (h, sizeof (gwhandle), GWNUM_VERSION)
/* The gwinit function has been superceeded by gwinit2.  By passing in the */
/* version number we can verify the caller used the same gwnum.h file as the */
/* one he eventually links with.  The sizeof (gwhandle) structure is used */
/* to verify he compiles with the same structure alignment options that */
/* were used when compiling gwnum.c.  For compatibility with existing code */
/* we delay reporting any compatibility problems until gwsetup is called. */
void gwinit2 (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	int	struct_size,	/* Size of the gwdata structure */
	const char *version_string);

/* There are three different setup routines.  The first, gwsetup, is for */
/* gwnum's primary use - support for fast operations modulo K*B^N+C. */
/* Smaller K and C values result in smaller FFT sizes and faster operations. */
/* Right now, if B<>2 defaults to the slower gwsetup_general_mod case. */

int gwsetup (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	double	k,		/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,	/* B in K*B^N+C. */
	unsigned long n,	/* N in K*B^N+C. Exponent to test. */
	signed long c);		/* C in K*B^N+C. Must be rel. prime to K. */

/* This setup routine is for operations modulo an arbitrary binary number. */
/* This is three times slower than the special forms above. */
/* The code will try to convert suitable k*2^n+c values into the faster */
/* gwsetup (gwdata,b,b,n,c) call above.  The caller would be better off */
/* not relying on this detection if at all possible. */

int gwsetup_general_mod (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	const uint32_t *array,	/* The modulus as an array of 32-bit values */
	uint32_t arraylen);	/* Number of values in the array */
int gwsetup_general_mod_64 (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	const uint64_t *array,	/* The modulus as an array of 64-bit values */
	uint64_t arraylen);	/* Number of values in the array */

/* This setup routine is for operations without a modulo. In essence, */
/* you are using gwnums as a general-purpose FFT multiply library. */

int gwsetup_without_mod (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	unsigned long n);	/* Maximum number of bits in OUTPUT numbers. */

/* Free all memory allocated by gwnum routines since gwsetup was called. */

void gwdone (
	gwhandle *gwdata);	/* Handle initialized by gwsetup */

/*---------------------------------------------------------------------+
|                    GWNUM OBSCURE GWSETUP OPTIONS                     |
+---------------------------------------------------------------------*/

/* Prior to calling one of the gwsetup routines, you can tell the library */
/* how many compute threads it can use to perform a multiply. */

#define gwset_num_threads(h,n)		((h)->num_threads = n)
#define gwget_num_threads(h)		((h)->num_threads)

/* Prior to calling one of the gwsetup routines, you can tell the library to use a hyperthread for memory prefetching. */
/* Only implemented for AVX-512 FFTs.  Caller must ensure the compute thread and prefetching hyperthread are set to use */
/* the same physical CPU core.  At present there are no known CPUs where this provides a benefit. */

#define gwset_hyperthread_prefetch(h)	((h)->hyperthread_prefetching = TRUE)
#define gwclear_hyperthread_prefetch(h)	((h)->hyperthread_prefetching = FALSE)

/* Specify a call back routine for the auxiliary threads to call when they */
/* are created.  This lets the user of the gwnum library set the thread */
/* priority and affinity as it sees fit.  You can also specify an arbitrary */
/* pointer to pass to the callback routine. */
/* The callback routine must be declared as follows: */
/*	void callback (int thread_num, int action, void *data) */
/* If you tell gwnum to use 4 threads, it will create 3 auxiliary threads */
/* and invoke the callback routine with thread_num = 1, 2, and 3. */
/* Action is 0 for thread starting and 1 for thread terminating. */
/* Action is 10 for prefetching hyperthread starting and 11 for prefetching hyperthread terminating. */

#define gwset_thread_callback(h,n)		((h)->thread_callback = n)
#define gwset_thread_callback_data(h,d)		((h)->thread_callback_data = d)

/* Prior to calling one of the gwsetup routines, you can have the library */
/* "play it safe" by reducing the maximum allowable bits per FFT data word. */
/* For example, the code normally tests a maximum of 22477 bits in a 1024 */
/* SSE2 FFT, or 21.95 bits per double.  If you set the safety margin to 0.5 */
/* then the code will only allow 21.45 bits per double, or a maximum of */
/* 21965 bits in a 1024 length FFT.  You can also use this option to */
/* "live dangerously" by increasing the maximum allowable bits per FFT */
/* data word - just set the safety margin to a negative value. */

#define gwset_safety_margin(h,m)	((h)->safety_margin = m)

/* The gwsetup routines need to know the maximum value that will be used */
/* in a call to gwsetmulbyconst.  By default this value is assumed to be 3, */
/* which is what you would use in a base-3 Fermat PRP test.  Gwsetup must */
/* switch to a generic modular reduction if k * mulbyconst or c * mulbyconst */
/* is too large.  Call this routine prior to calling gwsetup. */

#define gwset_maxmulbyconst(h,c)	((h)->maxmulbyconst = c)
#define gwsetmaxmulbyconst		gwset_maxmulbyconst

/* The gwsetup routines pick the fastest FFT implementation by default. */
/* Setting this option will cause gwsetup to give preference to FFT */
/* implementations that support the SUM(INPUTS) != SUM(OUTPUTS) error check. */
/* NOTE:  This error check is not available for k*b^n+c IBDWT FFTs when */
/* c is positive.  Setting this option will have no effect. */
/* NOTE: sum_inputs checking is only available in SSE2 FFTs and earlier. */

#define gwset_sum_inputs_checking(h,b) ((h)->sum_inputs_checking = (char) (b))

/* When doing a gwsetup_general_mod, the library prefers to use an */
/* integral number of bits per word (a rational FFT) because they are */
/* a little faster than irrational FFTs.  However, some moduli create */
/* non-random data when using rational FFTs.  For example, if we test */
/* (10^828809-1)/9 and put exactly 18 bits into each FFT word, then */
/* every FFT word in GW_MODULUS_FFT will contains the same value! */
/* Not exactly the random data the FFTs require for small roundoff errors. */
/* This routine takes a boolean to force use of the safer irrational FFTs. */

#define gwset_irrational_general_mod(h,b)  ((h)->use_irrational_general_mod = (char) (b))

/* Prior to calling one of the gwsetup routines, you can force the library */
/* to use a larger fft length than normal.  The input argument specifies */
/* how many FFT sizes larger than normal you would like.  You might use this */
/* routine if you are having roundoff errors using the normal FFT length. */

#define gwset_larger_fftlen_count(h,n)	((h)->larger_fftlen_count = n)

/* Prior to calling one of the gwsetup routines, you can set a minimum fft length. */
/* You might use this routine to select a larger FFT if you are having roundoff */
/* errors using the normal FFT length (or use gwset_larger_fftlen_count).  The library */
/* will use the first FFT meeting the minimum_fftlen criteria -- EVEN IF IT WOULD */
/* NOT ORDINARILY DO SO!!!  Set the FFT length below the default FFT length only if */
/* you know what you are doing!! */
#define gwset_minimum_fftlen(h,n)	((h)->minimum_fftlen = n)

/* Prior to calling one of the gwsetup routines, you can have the library */
/* use benchmark data stored in gwnum.txt to select the fastest */
/* FFT implementation.  This is the default behavior. */

#define gwset_use_benchmarks(h)		((h)->use_benchmarks = 1)
#define gwclear_use_benchmarks(h)	((h)->use_benchmarks = 0)

/* Set this if FFTs will use hyperthreading. This may affect selection of fastest */
/* FFT implementation.  By default, it is assumed hyperthreading will not be used. */

#define gwset_will_hyperthread(h,n)	((h)->will_hyperthread = n)
#define gwclear_will_hyperthread(h)	((h)->will_hyperthread = 0)

/* Set this if it is known how many cores will be used in total -- either by your */
/* program or multiple instances of your program.  By default, this value is */
/* the number of cores on the machine, which means the user of your program will */
/* keep all cores fully occupied with gwnum work.  This setting may affect */
/* selection of fastest FFT implementation and can be overriden in gwnum.txt. */

#define gwset_bench_cores(h,n)		((h)->bench_num_cores = n)

/* Set this if it is known how many independent gwnum FFTs will be active -- either by your */
/* program or multiple instances of your program.  Prime95 calls this "worker windows". */
/* By default, this value is the number of cores divided by number of threads.  This setting */
/* may affect selection of fastest FFT implementation and can be overriden in gwnum.txt. */

#define gwset_bench_workers(h,n)	((h)->bench_num_workers = n)

/* Set this if FFTs will always error check, will error check if near limit of FFT, or will */
/* not error check.  This setting may affect selection of fastest FFT implementation. */
/* By default, it is assumed round off error checking will not be used for every operation. */

#define gwset_will_error_check(h)		((h)->will_error_check = 1)
#define gwset_will_error_check_near_limit(h)	((h)->will_error_check = 2)
#define gwclear_will_error_check(h)		((h)->will_error_check = 0)

/* Prior to calling one of the gwsetup routines, you can have the library */
/* attempt to use large pages (2MB or 4MB on Intel architecture) rather than the */
/* standard 4KB pages.  This may improve performance by reducing TLB misses. */
/* It may have system-wide costs, as the OS may not page these to disk */
/* when not in use.  NOTE: Only the first gwalloc will return memory */
/* allocated using large pages. */

#define gwset_use_large_pages(h)	((h)->use_large_pages = 1)
#define gwclear_use_large_pages(h)	((h)->use_large_pages = 0)
#define gwget_use_large_pages(h)	((h)->use_large_pages)
#define gw_using_large_pages(h)		((h)->large_pages_ptr != NULL)

/* DEPRECATED, use gwset_minimum_fftlen instead. */
/* Prior to calling one of the gwsetup routines, you can force the library */
/* to use a specific fft length.  This should rarely (if ever) be used. */
/* I use it occasionally for benchmarking and/or checking round off errors */
/* at the FFT crossover points. */
/* Only choose a specific FFT size if you know what you are doing!! */
#define gwset_specific_fftlen(h,n)	((h)->minimum_fftlen = n)

/*---------------------------------------------------------------------+
|                     GWNUM MEMORY ALLOCATION ROUTINES                 |
+---------------------------------------------------------------------*/

/* Allocate memory for a gwnum */
gwnum gwalloc (
	gwhandle *gwdata);	/* Handle initialized by gwsetup */

/* Free a previously allocated gwnum */
void gwfree (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	val);		/* Gwnum to free */

/* Free all previously allocated gwnums */
void gwfreeall (
	gwhandle *gwdata);	/* Handle initialized by gwsetup */

/*---------------------------------------------------------------------+
|                        GWNUM CONVERSION ROUTINES                     |
+---------------------------------------------------------------------*/

/* Convert a double (must be an integer) to a gwnum */
void dbltogw (gwhandle *, double, gwnum);

/* Convert a binary value (array of 32-bit values) to a gwnum */
void binarytogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const uint32_t *array,	/* Array containing the binary value */
	uint32_t arraylen,	/* Length of the array */
	gwnum	n);		/* Destination gwnum */

/* Convert a binary value (array of 64-bit values) to a gwnum */
void binary64togw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const uint64_t *array,	/* Array containing the binary value */
	uint64_t arraylen,	/* Length of the array */
	gwnum	n);		/* Destination gwnum */

/* Convert a binary value (array of 32-bit or 64-bit values) to a gwnum. */
/* Check your C compiler specs to see if a long is 32 or 64 bits. */

void binarylongstogw (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	const unsigned long *array, /* Array containing the binary value */
	unsigned long arraylen,	/* Length of the array */
	gwnum	n);		/* Destination gwnum */

/* Convert a gwnum to a binary value (array of 32-bit values).  Returns */
/* the number of 32-bit values written to the array.  The array is NOT */
/* zero-padded.  Returns a negative number if an error occurs during the */
/* conversion.  An error can happen if the FFT data contains a NaN or */
/* infinity value. */
long gwtobinary (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	uint32_t *array,	/* Array to contain the binary value */
	uint32_t arraylen);	/* Maximum size of the array */

/* Convert a gwnum to a binary value (array of 64-bit values).  Returns */
/* the number of 64-bit values written to the array.  The array is NOT */
/* zero-padded.  Returns a negative number if an error occurs during the */
/* conversion.  An error can happen if the FFT data contains a NaN or */
/* infinity value. */
long gwtobinary64 (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	uint64_t *array,	/* Array to contain the binary value */
	uint64_t arraylen);	/* Maximum size of the array */

/* Convert a gwnum to a binary value (array of 32-bit or 64-bit values). */
/* Check your C compiler specs to see if a long is 32 or 64 bits. */
long gwtobinarylongs (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	n,		/* Source gwnum */
	unsigned long *array,	/* Array to contain the binary value */
	unsigned long arraylen);/* Maximum size of the array */

/* Generate a random number.  Can be useful for QA purposes. */
void gw_random_number (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	x);		/* Returned random number */

/*---------------------------------------------------------------------+
|                          GWNUM MATH OPERATIONS                       |
+---------------------------------------------------------------------*/

/* Macros to interface with assembly code */
/* The assembly routines are designed to provide a flexible way of */
/* multiplying two numbers.  If you will use a value in several multiplies */
/* you can perform the forward transform just once.  Furthermore, the */
/* multiply routines are tuned to allow one unnormalized addition prior */
/* to a multiply without introducing too much convolution error.  Thus: */
/* Legal:	gwaddquick (h, t1, t2); gwmul (h, t2, x); */
/* Legal:	gwfft (h, t1, t1); gwfft (h, t2, t2); gwfftadd (h, t1, t2); gwfftmul (h, t2, x); */
/* Not Legal:	gwaddquick (h, t1, t2); gwaddquick (h, y, x); gwmul (h, t2, x); */
/* Not Legal:	gwfft (h, t1, t1); gwfftadd (h, t1, t1); gwfftfftmul (h, t1, t1, t2); */

/* A brief description of each of the "gw" routines: */
/* gwswap	Quickly swaps two gw numbers */
/* gwcopy(s,d)	Copies gwnum s to d */
/* gwadd	Adds two numbers and normalizes result if necessary */
/* gwsub	Subtracts first number from second number and normalizes result if necessary */
/* gwadd3quick	Adds two numbers WITHOUT normalizing */
/* gwsub3quick	Subtracts second number from first WITHOUT normalizing */
/* gwadd3	Adds two numbers and normalizes them if necessary */
/* gwsub3	Subtracts second number from first number and normalizes result if necessary */
/* gwaddsub	Adds and subtracts 2 numbers (first+second and first-second) normalizes the results if necessary */
/* gwaddsub4	Like, gwaddsub but can store results in separate variables */
/* gwaddsub4quick Like, gwaddsub4 but will not do a normalize */
/* gwfft	Perform the forward Fourier transform on a number */
/* gwsquare	Multiplies a number by itself */
/* gwsquare2	Multiplies a number by itself, takes a source and destination */
/* gwsquare_carefully  Like gwsquare but uses a slower method that will */
/*		have a lower roundoff error even if input is non-random data */
/*		NOTE: Unlike gwsquare, input cannot have been partially FFTed */
/* gwsquare2_carefully  Like gwsquare_carefully but takes a source and destination */
/* gwmul(s,d)	Computes d=s*d.  NOTE: s is replaced by its FFT */
/* gwsafemul(s,d) Like gwmul but s is not replaced with its FFT */
/* gwfftmul(s,d) Computes d=s*d.  NOTE: s must have been previously FFTed */
/* gwfftfftmul(s1,s2,d) Computes d=s1*s2.  Both s1 and s2 must have been previously FFTed */
/* gwmul_carefully(s,d)  Like gwsafemul but uses a slower method that will */
/*		have a lower roundoff error even if input is non-random data */
/*		NOTE: Unlike gwsafemul, inputs cannot have been partially FFTed */

/* The routines below operate on numbers that have already been FFTed. */

/* gwfftadd	Adds two FFTed numbers */
/* gwfftsub	Subtracts first FFTed number from second FFTed number */
/* gwfftadd3	Adds two FFTed numbers */
/* gwfftsub3	Subtracts second FFTed number from first FFTed number */
/* gwfftaddsub	Adds and subtracts 2 FFTed numbers */
/* gwfftaddsub4	Like, gwfftaddsub but stores results in separate variables */

#define gwswap(s,d)	{gwnum t; t = s; s = d; d = t;}
#define gwsquare(h,s)	gwsquare2 (h,s,s)
#define gwsquare_carefully(h,s)	gwsquare2_carefully (h,s,s)
#define gwaddquick(h,s,d) gwadd3quick (h,s,d,d)
#define gwsubquick(h,s,d) gwsub3quick (h,d,s,d)
#define gwadd(h,s,d)	gwadd3 (h,s,d,d)
#define gwsub(h,s,d)	gwsub3 (h,d,s,d)
#define gwaddsub(h,a,b)	gwaddsub4 (h,a,b,a,b)
#define gwaddsubquick(h,a,b) gwaddsub4quick (h,a,b,a,b)
#define gwtouch(h,s)	gwcopy (h,s,s)
#define gwfftadd(h,s,d)	gwfftadd3 (h,s,d,d)
#define gwfftsub(h,s,d)	gwfftsub3 (h,d,s,d)
#define gwfftaddsub(h,a,b) gwfftaddsub4 (h,a,b,a,b)

/* Set the constant which the results of a multiplication should be */
/* multiplied by.  Use this macro in conjunction with the c argument of */
/* gwsetnormroutine. */

#define GWMULBYCONST_MAX	255		/* I think this is right */
void gwsetmulbyconst (gwhandle *gwdata, long s);

/* The multiplication code has two options that you can set using this */
/* macro.  The e argument tells the multiplication code whether or not */
/* it should perform round-off error checking - returning the maximum */
/* difference from an integer result in MAXERR.  The c argument tells the */
/* multiplication code whether or not it should multiply the result by */
/* a small constant. */

#define gwsetnormroutine(h,z,e,c) {(h)->NORMNUM=((c)?2:0)+((e)?1:0);}

/* If you know the result of a multiplication will be the input to another */
/* multiplication (but not gwsquare_carefully), then a small performance */
/* gain can be had in larger FFTs by doing some of the next forward FFT at */
/* the end of the multiplication.  Call this macro to tell the */
/* multiplication code whether or not it can start the forward FFT on */
/* the result. */

void gwstartnextfft (gwhandle *gwdata, int state);

void gwcopy (			/* Copy a gwnum */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d);		/* Dest */
void gwfft (			/* Forward FFT */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number */
	gwnum	d);		/* Destination (can overlap source) */
void gwsquare2 (		/* Square a number */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d);		/* Destination */
void gwmul (			/* Multiply source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number (changed to FFTed source!) */
	gwnum	d);		/* Source and destination */
void gwsafemul (		/* Multiply source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source number (not changed) */
	gwnum	d);		/* Source and destination */
void gwfftmul (			/* Multiply already FFTed source with dest */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Already FFTed source number */
	gwnum	d);		/* Non-FFTed source. Also destination */
void gwfftfftmul (		/* Multiply two already FFTed sources */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Already FFTed source number */
	gwnum	s2,		/* Already FFTed source number */
	gwnum	d);		/* Destination (can overlap sources) */
void gwadd3quick (		/* Add two numbers without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwsub3quick (		/* Compute s1 - s2 without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwaddsub4quick (		/* Add & sub two numbers without normalizing */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2);		/* Destination #2 */
void gwadd3 (			/* Add two numbers normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwsub3 (			/* Compute s1 - s2 normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwaddsub4 (		/* Add & sub two nums normalizing if needed */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2);		/* Destination #2 */
void gwfftadd3 (		/* Add two FFTed numbers */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwfftsub3 (		/* Compute FFTed s1 - FFTed s2 */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d);		/* Destination */
void gwfftaddsub4 (		/* Add & sub two FFTed numbers */
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s1,		/* Source #1 */
	gwnum	s2,		/* Source #2 */
	gwnum	d1,		/* Destination #1 */
	gwnum	d2);		/* Destination #2 */

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
	int	n);		/* Number of gwsquare calls to do carefully. */
				/* If n is -1, a default value is used */

/* Square or multiply numbers using a slower method that will have reduced */
/* round-off error on non-random input data.  Caller must make sure the */
/* input number has not been partially (via gwstartnextfft) or fully FFTed. */

void gwsquare2_carefully (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	d);		/* Destination */

void gwmul_carefully (
	gwhandle *gwdata,	/* Handle initialized by gwsetup */
	gwnum	s,		/* Source */
	gwnum	t);		/* Source and destination */


/* These routines can be used to add a constant to the result of a */
/* multiplication at virtually no cost.  Prime95 uses these routines to */
/* do the -2 operation in a Lucas-Lehmer test.  NOTE:  There are some */
/* number formats that cannot use these routines.  If abs(c) in k*b^n+c is 1, */
/* then gwsetaddin can be used.  To use gwsetaddinatpowerofb, k must also be 1. */
/* If you also use the mul-by-small-const normalization routine, the multiply */
/* is done after the addition. */

void gwsetaddin (gwhandle *, long);
void gwsetaddinatpowerofb (gwhandle *, long, unsigned long);

/* This routine adds a small value to a gwnum.  This lets us apply some */
/* optimizations that cannot be performed by general purpose gwadd */

#define GWSMALLADD_MAX		1125899906842624.0	/* 2^50 */
void gwsmalladd (gwhandle *gwdata, double addin, gwnum g);

/* This routine multiplies a gwnum by a small positive value.  This lets us apply some */
/* optimizations that cannot be performed by a full FFT multiplication. */

#define GWSMALLMUL_MAX		67108864.0		/* May allow more at a later date */
void gwsmallmul (gwhandle *gwdata, double mult, gwnum g);


/* DEPRECATED!!! These routines were deprecated because unlike all other gwnum routines */
/* the destination argument appeared before the source argument. */
#define gwaddsmall(h,g,a) gwsmalladd(h,a,g)	
#define gwmulsmall(h,g,m) gwsmallmul(h,m,g)
/* Replaced by better named gwsetaddinatpowerofb */
#define gwsetaddinatbit(h,v,b)	gwsetaddinatpowerofb(h,v,b)

/*-----------------------------------------------------------------+
|                      GWNUM COMPARISON ROUTINES                   |
+-----------------------------------------------------------------*/

/* Test if a gwnum is zero.  This routine was originally written by Jean Penne. */
/* It has not been adequately tested and MAY NOT BE BUG-FREE.  Use at your own risk! */
/* Returns TRUE if number is zero, FALSE if number is not zero, and a negative error */
/* code if a problem is found. */

int gwiszero (gwhandle *, gwnum);

/* Test two gwnums for equality.  Written by Jean Penne.  Uses the gwiszero routine */
/* which MAY NOT BE BUG-FREE.  Use this routine at your own risk! */
/* Returns TRUE if number is zero, FALSE if number is not zero, and a negative error */
/* code if a problem is found. */

int gwequal (gwhandle *, gwnum, gwnum);

/*---------------------------------------------------------------------+
|                      GWNUM ERROR-CHECKING ROUTINES                   |
+---------------------------------------------------------------------*/

#define gw_test_for_error(h)		((h)->GWERROR)
#define gw_test_illegal_sumout(h)	((h)->GWERROR & 1)
#define gw_test_mismatched_sums(h)	((h)->GWERROR & 2)
#define gwsuminp(h,g)			((g)[-2])
#define gwsumout(h,g)			((g)[-3])
#define gw_clear_error(h)		((h)->GWERROR = 0)

/* Get or clear the roundoff error.  Remember that if the roundoff error */
/* exceeds 0.5 then the FFT results will be wrong.  It is prudent to watch */
/* the roundoff error to make sure the roundoff error does not get close */
/* to 0.5. */

double gw_get_maxerr (gwhandle *gwdata);
void gw_clear_maxerr (gwhandle *gwdata);

/* Return TRUE if we are operating near the limit of this FFT length */
/* Input argument is the percentage to consider as near the limit. */
/* For example, if percent is 0.1 and the FFT can handle 20 bits per FFT */
/* data word, then if there are more than 19.98 bits per FFT data word */
/* this function will return TRUE. */

int gwnear_fft_limit (gwhandle *gwdata, double pct);

/*---------------------------------------------------------------------+
|                    GWNUM MISC. INFORMATION ROUTINES                  |
+---------------------------------------------------------------------*/

/* Map a gwerror code into human readable text */
void gwerror_text (gwhandle *gwdata, int error_code, char *buf, int buflen);

/* Return TRUE if this is a GPL'ed version of the GWNUM library. */
#define gwnum_is_gpl()		(0)

/* Return the FFT length being used */
#define gwfftlen(h)		((h)->FFTLEN)

/* Generate a human-readable description of the chosen FFT length and type */
void gwfft_description (gwhandle *, char *buf);

/* A human-readable string for the modulus currently in use */
#define gwmodulo_as_string(h)	((h)->GWSTRING_REP)

/* Get the number of threads gwnum can use */
#define gwget_num_threads(h)	((h)->num_threads)

/* Gwnum keeps a running count of the number of Fast Fourier transforms */
/* performed.  You can get and reset this counter. */

#define gw_get_fft_count(h)	((h)->fft_count)
#define gw_clear_fft_count(h)	((h)->fft_count = 0.0)

/* Get the amount of memory required for the gwnum's raw FFT data.  This */
/* does not include the GW_HEADER_SIZE bytes for the header or any pad */
/* bytes that might be allocated for alignment.  I see little need for */
/* a program to use this routine. */

unsigned long gwnum_datasize (gwhandle *);

/* Get the amount of memory likely to be allocated a gwnum.  This includes */
/* FFT data, headers, and pad bytes for alignment. */

unsigned long gwnum_size (gwhandle *);

/* Get the fixed amount of memory allocated during gwsetup.  Programs can */
/* use this and gwnum_size to determine working set size and act accordingly.*/
unsigned long gwmemused (gwhandle *);

/* Return TRUE if the gwnum value has been partially FFTed. */
#define gwnum_is_partially_ffted(h,g)	(((uint32_t *) g)[-7])

/*---------------------------------------------------------------------+
|                 ALTERNATIVE INTERFACES USING GIANTS                  |
+---------------------------------------------------------------------*/

/* The giants library from Dr. Richard Crandall, Perfectly Scientific, */
/* is used internally for a few infrequent operations.  It can optionally */
/* be used in the interfaces to convert between gwnum data type and binary. */
/* I do not recommend this.  There are many other faster and more robust */
/* libraries available. */

#include "giants.h"

/* Same as gwsetup_general_mod but uses giants instead of array of longs */
int gwsetup_general_mod_giant (
	gwhandle *gwdata,	/* Placeholder for gwnum global data */
	giant n);		/* The modulus */

/* Convert a giant to a gwnum */
void gianttogw (gwhandle *, giant, gwnum);

/* Convert a gwnum to a giant.  WARNING: Caller must allocate an array that */
/* is several words larger than the maximum result that can be returned. */
/* This is a gross kludge that lets gwtogiant use the giant for intermediate */
/* calculations.  Returns a negative number if an error occurs.  Returns */
/* zero on success. */
int gwtogiant (gwhandle *, gwnum, giant);

/*---------------------------------------------------------------------+
|          MISC. CONSTANTS YOU PROBABLY SHOULDN'T CARE ABOUT           |
+---------------------------------------------------------------------*/

/* The maximum value k * mulbyconst can be in a zero pad FFT.  Larger */
/* values must use generic modular reduction. */

#define MAX_ZEROPAD_K	2251799813685247.0	/* 51-bit k's are OK. */

/* The maximum value c * mulbyconst can be in a zero pad FFT.  Larger */
/* values must use generic modular reduction. */

#define MAX_ZEROPAD_C	8388607			/* 23-bit c's seem to work. */

/*---------------------------------------------------------------------+
|          SPECIAL ECM ROUTINE FOR GMP-ECM USING GWNUM LIBRARY         |
+---------------------------------------------------------------------*/

/* Return codes */

#define ES1_SUCCESS		0	/* Success, but no factor */
#define ES1_FACTOR_FOUND	1	/* Success, factor found */
#define ES1_CANNOT_DO_IT	2	/* This k,b,n,c cannot be handled */
#define ES1_MEMORY		3	/* Out of memory */
#define ES1_INTERRUPT		4	/* Execution interrupted */
#define ES1_CANNOT_DO_QUICKLY	5	/* Requires 3-multiply reduction */
#define ES1_HARDWARE_ERROR	6	/* An error was detected, most */
					/* likely a hardware error. */

/* Option codes */

#define ES1_DO_SLOW_CASE	0x1	/* Set this if ecmStage1 should do */
					/* slow 3-multiply reduction cases. */

/* INPUTS:

Input number (3 possibilities):

1) k,b,n,c set and num_being_factored_array = NULL.  k*b^n+c is factored.
2) k,b,n,c zero and num_being_factored_array set.  num_being_factored is
   worked on using generic 3-multiply reduction
3) k,b,n,c set and num_being_factored_array set.  num_being_factored is
   worked on - it must be a factor of k*b^n+c.

A_array, B1 are required

B1_done is optional.  Use it to resume a stage 1 calculation.

x_array, z_array is the starting point.  If z_array is not given, then
z is assumed to be one.

stop_check_proc is optional

options are defined above


   OUTPUTS:

On success:

   if z_array is NULL, then x_array is set to normalized point
   else x_array, z_array is set to the unnormalized point

On factor found:

   x_array is set to the factor found

On interrupt:

   B1_done is set to the last prime below B1 that was processed.
   If z_array != NULL (preferred) then x_array and z_array are set to the
current point.  The (x,z) point is not normalized because it will
be slow for large numbers.  This is unacceptable during system shutdown.
Caller must allocate x and z arrays large enough to hold any k*b^n+c value.
   If z_array == NULL, then a normalized x is returned. Caller must allocate
x array large enough to hold any value less than num_being_factored.

*/

int gwnum_ecmStage1_u32 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	uint32_t *num_being_factored_array, /* Number to factor */
	unsigned long num_being_factored_array_len,
	uint64_t B1,			/* Stage 1 bound */
	uint64_t *B1_done,		/* Stage 1 that is already done */
	uint32_t *A_array,		/* A - caller derives it from sigma */
	unsigned long A_array_len,
	uint32_t *x_array,		/* X value of point */
	unsigned long *x_array_len,
	uint32_t *z_array,		/* Z value of point */
	unsigned long *z_array_len,
	int	(*stop_check_proc)(int),/* Ptr to proc that returns TRUE */
					/* if user interrupts processing */
	unsigned long options);

int gwnum_ecmStage1_u64 (
	double	k,			/* K in K*B^N+C */
	unsigned long b,		/* B in K*B^N+C */
	unsigned long n,		/* N in K*B^N+C */
	signed long c,			/* C in K*B^N+C */
	uint64_t *num_being_factored_array, /* Number to factor */
	unsigned long num_being_factored_array_len,
	uint64_t B1,			/* Stage 1 bound */
	uint64_t *B1_done,		/* Stage 1 that is already done */
	uint64_t *A_array,		/* A - caller derives it from sigma */
	unsigned long A_array_len,
	uint64_t *x_array,		/* X value of point */
	unsigned long *x_array_len,
	uint64_t *z_array,		/* Z value of point */
	unsigned long *z_array_len,
	int	(*stop_check_proc)(int),/* Ptr to proc that returns TRUE */
					/* if user interrupts processing */
	unsigned long options);

/*---------------------------------------------------------------------+
|                             GWNUM INTERNALS                          |
+---------------------------------------------------------------------*/

/* This structure mimics a jmptable entry defined in the assembly code */
/* We use C code to read the entry and do lots of initialization. */

struct gwasm_jmptab {
	uint32_t max_exp;	/* Maximum exponent for this FFT len */
	uint32_t fftlen;	/* FFT length */
	float	timing;		/* Reference machine's time for a squaring */
	uint32_t flags;		/* Flags defined in mult.asm */
	void	*proc_ptr;	/* Ptr to assembly coded FFT routine */
	uint32_t mem_needed;	/* Memory needed */
	int32_t counts[8];
};

struct gwasm_alt_jmptab {	/* Used when pass 1 and pass 2 code is shared among FFT implementations */
	uint32_t max_exp;	/* Maximum exponent for this FFT len */
	uint32_t fftlen;	/* FFT length */
	float	timing;		/* Reference machine's time for a squaring */
	uint32_t flags;		/* Flags defined in mult.asm */
	void	*proc_ptr;	/* Ptr to assembly coded pass 1 FFT routine */
	void	*pass2_proc_ptr;/* Ptr to assembly coded pass 2 FFT routine */
	uint32_t mem_needed;	/* Memory needed */
	int32_t counts[8];
};

/* Structure for maintaining groups of blocks for each pass 1 thread to work on. */
/* Each thread wants to work on contiguous blocks for independent carry propagation. */

struct pass1_carry_sections {
	unsigned int start_block;	/* First block in section */
	unsigned int last_block;	/* Last block in section */
	unsigned int next_block;	/* First unassigned/unprocessed block */
	int	section_state;		/* Various states in processing this section -- see code */
	int	can_carry_into_next;	/* Flag indicating it is safe to propagate carries out of */
					/* the last block in this section into the next block */
	int	dependent_section;	/* Which section's carry out depends on this section finishing */
					/* processing of its first block to propagate carries into */
};

/* The FFT types currently implemented in assembly code */

#define FFT_TYPE_HOME_GROWN		0
#define FFT_TYPE_RADIX_4		1
#define FFT_TYPE_RADIX_4_DELAYED	2
#define FFT_TYPE_RADIX_4_DWPN		3	/* r4delay with partial normalization */

/* The gwhandle structure containing all of gwnum's "global" data. */

struct gwhandle_struct {

	/* Variables which affect gwsetup.  These are usually set by macros above. */
	double	safety_margin;		/* Reduce maximum allowable bits per FFT data word by this amount. */
	long	maxmulbyconst;		/* Gwsetup needs to know the maximum value the caller will use in */
					/* gwsetmulbyconst.  The default value is 3, commonly used */
					/* in a base-3 Fermat PRP test. */
	unsigned long minimum_fftlen;	/* Minimum fft length for gwsetup to use. */
	unsigned long num_threads;	/* Number of compute threads to use in multiply routines.  Default is obviously one. */
	char	hyperthread_prefetching; /* Set to true to launch a separate thread for prefetching.  Caller must set */
					/* affinity to make sure hyperthread and compute thread share the same physical core */
	char	larger_fftlen_count;	/* Force using larger FFT sizes.  This is a count of how many FFT sizes to "skip over". */
	char	sum_inputs_checking;	/* If possible, pick an FFT implementation that */
					/* supports the SUM(INPUTS) != SUM(OUTPUTS) error check. */
	char	force_general_mod;	/* Forces gwsetup_general_mod to not check for a k*2^n+c reduction */
	char	use_irrational_general_mod; /* Force using an irrational FFT when doing a general mod. */
					/* This is slower, but more immune to round off errors from */
					/* pathological bit patterns in the modulus. */
	char	use_large_pages;	/* FUTURE USE: Try to use 2MB/4MB pages */
	char	use_benchmarks;		/* Use benchmark data in gwnum.txt to select fastest FFT implementations */
	char	will_hyperthread;	/* Set if FFTs will use hyperthreading (affects select fastest FFT implementation) */
	char	will_error_check;	/* Set if FFTs will error check (affects select fastest FFT implementation) */
	char	unused_setup_flags[3];
	int	bench_num_cores;	/* Set to expected number of cores that will FFT (affects select fastest FFT implementation) */
	int	bench_num_workers;	/* Set to expected number of workers that will FFT (affects select fastest FFT implementation) */
	/* End of variables affecting gwsetup */

	double	k;			/* K in K*B^N+C */
	unsigned long b;		/* B in K*B^N+C */
	unsigned long n;		/* N in K*B^N+C */
	signed long c;			/* C in K*B^N+C */
	unsigned long FFTLEN;		/* The FFT size we are using */
	unsigned long PASS1_SIZE;	/* Number of real values FFTed in pass 1. */
	unsigned long PASS2_SIZE;	/* Number of complex values FFTed in pass 2. */
	int	cpu_flags;		/* Copy of CPU_FLAGS at time gwinit was called (just in case CPU_FLAGS changes) */
	char	ZERO_PADDED_FFT;	/* True if doing a zero pad FFT */
	char	ALL_COMPLEX_FFT;	/* True if using all-complex FFTs */
	char	RATIONAL_FFT;		/* True if bits per FFT word is integer */
	char	POSTFFT;		/* True if starting forward FFT on a result */
	char	GENERAL_MOD;		/* True if doing general-purpose mod as defined in gwsetup_general_mod */
	char	NO_PREFETCH_FFT;	/* True if this FFT does no prefetching */
	char	IN_PLACE_FFT;		/* True if this FFT is in-place (no scratch area) */
	char	UNUSED_CHARS[1];
	int	FFT_TYPE;		/* Home-grown, Radix-4, etc. */
	int	ARCH;			/* Architecture.  Which CPU type the FFT is optimized for. */
	void	(*GWPROCPTRS[16])(void*); /* Ptrs to assembly routines */
	giant	GW_MODULUS;		/* In the general purpose mod case, this is the number operations are modulo */
	gwnum	GW_MODULUS_FFT;		/* In the general purpose mod case, this is  the FFT of GW_MODULUS */
	gwnum	GW_RECIP_FFT;		/* FFT of shifted reciprocal of GW_MODULUS */
	unsigned long GW_ZEROWORDSLOW;	/* Count of words to zero during copy step of a general purpose mod */
	unsigned long GW_GEN_MOD_MAX;	/* Maximum number of words we can safely allow in a GENERAL_MOD number */
	unsigned long GW_GEN_MOD_MAX_OFFSET; /* Offset to the GW_GEN_MOD_MAX word */
	unsigned long NUM_B_PER_SMALL_WORD; /* Number of b's in a little word.  For the common case, b=2, this is the */
					/* number of bits in a little word. */
	double	avg_num_b_per_word;	/* Number of base b's in each fft word */
	double	bit_length;		/* Bit length of k*b^n */
	double	fft_max_bits_per_word;	/* Maximum bits per data word that this FFT size can support */
	long	FOURKBGAPSIZE;		/* Gap between 4KB blocks in pass 2 of a two-pass FFT. */
					/* Number of cache lines in a padded block in a one-pass FFT. */
	long	PASS2GAPSIZE;		/* Gap between blocks in pass 2 */
	unsigned long PASS1_CACHE_LINES; /* Cache lines grouped together in first pass of an FFT */
	unsigned long mem_needed;	/* Memory needed for sin/cos, weights, etc. */
	unsigned long SCRATCH_SIZE;	/* Size of the pass 1 scratch area */
	unsigned long EXTRA_BITS;	/* Number of unnormalized adds that can be safely performed */
	gwnum	GW_RANDOM;		/* A random number used in gwsquare_carefully. */
	unsigned long saved_copyz_n;	/* Used to reduce COPYZERO calculations */
	char	GWSTRING_REP[60];	/* The gwsetup modulo number as a string. */
	unsigned int NORMNUM;		/* The post-multiply normalize routine index */
	int	GWERROR;		/* Set if an error is detected */
	double	MAXDIFF;		/* Maximum allowable difference between sum of inputs and outputs */
	double	fft_count;		/* Count of forward and inverse FFTs */
	const struct gwasm_jmptab *jmptab; /* ASM jmptable pointer */
	void	*asm_data;		/* Memory allocated for ASM global data */
	void	*dd_data;		/* Memory allocated for gwdbldbl global data */
	ghandle	gdata;			/* Structure that allows sharing giants and gwnum memory allocations */
	double	*gwnum_memory;		/* Allocated memory */
	unsigned long GW_ALIGNMENT;	/* How to align allocated gwnums */
	unsigned long GW_ALIGNMENT_MOD; /* How to align allocated gwnums */
	gwnum	*gwnum_alloc;		/* Array of allocated gwnums */
	unsigned int gwnum_alloc_count; /* Count of allocated gwnums */
	unsigned int gwnum_alloc_array_size; /* Size of gwnum_alloc array */
	gwnum	*gwnum_free;		/* Array of available gwnums */
	unsigned int gwnum_free_count;	/* Count of available gwnums */
	size_t	GW_BIGBUF_SIZE;		/* Size of the optional buffer */
	char	*GW_BIGBUF;		/* Optional buffer to allocate gwnums in */
	void	*large_pages_ptr;	/* Pointer to the lage pages memory block we allocated. */
	void	*large_pages_gwnum;	/* Pointer to the one large pages gwnum */
	void	(*thread_callback)(int, int, void *); /* Auxiliary thread callback routine letting */
					/* the gwnum library user set auxiliary thread priority and affinity */
	void	*thread_callback_data;	/* User-supplied data to pass to the auxiliary thread callback routine */
	unsigned int num_active_threads; /* Count of the number of active auxiliary threads */
	gwmutex	thread_lock;		/* This mutex limits one thread at a time in critical sections. */
	gwevent	thread_work_to_do;	/* This event is set whenever the auxiliary threads have work to do. */
	gwevent	all_threads_done;	/* This event is set whenever the auxiliary threads are done and the */
					/* main thread can resume.  That is, it is set if and only if num_active_threads==0 */
	gwevent can_carry_into;		/* This event signals pass 1 sections that the block they are waiting on to carry */
					/* into may now be ready. */
	short	threads_must_exit;	/* Flag set to force all auxiliary threads to terminate */
	short	catch_straggler_threads;/* Flag set when auxiliary threads have finished their work */
	int	pass1_state;		/* Mainly used to keep track of what we are doing in pass 1 of an FFT.  See */
					/* pass1_get_next_block for details.  Also, 999 means we are in pass 2 of the FFT. */
	void	*pass1_var_data;	/* pass1 variable sin/cos/premultiplier/fudge/biglit data */
	void	*adjusted_pass2_premults; /* pass2_premults pointer adjusted for the fact the first block of real FFTs have */
					/* no premultipliers */
	unsigned long biglit_data_offset; /* Offset of the big/lit data in the pass 1 variable data */
	unsigned long pass1_var_data_size; /* Used to calculate address of pass 1 premultiplier data */
	unsigned long pass2_premult_block_size; /* Used to calculate address of pass 2 premultiplier data */
	unsigned long next_block;	/* Next block for threads to process */
	unsigned long num_pass1_blocks; /* Number of data blocks in pass 1 for threads to process */
	unsigned long num_pass2_blocks; /* Number of data blocks in pass 2 for threads to process */
	unsigned long num_postfft_blocks; /* Number of data blocks that must delay forward fft because POSTFFT is set. */
	gwthread *thread_ids;		/* Array of auxiliary thread ids */
	struct pass1_carry_sections *pass1_carry_sections; /* Array of pass1 sections for carry propagation */
	void	*multithread_op_data;	/* Data shared amongst add/sub/addsub/smallmul compute threads */
	uint32_t ASM_TIMERS[32];	/* Internal timers used by me to optimize code */
	int	bench_pick_nth_fft;	/* DO NOT set this variable.  Internal hack to force the FFT selection code to */
					/* pick the n-th possible implementation instead of the best one.  The prime95 */
					/* benchmarking code uses this to time every FFT implementation. */
	int	qa_pick_nth_fft;	/* DO NOT set this variable.  Internal hack to force the FFT selection code to */
					/* pick the n-th possible implementation instead of the best one.  The prime95 QA */
					/* code uses this to compare results from one FFT implementation to the (should */
					/* be identical) results of another FFT implementation. */
	int	qa_picked_nth_fft;	/* Internal hack returning which FFT was picked */
	int	square_carefully_count; /* Count of gwsquare calls to convert into gwsquare_carefully calls */
	double	ZPAD_COPY7_ADJUST[7];	/* Adjustments for copying the 7 words around the halfway point of a zero pad FFT. */
	double	ZPAD_0_6_ADJUST[7];	/* Adjustments for ZPAD0_6 in a r4dwpn FFT */
	unsigned long wpn_count;	/* Count of r4dwpn pass 1 blocks that use the same ttp/ttmp grp multipliers */
};

/* A psuedo declaration for our big numbers.  The actual pointers to */
/* these big numbers are to the data array.  The 96 bytes prior to the */
/* data contain: */
/* data-4:  integer containing number of unnormalized adds that have been */
/*	    done.  After a certain number of unnormalized adds, the next add */
/*	    must be normalized to avoid overflow errors during a multiply. */
/* data-8:  integer containing number of bytes in data area. Used by gwcopy. */
/* data-16: double containing the product of the two sums of the input FFT */
/*	    values. */
/* data-24: double containing the sum of the output FFT values.  These two */
/*	    values can be used as a sanity check when multiplying numbers. */
/*	    The two values should be "reasonably close" to one another. */
/* data-28: Flag indicating gwnum value has been partially FFTed. */
/* data-32: Pointer returned by malloc - used to free memory when done. */
/* data-88: Seven doubles (input FFT values near the halfway point */
/*	    when doing a zero-padded FFT). */
/* data-96: Eight unused bytes */
/* typedef struct { */
/*	char	pad[96];	   Used to track unnormalized add/sub */
/*				   and original address */
/*	double	data[512];	   The big number broken into chunks */
/*				   This array is variably sized. */
/* } *gwnum; */
#define GW_HEADER_SIZE	96	/* Number of data bytes before a gwnum ptr */

/* Some mis-named #defines that describe the maximum Mersenne number */
/* exponent that the gwnum routines can process. */

#define MAX_PRIME	79300000L	/* Maximum number of x87 bits */
#define MAX_PRIME_SSE2	595800000L	/* SSE2 bit limit */
#define MAX_PRIME_AVX	595700000L	/* AVX bit limit */
#define MAX_PRIME_FMA3	922600000L	/* FMA3 bit limit */
#define MAX_PRIME_AVX512 1169000000L	/* AVX-512 bit limit */
#define MAX_FFTLEN	4194304L	/* 4M FFT max for x87 */
#define MAX_FFTLEN_SSE2	33554432L	/* 32M FFT max for SSE2 */
#define MAX_FFTLEN_AVX	33554432L	/* 32M FFT max for AVX */
#define MAX_FFTLEN_FMA3	52428800L	/* 50M FFT max for FMA3 */
#define MAX_FFTLEN_AVX512 67108864L	/* 64M FFT max for AVX-512 */

/* Informational routines that can be called prior to gwsetup */
/* Many of these routines only work for k*b^n+c FFTs. */

unsigned long gwmap_to_fftlen (double, unsigned long, unsigned long, signed long);
unsigned long gwmap_with_cpu_flags_to_fftlen (int, double, unsigned long, unsigned long, signed long);
double gwmap_to_timing (double, unsigned long, unsigned long, signed long);
unsigned long gwmap_to_memused (double, unsigned long, unsigned long, signed long);
unsigned long gwmap_fftlen_to_max_exponent (unsigned long fftlen);
unsigned long gwmap_with_cpu_flags_fftlen_to_max_exponent (int, unsigned long fftlen);
unsigned long gwmap_to_estimated_size (double, unsigned long, unsigned long, signed long);
//int gwmap_to_fft_info (gwhandle *, double, unsigned long, unsigned long, signed long); /* DEPRECATED */

/* Generate a human-readable string for k*b^n+c */
void gw_as_string(char *buf, double k, unsigned long b, unsigned long n, signed long c);

/* Other routines used internally */

int gwinfo (gwhandle *, double, unsigned long, unsigned long, signed long);
double virtual_bits_per_word (gwhandle *);
unsigned long addr_offset (gwhandle *, unsigned long);
double *addr (gwhandle *, gwnum, unsigned long);
int get_fft_value (gwhandle *, gwnum, unsigned long, long *);
void set_fft_value (gwhandle *, gwnum, unsigned long, long);
int is_big_word (gwhandle *, unsigned long);
void bitaddr (gwhandle *, unsigned long, unsigned long *, unsigned long *);
void specialmodg (gwhandle *, giant);
#define gw_set_max_allocs(h,n)	if ((h)->gwnum_alloc==NULL) (h)->gwnum_alloc_array_size=n

/* Specialized routines that let the internal giants code share the free */
/* memory pool used by gwnums. */

void gwfree_temporarily (gwhandle *, gwnum);
void gwrealloc_temporarily (gwhandle *, gwnum);

/* Routines to share the memory of cached free gwnums with giants code. */
/* Used by prime95 to have the giants GCD code reuse the memory used */
/* during P-1 and ECM calculations. */

void *gwgiantalloc (void *);
void gwgiantfree (void *, void *);
void gwgiantdealloc (void *);

/* When debugging gwnum and giants, I sometimes write code that "cheats" */
/* by calling a routine that is part of prime95 rather than the gwnum and */
/* giants library.  Prime95 will set this routine pointer so that gwnum */
/* code can cheat while keeping the gwnum library interface clean. */

extern void (*OutputBothRoutine)(int, const char *);

/* These routines let me time many assembly language building blocks -- used */
/* when optimizing these building blocks. */

int gwtimeit (void *);
#define get_asm_timers(g) ((uint32_t *) &(g)->ASM_TIMERS)

#ifdef __cplusplus
}
#endif

#endif
