/**************************************************************
 *
 *  gwdbldbl.cpp
 *
 *  This file contains all the gwnum initialization routines that require
 *  extended precision floating point.  We want to initialize our sin/cos
 *  and FFT weight arrays with doubles that are as accurate as possible.
 *
 *  This is the only C++ routine in the gwnum library.  Since gwnum is
 *  a C based library, we declare all routines here as extern "C".
 *
 *  Copyright 2005-2018 Mersenne Research, Inc.  All rights reserved.
 *
 **************************************************************/

/* Include files */

#include "gwdbldbl.h"

/* Pick which doubledouble package we will use. */

#define QD
//#define KEITH_BRIGGS

/* Turn on #define that will disable extended precision floating point */
/* registers, as they wreak havoc with double-double library routines. */

#ifndef X86_64
#define x86
#endif

/* Use Hida, Li & Bailey's QD doubledouble C++ package. */

#ifdef QD
#include "dd.cc"
#endif

/* Use Keith Briggs' doubledouble C++ package.  I find the QD package a better choice. */

#ifdef KEITH_BRIGGS
#define DD_INLINE
#include "doubledouble.h"
#include "doubledouble.cc"
#include "math.cc"
#define	dd_real doubledouble
#define	_2pi TwoPi
#define	_log2 Log2
#endif

/* Epsilon value, 2^-250, should have an exact representation as a double */

#define epsilon 5.5271478752604445602472651921923E-76

/* Structure for "global" data that gwfft_weight_setup and passes */
/* to several gwdbldbl routines */
	
struct gwdbldbl_constants {
	dd_real	gw__b;
	dd_real	gw__logb;
	dd_real	gw__num_b_per_word;
	int	gw__c_is_one;
	dd_real gw__logb_abs_c_div_fftlen;
	dd_real gw__fftlen_inverse;
	dd_real gw__over_fftlen;
	dd_real gw__over_sqrt_fftlen;
	double	gwdbl__b;
	double	gwdbl__b_inverse;
	double	gwdbl__fractional_b_per_word;
	double	gwdbl__logb_abs_c_div_fftlen;
	unsigned long integral_b_per_word;
	unsigned long last_inv_sloppy_j;
	double	last_inv_sloppy_result;
	double	fast_inv_sloppy_multiplier;
	double	last_partial_sloppy_power[2];
	double	last_partial_sloppy_result[2];
	double	last_partial_inv_sloppy_power[2];
	double	last_partial_inv_sloppy_result[2];
};

/* Macro routines below use to type the cast untyped data pointer */
/* input argument. */

#define dd_data		((struct gwdbldbl_constants *) dd_data_arg)

/* Now write all the routines that use the dd_real package. */

/* This routine is used by me to analyze the rounding error due to converting from dd_real to double. */
/* Sometimes I have a choice of how to arrange sin/cos values in the building block macros. */
/* I use this routine to select an arrangement that minimizes roundoff errors. */

double dderr (dd_real x) {
	double	y = double (x);
	double	err = double (fabs (x - y));
	double	absx = double (fabs (x));
	double	pcterr = err / absx * pow (2, 53.0);
	double	rounderr = err * pow (2.0, 53.0 - ceil (log(absx)/log(2.0)));
	return (y);
}

/* Utility routine to compute many of the constants used by the assembly language code */

extern "C"
void gwasm_constants (
	double	*asm_values)
{
	dd_real arg, sine1, cosine1, sine2, cosine2, sine3, cosine3;
#define P951			asm_values[0]
#define P618			asm_values[1]
#define P309			asm_values[2]
#define M262			asm_values[3]
#define P588			asm_values[4]
#define M162			asm_values[5]
#define M809			asm_values[6]
#define M382			asm_values[7]
#define P866			asm_values[8]
#define P433			asm_values[9]
#define P577			asm_values[10]
#define P975			asm_values[11]
#define P434_P975		asm_values[12]
#define P180			asm_values[13]
#define P623			asm_values[14]
#define M358			asm_values[15]
#define P404			asm_values[16]
#define M223			asm_values[17]
#define M901			asm_values[18]
#define M691			asm_values[19]
#define P924			asm_values[20]
#define P383			asm_values[21]
#define P782			asm_values[22]
#define P434			asm_values[23]
#define P975_P434		asm_values[24]
#define P782_P434		asm_values[25]
#define P259			asm_values[26]
#define P966			asm_values[27]
#define P259_P707		asm_values[28]
#define P966_P707		asm_values[29]
#define P924_P383		asm_values[30]
#define P981_P195		asm_values[31]
#define P981			asm_values[32]
#define P195			asm_values[33]
#define P831_P556		asm_values[34]
#define P831			asm_values[35]
#define P556			asm_values[36]
#define P556_P195		asm_values[37]
#define P901_P975		asm_values[38]
#define P623_P975		asm_values[39]
#define P223_P975		asm_values[40]
#define P1_P975			asm_values[41]
#define P782_P975		asm_values[42]

/* Do some initial setup */

	x86_FIX

/* Initialize the five-complex sine/cosine data. */
/* NOTE: When computing cosine / sine, divide by the 64-bit sine not the */
/* extra precision sine since macros will multiply by the 64-bit sine. */

	arg = dd_real::_2pi / 5.0;		// 2*PI * 1 / 5
	sincos (arg, sine1, cosine1);

	arg = arg * 2.0;			// 2*PI * 2 / 5
	sincos (arg, sine2, cosine2);

	P951 = double (sine1);
	P618 = double (sine2 / P951);		// 0.588 / 0.951

	P309 = double (cosine1);
	M262 = double (cosine2 / P309);		// -0.809 / 0.309

	P588 = double (sine2);
	M162 = double (-sine1 / P588);		// -0.951 / 0.588

	M809 = double (cosine2);
	M382 = double (cosine1 / M809);		// 0.309 / -0.809

/* Initialize the six_reals sine-cosine data. */

	arg = dd_real::_2pi / 3.0;		// 2*PI / 3
	sine1 = sin (arg);			// Compute sine (0.866)

	P866 = double (sine1);
	P433 = double (sine1 * 0.5);		// 0.5 * P866
	P577 = double (dd_real (0.5) / sine1);	// 0.5 / 0.866

/* Initialize the seven-complex sine/cosine data. */

	arg = dd_real::_2pi / 7.0;		// 2*PI * 1 / 7
	sincos (arg, sine1, cosine1);		// cosine (0.623), sine (0.782)

	arg = arg * 2.0;			// 2*PI * 2 / 7
	sincos (arg, sine2, cosine2);		// cosine (-.223), sine (0.975)

	arg = arg * 1.5;			// 2*PI * 3 / 7
	sincos (arg, sine3, cosine3);		// cosine (-.901), sine (0.434)

	P782 = double (sine1);
	P975 = double (sine2);
	P434 = double (sine3);
	P623 = double (cosine1);
	M223 = double (cosine2);
	M901 = double (cosine3);

	P434_P975 = double (sine3 / P975);	// 0.434 / 0.975
	P180 = double (sine1 / P975 / P434_P975); // 0.782 / (0.975 * 0.445)
	M358 = double (cosine2 / P623);		// -0.223 / 0.623
	P404 = double (cosine3 / P623 / M358);	// -.901 / (.623 * -.358)
	M691 = double (cosine1 / M901);		// 0.623 / -0.901
	P782_P434 = double (sine1 / P434);	// 0.782 / 0.434
	P975_P434 = double (sine2 / P434);	// 0.975 / 0.434

	P782_P975 = double (sine1 / P975);	// Radix-7 constants for AVX-512
	P901_P975 = double (-cosine3 / P975);
	P623_P975 = double (cosine1 / P975);
	P223_P975 = double (-cosine2 / P975);
	P1_P975 = 1.0 / P975;

/* Initialize the 24-reals sine-cosine data. */

	arg = dd_real::_2pi / 24.0;		// 2*PI * 1 / 24
	sincos (arg, sine1, cosine1);		// cosine (0.966), sine (0.259)

	P259 = double (sine1);
	P966 = double (cosine1);

	P259_P707 = double (sine1 / sqrt (0.5));
	P966_P707 = double (cosine1 / sqrt (0.5));

/* Initialize the roots of -1 used by r4_four_complex_first_fft and zr8_eight_complex_first_fft. */

	arg = dd_real::_2pi / 16.0;		// 2*PI / 16
	sincos (arg, sine1, cosine1);		// cosine (0.924), sine (0.383)
	P924 = double (cosine1);
	P383 = double (sine1);
	P924_P383 = double (cosine1 / P383);

	arg = dd_real::_2pi / 32.0;		// 2*PI / 32
	sincos (arg, sine1, cosine1);		// cosine (0.981), sine (0.195)
	P981 = double (cosine1);
	P195 = double (sine1);
	P981_P195 = double (cosine1 / P195);

	arg = dd_real::_2pi / 32.0;		// 2*PI / 32
	arg = 3.0 * arg;			// 3 * 2*PI / 32
	sincos (arg, sine1, cosine1);		// cosine (0.831), sine (0.556)
	P831 = double (cosine1);
	P556 = double (sine1);
	P831_P556 = double (cosine1 / P556);
	P556_P195 = double (sine1 / P195);

	END_x86_FIX
}

// Utility routines to compute a sin/cos multipliers

extern "C"
void gwsincos (
	unsigned long x,
	unsigned long N,
	double	*results)
{
	gwsincos1by (x, N, results, 1);
}

extern "C"
void gwsincos1by (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr)
{
	gwsincos12345by (x, N, results, incr, 1);
}

extern "C"
void gwsincos12345by (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr,
	int	amt)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine3, cosine3, sine4, cosine4, sine5, cosine5, sine6, cosine6, sine7, cosine7;
	dd_real	sine8, cosine8;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	results[0] = sine;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine / results[0];
	results += incr + incr;
	if (amt == 1) goto done;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);
	results[0] = sine2;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine2 / results[0];
	results += incr + incr;
	if (amt == 2) goto done;

	sine3 = sine * cosine2 + sine2 * cosine;
	cosine3 = cosine * cosine2 - sine * sine2;
	results[0] = sine3;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine3 / results[0];
	results += incr + incr;
	if (amt == 3) goto done;

	sine4 = sine2 * cosine2 * 2.0;
	cosine4 = sqr (cosine2) - sqr (sine2);
	results[0] = sine4;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine4 / results[0];
	results += incr + incr;
	if (amt == 4) goto done;

	sine5 = sine * cosine4 + sine4 * cosine;
	cosine5 = cosine * cosine4 - sine * sine4;
	results[0] = sine5;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine5 / results[0];
	results += incr + incr;
	if (amt == 5) goto done;

	sine6 = sine2 * cosine4 + sine4 * cosine2;
	cosine6 = cosine2 * cosine4 - sine2 * sine4;
	results[0] = sine6;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine6 / results[0];
	results += incr + incr;
	if (amt == 6) goto done;

	sine7 = sine3 * cosine4 + sine4 * cosine3;
	cosine7 = cosine3 * cosine4 - sine3 * sine4;
	results[0] = sine7;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine7 / results[0];
	results += incr + incr;
	if (amt == 7) goto done;

	sine8 = sine4 * cosine4 + sine4 * cosine4;
	cosine8 = cosine4 * cosine4 - sine4 * sine4;
	results[0] = sine8;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine8 / results[0];
done:	;
	END_x86_FIX
}

extern "C"
void gwsincos13579by (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr,
	int	amt)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine3, cosine3, sine5, cosine5, sine7, cosine7, sine9, cosine9, sine11, cosine11;
	dd_real sine13, cosine13, sine15, cosine15;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	results[0] = sine;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine / results[0];
	results += incr + incr;
	if (amt == 1) goto done;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);

	sine3 = sine * cosine2 + sine2 * cosine;
	cosine3 = cosine * cosine2 - sine * sine2;

	results[0] = sine3;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine3 / results[0];
	results += incr + incr;
	if (amt == 2) goto done;

	sine5 = sine3 * cosine2 + sine2 * cosine3;
	cosine5 = cosine3 * cosine2 - sine3 * sine2;

	results[0] = sine5;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine5 / results[0];
	results += incr + incr;
	if (amt == 3) goto done;

	sine7 = sine5 * cosine2 + sine2 * cosine5;
	cosine7 = cosine5 * cosine2 - sine5 * sine2;

	results[0] = sine7;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine7 / results[0];
	results += incr + incr;
	if (amt == 4) goto done;

	sine9 = sine7 * cosine2 + sine2 * cosine7;
	cosine9 = cosine7 * cosine2 - sine7 * sine2;

	results[0] = sine9;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine9 / results[0];
	results += incr + incr;
	if (amt == 5) goto done;

	sine11 = sine9 * cosine2 + sine2 * cosine9;
	cosine11 = cosine9 * cosine2 - sine9 * sine2;

	results[0] = sine11;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine11 / results[0];
	results += incr + incr;
	if (amt == 6) goto done;

	sine13 = sine11 * cosine2 + sine2 * cosine11;
	cosine13 = cosine11 * cosine2 - sine11 * sine2;

	results[0] = sine13;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine13 / results[0];
	results += incr + incr;
	if (amt == 7) goto done;

	sine15 = sine13 * cosine2 + sine2 * cosine13;
	cosine15 = cosine13 * cosine2 - sine13 * sine2;

	results[0] = sine15;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine15 / results[0];
done:	;
	END_x86_FIX
}

extern "C"
void gwsincos159by (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr,
	int	amt)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine4, cosine4, sine5, cosine5, sine9, cosine9, sine13, cosine13;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	results[0] = sine;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine / results[0];
	results += incr + incr;
	if (amt == 1) goto done;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);

	sine4 = sine2 * cosine2 * 2.0;
	cosine4 = sqr (cosine2) - sqr (sine2);

	sine5 = sine * cosine4 + sine4 * cosine;
	cosine5 = cosine * cosine4 - sine * sine4;

	results[0] = sine5;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine5 / results[0];
	results += incr + incr;
	if (amt == 2) goto done;

	sine9 = sine5 * cosine4 + sine4 * cosine5;
	cosine9 = cosine5 * cosine4 - sine5 * sine4;

	results[0] = sine9;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine9 / results[0];
	results += incr + incr;
	if (amt == 3) goto done;

	sine13 = sine9 * cosine4 + sine4 * cosine9;
	cosine13 = cosine9 * cosine4 - sine9 * sine4;

	results[0] = sine13;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine13 / results[0];
done:	;
	END_x86_FIX
}

extern "C"
void gwsincos125by (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine4, cosine4, sine5, cosine5;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	results[0] = sine;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine / results[0];
	results += incr + incr;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);

	results[0] = sine2;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine2 / results[0];
	results += incr + incr;

	sine4 = sine2 * cosine2 * 2.0;
	cosine4 = sqr (cosine2) - sqr (sine2);

	sine5 = sine * cosine4 + sine4 * cosine;
	cosine5 = cosine * cosine4 - sine * sine4;

	results[0] = sine5;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine5 / results[0];
	END_x86_FIX
}

extern "C"
void gwsincos1234by_raw (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr,
	int	amt)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine3, cosine3, sine4, cosine4;

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	results[0] = sine;
	results[incr] = cosine;
	results += incr + incr;
	if (amt == 1) goto done;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);
	results[0] = sine2;
	results[incr] = cosine2;
	results += incr + incr;
	if (amt == 2) goto done;

	sine3 = sine * cosine2 + sine2 * cosine;
	cosine3 = cosine * cosine2 - sine * sine2;
	results[0] = sine3;
	results[incr] = cosine3;
	results += incr + incr;
	if (amt == 3) goto done;

	sine4 = sine2 * cosine2 * 2.0;
	cosine4 = sqr (cosine2) - sqr (sine2);
	results[0] = sine4;
	results[incr] = cosine4;
done:	;
	END_x86_FIX
}


extern "C"
void gwsincos1plusby (
	unsigned long x,
	unsigned long inc,
	unsigned long N,
	double	*results,
	int	incr,
	int	amt)
{
	dd_real twopi_over_N, sineinc, cosineinc, sine0, cosine0;
	int	num_output;

	x86_FIX
	twopi_over_N = dd_real::_2pi / (double) N;
	sincos (twopi_over_N * (double) inc, sineinc, cosineinc);
	sincos (twopi_over_N * (double) x, sine0, cosine0);

	for (num_output = 0; ; ) {
		dd_real last_sine;
		/* Output the sine and cosine/sine value while protecting against divide by zero */
		results[0] = sine0;
		results[0] += epsilon;
		results[incr] = cosine0 / results[0];
		results += incr + incr;
		num_output++;
		/* Break if we've output the requested amount */
		if (num_output == amt) break;
		/* Compute next cosine and sine value */
		last_sine = sine0;
		sine0 = sine0 * cosineinc + cosine0 * sineinc;
		cosine0 = cosine0 * cosineinc - last_sine * sineinc;
	}
	END_x86_FIX
}

/* Like gwsincos1plusby but only the cos/sin value is output */

extern "C"
void gwcos1plusby (
	unsigned long x,
	unsigned long inc,
	unsigned long N,
	double	*results,
	int	incr,
	int	amt)
{
	dd_real twopi_over_N, sineinc, cosineinc, sine0, cosine0;
	int	num_output;

	x86_FIX
	twopi_over_N = dd_real::_2pi / (double) N;
	sincos (twopi_over_N * (double) inc, sineinc, cosineinc);
	sincos (twopi_over_N * (double) x, sine0, cosine0);

	for (num_output = 0; ; ) {
		dd_real last_sine;
		/* Output the cosine / sine value while protecting against divide by zero */
		results[0] = cosine0 / (sine0 + epsilon);
		results += incr;
		num_output++;
		/* Break if we've output the requested amount */
		if (num_output == amt) break;
		/* Compute next cosine and sine value */
		last_sine = sine0;
		sine0 = sine0 * cosineinc + cosine0 * sineinc;
		cosine0 = cosine0 * cosineinc - last_sine * sineinc;
	}
	END_x86_FIX
}

/* Special version for 8-complex and 16-real macros that multiplies some sine values by SQRTHALF */

#ifdef TRY_SQRT2_TO_REDUCE_ROUNDOFF
extern "C"
void gwsincos1234by_sqrthalf (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine3, cosine3, sine4, cosine4;
	dd_real sqrthalf;

	x86_FIX
	sqrthalf = sqrt (dd_real (0.5));
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	results[0] = sine * sqrthalf;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine * sqrthalf / results[0];
	results += incr + incr;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);
	results[0] = sine2;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine2 / results[0];
	results += incr + incr;

	sine3 = sine * cosine2 + sine2 * cosine;
	cosine3 = cosine * cosine2 - sine * sine2;
	results[0] = sine3 * sqrthalf;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine3 * sqrthalf / results[0];
	results += incr + incr;

	sine4 = sine2 * cosine2 * 2.0;
	cosine4 = sqr (cosine2) - sqr (sine2);
	results[0] = sine4;
	results[0] += epsilon;		/* Protect against divide by zero */
	results[incr] = cosine4 / results[0];
	END_x86_FIX
}
#endif

/* Special version for 14-reals (and 7-complex) macro multiplies sine values by 0.975^(2/3) */

extern "C"
void gwsincos123by_special7 (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine3, cosine3, tmp;
	dd_real sine975, multiplier;

	arg1 = dd_real::_2pi * 2.0 / 7.0;		// 2*PI * 2 / 7
	sine975 = sin (arg1);				// sine (0.975)
	multiplier = exp (log (sine975) * 2.0 / 3.0);	// 0.975 ^ (2/3)

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	tmp = sine + epsilon;		/* Protect against divide by zero */
	results[0] = tmp * multiplier;
	results[incr] = cosine / tmp;
	results += incr + incr;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);
	tmp = sine2 + epsilon;		/* Protect against divide by zero */
	results[0] = tmp * multiplier;
	results[incr] = cosine2 / tmp;
	results += incr + incr;

	sine3 = sine * cosine2 + sine2 * cosine;
	cosine3 = cosine * cosine2 - sine * sine2;
	tmp = sine3 + epsilon;		/* Protect against divide by zero */
	results[0] = tmp * multiplier;
	results[incr] = cosine3 / tmp;
	END_x86_FIX
}

/* Special version for 14-reals (and 7-complex) macro multiplies sine values by 0.975^(2/3) */

extern "C"
void gwsincos135by_special7 (
	unsigned long x,
	unsigned long N,
	double	*results,
	int	incr)
{
	dd_real arg1, sine, cosine, sine2, cosine2, sine3, cosine3, sine5, cosine5, tmp;
	dd_real sine975, multiplier;

	arg1 = dd_real::_2pi * 2.0 / 7.0;		// 2*PI * 2 / 7
	sine975 = sin (arg1);				// sine (0.975)
	multiplier = exp (log (sine975) * 2.0 / 3.0);	// 0.975 ^ (2/3)

	x86_FIX
	arg1 = dd_real::_2pi * (double) x / (double) N;
	sincos (arg1, sine, cosine);
	tmp = sine + epsilon;		/* Protect against divide by zero */
	results[0] = tmp * multiplier;
	results[incr] = cosine / tmp;
	results += incr + incr;

	sine2 = sine * cosine * 2.0;
	cosine2 = sqr (cosine) - sqr (sine);

	sine3 = sine * cosine2 + sine2 * cosine;
	cosine3 = cosine * cosine2 - sine * sine2;
	tmp = sine3 + epsilon;		/* Protect against divide by zero */
	results[0] = tmp * multiplier;
	results[incr] = cosine3 / tmp;
	results += incr + incr;

	sine5 = sine3 * cosine2 + sine2 * cosine3;
	cosine5 = cosine3 * cosine2 - sine3 * sine2;
	tmp = sine5 + epsilon;		/* Protect against divide by zero */
	results[0] = tmp * multiplier;
	results[incr] = cosine5 / tmp;
	END_x86_FIX
}


//
// Utility routines to compute fft weights
//
// The FFT weight for the j-th FFT word doing a b^n+c weighted transform is
//	b ^ (ceil (j*n/FFTLEN) - j*n/FFTLEN)   *   abs(c) ^ j/FFTLEN
//
// NOTE:  Simply using the dd_real ceil function is not a good idea.  This is
// because when we calculate (j * (n/FFTLEN)) and the result is an exact integer,
// then if FFTLEN is not a power of 2, the dd_real result is exact integer +/- epsilon.
// If it is +epsilon then the ceil function returns the wrong value for us.
// There are numerous examples, one being testing M26000208 using 1440K FFT.
// We also cannot simply convert to a real value before applying the ceil function
// because sometimes we need more than 53 bits of precision as happens when
// testing 10223*2^29588045-1.
// Our solution is to subtract a smidge before applying the ceil function.
// Since j*n/FFTLEN is always less than FFTLEN and we don't envision handling
// FFTs above 64M, the integer part the input to ceil will be 26 bits or less.
// We think dd_real routines support 106 bits of precision, leaving a fractional
// part accurate to at least 2^-80 or about 8e-25.  Thus, if we define "smidge" as
// 1e-22 then smidge will swamp +/- epsilon and force ceil to give us the correct
// answer while still providing a LOT more than 53 bits of precision for handling
// j*n/FFTLEN values that are very, very close to an integer.
// An alternative solution would be to not precompute n/FFTLEN in our calculations
// and assume that the dd_real routines will return an exact integer when appropriate
// for ((j*n) / FFTLEN).

inline dd_real gwceil (dd_real val)
{
	return (ceil (val - 1.0e-22));
}

inline dd_real gwfloor (dd_real val)
{
	return (floor (val + 1.0e-22));
}

extern "C"
void *gwdbldbl_data_alloc (void)
{
	return (malloc (sizeof (struct gwdbldbl_constants)));
}

extern "C"
void gwfft_weight_setup (
	void	*dd_data_arg,
	int	zero_pad,
	double	k,
	unsigned long b,
	unsigned long n,
	signed long c,
	unsigned long fftlen)
{
	x86_FIX
	dd_data->gw__b = dd_real ((double) b);
	dd_data->gw__logb = log (dd_real ((double) b));
	dd_data->gw__fftlen_inverse = dd_real (1.0) / dd_real ((double) fftlen);
	if (zero_pad) {
		dd_data->gw__num_b_per_word = dd_real ((double) (n + n)) * dd_data->gw__fftlen_inverse;
		dd_data->gw__c_is_one = TRUE;
		dd_data->gw__over_fftlen = dd_real (2.0) * dd_data->gw__fftlen_inverse;
	} else {
		dd_data->gw__num_b_per_word = (dd_real ((double) n) + log (dd_real (k)) / dd_data->gw__logb) * dd_data->gw__fftlen_inverse;
		dd_data->gw__c_is_one = (abs ((int) c) == 1);
		dd_data->gw__logb_abs_c_div_fftlen = log (dd_real (abs ((int) c))) / dd_data->gw__logb * dd_data->gw__fftlen_inverse;
		dd_data->gw__over_fftlen = dd_real (k * 2.0) * dd_data->gw__fftlen_inverse;
	}
	dd_data->gw__over_sqrt_fftlen = sqrt (dd_data->gw__over_fftlen);
	dd_data->gwdbl__b = (double) b;
	dd_data->gwdbl__b_inverse = 1.0 / (double) b;
	dd_data->integral_b_per_word = (unsigned long) double (gwfloor (dd_data->gw__num_b_per_word));
	dd_data->gwdbl__fractional_b_per_word = (double) (dd_data->gw__num_b_per_word - gwfloor (dd_data->gw__num_b_per_word));
	dd_data->gwdbl__logb_abs_c_div_fftlen = (double) dd_data->gw__logb_abs_c_div_fftlen;
	dd_data->last_inv_sloppy_j = 0;
	dd_data->last_inv_sloppy_result = 1.0;
	dd_data->fast_inv_sloppy_multiplier = gwfft_weight_inverse (dd_data, 1);

	dd_data->last_partial_sloppy_power[0] = dd_data->last_partial_sloppy_power[1] = 999.0;
	dd_data->last_partial_inv_sloppy_power[0] = dd_data->last_partial_inv_sloppy_power[1] = 999.0;
	END_x86_FIX
}

extern "C"
double gwfft_weight (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	result = exp (dd_data->gw__logb * bpower);
	END_x86_FIX
	return (double (result));
}

// Like the above, but multiplies the weight by sqrt (2/FFTLEN).

extern "C"
double gwfft_weight_over_sqrt_fftlen (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	result = exp (dd_data->gw__logb * bpower) * dd_data->gw__over_sqrt_fftlen;
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.
// We cannot be very sloppy as these weights are used to set FFT data values
// when reading save files.  If we are too sloppy we get roundoff errors > 0.5.

extern "C"
double gwfft_weight_sloppy (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	END_x86_FIX
	return (pow (dd_data->gwdbl__b, double (bpower)));
}

// Compute the inverse of the fft weight

extern "C"
double gwfft_weight_inverse (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	result = exp (dd_data->gw__logb * -bpower);
	END_x86_FIX
	return (double (result));
}

// Compute the inverse of the fft weight squared

extern "C"
double gwfft_weight_inverse_squared (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	result = exp (dd_data->gw__logb * -2.0 * bpower);
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.
// We can be very sloppy as these weights are used to read FFT data values
// when writing save files.

extern "C"
double gwfft_weight_inverse_sloppy (
	void	*dd_data_arg,
	unsigned long j)
{
	double	temp, bpower, result;

// Our sequential sloppy optimizations won't work if abs(c) is not one.
// This is because the result is not in the range 0.5 to 1.0.

	if (! dd_data->gw__c_is_one) {
		temp = (double) j * dd_data->gwdbl__fractional_b_per_word;
		bpower = ceil (temp) - temp;
		if (bpower < 0.001 || bpower > 0.999) return (gwfft_weight_inverse (dd_data_arg, j));
		bpower += dd_data->gwdbl__logb_abs_c_div_fftlen * (double) j;
		return (pow (dd_data->gwdbl__b, - bpower));
	}

// Compute weight from previous weight, but don't do too many of
// these in a row as floating point roundoff errors will accumulate

	if (j == dd_data->last_inv_sloppy_j + 1 && (j & 0x7F)) {
		result = dd_data->last_inv_sloppy_result * dd_data->fast_inv_sloppy_multiplier;
		if (result <= dd_data->gwdbl__b_inverse) result = result * dd_data->gwdbl__b;

// Just to be safe, if result is at all close to the boundaries return
// the carefully computed weight.

		if (result < dd_data->gwdbl__b_inverse + 0.001 || result > 0.999)
			result = gwfft_weight_inverse (dd_data_arg, j);
	}

// Get an accurate value to start off a series of computing next weights

	else {
		result = gwfft_weight_inverse (dd_data_arg, j);
	}

// Save the result for faster sequential sloppy calls

	dd_data->last_inv_sloppy_j  = j;
	dd_data->last_inv_sloppy_result = result;
	return (result);
}

// This computes the inverse FFT weight multiplied by the appropriate constant
// to produce an integer during an FFT multiply's normalize stage.  This
// constant is 2/FFTLEN for a zero-padded FFT and k*2/FFTLEN for a
// non-zero-padded FFT.

extern "C"
double gwfft_weight_inverse_over_fftlen (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	result = exp (dd_data->gw__logb * -bpower) * dd_data->gw__over_fftlen;
	END_x86_FIX
	return (double (result));
}

// This computes the three FFT weights in one call.  It is faster than calling the above individually.

extern "C"
void gwfft_weights3 (
	void	*dd_data_arg,
	unsigned long j,
	double	*fft_weight,
	double	*fft_weight_inverse,
	double	*fft_weight_inverse_over_fftlen)
{
	dd_real temp, bpower, weight;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	weight = exp (dd_data->gw__logb * bpower);
	if (fft_weight != NULL) *fft_weight = double (weight);
	if (fft_weight_inverse != NULL) *fft_weight_inverse = double (1.0 / weight);
	if (fft_weight_inverse_over_fftlen != NULL) *fft_weight_inverse_over_fftlen = double (dd_data->gw__over_fftlen / weight);
	END_x86_FIX
}

// Returns logb(fft_weight).  This is used in determining the FFT weight
// fudge factor in two-pass FFTs.  This is much faster than computing the
// fft_weight because it eliminates a call to the double-double exp routine.

extern "C"
double gwfft_weight_exponent (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, bpowerdbl;
	dd_real temp, bpower;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * dd_data->gwdbl__fractional_b_per_word;
	bpowerdbl = ceil (tempdbl) - tempdbl;
	if (bpowerdbl > 0.001 && bpowerdbl < 0.999) return (bpowerdbl);

// We cannot be certain of the ceiling result, use doubledoubles to do the calculation

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	END_x86_FIX
	return (double (bpower));
}

// Funky version that does not factor in the weight for c != 1.  Required for building the AVX-512 XOR masks in the compressed fudge table,

extern "C"
double gwfft_weight_no_c (
	void	*dd_data_arg,
	unsigned long j)
{
	dd_real temp, bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	result = exp (dd_data->gw__logb * bpower);
	END_x86_FIX
	return (double (result));
}

//
// Utility routine to compute fft base for j-th fft word
//
// The FFT base for the j-th FFT word doing a b^n+c weighted transform is
//	ceil (j*n/FFTLEN)
// This routine returns ceil (j*n/FFTLEN) taking great care to return a
// value accurate to 53 bits.  This is important when j*n is really close to
// a multiple of FFTLEN (admittedly quite rare).  It would be really bad if
// rounding differences caused this routine to compute ceil (j*n/FFTLEN)
// differently than the weighting functions.
//

extern "C"
unsigned long gwfft_base (
	void	*dd_data_arg,
	unsigned long j)
{
	double	tempdbl, ceildbl, diffdbl;
	dd_real temp;
	unsigned long bpower;

// For speed, try this with plain old doubles first

	if (j == 0) return (0);
	tempdbl = (double) j * dd_data->gwdbl__fractional_b_per_word;
	ceildbl = ceil (tempdbl);
	diffdbl = ceildbl - tempdbl;
	if (diffdbl > 0.001 && diffdbl < 0.999)
		return ((unsigned long) ceildbl + j * dd_data->integral_b_per_word);

// We cannot be certain of the ceiling result, use doubledoubles to do the calculation

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = (int) gwceil (temp);
	END_x86_FIX
	return (bpower);
}

extern "C"
void gwsincos012by_weighted (
	void	*dd_data_arg,
	unsigned long x,
	unsigned long upper_x_incr,
	unsigned long N,
	unsigned long col,
	double	*results,
	int	incr)
{
	dd_real twopi_over_N, sine, cosine, sine2, cosine2;
	dd_real temp, bpower, weight, inv_weight;
	int	i;

	x86_FIX
	temp = (double) col * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) col;
	weight = exp (dd_data->gw__logb * bpower);
	inv_weight = dd_data->gw__over_fftlen / weight;

	// This hack causes SSE2 code to get duplicated weight and inv_weight values, while
	// AVX code does not (it uses the vbroadcastsd instruction to reduce memory usage)
	if (incr >= 4) {
		results[24] = weight;
		results[25] = inv_weight;
	} else {
		for (i = 0; i < incr; i++) {
			results[i] = weight;
			results[incr+i] = inv_weight;
		}
		results += 2*incr;
	}

	// Now output the weighted sin/cos data
	twopi_over_N = dd_real::_2pi / (double) N;
	for (i = 0; i < incr; i++) {
		sincos (twopi_over_N * (double) (x + i * upper_x_incr), sine, cosine);

		sine2 = sine * cosine * 2.0;
		cosine2 = sqr (cosine) - sqr (sine);

		sine = sine + epsilon;		// Hack to avoid divide-by-zero errors
		sine2 = sine2 + epsilon;

		// This hack outputs the weights in a different order for SSE2 and AVX.
		// We changed the AVX output order for better grouping of data in the cache lines.
		if (incr < 4) {
			results[i] = sine * weight;
			results[incr+i] = cosine / sine;
			results[2*incr+i] = sine * inv_weight;
			results[3*incr+i] = sine2 * weight;
			results[4*incr+i] = cosine2 / sine2;
			results[5*incr+i] = sine2 * inv_weight;
		} else {
			results[i] = cosine / sine;
			results[incr+i] = cosine2 / sine2;
			results[2*incr+i] = sine * weight;
			results[3*incr+i] = sine2 * weight;
			results[4*incr+i] = sine * inv_weight;
			results[5*incr+i] = sine2 * inv_weight;
		}
	}
	END_x86_FIX
}

extern "C"
void gwsincos15by_weighted (
	void	*dd_data_arg,
	unsigned long x,
	unsigned long upper_x_incr,
	unsigned long N,
	unsigned long col,
	double	*results,
	int	incr)
{
	dd_real twopi_over_N, sine, cosine, sine2, cosine2, sine4, cosine4, sine5, cosine5;
	dd_real temp, bpower, weight, inv_weight;
	int	i;

	x86_FIX
	temp = (double) col * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) col;
	weight = exp (dd_data->gw__logb * bpower);
	inv_weight = dd_data->gw__over_fftlen / weight;

	twopi_over_N = dd_real::_2pi / (double) N;
	for (i = 0; i < incr; i++) {
		sincos (twopi_over_N * (double) (x + i * upper_x_incr), sine, cosine);

		sine2 = sine * cosine * 2.0;
		cosine2 = sqr (cosine) - sqr (sine);

		sine4 = sine2 * cosine2 * 2.0;
		cosine4 = sqr (cosine2) - sqr (sine2);

		sine5 = sine * cosine4 + sine4 * cosine;
		cosine5 = cosine * cosine4 - sine * sine4;

		sine = sine + epsilon;		// Hack to avoid divide-by-zero errors
		sine5 = sine5 + epsilon;

		// This hack outputs the weights in a different order for SSE2 and AVX.
		// We changed the AVX output order for better grouping of data in the cache lines.
		if (incr < 4) {
			results[i] = sine * weight;
			results[incr+i] = cosine / sine;
			results[2*incr+i] = sine * inv_weight;
			results[3*incr+i] = sine5 * weight;
			results[4*incr+i] = cosine5 / sine5;
			results[5*incr+i] = sine5 * inv_weight;
		} else {
			results[i] = cosine / sine;
			results[incr+i] = cosine5 / sine5;
			results[2*incr+i] = sine * weight;
			results[3*incr+i] = sine5 * weight;
			results[4*incr+i] = sine * inv_weight;
			results[5*incr+i] = sine5 * inv_weight;
		}
	}
	END_x86_FIX
}

// Compute several column weights to later be applied to sine values of all-complex pre-multipliers.  One pass AVX-512 FFTs.

extern "C"
void gwfft_colweights (
	void	*dd_data_arg,
	void	*colweights_arg,	/* Area to output column weights */
	int	count)			/* Number of column weights to generate */
{
	unsigned long j;
	dd_real	*colweights;
	dd_real temp, bpower;

	x86_FIX
	colweights = (dd_real *) colweights_arg;
	for (j = 0; j < (unsigned long) count; j++) {
		temp = (double) j * dd_data->gw__num_b_per_word;
		bpower = gwceil (temp) - temp;
		if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
		*colweights++ = exp (dd_data->gw__logb * bpower);
	}
	END_x86_FIX
}

// Apply the column weights to the all-complex premultipliers.  One pass AVX-512 FFTs.

extern "C"
void gwsincos1plusby_colweighted (
	unsigned long x,
	unsigned long inc,
	unsigned long N,
	void	*colweights_arg,	/* Column weights generated by gwfft_colweights */
	double	*results,
	int	incr,
	int	amt)
{
	dd_real	*colweights;
	dd_real twopi_over_N, sineinc, cosineinc, sine0, cosine0;
	int	num_output;

	x86_FIX
	twopi_over_N = dd_real::_2pi / (double) N;
	sincos (twopi_over_N * (double) inc, sineinc, cosineinc);
	sincos (twopi_over_N * (double) x, sine0, cosine0);

	colweights = (dd_real *) colweights_arg;
	for (num_output = 0; ; colweights++) {
		dd_real last_sine;
		/* Output the sine * column weight and cosine/sine value while protecting against divide by zero */
		results[0] = sine0;
		results[0] += epsilon;
		results[incr] = cosine0 / results[0];
		results[0] = sine0 * *colweights;
		results[0] += epsilon;
		results += incr + incr;
		num_output++;
		/* Break if we've output the requested amount */
		if (num_output == amt) break;
		/* Compute next cosine and sine value */
		last_sine = sine0;
		sine0 = sine0 * cosineinc + cosine0 * sineinc;
		cosine0 = cosine0 * cosineinc - last_sine * sineinc;
	}
	END_x86_FIX
}

// This computes the two FFT weights and the two fudged weights in one call.

extern "C"
void gwfft_weights_fudged (
	void	*dd_data_arg,
	unsigned long j,
	unsigned long b,
	double	*fft_weight,
	double	*fft_weight_inverse,
	double	*fft_weight_over_b,
	double	*fft_weight_inverse_times_b)
{
	dd_real temp, bpower, weight, weight_inverse;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	weight = exp (dd_data->gw__logb * bpower);
	*fft_weight = double (weight);
	*fft_weight_over_b = double (weight / (double) b);
	weight_inverse = 1.0 / weight;
	*fft_weight_inverse = double (weight_inverse);
	*fft_weight_inverse_times_b = double (weight_inverse * (double) b);
	END_x86_FIX
}

/* Special version that merges group multiplier with the sine value from gwsincos1plus0123456789ABby. */
/* Used in AVX-512 all-omplex FFTs. */

extern "C"
void gwfft_weights_times_sine (
	void	*dd_data_arg,
	unsigned long j,		/* For computing group multiplier */
	unsigned long x,		/* For computing roots-of-minus-one sine value */
	unsigned long N,		/* For computing sine value */
	double	*fft_weight,		/* Returned group multiplier * sine */
	double	*fft_weight_inverse)	/* Returned 1/group multiplier * sine */
{
	dd_real temp, bpower, weight, sine0;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;
	weight = exp (dd_data->gw__logb * bpower);

	if (x == 0) {
		sine0 = epsilon;	/* Compensate for the way we calculated cos/sin */
	} else {
		sine0 = sin (dd_real::_2pi * (double) x / (double) N);
	}

	if (fft_weight != NULL) *fft_weight = double (weight * sine0);
	if (fft_weight_inverse != NULL) *fft_weight_inverse = double (sine0 / weight);

	END_x86_FIX
}

// Like the gwfft_weight routines but used for r4dwpn FFTs where part of the weight is applied
// during the FFT process.

extern "C"
double gwfft_partial_weight (
	void	*dd_data_arg,
	unsigned long j,
	unsigned long col)
{
	dd_real temp, j_bpower, col_bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	j_bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) j_bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;

	temp = (double) col * dd_data->gw__num_b_per_word;
	col_bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) col_bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) col;

	result = exp (dd_data->gw__logb * (j_bpower - col_bpower));
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.

extern "C"
double gwfft_partial_weight_sloppy (
	void	*dd_data_arg,
	unsigned long j,
	unsigned long col)
{
	dd_real temp, j_bpower, col_bpower;
	double	bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	j_bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) j_bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;

	temp = (double) col * dd_data->gw__num_b_per_word;
	col_bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) col_bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) col;
	END_x86_FIX

	bpower = double (j_bpower - col_bpower);
	if (bpower == dd_data->last_partial_sloppy_power[0]) return (dd_data->last_partial_sloppy_result[0]);
	if (bpower == dd_data->last_partial_sloppy_power[1]) return (dd_data->last_partial_sloppy_result[1]);
	result = pow (dd_data->gwdbl__b, bpower);
	dd_data->last_partial_sloppy_power[0] = dd_data->last_partial_sloppy_power[1];
	dd_data->last_partial_sloppy_result[0] = dd_data->last_partial_sloppy_result[1];
	dd_data->last_partial_sloppy_power[1] = bpower;
	dd_data->last_partial_sloppy_result[1] = result;
	return (result);
}

// Compute the inverse of the fft partial weight

extern "C"
double gwfft_partial_weight_inverse (
	void	*dd_data_arg,
	unsigned long j,
	unsigned long col)
{
	dd_real temp, j_bpower, col_bpower, result;

	x86_FIX
	temp = (double) j * dd_data->gw__num_b_per_word;
	j_bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) j_bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) j;

	temp = (double) col * dd_data->gw__num_b_per_word;
	col_bpower = gwceil (temp) - temp;
	if (! dd_data->gw__c_is_one) col_bpower += dd_data->gw__logb_abs_c_div_fftlen * (double) col;

	result = exp (dd_data->gw__logb * (col_bpower - j_bpower));
	END_x86_FIX
	return (double (result));
}

// Like the above, but faster and does not guarantee quite as much accuracy.
// We can be very sloppy as these weights are used to read FFT data values
// when writing save files.

extern "C"
double gwfft_partial_weight_inverse_sloppy (
	void	*dd_data_arg,
	unsigned long j,
	unsigned long col)
{
	double	temp, j_bpower, col_bpower, bpower, result;

	temp = (double) j * dd_data->gwdbl__fractional_b_per_word;
	j_bpower = ceil (temp) - temp;
	if (j_bpower < 0.001 || j_bpower > 0.999) {
		dd_real	dd_temp, dd_bpower;
		x86_FIX
		dd_temp = (double) j * dd_data->gw__num_b_per_word;
		dd_bpower = gwceil (dd_temp) - dd_temp;
		j_bpower = double (dd_bpower);
		END_x86_FIX
	}
	if (! dd_data->gw__c_is_one) j_bpower += dd_data->gwdbl__logb_abs_c_div_fftlen * (double) j;

	temp = (double) col * dd_data->gwdbl__fractional_b_per_word;
	col_bpower = ceil (temp) - temp;
	if (col_bpower < 0.001 || col_bpower > 0.999) {
		dd_real	dd_temp, dd_col_bpower;
		x86_FIX
		dd_temp = (double) col * dd_data->gw__num_b_per_word;
		dd_col_bpower = gwceil (dd_temp) - dd_temp;
		col_bpower = double (dd_col_bpower);
		END_x86_FIX
	}
	if (! dd_data->gw__c_is_one) col_bpower += dd_data->gwdbl__logb_abs_c_div_fftlen * (double) col;

	bpower = col_bpower - j_bpower;
	if (bpower == dd_data->last_partial_inv_sloppy_power[0]) return (dd_data->last_partial_inv_sloppy_result[0]);
	if (bpower == dd_data->last_partial_inv_sloppy_power[1]) return (dd_data->last_partial_inv_sloppy_result[1]);
	result = pow (dd_data->gwdbl__b, bpower);
	dd_data->last_partial_inv_sloppy_power[0] = dd_data->last_partial_inv_sloppy_power[1];
	dd_data->last_partial_inv_sloppy_result[0] = dd_data->last_partial_inv_sloppy_result[1];
	dd_data->last_partial_inv_sloppy_power[1] = bpower;
	dd_data->last_partial_inv_sloppy_result[1] = result;
	return (result);
}
