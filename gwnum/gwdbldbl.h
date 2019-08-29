/*----------------------------------------------------------------------
| gwdbldbl.h
|
| This file contains the headers for the gwnum helper routines that use
| extended-precision floats.
| 
|  Copyright 2005-2018 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWDBLDBL_H
#define _GWDBLDBL_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/* Include common definitions */

#include "gwcommon.h"

/* Extended precision helper routines */

void gwasm_constants (double *);
void gwsincos (unsigned long, unsigned long, double *);

#define gwsincos1by1(a,b,c)	gwsincos1by(a,b,c,1)
#define gwsincos1by2(a,b,c)	gwsincos1by(a,b,c,2)
#define gwsincos1by4(a,b,c)	gwsincos1by(a,b,c,4)
#define gwsincos1by8(a,b,c)	gwsincos1by(a,b,c,8)
void gwsincos1by (unsigned long, unsigned long, double *, int);

#define gwsincos12by1(a,b,c)		gwsincos12345by(a,b,c,1,2)
#define gwsincos12by2(a,b,c)		gwsincos12345by(a,b,c,2,2)
#define gwsincos12by4(a,b,c)		gwsincos12345by(a,b,c,4,2)
#define gwsincos12by8(a,b,c)		gwsincos12345by(a,b,c,8,2)
#define gwsincos123by1(a,b,c)		gwsincos12345by(a,b,c,1,3)
#define gwsincos123by2(a,b,c)		gwsincos12345by(a,b,c,2,3)
#define gwsincos123by4(a,b,c)		gwsincos12345by(a,b,c,4,3)
#define gwsincos123by8(a,b,c)		gwsincos12345by(a,b,c,8,3)
#define gwsincos1234by1(a,b,c)		gwsincos12345by(a,b,c,1,4)
#define gwsincos1234by2(a,b,c)		gwsincos12345by(a,b,c,2,4)
#define gwsincos1234by4(a,b,c)		gwsincos12345by(a,b,c,4,4)
#define gwsincos1234by8(a,b,c)		gwsincos12345by(a,b,c,8,4)
#define gwsincos12345by8(a,b,c)		gwsincos12345by(a,b,c,8,5)
#define gwsincos123456by8(a,b,c)	gwsincos12345by(a,b,c,8,6)
#define gwsincos1234567by2(a,b,c)	gwsincos12345by(a,b,c,2,7)
#define gwsincos1234567by4(a,b,c)	gwsincos12345by(a,b,c,4,7)
#define gwsincos1234567by8(a,b,c)	gwsincos12345by(a,b,c,8,7)
#define gwsincos12345678by8(a,b,c)	gwsincos12345by(a,b,c,8,8)
void gwsincos12345by (unsigned long, unsigned long, double *, int, int);

#define gwsincos13by1(a,b,c)		gwsincos13579by(a,b,c,1,2)
#define gwsincos13by2(a,b,c)		gwsincos13579by(a,b,c,2,2)
#define gwsincos13by4(a,b,c)		gwsincos13579by(a,b,c,4,2)
#define gwsincos13by8(a,b,c)		gwsincos13579by(a,b,c,8,2)
#define gwsincos135by1(a,b,c)		gwsincos13579by(a,b,c,1,3)
#define gwsincos135by8(a,b,c)		gwsincos13579by(a,b,c,8,3)
#define gwsincos1357by1(a,b,c)		gwsincos13579by(a,b,c,1,4)
#define gwsincos1357by8(a,b,c)		gwsincos13579by(a,b,c,8,4)
#define gwsincos13579by8(a,b,c)		gwsincos13579by(a,b,c,8,5)
#define gwsincos13579Bby8(a,b,c)	gwsincos13579by(a,b,c,8,6)
#define gwsincos13579BDFby8(a,b,c)	gwsincos13579by(a,b,c,8,8)
void gwsincos13579by (unsigned long, unsigned long, double *, int, int);

#define gwsincos15by1(a,b,c)		gwsincos159by(a,b,c,1,2)
#define gwsincos15by2(a,b,c)		gwsincos159by(a,b,c,2,2)
#define gwsincos15by4(a,b,c)		gwsincos159by(a,b,c,4,2)
#define gwsincos159Dby2(a,b,c)		gwsincos159by(a,b,c,2,4)
#define gwsincos159Dby4(a,b,c)		gwsincos159by(a,b,c,4,4)
#define gwsincos159Dby8(a,b,c)		gwsincos159by(a,b,c,8,4)
void gwsincos159by (unsigned long, unsigned long, double *, int, int);

#define gwsincos125by2(a,b,c)		gwsincos125by(a,b,c,2)
#define gwsincos125by4(a,b,c)		gwsincos125by(a,b,c,4)
#define gwsincos125by8(a,b,c)		gwsincos125by(a,b,c,8)
void gwsincos125by (unsigned long, unsigned long, double *, int);

#define gwsincos1by1_raw(a,b,c)		gwsincos1234by_raw(a,b,c,1,1)
#define gwsincos1by2_raw(a,b,c)		gwsincos1234by_raw(a,b,c,2,1)
#define gwsincos1by4_raw(a,b,c)		gwsincos1234by_raw(a,b,c,4,1)
#define gwsincos1by8_raw(a,b,c)		gwsincos1234by_raw(a,b,c,8,1)
#define gwsincos12by1_raw(a,b,c)	gwsincos1234by_raw(a,b,c,1,2)
#define gwsincos12by2_raw(a,b,c)	gwsincos1234by_raw(a,b,c,2,2)
#define gwsincos12by4_raw(a,b,c)	gwsincos1234by_raw(a,b,c,4,2)
#define gwsincos1234by2_raw(a,b,c)	gwsincos1234by_raw(a,b,c,2,4)
#define gwsincos1234by4_raw(a,b,c)	gwsincos1234by_raw(a,b,c,4,4)
#define gwsincos1234by8_raw(a,b,c)	gwsincos1234by_raw(a,b,c,8,4)
void gwsincos1234by_raw (unsigned long, unsigned long, double *, int, int);

#define gwsincos1plus0123by1(a,b,c,d)		gwsincos1plusby(a,b,c,d,1,4)
#define gwsincos1plus0123by2(a,b,c,d)		gwsincos1plusby(a,b,c,d,2,4)
#define gwsincos1plus0123by4(a,b,c,d)		gwsincos1plusby(a,b,c,d,4,4)
#define gwsincos1plus0123by8(a,b,c,d)		gwsincos1plusby(a,b,c,d,8,4)
#define gwsincos1plus01234by1(a,b,c,d)		gwsincos1plusby(a,b,c,d,1,5)
#define gwsincos1plus01234by2(a,b,c,d)		gwsincos1plusby(a,b,c,d,2,5)
#define gwsincos1plus01234by4(a,b,c,d)		gwsincos1plusby(a,b,c,d,4,5)
#define gwsincos1plus01234by8(a,b,c,d)		gwsincos1plusby(a,b,c,d,8,5)
#define gwsincos1plus01234567by2(a,b,c,d)	gwsincos1plusby(a,b,c,d,2,8)
#define gwsincos1plus01234567by4(a,b,c,d)	gwsincos1plusby(a,b,c,d,4,8)
#define gwsincos1plus01234567by8(a,b,c,d)	gwsincos1plusby(a,b,c,d,8,8)
#define gwsincos1plus0123456789ABby8(a,b,c,d)	gwsincos1plusby(a,b,c,d,8,12)
void gwsincos1plusby (unsigned long, unsigned long, unsigned long, double *, int, int);

#define gwcos1plus01234567by8(a,b,c,d)		gwcos1plusby(a,b,c,d,8,8)
#define gwcos1plus0123456789ABby8(a,b,c,d)	gwcos1plusby(a,b,c,d,8,12)
void gwcos1plusby (unsigned long, unsigned long, unsigned long, double *, int, int);

#define gwsincos012by2_weighted(a,b,c,d,e,f)	gwsincos012by_weighted(a,b,c,d,e,f,2)
#define gwsincos012by4_weighted(a,b,c,d,e,f)	gwsincos012by_weighted(a,b,c,d,e,f,4)
void gwsincos012by_weighted (void *, unsigned long, unsigned long, unsigned long, unsigned long, double *, int);

#define gwsincos15by2_weighted(a,b,c,d,e,f)	gwsincos15by_weighted(a,b,c,d,e,f,2)
#define gwsincos15by4_weighted(a,b,c,d,e,f)	gwsincos15by_weighted(a,b,c,d,e,f,4)
void gwsincos15by_weighted (void *, unsigned long, unsigned long, unsigned long, unsigned long, double *, int);

#define gwsincos1plus01234567by8_colweighted(a,b,c,d,e)	gwsincos1plusby_colweighted(a,b,c,d,e,8,8)
void gwsincos1plusby_colweighted (unsigned long, unsigned long, unsigned long, void *, double *, int, int);

#define gwsincos1234by8_sqrthalf(a,b,c)		gwsincos1234by_sqrthalf(a,b,c,8)
void gwsincos1234by_sqrthalf (unsigned long, unsigned long, double *, int);

#define gwsincos123by1_special7(a,b,c)		gwsincos123by_special7(a,b,c,1)
#define gwsincos123by8_special7(a,b,c)		gwsincos123by_special7(a,b,c,8)
void gwsincos123by_special7 (unsigned long, unsigned long, double *, int);

#define gwsincos135by1_special7(a,b,c)		gwsincos135by_special7(a,b,c,1)
#define gwsincos135by8_special7(a,b,c)		gwsincos135by_special7(a,b,c,8)
void gwsincos135by_special7 (unsigned long, unsigned long, double *, int);

void *gwdbldbl_data_alloc (void);
void gwfft_weight_setup (void *, int, double, unsigned long, unsigned long, signed long, unsigned long);
double gwfft_weight (void *, unsigned long);
double gwfft_weight_over_sqrt_fftlen (void *, unsigned long);
double gwfft_weight_sloppy (void *, unsigned long);
double gwfft_weight_inverse (void *, unsigned long);
double gwfft_weight_inverse_squared (void *, unsigned long);
double gwfft_weight_inverse_sloppy (void *, unsigned long);
double gwfft_weight_inverse_over_fftlen (void *, unsigned long);
void gwfft_weights3 (void *, unsigned long, double *, double *, double *);
double gwfft_weight_exponent (void *, unsigned long);
double gwfft_weight_no_c (void *, unsigned long);
unsigned long gwfft_base (void *, unsigned long);
void gwfft_weights_fudged (void	*, unsigned long, unsigned long, double *, double *, double *, double *);
void gwfft_weights_times_sine (void *, unsigned long, unsigned long, unsigned long, double *, double *);
double gwfft_partial_weight (void *, unsigned long, unsigned long);
double gwfft_partial_weight_sloppy (void *, unsigned long, unsigned long);
double gwfft_partial_weight_inverse (void *, unsigned long, unsigned long);
double gwfft_partial_weight_inverse_sloppy (void *, unsigned long, unsigned long);
void gwfft_colweights (void *, void *, int);

#ifdef __cplusplus
}
#endif

#endif
