/**************************************************************
 *
 *	ecm.c
 *
 *	ECM and P-1 factoring program
 *
 *	Original author:  Richard Crandall - www.perfsci.com
 *	Adapted to Mersenne numbers and optimized by George Woltman
 *	Further optimizations from Paul Zimmerman's GMP-ECM program
 *	Other important ideas courtesy of Peter Montgomery.
 *
 *	c. 1997 Perfectly Scientific, Inc.
 *	c. 1998-2020 Mersenne Research, Inc.
 *	All Rights Reserved.
 *
 *************************************************************/

/* IDEAS:
   The P-1 stage 2 bit array could also use a mod 30 scheme to cut the
	number of bits in the array nearly in half)
   We could allocate a bit array half the size and do the pairings as we fill
	the bit array.
   Have 16 PRAC values, and/or 16 PRAC values and +/-8 offset.  Then keep one
	   precomputed 4-bit or 8-bit value per B1 prime that points to
	   a precomputed optimal PRAC initial value.  There are 5.7 million
	   primes < 100M -- or maybe the cost of testing 256 different PRAC
	   combos isn't great.  Another possibility would be to remember only
	   those primes where a better solution is available
*/

/* Global variables */

int	QA_IN_PROGRESS = FALSE;
int	QA_TYPE = 0;
giant	QA_FACTOR = NULL;

/********************/
/* Utility routines */
/********************/

/* Test if N is a probable prime */
/* Compute i^(N-1) mod N for i = 3,5,7 */

int isProbablePrime (
	gwhandle *gwdata,
	giant	N)
{
	int	i, j, len, retval;
	gwnum	t1, t2;
	giant	x;

	if (isone (N)) return (TRUE);

	retval = TRUE;		/* Assume it is a probable prime */
	t1 = gwalloc (gwdata);
	len = bitlen (N);
	for (i = 3; retval && i <= 7; i += 2) {
		t2 = gwalloc (gwdata);
		dbltogw (gwdata, (double) 1.0, t1);
		dbltogw (gwdata, (double) i, t2);
		gwfft (gwdata, t2, t2);
		for (j = 1; j <= len; j++) {
			gwsquare (gwdata, t1);
			if (bitval (N, len-j)) gwfftmul (gwdata, t2, t1);
		}
		gwfree (gwdata, t2);
		x = popg (&gwdata->gdata, ((int) gwdata->bit_length >> 5) + 10);
		if (gwtogiant (gwdata, t1, x)) retval = FALSE;	/* Technically, prime status is unknown on an unexpected error */
		else {
			modgi (&gwdata->gdata, N, x);
			iaddg (-i, x);
			if (!isZero (x)) retval = FALSE;	/* Not a prime */
		}
		pushg (&gwdata->gdata, 1);
	}
	gwfree (gwdata, t1);
	return (retval);
}

/*************************************************/
/* ECM structures and setup/termination routines */
/*************************************************/

/* Data maintained during ECM process */

#define POOL_3MULT	2	/* Modinv algorithm that takes 3 multiplies */
#define POOL_N_SQUARED	4	/* Use O(N^2) multiplies modinv algorithm */

typedef struct {
	gwhandle gwdata;	/* GWNUM handle */
	int	thread_num;	/* Worker thread number */
	unsigned long D;	/* Stage 2 loop size */
	unsigned long E;	/* Suyama's power in stage 2 */
	gwnum	*nQx;		/* Array of data used in stage 2 */
	char	*pairings;	/* Bits used in determining if primes pair */
	gwnum	Qprevmx, Qprevmz, Qmx, Qmz, Q2Dxplus1, Q2Dxminus1, *mQx;
	unsigned int mQx_count;
	int	TWO_FFT_STAGE2;	/* Type of ECM stage 2 to execute */
	int	pool_type;	/* Modinv algorithm type to use */
	unsigned int pool_count;/* Count of pooled normalizes */
	int	pool_ffted;	/* TRUE if pooled values were pre-ffted */
	gwnum	pool_modinv_value;/* Value we will eventually do a modinv on */
	gwnum	*pool_values;	/* Array of values to normalize */
	gwnum	*poolz_values;	/* Array of z values we are normalize */
	unsigned long modinv_count; /* Stats - count of modinv calls */
	void	*sieve_info;
} ecmhandle;

/* Perform cleanup functions. */

void ecm_cleanup (
	ecmhandle *ecmdata)
{
	free (ecmdata->nQx);
	free (ecmdata->pool_values);
	free (ecmdata->poolz_values);
	free (ecmdata->mQx); 
	free (ecmdata->pairings);
	gwdone (&ecmdata->gwdata);
	end_sieve (ecmdata->sieve_info);	
	memset (ecmdata, 0, sizeof (ecmhandle));
}

/* Free some, but not all, memory in the ecmdata structure.  Keep */
/* memory that can be reused in next ECM curve. */

void ecm_partial_cleanup (
	ecmhandle *ecmdata)
{
	free (ecmdata->nQx); ecmdata->nQx = NULL;
	free (ecmdata->pool_values); ecmdata->pool_values = NULL;
	free (ecmdata->poolz_values); ecmdata->poolz_values = NULL;
	free (ecmdata->mQx); ecmdata->mQx = NULL;
	free (ecmdata->pairings); ecmdata->pairings = NULL;
	gwfreeall (&ecmdata->gwdata);
}

/**************************************************************
 *
 *	Functions
 *
 **************************************************************/

/* computes 2P=(x2:z2) from P=(x1:z1), uses the global variables Ad4 */

int ell_dbl (
	ecmhandle *ecmdata,
	const gwnum x1,
	const gwnum z1,
	gwnum	x2,
	gwnum	z2,
	const gwnum Ad4)
{					/* 10 FFTs */
	gwnum	t1, t3;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t3 = gwalloc (&ecmdata->gwdata);
	if (t3 == NULL) goto oom;
	gwaddsub4 (&ecmdata->gwdata, x1, z1, t1, x2);
	gwsquare (&ecmdata->gwdata, t1);	/* t1 = (x1 + z1)^2 */
	gwsquare (&ecmdata->gwdata, x2);	/* t2 = (x1 - z1)^2 (store in x2) */
	gwsub3 (&ecmdata->gwdata, t1, x2, t3);	/* t3 = t1 - t2 = 4 * x1 * z1 */
	gwfft (&ecmdata->gwdata, t3, t3);
	gwfft (&ecmdata->gwdata, x2, x2);
	gwfftadd3 (&ecmdata->gwdata, t3, x2, t1);/* Compute the fft of t1! */
	gwfftfftmul (&ecmdata->gwdata, Ad4, x2, x2);	/* x2 = t2 * Ad4 */
	gwfft (&ecmdata->gwdata, x2, x2);
	gwfftadd3 (&ecmdata->gwdata, x2, t3, z2);	/* z2 = (t2 * Ad4 + t3) */
	gwfftfftmul (&ecmdata->gwdata, t3, z2, z2);	/* z2 = z2 * t3 */
	gwfftfftmul (&ecmdata->gwdata, t1, x2, x2);	/* x2 = x2 * t1 */
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t3);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
   Assumes that Q-R=P or R-Q=P where P=(xdiff:zdiff). */

#ifdef ELL_ADD_USED
int ell_add (
	ecmhandle *ecmdata,
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum	xdiff,
	gwnum	zdiff,
	gwnum	x3,
	gwnum	z3)
{					/* 16 FFTs */
	gwnum	t1, t2, t3;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	t3 = gwalloc (&ecmdata->gwdata);
	if (t3 == NULL) goto oom;
	gwaddsub4 (&ecmdata->gwdata, x1, z1, t1, t2);
						/* t1 = (x1 + z1)(x2 - z2) */
						/* t2 = (x1 - z1)(x2 + z2) */
	gwsub3 (&ecmdata->gwdata, x2, z2, t3);
	gwmul (&ecmdata->gwdata, t3, t1);
	gwadd3 (&ecmdata->gwdata, x2, z2, t3);
	gwmul (&ecmdata->gwdata, t3, t2);
	gwaddsub (&ecmdata->gwdata, t2, t1);	/* x3 = (t2 + t1)^2 * zdiff */
	gwsquare (&ecmdata->gwdata, t2);
	gwmul (&ecmdata->gwdata, zdiff, t2);
	gwsquare (&ecmdata->gwdata, t1);	/* z3 = (t2 - t1)^2 * xdiff */
	gwmul (&ecmdata->gwdata, xdiff, t1);
	gwcopy (&ecmdata->gwdata, t2, x3);
	gwcopy (&ecmdata->gwdata, t1, z3);
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	gwfree (&ecmdata->gwdata, t3);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}
#endif

/* Like ell_add except that x1, z1, xdiff, and zdiff have been FFTed */
/* NOTE: x2 and z2 represent the FFTs of (x2+z2) and (x2-z2) respectively. */

#ifdef ELL_ADD_SPECIAL_USED
int ell_add_special (
	ecmhandle *ecmdata,
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum	xdiff,
	gwnum	zdiff,
	gwnum	x3,
	gwnum	z3)
{				/* 10 FFTs */
	gwnum	t1, t2;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	gwfftaddsub4 (&ecmdata->gwdata, x1, z1, t1, t2);
						/* t1 = (x1 + z1)(x2 - z2) */
						/* t2 = (x1 - z1)(x2 + z2) */
	gwfftfftmul (&ecmdata->gwdata, z2, t1, t1);
	gwfftfftmul (&ecmdata->gwdata, x2, t2, t2);
	gwaddsub (&ecmdata->gwdata, t2, t1);	/* x3 = (t2 + t1)^2 * zdiff */
	gwstartnextfft (&ecmdata->gwdata, TRUE);
	gwsquare (&ecmdata->gwdata, t2);
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftmul (&ecmdata->gwdata, zdiff, t2);
	gwstartnextfft (&ecmdata->gwdata, TRUE);
	gwsquare (&ecmdata->gwdata, t1);	/* z3 = (t2 - t1)^2 * xdiff */
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftmul (&ecmdata->gwdata, xdiff, t1);
	gwcopy (&ecmdata->gwdata, t2, x3);
	gwcopy (&ecmdata->gwdata, t1, z3);
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}
#endif

/* Like ell_add except that x1, z1, xdiff, and zdiff have been FFTed */
/* NOTE: x2 and z2 represent the FFTs of (x2+z2) and (x2-z2) respectively. */
/* ALSO: Like ell_add_special except x3 is returned in zdiff and z3 is returned in xdiff */
/* caller must then gwswap xdiff and zdiff.  This save two gwcopies. */

#ifdef ELL_ADD_SPECIAL2_USED
int ell_add_special2 (
	ecmhandle *ecmdata,
	const gwnum x1,
	const gwnum z1,
	const gwnum x2,
	const gwnum z2,
	gwnum	xdiff,
	gwnum	zdiff)
{				/* 10 FFTs */
	gwnum	t1, t2;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	gwfftaddsub4 (&ecmdata->gwdata, x1, z1, t1, t2);/* Calc (x1 + z1) and (z1 - z1) */
	gwfftfftmul (&ecmdata->gwdata, z2, t1, t1);	/* t1 = (x1 + z1)(x2 - z2) */
	gwfftfftmul (&ecmdata->gwdata, x2, t2, t2);	/* t2 = (x1 - z1)(x2 + z2) */
	gwaddsub (&ecmdata->gwdata, t2, t1);		/* Calc t2 + t1 and t2 - t1 */
	gwstartnextfft (&ecmdata->gwdata, TRUE);	/* x3 = (t2 + t1)^2 * zdiff */
	gwsquare (&ecmdata->gwdata, t2);
	gwfft (&ecmdata->gwdata, t2, t2);
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftfftmul (&ecmdata->gwdata, t2, zdiff, zdiff); /* Return x3 in zdiff */
	gwstartnextfft (&ecmdata->gwdata, TRUE);	/* z3 = (t2 - t1)^2 * xdiff */
	gwsquare (&ecmdata->gwdata, t1);
	gwfft (&ecmdata->gwdata, t1, t1);
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftfftmul (&ecmdata->gwdata, t1, xdiff, xdiff); /* Return z3 in xdiff */
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}
#endif

/* Like ell_add_special2 except that xdiff and zdiff will be freed and replaced with the new x3 and z3 values */
/* NOTE: x2 and z2 represent the FFTs of (x2+z2) and (x2-z2) respectively. */

int ell_add_special3 (
	ecmhandle *ecmdata,
	const gwnum x1,
	const gwnum z1,
	const gwnum x2,
	const gwnum z2,
	gwnum	*xdiff_and_returned_x3,
	gwnum	*zdiff_and_returned_z3)
{				/* 10 FFTs */
	gwnum	xdiff = *xdiff_and_returned_x3;
	gwnum	zdiff = *zdiff_and_returned_z3;
	gwnum	t1, t2;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	gwfftaddsub4 (&ecmdata->gwdata, x1, z1, t1, t2);/* Calc (x1 + z1) and (z1 - z1) */
	gwfftfftmul (&ecmdata->gwdata, z2, t1, t1);	/* t1 = (x1 + z1)(x2 - z2) */
	gwfftfftmul (&ecmdata->gwdata, x2, t2, t2);	/* t2 = (x1 - z1)(x2 + z2) */
	gwaddsub (&ecmdata->gwdata, t2, t1);		/* Calc t2 + t1 and t2 - t1 */
	gwstartnextfft (&ecmdata->gwdata, TRUE);	/* x3 = (t2 + t1)^2 * zdiff */
	gwsquare (&ecmdata->gwdata, t2);
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftmul (&ecmdata->gwdata, zdiff, t2);		/* Final x3 */
	gwstartnextfft (&ecmdata->gwdata, TRUE);	/* z3 = (t2 - t1)^2 * xdiff */
	gwsquare (&ecmdata->gwdata, t1);
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftmul (&ecmdata->gwdata, xdiff, t1);		/* Final z3 */
	gwfree (&ecmdata->gwdata, xdiff);
	gwfree (&ecmdata->gwdata, zdiff);
	*xdiff_and_returned_x3 = t2;
	*zdiff_and_returned_z3 = t1;
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* This routine is called prior to a series of many ell_add_fft and */
/* ell_dbl_fft calls.  The sequence ends by calling ell_add_fft_last. */
/* Note: We used to simply just FFT x1 and z1.  However, convolution error */
/* in computing (x1+z1)^2 and the like was too great.  Instead, we now */
/* save the FFTs of (x1+z1) and (x1-z1).  The multiplication by xdiff */
/* and zdiff is now more complicated, but convolution errors are reduced */
/* since only one argument of any multiply will involve a value that is */
/* the sum of two FFTs rather than computing a properly normalized sum */
/* and then taking the FFT. */

void ell_begin_fft (
	ecmhandle *ecmdata,
	const gwnum x1,
	const gwnum z1,
	gwnum	x2,
	gwnum	z2)
{
	gwaddsub4 (&ecmdata->gwdata, x1, z1, x2, z2);	/* x2 = x1 + z1, z2 = x1 - z1 */
	gwfft (&ecmdata->gwdata, x2, x2);
	gwfft (&ecmdata->gwdata, z2, z2);
}

/* Like ell_dbl, but the input arguments are FFTs of x1=x1+z1, z1=x1-z1 */
/* The output arguments are also FFTs of x2=x2+z2, z2=x2-z2 */

int ell_dbl_fft (
	ecmhandle *ecmdata,
	const gwnum x1,
	const gwnum z1,
	gwnum	x2,
	gwnum	z2,
	const gwnum Ad4)
{					/* 10 FFTs, 4 adds */
	gwnum	t1, t3;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t3 = gwalloc (&ecmdata->gwdata);
	if (t3 == NULL) goto oom;
	gwfftfftmul (&ecmdata->gwdata, x1, x1, t1);	/* t1 = (x1 + z1)^2 */
	gwfftfftmul (&ecmdata->gwdata, z1, z1, x2);	/* t2 = (x1 - z1)^2 (store in x2) */
	gwsub3 (&ecmdata->gwdata, t1, x2, t3);		/* t3 = t1 - t2 = 4 * x1 * z1 */
	gwfft (&ecmdata->gwdata, t3, t3);
	gwfft (&ecmdata->gwdata, x2, x2);
	gwfftadd3 (&ecmdata->gwdata, t3, x2, t1);	/* Compute fft of t1! */
	gwstartnextfft (&ecmdata->gwdata, TRUE);
	gwfftfftmul (&ecmdata->gwdata, Ad4, x2, x2);	/* x2 = t2 * Ad4 */
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfft (&ecmdata->gwdata, x2, x2);
	gwfftadd3 (&ecmdata->gwdata, x2, t3, z2);	/* z2 = (t2 * Ad4 + t3) * t3 */
	gwfftfftmul (&ecmdata->gwdata, t3, z2, z2);
	gwfftfftmul (&ecmdata->gwdata, t1, x2, x2);	/* x2 = x2 * t1 */
	gwaddsub (&ecmdata->gwdata, x2, z2);		/* x2 = x2 + z2, z2 = x2 - z2 */
	gwfft (&ecmdata->gwdata, x2, x2);
	gwfft (&ecmdata->gwdata, z2, z2);
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t3);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Like ell_add but input arguments are FFTs of x1=x1+z1, z1=x1-z1, */
/* x2=x2+z2, z2=x2-z2, xdiff=xdiff+zdiff, zdiff=xdiff-zdiff. */
/* The output arguments are also FFTs of x3=x3+z3, z3=x3-z3 */

int ell_add_fft (
	ecmhandle *ecmdata,
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum 	xdiff,
	gwnum 	zdiff,
	gwnum 	x3,
	gwnum 	z3)
{				/* 12 FFTs, 6 adds */
	gwnum	t1, t2;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	gwfftfftmul (&ecmdata->gwdata, x1, z2, t1);
					/* t1 = (x1 + z1)(x2 - z2) */
	gwfftfftmul (&ecmdata->gwdata, x2, z1, t2);
					/* t2 = (x1 - z1)(x2 + z2) */
	gwaddsub (&ecmdata->gwdata, t2, t1);
	gwstartnextfft (&ecmdata->gwdata, TRUE);
	gwsquare (&ecmdata->gwdata, t2);
					/* t2 = (t2 + t1)^2 (will become x3) */
	gwsquare (&ecmdata->gwdata, t1);
					/* t1 = (t2 - t1)^2 (will become z3) */
	gwstartnextfft (&ecmdata->gwdata, FALSE);
	gwfftaddsub4 (&ecmdata->gwdata, xdiff, zdiff, x3, z3);
					/* x3 = xdiff = (xdiff + zdiff) */
					/* z3 = zdiff = (xdiff - zdiff) */
	gwfftmul (&ecmdata->gwdata, z3, t2);
					/* t2 = t2 * zdiff (new x3) */
	gwfftmul (&ecmdata->gwdata, x3, t1);
					/* t1 = t1 * xdiff (new z3) */
	gwaddsub (&ecmdata->gwdata, t2, t1);
					/* t2 = x3 + z3, t1 = x3 - z3 */
	gwfft (&ecmdata->gwdata, t2, x3);
	gwfft (&ecmdata->gwdata, t1, z3);
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Like ell_add_fft but output arguments are not FFTed. */

int ell_add_fft_last (
	ecmhandle *ecmdata,
	gwnum 	x1,
	gwnum 	z1,
	gwnum 	x2,
	gwnum 	z2,
	gwnum 	xdiff,
	gwnum 	zdiff,
	gwnum 	x3,
	gwnum 	z3)
{				/* 10 FFTs, 6 adds */
	gwnum	t1, t2;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	gwfftfftmul (&ecmdata->gwdata, x1, z2, t1);
					/* t1 = (x1 + z1)(x2 - z2) */
	gwfftfftmul (&ecmdata->gwdata, x2, z1, t2);
					/* t2 = (x1 - z1)(x2 + z2) */
	if (xdiff != x3) {
		gwaddsub4 (&ecmdata->gwdata, t2, t1, x3, z3);
		gwstartnextfft (&ecmdata->gwdata, TRUE);
		gwsquare (&ecmdata->gwdata, x3);
					/* x3 = (t2 + t1)^2 */
		gwsquare (&ecmdata->gwdata, z3);
					/* z3 = (t2 - t1)^2 */
		gwstartnextfft (&ecmdata->gwdata, FALSE);
		gwfftaddsub4 (&ecmdata->gwdata, xdiff, zdiff, t1, t2);
					/* t1 = xdiff = (xdiff + zdiff) */
					/* t2 = zdiff = (xdiff - zdiff) */
		gwfftmul (&ecmdata->gwdata, t2, x3);
					/* x3 = x3 * zdiff */
		gwfftmul (&ecmdata->gwdata, t1, z3);
					/* z3 = z3 * xdiff */
	} else {
		gwaddsub (&ecmdata->gwdata, t2, t1);
		gwstartnextfft (&ecmdata->gwdata, TRUE);
		gwsquare (&ecmdata->gwdata, t2);
		gwfft (&ecmdata->gwdata, t2, t2);
		gwsquare (&ecmdata->gwdata, t1);
		gwfft (&ecmdata->gwdata, t1, t1);
		gwstartnextfft (&ecmdata->gwdata, FALSE);
		gwfftaddsub4 (&ecmdata->gwdata, xdiff, zdiff, z3, x3);
		gwfftfftmul (&ecmdata->gwdata, t2, x3, x3);
		gwfftfftmul (&ecmdata->gwdata, t1, z3, z3);
	}
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Perform an elliptic multiply using an algorithm developed by */
/* Peter Montgomery.  Basically, we try to find a near optimal */
/* Lucas chain of additions that generates the number we are */
/* multiplying by.  This minimizes the number of calls to ell_dbl */
/* and ell_add. */

/* The costing function assigns an ell_dbl call a cost of 10 and */
/* an ell_add call a cost of 12.  This cost estimates the number */
/* of forward and inverse transforms performed. */

#define swap(a,b)	{t=a;a=b;b=t;}

unsigned long lucas_cost (
	uint64_t n,
	double	inv_v)
{
	uint64_t d, e, t, dmod3, emod3;
	unsigned long c;

	c = 0;
	while (n != 1) {
	    d = (uint64_t) ((double) n * inv_v + 0.5); e = n - d;
	    d = d - e;

	    c += 12;

	    while (d != e) {
		if (d < e) {
			swap (d,e);
		}
		if (d <= e + (e >> 2)) {
			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
				t = d;
				d = (d+d-e)/3;
				e = (e+e-t)/3;
				c += 36;
				continue;
			}
			if (dmod3 == emod3 && (d&1) == (e&1)) {
				d = (d-e) >> 1;
				c += 22;
				continue;
			}
		}
		if (d <= (e << 2)) {
			d = d-e;
			c += 12;
		} else if ((d&1) == (e&1)) {
			d = (d-e) >> 1;
			c += 22;
		} else if ((d&1) == 0) {
			d = d >> 1;
			c += 22;
		} else if ((dmod3 = d%3) == 0) {
			d = d/3-e;
			c += 46;
		} else if (dmod3 == 3 - (emod3 = e%3)) {
			d = (d-e-e)/3;
			c += 46;
		} else if (dmod3 == emod3) {
			d = (d-e)/3;
			c += 46;
		} else {
			e = e >> 1;
			c += 22;
		}
	    }
	    c += 10;
	    n = d;
	}

	return (c);
}

int lucas_mul (
	ecmhandle *ecmdata,
	gwnum	xx,
	gwnum	zz,
	uint64_t n,
	double	inv_v,
	gwnum	Ad4)
{
	uint64_t d, e, t, dmod3, emod3;
	gwnum	xA, zA, xB, zB, xC, zC, xs, zs, xt, zt;
	int	stop_reason;

	xA = gwalloc (&ecmdata->gwdata);
	if (xA == NULL) goto oom;
	zA = gwalloc (&ecmdata->gwdata);
	if (zA == NULL) goto oom;
	xB = gwalloc (&ecmdata->gwdata);
	if (xB == NULL) goto oom;
	zB = gwalloc (&ecmdata->gwdata);
	if (zB == NULL) goto oom;
	xC = gwalloc (&ecmdata->gwdata);
	if (xC == NULL) goto oom;
	zC = gwalloc (&ecmdata->gwdata);
	if (zC == NULL) goto oom;
	xs = xx;
	zs = zz;
	xt = gwalloc (&ecmdata->gwdata);
	if (xt == NULL) goto oom;
	zt = gwalloc (&ecmdata->gwdata);
	if (zt == NULL) goto oom;

	while (n != 1) {
	    ell_begin_fft (ecmdata, xx, zz, xA, zA);		/* A */
	    stop_reason = ell_dbl_fft (ecmdata, xA, zA, xB, zB, Ad4);		/* B = 2*A */
	    if (stop_reason) return (stop_reason);
	    gwcopy (&ecmdata->gwdata, xA, xC);
	    gwcopy (&ecmdata->gwdata, zA, zC);			/* C = A */

	    d = (uint64_t) ((double) n * inv_v + 0.5); e = n - d;
	    d = d - e;

	    while (d != e) {
		if (d < e) {
			swap (d, e);
			gwswap (xA, xB); gwswap (zA, zB);
		}
		if (d <= e + (e >> 2)) {
			if ((dmod3 = d%3) == 3 - (emod3 = e%3)) {
				stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xs, zs);/* S = A+B */
				if (stop_reason) return (stop_reason);
				stop_reason = ell_add_fft (ecmdata, xA, zA, xs, zs, xB, zB, xt, zt);/* T = A+S */
				if (stop_reason) return (stop_reason);
				stop_reason = ell_add_fft (ecmdata, xs, zs, xB, zB, xA, zA, xB, zB);/* B = B+S */
				if (stop_reason) return (stop_reason);
				gwswap (xt, xA); gwswap (zt, zA);/* A = T */
				t = d;
				d = (d+d-e)/3;
				e = (e+e-t)/3;
				continue;
			}
			if (dmod3 == emod3 && (d&1) == (e&1)) {
				stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xB, zB);/* B = A+B */
				if (stop_reason) return (stop_reason);
				stop_reason = ell_dbl_fft (ecmdata, xA, zA, xA, zA, Ad4);	/* A = 2*A */
				if (stop_reason) return (stop_reason);
				d = (d-e) >> 1;
				continue;
			}
		}
		if (d <= (e << 2)) {
			stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xC, zC);/* B = A+B */
			if (stop_reason) return (stop_reason);
			gwswap (xB, xC); gwswap (zB, zC);	/* C = B */
			d = d-e;
		} else if ((d&1) == (e&1)) {
			stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xB, zB);/* B = A+B */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_dbl_fft (ecmdata, xA, zA, xA, zA, Ad4);	/* A = 2*A */
			if (stop_reason) return (stop_reason);
			d = (d-e) >> 1;
		} else if ((d&1) == 0) {
			stop_reason = ell_add_fft (ecmdata, xA, zA, xC, zC, xB, zB, xC, zC);/* C = A+C */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_dbl_fft (ecmdata, xA, zA, xA, zA, Ad4);	/* A = 2*A */
			if (stop_reason) return (stop_reason);
			d = d >> 1;
		} else if ((dmod3 = d%3) == 0) {
			stop_reason = ell_dbl_fft (ecmdata, xA, zA, xs, zs, Ad4);	/* S = 2*A */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xt, zt);/* T = A+B */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xs, zs, xA, zA, xA, zA, xA, zA);/* A = S+A */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xs, zs, xt, zt, xC, zC, xC, zC);/* B = S+T */
			if (stop_reason) return (stop_reason);
			gwswap (xB, xC); gwswap (zB, zC);	/* C = B */
			d = d/3-e;
		} else if (dmod3 == 3 - (emod3 = e%3)) {
			stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xs, zs);/* S = A+B */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xA, zA, xs, zs, xB, zB, xB, zB);/* B = A+S */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_dbl_fft (ecmdata, xA, zA, xs, zs, Ad4);	/* S = 2*A */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xs, zs, xA, zA, xA, zA, xA, zA);/* A = S+A */
			if (stop_reason) return (stop_reason);
			d = (d-e-e)/3;
		} else if (dmod3 == emod3) {
			stop_reason = ell_add_fft (ecmdata, xA, zA, xB, zB, xC, zC, xt, zt);/* T = A+B */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xA, zA, xC, zC, xB, zB, xC, zC);/* C = A+C */
			if (stop_reason) return (stop_reason);
			gwswap (xt, xB); gwswap (zt, zB);	/* B = T */
			stop_reason = ell_dbl_fft (ecmdata, xA, zA, xs, zs, Ad4);	/* S = 2*A */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xs, zs, xA, zA, xA, zA, xA, zA);/* A = S+A */
			if (stop_reason) return (stop_reason);
			d = (d-e)/3;
		} else {
			stop_reason = ell_add_fft (ecmdata, xB, zB, xC, zC, xA, zA, xC, zC);/* C = C-B */
			if (stop_reason) return (stop_reason);
			stop_reason = ell_dbl_fft (ecmdata, xB, zB, xB, zB, Ad4);	/* B = 2*B */
			if (stop_reason) return (stop_reason);
			e = e >> 1;
		}
	    }

	    stop_reason = ell_add_fft_last (ecmdata, xB, zB, xA, zA, xC, zC, xx, zz);	/* A = A+B */
	    if (stop_reason) return (stop_reason);

	    n = d;
	}
	gwfree (&ecmdata->gwdata, xA);
	gwfree (&ecmdata->gwdata, zA);
	gwfree (&ecmdata->gwdata, xB);
	gwfree (&ecmdata->gwdata, zB);
	gwfree (&ecmdata->gwdata, xC);
	gwfree (&ecmdata->gwdata, zC);
	gwfree (&ecmdata->gwdata, xt);
	gwfree (&ecmdata->gwdata, zt);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Multiplies the point (xx,zz) by n using a combination */
/* of ell_dbl and ell_add calls */

int bin_ell_mul (
	ecmhandle *ecmdata,
	gwnum	xx,
	gwnum	zz,
	uint64_t n,
	gwnum	Ad4)
{
	uint64_t c;
	unsigned long zeros;
	gwnum	xorg, zorg, xs, zs;
	int	stop_reason;

	xorg = gwalloc (&ecmdata->gwdata);
	if (xorg == NULL) goto oom;
	zorg = gwalloc (&ecmdata->gwdata);
	if (zorg == NULL) goto oom;
	xs = gwalloc (&ecmdata->gwdata);
	if (xs == NULL) goto oom;
	zs = gwalloc (&ecmdata->gwdata);
	if (zs == NULL) goto oom;

	for (zeros = 0; (n & 1) == 0; zeros++) n >>= 1;

	if (n > 1) {
		ell_begin_fft (ecmdata, xx, zz, xorg, zorg);

		c = 1; c <<= 63;
		while ((c&n) == 0) c >>= 1;
		c >>= 1;

		/* If the second bit is zero, we can save one ell_dbl call */

		if (c&n) {
			gwcopy (&ecmdata->gwdata, xorg, xx);
			gwcopy (&ecmdata->gwdata, zorg, zz);
			stop_reason = ell_dbl_fft (ecmdata, xx, zz, xs, zs, Ad4);
			if (stop_reason) return (stop_reason);
		} else {
			stop_reason = ell_dbl_fft (ecmdata, xorg, zorg, xx, zz, Ad4);
			if (stop_reason) return (stop_reason);
			stop_reason = ell_add_fft (ecmdata, xorg, zorg, xx, zz, xorg, zorg, xs, zs);
			if (stop_reason) return (stop_reason);
			c >>= 1;
		}

		/* Do the rest of the bits */

		do {
			if (c&n) {
				if (c == 1) {
					stop_reason = ell_add_fft_last (ecmdata, xs, zs, xx, zz,
									xorg, zorg, xx, zz);
					if (stop_reason) return (stop_reason);
				} else {
					stop_reason = ell_add_fft (ecmdata, xs, zs, xx, zz,
								   xorg, zorg, xx, zz);
					if (stop_reason) return (stop_reason);
					stop_reason = ell_dbl_fft (ecmdata, xs, zs, xs, zs, Ad4);
					if (stop_reason) return (stop_reason);
				}
			} else {
				stop_reason = ell_add_fft (ecmdata, xx, zz, xs, zs,
							   xorg, zorg, xs, zs);
				if (stop_reason) return (stop_reason);
				stop_reason = ell_dbl_fft (ecmdata, xx, zz, xx, zz, Ad4);
				if (stop_reason) return (stop_reason);
			}
			c >>= 1;
		} while (c);
	}

	gwfree (&ecmdata->gwdata, xorg); 
	gwfree (&ecmdata->gwdata, zorg); 
	gwfree (&ecmdata->gwdata, xs); 
	gwfree (&ecmdata->gwdata, zs); 

	while (zeros--) {
		stop_reason = ell_dbl (ecmdata, xx, zz, xx, zz, Ad4);
		if (stop_reason) return (stop_reason);
	}
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Try a series of Lucas chains to find the cheapest. */
/* First try v = (1+sqrt(5))/2, then (2+v)/(1+v), then (3+2*v)/(2+v), */
/* then (5+3*v)/(3+2*v), etc.  Finally, execute the cheapest. */
/* This is much faster than bin_ell_mul, but uses more memory. */

int ell_mul (
	ecmhandle *ecmdata,
	gwnum	xx,
	gwnum	zz,
	uint64_t n,
	gwnum	Ad4)
{
	unsigned long zeros;
	int	stop_reason;

	for (zeros = 0; (n & 1) == 0; zeros++) n >>= 1;

	if (n > 1) {
		unsigned long c, min;
		double	minv;

		minv = 0.6180339887498948;		/*v=(1+sqrt(5))/2*/
		min = lucas_cost (n, minv);

		c = lucas_cost (n, 0.7236067977499790);	/*(2+v)/(1+v)*/
		if (c < min) min = c, minv = 0.7236067977499790;

		c = lucas_cost (n, 0.5801787282954641);	/*(3+2*v)/(2+v)*/
		if (c < min) min = c, minv = 0.5801787282954641;

		c = lucas_cost (n, 0.6328398060887063);	/*(5+3*v)/(3+2*v)*/
		if (c < min) min = c, minv = 0.6328398060887063;

		c = lucas_cost (n, 0.6124299495094950);	/*(8+5*v)/(5+3*v)*/
		if (c < min) min = c, minv = 0.6124299495094950;

		c = lucas_cost (n, 0.6201819808074158);	/*(13+8*v)/(8+5*v)*/
		if (c < min) min = c, minv = 0.6201819808074158;

		c = lucas_cost (n, 0.6172146165344039);	/*(21+13*v)/(13+8*v)*/
		if (c < min) min = c, minv = 0.6172146165344039;

		c = lucas_cost (n, 0.6183471196562281);	/*(34+21*v)/(21+13*v)*/
		if (c < min) min = c, minv = 0.6183471196562281;

		c = lucas_cost (n, 0.6179144065288179);	/*(55+34*v)/(34+21*v)*/
		if (c < min) min = c, minv = 0.6179144065288179;

		c = lucas_cost (n, 0.6180796684698958);	/*(89+55*v)/(55+34*v)*/
		if (c < min) min = c, minv = 0.6180796684698958;

		stop_reason = lucas_mul (ecmdata, xx, zz, n, minv, Ad4);
		if (stop_reason) return (stop_reason);
	}
	while (zeros--) {
		stop_reason = ell_dbl (ecmdata, xx, zz, xx, zz, Ad4);
		if (stop_reason) return (stop_reason);
	}
	return (0);
}

/* Test if factor divides N, return TRUE if it does */

int testFactor (
	gwhandle *gwdata,
	struct work_unit *w,
	giant	f)		/* Factor to test */
{
	giant	tmp;
	int	divides_ok;

/* See if this is a valid factor */

	tmp = popg (&gwdata->gdata, f->sign + 5);	/* Allow room for mul by KARG */
	itog (w->b, tmp);
	powermod (tmp, w->n, f);
	dblmulg (w->k, tmp);
	iaddg (w->c, tmp);
	modgi (&gwdata->gdata, f, tmp);
	divides_ok = isZero (tmp);
	pushg (&gwdata->gdata, 1);
	if (!divides_ok) return (FALSE);

/* If QAing, see if we found the expected factor */

	if (QA_IN_PROGRESS) {
		tmp = popg (&gwdata->gdata, f->sign + 5);
		gtog (f, tmp);
		modg (QA_FACTOR, tmp);
		divides_ok = isZero (tmp);
		pushg (&gwdata->gdata, 1);
		if (!divides_ok) {
			char	buf[200];
			strcpy (buf, "ERROR: Factor not found. Expected ");
			gtoc (QA_FACTOR, buf+strlen(buf), 150);
			strcat (buf, "\n");
			OutputBoth (MAIN_THREAD_NUM, buf);
		}
	}

/* All done, return success */

	return (TRUE);
}

/* Set N, the number we are trying to factor */

int setN (
	gwhandle *gwdata,
	int	thread_num,
	struct work_unit *w,
	giant	*N)		/* k*b^n+c as a giant */
{
	unsigned long bits, p;
	char	buf[2500];

/* Create the binary representation of the number we are factoring */
/* Allocate 5 extra words to handle any possible k value. */

	bits = (unsigned long) (w->n * _log2 (w->b));
	*N = allocgiant ((bits >> 5) + 5);
	if (*N == NULL) return (OutOfMemory (thread_num));

/* This special code comes from Serge Batalov */

	if (IniGetInt (INI_FILE, "PhiExtensions", 0) &&
	    w->k == 1.0 && w->b == 2 && w->c == -1) {		/*=== this input means Phi(n,2) with n semiprime ===*/
		unsigned int i,k,q,knownSmallMers[] = {3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217,
						       4253, 4423, 9689, 9941, 11213, 19937, 999999999}; /* for now, just the cases where w->n = p * q, and 2^q-1 is prime */
		for (i=0; (q=knownSmallMers[i]) < w->n || q*q <= w->n; i++) if ((w->n%q) == 0) {
			giant tmp = allocgiant ((bits >> 5) + 5);
			if (!tmp) return (OutOfMemory (thread_num));
			ultog (1, tmp);
			ultog (1, *N);
			gshiftleft (w->n-w->n/q, *N);
			for (k=2; k < q; k++) {
				gshiftleft (w->n/q, tmp);
				addg (tmp, *N);
			}
			iaddg (1, *N);
			if (q != w->n/q) {
				ultog (w->b, tmp);
				power (tmp, w->n/q);
				iaddg (w->c, tmp);
				divg (tmp, *N);
			}
			free (tmp);
			/* w->minimum_fftlen = w->n; */ /*=== too late to do this here. Moved before gwsetup() --SB. */
			if (!w->known_factors || !strcmp (w->known_factors, "1")) {
				p = sprintf (buf, "M%lu", w->n/q); if(q != w->n/q) p += sprintf (buf+p, "/M%d", q);
				w->known_factors = (char *) malloc (p+1);
				memcpy (w->known_factors, buf, p+1);
			}
			return (0);
		}
        }

	if (IniGetInt (INI_FILE, "PhiExtensions", 0) &&
	    w->k == 1.0 && labs(w->c) == 1 && (w->n%3) == 0) {		/*=== this input means Phi(3,-b^(n/3)) ===*/
		giant	tmp = allocgiant ((bits >> 5) + 5);
		if (tmp == NULL) return (OutOfMemory (thread_num));
		ultog (w->b, tmp);
		power (tmp, w->n/3);
		gtog (tmp, *N);
		squareg (*N);
		if (w->c == 1) subg (tmp, *N); else addg (tmp, *N);
		iaddg (1, *N);
		free (tmp);
		/* w->minimum_fftlen = w->n; */ /*=== too late to do this here. Moved before gwsetup() --SB. */
		if (!w->known_factors) {
			p = sprintf (buf, "(%lu^%lu%+ld)", w->b, w->n/3, w->c);
			w->known_factors = (char *) malloc (p+1);
			memcpy (w->known_factors, buf, p+1);
		}
		return (0);
	}

/* Standard code for working on k*b^n+c */

	ultog (w->b, *N);
	power (*N, w->n);
	dblmulg (w->k, *N);
	iaddg (w->c, *N);

/* If we have a list of known factors then process it */

	if (w->known_factors != NULL) {
		char	*comma, *p;
		giant	tmp, f;

		tmp = allocgiant ((bits >> 5) + 5);
		if (tmp == NULL) return (OutOfMemory (thread_num));
		f = allocgiant ((bits >> 5) + 5);
		if (f == NULL) return (OutOfMemory (thread_num));

/* Process each factor */

		for (p = w->known_factors; ; ) {

/* Get the factor - raise error is it is less than or equal to one */

			ctog (p, f);
			if (gsign (f) < 1 || isone (f)) {
				char	msg[1000];
				strcpy (buf, p);
				buf[10] = 0;
				sprintf (msg, "Error parsing known factors of %s near: '%s'\n", gwmodulo_as_string (gwdata), buf);
				OutputStr (thread_num, msg);
				free (f);
				free (tmp);
				deleteWorkToDoLine (thread_num, w, FALSE);
				return (STOP_ABORT);
			}

/* Divide N by factor - then verify the factor */

			gtog (*N, tmp);
			divg (f, tmp);
			mulg (tmp, f);
			if (gcompg (f, *N)) {
				strcpy (buf, p);
				comma = strchr (buf, ',');
				if (comma != NULL) *comma = 0;
				sprintf (buf+strlen(buf), " does not divide %s\n", gwmodulo_as_string (gwdata));
				OutputBoth (thread_num, buf);
				free (f);
				free (tmp);
				deleteWorkToDoLine (thread_num, w, FALSE);
				return (STOP_ABORT);
			}
			gtog (tmp, *N);

/* Skip to next factor in list */

			comma = strchr (p, ',');
			if (comma == NULL) break;
			p = comma + 1;
		}
		free (f);
		free (tmp);
		return (0);
	}

/* Return success */

	return (0);
}

/* Do a GCD of the input value and N to see if a factor was found. */
/* The GCD is returned in factor iff a factor is found. */
/* This routine used to be interruptible and thus returns a stop_reason. */
/* Since switching to GMP's mpz code to implement the GCD this routine is no longer interruptible. */

int gcd (
	gwhandle *gwdata,
	int	thread_num,
	gwnum	gg,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found if any */
{
	giant	v;
	mpz_t	a, b;

/* Assume a factor will not be found */

	*factor = NULL;

/* Convert input number to binary */

	v = popg (&gwdata->gdata, ((int) gwdata->bit_length >> 5) + 10);
	if (v == NULL) goto oom;
	if (gwtogiant (gwdata, gg, v)) {	// On unexpected error, return no factor found
		pushg (&gwdata->gdata, 1);
		return (0);
	}

/* Do the GCD */

	mpz_init (a);
	mpz_init (b);
	gtompz (v, a);
	gtompz (N, b);
	pushg (&gwdata->gdata, 1);
	mpz_gcd (a, a, b);

/* If a factor was found, save it in FAC */

	if (mpz_cmp_ui (a, 1) && mpz_cmp (a, b)) {
		*factor = allocgiant ((int) divide_rounding_up (mpz_sizeinbase (a, 2), 32));
		if (*factor == NULL) goto oom;
		mpztog (a, *factor);
	}

/* Cleanup and return */

	mpz_clear (a);
	mpz_clear (b);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (thread_num));
}

/* Computes the modular inverse of a number.  This is done using the */
/* extended GCD algorithm.  If a factor is accidentally found, it is */
/* returned in factor.  Function returns stop_reason if it was */
/* interrupted by an escape. */

int ecm_modinv (
	ecmhandle *ecmdata,
	gwnum	b,
	giant	N,			/* Number we are factoring */
	giant	*factor)		/* Factor found, if any */
{
	giant	v;

/* Convert input number to binary */

	v = popg (&ecmdata->gwdata.gdata, ((int) ecmdata->gwdata.bit_length >> 5) + 10);
	if (v == NULL) goto oom;
	if (gwtogiant (&ecmdata->gwdata, b, v)) {
		// On unexpected, should-never-happen error, return out-of-memory for lack of a better error message
		goto oom;
	}

#ifdef MODINV_USING_GIANTS

	int	stop_reason;

/* Let the invg code use gwnum b's memory.  This code is slower, but at least it is interruptible. */
/* Compute 1/v mod N */

	gwfree_temporarily (&ecmdata->gwdata, b);
	stop_reason = invgi (&ecmdata->gwdata.gdata, ecmdata->thread_num, N, v);
	if (stop_reason == GIANT_OUT_OF_MEMORY)
		stop_reason = OutOfMemory (ecmdata->thread_num);
	gwrealloc_temporarily (&ecmdata->gwdata, b);
	if (stop_reason) {
		pushg (&ecmdata->gwdata.gdata, 1);
		return (stop_reason);
	}

/* If a factor was found, save it in FAC */

	if (v->sign < 0) {
		negg (v);
		*factor = allocgiant (v->sign);
		if (*factor == NULL) goto oom;
		gtog (v, *factor);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		*factor = NULL;
		gianttogw (&ecmdata->gwdata, v, b);
	}

/* Use the faster GMP library to do an extended GCD which gives us 1/v mod N */

#else
	{
	mpz_t	__v, __N, __gcd, __inv;

/* Do the extended GCD */

	mpz_init (__v);
	mpz_init (__N);
	mpz_init (__gcd);
	mpz_init (__inv);
	gtompz (v, __v);
	gtompz (N, __N);
	mpz_gcdext (__gcd, __inv, NULL, __v, __N);
	mpz_clear (__v);

/* If a factor was found (gcd != 1 && gcd != N), save it in FAC */

	if (mpz_cmp_ui (__gcd, 1) && mpz_cmp (__gcd, __N)) {
		*factor = allocgiant ((int) divide_rounding_up (mpz_sizeinbase (__gcd, 2), 32));
		if (*factor == NULL) goto oom;
		mpztog (__gcd, *factor);
	}

/* Otherwise, convert the inverse to FFT-ready form */

	else {
		*factor = NULL;
		if (mpz_sgn (__inv) < 0) mpz_add (__inv, __inv, __N);
		mpztog (__inv, v);
		gianttogw (&ecmdata->gwdata, v, b);
	}

/* Cleanup and return */

	mpz_clear (__gcd);
	mpz_clear (__inv);
	mpz_clear (__N);
	}
#endif

/* Clean up */

	pushg (&ecmdata->gwdata.gdata, 1);

/* Increment count and return */

	ecmdata->modinv_count++;
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}


/* Computes the modular inverse of an array of numbers */
/* Uses extra multiplications to make only one real modinv call */
/* Uses the simple formula 1/a = b * 1/ab, 1/b = a * 1/ab */
/* If we accidentally find a factor it is returned in factor. */
/* Return stop_reason if there is a user interrupt. */

int grouped_modinv (
	ecmhandle *ecmdata,
	gwnum	*b,
	unsigned int size,
	gwnum	*tmp,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	unsigned int i;
	gwnum	*orig_tmp;
	int	stop_reason;

/* Handle group of 1 as a special case */

	if (size == 1) return (ecm_modinv (ecmdata, *b, N, factor));

/* Handle an odd size */

	orig_tmp = tmp;
	if (size & 1) {
		gwswap (b[0], *tmp);
		tmp++;
	}

/* Multiply each pair of numbers */

	for (i = (size & 1); i < size; i += 2) {
		gwfft (&ecmdata->gwdata, b[i], b[i]);
		gwfft (&ecmdata->gwdata, b[i+1], b[i+1]);
		gwfftfftmul (&ecmdata->gwdata, b[i], b[i+1], *tmp);
		tmp++;
	}

/* Recurse */

	stop_reason = grouped_modinv (ecmdata, orig_tmp, (size+1) / 2, tmp, N, factor);
	if (stop_reason) return (stop_reason);
	if (*factor != NULL) return (0);

/* Handle an odd size */

	if (size & 1) {
		gwswap (b[0], *orig_tmp);
		orig_tmp++;
	}

/* Now perform multiplications on each pair to get the modular inverse */

	for (i = (size & 1); i < size; i += 2) {
		gwfft (&ecmdata->gwdata, *orig_tmp, *orig_tmp);
		gwfftfftmul (&ecmdata->gwdata, *orig_tmp, b[i], b[i]);
		gwfftfftmul (&ecmdata->gwdata, *orig_tmp, b[i+1], b[i+1]);
		gwswap (b[i], b[i+1]);
		orig_tmp++;
	}

/* All done, return */

	return (0);
}

/* Takes a point (a,b) and multiplies it by a value such that b will be one */
/* If we accidentally find a factor it is returned in factor.  Function */
/* returns stop_reason if it was interrupted by an escape. */

int normalize (
	ecmhandle *ecmdata,
	gwnum	a,
	gwnum	b,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	int	stop_reason;

/* Compute the modular inverse and scale up the first input value */

	stop_reason = ecm_modinv (ecmdata, b, N, factor);
	if (stop_reason) return (stop_reason);
	if (*factor != NULL) return (0);
	gwmul (&ecmdata->gwdata, b, a);
	return (0);
}

/* Adds a point (a,b) to the list of numbers that need normalizing. */
/* This is done in such a way as to minimize the amount of memory used. */

/* This is an interesting bit of code with a variety of algorithms */
/* available.  Assuming there are N pairs to normalize, then you can: */
/* 1) Use 3*N memory and use as few as 2 multiplies per pair. */
/* 2) Use 2*N memory and use 3 multiplies per pair. */
/* 3) Use N+log N memory and use O(log N) multiplies */
/* 4) Use N memory and use O(N^2) multiplies. */

int add_to_normalize_pool (
	ecmhandle *ecmdata,
	gwnum	a,
	gwnum	b,
	int	ffted)		/* TRUE if input arguments have been FFTed */
{

/* Switch off the type of pooling we are going to do */

	switch (ecmdata->pool_type) {

/* Implement algorithm 2 above */

	case POOL_3MULT:

/* If this is the first call allocate memory for the gwnum we use */

		if (ecmdata->pool_count == 0) {
			ecmdata->pool_modinv_value = gwalloc (&ecmdata->gwdata);
			if (ecmdata->pool_modinv_value == NULL) goto oom;
			gwcopy (&ecmdata->gwdata, b, ecmdata->pool_modinv_value);
			ecmdata->pool_ffted = ffted;
		}

/* Otherwise, multiply a by the accumulated b values */

		else if (ffted) {
			if (ecmdata->pool_count != 1)
				gwfft (&ecmdata->gwdata,
				       ecmdata->pool_modinv_value,
				       ecmdata->pool_modinv_value);
			gwfftfftmul (&ecmdata->gwdata,
				     ecmdata->pool_modinv_value, a, a);
			ecmdata->poolz_values[ecmdata->pool_count] = gwalloc (&ecmdata->gwdata);
			if (ecmdata->poolz_values[ecmdata->pool_count] == NULL) goto oom;
			gwcopy (&ecmdata->gwdata, b, ecmdata->poolz_values[ecmdata->pool_count]);
			gwfftfftmul (&ecmdata->gwdata,
				     ecmdata->poolz_values[ecmdata->pool_count],
				     ecmdata->pool_modinv_value,
				     ecmdata->pool_modinv_value);
		} else {
			gwmul (&ecmdata->gwdata, ecmdata->pool_modinv_value, a);
			ecmdata->poolz_values[ecmdata->pool_count] = gwalloc (&ecmdata->gwdata);
			if (ecmdata->poolz_values[ecmdata->pool_count] == NULL) goto oom;
			gwfft (&ecmdata->gwdata,
			       b, ecmdata->poolz_values[ecmdata->pool_count]);
			gwfftfftmul (&ecmdata->gwdata,
				     ecmdata->poolz_values[ecmdata->pool_count],
				     ecmdata->pool_modinv_value,
				     ecmdata->pool_modinv_value);
		}

/* Add a to array of values to normalize */

		ecmdata->pool_values[ecmdata->pool_count++] = a;
		break;

/* Implement algorithm 4 above */

	case POOL_N_SQUARED:

/* If this is the first call allocate memory for the gwnum we use */

		if (ecmdata->pool_count == 0) {
			ecmdata->pool_modinv_value = gwalloc (&ecmdata->gwdata);
			if (ecmdata->pool_modinv_value == NULL) goto oom;
			gwcopy (&ecmdata->gwdata, b, ecmdata->pool_modinv_value);
			ecmdata->pool_ffted = ffted;
		}

/* Otherwise, multiply a by the accumulated b values */
/* and multiply all previous a's by this b */

		else if (ffted) {
			unsigned int i;
			if (ecmdata->pool_count != 1)
				gwfft (&ecmdata->gwdata,
				       ecmdata->pool_modinv_value,
				       ecmdata->pool_modinv_value);
			gwfftfftmul (&ecmdata->gwdata,
				     ecmdata->pool_modinv_value, a, a);
			gwfftfftmul (&ecmdata->gwdata,
				     b, ecmdata->pool_modinv_value,
				     ecmdata->pool_modinv_value);
			for (i = 0; i < ecmdata->pool_count; i++)
				if (i == 0 && ecmdata->pool_ffted) {
					gwfftfftmul (&ecmdata->gwdata,
						     b, ecmdata->pool_values[i],
						     ecmdata->pool_values[i]);
					ecmdata->pool_ffted = FALSE;
				} else
					gwfftmul (&ecmdata->gwdata,
						  b, ecmdata->pool_values[i]);
		} else {
			unsigned int i;
			gwnum	tmp;
			gwmul (&ecmdata->gwdata, ecmdata->pool_modinv_value, a);
			tmp = gwalloc (&ecmdata->gwdata);
			if (tmp == NULL) goto oom;
			gwfft (&ecmdata->gwdata, b, tmp);
			gwfftfftmul (&ecmdata->gwdata, tmp, ecmdata->pool_modinv_value, ecmdata->pool_modinv_value);
			for (i = 0; i < ecmdata->pool_count; i++)
				if (i == 0 && ecmdata->pool_ffted) {
					gwfftfftmul (&ecmdata->gwdata, tmp, ecmdata->pool_values[i], ecmdata->pool_values[i]);
					ecmdata->pool_ffted = FALSE;
				} else
					gwfftmul (&ecmdata->gwdata, tmp, ecmdata->pool_values[i]);
			gwfree (&ecmdata->gwdata, tmp);
		}

/* Add a to array of values to normalize */

		ecmdata->pool_values[ecmdata->pool_count++] = a;
		break;
	}
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Takes each point from add_to_normalize_pool and normalizes it. */
/* If we accidentally find a factor, it is returned in factor. */
/* Return stop_reason if there is a user interrupt. */

int normalize_pool (
	ecmhandle *ecmdata,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	unsigned int i;
	int	stop_reason;

/* Compute the modular inverse */

	stop_reason = ecm_modinv (ecmdata, ecmdata->pool_modinv_value, N, factor);
	if (stop_reason) return (stop_reason);
	if (*factor != NULL) goto exit;

/* Now invert each value */
/* Switch off the type of pooling we are going to do */

	switch (ecmdata->pool_type) {

/* Implement algorithm 2 above */

	case POOL_3MULT:
		for (i = ecmdata->pool_count-1; ; i--) {
			if (i == 0 && ecmdata->pool_ffted) {
				gwfft (&ecmdata->gwdata,
				       ecmdata->pool_modinv_value,
				       ecmdata->pool_modinv_value);
				gwfftfftmul (&ecmdata->gwdata,
					     ecmdata->pool_modinv_value,
					     ecmdata->pool_values[i],
					     ecmdata->pool_values[i]);
			} else
				gwmul (&ecmdata->gwdata,
				       ecmdata->pool_modinv_value,
				       ecmdata->pool_values[i]);
			if (i == 0) break;
			gwfftfftmul (&ecmdata->gwdata,
				     ecmdata->poolz_values[i],
				     ecmdata->pool_modinv_value,
				     ecmdata->pool_modinv_value);
			gwfree (&ecmdata->gwdata, ecmdata->poolz_values[i]);
		}
		break;

/* Implement algorithm 4 above */

	case POOL_N_SQUARED:
		gwfft (&ecmdata->gwdata,
		       ecmdata->pool_modinv_value,
		       ecmdata->pool_modinv_value);
		for (i = 0; i < ecmdata->pool_count; i++)
			if (i == 0 && ecmdata->pool_ffted) {
				gwfftfftmul (&ecmdata->gwdata,
					     ecmdata->pool_modinv_value,
					     ecmdata->pool_values[i],
					     ecmdata->pool_values[i]);
			} else
				gwfftmul (&ecmdata->gwdata,
					  ecmdata->pool_modinv_value,
					  ecmdata->pool_values[i]);
		break;
	}

/* Cleanup and reinitialize */

exit:	ecmdata->pool_count = 0;
	gwfree (&ecmdata->gwdata, ecmdata->pool_modinv_value);
	return (0);
}


/* From R. P. Brent, priv. comm. 1996:
Let s > 5 be a pseudo-random seed (called $\sigma$ in the Tech. Report),

	u/v = (s^2 - 5)/(4s)

Then starting point is (x_1, y_1) where

	x_1 = (u/v)^3
and
	a = (v-u)^3(3u+v)/(4u^3 v) - 2
*/
int choose12 (
	ecmhandle *ecmdata,
	struct work_unit *w,
	gwnum 	x,
	gwnum 	z,
	double 	curve,
	gwnum	*Ad4,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	gwnum	xs, zs, t1, t2, t3;
	int	stop_reason;

	xs = gwalloc (&ecmdata->gwdata);
	if (xs == NULL) goto oom;
	zs = gwalloc (&ecmdata->gwdata);
	if (zs == NULL) goto oom;
	t1 = gwalloc (&ecmdata->gwdata);
	if (t1 == NULL) goto oom;
	t2 = gwalloc (&ecmdata->gwdata);
	if (t2 == NULL) goto oom;
	t3 = gwalloc (&ecmdata->gwdata);
	if (t3 == NULL) goto oom;

	dbltogw (&ecmdata->gwdata, curve, zs);
	gwcopy (&ecmdata->gwdata, zs, xs);
	gwsquare (&ecmdata->gwdata, xs);	/* s^2 */
	dbltogw (&ecmdata->gwdata, 5.0, t1);
	gwsub (&ecmdata->gwdata, t1, xs);	/* u = s^2 - 5 */
	dbltogw (&ecmdata->gwdata, 4.0, t1);
	gwmul (&ecmdata->gwdata, t1, zs);	/* v = 4*s */
	gwcopy (&ecmdata->gwdata, xs, x);
	gwsquare (&ecmdata->gwdata, x);
	gwsafemul (&ecmdata->gwdata, xs, x);	/* x = u^3 */
	gwcopy (&ecmdata->gwdata, zs, z);
	gwsquare (&ecmdata->gwdata, z);
	gwsafemul (&ecmdata->gwdata, zs, z);	/* z = v^3 */

	/* Now for A. */
	gwcopy (&ecmdata->gwdata, zs, t2);
	gwsub (&ecmdata->gwdata, xs, t2);
	gwcopy (&ecmdata->gwdata, t2, t3);
	gwsquare (&ecmdata->gwdata, t2);
	gwmul (&ecmdata->gwdata, t3, t2);	/* (v-u)^3 */
	gwcopy (&ecmdata->gwdata, xs, t3);
	gwadd (&ecmdata->gwdata, t3, t3);
	gwadd (&ecmdata->gwdata, xs, t3);
	gwadd (&ecmdata->gwdata, zs, t3);
	gwmul (&ecmdata->gwdata, t3, t2);	/* An = (v-u)^3 (3u+v) */
	gwcopy (&ecmdata->gwdata, zs, t3);
	gwsafemul (&ecmdata->gwdata, xs, t3);
	gwsquare (&ecmdata->gwdata, xs);
	gwsafemul (&ecmdata->gwdata, xs, t3);
	*Ad4 = gwalloc (&ecmdata->gwdata);
	if (*Ad4 == NULL) goto oom;
	dbltogw (&ecmdata->gwdata, 4.0, *Ad4);
	gwmul (&ecmdata->gwdata, t3, *Ad4);	/* An/Ad is now A + 2 */

	/* Normalize so that An is one */
	stop_reason = normalize (ecmdata, *Ad4, t2, N, factor);
	if (stop_reason) return (stop_reason);

	/* For extra speed, precompute Ad * 4 */
	dbltogw (&ecmdata->gwdata, 4.0, t1);
	gwmul (&ecmdata->gwdata, t1, *Ad4);

	/* Even more speed, save FFT of Ad4 */
	gwfft (&ecmdata->gwdata, *Ad4, *Ad4);

/* Clean up temporaries */

	gwfree (&ecmdata->gwdata, xs);
	gwfree (&ecmdata->gwdata, zs);
	gwfree (&ecmdata->gwdata, t1);
	gwfree (&ecmdata->gwdata, t2);
	gwfree (&ecmdata->gwdata, t3);
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}

/* Print message announcing the start of this curve */

void curve_start_msg (
	ecmhandle *ecmdata,
	int	thread_num,
	unsigned long curve,
	double	 sigma,
	uint64_t B,
	uint64_t C)
{
	char	buf[120];

	sprintf (buf, "%s ECM curve #%ld",
		 gwmodulo_as_string (&ecmdata->gwdata), curve);
	title (thread_num, buf);

	sprintf (buf,
		 "ECM on %s: curve #%ld with s=%.0f, B1=%.0f, B2=%.0f\n",
		 gwmodulo_as_string (&ecmdata->gwdata), curve, sigma,
		 (double) B, (double) C);
	OutputStr (thread_num, buf);
}

/* These routines manage the computing of Q^m in stage 2 */

int mQ_init (
	ecmhandle *ecmdata,
	gwnum	x,
	uint64_t m,
	gwnum	Q2Dx,
	gwnum	Ad4)
{
	int	stop_reason;

	ecmdata->Qprevmx = gwalloc (&ecmdata->gwdata);
	if (ecmdata->Qprevmx == NULL) goto oom;
	ecmdata->Qprevmz = gwalloc (&ecmdata->gwdata);
	if (ecmdata->Qprevmz == NULL) goto oom;
	ecmdata->Qmx = gwalloc (&ecmdata->gwdata);
	if (ecmdata->Qmx == NULL) goto oom;
	ecmdata->Qmz = gwalloc (&ecmdata->gwdata);
	if (ecmdata->Qmz == NULL) goto oom;
	gwcopy (&ecmdata->gwdata, x, ecmdata->Qprevmx);
	dbltogw (&ecmdata->gwdata, 1.0, ecmdata->Qprevmz);
	stop_reason = bin_ell_mul (ecmdata, ecmdata->Qprevmx, ecmdata->Qprevmz, m - 4*ecmdata->D, Ad4);
	if (stop_reason) return (stop_reason);
	gwfft (&ecmdata->gwdata, ecmdata->Qprevmx, ecmdata->Qprevmx);
	gwfft (&ecmdata->gwdata, ecmdata->Qprevmz, ecmdata->Qprevmz);
	gwcopy (&ecmdata->gwdata, x, ecmdata->Qmx);
	dbltogw (&ecmdata->gwdata, 1.0, ecmdata->Qmz);
	stop_reason = bin_ell_mul (ecmdata, ecmdata->Qmx, ecmdata->Qmz, m - 2*ecmdata->D, Ad4);
	if (stop_reason) return (stop_reason);
	gwfft (&ecmdata->gwdata, ecmdata->Qmx, ecmdata->Qmx);
	gwfft (&ecmdata->gwdata, ecmdata->Qmz, ecmdata->Qmz);

	/* There will be no more ell_dbl calls */
	gwfree (&ecmdata->gwdata, Ad4);

	/* Precompute the FFTs of Q2Dx+1 and Q2Dx-1 */
	ecmdata->Q2Dxplus1 = Q2Dx;
	ecmdata->Q2Dxminus1 = gwalloc (&ecmdata->gwdata);
	if (ecmdata->Q2Dxminus1 == NULL) goto oom;
	gwaddsmall (&ecmdata->gwdata, Q2Dx, -1);
	gwfft (&ecmdata->gwdata, Q2Dx, ecmdata->Q2Dxminus1);
	gwaddsmall (&ecmdata->gwdata, Q2Dx, 2);
	gwfft (&ecmdata->gwdata, Q2Dx, ecmdata->Q2Dxplus1);

	/* Init the arrays used in pooled normalizes of mQx values */
	if (ecmdata->TWO_FFT_STAGE2) {
		unsigned long i;
		for (i = 0; i < ecmdata->E; i++) {
			ecmdata->mQx[i] = gwalloc (&ecmdata->gwdata);
			if (ecmdata->mQx[i] == NULL) goto oom;
		}
		ecmdata->mQx_count = 0;
	}
	return (0);

/* Out of memory exit path */

oom:	return (OutOfMemory (ecmdata->thread_num));
}
int mQ_next (
	ecmhandle *ecmdata,
	gwnum	*retx,
	gwnum	*retz,
	giant	N,		/* Number we are factoring */
	giant	*factor)	/* Factor found, if any */
{
	int	stop_reason;

/* The non-normalized case - simple multiply the last Q^m value */
/* by Q^2D to get the next Q^m value */

	if (!ecmdata->TWO_FFT_STAGE2) {
		stop_reason = ell_add_special3 (ecmdata, ecmdata->Qmx, ecmdata->Qmz,
					        ecmdata->Q2Dxplus1, ecmdata->Q2Dxminus1,
						&ecmdata->Qprevmx, &ecmdata->Qprevmz);
		if (stop_reason) return (stop_reason);
		gwswap (ecmdata->Qmx, ecmdata->Qprevmx);
		gwswap (ecmdata->Qmz, ecmdata->Qprevmz);
		gwfft (&ecmdata->gwdata, ecmdata->Qmx, ecmdata->Qmx);
		gwfft (&ecmdata->gwdata, ecmdata->Qmz, ecmdata->Qmz);
		*retx = ecmdata->Qmx;
		*retz = ecmdata->Qmz;
		return (0);
	}

/* The normalized case - batch up a bunch of Q^m values and normalize */
/* them.  Then return them one at a time.  Obviously retz need not be */
/* returned since it is always one. */

	if (ecmdata->mQx_count == 0) {
		for ( ; ecmdata->mQx_count < ecmdata->E; ecmdata->mQx_count++) {
			stop_reason = ell_add_special3 (ecmdata, ecmdata->Qmx, ecmdata->Qmz,
						        ecmdata->Q2Dxplus1, ecmdata->Q2Dxminus1,
						        &ecmdata->Qprevmx, &ecmdata->Qprevmz);
			if (stop_reason) return (stop_reason);
			gwswap (ecmdata->Qmx, ecmdata->Qprevmx);
			gwswap (ecmdata->Qmz, ecmdata->Qprevmz);
			gwfft (&ecmdata->gwdata, ecmdata->Qmx, ecmdata->Qmx);
			gwfft (&ecmdata->gwdata, ecmdata->Qmz, ecmdata->Qmz);
			gwcopy (&ecmdata->gwdata, ecmdata->Qmx, ecmdata->mQx[ecmdata->mQx_count]);
			stop_reason = add_to_normalize_pool (ecmdata, ecmdata->mQx[ecmdata->mQx_count], ecmdata->Qmz, 1);
			if (stop_reason) return (stop_reason);
		}
		stop_reason = normalize_pool (ecmdata, N, factor);
		if (stop_reason) return (stop_reason);
		if (*factor != NULL) return (0);
	}
	*retx = ecmdata->mQx[ecmdata->E - ecmdata->mQx_count];
	gwfft (&ecmdata->gwdata, *retx, *retx);
	ecmdata->mQx_count--;
	return (0);
}
void mQ_term (
	ecmhandle *ecmdata)
{
	gwfree (&ecmdata->gwdata, ecmdata->Qprevmx);
	gwfree (&ecmdata->gwdata, ecmdata->Qprevmz);
	gwfree (&ecmdata->gwdata, ecmdata->Qmx);
	gwfree (&ecmdata->gwdata, ecmdata->Qmz);
	gwfree (&ecmdata->gwdata, ecmdata->Q2Dxminus1);
}

/* Record the amount of memory being used by this thread.  Until we get to */
/* stage 2, ECM uses 13 gwnums (see comments in the code). */

void ecm_stage1_memory_usage (
	int	thread_num,
	ecmhandle *ecmdata)
{
	set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&ecmdata->gwdata,  13));
}

/* Choose 4 FFT stage 2 of the 2 FFT stage 2.  Also choose a good */
/* value for D and a good algorithm for normalize_pool. */
/* We try to choose the above such that the number of multiplications */
/* are minimized, yet too much memory isn't used. */

int choose_stage2_plan (
	int	thread_num,
	ecmhandle *ecmdata,
	uint64_t B,			/* Stage 1 bound */
	uint64_t C)			/* Stage 2 bound */
{
	unsigned int memory;		/* MB of memory we can use */
	unsigned long numvals, d, e;
	unsigned long relprime, beste, bestnumvals;
	double	numprimes, numpairings, numsections, numgcdsections;
	double	cost, bestcost, density, gcd_cost;
	char	buf[120];
	int	stop_reason;

/* Get available memory.  We need 21 gwnums to do the smallest stage 2, plus */
/* another 0.5 gwnums of overhead in invgi (and gcd).  We assume 120 gwnums */
/* will allow us to do a reasonable efficient stage 2 implementation. */

replan:	stop_reason = avail_mem (thread_num,
				 cvt_gwnums_to_mem (&ecmdata->gwdata, 21.5),
				 cvt_gwnums_to_mem (&ecmdata->gwdata, 120),
				 &memory);
	if (stop_reason) return (stop_reason);
	if (memory < 8) memory = 8;

/* Define constants for the number of transforms for various operations */
/* The GCD cost is based on our timings and an Excel spreadsheet */

#define ELL_ADD_COST		12
#define N_SQUARED_POOL_COST	2
#define MULT3_POOL_COST		7
	gcd_cost = 861.0 * log ((double) ecmdata->gwdata.n) - 7775.0;
	if (gcd_cost < 100.0) gcd_cost = 100.0;

/* Reserve space for overhead needed by invgi and gcdg (about 1.5 gwnums) */
/* However, not all stage 2 implementations use gcds.  So only reserve space */
/* for the giant sin/cos data that is never freed once allocated (about 0.5 gwnums) */	

	memory = memory - (unsigned int) (gwnum_size (&ecmdata->gwdata) * 0.5 / 1048576.0);

/* Figure out how many gwnum values fit in our memory limit. */

	numvals = cvt_mem_to_gwnums (&ecmdata->gwdata, memory);
	ASSERTG (numvals >= 21);

/* If memory is really tight, then the 4 FFT - O(n^2) pooling is the */
/* most memory efficient ECM implementation.  This will be our default */
/* plan.  Note: D=30 (8 nQx values) requires 21 gwnums.  The next D */
/* value (60) requires 28 gwnums. */

	if (numvals < 28) {
		ecmdata->D = 30;
		ecmdata->E = 0;
		ecmdata->TWO_FFT_STAGE2 = FALSE;
		ecmdata->pool_type = POOL_N_SQUARED;
		bestnumvals = 21;
	}

/* Numprimes below C approximately equals C / (ln(C)-1) */
/* Compute numprimes between B and C */

	numprimes = ceil ((C / (log ((double) C) - 1)) - (B / (log ((double) B) - 1)));

/* Figure out the best value for E when using the O(N^2) pool method */

	beste = (unsigned long) sqrt (gcd_cost / N_SQUARED_POOL_COST) + 1;

/* Loop through various D values choosing the most cost effective one */

	bestcost = 1.0E99;
	d = ((unsigned long) sqrt ((double) (C-B)) / 2310 + 3) * 2310;
	for ( ; ; ) {
		if (d >= 2310) {
			relprime = d / 2310 * 480;
			density = 480.0 / 2310.0;
		} else if (d >= 210) {
			relprime = d / 210 * 48;
			density = 48.0 / 210.0;
		} else {
			relprime = d / 30 * 8;
			density = 8.0 / 30.0;
		}

/* Half the primes are eligible for pairing (numprimes / 2). */
/* The chance that a pairing occurs is numprimes / area.  Area would */
/* normally be C-B.  However, the relprime algorithm makes */
/* our primes much denser than that. */

		numpairings = ceil (
			(numprimes / 2.0 * numprimes / ((C-B) * density)));

/* There will be (C-B)/2D sections */

		numsections = ceil ((double) (C-B) / (double) (d+d));

/* Cost out the 4FFT stage 2 using this D			*/
/* The cost will be:						*/
/*	D/2 ell_add_ffts  + pool_cost (relprime) +		*/
/*	(C-B)/2D ell_add_specials + #primes*4			*/
/* The memory consumed will be:					*/
/*	13 + relprime gwnums if N^2 pooling			*/
/* or	13 + 2*relprime gwnums if 3N pooling			*/
/* Note that MQ_init requires B is at least 4 times D		*/

		if (B >= 4*d && 13 + relprime <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 1)) {
			cost = d/2 * ELL_ADD_COST +
			       relprime * relprime * N_SQUARED_POOL_COST +
			       numsections * ELL_ADD_COST +
			       (numprimes - numpairings) * 4;
			if (cost < bestcost) {
				ecmdata->TWO_FFT_STAGE2 = FALSE;
				ecmdata->pool_type = POOL_N_SQUARED;
				ecmdata->D = d;
				ecmdata->E = 0;
				bestcost = cost;
				bestnumvals = 13 + relprime; 
			}
		}
		if (B >= 4*d && 13 + relprime*2 <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 2)) {
			cost = d/2 * ELL_ADD_COST +
			       relprime * MULT3_POOL_COST +
			       numsections * ELL_ADD_COST +
			       (numprimes - numpairings) * 4;
			if (cost < bestcost) {
				ecmdata->TWO_FFT_STAGE2 = FALSE;
				ecmdata->pool_type = POOL_3MULT;
				ecmdata->D = d;
				ecmdata->E = 0;
				bestcost = cost;
				bestnumvals = 13 + relprime*2; 
			}
		}

/* Cost out the 2FFT stage 2 using this D			*/
/* The cost will be:						*/
/*	D/2 ell_add_ffts  + pool_cost (relprime) +		*/
/*	(C-B)/2D ell_add_specials + #primes*2 +			*/
/*	(C-B)/2D/E * pool_cost (e)				*/
/*	(C-B)/2D/E * gcd_cost					*/
/* The memory consumed will be:					*/
/*	13 + relprime + e gwnums if N^2 pooling			*/
/* or	13 + relprime + e + max (relprime, e) gwnums if 3N pooling */

		if (B >= 4*d && 13 + relprime <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 3)) {
			e = numvals - relprime - 13;
			if (e == 0) e = 1;
			if (e > beste) e = beste;
			numgcdsections = ceil (numsections / e);
			cost = d/2 * ELL_ADD_COST +
			       relprime * relprime * N_SQUARED_POOL_COST +
			       numsections * ELL_ADD_COST +
			       (numprimes - numpairings) * 2 +
			       numgcdsections * e * e * N_SQUARED_POOL_COST +
			       numgcdsections * gcd_cost;
			if (cost < bestcost) {
				ecmdata->TWO_FFT_STAGE2 = TRUE;
				ecmdata->pool_type = POOL_N_SQUARED;
				ecmdata->D = d;
				ecmdata->E = e;
				bestcost = cost;
				bestnumvals = 13 + relprime + e;
			}
		}
		if (B >= 4*d && 13 + relprime*2 <= numvals &&
		    (QA_TYPE == 0 || QA_TYPE == 4)) {
			e = numvals - 13 - relprime*2;
			if (e > relprime) e -= (e - relprime) / 2;
			if (e == 0) e = 1;
			numgcdsections = ceil (numsections / e);
			e = (unsigned long) ceil (numsections / numgcdsections);
			cost = d/2 * ELL_ADD_COST +
			       relprime * MULT3_POOL_COST +
			       numsections * (ELL_ADD_COST + 1) +
			       (numprimes - numpairings) * 2.0 +
			       numgcdsections * e * MULT3_POOL_COST +
			       numgcdsections * gcd_cost;
			if (cost < bestcost) {
				ecmdata->TWO_FFT_STAGE2 = TRUE;
				ecmdata->pool_type = POOL_3MULT;
				ecmdata->D = d;
				ecmdata->E = e;
				bestcost = cost;
				bestnumvals = 13 + relprime + e +
						(relprime > e ? relprime : e);
			}
		}

/* Cost out the next possible value of D */

		if (d > 2310) d = d - 2310;
		else if (d > 210) d = d - 210;
		else if (d > 30) d = d - 30;
		else break;
	}

/* Record the amount of memory this thread will be using in stage 2. */

	memory = cvt_gwnums_to_mem (&ecmdata->gwdata, bestnumvals);
	if (set_memory_usage (thread_num, MEM_VARIABLE_USAGE, memory))
		goto replan;

/* Output a useful message regarding memory usage */

	sprintf (buf, "Using %dMB of memory in stage 2.\n", memory);
	OutputStr (thread_num, buf);
	return (0);
}

/* Routines to create and read save files for an ECM factoring job */

#define ECM_MAGICNUM	0x1725bcd9
#define ECM_VERSION	1
#define ECM_STAGE1	0
#define ECM_STAGE2	1

void ecm_save (
	ecmhandle *ecmdata,
	writeSaveFileState *write_save_file_state,
	struct work_unit *w,
	int	stage,
	unsigned long curve,
	double	sigma,
	uint64_t B,
	uint64_t B_processed,
	uint64_t C_processed,
	gwnum	x,
	gwnum	gg)
{
	int	fd;
	unsigned long sum = 0;

/* Create the intermediate file */

	fd = openWriteSaveFile (write_save_file_state);
	if (fd < 0) return;

/* Write the file header. */

	if (! write_header (fd, ECM_MAGICNUM, ECM_VERSION, w)) goto writeerr;

/* Write the file data */

	if (! write_long (fd, stage, &sum)) goto writeerr;
	if (! write_long (fd, curve, &sum)) goto writeerr;
	if (! write_double (fd, sigma, NULL)) goto writeerr;
	if (! write_longlong (fd, B, &sum)) goto writeerr;
	if (! write_longlong (fd, B_processed, &sum)) goto writeerr;
	if (! write_longlong (fd, C_processed, &sum)) goto writeerr;

/* Write the data values */

	if (! write_gwnum (fd, &ecmdata->gwdata, x, &sum)) goto writeerr;
	if (! write_gwnum (fd, &ecmdata->gwdata, gg, &sum)) goto writeerr;

/* Write the checksum, we're done */

	if (! write_checksum (fd, sum)) goto writeerr;

	closeWriteSaveFile (write_save_file_state, fd);
	return;

/* An error occured.  Close and delete the current file. */

writeerr:
	deleteWriteSaveFile (write_save_file_state, fd);
}

/* Read a save file */

int old_ecm_restore (			/* For version 24 save files */
	ecmhandle *ecmdata,
	int	fd,
	int	*stage,
	unsigned long *curve,
	double	*sigma,
	uint64_t *B,
	uint64_t *B_processed,
	uint64_t *C_processed,
	gwnum	x,
	gwnum	gg)
{
	unsigned long magicnum, version;
	unsigned long tmp;
	unsigned long sum = 0, i;

/* Read the file header */

	_lseek (fd, 0, SEEK_SET);
	if (! read_long (fd, &magicnum, NULL)) return (FALSE);
	if (magicnum != 0x1a2b3cd4) return (FALSE);

	if (!read_long (fd, &version, NULL)) return (FALSE);
	if (version != 1) return (FALSE);

/* Read the file data */

	if (! read_long (fd, &tmp, &sum)) return (FALSE);
	*stage = (int) tmp;
	if (! read_long (fd, curve, &sum)) return (FALSE);
	if (! read_double (fd, sigma, NULL)) return (FALSE);
	if (! read_long (fd, &i, &sum)) return (FALSE);
	*B = i;
	if (! read_long (fd, &i, &sum)) return (FALSE);
	*B_processed = i;
	if (! read_long (fd, &i, &sum)) return (FALSE);
	*C_processed = i;

/* Read the values */

	if (! read_gwnum (fd, &ecmdata->gwdata, x, &sum)) return (FALSE);
	if (! read_gwnum (fd, &ecmdata->gwdata, gg, &sum)) return (FALSE);

/* Read and compare the checksum */

	if (! read_long (fd, &i, NULL)) return (FALSE);
	if (i != sum) return (FALSE);
	_close (fd);
	return (TRUE);
}

int ecm_restore (			/* For version 25 save files */
	ecmhandle *ecmdata,
	int	thread_num,
	char	*filename,
	struct work_unit *w,
	int	*stage,
	unsigned long *curve,
	double	*sigma,
	uint64_t *B,
	uint64_t *B_processed,
	uint64_t *C_processed,
	gwnum	x,
	gwnum	gg)
{
	int	fd;
	unsigned long version;
	unsigned long tmp;
	unsigned long sum = 0, filesum;

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto err;

/* Read the file header */

	if (! read_magicnum (fd, ECM_MAGICNUM)) {
		if (! old_ecm_restore (ecmdata, fd, stage, curve, sigma,
				       B, B_processed, C_processed, x, gg))
			goto readerr;
		return (TRUE);
	}
	if (! read_header (fd, &version, w, &filesum)) goto readerr;
	if (version != ECM_VERSION) goto readerr;

/* Read the file data */

	if (! read_long (fd, &tmp, &sum)) goto readerr;
	*stage = (int) tmp;
	if (! read_long (fd, curve, &sum)) goto readerr;
	if (! read_double (fd, sigma, NULL)) goto readerr;
	if (! read_longlong (fd, B, &sum)) goto readerr;
	if (! read_longlong (fd, B_processed, &sum)) goto readerr;
	if (! read_longlong (fd, C_processed, &sum)) goto readerr;

/* Read the values */

	if (! read_gwnum (fd, &ecmdata->gwdata, x, &sum)) goto readerr;
	if (! read_gwnum (fd, &ecmdata->gwdata, gg, &sum)) goto readerr;

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);
	return (TRUE);

/* An error occured.  Cleanup and return FALSE. */

readerr:
	_close (fd);
err:
	return (FALSE);
}


/**************************************************************
 *
 *	Main ECM Function
 *
 **************************************************************/

int ecm (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	ecmhandle ecmdata;
	uint64_t B;		/* Stage 1 bound */
	uint64_t C_start;	/* Stage 2 starting point (usually B) */
	uint64_t C;		/* Stage 2 ending point */
	uint64_t sieve_start, prime, m;
	unsigned long SQRT_B;
	double	sigma, last_output, last_output_t, one_over_B, one_over_C_minus_B;
	double	output_frequency, output_title_frequency;
	unsigned long i, j, curve, min_memory;
	readSaveFileState read_save_file_state;	/* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	char	filename[32], buf[255], JSONbuf[4000], fft_desc[200];
	int	res, stop_reason, stage, first_iter_msg;
	gwnum	x, z, t1, t2, gg;
	gwnum	Q2x, Q2z, Qiminus2x, Qiminus2z, Qdiffx, Qdiffz;
	giant	N;		/* Number being factored */
	giant	factor;		/* Factor found, if any */
	gwnum	Ad4 = NULL;
	int	msglen, continueECM, prpAfterEcmFactor;
	char	*str, *msg;
	double	timers[10];

/* Init local copies of B1 and B2 */

	B = (uint64_t) w->B1;
	C_start = (uint64_t) w->B2_start;
	C = (uint64_t) w->B2;

/* Clear pointers to allocated memory */

	memset (&ecmdata, 0, sizeof (ecmhandle));
	N = NULL;
	factor = NULL;
	str = NULL;
	msg = NULL;

/* Clear all timers */

restart:
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Time the giants squaring and multiply code in order to select the */
/* best crossover points.  This should only be done in the release code */
/* (optimized giants library). */

/*#define TIMING1*/
#ifdef TIMING1
if (w->n == 598) {
int i, j;
giant	x, y, z, a, m;
#define TESTSIZE	200
RDTSC_TIMING = 12;
x = allocgiant(TESTSIZE); y = allocgiant (2*TESTSIZE);
z = allocgiant (2*TESTSIZE), a = allocgiant (2*TESTSIZE);
m = allocgiant (TESTSIZE);
srand ((unsigned) time (NULL));
for (i = 0; i < TESTSIZE; i++) {
	x->n[i] = (rand () << 17) + rand ();
	m->n[i] = (rand () << 17) + rand ();
}
x->n[TESTSIZE-1] &= 0x00FFFFFF;
m->n[TESTSIZE-1] &= 0x00FFFFFF;
for (i = TESTSIZE; i >= 10; i--) {
	x->sign = i;
	m->sign = i;
	setmulmode (GRAMMAR_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, y);
		start_timer (timers, 0);
		if (B&1) mulg (m, y);
		else squareg (y);
		end_timer (timers, 0);
		if (timers[1] == 0 || timers[1] > timers[0]) timers[1] = timers[0];
		timers[0] = 0;
	}
	setmulmode (KARAT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, z);
		start_timer (timers, 0);
		if (B&1) mulg (m, z);
		else squareg (z);
		end_timer (timers, 0);
		if (timers[2] == 0 || timers[2] > timers[0]) timers[2] = timers[0];
		timers[0] = 0;
	}
	setmulmode (FFT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, a);
		start_timer (timers, 0);
		if (B&1) mulg (m, a);
		else squareg (a);
		end_timer (timers, 0);
		if (timers[3] == 0 || timers[3] > timers[0]) timers[3] = timers[0];
		timers[0] = 0;
	}
	sprintf (buf, "Size: %ld  G: ", i);
	print_timer (timers, 1, buf, TIMER_MS | TIMER_CLR);
	strcat (buf, ", K: ");
	print_timer (timers, 2, buf, TIMER_MS | TIMER_CLR);
	strcat (buf, ", F: ");
	print_timer (timers, 3, buf, TIMER_MS | TIMER_NL | TIMER_CLR);
	OutputBoth (thread_num, buf);
	if (gcompg (y, z) != 0)
		i--;
	if (gcompg (y, a) != 0)
		i--;
	Sleep (100);
}
return 0;
}

/* This code lets us time various giants FFT squarings and multiplies */

if (w->n == 601) {
int i, j;
giant	x, a, m;
#define TESTSIZE2	260000
RDTSC_TIMING = 12;
x = allocgiant(TESTSIZE2);
a = allocgiant (2*TESTSIZE2);
m = allocgiant (TESTSIZE2);
srand ((unsigned) time (NULL));
for (i = 0; i < TESTSIZE2; i++) {
	x->n[i] = (rand () << 17) + rand ();
	m->n[i] = (rand () << 17) + rand ();
}
x->n[TESTSIZE2-1] &= 0x00FFFFFF;
m->n[TESTSIZE2-1] &= 0x00FFFFFF;
for (i = 30; i < TESTSIZE2/2; i<<=1) {
	x->sign = i;
	m->sign = i;
	setmulmode (FFT_MUL);
	for (j = 0; j < 10; j++) {
		gtog (x, a);
		start_timer (timers, 0);
		if (B&1) mulg (m, a);
		else squareg (a);
		end_timer (timers, 0);
		if (timers[3] == 0 || timers[3] > timers[0]) timers[3] = timers[0];
		timers[0] = 0;
	}
	sprintf (buf, "Size: %ld  , F: ", i);
	print_timer (timers, 3, buf, TIMER_NL | TIMER_CLR | TIMER_MS);
	OutputStr (thread_num, buf);
	Sleep (100);
}
return 0;
}
#endif

/* Include timing code when building the debug version of prime95 */

#ifdef GDEBUG
if (w->n == 600) {
gwhandle gwdata;
void *workbuf;
int j, min_test, max_test, test, cnt, NUM_X87_TESTS, NUM_SSE2_TESTS, NUM_AVX_TESTS, NUM_AVX512_TESTS;
#define timeit(a,n,w) (((void**)a)[0]=w,((uint32_t*)a)[2]=n,gwtimeit(a))

gwinit (&gwdata); gwdata.cpu_flags &= ~CPU_AVX512F;
gwsetup (&gwdata, 1.0, 2, 10000000, -1);
workbuf = (void *) aligned_malloc (400000000, 4096);
memset (workbuf, 0, 400000000);
RDTSC_TIMING = 2;
min_test = IniGetInt (INI_FILE, "MinTest", 0);
max_test = IniGetInt (INI_FILE, "MaxTest", min_test);
NUM_X87_TESTS = timeit (gwdata.asm_data, -1, NULL);
NUM_SSE2_TESTS = timeit (gwdata.asm_data, -2, NULL);
NUM_AVX_TESTS = timeit (gwdata.asm_data, -3, NULL);
NUM_AVX512_TESTS = timeit (gwdata.asm_data, -4, NULL);
//SetThreadPriority (CURRENT_THREAD, THREAD_PRIORITY_TIME_CRITICAL);
for (j = 0; j < NUM_X87_TESTS + NUM_SSE2_TESTS + NUM_AVX_TESTS + NUM_AVX512_TESTS; j++) {
	cnt = 0;
	test = (j < NUM_X87_TESTS ? j :
		j < NUM_X87_TESTS + NUM_SSE2_TESTS ? 1000 + j - NUM_X87_TESTS :
		j < NUM_X87_TESTS + NUM_SSE2_TESTS + NUM_AVX_TESTS ? 2000 + j - NUM_X87_TESTS - NUM_SSE2_TESTS :
			3000 + j - NUM_X87_TESTS - NUM_SSE2_TESTS - NUM_AVX_TESTS);
	if (test == 3000 && CPU_FLAGS & CPU_AVX512F) {
		gwdone (&gwdata);
		gwinit (&gwdata);
		gwsetup (&gwdata, 1.0, 2, 10000000, -1);
	}
	if (min_test && (test < min_test || test > max_test)) continue;
	if (! (CPU_FLAGS & CPU_SSE2) && test >= 1000) break;
	if (! (CPU_FLAGS & CPU_AVX) && test >= 2000) break;
	if (! (CPU_FLAGS & CPU_AVX512F) && test >= 3000) break;
for (i = 1; i <= 50; i++) {
	start_timer (timers, 0);
	timeit (gwdata.asm_data, test, workbuf);
	end_timer (timers, 0);
	if (timers[1] == 0 || timers[1] > timers[0]) timers[1] = timers[0];
	if (i > 1 && timers[0] < 3.0 * timers[1]) {
		if (timers[0] > 1.5 * timers[1])
			i++;
		timers[2] += timers[0];
		cnt++;
	}
	timers[0] = 0;
}
sprintf (buf, "Test %d: ", test);
print_timer (timers, 1, buf, TIMER_CLR);
timers[2] /= cnt;
strcat (buf, ", avg: ");
print_timer (timers, 2, buf, TIMER_NL | TIMER_CLR);
OutputBoth (thread_num, buf);
}
aligned_free (workbuf);
gwdone (&gwdata);
if (min_test) exit (0);
return 0;
}
#endif


/* Init filename */

	tempFileName (w, filename);
	uniquifySaveFile (thread_num, filename);

/* Init the random number generator */

	srand ((unsigned) time (NULL));

/* MQ_init also requires that B is at least 120 (4 times the minimum D) */
/* Choose a default value for the second bound if none was specified */

	if (B < 120) {
		OutputStr (thread_num, "Using minimum bound #1 of 120\n");
		B = 120;
	}
	if (C == 0) C = B * 100;
	if (C <= B) C = B;

/* Set other constants */

	SQRT_B = (unsigned long) sqrt ((double) B);

/* Perform setup functions.  This includes decding how big an FFT to */
/* use, allocating memory, calling the FFT setup code, etc. */

/* Zero all data before beginning.  Init the thread number. */

	memset (&ecmdata, 0, sizeof (ecmhandle));
	ecmdata.thread_num = thread_num;

/* Setup the gwnum assembly code */

	gwinit (&ecmdata.gwdata);
	gwset_sum_inputs_checking (&ecmdata.gwdata, SUM_INPUTS_ERRCHK);
	if (IniGetInt (LOCALINI_FILE, "UseLargePages", 0)) gwset_use_large_pages (&ecmdata.gwdata);
	if (IniGetInt (INI_FILE, "HyperthreadPrefetch", 0)) gwset_hyperthread_prefetch (&ecmdata.gwdata);
	if (HYPERTHREAD_LL) {
		sp_info->normal_work_hyperthreads = IniGetInt (LOCALINI_FILE, "HyperthreadLLcount", CPU_HYPERTHREADS);
		gwset_will_hyperthread (&ecmdata.gwdata, sp_info->normal_work_hyperthreads);
	}
	gwset_bench_cores (&ecmdata.gwdata, NUM_CPUS);
	gwset_bench_workers (&ecmdata.gwdata, NUM_WORKER_THREADS);
	if (ERRCHK) gwset_will_error_check (&ecmdata.gwdata);
	gwset_num_threads (&ecmdata.gwdata, CORES_PER_TEST[thread_num] * sp_info->normal_work_hyperthreads);
	gwset_thread_callback (&ecmdata.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&ecmdata.gwdata, sp_info);
	gwset_safety_margin (&ecmdata.gwdata, IniGetFloat (INI_FILE, "ExtraSafetyMargin", 0.0));
	gwset_minimum_fftlen (&ecmdata.gwdata, w->minimum_fftlen);
	res = gwsetup (&ecmdata.gwdata, w->k, w->b, w->n, w->c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (thread_num, buf);
		return (STOP_FATAL_ERROR);
	}

/* A kludge so that the error checking code is not as strict. */

	ecmdata.gwdata.MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);

/* More random initializations */
	
	gwsetnormroutine (&ecmdata.gwdata, 0, ERRCHK, 0);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);
	first_iter_msg = TRUE;
	calc_output_frequencies (&ecmdata.gwdata, &output_frequency, &output_title_frequency);

/* Compute the number we are factoring */

	stop_reason = setN (&ecmdata.gwdata, thread_num, w, &N);
	if (stop_reason) {
		ecm_cleanup (&ecmdata);
		return (stop_reason);
	}

/* Optionally do a probable prime test */

	if (IniGetInt (INI_FILE, "ProbablePrimeTest", 0) && isProbablePrime (&ecmdata.gwdata, N)) {
		sprintf (buf, "%s is a probable prime\n", gwmodulo_as_string (&ecmdata.gwdata));
		OutputStr (thread_num, buf);
	}

/* Time various gwnum routines */

/*#define TIMING*/
#ifdef TIMING
{
unsigned long i, j, limit;
gwnum	n1, n2, n3;
n1 = gwalloc (&ecmdata.gwdata);
n2 = gwalloc (&ecmdata.gwdata);
n3 = gwalloc (&ecmdata.gwdata);
dbltogw (&ecmdata.gwdata, 283457283657.0, n2);
for (i = 1; i <= 50; i++) gwsquare (&ecmdata.gwdata, n2); /* gen random num */
gwcopy (&ecmdata.gwdata, n2, n3);
gwcopy (&ecmdata.gwdata, n2, n1);
if (n < 20000) limit = 100; else limit = 10;
for (i = 1; i <= limit; i++) {
	start_timer (timers, 0); gwsquare (&ecmdata.gwdata, n2); end_timer (0);
	start_timer (timers, 1); gwmul (&ecmdata.gwdata, n2, n3); end_timer (1);
	start_timer (timers, 2); gwfftmul (&ecmdata.gwdata, n2, n3); end_timer (2);
	start_timer (timers, 3); normalize (&ecmdata, n1, n3); end_timer (3);
	start_timer (timers, 4); gwfftfftmul (n2, n2, n2); end_timer (4);
	start_timer (timers, 5); gwfft (&ecmdata.gwdata, n2, n2); end_timer (5);
	start_timer (timers, 6); gwadd (&ecmdata.gwdata, n2, n2); end_timer (6);
	gwcopy (&ecmdata.gwdata, n1, n2);
}
strcpy (buf, "100 squares: ");
print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
strcpy (buf, "100 muls: ");
print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
strcpy (buf, "100 fftmuls: ");
print_timer (timers, 2, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
strcpy (buf, "100 ffts: ");
print_timer (timers, 5, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
strcpy (buf, "100 fftfftmuls: ");
print_timer (timers, 4, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
strcpy (buf, "100 normalizes: ");
print_timer (timers, 3, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
strcpy (buf, "100 adds: ");
print_timer (timers, 6, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
	start_timer (timers, 7);
	start_sieve (thread_num, 2, &ecmdata.sieve_info);
	for (i = 0; sieve (ecmdata.sieve_info) < 0xFFFFFFFF; i++);
	end_timer (timers, 7);
sprintf (buf, "Sieve: %ld primes found.  ", i);
print_timer (timers, 7, buf, TIMER_NL | TIMER_CLR);
OutputStr (thread_num, buf);
}		
#endif

/* Output a startup message */

	gwfft_description (&ecmdata.gwdata, fft_desc);
	sprintf (buf, "Using %s\n", fft_desc);
	OutputStr (thread_num, buf);

/* Check for a continuation file.  Limit number of backup files we try */
/* to read in case there is an error deleting bad save files. */

	readSaveFileStateInit (&read_save_file_state, thread_num, filename);
	writeSaveFileStateInit (&write_save_file_state, filename, 0);
	for ( ; ; ) {
		uint64_t save_B, save_B_processed, save_C_processed;

		if (! saveFileExists (&read_save_file_state)) {
			/* If there were save files, they are all bad.  Report a message */
			/* and temporarily abandon the work unit.  We do this in hopes that */
			/* we can successfully read one of the bad save files at a later time. */
			/* This sounds crazy, but has happened when OSes get in a funky state. */
			if (read_save_file_state.a_non_bad_save_file_existed) {
				OutputBoth (thread_num, ALLSAVEBAD_MSG);
				return (0);
			}
			/* No save files existed, start from scratch. */
			break;
		}

/* Allocate memory */

		x = gwalloc (&ecmdata.gwdata);
		if (x == NULL) goto oom;
		z = gwalloc (&ecmdata.gwdata);
		if (z == NULL) goto oom;
		gg = NULL;

/* Read in the save file.  If the save file is no good ecm_restore will have */
/* deleted it.  Loop trying to read a backup save file. */

		if (! ecm_restore (&ecmdata, thread_num, read_save_file_state.current_filename, w, &stage,
				   &curve, &sigma, &save_B, &save_B_processed,
				   &save_C_processed, x, z)) {
			gwfree (&ecmdata.gwdata, x);
			gwfree (&ecmdata.gwdata, z);
			/* Close and rename the bad save file */
			saveFileBad (&read_save_file_state);
			continue;
		}

/* Handle the case where we have a save file */
/* with a smaller bound #1 than the bound #1 we are presently working on. */
/* Restart the curve (and curve count) from scratch. */

		if (B > save_B) {
			gwfree (&ecmdata.gwdata, x);
			gwfree (&ecmdata.gwdata, z);
			break;
		}

/* Compute Ad4 from sigma */

		curve_start_msg (&ecmdata, thread_num, curve, sigma, B, C);
		t1 = gwalloc (&ecmdata.gwdata);
		if (t1 == NULL) goto oom;
		t2 = gwalloc (&ecmdata.gwdata);
		if (t2 == NULL) goto oom;
		stop_reason = choose12 (&ecmdata, w, t1, t2, sigma,
					&Ad4, N, &factor);
		if (stop_reason) goto exit;
		gwfree (&ecmdata.gwdata, t1);
		gwfree (&ecmdata.gwdata, t2);

/* Continue in the middle of stage 1 */

		if (stage == ECM_STAGE1) {
			sieve_start = save_B_processed + 1;
			goto restart1;
		}
		
/* Allocate more memory */

		gg = gwalloc (&ecmdata.gwdata);
		if (gg == NULL) goto oom;
		gwswap (z, gg);

/* We've finished stage 1, resume stage 2 */

		if (C > save_C_processed) {
			dbltogw (&ecmdata.gwdata, 1.0, z);
			stop_reason = start_sieve (thread_num, save_C_processed, &ecmdata.sieve_info);
			if (stop_reason) goto exit;
			prime = sieve (ecmdata.sieve_info);
			goto restart3;
		}
		
/* We've finished stage 2, but haven't done the GCD yet */

		goto restart4;
	}

/* Unless a save file indicates otherwise, we are testing our first curve */

	curve = 1;

/* Loop processing the requested number of ECM curves */

restart0:
	ecm_stage1_memory_usage (thread_num, &ecmdata);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Allocate memory */

	x = gwalloc (&ecmdata.gwdata);
	if (x == NULL) goto oom;
	z = gwalloc (&ecmdata.gwdata);
	if (z == NULL) goto oom;
	gg = NULL;

/* Choose curve with order divisible by 16 and choose a point (x/z) on */
/* said curve. */

	stage = 0;  /* In case we print out a factor found message! */
	do {
		uint32_t hi, lo;
		sigma = (rand () & 0x1F) * 65536.0 * 65536.0 * 65536.0;
		sigma += (rand () & 0xFFFF) * 65536.0 * 65536.0;
		if (CPU_FLAGS & CPU_RDTSC) rdtsc (&hi, &lo);
		sigma += lo ^ hi ^ ((unsigned long) rand () << 16);
	} while (sigma <= 5.0);
	if (w->curve > 5.0 && w->curve < 9007199254740992.0) sigma = w->curve;
	curve_start_msg (&ecmdata, thread_num, curve, sigma, B, C);
	stop_reason = choose12 (&ecmdata, w, x, z, sigma, &Ad4, N, &factor);
	if (stop_reason) goto exit;
	if (factor != NULL) goto bingo;
	sieve_start = 2;

/* The stage 1 restart point */

restart1:
	ecm_stage1_memory_usage (thread_num, &ecmdata);
	stage = 1;
	one_over_B = 1.0 / (double) B;
	sprintf (w->stage, "C%ldS1", curve);
	w->pct_complete = sieve_start * one_over_B;
	start_timer (timers, 0);
	stop_reason = start_sieve (thread_num, sieve_start, &ecmdata.sieve_info);
	if (stop_reason) goto exit;
	for ( ; ; ) {
		prime = sieve (ecmdata.sieve_info);
		if (prime > B) break;

/* Apply as many powers of prime as long as prime^n <= B */
/* MEMUSED: 3 gwnums (x, z, AD4) + 10 for ell_mul */

		stop_reason = ell_mul (&ecmdata, x, z, prime, Ad4);
		if (stop_reason) goto exit;
		if (prime <= SQRT_B) {
			uint64_t mult, max;
			mult = prime;
			max = B / prime;
			for ( ; ; ) {
				stop_reason = ell_mul (&ecmdata, x, z, prime, Ad4);
				if (stop_reason) goto exit;
				mult *= prime;
				if (mult > max) break;
			}
		}

/* Calculate stage 1 percent complete */

		w->pct_complete = prime * one_over_B;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&ecmdata.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s ECM curve %d stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&ecmdata.gwdata), (int) curve);
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&ecmdata.gwdata);
		}

/* Print a message every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&ecmdata.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s curve %d stage 1 at prime %.0f [%.*f%%].",
				 gwmodulo_as_string (&ecmdata.gwdata), (int) curve, (double) prime,
				 (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			if (prime != 2)
				OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&ecmdata.gwdata);
			first_iter_msg = FALSE;
		}

/* Check for errors */

		if (gw_test_for_error (&ecmdata.gwdata)) goto err;

/* Write a save file when the user interrupts the calculation and */
/* every DISK_WRITE_TIME minutes. */

		stop_reason = stopCheck (thread_num);
		if (stop_reason || testSaveFilesFlag (thread_num)) {
			ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE1, curve, sigma, B, prime, 0, x, z);
			if (stop_reason) goto exit;
		}
	}

/* Stage 1 complete */

	end_timer (timers, 0);
	sprintf (buf, "Stage 1 complete. %.0f transforms, %lu modular inverses. Time: ",
		 gw_get_fft_count (&ecmdata.gwdata), ecmdata.modinv_count);
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&ecmdata.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&ecmdata.gwdata);
	}

/* If we aren't doing a stage 2, then check to see if we found a factor. */
/* If we are doing a stage 2, then the stage 2 init will do this GCD for us. */

	if (C <= B) {
skip_stage_2:	start_timer (timers, 0);
		stop_reason = gcd (&ecmdata.gwdata, thread_num, z, N, &factor);
		if (stop_reason) {
			ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE1, curve, sigma, B, B, 0, x, z);
			goto exit;
		}
		end_timer (timers, 0);
		strcpy (buf, "Stage 1 GCD complete. Time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputStr (thread_num, buf);
		if (factor != NULL) goto bingo;

/* Alexander Kruppa wrote this code to normalize and output the x value */
/* along with N and sigma so that it can be used in Paul Zimmermann's */
/* superior GMP-ECM implementation of stage 2. */

		if (IniGetInt (INI_FILE, "GmpEcmHook", 0)) {
			char	*msg, *buf;
			int	msglen, i, leadingzeroes;
			giant	gx;
			char	*hex = "0123456789ABCDEF";

			stop_reason = normalize (&ecmdata, x, z, N, &factor);
			if (stop_reason) goto exit;

			gx = popg (&ecmdata.gwdata.gdata, ((int) ecmdata.gwdata.bit_length >> 5) + 10);
			if (gx == NULL) goto oom;
			if (gwtogiant (&ecmdata.gwdata, x, gx)) goto oom;  // Unexpected error, return oom for lack of a better error message
			modgi (&ecmdata.gwdata.gdata, N, gx);

			msglen = N->sign * 8 + 5;
			buf = (char *) malloc (msglen + msglen + 80);
			if (buf == NULL) goto oom;
			strcpy (buf, "N=0x");
			msg = buf + strlen (buf);
			leadingzeroes = 1; /* Still eat leading zeroes? */
			for (i = 0; i < N->sign * 8; i++) {
				char nibble = ( N->n[N->sign - 1 - i/8] >> ((7-i%8)*4)) & 0xF;
				if (nibble != 0) leadingzeroes = 0;
				if (!leadingzeroes) *msg++ = hex[nibble];
			}
			strcpy (msg, "; QX=0x");
			msg = msg + strlen (msg);
			leadingzeroes = 1;
			for (i = 0; i < gx->sign * 8; i++) {
				char nibble = ( gx->n[gx->sign - 1 - i/8] >> ((7-i%8)*4)) & 0xF;
				if (nibble != 0) leadingzeroes = 0;
				if (!leadingzeroes) *msg++ = hex[nibble];
			}
			strcpy (msg, "; SIGMA=");
			msg = msg + strlen (msg);
			sprintf (msg, "%.0f\n", sigma);
			writeResults (buf);
			free (buf);
			pushg (&ecmdata.gwdata.gdata, 1);
		}

/* Now do the next ECM curve */

		goto more_curves;
	}

/*
   Stage 2:  We support two types of stage 2's here.  One uses
   less memory and uses fewer extended GCDs, but is slower in accumulating
   each found prime.  Thanks to Richard Crandall and Paul Zimmermann
   for letting me liberally use their code and ideas here.
   x, z: coordinates of Q at the beginning of stage 2
*/

/* Make sure we will have enough memory to run stage 2 at some time */
/* We need at least 20 gwnums. */

restart3:
	min_memory = cvt_gwnums_to_mem (&ecmdata.gwdata, 20);
	if (max_mem (thread_num) < min_memory) {
		sprintf (buf, "Skipping stage 2 due to insufficient memory -- %ldMB needed.\n", min_memory);
		OutputStr (thread_num, buf);
		C = B;
		goto skip_stage_2;
	}

/* Initialize variables for second stage */
/* Our goal is to fill up the nQx array with Q^1, Q^3, Q^5, ... */
/* normalized with only one modular inverse call. */

	start_timer (timers, 0);
	sprintf (w->stage, "C%ldS2", curve);
	one_over_C_minus_B = 1.0 / (double) (C - B);
	w->pct_complete = 0.0;

/* Choose a good value for D.  One that reduces the number of */
/* multiplications, yet doesn't use too much memory. */

	stop_reason = choose_stage2_plan (thread_num, &ecmdata, B, C);
	if (stop_reason) {
		if (gg == NULL) {
			ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE1, curve, sigma, B, B, 0, x, z);
		}
		goto exit;
	}

/* Allocate more memory.  D/3 is enough for the nQx values and */
/* we need an additional D/3 or E*2 values for pooling in case we */
/* are using the POOL_3MULT algorithm */

	{
		unsigned long i, max;
		max = (ecmdata.D/3 > ecmdata.E*2) ? ecmdata.D/3 : ecmdata.E*2;
		gw_set_max_allocs (&ecmdata.gwdata, ecmdata.D/3 + max + 20);
		ecmdata.nQx = (gwnum *) malloc ((ecmdata.D/2) * sizeof (gwnum));
		if (ecmdata.nQx == NULL) goto oom;
		for (i = 0; i < ecmdata.D/2; i++) ecmdata.nQx[i] = NULL;
		ecmdata.pool_values = (gwnum *) malloc (max * sizeof (gwnum));
		if (ecmdata.pool_values == NULL) goto oom;
		ecmdata.poolz_values = (gwnum *) malloc (max * sizeof (gwnum));
		if (ecmdata.poolz_values == NULL) goto oom;
		ecmdata.mQx = (gwnum *) malloc (ecmdata.E * sizeof (gwnum));
		if (ecmdata.mQx == NULL) goto oom;
		ecmdata.pairings = (char *) malloc ((ecmdata.D + 15) >> 4);
		if (ecmdata.pairings == NULL) goto oom;
	}

/* Allocate memory for computing nQx values */
/* MEMUSED: 9 gwnums (x, z, AD4, 6 for nQx) */

	Q2x = gwalloc (&ecmdata.gwdata);
	if (Q2x == NULL) goto oom;
	Q2z = gwalloc (&ecmdata.gwdata);
	if (Q2z == NULL) goto oom;
	Qiminus2x = gwalloc (&ecmdata.gwdata);
	if (Qiminus2x == NULL) goto oom;
	Qiminus2z = gwalloc (&ecmdata.gwdata);
	if (Qiminus2z == NULL) goto oom;
	Qdiffx = gwalloc (&ecmdata.gwdata); Qdiffz = gwalloc (&ecmdata.gwdata);
	if (Qdiffx == NULL) goto oom;

/* Init values used in computing nQx.  We need Q^2, Q^1, and diff of Q^1. */
/* MEMUSED: 9 gwnums (x, z, AD4, 6 for computing nQx) + 2 temporaries */

	stop_reason = ell_dbl (&ecmdata, x, z, Q2x, Q2z, Ad4);
	if (stop_reason) goto exit;
	ell_begin_fft (&ecmdata, Q2x, Q2z, Q2x, Q2z);
	gwfft (&ecmdata.gwdata, x, x); gwfft (&ecmdata.gwdata, z, z);
	gwcopy (&ecmdata.gwdata, x, Qdiffx);
	gwcopy (&ecmdata.gwdata, z, Qdiffz);
	gwcopy (&ecmdata.gwdata, x, Qiminus2x);
	gwcopy (&ecmdata.gwdata, z, Qiminus2z);

/* Init the first nQx value with Q^1 */
/* MEMUSED: 9 gwnums (AD4, 6 for computing nQx, nQx[0], modinv_value) */

	ecmdata.nQx[0] = x;
	stop_reason = add_to_normalize_pool (&ecmdata, ecmdata.nQx[0], z, 1);
	if (stop_reason) goto exit;
	gwfree (&ecmdata.gwdata, z);

/* Compute the rest of the nQx values (Q^i for i >= 3) */
/* MEMUSED: 8 + nQx gwnums (AD4, 6 for computing nQx, nQx vals, modinv_val) */
/* MEMPEAK: 8 + nQx-1 + 2 for ell_add temporaries */

	for (i = 3; i < ecmdata.D; i = i + 2) {
		ell_add_special3 (&ecmdata, Qiminus2x, Qiminus2z, Q2x, Q2z, &Qdiffx, &Qdiffz);

		if (gw_test_for_error (&ecmdata.gwdata)) goto err;

		stop_reason = stopCheck (thread_num);
		if (stop_reason) {
			if (gg == NULL) {
				ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE1,
					  curve, sigma, B, B, 0, Qdiffx, Qdiffz);
			}
			goto exit;
		}

		gwfft (&ecmdata.gwdata, Qdiffx, Qdiffx);
		gwfft (&ecmdata.gwdata, Qdiffz, Qdiffz);
		if (relatively_prime (i, ecmdata.D)) {
			j = (i - 1) >> 1;
			ecmdata.nQx[j] = gwalloc (&ecmdata.gwdata);
			if (ecmdata.nQx[j] == NULL) goto oom;
			gwcopy (&ecmdata.gwdata, Qdiffx, ecmdata.nQx[j]);
			stop_reason = add_to_normalize_pool (&ecmdata, ecmdata.nQx[j], Qdiffz, 1);
			if (stop_reason) goto exit;
		}
		gwswap (Qdiffx, Qiminus2x); gwswap (Qdiffz, Qiminus2z);
	}

/* Now compute Q^2D.  This will be used in computing Q^m values. */
/* Qiminus2 is Q^(D-1) and Qdiff if Q^(D-3).  Add Q^(D-1) and Q^2 to */
/* get Q^(D+1).  Then add Q^(D-1) and Q^(D+1) to get Q^2D.  Store Q^2D */
/* in Q2x and Q2z.  Normalize them so we can free Q2z later on. */
/* MEMUSED: 8 + nQx gwnums (AD4, 6 for computing nQx, nQx vals, modinv_val) */
/* MEMPEAK: 8 + nQx + 2 for ell_add temporaries */

	ell_add_special3 (&ecmdata, Qiminus2x, Qiminus2z, Q2x, Q2z, &Qdiffx, &Qdiffz);
	gwfftaddsub (&ecmdata.gwdata, Q2x, Q2z); /* Recompute fft of Q2x,Q2z */
	ell_begin_fft (&ecmdata, Qdiffx, Qdiffz, Qdiffx, Qdiffz);
	ell_add_special3 (&ecmdata, Qiminus2x, Qiminus2z, Qdiffx, Qdiffz, &Q2x, &Q2z);
	gwfft (&ecmdata.gwdata, Q2x, Q2x); gwfft (&ecmdata.gwdata, Q2z, Q2z);
	stop_reason = add_to_normalize_pool (&ecmdata, Q2x, Q2z, 1);
	if (stop_reason) goto exit;

/* Free most of the memory used in computing nQx values */
/* Keep two values we could free in case the upcoming normalize is */
/* aborted and we need to write a save file. */
/* MEMUSED: 5 + nQx gwnums (AD4, Q2x, Qiminus2x&z, nQx values, modinv_value) */

	gwfree (&ecmdata.gwdata, Q2z);
	gwfree (&ecmdata.gwdata, Qdiffx); gwfree (&ecmdata.gwdata, Qdiffz);

/* Normalize all the nQx values */
/* MEMUSED: 4 + nQx gwnums (AD4, Q2x, Qiminus2x&z, nQx values) */
/* MEMPEAK: 4 + nQx + 3 for UV, normalize, and pooled_modinv temporaries */

	stop_reason = normalize_pool (&ecmdata, N, &factor);
	if (stop_reason) {
		if (gg == NULL) {
			dbltogw (&ecmdata.gwdata, 1.0, Q2x);
			gwfft (&ecmdata.gwdata, Q2x, Q2x);
			gwfftfftmul (&ecmdata.gwdata, Q2x, Qiminus2x, Qiminus2x);
			gwfftfftmul (&ecmdata.gwdata, Q2x, Qiminus2z, Qiminus2z);
			ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE1, curve,
				  sigma, B, B, 0, Qiminus2x, Qiminus2z);
		}
		goto exit;
	}

/* If we found a factor, we're done */

	if (factor != NULL) goto bingo;

/* Free rest of the memory used in computing nQx values */
/* MEMUSED: 2 + nQx gwnums (AD4, Q2x, nQx values) */

	gwfree (&ecmdata.gwdata, Qiminus2x);
	gwfree (&ecmdata.gwdata, Qiminus2z);

/* Init code that computes Q^m */
/* MEMUSED: 6 + nQx gwnums (6 for computing mQx, nQx values) */
/* MEMPEAK: 6 + nQx + 6 for bin_ell_mul temporaries */

	m = (prime / ecmdata.D + 1) * ecmdata.D;
	stop_reason = mQ_init (&ecmdata, ecmdata.nQx[0], m, Q2x, Ad4);
	if (stop_reason) goto exit;

/* Precompute the transforms of nQx */

	for (i = 0; i < ecmdata.D/2; i++)
		if (ecmdata.nQx[i] != NULL)
			gwfft (&ecmdata.gwdata, ecmdata.nQx[i], ecmdata.nQx[i]);

/* Now init the accumulator unless this value was read */
/* from a continuation file */
/* MEMUSED: 7 + nQx gwnums (6 for computing mQx, gg, nQx values) */

	if (gg == NULL) {
		gg = gwalloc (&ecmdata.gwdata);
		if (gg == NULL) goto oom;
		dbltogw (&ecmdata.gwdata, 1.0, gg);
	}

/* Initialization of stage 2 complete */

	end_timer (timers, 0);
	sprintf (buf, "Stage 2 init complete. %.0f transforms, %lu modular inverses. Time: ",
		 gw_get_fft_count (&ecmdata.gwdata), ecmdata.modinv_count);
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Now do stage 2 */
/* Accumulate (mQx - nQx)(mQz + nQz) - mQx mQz + nQx nQz.		*/
/* Since nQz = 1, we have (the 4 FFT per prime continuation)		*/
/*		== (mQx - nQx)(mQz + 1) - mQx mQz + nQx			*/
/*		== mQx mQz - nQx mQz + mQx - nQx - mQx mQz + nQx	*/
/*		== mQx - nQx mQz					*/
/* If mQz also = 1 (the 2 FFT per prime continuation) then we accumulate*/
/*		== mQx - nQx						*/

	start_timer (timers, 0);
	stage = 2;
	for ( ; C > m-ecmdata.D; m += ecmdata.D+ecmdata.D) {
		gwnum	mQx, mQz;

/* Compute next Q^m value */
/* MEMUSED: 7 + nQx + E gwnums (6 for computing mQx, gg, nQx and E values) */
/* MEMPEAK: 7 + nQx + E + 2 for ell_add temporaries */

		stop_reason = mQ_next (&ecmdata, &mQx, &mQz, N, &factor);
		if (stop_reason) {
			// In case stop_reason is out-of-memory, free some up
			// before calling ecm_save.
			mQ_term (&ecmdata);
			t1 = gwalloc (&ecmdata.gwdata);
			if (t1 == NULL) goto oom;
			dbltogw (&ecmdata.gwdata, 1.0, t1);
			gwmul (&ecmdata.gwdata, t1, gg);
			gwfftfftmul (&ecmdata.gwdata, t1, ecmdata.nQx[0], ecmdata.nQx[0]);
			ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE2, curve,
				  sigma, B, B, prime, ecmdata.nQx[0], gg);
			gwfree (&ecmdata.gwdata, t1);
			goto exit;
		}
		if (factor != NULL) goto bingo;
		memset (ecmdata.pairings, 0, (ecmdata.D + 15) >> 4);
		t1 = gwalloc (&ecmdata.gwdata);
		if (t1 == NULL) goto oom;

/* 2 FFT per prime continuation - deals with all normalized values */

		if (ecmdata.TWO_FFT_STAGE2) {
		    for ( ; ; prime = sieve (ecmdata.sieve_info)) {
			if (prime < m) {	/* Do the m-D to m range */
				i = (unsigned long) (m - prime) >> 1;
				bitset (ecmdata.pairings, i);
			} else if (prime < m+ecmdata.D) { /* Do the m to m+D range */
				i = (unsigned long) (prime - m) >> 1;
				if (bittst (ecmdata.pairings, i)) continue;
			} else
				break;
			gwfftsub3 (&ecmdata.gwdata, mQx, ecmdata.nQx[i], t1);
			gwstartnextfft (&ecmdata.gwdata, TRUE);
			gwfftmul (&ecmdata.gwdata, t1, gg);
			gwstartnextfft (&ecmdata.gwdata, FALSE);
		    }
		}

/* 4 FFT per prime continuation - deals with only nQx values normalized */

		else {
		    for ( ; ; prime = sieve (ecmdata.sieve_info)) {
			if (prime < m) {	/* Do the m-D to m range */
				i = (unsigned long) (m - prime) >> 1;
				bitset (ecmdata.pairings, i);
			} else if (prime < m+ecmdata.D) { /* Do the m to m+D range */
				i = (unsigned long) (prime - m) >> 1;
				if (bittst (ecmdata.pairings, i)) continue;
			} else
				break;
			gwstartnextfft (&ecmdata.gwdata, TRUE);
			gwfftfftmul (&ecmdata.gwdata, ecmdata.nQx[i], mQz, t1);
			gwstartnextfft (&ecmdata.gwdata, FALSE);
			gwfft (&ecmdata.gwdata, t1, t1);
			gwfftsub3 (&ecmdata.gwdata, mQx, t1, t1);
			gwstartnextfft (&ecmdata.gwdata, TRUE);
			gwfftmul (&ecmdata.gwdata, t1, gg);
			gwstartnextfft (&ecmdata.gwdata, FALSE);
		    }
		}
		gwfree (&ecmdata.gwdata, t1);

/* Calculate stage 2 percent complete */

		w->pct_complete = (prime - B) * one_over_C_minus_B;
		if (w->pct_complete > 1.0) w->pct_complete = 1.0;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&ecmdata.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s ECM curve %d stage 2",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&ecmdata.gwdata), (int) curve);
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&ecmdata.gwdata);
		}

/* Print a message every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&ecmdata.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s curve %d stage 2 at prime %.0f [%.*f%%].",
				 gwmodulo_as_string (&ecmdata.gwdata), (int) curve, (double) prime,
				 (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&ecmdata.gwdata);
			first_iter_msg = FALSE;
		}

/* Check for errors */

		if (gw_test_for_error (&ecmdata.gwdata)) goto err;

/* Write a save file when the user interrupts the calculation and */
/* every DISK_WRITE_TIME minutes. */

		stop_reason = stopCheck (thread_num);
		if (stop_reason || testSaveFilesFlag (thread_num)) {
			t1 = gwalloc (&ecmdata.gwdata);
			if (t1 == NULL) goto oom;
			dbltogw (&ecmdata.gwdata, 1.0, t1);
			gwmul (&ecmdata.gwdata, t1, gg);
			gwfftfftmul (&ecmdata.gwdata, t1, ecmdata.nQx[0], t1);
			ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE2, curve,
				  sigma, B, B, prime, t1, gg);
			gwfree (&ecmdata.gwdata, t1);
			if (stop_reason) goto exit;
		}
	}
	mQ_term (&ecmdata);
	t1 = gwalloc (&ecmdata.gwdata);
	if (t1 == NULL) goto oom;
	dbltogw (&ecmdata.gwdata, 1.0, t1);
	gwmul (&ecmdata.gwdata, t1, gg);
	gwfree (&ecmdata.gwdata, t1);

/* Stage 2 is complete */

	end_timer (timers, 0);
	sprintf (buf, "Stage 2 complete. %.0f transforms, %lu modular inverses. Time: ",
		 gw_get_fft_count (&ecmdata.gwdata), ecmdata.modinv_count);
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	last_output = last_output_t = ecmdata.modinv_count = 0;
	gw_clear_fft_count (&ecmdata.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&ecmdata.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&ecmdata.gwdata);
	}

/* See if we got lucky! */

restart4:
	ecm_stage1_memory_usage (thread_num, &ecmdata);
	sprintf (w->stage, "C%ldS2", curve);
	w->pct_complete = 1.0;
	start_timer (timers, 0);
	stop_reason = gcd (&ecmdata.gwdata, thread_num, gg, N, &factor);
	if (stop_reason) {
		ecm_save (&ecmdata, &write_save_file_state, w, ECM_STAGE2, curve, sigma, B, B, C, gg, gg);
		goto exit;
	}
	end_timer (timers, 0);
	strcpy (buf, "Stage 2 GCD complete. Time: ");
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	if (factor != NULL) goto bingo;

/* Do not loop if we are testing a specific curve */

more_curves:
	ecm_partial_cleanup (&ecmdata);
	if (w->curve < 5.0 && ++curve <= w->curves_to_do)
		goto restart0;

/* Output line to results file indicating the number of curves run */

	sprintf (buf, "%s completed %u ECM %s, B1=%.0f, B2=%.0f, Wh%d: %08lX\n",
		 gwmodulo_as_string (&ecmdata.gwdata), w->curves_to_do,
		 w->curves_to_do == 1 ? "curve" : "curves",
		 (double) B, (double) C, PORT, SEC5 (w->n, B, C));
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);
	writeResults (buf);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"NF", "exponent":45581713, "worktype":"ECM", "b1":50000, "b2":5000000, */
/* "curves":5, "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"NF\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"ECM\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%.0f, \"b2\":%.0f", (double) B, (double) C);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"curves\":%u", w->curves_to_do);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", ecmdata.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, B, C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send ECM completed message to the server.  Although don't do it for puny B1 values. */

	if (B >= 10000 || IniGetInt (INI_FILE, "SendAllFactorData", 0)) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		truncated_strcpy (pkt.message, sizeof (pkt.message), buf);
		pkt.result_type = PRIMENET_AR_ECM_NOFACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		pkt.B1 = (double) B;
		pkt.B2 = (double) C;
		pkt.curves = w->curves_to_do;
		pkt.fftlen = gwfftlen (&ecmdata.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Delete the save file */

	unlinkSaveFiles (&write_save_file_state);

/* Free memory and return */

	stop_reason = STOP_WORK_UNIT_COMPLETE;
exit:	ecm_cleanup (&ecmdata);
	free (N);
	free (factor);
	free (str);
	free (msg);
	return (stop_reason);

/* We've run out of memory.  Print error message and exit. */

oom:	stop_reason = OutOfMemory (thread_num);
	goto exit;

/* Print a message, we found a factor! */

bingo:	sprintf (buf, "ECM found a factor in curve #%ld, stage #%d\n", curve, stage);
	writeResults (buf);
	sprintf (buf, "Sigma=%.0f, B1=%.0f, B2=%.0f.\n", sigma, (double) B, (double) C);
	writeResults (buf);

/* Allocate memory for the string representation of the factor and for */
/* a message.  Convert the factor to a string. */ 

	msglen = factor->sign * 10 + 400;
	str = (char *) malloc (msglen);
	if (str == NULL) goto oom;
	msg = (char *) malloc (msglen);
	if (msg == NULL) goto oom;
	gtoc (factor, str, msglen);

/* Validate the factor we just found */

	if (!testFactor (&ecmdata.gwdata, w, factor)) {
		sprintf (msg, "ERROR: Bad factor for %s found: %s\n", gwmodulo_as_string (&ecmdata.gwdata), str);
		OutputBoth (thread_num, msg);
		OutputStr (thread_num, "Restarting ECM curve from scratch.\n");
		continueECM = TRUE;
		curve--;
		goto bad_factor_recovery;
	}

/* Output the validated factor */

	sprintf (msg, "%s has a factor: %s (ECM curve %d, B1=%.0f, B2=%.0f)\n",
		 gwmodulo_as_string (&ecmdata.gwdata), str, (int) curve, (double) B, (double) C);
	OutputStr (thread_num, msg);
	formatMsgForResultsFile (msg, w);
	writeResults (msg);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"F", "exponent":45581713, "worktype":"ECM", "factors":["430639100587696027847"], */
/* "b1":50000, "b2":5000000, "sigma":"123456789123456", "stage":2 */
/* "curves":5, "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"F\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"ECM\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"factors\":[\"%s\"]", str);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%.0f, \"b2\":%.0f", (double) B, (double) C);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"sigma\":%.0f", sigma);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"stage\":%d", stage);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"curves\":%lu", curve);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", ecmdata.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, B, C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* See if the cofactor is prime and set flag if we will be continuing ECM */

	continueECM = IniGetInt (INI_FILE, "ContinueECM", 0);
	prpAfterEcmFactor = IniGetInt (INI_FILE, "PRPAfterECMFactor", bitlen (N) < 100000);
	if (prpAfterEcmFactor || continueECM) divg (factor, N);
	if (prpAfterEcmFactor && isProbablePrime (&ecmdata.gwdata, N)) {
		OutputBoth (thread_num, "Cofactor is a probable prime!\n");
		continueECM = FALSE;
	}

/* Send assignment result to the server.  To avoid flooding the server */
/* with small factors from users needlessly redoing factoring work, make */
/* sure the factor is more than 50 bits or so. */

	if (strlen (str) >= 15 || IniGetInt (INI_FILE, "SendAllFactorData", 0)) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		truncated_strcpy (pkt.message, sizeof (pkt.message), msg);
		pkt.result_type = PRIMENET_AR_ECM_FACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		truncated_strcpy (pkt.factor, sizeof (pkt.factor), str);
		pkt.B1 = (double) B;
		pkt.B2 = (double) C;
		pkt.curves = curve;
		pkt.stage = stage;
		pkt.fftlen = gwfftlen (&ecmdata.gwdata);
		pkt.done = !continueECM;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);

/* If continuing ECM, subtract the curves we just reported from the */
/* worktodo count of curves to run.  Otherwise, delete all ECM entries */
/* for this number from the worktodo file. */

		if (continueECM) {
			unlinkSaveFiles (&write_save_file_state);
			w->curves_to_do -= curve;
			stop_reason = updateWorkToDoLine (thread_num, w);
			if (stop_reason) return (stop_reason);
			curve = 0;
		} else {
//bug - how to update worktodo such that all ECM's of this number are deleted???
		}
	}

/* Free memory */

bad_factor_recovery:
	free (str); str = NULL;
	free (msg); msg = NULL;
	free (factor); factor = NULL;

	clear_timer (timers, 0);

/* Since we found a factor, then we likely performed much fewer curves than */
/* expected.  Make sure we do not update the rolling average with */
/* this inaccurate data. */

	if (!continueECM) {
		unlinkSaveFiles (&write_save_file_state);
		stop_reason = STOP_WORK_UNIT_COMPLETE;
		invalidateNextRollingAverageUpdate ();
		goto exit;
	}

/* Do more curves despite finding a factor */

	goto more_curves;

/* Output a message saying we are restarting */

err:	OutputBoth (thread_num, "SUMOUT error occurred.\n");

/* Sleep five minutes before restarting */

	ecm_cleanup (&ecmdata);
	free (N); N = NULL;
	stop_reason = SleepFive (thread_num);
	if (stop_reason) return (stop_reason);

/* Restart from last save file */

	goto restart;
}

/* Read a file of ECM tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, sigma, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9991 to run the QA suite */

int ecm_QA (
	int	thread_num,
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	FILE	*fd;

/* Set the title */

	title (thread_num, "QA");

/* Open QA file */

	fd = fopen ("qa_ecm", "r");
	if (fd == NULL) {
		OutputStr (thread_num, "File named 'qa_ecm' could not be opened.\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop until the entire file is processed */

	QA_TYPE = 0;
	for ( ; ; ) {
		struct work_unit w;
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		double	sigma;
		int	stop_reason;

/* Read a line from the file */

		n = 0;
		(void) fscanf (fd, "%lf,%lu,%lu,%ld,%lf,%lu,%lu,%lu,%s\n",
			&k, &b, &n, &c, &sigma, &B1, &B2_start, &B2_end,
			fac_str);
		if (n == 0) break;

/* If b is 1, set QA_TYPE */

		if (b == 1) {
			QA_TYPE = c;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		QA_FACTOR = allocgiant ((int) strlen (fac_str));
		ctog (fac_str, QA_FACTOR);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Do the ECM */

		if (B2_start < B1) B2_start = B1;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_ECM;
		w.k = k;
		w.b = b;
		w.n = n;
		w.c = c;
		w.B1 = B1;
		w.B2_start = B2_start;
		w.B2 = B2_end;
		w.curves_to_do = 1;
		w.curve = sigma;
		QA_IN_PROGRESS = TRUE;
		stop_reason = ecm (0, sp_info, &w);
		QA_IN_PROGRESS = FALSE;
		free (QA_FACTOR);
		if (stop_reason != STOP_WORK_UNIT_COMPLETE) {
			fclose (fd);
			return (stop_reason);
		}
	}

/* Cleanup */

	fclose (fd);
	return (0);
}

/**************************************************************
 *
 *	P-1 Functions
 *
 **************************************************************/

/* Data maintained during P-1 process */

#define PM1_STAGE0	3	/* In stage 1, squaring small primes */
#define PM1_STAGE1	0	/* In stage 1, processing larger primes */
#define PM1_STAGE2	1	/* In stage 2 */
#define PM1_DONE	2	/* P-1 job complete */

typedef struct {
	gwhandle gwdata;	/* GWNUM handle */
	int	thread_num;	/* Worker thread number */
	unsigned long stage;	/* One of the 4 states listed above */
	unsigned long D;	/* Stage 2 loop size */
	unsigned long E;	/* Suyama's power in stage 2 */
	gwnum	*nQx;		/* Array of data used in stage 2 */
	gwnum	*eQx;		/* Array of data used in stage 2 of P-1 */
	void	*sieve_info;	/* Prime number sieve */
	unsigned long max_stage_0_prime; /* Maximum small prime that can be used in stage 0 exponent calculation. */
	uint64_t B_done;	/* We have completed calculating 3^e to this bound #1 */
	uint64_t B;		/* We are trying to increase bound #1 to this value */
	uint64_t C_done;	/* Bound #2 has been computed to this value */
	uint64_t C_start;	/* We are trying to increase bound #2 */
				/* from this starting point.  This is the */
				/* same as C_done except when using the */
				/* untested worktodo feature that allows */
				/* doing part of stage 2 on several machines */
	uint64_t C;		/* We are advancing bound #2 to this */
				/* value using the bit array.  Large B2 */
				/* values require us to break the bit */
				/* array into more than one chunk */
	unsigned long numrels;	/* Number of values relatively prime to D */
	unsigned long rels_done;/* In multi-pass processing of the bit */
				/* array, the number of relative primes */
				/* already processed */
	unsigned long rels_this_pass;
				/* In multi-pass processing of the bit */
				/* array, the number of relative primes we */
				/* are processing this pass. */
	char	*bitarray;	/* Bit array for primes between */
				/* bitarray_first_number and C */
	unsigned long bitarray_len;
				/* Number of bytes in the bit array */
	uint64_t bitarray_first_number; /* The number corresponding to the */
				/* first bit in the bit array. */
	unsigned long pairs_set;/* Number of pairs originally set in */
				/* bitarray */
	unsigned long pairs_done;/* Number of pairs completed */
	double	pct_mem_to_use;	/* If we get memory allocation errors, we */
				/* progressively try using less and less. */
} pm1handle;

/* Perform cleanup functions. */

void pm1_cleanup (
	pm1handle *pm1data)
{

/* Free memory */

	free (pm1data->nQx);
	free (pm1data->eQx);
	free (pm1data->bitarray);
	gwdone (&pm1data->gwdata);
	end_sieve (pm1data->sieve_info);
	memset (pm1data, 0, sizeof (pm1handle));
}

/* Raises number to the given power */

void pm1_mul (
	pm1handle *pm1data,
	gwnum	xx,
	gwnum	orig_xx_fft,
	uint64_t n)
{
	uint64_t c;

/* Find most significant bit and then ignore it */

	c = 1;
	c <<= 63;
	while ((c&n) == 0) c >>= 1;
	c >>= 1;

/* Handle the second most significant bit */

	gwstartnextfft (&pm1data->gwdata, c > 1);
	gwfft (&pm1data->gwdata, xx, orig_xx_fft);
	gwfftfftmul (&pm1data->gwdata, orig_xx_fft, orig_xx_fft, xx);
	if (c&n) gwfftmul (&pm1data->gwdata, orig_xx_fft, xx);
	c >>= 1;

/* Do the rest of the bits */

	while (c) {
		gwstartnextfft (&pm1data->gwdata, c > 1);
		gwsquare (&pm1data->gwdata, xx);
		if (c&n) gwfftmul (&pm1data->gwdata, orig_xx_fft, xx);
		c >>= 1;
	}
}

/* Code to init "finite differences" for computing successive */
/* values of x^(start+i*incr)^E */

int fd_init (
	pm1handle *pm1data,
	uint64_t start,
	unsigned long incr,
	gwnum	x)		/* Caller must pass in the FFT of x */
{
	unsigned long i, j;
	giant	p;

/* Treat each eQx[i] as a binary value and compute (start+i*incr)^e */

	for (i = 0; i <= pm1data->E; i++) {
		uint64_t val;
		p = allocgiant (pm1data->E * 2);
		if (p == NULL) goto oom;
		val = start + i * incr;
		ulltog (val, p);
		for (j = 2; j <= pm1data->E; j++) ullmulg (val, p);
		pm1data->eQx[i] = (gwnum) p;
	}		

/* Now do the finite differences */

	for (i = 1; i <= pm1data->E; i++) {
		for (j = pm1data->E; j >= i; j--) {
			subg ((giant) pm1data->eQx[j-1],
			      (giant) pm1data->eQx[j]);
		}
	}

/* Now compute each x^difference */

	for (i = 0; i <= pm1data->E; i++) {
		p = (giant) pm1data->eQx[i];
		pm1data->eQx[i] = gwalloc (&pm1data->gwdata);
		if (pm1data->eQx[i] == NULL) goto oom;

/* Test for easy cases */

		ASSERTG (!isZero (p));
		if (isone (p)) {
			gwcopy (&pm1data->gwdata, x, pm1data->eQx[i]);
		}

/* Find most significant bit and then ignore it */

		else {
			int	len;

			len = bitlen (p);
			len--;

/* Perform the first squaring using the already FFTed value of x */
/* Then process the second and remaining bits of p */

#ifndef SERVER_TESTING
			gwfftfftmul (&pm1data->gwdata, x, x, pm1data->eQx[i]);
			for ( ; ; ) {
				if (bitval (p, len-1))
					gwfftmul (&pm1data->gwdata, x, pm1data->eQx[i]);
				len--;
				if (len == 0) break;
				gwsquare (&pm1data->gwdata, pm1data->eQx[i]);
			}

/* FFT the final result */

			gwfft (&pm1data->gwdata, pm1data->eQx[i], pm1data->eQx[i]);
#endif
		}
		free (p);
	}
	return (0);

/* Out of memory exit path */

oom:	free (p);
	return (OutOfMemory (pm1data->thread_num));
}

/* Code to compute next x^(start+i*incr)^E value */
/* Value is returned in eQx[0] - already FFTed */

void fd_next (
	pm1handle *pm1data)
{
	unsigned long i;

	for (i = 0; i < pm1data->E; i++) {
#ifndef SERVER_TESTING
		gwfftfftmul (&pm1data->gwdata, pm1data->eQx[i], pm1data->eQx[i+1], pm1data->eQx[i]);
		gwfft (&pm1data->gwdata, pm1data->eQx[i], pm1data->eQx[i]);
#endif
	}
}

/* Terminate finite differences code */

void fd_term (
	pm1handle *pm1data)
{
	unsigned long i;

/* Free each eQx[i] */

	for (i = 0; i <= pm1data->E; i++)
		gwfree (&pm1data->gwdata, pm1data->eQx[i]);
}

/* Routines to create and read save files for a P-1 factoring job */

#define PM1_MAGICNUM	0x317a394b
//#define PM1_VERSION	2				/* Changed in 29.4 build 7 -- corrected calc_exp bug */
#define PM1_VERSION	3				/* Changed in 29.8 build 8.  Configurable calc_exp max exponent */

void pm1_save (
	pm1handle *pm1data,
	writeSaveFileState *write_save_file_state,
	struct work_unit *w,
	uint64_t processed,
	gwnum	x,
	gwnum	gg)
{
	int	fd;
	unsigned long sum = 0;

/* Create the intermediate file */

	fd = openWriteSaveFile (write_save_file_state);
	if (fd < 0) return;

/* Write the file header */

	if (!write_header (fd, PM1_MAGICNUM, PM1_VERSION, w)) goto writeerr;

/* Write the file data */

	if (! write_long (fd, pm1data->stage, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->max_stage_0_prime, &sum)) goto writeerr;
	if (! write_longlong (fd, pm1data->B_done, &sum)) goto writeerr;
	if (! write_longlong (fd, pm1data->B, &sum)) goto writeerr;
	if (! write_longlong (fd, pm1data->C_done, &sum)) goto writeerr;
	if (! write_longlong (fd, pm1data->C_start, &sum)) goto writeerr;
	if (! write_longlong (fd, pm1data->C, &sum)) goto writeerr;
	if (! write_longlong (fd, processed, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->D, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->E, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->rels_done, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->bitarray_len, &sum)) goto writeerr;
	if (! write_array (fd, pm1data->bitarray, pm1data->bitarray_len, &sum)) goto writeerr;
	if (! write_longlong (fd, pm1data->bitarray_first_number, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->pairs_set, &sum)) goto writeerr;
	if (! write_long (fd, pm1data->pairs_done, &sum)) goto writeerr;

/* Write the data values.  There are occasions where gg may be in a */
/* partially FFTed state.  If so, do a harmless squaring to convert gg */
/* to an integer. */

	if (! write_gwnum (fd, &pm1data->gwdata, x, &sum)) goto writeerr;
	if (gg != NULL) {
		if (gwnum_is_partially_ffted (&pm1data->gwdata, gg)) {
			gwstartnextfft (&pm1data->gwdata, FALSE);
			gwsquare (&pm1data->gwdata, gg);
		}
		if (! write_gwnum (fd, &pm1data->gwdata, gg, &sum)) goto writeerr;
	}

/* Write the checksum, we're done */

	if (! write_checksum (fd, sum)) goto writeerr;

	closeWriteSaveFile (write_save_file_state, fd);
	return;

/* An error occured.  Close and delete the current file. */

writeerr:
	deleteWriteSaveFile (write_save_file_state, fd);
}

/* Read a save file */

int pm1_restore (			/* For version 25 and later save files */
	int	thread_num,
	pm1handle *pm1data,
	char	*filename,
	struct work_unit *w,
	uint64_t *processed,
	gwnum	*x,
	gwnum	*gg)
{
	int	fd;
	unsigned long version;
	unsigned long sum = 0, filesum;

/* Open the intermediate file */

	fd = _open (filename, _O_BINARY | _O_RDONLY);
	if (fd < 0) goto err;

/* Read the file header */

	if (! read_magicnum (fd, PM1_MAGICNUM)) goto readerr;
	if (! read_header (fd, &version, w, &filesum)) goto readerr;
	if (version < 1 || version > PM1_VERSION) goto readerr;

/* Read the file data */

	if (! read_long (fd, &pm1data->stage, &sum)) goto readerr;
	if (version == 2) pm1data->max_stage_0_prime = 13333333;	// The hardwired value prior to version 29.8 build 8
	else if (! read_long (fd, &pm1data->max_stage_0_prime, &sum)) goto readerr;
	if (! read_longlong (fd, &pm1data->B_done, &sum)) goto readerr;
	if (! read_longlong (fd, &pm1data->B, &sum)) goto readerr;
	if (! read_longlong (fd, &pm1data->C_done, &sum)) goto readerr;
	if (! read_longlong (fd, &pm1data->C_start, &sum)) goto readerr;
	if (! read_longlong (fd, &pm1data->C, &sum)) goto readerr;
	if (! read_longlong (fd, processed, &sum)) goto readerr;
	if (! read_long (fd, &pm1data->D, &sum)) goto readerr;
	if (! read_long (fd, &pm1data->E, &sum)) goto readerr;
	if (! read_long (fd, &pm1data->rels_done, &sum)) goto readerr;
	if (! read_long (fd, &pm1data->bitarray_len, &sum)) goto readerr;
	if (pm1data->bitarray_len) {
		pm1data->bitarray = (char *) malloc (pm1data->bitarray_len);
		if (pm1data->bitarray == NULL) goto readerr;
		if (! read_array (fd, pm1data->bitarray, pm1data->bitarray_len, &sum)) goto readerr;
	}
	if (! read_longlong (fd, &pm1data->bitarray_first_number, &sum)) goto readerr;
	if (! read_long (fd, &pm1data->pairs_set, &sum)) goto readerr;
	if (! read_long (fd, &pm1data->pairs_done, &sum)) goto readerr;

/* Note version 29.4 build 7 changed the calc_exp algorithm which invalidates earlier save files that are in stage 0. */

	if (version == 1 && pm1data->stage == PM1_STAGE0) {
		OutputBoth (thread_num, "P-1 save file incompatible with this program version.  Restarting stage 1 from the beginning.\n");
		goto readerr;
	}

/* Read the values */

	*x = gwalloc (&pm1data->gwdata);
	if (*x == NULL) goto readerr;
	if (! read_gwnum (fd, &pm1data->gwdata, *x, &sum)) goto readerr;

	*gg = NULL;
	if (pm1data->stage == PM1_STAGE2) {
		*gg = gwalloc (&pm1data->gwdata);
		if (*gg == NULL) goto readerr;
		if (! read_gwnum (fd, &pm1data->gwdata, *gg, &sum)) goto readerr;
	}

/* Read and compare the checksum */

	if (filesum != sum) goto readerr;
	_close (fd);

/* All done */

	return (TRUE);

/* An error occured.  Cleanup and return. */

readerr:
	_close (fd);
err:
	return (FALSE);
}


/* Compute how many values we can allocate.  This function can calculate */
/* the value using either the maximum available memory or the currently */
/* available memory. */

int choose_pminus1_numvals (
	pm1handle *pm1data,
	int	use_max_mem,		/* True if calculation should use */
					/* maximum available memory */
	unsigned long *numvals)		/* Returned number of values we */
					/* can allocate. */
{
	unsigned int memory;		/* Available memory in MB */
	int	stop_reason;

/* Override numvals when QAing */

	if (QA_TYPE) {
		*numvals = QA_TYPE;
		return (0);
	}

/* We assume that 13 temporaries will provide us with a reasonable */
/* execution speed.  We must have a minimum of 5 temporaries. */

	if (use_max_mem)
		memory =  max_mem (pm1data->thread_num);
	else {
		unsigned int min_memory, desired_memory;
		min_memory = cvt_gwnums_to_mem (&pm1data->gwdata, 5);
		desired_memory = cvt_gwnums_to_mem (&pm1data->gwdata, 13);
		stop_reason = avail_mem (pm1data->thread_num, min_memory, desired_memory, &memory);
		if (stop_reason) return (stop_reason);

/* Factor in the multiplier that we set to less than 1.0 when we get unexpected */
/* memory allocation errors.  Make sure we can still allocate 5 temporaries. */

		memory = (unsigned int) (pm1data->pct_mem_to_use * (double) memory);
		if (memory < min_memory)
			return (avail_mem_not_sufficient (pm1data->thread_num, min_memory, desired_memory));
	}
	if (memory < 8) memory = 8;

/* Output a message telling us how much memory is available */

	if (!use_max_mem && NUM_WORKER_THREADS > 1) {
		char	buf[100];
		sprintf (buf, "Available memory is %dMB.\n", memory);
		OutputStr (pm1data->thread_num, buf);
	}

/* Compute the number of gwnum temporaries we can allocate. */

	*numvals = cvt_mem_to_gwnums (&pm1data->gwdata, memory);
	if (*numvals < 1) *numvals = 1;
	return (0);
}

/* Calculate the number of values relatively prime to D.  D is a multiple */
/* of 30, 210, or 2310. */

unsigned long calc_numrels (
	unsigned long d)
{
	if (d >= 2310) return (d / 2310 * 480);
	if (d >= 210) return (d / 210 * 48);
	return (d / 30 * 8);
}

/* Compute the cost of a particular P-1 stage 2 plan. */

double cost_pminus1_plan (
	uint64_t B,		/* Stage 2 start point */
	uint64_t C,		/* Stage 2 end point */
	unsigned long d,
	unsigned long e,	/* Suyama power */
	unsigned long numvals,	/* Temp gwnums available for stage 2 */
	int	using_t3)	/* True if t3 avoids a gwfftadd3 */
{
	double numprimes, numpairings, sets;
	unsigned long numrels, passes;
	double	cost;

/* Estimate the number of primes */

	numprimes = (double) C / (log ((double) C) - 1.0) -
		    (double) B / (log ((double) B) - 1.0);

/* Calculate the number of values relatively prime to D */

	numrels = calc_numrels (d);

/* Estimate the number of prime pairs */

	if (e >= 2)
		numpairings =
			numprimes / 2.0 *
			numprimes / ((double) (C-B) * (double) numrels / d);

/* Compute how many passes this will take to process all the numbers */
/* relatively prime to D. */

	passes = (unsigned long)
		ceil ((double) numrels / (numvals - (e+1) - using_t3));

/* Compute the nQx setup costs.  To calculate nQx values, one fd_init is */
/* required on each pass with an average start point of D/2.  A single */
/* fd_init does E+1 powerings of roughly E*log2(startpoint) bits each */
/* at a cost of about 1.5 multiplies per bit.  In other words, */
/* (E+1) * E*log2(startpoint) * 1.5. */

	cost = passes * (e+1) * e*_log2(d/2) * 1.5 +

/* Then there are the D/2 calls to fd_next at E multiplies each. */

		d/2 * e +

/* Compute the eQx setup costs.  To calculate eQx values, one fd_init is */
/* required on each pass with a start point of B. */

		passes * (e+1) * e*_log2(B) * 1.5;

/* If E=1 add the cost of (C-B)/D fd_next calls.  If E>=2, add the cost */
/* of (C-B)/(D+D) calls to fd_next (E multiplies). */

	if (e == 1)
		sets = ceil ((double) (C-B) / d);
	else
		sets = ceil ((double) (C-B) / (d+d));
	cost += passes * sets * e;

/* Finally, each prime pairing costs one multiply.  If E = 1, then there */
/* is no prime pairing.  If not using_t3 then each multiply costs more as */
/* there is an extra gwfftadd3 call. */

	if (e >= 2 && using_t3)
		cost += numprimes - numpairings;
	else if (e >= 2)
		cost += (numprimes - numpairings) * 1.1;
	else if (using_t3)
		cost += numprimes;
	else
		cost += numprimes * 1.1;

/* Return the resulting cost */

	return (cost);
}

/* Choose the best values for D and E.  One that reduces the number of */
/* multiplications, yet doesn't use too much memory. */

int choose_pminus1_plan (
	pm1handle *pm1data,
	struct work_unit *w)
{
	uint64_t B, C;
	unsigned long numvals, d, e, i;
	double	cost, best_cost;
	int	stop_reason;

/* Clear D and E in case we don't find any acceptable plan */

	pm1data->D = 0;
	pm1data->E = 0;

/* Handle case where there is no stage 2 */

	B = pm1data->C_start;		/* Stage 2 starting point */
	C = pm1data->C;			/* Stage 2 ending point */
	if (C <= B) return (0);

/* Calculate the number of temporaries we can use for nQx and eQx and t3. */
/* Base our decision on the maximum amount of memory available.  Also, */
/* the main loop uses one temp for gg, so subtract one from numvals. */

	stop_reason = choose_pminus1_numvals (pm1data, 1, &numvals);
	if (stop_reason) return (stop_reason);
	numvals--;

/* Try various values of D until we find the best one */

	best_cost = 1e99;
	for (d = 5 * 2310; d >= 30; ) {

/* Try various values of E and using_t3 until we find the best one */

		for (i = 0; i < 10; i++) {
			int	using_t3;

/* Calculate e and the using_t3 flag.  Don't cost out e values where we */
/* don't have enough memory. */

			if (i <= 1) e = 1;
			else if (i <= 3) e = 2;
			else if (i <= 5) e = 4;
			else if (i <= 7) e = 6;
			else e = 12;
			using_t3 = i & 1;
			if (numvals <= e + 1 + using_t3) break;

/* Calculate the cost of this stage 2 plan */

			cost = cost_pminus1_plan (B, C, d, e, numvals, using_t3);

/* Reward higher E values because they should find more factors */
/* These rewards are close to a complete guess.  Little studying has */
/* been done on how often higher e values will find a factor. */

			if (e == 4) cost *= .95;
			if (e == 6) cost *= .90;
			if (e == 12) cost *= .85;

/* Remember best cost and best d and e */

			if (cost < best_cost) {
				best_cost = cost;
				pm1data->D = d;
				pm1data->E = e;
			}
		}

/* Try next smaller value of d */

		if (d > 2310) d = d - 2310;
		else if (d > 210) d = d - 210;
		else d = d - 30;
	}

/* Return no error code */

	return (0);
}

/* Choose the best implementation of the pminus1 plan given the current */
/* memory settings.  We may decide to wait for more memory to be available. */
/* We may choose to use the t3 temporary variable. */

int choose_pminus1_implementation (
	pm1handle *pm1data,
	struct work_unit *w,
	int	*using_t3_result)	/* Should we use a temporary for t3 */
{
	uint64_t B, C;
	unsigned long numvals;
	int	using_t3;
	double	cost, best_cost;
	int	stop_reason;

/* Copy some pm1data variables for easier access */

	B = pm1data->C_start;		/* Stage 2 starting point */
	C = pm1data->C;			/* Stage 2 ending point */

/* Compute the number of values relatively prime to D */

	pm1data->numrels = calc_numrels (pm1data->D);

/* Calculate the number of temporaries we can use for nQx and eQx and t3. */
/* Base our decision on the maximum amount of memory available.  Also, */
/* the main loop uses one temp for gg, so subtract one from numvals. */

	stop_reason = choose_pminus1_numvals (pm1data, 0, &numvals);
	if (stop_reason) return (stop_reason);
	ASSERTG (numvals >= 5);
	numvals--;

/* If not much memory is available right now, try shrinking E so that */
/* we can make some progress now.  The 10 in the formula below is arbitrary. */

	if (pm1data->E > 2 && numvals < pm1data->E + 10) pm1data->E = 2;

/* Try with and without using_t3 to find the best cost. */

	best_cost = 1e99;
	for (using_t3 = 0; using_t3 < 2; using_t3++) {

/* If there isn't enough memory to run this scenario, then do not */
/* bother costing it out. */

		if (numvals <= pm1data->E + 1 + using_t3) break;

/* Calculate the cost of this stage 2 plan */

		cost = cost_pminus1_plan (B, C, pm1data->D, pm1data->E,
					  numvals, using_t3);

/* Remember best cost, rels_this_pass, and using_t3 */

		if (cost < best_cost) {
			best_cost = cost;
			*using_t3_result = using_t3;
			pm1data->rels_this_pass =
				numvals - (pm1data->E+1) - using_t3;
		}
	}

/* Adjust rels_this_pass down if it is too high */

	if (pm1data->rels_this_pass > pm1data->numrels - pm1data->rels_done)
		pm1data->rels_this_pass = pm1data->numrels - pm1data->rels_done;
	return (0);
}

/* Formula to convert a prime into its corresponding bit in the bitarray */

#define bitcvt(prime,pm1data)  ((prime - (pm1data)->bitarray_first_number) >> 1)

/* Fill the bit array in such a way that it maximizes prime pairings. */
/* This is really optimized for P-1 on big Mersenne numbers.  I say this */
/* because for smaller numbers, you are apt to use large B2 values and */
/* you get a big bit array allocated.  And if B2 is really large then the */
/* bit array must be created in chunks. */

int fill_pminus1_bitarray (
	pm1handle *pm1data)
{
	uint64_t adjusted_C_start, prime, clear, jprime, pair, m, first_m;
	unsigned long max_bitarray_size, stage2incr, i;
	unsigned long *j, *jset;
	unsigned long relp[] = {7,11,13,17,19,23,29,31,37,41,43,47,0};
	int	stop_reason;

/* Process stage 2 in chunks if the bit array will be really large. */
/* By default, the bit array is limited to 250MB. Remember each byte */
/* corresponds to 8 odd numbers which is a range of 16. */

	max_bitarray_size = IniGetInt (INI_FILE, "MaximumBitArraySize", 250);
	if (max_bitarray_size > 2000) max_bitarray_size = 2000;
	if ((pm1data->C - pm1data->C_start) / 1000000 / 16 > max_bitarray_size)
		pm1data->C = pm1data->C_start + max_bitarray_size * 1000000 * 16;

/* Make sure C_start is odd. */

	pm1data->C_start |= 1;

/* The bit array starts at the first multiple of D below C_start. */

	pm1data->bitarray_first_number =
		(pm1data->C_start / pm1data->D) * pm1data->D + 1;

/* If the range from C_start to C allows us to move all the smaller primes */
/* then set adjusted_C_start to the first number that cannot be moved to */
/* a higher spot in the bit array. */

	if (pm1data->D >= 2310)
		adjusted_C_start = pm1data->C / 13, jset = relp + 2;
	else if (pm1data->D >= 210)
		adjusted_C_start = pm1data->C / 11, jset = relp + 1;
	else
		adjusted_C_start = pm1data->C / 7, jset = relp;
	adjusted_C_start |= 1;

/* Allocate the bitarray, pad it so that the stage 2 bit testing loop */
/* does not examine unallocated memory */

	pm1data->bitarray_len = (unsigned long)
		(pm1data->C - pm1data->bitarray_first_number +
		 pm1data->D * 2 + 15) >> 4;
	pm1data->bitarray = (char *) malloc (pm1data->bitarray_len);
	if (pm1data->bitarray == NULL) {
		stop_reason = OutOfMemory (pm1data->thread_num);
errexit:	pm1data->bitarray_len = 0;
		free (pm1data->bitarray);
		pm1data->bitarray = NULL;
		return (stop_reason);
	}
	memset (pm1data->bitarray, 0, pm1data->bitarray_len);

/* Set one bit for each prime between C_start and C */

	stop_reason = start_sieve (pm1data->thread_num, pm1data->C_start, &pm1data->sieve_info);
	if (stop_reason) goto errexit;
	for (prime = sieve (pm1data->sieve_info);
	     prime <= pm1data->C;
	     prime = sieve (pm1data->sieve_info)) {
		bitset (pm1data->bitarray, bitcvt (prime, pm1data));
		stop_reason = stopCheck (pm1data->thread_num);
		if (stop_reason) goto errexit;
	}

/* Now "move" some of the primes around so that we both maximize pairings. */
/* We do this by moving prime to 13*prime or 17*prime, etc. (as long as */
/* the multiple of prime is also in the bit array). */

	stage2incr = (pm1data->E == 1) ? pm1data->D : pm1data->D + pm1data->D;
	first_m = (adjusted_C_start / pm1data->D + 1) * pm1data->D;
	for (prime = pm1data->C_start; prime < adjusted_C_start; prime+=2) {
		if (!bittst (pm1data->bitarray, bitcvt (prime, pm1data))) continue;
		clear = prime;
		for (j = jset; *j; j++) {
			jprime = *j * prime;
			if (jprime > pm1data->C) break;
			bitclr (pm1data->bitarray, bitcvt (clear, pm1data));
			bitset (pm1data->bitarray, bitcvt (jprime, pm1data));
			clear = jprime;
			if (jprime < adjusted_C_start) continue;

/* Test if jprime pairs up */

			if (pm1data->E == 1) break;
			m = (jprime - (first_m - pm1data->D)) / stage2incr *
				stage2incr + first_m;
			if (jprime < m) pair = m + (m - jprime);
			else pair = m - (jprime - m);
			if (bittst (pm1data->bitarray, bitcvt (pair, pm1data)))
				break;
		}
		stop_reason = stopCheck (pm1data->thread_num);
		if (stop_reason) goto errexit;
	}

/* Count the number of pairs.  This is used to calculate the stage 2 */
/* percent complete. */

	m = (adjusted_C_start < pm1data->C_start) ? pm1data->C_start : adjusted_C_start;
	m = (m / pm1data->D + 1) * pm1data->D;
	for (pm1data->pairs_set = 0; pm1data->C > m-pm1data->D; m += stage2incr) {
	    for (i = 1; i < pm1data->D; i += 2) {
		if (bittst (pm1data->bitarray, bitcvt (m - i, pm1data)))
			pm1data->pairs_set++;
		else if (pm1data->E > 1 &&
		         bittst (pm1data->bitarray, bitcvt (m + i, pm1data))) {
			bitset (pm1data->bitarray, bitcvt (m - i, pm1data));
			pm1data->pairs_set++;
		}
	    }
	}
	pm1data->pairs_done = 0;

/* All done */

	return (0);
}

/* Recursively compute exp used in initial 3^exp calculation of a P-1 */
/* factoring run.  Don't forget to include 2*n in exp when factoring */
/* Mersenne numbers since factors must be 1 mod 2n */

void calc_exp (
	pm1handle *pm1data,
	double	k,		/* K in K*B^N+C */
	unsigned long b,	/* B in K*B^N+C */
	unsigned long n,	/* N in K*B^N+C */
	signed long c,		/* C in K*B^N+C */
	mpz_t	g,		/* Variable to accumulate multiplied small primes */
	uint64_t B1,		/* P-1 stage 1 bound */
	uint64_t *p,		/* Variable to fetch next small prime into */
	unsigned long lower,
	unsigned long upper)
{
	unsigned long len;

/* Compute the number of result bits we are to calculate */

	len = upper - lower;

/* Use recursion to compute the exponent.  This will perform better */
/* because mpz_mul will be handling arguments of equal size. */

	if (len >= 1024) {
		mpz_t	x;
		calc_exp (pm1data, k, b, n, c, g, B1, p, lower, lower + (len >> 1));
		mpz_init (x);
		calc_exp (pm1data, k, b, n, c, x, B1, p, lower + (len >> 1), upper);
		mpz_mul (g, x, g);
		mpz_clear (x);
		return;
	}

/* For Mersenne numbers, 2^n-1, make sure we include 2n in the calculated exponent (since factors */
/* are of the form 2kn+1).  For generalized Fermat numbers, b^n+1 (n is a power of 2), make sure n */
/* is included in the calculated exponent as factors are of the form kn+1 (actually forum posters */
/* have pointed out that Fermat numbers should include 4n and generalized Fermat should include 2n). */
/* Heck, maybe other forms may also need n included, so just always include 2n -- it is very cheap. */

	if (lower == 0) mpz_set_ui (g, 2*n);
	else mpz_set_ui (g, 1);

/* Find all the primes in the range and use as many powers as possible */

	for ( ; *p <= B1 && mpz_sizeinbase (g, 2) < len; *p = sieve (pm1data->sieve_info)) {
		uint64_t val, max;
		val = *p;
		max = B1 / *p;
		while (val <= max) val *= *p;
		if (sizeof (unsigned int) == 4 && val > 0xFFFFFFFF) {
			mpz_t	mpz_val;
			mpz_init_set_d (mpz_val, (double) val);		/* Works for B1 up to 2^53 */
			mpz_mul (g, g, mpz_val);
			mpz_clear (mpz_val);
		} else
			mpz_mul_ui (g, g, (unsigned int) val);
	}
}

/* Main P-1 entry point */

int pminus1 (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	pm1handle pm1data;
	uint64_t B;		/* Stage 1 bound */
	uint64_t C_start;	/* Stage 2 starting point (usually B) */
	uint64_t C;		/* Stage 2 ending point */
	uint64_t processed;	/* Data read from save file */
	giant	N;		/* Number being factored */
	giant	factor;		/* Factor found, if any */
	mpz_t	exp;
	int	exp_initialized;
	uint64_t stage_0_limit, prime, m;
	unsigned long memused, SQRT_B;
	unsigned long numrels, first_rel, last_rel;
	unsigned long i, j, stage2incr, len, bit_number;
	unsigned long error_recovery_mode = 0;
	gwnum	x, gg, tmp_gwnum, t3;
	readSaveFileState read_save_file_state;	/* Manage savefile names during reading */
	writeSaveFileState write_save_file_state; /* Manage savefile names during writing */
	char	filename[32], buf[255], JSONbuf[4000], testnum[100];
	int	have_save_file;
	int	res, stop_reason, stage, saving, near_fft_limit, echk;
	double	one_over_len, one_over_B, one_pair_pct;
	double	base_pct_complete, last_output, last_output_t, last_output_r;
	double	allowable_maxerr, output_frequency, output_title_frequency;
	int	first_iter_msg;
	int	using_t3;	/* Indicates we are using the gwnum t3 */
				/* to avoid a gwfftadd3 in stage 2 */
	int	msglen;
	char	*str, *msg;
	double	timers[2];
	double	pct_mem_to_use;

/* Unless we get memory errors, use as much memory as we can */

	pct_mem_to_use = 1.0;

/* Clear pointers to allocated memory (so common error exit code knows */
/* what to free) */

	memset (&pm1data, 0, sizeof (pm1handle));
	N = NULL;
	factor = NULL;
	str = NULL;
	msg = NULL;
	exp_initialized = FALSE;

/* Init local copies of B1 and B2 */

	B = (uint64_t) w->B1;
	C_start = (uint64_t) w->B2_start;
	C = (uint64_t) w->B2;

/* Choose a default value for the second bound if none was specified */

	if (C == 0) C = B * 100;

/* Make sure C_start and C values make sense */

	if (C_start < B) C_start = B;
	if (C < B) C = B;

/* Output startup message, but only if work type is P-1.  Pfactor work */
/* type has already output a startup message. */

	gw_as_string (testnum, w->k, w->b, w->n, w->c);
	sprintf (buf, "%s P-1", testnum);
	title (thread_num, buf);
	if (w->work_type == WORK_PMINUS1) {
		if (C <= B)
			sprintf (buf, "P-1 on %s with B1=%.0f\n", testnum, (double) B);
		else
			sprintf (buf, "P-1 on %s with B1=%.0f, B2=%.0f\n", testnum, (double) B, (double) C);
		OutputStr (thread_num, buf);
		if (w->sieve_depth > 0.0) {
			double prob = guess_pminus1_probability (w);
			sprintf (buf, "Chance of finding a factor is an estimated %.3g%%\n", prob * 100.0);
			OutputStr (thread_num, buf);
		}
	}

/* Clear all timers */

restart:
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Init filename.  This is a little kludgy as we want to generate a P-1 */
/* save file that does not conflict with an LL or PRP save file name. */
/* Both save files can exist at the same time when stage 2 is delayed */
/* waiting for more memory. */

	tempFileName (w, filename);
	filename[0] = 'm';

/* Override silly bounds */

	if (B < 30) {
		OutputStr (thread_num, "Using minimum bound #1 of 30\n");
		B = 30;
	}
	if (C < B) C = B;

/* Perform setup functions.  This includes decding how big an FFT to */
/* use, allocating memory, calling the FFT setup code, etc. */

/* Zero all data before beginning.  Set thread number. */

	memset (&pm1data, 0, sizeof (pm1handle));
	pm1data.thread_num = thread_num;
	pm1data.pct_mem_to_use = pct_mem_to_use;

/* Setup the assembly code */

	gwinit (&pm1data.gwdata);
	gwset_sum_inputs_checking (&pm1data.gwdata, SUM_INPUTS_ERRCHK);
	if (IniGetInt (LOCALINI_FILE, "UseLargePages", 0)) gwset_use_large_pages (&pm1data.gwdata);
	if (IniGetInt (INI_FILE, "HyperthreadPrefetch", 0)) gwset_hyperthread_prefetch (&pm1data.gwdata);
	if (HYPERTHREAD_LL) {
		sp_info->normal_work_hyperthreads = IniGetInt (LOCALINI_FILE, "HyperthreadLLcount", CPU_HYPERTHREADS);
		gwset_will_hyperthread (&pm1data.gwdata, sp_info->normal_work_hyperthreads);
	}
	gwset_bench_cores (&pm1data.gwdata, NUM_CPUS);
	gwset_bench_workers (&pm1data.gwdata, NUM_WORKER_THREADS);
	if (ERRCHK) gwset_will_error_check (&pm1data.gwdata);
	else gwset_will_error_check_near_limit (&pm1data.gwdata);
	gwset_num_threads (&pm1data.gwdata, CORES_PER_TEST[thread_num] * sp_info->normal_work_hyperthreads);
	gwset_thread_callback (&pm1data.gwdata, SetAuxThreadPriority);
	gwset_thread_callback_data (&pm1data.gwdata, sp_info);
	gwset_safety_margin (&pm1data.gwdata, IniGetFloat (INI_FILE, "ExtraSafetyMargin", 0.0));
	gwset_minimum_fftlen (&pm1data.gwdata, w->minimum_fftlen);
	res = gwsetup (&pm1data.gwdata, w->k, w->b, w->n, w->c);
	if (res) {
		sprintf (buf, "Cannot initialize FFT code, errcode=%d\n", res);
		OutputBoth (thread_num, buf);
		return (STOP_FATAL_ERROR);
	}

/* A kludge so that the error checking code is not as strict. */

	pm1data.gwdata.MAXDIFF *= IniGetInt (INI_FILE, "MaxDiffMultiplier", 1);

/* More miscellaneous initializations */

	last_output = last_output_t = last_output_r = 0;
	gw_clear_fft_count (&pm1data.gwdata);
	first_iter_msg = TRUE;
	calc_output_frequencies (&pm1data.gwdata, &output_frequency, &output_title_frequency);

/* Output message about the FFT length chosen */

	{
		char	fft_desc[200];
		gwfft_description (&pm1data.gwdata, fft_desc);
		sprintf (buf, "Using %s\n", fft_desc);
		OutputStr (thread_num, buf);
	}

/* If we are near the maximum exponent this fft length can test, then we */
/* will roundoff check all multiplies */

	near_fft_limit = exponent_near_fft_limit (&pm1data.gwdata);
	gwsetnormroutine (&pm1data.gwdata, 0, ERRCHK || near_fft_limit, 0);

/* Figure out the maximum round-off error we will allow.  By default this is 27/64 when near the FFT limit and 26/64 otherwise. */
/* We've found that this default catches errors without raising too many spurious error messages.  We let the user override */
/* this default for user "Never Odd Or Even" who tests exponents well beyond an FFT's limit.  He does his error checking by */
/* running the first-test and double-check simultaneously. */

	allowable_maxerr = IniGetFloat (INI_FILE, "MaxRoundoffError", (float) (near_fft_limit ? 0.421875 : 0.40625));

/* Compute the number we are factoring */

	stop_reason = setN (&pm1data.gwdata, thread_num, w, &N);
	if (stop_reason) goto exit;

/* Check for a save file and read the save file.  If there is an error */
/* reading the file then restart the P-1 factoring job from scratch. */
/* Limit number of backup files we try */
/* to read in case there is an error deleting bad save files. */

	have_save_file = FALSE;
	readSaveFileStateInit (&read_save_file_state, thread_num, filename);
	writeSaveFileStateInit (&write_save_file_state, filename, 0);
	for ( ; ; ) {
		if (! saveFileExists (&read_save_file_state)) {
			/* If there were save files, they are all bad.  Report a message */
			/* and temporarily abandon the work unit.  We do this in hopes that */
			/* we can successfully read one of the bad save files at a later time. */
			/* This sounds crazy, but has happened when OSes get in a funky state. */
			if (read_save_file_state.a_non_bad_save_file_existed) {
				OutputBoth (thread_num, ALLSAVEBAD_MSG);
				return (0);
			}
			/* No save files existed, start from scratch. */
			break;
		}

		if (!pm1_restore (thread_num, &pm1data, read_save_file_state.current_filename, w, &processed, &x, &gg)) {
			/* Close and rename the bad save file */
			saveFileBad (&read_save_file_state);
			continue;
		}

		have_save_file = TRUE;
		break;
	}

/* Record the amount of memory being used by this thread.  Until we get to */
/* stage 2, P-1 uses as much memory as an LL test. */

	if (!have_save_file || pm1data.stage != PM1_STAGE2 || C <= pm1data.C_done)
		set_memory_usage (thread_num, 0, cvt_gwnums_to_mem (&pm1data.gwdata, 1));

/* Jump to the proper continuation point if processing a save file */

	if (have_save_file) {

/* Handle stage 0 save files.  If the B values do not match, then use */
/* the bound given in the save file -- this may well result in an */
/* increased execution time if the saved B is larger than the B passed */
/* in to this routine. */

		if (pm1data.stage == PM1_STAGE0) {
			bit_number = (unsigned long) processed;
			goto restart0;
		}

/* To avoid an infinite loop of repeatable roundoff errors, we square */
/* the value read in from the P-1 save file.  This won't affect our final */
/* results, but will change the FFT data. */

		if (error_recovery_mode) {
			gwstartnextfft (&pm1data.gwdata, FALSE);
			gwsetnormroutine (&pm1data.gwdata, 0, 0, 0);
			gwsquare_carefully (&pm1data.gwdata, x);
			pm1_save (&pm1data, &write_save_file_state, w, processed, x, gg);
			error_recovery_mode = 0;
		}

/* Handle stage 1 save files */

		if (pm1data.stage == PM1_STAGE1) {
			if (B <= processed) {
				pm1data.B_done = processed;
				pm1data.C_done = processed;
				pm1data.B = processed;
				goto restart2;
			}
			if (B < pm1data.B) pm1data.B = B;
			stop_reason = start_sieve (thread_num, processed + 1, &pm1data.sieve_info);
			if (stop_reason) goto exit;
			prime = sieve (pm1data.sieve_info);
			goto restart1;
		}

/* Handle stage 2 save files */

		if (pm1data.stage == PM1_STAGE2) {

/* Clear flag indicating we need to restart if the maximum amount of */
/* memory changes.  We cannot change P-1 bounds after we've picked our plan */

			clear_restart_if_max_memory_change (thread_num);

/* If B is larger than the one in the save file, then go back and */
/* do some more stage 1 processing.  Since this is very upsetting to */
/* an LL tester that has already begun stage 2 only do this for the */
/* non-LL tester. */

			if (B > pm1data.B_done && w->work_type == WORK_PMINUS1) {
				gwfree (&pm1data.gwdata, gg);
				goto more_B;
			}

/* If B is larger than the one in the save file, then use the one in the save */
/* file rather than discarding all the work done thusfar in stage 2. */

			if (B != pm1data.B_done) {
				B = pm1data.B_done;
				C_start = pm1data.B_done;
				sprintf (buf, "Ignoring suggested B1 value, using B1=%.0f from the save file\n", (double) B);
				OutputStr (thread_num, buf);
			}

/* If we've already done enough stage 2, go do the stage 2 GCD */

			if (C <= pm1data.C_done) {
				stage = 2;
				goto restart4;
			}

/* If we never really started stage 2, then do so now */

			if (pm1data.bitarray_len == 0) goto more_C;

/* If LL testing and bound #2 has changed then use the original bound #2. */
/* If explicit P-1 testing and bound #2 is larger in the save file then use the original bound #2. */
/* The user doing explicit P-1 testing that wants to discard the stage 2 work he has done thusfar */
/* and reduce the stage 2 bound must manually delete the save file. */

			if ((w->work_type != WORK_PMINUS1 && C != pm1data.C) ||
			    (w->work_type == WORK_PMINUS1 && C < pm1data.C)) {
				C = pm1data.C;
				sprintf (buf, "Ignoring suggested B2 value, using B2=%.0f from the save file\n", (double) C);
				OutputStr (thread_num, buf);
			}

/* Resume stage 2 */

			goto restart3b;
		}

/* Handle case where we have a completed save file (the PM1_DONE state) */

		if (B > pm1data.B_done) goto more_B;
		if (C > pm1data.C_done) goto restart3a;

/* Note: if C_start != B then the user is using the undocumented feature */
/* of doing stage 2 in pieces.  Assume he knows what he is doing */

		if (C_start != B) {
			pm1data.C_done = pm1data.B_done;
			goto restart3a;
		}

/* The save file indicates we've tested to these bounds already */

		sprintf (buf, "%s already tested to B1=%.0f and B2=%.0f.\n",
			 gwmodulo_as_string (&pm1data.gwdata),
			 (double) pm1data.B_done, (double) pm1data.C_done);
		OutputBoth (thread_num, buf);
		goto done;
	}

/* Start this P-1 run from scratch starting with x = 3 */

	bit_number = 0;
	x = gwalloc (&pm1data.gwdata);
	if (x == NULL) goto oom;
	dbltogw (&pm1data.gwdata, 3.0, x);
	pm1data.B_done = 0;
	pm1data.B = B;

/* Stage 0 pre-calculates an exponent that is the product of small primes.  Our default uses only small */
/* primes below 40,000,000 (roughly 60 million bits).  This is configurable starting in version 29.8 build 8. */

	pm1data.max_stage_0_prime = IniGetInt (INI_FILE, "MaxStage0Prime", 40000000);

/* First restart point.  Compute the big exponent (product of small primes).  Then compute 3^exponent. */
/* The exponent always contains 2*p.  We only use primes below max_stage_0_prime.  The rest of the */
/* exponentiation will be done one prime at a time in the second part of stage 1. */
/* This stage uses 2 transforms per exponent bit. */

restart0:
	strcpy (w->stage, "S1");
	w->pct_complete = 0.0;
	pm1data.stage = PM1_STAGE0;
	start_timer (timers, 0);
	start_timer (timers, 1);
	stop_reason = start_sieve (thread_num, 2, &pm1data.sieve_info);
	if (stop_reason) goto exit;
	prime = sieve (pm1data.sieve_info);
	stage_0_limit = (pm1data.B > pm1data.max_stage_0_prime) ? pm1data.max_stage_0_prime : pm1data.B;
	mpz_init (exp);  exp_initialized = TRUE;
	calc_exp (&pm1data, w->k, w->b, w->n, w->c, exp, pm1data.B, &prime, 0, (unsigned long) (stage_0_limit * 1.5));

/* Find number of bits, ignoring the most significant bit */

	len = (unsigned long) mpz_sizeinbase (exp, 2) - 1;
	one_over_len = 1.0 / (double) len;
	if (prime < B) one_over_len *= (double) prime / (double) B;

/* Now take the exponent and raise x to that power */

	gwsetmulbyconst (&pm1data.gwdata, 3);
	while (bit_number < len) {

/* To avoid an infinite loop of repeatable roundoff errors, carefully */
/* get us past the offending iteration. */

		if (error_recovery_mode && bit_number == error_recovery_mode) {
			gwstartnextfft (&pm1data.gwdata, FALSE);
			gwsetnormroutine (&pm1data.gwdata, 0, 0, mpz_tstbit (exp, len - bit_number - 1));
			gwsquare_carefully (&pm1data.gwdata, x);
			error_recovery_mode = 0;
			saving = TRUE;
		}

/* Set various flags.  They control whether error-checking or the next FFT can be started. */

		else {
			stop_reason = stopCheck (thread_num);
			saving = testSaveFilesFlag (thread_num);
			echk = stop_reason || saving || ERRCHK || near_fft_limit || ((bit_number & 127) == 64);

/* Either square x or square x and multiply it by three. */

#ifndef SERVER_TESTING
			gwstartnextfft (&pm1data.gwdata, !stop_reason && !saving && bit_number+1 != error_recovery_mode && bit_number+1 != len);
			gwsetnormroutine (&pm1data.gwdata, 0, echk, mpz_tstbit (exp, len - bit_number - 1));
			if (bit_number < 30) gwsquare_carefully (&pm1data.gwdata, x);
			else gwsquare (&pm1data.gwdata, x);
#endif
		}

/* Test for an error */

		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;
		bit_number++;

/* Calculate our stage 1 percentage complete */

		w->pct_complete = (double) bit_number * one_over_len;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s P-1 stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata));
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&pm1data.gwdata);
		}

/* Every N squarings, output a progress report */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			if (bit_number > 1)
				OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&pm1data.gwdata);
			first_iter_msg = FALSE;
		}

/* Every N squarings, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 && gw_get_fft_count (&pm1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && stop_reason)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.\n",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			writeResults (buf);
			last_output_r = gw_get_fft_count (&pm1data.gwdata);
		}

/* Check for escape and/or if its time to write a save file */

		if (stop_reason || saving) {
			pm1_save (&pm1data, &write_save_file_state, w, bit_number, x, NULL);
			if (stop_reason) goto exit;
		}
	}

/* If roundoff error recovery returned to restart0, but the roundoff error */
/* occurs after the above loop, then square the value and create a save file. */
/* This won't affect our final results, but will change the FFT data. */

	if (error_recovery_mode) {
		gwstartnextfft (&pm1data.gwdata, FALSE);
		gwsetnormroutine (&pm1data.gwdata, 0, 0, 0);
		gwsquare_carefully (&pm1data.gwdata, x);
		pm1_save (&pm1data, &write_save_file_state, w, bit_number, x, NULL);
		error_recovery_mode = 0;
	}

/* Do stage 0 cleanup */

	gwsetnormroutine (&pm1data.gwdata, 0, ERRCHK || near_fft_limit, 0);
	mpz_clear (exp), exp_initialized = FALSE;
	end_timer (timers, 0);
	end_timer (timers, 1);

/* This situation will probably never happen, but will handle it anyway */

	if (B > stage_0_limit && B < pm1data.B) pm1data.B = B;

/* Second restart point.  Do the larger primes of stage 1. */
/* This stage uses 2.5 transforms per exponent bit. */

restart1:
	one_over_B = 1.0 / (double) B;
	strcpy (w->stage, "S1");
	w->pct_complete = prime * one_over_B;
	start_timer (timers, 0);
	start_timer (timers, 1);
	pm1data.stage = PM1_STAGE1;
	tmp_gwnum = gwalloc (&pm1data.gwdata);
	if (tmp_gwnum == NULL) goto oom;
	SQRT_B = (unsigned long) sqrt ((double) pm1data.B);
	for ( ; prime <= pm1data.B; prime = sieve (pm1data.sieve_info)) {

/* Test for user interrupt, save files, and error checking */

		stop_reason = stopCheck (thread_num);
		saving = testSaveFilesFlag (thread_num);
		echk = stop_reason || saving || ERRCHK || near_fft_limit || ((prime & 127) == 127);
		gwsetnormroutine (&pm1data.gwdata, 0, echk, 0);

/* Apply as many powers of prime as long as prime^n <= B */

		if (prime > pm1data.B_done) {
			pm1_mul (&pm1data, x, tmp_gwnum, prime);
		}
		if (prime <= SQRT_B) {
			uint64_t mult, max;
			mult = prime;
			max = pm1data.B / prime;
			for ( ; ; ) {
				mult *= prime;
				if (mult > pm1data.B_done) {
					pm1_mul (&pm1data, x, tmp_gwnum, prime);
				}
				if (mult > max) break;
			}
		}

/* Test for an error */

		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;

/* Calculate our stage 1 percentage complete */

		w->pct_complete = (double) prime * one_over_B;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s P-1 stage 1",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata));
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&pm1data.gwdata);
		}

/* Every N primes, output a progress report */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&pm1data.gwdata);
			first_iter_msg = FALSE;
		}

/* Every N primes, output a progress report to the results file */

		if ((ITER_OUTPUT_RES != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && stop_reason)) {
			sprintf (buf, "%s stage 1 is %.*f%% complete.\n",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			writeResults (buf);
			last_output_r = gw_get_fft_count (&pm1data.gwdata);
		}

/* Check for escape and/or if its time to write a save file */

		if (stop_reason || saving) {
			pm1_save (&pm1data, &write_save_file_state, w, prime, x, NULL);
			if (stop_reason) goto exit;
		}
	}
	gwfree (&pm1data.gwdata, tmp_gwnum);
	pm1data.B_done = pm1data.B;
	pm1data.C_done = pm1data.B;
	end_timer (timers, 0);
	end_timer (timers, 1);

/* Check for the rare case where we need to do even more stage 1 */
/* This happens when a save file was created with a smaller bound #1 */
/* than the bound #1 passed into this routine */

	if (B > pm1data.B) {
more_B:		pm1data.B = B;
		stop_reason = start_sieve (thread_num, 2, &pm1data.sieve_info);
		if (stop_reason) goto exit;
		prime = sieve (pm1data.sieve_info);
		goto restart1;
	}

/* Stage 1 complete, print a message */

	sprintf (buf, "%s stage 1 complete. %.0f transforms. Time: ",
		 gwmodulo_as_string (&pm1data.gwdata),
		 gw_get_fft_count (&pm1data.gwdata));
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));
	last_output = last_output_t = last_output_r = 0;
	gw_clear_fft_count (&pm1data.gwdata);

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&pm1data.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&pm1data.gwdata);
	}

/* Check to see if we found a factor - do GCD (x-1, N) */

restart2:
	strcpy (w->stage, "S1");
	w->pct_complete = 1.0;
	if (C <= B ||
	    (!QA_IN_PROGRESS && IniGetInt (INI_FILE, "Stage1GCD", 1))) {
		if (w->work_type != WORK_PMINUS1)
			OutputStr (thread_num, "Starting stage 1 GCD - please be patient.\n");
		start_timer (timers, 0);
		gwaddsmall (&pm1data.gwdata, x, -1);
		stop_reason = gcd (&pm1data.gwdata, thread_num, x, N, &factor);
		gwaddsmall (&pm1data.gwdata, x, 1);
		if (stop_reason) {
			pm1_save (&pm1data, &write_save_file_state, w, B, x, NULL);
			goto exit;
		}
		end_timer (timers, 0);
		strcpy (buf, "Stage 1 GCD complete. Time: ");
		print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
		OutputStr (thread_num, buf);
		stage = 1;
		if (factor != NULL) goto bingo;
	}

/* Skip second stage if so requested */

	if (C <= B) goto msg_and_exit;

/*
   Stage 2:  Use ideas from Crandall, Zimmermann, and Montgomery on each
   prime below C.  This code is more efficient the more memory you can
   give it.
   x: value at the end of stage 1 
*/

/* Initialize variables for second stage.  We set gg to x-1 in case the */
/* user opted to skip the GCD after stage 1. */

restart3a:
	sprintf (buf, "%s P-1 stage 2 init", gwmodulo_as_string (&pm1data.gwdata));
	title (thread_num, buf);
	strcpy (w->stage, "S2");
	w->pct_complete = 0.0;
	gg = gwalloc (&pm1data.gwdata);
	if (gg == NULL) goto oom;
	gwcopy (&pm1data.gwdata, x, gg);
	gwaddsmall (&pm1data.gwdata, gg, -1);
	pm1data.stage = PM1_STAGE2;

/* Choose a good value for D and E - one that reduces the number of */
/* multiplications, yet doesn't use too much memory.  This plan is based */
/* on the maximum available memory.  Once the plan is selected it cannot */
/* be changed in multi-pass stage 2 runs. */

more_C:	pm1data.C_start = (C_start > pm1data.C_done) ? C_start : pm1data.C_done;
	pm1data.C = C;
	stop_reason = choose_pminus1_plan (&pm1data, w);
	if (stop_reason) {
		pm1_save (&pm1data, &write_save_file_state, w, 0, x, gg);
		goto exit;
	}
	if (pm1data.D == 0) { /* We'll never have enough memory for stage 2 */
		OutputStr (thread_num, "Insufficient memory to ever run stage 2.\n");
		C = pm1data.B_done;
		goto restart4;
	}
	stop_reason = fill_pminus1_bitarray (&pm1data);
	if (stop_reason) {
		pm1_save (&pm1data, &write_save_file_state, w, 0, x, gg);
		goto exit;
	}
	pm1data.rels_done = 0;

/* Clear flag indicating we need to restart if the maximum amount of */
/* memory changes.  We cannot change P-1 bounds after we've picked our plan */

	clear_restart_if_max_memory_change (thread_num);

/* Restart here when in the middle of stage 2 or */
/* move to the next pass of a multi-pass stage 2 run */

restart3b:
	sprintf (buf, "%s P-1 stage 2 init", gwmodulo_as_string (&pm1data.gwdata));
	title (thread_num, buf);
	stage = 2;
	strcpy (w->stage, "S2");
	base_pct_complete =
		(double) (pm1data.C_start - C_start) / (double) (C - C_start);
	one_pair_pct =
		1.0 / (double) pm1data.pairs_set *
		(double) (pm1data.C - pm1data.C_start) / (double) (C - C_start);
	w->pct_complete = base_pct_complete +
		(double) pm1data.pairs_done * one_pair_pct;

/* Choose the best plan implementation given the currently available memory. */
/* This implementation could be anything from "wait until we have more */
/* memory" to deciding whether using_t3 should be set. */

replan:	stop_reason = choose_pminus1_implementation (&pm1data, w, &using_t3);
	if (stop_reason) {
		pm1_save (&pm1data, &write_save_file_state, w, 0, x, gg);
		goto exit;
	}

/* Record the amount of memory we intend to use.  We use rels_this_pass */
/* gwnums in the NQx array, E+1 gwnums to calculate eQx values, one gwnum */
/* for gg, and an optional gwnum for t3. */

	memused = cvt_gwnums_to_mem (&pm1data.gwdata, pm1data.rels_this_pass + pm1data.E + 2 + using_t3);
	if (set_memory_usage (thread_num, MEM_VARIABLE_USAGE, memused)) goto replan;
	sprintf (buf,
		 "Using %luMB of memory.  Processing %lu relative primes (%lu of %lu already processed).\n",
		 memused, pm1data.rels_this_pass, pm1data.rels_done, pm1data.numrels);
	OutputStr (thread_num, buf);

/* Here is where we restart the next pass of a multi-pass stage 2 */

	start_timer (timers, 0);
	start_timer (timers, 1);

/* On first pass, allocate P-1 stage 2 memory */

	if (pm1data.nQx == NULL) {
		pm1data.nQx = (gwnum *)
			malloc ((pm1data.D >> 1) * sizeof (gwnum));
		if (pm1data.nQx == NULL) goto lowmem;
		pm1data.eQx = (gwnum *)
			malloc ((pm1data.E + 1) * sizeof (gwnum));
		if (pm1data.eQx == NULL) goto lowmem;
	}

/* Clear the nQx array for this pass */

	memset (pm1data.nQx, 0, (pm1data.D >> 1) * sizeof (gwnum));

/* Compute x^(1^e), x^(3^e), ..., x^((D-1)^e) */

	for (i = 1, j = 0; ; i += 2) {
		if (! relatively_prime (i, pm1data.D)) continue;
		if (++j > pm1data.rels_done) break;
	}
	first_rel = i;
	gwfft (&pm1data.gwdata, x, x);		/* fd_init requires fft of x */
	stop_reason = fd_init (&pm1data, i, 2, x);
	if (stop_reason) goto exit;
	for (numrels = 0; ; ) {			/* Compute x^(i^e) */
		if (relatively_prime (i, pm1data.D)) {
			j = (i - 1) >> 1;
			pm1data.nQx[j] = gwalloc (&pm1data.gwdata);
			if (pm1data.nQx[j] == NULL) {
				gwstartnextfft (&pm1data.gwdata, FALSE);
				gwfftfftmul (&pm1data.gwdata, x, x, x);	/* Unfft x for save */
				goto lowmem;
			}
			gwcopy (&pm1data.gwdata, pm1data.eQx[0], pm1data.nQx[j]);
			numrels++;
			last_rel = i;
		}
		i = i + 2;
		if (i >= pm1data.D) break;
		if (numrels == pm1data.rels_this_pass) break;
		fd_next (&pm1data);
		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;
		stop_reason = stopCheck (thread_num);
		if (stop_reason) {
			fd_term (&pm1data);
			gwstartnextfft (&pm1data.gwdata, FALSE);
			gwfftfftmul (&pm1data.gwdata, x, x, x);	/* Unfft x - generates x^2 */
			pm1_save (&pm1data, &write_save_file_state, w, 0, x, gg);
			goto exit;
		}
	}
	fd_term (&pm1data);

/* Compute m = CEIL(start/D)*D, the first group we work on in stage 2 */
/* For the count of paired primes to be accurate, this code must exactly mirror */
/* the calculation of adjusted_C_start and first_m in fill_pminus1_bitarray. */	

	if (pm1data.D >= 2310) m = pm1data.C / 13;
	else if (pm1data.D >= 210) m = pm1data.C / 11;
	else m = pm1data.C / 7;
	m |= 1;
	if (m < pm1data.C_start) m = pm1data.C_start;
	m = (m / pm1data.D + 1) * pm1data.D;
	stage2incr = (pm1data.E == 1) ? pm1data.D : pm1data.D + pm1data.D;

/* Scan the bit array until we find the first group with a bit set. */
/* When continuing from a save file there could be many groups that */
/* have already been completed. */

	for ( ; pm1data.C > m-pm1data.D; m += stage2incr) {
	    for (i = first_rel; i <= last_rel; i += 2) {
		if (pm1data.nQx[i>>1] == NULL) continue;
		if (bittst (pm1data.bitarray, bitcvt (m - i, &pm1data)))
			goto found_a_bit;
	    }
	}
found_a_bit:;

/* Initialize for computing successive x^(m^e) */

	fd_init (&pm1data, m, stage2incr, x);

/* Unfft x for use in save files.  Actually this generates x^2 which */
/* is just fine - no stage 2 factors will be missed (in fact it could */
/* find more factors) */

	gwstartnextfft (&pm1data.gwdata, FALSE);
	gwfftfftmul (&pm1data.gwdata, x, x, x);

/* Now touch all the nQx and eQx values so that when gg is used, x is */
/* swapped out rather than a value we will need in the near future. */
/* In other words, make the gwnum x the least-recently-used. */

	for (i = 0; i <= pm1data.E; i++)
		gwtouch (&pm1data.gwdata, pm1data.eQx[i]);
	for (i = first_rel; i < last_rel; i += 2) {
		j = i >> 1;
		if (pm1data.nQx[j] != NULL)
			gwtouch (&pm1data.gwdata, pm1data.nQx[j]);
	}

/* Stage 2 init complete, change the title */

	{
		sprintf (buf, "%.*f%% of %s P-1 stage 2 (using %dMB)",
			 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata), (int) memused);
		title (thread_num, buf);
	}

/* When E >= 2, we can do prime pairing and each loop iteration */
/* handles the range m-D to m+D.  When E = 1, each iteration handles */
/* the range m-D to m. */

	if (using_t3) {
		t3 = gwalloc (&pm1data.gwdata);
		if (t3 == NULL) goto lowmem;
	}
	for ( ; pm1data.C > m-pm1data.D; m += stage2incr) {
	    int	inner_loop_done = FALSE;
	    int	last_pass = (m - pm1data.D + stage2incr >= pm1data.C);
	    saving = testSaveFilesFlag (thread_num);

/* Test all the relprimes between m-D and m */

	    for (i = first_rel; ; i += 2) {

/* Move onto the next m value when we are done with all the relprimes */

		if (i > last_rel) {	/* Compute next x^(m^e) */
			if (!last_pass) fd_next (&pm1data);
			inner_loop_done = TRUE;
			stop_reason = stopCheck (thread_num);
			goto errchk;
		}

/* Skip this relprime if we aren't processing it this pass */ 

		j = i >> 1;
		if (pm1data.nQx[j] == NULL) continue;

/* Skip this relprime if neither m - i nor its pair m + i are set */
/* in the bitarray. */

		if (! bittst (pm1data.bitarray, bitcvt (m - i, &pm1data)))
			continue;

/* Mul this eQx - nQx value into gg */

		stop_reason = stopCheck (thread_num);
		gwstartnextfft (&pm1data.gwdata, !stop_reason && !saving);
#ifndef SERVER_TESTING
		if (using_t3) {
			gwfftsub3 (&pm1data.gwdata, pm1data.eQx[0], pm1data.nQx[j], t3);
			gwfftmul (&pm1data.gwdata, t3, gg);
		} else {
			gwfftsub3 (&pm1data.gwdata, pm1data.eQx[0], pm1data.nQx[j], pm1data.eQx[0]);
			gwfftmul (&pm1data.gwdata, pm1data.eQx[0], gg);
			gwfftadd3 (&pm1data.gwdata, pm1data.eQx[0], pm1data.nQx[j], pm1data.eQx[0]);
		}
#endif

/* Clear this bit or bits in case a save file is written. */
/* Calculate stage 2 percentage. */

		bitclr (pm1data.bitarray, bitcvt (m - i, &pm1data));
		if (pm1data.E >= 2)
			bitclr (pm1data.bitarray, bitcvt (m + i, &pm1data));
		pm1data.pairs_done++;
		w->pct_complete = base_pct_complete +
			(double) pm1data.pairs_done * one_pair_pct;

/* Test for errors */

errchk:		if (gw_test_for_error (&pm1data.gwdata) || gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) goto err;

/* Output the title every so often */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output_t + 2 * ITER_OUTPUT * output_title_frequency)) {
			sprintf (buf, "%.*f%% of %s P-1 stage 2 (using %dMB)",
				 (int) PRECISION, trunc_percent (w->pct_complete), gwmodulo_as_string (&pm1data.gwdata), (int) memused);
			title (thread_num, buf);
			last_output_t = gw_get_fft_count (&pm1data.gwdata);
		}

/* Write out a message every now and then */

		if (first_iter_msg ||
		    (ITER_OUTPUT != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output + 2 * ITER_OUTPUT * output_frequency)) {
			sprintf (buf, "%s stage 2 is %.*f%% complete.",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			end_timer (timers, 0);
			if (first_iter_msg) {
				strcat (buf, "\n");
				clear_timer (timers, 0);
			} else {
				strcat (buf, " Time: ");
				print_timer (timers, 0, buf, TIMER_NL | TIMER_OPT_CLR);
			}
			OutputStr (thread_num, buf);
			start_timer (timers, 0);
			last_output = gw_get_fft_count (&pm1data.gwdata);
			first_iter_msg = FALSE;
		}

/* Write out a message to the results file every now and then */

		if ((ITER_OUTPUT_RES != 999999999 &&
		     gw_get_fft_count (&pm1data.gwdata) >= last_output_r + 2 * ITER_OUTPUT_RES) ||
		    (NO_GUI && stop_reason)) {
			sprintf (buf, "%s stage 2 is %.*f%% complete.\n",
				 gwmodulo_as_string (&pm1data.gwdata), (int) PRECISION, trunc_percent (w->pct_complete));
			writeResults (buf);
			last_output_r = gw_get_fft_count (&pm1data.gwdata);
		}

/* Periodicly write a save file.  If we escaped, free eQx memory so */
/* that pm1_save can reuse it to convert x and gg to binary.  If we */
/* have been using t3 as a temporary, free that for the same reason. */
/* "Touch" gg so that in low memory situations, the reading in of x */
/* swaps out one of the eQx or nQx values rather than gg. */

		if (stop_reason || saving) {
			if (stop_reason) fd_term (&pm1data);
			if (using_t3) gwfree (&pm1data.gwdata, t3);
			gwtouch (&pm1data.gwdata, gg);
			pm1_save (&pm1data, &write_save_file_state, w, 0, x, gg);
			if (stop_reason) goto exit;
			saving = FALSE;
			if (using_t3) {
				t3 = gwalloc (&pm1data.gwdata);
				if (t3 == NULL) goto oom;
			}
		}

/* Leave inner loop to work on the next m value */

		if (inner_loop_done) break;
	    }
	}
	if (using_t3) gwfree (&pm1data.gwdata, t3);
	fd_term (&pm1data);

/* Free up the nQx values for the next pass */

	for (i = first_rel; i <= last_rel; i += 2) {
		j = i >> 1;
		if (pm1data.nQx[j] != NULL)
			gwfree (&pm1data.gwdata, pm1data.nQx[j]);
	}

/* Check to see if another pass is required */

	end_timer (timers, 0);
	end_timer (timers, 1);
	pm1data.rels_done += pm1data.rels_this_pass;
	if (pm1data.rels_done < pm1data.numrels) goto restart3b;
	free (pm1data.bitarray);
	pm1data.bitarray = NULL;
	pm1data.bitarray_len = 0;
	free (pm1data.nQx);
	pm1data.nQx = NULL;
	free (pm1data.eQx);
	pm1data.eQx = NULL;

/* Check for the rare cases where we need to do even more stage 2. */
/* This happens when a save file was created with a smaller bound #2 */
/* than the bound #2 passed into this routine.  This also happens when */
/* B2 is so large that we must create the bitarray in chunks. */

	pm1data.C_done = pm1data.C;
	if (C > pm1data.C_done) goto more_C;

/* Stage 2 is complete */

	sprintf (buf, "%s stage 2 complete. %.0f transforms. Time: ",
		 gwmodulo_as_string (&pm1data.gwdata),
		 gw_get_fft_count (&pm1data.gwdata));
	print_timer (timers, 1, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	clear_timers (timers, sizeof (timers) / sizeof (timers[0]));

/* Print out round off error */

	if (ERRCHK) {
		sprintf (buf, "Round off: %.10g\n", gw_get_maxerr (&pm1data.gwdata));
		OutputStr (thread_num, buf);
		gw_clear_maxerr (&pm1data.gwdata);
	}

/* Since we set gwstartnextfft above, we must do another harmless squaring */
/* here to make sure gg has not been partially FFTed.  We cannot convert gg */
/* to an integer for GCD if it has been partially FFTed. */

	if (gwnum_is_partially_ffted (&pm1data.gwdata, gg)) {
		gwstartnextfft (&pm1data.gwdata, FALSE);
		gwsquare (&pm1data.gwdata, gg);
	}

/* See if we got lucky! */

restart4:
	strcpy (w->stage, "S2");
	w->pct_complete = 1.0;
	if (w->work_type != WORK_PMINUS1)
		OutputStr (thread_num, "Starting stage 2 GCD - please be patient.\n");
	start_timer (timers, 0);
	stop_reason = gcd (&pm1data.gwdata, thread_num, gg, N, &factor);
	if (stop_reason) {
		pm1_save (&pm1data, &write_save_file_state, w, C, x, gg);
		goto exit;
	}
	pm1data.stage = PM1_DONE;
	end_timer (timers, 0);
	strcpy (buf, "Stage 2 GCD complete. Time: ");
	print_timer (timers, 0, buf, TIMER_NL | TIMER_CLR);
	OutputStr (thread_num, buf);
	if (factor != NULL) goto bingo;

/* Output line to results file indicating P-1 run */

msg_and_exit:
	sprintf (buf, "%s completed P-1, B1=%.0f", gwmodulo_as_string (&pm1data.gwdata), (double) B);
	if (C > B) {
		if (pm1data.E <= 2)
			sprintf (buf+strlen(buf), ", B2=%.0f", (double) C);
		else
			sprintf (buf+strlen(buf), ", B2=%.0f, E=%lu", (double) C, pm1data.E);
	}
	sprintf (buf+strlen(buf), ", Wh%d: %08lX\n", PORT, SEC5 (w->n, B, C));
	OutputStr (thread_num, buf);
	formatMsgForResultsFile (buf, w);
	writeResults (buf);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"NF", "exponent":45581713, "worktype":"P-1", "b1":50000, "b2":5000000, "brent-suyama":6, */
/* "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"NF\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"P-1\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%.0f", (double) B);
	if (C > B) {
		sprintf (JSONbuf+strlen(JSONbuf), ", \"b2\":%.0f", (double) C);
		if (pm1data.E > 2) sprintf (JSONbuf+strlen(JSONbuf), ", \"brent-suyama\":%lu", pm1data.E);
	}
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", pm1data.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, B, C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send P-1 completed message to the server.  Although don't do it for puny */
/* B1 values as this is just the user tinkering with P-1 factoring. */

	if (B >= 10000 || IniGetInt (INI_FILE, "SendAllFactorData", 0)) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, buf);
		pkt.result_type = PRIMENET_AR_P1_NOFACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		pkt.B1 = (double) B;
		pkt.B2 = (double) C;
		pkt.fftlen = gwfftlen (&pm1data.gwdata);
		pkt.done = (w->work_type == WORK_PMINUS1 || w->work_type == WORK_PFACTOR);
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* Create save file so that we can expand bound 1 or bound 2 at a later date */
/* If this is pre-factoring for an LL or PRP test, then delete the large */
/* save file. */

	if (w->work_type == WORK_PMINUS1 && IniGetInt (INI_FILE, "KeepPminus1SaveFiles", 1))
		pm1_save (&pm1data, &write_save_file_state, w, 0, x, NULL);
	else
		unlinkSaveFiles (&write_save_file_state);

/* Return stop code indicating success or work unit complete */ 

done:	if (w->work_type == WORK_PMINUS1 || w->work_type == WORK_PFACTOR)
		stop_reason = STOP_WORK_UNIT_COMPLETE;
	else {
		w->pminus1ed = 1;		// Flag to indicate LL test has completed P-1
		w->tests_saved = 0.0;		// Variable to indicate PRP test has completed P-1
		stop_reason = updateWorkToDoLine (thread_num, w);
		if (stop_reason) goto exit;
	}

/* Free memory and return */

exit:	pm1_cleanup (&pm1data);
	free (N);
	free (factor);
	free (str);
	free (msg);
	if (exp_initialized) mpz_clear (exp);
	return (stop_reason);

/* Low on memory, reduce memory settings and try again */

lowmem:	fd_term (&pm1data);
	pm1_save (&pm1data, &write_save_file_state, w, 0, x, gg);
	pm1_cleanup (&pm1data);
	free (N);
	N = NULL;
	OutputBoth (thread_num, "Memory allocation error.  Trying again using less memory.\n");
	pct_mem_to_use *= 0.8;
	goto restart;

/* We've run out of memory.  Print error message and exit. */

oom:	stop_reason = OutOfMemory (thread_num);
	goto exit;

/* Print a message if we found a factor! */

bingo:	if (stage == 1)
		sprintf (buf, "P-1 found a factor in stage #1, B1=%.0f.\n", (double) B);
	else if (pm1data.E <= 2)
		sprintf (buf, "P-1 found a factor in stage #2, B1=%.0f, B2=%.0f.\n", (double) B, (double) C);
	else
		sprintf (buf, "P-1 found a factor in stage #2, B1=%.0f, B2=%.0f, E=%lu.\n", (double) B, (double) C, pm1data.E);
	OutputBoth (thread_num, buf);

/* Allocate memory for the string representation of the factor and for */
/* a message.  Convert the factor to a string.  Allocate lots of extra space */
/* as formatMsgForResultsFile can append a lot of text. */	

	msglen = factor->sign * 10 + 400;
	str = (char *) malloc (msglen);
	if (str == NULL) {
		stop_reason = OutOfMemory (thread_num);
		goto exit;
	}
	msg = (char *) malloc (msglen);
	if (msg == NULL) {
		stop_reason = OutOfMemory (thread_num);
		goto exit;
	}
	gtoc (factor, str, msglen);

/* Validate the factor we just found */

	if (!testFactor (&pm1data.gwdata, w, factor)) {
		sprintf (msg, "ERROR: Bad factor for %s found: %s\n", gwmodulo_as_string (&pm1data.gwdata), str);
		OutputBoth (thread_num, msg);
		unlinkSaveFiles (&write_save_file_state);
		OutputStr (thread_num, "Restarting P-1 from scratch.\n");
		stop_reason = 0;
		goto error_restart;
	}

/* Output the validated factor */

	if (stage == 1)
		sprintf (msg, "%s has a factor: %s (P-1, B1=%.0f)\n",
			 gwmodulo_as_string (&pm1data.gwdata), str, (double) B);
	else if (pm1data.E <= 2)
		sprintf (msg, "%s has a factor: %s (P-1, B1=%.0f, B2=%.0f)\n",
			 gwmodulo_as_string (&pm1data.gwdata), str, (double) B, (double) C);
	else
		sprintf (msg, "%s has a factor: %s (P-1, B1=%.0f, B2=%.0f, E=%d)\n",
			 gwmodulo_as_string (&pm1data.gwdata), str, (double) B, (double) C, (int) pm1data.E);
	OutputStr (thread_num, msg);
	formatMsgForResultsFile (msg, w);
	writeResults (msg);

/* Format a JSON version of the result.  An example follows: */
/* {"status":"F", "exponent":45581713, "worktype":"P-1", "factors":["430639100587696027847"], */
/* "b1":50000, "b2":5000000, "brent-suyama":6, */
/* "fft-length":5120, "security-code":"39AB1238", */
/* "program":{"name":"prime95", "version":"29.5", "build":"8"}, "timestamp":"2019-01-15 23:28:16", */
/* "user":"gw_2", "cpu":"work_computer", "aid":"FF00AA00FF00AA00FF00AA00FF00AA00"} */

	strcpy (JSONbuf, "{\"status\":\"F\"");
	JSONaddExponent (JSONbuf, w);
	strcat (JSONbuf, ", \"worktype\":\"P-1\"");
	sprintf (JSONbuf+strlen(JSONbuf), ", \"factors\":[\"%s\"]", str);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"b1\":%.0f", (double) B);
	if (stage > 1) {
		sprintf (JSONbuf+strlen(JSONbuf), ", \"b2\":%.0f", (double) C);
		if (pm1data.E > 2) sprintf (JSONbuf+strlen(JSONbuf), ", \"brent-suyama\":%lu", pm1data.E);
	}
	sprintf (JSONbuf+strlen(JSONbuf), ", \"fft-length\":%lu", pm1data.gwdata.FFTLEN);
	sprintf (JSONbuf+strlen(JSONbuf), ", \"security-code\":\"%08lX\"", SEC5 (w->n, B, C));
	JSONaddProgramTimestamp (JSONbuf);
	JSONaddUserComputerAID (JSONbuf, w);
	strcat (JSONbuf, "}\n");
	if (IniGetInt (INI_FILE, "OutputJSON", 1)) writeResultsJSON (JSONbuf);

/* Send assignment result to the server.  To avoid flooding the server */
/* with small factors from users needlessly redoing factoring work, make */
/* sure the factor is more than 50 bits or so. */

	if (strlen (str) >= 15 || IniGetInt (INI_FILE, "SendAllFactorData", 0)) {
		struct primenetAssignmentResult pkt;
		memset (&pkt, 0, sizeof (pkt));
		strcpy (pkt.computer_guid, COMPUTER_GUID);
		strcpy (pkt.assignment_uid, w->assignment_uid);
		strcpy (pkt.message, msg);
		pkt.result_type = PRIMENET_AR_P1_FACTOR;
		pkt.k = w->k;
		pkt.b = w->b;
		pkt.n = w->n;
		pkt.c = w->c;
		truncated_strcpy (pkt.factor, sizeof (pkt.factor), str);
		pkt.B1 = (double) B;
		pkt.B2 = (double) (stage == 1 ? 0 : C);
		pkt.fftlen = gwfftlen (&pm1data.gwdata);
		pkt.done = TRUE;
		strcpy (pkt.JSONmessage, JSONbuf);
		spoolMessage (PRIMENET_ASSIGNMENT_RESULT, &pkt);
	}

/* If LL testing, free all save files -- including possible LL save files */

	if (w->work_type != WORK_PMINUS1 || !IniGetInt (INI_FILE, "KeepPminus1SaveFiles", 1)) {
		unlinkSaveFiles (&write_save_file_state);
		write_save_file_state.base_filename[0] = 'p';
		unlinkSaveFiles (&write_save_file_state);
	}

/* Otherwise create save file so that we can expand bound 1 or bound 2 */
/* at a later date. */

	else
		pm1_save (&pm1data, &write_save_file_state, w, 0, x, NULL);

/* Since we found a factor, then we may have performed less work than */
/* expected.  Make sure we do not update the rolling average with */
/* this inaccurate data. */

	invalidateNextRollingAverageUpdate ();

/* Remove the exponent from the worktodo.ini file */

	stop_reason = STOP_WORK_UNIT_COMPLETE;
	goto exit;

/* Output an error message saying we are restarting. */
/* Sleep five minutes before restarting from last save file. */

err:	if (gw_get_maxerr (&pm1data.gwdata) > allowable_maxerr) {
		sprintf (buf, "Possible roundoff error (%.8g), backtracking to last save file.\n", gw_get_maxerr (&pm1data.gwdata));
		OutputStr (thread_num, buf);
	} else {
		OutputBoth (thread_num, "SUMOUT error occurred.\n");
		stop_reason = SleepFive (thread_num);
		if (stop_reason) goto exit;
	}
	error_recovery_mode = bit_number ? bit_number : 1;
error_restart:
	pm1_cleanup (&pm1data);
	free (N), N = NULL;
	free (factor), factor = NULL;
	free (str), str = NULL;
	free (msg), msg = NULL;
	if (exp_initialized) mpz_clear (exp), exp_initialized = FALSE;
	goto restart;
}

/* Read a file of P-1 tests to run as part of a QA process */
/* The format of this file is: */
/*	k, n, c, B1, B2_start, B2_end, factor */
/* Use Advanced/Time 9992 to run the QA suite */

int pminus1_QA (
	int	thread_num,
	struct PriorityInfo *sp_info)	/* SetPriority information */
{
	FILE	*fd;

/* Set the title */

	title (thread_num, "QA");

/* Open QA file */

	fd = fopen ("qa_pm1", "r");
	if (fd == NULL) {
		OutputStr (thread_num, "File named 'qa_pm1' could not be opened.\n");
		return (STOP_FILE_IO_ERROR);
	}

/* Loop until the entire file is processed */

	QA_TYPE = 0;
	for ( ; ; ) {
		struct work_unit w;
		double	k;
		unsigned long b, n, B1, B2_start, B2_end;
		signed long c;
		char	fac_str[80];
		int	stop_reason;

/* Read a line from the file */

		n = 0;
		(void) fscanf (fd, "%lf,%lu,%lu,%ld,%lu,%lu,%lu,%s\n",
			&k, &b, &n, &c, &B1, &B2_start, &B2_end, fac_str);
		if (n == 0) break;

/* If p is 1, set QA_TYPE */

		if (n == 1) {
			QA_TYPE = c;
			continue;
		}

/* Convert the factor we expect to find into a "giant" type */

		QA_FACTOR = allocgiant ((int) strlen (fac_str));
		ctog (fac_str, QA_FACTOR);

/*test various num_tmps
test 4 (or more?) stage 2 code paths
print out each test case (all relevant data)*/

/* Do the P-1 */

		if (B2_start < B1) B2_start = B1;
		w.work_type = WORK_PMINUS1;
		w.k = k;
		w.b = b;
		w.n = n;
		w.c = c;
		w.B1 = B1;
		w.B2_start = B2_start;
		w.B2 = B2_end;
		QA_IN_PROGRESS = TRUE;
		stop_reason = pminus1 (0, sp_info, &w);
		QA_IN_PROGRESS = FALSE;
		free (QA_FACTOR);
		if (stop_reason) {
			fclose (fd);
			return (stop_reason);
		}
	}

/* Cleanup */

	fclose (fd);
	return (0);
}

/**************************************************************/
/* Routines to compute optimal and test to optimal P-1 bounds */
/**************************************************************/

struct global_pm1_cost_data {
	unsigned long n;				// Exponent being tested
	int	isMersenne;
	double	how_far_factored;			// How far the number has been trial factored
	double	gcd_cost;				// Cost (in squarings) of running GCD
	double	ll_testing_cost;			// Cost (in squarings) of running LL/PRP tests should we fail to find a factor
	unsigned long vals;				// Number of temporaries we can allocate in pass 2
};

struct cost_pm1_data {
	unsigned long B1;
	unsigned long B2;
	double	prob;
	double	pass1_squarings;
	double	pass2_squarings;
	double	savings;
};

/* For a given B1,B2 calculate the costs and savings */

void cost_pm1 (
	struct global_pm1_cost_data *g,
	struct cost_pm1_data *c)
{
	double	logB1, logB2, numprimes;

/* Compute how many squarings will be required in pass 2.  In the */
/* low-memory cases, assume choose_pminus1_plan will pick D = 210, E = 1 */
/* If more memory is available assume choose_pminus1_plan will pick */
/* D = 2310, E = 2.  This will provide an accurate enough cost for our */
/* purposes even if different D and E values are picked.  See */
/* choose_pminus1_plan for a description of the costs of P-1 stage 2. */

	logB1 = log ((double) c->B1);
	logB2 = log ((double) c->B2);
	numprimes = (unsigned long) (c->B2 / (logB2 - 1.0) - c->B1 / (logB1 - 1.0));
	if (c->B2 <= c->B1) {
		c->pass2_squarings = 0.0;
	} else if (g->vals <= 8) {		/* D = 210, E = 1, passes = 48/temps */
		unsigned long num_passes = (unsigned long) ceil (48.0 / (g->vals - 3));
		c->pass2_squarings = ceil ((c->B2 - c->B1) / 210.0) * num_passes;
		c->pass2_squarings += numprimes * 1.1;
	} else {
		unsigned long num_passes = (unsigned long) ceil (480.0 / (g->vals - 5));
		double	numpairings = (unsigned long) (numprimes / 2.0 * numprimes / ((c->B2 - c->B1) * 480.0/2310.0));
		c->pass2_squarings = 2400.0 + num_passes * 90.0; /* setup costs */
		c->pass2_squarings += ceil ((c->B2 - c->B1) / 4620.0) * 2.0 * num_passes;
		c->pass2_squarings += numprimes - numpairings;
	}

/* Pass 2 FFT multiplications seem to be at least 20% slower than */
/* the squarings in pass 1.  This is probably due to several factors. */
/* These include: better L2 cache usage and no calls to the faster */
/* gwsquare routine.  Nov, 2009:  On my Macbook Pro, with exponents */
/* around 45M and using 800MB memory, pass2 squarings are 40% slower. */

	c->pass2_squarings *= 1.35;

/* Calculate probability of finding a factor (courtesy of Mihai Preda) */

	c->prob = pm1prob (g->n, g->isMersenne, (unsigned int) g->how_far_factored, c->B1, c->B2);

/* Calculate our savings using this B1/B2.  Savings is success_probability * cost_of_LL_tests - cost_of_Pminus1. */

	c->savings =  c->prob * g->ll_testing_cost - (c->pass1_squarings + c->pass2_squarings + g->gcd_cost);
}

/* For a given B1, find the best B2 */

void bestB2 (
	struct global_pm1_cost_data *g,
	struct cost_pm1_data *c)
{
	struct cost_pm1_data best[3];

/* Estimate how many squarings will be required in pass 1 */

	c->pass1_squarings = ceil (1.44 * c->B1);

/* Look for the best B2 somewhere between 1*B1 and 100*B1 */

	best[0] = *c;
	best[0].B2 = c->B1;
	cost_pm1 (g, &best[0]);
	best[1] = *c;
	best[1].B2 = 50*c->B1;
	cost_pm1 (g, &best[1]);
	best[2] = *c;
	best[2].B2 = 100*c->B1;
	cost_pm1 (g, &best[2]);

/* Handle case where midpoint has worse savings than the start point */
/* The search code requires best[1] to be better than best[0] and best[2] */

	while (best[0].savings > best[1].savings) {
		best[2] = best[1];
		best[1].B2 = (best[0].B1 + best[2].B1) / 2;
		cost_pm1 (g, &best[1]);
	}

/* Handle case where midpoint has worse savings than the end point */
/* The search code requires best[1] to be better than best[0] and best[2] */

	while (best[2].savings > best[1].savings) {
		best[0] = best[1];
		best[1] = best[2];
		best[2].B2 = best[1].B2 * 2;
		cost_pm1 (g, &best[2]);
	}

/* Find the best B2.  We use a binary-like search to speed things up (new in version 30.3b3). */

	while (best[0].B2 != best[2].B2) {
		struct cost_pm1_data midpoint;

		ASSERTG (best[1].savings >= best[0].savings);
		ASSERTG (best[1].savings >= best[2].savings);

/* Work on the bigger of the lower section and upper section */

		if (best[1].B2 - best[0].B2 > best[2].B2 - best[1].B2) {	// Work on lower section
			// If B2's are close together or the savings difference is real small, then we've searched this section enough
			if (best[1].B2 - best[0].B2 < 1000 || best[1].savings - best[0].savings < 100.0) {
				best[0] = best[1];
				continue;
			}
			midpoint = *c;
			midpoint.B2 = (best[0].B2 + best[1].B2) / 2;
			cost_pm1 (g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			// If B2's are close together or the savings difference is real small, then we've searched this section enough
			if (best[2].B2 - best[1].B2 < 1000 || best[1].savings - best[2].savings < 100.0) {
				best[2] = best[1];
				continue;
			}
			midpoint = *c;
			midpoint.B2 = (best[1].B2 + best[2].B2) / 2;
			cost_pm1 (g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Return the best B2 we found */

	c->B2 = best[1].B2;
	c->prob = best[1].prob;
	c->pass2_squarings = best[1].pass2_squarings;
	c->savings = best[1].savings;
}

/* Analyze how well P-1 factoring will perform */

void guess_pminus1_bounds (
	int	thread_num,
	double	k,			/* K in K*B^N+C. Must be a positive integer. */
	unsigned long b,		/* B in K*B^N+C. Must be two. */
	unsigned long n,		/* N in K*B^N+C. Exponent to test. */
	signed long c,			/* C in K*B^N+C. */
	double	how_far_factored,	/* Bit depth of trial factoring */
	double	tests_saved,		/* 1 if doublecheck, 2 if first test */
	unsigned long *bound1,
	unsigned long *bound2,
	unsigned long *squarings,
	double	*success_rate)
{
	struct global_pm1_cost_data g;
	struct cost_pm1_data best[3];

/* Copy exponent, how_far_factored to global costing data */

	g.n = n;
	g.isMersenne = (k == 1.0 && b == 2 && c == -1);
	g.how_far_factored = how_far_factored;

/* Guard against wild tests_saved values.  Huge values will cause this routine */
/* to run for a very long time.  This shouldn't happen as auxiliaryWorkUnitInit */
/* now has the exact same test. */

	if (tests_saved > 10) tests_saved = 10;

/* Balance P-1 against 1 or 2 LL/PRP tests (actually more since we get a */
/* corrupt result reported some of the time). */

	g.ll_testing_cost = (tests_saved + 2 * ERROR_RATE) * n;

/* Precompute the cost of a GCD.  We used Excel to come up with the */
/* formula GCD is equivalent to 861 * Ln (p) - 7775 transforms. */
/* Since one squaring equals two transforms we get the formula below. */
/* NOTE: In version 22, the GCD speed has approximately doubled.  I've */
/* adjusted the formula accordingly. */

	g.gcd_cost = (430.5 * log ((double) n) - 3887.5) / 2.0;
	if (g.gcd_cost < 50.0) g.gcd_cost = 50.0;

/* Compute how many temporaries we can use given our memory constraints. */
/* Allow 1MB for code and data structures. */

	g.vals = cvt_mem_to_estimated_gwnums (max_mem (thread_num), k, b, n, c);
	if (g.vals < 1) g.vals = 1;

/* Find the best B1 somewhere between n/3300 and 250*(n/3300). */

	best[0].B1 = n / 3300;
	bestB2 (&g, &best[0]);
	best[1].B1 = 125 * best[0].B1;
	bestB2 (&g, &best[1]);
	best[2].B1 = 250 * best[0].B1;
	bestB2 (&g, &best[2]);

/* Handle case where midpoint has worse savings than the start point */
/* The search code requires best[1] to be better than best[0] and best[2] */

	while (best[0].savings > best[1].savings) {
		best[2] = best[1];
		best[1].B1 = (best[0].B1 + best[2].B1) / 2;
		bestB2 (&g, &best[1]);
	}

/* Handle case where midpoint has worse savings than the end point */
/* The search code requires best[1] to be better than best[0] and best[2] */

	while (best[2].savings > best[1].savings) {
		best[0] = best[1];
		best[1] = best[2];
		best[2].B1 = best[1].B1 * 2;
		bestB2 (&g, &best[2]);
	}

/* Find the best B1.  We use a binary-like search to speed things up (new in version 30.3b3). */

	while (best[0].B1 != best[2].B1) {
		struct cost_pm1_data midpoint;

		ASSERTG (best[1].savings >= best[0].savings);
		ASSERTG (best[1].savings >= best[2].savings);

/* Work on the bigger of the lower section and upper section */

		if (best[1].B1 - best[0].B1 > best[2].B1 - best[1].B1) {	// Work on lower section
			// If B1's are close together or the savings difference is real small, then we've searched this section enough
			if (best[1].B1 - best[0].B1 < 1000 || best[1].savings - best[0].savings < 100.0) {
				best[0] = best[1];
				continue;
			}
			midpoint.B1 = (best[0].B1 + best[1].B1) / 2;
			bestB2 (&g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new end point
				best[2] = best[1];
				best[1] = midpoint;
			} else {						// Create new start point
				best[0] = midpoint;
			}
		} else {							// Work on upper section
			// If B1's are close together or the savings difference is real small, then we've searched this section enough
			if (best[2].B1 - best[1].B1 < 1000 || best[1].savings - best[2].savings < 100.0) {
				best[2] = best[1];
				continue;
			}
			midpoint.B1 = (best[1].B1 + best[2].B1) / 2;
			bestB2 (&g, &midpoint);
			if (midpoint.savings > best[1].savings) {		// Make middle the new start point
				best[0] = best[1];
				best[1] = midpoint;
			} else {						// Create new end point
				best[2] = midpoint;
			}
		}
	}

/* Round up B1 and B2 to nearest 1000 -- just to look pretty */

	best[1].B1 = round_up_to_multiple_of (best[1].B1, 1000);
	best[1].B2 = round_up_to_multiple_of (best[1].B2, 1000);
	cost_pm1 (&g, &best[1]);
	
/* Return the final best choice */

	if (best[1].savings > 0.0) {
		*bound1 = best[1].B1;
		*bound2 = best[1].B2;
		*squarings = (unsigned long) (best[1].pass1_squarings + best[1].pass2_squarings + g.gcd_cost);
		*success_rate = best[1].prob;
	} else {
		*bound1 = 0;
		*bound2 = 0;
		*squarings = 0;
		*success_rate = 0.0;
	}
}

/* Determine the probability of P-1 finding a factor */

double guess_pminus1_probability (
	struct work_unit *w)
{

/* Return the Mihai estimated P-1 success probability */

	return (pm1prob (w->n, w->k == 1.0 && w->b == 2 && w->c == -1, (unsigned int) w->sieve_depth, w->B1, w->B2));
}

/* Do the P-1 factoring step prior to a Lucas-Lehmer test */
/* Similar to the main P-1 entry point, except bounds are not known */

int pfactor (
	int	thread_num,
	struct PriorityInfo *sp_info,	/* SetPriority information */
	struct work_unit *w)
{
	unsigned long bound1, bound2, squarings;
	double	prob;
	char	buf[120], testnum[120];
	int	stop_reason;

/* Choose the best FFT size */

	stop_reason = pick_fft_size (thread_num, w);
	if (stop_reason) return (stop_reason);

/* Make sure the first-time user runs a successful self-test. */
/* The one-hour self-test may have been useful when it was first introduced */
/* but I think it now does little to catch buggy machines (they eventually */
/* work OK for an hour) and does create user confusion and annoyance. */

#ifdef ONE_HOUR_SELF_TEST
	if (w->k == 1.0 && w->b == 2 && w->c == -1) {
		stop_reason = selfTest (thread_num, sp_info, w);
		if (stop_reason) return (stop_reason);
	}
#endif

/* Set flag indicating we need to restart if the maximum amount of */
/* memory changes (as opposed to the currently available memory!) */
/* If maximum memory changes we want to recompute the P-1 bounds. */

	set_restart_if_max_memory_change (thread_num);

/* Output a message that P-1 factoring is about to begin */

	gw_as_string (testnum, w->k, w->b, w->n, w->c);
	sprintf (buf, "Optimal P-1 factoring of %s using up to %luMB of memory.\n",
		 testnum, max_mem (thread_num));
	OutputStr (thread_num, buf);
	sprintf (buf, "Assuming no factors below 2^%.2g and %.2g primality test%s saved if a factor is found.\n",
		 w->sieve_depth, w->tests_saved,
		 w->tests_saved == 1.0 ? "" : "s");
	OutputStr (thread_num, buf);

/* Deduce the proper P-1 bounds */

	guess_pminus1_bounds (thread_num, w->k, w->b, w->n, w->c, w->sieve_depth,
			      w->tests_saved, &bound1, &bound2,
			      &squarings, &prob);
	if (bound1 == 0) {
		sprintf (buf, "%s does not need P-1 factoring.\n", testnum);
		OutputBoth (thread_num, buf);
		if (w->work_type == WORK_PFACTOR) {
			//bug - do we need to tell the server to cancel the
			//assignment?  In theory, server shouldn't ever send
			//this assignment out.
			return (STOP_WORK_UNIT_COMPLETE);
		} else {
			w->pminus1ed = 1;		// Flag to indicate LL test has completed P-1
			w->tests_saved = 0.0;		// Variable to indicate PRP test has completed P-1
			stop_reason = updateWorkToDoLine (thread_num, w);
			if (stop_reason) return (stop_reason);
			return (0);
		}
	}

/* Output a message that P-1 factoring is about to begin */

	sprintf (buf, "Optimal bounds are B1=%ld, B2=%ld\n", bound1, bound2);
	OutputStr (thread_num, buf);
	sprintf (buf, "Chance of finding a factor is an estimated %.3g%%\n", prob * 100.0);
	OutputStr (thread_num, buf);

/* Call the P-1 factoring code */

	w->B1 = bound1;
	w->B2_start = bound1;
	w->B2 = bound2;
	return (pminus1 (thread_num, sp_info, w));
}
