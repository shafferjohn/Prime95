// test2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "gwnum.h"

void foo (int x, char *str)
{
	printf ("%s", str);
}

#ifdef SIMPLE
int _tmain(int argc, _TCHAR* argv[])
{
	gwhandle gwdata;
	gwnum	x;

	OutputBothRoutine = foo;
	gwinit (&gwdata);
	gwsetup (&gwdata, 1.0, 2, 640, -1);

	x = gwalloc (&gwdata);
	dbltogw (&gwdata, 2.0, x);
	gwsetnormroutine (&gwdata, 0, 1, 0); /* Enable error checking */
	gwsquare (&gwdata, x);

	return 0;
}
#endif

/* Compare gwnum to a giant */

void compare (gwhandle *gwdata, gwnum x, giant g)
{
	giant	tmp;

	tmp = popg (&gwdata->gdata, ((int) gwdata->bit_length >> 5) + 13);
	gwtogiant (gwdata, x, tmp);
	if (gcompg (g, tmp) != 0)
		printf ("Test failed.\n");
	pushg (&gwdata->gdata, 1);
}

/* Generate random data to start the test */

void gen_data (gwhandle *gwdata, gwnum x, giant g)
{
	unsigned long i, len;
	int	seed;

/* Set and output seed so that we can re-generate the random data */

	seed = (int) time (NULL);
	srand (seed);
	printf ("Random seed is %d\n", seed);

/* Generate the random number */

	len = (((unsigned long) gwdata->bit_length) >> 5) + 1;
	for (i = 0; i < len; i++) {
		g->n[i] = ((unsigned long) rand() << 20) +
			  ((unsigned long) rand() << 10) +
			  (unsigned long) rand();
	}
	len = 1;  g->n[0] = 256;
	g->sign = len;
	specialmodg (gwdata, g);
	gianttogw (gwdata, g, x);
}

int _tmain(int argc, _TCHAR* argv[])
{
	gwhandle gwdatastruct;
	gwhandle *gwdata = &gwdatastruct;
	gwnum	x, x2, x3, x4;
	giant	g, g2, g3, g4;
	int	i, num_squarings = 50;
	double	diff, maxdiff = 0.0;
	char	buf[200];
	int	SQUARE_ONLY = 0, CHECK_OFTEN = 1;

	OutputBothRoutine = foo;
	gwinit (gwdata);
	gwsetup (gwdata, 1.0, 2, 640, -1);

/* Alloc and init numbers */

	x = gwalloc (gwdata);
	x2 = gwalloc (gwdata);
	x3 = gwalloc (gwdata);
	x4 = gwalloc (gwdata);
	g = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	g2 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	g3 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	g4 = popg (&gwdata->gdata, ((unsigned long) gwdata->bit_length >> 4) + 20);
	gen_data (gwdata, x, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);
	gwcopy (gwdata, x, x2); gtog (g, g2);

/* Test 50 squarings */	

	gwsetnormroutine (gwdata, 0, 1, 0);	/* Enable error checking */
	for (i = 0; i < num_squarings; i++) {

		/* Test POSTFFT sometimes */
		gwstartnextfft (gwdata, !CHECK_OFTEN && (i & 3) == 2);

		/* Test gwsetaddin without and with POSTFFT set */
		if ((i == 45 || i == 46) && abs (gwdata->c) == 1)
			gwsetaddin (gwdata, -31);

		/* Test several different ways to square a number */
		if (i % 50 >= 24 && i % 50 <= 27) {
			gwfft (gwdata, x, x);
			gwfftfftmul (gwdata, x, x, x);
		} else if (i % 50 >= 32 && i % 50 <= 35) {
			gwfft (gwdata, x, x3);
			gwfftmul (gwdata, x3, x);
		} else if (i % 50 >= 40 && i % 50 <= 43) {
			gwfft (gwdata, x, x3);
			gwcopy (gwdata, x3, x4);
			gwfftfftmul (gwdata, x3, x4, x);
		} else
			gwsquare (gwdata, x);

		/* Remember maximum difference */
		diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
		if (diff > maxdiff) maxdiff = diff;

		/* Square number (and do add-in) using giants code */
		squaregi (&gwdata->gdata, g);
		if ((i == 45 || i == 46) && abs (gwdata->c) == 1) {
			iaddg (-31, g);
			gwsetaddin (gwdata, 0);
		}
		specialmodg (gwdata, g);

		/* Compare results */
		if (CHECK_OFTEN) compare (gwdata, x, g);
	}
	if (SQUARE_ONLY) goto done;

	/* Report interim results */
	if (gwdata->MAXDIFF < 1e50)
		printf ("Squares complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n",
			gw_get_maxerr (gwdata), maxdiff, gwdata->MAXDIFF,
			(int) (gwdata->MAXDIFF / maxdiff));
	else
		printf ("Squares complete. MaxErr=%.10g\n", gw_get_maxerr (gwdata));

/* Test mul by const */

	gwsetmulbyconst (gwdata, 3);
	gwsetnormroutine (gwdata, 0, 1, 1);
	gwsquare (gwdata, x);
	gwsetnormroutine (gwdata, 0, 1, 0);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g); imulg (3, g);
	specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

	gwsetmulbyconst (gwdata, -3);
	gwsetnormroutine (gwdata, 0, 1, 1);
	gwsquare (gwdata, x);
	gwsetnormroutine (gwdata, 0, 1, 0);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g); imulg (-3, g);
	specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

/* Test square carefully */

	gwfree (gwdata, x3); gwfree (gwdata, x4);
	if (abs (gwdata->c) == 1) gwsetaddin (gwdata, -42);
	gwsquare_carefully (gwdata, x);
	gwfree (gwdata, gwdata->GW_RANDOM); gwdata->GW_RANDOM = NULL;
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g);
	if (abs (gwdata->c) == 1) { iaddg (-42, g); gwsetaddin (gwdata, 0); }
	specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

/* Test gwaddquick, gwsubquick */

	x3 = gwalloc (gwdata);
	x4 = gwalloc (gwdata);
	gwadd3quick (gwdata, x, x2, x3); gtog (g, g3); addg (g2, g3);
	gwsub3quick (gwdata, x, x2, x4); gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (gwdata, g3); compare (gwdata, x3, g3);
		specialmodg (gwdata, g4); compare (gwdata, x4, g4);
	}

/* Test gwadd and gwsub */

	gwadd (gwdata, x, x); gwadd (gwdata, x, x); gwadd (gwdata, x, x);
	imulg (8, g);
	gwsub (gwdata, x3, x); subg (g3, g);
	gwadd (gwdata, x4, x); addg (g4, g);
	gwadd3 (gwdata, x3, x4, x2); gtog (g3, g2); addg (g4, g2);
	gwsub3 (gwdata, x3, x, x4); gtog (g3, g4); subg (g, g4);
	if (CHECK_OFTEN) {
		specialmodg (gwdata, g); compare (gwdata, x, g);
		specialmodg (gwdata, g2); compare (gwdata, x2, g2);
		specialmodg (gwdata, g4); compare (gwdata, x4, g4);
	}
	gwadd (gwdata, x2, x); addg (g2, g);
	gwadd (gwdata, x4, x); addg (g4, g);

/* Test gwaddsub */

	gwaddsub (gwdata, x, x2);	// compute x+x2 and x-x2
	subg (g, g2);		// x2-x
	addg (g, g);		// x+x
	addg (g2, g);		// x+x2
	negg (g2);		// x-x2
	gwaddsub4 (gwdata, x, x2, x3, x4);	// compute x+x2 and x-x2
	gtog (g, g3); addg (g2, g3);
	gtog (g, g4); subg (g2, g4);
	if (CHECK_OFTEN) {
		specialmodg (gwdata, g); compare (gwdata, x, g);
		specialmodg (gwdata, g2); compare (gwdata, x2, g2);
		specialmodg (gwdata, g3); compare (gwdata, x3, g3);
		specialmodg (gwdata, g4); compare (gwdata, x4, g4);
	}
	gwadd (gwdata, x2, x); addg (g2, g);
	gwadd (gwdata, x3, x); addg (g3, g);
	gwadd (gwdata, x4, x); addg (g4, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

/* Test gwsmalladd and gwsmallmul */

	gwsmalladd (gwdata, GWSMALLADD_MAX, x);
	dbltog (GWSMALLADD_MAX, g4); addg (g4, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);
	i = rand() % 10 + 2;
	gwsmallmul (gwdata, i, x);
	ulmulg (i, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);
	gwsmallmul (gwdata, GWSMALLMUL_MAX-1.0, x);
	ulmulg ((unsigned long) (GWSMALLMUL_MAX-1.0), g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

/* Do some multiplies to make sure that the adds and subtracts above */
/* normalized properly. */

	gwfft (gwdata, x, x);
	gwfftfftmul (gwdata, x, x, x);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	squaregi (&gwdata->gdata, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

	gwfft (gwdata, x, x2); gwcopy (gwdata, x2, x); gwfftadd3 (gwdata, x, x2, x4);
	gwfftmul (gwdata, x4, x3); addg (g3, g3); mulgi (&gwdata->gdata, g, g3); specialmodg (gwdata, g3);
	diff = fabs (gwsuminp (gwdata, x3) - gwsumout (gwdata, x3));
	if (diff > maxdiff) maxdiff = diff;
	if (CHECK_OFTEN) compare (gwdata, x3, g3);
	gwfft (gwdata, x3, x4);
	gwfftfftmul (gwdata, x4, x2, x);
	diff = fabs (gwsuminp (gwdata, x) - gwsumout (gwdata, x));
	if (diff > maxdiff) maxdiff = diff;
	mulgi (&gwdata->gdata, g3, g); specialmodg (gwdata, g);
	if (CHECK_OFTEN) compare (gwdata, x, g);

/* Do the final compare */

done:	if (!CHECK_OFTEN) compare (gwdata, x, g);
	if (gwdata->GWERROR) printf ("GWERROR set during calculations.\n");

/* Print final stats */

	if (maxdiff > gwdata->MAXDIFF) printf ("Sumout failed during test.\n");
	if (gwdata->MAXDIFF < 1e50)
		printf ("Test complete. MaxErr=%.8g, SumoutDiff=%.8g/%.8g(%d to 1)\n",
			gw_get_maxerr (gwdata), maxdiff, gwdata->MAXDIFF,
			(int) (gwdata->MAXDIFF / maxdiff));
	else
		printf ("Test complete. MaxErr=%.10g\n", gw_get_maxerr (gwdata));
	printf ("\n");

	pushg (&gwdata->gdata, 4);

	return 0;
}

