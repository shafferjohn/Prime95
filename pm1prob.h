/* Copied with permission from https://github.com/preda/gpuowl/pm1 on 2020-08-11 */
/* Code courtesy of Mihai Preda */

// Returns the probability of PM1(B1,B2) success for a Mersenne 2^exponent -1 already TF'ed to factoredUpTo.
double pm1prob(unsigned exponent, int isMersenne, unsigned factoredUpTo, double B1, double B2);

// Dickman's "rho" function; rho(x) == F(1/x)
double rho(double x);
