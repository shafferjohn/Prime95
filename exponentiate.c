/*----------------------------------------------------------------------
| Copyright 2020 Mersenne Research, Inc.  All rights reserved
|
| Auxiliary routine for proof generator and proof verifier
+---------------------------------------------------------------------*/

// This version is a simple left-to-right powering

// Raise a gwnum to a power
void simple_exponentiate (gwhandle *gwdata, gwnum x, uint64_t power)
{
	uint64_t current_bit;
	int	first_squaring;
	gwnum	scratch;		// Scratch gwnum for code to use

	// Handle simple cases
	if (power == 0) { dbltogw (gwdata, 1.0, x); return; }
	if (power == 1) return;

	scratch = gwalloc (gwdata);
//BUG	if (scratch == NULL) goto oom;

	// Find the topmost bit in power
	for (current_bit = 0x8000000000000000ULL; power < current_bit; current_bit >>= 1);

	// Exponentiate
	first_squaring = 1;
	for (current_bit >>= 1; current_bit; current_bit >>= 1) {
		if (first_squaring) {
			gwfft (gwdata, x, scratch);
			gwfftfftmul (gwdata, scratch, scratch, x);
			first_squaring = 0;
		} else
			gwsquare (gwdata, x);
		if (power & current_bit)
			gwfftmul (gwdata, scratch, x);
	}

	// Free allocated memory
	gwfree (gwdata, scratch);
}


// Utility routine to help with bit testing

#define is_bit_set(array,bitnum)	((array)[(bitnum)>>6] & (1ULL << ((bitnum) & 63)))

// This version uses a 3-bit windowing scheme to reduce multiplies a little

// Raise a gwnum to a power
void exponentiate_array (gwhandle *gwdata, gwnum x, uint64_t *power, int arraylen)
{
	int	current_bit, top3, first_squaring_done;
	gwnum	x1, x3, x5, x7;		// Scratch gwnums for code to use
	int	x7_is_really_x2;

	// Handle simple cases.  We've not studied the best power to switch from simple to windowed exponentiate
	while (arraylen >= 1 && power[arraylen-1] == 0) arraylen--;
	if (arraylen <= 0) { simple_exponentiate (gwdata, x, 0); return; }
	if (arraylen == 1 && *power <= 64) { simple_exponentiate (gwdata, x, *power); return; }

	// Allocate space for temporaries
	x1 = gwalloc (gwdata);
	x3 = gwalloc (gwdata);
	x5 = gwalloc (gwdata);
	x7 = gwalloc (gwdata);
//BUG	if (x1 == NULL || x3 == NULL || x5 == NULL || x7 == NULL) goto oom;

	// Pre-compute x^1, x^3, x^5.  We delay computing x^7 until we are sure we will need it.
	gwstartnextfft (gwdata, 1);
	gwfft (gwdata, x, x1);
	gwfftfftmul(gwdata, x1, x1, x7); x7_is_really_x2 = 1; gwfft(gwdata, x7, x7);
	gwfftfftmul (gwdata, x7, x1, x3); gwfft (gwdata, x3, x3);
	gwfftfftmul (gwdata, x7, x3, x5); gwfft (gwdata, x5, x5);

	// Find the topmost bit in power
	current_bit = arraylen * 64 - 1;
	while (!is_bit_set (power, current_bit)) current_bit--;

	// Special case the first 3 bits of the exponent
	// 100:  Normal would be 1,2,4,8.  We do 1,2,3,5,8		Normal is better only if x3,x5 never used
	// 101:  Normal would be 1,2,4,5.  We do 1,2,3,5
	// 110:  Normal would be 1,2,3,6,12.  We do 1,2,3,5,7,12	Normal is better only if x5,x7 never used
	// 111:  Normal would be 1,2,3,6,7.  We do 1,2,3,5,7
	top3 = 4;
	if (is_bit_set (power, current_bit - 1)) top3 += 2;
	if (is_bit_set (power, current_bit - 2)) top3 += 1;
	if (top3 == 4) { gwfftfftmul (gwdata, x3, x5, x); }
	if (top3 == 5) { gwfftfftmul (gwdata, x5, x5, x); }
	if (top3 == 6) { gwfftfftmul (gwdata, x7, x5, x7); x7_is_really_x2 = 0; gwfft (gwdata, x7, x7); gwfftfftmul (gwdata, x5, x7, x); }
	if (top3 == 7) { gwfftfftmul (gwdata, x7, x5, x7); x7_is_really_x2 = 0; gwfft (gwdata, x7, x7); gwfftfftmul (gwdata, x7, x7, x); }
	first_squaring_done = 1;	// Set flag indicating we've already done the first squaring needed in exponentiate loop

	// Sliding window exponentiate
	for (current_bit -= 3; current_bit >= 0; current_bit--) {
		if (first_squaring_done)
			first_squaring_done = 0;
		else {
			gwstartnextfft (gwdata, current_bit || is_bit_set (power, current_bit));
			gwsquare (gwdata, x);
		}
		if (!is_bit_set (power, current_bit)) continue;
		// We have a window starting with a set bit, see if we want to apply x^1, x^3, x^5, or x^7
		if (current_bit >= 2 && is_bit_set (power, current_bit - 2)) {
			gwstartnextfft (gwdata, 1);
			gwsquare (gwdata, x);
			gwsquare (gwdata, x);
			if (is_bit_set (power, current_bit - 1)) {	// Apply x^7
				if (x7_is_really_x2) { gwfftfftmul (gwdata, x7, x5, x7); x7_is_really_x2 = 0; gwfft (gwdata, x7, x7); }
				gwstartnextfft (gwdata, current_bit != 2);
				gwfftmul (gwdata, x7, x);
			} else {				// Apply x^5
				gwstartnextfft (gwdata, current_bit != 2);
				gwfftmul (gwdata, x5, x);
			}
			current_bit -= 2;
		} else {
			if (current_bit >= 1 && is_bit_set (power, current_bit - 1)) {	// Apply x^3
				gwstartnextfft (gwdata, 1);
				gwsquare (gwdata, x);
				gwstartnextfft (gwdata, current_bit != 1);
				gwfftmul (gwdata, x3, x);
				current_bit -= 1;
			} else {				// Apply x^1
				gwstartnextfft (gwdata, current_bit != 0);
				gwfftmul (gwdata, x1, x);
			}
		}
	}

	// Free allocated memory
	gwfree (gwdata, x1);
	gwfree (gwdata, x3);
	gwfree (gwdata, x5);
	gwfree (gwdata, x7);
}

// Raise a gwnum to a power
void exponentiate (gwhandle *gwdata, gwnum x, uint64_t power)
{
	exponentiate_array (gwdata, x, &power, 1);
}

// Raise a gwnum to a mpz power
void exponentiate_mpz (gwhandle *gwdata, gwnum x, mpz_t power)
{
	uint64_t *array;
	size_t	len;

	array = (uint64_t *) malloc (divide_rounding_up (mpz_sizeinbase (power, 2), 64) * sizeof (uint64_t));
	mpz_export (array, &len, -1, sizeof (uint64_t), 0, 0, power);
	exponentiate_array (gwdata, x, array, (int) len);
	free (array);
}

