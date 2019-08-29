/*
 * src/bits.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines various routines to get / set bits of a IEEE floating point
 * number.  This used by the library for debugging purposes.
 */

#include <iostream>
#include <cmath>
#include <climits>
#include "bits.h"

#ifdef CRAY
#include <fp.h>
#endif

using namespace std;

#ifndef isinf
#define isinf(x) ((x) != 0.0 && (x) == (2.0*x))
#endif
#ifndef isnan
#define isnan(x) ((x) != (x))
#endif

int get_double_expn(double x) {
  if (x == 0.0)
    return INT_MIN;
  if (isinf(x) || isnan(x))
    return INT_MAX;

  double y = fabs(x);
  int i = 0;
  if (y < 1.0) {
    while (y < 1.0) {
      y *= 2.0;
      i++;
    }
    return -i;
  } else if (y >= 2.0) {
    while (y >= 2.0) {
      y *= 0.5;
      i++;
    }
    return i;
  }
  return 0;
}

#ifdef STREAMS

void print_double_info(double x) {
  if (isnan(x)) {
    cout << "NaN";
    return;
  }

  double sign = copysign(1.0, x);
  cout << (sign < 0.0 ? "- " : "+ ");
  if (isinf(x)) {
    cout << "Inf";
    return;
  }

  if (x == 0.0) {
    cout << "0";
    return;
  }

  /* Now that we handled NaNs, Infs, and Zeros, 
     x should be a normal number. */
  x = fabs(x);
  int expn = get_double_expn(x);
  double d = ldexp(1.0, expn);
  cout << expn << " ";
  while (x != 0.0) {
    if (x >= d) {
      x -= d;
      cout << '1';
    } else
      cout << '0';
    d *= 0.5;
  }
}

#endif
