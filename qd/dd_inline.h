/*
 * include/dd_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2018
 *
 * Contains small functions (suitable for inlining) in the double-double
 * arithmetic package.
 */
#ifndef _DD_INLINE_H_
#define _DD_INLINE_H_

#include "math.h"
#include "inline.h"

#ifdef NO_INLINE
#define inline
#endif


/*********** Additions ************/
/* double-double = double + double */
inline dd_real dd_real::add(double a, double b) {
  double s, e;
  s = two_sum(a, b, e);
  return dd_real(s, e);
}

/* double-double + double */
inline dd_real operator+(const dd_real &a, double b) {
  double s1, s2;
  s1 = two_sum(a.hi, b, s2);
  s2 += a.lo;
  s1 = quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double + double-double */
inline dd_real operator+(const dd_real &a, const dd_real &b) {
#ifndef ACCURATE

  /* This is the less accurate version ... obeys Cray-style
     error bound. */
  double s, e;

  s = two_sum(a.hi, b.hi, e);
  e += a.lo;
  e += b.lo;
  s = quick_two_sum(s, e, e);
  return dd_real(s, e);
#else

  /* This one satisfies IEEE style error bound, 
     due to K. Briggs and W. Kahan.                   */
  double s1, s2, t1, t2;

  s1 = two_sum(a.hi, b.hi, s2);
  t1 = two_sum(a.lo, b.lo, t2);
  s2 += t1;
  s1 = quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
#endif
}

/* double + double-double */
inline dd_real operator+(double a, const dd_real &b) {
  return (b + a);
}


/*********** Self-Additions ************/
/* double-double += double */
inline dd_real &dd_real::operator+=(double a) {
  double s1, s2;
  s1 = two_sum(hi, a, s2);
  s2 += lo;
  hi = quick_two_sum(s1, s2, lo);
  return *this;
}

/* double-double += double-double */
inline dd_real &dd_real::operator+=(const dd_real &a) {
#ifndef ACCURATE
  double s, e;
  s = two_sum(hi, a.hi, e);
  e += lo;
  e += a.lo;
  hi = quick_two_sum(s, e, lo);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = two_sum(hi, a.hi, s2);
  t1 = two_sum(lo, a.lo, t2);
  s2 += t1;
  s1 = quick_two_sum(s1, s2, s2);
  s2 += t2;
  hi = quick_two_sum(s1, s2, lo);
  return *this;
#endif
}

/*********** Subtractions ************/
/* double-double = double - double */
inline dd_real dd_real::sub(double a, double b) {
  double s, e;
  s = two_diff(a, b, e);
  return dd_real(s, e);
}

/* double-double - double */
inline dd_real operator-(const dd_real &a, double b) {
  double s1, s2;
  s1 = two_diff(a.hi, b, s2);
  s2 += a.lo;
  s1 = quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/* double-double - double-double */
inline dd_real operator-(const dd_real &a, const dd_real &b) {
#ifndef ACCURATE
  double s, e;
  s = two_diff(a.hi, b.hi, e);
  e += a.lo;
  e -= b.lo;
  s = quick_two_sum(s, e, e);
  return dd_real(s, e);
#else
  double s1, s2, t1, t2;
  s1 = two_diff(a.hi, b.hi, s2);
  t1 = two_diff(a.lo, b.lo, t2);
  s2 += t1;
  s1 = quick_two_sum(s1, s2, s2);
  s2 += t2;
  s1 = quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
#endif
}

/* double - double-double */
inline dd_real operator-(double a, const dd_real &b) {
  double s1, s2;
  s1 = two_diff(a, b.hi, s2);
  s2 -= b.lo;
  s1 = quick_two_sum(s1, s2, s2);
  return dd_real(s1, s2);
}

/*********** Self-Subtractions ************/
/* double-double -= double */
inline dd_real &dd_real::operator-=(double a) {
  double s1, s2;
  s1 = two_diff(hi, a, s2);
  s2 += lo;
  hi = quick_two_sum(s1, s2, lo);
  return *this;
}

/* double-double -= double-double */
inline dd_real &dd_real::operator-=(const dd_real &a) {
#ifndef ACCURATE
  double s, e;
  s = two_diff(hi, a.hi, e);
  e += lo;
  e -= a.lo;
  hi = quick_two_sum(s, e, lo);
  return *this;
#else
  double s1, s2, t1, t2;
  s1 = two_diff(hi, a.hi, s2);
  t1 = two_diff(lo, a.lo, t2);
  s2 += t1;
  s1 = quick_two_sum(s1, s2, s2);
  s2 += t2;
  hi = quick_two_sum(s1, s2, lo);
  return *this;
#endif
}

/*********** Unary Minus ***********/
inline dd_real dd_real::operator-() const {
  return dd_real(-hi, -lo);
}

/*********** Multiplications ************/
/* double-double = double * double */
inline dd_real dd_real::mul(double a, double b) {
  double p, e;
  p = two_prod(a, b, e);
  return dd_real(p, e);
}

/* double-double * (2.0 ^ exp) */
inline dd_real ldexp(const dd_real &a, int exp) {
  return dd_real(ldexp(a.hi, exp), ldexp(a.lo, exp));
}

/* double-double * double,  where double is a power of 2. */
inline dd_real mul_pwr2(const dd_real &a, double b) {
  return dd_real(a.hi * b, a.lo * b);
}

/* double-double * double */
inline dd_real operator*(const dd_real &a, double b) {
  double p1, p2;

  p1 = two_prod(a.hi, b, p2);
  p2 += (a.lo * b);
  p1 = quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double-double * double-double */
inline dd_real operator*(const dd_real &a, const dd_real &b) {
  double p1, p2;

  p1 = two_prod(a.hi, b.hi, p2);
  p2 += a.hi * b.lo;
  p2 += a.lo * b.hi;
  p1 = quick_two_sum(p1, p2, p2);
  return dd_real(p1, p2);
}

/* double * double-double */
inline dd_real operator*(double a, const dd_real &b) {
  return (b * a);
}

/*********** Self-Multiplications ************/
/* double-double *= double */
inline dd_real &dd_real::operator*=(double a) {
  double p1, p2;
  p1 = two_prod(hi, a, p2);
  p2 += lo * a;
  hi = quick_two_sum(p1, p2, lo);
  return *this;
}

/* double-double *= double-double */
inline dd_real &dd_real::operator*=(const dd_real &a) {
  double p1, p2;
  p1 = two_prod(hi, a.hi, p2);
  p2 += a.lo * hi;
  p2 += a.hi * lo;
  hi = quick_two_sum(p1, p2, lo);
  return *this;
}

/*********** Divisions ************/
inline dd_real dd_real::div(double a, double b) {
  double q1, q2;
  double p1, p2;
  double s, e;

  q1 = a / b;

  /* Compute  a - q1 * b */
  p1 = two_prod(q1, b, p2);
  s = two_diff(a, p1, e);
  e -= p2;

  /* get next approximation */
  q2 = (s + e) / b;

  s = quick_two_sum(q1, q2, e);

  return dd_real(s, e);
}

/* double-double / double */
inline dd_real operator/(const dd_real &a, double b) {

  double q1, q2;
  double p1, p2;
  double s, e;
  dd_real r;
  
  q1 = a.hi / b;   /* approximate quotient. */

  /* Compute  this - q1 * d */
  p1 = two_prod(q1, b, p2);
  s = two_diff(a.hi, p1, e);
  e += a.lo;
  e -= p2;
  
  /* get next approximation. */
  q2 = (s + e) / b;

  /* renormalize */
  r.hi = quick_two_sum(q1, q2, r.lo);

  return r;
}

/* double-double / double-double */
inline dd_real operator/(const dd_real &a, const dd_real &b) {
  double q1, q2;
  dd_real r;

  q1 = a.hi / b.hi;  /* approximate quotient */

#ifdef SLOPPY
  double s1, s2;

  /* compute  this - q1 * dd */
  r = b * q1;
  s1 = two_diff(a.hi, r.hi, s2);
  s2 -= r.lo;
  s2 += a.lo;

  /* get next approximation */
  q2 = (s1 + s2) / b.hi;

  /* renormalize */
  r.hi = quick_two_sum(q1, q2, r.lo);
  return r;
#else
  double q3;
  r = a - q1 * b;
  
  q2 = r.hi / b.hi;
  r -= (q2 * b);

  q3 = r.hi / b.hi;

  q1 = quick_two_sum(q1, q2, q2);
  r = dd_real(q1, q2) + q3;
  return r;
#endif
}

/* double / double-double */
inline dd_real operator/(double a, const dd_real &b) {
  return dd_real(a) / b;
}

inline dd_real inv(const dd_real &a) {
  return 1.0 / a;
}

/*********** Self-Divisions ************/
/* double-double /= double */
inline dd_real &dd_real::operator/=(double a) {
  *this = *this / a;
  return *this;
}

/* double-double /= double-double */
inline dd_real &dd_real::operator/=(const dd_real &a) {
  *this = *this / a;
  return *this;
}

/********** Remainder **********/
inline dd_real drem(const dd_real &a, const dd_real &b) {
  dd_real n = nint(a / b);
  return (a - n * b);
}

inline dd_real divrem(const dd_real &a, const dd_real &b, dd_real &r) {
  dd_real n = nint(a / b);
  r = a - n * b;
  return n;
}

/*********** Squaring **********/
inline dd_real sqr(const dd_real &a) {
  double p1, p2;
  double s1, s2;
  p1 = two_sqr(a.hi, p2);
  p2 += 2.0 * a.hi * a.lo;
  p2 += a.lo * a.lo;
  s1 = quick_two_sum(p1, p2, s2);
  return dd_real(s1, s2);
}

inline dd_real dd_real::sqr(double a) {
  double p1, p2;
  p1 = two_sqr(a, p2);
  return dd_real(p1, p2);
}


/********** Exponentiation **********/
inline dd_real dd_real::operator^(int n) {
  return npwr(*this, n);
}


/*********** Assignments ************/
/* double-double = double */
inline dd_real &dd_real::operator=(double a) {
  hi = a;
  lo = 0.0;
  return *this;
}

/*********** Equality Comparisons ************/
/* double-double == double */
inline bool operator==(const dd_real &a, double b) {
  return (a.hi == b && a.lo == 0.0);
}

/* double-double == double-double */
inline bool operator==(const dd_real &a, const dd_real &b) {
  return (a.hi == b.hi && a.lo == b.lo);
}

/* double == double-double */
inline bool operator==(double a, const dd_real &b) {
  return (a == b.hi && b.lo == 0.0);
}

/*********** Greater-Than Comparisons ************/
/* double-double > double */
inline bool operator>(const dd_real &a, double b) {
  return (a.hi > b || (a.hi == b && a.lo > 0.0));
}

/* double-double > double-double */
inline bool operator>(const dd_real &a, const dd_real &b) {
  return (a.hi > b.hi || (a.hi == b.hi && a.lo > b.lo));
}

/* double > double-double */
inline bool operator>(double a, const dd_real &b) {
  return (a > b.hi || (a == b.hi && b.lo < 0.0));
}

/*********** Less-Than Comparisons ************/
/* double-double < double */
inline bool operator<(const dd_real &a, double b) {
  return (a.hi < b || (a.hi == b && a.lo < 0.0));
}

/* double-double < double-double */
inline bool operator<(const dd_real &a, const dd_real &b) {
  return (a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo));
}

/* double < double-double */
inline bool operator<(double a, const dd_real &b) {
  return (a < b.hi || (a == b.hi && b.lo > 0.0));
}

/*********** Greater-Than-Or-Equal-To Comparisons ************/
/* double-double >= double */
inline bool operator>=(const dd_real &a, double b) {
  return (a.hi > b || (a.hi == b && a.lo >= 0.0));
}

/* double-double >= double-double */
inline bool operator>=(const dd_real &a, const dd_real &b) {
  return (a.hi > b.hi || (a.hi == b.hi && a.lo >= b.lo));
}

/* double >= double-double */
inline bool operator>=(double a, const dd_real &b) {
  return (b <= a);
}

/*********** Less-Than-Or-Equal-To Comparisons ************/
/* double-double <= double */
inline bool operator<=(const dd_real &a, double b) {
  return (a.hi < b || (a.hi == b && a.lo <= 0.0));
}

/* double-double <= double-double */
inline bool operator<=(const dd_real &a, const dd_real &b) {
  return (a.hi < b.hi || (a.hi == b.hi && a.lo <= b.lo));
}

/* double <= double-double */
inline bool operator<=(double a, const dd_real &b) {
  return (b >= a);
}

/*********** Not-Equal-To Comparisons ************/
/* double-double != double */
inline bool operator!=(const dd_real &a, double b) {
  return (a.hi != b || a.lo != 0.0);
}

/* double-double != double-double */
inline bool operator!=(const dd_real &a, const dd_real &b) {
  return (a.hi != b.hi || a.lo != b.lo);
}

/* double != double-double */
inline bool operator!=(double a, const dd_real &b) {
  return (a != b.hi || b.lo != 0.0);
}

/*********** Micellaneous ************/
/*  this == 0 */
inline bool dd_real::is_zero() const {
  return (hi == 0.0);
}

/*  this == 1 */
inline bool dd_real::is_one() const {
  return (hi == 1.0 && lo == 0.0);
}

/*  this > 0 */
inline bool dd_real::is_positive() const {
  return (hi > 0.0);
}

/* this < 0 */
inline bool dd_real::is_negative() const {
  return (hi < 0.0);
}

/* Absolute value */
inline dd_real fabs(const dd_real &a) {
  return (a.hi < 0.0) ? -a : a;
}

inline dd_real abs(const dd_real &a) {
  return fabs(a);
}

/* Round to Nearest integer */
inline dd_real nint(const dd_real &a) {
  double hi = nint(a.hi);
  double lo;

  if (hi == a.hi) {
    /* High word is an integer already.  Round the low word.*/
    lo = nint(a.lo);
    
    /* Renormalize. This is needed if hi = some integer, lo = 1/2.*/
    hi = quick_two_sum(hi, lo, lo);
  } else {
    /* High word is not an integer. */
    lo = 0.0;
    if (fabs(hi-a.hi) == 0.5 && a.lo < 0.0) {
      /* There is a tie in the high word, consult the low word 
	 to break the tie. */
      hi -= 1.0;      /* NOTE: This does not cause INEXACT. */
    }
  }

  return dd_real(hi, lo);
}

inline dd_real floor(const dd_real &a) {
  double hi = floor(a.hi);
  double lo = 0.0;

  if (hi == a.hi) {
    /* High word is integer already.  Round the low word. */
    lo = floor(a.lo);
    hi = quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

inline dd_real ceil(const dd_real &a) {
  double hi = ceil(a.hi);
  double lo = 0.0;

  if (hi == a.hi) {
    /* High word is integer already.  Round the low word. */
    lo = ceil(a.lo);
    hi = quick_two_sum(hi, lo, lo);
  }

  return dd_real(hi, lo);
}

inline dd_real aint(const dd_real &a) {
  return (a.hi >= 0.0) ? floor(a) : ceil(a);
}

/* Cast to double. */
inline dd_real::operator double() const {
  return hi;
}

/* Cast to int. */
inline dd_real::operator int() const {
  return (int) hi;
}

/* Random number generator */
inline dd_real dd_real::rand() {
  return ddrand();
}

#endif /* _DD_INLINE_H_ */
