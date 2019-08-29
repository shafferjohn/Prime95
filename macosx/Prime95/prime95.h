//
//  Prime95.h
//  Prime95
//
//  Created by George Woltman on 4/17/09.
//  Copyright 2009-2017 Mersenne Research, Inc. All rights reserved.
//
//  This file contains the include files needed to access any
//  of the common C routines
//

/* Some necessary definitions */

#define NO_GUI		0
#ifdef X86_64
#define PORT	10
#else
#define PORT	9
#endif

/* Include our common C header files */

#include <gmp.h>
#include <hwloc.h>
#include "common.h"
#include "cpuid.h"
#include "gwnum.h"
#include "gwini.h"
#include "gwbench.h"
#include "gwutil.h"
#include "commona.h"
#include "commonc.h"
#include "commonb.h"
#include "primenet.h"

/* Redefine some things with better equivalents on the Mac */

#undef ASSERTG
#define ASSERTG(a)	NSCAssert(a,@"Assertion raised")

/* Declare some functions */

void Sleep (long);
#define max(a,b)	((a)>(b)?(a):(b))
void BiggerFonts ();
void SmallerFonts ();

