/*----------------------------------------------------------------------
| gwcommon.h
|
| This file contains handy #defines that I use in all my projects
| 
|  Copyright 2005-2016 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWCOMMON_H
#define _GWCOMMON_H

#ifndef TRUE
#define TRUE	1
#endif
#ifndef FALSE
#define FALSE	0
#endif

/* In many cases where the C code is interfacing with the assembly code */
/* we must declare variables that are exactly 32-bits wide.  This is the */
/* portable way to do this, as the linux x86-64 C compiler defines the */
/* long data type as 64 bits.  We also use portable definitions for */
/* values that can be either an integer or a pointer. */

#ifdef _MSC_VER
typedef __int32			int32_t;
typedef unsigned __int32	uint32_t;
typedef __int64			int64_t;
typedef unsigned __int64	uint64_t;
#ifdef _WIN64
typedef __int64			intptr_t;
typedef unsigned __int64	uintptr_t;
#else
typedef int			intptr_t;
typedef unsigned int		uintptr_t;
#endif
#else
#include "inttypes.h"
#endif

/* Define the ASSERT macro I use while debugging */

#include <assert.h>
#ifdef GDEBUG
#define ASSERTG		assert
#else
#define ASSERTG(a)
#endif
#define GWASSERT  	assert

#endif
