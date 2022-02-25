/* Copyright 1995-2021 Mersenne Research, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Handy definitions */

#define _XOPEN_SOURCE 500
typedef int bool;
#include "common.h"
#include <strings.h>
#include <time.h>
#include <unistd.h>

/* Port number used in version numbers and result reporting. */

#ifdef __FreeBSD__
#ifdef X86_64
#define PORT	12
#else
#define PORT	6
#endif
#endif
#if defined (__EMX__) || defined (__IBMC__) || defined (__OS2__)
#define PORT	7
#endif
#ifdef __APPLE__
#ifdef X86_64
#define PORT	10
#else
#define PORT	9
#endif
#endif
#ifdef __HAIKU__
#define PORT	11
#endif

/* This controls whether we want to pause computation if the load average */
/* becomes too great.  This does not apply to OS/2. */

#if defined (__linux__) || defined (__FreeBSD__)
#define MPRIME_LOADAVG
#endif

/* The common include files */

#ifndef X86_64
#pragma pack(push)
#pragma pack(4)			// Hwloc library was (likely) built without -malign-double
#include <hwloc.h>
#pragma pack(pop)
#endif
#ifdef X86_64
#include <hwloc.h>
#endif

#include <gmp.h>

/*#define SERVER_TESTING*/
#include "cpuid.h"
#include "gwnum.h"
#include "gwbench.h"
#include "gwini.h"
#include "gwutil.h"
#include "commona.h"
#include "commonc.h"
#include "commonb.h"
#include "primenet.h"

/* Global variables */

extern int volatile KILL_MENUS;		/* TRUE if program should terminate */
extern int MENUING;			/* TRUE when main menu active */

/* Internal routines */

void sigterm_handler(int);
void main_menu (void);
void linuxContinue (char *, int, int);
void Sleep (long);
void test_user(void);
void test_welcome(void);
void test_status(void);

/* Routine definitions */

void rangeStatus (void);
void options_cpu (void);
