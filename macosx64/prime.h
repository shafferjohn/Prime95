/* Copyright 1995-2020 Mersenne Research, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Handy definitions */

#include "common.h"

/* Port number used in version numbers and result reporting. */

#ifdef __linux__
#ifdef X86_64
#define PORT	8
#else
#define PORT	2
#endif
#endif
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

/* Handle differences between Windows and Linux runtime libraries */

#define _commit(f)	fsync(f)
#define _open		open
#define _close		close
#define _read		read
#define _write		write
#define _lseek		lseek
#define _lseeki64	lseek
#define _fseeki64	fseek
#define _chsize_s	ftruncate
#define _unlink		unlink
#define _creat		creat
#define _chdir		chdir
#define closesocket	close
#define IsCharAlphaNumeric(c) isalnum(c)
#define _stricmp	strcasecmp

#ifndef __WATCOMC__
#define _O_APPEND	O_APPEND
#define _O_RDONLY	O_RDONLY
#define _O_WRONLY	O_WRONLY
#define _O_RDWR		O_RDWR
#define _O_CREAT	O_CREAT
#define _O_TRUNC	O_TRUNC
#define _O_BINARY 	0
#define _O_TEXT		0
#endif

/* Handle differences between Windows and OS/2 runtime libraries */

#ifdef __IBMC__
#define stricmp(x,y)  stricmp(x,y)
#define _commit(f)    /* no commit/fsync on OS/2 */
#define _ftime        _ftime
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

#include <time.h>
/*#define SERVER_TESTING*/
extern int NO_GUI;
#include "common.h"
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
