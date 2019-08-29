/*----------------------------------------------------------------------
| gwthread.h
|
| This file contains the headers and definitions that are used to implement
| multi-threading.  Specifically, threads, mutexes, and events.
|
| Gwnum math routines use these routines to implement multi-threaded
| multiplications.  However, you can use these routines without using
| the rest of the gwnum package if you so desire.
|
| NOTE:  MFC users may not be able to use the thread creation routines as
| because we call the C runtime library thread create rather than the
| MFC thread create routine. 
|
|  Copyright 2006 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

#ifndef _GWTHREAD_H
#define _GWTHREAD_H

/* This is a C library.  If used in a C++ program, don't let the C++ */
/* compiler mangle names. */

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
*                         Mutex and Events Routines                           *
******************************************************************************/

typedef void *gwmutex;			/* Definition of a mutex */

void gwmutex_init (gwmutex *mutex);	/* Mutex to init */
void gwmutex_lock (gwmutex *mutex);	/* Mutex to lock */
void gwmutex_unlock (gwmutex *mutex);	/* Mutex to unlock */
void gwmutex_destroy (gwmutex *mutex);	/* Mutex to destroy */

typedef void *gwevent;			/* Definition of an event handle */

void gwevent_init (gwevent *event);	/* Event to init */
#define GWEVENT_TIMED_OUT	1
#define GWEVENT_SIGNALED	2
int gwevent_wait (gwevent *event, int seconds);
void gwevent_signal (gwevent *event);	/* Event to signal */
void gwevent_reset (gwevent *event);	/* Event to reset */
void gwevent_destroy (gwevent *event);	/* Event to destroy */

/******************************************************************************
*                           Thread Routines                                   *
******************************************************************************/

typedef void *gwthread;			/* Definition of a thread handle */

void gwthread_create (gwthread *thread_id, void (*thread_proc)(void *), void *arg);
void gwthread_create_waitable (gwthread *thread_id, void (*thread_proc)(void *),	void *arg);
void gwthread_wait_for_exit (gwthread *thread_id);
void gwthread_kill (gwthread *thread_id);

#ifdef __cplusplus
}
#endif

#endif
