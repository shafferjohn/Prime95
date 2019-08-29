/*----------------------------------------------------------------------
| gwthread.c
|
| This file contains the C routines and global variables that are used
| to implement multi-threading, mutexes, and locking.
| 
|  Copyright 2006-2009 Mersenne Research, Inc.  All rights reserved.
+---------------------------------------------------------------------*/

/* Include files */

#include <stdlib.h>
#ifdef _WIN32
#include "windows.h"
#include <process.h>
#else
#include <sys/types.h>
#include <sys/time.h>
#include <errno.h>
#include <pthread.h>
#endif

#include "gwcommon.h"
#include "gwthread.h"
#include <memory.h>

/******************************************************************************
*                         Mutex and Events Routines                           *
******************************************************************************/

void gwmutex_init (
	gwmutex	*mutex)			/* Mutex to init */
{
#ifdef _WIN32
	*mutex = (gwmutex) malloc (sizeof (CRITICAL_SECTION));
//bug	if (*mutex == NULL) do something! ;
	InitializeCriticalSection ((LPCRITICAL_SECTION) *mutex);
#else
	*mutex = (gwmutex) malloc (sizeof (pthread_mutex_t));
//bug	if (*mutex == NULL) do something! ;
	pthread_mutex_init ((pthread_mutex_t *) *mutex, NULL);
#endif
}

void gwmutex_lock (
	gwmutex	*mutex)			/* Mutex to lock */
{
#ifdef _WIN32
	EnterCriticalSection ((LPCRITICAL_SECTION) *mutex);
#else
	pthread_mutex_lock ((pthread_mutex_t *) *mutex);
#endif
}

void gwmutex_unlock (
	gwmutex	*mutex)			/* Mutex to unlock */
{
#ifdef _WIN32
	LeaveCriticalSection ((LPCRITICAL_SECTION) *mutex);
#else
	pthread_mutex_unlock ((pthread_mutex_t *) *mutex);
#endif
}


void gwmutex_destroy (
	gwmutex	*mutex)			/* Mutex to destroy */
{
	if (mutex != NULL) {
#ifdef _WIN32
		DeleteCriticalSection ((LPCRITICAL_SECTION) *mutex);
#else
		pthread_mutex_destroy ((pthread_mutex_t *) *mutex);
#endif
		free (*mutex);
		*mutex = NULL;
	}
}

/* Data structures for implementing events */

#ifdef _WIN32
#else
struct event_data {
	pthread_mutex_t event_mutex;
	pthread_cond_t event_cond;
	int	signalled_state;
	int	threads_waiting;
	int	signalled_count;
};
#endif

/* Event routines */

void gwevent_init (
	gwevent	*event)			/* Event to init */
{
#ifdef _WIN32
 	* (HANDLE *) event = CreateEvent (
				NULL,		// No security stuff
				TRUE,		// Manual reset event
				FALSE,		// Initially not signaled
				NULL);		// No name
	ASSERTG (*(HANDLE *) event != 0);
#else
	struct event_data *e;

	e = (struct event_data *) malloc (sizeof (struct event_data));
//bug	if (e == NULL) do something! ;
	pthread_mutex_init (&e->event_mutex, NULL);
	pthread_cond_init (&e->event_cond, NULL);
	e->signalled_state = FALSE;
	e->threads_waiting = 0;
	e->signalled_count = 0;
	*event = (gwevent) e;
#endif
}

int gwevent_wait (
	gwevent	*event,			/* Event to wait on */
	int	seconds)		/* Seconds until timeout */
{
#ifdef _WIN32
	DWORD	rc;

	if (seconds == 0)
		rc = WaitForSingleObject (
			* (HANDLE *) event,	//handle to event
			INFINITE);		//wait forever
	else
		rc = WaitForSingleObject (
			* (HANDLE *) event,	//handle to event
			seconds * 1000);	//timeout value in ms.
	if (rc == WAIT_TIMEOUT) return (GWEVENT_TIMED_OUT);
	ASSERTG (rc == 0);
	return (GWEVENT_SIGNALED);
#else
	struct event_data *e;
	struct timeval now;
	struct timespec timeout;
	int	rc;

	e = (struct event_data *) *event;
	rc = 0;

	if (seconds) {
		gettimeofday (&now, NULL);
		timeout.tv_sec = now.tv_sec + seconds;
		timeout.tv_nsec = now.tv_usec * 1000;
	}

/* Obtain lock before reading event structure.  If we're in the signalled */
/* state skip the wait. */ 

	pthread_mutex_lock (&e->event_mutex);
	if (!e->signalled_state) {

/* Loop until we timeout or the signalled_count changes.  According to */
/* pthreads documentation the conditional wait can return even though no */
/* signal took place! */
		
		e->threads_waiting++;
		for ( ; ; ) {
			int	signalled_count;
			signalled_count = e->signalled_count;
			if (seconds) {
				rc = pthread_cond_timedwait (&e->event_cond,
							     &e->event_mutex,
							     &timeout);
				if (rc == ETIMEDOUT) break;
			} else
				rc = pthread_cond_wait (&e->event_cond,
							&e->event_mutex);
			if (signalled_count != e->signalled_count) break;
		}
		e->threads_waiting--;
	}
	pthread_mutex_unlock (&e->event_mutex);

/* Return timed-out code or event-signalled code */

	if (rc == ETIMEDOUT) return (GWEVENT_TIMED_OUT);
	ASSERTG (rc == 0);
	return (GWEVENT_SIGNALED);
#endif
}

void gwevent_signal (
	gwevent	*event)			/* Event to signal */
{
#ifdef _WIN32
	BOOL	rc;

	rc = SetEvent (* (HANDLE *) event);	//handle to event
	ASSERTG (rc);
#else
	struct event_data *e;

	e = (struct event_data *) *event;
	pthread_mutex_lock (&e->event_mutex);
	if (!e->signalled_state) {
		e->signalled_state = TRUE;
		e->signalled_count++;
		if (e->threads_waiting)
			pthread_cond_broadcast (&e->event_cond);
	}
	pthread_mutex_unlock (&e->event_mutex);
#endif
}

void gwevent_reset (
	gwevent	*event)			/* Event to reset */
{
#ifdef _WIN32
	BOOL	rc;

	rc = ResetEvent (* (HANDLE *) event);	//handle to event
	ASSERTG (rc);
#else
	struct event_data *e;

	e = (struct event_data *) *event;
	pthread_mutex_lock (&e->event_mutex);
	e->signalled_state = FALSE;
	pthread_mutex_unlock (&e->event_mutex);
#endif
}


void gwevent_destroy (
	gwevent	*event)			/* Event to destroy */
{
#ifdef _WIN32
	BOOL	rc;

	rc = CloseHandle (* (HANDLE *) event);	//handle to event
	ASSERTG (rc);
#else
	struct event_data *e;

	e = (struct event_data *) *event;
	pthread_mutex_destroy (&e->event_mutex);
	pthread_cond_destroy (&e->event_cond);
	free (e);
	*event = NULL;
#endif
}

/******************************************************************************
*                           Thread Routines                                   *
******************************************************************************/

/* Dummy structure and routine to pass to _beginthread and _beginthreadex */
/* so that we use the proper calling convention. */

struct threaddata {
	void	(*proc)(void *);/* Thread routine to call */
	void	*arg;		/* Argument to pass to thread routine */
};

#ifdef _WIN32
void __cdecl ThreadStarter (
	void	*arg)
{
	struct threaddata *td;

	td = (struct threaddata *) arg;
	(*td->proc)(td->arg);
	free (arg);
}

unsigned __stdcall ThreadExStarter (
	void	*arg)
{
	struct threaddata *td;

	td = (struct threaddata *) arg;
	(*td->proc)(td->arg);
	free (arg);
	return (0);
}
#else
void *ThreadStarter (
	void	*arg)
{
	struct threaddata *td;

	td = (struct threaddata *) arg;
	(*td->proc)(td->arg);
	free (arg);
	return (NULL);
}

#endif

/* Create a thread.  The thread ends when thread_proc returns (no cleanup */
/* necessary). There is no way for another thread to wait for the thread */
/* to end. */

void gwthread_create (
	gwthread *thread_id,
	void	(*thread_proc)(void *),
	void	*arg)
{
	struct threaddata *td;
	td = (struct threaddata *) malloc (sizeof (struct threaddata));
//bug - check for mem error
	td->proc = thread_proc;
	td->arg = arg;

#ifdef _WIN32
	* (HANDLE *) thread_id = (HANDLE)
		_beginthread (
			&ThreadStarter, // Routine to execute
			0,		// Standard stack size
			td);		// Routine's argument
#else
	{
	pthread_attr_t attr;
	pthread_t tid;

	pthread_attr_init (&attr);
	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_DETACHED);
	pthread_create (&tid, &attr, &ThreadStarter, td);
	pthread_attr_destroy (&attr);
	*thread_id = (gwthread) tid;
	}
#endif
}

/* Create a thread.  The thread_id remains valid after the thread_proc */
/* returns.  Another thread must wait for this thread by calling */
/* gwthread_wait_for_exit. */

void gwthread_create_waitable (
	gwthread *thread_id,
	void	(*thread_proc)(void *),
	void	*arg)
{
	unsigned int thread_num;	// Dummy argument
	struct threaddata *td;

	td = (struct threaddata *) malloc (sizeof (struct threaddata));
//bug - check for mem error
	td->proc = thread_proc;
	td->arg = arg;

#ifdef _WIN32
	* (HANDLE *) thread_id = (HANDLE)
		_beginthreadex (
			NULL,		// No security stuff
			0,		// Standard stack size
			&ThreadExStarter, // Routine to execute
			td,		// Routine's argument
			0,		// Start thread in running state
			&thread_num);	// Returned threadID - unused
#else
	{
	pthread_attr_t attr;
	pthread_t tid;

	pthread_attr_init (&attr);
	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
	pthread_create (&tid, &attr, &ThreadStarter, td);
	pthread_attr_destroy (&attr);
	*thread_id = (gwthread) tid;
	}
#endif
}

/* Wait for a thread to end then cleanup.  Thread_id is no longer valid */
/* when this routine returns. */

void gwthread_wait_for_exit (
	gwthread *thread_id)
{
#ifdef _WIN32
	WaitForSingleObject (* (HANDLE *) thread_id, INFINITE);
	CloseHandle (* (HANDLE *) thread_id);
#else
	void	*returned_status;
	pthread_join ((pthread_t) *thread_id, &returned_status);
#endif
}


void gwthread_kill (gwthread *thread_id)
{
/*BUG	call QueueUserAPC and catch/throw or setjmp/longjmp ????? */
/*BUG	or is this impossible to program cleanly ????? */
}
