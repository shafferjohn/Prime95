/* Many of the OS-specific routines are here.  They are in this common file */
/* so that they can be included in both the command-line mprime version as */
/* well as the Mac OS X GUI version. */

/* Copyright 1995-2016 Mersenne Research, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Do some work prior to launching worker threads */
/* Windows uses this to implement boot delay. */
/* On the Mac, we use this to disable OS X Mavericks power-saving heuristics */

#ifdef __APPLE__
#ifndef COMMAND_LINE_MPRIME
#include <Foundation/NSProcessInfo.h>
id	activity;
int	activity_is_set = FALSE;
#endif
#endif

void PreLaunchCallback (
	int	launch_type)
{
#ifdef __APPLE__
#ifndef COMMAND_LINE_MPRIME
	if (DEFEAT_POWER_SAVE) {
		activity = [[NSProcessInfo processInfo] beginActivityWithOptions:NSActivityUserInitiated reason:@"Prime Hunting!"];
		activity_is_set = TRUE;
	}
#endif
#endif
}

/* Do some work after worker threads have terminated */
/* NOTE: We can't use DEFEAT_POWER_SAVE in the if statement because the global */
/* variable can be changed in the preferences page before PostLaunchCallback is called. */

void PostLaunchCallback (
	int	launch_type)
{
#ifdef __APPLE__
#ifndef COMMAND_LINE_MPRIME
	if (activity_is_set) {
		[[NSProcessInfo processInfo] endActivity:activity];
		activity_is_set = FALSE;
	}
#endif
#endif
}

/* OSes that must poll for whether the ESC key was hit do it here. */
/* We use this opportunity to perform other miscellaneous tasks that */
/* can't be done any other way. */

void stopCheckCallback (
	int	thread_num)
{
}

/* Routine to get the current load average */

double get_load_average (void)
{
#if defined (__linux__)
	char	ldavgbuf[40];
	double	load_avg;
	int	fd, count;

	fd = open ("/proc/loadavg", O_RDONLY);
	if (fd == -1) return (-1.0);
	count = read (fd, ldavgbuf, 40);
	(void) close (fd);
	if (count <= 0) return (-1.0);
	count = sscanf (ldavgbuf, "%lf", &load_avg);
	if (count < 1) return (-1.0);
	return (load_avg);
#elif defined (__FreeBSD__) || defined (__APPLE__)
	double load[3];

	if (getloadavg (load, sizeof(load)/sizeof(load[0])) < 0) return (-1.0);
	return (load[0]);
#else
	return (-1.0);
#endif
}

/* Return TRUE if we are on battery power. */

int OnBattery (void)
{
#if defined (__APPLE__)
/* The following copyright applies to the battery detection code below.  I did modify */
/* substantially as it did far more than I needed. */

/* Copyright (c) 2003 Thomas Runge (coto@core.de)
 * Mach and Darwin specific code is Copyright (c) 2006 Eric Pooch (epooch@tenon.com)
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the author nor the names of its contributors
 *    may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
	CFTypeRef blob = IOPSCopyPowerSourcesInfo();
	CFArrayRef sources = IOPSCopyPowerSourcesList(blob);

	int i, acstat;
	CFDictionaryRef pSource = NULL;
	const void *psValue;

	acstat = TRUE;

	for(i = 0; i < CFArrayGetCount(sources); i++)
	{
		pSource = IOPSGetPowerSourceDescription(blob, CFArrayGetValueAtIndex(sources, i));
		if(!pSource) break;

		psValue = (CFStringRef)CFDictionaryGetValue(pSource, CFSTR(kIOPSNameKey));

		if (CFDictionaryGetValueIfPresent(pSource, CFSTR(kIOPSIsPresentKey), &psValue) && (CFBooleanGetValue(psValue) > 0))
		{
			psValue = (CFStringRef)CFDictionaryGetValue(pSource, CFSTR(kIOPSPowerSourceStateKey));

			if (CFStringCompare(psValue,CFSTR(kIOPSBatteryPowerValue),0)==kCFCompareEqualTo)
			{
				/* We are running on a battery power source. */
				acstat = FALSE;
			}
		}
	}

	CFRelease(blob);
	CFRelease(sources);

	return(!acstat);
#elif defined (__linux__)
	FILE	*fd;
	char	buf[180];
	int	ac_state;

	ac_state = -1;
	fd = fopen ("/sys/class/power_supply/AC/online", "r");
	if (fd != NULL) {
		fscanf (fd, "%d", &ac_state);
		fclose (fd);
		return (ac_state == 0);
	}
	fd = fopen ("/proc/acpi/battery/BAT0/state", "r");
	if (fd != NULL) {
		while (fgets (buf, sizeof (buf), fd) != NULL) {
			char	*p;
			p = strstr (buf, "charging state:");
			if (p == NULL) continue;
			if (strstr (p+14, "discharging") != NULL) ac_state = 0;
			else if (strstr (p+14, "charging") != NULL) ac_state = 1;
			else if (strstr (p+14, "charged") != NULL) ac_state = 1;
		}
		fclose (fd);
	}
	return (ac_state == 0);
#else
	return (FALSE);		/* Assume we're on AC power */
#endif
}

/* The current implementation comes courtesy of Tim Wood and Dennis Gregorovic */

unsigned long physical_memory (void)
{
#if defined (__APPLE__)
	int	mib[2];
	union {
		uint32_t ui32;
		uint64_t ui64;
	} value;
	size_t	len;

	mib[0] = CTL_HW;
	mib[1] = HW_MEMSIZE;
	len = sizeof (value);
	if (sysctl (mib, 2, &value, &len, NULL, 0) < 0)
		return (1024);		/* On error, guess 1GB */
	if (len == sizeof (uint32_t))
		return (value.ui32 >> 20);
	else
		return ((unsigned long) (value.ui64 >> 20));
#elif defined (__FreeBSD__)
	int	mib[2];
	union {
		uint32_t ui32;
		uint64_t ui64;
	} value;
	size_t	len;

	mib[0] = CTL_HW;
	mib[1] = HW_PHYSMEM;
	len = sizeof (value);
	if (sysctl (mib, 2, &value, &len, NULL, 0) < 0)
		return (1024);		/* On error, guess 1GB */
	if (len == sizeof (uint32_t))
		return (value.ui32 >> 20);
	else
		return ((unsigned long) (value.ui64 >> 20));
#elif defined (__HAIKU__)
	long phys_pages;
	long page_size;

	phys_pages = sysconf(_SC_PHYS_PAGES);
	page_size = sysconf(_SC_PAGE_SIZE);

	return ((unsigned long)((double)phys_pages * (double)page_size / 1048576.0));
#else
	struct sysinfo sys_info;

	if (sysinfo(&sys_info) != 0) return (1024);  /* Guess 1GB */

	return ((unsigned long)
		  ((double) sys_info.totalram *
		   (double) sys_info.mem_unit / 1048576.0));
#endif
}

/* Return a better guess for amount of memory to use in a torture test. */
/* Caller passes in its guess for amount of memory to use, but this routine */
/* can reduce that guess based on OS-specific code that looks at amount */
/* of available physical memory. */
/* This code was written by an anonymous GIMPS user. */

unsigned long GetSuggestedMemory (unsigned long nDesiredMemory)
{
	return (nDesiredMemory);
}

int getDefaultTimeFormat (void)
{
	return (2);
}

void Sleep (
	    long	ms) 
{
#ifdef __IBMC__
	DosSleep(ms);
#else
	usleep (ms * 1000);
#endif
}

/* Clear the array of active thread handles */

void clearThreadHandleArray (void)
{
}

/* Register a thread termination.  We remove the thread handle from the */
/* list of active worker threads. */

void registerThreadTermination (void)
{
}

/* When stopping or exiting we raise the priority of all worker threads */
/* so that they can terminate in a timely fashion even if there are other */
/* CPU bound tasks running. */

void raiseAllWorkerThreadPriority (void)
{
}


/* Set priority.  Map one (prime95's lowest priority) to 20 */
/* (linux's lowest priority).  Map eight (prime95's normal priority) to */
/* 0 (linux's normal priority). */

void setOsThreadPriority (
	int	priority)		/* Priority, 1=low, 9=high */
{
#ifdef __IBMC__
	DosSetPriority(PRTYS_PROCESS,
		       (priority < 6) ? PRTYC_IDLETIME : PRTYC_REGULAR,
		       (priority == 1 || priority == 6) ? PRTYD_MINIMUM :
		       (priority == 2 || priority == 7) ? -10 :
		       (priority == 3 || priority == 8) ? 0 :
		       (priority == 4 || priority == 9) ? 10 :
		       PRTYD_MAXIMUM,
		       0);
#endif
#ifdef __linux__
	int	linux_priority, errcode;
	pid_t	thread_id;

/* I couldn't get the standard syscall0 declaration of gettid to */
/* work in my Linux distro.  Use the direct system call instead. */

	thread_id = (pid_t) syscall (__NR_gettid);

/* Set priority.  Map one (prime95's lowest priority) to 19 */
/* (linux's lowest priority).  Map eight (prime95's normal priority) to */
/* 0 (linux's normal priority). */

	linux_priority = (8 - (int) priority) * 19 / 7;
	errcode = setpriority (PRIO_PROCESS, thread_id, linux_priority);
#endif

#if defined (__APPLE__) || defined (__FreeBSD__)
	static	int	default_policy = 0;
	static	int	default_priority = 999;
	static	int	min_priority = 0;
	static	int	max_priority = 0;
	struct sched_param sp;

/* My latest MacBook Pro running OS X Mavericks changes Prime95's priority to background. */
/* This wouldn't be a problem except that the scheduler feels free to invoke SpeedStep */
/* when background processes are using 100% of the CPU time.  Thus, we now change the */
/* task policy to "trick" the scheduler into thinking we are a foreground process. */
/* Alas, this trick doesn't seem to work.  Code left here in case it might be useful one day. */

#if defined (__APPLE__)
#include <mach/mach_init.h>
#include <mach/task_policy.h>
/* APPLE BUG - the SDK included in XCode 5.0.2 comments out the task_policy_set definition in mach/task_policy.h */
	kern_return_t	task_policy_set(
                                    task_t			task,
                                    task_policy_flavor_t	flavor,
                                    task_policy_t		policy_info,
                                    mach_msg_type_number_t	count);
 
	if (IniGetInt (INI_FILE, "SetForegroundPolicy", 0)) {
		struct task_category_policy tcatpolicy;
 
		tcatpolicy.role = TASK_FOREGROUND_APPLICATION;
		task_policy_set (mach_task_self (), TASK_CATEGORY_POLICY,
				 (thread_policy_t) &tcatpolicy, TASK_CATEGORY_POLICY_COUNT);
	}
#endif

/* Get the default thread priority when a thread is first launched */

	if (default_priority == 999) {
		memset (&sp, 0, sizeof(struct sched_param));
		if (pthread_getschedparam (pthread_self(), &default_policy, &sp) >= 0) {
			min_priority = sched_get_priority_min (default_policy);
			max_priority = sched_get_priority_max (default_policy);
			default_priority = sp.sched_priority;
		} else {
			default_policy = SCHED_RR;
			min_priority = sched_get_priority_min (default_policy);
			max_priority = sched_get_priority_max (default_policy);
			default_priority = min_priority + (max_priority - min_priority) / 2;
		}
	}

/* Map one (prime95's lowest priority) to min_priority (pthread's lowest priority). */
/* Map eight (prime95's normal priority) to pthread's default priority. */

	memset (&sp, 0, sizeof(struct sched_param));
	sp.sched_priority = min_priority + (priority - 1) * (default_priority - min_priority) / 7;
	pthread_setschedparam (pthread_self(), default_policy, &sp);
#endif
}

/* Init PrimeNet communication code, make sure an internet connection is active */
/* Return false if not connected to internet */

int LoadPrimeNet (void)
{

#if defined (__linux__)
	/* Open file that will hopefully tell us if we are connected to */
	/* the Internet.  There are four possible settings for RouteRequired */
	/* 0:	Always return TRUE */
	/* 1:   Use old version 19 code */
	/* 2:   Use new code supplied by Matthew Ashton. */
	/* 99:	Default.  Use case 2 above but if cannot open /proc/net/route*/
	/*	then assume you are connected (we probably do not have read */
	/*	permission or this is a funny Linux setup). */
	int	lines = 0;
	FILE	*fd;
	char	buffer[4096];
	int	RtReq = IniSectionGetInt (INI_FILE, "PrimeNet", "RouteRequired", 99);
	if (RtReq == 0) return (TRUE);
	fd = fopen("/proc/net/route","r");
	if (fd == NULL) return (RtReq == 99);
	/* We have a readable /proc/net/route file.  Use the new check */
	/* for an Internet connection written by Matthew Ashton. However, */
	/* we still support the old style check (just in case) by setting */
	/* RouteRequired to 1. */
	if (RtReq >= 2) {
		while (fgets(buffer, sizeof(buffer), fd)) {
			int dest;
			if(sscanf(buffer, "%*s %x", &dest) == 1 && dest == 0) {
				fclose (fd);
				return (TRUE);
			}
		}
	}
	/* The old code for testing an Internet connection is below */
	else {
		(void) fgets(buffer, 199, fd);
		(void) fgets(buffer, 199, fd);
		while (!feof(fd)) {
			if (strncmp(buffer, "lo", 2)) {
				fclose(fd);
				return TRUE;
			}
			(void) fgets(buffer, 199, fd);
		}
	}
	fclose(fd);
#elif defined (__FreeBSD__) || defined (__APPLE__) || defined (__WATCOMC__) || defined (__HAIKU__)
	/* The /proc/net/route test is only really meaningfulunder linux. */
	/* For other OSes, there doesn't seem to be any meaningful test to see whether the */
	/* computer is connected to the Internet at the time using a non- */
	/* invasive test (which wouldn't, say, activate diald or ppp or */
	/* something else */
	return TRUE;
#endif
	OutputStr (COMM_THREAD_NUM, "You are not connected to the Internet.\n");
	return FALSE;
}

/* Unload the PrimeNet communication code */

void UnloadPrimeNet (void)
{
}

/* Check if a program is currently running - not implemented for OS/2 */

void checkPauseListCallback (void)
{
#ifndef __OS2__
	FILE	*fd;
	char	buf[80];

	fd = popen ("ps -A -o ucomm", "r");
	if (fd != NULL) {
		while (fgets (buf, sizeof (buf), fd) != NULL) {
			int	len = strlen (buf);
			while (len && isspace (buf[len-1])) buf[--len] = 0;
			isInPauseList (buf);
		}
		pclose (fd);
	}
#endif
}

