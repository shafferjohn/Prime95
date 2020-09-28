/*----------------------------------------------------------------------
| Copyright 1995-2020 Mersenne Research, Inc.  All rights reserved
+---------------------------------------------------------------------*/

//#define SERVER_TESTING

/* These macros produce the security codes used for sending in results */
/* The security.h file is not required to compile a trial version of */
/* this program. */

#include "security.h"
#ifndef SEC1
#define SEC1(p)			0UL
#define SEC2(p,hi,lo,u,e)	0UL
#define SEC3(p)			0UL
#define SEC4(p)			0UL
#define SEC5(p,b1,b2)		0UL
#endif

/* Common routines */

#define ALL_WORKERS 8888	// Special value for first arg of LaunchWorkerThreads
int LaunchWorkerThreads (int, int);
int LaunchTortureTest (unsigned long, int);
int LaunchBench (int);
int LaunchAdvancedTime (unsigned long, unsigned long);

/* Callback routines */

#define LD_CONTINUE	0		/* Call primeContinue */
#define LD_TIME		1		/* Call primeTime */
#define LD_BENCH	2		/* Call primeBench */
#define LD_TORTURE	3		/* Call selfTest */
void PreLaunchCallback (int launch_type);	/* Implemented for each OS */
void PostLaunchCallback (int launch_type);	/* Implemented for each OS */
void stopCheckCallback (int thread_num);
void checkPauseListCallback (void);

/* Routines and structures to set thread priority.  Since we may have */
/* different rules doing normal work vs. benchmarking vs. torture testing, */
/* we send that info to the SetPriority routine to use as it sees fit. */

#define SET_PRIORITY_NORMAL_WORK	1
#define SET_PRIORITY_BENCHMARKING	2
#define SET_PRIORITY_TORTURE		3
#define SET_PRIORITY_TIME		4
#define SET_PRIORITY_QA			5
#define SET_PRIORITY_BUSY_LOOP		6
struct PriorityInfo {
 	int	type;			/* Type defined above */
	int	worker_num;		/* Worker number -- output informative messages to this worker's window */
	int	verbose_flag;		/* Output affinity messages to the worker window */
	int	aux_thread_num;		/* Set when gwnum launches auxiliary threads */
	int	aux_hyperthread;	/* Set when gwnum launches auxiliary hyperthreads */
	union {
		struct {		/* Normal work info */
			int	normal_work_hyperthreads;	/* Number of hyperthreads to be assigned to same core */
		};
		struct {		/* Torture test info */
			int	torture_num_workers;		/* Total number of torture test worker windows */
			int	torture_threads_per_test;	/* Number of threads per torture test (usually one) */
		};
		struct {		/* Advanced/Time info */
			int	time_hyperthreads;		/* Number of hyperthreads to be assigned to same core */
		};
		struct {		/* Benchmark info */
			int	bench_base_cpu_num;		/* First CPU core to set affinity to */
			int	bench_hyperthreads;		/* Number of hyperthreads to be assigned to same core */
		};
		struct {		/* Busy loop info */
			int	busy_loop_cpu;			/* CPU core to keep busy */
		};
	};
};
void SetPriority (struct PriorityInfo *);

/* Internal routines that do the real work */

void create_worker_windows (int num_threads);
void Launcher (void *arg);
void LauncherDispatch (void *arg);
int primeContinue (int);
void tortureTestDefaultSizes (int torture_type, int num_threads, int *minfft, int *maxfft);
int tortureTest (int, int);
int selfTest (int, struct PriorityInfo *, struct work_unit *);
int primeTime (int, unsigned long, unsigned long);
int primeBench (int, int);
int primeFactor (int, struct PriorityInfo *, struct work_unit *, unsigned int);
int prime (int, struct PriorityInfo *, struct work_unit *, int);
int prp (int, struct PriorityInfo *, struct work_unit *, int);
int cert (int, struct PriorityInfo *, struct work_unit *, int);
int ecm (int, struct PriorityInfo *, struct work_unit *);
int pminus1 (int, struct PriorityInfo *, struct work_unit *);
int pfactor (int, struct PriorityInfo *, struct work_unit *);
double guess_pminus1_probability (struct work_unit *w);
void autoBench (void);

/* Utility routines */

int isKnownMersennePrime (unsigned long);
void makestr (unsigned long, unsigned long, unsigned long, char *);

/* Stop routines */

/* Reasons returned for why we stopped processing a work unit */
#define STOP_ESCAPE		1	/* User hit escape key.  Stopping all worker threads. */
#define STOP_OUT_OF_MEM		2	/* Fatal out-of-memory error */
#define STOP_FILE_IO_ERROR	3	/* Fatal file I/O error */
#define STOP_FATAL_ERROR	4	/* Other fatal error */
#define STOP_ABORT		5	/* Abort current work unit */
#define STOP_WORKER		6	/* This worker thread was stopped or was never started */
#define STOP_PAUSE		7	/* This worker must pause because a program in the PauseWhileRunning is running. */
#define STOP_RETRY_LATER	8	/* Abort work unit, but try again later */
#define STOP_WORK_UNIT_COMPLETE	50	/* Work unit is done! */
#define STOP_PRIORITY_WORK	51	/* Priority work, restart thread */
#define STOP_BATTERY		52	/* On battery - pause */
#define STOP_AUTOBENCH		53	/* Stop worker for a little while to run auto-benchmarks */
#define STOP_REREAD_INI		100	/* Reread prime.ini because a during/else time period has changed */
#define STOP_RESTART		101	/* Important INI option changed */
#define STOP_MEM_CHANGED	102	/* Day/night memory change */
#define STOP_NOT_ENOUGH_MEM	103	/* Not enough memory for P-1 stage 2 */

EXTERNC int stopCheck (int);
void stop_workers_for_escape (void);
void stop_workers_for_restart (void);
void stop_workers_for_add_files (void);
void stop_workers_for_reread_ini (void);
#define RESTART_ALL			0xFFFF
#define RESTART_WORK_AVAILABLE		0x0001
#define RESTART_USER_START		0x0002
#define RESTART_END_PAUSE		0x0004
#define RESTART_MEM_WAIT		0x0008
#define RESTART_BATTERY			0x0010
#define RESTART_LOADAVG			0x0020
void restart_waiting_workers (int);
void restart_one_waiting_worker (int, int);
void stop_worker_for_abort (int);

/* Routines dealing with day/night memory settings */

void mem_settings_have_changed (void);
unsigned long max_mem (int);
int avail_mem (int, unsigned long, unsigned long, unsigned int *);
int set_memory_usage (int, int, unsigned long);

/* Handy macros to help in calling memory routines.  Macros tell */
/* us how many gwnums fit in given megabytes AND how many megabytes */
/* are used by a given number of gwnums.  We have to be real careful */
/* as machines with more than 4GB are becoming commonplace.  We assume */
/* code and other data will add 1MB or 2MB to our working set. */

#define cvt_mem_to_gwnums(g,m) ((unsigned long)(((double)(m)*1048576.0-1000000.0-(double)gwmemused(g))/(double)gwnum_size(g)))
#define cvt_gwnums_to_mem(g,n) ((unsigned long)(((double)gwmemused(g)+2000000.0+(double)(n)*(double)gwnum_size(g))/1048576.0))
#define cvt_mem_to_estimated_gwnums(m,k,b,n,c) ((unsigned long)(((double)(m)*1048576.0-1000000.0-(double)gwmap_to_memused(k,b,n,c))/(double)gwmap_to_estimated_size(k,b,n,c)))

/* battery routines */

void run_on_battery_changed (void);
int OnBattery (void);			/* Implemented for each OS */
void test_battery (void);

/* priority work routines */

void check_for_priority_work (void);
void stop_worker_for_advanced_test (int);

/* starting/stopping individual worker threads routines */

void start_one_worker (int);
void stop_one_worker (int);
unsigned int active_workers_count (void);

/* "pause while running" routines */

extern char **PAUSE_WHILE_RUNNING;	/* An array of program names that, */
					/* if running, prime95 should pause. */
extern int PAUSE_WHILE_RUNNING_FREQ;	/* How often prime95 should check */
					/* the pause-while-running list */
void read_pause_info (void);
void checkPauseWhileRunning (void);
void implement_pause (int thread_num);
int is_LowMemWhileRunning_active (int thread_num);
void isInPauseList (char *program_name);

/* Load average routines */

void read_load_average_info (void);
double get_load_average (void);		/* Implemented by each OS */
void checkLoadAverage (void);
void implement_loadavg (int thread_num);

/* throttle routines */

void start_throttle_timer (void);
void stop_throttle_timer (void);
int handleThrottleTimerEvent (void);
void implementThrottle (int thread_num);

/* Routines called by common routines */

void clearThreadHandleArray (void);
void setOsThreadPriority (int);
void registerThreadTermination (void);
void raiseAllWorkerThreadPriority (void);
void flashWindowAndBeep (void);
int primeSieveTest (int);
int setN (gwhandle *, int, struct work_unit *, giant *);
int ecm_QA (int, struct PriorityInfo *);
int pminus1_QA (int, struct PriorityInfo *);
int test_randomly (int, struct PriorityInfo *);
int test_all_impl (int, struct PriorityInfo *);

/* Messages */

#define BENCH_SPEED  "The CPU speed in Options/CPU may not be correct.  An incorrect value will result in inaccurate timings.  Are you sure this is the correct CPU speed value?"
