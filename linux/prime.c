/* Copyright 1995-2020 Mersenne Research, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Include files needed by all ports */
#include <ctype.h>
#include <fcntl.h>
#include <math.h>
#include <memory.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "prime.h"

/* Required Linux files */
#ifdef __linux__
#include <dirent.h>
#include <unistd.h>
#include <linux/unistd.h>
#include <asm/unistd.h>
#define __USE_GNU
#include <sched.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#endif

/* Required Mac OS X files */
#ifdef __APPLE__
#include <dirent.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <sys/time.h>
#include <CoreFoundation/CoreFoundation.h>
#include <IOKit/ps/IOPowerSources.h>
#include <IOKit/ps/IOPSKeys.h>
#endif

/* Required FreeBSD files */
#ifdef __FreeBSD__
#include <dirent.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <sys/param.h>
#include <sys/cpuset.h>
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#include <sys/time.h>
#endif

/* Required OS/2 header files */
#ifdef __IBMC__
#define INCL_DOS
#define INCL_DOSPROFILE
#include <os2.h>
#include <direct.h>
#include <io.h>
#include <process.h>
typedef int pid_t;
#include "dosqss.h"
#endif

/* Required Watcom C (is this defunct DOS port???) header files */
#ifdef __WATCOMC__
#include <dirent.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif

/* Required Haiku files */
#ifdef __HAIKU__
#include <dirent.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/time.h>
#endif

/* Globals */

int volatile KILL_MENUS = 0;
int NO_GUI = 1;
int VERBOSE = 0;
int MENUING = 0;
char pidfile[256];

/* Common code */

#include "gwutil.h"
#include "cJSON.h"
#include "cJSON.c"
#include "pm1prob.h"
#include "pm1prob.c"
#include "md5.c"
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "primenet.c"
#include "proof_upload.c"
#include "proof_getdata.c"
#include "gwtest.c"

/* Signal handlers */

void sigterm_handler(int signo)
{
	int	i;

	stop_workers_for_escape ();	// Gracefully stop any active worker threads
	for (i = 0; WORKER_THREADS_STOPPING && i < 100; i++) Sleep (50);  // Give the worker threads some time to stop gracefully
	if (signo != SIGINT) {
		KILL_MENUS = TRUE;	// Set flag so we exit the menus
		fclose (stdin);		// Makes fgets in menu.c return.  Thus, mprime will terminate rather than waiting for a menu choice.
	}
	(void)signal(signo, sigterm_handler);
}

/* Main entry point! */

int main (
	int	argc,
	char	*argv[])
{
	char	buf[256];
	int	named_ini_files = -1;
	int	contact_server = 0;
	int	torture_test = 0;
	int	i, nice_level;
	int	pid_fd;
	char	*p;

/* catch termination signals */

	(void)signal(SIGTERM, sigterm_handler);

	/* SIGINT handling.  See discussion in http://www.mersenneforum.org/showthread.php?t=21672 */
	if (signal(SIGINT, sigterm_handler) == SIG_DFL) {
		/* If we are run asynchronously, keep SIGINT ignored. */
		/* (e.g. `sh -c 'mprime & other_command'`) */
		(void)signal(SIGHUP, SIG_DFL);
	}

	/* SIGHUP handling.  See discussions in http://www.mersenneforum.org/showthread.php?t=21496 */
	/* and http://www.mersenneforum.org/showthread.php?t=21672 */
	if (signal(SIGHUP, sigterm_handler) == SIG_DFL) {
		/* If we are run through 'nohup', keep SIGHUP ignored. */
		(void)signal(SIGHUP, SIG_DFL);
	}

	/* This code was suggested for handling SIGPIPE.  Since 99.9% of mprime installs use */
	/* LIBCURL rather than the archaic sockets code, we ought to be OK with sigterm_handler. */
	/* This assumes LIBCURL does not depend on the SIGPIPE signal.  That said, perhaps it */
	/* is OK to stick with ignoring or using the default signal handler.  See the mersenneforum */
	/* threads mentioned above for more discussion. */
//#if !defined(SO_NOSIGPIPE) && !defined(MSG_NOSIGNAL)
//	(void)signal(SIGPIPE, SIG_IGN);
//#else
//	/* SIGPIPE might be generated when user pipes stdout to another */
//	/* program (e.g. `mprime | tee log.txt`), and that "tee" dies. */
//	/* Handle this to terminate when we're sure we never receive */
//	/* SIGPIPE from sockets. */
//	(void)signal(SIGPIPE, sigterm_handler);
//#endif

/* No buffering of output */

	setvbuf (stdout, NULL, _IONBF, 0);

/* Change to the executable's directory */
/* NOTE:  This only changes the working directory if the user typed */
/* in a full path to the executable (as opposed to finding it on the PATH) */

	strcpy (buf, argv[0]);
	p = strrchr (buf, '/');
	if (p != NULL) {
		*p = 0;
		(void) _chdir (buf);
	}

/* Default pid file is mprime.pid */

	strcpy (pidfile, "mprime.pid");

/* Initialize gwnum call back routines.  Using callback routines lets the */
/* gwnum library have a nice clean interface for users that do not need */
/* additional functionality that only prime95 uses. */

	StopCheckRoutine = stopCheck;
	OutputBothRoutine = OutputBoth;

/* Process command line switches */

	for (i = 1; i < argc; i++) {
		p = argv[i];

		if (*p++ != '-') break;
		switch (*p++) {

/* Accept a -A switch indicating an alternate set of INI files */
/* are to be used. */

		case 'A':
		case 'a':
			named_ini_files = 0;
			while (isspace (*p)) p++;
			while (isdigit (*p)) {
				named_ini_files = named_ini_files * 10 + (*p - '0');
				p++;
			}
			break;

/* -C - contact the server now, then exit */

		case 'C':
		case 'c':
			contact_server = 1;
			VERBOSE = TRUE;
			NO_GUI = FALSE;
			break;
			
/* -D - debug */

		case 'D':
		case 'd':
			VERBOSE = TRUE;
			NO_GUI = FALSE;
			break;

/* -H - help */

		case 'H':
		case 'h':
		case '?':
			goto usage;

/* -M - Menu */

		case 'M':
		case 'm':
			MENUING = 1;
			NO_GUI = FALSE;
			break;

/* -P - use a different pidfile */

		case 'P':
		case 'p':
			strcpy (pidfile, p);
			break; 

/* -S - status */

		case 'S':
		case 's':
			MENUING = 2;
			NO_GUI = FALSE;
			break;
		  

/* -T - Torture test */

		case 'T':
		case 't':
			torture_test = TRUE;
			break;

/* -V - version number */

		case 'V':
		case 'v':
			generate_application_string (buf);
			printf ("Mersenne Prime Test Program: %s\n", buf);
			return (0); 

/* -W - use a different working directory */

		case 'W':
		case 'w':
			(void) _chdir (p);
			break; 

/* Otherwise unknown switch */

		default:
			printf ("Invalid switch\n");
			goto usage;
		}
	}

/* Create the pidfile */

	pid_fd = _open (pidfile, _O_CREAT | _O_TRUNC | _O_WRONLY | _O_TEXT, CREATE_FILE_ACCESS);
	if (pid_fd >= 0) {
		sprintf (buf, "%d", (int) getpid ());
		_write (pid_fd, buf, strlen (buf));
		_close (pid_fd);
	}

/* Determine the names of the INI files, read them, do other initialization. */
/* Skip the comm code initialization if we are just displaying the status */
/* or running a torture test */

	nameAndReadIniFiles (named_ini_files);
	if (MENUING != 2 && !torture_test) initCommCode ();

/* If not running a torture test, set the program to nice priority. */
/* Technically, this is not necessary since worker threads are set to */
/* the lowest possible priority.  However, sysadmins might be alarmed */
/* to see a CPU intensive program not running at nice priority when */
/* executing a ps command. */

#if defined (__linux__) || defined (__APPLE__) || defined (__FreeBSD__)
	/* Linux/FreeBSD ranges from -20 to +19, lower values give more favorable scheduling */
	nice_level = IniGetInt (INI_FILE, "Nice", 10);
	if (!torture_test && nice_level) {
		setpriority (PRIO_PROCESS, 0, nice_level);
	}
#endif

/* If running the torture test, do so now. */

	if (torture_test) {
		int	num_threads;

		VERBOSE = TRUE;
		NO_GUI = FALSE;
		num_threads = IniGetInt (INI_FILE, "TortureThreads", NUM_CPUS * CPU_HYPERTHREADS);
		if (num_threads < 1) num_threads = 1;
		if (num_threads > MAX_NUM_WORKER_THREADS) num_threads = MAX_NUM_WORKER_THREADS;
 		LaunchTortureTest (num_threads, TRUE);
	}

/* If this is a stress tester, then turn on menuing. */

	else if (IniGetInt (INI_FILE, "StressTester", 99) == 1) {
		MENUING = 1;
		VERBOSE = TRUE;
		NO_GUI = FALSE;
		main_menu ();
	}

/* On first run, get user id before contacting server */
/* for a work assignment.  To make first time user more comfortable, we will */
/* display data to the screen, rather than running silently. */

	else if (IniGetInt (INI_FILE, "StressTester", 99) == 99) {
		VERBOSE = TRUE;
		NO_GUI = FALSE;
		test_welcome ();
	}

/* If we are to contact the server, do so now.  This option lets the */
/* user create a batch file that contacts the server at regular intervals */
/* or when the ISP is contacted, etc. */

	else if (contact_server) {
		do_manual_comm_now ();
		while (COMMUNICATION_THREAD) Sleep (50);
	}

/* Bring up the main menu */

	else if (MENUING == 1)
		main_menu ();
	else if (MENUING == 2)
		test_status();

/* Continue testing, return when worker threads exit. */

	else {
		linuxContinue ("Another mprime is already running!\n", ALL_WORKERS, TRUE);
	}

/* Write the worktodo file in case the WELL_BEHAVED_WORK flag caused us */
/* to delay writing the file. */

	writeWorkToDoFile (TRUE);

/* Delete the pidfile */

	_unlink (pidfile);

/* All done */

	return (0);

/* Invalid args message */

usage:	printf ("Usage: mprime [-cdhmstv] [-aN] [-wDIR] [-pPIDFILE]\n");
	printf ("-c\tContact the PrimeNet server, then exit.\n");
	printf ("-d\tPrint detailed information to stdout.\n");
	printf ("-h\tPrint this.\n");
	printf ("-m\tMenu to configure mprime.\n");
	printf ("-s\tDisplay status.\n");
	printf ("-t\tRun the torture test.\n");
	printf ("-v\tPrint the version number.\n");
	printf ("-aN\tUse an alternate set of INI and output files (obsolete).\n");
	printf ("-wDIR\tRun from a different working directory.\n");
	printf ("-pPIDFILE\tFilename for the PID file.  Default is mprime.pid.\n");
	printf ("\n");
	return (1);
}

/* Create an MDI output window for the thread -- unless we created */
/* one earlier and the user has not closed it */

void create_window (
	int	thread_num)
{
}

/* Set the title prefix for this MDI window - only called once */

void base_title (int thread_num, const char *str)
{
}

/* Put a title on the MDI window */

void title (int thread_num, const char *msg)
{
}

void flashWindowAndBeep (void)
{
	printf ("\007");
}

void RealOutputStr (int thread_num, const char *buf)
{
static	int	last_char_out_was_newline = TRUE;
	if (VERBOSE || MENUING) {
		if (last_char_out_was_newline && !(MERGE_WINDOWS & MERGE_NO_PREFIX)) {
			if (thread_num == MAIN_THREAD_NUM)
				printf ("[Main thread");
			else if (thread_num == COMM_THREAD_NUM)
				printf ("[Comm thread");
			else if (NUM_WORKER_THREADS == 1 && WORKER_THREADS_ACTIVE == 1)
				printf ("[Work thread");
			else
				printf ("[Worker #%d", thread_num+1);
			if (buf[0] == '[')
				printf (" %s", buf+1);
			else
				printf ("] %s", buf);
		} else
			printf ("%s", buf);
		last_char_out_was_newline = (buf[strlen(buf)-1] == '\n');
	}
}

void BlinkIcon (int thread_num, int x)
{
}

void ChangeIcon (int thread_num, int x)
{
}


/* This routine calls primeContinue unless there is another copy of mprime */
/* already running.  In that case, it outputs an optional error message. */

void linuxContinue (
	char	*error_message,
	int	thread_num,		/* Specific worker to launch or special value ALL_WORKERS */
	int	wait_flag)
{
#ifdef __linux__
#define PROCNAME	"/proc/%d/exe"
#endif
#ifdef __FreeBSD__
#define PROCNAME	"/proc/%d/file"
#endif
	pid_t	my_pid, running_pid;

/* Compare this process' ID and the pid from the INI file */

	my_pid = getpid ();
	IniFileReread (LOCALINI_FILE);
	running_pid = IniGetInt (LOCALINI_FILE, "Pid", 0);
	if (running_pid == 0 || my_pid == running_pid) goto ok;

#if defined (__APPLE__) || defined (__HAIKU__)
	goto ok;
#elif defined (__OS2__)

        {
            USHORT handle1 = 0, handle2 = 0;
            unsigned char buf[0x2000];
            if( !DosQuerySysState(0x01, 0, 0, 0, (PCHAR)buf, 0x2000) ) {
                PQPROCESS p = ((PQTOPLEVEL)buf)->procdata;
                while(p && p->rectype == 1) {
                    if( p->pid == running_pid ) handle1 = p->hndmod;
                    if( p->pid == my_pid ) handle2 = p->hndmod;
                    p = (PQPROCESS)(p->threads + p->threadcnt);
                }
                if( handle1 != handle2 ) goto ok;
            }
        }

#else

/* See if the two pids are running the same executable */

	{
	char	filename[30];
	int	fd;
	struct stat filedata;
	dev_t	dev1, dev2;
	ino_t	inode1, inode2;

	sprintf (filename, PROCNAME, my_pid);
	fd = _open (filename, _O_RDONLY);
	if (fd < 0) goto ok;
	fstat (fd, &filedata);
	dev1 = filedata.st_dev;
	inode1 = filedata.st_ino;
	_close (fd);
	sprintf (filename, PROCNAME, running_pid);
	fd = _open (filename, _O_RDONLY);
	if (fd < 0) goto ok;
	fstat (fd, &filedata);
	dev2 = filedata.st_dev;
	inode2 = filedata.st_ino;
	_close (fd);
	if (dev1 != dev2 || inode1 != inode2) goto ok;
	}
#endif

/* The two pids are running the same executable, raise an error and return */

	if (error_message != NULL) printf ("%s", error_message);
	return;

/* All is OK, save our pid, run, then delete our pid */

ok:	IniWriteInt (LOCALINI_FILE, "Pid", my_pid);
	LaunchWorkerThreads (thread_num, wait_flag);
	if (wait_flag) IniWriteInt (LOCALINI_FILE, "Pid", 0);
}

/* Implement the rest of the OS-specific routines */

#include "os_routines.c"

