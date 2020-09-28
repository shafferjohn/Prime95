/* Copyright 1995-2020 Mersenne Research, Inc. */
/* Author:  George Woltman */
/* Email: woltman@alum.mit.edu */

/* Include files */

#ifdef __IBMC__
#include <io.h>
#endif
#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "prime.h"

/* Forward declarations */

void sigterm_handler(int);

/* Global variables */

jmp_buf	MENU_JMPBUF;				// Jmpbuf to abort program if a signal is encountered while in the menus or dialogs.

/* Get line from the user (stdin) */

void get_line (
	char	*buf)
{
	int	len;
	buf[0] = 0;
	if (fgets (buf, 80, stdin) == NULL) sigterm_handler (SIGTERM);	// Treat EOF the same as a termination signal
	if (KILL_MENUS) longjmp (MENU_JMPBUF, 1); // If a signal occurred while waiting for user input, longjmp to exit menus
	len = strlen (buf);
	if (len > 0 && buf[len-1] == '\n') buf[len-1] = 0;
}

/* Get a number from the user */

unsigned long get_number (
	unsigned long dflt)
{
	char	line[80];
	unsigned long i;
	get_line (line);
	if (line[0] == 0) return (dflt);
	return (atol (line));
}

unsigned long get_number_or_quit (
	unsigned long dflt,
	unsigned long quit_value)
{
	char	line[80];
	unsigned long i;
	get_line (line);
	if (line[0] == 0) return (dflt);
	if (_stricmp (line, "q") == 0 || _stricmp (line, "quit") == 0) return (quit_value);
	return (atol (line));
}

/* Ask a Yes/No question */

void askYN (
	const char *str,
	int	*val)
{
	char	buf[80];
	printf ("%s (%s): ", str, *val ? "Y" : "N");
	get_line (buf);
	if (buf[0] == 0) return;
	*val = (buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask a number question */

void askNum (
	const char *str,
	unsigned long *val,
	unsigned long min,
	unsigned long max)
{
	char	buf[80];
	unsigned long newval;
	if (min && min >= max) {
		printf ("%s: %ld\n", str, min);
		*val = min;
		return;
	}
	printf ("%s (%ld): ", str, *val);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld: ", min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a number question */

void askInt (
	const char *str,
	long	*val,
	long	min,
	long	max)
{
	char	buf[80];
	long	newval;
	printf ("%s (%ld): ", str, *val);
	if (min && min >= max) {
		printf ("%s: %ld\n", str, min);
		*val = min;
		return;
	}
loop:	get_line (buf);
	if (buf[0] == 0) return;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld. ", min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a number question */

void askFloat (
	const char *str,
	float	*val,
	float	min,
	float	max)
{
	char	buf[80];
	double	newval;
	printf ("%s (%f): ", str, *val);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	newval = (float) atof (buf);
	if (min != 0.0 || max != 0.0) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %f and %f. ", min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a number question */

void askDbl (
	const char *str,
	double	*val,
	double	min,
	double	max)
{
	char	buf[80];
	double	newval;
	printf ("%s (%g): ", str, *val);
loop:	get_line (buf);
	if (buf[0] == 0) return;
	newval = atof (buf);
	if (min != 0.0 || max != 0.0) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %g and %g. ", min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a number question */

void askNumNoDflt (
	const char *str,
	unsigned long *val,
	unsigned long min,
	unsigned long max)
{
	char	buf[80];
	unsigned long newval;
	printf ("%s: ", str);
	if (min && min >= max) {
		printf ("%s: %ld\n", str, min);
		*val = min;
		return;
	}
loop:	get_line (buf);
	if (buf[0] == 0) goto loop;
	newval = atol (buf);
	if (min || max) {
		if (newval < min || newval > max) {
			printf ("Please enter a value between %ld and %ld. ", min, max);
			goto loop;
		}
	}
	*val = newval;
}

/* Ask a string question */

void askStr (
	const char *str,
	char	*val,
	unsigned long maxlen)
{
	char	buf[256], default_val[256];

	strcpy (default_val, val);
	for ( ; ; ) {
		if (default_val[0])
			printf ("%s (%s): ", str, default_val);
		else
			printf ("%s: ", str);
		get_line (buf);
		if (strlen (buf) <= maxlen) break;
		printf ("Maximum string length is %ld characters. ", maxlen);
		buf[maxlen] = 0;
		strcpy (default_val, buf);
	}
	if (buf[0] == 0) strcpy (val, default_val);
	else strcpy (val, buf);
}

/* Wait for user input - gives the user time to read the screen */

void askOK (void)
{
	char	str[80];
	printf ("\nHit enter to continue: ");
	get_line (str);
}

/* Ask user if he is satisfied with his dialog responses */

int askOkCancel (void)
{
	char	buf[80];
	printf ("\nAccept the answers above? (Y): ");
	get_line (buf);
	return (buf[0] == 0 || buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask user if he is satisfied with his dialog responses */

int askYesNo (
	char	dflt)
{
	char	buf[80];
	printf (" (%c): ", dflt);
	get_line (buf);
	if (buf[0] == 0) buf[0] = dflt;
	return (buf[0] == 'Y' || buf[0] == 'y');
}

/* Ask user if he is satisfied with his dialog responses */

int askYesNoCancel (
	char	dflt)
{
	char	buf[80];
	printf (" Y=Yes, N=No, C=Cancel (%c): ", dflt);
	get_line (buf);
	if (buf[0] == 0) buf[0] = dflt;
	return ((buf[0] == 'Y' || buf[0] == 'y') ? 0 :
		(buf[0] == 'N' || buf[0] == 'n') ? 1 : 2);
}

/* Output a long string with a max of 75 characters to a line */

void outputLongLine (
	const char *buf)
{
	char	line[80];
	const char *p;
	int	i, chars_to_output;

	for (p = buf; ; ) {
		chars_to_output = 75;			// Default split point if no natural break found
		for (i = 0; i < 75; i++) {
			line[i] = p[i];
			if (p[i] == 0) { chars_to_output = i; break; }
			if (p[i] == '\n') { chars_to_output = i + 1; break; }
			if (p[i] == ' ') { chars_to_output = i; }
			//if (p[i] == '.') { chars_to_output = i + 1; }		// Don't split at "mersenne.org"
			if (p[i] == ',') { chars_to_output = i + 1; }		// Do split "a,b,c"
		}
		line[chars_to_output] = 0;
		printf ("%s", line);
		p += chars_to_output;
		if (*p == 0) break;
		if (p[-1] != '\n') {			// Do not change any formatting after an explicit line break
			printf ("\n");			// Generate a line break
			while (*p == ' ') p++;		// Eliminate whitespace after a generated line break
		}
	}
}

/* Test/Primenet dialog */

void getProxyInfo (char *, unsigned short *, char *, char *);

void test_primenet (void)
{
	int	m_primenet, m_dialup, m_use_proxy;
	unsigned long m_proxy_port, m_debug;
	char	m_userid[21], m_compid[21], m_proxy_host[120];
	char	m_proxy_user[50], m_proxy_pwd[50], orig_proxy_pwd[50];
	unsigned short proxy_port;
	int	update_computer_info, primenet_debug;

	update_computer_info = FALSE;
	primenet_debug = IniSectionGetInt (INI_FILE, "PrimeNet", "Debug", 0);

	m_primenet = USE_PRIMENET;
	if (USERID[0] == 0)
		strcpy (m_userid, "ANONYMOUS");
	else
		strcpy (m_userid, USERID);
	strcpy (m_compid, COMPID);
	m_dialup = DIAL_UP;
	getProxyInfo (m_proxy_host, &proxy_port, m_proxy_user, m_proxy_pwd);
	m_proxy_port = proxy_port;
	strcpy (orig_proxy_pwd, m_proxy_pwd);
	m_debug = primenet_debug;

	askYN ("Use PrimeNet to get work and report results", &m_primenet);
	if (!m_primenet) goto done;

	outputLongLine ("\nCreate a user account at https://mersenne.org/update/ or you may join GIMPS anonymously but it is not recommended.  See the readme.txt file for details.\n");
	askStr ("Your user ID or \"ANONYMOUS\"", m_userid, 20);
	askStr ("Optional computer name", m_compid, 20);

	askYN ("Computer uses a dial-up connection to the Internet", &m_dialup);

	m_use_proxy = (m_proxy_host[0] != 0);
	askYN ("Use a proxy server", &m_use_proxy);
	if (m_use_proxy) {
		/* Sorry but we can only ask for up to 79 chars for this */
		/* field on the console. User that needs a longer host name */
		/* have to write it in INI directly. */
		askStr ("Proxy host name", m_proxy_host, 79);
		askNum ("Proxy port number", &m_proxy_port, 1, 65535);
		askStr ("Optional proxy user name", m_proxy_user, 49);
		askStr ("Optional proxy password", m_proxy_pwd, 49);
	} else {
		m_proxy_host[0] = 0;
	}
 
	askNum ("Output debug info to prime.log (0=none, 1=some, 2=too much)", &m_debug, 0, 2);

done:	if (askOkCancel ()) {
		DIAL_UP = m_dialup;
		IniWriteInt (INI_FILE, "DialUp", DIAL_UP);

		if (m_proxy_host[0] && m_proxy_port != 8080)
			sprintf (m_proxy_host + strlen (m_proxy_host), ":%lu", m_proxy_port);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyHost", m_proxy_host);
		if (m_proxy_host[0]) {
			IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyUser", m_proxy_user);
			if (strcmp (m_proxy_pwd, orig_proxy_pwd)) {
				IniSectionWriteString (INI_FILE, "PrimeNet",
					"ProxyPass", m_proxy_pwd);
				IniSectionWriteInt (INI_FILE, "PrimeNet",
					"ProxyMask", 0);
			}
		}
		if (m_debug != primenet_debug) {
			IniSectionWriteInt (INI_FILE, "PrimeNet", "Debug",
					    m_debug);
		}

		if (strcmp (USERID, m_userid) != 0) {
			strcpy (USERID, m_userid);
			sanitizeString (USERID);
			IniWriteString (INI_FILE, "V5UserID", USERID);
			update_computer_info = TRUE;
		}
		if (strcmp (COMPID, m_compid) != 0) {
			strcpy (COMPID, m_compid);
			sanitizeString (COMPID);
			IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
			update_computer_info = TRUE;
		}

		if (!USE_PRIMENET && m_primenet) {
			USE_PRIMENET = 1;
			create_window (COMM_THREAD_NUM);
			base_title (COMM_THREAD_NUM, "Communication thread");
			if (!STARTUP_IN_PROGRESS) set_comm_timers ();
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
			spoolExistingResultsFile ();
		} else if (USE_PRIMENET && !m_primenet) {
			USE_PRIMENET = 0;
			if (!STARTUP_IN_PROGRESS) set_comm_timers ();
		} else if (update_computer_info)
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);

		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET);
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	} else
		STARTUP_IN_PROGRESS = 0;
}

/* Test/Worker threads dialog */

int AreAllTheSame (
	unsigned long *array,
	int	arraylen)
{
	int	i;

	for (i = 1; i < arraylen; i++)
		if (array[i-1] != array[i]) return (FALSE);
	return (TRUE);
}

// In theory, the maximum number of workers should be number of logical cpus.
// However, local.txt could specify more worker threads than logical cpus (for
// example, when local.txt is copied from a dual-core to a single-core machine).
// We must let the user manipulate the options on these worker threads that
// don't have a CPU to run on.

unsigned int max_num_workers (void)
{
	if (NUM_WORKER_THREADS >= NUM_CPUS)
		return (NUM_WORKER_THREADS);
	else
		return (NUM_CPUS);
}

void test_worker_threads (void)
{
	unsigned long m_num_thread;
	unsigned long m_work_pref[MAX_NUM_WORKER_THREADS];
	unsigned long m_numcpus[MAX_NUM_WORKER_THREADS];
	int	i, cores_assigned, m_cert_work;

	m_num_thread = NUM_WORKER_THREADS;
	for (i = 0; i < NUM_WORKER_THREADS; i++) {
		m_work_pref[i] = WORK_PREFERENCE[i];
		m_numcpus[i] = CORES_PER_TEST[i];
	}
	for ( ; i < MAX_NUM_WORKER_THREADS; i++) {
		m_work_pref[i] = WORK_PREFERENCE[i];
		m_numcpus[i] = 0;
	}
	m_cert_work = IniGetInt (LOCALINI_FILE, "CertWork", 1);

	if (max_num_workers () <= 1 && !USE_PRIMENET) {
		outputLongLine ("This menu choice only makes sense if you've elected to use PrimeNet to get work and report results.\n");
		askOK ();
		return;
	}

again:	if (max_num_workers () > 1)
		askNum ("Number of workers to run", &m_num_thread, 1, max_num_workers ());

	if (USE_PRIMENET) {
		outputLongLine ("\nUse the following values to select a work type:\n  0 - Whatever makes the most sense\n 150 - First time prime tests\n  152 - World record sized numbers to prime test\n  151 - Double-check prime tests\n  2 - Trial factoring\n  4 - P-1 factoring\n  153 - 100 million digit numbers to prime test\n  160 - First time PRP on Mersenne cofactors\n  161 - Double-check PRP on Mersenne cofactors\n  5 - ECM for first factors of Mersenne numbers\n  8 - ECM on Mersenne cofactors\n  6 - ECM on Fermat numbers\n  1 - Trial factoring to low limits\n");
	}

	if (USE_PRIMENET || NUM_CPUS > 1) {
	    cores_assigned = 0;
	    for (i = 0; i < m_num_thread; i++) {
		if (m_num_thread > 1)
			printf ("\nOptions for worker #%d\n\n", i+1);
		else
			printf ("\n");

		if (USE_PRIMENET) {
			askNum ("Type of work to get", &m_work_pref[i], 0, 161);
			if (m_numcpus[i] < min_cores_for_work_pref (m_work_pref[i]))
				m_numcpus[i] = min_cores_for_work_pref (m_work_pref[i]);
		}

		if (NUM_CPUS > 1) {
			int	min_cores, max_cores;
			min_cores = min_cores_for_work_pref (m_work_pref[i]);
			// Max cores = num_cores_unassigned - num_workers_unconfigured + 1
			max_cores = ((int) NUM_CPUS - cores_assigned) - ((int) m_num_thread - i) + 1;
			if (max_cores < min_cores) max_cores = min_cores;
			if (m_numcpus[i] < min_cores) m_numcpus[i] = min_cores;
			if (m_numcpus[i] > max_cores) m_numcpus[i] = max_cores;
			askNum ("CPU cores to use (multithreading)", &m_numcpus[i], min_cores, max_cores);
		} else
			m_numcpus[i] = 1;

		cores_assigned += m_numcpus[i];
	    }
	}

	if (USE_PRIMENET) {
		askYN ("Get occasional proof certification work", &m_cert_work);
	}

/* Ask user if they are happy with their answers */

	if (askOkCancel ()) {
		int	restart = FALSE;
		int	new_options = FALSE;
		unsigned long i, total_num_cores;

/* If the user has selected 100M tests and per-worker temp disk is not enough for a power=8 proof, then do not permit it. */

		if (CPU_WORKER_DISK_SPACE < 12.0) {
			int	changed = FALSE;
			for (i = 0; i < m_num_thread; i++) {
				if (m_work_pref[i] == PRIMENET_WP_PRP_100M) {
					m_work_pref[i] = PRIMENET_WP_PRP_FIRST;
					changed = TRUE;
				}
			}
			if (changed) {
				outputLongLine (MSG_100M);
				askOK ();
			}
		}

/* If the user has selected first-time tests and per-worker temp disk is not enough for a power=6 proof, then warn the user. */

		if (CPU_WORKER_DISK_SPACE < 1.5) {
			int	warn = FALSE;
			for (i = 0; i < m_num_thread; i++) {
				if (m_work_pref[i] == PRIMENET_WP_PRP_FIRST || m_work_pref[i] == PRIMENET_WP_PRP_WORLD_RECORD) {
					warn = TRUE;
				}
			}
			if (warn) {
				outputLongLine (MSG_FIRST);
				askOK ();
			}
		}

/* If the user has allocated too many cores then raise a severe warning. */

		total_num_cores = 0;
		for (i = 0; i < m_num_thread; i++) total_num_cores += m_numcpus[i];
		if (total_num_cores > NUM_CPUS) {
			outputLongLine (MSG_THREADS);
			if (askYesNo ('Y')) goto again;
		}

/* If user changed the number of worker threads, then make the */
/* necessary changes.  Restart worker threads so that we are running */
/* the correct number of worker threads. */

		if (m_num_thread != NUM_WORKER_THREADS) {
			NUM_WORKER_THREADS = m_num_thread;
			IniWriteInt (LOCALINI_FILE, "WorkerThreads", NUM_WORKER_THREADS);
			new_options = TRUE;
			restart = TRUE;
		}

/* If the user changed any of the work preferences record it in the INI file */
/* and tell the server */

		if (AreAllTheSame (m_work_pref, m_num_thread)) {
			if (! PTOIsGlobalOption (WORK_PREFERENCE) || WORK_PREFERENCE[0] != m_work_pref[0]) {
				PTOSetAll (INI_FILE, "WorkPreference", NULL, WORK_PREFERENCE, m_work_pref[0]);
				new_options = TRUE;
			}
		} else {
			for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
				if (WORK_PREFERENCE[i] == m_work_pref[i]) continue;
				PTOSetOne (INI_FILE, "WorkPreference", NULL, WORK_PREFERENCE, i, m_work_pref[i]);
				new_options = TRUE;
			}
		}

/* If the user changed any of the cores_per_test record it in the INI file */

		if (AreAllTheSame (m_numcpus, m_num_thread))
			PTOSetAll (LOCALINI_FILE, "CoresPerTest", NULL, CORES_PER_TEST, m_numcpus[0]);
		else for (i = 0; i < (int) NUM_WORKER_THREADS; i++)
			PTOSetOne (LOCALINI_FILE, "CoresPerTest", NULL, CORES_PER_TEST, i, m_numcpus[i]);

/* Write the new CertWork setting */

		IniWriteInt (LOCALINI_FILE, "CertWork", m_cert_work);

/* Send new settings to the server */

		if (new_options) spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);

/* Restart worker threads with new options */

		if (restart) stop_workers_for_restart ();
	} else
		STARTUP_IN_PROGRESS = 0;
}

/* Output a status report for the range */

void test_status (void)
{
	char	buf[4000];

	rangeStatusMessage (buf, sizeof (buf));
	strcat (buf, "\n");
	outputLongLine (buf);
}

/* Start one or all workers */

void test_continue (void)
{
	unsigned long worker;
	int	thread_num;

	worker = 0;
	askNum ("Worker to start, 0=all", &worker, 0, WORKER_THREADS_ACTIVE > NUM_WORKER_THREADS ? WORKER_THREADS_ACTIVE : NUM_WORKER_THREADS);
	if (worker == 0) thread_num = ALL_WORKERS;
	else thread_num = worker - 1;
	linuxContinue ("Another mprime is running.\n", thread_num, FALSE);
}

/* Stop one or all workers */

void test_stop (void)
{
	unsigned long worker;

	worker = 0;
	askNum ("Worker to stop, 0=all", &worker, 0, WORKER_THREADS_ACTIVE);
	if (worker == 0) stop_workers_for_escape ();
	else stop_one_worker (worker - 1);
}

/* Start or stop one or all workers */

void test_continue_or_stop (void)
{
	outputLongLine ("Do you want to start some workers? ");
	if (askYesNo ('Y')) test_continue ();
	else test_stop ();
}

/* Advanced/Test dialog */

void advanced_test (void)
{
	unsigned long m_thread, m_p;
#define NOTPRIMEERR "This number is not prime, there is no need to test it.\n"

loop:	m_p = 0;

	m_thread = 1;
	if (NUM_WORKER_THREADS > 1)
		askNum ("Worker number", &m_thread, 1, NUM_WORKER_THREADS);

	askNumNoDflt ("Exponent to test", &m_p, MIN_PRIME,
		      CPU_FLAGS & CPU_FMA3 ? MAX_PRIME_FMA3 :
		      CPU_FLAGS & CPU_SSE2 ? MAX_PRIME_SSE2 : MAX_PRIME);

	if (askOkCancel ()) {
		struct work_unit w;
		if (! isPrime (m_p)) {
			printf (NOTPRIMEERR);
			goto loop;
		}
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_ADVANCEDTEST;
		w.k = 1.0;
		w.b = 2;
		w.n = m_p;
		w.c = -1;
		addWorkToDoLine (m_thread - 1, &w);
		if (WORKER_THREADS_ACTIVE)
			stop_worker_for_advanced_test (m_thread - 1);
		else
			linuxContinue ("\nWork added to worktodo.txt file.  Another mprime is running.\n", ALL_WORKERS, FALSE);
	}
}

/* Advanced/Time dialog */

void advanced_time (void)
{
	unsigned long m_p, m_iter;

	m_p = 10000000;
	m_iter = 10;

	askNum ("Exponent to time", &m_p, MIN_PRIME,
	        CPU_FLAGS & CPU_FMA3 ? MAX_PRIME_FMA3 :
	        CPU_FLAGS & CPU_SSE2 ? MAX_PRIME_SSE2 : MAX_PRIME);
	askNum ("Number of Iterations", &m_iter, 1, 1000);
	if (askOkCancel ()) {
		LaunchAdvancedTime (m_p, m_iter);
	}
}

/* Advanced/P-1 dialog */

void advanced_pminus1 (void)
{
	unsigned long m_thread, m_b, m_n, m_num_curves;
	long	m_c;
	double	m_k, m_bound1, m_bound2;

	m_k = 1.0;
	m_b = 2;
	m_n = 0;
	m_c = -1;
	m_bound1 = 50000.0;
	m_bound2 = 0.0;

	m_thread = 1;
	if (NUM_WORKER_THREADS > 1)
		askNum ("Worker number", &m_thread, 1, NUM_WORKER_THREADS);

	askDbl ("k in k*b^n+c", &m_k, 1.0, 1.0e15);
	askNum ("b in k*b^n+c", &m_b, 2, 1000000000);
	askNumNoDflt ("n in k*b^n+c", &m_n, 100, 600000000);
	askInt ("c in k*b^n+c", &m_c, -2000000000, 2000000000);
	askDbl ("Bound #1", &m_bound1, 100.0, 1.0e15);
	askDbl ("Bound #2", &m_bound2, 0.0, 1.0e15);

	if (askOkCancel ()) {
		struct work_unit w;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_PMINUS1;
		w.k = m_k;
		w.b = m_b;
		w.n = m_n;
		w.c = m_c;
		w.B1 = m_bound1;
		w.B2_start = 0;
		w.B2 = m_bound2;
		addWorkToDoLine (m_thread - 1, &w);
		if (!WORKER_THREADS_ACTIVE)
			linuxContinue ("\nWork added to worktodo.txt file.  Another mprime is running.\n", ALL_WORKERS, FALSE);
		askOK ();
	}
}

/* Advanced/ECM dialog */

void advanced_ecm (void)
{
	unsigned long m_thread, m_b, m_n, m_num_curves;
	long	m_c;
	double	m_k, m_bound1, m_bound2;

	m_k = 1.0;
	m_b = 2;
	m_n = 0;
	m_c = -1;
	m_bound1 = 50000.0;
	m_bound2 = 0.0;
	m_num_curves = 100;

	m_thread = 1;
	if (NUM_WORKER_THREADS > 1)
		askNum ("Worker number", &m_thread, 1, NUM_WORKER_THREADS);

	askDbl ("k in k*b^n+c", &m_k, 1.0, 1.0e15);
	askNum ("b in k*b^n+c", &m_b, 2, 1000000000);
	askNumNoDflt ("n in k*b^n+c", &m_n, 100, 600000000);
	askInt ("c in k*b^n+c", &m_c, -2000000000, 2000000000);
	askDbl ("Bound #1", &m_bound1, 100.0, 1.0e15);
	askDbl ("Bound #2", &m_bound2, 0.0, 1.0e15);
	askNum ("Curves to test", &m_num_curves, 1, 100000);

	if (askOkCancel ()) {
		struct work_unit w;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_ECM;
		w.k = m_k;
		w.b = m_b;
		w.n = m_n;
		w.c = m_c;
		w.B1 = m_bound1;
		w.B2_start = 0;
		w.B2 = m_bound2;
		w.curves_to_do = m_num_curves;
		w.curve = 0.0;
		addWorkToDoLine (m_thread - 1, &w);
		if (!WORKER_THREADS_ACTIVE)
			linuxContinue ("\nWork added to worktodo.txt file.  Another mprime is running.\n", ALL_WORKERS, FALSE);
		askOK ();
	}
}

/* Advanced/Manual Communication dialog */

void advanced_manualcomm (void)
{
	int	m_manual_comm, m_comm_now, m_new_dates;

	m_manual_comm = MANUAL_COMM;
	m_comm_now = 1;
	m_new_dates = 0;

	m_manual_comm = !m_manual_comm;
	askYN ("Contact PrimeNet server automatically", &m_manual_comm);
	m_manual_comm = !m_manual_comm;
	askYN ("Contact PrimeNet server now", &m_comm_now);
	askYN ("Send new expected completion dates to server", &m_new_dates);

	if (askOkCancel ()) {
		if ((MANUAL_COMM && !m_manual_comm) ||
		    (!MANUAL_COMM && m_manual_comm)) {
			MANUAL_COMM = m_manual_comm;
			IniWriteInt (INI_FILE, "ManualComm", MANUAL_COMM);
			set_comm_timers ();
		}
		if (m_new_dates) UpdateEndDates ();
		if (m_comm_now) do_manual_comm_now ();
	}
}

/* Advanced/Time dialog */

void advanced_unreserve (void)
{
	unsigned long m_p;

	m_p = 0;

	outputLongLine ("\nUse this only if you are sure you will not be finishing this exponent.  The exponent will be assigned to someone else.  It is not fair to them if you test an exponent assigned to someone else.\n");
	askNumNoDflt ("Exponent to unreserve", &m_p, 1000, 1000000000);
	if (askOkCancel ()) unreserve (m_p);
}

/* Advanced/Quit Gimps dialog */

void advanced_quit (void)
{

	if (!USE_PRIMENET) {
		outputLongLine (MANUAL_QUIT);
		if (askYesNo ('N')) {
			writeResults ("Quitting GIMPS.\n");
//bug - either delete file, or delete all work_units and write the file.
//bug			IniDeleteAllLines (WORKTODO_FILE);
			stop_workers_for_escape ();
		}
	} else {
		int	res;
		outputLongLine (PRIMENET_QUIT);
		res = askYesNoCancel ('C');
		if (res == 0) {
			OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS after current work completes.\n");
			IniWriteInt (INI_FILE, "NoMoreWork", 1);
			askOK ();
		}
		if (res == 1) {
			OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS immediately.\n");
			spoolMessage (MSG_QUIT_GIMPS, NULL);
			askOK ();
		}
	}
}

/* Options/CPU dialog */

void options_cpu (void)
{
	unsigned long m_hours;
	char buf[512];

	m_hours = CPU_HOURS;

	askNum ("Hours per day this program will run", &m_hours, 1, 24);

	getCpuDescription (buf, 0);
	printf ("\nCPU Information:\n%s\n", buf);

	if (askOkCancel ()) {
		if (CPU_HOURS != m_hours) {
			CPU_HOURS = m_hours;
			IniWriteInt (LOCALINI_FILE, "CPUHours", CPU_HOURS);
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
			delete_timed_event (TE_COMM_SERVER);
			UpdateEndDates ();
		}
	} else
		STARTUP_IN_PROGRESS = 0;
}

/* Options/Preferences dialog */

#define round_to_tenth(a)	((round((a) * 10.0)) / 10.0)

void options_resources (void)
{
	unsigned long m_download_mb, m_priority, m_cert_cpu;
	float	m_disk, m_upload_bandwidth, m_day_memory, m_night_memory, m_emergency_mem;
	char	m_upload_start[20], m_upload_end[20], m_temp_dir[512], m_archive_dir[512];
	char	m_start_time[13], m_end_time[13];
	int	can_upload, m_hyper_tf, m_hyper_ll, m_memory_editable;
	unsigned int day_memory, night_memory, day_start_time, day_end_time;

	outputLongLine ("Consult readme.txt prior to changing any of these settings.\n\n");

	m_disk = CPU_WORKER_DISK_SPACE;
	m_upload_bandwidth = IniSectionGetFloat (INI_FILE, "PrimeNet", "UploadRateLimit", 0.25);
	if (m_upload_bandwidth <= 0.0 || m_upload_bandwidth > 10000.0) m_upload_bandwidth = 10000.0;
	IniSectionGetString (INI_FILE, "PrimeNet", "UploadStartTime", m_upload_start, sizeof (m_upload_start), "00:00");
	if (strcmp (m_upload_start, "00:00") != 0) minutesToStr (strToMinutes (m_upload_start), m_upload_start);
	IniSectionGetString (INI_FILE, "PrimeNet", "UploadEndTime", m_upload_end, sizeof (m_upload_end), "24:00");
	if (strcmp (m_upload_end, "24:00") != 0) minutesToStr (strToMinutes (m_upload_end), m_upload_end);
	m_download_mb = IniSectionGetInt (INI_FILE, "PrimeNet", "DownloadDailyLimit", 40);
	IniGetString (LOCALINI_FILE, "ProofResiduesDir", m_temp_dir, sizeof (m_temp_dir), NULL);
	IniGetString (LOCALINI_FILE, "ProofArchiveDir", m_archive_dir, sizeof (m_archive_dir), NULL);
	m_memory_editable = read_memory_settings (&day_memory, &night_memory, &day_start_time, &day_end_time);
	m_day_memory = (float) round_to_tenth (day_memory / 1024.0);
	m_night_memory = (float) round_to_tenth (night_memory / 1024.0);
	minutesToStr (day_start_time, m_start_time);
	minutesToStr (day_end_time, m_end_time);
	m_emergency_mem = (float) round_to_tenth (IniGetInt (LOCALINI_FILE, "MaxEmergencyMemory", 1024) / 1024.0);
	m_priority = PRIORITY;
	m_cert_cpu = IniGetInt (LOCALINI_FILE, "CertDailyCPULimit", 10);
	m_hyper_tf = HYPERTHREAD_TF;
	m_hyper_ll = HYPERTHREAD_LL;
	can_upload = IniSectionGetInt (INI_FILE, "PrimeNet", "ProofUploads", 1);

	askFloat ("Temporary disk space limit in GB/worker", &m_disk, 0.0, 1000.0);
	if (can_upload) {
		askFloat ("Upload bandwidth limit in Mbps", &m_upload_bandwidth, 0.05, 10000.0);
		askStr ("Upload large files time period start", m_upload_start, 8);
		askStr ("Upload large files time period end", m_upload_end, 8);
	}
	askNum ("Download limit for certification work in MB/day", &m_download_mb, 0, 999999);

	outputLongLine ("Skip advanced resource settings");
	if (!askYesNo ('Y')) {
		float	max_emergency_mem = (float) (0.25 * physical_memory () / 1024.0);
		if (max_emergency_mem < 1.0) max_emergency_mem = 1.0;
		askStr ("Optional directory to hold large temporary files", m_temp_dir, 511);
		askStr ("Optional directory to hold archived proofs", m_archive_dir, 511);
		if (m_memory_editable) {
			float	max_mem;
			max_mem = (float) (0.9 * physical_memory () / 1024.0);
			askFloat ("Daytime P-1/ECM stage 2 memory in GB", &m_day_memory, 0.0, max_mem);
			askFloat ("Nighttime P-1/ECM stage 2 memory in GB", &m_night_memory, 0.0, max_mem);
			if (m_day_memory != m_night_memory) {
				askStr ("Daytime begins at", m_start_time, 12);
				askStr ("Daytime ends at", m_end_time, 12);
			}
		}
		askFloat ("Max emergency memory in GB/worker", &m_emergency_mem, 0.0, max_emergency_mem);
		askNum ("Priority -- 1 is highly recommended, see readme.txt", &m_priority, 1, 10);
		if (m_download_mb) {
			askNum ("Certification work limit in % of CPU time", &m_cert_cpu, 1, 100);
		}
		if (CPU_HYPERTHREADS > 1 && OS_CAN_SET_AFFINITY) {
			askYN ("Use hyperthreading for trial factoring (recommended)", &m_hyper_tf);
			askYN ("Use hyperthreading for PRP, LL, P-1, ECM (not recommended)", &m_hyper_ll);
		}
	}

	if (askOkCancel ()) {
		int	restart = FALSE;
		unsigned int new_day_start_time, new_day_end_time;

		// Raise a warning if uesr drops the temp disk space below the threshold for first time work.
		if (CPU_WORKER_DISK_SPACE >= 1.5 && m_disk < 1.5) {
			outputLongLine (MSG_DISK);
			askOK ();
		}
		CPU_WORKER_DISK_SPACE = m_disk;
		IniWriteFloat (LOCALINI_FILE, "WorkerDiskSpace", CPU_WORKER_DISK_SPACE);
		IniSectionWriteFloat (INI_FILE, "PrimeNet", "UploadRateLimit", m_upload_bandwidth);
		IniSectionWriteString (INI_FILE, "PrimeNet", "UploadStartTime", m_upload_start);
		IniSectionWriteString (INI_FILE, "PrimeNet", "UploadEndTime", m_upload_end);
		IniSectionWriteInt (INI_FILE, "PrimeNet", "DownloadDailyLimit", m_download_mb);
		IniWriteString (LOCALINI_FILE, "ProofResiduesDir", m_temp_dir);
		IniWriteString (LOCALINI_FILE, "ProofArchiveDir", m_archive_dir);

/* Save the new memory settings */

		new_day_start_time = strToMinutes (m_start_time);
		new_day_end_time = strToMinutes (m_end_time);
		if (day_memory != (int) (m_day_memory * 1024.0)  ||
		    night_memory != (int) (m_night_memory * 1024.0) ||
		    day_start_time != new_day_start_time ||
		    day_end_time != new_day_end_time) {
			write_memory_settings ((int) (m_day_memory * 1024.0), (int) (m_night_memory * 1024.0), new_day_start_time, new_day_end_time);
			mem_settings_have_changed ();
			spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
		}
		IniWriteInt (LOCALINI_FILE, "MaxEmergencyMemory", (long) (m_emergency_mem * 1024.0));

/* If user changed the priority of worker threads, then change the INI file. */
/* Restart worker threads so that they are running at the new priority. */

		if (PRIORITY != m_priority) {
			PRIORITY = m_priority;
			IniWriteInt (INI_FILE, "Priority", PRIORITY);
			restart = TRUE;
		}

/* Handle cert work CPU limit */

		IniWriteInt (LOCALINI_FILE, "CertDailyCPULimit", m_cert_cpu);

/* If user changed the hyperthreading options, then save the options to the INI file */

		if (m_hyper_tf != HYPERTHREAD_TF) {
			HYPERTHREAD_TF = m_hyper_tf;
			IniWriteInt (LOCALINI_FILE, "HyperthreadTF", HYPERTHREAD_TF);
			restart = TRUE;
		}
		if (m_hyper_ll != HYPERTHREAD_LL) {
			HYPERTHREAD_LL = m_hyper_ll;
			IniWriteInt (LOCALINI_FILE, "HyperthreadLL", HYPERTHREAD_LL);
			restart = TRUE;
		}

/* Restart worker threads with new options */

		if (restart) stop_workers_for_restart ();
	}
}

/* Options/Preferences dialog */

void options_preferences (void)
{
	unsigned long m_iter, m_r_iter, m_disk_write_time;
	unsigned long m_modem, m_retry, m_work, m_backup;
	float	m_end_dates;
	int	m_noise, m_battery;

	m_iter = ITER_OUTPUT;
	m_r_iter = ITER_OUTPUT_RES;
	m_disk_write_time = DISK_WRITE_TIME;
	m_modem = MODEM_RETRY_TIME;
	m_retry = NETWORK_RETRY_TIME;
	m_work = DAYS_OF_WORK;
	m_end_dates = DAYS_BETWEEN_CHECKINS;
	m_backup = NUM_BACKUP_FILES;
	m_noise = !SILENT_VICTORY;
	m_battery = RUN_ON_BATTERY;

	askNum ("Iterations between screen outputs", &m_iter, 1, 999999999);
	askNum ("Iterations between results file outputs",
		&m_r_iter, 10000, 999999999);
	askNum ("Minutes between disk writes", &m_disk_write_time, 10, 999999);
	if (USE_PRIMENET && DIAL_UP)
		askNum ("Minutes between modem retries", &m_modem, 1, 300);
	if (USE_PRIMENET)
		askNum ("Minutes between network retries", &m_retry, 1, 300);
	if (USE_PRIMENET)
		askNum ("Days of work to queue up", &m_work, 1, 90);
	if (USE_PRIMENET)
		askFloat ("Days between sending end dates", &m_end_dates, 0.125, 7);
	askNum ("Number of Backup Files", &m_backup, 1, 3);
	askYN ("Make noise if new Mersenne prime is found", &m_noise);
	askYN ("Run program even when using laptop battery power", &m_battery);

	if (askOkCancel ()) {
		ITER_OUTPUT = m_iter;
		ITER_OUTPUT_RES = m_r_iter;
		DISK_WRITE_TIME = m_disk_write_time;
		MODEM_RETRY_TIME = m_modem;
		NETWORK_RETRY_TIME = m_retry;
		DAYS_OF_WORK = m_work;
		DAYS_BETWEEN_CHECKINS = m_end_dates;
		NUM_BACKUP_FILES = m_backup;
		SILENT_VICTORY = !m_noise;
		if (RUN_ON_BATTERY != m_battery) {
			RUN_ON_BATTERY = m_battery;
			IniWriteInt (LOCALINI_FILE, "RunOnBattery", RUN_ON_BATTERY);
			run_on_battery_changed ();
		}
		IniWriteInt (INI_FILE, "OutputIterations", ITER_OUTPUT);
		IniWriteInt (INI_FILE, "ResultsFileIterations", ITER_OUTPUT_RES);
		IniWriteInt (INI_FILE, "DiskWriteTime", DISK_WRITE_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime", MODEM_RETRY_TIME);
		IniWriteInt (INI_FILE, "NetworkRetryTime2", NETWORK_RETRY_TIME);
		IniWriteInt (INI_FILE, "DaysOfWork", DAYS_OF_WORK);
		IniWriteFloat (INI_FILE, "DaysBetweenCheckins", DAYS_BETWEEN_CHECKINS);
		IniWriteInt (INI_FILE, "NumBackupFiles", NUM_BACKUP_FILES);
		IniWriteInt (INI_FILE, "SilentVictory", SILENT_VICTORY);
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
}

/* Options/Torture test dialog */

void options_torture (void)
{
	unsigned long m_thread, m_type, m_minfft, m_maxfft;
	unsigned long m_memory, m_timefft;
	unsigned long mem, blendmemory;
	int	m_custom, m_weak, m_avx512, m_fma3, m_avx, m_sse2;

	m_thread = NUM_CPUS * CPU_HYPERTHREADS;
	mem = physical_memory ();
	// New in 29.5, default to all but 3GB of memory for the large FFT and blend test.
	// On my Skylake-X linux machine mprime crashes if it uses all but 2.5GB (maybe due
	// to large pages allocated by a running production mprime).
	if (mem >= 5000) {
		blendmemory = GetSuggestedMemory (mem - 3000);
	} else if (mem >= 2000) {
		blendmemory = GetSuggestedMemory (mem - 512);
	} else if (mem >= 500) {
		blendmemory = GetSuggestedMemory (mem - 256);
	} else if (mem >= 200) {
		blendmemory = GetSuggestedMemory (mem / 2);
	} else {
		blendmemory = 8;
	}

	// Get number of threads to run.  It will affect our FFT size calculations.
	if (NUM_CPUS * CPU_HYPERTHREADS > 1)
		askNum ("Number of torture test threads to run", &m_thread, 1, NUM_CPUS * CPU_HYPERTHREADS);

	// Ask which torture test to run
	if (CPU_TOTAL_L4_CACHE_SIZE) {
		m_type = 5;
		outputLongLine ("Choose a type of torture test to run.\n  1 = Smallest FFTs (tests L1/L2 caches, high power/heat/CPU stress).\n  2 = Small FFTs (tests L1/L2/L3 caches, maximum power/heat/CPU stress).\n  3 = Medium FFTs (tests L1/L2/L3/L4 caches, high power/heat/CPU stress).\n  4 = Large FFTs (stresses memory controller and RAM).\n  5 = Blend (tests all of the above).\nBlend is the default.  NOTE: if you fail the blend test but pass the smaller FFT tests then your problem is likely bad memory or bad memory controller.\n");
		askNum ("Type of torture test to run", &m_type, 1, 5);
	} else if (CPU_TOTAL_L3_CACHE_SIZE) {
		m_type = 4;
		outputLongLine ("Choose a type of torture test to run.\n  1 = Smallest FFTs (tests L1/L2 caches, high power/heat/CPU stress).\n  2 = Small FFTs (tests L1/L2/L3 caches, maximum power/heat/CPU stress).\n  3 = Large FFTs (stresses memory controller and RAM).\n  4 = Blend (tests all of the above).\nBlend is the default.  NOTE: if you fail the blend test but pass the smaller FFT tests then your problem is likely bad memory or bad memory controller.\n");
		askNum ("Type of torture test to run", &m_type, 1, 4);
		if (m_type >= 3) m_type++;
	} else {
		m_type = 3;
		outputLongLine ("Choose a type of torture test to run.\n  1 = Smallest FFTs (tests L1/L2 caches, maximum power/heat/CPU stress).\n  2 = Large FFTs (stresses memory controller and RAM).\n  3 = Blend (tests all of the above).\nBlend is the default.  NOTE: if you fail the blend test but pass the smaller FFT tests then your problem is likely bad memory or bad memory controller.\n");
		askNum ("Type of torture test to run", &m_type, 1, 3);
		if (m_type >= 2) m_type += 2;
	}

	// Calculate default FFT sizes
	{
		int	minfft, maxfft;
		tortureTestDefaultSizes (m_type - 1, m_thread, &minfft, &maxfft);
		m_minfft = minfft, m_maxfft = maxfft;
	}
	if (m_minfft < 4) m_minfft = 4;
	if (m_maxfft < m_minfft) m_maxfft = m_minfft;
	if (m_maxfft > 32768) m_maxfft = 32768;

	// Assign other options
	m_memory = (m_type <= 3 ? 0 : blendmemory);
	m_timefft = (m_thread > NUM_CPUS ? 6 : 3);

	// Let user customize
	m_custom = FALSE;
	askYN ("Customize settings", &m_custom);
	if (m_custom) {
		askNum ("Min FFT size (in K)", &m_minfft, 4,
			CPU_FLAGS & CPU_AVX512F && !m_avx512 ? MAX_FFTLEN_AVX512 / 1024 :
			CPU_FLAGS & CPU_FMA3 && !m_fma3 ? MAX_FFTLEN_FMA3 / 1024 :
			CPU_FLAGS & CPU_SSE2 && !m_sse2 ? MAX_FFTLEN_SSE2 / 1024 :
							  MAX_FFTLEN / 1024);
		askNum ("Max FFT size (in K)", &m_maxfft, m_minfft,
			CPU_FLAGS & CPU_AVX512F && !m_avx512 ? MAX_FFTLEN_AVX512 / 1024 :
			CPU_FLAGS & CPU_FMA3 && !m_fma3 ? MAX_FFTLEN_FMA3 / 1024 :
			CPU_FLAGS & CPU_SSE2 && !m_sse2 ? MAX_FFTLEN_SSE2 / 1024 :
							  MAX_FFTLEN / 1024);
		askNum ("Memory to use (in MB, 0 = in-place FFTs)", &m_memory, 0, mem);
		askNum ("Time to run each FFT size (in minutes)", &m_timefft, 1, 60);
	}

	// Ask about running a less stressful torture test.  Stop on first negative answer.
	// X86-64 cannot go any lower than SSE2.
	m_weak = m_avx512 = m_fma3 = m_avx = m_sse2 = 0;
	askYN ("Run a weaker torture test (not recommended)", &m_weak);
	if (m_weak && CPU_FLAGS & CPU_AVX512F) askYN ("Disable AVX-512", &m_avx512), m_weak = m_avx512;
	if (m_weak && CPU_FLAGS & CPU_FMA3) askYN ("Disable AVX2 (fused multiply add)", &m_fma3), m_weak = m_fma3;
	if (m_weak && CPU_FLAGS & CPU_AVX) askYN ("Disable AVX", &m_avx), m_weak = m_avx;
#ifndef X86_64
	if (m_weak && CPU_FLAGS & CPU_SSE2) askYN ("Disable SSE2", &m_sse2), m_weak = m_sse2;
#endif

	if (askOkCancel ()) {
		IniWriteInt (INI_FILE, "MinTortureFFT", m_minfft);
		IniWriteInt (INI_FILE, "MaxTortureFFT", m_maxfft);
		IniWriteInt (INI_FILE, "TortureMem", m_memory);
		IniWriteInt (INI_FILE, "TortureTime", m_timefft);
		m_weak = m_avx512 * CPU_AVX512F + m_fma3 * CPU_FMA3 + m_avx * CPU_AVX + m_sse2 * CPU_SSE2;
		IniWriteInt (INI_FILE, "TortureWeak", m_weak);
		LaunchTortureTest (m_thread, TRUE);
	}
}

/* Options/Benchmark */

void options_benchmark (void)
{
	unsigned long m_bench_type, m_minFFT, m_maxFFT, m_bench_time;
	char	m_cores[512], m_workers[512];
	int	m_errchk, m_all_complex, m_limit_FFT_sizes, m_hyperthreading, m_all_FFT_impl;

	m_bench_type = 0;
	askNum ("Benchmark type (0 = Throughput, 1 = FFT timings, 2 = Trial factoring)", &m_bench_type, 0, 2);

	if (m_bench_type != 2) {
		printf ("\nFFTs to benchmark\n");
		m_minFFT = IniGetInt (INI_FILE, "MinBenchFFT", 2048);
		askNum ("Minimum FFT size (in K)", &m_minFFT, 1, 65536);
		m_maxFFT = IniGetInt (INI_FILE, "MaxBenchFFT", 8192);
		askNum ("Maximum FFT size (in K)", &m_maxFFT, m_minFFT, 65536);
		m_errchk = ERRCHK;			// IniGetInt (INI_FILE, "BenchErrorCheck", 0);
		askYN ("Benchmark with round-off checking enabled", &m_errchk);
		m_all_complex = 0;			// IniGetInt (INI_FILE, "BenchAllComplex", 0);
		askYN ("Benchmark all-complex FFTs (for LLR,PFGW,PRP users)", &m_all_complex);
		m_limit_FFT_sizes = 0;			// IniGetInt (INI_FILE, "OnlyBench5678", 1);
		if (m_minFFT != m_maxFFT) askYN ("Limit FFT sizes (mimic older benchmarking code)", &m_limit_FFT_sizes);
	}

	sprintf (m_cores, "%lu", NUM_CPUS);
	m_hyperthreading = IniGetInt (INI_FILE, "BenchHyperthreads", 1);
	if (NUM_CPUS > 1 || CPU_HYPERTHREADS > 1) {
		printf ("\nCPU cores to benchmark\n");
		if (NUM_CPUS > 1) askStr ("Number of CPU cores (comma separated list of ranges)", m_cores, 511);
		if (CPU_HYPERTHREADS > 1) askYN ("Benchmark hyperthreading", &m_hyperthreading);
	}

	if (m_bench_type == 0) {
		int	i, max_cores, vals[4], numvals;

		printf ("\nThroughput benchmark options\n");

		m_all_FFT_impl = IniGetInt (INI_FILE, "AllBench", 0);
		askYN ("Benchmark all FFT implementations to find best one for your machine", &m_all_FFT_impl);

		// To come up with a rational default for number of workers, we need to know the maximum number of
		// cores the benchmark will be running on.
		max_cores = 1;
		for (i = 2; i <= NUM_CPUS; i++) if (is_number_in_list (i, m_cores)) max_cores = i;
		// If testing all FFT implementations. then default to the current num_workers.
		numvals = 0;
		if (NUM_WORKER_THREADS <= max_cores) sorted_add_unique (vals, &numvals, NUM_WORKER_THREADS);
		else sorted_add_unique (vals, &numvals, max_cores);
		// Otherwise, assume user is trying to figure out how many workers to run and form a string
		// with the most common best values for number of workers: 1, num_threading_nodes, num_cores, num_workers
		if (!m_all_FFT_impl) {
			sorted_add_unique (vals, &numvals, 1);
			if (NUM_THREADING_NODES <= max_cores) sorted_add_unique (vals, &numvals, NUM_THREADING_NODES);
			sorted_add_unique (vals, &numvals, max_cores);
		}
		sprintf (m_workers, "%d", vals[0]);
		for (i = 1; i < numvals; i++) sprintf (m_workers + strlen (m_workers), ",%d", vals[i]);

		if (max_cores > 1) askStr ("Number of workers (comma separated list of ranges)", m_workers, 511);

		m_bench_time = IniGetInt (INI_FILE, "BenchTime", 15);
		askNum ("Time to run each benchmark (in seconds)", &m_bench_time, 5, 60);
	}

	if (askOkCancel ()) {
		if (m_bench_type != 2) {
			IniWriteInt (INI_FILE, "MinBenchFFT", m_minFFT);
			IniWriteInt (INI_FILE, "MaxBenchFFT", m_maxFFT);
			IniWriteInt (INI_FILE, "BenchErrorCheck", m_errchk);
			IniWriteInt (INI_FILE, "BenchAllComplex", m_all_complex ? 2 : 0);
			IniWriteInt (INI_FILE, "OnlyBench5678", m_limit_FFT_sizes);
		}
		IniWriteString (INI_FILE, "BenchCores", m_cores);
		IniWriteInt (INI_FILE, "BenchHyperthreads", m_hyperthreading);
		if (m_bench_type == 0) {
			IniWriteString (INI_FILE, "BenchWorkers", m_workers);
			IniWriteInt (INI_FILE, "AllBench", m_all_FFT_impl);
			IniWriteInt (INI_FILE, "BenchTime", m_bench_time);
		}
		LaunchBench (m_bench_type);
	}
}

/* Help/About */

void help_about (void)
{
	char	app_string[120];

	generate_application_string (app_string);
	printf ("GIMPS: Mersenne Prime Search\n");
	printf ("Web site: http://mersenne.org\n");
	printf ("%s\n", app_string);
	printf ("Copyright 1996-2020 Mersenne Research, Inc.\n");
	printf ("Author: George Woltman\n");
	printf ("Email:  woltman@alum.mit.edu\n");
	askOK ();
}

/* Help/About PrimeNet Server */

void help_about_server (void)
{
	char	*buildId;
	char	*serverName;
	char	*adminEmailAddr;

	pingServer ();
}

/* Welcome Information dialog */

void test_welcome (void)
{
	int	m_join = 1;

/* Set global flag indicating startup is in progress.  This will delay */
/* starting any communication with the server until the user has confirmed */
/* he wants to use primenet and he has selected his work preferences. */
	
	STARTUP_IN_PROGRESS = 1;

	outputLongLine ("\nWelcome to GIMPS, the hunt for huge prime numbers.  You will be asked a few simple questions and then the program will contact the primenet server to get some work for your computer.  Good luck!\n");
	outputLongLine ("\nAttention OVERCLOCKERS!!  Mprime has gained a reputation as a useful stress testing tool for people that enjoy pushing their hardware to the limit.  You are more than welcome to use this software for that purpose.  Please select the stress testing choice below to avoid interfering with the PrimeNet server.  Use the Options/Torture Test menu choice for your stress tests.  Also, read the stress.txt file.\n");
	outputLongLine ("\nIf you want to both join GIMPS and run stress tests, then Join GIMPS and answer the questions.  After the server gets some work for you, stop mprime, then run mprime -m and choose Options/Torture Test.\n\n");
	askYN ("Join Gimps? (Y=Yes, N=Just stress testing)", &m_join);
	if (m_join) {
		STRESS_TESTER = 0;
		IniWriteInt (INI_FILE, "StressTester", 0);
		USE_PRIMENET = 1;
		IniWriteInt (INI_FILE, "UsePrimenet", 1);
		test_primenet ();
		if (USE_PRIMENET && STARTUP_IN_PROGRESS) options_cpu ();
		if (STARTUP_IN_PROGRESS) options_resources ();
		if (USE_PRIMENET && STARTUP_IN_PROGRESS) test_worker_threads ();
		if (USE_PRIMENET && STARTUP_IN_PROGRESS) {
			STARTUP_IN_PROGRESS = 0;
			set_comm_timers ();
			linuxContinue (NULL, ALL_WORKERS, FALSE);
		} else
			STARTUP_IN_PROGRESS = 0;
	} else {
		STRESS_TESTER = 1;
		IniWriteInt (INI_FILE, "StressTester", 1);
		USE_PRIMENET = 0;
		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET = 0);
		STARTUP_IN_PROGRESS = 0;
		options_torture ();
	}
	main_menu ();
}

/* Display the main menu */

void main_menu (void)
{
	unsigned long choice;

	if (setjmp (MENU_JMPBUF)) return;		// Exit menus if a signal received while in the menus
	
	for ( ; ; ) {

	printf ("\t     Main Menu\n");
	printf ("\n");
	printf ("\t 1.  Test/Primenet\n");
	printf ("\t 2.  Test/Worker threads\n");
	printf ("\t 3.  Test/Status\n");
	if (WORKER_THREADS_ACTIVE && active_workers_count () < WORKER_THREADS_ACTIVE)
		printf ("\t 4.  Test/Continue or Stop\n");
	else if (!WORKER_THREADS_ACTIVE || WORKER_THREADS_STOPPING)
		printf ("\t 4.  Test/Continue\n");
	else
		printf ("\t 4.  Test/Stop\n");
	printf ("\t 5.  Test/Exit\n");
	printf ("\t 6.  Advanced/Test\n");
	printf ("\t 7.  Advanced/Time\n");
	printf ("\t 8.  Advanced/P-1\n");
	printf ("\t 9.  Advanced/ECM\n");
	printf ("\t10.  Advanced/Manual Communication\n");
	printf ("\t11.  Advanced/Unreserve Exponent\n");
	printf ("\t12.  Advanced/Quit Gimps\n");
	printf ("\t13.  Options/CPU\n");
	printf ("\t14.  Options/Resource Limits\n");
	printf ("\t15.  Options/Preferences\n");
	printf ("\t16.  Options/Torture Test\n");
	printf ("\t17.  Options/Benchmark\n");
	printf ("\t18.  Help/About\n");
	printf ("\t19.  Help/About PrimeNet Server\n");
	printf ("Your choice: ");
	choice = get_number_or_quit (0, 5);
	if (choice <= 0 || choice >= 20) {
		printf ("\n\t     Invalid choice\n\n");
		continue;
	}

/* Display the main menu and switch off the users choice */

	printf ("\n");
	switch (choice) {

/* Test/Primenet dialog */

	case 1:
		test_primenet ();
		break;

/* Test/User Information dialog */

	case 2:
		test_worker_threads ();
		break;

/* Test/Status message */

	case 3:
		test_status ();
		askOK ();
		break;

/* Test/Continue or Stop or Test/Continue or Test/Stop */

	case 4:
		if (WORKER_THREADS_ACTIVE && active_workers_count () < WORKER_THREADS_ACTIVE)
			test_continue_or_stop ();
		else if (NUM_WORKER_THREADS > 1 && active_workers_count () < WORKER_THREADS_ACTIVE - 1)
			test_continue ();
		else if (!WORKER_THREADS_ACTIVE || WORKER_THREADS_STOPPING) {
			while (WORKER_THREADS_STOPPING) Sleep (50);
			linuxContinue ("Another mprime is running.\n", ALL_WORKERS, FALSE);
		} else if (active_workers_count () > 1)
			test_stop ();
		else
			stop_workers_for_escape ();
		break;

/* Test/Exit */

	case 5:
		{
		int counter = 0;
		if (WORKER_THREADS_ACTIVE && !WORKER_THREADS_STOPPING)
			stop_workers_for_escape ();
		while (WORKER_THREADS_STOPPING) {
			if (counter++ % 100 == 0) printf ("Waiting for worker threads to stop.\n");
			Sleep (50);
		}
		}
		return;

/* Advanced/Test dialog */

	case 6:
		advanced_test ();
		break;

/* Advanced/Time dialog */

	case 7:
		advanced_time ();
		break;

/* Advanced/P-1 dialog */

	case 8:
		advanced_pminus1 ();
		break;

/* Advanced/ECM dialog */

	case 9:
		advanced_ecm ();
		break;

/* Advanced/Manual Communication dialog */

	case 10:
		advanced_manualcomm ();
		break;

/* Advanced/Unreserve exponent dialog */

	case 11:
		advanced_unreserve ();
		break;

/* Advanced/Quit Gimps dialog */

	case 12:
		advanced_quit ();
		break;

/* Options/CPU dialog */

	case 13:
		options_cpu ();
		break;

/* Options/Resource limits dialog */

	case 14:
		options_resources ();
		break;

/* Options/Preferences dialog */

	case 15:
		options_preferences ();
		break;

/* Options/Torture Test */

	case 16:
		options_torture ();
		askOK ();
		break;

/* Options/Benchmark Test */

	case 17:
		options_benchmark ();
		break;

/* Help/About */

	case 18:
		help_about ();
		break;

/* Help/About PrimeNet Server */

	case 19:
		help_about_server ();
		break;
	}

	}
}
