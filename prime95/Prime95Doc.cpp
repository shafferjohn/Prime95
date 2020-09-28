// Prime95Doc.cpp : implementation of the CPrime95Doc class
//
// Copyright 1995-2020 Mersenne Research, Inc.  All rights reserved
//

#include "stdafx.h"
#include "MainFrm.h"
#include "Prime95.h"
#include "Prime95Doc.h"
#include "Prime95View.h"

#include <direct.h>
#include "math.h"

#include "BenchmarkDlg.h"
#include "CpuDlg.h"
#include "EcmDlg.h"
#include "ManualCommDlg.h"
#include "Pminus1Dlg.h"
#include "PreferencesDlg.h"
#include "PrimeNetDlg.h"
#include "ResourcesDlg.h"
#include "StartDlg.h"
#include "StopDlg.h"
#include "TestDlg.h"
#include "TimeDlg.h"
#include "TortureDlg.h"
#include "UnreserveDlg.h"
#include "WelcomeDlg.h"
#include "WorkerDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc

IMPLEMENT_DYNCREATE(CPrime95Doc, CDocument)

BEGIN_MESSAGE_MAP(CPrime95Doc, CDocument)
	//{{AFX_MSG_MAP(CPrime95Doc)
	ON_COMMAND(IDM_PRIMENET, OnPrimenet)
	ON_UPDATE_COMMAND_UI(IDM_WORKER_THREADS, OnUpdateWorkerThreads)
	ON_COMMAND(IDM_WORKER_THREADS, OnWorkerThreads)
	ON_UPDATE_COMMAND_UI(IDM_CONTINUE_SWITCHER, OnUpdateContinueSwitcher)
	ON_COMMAND(IDM_CONTINUE_SWITCHER, OnContinueSwitcher)
	ON_COMMAND(IDM_CONTINUE, OnContinue)
	ON_UPDATE_COMMAND_UI(IDM_STOP_SWITCHER, OnUpdateStopSwitcher)
	ON_COMMAND(IDM_STOP_SWITCHER, OnStopSwitcher)
	ON_COMMAND(IDM_STOP, OnStop)
	ON_UPDATE_COMMAND_UI(IDM_ERRCHK, OnUpdateErrchk)
	ON_COMMAND(IDM_ERRCHK, OnErrchk)
	ON_COMMAND(IDM_CPU, OnCpu)
	ON_COMMAND(IDM_RESOURCES, OnResources)
	ON_COMMAND(IDM_PREFERENCES, OnPreferences)
	ON_UPDATE_COMMAND_UI(IDM_TEST, OnUpdateTest)
	ON_COMMAND(IDM_TEST, OnTest)
	ON_UPDATE_COMMAND_UI(IDM_TIME, OnUpdateTime)
	ON_COMMAND(IDM_TIME, OnTime)
	ON_COMMAND(IDM_STATUS, OnRangeStatus)
	ON_UPDATE_COMMAND_UI(ID_HELP_FINDER, OnUpdateHelpFinder)
	ON_COMMAND(IDM_TRAY, OnTray)
	ON_UPDATE_COMMAND_UI(IDM_TRAY, OnUpdateTray)
	ON_COMMAND(IDM_HIDE, OnHide)
	ON_UPDATE_COMMAND_UI(IDM_HIDE, OnUpdateHide)
	ON_COMMAND(IDM_TORTURE, OnTorture)
	ON_UPDATE_COMMAND_UI(IDM_TORTURE, OnUpdateTorture)
	ON_COMMAND(IDM_SERVER, OnServer)
	ON_UPDATE_COMMAND_UI(IDM_SERVER, OnUpdateServer)
	ON_COMMAND(IDM_QUIT, OnQuitGimps)
	ON_COMMAND(IDM_SERVICE, OnService)
	ON_UPDATE_COMMAND_UI(IDM_SERVICE, OnUpdateService)
	ON_COMMAND(IDM_MANUALCOMM, OnManualcomm)
	ON_UPDATE_COMMAND_UI(IDM_MANUALCOMM, OnUpdateManualcomm)
	ON_COMMAND(IDM_ECM, OnEcm)
	ON_UPDATE_COMMAND_UI(IDM_ECM, OnUpdateEcm)
	ON_COMMAND(IDM_PMINUS1, OnPminus1)
	ON_UPDATE_COMMAND_UI(IDM_PMINUS1, OnUpdatePminus1)
	ON_COMMAND(USR_WELCOME, OnWelcome)
	ON_COMMAND(USR_TORTURE, OnUsrTorture)
	ON_COMMAND(IDM_UNRESERVE, OnUnreserve)
	ON_UPDATE_COMMAND_UI(IDM_UNRESERVE, OnUpdateUnreserve)
	ON_UPDATE_COMMAND_UI(IDM_QUIT, OnUpdateQuit)
	ON_COMMAND(IDM_BENCHMARK, OnBenchmark)
	ON_UPDATE_COMMAND_UI(IDM_BENCHMARK, OnUpdateBenchmark)
	ON_COMMAND(IDM_MERGE_MAIN, OnMergeMain)
	ON_UPDATE_COMMAND_UI(IDM_MERGE_MAIN, OnUpdateMergeMain)
	ON_COMMAND(IDM_MERGE_COMM, OnMergeComm)
	ON_UPDATE_COMMAND_UI(IDM_MERGE_COMM, OnUpdateMergeComm)
	ON_COMMAND(IDM_MERGE_ALL, OnMergeAll)
	ON_UPDATE_COMMAND_UI(IDM_MERGE_ALL, OnUpdateMergeAll)
	ON_COMMAND(IDM_HELP_FORUM, OnForum)
	ON_COMMAND(IDM_HELP_WIKI, OnWiki)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc construction/destruction

CPrime95Doc::CPrime95Doc()
{
}

CPrime95Doc::~CPrime95Doc()
{
}

BOOL CPrime95Doc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc serialization

void CPrime95Doc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc diagnostics

#ifdef _DEBUG
void CPrime95Doc::AssertValid() const
{
	CDocument::AssertValid();
}

void CPrime95Doc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc commands

void CPrime95Doc::OnCloseDocument() 
{

// Set flag indicating we are exiting.  Needed because setting the icon
// while in the sleep 50 loop causes a hang.

	EXIT_IN_PROGRESS = 1;

// Stop background threads before exiting

	if (WORKER_THREADS_ACTIVE) {
		OnStop ();
		while (WORKER_THREADS_STOPPING) Sleep (50);
	}

// Remember the main window's size and position
// Note: INI_FILE may not be initialized when this instance
// was used to activate another instance

	CWinApp* pApp = AfxGetApp();
	WINDOWPLACEMENT wp;
	if (pApp->m_pMainWnd && INI_FILE[0]) {
		pApp->m_pMainWnd->GetWindowPlacement (&wp);
		IniWriteInt (INI_FILE, "Left", wp.rcNormalPosition.left);
		IniWriteInt (INI_FILE, "Top", wp.rcNormalPosition.top);
		IniWriteInt (INI_FILE, "Right", wp.rcNormalPosition.right);
		IniWriteInt (INI_FILE, "Bottom", wp.rcNormalPosition.bottom);
	}

// Free the networking library

	UnloadPrimeNet ();
	SaveViews();

/* Write the worktodo file in case the WELL_BEHAVED_WORK flag caused us */
/* to delay writing the file. */

	writeWorkToDoFile (TRUE);

// Finish closing

	CDocument::OnCloseDocument();
}

// Test menu

void getProxyInfo (char *, unsigned short *, char *, char *);

void CPrime95Doc::OnPrimenet() 
{
	int	update_computer_info, primenet_debug;

	update_computer_info = FALSE;
	primenet_debug = IniSectionGetInt (INI_FILE, "PrimeNet", "Debug", 0);

	PrimenetDlg dlg;
	char	szProxyHost[120], szProxyUser[50], szProxyPassword[50];
	unsigned short nProxyPort;

	dlg.m_primenet = USE_PRIMENET;
	if (strcmp (USERID, "ANONYMOUS") == 0)
		dlg.m_userid = "";
	else
		dlg.m_userid = USERID;
	dlg.m_compid = COMPID;
	dlg.m_dialup = DIAL_UP;
	dlg.m_debug = primenet_debug;
	getProxyInfo (szProxyHost, &nProxyPort, szProxyUser, szProxyPassword);
	if (szProxyHost[0]) {
		dlg.m_proxyhost = szProxyHost;
		dlg.m_proxyport = nProxyPort;
		dlg.m_proxyuser = szProxyUser;
		dlg.m_proxypassword = szProxyPassword;
	}
	if (dlg.DoModal () == IDOK) {
		DIAL_UP = dlg.m_dialup;
		IniWriteInt (INI_FILE, "DialUp", DIAL_UP);
		strcpy (szProxyHost, (const char *) dlg.m_proxyhost);
		if (szProxyHost[0] && dlg.m_proxyport != 8080)
			sprintf (szProxyHost + strlen (szProxyHost), ":%d", dlg.m_proxyport);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyHost", szProxyHost);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyUser", dlg.m_proxyuser);
		if (strcmp (szProxyPassword, dlg.m_proxypassword)) {
			IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyPass", dlg.m_proxypassword);
			IniSectionWriteInt (INI_FILE, "PrimeNet", "ProxyMask", 0);
		}
		if (!dlg.m_debug != !primenet_debug) {
			IniSectionWriteInt (INI_FILE, "PrimeNet", "Debug", dlg.m_debug ? 1 : 0);
		}

		if (dlg.m_userid[0] == 0)
			dlg.m_userid = "ANONYMOUS";

		if (strcmp (USERID, dlg.m_userid) != 0) {
			strcpy (USERID, (const char *) dlg.m_userid);
			sanitizeString (USERID);
			IniWriteString (INI_FILE, "V5UserID", USERID);
			update_computer_info = TRUE;
		}
		if (strcmp (COMPID, dlg.m_compid) != 0) {
			strcpy (COMPID, (const char *) dlg.m_compid);
			sanitizeString (COMPID);
			IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
			update_computer_info = TRUE;
		}
		if (!USE_PRIMENET && dlg.m_primenet) {
			USE_PRIMENET = 1;
			create_window (COMM_THREAD_NUM);
			base_title (COMM_THREAD_NUM, "Communication thread");
			if (!STARTUP_IN_PROGRESS) set_comm_timers ();
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
			spoolExistingResultsFile ();
		} else if (USE_PRIMENET && !dlg.m_primenet) {
			USE_PRIMENET = 0;
			if (!STARTUP_IN_PROGRESS) set_comm_timers ();
		} else if (update_computer_info)
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);

		IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET);

/* For historical reasons, this dialog box also does a Test/Continue */
/* when you are using primenet */

		if (!STARTUP_IN_PROGRESS && USE_PRIMENET) OnContinue ();
	} else
		STARTUP_IN_PROGRESS = 0;
}

void CPrime95Doc::OnUpdateQuit(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET || WORKTODO_COUNT);
}

void CPrime95Doc::OnQuitGimps() 
{
	if (!USE_PRIMENET) {
		if (AfxMessageBox (MANUAL_QUIT, MB_YESNO | MB_ICONQUESTION | MB_DEFBUTTON2) == IDYES) {
			OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS.\n");
//bug - either delete file, or delete all work_units and write the file.
//bug			IniDeleteAllLines (WORKTODO_FILE);
			stop_workers_for_escape ();
		}
	} else {
		int	res;
		res = AfxMessageBox (PRIMENET_QUIT, MB_YESNOCANCEL | MB_ICONQUESTION | MB_DEFBUTTON3);
		if (res == IDYES) {
			OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS after current work completes.\n");
			IniWriteInt (INI_FILE, "NoMoreWork", 1);
		}
		if (res == IDNO) {
			OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS immediately.\n");
			spoolMessage (MSG_QUIT_GIMPS, NULL);
		}
	}
}

void CPrime95Doc::OnUpdateWorkerThreads(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (1);
}

void CPrime95Doc::OnWorkerThreads() 
{
	CWorkerDlg dlg;
	int	i;

	dlg.m_num_thread = NUM_WORKER_THREADS;
	for (i = 0; i < MAX_NUM_WORKER_THREADS; i++) {
		dlg.m_work_pref[i] = WORK_PREFERENCE[i];
		dlg.m_numcpus[i] = CORES_PER_TEST[i];
	}
	dlg.m_cert_work = IniGetInt (LOCALINI_FILE, "CertWork", 1);

again:	if (dlg.DoModal () == IDOK) {
		int	restart = FALSE;
		int	new_options = FALSE;
		unsigned long i, total_num_cores;

/* If the user has selected 100M tests and per-worker temp disk is not enough for a power=8 proof, then do not permit it. */

		if (CPU_WORKER_DISK_SPACE < 12.0) {
			int	changed = FALSE;
			for (i = 0; i < dlg.m_num_thread; i++) {
				if (dlg.m_work_pref[i] == PRIMENET_WP_PRP_100M) {
					dlg.m_work_pref[i] = PRIMENET_WP_PRP_FIRST;
					changed = TRUE;
				}
			}
			if (changed)
				AfxMessageBox (MSG_100M, MB_ICONEXCLAMATION | MB_OK);
		}

/* If the user has selected first-time tests and per-worker temp disk is not enough for a power=6 proof, then warn the user. */

		if (CPU_WORKER_DISK_SPACE < 1.5) {
			int	warn = FALSE;
			for (i = 0; i < dlg.m_num_thread; i++) {
				if (dlg.m_work_pref[i] == PRIMENET_WP_PRP_FIRST || dlg.m_work_pref[i] == PRIMENET_WP_PRP_WORLD_RECORD) {
					warn = TRUE;
				}
			}
			if (warn)
				AfxMessageBox (MSG_FIRST, MB_ICONEXCLAMATION | MB_OK);
		}

/* If the user has allocated too many cores then raise a severe warning. */

		total_num_cores = 0;
		for (i = 0; i < dlg.m_num_thread; i++) total_num_cores += dlg.m_numcpus[i];
		if (total_num_cores > NUM_CPUS &&
		    AfxMessageBox (MSG_THREADS, MB_YESNO | MB_ICONQUESTION) == IDYES)
			goto again;

/* If user changed the number of worker threads, then make the */
/* necessary changes.  Restart worker threads so that we are running */
/* the correct number of worker threads. */

		if (dlg.m_num_thread != NUM_WORKER_THREADS) {
			NUM_WORKER_THREADS = dlg.m_num_thread;
			IniWriteInt (LOCALINI_FILE, "WorkerThreads", NUM_WORKER_THREADS);
			new_options = TRUE;
			restart = TRUE;
		}

/* If the user changed any of the work preferences record it in the INI file */
/* and tell the server */

		if (dlg.AreAllTheSame (dlg.m_work_pref)) {
			if (! PTOIsGlobalOption (WORK_PREFERENCE) || WORK_PREFERENCE[0] != dlg.m_work_pref[0]) {
				PTOSetAll (INI_FILE, "WorkPreference", NULL, WORK_PREFERENCE, dlg.m_work_pref[0]);
				new_options = TRUE;
			}
		} else {
			for (i = 0; i < (int) NUM_WORKER_THREADS; i++) {
				if (WORK_PREFERENCE[i] == dlg.m_work_pref[i]) continue;
				PTOSetOne (INI_FILE, "WorkPreference", NULL, WORK_PREFERENCE, i, dlg.m_work_pref[i]);
				new_options = TRUE;
			}
		}

/* If user changed any of the cores_per_test, then record it in the INI file */

		if (dlg.AreAllTheSame (dlg.m_numcpus))
			PTOSetAll (LOCALINI_FILE, "CoresPerTest", NULL, CORES_PER_TEST, dlg.m_numcpus[0]);
		else for (i = 0; i < (int) NUM_WORKER_THREADS; i++)
			PTOSetOne (LOCALINI_FILE, "CoresPerTest", NULL, CORES_PER_TEST, i, dlg.m_numcpus[i]);

/* Write the new CertWork setting */

		IniWriteInt (LOCALINI_FILE, "CertWork", dlg.m_cert_work);

/* Send new settings to the server */

		if (new_options) spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);

/* Restart worker threads with new options.  Since Windows must create */
/* worker windows in the main thread, the routine we call will do that */
/* if there are now more worker threads. */

		if (restart) stop_workers_for_restart ();
	} else
		STARTUP_IN_PROGRESS = 0;
}

void CPrime95Doc::OnRangeStatus() 
{
	rangeStatus ();
}

void CPrime95Doc::OnUpdateContinueSwitcher(CCmdUI* pCmdUI) 
{
	pCmdUI->SetText (((!WORKER_THREADS_ACTIVE && NUM_WORKER_THREADS > 1) ||
			  (WORKER_THREADS_ACTIVE && active_workers_count () != WORKER_THREADS_ACTIVE - 1)) ? "&Continue..." : "&Continue");
	pCmdUI->Enable ((!WORKER_THREADS_ACTIVE && (USE_PRIMENET || WORKTODO_COUNT)) ||
			(WORKER_THREADS_ACTIVE && active_workers_count () != WORKER_THREADS_ACTIVE));
}

void CPrime95Doc::OnContinueSwitcher() 
{
	if ((!WORKER_THREADS_ACTIVE && NUM_WORKER_THREADS > 1) ||
	    (WORKER_THREADS_ACTIVE && active_workers_count () != WORKER_THREADS_ACTIVE - 1)) {
		// Start the dialog box
		CStartDlg dlg;

		if (dlg.DoModal () == IDOK) {
			if (dlg.m_all_threads)
				OnContinue ();
			else
				LaunchWorkerThreads (dlg.m_thread-1, FALSE);
		}
	} else {
		// Start the thread
		OnContinue ();
	}
}

void CPrime95Doc::OnContinue() 
{
	// Start the threads
	LaunchWorkerThreads (ALL_WORKERS, FALSE);
}

void CPrime95Doc::OnUpdateStopSwitcher(CCmdUI* pCmdUI) 
{
	// Set text to "Stop..." if multiple workers or torture threads are running
	pCmdUI->SetText (active_workers_count () > 1 ? "St&op..." : "St&op");
	pCmdUI->Enable (WORKER_THREADS_ACTIVE && !WORKER_THREADS_STOPPING);
}

void CPrime95Doc::OnStopSwitcher() 
{
	if (active_workers_count () > 1) {
		// Start the dialog box
		CStopDlg dlg;

		if (dlg.DoModal () == IDOK) {
			if (dlg.m_all_threads)
				OnStop ();
			else
				stop_one_worker (dlg.m_thread-1);
		}
	} else {
		// Stop the one active thread
		OnStop ();
	}
}

void CPrime95Doc::OnStop() 
{
	stop_workers_for_escape ();
}

// Advanced Menu

void CPrime95Doc::OnUpdateTest(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (TRUE);
}

void CPrime95Doc::OnTest() 
{
	CTestDlg dlg;

	if (dlg.DoModal () == IDOK) {
		struct work_unit w;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_ADVANCEDTEST;
		w.k = 1.0;
		w.b = 2;
		w.n = dlg.m_p;
		w.c = -1;
		addWorkToDoLine (dlg.m_thread - 1, &w);
		if (WORKER_THREADS_ACTIVE)
			stop_worker_for_advanced_test (dlg.m_thread - 1);
		else
			OnContinue ();
	}
}

void CPrime95Doc::OnUpdateTime(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (!WORKER_THREADS_STOPPING);
}

void CPrime95Doc::OnTime() 
{
	CTimeDlg dlg;

	dlg.m_p = 38000000;
	dlg.m_iter = 10;
	if (dlg.DoModal () == IDOK) {
		LaunchAdvancedTime (dlg.m_p, dlg.m_iter);
	}
}

void CPrime95Doc::OnUpdatePminus1(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (TRUE);
}

void CPrime95Doc::OnPminus1() 
{
	CPminus1Dlg dlg;

	dlg.m_k = 1.0;
	dlg.m_b = 2;
	dlg.m_n = 1061;
	dlg.m_c = -1;
	dlg.m_bound1 = 1000000;
	if (dlg.DoModal () == IDOK) {
		struct work_unit w;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_PMINUS1;
		w.k = dlg.m_k;
		w.b = dlg.m_b;
		w.n = dlg.m_n;
		w.c = dlg.m_c;
		w.B1 = dlg.m_bound1;
		w.B2_start = 0;
		w.B2 = dlg.m_bound2;
		addWorkToDoLine (dlg.m_thread - 1, &w);

/* If worker threads are running, adding the work should have restarted */
/* threads waiting for work.  Otherwise, start the worker threads. */

		if (!WORKER_THREADS_ACTIVE) OnContinue ();
	}
}

void CPrime95Doc::OnUpdateEcm(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (TRUE);
}

void CPrime95Doc::OnEcm() 
{
	CEcmDlg dlg;

	dlg.m_k = 1.0;
	dlg.m_b = 2;
	dlg.m_n = 1061;
	dlg.m_c = -1;
	dlg.m_bound1 = 1000000;
	dlg.m_num_curves = 100;
	if (dlg.DoModal () == IDOK) {
		struct work_unit w;
		memset (&w, 0, sizeof (w));
		w.work_type = WORK_ECM;
		w.k = dlg.m_k;
		w.b = dlg.m_b;
		w.n = dlg.m_n;
		w.c = dlg.m_c;
		w.B1 = dlg.m_bound1;
		w.B2 = dlg.m_bound2;
		w.curves_to_do = dlg.m_num_curves;
		addWorkToDoLine (dlg.m_thread - 1, &w);

/* If worker threads are running, adding the work should have restarted */
/* threads waiting for work.  Otherwise, start the worker threads. */

		if (!WORKER_THREADS_ACTIVE) OnContinue ();
	}
}

void CPrime95Doc::OnUpdateErrchk(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck (ERRCHK);
	pCmdUI->Enable (1);
}

void CPrime95Doc::OnErrchk() 
{
	ERRCHK = !ERRCHK;
	IniWriteInt (INI_FILE, "ErrorCheck", ERRCHK);
}

void CPrime95Doc::OnUpdateManualcomm(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET);
}

void CPrime95Doc::OnManualcomm() 
{
	CManualCommDlg dlg;

	dlg.m_manual_comm = MANUAL_COMM;
	dlg.m_comm_now = 1;
	dlg.m_new_dates = 0;
	if (dlg.DoModal () == IDOK) {
		if ((MANUAL_COMM && !dlg.m_manual_comm) ||
		    (!MANUAL_COMM && dlg.m_manual_comm)) {
			MANUAL_COMM = dlg.m_manual_comm;
			IniWriteInt (INI_FILE, "ManualComm", MANUAL_COMM);
			set_comm_timers ();
		}
		if (dlg.m_new_dates) UpdateEndDates ();
		if (dlg.m_comm_now) do_manual_comm_now ();
	}
}

void CPrime95Doc::OnUpdateUnreserve(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET);
}

void CPrime95Doc::OnUnreserve() 
{
	CUnreserveDlg dlg;

	if (dlg.DoModal () == IDOK) unreserve (dlg.m_p);
}


// Options menu

void CPrime95Doc::OnCpu() 
{
	CCpuDlg dlg;
	char	buf[512];

	dlg.m_hours = CPU_HOURS;
	getCpuDescription (buf, 0);
	dlg.m_cpu_info = buf;
//again:
	if (dlg.DoModal () == IDOK) {

		if (CPU_HOURS != dlg.m_hours) {
			CPU_HOURS = dlg.m_hours;
			IniWriteInt (LOCALINI_FILE, "CPUHours", CPU_HOURS);
			ROLLING_AVERAGE = 1000;
			IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
			IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
			spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
			delete_timed_event (TE_COMM_SERVER);
			UpdateEndDates ();
		}

// Now that Primenet almost always hands out LL assignments that are P-1'ed,
// there is little reason to prompt user into allowing us to use more memory.
//		if (!IniGetInt (INI_FILE, "AskedAboutMemory", 0)) {
//			IniWriteInt (INI_FILE, "AskedAboutMemory", 1);
//			if (dlg.m_day_memory == 8 && dlg.m_night_memory == 8 &&
//			    AfxMessageBox (MSG_MEMORY, MB_YESNO | MB_ICONQUESTION) == IDYES)
//				goto again;
//		}
	} else
		STARTUP_IN_PROGRESS = 0;
}

void CPrime95Doc::OnResources() 
{
	CResourcesDlg dlg;
	char	timebuf[20];

	dlg.m_disk = CPU_WORKER_DISK_SPACE;
	dlg.m_upload_bandwidth = IniSectionGetFloat (INI_FILE, "PrimeNet", "UploadRateLimit", 0.25);
	if (dlg.m_upload_bandwidth <= 0.0 || dlg.m_upload_bandwidth > 10000.0) dlg.m_upload_bandwidth = 10000.0;
	IniSectionGetString (INI_FILE, "PrimeNet", "UploadStartTime", timebuf, sizeof (timebuf), "00:00");
	if (strcmp (timebuf, "00:00") != 0) minutesToStr (strToMinutes (timebuf), timebuf);
	dlg.m_upload_start = timebuf;
	IniSectionGetString (INI_FILE, "PrimeNet", "UploadEndTime", timebuf, sizeof (timebuf), "24:00");
	if (strcmp (timebuf, "24:00") != 0) minutesToStr (strToMinutes (timebuf), timebuf);
	dlg.m_upload_end = timebuf;
	dlg.m_download_mb = IniSectionGetInt (INI_FILE, "PrimeNet", "DownloadDailyLimit", 40);
	dlg.m_can_upload = IniSectionGetInt (INI_FILE, "PrimeNet", "ProofUploads", 1);
	if (dlg.DoModal () == IDOK) {
		// Raise a warning if uesr drops the temp disk space below the threshold for first time work.
		if (CPU_WORKER_DISK_SPACE >= 1.5 && dlg.m_disk < 1.5) {
			AfxMessageBox (MSG_DISK, MB_ICONEXCLAMATION | MB_OK);
		}
		CPU_WORKER_DISK_SPACE = dlg.m_disk;
		IniWriteFloat (LOCALINI_FILE, "WorkerDiskSpace", CPU_WORKER_DISK_SPACE);
		IniSectionWriteFloat (INI_FILE, "PrimeNet", "UploadRateLimit", dlg.m_upload_bandwidth);
		IniSectionWriteString (INI_FILE, "PrimeNet", "UploadStartTime", (const char *) dlg.m_upload_start);
		IniSectionWriteString (INI_FILE, "PrimeNet", "UploadEndTime", (const char *) dlg.m_upload_end);
		IniSectionWriteInt (INI_FILE, "PrimeNet", "DownloadDailyLimit", dlg.m_download_mb);
	}
}

void CPrime95Doc::OnPreferences() 
{
	CPreferencesDlg dlg;

	dlg.m_iter = ITER_OUTPUT;
	dlg.m_r_iter = ITER_OUTPUT_RES;
	dlg.m_disk_write_time = DISK_WRITE_TIME;
	dlg.m_modem = MODEM_RETRY_TIME;
	dlg.m_retry = NETWORK_RETRY_TIME;
	dlg.m_work = DAYS_OF_WORK;
	dlg.m_end_dates = DAYS_BETWEEN_CHECKINS;
	dlg.m_backup = NUM_BACKUP_FILES;
	dlg.m_noise = !SILENT_VICTORY;
	dlg.m_battery = RUN_ON_BATTERY;
	if (dlg.DoModal () == IDOK) {
		ITER_OUTPUT = dlg.m_iter;
		ITER_OUTPUT_RES = dlg.m_r_iter;
		DISK_WRITE_TIME = dlg.m_disk_write_time;
		MODEM_RETRY_TIME = dlg.m_modem;
		NETWORK_RETRY_TIME = dlg.m_retry;
		DAYS_OF_WORK = dlg.m_work;
		DAYS_BETWEEN_CHECKINS = dlg.m_end_dates;
		NUM_BACKUP_FILES = dlg.m_backup;
		SILENT_VICTORY = !dlg.m_noise;
		if (RUN_ON_BATTERY != dlg.m_battery) {
			RUN_ON_BATTERY = dlg.m_battery;
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

void CPrime95Doc::OnUpdateBenchmark(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (!WORKER_THREADS_STOPPING);
}

void CPrime95Doc::OnBenchmark() 
{
	CBenchmarkDlg dlg;
	char	default_cores_string[80], default_workers_string[80];
	int	i, vals[4], numvals;

	// Init FFT size dialog box entries
	dlg.m_minFFT = IniGetInt (INI_FILE, "MinBenchFFT", 2048);
	dlg.m_maxFFT = IniGetInt (INI_FILE, "MaxBenchFFT", 8192);
	dlg.m_errchk = ERRCHK;				// IniGetInt (INI_FILE, "BenchErrorCheck", 0);
	dlg.m_all_complex = 0;				// IniGetInt (INI_FILE, "BenchAllComplex", 0);
	dlg.m_limit_FFT_sizes = 0;			// IniGetInt (INI_FILE, "OnlyBench5678", 1);

	// Init CPU cores dialog box entries
	sprintf (default_cores_string, "%lu", NUM_CPUS);
	dlg.m_bench_cores = default_cores_string;
	dlg.m_hyperthreading = IniGetInt (INI_FILE, "BenchHyperthreads", 1);

	// Init throughput dialog box entries
	dlg.m_all_FFT_impl = IniGetInt (INI_FILE, "AllBench", 0);
	dlg.m_bench_time = IniGetInt (INI_FILE, "BenchTime", 15);
	// If testing all FFT implementations. then default to the current num_workers.
	// Otherwise, assume user is trying to figure out how many workers to run and form a string
	// with the most common best values for number of workers: 1, num_threading_nodes, num_cores, num_workers
	numvals = 0;
	sorted_add_unique (vals, &numvals, NUM_WORKER_THREADS);
	if (!dlg.m_all_FFT_impl) {
		sorted_add_unique (vals, &numvals, 1);
		sorted_add_unique (vals, &numvals, NUM_THREADING_NODES);
		sorted_add_unique (vals, &numvals, NUM_CPUS);
	}
	sprintf (default_workers_string, "%d", vals[0]);
	for (i = 1; i < numvals; i++) sprintf (default_workers_string + strlen (default_workers_string), ",%d", vals[i]);
	dlg.m_bench_workers = default_workers_string;

	if (dlg.DoModal () == IDOK) {
		if (dlg.m_bench_type != 2) {
			IniWriteInt (INI_FILE, "MinBenchFFT", dlg.m_minFFT);
			IniWriteInt (INI_FILE, "MaxBenchFFT", dlg.m_maxFFT);
			IniWriteInt (INI_FILE, "BenchErrorCheck", dlg.m_errchk);
			IniWriteInt (INI_FILE, "BenchAllComplex", dlg.m_all_complex ? 2 : 0);
			IniWriteInt (INI_FILE, "OnlyBench5678", dlg.m_limit_FFT_sizes);
		}
		IniWriteString (INI_FILE, "BenchCores", dlg.m_bench_cores);
		IniWriteInt (INI_FILE, "BenchHyperthreads", dlg.m_hyperthreading);
		if (dlg.m_bench_type == 0) {
			IniWriteString (INI_FILE, "BenchWorkers", dlg.m_bench_workers);
			IniWriteInt (INI_FILE, "AllBench", dlg.m_all_FFT_impl);
			IniWriteInt (INI_FILE, "BenchTime", dlg.m_bench_time);
		}
		LaunchBench (dlg.m_bench_type);
	}
}

void CPrime95Doc::OnUpdateTorture(CCmdUI* pCmdUI)
{
	pCmdUI->Enable (!WORKER_THREADS_STOPPING);
}

void CPrime95Doc::OnTorture() 
{
	CTortureDlg dlg;
	int	mem;

	dlg.m_minfft = 4;
	dlg.m_maxfft = (CPU_TOTAL_L4_CACHE_SIZE ? 32768 : 8192);
	dlg.m_thread = NUM_CPUS * CPU_HYPERTHREADS;
	mem = physical_memory ();
	// New in 29.5, raise the default memory used to all but 3GB on 64-bit machines.  Almost all machines today have
	// more than 5GB of memory installed.  If memory serves me, there may be issues in Win32 allocating more than 2GB.
#ifdef X86_64
	if (mem >= 5000) {
		dlg.m_blendmemory = GetSuggestedMemory (mem - 3000);
		dlg.m_in_place = FALSE;
	} else
#endif
	// These are the pre 29.5 memory defaults.
	if (mem >= 3000) {
		dlg.m_blendmemory = GetSuggestedMemory (2000);
		dlg.m_in_place = FALSE;
	} else if (mem >= 2000) {
		dlg.m_blendmemory = GetSuggestedMemory (1500);
		dlg.m_in_place = FALSE;
	} else if (mem >= 500) {
		dlg.m_blendmemory = GetSuggestedMemory (mem - 256);
		dlg.m_in_place = FALSE;
	} else if (mem >= 200) {
		dlg.m_blendmemory = GetSuggestedMemory (mem / 2);
		dlg.m_in_place = TRUE;
	} else {
		dlg.m_blendmemory = 8;
		dlg.m_in_place = TRUE;
	}
	dlg.m_memory = dlg.m_blendmemory;
	dlg.m_timefft = CPU_HYPERTHREADS * 3;
	if (dlg.DoModal () == IDOK) {
		int	m_weak;
		IniWriteInt (INI_FILE, "MinTortureFFT", dlg.m_minfft);
		IniWriteInt (INI_FILE, "MaxTortureFFT", dlg.m_maxfft);
		if (dlg.m_in_place) dlg.m_memory = 8;
		IniWriteInt (INI_FILE, "TortureMem", dlg.m_memory);
		IniWriteInt (INI_FILE, "TortureTime", dlg.m_timefft);
		m_weak = dlg.m_avx512 * CPU_AVX512F + dlg.m_fma3 * CPU_FMA3 + dlg.m_avx * CPU_AVX + dlg.m_sse2 * CPU_SSE2;
		IniWriteInt (INI_FILE, "TortureWeak", m_weak);
		LaunchTortureTest (dlg.m_thread, FALSE);
	}
}

void CPrime95Doc::OnUpdateTray(CCmdUI* pCmdUI) 
{
	pCmdUI->SetCheck (TRAY_ICON);
}

void CPrime95Doc::OnTray() 
{
	CPrime95App* pApp = (CPrime95App *) AfxGetApp();
	TRAY_ICON = ! TRAY_ICON;
	if (TRAY_ICON) {
		HIDE_ICON = 0;
		pApp->TrayMessage (NIM_ADD, NULL, 0);
		ChangeIcon (MAIN_THREAD_NUM, -1);
	} else {
		pApp->TrayMessage (NIM_DELETE, NULL, 0);
	}
	IniWriteInt (INI_FILE, "HideIcon", HIDE_ICON);
	IniWriteInt (INI_FILE, "TrayIcon", TRAY_ICON);
}

void CPrime95Doc::OnUpdateHide(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (1);
	pCmdUI->SetCheck (HIDE_ICON);
}

void CPrime95Doc::OnHide() 
{
	CPrime95App* pApp = (CPrime95App *) AfxGetApp();

	HIDE_ICON = ! HIDE_ICON;
	if (HIDE_ICON) {
		if (TRAY_ICON) pApp->TrayMessage (NIM_DELETE, NULL, 0);
		TRAY_ICON = 0;
	}
	IniWriteInt (INI_FILE, "HideIcon", HIDE_ICON);
	IniWriteInt (INI_FILE, "TrayIcon", TRAY_ICON);
}

// Check the menu item if start at logon is enabled

static int startup_at_logon_set;

void CPrime95Doc::OnUpdateService(CCmdUI* pCmdUI)
{
	char	pathname[512];
	DWORD	rc, dwtype;
	DWORD	len = sizeof(pathname);

/* See if there is an entry for prime95 */

	rc = RegGetValueA (HKEY_CURRENT_USER, "Software\\Microsoft\\Windows\\CurrentVersion\\Run", "Prime95", RRF_RT_REG_SZ, &dwtype, pathname, &len);
	startup_at_logon_set = (rc == ERROR_SUCCESS && pathname[0]);
	pCmdUI->Enable (1);
	pCmdUI->SetCheck (startup_at_logon_set);
}

void CPrime95Doc::OnService() 
{
	char	regkey[20];
	char	pathname[512];
	HKEY	hkey;
	DWORD	rc, disposition;

/* Get handle to the registry entry */

	if (RegCreateKeyEx (
			HKEY_CURRENT_USER,
			"Software\\Microsoft\\Windows\\CurrentVersion\\Run",
			0,
			NULL,
			REG_OPTION_NON_VOLATILE,
			KEY_ALL_ACCESS,
			NULL,
			&hkey,
			&disposition) != ERROR_SUCCESS) {
		OutputStr (MAIN_THREAD_NUM, "Can't access Run registry key.\n");
		return;
	}

/* Now create or delete an entry for prime95 */

	strcpy (regkey, "Prime95");
	startup_at_logon_set = !startup_at_logon_set;
	if (startup_at_logon_set) {
		GetModuleFileName (NULL, pathname, sizeof (pathname));		/* Get pathname of executable */
		rc = RegSetValueEx (hkey, regkey, 0, REG_SZ,
				(BYTE *) pathname, (DWORD) strlen (pathname) + 1);
		if (rc != ERROR_SUCCESS) {
			OutputStr (MAIN_THREAD_NUM, "Can't write Run registry entry.\n");
			return;
		}
	} else {
		rc = RegDeleteValue (hkey, regkey);
		if (rc != ERROR_SUCCESS && rc != ERROR_FILE_NOT_FOUND){
			OutputStr (MAIN_THREAD_NUM, "Can't delete Run registry entry.\n");
			return;
		}
	}
	RegCloseKey (hkey);
}

// Window menu

// This menu choice combines and uncombines the main and comm windows.

void CPrime95Doc::OnUpdateMergeMain(CCmdUI* pCmdUI)
{
	pCmdUI->Enable (1);
	pCmdUI->SetCheck (MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS);
}

void CPrime95Doc::OnMergeMain() 
{
	// If going from checked to unchecked state, create the main window
	if (MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS) {
		MERGE_WINDOWS &= ~MERGE_MAINCOMM_WINDOWS;
		MERGE_WINDOWS &= ~MERGE_MAIN_WINDOW;
		MERGE_WINDOWS &= ~MERGE_COMM_WINDOW;
		create_window (MAIN_THREAD_NUM);
		base_title (MAIN_THREAD_NUM, "Main thread");
	}
	// If going from unchecked to checked state, destroy the main window
	else {
		int	destroy;
		destroy = ! (MERGE_WINDOWS & MERGE_MAIN_WINDOW);
		if (destroy) destroy_window (MAIN_THREAD_NUM);
		MERGE_WINDOWS |= MERGE_MAINCOMM_WINDOWS;
		MERGE_WINDOWS &= ~MERGE_MAIN_WINDOW;
		MERGE_WINDOWS &= ~MERGE_COMM_WINDOW;
	}
	// In both checked and unchecked states we need to make sure the comm window exists
	create_window (COMM_THREAD_NUM);
	base_title (COMM_THREAD_NUM, "Communication thread");
	IniWriteInt (INI_FILE, "MergeWindows", MERGE_WINDOWS);
	PositionViews (TRUE);
}

// This menu choice, in checked state, combines the main and comm and 1st worker windows.
// In going to unchecked state we leave the main and comm windows combined.

void CPrime95Doc::OnUpdateMergeComm(CCmdUI* pCmdUI)
{
	pCmdUI->SetText (NUM_WORKER_THREADS == 1 ?
				"Merge Main && Comm && &Worker" :
			 MERGE_WINDOWS & MERGE_WORKER_WINDOWS ?
				"Merge Main && Comm && &Workers" :
				"Merge Main && Comm && 1st &Worker");
	pCmdUI->Enable (1);
	pCmdUI->SetCheck (MERGE_WINDOWS & MERGE_MAIN_WINDOW && MERGE_WINDOWS & MERGE_COMM_WINDOW);
}

void CPrime95Doc::OnMergeComm() 
{
	// If going from checked to unchecked state, create the comm window
	if (MERGE_WINDOWS & MERGE_MAIN_WINDOW && MERGE_WINDOWS & MERGE_COMM_WINDOW) {
		MERGE_WINDOWS |= MERGE_MAINCOMM_WINDOWS;
		MERGE_WINDOWS &= ~MERGE_MAIN_WINDOW;
		MERGE_WINDOWS &= ~MERGE_COMM_WINDOW;
		create_window (COMM_THREAD_NUM);
		base_title (COMM_THREAD_NUM, "Communication thread");
	}
	// If going from unchecked to checked state, destroy the main and comm windows
	else {
		int	destroy_main, destroy_comm;
		destroy_main = ! (MERGE_WINDOWS & MERGE_MAIN_WINDOW) && ! (MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS);
		destroy_comm = ! (MERGE_WINDOWS & MERGE_COMM_WINDOW);
		if (destroy_main) destroy_window (MAIN_THREAD_NUM);
		if (destroy_comm) destroy_window (COMM_THREAD_NUM);
		MERGE_WINDOWS &= ~MERGE_MAINCOMM_WINDOWS;
		MERGE_WINDOWS |= MERGE_MAIN_WINDOW;
		MERGE_WINDOWS |= MERGE_COMM_WINDOW;
	}
	IniWriteInt (INI_FILE, "MergeWindows", MERGE_WINDOWS);
	PositionViews (TRUE);
}

void CPrime95Doc::OnUpdateMergeAll(CCmdUI* pCmdUI)
{
//	pCmdUI->Enable (NUM_WORKER_THREADS > 1);
	pCmdUI->Enable (TRUE);					/* Stress testers may want to set this option */
	pCmdUI->SetCheck (MERGE_WINDOWS & MERGE_WORKER_WINDOWS);
}

void CPrime95Doc::OnMergeAll() 
{
	int	i;

	if (! (MERGE_WINDOWS & MERGE_WORKER_WINDOWS)) {
		for (i = 1; i < MAX_NUM_WORKER_THREADS; i++)
			destroy_window (i);
	}
	MERGE_WINDOWS ^= MERGE_WORKER_WINDOWS;
	IniWriteInt (INI_FILE, "MergeWindows", MERGE_WINDOWS);
	if (! (MERGE_WINDOWS & MERGE_WORKER_WINDOWS)) {
		create_worker_windows (NUM_WORKER_THREADS);
	}
	PositionViews (TRUE);
}

// Help menu

void CPrime95Doc::OnUpdateHelpFinder(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (TRUE);
}

void CPrime95Doc::OnForum() 
{
	CHyperLink dummy;
	dummy.GotoURL (_T("http://mersenneforum.org"), SW_SHOW);
}

void CPrime95Doc::OnWiki() 
{
	CHyperLink dummy;
	dummy.GotoURL (_T("https://rieselprime.de/ziki/Main_Page"), SW_SHOW);
}

void CPrime95Doc::OnUpdateServer(CCmdUI* pCmdUI) 
{
	pCmdUI->Enable (USE_PRIMENET);
}

void CPrime95Doc::OnServer() 
{
	pingServer ();
}

// Other internal commands

void CPrime95Doc::OnWelcome() 
{
	CWelcomeDlg dlg;

// Set global flag indicating startup is in progress.  This will delay
// starting any communication with the server until the user has confirmed
// he wants to use primenet and he has selected his work preferences.
	
	STARTUP_IN_PROGRESS = 1;

// After the welcome screen, install prime95 as an auto-start program
// and then go collect the user information.

	if (dlg.DoModal () == IDOK) {
		STRESS_TESTER = 0;
		IniWriteInt (INI_FILE, "StressTester", 0);
		USE_PRIMENET = 1;
		IniWriteInt (INI_FILE, "UsePrimenet", 1);
		OnPrimenet();
		if (USE_PRIMENET && STARTUP_IN_PROGRESS) OnCpu ();
		if (STARTUP_IN_PROGRESS) OnResources ();
		if (USE_PRIMENET && STARTUP_IN_PROGRESS) OnWorkerThreads ();
		if (USE_PRIMENET && STARTUP_IN_PROGRESS) {
			STARTUP_IN_PROGRESS = 0;
			set_comm_timers ();
			OnContinue ();
		} else
			STARTUP_IN_PROGRESS = 0;
	} else {
		STRESS_TESTER = 1;
		IniWriteInt (INI_FILE, "StressTester", 1);
		USE_PRIMENET = 0;
		IniWriteInt (INI_FILE, "UsePrimenet", 0);
		STARTUP_IN_PROGRESS = 0;
		OnTorture ();
	}
}

void CPrime95Doc::OnUsrTorture() 
{
	int	num_threads;

	num_threads = IniGetInt (INI_FILE, "TortureThreads", NUM_CPUS * CPU_HYPERTHREADS);
	LaunchTortureTest (num_threads, FALSE);
}

/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc private routines


void flashWindowAndBeep ()
{
	CWinApp* pApp = AfxGetApp();
	pApp->m_pMainWnd->FlashWindow (TRUE);
	MessageBeep (0xFFFFFFFF);
}


/////////////////////////////////////////////////////////////////////////////
// CPrime95Doc public routines

#include <ctype.h>
#include <dos.h>
#include <fcntl.h>
#include <io.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gwutil.h>

#ifdef X86_64
#define PORT	4
#else
#define PORT	1
#endif
#include "cJSON.h"
#include "cJSON.c"
#include "pm1prob.h"
#include "pm1prob.c"
#include "md5.c"
#include "comm95b.c"
#include "comm95c.c"
#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "primenet.c"
#include "proof_upload.c"
#include "proof_getdata.c"
#include "gwtest.c"

/* Do some work prior to launching worker threads */

void PreLaunchCallback (
	int	launch_type)
{

// Stall if we've just booted (within 5 minutes of Windows starting)

	if (GetTickCount () < 300000 && launch_type == LD_CONTINUE) {
		int	delay;
		delay = IniGetInt (INI_FILE, "BootDelay", 90);
		delay -= GetTickCount () / 1000;
		if (delay > 0) {
			char buf[50];
			sprintf (buf, "Waiting %d seconds for boot to complete.\n", delay);
			OutputStr (MAIN_THREAD_NUM, buf);
			Sleep (delay * 1000);
		}
	}
}

/* Do some work after worker threads have terminated */

void PostLaunchCallback (
	int	launch_type)
{
}

/* OSes that must poll for whether the ESC key was hit do it here. */
/* We use this opportunity to perform other miscellaneous tasks that */
/* can't be done any other way. */

void stopCheckCallback (
	int	thread_num)
{
	// If we couldn't add the icon to the task bar, then keep
	// trying until we finally succeed!
	if (WINDOWS95_TRAY_ADD) {
		CPrime95App *pApp = (CPrime95App *)AfxGetApp();
		pApp->TrayMessage (NIM_ADD, NULL, 0);
	}
}

// Output a status report for the range

void CPrime95Doc::rangeStatus ()
{
	char	buf[2000];

	rangeStatusMessage (buf, sizeof (buf));

	AfxMessageBox (buf, MB_ICONINFORMATION);
}
