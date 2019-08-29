// THIS CODE WAS LIBERALLY ADAPTED FROM A MICROSOFT SAMPLE NT SERVICE
//
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (C) 1993-1995  Microsoft Corporation.  All Rights Reserved.
//
//  MODULE:   service.c
//
//  PURPOSE:  Implements functions required by all services windows.
//
//  FUNCTIONS:
//    service_ctrl (DWORD dwCtrlCode);
//    service_main (DWORD dwArgc, LPTSTR *lpszArgv);
//
//  COMMENTS:
//
//  AUTHOR: Craig Link - Microsoft Developer Support
//


#include <windows.h>
#include <direct.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <tchar.h>
#include <afxres.h>
#include "resource.h"

// Variables used to communicate with the main prime95 code
char	NTSERVICENAME[32] = {0};	// name of the NT service
HWND	MAINFRAME_HWND = 0;		// Handle of main frame window
int	THREAD_ACTIVE = 0;		// TRUE if worker thread is active


// internal variables
SERVICE_STATUS		ssStatus;       // current status of the service
SERVICE_STATUS_HANDLE	sshStatusHandle;

// internal function prototypes
VOID WINAPI service_ctrl (DWORD dwCtrlCode);
VOID WINAPI service_main (DWORD dwArgc, LPTSTR *lpszArgv);
BOOL ReportStatusToSCMgr (DWORD dwCurrentState, DWORD dwWin32ExitCode, DWORD dwWaitHint);
void AddToMessageLog (LPTSTR lpszMsg);

// C runtime library's normal entry point
int WinMainCRTStartup (void);


// This is prime95's entry point.  It is called even before the C runtime
// library and MFC are initialized.  Our goal is to first determine if we are
// running as an NT service.  If so, we do some required NT service
// initialization (StartServiceCtrlDispatcher).  Then we fire up the
// regular C runtime library initialization code.

int MyWinMainCRTStartup (void)
{
	SERVICE_TABLE_ENTRY dispatchTable[] = {
		{ NTSERVICENAME, (LPSERVICE_MAIN_FUNCTION) service_main },
		{ NULL, NULL }
	};
	STARTUPINFO sui;

/* Is prime95 being run as an NT service or a normal program?  We determine */
/* this by looking at the USESHOWWINDOW flag.  If set this is not an NT */
/* service run.  Also, when started from an MS-DOS command prompt, the */
/* flags are zero.  Just jump to the normal C runtime library startup code. */

	GetStartupInfo (&sui);
	if (sui.dwFlags == 0 || (sui.dwFlags & STARTF_USESHOWWINDOW)) {
		return (WinMainCRTStartup ());
	}

/* Start the service thread */

	if (!StartServiceCtrlDispatcher (dispatchTable))
		AddToMessageLog ("StartServiceCtrlDispatcher failed.");

	return (0);
}


//  FUNCTION: service_main
//
//  PURPOSE: To perform actual initialization of the service
//
//  PARAMETERS:
//    dwArgc   - number of command line arguments
//    lpszArgv - array of command line arguments
//
//  RETURN VALUE:
//    none
//
//  COMMENTS:
//    This routine performs the service initialization and then calls
//    the user defined ServiceStart() routine to perform majority
//    of the work.
//
void WINAPI service_main (DWORD dwArgc, LPTSTR *lpszArgv)
{

// Copy the service name from the first lpszArgv
// Prime95 will use this to determine the named_ini_files value

	strcpy (NTSERVICENAME, lpszArgv[0]);

// register our service control handler:

	sshStatusHandle = RegisterServiceCtrlHandler (NTSERVICENAME, service_ctrl);
	if (!sshStatusHandle) goto cleanup;

// SERVICE_STATUS members that don't change in example

	ssStatus.dwServiceType = SERVICE_INTERACTIVE_PROCESS | SERVICE_WIN32_OWN_PROCESS;
	ssStatus.dwServiceSpecificExitCode = 0;

// report the status to the service control manager.

	if (! ReportStatusToSCMgr (
			SERVICE_RUNNING,	// service state
			NO_ERROR,		// code
			0))			// wait hint
		goto cleanup;

// Start the GUI service by calling the normal C runtime startup code

	WinMainCRTStartup ();
	return;

// try to report the stopped status to the service control manager.

cleanup:
	if (sshStatusHandle)
		(VOID) ReportStatusToSCMgr (SERVICE_STOPPED, 0, 0);
	return;
}


//  FUNCTION: service_ctrl
//
//  PURPOSE: This function is called by the SCM whenever
//           ControlService() is called on this service.
//
//  PARAMETERS:
//    dwCtrlCode - type of control requested
//
//  RETURN VALUE:
//    none
//
//  COMMENTS:
//
VOID WINAPI service_ctrl (DWORD dwCtrlCode)
{

// Handle the requested control code.

	switch (dwCtrlCode) {

// Stop the service.
// Note: Only the system can send SERVICE_CONTROL_SHUTDOWN to
// a service, otherwise the function ControlService will fail.
// USE WITH CAUTION:
// The SERVICE_CONTROL_SHUTDOWN control should only be processed by services
// that must absolutely clean up during shutdown, because there is an extremely
// limited time (about 20 seconds) available for service shutdown. After this
// time expires, system shutdown proceeds regardless of whether service
// shutdown is complete. If the service needs to take more time to shut down,
// it should send out STOP_PENDING status messages, along with a wait hint, so
// that the service controller knows how long to wait before reporting to the
// system that service shutdown is complete. For example, the server service
// needs to shut down so that network connections are not made when the system
// is in the shutdown state.

	case SERVICE_CONTROL_SHUTDOWN:
	case SERVICE_CONTROL_STOP:
		ssStatus.dwCurrentState = SERVICE_STOP_PENDING;
		ReportStatusToSCMgr (ssStatus.dwCurrentState, NO_ERROR, 0);
		if (THREAD_ACTIVE) {
			SendMessage (MAINFRAME_HWND, WM_COMMAND, IDM_STOP, 0);
			while (THREAD_ACTIVE) Sleep (50);
		}
		SendMessage (MAINFRAME_HWND, USR_SERVICE_STOP, 0, 0);
		ssStatus.dwCurrentState = SERVICE_STOPPED;
		break;

// Update the service status.

	case SERVICE_CONTROL_INTERROGATE:
		break;

// invalid control code

	default:
		break;
	}

	ReportStatusToSCMgr (ssStatus.dwCurrentState, NO_ERROR, 0);
}


//  FUNCTION: ReportStatusToSCMgr()
//
//  PURPOSE: Sets the current status of the service and
//           reports it to the Service Control Manager
//
//  PARAMETERS:
//    dwCurrentState - the state of the service
//    dwWin32ExitCode - error code to report
//    dwWaitHint - worst case estimate to next checkpoint
//
//  RETURN VALUE:
//    TRUE  - success
//    FALSE - failure
//
//  COMMENTS:
//
BOOL ReportStatusToSCMgr (
	DWORD	dwCurrentState,
	DWORD	dwWin32ExitCode,
	DWORD	dwWaitHint)
{
static	DWORD	dwCheckPoint = 1;
	BOOL	fResult = TRUE;

	if (dwCurrentState == SERVICE_START_PENDING)
		ssStatus.dwControlsAccepted = 0;
	else
		ssStatus.dwControlsAccepted = SERVICE_ACCEPT_STOP |
					      SERVICE_ACCEPT_SHUTDOWN;

	ssStatus.dwCurrentState = dwCurrentState;
	ssStatus.dwWin32ExitCode = dwWin32ExitCode;
	ssStatus.dwWaitHint = dwWaitHint;

	if (dwCurrentState == SERVICE_RUNNING ||
	    dwCurrentState == SERVICE_STOPPED)
		ssStatus.dwCheckPoint = 0;
	else
		ssStatus.dwCheckPoint = dwCheckPoint++;

// Report the status of the service to the service control manager.

	if (!(fResult = SetServiceStatus (sshStatusHandle, &ssStatus))) {
		AddToMessageLog ("SetServiceStatus");
	}

	return fResult;
}



//  FUNCTION: AddToMessageLog (LPTSTR lpszMsg)
//
//  PURPOSE: Allows any thread to log an error message
//
//  PARAMETERS:
//    lpszMsg - text for message
//
//  RETURN VALUE:
//    none
//
//  COMMENTS:
//
VOID AddToMessageLog (LPTSTR lpszMsg)
{
	TCHAR	szMsg[256];
	HANDLE	hEventSource;
	LPTSTR	lpszStrings[2];
	DWORD	dwErr;

	dwErr = GetLastError();

// Use event logging to log the error.

	hEventSource = RegisterEventSource (NULL, NTSERVICENAME);

	_stprintf (szMsg, "%s error: %d", NTSERVICENAME, dwErr);
	lpszStrings[0] = szMsg;
	lpszStrings[1] = lpszMsg;

	if (hEventSource != NULL) {
		ReportEvent (
			hEventSource,		// handle of event source
			EVENTLOG_ERROR_TYPE,	// event type
			0,			// event category
			0,			// event ID
			NULL,			// current user's SID
			2,			// strings in lpszStrings
			0,			// no bytes of raw data
			lpszStrings,		// array of error strings
			NULL);			// no raw data
		(VOID) DeregisterEventSource (hEventSource);
	}
}
