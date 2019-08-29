// MainFrm.cpp : implementation of the CMainFrame class
//

#include "stdafx.h"
#include "Prime95.h"

#include "MainFrm.h"
#include "Prime95Doc.h"
#include "Prime95View.h"

#include <winreg.h>
#include <pbt.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

const UINT WM_TASKBARCREATED = ::RegisterWindowMessage(_T("TaskbarCreated"));

/////////////////////////////////////////////////////////////////////////////
// CMainFrame

IMPLEMENT_DYNAMIC(CMainFrame, CMDIFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CMDIFrameWnd)
	//{{AFX_MSG_MAP(CMainFrame)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_ENDSESSION()
	ON_WM_ACTIVATEAPP()
	//}}AFX_MSG_MAP
	ON_WM_SYSCOMMAND()
	// Global help commands
	ON_COMMAND(ID_HELP_FINDER, CMDIFrameWnd::OnHelpFinder)
	ON_COMMAND(ID_HELP, CMDIFrameWnd::OnHelp)
	ON_COMMAND(ID_CONTEXT_HELP, CMDIFrameWnd::OnContextHelp)
	ON_COMMAND(ID_DEFAULT_HELP, CMDIFrameWnd::OnHelpFinder)
	ON_COMMAND(IDM_TRAY_OPEN, OnTrayOpenWindow)
	ON_COMMAND(IDM_STOP_CONTINUE, OnStopContinue)
	ON_COMMAND(ID_WINDOW_TILE_HORZ, OnTile)
	ON_COMMAND(ID_WINDOW_POSITION, OnPosition)
	ON_MESSAGE(USR_SERVICE_STOP, OnServiceStop)
	ON_MESSAGE(WM_POWERBROADCAST, OnPower)
	ON_MESSAGE(MYWM_TRAYMESSAGE, OnTrayMessage)
	ON_REGISTERED_MESSAGE(WM_TASKBARCREATED, OnTaskBarCreated)
END_MESSAGE_MAP()

static UINT indicators[] =
{
	ID_SEPARATOR,           // status line indicator
	ID_INDICATOR_CAPS,
	ID_INDICATOR_NUM,
	ID_INDICATOR_SCRL,
};

/////////////////////////////////////////////////////////////////////////////
// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
	// TODO: add member initialization code here
}

CMainFrame::~CMainFrame()
{
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

//	cs.style &= ~FWS_ADDTOTITLE;	// We'll control window titles!
	return CMDIFrameWnd::PreCreateWindow(cs);
}

void CMainFrame::OnUpdateFrameTitle(BOOL bAddToTitle)
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CMDIFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;

	if (!m_wndStatusBar.Create(this) ||
		!m_wndStatusBar.SetIndicators(indicators,
		  sizeof(indicators)/sizeof(UINT)))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}

	MAINFRAME_HWND = m_hWnd;

	return 0;
}

/////////////////////////////////////////////////////////////////////////////
// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CMDIFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CMDIFrameWnd::Dump(dc);
}

#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CMainFrame message handlers
LRESULT CMainFrame::OnPower(WPARAM uID, LPARAM uMouseMsg)
{
// We used to check battery status only after we got this message.  However,
// on my laptop prime95 sometimes would not restart when AC power resumed.
// Now we'll test using a timer as well as this Windows message (so we don't
// have a delay in shutting threads down when we go on battery)
	if (uID == PBT_APMPOWERSTATUSCHANGE) {
		if (is_timed_event_active (TE_BATTERY_CHECK)) test_battery ();
	}
	return 0;
}

LRESULT CMainFrame::OnTrayMessage(WPARAM uID, LPARAM uMouseMsg)
{
	if (uID == 352 && uMouseMsg == WM_LBUTTONDBLCLK) {
		if (IsWindowVisible())
			ShowWindow(FALSE);	// hide it
		else {
			ShowWindow(TRUE);	// show it
			SetForegroundWindow();
		}
	}

	if (uID == 352 && uMouseMsg == WM_RBUTTONUP) {
		POINT	point;
		HMENU	hMenu, hSubMenu;
		MENUITEMINFO m;

		GetCursorPos (&point);		// Get cursor location
		hMenu = LoadMenu (AfxGetApp()->m_hInstance,
				  MAKEINTRESOURCE (IDR_TRAYMENU));
		hSubMenu = GetSubMenu (hMenu, 0);
		m.cbSize = sizeof (m);
		if (!IsIconic ()) {
			m.fMask = MIIM_STRING;
			m.fType = MFT_STRING;
			m.dwTypeData = "Minimize Window";
			m.cch = (UINT) strlen (m.dwTypeData);
			SetMenuItemInfo (hSubMenu, IDM_TRAY_OPEN, FALSE, &m);
		}
		if (!WORKER_THREADS_ACTIVE) {
			m.fMask = MIIM_STRING;
			m.fType = MFT_STRING;
			m.dwTypeData = "Continue";
			m.cch = (UINT) strlen (m.dwTypeData);
			SetMenuItemInfo (hSubMenu, IDM_STOP_CONTINUE, FALSE, &m);
		}
		if ((WORKER_THREADS_ACTIVE && WORKER_THREADS_STOPPING) ||
		    (!WORKER_THREADS_ACTIVE &&
		     !USE_PRIMENET && WORKTODO_COUNT == 0)) {
			m.fMask = MIIM_STATE;
			m.fState = MFS_DISABLED;
			SetMenuItemInfo (hSubMenu, IDM_STOP_CONTINUE, FALSE, &m);
		}
		SetMenuDefaultItem (hSubMenu, IDM_TRAY_OPEN, FALSE);
		SetForegroundWindow ();	// Per KB Article Q135788
		TrackPopupMenu (hSubMenu,
			TPM_LEFTBUTTON | TPM_RIGHTBUTTON | TPM_LEFTALIGN,
			point.x, point.y, 0, m_hWnd, NULL);
		PostMessage (WM_NULL, 0, 0);// Per KB Article Q135788
		DestroyMenu (hMenu);
	}

	return 0;
}

LRESULT CMainFrame::OnTaskBarCreated (WPARAM wp, LPARAM lp)
{
	if (!HIDE_ICON)
		((CPrime95App *)AfxGetApp())->TrayMessage (NIM_ADD, "Prime95", 0);
	return 0;
}

void CMainFrame::OnTrayOpenWindow()
{
	if (IsIconic ())
		PostMessage (WM_SYSCOMMAND, SC_RESTORE, 0L);
	else
		PostMessage (WM_SYSCOMMAND, SC_MINIMIZE, 0L);
}

void CMainFrame::OnStopContinue()
{
	if (WORKER_THREADS_ACTIVE)
		PostMessage (WM_COMMAND, IDM_STOP, 0);
	else
		PostMessage (WM_COMMAND, IDM_CONTINUE, 0);
}

LRESULT CMainFrame::OnServiceStop (WPARAM wParam, LPARAM lParam)
{
	if (TRAY_ICON) ((CPrime95App *)AfxGetApp())->TrayMessage (NIM_DELETE, NULL, 0);
	ShowWindow (FALSE);
	CloseHandle (g_hMutexInst);
	return (0);
}

BOOL CMainFrame::DestroyWindow() 
{
	if (TRAY_ICON) ((CPrime95App *)AfxGetApp())->TrayMessage (NIM_DELETE, NULL, 0);
	return CMDIFrameWnd::DestroyWindow();
}

void CMainFrame::OnSize(UINT nType, int cx, int cy) 
{
static	int first_sizing = 0;

	if (nType == SIZE_MINIMIZED) {
		if (TRAY_ICON) ShowWindow (FALSE);// hide it
		if (HIDE_ICON) ShowWindow (FALSE);
	}
	CMDIFrameWnd::OnSize(nType, cx, cy);

// On Windows XP, we need to tile on the second message, while on 64-bit
// windows tiling on the first size message works.  Go figure.

	if (nType != SIZE_MINIMIZED && first_sizing < 2) {
		int id;
		char rgch[80];

		rgch[0] = 0;
		IniGetString (INI_FILE, "W0", rgch, sizeof(rgch), NULL);

		first_sizing++;
		id = (0 == rgch[0]) ? ID_WINDOW_TILE_HORZ : ID_WINDOW_POSITION;
		PostMessage (WM_COMMAND, id, 0);
	}
}

LRESULT CMainFrame::WindowProc (UINT message, WPARAM wParam, LPARAM lParam)
{
	WM_ENDSESSION_LPARAM = (LONG) lParam;
	return CMDIFrameWnd::WindowProc(message, wParam, lParam);
}


// Special WM_ENDSESSION code from knowledge base article Q164166

static const TCHAR szAfxOldWndProc[] = _T("AfxOldWndProc423");
BOOL CALLBACK EnumProc (HWND hWnd, LPARAM lParam)
{
	//check for property and unsubclass if necessary
	WNDPROC oldWndProc = (WNDPROC)::GetProp (hWnd, szAfxOldWndProc);
	if (oldWndProc!=NULL) {
		SetWindowLongPtr (hWnd, GWLP_WNDPROC, (LONG_PTR) oldWndProc);
		RemoveProp (hWnd, szAfxOldWndProc);
	}
	return TRUE;
}

void CMainFrame::OnEndSession (BOOL bEnding)
{

// If we are running as a service, then do not end the program at logoff
// Instead, just remove the icon from the system tray.

	if (NTSERVICENAME[0] ||
	    (WINDOWS95_SERVICE && WM_ENDSESSION_LPARAM && isWindows95 ())) {
		if (TRAY_ICON) {
			((CPrime95App *)AfxGetApp())->TrayMessage (NIM_DELETE, NULL, 0);
			WINDOWS95_TRAY_ADD = 1;
		}

// In addition a Windows NT service must take special actions.  MFC was
// not designed to be used in a NT service as it uses Global Atoms which are
// cleared at logoff.  This fix from knowledge base article Q164166 seems
// to fix the problem.

		if (NTSERVICENAME[0]) {
			DWORD	dwProcessId, dwThreadId;
			dwThreadId = GetWindowThreadProcessId (m_hWnd, &dwProcessId);
			EnumThreadWindows (dwThreadId, EnumProc, (LPARAM) dwThreadId);
		}
	}

// If we aren't running as a service, just do normal processing

	else
		CMDIFrameWnd::OnEndSession (bEnding);
}

/* Return true if user has permission to create and delete services */

#include <winsvc.h>
int canModifyServices ()
{
	SC_HANDLE schSCManager;

/* All Windows 9x users have permission */

	if (isWindows95 ()) return (TRUE);

/* Vista won't let services interact with the desktop */

	if (isWindowsVista ()) return (FALSE);

/* See if Windows NT user can open service control manager */

	schSCManager = OpenSCManager (
		NULL,		// machine (NULL == local)
		NULL,		// database (NULL == default)
		SC_MANAGER_ALL_ACCESS);	// access required
	if (! schSCManager) return (FALSE);
	CloseServiceHandle (schSCManager);
	return (TRUE);
}

/* Handle all the details for running prime95 as a service */

#define RSP_SIMPLE_SERVICE	1
#define RSP_UNREGISTER_SERVICE	0
void Service95 ()
{
	char	pathname[256];
	char	regkey[20];		/* Win9x registry name */
	SC_HANDLE schSCManager = 0;
	SC_HANDLE schService = 0;
	HKEY	hkey = 0;
	DWORD	rc, disposition;

/* In Windows 95/98/Me, call RegisterServiceProcess in the Kernel */
/* This will prevent prime95 from terminating on logoff. */

	if (isWindows95 ()) {
		HMODULE	hlib;
		DWORD (__stdcall *proc)(DWORD, DWORD);

		hlib = LoadLibrary ("KERNEL32.DLL");
		if (!hlib) {
			OutputStr (MAIN_THREAD_NUM, "Unable to load KERNEL32.DLL\n");
			goto done;
		}
		proc = (DWORD (__stdcall *)(DWORD, DWORD))
			GetProcAddress (hlib, "RegisterServiceProcess");
		if (proc == NULL)
			OutputStr (MAIN_THREAD_NUM, "Unable to find RegisterServiceProcess\n");
		else {
			if (WINDOWS95_SERVICE)
				rc = (*proc) (NULL, RSP_SIMPLE_SERVICE);
			else
				rc = (*proc) (NULL, RSP_UNREGISTER_SERVICE);
			if (!rc)
				OutputStr (MAIN_THREAD_NUM, "RegisterServiceProcess failed\n");
		}
		FreeLibrary (hlib);
	}

/* Now we deal with making the registry entries correct for proper starting */
/* or not starting of the service. */

/* Get pathname of executable */

	GetModuleFileName (NULL, pathname, sizeof (pathname));

/* In Win95/98/Me, we create a registry entry for each -A command line value */
/* We used to do this in WinNT/2000/XP, but now we just delete these old */
/* registry entries. */

	if (WINDOWS95_A_SWITCH < 0)
		strcpy (regkey, "Prime95");
	else
		sprintf (regkey, "Prime95-%d", WINDOWS95_A_SWITCH);

// For Windows 95/98/Me create a RunServices entry

	if (isWindows95 ()) {
		if (RegCreateKeyEx (
				HKEY_LOCAL_MACHINE,
				"Software\\Microsoft\\Windows\\CurrentVersion\\RunServices",
				0,
				NULL,
				REG_OPTION_NON_VOLATILE,
				KEY_ALL_ACCESS,
				NULL,
				&hkey,
				&disposition) != ERROR_SUCCESS) {
			OutputStr (MAIN_THREAD_NUM, "Can't create registry key.\n");
			goto done;
		}

/* Now create or delete an entry for prime95 */

		if (WINDOWS95_SERVICE) {
			if (WINDOWS95_A_SWITCH >= 0) {
				char	append[20];
				sprintf (append, " -A%d", WINDOWS95_A_SWITCH);
				strcat (pathname, append);
			}
			rc = RegSetValueEx (hkey, regkey, 0, REG_SZ,
				(BYTE *) pathname, (DWORD) strlen (pathname) + 1);
			if (rc != ERROR_SUCCESS) {
				OutputStr (MAIN_THREAD_NUM, "Can't write registry value.\n");
				goto done;
			}
		} else {
			rc = RegDeleteValue (hkey, regkey);
			if (rc != ERROR_SUCCESS && rc != ERROR_FILE_NOT_FOUND){
				OutputStr (MAIN_THREAD_NUM, "Can't delete registry entry.\n");
				goto done;
			}
		}
		RegCloseKey (hkey);
		hkey = 0;
	}

// For Windows NT/2000/XP we call the service control manager to maintain the
// services database.  If we don't have administrator privileges (can't open
// the service control manager) then simply create a registry entry to start
// the program at logon.  So, attempt to open the service control manager on
// NULL = local machine, NULL = default database, all access required.
// Also Vista won't let services interact with the desktop, so we can't run
// as a service.

	else if (isWindowsVista () ||
		 ! (schSCManager = OpenSCManager (NULL, NULL,
						  SC_MANAGER_ALL_ACCESS))) {
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
			OutputStr (MAIN_THREAD_NUM, "Can't create registry key.\n");
			goto done;
		}

/* Now create or delete an entry for prime95 */

		if (WINDOWS95_SERVICE) {
			if (WINDOWS95_A_SWITCH >= 0) {
				char	append[20];
				sprintf (append, " -A%d", WINDOWS95_A_SWITCH);
				strcat (pathname, append);
			}
			rc = RegSetValueEx (hkey, regkey, 0, REG_SZ,
				(BYTE *) pathname, (DWORD) strlen (pathname) + 1);
			if (rc != ERROR_SUCCESS) {
				OutputStr (MAIN_THREAD_NUM, "Can't write registry value.\n");
				goto done;
			}
		} else {
			rc = RegDeleteValue (hkey, regkey);
			if (rc != ERROR_SUCCESS && rc != ERROR_FILE_NOT_FOUND){
				OutputStr (MAIN_THREAD_NUM, "Can't delete registry entry.\n");
				goto done;
			}
		}
		RegCloseKey (hkey);
		hkey = 0;
	}

// Make the necessary NT/2000/XP service control changes

	else {
		char	servicename[80];
		char	displayname[80];

// Create the Windows NT service name and display name

		IniGetString (LOCALINI_FILE, "ServiceName", servicename,
			      sizeof (servicename), NULL);
		if (servicename[0] == 0) {
			if (WINDOWS95_A_SWITCH < 0)
				strcpy (servicename, "Prime95 Service");
			else
				sprintf (servicename, "Prime95 Service-%d",
					 WINDOWS95_A_SWITCH);
		}

		IniGetString (LOCALINI_FILE, "DisplayName", displayname,
			      sizeof (displayname), servicename);

// Create the service entry

		if (WINDOWS95_SERVICE) {
			schService = CreateService (
				schSCManager,		// SCManager database
				servicename,		// name of service
				displayname,		// display name
				SERVICE_ALL_ACCESS,	// desired access
				SERVICE_INTERACTIVE_PROCESS |
				SERVICE_WIN32_OWN_PROCESS,  // service type
				SERVICE_AUTO_START,	// start type
				SERVICE_ERROR_NORMAL,	// error control type
				pathname,		// service's binary
				NULL,			// no load ordering
				NULL,			// no tag identifier
				NULL,			// no dependencies
				NULL,			// LocalSystem account
				NULL);			// no password
			if (!schService) {
				if (GetLastError () != ERROR_SERVICE_EXISTS)
					OutputStr (MAIN_THREAD_NUM, "Error creating service.\n");
				goto done;
			}

// Set description for Win2K and later

			if (isWindows2000 ()) {
				SERVICE_DESCRIPTION svc_desc;
				svc_desc.lpDescription =
					"GIMPS client to find large prime numbers";
				ChangeServiceConfig2 (
					schService,
					SERVICE_CONFIG_DESCRIPTION,
					&svc_desc);
			}
		}

// Remove the service entry

		else {
			schService = OpenService (
					schSCManager,
					servicename,
					SERVICE_ALL_ACCESS);
			if (!schService) {
				if (GetLastError () != ERROR_SERVICE_DOES_NOT_EXIST)
					OutputStr (MAIN_THREAD_NUM, "Error opening service.\n");
				goto done;
			}
			if (! DeleteService (schService)) {
				OutputStr (MAIN_THREAD_NUM, "Error deleting service.\n");
				goto done;
			}
		}

/* Delete the old-style (version 21) Run entry for the current user */

		if (RegCreateKeyEx (
				HKEY_CURRENT_USER,
				"Software\\Microsoft\\Windows\\CurrentVersion\\Run",
				0,
				NULL,
				REG_OPTION_NON_VOLATILE,
				KEY_ALL_ACCESS,
				NULL,
				&hkey,
				&disposition) == ERROR_SUCCESS) {
			RegDeleteValue (hkey, regkey);
			RegCloseKey (hkey);
			hkey = 0;
		}
	}

/* For Windows NT/2000/XP delete the old-style (version 21) Run entry for */
/* the local machine */

	if (! isWindows95 ()) {
		if (RegCreateKeyEx (
				HKEY_LOCAL_MACHINE,
				"Software\\Microsoft\\Windows\\CurrentVersion\\Run",
				0,
				NULL,
				REG_OPTION_NON_VOLATILE,
				KEY_ALL_ACCESS,
				NULL,
				&hkey,
				&disposition) == ERROR_SUCCESS) {
			RegDeleteValue (hkey, regkey);
			RegCloseKey (hkey);
			hkey = 0;
		}
	}

/* Now delete any shortcuts to the program.  We do this so that as users */
/* upgrade from version 20 and try this menu choice they do not end up with */
/* both a registry entry and a StartUp menu shortcut. */

	if (WINDOWS95_SERVICE) {
		char	buf[256];
		DWORD	type;
		DWORD	bufsize = sizeof (buf);

		if (RegCreateKeyEx (
				HKEY_CURRENT_USER,
				"Software\\Microsoft\\Windows\\CurrentVersion\\Explorer\\Shell Folders",
				0,
				NULL,
				REG_OPTION_NON_VOLATILE,
				KEY_ALL_ACCESS,
				NULL,
				&hkey,
				&disposition) != ERROR_SUCCESS)
			goto done;
		if (RegQueryValueEx (hkey, "Startup", NULL, &type,
				(BYTE *) buf, &bufsize) == ERROR_SUCCESS &&
		    type == REG_SZ) {
			strcat (buf, "\\prime95.lnk");
			_unlink (buf);
		}
	}

// Cleanup and return

done:	if (schService) CloseServiceHandle (schService);
	if (schSCManager) CloseServiceHandle (schSCManager);
	if (hkey) RegCloseKey (hkey);
}

void CMainFrame::OnActivateApp(BOOL bActive, DWORD hTask) 
{
	CMDIFrameWnd::OnActivateApp(bActive, hTask);
}

// Override the Upper Right X to do a minimize instead of a close.
// This has become common practice for tray applications.

void CMainFrame::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == SC_CLOSE &&
	    ! IniGetInt (INI_FILE, "ExitOnX", 0) &&
	    (TRAY_ICON || HIDE_ICON)) {
		SendMessage (WM_SYSCOMMAND, SC_MINIMIZE, 0);
		return;
	}
	CMDIFrameWnd::OnSysCommand (nID, lParam);
}

// Override the MFC implementation of tile MDI windows horizontally.

void CMainFrame::OnTile()
{

// Call default implementation to de-maximize windows if necessary

	MDITile (MDITILE_HORIZONTAL);

// Now call our implementation of horizontal tiling

	PositionViews (TRUE);
}

void CMainFrame::OnPosition()
{

// Call default implementation to de-maximize windows if necessary

	MDITile (MDITILE_HORIZONTAL);

// Now call our implementation of horizontal tiling

	PositionViews (FALSE);
}

