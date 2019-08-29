//
//  AppController.m
//  Prime95
//
//  Created by George Woltman on 4/17/09.
//  Copyright 2009-2019 Mersenne Research, Inc. All rights reserved.
//

#import "AppController.h"
#import "AboutController.h"
#import "BenchmarkController.h"
#import "ContinueController.h"
#import "CPUController.h"
#import "ECMController.h"
#import "ManualCommunicationController.h"
#import "Pminus1Controller.h"
#import "PreferencesController.h"
#import "PrimeNetController.h"
#import "PRPController.h"
#import "QuitGIMPSController.h"
#import "StatusController.h"
#import "StopController.h"
#import "TestController.h"
#import "TimeController.h"
#import "TortureTestController.h"
#import "UnreserveController.h"
#import "WorkerWindowsController.h"
#include "prime95.h"


/* Required Mac OS X files */
#ifdef __APPLE__
#include <dirent.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/sysctl.h>
#include <sys/time.h>
#include <sys/timeb.h>
#define PTHREAD_MIN_PRIORITY 0		/* Missing #defines from pthreads.h */
#define PTHREAD_MAX_PRIORITY 31		/* Missing #defines from pthreads.h */
#include <CoreFoundation/CoreFoundation.h>
#include <IOKit/ps/IOPowerSources.h>
#include <IOKit/ps/IOPSKeys.h>
#endif

/* Handle differences between Windows and Linux runtime libraries */

#define _commit(f)	fsync(f)
#define _open		open
#define _close		close
#define _read		read
#define _write		write
#define _lseek		lseek
#define _unlink		unlink
#define _creat		creat
#define _chdir		chdir
#define closesocket	close
#define IsCharAlphaNumeric(c) isalnum(c)
#define _stricmp	strcasecmp
#define _timeb		timeb
#define _ftime		ftime
#define _O_APPEND	O_APPEND
#define _O_RDONLY	O_RDONLY
#define _O_WRONLY	O_WRONLY
#define _O_RDWR		O_RDWR
#define _O_CREAT	O_CREAT
#define _O_TRUNC	O_TRUNC
#define _O_BINARY 	0
#define _O_TEXT		0

AppController *myAppController;			// Global variable to allow access to this object

@implementation AppController

- (void) awakeFromNib
{
	NSString *workingDir;

/* Init the global variable */

	myAppController = self;

/* Unlike other ports we do not change to the executable's directory.  This is because */
/* that would place our temporary files within the application bundle.  They would then */
/* be deleted whenever the user upgrades to a new version.  Instead, we'll make our */
/* default working directory ~/Prime95.  This default can be overridden using the plist */
/* editor (we write the default there so it is easy to find and overwrite by the user) */
/* or by specifying -WorkingDirectory "some_path" as a command line argument. */

/* SoB doesn't like this approach as it makes it hard for them to package up an executable */
/* with a default prime.txt to contact their server using the MersenneIP= option. */
/* Our workaround: if prime.txt exists in the application bundle we copy it to the new */
/* working directory if no prime.txt file exists there. */

	/* Change the working directory */
	workingDir = [[NSUserDefaults standardUserDefaults] stringForKey:@"WorkingDirectory"];
	if (workingDir == nil) {
		workingDir = @"~/Prime95";
		[[NSUserDefaults standardUserDefaults] setObject:workingDir forKey:@"WorkingDirectory"];
	}
	workingDir = [workingDir stringByExpandingTildeInPath];
	[[NSFileManager defaultManager] createDirectoryAtPath:workingDir withIntermediateDirectories: TRUE attributes:nil error:nil];
        [[NSFileManager defaultManager] changeCurrentDirectoryPath:workingDir];

	/* If no prime.txt file exists, optionally copy one from the application bundle (for SoB) */
	if (! fileExists ("prime.txt")) {
		char	filename[1025];
		int	infd;

		// First try getting the prime.txt from within the bundle
		sprintf (filename, "%s/prime.txt",
			 [[[[[NSBundle mainBundle] executablePath] stringByDeletingLastPathComponent] stringByDeletingLastPathComponent] UTF8String]);

		// Next try getting the prime.txt from the same directory as the app
		if (! fileExists (filename))
			sprintf (filename, "%s/prime.txt",
				 [[[[NSBundle mainBundle] bundlePath] stringByDeletingLastPathComponent] UTF8String]);

		// Open and copy the prime.txt file
		infd = _open (filename, _O_TEXT | _O_RDONLY);
		if (infd >= 0) {
			int	outfd;

			outfd = _open ("prime.txt", _O_TEXT | _O_WRONLY | _O_CREAT, CREATE_FILE_ACCESS);
			if (outfd >= 0) {
				char	*buf;
				int	buflen;

				buf = (char *) malloc (100000);
				if (buf != NULL) {
					buflen = _read (infd, buf, 100000);
					(void) _write (outfd, buf, buflen);
					free (buf);
				}
				_close (outfd);
			}
			_close (infd);
		}
	}

/* Initialize gwnum call back routines.  Using callback routines lets the */
/* gwnum library have a nice clean interface for users that do not need */
/* additional functionality that only prime95 uses. */

	StopCheckRoutine = stopCheck;
	OutputBothRoutine = OutputBoth;

/* Name and read the INI files.  Perform some other startup initializations. */

	nameAndReadIniFiles (-1);
	initCommCode ();

// On first run step user through the primenet dialog boxes.

	if (STRESS_TESTER == 99) {

// Set global flag indicating startup is in progress.  This will delay
// starting any communication with the server until the user has confirmed
// he wants to use primenet and he has selected his work preferences.

		STARTUP_IN_PROGRESS = 1;

// Clear flag saying this is the first run

		STRESS_TESTER = 0;
		IniWriteInt (INI_FILE, "StressTester", 0);

// Make using primenet the default

		USE_PRIMENET = 1;
		IniWriteInt (INI_FILE, "UsePrimenet", 1);

// Go collect the user information.
		
		[self testPrimeNet:nil];
	}

// Auto-continue if there is any work to do.

	else if (USE_PRIMENET || WELL_BEHAVED_WORK || WORKTODO_COUNT) {
		LaunchWorkerThreads (ALL_WORKERS, FALSE);
	}
}

// Intercept shutdown messages and save our state

- (NSApplicationTerminateReply)applicationShouldTerminate:(NSApplication *)app {

// Stop background threads before exiting

	if (WORKER_THREADS_ACTIVE) {
		stop_workers_for_escape ();
		while (WORKER_THREADS_STOPPING) Sleep (50);
	}

/* Write the worktodo file in case the WELL_BEHAVED_WORK flag caused us */
/* to delay writing the file. */

	writeWorkToDoFile (TRUE);

// Finish exiting

	return NSTerminateNow;
}

- (BOOL)validateMenuItem:(NSMenuItem *)item
{
	SEL theAction = [item action];

	if (theAction == @selector(showAboutPrimenet:)) {
		if (USE_PRIMENET) return YES;
		return NO;
	}

	if (theAction == @selector(testContinue:)) {
		[item setTitle:((!WORKER_THREADS_ACTIVE && NUM_WORKER_THREADS > 1) ||
				(WORKER_THREADS_ACTIVE && active_workers_count () != WORKER_THREADS_ACTIVE - 1) ?
				@"Continue..." : @"Continue")];
		if ((!WORKER_THREADS_ACTIVE && (USE_PRIMENET || WORKTODO_COUNT)) ||
		    (WORKER_THREADS_ACTIVE && active_workers_count () != WORKER_THREADS_ACTIVE)) return YES;
		return NO;
	}

	if (theAction == @selector(testStop:)) {
		[item setTitle:(active_workers_count () > 1 ? @"Stop..." : @"Stop")];
		if (WORKER_THREADS_ACTIVE && !WORKER_THREADS_STOPPING) return YES;
		return NO;
	}

	if (theAction == @selector(advancedManualCommunication:)) {
		if (USE_PRIMENET) return YES;
		return NO;
	}

	if (theAction == @selector(toggleSuminpErrorChecking:)) {
		if (SUM_INPUTS_ERRCHK) [item setState:NSOnState];
		else [item setState:NSOffState];
		return YES;
	}

	if (theAction == @selector(toggleErrorChecking:)) {
		if (ERRCHK) [item setState:NSOnState];
		else [item setState:NSOffState];
		return YES;
	}

	if (theAction == @selector(advancedUnreserve:)) {
		if (USE_PRIMENET) return YES;
		return NO;
	}

	if (theAction == @selector(advancedQuitGIMPS:)) {
		if (USE_PRIMENET || WORKTODO_COUNT) return YES;
		return NO;
	}

	if (theAction == @selector(toggleMergeMainComm:)) {
		if (MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS) [item setState:NSOnState];
		else [item setState:NSOffState];
		return YES;
	}

	if (theAction == @selector(toggleMergeMainCommWorker:)) {
		[item setTitle:
			(NUM_WORKER_THREADS == 1 ?
				@"Merge Main & Comm & Worker" :
			 MERGE_WINDOWS & MERGE_WORKER_WINDOWS ?
				@"Merge Main & Comm & Workers" :
				@"Merge Main & Comm & 1st Worker")];
		if (MERGE_WINDOWS & MERGE_MAIN_WINDOW && MERGE_WINDOWS & MERGE_COMM_WINDOW) [item setState:NSOnState];
		else [item setState:NSOffState];
		return YES;
	}

	if (theAction == @selector(toggleMergeAllWorkers:)) {
		if (MERGE_WINDOWS & MERGE_WORKER_WINDOWS) [item setState:NSOnState];
		else [item setState:NSOffState];
		if (NUM_WORKER_THREADS > 1) return YES;
		return NO;
	}

	if (theAction == @selector(advancedTime:)) {
		if (!WORKER_THREADS_STOPPING) return YES;
		return NO;
	}

	if (theAction == @selector(optionsTortureTest:)) {
		if (!WORKER_THREADS_STOPPING) return YES;
		return NO;
	}

	if (theAction == @selector(optionsBenchmark:)) {
		if (!WORKER_THREADS_STOPPING) return YES;
		return NO;
	}

	return YES;
}

// This menu choice displays the about box

- (IBAction)showAboutPanel:(id)sender
{
	if (!aboutController) aboutController = [[AboutController alloc] init];
	else [aboutController reInit];
	[aboutController showWindow:self];
}

// This menu choice displays the about primenet box

- (IBAction)showAboutPrimenet:(id)sender
{
	pingServer ();
}

// This menu choice displays the PrimeNet dialog

- (IBAction)testPrimeNet:(id)sender
{
	if (!primeNetController) primeNetController = [[PrimeNetController alloc] init];
	else [primeNetController reInit];
	[primeNetController showWindow:self];
}

// This menu choice displays the Worker Windows dialog

- (IBAction)testWorkerWindows:(id)sender
{
	if (!workerWindowsController) workerWindowsController = [[WorkerWindowsController alloc] init];
	else [workerWindowsController reInit];
	[workerWindowsController showWindow:self];
}

// This menu choice displays the work queue status

- (IBAction)testStatus:(id)sender
{
	if (!statusController) statusController = [[StatusController alloc] init];
	else [statusController reInit];
	[statusController showWindow:self];
}

// This menu choice quits the application

- (IBAction)shutDown:(id)sender
{

// Stop background threads before exiting

	if (WORKER_THREADS_ACTIVE) {
		stop_workers_for_escape ();
		while (WORKER_THREADS_STOPPING) Sleep (50);
	}

// Do the normal termination

	[NSApp terminate:self];
}

// This menu choice starts one or all workers

- (IBAction)testContinue:(id)sender
{
	if (NUM_WORKER_THREADS > 1 &&
	    active_workers_count () != WORKER_THREADS_ACTIVE - 1) {
		// Start the dialog box
		if (!continueController) continueController = [[ContinueController alloc] init];
		else [continueController reInit];
		[continueController showWindow:self];
	} else {
		// Start the thread
		LaunchWorkerThreads (ALL_WORKERS, FALSE);
	}
}

// This menu choice stops one or all workers

- (IBAction)testStop:(id)sender
{
	if (NUM_WORKER_THREADS > 1 && active_workers_count () != 1) {
		if (!stopController) stopController = [[StopController alloc] init];
		else [stopController reInit];
		[stopController showWindow:self];
	} else {
		// Stop the thread
		stop_workers_for_escape ();
	}
}

// This menu choice makes font bigger

- (IBAction)editBigger:(id)sender
{
	BiggerFonts ();
}

// This menu choice makes font smaller

- (IBAction)editSmaller:(id)sender
{
	SmallerFonts ();
}

// This menu choice LL tests one exponent

- (IBAction)advancedTest:(id)sender
{
	// Start the dialog box
	if (!testController) testController = [[TestController alloc] init];
	else [testController reInit];
	[testController showWindow:self];
}

// This menu choice times one exponent

- (IBAction)advancedTime:(id)sender
{
	// Start the dialog box
	if (!timeController) timeController = [[TimeController alloc] init];
	else [timeController reInit];
	[timeController showWindow:self];
}

// This menu choice does P-1 on a number

- (IBAction)advancedPminus1:(id)sender
{
	// Start the dialog box
	if (!pminus1Controller) pminus1Controller = [[Pminus1Controller alloc] init];
	else [pminus1Controller reInit];
	[pminus1Controller showWindow:self];
}

// This menu choice does ECM on a number

- (IBAction)advancedECM:(id)sender
{
	// Start the dialog box
	if (!ecmController) ecmController = [[ECMController alloc] init];
	else [ecmController reInit];
	[ecmController showWindow:self];
}

// This menu choice does PRP on a number

- (IBAction)advancedPRP:(id)sender
{
	// Start the dialog box
	if (!prpController) prpController = [[PRPController alloc] init];
	else [prpController reInit];
	[prpController showWindow:self];
}

// This menu choice manually communicates with the server

- (IBAction)advancedManualCommunication:(id)sender
{
	// Start the dialog box
	if (!manualCommunicationController) manualCommunicationController = [[ManualCommunicationController alloc] init];
	else [manualCommunicationController reInit];
	[manualCommunicationController showWindow:self];
}

// This menu choice toggles the SUM(INPUTS) != SUM(OUTPUTS) error checking

- (IBAction)toggleSuminpErrorChecking:(id)sender
{
	SUM_INPUTS_ERRCHK = !SUM_INPUTS_ERRCHK;
	IniWriteInt (INI_FILE, "SumInputsErrorCheck", SUM_INPUTS_ERRCHK);
}

// This menu choice toggles the roundoff error checking

- (IBAction)toggleErrorChecking:(id)sender
{
	ERRCHK = !ERRCHK;
	IniWriteInt (INI_FILE, "ErrorCheck", ERRCHK);
}

// This menu choice allows quitting GIMPS

- (IBAction)advancedQuitGIMPS:(id)sender
{
	// Start the dialog box
	if (!quitGIMPSController) quitGIMPSController = [[QuitGIMPSController alloc] init];
	else [quitGIMPSController reInit];
	[quitGIMPSController showWindow:self];
}

// This menu choice allows unreserving a number

- (IBAction)advancedUnreserve:(id)sender
{
	// Start the dialog box
	if (!unreserveController) unreserveController = [[UnreserveController alloc] init];
	else [unreserveController reInit];
	[unreserveController showWindow:self];
}

// This menu choice displays CPU info and settings

- (IBAction)optionsCPU:(id)sender
{
	// Start the dialog box
	if (!cpuController) cpuController = [[CPUController alloc] init];
	else [cpuController reInit];
	[cpuController showWindow:self];
}

// This menu choice displays prefernces panel

- (IBAction)optionsPreferences:(id)sender
{
	// Start the dialog box
	if (!preferencesController) preferencesController = [[PreferencesController alloc] init];
	else [preferencesController reInit];
	[preferencesController showWindow:self];
}

// This menu choice runs the torture test

- (IBAction)optionsTortureTest:(id)sender
{
	// Start the dialog box
	if (!tortureTestController) tortureTestController = [[TortureTestController alloc] init];
	else [tortureTestController reInit];
	[tortureTestController showWindow:self];
}

// This menu choice starts a benchmark

- (IBAction)optionsBenchmark:(id)sender
{
	// Start the dialog box
	if (!benchmarkController) benchmarkController = [[BenchmarkController alloc] init];
	else [benchmarkController reInit];
	[benchmarkController showWindow:self];
}

// This menu choice combines and uncombines the main and comm windows.

- (IBAction)toggleMergeMainComm:(id)sender
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
}

// This menu choice, in checked state, combines the main and comm and 1st worker windows.
// In going to unchecked state we leave the main and comm windows combined.

- (IBAction)toggleMergeMainCommWorker:(id)sender
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
}

// This menu choice, in checked state, combines the all worker windows into one window.

- (IBAction)toggleMergeAllWorkers:(id)sender
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
}

- (IBAction)helpMersenneForum:(id)sender
{
	NSURL *url = [NSURL URLWithString:@"http://mersenneforum.org"];
	[[NSWorkspace sharedWorkspace] openURL:url];
}

- (IBAction)helpMersenneWiki:(id)sender
{
	NSURL *url = [NSURL URLWithString:@"https://rieselprime.de/ziki/Main_Page"];
	[[NSWorkspace sharedWorkspace] openURL:url];
}

@end


/* Include our common C source files */

#include "commona.c"
#include "commonb.c"
#include "commonc.c"
#include "ecm.c"
#include "primenet.c"
#include "gwtest.c"

/* Implement the OS-specific routines */

#include "os_routines.c"

