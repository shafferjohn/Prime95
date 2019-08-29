//
//  WorkerWindowsController.m
//  Prime95
//
//  Created by George Woltman on 4/26/09.
//  Copyright 2009-2017 Mersenne Research, Inc. All rights reserved.
//

#import "WorkerWindowsController.h"
#import "AppController.h"
#include "prime95.h"

int map_work_pref_to_sel (
	int	work_pref)
{
	switch (work_pref) {
	case PRIMENET_WP_WHATEVER:
		return (0);
	case PRIMENET_WP_LL_FIRST:
		return (1);
	case PRIMENET_WP_LL_WORLD_RECORD:
		return (2);
	case PRIMENET_WP_LL_DBLCHK:
		return (3);
	case PRIMENET_WP_PRP_FIRST:
		return (4);
	case PRIMENET_WP_PRP_WORLD_RECORD:
		return (5);
	case PRIMENET_WP_PRP_DBLCHK:
		return (6);
	case PRIMENET_WP_FACTOR:
		return (7);
	case PRIMENET_WP_PFACTOR:
		return (8);
	case PRIMENET_WP_PRP_100M:
		return (9);
	case PRIMENET_WP_LL_100M:
		return (10);
	case PRIMENET_WP_PRP_COFACTOR:
		return (11);
	case PRIMENET_WP_PRP_COFACTOR_DBLCHK:
		return (12);
	case PRIMENET_WP_ECM_SMALL:
		return (13);
	case PRIMENET_WP_ECM_COFACTOR:
		return (14);
	case PRIMENET_WP_ECM_FERMAT:
		return (15);
	case PRIMENET_WP_FACTOR_LMH:
		return (16);
	default:
		return (17);
	}
}

int map_sel_to_work_pref (
	int	sel)
{
	switch (sel) {
	case 0:
		return (PRIMENET_WP_WHATEVER);
	case 1:
		return (PRIMENET_WP_LL_FIRST);
	case 2:
		return (PRIMENET_WP_LL_WORLD_RECORD);
	case 3:
		return (PRIMENET_WP_LL_DBLCHK);
	case 4:
		return (PRIMENET_WP_PRP_FIRST);
	case 5:
		return (PRIMENET_WP_PRP_WORLD_RECORD);
	case 6:
		return (PRIMENET_WP_PRP_DBLCHK);
	case 7:
		return (PRIMENET_WP_FACTOR);
	case 8:
		return (PRIMENET_WP_PFACTOR);
	case 9:
		return (PRIMENET_WP_PRP_100M);
	case 10:
		return (PRIMENET_WP_LL_100M);
	case 11:
		return (PRIMENET_WP_PRP_COFACTOR);
	case 12:
		return (PRIMENET_WP_PRP_COFACTOR_DBLCHK);
	case 13:
		return (PRIMENET_WP_ECM_SMALL);
	case 14:
		return (PRIMENET_WP_ECM_COFACTOR);
	case 15:
		return (PRIMENET_WP_ECM_FERMAT);
	case 16:
		return (PRIMENET_WP_FACTOR_LMH);
	}
	return (-1);
}

int AreAllTheSame (
	int	*array,
	int	len)
{
	int	i;

	for (i = 1; i < len; i++)
		if (array[i-1] != array[i]) return (FALSE);
	return (TRUE);
}



@implementation WorkerWindowsController

- (id)init
{
	if (![super initWithWindowNibName:@"WorkerWindows"]) return nil;
	workerData = [[NSMutableArray alloc] init];
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	int	i;

// In theory, the maximum number of workers should be the number of cpu cores.
// However, local.ini could specify more worker threads than cpu cores (for
// example, when local.ini is copied from a dual-core to a single-core machine).
// We must let the user manipulate the options on these worker threads that
// don't have a CPU to run on.

	[self setNumWorkersMax:max (NUM_WORKER_THREADS, NUM_CPUS)];
	[self setNumWorkersEnabled:(numWorkersMax > 1)];

// delete old rows
	
	while ([workerData count])
		[workerDataArrayController removeObjectAtArrangedObjectIndex:0];

// create and configure a row for each worker

	numWorkers = NUM_WORKER_THREADS;
	for (i = 0; i < numWorkers; i++) {
		WorkerData *newRow = [[WorkerData alloc] init];
		[newRow setWorkerNumber:[[NSString alloc] initWithFormat:@"Worker #%d", i+1]];
		[newRow setTypeOfWork:map_work_pref_to_sel(WORK_PREFERENCE[i])];
		[newRow setMultithreading:CORES_PER_TEST[i]];
		[newRow setMultithreadingMax:NUM_CPUS];
		[newRow setMultithreadingEnabled:(CORES_PER_TEST[i] > 1 || numWorkers < numWorkersMax)];
		[workerDataArrayController addObject:newRow];
		[newRow release];
	}

// finish off the initialization

	[self setNumWorkers:NUM_WORKER_THREADS];
	[self setPriority:PRIORITY];
}

/* Remember and clear the START_UP_IN_PROGRESS flag.  Why?  This routine */
/* is our only chance to capture close/cancel, but is also called if we hit apply. */
/* Thus, process the request as if this is a cancel and then undo the cancel in the ok method. */

- (void)windowWillClose:(NSNotification *)notification
{
	startupInProgress = STARTUP_IN_PROGRESS;
	STARTUP_IN_PROGRESS = 0;
}

- (int)numWorkers
{
	return numWorkers;
}

- (void)setNumWorkers:(int) _value
{
	int	i;

	if (numWorkers < _value) {
		for (i = numWorkers; i < _value; i++) {	// create and configure a row for each new worker
			WorkerData *newRow = [[WorkerData alloc] init];
			[newRow setWorkerNumber:[[NSString alloc] initWithFormat:@"Worker #%d", i+1]];
			[newRow setTypeOfWork:map_work_pref_to_sel(PRIMENET_WP_WHATEVER)];
			[newRow setMultithreading:1];
			[newRow setMultithreadingMax:NUM_CPUS];
			[newRow setMultithreadingEnabled:(_value < numWorkersMax)];
			[workerDataArrayController addObject:newRow];
			[newRow release];
		}
	} else if (numWorkers > _value) {
		for (i = numWorkers; i > _value; i--)	// delete rows
			[workerDataArrayController removeObjectAtArrangedObjectIndex:i-1];
	}
	numWorkers = _value;
}

- (IBAction)ok:(id)sender
{
	int	i, tot_cores, work_prefs[MAX_NUM_WORKER_THREADS];
	int	num_cpus[MAX_NUM_WORKER_THREADS];
	int	restart = FALSE;
	int	new_options = FALSE;

	[[self window] makeFirstResponder:nil];			// End any active text field edits

/* Make sure user has not allocated too many cores */

	for (i = 0, tot_cores = 0; i < numWorkers; i++) {	// examine each worker row
		WorkerData *row = [workerData objectAtIndex:i];
		tot_cores += [row multithreading];
	}
	if (tot_cores > NUM_CPUS) {
		NSAlert *alert = [[NSAlert alloc] init];
		[alert addButtonWithTitle:@"OK"];
		[alert setMessageText:@"There are not enough CPU cores to run all the workers.  Reduce the number of workers or CPU cores to use."];
		[alert setAlertStyle:NSWarningAlertStyle];
		[alert runModal];
		[alert release];
		return;
	}

/* If user changed the number of worker threads, then make the */
/* necessary changes.  Restart worker threads so that we are running */
/* the correct number of worker threads. */

	if (NUM_WORKER_THREADS != numWorkers) {
		NUM_WORKER_THREADS = numWorkers;
		IniWriteInt (LOCALINI_FILE, "WorkerThreads", NUM_WORKER_THREADS);
		new_options = TRUE;
		restart = TRUE;
	}

/* If user changed the priority of worker threads, then change */
/* the INI file.  Restart worker threads so that they are running at */
/* the new priority. */

	if (PRIORITY != priority) {
		PRIORITY = priority;
		IniWriteInt (INI_FILE, "Priority", PRIORITY);
		new_options = TRUE;
		restart = TRUE;
	}

/* If the user changed any of the work preferences record it in the INI file */
/* and tell the server */

	for (i = 0; i < numWorkers; i++) {	// examine each worker row
		WorkerData *row = [workerData objectAtIndex:i];
		work_prefs[i] = map_sel_to_work_pref ([row typeOfWork]);
		num_cpus[i] = [row multithreading];
	}
	if (AreAllTheSame (work_prefs, numWorkers)) {
		if (! PTOIsGlobalOption (WORK_PREFERENCE) || WORK_PREFERENCE[0] != work_prefs[0]) {
			PTOSetAll (INI_FILE, "WorkPreference", NULL, WORK_PREFERENCE, work_prefs[0]);
			new_options = TRUE;
		}
	} else {
		for (i = 0; i < numWorkers; i++) {
			if (WORK_PREFERENCE[i] == work_prefs[i]) continue;
			PTOSetOne (INI_FILE, "WorkPreference", NULL, WORK_PREFERENCE, i, work_prefs[i]);
			new_options = TRUE;
		}
	}

/* If the user changed any of the cores_per_test record it in the INI file */

	if (AreAllTheSame (num_cpus, numWorkers))
		PTOSetAll (LOCALINI_FILE, "CoresPerTest", NULL, CORES_PER_TEST, num_cpus[0]);
	else for (i = 0; i < numWorkers; i++)
		PTOSetOne (LOCALINI_FILE, "CoresPerTest", NULL, CORES_PER_TEST, i, num_cpus[i]);

/* Spool message if any options changed, restart workers if necessary */

	if (new_options) spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	if (restart) stop_workers_for_restart ();

/* Close the window */

	[[self window] performClose:self];

/* Finish off the startup sequence */

	if (startupInProgress) {
		set_comm_timers ();
		LaunchWorkerThreads (ALL_WORKERS, FALSE);
	}
}

@synthesize numWorkersMax;
@synthesize numWorkersEnabled;
@synthesize priority;
@synthesize workerData;

@end


@implementation WorkerData

- (int)typeOfWork
{
	return typeOfWork;
}

- (void)setTypeOfWork:(int) _value
{
	int	min_cores;

	min_cores = min_cores_for_work_pref (map_sel_to_work_pref (_value));
	if (multithreading < min_cores) [self setMultithreading:min_cores];
	typeOfWork = _value;
}

- (int)multithreading
{
	return multithreading;
}

- (void)setMultithreading:(int) _value
{
	int	min_cores;

	min_cores = min_cores_for_work_pref (map_sel_to_work_pref (typeOfWork));
	if (_value < min_cores) _value = min_cores;
	multithreading = _value;
}

@synthesize workerNumber;
@synthesize multithreadingMax;
@synthesize multithreadingEnabled;

@end
