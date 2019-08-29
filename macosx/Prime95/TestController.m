//
//  TestController.m
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009-2010 Mersenne Research, Inc. All rights reserved.
//

#import "TestController.h"
#include "prime95.h"

@implementation TestController

- (id)init
{
	if (![super initWithWindowNibName:@"Test"]) return nil;
	exponentToTest = 5;
	exponentToTestMax = MAX_PRIME_SSE2;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setWorkerNumber:1];
	[self setWorkerNumberMax:NUM_WORKER_THREADS];
	[self setWorkerNumberEnabled:(NUM_WORKER_THREADS > 1)];
}

@synthesize workerNumber;
@synthesize workerNumberMax;
@synthesize workerNumberEnabled;
@synthesize exponentToTest;
@synthesize exponentToTestMax;

- (IBAction)ok:(id)sender
{
	struct work_unit w;

	[[self window] makeFirstResponder:nil];			// End any active text field edits

	memset (&w, 0, sizeof (w));
	w.work_type = WORK_ADVANCEDTEST;
	w.k = 1.0;
	w.b = 2;
	w.n = exponentToTest;
	w.c = -1;
	addWorkToDoLine (workerNumber - 1, &w);
	if (WORKER_THREADS_ACTIVE)
		stop_worker_for_advanced_test (workerNumber - 1);
	else
		LaunchWorkerThreads (ALL_WORKERS, FALSE);
	[[self window] performClose:self];
}

@end
