//
//  ECMController.m
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009-2019 Mersenne Research, Inc. All rights reserved.
//

#import "ECMController.h"
#include "prime95.h"

@implementation ECMController

- (id)init
{
	if (![super initWithWindowNibName:@"ECM"]) return nil;
	k = 1.0;
	b = 2;
	n = 1277;
	nMax = MAX_PRIME_SSE2;
	c = -1;
	bound1 = 850.0e6;
	bound2 = 0.0;
	numberOfCurves = 10;

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
@synthesize k;
@synthesize b;
@synthesize n;
@synthesize nMax;
@synthesize c;
@synthesize bound1;
@synthesize bound2;
@synthesize numberOfCurves;

- (IBAction)ok:(id)sender
{
	struct work_unit w;

	[[self window] makeFirstResponder:nil];			// End any active text field edits

	memset (&w, 0, sizeof (w));
	w.work_type = WORK_ECM;
	w.k = k;
	w.b = b;
	w.n = n;
	w.c = c;
	w.B1 = bound1;
	w.B2_start = 0;
	w.B2 = bound2;
	w.curves_to_do = numberOfCurves;
	addWorkToDoLine (workerNumber - 1, &w);

	if (! WORKER_THREADS_ACTIVE)
		LaunchWorkerThreads (ALL_WORKERS, FALSE);

	[[self window] performClose:self];
}

@end
