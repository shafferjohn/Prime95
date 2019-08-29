//
//  Pminus1Controller.m
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009-2010 Mersenne Research, Inc. All rights reserved.
//

#import "Pminus1Controller.h"
#include "prime95.h"

@implementation Pminus1Controller

- (id)init
{
	if (![super initWithWindowNibName:@"Pminus1"]) return nil;
	k = 1.0;
	b = 2;
	n = 1061;
	nMax = MAX_PRIME_SSE2;
	c = -1;
	bound1 = 1000000.0;
	bound2 = 0.0;

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

- (IBAction)ok:(id)sender
{
	struct work_unit w;

	[[self window] makeFirstResponder:nil];			// End any active text field edits

	memset (&w, 0, sizeof (w));
	w.work_type = WORK_PMINUS1;
	w.k = k;
	w.b = b;
	w.n = n;
	w.c = c;
	w.B1 = bound1;
	w.B2_start = 0;
	w.B2 = bound2;
	addWorkToDoLine (workerNumber - 1, &w);

	if (! WORKER_THREADS_ACTIVE)
		LaunchWorkerThreads (ALL_WORKERS, FALSE);

	[[self window] performClose:self];
}

@end
