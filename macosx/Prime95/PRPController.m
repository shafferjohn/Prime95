//
//  PRPController.m
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009-2010 Mersenne Research, Inc. All rights reserved.
//

#import "PRPController.h"
#include "prime95.h"

@implementation PRPController

- (id)init
{
	if (![super initWithWindowNibName:@"PRP"]) return nil;
	k = 4605.0;
	b = 2;
	n = 3313;
	nMax = MAX_PRIME_SSE2;
	c = 1;

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

- (IBAction)ok:(id)sender
{
	struct work_unit w;

	[[self window] makeFirstResponder:nil];			// End any active text field edits

	memset (&w, 0, sizeof (w));
	w.work_type = WORK_PRP;
	w.k = k;
	w.b = b;
	w.n = n;
	w.c = c;
	addWorkToDoLine (workerNumber - 1, &w);

	if (! WORKER_THREADS_ACTIVE)
		LaunchWorkerThreads (ALL_WORKERS, FALSE);

	[[self window] performClose:self];
}

@end
