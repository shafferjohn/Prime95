//
//  ContinueController.m
//  Prime95
//
//  Created by George Woltman on 4/19/09.
//  Copyright 2009-2019 Mersenne Research, Inc. All rights reserved.
//

#import "ContinueController.h"
#include "prime95.h"

@implementation ContinueController

- (id)init
{
	if (![super initWithWindowNibName:@"Continue"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setStartAllWorkers:YES];
	[self setWorkerNumber:1];
	[self setWorkerNumberMax:WORKER_THREADS_ACTIVE > NUM_WORKER_THREADS ? WORKER_THREADS_ACTIVE : NUM_WORKER_THREADS];
}

@synthesize startAllWorkers;
@synthesize workerNumber;
@synthesize workerNumberMax;

- (IBAction)ok:(id)sender
{
	[[self window] makeFirstResponder:nil];			// End any active text field edits

	if (startAllWorkers)
		LaunchWorkerThreads (ALL_WORKERS, FALSE);
	else
		LaunchWorkerThreads (workerNumber - 1, FALSE);
	[[self window] performClose:self];
}

@end
