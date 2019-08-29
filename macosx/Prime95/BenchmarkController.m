//
//  BenchmarkController.m
//  Prime95
//
//  Created by George Woltman on 1/25/17.
//  Copyright 2017-2019 Mersenne Research, Inc. All rights reserved.
//

#import "BenchmarkController.h"
#include "prime95.h"

@implementation BenchmarkController

- (id)init
{
	if (![super initWithWindowNibName:@"Benchmark"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	char	default_cores_string[80], default_workers_string[80];
	int	i, vals[4], numvals;

	[self setBenchType:0];
	[self setMinFFT:IniGetInt (INI_FILE, "MinBenchFFT", 2048)];
	[self setMaxFFT:IniGetInt (INI_FILE, "MaxBenchFFT", 8192)];
	[self setMinBenchableFFT:4];
	[self setMaxBenchableFFT:51200];
	[self setErrchk:ERRCHK];			//[self setAllComplex:IniGetInt (INI_FILE, "BenchErrorCheck", 0)];
	[self setAllComplex:0];				//[self setAllComplex:IniGetInt (INI_FILE, "BenchComplex", 0)];
	[self setLimitFFTSizes:0];			//[self setLimitFFTSizes:IniGetInt (INI_FILE, "OnlyBench5678", 1)];
	sprintf (default_cores_string, "%lu", NUM_CPUS);
	[self setBenchCores:[[NSString alloc] initWithFormat:@"%s", default_cores_string]];
	[self setHyperthreading:IniGetInt (INI_FILE, "BenchHyperthreads", 1)];

	// Init throughput dialog box entries
	[self setAllFFTImpl:IniGetInt (INI_FILE, "AllBench", 0)];
	[self setBenchTime:IniGetInt (INI_FILE, "BenchTime", 15)];
	[self setMinBenchTime:5];
	[self setMaxBenchTime:60];
	// If testing all FFT implementations. then default to the current num_workers.
	// Otherwise, assume user is trying to figure out how many workers to run and form a string
	// with the most common best values for number of workers: 1, num_threading_nodes, num_cores, num_workers
	numvals = 0;
	sorted_add_unique (vals, &numvals, NUM_WORKER_THREADS);
	if (!allFFTImpl) {
		sorted_add_unique (vals, &numvals, 1);
		sorted_add_unique (vals, &numvals, NUM_THREADING_NODES);
		sorted_add_unique (vals, &numvals, NUM_CPUS);
	}
	sprintf (default_workers_string, "%d", vals[0]);
	for (i = 1; i < numvals; i++) sprintf (default_workers_string + strlen (default_workers_string), ",%d", vals[i]);
	[self setBenchWorkers:[[NSString alloc] initWithFormat:@"%s", default_workers_string]];

	[self setEnableds];
}

- (void)setEnableds
{
	[self setMinmaxFFTEnabled: (benchType != 2)];
	[self setLimitFFTSizesEnabled: (benchType != 2 && minFFT != maxFFT && !allFFTImpl)];
	[self setBenchCoresEnabled: (NUM_CPUS > 1)];
	[self setHyperthreadingEnabled: (CPU_HYPERTHREADS > 1)];
	[self setBenchWorkersEnabled: (benchType == 0 && NUM_CPUS > 1)];
	[self setAllFFTImplEnabled: (benchType == 0)];
	[self setBenchTimeEnabled: (benchType == 0)];
}

- (int)benchType
{
	return benchType;
}

- (void)setBenchType:(int) _value
{
	benchType = _value;
	[self setEnableds];
}

- (int)minFFT
{
	return minFFT;
}

- (void)setMinFFT:(int) _value
{
	minFFT = _value;
	[self setEnableds];
}

- (int)maxFFT
{
	return maxFFT;
}

- (void)setMaxFFT:(int) _value
{
	maxFFT = _value;
	[self setEnableds];
}

- (int)allFFTImpl
{
	return allFFTImpl;
}

- (void)setAllFFTImpl:(int) _value
{
	allFFTImpl = _value;
	if (allFFTImpl) [self setLimitFFTSizes: FALSE];
	[self setEnableds];
}

@synthesize minmaxFFTEnabled;
@synthesize minBenchableFFT;
@synthesize maxBenchableFFT;
@synthesize errchk;
@synthesize allComplex;
@synthesize limitFFTSizes;
@synthesize limitFFTSizesEnabled;
@synthesize benchCores;
@synthesize benchCoresEnabled;
@synthesize hyperthreading;
@synthesize hyperthreadingEnabled;
@synthesize benchWorkers;
@synthesize benchWorkersEnabled;
@synthesize allFFTImplEnabled;
@synthesize benchTime;
@synthesize benchTimeEnabled;
@synthesize minBenchTime;
@synthesize maxBenchTime;

- (IBAction)ok:(id)sender
{
	[[self window] makeFirstResponder:nil];			// End any active text field edits

	if (benchType != 2) {
		IniWriteInt (INI_FILE, "MinBenchFFT", minFFT);
		IniWriteInt (INI_FILE, "MaxBenchFFT", maxFFT);
		IniWriteInt (INI_FILE, "BenchErrorCheck", errchk);
		IniWriteInt (INI_FILE, "BenchAllComplex", allComplex ? 2 : 0);
		IniWriteInt (INI_FILE, "OnlyBench5678", limitFFTSizes);
	}
	IniWriteString (INI_FILE, "BenchCores", [benchCores UTF8String]);
	IniWriteInt (INI_FILE, "BenchHyperthreads", hyperthreading);
	if (benchType == 0) {
		IniWriteString (INI_FILE, "BenchWorkers", [benchWorkers UTF8String]);
		IniWriteInt (INI_FILE, "AllBench", allFFTImpl);
		IniWriteInt (INI_FILE, "BenchTime", benchTime);
	}
	LaunchBench (benchType);

	[[self window] performClose:self];
}

@end
