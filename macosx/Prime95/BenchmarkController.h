//
//  BenchmarkController.h
//  Prime95
//
//  Created by George Woltman on 1/25/17.
//  Copyright 2017 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface BenchmarkController : NSWindowController {
	int	benchType;
	int	minFFT;
	int	maxFFT;
	int	minmaxFFTEnabled;
	int	minBenchableFFT;
	int	maxBenchableFFT;
	int	errchk;
	int	allComplex;
	int	limitFFTSizes;
	int	limitFFTSizesEnabled;
	NSString *benchCores;
	int	benchCoresEnabled;
	int	hyperthreading;
	int	hyperthreadingEnabled;
	NSString *benchWorkers;
	int	benchWorkersEnabled;
	int	allFFTImpl;
	int	allFFTImplEnabled;
	int	benchTime;
	int	benchTimeEnabled;
	int	minBenchTime;
	int	maxBenchTime;
}

@property(readwrite, assign) int benchType;
@property(readwrite, assign) int minFFT;
@property(readwrite, assign) int maxFFT;
@property(readwrite, assign) int minmaxFFTEnabled;
@property(readwrite, assign) int minBenchableFFT;
@property(readwrite, assign) int maxBenchableFFT;
@property(readwrite, assign) int errchk;
@property(readwrite, assign) int allComplex;
@property(readwrite, assign) int limitFFTSizes;
@property(readwrite, assign) int limitFFTSizesEnabled;
@property(readwrite, retain) NSString *benchCores;
@property(readwrite, assign) int benchCoresEnabled;
@property(readwrite, assign) int hyperthreading;
@property(readwrite, assign) int hyperthreadingEnabled;
@property(readwrite, retain) NSString *benchWorkers;
@property(readwrite, assign) int benchWorkersEnabled;
@property(readwrite, assign) int allFFTImpl;
@property(readwrite, assign) int allFFTImplEnabled;
@property(readwrite, assign) int benchTime;
@property(readwrite, assign) int benchTimeEnabled;
@property(readwrite, assign) int minBenchTime;
@property(readwrite, assign) int maxBenchTime;

	- (void)reInit;
	- (void)setEnableds;
- (IBAction)ok:(id)sender;

@end
