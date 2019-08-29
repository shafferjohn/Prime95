//
//  TortureTestController.h
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009-2019 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface TortureTestController : NSWindowController {
	int	tortureType;
	int	numberOfThreads;
	int	numberOfThreadsMin;
	int	numberOfThreadsMax;
	int	numberOfThreadsEnabled;
        int     smallFFTsEnabled;
        int     mediumFFTsEnabled;
	int	customSettingsEnabled;
	int	customMemoryEnabled;
	int	minFFTSize;
	int	maxFFTSize;
	int	runFFTsInPlace;
	int	memoryToUse;
	int	timeToRunEachFFT;
	int	blendMemory;
	int	disableAVX512;
	int	disableFMA3;
	int	disableAVX;
	int	disableAVX512Enabled;
	int	disableFMA3Enabled;
	int	disableAVXEnabled;
}

@property(readwrite, assign) int tortureType;
@property(readwrite, assign) int numberOfThreads;
@property(readwrite, assign) int numberOfThreadsMin;
@property(readwrite, assign) int numberOfThreadsMax;
@property(readwrite, assign) int numberOfThreadsEnabled;
@property(readwrite, assign) int smallFFTsEnabled;
@property(readwrite, assign) int mediumFFTsEnabled;
@property(readwrite, assign) int customSettingsEnabled;
@property(readwrite, assign) int customMemoryEnabled;
@property(readwrite, assign) int minFFTSize;
@property(readwrite, assign) int maxFFTSize;
@property(readwrite, assign) int runFFTsInPlace;
@property(readwrite, assign) int memoryToUse;
@property(readwrite, assign) int timeToRunEachFFT;
@property(readwrite, assign) int disableAVX512;
@property(readwrite, assign) int disableFMA3;
@property(readwrite, assign) int disableAVX;
@property(readwrite, assign) int disableAVX512Enabled;
@property(readwrite, assign) int disableFMA3Enabled;
@property(readwrite, assign) int disableAVXEnabled;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
