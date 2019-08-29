//
//  TortureTestController.m
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009-2019 Mersenne Research, Inc. All rights reserved.
//

#import "TortureTestController.h"
#include "prime95.h"

@implementation TortureTestController

- (id)init
{
	if (![super initWithWindowNibName:@"TortureTest"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	int	mem, in_place_fft;

	[self setNumberOfThreads:NUM_CPUS * CPU_HYPERTHREADS];
	[self setNumberOfThreadsMin:1];
	[self setNumberOfThreadsMax:NUM_CPUS * CPU_HYPERTHREADS];
	[self setNumberOfThreadsEnabled:(NUM_CPUS * CPU_HYPERTHREADS > 1)];
	[self setTortureType:4];
	[self setSmallFFTsEnabled:(CPU_TOTAL_L3_CACHE_SIZE > 0)];
	[self setMediumFFTsEnabled:(CPU_TOTAL_L4_CACHE_SIZE > 0)];
	[self setCustomSettingsEnabled:NO];
	[self setCustomMemoryEnabled:NO];
	[self setMinFFTSize:4];
	[self setMaxFFTSize:(CPU_TOTAL_L4_CACHE_SIZE ? 32768 : 8192)];
	[self setDisableAVX512:!(CPU_FLAGS & CPU_AVX512F)];
	[self setDisableFMA3:!(CPU_FLAGS & CPU_FMA3)];
	[self setDisableAVX:!(CPU_FLAGS & CPU_AVX)];

	mem = physical_memory ();
	// New in 29.5 default to all but 2.5GB
	if (mem >= 4500) {
		blendMemory = GetSuggestedMemory (mem - 2500);
		in_place_fft = FALSE;
	} else if (mem >= 2000) {
		blendMemory = GetSuggestedMemory (1600);
		in_place_fft = FALSE;
	} else if (mem >= 500) {
		blendMemory = GetSuggestedMemory (mem - 256);
		in_place_fft = FALSE;
	} else if (mem >= 200) {
		blendMemory = GetSuggestedMemory (mem / 2);
		in_place_fft = TRUE;
	} else {
		blendMemory = 8;
		in_place_fft = TRUE;
	}
	[self setRunFFTsInPlace:in_place_fft];
	[self setMemoryToUse:blendMemory];
	[self setTimeToRunEachFFT:CPU_HYPERTHREADS * 3];
}

- (int)numberOfThreads
{
	return numberOfThreads;
}

- (void)setNumberOfThreads:(int) _value
{
	numberOfThreads = _value;
	[self setTortureType:tortureType];		// Set FFT lengths for new numberOfThreads
}

- (int)tortureType
{
	return tortureType;
}

- (void)setTortureType:(int) _value
{
	if (_value != 5) {				// Not custom
		int	minfft, maxfft;

		// Calculate default FFT sizes
		tortureTestDefaultSizes (_value, numberOfThreads, &minfft, &maxfft);
		if (minfft < 4) minfft = 4;
		if (maxfft < minfft) maxfft = minfft;
		if (maxfft > 32768) maxfft = 32768;
		[self setMinFFTSize:minfft];
		[self setMaxFFTSize:maxfft];

		// Assign other options
		[self setRunFFTsInPlace:(_value <= 2)];	// TRUE for L2/L3/L4 cache
		[self setMemoryToUse:(runFFTsInPlace ? 0 : blendMemory)];
		[self setTimeToRunEachFFT:(numberOfThreads > NUM_CPUS ? 6 : 3)];
		[self setCustomSettingsEnabled:NO];
		[self setCustomMemoryEnabled:NO];
	}
	else {						// Custom
		[self setCustomSettingsEnabled:YES];
		[self setCustomMemoryEnabled:!runFFTsInPlace];
	}

	tortureType = _value;
}

- (int)runFFTsInPlace
{
	return runFFTsInPlace;
}

- (void)setRunFFTsInPlace:(int) _value
{
	[self setCustomMemoryEnabled:(customMemoryEnabled && !_value)];

	runFFTsInPlace = _value;
}

- (int)disableAVX512
{
	return disableAVX512;
}

- (void)setDisableAVX512:(int) _value
{
	disableAVX512 = _value;
	[self setDisableAVX512Enabled:(CPU_FLAGS & CPU_AVX512F) && !disableFMA3];
	[self setDisableFMA3Enabled:(CPU_FLAGS & CPU_FMA3) && disableAVX512 && !disableAVX];
	[self setDisableAVXEnabled:(CPU_FLAGS & CPU_AVX) && disableFMA3];
}

- (int)disableFMA3
{
	return disableFMA3;
}

- (void)setDisableFMA3:(int) _value
{
	disableFMA3 = _value;
	[self setDisableAVX512Enabled:(CPU_FLAGS & CPU_AVX512F) && !disableFMA3];
	[self setDisableFMA3Enabled:(CPU_FLAGS & CPU_FMA3) && disableAVX512 && !disableAVX];
	[self setDisableAVXEnabled:(CPU_FLAGS & CPU_AVX) && disableFMA3];
}

- (int)disableAVX
{
	return disableAVX;
}

- (void)setDisableAVX:(int) _value
{
	disableAVX = _value;
	[self setDisableAVX512Enabled:(CPU_FLAGS & CPU_AVX512F) && !disableFMA3];
	[self setDisableFMA3Enabled:(CPU_FLAGS & CPU_FMA3) && disableAVX512 && !disableAVX];
	[self setDisableAVXEnabled:(CPU_FLAGS & CPU_AVX) && disableFMA3];
}

@synthesize numberOfThreadsMin;
@synthesize numberOfThreadsMax;
@synthesize numberOfThreadsEnabled;
@synthesize smallFFTsEnabled;
@synthesize mediumFFTsEnabled;
@synthesize customSettingsEnabled;
@synthesize customMemoryEnabled;
@synthesize minFFTSize;
@synthesize maxFFTSize;
@synthesize memoryToUse;
@synthesize timeToRunEachFFT;
@synthesize disableAVX512Enabled;
@synthesize disableFMA3Enabled;
@synthesize disableAVXEnabled;

- (IBAction)ok:(id)sender
{
	int	weak;

	[[self window] makeFirstResponder:nil];			// End any active text field edits

	IniWriteInt (INI_FILE, "MinTortureFFT", minFFTSize);
	IniWriteInt (INI_FILE, "MaxTortureFFT", maxFFTSize);
	if (runFFTsInPlace) memoryToUse = 8;
	IniWriteInt (INI_FILE, "TortureMem", memoryToUse);
	IniWriteInt (INI_FILE, "TortureTime", timeToRunEachFFT);
	weak = (disableAVX512 ? CPU_AVX512F : 0) + (disableFMA3 ? CPU_FMA3 : 0) + (disableAVX ? CPU_AVX : 0);
	IniWriteInt (INI_FILE, "TortureWeak", weak);
	LaunchTortureTest (numberOfThreads, FALSE);

	[[self window] performClose:self];
}

@end
