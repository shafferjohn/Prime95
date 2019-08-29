//
//  PreferencesController.h
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009-2014 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface PreferencesController : NSWindowController {
	int	iterationsOutput;
	int	iterationsResultsFile;
	int	writeSaveFileMinutes;
	int	modemRetryMinutes;
	int	modemRetryMinutesEnabled;
	int	primenetOptionsEnabled;
	int	networkRetryMinutes;
	int	daysOfWork;
	int	daysBetweenEndDates;
	int	numberOfBackupFiles;
	int	makeNoise;
	int	runOnBattery;
        int     defeatPowerSave;
}

@property(readwrite, assign) int iterationsOutput;
@property(readwrite, assign) int iterationsResultsFile;
@property(readwrite, assign) int writeSaveFileMinutes;
@property(readwrite, assign) int modemRetryMinutes;
@property(readwrite, assign) int modemRetryMinutesEnabled;
@property(readwrite, assign) int primenetOptionsEnabled;
@property(readwrite, assign) int networkRetryMinutes;
@property(readwrite, assign) int daysOfWork;
@property(readwrite, assign) int daysBetweenEndDates;
@property(readwrite, assign) int numberOfBackupFiles;
@property(readwrite, assign) int makeNoise;
@property(readwrite, assign) int runOnBattery;
@property(readwrite, assign) int defeatPowerSave;

- (void)reInit;


@end
