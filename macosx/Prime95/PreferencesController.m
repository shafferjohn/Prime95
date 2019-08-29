//
//  PreferencesController.m
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009-2014 Mersenne Research, Inc. All rights reserved.
//

#import "PreferencesController.h"
#include "prime95.h"

@implementation PreferencesController

- (id)init
{
	if (![super initWithWindowNibName:@"Preferences"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setIterationsOutput:ITER_OUTPUT];
	[self setIterationsResultsFile:ITER_OUTPUT_RES];
	[self setWriteSaveFileMinutes:DISK_WRITE_TIME];
	[self setModemRetryMinutes:MODEM_RETRY_TIME];
	[self setModemRetryMinutesEnabled:(USE_PRIMENET && DIAL_UP)];
	[self setPrimenetOptionsEnabled:USE_PRIMENET];
	[self setNetworkRetryMinutes:NETWORK_RETRY_TIME];
	[self setDaysOfWork:DAYS_OF_WORK];
	[self setDaysBetweenEndDates:DAYS_BETWEEN_CHECKINS];
	[self setNumberOfBackupFiles:NUM_BACKUP_FILES];
	[self setMakeNoise:!SILENT_VICTORY];
	[self setRunOnBattery:RUN_ON_BATTERY];
	[self setDefeatPowerSave:DEFEAT_POWER_SAVE];
}

- (int)iterationsOutput
{
	return iterationsOutput;
}

- (void)setIterationsOutput:(int) _value
{
	if (ITER_OUTPUT != _value) {
		ITER_OUTPUT = _value;
		IniWriteInt (INI_FILE, "OutputIterations", ITER_OUTPUT);
	}
	iterationsOutput = _value;
}

- (int)iterationsResultsFile
{
	return iterationsResultsFile;
}

- (void)setIterationsResultsFile:(int) _value
{
	if (ITER_OUTPUT_RES != _value) {
		ITER_OUTPUT_RES = _value;
		IniWriteInt (INI_FILE, "ResultsFileIterations", ITER_OUTPUT_RES);
	}
	iterationsResultsFile = _value;
}

- (int)writeSaveFileMinutes
{
	return writeSaveFileMinutes;
}

- (void)setWriteSaveFileMinutes:(int) _value
{
	if (DISK_WRITE_TIME != _value) {
		DISK_WRITE_TIME = _value;
		IniWriteInt (INI_FILE, "DiskWriteTime", DISK_WRITE_TIME);
	}
	writeSaveFileMinutes = _value;
}

- (int)modemRetryMinutes
{
	return modemRetryMinutes;
}

- (void)setModemRetryMinutes:(int) _value
{
	if (MODEM_RETRY_TIME != _value) {
		MODEM_RETRY_TIME = _value;
		IniWriteInt (INI_FILE, "NetworkRetryTime", MODEM_RETRY_TIME);
	}
	modemRetryMinutes = _value;
}

- (int)networkRetryMinutes
{
	return networkRetryMinutes;
}

- (void)setNetworkRetryMinutes:(int) _value
{
	if (NETWORK_RETRY_TIME != _value) {
		NETWORK_RETRY_TIME = _value;
		IniWriteInt (INI_FILE, "NetworkRetryTime2", NETWORK_RETRY_TIME);
	}
	networkRetryMinutes = _value;
}

- (int)daysOfWork
{
	return daysOfWork;
}

- (void)setDaysOfWork:(int) _value
{
	if (DAYS_OF_WORK != _value) {
		DAYS_OF_WORK = _value;
		IniWriteInt (INI_FILE, "DaysOfWork", DAYS_OF_WORK);
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
	daysOfWork = _value;
}

- (int)daysBetweenEndDates
{
	return daysBetweenEndDates;
}

- (void)setDaysBetweenEndDates:(int) _value
{
	if (DAYS_BETWEEN_CHECKINS != _value) {
		DAYS_BETWEEN_CHECKINS = _value;
		IniWriteInt (INI_FILE, "DaysBetweenCheckins", (int) DAYS_BETWEEN_CHECKINS);
	}
	daysBetweenEndDates = _value;
}

- (int)numberOfBackupFiles
{
	return numberOfBackupFiles;
}

- (void)setNumberOfBackupFiles:(int) _value
{
	if (NUM_BACKUP_FILES != _value) {
		NUM_BACKUP_FILES = _value;
		IniWriteInt (INI_FILE, "NumBackupFiles", NUM_BACKUP_FILES);
	}
	numberOfBackupFiles = _value;
}

- (int)makeNoise
{
	return makeNoise;
}

- (void)setMakeNoise:(int) _value
{
	if (SILENT_VICTORY != !_value) {
		SILENT_VICTORY = !_value;
		IniWriteInt (INI_FILE, "SilentVictory", SILENT_VICTORY);
	}
	makeNoise = _value;
}

- (int)runOnBattery
{
	return runOnBattery;
}

- (void)setRunOnBattery:(int) _value
{
	if (RUN_ON_BATTERY != _value) {
		RUN_ON_BATTERY = _value;
		IniWriteInt (LOCALINI_FILE, "RunOnBattery", RUN_ON_BATTERY);
		run_on_battery_changed ();
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
	runOnBattery = _value;
}

- (int)defeatPowerSave
{
	return defeatPowerSave;
}

- (void)setDefeatPowerSave:(int) _value
{
	if (DEFEAT_POWER_SAVE != _value) {
		DEFEAT_POWER_SAVE = _value;
		IniWriteInt (LOCALINI_FILE, "DefeatPowerSave", DEFEAT_POWER_SAVE);
		stop_workers_for_restart();
	}
	defeatPowerSave = _value;
}

@synthesize modemRetryMinutesEnabled;
@synthesize primenetOptionsEnabled;

@end
