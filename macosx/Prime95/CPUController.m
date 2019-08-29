//
//  CPUController.m
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "CPUController.h"
#import "AppController.h"
#include "prime95.h"

@implementation CPUController

- (id)init
{
	if (![super initWithWindowNibName:@"CPU"]) return nil;
	memoryMax = 0.9 * physical_memory ();
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	int	editable;
	char	buf[512];

	[self setHoursPerDay:CPU_HOURS];

	editable = read_memory_settings (&origDayMemory, &origNightMemory, &origDayStartTime, &origDayEndTime);
	[self setDayMemory:origDayMemory];
	[self setNightMemory:origNightMemory];
	[self setMemoryEnabled:editable];

	minutesToStr (origDayStartTime, buf);
	[self setDayStartTime:[[NSString alloc] initWithFormat:@"%s", buf]];
	minutesToStr (origDayEndTime, buf);
	[self setDayEndTime:[[NSString alloc] initWithFormat:@"%s", buf]];

	getCpuDescription (buf, 0);
	[textField setStringValue:[[NSString alloc] initWithFormat:@"%s", buf]];
}

/* Start next dialog in startup chain */

- (void)windowWillClose:(NSNotification *)notification
{
	if (STARTUP_IN_PROGRESS)
		[myAppController testWorkerWindows:nil];
}

- (int)hoursPerDay
{
	return hoursPerDay;
}

- (void)setHoursPerDay:(int) _value
{
	if (CPU_HOURS != _value) {
		CPU_HOURS = _value;
		IniWriteInt (LOCALINI_FILE, "CPUHours", CPU_HOURS);
		ROLLING_AVERAGE = 1000;
		IniWriteInt (LOCALINI_FILE, "RollingAverage", 1000);
		IniWriteInt (LOCALINI_FILE, "RollingStartTime", 0);
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		delete_timed_event (TE_COMM_SERVER);
		UpdateEndDates ();
	}
	hoursPerDay = _value;
}

- (int)dayMemory
{
	return dayMemory;
}

- (void)setDayMemory:(int) _value
{
	if (origDayMemory != (unsigned int) _value) {
		origDayMemory = _value;
		write_memory_settings (origDayMemory, origNightMemory, origDayStartTime, origDayEndTime);
		mem_settings_have_changed ();
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
	dayMemory = _value;
}

- (int)nightMemory
{
	return nightMemory;
}

- (void)setNightMemory:(int) _value
{
	if (origNightMemory != (unsigned int) _value) {
		origNightMemory = _value;
		write_memory_settings (origDayMemory, origNightMemory, origDayStartTime, origDayEndTime);
		mem_settings_have_changed ();
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
	nightMemory = _value;
}

- (NSString *)dayStartTime
{
	return dayStartTime;
}

- (void)setDayStartTime:(NSString *) _value
{
	unsigned int new_day_start_time;

	new_day_start_time = strToMinutes ([_value UTF8String]);
	if (origDayStartTime != new_day_start_time) {
		origDayStartTime = new_day_start_time;
		write_memory_settings (origDayMemory, origNightMemory, origDayStartTime, origDayEndTime);
		mem_settings_have_changed ();
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
	[_value retain];
	[dayStartTime release];
	dayStartTime = _value;
}

- (NSString *)dayEndTime
{
	return dayEndTime;
}

- (void)setDayEndTime:(NSString *) _value
{
	unsigned int new_day_end_time;

	new_day_end_time = strToMinutes ([_value UTF8String]);
	if (origDayEndTime != new_day_end_time) {
		origDayEndTime = new_day_end_time;
		write_memory_settings (origDayMemory, origNightMemory, origDayStartTime, origDayEndTime);
		mem_settings_have_changed ();
		spoolMessage (PRIMENET_PROGRAM_OPTIONS, NULL);
	}
	[_value retain];
	[dayEndTime release];
	dayEndTime = _value;
}

@synthesize memoryMax;
@synthesize memoryEnabled;

@end
