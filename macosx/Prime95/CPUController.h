//
//  CPUController.h
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface CPUController : NSWindowController {
	IBOutlet NSTextField *textField;
	int	hoursPerDay;
	int	dayMemory;
	int	nightMemory;
	int	memoryMax;
	int	memoryEnabled;
	NSString *dayStartTime;
	NSString *dayEndTime;
	unsigned int origDayMemory;
	unsigned int origNightMemory;
	unsigned int origDayStartTime;
	unsigned int origDayEndTime;
}

@property(readwrite, assign) int hoursPerDay;
@property(readwrite, assign) int dayMemory;
@property(readwrite, assign) int nightMemory;
@property(readonly) int memoryMax;
@property(readwrite, assign) int memoryEnabled;
@property(readwrite, retain) NSString *dayStartTime;
@property(readwrite, retain) NSString *dayEndTime;

- (void)reInit;


@end
