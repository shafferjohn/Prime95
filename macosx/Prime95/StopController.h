//
//  StopController.h
//  Prime95
//
//  Created by George Woltman on 4/19/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface StopController : NSWindowController {
	int	stopAllWorkers;
	int	workerNumber;
	int	workerNumberMax;

}

@property(readwrite, assign) int stopAllWorkers;
@property(readwrite, assign) int workerNumber;
@property(readwrite, assign) int workerNumberMax;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
