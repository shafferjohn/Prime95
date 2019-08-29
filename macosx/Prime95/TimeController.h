//
//  TimeController.h
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface TimeController : NSWindowController {

	int	exponentToTime;
	int	exponentToTimeMax;
	int	numberOfIterations;
}

@property(readwrite, assign) int exponentToTime;
@property(readonly) int exponentToTimeMax;
@property(readwrite, assign) int numberOfIterations;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
