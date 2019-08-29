//
//  TestController.h
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface TestController : NSWindowController {
	int	workerNumber;
	int	workerNumberMax;
	int	workerNumberEnabled;
	int	exponentToTest;
	int	exponentToTestMax;
}

@property(readwrite, assign) int workerNumber;
@property(readwrite, assign) int workerNumberMax;
@property(readwrite, assign) int workerNumberEnabled;
@property(readwrite, assign) int exponentToTest;
@property(readonly) int exponentToTestMax;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
