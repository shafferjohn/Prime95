//
//  ECMController.h
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface ECMController : NSWindowController {
	int	workerNumber;
	int	workerNumberMax;
	int	workerNumberEnabled;
	double	k;
	int	b;
	int	n;
	int	nMax;
	int	c;
	double	bound1;
	double	bound2;
	int	numberOfCurves;
}

@property(readwrite, assign) int workerNumber;
@property(readwrite, assign) int workerNumberMax;
@property(readwrite, assign) int workerNumberEnabled;
@property(readwrite, assign) double k;
@property(readwrite, assign) int b;
@property(readwrite, assign) int n;
@property(readonly) int nMax;
@property(readwrite, assign) int c;
@property(readwrite, assign) double bound1;
@property(readwrite, assign) double bound2;
@property(readwrite, assign) int numberOfCurves;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
