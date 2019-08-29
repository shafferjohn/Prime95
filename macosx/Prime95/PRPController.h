//
//  PRPController.h
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface PRPController : NSWindowController {

	int	workerNumber;
	int	workerNumberMax;
	int	workerNumberEnabled;
	double	k;
	int	b;
	int	n;
	int	nMax;
	int	c;
}

@property(readwrite, assign) int workerNumber;
@property(readwrite, assign) int workerNumberMax;
@property(readwrite, assign) int workerNumberEnabled;
@property(readwrite, assign) double k;
@property(readwrite, assign) int b;
@property(readwrite, assign) int n;
@property(readonly) int nMax;
@property(readwrite, assign) int c;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
