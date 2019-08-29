//
//  UnreserveController.h
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface UnreserveController : NSWindowController {

	int	exponent;
}

@property(readwrite, assign) int exponent;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
