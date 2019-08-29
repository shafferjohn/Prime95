//
//  ManualCommunicationController.h
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface ManualCommunicationController : NSWindowController {
	int	manualContact;
}

@property(readwrite, assign) int manualContact;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end
