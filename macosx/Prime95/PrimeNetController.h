//
//  PrimeNetController.h
//  Prime95
//
//  Created by George Woltman on 4/26/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>
@class ConnectionController;

@interface PrimeNetController : NSWindowController {
	int	usePrimeNet;
	NSString *userID;
	NSString *computerName;
	ConnectionController *connectionController;
	int	startupInProgress;
}

@property(readwrite, assign) int usePrimeNet;
@property(readwrite, retain) NSString *userID;
@property(readwrite, retain) NSString *computerName;

- (void)reInit;
- (IBAction)connection:(id)sender;
- (IBAction)ok:(id)sender;
- (IBAction)linkToMersenneOrg:(id)sender;

@end
