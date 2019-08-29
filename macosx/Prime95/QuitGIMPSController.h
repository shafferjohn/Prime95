//
//  QuitGIMPSController.h
//  Prime95
//
//  Created by George Woltman on 5/1/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface QuitGIMPSController : NSWindowController {
	NSString *explanation;
}

@property(readwrite, retain) NSString *explanation;

- (void)reInit;
- (IBAction)quitLater:(id)sender;
- (IBAction)quitNow:(id)sender;

@end
