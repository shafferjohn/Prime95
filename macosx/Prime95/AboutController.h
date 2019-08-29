//
//  AboutController.h
//  Prime95
//
//  Created by George Woltman on 4/17/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface AboutController : NSWindowController {
	IBOutlet NSTextField *textField;
}

- (void)reInit;
- (IBAction)linkToMersenneOrg:(id)sender;

@end
