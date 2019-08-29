//
//  WindowController.h
//  Prime95
//
//  Created by George Woltman on 4/18/09.
//  Copyright 2009 Merenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface WindowController : NSWindowController {
	IBOutlet NSTextView *textView;
	NSString *baseTitle;
}

- (void)setFontSize:(int)newSize;

@end
