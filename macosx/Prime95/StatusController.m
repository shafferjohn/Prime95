//
//  StatusController.m
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "StatusController.h"
#include "prime95.h"

@implementation StatusController

- (id)init
{
	if (![super initWithWindowNibName:@"Status"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	char	buf[20000];
	NSRange endRange;

	rangeStatusMessage (buf, sizeof (buf));

	endRange.location = 0;
	endRange.length = [[textView textStorage] length];
	[textView replaceCharactersInRange:endRange withString:[[NSString alloc] initWithFormat:@"%s", buf]];
}

@end
