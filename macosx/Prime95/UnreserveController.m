//
//  UnreserveController.m
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009-2010 Mersenne Research, Inc. All rights reserved.
//

#import "UnreserveController.h"
#include "prime95.h"

@implementation UnreserveController

- (id)init
{
	if (![super initWithWindowNibName:@"Unreserve"]) return nil;
	exponent = 50000000;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
}

@synthesize exponent;

- (IBAction)ok:(id)sender
{
	[[self window] makeFirstResponder:nil];			// End any active text field edits

	unreserve (exponent);
	[[self window] performClose:self];
}

@end
