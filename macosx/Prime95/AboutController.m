//
//  AboutController.m
//  Prime95
//
//  Created by George Woltman on 4/17/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "AboutController.h"
#include "prime95.h"

@implementation AboutController

- (id)init
{
	if (![super initWithWindowNibName:@"About"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	char	app_string[100];
	generate_application_string (app_string);
	[textField setStringValue:[[NSString alloc] initWithFormat:@"%s", app_string]];
}

- (IBAction)linkToMersenneOrg:(id)sender
{
	NSURL *url = [NSURL URLWithString:@"http://mersenne.org"];
	[[NSWorkspace sharedWorkspace] openURL:url];
}

@end
