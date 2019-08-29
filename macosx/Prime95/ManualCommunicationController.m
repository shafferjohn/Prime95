//
//  ManualCommunicationController.m
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009-2010 Mersenne Research, Inc. All rights reserved.
//

#import "ManualCommunicationController.h"
#include "prime95.h"

@implementation ManualCommunicationController

- (id)init
{
	if (![super initWithWindowNibName:@"ManualCommunication"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setManualContact:MANUAL_COMM];
}

- (int)manualContact
{
	return manualContact;
}

- (void)setManualContact:(int) _value
{
	if ((MANUAL_COMM && !_value) || (!MANUAL_COMM && _value)) {
		MANUAL_COMM = _value;
		IniWriteInt (INI_FILE, "ManualComm", MANUAL_COMM);
		set_comm_timers ();
	}
	manualContact = _value;
}

- (IBAction)ok:(id)sender
{
	[[self window] makeFirstResponder:nil];			// End any active text field edits

	UpdateEndDates ();
	do_manual_comm_now ();

	[[self window] performClose:self];
}

@end
