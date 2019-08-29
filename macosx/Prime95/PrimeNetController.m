//
//  PrimeNetController.m
//  Prime95
//
//  Created by George Woltman on 4/26/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "PrimeNetController.h"
#import "AppController.h"
#import "ConnectionController.h"
#include "prime95.h"

@implementation PrimeNetController

- (id)init
{
	if (![super initWithWindowNibName:@"PrimeNet"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setUsePrimeNet:USE_PRIMENET];

	if (strcmp (USERID, "ANONYMOUS") == 0)
		[self setUserID:@""];
	else
		[self setUserID:[[NSString alloc] initWithFormat:@"%s", USERID]];

	[self setComputerName:[[NSString alloc] initWithFormat:@"%s", COMPID]];
}

/* Remember and clear the START_UP_IN_PROGRESS flag.  Why?  This routine */
/* is our only chance to capture close/cancel, but is also called if we hit apply. */
/* Thus, process the request as if this is a cancel and then undo the cancel in the ok method. */

- (void)windowWillClose:(NSNotification *)notification
{
	startupInProgress = STARTUP_IN_PROGRESS;
	STARTUP_IN_PROGRESS = 0;
}

@synthesize usePrimeNet;
@synthesize userID;
@synthesize computerName;

- (IBAction)connection:(id)sender
{
	if (!connectionController) connectionController = [[ConnectionController alloc] init];
	else [connectionController reInit];
	[connectionController showWindow:self];
}

- (IBAction)ok:(id)sender
{
	const char *m_userid, *m_compid;
	int	update_computer_info = FALSE;

// Annoyingly, we have to force any in-progress editing to end.
// We do this by changing the first responder.

// Save the current first responder

	id oldFirstResponder = [[self window] firstResponder];

// Gracefully end all editing in our window (from Erik Buck).
// This will cause the user's changes to be committed.
	
	if ([[self window] makeFirstResponder:[self window]]) {
		// All editing is now ended and delegate messages sent etc.
	} else {
		// For some reason the text object being edited will
		// not resign first responder status so force an
		/// end to editing anyway
		[[self window] endEditingFor:nil];
	}

// Process the dialog box contents

	m_userid = [userID UTF8String];
	m_compid = [computerName UTF8String];
	if (m_userid == NULL || m_userid[0] == 0)
		m_userid = "ANONYMOUS";
	if (m_compid == NULL)
		m_compid = "";

	if (strcmp (USERID, m_userid) != 0) {
		strcpy (USERID, m_userid);
		sanitizeString (USERID);
		IniWriteString (INI_FILE, "V5UserID", USERID);
		update_computer_info = TRUE;
	}
	if (strcmp (COMPID, m_compid) != 0) {
		strcpy (COMPID, m_compid);
		sanitizeString (COMPID);
		IniWriteString (LOCALINI_FILE, "ComputerID", COMPID);
		update_computer_info = TRUE;
	}
	if (!USE_PRIMENET && usePrimeNet) {
		USE_PRIMENET = 1;
		create_window (COMM_THREAD_NUM);
		base_title (COMM_THREAD_NUM, "Communication thread");
		if (!STARTUP_IN_PROGRESS) set_comm_timers ();
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);
		spoolExistingResultsFile ();
	} else if (USE_PRIMENET && !usePrimeNet) {
		USE_PRIMENET = 0;
		if (!STARTUP_IN_PROGRESS) set_comm_timers ();
	} else if (update_computer_info)
		spoolMessage (PRIMENET_UPDATE_COMPUTER_INFO, NULL);

	IniWriteInt (INI_FILE, "UsePrimenet", USE_PRIMENET);

/* For historical reasons, this dialog box also does a Test/Continue */
/* when you are using primenet */

	if (!STARTUP_IN_PROGRESS && USE_PRIMENET)
		LaunchWorkerThreads (ALL_WORKERS, FALSE);

// If we had a first responder before, restore it

	if (oldFirstResponder != nil) {
		[[self window] makeFirstResponder:oldFirstResponder];
	}

/* Close this window */

	[[self window] performClose:self];

/* Start next dialog in startup chain */

	if (startupInProgress && USE_PRIMENET) {
		STARTUP_IN_PROGRESS = 1;
		[myAppController optionsCPU:nil];
	}
}

- (IBAction)linkToMersenneOrg:(id)sender
{
	NSURL *url = [NSURL URLWithString:@"http://v5www.mersenne.org/update/"];
	[[NSWorkspace sharedWorkspace] openURL:url];
}

@end
