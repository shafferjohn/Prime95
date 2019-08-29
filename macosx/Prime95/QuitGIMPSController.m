//
//  QuitGIMPSController.m
//  Prime95
//
//  Created by George Woltman on 5/1/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "QuitGIMPSController.h"
#include "prime95.h"

@implementation QuitGIMPSController

- (id)init
{
	if (![super initWithWindowNibName:@"QuitGIMPS"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

#define MAC_MANUAL_QUIT	"You have elected to remove this computer from the Great Internet Mersenne Prime Search.  Other computers using this user ID will not be affected.\n\nPlease send the file results.txt to woltman@alum.mit.edu.\n\nAre you sure you want to do this?"
#define MAC_PRIMENET_QUIT "You have elected to remove this computer from the Great Internet Mersenne Prime Search.  Other computers using this user ID will not be affected.\n\nPlease make sure your results have been successfully sent to the server before uninstalling the program. If in doubt, you can send the results.txt file to woltman@alum.mit.edu.\n\nYou can either complete your current assignment or you can quit GIMPS immediately."

- (void)reInit
{
	[self setExplanation:[[NSString alloc] initWithFormat:@"%s", !USE_PRIMENET ? MAC_MANUAL_QUIT : MAC_PRIMENET_QUIT]];
}

@synthesize explanation;

- (IBAction)quitLater:(id)sender
{
	OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS after current work completes.\n");
	IniWriteInt (INI_FILE, "NoMoreWork", 1);
	[[self window] performClose:self];
}

- (IBAction)quitNow:(id)sender
{
	OutputBoth (MAIN_THREAD_NUM, "Quitting GIMPS immediately.\n");
	stop_workers_for_escape ();
	if (USE_PRIMENET) spoolMessage (MSG_QUIT_GIMPS, NULL);
	[[self window] performClose:self];
}

@end
