//
//  WindowController.m
//  Prime95
//
//  Created by George Woltman on 4/18/09.
//  Copyright 2009-2015 Mersenne Research, Inc. All rights reserved.
//

#import "WindowController.h"
#include "prime95.h"

gwmutex	VIEW_MUTEX;		/* Lock for accessing Views Array */
int	VIEW_MUTEX_INITIALIZED = 0;
#define MAX_VIEWS	(MAX_NUM_WORKER_THREADS+2)	/* Main_thread, comm_thread and each worker thread. */
WindowController *Views[MAX_VIEWS] = {0};
char	ThreadTitles[MAX_VIEWS][80] = {0};
double	currentFontSize = 0.0;


@implementation WindowController

- (id)init
{
	if (![super initWithWindowNibName:@"Window"]) return nil;
	baseTitle = nil;
	[self showWindow:nil];
	return self;
}

- (void)dealloc
{
	[[self window] performClose:nil];
	[baseTitle release];
	[super dealloc];
}

- (void)setBaseTitle:(NSString *)newBaseTitle
{

// Save the string for later use

	[newBaseTitle retain];
	[baseTitle release];
	baseTitle = newBaseTitle;
}

// Set the window title.  Cocoa is not thread-safe so setting the
// title must be done by the main thread.

- (void)setTitle:(NSString *)title
{
	NSString *combinedTitle;
	if ([baseTitle length] && [title length])
		combinedTitle = [NSString stringWithFormat:@"%@ - %@", baseTitle, title];
	else if ([title length])
		combinedTitle = title;
	else
		combinedTitle = baseTitle;
	[[self window] performSelectorOnMainThread:@selector(setTitle:) withObject:combinedTitle waitUntilDone:NO];
}

- (void)outputStr:(NSString *)textToAdd
{
	NSRange range;

// If length > 250000, truncate the view

	if ([[textView textStorage] length] >= 250000) {
		range.location = 0;
		range.length = 5000;
		[textView replaceCharactersInRange:range withString:@"..."];
	}

// Output the text and scroll it into view (replace dummy newline written by create_window)

	if ([[textView textStorage] length] == 1) {
		range.location = 0;
		range.length = 1;
	} else {
		range.location = [[textView textStorage] length];
		range.length = 0;
	}

	[textView replaceCharactersInRange:range withString:textToAdd];

	range.length = [textToAdd length];
	[textView scrollRangeToVisible:range];
}

// Change the font size

- (void)setFontSize:(int)newSize
{
	NSFont *existingFont = [textView font];
	[textView setFont:[[NSFontManager sharedFontManager] convertFont:existingFont toSize:newSize]];
}

// Get the font size from the text view, then check if there is a different one specified in
// the INI file.

- (void)getFontSize
{
	NSFont *existingFont = [textView font];
	currentFontSize = [existingFont pointSize];
	currentFontSize = IniGetInt (INI_FILE, "FontSize", (int) currentFontSize);
}

@end


/* Implement the OS-specific routines that deal with displaying data */


/* Internal routine to convert thread number into it's corresponding */
/* view number.  This routine handles merging windows together. */

int map_thread_num_to_view_num (
	int	thread_num)
{

/* Merge main window into first worker window if so requested */

	if (thread_num == MAIN_THREAD_NUM) {
		if (MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS) return (1);
		if (MERGE_WINDOWS & MERGE_MAIN_WINDOW) return (2);
	}

/* Merge comm window into first worker window if so requested */

	if (thread_num == COMM_THREAD_NUM && MERGE_WINDOWS & MERGE_COMM_WINDOW)
		return (2);

/* Merge all worker windows into first worker window if so requested */

	if (thread_num >= 0 && MERGE_WINDOWS & MERGE_WORKER_WINDOWS)
		return (2);

/* No merging.  Add 2 to map thread numbers (-2...32) into */
/* view numbers (0..34) */

	return (thread_num+2);
}

/* Internal routine to convert thread number into it's corresponding */
/* windowController */

WindowController *map_thread_num_to_view (
	int	thread_num)
{
	return (Views[map_thread_num_to_view_num (thread_num)]);
}

/* Create an output window for the thread -- unless we created */
/* one earlier and the user has not closed it */

void create_window (
	int	thread_num)
{
	int	viewNum;

// Init our view mutexes

	if (!VIEW_MUTEX_INITIALIZED) {
		gwmutex_init (&VIEW_MUTEX);
		VIEW_MUTEX_INITIALIZED = 1;
	}

// Create the window if it hasn't been created yet

	gwmutex_lock (&VIEW_MUTEX);
	viewNum = map_thread_num_to_view_num (thread_num);
	if (Views[viewNum] == NULL) {
		static int firstWindow = TRUE;
		static NSPoint point;

		Views[viewNum] = [[WindowController alloc] init];

		// Cascade windows (the auto-cascade feature doesn't seem to work)
		[Views[viewNum] setShouldCascadeWindows:NO];
		if (firstWindow) {
			NSRect windowFrame = [[Views[viewNum] window] frame];
			point = NSMakePoint(NSMinX(windowFrame), NSMaxY(windowFrame));
			firstWindow = FALSE;
		}
		point = [[Views[viewNum] window] cascadeTopLeftFromPoint:point];

		// Load saved window coordinates
		[Views[viewNum] setWindowFrameAutosaveName:[NSString stringWithFormat:@"Window%d", viewNum]];

		// Add some text (to create font) before getting or setting font size
		[Views[viewNum] outputStr:@"\n"];

		// Set initial font size (use the default or one from the INI file)
		if (currentFontSize == 0.0) [Views[viewNum] getFontSize];
		[Views[viewNum] setFontSize:currentFontSize];
	}
	gwmutex_unlock (&VIEW_MUTEX);
}

/* Destroy an MDI output window for the thread */

void destroy_window (
	int	thread_num)
{
	int	viewNum;

	gwmutex_lock (&VIEW_MUTEX);
	viewNum = map_thread_num_to_view_num (thread_num);
	if (Views[viewNum] != NULL) {
		[Views[viewNum] release];
		Views[viewNum] = NULL;
	}
	gwmutex_unlock (&VIEW_MUTEX);
}

/* Set the title prefix for this window - only called once */

void base_title (
	int	thread_num,
	const char *str)
{
	NSAutoreleasePool *pool;
	WindowController *view;

// Create an autorelease pool.  This routine is probably called from
// a worker thread where an autorelease pool is not in place.

	pool = [[NSAutoreleasePool alloc] init];

/* When merging windows, apply some arbitrary rules to decide which */
/* base title to use. */

	if (MERGE_WINDOWS & MERGE_MAIN_WINDOW &&
	    MERGE_WINDOWS & MERGE_COMM_WINDOW &&
	    (MERGE_WINDOWS & MERGE_WORKER_WINDOWS || NUM_WORKER_THREADS == 1))
		str = "";
	else if ((thread_num == MAIN_THREAD_NUM &&
		  MERGE_WINDOWS & MERGE_MAIN_WINDOW) ||
		 (thread_num == COMM_THREAD_NUM &&
		  MERGE_WINDOWS & MERGE_COMM_WINDOW) ||
		 (thread_num >= 0 &&
		  MERGE_WINDOWS & MERGE_WORKER_WINDOWS)) {
		if (MERGE_WINDOWS & MERGE_WORKER_WINDOWS &&
		    NUM_WORKER_THREADS > 1)
			str = "Workers";
		else
			str = "Worker";
	}

/* Pass the base title on to the view */

	gwmutex_lock (&VIEW_MUTEX);
	view = map_thread_num_to_view (thread_num);
	if (view != NULL) [view setBaseTitle:[NSString stringWithFormat:@"%s", str]];
	gwmutex_unlock (&VIEW_MUTEX);

// Free resources

	[pool drain];
}

/* Put a title on the window */

void title (
	int	thread_num,
	const char *str)
{
	NSAutoreleasePool *pool;
	WindowController *view;
	char	merged_title[160];

	if (thread_num == COMM_THREAD_NUM && MERGE_WINDOWS & MERGE_COMM_WINDOW)
		return;

// Create an autorelease pool.  This routine is probably called from
// a worker thread where an autorelease pool is not in place.

	pool = [[NSAutoreleasePool alloc] init];

/* When merging windows, apply some arbitrary rules to decide which */
/* what the title should be. */

	strcpy (ThreadTitles[thread_num+2], str);
	if (thread_num == MAIN_THREAD_NUM && MERGE_WINDOWS & MERGE_MAIN_WINDOW)
		thread_num = 0;
	if (thread_num >= 0 && MERGE_WINDOWS & MERGE_WORKER_WINDOWS) {
		int	i;
		merged_title[0] = 0;
		for (i = 2; i < MAX_VIEWS; i++) {
			if (merged_title[0] && ThreadTitles[i][0])
				strcat (merged_title, ",");
			strcat (merged_title, ThreadTitles[i]);
			if (strlen (merged_title) > 80) {
				merged_title[77] = 0;
				strcat (merged_title, "...");
				break;
			}
		}
		str = merged_title[0] ? merged_title : "Not running";
	}

/* Pass the base title on to the view */

	gwmutex_lock (&VIEW_MUTEX);
	view = map_thread_num_to_view (thread_num);
	if (view != NULL) [view setTitle:[NSString stringWithFormat:@"%s", str]];
	gwmutex_unlock (&VIEW_MUTEX);

// Free resources

	[pool drain];
}

void flashWindowAndBeep (void)
{
	printf ("\007");
}

void RealOutputStr (int thread_num, const char *str)
{
static	int	partial_line_output[MAX_VIEWS] = {FALSE};
	NSAutoreleasePool *pool;
	WindowController *view;
	NSString *textToAdd;

/* My latest MacBook Pro running OS X Mavericks changes Prime95's priority to near zero. */
/* This wouldn't be a problem except that the scheduler feels free to invoke SpeedStep */
/* when low priority processes are using 100% of the CPU time.  This hack changes the */
/* priority frequently to "trick" the scheduler into keeping the CPU running full speed. */
/* Alas, this trick does not seem to work.  ode left here in case it might be useful one day. */

	if (thread_num != MAIN_THREAD_NUM &&
	    thread_num != COMM_THREAD_NUM &&
	    IniGetInt (INI_FILE, "FrequentSetPriority", 0)) {
		setOsThreadPriority (PRIORITY);
	}

// Create an autorelease pool.  This routine is probably called from
// a worker thread where an autorelease pool is not in place.

	pool = [[NSAutoreleasePool alloc] init];

// Find the view (window) to output to

	gwmutex_lock (&VIEW_MUTEX);
	view = map_thread_num_to_view (thread_num);

// Shouldn't happen, but catch it just in case

	if (view == NULL) goto done;

// When merging windows output a prefix so user knows which thread is
// responsible for this line of output

	if (!(MERGE_WINDOWS & MERGE_NO_PREFIX) &&
	    !partial_line_output[thread_num+2] &&
	    ((thread_num == MAIN_THREAD_NUM &&
	      MERGE_WINDOWS & MERGE_MAIN_WINDOW) ||
	     (thread_num == MAIN_THREAD_NUM &&
	      MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS) ||
	     (thread_num == COMM_THREAD_NUM &&
	      MERGE_WINDOWS & MERGE_COMM_WINDOW) ||
	     (thread_num == COMM_THREAD_NUM &&
	      MERGE_WINDOWS & MERGE_MAINCOMM_WINDOWS) ||
	     (thread_num == 0 &&
	      MERGE_WINDOWS & (MERGE_MAIN_WINDOW | MERGE_COMM_WINDOW)) ||
	     (thread_num == 0 &&
	      MERGE_WINDOWS & MERGE_WORKER_WINDOWS &&
	      NUM_WORKER_THREADS > 1) ||
	     (thread_num >= 1 &&
	      MERGE_WINDOWS & MERGE_WORKER_WINDOWS))) {
		char	prefix[50];
		if (thread_num == MAIN_THREAD_NUM)
			strcpy (prefix, "[Main thread");
		else if (thread_num == COMM_THREAD_NUM)
			strcpy (prefix, "[Comm thread");
		else if (! (MERGE_WINDOWS & MERGE_WORKER_WINDOWS) ||
			 NUM_WORKER_THREADS == 1)
			strcpy (prefix, "[Work thread");
		else
			sprintf (prefix, "[Worker #%d", thread_num+1);
		if (str[0] == '[') {
			strcat (prefix, " ");
			str++;
		} else
			strcat (prefix, "] ");
		textToAdd = [NSString stringWithFormat:@"%s%s", prefix, str];
	}

// No prefix is required.  Simply output the text to the view.

	else
		textToAdd = [NSString stringWithFormat:@"%s", str];

// Output the text

	[view performSelectorOnMainThread:@selector(outputStr:) withObject:textToAdd waitUntilDone:NO];

// Remember if we are in the middle of outputting a line

	partial_line_output[thread_num+2] = (str[strlen(str)-1] != '\n');

// Free resources and return

done:	gwmutex_unlock (&VIEW_MUTEX);
	[pool drain];
}

void BlinkIcon (int thread_num, int x)
{
}

void ChangeIcon (int thread_num, int x)
{
}

void BiggerFonts ()
{
	int	i;

	gwmutex_lock (&VIEW_MUTEX);
	currentFontSize *= 1.1;
	IniWriteInt (INI_FILE, "FontSize", (int) currentFontSize);
	for (i = 0; i < MAX_VIEWS; i++) {
		if (Views[i] != NULL) [Views[i] setFontSize:currentFontSize];
	}
	gwmutex_unlock (&VIEW_MUTEX);
}

void SmallerFonts ()
{
	int	i;

	gwmutex_lock (&VIEW_MUTEX);
	currentFontSize /= 1.1;
	IniWriteInt (INI_FILE, "FontSize", (int) currentFontSize);
	for (i = 0; i < MAX_VIEWS; i++) {
		if (Views[i] != NULL) [Views[i] setFontSize:currentFontSize];
	}
	gwmutex_unlock (&VIEW_MUTEX);
}
