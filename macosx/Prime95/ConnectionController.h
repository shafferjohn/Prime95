//
//  ConnectionController.h
//  Prime95
//
//  Created by George Woltman on 4/26/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface ConnectionController : NSWindowController {
	NSString *hostName;
	int	proxyEnabled;
	int	portNumber;
	NSString *userName;
	NSString *password;
	int	debug;
	char	szProxyHost[120], szProxyUser[50], szProxyPassword[50];
	unsigned short nProxyPort;
}

@property(readwrite, retain) NSString *hostName;
@property(readwrite, assign) int proxyEnabled;
@property(readwrite, assign) int portNumber;
@property(readwrite, retain) NSString *userName;
@property(readwrite, retain) NSString *password;
@property(readwrite, assign) int debug;

- (void)reInit;

@end
