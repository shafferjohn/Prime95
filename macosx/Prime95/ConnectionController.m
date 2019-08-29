//
//  ConnectionController.m
//  Prime95
//
//  Created by George Woltman on 4/26/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "ConnectionController.h"
#include "prime95.h"
void getProxyInfo (char *, unsigned short *, char *, char *);

@implementation ConnectionController

- (id)init
{
	if (![super initWithWindowNibName:@"Connection"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	int	primenet_debug;

	nProxyPort = 8080;
	getProxyInfo (szProxyHost, &nProxyPort, szProxyUser, szProxyPassword);
	primenet_debug = IniSectionGetInt (INI_FILE, "PrimeNet", "Debug", 0);

	[self setHostName:[[NSString alloc] initWithFormat:@"%s", szProxyHost]];
	[self setProxyEnabled:(szProxyHost[0] != 0)];
	[self setPortNumber:nProxyPort];
	[self setUserName:[[NSString alloc] initWithFormat:@"%s", szProxyUser]];
	[self setPassword:[[NSString alloc] initWithFormat:@"%s", szProxyPassword]];
	[self setDebug:primenet_debug];
}

- (NSString *)hostName
{
	return hostName;
}

- (void)setHostName:(NSString *) _value
{
	const char *new_host_name;

	new_host_name = [_value UTF8String];
	if (new_host_name == NULL) new_host_name = "";

	if (strcmp (szProxyHost, new_host_name) != 0) {
		char	iniProxyHost[120];

		strcpy (szProxyHost, new_host_name);
		[self setProxyEnabled:(szProxyHost[0] != 0)];

		strcpy (iniProxyHost, szProxyHost);
		if (iniProxyHost[0] && nProxyPort != 8080)
			sprintf (iniProxyHost + strlen (iniProxyHost), ":%d", nProxyPort);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyHost", iniProxyHost);
	}

	[_value retain];
	[hostName release];
	hostName = _value;
}

- (int)portNumber
{
	return portNumber;
}

- (void)setPortNumber:(int) _value
{
	if (nProxyPort != _value) {
		char	iniProxyHost[120];

		nProxyPort = _value;

		strcpy (iniProxyHost, szProxyHost);
		if (iniProxyHost[0] && nProxyPort != 8080)
			sprintf (iniProxyHost + strlen (iniProxyHost), ":%d", nProxyPort);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyHost", iniProxyHost);
	}

	portNumber = _value;
}

- (NSString *)userName
{
	return userName;
}

- (void)setUserName:(NSString *) _value
{
	const char *new_user_name;

	new_user_name = [_value UTF8String];
	if (new_user_name == NULL) new_user_name = "";

	if (strcmp (szProxyUser, new_user_name) != 0) {
		strcpy (szProxyUser, new_user_name);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyUser", szProxyUser);
	}

	[_value retain];
	[userName release];
	userName = _value;
}

- (NSString *)password
{
	return password;
}

- (void)setPassword:(NSString *) _value
{
	const char *new_password;

	new_password = [_value UTF8String];
	if (new_password == NULL) new_password = "";

	if (strcmp (szProxyPassword, new_password) != 0) {
		strcpy (szProxyPassword, new_password);
		IniSectionWriteString (INI_FILE, "PrimeNet", "ProxyPass", szProxyPassword);
		IniSectionWriteInt (INI_FILE, "PrimeNet", "ProxyMask", 0);
	}

	[_value retain];
	[password release];
	password = _value;
}

- (int)debug
{
	return debug;
}

- (void)setDebug:(int) _value
{
	if (IniSectionGetInt (INI_FILE, "PrimeNet", "Debug", 0) != _value) {
		IniSectionWriteInt (INI_FILE, "PrimeNet", "Debug", _value);
	}

	debug = _value;
}

@synthesize proxyEnabled;

@end
