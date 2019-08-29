/*
 * Primenet communication routines for all operating systems
 * Uses sockets and HTTP
 */ 

/*
// THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A
// PARTICULAR PURPOSE.
//
// Copyright (c) 1997-2019 Mersenne Research, Inc. All Rights Reserved.
//
//  MODULE:   primenet.c
//
//  PURPOSE:  Implements PrimeNet Version 4 and 5 API as HTTP network client
//
//  AUTHOR:   Peter Hunter, on the basis of work by Scott Kurowski (v3 API)
//            Michiel van Loon, OS/2 adaptations 
//            Kurowski 5/1998, 4.0 API support for MPrime 16.x
//            Kurowski 9/1999, 4.0 API changes for MPrime 19.x
//	      Woltman 1/2002, Windows support and bug fixes
//	      Woltman 10/2005, Version 5 API support, CURL library
//	      Woltman 9/2017, Added interim residues to AP msg, added PRP support
//	      Woltman 9/2019, Added PRP dblchk support
//
//  ASSUMPTIONS: 1. less than 4k of data is sent or received per call
//               2. HTTP/1.1
//               3. PrimeNet Version 5 or later API on server and client
*/

/* Linux defines, adapted for OS/2, FreeBSD, and Windows */

#define CURL_STATICLIB

#ifdef __WATCOMC__
#include <types.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>

#include "secure5.c"

#ifdef _WINDOWS_
#include <curl.h>
#include <winsock.h>
#else
#include <curl/curl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <errno.h>
#include <arpa/inet.h>
#ifndef __IBMC__
#include <unistd.h>
#endif
typedef int SOCKET;
#endif

#if defined (__EMX__) && ! defined (AOUT)
#define htonl(x) (lswap(x))
#define ntohl(x) (lswap(x))
#define htons(x) (bswap(x))
#define ntohs(x) (bswap(x))
unsigned short bswap(unsigned short);
#endif

static const char iniSection[] = "PrimeNet";
static const char hx[] = "0123456789ABCDEF";
static const char szSITEstr[] = "v5.mersenne.org";	/* PrimeNet Server's home domain */
#define nHostPort 80					/* Internet PrimeNet port */
static const char szFILE[] = "/v5server/?";		/* HTTP GET string */

#define PROXY_HOST_BUFSIZE	120
#define PROXY_USER_BUFSIZE	50
#define PROXY_PASSWORD_BUFSIZE	50

/* implement the missing inet_aton call */

#if defined (_WINDOWS_) || defined (__EMX__)
int inet_aton (char *cp, struct in_addr *inp)
{
	u_long temp;
 
	temp = inet_addr(cp);
	if (temp == -1) return (0);
	inp->s_addr = temp;
	return (1);
}
#endif

/* Routine to get the error number after a failed sockets call */

int getLastSocketError (void)
{
#ifdef _WINDOWS_
	return (WSAGetLastError ());
#else
	return (errno);
#endif
}

/* simple password de/scrambler */

char SCRAMBLE_STRING[] = "/cgi-bin/pnHttp.exe";

void scramble (char *s)
{
	char	out[100];
	char	*p = s, *z = out;
	unsigned int i, c = (unsigned int) strlen (SCRAMBLE_STRING);

	for (i = 0; i < strlen (s); i++) {
		int b = (unsigned char) *p++ ^ SCRAMBLE_STRING[i % c];
		*z++ = hx[b >> 4];
		*z++ = hx[b % 16];
	}
	*z = 0;
	strcpy (s, out);
}

void unscramble (char *s)
{
	char	out[50];
	char	*q = s, *z = out;
	unsigned int i, c = (unsigned int) strlen (SCRAMBLE_STRING);

	for (i = 0; i < strlen (s) >> 1; i++) {
		*z = (char) (strchr (hx, *q++) - hx) * 16;
		*z += (char) (strchr (hx, *q++) - hx);
		*z++ ^= SCRAMBLE_STRING[i % c];
	}
	*z = 0;
	strcpy (s, out);
}

/* base64 encode for basic proxy passwords */

static const char encode[] = {
  'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
  'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
  'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
  'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
  'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
  'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
  'w', 'x', 'y', 'z', '0', '1', '2', '3',
  '4', '5', '6', '7', '8', '9', '+', '/'
};

/* Base64-encode a null-terminated string. */

void encode64 (
	char	*data)
{
	char	*s, *end, *buf;
	unsigned int x;
	size_t	length;
	int	i, j;
	char	temp[400];

	length = strlen (data);
	if (length == 0) return;

	end = data + length - 3;

	buf = temp;
	for (s = data; s < end; ) {
		x = *s++ << 24;
		x |= *s++ << 16;
		x |= *s++ << 8;

		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
		x <<= 6;
		*buf++ = encode[x >> 26];
	}
	end += 3;

	x = 0;
	for (i = 0; s < end; i++)
		x |= *s++ << (24 - 8 * i);

	for (j = 0; j < 4; j++) {
		if (8 * i >= 6 * j) {
			*buf++ = encode [x >> 26];
			x <<= 6;
		} else {
			*buf++ = '=';
		}
	}

	*buf = 0;

	strcpy (data, temp);
}


/* Get proxy information from INI file */

void getProxyInfo (
	char	*szProxyHost,
	unsigned short *nProxyPort,
	char	*szProxyUser,
	char	*szProxyPassword)
{
	char	*colon;

/* Initialize return variables in case we return early */

	*nProxyPort = 8080;
	*szProxyUser = 0;
	*szProxyPassword = 0;

/* Get the host name of the optional proxy server.  If using a proxy */
/* server strip the optional http:// prefix. */

	IniSectionGetString (INI_FILE, iniSection, "ProxyHost", szProxyHost, PROXY_HOST_BUFSIZE, NULL);
	if (szProxyHost[0] == 0) return;

	if ((szProxyHost[0] == 'H' || szProxyHost[0] == 'h') &&
	    (szProxyHost[1] == 'T' || szProxyHost[1] == 't') &&
	    (szProxyHost[2] == 'T' || szProxyHost[2] == 't') &&
	    (szProxyHost[3] == 'P' || szProxyHost[3] == 'p') &&
	    szProxyHost[4] == ':' && szProxyHost[5] == '/' &&
	    szProxyHost[6] == '/')
		safe_strcpy (szProxyHost, szProxyHost + 7);

/* Get optional port number */

	if ((colon = strchr (szProxyHost, ':'))) {
		*nProxyPort = (unsigned short) atoi (colon + 1);
		*colon = 0;
	} else
		*nProxyPort = 8080;

/* Secure proxy - get username and password to negotiate access */

	IniSectionGetString (INI_FILE, iniSection, "ProxyUser", szProxyUser, PROXY_USER_BUFSIZE, NULL);
	IniSectionGetString (INI_FILE, iniSection, "ProxyPass", szProxyPassword, PROXY_PASSWORD_BUFSIZE, NULL);

/* Scramble or unscramble the password as necessary */

	if (!IniSectionGetInt (INI_FILE, iniSection, "ProxyMask", 0)) {
		scramble (szProxyPassword);
		IniSectionWriteString (INI_FILE, iniSection, "ProxyPass", szProxyPassword);
		IniSectionWriteInt (INI_FILE, iniSection, "ProxyMask", 1);
	}
	unscramble (szProxyPassword);
}

/*///////////////////////////////////////////////////////////////////////////
//
// HTTP GET procedure (Sockets Implementation)
//
///////////////////////////////////////////////////////////////////////////*/

/*
// pnHttpServer: Uses GET to send a formatted HTTP argument string
//               and downloads the server result page
*/

int pnHttpServer (char *pbuf, unsigned cbuf, char* postargs)
{
	char	szBuffer[1000];
	char	szURL[4096];			/* URL assembly buffer */
	char	buf[4096];
	struct in_addr defaddr;
	struct hostent *hp, def;
	struct sockaddr_in sn;
	SOCKET	s;
	int	res, debug, url_format;
	int	timeout;
	char	*alist[1];
	unsigned int count;
	char	szProxyHost[PROXY_HOST_BUFSIZE];
	char	szProxyUser[PROXY_USER_BUFSIZE];
	char	szProxyPassword[PROXY_PASSWORD_BUFSIZE];
	char	szSITE[80];
	char	*con_host;
	char	szOtherGetInfo[256];
	unsigned short nProxyPort, con_port;

/* Get debug logging and URL format flags */

	debug = IniSectionGetInt (INI_FILE, iniSection, "Debug", 0);
	url_format = IniSectionGetInt (INI_FILE, iniSection, "UseFullURL", 2);
 
/* Use MersenneIP setting so that SeventeenOrBust can use the */
/* UseCurl=0 option. */

	IniSectionGetString (INI_FILE, iniSection, "MersenneIP", szSITE, sizeof (szSITE), szSITEstr);

/* Get information about the optional proxy server */

	getProxyInfo (szProxyHost, &nProxyPort, szProxyUser, szProxyPassword);
	if (szProxyHost[0]) {
		con_host = szProxyHost;
		con_port = nProxyPort;
	}

/* No proxy server - use default site (mersenne.org) */

	else {
		con_host = szSITE;
		con_port = nHostPort;
	}

/* Output debug info */

redirect:
	if (debug) {
		sprintf (buf, "host = %s, port = %d\n", con_host, con_port);
		LogMsg (buf);
	}

/* Convert host name into an IP address */

	hp = gethostbyname (con_host);
	if (!hp) {
		if (!inet_aton (con_host, &defaddr)) {
			sprintf (buf, "Unknown host: %s\n", con_host);
			if (debug) LogMsg (buf);
			else OutputStr (COMM_THREAD_NUM, buf);
			return (PRIMENET_ERROR_CONNECT_FAILED);
		}
		alist[0] = (char *) &defaddr;
		def.h_name = szSITE;
		def.h_addr_list = alist;
		def.h_length = sizeof (struct in_addr);
		def.h_addrtype = AF_INET;
		def.h_aliases = 0;
		hp = &def;
	}

	if (debug) {
		sprintf (buf, "IP-addr = %s\n", inet_ntoa (*(struct in_addr *)hp->h_addr));
		LogMsg (buf);
	}

	memset (&sn, 0, sizeof (sn));
	sn.sin_family = hp->h_addrtype;
	if (hp->h_length > (int) sizeof (sn.sin_addr)) {
		hp->h_length = sizeof (sn.sin_addr);
	}
	memcpy (&sn.sin_addr, hp->h_addr, hp->h_length);
	sn.sin_port = htons (con_port);

/* Create a socket and connect to server */

rel_url:
	if ((s = socket (hp->h_addrtype, SOCK_STREAM, 0)) < 0) {
		sprintf (buf, "Error in socket call: %d\n", getLastSocketError ());
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		return (PRIMENET_ERROR_CONNECT_FAILED);
	}

	if (connect (s, (struct sockaddr *) &sn, sizeof (sn)) < 0) {
		sprintf (buf, "Error in connect call: %d\n", getLastSocketError ());
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		closesocket (s);
		return (PRIMENET_ERROR_CONNECT_FAILED);
	}

/* Prevent SIGPIPE signals in OS X, BSD and other environments */

#ifdef SO_NOSIGPIPE
	{
		int	i = 1;
		res = setsockopt (s, SOL_SOCKET, SO_NOSIGPIPE, &i, sizeof(i));
		if (res < 0) {
			if (debug) {
				sprintf (buf, "Error in NOSIGPIPE call: %d\n", getLastSocketError ());
				LogMsg (buf);
			}
		}
	}
#endif

/* GET method, data follows ? in URL */

	strcpy (szURL, "GET ");
	if (*szProxyHost || url_format == 1) {
		strcat (szURL, "http://");
		strcat (szURL, szSITE);
		if (IniSectionGetInt (INI_FILE, iniSection, "SendPortNumber", 1))
			sprintf (szURL + strlen (szURL), ":%d", nHostPort);
	}
	strcat (szURL, szFILE);
	strcat (szURL, postargs);
	strcat (szURL, " HTTP/1.1\r\n");

/* Append host: header */

	if (!*szProxyHost)
		sprintf (szURL+strlen(szURL), "Host: %s\r\n", szSITE);

/* Persistent connections are not supported */

	strcat (szURL, "Connection: close\r\n");

/* Append other GET info */

	IniSectionGetString (INI_FILE, iniSection, "OtherGetInfo", szOtherGetInfo,
		sizeof (szOtherGetInfo), NULL);
	if (*szOtherGetInfo) {
		strcat (szURL, szOtherGetInfo);
		strcat (szURL, "\r\n");
	}

/* Append proxy authorization here */

	if (*szProxyHost && *szProxyUser) {
		char	buf[200];
		strcat (szURL, "Proxy-Authorization: Basic ");
		sprintf (buf, "%s:%s", szProxyUser, szProxyPassword);
		encode64 (buf);
		strcat (szURL, buf);
		strcat (szURL, "\r\n");
	}

	strcat (szURL, "\r\n");
	if (debug) LogMsg (szURL);

/* Send the URL request */

	timeout = 90000;		/* 90 seconds */
	res = setsockopt (s, SOL_SOCKET, SO_SNDTIMEO, (char *) &timeout, sizeof (timeout));
	if (res < 0) {
		if (debug) {
			sprintf (buf, "Error in send timeout call: %d\n", getLastSocketError ());
			LogMsg (buf);
		}
	}
	res = send (s, szURL, (int) strlen (szURL), 0
#ifdef MSG_NOSIGNAL /* Prevent SIGPIPE in Linux (and others) */
		    | MSG_NOSIGNAL
#endif
		    );
	if (res < 0) {
		sprintf (buf, "Error in send call: %d\n", getLastSocketError ());
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		closesocket (s);
		if (url_format == 2 && *szProxyHost == 0) {
			if (debug) LogMsg ("Trying full URL\n");
			url_format = 1;
			goto rel_url;
		}
		return (PRIMENET_ERROR_SEND_FAILED);
	}

/* Now accumulate the response */

	timeout = 90000;		/* 90 seconds */
	res = setsockopt (s, SOL_SOCKET, SO_RCVTIMEO, (char *) &timeout, sizeof (timeout));
	if (res < 0) {
		if (debug) {
			sprintf (buf, "Error in receive timeout call: %d\n", getLastSocketError ());
			LogMsg (buf);
		}
	}
	*pbuf = 0; count = 1;
	while (count < cbuf) {
		res = recv (s, szBuffer, 999, 0);
		if (res < 0) {
			sprintf (buf, "Error in recv call: %d\n", getLastSocketError ());
			if (debug) LogMsg (buf);
			else OutputStr (COMM_THREAD_NUM, buf);
			closesocket (s);
			if (url_format == 2 && *szProxyHost == 0) {
				if (debug) LogMsg ("Trying full URL\n");
				url_format = 1;
				goto rel_url;
			}
			return (PRIMENET_ERROR_RECV_FAILED);
		}
		if (res == 0) break;
		szBuffer[res] = 0;
		if (debug) {
			sprintf (buf, "RECV: %s\n", szBuffer);
			LogMsg (buf);
		}
		if ((count + strlen(szBuffer)) >= cbuf)
			szBuffer[cbuf - count] = 0;
		strcat (pbuf, szBuffer);
		count += (unsigned int) strlen (szBuffer);
	}

	closesocket (s);

/* pbuf + 9 is where the message code following HTTP/1.1 starts */

	if (count <= 10) res = -1;
	else res = atoi (pbuf + 9);

/* Some proxy servers can redirect us to another host.  These are */
/* the 300 series of error codes.  This code probably isn't right. */
/* We can improve it as people find problems. */

	if (res >= 300 && res <= 399) {
		char	*location, *colon;
		location = strstr (pbuf, "Location:");
		if (location != NULL) {
			if (debug) LogMsg ("Attempting a redirect.\n");
			location += 9;
			while (isspace (*location)) location++;
			strcpy (szProxyHost, location);

/* Parse the redirection address */

			if ((location[0] == 'H' || location[0] == 'h') &&
			    (location[1] == 'T' || location[1] == 't') &&
			    (location[2] == 'T' || location[2] == 't') &&
			    (location[3] == 'P' || location[3] == 'p') &&
			    location[4] == ':' && location[5] == '/' &&
			    location[6] == '/')
				safe_strcpy (location, location + 7);
			con_host = location;

/* Get optional port number */

			if ((colon = strchr (location, ':')) != NULL) {
				con_port = (unsigned short) atoi (colon + 1);
				*colon = 0;
			} else
				con_port = 80;
			goto redirect;
		}
	}

/* Any return code other than 200 is an error.  We've had problems using */
/* both full URLs and relative URLs.  Thus, our default behavior is to try */
/* both before giving up. */

	if (res != 200) {
		if (url_format == 2 && *szProxyHost == 0) {
			if (debug) LogMsg ("Trying full URL\n");
			url_format = 1;
			goto rel_url;
		}
		strcpy (buf, "HTTP return code is not 200\n");
		if (debug) LogMsg (buf);
		else OutputStr (COMM_THREAD_NUM, buf);
		return (PRIMENET_ERROR_SERVER_UNSPEC);
	}

	return (PRIMENET_NO_ERROR);
}


/*///////////////////////////////////////////////////////////////////////////
//
// HTTP GET procedure (cURL Implementation)
//
///////////////////////////////////////////////////////////////////////////*/

/* This callback routine assembles the server's response to our request */

size_t WriteMemoryCallback (
	void	*ptr,
	size_t	size,
	size_t	nmemb,
	void	*data)
{
	size_t realsize = size * nmemb;
	char	*buf = (char *) data;
	int	buflen = (int) strlen (buf);

/* Truncate response to fit in a 4096 byte buffer */

	if (buflen + realsize <= 4095) {
		memcpy (buf + buflen, ptr, realsize);
		buf[buflen + realsize] = 0;
	} else {
		memcpy (buf + buflen, ptr, 4095 - buflen);
		buf[4095] = 0;
	}
	return realsize;
}

/* This callback is for dumping out cURL debug information */

int curl_trace (
	CURL	*handle,
	curl_infotype type,
	char	*data,
	size_t	size,
	void	*userp)
{
	char	*text;
	char	buf[4096];
	int	len;

	switch (type) {
	case CURLINFO_TEXT:
		text = "== Info";
		break;
	case CURLINFO_HEADER_OUT:
		text = "=> Send header";
		break;
	case CURLINFO_DATA_OUT:
		text = "=> Send data";
		break;
	case CURLINFO_HEADER_IN:
		text = "<= Recv header";
		break;
	case CURLINFO_DATA_IN:
		text = "<= Recv data";
		break;
	case CURLINFO_SSL_DATA_IN:
		text = "<= Recv SSL data";
		break;
	case CURLINFO_SSL_DATA_OUT:
		text = "<= Send SSL data";
		break;
	default: /* in case a new one is introduced to shock us */
		return 0;
	}

/* Output the data. */

	strcpy (buf, text);
	strcat (buf, ": ");
	len = (int) strlen (buf);
	memcpy (buf + len, data, size);
	if (data[size-1] != '\n') {
		buf[len + size] = '\n';
		buf[len + size + 1] = 0;
	} else
		buf[len + size] = 0;
	LogMsg (buf);

	return 0;
}

/*
// pnHttpServerCURL: Uses GET to send a formatted HTTP argument string
//               and downloads the server result page
*/

int pnHttpServerCURL (char *pbuf, unsigned cbuf, char* postargs)
{
	CURL	*curl;
	CURLcode res;
	int	try_proxy, debug;
	char	szSITE[120];
	char	url[4096], buf[4150], errbuf[CURL_ERROR_SIZE];
	char	szProxyHost[PROXY_HOST_BUFSIZE];
	char	szProxyUser[PROXY_USER_BUFSIZE];
	char	szProxyPassword[PROXY_PASSWORD_BUFSIZE];
	unsigned short nProxyPort;

/* Get debug logging flag */

	debug = IniSectionGetInt (INI_FILE, iniSection, "Debug", 0);
 
/* Loop to try with proxy, then after a failure without proxy.  Ixfd64 requested this */
/* feature here:  https://www.mersenneforum.org/showpost.php?p=505557&postcount=415 */

	for (try_proxy = 1; ; try_proxy = 0) {

/* Init the cURL structures */

		curl = curl_easy_init ();
		if (curl == NULL) return (PRIMENET_ERROR_CURL_INIT);
		curl_easy_setopt (curl, CURLOPT_NOPROGRESS, 1);

/* Give curl library the HTTP string to send */

		strcpy (url, "http://");
		IniSectionGetString (INI_FILE, iniSection, "MersenneIP", szSITE, sizeof (szSITE), szSITEstr);
		strcat (url, szSITE);
		if (IniSectionGetInt (INI_FILE, iniSection, "SendPortNumber", 0))
			sprintf (url + strlen (url), ":%d", nHostPort);
		strcat (url, szFILE);
		strcat (url, postargs);
		curl_easy_setopt (curl, CURLOPT_URL, url);
		if (debug) {
			sprintf (buf, "URL: %s\n", url);
			LogMsg (buf);
		}

/* Get information about the optional proxy server */

		if (try_proxy) {
			getProxyInfo (szProxyHost, &nProxyPort, szProxyUser, szProxyPassword);
			if (szProxyHost[0]) {
				curl_easy_setopt (curl, CURLOPT_PROXY, szProxyHost);
				curl_easy_setopt (curl, CURLOPT_PROXYPORT, (long) nProxyPort);
//bug?				curl_easy_setopt (curl, CURLOPT_PROXYTYPE, ???);
				if (szProxyUser[0]) {
					sprintf (buf, "%s:%s", szProxyUser, szProxyPassword);
					curl_easy_setopt (curl, CURLOPT_PROXYUSERPWD, buf);
					curl_easy_setopt (curl, CURLOPT_PROXYAUTH, CURLAUTH_ANY);
				}
			}
		}

/* Setup function to receive the response */

		curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
		curl_easy_setopt (curl, CURLOPT_WRITEDATA, (void *) pbuf);
		pbuf[0] = 0;

/* Output verbose debug information */

		if (debug >= 2) {
			curl_easy_setopt (curl, CURLOPT_DEBUGFUNCTION, curl_trace);
			curl_easy_setopt (curl, CURLOPT_DEBUGDATA, NULL);
			/* the DEBUGFUNCTION has no effect until we enable VERBOSE */
			curl_easy_setopt (curl, CURLOPT_VERBOSE, 1);
		}

/* Send the URL request */

		curl_easy_setopt (curl, CURLOPT_FOLLOWLOCATION, 1);
		curl_easy_setopt (curl, CURLOPT_NOSIGNAL, 1);
		curl_easy_setopt (curl, CURLOPT_CONNECTTIMEOUT, 180);
		curl_easy_setopt (curl, CURLOPT_TIMEOUT, 180);
		curl_easy_setopt (curl, CURLOPT_ERRORBUFFER, errbuf);
		res = curl_easy_perform (curl);
		if (res != CURLE_OK) {
			sprintf (buf, "CURL library error: %s\n", errbuf);
			LogMsg (buf);
			OutputStr (COMM_THREAD_NUM, buf);
			curl_easy_cleanup (curl);
			// By default, try again without using a proxy server
			if (try_proxy && szProxyHost[0] && IniSectionGetInt (INI_FILE, iniSection, "TryNoProxyAfterProxyFailure", 1))
				continue;
			return (PRIMENET_ERROR_CURL_PERFORM);
		}

/* Success, break out of loop */

		break;
	}

/* Cleanup */

	curl_easy_cleanup (curl);

/* Log a debug message */

	if (debug) {
		LogMsg ("RESPONSE:\n");
		LogMsg (pbuf);
	}

/* Return success */

	return (PRIMENET_NO_ERROR);
}



/*///////////////////////////////////////////////////////////////////////////
//
// HTTP GET argument formatting procedures
//
////////////////////////////////////////////////////////////////////////////*/

/* armor parameter control chars as hex codes for transport */

#define ARMOR_CHARS		"&+%\r\n"

char *armor (char *d, char *s)
{

/* & is token delimiter, '+' is space char */

	while (*s) {
		if (strchr (ARMOR_CHARS, *s)) {	
			*d++ = '%';	/* convert chars to %nn hex codes */
			*d++ = hx[(*s) / 16];
			*d++ = hx[(*s) % 16];
		} else if (*s == ' ')	/* convert spaces to '+' */
			*d++ = '+';
		else *d++ = *s;		/* copy normal character */
		s++;
	}
	*d = 0;
	return (d);
}

/*
// format_args: format a HTTP argument string from a PrimeNet v5 packet
*/

int format_args (char* args, short operation, void* pkt)
{
	char	*p;

/* Format the common header */

	sprintf (args, "v=%.2f&px=GIMPS", PRIMENET_TRANSACTION_API_VERSION);
	p = args + strlen (args);

/* Format the message dependent args */

	switch (operation) {
	case PRIMENET_UPDATE_COMPUTER_INFO:	/* update computer info */
		{
		struct primenetUpdateComputerInfo *z;

		z = (struct primenetUpdateComputerInfo *) pkt;

//		if (!_stricmp (z->user_id, "ANONYMOUS"))
//			strcpy (z->user_id, "admin_user_anon");

		strcpy (p, "&t=uc&g=");
		p = armor (p + strlen (p), z->computer_guid);
		strcpy (p, "&hg=");
		p = armor (p + strlen (p), z->hardware_guid);
		strcpy (p, "&wg=");
		p = armor (p + strlen (p), z->windows_guid);
		strcpy (p, "&a=");
		p = armor (p + strlen (p), z->application);
		strcpy (p, "&c=");
		p = armor (p + strlen (p), z->cpu_model);
		strcpy (p, "&f=");
		p = armor (p + strlen (p), z->cpu_features);
		sprintf (p, "&L1=%d&L2=%d&np=%d&hp=%d&m=%d&s=%d&h=%d&r=%d",
			 z->L1_cache_size, z->L2_cache_size, z->num_cpus,
			 z->num_hyperthread, z->mem_installed, z->cpu_speed,
			 z->hours_per_day, z->rolling_average);
		p += strlen (p);
		if (z->L3_cache_size > 0) {
			sprintf (p, "&L3=%d", z->L3_cache_size);
			p += strlen (p);
		}
		if (z->user_id[0]) {
			strcpy (p, "&u=");
			p = armor (p + strlen (p), z->user_id);
		}
		if (z->computer_name[0]) {
			strcpy (p, "&cn=");
			p = armor (p + strlen (p), z->computer_name);
		}
		break;
		}
	case PRIMENET_PROGRAM_OPTIONS:
		{
		struct primenetProgramOptions *z;

		z = (struct primenetProgramOptions *) pkt;
		strcpy (p, "&t=po&g=");
		p = armor (p + strlen (p), z->computer_guid);
		if (z->cpu_num != -1) {
			sprintf (p, "&c=%d", z->cpu_num);
			p += strlen (p);
		}
		if (z->num_workers != -1) {
			sprintf (p, "&nw=%d", z->num_workers);
			p += strlen (p);
		}
		if (z->work_preference != -1) {
			sprintf (p, "&w=%d", z->work_preference);
			p += strlen (p);
		}
		if (z->priority != -1) {
			sprintf (p, "&Priority=%d", z->priority);
			p += strlen (p);
		}
		if (z->daysOfWork != -1) {
			sprintf (p, "&DaysOfWork=%d", z->daysOfWork);
			p += strlen (p);
		}
		if (z->dayMemory != -1) {
			sprintf (p, "&DayMemory=%d", z->dayMemory);
			p += strlen (p);
		}
		if (z->nightMemory != -1) {
			sprintf (p, "&NightMemory=%d", z->nightMemory);
			p += strlen (p);
		}
		if (z->dayStartTime != -1) {
			sprintf (p, "&DayStartTime=%d", z->dayStartTime);
			p += strlen (p);
		}
		if (z->nightStartTime != -1) {
			sprintf (p, "&NightStartTime=%d", z->nightStartTime);
			p += strlen (p);
		}
		if (z->runOnBattery != -1) {
			sprintf (p, "&RunOnBattery=%d", z->runOnBattery);
			p += strlen (p);
		}
		break;
		}
	case PRIMENET_REGISTER_ASSIGNMENT:	/* register assignment */
		{
		struct primenetRegisterAssignment *z;

		z = (struct primenetRegisterAssignment *) pkt;
		strcpy (p, "&t=ra&g=");
		p = armor (p + strlen (p), z->computer_guid);
		sprintf (p, "&c=%d&w=%d", z->cpu_num, z->work_type);
		p = p + strlen (p);
		if (z->work_type == PRIMENET_WORK_TYPE_FACTOR) {
			sprintf (p, "&n=%d&sf=%g",
				 z->n, z->how_far_factored);
			p = p + strlen (p);
			if (z->factor_to != 0.0) {
				sprintf (p, "&ef=%g", z->factor_to);
				p = p + strlen (p);
			}
		}
		if (z->work_type == PRIMENET_WORK_TYPE_PFACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&sf=%g&saved=%g",
				 z->k, z->b, z->n, z->c, z->how_far_factored,
				 z->tests_saved);
			p = p + strlen (p);
		}
		if (z->work_type == PRIMENET_WORK_TYPE_FIRST_LL ||
		    z->work_type == PRIMENET_WORK_TYPE_DBLCHK) {
			sprintf (p, "&n=%d&sf=%g&p1=%d",
				 z->n, z->how_far_factored,
				 z->has_been_pminus1ed);
			p = p + strlen (p);
		}
		if (z->work_type == PRIMENET_WORK_TYPE_PMINUS1) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p = p + strlen (p);
			if (z->B2 != 0.0) {
				sprintf (p, "&B2=%.0f", z->B2);
				p = p + strlen (p);
			}
		}
		if (z->work_type == PRIMENET_WORK_TYPE_ECM) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p = p + strlen (p);
			if (z->B2 != 0.0) {
				sprintf (p, "&B2=%.0f", z->B2);
				p = p + strlen (p);
			}
			sprintf (p, "&CR=%d", z->curves);
			p = p + strlen (p);
		}
		if (z->work_type == PRIMENET_WORK_TYPE_PRP) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&C=%d&sf=%g&saved=%g",
				 z->k, z->b, z->n, z->c, z->how_far_factored,
				 z->tests_saved);
			p = p + strlen (p);
		}
		break;
		}
	case PRIMENET_GET_ASSIGNMENT:		/* get assignment */
		{
		struct primenetGetAssignment *z;

		z = (struct primenetGetAssignment *) pkt;
		strcpy (p, "&t=ga&g=");
		p = armor (p + strlen (p), z->computer_guid);
		sprintf (p, "&c=%d", z->cpu_num);
		p += strlen (p);
		break;
		}
	case PRIMENET_ASSIGNMENT_PROGRESS:
		{
		struct primenetAssignmentProgress *z;

		z = (struct primenetAssignmentProgress *) pkt;
		strcpy (p, "&t=ap&g=");
		p = armor (p + strlen (p), z->computer_guid);
		strcpy (p, "&k=");
		p = armor (p + strlen (p), z->assignment_uid);
		if (z->stage[0]) {
			strcpy (p, "&stage=");
			p = armor (p + strlen (p), z->stage);
		}
		/* Server does not like a pcercent complete of 100%. */
		/* Just in case caller passes that value in, convert it */
		sprintf (p, "&c=%lu&p=%.4f&d=%lu&e=%lu",
			 (unsigned long) z->cpu_num,
			 z->pct_complete < 99.9999 ? z->pct_complete : 99.9999,
			 (unsigned long) z->next_update, (unsigned long) z->end_date);
		p += strlen (p);
		if (z->fftlen) {
			sprintf (p, "&fftlen=%d", z->fftlen);
			p += strlen (p);
		}
		if (z->iteration) {
			sprintf (p, "&iteration=%d", z->iteration);
			p += strlen (p);
			strcpy (p, "&res64=");
			p = armor (p + strlen (p), z->residue);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		break;
		}
	case PRIMENET_ASSIGNMENT_RESULT:
		{
		struct primenetAssignmentResult *z;

		z = (struct primenetAssignmentResult *) pkt;
		strcpy (p, "&t=ar&g=");
		p = armor (p + strlen (p), z->computer_guid);
		if (z->assignment_uid[0]) {
			strcpy (p, "&k=");
			p = armor (p + strlen (p), z->assignment_uid);
		} else {
			strcpy (p, "&k=0");
			p += strlen (p);
		}
		if (z->JSONmessage[0]) {
			strcpy (p, "&m=");
			p = armor (p + strlen (p), z->JSONmessage);
		}
		else if (z->message[0]) {
			strcpy (p, "&m=");
			p = armor (p + strlen (p), z->message);
		}
		sprintf (p, "&r=%d&d=%d", z->result_type, z->done);
		p += strlen (p);
		if (z->result_type == PRIMENET_AR_LL_RESULT) {
			sprintf (p, "&n=%d&sc=%d", z->n, z->shift_count);
			p += strlen (p);
			strcpy (p, "&rd=");
			p = armor (p + strlen (p), z->residue);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		if (z->result_type == PRIMENET_AR_LL_PRIME) {
			sprintf (p, "&n=%d&sc=%d", z->n, z->shift_count);
			p += strlen (p);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
		}
		if (z->result_type == PRIMENET_AR_TF_FACTOR) {
			sprintf (p, "&n=%d&sf=%g", z->n, z->start_bits);
			p += strlen (p);
			strcpy (p, "&f=");
			p = armor (p + strlen (p), z->factor);
		}
		if (z->result_type == PRIMENET_AR_P1_FACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
			strcpy (p, "&f=");
			p = armor (p + strlen (p), z->factor);
		}
		if (z->result_type == PRIMENET_AR_ECM_FACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&CR=%d&B1=%.0f&stage=%d",
				 z->k, z->b, z->n, z->c, z->curves, z->B1, z->stage);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
			strcpy (p, "&f=");
			p = armor (p + strlen (p), z->factor);
		}

		if (z->result_type == PRIMENET_AR_TF_NOFACTOR) {
			sprintf (p, "&n=%d&sf=%g&ef=%g",
				 z->n, z->start_bits, z->end_bits);
			p += strlen (p);
		}
		if (z->result_type == PRIMENET_AR_P1_NOFACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->B1);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
		}
		if (z->result_type == PRIMENET_AR_ECM_NOFACTOR) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d&CR=%d&B1=%.0f",
				 z->k, z->b, z->n, z->c, z->curves, z->B1);
			p += strlen (p);
			if (z->B2) {
				sprintf (p, "&B2=%.0f", z->B2);
				p += strlen (p);
			}
		}
		if (z->result_type == PRIMENET_AR_PRP_RESULT) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d", z->k, z->b, z->n, z->c);
			p += strlen (p);
			strcpy (p, "&rd=");
			p = armor (p + strlen (p), z->residue);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
			if (z->num_known_factors) sprintf (p, "&nkf=%d", z->num_known_factors), p += strlen (p);
			if (z->prp_base) sprintf (p, "&base=%d", z->prp_base), p += strlen (p);
			if (z->prp_residue_type) sprintf (p, "&rt=%d", z->prp_residue_type), p += strlen (p);
			if (z->shift_count) sprintf (p, "&sc=%d", z->shift_count), p += strlen (p);
			if (z->gerbicz) strcpy (p, "&gbz=1"), p += strlen (p);
		}
		if (z->result_type == PRIMENET_AR_PRP_PRIME) {
			sprintf (p, "&A=%.0f&b=%d&n=%d&c=%d", z->k, z->b, z->n, z->c);
			p += strlen (p);
			strcpy (p, "&ec=");
			p = armor (p + strlen (p), z->error_count);
			if (z->num_known_factors) sprintf (p, "&nkf=%d", z->num_known_factors), p += strlen (p);
			if (z->prp_base) sprintf (p, "&base=%d", z->prp_base), p += strlen (p);
			if (z->shift_count) sprintf (p, "&sc=%d", z->shift_count), p += strlen (p);
			if (z->gerbicz) strcpy (p, "&gbz=1"), p += strlen (p);
		}
		if (z->fftlen) {
			sprintf (p, "&fftlen=%d", z->fftlen);
			p += strlen (p);
		}
//bug - can we support a 0 result_type that only sends a message.
//might need this for sending old results file
		break;
		}
	case PRIMENET_ASSIGNMENT_UNRESERVE:	/* assignment unreserve */
		{
		struct primenetAssignmentUnreserve *z;

		z = (struct primenetAssignmentUnreserve *) pkt;
		strcpy (p, "&t=au&g=");
		p = armor (p + strlen (p), z->computer_guid);
		strcpy (p, "&k=");
		p = armor (p + strlen (p), z->assignment_uid);
		break;
		}
	case PRIMENET_BENCHMARK_DATA:
		{
		struct primenetBenchmarkData *z;
		unsigned int i;

		z = (struct primenetBenchmarkData *) pkt;
		strcpy (p, "&t=bd&g=");
		p = armor (p + strlen (p), z->computer_guid);
		if (z->user_comment[0]) {
			strcpy (p, "&c=");
			p = armor (p + strlen (p), z->user_comment);
		}
		for (i = 0; i < z->num_data_points; i++) {
			*p++ = '&';
			p = armor (p, z->data_points[i].bench);
			sprintf (p, "=%f", z->data_points[i].timing);
			p += strlen (p);
		}
		break;
		}
	case PRIMENET_PING_SERVER:
		{
		struct primenetPingServer *z;

		z = (struct primenetPingServer *) pkt;
		sprintf (p, "&t=ps&q=%d", z->ping_type);
		p += strlen (p);
		break;
		}
	}

/* Append the security string */

#ifdef _V5_SECURITY_MODULE_PRESENT_
	{
		char	p_v5key[33];
		make_v5_client_key (p_v5key, COMPUTER_GUID);
		secure_v5_url (args, p_v5key);
	}
#endif
	return (0);
}


/*////////////////////////////////////////////////////////////////////////////
//
// HTTP downloaded response page parsing procedures
//
/////////////////////////////////////////////////////////////////////////////*/

/* skip over the token name and point to the data string */

char* skip_token (char *s)
{
	while (*s && *s != '=') s++;
	if (*s == '=') s++;
	return (s);
}


/* copy the data string up to the next '\r' delimiter character */

char* copy_value (char *buf, char *s)
{
	while (*s && *s != '\r') *buf++ = *s++;
	if (*s == '\r') s++;
	*buf = 0;
	return (s);
}

/* parse various tokens from the response */

char *find_id (
	char	*buf,
	char	*id)
{
	unsigned int idlen;
	idlen = (unsigned int) strlen (id);
	while (*buf) {
		if (memcmp (buf, id, idlen) == 0 && buf[idlen] == '=')
			return (buf + idlen + 1);
		while (*buf && *buf != '>' && *buf != '\n') buf++;
		while (*buf && (*buf == '>' || *buf == '\r' || *buf == '\n' || *buf == ' ')) buf++;
	}
	return (NULL);
}

int parse_string (
	char	*buf,
	char	*id,
	char	*result_buf,
	unsigned int result_buflen)
{
	unsigned int i;
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	for  (i = 0; i < result_buflen-1 && buf[i] && buf[i] != '\n'; i++)
		result_buf[i] = (buf[i] == '\r' ? '\n' : buf[i]);
	result_buf[i] = 0;
//bug -raise error if string too long
	return (TRUE);
}

int parse_multiline_string (
	char	*buf,
	char	*id,
	char	*result_buf,
	unsigned int result_buflen)
{
	unsigned int i;
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	for  (i = 0; i < result_buflen-1 && buf[i]; i++) {
		if (buf[i] == '\n') {
			char	*equal_sign, *next_newline;
			equal_sign = strchr (buf+i+1, '=');
			next_newline = strchr (buf+i+1, '\n');
			if (equal_sign < next_newline ||
			    (equal_sign != NULL && next_newline == NULL) ||
			    (equal_sign == NULL && next_newline == NULL && buf[i+1] == 0))
				break;
		}
		result_buf[i] = (buf[i] == '\r' ? '\n' : buf[i]);
	}
	result_buf[i] = 0;
//bug -raise error if string too long
	return (TRUE);
}

int parse_int (
	char	*buf,
	char	*id,
	int32_t	*result)
{
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	if (buf[0] != '-' && (buf[0] < '0' || buf[0] > '9')) return (FALSE);
	*result = atoi (buf);
//bug -raise error if not integer
	return (TRUE);
}

int parse_uint (
	char	*buf,
	char	*id,
	uint32_t *result)
{
	buf = find_id (buf, id);
	if (buf == NULL) return (FALSE);
	if (buf[0] < '0' || buf[0] > '9') return (FALSE);
	*result = atoi (buf);
//bug -raise error if not integer
	return (TRUE);
}

void parse_double (
	char	*buf,
	char	*id,
	double	*result)
{
	buf = find_id (buf, id);
	if (buf == NULL) return;
	*result = atof (buf);
//bug -raise error if not number
}
	


/*
// parse_page: reads the server response page tokens and values
//             and converts these back into a C structure
*/

int parse_page (char *response_buf, short operation, void *pkt)

{
	char	*s;
	char	buf[400], errtxt[200];
	int32_t	res;

/* Get result code, which is always first */

	s = response_buf;
	if (!parse_int (s, "pnErrorResult", &res)) {
		LogMsg ("PnErrorResult value missing.  Full response was:\n");
		LogMsg (response_buf);

		/* Look for PHP timeout response and convert it to */
		/* a server busy error code. */

		if (strstr (s, "execution time") != NULL &&
		    strstr (s, "exceeded") != NULL)
			return (PRIMENET_ERROR_SERVER_BUSY);
			
		return (PRIMENET_ERROR_PNERRORRESULT);
	}
	if (!parse_multiline_string (s, "pnErrorDetail", errtxt, sizeof (errtxt))) {
		LogMsg ("PnErrorDetail string missing\n");
		return (PRIMENET_ERROR_PNERRORDETAIL);
	}

/* If result is non-zero print out an error message */

	if (res) {
		char	buf[2000];
		char	*resmsg;

/* Convert the error number to text */

		switch (res) {
		case PRIMENET_ERROR_SERVER_BUSY:
			resmsg = "Server busy";
			break;
		case PRIMENET_ERROR_INVALID_VERSION:
			resmsg = "Invalid version";
			break;
		case PRIMENET_ERROR_INVALID_TRANSACTION:
			resmsg = "Invalid transaction";
			break;
		case PRIMENET_ERROR_INVALID_PARAMETER:
			resmsg = "Invalid parameter";
			break;
		case PRIMENET_ERROR_ACCESS_DENIED:
			resmsg = "Access denied";
			break;
		case PRIMENET_ERROR_DATABASE_CORRUPT:
			resmsg = "Server database malfunction";
			break;
		case PRIMENET_ERROR_DATABASE_FULL_OR_BROKEN:
			resmsg = "Server database full or broken";
			break;
		case PRIMENET_ERROR_INVALID_USER:
			resmsg = "Invalid user";
			break;
		case PRIMENET_ERROR_UNREGISTERED_CPU:
			resmsg = "CPU not registered";
			break;
		case PRIMENET_ERROR_OBSOLETE_CLIENT:
			resmsg = "Obsolete client, please upgrade";
			break;
		case PRIMENET_ERROR_STALE_CPU_INFO:
			resmsg = "Stale cpu info";
			break;
		case PRIMENET_ERROR_CPU_IDENTITY_MISMATCH:
			resmsg = "CPU identity mismatch";
			break;
		case PRIMENET_ERROR_CPU_CONFIGURATION_MISMATCH:
			resmsg = "CPU configuration mismatch";
			break;
		case PRIMENET_ERROR_NO_ASSIGNMENT:
			resmsg = "No assignment";
			break;
		case PRIMENET_ERROR_INVALID_ASSIGNMENT_KEY:
			resmsg = "Invalid assignment key";
			break;
		case PRIMENET_ERROR_INVALID_ASSIGNMENT_TYPE:
			resmsg = "Invalid assignment type";
			break;
		case PRIMENET_ERROR_INVALID_RESULT_TYPE:
			resmsg = "Invalid result type";
			break;
		default:
			resmsg = "Unknown error code";
			break;
		}

/* Print out the error code, text, and details */
		
		sprintf (buf, "PrimeNet error %d: %s\n", res, resmsg);
		LogMsg (buf);
		sprintf (buf, "%s\n", errtxt);
		LogMsg (buf);
	}

/* If there was no error code but there was some error text, then print */
/* the error text. */

	else if (strcmp (errtxt, "SUCCESS")) {
		LogMsg ("PrimeNet success code with additional info:\n");
		sprintf (buf, "%s\n", errtxt);
		LogMsg (buf);
	}

/* Parse remaining response parameters */

	switch (operation) {
	case PRIMENET_UPDATE_COMPUTER_INFO:	/* update computer info */
		{
		struct primenetUpdateComputerInfo *z;

		z = (struct primenetUpdateComputerInfo *) pkt;
		parse_string (s, "g", z->computer_guid, sizeof (z->computer_guid));
		parse_string (s, "u", z->user_id, sizeof (z->user_id));
		parse_string (s, "un", z->user_name, sizeof (z->user_name));
		parse_string (s, "cn", z->computer_name, sizeof (z->computer_name));
		parse_uint (s, "od", &z->options_counter);

		if (!strcmp (z->user_id, "admin_user_anon"))
			strcpy (z->user_id, "ANONYMOUS");

		break;
		}
	case PRIMENET_PROGRAM_OPTIONS:
		{
		struct primenetProgramOptions *z;

		z = (struct primenetProgramOptions *) pkt;
		z->num_workers = -1;
		parse_int (s, "nw", &z->num_workers);
		z->work_preference = -1;
		parse_int (s, "w", &z->work_preference);
		z->priority = -1;
		parse_int (s, "Priority", &z->priority);
		z->daysOfWork = -1;
		parse_int (s, "DaysOfWork", &z->daysOfWork);
		z->dayMemory = -1;
		parse_int (s, "DayMemory", &z->dayMemory);
		z->nightMemory = -1;
		parse_int (s, "NightMemory", &z->nightMemory);
		z->dayStartTime = -1;
		parse_int (s, "DayStartTime", &z->dayStartTime);
		z->nightStartTime = -1;
		parse_int (s, "NightStartTime", &z->nightStartTime);
		z->runOnBattery = -1;
		parse_int (s, "RunOnBattery", &z->runOnBattery);
		parse_uint (s, "od", &z->options_counter);
		break;
		}
	case PRIMENET_GET_ASSIGNMENT:
		{
		struct primenetGetAssignment *z;

		z = (struct primenetGetAssignment *) pkt;
		parse_string (s, "k", z->assignment_uid, sizeof (z->assignment_uid));
		parse_uint (s, "w", &z->work_type);
		parse_double (s, "A", &z->k);
		parse_uint (s, "b", &z->b);
		parse_uint (s, "n", &z->n);
		parse_int (s, "c", &z->c);
		parse_uint (s, "p1", &z->has_been_pminus1ed);
		parse_double (s, "sf", &z->how_far_factored);
		parse_double (s, "ef", &z->factor_to);
		parse_double (s, "B1", &z->B1);
		parse_double (s, "B2", &z->B2);
		parse_uint (s, "CR", &z->curves);
		parse_double (s, "saved", &z->tests_saved);
		parse_uint (s, "base", &z->prp_base);
		parse_uint (s, "rt", &z->prp_residue_type);
		parse_uint (s, "dc", &z->prp_dblchk);
		parse_string (s, "kf", z->known_factors, sizeof (z->known_factors));
		break;
		}
	case PRIMENET_REGISTER_ASSIGNMENT:
		{
		struct primenetRegisterAssignment *z;

		z = (struct primenetRegisterAssignment *) pkt;
		parse_string (s, "k", z->assignment_uid, sizeof (z->assignment_uid));
		break;
		}
	case PRIMENET_ASSIGNMENT_PROGRESS:
		break;
	case PRIMENET_ASSIGNMENT_RESULT:
		break;
	case PRIMENET_ASSIGNMENT_UNRESERVE:
		break;
	case PRIMENET_BENCHMARK_DATA:
		break;
	case PRIMENET_PING_SERVER:
		{
		struct primenetPingServer *z;

		z = (struct primenetPingServer *) pkt;
		parse_string (s, "r", z->ping_response, sizeof (z->ping_response));
		break;
		}
	}
	return (res);
}

/*
// LoadPrimenet: call from main thread at startup
*/

void LoadPrimenet ()
{

/* The cURL documentation strongly recommends calling curl_global_init */
/* from the main thread rather than letting curl_easy_init do the */
/* initialization when other threads are running. */

	if (IniSectionGetInt (INI_FILE, iniSection, "UseCURL", 1))
		curl_global_init (CURL_GLOBAL_ALL);
}


/*
// Primenet: main interface to Prime95.exe
*/

int PRIMENET (short operation, void *pkt)
{
	int	status;
	char args[4096];		/* formatted arguments buffer */
	char pbuf[4096];		/* return page buffer */

/* Assemble GET/POST arguments */

	status = format_args (args, operation, pkt);
	if (status) return (status);

/* Send arguments, read back resulting page */

	if (IniSectionGetInt (INI_FILE, iniSection, "UseCURL", 1))
		status = pnHttpServerCURL (pbuf, sizeof (pbuf), args);
	else
		status = pnHttpServer (pbuf, sizeof (pbuf), args);
	if (status) return (status);

/* Extract results from returned page into packet */

	return (parse_page (pbuf, operation, pkt));
}

