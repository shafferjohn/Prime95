; Copyright 2011-2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the small AVX-512 FFTs.
;

	TITLE   setup

zfft_type TEXTEQU <r4dwpn>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE zonepass.mac
INCLUDE zr4.mac
INCLUDE zr4dwpnpass1sc.mac

_TEXT SEGMENT

;; Implement the small one pass FFTs

;PREFETCHING = 0   ;; this turns on the _ip suffix

;; BUG??? move to zr4dwpn.asm????

buildfor SKX,   zpass1gen 2,0,,1		; Build real wrapper

IFDEF NOT_YET_IMPL
buildfor SKX,	zonepass 128, 0
buildfor SKX,	zonepass 640, 0
buildfor SKX,	zonepass 768, 0
buildfor SKX,	zonepass 896, 0
buildfor SKX,	zonepass 160, 0
buildfor SKX,	zonepass 192, 0
buildfor SKX,	zonepass 256, 0
buildfor SKX,	zonepass 320, 0
buildfor SKX,	zonepass 384, 0
buildfor SKX,	zonepass 512, 0
buildfor SKX,	zonepass 1K, 0
buildfor SKX,	zonepass 1280, 0
buildfor SKX,	zonepass 1536, 0
buildfor SKX,	zonepass 2K, 0
buildfor SKX,	zonepass 2560, 0
buildfor SKX,	zonepass 3K, 0
buildfor SKX,	zonepass 4K, 0
buildfor SKX,	zonepass 5K, 0
buildfor     ,	zonepass 6K, 0
buildfor SKX,	zonepass 8K, 0
buildfor     ,	zonepass 10K, 0
buildfor     ,	zonepass 12K, 0
buildfor     ,	zonepass 16K, 0
buildfor     ,	zonepass 20K, 0
buildfor     ,	zonepass 18K, 0
buildfor     ,	zonepass 24K, 0
buildfor     ,	zonepass 32K, 0
ENDIF

;buildfor SKX,   zpass1gen 2,1,,1		; Build complex wrapper

IFDEF NOT_YET_IMPL
buildfor SKX,	zonepass 128, 1
buildfor SKX,	zonepass 640, 1
buildfor SKX,	zonepass 768, 1
buildfor SKX,	zonepass 896, 1
buildfor SKX,	zonepass 1K, 1
buildfor SKX,	zonepass 8K, 1
buildfor SKX,	zonepass 160, 1
buildfor SKX,	zonepass 192, 1
buildfor SKX,	zonepass 256, 1
buildfor SKX,	zonepass 320, 1
buildfor SKX,	zonepass 384, 1
buildfor SKX,	zonepass 512, 1
buildfor SKX,	zonepass 1280, 1
buildfor SKX,	zonepass 1536, 1
buildfor SKX,	zonepass 2K, 1
buildfor SKX,	zonepass 2560, 1
buildfor SKX,	zonepass 3K, 1
buildfor SKX,	zonepass 4K, 1
buildfor SKX,	zonepass 5K, 1
buildfor SKX,	zonepass 6K, 1
buildfor SKX,	zonepass 10K, 1
buildfor SKX,	zonepass 12K, 1
buildfor SKX,	zonepass 16K, 1
buildfor     ,	zonepass 20K, 1
buildfor     ,	zonepass 18K, 1
buildfor     ,	zonepass 24K, 1
buildfor     ,	zonepass 32K, 1
ENDIF

_TEXT	ENDS
END
