; Copyright 2011-2019 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the AVX-512 radix-4/8 DJB FFT with delayed sin/cos multiplies and partial normalization.
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
INCLUDE zr4dwpnpass2.mac

_TEXT SEGMENT

;; Generate real shared pass 1 routines optimized for this architecture

buildfor SKX, zpass1gen 2,0,,1		; Build real one-pass FFT wrapper

build421 SKX, SKX, SKX,	zpass1gen 128,0,
build421 SKX, SKX, SKX,	zpass1gen 192,0,
build421 SKX, SKX, SKX,	zpass1gen 640,0,
build421 SKX, SKX, SKX,	zpass1gen 768,0,
build421 SKX, SKX, SKX,	zpass1gen 896,0,
build421 SKX, SKX, SKX,	zpass1gen 960,0,
build421 SKX, SKX, SKX,	zpass1gen 1024,0,
build421 SKX, SKX, SKX,	zpass1gen 1152,0,
build421 SKX, SKX, SKX,	zpass1gen 1280,0,
build421 SKX, SKX, SKX,	zpass1gen 1344,0,
build421 SKX, SKX, SKX,	zpass1gen 1536,0,
build421 SKX, SKX, SKX,	zpass1gen 1920,0,
build421 SKX, SKX, SKX,	zpass1gen 2048,0,
build421 SKX, SKX, SKX,	zpass1gen 2304,0,
build421 SKX, SKX, SKX,	zpass1gen 3072,0,

;; Generate all-complex shared pass 1 routines optimized for this architecture

buildfor SKX, zpass1gen 2,1,,1		; Build complex one-pass FFT wrapper

build421 SKX, SKX, SKX,	zpass1gen 128,1,
build421 SKX, SKX, SKX,	zpass1gen 192,1,
build421 SKX, SKX, SKX,	zpass1gen 640,1,
build421 SKX, SKX, SKX,	zpass1gen 768,1,
build421 SKX, SKX, SKX,	zpass1gen 896,1,
build421 SKX, SKX, SKX,	zpass1gen 960,1,
build421 SKX, SKX, SKX,	zpass1gen 1024,1,
build421 SKX, SKX, SKX,	zpass1gen 1152,1,
build421 SKX, SKX, SKX,	zpass1gen 1280,1,
build421 SKX, SKX, SKX,	zpass1gen 1344,1,
build421 SKX, SKX, SKX,	zpass1gen 1536,1,
build421 SKX, SKX, SKX,	zpass1gen 1920,1,
build421 SKX, SKX, SKX,	zpass1gen 2048,1,
build421 SKX, SKX, SKX,	zpass1gen 2304,1,
build421 SKX, SKX, SKX,	zpass1gen 3072,1,

_TEXT	ENDS
END
