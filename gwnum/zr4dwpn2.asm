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
INCLUDE zr4.mac
INCLUDE zr4dwpnpass1sc.mac
INCLUDE zr4dwpnpass2.mac

_TEXT SEGMENT

;; Generate real and complex pass 2 routines optimized for this architecture

buildfor SKX,	zpass2gen 64
buildfor SKX,	zpass2gen 320
buildfor SKX,	zpass2gen 384
buildfor SKX,	zpass2gen 448
buildfor SKX,	zpass2gen 512
buildfor SKX,	zpass2gen 640
buildfor SKX,	zpass2gen 768
buildfor SKX,	zpass2gen 1024
buildfor SKX,	zpass2gen 1600
buildfor SKX,	zpass2gen 1920
buildfor SKX,	zpass2gen 2240
buildfor SKX,	zpass2gen 2304
buildfor SKX,	zpass2gen 2560
buildfor SKX,	zpass2gen 2688
buildfor SKX,	zpass2gen 3072
buildfor SKX,	zpass2gen 3136
buildfor SKX,	zpass2gen 3200
buildfor SKX,	zpass2gen 3584
buildfor SKX,	zpass2gen 3840
buildfor SKX,	zpass2gen 4096
buildfor SKX,	zpass2gen 4480
buildfor SKX,	zpass2gen 4608
buildfor SKX,	zpass2gen 5120
buildfor SKX,	zpass2gen 5376
buildfor SKX,	zpass2gen 6144
buildfor SKX,	zpass2gen 6400
buildfor SKX,	zpass2gen 7168
buildfor SKX,	zpass2gen 7680
buildfor SKX,	zpass2gen 8192
buildfor SKX,	zpass2gen 9216
buildfor SKX,	zpass2gen 10240
buildfor SKX,	zpass2gen 12288
buildfor SKX,	zpass2gen 12800
buildfor SKX,	zpass2gen 15360
buildfor SKX,	zpass2gen 16384
buildfor SKX,	zpass2gen 17920
buildfor SKX,	zpass2gen 18432
buildfor SKX,	zpass2gen 20480
buildfor SKX,	zpass2gen 21504
buildfor SKX,	zpass2gen 24576
buildfor SKX,	zpass2gen 25088
buildfor SKX,	zpass2gen 28672
buildfor SKX,	zpass2gen 30720
buildfor SKX,	zpass2gen 32768

_TEXT	ENDS
END
