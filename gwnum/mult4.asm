; Copyright 1995-2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that uses a scratch area in pass 1.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

	.686
	.XMM
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE	lucas.mac
INCLUDE pfa.mac
INCLUDE mult.mac
INCLUDE pass1.mac
INCLUDE fft4.mac
INCLUDE memory.mac

EXTRN gw_finish_fft:PROC
EXTRN gw_carries:PROC
EXTRN gw_finish_mult:PROC

EXTRNP	pass2_12_levels

_TEXT SEGMENT

	flat_distances

;; All the FFT routines for each FFT length

;; Distance between two pass 2 data blocks.  Pass 2 does 12 FFT levels.
;; 2^12 complex values = 2^13 doubles = 64KB.

blkdst = (16*(4096+64)+64)

;	fft	80K
;	fft	96K
;	fft	96Kp
;	fft	112K
;	fft	128K
;	fft	128Kp
;	fft	160K
;	fft	192K
;	fft	192Kp
;	fft	224K
;	fft	256K
;	fft	256Kp
;	fft	320K
;	fft	384K
;	fft	384Kp
;	fft	448K
;	fft	512K
;	fft	512Kp

INCLUDE pass1scr.mac

;	fft	640K
;	fft	768K
;	fft	768Kp
;	fft	896K
;	fftclm	1024K, 2
;	fft	1024Kp
	fft	1280K
	fft	1536K
	fft	1536Kp
	fft	1792K
	fftclm	2048K, 2
	fft	2048Kp
	fftclm	2560K, 2
	fftclm	3072K, 2
	fftclm	3072Kp, 2
	fftclm	3584K, 2
	fftclm	4096K, 2
	fftclm	4096Kp, 2

_TEXT	ENDS
END
