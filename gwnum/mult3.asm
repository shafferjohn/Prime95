; Copyright 1995-2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This code handles FFTs that use two passes with 10 levels done on the
; second pass.
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
INCLUDE fft3.mac
INCLUDE memory.mac

EXTRN gw_finish_fft:PROC
EXTRN gw_carries:PROC
EXTRN gw_finish_mult:PROC

EXTRNP	pass2_10_levels

_TEXT SEGMENT

	flat_distances

;; Distance between two pass 2 data blocks.  Pass 2 does 10 FFT levels.
;; 2^10 complex values = 2^11 doubles = 16KB.

blkdst = (4*(4096+64)+64)

;; All the FFT routines for each FFT length.

IFNDEF PFETCH
;	fft	20K
;	fft	24K
;	fft	24Kp
;	fft	28K
;	fft	32K
;	fft	32Kp
ENDIF
;	fft	40K
;	fft	48K
;	fft	48Kp
;	fft	56K
;	fft	64K
;	fft	64Kp
;	fft	80K
	fft	96K
	fft	96Kp
	fft	112K
	fft	128K
	fft	128Kp

INCLUDE pass1scr.mac

	fft	160K
	fft	192K
	fft	192Kp
	fft	224K
	fft	256K
	fft	256Kp
	fft	320K
	fft	384K
	fft	384Kp
	fft	448K
	fft	512K
	fft	512Kp
	fft	640K
	fft	768K
	fft	768Kp
	fft	896K
	fftclm	1024K, 2
	fft	1024Kp

_TEXT	ENDS
END
