; Copyright 1995-2009 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use a simple linear memory model and
; the simplified normalization code.  FFT sizes from 32 doubles to
; 128 doubles are supported.
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

_TEXT SEGMENT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE	lucas.mac
INCLUDE pfa.mac
INCLUDE mult.mac
INCLUDE fft1.mac
INCLUDE memory.mac

	flat_distances

;; All the FFT routines for each FFT length

	fft	32
	fft	32p
	fft	40
	fft	48
	fft	48p
	fft	56
	fft	64
	fft	64p
	fft	80
	fft	96
	fft	96p
	fft	112
	fft	128
	fft	128p
	fft	160
	fft	192
	fft	192p
	fft	224
	fft	256
	fft	256p
	fft	320
	fft	384
	fft	384p
	fft	448
	fft	512
	fft	512p
	fft	640
	fft	768
	fft	768p
	fft	896
	fft	1024
	fft	1024p
	fft	1280
	fft	1536
	fft	1536p
	fft	1792
	fft	2048
	fft	2048p
	fft	2560
	fft	3072
	fft	3072p
	fft	3584
	fft	4096
	fft	4096p

_TEXT	ENDS
END
