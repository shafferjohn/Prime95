; Copyright 2011-2017 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the small AVX FFTs.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

yfft_type TEXTEQU <r4>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE ymult.mac
INCLUDE yonepass.mac
INCLUDE yr4.mac

_TEXT SEGMENT

;; Implement the small one pass FFTs

PREFETCHING = 0

buildfor CORE + FMA3_64,	yonepass 32, 0
buildfor CORE + FMA3_64,	yonepass 64, 0
buildfor CORE + FMA3_64,	yonepass 96, 0
buildfor CORE + FMA3_64,	yonepass 128, 0
buildfor CORE + FMA3_64,	yonepass 160, 0
buildfor CORE + FMA3_64,	yonepass 192, 0
buildfor CORE + FMA3_64,	yonepass 256, 0
buildfor CORE + FMA3_64,	yonepass 320, 0
buildfor CORE + FMA3_64,	yonepass 384, 0
buildfor CORE + FMA3_64,	yonepass 512, 0
buildfor CORE + FMA3_64,	yonepass 640, 0
buildfor CORE + FMA3_64,	yonepass 768, 0
buildfor CORE + FMA3_64,	yonepass 1K, 0
buildfor CORE + FMA3_64,	yonepass 1280, 0
buildfor CORE + FMA3_64,	yonepass 1536, 0
buildfor CORE + FMA3_64,	yonepass 2K, 0
buildfor CORE + FMA3_64,	yonepass 2560, 0
buildfor CORE + FMA3_64,	yonepass 3K, 0
buildfor CORE + FMA3_64,	yonepass 4K, 0
buildfor CORE + FMA3_64,	yonepass 5K, 0
buildfor               ,	yonepass 6K, 0
buildfor               ,	yonepass 8K, 0
buildfor               ,	yonepass 10K, 0
buildfor               ,	yonepass 12K, 0
buildfor        FMA3_64,	yonepass 16K, 0
buildfor               ,	yonepass 18K, 0
buildfor        FMA3_64,	yonepass 20K, 0
buildfor               ,	yonepass 24K, 0
buildfor               ,	yonepass 32K, 0

buildfor CORE + FMA3_64,	yonepass 32, 1
buildfor CORE + FMA3_64,	yonepass 64, 1
buildfor CORE + FMA3_64,	yonepass 96, 1
buildfor CORE + FMA3_64,	yonepass 128, 1
buildfor CORE + FMA3_64,	yonepass 160, 1
buildfor CORE + FMA3_64,	yonepass 192, 1
buildfor CORE + FMA3_64,	yonepass 256, 1
buildfor CORE + FMA3_64,	yonepass 320, 1
buildfor CORE + FMA3_64,	yonepass 384, 1
buildfor CORE + FMA3_64,	yonepass 512, 1
buildfor CORE + FMA3_64,	yonepass 640, 1
buildfor CORE + FMA3_64,	yonepass 768, 1
buildfor CORE + FMA3_64,	yonepass 1K, 1
buildfor CORE + FMA3_64,	yonepass 1280, 1
buildfor CORE + FMA3_64,	yonepass 1536, 1
buildfor CORE + FMA3_64,	yonepass 2K, 1
buildfor CORE + FMA3_64,	yonepass 2560, 1
buildfor CORE + FMA3_64,	yonepass 3K, 1
buildfor CORE + FMA3_64,	yonepass 4K, 1
buildfor CORE + FMA3_64,	yonepass 5K, 1
buildfor CORE_32 + FMA3_64,	yonepass 6K, 1
buildfor CORE_32 + FMA3_64,	yonepass 8K, 1
buildfor CORE_32 + FMA3_64,	yonepass 10K, 1
buildfor CORE_32 + FMA3_64,	yonepass 12K, 1
buildfor CORE_32 + FMA3_64,	yonepass 16K, 1
buildfor               ,	yonepass 20K, 1
buildfor               ,	yonepass 18K, 1
buildfor               ,	yonepass 24K, 1
buildfor               ,	yonepass 32K, 1

_TEXT	ENDS
END
