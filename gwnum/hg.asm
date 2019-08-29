; Copyright 2001-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4 home-grown FFTs.  These were developed over many years
; and used up until version 25 of the gwnum library.  The one-pass FFTs and a few
; of the smaller two-pass FFTs are still in use today.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <hg>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE hg.mac
INCLUDE hgonepass.mac
INCLUDE hgpass1.mac
INCLUDE hgpass1sc.mac
INCLUDE hgpass2.mac

EXTRN	xgw_carries:PROC

_TEXT SEGMENT

;; Implement the small one pass FFTs, for now only assemble one-pass FFTs for the BLEND architecture

IF @INSTR(,%xarch,<BLEND>) NE 0

PREFETCHING = 0

	xonepass 32, 0
	xonepass 48, 0
	xonepass 64, 0
	xonepass 80, 0
	xonepass 96, 0
	xonepass 112, 0
	xonepass 128, 0
	xonepass 160, 0
	xonepass 192, 0
	xonepass 224, 0
	xonepass 256, 0
	xonepass 320, 0
	xonepass 384, 0
	xonepass 448, 0
	xonepass 512, 0
	xonepass 640, 0
	xonepass 768, 0
	xonepass 896, 0
	xonepass 1024, 0
	xonepass 1280, 0
	xonepass 1536, 0
	xonepass 1792, 0
	xonepass 2048, 0
	xonepass 2560, 0
	xonepass 3072, 0
	xonepass 3584, 0
	xonepass 4096, 0
	xonepass 5120, 0
	xonepass 6144, 0
	xonepass 7168, 0
	xonepass 8192, 0

	xonepass 32, 1
	xonepass 48, 1
	xonepass 64, 1
	xonepass 96, 1
	xonepass 128, 1
	xonepass 192, 1
	xonepass 256, 1
	xonepass 384, 1
	xonepass 512, 1
	xonepass 768, 1
	xonepass 1024, 1
	xonepass 1536, 1
	xonepass 2048, 1
	xonepass 3072, 1
	xonepass 4096, 1
	xonepass 6144, 1
	xonepass 8192, 1

ELSE

;; Generate pass 2 routines optmized for this architecture

buildfor CORE    + P4    + P4TP + K8    + K10,	xpass2gen 8
buildfor CORE    + P4    + P4TP + K8         ,	xpass2gen 10
buildfor CORE    + P4    + P4TP + K8_32      ,	xpass2gen 11
buildfor CORE_64 + P4    + P4TP + K8_32      ,	xpass2gen 12
buildfor                 + P4TP              ,	xpass2gen 13

;; Routines for many FFT sizes

buildfor CORE    + P4    + P4TP + K8    + K10,	hg_pass1levels6pfa5 10K, 8, 4
buildfor CORE    + P4    + P4TP + K8    + K10,	hg_pass1levels6pfa6 12K, 8, 4
buildfor CORE_64 + P4    + P4TP + K8_32      ,	hg_pass1levels6pfa7 14K, 8, 4
buildfor         + P4_32                     ,	hg_pass1levels6pfa8 16K, 8, 4
buildfor         + P4_32                     ,	hg_pass1levels7pfa5 20K, 8, 4
buildfor         + P4    + P4TP              ,	hg_pass1levels7pfa6 24K, 8, 4
buildfor CORE    + P4_32 + P4TP              ,	hg_pass1levels7pfa7 28K, 8, 4
buildfor         + P4                        ,	hg_pass1levels7pfa8 32K, 8, 4

buildfor CORE_64 + P4    + P4TP + K8    + K10,	hg_pass1levels6complex3	12K, 8, 4
buildfor                                     ,	hg_pass1levels6complex4 16K, 8, 4
buildfor         + P4                        ,	hg_pass1levels7complex3 24K, 8, 4
buildfor         + P4                        ,	hg_pass1levels7complex4 32K, 8, 4

buildfor                 + P4TP              ,	hg_pass1sclevels8pfa5 160K, 10, 4
buildfor                                     ,	hg_pass1sclevels8pfa5 160K, 10, 2
buildfor                 + P4TP              ,	hg_pass1sclevels8pfa6 192K, 10, 4
buildfor CORE            + P4TP              ,	hg_pass1sclevels8pfa6 192K, 10, 2
buildfor         + P4_32                     ,	hg_pass1sclevels8pfa6 192K, 10, 1
buildfor                 + P4TP              ,	hg_pass1sclevels8pfa7 224K, 10, 4
buildfor                 + P4TP + K8         ,	hg_pass1sclevels8pfa7 224K, 10, 2
buildfor                 + P4TP              ,	hg_pass1sclevels8pfa8 256K, 10, 4
buildfor                                     ,	hg_pass1sclevels8pfa8 256K, 10, 2
buildfor         + P4_32 + P4TP              ,	hg_pass1sclevels9pfa5 320K, 10, 4
buildfor                                     ,	hg_pass1sclevels9pfa5 320K, 10, 2
buildfor                        + K8_32      ,	hg_pass1sclevels9pfa5 320K, 10, 1
buildfor         + P4_32 + P4TP              ,	hg_pass1sclevels9pfa6 384K, 10, 4
buildfor                                     ,	hg_pass1sclevels9pfa6 384K, 10, 2
buildfor                                     ,	hg_pass1sclevels9pfa6 384K, 10, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa7 448K, 10, 4
buildfor                                     ,	hg_pass1sclevels9pfa7 448K, 10, 2
buildfor                                     ,	hg_pass1sclevels9pfa7 448K, 10, 1
buildfor         + P4_32 + P4TP              ,	hg_pass1sclevels9pfa8 512K, 10, 4
buildfor                                     ,	hg_pass1sclevels9pfa8 512K, 10, 2
buildfor                                     ,	hg_pass1sclevels9pfa8 512K, 10, 1

buildfor                                     ,	hg_pass1sclevels8complex3 192K, 10, 4
buildfor                                     ,	hg_pass1sclevels8complex3 192K, 10, 2
buildfor                 + P4TP              ,	hg_pass1sclevels8complex3 192K, 10, 1
buildfor                                     ,	hg_pass1sclevels8complex4 256K, 10, 4
buildfor                                     ,	hg_pass1sclevels8complex4 256K, 10, 2
buildfor                                     ,	hg_pass1sclevels8complex4 256K, 10, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9complex3 384K, 10, 4
buildfor                                     ,	hg_pass1sclevels9complex3 384K, 10, 2
buildfor                                     ,	hg_pass1sclevels9complex3 384K, 10, 1
buildfor                                     ,	hg_pass1sclevels9complex4 512K, 10, 4
buildfor                                     ,	hg_pass1sclevels9complex4 512K, 10, 2
buildfor                                     ,	hg_pass1sclevels9complex4 512K, 10, 1

buildfor         + P4_32 + P4TP              ,	hg_pass1levels5pfa5 40K, 11, 4
buildfor CORE            + P4TP              ,	hg_pass1levels5pfa6 48K, 11, 4
buildfor                 + P4TP              ,	hg_pass1levels5pfa7 56K, 11, 4
buildfor                 + P4TP              ,	hg_pass1levels5pfa8 64K, 11, 4
buildfor                                     ,	hg_pass1levels6pfa5 80K, 11, 4
buildfor                                     ,	hg_pass1levels6pfa6 96K, 11, 4
buildfor                                     ,	hg_pass1levels6pfa7 112K, 11, 4
buildfor                                     ,	hg_pass1levels6pfa8 128K, 11, 4
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa5 640K, 11, 4
buildfor                        + K8_32      ,	hg_pass1sclevels9pfa5 640K, 11, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa6 768K, 11, 4
buildfor                                     ,	hg_pass1sclevels9pfa6 768K, 11, 2
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa7 896K, 11, 4
buildfor                                     ,	hg_pass1sclevels9pfa7 896K, 11, 2
buildfor         + P4_32 + P4TP              ,	hg_pass1sclevels9pfa8 1024K, 11, 4
buildfor                                     ,	hg_pass1sclevels9pfa8 1024K, 11, 2
buildfor                                     ,	hg_pass1sclevels10pfa5 1280K, 11, 2
buildfor                                     ,	hg_pass1sclevels10pfa6 1536K, 11, 1
buildfor                                     ,	hg_pass1sclevels10pfa7 1792K, 11, 1

buildfor                 + P4TP              ,	hg_pass1levels5complex3 48K, 11, 4
buildfor                                     ,	hg_pass1levels5complex4 64K, 11, 4
buildfor                 + P4TP              ,	hg_pass1levels6complex3 96K, 11, 4
buildfor                                     ,	hg_pass1levels6complex4	128K, 11, 4
buildfor                                     ,	hg_pass1sclevels9complex3 768K, 11, 4
buildfor                                     ,	hg_pass1sclevels9complex3 768K, 11, 2
buildfor                                     ,	hg_pass1sclevels9complex4 1024K, 11, 4
buildfor                                     ,	hg_pass1sclevels9complex4 1024K, 11, 2
buildfor                                     ,	hg_pass1sclevels10complex3 1536K, 11, 1

buildfor                                     ,	hg_pass1sclevels8pfa5 640K, 12, 2
buildfor         + P4                        ,	hg_pass1sclevels8pfa6 768K, 12, 2
buildfor                        + K8_32      ,	hg_pass1sclevels8pfa7 896K, 12, 2
buildfor                                     ,	hg_pass1sclevels8pfa7 896K, 12, 1
buildfor                                     ,	hg_pass1sclevels8pfa8 1024K, 12, 4
buildfor                                     ,	hg_pass1sclevels8pfa8 1024K, 12, 2
buildfor                                     ,	hg_pass1sclevels8pfa8 1024K, 12, 1
buildfor         + P4_32 + P4TP              ,	hg_pass1sclevels9pfa5 1280K, 12, 4
buildfor                        + K8_32      ,	hg_pass1sclevels9pfa5 1280K, 12, 1
buildfor CORE_64 + P4_32 + P4TP              ,	hg_pass1sclevels9pfa6 1536K, 12, 4
buildfor                                     ,	hg_pass1sclevels9pfa6 1536K, 12, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa7 1792K, 12, 4
buildfor                        + K8_32      ,	hg_pass1sclevels9pfa7 1792K, 12, 1
buildfor         + P4    + P4TP              ,	hg_pass1sclevels9pfa8 2048K, 12, 4
buildfor                                     ,	hg_pass1sclevels9pfa8 2048K, 12, 2
buildfor                                     ,	hg_pass1sclevels9pfa8 2048K, 12, 1
buildfor                 + P4TP              ,	hg_pass1sclevels10pfa5 2560K, 12, 4
buildfor                                     ,	hg_pass1sclevels10pfa5 2560K, 12, 1
buildfor CORE_64 + P4_32                     ,	hg_pass1sclevels10pfa6 3072K, 12, 4
buildfor                                     ,	hg_pass1sclevels10pfa6 3072K, 12, 1
buildfor                                     ,	hg_pass1sclevels10pfa7 3584K, 12, 4
buildfor                                     ,	hg_pass1sclevels10pfa7 3584K, 12, 2
buildfor                                     ,	hg_pass1sclevels10pfa7 3584K, 12, 1
buildfor         + P4_32                     ,	hg_pass1sclevels10pfa8 4096K, 12, 4
buildfor                                     ,	hg_pass1sclevels10pfa8 4096K, 12, 2
buildfor                                     ,	hg_pass1sclevels10pfa8 4096K, 12, 1
buildfor         + P4_32        + K8_32      ,	hg_pass1sclevels11pfa5 5M, 12, 4
buildfor                                     ,	hg_pass1sclevels11pfa5 5M, 12, 1
buildfor                                     ,	hg_pass1sclevels11pfa6 6M, 12, 4
buildfor                                     ,	hg_pass1sclevels11pfa6 6M, 12, 2
buildfor                                     ,	hg_pass1sclevels11pfa6 6M, 12, 1
buildfor                                     ,	hg_pass1sclevels11pfa7 7M, 12, 4
buildfor                                     ,	hg_pass1sclevels11pfa7 7M, 12, 2
buildfor                                     ,	hg_pass1sclevels11pfa7 7M, 12, 1
buildfor                                     ,	hg_pass1sclevels11pfa8 8M, 12, 4
buildfor                                     ,	hg_pass1sclevels11pfa8 8M, 12, 2
buildfor                                     ,	hg_pass1sclevels11pfa8 8M, 12, 1
   
buildfor                                     ,	hg_pass1sclevels8complex3 768K, 12, 2
buildfor                                     ,	hg_pass1sclevels8complex4 1024K, 12, 4
buildfor                                     ,	hg_pass1sclevels8complex4 1024K, 12, 1
buildfor                                     ,	hg_pass1sclevels9complex3 1536K, 12, 4
buildfor                                     ,	hg_pass1sclevels9complex3 1536K, 12, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9complex4 2048K, 12, 4
buildfor                                     ,	hg_pass1sclevels9complex4 2048K, 12, 2
buildfor                                     ,	hg_pass1sclevels9complex4 2048K, 12, 1
buildfor                                     ,	hg_pass1sclevels10complex3 3072K, 12, 4
buildfor                                     ,	hg_pass1sclevels10complex3 3072K, 12, 2
buildfor                                     ,	hg_pass1sclevels10complex3 3072K, 12, 1
buildfor                                     ,	hg_pass1sclevels10complex4 4096K, 12, 4
buildfor                                     ,	hg_pass1sclevels10complex4 4096K, 12, 2
buildfor                                     ,	hg_pass1sclevels10complex4 4096K, 12, 1
buildfor                                     ,	hg_pass1sclevels11complex3 6M, 12, 4
buildfor                                     ,	hg_pass1sclevels11complex3 6M, 12, 2
buildfor                                     ,	hg_pass1sclevels11complex3 6M, 12, 1
buildfor                                     ,	hg_pass1sclevels11complex4 8M, 12, 4
buildfor                                     ,	hg_pass1sclevels11complex4 8M, 12, 2
buildfor                                     ,	hg_pass1sclevels11complex4 8M, 12, 1

buildfor                                     ,	hg_pass1sclevels9pfa5 2560K, 13, 1
buildfor                                     ,	hg_pass1sclevels9pfa6 3072K, 13, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa7 3584K, 13, 1
buildfor                 + P4TP              ,	hg_pass1sclevels9pfa8 4096K, 13, 1
buildfor                                     ,	hg_pass1sclevels11pfa5 10M, 13, 4
buildfor                                     ,	hg_pass1sclevels11pfa6 12M, 13, 4
buildfor                                     ,	hg_pass1sclevels11pfa7 14M, 13, 4
buildfor                                     ,	hg_pass1sclevels11pfa8 16M, 13, 4
buildfor                                     ,	hg_pass1sclevels12pfa5 20M, 13, 2
buildfor                                     ,	hg_pass1sclevels12pfa6 24M, 13, 2
buildfor                                     ,	hg_pass1sclevels12pfa7 28M, 13, 2
buildfor                                     ,	hg_pass1sclevels12pfa8 32M, 13, 2

buildfor                                     ,	hg_pass1sclevels9complex3 3072K, 13, 1
buildfor                                     ,	hg_pass1sclevels9complex4 4096K, 13, 1
buildfor                                     ,	hg_pass1sclevels11complex3 12M, 13, 4
buildfor                                     ,	hg_pass1sclevels11complex4 16M, 13, 4
buildfor                                     ,	hg_pass1sclevels12complex3 24M, 13, 2
buildfor                                     ,	hg_pass1sclevels12complex4 32M, 13, 2

ENDIF

_TEXT	ENDS
END
