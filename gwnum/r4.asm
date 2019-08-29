; Copyright 2009-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble traditional radix-4 FFTs optimized for various architectures.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <r4>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4.mac
INCLUDE r4pass1sc.mac
INCLUDE r4pass2.mac

EXTRN	xgw_carries:PROC

_TEXT SEGMENT

;; Generate pass 2 routines optmized for this architecture

buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 8
buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 768
buildfor CORE    + P4    + P4TP         + K10   ,	xpass2gen 10
buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 1280
buildfor CORE    + P4_32 + P4TP + K8    + K10   ,	xpass2gen 1536

;; Routines for many FFT sizes

buildfor                                        ,	r4_pass1sc32 8K, 8, 4
buildfor         + P4_64 + P4TP                 ,	r4_pass1sc64 16K, 8, 4
buildfor CORE_64 + P4_32 + P4TP                 ,	r4_pass1sc80 20K, 8, 4
buildfor CORE_64 + P4_32                        ,	r4_pass1sc96 24K, 8, 4
buildfor         + P4                           ,	r4_pass1sc112 28K, 8, 4
buildfor CORE_32                        + K10   ,	r4_pass1sc112 28K, 8, 2
buildfor         + P4                           ,	r4_pass1sc128 32K, 8, 4
buildfor CORE_64 + P4                           ,	r4_pass1sc224 56K, 8, 4
buildfor                                + K10   ,	r4_pass1sc224 56K, 8, 2
buildfor                        + K8            ,	r4_pass1sc224 56K, 8, 1
buildfor                                        ,	r4_pass1sc256 64K, 8, 4
buildfor                 + P4TP                 ,	r4_pass1sc32 24K, 768, 4
buildfor         + P4           + K8    + K10   ,	r4_pass1sc64 48K, 768, 4
buildfor                                        ,	r4_pass1sc80 60K, 768, 4
buildfor                        + K8_32         ,	r4_pass1sc80 60K, 768, 2
buildfor                                        ,	r4_pass1sc80 60K, 768, 1
buildfor         + P4_32 + P4TP         + K10   ,	r4_pass1sc96 72K, 768, 4
buildfor CORE                   + K8            ,	r4_pass1sc96 72K, 768, 2
buildfor                                        ,	r4_pass1sc96 72K, 768, 1
buildfor CORE    + P4    + P4TP                 ,	r4_pass1sc112 84K, 768, 4
buildfor                                + K10   ,	r4_pass1sc112 84K, 768, 2
buildfor                        + K8            ,	r4_pass1sc112 84K, 768, 1
buildfor                                        ,	r4_pass1sc128 96K, 768, 4
buildfor                                        ,	r4_pass1sc256 192K, 768, 4
buildfor                 + P4TP                 ,	r4_pass1sc32 32K, 10, 4
buildfor                                        ,	r4_pass1sc64 64K, 10, 4
buildfor         + P4                   + K10_32,	r4_pass1sc80 80K, 10, 4
buildfor                                        ,	r4_pass1sc96 96K, 10, 4
buildfor CORE_64 + P4_32                + K10   ,	r4_pass1sc112 112K, 10, 4
buildfor CORE    + P4_32 + P4TP + K8    + K10   ,	r4_pass1sc32 40K, 1280, 4
buildfor                                        ,	r4_pass1sc64 80K, 1280, 4
buildfor                                        ,	r4_pass1sc96 120K, 1280, 4
buildfor CORE    + P4_32 + P4TP                 ,	r4_pass1sc96 144K, 1536, 4
buildfor CORE_32                + K8    + K10   ,	r4_pass1sc96 144K, 1536, 2
buildfor                                        ,	r4_pass1sc96 144K, 1536, 1

buildfor                 + P4TP                 ,	r4_pass1sc32ac 8K, 8, 4
buildfor         + P4_32                        ,	r4_pass1sc64ac 16K, 8, 4
buildfor         + P4_32 + P4TP                 ,	r4_pass1sc96ac 24K, 8, 4
buildfor                                        ,	r4_pass1sc128ac 32K, 8, 4
buildfor                                        ,	r4_pass1sc256ac 64K, 8, 4
buildfor         + P4_32                        ,	r4_pass1sc32ac 24K, 768, 4
buildfor CORE_64 + P4_32        + K8    + K10   ,	r4_pass1sc64ac 48K, 768, 4
buildfor CORE_64 + P4_32 + P4TP         + K10   ,	r4_pass1sc96ac 72K, 768, 4
buildfor CORE                   + K8            ,	r4_pass1sc96ac 72K, 768, 2
buildfor                                        ,	r4_pass1sc128ac 96K, 768, 4
buildfor                                        ,	r4_pass1sc256ac 192K, 768, 4
buildfor                                        ,	r4_pass1sc32ac 32K, 10, 4
buildfor                                        ,	r4_pass1sc64ac 64K, 10, 4
buildfor                                        ,	r4_pass1sc96ac 96K, 10, 4
buildfor CORE            + P4TP + K8    + K10   ,	r4_pass1sc32ac 40K, 1280, 4
buildfor         + P4    + P4TP         + K10   ,	r4_pass1sc64ac 80K, 1280, 4
buildfor CORE                   + K8            ,	r4_pass1sc64ac 80K, 1280, 2
buildfor                                        ,	r4_pass1sc96ac 120K, 1280, 4
buildfor CORE    + P4_32 + P4TP                 ,	r4_pass1sc96ac 144K, 1536, 4
buildfor                        + K8    + K10_32,	r4_pass1sc96ac 144K, 1536, 2

;; Some small FFTs fit completely in a large L2 cache and thus do
;; not need to be prefetched.

IF TLB_PRIMING EQ 0

PREFETCHING = 0

buildfor CORE    + P4           + K8    + K10   ,	xpass2gen 8
buildfor CORE    + P4_64        + K8_64 + K10_64,	xpass2gen 768
buildfor CORE_32                                ,	xpass2gen 10
buildfor CORE    + P4                           ,	xpass2gen 1280
buildfor                                        ,	xpass2gen 1536

buildfor CORE    + P4           + K8    + K10   ,	r4_pass1sc32 8K, 8, 4
buildfor CORE    + P4           + K8    + K10   ,	r4_pass1sc64 16K, 8, 4
buildfor CORE    + P4           + K8    + K10   ,	r4_pass1sc80 20K, 8, 4
buildfor CORE    + P4_32        + K8_32 + K10   ,	r4_pass1sc96 24K, 8, 4
buildfor CORE_32                + K8_32         ,	r4_pass1sc112 28K, 8, 4
buildfor                        + K8_64         ,	r4_pass1sc112 28K, 8, 2
buildfor CORE_64                                ,	r4_pass1sc128 32K, 8, 4
buildfor         + P4_32                        ,	r4_pass1sc224 56K, 8, 4
buildfor                                        ,	r4_pass1sc256 64K, 8, 4
buildfor                        + K8_64         ,	r4_pass1sc32 24K, 768, 4
buildfor CORE                                   ,	r4_pass1sc64 48K, 768, 4
buildfor                                        ,	r4_pass1sc80 60K, 768, 4
buildfor CORE                                   ,	r4_pass1sc96 72K, 768, 4
buildfor CORE                                   ,	r4_pass1sc112 84K, 768, 4
buildfor CORE_32                                ,	r4_pass1sc128 96K, 768, 4
buildfor                                        ,	r4_pass1sc32 32K, 10, 4
buildfor                                        ,	r4_pass1sc64 64K, 10, 4
buildfor                                        ,	r4_pass1sc80 80K, 10, 4
buildfor                                        ,	r4_pass1sc96 96K, 10, 4
buildfor CORE_32                                ,	r4_pass1sc112 112K, 10, 4
buildfor CORE    + P4_64                        ,	r4_pass1sc32 40K, 1280, 4
buildfor                                        ,	r4_pass1sc64 80K, 1280, 4
buildfor                                        ,	r4_pass1sc96 120K, 1280, 4

buildfor CORE    + P4           + K8    + K10   ,	r4_pass1sc32ac 8K, 8, 4
buildfor CORE    + P4_32        + K8    + K10   ,	r4_pass1sc64ac 16K, 8, 4
buildfor CORE_64 + P4_32        + K8_32 + K10_32,	r4_pass1sc96ac 24K, 8, 4
buildfor                                        ,	r4_pass1sc128ac 32K, 8, 4
buildfor                                        ,	r4_pass1sc256ac 64K, 8, 4
buildfor CORE                   + K8_64 + K10_64,	r4_pass1sc32ac 24K, 768, 4
buildfor CORE    + P4_64                        ,	r4_pass1sc64ac 48K, 768, 4
buildfor CORE                                   ,	r4_pass1sc96ac 72K, 768, 4
buildfor                                        ,	r4_pass1sc128ac 96K, 768, 4
buildfor                                        ,	r4_pass1sc32ac 32K, 10, 4
buildfor                                        ,	r4_pass1sc64ac 64K, 10, 4
buildfor                                        ,	r4_pass1sc96ac 96K, 10, 4
buildfor CORE    + P4                           ,	r4_pass1sc32ac 40K, 1280, 4
buildfor CORE                                   ,	r4_pass1sc64ac 80K, 1280, 4
buildfor                                        ,	r4_pass1sc96ac 120K, 1280, 4

ENDIF

_TEXT	ENDS
END
