; Copyright 2009-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the radix-4/8 DJB FFT with delayed sin/cos multiplies.
; This is identical to the radix-4/8 DJB split premultiplier FFT
; except that the first levels use common sin/cos data and the last levels use
; extra complex multiplies to repair the introduced error.  Saves lots of sin/cos memory.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

xfft_type TEXTEQU <r4delay>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE r4.mac
INCLUDE r4delaypass1sc.mac
INCLUDE r4delaypass2.mac

EXTRN	xgw_carries:PROC

_TEXT SEGMENT

;; Generate pass 2 routines

buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 8
buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 768
buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 10
buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 1280
buildfor CORE    + P4    + P4TP + K8    + K10   ,	xpass2gen 1536
buildfor CORE    + P4    + P4TP         + K10_32,	xpass2gen 11
buildfor CORE            + P4TP         + K10_32,	xpass2gen 2560

;; Routines for many FFT sizes

	; The 8 levels variants

buildfor                                        ,	r4delay_pass1sc128 32K, 8, 4
buildfor CORE_64 + P4_32        + K8            ,	r4delay_pass1sc128 32K, 8, 2
buildfor CORE_64 + P4                           ,	r4delay_pass1sc256 64K, 8, 4
buildfor CORE_32                        + K10   ,	r4delay_pass1sc256 64K, 8, 2
buildfor CORE_64                + K8            ,	r4delay_pass1sc256 64K, 8, 1
buildfor CORE_32 + P4_64                        ,	r4delay_pass1sc320 80K, 8, 4
buildfor                 + P4TP + K8_32 + K10_64,	r4delay_pass1sc320 80K, 8, 2
buildfor CORE                   + K8_64         ,	r4delay_pass1sc320 80K, 8, 1
buildfor CORE                                   ,	r4delay_pass1sc384 96K, 8, 4
buildfor CORE                                   ,	r4delay_pass1sc448 112K, 8, 4
buildfor                 + P4TP + K8            ,	r4delay_pass1sc448 112K, 8, 2
buildfor                                        ,	r4delay_pass1sc448 112K, 8, 1
buildfor CORE                                   ,	r4delay_pass1sc512 128K, 8, 4
buildfor                                        ,	r4delay_pass1sc512 128K, 8, 2
buildfor                                        ,	r4delay_pass1sc512 128K, 8, 1
buildfor CORE    + P4                           ,	r4delay_pass1sc640 160K, 8, 4
buildfor         + P4_64                        ,	r4delay_pass1sc768 192K, 8, 4
buildfor CORE    + P4                   + K10_32,	r4delay_pass1sc896 224K, 8, 4
buildfor                                + K10_64,	r4delay_pass1sc896 224K, 8, 2
buildfor                                        ,	r4delay_pass1sc896 224K, 8, 1
buildfor CORE                                   ,	r4delay_pass1sc1024 256K, 8, 4

	; The 10 levels variants (768, 1024, 1280)

buildfor CORE_64         + P4TP         + K10   ,	r4delay_pass1sc128 96K, 768, 4
buildfor CORE                   + K8_64         ,	r4delay_pass1sc128 96K, 768, 2
buildfor                        + K8_32         ,	r4delay_pass1sc128 96K, 768, 1
buildfor CORE_64 + P4_32 + P4TP         + K10   ,	r4delay_pass1sc128 128K, 10, 4
buildfor CORE                   + K8_64         ,	r4delay_pass1sc128 128K, 10, 2
buildfor                        + K8_32         ,	r4delay_pass1sc128 128K, 10, 1
buildfor                                + K10   ,	r4delay_pass1sc128 160K, 1280, 4
buildfor CORE                   + K8            ,	r4delay_pass1sc128 160K, 1280, 2
buildfor                                        ,	r4delay_pass1sc128 160K, 1280, 1

buildfor CORE_64 + P4_32                + K10   ,	r4delay_pass1sc256 192K, 768, 4
buildfor                                        ,	r4delay_pass1sc256 192K, 768, 2
buildfor                        + K8            ,	r4delay_pass1sc256 192K, 768, 1
buildfor         + P4    + P4TP                 ,	r4delay_pass1sc256 256K, 10, 4
buildfor CORE                   + K8            ,	r4delay_pass1sc256 256K, 10, 2
buildfor                                + K10   ,	r4delay_pass1sc256 256K, 10, 1
buildfor         + P4_64                        ,	r4delay_pass1sc256 320K, 1280, 4
buildfor                                        ,	r4delay_pass1sc256 320K, 1280, 2
buildfor                        + K8_64         ,	r4delay_pass1sc256 320K, 1280, 1

buildfor CORE    + P4                   + K10_64,	r4delay_pass1sc320 240K, 768, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc320 240K, 768, 2
buildfor                        + K8            ,	r4delay_pass1sc320 240K, 768, 1
buildfor CORE    + P4                           ,	r4delay_pass1sc320 320K, 10, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc320 320K, 10, 2
buildfor CORE                           + K10   ,	r4delay_pass1sc320 320K, 10, 1
buildfor CORE_64 + P4    + P4TP                 ,	r4delay_pass1sc320 400K, 1280, 4
buildfor CORE            + P4TP                 ,	r4delay_pass1sc320 400K, 1280, 2
buildfor CORE_32                + K8    + K10   ,	r4delay_pass1sc320 400K, 1280, 1

buildfor CORE    + P4                   + K10   ,	r4delay_pass1sc384 288K, 768, 4
buildfor         + P4_32 + P4TP + K8            ,	r4delay_pass1sc384 288K, 768, 2
buildfor                                        ,	r4delay_pass1sc384 288K, 768, 1
buildfor CORE_64                                ,	r4delay_pass1sc384 384K, 10, 4
buildfor                                + K10_32,	r4delay_pass1sc384 384K, 10, 2
buildfor CORE_64 + P4_32        + K8    + K10_64,	r4delay_pass1sc384 384K, 10, 1
buildfor                                        ,	r4delay_pass1sc384 480K, 1280, 4
buildfor                                        ,	r4delay_pass1sc384 480K, 1280, 2
buildfor CORE                   + K8    + K10_64,	r4delay_pass1sc384 480K, 1280, 1

buildfor CORE    + P4                   + K10   ,	r4delay_pass1sc448 336K, 768, 4
buildfor                 + P4TP + K8            ,	r4delay_pass1sc448 336K, 768, 2
buildfor                                        ,	r4delay_pass1sc448 336K, 768, 1
buildfor CORE    + P4                   + K10   ,	r4delay_pass1sc448 448K, 10, 4
buildfor CORE            + P4TP + K8_64         ,	r4delay_pass1sc448 448K, 10, 2
buildfor CORE_64                + K8_32         ,	r4delay_pass1sc448 448K, 10, 1
buildfor CORE    + P4    + P4TP         + K10_64,	r4delay_pass1sc448 560K, 1280, 4
buildfor CORE_32         + P4TP                 ,	r4delay_pass1sc448 560K, 1280, 2

buildfor         + P4_32                        ,	r4delay_pass1sc512 384K, 768, 4
buildfor                                        ,	r4delay_pass1sc512 384K, 768, 2
buildfor                                        ,	r4delay_pass1sc512 384K, 768, 1
buildfor CORE                           + K10   ,	r4delay_pass1sc512 512K, 10, 4
buildfor CORE                   + K8            ,	r4delay_pass1sc512 512K, 10, 2
buildfor                                        ,	r4delay_pass1sc512 512K, 10, 1
buildfor                        + K8_64         ,	r4delay_pass1sc512 640K, 1280, 4

buildfor CORE    + P4_32                        ,	r4delay_pass1sc640 480K, 768, 4
buildfor CORE    + P4                   + K10_64,	r4delay_pass1sc640 640K, 10, 4
buildfor                                        ,	r4delay_pass1sc640 640K, 10, 2
buildfor CORE_32                                ,	r4delay_pass1sc640 640K, 10, 1
buildfor CORE    + P4                   + K10_64,	r4delay_pass1sc640 800K, 1280, 4
buildfor                                        ,	r4delay_pass1sc640 800K, 1280, 2
buildfor CORE                   + K8            ,	r4delay_pass1sc640 800K, 1280, 1

buildfor CORE    + P4_32                        ,	r4delay_pass1sc768 576K, 768, 4
buildfor                                        ,	r4delay_pass1sc768 576K, 768, 2
buildfor                                        ,	r4delay_pass1sc768 576K, 768, 1
buildfor CORE_32 + P4_32                        ,	r4delay_pass1sc768 768K, 10, 4
buildfor                                        ,	r4delay_pass1sc768 768K, 10, 2
buildfor                                        ,	r4delay_pass1sc768 768K, 10, 1
buildfor                                        ,	r4delay_pass1sc768 960K, 1280, 4

buildfor CORE    + P4                           ,	r4delay_pass1sc896 672K, 768, 4
buildfor                                        ,	r4delay_pass1sc896 672K, 768, 2
buildfor                                        ,	r4delay_pass1sc896 672K, 768, 1
buildfor CORE    + P4                   + K10_64,	r4delay_pass1sc896 896K, 10, 4
buildfor                                        ,	r4delay_pass1sc896 896K, 10, 2
buildfor                        + K8_64         ,	r4delay_pass1sc896 896K, 10, 1
buildfor                                + K10_64,	r4delay_pass1sc896 1120K, 1280, 4

buildfor CORE_64                                ,	r4delay_pass1sc1024 768K, 768, 4
buildfor                                        ,	r4delay_pass1sc1024 768K, 768, 2
buildfor                                        ,	r4delay_pass1sc1024 768K, 768, 1
buildfor CORE                                   ,	r4delay_pass1sc1024 1M, 10, 4
buildfor                                        ,	r4delay_pass1sc1024 1280K, 1280, 4

buildfor CORE                                   ,	r4delay_pass1sc1280 960K, 768, 4
buildfor                                        ,	r4delay_pass1sc1280 1280K, 10, 4
buildfor                                        ,	r4delay_pass1sc1280 1600K, 1280, 4

buildfor                                        ,	r4delay_pass1sc1536 1536K, 10, 4

buildfor CORE_64                                ,	r4delay_pass1sc1792 1792K, 10, 4

buildfor                                        ,	r4delay_pass1sc2048 2M, 10, 4

buildfor                                        ,	r4delay_pass1sc2560 2560K, 10, 4

buildfor                                        ,	r4delay_pass1sc3072 3M, 10, 4

buildfor CORE_32                                ,	r4delay_pass1sc3584 3584K, 10, 4

buildfor                                        ,	r4delay_pass1sc4096 4M, 10, 4
buildfor                                        ,	r4delay_pass1sc4096 4M, 10, 2

	; The 11 levels variants (1536, 2048, 2560)

buildfor CORE_64 + P4    + P4TP                 ,	r4delay_pass1sc256 384K, 1536, 4
buildfor CORE_64 + P4    + P4TP                 ,	r4delay_pass1sc256 512K, 11, 4
buildfor                                        ,	r4delay_pass1sc256 640K, 2560, 4

buildfor CORE                                   ,	r4delay_pass1sc320 480K, 1536, 4
buildfor                 + P4TP         + K10_32,	r4delay_pass1sc320 480K, 1536, 2
buildfor         + P4_32                        ,	r4delay_pass1sc320 640K, 11, 4
buildfor CORE_64         + P4TP         + K10_32,	r4delay_pass1sc320 640K, 11, 2
buildfor CORE            + P4TP         + K10_32,	r4delay_pass1sc320 800K, 2560, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc320 800K, 2560, 2

buildfor CORE_32         + P4TP         + K10_32,	r4delay_pass1sc384 576K, 1536, 4
buildfor CORE            + P4TP         + K10_64,	r4delay_pass1sc384 576K, 1536, 2
buildfor                                        ,	r4delay_pass1sc384 768K, 11, 4
buildfor                                        ,	r4delay_pass1sc384 960K, 2560, 4

buildfor CORE                           + K10_32,	r4delay_pass1sc448 672K, 1536, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc448 672K, 1536, 2
buildfor CORE                   + K8    + K10_64,	r4delay_pass1sc448 672K, 1536, 1
buildfor CORE_32                        + K10_32,	r4delay_pass1sc448 896K, 11, 4
buildfor CORE            + P4TP                 ,	r4delay_pass1sc448 896K, 11, 2
buildfor CORE            + P4TP                 ,	r4delay_pass1sc448 1120K, 2560, 4
buildfor CORE            + P4TP                 ,	r4delay_pass1sc448 1120K, 2560, 2
buildfor                                        ,	r4delay_pass1sc448 1120K, 2560, 1

buildfor                                + K10   ,	r4delay_pass1sc512 768K, 1536, 4
buildfor         + P4_64                + K10_32,	r4delay_pass1sc512 1M, 11, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc512 1M, 11, 2
buildfor                                        ,	r4delay_pass1sc512 1280K, 2560, 4
buildfor                                        ,	r4delay_pass1sc512 1280K, 2560, 2

buildfor CORE    + P4_32                + K10   ,	r4delay_pass1sc640 960K, 1536, 4
buildfor CORE_32 + P4_64                        ,	r4delay_pass1sc640 1280K, 11, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc640 1280K, 11, 2
buildfor CORE                                   ,	r4delay_pass1sc640 1600K, 2560, 4

buildfor                                        ,	r4delay_pass1sc768 1152K, 1536, 4
buildfor CORE_32                                ,	r4delay_pass1sc768 1536K, 11, 4
buildfor                                        ,	r4delay_pass1sc768 1536K, 11, 2
buildfor                                        ,	r4delay_pass1sc768 1920K, 2560, 4

buildfor CORE_64 + P4_32                        ,	r4delay_pass1sc896 1344K, 1536, 4
buildfor CORE_64 + P4_32                        ,	r4delay_pass1sc896 1792K, 11, 4
buildfor CORE                                   ,	r4delay_pass1sc896 2240K, 2560, 4

buildfor                                        ,	r4delay_pass1sc1024 1536K, 1536, 4
buildfor CORE                                   ,	r4delay_pass1sc1024 2M, 11, 4
buildfor CORE_64                                ,	r4delay_pass1sc1024 2560K, 2560, 4

buildfor                                + K10_32,	r4delay_pass1sc1280 1920K, 1536, 4
buildfor                                        ,	r4delay_pass1sc1280 2560K, 11, 4
buildfor                                        ,	r4delay_pass1sc1280 3200K, 2560, 4

buildfor                                        ,	r4delay_pass1sc1536 3M, 11, 4

buildfor                                        ,	r4delay_pass1sc1792 3584K, 11, 4

buildfor                                        ,	r4delay_pass1sc2048 4M, 11, 4

buildfor                                        ,	r4delay_pass1sc2560 5M, 11, 4

buildfor                                        ,	r4delay_pass1sc3072 6M, 11, 4

buildfor                                        ,	r4delay_pass1sc3584 7M, 11, 4

buildfor                                        ,	r4delay_pass1sc4096 8M, 11, 4
buildfor                                        ,	r4delay_pass1sc4096 8M, 11, 2

	; The all-complex 8 levels variants

buildfor                                        ,	r4delay_pass1sc128ac 32K, 8, 4
buildfor                                        ,	r4delay_pass1sc128ac 32K, 8, 2
buildfor                                        ,	r4delay_pass1sc128ac 32K, 8, 1
buildfor                                        ,	r4delay_pass1sc256ac 64K, 8, 4
buildfor                                        ,	r4delay_pass1sc256ac 64K, 8, 2
buildfor                                        ,	r4delay_pass1sc256ac 64K, 8, 1
buildfor                                        ,	r4delay_pass1sc384ac 96K, 8, 4
buildfor                                        ,	r4delay_pass1sc512ac 128K, 8, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc640ac 160K, 8, 4
buildfor                                        ,	r4delay_pass1sc768ac 192K, 8, 4
buildfor                                        ,	r4delay_pass1sc1024ac 256K, 8, 4
buildfor                                        ,	r4delay_pass1sc1280ac 320K, 8, 4

	; The all-complex 10 levels variants (768, 1024, 1280)

buildfor                 + P4TP                 ,	r4delay_pass1sc128ac 96K, 768, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc128ac 128K, 10, 4
buildfor                                        ,	r4delay_pass1sc128ac 160K, 1280, 4

buildfor                                        ,	r4delay_pass1sc256ac 192K, 768, 4
buildfor                                        ,	r4delay_pass1sc256ac 256K, 10, 4
buildfor                                        ,	r4delay_pass1sc256ac 320K, 1280, 4

buildfor                 + P4TP                 ,	r4delay_pass1sc384ac 288K, 768, 4
buildfor                                        ,	r4delay_pass1sc384ac 384K, 10, 4
buildfor                                        ,	r4delay_pass1sc384ac 480K, 1280, 4

buildfor                                        ,	r4delay_pass1sc512ac 384K, 768, 4
buildfor                                        ,	r4delay_pass1sc512ac 512K, 10, 4
buildfor                                        ,	r4delay_pass1sc512ac 640K, 1280, 4
buildfor                                        ,	r4delay_pass1sc512ac 640K, 1280, 2

buildfor                                        ,	r4delay_pass1sc640ac 480K, 768, 4
buildfor                                        ,	r4delay_pass1sc640ac 640K, 10, 4
buildfor                                        ,	r4delay_pass1sc640ac 800K, 1280, 4
buildfor                                        ,	r4delay_pass1sc640ac 800K, 1280, 2

buildfor                                        ,	r4delay_pass1sc768ac 576K, 768, 4
buildfor                                        ,	r4delay_pass1sc768ac 768K, 10, 4
buildfor                                        ,	r4delay_pass1sc768ac 960K, 1280, 4

buildfor                                        ,	r4delay_pass1sc1024ac 768K, 768, 4
buildfor                                        ,	r4delay_pass1sc1024ac 768K, 768, 2
buildfor                                        ,	r4delay_pass1sc1024ac 1M, 10, 4
buildfor                                        ,	r4delay_pass1sc1024ac 1280K, 1280, 4

buildfor                                        ,	r4delay_pass1sc1280ac 960K, 768, 4
buildfor                                        ,	r4delay_pass1sc1280ac 1280K, 10, 4
buildfor                                        ,	r4delay_pass1sc1280ac 1600K, 1280, 4

buildfor                                        ,	r4delay_pass1sc1536ac 1536K, 10, 4

buildfor                                        ,	r4delay_pass1sc2048ac 2M, 10, 4

buildfor                                        ,	r4delay_pass1sc2560ac 2560K, 10, 4

buildfor                                        ,	r4delay_pass1sc3072ac 3M, 10, 4

buildfor                                        ,	r4delay_pass1sc4096ac 4M, 10, 4

buildfor                                        ,	r4delay_pass1sc5120ac 5M, 10, 4

	; The all-complex 11 levels variants (1536, 2048, 2560)

buildfor                                        ,	r4delay_pass1sc256ac 384K, 1536, 4
buildfor                                        ,	r4delay_pass1sc256ac 512K, 11, 4
buildfor                                        ,	r4delay_pass1sc256ac 640K, 2560, 4

buildfor                 + P4TP                 ,	r4delay_pass1sc384ac 576K, 1536, 4
buildfor                                        ,	r4delay_pass1sc384ac 768K, 11, 4
buildfor                                        ,	r4delay_pass1sc384ac 768K, 11, 2
buildfor                                        ,	r4delay_pass1sc384ac 960K, 2560, 4
buildfor                                        ,	r4delay_pass1sc384ac 960K, 2560, 2

buildfor                                        ,	r4delay_pass1sc512ac 768K, 1536, 4
buildfor                                        ,	r4delay_pass1sc512ac 768K, 1536, 2
buildfor                                        ,	r4delay_pass1sc512ac 1M, 11, 4
buildfor                                        ,	r4delay_pass1sc512ac 1M, 11, 2
buildfor                                        ,	r4delay_pass1sc512ac 1280K, 2560, 4
buildfor                                        ,	r4delay_pass1sc512ac 1280K, 2560, 2

buildfor                                        ,	r4delay_pass1sc640ac 960K, 1536, 4
buildfor                                        ,	r4delay_pass1sc640ac 960K, 1536, 2
buildfor                                        ,	r4delay_pass1sc640ac 1280K, 11, 4
buildfor                                        ,	r4delay_pass1sc640ac 1280K, 11, 2
buildfor                                        ,	r4delay_pass1sc640ac 1600K, 2560, 4
buildfor                                        ,	r4delay_pass1sc640ac 1600K, 2560, 2

buildfor                                        ,	r4delay_pass1sc768ac 1152K, 1536, 4
buildfor                                        ,	r4delay_pass1sc768ac 1536K, 11, 4
buildfor                                        ,	r4delay_pass1sc768ac 1536K, 11, 2
buildfor                                        ,	r4delay_pass1sc768ac 1920K, 2560, 4

buildfor                                        ,	r4delay_pass1sc1024ac 1536K, 1536, 4
buildfor                                        ,	r4delay_pass1sc1024ac 2M, 11, 4
buildfor                                        ,	r4delay_pass1sc1024ac 2560K, 2560, 4

buildfor                                        ,	r4delay_pass1sc1280ac 1920K, 1536, 4
buildfor                                        ,	r4delay_pass1sc1280ac 2560K, 11, 4
buildfor                                        ,	r4delay_pass1sc1280ac 3200K, 2560, 4

;; Some small FFTs fit completely in a large L2 cache and thus do
;; not need to be prefetched.

IF TLB_PRIMING EQ 0

PREFETCHING = 0

buildfor CORE    + P4_32                + K10   ,	xpass2gen 8
buildfor CORE_64                                ,	xpass2gen 768
buildfor CORE_64                                ,	xpass2gen 10
buildfor                                        ,	xpass2gen 1280

buildfor CORE_32                        + K10   ,	r4delay_pass1sc128 32K, 8, 4
buildfor                                        ,	r4delay_pass1sc128 32K, 8, 2
buildfor CORE_64 + P4_32                        ,	r4delay_pass1sc256 64K, 8, 4
buildfor                                        ,	r4delay_pass1sc256 64K, 8, 2
buildfor                                        ,	r4delay_pass1sc256 64K, 8, 1
buildfor CORE                                   ,	r4delay_pass1sc320 80K, 8, 4
buildfor                                        ,	r4delay_pass1sc384 96K, 8, 4
buildfor CORE_64                                ,	r4delay_pass1sc448 112K, 8, 4
buildfor CORE_32                                ,	r4delay_pass1sc512 128K, 8, 4
buildfor CORE_64                                ,	r4delay_pass1sc128 96K, 768, 4
buildfor CORE_64                                ,	r4delay_pass1sc128 128K, 10, 4

buildfor                                        ,	r4delay_pass1sc128ac 32K, 8, 4
buildfor                                        ,	r4delay_pass1sc128ac 32K, 8, 2
buildfor                                        ,	r4delay_pass1sc128ac 32K, 8, 1
buildfor                                        ,	r4delay_pass1sc256ac 64K, 8, 4
buildfor                                        ,	r4delay_pass1sc256ac 64K, 8, 2
buildfor                                        ,	r4delay_pass1sc256ac 64K, 8, 1
buildfor                                        ,	r4delay_pass1sc384ac 96K, 8, 4
buildfor                                        ,	r4delay_pass1sc512ac 128K, 8, 4
buildfor                                        ,	r4delay_pass1sc128ac 96K, 768, 4
buildfor                                        ,	r4delay_pass1sc128ac 128K, 10, 4

ENDIF

_TEXT	ENDS
END
