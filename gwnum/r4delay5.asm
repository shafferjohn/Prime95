; Copyright 2009-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These FFTs were split out of r4delay.asm because of MASM limitations.
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

;; Generate pass 2 routines optimized for this architecture

buildfor CORE    + P4    + P4TP         + K10   ,	xpass2gen 4608
buildfor CORE    + P4    + P4TP         + K10   ,	xpass2gen 6144
buildfor CORE    + P4                   + K10_64,	xpass2gen 7680
buildfor CORE    + P4    + P4TP         + K10   ,	xpass2gen 13
buildfor CORE    + P4                   + K10   ,	xpass2gen 10240
buildfor CORE    + P4                           ,	xpass2gen 12800

;; Routines for many FFT sizes

	; The 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

buildfor CORE_64 + P4    + P4TP         + K10_32,	r4delay_pass1sc256 1152K, 4608, 4
buildfor CORE_32                                ,	r4delay_pass1sc256 1152K, 4608, 2
buildfor                                        ,	r4delay_pass1sc256 1536K, 6144, 4
buildfor                                        ,	r4delay_pass1sc256 1920K, 7680, 4
buildfor                                        ,	r4delay_pass1sc256 2M, 13, 4
buildfor                                        ,	r4delay_pass1sc256 2560K, 10240, 4
buildfor                                        ,	r4delay_pass1sc256 3200K, 12800, 4

buildfor CORE_64 + P4                   + K10_32,	r4delay_pass1sc320 1440K, 4608, 4
buildfor CORE_32         + P4TP                 ,	r4delay_pass1sc320 1440K, 4608, 2
buildfor         + P4_32                        ,	r4delay_pass1sc320 1920K, 6144, 4
buildfor         + P4_32                        ,	r4delay_pass1sc320 2400K, 7680, 4
buildfor                                        ,	r4delay_pass1sc320 2560K, 13, 4
buildfor                                        ,	r4delay_pass1sc320 3200K, 10240, 4
buildfor                                        ,	r4delay_pass1sc320 4000K, 12800, 4

buildfor CORE_64                        + K10   ,	r4delay_pass1sc384 1728K, 4608, 4
buildfor                                        ,	r4delay_pass1sc384 2304K, 6144, 4
buildfor                                        ,	r4delay_pass1sc384 2880K, 7680, 4
buildfor                                        ,	r4delay_pass1sc384 3M, 13, 4
buildfor                                        ,	r4delay_pass1sc384 3840K, 10240, 4
buildfor                                        ,	r4delay_pass1sc384 4800K, 12800, 4

buildfor                                        ,	r4delay_pass1sc448 2016K, 4608, 4
buildfor         + P4_32                + K10_32,	r4delay_pass1sc448 2688K, 6144, 4
buildfor                                        ,	r4delay_pass1sc448 2688K, 6144, 2
buildfor         + P4_32                        ,	r4delay_pass1sc448 3360K, 7680, 4
buildfor                                        ,	r4delay_pass1sc448 3584K, 13, 4
buildfor                                        ,	r4delay_pass1sc448 4480K, 10240, 4
buildfor         + P4_32                        ,	r4delay_pass1sc448 5600K, 12800, 4

buildfor         + P4_32                + K10   ,	r4delay_pass1sc512 2304K, 4608, 4
buildfor CORE_32         + P4TP                 ,	r4delay_pass1sc512 2304K, 4608, 2
buildfor                                + K10_32,	r4delay_pass1sc512 3M, 6144, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc512 3M, 6144, 2
buildfor                                        ,	r4delay_pass1sc512 3840K, 7680, 4
buildfor                                        ,	r4delay_pass1sc512 4M, 13, 4
buildfor                                        ,	r4delay_pass1sc512 4M, 13, 2
buildfor                                        ,	r4delay_pass1sc512 4M, 13, 1
buildfor                                        ,	r4delay_pass1sc512 5M, 10240, 4
buildfor                                        ,	r4delay_pass1sc512 6400K, 12800, 4

buildfor         + P4                   + K10   ,	r4delay_pass1sc640 2880K, 4608, 4
buildfor CORE            + P4TP                 ,	r4delay_pass1sc640 2880K, 4608, 2
buildfor         + P4    + P4TP         + K10_32,	r4delay_pass1sc640 3840K, 6144, 4
buildfor CORE    + P4                           ,	r4delay_pass1sc640 4800K, 7680, 4
buildfor         + P4                           ,	r4delay_pass1sc640 5M, 13, 4
buildfor CORE_64                                ,	r4delay_pass1sc640 6400K, 10240, 4
buildfor CORE_64                                ,	r4delay_pass1sc640 8000K, 12800, 4

buildfor CORE_64                        + K10   ,	r4delay_pass1sc768 3456K, 4608, 4
buildfor                                + K10   ,	r4delay_pass1sc768 4608K, 6144, 4
buildfor CORE_64                                ,	r4delay_pass1sc768 5760K, 7680, 4
buildfor CORE_64                                ,	r4delay_pass1sc768 6M, 13, 4
buildfor                                        ,	r4delay_pass1sc768 7680K, 10240, 4
buildfor                                        ,	r4delay_pass1sc768 9600K, 12800, 4

buildfor CORE    + P4_32                        ,	r4delay_pass1sc896 4032K, 4608, 4
buildfor CORE    + P4                   + K10   ,	r4delay_pass1sc896 5376K, 6144, 4
buildfor CORE                           + K10_64,	r4delay_pass1sc896 6720K, 7680, 4
buildfor CORE_32                                ,	r4delay_pass1sc896 7M, 13, 4
buildfor CORE                                   ,	r4delay_pass1sc896 8960K, 10240, 4
buildfor CORE                                   ,	r4delay_pass1sc896 11200K, 12800, 4

buildfor CORE    + P4                           ,	r4delay_pass1sc1024 4608K, 4608, 4
buildfor CORE_32 + P4                   + K10_32,	r4delay_pass1sc1024 6M, 6144, 4
buildfor CORE_64                                ,	r4delay_pass1sc1024 7680K, 7680, 4
buildfor CORE                           + K10   ,	r4delay_pass1sc1024 8M, 13, 4
buildfor CORE_64                                ,	r4delay_pass1sc1024 10M, 10240, 4
buildfor CORE                                   ,	r4delay_pass1sc1024 12800K, 12800, 4

buildfor                                + K10   ,	r4delay_pass1sc1280 5760K, 4608, 4
buildfor                                + K10_32,	r4delay_pass1sc1280 7680K, 6144, 4
buildfor                                        ,	r4delay_pass1sc1280 9600K, 7680, 4
buildfor                                + K10   ,	r4delay_pass1sc1280 10M, 13, 4
buildfor                                + K10   ,	r4delay_pass1sc1280 12800K, 10240, 4
buildfor                                        ,	r4delay_pass1sc1280 16000K, 12800, 4

buildfor CORE_64                        + K10_32,	r4delay_pass1sc1536 6912K, 4608, 4
buildfor                                + K10_32,	r4delay_pass1sc1536 9M, 6144, 4
buildfor                                        ,	r4delay_pass1sc1536 11520K, 7680, 4
buildfor                                        ,	r4delay_pass1sc1536 12M, 13, 4
buildfor                                        ,	r4delay_pass1sc1536 15M, 10240, 4
buildfor                                        ,	r4delay_pass1sc1536 19200K, 12800, 4

buildfor                                        ,	r4delay_pass1sc1792 8064K, 4608, 4
buildfor                                        ,	r4delay_pass1sc1792 10752K, 6144, 4
buildfor                                        ,	r4delay_pass1sc1792 13440K, 7680, 4
buildfor                                + K10   ,	r4delay_pass1sc1792 14M, 13, 4
buildfor                                        ,	r4delay_pass1sc1792 17920K, 10240, 4
buildfor                                        ,	r4delay_pass1sc1792 22400K, 12800, 4

buildfor                                        ,	r4delay_pass1sc2048 9M, 4608, 4
buildfor                                        ,	r4delay_pass1sc2048 12M, 6144, 4
buildfor                                        ,	r4delay_pass1sc2048 15M, 7680, 4
buildfor                                        ,	r4delay_pass1sc2048 16M, 13, 4
buildfor                                + K10   ,	r4delay_pass1sc2048 16M, 13, 2
buildfor                                        ,	r4delay_pass1sc2048 16M, 13, 1
buildfor                                        ,	r4delay_pass1sc2048 20M, 10240, 4
buildfor                                + K10   ,	r4delay_pass1sc2048 20M, 10240, 2
buildfor                                        ,	r4delay_pass1sc2048 25M, 12800, 4

buildfor CORE                                   ,	r4delay_pass1sc2560 11520K, 4608, 4
buildfor                                        ,	r4delay_pass1sc2560 15M, 6144, 4
buildfor                                        ,	r4delay_pass1sc2560 19200K, 7680, 4
buildfor CORE                                   ,	r4delay_pass1sc2560 20M, 13, 4
buildfor CORE_64                                ,	r4delay_pass1sc2560 25M, 10240, 4
buildfor                                        ,	r4delay_pass1sc2560 32000K, 12800, 4

buildfor CORE_64                                ,	r4delay_pass1sc3072 13824K, 4608, 4
buildfor                                        ,	r4delay_pass1sc3072 18M, 6144, 4
buildfor                                        ,	r4delay_pass1sc3072 23040K, 7680, 4
buildfor CORE                                   ,	r4delay_pass1sc3072 24M, 13, 4
buildfor                                        ,	r4delay_pass1sc3072 30M, 10240, 4
;buildfor                                       ,	r4delay_pass1sc3072 38400K, 12800, 4

buildfor                                        ,	r4delay_pass1sc3584 16128K, 4608, 4
buildfor CORE                                   ,	r4delay_pass1sc3584 21M, 6144, 4
buildfor CORE                                   ,	r4delay_pass1sc3584 26880K, 7680, 4
buildfor CORE                                   ,	r4delay_pass1sc3584 28M, 13, 4
;buildfor                                       ,	r4delay_pass1sc3584 35840K, 10240, 4
;buildfor                                       ,	r4delay_pass1sc3584 44800K, 12800, 4

buildfor CORE                                   ,	r4delay_pass1sc4096 18M, 4608, 4
buildfor                                        ,	r4delay_pass1sc4096 18M, 4608, 2
buildfor CORE                                   ,	r4delay_pass1sc4096 24M, 6144, 4
buildfor CORE_64                                ,	r4delay_pass1sc4096 30M, 7680, 4
buildfor CORE                                   ,	r4delay_pass1sc4096 32M, 13, 4
buildfor                                        ,	r4delay_pass1sc4096 32M, 13, 2
buildfor                                        ,	r4delay_pass1sc4096 32M, 13, 1
;buildfor                                       ,	r4delay_pass1sc4096 40M, 10240, 4
;buildfor                                       ,	r4delay_pass1sc4096 50M, 12800, 4

	; The all-complex 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

buildfor                                        ,	r4delay_pass1sc256ac 1152K, 4608, 4
buildfor                                        ,	r4delay_pass1sc256ac 1536K, 6144, 4
buildfor                                        ,	r4delay_pass1sc256ac 1920K, 7680, 4
buildfor                                        ,	r4delay_pass1sc256ac 2M, 13, 4
buildfor                                        ,	r4delay_pass1sc256ac 2560K, 10240, 4
buildfor                                        ,	r4delay_pass1sc256ac 3200K, 12800, 4

buildfor                                        ,	r4delay_pass1sc384ac 1728K, 4608, 4
buildfor                                        ,	r4delay_pass1sc384ac 2304K, 6144, 4
buildfor                                        ,	r4delay_pass1sc384ac 2880K, 7680, 4
buildfor                                        ,	r4delay_pass1sc384ac 3M, 13, 4
buildfor                                        ,	r4delay_pass1sc384ac 3840K, 10240, 4
buildfor                                        ,	r4delay_pass1sc384ac 4800K, 12800, 4

buildfor                                        ,	r4delay_pass1sc512ac 2304K, 4608, 4
buildfor                 + P4TP                 ,	r4delay_pass1sc512ac 2304K, 4608, 2
buildfor                 + P4TP                 ,	r4delay_pass1sc512ac 3M, 6144, 4
buildfor                                        ,	r4delay_pass1sc512ac 3M, 6144, 2
buildfor                                        ,	r4delay_pass1sc512ac 3840K, 7680, 4
buildfor                                        ,	r4delay_pass1sc512ac 4M, 13, 4
buildfor                                        ,	r4delay_pass1sc512ac 4M, 13, 2
buildfor                                        ,	r4delay_pass1sc512ac 4M, 13, 1
buildfor                                        ,	r4delay_pass1sc512ac 5M, 10240, 4
buildfor                                        ,	r4delay_pass1sc512ac 6400K, 12800, 4

buildfor                                        ,	r4delay_pass1sc640ac 2880K, 4608, 4
buildfor                                        ,	r4delay_pass1sc640ac 2880K, 4608, 2
buildfor                 + P4TP                 ,	r4delay_pass1sc640ac 3840K, 6144, 4
buildfor                                        ,	r4delay_pass1sc640ac 4800K, 7680, 4
buildfor                                        ,	r4delay_pass1sc640ac 5M, 13, 4
buildfor                                        ,	r4delay_pass1sc640ac 6400K, 10240, 4
buildfor                                        ,	r4delay_pass1sc640ac 8000K, 12800, 4

buildfor                                        ,	r4delay_pass1sc768ac 3456K, 4608, 4
buildfor                                        ,	r4delay_pass1sc768ac 4608K, 6144, 4
buildfor                                        ,	r4delay_pass1sc768ac 5760K, 7680, 4
buildfor                                        ,	r4delay_pass1sc768ac 6M, 13, 4
buildfor                                        ,	r4delay_pass1sc768ac 7680K, 10240, 4
buildfor                                        ,	r4delay_pass1sc768ac 9600K, 12800, 4

buildfor                                        ,	r4delay_pass1sc1024ac 4608K, 4608, 4
buildfor                                        ,	r4delay_pass1sc1024ac 6M, 6144, 4
buildfor                                        ,	r4delay_pass1sc1024ac 7680K, 7680, 4
buildfor                                        ,	r4delay_pass1sc1024ac 8M, 13, 4
buildfor                                        ,	r4delay_pass1sc1024ac 10M, 10240, 4
buildfor                                        ,	r4delay_pass1sc1024ac 12800K, 12800, 4

buildfor                                        ,	r4delay_pass1sc1280ac 5760K, 4608, 4
buildfor                                        ,	r4delay_pass1sc1280ac 7680K, 6144, 4
buildfor                                        ,	r4delay_pass1sc1280ac 9600K, 7680, 4
buildfor                                        ,	r4delay_pass1sc1280ac 10M, 13, 4
buildfor                                        ,	r4delay_pass1sc1280ac 12800K, 10240, 4
buildfor                                        ,	r4delay_pass1sc1280ac 16000K, 12800, 4

buildfor                                        ,	r4delay_pass1sc1536ac 6912K, 4608, 4
buildfor                                        ,	r4delay_pass1sc1536ac 9M, 6144, 4
buildfor                                        ,	r4delay_pass1sc1536ac 11520K, 7680, 4
buildfor                                        ,	r4delay_pass1sc1536ac 12M, 13, 4
buildfor                                        ,	r4delay_pass1sc1536ac 15M, 10240, 4
buildfor                                        ,	r4delay_pass1sc1536ac 19200K, 12800, 4

buildfor                                        ,	r4delay_pass1sc2048ac 9M, 4608, 4
buildfor                                        ,	r4delay_pass1sc2048ac 12M, 6144, 4
buildfor                                        ,	r4delay_pass1sc2048ac 15M, 7680, 4
buildfor                                        ,	r4delay_pass1sc2048ac 16M, 13, 4
buildfor                                        ,	r4delay_pass1sc2048ac 16M, 13, 2
buildfor                                        ,	r4delay_pass1sc2048ac 16M, 13, 1
buildfor                                        ,	r4delay_pass1sc2048ac 20M, 10240, 4
buildfor                                        ,	r4delay_pass1sc2048ac 20M, 10240, 2
buildfor                                        ,	r4delay_pass1sc2048ac 25M, 12800, 4

buildfor                                        ,	r4delay_pass1sc2560ac 11520K, 4608, 4
buildfor                                        ,	r4delay_pass1sc2560ac 15M, 6144, 4
buildfor                                        ,	r4delay_pass1sc2560ac 19200K, 7680, 4
buildfor                                        ,	r4delay_pass1sc2560ac 20M, 13, 4
buildfor                                        ,	r4delay_pass1sc2560ac 25M, 10240, 4
buildfor                                        ,	r4delay_pass1sc2560ac 32000K, 12800, 4

buildfor                                        ,	r4delay_pass1sc3072ac 13824K, 4608, 4
buildfor                                        ,	r4delay_pass1sc3072ac 18M, 6144, 4
buildfor                                        ,	r4delay_pass1sc3072ac 23040K, 7680, 4
buildfor                                        ,	r4delay_pass1sc3072ac 24M, 13, 4
buildfor                                        ,	r4delay_pass1sc3072ac 30M, 10240, 4
;buildfor                                       ,	r4delay_pass1sc3072ac 38400K, 12800, 4

buildfor                                        ,	r4delay_pass1sc4096ac 18M, 4608, 4
buildfor                                        ,	r4delay_pass1sc4096ac 24M, 6144, 4
buildfor                                        ,	r4delay_pass1sc4096ac 30M, 7680, 4
buildfor                                        ,	r4delay_pass1sc4096ac 32M, 13, 4
buildfor                                        ,	r4delay_pass1sc4096ac 32M, 13, 2
buildfor                                        ,	r4delay_pass1sc4096ac 32M, 13, 1
;buildfor                                       ,	r4delay_pass1sc4096ac 40M, 10240, 4
;buildfor                                       ,	r4delay_pass1sc4096ac 50M, 12800, 4

buildfor                                        ,	r4delay_pass1sc5120ac 23040K, 4608, 4
buildfor                                        ,	r4delay_pass1sc5120ac 30M, 6144, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 38400K, 7680, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 40M, 13, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 50M, 10240, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 64000K, 12800, 4

;; Generate pass 2 routines optimized for this architecture

buildfor CORE    + P4           + K8    + K10   ,	xpass2gen 9216
buildfor CORE                   + K8    + K10   ,	xpass2gen 12288
buildfor CORE                           + K10   ,	xpass2gen 15360
buildfor CORE                           + K10   ,	xpass2gen 14
buildfor CORE_32                                ,	xpass2gen 20480
buildfor CORE                           + K10   ,	xpass2gen 25600

;; Routines for many FFT sizes

	; The 14 levels variants (9216, 12288, 15360, 16384, 20480, 25600)

buildfor         + P4_64                        ,	r4delay_pass1sc256 2304K, 9216, 4
buildfor                                        ,	r4delay_pass1sc256 3M, 12288, 4
buildfor                                        ,	r4delay_pass1sc256 3840K, 15360, 4
buildfor                                        ,	r4delay_pass1sc256 4M, 14, 4
buildfor                                        ,	r4delay_pass1sc256 5M, 20480, 4
buildfor                                        ,	r4delay_pass1sc256 6400K, 25600, 4

buildfor CORE                                   ,	r4delay_pass1sc320 2880K, 9216, 4
buildfor CORE_32                                ,	r4delay_pass1sc320 3840K, 12288, 4
buildfor                                        ,	r4delay_pass1sc320 4800K, 15360, 4
buildfor CORE_64                                ,	r4delay_pass1sc320 5M, 14, 4
buildfor                                        ,	r4delay_pass1sc320 6400K, 20480, 4
buildfor                                        ,	r4delay_pass1sc320 8000K, 25600, 4

buildfor CORE_64                                ,	r4delay_pass1sc384 3456K, 9216, 4
buildfor                                        ,	r4delay_pass1sc384 3456K, 9216, 2
buildfor                        + K8            ,	r4delay_pass1sc384 3456K, 9216, 1
buildfor                                        ,	r4delay_pass1sc384 4608K, 12288, 4
buildfor                                        ,	r4delay_pass1sc384 5760K, 15360, 4
buildfor                                        ,	r4delay_pass1sc384 6M, 14, 4
buildfor                                        ,	r4delay_pass1sc384 7680K, 20480, 4
buildfor                                        ,	r4delay_pass1sc384 9600K, 25600, 4

buildfor                                        ,	r4delay_pass1sc448 4032K, 9216, 4
buildfor                        + K8            ,	r4delay_pass1sc448 5376K, 12288, 4
buildfor                                        ,	r4delay_pass1sc448 6720K, 15360, 4
buildfor                                        ,	r4delay_pass1sc448 7M, 14, 4
buildfor                                        ,	r4delay_pass1sc448 8960K, 20480, 4
buildfor                                        ,	r4delay_pass1sc448 11200K, 25600, 4

buildfor CORE_64                                ,	r4delay_pass1sc512 4608K, 9216, 4
buildfor                                        ,	r4delay_pass1sc512 6M, 12288, 4
buildfor                                        ,	r4delay_pass1sc512 7680K, 15360, 4
buildfor CORE_64                                ,	r4delay_pass1sc512 8M, 14, 4
buildfor                                        ,	r4delay_pass1sc512 10M, 20480, 4
buildfor                                        ,	r4delay_pass1sc512 12800K, 25600, 4

buildfor CORE_64 + P4                           ,	r4delay_pass1sc640 5760K, 9216, 4
buildfor CORE_32                                ,	r4delay_pass1sc640 7680K, 12288, 4
buildfor CORE                                   ,	r4delay_pass1sc640 9600K, 15360, 4
buildfor CORE                                   ,	r4delay_pass1sc640 10M, 14, 4
buildfor                                        ,	r4delay_pass1sc640 12800K, 20480, 4
buildfor CORE_64                                ,	r4delay_pass1sc640 16000K, 25600, 4

buildfor CORE_64                        + K10_64,	r4delay_pass1sc768 6912K, 9216, 4
buildfor                                        ,	r4delay_pass1sc768 9M, 12288, 4
buildfor                                        ,	r4delay_pass1sc768 11520K, 15360, 4
buildfor CORE_64                                ,	r4delay_pass1sc768 12M, 14, 4
buildfor                                        ,	r4delay_pass1sc768 15M, 20480, 4
buildfor CORE                                   ,	r4delay_pass1sc768 19200K, 25600, 4

buildfor                                        ,	r4delay_pass1sc896 8064K, 9216, 4
buildfor CORE_32                        + K10   ,	r4delay_pass1sc896 10752K, 12288, 4
buildfor CORE_32                        + K10   ,	r4delay_pass1sc896 13440K, 15360, 4
buildfor                                        ,	r4delay_pass1sc896 14M, 14, 4
buildfor                                        ,	r4delay_pass1sc896 17920K, 20480, 4
buildfor CORE                           + K10   ,	r4delay_pass1sc896 22400K, 25600, 4

buildfor CORE                           + K10_64,	r4delay_pass1sc1024 9M, 9216, 4
buildfor CORE_32                        + K10   ,	r4delay_pass1sc1024 12M, 12288, 4
buildfor CORE                                   ,	r4delay_pass1sc1024 15M, 15360, 4
buildfor CORE                                   ,	r4delay_pass1sc1024 16M, 14, 4
buildfor                                        ,	r4delay_pass1sc1024 20M, 20480, 4
buildfor CORE                           + K10   ,	r4delay_pass1sc1024 25M, 25600, 4

buildfor CORE_64                        + K10   ,	r4delay_pass1sc1280 11520K, 9216, 4
buildfor                                + K10   ,	r4delay_pass1sc1280 15M, 12288, 4
buildfor                                + K10   ,	r4delay_pass1sc1280 19200K, 15360, 4
buildfor CORE_32                                ,	r4delay_pass1sc1280 20M, 14, 4
buildfor CORE_32                                ,	r4delay_pass1sc1280 25M, 20480, 4
buildfor CORE_64                        + K10   ,	r4delay_pass1sc1280 32000K, 25600, 4

buildfor CORE_64                        + K10   ,	r4delay_pass1sc1536 13824K, 9216, 4
buildfor                                + K10   ,	r4delay_pass1sc1536 18M, 12288, 4
buildfor CORE_64                        + K10   ,	r4delay_pass1sc1536 23040K, 15360, 4
buildfor                                        ,	r4delay_pass1sc1536 24M, 14, 4
buildfor                                        ,	r4delay_pass1sc1536 30M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc1536 38400K, 25600, 4

buildfor                                        ,	r4delay_pass1sc1792 16128K, 9216, 4
buildfor CORE_32                        + K10   ,	r4delay_pass1sc1792 21M, 12288, 4
buildfor CORE_32                        + K10   ,	r4delay_pass1sc1792 26880K, 15360, 4
buildfor CORE                           + K10   ,	r4delay_pass1sc1792 28M, 14, 4
;buildfor                                       ,	r4delay_pass1sc1792 35840K, 20480, 4
;buildfor                                       ,	r4delay_pass1sc1792 44800K, 25600, 4

buildfor CORE                                   ,	r4delay_pass1sc2048 18M, 9216, 4
buildfor                                + K10   ,	r4delay_pass1sc2048 24M, 12288, 4
buildfor                                        ,	r4delay_pass1sc2048 30M, 15360, 4
buildfor CORE                           + K10   ,	r4delay_pass1sc2048 32M, 14, 4
;buildfor                                       ,	r4delay_pass1sc2048 40M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc2048 50M, 25600, 4

buildfor CORE                                   ,	r4delay_pass1sc2560 23040K, 9216, 4
buildfor CORE                           + K10   ,	r4delay_pass1sc2560 30M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc2560 38400K, 15360, 4
;buildfor                                       ,	r4delay_pass1sc2560 40M, 14, 4
;buildfor                                       ,	r4delay_pass1sc2560 50M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc2560 64000K, 25600, 4

buildfor CORE                           + K10   ,	r4delay_pass1sc3072 27M, 9216, 4
;buildfor                                       ,	r4delay_pass1sc3072 36M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc3072 45M, 15360, 4
;buildfor                                       ,	r4delay_pass1sc3072 48M, 14, 4
;buildfor                                       ,	r4delay_pass1sc3072 60M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc3072 75M, 25600, 4

buildfor                                        ,	r4delay_pass1sc3584 32256K, 9216, 4
;buildfor                                       ,	r4delay_pass1sc3584 42M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc3584 53760K, 15360, 4
;buildfor                                       ,	r4delay_pass1sc3584 56M, 14, 4
;buildfor                                       ,	r4delay_pass1sc3584 70M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc3584 89600K, 25600, 4

;buildfor                                       ,	r4delay_pass1sc4096 36M, 9216, 4
;buildfor                                       ,	r4delay_pass1sc4096 36M, 9216, 4
;buildfor                                       ,	r4delay_pass1sc4096 48M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc4096 60M, 15360, 4
;buildfor                                       ,	r4delay_pass1sc4096 64M, 14, 4
;buildfor                                       ,	r4delay_pass1sc4096 80M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc4096 100M, 25600, 4

	; The all-complex 14 levels variants (9216, 12288, 15360, 16384, 20480, 25600)

buildfor                                        ,	r4delay_pass1sc256ac 2304K, 9216, 4
buildfor                                        ,	r4delay_pass1sc256ac 3M, 12288, 4
buildfor                                        ,	r4delay_pass1sc256ac 3840K, 15360, 4
buildfor                                        ,	r4delay_pass1sc256ac 4M, 14, 4
buildfor                                        ,	r4delay_pass1sc256ac 5M, 20480, 4
buildfor                                        ,	r4delay_pass1sc256ac 6400K, 25600, 4

buildfor                                        ,	r4delay_pass1sc384ac 3456K, 9216, 4
buildfor                                        ,	r4delay_pass1sc384ac 4608K, 12288, 4
buildfor                                        ,	r4delay_pass1sc384ac 5760K, 15360, 4
buildfor                                        ,	r4delay_pass1sc384ac 6M, 14, 4
buildfor                                        ,	r4delay_pass1sc384ac 7680K, 20480, 4
buildfor                                        ,	r4delay_pass1sc384ac 9600K, 25600, 4

buildfor                                        ,	r4delay_pass1sc512ac 4608K, 9216, 4
buildfor                                        ,	r4delay_pass1sc512ac 6M, 12288, 4
buildfor                                        ,	r4delay_pass1sc512ac 7680K, 15360, 4
buildfor                                        ,	r4delay_pass1sc512ac 8M, 14, 4
buildfor                                        ,	r4delay_pass1sc512ac 10M, 20480, 4
buildfor                                        ,	r4delay_pass1sc512ac 12800K, 25600, 4

buildfor                                        ,	r4delay_pass1sc640ac 5760K, 9216, 4
buildfor                                        ,	r4delay_pass1sc640ac 7680K, 12288, 4
buildfor                                        ,	r4delay_pass1sc640ac 9600K, 15360, 4
buildfor                                        ,	r4delay_pass1sc640ac 10M, 14, 4
buildfor                                        ,	r4delay_pass1sc640ac 12800K, 20480, 4
buildfor                                        ,	r4delay_pass1sc640ac 16000K, 25600, 4

buildfor                                        ,	r4delay_pass1sc768ac 6912K, 9216, 4
buildfor                                        ,	r4delay_pass1sc768ac 9M, 12288, 4
buildfor                                        ,	r4delay_pass1sc768ac 11520K, 15360, 4
buildfor                                        ,	r4delay_pass1sc768ac 12M, 14, 4
buildfor                                        ,	r4delay_pass1sc768ac 15M, 20480, 4
buildfor                                        ,	r4delay_pass1sc768ac 19200K, 25600, 4

buildfor                                        ,	r4delay_pass1sc1024ac 9M, 9216, 4
buildfor                                        ,	r4delay_pass1sc1024ac 12M, 12288, 4
buildfor                                        ,	r4delay_pass1sc1024ac 15M, 15360, 4
buildfor                                        ,	r4delay_pass1sc1024ac 16M, 14, 4
buildfor                                        ,	r4delay_pass1sc1024ac 20M, 20480, 4
buildfor                                        ,	r4delay_pass1sc1024ac 25M, 25600, 4

buildfor                                        ,	r4delay_pass1sc1280ac 11520K, 9216, 4
buildfor                                        ,	r4delay_pass1sc1280ac 15M, 12288, 4
buildfor                                        ,	r4delay_pass1sc1280ac 19200K, 15360, 4
buildfor                                        ,	r4delay_pass1sc1280ac 20M, 14, 4
buildfor                                        ,	r4delay_pass1sc1280ac 25M, 20480, 4
buildfor                                        ,	r4delay_pass1sc1280ac 32000K, 25600, 4

buildfor                                        ,	r4delay_pass1sc1536ac 13824K, 9216, 4
buildfor                                        ,	r4delay_pass1sc1536ac 18M, 12288, 4
buildfor                                        ,	r4delay_pass1sc1536ac 23040K, 15360, 4
buildfor                                        ,	r4delay_pass1sc1536ac 24M, 14, 4
buildfor                                        ,	r4delay_pass1sc1536ac 30M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc1536ac 38400K, 25600, 4

buildfor                                        ,	r4delay_pass1sc2048ac 18M, 9216, 4
buildfor                                        ,	r4delay_pass1sc2048ac 24M, 12288, 4
buildfor                                        ,	r4delay_pass1sc2048ac 30M, 15360, 4
buildfor                                        ,	r4delay_pass1sc2048ac 32M, 14, 4
;buildfor                                       ,	r4delay_pass1sc2048ac 40M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc2048ac 50M, 25600, 4

buildfor                                        ,	r4delay_pass1sc2560ac 23040K, 9216, 4
buildfor                                        ,	r4delay_pass1sc2560ac 30M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc2560ac 38400K, 15360, 4
;buildfor                                       ,	r4delay_pass1sc2560ac 40M, 14, 4
;buildfor                                       ,	r4delay_pass1sc2560ac 50M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc2560ac 64000K, 25600, 4

buildfor                                        ,	r4delay_pass1sc3072ac 27M, 9216, 4
;buildfor                                       ,	r4delay_pass1sc3072ac 36M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc3072ac 45M, 15360, 4
;buildfor                                       ,	r4delay_pass1sc3072ac 48M, 14, 4
;buildfor                                       ,	r4delay_pass1sc3072ac 60M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc3072ac 75M, 25600, 4

;buildfor                                       ,	r4delay_pass1sc4096ac 36M, 9216, 4
;buildfor                                       ,	r4delay_pass1sc4096ac 48M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc4096ac 60M, 15360, 4
;buildfor                                       ,	r4delay_pass1sc4096ac 64M, 14, 4
;buildfor                                       ,	r4delay_pass1sc4096ac 80M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc4096ac 100M, 25600, 4

;buildfor                                       ,	r4delay_pass1sc5120ac 45M, 9216, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 60M, 12288, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 75M, 15360, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 80M, 14, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 100M, 20480, 4
;buildfor                                       ,	r4delay_pass1sc5120ac 125M, 25600, 4

_TEXT	ENDS
END
