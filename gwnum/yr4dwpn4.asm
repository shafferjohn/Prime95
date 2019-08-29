; Copyright 2011-2017 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These FFTs were split out of yr4dwpn5.asm because of HJWASM limitations.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

yfft_type TEXTEQU <r4dwpn>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE ymult.mac
INCLUDE yr4.mac
INCLUDE yr4dwpnpass1sc.mac
INCLUDE yr4dwpnpass2.mac

_TEXT SEGMENT

;; Generate pass 2 routines optimized for this architecture

buildfor CORE    + FMA3_64,	ypass2gen 4608
buildfor CORE    + FMA3_64,	ypass2gen 6144
buildfor CORE    + FMA3_64,	ypass2gen 7680
buildfor CORE    + FMA3_64,	ypass2gen 13
buildfor CORE    + FMA3_64,	ypass2gen 10240
buildfor CORE    + FMA3_64,	ypass2gen 12800

;; Routines for many FFT sizes

	; The 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

IFNDEF X86_64
build421 ,			,			,			yr4dwpn_pass1sc128 576K, 4608
build421 ,			,			,			yr4dwpn_pass1sc128 768K, 6144
build421 ,			,			,			yr4dwpn_pass1sc128 960K, 7680
build421 ,			,			,			yr4dwpn_pass1sc128 1M, 13
build421 ,			,			,			yr4dwpn_pass1sc128 1280K, 10240
build421 ,			,			,			yr4dwpn_pass1sc128 1600K, 12800

build421 ,			,			,			yr4dwpn_pass1sc256 1152K, 4608
build421 ,			,			,			yr4dwpn_pass1sc256 1536K, 6144
build421 ,			,			,			yr4dwpn_pass1sc256 1920K, 7680
build421 ,			,			,			yr4dwpn_pass1sc256 2M, 13
build421 ,			,			,			yr4dwpn_pass1sc256 2560K, 10240
build421 ,			,			CORE_64,		yr4dwpn_pass1sc256 3200K, 12800

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc320 1440K, 4608
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc320 1920K, 6144
build421 CORE_64 + FMA3_64,	CORE_64 + FMA3_64,	CORE_64 + FMA3_64,	yr4dwpn_pass1sc320 2400K, 7680
build421 ,			,			,			yr4dwpn_pass1sc320 2560K, 13
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc320 3200K, 10240
build421 ,			,			,			yr4dwpn_pass1sc320 4000K, 12800

build421 ,			,			CORE_32 + FMA3_64,	yr4dwpn_pass1sc384 1728K, 4608
build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc384 2304K, 6144
build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc384 2880K, 7680
build421 ,			,			,			yr4dwpn_pass1sc384 3M, 13
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc384 3840K, 10240
build421 ,			,			,			yr4dwpn_pass1sc384 4800K, 12800

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc448 2016K, 4608
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc448 2688K, 6144
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc448 3360K, 7680
build421 ,			CORE_32,		,			yr4dwpn_pass1sc448 3584K, 13
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc448 4480K, 10240
build421 ,			,			,			yr4dwpn_pass1sc448 5600K, 12800

build421 ,			,			,			yr4dwpn_pass1sc512 2304K, 4608
build421 ,			,			,			yr4dwpn_pass1sc512 3M, 6144
build421 ,			,			,			yr4dwpn_pass1sc512 3840K, 7680
build421 ,			CORE_32,		,			yr4dwpn_pass1sc512 4M, 13
build421 ,			,			,			yr4dwpn_pass1sc512 5M, 10240
build421 ,			,			,			yr4dwpn_pass1sc512 6400K, 12800

build421 ,			,			,			yr4dwpn_pass1sc640 2880K, 4608
build421 ,			,			,			yr4dwpn_pass1sc640 3840K, 6144
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc640 4800K, 7680
build421 ,			CORE_32,		,			yr4dwpn_pass1sc640 5M, 13
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc640 6400K, 10240
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc640 8000K, 12800

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc768 3456K, 4608
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc768 4608K, 6144
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc768 5760K, 7680
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc768 6M, 13
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc768 7680K, 10240
build421 ,			,			,			yr4dwpn_pass1sc768 9600K, 12800

build421 ,			,			,			yr4dwpn_pass1sc896 4032K, 4608
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc896 5376K, 6144
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc896 6720K, 7680
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc896 7M, 13
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc896 8960K, 10240
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc896 11200K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1024 4608K, 4608
build421 ,			,			,			yr4dwpn_pass1sc1024 6M, 6144
build421 ,			,			,			yr4dwpn_pass1sc1024 7680K, 7680
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc1024 8M, 13
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc1024 10M, 10240
build421 ,			,			,			yr4dwpn_pass1sc1024 12800K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1280 5760K, 4608
build421 ,			,			,			yr4dwpn_pass1sc1280 7680K, 6144
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1280 9600K, 7680
build421 ,			,			,			yr4dwpn_pass1sc1280 10M, 13
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1280 12800K, 10240
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1280 16000K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1536 6912K, 4608
build421 ,			,			,			yr4dwpn_pass1sc1536 9M, 6144
build421 ,			,			,			yr4dwpn_pass1sc1536 11520K, 7680
build421 CORE_64 + FMA3_64,	,			,			yr4dwpn_pass1sc1536 12M, 13
build421 CORE + FMA3_64,	,			,			yr4dwpn_pass1sc1536 15M, 10240
build421 ,			,			,			yr4dwpn_pass1sc1536 19200K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1792 8064K, 4608
build421 ,			,			,			yr4dwpn_pass1sc1792 10752K, 6144
build421 ,			,			,			yr4dwpn_pass1sc1792 13440K, 7680
build421 CORE,			,			,			yr4dwpn_pass1sc1792 14M, 13
build421 CORE,			,			,			yr4dwpn_pass1sc1792 17920K, 10240
build421 CORE_64 + FMA3_64,	,			,			yr4dwpn_pass1sc1792 22400K, 12800

build421 ,			,			,			yr4dwpn_pass1sc2048 9M, 4608
build421 CORE_32,		,			,			yr4dwpn_pass1sc2048 12M, 6144
build421 ,			,			,			yr4dwpn_pass1sc2048 15M, 7680
build421 CORE_32,		FMA3_64,		,			yr4dwpn_pass1sc2048 16M, 13
build421 CORE_32,		,			,			yr4dwpn_pass1sc2048 20M, 10240
build421 ,			,			,			yr4dwpn_pass1sc2048 25M, 12800

;build421 ,			,			,			yr4dwpn_pass1sc2560 11520K, 4608
;build421 ,			,			,			yr4dwpn_pass1sc2560 15M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc2560 19200K, 7680
;build421 ,			,			,			yr4dwpn_pass1sc2560 20M, 13
;build421 ,			,			,			yr4dwpn_pass1sc2560 25M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc2560 32000K, 12800

;build421 ,			,			,			yr4dwpn_pass1sc3072 13824K, 4608
;build421 ,			,			,			yr4dwpn_pass1sc3072 18M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc3072 23040K, 7680
;build421 ,			,			,			yr4dwpn_pass1sc3072 24M, 13
;build421 ,			,			,			yr4dwpn_pass1sc3072 30M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc3072 38400K, 12800

;build421 ,			,			,			yr4dwpn_pass1sc3584 16128K, 4608
;build421 ,			,			,			yr4dwpn_pass1sc3584 21M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc3584 26880K, 7680
;build421 ,			,			,			yr4dwpn_pass1sc3584 28M, 13
;build421 ,			,			,			yr4dwpn_pass1sc3584 35840K, 10240
;build421 ,			,			,			yr4dwpn_pass1sc3584 44800K, 12800

;build421 ,			,			,			yr4dwpn_pass1sc4096 18M, 4608
;build421 ,			,			,			yr4dwpn_pass1sc4096 24M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc4096 30M, 7680
;build421 ,			,			,			yr4dwpn_pass1sc4096 32M, 13
;build421 ,			,			,			yr4dwpn_pass1sc4096 40M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc4096 50M, 12800
ENDIF

	; The all-complex 13 levels variants (4608, 6144, 7680, 8192, 10240, 12800)

IFNDEF X86_64
build421 ,			,			,			yr4dwpn_pass1sc128ac 576K, 4608
build421 ,			,			,			yr4dwpn_pass1sc128ac 768K, 6144
build421 ,			,			,			yr4dwpn_pass1sc128ac 960K, 7680
build421 ,			,			,			yr4dwpn_pass1sc128ac 1M, 13
build421 ,			,			,			yr4dwpn_pass1sc128ac 1280K, 10240
build421 ,			,			,			yr4dwpn_pass1sc128ac 1600K, 12800

build421 ,			,			,			yr4dwpn_pass1sc256ac 1152K, 4608
build421 ,			,			,			yr4dwpn_pass1sc256ac 1536K, 6144
build421 ,			,			,			yr4dwpn_pass1sc256ac 1920K, 7680
build421 ,			,			,			yr4dwpn_pass1sc256ac 2M, 13
build421 ,			,			,			yr4dwpn_pass1sc256ac 2560K, 10240
build421 ,			,			,			yr4dwpn_pass1sc256ac 3200K, 12800

build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc384ac 1728K, 4608
build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc384ac 2304K, 6144
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc384ac 2880K, 7680
build421 ,			,			,			yr4dwpn_pass1sc384ac 3M, 13
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc384ac 3840K, 10240
build421 ,			,			,			yr4dwpn_pass1sc384ac 4800K, 12800

build421 ,			,			,			yr4dwpn_pass1sc512ac 2304K, 4608
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc512ac 3M, 6144
build421 ,			CORE,			,			yr4dwpn_pass1sc512ac 3840K, 7680
build421 CORE_64,		CORE_32,		,			yr4dwpn_pass1sc512ac 4M, 13
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc512ac 5M, 10240
build421 ,			CORE,			,			yr4dwpn_pass1sc512ac 6400K, 12800

build421 ,			CORE_32,		,			yr4dwpn_pass1sc640ac 2880K, 4608
build421 ,			,			,			yr4dwpn_pass1sc640ac 3840K, 6144
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc640ac 4800K, 7680
build421 ,			CORE_32,		,			yr4dwpn_pass1sc640ac 5M, 13
build421 FMA3_64,		,			,			yr4dwpn_pass1sc640ac 6400K, 10240
build421 CORE_64,		,			,			yr4dwpn_pass1sc640ac 8000K, 12800

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc768ac 3456K, 4608
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc768ac 4608K, 6144
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc768ac 5760K, 7680
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768ac 6M, 13
build421 ,			,			,			yr4dwpn_pass1sc768ac 7680K, 10240
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768ac 9600K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1024ac 4608K, 4608
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1024ac 6M, 6144
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1024ac 7680K, 7680
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc1024ac 8M, 13
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1024ac 10M, 10240
build421 ,			,			,			yr4dwpn_pass1sc1024ac 12800K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1280ac 5760K, 4608
build421 ,			,			,			yr4dwpn_pass1sc1280ac 7680K, 6144
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1280ac 9600K, 7680
build421 ,			,			,			yr4dwpn_pass1sc1280ac 10M, 13
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1280ac 12800K, 10240
build421 ,			,			,			yr4dwpn_pass1sc1280ac 16000K, 12800

build421 ,			,			,			yr4dwpn_pass1sc1536ac 6912K, 4608
build421 ,			,			,			yr4dwpn_pass1sc1536ac 9M, 6144
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1536ac 11520K, 7680
build421 CORE_32,		FMA3_64,		,			yr4dwpn_pass1sc1536ac 12M, 13
build421 CORE_32 + FMA3_64,	,			,			yr4dwpn_pass1sc1536ac 15M, 10240
build421 FMA3_64,		,			,			yr4dwpn_pass1sc1536ac 19200K, 12800

build421 ,			,			,			yr4dwpn_pass1sc2048ac 9M, 4608
build421 ,			,			,			yr4dwpn_pass1sc2048ac 12M, 6144
build421 ,			,			,			yr4dwpn_pass1sc2048ac 15M, 7680
build421 CORE_32,		FMA3_64,		,			yr4dwpn_pass1sc2048ac 16M, 13
build421 CORE,			,			,			yr4dwpn_pass1sc2048ac 20M, 10240
build421 ,			,			,			yr4dwpn_pass1sc2048ac 25M, 12800

;build421 ,			,			,			yr4dwpn_pass1sc2560ac 11520K, 4608
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 15M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 19200K, 7680
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 20M, 13
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 25M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 32000K, 12800

;build421 ,			,			,			yr4dwpn_pass1sc3072ac 13824K, 4608
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 18M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 23040K, 7680
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 24M, 13
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 30M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 38400K, 12800

;build421 ,			,			,			yr4dwpn_pass1sc4096ac 18M, 4608
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 24M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 30M, 7680
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 32M, 13
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 40M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 50M, 12800

;build421 ,			,			,			yr4dwpn_pass1sc5120ac 23040K, 4608
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 30M, 6144
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 38400K, 7680
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 40M, 13
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 50M, 10240
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 64000K, 12800
ENDIF

_TEXT	ENDS
END
