; Copyright 2011-2017 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These FFTs were split out of yr4dwpn.asm because of MASM limitations.
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

buildfor CORE    + FMA3_64,	ypass2gen 2304
buildfor CORE    + FMA3_64,	ypass2gen 3072
buildfor CORE    + FMA3_64,	ypass2gen 3840
buildfor CORE    + FMA3_64,	ypass2gen 12
buildfor CORE    + FMA3_64,	ypass2gen 5120
buildfor CORE    + FMA3_64,	ypass2gen 6400

;; Routines for many FFT sizes

	; The 12 levels variants (2304, 3072, 3840, 4096, 5120, 6400)

IFNDEF X86_64
build421 ,			,			,			yr4dwpn_pass1sc128 288K, 2304
build421 ,			,			,			yr4dwpn_pass1sc128 384K, 3072
build421 ,			,			,			yr4dwpn_pass1sc128 480K, 3840
build421 ,			CORE_32,		,			yr4dwpn_pass1sc128 512K, 12
build421 ,			,			,			yr4dwpn_pass1sc128 640K, 5120
build421 ,			,			,			yr4dwpn_pass1sc128 800K, 6400

build421 ,			,			,			yr4dwpn_pass1sc256 576K, 2304
build421 ,			,			,			yr4dwpn_pass1sc256 768K, 3072
build421 ,			,			,			yr4dwpn_pass1sc256 960K, 3840
build421 CORE_64,		FMA3_64,		,			yr4dwpn_pass1sc256 1M, 12
build421 ,			,			,			yr4dwpn_pass1sc256 1280K, 5120
build421 ,			,			,			yr4dwpn_pass1sc256 1600K, 6400

build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc320 720K, 2304
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc320 960K, 3072
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc320 1200K, 3840
build421 CORE_64,		CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc320 1280K, 12
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc320 1600K, 5120
build421 ,			,			,			yr4dwpn_pass1sc320 2000K, 6400

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc384 864K, 2304
build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc384 1152K, 3072
build421 ,			CORE_32,		,			yr4dwpn_pass1sc384 1440K, 3840
build421 CORE_64,		CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc384 1536K, 12
build421 ,			CORE_32,		,			yr4dwpn_pass1sc384 1920K, 5120
build421 CORE_64 + FMA3_64,	CORE + FMA3_64,		CORE_64 + FMA3_64,	yr4dwpn_pass1sc384 2400K, 6400

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc448 1008K, 2304
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc448 1344K, 3072
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc448 1680K, 3840
build421 CORE_64,		CORE_32,		FMA3_64,		yr4dwpn_pass1sc448 1792K, 12
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc448 2240K, 5120
build421 ,			CORE_32,		,			yr4dwpn_pass1sc448 2800K, 6400

build421 ,			,			,			yr4dwpn_pass1sc512 1152K, 2304
build421 ,			,			,			yr4dwpn_pass1sc512 1536K, 3072
build421 ,			,			,			yr4dwpn_pass1sc512 1920K, 3840
build421 CORE_64,		CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc512 2M, 12
build421 ,			,			,			yr4dwpn_pass1sc512 2560K, 5120
build421 ,			,			,			yr4dwpn_pass1sc512 3200K, 6400

build421 ,			CORE_64,		,			yr4dwpn_pass1sc640 1440K, 2304
build421 ,			,			,			yr4dwpn_pass1sc640 1920K, 3072
build421 CORE_64 + FMA3_64,	CORE_64 + FMA3_64,	CORE_64 + FMA3_64,	yr4dwpn_pass1sc640 2400K, 3840
build421 ,			CORE_32,		CORE_64 + FMA3_64,	yr4dwpn_pass1sc640 2560K, 12
build421 ,			CORE_32,		,			yr4dwpn_pass1sc640 3200K, 5120
build421 ,			,			,			yr4dwpn_pass1sc640 4000K, 6400

build421 ,			CORE_64,		,			yr4dwpn_pass1sc768 1728K, 2304
build421 ,			,			,			yr4dwpn_pass1sc768 2304K, 3072
build421 ,			,			,			yr4dwpn_pass1sc768 2880K, 3840
build421 ,			CORE_32,		CORE_64 + FMA3_64,	yr4dwpn_pass1sc768 3M, 12
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768 3840K, 5120
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768 4800K, 6400

build421 ,			,			,			yr4dwpn_pass1sc896 2016K, 2304
build421 ,			,			,			yr4dwpn_pass1sc896 2688K, 3072
build421 ,			,			,			yr4dwpn_pass1sc896 3360K, 3840
build421 ,			,			CORE_64 + FMA3_64,	yr4dwpn_pass1sc896 3584K, 12
build421 ,			CORE_32,		,			yr4dwpn_pass1sc896 4480K, 5120
build421 ,			CORE_32,		,			yr4dwpn_pass1sc896 5600K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1024 2304K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1024 3M, 3072
build421 ,			,			,			yr4dwpn_pass1sc1024 3840K, 3840
build421 ,			,			FMA3_64,		yr4dwpn_pass1sc1024 4M, 12
build421 ,			,			,			yr4dwpn_pass1sc1024 5M, 5120
build421 ,			,			,			yr4dwpn_pass1sc1024 6400K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1280 2880K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1280 3840K, 3072
build421 ,			,			,			yr4dwpn_pass1sc1280 4800K, 3840
build421 ,			,			FMA3_64,		yr4dwpn_pass1sc1280 5M, 12
build421 ,			,			,			yr4dwpn_pass1sc1280 6400K, 5120
build421 ,			,			,			yr4dwpn_pass1sc1280 8000K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1536 3456K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1536 4608K, 3072
build421 ,			,			,			yr4dwpn_pass1sc1536 5760K, 3840
build421 ,			,			,			yr4dwpn_pass1sc1536 6M, 12
build421 ,			,			,			yr4dwpn_pass1sc1536 7680K, 5120
build421 ,			,			,			yr4dwpn_pass1sc1536 9600K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1792 4032K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1792 5376K, 3072
build421 ,			,			,			yr4dwpn_pass1sc1792 6720K, 3840
build421 ,			,			,			yr4dwpn_pass1sc1792 7M, 12
build421 ,			,			,			yr4dwpn_pass1sc1792 8960K, 5120
build421 ,			,			,			yr4dwpn_pass1sc1792 11200K, 6400

build421 ,			,			,			yr4dwpn_pass1sc2048 4608K, 2304
build421 ,			,			,			yr4dwpn_pass1sc2048 6M, 3072
build421 ,			,			,			yr4dwpn_pass1sc2048 7680K, 3840
build421 ,			,			,			yr4dwpn_pass1sc2048 8M, 12
build421 ,			,			,			yr4dwpn_pass1sc2048 10M, 5120
build421 ,			,			,			yr4dwpn_pass1sc2048 12800K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc2560 5760K, 2304
;build421 ,			,			,			yr4dwpn_pass1sc2560 7680K, 3072
;build421 ,			,			,			yr4dwpn_pass1sc2560 9600K, 3840
;build421 ,			,			,			yr4dwpn_pass1sc2560 10M, 12
;build421 ,			,			,			yr4dwpn_pass1sc2560 12800K, 5120
;build421 ,			,			,			yr4dwpn_pass1sc2560 16000K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc3072 6912K, 2304
;build421 ,			,			,			yr4dwpn_pass1sc3072 9M, 3072
;build421 ,			,			,			yr4dwpn_pass1sc3072 11520K, 3840
;build421 ,			,			,			yr4dwpn_pass1sc3072 12M, 12
;build421 ,			,			,			yr4dwpn_pass1sc3072 15M, 5120
;build421 ,			,			,			yr4dwpn_pass1sc3072 19200K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc3584 8064K, 2304
;build421 ,			,			,			yr4dwpn_pass1sc3584 10752K, 3072
;build421 ,			,			,			yr4dwpn_pass1sc3584 13440K, 3840
;build421 ,			,			,			yr4dwpn_pass1sc3584 14M, 12
;build421 ,			,			,			yr4dwpn_pass1sc3584 17920K, 5120
;build421 ,			,			,			yr4dwpn_pass1sc3584 22400K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc4096 9M, 2304
;build421 ,			,			,			yr4dwpn_pass1sc4096 12M, 3072
;build421 ,			,			,			yr4dwpn_pass1sc4096 15M, 3840
;build421 ,			,			,			yr4dwpn_pass1sc4096 16M, 12
;build421 ,			,			,			yr4dwpn_pass1sc4096 20M, 5120
;build421 ,			,			,			yr4dwpn_pass1sc4096 25M, 6400
ENDIF

	; The all-complex 12 levels variants (2304, 3072, 3840, 4096, 5120, 6400)

IFNDEF X86_64
build421 ,			,			,			yr4dwpn_pass1sc128ac 288K, 2304
build421 ,			,			,			yr4dwpn_pass1sc128ac 384K, 3072
build421 ,			,			,			yr4dwpn_pass1sc128ac 480K, 3840
build421 CORE_64,		,			,			yr4dwpn_pass1sc128ac 512K, 12
build421 ,			,			,			yr4dwpn_pass1sc128ac 640K, 5120
build421 ,			,			,			yr4dwpn_pass1sc128ac 800K, 6400

build421 ,			,			,			yr4dwpn_pass1sc256ac 576K, 2304
build421 ,			,			,			yr4dwpn_pass1sc256ac 768K, 3072
build421 ,			,			,			yr4dwpn_pass1sc256ac 960K, 3840
build421 CORE_64,		FMA3_64,		,			yr4dwpn_pass1sc256ac 1M, 12
build421 ,			,			,			yr4dwpn_pass1sc256ac 1280K, 5120
build421 ,			,			,			yr4dwpn_pass1sc256ac 1600K, 6400

build421 ,			,			CORE    + FMA3_64,	yr4dwpn_pass1sc384ac 864K, 2304
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc384ac 1152K, 3072
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc384ac 1440K, 3840
build421 CORE_64,		CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc384ac 1536K, 12
build421 ,			CORE_32,		,			yr4dwpn_pass1sc384ac 1920K, 5120
build421 ,			CORE,			,			yr4dwpn_pass1sc384ac 2400K, 6400

build421 ,			CORE_64,		,			yr4dwpn_pass1sc512ac 1152K, 2304
build421 ,			,			,			yr4dwpn_pass1sc512ac 1536K, 3072
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc512ac 1920K, 3840
build421 CORE_64,		CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc512ac 2M, 12
build421 ,			CORE_32,		,			yr4dwpn_pass1sc512ac 2560K, 5120
build421 ,			CORE,			,			yr4dwpn_pass1sc512ac 3200K, 6400

build421 ,			CORE_32,		,			yr4dwpn_pass1sc640ac 1440K, 2304
build421 ,			,			,			yr4dwpn_pass1sc640ac 1920K, 3072
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc640ac 2400K, 3840
build421 CORE_64,		FMA3_64,		,			yr4dwpn_pass1sc640ac 2560K, 12
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc640ac 3200K, 5120
build421 ,			,			,			yr4dwpn_pass1sc640ac 4000K, 6400

build421 ,			,			,			yr4dwpn_pass1sc768ac 1728K, 2304
build421 ,			,			,			yr4dwpn_pass1sc768ac 2304K, 3072
build421 ,			,			,			yr4dwpn_pass1sc768ac 2880K, 3840
build421 ,			,			,			yr4dwpn_pass1sc768ac 3M, 12
build421 ,			,			,			yr4dwpn_pass1sc768ac 3840K, 5120
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768ac 4800K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1024ac 2304K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1024ac 3M, 3072
build421 ,			,			,			yr4dwpn_pass1sc1024ac 3840K, 3840
build421 ,			,			FMA3_64,		yr4dwpn_pass1sc1024ac 4M, 12
build421 ,			,			,			yr4dwpn_pass1sc1024ac 5M, 5120
build421 ,			,			,			yr4dwpn_pass1sc1024ac 6400K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1280ac 2880K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1280ac 3840K, 3072
build421 ,			,			,			yr4dwpn_pass1sc1280ac 4800K, 3840
build421 ,			,			,			yr4dwpn_pass1sc1280ac 5M, 12
build421 ,			,			,			yr4dwpn_pass1sc1280ac 6400K, 5120
build421 ,			,			,			yr4dwpn_pass1sc1280ac 8000K, 6400

build421 ,			,			,			yr4dwpn_pass1sc1536ac 3456K, 2304
build421 ,			,			,			yr4dwpn_pass1sc1536ac 4608K, 3072
build421 ,			,			,			yr4dwpn_pass1sc1536ac 5760K, 3840
build421 ,			,			,			yr4dwpn_pass1sc1536ac 6M, 12
build421 ,			,			,			yr4dwpn_pass1sc1536ac 7680K, 5120
build421 ,			,			,			yr4dwpn_pass1sc1536ac 9600K, 6400

build421 ,			,			,			yr4dwpn_pass1sc2048ac 4608K, 2304
build421 ,			,			,			yr4dwpn_pass1sc2048ac 6M, 3072
build421 ,			,			,			yr4dwpn_pass1sc2048ac 7680K, 3840
build421 ,			,			,			yr4dwpn_pass1sc2048ac 8M, 12
build421 CORE_32,		,			,			yr4dwpn_pass1sc2048ac 10M, 5120
build421 ,			,			,			yr4dwpn_pass1sc2048ac 12800K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc2560ac 5760K, 2304
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 7680K, 3072
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 9600K, 3840
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 10M, 12
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 12800K, 5120
;build421 ,			,			,			yr4dwpn_pass1sc2560ac 16000K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc3072ac 6912K, 2304
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 9M, 3072
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 11520K, 3840
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 12M, 12
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 15M, 5120
;build421 ,			,			,			yr4dwpn_pass1sc3072ac 19200K, 6400

;build421 ,			,			,			yr4dwpn_pass1sc4096ac 9M, 2304
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 12M, 3072
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 15M, 3840
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 16M, 12
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 20M, 5120
;build421 ,			,			,			yr4dwpn_pass1sc4096ac 25M, 6400

;build421 ,			,			,			yr4dwpn_pass1sc5120ac 11520K, 2304
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 15M, 3072
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 19200K, 3840
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 20M, 12
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 25M, 5120
;build421 ,			,			,			yr4dwpn_pass1sc5120ac 32000K, 6400
ENDIF

_TEXT	ENDS
END
