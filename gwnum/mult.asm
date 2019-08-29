; Copyright 1995-2019 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the setup, common routines, and global variables
; for the various discrete-weighted transforms
;

	TITLE   gwdata

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac

VERSION_NUMBER = 2908		;; Version 29.8

;
; Global variables needed by FFT setup code
;

IFNDEF X86_64
_GWDATA SEGMENT PAGE PUBLIC 'DATA'
ELSE
_GWDATA SEGMENT PAGE
ENDIF

; Jmptable definitions

; The flags word of a procedure entry is divided as follows:
;	80000000h		always on
;	2 SHL 26		(no prefetching - not used by gwnum)
;	1 SHL 26		(in_place)
;	fft_type SHL 21		(hg=0, r4=1, r4delay=2, r4dwpn=3)
;	arch SHL 17		(blend=0,p3=1,p4=2,core=3,fma3=4,k8=5,etc.)
;	best_impl_for SHL 13	(CORE2=0,I7=1,P4_1024=2,etc.)
;	clm SHL 9		(1,2,4,8)
;	pass2size_over_64	(many valid values)

; Machine types.  FFTs often have several different implementations.  For the most common CPUs,
; we run benchmarks to determine which FFT implementation is the fastest.  A jmptable entry is
; built for the fastest implementation.  At runtime, the GWNUM setup code figures out what CPU
; it is running on and selects the best matching FFT implementation.

; The jmptable uses a bit number to encode the "best implementation for machine type" value.

	BIF_CORE2   = 0			; Core 2 CPUs, big L2 caches
	BIF_CORE2_512 = 1		; Core 2 Celerons, 512K L2 cache, no L3 cache
	BIF_I7	    = 2			; Core i3/5/7 CPUs, 256K L2, big L3 caches
	BIF_FMA3    = 3			; Core i3/5/7 CPUs, Haswell architecture with FMA3 support
	BIF_P4_1024 = 4			; Pentium 4, 1MB cache
	BIF_P4TP_512  = 5		; Pentium 4, 512K cache (did not support EMT64)
	BIF_P4TP_256  = 6		; Pentium 4, 256K cache (did not support EMT64)
	BIF_SKX	    = 7			; Intel CPUs supporting AVX-512 (Skylake-X)
	BIF_K8	    = 8			; AMD K8 CPUs
	BIF_K10	    = 9			; AMD K10 CPUs
	BIF_RYZEN   = 11		; AMD Ryzen CPUs

; The PRCENTRY macro uses a bit mask to encode the "best implementation for machine type" values.
; The PRCENTRY macro also lets us specify different "best implementation for machine type" values
; for 32-bit and 64-bit operating systems.

	CORE2	= 00000001h		; Core 2 CPUs, big L2 caches
	CORE2_512 = 00000002h		; Core 2 Celerons, 512K L2 cache, no L2 cache
	I7	= 00000004h		; Core i3/5/7 CPUs, 256K L2, big L3 caches
	FMA3	= 00000008h		; Core i3/5/7 CPUs, Haswell architecture with FMA3 support
	P4_1024	= 00000010h		; Pentium 4, 1MB cache
	P4TP_512 = 00000020h		; Pentium 4, 512K cache (did not support EMT64)
	P4TP_256 = 00000040h		; Pentium 4, 256K cache (did not support EMT64)
	SKX	= 00000080h		; Intel CPUs supporting AVX-512 (Skylake-X)
	K8	= 00000100h		; AMD K8 CPUs
	K10	= 00000200h		; AMD K10 CPUs
	RYZEN	= 00000800h		; AMD Ryzen CPUs

	MAYBE_P4 = 01000000h		; This _P4 implementation may be best on some CPUs
	MAYBE_CORE = 02000000h		; This _CORE implementation may be best on some CPUs
	MAYBE_K8 = 08000000h		; This _K8 implementation may be best on some CPUs
	MAYBE_K10 = 10000000h		; This _K10 implementation may be best on some CPUs

; Apple products do not contain AMD CPUs, nor small-cache Pentium 4s.  They don't contain
; large cache Pentium 4s either, but gwnum treats Core CPUs as large-cache Pentium 4s.
; Since the Core CPUs Apple used didn't support EMT64, we don't need to support large
; cache Pentium 4's in 64-bit mode.  Apple doesn't use small cache Core 2 Celerons.

IFDEF APPLE
	K8 = 0
	K10 = 0
	RYZEN = 0
	P4TP_512 = 0
	P4TP_256 = 0
IFDEF X86_64
	P4_1024	= 0
ENDIF
	CORE2_512 = 0
ENDIF

; Define constants that let us specify different best FFT implementations for
; 32-bit and 64-bit CPUs.

IFDEF X86_64
	CORE2_32 = 0
	CORE2_64 = CORE2
	CORE2_512_32 = 0
	CORE2_512_64 = CORE2_512
	I7_32 = 0
	I7_64 = I7
	FMA3_32 = 0
	FMA3_64 = FMA3
	P4_1024_32 = 0
	P4_1024_64 = P4_1024
	P4TP_512 = 0
	P4TP_256 = 0
	K8_32 = 0
	K8_64 = K8
	K10_32 = 0
	K10_64 = K10
	RYZEN_32 = 0
	RYZEN_64 = RYZEN
	MAYBE_P4_32 = 0
	MAYBE_P4_64 = MAYBE_P4
	MAYBE_CORE_32 = 0
	MAYBE_CORE_64 = MAYBE_CORE
	MAYBE_K8_32 = 0
	MAYBE_K8_64 = MAYBE_K8
	MAYBE_K10_32 = 0
	MAYBE_K10_64 = MAYBE_K10
ELSE
	CORE2_32 = CORE2
	CORE2_64 = 0
	CORE2_512_32 = CORE2_512
	CORE2_512_64 = 0
	I7_32 = I7
	I7_64 = 0
	FMA3_32 = FMA3
	FMA3_64 = 0
	P4_1024_32 = P4_1024
	P4_1024_64 = 0
	K8_32 = K8
	K8_64 = 0
	K10_32 = K10
	K10_64 = 0
	RYZEN_32 = RYZEN
	RYZEN_64 = 0
	MAYBE_P4_32 = MAYBE_P4
	MAYBE_P4_64 = 0
	MAYBE_CORE_32 = MAYBE_CORE
	MAYBE_CORE_64 = 0
	MAYBE_K8_32 = MAYBE_K8
	MAYBE_K8_64 = 0
	MAYBE_K10_32 = MAYBE_K10
	MAYBE_K10_64 = 0
ENDIF

; Macros to build jmptable entries

PRCSTRT MACRO max_exp, fftlen, speed
	DD	max_exp, fftlen, speed
	saved_fftlen = fftlen			; Remember for use in PRCENTRY macro
	ENDM

X87PRC MACRO procname, mem_needed, p2larg, clmarg
	EXTRN	procname:PROC
	IF (@INSTR (,procname,PPRO)) NE 0
		arch = 0
	ELSEIF (@INSTR (,procname,P3)) NE 0
		arch = 1
	ELSE
		bad_x87_arch
	ENDIF
	IFB <clmarg>
		clm = 0
	ELSE
		clm = clmarg
	ENDIF
	IFB <p2larg>
		pass2_size_over_64 = 0
	ELSE
		pass2_size_over_64 = (1 SHL p2larg) SHR 6
	ENDIF
	flags = 0
	;; Kludge - if mem_needed is below 400000 then pass 2 is in place
	IF (mem_needed LT 400000 AND pass2_size_over_64 NE 0)
		flags = 1
	ENDIF
	DD	80000000h + flags SHL 26 + arch SHL 17 + clm SHL 9 + pass2_size_over_64
	DP	OFFSET procname
	DD	mem_needed
	ENDM

; Used for small one-pass FFTs.  At present there is only one FFT implementation.

PRCENTRY MACRO procname, mem_needed, best_impl_for
	IFB <best_impl_for>
		PRCENTRY3 procname, mem_needed, 0
	ELSE
		PRCENTRY3 procname, mem_needed, best_impl_for
	ENDIF
	ENDM

; Used for FFTs with multiple implementations.  This macro generates a table entry
; for each CPU architecture where this is the fastest FFT implementation.

PRCENTRY2 MACRO procname, mem_needed, best_impl_for
	;; Copy input argument
	IFB <best_impl_for>
		bif = 0
	ELSE
		bif = best_impl_for
	ENDIF
	;; Hack to force output of jmptable entries.  Used when we are benchmarking
	;; all possible FFT implementations to determine best implementation
	IFDEF IMPL_ALL_CORE
		bif = bif OR MAYBE_CORE
	ENDIF
	IF (@INSTR (,&procname,_op)) EQ 0
	IFDEF IMPL_ALL_P4
		bif = bif OR P4_1024
	ENDIF
	IFDEF IMPL_ALL_P4TP
		;; We do not generate a TLB priming version when we are not prefetching
		IF (@INSTR (,procname,_np_)) EQ 0
			bif = bif OR P4TP_512
		ENDIF
	ENDIF
	IFDEF IMPL_ALL_K8
		bif = bif OR MAYBE_K8
	ENDIF
	IFDEF IMPL_ALL_K10
		bif = bif OR MAYBE_K10
	ENDIF
	ENDIF
	;; We arbitrarily decide that older architectures are ill-equiped to run large FFTs.
	;; Thus, we make sure that bif bit is cleared
	IF (saved_fftlen GT 4194304) AND ((bif AND P4TP_256) NE 0)
		bif = bif - P4TP_256
	ENDIF
	IF (saved_fftlen GT 4194304) AND ((bif AND P4TP_512) NE 0)
		bif = bif - P4TP_512
	ENDIF
	IF (saved_fftlen GT 6291456) AND ((bif AND P4_1024) NE 0)
		bif = bif - P4_1024
	ENDIF
	IF (saved_fftlen GT 6291456) AND ((bif AND K8) NE 0)
		bif = bif - K8
	ENDIF
	;; Output one jmptable entry for each CPU where this is the best FFT implementation.
	IF (bif AND CORE2) NE 0
		PRCENTRY3 @CATSTR(procname, <_CORE>), mem_needed, CORE2
		bif = bif AND NOT MAYBE_CORE
	ENDIF
	IF (bif AND CORE2_512) NE 0
		PRCENTRY3 @CATSTR(procname, <_CORE>), mem_needed, CORE2_512
		bif = bif AND NOT MAYBE_CORE
	ENDIF
	IF (bif AND I7) NE 0
		PRCENTRY3 @CATSTR(procname, <_CORE>), mem_needed, I7
		bif = bif AND NOT MAYBE_CORE
	ENDIF
	IF (bif AND FMA3) NE 0
		PRCENTRY3 @CATSTR(procname, <_FMA3>), mem_needed, FMA3
	ENDIF
	IF (bif AND RYZEN) NE 0
		PRCENTRY3 @CATSTR(procname, <_FMA3>), mem_needed, RYZEN
	ENDIF
	IF (bif AND SKX) NE 0
		PRCENTRY3 @CATSTR(procname, <_SKX>), mem_needed, SKX
	ENDIF
	IF (bif AND P4_1024) NE 0
		PRCENTRY3 @CATSTR(procname, <_P4>), mem_needed, P4_1024
		bif = bif AND NOT MAYBE_P4
	ENDIF
	IF (bif AND P4TP_512) NE 0
		IF (@INSTR (,procname,_np_)) EQ 0
			PRCENTRY3 @CATSTR(procname, <_P4TP>), mem_needed, P4TP_512
		ELSE
			PRCENTRY3 @CATSTR(procname, <_P4>), mem_needed, P4TP_512
		ENDIF
	ENDIF
	IF (bif AND P4TP_256) NE 0
		IF (@INSTR (,procname,_np_)) EQ 0
			PRCENTRY3 @CATSTR(procname, <_P4TP>), mem_needed, P4TP_256
		ELSE
			PRCENTRY3 @CATSTR(procname, <_P4>), mem_needed, P4TP_256
		ENDIF
	ENDIF
	IF (bif AND K8) NE 0
		PRCENTRY3 @CATSTR(procname, <_K8>), mem_needed, K8
		bif = bif AND NOT MAYBE_K8
	ENDIF
	IF (bif AND K10) NE 0
		PRCENTRY3 @CATSTR(procname, <_K10>), mem_needed, K10
		bif = bif AND NOT MAYBE_K10
	ENDIF
	IF (bif AND MAYBE_P4) NE 0
		PRCENTRY3 @CATSTR(procname, <_P4>), mem_needed, 0
	ENDIF
	IF (bif AND MAYBE_CORE) NE 0
		PRCENTRY3 @CATSTR(procname, <_CORE>), mem_needed, 0
	ENDIF
	IF (bif AND MAYBE_K8) NE 0
		PRCENTRY3 @CATSTR(procname, <_K8>), mem_needed, 0
	ENDIF
	IF (bif AND MAYBE_K10) NE 0
		PRCENTRY3 @CATSTR(procname, <_K10>), mem_needed, 0
	ENDIF
	ENDM

;; Macro to create an entry for one CPU that is optimized for a different CPU type.
;; For example, on a P4, some CORE optimized implementations are faster than the
;; P4 optimized implementations.
PRCENTRY2A MACRO procname, mem_needed, best_impl_for
	;; When building a special version for timing all implementations,
	;; don't bother creating these specialized entries.
	IFDEF IMPL_ALL_CORE
		exitm
	ENDIF
	IFDEF IMPL_ALL_P4
		exitm
	ENDIF
	IFDEF IMPL_ALL_P4TP
		exitm
	ENDIF
	IFDEF IMPL_ALL_K8
		exitm
	ENDIF
	IFDEF IMPL_ALL_K10
		exitm
	ENDIF
	;; Copy input argument
	bif = best_impl_for
	;; Output one jmptable entry for each CPU where this is the best FFT implementation.
	IF (bif AND CORE2) NE 0
		PRCENTRY3 procname, mem_needed, CORE2
	ENDIF
	IF (bif AND CORE2_512) NE 0
		PRCENTRY3 procname, mem_needed, CORE2_512
	ENDIF
	IF (bif AND I7) NE 0
		PRCENTRY3 procname, mem_needed, I7
	ENDIF
	IF (bif AND FMA3) NE 0
		PRCENTRY3 procname, mem_needed, FMA3
	ENDIF
	IF (bif AND SKX) NE 0
		PRCENTRY3 procname, mem_needed, SKX
	ENDIF
	IF (bif AND P4_1024) NE 0
		PRCENTRY3 procname, mem_needed, P4_1024
	ENDIF
	IF (bif AND P4TP_512) NE 0
		PRCENTRY3 procname, mem_needed, P4TP_512
	ENDIF
	IF (bif AND P4TP_256) NE 0
		PRCENTRY3 procname, mem_needed, P4TP_256
	ENDIF
	IF (bif AND K8) NE 0
		PRCENTRY3 procname, mem_needed, K8
	ENDIF
	IF (bif AND K10) NE 0
		PRCENTRY3 procname, mem_needed, K10
	ENDIF
	ENDM

PRCENTRY3 MACRO procname, mem_needed, bif
	EXTRN	&procname:PROC
	flags = 0
	IF (@INSTR (,&procname,_np_)) NE 0
		flags = flags + 2
	ENDIF
	IF (@INSTR (,&procname,_ip_)) NE 0
		flags = flags + 1
	ENDIF
	IF (@INSTR (,&procname,_hg_)) NE 0
		fft_type = 0
	ELSEIF (@INSTR (,&procname,_r4_)) NE 0
		fft_type = 1
	ELSEIF (@INSTR (,&procname,_r4delay_)) NE 0
		fft_type = 2
	ELSEIF (@INSTR (,&procname,_r4dwpn_)) NE 0
		fft_type = 3
	ELSE
		bad_sse2_fft_type
	ENDIF
	IF (@INSTR (,&procname,_P4TP)) NE 0
		arch = 1
	ELSEIF (@INSTR (,&procname,_P4)) NE 0
		arch = 2
	ELSEIF (@INSTR (,&procname,_CORE)) NE 0
		arch = 3
	ELSEIF (@INSTR (,&procname,_FMA3)) NE 0
		arch = 4
	ELSEIF (@INSTR (,&procname,_K8)) NE 0
		arch = 5
	ELSEIF (@INSTR (,&procname,_K10)) NE 0
		arch = 6
	ELSEIF (@INSTR (,&procname,_SKX)) NE 0
		arch = 8
	ELSEIF (@INSTR (,&procname,_BLEND)) NE 0
		arch = 0
	ELSE
		bad_sse2_arch
	ENDIF
	IF (@INSTR (,&procname,_1_)) NE 0
		clm = 1
	ELSEIF (@INSTR (,&procname,_2_)) NE 0
		clm = 2
	ELSEIF (@INSTR (,&procname,_4_)) NE 0
		clm = 4
	ELSEIF (@INSTR (,&procname,_8_)) NE 0
		clm = 8
	ELSE
		clm = 0
	ENDIF
	IF (@INSTR (,&procname,_op_)) NE 0
		pass2_size_over_64 = 0
	ELSEIF (@INSTR (,&procname,_48_)) NE 0
		pass2_size_over_64 = 511
	ELSEIF (@INSTR (,&procname,_6_)) NE 0
		pass2_size_over_64 = 64 SHR 6
	ELSEIF (@INSTR (,&procname,_80_)) NE 0
		pass2_size_over_64 = 510
	ELSEIF (@INSTR (,&procname,_192_)) NE 0
		pass2_size_over_64 = 192 SHR 6
	ELSEIF (@INSTR (,&procname,_320_)) NE 0
		pass2_size_over_64 = 320 SHR 6
	ELSEIF (@INSTR (,&procname,_768_)) NE 0
		pass2_size_over_64 = 768 SHR 6
	ELSEIF (@INSTR (,&procname,_10_)) NE 0
		pass2_size_over_64 = 1024 SHR 6
	ELSEIF (@INSTR (,&procname,_1280_)) NE 0
		pass2_size_over_64 = 1280 SHR 6
	ELSEIF (@INSTR (,&procname,_1536_)) NE 0
		pass2_size_over_64 = 1536 SHR 6
	ELSEIF (@INSTR (,&procname,_11_)) NE 0
		pass2_size_over_64 = 2048 SHR 6
	ELSEIF (@INSTR (,&procname,_2304_)) NE 0
		pass2_size_over_64 = 2304 SHR 6
	ELSEIF (@INSTR (,&procname,_2560_)) NE 0
		pass2_size_over_64 = 2560 SHR 6
	ELSEIF (@INSTR (,&procname,_3072_)) NE 0
		pass2_size_over_64 = 3072 SHR 6
	ELSEIF (@INSTR (,&procname,_3840_)) NE 0
		pass2_size_over_64 = 3840 SHR 6
	ELSEIF (@INSTR (,&procname,_12_)) NE 0
		pass2_size_over_64 = 4096 SHR 6
	ELSEIF (@INSTR (,&procname,_4608_)) NE 0
		pass2_size_over_64 = 4608 SHR 6
	ELSEIF (@INSTR (,&procname,_5120_)) NE 0
		pass2_size_over_64 = 5120 SHR 6
	ELSEIF (@INSTR (,&procname,_6144_)) NE 0
		pass2_size_over_64 = 6144 SHR 6
	ELSEIF (@INSTR (,&procname,_6400_)) NE 0
		pass2_size_over_64 = 6400 SHR 6
	ELSEIF (@INSTR (,&procname,_7680_)) NE 0
		pass2_size_over_64 = 7680 SHR 6
	ELSEIF (@INSTR (,&procname,_13_)) NE 0
		pass2_size_over_64 = 8192 SHR 6
	ELSEIF (@INSTR (,&procname,_9216_)) NE 0
		pass2_size_over_64 = 9216 SHR 6
	ELSEIF (@INSTR (,&procname,_10240_)) NE 0
		pass2_size_over_64 = 10240 SHR 6
	ELSEIF (@INSTR (,&procname,_12288_)) NE 0
		pass2_size_over_64 = 12288 SHR 6
	ELSEIF (@INSTR (,&procname,_12800_)) NE 0
		pass2_size_over_64 = 12800 SHR 6
	ELSEIF (@INSTR (,&procname,_15360_)) NE 0
		pass2_size_over_64 = 15360 SHR 6
	ELSEIF (@INSTR (,&procname,_14_)) NE 0
		pass2_size_over_64 = 16384 SHR 6
	ELSEIF (@INSTR (,&procname,_20480_)) NE 0
		pass2_size_over_64 = 20480 SHR 6
	ELSEIF (@INSTR (,&procname,_25600_)) NE 0
		pass2_size_over_64 = 25600 SHR 6
	ELSEIF (@INSTR (,&procname,_15_)) NE 0
		pass2_size_over_64 = 509
	ELSEIF (@INSTR (,&procname,_8_)) NE 0
		pass2_size_over_64 = 256 SHR 6
	ELSE
		unknown_pass2_size
	ENDIF
	IF bif EQ CORE2
		best_impl_for = BIF_CORE2
	ELSEIF bif EQ CORE2_512
		best_impl_for = BIF_CORE2_512
	ELSEIF bif EQ I7
		best_impl_for = BIF_I7
	ELSEIF bif EQ FMA3
		best_impl_for = BIF_FMA3
	ELSEIF bif EQ SKX
		best_impl_for = BIF_SKX
	ELSEIF bif EQ P4_1024
		best_impl_for = BIF_P4_1024
	ELSEIF bif EQ P4TP_512
		best_impl_for = BIF_P4TP_512
	ELSEIF bif EQ P4TP_256
		best_impl_for = BIF_P4TP_256
	ELSEIF bif EQ K8
		best_impl_for = BIF_K8
	ELSEIF bif EQ K10
		best_impl_for = BIF_K10
	ELSEIF bif EQ RYZEN
		best_impl_for = BIF_RYZEN
	ELSE
		best_impl_for = 0
	ENDIF
	DD	80000000h + flags SHL 26 + fft_type SHL 21 + arch SHL 17 + best_impl_for SHL 13 + clm SHL 9 + pass2_size_over_64
	DP	OFFSET procname
	DD	mem_needed
	ENDM

; Used for 64-bit AVX FFTs with shared code for pass 1 and for pass 2.  This macro generates a table entry
; for all possible FFT implemenations so that gwnum can use prior benchmarks to chose the fastest.
; In 32-bit mode, convert to an old-style PRCENTRY2 call

YPRCENTRY421 MACRO procprefix, fftlen, procsuffix, p2size, mem_needed, best_impl_for4, best_impl_for2, best_impl_for1
	YPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 4, mem_needed, best_impl_for4
	YPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 2, mem_needed, best_impl_for2
	YPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 1, mem_needed, best_impl_for1
	ENDM

YPRCENTRY2 MACRO procprefix, fftlen, procsuffix, p2size, clm, mem_needed, best_impl_for
	IFNDEF X86_64
		yproc32name CATSTR <&procprefix>,<_>,<&fftlen>,<&procsuffix>,<_>,%p2size,<_>,%clm
		PRCENTRY2 %yproc32name, mem_needed, best_impl_for
		exitm
	ENDIF
	;; Copy input argument
	IFB <best_impl_for>
		bif = 0
	ELSE
		bif = best_impl_for
	ENDIF
	;; Output one jmptable entry for each CPU where this is the best FFT implementation.
	IF (bif AND I7) NE 0
		YPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_CORE>, mem_needed, I7
	ENDIF
	IF (bif AND FMA3) NE 0
		YPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_FMA3>, mem_needed, FMA3
	ENDIF
	IF (bif AND RYZEN) NE 0
		YPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_FMA3>, mem_needed, RYZEN
	ENDIF
	;; Output a jmptable entry even if this implementation is not best
	IF (bif AND I7) EQ 0
		YPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_CORE>, mem_needed, 0
	ENDIF
	IF (bif AND (FMA3 OR RYZEN)) EQ 0
		YPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_FMA3>, mem_needed, 0
	ENDIF
	ENDM

YPRCENTRY3 MACRO procprefix, fftlen, procsuffix, p2size, clm, procarch, mem_needed, bif
	flags = 0
	IFNB <procsuffix>
	IF (@INSTR (,&procsuffix,_op)) NE 0
		one_pass_ffts_not_supported_here
	ENDIF
	IF (@INSTR (,&procsuffix,_np)) NE 0
		flags = flags + 2
	ENDIF
	IF (@INSTR (,&procsuffix,_ip)) NE 0
		flags = flags + 1
	ENDIF
	ENDIF
	IF (@INSTR (,&procprefix,_r4dwpn)) NE 0
		fft_type = 3
	ELSE
		bad_avx_fft_type
	ENDIF
	IF (@INSTR (,&procarch,_CORE)) NE 0
		arch = 3
	ELSEIF (@INSTR (,&procarch,_FMA3)) NE 0
		arch = 4
	ELSE
		bad_avx_arch
	ENDIF
	IF p2size LT 24
		real_p2size = 1 SHL p2size
	ELSE
		real_p2size = p2size
	ENDIF
	IF real_p2size EQ 48
		pass2_size_over_64 = 511
	ELSEIF real_p2size EQ 80
		pass2_size_over_64 = 510
	ELSE
		pass2_size_over_64 = real_p2size SHR 6
	ENDIF
	IF bif EQ I7
		best_impl_for = BIF_I7
	ELSEIF bif EQ FMA3
		best_impl_for = BIF_FMA3
	ELSEIF bif EQ RYZEN
		best_impl_for = BIF_RYZEN
	ELSE
		best_impl_for = 0
	ENDIF
	IF (@INSTR (,&fftlen,K)) NE 0
		p1size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 / real_p2size
	ELSE
		p1size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 * 1024 / real_p2size
	ENDIF
	yprocpass1name CATSTR <&procprefix>,<_>,%p1size,<&procsuffix>,<_>,%clm,<&procarch>
	IF p2size LT 24
		yprocpass2name CATSTR <ypass2_>,@SUBSTR(&procprefix,6),<_>,%p2size,<_levels>,<&procarch>
	ELSE
		yprocpass2name CATSTR <ypass2_>,@SUBSTR(&procprefix,6),<_>,%p2size,<&procarch>
	ENDIF
	EXTRN	yprocpass1name:PROC
	EXTRN	yprocpass2name:PROC
	DD	80000000h + 40000000h + flags SHL 26 + fft_type SHL 21 + arch SHL 17 + best_impl_for SHL 13 + clm SHL 9 + pass2_size_over_64
	DP	OFFSET yprocpass1name, OFFSET yprocpass2name
	DD	mem_needed
	ENDM

; Used for so called "one-pass" AVX-512 FFTs that use a wrapper pass 1 and a standard pass 2.

ZPRCENTRY1 MACRO procprefix, fftlen, procsuffix, mem_needed, best_impl_for1
	IF (@INSTR (,&fftlen,K)) NE 0
		p2size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 / 2
	ELSE
		p2size = fftlen / 2
	ENDIF
	p1size = 2
	ZPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 1, mem_needed, best_impl_for1
	ENDM

; Used for AVX-512 FFTs with shared code for pass 1 and for pass 2.  This macro generates a table entry
; for all possible FFT implemenations so that gwnum can use prior benchmarks to chose the fastest.

ZPRCENTRY421 MACRO procprefix, fftlen, procsuffix, p2size, mem_needed, best_impl_for4, best_impl_for2, best_impl_for1
	IF (@INSTR (,&fftlen,K)) NE 0
		p1size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 / p2size
	ELSE
		p1size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 * 1024 / p2size
	ENDIF
	;; Target a 1MB L2 cache size, these IFs must match zpass1gen's IFs in zmult.mac
	IF p1size LT 1024
	ZPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 4, mem_needed, best_impl_for4
	ENDIF
	IF p1size LT 2048
	ZPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 2, mem_needed, best_impl_for2
	ENDIF
	ZPRCENTRY2 procprefix, fftlen, procsuffix, p2size, 1, mem_needed, best_impl_for1
	ENDM

ZPRCENTRY2 MACRO procprefix, fftlen, procsuffix, p2size, clm, mem_needed, best_impl_for
	;; Copy input argument
	IFB <best_impl_for>
		bif = 0
	ELSE
		bif = best_impl_for
	ENDIF
	;; Output one jmptable entry for each CPU where this is the best FFT implementation.
	IF (bif AND SKX) NE 0
		ZPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_SKX>, mem_needed, SKX
	ENDIF
	;; Output a jmptable entry even if this implementation is not best
	IF (bif AND SKX) EQ 0
		ZPRCENTRY3 procprefix, fftlen, procsuffix, p2size, clm, <_SKX>, mem_needed, 0
	ENDIF
	ENDM

ZPRCENTRY3 MACRO procprefix, fftlen, procsuffix, p2size, clm, procarch, mem_needed, bif
	flags = 0
	IFNB <procsuffix>
	IF (@INSTR (,&procsuffix,_op)) NE 0
		one_pass_ffts_not_supported_here
	ENDIF
	IF (@INSTR (,&procsuffix,_np)) NE 0
		no_prefetch_ffts_not_supported_here
		flags = flags + 2
	ENDIF
	IF (@INSTR (,&procsuffix,_ip)) NE 0
		in_place_ffts_not_supported_here
		flags = flags + 1
	ENDIF
	ENDIF
	IF (@INSTR (,&procprefix,_r4dwpn)) NE 0
		fft_type = 3
	ELSE
		bad_avx512_fft_type
	ENDIF
	IF (@INSTR (,&procarch,_SKX)) NE 0
		arch = 8
	ELSE
		bad_avx512_arch
	ENDIF
	IF p2size EQ 32768
		pass2_size_over_64 = 509
	ELSE
		pass2_size_over_64 = p2size SHR 6
	ENDIF
	IF bif EQ SKX
		best_impl_for = BIF_SKX
	ELSEIF bif EQ RYZEN
		best_impl_for = BIF_RYZEN
	ELSE
		best_impl_for = 0
	ENDIF
	IF (@INSTR (,&fftlen,K)) NE 0
		p1size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 / p2size
	ELSEIF (@INSTR (,&fftlen,M)) NE 0
		p1size = @SUBSTR (fftlen, 1, @SIZESTR (fftlen) - 1) * 1024 * 1024 / p2size
	ELSE
		p1size = fftlen / p2size
	ENDIF
	zprocpass1name CATSTR <&procprefix>,<_>,%p1size,<&procsuffix>,<_>,%clm,<&procarch>
	zprocpass2name CATSTR <zpass2_>,@SUBSTR(&procprefix,6),<_>,%p2size,<&procarch>
	EXTRN	zprocpass1name:PROC
	EXTRN	zprocpass2name:PROC
	DD	80000000h + 40000000h + flags SHL 26 + fft_type SHL 21 + arch SHL 17 + best_impl_for SHL 13 + clm SHL 9 + pass2_size_over_64
	DP	OFFSET zprocpass1name, OFFSET zprocpass2name
	DD	mem_needed
	ENDM

;; jmptable entries for x87 FFTs

IFNDEF	X86_64
jmptable DD	755,	32,	0.0000036
 	X87PRC		fft32PPRO, 672
	DD			3, 1, 1, 1, 0
	DD	939,	40,	0.0000057
	X87PRC		fft40PPRO, 948
	DD			4, 1, 1, 1, 0
	DD	1113,	48,	0.0000065
	X87PRC		fft48PPRO, 1128
	DD			5, 1, 1, 1, 0
	DD	1303,	56,	0.0000084
	X87PRC		fft56PPRO, 1356
	DD			6, 1, 1, 1, 0
	DD	1499,	64,	0.0000083
	X87PRC		fft64PPRO, 1392
	DD			7, 1, 1, 1, 0
	DD	1857,	80,	0.0000121
	X87PRC		fft80PPRO, 1848
	DD			9, 2, 1, 1, 0
	DD	2211,	96,	0.0000141
	X87PRC		fft96PPRO, 2208
	DD			11, 2, 1, 1, 0
	DD	2585,	112,	0.0000179
	X87PRC		fft112PPRO, 2616
	DD			13, 3, 1, 1, 0
	DD	2953,	128,	0.0000178
	X87PRC		fft128PPRO, 2832
	DD			15, 3, 1, 1, 0
	DD	3663,	160,	0.0000296
	X87PRC		fft160PPRO, 3840
	DD			19, 4, 2, 1, 0
	DD	4359,	192,	0.000035
	X87PRC		fft192PPRO, 4608
	DD			23, 5, 2, 1, 0
	DD	5093,	224,	0.000045
	X87PRC		fft224PPRO, 5424
	DD			27, 6, 3, 1, 0
	DD	5833,	256,	0.000045
	X87PRC		fft256PPRO, 5712
	DD			31, 7, 3, 1, 0
	DD	7243,	320,	0.000062
	X87PRC		fft320PPRO, 7680
	DD			39, 9, 2, 1, 0
	DD	8639,	384,	0.000075
	X87PRC		fft384PPRO, 9216
	DD			47, 11, 2, 1, 0
	DD	10085,	448,	0.000093
	X87PRC		fft448PPRO, 10800
	DD			55, 13, 3, 1, 0
	DD	11537,	512,	0.000097
	X87PRC		fft512PPRO, 11472
	DD			63, 15, 3, 1, 0
	DD	14301,	640,	0.000140
	X87PRC		fft640PPRO, 15552
	DD			79, 19, 4, 2, 0
	DD	17047,	768,	0.000167
	X87PRC		fft768PPRO, 18672
	DD			95, 23, 5, 2, 0
	DD	19881,	896,	0.000210
	X87PRC		fft896PPRO, 21840
	DD			111, 27, 6, 3, 0
	DD	22799,	1024,	0.000218
	X87PRC		fft1024PPRO, 22992
	DD			127, 31, 7, 3, 0
	DD	28295,	1280,	0.000302
	X87PRC		fft1280PPRO, 31152
	DD			159, 39, 9, 2, 0
	DD	33761,	1536,	0.000365
	X87PRC		fft1536PPRO, 37392
	DD			191, 47, 11, 2, 0
	DD	39411,	1792,	0.000456
	X87PRC		fft1792PPRO, 43680
	DD			223, 55, 13, 3, 0
	DD	45061,	2048,	0.000490
	X87PRC		fft2048PPRO, 46032
	DD			255, 63, 15, 3, 0
	DD	55825,	2560,	0.000708
	X87PRC		fft2560PPRO, 62544
	DD			319, 79, 19, 4, 2, 0
	DD	66519,	3072,	0.000851
	X87PRC		fft3072PPRO, 75072
	DD			383, 95, 23, 5, 2, 0
	DD	77599,	3584,	0.00107
	X87PRC		fft3584PPRO, 87648
	DD			447, 111, 27, 6, 3, 0
	DD	89047,	4096,	0.00113
	X87PRC		fft4096PPRO, 92112
	DD			511, 127, 31, 7, 3, 0
	DD	110400,	5120,	0.00152
	X87PRC		fft5120PPRO, 24848, 8, 2
	DD			0
	DD	131100,	6144,	0.00191
	X87PRC		fft6144PPRO, 28016, 8, 2
	DD			0
	DD	152800,	7168,	0.00226
	X87PRC		fft7168PPRO, 31232, 8, 2
	DD			0
	DD	175300,	8192,	0.00242
	X87PRC		fft8192PPRO, 34400, 8, 2
	DD			0
	DD	217700,	10240,	0.00333
	X87PRC		fft10KPPRO, 40880, 8, 2
	DD			0
	DD	258200,	12288,	0.00397
	X87PRC		fft12KPPRO, 47264, 8, 2
	DD			0
	DD	301400,	14336,	0.00488
	X87PRC		fft14KPPRO, 53696, 8, 2
	DD			0
	DD	346100,	16384,	0.00522
	X87PRC		fft16KPPRO, 59936, 8, 2
	DD			0
	DD	430300,	20480,	0.00692
	X87PRC		fft20KPPRO, 72800, 8, 2
	DD			0
	DD	511600,	24576,	0.00826
	X87PRC		fft24KPPRO, 85568, 8, 2
	DD			0
	DD	596100,	28672,	0.0101
	X87PRC		fft28KPPRO, 98384, 8, 2
	DD			0
	DD	683700,	32768,	0.0109
	X87PRC		fft32KPPRO, 111008, 8, 2
	DD			0
	DD	848800,	40960,	0.0151
	X87PRC		fft40KP3, 136832, 8, 2
	X87PRC		fft40KPPRO, 136832, 8, 2
	DD			0
	DD	1009000, 49152,	0.0184
	X87PRC		fft48KP3, 162416, 8, 2
	X87PRC		fft48KPPRO, 162416, 8, 2
	DD			0
	DD	1177000, 57344,	0.0227
	X87PRC		fft56KP3, 188048, 8, 2
	X87PRC		fft56KPPRO, 188048, 8, 2
	DD			0
	DD	1350000, 65536,	0.0252
	X87PRC		fft64KP3, 213152, 8, 1
	X87PRC		fft64KPPRO, 213152, 8, 1
	DD			0
	DD	1678000, 81920, 0.0360
	X87PRC		fft80KP3, 264752, 8, 1
	X87PRC		fft80KPPRO, 264752, 8, 1
	DD			0
	DD	1994000, 98304, 0.0445
	X87PRC		fft96KP3, 297536, 10, 1
	X87PRC		fft96KPPRO, 297536, 10, 1
	DD			0
	DD	2324000, 114688, 0.0548
	X87PRC		fft112KP3, 341072, 10, 1
	X87PRC		fft112KPPRO, 341072, 10, 1
	DD			0
	DD	2664000, 131072, 0.0604
	X87PRC		fft128KP3, 384416, 10, 1
	X87PRC		fft128KPPRO, 384416, 10, 1
	DD			0
	DD	3310000, 163840, 0.0830
	X87PRC		fft160KP3, 471680, 10, 1
	X87PRC		fft160KPPRO, 471680, 10, 1
	DD			0
	DD	3933000, 196608, 0.0982
	X87PRC		fft192KP3, 558704, 10, 1
	X87PRC		fft192KPPRO, 558704, 10, 1
	DD			0
	DD	4593000, 229376, 0.1193
	X87PRC		fft224KP3, 645776, 10, 1
	X87PRC		fft224KPPRO, 645776, 10, 1
	DD			0
	DD	5264000, 262144, 0.1316
	X87PRC		fft256KP3, 732320, 10, 1
	X87PRC		fft256KPPRO, 732320, 10, 1
	DD			0
	DD	6545000, 327680, 0.1726
	X87PRC		fft320KP3, 906800, 10, 1
	X87PRC		fft320KPPRO, 906800, 10, 1
	DD			0
	DD	7772000, 393216, 0.2107
	X87PRC		fft384KP3, 1080848, 10, 1
	X87PRC		fft384KPPRO, 1080848, 10, 1
	DD			0
	DD	9071000, 458752, 0.2520
	X87PRC		fft448KP3, 1254944, 10, 1
	X87PRC		fft448KPPRO, 1254944, 10, 1
	DD			0
	DD	10380000, 524288, 0.2808
	X87PRC		fft512KP3, 1428128, 10, 1
	X87PRC		fft512KPPRO, 1428128, 10, 1
	DD			0
	DD	12890000, 655360, 0.372
	X87PRC		fft640KP3, 1777232, 10, 1
	X87PRC		fft640KPPRO, 1777232, 10, 1
	DD			0
	DD	15310000, 786432, 0.453
	X87PRC		fft768KP3, 2125376, 10, 1
	X87PRC		fft768KPPRO, 2125376, 10, 1
	DD			0
	DD	17890000, 917504, 0.536
	X87PRC		fft896KP3, 2473568, 10, 1
	X87PRC		fft896KPPRO, 2473568, 10, 1
	DD			0
	DD	20460000, 1048576, 0.600
	X87PRC		fft1024K2P3, 2819744, 10, 1
	X87PRC		fft1024K2PPRO, 2819744, 10, 1
	DD			0
	DD	25390000, 1310720, 0.776
	X87PRC		fft1280KP3, 3474992, 12, 1
	X87PRC		fft1280KPPRO, 3474992, 12, 1
	DD			0
	DD	30190000, 1572864, 0.934
	X87PRC		fft1536KP3, 4140560, 12, 1
	X87PRC		fft1536KPPRO, 4140560, 12, 1
	DD			0
	DD	35200000, 1835008, 1.113
	X87PRC		fft1792KP3, 4806176, 12, 1
	X87PRC		fft1792KPPRO, 4806176, 12, 1
	DD			0
	DD	40300000, 2097152, 1.226
	X87PRC		fft2048K2P3, 5470880, 12, 1
	X87PRC		fft2048K2PPRO, 5470880, 12, 1
	DD			0
	DD	50020000, 2621440, 1.636
	X87PRC		fft2560K2P3, 6803024, 12, 1
	X87PRC		fft2560K2PPRO, 6803024, 12, 1
	DD			0
	DD	59510000, 3145728, 1.990
	X87PRC		fft3072K2P3, 8134208, 12, 1
	X87PRC		fft3072K2PPRO, 8134208, 12, 1
	DD			0
	DD	69360000, 3670016, 2.380
	X87PRC		fft3584K2P3, 9465440, 12, 1
	X87PRC		fft3584K2PPRO, 9465440, 12, 1
	DD			0
	DD	79370000, 4194304, 2.604
	X87PRC		fft4096K2P3, 10794656, 12, 1
	X87PRC		fft4096K2PPRO, 10794656, 12, 1
	DD			0
	DD	0
jmptablep DD	755,	32,	0.000004
	X87PRC		fft32pPPRO, 976
	DD			4, 1, 1, 1, 0
	DD	1111,	48,	0.000007
	X87PRC		fft48pPPRO, 1608
	DD			6, 3, 1, 1, 0
	DD	1485,	64,	0.000010
	X87PRC		fft64pPPRO, 1952
	DD			8, 4, 1, 1, 0
	DD	2199,	96,	0.0000141
	X87PRC		fft96pPPRO, 3072
	DD			12, 3, 1, 1, 0
	DD	2947,	128,	0.000021
	X87PRC		fft128pPPRO, 3904
	DD			16, 4, 1, 1, 1, 0
	DD	4345,	192,	0.000035
	X87PRC		fft192pPPRO, 6288
	DD			24, 6, 3, 1, 0
	DD	5817,	256,	0.000051
	X87PRC		fft256pPPRO, 7808
	DD			32, 8, 4, 1, 1, 0
	DD	8607,	384,	0.000075
	X87PRC		fft384pPPRO, 12432
	DD			48, 12, 3, 1, 0
	DD	11515,	512,	0.000106
	X87PRC		fft512pPPRO, 15616
	DD			64, 16, 4, 1, 1, 0
	DD	17001,	768,	0.000167
	X87PRC		fft768pPPRO, 25008
	DD			96, 24, 6, 3, 1, 0
	DD	22701,	1024,	0.000249
	X87PRC		fft1024pPPRO, 31232
	DD			128, 32, 8, 4, 1, 0
	DD	33569,	1536,	0.000365
	X87PRC		fft1536pPPRO, 49872
	DD			192, 48, 12, 3, 0
	DD	44951,	2048,	0.000582
	X87PRC		fft2048pPPRO, 62464
	DD			256, 64, 16, 4, 1, 0
	DD	66319,	3072,	0.000851
	X87PRC		fft3072pPPRO, 99888
	DD			384, 96, 24, 6, 3, 0
	DD	88747,	4096,	0.00135
	X87PRC		fft4096pPPRO, 124928
	DD			512, 128, 32, 8, 4, 0
	DD	130600,	6144,	0.00191
	X87PRC		fft6144pPPRO, 26512, 8, 2
	DD			0
	DD	174000,	8192,	0.00284
	X87PRC		fft8192pPPRO, 32960, 8, 2
	DD			0
	DD	257700,	12288,	0.00397
	X87PRC		fft12KpPPRO, 46000, 8, 2
	DD			0
	DD	344700,	16384,	0.00588
	X87PRC		fft16KpPPRO, 58752, 8, 2
	DD			0
	DD	508600,	24576,	0.00826
	X87PRC		fft24KpPPRO, 84688, 8, 2
	DD			0
	DD	679400,	32768,	0.01299
	X87PRC		fft32KpPPRO, 110336, 8, 2
	DD			0
	DD	1006000, 49152,	0.0184
	X87PRC		fft48KpP3, 162352, 8, 2
	X87PRC		fft48KpPPRO, 162352, 8, 2
	DD			0
	DD	1345000, 65536,	0.03283
	X87PRC		fft64KpP3, 213504, 8, 1
	X87PRC		fft64KpPPRO, 213504, 8, 1
	DD			0
	DD	1983000, 98304, 0.0445
	X87PRC		fft96KpP3, 290512, 10, 1
	X87PRC		fft96KpPPRO, 290512, 10, 1
	DD			0
	DD	2652000, 131072, 0.0719
	X87PRC		fft128KpP3, 377600, 10, 1
	X87PRC		fft128KpPPRO, 377600, 10, 1
	DD			0
	DD	3924000, 196608, 0.0982
	X87PRC		fft192KpP3, 552496, 10, 1
	X87PRC		fft192KpPPRO, 552496, 10, 1
	DD			0
	DD	5242000, 262144, 0.155
	X87PRC		fft256KpP3, 726528, 10, 1
	X87PRC		fft256KpPPRO, 726528, 10, 1
	DD			0
	DD	7733000, 393216, 0.2107
	X87PRC		fft384KpP3, 1076176, 10, 1
	X87PRC		fft384KpPPRO, 1076176, 10, 1
	DD			0
	DD	10320000, 524288, 0.322
	X87PRC		fft512KpP3, 1424384, 10, 1
	X87PRC		fft512KpPPRO, 1424384, 10, 1
	DD			0
	DD	15260000, 786432, 0.453
	X87PRC		fft768KpP3, 2123824, 10, 1
	X87PRC		fft768KpPPRO, 2123824, 10, 1
	DD			0
	DD	20360000, 1048576, 0.681
	X87PRC		fft1024KpP3, 2820096, 10, 1
	X87PRC		fft1024KpPPRO, 2820096, 10, 1
	DD			0
	DD	30070000, 1572864, 0.934
	X87PRC		fft1536KpP3, 4111312, 12, 1
	X87PRC		fft1536KpPPRO, 4111312, 12, 1
	DD			0
	DD	40110000, 2097152, 1.380
	X87PRC		fft2048KpP3, 5442560, 12, 1
	X87PRC		fft2048KpPPRO, 5442560, 12, 1
	DD			0
	DD	59360000, 3145728, 1.990
	X87PRC		fft3072Kp2P3, 8108080, 12, 1
	X87PRC		fft3072Kp2PPRO, 8108080, 12, 1
	DD			0
	DD	79100000, 4194304, 2.919
	X87PRC		fft4096Kp2P3, 10770432, 12, 1
	X87PRC		fft4096Kp2PPRO, 10770432, 12, 1
	DD			0
	DD	0
ENDIF

;; Jump tables for the Pentium 4 SSE2 optimized code

xjmptable DD	0
	org	$-4
	PRCSTRT	743,	32,	0.00000111
	PRCENTRY		xfft_hg_32_op_BLEND, 896
	DD			4, 4
	DD			1, 1, 1, 1, 1, 0
	PRCSTRT	1099,	48,	0.00000144
	PRCENTRY		xfft_hg_48_op_BLEND, 1408 
	DD			6, 6
	DD			2, 1, 1, 1, 1, 0
	PRCSTRT	1469,	64,	0.00000178
	PRCENTRY		xfft_hg_64_op_BLEND, 1920 
	DD			8, 8
	DD			3, 1, 1, 1, 1, 0
	PRCSTRT	1827,	80,	0.00000222
	PRCENTRY		xfft_hg_80_op_BLEND, 2176 
	DD			10, 8*2048+2
	DD			4, 2, 1, 1, 1, 0
	PRCSTRT	2179,	96,	0.00000259
	PRCENTRY		xfft_hg_96_op_BLEND, 2432 
	DD			12, 12
	DD			5, 2, 1, 1, 1, 0
	PRCSTRT	2539,	112,	0.00000311
	PRCENTRY		xfft_hg_112_op_BLEND, 2944
	DD			14, (8*2048+4)*2048+2
	DD			6, 3, 1, 1, 1, 0
	PRCSTRT	2905,	128,	0.00000319
	PRCENTRY		xfft_hg_128_op_BLEND, 3328
	DD			16, 16
	DD			7, 3, 1, 1, 1, 0
	PRCSTRT	3613,	160,	0.00000450
	PRCENTRY		xfft_hg_160_op_BLEND, 4736
	DD			20, 16*2048+4
	DD			9, 9, 2, 1, 4, 0
	PRCSTRT	4311,	192,	0.00000542
	PRCENTRY		xfft_hg_192_op_BLEND, 5632
	DD			24, 24
	DD			11, 11, 2, 1, 5, 0
	PRCSTRT	5029,	224,	0.00000663
	PRCENTRY		xfft_hg_224_op_BLEND, 6656
	DD			28, (16*2048+8)*2048+4
	DD			13, 13, 3, 1, 6, 0
	PRCSTRT	5755,	256,	0.00000691
	PRCENTRY		xfft_hg_256_op_BLEND, 7296
	DD			32, 32
	DD			15, 15, 3, 1, 7, 0
	PRCSTRT	7149,	320,	0.00000928
	PRCENTRY		xfft_hg_320_op_BLEND, 8448
	DD			40, 32*2048+8
	DD			19, 9, 2, 4, 1, 0
	PRCSTRT	8527,	384,	0.0000111
	PRCENTRY		xfft_hg_384_op_BLEND, 9984
	DD			48, 48
	DD			23, 11, 2, 5, 1, 0
	PRCSTRT	9933,	448,	0.0000133
	PRCENTRY		xfft_hg_448_op_BLEND, 11648
	DD			56, (32*2048+16)*2048+8
	DD			27, 13, 3, 6, 1, 0
	PRCSTRT	11359,	512,	0.0000143
	PRCENTRY		xfft_hg_512_op_BLEND, 13056
	DD			64, 64
	DD			31, 15, 3, 7, 1, 0
	PRCSTRT	14119,	640,	0.0000215
	PRCENTRY		xfft_hg_640_op_BLEND, 17408
	DD			80, 64*2048+16
	DD			39, 19, 9, 9, 4*256+2, 0
	PRCSTRT	16839,	768,	0.0000260
	PRCENTRY		xfft_hg_768_op_BLEND, 20736
	DD			96, 96
	DD			47, 23, 11, 11, 5*256+2, 0
	PRCSTRT	19639,	896,	0.0000321
	PRCENTRY		xfft_hg_896_op_BLEND, 24448
	DD			112, (64*2048+32)*2048+16
	DD			55, 27, 13, 13, 6*256+3, 0
	PRCSTRT	22477,	1024,	0.0000349
	PRCENTRY		xfft_hg_1024_op_BLEND, 26112
	DD			128, 128
	DD			63, 31, 15, 15, 7*256+3, 0
	PRCSTRT	27899,	1280,	0.0000494
	PRCENTRY		xfft_hg_1280_op_BLEND, 33664
	DD			160, 128*2048+32
	DD			79, 39, 9, 4*256+19, 2, 0
	PRCSTRT	33289,	1536,	0.0000601
	PRCENTRY		xfft_hg_1536_op_BLEND, 40320
	DD			192, 192
	DD			95, 47, 11, 5*256+23, 2, 0
	PRCSTRT	38799,	1792,	0.0000719
	PRCENTRY		xfft_hg_1792_op_BLEND, 47232
	DD			224, (128*2048+64)*2048+32
	DD			111, 55, 13, 6*256+27, 3, 0
	PRCSTRT	44339,	2048,	0.0000773
	PRCENTRY		xfft_hg_2048_op_BLEND, 52224
	DD			256, 256
	DD			127, 63, 15, 7*256+31, 3, 0
	PRCSTRT	55099,	2560,	0.000111
	PRCENTRY		xfft_hg_2560_op_BLEND, 68096
	DD			320, 256*2048+64
	DD			159, 79, 9*256+19, 9*256+39, 4*256+2, 0
	PRCSTRT	65729,	3072,	0.000131
	PRCENTRY		xfft_hg_3072_op_BLEND, 81792
	DD			384, 384
	DD			191, 95, 11*256+23, 11*256+47, 5*256+2, 0
	PRCSTRT	76559,	3584,	0.000165
	PRCENTRY		xfft_hg_3584_op_BLEND, 95488
	DD			448, (256*2048+128)*2048+64
	DD			223, 111, 13*256+27, 13*256+55, 6*256+3, 0
	PRCSTRT	87549,	4096,	0.000163
	PRCENTRY		xfft_hg_4096_op_BLEND, 104448
	DD			512, 512
	DD			255, 127, 15*256+31, 15*256+63, 7*256+3, 0
	PRCSTRT	108800,	5120,	0.000215
	PRCENTRY		xfft_hg_5120_op_BLEND, 135296
	DD			640, 512*2048+128
	DD			319, 159, 9*256+39, 19*256+79, 4*256+2, 0
	PRCSTRT	129900,	6144,	0.000276
	PRCENTRY		xfft_hg_6144_op_BLEND, 162432
	DD			768, 768
	DD			383, 191, 11*256+47, 23*256+95, 5*256+2, 0
	PRCSTRT	151300,	7168,	0.000374
	PRCENTRY2A		xfft_hg_7168_op_BLEND, 189568, I7_32 + P4_1024_32
	DD			896, (512*2048+256)*2048+128
	DD			447, 223, 13*256+55, 27*256+111, 6*256+3, 0
	PRCSTRT	172800,	8192,	0.000398
	PRCENTRY2		xfft_r4_8K_8_4, 81408,
	PRCENTRY2		xfft_r4_8K_np_8_4, 81408, P4_1024 + P4TP_512 + P4TP_256 + I7_32 + CORE2 + K8 + K10
	PRCENTRY2A		xfft_r4_8K_np_8_4_P4, 81408, CORE2_512
	PRCENTRY2A		xfft_hg_8192_op_BLEND, 208896, I7_64
	DD			1024, 1024
	DD			511, 255, 15*256+63, 31*256+127, 7*256+3, 0
	PRCSTRT	214400,	10240,	0.000470
	PRCENTRY2		xfft_hg_10K_ip_8_4, 59904, P4_1024 + P4TP_512 + P4TP_256 + I7_64 + CORE2 + K8 + K10
	PRCENTRY2A		xfft_hg_10K_ip_8_4_P4, 59904, I7_32 + CORE2_512
	DD			0
	PRCSTRT	256000,	12288,	0.000590
	PRCENTRY2		xfft_hg_12K_ip_8_4, 69632, P4_1024 + P4TP_512 + P4TP_256 + I7_64 + CORE2 + K8 + K10
	PRCENTRY2A		xfft_hg_12K_ip_8_4_P4, 69632, I7_32 + CORE2_512
	DD			0
	PRCSTRT	297700,	14336,	0.000716
	PRCENTRY2		xfft_hg_14K_ip_8_4, 79488, P4_1024 + P4TP_256 + I7_64 + K8_32
	PRCENTRY2A		xfft_hg_14K_ip_8_4_P4, 79488, CORE2_512_64
	DD			0
	PRCSTRT	340600,	16384,	0.000787
	PRCENTRY2		xfft_r4_16K_8_4, 168960, P4TP_256
	PRCENTRY2A		xfft_r4_16K_8_4_P4, 168960, CORE2_512_64
	PRCENTRY2		xfft_r4_16K_np_8_4, 168960, P4_1024_64 + P4TP_512 + I7 + CORE2 + K8 + K10
	PRCENTRY2A		xfft_r4_16K_np_8_4_P4, 168960, CORE2_512_32
	PRCENTRY2		xfft_hg_16K_ip_8_4, 89088, P4_1024_32
	DD			0
	PRCSTRT	424300,	20480,	0.00103
	PRCENTRY2		xfft_r4_20K_8_4, 188160, P4TP_256 + CORE2_512_64
	PRCENTRY2A		xfft_r4_20K_8_4_P4, 188160, I7_32
	PRCENTRY2		xfft_r4_20K_np_8_4, 188160, P4_1024_64 + P4TP_512 + I7_64 + CORE2 + K8 + K10
	PRCENTRY2A		xfft_r4_20K_np_8_4_P4, 188160, CORE2_512_32
	PRCENTRY2		xfft_hg_20K_ip_8_4, 107904, P4_1024_32
	DD			0
	PRCSTRT	506900,	24576,	0.00132
	PRCENTRY2		xfft_r4_24K_768_4, 241152, P4TP_256
	PRCENTRY2		xfft_r4_24K_np_768_4, 241152, K8_64
	PRCENTRY2		xfft_r4_24K_8_4, 223744, I7_64 + CORE2_512_64
	PRCENTRY2A		xfft_r4_24K_8_4_P4, 223744, I7_32
	PRCENTRY2		xfft_r4_24K_np_8_4, 223744, P4TP_512 + CORE2 + K8_32 + K10
	PRCENTRY2A		xfft_r4_24K_np_8_4_P4, 223744, CORE2_512_32
	PRCENTRY2		xfft_hg_24K_ip_8_4, 127232, P4_1024
	DD			0
	PRCSTRT	590600,	28672,	0.00156
	PRCENTRY2		xfft_r4_28K_8_4, 259328
	PRCENTRY2A		xfft_r4_28K_8_4_P4, 259328, I7
	PRCENTRY2		xfft_r4_28K_8_2, 259328, CORE2_512_32 + K10
	PRCENTRY2		xfft_r4_28K_np_8_4, 259328, CORE2_32 + K8_32
	PRCENTRY2		xfft_r4_28K_np_8_2, 259328, K8_64
	PRCENTRY2		xfft_hg_28K_ip_8_4, 146688, P4TP_256
	PRCENTRY2A		xfft_hg_28K_ip_8_4_CORE, 146688, P4_1024
	PRCENTRY2A		xfft_hg_28K_ip_8_4_P4, 146688, P4TP_512
	DD			0
	PRCSTRT	673100,	32768,	0.00175
	PRCENTRY2		xfft_r4dwpn_32K_8_4, 211456, P4TP_512 + P4TP_256
	PRCENTRY2A		xfft_r4dwpn_32K_8_4_P4, 211456, I7
	PRCENTRY2		xfft_r4dwpn_32K_np_8_4, 211456, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_32K_8_2, 211456, CORE2_512_64 + K8
	PRCENTRY2A		xfft_r4dwpn_32K_8_2_P4, 211456, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_32K_np_8_2, 211456
	PRCENTRY2		xfft_r4delay_32K_8_4, 173568
	PRCENTRY2		xfft_r4delay_32K_np_8_4, 173568, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_32K_8_2, 173568, CORE2_512_64 + K8
	PRCENTRY2A		xfft_r4delay_32K_8_2_P4, 173568, CORE2_512_32
	PRCENTRY2		xfft_r4delay_32K_np_8_2, 173568
	PRCENTRY2		xfft_r4_32K_10_4, 325120, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_r4_32K_np_10_4, 325120
	PRCENTRY2		xfft_r4_32K_8_4, 311296
	PRCENTRY2A		xfft_r4_32K_8_4_P4, 311296, I7
	PRCENTRY2		xfft_r4_32K_np_8_4, 311296, CORE2_64
	PRCENTRY2		xfft_hg_32K_ip_8_4, 165888, P4_1024
	DD			0
	PRCSTRT	836800,	40960,	0.00225
	PRCENTRY2		xfft_r4_40K_1280_4, 400896, P4TP_256 + CORE2_512 + K8 + K10
	PRCENTRY2A		xfft_r4_40K_1280_4_P4, 400896, I7_32
	PRCENTRY2		xfft_r4_40K_np_1280_4, 400896, P4_1024_64 + CORE2
	PRCENTRY2A		xfft_r4_40K_np_1280_4_P4, 400896, I7_64
	PRCENTRY2		xfft_hg_40K_ip_11_4, 189696, P4_1024_32 + P4TP_512
	DD			0
	PRCSTRT	999900, 49152,	0.00279
	PRCENTRY2		xfft_r4_48K_768_4, 500736, P4_1024 + K8 + K10
	PRCENTRY2		xfft_r4_48K_np_768_4, 500736, I7 + CORE2
	PRCENTRY2		xfft_hg_48K_ip_11_4, 206208, P4TP_512 + P4TP_256 + CORE2_512
	DD			0
	PRCSTRT	1164000, 57344,	0.00327
	PRCENTRY2		xfft_r4_56K_8_4, 524800, P4_1024_32
	PRCENTRY2A		xfft_r4_56K_8_4_P4, 524800, I7
	PRCENTRY2A		xfft_r4_56K_8_4_CORE, 524800, P4_1024_64
	PRCENTRY2		xfft_r4_56K_np_8_4, 524800
	PRCENTRY2A		xfft_r4_56K_np_8_4_P4, 524800, CORE2_32
	PRCENTRY2		xfft_r4_56K_8_2, 524800, K10
	PRCENTRY2		xfft_r4_56K_8_1, 524800, K8
	PRCENTRY2		xfft_hg_56K_ip_11_4, 222976, P4TP_512 + P4TP_256
	DD			0
	PRCSTRT	1245000, 61440,	0.00344
	PRCENTRY2		xfft_r4_60K_768_4, 556800
	PRCENTRY2		xfft_r4_60K_np_768_4, 556800
	PRCENTRY2		xfft_r4_60K_768_2, 556800, K8_32
	PRCENTRY2		xfft_r4_60K_768_1, 556800
	DD			0
	PRCSTRT	1327000, 65536,	0.00367
	PRCENTRY2		xfft_r4dwpn_64K_8_4, 411648, P4_1024_32 + P4TP_256
	PRCENTRY2A		xfft_r4dwpn_64K_8_4_P4, 411648, I7
	PRCENTRY2A		xfft_r4dwpn_64K_8_4_CORE, 411648, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_64K_np_8_4, 411648, CORE2_64
	PRCENTRY2A		xfft_r4dwpn_64K_np_8_4_P4, 411648, CORE2_32
	PRCENTRY2		xfft_r4dwpn_64K_8_2, 411648, K10
	PRCENTRY2		xfft_r4dwpn_64K_np_8_2, 411648
	PRCENTRY2		xfft_r4dwpn_64K_8_1, 411648, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_64K_np_8_1, 411648
	PRCENTRY2		xfft_r4delay_64K_8_4, 336896, P4_1024_32
	PRCENTRY2A		xfft_r4delay_64K_8_4_P4, 336896, I7
	PRCENTRY2A		xfft_r4delay_64K_8_4_CORE, 336896, P4_1024_64
	PRCENTRY2		xfft_r4delay_64K_np_8_4, 336896, CORE2_64
	PRCENTRY2A		xfft_r4delay_64K_np_8_4_P4, 336896, CORE2_32
	PRCENTRY2		xfft_r4delay_64K_8_2, 336896, CORE2_512_32 + K10
	PRCENTRY2		xfft_r4delay_64K_np_8_2, 336896
	PRCENTRY2		xfft_r4delay_64K_8_1, 336896, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4delay_64K_np_8_1, 336896
	PRCENTRY2		xfft_r4_64K_10_4, 670720
	PRCENTRY2		xfft_r4_64K_np_10_4, 670720
	PRCENTRY2		xfft_r4_64K_8_4, 628736
	PRCENTRY2		xfft_r4_64K_np_8_4, 628736
	PRCENTRY2		xfft_hg_64K_ip_11_4, 239488, P4TP_512 + P4TP_256
	DD			0
	PRCSTRT	1486000, 73728,	0.00422
	PRCENTRY2		xfft_r4_72K_768_4, 662016, P4_1024_32 + P4TP_256 + K10
	PRCENTRY2A		xfft_r4_72K_768_4_P4, 662016, I7_32
	PRCENTRY2		xfft_r4_72K_np_768_4, 662016, I7_64 + CORE2
	PRCENTRY2		xfft_r4_72K_768_2, 662016, CORE2_512 + K8
	PRCENTRY2		xfft_r4_72K_768_1, 662016
	DD			0
	PRCSTRT	1653000, 81920, 0.00474
	PRCENTRY2		xfft_r4dwpn_80K_8_4, 335872, P4_1024_64 + K10_32
	PRCENTRY2A		xfft_r4dwpn_80K_8_4_P4, 335872, I7
	PRCENTRY2A		xfft_r4dwpn_80K_8_4_CORE, 335872, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_80K_np_8_4, 335872, CORE2
	PRCENTRY2		xfft_r4dwpn_80K_8_2, 335872, P4TP_512 + P4TP_256 + K10_64
	PRCENTRY2		xfft_r4dwpn_80K_8_1, 335872, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_80K_8_4, 308736, P4_1024_64
	PRCENTRY2A		xfft_r4delay_80K_8_4_CORE, 308736, P4_1024_32
	PRCENTRY2		xfft_r4delay_80K_np_8_4, 308736, CORE2
	PRCENTRY2		xfft_r4delay_80K_8_2, 308736, P4TP_512 + P4TP_256 + K8_32 + K10_64
	PRCENTRY2		xfft_r4delay_80K_8_1, 308736, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4_80K_1280_4, 832512
	PRCENTRY2		xfft_r4_80K_np_1280_4, 832512
	PRCENTRY2		xfft_r4_80K_10_4, 745216, K10_32
	PRCENTRY2A		xfft_r4_80K_10_4_P4, 745216, I7
	PRCENTRY2		xfft_r4_80K_np_10_4, 745216
	PRCENTRY2		xfft_hg_80K_ip_11_4, 273408
	DD			0
	PRCSTRT	1731000, 86016, 0.00503
	PRCENTRY2		xfft_r4_84K_768_4, 767232, P4_1024_64 + P4TP_512 + P4TP_256 + I7_64
	PRCENTRY2A		xfft_r4_84K_768_4_P4, 767232, I7_32
	PRCENTRY2A		xfft_r4_84K_768_4_CORE, 767232, P4_1024_32
	PRCENTRY2		xfft_r4_84K_np_768_4, 767232, CORE2
	PRCENTRY2		xfft_r4_84K_768_2, 767232, K10
	PRCENTRY2		xfft_r4_84K_768_1, 767232, K8
	DD			0
	PRCSTRT	1975000, 98304, 0.00584
	PRCENTRY2		xfft_r4dwpn_96K_768_4, 616960, P4TP_512 + I7_64 + K10
	PRCENTRY2A		xfft_r4dwpn_96K_768_4_P4, 616960, I7_32
	PRCENTRY2		xfft_r4dwpn_96K_np_768_4, 616960, CORE2
	PRCENTRY2		xfft_r4dwpn_96K_768_2, 616960, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_96K_768_1, 616960
	PRCENTRY2		xfft_r4dwpn_96K_8_4, 376320
	PRCENTRY2A		xfft_r4dwpn_96K_8_4_CORE, 376320, P4_1024
	PRCENTRY2		xfft_r4dwpn_96K_np_8_4, 376320
	PRCENTRY2		xfft_r4delay_96K_768_4, 505344, P4TP_512 + P4TP_256 + I7_64 + K10
	PRCENTRY2		xfft_r4delay_96K_np_768_4, 505344, CORE2_64
	PRCENTRY2		xfft_r4delay_96K_768_2, 505344, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4delay_96K_768_1, 505344, K8_32
	PRCENTRY2		xfft_r4delay_96K_8_4, 352768
	PRCENTRY2A		xfft_r4delay_96K_8_4_CORE, 352768, P4_1024
	PRCENTRY2		xfft_r4delay_96K_np_8_4, 352768
	PRCENTRY2		xfft_r4_96K_10_4, 885248
	PRCENTRY2		xfft_r4_96K_np_10_4, 885248
	PRCENTRY2		xfft_r4_96K_768_4, 921600
	PRCENTRY2		xfft_r4_96K_np_768_4, 921600, I7_32 + CORE2_32
	PRCENTRY2		xfft_hg_96K_ip_11_4, 306688
	DD			0
	PRCSTRT	2300000, 114688, 0.00693
	PRCENTRY2		xfft_r4dwpn_112K_8_4, 418816, K10_64
	PRCENTRY2A		xfft_r4dwpn_112K_8_4_P4, 418816, I7
	PRCENTRY2A		xfft_r4dwpn_112K_8_4_CORE, 418816, P4_1024
	PRCENTRY2		xfft_r4dwpn_112K_np_8_4, 418816, CORE2
	PRCENTRY2		xfft_r4dwpn_112K_8_2, 418816, P4TP_512 + P4TP_256 + K8
	PRCENTRY2		xfft_r4dwpn_112K_8_1, 418816, CORE2_512_64
	PRCENTRY2		xfft_r4delay_112K_8_4, 398848
	PRCENTRY2A		xfft_r4delay_112K_8_4_CORE, 398848, P4_1024
	PRCENTRY2		xfft_r4delay_112K_np_8_4, 398848, CORE2_64
	PRCENTRY2		xfft_r4delay_112K_8_2, 398848, P4TP_512 + P4TP_256 + K8
	PRCENTRY2		xfft_r4delay_112K_8_1, 398848
	PRCENTRY2		xfft_r4_112K_10_4, 1025280, I7_64 + CORE2_512_64 + K10
	PRCENTRY2A		xfft_r4_112K_10_4_P4, 1025280, I7_32
	PRCENTRY2		xfft_r4_112K_np_10_4, 1025280, CORE2_32
	PRCENTRY2		xfft_hg_112K_ip_11_4, 340096
	DD			0
	PRCSTRT	2457000, 122880, 0.00730
	PRCENTRY2		xfft_r4_120K_1280_4, 1100288
	PRCENTRY2		xfft_r4_120K_np_1280_4, 1100288
	DD			0
	PRCSTRT	2622000, 131072, 0.00779
	PRCENTRY2		xfft_r4dwpn_128K_10_4, 823808, I7_64 + K10
	PRCENTRY2A		xfft_r4dwpn_128K_10_4_P4, 823808, I7_32
	PRCENTRY2		xfft_r4dwpn_128K_np_10_4, 823808, CORE2
	PRCENTRY2		xfft_r4dwpn_128K_10_2, 823808, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_128K_10_1, 823808, K8_32
	PRCENTRY2		xfft_r4dwpn_128K_8_4, 460800
	PRCENTRY2A		xfft_r4dwpn_128K_8_4_CORE, 460800, P4_1024
	PRCENTRY2		xfft_r4dwpn_128K_np_8_4, 460800
	PRCENTRY2		xfft_r4dwpn_128K_8_2, 460800
	PRCENTRY2		xfft_r4dwpn_128K_8_1, 460800
	PRCENTRY2		xfft_r4delay_128K_10_4, 675328, P4TP_512 + P4TP_256 + I7_64 + K10
	PRCENTRY2A		xfft_r4delay_128K_10_4_P4, 675328, I7_32
	PRCENTRY2		xfft_r4delay_128K_np_10_4, 675328, CORE2_64
	PRCENTRY2		xfft_r4delay_128K_10_2, 675328, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4delay_128K_10_1, 675328, K8_32
	PRCENTRY2		xfft_r4delay_128K_8_4, 444416
	PRCENTRY2A		xfft_r4delay_128K_8_4_CORE, 444416, P4_1024
	PRCENTRY2		xfft_r4delay_128K_np_8_4, 444416, CORE2_32
	PRCENTRY2		xfft_r4delay_128K_8_2, 444416
	PRCENTRY2		xfft_r4delay_128K_8_1, 444416
	PRCENTRY2		xfft_hg_128K_ip_11_4, 373248
	DD			0
	PRCSTRT	2937000, 147456, 0.00870
	PRCENTRY2		xfft_r4_144K_1536_4, 1319424, P4TP_256 + I7_64
	PRCENTRY2A		xfft_r4_144K_1536_4_P4, 1319424, CORE2_32 + I7_32
	PRCENTRY2A		xfft_r4_144K_1536_4_CORE, 1319424, P4_1024
	PRCENTRY2		xfft_r4_144K_1536_2, 1319424, CORE2_512_32 + K8 + K10
	PRCENTRY2		xfft_r4_144K_1536_1, 1319424
	DD			0
	PRCSTRT	3260000, 163840, 0.00914
	PRCENTRY2		xfft_r4dwpn_160K_1280_4, 1022464, P4TP_256 + K10_32
	PRCENTRY2		xfft_r4dwpn_160K_1280_2, 1022464, CORE2_512 + K8 + K10_64
	PRCENTRY2		xfft_r4dwpn_160K_1280_1, 1022464
	PRCENTRY2		xfft_r4dwpn_160K_8_4, 652800, CORE2_64
	PRCENTRY2A		xfft_r4dwpn_160K_8_4_P4, 652800, CORE2_32 + I7
	PRCENTRY2A		xfft_r4dwpn_160K_8_4_CORE, 652800, P4_1024
	PRCENTRY2		xfft_r4delay_160K_1280_4, 837120, K10
	PRCENTRY2		xfft_r4delay_160K_1280_2, 837120, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_160K_1280_1, 837120
	PRCENTRY2		xfft_r4delay_160K_8_4, 607232, CORE2_64
	PRCENTRY2A		xfft_r4delay_160K_8_4_P4, 607232, CORE2_32 + I7
	PRCENTRY2A		xfft_r4delay_160K_8_4_CORE, 607232, P4_1024
	PRCENTRY2		xfft_hg_160K_10_4, 471424, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_hg_160K_10_2, 471424
	DD			0
	PRCSTRT	3897000, 196608, 0.0114
	PRCENTRY2		xfft_r4dwpn_192K_768_4, 1210368, P4TP_256 + I7_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4dwpn_192K_768_4_P4, 1210368, CORE2_32 + I7_32
	PRCENTRY2		xfft_r4dwpn_192K_768_2, 1210368, K8_64
	PRCENTRY2		xfft_r4dwpn_192K_768_1, 1210368, CORE2_512 + K8_32
	PRCENTRY2		xfft_r4dwpn_192K_8_4, 731136
	PRCENTRY2A		xfft_r4dwpn_192K_8_4_CORE, 731136, P4_1024
	PRCENTRY2		xfft_r4delay_192K_768_4, 979968, I7_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4delay_192K_768_4_P4, 979968, CORE2_32 + I7_32
	PRCENTRY2		xfft_r4delay_192K_768_2, 979968
	PRCENTRY2		xfft_r4delay_192K_768_1, 979968, K8
	PRCENTRY2		xfft_r4delay_192K_8_4, 695296, P4_1024_64
	PRCENTRY2		xfft_r4_192K_768_4, 1861632
	PRCENTRY2		xfft_hg_192K_10_4, 555392, P4TP_512
	PRCENTRY2		xfft_hg_192K_10_2, 555392, P4TP_256 + CORE2_512
	PRCENTRY2		xfft_hg_192K_10_1, 555392, P4_1024_32
	DD			0
	PRCSTRT	4538000, 229376, 0.0134
	PRCENTRY2		xfft_r4dwpn_224K_8_4, 813568, CORE2_64 + K10
	PRCENTRY2A		xfft_r4dwpn_224K_8_4_P4, 813568, CORE2_32 + I7
	PRCENTRY2A		xfft_r4dwpn_224K_8_4_CORE, 813568, P4_1024
	PRCENTRY2		xfft_r4dwpn_224K_8_2, 813568
	PRCENTRY2		xfft_r4dwpn_224K_8_1, 813568, K8
	PRCENTRY2		xfft_r4delay_224K_8_4, 787456, CORE2_64 + K10_32
	PRCENTRY2A		xfft_r4delay_224K_8_4_P4, 787456, CORE2_32 + I7
	PRCENTRY2A		xfft_r4delay_224K_8_4_CORE, 787456, P4_1024
	PRCENTRY2		xfft_r4delay_224K_8_2, 787456, K10_64
	PRCENTRY2		xfft_r4delay_224K_8_1, 787456
	PRCENTRY2		xfft_hg_224K_10_4, 639616, P4TP_512
	PRCENTRY2		xfft_hg_224K_10_2, 639616, P4TP_256 + K8
	DD			0
	PRCSTRT	4846000, 245760, 0.0141
	PRCENTRY2		xfft_r4dwpn_240K_768_4, 962560, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2A		xfft_r4dwpn_240K_768_4_P4, 962560, CORE2_32
	PRCENTRY2A		xfft_r4dwpn_240K_768_4_CORE, 962560, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_240K_768_2, 962560, P4TP_256 + K8_64
	PRCENTRY2		xfft_r4dwpn_240K_768_1, 962560, K8_32
	PRCENTRY2		xfft_r4delay_240K_768_4, 886272, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2A		xfft_r4delay_240K_768_4_P4, 886272, CORE2_32
	PRCENTRY2A		xfft_r4delay_240K_768_4_CORE, 886272, P4_1024_32
	PRCENTRY2		xfft_r4delay_240K_768_2, 886272, P4TP_256
	PRCENTRY2		xfft_r4delay_240K_768_1, 886272, K8
	DD			0
	PRCSTRT	5168000, 262144, 0.0150
	PRCENTRY2		xfft_r4dwpn_256K_10_4, 1613824, P4TP_256 + I7_64 + K10_64
	PRCENTRY2A		xfft_r4dwpn_256K_10_4_P4, 1613824, I7_32
	PRCENTRY2		xfft_r4dwpn_256K_10_2, 1613824, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_256K_10_1, 1613824, K8_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_256K_8_4, 894976, CORE2
	PRCENTRY2A		xfft_r4dwpn_256K_8_4_CORE, 894976, P4_1024
	PRCENTRY2		xfft_r4delay_256K_10_4, 1305600, P4TP_256
	PRCENTRY2A		xfft_r4delay_256K_10_4_P4, 1305600, I7
	PRCENTRY2		xfft_r4delay_256K_10_2, 1305600, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_256K_10_1, 1305600, K10
	PRCENTRY2		xfft_r4delay_256K_8_4, 878592, CORE2
	PRCENTRY2A		xfft_r4delay_256K_8_4_CORE, 878592, P4_1024
	PRCENTRY2		xfft_hg_256K_10_4, 721920, P4TP_512
	PRCENTRY2		xfft_hg_256K_10_2, 721920
	DD			0
	PRCSTRT	5790000, 294912, 0.0172
	PRCENTRY2		xfft_r4dwpn_288K_768_4, 1076736, P4_1024 + I7_64 + CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_288K_768_2, 1076736, P4TP_256 + K8
	PRCENTRY2A		xfft_r4dwpn_288K_768_2_P4, 1076736, I7_32
	PRCENTRY2		xfft_r4dwpn_288K_768_1, 1076736
	PRCENTRY2		xfft_r4delay_288K_768_4, 1012224, P4_1024 + I7_64 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_288K_768_2, 1012224, P4TP_256 + K8
	PRCENTRY2A		xfft_r4delay_288K_768_2_P4, 1012224, I7_32
	PRCENTRY2		xfft_r4delay_288K_768_1, 1012224
	DD			0
	PRCSTRT	6435000, 327680, 0.0192
	PRCENTRY2		xfft_r4dwpn_320K_1280_4, 2009088
	PRCENTRY2		xfft_r4dwpn_320K_1280_2, 2009088
	PRCENTRY2		xfft_r4dwpn_320K_1280_1, 2009088, K8_64
	PRCENTRY2A		xfft_r4dwpn_320K_1280_4_P4, 2009088, I7_64
	PRCENTRY2		xfft_r4dwpn_320K_10_4, 1280000, P4_1024_64 + CORE2
	PRCENTRY2A		xfft_r4dwpn_320K_10_4_P4, 1280000, I7_32
	PRCENTRY2		xfft_r4dwpn_320K_10_2, 1280000, P4TP_256 + CORE2_512
	PRCENTRY2		xfft_r4dwpn_320K_10_1, 1280000, K8_32 + K10
	PRCENTRY2		xfft_r4delay_320K_1280_4, 1623040
	PRCENTRY2A		xfft_r4delay_320K_1280_4_P4, 1623040, I7_64
	PRCENTRY2		xfft_r4delay_320K_1280_2, 1623040
	PRCENTRY2		xfft_r4delay_320K_1280_1, 1623040, K8_64
	PRCENTRY2		xfft_r4delay_320K_10_4, 1179136, P4_1024_64 + CORE2
	PRCENTRY2A		xfft_r4delay_320K_10_4_P4, 1179136, I7_32
	PRCENTRY2		xfft_r4delay_320K_10_2, 1179136, P4TP_256
	PRCENTRY2		xfft_r4delay_320K_10_1, 1179136, CORE2_512 + K10
	PRCENTRY2		xfft_hg_320K_10_4, 890624, P4_1024_32 + P4TP_512
	PRCENTRY2		xfft_hg_320K_10_2, 890624
	PRCENTRY2		xfft_hg_320K_10_1, 890624, K8_32
	DD			0
	PRCSTRT	6749000, 344064, 0.0205
	PRCENTRY2		xfft_r4dwpn_336K_768_4, 1192960, P4_1024_64 + I7_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4dwpn_336K_768_4_P4, 1192960, CORE2_32 + I7_32
	PRCENTRY2A		xfft_r4dwpn_336K_768_4_CORE, 1192960, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_336K_768_2, 1192960, P4TP_512 + P4TP_256 + K8
	PRCENTRY2		xfft_r4dwpn_336K_768_1, 1192960
	PRCENTRY2		xfft_r4delay_336K_768_4, 1140224, P4_1024_64 + I7_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4delay_336K_768_4_P4, 1140224, CORE2_32 + I7_32
	PRCENTRY2A		xfft_r4delay_336K_768_4_CORE, 1140224, P4_1024_32
	PRCENTRY2		xfft_r4delay_336K_768_2, 1140224, P4TP_512 + P4TP_256 + K8
	PRCENTRY2		xfft_r4delay_336K_768_1, 1140224
	DD			0
	PRCSTRT	7692000, 393216, 0.0238
	PRCENTRY2		xfft_r4dwpn_384K_1536_4, 2408448, P4_1024_64 + P4TP_256 + I7_64
	PRCENTRY2A		xfft_r4dwpn_384K_1536_4_P4, 2408448, I7_32
	PRCENTRY2A		xfft_r4dwpn_384K_1536_4_CORE, 2408448, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_384K_10_4, 1431040, K10
	PRCENTRY2		xfft_r4dwpn_384K_10_2, 1431040, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_384K_10_1, 1431040, K8_32
	PRCENTRY2A		xfft_r4dwpn_384K_10_1_P4, 1431040, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_384K_768_4, 1308672, CORE2
	PRCENTRY2		xfft_r4dwpn_384K_768_2, 1308672
	PRCENTRY2		xfft_r4dwpn_384K_768_1, 1308672
	PRCENTRY2		xfft_r4delay_384K_1536_4, 1944576, P4_1024_64 + P4TP_256 + I7_64
	PRCENTRY2A		xfft_r4delay_384K_1536_4_P4, 1944576, I7_32
	PRCENTRY2		xfft_r4delay_384K_10_4, 1346048, CORE2_64
	PRCENTRY2		xfft_r4delay_384K_10_2, 1346048, K10_32
	PRCENTRY2		xfft_r4delay_384K_10_1, 1346048, CORE2_512_64 + K8 + K10_64
	PRCENTRY2A		xfft_r4delay_384K_10_1_P4, 1346048, CORE2_512_32
	PRCENTRY2		xfft_r4delay_384K_768_4, 1267712
	PRCENTRY2A		xfft_r4delay_384K_768_4_P4, 1267712, CORE2_32
	PRCENTRY2		xfft_r4delay_384K_768_2, 1267712
	PRCENTRY2		xfft_r4delay_384K_768_1, 1267712
	PRCENTRY2		xfft_hg_384K_10_4, 1058432, P4_1024_32 + P4TP_512
	PRCENTRY2		xfft_hg_384K_10_2, 1058432
	PRCENTRY2		xfft_hg_384K_10_1, 1058432
	DD			0
	PRCSTRT	8012000, 409600, 0.0248
	PRCENTRY2		xfft_r4dwpn_400K_1280_4, 1589248, P4_1024_32 + P4TP_512 + CORE2
	PRCENTRY2A		xfft_r4dwpn_400K_1280_4_P4, 1589248, I7
	PRCENTRY2A		xfft_r4dwpn_400K_1280_4_CORE, 1589248, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_400K_1280_2, 1589248, P4TP_256 + CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_400K_1280_1, 1589248, K8 + K10
	PRCENTRY2		xfft_r4delay_400K_1280_4, 1463808, P4_1024_32 + P4TP_512 + CORE2_64
	PRCENTRY2A		xfft_r4delay_400K_1280_4_P4, 1463808, I7
	PRCENTRY2A		xfft_r4delay_400K_1280_4_CORE, 1463808, P4_1024_64
	PRCENTRY2		xfft_r4delay_400K_1280_2, 1463808, P4TP_256 + CORE2_32 + CORE2_512_64
	PRCENTRY2		xfft_r4delay_400K_1280_1, 1463808, CORE2_512_32 + K8 + K10
	DD			0
	PRCSTRT	8958000, 458752, 0.0283
	PRCENTRY2		xfft_r4dwpn_448K_10_4, 1584128, P4_1024_64 + CORE2 + K10
	PRCENTRY2A		xfft_r4dwpn_448K_10_4_P4, 1584128, I7_64
	PRCENTRY2A		xfft_r4dwpn_448K_10_4_CORE, 1584128, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_448K_10_2, 1584128, K8
	PRCENTRY2A		xfft_r4dwpn_448K_10_2_P4, 1584128, I7_32
	PRCENTRY2		xfft_r4dwpn_448K_10_1, 1584128, CORE2_512
	PRCENTRY2		xfft_r4delay_448K_10_4, 1515008, P4_1024_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4delay_448K_10_4_P4, 1515008, CORE2_32 + I7_32
	PRCENTRY2A		xfft_r4delay_448K_10_4_CORE, 1515008, P4_1024_32
	PRCENTRY2		xfft_r4delay_448K_10_2, 1515008, P4TP_256 + CORE2_512_32 + I7_64 + K8_64
	PRCENTRY2		xfft_r4delay_448K_10_1, 1515008, CORE2_512_64 + K8_32
	PRCENTRY2		xfft_hg_448K_10_4, 1226496, P4TP_512
	PRCENTRY2		xfft_hg_448K_10_2, 1226496
	PRCENTRY2		xfft_hg_448K_10_1, 1226496
	DD			0
	PRCSTRT	9567000, 491520, 0.0298
	PRCENTRY2		xfft_r4dwpn_480K_1536_4, 1902592, I7_64
	PRCENTRY2A		xfft_r4dwpn_480K_1536_4_CORE, 1902592, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_480K_1536_2, 1902592, P4TP_512 + P4TP_256 + CORE2_512_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_480K_1280_4, 1777152, K10_64
	PRCENTRY2		xfft_r4dwpn_480K_1280_2, 1777152, K8_64
	PRCENTRY2		xfft_r4dwpn_480K_1280_1, 1777152, CORE2_512_32 + K8_32
	PRCENTRY2		xfft_r4dwpn_480K_768_4, 1893888, P4_1024_64 + CORE2
	PRCENTRY2A		xfft_r4dwpn_480K_768_4_P4, 1893888, I7_32
	PRCENTRY2		xfft_r4delay_480K_1536_4, 1752576, I7_64
	PRCENTRY2A		xfft_r4delay_480K_1536_4_CORE, 1752576, P4_1024
	PRCENTRY2		xfft_r4delay_480K_1536_2, 1752576, P4TP_512 + P4TP_256 + K10_32
	PRCENTRY2		xfft_r4delay_480K_1280_4, 1671680
	PRCENTRY2		xfft_r4delay_480K_1280_2, 1671680
	PRCENTRY2		xfft_r4delay_480K_1280_1, 1671680, CORE2_512 + K8 + K10_64
	PRCENTRY2		xfft_r4delay_480K_768_4, 1741824, CORE2
	PRCENTRY2A		xfft_r4delay_480K_768_4_P4, 1741824, I7_32
	DD			0
	PRCSTRT	10180000, 524288, 0.0319
	PRCENTRY2		xfft_r4dwpn_512K_11_4, 3215360, P4TP_256
	PRCENTRY2A		xfft_r4dwpn_512K_11_4_P4, 3215360, I7
	PRCENTRY2A		xfft_r4dwpn_512K_11_4_CORE, 3215360, P4_1024
	PRCENTRY2		xfft_r4dwpn_512K_10_4, 1736704, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_512K_10_2, 1736704, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_512K_10_1, 1736704
	PRCENTRY2A		xfft_r4dwpn_512K_10_1_P4, 1736704, CORE2_512_32
	PRCENTRY2		xfft_r4delay_512K_11_4, 2595840, P4TP_256
	PRCENTRY2A		xfft_r4delay_512K_11_4_P4, 2595840, I7
	PRCENTRY2A		xfft_r4delay_512K_11_4_CORE, 2595840, P4_1024_64
	PRCENTRY2		xfft_r4delay_512K_10_4, 1683456, CORE2 + K10
	PRCENTRY2		xfft_r4delay_512K_10_2, 1683456, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_512K_10_1, 1683456
	PRCENTRY2		xfft_hg_512K_10_4, 1392640, P4_1024_32 + P4TP_512
	PRCENTRY2		xfft_hg_512K_10_2, 1392640
	PRCENTRY2		xfft_hg_512K_10_1, 1392640
	DD			0
	PRCSTRT	11130000, 573440, 0.0353
	PRCENTRY2		xfft_r4dwpn_560K_1280_4, 1967104, P4TP_512 + CORE2 + K10_64
	PRCENTRY2A		xfft_r4dwpn_560K_1280_4_P4, 1967104, I7_64
	PRCENTRY2A		xfft_r4dwpn_560K_1280_4_CORE, 1967104, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_560K_1280_2, 1967104, P4TP_256
	PRCENTRY2A		xfft_r4dwpn_560K_1280_2_P4, 1967104, I7_32
	PRCENTRY2		xfft_r4delay_560K_1280_4, 1881600, P4_1024_32 + P4TP_512 + CORE2 + K10_64
	PRCENTRY2A		xfft_r4delay_560K_1280_4_P4, 1881600, I7_64
	PRCENTRY2A		xfft_r4delay_560K_1280_4_CORE, 1881600, P4_1024_64
	PRCENTRY2		xfft_r4delay_560K_1280_2, 1881600, P4TP_256 + I7_32
	DD			0
	PRCSTRT	11430000, 589824, 0.0365
	PRCENTRY2		xfft_r4dwpn_576K_2304_4, 3594240, I7_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_576K_2304_2, 3594240, K8_32
	PRCENTRY2		xfft_r4dwpn_576K_2304_1, 3594240
	PRCENTRY2		xfft_r4dwpn_576K_1536_4, 2127360, P4_1024_64 + P4TP_512 + K10_64
	PRCENTRY2A		xfft_r4dwpn_576K_1536_4_CORE, 2127360, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_576K_1536_2, 2127360, P4TP_256 + CORE2_512 + K10_32
	PRCENTRY2		xfft_r4dwpn_576K_768_4, 2119680, CORE2
	PRCENTRY2A		xfft_r4dwpn_576K_768_4_P4, 2119680, I7_32
	PRCENTRY2		xfft_r4dwpn_576K_768_2, 2119680
	PRCENTRY2		xfft_r4dwpn_576K_768_1, 2119680
	PRCENTRY2		xfft_r4delay_576K_2304_4, 2896896, P4_1024_64 + I7_64 + K8_64
	PRCENTRY2		xfft_r4delay_576K_2304_2, 2896896, K8_32
	PRCENTRY2		xfft_r4delay_576K_2304_1, 2896896
	PRCENTRY2		xfft_r4delay_576K_1536_4, 2001408, P4TP_512 + K10_32
	PRCENTRY2A		xfft_r4delay_576K_1536_4_CORE, 2001408, P4_1024_32
	PRCENTRY2		xfft_r4delay_576K_1536_2, 2001408, P4TP_256 + CORE2_512 + K10_64
	PRCENTRY2		xfft_r4delay_576K_768_4, 1993728, CORE2
	PRCENTRY2A		xfft_r4delay_576K_768_4_P4, 1993728, I7_32
	PRCENTRY2		xfft_r4delay_576K_768_2, 1993728
	PRCENTRY2		xfft_r4delay_576K_768_1, 1993728
	DD			0
	PRCSTRT	12680000, 655360, 0.0410
	PRCENTRY2		xfft_r4dwpn_640K_2560_4, 4005888
	PRCENTRY2		xfft_r4dwpn_640K_11_4, 2537472, P4_1024_32 + P4TP_512
	PRCENTRY2		xfft_r4dwpn_640K_11_2, 2537472, P4TP_256
	PRCENTRY2		xfft_r4dwpn_640K_1280_4, 2156544
	PRCENTRY2		xfft_r4dwpn_640K_10_4, 2518528, P4_1024_64 + CORE2 + K10_64
	PRCENTRY2A		xfft_r4dwpn_640K_10_4_P4, 2518528, I7
	PRCENTRY2		xfft_r4dwpn_640K_10_2, 2518528
	PRCENTRY2		xfft_r4dwpn_640K_10_1, 2518528, CORE2_512 + K8 + K10_32
	PRCENTRY2		xfft_r4delay_640K_2560_4, 3230720
	PRCENTRY2		xfft_r4delay_640K_11_4, 2338304, P4_1024_32
	PRCENTRY2		xfft_r4delay_640K_11_2, 2338304, P4TP_256 + CORE2_512_64 + K10_32
	PRCENTRY2		xfft_r4delay_640K_1280_4, 2091008, K8_64
	PRCENTRY2		xfft_r4delay_640K_10_4, 2313216, P4_1024_64 + CORE2 + K10_64
	PRCENTRY2A		xfft_r4delay_640K_10_4_P4, 2313216, I7
	PRCENTRY2		xfft_r4delay_640K_10_2, 2313216
	PRCENTRY2		xfft_r4delay_640K_10_1, 2313216, CORE2_512_32
	PRCENTRY2		xfft_hg_640K_12_2, 1214848
	PRCENTRY2		xfft_hg_640K_11_4, 1437440, P4TP_512
	PRCENTRY2		xfft_hg_640K_11_1, 1437440, K8_32
	DD			0
	PRCSTRT	13300000, 688128, 0.0434
	PRCENTRY2		xfft_r4dwpn_672K_1536_4, 2354176, I7_64 + K10
	PRCENTRY2A		xfft_r4dwpn_672K_1536_4_CORE, 2354176, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_672K_1536_2, 2354176, P4TP_512 + P4TP_256 + K8_32
	PRCENTRY2		xfft_r4dwpn_672K_1536_1, 2354176, CORE2_512
	PRCENTRY2		xfft_r4dwpn_672K_768_4, 2349568, P4_1024_64 + CORE2
	PRCENTRY2A		xfft_r4dwpn_672K_768_4_P4, 2349568, I7_32
	PRCENTRY2		xfft_r4dwpn_672K_768_2, 2349568, K8_64
	PRCENTRY2		xfft_r4dwpn_672K_768_1, 2349568
	PRCENTRY2		xfft_r4delay_672K_1536_4, 2252288, I7_64 + K10_32
	PRCENTRY2A		xfft_r4delay_672K_1536_4_CORE, 2252288, P4_1024_32
	PRCENTRY2		xfft_r4delay_672K_1536_2, 2252288, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_r4delay_672K_1536_1, 2252288, CORE2_512 + K8 + K10_64
	PRCENTRY2		xfft_r4delay_672K_768_4, 2249728, P4_1024_64 + CORE2
	PRCENTRY2A		xfft_r4delay_672K_768_4_P4, 2249728, I7_32
	PRCENTRY2		xfft_r4delay_672K_768_2, 2249728
	PRCENTRY2		xfft_r4delay_672K_768_1, 2249728
	DD			0
	PRCSTRT	14180000, 737280, 0.0471
	PRCENTRY2		xfft_r4dwpn_720K_2304_4, 2830336, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_720K_2304_2, 2830336, P4TP_256
	PRCENTRY2		xfft_r4dwpn_720K_2304_1, 2830336
	PRCENTRY2		xfft_r4delay_720K_2304_4, 2606592, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_720K_2304_2, 2606592, P4TP_256
	PRCENTRY2		xfft_r4delay_720K_2304_1, 2606592
	DD			0
	PRCSTRT	15160000, 786432, 0.0507
	PRCENTRY2		xfft_r4dwpn_768K_3072_4, 4808704, I7_64
	PRCENTRY2		xfft_r4dwpn_768K_3072_2, 4808704, P4TP_256 + I7_32 + CORE2_512
	PRCENTRY2		xfft_r4dwpn_768K_3072_1, 4808704, K8
	PRCENTRY2		xfft_r4dwpn_768K_11_4, 2835968
	PRCENTRY2		xfft_r4dwpn_768K_1536_4, 2580480, K10_32
	PRCENTRY2		xfft_r4dwpn_768K_10_4, 2818048, K10_64
	PRCENTRY2		xfft_r4dwpn_768K_10_2, 2818048
	PRCENTRY2		xfft_r4dwpn_768K_10_1, 2818048
	PRCENTRY2		xfft_r4dwpn_768K_768_4, 2578432, CORE2
	PRCENTRY2		xfft_r4dwpn_768K_768_2, 2578432
	PRCENTRY2		xfft_r4dwpn_768K_768_1, 2578432
	PRCENTRY2		xfft_r4delay_768K_3072_4, 3877888, I7_64
	PRCENTRY2		xfft_r4delay_768K_3072_2, 3877888, CORE2_512
	PRCENTRY2		xfft_r4delay_768K_3072_1, 3877888, K8
	PRCENTRY2		xfft_r4delay_768K_11_4, 2669056
	PRCENTRY2		xfft_r4delay_768K_1536_4, 2502656, K10
	PRCENTRY2		xfft_r4delay_768K_10_4, 2647040, CORE2_32
	PRCENTRY2A		xfft_r4delay_768K_10_4_P4, 2647040, I7_32
	PRCENTRY2		xfft_r4delay_768K_10_2, 2647040
	PRCENTRY2		xfft_r4delay_768K_10_1, 2647040
	PRCENTRY2		xfft_r4delay_768K_768_4, 2504704, CORE2_64
	PRCENTRY2		xfft_r4delay_768K_768_2, 2504704
	PRCENTRY2		xfft_r4delay_768K_768_1, 2504704
	PRCENTRY2		xfft_hg_768K_12_2, 1413504, P4_1024
	PRCENTRY2		xfft_hg_768K_11_4, 1703552, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_hg_768K_11_2, 1703552
	DD			0
	PRCSTRT	15770000, 819200, 0.0528
	PRCENTRY2		xfft_r4dwpn_800K_2560_4, 3155968
	PRCENTRY2A		xfft_r4dwpn_800K_2560_4_CORE, 3155968, P4_1024
	PRCENTRY2		xfft_r4dwpn_800K_2560_2, 3155968, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_r4dwpn_800K_1280_4, 3134976, I7_32 + CORE2 + K10
	PRCENTRY2A		xfft_r4dwpn_800K_1280_4_P4, 3134976, I7_64
	PRCENTRY2		xfft_r4dwpn_800K_1280_2, 3134976
	PRCENTRY2		xfft_r4dwpn_800K_1280_1, 3134976, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_800K_2560_4, 2907648, P4TP_512 + K10_32
	PRCENTRY2A		xfft_r4delay_800K_2560_4_CORE, 2907648, P4_1024
	PRCENTRY2		xfft_r4delay_800K_2560_2, 2907648, P4TP_256
	PRCENTRY2		xfft_r4delay_800K_1280_4, 2876416, CORE2 + K10_64
	PRCENTRY2A		xfft_r4delay_800K_1280_4_P4, 2876416, I7
	PRCENTRY2		xfft_r4delay_800K_1280_2, 2876416
	PRCENTRY2		xfft_r4delay_800K_1280_1, 2876416, CORE2_512 + K8
	DD			0
	PRCSTRT	16940000, 884736, 0.0582
	PRCENTRY2		xfft_r4dwpn_864K_2304_4, 3165696, K10
	PRCENTRY2		xfft_r4dwpn_864K_2304_2, 3165696, K8
	PRCENTRY2		xfft_r4dwpn_864K_2304_1, 3165696
	PRCENTRY2		xfft_r4delay_864K_2304_4, 2978304, K8 + K10
	PRCENTRY2		xfft_r4delay_864K_2304_2, 2978304
	PRCENTRY2		xfft_r4delay_864K_2304_1, 2978304
	DD			0
	PRCSTRT	17640000, 917504, 0.0607
	PRCENTRY2		xfft_r4dwpn_896K_11_4, 3136512
	PRCENTRY2		xfft_r4dwpn_896K_11_2, 3136512, P4TP_256 + CORE2_512
	PRCENTRY2		xfft_r4dwpn_896K_10_4, 3121664, P4_1024_64 + CORE2 + K10
	PRCENTRY2A		xfft_r4dwpn_896K_10_4_P4, 3121664, I7
	PRCENTRY2A		xfft_r4dwpn_896K_10_4_CORE, 3121664, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_896K_10_2, 3121664
	PRCENTRY2		xfft_r4dwpn_896K_10_1, 3121664, K8
	PRCENTRY2		xfft_r4delay_896K_11_4, 3001856, K10_32
	PRCENTRY2A		xfft_r4delay_896K_11_4_CORE, 3001856, P4_1024_32
	PRCENTRY2		xfft_r4delay_896K_11_2, 3001856, P4TP_256 + CORE2_512
	PRCENTRY2		xfft_r4delay_896K_10_4, 2984960, P4_1024_64 + CORE2 + K10_64
	PRCENTRY2A		xfft_r4delay_896K_10_4_P4, 2984960, I7
	PRCENTRY2		xfft_r4delay_896K_10_2, 2984960
	PRCENTRY2		xfft_r4delay_896K_10_1, 2984960, K8_64
	PRCENTRY2		xfft_hg_896K_12_2, 1612416, K8_32
	PRCENTRY2		xfft_hg_896K_12_1, 1612416
	PRCENTRY2		xfft_hg_896K_11_4, 1969920, P4TP_512
	PRCENTRY2		xfft_hg_896K_11_2, 1969920
	DD			0
	PRCSTRT	18800000, 983040, 0.0642
	PRCENTRY2		xfft_r4dwpn_960K_3840_4, 5994496, I7_64
	PRCENTRY2		xfft_r4dwpn_960K_3840_2, 5994496, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_960K_3840_1, 5994496, K8_32
	PRCENTRY2		xfft_r4dwpn_960K_3072_4, 3786752
	PRCENTRY2		xfft_r4dwpn_960K_3072_2, 3786752, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_r4dwpn_960K_2560_4, 3528192
	PRCENTRY2		xfft_r4dwpn_960K_1536_4, 3755520, P4_1024_64 + I7_32 + K10
	PRCENTRY2A		xfft_r4dwpn_960K_1536_4_CORE, 3755520, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_960K_1280_4, 3508224, CORE2_32
	PRCENTRY2		xfft_r4dwpn_960K_768_4, 2000896, CORE2_64
	PRCENTRY2		xfft_r4delay_960K_3840_4, 4830208, I7_64
	PRCENTRY2		xfft_r4delay_960K_3840_2, 4830208, CORE2_512
	PRCENTRY2		xfft_r4delay_960K_3840_1, 4830208, K8
	PRCENTRY2		xfft_r4delay_960K_3072_4, 3489280
	PRCENTRY2		xfft_r4delay_960K_3072_2, 3489280, P4TP_512 + P4TP_256
	PRCENTRY2		xfft_r4delay_960K_2560_4, 3320320
	PRCENTRY2		xfft_r4delay_960K_1536_4, 3443712, K10
	PRCENTRY2A		xfft_r4delay_960K_1536_4_P4, 3443712, I7_32
	PRCENTRY2A		xfft_r4delay_960K_1536_4_CORE, 3443712, P4_1024
	PRCENTRY2		xfft_r4delay_960K_1280_4, 3292160
	PRCENTRY2		xfft_r4delay_960K_768_4, 2101248, CORE2
	DD			0
	PRCSTRT	19740000, 1032192, 0.0667
	PRCENTRY2		xfft_r4dwpn_1008K_2304_4, 3503104
	PRCENTRY2		xfft_r4delay_1008K_2304_4, 3352064
	DD			0
	PRCSTRT	20090000, 1048576, 0.0676
	PRCENTRY2		xfft_r4dwpn_1M_12_4, 6406144, P4_1024 + I7 + K10_64
	PRCENTRY2		xfft_r4dwpn_1M_12_2, 6406144, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_1M_12_1, 6406144, K8_32
	PRCENTRY2		xfft_r4dwpn_1M_11_4, 3436544, K10_32
	PRCENTRY2		xfft_r4dwpn_1M_11_2, 3436544, P4TP_256
	PRCENTRY2		xfft_r4dwpn_1M_10_4, 3424256, CORE2
	PRCENTRY2		xfft_r4delay_1M_12_4, 5164032, P4_1024_64 + I7 + K10_64
	PRCENTRY2		xfft_r4delay_1M_12_2, 5164032, CORE2_512
	PRCENTRY2		xfft_r4delay_1M_12_1, 5164032, K8
	PRCENTRY2		xfft_r4delay_1M_11_4, 3334144, K10_32
	PRCENTRY2		xfft_r4delay_1M_11_2, 3334144, P4TP_256
	PRCENTRY2		xfft_r4delay_1M_10_4, 3321856, CORE2
	PRCENTRY2		xfft_hg_1024K_12_4, 1809408
	PRCENTRY2		xfft_hg_1024K_12_2, 1809408
	PRCENTRY2		xfft_hg_1024K_12_1, 1809408
	PRCENTRY2		xfft_hg_1024K_11_4, 2234368, P4_1024_32 + P4TP_512
	PRCENTRY2		xfft_hg_1024K_11_2, 2234368
	DD			0
	PRCSTRT	21930000, 1146880, 0.0761
	PRCENTRY2		xfft_r4dwpn_1120K_2560_4, 3902464, P4_1024_64 + I7_64 + CORE2_64
	PRCENTRY2A		xfft_r4dwpn_1120K_2560_4_CORE, 3902464, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_1120K_2560_2, 3902464, P4TP_512 + P4TP_256 + I7_32
	PRCENTRY2		xfft_r4dwpn_1120K_2560_1, 3902464
	PRCENTRY2		xfft_r4dwpn_1120K_1280_4, 3885568, K10_64
	PRCENTRY2		xfft_r4delay_1120K_2560_4, 3735040, P4TP_512 + CORE2_64
	PRCENTRY2A		xfft_r4delay_1120K_2560_4_CORE, 3735040, P4_1024
	PRCENTRY2		xfft_r4delay_1120K_2560_2, 3735040, P4TP_256 + I7
	PRCENTRY2		xfft_r4delay_1120K_2560_1, 3735040
	PRCENTRY2		xfft_r4delay_1120K_1280_4, 3712000, K10_64
	DD			0
	PRCSTRT	22510000, 1179648, 0.0784
	PRCENTRY2		xfft_r4dwpn_1152K_4608_4, 7176192, P4_1024 + I7_64
	PRCENTRY2		xfft_r4dwpn_1152K_4608_2, 7176192, P4TP_512 + I7_32
	PRCENTRY2		xfft_r4dwpn_1152K_3072_4, 4232704, CORE2_64 + K10
	PRCENTRY2		xfft_r4dwpn_1152K_3072_2, 4232704, P4TP_256 + CORE2_32 + CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_1152K_3072_1, 4232704, CORE2_512_32 + K8
	PRCENTRY2		xfft_r4dwpn_1152K_2304_4, 3840000
	PRCENTRY2		xfft_r4dwpn_1152K_2304_2, 3840000
	PRCENTRY2		xfft_r4dwpn_1152K_2304_1, 3840000
	PRCENTRY2		xfft_r4delay_1152K_4608_4, 5778432, P4_1024 + P4TP_512 + I7_64 + K10_32
	PRCENTRY2		xfft_r4delay_1152K_4608_2, 5778432, I7_32
	PRCENTRY2		xfft_r4delay_1152K_3072_4, 3983872, CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_1152K_3072_2, 3983872, P4TP_256 + CORE2_32 + CORE2_512_32
	PRCENTRY2		xfft_r4delay_1152K_3072_1, 3983872, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4delay_1152K_2304_4, 3725312
	PRCENTRY2		xfft_r4delay_1152K_2304_2, 3725312
	PRCENTRY2		xfft_r4delay_1152K_2304_1, 3725312
	DD			0
	PRCSTRT	23440000, 1228800, 0.0832
	PRCENTRY2		xfft_r4dwpn_1200K_3840_4, 4714496, P4_1024_64 + I7_64 + CORE2 + CORE2_512_32 + K10
	PRCENTRY2		xfft_r4dwpn_1200K_3840_2, 4714496, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_1200K_3840_1, 4714496, K8_32
	PRCENTRY2		xfft_r4delay_1200K_3840_4, 4343296, P4_1024_64 + I7_64 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_1200K_3840_2, 4343296, CORE2_512_64
	PRCENTRY2		xfft_r4delay_1200K_3840_1, 4343296, CORE2_512_32 + K8
	DD			0
	PRCSTRT	24980000, 1310720, 0.0892
	PRCENTRY2		xfft_r4dwpn_1280K_5120_4, 8003584
	PRCENTRY2		xfft_r4dwpn_1280K_5120_2, 8003584, K8_64
	PRCENTRY2		xfft_r4dwpn_1280K_5120_1, 8003584
	PRCENTRY2		xfft_r4dwpn_1280K_12_4, 5040128, CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_1280K_12_2, 5040128, P4TP_512 + CORE2_512
	PRCENTRY2		xfft_r4dwpn_1280K_12_1, 5040128, K8_32
	PRCENTRY2		xfft_r4dwpn_1280K_2560_4, 4276224
	PRCENTRY2		xfft_r4dwpn_1280K_2560_2, 4276224
	PRCENTRY2		xfft_r4dwpn_1280K_11_4, 5004800, P4_1024_64 + I7_32 + K10_32
	PRCENTRY2A		xfft_r4dwpn_1280K_11_4_CORE, 5004800, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_1280K_11_2, 5004800, P4TP_256 + I7_64
	PRCENTRY2		xfft_r4dwpn_1280K_1280_4, 4261888, CORE2_32
	PRCENTRY2		xfft_r4dwpn_1280K_10_4, 2641920
	PRCENTRY2		xfft_r4delay_1280K_5120_4, 6450176, I7_64
	PRCENTRY2		xfft_r4delay_1280K_5120_2, 6450176
	PRCENTRY2		xfft_r4delay_1280K_5120_1, 6450176, K8_64
	PRCENTRY2		xfft_r4delay_1280K_12_4, 4644352, CORE2 + K10
	PRCENTRY2		xfft_r4delay_1280K_12_2, 4644352, CORE2_512_64
	PRCENTRY2		xfft_r4delay_1280K_12_1, 4644352, CORE2_512_32
	PRCENTRY2		xfft_r4delay_1280K_2560_4, 4149248
	PRCENTRY2		xfft_r4delay_1280K_2560_2, 4149248
	PRCENTRY2		xfft_r4delay_1280K_11_4, 4586496, P4_1024_64 + I7_32
	PRCENTRY2		xfft_r4delay_1280K_11_2, 4586496, P4TP_256
	PRCENTRY2		xfft_r4delay_1280K_1280_4, 4130816
	PRCENTRY2		xfft_r4delay_1280K_10_4, 2779136
	PRCENTRY2		xfft_hg_1280K_12_4, 2207488, P4_1024_32 + P4TP_512
	PRCENTRY2		xfft_hg_1280K_12_1, 2207488, K8_32
	PRCENTRY2		xfft_hg_1280K_11_2, 2769152
	DD			0
	PRCSTRT	26200000, 1376256, 0.0952
	PRCENTRY2		xfft_r4dwpn_1344K_3072_4, 4680704, CORE2_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_1344K_3072_2, 4680704, P4TP_512 + P4TP_256 + I7_32
	PRCENTRY2		xfft_r4dwpn_1344K_3072_1, 4680704, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_1344K_1536_4, 4653568, P4_1024_64 + I7_64 + CORE2_32 + K10_64
	PRCENTRY2A		xfft_r4dwpn_1344K_1536_4_CORE, 4653568, P4_1024_32
	PRCENTRY2		xfft_r4delay_1344K_3072_4, 4480512, CORE2 + K10_32
	PRCENTRY2		xfft_r4delay_1344K_3072_2, 4480512, P4TP_512 + P4TP_256 + I7 + CORE2_512
	PRCENTRY2		xfft_r4delay_1344K_3072_1, 4480512, K8 + K10_64
	PRCENTRY2		xfft_r4delay_1344K_1536_4, 4443136, P4_1024_32
	PRCENTRY2A		xfft_r4delay_1344K_1536_4_CORE, 4443136, P4_1024_64
	DD			0
	PRCSTRT	27990000, 1474560, 0.104
	PRCENTRY2		xfft_r4dwpn_1440K_4608_4, 5638144, P4_1024 + P4TP_512 + I7_64
	PRCENTRY2		xfft_r4dwpn_1440K_4608_2, 5638144, P4TP_256 + I7_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_1440K_3840_4, 5271040, K10_64
	PRCENTRY2		xfft_r4dwpn_1440K_3840_2, 5271040, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_1440K_3840_1, 5271040, K8_32
	PRCENTRY2		xfft_r4dwpn_1440K_2304_4, 5604864, CORE2
	PRCENTRY2		xfft_r4delay_1440K_4608_4, 5193216, P4_1024 + I7_64 + K10_32
	PRCENTRY2		xfft_r4delay_1440K_4608_2, 5193216, P4TP_512 + P4TP_256 + I7_32
	PRCENTRY2		xfft_r4delay_1440K_3840_4, 4960768, K10_64
	PRCENTRY2		xfft_r4delay_1440K_3840_2, 4960768
	PRCENTRY2		xfft_r4delay_1440K_3840_1, 4960768, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4delay_1440K_2304_4, 5133312, CORE2
	DD			0
	PRCSTRT	29870000, 1572864, 0.113
	PRCENTRY2		xfft_r4dwpn_1536K_6144_4, 9605120
	PRCENTRY2		xfft_r4dwpn_1536K_12_4, 5633536, CORE2_32 + K10_64
	PRCENTRY2		xfft_r4dwpn_1536K_12_2, 5633536, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_1536K_12_1, 5633536, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_1536K_3072_4, 5128192, I7_64 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_1536K_3072_2, 5128192, P4TP_256
	PRCENTRY2		xfft_r4dwpn_1536K_11_4, 5599232, I7_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_1536K_11_2, 5599232
	PRCENTRY2		xfft_r4dwpn_1536K_1536_4, 5103616
	PRCENTRY2		xfft_r4dwpn_1536K_10_4, 2596864
	PRCENTRY2		xfft_r4delay_1536K_6144_4, 7740416
	PRCENTRY2		xfft_r4delay_1536K_12_4, 5302784, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_1536K_12_2, 5302784, CORE2_512_32 + K8_64
	PRCENTRY2		xfft_r4delay_1536K_12_1, 5302784, CORE2_512_64 + K8_32
	PRCENTRY2		xfft_r4delay_1536K_3072_4, 4976640, CORE2_64
	PRCENTRY2		xfft_r4delay_1536K_3072_2, 4976640, P4TP_256 + I7_64
	PRCENTRY2		xfft_r4delay_1536K_11_4, 5248000, I7_32
	PRCENTRY2		xfft_r4delay_1536K_11_2, 5248000
	PRCENTRY2		xfft_r4delay_1536K_1536_4, 4943872
	PRCENTRY2		xfft_r4delay_1536K_10_4, 2797568
	PRCENTRY2		xfft_hg_1536K_12_4, 2604672, P4_1024_32 + P4TP_512
	PRCENTRY2A		xfft_hg_1536K_12_4_CORE, 2604672, P4_1024_64
	PRCENTRY2		xfft_hg_1536K_12_1, 2604672
	PRCENTRY2		xfft_hg_1536K_11_1, 3301760
	DD			0
	PRCSTRT	31060000, 1638400, 0.119
	PRCENTRY2		xfft_r4dwpn_1600K_5120_4, 6293504, P4TP_512 + I7_64 + CORE2 + CORE2_512_32 + K10
	PRCENTRY2		xfft_r4dwpn_1600K_5120_2, 6293504, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_1600K_5120_1, 6293504, K8_32
	PRCENTRY2		xfft_r4dwpn_1600K_2560_4, 6237696, P4_1024_64 + I7_32
	PRCENTRY2A		xfft_r4dwpn_1600K_2560_4_CORE, 6237696, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_1600K_1280_4, 3274752
	PRCENTRY2		xfft_r4delay_1600K_5120_4, 5799424, P4TP_512 + I7_64 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_1600K_5120_2, 5799424, CORE2_512_64
	PRCENTRY2		xfft_r4delay_1600K_5120_1, 5799424, CORE2_512_32 + K8
	PRCENTRY2		xfft_r4delay_1600K_2560_4, 5712896, I7_32
	PRCENTRY2A		xfft_r4delay_1600K_2560_4_CORE, 5712896, P4_1024
	PRCENTRY2		xfft_r4delay_1600K_1280_4, 3448832
	DD			0
	PRCSTRT	32590000, 1720320, 0.125
	PRCENTRY2		xfft_r4dwpn_1680K_3840_4, 5829632, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_1680K_3840_4, 5580288, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	DD			0
	PRCSTRT	33450000, 1769472, 0.129
	PRCENTRY2		xfft_r4dwpn_1728K_4608_4, 6305280, K10
	PRCENTRY2		xfft_r4dwpn_1728K_2304_4, 6273024, I7_64 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_1728K_2304_2, 6273024, K8
	PRCENTRY2		xfft_r4dwpn_1728K_2304_1, 6273024
	PRCENTRY2		xfft_r4delay_1728K_4608_4, 5933568, I7_64 + K10
	PRCENTRY2		xfft_r4delay_1728K_2304_4, 5876736, CORE2_64 + K8
	PRCENTRY2		xfft_r4delay_1728K_2304_2, 5876736
	PRCENTRY2		xfft_r4delay_1728K_2304_1, 5876736
	DD			0
	PRCSTRT	34750000, 1835008, 0.135
	PRCENTRY2		xfft_r4dwpn_1792K_12_4, 6228992, K10
	PRCENTRY2		xfft_r4dwpn_1792K_12_2, 6228992, I7_32 + K8
	PRCENTRY2		xfft_r4dwpn_1792K_12_1, 6228992, CORE2_512
	PRCENTRY2		xfft_r4dwpn_1792K_11_4, 6197760, I7_64 + CORE2_64
	PRCENTRY2A		xfft_r4dwpn_1792K_11_4_CORE, 6197760, P4_1024
	PRCENTRY2		xfft_r4dwpn_1792K_10_4, 3608064, CORE2_32
	PRCENTRY2		xfft_r4delay_1792K_12_4, 5963264, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_1792K_12_2, 5963264, I7_32 + CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4delay_1792K_12_1, 5963264, CORE2_512_32
	PRCENTRY2		xfft_r4delay_1792K_11_4, 5913600, P4_1024_32 + I7_64
	PRCENTRY2A		xfft_r4delay_1792K_11_4_CORE, 5913600, P4_1024_64
	PRCENTRY2		xfft_r4delay_1792K_10_4, 3336192, CORE2_64
	PRCENTRY2		xfft_hg_1792K_12_4, 3002112, P4TP_512
	PRCENTRY2		xfft_hg_1792K_12_1, 3002112, K8_32
	PRCENTRY2		xfft_hg_1792K_11_1, 3834368
	DD			0
	PRCSTRT	37130000, 1966080, 0.145
	PRCENTRY2		xfft_r4dwpn_1920K_7680_4, 11976704
	PRCENTRY2		xfft_r4dwpn_1920K_6144_4, 7550976, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_1920K_5120_4, 7034368
	PRCENTRY2		xfft_r4dwpn_1920K_5120_2, 7034368, K8
	PRCENTRY2		xfft_r4dwpn_1920K_5120_1, 7034368, CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_1920K_3840_4, 6387712
	PRCENTRY2		xfft_r4dwpn_1920K_3840_2, 6387712
	PRCENTRY2		xfft_r4dwpn_1920K_3840_1, 6387712
	PRCENTRY2		xfft_r4dwpn_1920K_3072_4, 7482880, P4_1024_64 + I7 + CORE2 + K10_64
	PRCENTRY2		xfft_r4dwpn_1920K_2560_4, 6979584
	PRCENTRY2		xfft_r4dwpn_1920K_1536_4, 3911680, K10_32
	PRCENTRY2		xfft_r4delay_1920K_7680_4, 9645056
	PRCENTRY2		xfft_r4delay_1920K_6144_4, 6958592, P4_1024_32
	PRCENTRY2		xfft_r4delay_1920K_5120_4, 6621696
	PRCENTRY2		xfft_r4delay_1920K_5120_2, 6621696
	PRCENTRY2		xfft_r4delay_1920K_5120_1, 6621696, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4delay_1920K_3840_4, 6199296, K8_32
	PRCENTRY2		xfft_r4delay_1920K_3840_2, 6199296
	PRCENTRY2		xfft_r4delay_1920K_3840_1, 6199296
	PRCENTRY2		xfft_r4delay_1920K_3072_4, 6851584, P4_1024_64 + I7 + CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_1920K_2560_4, 6538240
	PRCENTRY2		xfft_r4delay_1920K_1536_4, 4122624, K10_32
	DD			0
	PRCSTRT	38660000, 2048000, 0.151
	PRCENTRY2		xfft_r4dwpn_2000K_6400_4, 7831552, P4_1024_64
	PRCENTRY2		xfft_r4delay_2000K_6400_4, 7214592, P4_1024_64
	DD			0
	PRCSTRT	38880000, 2064384, 0.153
	PRCENTRY2		xfft_r4dwpn_2016K_4608_4, 6974464
	PRCENTRY2		xfft_r4dwpn_2016K_2304_4, 6945280
	PRCENTRY2		xfft_r4delay_2016K_4608_4, 6675968
	PRCENTRY2		xfft_r4delay_2016K_2304_4, 6624256
	DD			0
	PRCSTRT	39530000, 2097152, 0.155
	PRCENTRY2		xfft_r4dwpn_2M_13_4, 12800000, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_2M_12_4, 6823936, CORE2 + K8_64 + K10
	PRCENTRY2		xfft_r4dwpn_2M_12_2, 6823936, CORE2_512 + K8_32
	PRCENTRY2		xfft_r4dwpn_2M_12_1, 6823936
	PRCENTRY2		xfft_r4dwpn_2M_11_4, 6795264, I7
	PRCENTRY2		xfft_r4dwpn_2M_10_4, 3817472
	PRCENTRY2		xfft_r4delay_2M_13_4, 10312704
	PRCENTRY2		xfft_r4delay_2M_12_4, 6623232, CORE2 + K8 + K10
	PRCENTRY2		xfft_r4delay_2M_12_2, 6623232, CORE2_512
	PRCENTRY2		xfft_r4delay_2M_12_1, 6623232
	PRCENTRY2		xfft_r4delay_2M_11_4, 6578176, I7
	PRCENTRY2		xfft_r4delay_2M_10_4, 3616768
	PRCENTRY2		xfft_hg_2048K_12_4, 3397632, P4_1024 + P4TP_512
	PRCENTRY2		xfft_hg_2048K_12_2, 3397632
	PRCENTRY2		xfft_hg_2048K_12_1, 3397632
	DD			0
	PRCSTRT	43220000, 2293760, 0.175
	PRCENTRY2		xfft_r4dwpn_2240K_5120_4, 7777280, I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_2240K_5120_2, 7777280, P4TP_512 + I7_32
	PRCENTRY2		xfft_r4dwpn_2240K_5120_1, 7777280
	PRCENTRY2		xfft_r4dwpn_2240K_2560_4, 7725568
	PRCENTRY2A		xfft_r4dwpn_2240K_2560_4_CORE, 7725568, P4_1024
	PRCENTRY2		xfft_r4delay_2240K_5120_4, 7446016, I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_2240K_5120_2, 7446016, P4TP_512 + I7_32
	PRCENTRY2		xfft_r4delay_2240K_5120_1, 7446016
	PRCENTRY2		xfft_r4delay_2240K_2560_4, 7367680
	PRCENTRY2A		xfft_r4delay_2240K_2560_4_CORE, 7367680, P4_1024
	DD			0
	PRCSTRT	44380000, 2359296, 0.181
	PRCENTRY2		xfft_r4dwpn_2304K_9216_4, 14376960, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_2304K_6144_4, 8439296
	PRCENTRY2		xfft_r4dwpn_2304K_4608_4, 7643136, P4_1024_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_2304K_4608_2, 7643136, P4TP_512
	PRCENTRY2		xfft_r4dwpn_2304K_3072_4, 8372224, I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_2304K_3072_2, 8372224
	PRCENTRY2		xfft_r4dwpn_2304K_3072_1, 8372224, CORE2_512 + K8_32
	PRCENTRY2		xfft_r4dwpn_2304K_2304_4, 7616512, CORE2_32 + K8_64
	PRCENTRY2		xfft_r4dwpn_2304K_2304_2, 7616512
	PRCENTRY2		xfft_r4dwpn_2304K_2304_1, 7616512
	PRCENTRY2		xfft_r4delay_2304K_9216_4, 11578368, P4_1024_64
	PRCENTRY2		xfft_r4delay_2304K_6144_4, 7944704
	PRCENTRY2		xfft_r4delay_2304K_4608_4, 7417856, P4_1024_32 + K10
	PRCENTRY2		xfft_r4delay_2304K_4608_2, 7417856, P4TP_512 + CORE2_512_32
	PRCENTRY2		xfft_r4delay_2304K_3072_4, 7840768, I7 + CORE2
	PRCENTRY2		xfft_r4delay_2304K_3072_2, 7840768
	PRCENTRY2		xfft_r4delay_2304K_3072_1, 7840768, CORE2_512_64 + K8_32
	PRCENTRY2		xfft_r4delay_2304K_2304_4, 7370752, K8_64
	PRCENTRY2		xfft_r4delay_2304K_2304_2, 7370752
	PRCENTRY2		xfft_r4delay_2304K_2304_1, 7370752
	DD			0
	PRCSTRT	46170000, 2457600, 0.191
	PRCENTRY2		xfft_r4dwpn_2400K_7680_4, 9406464, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_2400K_6400_4, 8756736
	PRCENTRY2		xfft_r4dwpn_2400K_6400_2, 8756736
	PRCENTRY2		xfft_r4dwpn_2400K_6400_1, 8756736, K8_64
	PRCENTRY2		xfft_r4dwpn_2400K_3840_4, 9332224, P4_1024_64 + I7_64 + CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_2400K_7680_4, 8666624, P4_1024_32
	PRCENTRY2		xfft_r4delay_2400K_6400_4, 8241664
	PRCENTRY2		xfft_r4delay_2400K_6400_2, 8241664, K8_64
	PRCENTRY2		xfft_r4delay_2400K_6400_1, 8241664
	PRCENTRY2		xfft_r4delay_2400K_3840_4, 8541184, P4_1024_64 + I7_64 + CORE2 + K10_64
	DD			0
	PRCSTRT	49250000, 2621440, 0.204
	PRCENTRY2		xfft_r4dwpn_2560K_10240_4, 15994880
	PRCENTRY2		xfft_r4dwpn_2560K_13_4, 10057728
	PRCENTRY2		xfft_r4dwpn_2560K_12_4, 9965056, P4_1024 + I7 + CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_2560K_12_2, 9965056, P4TP_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_2560K_12_1, 9965056, CORE2_512 + K8_32
	PRCENTRY2		xfft_r4dwpn_2560K_2560_4, 8470528
	PRCENTRY2		xfft_r4dwpn_2560K_11_4, 5193728
	PRCENTRY2		xfft_r4dwpn_2560K_10_4, 5203968
	PRCENTRY2		xfft_r4delay_2560K_10240_4, 12884992
	PRCENTRY2		xfft_r4delay_2560K_13_4, 9268736
	PRCENTRY2		xfft_r4delay_2560K_12_4, 9120768, P4_1024 + I7_32 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_2560K_12_2, 9120768
	PRCENTRY2		xfft_r4delay_2560K_12_1, 9120768, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_2560K_2560_4, 8196096, I7_64
	PRCENTRY2		xfft_r4delay_2560K_11_4, 5478400
	PRCENTRY2		xfft_r4delay_2560K_10_4, 5513216
	PRCENTRY2		xfft_hg_2560K_13_1, 3747584
	PRCENTRY2		xfft_hg_2560K_12_4, 4194560, P4TP_512
	PRCENTRY2		xfft_hg_2560K_12_1, 4194560
	DD			0
	PRCSTRT	51660000, 2752512, 0.218
	PRCENTRY2		xfft_r4dwpn_2688K_6144_4, 9329664
	PRCENTRY2		xfft_r4dwpn_2688K_6144_2, 9329664
	PRCENTRY2		xfft_r4dwpn_2688K_3072_4, 9265664, P4_1024 + I7_64 + CORE2_64 + K10
	PRCENTRY2		xfft_r4dwpn_2688K_3072_2, 9265664, P4TP_512 + I7_32 + CORE2_32
	PRCENTRY2		xfft_r4dwpn_2688K_3072_1, 9265664, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_2688K_6144_4, 8932864, P4_1024_32 + K10_32
	PRCENTRY2		xfft_r4delay_2688K_6144_2, 8932864
	PRCENTRY2		xfft_r4delay_2688K_3072_4, 8834048, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_2688K_3072_2, 8834048, P4TP_512 + I7_32 + CORE2_32
	PRCENTRY2		xfft_r4delay_2688K_3072_1, 8834048, CORE2_512 + K8
	DD			0
	PRCSTRT	53740000, 2867200, 0.232
	PRCENTRY2		xfft_r4dwpn_2800K_6400_4, 9683968, I7_64 + CORE2 + K10_64
	PRCENTRY2		xfft_r4dwpn_2800K_6400_2, 9683968, CORE2_512_64
	PRCENTRY2		xfft_r4delay_2800K_6400_4, 9270784, CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_2800K_6400_2, 9270784, CORE2_512_64 + I7_64
	DD			0
	PRCSTRT	55180000, 2949120, 0.240
	PRCENTRY2		xfft_r4dwpn_2880K_9216_4, 11290624, CORE2_64
	PRCENTRY2		xfft_r4dwpn_2880K_7680_4, 10515968
	PRCENTRY2		xfft_r4dwpn_2880K_4608_4, 11177472, P4_1024 + K10
	PRCENTRY2		xfft_r4dwpn_2880K_4608_2, 11177472, P4TP_512 + CORE2_512_64 + I7
	PRCENTRY2		xfft_r4dwpn_2880K_3840_4, 10442752
	PRCENTRY2		xfft_r4dwpn_2880K_3840_2, 10442752, K8_64
	PRCENTRY2		xfft_r4dwpn_2880K_3840_1, 10442752
	PRCENTRY2		xfft_r4dwpn_2880K_2304_4, 5810176, CORE2_32
	PRCENTRY2		xfft_r4delay_2880K_9216_4, 10403328, CORE2
	PRCENTRY2		xfft_r4delay_2880K_7680_4, 9898496
	PRCENTRY2		xfft_r4delay_2880K_4608_4, 10226688, P4_1024 + K10
	PRCENTRY2		xfft_r4delay_2880K_4608_2, 10226688, P4TP_512 + CORE2_512_64 + I7
	PRCENTRY2		xfft_r4delay_2880K_3840_4, 9776128, K8_64
	PRCENTRY2		xfft_r4delay_2880K_3840_2, 9776128
	PRCENTRY2		xfft_r4delay_2880K_3840_1, 9776128
	PRCENTRY2		xfft_r4delay_2880K_2304_4, 6131712
	DD			0
	PRCSTRT	58850000, 3145728, 0.259
	PRCENTRY2		xfft_r4dwpn_3M_12288_4, 19185664
	PRCENTRY2		xfft_r4dwpn_3M_13_4, 11240960
	PRCENTRY2		xfft_r4dwpn_3M_6144_4, 10219520
	PRCENTRY2		xfft_r4dwpn_3M_6144_2, 10219520, P4TP_512
	PRCENTRY2		xfft_r4dwpn_3M_12_4, 11149312, P4_1024_64 + I7_32 + K10
	PRCENTRY2		xfft_r4dwpn_3M_12_2, 11149312
	PRCENTRY2		xfft_r4dwpn_3M_12_1, 11149312, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_3M_3072_4, 10158080, I7_64 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_3M_11_4, 5083136
	PRCENTRY2		xfft_r4dwpn_3M_10_4, 5103616, CORE2_32
	PRCENTRY2		xfft_r4delay_3M_12288_4, 15453184
	PRCENTRY2		xfft_r4delay_3M_13_4, 10582528
	PRCENTRY2		xfft_r4delay_3M_6144_4, 9920512, K10_32
	PRCENTRY2		xfft_r4delay_3M_6144_2, 9920512, P4TP_512
	PRCENTRY2		xfft_r4delay_3M_12_4, 10437632, I7_32 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_3M_12_2, 10437632, CORE2_512_32 + K8_64
	PRCENTRY2		xfft_r4delay_3M_12_1, 10437632, CORE2_512_64 + K8_32
	PRCENTRY2		xfft_r4delay_3M_3072_4, 9826304, I7_64 + CORE2_32
	PRCENTRY2		xfft_r4delay_3M_11_4, 5496832
	PRCENTRY2		xfft_r4delay_3M_10_4, 5550080
	PRCENTRY2		xfft_hg_3072K_13_1, 4406912
	PRCENTRY2		xfft_hg_3072K_12_4, 4989312, P4_1024_32
	PRCENTRY2A		xfft_hg_3072K_12_4_CORE, 4989312, P4_1024_64
	PRCENTRY2		xfft_hg_3072K_12_1, 4989312
	DD			0
	PRCSTRT	61220000, 3276800, 0.275
	PRCENTRY2		xfft_r4dwpn_3200K_12800_4, 19931136
	PRCENTRY2		xfft_r4dwpn_3200K_10240_4, 12564480
	PRCENTRY2		xfft_r4dwpn_3200K_5120_4, 12447232, P4_1024 + CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_3200K_5120_2, 12447232, P4TP_512 + I7
	PRCENTRY2		xfft_r4dwpn_3200K_5120_1, 12447232, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_3200K_2560_4, 6459392
	PRCENTRY2		xfft_r4delay_3200K_12800_4, 16043008
	PRCENTRY2		xfft_r4delay_3200K_10240_4, 11578880
	PRCENTRY2		xfft_r4delay_3200K_5120_4, 11389952, P4_1024 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_3200K_5120_2, 11389952, P4TP_512 + I7
	PRCENTRY2		xfft_r4delay_3200K_5120_1, 11389952, CORE2_512 + K8
	PRCENTRY2		xfft_r4delay_3200K_2560_4, 6817792
	DD			0
	PRCSTRT	64230000, 3440640, 0.294
	PRCENTRY2		xfft_r4dwpn_3360K_7680_4, 11627520, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_3360K_3840_4, 11557376, P4_1024_64 + I7_64 + CORE2
	PRCENTRY2		xfft_r4delay_3360K_7680_4, 11132416, P4_1024_32
	PRCENTRY2		xfft_r4delay_3360K_3840_4, 11015168, P4_1024_64 + I7_64 + CORE2
	DD			0
	PRCSTRT	65950000, 3538944, 0.305
	PRCENTRY2		xfft_r4dwpn_3456K_9216_4, 12621312, CORE2_64
	PRCENTRY2		xfft_r4dwpn_3456K_9216_2, 12621312
	PRCENTRY2		xfft_r4dwpn_3456K_9216_1, 12621312, K8_32
	PRCENTRY2		xfft_r4dwpn_3456K_4608_4, 12509184, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_3456K_2304_4, 5683200, K8_64
	PRCENTRY2		xfft_r4dwpn_3456K_2304_2, 5683200
	PRCENTRY2		xfft_r4delay_3456K_9216_4, 11880960, CORE2_64
	PRCENTRY2		xfft_r4delay_3456K_9216_2, 11880960
	PRCENTRY2		xfft_r4delay_3456K_9216_1, 11880960, K8
	PRCENTRY2		xfft_r4delay_3456K_4608_4, 11707392, I7_64 + K10
	PRCENTRY2		xfft_r4delay_3456K_2304_4, 6150144
	PRCENTRY2		xfft_r4delay_3456K_2304_2, 6150144
	DD			0
	PRCSTRT	68490000, 3670016, 0.323
	PRCENTRY2		xfft_r4dwpn_3584K_13_4, 12426240
	PRCENTRY2		xfft_r4dwpn_3584K_12_4, 12337664, P4_1024 + I7 + K10
	PRCENTRY2		xfft_r4dwpn_3584K_12_2, 12337664, CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_3584K_12_1, 12337664, CORE2_512_32 + K8
	PRCENTRY2		xfft_r4dwpn_3584K_11_4, 7142912
	PRCENTRY2		xfft_r4dwpn_3584K_10_4, 6043648, CORE2
	PRCENTRY2		xfft_r4delay_3584K_13_4, 11898368
	PRCENTRY2		xfft_r4delay_3584K_12_4, 11758592, P4_1024 + I7 + CORE2_64 + K10
	PRCENTRY2		xfft_r4delay_3584K_12_2, 11758592, CORE2_512_64
	PRCENTRY2		xfft_r4delay_3584K_12_1, 11758592, CORE2_512_32 + K8
	PRCENTRY2		xfft_r4delay_3584K_11_4, 6559744
	PRCENTRY2		xfft_r4delay_3584K_10_4, 6627328, CORE2_32
	PRCENTRY2		xfft_hg_3584K_13_1, 5066496, P4TP_512
	PRCENTRY2		xfft_hg_3584K_12_4, 5784064
	PRCENTRY2		xfft_hg_3584K_12_2, 5784064
	PRCENTRY2		xfft_hg_3584K_12_1, 5784064
	DD			0
	PRCSTRT	73180000, 3932160, 0.352
	PRCENTRY2		xfft_r4dwpn_3840K_15360_4, 23977984
	PRCENTRY2		xfft_r4dwpn_3840K_12288_4, 15067136
	PRCENTRY2		xfft_r4dwpn_3840K_10240_4, 14042624
	PRCENTRY2		xfft_r4dwpn_3840K_7680_4, 12738560
	PRCENTRY2		xfft_r4dwpn_3840K_6144_4, 14933504, P4_1024 + P4TP_512 + K10_32
	PRCENTRY2		xfft_r4dwpn_3840K_5120_4, 13926400, I7 + K10_64
	PRCENTRY2		xfft_r4dwpn_3840K_5120_2, 13926400
	PRCENTRY2		xfft_r4dwpn_3840K_5120_1, 13926400, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_3840K_3840_4, 12670976, CORE2
	PRCENTRY2		xfft_r4dwpn_3840K_3840_2, 12670976
	PRCENTRY2		xfft_r4dwpn_3840K_3840_1, 12670976, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_3840K_3072_4, 7737344
	PRCENTRY2		xfft_r4delay_3840K_15360_4, 19311616
	PRCENTRY2		xfft_r4delay_3840K_12288_4, 13884928, CORE2_32
	PRCENTRY2		xfft_r4delay_3840K_10240_4, 13220352
	PRCENTRY2		xfft_r4delay_3840K_7680_4, 12365824
	PRCENTRY2		xfft_r4delay_3840K_6144_4, 13663232, P4_1024 + P4TP_512 + K10_32
	PRCENTRY2		xfft_r4delay_3840K_5120_4, 13034496, I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_3840K_5120_2, 13034496, K8_64
	PRCENTRY2		xfft_r4delay_3840K_5120_1, 13034496, CORE2_512_64 + K8_32
	PRCENTRY2		xfft_r4delay_3840K_3840_4, 12253184
	PRCENTRY2		xfft_r4delay_3840K_3840_2, 12253184
	PRCENTRY2		xfft_r4delay_3840K_3840_1, 12253184, CORE2_512_32
	PRCENTRY2		xfft_r4delay_3840K_3072_4, 8169472
	DD			0
	PRCSTRT	76210000, 4096000, 0.371
	PRCENTRY2		xfft_r4dwpn_4000K_12800_4, 15640576
	PRCENTRY2		xfft_r4dwpn_4000K_6400_4, 15521280, P4_1024_64
	PRCENTRY2		xfft_r4delay_4000K_12800_4, 14409216
	PRCENTRY2		xfft_r4delay_4000K_6400_4, 14197760, P4_1024_64
	DD			0
	PRCSTRT	76790000, 4128768, 0.375
	PRCENTRY2		xfft_r4dwpn_4032K_9216_4, 13954048
	PRCENTRY2		xfft_r4dwpn_4032K_4608_4, 13844992, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_4032K_2304_4, 8005120
	PRCENTRY2		xfft_r4delay_4032K_9216_4, 13360640
	PRCENTRY2		xfft_r4delay_4032K_4608_4, 13192192, P4_1024_32
	PRCENTRY2		xfft_r4delay_4032K_2304_4, 7344128
	DD			0
	PRCSTRT	77990000, 4194304, 0.382
	PRCENTRY2		xfft_r4dwpn_4M_14_4, 25575424
	PRCENTRY2		xfft_r4dwpn_4M_13_4, 13611008
	PRCENTRY2		xfft_r4dwpn_4M_13_2, 13611008, P4TP_512
	PRCENTRY2		xfft_r4dwpn_4M_13_1, 13611008
	PRCENTRY2		xfft_r4dwpn_4M_12_4, 13524992, P4_1024_64 + I7 + CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_4M_12_2, 13524992, K8
	PRCENTRY2		xfft_r4dwpn_4M_12_1, 13524992, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_4M_11_4, 7548928
	PRCENTRY2		xfft_r4dwpn_4M_10_4, 6467584
	PRCENTRY2		xfft_r4dwpn_4M_10_2, 6467584
	PRCENTRY2		xfft_r4delay_4M_14_4, 20597760
	PRCENTRY2		xfft_r4delay_4M_13_4, 13213696
	PRCENTRY2		xfft_r4delay_4M_13_2, 13213696
	PRCENTRY2		xfft_r4delay_4M_13_1, 13213696
	PRCENTRY2		xfft_r4delay_4M_12_4, 13078528, P4_1024_64 + I7_64 + CORE2 + K8 + K10_32
	PRCENTRY2		xfft_r4delay_4M_12_2, 13078528, I7_32 + K10_64
	PRCENTRY2		xfft_r4delay_4M_12_1, 13078528, CORE2_512
	PRCENTRY2		xfft_r4delay_4M_11_4, 7102464
	PRCENTRY2		xfft_r4delay_4M_10_4, 7188480
	PRCENTRY2		xfft_r4delay_4M_10_2, 7188480
	PRCENTRY2		xfft_hg_4096K_13_1, 5724160, P4TP_512
	PRCENTRY2		xfft_hg_4096K_12_4, 6574080, P4_1024_32
	PRCENTRY2		xfft_hg_4096K_12_2, 6574080
	PRCENTRY2		xfft_hg_4096K_12_1, 6574080
	DD			0
	PRCSTRT	85200000, 4587520, 0.436
	PRCENTRY2		xfft_r4dwpn_4480K_10240_4, 15522816
	PRCENTRY2		xfft_r4dwpn_4480K_5120_4, 15409664, P4_1024 + I7 + CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_4480K_10240_4, 14863872
	PRCENTRY2		xfft_r4delay_4480K_5120_4, 14683136, P4_1024 + I7 + CORE2 + K10_64
	DD			0
	PRCSTRT	87400000, 4718592, 0.454
	PRCENTRY2		xfft_r4dwpn_4608K_12288_4, 16840192
	PRCENTRY2		xfft_r4dwpn_4608K_9216_4, 15286272, CORE2_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_4608K_6144_4, 16707584, K10
	PRCENTRY2		xfft_r4dwpn_4608K_4608_4, 15179776, P4_1024 + I7
	PRCENTRY2		xfft_r4dwpn_4608K_3072_4, 7561216, K8_32
	PRCENTRY2		xfft_r4dwpn_4608K_2304_4, 8460288, CORE2_32
	PRCENTRY2		xfft_r4delay_4608K_12288_4, 15854080
	PRCENTRY2		xfft_r4delay_4608K_9216_4, 14839808, CORE2_64
	PRCENTRY2		xfft_r4delay_4608K_6144_4, 15635456, K10
	PRCENTRY2		xfft_r4delay_4608K_4608_4, 14675968, P4_1024 + I7
	PRCENTRY2		xfft_r4delay_4608K_3072_4, 8187904, K8
	PRCENTRY2		xfft_r4delay_4608K_2304_4, 7952384, CORE2_32
	DD			0
	PRCSTRT	91020000, 4915200, 0.480
	PRCENTRY2		xfft_r4dwpn_4800K_15360_4, 18827264
	PRCENTRY2		xfft_r4dwpn_4800K_12800_4, 17487360
	PRCENTRY2		xfft_r4dwpn_4800K_7680_4, 18632192, P4_1024 + I7_32
	PRCENTRY2		xfft_r4dwpn_4800K_6400_4, 17369088, I7_64 + CORE2_64 + K8_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_4800K_3840_4, 9635840, CORE2_32 + K8_32
	PRCENTRY2		xfft_r4delay_4800K_15360_4, 17350144
	PRCENTRY2		xfft_r4delay_4800K_12800_4, 16460288
	PRCENTRY2		xfft_r4delay_4800K_7680_4, 17042432, P4_1024 + I7_32 + CORE2_64
	PRCENTRY2		xfft_r4delay_4800K_6400_4, 16251904, I7_64 + CORE2_32 + K8_64 + K10_64
	PRCENTRY2		xfft_r4delay_4800K_3840_4, 10178560, K8_32
	DD			0
	PRCSTRT	97020000, 5242880, 0.485
	PRCENTRY2		xfft_r4dwpn_5M_20480_4, 31965184
	PRCENTRY2		xfft_r4dwpn_5M_14_4, 20080640
	PRCENTRY2		xfft_r4dwpn_5M_10240_4, 17002496
	PRCENTRY2		xfft_r4dwpn_5M_13_4, 19897856, P4_1024_64 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_5M_5120_4, 16891904, I7 + K8
	PRCENTRY2		xfft_r4dwpn_5M_12_4, 10285056, CORE2_32 + K10
	PRCENTRY2		xfft_r4dwpn_5M_11_4, 10278912
	PRCENTRY2		xfft_r4delay_5M_20480_4, 25742336
	PRCENTRY2		xfft_r4delay_5M_14_4, 18505216, CORE2_64
	PRCENTRY2		xfft_r4delay_5M_10240_4, 16506880
	PRCENTRY2		xfft_r4delay_5M_13_4, 18201600, P4_1024_64
	PRCENTRY2		xfft_r4delay_5M_5120_4, 16330752, I7 + K8_64
	PRCENTRY2		xfft_r4delay_5M_12_4, 10864640, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_5M_11_4, 10866688
	PRCENTRY2		xfft_hg_5M_12_4, 8167040, P4_1024_32 + K8_32
	PRCENTRY2		xfft_hg_5M_12_1, 8167040
	DD			0
	PRCSTRT	101700000, 5505024, 0.561
	PRCENTRY2		xfft_r4dwpn_5376K_12288_4, 18615296, K8_64
	PRCENTRY2		xfft_r4dwpn_5376K_6144_4, 18485760, P4_1024 + I7 + K8_32 + K10
	PRCENTRY2		xfft_r4dwpn_5376K_3072_4, 10669568, CORE2
	PRCENTRY2		xfft_r4delay_5376K_12288_4, 17825280, K8
	PRCENTRY2		xfft_r4delay_5376K_6144_4, 17611776, P4_1024 + I7 + K10
	PRCENTRY2		xfft_r4delay_5376K_3072_4, 9775104, CORE2
	DD			0
	PRCSTRT	105900000, 5734400, 0.592
	PRCENTRY2		xfft_r4dwpn_5600K_12800_4, 19336192
	PRCENTRY2		xfft_r4dwpn_5600K_6400_4, 19220992, P4_1024_64 + I7 + CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_5600K_12800_4, 18513408, P4_1024_32
	PRCENTRY2		xfft_r4delay_5600K_6400_4, 18310144, P4_1024_64 + I7 + CORE2 + K10_64
	DD			0
	PRCSTRT	108700000, 5898240, 0.614
	PRCENTRY2		xfft_r4dwpn_5760K_15360_4, 21042688
	PRCENTRY2		xfft_r4dwpn_5760K_9216_4, 22359552, P4_1024 + I7_32 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_5760K_7680_4, 20848640
	PRCENTRY2		xfft_r4dwpn_5760K_4608_4, 11530240, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_5760K_3840_4, 9410560, K8
	PRCENTRY2		xfft_r4dwpn_5760K_2304_4, 11526144, CORE2_32
	PRCENTRY2		xfft_r4delay_5760K_15360_4, 19810816
	PRCENTRY2		xfft_r4delay_5760K_9216_4, 20450304, P4_1024 + CORE2_64
	PRCENTRY2		xfft_r4delay_5760K_7680_4, 19506176, I7_64
	PRCENTRY2		xfft_r4delay_5760K_4608_4, 12183552, K10
	PRCENTRY2		xfft_r4delay_5760K_3840_4, 10196992, K8
	PRCENTRY2		xfft_r4delay_5760K_2304_4, 12183552, I7_32 + CORE2_32
	DD			0
	PRCSTRT	115800000, 6291456, 0.668
	PRCENTRY2		xfft_r4dwpn_6M_14_4, 22443520
	PRCENTRY2		xfft_r4dwpn_6M_12288_4, 20389888
	PRCENTRY2		xfft_r4dwpn_6M_13_4, 22261760
	PRCENTRY2		xfft_r4dwpn_6M_6144_4, 20262912, P4_1024 + I7 + CORE2_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_6M_12_4, 10043392, CORE2_32 + K8 + K10_64
	PRCENTRY2		xfft_r4dwpn_6M_3072_4, 11272192
	PRCENTRY2		xfft_r4dwpn_6M_11_4, 10047488
	PRCENTRY2		xfft_r4delay_6M_14_4, 21129728
	PRCENTRY2		xfft_r4delay_6M_12288_4, 19795968
	PRCENTRY2		xfft_r4delay_6M_13_4, 20829184, CORE2_64
	PRCENTRY2		xfft_r4delay_6M_6144_4, 19587072, P4_1024 + I7_32 + K10_32
	PRCENTRY2		xfft_r4delay_6M_12_4, 10883072, I7_64 + CORE2_32 + K8 + K10_64
	PRCENTRY2		xfft_r4delay_6M_3072_4, 10579968
	PRCENTRY2		xfft_r4delay_6M_11_4, 10903552
	PRCENTRY2		xfft_hg_6M_12_4, 9756288
	PRCENTRY2		xfft_hg_6M_12_2, 9756288
	PRCENTRY2		xfft_hg_6M_12_1, 9756288
	DD			0
	PRCSTRT	120500000, 6553600, 0.715
	PRCENTRY2		xfft_r4dwpn_6400K_25600_4, 39919616
	PRCENTRY2		xfft_r4dwpn_6400K_20480_4, 25094144
	PRCENTRY2		xfft_r4dwpn_6400K_12800_4, 21184512
	PRCENTRY2		xfft_r4dwpn_6400K_10240_4, 24862208, CORE2_32
	PRCENTRY2		xfft_r4dwpn_6400K_6400_4, 21071872, I7 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_6400K_5120_4, 12832768, K10
	PRCENTRY2		xfft_r4delay_6400K_25600_4, 32140288
	PRCENTRY2		xfft_r4delay_6400K_20480_4, 23125504
	PRCENTRY2		xfft_r4delay_6400K_12800_4, 20566016
	PRCENTRY2		xfft_r4delay_6400K_10240_4, 22739968, CORE2_64
	PRCENTRY2		xfft_r4delay_6400K_6400_4, 20367360, I7_32 + CORE2_32
	PRCENTRY2		xfft_r4delay_6400K_5120_4, 13559808, I7_64 + K10
	DD			0
	PRCSTRT	126400000, 6881280, 0.774
	PRCENTRY2		xfft_r4dwpn_6720K_15360_4, 23260160
	PRCENTRY2		xfft_r4dwpn_6720K_7680_4, 23069184, I7 + K10_64
	PRCENTRY2		xfft_r4dwpn_6720K_3840_4, 13305344, CORE2
	PRCENTRY2		xfft_r4delay_6720K_15360_4, 22273536
	PRCENTRY2		xfft_r4delay_6720K_7680_4, 21974016, I7 + K10_64
	PRCENTRY2		xfft_r4delay_6720K_3840_4, 12177408, CORE2
	DD			0
	PRCSTRT	129700000, 7077888, 0.809
	PRCENTRY2		xfft_r4dwpn_6912K_9216_4, 25018368, I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_6912K_4608_4, 11255808, K10_32
	PRCENTRY2		xfft_r4dwpn_6912K_2304_4, 11261952, CORE2_32
	PRCENTRY2		xfft_r4delay_6912K_9216_4, 23405568, CORE2_64 + K10_64
	PRCENTRY2		xfft_r4delay_6912K_4608_4, 12201984, I7_64 + K10_32
	PRCENTRY2		xfft_r4delay_6912K_2304_4, 12220416, CORE2_32
	DD			0
	PRCSTRT	134800000, 7340032, 0.886
	PRCENTRY2		xfft_r4dwpn_7M_14_4, 24808448
	PRCENTRY2		xfft_r4dwpn_7M_13_4, 24629760, I7_32
	PRCENTRY2		xfft_r4dwpn_7M_12_4, 14200320, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_7M_11_4, 11905024, I7_64
	PRCENTRY2		xfft_r4delay_7M_14_4, 23756288
	PRCENTRY2		xfft_r4delay_7M_13_4, 23460864, I7_32
	PRCENTRY2		xfft_r4delay_7M_12_4, 12994560, I7_64 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_7M_11_4, 13029376
	PRCENTRY2		xfft_hg_7M_12_4, 11345536
	PRCENTRY2		xfft_hg_7M_12_2, 11345536
	PRCENTRY2		xfft_hg_7M_12_1, 11345536
	DD			0
	PRCSTRT	144000000, 7864320, 0.949
	PRCENTRY2		xfft_r4dwpn_7680K_20480_4, 28046848
	PRCENTRY2		xfft_r4dwpn_7680K_15360_4, 25477120
	PRCENTRY2		xfft_r4dwpn_7680K_12288_4, 29822464
	PRCENTRY2		xfft_r4dwpn_7680K_10240_4, 27815936
	PRCENTRY2		xfft_r4dwpn_7680K_7680_4, 25288704, I7
	PRCENTRY2		xfft_r4dwpn_7680K_6144_4, 15384576, K10
	PRCENTRY2		xfft_r4dwpn_7680K_5120_4, 12525568
	PRCENTRY2		xfft_r4dwpn_7680K_3840_4, 14055424
	PRCENTRY2		xfft_r4dwpn_7680K_3072_4, 15345664, CORE2
	PRCENTRY2		xfft_r4delay_7680K_20480_4, 26405376
	PRCENTRY2		xfft_r4delay_7680K_15360_4, 24735744
	PRCENTRY2		xfft_r4delay_7680K_12288_4, 27274240, CORE2_32
	PRCENTRY2		xfft_r4delay_7680K_10240_4, 26022912
	PRCENTRY2		xfft_r4delay_7680K_7680_4, 24440832, I7_64
	PRCENTRY2		xfft_r4delay_7680K_6144_4, 16259072, K10_32
	PRCENTRY2		xfft_r4delay_7680K_5120_4, 13578240, K10_64
	PRCENTRY2		xfft_r4delay_7680K_3840_4, 13178880
	PRCENTRY2		xfft_r4delay_7680K_3072_4, 16211968, I7_32 + CORE2_64
	DD			0
	PRCSTRT	149900000, 8192000, 1.007
	PRCENTRY2		xfft_r4dwpn_8000K_25600_4, 31328256
	PRCENTRY2		xfft_r4dwpn_8000K_12800_4, 31010304, I7_64
	PRCENTRY2		xfft_r4dwpn_8000K_6400_4, 15988736, K10_64
	PRCENTRY2		xfft_r4delay_8000K_25600_4, 28868096
	PRCENTRY2		xfft_r4delay_8000K_12800_4, 28355584, I7_64
	PRCENTRY2		xfft_r4delay_8000K_6400_4, 16900096, K10_64
	DD			0
	PRCSTRT	151000000, 8257536, 1.019
	PRCENTRY2		xfft_r4dwpn_8064K_9216_4, 27681280
	PRCENTRY2		xfft_r4dwpn_8064K_4608_4, 15937024
	PRCENTRY2		xfft_r4dwpn_8064K_2304_4, 13348864
	PRCENTRY2		xfft_r4delay_8064K_9216_4, 26364928
	PRCENTRY2		xfft_r4delay_8064K_4608_4, 14575616
	PRCENTRY2		xfft_r4delay_8064K_2304_4, 14608384
	DD			0
	PRCSTRT	153300000, 8388608, 1.042
	PRCENTRY2		xfft_r4dwpn_8M_14_4, 27172864, CORE2
	PRCENTRY2		xfft_r4dwpn_8M_13_4, 26996736, I7 + K10
	PRCENTRY2		xfft_r4dwpn_8M_12_4, 14999552
	PRCENTRY2		xfft_r4dwpn_8M_11_4, 12722176
	PRCENTRY2		xfft_r4dwpn_8M_11_2, 12722176
	PRCENTRY2		xfft_r4delay_8M_14_4, 26382336, CORE2_64
	PRCENTRY2		xfft_r4delay_8M_13_4, 26091520, I7 + K10
	PRCENTRY2		xfft_r4delay_8M_12_4, 14061568, CORE2_32
	PRCENTRY2		xfft_r4delay_8M_11_4, 14114816
	PRCENTRY2		xfft_r4delay_8M_11_2, 14114816
	PRCENTRY2		xfft_hg_8M_12_4, 12926976
	PRCENTRY2		xfft_hg_8M_12_2, 12926976
	PRCENTRY2		xfft_hg_8M_12_1, 12926976
	DD			0
	PRCSTRT	167600000, 9175040, 1.110
	PRCENTRY2		xfft_r4dwpn_8960K_20480_4, 31001600
	PRCENTRY2		xfft_r4dwpn_8960K_10240_4, 30773760, I7
	PRCENTRY2		xfft_r4dwpn_8960K_5120_4, 17731072, CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_8960K_20480_4, 29687296
	PRCENTRY2		xfft_r4delay_8960K_10240_4, 29309952, I7
	PRCENTRY2		xfft_r4delay_8960K_5120_4, 16214016, CORE2 + K10_64
	DD			0
	PRCSTRT	172200000, 9437184, 1.131
	PRCENTRY2		xfft_r4dwpn_9M_12288_4, 33366016, K10_64
	PRCENTRY2		xfft_r4dwpn_9M_9216_4, 30343168, I7 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_9M_6144_4, 15011840, K10_32
	PRCENTRY2		xfft_r4dwpn_9M_4608_4, 16834560
	PRCENTRY2		xfft_r4dwpn_9M_3072_4, 14983168, CORE2_32
	PRCENTRY2		xfft_r4dwpn_9M_2304_4, 14264320
	PRCENTRY2		xfft_r4delay_9M_12288_4, 31212544
	PRCENTRY2		xfft_r4delay_9M_9216_4, 29323264, I7 + CORE2 + K10_64
	PRCENTRY2		xfft_r4delay_9M_6144_4, 16277504, K10_32
	PRCENTRY2		xfft_r4delay_9M_4608_4, 15773696
	PRCENTRY2		xfft_r4delay_9M_3072_4, 16248832
	PRCENTRY2		xfft_r4delay_9M_2304_4, 15824896
	DD			0
	PRCSTRT	179100000, 9830400, 1.166
	PRCENTRY2		xfft_r4dwpn_9600K_25600_4, 35018240
	PRCENTRY2		xfft_r4dwpn_9600K_15360_4, 37268992, I7_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_9600K_12800_4, 34701312, I7_64
	PRCENTRY2		xfft_r4dwpn_9600K_7680_4, 19181568
	PRCENTRY2		xfft_r4dwpn_9600K_6400_4, 15599616, K10_64
	PRCENTRY2		xfft_r4dwpn_9600K_3840_4, 19136512, CORE2
	PRCENTRY2		xfft_r4delay_9600K_25600_4, 32967168
	PRCENTRY2		xfft_r4delay_9600K_15360_4, 34081792, I7_32 + CORE2_64
	PRCENTRY2		xfft_r4delay_9600K_12800_4, 32457728
	PRCENTRY2		xfft_r4delay_9600K_7680_4, 20277248
	PRCENTRY2		xfft_r4delay_9600K_6400_4, 16918528, K10
	PRCENTRY2		xfft_r4delay_9600K_3840_4, 20211712, I7_64 + CORE2_32
	DD			0
	PRCSTRT	190600000, 10485760, 1.221
	PRCENTRY2		xfft_r4dwpn_10M_20480_4, 33955840
	PRCENTRY2		xfft_r4dwpn_10M_14_4, 39751168, CORE2
	PRCENTRY2		xfft_r4dwpn_10M_10240_4, 33730560, I7_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_10M_13_4, 20480000, K10_64
	PRCENTRY2		xfft_r4dwpn_10M_12_4, 20416512, I7_32
	PRCENTRY2		xfft_r4dwpn_10M_5120_4, 18726912
	PRCENTRY2		xfft_r4delay_10M_20480_4, 32968704
	PRCENTRY2		xfft_r4delay_10M_14_4, 36350976, CORE2
	PRCENTRY2		xfft_r4delay_10M_10240_4, 32595968, I7_64
	PRCENTRY2		xfft_r4delay_10M_13_4, 21649408, K10
	PRCENTRY2		xfft_r4delay_10M_12_4, 21561344, I7_32
	PRCENTRY2		xfft_r4delay_10M_5120_4, 17543168
	PRCENTRY2		xfft_hg_10M_13_4, 13639296
	DD			0
	PRCSTRT	200300000, 11010048, 1.266
	PRCENTRY2		xfft_r4dwpn_10752K_12288_4, 36913664, K10
	PRCENTRY2		xfft_r4dwpn_10752K_6144_4, 21265920
	PRCENTRY2		xfft_r4dwpn_10752K_3072_4, 17758208, I7 + CORE2
	PRCENTRY2		xfft_r4delay_10752K_12288_4, 35154944, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_10752K_6144_4, 19437568
	PRCENTRY2		xfft_r4delay_10752K_3072_4, 19423232, I7 + CORE2_64
	DD			0
	PRCSTRT	208300000, 11468800, 1.305
	PRCENTRY2		xfft_r4dwpn_11200K_25600_4, 38710272
	PRCENTRY2		xfft_r4dwpn_11200K_12800_4, 38396416, I7
	PRCENTRY2		xfft_r4dwpn_11200K_6400_4, 22115840, CORE2
	PRCENTRY2		xfft_r4delay_11200K_25600_4, 37068288
	PRCENTRY2		xfft_r4delay_11200K_12800_4, 36563968, I7
	PRCENTRY2		xfft_r4delay_11200K_6400_4, 20209664, CORE2
	DD			0
	PRCSTRT	214000000, 11796480, 1.333
	PRCENTRY2		xfft_r4dwpn_11520K_15360_4, 41697280
	PRCENTRY2		xfft_r4dwpn_11520K_9216_4, 23007232, CORE2_64 + K10
	PRCENTRY2		xfft_r4dwpn_11520K_7680_4, 18710528
	PRCENTRY2		xfft_r4dwpn_11520K_4608_4, 22923264, I7
	PRCENTRY2		xfft_r4dwpn_11520K_3840_4, 18675712, CORE2_32
	PRCENTRY2		xfft_r4delay_11520K_15360_4, 39003136
	PRCENTRY2		xfft_r4delay_11520K_9216_4, 24324096, CORE2_64 + K10
	PRCENTRY2		xfft_r4delay_11520K_7680_4, 20295680
	PRCENTRY2		xfft_r4delay_11520K_4608_4, 24207360, I7 + CORE2_32
	PRCENTRY2		xfft_r4delay_11520K_3840_4, 20248576
	DD			0
	PRCSTRT	228400000, 12582912, 1.400
	PRCENTRY2		xfft_r4dwpn_12M_14_4, 44474368, I7_64 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_12M_12288_4, 40460288, I7_32 + K10
	PRCENTRY2		xfft_r4dwpn_12M_13_4, 19976192
	PRCENTRY2		xfft_r4dwpn_12M_6144_4, 22458368
	PRCENTRY2		xfft_r4dwpn_12M_12_4, 19922944, CORE2_32
	PRCENTRY2		xfft_r4dwpn_12M_3072_4, 18968576
	PRCENTRY2		xfft_r4delay_12M_14_4, 41600000, I7_64 + CORE2_64
	PRCENTRY2		xfft_r4delay_12M_12288_4, 39096320, I7_32 + K10
	PRCENTRY2		xfft_r4delay_12M_13_4, 21667840
	PRCENTRY2		xfft_r4delay_12M_6144_4, 21028864
	PRCENTRY2		xfft_r4delay_12M_12_4, 21598208, CORE2_32
	PRCENTRY2		xfft_r4delay_12M_3072_4, 21032960
	PRCENTRY2		xfft_hg_12M_13_4, 16277120
	DD			0
	PRCSTRT	237100000, 13107200, 1.488
	PRCENTRY2		xfft_r4dwpn_12800K_25600_4, 42401792
	PRCENTRY2		xfft_r4dwpn_12800K_20480_4, 49679872
	PRCENTRY2		xfft_r4dwpn_12800K_12800_4, 42090496, I7
	PRCENTRY2		xfft_r4dwpn_12800K_10240_4, 25575424, K10
	PRCENTRY2		xfft_r4dwpn_12800K_6400_4, 23357440
	PRCENTRY2		xfft_r4dwpn_12800K_5120_4, 25487360, CORE2
	PRCENTRY2		xfft_r4delay_12800K_25600_4, 41168896
	PRCENTRY2		xfft_r4delay_12800K_20480_4, 45427712
	PRCENTRY2		xfft_r4delay_12800K_12800_4, 40669184, I7
	PRCENTRY2		xfft_r4delay_12800K_10240_4, 27039744, K10
	PRCENTRY2		xfft_r4delay_12800K_6400_4, 21866496
	PRCENTRY2		xfft_r4delay_12800K_5120_4, 26910720, CORE2
	DD			0
	PRCSTRT	249100000, 13762560, 1.597
	PRCENTRY2		xfft_r4dwpn_13440K_15360_4, 46129664, I7_32 + K10
	PRCENTRY2		xfft_r4dwpn_13440K_7680_4, 26537472
	PRCENTRY2		xfft_r4dwpn_13440K_3840_4, 22138880, I7_64 + CORE2
	PRCENTRY2		xfft_r4delay_13440K_15360_4, 43928576, I7_32 + K10
	PRCENTRY2		xfft_r4delay_13440K_7680_4, 24242176
	PRCENTRY2		xfft_r4delay_13440K_3840_4, 24209408, I7_64 + CORE2
	DD			0
	PRCSTRT	255300000, 14155776, 1.663
	PRCENTRY2		xfft_r4dwpn_13824K_9216_4, 22437888, K10
	PRCENTRY2		xfft_r4dwpn_13824K_4608_4, 22364160, I7_64 + CORE2_64
	PRCENTRY2		xfft_r4delay_13824K_9216_4, 24342528, CORE2_64 + K10
	PRCENTRY2		xfft_r4delay_13824K_4608_4, 24244224, I7_64
	DD			0
	PRCSTRT	265600000, 14680064, 1.700
	PRCENTRY2		xfft_r4dwpn_14M_14_4, 49201664
	PRCENTRY2		xfft_r4dwpn_14M_13_4, 28327424, K10
	PRCENTRY2		xfft_r4dwpn_14M_12_4, 23615488, I7 + CORE2
	PRCENTRY2		xfft_r4delay_14M_14_4, 46853120
	PRCENTRY2		xfft_r4delay_14M_13_4, 25876480, K10
	PRCENTRY2		xfft_r4delay_14M_12_4, 25821184, I7 + CORE2
	PRCENTRY2		xfft_hg_14M_13_4, 18914944
	DD			0
	PRCSTRT	283900000, 15728640, 1.968
	PRCENTRY2		xfft_r4dwpn_15M_20480_4, 55582720
	PRCENTRY2		xfft_r4dwpn_15M_15360_4, 50561024, I7
	PRCENTRY2		xfft_r4dwpn_15M_12288_4, 30666752, CORE2_32 + K10
	PRCENTRY2		xfft_r4dwpn_15M_10240_4, 24940544
	PRCENTRY2		xfft_r4dwpn_15M_7680_4, 28024832
	PRCENTRY2		xfft_r4dwpn_15M_6144_4, 30562304
	PRCENTRY2		xfft_r4dwpn_15M_5120_4, 24862720, CORE2_64
	PRCENTRY2		xfft_r4dwpn_15M_3840_4, 23644160
	PRCENTRY2		xfft_r4delay_15M_20480_4, 51987456
	PRCENTRY2		xfft_r4delay_15M_15360_4, 48852992, I7
	PRCENTRY2		xfft_r4delay_15M_12288_4, 32425984, K10
	PRCENTRY2		xfft_r4delay_15M_10240_4, 27058176
	PRCENTRY2		xfft_r4delay_15M_7680_4, 26226688
	PRCENTRY2		xfft_r4delay_15M_6144_4, 32264192
	PRCENTRY2		xfft_r4delay_15M_5120_4, 26947584, CORE2
	PRCENTRY2		xfft_r4delay_15M_3840_4, 26212352
	DD			0
	PRCSTRT	295500000, 16384000, 2.035
	PRCENTRY2		xfft_r4dwpn_16000K_25600_4, 62057984, I7_64
	PRCENTRY2		xfft_r4dwpn_16000K_12800_4, 31887360
	PRCENTRY2		xfft_r4dwpn_16000K_6400_4, 31797248
	PRCENTRY2		xfft_r4delay_16000K_25600_4, 56740864, I7_64
	PRCENTRY2		xfft_r4delay_16000K_12800_4, 33720320
	PRCENTRY2		xfft_r4delay_16000K_6400_4, 33568768
	DD			0
	PRCSTRT	297800000, 16515072, 2.057
	PRCENTRY2		xfft_r4dwpn_16128K_9216_4, 31837696
	PRCENTRY2		xfft_r4dwpn_16128K_4608_4, 26515456
	PRCENTRY2		xfft_r4delay_16128K_9216_4, 29075456
	PRCENTRY2		xfft_r4delay_16128K_4608_4, 28991488
	DD			0
	PRCSTRT	302400000, 16777216, 2.100
	PRCENTRY2		xfft_r4dwpn_16M_14_4, 53927936, CORE2
	PRCENTRY2		xfft_r4dwpn_16M_13_4, 29913088
	PRCENTRY2		xfft_r4dwpn_16M_13_2, 29913088, K10
	PRCENTRY2		xfft_r4dwpn_16M_13_1, 29913088
	PRCENTRY2		xfft_r4dwpn_16M_12_4, 25219072, I7
	PRCENTRY2		xfft_r4dwpn_16M_12_2, 25219072
	PRCENTRY2		xfft_r4dwpn_16M_12_1, 25219072
	PRCENTRY2		xfft_r4delay_16M_14_4, 52105216, CORE2
	PRCENTRY2		xfft_r4delay_16M_13_4, 27992064
	PRCENTRY2		xfft_r4delay_16M_13_2, 27992064, K10
	PRCENTRY2		xfft_r4delay_16M_13_1, 27992064
	PRCENTRY2		xfft_r4delay_16M_12_4, 27955200, I7
	PRCENTRY2		xfft_r4delay_16M_12_2, 27955200
	PRCENTRY2		xfft_r4delay_16M_12_1, 27955200
	PRCENTRY2		xfft_hg_16M_13_4, 21544960
	DD			0
	PRCSTRT	329800000, 18350080, 2.325
	PRCENTRY2		xfft_r4dwpn_17920K_20480_4, 61489664
	PRCENTRY2		xfft_r4dwpn_17920K_10240_4, 35388928
	PRCENTRY2		xfft_r4dwpn_17920K_5120_4, 29472768, I7 + CORE2
	PRCENTRY2		xfft_r4delay_17920K_20480_4, 58551296
	PRCENTRY2		xfft_r4delay_17920K_10240_4, 32315392
	PRCENTRY2		xfft_r4delay_17920K_5120_4, 32219136, I7 + CORE2
	DD			0
	PRCSTRT	339200000, 18874368, 2.400
	PRCENTRY2		xfft_r4dwpn_18M_12288_4, 29900800, K10
	PRCENTRY2		xfft_r4dwpn_18M_9216_4, 33619968, CORE2
	PRCENTRY2		xfft_r4dwpn_18M_6144_4, 29806592
	PRCENTRY2		xfft_r4dwpn_18M_4608_4, 28315648, I7
	PRCENTRY2		xfft_r4dwpn_18M_4608_2, 28315648
	PRCENTRY2		xfft_r4delay_18M_12288_4, 32444416, K10
	PRCENTRY2		xfft_r4delay_18M_9216_4, 31453184, CORE2
	PRCENTRY2		xfft_r4delay_18M_6144_4, 32301056
	PRCENTRY2		xfft_r4delay_18M_4608_4, 31387648, I7
	PRCENTRY2		xfft_r4delay_18M_4608_2, 31387648
	DD			0
	PRCSTRT	352900000, 19660800, 2.513
	PRCENTRY2		xfft_r4dwpn_19200K_25600_4, 69435392, I7
	PRCENTRY2		xfft_r4dwpn_19200K_15360_4, 38309888, K10
	PRCENTRY2		xfft_r4dwpn_19200K_12800_4, 31088640
	PRCENTRY2		xfft_r4dwpn_19200K_7680_4, 38144000
	PRCENTRY2		xfft_r4dwpn_19200K_6400_4, 31008768, CORE2
	PRCENTRY2		xfft_r4delay_19200K_25600_4, 64939008, I7
	PRCENTRY2		xfft_r4delay_19200K_15360_4, 40511488, K10
	PRCENTRY2		xfft_r4delay_19200K_12800_4, 33738752
	PRCENTRY2		xfft_r4delay_19200K_7680_4, 40263680
	PRCENTRY2		xfft_r4delay_19200K_6400_4, 33605632, CORE2
	DD			0
	PRCSTRT	375900000, 20971520, 2.700
	PRCENTRY2		xfft_r4dwpn_20M_20480_4, 67395584
	PRCENTRY2		xfft_r4dwpn_20M_14_4, 40857600, CORE2_64
	PRCENTRY2		xfft_r4dwpn_20M_10240_4, 37367808
	PRCENTRY2		xfft_r4dwpn_20M_10240_2, 37367808, K10
	PRCENTRY2		xfft_r4dwpn_20M_13_4, 40704000, I7_32 + CORE2_32
	PRCENTRY2		xfft_r4dwpn_20M_5120_4, 31469568, I7_64
	PRCENTRY2		xfft_r4delay_20M_20480_4, 65114112
	PRCENTRY2		xfft_r4delay_20M_14_4, 43206656, CORE2_32
	PRCENTRY2		xfft_r4delay_20M_10240_4, 34955264
	PRCENTRY2		xfft_r4delay_20M_10240_2, 34955264, K10
	PRCENTRY2		xfft_r4delay_20M_13_4, 42962944, I7_32 + CORE2_64
	PRCENTRY2		xfft_r4delay_20M_5120_4, 34877440, I7_64
	PRCENTRY2		xfft_hg_20M_13_2, 26828928
	DD			0
	PRCSTRT	394700000, 22020096, 2.850
	PRCENTRY2		xfft_r4dwpn_21M_12288_4, 42446336, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_21M_6144_4, 35334144, I7
	PRCENTRY2		xfft_r4delay_21M_12288_4, 38750208, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_21M_6144_4, 38621184, I7 + CORE2_64
	DD			0
	PRCSTRT	410600000, 22937600, 2.982
	PRCENTRY2		xfft_r4dwpn_22400K_25600_4, 76816896, I7 + K10
	PRCENTRY2		xfft_r4dwpn_22400K_12800_4, 44158464
	PRCENTRY2		xfft_r4dwpn_22400K_6400_4, 36765696, CORE2
	PRCENTRY2		xfft_r4delay_22400K_25600_4, 73141248, I7 + K10
	PRCENTRY2		xfft_r4delay_22400K_12800_4, 40306688
	PRCENTRY2		xfft_r4delay_22400K_6400_4, 40187904, CORE2
	DD			0
	PRCSTRT	421100000, 23592960, 3.075
	PRCENTRY2		xfft_r4dwpn_23040K_15360_4, 37347328, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_23040K_9216_4, 45754368, I7_32 + CORE2
	PRCENTRY2		xfft_r4dwpn_23040K_7680_4, 37191680
	PRCENTRY2		xfft_r4delay_23040K_15360_4, 40529920, I7_64 + K10
	PRCENTRY2		xfft_r4delay_23040K_9216_4, 48291840, I7_32 + CORE2
	PRCENTRY2		xfft_r4delay_23040K_7680_4, 40300544
	DD			0
	PRCSTRT	449200000, 25165824, 3.300
	PRCENTRY2		xfft_r4dwpn_24M_14_4, 39829504
	PRCENTRY2		xfft_r4dwpn_24M_12288_4, 44818432, K10
	PRCENTRY2		xfft_r4dwpn_24M_13_4, 39686144, CORE2
	PRCENTRY2		xfft_r4dwpn_24M_6144_4, 37724160, I7
	PRCENTRY2		xfft_r4delay_24M_14_4, 43225088
	PRCENTRY2		xfft_r4delay_24M_12288_4, 41914368, K10
	PRCENTRY2		xfft_r4delay_24M_13_4, 42999808, CORE2
	PRCENTRY2		xfft_r4delay_24M_6144_4, 41803776, I7
	PRCENTRY2		xfft_hg_24M_13_2, 32104704
	DD			0
	PRCSTRT	467500000, 26214400, 3.450
	PRCENTRY2		xfft_r4dwpn_25M_25600_4, 84197376, I7 + K10
	PRCENTRY2		xfft_r4dwpn_25M_20480_4, 51048448
	PRCENTRY2		xfft_r4dwpn_25M_12800_4, 46628864
	PRCENTRY2		xfft_r4dwpn_25M_10240_4, 50845696, CORE2
	PRCENTRY2		xfft_r4dwpn_25M_6400_4, 39254016
	PRCENTRY2		xfft_r4delay_25M_25600_4, 81342464, I7 + K10
	PRCENTRY2		xfft_r4delay_25M_20480_4, 53987328, CORE2_32
	PRCENTRY2		xfft_r4delay_25M_12800_4, 43601920
	PRCENTRY2		xfft_r4delay_25M_10240_4, 53661696, CORE2_64
	PRCENTRY2		xfft_r4delay_25M_6400_4, 43501568
	DD			0
	PRCSTRT	490100000, 27525120, 3.638
	PRCENTRY2		xfft_r4dwpn_26880K_15360_4, 53038592, K10
	PRCENTRY2		xfft_r4dwpn_26880K_7680_4, 44095488, I7 + CORE2
	PRCENTRY2		xfft_r4delay_26880K_15360_4, 48408576, CORE2_32 + K10
	PRCENTRY2		xfft_r4delay_26880K_7680_4, 48193536, I7 + CORE2_64
	DD			0
	PRCSTRT	503100000, 28311552, 3.750
	PRCENTRY2		xfft_r4dwpn_27M_9216_4, 44605440, I7_64 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_27M_9216_4, 48328704, I7_64 + CORE2 + K10
	DD			0
	PRCSTRT	522800000, 29360128, 3.900
	PRCENTRY2		xfft_r4dwpn_28M_14_4, 56569344, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_28M_13_4, 47048704, I7
	PRCENTRY2		xfft_r4delay_28M_14_4, 51628032, CORE2 + K10
	PRCENTRY2		xfft_r4delay_28M_13_4, 51417088, I7
	PRCENTRY2		xfft_hg_28M_13_2, 37380480
	DD			0
	PRCSTRT	558300000, 31457280, 4.200
	PRCENTRY2		xfft_r4dwpn_30M_20480_4, 49758208
	PRCENTRY2		xfft_r4dwpn_30M_15360_4, 56000512, K10_64
	PRCENTRY2		xfft_r4dwpn_30M_12288_4, 60983296, I7_32 + CORE2 + K10_32
	PRCENTRY2		xfft_r4dwpn_30M_10240_4, 49565696
	PRCENTRY2		xfft_r4dwpn_30M_7680_4, 47075328, I7_64
	PRCENTRY2		xfft_r4delay_30M_20480_4, 54005760
	PRCENTRY2		xfft_r4delay_30M_15360_4, 52359168
	PRCENTRY2		xfft_r4delay_30M_12288_4, 64356352, I7_32 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_30M_10240_4, 53698560
	PRCENTRY2		xfft_r4delay_30M_7680_4, 52162560, I7_64
	DD			0
	PRCSTRT	581500000, 32768000, 4.388
	PRCENTRY2		xfft_r4dwpn_32000K_25600_4, 63754240, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_32000K_12800_4, 63465472
	PRCENTRY2		xfft_r4delay_32000K_25600_4, 67430400, I7_64 + K10
	PRCENTRY2		xfft_r4delay_32000K_12800_4, 66977792
	DD			0
	PRCSTRT	585300000, 33030144, 4.425
	PRCENTRY2		xfft_r4dwpn_32256K_9216_4, 52885504
	PRCENTRY2		xfft_r4delay_32256K_9216_4, 57794560
	DD			0
	PRCSTRT	595800000, 33554432, 4.500
	PRCENTRY2		xfft_r4dwpn_32M_14_4, 59727872, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_32M_13_4, 50225152, I7
	PRCENTRY2		xfft_r4dwpn_32M_13_2, 50225152
	PRCENTRY2		xfft_r4dwpn_32M_13_1, 50225152
	PRCENTRY2		xfft_r4delay_32M_14_4, 55840768, CORE2 + K10
	PRCENTRY2		xfft_r4delay_32M_13_4, 55648256, I7
	PRCENTRY2		xfft_r4delay_32M_13_2, 55648256
	PRCENTRY2		xfft_r4delay_32M_13_1, 55648256
	PRCENTRY2		xfft_hg_32M_13_2, 42639360
	DD			0
	DD	0

xjmptablep DD	0
	org	$-4
	PRCSTRT	739,	32,	0.00000111
	PRCENTRY		xfft_hg_32_op_ac_BLEND, 1152
	DD			4, 4
	DD			2, 1, 1, 1, 1, 0
	PRCSTRT	1095,	48,	0.00000144
	PRCENTRY		xfft_hg_48_op_ac_BLEND, 1920
	DD			6, 6
	DD			3, 1, 1, 1, 1, 0
	PRCSTRT	1465,	64,	0.00000178
	PRCENTRY		xfft_hg_64_op_ac_BLEND, 2432
	DD			8, 8
	DD			4, 1, 1, 1, 1, 0
	PRCSTRT	2173,	96,	0.00000259
	PRCENTRY		xfft_hg_96_op_ac_BLEND, 3328
	DD			12, 12
	DD			6, 3, 1, 1, 1, 0
	PRCSTRT	2897,	128,	0.00000319
	PRCENTRY		xfft_hg_128_op_ac_BLEND, 4352
	DD			16, 16
	DD			8, 4, 1, 2, 1, 0
	PRCSTRT	4295,	192,	0.00000542
	PRCENTRY		xfft_hg_192_op_ac_BLEND, 6784
	DD			24, 24
	DD			12, 6, 3, 3, 1, 0
	PRCSTRT	5729,	256,	0.00000691
	PRCENTRY		xfft_hg_256_op_ac_BLEND, 8576
	DD			32, 32
	DD			16, 8, 4, 4, 1, 0
	PRCSTRT	8493,	384,	0.0000111
	PRCENTRY		xfft_hg_384_op_ac_BLEND, 13312
	DD			48, 48
	DD			24, 12, 3, 6, 1, 0
	PRCSTRT	11319,	512,	0.0000143
	PRCENTRY		xfft_hg_512_op_ac_BLEND, 17152
	DD			64, 64
	DD			32, 16, 4, 2*256+8, 1, 0
	PRCSTRT	16779,	768,	0.0000260
	PRCENTRY		xfft_hg_768_op_ac_BLEND, 26624
	DD			96, 96
	DD			48, 24, 3*256+6, 3*256+12, 1, 0
	PRCSTRT	22381,	1024,	0.0000349
	PRCENTRY		xfft_hg_1024_op_ac_BLEND, 34304
	DD			128, 128
	DD			64, 32, 4*256+8, 4*256+16, 1, 0
	PRCSTRT	33189,	1536,	0.0000601
	PRCENTRY		xfft_hg_1536_op_ac_BLEND, 52992
	DD			192, 192
	DD			96, 48, 12, 6*256+24, 3, 0
	PRCSTRT	44221,	2048,	0.0000773
	PRCENTRY		xfft_hg_2048_op_ac_BLEND, 68608
	DD			256, 256
	DD			128, 64, 16, 8*256+32, 2*256+4, 0
	PRCSTRT	65519,	3072,	0.000131
	PRCENTRY		xfft_hg_3072_op_ac_BLEND, 106112
	DD			384, 384
	DD			192, 96, 6*256+24, 12*256+48, 3*256+3, 0
	PRCSTRT	87271,	4096,	0.000172
	PRCENTRY		xfft_hg_4096_op_ac_BLEND, 137216
	DD			512, 512
	DD			256, 128, 8*256+32, 16*256+64, 4*256+4, 0
	PRCSTRT	129600,	6144,	0.000291
	PRCENTRY		xfft_hg_6144_op_ac_BLEND, 211968
	DD			768, 768
	DD			384, 192, 12*256+48, 24*256+96, 6*256+3, 0
	PRCSTRT	172800,	8192,	0.000395
	PRCENTRY2		xfft_r4_8K_ac_8_4, 87552, P4TP_256
	PRCENTRY2		xfft_r4_8K_ac_np_8_4, 87552, P4_1024 + P4TP_512 + I7 + CORE2 + CORE2_512_64 + K8 + K10
	PRCENTRY2A		xfft_r4_8K_ac_np_8_4_P4, 87552, CORE2_512_32
IFDEF IMPL_ALL_CORE
	PRCENTRY2A		xfft_hg_8192_op_ac_BLEND, 274432, CORE2	;; Only used when testing all possible FFT implementations
	DD			1024, 1024
	DD			512, 256, 16*256+64, 32*256+128, 8*256+4, 0
ELSE
	DD			0
ENDIF
	PRCSTRT	256200,	12288,	0.000626
	PRCENTRY2		xfft_hg_12K_ac_ip_8_4, 64896, P4_1024 + P4TP_512 + P4TP_256 + CORE2_64 + K8 + K10
	PRCENTRY2A		xfft_hg_12K_ac_ip_8_4_P4, 64896, CORE2_32 + CORE2_512 + I7
	DD			0
	PRCSTRT	340800,	16384,	0.000857
	PRCENTRY2		xfft_r4_16K_ac_8_4, 183296
	PRCENTRY2A		xfft_r4_16K_ac_8_4_P4, 183296, CORE2_512_32
	PRCENTRY2		xfft_r4_16K_ac_np_8_4, 183296, P4TP_512 + P4TP_256 + I7 + CORE2 + CORE2_512_64 + K8 + K10
	PRCENTRY2A		xfft_r4_16K_ac_np_8_4_CORE, 183296, P4_1024
	PRCENTRY2		xfft_hg_16K_ac_ip_8_4, 84224
	DD			0
	PRCSTRT	505700,	24576,	0.00135
	PRCENTRY2		xfft_r4_24K_ac_768_4, 259584
	PRCENTRY2A		xfft_r4_24K_ac_768_4_P4, 259584, CORE2_512_32
	PRCENTRY2		xfft_r4_24K_ac_np_768_4, 259584, I7_64 + CORE2 + K8_64 + K10_64
	PRCENTRY2		xfft_r4_24K_ac_8_4, 246272, P4TP_256
	PRCENTRY2A		xfft_r4_24K_ac_8_4_P4, 246272, I7_32
	PRCENTRY2		xfft_r4_24K_ac_np_8_4, 246272, P4TP_512 + CORE2_512_64 + K8_32 + K10_32
	PRCENTRY2		xfft_hg_24K_ac_ip_8_4, 123904, P4_1024
	DD			0
	PRCSTRT	672300,	32768,	0.00179
	PRCENTRY2		xfft_r4dwpn_32K_ac_8_4, 153600, P4TP_256
	PRCENTRY2A		xfft_r4dwpn_32K_ac_8_4_P4, 153600, I7
	PRCENTRY2		xfft_r4dwpn_32K_ac_np_8_4, 153600, P4TP_512 + CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_32K_ac_8_2, 153600, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_32K_ac_np_8_2, 153600
	PRCENTRY2A		xfft_r4dwpn_32K_ac_np_8_2_P4, 153600, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_32K_ac_8_1, 153600
	PRCENTRY2		xfft_r4dwpn_32K_ac_np_8_1, 153600
	PRCENTRY2		xfft_r4delay_32K_ac_8_4, 132096
	PRCENTRY2		xfft_r4delay_32K_ac_np_8_4, 132096
	PRCENTRY2		xfft_r4delay_32K_ac_8_2, 132096
	PRCENTRY2		xfft_r4delay_32K_ac_np_8_2, 132096
	PRCENTRY2		xfft_r4delay_32K_ac_8_1, 132096
	PRCENTRY2		xfft_r4delay_32K_ac_np_8_1, 132096
	PRCENTRY2		xfft_r4_32K_ac_10_4, 349696
	PRCENTRY2		xfft_r4_32K_ac_np_10_4, 349696
	PRCENTRY2		xfft_r4_32K_ac_8_4, 333824
	PRCENTRY2		xfft_r4_32K_ac_np_8_4, 333824
	PRCENTRY2		xfft_hg_32K_ac_ip_8_4, 162816, P4_1024
	DD			0
	PRCSTRT	836000, 40960,	0.00236
	PRCENTRY2		xfft_r4_40K_ac_1280_4, 431616, P4TP_512 + P4TP_256 + CORE2_512 + K8 + K10
	PRCENTRY2		xfft_r4_40K_ac_np_1280_4, 431616, P4_1024 + I7 + CORE2
	DD			0
	PRCSTRT	999000, 49152,	0.00303
	PRCENTRY2		xfft_r4_48K_ac_768_4, 543744, CORE2_512_64 + K8 + K10
	PRCENTRY2A		xfft_r4_48K_ac_768_4_P4, 543744, CORE2_512_32
	PRCENTRY2		xfft_r4_48K_ac_np_768_4, 543744, P4_1024_64 + I7 + CORE2
	PRCENTRY2A		xfft_r4_48K_ac_np_768_4_CORE, 543744, P4_1024_32
	PRCENTRY2		xfft_hg_48K_ac_ip_11_4, 145280, P4TP_512 + P4TP_256
	DD			0
	PRCSTRT	1327000, 65536,	0.00404
	PRCENTRY2		xfft_r4dwpn_64K_ac_8_4, 297984, P4_1024_32 + P4TP_512 + P4TP_256
	PRCENTRY2A		xfft_r4dwpn_64K_ac_8_4_P4, 297984, I7
	PRCENTRY2A		xfft_r4dwpn_64K_ac_8_4_CORE, 297984, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_64K_ac_np_8_4, 297984, CORE2
	PRCENTRY2		xfft_r4dwpn_64K_ac_8_2, 297984, K10
	PRCENTRY2		xfft_r4dwpn_64K_ac_np_8_2, 297984
	PRCENTRY2		xfft_r4dwpn_64K_ac_8_1, 297984, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_64K_ac_np_8_1, 297984
	PRCENTRY2		xfft_r4delay_64K_ac_8_4, 256000
	PRCENTRY2		xfft_r4delay_64K_ac_np_8_4, 256000
	PRCENTRY2		xfft_r4delay_64K_ac_8_2, 256000
	PRCENTRY2		xfft_r4delay_64K_ac_np_8_2, 256000
	PRCENTRY2		xfft_r4delay_64K_ac_8_1, 256000
	PRCENTRY2		xfft_r4delay_64K_ac_np_8_1, 256000
	PRCENTRY2		xfft_r4_64K_ac_10_4, 728064
	PRCENTRY2		xfft_r4_64K_ac_np_10_4, 728064
	PRCENTRY2		xfft_r4_64K_ac_8_4, 675840
	PRCENTRY2		xfft_r4_64K_ac_np_8_4, 675840
	PRCENTRY2		xfft_hg_64K_ac_ip_11_4, 178560
	DD			0
	PRCSTRT	1482000, 73728,	0.00463
	PRCENTRY2		xfft_r4_72K_ac_768_4, 729600, P4_1024_32 + P4TP_512 + P4TP_256 + K10
	PRCENTRY2A		xfft_r4_72K_ac_768_4_CORE, 729600, P4_1024_64
	PRCENTRY2		xfft_r4_72K_ac_np_768_4, 729600, I7 + CORE2
	PRCENTRY2		xfft_r4_72K_ac_768_2, 729600, CORE2_512 + K8
	DD			0
	PRCSTRT	1650000, 81920,	0.00521
	PRCENTRY2		xfft_r4_80K_ac_1280_4, 904192, P4_1024 + P4TP_512 + P4TP_256 + K10
	PRCENTRY2		xfft_r4_80K_ac_np_1280_4, 904192, I7 + CORE2
	PRCENTRY2		xfft_r4_80K_ac_1280_2, 904192, CORE2_512 + K8
	DD			0
	PRCSTRT	1970000, 98304, 0.00644
	PRCENTRY2		xfft_r4dwpn_96K_ac_768_4, 440320, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_96K_ac_np_768_4, 440320, I7_32 + CORE2
	PRCENTRY2		xfft_r4dwpn_96K_ac_768_2, 440320, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_96K_ac_8_4, 256000
	PRCENTRY2A		xfft_r4dwpn_96K_ac_8_4_CORE, 256000, P4_1024
	PRCENTRY2		xfft_r4dwpn_96K_ac_np_8_4, 256000
	PRCENTRY2		xfft_r4delay_96K_ac_768_4, 377856, P4TP_256
	PRCENTRY2		xfft_r4delay_96K_ac_np_768_4, 377856
	PRCENTRY2		xfft_r4delay_96K_ac_8_4, 248832
	PRCENTRY2		xfft_r4delay_96K_ac_np_8_4, 248832
	PRCENTRY2		xfft_r4_96K_ac_10_4, 975360
	PRCENTRY2		xfft_r4_96K_ac_np_10_4, 975360
	PRCENTRY2		xfft_r4_96K_ac_768_4, 989184
	PRCENTRY2		xfft_r4_96K_ac_np_768_4, 989184
	PRCENTRY2		xfft_hg_96K_ac_ip_11_4, 245632, P4TP_512
	DD			0
	PRCSTRT	2448000, 122880, 0.00815
	PRCENTRY2		xfft_r4_120K_ac_1280_4, 1212928
	PRCENTRY2		xfft_r4_120K_ac_np_1280_4, 1212928
	DD			0
	PRCSTRT	2618000, 131072, 0.00871
	PRCENTRY2		xfft_r4dwpn_128K_ac_10_4, 587776, P4_1024_64 + I7_64 + K10
	PRCENTRY2A		xfft_r4dwpn_128K_ac_10_4_P4, 587776, I7_32
	PRCENTRY2A		xfft_r4dwpn_128K_ac_10_4_CORE, 587776, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_128K_ac_np_10_4, 587776, CORE2
	PRCENTRY2		xfft_r4dwpn_128K_ac_10_2, 587776, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_128K_ac_8_4, 339968, P4TP_512
	PRCENTRY2		xfft_r4dwpn_128K_ac_8_2, 339968
	PRCENTRY2		xfft_r4dwpn_128K_ac_np_8_4, 339968
	PRCENTRY2		xfft_r4delay_128K_ac_10_4, 504832, P4TP_256
	PRCENTRY2		xfft_r4delay_128K_ac_np_10_4, 504832
	PRCENTRY2		xfft_r4delay_128K_ac_8_4, 339968
	PRCENTRY2		xfft_r4delay_128K_ac_np_8_4, 339968
	PRCENTRY2		xfft_hg_128K_ac_ip_11_4, 312064
	DD			0
	PRCSTRT	2930000, 147456, 0.00994
	PRCENTRY2		xfft_r4_144K_ac_1536_4, 1454592, P4TP_512 + P4TP_256 + I7_64
	PRCENTRY2A		xfft_r4_144K_ac_1536_4_P4, 1454592, CORE2_32 + I7_32
	PRCENTRY2A		xfft_r4_144K_ac_1536_4_CORE, 1454592, P4_1024
	PRCENTRY2		xfft_r4_144K_ac_1536_2, 1454592, K8 + K10_32
	DD			0
	PRCSTRT	3255000, 163840, 0.0112
	PRCENTRY2		xfft_r4dwpn_160K_ac_1280_4, 727040, P4TP_256 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4dwpn_160K_ac_1280_4_P4, 727040, CORE2_32 + I7
	PRCENTRY2A		xfft_r4dwpn_160K_ac_1280_4_CORE, 727040, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_160K_ac_1280_2, 727040, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_160K_ac_8_4, 358400
	PRCENTRY2A		xfft_r4dwpn_160K_ac_8_4_CORE, 358400, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_160K_ac_8_2, 358400
	PRCENTRY2		xfft_r4delay_160K_ac_1280_4, 623616
	PRCENTRY2		xfft_r4delay_160K_ac_8_4, 365568, P4TP_512
	DD			0
	PRCSTRT	3891000, 196608, 0.0136
	PRCENTRY2		xfft_r4dwpn_192K_ac_768_4, 863232, P4TP_256 + I7_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4dwpn_192K_ac_768_4_P4, 863232, CORE2_32 + I7_32
	PRCENTRY2		xfft_r4dwpn_192K_ac_768_2, 863232, K8_64
	PRCENTRY2		xfft_r4dwpn_192K_ac_768_1, 863232, CORE2_512 + K8_32
	PRCENTRY2		xfft_r4dwpn_192K_ac_8_4, 492544
	PRCENTRY2A		xfft_r4dwpn_192K_ac_8_4_CORE, 492544, P4_1024
	PRCENTRY2		xfft_r4delay_192K_ac_768_4, 731136
	PRCENTRY2		xfft_r4delay_192K_ac_8_4, 489472
	PRCENTRY2		xfft_r4_192K_ac_768_4, 2002944
	PRCENTRY2		xfft_hg_192K_ac_10_4, 529536
	PRCENTRY2		xfft_hg_192K_ac_10_2, 529536
	PRCENTRY2		xfft_hg_192K_ac_10_1, 529536, P4TP_512
	DD			0
	PRCSTRT	5167000, 262144, 0.0179
	PRCENTRY2		xfft_r4dwpn_256K_ac_10_4, 1149952, P4_1024_32 + P4TP_512 + P4TP_256 + CORE2_64 + K10_64
	PRCENTRY2A		xfft_r4dwpn_256K_ac_10_4_P4, 1149952, I7
	PRCENTRY2		xfft_r4dwpn_256K_ac_10_2, 1149952, CORE2_512 + K8_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_256K_ac_10_1, 1149952, K8_32
	PRCENTRY2		xfft_r4dwpn_256K_ac_8_4, 655360
	PRCENTRY2A		xfft_r4dwpn_256K_ac_8_4_P4, 655360, CORE2_32
	PRCENTRY2A		xfft_r4dwpn_256K_ac_8_4_CORE, 655360, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_256K_ac_8_2, 655360
	PRCENTRY2		xfft_r4delay_256K_ac_10_4, 972800
	PRCENTRY2		xfft_r4delay_256K_ac_8_4, 671744
	PRCENTRY2		xfft_hg_256K_ac_10_4, 697344
	PRCENTRY2		xfft_hg_256K_ac_10_2, 697344
	PRCENTRY2		xfft_hg_256K_ac_10_1, 697344
	DD			0
	PRCSTRT	5773000, 294912, 0.0203
	PRCENTRY2		xfft_r4dwpn_288K_ac_768_4, 706560, P4_1024_64 + CORE2_64 + K10
	PRCENTRY2A		xfft_r4dwpn_288K_ac_768_4_P4, 706560, CORE2_32 + I7_32
	PRCENTRY2A		xfft_r4dwpn_288K_ac_768_4_CORE, 706560, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_288K_ac_768_2, 706560, P4TP_256 + I7_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_288K_ac_768_1, 706560, K8_32
	PRCENTRY2		xfft_r4delay_288K_ac_768_4, 691200, P4TP_512
	DD			0
	PRCSTRT	6432000, 327680, 0.0226
	PRCENTRY2		xfft_r4dwpn_320K_ac_1280_4, 1428480, P4TP_512 + P4TP_256 + CORE2_64 + K10_64
	PRCENTRY2A		xfft_r4dwpn_320K_ac_1280_4_P4, 1428480, CORE2_32 + I7
	PRCENTRY2A		xfft_r4dwpn_320K_ac_1280_4_CORE, 1428480, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_320K_ac_1280_2, 1428480, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_320K_ac_1280_1, 1428480, K8_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_320K_ac_8_4, 687104
	PRCENTRY2A		xfft_r4dwpn_320K_ac_8_4_CORE, 687104, P4_1024_32
	PRCENTRY2		xfft_r4delay_320K_ac_1280_4, 1206272
	PRCENTRY2		xfft_r4delay_320K_ac_8_4, 722944
	DD			0
	PRCSTRT	7689000, 393216, 0.0273
	PRCENTRY2		xfft_r4dwpn_384K_ac_1536_4, 1711104, P4TP_256 + I7_64
	PRCENTRY2A		xfft_r4dwpn_384K_ac_1536_4_P4, 1711104, I7_32
	PRCENTRY2A		xfft_r4dwpn_384K_ac_1536_4_CORE, 1711104, P4_1024
	PRCENTRY2		xfft_r4dwpn_384K_ac_10_4, 935936, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_384K_ac_10_2, 935936, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_384K_ac_10_1, 935936, K8_32
	PRCENTRY2A		xfft_r4dwpn_384K_ac_10_1_P4, 935936, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_384K_ac_768_4, 937984
	PRCENTRY2		xfft_r4dwpn_384K_ac_768_2, 937984
	PRCENTRY2		xfft_r4dwpn_384K_ac_768_1, 937984
	PRCENTRY2		xfft_r4delay_384K_ac_1536_4, 1443840
	PRCENTRY2		xfft_r4delay_384K_ac_10_4, 916480
	PRCENTRY2		xfft_r4delay_384K_ac_768_4, 946176
	PRCENTRY2		xfft_hg_384K_ac_10_4, 1036288, P4TP_512
	PRCENTRY2		xfft_hg_384K_ac_10_2, 1036288
	PRCENTRY2		xfft_hg_384K_ac_10_1, 1036288
	DD			0
	PRCSTRT	9545000, 491520, 0.0345
	PRCENTRY2		xfft_r4dwpn_480K_ac_1280_4, 1157120, I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_480K_ac_1280_2, 1157120
	PRCENTRY2		xfft_r4dwpn_480K_ac_1280_1, 1157120, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_480K_ac_768_4, 972800
	PRCENTRY2		xfft_r4delay_480K_ac_1280_4, 1133568
	PRCENTRY2		xfft_r4delay_480K_ac_768_4, 1004544
	DD			0
	PRCSTRT	10190000, 524288, 0.0368
	PRCENTRY2		xfft_r4dwpn_512K_ac_11_4, 2284544, P4_1024_32 + P4TP_512 + P4TP_256 + CORE2_512_64
	PRCENTRY2A		xfft_r4dwpn_512K_ac_11_4_P4, 2284544, I7
	PRCENTRY2A		xfft_r4dwpn_512K_ac_11_4_CORE, 2284544, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_512K_ac_10_4, 1241088, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_512K_ac_10_2, 1241088, K8
	PRCENTRY2		xfft_r4dwpn_512K_ac_10_1, 1241088
	PRCENTRY2A		xfft_r4dwpn_512K_ac_10_1_P4, 1241088, CORE2_512_32
	PRCENTRY2		xfft_r4delay_512K_ac_11_4, 1927168
	PRCENTRY2		xfft_r4delay_512K_ac_10_4, 1253376
	PRCENTRY2		xfft_hg_512K_ac_10_4, 1372160
	PRCENTRY2		xfft_hg_512K_ac_10_2, 1372160
	PRCENTRY2		xfft_hg_512K_ac_10_1, 1372160
	DD			0
	PRCSTRT	11400000, 589824, 0.0418
	PRCENTRY2		xfft_r4dwpn_576K_ac_2304_4, 2550784, P4_1024 + I7 + K8_32
	PRCENTRY2		xfft_r4dwpn_576K_ac_2304_2, 2550784
	PRCENTRY2		xfft_r4dwpn_576K_ac_2304_1, 2550784
	PRCENTRY2		xfft_r4dwpn_576K_ac_1536_4, 1382400, P4TP_256 + CORE2_512 + K10
	PRCENTRY2		xfft_r4dwpn_576K_ac_768_4, 1385472, CORE2
	PRCENTRY2		xfft_r4dwpn_576K_ac_768_2, 1385472, K8_64
	PRCENTRY2		xfft_r4dwpn_576K_ac_768_1, 1385472
	PRCENTRY2		xfft_r4delay_576K_ac_2304_4, 2148352
	PRCENTRY2		xfft_r4delay_576K_ac_1536_4, 1354752, P4TP_512
	PRCENTRY2		xfft_r4delay_576K_ac_768_4, 1357824
	DD			0
	PRCSTRT	12660000, 655360, 0.0467
	PRCENTRY2		xfft_r4dwpn_640K_ac_2560_4, 2841600, P4TP_512 + P4TP_256 + CORE2_512_64
	PRCENTRY2A		xfft_r4dwpn_640K_ac_2560_4_P4, 2841600, I7
	PRCENTRY2A		xfft_r4dwpn_640K_ac_2560_4_CORE, 2841600, P4_1024
	PRCENTRY2		xfft_r4dwpn_640K_ac_1280_4, 1536000, CORE2 + K8_64 + K10
	PRCENTRY2		xfft_r4dwpn_640K_ac_1280_2, 1536000, K8_32
	PRCENTRY2		xfft_r4dwpn_640K_ac_1280_1, 1536000, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_640K_ac_10_4, 1284096
	PRCENTRY2		xfft_r4delay_640K_ac_2560_4, 2394112
	PRCENTRY2		xfft_r4delay_640K_ac_1280_4, 1552384
	PRCENTRY2		xfft_r4delay_640K_ac_1280_2, 1552384
	PRCENTRY2		xfft_r4delay_640K_ac_10_4, 1328128
	DD			0
	PRCSTRT	15130000, 786432, 0.0566
	PRCENTRY2		xfft_r4dwpn_768K_ac_3072_4, 3410944, P4_1024 + I7_64
	PRCENTRY2		xfft_r4dwpn_768K_ac_11_4, 1841152, P4TP_512
	PRCENTRY2		xfft_r4dwpn_768K_ac_11_2, 1841152, P4TP_256 + CORE2_512
	PRCENTRY2		xfft_r4dwpn_768K_ac_1536_4, 1835008, K10_32
	PRCENTRY2		xfft_r4dwpn_768K_ac_1536_2, 1835008
	PRCENTRY2		xfft_r4dwpn_768K_ac_10_4, 1836032, CORE2 + K10_64
	PRCENTRY2A		xfft_r4dwpn_768K_ac_10_4_P4, 1836032, I7_32
	PRCENTRY2		xfft_r4dwpn_768K_ac_10_2, 1836032
	PRCENTRY2		xfft_r4dwpn_768K_ac_10_1, 1836032, K8
	PRCENTRY2		xfft_r4dwpn_768K_ac_768_4, 1843200
	PRCENTRY2		xfft_r4delay_768K_ac_3072_4, 2873344
	PRCENTRY2		xfft_r4delay_768K_ac_11_4, 1805312
	PRCENTRY2		xfft_r4delay_768K_ac_11_2, 1805312
	PRCENTRY2		xfft_r4delay_768K_ac_1536_4, 1855488
	PRCENTRY2		xfft_r4delay_768K_ac_1536_2, 1855488
	PRCENTRY2		xfft_r4delay_768K_ac_10_4, 1796096
	PRCENTRY2		xfft_r4delay_768K_ac_768_4, 1867776
	PRCENTRY2		xfft_hg_768K_ac_12_2, 1285248
	PRCENTRY2		xfft_hg_768K_ac_11_4, 1648640
	PRCENTRY2		xfft_hg_768K_ac_11_2, 1648640
	DD			0
	PRCSTRT	15740000, 819200, 0.0591
	PRCENTRY2		xfft_r4dwpn_800K_ac_1280_4, 1587200, P4_1024 + CORE2 + K10
	PRCENTRY2A		xfft_r4dwpn_800K_ac_1280_4_P4, 1587200, I7
	PRCENTRY2		xfft_r4dwpn_800K_ac_1280_2, 1587200, P4TP_512 + P4TP_256 + K8_64
	PRCENTRY2		xfft_r4dwpn_800K_ac_1280_1, 1587200, CORE2_512 + K8_32
	PRCENTRY2		xfft_r4delay_800K_ac_1280_4, 1643520
	PRCENTRY2		xfft_r4delay_800K_ac_1280_2, 1643520
	DD			0
	PRCSTRT	16890000, 884736, 0.0639
	PRCENTRY2		xfft_r4dwpn_864K_ac_2304_4, 2050048, P4_1024 + I7 + CORE2 + K8 + K10
	PRCENTRY2		xfft_r4dwpn_864K_ac_2304_2, 2050048, P4TP_256
	PRCENTRY2		xfft_r4dwpn_864K_ac_2304_1, 2050048
	PRCENTRY2		xfft_r4delay_864K_ac_2304_4, 2010112
	PRCENTRY2		xfft_r4delay_864K_ac_2304_2, 2010112
	DD			0
	PRCSTRT	18820000, 983040, 0.0713
	PRCENTRY2		xfft_r4dwpn_960K_ac_3840_4, 4246528, P4_1024_64 + I7
	PRCENTRY2		xfft_r4dwpn_960K_ac_3840_2, 4246528, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_960K_ac_3840_1, 4246528, K8_32
	PRCENTRY2		xfft_r4dwpn_960K_ac_2560_4, 2283520
	PRCENTRY2		xfft_r4dwpn_960K_ac_2560_2, 2283520, P4TP_256
	PRCENTRY2		xfft_r4dwpn_960K_ac_1536_4, 1894400, K10_32
	PRCENTRY2A		xfft_r4dwpn_960K_ac_1536_4_CORE, 1894400, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_960K_ac_1536_2, 1894400
	PRCENTRY2		xfft_r4dwpn_960K_ac_1280_4, 2278400, CORE2 + K10_64
	PRCENTRY2		xfft_r4dwpn_960K_ac_1280_2, 2278400
	PRCENTRY2		xfft_r4dwpn_960K_ac_1280_1, 2278400
	PRCENTRY2		xfft_r4dwpn_960K_ac_768_4, 1907712
	PRCENTRY2		xfft_r4delay_960K_ac_3840_4, 3573760
	PRCENTRY2		xfft_r4delay_960K_ac_2560_4, 2239488
	PRCENTRY2		xfft_r4delay_960K_ac_2560_2, 2239488
	PRCENTRY2		xfft_r4delay_960K_ac_1536_4, 1963008
	PRCENTRY2		xfft_r4delay_960K_ac_1536_2, 1963008
	PRCENTRY2		xfft_r4delay_960K_ac_1280_4, 2226176
	PRCENTRY2		xfft_r4delay_960K_ac_768_4, 1984512
	DD			0
	PRCSTRT	20080000, 1048576, 0.0761
	PRCENTRY2		xfft_r4dwpn_1M_ac_12_4, 4541440, P4_1024 + P4TP_512 + I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_1M_ac_12_2, 4541440, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_1M_ac_12_1, 4541440, K8_32
	PRCENTRY2		xfft_r4dwpn_1M_ac_11_4, 2441216, K10_32
	PRCENTRY2		xfft_r4dwpn_1M_ac_11_2, 2441216, P4TP_256
	PRCENTRY2		xfft_r4dwpn_1M_ac_10_4, 2441216, CORE2_32
	PRCENTRY2		xfft_r4delay_1M_ac_12_4, 3823616
	PRCENTRY2		xfft_r4delay_1M_ac_11_4, 2469888
	PRCENTRY2		xfft_r4delay_1M_ac_11_2, 2469888
	PRCENTRY2		xfft_r4delay_1M_ac_10_4, 2469888
	PRCENTRY2		xfft_hg_1024K_ac_12_4, 1682432
	PRCENTRY2		xfft_hg_1024K_ac_12_1, 1682432
	PRCENTRY2		xfft_hg_1024K_ac_11_4, 2181120
	PRCENTRY2		xfft_hg_1024K_ac_11_2, 2181120
	DD			0
	PRCSTRT	22510000, 1179648, 0.088
	PRCENTRY2		xfft_r4dwpn_1152K_ac_4608_4, 5086208, P4_1024 + P4TP_512 + I7
	PRCENTRY2		xfft_r4dwpn_1152K_ac_3072_4, 2738176, CORE2_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_1152K_ac_3072_2, 2738176, P4TP_256 + CORE2_32 + CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_1152K_ac_3072_1, 2738176, K8_32
	PRCENTRY2		xfft_r4dwpn_1152K_ac_2304_4, 2723840
	PRCENTRY2		xfft_r4dwpn_1152K_ac_1536_4, 2724864, K10_64
	PRCENTRY2		xfft_r4delay_1152K_ac_4608_4, 4278272
	PRCENTRY2		xfft_r4delay_1152K_ac_3072_4, 2685952
	PRCENTRY2		xfft_r4delay_1152K_ac_3072_2, 2685952
	PRCENTRY2		xfft_r4delay_1152K_ac_2304_4, 2756608
	PRCENTRY2		xfft_r4delay_1152K_ac_1536_4, 2660352
	DD			0
	PRCSTRT	24980000, 1310720, 0.101
	PRCENTRY2		xfft_r4dwpn_1280K_ac_5120_4, 5671936, P4_1024 + I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_1280K_ac_5120_2, 5671936, CORE2_512
	PRCENTRY2		xfft_r4dwpn_1280K_ac_5120_1, 5671936, K8
	PRCENTRY2		xfft_r4dwpn_1280K_ac_2560_4, 3031040, K10_32
	PRCENTRY2		xfft_r4dwpn_1280K_ac_2560_2, 3031040, P4TP_256
	PRCENTRY2		xfft_r4dwpn_1280K_ac_11_4, 2516992, P4TP_512
	PRCENTRY2		xfft_r4dwpn_1280K_ac_11_2, 2516992
	PRCENTRY2		xfft_r4dwpn_1280K_ac_1280_4, 3031040, CORE2_32
	PRCENTRY2		xfft_r4dwpn_1280K_ac_10_4, 2522112
	PRCENTRY2		xfft_r4delay_1280K_ac_5120_4, 4773888
	PRCENTRY2		xfft_r4delay_1280K_ac_2560_4, 3067904
	PRCENTRY2		xfft_r4delay_1280K_ac_2560_2, 3067904
	PRCENTRY2		xfft_r4delay_1280K_ac_11_4, 2610176
	PRCENTRY2		xfft_r4delay_1280K_ac_11_2, 2610176
	PRCENTRY2		xfft_r4delay_1280K_ac_1280_4, 3063808
	PRCENTRY2		xfft_r4delay_1280K_ac_10_4, 2619392
	DD			0
	PRCSTRT	27930000, 1474560, 0.116
	PRCENTRY2		xfft_r4dwpn_1440K_ac_3840_4, 3401728, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_1440K_ac_3840_2, 3401728, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_1440K_ac_3840_1, 3401728
	PRCENTRY2		xfft_r4dwpn_1440K_ac_2304_4, 2807808
	PRCENTRY2		xfft_r4delay_1440K_ac_3840_4, 3337216
	PRCENTRY2		xfft_r4delay_1440K_ac_3840_2, 3337216, P4TP_256
	PRCENTRY2		xfft_r4delay_1440K_ac_2304_4, 2913280
	DD			0
	PRCSTRT	29800000, 1572864, 0.125
	PRCENTRY2		xfft_r4dwpn_1536K_ac_6144_4, 6806528, P4_1024
	PRCENTRY2		xfft_r4dwpn_1536K_ac_12_4, 3639296, CORE2_32 + K10
	PRCENTRY2		xfft_r4dwpn_1536K_ac_12_2, 3639296, CORE2_512_32 + K8
	PRCENTRY2		xfft_r4dwpn_1536K_ac_12_1, 3639296, CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_1536K_ac_3072_4, 3633152, CORE2_64
	PRCENTRY2		xfft_r4dwpn_1536K_ac_3072_2, 3633152, P4TP_256
	PRCENTRY2		xfft_r4dwpn_1536K_ac_11_4, 3625984, I7
	PRCENTRY2		xfft_r4dwpn_1536K_ac_11_2, 3625984, P4TP_512
	PRCENTRY2		xfft_r4dwpn_1536K_ac_1536_4, 3624960
	PRCENTRY2		xfft_r4dwpn_1536K_ac_10_4, 2112512
	PRCENTRY2		xfft_r4delay_1536K_ac_6144_4, 5728256
	PRCENTRY2		xfft_r4delay_1536K_ac_12_4, 3570688
	PRCENTRY2		xfft_r4delay_1536K_ac_3072_4, 3678208
	PRCENTRY2		xfft_r4delay_1536K_ac_3072_2, 3678208
	PRCENTRY2		xfft_r4delay_1536K_ac_11_4, 3536896
	PRCENTRY2		xfft_r4delay_1536K_ac_11_2, 3536896
	PRCENTRY2		xfft_r4delay_1536K_ac_1536_4, 3661824
	PRCENTRY2		xfft_r4delay_1536K_ac_10_4, 2378752
	PRCENTRY2		xfft_hg_1536K_ac_12_4, 2480128
	PRCENTRY2		xfft_hg_1536K_ac_12_1, 2480128
	PRCENTRY2		xfft_hg_1536K_ac_11_1, 3252352
	DD			0
	PRCSTRT	31060000, 1638400, 0.131
	PRCENTRY2		xfft_r4dwpn_1600K_ac_6400_4, 7064576, P4_1024_64 + I7 + CORE2
	PRCENTRY2		xfft_r4dwpn_1600K_ac_6400_2, 7064576, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_1600K_ac_6400_1, 7064576, K8_32
	PRCENTRY2		xfft_r4dwpn_1600K_ac_2560_4, 3123200, K10_32
	PRCENTRY2A		xfft_r4dwpn_1600K_ac_2560_4_CORE, 3123200, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_1600K_ac_2560_2, 3123200, P4TP_512
	PRCENTRY2		xfft_r4dwpn_1600K_ac_1280_4, 3128320, K10_64
	PRCENTRY2		xfft_r4delay_1600K_ac_6400_4, 5941248
	PRCENTRY2		xfft_r4delay_1600K_ac_2560_4, 3240960
	PRCENTRY2		xfft_r4delay_1600K_ac_2560_2, 3240960
	PRCENTRY2		xfft_r4delay_1600K_ac_1280_4, 3246080
	DD			0
	PRCSTRT	33350000, 1769472, 0.142
	PRCENTRY2		xfft_r4dwpn_1728K_ac_4608_4, 4069376, P4_1024_32 + P4TP_512 + I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_1728K_ac_2304_4, 4056064, P4_1024_64 + I7_32 + CORE2 + K8_64
	PRCENTRY2		xfft_r4dwpn_1728K_ac_2304_2, 4056064, K8_32
	PRCENTRY2		xfft_r4dwpn_1728K_ac_2304_1, 4056064
	PRCENTRY2		xfft_r4delay_1728K_ac_4608_4, 3992576
	PRCENTRY2		xfft_r4delay_1728K_ac_2304_4, 3954688
	DD			0
	PRCSTRT	37110000, 1966080, 0.158
	PRCENTRY2		xfft_r4dwpn_1920K_ac_7680_4, 8477696, P4_1024
	PRCENTRY2		xfft_r4dwpn_1920K_ac_5120_4, 4540416, I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_1920K_ac_5120_2, 4540416, K8_64
	PRCENTRY2		xfft_r4dwpn_1920K_ac_5120_1, 4540416, CORE2_512_64 + K8_32
	PRCENTRY2		xfft_r4dwpn_1920K_ac_3840_4, 4517888, K10_32
	PRCENTRY2		xfft_r4dwpn_1920K_ac_3072_4, 3741696
	PRCENTRY2		xfft_r4dwpn_1920K_ac_3072_2, 3741696
	PRCENTRY2		xfft_r4dwpn_1920K_ac_2560_4, 4510720, I7_32
	PRCENTRY2		xfft_r4dwpn_1920K_ac_1536_4, 3738624
	PRCENTRY2		xfft_r4delay_1920K_ac_7680_4, 7129088
	PRCENTRY2		xfft_r4delay_1920K_ac_5120_4, 4455424
	PRCENTRY2		xfft_r4delay_1920K_ac_3840_4, 4575232
	PRCENTRY2		xfft_r4delay_1920K_ac_3072_4, 3884032
	PRCENTRY2		xfft_r4delay_1920K_ac_3072_2, 3884032
	PRCENTRY2		xfft_r4delay_1920K_ac_2560_4, 4397056
	PRCENTRY2		xfft_r4delay_1920K_ac_1536_4, 3876864
	DD			0
	PRCSTRT	39550000, 2097152, 0.169
	PRCENTRY2		xfft_r4dwpn_2M_ac_13_4, 9067520, P4_1024
	PRCENTRY2		xfft_r4dwpn_2M_ac_12_4, 4829184, CORE2 + K8 + K10
	PRCENTRY2		xfft_r4dwpn_2M_ac_12_2, 4829184, CORE2_512
	PRCENTRY2		xfft_r4dwpn_2M_ac_12_1, 4829184
	PRCENTRY2		xfft_r4dwpn_2M_ac_11_4, 4820992, I7
	PRCENTRY2		xfft_r4dwpn_2M_ac_10_4, 2879488
	PRCENTRY2		xfft_r4delay_2M_ac_13_4, 7628800
	PRCENTRY2		xfft_r4delay_2M_ac_12_4, 4890624
	PRCENTRY2		xfft_r4delay_2M_ac_12_2, 4890624
	PRCENTRY2		xfft_r4delay_2M_ac_12_1, 4890624
	PRCENTRY2		xfft_r4delay_2M_ac_11_4, 4866048
	PRCENTRY2		xfft_r4delay_2M_ac_10_4, 2940928
	PRCENTRY2		xfft_hg_2048K_ac_12_4, 3274752, P4TP_512
	PRCENTRY2		xfft_hg_2048K_ac_12_2, 3274752
	PRCENTRY2		xfft_hg_2048K_ac_12_1, 3274752
	DD			0
	PRCSTRT	44270000, 2359296, 0.199
	PRCENTRY2		xfft_r4dwpn_2304K_ac_9216_4, 10193920, P4_1024
	PRCENTRY2		xfft_r4dwpn_2304K_ac_6144_4, 5445632
	PRCENTRY2		xfft_r4dwpn_2304K_ac_4608_4, 5406720, K10_32
	PRCENTRY2		xfft_r4dwpn_2304K_ac_4608_2, 5406720
	PRCENTRY2		xfft_r4dwpn_2304K_ac_3072_4, 5407744, I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_2304K_ac_3072_2, 5407744, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_2304K_ac_3072_1, 5407744, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_2304K_ac_2304_4, 5398528, CORE2_32
	PRCENTRY2		xfft_r4delay_2304K_ac_9216_4, 8574976
	PRCENTRY2		xfft_r4delay_2304K_ac_6144_4, 5344256
	PRCENTRY2		xfft_r4delay_2304K_ac_4608_4, 5476352
	PRCENTRY2		xfft_r4delay_2304K_ac_4608_2, 5476352, P4TP_512
	PRCENTRY2		xfft_r4delay_2304K_ac_3072_4, 5269504
	PRCENTRY2		xfft_r4delay_2304K_ac_2304_4, 5447680
	DD			0
	PRCSTRT	46050000, 2457600, 0.211
	PRCENTRY2		xfft_r4dwpn_2400K_ac_6400_4, 5646336, P4_1024_64 + I7_64 + CORE2
	PRCENTRY2		xfft_r4dwpn_2400K_ac_6400_2, 5646336, CORE2_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_2400K_ac_6400_1, 5646336, K8_32
	PRCENTRY2		xfft_r4dwpn_2400K_ac_3840_4, 4651008, K10_64
	PRCENTRY2		xfft_r4dwpn_2400K_ac_3840_2, 4651008
	PRCENTRY2		xfft_r4delay_2400K_ac_6400_4, 5540864
	PRCENTRY2		xfft_r4delay_2400K_ac_3840_4, 4830208
	PRCENTRY2		xfft_r4delay_2400K_ac_3840_2, 4830208
	DD			0
	PRCSTRT	49170000, 2621440, 0.229
	PRCENTRY2		xfft_r4dwpn_2560K_ac_10240_4, 11328512, P4_1024
	PRCENTRY2		xfft_r4dwpn_2560K_ac_5120_4, 6025216, P4TP_512 + CORE2_64 + K10
	PRCENTRY2		xfft_r4dwpn_2560K_ac_5120_2, 6025216, CORE2_512 + K8
	PRCENTRY2		xfft_r4dwpn_2560K_ac_5120_1, 6025216
	PRCENTRY2		xfft_r4dwpn_2560K_ac_12_4, 4970496, CORE2_32
	PRCENTRY2		xfft_r4dwpn_2560K_ac_2560_4, 6000640
	PRCENTRY2A		xfft_r4dwpn_2560K_ac_2560_4_P4, 6000640, I7_32
	PRCENTRY2		xfft_r4dwpn_2560K_ac_11_4, 4967424, I7_64
	PRCENTRY2		xfft_r4dwpn_2560K_ac_10_4, 2720768
	PRCENTRY2		xfft_r4delay_2560K_ac_10240_4, 9529344
	PRCENTRY2		xfft_r4delay_2560K_ac_5120_4, 6103040
	PRCENTRY2		xfft_r4delay_2560K_ac_5120_2, 6103040
	PRCENTRY2		xfft_r4delay_2560K_ac_12_4, 5161984
	PRCENTRY2		xfft_r4delay_2560K_ac_2560_4, 6053888
	PRCENTRY2		xfft_r4delay_2560K_ac_11_4, 5146624
	PRCENTRY2		xfft_r4delay_2560K_ac_10_4, 3240960
	DD			0
	PRCSTRT	55050000, 2949120, 0.267
	PRCENTRY2		xfft_r4dwpn_2880K_ac_7680_4, 6772736
	PRCENTRY2		xfft_r4dwpn_2880K_ac_4608_4, 5564416, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_2880K_ac_4608_2, 5564416, P4TP_512
	PRCENTRY2		xfft_r4dwpn_2880K_ac_3840_4, 6734848, P4_1024_64 + I7_64 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_2880K_ac_3840_2, 6734848, K8
	PRCENTRY2		xfft_r4dwpn_2880K_ac_3840_1, 6734848, CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_2880K_ac_2304_4, 5561344
	PRCENTRY2		xfft_r4delay_2880K_ac_7680_4, 6646784
	PRCENTRY2		xfft_r4delay_2880K_ac_4608_4, 5780480
	PRCENTRY2		xfft_r4delay_2880K_ac_4608_2, 5780480
	PRCENTRY2		xfft_r4delay_2880K_ac_3840_4, 6559744
	PRCENTRY2		xfft_r4delay_2880K_ac_2304_4, 5761024
	DD			0
	PRCSTRT	58760000, 3145728, 0.289
	PRCENTRY2		xfft_r4dwpn_3M_ac_12288_4, 13585408
	PRCENTRY2		xfft_r4dwpn_3M_ac_13_4, 7247872
	PRCENTRY2		xfft_r4dwpn_3M_ac_6144_4, 7225344
	PRCENTRY2		xfft_r4dwpn_3M_ac_6144_2, 7225344
	PRCENTRY2		xfft_r4dwpn_3M_ac_12_4, 7193600, P4_1024 + I7 + CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_3M_ac_12_2, 7193600, CORE2_512_32
	PRCENTRY2		xfft_r4dwpn_3M_ac_12_1, 7193600, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_3M_ac_3072_4, 7192576
	PRCENTRY2		xfft_r4dwpn_3M_ac_10_4, 4143104
	PRCENTRY2		xfft_r4delay_3M_ac_12288_4, 11425792
	PRCENTRY2		xfft_r4delay_3M_ac_13_4, 7113728
	PRCENTRY2		xfft_r4delay_3M_ac_6144_4, 7319552, P4TP_512
	PRCENTRY2		xfft_r4delay_3M_ac_6144_2, 7319552
	PRCENTRY2		xfft_r4delay_3M_ac_12_4, 7006208
	PRCENTRY2		xfft_r4delay_3M_ac_3072_4, 7254016
	PRCENTRY2		xfft_r4delay_3M_ac_10_4, 4720640
	PRCENTRY2		xfft_hg_3072K_ac_13_1, 4143104
	PRCENTRY2		xfft_hg_3072K_ac_12_4, 4870272
	PRCENTRY2		xfft_hg_3072K_ac_12_2, 4870272
	PRCENTRY2		xfft_hg_3072K_ac_12_1, 4870272
	DD			0
	PRCSTRT	61130000, 3276800, 0.306
	PRCENTRY2		xfft_r4dwpn_3200K_ac_12800_4, 14113792, P4_1024
	PRCENTRY2		xfft_r4dwpn_3200K_ac_6400_4, 7499776, CORE2 + CORE2_512 + K8_32
	PRCENTRY2		xfft_r4dwpn_3200K_ac_5120_4, 6199296, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_3200K_ac_5120_2, 6199296, P4TP_512 + K8_64
	PRCENTRY2		xfft_r4dwpn_3200K_ac_5120_1, 6199296
	PRCENTRY2		xfft_r4dwpn_3200K_ac_2560_4, 6179840
	PRCENTRY2A		xfft_r4dwpn_3200K_ac_2560_4_P4, 6179840, I7_32
	PRCENTRY2		xfft_r4delay_3200K_ac_12800_4, 11864064
	PRCENTRY2		xfft_r4delay_3200K_ac_6400_4, 7598080
	PRCENTRY2		xfft_r4delay_3200K_ac_5120_4, 6439936
	PRCENTRY2		xfft_r4delay_3200K_ac_5120_2, 6439936
	PRCENTRY2		xfft_r4delay_3200K_ac_2560_4, 6400000
	DD			0
	PRCSTRT	65820000, 3538944, 0.340
	PRCENTRY2		xfft_r4dwpn_3456K_ac_9216_4, 8144896
	PRCENTRY2		xfft_r4dwpn_3456K_ac_9216_2, 8144896, CORE2_512_64
	PRCENTRY2		xfft_r4dwpn_3456K_ac_9216_1, 8144896, K8_32
	PRCENTRY2		xfft_r4dwpn_3456K_ac_4608_4, 8066048, P4_1024 + P4TP_512 + I7 + CORE2_64 + K10
	PRCENTRY2		xfft_r4dwpn_3456K_ac_2304_4, 4578304, CORE2_32 + K8_64
	PRCENTRY2		xfft_r4delay_3456K_ac_9216_4, 7994368
	PRCENTRY2		xfft_r4delay_3456K_ac_4608_4, 7854080
	PRCENTRY2		xfft_r4delay_3456K_ac_2304_4, 5192704
	DD			0
	PRCSTRT	73080000, 3932160, 0.391
	PRCENTRY2		xfft_r4dwpn_3840K_ac_15360_4, 16976896
	PRCENTRY2		xfft_r4dwpn_3840K_ac_10240_4, 9050112
	PRCENTRY2		xfft_r4dwpn_3840K_ac_7680_4, 8994816
	PRCENTRY2		xfft_r4dwpn_3840K_ac_6144_4, 7432192, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_3840K_ac_5120_4, 8979456, P4_1024_64 + I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_3840K_ac_5120_2, 8979456
	PRCENTRY2		xfft_r4dwpn_3840K_ac_5120_1, 8979456, CORE2_512_64 + K8
	PRCENTRY2		xfft_r4dwpn_3840K_ac_3840_4, 8962048, CORE2_32
	PRCENTRY2		xfft_r4dwpn_3840K_ac_3072_4, 7404544
	PRCENTRY2		xfft_r4delay_3840K_ac_15360_4, 14276608
	PRCENTRY2		xfft_r4delay_3840K_ac_10240_4, 8883200
	PRCENTRY2		xfft_r4delay_3840K_ac_7680_4, 9113600
	PRCENTRY2		xfft_r4delay_3840K_ac_6144_4, 7721984, P4TP_512
	PRCENTRY2		xfft_r4delay_3840K_ac_5120_4, 8742912
	PRCENTRY2		xfft_r4delay_3840K_ac_3840_4, 9035776
	PRCENTRY2		xfft_r4delay_3840K_ac_3072_4, 7665664
	DD			0
	PRCSTRT	75890000, 4096000, 0.412
	PRCENTRY2		xfft_r4dwpn_4000K_ac_6400_4, 7714816
	PRCENTRY2		xfft_r4delay_4000K_ac_6400_4, 8016896
	DD			0
	PRCSTRT	77950000, 4194304, 0.425
	PRCENTRY2		xfft_r4dwpn_4M_ac_14_4, 18107392, CORE2_64
	PRCENTRY2		xfft_r4dwpn_4M_ac_13_4, 9617408, P4TP_512
	PRCENTRY2		xfft_r4dwpn_4M_ac_13_2, 9617408
	PRCENTRY2		xfft_r4dwpn_4M_ac_13_1, 9617408
	PRCENTRY2		xfft_r4dwpn_4M_ac_12_4, 9568256, P4_1024 + I7 + CORE2_32 + K10
	PRCENTRY2		xfft_r4dwpn_4M_ac_12_2, 9568256, CORE2_512_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_4M_ac_12_1, 9568256, CORE2_512_32 + K8_32
	PRCENTRY2		xfft_r4dwpn_4M_ac_10_4, 4993024
	PRCENTRY2		xfft_r4delay_4M_ac_14_4, 15226880
	PRCENTRY2		xfft_r4delay_4M_ac_13_4, 9744384
	PRCENTRY2		xfft_r4delay_4M_ac_13_2, 9744384
	PRCENTRY2		xfft_r4delay_4M_ac_13_1, 9744384
	PRCENTRY2		xfft_r4delay_4M_ac_12_4, 9646080
	PRCENTRY2		xfft_r4delay_4M_ac_12_2, 9646080
	PRCENTRY2		xfft_r4delay_4M_ac_12_1, 9646080
	PRCENTRY2		xfft_r4delay_4M_ac_10_4, 5844992
	PRCENTRY2		xfft_hg_4096K_ac_13_1, 5462016
	PRCENTRY2		xfft_hg_4096K_ac_12_4, 6459392
	PRCENTRY2		xfft_hg_4096K_ac_12_2, 6459392
	PRCENTRY2		xfft_hg_4096K_ac_12_1, 6459392
	DD			0
	PRCSTRT	87170000, 4718592, 0.486
	PRCENTRY2		xfft_r4dwpn_4608K_ac_12288_4, 10848256
	PRCENTRY2		xfft_r4dwpn_4608K_ac_9216_4, 10809344, CORE2_64 + K8_64
	PRCENTRY2		xfft_r4dwpn_4608K_ac_6144_4, 10769408, P4_1024_32 + K10_32
	PRCENTRY2		xfft_r4dwpn_4608K_ac_4608_4, 10735616, P4_1024_64 + I7 + K8_32 + K10_64
	PRCENTRY2		xfft_r4dwpn_4608K_ac_3072_4, 6077440
	PRCENTRY2		xfft_r4dwpn_4608K_ac_2304_4, 6328320, CORE2_32
	PRCENTRY2		xfft_r4delay_4608K_ac_12288_4, 10648576
	PRCENTRY2		xfft_r4delay_4608K_ac_9216_4, 10952704
	PRCENTRY2		xfft_r4delay_4608K_ac_6144_4, 10483712
	PRCENTRY2		xfft_r4delay_4608K_ac_4608_4, 10821632
	PRCENTRY2		xfft_r4delay_4608K_ac_3072_4, 6900736
	PRCENTRY2		xfft_r4delay_4608K_ac_2304_4, 6410240
	DD			0
	PRCSTRT	90920000, 4915200, 0.509
	PRCENTRY2		xfft_r4dwpn_4800K_ac_12800_4, 11261952
	PRCENTRY2		xfft_r4dwpn_4800K_ac_7680_4, 9250816, P4_1024_32
	PRCENTRY2		xfft_r4dwpn_4800K_ac_6400_4, 11191296, P4_1024_64 + I7 + CORE2_64 + K8 + K10_64
	PRCENTRY2		xfft_r4dwpn_4800K_ac_3840_4, 9223168, CORE2_32
	PRCENTRY2		xfft_r4delay_4800K_ac_12800_4, 11054080
	PRCENTRY2		xfft_r4delay_4800K_ac_7680_4, 9614336
	PRCENTRY2		xfft_r4delay_4800K_ac_6400_4, 10893312
	PRCENTRY2		xfft_r4delay_4800K_ac_3840_4, 9545728
	DD			0
	PRCSTRT	96940000, 5242880, 0.550
	PRCENTRY2		xfft_r4dwpn_5M_ac_20480_4, 22629376
	PRCENTRY2		xfft_r4dwpn_5M_ac_10240_4, 12009472
	PRCENTRY2		xfft_r4dwpn_5M_ac_13_4, 9889792
	PRCENTRY2		xfft_r4dwpn_5M_ac_5120_4, 11943936, I7 + CORE2 + K8 + K10_32
	PRCENTRY2		xfft_r4dwpn_5M_ac_12_4, 9845760, P4_1024 + K10_64
	PRCENTRY2		xfft_r4dwpn_5M_ac_10_4, 5318656
	PRCENTRY2		xfft_r4delay_5M_ac_20480_4, 19027968
	PRCENTRY2		xfft_r4delay_5M_ac_10240_4, 12169216
	PRCENTRY2		xfft_r4delay_5M_ac_13_4, 10277888
	PRCENTRY2		xfft_r4delay_5M_ac_5120_4, 12038144
	PRCENTRY2		xfft_r4delay_5M_ac_12_4, 10188800
	PRCENTRY2		xfft_r4delay_5M_ac_10_4, 6445056
	DD			0
	PRCSTRT	108600000, 5898240, 0.623
	PRCENTRY2		xfft_r4dwpn_5760K_ac_15360_4, 13551616
	PRCENTRY2		xfft_r4dwpn_5760K_ac_9216_4, 11114496
	PRCENTRY2		xfft_r4dwpn_5760K_ac_7680_4, 13423616, I7_32 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_5760K_ac_4608_4, 11045888, P4_1024 + I7_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_5760K_ac_3840_4, 7552000, K8
	PRCENTRY2		xfft_r4dwpn_5760K_ac_2304_4, 5841920, CORE2_32
	PRCENTRY2		xfft_r4delay_5760K_ac_15360_4, 13302784
	PRCENTRY2		xfft_r4delay_5760K_ac_9216_4, 11551744
	PRCENTRY2		xfft_r4delay_5760K_ac_7680_4, 13064192
	PRCENTRY2		xfft_r4delay_5760K_ac_4608_4, 11429888
	PRCENTRY2		xfft_r4delay_5760K_ac_3840_4, 8584192
	PRCENTRY2		xfft_r4delay_5760K_ac_2304_4, 7037952
	DD			0
	PRCSTRT	115800000, 6291456, 0.668
	PRCENTRY2		xfft_r4dwpn_6M_ac_14_4, 14452736, CORE2_64
	PRCENTRY2		xfft_r4dwpn_6M_ac_12288_4, 14397440, CORE2_32
	PRCENTRY2		xfft_r4dwpn_6M_ac_13_4, 14341120, P4_1024_64
	PRCENTRY2		xfft_r4dwpn_6M_ac_6144_4, 14323712, P4_1024_32 + I7 + K10
	PRCENTRY2		xfft_r4dwpn_6M_ac_12_4, 8059904, K8
	PRCENTRY2		xfft_r4dwpn_6M_ac_3072_4, 8417280
	PRCENTRY2		xfft_r4delay_6M_ac_14_4, 14187520
	PRCENTRY2		xfft_r4delay_6M_ac_12288_4, 14589952
	PRCENTRY2		xfft_r4delay_6M_ac_13_4, 13957120
	PRCENTRY2		xfft_r4delay_6M_ac_6144_4, 14434304
	PRCENTRY2		xfft_r4delay_6M_ac_12_4, 9161728
	PRCENTRY2		xfft_r4delay_6M_ac_3072_4, 8511488
	PRCENTRY2		xfft_hg_6M_ac_12_4, 9650176
	PRCENTRY2		xfft_hg_6M_ac_12_2, 9650176
	PRCENTRY2		xfft_hg_6M_ac_12_1, 9650176
	DD			0
	PRCSTRT	120500000, 6553600, 0.715
	PRCENTRY2		xfft_r4dwpn_6400K_ac_25600_4, 28281856
	PRCENTRY2		xfft_r4dwpn_6400K_ac_12800_4, 14958592
	PRCENTRY2		xfft_r4dwpn_6400K_ac_10240_4, 12347392
	PRCENTRY2		xfft_r4dwpn_6400K_ac_6400_4, 14893056, CORE2
	PRCENTRY2		xfft_r4dwpn_6400K_ac_5120_4, 12286976, I7 + K10
	PRCENTRY2		xfft_r4delay_6400K_ac_25600_4, 23779328
	PRCENTRY2		xfft_r4delay_6400K_ac_12800_4, 15159296
	PRCENTRY2		xfft_r4delay_6400K_ac_10240_4, 12833792
	PRCENTRY2		xfft_r4delay_6400K_ac_6400_4, 15007744
	PRCENTRY2		xfft_r4delay_6400K_ac_5120_4, 12711936
	DD			0
	PRCSTRT	129400000, 7077888, 0.760
	PRCENTRY2		xfft_r4dwpn_6912K_ac_9216_4, 16122880, I7 + CORE2_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_6912K_ac_4608_4, 9030656, K10_32
	PRCENTRY2		xfft_r4dwpn_6912K_ac_2304_4, 9066496, CORE2_32
	PRCENTRY2		xfft_r4delay_6912K_ac_9216_4, 15689728
	PRCENTRY2		xfft_r4delay_6912K_ac_4608_4, 10271744
	PRCENTRY2		xfft_r4delay_6912K_ac_2304_4, 10319872
	DD			0
	PRCSTRT	143700000, 7864320, 0.851
	PRCENTRY2		xfft_r4dwpn_7680K_ac_20480_4, 18057216
	PRCENTRY2		xfft_r4dwpn_7680K_ac_15360_4, 17985536
	PRCENTRY2		xfft_r4dwpn_7680K_ac_12288_4, 14800896
	PRCENTRY2		xfft_r4dwpn_7680K_ac_10240_4, 17912832, CORE2_64
	PRCENTRY2		xfft_r4dwpn_7680K_ac_7680_4, 17862656, I7
	PRCENTRY2		xfft_r4dwpn_7680K_ac_6144_4, 14732288, K10_32
	PRCENTRY2		xfft_r4dwpn_7680K_ac_5120_4, 10042368, K10_64
	PRCENTRY2		xfft_r4dwpn_7680K_ac_3840_4, 10481664, CORE2_32
	PRCENTRY2		xfft_r4dwpn_7680K_ac_3072_4, 7734272
	PRCENTRY2		xfft_r4delay_7680K_ac_20480_4, 17726464
	PRCENTRY2		xfft_r4delay_7680K_ac_15360_4, 18227200
	PRCENTRY2		xfft_r4delay_7680K_ac_12288_4, 15385600
	PRCENTRY2		xfft_r4delay_7680K_ac_10240_4, 17430528
	PRCENTRY2		xfft_r4delay_7680K_ac_7680_4, 17997824
	PRCENTRY2		xfft_r4delay_7680K_ac_6144_4, 15239168
	PRCENTRY2		xfft_r4delay_7680K_ac_5120_4, 11422720
	PRCENTRY2		xfft_r4delay_7680K_ac_3840_4, 10588160
	PRCENTRY2		xfft_r4delay_7680K_ac_3072_4, 9335808
	DD			0
	PRCSTRT	149700000, 8192000, 1.007
	PRCENTRY2		xfft_r4dwpn_8000K_ac_12800_4, 15378432
	PRCENTRY2		xfft_r4dwpn_8000K_ac_6400_4, 15318016, I7_64 + K10_64
	PRCENTRY2		xfft_r4delay_8000K_ac_12800_4, 15987712
	PRCENTRY2		xfft_r4delay_8000K_ac_6400_4, 15845376
	DD			0
	PRCSTRT	153400000, 8388608, 1.042
	PRCENTRY2		xfft_r4dwpn_8M_ac_14_4, 19181568, CORE2_64
	PRCENTRY2		xfft_r4dwpn_8M_ac_13_4, 19075072, I7 + K10
	PRCENTRY2		xfft_r4dwpn_8M_ac_12_4, 11186176, CORE2_32
	PRCENTRY2		xfft_r4delay_8M_ac_14_4, 19439616
	PRCENTRY2		xfft_r4delay_8M_ac_13_4, 19218432
	PRCENTRY2		xfft_r4delay_8M_ac_12_4, 11296768
	PRCENTRY2		xfft_hg_8M_ac_12_4, 12828672
	PRCENTRY2		xfft_hg_8M_ac_12_2, 12828672
	PRCENTRY2		xfft_hg_8M_ac_12_1, 12828672
	DD			0
	PRCSTRT	171700000, 9437184, 1.131
	PRCENTRY2		xfft_r4dwpn_9M_ac_12288_4, 21480448, K10
	PRCENTRY2		xfft_r4dwpn_9M_ac_9216_4, 21446656, I7 + CORE2
	PRCENTRY2		xfft_r4dwpn_9M_ac_6144_4, 12028928
	PRCENTRY2		xfft_r4dwpn_9M_ac_4608_4, 12550144
	PRCENTRY2		xfft_r4dwpn_9M_ac_3072_4, 12040192
	PRCENTRY2		xfft_r4dwpn_9M_ac_2304_4, 10899456
	PRCENTRY2		xfft_r4delay_9M_ac_12288_4, 20899840
	PRCENTRY2		xfft_r4delay_9M_ac_9216_4, 21606400
	PRCENTRY2		xfft_r4delay_9M_ac_6144_4, 13687808
	PRCENTRY2		xfft_r4delay_9M_ac_4608_4, 12668928
	PRCENTRY2		xfft_r4delay_9M_ac_3072_4, 13699072
	PRCENTRY2		xfft_r4delay_9M_ac_2304_4, 12754944
	DD			0
	PRCSTRT	178900000, 9830400, 1.166
	PRCENTRY2		xfft_r4dwpn_9600K_ac_25600_4, 22562816
	PRCENTRY2		xfft_r4dwpn_9600K_ac_15360_4, 18487296
	PRCENTRY2		xfft_r4dwpn_9600K_ac_12800_4, 22336512, I7 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_9600K_ac_7680_4, 18369536
	PRCENTRY2		xfft_r4dwpn_9600K_ac_6400_4, 12499968, K10_64
	PRCENTRY2		xfft_r4dwpn_9600K_ac_3840_4, 9602048, CORE2_32
	PRCENTRY2		xfft_r4delay_9600K_ac_25600_4, 22150144
	PRCENTRY2		xfft_r4delay_9600K_ac_15360_4, 19219456
	PRCENTRY2		xfft_r4delay_9600K_ac_12800_4, 21731328
	PRCENTRY2		xfft_r4delay_9600K_ac_7680_4, 18999296
	PRCENTRY2		xfft_r4delay_9600K_ac_6400_4, 14228480
	PRCENTRY2		xfft_r4delay_9600K_ac_3840_4, 11609088
	DD			0
	PRCSTRT	190700000, 10485760, 1.221
	PRCENTRY2		xfft_r4dwpn_10M_ac_20480_4, 23965696
	PRCENTRY2		xfft_r4dwpn_10M_ac_14_4, 19716096, CORE2_32
	PRCENTRY2		xfft_r4dwpn_10M_ac_10240_4, 23826432, I7 + CORE2_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_10M_ac_13_4, 19614720, K10_64
	PRCENTRY2		xfft_r4dwpn_10M_ac_5120_4, 13955072
	PRCENTRY2		xfft_r4dwpn_10M_ac_12_4, 10241024
	PRCENTRY2		xfft_r4delay_10M_ac_20480_4, 24289280
	PRCENTRY2		xfft_r4delay_10M_ac_14_4, 20497408
	PRCENTRY2		xfft_r4delay_10M_ac_10240_4, 24002560
	PRCENTRY2		xfft_r4delay_10M_ac_13_4, 20285440
	PRCENTRY2		xfft_r4delay_10M_ac_5120_4, 14082048
	PRCENTRY2		xfft_r4delay_10M_ac_12_4, 12383232
	DD			0
	PRCSTRT	213600000, 11796480, 1.330
	PRCENTRY2		xfft_r4dwpn_11520K_ac_15360_4, 26838016, I7 + CORE2_64
	PRCENTRY2		xfft_r4dwpn_11520K_ac_9216_4, 22051840, K10_64
	PRCENTRY2		xfft_r4dwpn_11520K_ac_7680_4, 14978048
	PRCENTRY2		xfft_r4dwpn_11520K_ac_4608_4, 11473920
	PRCENTRY2		xfft_r4dwpn_11520K_ac_3840_4, 14989312, CORE2_32
	PRCENTRY2		xfft_r4dwpn_11520K_ac_2304_4, 11552768
	PRCENTRY2		xfft_r4delay_11520K_ac_15360_4, 26109952
	PRCENTRY2		xfft_r4delay_11520K_ac_9216_4, 22804480
	PRCENTRY2		xfft_r4delay_11520K_ac_7680_4, 17054720
	PRCENTRY2		xfft_r4delay_11520K_ac_4608_4, 13886464
	PRCENTRY2		xfft_r4delay_11520K_ac_3840_4, 17053696
	PRCENTRY2		xfft_r4delay_11520K_ac_2304_4, 14010368
	DD			0
	PRCSTRT	228000000, 12582912, 1.400
	PRCENTRY2		xfft_r4dwpn_12M_ac_14_4, 28623872, I7_64 + CORE2
	PRCENTRY2		xfft_r4dwpn_12M_ac_12288_4, 28573696, I7_32 + K10
	PRCENTRY2		xfft_r4dwpn_12M_ac_13_4, 15993856
	PRCENTRY2		xfft_r4dwpn_12M_ac_6144_4, 16728064
	PRCENTRY2		xfft_r4dwpn_12M_ac_12_4, 15988736
	PRCENTRY2		xfft_r4dwpn_12M_ac_3072_4, 14462976
	PRCENTRY2		xfft_r4delay_12M_ac_14_4, 27846656
	PRCENTRY2		xfft_r4delay_12M_ac_12288_4, 28782592
	PRCENTRY2		xfft_r4delay_12M_ac_13_4, 18209792
	PRCENTRY2		xfft_r4delay_12M_ac_6144_4, 16871424
	PRCENTRY2		xfft_r4delay_12M_ac_12_4, 18188288
	PRCENTRY2		xfft_r4delay_12M_ac_3072_4, 16920576
	PRCENTRY2		xfft_hg_12M_ac_13_4, 16031744
	DD			0
	PRCSTRT	237200000, 13107200, 1.488
	PRCENTRY2		xfft_r4dwpn_12800K_ac_25600_4, 29945856, CORE2_64
	PRCENTRY2		xfft_r4dwpn_12800K_ac_20480_4, 24631296
	PRCENTRY2		xfft_r4dwpn_12800K_ac_12800_4, 29724672, I7
	PRCENTRY2		xfft_r4dwpn_12800K_ac_10240_4, 24497152, K10
	PRCENTRY2		xfft_r4dwpn_12800K_ac_6400_4, 17395712, CORE2_32
	PRCENTRY2		xfft_r4dwpn_12800K_ac_5120_4, 12747776
	PRCENTRY2		xfft_r4delay_12800K_ac_25600_4, 30351360
	PRCENTRY2		xfft_r4delay_12800K_ac_20480_4, 25609216
	PRCENTRY2		xfft_r4delay_12800K_ac_12800_4, 29941760
	PRCENTRY2		xfft_r4delay_12800K_ac_10240_4, 25331712
	PRCENTRY2		xfft_r4delay_12800K_ac_6400_4, 17543168
	PRCENTRY2		xfft_r4delay_12800K_ac_5120_4, 15430656
	DD			0
	PRCSTRT	254500000, 14155776, 1.663
	PRCENTRY2		xfft_r4dwpn_13824K_ac_9216_4, 17972224, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_13824K_ac_4608_4, 17942528, I7
	PRCENTRY2		xfft_r4delay_13824K_ac_9216_4, 20466688
	PRCENTRY2		xfft_r4delay_13824K_ac_4608_4, 20412416
	DD			0
	PRCSTRT	283500000, 15728640, 1.875
	PRCENTRY2		xfft_r4dwpn_15M_ac_20480_4, 35767296
	PRCENTRY2		xfft_r4dwpn_15M_ac_15360_4, 35700736, I7
	PRCENTRY2		xfft_r4dwpn_15M_ac_12288_4, 29375488, K10
	PRCENTRY2		xfft_r4dwpn_15M_ac_10240_4, 19958784
	PRCENTRY2		xfft_r4dwpn_15M_ac_7680_4, 20856832
	PRCENTRY2		xfft_r4dwpn_15M_ac_6144_4, 15258624
	PRCENTRY2		xfft_r4dwpn_15M_ac_5120_4, 19937280, CORE2
	PRCENTRY2		xfft_r4dwpn_15M_ac_3840_4, 18001920
	PRCENTRY2		xfft_r4dwpn_15M_ac_3072_4, 15312896
	PRCENTRY2		xfft_r4delay_15M_ac_20480_4, 34793472
	PRCENTRY2		xfft_r4delay_15M_ac_15360_4, 35958784
	PRCENTRY2		xfft_r4delay_15M_ac_12288_4, 30373888
	PRCENTRY2		xfft_r4delay_15M_ac_10240_4, 22731776
	PRCENTRY2		xfft_r4delay_15M_ac_7680_4, 21024768
	PRCENTRY2		xfft_r4delay_15M_ac_6144_4, 18482176
	PRCENTRY2		xfft_r4delay_15M_ac_5120_4, 22677504
	PRCENTRY2		xfft_r4delay_15M_ac_3840_4, 21061632
	PRCENTRY2		xfft_r4delay_15M_ac_3072_4, 18569216
	DD			0
	PRCSTRT	295000000, 16384000, 2.078
	PRCENTRY2		xfft_r4dwpn_16000K_ac_25600_4, 30775296
	PRCENTRY2		xfft_r4dwpn_16000K_ac_12800_4, 30559232, I7_64 + K10_64
	PRCENTRY2		xfft_r4dwpn_16000K_ac_6400_4, 15860736
	PRCENTRY2		xfft_r4delay_16000K_ac_25600_4, 31998976
	PRCENTRY2		xfft_r4delay_16000K_ac_12800_4, 31598592
	PRCENTRY2		xfft_r4delay_16000K_ac_6400_4, 19219456
	DD			0
	PRCSTRT	302400000, 16777216, 2.100
	PRCENTRY2		xfft_r4dwpn_16M_ac_14_4, 38076416, CORE2 + K10_64
	PRCENTRY2		xfft_r4dwpn_16M_ac_13_4, 22265856
	PRCENTRY2		xfft_r4dwpn_16M_ac_13_2, 22265856, K10_32
	PRCENTRY2		xfft_r4dwpn_16M_ac_13_1, 22265856
	PRCENTRY2		xfft_r4dwpn_16M_ac_12_4, 19197952, I7
	PRCENTRY2		xfft_r4dwpn_16M_ac_12_2, 19197952
	PRCENTRY2		xfft_r4dwpn_16M_ac_12_1, 19197952
	PRCENTRY2		xfft_r4delay_16M_ac_14_4, 38350848
	PRCENTRY2		xfft_r4delay_16M_ac_13_4, 22441984
	PRCENTRY2		xfft_r4delay_16M_ac_13_2, 22441984
	PRCENTRY2		xfft_r4delay_16M_ac_13_1, 22441984
	PRCENTRY2		xfft_r4delay_16M_ac_12_4, 22458368
	PRCENTRY2		xfft_r4delay_16M_ac_12_2, 22458368
	PRCENTRY2		xfft_r4delay_16M_ac_12_1, 22458368
	PRCENTRY2		xfft_hg_16M_ac_13_4, 21307392
	DD			0
	PRCSTRT	338900000, 18874368, 2.400
	PRCENTRY2		xfft_r4dwpn_18M_ac_12288_4, 23919616, K10
	PRCENTRY2		xfft_r4dwpn_18M_ac_9216_4, 25030656, CORE2
	PRCENTRY2		xfft_r4dwpn_18M_ac_6144_4, 23889920
	PRCENTRY2		xfft_r4dwpn_18M_ac_4608_4, 21544960, I7
	PRCENTRY2		xfft_r4delay_18M_ac_12288_4, 27249664
	PRCENTRY2		xfft_r4delay_18M_ac_9216_4, 25223168
	PRCENTRY2		xfft_r4delay_18M_ac_6144_4, 27170816
	PRCENTRY2		xfft_r4delay_18M_ac_4608_4, 25206784
	DD			0
	PRCSTRT	352400000, 19660800, 2.469
	PRCENTRY2		xfft_r4dwpn_19200K_ac_25600_4, 44696576, I7 + CORE2_64 + K10_32
	PRCENTRY2		xfft_r4dwpn_19200K_ac_15360_4, 36699136, K10_64
	PRCENTRY2		xfft_r4dwpn_19200K_ac_12800_4, 24873984
	PRCENTRY2		xfft_r4dwpn_19200K_ac_7680_4, 18994176
	PRCENTRY2		xfft_r4dwpn_19200K_ac_6400_4, 24852480, CORE2_32
	PRCENTRY2		xfft_r4dwpn_19200K_ac_3840_4, 19048448
	PRCENTRY2		xfft_r4delay_19200K_ac_25600_4, 43476992
	PRCENTRY2		xfft_r4delay_19200K_ac_15360_4, 37943296
	PRCENTRY2		xfft_r4delay_19200K_ac_12800_4, 28343296
	PRCENTRY2		xfft_r4delay_19200K_ac_7680_4, 23028736
	PRCENTRY2		xfft_r4delay_19200K_ac_6400_4, 28268544
	PRCENTRY2		xfft_r4delay_19200K_ac_3840_4, 23103488
	DD			0
	PRCSTRT	375900000, 20971520, 2.700
	PRCENTRY2		xfft_r4dwpn_20M_ac_20480_4, 47579136, K10_64
	PRCENTRY2		xfft_r4dwpn_20M_ac_14_4, 39140352, CORE2_64
	PRCENTRY2		xfft_r4dwpn_20M_ac_10240_4, 27803648
	PRCENTRY2		xfft_r4dwpn_20M_ac_10240_2, 27803648, K10_32
	PRCENTRY2		xfft_r4dwpn_20M_ac_13_4, 20272128
	PRCENTRY2		xfft_r4dwpn_20M_ac_5120_4, 23932928, I7_32
	PRCENTRY2		xfft_r4dwpn_20M_ac_12_4, 20310016, I7_64 + CORE2_32
	PRCENTRY2		xfft_r4delay_20M_ac_20480_4, 47919104
	PRCENTRY2		xfft_r4delay_20M_ac_14_4, 40466432
	PRCENTRY2		xfft_r4delay_20M_ac_10240_4, 28012544
	PRCENTRY2		xfft_r4delay_20M_ac_10240_2, 28012544
	PRCENTRY2		xfft_r4delay_20M_ac_13_4, 24577024
	PRCENTRY2		xfft_r4delay_20M_ac_5120_4, 27996160
	PRCENTRY2		xfft_r4delay_20M_ac_12_4, 24631296
	DD			0
	PRCSTRT	419600000, 23592960, 3.075
	PRCENTRY2		xfft_r4dwpn_23040K_ac_15360_4, 29867008, K10
	PRCENTRY2		xfft_r4dwpn_23040K_ac_9216_4, 22774784, CORE2
	PRCENTRY2		xfft_r4dwpn_23040K_ac_7680_4, 29788160
	PRCENTRY2		xfft_r4dwpn_23040K_ac_4608_4, 22788096, I7
	PRCENTRY2		xfft_r4delay_23040K_ac_15360_4, 34032640
	PRCENTRY2		xfft_r4delay_23040K_ac_9216_4, 27620352
	PRCENTRY2		xfft_r4delay_23040K_ac_7680_4, 33880064
	PRCENTRY2		xfft_r4delay_23040K_ac_4608_4, 27641856
	DD			0
	PRCSTRT	449200000, 25165824, 3.300
	PRCENTRY2		xfft_r4dwpn_24M_ac_14_4, 31849472, CORE2
	PRCENTRY2		xfft_r4dwpn_24M_ac_12288_4, 33337344, K10
	PRCENTRY2		xfft_r4dwpn_24M_ac_13_4, 31787008
	PRCENTRY2		xfft_r4dwpn_24M_ac_6144_4, 28672000, I7_32
	PRCENTRY2		xfft_r4delay_24M_ac_14_4, 36293632
	PRCENTRY2		xfft_r4delay_24M_ac_12288_4, 33579008
	PRCENTRY2		xfft_r4delay_24M_ac_13_4, 36149248
	PRCENTRY2		xfft_r4delay_24M_ac_6144_4, 33538048
	PRCENTRY2		xfft_hg_24M_ac_13_2, 31883392
	DD			0
	PRCSTRT	467200000, 26214400, 3.450
	PRCENTRY2		xfft_r4dwpn_25M_ac_25600_4, 59457536, I7 + CORE2_64 + K10
	PRCENTRY2		xfft_r4dwpn_25M_ac_20480_4, 48905216, CORE2_32
	PRCENTRY2		xfft_r4dwpn_25M_ac_12800_4, 34684928
	PRCENTRY2		xfft_r4dwpn_25M_ac_10240_4, 25285632
	PRCENTRY2		xfft_r4dwpn_25M_ac_6400_4, 29831168
	PRCENTRY2		xfft_r4dwpn_25M_ac_5120_4, 25307136
	PRCENTRY2		xfft_r4delay_25M_ac_25600_4, 59879424
	PRCENTRY2		xfft_r4delay_25M_ac_20480_4, 50558976
	PRCENTRY2		xfft_r4delay_25M_ac_12800_4, 34934784
	PRCENTRY2		xfft_r4delay_25M_ac_10240_4, 30671872
	PRCENTRY2		xfft_r4delay_25M_ac_6400_4, 34897920
	PRCENTRY2		xfft_r4delay_25M_ac_5120_4, 30693376
	DD			0
	PRCSTRT	501700000, 28311552, 3.750
	PRCENTRY2		xfft_r4dwpn_27M_ac_9216_4, 35731456, I7 + CORE2 + K10
	PRCENTRY2		xfft_r4delay_27M_ac_9216_4, 40634368
	DD			0
	PRCSTRT	557100000, 31457280, 4.200
	PRCENTRY2		xfft_r4dwpn_30M_ac_20480_4, 39779328, K10_32
	PRCENTRY2		xfft_r4dwpn_30M_ac_15360_4, 41644032, K10_64
	PRCENTRY2		xfft_r4dwpn_30M_ac_12288_4, 30295040, CORE2_32
	PRCENTRY2		xfft_r4dwpn_30M_ac_10240_4, 39684096, CORE2_64
	PRCENTRY2		xfft_r4dwpn_30M_ac_7680_4, 35749888
	PRCENTRY2		xfft_r4dwpn_30M_ac_6144_4, 30308352, I7
	PRCENTRY2		xfft_r4delay_30M_ac_20480_4, 45337600
	PRCENTRY2		xfft_r4delay_30M_ac_15360_4, 41934848
	PRCENTRY2		xfft_r4delay_30M_ac_12288_4, 36762624
	PRCENTRY2		xfft_r4delay_30M_ac_10240_4, 45127680
	PRCENTRY2		xfft_r4delay_30M_ac_7680_4, 41820160
	PRCENTRY2		xfft_r4delay_30M_ac_6144_4, 36759552
	DD			0
	PRCSTRT	580600000, 32768000, 4.388
	PRCENTRY2		xfft_r4dwpn_32000K_ac_25600_4, 61111296, I7_64 + K10
	PRCENTRY2		xfft_r4dwpn_32000K_ac_12800_4, 31511552
	PRCENTRY2		xfft_r4dwpn_32000K_ac_6400_4, 31533056
	PRCENTRY2		xfft_r4delay_32000K_ac_25600_4, 63174656
	PRCENTRY2		xfft_r4delay_32000K_ac_12800_4, 38249472
	PRCENTRY2		xfft_r4delay_32000K_ac_6400_4, 38250496
	DD			0
	PRCSTRT	594800000, 33554432, 4.500
	PRCENTRY2		xfft_r4dwpn_32M_ac_14_4, 44412928, CORE2 + K10
	PRCENTRY2		xfft_r4dwpn_32M_ac_13_4, 38141952, I7
	PRCENTRY2		xfft_r4dwpn_32M_ac_13_2, 38141952
	PRCENTRY2		xfft_r4dwpn_32M_ac_13_1, 38141952
	PRCENTRY2		xfft_r4delay_32M_ac_14_4, 44720128
	PRCENTRY2		xfft_r4delay_32M_ac_13_4, 44613632
	PRCENTRY2		xfft_r4delay_32M_ac_13_2, 44613632
	PRCENTRY2		xfft_r4delay_32M_ac_13_1, 44613632
	PRCENTRY2		xfft_hg_32M_ac_13_2, 42434560
	DD			0
	DD	0

;; Jump tables for the AVX & FMA3 optimized code

yjmptable DD	0
	org	$-4
	PRCSTRT	743,	32,	0.00000015
	PRCENTRY2		yfft_r4_32_op, 712, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	1469,	64,	0.0000002
	PRCENTRY2		yfft_r4_64_op, 1424, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	2179,	96,	0.0000003
	PRCENTRY2		yfft_r4_96_op, 2648, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	2905,	128,	0.0000004
	PRCENTRY2		yfft_r4_128_op, 3616, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	3613,	160,	0.0000005
	PRCENTRY2		yfft_r4_160_op, 4584, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	4311,	192,	0.0000006
	PRCENTRY2		yfft_r4_192_op, 5296, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	5755,	256,	0.0000007
	PRCENTRY2		yfft_r4_256_op, 7232, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	7149,	320,	0.00000093
	PRCENTRY2		yfft_r4_320_op, 9168, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	8527,	384,	0.000001
	PRCENTRY2		yfft_r4_384_op, 11360, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	11359,	512,	0.000002
	PRCENTRY2		yfft_r4_512_op, 14464, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	14119,	640,	0.0000025
	PRCENTRY2		yfft_r4_640_op, 19104, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	16839,	768,	0.000003
	PRCENTRY2		yfft_r4_768_op, 22720, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	22473,	1024,	0.000004
	PRCENTRY2		yfft_r4_1K_op, 28928, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	27913,	1280,	0.000005
	PRCENTRY2		yfft_r4_1280_op, 38208, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	33375,	1536,	0.000006
	PRCENTRY2		yfft_r4_1536_op, 45440, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	44363,	2048,	0.000008
	PRCENTRY2		yfft_r4_2K_op, 57856, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	55153,	2560,	0.000011
	PRCENTRY2		yfft_r4_2560_op, 76416, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	65887,	3072,	0.000013
	PRCENTRY2		yfft_r4_3K_op, 90880, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	87477,	4096,	0.000017
	PRCENTRY2		yfft_r4_4K_op, 115712, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	108800,	5120,	0.000022
	PRCENTRY2		yfft_r4_5K_op, 152832, I7 + FMA3_64 + RYZEN_64
	DD			0
	PRCSTRT	130000,	6144,	0.000025
	PRCENTRY2		yfft_r4_6K_op, 181760
	YPRCENTRY421		yfft_r4dwpn,6K,,48,27712,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	172700,	8192,	0.000033
	PRCENTRY2		yfft_r4_8K_op, 231424
	YPRCENTRY421		yfft_r4dwpn,8K,,6,35968,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	214600,	10240,	0.000043
	PRCENTRY2		yfft_r4_10K_op, 305664
	YPRCENTRY421		yfft_r4dwpn,10K,,80,43072,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	256400,	12288,	0.000055
	PRCENTRY2		yfft_r4_12K_op, 363520
	YPRCENTRY421		yfft_r4dwpn,12K,,48,36928,,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	319200,	15360,	0.000073
	YPRCENTRY421		yfft_r4dwpn,15K,,48,41536,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	340700,	16384,	0.000072
	PRCENTRY2		yfft_r4_16K_op, 462848, RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,16K,,6,46720,,I7_32+FMA3_64,I7_64
	DD			0
	PRCSTRT	381700,	18432,	0.000086
	PRCENTRY2		yfft_r4_18K_op, 520704
	YPRCENTRY421		yfft_r4dwpn,18K,,48,46144,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	424300,	20480,	0.000089
	PRCENTRY2		yfft_r4_20K_op, 611328, RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,20K,,80,55360
	YPRCENTRY421		yfft_r4dwpn,20K,,6,52096,FMA3_64,,I7
	DD			0
	PRCSTRT	444400,	21504,	0.000108
	YPRCENTRY421		yfft_r4dwpn,21K,,48,50752,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	506800,	24576,	0.000105
	PRCENTRY2		yfft_r4_24K_op, 727040
	YPRCENTRY421		yfft_r4dwpn,24K,,192,97408,,I7
	YPRCENTRY421		yfft_r4dwpn,24K,,6,57472,,FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,24K,,48,62272
	DD			0
	PRCSTRT	527800,	25600,	0.000123
	YPRCENTRY421		yfft_r4dwpn,25K,,80,61504,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	590400,	28672,	0.000141
	YPRCENTRY421		yfft_r4dwpn,28K,,6,62848,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	630600,	30720,	0.000147
	YPRCENTRY421		yfft_r4dwpn,30K,,80,67648,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,30K,,48,68672,,FMA3_64
	DD			0
	PRCSTRT	672900,	32768,	0.000133
	PRCENTRY2		yfft_r4_32K_op, 925696
	YPRCENTRY421		yfft_r4dwpn,32K,,8,130816,,I7
	YPRCENTRY421		yfft_r4dwpn,32K,,6,79232,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	734400,	35840,	0.000179
	YPRCENTRY421		yfft_r4dwpn,35K,,80,73792,,RYZEN_64,I7
	DD			0
	PRCSTRT	754400,	36864,	0.000183
	YPRCENTRY421		yfft_r4dwpn,36K,,48,75072,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	836300,	40960,	0.000183
	YPRCENTRY421		yfft_r4dwpn,40K,,320,158848,,I7_32
	YPRCENTRY421		yfft_r4dwpn,40K,,80,95040,,RYZEN_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,40K,,6,86656,,FMA3_64
	DD			0
	PRCSTRT	1000000, 49152,	0.000232
	YPRCENTRY421		yfft_r4dwpn,48K,,192,120448,,,I7
	YPRCENTRY421		yfft_r4dwpn,48K,,6,94080,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	1041000, 51200,	0.000260
	YPRCENTRY421		yfft_r4dwpn,50K,,80,103488,,FMA3_64,I7+RYZEN_64
	DD			0
	PRCSTRT	1244000, 61440,	0.000301
	YPRCENTRY421		yfft_r4dwpn,60K,,192,131968,,I7
	YPRCENTRY421		yfft_r4dwpn,60K,,80,111936,,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	1327000, 65536,	0.000300
	YPRCENTRY421		yfft_r4dwpn,64K,,8,160000,,I7_32+FMA3_64+RYZEN_64,I7_64
	DD			0
	PRCSTRT	1487000, 73728,	0.000367
	YPRCENTRY421		yfft_r4dwpn,72K,,192,143488,FMA3_64,RYZEN_64,I7
	DD			0
	PRCSTRT	1652000, 81920, 0.000386
	YPRCENTRY421		yfft_r4dwpn,80K,,320,194176
	YPRCENTRY421		yfft_r4dwpn,80K,,8,174592,FMA3_64,RYZEN_64,I7
	DD			0
	PRCSTRT	1730000, 86016, 0.000443
	YPRCENTRY421		yfft_r4dwpn,84K,,192,155008,FMA3_64,,I7+RYZEN_64
	DD			0
	PRCSTRT	1975000, 98304, 0.000458
	YPRCENTRY421		yfft_r4dwpn,96K,,768,376576
	YPRCENTRY421		yfft_r4dwpn,96K,,8,189184,FMA3_64,I7_64+RYZEN_64,I7_32
	YPRCENTRY421		yfft_r4dwpn,96K,,192,210304
	DD			0
	PRCSTRT	2055000, 102400, 0.000508
	YPRCENTRY421		yfft_r4dwpn,100K,,320,211840,FMA3_64,I7+RYZEN_64
	DD			0
	PRCSTRT	2299000, 114688, 0.000567
	YPRCENTRY421		yfft_r4dwpn,112K,,8,203776,,I7_64+FMA3_64+RYZEN_64,I7_32
	DD			0
	PRCSTRT	2457000, 122880, 0.000639
	YPRCENTRY421		yfft_r4dwpn,120K,,320,229504,FMA3_64,I7_64,I7_32
	YPRCENTRY421		yfft_r4dwpn,120K,,192,225920,,RYZEN_64
	DD			0
	PRCSTRT	2619000, 131072, 0.000594
	YPRCENTRY421		yfft_r4dwpn,128K,,10,510208,,I7_32
	YPRCENTRY421		yfft_r4dwpn,128K,,8,278528,,I7_64+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	2855000, 143360, 0.000756
	YPRCENTRY421		yfft_r4dwpn,140K,,320,247168,FMA3_64,,I7
	DD			0
	PRCSTRT	2936000, 147456, 0.000801
	YPRCENTRY421		yfft_r4dwpn,144K,,192,241536,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	3258000, 163840, 0.000805
	YPRCENTRY421		yfft_r4dwpn,160K,,1280,622336,,I7_32
	YPRCENTRY421		yfft_r4dwpn,160K,,320,341376
	YPRCENTRY421		yfft_r4dwpn,160K,,8,298240,,I7_64+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	3412000, 172032, 0.000959
	YPRCENTRY421		yfft_r4dwpn,168K,,192,257152,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	3895000, 196608, 0.000984
	YPRCENTRY421		yfft_r4dwpn,192K,,1536,748800
	YPRCENTRY421		yfft_r4dwpn,192K,,768,454912
	YPRCENTRY421		yfft_r4dwpn,192K,,192,273280
	YPRCENTRY421		yfft_r4dwpn,192K,,8,317952,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	4052000, 204800, 0.001138
	YPRCENTRY421		yfft_r4dwpn,200K,,320,365184,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	4532000, 229376, 0.001229
	YPRCENTRY421		yfft_r4dwpn,224K,,8,337664,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	4842000, 245760, 0.001300
	YPRCENTRY421		yfft_r4dwpn,240K,,768,494080,,I7,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,240K,,320,388992,,RYZEN_64
	DD			0
	PRCSTRT	5164000, 262144, 0.001306
	YPRCENTRY421		yfft_r4dwpn,256K,,11,1012992
	YPRCENTRY421		yfft_r4dwpn,256K,,10,613120,,I7_32+FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,256K,,8,357888,,RYZEN_64
	DD			0
	PRCSTRT	5636000, 286720, 0.001668
	YPRCENTRY421		yfft_r4dwpn,280K,,320,412800,,RYZEN_64
	DD			0
	PRCSTRT	5786000, 294912, 0.001575
	YPRCENTRY421		yfft_r4dwpn,288K,,2304,1113856
	YPRCENTRY421		yfft_r4dwpn,288K,,768,533248,FMA3_64,I7+RYZEN_64
	DD			0
	PRCSTRT	6428000, 327680, 0.001679
	YPRCENTRY421		yfft_r4dwpn,320K,,2560,1240320
	YPRCENTRY421		yfft_r4dwpn,320K,,1280,749824
	YPRCENTRY421		yfft_r4dwpn,320K,,10,664576,FMA3_64,I7,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,320K,,320,437120
	DD			0
	PRCSTRT	6735000, 344064, 0.001969
	YPRCENTRY421		yfft_r4dwpn,336K,,768,572416,FMA3_64,I7,RYZEN_64
	DD			0
	PRCSTRT	7685000, 393216, 0.002047
	YPRCENTRY421		yfft_r4dwpn,384K,,3072,1493248
	YPRCENTRY421		yfft_r4dwpn,384K,,1536,900864
	YPRCENTRY421		yfft_r4dwpn,384K,,10,716032,FMA3_64,I7,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,384K,,768,802816
	DD			0
	PRCSTRT	7998000, 409600, 0.002250
	YPRCENTRY421		yfft_r4dwpn,400K,,1280,813568,FMA3_64,I7,RYZEN_64
	DD			0
	PRCSTRT	8943000, 458752, 0.002446
	YPRCENTRY421		yfft_r4dwpn,448K,,10,767488,FMA3_64,I7+RYZEN_64
	DD			0
	PRCSTRT	9557000, 491520, 0.002709
	YPRCENTRY421		yfft_r4dwpn,480K,,3840,1851136
	YPRCENTRY421		yfft_r4dwpn,480K,,1536,976896
	YPRCENTRY421		yfft_r4dwpn,480K,,1280,877312,FMA3_64,I7,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,480K,,768,855296
	DD			0
	PRCSTRT	10180000, 524288, 0.002712
	YPRCENTRY421		yfft_r4dwpn,512K,,12,2021632,,I7_32
	YPRCENTRY421		yfft_r4dwpn,512K,,11,1214208,,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,512K,,10,1075712,,,RYZEN_64
	DD			0
	PRCSTRT	11120000, 573440, 0.003284
	YPRCENTRY421		yfft_r4dwpn,560K,,1280,941056,FMA3_64,I7,RYZEN_64
	DD			0
	PRCSTRT	11420000, 589824, 0.003313
	YPRCENTRY421		yfft_r4dwpn,576K,,4608,2223360
	YPRCENTRY421		yfft_r4dwpn,576K,,2304,1339648
	YPRCENTRY421		yfft_r4dwpn,576K,,1536,1052928,FMA3_64,I7_32,I7_64
	YPRCENTRY421		yfft_r4dwpn,576K,,768,907776
	DD			0
	PRCSTRT	12670000, 655360, 0.003597
	YPRCENTRY421		yfft_r4dwpn,640K,,5120,2476288
	YPRCENTRY421		yfft_r4dwpn,640K,,2560,1490688
	YPRCENTRY421		yfft_r4dwpn,640K,,11,1314816,FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,640K,,1280,1327104
	YPRCENTRY421		yfft_r4dwpn,640K,,10,1144576,,I7_32,RYZEN_64
	DD			0
	PRCSTRT	13290000, 688128, 0.004010
	YPRCENTRY421		yfft_r4dwpn,672K,,1536,1128960,FMA3_64,I7
	YPRCENTRY421		yfft_r4dwpn,672K,,768,960256,,,RYZEN_64
	DD			0
	PRCSTRT	14170000, 737280, 0.004431
	YPRCENTRY421		yfft_r4dwpn,720K,,2304,1452544,FMA3_64,I7_64
	DD			0
	PRCSTRT	15150000, 786432, 0.004390
	YPRCENTRY421		yfft_r4dwpn,768K,,6144,2979072
	YPRCENTRY421		yfft_r4dwpn,768K,,3072,1792768
	YPRCENTRY421		yfft_r4dwpn,768K,,11,1415424,,FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,768K,,1536,1592832
	YPRCENTRY421		yfft_r4dwpn,768K,,10,1213440,,I7_32,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,768K,,768,1013248
	DD			0
	PRCSTRT	15760000, 819200, 0.004791
	YPRCENTRY421		yfft_r4dwpn,800K,,6400,3079936
	YPRCENTRY421		yfft_r4dwpn,800K,,2560,1615872,,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,800K,,1280,1412352,,I7_32,RYZEN_64
	DD			0
	PRCSTRT	16930000, 884736, 0.005290
	YPRCENTRY421		yfft_r4dwpn,864K,,2304,1565440,,FMA3_64,I7
	DD			0
	PRCSTRT	17620000, 917504, 0.005404
	YPRCENTRY421		yfft_r4dwpn,896K,,11,1516032,,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,896K,,10,1282304,,I7_32,RYZEN_64
	DD			0
	PRCSTRT	18820000, 983040, 0.005806
	YPRCENTRY421		yfft_r4dwpn,960K,,7680,3697920
	YPRCENTRY421		yfft_r4dwpn,960K,,3840,2224384
	YPRCENTRY421		yfft_r4dwpn,960K,,3072,1942528
	YPRCENTRY421		yfft_r4dwpn,960K,,2560,1741056,,I7_32+FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,960K,,1536,1694464
	YPRCENTRY421		yfft_r4dwpn,960K,,1280,1497600,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,960K,,768,1081856
	DD			0
	PRCSTRT	19710000, 1032192, 0.006539
	YPRCENTRY421		yfft_r4dwpn,1008K,,2304,1678336
	DD			0
	PRCSTRT	20080000, 1048576, 0.005906
	YPRCENTRY421		yfft_r4dwpn,1M,,13,4035840
	YPRCENTRY421		yfft_r4dwpn,1M,,12,2419456,I7_64,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1M,,11,2135552
	YPRCENTRY421		yfft_r4dwpn,1M,,10,1351680,,I7_32,RYZEN_64
	DD			0
	PRCSTRT	21910000, 1146880, 0.007069
	YPRCENTRY421		yfft_r4dwpn,1120K,,2560,1866240,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1120K,,1280,1582848,,I7_32,RYZEN_64
	DD			0
	PRCSTRT	22500000, 1179648, 0.007022
	YPRCENTRY421		yfft_r4dwpn,1152K,,9216,4442368
	YPRCENTRY421		yfft_r4dwpn,1152K,,4608,2670336
	YPRCENTRY421		yfft_r4dwpn,1152K,,3072,2092288,,I7_32+FMA3_64+RYZEN_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,1152K,,2304,2375680
	YPRCENTRY421		yfft_r4dwpn,1152K,,1536,1796096
	DD			0
	PRCSTRT	23430000, 1228800, 0.007757
	YPRCENTRY421		yfft_r4dwpn,1200K,,3840,2411008,,I7_32+FMA3_64
	DD			0
	PRCSTRT	24990000, 1310720, 0.007707
	YPRCENTRY421		yfft_r4dwpn,1280K,,10240,4945152
	YPRCENTRY421		yfft_r4dwpn,1280K,,5120,2972416
	YPRCENTRY421		yfft_r4dwpn,1280K,,12,2618368,I7_64,I7_32+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1280K,,2560,2641408
	YPRCENTRY421		yfft_r4dwpn,1280K,,11,2269952
	YPRCENTRY421		yfft_r4dwpn,1280K,,1280,1668608
	YPRCENTRY421		yfft_r4dwpn,1280K,,10,1436672,,,RYZEN_64
	DD			0
	PRCSTRT	26190000, 1376256, 0.008592
	YPRCENTRY421		yfft_r4dwpn,1344K,,3072,2242048,,I7_64+FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1344K,,1536,1897728,,I7_32
	DD			0
	PRCSTRT	27990000, 1474560, 0.009242
	YPRCENTRY421		yfft_r4dwpn,1440K,,4608,2893824,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1440K,,3840,2597632,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1440K,,2304,2526464,,I7_64
	DD			0
	PRCSTRT	29870000, 1572864, 0.009264
	YPRCENTRY421		yfft_r4dwpn,1536K,,12288,5953792
	YPRCENTRY421		yfft_r4dwpn,1536K,,6144,3573504
	YPRCENTRY421		yfft_r4dwpn,1536K,,12,2817280,I7_64,I7_32+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1536K,,3072,3172864
	YPRCENTRY421		yfft_r4dwpn,1536K,,11,2404352
	YPRCENTRY421		yfft_r4dwpn,1536K,,1536,1999872
	YPRCENTRY421		yfft_r4dwpn,1536K,,10,1642496,,,RYZEN_64
	DD			0
	PRCSTRT	31050000, 1638400, 0.010048
	YPRCENTRY421		yfft_r4dwpn,1600K,,12800,6155520
	YPRCENTRY421		yfft_r4dwpn,1600K,,6400,3698944
	YPRCENTRY421		yfft_r4dwpn,1600K,,5120,3220480,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1600K,,2560,2808576,,I7_64,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1600K,,1280,1769984
	DD			0
	PRCSTRT	32560000, 1720320, 0.011459
	YPRCENTRY421		yfft_r4dwpn,1680K,,3840,2784256,,I7_64+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	33440000, 1769472, 0.011102
	YPRCENTRY421		yfft_r4dwpn,1728K,,4608,3117312,,,I7_32
	YPRCENTRY421		yfft_r4dwpn,1728K,,2304,2677248,,I7_64,FMA3_64
	DD			0
	PRCSTRT	34740000, 1835008, 0.011146
	YPRCENTRY421		yfft_r4dwpn,1792K,,12,3016192,I7_64,I7_32,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1792K,,11,2538752
	YPRCENTRY421		yfft_r4dwpn,1792K,,10,1848320,,,RYZEN_64
	DD			0
	PRCSTRT	37120000, 1966080, 0.012046
	YPRCENTRY421		yfft_r4dwpn,1920K,,15360,7391488
	YPRCENTRY421		yfft_r4dwpn,1920K,,7680,4439808
	YPRCENTRY421		yfft_r4dwpn,1920K,,6144,3870720,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1920K,,5120,3468544,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1920K,,3840,3948544
	YPRCENTRY421		yfft_r4dwpn,1920K,,3072,3372800
	YPRCENTRY421		yfft_r4dwpn,1920K,,2560,2975744,,I7_64
	YPRCENTRY421		yfft_r4dwpn,1920K,,1536,2117632
	DD			0
	PRCSTRT	38650000, 2048000, 0.013344
	YPRCENTRY421		yfft_r4dwpn,2000K,,6400,4008448
	DD			0
	PRCSTRT	38900000, 2064384, 0.013399
	YPRCENTRY421		yfft_r4dwpn,2016K,,4608,3340800,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2016K,,2304,2828032
	DD			0
	PRCSTRT	39560000, 2097152, 0.012290
	YPRCENTRY421		yfft_r4dwpn,2M,,14,8067328
	YPRCENTRY421		yfft_r4dwpn,2M,,13,4826880
	YPRCENTRY421		yfft_r4dwpn,2M,,12,4258304,I7_64,I7_32,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2M,,11,2673664
	YPRCENTRY421		yfft_r4dwpn,2M,,10,2482176,,,RYZEN_64
	DD			0
	PRCSTRT	43180000, 2293760, 0.014710
	YPRCENTRY421		yfft_r4dwpn,2240K,,5120,3716608
	YPRCENTRY421		yfft_r4dwpn,2240K,,2560,3142912,,,RYZEN_64
	DD			0
	PRCSTRT	44390000, 2359296, 0.014635
	YPRCENTRY421		yfft_r4dwpn,2304K,,9216,5331712
	YPRCENTRY421		yfft_r4dwpn,2304K,,6144,4167936,,I7_32+FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,2304K,,4608,4738560
	YPRCENTRY421		yfft_r4dwpn,2304K,,3072,3572736,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,2304K,,2304,2979328
	DD			0
	PRCSTRT	46180000, 2457600, 0.015855
	YPRCENTRY421		yfft_r4dwpn,2400K,,7680,4810752,,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2400K,,6400,4317952,,I7_32
	YPRCENTRY421		yfft_r4dwpn,2400K,,3840,4197632,,,RYZEN_64
	DD			0
	PRCSTRT	49250000, 2621440, 0.016017
	YPRCENTRY421		yfft_r4dwpn,2560K,,20480,9885952
	YPRCENTRY421		yfft_r4dwpn,2560K,,10240,5932800
	YPRCENTRY421		yfft_r4dwpn,2560K,,13,5222400
	YPRCENTRY421		yfft_r4dwpn,2560K,,5120,5270016
	YPRCENTRY421		yfft_r4dwpn,2560K,,12,4523776,,I7_32,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2560K,,2560,3310592
	YPRCENTRY421		yfft_r4dwpn,2560K,,11,2824192,,,RYZEN_64
	DD			0
	PRCSTRT	51630000, 2752512, 0.017883
	YPRCENTRY421		yfft_r4dwpn,2688K,,6144,4465152,,I7+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2688K,,3072,3772672,,,RYZEN_64
	DD			0
	PRCSTRT	53710000, 2867200, 0.019289
	YPRCENTRY421		yfft_r4dwpn,2800K,,6400,4627456,,I7_32,RYZEN_64
	DD			0
	PRCSTRT	55190000, 2949120, 0.019245
	YPRCENTRY421		yfft_r4dwpn,2880K,,9216,5776384
	YPRCENTRY421		yfft_r4dwpn,2880K,,7680,5181696,,I7_32+FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,2880K,,4608,5036800
	YPRCENTRY421		yfft_r4dwpn,2880K,,3840,4446720,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,2880K,,2304,3146240
	DD			0
	PRCSTRT	58870000, 3145728, 0.019430
	YPRCENTRY421		yfft_r4dwpn,3M,,12288,7138048
	YPRCENTRY421		yfft_r4dwpn,3M,,13,5617920
	YPRCENTRY421		yfft_r4dwpn,3M,,6144,6329856
	YPRCENTRY421		yfft_r4dwpn,3M,,12,4789248,,I7_32+RYZEN_64,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3M,,3072,3973120
	YPRCENTRY421		yfft_r4dwpn,3M,,11,3226624
	DD			0
	PRCSTRT	61200000, 3276800, 0.020978
	YPRCENTRY421		yfft_r4dwpn,3200K,,25600,12306688
	YPRCENTRY421		yfft_r4dwpn,3200K,,12800,7388928,,,I7_64
	YPRCENTRY421		yfft_r4dwpn,3200K,,10240,6426624,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3200K,,6400,6569984
	YPRCENTRY421		yfft_r4dwpn,3200K,,5120,5601024,,I7_32
	YPRCENTRY421		yfft_r4dwpn,3200K,,2560,3493888,,,RYZEN_64
	DD			0
	PRCSTRT	64210000, 3440640, 0.023358
	YPRCENTRY421		yfft_r4dwpn,3360K,,7680,5552640,,I7_64
	YPRCENTRY421		yfft_r4dwpn,3360K,,3840,4695808,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	65860000, 3538944, 0.023145
	YPRCENTRY421		yfft_r4dwpn,3456K,,9216,6221056,,I7
	YPRCENTRY421		yfft_r4dwpn,3456K,,4608,5335040,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3456K,,2304,3597824
	DD			0
	PRCSTRT	68480000, 3670016, 0.024079
	YPRCENTRY421		yfft_r4dwpn,3584K,,13,6013440,,I7_32
	YPRCENTRY421		yfft_r4dwpn,3584K,,12,5054720,,RYZEN_64,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3584K,,11,3629056
	DD			0
	PRCSTRT	73160000, 3932160, 0.025515
	YPRCENTRY421		yfft_r4dwpn,3840K,,15360,8870656,,I7_64
	YPRCENTRY421		yfft_r4dwpn,3840K,,12288,7730176
	YPRCENTRY421		yfft_r4dwpn,3840K,,10240,6920448
	YPRCENTRY421		yfft_r4dwpn,3840K,,7680,7884288
	YPRCENTRY421		yfft_r4dwpn,3840K,,6144,6726400,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3840K,,5120,5932032,,I7_32
	YPRCENTRY421		yfft_r4dwpn,3840K,,3840,4945408
	YPRCENTRY421		yfft_r4dwpn,3840K,,3072,4189184,,,RYZEN_64
	DD			0
	PRCSTRT	76180000, 4096000, 0.028093
	YPRCENTRY421		yfft_r4dwpn,4000K,,12800,8005632
	YPRCENTRY421		yfft_r4dwpn,4000K,,6400,6982912
	DD			0
	PRCSTRT	76660000, 4128768, 0.027915
	YPRCENTRY421		yfft_r4dwpn,4032K,,9216,6665728,,I7_64
	YPRCENTRY421		yfft_r4dwpn,4032K,,4608,5633280
	YPRCENTRY421		yfft_r4dwpn,4032K,,2304,4049408
	DD			0
	PRCSTRT	78040000, 4194304, 0.026160
	YPRCENTRY421		yfft_r4dwpn,4M,,14,9644800,I7_64
	YPRCENTRY421		yfft_r4dwpn,4M,,13,8500736,,I7_32
	YPRCENTRY421		yfft_r4dwpn,4M,,12,5320704,,RYZEN_64,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,4M,,11,4918272
	DD			0
	PRCSTRT	85170000, 4587520, 0.031930
	YPRCENTRY421		yfft_r4dwpn,4480K,,10240,7414272,RYZEN_64,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,4480K,,5120,6263040,,I7_32
	DD			0
	PRCSTRT	87450000, 4718592, 0.030538
	YPRCENTRY421		yfft_r4dwpn,4608K,,12288,8322304,,I7_64
	YPRCENTRY421		yfft_r4dwpn,4608K,,9216,9464320
	YPRCENTRY421		yfft_r4dwpn,4608K,,6144,7122944,,I7_32,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,4608K,,4608,5932032
	YPRCENTRY421		yfft_r4dwpn,4608K,,3072,4788224,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,4608K,,2304,5502464
	DD			0
	PRCSTRT	91050000, 4915200, 0.033479
	YPRCENTRY421		yfft_r4dwpn,4800K,,15360,9610240
	YPRCENTRY421		yfft_r4dwpn,4800K,,12800,8622336
	YPRCENTRY421		yfft_r4dwpn,4800K,,7680,8379136,,I7_64+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,4800K,,6400,7395840,,I7_32
	YPRCENTRY421		yfft_r4dwpn,4800K,,3840,5210624,,,RYZEN_64
	DD			0
	PRCSTRT	97090000, 5242880, 0.034976
	YPRCENTRY421		yfft_r4dwpn,5M,,20480,11856640
	YPRCENTRY421		yfft_r4dwpn,5M,,14,10433536,I7_64
	YPRCENTRY421		yfft_r4dwpn,5M,,10240,10524160
	YPRCENTRY421		yfft_r4dwpn,5M,,13,9028352,,I7_32
	YPRCENTRY421		yfft_r4dwpn,5M,,5120,6594560
	YPRCENTRY421		yfft_r4dwpn,5M,,12,5602304,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	101700000, 5505024, 0.038004
	YPRCENTRY421		yfft_r4dwpn,5376K,,12288,8914432,,I7_64
	YPRCENTRY421		yfft_r4dwpn,5376K,,6144,7519488,,I7_32+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,5376K,,3072,5387264,,,RYZEN_64
	DD			0
	PRCSTRT	105900000, 5734400, 0.040691
	YPRCENTRY421		yfft_r4dwpn,5600K,,12800,9239040
	YPRCENTRY421		yfft_r4dwpn,5600K,,6400,7808768,,I7_32+RYZEN_64
	DD			0
	PRCSTRT	108700000, 5898240, 0.040324
	YPRCENTRY421		yfft_r4dwpn,5760K,,15360,10349824,,I7_64,I7_32
	YPRCENTRY421		yfft_r4dwpn,5760K,,9216,10057472
	YPRCENTRY421		yfft_r4dwpn,5760K,,7680,8873984,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,5760K,,4608,6246400
	YPRCENTRY421		yfft_r4dwpn,5760K,,3840,5957120,,FMA3_64
	DD			0
	PRCSTRT	115900000, 6291456, 0.041752
	YPRCENTRY421		yfft_r4dwpn,6M,,14,11222272,I7_64
	YPRCENTRY421		yfft_r4dwpn,6M,,12288,12646912
	YPRCENTRY421		yfft_r4dwpn,6M,,13,9555968,,I7_32
	YPRCENTRY421		yfft_r4dwpn,6M,,6144,7916544
	YPRCENTRY421		yfft_r4dwpn,6M,,12,6397952,,FMA3_64,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,6M,,3072,7331840
	DD			0
	PRCSTRT	120600000, 6553600, 0.044473
	YPRCENTRY421		yfft_r4dwpn,6400K,,25600,14768896
	YPRCENTRY421		yfft_r4dwpn,6400K,,20480,12841984
	YPRCENTRY421		yfft_r4dwpn,6400K,,12800,13127168
	YPRCENTRY421		yfft_r4dwpn,6400K,,10240,11182848,,I7+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,6400K,,6400,8222208
	YPRCENTRY421		yfft_r4dwpn,6400K,,5120,6941696,,,RYZEN_64
	DD			0
	PRCSTRT	126400000, 6881280, 0.048606
	YPRCENTRY421		yfft_r4dwpn,6720K,,15360,11089408,,I7
	YPRCENTRY421		yfft_r4dwpn,6720K,,7680,9368832,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,6720K,,3840,6703616,,,FMA3_64
	DD			0
	PRCSTRT	129800000, 7077888, 0.050847
	YPRCENTRY421		yfft_r4dwpn,6912K,,9216,10650624,,I7
	YPRCENTRY421		yfft_r4dwpn,6912K,,4608,7140352,,FMA3_64
	DD			0
	PRCSTRT	134900000, 7340032, 0.052097
	YPRCENTRY421		yfft_r4dwpn,7M,,14,12011008,I7_64
	YPRCENTRY421		yfft_r4dwpn,7M,,13,10083584,,I7_32
	YPRCENTRY421		yfft_r4dwpn,7M,,12,7193600,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	144100000, 7864320, 0.053362
	YPRCENTRY421		yfft_r4dwpn,7680K,,20480,13827328
	YPRCENTRY421		yfft_r4dwpn,7680K,,15360,15755776
	YPRCENTRY421		yfft_r4dwpn,7680K,,12288,13436672
	YPRCENTRY421		yfft_r4dwpn,7680K,,10240,11841536,RYZEN_64,I7
	YPRCENTRY421		yfft_r4dwpn,7680K,,7680,9864192
	YPRCENTRY421		yfft_r4dwpn,7680K,,6144,8329216,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,7680K,,5120,7933952
	YPRCENTRY421		yfft_r4dwpn,7680K,,3840,9139712
	DD			0
	PRCSTRT	149900000, 8192000, 0.059014
	YPRCENTRY421		yfft_r4dwpn,8000K,,25600,16000000
	YPRCENTRY421		yfft_r4dwpn,8000K,,12800,13949696,,I7_64
	YPRCENTRY421		yfft_r4dwpn,8000K,,6400,8651264,,,RYZEN_64
	DD			0
	PRCSTRT	151000000, 8257536, 0.063072
	YPRCENTRY421		yfft_r4dwpn,8064K,,9216,11243776
	YPRCENTRY421		yfft_r4dwpn,8064K,,4608,8034304,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	153500000, 8388608, 0.058218
	YPRCENTRY421		yfft_r4dwpn,8M,,14,16988672,I7_64
	YPRCENTRY421		yfft_r4dwpn,8M,,13,10611712,,I7_32
	YPRCENTRY421		yfft_r4dwpn,8M,,12,9793536,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	167700000, 9175040, 0.068512
	YPRCENTRY421		yfft_r4dwpn,8960K,,20480,14812672,,I7_64
	YPRCENTRY421		yfft_r4dwpn,8960K,,10240,12500224,RYZEN_64,I7_32+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,8960K,,5120,8926208
	DD			0
	PRCSTRT	172300000, 9437184, 0.067985
	YPRCENTRY421		yfft_r4dwpn,9M,,12288,14226432,RYZEN_64,I7
	YPRCENTRY421		yfft_r4dwpn,9M,,9216,11837440
	YPRCENTRY421		yfft_r4dwpn,9M,,6144,9518080,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,9M,,4608,10961920
	DD			0
	PRCSTRT	179200000, 9830400, 0.070403
	YPRCENTRY421		yfft_r4dwpn,9600K,,25600,17231104,,,I7_32
	YPRCENTRY421		yfft_r4dwpn,9600K,,15360,16742144,,I7_64
	YPRCENTRY421		yfft_r4dwpn,9600K,,12800,14772224
	YPRCENTRY421		yfft_r4dwpn,9600K,,7680,10375168,RYZEN_64,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,9600K,,6400,9889280
	DD			0
	PRCSTRT	190700000, 10485760, 0.076827
	YPRCENTRY421		yfft_r4dwpn,10M,,20480,21035520
	YPRCENTRY421		yfft_r4dwpn,10M,,14,18040576
	YPRCENTRY421		yfft_r4dwpn,10M,,10240,13159424,RYZEN_64,I7,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,10M,,13,11155456
	YPRCENTRY421		yfft_r4dwpn,10M,,5120,12181504
	DD			0
	PRCSTRT	200300000, 11010048, 0.083177
	YPRCENTRY421		yfft_r4dwpn,10752K,,12288,15016192,,I7
	YPRCENTRY421		yfft_r4dwpn,10752K,,6144,10706944,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	208400000, 11468800, 0.085770
	YPRCENTRY421		yfft_r4dwpn,11200K,,25600,18462208
	YPRCENTRY421		yfft_r4dwpn,11200K,,12800,15594752,,I7_32
	YPRCENTRY421		yfft_r4dwpn,11200K,,6400,11127296,,,RYZEN_64
	DD			0
	PRCSTRT	214100000, 11796480, 0.086096
	YPRCENTRY421		yfft_r4dwpn,11520K,,15360,17728512,,I7
	YPRCENTRY421		yfft_r4dwpn,11520K,,9216,12446720
	YPRCENTRY421		yfft_r4dwpn,11520K,,7680,11858944,RYZEN_64,,FMA3_64
	DD			0
	PRCSTRT	228400000, 12582912, 0.094441
	YPRCENTRY421		yfft_r4dwpn,12M,,14,19092480
	YPRCENTRY421		yfft_r4dwpn,12M,,12288,15806464,,RYZEN_64,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,12M,,13,12737536,I7_64
	YPRCENTRY421		yfft_r4dwpn,12M,,6144,14617600,I7_32
	DD			0
	PRCSTRT	237400000, 13107200, 0.094766
	YPRCENTRY421		yfft_r4dwpn,12800K,,25600,26241536,,I7
	YPRCENTRY421		yfft_r4dwpn,12800K,,20480,22349568
	YPRCENTRY421		yfft_r4dwpn,12800K,,12800,16417792
	YPRCENTRY421		yfft_r4dwpn,12800K,,10240,13834240,RYZEN_64,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,12800K,,6400,15201792
	DD			0
	PRCSTRT	249100000, 13762560, 0.106540
	YPRCENTRY421		yfft_r4dwpn,13440K,,15360,18714880,,I7
	YPRCENTRY421		yfft_r4dwpn,13440K,,7680,13342720,RYZEN_64,,FMA3_64
	DD			0
	PRCSTRT	255400000, 14155776, 0.111224
	YPRCENTRY421		yfft_r4dwpn,13824K,,9216,14225408,I7+RYZEN_64,,FMA3_64
	DD			0
	PRCSTRT	265600000, 14680064, 0.114093
	YPRCENTRY421		yfft_r4dwpn,14M,,14,20144384
	YPRCENTRY421		yfft_r4dwpn,14M,,13,14319616,I7+RYZEN_64,,FMA3_64
	DD			0
	PRCSTRT	284000000, 15728640, 0.120880
	YPRCENTRY421		yfft_r4dwpn,15M,,20480,23663616
	YPRCENTRY421		yfft_r4dwpn,15M,,15360,19701760
	YPRCENTRY421		yfft_r4dwpn,15M,,12288,16612352
	YPRCENTRY421		yfft_r4dwpn,15M,,10240,15809536,I7+RYZEN_64,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,15M,,7680,18236416
	DD			0
	PRCSTRT	295500000, 16384000, 0.127481
	YPRCENTRY421		yfft_r4dwpn,16000K,,25600,27883264,,I7
	YPRCENTRY421		yfft_r4dwpn,16000K,,12800,17256448,,RYZEN_64,FMA3_64
	DD			0
	PRCSTRT	297400000, 16515072, 0.132344
	YPRCENTRY421		yfft_r4dwpn,16128K,,9216,16004096,I7_64+RYZEN_64
	DD			0
	PRCSTRT	302600000, 16777216, 0.129355
	YPRCENTRY421		yfft_r4dwpn,16M,,14,21196800,I7_64
	YPRCENTRY421		yfft_r4dwpn,16M,,13,19540992,I7_32+RYZEN_64,,FMA3_64
	DD			0
	PRCSTRT	330200000, 18350080, 0.147854
	YPRCENTRY421		yfft_r4dwpn,17920K,,20480,24977664
	YPRCENTRY421		yfft_r4dwpn,17920K,,10240,17784832,I7,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	339000000, 18874368, 0.149016
	YPRCENTRY421		yfft_r4dwpn,18M,,12288,18980864,I7_64,FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,18M,,9216,21880832,I7_32
	DD			0
	PRCSTRT	353100000, 19660800, 0.153223
	YPRCENTRY421		yfft_r4dwpn,19200K,,25600,29524992,,I7
	YPRCENTRY421		yfft_r4dwpn,19200K,,15360,20704256
	YPRCENTRY421		yfft_r4dwpn,19200K,,12800,19723264,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	376200000, 20971520, 0.168125
	YPRCENTRY421		yfft_r4dwpn,20M,,20480,26292224
	YPRCENTRY421		yfft_r4dwpn,20M,,14,22264832,I7_64
	YPRCENTRY421		yfft_r4dwpn,20M,,10240,24316928,I7_32,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	394400000, 22020096, 0.179899
	YPRCENTRY421		yfft_r4dwpn,21M,,12288,21349376,I7,RYZEN_64,FMA3_64
	DD			0
	PRCSTRT	410600000, 22937600, 0.186258
	YPRCENTRY421		yfft_r4dwpn,22400K,,25600,31166720,,I7_32
	YPRCENTRY421		yfft_r4dwpn,22400K,,12800,22190080,I7_64,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	421200000, 23592960, 0.194140
	YPRCENTRY421		yfft_r4dwpn,23040K,,15360,23662592,I7,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	449500000, 25165824, 0.203308
	YPRCENTRY421		yfft_r4dwpn,24M,,14,25419776,I7_64
	YPRCENTRY421		yfft_r4dwpn,24M,,12288,29192192,I7_32,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	467500000, 26214400, 0.208601
	YPRCENTRY421		yfft_r4dwpn,25M,,25600,32808960,I7
	YPRCENTRY421		yfft_r4dwpn,25M,,20480,27622400
	YPRCENTRY421		yfft_r4dwpn,25M,,12800,30360576,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	489900000, 27525120, 0.229014
	YPRCENTRY421		yfft_r4dwpn,26880K,,15360,26620928,I7,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	522900000, 29360128, 0.247130
	YPRCENTRY421		yfft_r4dwpn,28M,,14,28574720,I7,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	558500000, 31457280, 0.258182
	YPRCENTRY421		yfft_r4dwpn,30M,,20480,31563776,,FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,30M,,15360,36429824,I7
	DD			0
	PRCSTRT	581300000, 32768000, 0.271767
	YPRCENTRY421		yfft_r4dwpn,32000K,,25600,34466816,I7,RYZEN_64
	DD			0
	PRCSTRT	595700000, 33554432, 0.279123
	YPRCENTRY421		yfft_r4dwpn,32M,,14,39038976,I7,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	649800000, 36700160, 0.306
	YPRCENTRY421		yfft_r4dwpn,35M,,20480,35505152,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	694900000, 39321600, 0.328
	YPRCENTRY421		yfft_r4dwpn,38400K,,25600,39391232,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	740200000, 41943040, 0.350
	YPRCENTRY421		yfft_r4dwpn,40M,,20480,48590848,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	808200000, 45875200, 0.390
	YPRCENTRY421		yfft_r4dwpn,44800K,,25600,44315648,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	920800000, 52428800, 0.440
	YPRCENTRY421		yfft_r4dwpn,50M,,25600,60678144,,,FMA3_64+RYZEN_64
	DD			0
	DD	0

yjmptablep DD	0
	org	$-4
	PRCSTRT	739,	32,	0.00000015
	PRCENTRY2		yfft_r4_32_op_ac, 776, I7 + FMA3_64
	DD			0
	PRCSTRT	1465,	64,	0.0000002
	PRCENTRY2		yfft_r4_64_op_ac, 1552, I7 + FMA3_64
	DD			0
	PRCSTRT	2173,	96,	0.0000003
	PRCENTRY2		yfft_r4_96_op_ac, 2392, I7 + FMA3_64
	DD			0
	PRCSTRT	2897,	128,	0.00000035
	PRCENTRY2		yfft_r4_128_op_ac, 3232, I7 + FMA3_64
	DD			0
	PRCSTRT	3605,	160,	0.0000004
	PRCENTRY2		yfft_r4_160_op_ac, 4008, I7 + FMA3_64
	DD			0
	PRCSTRT	4295,	192,	0.0000006
	PRCENTRY2		yfft_r4_192_op_ac, 4784, I7 + FMA3_64
	DD			0
	PRCSTRT	5729,	256,	0.0000007
	PRCENTRY2		yfft_r4_256_op_ac, 6464, I7 + FMA3_64
	DD			0
	PRCSTRT	7111,	320,	0.0000010
	PRCENTRY2		yfft_r4_320_op_ac, 8016, I7 + FMA3_64
	DD			0
	PRCSTRT	8493,	384,	0.0000012
	PRCENTRY2		yfft_r4_384_op_ac, 9696, I7 + FMA3_64
	DD			0
	PRCSTRT	11319,	512,	0.0000015
	PRCENTRY2		yfft_r4_512_op_ac, 12928, I7 + FMA3_64
	DD			0
	PRCSTRT	14049,	640,	0.0000021
	PRCENTRY2		yfft_r4_640_op_ac, 16160, I7 + FMA3_64
	DD			0
	PRCSTRT	16779,	768,	0.0000028
	PRCENTRY2		yfft_r4_768_op_ac, 19392, I7 + FMA3_64
	DD			0
	PRCSTRT	22431,	1024,	0.0000036
	PRCENTRY2		yfft_r4_1K_op_ac, 25856, I7 + FMA3_64
	DD			0
	PRCSTRT	27863,	1280,	0.0000050
	PRCENTRY2		yfft_r4_1280_op_ac, 32320, I7 + FMA3_64
	DD			0
	PRCSTRT	33313,	1536,	0.0000068
	PRCENTRY2		yfft_r4_1536_op_ac, 38784, I7 + FMA3_64
	DD			0
	PRCSTRT	44303,	2048,	0.0000083
	PRCENTRY2		yfft_r4_2K_op_ac, 51712, I7 + FMA3_64
	DD			0
	PRCSTRT	55043,	2560,	0.000012
	PRCENTRY2		yfft_r4_2560_op_ac, 64640, I7 + FMA3_64
	DD			0
	PRCSTRT	65767,	3072,	0.000015
	PRCENTRY2		yfft_r4_3K_op_ac, 77568, I7 + FMA3_64
	DD			0
	PRCSTRT	87375,	4096,	0.000019
	PRCENTRY2		yfft_r4_4K_op_ac, 103424, I7 + FMA3_64
	DD			0
	PRCSTRT	108600,	5120,	0.000025
	PRCENTRY2		yfft_r4_5K_op_ac, 129280, I7 + FMA3_64
	DD			0
	PRCSTRT	129700,	6144,	0.000025
	PRCENTRY2		yfft_r4_6K_op_ac, 155136, I7_32
	YPRCENTRY421		yfft_r4dwpn,6K,_ac,48,20992,,I7_64+FMA3_64
	DD			0
	PRCSTRT	172600,	8192,	0.000033
	PRCENTRY2		yfft_r4_8K_op_ac, 206848, I7_32
	YPRCENTRY421		yfft_r4dwpn,8K,_ac,6,26240,,I7_64+FMA3_64
	DD			0
	PRCSTRT	214700,	10240,	0.000046
	PRCENTRY2		yfft_r4_10K_op_ac, 258560, I7_32
	YPRCENTRY421		yfft_r4dwpn,10K,_ac,80,31488,,I7_64+FMA3_64
	DD			0
	PRCSTRT	256500,	12288,	0.000054
	PRCENTRY2		yfft_r4_12K_op_ac, 310272, I7_32
	YPRCENTRY421		yfft_r4dwpn,12K,_ac,48,35840,,I7_64+FMA3_64
	DD			0
	PRCSTRT	340400,	16384,	0.000073
	PRCENTRY2		yfft_r4_16K_op_ac, 413696, I7_32
	YPRCENTRY421		yfft_r4dwpn,16K,_ac,6,45184,,I7_64+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	380200,	18432,	0.000088
	PRCENTRY2		yfft_r4_18K_op_ac, 461312
	YPRCENTRY421		yfft_r4dwpn,18K,_ac,48,40704,,I7_64+FMA3_64+RYZEN_64,I7_32
	DD			0
	PRCSTRT	423300,	20480,	0.000095
	PRCENTRY2		yfft_r4_20K_op_ac, 517120
	YPRCENTRY421		yfft_r4dwpn,20K,_ac,80,54528,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	505500,	24576,	0.000109
	PRCENTRY2		yfft_r4_24K_op_ac, 620544
	YPRCENTRY421		yfft_r4dwpn,24K,_ac,192,68224,RYZEN_64,I7
	YPRCENTRY421		yfft_r4dwpn,24K,_ac,6,49024
	YPRCENTRY421		yfft_r4dwpn,24K,_ac,48,50688,,FMA3_64
	DD			0
	PRCSTRT	628900,	30720,	0.000152
	YPRCENTRY421		yfft_r4dwpn,30K,_ac,80,57344,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,30K,_ac,48,60416,,FMA3_64
	DD			0
	PRCSTRT	672100,	32768,	0.000140
	PRCENTRY2		yfft_r4_32K_op_ac,827392
	YPRCENTRY421		yfft_r4dwpn,32K,_ac,8,89600,,I7+FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,32K,_ac,6,60544
	DD			0
	PRCSTRT	751200,	36864,	0.000180
	YPRCENTRY421		yfft_r4dwpn,36K,_ac,48,64000,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	835600, 40960,	0.000191
	YPRCENTRY421		yfft_r4dwpn,40K,_ac,320,110208,RYZEN_64,I7
	YPRCENTRY421		yfft_r4dwpn,40K,_ac,80,70400,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,40K,_ac,6,71808
	DD			0
	PRCSTRT	998700, 49152,	0.000240
	YPRCENTRY421		yfft_r4dwpn,48K,_ac,192,119936,,I7+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,48K,_ac,6,77440,,FMA3_64
	DD			0
	PRCSTRT	1039000, 51200,	0.000272
	YPRCENTRY421		yfft_r4dwpn,50K,_ac,80,83200,FMA3_64,I7_64+RYZEN_64,I7_32
	DD			0
	PRCSTRT	1241000, 61440,	0.000333
	YPRCENTRY421		yfft_r4dwpn,60K,_ac,80,90880,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	1326000, 65536,	0.000314
	YPRCENTRY421		yfft_r4dwpn,64K,_ac,8,157696,,I7_64+FMA3_64+RYZEN_64,I7_32
	DD			0
	PRCSTRT	1481000, 73728,	0.000374
	YPRCENTRY421		yfft_r4dwpn,72K,_ac,192,115584,FMA3_64,I7+RYZEN_64
	DD			0
	PRCSTRT	1649000, 81920,	0.000413
	YPRCENTRY421		yfft_r4dwpn,80K,_ac,320,194688,FMA3_64,I7+RYZEN_64
	DD			0
	PRCSTRT	1969000, 98304, 0.000479
	YPRCENTRY421		yfft_r4dwpn,96K,_ac,768,257536,,I7_64
	YPRCENTRY421		yfft_r4dwpn,96K,_ac,8,149248,FMA3_64,,I7_32
	YPRCENTRY421		yfft_r4dwpn,96K,_ac,192,139392,,RYZEN_64
	DD			0
	PRCSTRT	2448000, 122880, 0.000652
	YPRCENTRY421		yfft_r4dwpn,120K,_ac,320,182144,FMA3_64,I7_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,120K,_ac,192,162944
	DD			0
	PRCSTRT	2615000, 131072, 0.000622
	YPRCENTRY421		yfft_r4dwpn,128K,_ac,10,343040,,I7
	YPRCENTRY421		yfft_r4dwpn,128K,_ac,8,179200,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	2927000, 147456, 0.000841
	YPRCENTRY421		yfft_r4dwpn,144K,_ac,192,184960,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	3252000, 163840, 0.000839
	YPRCENTRY421		yfft_r4dwpn,160K,_ac,1280,425472,,I7_64
	YPRCENTRY421		yfft_r4dwpn,160K,_ac,320,218240,,RYZEN_64,I7_32
	YPRCENTRY421		yfft_r4dwpn,160K,_ac,8,208896,FMA3_64
	DD	0
	PRCSTRT	3886000, 196608, 0.001030
	YPRCENTRY421		yfft_r4dwpn,192K,_ac,1536,509952,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,192K,_ac,768,456704,,I7_64
	YPRCENTRY421		yfft_r4dwpn,192K,_ac,8,239104,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,192K,_ac,192,217728
	DD	0
	PRCSTRT	4045000, 204800, 0.001191
	YPRCENTRY421		yfft_r4dwpn,200K,_ac,320,254080,FMA3_64,I7_64+RYZEN_64,I7_32
	DD	0
	PRCSTRT	4835000, 245760, 0.001449
	YPRCENTRY421		yfft_r4dwpn,240K,_ac,320,292480
	YPRCENTRY421		yfft_r4dwpn,240K,_ac,192,249984,,I7_64+FMA3_64+RYZEN_64
	DD	0
	PRCSTRT	5158000, 262144, 0.001365
	YPRCENTRY421		yfft_r4dwpn,256K,_ac,11,677888,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,256K,_ac,10,607744,,I7
	YPRCENTRY421		yfft_r4dwpn,256K,_ac,8,280064,,,RYZEN_64
	DD	0
	PRCSTRT	5765000, 294912, 0.001632
	YPRCENTRY421		yfft_r4dwpn,288K,_ac,2304,761344
	YPRCENTRY421		yfft_r4dwpn,288K,_ac,768,415488,RYZEN_64,I7+FMA3_64
	DD			0
	PRCSTRT	6418000, 327680, 0.001802
	YPRCENTRY421		yfft_r4dwpn,320K,_ac,2560,845824,FMA3_64,I7_32
	YPRCENTRY421		yfft_r4dwpn,320K,_ac,1280,755712,,I7_64
	YPRCENTRY421		yfft_r4dwpn,320K,_ac,320,341632
	YPRCENTRY421		yfft_r4dwpn,320K,_ac,8,320512,,,RYZEN_64
	DD			0
	PRCSTRT	7665000, 393216, 0.002134
	YPRCENTRY421		yfft_r4dwpn,384K,_ac,3072,1014784
	YPRCENTRY421		yfft_r4dwpn,384K,_ac,1536,905728
	YPRCENTRY421		yfft_r4dwpn,384K,_ac,10,550144,FMA3_64,I7,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,384K,_ac,768,494592
	DD			0
	PRCSTRT	7982000, 409600, 0.002654
	YPRCENTRY421		yfft_r4dwpn,400K,_ac,320,390272,,I7+RYZEN_64
	DD			0
	PRCSTRT	9526000, 491520, 0.002815
	YPRCENTRY421		yfft_r4dwpn,480K,_ac,3840,1265152
	YPRCENTRY421		yfft_r4dwpn,480K,_ac,1280,681728,FMA3_64,I7
	YPRCENTRY421		yfft_r4dwpn,480K,_ac,768,573440,,,RYZEN_64
	DD			0
	PRCSTRT	10170000, 524288, 0.002819
	YPRCENTRY421		yfft_r4dwpn,512K,_ac,12,1350656,I7_64
	YPRCENTRY421		yfft_r4dwpn,512K,_ac,11,1204736
	YPRCENTRY421		yfft_r4dwpn,512K,_ac,10,653824,,I7_32+FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	11380000, 589824, 0.003433
	YPRCENTRY421		yfft_r4dwpn,576K,_ac,4608,1517568
	YPRCENTRY421		yfft_r4dwpn,576K,_ac,2304,1353728
	YPRCENTRY421		yfft_r4dwpn,576K,_ac,1536,815360,FMA3_64,I7
	YPRCENTRY421		yfft_r4dwpn,576K,_ac,768,669184,,,RYZEN_64
	DD			0
	PRCSTRT	12640000, 655360, 0.003735
	YPRCENTRY421		yfft_r4dwpn,640K,_ac,5120,1686528
	YPRCENTRY421		yfft_r4dwpn,640K,_ac,2560,1503744
	YPRCENTRY421		yfft_r4dwpn,640K,_ac,1280,809984,,I7+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,640K,_ac,10,757248,,FMA3_64
	DD			0
	PRCSTRT	15110000, 786432, 0.004558
	YPRCENTRY421		yfft_r4dwpn,768K,_ac,6144,2021376
	YPRCENTRY421		yfft_r4dwpn,768K,_ac,3072,1803776
	YPRCENTRY421		yfft_r4dwpn,768K,_ac,11,1081600,FMA3_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,768K,_ac,1536,968192
	YPRCENTRY421		yfft_r4dwpn,768K,_ac,10,885760,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,768K,_ac,768,775680
	DD			0
	PRCSTRT	15730000, 819200, 0.004927
	YPRCENTRY421		yfft_r4dwpn,800K,_ac,6400,2104832
	YPRCENTRY421		yfft_r4dwpn,800K,_ac,1280,937984,,I7+FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	16870000, 884736, 0.005499
	YPRCENTRY421		yfft_r4dwpn,864K,_ac,2304,1214208,,FMA3_64+RYZEN_64,I7
	DD			0
	PRCSTRT	18800000, 983040, 0.006015
	YPRCENTRY421		yfft_r4dwpn,960K,_ac,7680,2525184
	YPRCENTRY421		yfft_r4dwpn,960K,_ac,3840,2250752
	YPRCENTRY421		yfft_r4dwpn,960K,_ac,2560,1347840,,I7+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,960K,_ac,1536,1120768
	YPRCENTRY421		yfft_r4dwpn,960K,_ac,1280,1099264,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,960K,_ac,768,881664
	DD			0
	PRCSTRT	20060000, 1048576, 0.006127
	YPRCENTRY421		yfft_r4dwpn,1M,_ac,13,2693120
	YPRCENTRY421		yfft_r4dwpn,1M,_ac,12,2401792,I7_64
	YPRCENTRY421		yfft_r4dwpn,1M,_ac,11,1283584,,I7_32+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1M,_ac,10,1025024,,,RYZEN_64
	DD			0
	PRCSTRT	22440000, 1179648, 0.007165
	YPRCENTRY421		yfft_r4dwpn,1152K,_ac,9216,3030016
	YPRCENTRY421		yfft_r4dwpn,1152K,_ac,4608,2699776
	YPRCENTRY421		yfft_r4dwpn,1152K,_ac,3072,1615104,,I7_32+FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1152K,_ac,2304,1440768,,I7_64
	YPRCENTRY421		yfft_r4dwpn,1152K,_ac,1536,1314816
	DD			0
	PRCSTRT	24950000, 1310720, 0.007975
	YPRCENTRY421		yfft_r4dwpn,1280K,_ac,10240,3364864
	YPRCENTRY421		yfft_r4dwpn,1280K,_ac,5120,2999808
	YPRCENTRY421		yfft_r4dwpn,1280K,_ac,2560,1598976,,I7+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1280K,_ac,11,1485312
	YPRCENTRY421		yfft_r4dwpn,1280K,_ac,1280,1271296,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1280K,_ac,10,1163776
	DD			0
	PRCSTRT	27890000, 1474560, 0.009692
	YPRCENTRY421		yfft_r4dwpn,1440K,_ac,3840,2012928,,I7_64+FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1440K,_ac,2304,1667072,,I7_32
	DD			0
	PRCSTRT	29770000, 1572864, 0.009570
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,12288,4037632
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,6144,3596800
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,12,2147584,I7_64,I7_32+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,3072,1915392
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,11,1744896
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,1536,1519616
	YPRCENTRY421		yfft_r4dwpn,1536K,_ac,10,1220096,,,RYZEN_64
	DD			0
	PRCSTRT	31020000, 1638400, 0.010447
	YPRCENTRY421		yfft_r4dwpn,1600K,_ac,12800,4204544
	YPRCENTRY421		yfft_r4dwpn,1600K,_ac,6400,3745792
	YPRCENTRY421		yfft_r4dwpn,1600K,_ac,2560,1849856,,I7+FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1600K,_ac,1280,1442816,,,RYZEN_64
	DD			0
	PRCSTRT	33320000, 1769472, 0.011461
	YPRCENTRY421		yfft_r4dwpn,1728K,_ac,4608,2412800,,I7_32+FMA3_64+RYZEN_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,1728K,_ac,2304,1959424
	DD			0
	PRCSTRT	37040000, 1966080, 0.012478
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,15360,5045248
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,7680,4493824
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,5120,2680064,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,3840,2386944,,I7_64
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,3072,2215424,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,2560,2174976
	YPRCENTRY421		yfft_r4dwpn,1920K,_ac,1536,1723904
	DD			0
	PRCSTRT	39520000, 2097152, 0.012734
	YPRCENTRY421		yfft_r4dwpn,2M,_ac,14,5381120
	YPRCENTRY421		yfft_r4dwpn,2M,_ac,13,4792832
	YPRCENTRY421		yfft_r4dwpn,2M,_ac,12,2546176,I7_64,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,2M,_ac,11,2015232,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2M,_ac,10,1800192
	DD			0
	PRCSTRT	44230000, 2359296, 0.015112
	YPRCENTRY421		yfft_r4dwpn,2304K,_ac,9216,5391872
	YPRCENTRY421		yfft_r4dwpn,2304K,_ac,6144,3211520,,I7_32+FMA3_64+RYZEN_64,I7_64
	YPRCENTRY421		yfft_r4dwpn,2304K,_ac,4608,2860544
	YPRCENTRY421		yfft_r4dwpn,2304K,_ac,3072,2606080
	YPRCENTRY421		yfft_r4dwpn,2304K,_ac,2304,2262528
	DD			0
	PRCSTRT	46010000, 2457600, 0.016488
	YPRCENTRY421		yfft_r4dwpn,2400K,_ac,6400,3344128,,I7+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,2400K,_ac,3840,2760704,,,FMA3_64
	DD			0
	PRCSTRT	49160000, 2621440, 0.016596
	YPRCENTRY421		yfft_r4dwpn,2560K,_ac,20480,6724608
	YPRCENTRY421		yfft_r4dwpn,2560K,_ac,10240,5988864
	YPRCENTRY421		yfft_r4dwpn,2560K,_ac,5120,3176960,,I7_32+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,2560K,_ac,12,2944512,I7_64,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2560K,_ac,2560,2510848
	YPRCENTRY421		yfft_r4dwpn,2560K,_ac,11,2285056
	DD			0
	PRCSTRT	55020000, 2949120, 0.020031
	YPRCENTRY421		yfft_r4dwpn,2880K,_ac,7680,4010240,,I7_64
	YPRCENTRY421		yfft_r4dwpn,2880K,_ac,4608,3308032,,I7_32,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,2880K,_ac,3840,3249664,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,2880K,_ac,2304,2565120
	DD			0
	PRCSTRT	58730000, 3145728, 0.020198
	YPRCENTRY421		yfft_r4dwpn,3M,_ac,12288,7185920
	YPRCENTRY421		yfft_r4dwpn,3M,_ac,13,4276480
	YPRCENTRY421		yfft_r4dwpn,3M,_ac,6144,3806720,,I7
	YPRCENTRY421		yfft_r4dwpn,3M,_ac,12,3466240,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3M,_ac,3072,3007488
	DD			0
	PRCSTRT	61110000, 3276800, 0.021621
	YPRCENTRY421		yfft_r4dwpn,3200K,_ac,25600,8403968
	YPRCENTRY421		yfft_r4dwpn,3200K,_ac,12800,7483904
	YPRCENTRY421		yfft_r4dwpn,3200K,_ac,6400,3963904,,I7,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,3200K,_ac,5120,3673600
	YPRCENTRY421		yfft_r4dwpn,3200K,_ac,2560,2846208,,,FMA3_64
	DD			0
	PRCSTRT	65640000, 3538944, 0.023780
	YPRCENTRY421		yfft_r4dwpn,3456K,_ac,9216,4809984,,I7
	YPRCENTRY421		yfft_r4dwpn,3456K,_ac,4608,3895296,,,FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,3456K,_ac,2304,2662400
	DD			0
	PRCSTRT	73070000, 3932160, 0.026437
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,15360,8979968
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,10240,5341440
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,7680,4752896,,I7
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,6144,4401664,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,5120,4326400,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,3840,3749376
	YPRCENTRY421		yfft_r4dwpn,3840K,_ac,3072,3408384
	DD			0
	PRCSTRT	76060000, 4096000, 0.029047
	YPRCENTRY421		yfft_r4dwpn,4000K,_ac,6400,4583424,,,RYZEN_64
	DD			0
	PRCSTRT	77950000, 4194304, 0.027172
	YPRCENTRY421		yfft_r4dwpn,4M,_ac,14,9577984
	YPRCENTRY421		yfft_r4dwpn,4M,_ac,13,5068288,I7_64,I7_32
	YPRCENTRY421		yfft_r4dwpn,4M,_ac,12,3998720,,RYZEN_64,FMA3_64
	DD			0
	PRCSTRT	87120000, 4718592, 0.031494
	YPRCENTRY421		yfft_r4dwpn,4608K,_ac,12288,6407424
	YPRCENTRY421		yfft_r4dwpn,4608K,_ac,9216,5700096,,I7
	YPRCENTRY421		yfft_r4dwpn,4608K,_ac,6144,5185536
	YPRCENTRY421		yfft_r4dwpn,4608K,_ac,4608,4493312,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,4608K,_ac,3072,3530240,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,4608K,_ac,2304,3979776
	DD			0
	PRCSTRT	90930000, 4915200, 0.034693
	YPRCENTRY421		yfft_r4dwpn,4800K,_ac,12800,6672640
	YPRCENTRY421		yfft_r4dwpn,4800K,_ac,7680,5495296,,I7_64,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,4800K,_ac,6400,5400064,,I7_32,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,4800K,_ac,3840,4248576
	DD			0
	PRCSTRT	96940000, 5242880, 0.035846
	YPRCENTRY421		yfft_r4dwpn,5M,_ac,20480,11970048
	YPRCENTRY421		yfft_r4dwpn,5M,_ac,10240,6329856,,I7_64
	YPRCENTRY421		yfft_r4dwpn,5M,_ac,13,5859840,,I7_32
	YPRCENTRY421		yfft_r4dwpn,5M,_ac,5120,4989952,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,5M,_ac,12,4530688,,,FMA3_64
	DD			0
	PRCSTRT	108500000, 5898240, 0.041629
	YPRCENTRY421		yfft_r4dwpn,5760K,_ac,15360,8004864
	YPRCENTRY421		yfft_r4dwpn,5760K,_ac,9216,6589952,I7
	YPRCENTRY421		yfft_r4dwpn,5760K,_ac,7680,6475776
	YPRCENTRY421		yfft_r4dwpn,5760K,_ac,4608,5090816
	YPRCENTRY421		yfft_r4dwpn,5760K,_ac,3840,4395008,,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	115600000, 6291456, 0.043121
	YPRCENTRY421		yfft_r4dwpn,6M,_ac,14,8537344,I7_64
	YPRCENTRY421		yfft_r4dwpn,6M,_ac,12288,7592448
	YPRCENTRY421		yfft_r4dwpn,6M,_ac,13,6905856,,I7_32
	YPRCENTRY421		yfft_r4dwpn,6M,_ac,6144,5980160
	YPRCENTRY421		yfft_r4dwpn,6M,_ac,12,4685312,,FMA3_64,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,6M,_ac,3072,5289984
	DD			0
	PRCSTRT	120500000, 6553600, 0.045826
	YPRCENTRY421		yfft_r4dwpn,6400K,_ac,25600,14960128
	YPRCENTRY421		yfft_r4dwpn,6400K,_ac,12800,7906816,,I7
	YPRCENTRY421		yfft_r4dwpn,6400K,_ac,10240,7318016,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,6400K,_ac,6400,6227456,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,6400K,_ac,5120,5652992
	DD			0
	PRCSTRT	129300000, 7077888, 0.051987
	YPRCENTRY421		yfft_r4dwpn,6912K,_ac,9216,7767040,,I7
	YPRCENTRY421		yfft_r4dwpn,6912K,_ac,4608,5261824,,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	143600000, 7864320, 0.054932
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,20480,10667264
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,15360,9484800,,I7
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,12288,8777216
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,10240,8626176
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,7680,7467008,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,6144,6774272,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,5120,5840384
	YPRCENTRY421		yfft_r4dwpn,7680K,_ac,3840,6597120
	DD			0
	PRCSTRT	149700000, 8192000, 0.062103
	YPRCENTRY421		yfft_r4dwpn,8000K,_ac,12800,9140736,I7_64
	YPRCENTRY421		yfft_r4dwpn,8000K,_ac,6400,7054336
	DD			0
	PRCSTRT	153400000, 8388608, 0.059916
	YPRCENTRY421		yfft_r4dwpn,8M,_ac,14,10115584,I7_64
	YPRCENTRY421		yfft_r4dwpn,8M,_ac,13,7962624,,I7_32
	YPRCENTRY421		yfft_r4dwpn,8M,_ac,12,7034880,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	171700000, 9437184, 0.069579
	YPRCENTRY421		yfft_r4dwpn,9M,_ac,12288,10347520,,I7_32,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,9M,_ac,9216,8954880
	YPRCENTRY421		yfft_r4dwpn,9M,_ac,6144,6994432,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,9M,_ac,4608,7906304
	DD			0
	PRCSTRT	178800000, 9830400, 0.072967
	YPRCENTRY421		yfft_r4dwpn,9600K,_ac,25600,13329664
	YPRCENTRY421		yfft_r4dwpn,9600K,_ac,15360,10964480,I7_64
	YPRCENTRY421		yfft_r4dwpn,9600K,_ac,12800,10776576,,I7_32
	YPRCENTRY421		yfft_r4dwpn,9600K,_ac,7680,8457728,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,9600K,_ac,6400,7282688,,,RYZEN_64
	DD			0
	PRCSTRT	190600000, 10485760, 0.078834
	YPRCENTRY421		yfft_r4dwpn,10M,_ac,20480,12638720,,I7_64
	YPRCENTRY421		yfft_r4dwpn,10M,_ac,14,11693568
	YPRCENTRY421		yfft_r4dwpn,10M,_ac,10240,9945088,,,FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,10M,_ac,13,9018880
	YPRCENTRY421		yfft_r4dwpn,10M,_ac,5120,8779776,I7_32
	DD			0
	PRCSTRT	213500000, 11796480, 0.088274
	YPRCENTRY421		yfft_r4dwpn,11520K,_ac,15360,12928000,,I7
	YPRCENTRY421		yfft_r4dwpn,11520K,_ac,9216,10142208
	YPRCENTRY421		yfft_r4dwpn,11520K,_ac,7680,8727040,,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	227900000, 12582912, 0.096630
	YPRCENTRY421		yfft_r4dwpn,12M,_ac,14,13788160
	YPRCENTRY421		yfft_r4dwpn,12M,_ac,12288,11928576,,,FMA3_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,12M,_ac,13,9304576,I7_32
	YPRCENTRY421		yfft_r4dwpn,12M,_ac,6144,10523648
	DD			0
	PRCSTRT	237100000, 13107200, 0.097776
	YPRCENTRY421		yfft_r4dwpn,12800K,_ac,25600,15792640,,I7
	YPRCENTRY421		yfft_r4dwpn,12800K,_ac,20480,14609920
	YPRCENTRY421		yfft_r4dwpn,12800K,_ac,12800,12423168,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,12800K,_ac,10240,11263488,,,FMA3_64
	YPRCENTRY421		yfft_r4dwpn,12800K,_ac,6400,10959360
	DD			0
	PRCSTRT	254500000, 14155776, 0.115310
	YPRCENTRY421		yfft_r4dwpn,13824K,_ac,9216,10460672,I7,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	283500000, 15728640, 0.128292
	YPRCENTRY421		yfft_r4dwpn,15M,_ac,20480,17228800
	YPRCENTRY421		yfft_r4dwpn,15M,_ac,12288,13509120
	YPRCENTRY421		yfft_r4dwpn,15M,_ac,10240,11614720,I7_32,FMA3_64,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,15M,_ac,7680,13140992
	DD			0
	PRCSTRT	295100000, 16384000, 0.130344
	YPRCENTRY421		yfft_r4dwpn,16000K,_ac,25600,18255360,I7
	YPRCENTRY421		yfft_r4dwpn,16000K,_ac,12800,14069248,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	302400000, 16777216, 0.133010
	YPRCENTRY421		yfft_r4dwpn,16M,_ac,14,15893504,I7_64+RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,16M,_ac,13,14013440,I7_32,FMA3_64
	DD			0
	PRCSTRT	338200000, 18874368, 0.153546
	YPRCENTRY421		yfft_r4dwpn,18M,_ac,12288,13925888,I7_64,FMA3_64,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,18M,_ac,9216,15759360,I7_32
	DD			0
	PRCSTRT	352300000, 19660800, 0.156971
	YPRCENTRY421		yfft_r4dwpn,19200K,_ac,25600,21529600,,I7
	YPRCENTRY421		yfft_r4dwpn,19200K,_ac,15360,16876032
	YPRCENTRY421		yfft_r4dwpn,19200K,_ac,12800,14502400,,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	375800000, 20971520, 0.172761
	YPRCENTRY421		yfft_r4dwpn,20M,_ac,20480,19858432
	YPRCENTRY421		yfft_r4dwpn,20M,_ac,14,17998336
	YPRCENTRY421		yfft_r4dwpn,20M,_ac,10240,17503232,I7,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	419600000, 23592960, 0.198884
	YPRCENTRY421		yfft_r4dwpn,23040K,_ac,15360,17391104,I7,FMA3_64,RYZEN_64
	DD			0
	PRCSTRT	448500000, 25165824, 0.209253
	YPRCENTRY421		yfft_r4dwpn,24M,_ac,12288,20994048,I7_32,FMA3_64,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,24M,_ac,14,18546176,I7_64
	DD			0
	PRCSTRT	467100000, 26214400, 0.212801
	YPRCENTRY421		yfft_r4dwpn,25M,_ac,25600,24814592,I7
	YPRCENTRY421		yfft_r4dwpn,25M,_ac,20480,22487552
	YPRCENTRY421		yfft_r4dwpn,25M,_ac,12800,21865472,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	557000000, 31457280, 0.265964
	YPRCENTRY421		yfft_r4dwpn,30M,_ac,20480,23166464,,,RYZEN_64
	YPRCENTRY421		yfft_r4dwpn,30M,_ac,15360,26228736,I7,FMA3_64
	DD			0
	PRCSTRT	580500000, 32768000, 0.282467
	YPRCENTRY421		yfft_r4dwpn,32000K,_ac,25600,28099072,I7_64,,RYZEN_64
	DD			0
	PRCSTRT	594700000, 33554432, 0.287437
	YPRCENTRY421		yfft_r4dwpn,32M,_ac,14,27973632,I7,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	693000000, 39321600, 0.330
	YPRCENTRY421		yfft_r4dwpn,38400K,_ac,25600,28941824,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	740000000, 41943040, 0.360
	YPRCENTRY421		yfft_r4dwpn,40M,_ac,20480,34953216,,,FMA3_64+RYZEN_64
	DD			0
	PRCSTRT	920000000, 52428800, 0.450
	YPRCENTRY421		yfft_r4dwpn,50M,_ac,25600,43677696,,,FMA3_64+RYZEN_64
	DD			0
	DD	0

;; Jump tables for the AVX-512 optimized code

IFDEF X86_64
zjmptable DD	0
	org	$-4
	PRCSTRT	22469,	1024,	0.000004
	ZPRCENTRY1		zfft_r4dwpn,1K,,13312,SKX
	DD			0
	PRCSTRT	33423,	1536,	0.000006
	ZPRCENTRY1		zfft_r4dwpn,1536,,19008,SKX
	DD			0
	PRCSTRT	44307,	2048,	0.000008
	ZPRCENTRY1		zfft_r4dwpn,2K,,24704,SKX
	DD			0
	PRCSTRT	98407,	4608,	0.000020
	ZPRCENTRY1		zfft_r4dwpn,4608,,59328,SKX
	DD			0
	PRCSTRT	108700,	5120,	0.000022
	ZPRCENTRY1		zfft_r4dwpn,5K,,58880,SKX
	DD			0
	PRCSTRT	130200,	6144,	0.000025
	ZPRCENTRY1		zfft_r4dwpn,6K,,78464,SKX
	DD			0
	PRCSTRT	151300,	7168,	0.000029
	ZPRCENTRY1		zfft_r4dwpn,7K,,81664,SKX
	DD			0
	PRCSTRT	162400,	7680,	0.000031
	ZPRCENTRY1		zfft_r4dwpn,7680,,87360,SKX
	DD			0
	PRCSTRT	172600,	8192,	0.000033
	ZPRCENTRY1		zfft_r4dwpn,8K,,101248
	ZPRCENTRY421		zfft_r4dwpn,8K,,64,24192,,,SKX
	DD			0
	PRCSTRT	194300,	9216,	0.000038
	ZPRCENTRY1		zfft_r4dwpn,9K,,110592,SKX
	DD			0
	PRCSTRT	214800,	10240,	0.000043
	ZPRCENTRY1		zfft_r4dwpn,10K,,115840,SKX
	DD			0
	PRCSTRT	225600,	10752,	0.000046
	ZPRCENTRY1		zfft_r4dwpn,10752,,121536,SKX
	DD			0
	PRCSTRT	256500,	12288,	0.000055
	ZPRCENTRY1		zfft_r4dwpn,12K,,146816
	ZPRCENTRY421		zfft_r4dwpn,12K,,64,30336,,,SKX
	DD			0
	PRCSTRT	267400,	12800,	0.000057
	ZPRCENTRY1		zfft_r4dwpn,12800,,154560,SKX
	DD			0
	PRCSTRT	298600,	14336,	0.000065
	ZPRCENTRY1		zfft_r4dwpn,14K,,161408				;,SKX
	DD			0
	PRCSTRT	320200,	15360,	0.000073
	ZPRCENTRY1		zfft_r4dwpn,15K,,185088				;,SKX
	DD			0
	PRCSTRT	340000,	16384,	0.000072
	ZPRCENTRY1		zfft_r4dwpn,16K,,192384,SKX
	DD			0
	PRCSTRT	382900,	18432,	0.000086
	ZPRCENTRY1		zfft_r4dwpn,18K,,219264,SKX
	DD			0
	PRCSTRT	424300,	20480,	0.000089
	ZPRCENTRY1		zfft_r4dwpn,20K,,246144,SKX
	DD			0
	PRCSTRT	507400,	24576,	0.000105
	ZPRCENTRY1		zfft_r4dwpn,24K,,291712,SKX
	DD			0
	PRCSTRT	527800,	25600,	0.000108
	ZPRCENTRY1		zfft_r4dwpn,25K,,286720,SKX
	DD			0
	PRCSTRT	632000,	30720,	0.000121
	ZPRCENTRY1		zfft_r4dwpn,30K,,351872				;,SKX
	DD			0
	PRCSTRT	672200,	32768,	0.000133
	ZPRCENTRY1		zfft_r4dwpn,32K,,382848,SKX
	DD			0
	PRCSTRT	734100,	35840,	0.000140
	ZPRCENTRY1		zfft_r4dwpn,35K,,400640				;,SKX
	DD			0
	PRCSTRT	756900,	36864,	0.000145
	ZPRCENTRY1		zfft_r4dwpn,36K,,469376				;,SKX
	DD			0
	PRCSTRT	835300,	40960,	0.000183
	ZPRCENTRY1		zfft_r4dwpn,40K,,465792
	ZPRCENTRY421		zfft_r4dwpn,40K,,320,95872,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,40K,,64,73856
	DD			0
	PRCSTRT	878800, 43008,	0.000191
	ZPRCENTRY1		zfft_r4dwpn,42K,,488576				;,SKX
	DD			0
	PRCSTRT	1001000, 49152,	0.000232
	ZPRCENTRY1		zfft_r4dwpn,48K,,622464
	ZPRCENTRY421		zfft_r4dwpn,48K,,384,114816,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,48K,,64,94208
	DD			0
	PRCSTRT	1161000, 57344,	0.000270
	ZPRCENTRY421		zfft_r4dwpn,56K,,448,131712,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,56K,,64,98688
	DD			0
	PRCSTRT	1244000, 61440,	0.000301
	ZPRCENTRY421		zfft_r4dwpn,60K,,320,118400,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,60K,,64,88704
	DD			0
	PRCSTRT	1324000, 65536,	0.000300
	ZPRCENTRY421		zfft_r4dwpn,64K,,512,150656,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,64K,,64,119040
	DD			0
	PRCSTRT	1490000, 73728,	0.000367
	ZPRCENTRY421		zfft_r4dwpn,72K,,384,141440,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,72K,,64,111616
	DD			0
	PRCSTRT	1649000, 81920, 0.000386
	ZPRCENTRY421		zfft_r4dwpn,80K,,640,186496,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,80K,,64,143872
	DD			0
	PRCSTRT	1730000, 86016, 0.000443
	ZPRCENTRY421		zfft_r4dwpn,84K,,448,162432,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,84K,,64,118144
	DD			0
	PRCSTRT	1976000, 98304, 0.000458
	ZPRCENTRY421		zfft_r4dwpn,96K,,768,222336,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,96K,,512,185472
	ZPRCENTRY421		zfft_r4dwpn,96K,,64,168704
	DD			0
	PRCSTRT	2457000, 122880, 0.000639
	ZPRCENTRY421		zfft_r4dwpn,120K,,640,229504,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,120K,,64,169984
	DD			0
	PRCSTRT	2616000, 131072, 0.000594
	ZPRCENTRY421		zfft_r4dwpn,128K,,1024,294016,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,128K,,64,218368
	DD			0
	PRCSTRT	2941000, 147456, 0.000801
	ZPRCENTRY421		zfft_r4dwpn,144K,,768,273536,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,144K,,64,198912
	DD			0
	PRCSTRT	3893000, 196608, 0.000984
	ZPRCENTRY421		zfft_r4dwpn,192K,,1024,361600,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,192K,,64,257280
	DD			0
	PRCSTRT	4043000, 204800, 0.001138
	ZPRCENTRY421		zfft_r4dwpn,200K,,1600,454272
	ZPRCENTRY421		zfft_r4dwpn,200K,,320,280704,,,SKX
	DD			0
	PRCSTRT	4843000, 245760, 0.001300
	ZPRCENTRY421		zfft_r4dwpn,240K,,1920,544896
	ZPRCENTRY421		zfft_r4dwpn,240K,,384,333440,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,240K,,320,366592
	DD			0
	PRCSTRT	5623000, 286720, 0.001668
	ZPRCENTRY421		zfft_r4dwpn,280K,,2240,633472,,SKX
	ZPRCENTRY421		zfft_r4dwpn,280K,,448,384128
	ZPRCENTRY421		zfft_r4dwpn,280K,,320,373120
	DD			0
	PRCSTRT	5798000, 294912, 0.001575
	ZPRCENTRY421		zfft_r4dwpn,288K,,2304,658560,,SKX
	ZPRCENTRY421		zfft_r4dwpn,288K,,384,435712
	DD			0
	PRCSTRT	6022000, 307200, 0.001679
	ZPRCENTRY421		zfft_r4dwpn,300K,,1600,558720
	ZPRCENTRY421		zfft_r4dwpn,300K,,320,313984,,,SKX
	DD			0
	PRCSTRT	6408000, 327680, 0.001679
	ZPRCENTRY421		zfft_r4dwpn,320K,,2560,724096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,320K,,512,436864
	ZPRCENTRY421		zfft_r4dwpn,320K,,320,459008
	DD			0
	PRCSTRT	6730000, 344064, 0.001969
	ZPRCENTRY421		zfft_r4dwpn,336K,,2688,759936
	ZPRCENTRY421		zfft_r4dwpn,336K,,448,502784
	ZPRCENTRY421		zfft_r4dwpn,336K,,384,442752,,,SKX
	DD			0
	PRCSTRT	7219000, 368640, 0.001969
	ZPRCENTRY421		zfft_r4dwpn,360K,,1920,669824
	ZPRCENTRY421		zfft_r4dwpn,360K,,384,371328,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,360K,,320,404480
	DD			0
	PRCSTRT	7682000, 393216, 0.002047
	ZPRCENTRY421		zfft_r4dwpn,384K,,3072,875648
	ZPRCENTRY421		zfft_r4dwpn,384K,,512,571904
	ZPRCENTRY421		zfft_r4dwpn,384K,,384,545024,,,SKX
	DD			0
	PRCSTRT	7815000, 401408, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,392K,,3136,884352
	ZPRCENTRY421		zfft_r4dwpn,392K,,448,510336,,,SKX
	DD			0
	PRCSTRT	7985000, 409600, 0.002250
	ZPRCENTRY421		zfft_r4dwpn,400K,,3200,903296
	ZPRCENTRY421		zfft_r4dwpn,400K,,640,540288,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,400K,,320,551424
	DD			0
	PRCSTRT	8384000, 430080, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,420K,,2240,778880
	ZPRCENTRY421		zfft_r4dwpn,420K,,448,426624,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,420K,,320,413056
	DD			0
	PRCSTRT	8640000, 442368, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,432K,,2304,808064
	ZPRCENTRY421		zfft_r4dwpn,432K,,384,478720,,,SKX
	DD			0
	PRCSTRT	8907000, 458752, 0.002446
	ZPRCENTRY421		zfft_r4dwpn,448K,,3584,1010816
	ZPRCENTRY421		zfft_r4dwpn,448K,,512,579968
	ZPRCENTRY421		zfft_r4dwpn,448K,,448,628992,,,SKX
	DD			0
	PRCSTRT	9547000, 491520, 0.002709
	ZPRCENTRY421		zfft_r4dwpn,480K,,3840,1082496
	ZPRCENTRY421		zfft_r4dwpn,480K,,2560,889984
	ZPRCENTRY421		zfft_r4dwpn,480K,,768,643712
	ZPRCENTRY421		zfft_r4dwpn,480K,,640,708096
	ZPRCENTRY421		zfft_r4dwpn,480K,,512,483968,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,480K,,384,654336
	ZPRCENTRY421		zfft_r4dwpn,480K,,320,643840
	DD			0
	PRCSTRT	10020000, 516096, 0.002712
	ZPRCENTRY421		zfft_r4dwpn,504K,,2688,934016
	ZPRCENTRY421		zfft_r4dwpn,504K,,448,550912
	ZPRCENTRY421		zfft_r4dwpn,504K,,384,487808,,,SKX
	DD			0
	PRCSTRT	10150000, 524288, 0.002712
	ZPRCENTRY421		zfft_r4dwpn,512K,,4096,1162368
	ZPRCENTRY421		zfft_r4dwpn,512K,,512,715008,,,SKX
	DD			0
	PRCSTRT	11080000, 573440, 0.003284
	ZPRCENTRY421		zfft_r4dwpn,560K,,4480,1261696
	ZPRCENTRY421		zfft_r4dwpn,560K,,640,717184,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,560K,,448,755200
	DD			0
	PRCSTRT	11430000, 589824, 0.003313
	ZPRCENTRY421		zfft_r4dwpn,576K,,4608,1303680
	ZPRCENTRY421		zfft_r4dwpn,576K,,3072,1074304
	ZPRCENTRY421		zfft_r4dwpn,576K,,768,844288
	ZPRCENTRY421		zfft_r4dwpn,576K,,512,625152,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,576K,,384,763648
	DD			0
	PRCSTRT	11630000, 602112, 0.003313
	ZPRCENTRY421		zfft_r4dwpn,588K,,3136,1087104
	ZPRCENTRY421		zfft_r4dwpn,588K,,448,560512,,,SKX
	DD			0
	PRCSTRT	11880000, 614400, 0.003381
	ZPRCENTRY421		zfft_r4dwpn,600K,,3200,1110144
	ZPRCENTRY421		zfft_r4dwpn,600K,,640,596608,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,600K,,320,600064
	DD			0
	PRCSTRT	12640000, 655360, 0.003597
	ZPRCENTRY421		zfft_r4dwpn,640K,,5120,1440896
	ZPRCENTRY421		zfft_r4dwpn,640K,,1024,850560
	ZPRCENTRY421		zfft_r4dwpn,640K,,640,884992,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,640K,,512,858112
	ZPRCENTRY421		zfft_r4dwpn,640K,,320,828672
	DD			0
	PRCSTRT	13280000, 688128, 0.004010
	ZPRCENTRY421		zfft_r4dwpn,672K,,5376,1512576
	ZPRCENTRY421		zfft_r4dwpn,672K,,3584,1242240
	ZPRCENTRY421		zfft_r4dwpn,672K,,768,854400
	ZPRCENTRY421		zfft_r4dwpn,672K,,512,635264,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,672K,,448,881408
	DD			0
	PRCSTRT	14210000, 737280, 0.004431
	ZPRCENTRY421		zfft_r4dwpn,720K,,3840,1330304
	ZPRCENTRY421		zfft_r4dwpn,720K,,768,709248
	ZPRCENTRY421		zfft_r4dwpn,720K,,640,771584,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,720K,,384,708608
	ZPRCENTRY421		zfft_r4dwpn,720K,,320,696576
	DD			0
	PRCSTRT	15140000, 786432, 0.004390
	ZPRCENTRY421		zfft_r4dwpn,768K,,6144,1735808,,SKX
	ZPRCENTRY421		zfft_r4dwpn,768K,,4096,1426560
	ZPRCENTRY421		zfft_r4dwpn,768K,,1024,1116672
	ZPRCENTRY421		zfft_r4dwpn,768K,,768,1054976
	ZPRCENTRY421		zfft_r4dwpn,768K,,512,1001216
	ZPRCENTRY421		zfft_r4dwpn,768K,,384,982272
	DD			0
	PRCSTRT	15730000, 819200, 0.004791
	ZPRCENTRY421		zfft_r4dwpn,800K,,6400,1809536,,SKX
	ZPRCENTRY421		zfft_r4dwpn,800K,,640,1061888
	DD			0
	PRCSTRT	16500000, 860160, 0.005143
	ZPRCENTRY421		zfft_r4dwpn,840K,,4480,1550464
	ZPRCENTRY421		zfft_r4dwpn,840K,,640,782720,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,840K,,448,815104
	DD			0
	PRCSTRT	17000000, 884736, 0.005290
	ZPRCENTRY421		zfft_r4dwpn,864K,,4608,1600640,SKX
	ZPRCENTRY421		zfft_r4dwpn,864K,,768,918016
	ZPRCENTRY421		zfft_r4dwpn,864K,,384,822016
	DD			0
	PRCSTRT	17570000, 917504, 0.005404
	ZPRCENTRY421		zfft_r4dwpn,896K,,7168,2014336
	ZPRCENTRY421		zfft_r4dwpn,896K,,1024,1128832
	ZPRCENTRY421		zfft_r4dwpn,896K,,448,1133824,,,SKX
	DD			0
	PRCSTRT	18820000, 983040, 0.005806
	ZPRCENTRY421		zfft_r4dwpn,960K,,7680,2169984
	ZPRCENTRY421		zfft_r4dwpn,960K,,5120,1770624
	ZPRCENTRY421		zfft_r4dwpn,960K,,1024,934528
	ZPRCENTRY421		zfft_r4dwpn,960K,,768,1265664
	ZPRCENTRY421		zfft_r4dwpn,960K,,640,1238784,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,960K,,512,923648
	ZPRCENTRY421		zfft_r4dwpn,960K,,320,892160
	DD			0
	PRCSTRT	19540000, 1024000, 0.006539
	ZPRCENTRY421		zfft_r4dwpn,1000K,,1600,1314944,,SKX
	DD			0
	PRCSTRT	19760000, 1032192, 0.006539
	ZPRCENTRY421		zfft_r4dwpn,1008K,,5376,1858688,SKX
	ZPRCENTRY421		zfft_r4dwpn,1008K,,768,930176
	ZPRCENTRY421		zfft_r4dwpn,1008K,,448,945408
	DD			0
	PRCSTRT	20030000, 1048576, 0.005906
	ZPRCENTRY421		zfft_r4dwpn,1M,,8192,2309248,SKX
	ZPRCENTRY421		zfft_r4dwpn,1M,,1024,1394944
	ZPRCENTRY421		zfft_r4dwpn,1M,,512,1287424
	DD			0
	PRCSTRT	22540000, 1179648, 0.007022
	ZPRCENTRY421		zfft_r4dwpn,1152K,,9216,2600064
	ZPRCENTRY421		zfft_r4dwpn,1152K,,6144,2131072,SKX
	ZPRCENTRY421		zfft_r4dwpn,1152K,,1024,1210880
	ZPRCENTRY421		zfft_r4dwpn,1152K,,768,1476352
	ZPRCENTRY421		zfft_r4dwpn,1152K,,512,1070848
	ZPRCENTRY421		zfft_r4dwpn,1152K,,384,1051904
	DD			0
	PRCSTRT	23440000, 1228800, 0.007757
	ZPRCENTRY421		zfft_r4dwpn,1200K,,6400,2221184
	ZPRCENTRY421		zfft_r4dwpn,1200K,,1920,1574528,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1200K,,1600,1728512
	ZPRCENTRY421		zfft_r4dwpn,1200K,,640,1138688
	DD			0
	PRCSTRT	24940000, 1310720, 0.007707
	ZPRCENTRY421		zfft_r4dwpn,1280K,,10240,2890880,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1280K,,1024,1673216
	ZPRCENTRY421		zfft_r4dwpn,1280K,,640,1592576
	DD			0
	PRCSTRT	26170000, 1376256, 0.008592
	ZPRCENTRY421		zfft_r4dwpn,1344K,,7168,2475136,SKX
	ZPRCENTRY421		zfft_r4dwpn,1344K,,1024,1225088
	ZPRCENTRY421		zfft_r4dwpn,1344K,,448,1209600
	DD			0
	PRCSTRT	27190000, 1433600, 0.009242
	ZPRCENTRY421		zfft_r4dwpn,1400K,,2240,1832064
	ZPRCENTRY421		zfft_r4dwpn,1400K,,1600,1745280,,SKX
	DD			0
	PRCSTRT	28050000, 1474560, 0.009249
	ZPRCENTRY421		zfft_r4dwpn,1440K,,7680,2663552,SKX
	ZPRCENTRY421		zfft_r4dwpn,1440K,,2304,1890944
	ZPRCENTRY421		zfft_r4dwpn,1440K,,1920,2070016
	ZPRCENTRY421		zfft_r4dwpn,1440K,,768,1353728
	ZPRCENTRY421		zfft_r4dwpn,1440K,,640,1319680
	DD			0
	PRCSTRT	29140000, 1536000, 0.009252
	ZPRCENTRY421		zfft_r4dwpn,1500K,,1600,1440384,,SKX
	DD			0
	PRCSTRT	29880000, 1572864, 0.009264
	ZPRCENTRY421		zfft_r4dwpn,1536K,,12288,3464320,SKX
	ZPRCENTRY421		zfft_r4dwpn,1536K,,8192,2835584
	ZPRCENTRY421		zfft_r4dwpn,1536K,,1024,1951488
	ZPRCENTRY421		zfft_r4dwpn,1536K,,768,1897728
	ZPRCENTRY421		zfft_r4dwpn,1536K,,512,1369344
	DD			0
	PRCSTRT	30990000, 1638400, 0.010048
	ZPRCENTRY421		zfft_r4dwpn,1600K,,12800,3591296
	ZPRCENTRY421		zfft_r4dwpn,1600K,,2560,2091648
	ZPRCENTRY421		zfft_r4dwpn,1600K,,1600,2158848,,,SKX
	DD			0
	PRCSTRT	32550000, 1720320, 0.011459
	ZPRCENTRY421		zfft_r4dwpn,1680K,,2688,2195072
	ZPRCENTRY421		zfft_r4dwpn,1680K,,2240,2409472
	ZPRCENTRY421		zfft_r4dwpn,1680K,,1920,2089344,,SKX
	DD			0
	PRCSTRT	33540000, 1769472, 0.011102
	ZPRCENTRY421		zfft_r4dwpn,1728K,,9216,3191936
	ZPRCENTRY421		zfft_r4dwpn,1728K,,2304,2484736,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1728K,,768,1568512
	DD			0
	PRCSTRT	34870000, 1843200, 0.0115
	ZPRCENTRY421		zfft_r4dwpn,1800K,,1920,1723008,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1800K,,1600,1868800
	DD			0
	PRCSTRT	37120000, 1966080, 0.012046
	ZPRCENTRY421		zfft_r4dwpn,1920K,,15360,4316288,SKX
	ZPRCENTRY421		zfft_r4dwpn,1920K,,10240,3548288
	ZPRCENTRY421		zfft_r4dwpn,1920K,,3072,2513536
	ZPRCENTRY421		zfft_r4dwpn,1920K,,2560,2750976
	ZPRCENTRY421		zfft_r4dwpn,1920K,,1920,2584832
	ZPRCENTRY421		zfft_r4dwpn,1920K,,1024,1783808
	ZPRCENTRY421		zfft_r4dwpn,1920K,,640,1686784
	DD			0
	PRCSTRT	37760000, 2007040, 0.012046
	ZPRCENTRY421		zfft_r4dwpn,1960K,,3136,2556032
	ZPRCENTRY421		zfft_r4dwpn,1960K,,2240,2431360,,SKX
	DD			0
	PRCSTRT	38580000, 2048000, 0.013344
	ZPRCENTRY421		zfft_r4dwpn,2000K,,3200,2608768
	ZPRCENTRY421		zfft_r4dwpn,2000K,,1600,2589184			;,,,SKX
	DD			0
	PRCSTRT	38930000, 2064384, 0.013399
	ZPRCENTRY421		zfft_r4dwpn,2016K,,2688,2887168
	ZPRCENTRY421		zfft_r4dwpn,2016K,,2304,2507136			;,,,SKX
	DD			0
	PRCSTRT	39490000, 2097152, 0.012290
	ZPRCENTRY421		zfft_r4dwpn,2M,,16384,4611200,SKX
	ZPRCENTRY421		zfft_r4dwpn,2M,,1024,2508032
	DD			0
	PRCSTRT	40440000, 2150400, 0.0135
	ZPRCENTRY421		zfft_r4dwpn,2100K,,2240,2003584,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2100K,,1600,1887616
	DD			0
	PRCSTRT	41740000, 2211840, 0.014
	ZPRCENTRY421		zfft_r4dwpn,2160K,,2304,2067072
	ZPRCENTRY421		zfft_r4dwpn,2160K,,1920,2235904,,SKX
	DD			0
	PRCSTRT	43050000, 2293760, 0.014710
	ZPRCENTRY421		zfft_r4dwpn,2240K,,17920,5024896,SKX
	ZPRCENTRY421		zfft_r4dwpn,2240K,,3584,2919040
	ZPRCENTRY421		zfft_r4dwpn,2240K,,2560,2775424
	ZPRCENTRY421		zfft_r4dwpn,2240K,,2240,3008768
	DD			0
	PRCSTRT	44420000, 2359296, 0.014635
	ZPRCENTRY421		zfft_r4dwpn,2304K,,18432,5225600,SKX
	ZPRCENTRY421		zfft_r4dwpn,2304K,,12288,4252800
	ZPRCENTRY421		zfft_r4dwpn,2304K,,3072,3303936
	ZPRCENTRY421		zfft_r4dwpn,2304K,,2304,3100928
	ZPRCENTRY421		zfft_r4dwpn,2304K,,1024,2066176
	ZPRCENTRY421		zfft_r4dwpn,2304K,,768,2004224
	DD			0
	PRCSTRT	45170000, 2408448, 0.015855
	ZPRCENTRY421		zfft_r4dwpn,2352K,,3136,3362816
	ZPRCENTRY421		zfft_r4dwpn,2352K,,2688,2912640			;,,,SKX
	DD			0
	PRCSTRT	46170000, 2457600, 0.015855
	ZPRCENTRY421		zfft_r4dwpn,2400K,,12800,4412544
	ZPRCENTRY421		zfft_r4dwpn,2400K,,3840,3125888
	ZPRCENTRY421		zfft_r4dwpn,2400K,,3200,3431936
	ZPRCENTRY421		zfft_r4dwpn,2400K,,2560,2286208
	ZPRCENTRY421		zfft_r4dwpn,2400K,,1920,3099648,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2400K,,1600,3019520
	DD			0
	PRCSTRT	48440000, 2580480, 0.0159
	ZPRCENTRY421		zfft_r4dwpn,2520K,,2688,2398848
	ZPRCENTRY421		zfft_r4dwpn,2520K,,2240,2600960
	ZPRCENTRY421		zfft_r4dwpn,2520K,,1920,2257280,,,SKX
	DD			0
	PRCSTRT	49080000, 2621440, 0.016017
	ZPRCENTRY421		zfft_r4dwpn,2560K,,20480,5749888
	ZPRCENTRY421		zfft_r4dwpn,2560K,,4096,3340928
	ZPRCENTRY421		zfft_r4dwpn,2560K,,2560,3434752,,SKX
	DD			0
	PRCSTRT	49860000, 2654208, 0.0165
	ZPRCENTRY421		zfft_r4dwpn,2592K,,2304,2681344,,SKX
	DD			0
	PRCSTRT	51560000, 2752512, 0.017883
	ZPRCENTRY421		zfft_r4dwpn,2688K,,21504,6036608
	ZPRCENTRY421		zfft_r4dwpn,2688K,,3584,3840512
	ZPRCENTRY421		zfft_r4dwpn,2688K,,3072,3332480
	ZPRCENTRY421		zfft_r4dwpn,2688K,,2688,3604736,,SKX
	DD			0
	PRCSTRT	52420000, 2809856, 0.017883
	ZPRCENTRY421		zfft_r4dwpn,2744K,,3136,3391872			;,,SKX
	DD			0
	PRCSTRT	53620000, 2867200, 0.019289
	ZPRCENTRY421		zfft_r4dwpn,2800K,,4480,3643008
	ZPRCENTRY421		zfft_r4dwpn,2800K,,3200,3461504
	ZPRCENTRY421		zfft_r4dwpn,2800K,,2240,3608064			;,,SKX
	DD			0
	PRCSTRT	55240000, 2949120, 0.019245
	ZPRCENTRY421		zfft_r4dwpn,2880K,,15360,5301376,SKX
	ZPRCENTRY421		zfft_r4dwpn,2880K,,4608,3752576
	ZPRCENTRY421		zfft_r4dwpn,2880K,,3840,4112896
	ZPRCENTRY421		zfft_r4dwpn,2880K,,3072,2744960
	ZPRCENTRY421		zfft_r4dwpn,2880K,,2560,2968064
	ZPRCENTRY421		zfft_r4dwpn,2880K,,2304,3717120
	ZPRCENTRY421		zfft_r4dwpn,2880K,,1920,3614464
	DD			0
	PRCSTRT	56200000, 3010560, 0.01927
	ZPRCENTRY421		zfft_r4dwpn,2940K,,3136,2792064
	ZPRCENTRY421		zfft_r4dwpn,2940K,,2240,2624896,,SKX
	DD			0
	PRCSTRT	57410000, 3072000, 0.01927
	ZPRCENTRY421		zfft_r4dwpn,3000K,,3200,2849408
	ZPRCENTRY421		zfft_r4dwpn,3000K,,1600,2750464,,,SKX
	DD			0
	PRCSTRT	57970000, 3096576, 0.0193
	ZPRCENTRY421		zfft_r4dwpn,3024K,,2688,3114496
	ZPRCENTRY421		zfft_r4dwpn,3024K,,2304,2705792			;,,,SKX
	DD			0
	PRCSTRT	58810000, 3145728, 0.019430
	ZPRCENTRY421		zfft_r4dwpn,3M,,24576,6962304
	ZPRCENTRY421		zfft_r4dwpn,3M,,16384,5661824,SKX
	ZPRCENTRY421		zfft_r4dwpn,3M,,4096,4393472
	ZPRCENTRY421		zfft_r4dwpn,3M,,3072,4122880
	ZPRCENTRY421		zfft_r4dwpn,3M,,1024,2639104
	DD			0
	PRCSTRT	59720000, 3211264, 0.024079
	ZPRCENTRY421		zfft_r4dwpn,3136K,,25088,7031936,SKX
	ZPRCENTRY421		zfft_r4dwpn,3136K,,3584,3873152
	ZPRCENTRY421		zfft_r4dwpn,3136K,,3136,4198656
	DD			0
	PRCSTRT	61070000, 3276800, 0.020978
	ZPRCENTRY421		zfft_r4dwpn,3200K,,5120,4160128
	ZPRCENTRY421		zfft_r4dwpn,3200K,,3200,4284672
	ZPRCENTRY421		zfft_r4dwpn,3200K,,2560,4118528
	ZPRCENTRY421		zfft_r4dwpn,3200K,,1600,3880192,,,SKX
	DD			0
	PRCSTRT	64090000, 3440640, 0.023358
	ZPRCENTRY421		zfft_r4dwpn,3360K,,17920,6173824,SKX
	ZPRCENTRY421		zfft_r4dwpn,3360K,,5376,4366976
	ZPRCENTRY421		zfft_r4dwpn,3360K,,4480,4793856
	ZPRCENTRY421		zfft_r4dwpn,3360K,,3840,4147584
	ZPRCENTRY421		zfft_r4dwpn,3360K,,3584,3187328
	ZPRCENTRY421		zfft_r4dwpn,3360K,,2688,4322304
	ZPRCENTRY421		zfft_r4dwpn,3360K,,2560,2994560
	ZPRCENTRY421		zfft_r4dwpn,3360K,,2240,4207360
	DD			0
	PRCSTRT	66110000, 3538944, 0.023145
	ZPRCENTRY421		zfft_r4dwpn,3456K,,18432,6407296,SKX
	ZPRCENTRY421		zfft_r4dwpn,3456K,,4608,4936192
	ZPRCENTRY421		zfft_r4dwpn,3456K,,3072,3561984
	ZPRCENTRY421		zfft_r4dwpn,3456K,,2304,4333312
	DD			0
	PRCSTRT	67300000, 3612672, 0.024
	ZPRCENTRY421		zfft_r4dwpn,3528K,,3136,3625984
	ZPRCENTRY421		zfft_r4dwpn,3528K,,2688,3142016			;,,SKX
	DD			0
	PRCSTRT	68160000, 3670016, 0.024079
	ZPRCENTRY421		zfft_r4dwpn,3584K,,28672,8043648
	ZPRCENTRY421		zfft_r4dwpn,3584K,,4096,4430208
	ZPRCENTRY421		zfft_r4dwpn,3584K,,3584,4794624			;,,,SKX
	DD			0
	PRCSTRT	68730000, 3686400, 0.0242
	ZPRCENTRY421		zfft_r4dwpn,3600K,,3840,3412608
	ZPRCENTRY421		zfft_r4dwpn,3600K,,3200,3700224
	ZPRCENTRY421		zfft_r4dwpn,3600K,,1920,3289088
	ZPRCENTRY421		zfft_r4dwpn,3600K,,1600,3184896,,,SKX
	DD			0
	PRCSTRT	73200000, 3932160, 0.025515
	ZPRCENTRY421		zfft_r4dwpn,3840K,,30720,8617088
	ZPRCENTRY421		zfft_r4dwpn,3840K,,20480,7062656
	ZPRCENTRY421		zfft_r4dwpn,3840K,,6144,4995712
	ZPRCENTRY421		zfft_r4dwpn,3840K,,5120,5474816
	ZPRCENTRY421		zfft_r4dwpn,3840K,,4096,3646080
	ZPRCENTRY421		zfft_r4dwpn,3840K,,3840,5134592
	ZPRCENTRY421		zfft_r4dwpn,3840K,,3072,4941824
	ZPRCENTRY421		zfft_r4dwpn,3840K,,2560,4802304
	ZPRCENTRY421		zfft_r4dwpn,3840K,,1920,4644096,,,SKX
	DD			0
	PRCSTRT	74510000, 4014080, 0.02605
	ZPRCENTRY421		zfft_r4dwpn,3920K,,4480,4833664
	ZPRCENTRY421		zfft_r4dwpn,3920K,,3136,5034496,,SKX
	DD			0
	PRCSTRT	76100000, 4096000, 0.028093
	ZPRCENTRY421		zfft_r4dwpn,4000K,,6400,5204608
	ZPRCENTRY421		zfft_r4dwpn,4000K,,3200,5137408			;,,,SKX
	DD			0
	PRCSTRT	76760000, 4128768, 0.027915
	ZPRCENTRY421		zfft_r4dwpn,4032K,,21504,7414912
	ZPRCENTRY421		zfft_r4dwpn,4032K,,5376,5747200,,SKX
	ZPRCENTRY421		zfft_r4dwpn,4032K,,4608,4977024
	ZPRCENTRY421		zfft_r4dwpn,4032K,,3584,4139520
	ZPRCENTRY421		zfft_r4dwpn,4032K,,3072,3592576
	ZPRCENTRY421		zfft_r4dwpn,4032K,,2688,5039872
	DD			0
	PRCSTRT	77720000, 4194304, 0.026160
	ZPRCENTRY421		zfft_r4dwpn,4M,,32768,9256064
	ZPRCENTRY421		zfft_r4dwpn,4M,,4096,5482752			;,,SKX
	DD			0
	PRCSTRT	78130000, 4214784, 0.027
	ZPRCENTRY421		zfft_r4dwpn,4116K,,3136,3657088			;,,,SKX
	DD			0
	PRCSTRT	79860000, 4300800, 0.027
	ZPRCENTRY421		zfft_r4dwpn,4200K,,4480,3975808
	ZPRCENTRY421		zfft_r4dwpn,4200K,,3200,3731840
	ZPRCENTRY421		zfft_r4dwpn,4200K,,2240,3825664,,,SKX
	DD			0
	PRCSTRT	82310000, 4423680, 0.028
	ZPRCENTRY421		zfft_r4dwpn,4320K,,4608,4094592
	ZPRCENTRY421		zfft_r4dwpn,4320K,,3840,4432384
	ZPRCENTRY421		zfft_r4dwpn,4320K,,2304,3940352
	ZPRCENTRY421		zfft_r4dwpn,4320K,,1920,3808000,,,SKX
	DD			0
	PRCSTRT	84980000, 4587520, 0.030330
	ZPRCENTRY421		zfft_r4dwpn,4480K,,7168,5814912
	ZPRCENTRY421		zfft_r4dwpn,4480K,,5120,5519744
	ZPRCENTRY421		zfft_r4dwpn,4480K,,4480,5984512
	ZPRCENTRY421		zfft_r4dwpn,4480K,,3584,5748736
	ZPRCENTRY421		zfft_r4dwpn,4480K,,2240,5405952,,,SKX
	DD			0
	PRCSTRT	87510000, 4718592, 0.030538
	ZPRCENTRY421		zfft_r4dwpn,4608K,,24576,8537216,SKX
	ZPRCENTRY421		zfft_r4dwpn,4608K,,6144,6572544
	ZPRCENTRY421		zfft_r4dwpn,4608K,,4608,6160640
	ZPRCENTRY421		zfft_r4dwpn,4608K,,4096,4733440
	ZPRCENTRY421		zfft_r4dwpn,4608K,,3072,5760768
	ZPRCENTRY421		zfft_r4dwpn,4608K,,2304,5565696
	DD			0
	PRCSTRT	89040000, 4816896, 0.033388
	ZPRCENTRY421		zfft_r4dwpn,4704K,,25088,8639616,SKX
	ZPRCENTRY421		zfft_r4dwpn,4704K,,5376,5794176
	ZPRCENTRY421		zfft_r4dwpn,4704K,,3584,4174208
	ZPRCENTRY421		zfft_r4dwpn,4704K,,3136,5870336
	DD			0
	PRCSTRT	91010000, 4915200, 0.033479
	ZPRCENTRY421		zfft_r4dwpn,4800K,,7680,6240896
	ZPRCENTRY421		zfft_r4dwpn,4800K,,6400,6846976
	ZPRCENTRY421		zfft_r4dwpn,4800K,,5120,4539008
	ZPRCENTRY421		zfft_r4dwpn,4800K,,3840,6156288
	ZPRCENTRY421		zfft_r4dwpn,4800K,,3200,5990144
	ZPRCENTRY421		zfft_r4dwpn,4800K,,2560,4364288
	ZPRCENTRY421		zfft_r4dwpn,4800K,,1600,4066560,,,SKX
	DD			0
	PRCSTRT	95570000, 5160960, 0.034
	ZPRCENTRY421		zfft_r4dwpn,5040K,,5376,4764288,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5040K,,4480,5164544
	ZPRCENTRY421		zfft_r4dwpn,5040K,,3840,4469120
	ZPRCENTRY421		zfft_r4dwpn,5040K,,2688,4579328
	ZPRCENTRY421		zfft_r4dwpn,5040K,,2240,4429056
	DD			0
	PRCSTRT	96890000, 5242880, 0.034976
	ZPRCENTRY421		zfft_r4dwpn,5M,,8192,6650496
	ZPRCENTRY421		zfft_r4dwpn,5M,,5120,6834432,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5M,,4096,6572032
	ZPRCENTRY421		zfft_r4dwpn,5M,,2560,6169856
	DD			0
	PRCSTRT	98320000, 5308416, 0.036
	ZPRCENTRY421		zfft_r4dwpn,5184K,,4608,5317120
	ZPRCENTRY421		zfft_r4dwpn,5184K,,2304,4560640,,,SKX
	DD			0
	PRCSTRT	101500000, 5505024, 0.038004
	ZPRCENTRY421		zfft_r4dwpn,5376K,,28672,9880704
	ZPRCENTRY421		zfft_r4dwpn,5376K,,7168,7653888
	ZPRCENTRY421		zfft_r4dwpn,5376K,,6144,6625664
	ZPRCENTRY421		zfft_r4dwpn,5376K,,5376,7174400,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5376K,,4096,4772224
	ZPRCENTRY421		zfft_r4dwpn,5376K,,3584,6702848
	ZPRCENTRY421		zfft_r4dwpn,5376K,,2688,6475008
	DD			0
	PRCSTRT	105700000, 5734400, 0.0395
	ZPRCENTRY421		zfft_r4dwpn,5600K,,6400,6902144
	ZPRCENTRY421		zfft_r4dwpn,5600K,,4480,7175168			;,,,SKX
	DD			0
	PRCSTRT	108900000, 5898240, 0.040324
	ZPRCENTRY421		zfft_r4dwpn,5760K,,30720,10585216
	ZPRCENTRY421		zfft_r4dwpn,5760K,,9216,7481984
	ZPRCENTRY421		zfft_r4dwpn,5760K,,7680,8210944
	ZPRCENTRY421		zfft_r4dwpn,5760K,,6144,5448320
	ZPRCENTRY421		zfft_r4dwpn,5760K,,5120,5896704
	ZPRCENTRY421		zfft_r4dwpn,5760K,,4608,7385088
	ZPRCENTRY421		zfft_r4dwpn,5760K,,3840,7177984
	ZPRCENTRY421		zfft_r4dwpn,5760K,,3072,5232640
	ZPRCENTRY421		zfft_r4dwpn,5760K,,2560,5052160
	ZPRCENTRY421		zfft_r4dwpn,5760K,,1920,4861184,,,SKX
	DD			0
	PRCSTRT	110800000, 6021120, 0.040691
	ZPRCENTRY421		zfft_r4dwpn,5880K,,4480,5206400
	ZPRCENTRY421		zfft_r4dwpn,5880K,,3136,5330944			;,,,SKX
	DD			0
	PRCSTRT	113200000, 6144000, 0.0407
	ZPRCENTRY421		zfft_r4dwpn,6000K,,6400,5675648
	ZPRCENTRY421		zfft_r4dwpn,6000K,,3200,5439488			;,,,SKX
	DD			0
	PRCSTRT	114300000, 6193152, 0.041
	ZPRCENTRY421		zfft_r4dwpn,6048K,,5376,6189568,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6048K,,4608,5360000
	ZPRCENTRY421		zfft_r4dwpn,6048K,,2688,5300992
	DD			0
	PRCSTRT	115800000, 6291456, 0.041752
	ZPRCENTRY421		zfft_r4dwpn,6M,,32768,11355264
	ZPRCENTRY421		zfft_r4dwpn,6M,,8192,8751616
	ZPRCENTRY421		zfft_r4dwpn,6M,,6144,8202496,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6M,,4096,7661312
	ZPRCENTRY421		zfft_r4dwpn,6M,,3072,7398656
	DD			0
	PRCSTRT	117800000, 6422528, 0.042000
	ZPRCENTRY421		zfft_r4dwpn,6272K,,7168,7715200,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6272K,,3136,7542016
	DD			0
	PRCSTRT	120400000, 6553600, 0.044473
	ZPRCENTRY421		zfft_r4dwpn,6400K,,10240,8313472,SKX
	ZPRCENTRY421		zfft_r4dwpn,6400K,,6400,8544512
	ZPRCENTRY421		zfft_r4dwpn,6400K,,5120,8194048
	ZPRCENTRY421		zfft_r4dwpn,6400K,,3200,7695616
	DD			0
	PRCSTRT	126400000, 6881280, 0.048606
	ZPRCENTRY421		zfft_r4dwpn,6720K,,7680,8276352
	ZPRCENTRY421		zfft_r4dwpn,6720K,,7168,6341248
	ZPRCENTRY421		zfft_r4dwpn,6720K,,5376,8601600,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6720K,,5120,5943680
	ZPRCENTRY421		zfft_r4dwpn,6720K,,4480,8365824
	ZPRCENTRY421		zfft_r4dwpn,6720K,,3584,6084608
	ZPRCENTRY421		zfft_r4dwpn,6720K,,2240,5653760
	DD			0
	PRCSTRT	130200000, 7077888, 0.050847
	ZPRCENTRY421		zfft_r4dwpn,6912K,,9216,9845248
	ZPRCENTRY421		zfft_r4dwpn,6912K,,6144,7076352
	ZPRCENTRY421		zfft_r4dwpn,6912K,,4608,8609536			;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6912K,,3072,6055680
	ZPRCENTRY421		zfft_r4dwpn,6912K,,2304,5819648
	DD			0
	PRCSTRT	132600000, 7225344, 0.051
	ZPRCENTRY421		zfft_r4dwpn,7056K,,5376,6238592,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7056K,,3136,6170880
	DD			0
	PRCSTRT	134300000, 7340032, 0.052097
	ZPRCENTRY421		zfft_r4dwpn,7M,,8192,8821120
	ZPRCENTRY421		zfft_r4dwpn,7M,,7168,9554176,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7M,,3584,8611072
	DD			0
	PRCSTRT	135400000, 7372800, 0.051
	ZPRCENTRY421		zfft_r4dwpn,7200K,,7680,6804096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7200K,,6400,7371264
	ZPRCENTRY421		zfft_r4dwpn,7200K,,3840,6514688
	ZPRCENTRY421		zfft_r4dwpn,7200K,,3200,6296320
	DD			0
	PRCSTRT	144200000, 7864320, 0.053362
	ZPRCENTRY421		zfft_r4dwpn,7680K,,12288,9968256
	ZPRCENTRY421		zfft_r4dwpn,7680K,,10240,10938880,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7680K,,8192,7250560
	ZPRCENTRY421		zfft_r4dwpn,7680K,,7680,10246400
	ZPRCENTRY421		zfft_r4dwpn,7680K,,6144,9832448
	ZPRCENTRY421		zfft_r4dwpn,7680K,,5120,9553664
	ZPRCENTRY421		zfft_r4dwpn,7680K,,4096,6952960
	ZPRCENTRY421		zfft_r4dwpn,7680K,,3840,9221376
	ZPRCENTRY421		zfft_r4dwpn,7680K,,2560,6448384
	DD			0
	PRCSTRT	149800000, 8192000, 0.059014
	ZPRCENTRY421		zfft_r4dwpn,8000K,,12800,10365568
	ZPRCENTRY421		zfft_r4dwpn,8000K,,6400,10242048		;,,,SKX
	DD			0
	PRCSTRT	151300000, 8257536, 0.063072
	ZPRCENTRY421		zfft_r4dwpn,8064K,,9216,9922944
	ZPRCENTRY421		zfft_r4dwpn,8064K,,7168,8239616
	ZPRCENTRY421		zfft_r4dwpn,8064K,,6144,7131520
	ZPRCENTRY421		zfft_r4dwpn,8064K,,5376,10028800,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8064K,,3584,7042816
	ZPRCENTRY421		zfft_r4dwpn,8064K,,2688,6765824
	DD			0
	PRCSTRT	153000000, 8388608, 0.064
	ZPRCENTRY421		zfft_r4dwpn,8M,,8192,10922240
	ZPRCENTRY421		zfft_r4dwpn,8M,,4096,9839872			;,,,SKX
	DD			0
	PRCSTRT	157000000, 8601600, 0.066
	ZPRCENTRY421		zfft_r4dwpn,8400K,,6400,7428480,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8400K,,4480,7589888
	DD			0
	PRCSTRT	162000000, 8847360, 0.067
	ZPRCENTRY421		zfft_r4dwpn,8640K,,9216,8155776
	ZPRCENTRY421		zfft_r4dwpn,8640K,,7680,8837632,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8640K,,4608,7811072
	ZPRCENTRY421		zfft_r4dwpn,8640K,,3840,7540480
	DD			0
	PRCSTRT	167100000, 9175040, 0.068512
	ZPRCENTRY421		zfft_r4dwpn,8960K,,10240,11024768,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8960K,,7168,11454464
	ZPRCENTRY421		zfft_r4dwpn,8960K,,4480,10747136
	DD			0
	PRCSTRT	172500000, 9437184, 0.067985
	ZPRCENTRY421		zfft_r4dwpn,9M,,12288,13117952
	ZPRCENTRY421		zfft_r4dwpn,9M,,9216,12286208
	ZPRCENTRY421		zfft_r4dwpn,9M,,8192,9419264
	ZPRCENTRY421		zfft_r4dwpn,9M,,6144,11462400
	ZPRCENTRY421		zfft_r4dwpn,9M,,4608,11058432			;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,9M,,4096,8046336
	ZPRCENTRY421		zfft_r4dwpn,9M,,3072,7726336
	DD			0
	PRCSTRT	175200000, 9633792, 0.068994
	ZPRCENTRY421		zfft_r4dwpn,9408K,,7168,8302976			;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,9408K,,3136,7875840
	DD			0
	PRCSTRT	179000000, 9830400, 0.070403
	ZPRCENTRY421		zfft_r4dwpn,9600K,,15360,12442240,SKX
	ZPRCENTRY421		zfft_r4dwpn,9600K,,12800,13646336
	ZPRCENTRY421		zfft_r4dwpn,9600K,,10240,9060992
	ZPRCENTRY421		zfft_r4dwpn,9600K,,7680,12281856
	ZPRCENTRY421		zfft_r4dwpn,9600K,,6400,11939584
	ZPRCENTRY421		zfft_r4dwpn,9600K,,5120,8665088
	ZPRCENTRY421		zfft_r4dwpn,9600K,,3200,8035584
	DD			0
	PRCSTRT	188000000, 10321920, 0.075
	ZPRCENTRY421		zfft_r4dwpn,10080K,,7680,8905088
	ZPRCENTRY421		zfft_r4dwpn,10080K,,5376,9095168		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,10080K,,4480,8784640
	DD			0
	PRCSTRT	190500000, 10485760, 0.076827
	ZPRCENTRY421		zfft_r4dwpn,10M,,16384,13277824,SKX
	ZPRCENTRY421		zfft_r4dwpn,10M,,10240,13650176
	ZPRCENTRY421		zfft_r4dwpn,10M,,8192,13092864
	ZPRCENTRY421		zfft_r4dwpn,10M,,5120,12272896
	DD			0
	PRCSTRT	193600000, 10616832, 0.077
	ZPRCENTRY421		zfft_r4dwpn,10368K,,9216,10594816,,SKX
	ZPRCENTRY421		zfft_r4dwpn,10368K,,4608,9039616
	DD			0
	PRCSTRT	200100000, 11010048, 0.083177
	ZPRCENTRY421		zfft_r4dwpn,10752K,,12288,13220224
	ZPRCENTRY421		zfft_r4dwpn,10752K,,8192,9490816
	ZPRCENTRY421		zfft_r4dwpn,10752K,,7168,13354752		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,10752K,,5376,12883200
	ZPRCENTRY421		zfft_r4dwpn,10752K,,3584,8987904
	DD			0
	PRCSTRT	207600000, 11468800, 0.085770
	ZPRCENTRY421		zfft_r4dwpn,11200K,,17920,14502528,SKX
	ZPRCENTRY421		zfft_r4dwpn,11200K,,12800,13752704
	DD			0
	PRCSTRT	214300000, 11796480, 0.086096
	ZPRCENTRY421		zfft_r4dwpn,11520K,,18432,14973568
	ZPRCENTRY421		zfft_r4dwpn,11520K,,15360,16378368,,SKX
	ZPRCENTRY421		zfft_r4dwpn,11520K,,12288,10863232
	ZPRCENTRY421		zfft_r4dwpn,11520K,,10240,11770368
	ZPRCENTRY421		zfft_r4dwpn,11520K,,9216,14727168
	ZPRCENTRY421		zfft_r4dwpn,11520K,,7680,14317312
	ZPRCENTRY421		zfft_r4dwpn,11520K,,6144,10393600
	ZPRCENTRY421		zfft_r4dwpn,11520K,,5120,10028800
	ZPRCENTRY421		zfft_r4dwpn,11520K,,3840,9622784
	DD			0
	PRCSTRT	223000000, 12288000, 0.0917
	ZPRCENTRY421		zfft_r4dwpn,12000K,,12800,11297408
	ZPRCENTRY421		zfft_r4dwpn,12000K,,6400,10825728		;,,,SKX
	DD			0
	PRCSTRT	225100000, 12386304, 0.092
	ZPRCENTRY421		zfft_r4dwpn,12096K,,9216,10674560		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,12096K,,5376,10526464
	DD			0
	PRCSTRT	228300000, 12582912, 0.094441
	ZPRCENTRY421		zfft_r4dwpn,12M,,16384,17476096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,12M,,12288,16369920
	ZPRCENTRY421		zfft_r4dwpn,12M,,8192,15263488
	ZPRCENTRY421		zfft_r4dwpn,12M,,6144,14722304
	ZPRCENTRY421		zfft_r4dwpn,12M,,4096,10265856
	DD			0
	PRCSTRT	237200000, 13107200, 0.094766
	ZPRCENTRY421		zfft_r4dwpn,12800K,,20480,16579200
	ZPRCENTRY421		zfft_r4dwpn,12800K,,12800,17033472
	ZPRCENTRY421		zfft_r4dwpn,12800K,,10240,16361472,,SKX
	ZPRCENTRY421		zfft_r4dwpn,12800K,,6400,15334656
	DD			0
	PRCSTRT	248700000, 13762560, 0.106540
	ZPRCENTRY421		zfft_r4dwpn,13440K,,21504,17406592
	ZPRCENTRY421		zfft_r4dwpn,13440K,,17920,19094016
	ZPRCENTRY421		zfft_r4dwpn,13440K,,15360,16505216,,SKX
	ZPRCENTRY421		zfft_r4dwpn,13440K,,10240,11858304
	ZPRCENTRY421		zfft_r4dwpn,13440K,,7168,12105728
	ZPRCENTRY421		zfft_r4dwpn,13440K,,4480,11209984
	DD			0
	PRCSTRT	256600000, 14155776, 0.111224
	ZPRCENTRY421		zfft_r4dwpn,13824K,,18432,19696128		;,,SKX
	ZPRCENTRY421		zfft_r4dwpn,13824K,,12288,14113280
	ZPRCENTRY421		zfft_r4dwpn,13824K,,9216,17168128
	ZPRCENTRY421		zfft_r4dwpn,13824K,,6144,12027648
	ZPRCENTRY421		zfft_r4dwpn,13824K,,4608,11533568
	DD			0
	PRCSTRT	264700000, 14680064, 0.114093
	ZPRCENTRY421		zfft_r4dwpn,14M,,16384,17611136			;,,SKX
	ZPRCENTRY421		zfft_r4dwpn,14M,,7168,17155328
	DD			0
	PRCSTRT	266700000, 14745600, 0.116
	ZPRCENTRY421		zfft_r4dwpn,14400K,,15360,13558400,,SKX
	ZPRCENTRY421		zfft_r4dwpn,14400K,,12800,14682624
	ZPRCENTRY421		zfft_r4dwpn,14400K,,7680,12978176
	ZPRCENTRY421		zfft_r4dwpn,14400K,,6400,12527360
	DD			0
	PRCSTRT	283700000, 15728640, 0.120880
	ZPRCENTRY421		zfft_r4dwpn,15M,,24576,19954304
	ZPRCENTRY421		zfft_r4dwpn,15M,,20480,21826048
	ZPRCENTRY421		zfft_r4dwpn,15M,,16384,14467712
	ZPRCENTRY421		zfft_r4dwpn,15M,,15360,20441344,,SKX
	ZPRCENTRY421		zfft_r4dwpn,15M,,12288,19621888
	ZPRCENTRY421		zfft_r4dwpn,15M,,10240,19072768
	ZPRCENTRY421		zfft_r4dwpn,15M,,8192,13834240
	ZPRCENTRY421		zfft_r4dwpn,15M,,7680,18388224
	ZPRCENTRY421		zfft_r4dwpn,15M,,5120,12797184
	DD			0
	PRCSTRT	288500000, 16056320, 0.120880
	ZPRCENTRY421		zfft_r4dwpn,15680K,,25088,20294272,SKX
	ZPRCENTRY421		zfft_r4dwpn,15680K,,17920,19241344
	DD			0
	PRCSTRT	295000000, 16384000, 0.127481
	ZPRCENTRY421		zfft_r4dwpn,16000K,,12800,20420608		;,,,SKX
	DD			0
	PRCSTRT	297800000, 16515072, 0.132344
	ZPRCENTRY421		zfft_r4dwpn,16128K,,21504,22915584
	ZPRCENTRY421		zfft_r4dwpn,16128K,,18432,19847552,,SKX
	ZPRCENTRY421		zfft_r4dwpn,16128K,,12288,14217600
	ZPRCENTRY421		zfft_r4dwpn,16128K,,7168,14010112
	ZPRCENTRY421		zfft_r4dwpn,16128K,,5376,13432064
	DD			0
	PRCSTRT	302000000, 16777216, 0.129355
	ZPRCENTRY421		zfft_r4dwpn,16M,,16384,21809408,,SKX
	ZPRCENTRY421		zfft_r4dwpn,16M,,8192,19604736
	DD			0
	PRCSTRT	309100000, 17203200, 0.129355
	ZPRCENTRY421		zfft_r4dwpn,16800K,,17920,15803008,,SKX
	ZPRCENTRY421		zfft_r4dwpn,16800K,,12800,14791040
	DD			0
	PRCSTRT	318900000, 17694720, 0.141
	ZPRCENTRY421		zfft_r4dwpn,17280K,,18432,16310912
	ZPRCENTRY421		zfft_r4dwpn,17280K,,15360,17619456,,SKX
	ZPRCENTRY421		zfft_r4dwpn,17280K,,9216,15558656
	ZPRCENTRY421		zfft_r4dwpn,17280K,,7680,15017728
	DD			0
	PRCSTRT	328900000, 18350080, 0.147854
	ZPRCENTRY421		zfft_r4dwpn,17920K,,28672,23198336
	ZPRCENTRY421		zfft_r4dwpn,17920K,,20480,21993856
	ZPRCENTRY421		zfft_r4dwpn,17920K,,17920,23832832,,SKX
	DD			0
	PRCSTRT	339700000, 18874368, 0.149016
	ZPRCENTRY421		zfft_r4dwpn,18M,,24576,26249728
	ZPRCENTRY421		zfft_r4dwpn,18M,,18432,24570112
	ZPRCENTRY421		zfft_r4dwpn,18M,,16384,18799104,,SKX
	ZPRCENTRY421		zfft_r4dwpn,18M,,12288,22873856
	ZPRCENTRY421		zfft_r4dwpn,18M,,9216,22050048
	ZPRCENTRY421		zfft_r4dwpn,18M,,8192,16008960
	ZPRCENTRY421		zfft_r4dwpn,18M,,6144,15344896
	DD			0
	PRCSTRT	345300000, 19267584, 0.149016
	ZPRCENTRY421		zfft_r4dwpn,18816K,,25088,26720768,SKX
	ZPRCENTRY421		zfft_r4dwpn,18816K,,21504,23091584
	DD			0
	PRCSTRT	353000000, 19660800, 0.153223
	ZPRCENTRY421		zfft_r4dwpn,19200K,,30720,24853120
	ZPRCENTRY421		zfft_r4dwpn,19200K,,20480,18064000
	ZPRCENTRY421		zfft_r4dwpn,19200K,,15360,24504320,,SKX
	ZPRCENTRY421		zfft_r4dwpn,19200K,,12800,23807744
	ZPRCENTRY421		zfft_r4dwpn,19200K,,10240,17283072
	ZPRCENTRY421		zfft_r4dwpn,19200K,,6400,15981824
	DD			0
	PRCSTRT	370200000, 20643840, 0.168
	ZPRCENTRY421		zfft_r4dwpn,20160K,,21504,18965120
	ZPRCENTRY421		zfft_r4dwpn,20160K,,17920,20539904
	ZPRCENTRY421		zfft_r4dwpn,20160K,,15360,17748352,,SKX
	DD			0
	PRCSTRT	375600000, 20971520, 0.168125
	ZPRCENTRY421		zfft_r4dwpn,20M,,32768,26573440
	ZPRCENTRY421		zfft_r4dwpn,20M,,20480,27240704
	ZPRCENTRY421		zfft_r4dwpn,20M,,16384,26142720,,SKX
	ZPRCENTRY421		zfft_r4dwpn,20M,,10240,24495360
	DD			0
	PRCSTRT	381400000, 21233664, 0.171
	ZPRCENTRY421		zfft_r4dwpn,20736K,,18432,21182976,,SKX
	ZPRCENTRY421		zfft_r4dwpn,20736K,,9216,18003712
	DD			0
	PRCSTRT	393900000, 22020096, 0.179899
	ZPRCENTRY421		zfft_r4dwpn,21M,,28672,30542336
	ZPRCENTRY421		zfft_r4dwpn,21M,,24576,26450304
	ZPRCENTRY421		zfft_r4dwpn,21M,,21504,28600576
	ZPRCENTRY421		zfft_r4dwpn,21M,,16384,18936192,,SKX
	ZPRCENTRY421		zfft_r4dwpn,21M,,7168,17876224
	DD			0
	PRCSTRT	400200000, 22478848, 0.194140
	ZPRCENTRY421		zfft_r4dwpn,21952K,,25088,26925440,,SKX
	DD			0
	PRCSTRT	409300000, 22937600, 0.186258
	ZPRCENTRY421		zfft_r4dwpn,22400K,,17920,28571648,,SKX
	DD			0
	PRCSTRT	422000000, 23592960, 0.194140
	ZPRCENTRY421		zfft_r4dwpn,23040K,,30720,32721408
	ZPRCENTRY421		zfft_r4dwpn,23040K,,24576,21734016
	ZPRCENTRY421		zfft_r4dwpn,23040K,,20480,23476736
	ZPRCENTRY421		zfft_r4dwpn,23040K,,18432,29444096
	ZPRCENTRY421		zfft_r4dwpn,23040K,,15360,28567296		;,,SKX
	ZPRCENTRY421		zfft_r4dwpn,23040K,,12288,20723712
	ZPRCENTRY421		zfft_r4dwpn,23040K,,10240,19998464
	ZPRCENTRY421		zfft_r4dwpn,23040K,,7680,19158272
	DD			0
	PRCSTRT	430000000, 24084480, 0.199
	ZPRCENTRY421		zfft_r4dwpn,23520K,,25088,22110848
	ZPRCENTRY421		zfft_r4dwpn,23520K,,17920,20689280,,SKX
	DD			0
	PRCSTRT	438600000, 24576000, 0.200
	ZPRCENTRY421		zfft_r4dwpn,24000K,,12800,21567488		;,,,SKX
	DD			0
	PRCSTRT	442700000, 24772608, 0.201
	ZPRCENTRY421		zfft_r4dwpn,24192K,,21504,24648192
	ZPRCENTRY421		zfft_r4dwpn,24192K,,18432,21336448,,SKX
	DD			0
	PRCSTRT	448400000, 25165824, 0.203308
	ZPRCENTRY421		zfft_r4dwpn,24M,,32768,34966016
	ZPRCENTRY421		zfft_r4dwpn,24M,,24576,32745728,,SKX
	ZPRCENTRY421		zfft_r4dwpn,24M,,16384,30476032
	ZPRCENTRY421		zfft_r4dwpn,24M,,12288,29377792
	ZPRCENTRY421		zfft_r4dwpn,24M,,8192,20423936
	DD			0
	PRCSTRT	456100000, 25690112, 0.209253
	ZPRCENTRY421		zfft_r4dwpn,25088K,,28672,30775680
	ZPRCENTRY421		zfft_r4dwpn,25088K,,25088,33351936,,SKX
	DD			0
	PRCSTRT	466200000, 26214400, 0.208601
	ZPRCENTRY421		zfft_r4dwpn,25M,,20480,32655360,,SKX
	ZPRCENTRY421		zfft_r4dwpn,25M,,12800,30582016
	DD			0
	PRCSTRT	489500000, 27525120, 0.229014
	ZPRCENTRY421		zfft_r4dwpn,26880K,,30720,32971136
	ZPRCENTRY421		zfft_r4dwpn,26880K,,28672,25272960
	ZPRCENTRY421		zfft_r4dwpn,26880K,,21504,34285568
	ZPRCENTRY421		zfft_r4dwpn,26880K,,20480,23646592
	ZPRCENTRY421		zfft_r4dwpn,26880K,,17920,33310464,,SKX
	DD			0
	PRCSTRT	504600000, 28311552, 0.23830
	ZPRCENTRY421		zfft_r4dwpn,27M,,24576,28228096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,27M,,18432,34318080
	ZPRCENTRY421		zfft_r4dwpn,27M,,12288,23979776
	ZPRCENTRY421		zfft_r4dwpn,27M,,9216,22967552
	DD			0
	PRCSTRT	513500000, 28901376, 0.245
	ZPRCENTRY421		zfft_r4dwpn,28224K,,25088,28740096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,28224K,,21504,24826240
	DD			0
	PRCSTRT	519800000, 29360128, 0.247130
	ZPRCENTRY421		zfft_r4dwpn,28M,,32768,35232128
	ZPRCENTRY421		zfft_r4dwpn,28M,,28672,38119680			;,,,SKX
	DD			0
	PRCSTRT	524500000, 29491200, 0.249
	ZPRCENTRY421		zfft_r4dwpn,28800K,,30720,27075200
	ZPRCENTRY421		zfft_r4dwpn,28800K,,15360,25876480,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,28800K,,12800,24958720
	DD			0
	PRCSTRT	557900000, 31457280, 0.258182
	ZPRCENTRY421		zfft_r4dwpn,30M,,30720,40839424
	ZPRCENTRY421		zfft_r4dwpn,30M,,32768,28942976
	ZPRCENTRY421		zfft_r4dwpn,30M,,24576,39241728,,SKX
	ZPRCENTRY421		zfft_r4dwpn,30M,,20480,38070016
	ZPRCENTRY421		zfft_r4dwpn,30M,,16384,27604992
	ZPRCENTRY421		zfft_r4dwpn,30M,,15360,36693248
	ZPRCENTRY421		zfft_r4dwpn,30M,,10240,25511168
	DD			0
	PRCSTRT	567400000, 32112640, 0.2700
	ZPRCENTRY421		zfft_r4dwpn,31360K,,25088,39983104,,SKX
	DD			0
	PRCSTRT	585200000, 33030144, 0.2747
	ZPRCENTRY421		zfft_r4dwpn,32256K,,28672,32848384
	ZPRCENTRY421		zfft_r4dwpn,32256K,,24576,28430720,,SKX
	ZPRCENTRY421		zfft_r4dwpn,32256K,,21504,39970560
	DD			0
	PRCSTRT	594100000, 33554432, 0.279123
	ZPRCENTRY421		zfft_r4dwpn,32M,,32768,43624704
	ZPRCENTRY421		zfft_r4dwpn,32M,,16384,39142656			;,,,SKX
	DD			0
	PRCSTRT	594700000, 33718272, 0.280
	ZPRCENTRY421		zfft_r4dwpn,32928K,,25088,28946816,,SKX
	DD			0
	PRCSTRT	608300000, 34406400, 0.285
	ZPRCENTRY421		zfft_r4dwpn,33600K,,17920,30169088		;,,,SKX
	DD			0
	PRCSTRT	627800000, 35389440, 0.290
	ZPRCENTRY421		zfft_r4dwpn,34560K,,30720,35191296
	ZPRCENTRY421		zfft_r4dwpn,34560K,,18432,31086592
	ZPRCENTRY421		zfft_r4dwpn,34560K,,15360,29943552,,,SKX
	DD			0
	PRCSTRT	647200000, 36700160, 0.306
	ZPRCENTRY421		zfft_r4dwpn,35M,,28672,45697024
	ZPRCENTRY421		zfft_r4dwpn,35M,,17920,42788096			;,,,SKX
	DD			0
	PRCSTRT	668000000, 37748736, 0.310
	ZPRCENTRY421		zfft_r4dwpn,36M,,32768,37599744
	ZPRCENTRY421		zfft_r4dwpn,36M,,24576,45737728,,SKX
	ZPRCENTRY421		zfft_r4dwpn,36M,,18432,44066048
	ZPRCENTRY421		zfft_r4dwpn,36M,,16384,31942400
	ZPRCENTRY421		zfft_r4dwpn,36M,,12288,30590208
	DD			0
	PRCSTRT	679200000, 38535168, 0.328
	ZPRCENTRY421		zfft_r4dwpn,37632K,,28672,33083776
	ZPRCENTRY421		zfft_r4dwpn,37632K,,25088,46614272,,SKX
	DD			0
	PRCSTRT	693600000, 39321600, 0.328
	ZPRCENTRY421		zfft_r4dwpn,38400K,,30720,48957440
	ZPRCENTRY421		zfft_r4dwpn,38400K,,20480,34478080,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,38400K,,12800,31843584
	DD			0
	PRCSTRT	728500000, 41287680, 0.340
	ZPRCENTRY421		zfft_r4dwpn,40320K,,30720,35443072
	ZPRCENTRY421		zfft_r4dwpn,40320K,,21504,36198400
	ZPRCENTRY421		zfft_r4dwpn,40320K,,17920,34912000,,,SKX
	DD			0
	PRCSTRT	737700000, 41943040, 0.350
	ZPRCENTRY421		zfft_r4dwpn,40M,,32768,52283392
	ZPRCENTRY421		zfft_r4dwpn,40M,,20480,48899328,,,SKX
	DD			0
	PRCSTRT	751000000, 42467328, 0.360
	ZPRCENTRY421		zfft_r4dwpn,41472K,,18432,35964672,,,SKX
	DD			0
	PRCSTRT	775000000, 44040192, 0.370
	ZPRCENTRY421		zfft_r4dwpn,42M,,32768,37867904
	ZPRCENTRY421		zfft_r4dwpn,42M,,28672,53274368
	ZPRCENTRY421		zfft_r4dwpn,42M,,21504,51340544,,,SKX
	DD			0
	PRCSTRT	831500000, 47185920, 0.395
	ZPRCENTRY421		zfft_r4dwpn,45M,,30720,57075456
	ZPRCENTRY421		zfft_r4dwpn,45M,,24576,41424896
	ZPRCENTRY421		zfft_r4dwpn,45M,,20480,39896832
	ZPRCENTRY421		zfft_r4dwpn,45M,,15360,38200576,,,SKX
	DD			0
	PRCSTRT	845200000, 48168960, 0.400
	ZPRCENTRY421		zfft_r4dwpn,47040K,,25088,42211328,,,SKX
	DD			0
	PRCSTRT	872700000, 49545216, 0.410
	ZPRCENTRY421		zfft_r4dwpn,48384K,,21504,41887488,,,SKX
	DD			0
	PRCSTRT	884000000, 50331648, 0.420
	ZPRCENTRY421		zfft_r4dwpn,48M,,32768,60942080
	ZPRCENTRY421		zfft_r4dwpn,48M,,24576,58729728,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,48M,,16384,40748288
	DD			0
	PRCSTRT	898800000, 51380224, 0.430
	ZPRCENTRY421		zfft_r4dwpn,49M,,25088,59876608,,,SKX
	DD			0
	PRCSTRT	964700000, 55050240, 0.460
	ZPRCENTRY421		zfft_r4dwpn,53760K,,28672,48240640
	ZPRCENTRY421		zfft_r4dwpn,53760K,,17920,44541184,,,SKX
	DD			0
	PRCSTRT	995600000, 56623104, 0.480
	ZPRCENTRY421		zfft_r4dwpn,55296K,,24576,47924992,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,55296K,,18432,45868288
	DD			0
	PRCSTRT	1013000000, 57802752, 0.520
	ZPRCENTRY421		zfft_r4dwpn,56448K,,25088,48846592,,,SKX
	DD			0
	PRCSTRT	1025000000, 58720256, 0.530
	ZPRCENTRY421		zfft_r4dwpn,56M,,28672,68429056,,,SKX
	DD			0
	PRCSTRT	1036000000, 58982400, 0.540
	ZPRCENTRY421		zfft_r4dwpn,57600K,,30720,51681280		;,,,SKX
	DD			0
	PRCSTRT	1100000000, 62914560, 0.580
	ZPRCENTRY421		zfft_r4dwpn,60M,,32768,55187456
	ZPRCENTRY421		zfft_r4dwpn,60M,,30720,73311488
	ZPRCENTRY421		zfft_r4dwpn,60M,,20480,50898176,,,SKX
	DD			0
	PRCSTRT	1169000000, 67108864, 0.640
	ZPRCENTRY421		zfft_r4dwpn,64M,,32768,78259456,,,SKX
	DD			0
	DD	0

zjmptablep DD	0
	org	$-4
	PRCSTRT	22397,	1024,	0.0000036
	ZPRCENTRY1		zfft_r4dwpn,1K,_ac,16832,SKX
	DD			0
	PRCSTRT	33371,	1536,	0.0000068
	ZPRCENTRY1		zfft_r4dwpn,1536,_ac,24576,SKX
	DD			0
	PRCSTRT	44269,	2048,	0.0000083
	ZPRCENTRY1		zfft_r4dwpn,2K,_ac,32320,SKX
	DD			0
	PRCSTRT	98500,	4608,	0.000019
	ZPRCENTRY1		zfft_r4dwpn,4608,_ac,74112,SKX
	DD			0
	PRCSTRT	109000,	5120,	0.000021
	ZPRCENTRY1		zfft_r4dwpn,5K,_ac,78784,SKX
	DD			0
	PRCSTRT	130000,	6144,	0.000025
	ZPRCENTRY1		zfft_r4dwpn,6K,_ac,98368,SKX
	DD			0
	PRCSTRT	151000,	7168,	0.000030
	ZPRCENTRY1		zfft_r4dwpn,7K,_ac,109760,SKX
	DD			0
	PRCSTRT	162100,	7680,	0.000032
	ZPRCENTRY1		zfft_r4dwpn,7680,_ac,117504,SKX
	DD			0
	PRCSTRT	172300,	8192,	0.000033
	ZPRCENTRY1		zfft_r4dwpn,8K,_ac,129344
	ZPRCENTRY421		zfft_r4dwpn,8K,_ac,64,20736,,,SKX
	DD			0
	PRCSTRT	194100,	9216,	0.000038
	ZPRCENTRY1		zfft_r4dwpn,9K,_ac,143808			;,SKX
	DD			0
	PRCSTRT	214800,	10240,	0.000046
	ZPRCENTRY1		zfft_r4dwpn,10K,_ac,156224,SKX
	DD			0
	PRCSTRT	225500,	10752,	0.000048
	ZPRCENTRY1		zfft_r4dwpn,10752,_ac,163968,SKX
	DD			0
	PRCSTRT	256900,	12288,	0.000054
	ZPRCENTRY1		zfft_r4dwpn,12K,_ac,191296,
	ZPRCENTRY421		zfft_r4dwpn,12K,_ac,64,26880,,,SKX
	DD			0
	PRCSTRT	267200,	12800,	0.000056
	ZPRCENTRY1		zfft_r4dwpn,12800,_ac,200064,SKX
	DD			0
	PRCSTRT	298300,	14336,	0.000063
	ZPRCENTRY1		zfft_r4dwpn,14K,_ac,218176,SKX
	DD			0
	PRCSTRT	320300,	15360,	0.000068
	ZPRCENTRY1		zfft_r4dwpn,15K,_ac,239808,SKX
	DD			0
	PRCSTRT	340000,	16384,	0.000073
	ZPRCENTRY1		zfft_r4dwpn,16K,_ac,253248,SKX
	DD			0
	PRCSTRT	382100,	18432,	0.000088
	ZPRCENTRY1		zfft_r4dwpn,18K,_ac,286272,SKX
	DD			0
	PRCSTRT	423400,	20480,	0.000095
	ZPRCENTRY1		zfft_r4dwpn,20K,_ac,319296,SKX
	DD			0
	PRCSTRT	507100,	24576,	0.000109
	ZPRCENTRY1		zfft_r4dwpn,24K,_ac,381248,SKX
	DD			0
	PRCSTRT	526700,	25600,	0.000113
	ZPRCENTRY1		zfft_r4dwpn,25K,_ac,388544,SKX
	DD			0
	PRCSTRT	630800,	30720,	0.000130
	ZPRCENTRY1		zfft_r4dwpn,30K,_ac,470080			;,SKX
	DD			0
	PRCSTRT	671900,	32768,	0.000140
	ZPRCENTRY1		zfft_r4dwpn,32K,_ac,505152,SKX
	DD			0
	PRCSTRT	732200,	35840,	0.000150
	ZPRCENTRY1		zfft_r4dwpn,35K,_ac,543424			;,SKX
	DD			0
	PRCSTRT	755800,	36864,	0.000155
	ZPRCENTRY1		zfft_r4dwpn,36K,_ac,587584			;,SKX
	DD			0
	PRCSTRT	835000, 40960,	0.000191
	ZPRCENTRY1		zfft_r4dwpn,40K,_ac,624960
	ZPRCENTRY421		zfft_r4dwpn,40K,_ac,320,78080,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,40K,_ac,64,54528
	DD			0
	PRCSTRT	876400, 43008,	0.000191
	ZPRCENTRY1		zfft_r4dwpn,42K,_ac,655936			;,SKX
	DD			0
	PRCSTRT	1000000, 49152,	0.000240
	ZPRCENTRY1		zfft_r4dwpn,48K,_ac,781632
	ZPRCENTRY421		zfft_r4dwpn,48K,_ac,384,92928,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,48K,_ac,64,66816
	DD			0
	PRCSTRT	1161000, 57344, 0.000280
	ZPRCENTRY421		zfft_r4dwpn,56K,_ac,448,106752,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,56K,_ac,64,71424
	DD			0
	PRCSTRT	1243000, 61440,	0.000333
	ZPRCENTRY421		zfft_r4dwpn,60K,_ac,320,100608,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,60K,_ac,64,69376
	DD			0
	PRCSTRT	1324000, 65536,	0.000314
	ZPRCENTRY421		zfft_r4dwpn,64K,_ac,512,121600,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,64K,_ac,64,83712
	DD			0
	PRCSTRT	1488000, 73728,	0.000374
	ZPRCENTRY421		zfft_r4dwpn,72K,_ac,384,119552,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,72K,_ac,64,84224
	DD			0
	PRCSTRT	1648000, 81920,	0.000413
	ZPRCENTRY421		zfft_r4dwpn,80K,_ac,640,150272,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,80K,_ac,64,100608
	DD			0
	PRCSTRT	1731000, 86016, 0.000443
	ZPRCENTRY421		zfft_r4dwpn,84K,_ac,448,137472,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,84K,_ac,64,90880
	DD			0
	PRCSTRT	1974000, 98304, 0.000479
	ZPRCENTRY421		zfft_r4dwpn,96K,_ac,768,178944,,SKX
	ZPRCENTRY421		zfft_r4dwpn,96K,_ac,512,156416
	ZPRCENTRY421		zfft_r4dwpn,96K,_ac,64,117504
	DD			0
	PRCSTRT	2455000, 122880, 0.000652
	ZPRCENTRY421		zfft_r4dwpn,120K,_ac,640,193280,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,120K,_ac,64,126720
	DD			0
	PRCSTRT	2615000, 131072, 0.000622
	ZPRCENTRY421		zfft_r4dwpn,128K,_ac,1024,236288,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,128K,_ac,64,151296
	DD			0
	PRCSTRT	2939000, 147456, 0.000841
	ZPRCENTRY421		zfft_r4dwpn,144K,_ac,768,230144,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,144K,_ac,64,147712
	DD			0
	PRCSTRT	3894000, 196608, 0.001030
	ZPRCENTRY421		zfft_r4dwpn,192K,_ac,1024,303872,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,192K,_ac,64,190208
	DD	0
	PRCSTRT	4037000, 204800, 0.001191
	ZPRCENTRY421		zfft_r4dwpn,200K,_ac,1600,364800
	ZPRCENTRY421		zfft_r4dwpn,200K,_ac,320,181504,,,SKX
	DD	0
	PRCSTRT	4839000, 245760, 0.001449
	ZPRCENTRY421		zfft_r4dwpn,240K,_ac,1920,436992
	ZPRCENTRY421		zfft_r4dwpn,240K,_ac,384,213760,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,240K,_ac,320,226560
	DD	0
	PRCSTRT	5618000, 286720, 0.001668
	ZPRCENTRY421		zfft_r4dwpn,280K,_ac,2240,508160
	ZPRCENTRY421		zfft_r4dwpn,280K,_ac,448,244992,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,280K,_ac,320,233216
	DD			0
	PRCSTRT	5789000, 294912, 0.001632
	ZPRCENTRY421		zfft_r4dwpn,288K,_ac,2304,526080
	ZPRCENTRY421		zfft_r4dwpn,288K,_ac,384,267008,,,SKX
	DD			0
	PRCSTRT	6024000, 307200, 0.001679
	ZPRCENTRY421		zfft_r4dwpn,300K,_ac,1600,469248
	ZPRCENTRY421		zfft_r4dwpn,300K,_ac,320,214784,,,SKX
	DD			0
	PRCSTRT	6403000, 327680, 0.001802
	ZPRCENTRY421		zfft_r4dwpn,320K,_ac,2560,580352
	ZPRCENTRY421		zfft_r4dwpn,320K,_ac,512,277248,,SKX
	ZPRCENTRY421		zfft_r4dwpn,320K,_ac,320,278272
	DD			0
	PRCSTRT	6727000, 344064, 0.001969
	ZPRCENTRY421		zfft_r4dwpn,336K,_ac,2688,609024
	ZPRCENTRY421		zfft_r4dwpn,336K,_ac,448,306432
	ZPRCENTRY421		zfft_r4dwpn,336K,_ac,384,274176,,,SKX
	DD			0
	PRCSTRT	7219000, 368640, 0.001969
	ZPRCENTRY421		zfft_r4dwpn,360K,_ac,1920,561920
	ZPRCENTRY421		zfft_r4dwpn,360K,_ac,384,251648,,SKX
	ZPRCENTRY421		zfft_r4dwpn,360K,_ac,320,264448
	DD			0
	PRCSTRT	7672000, 393216, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,384K,_ac,3072,699136
	ZPRCENTRY421		zfft_r4dwpn,384K,_ac,512,346880,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,384K,_ac,384,327424
	DD			0
	PRCSTRT	7810000, 401408, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,392K,_ac,3136,708864
	ZPRCENTRY421		zfft_r4dwpn,392K,_ac,448,314112,,,SKX
	DD			0
	PRCSTRT	7975000, 409600, 0.002654
	ZPRCENTRY421		zfft_r4dwpn,400K,_ac,3200,723712
	ZPRCENTRY421		zfft_r4dwpn,400K,_ac,640,340736
	ZPRCENTRY421		zfft_r4dwpn,400K,_ac,320,329984,,,SKX
	DD			0
	PRCSTRT	8381000, 430080, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,420K,_ac,2240,653568
	ZPRCENTRY421		zfft_r4dwpn,420K,_ac,448,287488
	ZPRCENTRY421		zfft_r4dwpn,420K,_ac,320,273152,,,SKX
	DD			0
	PRCSTRT	8636000, 442368, 0.002134
	ZPRCENTRY421		zfft_r4dwpn,432K,_ac,2304,675584
	ZPRCENTRY421		zfft_r4dwpn,432K,_ac,384,310016,,,SKX
	DD			0
	PRCSTRT	8909000, 458752, 0.002446
	ZPRCENTRY421		zfft_r4dwpn,448K,_ac,3584,809728
	ZPRCENTRY421		zfft_r4dwpn,448K,_ac,512,355072
	ZPRCENTRY421		zfft_r4dwpn,448K,_ac,448,375552,,,SKX
	DD			0
	PRCSTRT	9543000, 491520, 0.002815
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,3840,867072
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,2560,746240
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,768,404224
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,640,426752,,SKX
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,512,324352
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,384,387840
	ZPRCENTRY421		zfft_r4dwpn,480K,_ac,320,381696
	DD			0
	PRCSTRT	10020000, 516096, 0.002712
	ZPRCENTRY421		zfft_r4dwpn,504K,_ac,2688,783104
	ZPRCENTRY421		zfft_r4dwpn,504K,_ac,448,354560,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,504K,_ac,384,319232
	DD			0
	PRCSTRT	10140000, 524288, 0.002819
	ZPRCENTRY421		zfft_r4dwpn,512K,_ac,4096,928512
	ZPRCENTRY421		zfft_r4dwpn,512K,_ac,512,424704,,SKX
	DD			0
	PRCSTRT	11080000, 573440, 0.003284
	ZPRCENTRY421		zfft_r4dwpn,560K,_ac,4480,1010432
	ZPRCENTRY421		zfft_r4dwpn,560K,_ac,640,435968
	ZPRCENTRY421		zfft_r4dwpn,560K,_ac,448,444672,,,SKX
	DD			0
	PRCSTRT	11420000, 589824, 0.003433
	ZPRCENTRY421		zfft_r4dwpn,576K,_ac,4608,1042176
	ZPRCENTRY421		zfft_r4dwpn,576K,_ac,3072,897792
	ZPRCENTRY421		zfft_r4dwpn,576K,_ac,768,506624
	ZPRCENTRY421		zfft_r4dwpn,576K,_ac,512,400128
	ZPRCENTRY421		zfft_r4dwpn,576K,_ac,384,448256,,,SKX
	DD			0
	PRCSTRT	11640000, 602112, 0.003313
	ZPRCENTRY421		zfft_r4dwpn,588K,_ac,3136,911616
	ZPRCENTRY421		zfft_r4dwpn,588K,_ac,448,364288,,,SKX
	DD			0
	PRCSTRT	11880000, 614400, 0.003381
	ZPRCENTRY421		zfft_r4dwpn,600K,_ac,3200,930560
	ZPRCENTRY421		zfft_r4dwpn,600K,_ac,640,397056,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,600K,_ac,320,378624
	DD			0
	PRCSTRT	12620000, 655360, 0.003735
	ZPRCENTRY421		zfft_r4dwpn,640K,_ac,5120,1153792
	ZPRCENTRY421		zfft_r4dwpn,640K,_ac,1024,531200,,SKX
	ZPRCENTRY421		zfft_r4dwpn,640K,_ac,640,521984
	ZPRCENTRY421		zfft_r4dwpn,640K,_ac,512,502528
	ZPRCENTRY421		zfft_r4dwpn,640K,_ac,320,485120
	DD			0
	PRCSTRT	13270000, 688128, 0.004010
	ZPRCENTRY421		zfft_r4dwpn,672K,_ac,5376,1211136
	ZPRCENTRY421		zfft_r4dwpn,672K,_ac,3584,1041152
	ZPRCENTRY421		zfft_r4dwpn,672K,_ac,768,516864
	ZPRCENTRY421		zfft_r4dwpn,672K,_ac,512,410368
	ZPRCENTRY421		zfft_r4dwpn,672K,_ac,448,513792,,,SKX
	DD			0
	PRCSTRT	14210000, 737280, 0.004431
	ZPRCENTRY421		zfft_r4dwpn,720K,_ac,3840,1114880
	ZPRCENTRY421		zfft_r4dwpn,720K,_ac,768,469760
	ZPRCENTRY421		zfft_r4dwpn,720K,_ac,640,490240,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,720K,_ac,384,442112
	ZPRCENTRY421		zfft_r4dwpn,720K,_ac,320,434432
	DD			0
	PRCSTRT	15120000, 786432, 0.004558
	ZPRCENTRY421		zfft_r4dwpn,768K,_ac,6144,1387264
	ZPRCENTRY421		zfft_r4dwpn,768K,_ac,4096,1192704
	ZPRCENTRY421		zfft_r4dwpn,768K,_ac,1024,666368,,SKX
	ZPRCENTRY421		zfft_r4dwpn,768K,_ac,768,619264
	ZPRCENTRY421		zfft_r4dwpn,768K,_ac,512,580352
	ZPRCENTRY421		zfft_r4dwpn,768K,_ac,384,569088
	DD			0
	PRCSTRT	15710000, 819200, 0.004927
	ZPRCENTRY421		zfft_r4dwpn,800K,_ac,6400,1445632
	ZPRCENTRY421		zfft_r4dwpn,800K,_ac,640,617216,,,SKX
	DD			0
	PRCSTRT	16490000, 860160, 0.005143
	ZPRCENTRY421		zfft_r4dwpn,840K,_ac,4480,1299200
	ZPRCENTRY421		zfft_r4dwpn,840K,_ac,640,501504,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,840K,_ac,448,504576
	DD			0
	PRCSTRT	16980000, 884736, 0.005499
	ZPRCENTRY421		zfft_r4dwpn,864K,_ac,4608,1339136
	ZPRCENTRY421		zfft_r4dwpn,864K,_ac,768,580352
	ZPRCENTRY421		zfft_r4dwpn,864K,_ac,384,506624,,,SKX
	DD			0
	PRCSTRT	17560000, 917504, 0.005404
	ZPRCENTRY421		zfft_r4dwpn,896K,_ac,7168,1612544
	ZPRCENTRY421		zfft_r4dwpn,896K,_ac,1024,678656
	ZPRCENTRY421		zfft_r4dwpn,896K,_ac,448,652032,,,SKX
	DD			0
	PRCSTRT	18830000, 983040, 0.006015
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,7680,1733376
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,5120,1483520
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,1024,615168,,SKX
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,768,731904
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,640,712448
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,512,568064
	ZPRCENTRY421		zfft_r4dwpn,960K,_ac,320,548608
	DD			0
	PRCSTRT	19530000, 1024000, 0.006539
	ZPRCENTRY421		zfft_r4dwpn,1000K,_ac,1600,816384,,SKX
	DD			0
	PRCSTRT	19770000, 1032192, 0.006539
	ZPRCENTRY421		zfft_r4dwpn,1008K,_ac,5376,1557248
	ZPRCENTRY421		zfft_r4dwpn,1008K,_ac,768,592640
	ZPRCENTRY421		zfft_r4dwpn,1008K,_ac,448,577792,,,SKX
	DD			0
	PRCSTRT	20000000, 1048576, 0.006127
	ZPRCENTRY421		zfft_r4dwpn,1M,_ac,8192,1846016
	ZPRCENTRY421		zfft_r4dwpn,1M,_ac,1024,813824,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1M,_ac,512,736000
	DD			0
	PRCSTRT	22550000, 1179648, 0.007165
	ZPRCENTRY421		zfft_r4dwpn,1152K,_ac,9216,2077440
	ZPRCENTRY421		zfft_r4dwpn,1152K,_ac,6144,1782528
	ZPRCENTRY421		zfft_r4dwpn,1152K,_ac,1024,760576,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1152K,_ac,768,844544
	ZPRCENTRY421		zfft_r4dwpn,1152K,_ac,512,649984
	ZPRCENTRY421		zfft_r4dwpn,1152K,_ac,384,638720
	DD			0
	PRCSTRT	23420000, 1228800, 0.007757
	ZPRCENTRY421		zfft_r4dwpn,1200K,_ac,6400,1857280
	ZPRCENTRY421		zfft_r4dwpn,1200K,_ac,1920,975616,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1200K,_ac,1600,1025280
	ZPRCENTRY421		zfft_r4dwpn,1200K,_ac,640,694016
	DD			0
	PRCSTRT	24920000, 1310720, 0.007975
	ZPRCENTRY421		zfft_r4dwpn,1280K,_ac,10240,2308864
	ZPRCENTRY421		zfft_r4dwpn,1280K,_ac,1024,961280
	ZPRCENTRY421		zfft_r4dwpn,1280K,_ac,640,902912,,,SKX
	DD			0
	PRCSTRT	26190000, 1376256, 0.008592
	ZPRCENTRY421		zfft_r4dwpn,1344K,_ac,7168,2073344,SKX
	ZPRCENTRY421		zfft_r4dwpn,1344K,_ac,1024,774912
	ZPRCENTRY421		zfft_r4dwpn,1344K,_ac,448,727808
	DD			0
	PRCSTRT	27180000, 1433600, 0.009242
	ZPRCENTRY421		zfft_r4dwpn,1400K,_ac,2240,1133824,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1400K,_ac,1600,1042176
	DD			0
	PRCSTRT	28010000, 1474560, 0.009692
	ZPRCENTRY421		zfft_r4dwpn,1440K,_ac,7680,2226944
	ZPRCENTRY421		zfft_r4dwpn,1440K,_ac,2304,1169152
	ZPRCENTRY421		zfft_r4dwpn,1440K,_ac,1920,1225472,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1440K,_ac,768,819968
	ZPRCENTRY421		zfft_r4dwpn,1440K,_ac,640,793344
	DD			0
	PRCSTRT	29140000, 1536000, 0.009700
	ZPRCENTRY421		zfft_r4dwpn,1500K,_ac,1600,941824,,SKX
	DD			0
	PRCSTRT	29840000, 1572864, 0.009770
	ZPRCENTRY421		zfft_r4dwpn,1536K,_ac,12288,2767616
	ZPRCENTRY421		zfft_r4dwpn,1536K,_ac,8192,2372352,SKX
	ZPRCENTRY421		zfft_r4dwpn,1536K,_ac,1024,1108736
	ZPRCENTRY421		zfft_r4dwpn,1536K,_ac,768,1069824
	ZPRCENTRY421		zfft_r4dwpn,1536K,_ac,512,817920
	DD			0
	PRCSTRT	30960000, 1638400, 0.010447
	ZPRCENTRY421		zfft_r4dwpn,1600K,_ac,12800,2874112
	ZPRCENTRY421		zfft_r4dwpn,1600K,_ac,2560,1293056
	ZPRCENTRY421		zfft_r4dwpn,1600K,_ac,1600,1251072,,SKX
	DD			0
	PRCSTRT	32520000, 1720320, 0.011459
	ZPRCENTRY421		zfft_r4dwpn,1680K,_ac,2688,1356544
	ZPRCENTRY421		zfft_r4dwpn,1680K,_ac,2240,1424640
	ZPRCENTRY421		zfft_r4dwpn,1680K,_ac,1920,1244928,,SKX
	DD			0
	PRCSTRT	33480000, 1769472, 0.011461
	ZPRCENTRY421		zfft_r4dwpn,1728K,_ac,9216,2669312
	ZPRCENTRY421		zfft_r4dwpn,1728K,_ac,2304,1468160,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1728K,_ac,768,936704
	DD			0
	PRCSTRT	34870000, 1843200, 0.012
	ZPRCENTRY421		zfft_r4dwpn,1800K,_ac,1920,1124096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,1800K,_ac,1600,1165568
	DD			0
	PRCSTRT	37110000, 1966080, 0.012478
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,15360,3451648,SKX
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,10240,2966272
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,3072,1551104
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,2560,1624832
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,1920,1494784
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,1024,1071872
	ZPRCENTRY421		zfft_r4dwpn,1920K,_ac,640,997120
	DD			0
	PRCSTRT	37740000, 2007040, 0.012046
	ZPRCENTRY421		zfft_r4dwpn,1960K,_ac,3136,1578240
	ZPRCENTRY421		zfft_r4dwpn,1960K,_ac,2240,1446656,,SKX
	DD			0
	PRCSTRT	38530000, 2048000, 0.013344
	ZPRCENTRY421		zfft_r4dwpn,2000K,_ac,3200,1610496
	ZPRCENTRY421		zfft_r4dwpn,2000K,_ac,1600,1476864,,,SKX
	DD			0
	PRCSTRT	38910000, 2064384, 0.013399
	ZPRCENTRY421		zfft_r4dwpn,2016K,_ac,2688,1704704		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2016K,_ac,2304,1490688
	DD			0
	PRCSTRT	39490000, 2097152, 0.012734
	ZPRCENTRY421		zfft_r4dwpn,2M,_ac,16384,3685120,SKX
	ZPRCENTRY421		zfft_r4dwpn,2M,_ac,1024,1403648
	DD			0
	PRCSTRT	40460000, 2150400, 0.0135
	ZPRCENTRY421		zfft_r4dwpn,2100K,_ac,2240,1305344
	ZPRCENTRY421		zfft_r4dwpn,2100K,_ac,1600,1184512,,,SKX
	DD			0
	PRCSTRT	41710000, 2211840, 0.014
	ZPRCENTRY421		zfft_r4dwpn,2160K,_ac,2304,1345280
	ZPRCENTRY421		zfft_r4dwpn,2160K,_ac,1920,1391360,,SKX
	DD			0
	PRCSTRT	43030000, 2293760, 0.014710
	ZPRCENTRY421		zfft_r4dwpn,2240K,_ac,17920,4020992,SKX
	ZPRCENTRY421		zfft_r4dwpn,2240K,_ac,3584,1800960
	ZPRCENTRY421		zfft_r4dwpn,2240K,_ac,2560,1649408
	ZPRCENTRY421		zfft_r4dwpn,2240K,_ac,2240,1737472
	DD			0
	PRCSTRT	44350000, 2359296, 0.015112
	ZPRCENTRY421		zfft_r4dwpn,2304K,_ac,18432,4164352
	ZPRCENTRY421		zfft_r4dwpn,2304K,_ac,12288,3556096
	ZPRCENTRY421		zfft_r4dwpn,2304K,_ac,3072,1948416
	ZPRCENTRY421		zfft_r4dwpn,2304K,_ac,2304,1789696,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2304K,_ac,1024,1223424
	ZPRCENTRY421		zfft_r4dwpn,2304K,_ac,768,1176320
	DD			0
	PRCSTRT	45130000, 2408448, 0.015855
	ZPRCENTRY421		zfft_r4dwpn,2352K,_ac,3136,1983744
	ZPRCENTRY421		zfft_r4dwpn,2352K,_ac,2688,1730304		;,,,SKX
	DD			0
	PRCSTRT	46120000, 2457600, 0.016488
	ZPRCENTRY421		zfft_r4dwpn,2400K,_ac,12800,3695360
	ZPRCENTRY421		zfft_r4dwpn,2400K,_ac,3840,1927936
	ZPRCENTRY421		zfft_r4dwpn,2400K,_ac,3200,2024192
	ZPRCENTRY421		zfft_r4dwpn,2400K,_ac,2560,1487616
	ZPRCENTRY421		zfft_r4dwpn,2400K,_ac,1920,1764096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2400K,_ac,1600,1702656
	DD			0
	PRCSTRT	48450000, 2580480, 0.0165
	ZPRCENTRY421		zfft_r4dwpn,2520K,_ac,2688,1560320
	ZPRCENTRY421		zfft_r4dwpn,2520K,_ac,2240,1616128
	ZPRCENTRY421		zfft_r4dwpn,2520K,_ac,1920,1412864,,,SKX
	DD			0
	PRCSTRT	48990000, 2621440, 0.016596
	ZPRCENTRY421		zfft_r4dwpn,2560K,_ac,20480,4598528
	ZPRCENTRY421		zfft_r4dwpn,2560K,_ac,4096,2059008
	ZPRCENTRY421		zfft_r4dwpn,2560K,_ac,2560,1981184,,SKX
	DD			0
	PRCSTRT	49790000, 2654208, 0.0165
	ZPRCENTRY421		zfft_r4dwpn,2592K,_ac,2304,1664768,,SKX
	DD			0
	PRCSTRT	51460000, 2752512, 0.017883
	ZPRCENTRY421		zfft_r4dwpn,2688K,_ac,21504,4827904
	ZPRCENTRY421		zfft_r4dwpn,2688K,_ac,3584,2263808,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2688K,_ac,3072,1977088
	ZPRCENTRY421		zfft_r4dwpn,2688K,_ac,2688,2078464
	DD			0
	PRCSTRT	52390000, 2809856, 0.017883
	ZPRCENTRY421		zfft_r4dwpn,2744K,_ac,3136,2012928		;,,,SKX
	DD			0
	PRCSTRT	53510000, 2867200, 0.019289
	ZPRCENTRY421		zfft_r4dwpn,2800K,_ac,4480,2245376
	ZPRCENTRY421		zfft_r4dwpn,2800K,_ac,3200,2053888		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,2800K,_ac,2240,2050304
	DD			0
	PRCSTRT	55260000, 2949120, 0.020031
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,15360,4436736,SKX
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,4608,2311936
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,3840,2423552
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,3072,1782528
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,2560,1841920
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,2304,2111232
	ZPRCENTRY421		zfft_r4dwpn,2880K,_ac,1920,2033408
	DD			0
	PRCSTRT	56240000, 3010560, 0.02007
	ZPRCENTRY421		zfft_r4dwpn,2940K,_ac,3136,1814272
	ZPRCENTRY421		zfft_r4dwpn,2940K,_ac,2240,1640192,,SKX
	DD			0
	PRCSTRT	57430000, 3072000, 0.02008
	ZPRCENTRY421		zfft_r4dwpn,3000K,_ac,3200,1851136
	ZPRCENTRY421		zfft_r4dwpn,3000K,_ac,1600,1638144,,,SKX
	DD			0
	PRCSTRT	57960000, 3096576, 0.0201
	ZPRCENTRY421		zfft_r4dwpn,3024K,_ac,2688,1932032
	ZPRCENTRY421		zfft_r4dwpn,3024K,_ac,2304,1689344		;,,,SKX
	DD			0
	PRCSTRT	58830000, 3145728, 0.020198
	ZPRCENTRY421		zfft_r4dwpn,3M,_ac,24576,5548800
	ZPRCENTRY421		zfft_r4dwpn,3M,_ac,16384,4735744,SKX
	ZPRCENTRY421		zfft_r4dwpn,3M,_ac,4096,2587392
	ZPRCENTRY421		zfft_r4dwpn,3M,_ac,3072,2374400
	ZPRCENTRY421		zfft_r4dwpn,3M,_ac,1024,1534720
	DD			0
	PRCSTRT	59720000, 3211264, 0.024079
	ZPRCENTRY421		zfft_r4dwpn,3136K,_ac,25088,5626624,SKX
	ZPRCENTRY421		zfft_r4dwpn,3136K,_ac,3584,2296576
	ZPRCENTRY421		zfft_r4dwpn,3136K,_ac,3136,2418432
	DD			0
	PRCSTRT	61010000, 3276800, 0.021621
	ZPRCENTRY421		zfft_r4dwpn,3200K,_ac,5120,2562816
	ZPRCENTRY421		zfft_r4dwpn,3200K,_ac,3200,2467584
	ZPRCENTRY421		zfft_r4dwpn,3200K,_ac,2560,2337536
	ZPRCENTRY421		zfft_r4dwpn,3200K,_ac,1600,2154240,,,SKX
	DD			0
	PRCSTRT	64140000, 3440640, 0.023358
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,17920,5169920,SKX
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,5376,2689792
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,4480,2822912
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,3840,2458368
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,3584,2069248
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,2688,2452224
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,2560,1868544
	ZPRCENTRY421		zfft_r4dwpn,3360K,_ac,2240,2363136
	DD			0
	PRCSTRT	66070000, 3538944, 0.023780
	ZPRCENTRY421		zfft_r4dwpn,3456K,_ac,18432,5346048,SKX
	ZPRCENTRY421		zfft_r4dwpn,3456K,_ac,4608,2905856
	ZPRCENTRY421		zfft_r4dwpn,3456K,_ac,3072,2206464
	ZPRCENTRY421		zfft_r4dwpn,3456K,_ac,2304,2432768
	DD			0
	PRCSTRT	67340000, 3612672, 0.024
	ZPRCENTRY421		zfft_r4dwpn,3528K,_ac,3136,2246912
	ZPRCENTRY421		zfft_r4dwpn,3528K,_ac,2688,1959680		;,,,SKX
	DD			0
	PRCSTRT	68060000, 3670016, 0.024079
	ZPRCENTRY421		zfft_r4dwpn,3584K,_ac,28672,6433536
	ZPRCENTRY421		zfft_r4dwpn,3584K,_ac,4096,2624256		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,3584K,_ac,3584,2759424
	DD			0
	PRCSTRT	68740000, 3686400, 0.025
	ZPRCENTRY421		zfft_r4dwpn,3600K,_ac,3840,2214656
	ZPRCENTRY421		zfft_r4dwpn,3600K,_ac,3200,2292480
	ZPRCENTRY421		zfft_r4dwpn,3600K,_ac,1920,1953536
	ZPRCENTRY421		zfft_r4dwpn,3600K,_ac,1600,1868032,,,SKX
	DD			0
	PRCSTRT	73140000, 3932160, 0.026437
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,30720,6892288
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,20480,5911296
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,6144,3074816
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,5120,3222272
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,4096,2364160
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,3840,2953984
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,3072,2800384
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,2560,2693888
	ZPRCENTRY421		zfft_r4dwpn,3840K,_ac,1920,2572032,,,SKX
	DD			0
	PRCSTRT	74430000, 4014080, 0.02605
	ZPRCENTRY421		zfft_r4dwpn,3920K,_ac,4480,2862848
	ZPRCENTRY421		zfft_r4dwpn,3920K,_ac,3136,2853120,,SKX
	DD			0
	PRCSTRT	76010000, 4096000, 0.029047
	ZPRCENTRY421		zfft_r4dwpn,4000K,_ac,6400,3202816
	ZPRCENTRY421		zfft_r4dwpn,4000K,_ac,3200,2910976		;,,,SKX
	DD			0
	PRCSTRT	76780000, 4128768, 0.027915
	ZPRCENTRY421		zfft_r4dwpn,4032K,_ac,21504,6206208
	ZPRCENTRY421		zfft_r4dwpn,4032K,_ac,5376,3382016
	ZPRCENTRY421		zfft_r4dwpn,4032K,_ac,4608,2946816
	ZPRCENTRY421		zfft_r4dwpn,4032K,_ac,3584,2562816,,SKX
	ZPRCENTRY421		zfft_r4dwpn,4032K,_ac,3072,2237184
	ZPRCENTRY421		zfft_r4dwpn,4032K,_ac,2688,2825984
	DD			0
	PRCSTRT	77630000, 4194304, 0.027172
	ZPRCENTRY421		zfft_r4dwpn,4M,_ac,32768,7383808
	ZPRCENTRY421		zfft_r4dwpn,4M,_ac,4096,3152640			;,,,SKX
	DD			0
	PRCSTRT	78180000, 4214784, 0.0275
	ZPRCENTRY421		zfft_r4dwpn,4116K,_ac,3136,2278144		;,,,SKX
	DD			0
	PRCSTRT	79810000, 4300800, 0.028
	ZPRCENTRY421		zfft_r4dwpn,4200K,_ac,4480,2578176
	ZPRCENTRY421		zfft_r4dwpn,4200K,_ac,3200,2324224
	ZPRCENTRY421		zfft_r4dwpn,4200K,_ac,2240,2267904,,,SKX
	DD			0
	PRCSTRT	82260000, 4423680, 0.029
	ZPRCENTRY421		zfft_r4dwpn,4320K,_ac,4608,2653952
	ZPRCENTRY421		zfft_r4dwpn,4320K,_ac,3840,2743040
	ZPRCENTRY421		zfft_r4dwpn,4320K,_ac,2304,2334464
	ZPRCENTRY421		zfft_r4dwpn,4320K,_ac,1920,2226944,,,SKX
	DD			0
	PRCSTRT	84940000, 4587520, 0.031930
	ZPRCENTRY421		zfft_r4dwpn,4480K,_ac,7168,3578624
	ZPRCENTRY421		zfft_r4dwpn,4480K,_ac,5120,3267328
	ZPRCENTRY421		zfft_r4dwpn,4480K,_ac,4480,3440384
	ZPRCENTRY421		zfft_r4dwpn,4480K,_ac,3584,3255040
	ZPRCENTRY421		zfft_r4dwpn,4480K,_ac,2240,2988800,,,SKX
	DD			0
	PRCSTRT	87470000, 4718592, 0.031494
	ZPRCENTRY421		zfft_r4dwpn,4608K,_ac,24576,7123712,SKX
	ZPRCENTRY421		zfft_r4dwpn,4608K,_ac,6144,3865344
	ZPRCENTRY421		zfft_r4dwpn,4608K,_ac,4608,3540736
	ZPRCENTRY421		zfft_r4dwpn,4608K,_ac,4096,2927360
	ZPRCENTRY421		zfft_r4dwpn,4608K,_ac,3072,3226368
	ZPRCENTRY421		zfft_r4dwpn,4608K,_ac,2304,3075840
	DD			0
	PRCSTRT	89110000, 4816896, 0.033388
	ZPRCENTRY421		zfft_r4dwpn,4704K,_ac,25088,7234304,SKX
	ZPRCENTRY421		zfft_r4dwpn,4704K,_ac,5376,3429120
	ZPRCENTRY421		zfft_r4dwpn,4704K,_ac,3584,2597632
	ZPRCENTRY421		zfft_r4dwpn,4704K,_ac,3136,3287808
	DD			0
	PRCSTRT	91020000, 4915200, 0.034693
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,7680,3838720
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,6400,4026112
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,5120,2941696
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,3840,3484416
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,3200,3354368
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,2560,2583296
	ZPRCENTRY421		zfft_r4dwpn,4800K,_ac,1600,2340608,,,SKX
	DD			0
	PRCSTRT	95590000, 5160960, 0.035
	ZPRCENTRY421		zfft_r4dwpn,5040K,_ac,5376,3087104,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5040K,_ac,4480,3193600
	ZPRCENTRY421		zfft_r4dwpn,5040K,_ac,3840,2779904
	ZPRCENTRY421		zfft_r4dwpn,5040K,_ac,2688,2709248
	ZPRCENTRY421		zfft_r4dwpn,5040K,_ac,2240,2584832
	DD			0
	PRCSTRT	96750000, 5242880, 0.035846
	ZPRCENTRY421		zfft_r4dwpn,5M,_ac,8192,4090624
	ZPRCENTRY421		zfft_r4dwpn,5M,_ac,5120,3926784,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5M,_ac,4096,3717888
	ZPRCENTRY421		zfft_r4dwpn,5M,_ac,2560,3406592
	DD			0
	PRCSTRT	98220000, 5308416, 0.036
	ZPRCENTRY421		zfft_r4dwpn,5184K,_ac,4608,3286784
	ZPRCENTRY421		zfft_r4dwpn,5184K,_ac,2304,2660096,,,SKX
	DD			0
	PRCSTRT	101500000, 5505024, 0.038004
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,28672,8270592
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,7168,4500224
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,6144,3918592
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,5376,4121344,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,4096,2966272
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,3584,3750656
	ZPRCENTRY421		zfft_r4dwpn,5376K,_ac,2688,3573504
	DD			0
	PRCSTRT	105500000, 5734400, 0.040691
	ZPRCENTRY421		zfft_r4dwpn,5600K,_ac,6400,4081408
	ZPRCENTRY421		zfft_r4dwpn,5600K,_ac,4480,4057856		;,,,SKX
	DD			0
	PRCSTRT	108900000, 5898240, 0.041629
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,30720,8860416
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,9216,4600576
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,7680,4825856
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,6144,3527424
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,5120,3644160
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,4608,4175616
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,3840,4014848
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,3072,3091200
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,2560,2943744
	ZPRCENTRY421		zfft_r4dwpn,5760K,_ac,1920,2789120,,,SKX
	DD			0
	PRCSTRT	110700000, 6021120, 0.0418
	ZPRCENTRY421		zfft_r4dwpn,5880K,_ac,4480,3235584		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,5880K,_ac,3136,3149568
	DD			0
	PRCSTRT	113200000, 6144000, 0.042
	ZPRCENTRY421		zfft_r4dwpn,6000K,_ac,6400,3673856		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6000K,_ac,3200,3213056
	DD			0
	PRCSTRT	114200000, 6193152, 0.043
	ZPRCENTRY421		zfft_r4dwpn,6048K,_ac,5376,3824384,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6048K,_ac,4608,3329792
	ZPRCENTRY421		zfft_r4dwpn,6048K,_ac,2688,3087104
	DD			0
	PRCSTRT	115700000, 6291456, 0.043121
	ZPRCENTRY421		zfft_r4dwpn,6M,_ac,32768,9483008
	ZPRCENTRY421		zfft_r4dwpn,6M,_ac,8192,5143296
	ZPRCENTRY421		zfft_r4dwpn,6M,_ac,6144,4709120,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6M,_ac,4096,4283136
	ZPRCENTRY421		zfft_r4dwpn,6M,_ac,3072,4078336
	DD			0
	PRCSTRT	117700000, 6422528, 0.042000
	ZPRCENTRY421		zfft_r4dwpn,6272K,_ac,7168,4561664,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6272K,_ac,3136,4157184
	DD			0
	PRCSTRT	120300000, 6553600, 0.045826
	ZPRCENTRY421		zfft_r4dwpn,6400K,_ac,10240,5110528,SKX
	ZPRCENTRY421		zfft_r4dwpn,6400K,_ac,6400,4904704
	ZPRCENTRY421		zfft_r4dwpn,6400K,_ac,5120,4631296
	ZPRCENTRY421		zfft_r4dwpn,6400K,_ac,3200,4241152
	DD			0
	PRCSTRT	126300000, 6881280, 0.048606
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,7680,4891392
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,7168,4104960
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,5376,4860672,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,5120,3691264
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,4480,4675328
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,3584,3590912
	ZPRCENTRY421		zfft_r4dwpn,6720K,_ac,2240,3236608
	DD			0
	PRCSTRT	130100000, 7077888, 0.051987
	ZPRCENTRY421		zfft_r4dwpn,6912K,_ac,9216,5784320
	ZPRCENTRY421		zfft_r4dwpn,6912K,_ac,6144,4369152
	ZPRCENTRY421		zfft_r4dwpn,6912K,_ac,4608,4810496
	ZPRCENTRY421		zfft_r4dwpn,6912K,_ac,3072,3521280		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,6912K,_ac,2304,3329792
	DD			0
	PRCSTRT	132700000, 7225344, 0.052
	ZPRCENTRY421		zfft_r4dwpn,7056K,_ac,5376,3873536,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7056K,_ac,3136,3588352
	DD			0
	PRCSTRT	134200000, 7340032, 0.052097
	ZPRCENTRY421		zfft_r4dwpn,7M,_ac,8192,5212928
	ZPRCENTRY421		zfft_r4dwpn,7M,_ac,7168,5483264,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7M,_ac,3584,4741888
	DD			0
	PRCSTRT	135400000, 7372800, 0.053
	ZPRCENTRY421		zfft_r4dwpn,7200K,_ac,7680,4401920,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7200K,_ac,6400,4550400
	ZPRCENTRY421		zfft_r4dwpn,7200K,_ac,3840,3842816
	ZPRCENTRY421		zfft_r4dwpn,7200K,_ac,3200,3660544
	DD			0
	PRCSTRT	143900000, 7864320, 0.054932
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,12288,6126336
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,10240,6425344,,SKX
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,8192,4690688
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,7680,5878528
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,6144,5552896
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,5120,5335808
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,4096,4098816
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,3840,5075712
	ZPRCENTRY421		zfft_r4dwpn,7680K,_ac,2560,3685120
	DD			0
	PRCSTRT	149600000, 8192000, 0.062103
	ZPRCENTRY421		zfft_r4dwpn,8000K,_ac,12800,6372096
	ZPRCENTRY421		zfft_r4dwpn,8000K,_ac,6400,5783296		;,,,SKX
	DD			0
	PRCSTRT	151200000, 8257536, 0.063072
	ZPRCENTRY421		zfft_r4dwpn,8064K,_ac,9216,5862144
	ZPRCENTRY421		zfft_r4dwpn,8064K,_ac,7168,5085952
	ZPRCENTRY421		zfft_r4dwpn,8064K,_ac,6144,4424448
	ZPRCENTRY421		zfft_r4dwpn,8064K,_ac,5376,5600000,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8064K,_ac,3584,4090624
	ZPRCENTRY421		zfft_r4dwpn,8064K,_ac,2688,3864320
	DD			0
	PRCSTRT	152900000, 8388608, 0.059916
	ZPRCENTRY421		zfft_r4dwpn,8M,_ac,8192,6265600
	ZPRCENTRY421		zfft_r4dwpn,8M,_ac,4096,5413632			;,,,SKX
	DD			0
	PRCSTRT	157100000, 8601600, 0.063
	ZPRCENTRY421		zfft_r4dwpn,8400K,_ac,6400,4607744,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8400K,_ac,4480,4472576
	DD			0
	PRCSTRT	161900000, 8847360, 0.064
	ZPRCENTRY421		zfft_r4dwpn,8640K,_ac,9216,5274368
	ZPRCENTRY421		zfft_r4dwpn,8640K,_ac,7680,5452544,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8640K,_ac,4608,4601600
	ZPRCENTRY421		zfft_r4dwpn,8640K,_ac,3840,4377344
	DD			0
	PRCSTRT	167000000, 9175040, 0.068512
	ZPRCENTRY421		zfft_r4dwpn,8960K,_ac,10240,6511360,,SKX
	ZPRCENTRY421		zfft_r4dwpn,8960K,_ac,7168,6466304
	ZPRCENTRY421		zfft_r4dwpn,8960K,_ac,4480,5910272
	DD			0
	PRCSTRT	172200000, 9437184, 0.069579
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,12288,7703296
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,9216,7045888
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,8192,5810944
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,6144,6396672
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,4608,6080256
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,4096,4668160			;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,9M,_ac,3072,4406016
	DD			0
	PRCSTRT	175300000, 9633792, 0.071508
	ZPRCENTRY421		zfft_r4dwpn,9408K,_ac,7168,5149440		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,9408K,_ac,3136,4491008
	DD			0
	PRCSTRT	178900000, 9830400, 0.072967
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,15360,7645952,SKX
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,12800,8014592
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,10240,5858048
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,7680,6931200
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,6400,6661888
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,5120,5102336
	ZPRCENTRY421		zfft_r4dwpn,9600K,_ac,3200,4581120
	DD			0
	PRCSTRT	188000000, 10321920, 0.077
	ZPRCENTRY421		zfft_r4dwpn,10080K,_ac,7680,5520128
	ZPRCENTRY421		zfft_r4dwpn,10080K,_ac,5376,5354240		;,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,10080K,_ac,4480,5094144
	DD			0
	PRCSTRT	190300000, 10485760, 0.078834
	ZPRCENTRY421		zfft_r4dwpn,10M,_ac,16384,8157952,SKX
	ZPRCENTRY421		zfft_r4dwpn,10M,_ac,10240,7826176
	ZPRCENTRY421		zfft_r4dwpn,10M,_ac,8192,7387904
	ZPRCENTRY421		zfft_r4dwpn,10M,_ac,5120,6744832
	DD			0
	PRCSTRT	193300000, 10616832, 0.079
	ZPRCENTRY421		zfft_r4dwpn,10368K,_ac,9216,6533888,,SKX
	ZPRCENTRY421		zfft_r4dwpn,10368K,_ac,4608,5240576
	DD			0
	PRCSTRT	200000000, 11010048, 0.083177
	ZPRCENTRY421		zfft_r4dwpn,10752K,_ac,12288,7805696
	ZPRCENTRY421		zfft_r4dwpn,10752K,_ac,8192,5882624
	ZPRCENTRY421		zfft_r4dwpn,10752K,_ac,7168,7449344,,SKX
	ZPRCENTRY421		zfft_r4dwpn,10752K,_ac,5376,7078656
	ZPRCENTRY421		zfft_r4dwpn,10752K,_ac,3584,5118720
	DD			0
	PRCSTRT	207400000, 11468800, 0.085770
	ZPRCENTRY421		zfft_r4dwpn,11200K,_ac,17920,8911616,SKX
	ZPRCENTRY421		zfft_r4dwpn,11200K,_ac,12800,8121088
	DD			0
	PRCSTRT	214100000, 11796480, 0.088274
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,18432,9194240
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,15360,9616128,,SKX
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,12288,7021312
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,10240,7256832
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,9216,8307456
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,7680,7983872
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,6144,6114048
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,5120,5810944
	ZPRCENTRY421		zfft_r4dwpn,11520K,_ac,3840,5477120
	DD			0
	PRCSTRT	223100000, 12288000, 0.094
	ZPRCENTRY421		zfft_r4dwpn,12000K,_ac,12800,7303936
	ZPRCENTRY421		zfft_r4dwpn,12000K,_ac,6400,6366976		;,,,SKX
	DD			0
	PRCSTRT	225100000, 12386304, 0.095
	ZPRCENTRY421		zfft_r4dwpn,12096K,_ac,9216,6613760
	ZPRCENTRY421		zfft_r4dwpn,12096K,_ac,5376,6097664		;,,,SKX
	DD			0
	PRCSTRT	228000000, 12582912, 0.096630
	ZPRCENTRY421		zfft_r4dwpn,12M,_ac,16384,10259200,,SKX
	ZPRCENTRY421		zfft_r4dwpn,12M,_ac,12288,9382656
	ZPRCENTRY421		zfft_r4dwpn,12M,_ac,8192,8510208
	ZPRCENTRY421		zfft_r4dwpn,12M,_ac,6144,8084224
	ZPRCENTRY421		zfft_r4dwpn,12M,_ac,4096,5839616
	DD			0
	PRCSTRT	237000000, 13107200, 0.097776
	ZPRCENTRY421		zfft_r4dwpn,12800K,_ac,20480,10185472
	ZPRCENTRY421		zfft_r4dwpn,12800K,_ac,12800,9763584
	ZPRCENTRY421		zfft_r4dwpn,12800K,_ac,10240,9227008,,SKX
	ZPRCENTRY421		zfft_r4dwpn,12800K,_ac,6400,8419072
	DD			0
	PRCSTRT	248500000, 13762560, 0.106540
	ZPRCENTRY421		zfft_r4dwpn,13440K,_ac,21504,10693376
	ZPRCENTRY421		zfft_r4dwpn,13440K,_ac,17920,11209472
	ZPRCENTRY421		zfft_r4dwpn,13440K,_ac,15360,9743104,,SKX
	ZPRCENTRY421		zfft_r4dwpn,13440K,_ac,10240,7344896
	ZPRCENTRY421		zfft_r4dwpn,13440K,_ac,7168,7117568
	ZPRCENTRY421		zfft_r4dwpn,13440K,_ac,4480,6373120
	DD			0
	PRCSTRT	255900000, 14155776, 0.115310
	ZPRCENTRY421		zfft_r4dwpn,13824K,_ac,18432,11557632,SKX
	ZPRCENTRY421		zfft_r4dwpn,13824K,_ac,12288,8698624
	ZPRCENTRY421		zfft_r4dwpn,13824K,_ac,9216,9569024
	ZPRCENTRY421		zfft_r4dwpn,13824K,_ac,6144,6961920
	ZPRCENTRY421		zfft_r4dwpn,13824K,_ac,4608,6555392
	DD			0
	PRCSTRT	264500000, 14680064, 0.114093
	ZPRCENTRY421		zfft_r4dwpn,14M,_ac,16384,10394368,,SKX
	ZPRCENTRY421		zfft_r4dwpn,14M,_ac,7168,9415424
	DD			0
	PRCSTRT	266700000, 14745600, 0.118
	ZPRCENTRY421		zfft_r4dwpn,14400K,_ac,15360,8762112,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,14400K,_ac,12800,9050880
	ZPRCENTRY421		zfft_r4dwpn,14400K,_ac,7680,7627520
	ZPRCENTRY421		zfft_r4dwpn,14400K,_ac,6400,7249664
	DD			0
	PRCSTRT	283400000, 15728640, 0.128292
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,24576,12249856
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,20480,12811008
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,16384,9347840
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,15360,11713280,,SKX
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,12288,11062016
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,10240,10627840
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,8192,8129280
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,7680,10089216
	ZPRCENTRY421		zfft_r4dwpn,15M,_ac,5120,7269120
	DD			0
	PRCSTRT	288300000, 16056320, 0.120880
	ZPRCENTRY421		zfft_r4dwpn,15680K,_ac,25088,12466944,SKX
	ZPRCENTRY421		zfft_r4dwpn,15680K,_ac,17920,11356928
	DD			0
	PRCSTRT	294700000, 16384000, 0.130344
	ZPRCENTRY421		zfft_r4dwpn,16000K,_ac,12800,11512576		;,,,SKX
	DD			0
	PRCSTRT	297600000, 16515072, 0.132344
	ZPRCENTRY421		zfft_r4dwpn,16128K,_ac,21504,13449984
	ZPRCENTRY421		zfft_r4dwpn,16128K,_ac,18432,11709184,,SKX
	ZPRCENTRY421		zfft_r4dwpn,16128K,_ac,12288,8803072
	ZPRCENTRY421		zfft_r4dwpn,16128K,_ac,7168,8104704
	ZPRCENTRY421		zfft_r4dwpn,16128K,_ac,5376,7627520
	DD			0
	PRCSTRT	301700000, 16777216, 0.133010
	ZPRCENTRY421		zfft_r4dwpn,16M,_ac,16384,12495616,,SKX
	ZPRCENTRY421		zfft_r4dwpn,16M,_ac,8192,10754816
	DD			0
	PRCSTRT	309200000, 17203200, 0.139
	ZPRCENTRY421		zfft_r4dwpn,16800K,_ac,17920,10212096,,SKX
	ZPRCENTRY421		zfft_r4dwpn,16800K,_ac,12800,9159424
	DD			0
	PRCSTRT	318800000, 17694720, 0.141
	ZPRCENTRY421		zfft_r4dwpn,17280K,_ac,18432,10531584
	ZPRCENTRY421		zfft_r4dwpn,17280K,_ac,15360,10857216,,SKX
	ZPRCENTRY421		zfft_r4dwpn,17280K,_ac,9216,9138944
	ZPRCENTRY421		zfft_r4dwpn,17280K,_ac,7680,8684288
	DD			0
	PRCSTRT	328500000, 18350080, 0.147854
	ZPRCENTRY421		zfft_r4dwpn,17920K,_ac,28672,14248704
	ZPRCENTRY421		zfft_r4dwpn,17920K,_ac,20480,12978944
	ZPRCENTRY421		zfft_r4dwpn,17920K,_ac,17920,13654784,,SKX
	DD			0
	PRCSTRT	339700000, 18874368, 0.153546
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,24576,15399680
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,18432,14072576
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,16384,11582208,,SKX
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,12288,12741376
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,9216,12092160
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,8192,9255680
	ZPRCENTRY421		zfft_r4dwpn,18M,_ac,6144,8706816
	DD			0
	PRCSTRT	345000000, 19267584, 0.149016
	ZPRCENTRY421		zfft_r4dwpn,18816K,_ac,25088,15682304
	ZPRCENTRY421		zfft_r4dwpn,18816K,_ac,21504,13626112		;,,,SKX
	DD			0
	PRCSTRT	352500000, 19660800, 0.156971
	ZPRCENTRY421		zfft_r4dwpn,19200K,_ac,30720,15264512
	ZPRCENTRY421		zfft_r4dwpn,19200K,_ac,20480,11670272
	ZPRCENTRY421		zfft_r4dwpn,19200K,_ac,15360,13810432,,SKX
	ZPRCENTRY421		zfft_r4dwpn,19200K,_ac,12800,13261568
	ZPRCENTRY421		zfft_r4dwpn,19200K,_ac,10240,10148608
	ZPRCENTRY421		zfft_r4dwpn,19200K,_ac,6400,9066240
	DD			0
	PRCSTRT	370300000, 20643840, 0.172
	ZPRCENTRY421		zfft_r4dwpn,20160K,_ac,21504,12251904
	ZPRCENTRY421		zfft_r4dwpn,20160K,_ac,17920,12655360
	ZPRCENTRY421		zfft_r4dwpn,20160K,_ac,15360,10986240,,SKX
	DD			0
	PRCSTRT	375200000, 20971520, 0.172761
	ZPRCENTRY421		zfft_r4dwpn,20M,_ac,32768,16313088
	ZPRCENTRY421		zfft_r4dwpn,20M,_ac,20480,15604480
	ZPRCENTRY421		zfft_r4dwpn,20M,_ac,16384,14732032,,SKX
	ZPRCENTRY421		zfft_r4dwpn,20M,_ac,10240,13429504
	DD			0
	PRCSTRT	381100000, 21233664, 0.174
	ZPRCENTRY421		zfft_r4dwpn,20736K,_ac,18432,13044480,,SKX
	ZPRCENTRY421		zfft_r4dwpn,20736K,_ac,9216,10404608
	DD			0
	PRCSTRT	394100000, 22020096, 0.179899
	ZPRCENTRY421		zfft_r4dwpn,21M,_ac,28672,17922816
	ZPRCENTRY421		zfft_r4dwpn,21M,_ac,24576,15600384
	ZPRCENTRY421		zfft_r4dwpn,21M,_ac,21504,16382720
	ZPRCENTRY421		zfft_r4dwpn,21M,_ac,16384,11719424,,SKX
	ZPRCENTRY421		zfft_r4dwpn,21M,_ac,7168,10136320
	DD			0
	PRCSTRT	400000000, 22478848, 0.194140
	ZPRCENTRY421		zfft_r4dwpn,21952K,_ac,25088,15887104,,SKX
	DD			0
	PRCSTRT	408800000, 22937600, 0.186258
	ZPRCENTRY421		zfft_r4dwpn,22400K,_ac,17920,16100096,,SKX
	DD			0
	PRCSTRT	421300000, 23592960, 0.198884
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,30720,19200768
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,24576,14029568
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,20480,14461696
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,18432,16587520
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,15360,15907584,,SKX
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,12288,12163840
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,10240,11553536
	ZPRCENTRY421		zfft_r4dwpn,23040K,_ac,7680,10859264
	DD			0
	PRCSTRT	430200000, 24084480, 0.199
	ZPRCENTRY421		zfft_r4dwpn,23520K,_ac,25088,14283520
	ZPRCENTRY421		zfft_r4dwpn,23520K,_ac,17920,12804864,,SKX
	DD			0
	PRCSTRT	438700000, 24576000, 0.200
	ZPRCENTRY421		zfft_r4dwpn,24000K,_ac,12800,12659456		;,,,SKX
	DD			0
	PRCSTRT	442800000, 24772608, 0.209253
	ZPRCENTRY421		zfft_r4dwpn,24192K,_ac,21504,15182592
	ZPRCENTRY421		zfft_r4dwpn,24192K,_ac,18432,13198080,,SKX
	DD			0
	PRCSTRT	447900000, 25165824, 0.209253
	ZPRCENTRY421		zfft_r4dwpn,24M,_ac,32768,20511488
	ZPRCENTRY421		zfft_r4dwpn,24M,_ac,24576,18750208,,SKX
	ZPRCENTRY421		zfft_r4dwpn,24M,_ac,16384,16968448
	ZPRCENTRY421		zfft_r4dwpn,24M,_ac,12288,16100096
	ZPRCENTRY421		zfft_r4dwpn,24M,_ac,8192,11574016
	DD			0
	PRCSTRT	455400000, 25690112, 0.209253
	ZPRCENTRY421		zfft_r4dwpn,25088K,_ac,28672,18156288
	ZPRCENTRY421		zfft_r4dwpn,25088K,_ac,25088,19102464,,SKX
	DD			0
	PRCSTRT	465700000, 26214400, 0.212801
	ZPRCENTRY421		zfft_r4dwpn,25M,_ac,20480,18397952,,SKX
	ZPRCENTRY421		zfft_r4dwpn,25M,_ac,12800,16759552
	DD			0
	PRCSTRT	489200000, 27525120, 0.229014
	ZPRCENTRY421		zfft_r4dwpn,26880K,_ac,30720,19450624
	ZPRCENTRY421		zfft_r4dwpn,26880K,_ac,28672,16323328
	ZPRCENTRY421		zfft_r4dwpn,26880K,_ac,21504,19315456
	ZPRCENTRY421		zfft_r4dwpn,26880K,_ac,20480,14631680
	ZPRCENTRY421		zfft_r4dwpn,26880K,_ac,17920,18545408,,SKX
	DD			0
	PRCSTRT	504800000, 28311552, 0.23830
	ZPRCENTRY421		zfft_r4dwpn,27M,_ac,24576,17378048,,SKX
	ZPRCENTRY421		zfft_r4dwpn,27M,_ac,18432,19102464
	ZPRCENTRY421		zfft_r4dwpn,27M,_ac,12288,13847296
	ZPRCENTRY421		zfft_r4dwpn,27M,_ac,9216,13009664
	DD			0
	PRCSTRT	513700000, 28901376, 0.251
	ZPRCENTRY421		zfft_r4dwpn,28224K,_ac,25088,17701632,,SKX
	ZPRCENTRY421		zfft_r4dwpn,28224K,_ac,21504,15360768
	DD			0
	PRCSTRT	519100000, 29360128, 0.247130
	ZPRCENTRY421		zfft_r4dwpn,28M,_ac,32768,20777728
	ZPRCENTRY421		zfft_r4dwpn,28M,_ac,28672,21830400		;,,,SKX
	DD			0
	PRCSTRT	524600000, 29491200, 0.251
	ZPRCENTRY421		zfft_r4dwpn,28800K,_ac,30720,17486592
	ZPRCENTRY421		zfft_r4dwpn,28800K,_ac,15360,15182592,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,28800K,_ac,12800,14412544
	DD			0
	PRCSTRT	557000000, 31457280, 0.265964
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,32768,18682624
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,30720,23386880
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,24576,22100736,,SKX
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,20480,21191424
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,16384,16194304
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,15360,20101888
	ZPRCENTRY421		zfft_r4dwpn,30M,_ac,10240,14445312
	DD			0
	PRCSTRT	566800000, 32112640, 0.306
	ZPRCENTRY421		zfft_r4dwpn,31360K,_ac,25088,22522624,,SKX
	DD			0
	PRCSTRT	585600000, 33030144, 0.2747
	ZPRCENTRY421		zfft_r4dwpn,32256K,_ac,28672,20228864
	ZPRCENTRY421		zfft_r4dwpn,32256K,_ac,24576,17580800,,SKX
	ZPRCENTRY421		zfft_r4dwpn,32256K,_ac,21504,22248192
	DD			0
	PRCSTRT	593400000, 33554432, 0.287437
	ZPRCENTRY421		zfft_r4dwpn,32M,_ac,32768,24976128
	ZPRCENTRY421		zfft_r4dwpn,32M,_ac,16384,21441280		;,,,SKX
	DD			0
	PRCSTRT	595100000, 33718272, 0.290
	ZPRCENTRY421		zfft_r4dwpn,32928K,_ac,25088,17908480,,SKX
	DD			0
	PRCSTRT	608500000, 34406400, 0.295
	ZPRCENTRY421		zfft_r4dwpn,33600K,_ac,17920,17697536		;,,,SKX
	DD			0
	PRCSTRT	627800000, 35389440, 0.300
	ZPRCENTRY421		zfft_r4dwpn,34560K,_ac,30720,21670656
	ZPRCENTRY421		zfft_r4dwpn,34560K,_ac,18432,18230016
	ZPRCENTRY421		zfft_r4dwpn,34560K,_ac,15360,17283840,,,SKX
	DD			0
	PRCSTRT	646700000, 36700160, 0.306
	ZPRCENTRY421		zfft_r4dwpn,35M,_ac,28672,25737984
	ZPRCENTRY421		zfft_r4dwpn,35M,_ac,17920,23436032		;,,,SKX
	DD			0
	PRCSTRT	667600000, 37748736, 0.306
	ZPRCENTRY421		zfft_r4dwpn,36M,_ac,32768,23145216
	ZPRCENTRY421		zfft_r4dwpn,36M,_ac,24576,25451264,,SKX
	ZPRCENTRY421		zfft_r4dwpn,36M,_ac,18432,24132352
	ZPRCENTRY421		zfft_r4dwpn,36M,_ac,16384,18434816
	ZPRCENTRY421		zfft_r4dwpn,36M,_ac,12288,17312512
	DD			0
	PRCSTRT	678800000, 38535168, 0.328
	ZPRCENTRY421		zfft_r4dwpn,37632K,_ac,28672,20464384
	ZPRCENTRY421		zfft_r4dwpn,37632K,_ac,25088,25942784,,SKX
	DD			0
	PRCSTRT	693600000, 39321600, 0.330
	ZPRCENTRY421		zfft_r4dwpn,38400K,_ac,30720,27572992
	ZPRCENTRY421		zfft_r4dwpn,38400K,_ac,20480,20220672,,,SKX
	ZPRCENTRY421		zfft_r4dwpn,38400K,_ac,12800,18021120
	DD			0
	PRCSTRT	728800000, 41287680, 0.350
	ZPRCENTRY421		zfft_r4dwpn,40320K,_ac,30720,21922560
	ZPRCENTRY421		zfft_r4dwpn,40320K,_ac,21504,21228288
	ZPRCENTRY421		zfft_r4dwpn,40320K,_ac,17920,20146944,,,SKX
	DD			0
	PRCSTRT	737200000, 41943040, 0.360
	ZPRCENTRY421		zfft_r4dwpn,40M,_ac,32768,29440768
	ZPRCENTRY421		zfft_r4dwpn,40M,_ac,20480,26778368,,,SKX
	DD			0
	PRCSTRT	750300000, 42467328, 0.365
	ZPRCENTRY421		zfft_r4dwpn,41472K,_ac,18432,20749056,,,SKX
	DD			0
	PRCSTRT	769900000, 44040192, 0.370
	ZPRCENTRY421		zfft_r4dwpn,42M,_ac,32768,23413504
	ZPRCENTRY421		zfft_r4dwpn,42M,_ac,28672,29645568
	ZPRCENTRY421		zfft_r4dwpn,42M,_ac,21504,28113664,,,SKX
	DD			0
	PRCSTRT	831700000, 47185920, 0.395
	ZPRCENTRY421		zfft_r4dwpn,45M,_ac,30720,31759104
	ZPRCENTRY421		zfft_r4dwpn,45M,_ac,24576,24283904
	ZPRCENTRY421		zfft_r4dwpn,45M,_ac,20480,23018240
	ZPRCENTRY421		zfft_r4dwpn,45M,_ac,15360,21609216,,,SKX
	DD			0
	PRCSTRT	845600000, 48168960, 0.400
	ZPRCENTRY421		zfft_r4dwpn,47040K,_ac,25088,24750848,,,SKX
	DD			0
	PRCSTRT	873000000, 49545216, 0.420
	ZPRCENTRY421		zfft_r4dwpn,48384K,_ac,21504,24165120,,,SKX
	DD			0
	PRCSTRT	885900000, 50331648, 0.430
	ZPRCENTRY421		zfft_r4dwpn,48M,_ac,32768,33905408
	ZPRCENTRY421		zfft_r4dwpn,48M,_ac,24576,32152320
	ZPRCENTRY421		zfft_r4dwpn,48M,_ac,16384,23046912,,,SKX
	DD			0
	PRCSTRT	898300000, 51380224, 0.470
	ZPRCENTRY421		zfft_r4dwpn,49M,_ac,25088,32783104,,,SKX
	DD			0
	PRCSTRT	965200000, 55050240, 0.500
	ZPRCENTRY421		zfft_r4dwpn,53760K,_ac,28672,28281600
	ZPRCENTRY421		zfft_r4dwpn,53760K,_ac,17920,25189120,,,SKX
	DD			0
	PRCSTRT	996200000, 56623104, 0.480
	ZPRCENTRY421		zfft_r4dwpn,55296K,_ac,24576,27638528
	ZPRCENTRY421		zfft_r4dwpn,55296K,_ac,18432,25934592,,,SKX
	DD			0
	PRCSTRT	1013000000, 57802752, 0.520
	ZPRCENTRY421		zfft_r4dwpn,56448K,_ac,25088,28175104,,,SKX
	DD			0
	PRCSTRT	1024000000, 58720256, 0.530
	ZPRCENTRY421		zfft_r4dwpn,56M,_ac,28672,37460736,,,SKX
	DD			0
	PRCSTRT	1036000000, 58982400, 0.540
	ZPRCENTRY421		zfft_r4dwpn,57600K,_ac,30720,30296832,,,SKX
	DD			0
	PRCSTRT	1101000000, 62914560, 0.590
	ZPRCENTRY421		zfft_r4dwpn,60M,_ac,32768,32344832
	ZPRCENTRY421		zfft_r4dwpn,60M,_ac,30720,40131328
	ZPRCENTRY421		zfft_r4dwpn,60M,_ac,20480,28777216,,,SKX
	DD			0
	PRCSTRT	1168000000, 67108864, 0.640
	ZPRCENTRY421		zfft_r4dwpn,64M,_ac,32768,42834688,,,SKX
	DD			0
	DD	0
ENDIF

	;; Align so that other GWDATA areas are also aligned on a cache line
	align 128
_GWDATA ENDS

;; FFT setup routines

_TEXT SEGMENT

; gwinfo1 (resptr)
;	Return address of jmp tables for C code to examine
; Windows 32-bit (_gwinfo1)
; Linux 32-bit (gwinfo1)
;	Parameter resptr = [esp+4]
; Windows 64-bit (gwinfo1) - leaf routine, no unwind info necessary
;	Parameter resptr = rcx
; Linux 64-bit (gwinfo1)
;	Parameter resptr = rdi

PROCL	gwinfo1
	IFNDEF X86_64
	mov	ecx, [esp+4]		; Address of data struct to return info
	ENDIF
	IFDEF LINUX64
	mov	rcx, rdi		; Address of data struct to return info
	ENDIF
	mov	rax, OFFSET xjmptable	; SSE2 mersenne mod FFTs
	mov	[rcx+0*SZPTR], rax
	mov	rax, OFFSET xjmptablep	; SSE2 2^N+1 mod FFTs
	mov	[rcx+1*SZPTR], rax
	IFNDEF X86_64
	mov	rax, OFFSET jmptable	; x86 mersenne mod FFTs
	mov	[rcx+2*SZPTR], rax
	mov	rax, OFFSET jmptablep	; x86 2^N+1 mod FFTs
	mov	[rcx+3*SZPTR], rax
	ENDIF
	mov	rax, OFFSET yjmptable	; AVX mersenne mod FFTs
	mov	[rcx+4*SZPTR], rax
	mov	rax, OFFSET yjmptablep	; AVX 2^N+1 mod FFTs
	mov	[rcx+5*SZPTR], rax
	IFDEF X86_64
	mov	rax, OFFSET zjmptable	; AVX-512 mersenne mod FFTs
	mov	[rcx+6*SZPTR], rax
	mov	rax, OFFSET zjmptablep	; AVX-512 2^N+1 mod FFTs
	mov	[rcx+7*SZPTR], rax
	ENDIF
	mov	eax, VERSION_NUMBER
	mov	[rcx+8*SZPTR], eax
	ret
gwinfo1 ENDP

_TEXT	ENDS
END
