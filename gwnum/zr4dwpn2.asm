; Copyright 2011-2021 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the AVX-512 radix-4/8 DJB FFT with delayed sin/cos multiplies and partial normalization.
;

	TITLE   setup

zfft_type TEXTEQU <r4dwpn>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE zr4.mac
INCLUDE zr4dwpnpass1sc.mac
INCLUDE zr4dwpnpass2.mac

;;EXTRN zrfoo:PROC

_TEXT SEGMENT

PROCFP	zreal_fft_final
	int_prolog 0,0,0
	zr64_hundredtwentyeight_real_fft_final_preload
	zr64_hundredtwentyeight_real_fft_final rsi, 8*128, 1*128, 2*128, 4*128, rdx, 0*ZMM_SCD4, rdi, 0*ZMM_SCD4, 1
	int_epilog 0,0,0
ENDPP	zreal_fft_final 

PROCFP	zreal_square
	int_prolog 0,0,0
	zr64_hundredtwentyeight_real_with_square_preload
	zr64_hundredtwentyeight_real_with_square rsi, 8*128, 1*128, 2*128, 4*128, rdx, 0*ZMM_SCD4, rdi, 0*ZMM_SCD4, 1
	int_epilog 0,0,0
ENDPP	zreal_square 

PROCFP	zreal_mult
	int_prolog 0,0,0
	zr64_hundredtwentyeight_real_with_mult_preload
	zr64_hundredtwentyeight_real_with_mult rsi, 8*128, 1*128, 2*128, 4*128, rdx, 0*ZMM_SCD4, rdi, 0*ZMM_SCD4, 1
	int_epilog 0,0,0
ENDPP	zreal_mult 

PROCFP	zreal_mulf
	int_prolog 0,0,0
	zr64f_hundredtwentyeight_real_with_mulf_preload
	zr64f_hundredtwentyeight_real_with_mulf rsi, 8*128, 1*128, 2*128, 4*128, rdx, 0*ZMM_SCD4, rdi, 0*ZMM_SCD4, 1
	int_epilog 0,0,0
ENDPP	zreal_mulf 

PROCF	zcomplex_square_opcode
	int_prolog 0,0,0
	zr64_64c_square_opcode rsi, 1*128, 2*128, 4*128
	int_epilog 0,0,0
zcomplex_square_opcode ENDP

PROCF	zcomplex_mult_opcode
	int_prolog 0,0,0
	zr64_64c_mult_opcode rsi, 1*128, 2*128, 4*128
	int_epilog 0,0,0
zcomplex_mult_opcode ENDP

PROCF	zcomplex_mulf_opcode
	int_prolog 0,0,0
	zr64_64c_mulf_opcode rsi, 1*128, 2*128, 4*128
	int_epilog 0,0,0
zcomplex_mulf_opcode ENDP

;; Generate real and complex pass 2 routines optimized for this architecture

buildfor SKX,	zpass2gen 64
buildfor SKX,	zpass2gen 320
buildfor SKX,	zpass2gen 384
buildfor SKX,	zpass2gen 448
buildfor SKX,	zpass2gen 512
buildfor SKX,	zpass2gen 640
buildfor SKX,	zpass2gen 768
buildfor SKX,	zpass2gen 1024
buildfor SKX,	zpass2gen 1600
buildfor SKX,	zpass2gen 1920
buildfor SKX,	zpass2gen 2240
buildfor SKX,	zpass2gen 2304
buildfor SKX,	zpass2gen 2560
buildfor SKX,	zpass2gen 2688
buildfor SKX,	zpass2gen 3072
buildfor SKX,	zpass2gen 3136
buildfor SKX,	zpass2gen 3200
buildfor SKX,	zpass2gen 3584
buildfor SKX,	zpass2gen 3840
buildfor SKX,	zpass2gen 4096
buildfor SKX,	zpass2gen 4480
buildfor SKX,	zpass2gen 4608
buildfor SKX,	zpass2gen 5120
buildfor SKX,	zpass2gen 5376
buildfor SKX,	zpass2gen 6144
buildfor SKX,	zpass2gen 6400
buildfor SKX,	zpass2gen 7168
buildfor SKX,	zpass2gen 7680
buildfor SKX,	zpass2gen 8192
buildfor SKX,	zpass2gen 9216
buildfor SKX,	zpass2gen 10240
buildfor SKX,	zpass2gen 12288
buildfor SKX,	zpass2gen 12800
buildfor SKX,	zpass2gen 15360
buildfor SKX,	zpass2gen 16384
buildfor SKX,	zpass2gen 17920
buildfor SKX,	zpass2gen 18432
buildfor SKX,	zpass2gen 20480
buildfor SKX,	zpass2gen 21504
buildfor SKX,	zpass2gen 24576
buildfor SKX,	zpass2gen 25088
buildfor SKX,	zpass2gen 28672
buildfor SKX,	zpass2gen 30720
buildfor SKX,	zpass2gen 32768

_TEXT	ENDS
END
