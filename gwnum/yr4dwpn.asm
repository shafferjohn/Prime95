; Copyright 2011-2020 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the AVX radix-4/8 DJB FFT with delayed sin/cos multiplies and partial normalization.
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

;; To reduce executable size, we call common code to do the very few real FFT words in the last FFT level

IF (@INSTR(,%yarch,<CORE>) NE 0)			;; Only build these once -- contains no FMA3-specific code

IFNDEF X86_64
PROCF	yreal_fft_final1a				;; One pass version has no current push_amt
	int_prolog 0,0,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 4*64, 64, 2*64
	int_epilog 0,0,0
yreal_fft_final1a ENDP
ENDIF
PROCF	yreal_fft_final1
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 4*64, 64, 2*64
	int_epilog 0,16,0
yreal_fft_final1 ENDP
PROCF	yreal_fft_final2
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 64, 4*64, 8*64
	int_epilog 0,16,0
yreal_fft_final2 ENDP
PROCF	yreal_fft_final3
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 64, 16*64, 32*64
	int_epilog 0,16,0
yreal_fft_final3 ENDP
PROCF	yreal_fft_final4
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 64, (64*64+64), 2*(64*64+64)
	int_epilog 0,16,0
yreal_fft_final4 ENDP
PROCF	yreal_fft_final4a
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 64, (64*64+128), 2*(64*64+128)
	int_epilog 0,16,0
yreal_fft_final4a ENDP
PROCF	yreal_fft_final5
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 64, 4*(64*64+64), 8*(64*64+64)
	int_epilog 0,16,0
yreal_fft_final5 ENDP
PROCF	yreal_fft_final6
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_fft_final rsi, 64, 16*(64*64+64), 32*(64*64+64)
	int_epilog 0,16,0
yreal_fft_final6 ENDP

IFNDEF X86_64
PROCF	yreal_square1a				;; One pass version has no current push_amt
	int_prolog 0,0,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 4*64, 64, 2*64
	int_epilog 0,0,0
yreal_square1a ENDP
ENDIF
PROCF	yreal_square1
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 4*64, 64, 2*64
	int_epilog 0,16,0
yreal_square1 ENDP
PROCF	yreal_square2
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 64, 4*64, 8*64
	int_epilog 0,16,0
yreal_square2 ENDP
PROCF	yreal_square3
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 64, 16*64, 32*64
	int_epilog 0,16,0
yreal_square3 ENDP
PROCF	yreal_square4
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 64, (64*64+64), 2*(64*64+64)
	int_epilog 0,16,0
yreal_square4 ENDP
PROCF	yreal_square4a
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 64, (64*64+128), 2*(64*64+128)
	int_epilog 0,16,0
yreal_square4a ENDP
PROCF	yreal_square5
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 64, 4*(64*64+64), 8*(64*64+64)
	int_epilog 0,16,0
yreal_square5 ENDP
PROCF	yreal_square6
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_square rsi, 64, 16*(64*64+64), 32*(64*64+64)
	int_epilog 0,16,0
yreal_square6 ENDP

IFNDEF X86_64
PROCF	yreal_mult1a				;; One pass version has no current push_amt
	int_prolog 0,0,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 4*64, 64, 2*64
	int_epilog 0,0,0
yreal_mult1a ENDP
ENDIF
PROCF	yreal_mult1
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 4*64, 64, 2*64
	int_epilog 0,16,0
yreal_mult1 ENDP
PROCF	yreal_mult2
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 64, 4*64, 8*64
	int_epilog 0,16,0
yreal_mult2 ENDP
PROCF	yreal_mult3
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 64, 16*64, 32*64
	int_epilog 0,16,0
yreal_mult3 ENDP
PROCF	yreal_mult4
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 64, (64*64+64), 2*(64*64+64)
	int_epilog 0,16,0
yreal_mult4 ENDP
PROCF	yreal_mult4a
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 64, (64*64+128), 2*(64*64+128)
	int_epilog 0,16,0
yreal_mult4a ENDP
PROCF	yreal_mult5
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 64, 4*(64*64+64), 8*(64*64+64)
	int_epilog 0,16,0
yreal_mult5 ENDP
PROCF	yreal_mult6
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mult rsi, 64, 16*(64*64+64), 32*(64*64+64)
	int_epilog 0,16,0
yreal_mult6 ENDP

IFNDEF X86_64
PROCF	yreal_mulf1a				;; One pass version has no current push_amt
	int_prolog 0,0,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 4*64, 64, 2*64
	int_epilog 0,0,0
yreal_mulf1a ENDP
ENDIF
PROCF	yreal_mulf1
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 4*64, 64, 2*64
	int_epilog 0,16,0
yreal_mulf1 ENDP
PROCF	yreal_mulf2
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 64, 4*64, 8*64
	int_epilog 0,16,0
yreal_mulf2 ENDP
PROCF	yreal_mulf3
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 64, 16*64, 32*64
	int_epilog 0,16,0
yreal_mulf3 ENDP
PROCF	yreal_mulf4
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 64, (64*64+64), 2*(64*64+64)
	int_epilog 0,16,0
yreal_mulf4 ENDP
PROCF	yreal_mulf4a
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 64, (64*64+128), 2*(64*64+128)
	int_epilog 0,16,0
yreal_mulf4a ENDP
PROCF	yreal_mulf5
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 64, 4*(64*64+64), 8*(64*64+64)
	int_epilog 0,16,0
yreal_mulf5 ENDP
PROCF	yreal_mulf6
	int_prolog 0,16,0
	yr4_4cl_eight_reals_four_complex_with_mulf rsi, 64, 16*(64*64+64), 32*(64*64+64)
	int_epilog 0,16,0
yreal_mulf6 ENDP


PROCF	y8real_fft_final1
	int_prolog 0,0,0
	yr8_8cl_sixteen_reals_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	int_epilog 0,0,0
y8real_fft_final1 ENDP
PROCF	y8real_fft_final2
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_fft_final rsi, 2*64, 32*64, (64*64+128), 64
	int_epilog 0,16,0
y8real_fft_final2 ENDP
PROCF	y8real_fft_final3
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_fft_final rsi, 2*64, 2*(64*64+64), 4*(64*64+64), 64
	int_epilog 0,16,0
y8real_fft_final3 ENDP
PROCF	y8real_fft_final4
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_fft_final rsi, 2*64, 8*(64*64+64), 16*(64*64+64), 64
	int_epilog 0,16,0
y8real_fft_final4 ENDP

PROCF	y8real_square1
	int_prolog 0,0,0
	yr8_8cl_sixteen_reals_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	int_epilog 0,0,0
y8real_square1 ENDP
PROCF	y8real_square2
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_square rsi, 2*64, 32*64, (64*64+128), 64
	int_epilog 0,16,0
y8real_square2 ENDP
PROCF	y8real_square3
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_square rsi, 2*64, 2*(64*64+64), 4*(64*64+64), 64
	int_epilog 0,16,0
y8real_square3 ENDP
PROCF	y8real_square4
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_square rsi, 2*64, 8*(64*64+64), 16*(64*64+64), 64
	int_epilog 0,16,0
y8real_square4 ENDP

PROCF	y8real_mult1
	int_prolog 0,0,0
	yr8_8cl_sixteen_reals_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	int_epilog 0,0,0
y8real_mult1 ENDP
PROCF	y8real_mult2
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_mult rsi, 2*64, 32*64, (64*64+128), 64
	int_epilog 0,16,0
y8real_mult2 ENDP
PROCF	y8real_mult3
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_mult rsi, 2*64, 2*(64*64+64), 4*(64*64+64), 64
	int_epilog 0,16,0
y8real_mult3 ENDP
PROCF	y8real_mult4
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_mult rsi, 2*64, 8*(64*64+64), 16*(64*64+64), 64
	int_epilog 0,16,0
y8real_mult4 ENDP

PROCF	y8real_mulf1
	int_prolog 0,0,0
	yr8_8cl_sixteen_reals_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64
	int_epilog 0,0,0
y8real_mulf1 ENDP
PROCF	y8real_mulf2
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_mulf rsi, 2*64, 32*64, (64*64+128), 64
	int_epilog 0,16,0
y8real_mulf2 ENDP
PROCF	y8real_mulf3
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_mulf rsi, 2*64, 2*(64*64+64), 4*(64*64+64), 64
	int_epilog 0,16,0
y8real_mulf3 ENDP
PROCF	y8real_mulf4
	int_prolog 0,16,0
	yr8_8cl_sixteen_reals_eight_complex_with_mulf rsi, 2*64, 8*(64*64+64), 16*(64*64+64), 64
	int_epilog 0,16,0
y8real_mulf4 ENDP

;; Generate code to support gwaddmul4 and gwsubmul4 with radix-8 last step.  There aren't FMA3 versions of these routines

IFDEF X86_64
PROCF	y8complex_mult_opcode1
	int_prolog 0,16,0
	yr8_8c_mult_opcode rsi, 64, 2*64, 4*64
	int_epilog 0,16,0
y8complex_mult_opcode1 ENDP
PROCF	y8complex_mult_opcode2
	int_prolog 0,16,0
	yr8_8c_mult_opcode rsi, 32*64, (64*64+128), 64
	int_epilog 0,16,0
y8complex_mult_opcode2 ENDP
PROCF	y8complex_mult_opcode3
	int_prolog 0,16,0
	yr8_8c_mult_opcode rsi, 2*(64*64+64), 4*(64*64+64), 64
	int_epilog 0,16,0
y8complex_mult_opcode3 ENDP
PROCF	y8complex_mult_opcode4
	int_prolog 0,16,0
	yr8_8c_mult_opcode rsi, 8*(64*64+64), 16*(64*64+64), 64
	int_epilog 0,16,0
y8complex_mult_opcode4 ENDP

PROCF	y8complex_mulf_opcode1
	int_prolog 0,16,0
	yr8_8c_mulf_opcode rsi, 64, 2*64, 4*64
	int_epilog 0,16,0
y8complex_mulf_opcode1 ENDP
PROCF	y8complex_mulf_opcode2
	int_prolog 0,16,0
	yr8_8c_mulf_opcode rsi, 32*64, (64*64+128), 64
	int_epilog 0,16,0
y8complex_mulf_opcode2 ENDP
PROCF	y8complex_mulf_opcode3
	int_prolog 0,16,0
	yr8_8c_mulf_opcode rsi, 2*(64*64+64), 4*(64*64+64), 64
	int_epilog 0,16,0
y8complex_mulf_opcode3 ENDP
PROCF	y8complex_mulf_opcode4
	int_prolog 0,16,0
	yr8_8c_mulf_opcode rsi, 8*(64*64+64), 16*(64*64+64), 64
	int_epilog 0,16,0
y8complex_mulf_opcode4 ENDP
ENDIF

ELSE

EXTRN	yreal_fft_final1:PROC
EXTRN	yreal_fft_final2:PROC
EXTRN	yreal_fft_final3:PROC
EXTRN	yreal_fft_final4:PROC
EXTRN	yreal_fft_final4a:PROC
EXTRN	yreal_fft_final5:PROC
EXTRN	yreal_fft_final6:PROC
EXTRN	yreal_square1:PROC
EXTRN	yreal_square2:PROC
EXTRN	yreal_square3:PROC
EXTRN	yreal_square4:PROC
EXTRN	yreal_square4a:PROC
EXTRN	yreal_square5:PROC
EXTRN	yreal_square6:PROC
EXTRN	yreal_mult1:PROC
EXTRN	yreal_mult2:PROC
EXTRN	yreal_mult3:PROC
EXTRN	yreal_mult4:PROC
EXTRN	yreal_mult4a:PROC
EXTRN	yreal_mult5:PROC
EXTRN	yreal_mult6:PROC
EXTRN	yreal_mulf1:PROC
EXTRN	yreal_mulf2:PROC
EXTRN	yreal_mulf3:PROC
EXTRN	yreal_mulf4:PROC
EXTRN	yreal_mulf4a:PROC
EXTRN	yreal_mulf5:PROC
EXTRN	yreal_mulf6:PROC
EXTRN	y8real_fft_final1:PROC
EXTRN	y8real_fft_final2:PROC
EXTRN	y8real_fft_final3:PROC
EXTRN	y8real_fft_final4:PROC
EXTRN	y8real_square1:PROC
EXTRN	y8real_square2:PROC
EXTRN	y8real_square3:PROC
EXTRN	y8real_square4:PROC
EXTRN	y8real_mult1:PROC
EXTRN	y8real_mult2:PROC
EXTRN	y8real_mult3:PROC
EXTRN	y8real_mult4:PROC
EXTRN	y8real_mulf1:PROC
EXTRN	y8real_mulf2:PROC
EXTRN	y8real_mulf3:PROC
EXTRN	y8real_mulf4:PROC

IFDEF X86_64
EXTRN	y8complex_mult_opcode1:PROC
EXTRN	y8complex_mult_opcode2:PROC
EXTRN	y8complex_mult_opcode3:PROC
EXTRN	y8complex_mult_opcode4:PROC
EXTRN	y8complex_mulf_opcode1:PROC
EXTRN	y8complex_mulf_opcode2:PROC
EXTRN	y8complex_mulf_opcode3:PROC
EXTRN	y8complex_mulf_opcode4:PROC
ENDIF

ENDIF

;; Generate the AVX and FMA3 versions of the radix-4 opcode routines

IFDEF X86_64
PROCFP	ycomplex_mult_opcode1
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, 64, 2*64
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode1
PROCFP	ycomplex_mult_opcode2
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, 4*64, 8*64
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode2
PROCFP	ycomplex_mult_opcode3
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, 16*64, 32*64
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode3
PROCFP	ycomplex_mult_opcode4
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, (64*64+64), 2*(64*64+64)
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode4
PROCFP	ycomplex_mult_opcode4a
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, (64*64+128), 2*(64*64+128)
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode4a
PROCFP	ycomplex_mult_opcode5
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, 4*(64*64+64), 8*(64*64+64)
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode5
PROCFP	ycomplex_mult_opcode6
	int_prolog 0,16,0
	yr4_4c_mult_opcode rsi, 16*(64*64+64), 32*(64*64+64)
	int_epilog 0,16,0
ENDPP	ycomplex_mult_opcode6

PROCFP	ycomplex_mulf_opcode1
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, 64, 2*64
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode1
PROCFP	ycomplex_mulf_opcode2
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, 4*64, 8*64
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode2
PROCFP	ycomplex_mulf_opcode3
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, 16*64, 32*64
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode3
PROCFP	ycomplex_mulf_opcode4
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, (64*64+64), 2*(64*64+64)
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode4
PROCFP	ycomplex_mulf_opcode4a
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, (64*64+128), 2*(64*64+128)
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode4a
PROCFP	ycomplex_mulf_opcode5
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, 4*(64*64+64), 8*(64*64+64)
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode5
PROCFP	ycomplex_mulf_opcode6
	int_prolog 0,16,0
	yr4_4c_mulf_opcode rsi, 16*(64*64+64), 32*(64*64+64)
	int_epilog 0,16,0
ENDPP	ycomplex_mulf_opcode6
ENDIF

;; Generate shared pass 1 routines optimized for this architecture

IFDEF X86_64
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 128,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 256,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 320,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 384,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 448,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 512,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 640,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 768,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 896,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1024,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1280,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1536,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1792,0,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 2048,0,
ENDIF

;; Generate pass 2 routines optimized for this architecture

buildfor CORE    + FMA3_64,	ypass2gen 48
buildfor CORE    + FMA3_64,	ypass2gen 6
buildfor CORE    + FMA3_64,	ypass2gen 80
buildfor CORE    + FMA3_64,	ypass2gen 192
buildfor CORE    + FMA3_64,	ypass2gen 8
buildfor CORE    + FMA3_64,	ypass2gen 320
buildfor CORE    + FMA3_64,	ypass2gen 768
buildfor CORE    + FMA3_64,	ypass2gen 10
buildfor CORE    + FMA3_64,	ypass2gen 1280
buildfor CORE    + FMA3_64,	ypass2gen 1536
buildfor CORE    + FMA3_64,	ypass2gen 11
buildfor CORE    + FMA3_64,	ypass2gen 2560

;; Routines for many FFT sizes

	; The 6 levels variants (48, 64, 80)

IFNDEF X86_64
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc128 6K, 48
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc128 8K, 6
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc128 10K, 80

build421 ,			,			CORE + FMA3_64,		yr4dwpn_pass1sc256 12K, 48
build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc256 16K, 6
build421 ,			,			,			yr4dwpn_pass1sc256 20K, 80

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc320 15K, 48
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc320 20K, 6
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc320 25K, 80

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc384 18K, 48
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc384 24K, 6
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc384 30K, 80

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc448 21K, 48
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc448 28K, 6
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc448 35K, 80

build421 ,			,			,			yr4dwpn_pass1sc512 24K, 48
build421 ,			,			,			yr4dwpn_pass1sc512 32K, 6
build421 ,			,			CORE_64 + FMA3_64,	yr4dwpn_pass1sc512 40K, 80

build421 ,			,			,			yr4dwpn_pass1sc640 30K, 48
build421 ,			,			,			yr4dwpn_pass1sc640 40K, 6
build421 ,			,			CORE + FMA3_64,		yr4dwpn_pass1sc640 50K, 80

build421 ,			,			CORE + FMA3_64,		yr4dwpn_pass1sc768 36K, 48
build421 ,			,			,			yr4dwpn_pass1sc768 48K, 6
build421 ,			,			,			yr4dwpn_pass1sc768 60K, 80
ENDIF

	; The 8 levels variants (192, 256, 320)

IFNDEF X86_64
build421 ,			CORE,			,			yr4dwpn_pass1sc128 24K, 192
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc128 32K, 8
build421 ,			CORE_32,		,			yr4dwpn_pass1sc128 40K, 320

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc256 48K, 192
build421 ,			CORE_32,		CORE_64 + FMA3_64,	yr4dwpn_pass1sc256 64K, 8
build421 ,			,			,			yr4dwpn_pass1sc256 80K, 320

build421 ,			CORE + FMA3_64,		CORE_64,		yr4dwpn_pass1sc320 60K, 192
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc320 80K, 8
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc320 100K, 320

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc384 72K, 192
build421 ,			CORE_64 + FMA3_64,	CORE_32,		yr4dwpn_pass1sc384 96K, 8
build421 ,			CORE_64 + FMA3_64,	CORE_32,		yr4dwpn_pass1sc384 120K, 320

build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc448 84K, 192
build421 ,			CORE_64 + FMA3_64,	CORE_32,		yr4dwpn_pass1sc448 112K, 8
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc448 140K, 320

build421 ,			,			,			yr4dwpn_pass1sc512 96K, 192
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc512 128K, 8
build421 ,			,			,			yr4dwpn_pass1sc512 160K, 320

build421 ,			,			,			yr4dwpn_pass1sc640 120K, 192
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc640 160K, 8
build421 ,			FMA3_64,		CORE,			yr4dwpn_pass1sc640 200K, 320

build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc768 144K, 192
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc768 192K, 8
build421 ,			,			,			yr4dwpn_pass1sc768 240K, 320

build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc896 168K, 192
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc896 224K, 8
build421 ,			,			,			yr4dwpn_pass1sc896 280K, 320

build421 ,			,			,			yr4dwpn_pass1sc1024 192K, 192
build421 ,			,			,			yr4dwpn_pass1sc1024 256K, 8
build421 ,			,			,			yr4dwpn_pass1sc1024 320K, 320

	; The 10 levels variants (768, 1024, 1280)

build421 ,			,			,			yr4dwpn_pass1sc128 96K, 768
build421 ,			CORE_32,		,			yr4dwpn_pass1sc128 128K, 10
build421 ,			CORE_32,		,			yr4dwpn_pass1sc128 160K, 1280

build421 ,			,			,			yr4dwpn_pass1sc256 192K, 768
build421 ,			CORE_32,		CORE_64,		yr4dwpn_pass1sc256 256K, 10
build421 ,			,			,			yr4dwpn_pass1sc256 320K, 1280

build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc320 240K, 768
build421 FMA3_64,		CORE,			,			yr4dwpn_pass1sc320 320K, 10
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc320 400K, 1280

build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc384 288K, 768
build421 ,			CORE,			,			yr4dwpn_pass1sc384 384K, 10
build421 ,			CORE,			,			yr4dwpn_pass1sc384 480K, 1280

build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc448 336K, 768
build421 FMA3_64,		CORE,			,			yr4dwpn_pass1sc448 448K, 10
build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc448 560K, 1280

build421 ,			,			,			yr4dwpn_pass1sc512 384K, 768
build421 ,			,			,			yr4dwpn_pass1sc512 512K, 10
build421 ,			,			,			yr4dwpn_pass1sc512 640K, 1280

build421 ,			,			,			yr4dwpn_pass1sc640 480K, 768
build421 ,			CORE_32,		,			yr4dwpn_pass1sc640 640K, 10
build421 ,			CORE_32,		,			yr4dwpn_pass1sc640 800K, 1280

build421 ,			,			,			yr4dwpn_pass1sc768 576K, 768
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768 768K, 10
build421 ,			,			,			yr4dwpn_pass1sc768 960K, 1280

build421 ,			,			,			yr4dwpn_pass1sc896 672K, 768
build421 ,			CORE_32,		,			yr4dwpn_pass1sc896 896K, 10
build421 ,			CORE_32,		,			yr4dwpn_pass1sc896 1120K, 1280

build421 ,			,			,			yr4dwpn_pass1sc1024 768K, 768
build421 ,			CORE_32,		,			yr4dwpn_pass1sc1024 1M, 10
build421 ,			,			,			yr4dwpn_pass1sc1024 1280K, 1280

build421 ,			,			,			yr4dwpn_pass1sc1280 960K, 768
build421 ,			,			,			yr4dwpn_pass1sc1280 1280K, 10
build421 ,			,			,			yr4dwpn_pass1sc1280 1600K, 1280

build421 ,			,			,			yr4dwpn_pass1sc1536 1536K, 10

build421 ,			,			,			yr4dwpn_pass1sc1792 1792K, 10

build421 ,			,			,			yr4dwpn_pass1sc2048 2M, 10

;build421 ,			,			,			yr4dwpn_pass1sc2560 2560K, 10

;build421 ,			,			,			yr4dwpn_pass1sc3072 3M, 10

;build421 ,			,			,			yr4dwpn_pass1sc3584 3584K, 10

;build421 ,			,			,			yr4dwpn_pass1sc4096 4M, 10
ENDIF

	; The 11 levels variants (1536, 2048, 2560)

IFNDEF X86_64
build421 ,			,			,			yr4dwpn_pass1sc128 192K, 1536
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc128 256K, 11
build421 ,			,			,			yr4dwpn_pass1sc128 320K, 2560

build421 ,			,			FMA3_64,		yr4dwpn_pass1sc256 384K, 1536
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc256 512K, 11
build421 ,			,			,			yr4dwpn_pass1sc256 640K, 2560

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc320 480K, 1536
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc320 640K, 11
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc320 800K, 2560

build421 ,			CORE_32 + FMA3_64,	CORE_64,		yr4dwpn_pass1sc384 576K, 1536
build421 ,			FMA3_64,		CORE_64,		yr4dwpn_pass1sc384 768K, 11
build421 ,			CORE_32,		CORE_64,		yr4dwpn_pass1sc384 960K, 2560

build421 ,			CORE + FMA3_64,		,			yr4dwpn_pass1sc448 672K, 1536
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc448 896K, 11
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc448 1120K, 2560

build421 ,			,			,			yr4dwpn_pass1sc512 768K, 1536
build421 ,			,			,			yr4dwpn_pass1sc512 1M, 11
build421 ,			,			,			yr4dwpn_pass1sc512 1280K, 2560

build421 ,			,			,			yr4dwpn_pass1sc640 960K, 1536
build421 ,			,			,			yr4dwpn_pass1sc640 1280K, 11
build421 ,			CORE_64,		,			yr4dwpn_pass1sc640 1600K, 2560

build421 ,			,			,			yr4dwpn_pass1sc768 1152K, 1536
build421 ,			,			,			yr4dwpn_pass1sc768 1536K, 11
build421 ,			CORE_64,		,			yr4dwpn_pass1sc768 1920K, 2560

build421 ,			CORE_32,		,			yr4dwpn_pass1sc896 1344K, 1536
build421 ,			,			,			yr4dwpn_pass1sc896 1792K, 11
build421 ,			,			,			yr4dwpn_pass1sc896 2240K, 2560

build421 ,			,			,			yr4dwpn_pass1sc1024 1536K, 1536
build421 ,			,			,			yr4dwpn_pass1sc1024 2M, 11
build421 ,			,			,			yr4dwpn_pass1sc1024 2560K, 2560

build421 ,			,			,			yr4dwpn_pass1sc1280 1920K, 1536
build421 ,			,			,			yr4dwpn_pass1sc1280 2560K, 11
build421 ,			,			,			yr4dwpn_pass1sc1280 3200K, 2560

build421 ,			,			,			yr4dwpn_pass1sc1536 3M, 11

build421 ,			,			,			yr4dwpn_pass1sc1792 3584K, 11

build421 ,			,			,			yr4dwpn_pass1sc2048 4M, 11

;build421 ,			,			,			yr4dwpn_pass1sc2560 5M, 11

;build421 ,			,			,			yr4dwpn_pass1sc3072 6M, 11

;build421 ,			,			,			yr4dwpn_pass1sc3584 7M, 11

;build421 ,			,			,			yr4dwpn_pass1sc4096 8M, 11
ENDIF

;; Generate shared pass 1 all-complex routines optimized for this architecture

IFDEF X86_64
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 128,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 256,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 384,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 512,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 640,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 768,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1024,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1280,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 1536,1,
build421 CORE    + FMA3_64, CORE    + FMA3_64, CORE    + FMA3_64,	ypass1gen 2048,1,
ENDIF

	; The all-complex 6 levels variants (48, 64, 80)

IFNDEF X86_64
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc128ac 6K, 48
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc128ac 8K, 6
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc128ac 10K, 80

build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc256ac 12K, 48
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc256ac 16K, 6
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc256ac 20K, 80

build421 ,			CORE_64 + FMA3_64,	CORE_32,		yr4dwpn_pass1sc384ac 18K, 48
build421 ,			,			,			yr4dwpn_pass1sc384ac 24K, 6
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc384ac 30K, 80

build421 ,			,			,			yr4dwpn_pass1sc512ac 24K, 48
build421 ,			,			,			yr4dwpn_pass1sc512ac 32K, 6
build421 ,			,			,			yr4dwpn_pass1sc512ac 40K, 80

build421 ,			,			,			yr4dwpn_pass1sc640ac 30K, 48
build421 ,			,			,			yr4dwpn_pass1sc640ac 40K, 6
build421 ,			CORE_64 + FMA3_64,	CORE_32,		yr4dwpn_pass1sc640ac 50K, 80

build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc768ac 36K, 48
build421 ,			,			,			yr4dwpn_pass1sc768ac 48K, 6
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc768ac 60K, 80
ENDIF

	; The all-complex 8 levels variants (192, 256, 320)

IFNDEF X86_64
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc128ac 24K, 192
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc128ac 32K, 8
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc128ac 40K, 320

build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc256ac 48K, 192
build421 FMA3_64,		CORE_64,		CORE_32,		yr4dwpn_pass1sc256ac 64K, 8
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc256ac 80K, 320

build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc384ac 72K, 192
build421 ,			,			CORE_32,		yr4dwpn_pass1sc384ac 96K, 8
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc384ac 120K, 320

build421 ,			,			,			yr4dwpn_pass1sc512ac 96K, 192
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc512ac 128K, 8
build421 ,			,			CORE_32,		yr4dwpn_pass1sc512ac 160K, 320

build421 ,			,			,			yr4dwpn_pass1sc640ac 120K, 192
build421 ,			,			,			yr4dwpn_pass1sc640ac 160K, 8
build421 ,			CORE_64 + FMA3_64,	CORE_32,		yr4dwpn_pass1sc640ac 200K, 320

build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc768ac 144K, 192
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768ac 192K, 8
build421 ,			,			,			yr4dwpn_pass1sc768ac 240K, 320

build421 ,			,			,			yr4dwpn_pass1sc1024ac 192K, 192
build421 ,			,			,			yr4dwpn_pass1sc1024ac 256K, 8
build421 ,			,			,			yr4dwpn_pass1sc1024ac 320K, 320

build421 ,			CORE_64,		,			yr4dwpn_pass1sc1280ac 240K, 192
build421 ,			,			,			yr4dwpn_pass1sc1280ac 320K, 8
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc1280ac 400K, 320
ENDIF

	; The all-complex 10 levels variants (768, 1024, 1280)

IFNDEF X86_64
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc128ac 96K, 768
build421 ,			CORE,			,			yr4dwpn_pass1sc128ac 128K, 10
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc128ac 160K, 1280

build421 ,			CORE_64,		,			yr4dwpn_pass1sc256ac 192K, 768
build421 ,			CORE,			,			yr4dwpn_pass1sc256ac 256K, 10
build421 ,			CORE_64,		,			yr4dwpn_pass1sc256ac 320K, 1280

build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc384ac 288K, 768
build421 ,			CORE,			,			yr4dwpn_pass1sc384ac 384K, 10
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc384ac 480K, 1280

build421 ,			FMA3_64,		,			yr4dwpn_pass1sc512ac 384K, 768
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc512ac 512K, 10
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc512ac 640K, 1280

build421 ,			,			,			yr4dwpn_pass1sc640ac 480K, 768
build421 ,			,			,			yr4dwpn_pass1sc640ac 640K, 10
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc640ac 800K, 1280

build421 ,			,			,			yr4dwpn_pass1sc768ac 576K, 768
build421 ,			CORE_32,		,			yr4dwpn_pass1sc768ac 768K, 10
build421 ,			,			,			yr4dwpn_pass1sc768ac 960K, 1280

build421 ,			,			,			yr4dwpn_pass1sc1024ac 768K, 768
build421 ,			,			,			yr4dwpn_pass1sc1024ac 1M, 10
build421 ,			,			,			yr4dwpn_pass1sc1024ac 1280K, 1280

build421 ,			,			,			yr4dwpn_pass1sc1280ac 960K, 768
build421 ,			,			,			yr4dwpn_pass1sc1280ac 1280K, 10
build421 ,			,			,			yr4dwpn_pass1sc1280ac 1600K, 1280

build421 ,			,			,			yr4dwpn_pass1sc1536ac 1536K, 10

build421 ,			,			,			yr4dwpn_pass1sc2048ac 2M, 10
ENDIF

	; The all-complex 11 levels variants (1536, 2048, 2560)

IFNDEF X86_64
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc128ac 192K, 1536
build421 ,			FMA3_64,		,			yr4dwpn_pass1sc128ac 256K, 11
build421 ,			CORE_32 + FMA3_64,	,			yr4dwpn_pass1sc128ac 320K, 2560

build421 ,			,			,			yr4dwpn_pass1sc256ac 384K, 1536
build421 ,			,			,			yr4dwpn_pass1sc256ac 512K, 11
build421 ,			,			,			yr4dwpn_pass1sc256ac 640K, 2560

build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc384ac 576K, 1536
build421 ,			CORE_64 + FMA3_64,	,			yr4dwpn_pass1sc384ac 768K, 11
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc384ac 960K, 2560

build421 ,			,			,			yr4dwpn_pass1sc512ac 768K, 1536
build421 ,			CORE_32,		,			yr4dwpn_pass1sc512ac 1M, 11
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc512ac 1280K, 2560

build421 ,			,			,			yr4dwpn_pass1sc640ac 960K, 1536
build421 ,			,			,			yr4dwpn_pass1sc640ac 1280K, 11
build421 ,			CORE    + FMA3_64,	,			yr4dwpn_pass1sc640ac 1600K, 2560

build421 ,			,			,			yr4dwpn_pass1sc768ac 1152K, 1536
build421 ,			,			,			yr4dwpn_pass1sc768ac 1536K, 11
build421 ,			,			,			yr4dwpn_pass1sc768ac 1920K, 2560

build421 ,			,			,			yr4dwpn_pass1sc1024ac 1536K, 1536
build421 ,			,			,			yr4dwpn_pass1sc1024ac 2M, 11
build421 ,			,			,			yr4dwpn_pass1sc1024ac 2560K, 2560

build421 ,			,			,			yr4dwpn_pass1sc1280ac 1920K, 1536
build421 ,			,			,			yr4dwpn_pass1sc1280ac 2560K, 11
build421 ,			,			,			yr4dwpn_pass1sc1280ac 3200K, 2560
ENDIF

;; Some small FFTs fit completely in a large L2 cache and thus do
;; not need to be prefetched.

IF TLB_PRIMING EQ 0

PREFETCHING = 0

;buildfor                  ,	ypass2gen 8
;buildfor                  ,	ypass2gen 768
;buildfor                  ,	ypass2gen 10
;buildfor                  ,	ypass2gen 1280

;build421 ,			,			,			yr4dwpn_pass1sc128 32K, 8
;build421 ,			,			,			yr4dwpn_pass1sc256 64K, 8
;build421 ,			,			,			yr4dwpn_pass1sc320 80K, 8
;build421 ,			,			,			yr4dwpn_pass1sc384 96K, 8
;build421 ,			,			,			yr4dwpn_pass1sc448 112K, 8
;build421 ,			,			,			yr4dwpn_pass1sc512 128K, 8
;build421 ,			,			,			yr4dwpn_pass1sc128 96K, 768
;build421 ,			,			,			yr4dwpn_pass1sc128 128K, 10

;build421 ,			,			,			yr4dwpn_pass1sc128ac 32K, 8
;build421 ,			,			,			yr4dwpn_pass1sc256ac 64K, 8
;build421 ,			,			,			yr4dwpn_pass1sc384ac 96K, 8
;build421 ,			,			,			yr4dwpn_pass1sc512ac 128K, 8
;build421 ,			,			,			yr4dwpn_pass1sc128ac 96K, 768
;build421 ,			,			,			yr4dwpn_pass1sc128ac 128K, 10

ENDIF

_TEXT	ENDS
END
