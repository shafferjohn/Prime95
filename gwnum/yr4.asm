; Copyright 2011-2020 Mersenne Research, Inc.  All rights reserved
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

IFNDEF X86_64
EXTRN	yreal_fft_final1a:PROC
ENDIF
EXTRN	yreal_fft_final1:PROC
EXTRN	yreal_fft_final2:PROC
EXTRN	yreal_fft_final3:PROC
EXTRN	yreal_fft_final4:PROC
EXTRN	yreal_fft_final4a:PROC
EXTRN	yreal_fft_final5:PROC
EXTRN	yreal_fft_final6:PROC
IFNDEF X86_64
EXTRN	yreal_square1a:PROC
ENDIF
EXTRN	yreal_square1:PROC
EXTRN	yreal_square2:PROC
EXTRN	yreal_square3:PROC
EXTRN	yreal_square4:PROC
EXTRN	yreal_square4a:PROC
EXTRN	yreal_square5:PROC
EXTRN	yreal_square6:PROC
IFNDEF X86_64
EXTRN	yreal_mult1a:PROC
ENDIF
EXTRN	yreal_mult1:PROC
EXTRN	yreal_mult2:PROC
EXTRN	yreal_mult3:PROC
EXTRN	yreal_mult4:PROC
EXTRN	yreal_mult4a:PROC
EXTRN	yreal_mult5:PROC
EXTRN	yreal_mult6:PROC
IFNDEF X86_64
EXTRN	yreal_mulf1a:PROC
ENDIF
EXTRN	yreal_mulf1:PROC
EXTRN	yreal_mulf2:PROC
EXTRN	yreal_mulf3:PROC
EXTRN	yreal_mulf4:PROC
EXTRN	yreal_mulf4a:PROC
EXTRN	yreal_mulf5:PROC
EXTRN	yreal_mulf6:PROC
IFDEF X86_64
EXTRNP	ycomplex_mult_opcode1
EXTRNP	ycomplex_mult_opcode2
EXTRNP	ycomplex_mult_opcode3
EXTRNP	ycomplex_mult_opcode4
EXTRNP	ycomplex_mult_opcode4a
EXTRNP	ycomplex_mult_opcode5
EXTRNP	ycomplex_mult_opcode6
EXTRNP	ycomplex_mulf_opcode1
EXTRNP	ycomplex_mulf_opcode2
EXTRNP	ycomplex_mulf_opcode3
EXTRNP	ycomplex_mulf_opcode4
EXTRNP	ycomplex_mulf_opcode4a
EXTRNP	ycomplex_mulf_opcode5
EXTRNP	ycomplex_mulf_opcode6
ENDIF
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
