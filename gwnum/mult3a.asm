; Copyright 1995-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the 10 level 2nd pass for FFTs.
;

	TITLE   setup

	.686
	.XMM
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE	lucas.mac
INCLUDE mult.mac
INCLUDE pass2.mac
INCLUDE memory.mac

_TEXT SEGMENT

	flat_distances

;;*****************************************************
;; Macros to do the code common to all pass 2 routines
;;*****************************************************

;; Entry point for real FFTs: do one real and count1 complex blocks, also
;; entry point for all-complex FFTs: do count1 complex blocks

pass2_entry MACRO complex_start
	int_prolog 0,0,0
	cmp	ALL_COMPLEX_FFT, 1	;; Test if there is a real-data block
	je	complex_start		;; Jump to process all-complex blocks
	ENDM

;; Exit code for pass2 routines

pass2_exit MACRO complex_start
	mov	esi, DESTARG		;; Restore source pointer
	sub	ebx, ebx		;; Clear DIST_TO_FFTSRCARG
	mov	DIST_TO_FFTSRCARG, ebx
	int_epilog 0,0,0
	ENDM



;; Routines to do the last 8 levels in a two-pass FFT

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels
;; 2^8 complex values = 2^9 doubles = 4KB.

blkdst	=	(4096+64+64)

;; Do the last 8 levels of a two pass FFT

PROCFP	pass2_8_levels
	pass2_entry p8cs
	pass2_eight_levels_real
p8cs:	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p8lp:	pass2_eight_levels_complex
	dec	ecx
	jnz	p8lp
	pass2_exit
ENDPP	pass2_8_levels


;; Routines to do the last 10 levels in a two-pass FFT

;; Distance between two pass 2 data blocks.  Pass 2 does 10 FFT levels
;; 2^10 complex values = 2^11 doubles = 16KB.

blkdst	=	(4*(4096+64)+64)

;; Do the last 10 levels of a two pass FFT

PROCFP	pass2_10_levels
	pass2_entry p10cs
	pass2_ten_levels_real
p10cs:	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p10lp:	pass2_ten_levels_complex
	dec	ecx
	jnz	p10lp
	pass2_exit
ENDPP	pass2_10_levels


;; Routines to do the last 12 levels in a two-pass FFT

;; Distance between two pass 2 data blocks.  Pass 2 does 12 FFT levels
;; 2^12 complex values = 2^13 doubles = 64KB.

blkdst	=	(16*(4096+64)+64)

;; Do the last 12 levels of a two pass FFT

PROCFP	pass2_12_levels
	pass2_entry p12cs
	pass2_twelve_levels_real
p12cs:	mov	ecx, count1		; Number of complex iterations
	mov	edx, pass2_premults	; Address of the group multipliers
p12lp:	pass2_twelve_levels_complex
	dec	ecx
	jnz	p12lp
	pass2_exit
ENDPP	pass2_12_levels


_TEXT	ENDS
END
