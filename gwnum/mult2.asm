; Copyright 1995-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code handles FFTs that use two passes with 8 levels done on the
; second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

	.686
	.XMM
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE	lucas.mac
INCLUDE pfa.mac
INCLUDE mult.mac
INCLUDE pass1.mac
INCLUDE fft2.mac
INCLUDE memory.mac
INCLUDE normal.mac

IFNDEF PFETCH
PUBLIC gw_finish_fft
PUBLIC gw_carries
PUBLIC gw_finish_mult
ELSE
EXTRN gw_finish_fft:PROC
EXTRN gw_carries:PROC
EXTRN gw_finish_mult:PROC
ENDIF

EXTRNP	pass2_8_levels

_TEXT SEGMENT

	flat_distances

;; Distance between two pass 2 data blocks.  Pass 2 does 8 FFT levels.
;; 2^8 complex values = 2^9 doubles = 4KB.

blkdst	EQU	(4096+64+64)

;; All the FFT routines for each FFT length.  We don't implement prefetching
;; versions for some of the smaller FFTs.

IFNDEF PFETCH
	fft	5120
	fft	6144
	fft	6144p
	fft	7168
	fft	8192
	fft	8192p
	fft	10K
	fft	12K
	fft	12Kp
	fft	14K
	fft	16K
	fft	16Kp
	fft	20K
	fft	24K
	fft	24Kp
	fft	28K
	fft	32K
	fft	32Kp
ENDIF

;INCLUDE pass1scr.mac

	fft	40K
	fft	48K
	fft	48Kp
	fft	56K
	fft	64K
	fft	64Kp
	fft	80K
;	fft	96K
;	fft	96Kp
;	fft	112K
;	fft	128K
;	fft	128Kp
;	fft	160K
;	fft	192K
;	fft	192Kp
;	fft	224K
;	fft	256K
;	fft	256Kp

IFNDEF PFETCH

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

loopcount1	EQU	DPTR [rsp+first_local]

PROCF	gw_carries
	int_prolog 4,0,0
	cmp	ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	gw_carries_zpad
	mov	esi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute end of carries addr
	shl	ebx, 4			; Two 8-byte doubles per section
	add	ebx, esi
	norm012_2d_part1
	mov	ebp, DESTARG		; Addr of the FFT data
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	mov	ecx, addcount1		; Load section count
	sub	eax, eax		; Clear big/little flag
ilp1:	mov	ebx, norm_col_mults	; Addr of the column multipliers
	norm012_2d			; Split carries for one cache line
	mov	ebx, cache_line_multiplier ; Cache lines in each pass1 loop
	lea	esi, [esi+2*8]		; Next carries pointer
	add	ebp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short iskip		; for rational FFTs
	lea	edx, [edx+2*16]		; Next group multiplier
	lea	edi, [edi+ebx*2]	; Next big/little flags pointer
iskip:	dec	ecx			; Test loop counter
	jnz	ilp1			; Next carry row in section
	fcompp				; Pop last 2 carries
	jmp	cdn			; Jump to common exit code

gw_carries_zpad:
	mov	esi, carries		; Addr of the carries
	mov	ebp, addcount1		; Compute end of carries addr
	shl	ebp, 4
	add	ebp, esi
	mov	esi, DESTARG		; Addr of the FFT data
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	mov	ebx, norm_col_mults	; Addr of the group multipliers
	sub	eax, eax		; Clear big/little flag
	cmp	const_fft, 0		; Call correct part1 macro
	je	c2a			; Jump if not const
	norm012_2d_zpad_part1 exec
	jmp	c2b
c2a:	norm012_2d_zpad_part1 noexec
c2b:	mov	ebp, carries		; Addr of the carries
	mov	edx, norm_grp_mults	; Addr of the group multipliers
	mov	ecx, addcount1		; Load counter
	mov	loopcount1, ecx		; Save counter
zlp1:	mov	ebx, norm_col_mults	; Addr of the column multipliers
	cmp	const_fft, 0		; Call correct zpad macro
	je	c2c			; Jump if not const
	norm012_2d_zpad exec		; Split carries for one cache line
	jmp	c2d
c2c:	norm012_2d_zpad noexec		; Split carries for one cache line
c2d:	mov	ebx, cache_line_multiplier ; Cache lines in each pass1 loop
	lea	ebp, [ebp+2*8]		; Next carries pointer
	add	esi, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short zskip		; for rational FFTs
	lea	edx, [edx+2*16]		; Next group multiplier
	lea	edi, [edi+ebx*2]	; Next big/little flags pointer
zskip:	dec	loopcount1		; Test loop counter
	jnz	zlp1			; Next carry row in section
	fcompp				; Pop last 2 carries

cdn:	int_epilog 4,0,0
gw_carries ENDP


; Common code to finish off the two-pass FFTs.  The Windows 64-bit ABI
; frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.
; Of course, this code won't be linked into a 64-bit Windows executable,
; but we include the dummy prolog to be consistent.

PROCF	__common_2pass_fft_exit_code

	;; Create a dummy prolog
	ad_prolog 0,0,rbx,rbp,rsi,rdi

; Common code to finish the FFT by restoring state and returning.

gw_finish_fft:
	mov	DWORD PTR [esi-28], 3	; Set has-been-FFTed flags
	fft_1_ret

; Common code to finish up multiplies

gw_finish_mult:

; Set FFT-started flag

	mov	esi, DESTARG		; Addr of FFT data
	mov	al, POSTFFT		; Set FFT started flag
	mov	BYTE PTR [esi-28], al

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	fld	SUMOUT
	fmul	ttmp_ff_inv
	fstp	QWORD PTR [esi-24]	; Save sum of FFT outputs

; Return

	ad_epilog 0,0,rbx,rbp,rsi,rdi

__common_2pass_fft_exit_code ENDP

ENDIF

_TEXT	ENDS
END
