; Copyright 1995-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the 8 level 2nd pass for FFTs.
;

	TITLE   setup

	.686
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE mult.mac
INCLUDE memory.mac
INCLUDE normal.mac

_TEXT SEGMENT

	flat_distances

;; Routines to do the normalization after a multiply

;; When doing zero-padded FFTs, the multiplied 7 words around the halfway point
;; must be subtracted from the bottom of the FFT.  This must be done before
;; normalization multiplies the FFT data by k.  This macro does that.

sub_7_words MACRO
	LOCAL	nozpad, zlp, done7
	lea	ebp, ZPAD6
	cmp	zpad_addr, ebp		;; Have we subtracted all 7 words?
	jg	short nozpad		;; Yes, skip this code
	mov	edi, zpad_addr		;; Addr of next zpad element to process
	mov	eax, cache_line_multiplier ;; Load loop counter
	add	eax, eax		;; Two values per cache line
	mov	ecx, esi		;; Copy source ptr (we preserve esi)
zlp:	fld	QWORD PTR [ecx]		;; Load FFT word
	fld	QWORD PTR [edi]		;; Load ZPAD data
	fmul	NORM012_FF		;; Scale by FFTLEN/2
	fsub	st(1), st
	faddp	st(2), st		;; Adjust sumout
	fstp	QWORD PTR [ecx]		;; Store FFT word
	lea	ecx, [ecx+dist1]	;; Bump pointers
	lea	edi, [edi+8]
	cmp	edi, ebp		;; Only subtract 7 zpad words
	jg	short done7
	dec	eax			;; Iterate 2*clm times
	jnz	short zlp		;; Loop if necessary
done7:	mov	zpad_addr, edi
nozpad:
	ENDM


; Macro to loop through all the FFT values and apply the proper normalization
; routine.  Used whenever we are using an irrational-base FFT.

saved_ptr	EQU	DPTR [rsp+first_local]
loopcount2	EQU	DPTR [rsp+first_local+4]
loopcount3	EQU	DPTR [rsp+first_local+8]
saved_edx	EQU	DPTR [rsp+first_local+12]
saved_esi	EQU	DPTR [rsp+first_local+16]

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp0, ilp1, ilpdn, done
	PROCFLP	lab
	int_prolog 20,0,0,rcx,rbp
	mov	saved_edx, edx		;; Save registers for top_carry_adjust
	mov	saved_esi, esi
no zero	mov	zero_fft, 0		;; Set flag saying not zero upper half
zero	mov	zero_fft, 1		;; Set flag saying zero upper half
echk	fld	MAXERR			;; Load maximum error
	fld	SUMOUT			;; Load SUMOUT
no zero	cmp	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	fld	QWORD PTR [esi][edi]	;; Get the value
no zero	fadd	ADDIN_VALUE		;; Add in the requested value
no zero	fstp	QWORD PTR [esi][edi]	;; Save the new value
no zero	fsub	ADDIN_VALUE		;; Do not include addin in sumout
noadd:	mov	edx, norm_grp_mults	;; Addr of the group multipliers
	mov	ebp, carries		;; Addr of the carries
	mov	edi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	eax, eax		;; Clear big/lit flags
ilp0:	mov	ecx, cache_line_multiplier ;; Load inner loop counter
	mov	ebx, norm_ptr2		;; Load column multipliers ptr
	fld	QWORD PTR [ebp+0*8]	;; Load carry1
	fadd	BIGVAL			;; c1 = c1 + rounding constant
	fld	QWORD PTR [ebp+1*8]	;; Load carry2
	fadd	BIGVAL			;; c2 = c2 + rounding constant
	mov	saved_ptr, esi		;; Save FFT data ptr
ilp1:	norm_2d ttp, zero, echk, const	;; Normalize 4 vals (2 grps, 2 cols)
	lea	esi, [esi+4*8]		;; Next FFT source
ttp	lea	ebx, [ebx+2*16]		;; Next column multipliers
ttp	lea	edi, [edi+2]		;; Next big/little flags
	dec	ecx			;; Test clm loop counter
	jnz	ilp1			;; Loop til done
	mov	esi, saved_ptr		;; Restore FFT data ptr
	fsub	BIGVAL			;; c2 = c2 - rounding constant
	fstp	QWORD PTR [ebp+1*8]	;; Save carry2
	fsub	BIGVAL			;; c1 = c1 - rounding constant
	fstp	QWORD PTR [ebp+0*8]	;; Save carry1
	add	esi, normblkdst		;; Next source pointer
	lea	ebp, [ebp+2*8]		;; Next set of carries
ttp	lea	edx, [edx+2*16]		;; Next set of 2 group multipliers
	dec	loopcount2		;; Test outer loop counter
	jz	short ilpdn		;; Iterate
	add	loopcount3, 80000000h/16;; 32 iterations
	jnc	ilp0
	add	esi, normblkdst8	;; Add 64 pad bytes every 32 clmblkdsts
	jmp	ilp0			;; Iterate
ilpdn:	fstp	SUMOUT			;; Save SUMOUT
echk	fstp	MAXERR			;; Save maximum error
ttp	mov	norm_ptr1, edi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, ebx		;; Save column multipliers ptr

	;; Handle adjusting the carry out of the topmost FFT word
	mov	edx, saved_edx		;; Restore edx (pass1 loop counter)
	mov	esi, saved_esi		;; Restore esi (FFT data pointer)
	cmp	edx, 65536+256		;; Check for last iteration
	jne	done			;; Top carry may require adjusting
	norm_top_carry_2d		;; Adjust carry for k > 1

done:	int_epilog 20,0,0,rcx,rbp
	ENDPP	lab
	ENDM

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp0, ilp1, ilpdn
	PROCFLP	lab
	int_prolog 12,0,0,rcx,rdx,rsi,rbp
no const mov	const_fft, 0		;; Set flag saying not mul-by-const
const	mov	const_fft, 1		;; Set flag saying mul-by-const
echk	fld	MAXERR			;; Load maximum error
	fld	SUMOUT			;; Load SUMOUT
	sub_7_words
	mov	edx, norm_grp_mults	;; Addr of the group multipliers
	mov	ebp, carries		;; Addr of the carries
	mov	edi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	eax, eax		;; Clear big/lit flags
ilp0:	mov	ecx, cache_line_multiplier ;; Load inner loop counter
	mov	ebx, norm_ptr2		;; Load column multipliers ptr
	fld	QWORD PTR [ebp+0*8]	;; Load carry1
	fadd	BIGVAL			;; c1 = c1 + rounding constant
	fld	QWORD PTR [ebp+1*8]	;; Load carry2
	mov	saved_ptr, esi		;; Save FFT data ptr
ilp1:	norm_2d_zpad ttp, echk, const	;; Normalize 4 vals (2 grps, 2 cols)
	lea	esi, [esi+4*8]		;; Next FFT source
ttp	lea	ebx, [ebx+2*16]		;; Next column multipliers
ttp	lea	edi, [edi+2]		;; Next big/little flags
	dec	ecx			;; Test clm loop counter
	jnz	ilp1			;; Loop til done
	mov	esi, saved_ptr		;; Restore FFT data ptr
	fstp	QWORD PTR [ebp+1*8]	;; Save carry2
	fsub	BIGVAL			;; c1 = c1 - rounding constant
	fstp	QWORD PTR [ebp+0*8]	;; Save carry1
	add	esi, normblkdst		;; Next source pointer
	lea	ebp, [ebp+2*8]		;; Next set of carries
ttp	lea	edx, [edx+2*16]		;; Next set of 2 group multipliers
	dec	loopcount2		;; Test outer loop counter
	jz	short ilpdn		;; Iterate
	add	loopcount3, 80000000h/16;; 32 iterations
	jnc	ilp0
	add	esi, normblkdst8	;; Add 64 pad bytes every 32 clmblkdsts
	jmp	ilp0			;; Iterate
ilpdn:	fstp	SUMOUT			;; Save SUMOUT
echk	fstp	MAXERR			;; Save maximum error
ttp	mov	norm_ptr1, edi		;; Save big/little flags array ptr
ttp	mov	norm_ptr2, ebx		;; Save column multipliers ptr
	int_epilog 12,0,0,rcx,rdx,rsi,rbp
	ENDPP	lab
	ENDM

; The 16 different normalization routines

	inorm	r2, noexec, noexec, noexec, noexec
	inorm	r2e, noexec, noexec, exec, noexec
	inorm	r2c, noexec, noexec, noexec, exec
	inorm	r2ec, noexec, noexec, exec, exec
	inorm	r2z, noexec, exec, noexec, noexec
	inorm	r2ze, noexec, exec, exec, noexec
	inorm	i2, exec, noexec, noexec, noexec
	inorm	i2e, exec, noexec, exec, noexec
	inorm	i2c, exec, noexec, noexec, exec
	inorm	i2ec, exec, noexec, exec, exec
	inorm	i2z, exec, exec, noexec, noexec
	inorm	i2ze, exec, exec, exec, noexec

	zpnorm	r2zp, noexec, noexec, noexec
	zpnorm	r2zpe, noexec, exec, noexec
	zpnorm	r2zpc, noexec, noexec, exec
	zpnorm	r2zpec, noexec, exec, exec
	zpnorm	i2zp, exec, noexec, noexec
	zpnorm	i2zpe, exec, exec, noexec
	zpnorm	i2zpc, exec, noexec, exec
	zpnorm	i2zpec, exec, exec, exec

_TEXT	ENDS
END
