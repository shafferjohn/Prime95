; Copyright 2001-2014 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes between than 5K and 128K doubles are supported.
; This code does two passes, 8 levels on the second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac

PUBLIC	xgw_carries

_TEXT SEGMENT

;;*************************************************************************
;; Routine for auxiliary threads to call to start processing pass 1 blocks
;;*************************************************************************

; pass1_aux_entry_point ()
; Entry point for auxiliary threads to do process blocks in pass 1
; Windows 32-bit and Linux 32-bit
;	Parameter asm_data = [esp+4]
; Windows 64-bit
;	Parameter asm_data = rcx
; Linux 64-bit
;	Parameter asm_data = rdi

PROCFL pass1_aux_entry_point
	ad_prolog 0,1,rbx,rbp,rsi,rdi,r12,r13,r14,r15,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

;; Jump to common pass 1 code.  We must jump rather than call because the
;; pass 1 code operates with a push_amt of zero.

	mov	rax, THREAD_WORK_ROUTINE ; Pass 1 entry point
	jmp	rax			; Go process data blocks
pass1_aux_entry_point ENDP

;;*************************************************************************
;; Routine for auxiliary threads to call to start processing pass 2 blocks
;;*************************************************************************

; pass2_aux_entry_point ()
; Entry point for auxiliary threads to process blocks in pass 2
; Windows 32-bit and Linux 32-bit
;	Parameter asm_data = [esp+4]
; Windows 64-bit
;	Parameter asm_data = rcx
; Linux 64-bit
;	Parameter asm_data = rdi

PROCFL	pass2_aux_entry_point
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

;; Call intermediate routine to mimic main thread's prologs.  We must do this
;; so that common epilog code performs properly.

	call	internal_pass2_aux_entry_point

;; Return to C code

	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15
pass2_aux_entry_point ENDP

; Intermediate routine to call to set up an identical prolog as the
; main thread's processing of pass 2 data blocks.

PROCF	internal_pass2_aux_entry_point
	int_prolog 0,0,1

;; Jump to common code handling pass 2 data blocks.

	mov	rax, THREAD_WORK_ROUTINE ; Pass 2 entry point
	jmp	rax			; Go process data blocks

internal_pass2_aux_entry_point ENDP


;;*****************************************
;; Routine for finishing off a two-pass FFT
;;*****************************************

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

loopcount1	EQU	DPTR [rsp+first_local]

PROCF	xgw_carries
	int_prolog 4,0,0
	mov	rsi, carries		; Addr of the carries
	cmp	ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	xgw_carries_zpad
	xnorm012_2d_part1
	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	rbx, norm_ptr2		; Addr of the column multipliers
	mov	eax, addcount1		; Load count of cache lines in carries array
	mov	loopcount1, eax		; Save for later
	cmp	B_IS_2, 0		; Is b = 2?
	jne	ilp1			; Yes, do simpler roundings
nb2ilp1:xnorm012_2d noexec		; Split carries for one cache line
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short nb2iskip		; for rational FFTs
	bump	rdx, 128		; Next group multiplier
	lea	rdi, [rdi+rax*4]	; Next big/little flags pointer
nb2iskip:sub	loopcount1, 1		; Test loop counter
	jnz	nb2ilp1			; Next carry row
	jmp	cdn			; Done
ilp1:	xnorm012_2d exec		; Split carries for one cache line
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short iskip		; for rational FFTs
	bump	rdx, 128		; Next group multiplier
	lea	rdi, [rdi+rax*4]	; Next big/little flags pointer
iskip:	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1			; Next carry row
	jmp	cdn			; Done

xgw_carries_zpad:
	xnorm012_2d_zpad_part1
	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	rbx, norm_ptr2		; Addr of the column multipliers
	mov	eax, addcount1		; Load count of cache lines in carries array
	mov	loopcount1, eax		; Save for later
zlp1:	xnorm012_2d_zpad		; Split carries for one cache line
c2d:	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short zskip		; for rational FFTs
	bump	rdx, 128		; Next group multiplier
	lea	rdi, [rdi+rax*4]	; Next big/little flags pointer
zskip:	sub	loopcount1, 1		; Test loop counter
	jnz	zlp1			; Next carry row in section
cdn:	int_epilog 4,0,0
xgw_carries ENDP


_TEXT	ENDS
END
