; Copyright 2001-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routine implement some common cleanup code for r4dwpn (r4delay with partial normalization) FFTs
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

PUBLIC	xgw_carries_wpn

_TEXT SEGMENT

;;*****************************************
;; Routine for finishing off a r4dwpn FFT
;;*****************************************

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

biglit_incr	EQU	PPTR [rsp+first_local]
grp_incr	EQU	PPTR [rsp+first_local+SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]

PROCF	xgw_carries_wpn
	int_prolog 2*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	edx, 4*XMM_GMD		; Compute grp increment
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	je	short iskip		; for rational FFTs
	sub	rax, rax		; Zero biglit_incr
	sub	rdx, rdx		; Zero grp_incr
iskip:	mov	biglit_incr, rax	; Save computed increments
	mov	grp_incr, rdx
	mov	rsi, carries		; Addr of the carries
	cmp	ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	xgw_carries_zpad
	xnorm012_wpn_part1
	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1:	cmp	B_IS_2, 0		; Is b = 2?
	jne	b2			; Yes, do simpler roundings
	xnorm012_wpn noexec		; Split carries for one cache line
	jmp	nb2			; Rejoin common code
b2:	xnorm012_wpn exec		; Split carries for one cache line
nb2:	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1
	add	rdx, grp_incr		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0			; Next carry row
	jmp	cdn			; Jump to common exit code

xgw_carries_zpad:
	xnorm012_wpn_zpad_part1
	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
zlp0:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
zlp1:	xnorm012_wpn_zpad		; Split carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	zlp1
	add	rdx, grp_incr		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	zlp0			; Next carry row

cdn:	int_epilog 2*SZPTR+8,0,0
xgw_carries_wpn ENDP


_TEXT	ENDS
END
