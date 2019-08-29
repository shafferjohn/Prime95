; Copyright 2011-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routines implement some common cleanup code for yr4dwpn (r4delay with partial normalization) FFTs
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE ymult.mac
INCLUDE ynormal.mac

_TEXT SEGMENT

;;*****************************************
;; Routine for finishing off a r4dwpn FFT
;;*****************************************

;; Multiple versions of the carry propagation routines

biglit_incr	EQU	PPTR [rsp+first_local]
grp_incr	EQU	PPTR [rsp+first_local+SZPTR]
saveptr1	EQU	PPTR [rsp+first_local+2*SZPTR]
saveptr2	EQU	PPTR [rsp+first_local+3*SZPTR]
saveptr3	EQU	PPTR [rsp+first_local+4*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+5*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+5*SZPTR+4]

;; Base 2, irrational, no zero padding
PROCFL	ygw_carries_wpn3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_part1 exec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1:	ynorm012_wpn exec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1
	bump	rdx, 2*YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpn3 ENDP


;; Base 2, rational, no zero padding
PROCFL	ygw_carries_wpnr3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_part1 exec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1:	ynorm012_wpn noexec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnr3 ENDP


;; Base 2, irrational, zero padding
PROCFL	ygw_carries_wpnzp3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_zpad_part1 exec, exec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0zp:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1zp:	ynorm012_wpn_zpad exec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1zp
	bump	rdx, YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0zp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnzp3 ENDP


;; Base 2, rational, zero padding
PROCFL	ygw_carries_wpnrzp3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_zpad_part1 noexec, exec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1zp:	ynorm012_wpn_zpad noexec, exec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1zp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnrzp3 ENDP


;; Not base 2, irrational, no zero padding
PROCFL	ygw_carries_wpnn3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_part1 noexec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0n:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1n:	ynorm012_wpn exec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1n
	bump	rdx, 2*YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0n			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnn3 ENDP


;; Not base 2, rational, no zero padding
PROCFL	ygw_carries_wpnnr3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_part1 noexec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1n:	ynorm012_wpn noexec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1n			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnnr3 ENDP


;; Not base 2, irrational, zero padding
PROCFL	ygw_carries_wpnnzp3
	int_prolog 5*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	biglit_incr, rax	; Save biglit increment

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_zpad_part1 exec, noexec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0nzp: mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1nzp: ynorm012_wpn_zpad exec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1nzp
	bump	rdx, YMM_GMD		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0nzp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnnzp3 ENDP


;; Not base 2, rational, zero padding
PROCFL	ygw_carries_wpnnrzp3
	int_prolog 5*SZPTR+8,0,0

	mov	rsi, carries		; Addr of the carries
	ynorm012_wpn_zpad_part1 noexec, noexec

	mov	rbp, DATA_ADDR		; Addr of the FFT data
	mov	eax, addcount1		; Load count of pass 1 blocks
	mov	loopcount2, eax		; Save count
rlp1nzp: ynorm012_wpn_zpad noexec, noexec, saveptr1, saveptr2, saveptr3 ; Propagate carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	sub	loopcount2, 1		; Test loop counter
	jnz	rlp1nzp			; Next carry row

	int_epilog 5*SZPTR+8,0,0
ygw_carries_wpnnrzp3 ENDP


_TEXT	ENDS
END
