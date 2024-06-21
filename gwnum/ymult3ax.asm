; Copyright 2011-2023 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu

;
; Additional routines used with yr4dwpn (r4delay with partial normalization) FFTs,
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

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwyaddq3
	ad_prolog 0,0,rbx,rsi,rdi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	ebx, 4			; Four pass 2 blks in one pass 1 block
qadd0:	mov	eax, normval4		; Load count of 4KB chunks in a block
qadd1:	mov	edi, normval1		; Load count of clms in 4KB
	imul	edi, cache_line_multiplier ; Compute cache lines in 4KB chunk
	shr	edi, 2
qaddlp:	vmovapd	ymm0, [rdx]		; Load second number
	vaddpd	ymm0, ymm0, [rcx]	; Add in first number
	vmovapd	ymm1, [rdx+32]		; Load second number
	vaddpd	ymm1, ymm1, [rcx+32]	; Add in first number
	ystore	[rsi], ymm0		; Save result
	ystore	[rsi+32], ymm1		; Save result
	bump	rcx, 64			; Next source
	bump	rdx, 64			; Next source
	bump	rsi, 64			; Next dest
	dec	rdi			; Test for end of 4KB chunk
	jnz	short qaddlp		; Loop if necessary
	add	rcx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rsi, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	dec	rax			; Check middle loop counter
	jnz	short qadd1		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	qadd0			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwyaddq3 ENDP

;;
;; Add two numbers with carry propagation (eight different versions)
;;

dist_to_dest	EQU	PPTR [rsp+first_local+0*SZPTR]
loopcount3	EQU	DPTR [rsp+first_local+1*SZPTR]
loopcount4	EQU	DPTR [rsp+first_local+1*SZPTR+4]
loopcount5	EQU	DPTR [rsp+first_local+1*SZPTR+8]

	; Base 2, irrational, not zero-padded
PROCFL	gwyadd3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rax, DESTARG			; Address of destination
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2:	ynorm_op_wpn vaddpd, exec, exec, dist_to_dest	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyadd3 ENDP


	; Base 2, rational, not zero-padded
PROCFL	gwyaddr3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rdi, DESTARG			; Address of destination
	sub	rdi, rsi			; Calculate distance from first number to destination
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0r:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1r:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2r:	ynorm_op_wpn vaddpd, noexec, exec, rdi	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2r
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2r
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1r				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0r

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddr3 ENDP


	; Not base 2, irrational, not zero-padded
PROCFL	gwyaddn3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rax, DESTARG			; Address of destination
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0n:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1n:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2n:	ynorm_op_wpn vaddpd, exec, noexec, dist_to_dest	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2n
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2n
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1n				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0n

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddn3 ENDP


	; Not base 2, rational, not zero-padded
PROCFL	gwyaddnr3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rdi, DESTARG			; Address of destination
	sub	rdi, rsi			; Calculate distance from first number to destination
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0nr:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1nr:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2nr:	ynorm_op_wpn vaddpd, noexec, noexec, rdi ; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2nr
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2nr
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1nr				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0nr

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddnr3 ENDP


	; Base 2, irrational, zero-padded
PROCFL	gwyaddzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0zp:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1zp:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2zp:	ynorm_op_wpn_zpad vaddpd, exec, exec	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2zp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0zp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddzp3 ENDP


	; Base 2, rational, zero-padded
PROCFL	gwyaddrzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0rzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1rzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2rzp: ynorm_op_wpn_zpad vaddpd, noexec, exec	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2rzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2rzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1rzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0rzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddrzp3 ENDP


	; Not base 2, irrational, zero-padded
PROCFL	gwyaddnzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0nzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1nzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2nzp: ynorm_op_wpn_zpad vaddpd, exec, noexec	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2nzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2nzp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1nzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0nzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddnzp3 ENDP


	; Not base 2, rational, zero-padded
PROCFL	gwyaddnrzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
add0nrzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
add1nrzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
add2nrzp: ynorm_op_wpn_zpad vaddpd, noexec, noexec ; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	add2nrzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	add2nrzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	add1nrzp			; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	add0nrzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwyaddnrzp3 ENDP


;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwysubq3
	ad_prolog 0,0,rbx,rsi,rdi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	ebx, 4			; Four pass 2 blks in one pass 1 block
qsub0:	mov	eax, normval4		; Load count of 4KB chunks in a block
qsub1:	mov	edi, normval1		; Load count of clms in 4KB
	imul	edi, cache_line_multiplier ; Compute cache lines in 4KB chunk
	shr	edi, 2
qsublp:	vmovapd	ymm0, [rdx]		; Load second number
	vsubpd	ymm0, ymm0, [rcx]	; Subtract first number
	vmovapd	ymm1, [rdx+32]		; Load second number
	vsubpd	ymm1, ymm1, [rcx+32]	; Subtract first number
	ystore	[rsi], ymm0		; Save result
	ystore	[rsi+32], ymm1		; Save result
	bump	rcx, 64			; Next source
	bump	rdx, 64			; Next source
	bump	rsi, 64			; Next dest
	dec	rdi			; Test for end of 4KB chunk
	jnz	short qsublp		; Loop if necessary
	add	rcx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rsi, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	dec	rax			; Check middle loop counter
	jnz	short qsub1		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	qsub0			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwysubq3 ENDP

;;
;; Subtract two numbers with carry propagation (eight different versions)
;;

dist_to_dest	EQU	PPTR [rsp+first_local+0*SZPTR]
loopcount3	EQU	DPTR [rsp+first_local+1*SZPTR]
loopcount4	EQU	DPTR [rsp+first_local+1*SZPTR+4]
loopcount5	EQU	DPTR [rsp+first_local+1*SZPTR+8]

	; Base 2, irrational, not zero-padded
PROCFL	gwysub3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rax, DESTARG			; Address of destination
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2:	ynorm_op_wpn vsubpd, exec, exec, dist_to_dest	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysub3 ENDP


	; Base 2, rational, not zero-padded
PROCFL	gwysubr3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rdi, DESTARG			; Address of destination
	sub	rdi, rsi			; Calculate distance from first number to destination
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0r:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1r:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2r:	ynorm_op_wpn vsubpd, noexec, exec, rdi	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2r
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2r
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1r				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0r

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubr3 ENDP


	; Not base 2, irrational, not zero-padded
PROCFL	gwysubn3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rax, DESTARG			; Address of destination
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0n:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1n:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2n:	ynorm_op_wpn vsubpd, exec, noexec, dist_to_dest	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2n
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2n
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1n				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0n

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubn3 ENDP


	; Not base 2, rational, not zero-padded
PROCFL	gwysubnr3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rdi, DESTARG			; Address of destination
	sub	rdi, rsi			; Calculate distance from first number to destination
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0nr:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1nr:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2nr:	ynorm_op_wpn vsubpd, noexec, noexec, rdi ; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2nr
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2nr
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1nr				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0nr

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubnr3 ENDP


	; Base 2, irrational, zero-padded
PROCFL	gwysubzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0zp:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1zp:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2zp:	ynorm_op_wpn_zpad vsubpd, exec, exec	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2zp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0zp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubzp3 ENDP


	; Base 2, rational, zero-padded
PROCFL	gwysubrzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0rzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1rzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2rzp: ynorm_op_wpn_zpad vsubpd, noexec, exec	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2rzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2rzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1rzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0rzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubrzp3 ENDP


	; Not base 2, irrational, zero-padded
PROCFL	gwysubnzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0nzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1nzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2nzp: ynorm_op_wpn_zpad vsubpd, exec, noexec	; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2nzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2nzp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1nzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0nzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubnzp3 ENDP


	; Not base 2, rational, zero-padded
PROCFL	gwysubnrzp3
	ad_prolog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
sub0nrzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
sub1nrzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
sub2nrzp: ynorm_op_wpn_zpad vsubpd, noexec, noexec ; Add and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	sub2nrzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	sub2nrzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	sub1nrzp			; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	sub0nrzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 1*SZPTR+12,0,rbx,rbp,rsi,rdi
gwysubnrzp3 ENDP


;;
;; Add and subtract two numbers without carry propagation
;;

PROCFL	gwyaddsubq3
	ad_prolog 0,0,rbx,rbp,rsi,rdi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	ebx, 4			; Four pass 2 blks in one pass 1 block
qaddsub0:mov	eax, normval4		; Load count of 4KB chunks in a block
qaddsub1:mov	edi, normval1		; Load count of clms in 4KB
	imul	edi, cache_line_multiplier ; Compute cache lines in 4KB chunk
	shr	edi, 2
qaddsublp:
	vmovapd	ymm0, [rcx]		; Load first number
	vsubpd	ymm1, ymm0, [rdx]	; Subtract out second number
	vaddpd	ymm0, ymm0, [rdx]	; Add in second number
	vmovapd	ymm2, [rcx+32]		; Load first number
	vsubpd	ymm3, ymm2, [rdx+32]	; Subtract out second number
	vaddpd	ymm2, ymm2, [rdx+32]	; Add in second number
	ystore	[rsi], ymm0		; Save result
	ystore	[rbp], ymm1		; Save result
	ystore	[rsi+32], ymm2		; Save result
	ystore	[rbp+32], ymm3		; Save result
	bump	rcx, 64			; Next source
	bump	rdx, 64			; Next source
	bump	rsi, 64			; Next dest
	bump	rbp, 64			; Next dest
	dec	rdi			; Test for end of 4KB chunk
	jnz	short qaddsublp		; Loop if necessary
	add	rcx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rsi, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rbp, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	dec	rax			; Check middle loop counter
	jnz	qaddsub1		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	qaddsub0		; Loop if necessary
	ad_epilog 0,0,rbx,rbp,rsi,rdi
gwyaddsubq3 ENDP

;;
;; Add and subtract two numbers with carry propagation (eight different versions)
;;

dist_to_dest1	EQU	PPTR [rsp+first_local+0*SZPTR]
dist_to_dest2	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount3	EQU	DPTR [rsp+first_local+2*SZPTR]
loopcount4	EQU	DPTR [rsp+first_local+2*SZPTR+4]
loopcount5	EQU	DPTR [rsp+first_local+2*SZPTR+8]

	; Base 2, irrational, not zero-padded
PROCFL	gwyaddsub3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rax, DESTARG			; Address of destination
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest1, rax		; Save distance to dest
	mov	rax, DEST2ARG		  	; Address of destination #2
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest2, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vmovapd	ymm2, YMM_BIGVAL		; Init 4 carry registers
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2:	ynorm_addsub_wpn exec, exec, dist_to_dest1, dist_to_dest2 ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsub3 ENDP


	; Base 2, rational, not zero-padded
PROCFL	gwyaddsubr3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rdi, DESTARG			; Address of destination
	sub	rdi, rsi			; Calculate distance from first number to destination
	mov	rbp, DEST2ARG		  	; Address of destination #2
	sub	rbp, rsi			; Calculate distance from first number to destination
	vmovapd	ymm2, YMM_BIGVAL		; Init 4 carry registers
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0r:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1r:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2r:	ynorm_addsub_wpn noexec, exec, rdi, rbp	; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2r
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2r
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1r				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0r

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubr3 ENDP


	; Not base 2, irrational, not zero-padded
PROCFL	gwyaddsubn3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rax, DESTARG			; Address of destination
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest1, rax		; Save distance to dest
	mov	rax, DEST2ARG		  	; Address of destination #2
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest2, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vxorpd	ymm2, ymm2, ymm2		; Init 4 carry registers
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0n:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1n:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2n:	ynorm_addsub_wpn exec, noexec, dist_to_dest1, dist_to_dest2 ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2n
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2n
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1n				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0n

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubn3 ENDP


	; Not base 2, rational, not zero-padded
PROCFL	gwyaddsubnr3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rdi, DESTARG			; Address of destination
	sub	rdi, rsi			; Calculate distance from first number to destination
	mov	rbp, DEST2ARG		  	; Address of destination #2
	sub	rbp, rsi			; Calculate distance from first number to destination
	vxorpd	ymm2, ymm2, ymm2		; Init 4 carry registers
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0nr:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1nr:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2nr:	ynorm_addsub_wpn noexec, noexec, rdi, rbp ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2nr
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2nr
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1nr				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0nr

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubnr3 ENDP


	; Base 2, irrational, zero-padded
PROCFL	gwyaddsubzp3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rax, DEST2ARG		  	; Address of destination #2
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest2, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vmovapd	ymm2, YMM_BIGVAL		; Init 4 carry registers
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0zp:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1zp:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2zp:	ynorm_addsub_wpn_zpad exec, exec, dist_to_dest2 ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2zp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0zp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2_zpad		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubzp3 ENDP


	; Base 2, rational, zero-padded
PROCFL	gwyaddsubrzp3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rdi, DEST2ARG		  	; Address of destination #2
	sub	rdi, rsi			; Calculate distance from first number to destination
	vmovapd	ymm2, YMM_BIGVAL		; Init 4 carry registers
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0rzp:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1rzp:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2rzp:	ynorm_addsub_wpn_zpad noexec, exec, rdi ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2rzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2rzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1rzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0rzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2_zpad		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubrzp3 ENDP


	; Not base 2, irrational, zero-padded
PROCFL	gwyaddsubnzp3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rax, DEST2ARG		  	; Address of destination #2
	sub	rax, rsi			; Calculate distance from first number to destination
	mov	dist_to_dest2, rax		; Save distance to dest
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vxorpd	ymm2, ymm2, ymm2		; Init 4 carry registers
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0nzp:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1nzp:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2nzp:	ynorm_addsub_wpn_zpad exec, noexec, dist_to_dest2 ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2nzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2nzp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1nzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0nzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2_zpad		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubnzp3 ENDP


	; Not base 2, rational, zero-padded
PROCFL	gwyaddsubnrzp3
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, SRCARG			; Address of first number
	mov	rdx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	sub	rbx, rsi			; Calculate distance from first number to destination
	mov	rdi, DEST2ARG		  	; Address of destination #2
	sub	rdi, rsi			; Calculate distance from first number to destination
	vxorpd	ymm2, ymm2, ymm2		; Init 4 carry registers
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
as0nrzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
as1nrzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
as2nrzp: ynorm_addsub_wpn_zpad noexec, noexec, rdi ; Add & subtract and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	rdx, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	as2nrzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdx, pass1blkdst
	bump	rsi, 64
	bump	rdx, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	as2nrzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	as1nrzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	as0nrzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries2_zpad		; Save carries in the carries array

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubnrzp3 ENDP


;;
;; Copy one number and zero some low order words.
;;

PROCFL	gwycopyzero3
	ad_prolog 0,0,rbx,rsi,rdi,rbp
	mov	rsi, SRCARG		; Address of first number
	mov	rdi, DESTARG		; Address of destination
	sub	ecx, ecx		; Offset to compare to COPYZERO

	mov	al, -1			; Create 4 masks for copying values
	mov	BYTE PTR YMM_TMP4[7], al ; Create the copy all four values mask
	mov	BYTE PTR YMM_TMP4[15], al
	mov	BYTE PTR YMM_TMP4[23], al
	mov	BYTE PTR YMM_TMP4[31], al
	mov	BYTE PTR YMM_TMP3[7], cl ; Create the copy three values mask
	mov	BYTE PTR YMM_TMP3[15], al
	mov	BYTE PTR YMM_TMP3[23], al
	mov	BYTE PTR YMM_TMP3[31], al
	mov	BYTE PTR YMM_TMP2[7], cl ; Create the copy two values mask
	mov	BYTE PTR YMM_TMP2[15], cl
	mov	BYTE PTR YMM_TMP2[23], al
	mov	BYTE PTR YMM_TMP2[31], al
	mov	BYTE PTR YMM_TMP1[7], cl ; Create the copy one value mask
	mov	BYTE PTR YMM_TMP1[15], cl
	mov	BYTE PTR YMM_TMP1[23], cl
	mov	BYTE PTR YMM_TMP1[31], al
	vxorpd	ymm1, ymm1, ymm1	; Start with the copy zero values mask

	mov	ebx, addcount1		; Load pass 1 blk count
cz1:	mov	ebp, normval4		; Load count of 4KB chunks in a block
cz2:	mov	edx, normval1		; Load count of clms in 4KB
cz3:	mov	eax, cache_line_multiplier ; Compute cache lines in 4KB chunk
cz4:	ycopyzero
	add	rsi, pass2blkdst	; Next src ptr
	add	rdi, pass2blkdst	; Next dest ptr
	add	rcx, pass2blkdst	; Next comparison offset
	add	al, 256/4		; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	short cz4
	sub	rsi, pass1blkdst	; Adjust source pointers
	sub	rdi, pass1blkdst
	sub	rcx, pass1blkdst
	bump	rsi, 64			; Next source
	bump	rdi, 64			; Next dest
	bump	rcx, 64			; Next compare offset
	sub	al, 4			; Loop clm times
	jnz	cz4			; Loop if necessary
	dec	rdx			; Test for end of 4KB chunk
	jnz	cz3			; Loop if necessary
	add	rsi, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rdi, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	add	rcx, fourKBgapsize	; Skip 64 to 192 bytes every 4KB
	dec	rbp			; Test loop counter
	jnz	cz2			; Loop if necessary
	add	rsi, pass2gapsize
	add	rdi, pass2gapsize
	add	rcx, pass2gapsize
	sub	rsi, pass2blkdst
	sub	rdi, pass2blkdst
	sub	rcx, pass2blkdst
	add	rsi, pass1blkdst
	add	rdi, pass1blkdst
	add	rcx, pass1blkdst
	dec	rbx			; Test loop counter
	jnz	cz1			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi,rbp
gwycopyzero3 ENDP

;;
;; Multiply a number by a small value with carry propagation (eight different versions)
;;

loopcount3	EQU	DPTR [rsp+first_local+0*SZPTR]
loopcount4	EQU	DPTR [rsp+first_local+0*SZPTR+4]
loopcount5	EQU	DPTR [rsp+first_local+0*SZPTR+8]

	; Base 2, irrational version, not zero-padded
PROCFL	gwymuls3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2:	ynorm_smallmul_wpn exec, exec		; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymuls3 ENDP


	; Base 2, rational version, not zero-padded
PROCFL	gwymulsr3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0r:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1r:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2r:	ynorm_smallmul_wpn noexec, exec		; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2r
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2r
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1r				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0r

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulsr3 ENDP


	; Not base 2, irrational version, not zero-padded
PROCFL	gwymulsn3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0n:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1n:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2n:	ynorm_smallmul_wpn exec, noexec		; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2n
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2n
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1n				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0n

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulsn3 ENDP


	; Not base 2, rational version, not zero-padded
PROCFL	gwymulsnr3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0nr:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1nr:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2nr:	ynorm_smallmul_wpn noexec, noexec	; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2nr
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2nr
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1nr				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0nr

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulsnr3 ENDP


	; Base 2, irrational version, zero-padded
PROCFL	gwymulszp3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0zp:	mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1zp:	mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2zp:	ynorm_smallmul_wpn_zpad exec, exec	; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2zp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0zp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulszp3 ENDP


	; Base 2, rational version, zero-padded
PROCFL	gwymulsrzp3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vmovapd	ymm2, YMM_BIGVAL		; Init 2 carry registers
	vmovapd	ymm3, ymm2

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0rzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1rzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2rzp: ynorm_smallmul_wpn_zpad noexec, exec	; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2rzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2rzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1rzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0rzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulsrzp3 ENDP


	; Non base 2, irrational version, zero-padded
PROCFL	gwymulsnzp3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	mov	rbp, PREMULT_PREFETCH		; Addr of the group multipliers
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0nzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1nzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2nzp: ynorm_smallmul_wpn_zpad exec, noexec	; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	bump	rdi, 2				; Next flags ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2nzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2nzp
	add	rdi, normval2			; Adjust ptr to little/big flags
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1nzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0nzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulsnzp3 ENDP


	; Non base 2, rational version, zero-padded
PROCFL	gwymulsnrzp3
	ad_prolog 12,0,rbx,rbp,rsi,rdi,xmm6
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm6, DBLARG		; Load small multiplier value
	vxorpd	ymm2, ymm2, ymm2		; Init 2 carry registers
	vxorpd	ymm3, ymm3, ymm3

	;; Do a pass 1 block

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
	mov	loopcount3, eax
mul0nrzp: mov	eax, normval1			; Load count of clms in 4KB
	mov	loopcount4, eax
mul1nrzp: mov	eax, cache_line_multiplier	; Load inner loop count
	mov	loopcount5, eax			; Save inner loop count
mul2nrzp: ynorm_smallmul_wpn_zpad noexec, noexec ; Multiply and normalize 8 values
	add	rsi, pass2blkdst		; Next src ptr
	add	BYTE PTR loopcount5, 256/4	; Loop 4 times (4 pass 2 blocks in each pass 1 block)
	jnc	mul2nrzp
	sub	rsi, pass1blkdst		; Adjust source pointers
	bump	rsi, 64
	sub	loopcount5, 4			; Loop clm times
	jnz	mul2nrzp
	dec	loopcount4			; Loop until 4KB processed
	jnz	mul1nrzp			; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	loopcount3			; Loop until pass 1 block completed
	jnz	mul0nrzp

	;; Pass 1 block done, save carries

	ynorm_op_wpn_save_carries_zpad		; Save carries in the carries array

	ad_epilog 12,0,rbx,rbp,rsi,rdi,xmm6
gwymulsnrzp3 ENDP

;;
;; Do final carry propagation for add/sub/addsub/smallmul operations
;;

PROCFL	gwy3_apply_carries
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rax, YMM_CARRIES_ROUTINE
	call	rax
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwy3_apply_carries ENDP

_TEXT	ENDS
END
