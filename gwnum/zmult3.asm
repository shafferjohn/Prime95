; Copyright 2011-2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routines implement some common cleanup code for AVX-512 FFTs
;

	TITLE   setup

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE znormal.mac
INCLUDE znormal_zpad.mac

_TEXT SEGMENT

;;***********************************************************
;; Routines for processing carries out of a pass 1 data block
;;***********************************************************

;; Irrational, no zero padding, carries from a squaring/multiply/add/sub/addsub/smallmul operation
PROCFL	zgw_carries_wpn3
	int_prolog 0,0,0

	cmp	THIS_BLOCK, 0		; Are we carrying into the first data block?
	jne	i3_not_block0_a		; If not, skip wrapping carries and negating the last carry
	zrotate_carries_array		; Rotate the entire carries array
	zprocess_last_two_carries exec	; Handle the last two carries
i3_not_block0_a:

	mov	rbp, carries		; Addr of the carries
	mov	rsi, DATA_ADDR		; Addr of the FFT data
	mov	r14, pass1blkdst	; Distance to FFT source #2
	lea	r13, [rsi+2*r14]	; FFT source #3
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	ecx, addcount1		; Load count of pass 1 blocks
	mov	r12, compressed_biglits	; Load address of the compressed biglit table
	sub	rdx, rdx		; Clear register used to load compressed biglit index

;; bug - precalculate this?
mov	r15d, cache_line_multiplier; Cache lines in each pass1 loop
;;shl	r15, 3			; Compute biglit increment

	zadd_carry_rows_preload exec
ilp1:	zadd_carry_rows exec		; Propagate 4 low/high carry pairs
	bump	rbp, 4*128		; Next carries pointer
	lea	rsi, [rsi+4*r14]	; Next FFT source #1
	lea	r13, [r13+4*r14]	; Next FFT source #3
	add	rsi, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	r13, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	rdi, r15		; Next big/little flags pointer
	sub	ecx, 4			; Test loop counter
	jnz	ilp1
	int_epilog 0,0,0
zgw_carries_wpn3 ENDP


;; Rational, no zero padding, carries from a squaring/multiply/add/sub/addsub/smallmul operation
PROCFL	zgw_carries_wpnr3
	int_prolog 0,0,0

	cmp	THIS_BLOCK, 0		; Are we carrying into the first data block?
	jne	r3_not_block0_a		; If not, skip wrapping carries and negating the last carry
	zrotate_carries_array		; Rotate the entire carries array
	zprocess_last_two_carries noexec ; Handle the last two carries
r3_not_block0_a:

	mov	rbp, carries		; Addr of the carries
	mov	rsi, DATA_ADDR		; Addr of the FFT data
	mov	r14, pass1blkdst	; Distance to FFT source #2
	lea	r13, [rsi+2*r14]	; FFT source #3
	mov	ecx, addcount1		; Load count of pass 1 blocks

	zadd_carry_rows_preload noexec
rlp1:	zadd_carry_rows noexec		; Propagate 4 low/high carry pairs
	bump	rbp, 4*128		; Next carries pointer
	lea	rsi, [rsi+4*r14]	; Next FFT source #1
	lea	r13, [r13+4*r14]	; Next FFT source #3
	add	rsi, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	r13, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	sub	ecx, 4			; Test loop counter
	jnz	rlp1			; Next carry row
	int_epilog 0,0,0
zgw_carries_wpnr3 ENDP


;; Irrational, zero padding, carries from a squaring/multiply operation
PROCFL	zgw_carries_wpnzp3
	int_prolog 0,0,0

	cmp	THIS_BLOCK, 0		; Are we carrying into the first data block?
	jne	izp3_not_block0_a	; If not, skip wrapping carries and negating the last carry
	zrotate_carries_array		; Rotate the entire carries array
	zprocess_last_two_carries_zpad exec ; Handle the last two carries
izp3_not_block0_a:

	mov	rbp, carries		; Addr of the carries
	mov	rsi, DATA_ADDR		; Addr of the FFT data
	mov	r14, pass1blkdst	; Distance to FFT source #2
	lea	r13, [rsi+2*r14]	; FFT source #3
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	ecx, addcount1		; Load count of pass 1 blocks
	mov	r12, compressed_biglits	; Load address of the compressed biglit table
	sub	rdx, rdx		; Clear register used to load compressed biglit index

;; bug - precalculate this?
mov	r15d, cache_line_multiplier; Cache lines in each pass1 loop
shr	r15, 1			; Compute biglit increment

	zadd_carry_rows_zpad_preload exec
izplp1: zadd_carry_rows_zpad exec	; Propagate 4 low/high carry pairs
	bump	rbp, 4*128		; Next carries pointer
	lea	rsi, [rsi+4*r14]	; Next FFT source #1
	lea	r13, [r13+4*r14]	; Next FFT source #3
	add	rsi, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	r13, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	rdi, r15		; Next big/little flags pointer
	sub	ecx, 4			; Test loop counter
	jnz	izplp1

;; Do the final step in zero-pad normalization:  Make sure the top 4 words in the high half of the FFT are less than k.

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	izp3_not_block0_b	;; If not, skip wrapping carries and negating the last carry
	zpad_prep_final_div_by_k_after_multiply exec
	zpad_final_div_by_k exec
izp3_not_block0_b:

	int_epilog 0,0,0
zgw_carries_wpnzp3 ENDP


;; Rational, zero padding, carries from a squaring/multiply operation
PROCFL	zgw_carries_wpnrzp3
	int_prolog 0,0,0

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	rzp3_not_block0_a	;; If not, skip wrapping carries and negating the last carry
	zrotate_carries_array		;; Rotate the entire carries array
	zprocess_last_two_carries_zpad noexec ; Handle the last two carries
rzp3_not_block0_a:

	mov	rbp, carries		; Addr of the carries
	mov	rsi, DATA_ADDR		; Addr of the FFT data
	mov	r14, pass1blkdst	; Distance to FFT source #2
	lea	r13, [rsi+2*r14]	; FFT source #3
	mov	ecx, addcount1		; Load count of pass 1 blocks

	zadd_carry_rows_zpad_preload noexec
rzplp1:	zadd_carry_rows_zpad noexec	; Propagate 4 low/high carry pairs
	bump	rbp, 4*128		; Next carries pointer
	lea	rsi, [rsi+4*r14]	; Next FFT source #1
	lea	r13, [r13+4*r14]	; Next FFT source #3
	add	rsi, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	r13, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	sub	ecx, 4			; Test loop counter
	jnz	rzplp1			; Next carry row

;; Do the final step in zero-pad normalization:  Make sure the top 4 words in the high half of the FFT are less than k.

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	rzp3_not_block0_b	;; If not, skip wrapping carries and negating the last carry
	zpad_prep_final_div_by_k_after_multiply noexec
	zpad_final_div_by_k noexec
rzp3_not_block0_b:

	int_epilog 0,0,0
zgw_carries_wpnrzp3 ENDP


;; Irrational, zero padding, carries from an add/sub/addsub/smallmul operation
PROCFL	zgw_carries_op_wpnzp3
	int_prolog 0,0,0

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	iopzp3_not_block0_a	;; If not, skip wrapping carries and negating the last carry
	zrotate_carries_array		;; Rotate the entire carries array
	zprocess_last_two_carries_op_zpad exec ;; Handle the last two carries
iopzp3_not_block0_a:

	mov	rbp, carries		; Addr of the carries
	mov	rsi, DATA_ADDR		; Addr of the FFT data
	mov	r14, pass1blkdst	; Distance to FFT source #2
	lea	r13, [rsi+2*r14]	; FFT source #3
	mov	rdi, norm_ptr1		; Addr of the big/little flags array
	mov	ecx, addcount1		; Load count of pass 1 blocks
	mov	r12, compressed_biglits	; Load address of the compressed biglit table
	sub	rdx, rdx		; Clear register used to load compressed biglit index

;; bug - precalculate this?
mov	r15d, cache_line_multiplier; Cache lines in each pass1 loop
shr	r15, 1			; Compute biglit increment

	zadd_carry_rows_op_zpad_preload exec
iopzplp1: zadd_carry_rows_op_zpad exec	; Propagate 4 low/high carry pairs
	bump	rbp, 4*128		; Next carries pointer
	lea	rsi, [rsi+4*r14]	; Next FFT source #1
	lea	r13, [r13+4*r14]	; Next FFT source #3
	add	rsi, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	r13, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	rdi, r15		; Next big/little flags pointer
	sub	ecx, 4			; Test loop counter
	jnz	iopzplp1

;; Do the final step in zero-pad normalization:  Make sure the top 4 words in the high half of the FFT are less than k.

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	iopzp3_not_block0_b	;; If not, skip wrapping carries and negating the last carry
	zpad_prep_final_div_by_k_after_op exec
	zpad_final_div_by_k exec
iopzp3_not_block0_b:

	int_epilog 0,0,0
zgw_carries_op_wpnzp3 ENDP


;; Rational, zero padding, carries from an add/sub/addsub/smallmul operation
PROCFL	zgw_carries_op_wpnrzp3
	int_prolog 0,0,0

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	ropzp3_not_block0_a	;; If not, skip wrapping carries and negating the last carry
	zrotate_carries_array		;; Rotate the entire carries array
	zprocess_last_two_carries_op_zpad noexec ;; Handle the last two carries
ropzp3_not_block0_a:

	mov	rbp, carries		; Addr of the carries
	mov	rsi, DATA_ADDR		; Addr of the FFT data
	mov	r14, pass1blkdst	; Distance to FFT source #2
	lea	r13, [rsi+2*r14]	; FFT source #3
	mov	ecx, addcount1		; Load count of pass 1 blocks

	zadd_carry_rows_op_zpad_preload noexec
ropzplp1: zadd_carry_rows_op_zpad noexec ; Propagate 4 low/high carry pairs
	bump	rbp, 4*128		; Next carries pointer
	lea	rsi, [rsi+4*r14]	; Next FFT source #1
	lea	r13, [r13+4*r14]	; Next FFT source #3
	add	rsi, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	add	r13, normblkdst4	; Add 0 or 64 every 4KB for one-pass FFTs
	sub	ecx, 4			; Test loop counter
	jnz	ropzplp1		; Next carry row

;; Do the final step in zero-pad normalization:  Make sure the top 4 words in the high half of the FFT are less than k.

	cmp	THIS_BLOCK, 0		;; Are we carrying into the first data block?
	jne	ropzp3_not_block0_b	;; If not, skip wrapping carries and negating the last carry
	zpad_prep_final_div_by_k_after_op noexec
	zpad_final_div_by_k noexec
ropzp3_not_block0_b:

	int_epilog 0,0,0
zgw_carries_op_wpnrzp3 ENDP

_TEXT	ENDS
END
