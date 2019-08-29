; Copyright 2011-2019 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu

;
; Additional routines used with zr4dwpn FFTs
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

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwzaddq3
	ad_prolog 0,0,rbx,rsi,rdi,r8
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rdi, fourKBgapsize	; Preload 4KB gap distance
	sub	rdi, normblkdst4	; One pass FFTs do not apply a 4KB gap distance
	mov	ebx, 8			; Eight pass 2 blks in one pass 1 block
qadd0:	mov	eax, normval4		; Load count of 4KB chunks in a pass 2 block
qadd1:	mov	r8d, normval1		; Load count of loword/highword pairs in 4KB chunk
qaddlp:	vmovapd	zmm0, [rdx]		; Load second number
	vmovapd	zmm1, [rdx+64]		; Load second number
	vaddpd	zmm0, zmm0, [rcx]	; Add in first number
	vaddpd	zmm1, zmm1, [rcx+64]	; Add in first number
;; bug -- unroll this a bit, use common register naming rsi=src, rbx = dest
	zstore	[rsi], zmm0		; Save result
	zstore	[rsi+64], zmm1		; Save result
	bump	rcx, 128		; Next source
	bump	rdx, 128		; Next source
	bump	rsi, 128		; Next dest
	dec	r8			; Test for end of 4KB chunk
	jnz	short qaddlp		; Loop if necessary
	add	rcx, rdi		; Skip 64 to 192 bytes every 4KB
	add	rdx, rdi		; Skip 64 to 192 bytes every 4KB
	add	rsi, rdi		; Skip 64 to 192 bytes every 4KB
	dec	rax			; Check middle loop counter
	jnz	short qadd1		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	qadd0			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi,r8
gwzaddq3 ENDP

;;
;; Add two numbers with carry propagation (four different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzadd3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
add0:	mov	r8, normval2			; Load count of clms in 4KB chunk
add1:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
add2:	znorm_op_wpn vaddpd, exec		; Add and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	bump	rdi, 1				; Next big/lit ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	add2
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust pointers
	sub	rcx, r9
	sub	rbx, r9
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	add2
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	add1				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	dec	rax				; Loop until pass 1 block completed
	jnz	add0

	znorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzadd3 ENDP


	; Rational, not zero-padded
PROCFL	gwzaddr3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_preload noexec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
add0r:	mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
add2r:	znorm_op_wpn vaddpd, noexec		; Add and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	add2r
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust pointers
	sub	rcx, r9
	sub	rbx, r9
	dec	r8				; Loop until 4KB processed
	jnz	add2r				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	add0r

	znorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddr3 ENDP


	; Irrational, zero-padded
PROCFL	gwzaddzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_zpad_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index
	sub	r9, r9				; Clear register used to alternate the increment of rdi

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
add0zp:	mov	r8, normval2			; Load count of clms in 4KB chunk
add1zp:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
add2zp:	znorm_op_wpn_zpad vaddpd, exec		; Add and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rcx, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	add	rdi, r9				; Increment (or not) the biglit table pointer
	xor	r9, 1				; Alternate the increment for rdi
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	add2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rcx, pass1blkdst
	sub	rbx, pass1blkdst
	bump	rsi, 128
	bump	rcx, 128
	bump	rbx, 128
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	add2zp
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	add1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rcx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	add0zp

	znorm_op_wpn_zpad_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddzp3 ENDP


	; Rational, zero-padded
PROCFL	gwzaddrzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_zpad_preload noexec	; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
add0rzp:mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
add2rzp:znorm_op_wpn_zpad vaddpd, noexec	; Add and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	add2rzp
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	dec	r8				; Loop until 4KB processed
	jnz	add2rzp				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	add0rzp

	znorm_op_wpn_zpad_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddrzp3 ENDP


;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwzsubq3
	ad_prolog 0,0,rbx,rsi,rdi,r8
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rdi, fourKBgapsize	; Preload 4KB gap distance
	sub	rdi, normblkdst4	; One pass FFTs do not apply a 4KB gap distance
	mov	ebx, 8			; Eight pass 2 blks in one pass 1 block
qsub0:	mov	eax, normval4		; Load count of 4KB chunks in a pass 2 block
qsub1:	mov	r8d, normval1		; Load count of loword/highword pairs in 4KB chunk
qsublp:	vmovapd	zmm0, [rdx]		; Load second number
	vmovapd	zmm1, [rdx+64]		; Load second number
	vsubpd	zmm0, zmm0, [rcx]	; Subtract first number
	vsubpd	zmm1, zmm1, [rcx+64]	; Subtract first number
	zstore	[rsi], zmm0		; Save result
	zstore	[rsi+64], zmm1		; Save result
	bump	rcx, 128		; Next source
	bump	rdx, 128		; Next source
	bump	rsi, 128		; Next dest
	dec	r8			; Test for end of 4KB chunk
	jnz	short qsublp		; Loop if necessary
	add	rcx, rdi		; Skip 64 to 192 bytes every 4KB
	add	rdx, rdi		; Skip 64 to 192 bytes every 4KB
	add	rsi, rdi		; Skip 64 to 192 bytes every 4KB
	dec	rax			; Check middle loop counter
	jnz	short qsub1		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	qsub0			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi,r8
gwzsubq3 ENDP

;;
;; Subtract two numbers with carry propagation (four different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzsub3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
sub0:	mov	r8, normval2			; Load count of clms in 4KB chunk
sub1:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
sub2:	znorm_op_wpn vsubpd, exec		; Subtract and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	bump	rdi, 1				; Next big/lit ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	sub2
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	sub2
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	sub1				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	sub0

	znorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsub3 ENDP


	; Rational, not zero-padded
PROCFL	gwzsubr3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_preload noexec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
sub0r:	mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
sub2r:	znorm_op_wpn vsubpd, noexec		; Subtract and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	sub2r
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	dec	r8				; Loop until 4KB processed
	jnz	sub2r				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	sub0r

	znorm_op_wpn_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsubr3 ENDP


	; Irrational, zero-padded
PROCFL	gwzsubzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_zpad_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index
	sub	r9, r9				; Clear register used to alternate the increment of rdi

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
sub0zp:	mov	r8, normval2			; Load count of clms in 4KB chunk
sub1zp:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
sub2zp:	znorm_op_wpn_zpad vsubpd, exec		; Subtract and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rcx, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	add	rdi, r9				; Increment (or not) the biglit table pointer
	xor	r9, 1				; Alternate the increment for rdi
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	sub2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rcx, pass1blkdst
	sub	rbx, pass1blkdst
	bump	rsi, 128
	bump	rcx, 128
	bump	rbx, 128
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	sub2zp
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	sub1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rcx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	sub0zp

	znorm_op_wpn_zpad_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsubzp3 ENDP


	; Rational, zero-padded
PROCFL	gwzsubrzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_wpn_zpad_preload noexec	; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
sub0rzp:mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
sub2rzp:znorm_op_wpn_zpad vsubpd, noexec	; Subtract and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	sub2rzp
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	dec	r8				; Loop until 4KB processed
	jnz	sub2rzp				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	sub0rzp

	znorm_op_wpn_zpad_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzsubrzp3 ENDP


;;
;; Add and subtract two numbers without carry propagation
;;

PROCFL	gwzaddsubq3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	rdi, fourKBgapsize	; Preload 4KB gap distance
	sub	rdi, normblkdst4	; One pass FFTs do not apply a 4KB gap distance
	mov	ebx, 8			; Eight pass 2 blks in one pass 1 block
qaddsub0:mov	eax, normval4		; Load count of 4KB chunks in a pass 2 block
qaddsub1:mov	r8d, normval1		; Load count of loword/highword pairs in 4KB chunk
qaddsublp:
	vmovapd	zmm1, [rcx]		; Load first number
	vmovapd	zmm3, [rcx+64]		; Load first number
	vaddpd	zmm0, zmm1, [rdx]	; Add in second number
	vsubpd	zmm1, zmm1, [rdx]	; Subtract out second number
	vaddpd	zmm2, zmm3, [rdx+64]	; Add in second number
	vsubpd	zmm3, zmm3, [rdx+64]	; Subtract out second number
	zstore	[rsi], zmm0		; Save result
	zstore	[rbp], zmm1		; Save result
	zstore	[rsi+64], zmm2		; Save result
	zstore	[rbp+64], zmm3		; Save result
	bump	rcx, 128		; Next source
	bump	rdx, 128		; Next source
	bump	rsi, 128		; Next dest
	bump	rbp, 128		; Next dest
	dec	r8			; Test for end of 4KB chunk
	jnz	short qaddsublp		; Loop if necessary
	add	rcx, rdi		; Skip 64 to 192 bytes every 4KB
	add	rdx, rdi		; Skip 64 to 192 bytes every 4KB
	add	rsi, rdi		; Skip 64 to 192 bytes every 4KB
	add	rbp, rdi		; Skip 64 to 192 bytes every 4KB
	dec	rax			; Check middle loop counter
	jnz	qaddsub1		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	qaddsub0		; Loop if necessary
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8
gwzaddsubq3 ENDP

;;
;; Add and subtract two numbers with carry propagation (four different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzaddsub3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG		  	; Address of destination #2

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_wpn_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
as0:	mov	r8, normval2			; Load count of clms in 4KB chunk
as1:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
as2:	znorm_addsub_wpn exec			; Add & subtract and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	rbp, r9				; Next dest ptr
	bump	rdi, 1				; Next big/lit ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	as2
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	sub	rbp, r9
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	as2
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	as1				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbp, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	as0

	znorm_addsub_wpn_save_carries		; Save carries in the two carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsub3 ENDP


	; Rational, not zero-padded
PROCFL	gwzaddsubr3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG		  	; Address of destination #2

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_wpn_preload noexec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
as0r:	mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
as2r:	znorm_addsub_wpn noexec			; Add & subtract and normalize four low/high FFT word pair (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	rbp, r9				; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	as2r
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	sub	rbp, r9
	dec	r8				; Loop until 4KB processed
	jnz	as2r				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbp, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	as0r

	znorm_addsub_wpn_save_carries		; Save carries in the two carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsubr3 ENDP


	; Irrational, zero-padded
PROCFL	gwzaddsubzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG		  	; Address of destination #2

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_wpn_zpad_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index
	sub	r9, r9				; Clear register used to alternate the increment of rdi

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
as0zp:	mov	r8, normval2			; Load count of clms in 4KB chunk
as1zp:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
as2zp:	znorm_addsub_wpn_zpad exec		; Add & subtract and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rcx, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	add	rbp, pass2blkdst		; Next dest ptr
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	add	rdi, r9				; Increment (or not) the biglit table pointer
	xor	r9, 1				; Alternate the increment for rdi
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	as2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rcx, pass1blkdst
	sub	rbx, pass1blkdst
	sub	rbp, pass1blkdst
	bump	rsi, 128
	bump	rcx, 128
	bump	rbx, 128
	bump	rbp, 128
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	as2zp
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	as1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rcx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbp, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	as0zp

	znorm_addsub_wpn_zpad_save_carries	; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsubzp3 ENDP


	; Rational, zero-padded
PROCFL	gwzaddsubrzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of first number
	mov	rcx, SRC2ARG			; Address of second number
	mov	rbx, DESTARG			; Address of destination
	mov	rbp, DEST2ARG		  	; Address of destination #2

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_wpn_zpad_preload noexec	; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
as0rzp:	mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
as2rzp:	znorm_addsub_wpn_zpad noexec		; Add & subtract and normalize four low/high FFT word pairs (64 values)
	mov	r9, pass2blkdst
	add	rsi, r9				; Next src ptr
	add	rcx, r9				; Next src ptr
	add	rbx, r9				; Next dest ptr
	add	rbp, r9				; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	as2rzp
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rcx, r9
	sub	rbx, r9
	sub	rbp, r9
	dec	r8				; Loop until 4KB processed
	jnz	as2rzp				; Loop til done
	mov	r9, fourKBgapsize
	add	rsi, r9				; Skip 64 to 192 bytes every 4KB
	add	rcx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbx, r9				; Skip 64 to 192 bytes every 4KB
	add	rbp, r9				; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	as0rzp

	znorm_addsub_wpn_zpad_save_carries	; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzaddsubrzp3 ENDP

;;
;; Add in a small number with carry propagation (four different versions)
;;

	; Irrational version
PROCFL	gwzadds3
	ad_prolog 0,0,rbx,rsi
	mov	rsi, DESTARG		; Address of destination
	vmovsd	xmm5, DBLARG		; Small addin value
	znorm_common_op_wpn_preload exec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smalladd_wpn_preload exec
	znorm_smalladd_wpn exec
	ad_epilog 0,0,rbx,rsi
gwzadds3 ENDP

	; Rational version
PROCFL	gwzaddsr3
	ad_prolog 0,0,rbx,rsi
	mov	rsi, DESTARG		; Address of destination
	vmovsd	xmm5, DBLARG		; Small addin value
	znorm_common_op_wpn_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smalladd_wpn_preload noexec
	znorm_smalladd_wpn noexec
	ad_epilog 0,0,rbx,rsi
gwzaddsr3 ENDP

;;
;; Multiply a number by a small value with carry propagation (four different versions)
;;

	; Irrational version, not zero-padded
PROCFL	gwzmuls3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of destination
	mov	rbx, DESTARG			; Address of destination

	vbroadcastsd zmm31, DBLARG		; Load small multiplier value

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_wpn_preload exec		; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
mul0:	mov	r8, normval2			; Load count of clms in 4KB chunk
mul1:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
mul2:	znorm_smallmul_wpn exec			; Multiply and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	bump	rdi, 1				; Next big/lit ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	mul2
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rbx, r9
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	mul2
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	mul1				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	mul0

	znorm_smallmul_wpn_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmuls3 ENDP


	; Rational version, not zero-padded
PROCFL	gwzmulsr3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of destination
	mov	rbx, DESTARG			; Address of destination

	vbroadcastsd zmm31, DBLARG		; Load small multiplier value

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_wpn_preload noexec	; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
mul0r:	mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
mul2r:	znorm_smallmul_wpn noexec		; Multiply and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	mul2r
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rbx, r9
	dec	r8				; Loop until 4KB processed
	jnz	mul2r				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	mul0r

	znorm_smallmul_wpn_save_carries		; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmulsr3 ENDP


	; Irrational version, zero-padded
PROCFL	gwzmulszp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of destination
	mov	rbx, DESTARG			; Address of destination

	vbroadcastsd zmm31, DBLARG		; Load small multiplier value

	znorm_common_op_wpn_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_wpn_zpad_preload exec	; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4
	mov	rdi, norm_ptr1			; Addr of the big/little flags array
	mov	r12, compressed_biglits		; Get pointer to compressed biglit table
	sub	rdx, rdx			; Clear register used to load compressed biglit index
	sub	r9, r9				; Clear register used to alternate the increment of rdi

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
mul0zp:	mov	r8, normval2			; Load count of clms in 4KB chunk
mul1zp:	mov	r10d, cache_line_multiplier	; Load clm*8.  This affects how we step through the big/lit array.
mul2zp:	znorm_smallmul_wpn_zpad exec		; Multiply and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	xor	r12, 4				; Bump or unbump pointer into the compressed biglit table
	add	rdi, r9				; Increment (or not) the biglit table pointer
	xor	r9, 1				; Alternate the increment for rdi
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	mul2zp
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rbx, pass1blkdst
	bump	rsi, 128
	bump	rbx, 128
	sub	r10, 8				; We just finished 8 of the cache_line_multiplier
	jnz	mul2zp
	add	rdi, normval3			; Next big/lit ptr
	dec	r8				; Loop until 4KB processed
	jnz	mul1zp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	mul0zp

	znorm_smallmul_wpn_zpad_save_carries	; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmulszp3 ENDP


	; Rational version, zero-padded
PROCFL	gwzmulsrzp3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rsi, SRCARG			; Address of destination
	mov	rbx, DESTARG			; Address of destination

	vbroadcastsd zmm31, DBLARG		; Load small multiplier value

	znorm_common_op_wpn_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_wpn_zpad_preload noexec	; Preload constants, init carry registers

	mov	r13, pass1blkdst		; Distance to FFT source #2
	lea	r14, [r13+r13]			; Distance to FFT source #3
	lea	r15, [r13+r14]			; Distance to FFT source #4

	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
mul0rzp:mov	r8d, normval1			; Load count of loword/highword pairs in 4KB chunk
mul2rzp:znorm_smallmul_wpn_zpad noexec		; Multiply and normalize four low/high FFT word pairs (64 values)
	add	rsi, pass2blkdst		; Next src ptr
	add	rbx, pass2blkdst		; Next dest ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	mul2rzp
	lea	r9, [r13-128]			; pass1blkdst - 128
	sub	rsi, r9				; Adjust source pointers
	sub	rbx, r9
	dec	r8				; Loop until 4KB processed
	jnz	mul2rzp				; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rbx, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	mul0rzp

	znorm_smallmul_wpn_zpad_save_carries	; Save carries in the carries array

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwzmulsrzp3 ENDP

;;
;; Do final carry propagation for add/sub/addsub/smallmul operations
;;

	; Irrational version
PROCFL	gwz3_apply_carries
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
	mov	rax, ZMM_OP_CARRIES_ROUTINE
	call	rax
	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r13,r14,r15,zmm6,zmm7,zmm8,zmm9,zmm10,zmm11,zmm12,zmm13,zmm14,zmm15
gwz3_apply_carries ENDP


;;
;; Copy a number zeroing some low order words
;;

PROCFL	gwzcopyzero3
	ad_prolog 0,0,rbx,rsi,rdi
	mov	rsi, SRCARG			; Address of first number
	mov	rdi, DESTARG			; Address of destination
	mov	rdx, SRC2ARG			; Load address of masks and counts computed in C code
	sub	ecx, ecx			; Count of iterations before new masks are needed
	mov	eax, normval4			; Load count of 4KB chunks in a pass 2 block
cpzlp0:	mov	ebx, normval1			; Load count of loword/highword pairs in 4KB chunk
cpzlp2:	and	ecx, ecx			; Time to load new masks?
	jnz	short noload			; No, wait for count to reach zero
	kmovw	k1, [rdx+0]			; Load 2 new masks
	kshiftrw k2, k1, 8
	mov	ecx, [rdx+4]			; Count of iterations these masks are valid
	bump	rdx, 8				; Move onto next set of masks
noload:	vmovapd	zmm0 {k1}{z}, [rsi]		; Load source data
	vmovapd	zmm2 {k2}{z}, [rsi+64]
	zstore	[rdi], zmm0			; Store destination data
	zstore	[rdi+64], zmm2
	bump	rcx, -1				; Decrement mask count
	add	rsi, pass2blkdst		; Next src ptr
	add	rdi, pass2blkdst		; Next src ptr
	add	eax, 0x80000000/4		; Loop 8 times (8 pass 2 blocks in each pass 1 block)
	jnc	short cpzlp2
	sub	rsi, pass1blkdst		; Adjust source pointers
	sub	rdi, pass1blkdst
	bump	rsi, 128
	bump	rdi, 128
	dec	rbx				; Loop until 4KB processed
	jnz	short cpzlp2			; Loop til done
	add	rsi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	add	rdi, fourKBgapsize		; Skip 64 to 192 bytes every 4KB
	dec	eax				; Loop until pass 1 block completed
	jnz	cpzlp0
	ad_epilog 0,0,rbx,rsi,rdi
gwzcopyzero3 ENDP


_TEXT	ENDS
END
