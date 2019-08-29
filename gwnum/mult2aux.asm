; Copyright 1995-2011 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

	TITLE   setup

	.686
	.XMM
	.MODEL	FLAT

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE mult.mac
INCLUDE memory.mac
INCLUDE normal.mac

_TEXT SEGMENT

	flat_distances

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwaddq2
	ad_prolog 0,0,rbx,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
uadd0:	mov	eax, normval4		; Load count of 4KB pages in a block
uaddlp:	fld	QWORD PTR [ecx]		; Load first number
	fadd	QWORD PTR [edx]		; Add in second number
	fld	QWORD PTR [ecx+8]	; Load first number
	fadd	QWORD PTR [edx+8]	; Add in second number
	fld	QWORD PTR [ecx+16]	; Load first number
	fadd	QWORD PTR [edx+16]	; Add in second number
	fxch	st(1)
	fld	QWORD PTR [ecx+24]	; Load first number
	fadd	QWORD PTR [edx+24]	; Add in second number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+32]		; Bump source pointer
	lea	edx, [edx+32]		; Bump source pointer
	lea	esi, [esi+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short uaddlp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short uaddlp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short uadd0		; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwaddq2	ENDP


;;
;; Add two numbers with carry propagation
;;

loopcount1	EQU	DPTR [rsp+first_local+12]
loopcount2	EQU	DPTR [rsp+first_local+8]
loopcount3	EQU	DPTR [rsp+first_local+4]
loopcount4	EQU	DPTR [rsp+first_local]

PROCFL	gwadd2
	ad_prolog 16,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	eax, addcount1		; Load block count
	mov	loopcount1, eax
ablk:	mov	eax, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, eax
	mov	ebx, norm_col_mults	; Addr of the column multipliers
iadd0:	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	sub	eax, eax		; Clear big/lit flag
iadd2:	norm_op_2d fadd			; Add and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	iadd2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	iadd1			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	iadd0			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	ablk 			; Loop til done

	;; All blocks done

	mov	esi, DESTARG		; Addr of FFT data
	mov	ebp, norm_grp_mults	; Addr of the group multipliers
	norm_op_2d_cleanup		; Add 2 carries to start of fft

	ad_epilog 16,0,rbx,rbp,rsi,rdi
gwadd2	ENDP


;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwsubq2
	ad_prolog 0,0,rbx,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
usub0:	mov	eax, normval4		; Load count of 4KB pages in a block
usublp:	fld	QWORD PTR [edx]		; Load second number
	fsub	QWORD PTR [ecx]		; Subtract first number
	fld	QWORD PTR [edx+8]	; Load second number
	fsub	QWORD PTR [ecx+8]	; Subtract first number
	fld	QWORD PTR [edx+16]	; Load second number
	fsub	QWORD PTR [ecx+16]	; Subtract first number
	fxch	st(1)
	fld	QWORD PTR [edx+24]	; Load second number
	fsub	QWORD PTR [ecx+24]	; Subtract first number
	fxch	st(3)
	fstp	QWORD PTR [esi]		; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+32]		; Bump source pointer
	lea	edx, [edx+32]		; Bump source pointer
	lea	esi, [esi+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short usublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short usublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short usub0		; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwsubq2	ENDP

;;
;; Subtract two numbers with carry propagation
;;

loopcount1	EQU	DPTR [rsp+first_local+12]
loopcount2	EQU	DPTR [rsp+first_local+8]
loopcount3	EQU	DPTR [rsp+first_local+4]
loopcount4	EQU	DPTR [rsp+first_local]

PROCFL	gwsub2
	ad_prolog 16,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	eax, addcount1		; Load block count
	mov	loopcount1, eax

sblk:	mov	eax, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, eax
	mov	ebx, norm_col_mults	; Addr of the column multipliers
isub0:	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	sub	eax, eax		; Clear big/lit flag
isub2:	norm_op_2d fsub			; Add and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	isub2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	isub1			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	isub0			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	sblk 			; Loop til done

	;; All blocks done

	mov	esi, DESTARG		; Addr of FFT data
	mov	ebp, norm_grp_mults	; Addr of the group multipliers
	norm_op_2d_cleanup		; Add 2 carries to start of fft

	ad_epilog 16,0,rbx,rbp,rsi,rdi
gwsub2	ENDP

;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwaddsubq2
	ad_prolog 0,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, DEST2ARG	  	; Address of destination #2
	mov	ebx, addcount1		; Load blk count
uaddsub0:
	mov	eax, normval4		; Load count of 4KB pages in a block
uaddsublp:
	fld	QWORD PTR [ecx]		; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx]		; Add in second number
	fxch	st(1)			; S0,A0
	fsub	QWORD PTR [edx]		; Subtract out second number
	fld	QWORD PTR [ecx+8]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx+8]	; Add in second number
	fxch	st(1)			; S1,A1,S0,A0
	fsub	QWORD PTR [edx+8]	; Subtract out second number
	fld	QWORD PTR [ecx+16]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx+16]	; Add in second number
	fxch	st(1)			; S2,A2,S1,A1,S0,A0
	fsub	QWORD PTR [edx+16]	; Subtract out second number
	fld	QWORD PTR [ecx+24]	; Load first number
	fld	st(0)			; Dup first number
	fadd	QWORD PTR [edx+24]	; Add in second number
	fxch	st(7)			; A0,S3,S2,A2,S1,A1,S0,A3
	fstp	QWORD PTR [esi]		; Save result
	fsub	QWORD PTR [edx+24]	; Subtract out second number
	fxch	st(5)			; S0,S2,A2,S1,A1,S3,A3
	fstp	QWORD PTR [ebp]		; Save result
	fstp	QWORD PTR [ebp+16]	; Save result
	fstp	QWORD PTR [esi+16]	; Save result
	fstp	QWORD PTR [ebp+8]	; Save result
	fstp	QWORD PTR [esi+8]	; Save result
	fstp	QWORD PTR [ebp+24]	; Save result
	fstp	QWORD PTR [esi+24]	; Save result
	lea	ecx, [ecx+32]		; Bump source pointer
	lea	edx, [edx+32]		; Bump source pointer
	lea	esi, [esi+32]		; Bump dest pointer
	lea	ebp, [ebp+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short uaddsublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	lea	ebp, [ebp+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short uaddsublp		; Loop if necessary
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	uaddsub0		; Loop if necessary
	ad_epilog 0,0,rbx,rbp,rsi,rdi
gwaddsubq2 ENDP

;;
;; Add and subtract two numbers with carry propagation
;;

loopcount1	EQU	DPTR [rsp+first_local]
loopcount2	EQU	DPTR [rsp+first_local+4]
loopcount3	EQU	DPTR [rsp+first_local+8]
loopcount4	EQU	DPTR [rsp+first_local+12]
save_reg	EQU	DPTR [rsp+first_local+16]

PROCFL	gwaddsub2
	ad_prolog 20,0,rbx,rbp,rsi,rdi
	mov	ecx, SRCARG		; Address of first number
	mov	edx, SRC2ARG		; Address of second number
	mov	esi, DESTARG		; Address of destination
	mov	ebp, DEST2ARG	  	; Address of destination #2
	mov	eax, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	fld	BIGVAL
	fld	BIGVAL
	mov	ebx, addcount1		; Load block count
	mov	loopcount1, ebx
asblk:	mov	ebx, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, ebx
	mov	ebx, norm_col_mults	; Addr of the column multipliers
ias0:	mov	save_reg, eax
	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
	mov	eax, save_reg
ias1:	mov	save_reg, eax
	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	mov	eax, save_reg
ias2:	norm_addsub_2d			; Add and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	ias2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	ias1			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every 4KB
	lea	edx, [edx+64]		; Skip 64 bytes every 4KB
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	lea	ebp, [ebp+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	ias0			; Loop til done
	lea	ecx, [ecx+64]		; Skip 64 bytes every blk
	lea	edx, [edx+64]		; Skip 64 bytes every blk
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+64]		; Skip 64 bytes every blk
	lea	eax, [eax+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	asblk 			; Loop til done

	;; All blocks done

	mov	esi, DESTARG		; Addr of add FFT data
	mov	ebp, DEST2ARG	  	; Addr of sub FFT data
	mov	ebx, norm_grp_mults	; Addr of the group multipliers
	norm_addsub_2d_cleanup		; Add 2 carries to start of fft

	ad_epilog 20,0,rbx,rbp,rsi,rdi
gwaddsub2 ENDP

;;
;; Copy one number and zero some low order words.
;;

PROCFL	gwcopyzero2
	ad_prolog 0,0,rbx,rsi,rdi
	mov	esi, SRCARG		; Address of first number
	mov	edi, DESTARG		; Address of destination
	mov	ecx, NUMARG		; Count of words to zero
	mov	ebx, addcount1		; Load blk count
zlp0:	mov	eax, normval4		; Load count of 4KB pages in a block
zlp:	dec	ecx			; Decrement count of words to zero
	js	short c2		; If negative we are now copying
	fldz
	dec	ecx			; Decrement count of words to zero
	js	short c1		; If negative we are now copying
	fldz
	jmp	short c0
c2:	fld	QWORD PTR [esi]		; Load number
c1:	fld	QWORD PTR [esi+16]	; Load number
c0:	fstp	QWORD PTR [edi+16]	; Save result
	fstp	QWORD PTR [edi]		; Save result
	fld	QWORD PTR [esi+8]	; Copy high word
	fstp	QWORD PTR [edi+8]
	fld	QWORD PTR [esi+24]	; Copy high word
	fstp	QWORD PTR [edi+24]
	lea	esi, [esi+32]		; Bump source pointer
	lea	edi, [edi+32]		; Bump dest pointer
	add	eax, 80000000h/64	; 128 cache lines in a 4KB page
	jnc	short zlp		; Loop if necessary
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	lea	edi, [edi+64]		; Skip 64 bytes every 4KB
	dec	eax			; Check middle loop counter
	jnz	short zlp		; Loop if necessary
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	edi, [edi+64]		; Skip 64 bytes every blk
	dec	ebx			; Check outer/restore inner count
	jnz	short zlp0		; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwcopyzero2 ENDP


;;
;; Mul by a small value with carry propagation
;;

loopcount1	EQU	DPTR [rsp+first_local+12]
loopcount2	EQU	DPTR [rsp+first_local+8]
loopcount3	EQU	DPTR [rsp+first_local+4]
loopcount4	EQU	DPTR [rsp+first_local]

PROCFL	gwmuls2
	ad_prolog 16,0,rbx,rbp,rsi,rdi
	mov	esi, DESTARG		; Address of destination
	fld	DBLARG			; Load small value
	fmul	NORM012_FF		; Mul by two-to-minus-phi fudge
	fstp	TMP5			; Save multiplier
	mov	ebp, norm_grp_mults	; Address of group multipliers
	mov	edi, norm_biglit_array	; Addr of the big/little flags array
	fld	BIGVAL			; Start process with no carry
	fld	BIGVAL
	mov	eax, addcount1		; Load block count
	mov	loopcount1, eax
mblk:	mov	eax, normval4		; Load outer count (4KB pages in a blk)
	mov	loopcount2, eax
	mov	ebx, norm_col_mults	; Addr of the column multipliers
imul0:	mov	eax, normval1		; Load middle count (clms in 4KB page)
	mov	loopcount3, eax		; Save middle count
imul1:	mov	eax, cache_line_multiplier ; Load inner loop count (clm)
	mov	loopcount4, eax		; Save inner loop count
	sub	eax, eax		; Clear big/lit flag
imul2:	norm_smallmul_2d		; Mul and normalize 4 values
	sub	loopcount4, 1		; Decrement inner loop counter
	jnz	imul2 			; Loop til done
	add	edi, normval2		; Adjust ptr to little/big flags
	sub	loopcount3, 1		; Decrement middle loop counter
	jnz	imul1			; Loop til done
	lea	esi, [esi+64]		; Skip 64 bytes every 4KB
	sub	loopcount2, 1		; Decrement outer loop counter
	jnz	imul0			; Loop til done
	lea	esi, [esi+64]		; Skip 64 bytes every blk
	lea	ebp, [ebp+2*16]		; Next set of 2 group multipliers
	add	edi, normval3		; Adjust little/big flags ptr
	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	mblk 			; Loop til done

	;; All blocks done

	mov	esi, DESTARG		; Addr of FFT data
	mov	ebp, norm_grp_mults	; Addr of the group multipliers
	norm_op_2d_cleanup		; Add 2 carries to start of fft

	ad_epilog 16,0,rbx,rbp,rsi,rdi
gwmuls2	ENDP


_TEXT	ENDS
END
