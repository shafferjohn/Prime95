; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
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
INCLUDE memory.mac
INCLUDE xnormal.mac

; Internal routine to add in the two wraparound carries
; I'd like to make this a subroutine, but it is too difficult
; to get push_amt correct.

final_carries_2 MACRO
	LOCAL	b2c, zpc, b2zpc, c2dn

	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rbx, norm_col_mults	; Addr of the column multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array

	cmp	ZERO_PADDED_FFT, 0	;; Zero-padded FFT?
	jne	zpc			;; Yes, do special zpad carry

	xnorm_top_carry_cmn rsi, xmm7, 2
	sub	rax, rax		; Trashed by xnorm_top_carry_cmn
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2c			; yes, do simpler rounding
	xnorm_smallmul_2d_fft noexec	; Add 2 carries to start of fft
	jmp	c2dn
b2c:	xnorm_smallmul_2d_fft exec	; Add 2 carries to start of fft
	jmp	c2dn

zpc:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2zpc			; yes, do simpler rounding
	xnorm_smallmul_2d_fft_zpad noexec ; Do the special zpad carry not base 2
	jmp	c2dn
b2zpc:	xnorm_smallmul_2d_fft_zpad exec	; Do the special zpad carry base 2
c2dn:
	ENDM

_TEXT SEGMENT

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwxaddq2
	ad_prolog 0,0,rbx,rsi,rdi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
uadd0:	mov	eax, normval4		; Load count of 8KB chunks in a block
uaddlp:	xload	xmm0, [rdx]		; Load second number
	addpd	xmm0, [rcx]		; Add in first number
	xload	xmm1, [rdx+16]		; Load second number
	addpd	xmm1, [rcx+16]		; Add in first number
	xload	xmm2, [rdx+32]		; Load second number
	addpd	xmm2, [rcx+32]		; Add in first number
	xload	xmm3, [rdx+48]		; Load second number
	addpd	xmm3, [rcx+48]		; Add in first number
	xstore	[rsi], xmm0		; Save result
	xstore	[rsi+16], xmm1		; Save result
	xstore	[rsi+32], xmm2		; Save result
	xstore	[rsi+48], xmm3		; Save result
	bump	rcx, 64			; Next source
	bump	rdx, 64			; Next source
	bump	rsi, 64			; Next dest
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	short uaddlp		; Loop if necessary
	bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	rax			; Check middle loop counter
	jnz	short uaddlp		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	uadd0			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwxaddq2 ENDP

;;
;; Add two numbers with carry propagation
;;

saved_reg	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+SZPTR+16]

PROCFL	gwxadd2
	ad_prolog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xload	xmm7, XMM_BIGVAL	; Init 4 carries
	xstore	XMM_TMP1, xmm7
	xstore	XMM_TMP2, xmm7
	xstore	XMM_TMP3, xmm7
	xstore	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

asec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

ablk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rbx, norm_col_mults	; Addr of the column multipliers

add0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2add			; yes, do simpler rounding
nb2add0:cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2iadd0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
nb2radd1:
	xnorm_op_2d addpd, noexec, noexec, saved_reg ; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2radd1		; Loop til done
	jmp	achunkdn		; Jump to chunk done code
nb2iadd0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2iadd1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax	; Save inner loop count
nb2iadd2:xnorm_op_2d addpd, exec, noexec, saved_reg ; Add and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2iadd2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2iadd1		; Loop til done
	jmp	achunkdn		; Jump to chunk done code
b2add:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	iadd0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
radd1:	xnorm_op_2d addpd, noexec, exec, saved_reg ; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	radd1 			; Loop til done
	jmp	achunkdn		; Jump to chunk done code
iadd0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax	; Save inner loop count
iadd2:	xnorm_op_2d addpd, exec, exec, saved_reg ; Add and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	iadd2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	iadd1			; Loop til done

	;; Chunk done

achunkdn:bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	add0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest

	;; Block done

ablkdn:	sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	add	rsi, pass1blkdst
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	askip2			; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 128		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
askip2:	dec	loopcount2		; Decrement outer loop counter
	jnz	ablk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	asec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_2			; Add the carries back in

	ad_epilog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadd2 ENDP

;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwxsubq2
	ad_prolog 0,0,rbx,rsi,rdi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	ebx, addcount1		; Load blk count
usub0:	mov	eax, normval4		; Load count of 8KB chunks in a block
usublp:	xload	xmm0, [rdx]		; Load second number
	subpd	xmm0, [rcx]		; Subtract first number
	xload	xmm1, [rdx+16]		; Load second number
	subpd	xmm1, [rcx+16]		; Subtract first number
	xload	xmm2, [rdx+32]		; Load second number
	subpd	xmm2, [rcx+32]		; Subtract first number
	xload	xmm3, [rdx+48]		; Load second number
	subpd	xmm3, [rcx+48]		; Subtract first number
	xstore	[rsi], xmm0		; Save result
	xstore	[rsi+16], xmm1		; Save result
	xstore	[rsi+32], xmm2		; Save result
	xstore	[rsi+48], xmm3		; Save result
	bump	rcx, 64			; Next source
	bump	rdx, 64			; Next source
	bump	rsi, 64			; Next dest
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	short usublp		; Loop if necessary
	bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	rax			; Check middle loop counter
	jnz	short usublp		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	usub0			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwxsubq2 ENDP

;;
;; Subtract two numbers with carry propagation
;;

saved_reg	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+SZPTR+16]

PROCFL	gwxsub2
	ad_prolog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xload	xmm7, XMM_BIGVAL	; Init 4 carries
	xstore	XMM_TMP1, xmm7
	xstore	XMM_TMP2, xmm7
	xstore	XMM_TMP3, xmm7
	xstore	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

ssec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address

	;; Do a block

sblk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rbx, norm_col_mults	; Addr of the column multipliers

sub0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2sub			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2isub0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
nb2rsub1:
	xnorm_op_2d subpd, noexec, noexec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2rsub1		; Loop til done
	jmp	schunkdn		; Jump to chunk done code
nb2isub0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2isub1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
nb2isub2:xnorm_op_2d subpd, exec, noexec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2isub2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2isub1		; Loop til done
	jmp	schunkdn		; Jump to chunk done code
b2sub:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	isub0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
rsub1:	xnorm_op_2d subpd, noexec, exec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rsub1 			; Loop til done
	jmp	schunkdn		; Jump to chunk done code
isub0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
isub2:	xnorm_op_2d subpd, exec, exec, saved_reg ; Subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	isub2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	isub1			; Loop til done

	;; Chunk done

schunkdn:bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	sub0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest

	;; Block done

sblkdn:	sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	add	rsi, pass1blkdst
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	sskip2			; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 128		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
sskip2:	dec	loopcount2		; Decrement outer loop counter
	jnz	sblk 			; Loop til done

	;; Section done

	mov	rax, norm_ptr1		; Reload section start ptr
	mov	rbx, norm_ptr2		; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	ssec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_2			; Add the carries back in

	ad_epilog SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxsub2 ENDP

;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwxaddsubq2
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	ebx, addcount1		; Load blk count
uaddsub0:mov	eax, normval4		; Load count of 8KB chunks in a block
uaddsublp:
	xload	xmm0, [rcx]		; Load first number
	xcopy	xmm1, xmm0		; Dup first number
	addpd	xmm0, [rdx]		; Add in second number
	subpd	xmm1, [rdx]		; Subtract out second number
	xload	xmm2, [rcx+16]		; Load first number
	xcopy	xmm3, xmm2		; Dup first number
	addpd	xmm2, [rdx+16]		; Add in second number
	subpd	xmm3, [rdx+16]		; Subtract out second number
	xload	xmm4, [rcx+32]		; Load first number
	xcopy	xmm5, xmm4		; Dup first number
	addpd	xmm4, [rdx+32]		; Add in second number
	subpd	xmm5, [rdx+32]		; Subtract out second number
	xload	xmm6, [rcx+48]		; Load first number
	xcopy	xmm7, xmm6		; Dup first number
	addpd	xmm6, [rdx+48]		; Add in second number
	subpd	xmm7, [rdx+48]		; Subtract out second number
	xstore	[rsi], xmm0		; Save result
	xstore	[rbp], xmm1		; Save result
	xstore	[rsi+16], xmm2		; Save result
	xstore	[rbp+16], xmm3		; Save result
	xstore	[rsi+32], xmm4		; Save result
	xstore	[rbp+32], xmm5		; Save result
	xstore	[rsi+48], xmm6		; Save result
	xstore	[rbp+48], xmm7		; Save result
	bump	rcx, 64			; Next source
	bump	rdx, 64			; Next source
	bump	rsi, 64			; Next dest
	bump	rbp, 64			; Next dest
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	uaddsublp		; Loop if necessary
	bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	bump	rbp, 128		; Skip 128 bytes every 8KB
	dec	rax			; Check middle loop counter
	jnz	uaddsublp		; Loop if necessary
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest
	dec	rbx			; Check loop counter
	jnz	uaddsub0		; Loop if necessary
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxaddsubq2 ENDP

;;
;; Add and subtract two numbers with carry propagation
;;

saved_dest1_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_dest2_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_grp_ptr	EQU	PPTR [rsp+first_local+2*SZPTR]
saved_reg	EQU	PPTR [rsp+first_local+3*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+4*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+4*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+4*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+4*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+4*SZPTR+16]

PROCFL	gwxaddsub2
	ad_prolog 4*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	rbx, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xload	xmm7, XMM_BIGVAL	; Init 8 carries
	xstore	XMM_TMP1, xmm7
	xstore	XMM_TMP2, xmm7
	xstore	XMM_TMP3, xmm7
	xstore	XMM_TMP4, xmm7
	xstore	XMM_TMP5, xmm7
	xstore	XMM_TMP6, xmm7
	xstore	XMM_TMP7, xmm7
	xstore	XMM_TMP8, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

assec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	saved_dest1_ptr, rsi	; Save section start address
	mov	saved_dest2_ptr, rbp	; Save section start address
	mov	saved_grp_ptr, rbx	; Save section group ttp address

	;; Do a block

asblk:	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rax, norm_col_mults	; Addr of the column multipliers

as0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2as			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2ias0	   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in an 8KB chunk
nb2ras1: xnorm_addsub_2d noexec, noexec, saved_reg ; Add & sub and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2ras1			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
nb2ias0:mov	saved_reg, rax
	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
	mov	rax, saved_reg
nb2ias1:mov	saved_reg, rax
	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	mov	rax, saved_reg
nb2ias2:xnorm_addsub_2d exec, noexec, saved_reg	; Add & subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2ias2			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2ias1			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
b2as:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	ias0	   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in an 8KB chunk
ras1:	xnorm_addsub_2d noexec, exec, saved_reg ; Add & sub and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	ras1 			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
ias0:	mov	saved_reg, rax
	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
	mov	rax, saved_reg
ias1:	mov	saved_reg, rax
	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
	mov	rax, saved_reg
ias2:	xnorm_addsub_2d exec, exec, saved_reg ; Add & subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	ias2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	ias1			; Loop til done

	;; Chunk done

aschunkdn:bump	rcx, 128		; Skip 128 bytes every 8KB
	bump	rdx, 128		; Skip 128 bytes every 8KB
	bump	rsi, 128		; Skip 128 bytes every 8KB
	bump	rbp, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	as0
	add	rcx, pass2gapsize	; Next source
	add	rdx, pass2gapsize	; Next source
	add	rsi, pass2gapsize	; Next dest
	add	rbp, pass2gapsize	; Next dest

	;; Block done

asblkdn:sub	rsi, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rsi, rbx, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	sub	rbp, pass1blkdst	; Restore start of block ptr
					; Add 4 carries to start of block
	xnorm_op_2d_blk rbp, rbx, XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	add	rsi, pass1blkdst
	add	rbp, pass1blkdst
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	asskip2			; Yes, skip bumping ttp/biglit ptrs
	bump	rbx, 128		; Next set of group multipliers
	add	rdi, normval3		; Adjust little/big flags ptr
asskip2:dec	loopcount2		; Decrement outer loop counter
	jnz	asblk 			; Loop til done

	;; Section done

	mov	saved_reg, rbx
	mov	rax, saved_dest1_ptr	; Reload section start ptr for dest #1
	mov	rbx, saved_grp_ptr	; Reload section group ttp ptr
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4
	mov	rax, saved_dest2_ptr	; Reload section start ptr for dest #2
					; Add 2 carries to start of section
	xnorm_op_2d_sec XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8
	mov	rbx, saved_reg
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	assec			; Do next section if necessary

	;; All sections done

	xload	xmm6, XMM_TMP5		; Load non-wraparound carry
	xstore	XMM_TMP8, xmm6		; Save carry, final_carries_2 destroys XMM1-6

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_2			; Add the carries back in

	mov	rsi, DEST2ARG		; Addr of FFT data
	xload	xmm6, XMM_TMP8		; Load non-wraparound carry
	xload	xmm7, XMM_TMP7		; Load wraparound carry
	final_carries_2			; Add the carries back in

	ad_epilog 4*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxaddsub2 ENDP

;;
;; Copy one number and zero some low order words.
;;

PROCFL	gwxcopyzero2
	ad_prolog 0,0,rbx,rsi,rdi
	mov	rsi, SRCARG		; Address of first number
	mov	rdi, DESTARG		; Address of destination
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	ebx, addcount1		; Get number of blocks
cz1:	mov	eax, normval4		; Load count of 8KB chunks in a block
cz2:	xcopyzero
	bump	rsi, 64			; Next source
	bump	rdi, 64			; Next dest
	bump	rcx, 64			; Next compare offset
	add	eax, 80000000h/64	; 128 cache lines in a 8KB chunk
	jnc	short cz2		; Loop if necessary
	bump	rsi, 128		; Skip 128 bytes every 8KB
	bump	rdi, 128		; Skip 128 bytes every 8KB
	bump	rcx, 128		; Skip 128 bytes every 8KB
	dec	rax			; Test loop counter
	jnz	cz2			; Loop if necessary
	add	rsi, pass2gapsize
	add	rdi, pass2gapsize
	add	rcx, pass2gapsize
	dec	rbx			; Test loop counter
	jnz	cz1			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi
gwxcopyzero2 ENDP

;;
;; Add in a small number with carry propagation
;;

PROCFL	gwxadds2
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	movsd	xmm7, DBLARG		; Small addin value

	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rbx, norm_col_mults	; Addr of the column multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2addsm			; yes, do simpler rounding
	xnorm_smalladd_2d noexec	; Similar to add last carry code
	jmp	addsmdn
b2addsm:xnorm_smalladd_2d exec		; Similar to add last carry code
addsmdn:
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadds2 ENDP

;;
;; Multiply a number by a small value with carry propagation
;;

saved_sec_biglit EQU	PPTR [rsp+first_local+0*SZPTR]
saved_blk_start EQU	PPTR [rsp+first_local+1*SZPTR]
saved_blk_biglit EQU	PPTR [rsp+first_local+2*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+3*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+3*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+3*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+3*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+3*SZPTR+16]

PROCFL	gwxmuls2
	ad_prolog 3*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movlpd	xmm0, DBLARG		; Small multiplier value
	movhpd	xmm0, DBLARG		; Small multiplier value
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	jne	skip			; No, skip mul by two-to-phi fudge factor
	mulpd	xmm0, XMM_NORM012_FF	; Mul by FFTLEN/2
skip:	xstore	XMM_TMP5, xmm0		; Save small value * FFTLEN/2
	xload	xmm7, XMM_BIGVAL	; Init 4 carries
	xstore	XMM_TMP1, xmm7
	xstore	XMM_TMP2, xmm7
	xstore	XMM_TMP3, xmm7
	xstore	XMM_TMP4, xmm7
	mov	eax, count3		; Load 3 section counts

	;; Do a section

msec:	mov	loopcount1, eax		; Save section counts
	and	eax, 07FFh		; Form block count for this section
	mov	loopcount2, eax		; Save count of blocks in this section
	mov	norm_ptr1, rsi		; Save section start address
	mov	norm_ptr2, rbp		; Save section group ttp address
	mov	saved_sec_biglit, rdi	; Save section biglit ptr

	;; Do a block

mblk:	mov	saved_blk_start, rsi
	mov	saved_blk_biglit, rdi
	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax
	mov	rbx, norm_col_mults	; Addr of the column multipliers

mul0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mul			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2imul0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
nb2rmul1:xnorm_smallmul_2d noexec, noexec ; Mul and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2rmul1		; Loop til done
	jmp	mchunkdn		; Jump to chunk done code
nb2imul0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2imul1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
nb2imul2:xnorm_smallmul_2d exec, noexec ; Mul and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2imul2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2imul1		; Loop til done
	jmp	mchunkdn		; Jump to chunk done code

b2mul:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	imul0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
rmul1:	xnorm_smallmul_2d noexec, exec	; Mul and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rmul1 			; Loop til done
	jmp	mchunkdn		; Jump to chunk done code
imul0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
imul1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
imul2:	xnorm_smallmul_2d exec, exec	; Mul and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	imul2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	imul1			; Loop til done

	;; Chunk done

mchunkdn:bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	mul0

	;; Block done

	xchg	rsi, saved_blk_start	; Restore/save block start ptr
	xchg	rdi, saved_blk_biglit	; Restore/save block biglit ptr
	mov	rbx, norm_col_mults	; Addr of the column multipliers
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mblk			; yes, do simpler rounding
	xnorm_smallmul_2d_blk noexec	; Add 4 carries to start of block
	jmp	mblkdn
b2mblk:	xnorm_smallmul_2d_blk exec	; Add 4 carries to start of block
mblkdn:	mov	rsi, saved_blk_start	; Restore start ptr
	mov	rdi, saved_blk_biglit	; Restore biglit ptr
	add	rsi, pass2gapsize	; Next dest
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	mskip2			; Yes, skip bumping ttp/biglit ptrs
	add	rdi, normval3		; Adjust little/big flags ptr
	bump	rbp, 128		; Next set of group multipliers
mskip2:	dec	loopcount2		; Decrement outer loop counter
	jnz	mblk 			; Loop til section done

	;; Section done

	xchg	rsi, norm_ptr1		; Restore/save section start ptr
	xchg	rbp, norm_ptr2		; Restore/save section group ttp ptr
	xchg	rdi, saved_sec_biglit	; Restore/save section biglit ptr
	mov	rbx, norm_col_mults	; Addr of the column multipliers
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2msec			; yes, do simpler rounding
	xnorm_smallmul_2d_sec noexec	; Add 2 carries to start of section
	jmp	msecdn
b2msec:	xnorm_smallmul_2d_sec exec	; Add 2 carries to start of section
msecdn:	mov	rsi, norm_ptr1		; Restore next section start ptr
	mov	rbp, norm_ptr2		; Restore next section group ttp ptr
	mov	rdi, saved_sec_biglit	; Restore/save section biglit ptr
	mov	eax, loopcount1		; Get list of counts
	shr	eax, 11			; Move counts list along
	jnz	msec			; Do next section if necessary

	;; All sections done

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_2			; Add the carries back in

	ad_epilog 3*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxmuls2 ENDP

_TEXT	ENDS
END
