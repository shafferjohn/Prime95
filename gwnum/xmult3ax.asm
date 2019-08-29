; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu

;
; Additional routines used with r4dwpn (r4delay with partial normalization) FFTs,
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

final_carries_3 MACRO
	LOCAL	b2c, zpc, b2zpc, c3dn

	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rbx, norm_col_mults	; Addr of the column multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array

	cmp	ZERO_PADDED_FFT, 0	;; Zero-padded FFT?
	jne	zpc			;; Yes, do special zpad carry

	xnorm_top_carry_cmn rsi, xmm7, 2
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2c			; yes, do simpler rounding
	xnorm_smallmul_wpn_fft noexec	; Add 2 carries to start of fft
	jmp	c3dn
b2c:	xnorm_smallmul_wpn_fft exec	; Add 2 carries to start of fft
	jmp	c3dn

zpc:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2zpc			; yes, do simpler rounding
	xnorm_smallmul_wpn_fft_zpad noexec ; Do the special zpad carry not base 2
	jmp	c3dn
b2zpc:	xnorm_smallmul_wpn_fft_zpad exec ; Do the special zpad carry base 2
c3dn:
	ENDM


_TEXT SEGMENT

;;
;; Add two numbers with carry propagation
;;

saved_blk_start EQU	PPTR [rsp+first_local+0*SZPTR]
saved_blk_biglit EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+2*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+2*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+2*SZPTR+16]

PROCFL	gwxadd3
	ad_prolog 2*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
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

	;; Do a block

	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount1, eax		; Save count
ablk0:	mov	eax, count2		; Load wpn count
	mov	loopcount2, eax		; Save count
ablk:	mov	saved_blk_start, rsi
	mov	saved_blk_biglit, rdi
	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

add0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2add			; yes, do simpler rounding
nb2add0:cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2iadd0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
nb2radd1:
	xnorm_op_wpn addpd, noexec, noexec ; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2radd1		; Loop til done
	jmp	achunkdn		; Jump to chunk done code
nb2iadd0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2iadd1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
nb2iadd2:xnorm_op_wpn addpd, exec, noexec ; Add and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2iadd2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2iadd1		; Loop til done
	jmp	achunkdn		; Jump to chunk done code
b2add:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	iadd0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
radd1:	xnorm_op_wpn addpd, noexec, exec ; Add and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	radd1 			; Loop til done
	jmp	achunkdn		; Jump to chunk done code
iadd0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
iadd1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
iadd2:	xnorm_op_wpn addpd, exec, exec	; Add and normalize 8 values
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

ablkdn:	xchg	rsi, saved_blk_start	; Restore/save block start ptr
	xchg	rdi, saved_blk_biglit	; Restore/save block biglit ptr
	xnorm_op_wpn_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4 ; Add 4 carries to start of block
	mov	rsi, saved_blk_start	; Restore start ptr
	mov	rdi, saved_blk_biglit	; Restore biglit ptr
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short askip2		; Yes, skip bumping ttp/biglit ptrs
	add	rdi, normval3		; Adjust little/big flags ptr
askip2:	sub	loopcount2, 1		; Test loop counter
	jnz	ablk
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short askip3		; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 4*XMM_GMD		; Next set of group multipliers
askip3:	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	ablk0 			; Loop til done

	;; All blocks done

	mov	rax, DESTARG		; Address of destination
	mov	rbx, norm_grp_mults	; Group ttp ptr
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xnorm_op_wpn_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4 ; Add 2 carries to start of section

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_3			; Add the carries back in

	ad_epilog 2*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadd3 ENDP

;;
;; Subtract two numbers with carry propagation
;;

saved_blk_start EQU	PPTR [rsp+first_local+0*SZPTR]
saved_blk_biglit EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+2*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+2*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+2*SZPTR+16]

PROCFL	gwxsub3
	ad_prolog 2*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
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

	;; Do a block

	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount1, eax		; Save count
sblk0:	mov	eax, count2		; Load wpn count
	mov	loopcount2, eax		; Save count
sblk:	mov	saved_blk_start, rsi
	mov	saved_blk_biglit, rdi
	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

sub0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2sub			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2isub0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
nb2rsub1:
	xnorm_op_wpn subpd, noexec, noexec ; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2rsub1		; Loop til done
	jmp	schunkdn		; Jump to chunk done code
nb2isub0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2isub1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
nb2isub2:xnorm_op_wpn subpd, exec, noexec ; Subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2isub2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2isub1		; Loop til done
	jmp	schunkdn		; Jump to chunk done code
b2sub:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	isub0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
rsub1:	xnorm_op_wpn subpd, noexec, exec ; Subtract and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rsub1 			; Loop til done
	jmp	schunkdn		; Jump to chunk done code
isub0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
isub1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
isub2:	xnorm_op_wpn subpd, exec, exec	; Subtract and normalize 8 values
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

sblkdn:	xchg	rsi, saved_blk_start	; Restore/save block start ptr
	xchg	rdi, saved_blk_biglit	; Restore/save block biglit ptr
	xnorm_op_wpn_blk rsi, rbp, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4 ; Add 4 carries to start of block
	mov	rsi, saved_blk_start	; Restore start ptr
	mov	rdi, saved_blk_biglit	; Restore biglit ptr
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short sskip2		; Yes, skip bumping ttp/biglit ptrs
	add	rdi, normval3		; Adjust little/big flags ptr
sskip2:	sub	loopcount2, 1		; Test loop counter
	jnz	sblk
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short sskip3		; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 4*XMM_GMD		; Next set of group multipliers
sskip3:	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	sblk0 			; Loop til done

	;; All blocks done

	mov	rax, DESTARG		; Reload FFT data ptr
	mov	rbx, norm_grp_mults	; Reload group ttp ptr
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xnorm_op_wpn_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4 ; Add 2 carries to start of section

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_3			; Add the carries back in

	ad_epilog 2*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxsub3 ENDP

;;
;; Add and subtract two numbers with carry propagation
;;

saved_blk_start EQU	PPTR [rsp+first_local+0*SZPTR]
saved_blk_start2 EQU	PPTR [rsp+first_local+1*SZPTR]
saved_blk_biglit EQU	PPTR [rsp+first_local+2*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+3*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+3*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+3*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+3*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+3*SZPTR+16]

PROCFL	gwxaddsub3
	ad_prolog 3*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
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

	;; Do a block

	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount1, eax		; Save count
asblk0:	mov	eax, count2		; Load wpn count
	mov	loopcount2, eax		; Save count
asblk:	mov	saved_blk_start, rsi
	mov	saved_blk_start2, rbp
	mov	saved_blk_biglit, rdi
	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

as0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2as			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2ias0	   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in an 8KB chunk
nb2ras1: xnorm_addsub_wpn noexec, noexec ; Add & sub and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2ras1			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
nb2ias0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2ias1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
nb2ias2:xnorm_addsub_wpn exec, noexec	; Add & subtract and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2ias2			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2ias1			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
b2as:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	ias0	   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in an 8KB chunk
ras1:	xnorm_addsub_wpn noexec, exec	; Add & sub and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	ras1 			; Loop til done
	jmp	aschunkdn		; Jump to chunk done code
ias0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
ias1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
ias2:	xnorm_addsub_wpn exec, exec	; Add & subtract and normalize 8 values
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

asblkdn:xchg	rsi, saved_blk_start	; Restore/save block start ptr
	xchg	rbp, saved_blk_start2	; Restore/save block start ptr
	xchg	rdi, saved_blk_biglit	; Restore/save block biglit ptr
	xnorm_op_wpn_blk rsi, rbx, XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4 ; Add 4 carries to start of block
	xnorm_op_wpn_blk rbp, rbx, XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8 ; Add 4 carries to start of block
	mov	rsi, saved_blk_start	; Restore blk ptr
	mov	rbp, saved_blk_start2	; Restore blk ptr
	mov	rdi, saved_blk_biglit	; Restore biglit ptr
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short asskip2		; Yes, skip bumping ttp/biglit ptrs
	add	rdi, normval3		; Adjust little/big flags ptr
asskip2:sub	loopcount2, 1		; Test loop counter
	jnz	asblk
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short asskip3		; Yes, skip bumping ttp/biglit ptrs
	bump	rbx, 4*XMM_GMD		; Next set of group multipliers
asskip3:sub	loopcount1, 1		; Decrement outer loop counter
	jnz	asblk0 			; Loop til done

	;; All blocks done

	mov	rax, DESTARG		; Reload ptr for dest #1
	mov	rbx, norm_grp_mults	; Reload group ttp ptr
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xnorm_op_wpn_sec XMM_TMP1, XMM_TMP2, XMM_TMP3, XMM_TMP4 ; Add 2 carries to start of section
	mov	rax, DEST2ARG		; Reload ptr for dest #2
	xnorm_op_wpn_sec XMM_TMP5, XMM_TMP6, XMM_TMP7, XMM_TMP8 ; Add 2 carries to start of section

	;; All sections done

	xload	xmm6, XMM_TMP5		; Load non-wraparound carry
	xstore	XMM_TMP8, xmm6		; Save carry, final_carries_3 destroys XMM1-6

	mov	rsi, DESTARG		; Addr of FFT data
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_3			; Add the carries back in

	mov	rsi, DEST2ARG		; Addr of FFT data
	xload	xmm6, XMM_TMP8		; Load non-wraparound carry
	xload	xmm7, XMM_TMP7		; Load wraparound carry
	final_carries_3			; Add the carries back in

	ad_epilog 3*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxaddsub3 ENDP

;;
;; Add in a small number with carry propagation
;;

PROCFL	gwxadds3
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	movsd	xmm7, DBLARG		; Small addin value

	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2addsm			; yes, do simpler rounding
	xnorm_smalladd_wpn noexec	; Similar to add last carry code
	jmp	addsmdn
b2addsm:xnorm_smalladd_wpn exec		; Similar to add last carry code
addsmdn:
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadds3 ENDP

;;
;; Multiply a number by a small value with carry propagation
;;

saved_blk_start EQU	PPTR [rsp+first_local+0*SZPTR]
saved_blk_biglit EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+2*SZPTR+8]
loopcount4	EQU	DPTR [rsp+first_local+2*SZPTR+12]
loopcount5	EQU	DPTR [rsp+first_local+2*SZPTR+16]

PROCFL	gwxmuls3
	ad_prolog 2*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, norm_grp_mults	; Addr of the group multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	movlpd	xmm0, DBLARG		; Small multiplier value
	movhpd	xmm0, DBLARG		; Small multiplier value
	xstore	XMM_TMP5, xmm0		; Save small value
	xload	xmm7, XMM_BIGVAL	; Init 4 carries
	xstore	XMM_TMP1, xmm7
	xstore	XMM_TMP2, xmm7
	xstore	XMM_TMP3, xmm7
	xstore	XMM_TMP4, xmm7

	;; Do a block

	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount1, eax		; Save count
mblk0:	mov	eax, count2		; Load wpn count
	mov	loopcount2, eax		; Save count
mblk:	mov	saved_blk_start, rsi
	mov	saved_blk_biglit, rdi
	mov	eax, normval4		; Load count of 8KB chunks in a blk
	mov	loopcount3, eax

mul0:	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mul			; yes, do simpler rounding
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2imul0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
nb2rmul1:xnorm_smallmul_wpn noexec, noexec ; Mul and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	nb2rmul1		; Loop til done
	jmp	mchunkdn		; Jump to chunk done code
nb2imul0:mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
nb2imul1:mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
nb2imul2:xnorm_smallmul_wpn exec, noexec ; Mul and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	nb2imul2		; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	nb2imul1		; Loop til done
	jmp	mchunkdn		; Jump to chunk done code

b2mul:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	imul0   		; Yes, use two-to-phi multipliers
	mov	loopcount4, 128		; Cache lines in 8KB chunk
rmul1:	xnorm_smallmul_wpn noexec, exec	; Mul and normalize 8 values
	dec	loopcount4		; Decrement inner loop counter
	jnz	rmul1 			; Loop til done
	jmp	mchunkdn		; Jump to chunk done code
imul0:	mov	eax, normval1		; Load count of clms in 8KB chunk
	mov	loopcount4, eax
imul1:	mov	eax, cache_line_multiplier ; Load inner loop count
	mov	loopcount5, eax		; Save inner loop count
imul2:	xnorm_smallmul_wpn exec, exec	; Mul and normalize 8 values
	dec	loopcount5		; Decrement inner loop counter
	jnz	imul2 			; Loop til done
	add	rdi, normval2		; Adjust ptr to little/big flags
	dec	loopcount4		; Decrement middle loop counter
	jnz	imul1			; Loop til done

	;; Chunk done

mchunkdn:bump	rsi, 128		; Skip 128 bytes every 8KB
	dec	loopcount3		; Test loop counter
	jnz	mul0
	add	rsi, pass2gapsize	; Next dest

	;; Block done

	xchg	rsi, saved_blk_start	; Restore/save block start ptr
	xchg	rdi, saved_blk_biglit	; Restore/save block biglit ptr
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2mblk			; yes, do simpler rounding
	xnorm_smallmul_wpn_blk noexec	; Add 4 carries to start of block
	jmp	mblkdn
b2mblk:	xnorm_smallmul_wpn_blk exec	; Add 4 carries to start of block
mblkdn:	mov	rsi, saved_blk_start	; Restore start ptr
	mov	rdi, saved_blk_biglit	; Restore biglit ptr
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short mskip2		; Yes, skip bumping ttp/biglit ptrs
	add	rdi, normval3		; Adjust little/big flags ptr
mskip2:	sub	loopcount2, 1		; Test loop counter
	jnz	mblk
	cmp	RATIONAL_FFT, 0		; Rational FFT?
	jne	short mskip3		; Yes, skip bumping ttp/biglit ptrs
	bump	rbp, 4*XMM_GMD		; Next set of group multipliers
mskip3:	sub	loopcount1, 1		; Decrement outer loop counter
	jnz	mblk0 			; Loop til section done

	;; Section done

	mov	rsi, DESTARG		; Restore data ptr
	mov	rbp, norm_grp_mults	; Restore group ttp ptr
	mov	rdi, norm_biglit_array	; Restore biglit ptr
	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2msec			; yes, do simpler rounding
	xnorm_smallmul_wpn_sec noexec	; Add 2 carries to start of section
	jmp	msecdn
b2msec:	xnorm_smallmul_wpn_sec exec	; Add 2 carries to start of section
msecdn:

	;; All sections done

	mov	rsi, DESTARG		; Restore data ptr
	xload	xmm6, XMM_TMP1		; Load non-wraparound carry
	xload	xmm7, XMM_TMP3		; Load wraparound carry
	final_carries_3			; Add the carries back in

	ad_epilog 2*SZPTR+20,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxmuls3 ENDP

_TEXT	ENDS
END
