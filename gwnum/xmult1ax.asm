; Copyright 2001-2016 Mersenne Research, Inc.  All rights reserved
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
INCLUDE memory.mac
INCLUDE xnormal.mac

; Internal routine to add in the wraparound carry
; I'd like to make this a subroutine, but it is too difficult
; to get push_amt correct.

final_carries_1 MACRO
	LOCAL	b2c, c1dn

	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array

	cmp	B_IS_2, 0		; Is this base 2?
	jne	b2c			; yes, do simpler rounding
	xnorm_smallmul_1d_cleanup noexec ; Non-base 2 carry propagation
	jmp	c1dn
b2c:	xnorm_smallmul_1d_cleanup exec	; Base 2 carry propagation
c1dn:
	ENDM

_TEXT SEGMENT

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwxaddq1
	ad_prolog 0,0,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
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
	sub	eax, 1			; Check loop counter
	jnz	short uaddlp		; Loop if necessary
	ad_epilog 0,0,rsi
gwxaddq1 ENDP

;;
;; Add two numbers with carry propagation
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]

PROCFL	gwxadd1
	ad_prolog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	xload	xmm2, XMM_BIGVAL	; Start process with no carry
	xcopy	xmm3, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nadd0:	mov	loopcount1, eax		; Save loop counter
	and	eax, 07FFh		; Grab 11 bits of the counter
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for xnorm_op_1d_mid
	mov	saved_col_ptr, rbp	; remember rbp for xnorm_op_1d_mid
	sub	rax, rax		; Clear big/lit flags
	sub	rbx, rbx

	cmp	B_IS_2, 0		; Is b = 2?
	jne	b2add			; Yes, do the easier normalization
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2naddlp		; Yes, use two-to-phi multipliers
nb2raddlp:
	xnorm_op_1d addpd, noexec, noexec ; Add and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nb2raddlp		; Loop til done
	jmp	nadddn			; Loop til done
nb2naddlp:
	xnorm_op_1d addpd, exec, noexec ; Add and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nb2naddlp		; Loop til done
	jmp	nadddn			; Loop til done
b2add:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	naddlp			; Yes, use two-to-phi multipliers
raddlp:	xnorm_op_1d addpd, noexec, exec	; Add and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	raddlp			; Loop til done
	jmp	nadddn			; Loop til done
naddlp:	xnorm_op_1d addpd, exec, exec	; Add and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	naddlp			; Loop til done

nadddn:	mov	rbx, saved_col_ptr	; Restore multipliers pointer
	mov	rax, saved_dest_ptr	; Restore dest pointer
	xnorm_op_1d_mid_cleanup		; Rotate carries and add in carries
	mov	eax, loopcount1		; Restore loop counter
	shr	eax, 11			; Get next loop amount
	jnz	nadd0

	mov	rsi, DESTARG		; Address of result
	final_carries_1

	ad_epilog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadd1 ENDP

;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwxsubq1
	ad_prolog 0,0,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
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
	sub	eax, 1			; Check loop counter
	jnz	short usublp		; Loop if necessary
	ad_epilog 0,0,rsi
gwxsubq1 ENDP

;;
;; Subtract two numbers with carry propagation
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]

PROCFL	gwxsub1
	ad_prolog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	xload	xmm2, XMM_BIGVAL	; Start process with no carry
	xcopy	xmm3, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nsub0:	mov	loopcount1, eax		; Save loop counter
	and	eax, 07FFh		; Grab 11 bits of the counter
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for xnorm_op_1d_mid
	mov	saved_col_ptr, rbp	; remember rbp for xnorm_op_1d_mid
	sub	rax, rax		; Clear big/lit flag
	sub	rbx, rbx

	cmp	B_IS_2, 0		; Is b = 2?
	jne	b2sub			; Yes, do the easier normalization
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2nsublp		; Yes, use two-to-phi multipliers
nb2rsublp:
	xnorm_op_1d subpd, noexec, noexec ; Subtract and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nb2rsublp		; Loop til done
	jmp	nsubdn			; Jump if done
nb2nsublp:
	xnorm_op_1d subpd, exec, noexec	; Subtract and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nb2nsublp		; Loop til done
	jmp	nsubdn			; Jump if done
b2sub:	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nsublp			; Yes, use two-to-phi multipliers
rsublp:	xnorm_op_1d subpd, noexec, exec	; Subtract and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	rsublp			; Loop til done
	jmp	nsubdn			; Jump if done
nsublp:	xnorm_op_1d subpd, exec, exec	; Subtract and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nsublp			; Loop til done

nsubdn:	mov	rbx, saved_col_ptr	; Restore multipliers pointer
	mov	rax, saved_dest_ptr	; Restore dest pointer
	xnorm_op_1d_mid_cleanup		; Rotate carries and add in carries
	mov	eax, loopcount1		; Restore loop counter
	shr	eax, 11			; Get next loop amount
	jnz	nsub0

	mov	rsi, DESTARG		; Address of result
	final_carries_1

	ad_epilog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxsub1 ENDP

;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwxaddsubq1
	ad_prolog 0,0,rbp,rsi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	eax, addcount1		; Load loop counter
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
	sub	eax, 1			; Check loop counter
	jnz	uaddsublp		; Loop if necessary
	ad_epilog 0,0,rbp,rsi,xmm6,xmm7
gwxaddsubq1 ENDP

;;
;; Add and subtract two numbers with carry propagation
;;

saved_dest1_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_dest2_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+2*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+3*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+3*SZPTR+4]

PROCFL	gwxaddsub1
	ad_prolog 3*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	xload	xmm2, XMM_BIGVAL	; Start process with no carry
	xcopy	xmm3, xmm2
	xcopy	xmm6, xmm2
	xcopy	xmm7, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	rbx, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
naddsub0:
	mov	loopcount1, eax		; Save loop counter
	and	eax, 07FFh		; Grab 11 bits of the counter
	mov	loopcount2, eax
	mov	saved_col_ptr, rbx	; remember rbx for xnorm_op_1d_mid
	mov	saved_dest2_ptr, rbp	; remember dest #2 for xnorm_op_1d_mid
	mov	saved_dest1_ptr, rsi	; remember dest #1 for xnorm_op_1d_mid
	sub	rax, rax		; Clear big/lit flag

	cmp	B_IS_2, 0		; Is b = 2?
	jne	b2addsub		; Yes, do the easier normalization
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2naddsublp		; Yes, use two-to-phi multipliers
nb2raddsublp:
	xnorm_addsub_1d noexec, noexec	; Add/sub and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nb2raddsublp		; Loop til done
	jmp	naddsubdn		; Jump if done
nb2naddsublp:
	xnorm_addsub_1d exec, noexec	; Add/sub and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nb2naddsublp		; Loop til done
	jmp	naddsubdn		; Jump if done
b2addsub:
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	naddsublp		; Yes, use two-to-phi multipliers
raddsublp:
	xnorm_addsub_1d noexec, exec	; Add/sub and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	raddsublp		; Loop til done
	jmp	naddsubdn		; Jump if done
naddsublp:
	xnorm_addsub_1d exec, exec	; Add/sub and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	naddsublp		; Loop til done

naddsubdn:
	xchg	rbx, saved_col_ptr	; Save/Restore multipliers pointer
					; Rotate carries and add in carries
	xnorm_addsub_1d_mid_cleanup saved_dest1_ptr, saved_dest2_ptr
	mov	rbx, saved_col_ptr	; Restore multipliers pointer
	mov	eax, loopcount1		; Restore loop counter
	shr	eax, 11			; Get next loop amount
	jnz	naddsub0

	xstore	XMM_TMP7, xmm6		; Save carry
	xstore	XMM_TMP8, xmm7		; Save carry

	mov	rsi, DESTARG		; Address of result #1
	final_carries_1

	xload	xmm2, XMM_TMP7		; Load carry
	xload	xmm3, XMM_TMP8		; Load carry
	mov	rsi, DEST2ARG		; Address of result #2
	final_carries_1

	ad_epilog 3*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxaddsub1 ENDP

;;
;; Copy one number and zero some low order words.
;;

PROCFL	gwxcopyzero1
	ad_prolog 0,0,rsi,rdi
	mov	rsi, SRCARG		; Address of first number
	mov	rdi, DESTARG		; Address of destination
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	eax, addcount1		; Load loop counter
cz1:	xcopyzero			; Copy/zero 8 values
	bump	rsi, 64			; Next source
	bump	rdi, 64			; Next dest
	bump	rcx, 64			; Next compare offset
	sub	eax, 1			; Test loop counter
	jnz	cz1			; Loop if necessary
	ad_epilog 0,0,rsi,rdi
gwxcopyzero1 ENDP

;;
;; Add in a small number with carry propagation
;;

PROCFL	gwxadds1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	movsd	xmm7, DBLARG		; Small addin value

	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	sub	rax, rax		; Clear big/lit flag
	cmp	B_IS_2, 0		; Is b = 2?
	jne	b2adds			; yes, do simpler normalization
	xnorm_smalladd_1d noexec
	jmp	addsdn
b2adds:	xnorm_smalladd_1d exec
addsdn:
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadds1 ENDP

;;
;; Multiply a number by a small value
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_biglit_ptr EQU	PPTR [rsp+first_local+2*SZPTR]

PROCFL	gwxmuls1
	ad_prolog 3*SZPTR,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	movlpd	xmm7, DBLARG		; Small multiplier value
	movhpd	xmm7, DBLARG		; Small multiplier value
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	jne	skip			; No, skip mul by two-to-phi fudge factor
	mulpd	xmm7, XMM_NORM012_FF	; Mul by FFTLEN/2
skip:	xload	xmm2, XMM_BIGVAL	; Start process with no carry
	xcopy	xmm3, xmm2
	mov	edx, normcount1		; Load loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	saved_dest_ptr, rsi	; remember rsi for xnorm_smallmul_1d_mid
	mov	saved_col_ptr, rbp	; remember rbp for xnorm_smallmul_1d_mid
	mov	saved_biglit_ptr, rdi	; remember rdi for xnorm_smallmul_1d_mid
	sub	rax, rax		; Clear big/lit flags
	sub	rcx, rcx

	cmp	B_IS_2, 0		; Is b = 2?
	jne	nmul0			; yes, do simpler normalization

nb2nmul0:
	mov	ebx, edx		; Save loop counter
	and	ebx, 07FFh		; Grab 11 bits of the counter
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nb2nmullp		; Yes, use two-to-phi multipliers
nb2rmullp:
	xnorm_smallmul_1d noexec, noexec ; Mul and normalize 8 values
	sub	ebx, 1			; Decrement loop counter
	jnz	nb2rmullp		; Loop til done
	jmp	nb2nmuldn		; Rejoin common code
nb2nmullp:
	xnorm_smallmul_1d exec, noexec	; Mul and normalize 8 values
	sub	ebx, 1			; Decrement loop counter
	jnz	nb2nmullp		; Loop til done
nb2nmuldn:
	xchg	rsi, saved_dest_ptr	; Restore dest pointer
	xchg	rdi, saved_biglit_ptr	; Restore biglit pointer
	xchg	rbp, saved_col_ptr	; Restore multipliers pointer
	xnorm_smallmul_1d_mid_cleanup noexec ; Rotate carries and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	rdi, saved_biglit_ptr	; Restore biglit pointer
	mov	rbp, saved_col_ptr	; Restore multipliers pointer
	shr	edx, 11			; Get next loop amount
	jnz	nb2nmul0
	jmp	mulsdn

nmul0:	mov	ebx, edx		; Save loop counter
	and	ebx, 07FFh		; Grab 11 bits of the counter
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nmullp			; Yes, use two-to-phi multipliers
rmullp:	xnorm_smallmul_1d noexec, exec	; Mul and normalize 8 values
	sub	ebx, 1			; Decrement loop counter
	jnz	rmullp			; Loop til done
	jmp	nmuldn			; Rejoin common code
nmullp:	xnorm_smallmul_1d exec, exec	; Mul and normalize 8 values
	sub	ebx, 1			; Decrement loop counter
	jnz	nmullp			; Loop til done
nmuldn:	xchg	rsi, saved_dest_ptr	; Restore dest pointer
	xchg	rdi, saved_biglit_ptr	; Restore biglit pointer
	xchg	rbp, saved_col_ptr	; Restore multipliers pointer
	xnorm_smallmul_1d_mid_cleanup exec ; Rotate carries and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	rdi, saved_biglit_ptr	; Restore biglit pointer
	mov	rbp, saved_col_ptr	; Restore multipliers pointer
	shr	edx, 11			; Get next loop amount
	jnz	nmul0

mulsdn:	mov	rsi, DESTARG		; Address of result
	final_carries_1

	ad_epilog 3*SZPTR,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxmuls1 ENDP

;;
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_reg1	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_reg2	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_reg3	EQU	PPTR [rsp+first_local+2*SZPTR]

inorm	MACRO	lab, ttp, zero, echk, const, base2, sse4
	LOCAL	ilp0, ilp1
	PROCFLP	lab
	int_prolog 3*SZPTR,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	movsd	xmm0, Q [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, ADDIN_VALUE	;; Add in the requested value
no zero	movsd	Q [rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, ADDIN_VALUE	;; Do not include addin in sumout
	xload	xmm2, XMM_BIGVAL	;; Start process with no carry
	xcopy	xmm3, xmm2
	movlpd	xmm6, MAXERR		;; Current maximum error
	movhpd	xmm6, MAXERR
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
	mov	edx, normcount1		;; Load loop counter
	mov	saved_reg3, rsi		;; remember esi for xnorm012_1d_mid
ttp	mov	saved_reg2, rdi		;; remember edi for xnorm012_1d_mid
ttp	mov	saved_reg1, rbp		;; remember ebp for xnorm012_1d_mid

ilp0:	mov	ebx, edx		;; Load loop counter
	and	ebx, 07FFh		;; Grab 11 bits of the counter
ilp1:	xnorm_1d ttp, zero, echk, const, base2, sse4 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rbp, 128		;; Next set of 8 multipliers
ttp	bump	rdi, 4			;; Next big/little flags
	sub	ebx, 1			;; Test loop counter
	jnz	ilp1			;; Loop til done

	xchg	rsi, saved_reg3		;; Restore FFT data addr
ttp	xchg	rdi, saved_reg2		;; Restore big/lit pointer
ttp	xchg	rbp, saved_reg1		;; Restore ttp pointer
	xnorm012_1d_mid ttp, zero, base2 ;; Rotate carries and add in carries
	mov	rsi, saved_reg3		;; Restore FFT data addr
ttp	mov	rdi, saved_reg2		;; Restore big/lit pointer
ttp	mov	rbp, saved_reg1		;; Restore ttp pointer

	shr	edx, 11			;; Get next loop amount
	jnz	ilp0
no base2 jmp	non2dn			;; Go to non-base2 end code
zero	jmp	zdn			;; Go to zero upper half end code
base2 no zero jmp idn			;; Go to normal end code
	ENDPP lab
	ENDM

zpnorm	MACRO	lab, ttp, echk, const, base2, sse4, khi, c1, cm1
	LOCAL	ilp0, ilp1
	PROCFLP	lab
	int_prolog 3*SZPTR,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
	xload	xmm2, XMM_BIGVAL	;; Start process with no carry
	subpd	xmm3, xmm3
	movlpd	xmm6, MAXERR		;; Current maximum error
	movhpd	xmm6, MAXERR
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	sub	rax, rax		;; Clear big/lit flag
	mov	ebx, normcount1		;; Load loop counter
ilp0:	mov	saved_reg1, rbx		;; Save loop counter
	and	ebx, 07FFh		;; Grab 11 bits of the counter
	mov	saved_reg2, rdi		;; remember edi for xnorm012_1d_mid
	mov	saved_reg3, rsi		;; remember esi for xnorm012_1d_mid
	mov	rdx, rbp		;; remember ebp for xnorm012_1d_mid
ilp1:	xnorm_1d_zpad ttp, echk, const, base2, sse4, khi, c1, cm1 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rbp, 128		;; Next set of 8 multipliers
ttp	bump	rdi, 4			;; Next big/little flags
	sub	ebx, 1			;; Test loop counter
	jnz	ilp1			;; Loop til done
	mov	rbx, saved_reg3		;; Restore FFT data addr
	xchg	rdi, saved_reg2		;; Restore big/lit pointer
	xnorm012_1d_mid_zpad const, base2 ;; Rotate carries and add in carries
	mov	rdi, saved_reg2		;; Restore big/lit pointer
	mov	rbx, saved_reg1		;; Restore loop counter
	shr	ebx, 11			;; Get next loop amount
	jnz	ilp0
no base2 const jmp non2zpcdn		;; Go to zero padded FFT end code
no base2 no const jmp non2zpdn		;; Go to non-base2 end code
base2 const jmp	zpcdn			;; Go to zero padded FFT end code
base2 no const jmp zpdn			;; Go to zero padded FFT end code
	ENDPP lab
	ENDM

; The many different normalization routines.  One for each valid combination of
; rational/irrational, zeroing/no zeroing, error check/no error check,
; mul by const/no mul by const, base2 / other than base 2

	inorm	xr1, noexec, noexec, noexec, noexec, exec, noexec
	inorm	xr1e, noexec, noexec, exec, noexec, exec, noexec
	inorm	xr1c, noexec, noexec, noexec, exec, exec, noexec
	inorm	xr1ec, noexec, noexec, exec, exec, exec, noexec
	inorm	xr1z, noexec, exec, noexec, noexec, exec, noexec
	inorm	xr1ze, noexec, exec, exec, noexec, exec, noexec
	inorm	xi1, exec, noexec, noexec, noexec, exec, noexec
	inorm	xi1e, exec, noexec, exec, noexec, exec, noexec
	inorm	xi1c, exec, noexec, noexec, exec, exec, noexec
	inorm	xi1ec, exec, noexec, exec, exec, exec, noexec
	inorm	xi1z, exec, exec, noexec, noexec, exec, noexec
	inorm	xi1ze, exec, exec, exec, noexec, exec, noexec
	zpnorm	xr1zp, noexec, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xr1zpc1, noexec, noexec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xr1zpcm1, noexec, noexec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xr1zpe, noexec, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xr1zpec1, noexec, exec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xr1zpecm1, noexec, exec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xr1zpc, noexec, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xr1zpec, noexec, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xi1zp, exec, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xi1zpc1, exec, noexec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xi1zpcm1, exec, noexec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xi1zpe, exec, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xi1zpec1, exec, exec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xi1zpecm1, exec, exec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xi1zpc, exec, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xi1zpec, exec, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xr1zpk, noexec, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr1zpkc1, noexec, noexec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xr1zpkcm1, noexec, noexec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xr1zpek, noexec, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr1zpekc1, noexec, exec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xr1zpekcm1, noexec, exec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xr1zpck, noexec, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr1zpeck, noexec, exec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpk, exec, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpkc1, exec, noexec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xi1zpkcm1, exec, noexec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xi1zpek, exec, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpekc1, exec, exec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xi1zpekcm1, exec, exec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xi1zpck, exec, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpeck, exec, exec, exec, exec, noexec, noexec, noexec, noexec

	inorm	xr1b, noexec, noexec, noexec, noexec, noexec, noexec
	inorm	xr1eb, noexec, noexec, exec, noexec, noexec, noexec
	inorm	xr1cb, noexec, noexec, noexec, exec, noexec, noexec
	inorm	xr1ecb, noexec, noexec, exec, exec, noexec, noexec
	inorm	xi1b, exec, noexec, noexec, noexec, noexec, noexec
	inorm	xi1eb, exec, noexec, exec, noexec, noexec, noexec
	inorm	xi1cb, exec, noexec, noexec, exec, noexec, noexec
	inorm	xi1ecb, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr1zpb, noexec, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr1zpbc1, noexec, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xr1zpbcm1, noexec, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xr1zpeb, noexec, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr1zpebc1, noexec, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xr1zpebcm1, noexec, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xr1zpcb, noexec, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr1zpecb, noexec, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi1zpb, exec, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi1zpbc1, exec, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xi1zpbcm1, exec, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xi1zpeb, exec, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi1zpebc1, exec, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xi1zpebcm1, exec, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xi1zpcb, exec, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi1zpecb, exec, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr1zpbk, noexec, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr1zpbkc1, noexec, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xr1zpbkcm1, noexec, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xr1zpebk, noexec, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr1zpebkc1, noexec, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xr1zpebkcm1, noexec, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xr1zpcbk, noexec, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr1zpecbk, noexec, exec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpbk, exec, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpbkc1, exec, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xi1zpbkcm1, exec, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xi1zpebk, exec, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpebkc1, exec, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xi1zpebkcm1, exec, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xi1zpcbk, exec, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi1zpecbk, exec, exec, exec, noexec, noexec, noexec, noexec, noexec

	inorm	xr1s4, noexec, noexec, noexec, noexec, exec, exec
	inorm	xr1es4, noexec, noexec, exec, noexec, exec, exec
	inorm	xr1cs4, noexec, noexec, noexec, exec, exec, exec
	inorm	xr1ecs4, noexec, noexec, exec, exec, exec, exec
	inorm	xi1s4, exec, noexec, noexec, noexec, exec, exec
	inorm	xi1es4, exec, noexec, exec, noexec, exec, exec
	inorm	xi1cs4, exec, noexec, noexec, exec, exec, exec
	inorm	xi1ecs4, exec, noexec, exec, exec, exec, exec
	zpnorm	xr1zps4, noexec, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xr1zps4c1, noexec, noexec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xr1zps4cm1, noexec, noexec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xr1zpes4, noexec, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xr1zpes4c1, noexec, exec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xr1zpes4cm1, noexec, exec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xr1zpcs4, noexec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xr1zpecs4, noexec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xi1zps4, exec, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xi1zps4c1, exec, noexec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xi1zps4cm1, exec, noexec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xi1zpes4, exec, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xi1zpes4c1, exec, exec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xi1zpes4cm1, exec, exec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xi1zpcs4, exec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xi1zpecs4, exec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xr1zps4k, noexec, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xr1zps4kc1, noexec, noexec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xr1zps4kcm1, noexec, noexec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xr1zpes4k, noexec, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xr1zpes4kc1, noexec, exec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xr1zpes4kcm1, noexec, exec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xr1zpcs4k, noexec, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xr1zpecs4k, noexec, exec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xi1zps4k, exec, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xi1zps4kc1, exec, noexec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xi1zps4kcm1, exec, noexec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xi1zpes4k, exec, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xi1zpes4kc1, exec, exec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xi1zpes4kcm1, exec, exec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xi1zpcs4k, exec, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xi1zpecs4k, exec, exec, exec, exec, exec, noexec, noexec, noexec

	inorm	xr1bs4, noexec, noexec, noexec, noexec, noexec, exec
	inorm	xr1ebs4, noexec, noexec, exec, noexec, noexec, exec
	inorm	xr1cbs4, noexec, noexec, noexec, exec, noexec, exec
	inorm	xr1ecbs4, noexec, noexec, exec, exec, noexec, exec
	inorm	xi1bs4, exec, noexec, noexec, noexec, noexec, exec
	inorm	xi1ebs4, exec, noexec, exec, noexec, noexec, exec
	inorm	xi1cbs4, exec, noexec, noexec, exec, noexec, exec
	inorm	xi1ecbs4, exec, noexec, exec, exec, noexec, exec
	zpnorm	xr1zpbs4, noexec, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xr1zpbs4c1, noexec, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xr1zpbs4cm1, noexec, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xr1zpebs4, noexec, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xr1zpebs4c1, noexec, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xr1zpebs4cm1, noexec, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xr1zpcbs4, noexec, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr1zpecbs4, noexec, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xi1zpbs4, exec, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xi1zpbs4c1, exec, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xi1zpbs4cm1, exec, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xi1zpebs4, exec, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xi1zpebs4c1, exec, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xi1zpebs4cm1, exec, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xi1zpcbs4, exec, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xi1zpecbs4, exec, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr1zpbs4k, noexec, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr1zpbs4kc1, noexec, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xr1zpbs4kcm1, noexec, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xr1zpebs4k, noexec, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr1zpebs4kc1, noexec, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xr1zpebs4kcm1, noexec, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xr1zpcbs4k, noexec, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr1zpecbs4k, noexec, exec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi1zpbs4k, exec, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi1zpbs4kc1, exec, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xi1zpbs4kcm1, exec, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xi1zpebs4k, exec, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi1zpebs4kc1, exec, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xi1zpebs4kcm1, exec, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xi1zpcbs4k, exec, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi1zpecbs4k, exec, exec, exec, noexec, exec, noexec, noexec, noexec

; Common code to finish off the one-pass FFTs normalization.  The
; Windows 64-bit ABI frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.

PROCF	__common_xnorm1_end_code

	;; Dummy prolog to match normalization code
	int_prolog 3*SZPTR,0,0

; Finish off the normalization process by adding any carry to first values.
; Handle both the with and without two-to-phi array cases.

non2zpdn:mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d_zpad noexec, noexec	; Add in carries
	jmp	cmnend			; All done, go cleanup

non2zpcdn:mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d_zpad exec, noexec	; Add in carries
	jmp	cmnend			; All done, go cleanup

zpcdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d_zpad exec, exec	; Add in carries
	jmp	cmnend			; All done, go cleanup

zpdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d_zpad noexec, exec	; Add in carries
	jmp	cmnend			; All done, go cleanup

non2dn:	mov	rsi, DESTARG		; Address of squared number
	xnorm_top_carry_1d		; Adjust top carry when k > 1
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d noexec, noexec	; Add in carries
	jmp	cmnend			; All done, go cleanup

zdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d exec, exec		; Add in carries
	jmp	cmnend			; All done, go cleanup

idn:	mov	rsi, DESTARG		; Address of squared number
	xnorm_top_carry_1d		; Adjust top carry when k > 1
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d noexec, exec	; Add in carries

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

cmnend:	mov	rsi, DESTARG		; Address of squared number
	xstore	XMM_TMP1, xmm7		; Add together the two partial sumouts
	addsd	xmm7, Q XMM_TMP1+8
	mulsd	xmm7, ttmp_ff_inv
	movsd	Q [rsi-24], xmm7	; Save sum of FFT outputs
	xstore	XMM_TMP1, xmm6		; Compute new maximum error
	maxsd	xmm6, Q XMM_TMP1+8
	movsd	MAXERR, xmm6

; Return

	int_epilog 3*SZPTR,0,0
__common_xnorm1_end_code ENDP

_TEXT	ENDS
END
