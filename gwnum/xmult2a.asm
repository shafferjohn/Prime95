; Copyright 2001-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8 levels of two pass FFTs that use the
; SSE2 instructions.
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

_TEXT SEGMENT

;;
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_rsi	EQU	PPTR [rsp+first_local]
IFDEF X86_64
loopcount1	EQU	r10
ELSE
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
ENDIF
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]

inorm	MACRO	lab, ttp, zero, echk, const, base2, sse4
	LOCAL	noadd, setlp, ilp0, ilp1, ilexit, done
	PROCFLP	lab
	int_prolog SZPTR+12,0,0
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
no zero	mov	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	edx, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	movsd	xmm0, Q [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, ADDIN_VALUE	;; Add in the requested value
no zero	movsd	Q [rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, ADDIN_VALUE	;; Do not include addin in sumout
noadd:	mov	saved_rsi, rsi		;; Save for xtop_carry_adjust
	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ttp	mov	eax, cache_line_multiplier ;; Load inner loop counter
	lea	rdi, XMM_COL_MULTS[128]	;; Load col mult scratch area
setlp:	xnorm_2d_setup ttp
ttp	bump	rdi, 512		;; Next scratch area section
ttp	bump	rbx, 32			;; Next column multiplier
ttp	sub	al, 1			;; Each cache line has its own col mult
ttp	jnz	setlp

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
ttp	movzx	rax, BYTE PTR [rdi+0]	;; Load big vs. little flags
ttp	movzx	rcx, BYTE PTR [rdi+1]	;; Load big vs. little flags
no ttp	sub	rax, rax
no ttp	sub	rcx, rcx
IFDEF X86_64
ttp	movzx	r8, BYTE PTR [rdi+2]	;; Load big vs. little flags
ttp	movzx	r9, BYTE PTR [rdi+3]	;; Load big vs. little flags
no ttp	sub	r8, r8
no ttp	sub	r9, r9
ENDIF
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, rbx		;; Save loop counter
	lea	rbx, XMM_COL_MULTS	;; Load col mult scratch area
	L2prefetch128 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_2d ttp, zero, echk, const, base2, sse4 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rbx, 512		;; Next column multipliers
ttp	bump	rdi, 4			;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	bump	rbp, 64			;; Next set of carries
ttp	bump	rdx, 128		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error

	; Handle adjusting the carry out of the topmost FFT word

	mov	eax, THIS_BLOCK		;; Check for processing last block
	cmp	eax, LAST_PASS1_BLOCK
	jne	done			;; Jump if not last block
	mov	rsi, saved_rsi		;; Restore FFT data ptr
	xnorm_top_carry			;; Adjust carry if k > 1

done:	int_epilog SZPTR+12,0,0
	ENDPP	lab
	ENDM

IFDEF X86_64
loopcount1z	EQU	r10
ELSE
loopcount1z	EQU	DPTR [rsp+first_local]
ENDIF
loopcount2z	EQU	DPTR [rsp+first_local+4]
loopcount3z	EQU	DPTR [rsp+first_local+8]

zpnorm	MACRO	lab, ttp, echk, const, base2, sse4, khi, c1, cm1
	LOCAL	setlp, ilp0, ilp1, ilexit
	PROCFLP	lab
	int_prolog 12,0,0
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
echk	movapd	xmm6, XMM_MAXERR	;; Load maximum error

	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ttp	mov	eax, cache_line_multiplier ;; Load inner loop counter
	lea	rdi, XMM_COL_MULTS[128]	;; Load col mult scratch area
setlp:	xnorm_2d_setup ttp
ttp	bump	rdi, 512		;; Next scratch area section
ttp	bump	rbx, 32			;; Next column multiplier
ttp	sub	al, 1			;; Each cache line has its own col mult
ttp	jnz	setlp

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2z, eax	;; Save loop counter
	mov	loopcount3z, 0		;; Clear outermost loop counter
ttp	movzx	rax, BYTE PTR [rdi+0]	;; Load big vs. little flags
no ttp	sub	rax, rax
no ttp	sub	rcx, rcx
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1z, rbx	;; Save loop counter
	lea	rbx, XMM_COL_MULTS	;; Load col mult scratch area
IFDEF X86_64
	xload	xmm4, [rbp+0*16]	;; Preload carries
	xload	xmm12, [rbp+1*16]
	xload	xmm3, [rbp+2*16]
	xload	xmm11, [rbp+3*16]
ENDIF
	L2prefetch128 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_2d_zpad ttp, echk, const, base2, sse4, khi, c1, cm1 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rbx, 512		;; Next column multipliers
ttp	bump	rdi, 4			;; Next big/little flags
	sub	loopcount1z, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
IFDEF X86_64
	xstore	[rbp+0*16], xmm4	;; Store carries
	xstore	[rbp+1*16], xmm12
	xstore	[rbp+2*16], xmm3
	xstore	[rbp+3*16], xmm11
ENDIF
	bump	rbp, 64			;; Next set of carries
ttp	bump	rdx, 128		;; Next set of 8 group multipliers
	sub	loopcount2z, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3z, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
echk	movapd	XMM_MAXERR, xmm6	;; Save maximum error
	int_epilog 12,0,0
	ENDPP	lab
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	xr2, noexec, noexec, noexec, noexec, exec, noexec
	inorm	xr2e, noexec, noexec, exec, noexec, exec, noexec
	inorm	xr2c, noexec, noexec, noexec, exec, exec, noexec
	inorm	xr2ec, noexec, noexec, exec, exec, exec, noexec
	inorm	xr2z, noexec, exec, noexec, noexec, exec, noexec
	inorm	xr2ze, noexec, exec, exec, noexec, exec, noexec
	inorm	xi2, exec, noexec, noexec, noexec, exec, noexec
	inorm	xi2e, exec, noexec, exec, noexec, exec, noexec
	inorm	xi2c, exec, noexec, noexec, exec, exec, noexec
	inorm	xi2ec, exec, noexec, exec, exec, exec, noexec
	inorm	xi2z, exec, exec, noexec, noexec, exec, noexec
	inorm	xi2ze, exec, exec, exec, noexec, exec, noexec
	zpnorm	xr2zp, noexec, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xr2zpc1, noexec, noexec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xr2zpcm1, noexec, noexec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xr2zpe, noexec, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xr2zpec1, noexec, exec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xr2zpecm1, noexec, exec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xr2zpc, noexec, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xr2zpec, noexec, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xi2zp, exec, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xi2zpc1, exec, noexec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xi2zpcm1, exec, noexec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xi2zpe, exec, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	xi2zpec1, exec, exec, noexec, exec, noexec, exec, exec, noexec
	zpnorm	xi2zpecm1, exec, exec, noexec, exec, noexec, exec, noexec, exec
	zpnorm	xi2zpc, exec, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xi2zpec, exec, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	xr2zpk, noexec, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr2zpkc1, noexec, noexec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xr2zpkcm1, noexec, noexec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xr2zpek, noexec, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr2zpekc1, noexec, exec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xr2zpekcm1, noexec, exec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xr2zpck, noexec, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xr2zpeck, noexec, exec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpk, exec, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpkc1, exec, noexec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xi2zpkcm1, exec, noexec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xi2zpek, exec, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpekc1, exec, exec, noexec, exec, noexec, noexec, exec, noexec
	zpnorm	xi2zpekcm1, exec, exec, noexec, exec, noexec, noexec, noexec, exec
	zpnorm	xi2zpck, exec, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpeck, exec, exec, exec, exec, noexec, noexec, noexec, noexec

	inorm	xr2b, noexec, noexec, noexec, noexec, noexec, noexec
	inorm	xr2eb, noexec, noexec, exec, noexec, noexec, noexec
	inorm	xr2cb, noexec, noexec, noexec, exec, noexec, noexec
	inorm	xr2ecb, noexec, noexec, exec, exec, noexec, noexec
	inorm	xi2b, exec, noexec, noexec, noexec, noexec, noexec
	inorm	xi2eb, exec, noexec, exec, noexec, noexec, noexec
	inorm	xi2cb, exec, noexec, noexec, exec, noexec, noexec
	inorm	xi2ecb, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr2zpb, noexec, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr2zpbc1, noexec, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xr2zpbcm1, noexec, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xr2zpeb, noexec, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr2zpebc1, noexec, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xr2zpebcm1, noexec, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xr2zpcb, noexec, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr2zpecb, noexec, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi2zpb, exec, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi2zpbc1, exec, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xi2zpbcm1, exec, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xi2zpeb, exec, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi2zpebc1, exec, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	xi2zpebcm1, exec, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	xi2zpcb, exec, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xi2zpecb, exec, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	xr2zpbk, noexec, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr2zpbkc1, noexec, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xr2zpbkcm1, noexec, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xr2zpebk, noexec, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr2zpebkc1, noexec, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xr2zpebkcm1, noexec, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xr2zpcbk, noexec, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xr2zpecbk, noexec, exec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpbk, exec, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpbkc1, exec, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xi2zpbkcm1, exec, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xi2zpebk, exec, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpebkc1, exec, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	xi2zpebkcm1, exec, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	xi2zpcbk, exec, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	xi2zpecbk, exec, exec, exec, noexec, noexec, noexec, noexec, noexec

	inorm	xr2s4, noexec, noexec, noexec, noexec, exec, exec
	inorm	xr2es4, noexec, noexec, exec, noexec, exec, exec
	inorm	xr2cs4, noexec, noexec, noexec, exec, exec, exec
	inorm	xr2ecs4, noexec, noexec, exec, exec, exec, exec
	inorm	xi2s4, exec, noexec, noexec, noexec, exec, exec
	inorm	xi2es4, exec, noexec, exec, noexec, exec, exec
	inorm	xi2cs4, exec, noexec, noexec, exec, exec, exec
	inorm	xi2ecs4, exec, noexec, exec, exec, exec, exec
	zpnorm	xr2zps4, noexec, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xr2zps4c1, noexec, noexec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xr2zps4cm1, noexec, noexec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xr2zpes4, noexec, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xr2zpes4c1, noexec, exec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xr2zpes4cm1, noexec, exec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xr2zpcs4, noexec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xr2zpecs4, noexec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xi2zps4, exec, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xi2zps4c1, exec, noexec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xi2zps4cm1, exec, noexec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xi2zpes4, exec, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	xi2zpes4c1, exec, exec, noexec, exec, exec, exec, exec, noexec
	zpnorm	xi2zpes4cm1, exec, exec, noexec, exec, exec, exec, noexec, exec
	zpnorm	xi2zpcs4, exec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xi2zpecs4, exec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	xr2zps4k, noexec, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xr2zps4kc1, noexec, noexec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xr2zps4kcm1, noexec, noexec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xr2zpes4k, noexec, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xr2zpes4kc1, noexec, exec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xr2zpes4kcm1, noexec, exec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xr2zpcs4k, noexec, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xr2zpecs4k, noexec, exec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xi2zps4k, exec, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xi2zps4kc1, exec, noexec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xi2zps4kcm1, exec, noexec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xi2zpes4k, exec, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	xi2zpes4kc1, exec, exec, noexec, exec, exec, noexec, exec, noexec
	zpnorm	xi2zpes4kcm1, exec, exec, noexec, exec, exec, noexec, noexec, exec
	zpnorm	xi2zpcs4k, exec, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	xi2zpecs4k, exec, exec, exec, exec, exec, noexec, noexec, noexec

	inorm	xr2bs4, noexec, noexec, noexec, noexec, noexec, exec
	inorm	xr2ebs4, noexec, noexec, exec, noexec, noexec, exec
	inorm	xr2cbs4, noexec, noexec, noexec, exec, noexec, exec
	inorm	xr2ecbs4, noexec, noexec, exec, exec, noexec, exec
	inorm	xi2bs4, exec, noexec, noexec, noexec, noexec, exec
	inorm	xi2ebs4, exec, noexec, exec, noexec, noexec, exec
	inorm	xi2cbs4, exec, noexec, noexec, exec, noexec, exec
	inorm	xi2ecbs4, exec, noexec, exec, exec, noexec, exec
	zpnorm	xr2zpbs4, noexec, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xr2zpbs4c1, noexec, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xr2zpbs4cm1, noexec, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xr2zpebs4, noexec, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xr2zpebs4c1, noexec, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xr2zpebs4cm1, noexec, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xr2zpcbs4, noexec, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr2zpecbs4, noexec, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xi2zpbs4, exec, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xi2zpbs4c1, exec, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xi2zpbs4cm1, exec, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xi2zpebs4, exec, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	xi2zpebs4c1, exec, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	xi2zpebs4cm1, exec, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	xi2zpcbs4, exec, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xi2zpecbs4, exec, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	xr2zpbs4k, noexec, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr2zpbs4kc1, noexec, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xr2zpbs4kcm1, noexec, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xr2zpebs4k, noexec, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr2zpebs4kc1, noexec, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xr2zpebs4kcm1, noexec, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xr2zpcbs4k, noexec, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xr2zpecbs4k, noexec, exec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi2zpbs4k, exec, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi2zpbs4kc1, exec, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xi2zpbs4kcm1, exec, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xi2zpebs4k, exec, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi2zpebs4kc1, exec, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	xi2zpebs4kcm1, exec, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	xi2zpcbs4k, exec, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	xi2zpecbs4k, exec, exec, exec, noexec, exec, noexec, noexec, noexec

_TEXT	ENDS
END
