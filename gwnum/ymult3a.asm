; Copyright 2011-2017 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routines implement the AVX version of normalization part of a r4dwpn (radix-4 delayed with partial normalization) FFT.
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
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

IFNDEF X86_64

saved_rsi	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]
blk8_counter	EQU	BYTE PTR [rsp+first_local+SZPTR+12]

inorm	MACRO	lab, ttp, zero, echk, const, base2
	LOCAL	noadd, ilp0, ilp1, ilp2, not8, done
	PROCFLP	lab
	int_prolog SZPTR+16,0,0
echk	vmovapd	ymm6, YMM_MAXERR	;; Load maximum error
no zero	mov	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	edx, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	vmovsd	xmm0, ADDIN_VALUE	;; Get the requested add-in value
no zero	vaddsd	xmm0, xmm0, Q [rsi][rdi] ;; Add in the FFT value
no zero	vmovsd	Q [rsi][rdi], xmm0	;; Save the new value
noadd:	mov	saved_rsi, rsi		;; Save for top_carry_adjust

	ynorm_wpn_preload ttp, base2, zero, echk, const

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little & fudge flags array ptr
	mov	blk8_counter, 0		;; Clear counter
	mov	eax, count3		;; Load count of grp multipliers
	mov	loopcount3, eax		;; Save loop counter
ttp	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
ilp0:	mov	eax, count2		;; Load wpn count
	mov	loopcount2, eax		;; Save count
ilp1:	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, rax		;; Save loop counter
	vmovapd ymm2, [rbp+0*32]	;; Load carries
	vmovapd ymm3, [rbp+1*32]
ilp2:	ynorm_wpn ttp, base2, zero, echk, const ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rdi, 2			;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp2			;; Loop til done
	ystore	[rbp+0*32], ymm2	;; Save carries
	ystore	[rbp+1*32], ymm3
	add	rsi, normblkdst		;; Add 0 or 64 every clmblkdst
	bump	rbp, 64			;; Next set of carries
	add	blk8_counter, 80h/4	;; Test for a multiple of 8 blocks
	jnc	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	loopcount2, 1		;; Test loop counter
	jnz	ilp1			;; Iterate
ttp	bump	rdx, 2*YMM_GMD		;; Next set of group multipliers
	sub	loopcount3, 1		;; Test outer loop counter
	jnz	ilp0			;; Iterate

echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error

	; Handle adjusting the carry out of the topmost FFT word

	mov	eax, THIS_BLOCK		;; Check for processing last block
	cmp	eax, LAST_PASS1_BLOCK
;; BUG - should we jump to 4 common top carry propagate end codes?  does it save much?  there are a lot of inorm variants!
	jne	done			;; Jump if not last block
	mov	rsi, saved_rsi		;; Restore FFT data ptr
;; BUG - isn't rbp pointing just past last carry?? isn't last carry in ymm5 now? ynorm_top_carry_wpn doesn't use this info
ttp	ynorm_top_carry_wpn ttp, base2	;; Adjust carry if k > 1

done:	int_epilog SZPTR+16,0,0
	ENDPP	lab
	ENDM

ELSE

IFDEF YIMPL_WPN1_FFTS

saved_rsi	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]
blk8_counter	EQU	BYTE PTR [rsp+first_local+SZPTR+12]

inorm	MACRO	lab, ttp, zero, echk, const, base2
	LOCAL	noadd, noinc, ilp0, ilp1, ilp2, not8, done
	PROCFLP	lab
	int_prolog SZPTR+16,0,0
echk	vmovapd	ymm6, YMM_MAXERR	;; Load maximum error
no zero	mov	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	edx, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	vmovsd	xmm0, ADDIN_VALUE	;; Get the requested add-in value
no zero	vaddsd	xmm0, xmm0, Q [rsi][rdi] ;; Add in the FFT value
no zero	vmovsd	Q [rsi][rdi], xmm0	;; Save the new value
noadd:	mov	saved_rsi, rsi		;; Save for top_carry_adjust

	ynorm_wpn_preload ttp, base2, zero, echk, const

	mov	r12, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
ttp	mov	rdi, norm_ptr1		;; Load big/little & fudge flags array ptr
	mov	blk8_counter, 0		;; Clear counter
	mov	eax, count5		;; Load count of grp multipliers (divided by 2)
	mov	loopcount3, eax		;; Save loop counter
ilp0:	mov	eax, count4		;; Load wpn count (divided by 2)
	mov	loopcount2, eax		;; Save count
	mov	r15, r12		;; Calc 2nd group multiplier pointer
ttp	cmp	count2, 1		;; If count is one, 2nd group multiplier must point to the next group
ttp	jne	short noinc		;; If count2 is more than one, 2nd group multiplier is same as 1st group multiplier
ttp	bump	r15, 2*YMM_GMD		;; Bump 2nd group multiplier pointer
noinc:
ttp	lea	r9, [r15+2*YMM_GMD]	;; Prefetch pointer for group multipliers
ilp1:	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, eax		;; Save loop counter
ttp	lea	r14, [rdi+2*rax]	;; Calc 2nd big/lit ptr
	shl	rax, 6			;; Calc 2nd src pointer (rsi + clm * 64 + normblkdst)
	lea	r13, [rsi+rax]
	add	r13, normblkdst
;; OPTIMIZATION - old wpn didn't reload rbx every clm block.  Rearranging big/lit data could achieve this as well (and can save the r14 register in ynorm_wpn)
ttp	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
ttp	movzx	rcx, WORD PTR [r14]	;; Preload 4 big vs. little & fudge flags
	vmovapd ymm2, [rbp+0*32]	;; Load carries
	vmovapd ymm3, [rbp+1*32]
	vmovapd ymm9, [rbp+2*32]
	vmovapd ymm10, [rbp+3*32]
ilp2:	ynorm_wpn ttp, base2, zero, echk, const ;; Normalize 2 sets of 8 values
	bump	rsi, 64			;; Next cache line
	bump	r13, 64			;; Next cache line
ttp	bump	rdi, 2			;; Next big/little flags
ttp	bump	r14, 2			;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp2			;; Loop til done
	ystore	[rbp+0*32], ymm2	;; Save carries
	ystore	[rbp+1*32], ymm3
	ystore	[rbp+2*32], ymm9
	ystore	[rbp+3*32], ymm10
	bump	rbp, 4*32		;; Next set of carries
ttp	mov	rdi, r14		;; Calculate address of next big/lit ptr
	mov	rsi, r13		;; Calculate address of next source
	add	rsi, normblkdst		;; Add 0 or 64 every clmblkdst
	add	blk8_counter, 80h/2	;; Test for a multiple of 8 blocks
	jnc	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	loopcount2, 1		;; Test loop counter
	jnz	ilp1			;; Iterate
ttp	lea	r12, [r15+2*YMM_GMD]	;; Next set of group multipliers
	sub	loopcount3, 1		;; Test outer loop counter
	jnz	ilp0			;; Iterate

echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error

	; Handle adjusting the carry out of the topmost FFT word

	mov	eax, THIS_BLOCK		;; Check for processing last block
	cmp	eax, LAST_PASS1_BLOCK
;; BUG - should we jump to 4 common top carry propagate end codes?  does it save much?  there are a lot of inorm variants!
	jne	done			;; Jump if not last block
	mov	rsi, saved_rsi		;; Restore FFT data ptr
;; BUG - isn't rbp pointing just past last carry?? isn't last carry in ymm5 now? ynorm_top_carry_wpn doesn't use this info
	ynorm_top_carry_wpn ttp, base2	;; Adjust carry if k > 1

done:	int_epilog SZPTR+16,0,0
	ENDPP	lab
	ENDM
ENDIF

;; In wpn4 FFTs we know count2 will be even.  Consequently, we only need one group pointer in ynorm_wpn

IFDEF YIMPL_WPN4_FFTS

saved_rsi	EQU	PPTR [rsp+first_local]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+0]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+4]
blk8_counter	EQU	BYTE PTR [rsp+first_local+SZPTR+8]

inorm	MACRO	lab, ttp, zero, echk, const, base2
	LOCAL	noadd, ilp0, ilp1, ilp2, not8, done
	PROCFLP	lab
	int_prolog SZPTR+12,0,0
echk	vmovapd	ymm6, YMM_MAXERR	;; Load maximum error
no zero	mov	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	edx, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	vmovsd	xmm0, ADDIN_VALUE	;; Get the requested add-in value
no zero	vaddsd	xmm0, xmm0, Q [rsi][rdi] ;; Add in the FFT value
no zero	vmovsd	Q [rsi][rdi], xmm0	;; Save the new value
noadd:	mov	saved_rsi, rsi		;; Save for top_carry_adjust

	ynorm_wpn_preload ttp, base2, zero, echk, const

	mov	r12, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
ttp	mov	rdi, norm_ptr1		;; Load big/little & fudge flags array ptr
	mov	blk8_counter, 0		;; Clear counter
	mov	eax, count5		;; Load count of grp multipliers (divided by 2)
	mov	loopcount3, eax		;; Save loop counter
ilp0:	mov	eax, count4		;; Load wpn count (divided by 2)
	mov	loopcount2, eax		;; Save count
ttp	lea	r9, [r12+2*YMM_GMD]	;; Prefetch pointer for group multipliers
ilp1:	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	r15, rax		;; Save loop counter
ttp	lea	r14, [rdi+2*rax]	;; Calc 2nd big/lit ptr
	shl	rax, 6			;; Calc 2nd src pointer (rsi + clm * 64 + normblkdst)
	lea	r13, [rsi+rax]
	add	r13, normblkdst
;; OPTIMIZATION - old wpn didn't reload rbx every clm block.  Rearranging big/lit data could achieve this as well (and can save the r14 register in ynorm_wpn)
ttp	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
ttp	movzx	rcx, WORD PTR [r14]	;; Preload 4 big vs. little & fudge flags
	vmovapd ymm2, [rbp+0*32]	;; Load carries
	vmovapd ymm3, [rbp+1*32]
	vmovapd ymm9, [rbp+2*32]
	vmovapd ymm10, [rbp+3*32]
ilp2:	ynorm_wpn ttp, base2, zero, echk, const ;; Normalize 2 sets of 8 values
	bump	rsi, 64			;; Next cache line
	bump	r13, 64			;; Next cache line
ttp	bump	rdi, 2			;; Next big/little flags
ttp	bump	r14, 2			;; Next big/little flags
	sub	r15, 1			;; Test loop counter
	jnz	ilp2			;; Loop til done
	ystore	[rbp+0*32], ymm2	;; Save carries
	ystore	[rbp+1*32], ymm3
	ystore	[rbp+2*32], ymm9
	ystore	[rbp+3*32], ymm10
	bump	rbp, 4*32		;; Next set of carries
ttp	mov	rdi, r14		;; Calculate address of next big/lit ptr
	mov	rsi, r13		;; Calculate address of next source
	add	rsi, normblkdst		;; Add 0 or 64 every clmblkdst
	add	blk8_counter, 80h/2	;; Test for a multiple of 8 blocks
	jnc	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	loopcount2, 1		;; Test loop counter
	jnz	ilp1			;; Iterate
ttp	lea	r12, [r12+2*YMM_GMD]	;; Next set of group multipliers
	sub	loopcount3, 1		;; Test outer loop counter
	jnz	ilp0			;; Iterate

echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error

	; Handle adjusting the carry out of the topmost FFT word

	mov	eax, THIS_BLOCK		;; Check for processing last block
	cmp	eax, LAST_PASS1_BLOCK
;; BUG - should we jump to 4 common top carry propagate end codes?  does it save much?  there are a lot of inorm variants!
	jne	done			;; Jump if not last block
	mov	rsi, saved_rsi		;; Restore FFT data ptr
;; BUG - isn't rbp pointing just past last carry?? isn't last carry in ymm5 now? ynorm_top_carry_wpn doesn't use this info
	ynorm_top_carry_wpn ttp, base2	;; Adjust carry if k > 1

done:	int_epilog SZPTR+12,0,0
	ENDPP	lab
	ENDM
ENDIF

ENDIF


IFNDEF X86_64

loopcount2z	EQU	DPTR [rsp+first_local]
loopcount3z	EQU	DPTR [rsp+first_local+4]
blk8_counterz	EQU	BYTE PTR [rsp+first_local+8]

zpnorm	MACRO	lab, ttp, echk, const, base2, khi, c1, cm1
	LOCAL	ilp0, ilp1, ilp2, not8
	PROCFLP	lab
	int_prolog 12,0,0
echk	vmovapd	ymm6, YMM_MAXERR	;; Load maximum error

	ynorm_wpn_zpad_preload ttp, base2, echk, const, khi, c1, cm1

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	blk8_counterz, 0	;; Clear counter
	mov	eax, count3		;; Load count of grp multipliers
	mov	loopcount3z, eax	;; Save loop counter
ttp	movzx	rbx, WORD PTR [rdi]	;; Preload big vs. little & fudge flags
ilp0:	mov	eax, count2		;; Load wpn count
	mov	loopcount2z, eax	;; Save loop counter
ilp1:	mov	ecx, cache_line_multiplier ;; Load inner loop counter
	ystore	ymm2, [rbp+0*32]	;; Preload carries
	ystore	ymm3, [rbp+1*32]
ilp2:	ynorm_wpn_zpad ttp, base2, echk, const, khi, c1, cm1 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rdi, 2			;; Next big/little flags
	sub	rcx, 1			;; Test loop counter
	jnz	ilp2			;; Loop til done
	ystore	[rbp+0*32], ymm2	;; Store carries
	ystore	[rbp+1*32], ymm3
	add	rsi, normblkdst		;; Add 0 or 64 every clmblkdst
	bump	rbp, 64			;; Next set of carries
	add	blk8_counterz, 80h/4	;; Test for a multiple of 8 blocks
	jnc	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	loopcount2z, 1		;; Test loop counter
	jnz	ilp1			;; Iterate
ttp	bump	rdx, YMM_GMD		;; Next set of group multipliers
	sub	loopcount3z, 1		;; Test outer loop counter
	jnz	ilp0			;; Iterate

echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error
	int_epilog 12,0,0
	ENDPP	lab
	ENDM

ELSE

IFDEF YIMPL_WPN1_FFTS

loopcount2z	EQU	DPTR [rsp+first_local+0]

zpnorm	MACRO	lab, ttp, echk, const, base2, khi, c1, cm1
	LOCAL	noinc, ilp0, ilp1, ilp2, not8
	PROCFLP	lab
	int_prolog 4,0,0
echk	vmovapd	ymm6, YMM_MAXERR	;; Load maximum error

	ynorm_wpn_zpad_preload ttp, base2, echk, const, khi, c1, cm1

	mov	r12, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, count5		;; Load count of grp multipliers (divided by 2)
	mov	r10, rax		;; Save loop counter
ilp0:	mov	eax, count4		;; Load wpn count (divided by 2)
	mov	loopcount2z, eax	;; Save loop counter
	mov	r15, r12		;; Calc 2nd group multiplier pointer
ttp	cmp	count2, 1		;; If count is one, 2nd group multiplier must point to the next group
ttp	jne	short noinc		;; If count2 is more than one, 2nd group multiplier is same as 1st group multiplier
ttp	bump	r15, YMM_GMD		;; Bump 2nd group multiplier pointer
noinc:
ttp	lea	r9, [r15+YMM_GMD]	;; Prefetch pointer for group multipliers
ilp1:	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	r8, rax			;; Save loop counter
ttp	lea	r14, [rdi+2*rax]	;; Calc 2nd big/lit ptr
	shl	rax, 6			;; Calc 2nd src pointer (rsi + clm * 64 + normblkdst)
	lea	r13, [rsi+rax]
	add	r13, normblkdst
;; OPTIMIZATION - old wpn didn't reload rbx every clm block.  Rearranging big/lit data could achieve this as well (and can save the r14 register in ynorm_wpn)
ttp	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
ttp	movzx	rcx, WORD PTR [r14]	;; Preload 4 big vs. little & fudge flags
	vmovapd	ymm2, [rbp+0*32]	;; Preload carries
	vmovapd	ymm3, [rbp+1*32]
	vmovapd ymm9, [rbp+2*32]
	vmovapd ymm10, [rbp+3*32]
;; OPTIMZATION - always store carries 3/10 in the +BIGVAL format
IF (@INSTR(,%yarch,<FMA3>) NE 0)
	vaddpd	ymm3, ymm3, ymm12	;; Add in YMM_BIGVAL
	vaddpd	ymm10, ymm10, ymm12
ENDIF
ilp2:	ynorm_wpn_zpad ttp, base2, echk, const, khi, c1, cm1 ;; Normalize 2 sets of 8 values
	bump	rsi, 64			;; Next cache line
	bump	r13, 64			;; Next cache line
ttp	bump	rdi, 2			;; Next big/little flags
ttp	bump	r14, 2			;; Next big/little flags
	sub	r8, 1			;; Test loop counter
	jnz	ilp2			;; Loop til done
;; OPTIMZATION - always store carries 3/10 in the +BIGVAL format
IF (@INSTR(,%yarch,<FMA3>) NE 0)
	vsubpd	ymm3, ymm3, ymm12	;; Subtract out YMM_BIGVAL
	vsubpd	ymm10, ymm10, ymm12
ENDIF
	ystore	[rbp+0*32], ymm2	;; Store carries
	ystore	[rbp+1*32], ymm3
	ystore	[rbp+2*32], ymm9
	ystore	[rbp+3*32], ymm10
	bump	rbp, 4*32		;; Next set of carries
ttp	mov	rdi, r14		;; Calculate address of next big/lit ptr
	mov	rsi, r13		;; Calculate address of next source
	add	rsi, normblkdst		;; Add 0 or 64 every clmblkdst
	add	r10w, 8000h/2		;; Test for a multiple of 8 blocks
	jnc	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	loopcount2z, 1		;; Test loop counter
	jnz	ilp1			;; Iterate
ttp	lea	r12, [r15+YMM_GMD]	;; Next set of group multipliers
	sub	r10w, 1			;; Test outer loop counter
	jnz	ilp0			;; Iterate

echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error
	int_epilog 4,0,0
	ENDPP	lab
	ENDM

ENDIF

;; In wpn4 FFTs we know count2 will be even.  Consequently, we only need one group pointer in ynorm_wpn

IFDEF YIMPL_WPN4_FFTS

zpnorm	MACRO	lab, ttp, echk, const, base2, khi, c1, cm1
	LOCAL	ilp0, ilp1, ilp2, not8
	PROCFLP	lab
	int_prolog 0,0,0
echk	vmovapd	ymm6, YMM_MAXERR	;; Load maximum error

	ynorm_wpn_zpad_preload ttp, base2, echk, const, khi, c1, cm1

	mov	r12, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, count5		;; Load count of grp multipliers (divided by 2)
	mov	r10, rax		;; Save loop counter
ilp0:	mov	eax, count4		;; Load wpn count (divided by 2)
	mov	r15, rax		;; Save loop counter
ttp	lea	r9, [r12+YMM_GMD]	;; Prefetch pointer for group multipliers
ilp1:	mov	eax, cache_line_multiplier ;; Load inner loop counter
	mov	r8, rax			;; Save loop counter
ttp	lea	r14, [rdi+2*rax]	;; Calc 2nd big/lit ptr
	shl	rax, 6			;; Calc 2nd src pointer (rsi + clm * 64 + normblkdst)
	lea	r13, [rsi+rax]
	add	r13, normblkdst
;; OPTIMIZATION - old wpn didn't reload rbx every clm block.  Rearranging big/lit data could achieve this as well (and can save the r14 register in ynorm_wpn)
ttp	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
ttp	movzx	rcx, WORD PTR [r14]	;; Preload 4 big vs. little & fudge flags
	vmovapd	ymm2, [rbp+0*32]	;; Preload carries
	vmovapd	ymm3, [rbp+1*32]
	vmovapd ymm9, [rbp+2*32]
	vmovapd ymm10, [rbp+3*32]
;; OPTIMZATION - always store carries 3/10 in the +BIGVAL format
IF (@INSTR(,%yarch,<FMA3>) NE 0)
	vaddpd	ymm3, ymm3, ymm12	;; Add in YMM_BIGVAL
	vaddpd	ymm10, ymm10, ymm12
ENDIF
ilp2:	ynorm_wpn_zpad ttp, base2, echk, const, khi, c1, cm1 ;; Normalize 2 sets of 8 values
	bump	rsi, 64			;; Next cache line
	bump	r13, 64			;; Next cache line
ttp	bump	rdi, 2			;; Next big/little flags
ttp	bump	r14, 2			;; Next big/little flags
	sub	r8, 1			;; Test loop counter
	jnz	ilp2			;; Loop til done
;; OPTIMZATION - always store carries 3/10 in the +BIGVAL format
IF (@INSTR(,%yarch,<FMA3>) NE 0)
	vsubpd	ymm3, ymm3, ymm12	;; Subtract out YMM_BIGVAL
	vsubpd	ymm10, ymm10, ymm12
ENDIF
	ystore	[rbp+0*32], ymm2	;; Store carries
	ystore	[rbp+1*32], ymm3
	ystore	[rbp+2*32], ymm9
	ystore	[rbp+3*32], ymm10
	bump	rbp, 4*32		;; Next set of carries
ttp	mov	rdi, r14		;; Calculate address of next big/lit ptr
	mov	rsi, r13		;; Calculate address of next source
	add	rsi, normblkdst		;; Add 0 or 64 every clmblkdst
	add	r10w, 8000h/2		;; Test for a multiple of 8 blocks
	jnc	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	r15, 1			;; Test loop counter
	jnz	ilp1			;; Iterate
ttp	lea	r12, [r12+YMM_GMD]	;; Next set of group multipliers
	sub	r10w, 1			;; Test outer loop counter
	jnz	ilp0			;; Iterate

echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error
	int_epilog 0,0,0
	ENDPP	lab
	ENDM

ENDIF

ENDIF

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	yr3, noexec, noexec, noexec, noexec, exec
	inorm	yr3e, noexec, noexec, exec, noexec, exec
	inorm	yr3c, noexec, noexec, noexec, exec, exec
	inorm	yr3ec, noexec, noexec, exec, exec, exec
	inorm	yr3z, noexec, exec, noexec, noexec, exec
	inorm	yr3ze, noexec, exec, exec, noexec, exec
	inorm	yi3, exec, noexec, noexec, noexec, exec
	inorm	yi3e, exec, noexec, exec, noexec, exec
	inorm	yi3c, exec, noexec, noexec, exec, exec
	inorm	yi3ec, exec, noexec, exec, exec, exec
	inorm	yi3z, exec, exec, noexec, noexec, exec
	inorm	yi3ze, exec, exec, exec, noexec, exec
	zpnorm	yr3zp, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	yr3zpc1, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	yr3zpcm1, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	yr3zpe, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	yr3zpec1, noexec, exec, noexec, exec, exec, exec, noexec
	zpnorm	yr3zpecm1, noexec, exec, noexec, exec, exec, noexec, exec
	zpnorm	yr3zpc, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	yr3zpec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	yi3zp, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	yi3zpc1, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	yi3zpcm1, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	yi3zpe, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	yi3zpec1, exec, exec, noexec, exec, exec, exec, noexec
	zpnorm	yi3zpecm1, exec, exec, noexec, exec, exec, noexec, exec
	zpnorm	yi3zpc, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	yi3zpec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	yr3zpk, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	yr3zpkc1, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	yr3zpkcm1, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	yr3zpek, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	yr3zpekc1, noexec, exec, noexec, exec, noexec, exec, noexec
	zpnorm	yr3zpekcm1, noexec, exec, noexec, exec, noexec, noexec, exec
	zpnorm	yr3zpck, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	yr3zpeck, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	yi3zpk, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	yi3zpkc1, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	yi3zpkcm1, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	yi3zpek, exec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	yi3zpekc1, exec, exec, noexec, exec, noexec, exec, noexec
	zpnorm	yi3zpekcm1, exec, exec, noexec, exec, noexec, noexec, exec
	zpnorm	yi3zpck, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	yi3zpeck, exec, exec, exec, exec, noexec, noexec, noexec

	inorm	yr3b, noexec, noexec, noexec, noexec, noexec
	inorm	yr3eb, noexec, noexec, exec, noexec, noexec
	inorm	yr3cb, noexec, noexec, noexec, exec, noexec
	inorm	yr3ecb, noexec, noexec, exec, exec, noexec
	inorm	yi3b, exec, noexec, noexec, noexec, noexec
	inorm	yi3eb, exec, noexec, exec, noexec, noexec
	inorm	yi3cb, exec, noexec, noexec, exec, noexec
	inorm	yi3ecb, exec, noexec, exec, exec, noexec
	zpnorm	yr3zpb, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	yr3zpbc1, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	yr3zpbcm1, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	yr3zpeb, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	yr3zpebc1, noexec, exec, noexec, noexec, exec, exec, noexec
	zpnorm	yr3zpebcm1, noexec, exec, noexec, noexec, exec, noexec, exec
	zpnorm	yr3zpcb, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	yr3zpecb, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	yi3zpb, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	yi3zpbc1, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	yi3zpbcm1, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	yi3zpeb, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	yi3zpebc1, exec, exec, noexec, noexec, exec, exec, noexec
	zpnorm	yi3zpebcm1, exec, exec, noexec, noexec, exec, noexec, exec
	zpnorm	yi3zpcb, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	yi3zpecb, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	yr3zpbk, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yr3zpbkc1, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	yr3zpbkcm1, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	yr3zpebk, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yr3zpebkc1, noexec, exec, noexec, noexec, noexec, exec, noexec
	zpnorm	yr3zpebkcm1, noexec, exec, noexec, noexec, noexec, noexec, exec
	zpnorm	yr3zpcbk, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	yr3zpecbk, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	yi3zpbk, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yi3zpbkc1, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	yi3zpbkcm1, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	yi3zpebk, exec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yi3zpebkc1, exec, exec, noexec, noexec, noexec, exec, noexec
	zpnorm	yi3zpebkcm1, exec, exec, noexec, noexec, noexec, noexec, exec
	zpnorm	yi3zpcbk, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	yi3zpecbk, exec, exec, exec, noexec, noexec, noexec, noexec

_TEXT	ENDS
END
