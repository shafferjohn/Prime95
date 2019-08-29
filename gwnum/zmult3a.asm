; Copyright 2011-2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routines implement the AVX version of normalization part of a r4dwpn (radix-4 delayed with partial normalization) FFT.
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
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization routine.

;; On input (see pass1_normalize macro in zmult.mac):
;; rsi = FFT data (scratch area for two-pass FFTs, DESTARG for one-pass FFTs)
;; rdx = count of clmblks to normalize
;; rcx = address to do clm prefetching
;; rbx = clm-1 (mask used to detect every clm-th iteration where rcx must bump by a large stride)
;; r12 = large stride increment (which equals -clm*128 + large stride)
;; rdi = pointer to big/lit flags
;; r10 = address of inverse group multipliers

;; During computation registers used thusly:
;; rax = scratch register
;; rsi = 1st ptr to FFT data area
;; rdi = big/lit/fudge flags array
;; rbp = carries array
;; r15 = clm loop counter
;; r8 = compressed biglit table pointer
;; r9 = alternating zero and one for incrementing rdi during a zero pad normalize
;; rdx = outer loop counter (which gets pushed/popped), then register for loading compressed biglit index
;; r13 = 3rd ptr to FFT data area
;; r14 = distance between 1st and 2nd FFT data area (as well as 3rd and 4th)

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, ilp1, ilp2, not8, done
	PROCFLP	lab
	int_prolog 0,0,0

no zero	mov	eax, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	eax, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	eax, ADDIN_OFFSET	;; Get address to add value into
no zero	vmovsd	xmm0, ADDIN_VALUE	;; Get the requested add-in value
no zero	vaddsd	xmm0, xmm0, Q [rsi][rax] ;; Add in the FFT value
no zero	vmovsd	Q [rsi][rax], xmm0	;; Save the new value
noadd:	

echk	vbroadcastsd zmm31, MAXERR	;; Load maximum error
	znorm_wpn_preload ttp, zero, echk, const

	mov	rbp, carries		;; Addr of the carries

;; BUG - why calculate this every time?
	;; Calculate distance between 1st and 2nd FFT data areas (cache_line_multiplier * 128 + normblkdst)
	mov	r14d, cache_line_multiplier
	shl	r14, 7
	add	r14, normblkdst

ttp	mov	r8, compressed_biglits	;; Get pointer to compressed biglit table
ilp1:	push	rdx			;; Save outer loop counter
	mov	r15d, cache_line_multiplier ;; Load inner loop counter (8*clm)
	lea	r13, [rsi+r14*2]	;; Calc 3rd src pointer
;; BUG/OPT?? use rax register instead of rdx to save a push/pop?
ttp	sub	rdx, rdx		;; Clear register used to load compressed biglit index

	vmovapd zmm0, [rbp+0*128]	;; Load low carries
	vmovapd zmm1, [rbp+1*128]
	vmovapd zmm2, [rbp+2*128]
	vmovapd zmm3, [rbp+3*128]
	vmovapd zmm4, [rbp+0*128+64]	;; Load high carries
	vmovapd zmm5, [rbp+1*128+64]
	vmovapd zmm6, [rbp+2*128+64]
	vmovapd zmm7, [rbp+3*128+64]

const ttp vsubpd zmm0, zmm0, zmm30	;; Remove RNDVAL from irrational mul-by-small-const normalization
const ttp vsubpd zmm1, zmm1, zmm30
const ttp vsubpd zmm2, zmm2, zmm30
const ttp vsubpd zmm3, zmm3, zmm30
const ttp vsubpd zmm4, zmm4, zmm30
const ttp vsubpd zmm5, zmm5, zmm30
const ttp vsubpd zmm6, zmm6, zmm30
const ttp vsubpd zmm7, zmm7, zmm30

ilp2:	znorm_wpn ttp, zero, echk, const ;; Normalize 64 values (4 simultaneous clmblks)
	bump	rsi, 128		;; Next 1st source pointer
	bump	r13, 128		;; Next 3rd source pointer
ttp	bump	rdi, 1			;; Next big/little flags
	sub	r15, 1			;; Test clm loop counter

	prefetchwt1 [rcx]		;; Prefetch with write intent FFT data for next pass 1 block
	bump	rcx, 64			;; Tentatively compute next cache line to prefetch
	lea	rax, [rcx+r12]		;; Tentatively compute next cache line to prefetch if it is time for a large stride
	test	r15, rbx		;; Compare loop counter to the clm-1 mask
	cmovz	rcx, rax		;; Every clm-th iteration use a large stride

	test	r15, r15		;; Are we done with the 8*clm inner loop?
	jnz	ilp2			;; Loop til done

const ttp vaddpd zmm0, zmm0, zmm30	;; Reapply RNDVAL for irrational mul-by-small-const normalization
const ttp vaddpd zmm1, zmm1, zmm30
const ttp vaddpd zmm2, zmm2, zmm30
const ttp vaddpd zmm3, zmm3, zmm30
const ttp vaddpd zmm4, zmm4, zmm30
const ttp vaddpd zmm5, zmm5, zmm30
const ttp vaddpd zmm6, zmm6, zmm30
const ttp vaddpd zmm7, zmm7, zmm30

	zstore	[rbp+0*128], zmm0	;; Save low carries
	zstore	[rbp+1*128], zmm1
	zstore	[rbp+2*128], zmm2
	zstore	[rbp+3*128], zmm3
	zstore	[rbp+0*128+64], zmm4	;; Save high carries
	zstore	[rbp+1*128+64], zmm5
	zstore	[rbp+2*128+64], zmm6
	zstore	[rbp+3*128+64], zmm7
	bump	rbp, 4*128		;; Next set of carries

ttp	bump	r10, 4*128		;; Next set of inverse group multipliers
no ttp	bump	r10, 4*64		;; Next set of inverse group multipliers

	add	r13, normblkdst		;; Add 0 or 64 every clmblkdst.  r13 is now the just completed 4th source ptr
	lea	rsi, [r13+r14]		;; Calculate next 1st source ptr 
	add	rsi, normblkdst4	;; Add 0 or 64 every 4KB for a one-pass FFT

;;BUG - change/alias name of normblkdst to clmblkdst_pad?
;;BUG - change/alias name of normblkdst8 to clmblkdst8_pad?

	pop	rdx			;; Restore outer loop counter
	test	rdx, 4			;; Test for a multiple of 8 blocks
	jz	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	rdx, 4			;; Test loop counter of clmblks
	jnz	ilp1			;; Iterate

echk	zcollapse_maxerr		;; Collapse and save MAXERR

	; Handle adjusting the carry out of the topmost FFT word

;; BUG where in C code do we set THIS_BLOCK and LAST_PASS1_BLOCK for one-pass??

ttp	mov	eax, THIS_BLOCK		;; Check for processing last block
ttp	cmp	eax, LAST_PASS1_BLOCK
ttp	jne	done			;; Jump if not last block
;; BUG/OPT - isn't rbp pointing just past last carry?? isn't last carry in zmm7 now? znorm_top_carry_wpn doesn't use this info
ttp	znorm_top_carry_wpn		;; Adjust carry in xmm7 if k > 1

done:	int_epilog 0,0,0
	ENDPP	lab
	ENDM


zpnorm	MACRO	lab, ttp, echk, const, khi
	LOCAL	ilp1, ilp2, not8
	PROCFLP	lab
	int_prolog 0,0,0

echk	vbroadcastsd zmm31, MAXERR	;; Load maximum error
	znorm_wpn_zpad_preload ttp, echk, const, khi

	mov	rbp, carries		;; Addr of the carries

	;; BUG - why calculate this every time?
	;; Calculate distance between 1st and 2nd FFT data areas (cache_line_multiplier * 128 + normblkdst)
	mov	r14d, cache_line_multiplier
	shl	r14, 7
	add	r14, normblkdst

ttp	mov	r8, compressed_biglits	;; Get pointer to compressed biglit table
ttp	sub	r9, r9			;; Clear register used to alternate the increment of rdi
ilp1:	push	rdx			;; Save outer loop counter
	mov	r15d, cache_line_multiplier ;; Load inner loop counter
	lea	r13, [rsi+r14*2]	;; Calc 3rd src pointer
ttp	sub	rdx, rdx		;; Clear register used to load compressed biglit index

	vmovapd zmm0, [rbp+0*128]	;; Load low carries
	vmovapd zmm1, [rbp+1*128]
	vmovapd zmm2, [rbp+2*128]
	vmovapd zmm3, [rbp+3*128]
	vmovapd zmm4, [rbp+0*128+64]	;; Load high carries
	vmovapd zmm5, [rbp+1*128+64]
	vmovapd zmm6, [rbp+2*128+64]
	vmovapd zmm7, [rbp+3*128+64]
ttp	vsubpd	zmm4, zmm4, zmm30	;; Remove RNDVAL from high carries for irrational normalization
ttp	vsubpd	zmm5, zmm5, zmm30
ttp	vsubpd	zmm6, zmm6, zmm30
ttp	vsubpd	zmm7, zmm7, zmm30

ilp2:	znorm_wpn_zpad ttp, echk, const, khi ;; Normalize 64 values (4 simultaneous clmblks)
	bump	rsi, 128		;; Next 1st source pointer
	bump	r13, 128		;; Next 3rd source pointer
	sub	r15, 1			;; Test clm loop counter

ttp	xor	r8, 4			;; Bump or unbump pointer into the compressed biglit table
ttp	add	rdi, r9			;; Increment (or not) the biglit table pointer
ttp	xor	r9, 1			;; Alternate the increment for rdi

	prefetchwt1 [rcx]		;; Prefetch with write intent FFT data for next pass 1 block
	bump	rcx, 64			;; Tentatively compute next cache line to prefetch
	lea	rax, [rcx+r12]		;; Tentatively compute next cache line to prefetch if it is time for a large stride
	test	r15, rbx		;; Compare loop counter to the clm-1 mask
	cmovz	rcx, rax		;; Every clm-th iteration use a large stride

	test	r15, r15		;; Are we done with the 8*clm inner loop?
	jnz	ilp2			;; Loop til done

ttp	vaddpd	zmm4, zmm4, zmm30	;; Restore RNDVAL to high carries for irrational normalization
ttp	vaddpd	zmm5, zmm5, zmm30
ttp	vaddpd	zmm6, zmm6, zmm30
ttp	vaddpd	zmm7, zmm7, zmm30
	zstore	[rbp+0*128], zmm0	;; Save low carries
	zstore	[rbp+1*128], zmm1
	zstore	[rbp+2*128], zmm2
	zstore	[rbp+3*128], zmm3
	zstore	[rbp+0*128+64], zmm4	;; Save high carries
	zstore	[rbp+1*128+64], zmm5
	zstore	[rbp+2*128+64], zmm6
	zstore	[rbp+3*128+64], zmm7
	bump	rbp, 4*128		;; Next set of carries

ttp	bump	r10, 4*64		;; Next set of inverse group multipliers

	add	r13, normblkdst		;; Add 0 or 64 every clmblkdst.  r13 is now the just completed 4th source ptr
	lea	rsi, [r13+r14]		;; Calculate next 1st source ptr 
	add	rsi, normblkdst4	;; Add 0 or 64 every 4KB for a one-pass FFT

	pop	rdx			;; Restore outer loop counter
	test	rdx, 4			;; Test for a multiple of 8 blocks
	jz	short not8
	add	rsi, normblkdst8	;; Add 64 or -64 every 8 clmblkdsts
not8:	sub	rdx, 4			;; Test loop counter of clmblks
	jnz	ilp1			;; Iterate

echk	zcollapse_maxerr		;; Collapse and save MAXERR

	int_epilog 0,0,0
	ENDPP	lab
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	zr3, noexec, noexec, noexec, noexec
	inorm	zr3e, noexec, noexec, exec, noexec
	inorm	zr3c, noexec, noexec, noexec, exec
	inorm	zr3ec, noexec, noexec, exec, exec
	inorm	zr3z, noexec, exec, noexec, noexec
	inorm	zr3ze, noexec, exec, exec, noexec
	inorm	zi3, exec, noexec, noexec, noexec
	inorm	zi3e, exec, noexec, exec, noexec
	inorm	zi3c, exec, noexec, noexec, exec
	inorm	zi3ec, exec, noexec, exec, exec
	inorm	zi3z, exec, exec, noexec, noexec
	inorm	zi3ze, exec, exec, exec, noexec
	zpnorm	zr3zp, noexec, noexec, noexec, exec
	zpnorm	zr3zpe, noexec, exec, noexec, exec
	zpnorm	zr3zpc, noexec, noexec, exec, exec
	zpnorm	zr3zpec, noexec, exec, exec, exec
	zpnorm	zi3zp, exec, noexec, noexec, exec
	zpnorm	zi3zpe, exec, exec, noexec, exec
	zpnorm	zi3zpc, exec, noexec, exec, exec
	zpnorm	zi3zpec, exec, exec, exec, exec
	zpnorm	zr3zpk, noexec, noexec, noexec, noexec
	zpnorm	zr3zpek, noexec, exec, noexec, noexec
	zpnorm	zr3zpck, noexec, noexec, exec, noexec
	zpnorm	zr3zpeck, noexec, exec, exec, noexec
	zpnorm	zi3zpk, exec, noexec, noexec, noexec
	zpnorm	zi3zpek, exec, exec, noexec, noexec
	zpnorm	zi3zpck, exec, noexec, exec, noexec
	zpnorm	zi3zpeck, exec, exec, exec, noexec

_TEXT	ENDS
END
