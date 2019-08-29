; Copyright 2011-2017 Mersenne Research, Inc.  All rights reserved
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
INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE memory.mac
INCLUDE ynormal.mac

_TEXT SEGMENT

;; Only assemble the add/sub/etc routines when compiled with the CORE architecture
;; Maybe someday we can look into whether or not FMA versions would be faster

IF (@INSTR(,%yarch,<CORE>) NE 0)

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwyaddq1
	ad_prolog 0,0,rbx,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
uaddlp2:mov	ebx, count2		; Load loop counter
uaddlp:	vmovapd	ymm0, [rdx]		; Load second number
	vaddpd	ymm0, ymm0, [rcx]	; Add in first number
	vmovapd	ymm1, [rdx+32]		; Load second number
	vaddpd	ymm1, ymm1, [rcx+32]	; Add in first number
	vmovapd	ymm2, [rdx+64]		; Load second number
	vaddpd	ymm2, ymm2, [rcx+64]	; Add in first number
	vmovapd	ymm3, [rdx+96]		; Load second number
	vaddpd	ymm3, ymm3, [rcx+96]	; Add in first number
	ystore	[rsi], ymm0		; Save result
	ystore	[rsi+32], ymm1		; Save result
	ystore	[rsi+64], ymm2		; Save result
	ystore	[rsi+96], ymm3		; Save result
	bump	rcx, 128		; Next source
	bump	rdx, 128		; Next source
	bump	rsi, 128		; Next dest
	sub	ebx, 1			; Check loop counter
	jnz	short uaddlp		; Loop if necessary
	bump	rcx, 64			; Pad 64 bytes
	bump	rdx, 64			; Pad 64 bytes
	bump	rsi, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	short uaddlp2		; Loop if necessary
	ad_epilog 0,0,rbx,rsi
gwyaddq1 ENDP

;;
;; Add two numbers with carry propagation (eight different versions)
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_ttp_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+2*SZPTR+8]

	; Base-2, irrational, not zero-padded
PROCFL	gwyadd1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
b2addsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_mid_cleanup
b2addlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2addlp:
	ynorm_op_1d vaddpd, exec, exec	; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2addlp			; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2addsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2addlp2		; Loop til done
b2addsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_mid_cleanup exec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2addsec

	ynorm_op_1d_cleanup exec, exec, DESTARG	; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyadd1 ENDP

	; Base-2, rational, not zero-padded
PROCFL	gwyaddr1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
b2raddsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
b2raddlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2raddlp:
	ynorm_op_1d vaddpd, noexec, exec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short b2raddlp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2raddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2raddlp2		; Loop til done
b2raddsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_mid_cleanup noexec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2raddsec

	ynorm_op_1d_cleanup noexec, exec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddr1 ENDP

	; Base-2, irrational, zero-padded
PROCFL	gwyaddzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
b2zpaddsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_zpad_mid_cleanup
b2zpaddlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2zpaddlp:
	ynorm_op_1d_zpad vaddpd, exec, exec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2zpaddlp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2zpaddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2zpaddlp2		; Loop til done
b2zpaddsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_zpad_mid_cleanup exec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2zpaddsec

	ynorm_op_1d_zpad_cleanup exec, exec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddzp1 ENDP

	; Base-2, rational, zero-padded
PROCFL	gwyaddrzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
b2rzpaddsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
b2rzpaddlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2rzpaddlp:
	ynorm_op_1d_zpad vaddpd, noexec, exec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short b2rzpaddlp	; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2rzpaddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2rzpaddlp2		; Loop til done
b2rzpaddsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_zpad_mid_cleanup noexec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2rzpaddsec

	ynorm_op_1d_zpad_cleanup noexec, exec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddrzp1 ENDP

	; Not base-2, irrational, not zero-padded
PROCFL	gwyaddn1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nb2addsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_mid_cleanup
nb2addlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2addlp:
	ynorm_op_1d vaddpd, exec, noexec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2addlp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2addsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2addlp2		; Loop til done
nb2addsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_mid_cleanup exec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2addsec

	ynorm_op_1d_cleanup exec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddn1 ENDP

	; Not base-2, rational, not zero-padded
PROCFL	gwyaddnr1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
nb2raddsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
nb2raddlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2raddlp:
	ynorm_op_1d vaddpd, noexec, noexec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short nb2raddlp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2raddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2raddlp2		; Loop til done
nb2raddsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_mid_cleanup noexec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2raddsec

	ynorm_op_1d_cleanup noexec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddnr1 ENDP

	; Not base-2, irrational, zero-padded
PROCFL	gwyaddnzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nb2zpaddsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_zpad_mid_cleanup
nb2zpaddlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2zpaddlp:
	ynorm_op_1d_zpad vaddpd, exec, noexec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2zpaddlp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2zpaddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2zpaddlp2		; Loop til done
nb2zpaddsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_zpad_mid_cleanup exec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2zpaddsec

	ynorm_op_1d_zpad_cleanup exec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddnzp1 ENDP

	; Not base-2, rational, zero-padded
PROCFL	gwyaddnrzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
nb2rzpaddsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
nb2rzpaddlp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2rzpaddlp:
	ynorm_op_1d_zpad vaddpd, noexec, noexec ; Add and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short nb2rzpaddlp	; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2rzpaddsecdn; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2rzpaddlp2		; Loop til done
nb2rzpaddsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_zpad_mid_cleanup noexec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2rzpaddsec

	ynorm_op_1d_zpad_cleanup noexec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddnrzp1 ENDP

;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwysubq1
	ad_prolog 0,0,rbx,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
usublp2:mov	ebx, count2		; Load loop counter
usublp:	vmovapd	ymm0, [rdx]		; Load second number
	vsubpd	ymm0, ymm0, [rcx]	; Subtract first number
	vmovapd	ymm1, [rdx+32]		; Load second number
	vsubpd	ymm1, ymm1, [rcx+32]	; Subtract first number
	vmovapd	ymm2, [rdx+64]		; Load second number
	vsubpd	ymm2, ymm2, [rcx+64]	; Subtract first number
	vmovapd	ymm3, [rdx+96]		; Load second number
	vsubpd	ymm3, ymm3, [rcx+96]	; Subtract first number
	ystore	[rsi], ymm0		; Save result
	ystore	[rsi+32], ymm1		; Save result
	ystore	[rsi+64], ymm2		; Save result
	ystore	[rsi+96], ymm3		; Save result
	bump	rcx, 128		; Next source
	bump	rdx, 128		; Next source
	bump	rsi, 128		; Next dest
	sub	ebx, 1			; Check loop counter
	jnz	short usublp		; Loop if necessary
	bump	rcx, 64			; Pad 64 bytes
	bump	rdx, 64			; Pad 64 bytes
	bump	rsi, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	short usublp2		; Loop if necessary
	ad_epilog 0,0,rbx,rsi
gwysubq1 ENDP

;;
;; Subtract two numbers with carry propagation (eight different versions)
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_ttp_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+2*SZPTR+8]

	; Base-2, irrational, not zero-padded
PROCFL	gwysub1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
b2subsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_mid_cleanup
b2sublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2sublp:
	ynorm_op_1d vsubpd, exec, exec	; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2sublp			; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2subsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2sublp2		; Loop til done
b2subsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_mid_cleanup exec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2subsec

	ynorm_op_1d_cleanup exec, exec, DESTARG	; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysub1 ENDP

	; Base-2, rational, not zero-padded
PROCFL	gwysubr1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
b2rsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
b2rsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2rsublp:
	ynorm_op_1d vsubpd, noexec, exec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short b2rsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2rsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2rsublp2		; Loop til done
b2rsubsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_mid_cleanup noexec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2rsubsec

	ynorm_op_1d_cleanup noexec, exec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubr1 ENDP

	; Base-2, irrational, zero-padded
PROCFL	gwysubzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
b2zpsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_zpad_mid_cleanup
b2zpsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2zpsublp:
	ynorm_op_1d_zpad vsubpd, exec, exec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2zpsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2zpsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2zpsublp2		; Loop til done
b2zpsubsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_zpad_mid_cleanup exec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2zpsubsec

	ynorm_op_1d_zpad_cleanup exec, exec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubzp1 ENDP

	; Base-2, rational, zero-padded
PROCFL	gwysubrzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
b2rzpsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
b2rzpsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2rzpsublp:
	ynorm_op_1d_zpad vsubpd, noexec, exec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short b2rzpsublp	; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2rzpsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	b2rzpsublp2		; Loop til done
b2rzpsubsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_zpad_mid_cleanup noexec, exec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2rzpsubsec

	ynorm_op_1d_zpad_cleanup noexec, exec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubrzp1 ENDP

	; Not base-2, irrational, not zero-padded
PROCFL	gwysubn1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nb2subsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_mid_cleanup
nb2sublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2sublp:
	ynorm_op_1d vsubpd, exec, noexec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2sublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2subsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2sublp2		; Loop til done
nb2subsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_mid_cleanup exec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2subsec

	ynorm_op_1d_cleanup exec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubn1 ENDP

	; Not base-2, rational, not zero-padded
PROCFL	gwysubnr1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
nb2rsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_mid_cleanup
nb2rsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2rsublp:
	ynorm_op_1d vsubpd, noexec, noexec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short nb2rsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2rsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2rsublp2		; Loop til done
nb2rsubsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_mid_cleanup noexec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2rsubsec

	ynorm_op_1d_cleanup noexec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubnr1 ENDP

	; Not base-2, irrational, zero-padded
PROCFL	gwysubnzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nb2zpsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for ynorm_op_1d_zpad_mid_cleanup
nb2zpsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2zpsublp:
	ynorm_op_1d_zpad vsubpd, exec, noexec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2zpsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2zpsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2zpsublp2		; Loop til done
nb2zpsubsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	ynorm_op_1d_zpad_mid_cleanup exec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2zpsubsec

	ynorm_op_1d_zpad_cleanup exec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubnzp1 ENDP

	; Not base-2, rational, zero-padded
PROCFL	gwysubnrzp1
	ad_prolog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
nb2rzpsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for ynorm_op_1d_zpad_mid_cleanup
nb2rzpsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2rzpsublp:
	ynorm_op_1d_zpad vsubpd, noexec, noexec ; Subtract and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	short nb2rzpsublp	; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2rzpsubsecdn; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	nb2rzpsublp2		; Loop til done
nb2rzpsubsecdn:

	mov	rax, saved_dest_ptr	; Restore dest pointer
	ynorm_op_1d_zpad_mid_cleanup noexec, noexec ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2rzpsubsec

	ynorm_op_1d_zpad_cleanup noexec, noexec, DESTARG ; Do final carry cleanup

	ad_epilog 2*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwysubnrzp1 ENDP

;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwyaddsubq1
	ad_prolog 0,0,rbx,rbp,rsi,ymm6,ymm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	eax, addcount1		; Load loop counter
uaddsublp2:
	mov	ebx, count2		; Load loop counter
uaddsublp:
	vmovapd	ymm1, [rcx]		; Load first number
	vaddpd	ymm0, ymm1, [rdx]	; Add in second number
	vsubpd	ymm1, ymm1, [rdx]	; Subtract out second number
	vmovapd	ymm3, [rcx+32]		; Load first number
	vaddpd	ymm2, ymm3, [rdx+32]	; Add in second number
	vsubpd	ymm3, ymm3, [rdx+32]	; Subtract out second number
	vmovapd	ymm5, [rcx+64]		; Load first number
	vaddpd	ymm4, ymm5, [rdx+64]	; Add in second number
	vsubpd	ymm5, ymm5, [rdx+64]	; Subtract out second number
	vmovapd	ymm7, [rcx+96]		; Load first number
	vaddpd	ymm6, ymm7, [rdx+96]	; Add in second number
	vsubpd	ymm7, ymm7, [rdx+96]	; Subtract out second number
	ystore	[rsi], ymm0		; Save result
	ystore	[rbp], ymm1		; Save result
	ystore	[rsi+32], ymm2		; Save result
	ystore	[rbp+32], ymm3		; Save result
	ystore	[rsi+64], ymm4		; Save result
	ystore	[rbp+64], ymm5		; Save result
	ystore	[rsi+96], ymm6		; Save result
	ystore	[rbp+96], ymm7		; Save result
	bump	rcx, 128		; Next source
	bump	rdx, 128		; Next source
	bump	rsi, 128		; Next dest
	bump	rbp, 128		; Next dest
	sub	ebx, 1			; Check loop counter
	jnz	uaddsublp		; Loop if necessary
	bump	rcx, 64			; Pad 64 bytes
	bump	rdx, 64			; Pad 64 bytes
	bump	rsi, 64			; Pad 64 bytes
	bump	rbp, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	uaddsublp2		; Loop if necessary
	ad_epilog 0,0,rbx,rbp,rsi,ymm6,ymm7
gwyaddsubq1 ENDP

;;
;; Add and subtract two numbers with carry propagation (eight different versions)
;;

saved_dest1_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_dest2_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_ttp_ptr	EQU	PPTR [rsp+first_local+2*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+3*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+3*SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+3*SZPTR+8]

	; Base-2, irrational, not zero-padded
PROCFL	gwyaddsub1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbx, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
b2addsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_ttp_ptr, rbx	; remember rbx for ynorm_addsub_1d_mid_cleanup
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_mid_cleanup
b2addsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2addsublp:
	ynorm_addsub_1d exec, exec	; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2addsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2addsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	b2addsublp2		; Loop til done
b2addsubsecdn:

	xchg	rbx, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	ynorm_addsub_1d_mid_cleanup exec, exec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	b2addsubsec

	ynorm_op_1d_cleanup exec, exec, DESTARG	; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_cleanup exec, exec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsub1 ENDP

	; Base-2, rational, not zero-padded
PROCFL	gwyaddsubr1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
b2raddsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_mid_cleanup
b2raddsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2raddsublp:
	ynorm_addsub_1d noexec, exec	; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2raddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2raddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	b2raddsublp2		; Loop til done
b2raddsubsecdn:

	ynorm_addsub_1d_mid_cleanup noexec, exec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2raddsubsec

	ynorm_op_1d_cleanup noexec, exec, DESTARG ; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_cleanup noexec, exec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubr1 ENDP

	; Base-2, irrational, zero-padded
PROCFL	gwyaddsubzp1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbx, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
b2zpaddsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_ttp_ptr, rbx	; remember rbx for ynorm_addsub_1d_zpad_mid_cleanup
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_zpad_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_zpad_mid_cleanup
b2zpaddsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2zpaddsublp:
	ynorm_addsub_1d_zpad exec, exec	; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2zpaddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2zpaddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	b2zpaddsublp2		; Loop til done
b2zpaddsubsecdn:

	xchg	rbx, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	ynorm_addsub_1d_zpad_mid_cleanup exec, exec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	b2zpaddsubsec

	ynorm_op_1d_zpad_cleanup exec, exec, DESTARG	; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_zpad_cleanup exec, exec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubzp1 ENDP

	; Base-2, rational, zero-padded
PROCFL	gwyaddsubrzp1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vmovapd	ymm2, YMM_BIGVAL	; Start process with no carry
	vmovapd	ymm3, ymm2
	vmovapd	ymm6, ymm2
	vmovapd	ymm7, ymm2
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
b2rzpaddsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_zpad_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_zpad_mid_cleanup
b2rzpaddsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
b2rzpaddsublp:
	ynorm_addsub_1d_zpad noexec, exec ; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	b2rzpaddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short b2rzpaddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	b2rzpaddsublp2		; Loop til done
b2rzpaddsubsecdn:

	ynorm_addsub_1d_zpad_mid_cleanup noexec, exec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	b2rzpaddsubsec

	ynorm_op_1d_zpad_cleanup noexec, exec, DESTARG ; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_zpad_cleanup noexec, exec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubrzp1 ENDP

	; Non base-2, irrational, not zero-padded
PROCFL	gwyaddsubn1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbx, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nb2addsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_ttp_ptr, rbx	; remember rbx for ynorm_addsub_1d_mid_cleanup
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_mid_cleanup
nb2addsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2addsublp:
	ynorm_addsub_1d exec, noexec	; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2addsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2addsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	nb2addsublp2		; Loop til done
nb2addsubsecdn:

	xchg	rbx, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	ynorm_addsub_1d_mid_cleanup exec, noexec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	nb2addsubsec

	ynorm_op_1d_cleanup exec, noexec, DESTARG ; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_cleanup exec, noexec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubn1 ENDP

	; Non base-2, rational, not zero-padded
PROCFL	gwyaddsubnr1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
nb2raddsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_mid_cleanup
nb2raddsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2raddsublp:
	ynorm_addsub_1d noexec, noexec	; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2raddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2raddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	nb2raddsublp2		; Loop til done
nb2raddsubsecdn:

	ynorm_addsub_1d_mid_cleanup noexec, noexec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2raddsubsec

	ynorm_op_1d_cleanup noexec, noexec, DESTARG ; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_cleanup noexec, noexec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubnr1 ENDP

	; Non base-2, irrational, zero-padded
PROCFL	gwyaddsubnzp1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
	mov	rbx, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nb2zpaddsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_ttp_ptr, rbx	; remember rbx for ynorm_addsub_1d_zpad_mid_cleanup
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_zpad_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_zpad_mid_cleanup
nb2zpaddsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2zpaddsublp:
	ynorm_addsub_1d_zpad exec, noexec ; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2zpaddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2zpaddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	nb2zpaddsublp2		; Loop til done
nb2zpaddsubsecdn:

	xchg	rbx, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	ynorm_addsub_1d_zpad_mid_cleanup exec, noexec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	mov	rbx, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	nb2zpaddsubsec

	ynorm_op_1d_zpad_cleanup exec, noexec, DESTARG ; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_zpad_cleanup exec, noexec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubnzp1 ENDP

	; Non base-2, rational, zero-padded
PROCFL	gwyaddsubnrzp1
	ad_prolog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	vxorpd	ymm2, ymm2, ymm2	; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	vxorpd	ymm6, ymm6, ymm6
	vxorpd	ymm7, ymm7, ymm7
	mov	eax, 4			; 4 sections
	mov	loopcount1, eax		; Save loop counter
nb2rzpaddsubsec:
	mov	eax, normcount1		; Get count of cache lines in a section
	mov	loopcount2, eax
	mov	saved_dest2_ptr, rbp	; remember dest #2 for ynorm_addsub_1d_zpad_mid_cleanup
	mov	saved_dest1_ptr, rsi	; remember dest #1 for ynorm_addsub_1d_zpad_mid_cleanup
nb2rzpaddsublp2:
	mov	eax, count1		; Get count of cache lines before a padding occurs
	mov	loopcount3, eax
nb2rzpaddsublp:
	ynorm_addsub_1d_zpad noexec, noexec ; Add/sub and normalize 8 values
	sub	loopcount3, 1		; Decrement loop counter
	jnz	nb2rzpaddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short nb2rzpaddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	rbp, [rbp+64]		; Pad by 64 bytes
	jnz	nb2rzpaddsublp2		; Loop til done
nb2rzpaddsubsecdn:

	ynorm_addsub_1d_zpad_mid_cleanup noexec, noexec, saved_dest1_ptr, saved_dest2_ptr ; Rotate and add in carries
	sub	loopcount1, 1		; Test section counter
	jnz	nb2rzpaddsubsec

	ynorm_op_1d_zpad_cleanup noexec, noexec, DESTARG ; Do final carry cleanup of result #1

	vmovapd	ymm2, ymm6		; Load carry
	vmovapd	ymm3, ymm7		; Load carry
	ynorm_op_1d_zpad_cleanup noexec, noexec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 3*SZPTR+12,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsubnrzp1 ENDP

;;
;; Copy one number and zero some low order words.
;;

PROCFL	gwycopyzero1
	ad_prolog 0,0,rsi,rdi
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

	mov	eax, addcount1		; Load loop counter
cz2:	mov	edx, count2		; Load loop counter
	add	rdx, rdx
cz1:	ycopyzero			; Copy/zero 8 values
	bump	rsi, 64			; Next source
	bump	rdi, 64			; Next dest
	bump	rcx, 64			; Next compare offset
	sub	edx, 1			; Check loop counter
	jnz	short cz1		; Loop if necessary
	bump	rsi, 64			; Pad 64 bytes
	bump	rdi, 64			; Pad 64 bytes
	bump	rcx, 64			; Pad 64 bytes
	sub	eax, 1			; Test loop counter
	jnz	cz2			; Loop if necessary
	ad_epilog 0,0,rsi,rdi
gwycopyzero1 ENDP

;;
;; Add in a small number with carry propagation (two different versions)
;;

	; Base-2 version
PROCFL	gwyadds1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	vmovsd	xmm7, DBLARG		; Small addin value
	ynorm_smalladd_1d exec
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyadds1 ENDP

	; Non base-2 version
PROCFL	gwyaddsn1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG		; Address of destination
	vmovsd	xmm7, DBLARG		; Small addin value
	ynorm_smalladd_1d noexec
	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwyaddsn1 ENDP

;;
;; Multiply a number by a small value (eight versions)
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_ttp_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_biglit_ptr EQU	PPTR [rsp+first_local+2*SZPTR]
saved_reg1	EQU	PPTR [rsp+first_local+3*SZPTR]
saved_reg2	EQU	PPTR [rsp+first_local+4*SZPTR]
saved_reg3	EQU	PPTR [rsp+first_local+5*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+6*SZPTR]
loopcount3	EQU	DPTR [rsp+first_local+6*SZPTR+4]

	; Base-2, irrational, not zero-padded
PROCFL	gwymuls1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vmulpd	ymm7, ymm7, YMM_NORM012_FF	; Mul by FFTLEN/2
	vmovapd	ymm2, YMM_BIGVAL		; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
b2mulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for ynorm_smallmul_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
b2mullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
b2mullp:
	ynorm_smallmul_1d exec, exec		; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	b2mullp				; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short b2mulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	b2mullp2			; Loop til done
b2mulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	mov	saved_reg2, rdi			; Save big/lit pointer
	mov	saved_reg3, rbp			; Save ttp pointer
	ynorm_smallmul_1d_mid_cleanup exec, exec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	mov	rdi, saved_reg2			; Restore big/lit pointer
	mov	rbp, saved_reg3			; Restore ttp pointer
	sub	loopcount1, 1			; Test section counter
	jnz	b2mulsec			; Do another section

	ynorm_smallmul_1d_cleanup exec, exec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymuls1 ENDP

	; Base-2, rational, not zero-padded
PROCFL	gwymulsr1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vmovapd	ymm2, YMM_BIGVAL		; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
b2rmulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
b2rmullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
b2rmullp:
	ynorm_smallmul_1d noexec, exec		; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	b2rmullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short b2rmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	b2rmullp2			; Loop til done
b2rmulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	ynorm_smallmul_1d_mid_cleanup noexec, exec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	sub	loopcount1, 1			; Test section counter
	jnz	b2rmulsec			; Do another section

	ynorm_smallmul_1d_cleanup noexec, exec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulsr1 ENDP

	; Base-2, irrational, zero-padded
PROCFL	gwymulszp1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vmulpd	ymm7, ymm7, YMM_NORM012_FF	; Mul by FFTLEN/2
	vmovapd	ymm2, YMM_BIGVAL		; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
b2zpmulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for ynorm_smallmul_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
b2zpmullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
b2zpmullp:
	ynorm_smallmul_1d_zpad exec, exec	; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	b2zpmullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short b2zpmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	b2zpmullp2			; Loop til done
b2zpmulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	mov	saved_reg2, rdi			; Save big/lit pointer
	mov	saved_reg3, rbp			; Save ttp pointer
	ynorm_smallmul_1d_zpad_mid_cleanup exec, exec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	mov	rdi, saved_reg2			; Restore big/lit pointer
	mov	rbp, saved_reg3			; Restore ttp pointer
	sub	loopcount1, 1			; Test section counter
	jnz	b2zpmulsec			; Do another section

	ynorm_smallmul_1d_zpad_cleanup exec, exec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulszp1 ENDP

	; Base-2, rational, zero-padded
PROCFL	gwymulsrzp1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vmovapd	ymm2, YMM_BIGVAL		; Start process with no carry
	vmovapd	ymm3, ymm2
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
b2rzpmulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
b2rzpmullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
b2rzpmullp:
	ynorm_smallmul_1d_zpad noexec, exec	; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	b2rzpmullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short b2rzpmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	b2rzpmullp2			; Loop til done
b2rzpmulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	ynorm_smallmul_1d_zpad_mid_cleanup noexec, exec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	sub	loopcount1, 1			; Test section counter
	jnz	b2rzpmulsec			; Do another section

	ynorm_smallmul_1d_zpad_cleanup noexec, exec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulsrzp1 ENDP

	; Non base-2, irrational, not zero-padded
PROCFL	gwymulsn1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vmulpd	ymm7, ymm7, YMM_NORM012_FF	; Mul by FFTLEN/2
	vxorpd	ymm2, ymm2, ymm2		; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
nb2mulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for ynorm_smallmul_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
nb2mullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
nb2mullp:
	ynorm_smallmul_1d exec, noexec		; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	nb2mullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short nb2mulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	nb2mullp2			; Loop til done
nb2mulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	mov	saved_reg2, rdi			; Save big/lit pointer
	mov	saved_reg3, rbp			; Save ttp pointer
	ynorm_smallmul_1d_mid_cleanup exec, noexec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	mov	rdi, saved_reg2			; Restore big/lit pointer
	mov	rbp, saved_reg3			; Restore ttp pointer
	sub	loopcount1, 1			; Test section counter
	jnz	nb2mulsec			; Do another section

	ynorm_smallmul_1d_cleanup exec, noexec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulsn1 ENDP

	; Non base-2, rational, not zero-padded
PROCFL	gwymulsnr1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vxorpd	ymm2, ymm2, ymm2		; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
nb2rmulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
nb2rmullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
nb2rmullp:
	ynorm_smallmul_1d noexec, noexec	; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	nb2rmullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short nb2rmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	nb2rmullp2			; Loop til done
nb2rmulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	ynorm_smallmul_1d_mid_cleanup noexec, noexec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	sub	loopcount1, 1			; Test section counter
	jnz	nb2rmulsec			; Do another section

	ynorm_smallmul_1d_cleanup noexec, noexec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulsnr1 ENDP

	; Non base-2, irrational, zero-padded
PROCFL	gwymulsnzp1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vmulpd	ymm7, ymm7, YMM_NORM012_FF	; Mul by FFTLEN/2
	vxorpd	ymm2, ymm2, ymm2		; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
nb2zpmulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for ynorm_smallmul_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
nb2zpmullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
nb2zpmullp:
	ynorm_smallmul_1d_zpad exec, noexec	; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	nb2zpmullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short nb2zpmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	nb2zpmullp2			; Loop til done
nb2zpmulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	mov	saved_reg2, rdi			; Save big/lit pointer
	mov	saved_reg3, rbp			; Save ttp pointer
	ynorm_smallmul_1d_zpad_mid_cleanup exec, noexec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	mov	rdi, saved_reg2			; Restore big/lit pointer
	mov	rbp, saved_reg3			; Restore ttp pointer
	sub	loopcount1, 1			; Test section counter
	jnz	nb2zpmulsec			; Do another section

	ynorm_smallmul_1d_zpad_cleanup exec, noexec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulsnzp1 ENDP

	; Non base-2, rational, zero-padded
PROCFL	gwymulsnrzp1
	ad_prolog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rsi, DESTARG			; Address of destination
	vbroadcastsd ymm7, DBLARG		; Load small multiplier value
	vxorpd	ymm2, ymm2, ymm2		; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
	mov	eax, 4				; 4 sections
	mov	loopcount1, eax			; Save counter
nb2rzpmulsec:
	mov	saved_dest_ptr, rsi		; remember rsi for ynorm_smallmul_1d_mid_cleanup
	mov	eax, normcount1			; Load count of cache lines in a section
	mov	loopcount3, eax
nb2rzpmullp2:
	mov	ebx, count1			; Get count of cache lines before a padding occurs
nb2rzpmullp:
	ynorm_smallmul_1d_zpad noexec, noexec	; Mul and normalize 8 values
	sub	ebx, 1				; Decrement loop counter
	jnz	nb2rzpmullp			; Loop til done
	sub	loopcount3, 1			; Decrement loop counter
	js	short nb2rzpmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	nb2rzpmullp2			; Loop til done
nb2rzpmulsecdn:

	mov	saved_reg1, rsi			; Save FFT data addr
	ynorm_smallmul_1d_zpad_mid_cleanup noexec, noexec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr ; Rotate carries and add in carries
	mov	rsi, saved_reg1			; Restore FFT data addr
	sub	loopcount1, 1			; Test section counter
	jnz	nb2rzpmulsec			; Do another section

	ynorm_smallmul_1d_zpad_cleanup noexec, noexec, DESTARG ; Do final carry propagations

	ad_epilog 6*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwymulsnrzp1 ENDP

ENDIF

;;
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_reg1		EQU	PPTR [rsp+first_local+0*SZPTR]
saved_reg2		EQU	PPTR [rsp+first_local+1*SZPTR]
saved_reg3		EQU	PPTR [rsp+first_local+2*SZPTR]
section_srcptr		EQU	PPTR [rsp+first_local+3*SZPTR]
section_biglitptr	EQU	PPTR [rsp+first_local+4*SZPTR]
section_ttpptr		EQU	PPTR [rsp+first_local+5*SZPTR]
loopcount1		EQU	DPTR [rsp+first_local+6*SZPTR]
loopcount2		EQU	DPTR [rsp+first_local+6*SZPTR+4]

inorm	MACRO	lab, ttp, zero, echk, const, base2
	LOCAL	ilp0, ilp1, ilp2, ilpsecdn
	PROCFLP	lab
	int_prolog 6*SZPTR+8,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	vmovsd	xmm0, ADDIN_VALUE	;; Get the addin value
no zero	vsubpd	ymm7, ymm7, ymm0	;; Do not include addin in sumout
no zero	vaddsd	xmm0, xmm0, Q [rsi][rdi] ;; Add in the FFT value
no zero	vmovsd	Q [rsi][rdi], xmm0	;; Save the new value
base2	vmovapd	ymm2, YMM_BIGVAL	;; Start process with no carry
no base2 vxorpd	ymm2, ymm2, ymm2	;; Start process with no carry
	vmovapd	ymm3, ymm2
echk	vxorpd	ymm6, ymm6, ymm6	;; Clear maximum error
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	ynorm_1d_preload ttp, base2, zero, echk, const ;; Preload useful constants
	mov	eax, 4			;; Load loop counter (4 sections)
	mov	loopcount1, eax
ilp0:	mov	section_srcptr, rsi	;; remember rsi for ynorm_1d_mid_cleanup
ttp	mov	section_biglitptr, rdi	;; remember rdi for ynorm_1d_mid_cleanup
ttp	mov	section_ttpptr, rbp	;; remember rbp for ynorm_1d_mid_cleanup
	mov	eax, normcount1		;; Load loop counter
	mov	loopcount2, eax
ilp2:	mov	ebx, count1		;; Get count of cache lines before a padding occurs
ilp1:	ynorm_1d ttp, base2, zero, echk, const ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rbp, 128		;; Next set of 8 multipliers
ttp	bump	rdi, 2			;; Next big/little flags
	sub	ebx, 1			;; Test loop counter
	jnz	ilp1			;; Loop til done
	sub	loopcount2, 1		;; Test loop counter
	js	short ilpsecdn		;; Loop til done
	lea	rsi, [rsi+64]		;; Pad 64 bytes, preserve flags
	jnz	ilp2			;; Loop til done
ilpsecdn:
	mov	saved_reg1, rsi		;; Save FFT data addr
ttp	mov	saved_reg2, rdi		;; Save big/lit pointer
ttp	mov	saved_reg3, rbp		;; Save ttp pointer
	ynorm_1d_mid_cleanup ttp, base2, zero, section_srcptr, section_biglitptr, section_ttpptr ;; Rotate carries and add in carries
	mov	rsi, saved_reg1		;; Restore FFT data addr
ttp	mov	rdi, saved_reg2		;; Restore big/lit pointer
ttp	mov	rbp, saved_reg3		;; Restore ttp pointer
	sub	loopcount1, 1		;; Test section counter
	jnz	ilp0
echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error
zero ttp		jmp	zdn	;; Go to zero upper half irrational end code
zero no ttp		jmp	zrdn	;; Go to zero upper half rational end code
no base2 ttp		jmp	nb2dn	;; Go to non-base2 irrational end code
no base2 no ttp		jmp	nb2rdn	;; Go to non-base2 rational end code
base2 no zero ttp	jmp	b2dn	;; Go to base2 irrational end code
base2 no zero no ttp	jmp	b2rdn	;; Go to base2 rational end code
	ENDPP lab
	ENDM

zpnorm	MACRO	lab, ttp, echk, const, base2, khi, c1, cm1
	LOCAL	ilp0, ilp1, ilp2, ilpsecdn
	PROCFLP	lab
	int_prolog 6*SZPTR+8,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
base2	vmovapd	ymm2, YMM_BIGVAL	;; Start process with no carry
no base2 vxorpd	ymm2, ymm2, ymm2	;; Start process with no carry
	vxorpd	ymm3, ymm3, ymm3
echk	vxorpd	ymm6, ymm6, ymm6	;; Clear maximum error
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	ynorm_1d_zpad_preload ttp, base2, echk, const, khi, c1, cm1 ;; Preload useful constants
	mov	eax, 4			;; Do 4 sections
	mov	loopcount1, eax		;; Save loop counter
ilp0:	mov	section_srcptr, rsi	;; remember rsi for ynorm_1d_mid_cleanup
ttp	mov	section_biglitptr, rdi	;; remember rdi for ynorm_1d_mid_cleanup
ttp	mov	section_ttpptr, rbp	;; remember rbp for ynorm_1d_mid_cleanup
	mov	eax, normcount1		;; Get count of padded groups in a section
	mov	loopcount2, eax
ilp2:	mov	ebx, count1		;; Get count of cache lines before a padding occurs
ilp1:	ynorm_1d_zpad ttp, base2, echk, const, khi, c1, cm1 ;; Normalize 8 values
	bump	rsi, 64			;; Next cache line
ttp	bump	rbp, 64			;; Next set of 8 multipliers
ttp	bump	rdi, 1			;; Next big/little flags
	sub	ebx, 1			;; Test loop counter
	jnz	ilp1			;; Loop til done
	sub	loopcount2, 1		;; Test loop counter
	js	short ilpsecdn		;; Loop til done
	lea	rsi, [rsi+64]		;; Pad 64 bytes, preserve flags
	jnz	ilp2			;; Loop til done
ilpsecdn:
	mov	saved_reg1, rsi		;; Save FFT data addr
ttp	mov	saved_reg2, rdi		;; Save big/lit pointer
ttp	mov	saved_reg3, rbp		;; Save ttp pointer
	ynorm_1d_zpad_mid_cleanup ttp, base2, const, section_srcptr, section_biglitptr, section_ttpptr ;; Rotate carries and add in carries
	mov	rsi, saved_reg1		;; Restore FFT data addr
ttp	mov	rdi, saved_reg2		;; Restore big/lit pointer
ttp	mov	rbp, saved_reg3		;; Restore ttp pointer
	sub	loopcount1, 1		;; Test section counter
	jnz	ilp0
echk	ystore	YMM_MAXERR, ymm6	;; Save maximum error
no base2 const jmp nb2zpcdn		;; Go to zero padded FFT end code
no base2 no const jmp nb2zpdn		;; Go to non-base2 end code
base2 const jmp	zpcdn			;; Go to zero padded FFT end code
base2 no const jmp zpdn			;; Go to zero padded FFT end code
	ENDPP lab
	ENDM

; The many different normalization routines.  One for each valid combination of
; rational/irrational, zeroing/no zeroing, error check/no error check,
; mul by const/no mul by const, base2 / other than base 2

	inorm	yr1, noexec, noexec, noexec, noexec, exec
	inorm	yr1e, noexec, noexec, exec, noexec, exec
	inorm	yr1c, noexec, noexec, noexec, exec, exec
	inorm	yr1ec, noexec, noexec, exec, exec, exec
	inorm	yr1z, noexec, exec, noexec, noexec, exec
	inorm	yr1ze, noexec, exec, exec, noexec, exec
	inorm	yi1, exec, noexec, noexec, noexec, exec
	inorm	yi1e, exec, noexec, exec, noexec, exec
	inorm	yi1c, exec, noexec, noexec, exec, exec
	inorm	yi1ec, exec, noexec, exec, exec, exec
	inorm	yi1z, exec, exec, noexec, noexec, exec
	inorm	yi1ze, exec, exec, exec, noexec, exec
	zpnorm	yr1zp, noexec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	yr1zpc1, noexec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	yr1zpcm1, noexec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	yr1zpe, noexec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	yr1zpec1, noexec, exec, noexec, exec, exec, exec, noexec
	zpnorm	yr1zpecm1, noexec, exec, noexec, exec, exec, noexec, exec
	zpnorm	yr1zpc, noexec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	yr1zpec, noexec, exec, exec, exec, exec, noexec, noexec
	zpnorm	yi1zp, exec, noexec, noexec, exec, exec, noexec, noexec
	zpnorm	yi1zpc1, exec, noexec, noexec, exec, exec, exec, noexec
	zpnorm	yi1zpcm1, exec, noexec, noexec, exec, exec, noexec, exec
	zpnorm	yi1zpe, exec, exec, noexec, exec, exec, noexec, noexec
	zpnorm	yi1zpec1, exec, exec, noexec, exec, exec, exec, noexec
	zpnorm	yi1zpecm1, exec, exec, noexec, exec, exec, noexec, exec
	zpnorm	yi1zpc, exec, noexec, exec, exec, exec, noexec, noexec
	zpnorm	yi1zpec, exec, exec, exec, exec, exec, noexec, noexec
	zpnorm	yr1zpk, noexec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	yr1zpkc1, noexec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	yr1zpkcm1, noexec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	yr1zpek, noexec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	yr1zpekc1, noexec, exec, noexec, exec, noexec, exec, noexec
	zpnorm	yr1zpekcm1, noexec, exec, noexec, exec, noexec, noexec, exec
	zpnorm	yr1zpck, noexec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	yr1zpeck, noexec, exec, exec, exec, noexec, noexec, noexec
	zpnorm	yi1zpk, exec, noexec, noexec, exec, noexec, noexec, noexec
	zpnorm	yi1zpkc1, exec, noexec, noexec, exec, noexec, exec, noexec
	zpnorm	yi1zpkcm1, exec, noexec, noexec, exec, noexec, noexec, exec
	zpnorm	yi1zpek, exec, exec, noexec, exec, noexec, noexec, noexec
	zpnorm	yi1zpekc1, exec, exec, noexec, exec, noexec, exec, noexec
	zpnorm	yi1zpekcm1, exec, exec, noexec, exec, noexec, noexec, exec
	zpnorm	yi1zpck, exec, noexec, exec, exec, noexec, noexec, noexec
	zpnorm	yi1zpeck, exec, exec, exec, exec, noexec, noexec, noexec

	inorm	yr1b, noexec, noexec, noexec, noexec, noexec
	inorm	yr1eb, noexec, noexec, exec, noexec, noexec
	inorm	yr1cb, noexec, noexec, noexec, exec, noexec
	inorm	yr1ecb, noexec, noexec, exec, exec, noexec
	inorm	yi1b, exec, noexec, noexec, noexec, noexec
	inorm	yi1eb, exec, noexec, exec, noexec, noexec
	inorm	yi1cb, exec, noexec, noexec, exec, noexec
	inorm	yi1ecb, exec, noexec, exec, exec, noexec
	zpnorm	yr1zpb, noexec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	yr1zpbc1, noexec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	yr1zpbcm1, noexec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	yr1zpeb, noexec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	yr1zpebc1, noexec, exec, noexec, noexec, exec, exec, noexec
	zpnorm	yr1zpebcm1, noexec, exec, noexec, noexec, exec, noexec, exec
	zpnorm	yr1zpcb, noexec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	yr1zpecb, noexec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	yi1zpb, exec, noexec, noexec, noexec, exec, noexec, noexec
	zpnorm	yi1zpbc1, exec, noexec, noexec, noexec, exec, exec, noexec
	zpnorm	yi1zpbcm1, exec, noexec, noexec, noexec, exec, noexec, exec
	zpnorm	yi1zpeb, exec, exec, noexec, noexec, exec, noexec, noexec
	zpnorm	yi1zpebc1, exec, exec, noexec, noexec, exec, exec, noexec
	zpnorm	yi1zpebcm1, exec, exec, noexec, noexec, exec, noexec, exec
	zpnorm	yi1zpcb, exec, noexec, exec, noexec, exec, noexec, noexec
	zpnorm	yi1zpecb, exec, exec, exec, noexec, exec, noexec, noexec
	zpnorm	yr1zpbk, noexec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yr1zpbkc1, noexec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	yr1zpbkcm1, noexec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	yr1zpebk, noexec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yr1zpebkc1, noexec, exec, noexec, noexec, noexec, exec, noexec
	zpnorm	yr1zpebkcm1, noexec, exec, noexec, noexec, noexec, noexec, exec
	zpnorm	yr1zpcbk, noexec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	yr1zpecbk, noexec, exec, exec, noexec, noexec, noexec, noexec
	zpnorm	yi1zpbk, exec, noexec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yi1zpbkc1, exec, noexec, noexec, noexec, noexec, exec, noexec
	zpnorm	yi1zpbkcm1, exec, noexec, noexec, noexec, noexec, noexec, exec
	zpnorm	yi1zpebk, exec, exec, noexec, noexec, noexec, noexec, noexec
	zpnorm	yi1zpebkc1, exec, exec, noexec, noexec, noexec, exec, noexec
	zpnorm	yi1zpebkcm1, exec, exec, noexec, noexec, noexec, noexec, exec
	zpnorm	yi1zpcbk, exec, noexec, exec, noexec, noexec, noexec, noexec
	zpnorm	yi1zpecbk, exec, exec, exec, noexec, noexec, noexec, noexec

; Common code to finish off the one-pass FFTs normalization.  The
; Windows 64-bit ABI frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.

PROCFP	__common_ynorm1_end_code

	;; Dummy prolog to match normalization code
	int_prolog 6*SZPTR+8,0,0

; Finish off the normalization process by adding any carry to first values.
; Handle both the with and without two-to-phi array cases.

nb2zpdn:ynorm_1d_zpad_cleanup noexec, noexec ; Add in carries
	jmp	cmnend			; All done, go cleanup

nb2zpcdn:ynorm_1d_zpad_cleanup exec, noexec ; Add in carries
	jmp	cmnend			; All done, go cleanup

zpcdn:	ynorm_1d_zpad_cleanup exec, exec ; Add in carries
	jmp	cmnend			; All done, go cleanup

zpdn:	ynorm_1d_zpad_cleanup noexec, exec ; Add in carries
	jmp	cmnend			; All done, go cleanup

zdn:	mov	rsi, DESTARG		; Address of squared number
	ynorm_1d_cleanup exec, exec, exec ; Add in carries
	jmp	cmnend			; All done, go cleanup

zrdn:	mov	rsi, DESTARG		; Address of squared number
	ynorm_1d_cleanup noexec, exec, exec ; Add in carries
	jmp	cmnend			; All done, go cleanup

nb2dn:	mov	rsi, DESTARG		; Address of squared number
	ynorm_top_carry_1d exec, noexec	; Adjust top carry when k > 1
	ynorm_1d_cleanup exec, noexec, noexec ; Add in carries
	jmp	cmnend			; All done, go cleanup

nb2rdn:	mov	rsi, DESTARG		; Address of squared number
	ynorm_1d_cleanup noexec, noexec, noexec ; Add in carries
	jmp	cmnend			; All done, go cleanup

b2dn:	mov	rsi, DESTARG		; Address of squared number
	ynorm_top_carry_1d exec, exec	; Adjust top carry when k > 1
	ynorm_1d_cleanup exec, exec, noexec ; Add in carries
	jmp	cmnend			; All done, go cleanup

b2rdn:	mov	rsi, DESTARG		; Address of squared number
	ynorm_1d_cleanup noexec, exec, noexec ; Add in carries

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

cmnend:	mov	rsi, DESTARG		; Address of squared number
	ystore	YMM_TMP1, ymm7		; Add together the four partial sumouts
	vaddsd	xmm7, xmm7, Q YMM_TMP1[8]
	vaddsd	xmm7, xmm7, Q YMM_TMP1[16]
	vaddsd	xmm7, xmm7, Q YMM_TMP1[24]
	vmulsd	xmm7, xmm7, ttmp_ff_inv
	vmovsd	Q [rsi-24], xmm7	; Save sum of FFT outputs
	vmovsd	xmm6, MAXERR		; Compute new maximum error
	vmaxsd	xmm6, xmm6, Q YMM_MAXERR[0]
	vmaxsd	xmm6, xmm6, Q YMM_MAXERR[8]
	vmaxsd	xmm6, xmm6, Q YMM_MAXERR[16]
	vmaxsd	xmm6, xmm6, Q YMM_MAXERR[24]
	vmovsd	MAXERR, xmm6

; Return

	int_epilog 6*SZPTR+8,0,0
	ENDPP __common_ynorm1_end_code 

_TEXT	ENDS
END
