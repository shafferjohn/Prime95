; Copyright 2011-2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;

	TITLE   setup

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE znormal.mac
INCLUDE znormal_zpad.mac

_TEXT SEGMENT

;; Only assemble the add/sub/etc routines when compiled with the CORE architecture

IF (@INSTR(,%zarch,<SKX>) NE 0)

;; General register layout for routines below:

;; eax,ebx,r8	loop count registers
;; rcx, rdx	source pointers
;; rsi		dest pointer
;; r14		dest pointer #2 (for addsub)
;; r12		distance to second independent set of FFT data to normalize
;; rdi		big/lit flags ptr
;; rbp		two-to-phi data ptr
;; r9		saved dest ptr at section start
;; r10		saved two-to-phi ptr at section start
;; r15		saved dest ptr #2 at section start (for addsub)
;; r15		saved big/lit ptr at section start (for smallmul)

saved_dest_ptr	EQU	PPTR r9
saved_ttp_ptr	EQU	PPTR r10
saved_dest2_ptr	EQU	PPTR r15
saved_biglit_ptr EQU	PPTR r15
loopcount1	EQU	DPTR eax
loopcount2	EQU	DPTR ebx
loopcount3	EQU	DPTR r8d

;;
;; Add two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwzaddq1
	ad_prolog 0,0,rbx,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, ZMM_QUICK1_PAD_GRPS ; Number of paddings
uaddlp2:mov	ebx, ZMM_QUICK1_CL_TIL_PAD ; Cache lines until a padding occurs
uaddlp:	vmovapd	zmm0, [rdx+0*64]	; Load second number
	vaddpd	zmm0, zmm0, [rcx+0*64]	; Add in first number
	vmovapd	zmm1, [rdx+1*64]	; Load second number
	vaddpd	zmm1, zmm1, [rcx+1*64]	; Add in first number
	vmovapd	zmm2, [rdx+2*64]	; Load second number
	vaddpd	zmm2, zmm2, [rcx+2*64]	; Add in first number
	vmovapd	zmm3, [rdx+3*64]	; Load second number
	vaddpd	zmm3, zmm3, [rcx+3*64]	; Add in first number
	zstore	[rsi+0*64], zmm0	; Save result
	zstore	[rsi+1*64], zmm1	; Save result
	zstore	[rsi+2*64], zmm2	; Save result
	zstore	[rsi+3*64], zmm3	; Save result
	bump	rcx, 4*64		; Next source
	bump	rdx, 4*64		; Next source
	bump	rsi, 4*64		; Next dest
	sub	ebx, 4			; Completed 4 cache lines
	jnz	short uaddlp		; Loop if necessary
	bump	rcx, 64			; Pad 64 bytes
	bump	rdx, 64			; Pad 64 bytes
	bump	rsi, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	uaddlp2			; Loop if necessary
	ad_epilog 0,0,rbx,rsi
gwzaddq1 ENDP

;;
;; Add two numbers with carry propagation (eight different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzadd1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_preload exec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
addsec:	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_op_1d_mid_cleanup
addlp2:	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
addlp:	znorm_op_1d vaddpd, exec	; Add and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	addlp			; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short addsecdn		; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	addlp2			; Loop til done
addsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	rbp, saved_ttp_ptr	; Save/Restore multipliers pointer
	znorm_op_1d_mid_cleanup exec	; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	addsec

	mov	rsi, DESTARG		; Address of FFT data
	mov	rbp, norm_col_mults	; Address of the ttp/ttmp multipliers
	znorm_top_carry_1d		; Do a standard top carry
	znorm_op_1d_cleanup exec	; Do final carry cleanup

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
gwzadd1 ENDP

	; Rational, not zero-padded
PROCFL	gwzaddr1
	ad_prolog 0,0,rbx,rsi,r8,r9,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_preload noexec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
raddsec:mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
raddlp2:mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
raddlp:	znorm_op_1d vaddpd, noexec	; Add and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	raddlp			; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short raddsecdn		; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	raddlp2			; Loop til done
raddsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	znorm_op_1d_mid_cleanup noexec	; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	sub	loopcount1, 1		; Test section counter
	jnz	raddsec

	mov	rsi, DESTARG		; Address of FFT data
	znorm_op_1d_cleanup noexec	; Do final carry cleanup

	ad_epilog 0,0,rbx,rsi,r8,r9,r12,zmm6,zmm7
gwzaddr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzaddzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_zpad_preload exec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
zpaddsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_op_1d_mid_cleanup
zpaddlp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
zpaddlp:
	znorm_op_1d_zpad vaddpd, exec	; Add and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	zpaddlp			; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short zpaddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	zpaddlp2		; Loop til done
zpaddsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	rbp, saved_ttp_ptr	; Save/Restore multipliers pointer
	znorm_op_1d_zpad_mid_cleanup exec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	zpaddsec

	znorm_op_1d_zpad_final_cleanup_preload exec
	znorm_op_1d_zpad_final_cleanup exec, DESTARG ; Add in carry that does not wrap around
	znorm_op_1d_zpad_cleanup_preload exec
	znorm_op_1d_zpad_cleanup exec, DESTARG ; Do final carry cleanup

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
gwzaddzp1 ENDP

	; Rational, zero-padded
PROCFL	gwzaddrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_zpad_preload noexec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
rzpaddsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
rzpaddlp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
rzpaddlp:
	znorm_op_1d_zpad vaddpd, noexec	; Add and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	rzpaddlp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short rzpaddsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	rzpaddlp2		; Loop til done
rzpaddsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	znorm_op_1d_zpad_mid_cleanup noexec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	sub	loopcount1, 1		; Test section counter
	jnz	rzpaddsec

	znorm_op_1d_zpad_final_cleanup_preload noexec
	znorm_op_1d_zpad_final_cleanup noexec, DESTARG ; Add in carry that does not wrap around
	znorm_op_1d_zpad_cleanup_preload noexec
	znorm_op_1d_zpad_cleanup noexec, DESTARG ; Do final carry cleanup

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
gwzaddrzp1 ENDP


;;
;; Subtract two numbers without carry propagation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCFL	gwzsubq1
	ad_prolog 0,0,rbx,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, ZMM_QUICK1_PAD_GRPS ; Number of paddings
usublp2:mov	ebx, ZMM_QUICK1_CL_TIL_PAD ; Cache lines until a padding occurs
usublp:	vmovapd	zmm0, [rdx+0*64]	; Load second number
	vsubpd	zmm0, zmm0, [rcx+0*64]	; Subtract first number
	vmovapd	zmm1, [rdx+1*64]	; Load second number
	vsubpd	zmm1, zmm1, [rcx+1*64]	; Subtract first number
	vmovapd	zmm2, [rdx+2*64]	; Load second number
	vsubpd	zmm2, zmm2, [rcx+2*64]	; Subtract first number
	vmovapd	zmm3, [rdx+3*64]	; Load second number
	vsubpd	zmm3, zmm3, [rcx+3*64]	; Subtract first number
	zstore	[rsi+0*64], zmm0	; Save result
	zstore	[rsi+1*64], zmm1	; Save result
	zstore	[rsi+2*64], zmm2	; Save result
	zstore	[rsi+3*64], zmm3	; Save result
	bump	rcx, 4*64		; Next source
	bump	rdx, 4*64		; Next source
	bump	rsi, 4*64		; Next dest
	sub	ebx, 4			; Completed 4 cache lines
	jnz	short usublp		; Loop if necessary
	bump	rcx, 64			; Pad 64 bytes
	bump	rdx, 64			; Pad 64 bytes
	bump	rsi, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	usublp2			; Loop if necessary
	ad_epilog 0,0,rbx,rsi
gwzsubq1 ENDP

;;
;; Subtract two numbers with carry propagation (eight different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzsub1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_preload exec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
subsec:	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_op_1d_mid_cleanup
sublp2:	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
sublp:	znorm_op_1d vsubpd, exec	; Subtract and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	sublp			; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short subsecdn		; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	sublp2			; Loop til done
subsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	rbp, saved_ttp_ptr	; Save/Restore multipliers pointer
	znorm_op_1d_mid_cleanup exec	; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	subsec

	mov	rsi, DESTARG		; Address of FFT data
	mov	rbp, norm_col_mults	; Address of the ttp/ttmp multipliers
	znorm_top_carry_1d		; Do a standard top carry
	znorm_op_1d_cleanup exec	; Do final carry cleanup

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
gwzsub1 ENDP

	; Rational, not zero-padded
PROCFL	gwzsubr1
	ad_prolog 0,0,rbx,rsi,r8,r9,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_preload noexec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
rsubsec:mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
rsublp2:mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
rsublp:	znorm_op_1d vsubpd, noexec	; Subtract and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	rsublp			; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short rsubsecdn		; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	rsublp2			; Loop til done
rsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	znorm_op_1d_mid_cleanup noexec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	sub	loopcount1, 1		; Test section counter
	jnz	rsubsec

	mov	rsi, DESTARG		; Address of FFT data
	znorm_op_1d_cleanup noexec	; Do final carry cleanup

	ad_epilog 0,0,rbx,rsi,r8,r9,r12,zmm6,zmm7
gwzsubr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzsubzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_zpad_preload exec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
zpsubsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_op_1d_mid_cleanup
zpsublp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
zpsublp:
	znorm_op_1d_zpad vsubpd, exec	; Subtract and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	zpsublp			; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short zpsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	zpsublp2		; Loop til done
zpsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	rbp, saved_ttp_ptr	; Save/Restore multipliers pointer
	znorm_op_1d_zpad_mid_cleanup exec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	zpsubsec

	znorm_op_1d_zpad_final_cleanup_preload exec
	znorm_op_1d_zpad_final_cleanup exec, DESTARG ; Add in carry that does not wrap around
	znorm_op_1d_zpad_cleanup_preload exec
	znorm_op_1d_zpad_cleanup exec, DESTARG ; Do final carry cleanup

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
gwzsubzp1 ENDP

	; Rational, zero-padded
PROCFL	gwzsubrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_op_1d_zpad_preload noexec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
rzpsubsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember rsi for znorm_op_1d_mid_cleanup
rzpsublp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
rzpsublp:
	znorm_op_1d_zpad vsubpd, noexec ; Subtract and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	rzpsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short rzpsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	jnz	rzpsublp2		; Loop til done
rzpsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	znorm_op_1d_zpad_mid_cleanup noexec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	sub	loopcount1, 1		; Test section counter
	jnz	rzpsubsec

	znorm_op_1d_zpad_final_cleanup_preload noexec
	znorm_op_1d_zpad_final_cleanup noexec, DESTARG ; Add in carry that does not wrap around
	znorm_op_1d_zpad_cleanup_preload noexec
	znorm_op_1d_zpad_cleanup noexec, DESTARG ; Do final carry cleanup

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,zmm6,zmm7
gwzsubrzp1 ENDP


;;
;; Add and subtract two numbers without carry propagation.
;;

PROCFL	gwzaddsubq1
	ad_prolog 0,0,rbx,rsi,r14,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	r14, DEST2ARG	  	; Address of destination #2
	mov	eax, ZMM_QUICK1_PAD_GRPS ; Number of paddings
uaddsublp2:
	mov	ebx, ZMM_QUICK1_CL_TIL_PAD ; Cache lines until a padding occurs
uaddsublp:
	vmovapd	zmm1, [rcx+0*64]	; Load first number
	vaddpd	zmm0, zmm1, [rdx+0*64]	; Add in second number
	vsubpd	zmm1, zmm1, [rdx+0*64]	; Subtract out second number
	vmovapd	zmm3, [rcx+1*64]	; Load first number
	vaddpd	zmm2, zmm3, [rdx+1*64]	; Add in second number
	vsubpd	zmm3, zmm3, [rdx+1*64]	; Subtract out second number
	vmovapd	zmm5, [rcx+2*64]	; Load first number
	vaddpd	zmm4, zmm5, [rdx+2*64]	; Add in second number
	vsubpd	zmm5, zmm5, [rdx+2*64]	; Subtract out second number
	vmovapd	zmm7, [rcx+3*64]	; Load first number
	vaddpd	zmm6, zmm7, [rdx+3*64]	; Add in second number
	vsubpd	zmm7, zmm7, [rdx+3*64]	; Subtract out second number
	zstore	[rsi+0*64], zmm0	; Save result
	zstore	[r14+0*64], zmm1	; Save result
	zstore	[rsi+1*64], zmm2	; Save result
	zstore	[r14+1*64], zmm3	; Save result
	zstore	[rsi+2*64], zmm4	; Save result
	zstore	[r14+2*64], zmm5	; Save result
	zstore	[rsi+3*64], zmm6	; Save result
	zstore	[r14+3*64], zmm7	; Save result
	bump	rcx, 4*64		; Next source
	bump	rdx, 4*64		; Next source
	bump	rsi, 4*64		; Next dest
	bump	r14, 4*64		; Next dest
	sub	ebx, 4			; Completed 4 cache lines
	jnz	uaddsublp		; Loop if necessary
	bump	rcx, 64			; Pad 64 bytes
	bump	rdx, 64			; Pad 64 bytes
	bump	rsi, 64			; Pad 64 bytes
	bump	r14, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	uaddsublp2		; Loop if necessary
	ad_epilog 0,0,rbx,rsi,r14,zmm6,zmm7
gwzaddsubq1 ENDP

;;
;; Add and subtract two numbers with carry propagation (eight different versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzaddsub1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r14,r15,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of add destination
	mov	r14, DEST2ARG	  	; Address of sub destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_1d_preload exec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
addsubsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_addsub_1d_mid_cleanup
	mov	saved_dest_ptr, rsi	; remember dest #1 for znorm_addsub_1d_mid_cleanup
	mov	saved_dest2_ptr, r14	; remember dest #2 for znorm_addsub_1d_mid_cleanup
addsublp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
addsublp:
	znorm_addsub_1d exec		; Add/sub and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	addsublp		; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short addsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	r14, [r14+64]		; Pad by 64 bytes
	jnz	addsublp2		; Loop til done
addsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	r14, saved_dest2_ptr	; Save/Restore dest2 pointer
	xchg	rbp, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	znorm_addsub_1d_mid_cleanup exec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	r14, saved_dest2_ptr	; Restore dest2 pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	addsubsec

	mov	rsi, DESTARG		; Address of add FFT data
	mov	rbp, norm_col_mults	; Address of the ttp/ttmp multipliers
	znorm_top_carry_1d		; Do a standard top carry
	znorm_op_1d_cleanup exec	; Do final carry cleanup of add result
	vmovapd	zmm0, zmm4		; Copy carries
	vmovapd	zmm1, zmm5
	vmovapd	zmm2, zmm6
	vmovapd	zmm3, zmm7
	mov	rsi, DEST2ARG		; Address of sub FFT data
	znorm_top_carry_1d		; Do a standard top carry
	znorm_op_1d_cleanup exec	; Do final carry cleanup of sub result

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r14,r15,zmm6,zmm7
gwzaddsub1 ENDP

	; Rational, not zero-padded
PROCFL	gwzaddsubr1
	ad_prolog 0,0,rbx,rsi,r8,r9,r12,r14,r15,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of add destination
	mov	r14, DEST2ARG	  	; Address of sub destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_1d_preload noexec	; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
raddsubsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi	; remember dest #1 for znorm_addsub_1d_mid_cleanup
	mov	saved_dest2_ptr, r14	; remember dest #2 for znorm_addsub_1d_mid_cleanup
raddsublp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
raddsublp:
	znorm_addsub_1d noexec		; Add/sub and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	raddsublp		; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short raddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	r14, [r14+64]		; Pad by 64 bytes
	jnz	raddsublp2		; Loop til done
raddsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	r14, saved_dest2_ptr	; Save/Restore dest2 pointer
	znorm_addsub_1d_mid_cleanup noexec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	r14, saved_dest2_ptr	; Restore dest2 pointer
	sub	loopcount1, 1		; Test section counter
	jnz	raddsubsec

	mov	rsi, DESTARG		; Address of add FFT data
	znorm_op_1d_cleanup noexec	; Do final carry cleanup of add result
	vmovapd	zmm0, zmm4		; Copy carries
	vmovapd	zmm1, zmm5
	vmovapd	zmm2, zmm6
	vmovapd	zmm3, zmm7
	mov	rsi, DEST2ARG		; Address of sub FFT data
	znorm_op_1d_cleanup noexec	; Do final carry cleanup of sub result

	ad_epilog 0,0,rbx,rsi,r8,r9,r12,r14,r15,zmm6,zmm7
gwzaddsubr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzaddsubzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r14,r15,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of add destination
	mov	r14, DEST2ARG	  	; Address of sub destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_1d_zpad_preload exec ; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
zpaddsubsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_addsub_1d_mid_cleanup
	mov	saved_dest_ptr, rsi	; remember dest #1 for znorm_addsub_1d_mid_cleanup
	mov	saved_dest2_ptr, r14	; remember dest #2 for znorm_addsub_1d_mid_cleanup
zpaddsublp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
zpaddsublp:
	znorm_addsub_1d_zpad exec	; Add/sub and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	zpaddsublp		; Loop til done
	sub	loopcount2, 1		; Test count of pad groups (may be zero!)
	js	short zpaddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	r14, [r14+64]		; Pad by 64 bytes
	jnz	zpaddsublp2		; Loop til done
zpaddsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	r14, saved_dest2_ptr	; Save/Restore dest2 pointer
	xchg	rbp, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	znorm_addsub_1d_zpad_mid_cleanup exec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	r14, saved_dest2_ptr	; Restore dest2 pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	zpaddsubsec

	znorm_addsub_1d_zpad_final_cleanup_preload exec
	znorm_addsub_1d_zpad_final_cleanup exec, DESTARG, DEST2ARG ; Add in carry that does not wrap around

	vmovsd	Q ZMM_TMPS+0*8, xmm5	; Save carry
	znorm_op_1d_zpad_cleanup_preload exec
	znorm_op_1d_zpad_cleanup exec, DESTARG ; Do final carry cleanup of result #1
	vmovsd	xmm1, Q ZMM_TMPS+0*8	; Restore carry
	znorm_op_1d_zpad_cleanup exec, DEST2ARG	; Do final carry cleanup of result #2

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r14,r15,zmm6,zmm7
gwzaddsubzp1 ENDP

	; Rational, zero-padded
PROCFL	gwzaddsubrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r14,r15,zmm6,zmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of add destination
	mov	r14, DEST2ARG	  	; Address of sub destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_addsub_1d_zpad_preload noexec ; Preload constants
	mov	loopcount1, ZMM_NORM1_SECTIONS ; Load loop counter (4 sections)
	mov	rbp, norm_col_mults	; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
rzpaddsubsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ; Number of subsections (padding groups)
	mov	saved_ttp_ptr, rbp	; remember rbp for znorm_addsub_1d_mid_cleanup
	mov	saved_dest_ptr, rsi	; remember dest #1 for znorm_addsub_1d_mid_cleanup
	mov	saved_dest2_ptr, r14	; remember dest #2 for znorm_addsub_1d_mid_cleanup
rzpaddsublp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
rzpaddsublp:
	znorm_addsub_1d_zpad noexec	; Add/sub and normalize 4 cache lines
	sub	loopcount3, 4		; Test cache line loop counter
	jnz	rzpaddsublp		; Loop til done
	sub	loopcount2, 1		; Decrement loop counter
	js	short rzpaddsubsecdn	; Loop til done
	lea	rcx, [rcx+64]		; Pad by 64 bytes
	lea	rdx, [rdx+64]		; Pad by 64 bytes
	lea	rsi, [rsi+64]		; Pad by 64 bytes
	lea	r14, [r14+64]		; Pad by 64 bytes
	jnz	rzpaddsublp2		; Loop til done
rzpaddsubsecdn:
	xchg	rsi, saved_dest_ptr	; Save/Restore dest pointer
	xchg	r14, saved_dest2_ptr	; Save/Restore dest2 pointer
	xchg	rbp, saved_ttp_ptr	; Save/restore ttp/ttmp multipliers pointer
	znorm_addsub_1d_zpad_mid_cleanup noexec ; Rotate and add in carries
	mov	rsi, saved_dest_ptr	; Restore dest pointer
	mov	r14, saved_dest2_ptr	; Restore dest2 pointer
	mov	rbp, saved_ttp_ptr	; Restore multipliers pointer
	sub	loopcount1, 1		; Test section counter
	jnz	rzpaddsubsec

	znorm_addsub_1d_zpad_final_cleanup_preload noexec
	znorm_addsub_1d_zpad_final_cleanup noexec, DESTARG, DEST2ARG ; Add in carry that does not wrap around

	znorm_op_1d_zpad_cleanup_preload noexec
	vmovsd	Q ZMM_TMPS+0*8, xmm5	; Save carry
	znorm_op_1d_zpad_cleanup noexec, DESTARG ; Do final carry cleanup of result #1
	vmovsd	xmm1, Q ZMM_TMPS+0*8	; Restore carry
	znorm_op_1d_zpad_cleanup noexec, DEST2ARG ; Do final carry cleanup of result #2

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r14,r15,zmm6,zmm7
gwzaddsubrzp1 ENDP


;;
;; Copy one number zeroing some low order words.
;;

PROCFL	gwzcopyzero1
	ad_prolog 0,0,rbx,rsi,rdi,r12
	mov	rsi, SRCARG		; Address of first number
	mov	rdi, DESTARG		; Address of destination
	mov	r12, ZMM_QUARTER_DIST	; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	sub	ecx, ecx		; Count of 4-cache-line iterations before new masks are needed
	lea	rdx, COPYZERO		; Load address of offsets and masks computed in C code
	mov	eax, ZMM_QUICK1_PAD_GRPS ; Number of paddings
cpzlp2:	mov	ebx, ZMM_QUICK1_CL_TIL_PAD ; Cache lines until a padding occurs
cpzlp:	and	ecx, ecx		; Time to load new masks?
	jnz	short noload		; No, wait for offsets to match
	kmovw	k1, [rdx+0]		; Load 4 new masks
	kshiftrw k2, k1, 8
	kmovw	k3, [rdx+2]
	kshiftrw k4, k3, 8
	mov	ecx, [rdx+4]		; Count of 4-cache-line iterations these masks are valid
	bump	rdx, 8			; Move onto next set of masks
noload:	vmovapd	zmm0 {k1}{z}, [rsi]	; Load source data
	vmovapd	zmm1 {k2}{z}, [rsi+r12]
	vmovapd	zmm2 {k3}{z}, [rsi+64]
	vmovapd	zmm3 {k4}{z}, [rsi+r12+64]
	zstore	[rdi], zmm0		; Store destination data
	zstore	[rdi+r12], zmm1
	zstore	[rdi+64], zmm2
	zstore	[rdi+r12+64], zmm3
	bump	rcx, -1			; Decrement mask count
	bump	rsi, 2*64		; Next source
	bump	rdi, 2*64		; Next dest
	sub	ebx, 4			; Completed 4 cache lines
	jnz	short cpzlp		; Loop if necessary
	bump	rsi, 64			; Pad 64 bytes
	bump	rdi, 64			; Pad 64 bytes
	sub	eax, 1			; Check loop counter
	jnz	cpzlp2			; Loop if necessary
	ad_epilog 0,0,rbx,rsi,rdi,r12
gwzcopyzero1 ENDP

;;
;; Add in a small number with carry propagation
;;

PROCFL	gwzadds1
	ad_prolog 0,0,rbx,rbp,rsi
	mov	rsi, DESTARG		; Address of destination
	vmovsd	xmm5, DBLARG		; Small addin value
	mov	rbp, norm_col_mults	; Address of the two-to-phi multipliers
	znorm_common_op_1d_preload exec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smalladd_1d_preload exec
	znorm_smalladd_1d exec
	ad_epilog 0,0,rbx,rbp,rsi
gwzadds1 ENDP

PROCFL	gwzaddsr1
	ad_prolog 0,0,rbx,rsi
	mov	rsi, DESTARG		; Address of destination
	vmovsd	xmm5, DBLARG		; Small addin value
	znorm_common_op_1d_preload noexec ; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smalladd_1d_preload noexec
	znorm_smalladd_1d noexec
	ad_epilog 0,0,rbx,rsi
gwzaddsr1 ENDP

;;
;; Multiply a number by a small value (four versions)
;;

	; Irrational, not zero-padded
PROCFL	gwzmuls1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r15,zmm6,zmm7
	mov	rsi, DESTARG			; Address of destination
	mov	r12, ZMM_QUARTER_DIST		; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	znorm_common_op_1d_preload exec		; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_1d_preload exec		; Preload constants, init carry to zero
	mov	loopcount1, ZMM_NORM1_SECTIONS	; Load loop counter (4 sections)
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
mulsec:	mov	loopcount2, ZMM_NORM1_PAD_GRPS	; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi		; remember rsi for znorm_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for znorm_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for znorm_1d_mid_cleanup
mullp2:	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
mullp:	znorm_smallmul_1d exec			; Mul and normalize 4 cache lines
	sub	loopcount3, 4			; Test cache line loop counter
	jnz	mullp				; Loop til done
	sub	loopcount2, 1			; Test count of pad groups (may be zero!)
	js	short mulsecdn			; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	mullp2				; Loop til done
mulsecdn:
	push	rsi				; Save src/dest pointer
	push	rbp		 		; Save multipliers pointer
	push	rdi				; Save big/lit pointer
	znorm_smallmul_1d_mid_cleanup exec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr, loopcount2, loopcount3 ; Rotate carries and add in carries
	pop	rdi				; Restore big/lit pointer
	pop	rbp				; Restore multipliers pointer
	pop	rsi				; Restore src/dest pointer
	sub	loopcount1, 1			; Test section counter
	jnz	mulsec				; Do another section

	mov	rsi, DESTARG			; Address of FFT data
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
	; znorm_1d_cleanup_preload		; Not needed - same as znorm_smallmul_1d_preload
	znorm_top_carry_1d			; Adjust top carry when k > 1
	znorm_1d_cleanup exec, noexec		; Do final carry propagations

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r15,zmm6,zmm7
gwzmuls1 ENDP

	; Rational, not zero-padded
PROCFL	gwzmulsr1
	ad_prolog 0,0,rbx,rsi,r8,r9,r12,zmm6,zmm7
	mov	rsi, DESTARG			; Address of destination
	mov	r12, ZMM_QUARTER_DIST		; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	znorm_common_op_1d_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_1d_preload noexec	; Preload constants, init carry to zero
	mov	loopcount1, ZMM_NORM1_SECTIONS	; Load loop counter (4 sections)
rmulsec:mov	loopcount2, ZMM_NORM1_PAD_GRPS	; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi		; remember rsi for znorm_1d_mid_cleanup
rmullp2:mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
rmullp:	znorm_smallmul_1d noexec		; Mul and normalize 4 cache lines
	sub	loopcount3, 4			; Test cache line loop counter
	jnz	rmullp				; Loop til done
	sub	loopcount2, 1			; Test count of pad groups (may be zero!)
	js	short rmulsecdn			; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	rmullp2				; Loop til done
rmulsecdn:
	push	rsi				; Save src/dest pointer
	znorm_smallmul_1d_mid_cleanup noexec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr, loopcount2, loopcount3 ; Rotate carries and add in carries
	pop	rsi				; Restore src/dest pointer
	sub	loopcount1, 1			; Test section counter
	jnz	rmulsec				; Do another section

	mov	rsi, DESTARG			; Address of FFT data
	; znorm_1d_cleanup_preload		; Not needed - same as znorm_smallmul_1d_preload
	znorm_1d_cleanup noexec, noexec		; Do final carry propagations

	ad_epilog 0,0,rbx,rsi,r8,r9,r12,zmm6,zmm7
gwzmulsr1 ENDP

	; Irrational, zero-padded
PROCFL	gwzmulszp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r15,zmm6,zmm7
	mov	rsi, DESTARG			; Address of destination
	mov	r12, ZMM_QUARTER_DIST		; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	znorm_common_op_1d_preload exec		; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_1d_zpad_preload exec	; Preload constants, init carry to zero
	mov	loopcount1, ZMM_NORM1_SECTIONS	; Load loop counter (4 sections)
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
zpmulsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS	; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi		; remember rsi for znorm_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for znorm_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for znorm_1d_mid_cleanup
zpmullp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
zpmullp:
	znorm_smallmul_1d_zpad exec		; Mul and normalize 4 cache lines
	sub	loopcount3, 4			; Test cache line loop counter
	jnz	zpmullp				; Loop til done
	sub	loopcount2, 1			; Test count of pad groups (may be zero!)
	js	short zpmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	zpmullp2			; Loop til done
zpmulsecdn:
	push	rsi				; Save src/dest pointer
	push	rbp		 		; Save multipliers pointer
	push	rdi				; Save big/lit pointer
	znorm_smallmul_1d_zpad_mid_cleanup exec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr, loopcount2, loopcount3 ; Rotate carries and add in carries
	pop	rdi				; Restore big/lit pointer
	pop	rbp				; Restore multipliers pointer
	pop	rsi				; Restore src/dest pointer
	sub	loopcount1, 1			; Test section counter
	jnz	zpmulsec			; Do another section

	znorm_smallmul_1d_zpad_final_cleanup_preload exec
	znorm_smallmul_1d_zpad_final_cleanup exec, DESTARG ; Add in carry that does not wrap around

	znorm_op_1d_zpad_cleanup_preload exec
	znorm_op_1d_zpad_cleanup exec, DESTARG	; Do final carry propagations

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r15,zmm6,zmm7
gwzmulszp1 ENDP

	; Rational, zero-padded
PROCFL	gwzmulsrzp1
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r15,zmm6,zmm7
	mov	rsi, DESTARG			; Address of destination
	mov	r12, ZMM_QUARTER_DIST		; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	vbroadcastsd zmm31, DBLARG		; Load small multiplier value
	znorm_common_op_1d_preload noexec	; Preload constants common to add,sub,addsub,smalladd,smallmul,with and without zpad
	znorm_smallmul_1d_zpad_preload noexec	; Preload constants, init carry to zero
	mov	loopcount1, ZMM_NORM1_SECTIONS	; Load loop counter (4 sections)
	mov	rbp, norm_col_mults		; Address of the ttp/ttmp multipliers
	mov	rdi, norm_biglit_array		; Addr of the big/little flags array
rzpmulsec:
	mov	loopcount2, ZMM_NORM1_PAD_GRPS	; Number of subsections (padding groups)
	mov	saved_dest_ptr, rsi		; remember rsi for znorm_1d_mid_cleanup
	mov	saved_ttp_ptr, rbp		; remember rbp for znorm_1d_mid_cleanup
	mov	saved_biglit_ptr, rdi		; remember rdi for znorm_1d_mid_cleanup
rzpmullp2:
	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ; Count of cache lines before a padding occurs or subsection ends
rzpmullp:
	znorm_smallmul_1d_zpad noexec		; Mul and normalize 4 cache lines
	sub	loopcount3, 4			; Test cache line loop counter
	jnz	rzpmullp			; Loop til done
	sub	loopcount2, 1			; Test count of pad groups (may be zero!)
	js	short rzpmulsecdn		; Loop til done
	lea	rsi, [rsi+64]			; Pad 64 bytes
	jnz	rzpmullp2			; Loop til done
rzpmulsecdn:
	push	rsi				; Save src/dest pointer
	push	rbp		 		; Save multipliers pointer
	push	rdi				; Save big/lit pointer
	znorm_smallmul_1d_zpad_mid_cleanup noexec, saved_dest_ptr, saved_biglit_ptr, saved_ttp_ptr, loopcount2, loopcount3 ; Rotate carries and add in carries
	pop	rdi				; Restore big/lit pointer
	pop	rbp				; Restore multipliers pointer
	pop	rsi				; Restore src/dest pointer
	sub	loopcount1, 1			; Test section counter
	jnz	rzpmulsec			; Do another section

	znorm_smallmul_1d_zpad_final_cleanup_preload noexec
	znorm_smallmul_1d_zpad_final_cleanup noexec, DESTARG ; Add in carry that does not wrap around

	znorm_op_1d_zpad_cleanup_preload noexec
	znorm_op_1d_zpad_cleanup noexec, DESTARG ; Do final carry propagations

	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r12,r15,zmm6,zmm7
gwzmulsrzp1 ENDP

ENDIF 

;;
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_reg1		EQU	r8
saved_reg2		EQU	r9
saved_reg3		EQU	r10
section_srcptr		EQU	r15
section_biglitptr	EQU	r14
section_ttpptr		EQU	r13
loopcount1		EQU	eax
loopcount2		EQU	ebx
loopcount3		EQU	ecx

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	ilp0, ilp1, ilp2, ilpsecdn
	PROCFLP	lab
	int_prolog 0,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
	mov	r12, ZMM_QUARTER_DIST	;; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	vmovsd	xmm0, ADDIN_VALUE	;; Get the addin value
no zero	vsubsd	xmm31, xmm31, xmm0	;; Do not include addin in sumout
no zero	vaddsd	xmm0, xmm0, Q [rsi][rdi] ;; Add in the FFT value
no zero	vmovsd	Q [rsi][rdi], xmm0	;; Save the new value
	vbroadcastsd zmm0, ZMM_RNDVAL	;; Start process with no carry
	vmovapd	zmm1, zmm0
	vmovapd	zmm2, zmm0
	vmovapd	zmm3, zmm0
echk	vpxorq	zmm27, zmm27, zmm27	;; Clear 4 maximum error registers

ttp	mov	rbp, norm_col_mults	;; Addr of the multipliers
ttp	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	mov	loopcount1, ZMM_NORM1_SECTIONS ;; Load loop counter (4 sections)
ilp0:	znorm_1d_preload ttp, zero, echk, const ;; Preload useful constants
	vpxorq	zmm30, zmm30, zmm30	;; Clear 3 auxilary sumout registers (zmm31 was initialized in zsub_7_words)
	vmovapd	zmm29, zmm30
	vmovapd	zmm28, zmm30
echk	vmovapd	zmm26, zmm30		;; Clear 3 auxilary maximum error registers
echk	vmovapd	zmm25, zmm30
echk	vmovapd	zmm24, zmm30

				;; BUG/OPT - use half as many sumout/maxerr registers as we do in zpad code

	mov	section_srcptr, rsi	;; remember rsi for znorm_1d_mid_cleanup
ttp	mov	section_biglitptr, rdi	;; remember rdi for znorm_1d_mid_cleanup
ttp	mov	section_ttpptr, rbp	;; remember rbp for znorm_1d_mid_cleanup
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ;; Number of subsections (padding groups)
ilp2:	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ;; Count of cache lines before a padding occurs or subsection ends
ilp1:	znorm_1d ttp, zero, echk, const ;; Normalize 32 values (4 cache lines)
	bump	rsi, 2*64		;; Next set of 4 cache lines (2 rsi cache lines, 2 rsi+r12 cache lines)
ttp	bump	rbp, 4*128		;; Next set of multipliers
ttp	bump	rdi, 4*1		;; Next set of big/little flags
	sub	loopcount3, 4		;; Test cache line loop counter
	jnz	ilp1			;; Loop til done
	sub	loopcount2, 1		;; Test count of pad groups (may be zero!)
	js	short ilpsecdn		;; Branch if not padding
	lea	rsi, [rsi+64]		;; Pad 64 bytes, preserve flags
	jnz	ilp2			;; Loop til done
ilpsecdn:
	vaddpd	zmm31, zmm31, zmm30	;; Coalesce sumouts
	vaddpd	zmm29, zmm29, zmm28
	vaddpd	zmm31, zmm31, zmm29
echk	vmaxpd	zmm27, zmm27, zmm26	;; Coalesce maximum error
echk	vmaxpd	zmm25, zmm25, zmm24
echk	vmaxpd	zmm27, zmm27, zmm25
echk	zstore	ZMM_MAXERR, zmm27	;; Save maximum error	BUG/OPT - use only 2 sumout and maxerr registers (like zpad case)
								;; gives us 4 more free registers and makes sharing cleanup code easier
								;; (constants in the same registers)
	mov	saved_reg1, rsi		;; Save FFT data addr
ttp	mov	saved_reg2, rdi		;; Save big/lit pointer
ttp	mov	saved_reg3, rbp		;; Save ttp pointer
	znorm_1d_mid_cleanup_preload ttp
	znorm_1d_mid_cleanup ttp, zero, section_srcptr, section_biglitptr, section_ttpptr, loopcount2, loopcount3 ;; Rotate carries and add in carries
	mov	rsi, saved_reg1		;; Restore FFT data addr
ttp	mov	rdi, saved_reg2		;; Restore big/lit pointer
ttp	mov	rbp, saved_reg3		;; Restore ttp pointer
echk	vmovapd	zmm27, ZMM_MAXERR	;; Restore maximum error BUG/OPT - leave maxerr in zmm29 like zpad code does

	sub	loopcount1, 1		;; Test section counter
	jnz	ilp0
;;echk	zstore	ZMM_MAXERR, zmm27	;; Save maximum error		BUG/OPT - leave maxerr in zmm27 for common code
									;; Better yet, ditch common code and do sumout / maxerr here
zero ttp	jmp	zdn		;; Go to zero upper half irrational end code
zero no ttp	jmp	zrdn		;; Go to zero upper half rational end code
no zero ttp	jmp	dn		;; Go to irrational end code
no zero no ttp	jmp	rdn		;; Go to rational end code
	ENDPP lab
	ENDM

zpnorm	MACRO	lab, ttp, echk, const, khi
	LOCAL	ilp0, ilp1, ilp2, ilpsecdn
	PROCFLP	lab
	int_prolog 0,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
	mov	r12, ZMM_QUARTER_DIST	;; Distance to FFTLEN/4 FFT word -- independent set of FFT data to work on
	vbroadcastsd zmm0, ZMM_RNDVAL	;; Start process with no carry
	vmovapd	zmm1, zmm0
	vmovapd	zmm2, zmm0
	vmovapd	zmm3, zmm0
echk	vpxorq	zmm29, zmm29, zmm29	;; Clear maximum error
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	mov	loopcount1, ZMM_NORM1_SECTIONS ;; Load loop counter (4 sections)
ilp0:	znorm_1d_zpad_preload ttp, echk, const, khi ;; Preload useful constants
	vpxorq	zmm30, zmm30, zmm30	;; Clear auxilary sumout register (zmm31 was initialized in zsub_7_words)
echk	vmovapd	zmm28, zmm29		;; Init (or clear) auxilary maximum error register
	mov	section_srcptr, rsi	;; remember rsi for znorm_1d_mid_cleanup
ttp	mov	section_biglitptr, rdi	;; remember rdi for znorm_1d_mid_cleanup
ttp	mov	section_ttpptr, rbp	;; remember rbp for znorm_1d_mid_cleanup
	mov	loopcount2, ZMM_NORM1_PAD_GRPS ;; Number of subsections (padding groups)
ilp2:	mov	loopcount3, ZMM_NORM1_CL_TIL_PAD ;; Count of cache lines before a padding occurs or subsection ends
ilp1:	znorm_1d_zpad ttp, echk, const, khi ;; Normalize 32 values (4 cache lines)
	bump	rsi, 2*64		;; Next set of 4 cache lines (2 rsi cache lines, 2 rsi+r12 cache lines)
ttp	bump	rbp, 2*128		;; Next set of multipliers
ttp	bump	rdi, 2*1		;; Next big/little flags
	sub	loopcount3, 4		;; Test cache line loop counter
	jnz	ilp1			;; Loop til done
	sub	loopcount2, 1		;; Test count of pad groups (may be zero!)
	js	short ilpsecdn		;; Loop til done
	lea	rsi, [rsi+64]		;; Pad 64 bytes, preserve flags
	jnz	ilp2			;; Loop til done
ilpsecdn:
	vaddpd	zmm31, zmm31, zmm30	;; Coalesce sumouts
echk	vmaxpd	zmm29, zmm29, zmm28	;; Coalesce maximum error
	mov	saved_reg1, rsi		;; Save FFT data addr
ttp	mov	saved_reg2, rdi		;; Save big/lit pointer
ttp	mov	saved_reg3, rbp		;; Save ttp pointer
	znorm_1d_zpad_mid_cleanup_preload ttp, const, khi
	znorm_1d_zpad_mid_cleanup ttp, const, khi, section_srcptr, section_biglitptr, section_ttpptr, loopcount2, loopcount3 ;; Rotate carries and add in carries
	mov	rsi, saved_reg1		;; Restore FFT data addr
ttp	mov	rdi, saved_reg2		;; Restore big/lit pointer
ttp	mov	rbp, saved_reg3		;; Restore ttp pointer
	sub	loopcount1, 1		;; Test section counter
	jnz	ilp0
	znorm_1d_zpad_final_cleanup_preload ttp, const, khi
	znorm_1d_zpad_final_cleanup ttp, const, khi ;; Apply final carries that do not wraparound
echk	zstore	ZMM_MAXERR, zmm29	;; Save maximum error		BUG/OPT - leave maxerr in zmm29 for common code
									;; Better yet, ditch common code and do sumout / maxerr here
ttp no const jmp zpdn			;; Go to zero padded FFT end code
ttp const jmp	zpcdn			;; Go to zero padded FFT end code
no ttp no const jmp zprdn		;; Go to zero padded FFT end code
no ttp const jmp zprcdn			;; Go to zero padded FFT end code
	ENDPP lab
	ENDM

; The many different normalization routines.  One for each valid combination of
; rational/irrational, zeroing/no zeroing, error check/no error check,
; mul by const/no mul by const

	inorm	zr1, noexec, noexec, noexec, noexec
	inorm	zr1e, noexec, noexec, exec, noexec
	inorm	zr1c, noexec, noexec, noexec, exec
	inorm	zr1ec, noexec, noexec, exec, exec
	inorm	zr1z, noexec, exec, noexec, noexec
	inorm	zr1ze, noexec, exec, exec, noexec
	inorm	zi1, exec, noexec, noexec, noexec
	inorm	zi1e, exec, noexec, exec, noexec
	inorm	zi1c, exec, noexec, noexec, exec
	inorm	zi1ec, exec, noexec, exec, exec
	inorm	zi1z, exec, exec, noexec, noexec
	inorm	zi1ze, exec, exec, exec, noexec
	zpnorm	zr1zp, noexec, noexec, noexec, exec
	zpnorm	zr1zpe, noexec, exec, noexec, exec
	zpnorm	zr1zpc, noexec, noexec, exec, exec
	zpnorm	zr1zpec, noexec, exec, exec, exec
	zpnorm	zi1zp, exec, noexec, noexec, exec
	zpnorm	zi1zpe, exec, exec, noexec, exec
	zpnorm	zi1zpc, exec, noexec, exec, exec
	zpnorm	zi1zpec, exec, exec, exec, exec
	zpnorm	zr1zpk, noexec, noexec, noexec, noexec
	zpnorm	zr1zpek, noexec, exec, noexec, noexec
	zpnorm	zr1zpck, noexec, noexec, exec, noexec
	zpnorm	zr1zpeck, noexec, exec, exec, noexec
	zpnorm	zi1zpk, exec, noexec, noexec, noexec
	zpnorm	zi1zpek, exec, exec, noexec, noexec
	zpnorm	zi1zpck, exec, noexec, exec, noexec
	zpnorm	zi1zpeck, exec, exec, exec, noexec

; Common code to finish off the one-pass FFTs normalization.  The
; Windows 64-bit ABI frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.

PROCFP	__common_znorm1_end_code

	;; Dummy prolog to match normalization code
	int_prolog 0,0,0

; Finish off the normalization process by adding any carry to first values.
; Handle both the with and without two-to-phi array cases.

dn:	mov	rsi, DESTARG		; Address of squared number
	mov	rbp, norm_col_mults	; Address of the two-to-phi multipliers
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	znorm_1d_cleanup_preload exec	; Preloads for znorm_top_carry_1d too!
	znorm_top_carry_1d		; Adjust top carry when k > 1
	znorm_1d_cleanup exec, noexec	; Add in carries
	jmp	cmnend			; All done, go cleanup

rdn:	mov	rsi, DESTARG		; Address of squared number
	znorm_1d_cleanup_preload noexec
	znorm_1d_cleanup noexec, noexec ; Add in carries
	jmp	cmnend			; All done, go cleanup

zdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rbp, norm_col_mults	; Address of the two-to-phi multipliers
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	znorm_1d_cleanup_preload exec
	znorm_1d_cleanup exec, exec	; Add in carries
	jmp	cmnend			; All done, go cleanup

zrdn:	mov	rsi, DESTARG		; Address of squared number
	znorm_1d_cleanup_preload noexec
	znorm_1d_cleanup noexec, exec	; Add in carries
	jmp	cmnend			; All done, go cleanup

zpdn:	znorm_1d_zpad_cleanup_preload exec, noexec
	znorm_1d_zpad_cleanup exec, noexec	; Add in carries
	jmp	cmnend				; All done, go cleanup

zpcdn:	znorm_1d_zpad_cleanup_preload exec, exec
	znorm_1d_zpad_cleanup exec, exec	; Add in carries
	jmp	cmnend				; All done, go cleanup

zprdn:	znorm_1d_zpad_cleanup_preload noexec, noexec
	znorm_1d_zpad_cleanup noexec, noexec	; Add in carries
	jmp	cmnend				; All done, go cleanup

zprcdn:	znorm_1d_zpad_cleanup_preload noexec, exec
	znorm_1d_zpad_cleanup noexec, exec	; Add in carries
;	jmp	cmnend				; All done, go cleanup

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

cmnend:	mov	rsi, DESTARG			; Address of squared number
	vshuff64x2 zmm1, zmm31, zmm31, 00001011b ; Move top 256-bits to bottom
	vaddpd	zmm0, zmm31, zmm1		; We now have just 4 sumout values
	vshuff64x2 zmm1, zmm0, zmm0, 00000001b	; Move top 128-bits (of the 256-bits) to bottom
	vaddpd	zmm0, zmm0, zmm1		; We now have just 2 sumout values
	vshufpd	zmm1, zmm0, zmm0, 1		; Move top 64-bits (of the 128-bits) to bottom
	vaddsd	xmm0, xmm0, xmm1
	vmulsd	xmm0, xmm0, ttmp_ff_inv
	vmovsd	Q [rsi-24], xmm0		; Save sum of FFT outputs

; Collapse the maximum error into one double

	vmovapd	zmm29, ZMM_MAXERR
	vshuff64x2 zmm1, zmm29, zmm29, 00001011b ; Move top 256-bits to bottom
	vmaxpd	zmm0, zmm29, zmm1		; We now have just 4 maxerr values
	vshuff64x2 zmm1, zmm0, zmm0, 00000001b	; Move top 128-bits (of the 256-bits) to bottom
	vmaxpd	zmm0, zmm0, zmm1		; We now have just 2 maxerr values
	vshufpd	zmm1, zmm0, zmm0, 1		; Move top 64-bits (of the 128-bits) to bottom
	vmaxsd	xmm0, xmm0, xmm1
	vmaxsd	xmm0, xmm0, MAXERR		; Compute new maximum error
	vmovsd	MAXERR, xmm0

; Return

	int_epilog 0,0,0
	ENDPP __common_znorm1_end_code 

_TEXT	ENDS
END
