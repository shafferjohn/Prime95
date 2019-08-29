
	TITLE   cpuidhlp

OPTION	EVEX:1

_TEXT SEGMENT

Q	EQU	<QWORD PTR>
SZPTR	EQU	8
YPTR	EQU	<YMMWORD PTR>
ZPTR	EQU	<ZMMWORD PTR>
AD_BASE		EQU	<r11>
UNION_BASE		EQU	AD_BASE+48*SZPTR+576
ZMM_MAXERR		EQU	ZPTR [UNION_BASE+16*SZPTR+56*8]
YMM_TMP6		EQU	YPTR [UNION_BASE+16*SZPTR+564*8]


zfmaddpd MACRO dest, mulval1, mulval2, addval
	LOCAL	destreglen
	destreglen = @SIZESTR(<&dest>)
	IF @INSTR(,@substr(<&dest>, 1, destreglen),<{>) NE 0
		destreglen = @INSTR(,@substr(<&dest>, 1, destreglen),<{>) - 1
	ENDIF
	IF @INSTR(,@substr(<&dest>, 1, destreglen),< >) NE 0
		destreglen = @INSTR(,@substr(<&dest>, 1, destreglen),< >) - 1
	ENDIF
	IFIDNI @substr(<&dest>, 1, destreglen), <&addval>
		IF (OPATTR (mulval1)) AND 10000b
			vfmadd231pd dest, mulval1, mulval2
		ELSE
			vfmadd231pd dest, mulval2, mulval1
		ENDIF
	ELSE
	IFIDNI @substr(<&dest>, 1, destreglen), <&mulval1>
		IF (OPATTR (mulval2)) AND 10000b
			vfmadd213pd dest, mulval2, addval
		ELSE
			vfmadd132pd dest, addval, mulval2
		ENDIF
	ELSE
	IFIDNI @substr(<&dest>, 1, destreglen), <&mulval2>
		IF (OPATTR (mulval1)) AND 10000b
			vfmadd213pd dest, mulval1, addval
		ELSE
			vfmadd132pd dest, addval, mulval1
		ENDIF
	ELSE
	IF (OPATTR (mulval2)) AND 10000b
		vmovapd dest, mulval2
		IF (OPATTR (mulval1)) AND 10000b
			vfmadd213pd dest, mulval1, addval
		ELSE
			vfmadd132pd dest, addval, mulval1
		ENDIF
	ELSE
		vmovapd dest, mulval1
		vfmadd132pd dest, addval, mulval2
	ENDIF
	ENDIF
	ENDIF
	ENDIF
	ENDM


fpu_init	proc
;;	mov	rax, 1

	vmovsd	xmm7, [r11+88]


	zfmaddpd zmm0, zmm1, zmm2, zmm3
	zfmaddpd zmm0{k1}, zmm1, zmm2, zmm3
	zfmaddpd zmm0 {k1}, zmm1, zmm2, zmm3
	zfmaddpd zmm0, zmm1, zmm2, zmm0
	zfmaddpd zmm0 , zmm1, zmm2, zmm0
	zfmaddpd zmm0{k1}, zmm1, zmm2, zmm0
	zfmaddpd zmm0 {k1}, zmm1, zmm2, zmm0
	zfmaddpd zmm0{k1} , zmm1, zmm2, zmm0

vmovsd xmm0, [r11+16*8+48]
mov ecx, 01010101b

vmovsd xmm10, xmm12, xmm12

	vsubsd	xmm3, xmm3, xmm30
	vmovsd	xmm3, [r11+320]
	vmovsd	xmm3, Q YMM_TMP6

	vmovsd	xmm0, QWORD PTR [r11+8]
	vsubsd	xmm31, xmm31, xmm0
	vaddsd	xmm0, xmm0, QWORD PTR [rsi][rdi]
	vmovsd	QWORD PTR [rsi][rdi], xmm0
mov	BYTE PTR YMMWORD PTR [r11+56], al
;vmaxsd	xmm6, xmm6, Q ZMM_MAXERR+8

	vmovdqa64 zmm31, [r11+32]

vmaxsd xmm6, xmm6, [rdi]+8

;vscatterdpd [srcreg+ymm1+8], zmm14
;vgatherdpd zmm0, [rsi+ymm1*1+0], zmm2

vgatherqpd zmm0{k1}, [rsi+zmm31*1+0]
vscatterqpd [rsi+zmm31*1+0] {k1}, zmm0
;	kmovw	k1, 01010101b
	kmovw	k1, eax
	kmovw	k1, ebx
	kmovw	k1, ecx
	kmovw	k1, edx
	kmovw	k1, esi
	kmovw	k1, edi
	kmovw	k1, ebp
;	kmovw	k1, r8w
	kmovw	k1, r8d
;	kmovw	k1, r8
;	kmovw	k1, r15w
	kmovw	k1, r15d
;	kmovw	k1, r15
	vmovsd  xmm6, QWORD PTR [r11+16*8+16]
	vblendmpd zmm16 {k1}, zmm22, zmm23		;; Create (RNDVAL * base - RNDVAL) constant used new value calculations
	vblendmpd zmm16 {k1}, zmm22, zmm31		;; Create (RNDVAL * base - RNDVAL) constant used new value calculations
	vblendmpd zmm16 {k1}, zmm22, zmm16		;; Create (RNDVAL * base - RNDVAL) constant used new value calculations
	vblendmpd zmm0{k1}, zmm2, zmm1
	vblendmpd zmm0{k1}, zmm2, zmm1
	ret
fpu_init ENDP

_TEXT	ENDS
END
