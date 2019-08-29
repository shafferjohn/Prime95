; Copyright 1995-2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Time low-level operations for optimizing macros and comparing CPUs
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE memory.mac
INCLUDE pfa.mac
INCLUDE lucas.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE hg.mac
INCLUDE r4.mac

IFDEF X86_64
X87_CASES	EQU	0
ELSE
X87_CASES	EQU	13
ENDIF
SSE2_CASES	EQU	216
AVX_CASES	EQU	300
IFDEF X86_64
AVX512_CASES	EQU	140
ELSE
AVX512_CASES	EQU	0
ENDIF
IACA_TEST_CASE	EQU	999

loopent	MACRO	y,z		; Create a entry in the loop entry table
	DP	&y&z
	ENDM
looptab	MACRO	y, cnt		; Create the loop entry table
	x = 0
	REPT	cnt
	loopent	y, %x
	x = x + 1
	ENDM
	ENDM

;; Macros from p4notes.doc

read1	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
	cnt = mem/64		        ;; Read 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:	mov     eax, [rsi]		;; Read one cache line
	mov     eax, [rsi+4]		;; 4 bytes at a time
	mov     eax, [rsi+8]
	mov     eax, [rsi+12]
	mov     eax, [rsi+16]
	mov     eax, [rsi+20]
	mov     eax, [rsi+24]
	mov     eax, [rsi+28]
	mov     eax, [rsi+32]
	mov     eax, [rsi+36]
	mov     eax, [rsi+40]
	mov     eax, [rsi+44]
	mov     eax, [rsi+48]
	mov     eax, [rsi+52]
	mov     eax, [rsi+56]
	mov     eax, [rsi+60]
	lea     rsi, [rsi+64]			; Next cache line
	sub     ecx, 1
	jnz     loop1
	lea     rsi, [rsi-mem]		; Restore esi
	dec	edx
	jnz	loop2
	ENDM

read2	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm1, [rsi]			; Read one cache line
        movapd  xmm1, [rsi+16]		; 16 bytes at a time
        movapd  xmm1, [rsi+32]
        movapd  xmm1, [rsi+48]
        lea     rsi, [rsi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

read4	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
        cnt = mem/(4*64)                ;; 4 cache lines per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  vmovapd  ymm0, [rsi]		;; Read cache lines
        vmovapd  ymm1, [rsi+32]
        vmovapd  ymm2, [rsi+64]
        vmovapd  ymm3, [rsi+96]
        vmovapd  ymm4, [rsi+128]
        vmovapd  ymm5, [rsi+160]
        vmovapd  ymm6, [rsi+192]
        vmovapd  ymm7, [rsi+224]
        add     rsi, 4*64		;; Next cache lines
        sub     ecx, 1
        jnz     loop1
        sub     rsi, cnt*4*64
	dec	edx
	jnz	loop2
	ENDM

read8	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
        cnt = mem/(4*64)                ;; 4 cache lines per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm1, [rsi+64]
        vmovapd  zmm2, [rsi+128]
        vmovapd  zmm3, [rsi+192]
        add     rsi, 4*64		;; Next cache lines
        sub     ecx, 1
        jnz     loop1
        sub     rsi, cnt*4*64
	dec	edx
	jnz	loop2
	ENDM

write1	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  [rsi], xmm1			; Write one cache line
        movapd  [rsi+16], xmm1		; 16 bytes at a time
        movapd  [rsi+32], xmm1
        movapd  [rsi+48], xmm1
        lea     rsi, [rsi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

write2	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/(4*128)               ; 128 bytes per iteration
        dist = 64
	mov	edx, c
	sub     ebx, ebx
loop2:	mov     ecx, cnt
loop1:  movapd  [rsi+0*dist], xmm1      ; Write 8 cache lines
        movapd  [rsi+1*dist], xmm1
        movapd  [rsi+2*dist], xmm1
        movapd  [rsi+3*dist], xmm1
        movapd  [rsi+4*dist], xmm1
        movapd  [rsi+5*dist], xmm1
        movapd  [rsi+6*dist], xmm1
        movapd  [rsi+7*dist], xmm1
        lea     rsi, [rsi+16]           ; Same cache lines
        add     bl, 256/4			; 4 inner loop iterations
        jnc     loop1
        lea     rsi, [rsi-4*16+8*dist]  ; Next set of 8 cache lines
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

write4	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/(4*64)                ;; 4 cache lines per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  vmovapd  [rsi], ymm0		;; Write cache lines
        vmovapd  [rsi+32], ymm1
        vmovapd  [rsi+64], ymm2
        vmovapd  [rsi+96], ymm3
        vmovapd  [rsi+128], ymm4
        vmovapd  [rsi+160], ymm5
        vmovapd  [rsi+192], ymm6
        vmovapd  [rsi+224], ymm7
        add     rsi, 4*64		;; Next cache lines
        sub     ecx, 1
        jnz     loop1
        sub     rsi, cnt*4*64
	dec	edx
	jnz	loop2
	ENDM

write8	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/(4*64)                ;; 4 cache lines per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  vmovapd  [rsi], zmm0		;; Write cache lines
        vmovapd  [rsi+64], zmm1
        vmovapd  [rsi+128], zmm2
        vmovapd  [rsi+192], zmm3
        add     rsi, 4*64		;; Next cache lines
        sub     ecx, 1
        jnz     loop1
        sub     rsi, cnt*4*64
	dec	edx
	jnz	loop2
	ENDM

readwrite1 MACRO mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm0, [rsi]             ; Read one cache line
        movapd  xmm1, [rsi+16]
        movapd  xmm2, [rsi+32]
        movapd  xmm3, [rsi+48]
        subpd   xmm0, xmm0              ; Operate on the data
        pxor    xmm1, xmm1
        subpd   xmm2, xmm2
        pxor    xmm3, xmm3
        movapd  [rsi], xmm0             ; Write the cache line
        movapd  [rsi+16], xmm1
        movapd  [rsi+32], xmm2
        movapd  [rsi+48], xmm3
        lea     rsi, [rsi+64]           ; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

readwrite4 MACRO mem, c			;; Bytes to read/write
	LOCAL loop1, loop2
        cnt = mem/(4*64)                ;; 4 cache lines per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  vmovapd  ymm0, [rsi]		;; Read cache lines
        vmovapd  ymm1, [rsi+32]
        vmovapd  ymm2, [rsi+64]
        vmovapd  ymm3, [rsi+96]
        vmovapd  ymm4, [rsi+128]
        vmovapd  ymm5, [rsi+160]
        vmovapd  ymm6, [rsi+192]
        vmovapd  ymm7, [rsi+224]
	vmovapd  [rsi], ymm0		;; Write cache lines
        vmovapd  [rsi+32], ymm1
        vmovapd  [rsi+64], ymm2
        vmovapd  [rsi+96], ymm3
        vmovapd  [rsi+128], ymm4
        vmovapd  [rsi+160], ymm5
        vmovapd  [rsi+192], ymm6
        vmovapd  [rsi+224], ymm7
        add     rsi, 4*64		;; Next cache lines
        sub     ecx, 1
        jnz     loop1
        sub     rsi, cnt*4*64
	dec	edx
	jnz	loop2
	ENDM

readwrite8 MACRO mem, c			;; Bytes to read/write
	LOCAL loop1, loop2
        cnt = mem/(4*64)                ;; 4 cache lines per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm2, [rsi+64]
        vmovapd  zmm4, [rsi+128]
        vmovapd  zmm6, [rsi+192]
	vmovapd  [rsi], zmm0		;; Write cache lines
        vmovapd  [rsi+64], zmm2
        vmovapd  [rsi+128], zmm4
        vmovapd  [rsi+192], zmm6
        add     rsi, 4*64		;; Next cache lines
        sub     ecx, 1
        jnz     loop1
        sub     rsi, cnt*4*64		;; Restore source pointer
	dec	edx
	jnz	loop2
	ENDM

clocks50 MACRO
	REPEAT 50
	vpandq zmm0, zmm0, zmm1
	ENDM
	ENDM

sequential_readwrite8 MACRO mem, c		;; Bytes to read/write
	LOCAL loop1, loop2, loop3
        cnt = mem/(4*64)                ;; 4 cache lines per iteration -- equivalent to clm=2
	cnt1 = 1280 / 2			;; Assumes pass 1 size of 1280, which equals 640 cache line pairs
        cnt2 = cnt / cnt1		;; This is number of pass 1 blocks
	dist = 7680*16 + 29*64		;; Assumes pass 2 does 7680 * 16 bytes + 29*64 pad bytes
	mov	edx, c
loop3:	mov     ecx, cnt1
loop2:	mov     ebx, cnt2
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm2, [rsi+64]
        vmovapd  zmm4, [rsi+128]
        vmovapd  zmm6, [rsi+192]
prefetcht0 [rsi+4*64+0]
prefetcht0 [rsi+4*64+64]
;	clocks50
prefetcht0 [rsi+4*64+128]
prefetcht0 [rsi+4*64+192]
	clocks50
	vmovapd  [rsi], zmm0		;; Write cache lines
        vmovapd  [rsi+64], zmm2
        vmovapd  [rsi+128], zmm4
        vmovapd  [rsi+192], zmm6
;	vmovntpd  [rsi], zmm0		;; Write cache lines
;        vmovntpd  [rsi+64], zmm2
;        vmovntpd  [rsi+128], zmm4
;        vmovntpd  [rsi+192], zmm6
        add     rsi, 4*64		;; Next cache lines
        sub     ebx, 1
        jnz     loop1
        sub     ecx, 1
        jnz     loop2
        sub     rsi, cnt2*cnt1*4*64	;; Restore source pointer
	dec	edx
	jnz	loop3
	ENDM

scatter_readwrite8 MACRO mem, c		;; Bytes to read/write
	LOCAL loop1, loop2, loop3
        cnt = mem/(4*64)                ;; 4 cache lines per iteration -- equivalent to clm=2
	cnt1 = 1280 / 2			;; Assumes pass 1 size of 1280, which equals 640 cache line pairs
        cnt2 = cnt / cnt1		;; This is number of pass 1 blocks
	dist = 7680*16 + 29*64		;; Assumes pass 2 does 7680 * 16 bytes + 29*64 pad bytes
	mov	edx, c
loop3:	mov     ecx, cnt2
loop2:	mov     ebx, cnt1
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm2, [rsi+64]
        vmovapd  zmm4, [rsi+128]
        vmovapd  zmm6, [rsi+192]
prefetcht0 [rsi+dist+0]
prefetcht0 [rsi+dist+64]
;	clocks50
prefetcht0 [rsi+dist+128]
prefetcht0 [rsi+dist+192]
	clocks50
	vmovapd  [rsi], zmm0		;; Write cache lines
        vmovapd  [rsi+64], zmm2
        vmovapd  [rsi+128], zmm4
        vmovapd  [rsi+192], zmm6
;	vmovntpd  [rsi], zmm0		;; Write cache lines
;        vmovntpd  [rsi+64], zmm2
;        vmovntpd  [rsi+128], zmm4
;        vmovntpd  [rsi+192], zmm6
        add     rsi, dist		;; Next cache lines
        sub     ebx, 1
        jnz     loop1
        add     rsi, -cnt1*dist+4*64	;; Next cache lines
        sub     ecx, 1
        jnz     loop2
        sub     rsi, cnt2*4*64		;; Restore source pointer
	dec	edx
	jnz	loop3
	ENDM

;;newer version simulating half as many large strides
scatter_readwrite8 MACRO mem, c		;; Bytes to read/write
	LOCAL loop1, loop2, loop3
        cnt = mem/(8*64)                ;; 4 cache lines per iteration -- equivalent to clm=4
	cnt1 = 1280 / 2			;; Assumes pass 1 size of 1280, which equals 640 cache line pairs
        cnt2 = cnt / cnt1		;; This is number of pass 1 blocks
	dist = 7680*16 + 29*64		;; Assumes pass 2 does 7680 * 16 bytes + 29*64 pad bytes
	mov	edx, c
loop3:	mov     ecx, cnt2
loop2:	mov     ebx, cnt1
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm2, [rsi+64]
        vmovapd  zmm4, [rsi+128]
        vmovapd  zmm6, [rsi+192]
prefetcht0 [rsi+dist+0]
prefetcht0 [rsi+dist+64]
prefetcht0 [rsi+dist+128]
prefetcht0 [rsi+dist+192]
	clocks50
	vmovapd  [rsi], zmm0		;; Write cache lines
        vmovapd  [rsi+64], zmm2
        vmovapd  [rsi+128], zmm4
        vmovapd  [rsi+192], zmm6

	vmovapd  zmm0, [rsi+256]		;; Read cache lines
        vmovapd  zmm2, [rsi+256+64]
        vmovapd  zmm4, [rsi+256+128]
        vmovapd  zmm6, [rsi+256+192]
prefetcht0 [rsi+dist+256+0]
prefetcht0 [rsi+dist+256+64]
prefetcht0 [rsi+dist+256+128]
prefetcht0 [rsi+dist+256+192]
	clocks50
	vmovapd  [rsi+256], zmm0		;; Write cache lines
        vmovapd  [rsi+256+64], zmm2
        vmovapd  [rsi+256+128], zmm4
        vmovapd  [rsi+256+192], zmm6

        add     rsi, dist		;; Next cache lines
        sub     ebx, 1
        jnz     loop1
        add     rsi, -cnt1*dist+8*64	;; Next cache lines
        sub     ecx, 1
        jnz     loop2
        sub     rsi, cnt2*8*64		;; Restore source pointer
	dec	edx
	jnz	loop3
	ENDM

sequential_read_scatter_write8 MACRO mem, c ;; Bytes to read/write
	LOCAL loop1, loop2, loop3
        cnt = mem/(4*64)                ;; 4 cache lines per iteration -- equivalent to clm=2
	cnt1 = 1280 / 2			;; Assumes pass 1 size of 1280, which equals 640 cache line pairs
        cnt2 = cnt / cnt1		;; This is number of pass 1 blocks
	dist = 7680*16 + 29*64		;; Assumes pass 2 does 7680 * 16 bytes + 29*64 pad bytes
	lea	rdi, [rsi+96*1024*1024]	;; write to separate mem address 96MB away
	mov	edx, c
loop3:	mov     ecx, cnt2
loop2:	mov     ebx, cnt1
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm2, [rsi+64]
        vmovapd  zmm4, [rsi+128]
        vmovapd  zmm6, [rsi+192]
prefetcht0 [rsi+4*64+0]
prefetcht0 [rsi+4*64+64]
;	clocks50
prefetcht0 [rsi+4*64+128]
prefetcht0 [rsi+4*64+192]
	clocks50
;	vmovapd  [rdi], zmm0		;; Write cache lines
;        vmovapd  [rdi+64], zmm2
;        vmovapd  [rdi+128], zmm4
;        vmovapd  [rdi+192], zmm6
	vmovntpd  [rdi], zmm0		;; Write cache lines
        vmovntpd  [rdi+64], zmm2
        vmovntpd  [rdi+128], zmm4
        vmovntpd  [rdi+192], zmm6
        add     rsi, 4*64		;; Next cache lines
        add     rdi, dist		;; Next cache lines
        sub     ebx, 1
        jnz     loop1
        add     rdi, -cnt1*dist+4*64	;; Next cache lines
        sub     ecx, 1
        jnz     loop2
        sub     rsi, cnt1*cnt2*4*64	;; Restore source pointer
        sub     rdi, cnt2*4*64		;; Restore dest pointer
	dec	edx
	jnz	loop3
	ENDM

scatter_read_sequential_write8 MACRO mem, c ;; Bytes to read/write
	LOCAL loop1, loop2, loop3
        cnt = mem/(4*64)                ;; 4 cache lines per iteration -- equivalent to clm=2
	cnt1 = 1280 / 2			;; Assumes pass 1 size of 1280, which equals 640 cache line pairs
        cnt2 = cnt / cnt1		;; This is number of pass 1 blocks
	dist = 7680*16 + 29*64		;; Assumes pass 2 does 7680 * 16 bytes + 29*64 pad bytes
	lea	rdi, [rsi+96*1024*1024]	;; write to separate mem address 96MB away
	mov	edx, c
loop3:	mov     ecx, cnt2
loop2:	mov     ebx, cnt1
loop1:  vmovapd  zmm0, [rsi]		;; Read cache lines
        vmovapd  zmm2, [rsi+64]
        vmovapd  zmm4, [rsi+128]
        vmovapd  zmm6, [rsi+192]
prefetcht0 [rsi+dist+0]
prefetcht0 [rsi+dist+64]
;	clocks50
prefetcht0 [rsi+dist+128]
prefetcht0 [rsi+dist+192]
	clocks50
;	vmovapd  [rdi], zmm0		;; Write cache lines
;        vmovapd  [rdi+64], zmm2
;        vmovapd  [rdi+128], zmm4
;        vmovapd  [rdi+192], zmm6
	vmovntpd  [rdi], zmm0		;; Write cache lines
        vmovntpd  [rdi+64], zmm2
        vmovntpd  [rdi+128], zmm4
        vmovntpd  [rdi+192], zmm6
        add     rsi, dist		;; Next cache lines
        add     rdi, 4*64		;; Next cache lines
        sub     ebx, 1
        jnz     loop1
        add     rsi, -cnt1*dist+4*64	;; Next cache lines
        sub     ecx, 1
        jnz     loop2
        sub     rsi, cnt2*4*64		;; Restore source pointer
        sub     rdi, cnt1*cnt2*4*64	;; Restore dest pointer
	dec	edx
	jnz	loop3
	ENDM

x4cl_empty MACRO srcreg,srcinc,d1,d2,screg,scoff
	xload	xmm0, [srcreg+0]
	xload	xmm1, [srcreg+32]
	xload	xmm2, [srcreg+d1+0]
	xload	xmm3, [srcreg+d1+32]
	xload	xmm4, [srcreg+d2+0]
	xload	xmm5, [srcreg+d2+32]
	xload	xmm6, [srcreg+d2+d1+0]
	xload	xmm7, [srcreg+d2+d1+32]
	xstore	[srcreg+0], xmm0
	xstore	[srcreg+32], xmm1
	xload	xmm0, [srcreg+0+16]
	xload	xmm1, [srcreg+32+16]
	xstore	[srcreg+16], xmm2
	xstore	[srcreg+48], xmm3
	xload	xmm2, [srcreg+d1+0+16]
	xload	xmm3, [srcreg+d1+32+16]
	xstore	[srcreg+d1], xmm4
	xstore	[srcreg+d1+16], xmm5
	xstore	[srcreg+d1+32], xmm6
	xstore	[srcreg+d1+48], xmm7
	xload	xmm4, [srcreg+d2+0+16]
	xload	xmm5, [srcreg+d2+32+16]
	xload	xmm6, [srcreg+d2+d1+0+16]
	xload	xmm7, [srcreg+d2+d1+32+16]
	xstore	[srcreg+d2], xmm0
	xstore	[srcreg+d2+16], xmm1
	xstore	[srcreg+d2+32], xmm2
	xstore	[srcreg+d2+48], xmm3
	xstore	[srcreg+d2+d1], xmm4
	xstore	[srcreg+d2+d1+16], xmm5
	xstore	[srcreg+d2+d1+32], xmm6
	xstore	[srcreg+d2+d1+48], xmm7
	bump	srcreg, srcinc
	ENDM

g4cl_empty MACRO srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg,scoff
	xload	xmm0, [srcreg+0]
	xload	xmm1, [srcreg+32]
	xload	xmm2, [srcreg+d1+0]
	xload	xmm3, [srcreg+d1+32]
	xload	xmm4, [srcreg+d2+0]
	xload	xmm5, [srcreg+d2+32]
	xload	xmm6, [srcreg+d2+d1+0]
	xload	xmm7, [srcreg+d2+d1+32]
	xstore	[dstreg+0], xmm0
	xstore	[dstreg+32], xmm1
	xload	xmm0, [srcreg+0+16]
	xload	xmm1, [srcreg+32+16]
	xstore	[dstreg+16], xmm2
	xstore	[dstreg+48], xmm3
	xload	xmm2, [srcreg+d1+0+16]
	xload	xmm3, [srcreg+d1+32+16]
	xstore	[dstreg+e1], xmm4
	xstore	[dstreg+e1+16], xmm5
	xstore	[dstreg+e1+32], xmm6
	xstore	[dstreg+e1+48], xmm7
	xload	xmm4, [srcreg+d2+0+16]
	xload	xmm5, [srcreg+d2+32+16]
	xload	xmm6, [srcreg+d2+d1+0+16]
	xload	xmm7, [srcreg+d2+d1+32+16]
	bump	srcreg, srcinc
	xstore	[dstreg+e2], xmm0
	xstore	[dstreg+e2+16], xmm1
	xstore	[dstreg+e2+32], xmm2
	xstore	[dstreg+e2+48], xmm3
	xstore	[dstreg+e2+e1], xmm4
	xstore	[dstreg+e2+e1+16], xmm5
	xstore	[dstreg+e2+e1+32], xmm6
	xstore	[dstreg+e2+e1+48], xmm7
	bump	dstreg, dstinc
	ENDM

g4cl_empty_nt MACRO srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg,scoff
	xload	xmm0, [srcreg+0]
	xload	xmm1, [srcreg+32]
	xload	xmm2, [srcreg+d1+0]
	xload	xmm3, [srcreg+d1+32]
	xload	xmm4, [srcreg+d2+0]
	xload	xmm5, [srcreg+d2+32]
	xload	xmm6, [srcreg+d2+d1+0]
	xload	xmm7, [srcreg+d2+d1+32]
	movntpd	[dstreg+0], xmm0
	movntpd	[dstreg+32], xmm1
	xload	xmm0, [srcreg+0+16]
	xload	xmm1, [srcreg+32+16]
	movntpd	[dstreg+16], xmm2
	movntpd	[dstreg+48], xmm3
	xload	xmm2, [srcreg+d1+0+16]
	xload	xmm3, [srcreg+d1+32+16]
	movntpd	[dstreg+e1], xmm4
	movntpd	[dstreg+e1+16], xmm5
	movntpd	[dstreg+e1+32], xmm6
	movntpd	[dstreg+e1+48], xmm7
	xload	xmm4, [srcreg+d2+0+16]
	xload	xmm5, [srcreg+d2+32+16]
	xload	xmm6, [srcreg+d2+d1+0+16]
	xload	xmm7, [srcreg+d2+d1+32+16]
	bump	srcreg, srcinc
	movntpd	[dstreg+e2], xmm0
	movntpd	[dstreg+e2+16], xmm1
	movntpd	[dstreg+e2+32], xmm2
	movntpd	[dstreg+e2+48], xmm3
	movntpd	[dstreg+e2+e1], xmm4
	movntpd	[dstreg+e2+e1+16], xmm5
	movntpd	[dstreg+e2+e1+32], xmm6
	movntpd	[dstreg+e2+e1+48], xmm7
	bump	dstreg, dstinc
	ENDM


;; Time one of the basic FFT building blocks

x87mac	MACRO	memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	mov	eax, outer_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr
	lea	rsi, [rdi+4096]
	mov	ecx, inner_iters
ss0b:	disp	&ops
	lea	rsi, [rsi+memused]	;; Next source pointer
	lea	rdi, [rdi+SCD]		;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	ENDM

sse2mac MACRO	lab, memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 0			;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+262144+4096+64] ;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
ss0b:	&ops
IF memused NE 192
	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
ELSE
	lea	rdi, [rdi+XMM_SCD1]	;; Next sine/cosine pointer
ENDIF
	dec	ecx
	jnz	ss0b
	mov	ecx, inner_iters
	dec	eax
	jnz	ss0a
	jmp	exit
	ENDM
sse2macbx MACRO	lab, memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 262144+128		;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+262144+4096+64] ;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
ss0b:	&ops
IF memused NE 192
	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
ELSE
	lea	rdi, [rdi+XMM_SCD1]	;; Next sine/cosine pointer
ENDIF
	dec	ecx
	jnz	ss0b
	mov	ecx, inner_iters
	dec	eax
	jnz	ss0a
	jmp	exit
	ENDM

avxmac MACRO	memused, memarea, rdi_incr, ops:vararg
	LOCAL	avxlabel
	avxlabel CATSTR <avcase>,%avx_case_num
	avx_case_num = avx_case_num + 1
	avxmac1 avxlabel, memused, memarea, rdi_incr, ops
	ENDM
avxmac1 MACRO	lab, memused, memarea, rdi_incr, ops:vararg
	LOCAL	av00, av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 0			;; Offset for some avx macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf avx macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+262144+4096+64] ;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
av0b:	&ops
	bump	rdi, rdi_incr		;; Next sine/cosine pointer
	dec	ecx
	jnz	av0b
	mov	ecx, inner_iters
	dec	eax
	jnz	av0a
	jmp	exit
	ENDM

avxmacbx MACRO	memused, memarea, rdi_incr, ops:vararg
	LOCAL	avxlabel
	avxlabel CATSTR <avcase>,%avx_case_num
	avx_case_num = avx_case_num + 1
	avxmacbx1 avxlabel, memused, memarea, rdi_incr, ops
	ENDM
avxmacbx1 MACRO	lab, memused, memarea, rdi_incr, ops:vararg
	LOCAL	av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 262144+128		;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+262144+4096+64] ;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
av0b:	&ops
	bump	rdi, rdi_incr		;; Next sine/cosine pointer
	dec	ecx
	jnz	av0b
	mov	ecx, inner_iters
	dec	eax
	jnz	av0a
	jmp	exit
	ENDM

ynormmac MACRO	memused, memarea, ops:vararg
	LOCAL	avxlabel
	avxlabel CATSTR <avcase>,%avx_case_num
	avx_case_num = avx_case_num + 1
	ynormmac1 avxlabel, memused, memarea, ops
	ENDM
ynormmac1 MACRO	lab, memused, memarea, ops:vararg
	LOCAL	av00, av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	edx, outer_iters
	mov	ebx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rsi, SRCARG		;; Reload work buf addr (FFT data)
	lea	rbp, [rsi+524288+256]	;; Addr of the multipliers
;	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	lea	rdi, [rsi+5000000]
av0b:	&ops
	bump	rsi, 64			;; Next cache line
	bump	rbp, 128		;; Next set of 8 multipliers
	bump	rdi, 2			;; Next big/little flags
	dec	ebx
	jnz	av0b
	mov	ebx, inner_iters
	dec	edx
	jnz	av0a
	jmp	exit
	ENDM

ynormwpnmac MACRO memused, memarea, ops:vararg
	LOCAL	avxlabel
	avxlabel CATSTR <avcase>,%avx_case_num
	avx_case_num = avx_case_num + 1
	ynormwpnmac1 avxlabel, memused, memarea, ops
	ENDM
ynormwpnmac1 MACRO lab, memused, memarea, ops:vararg
	LOCAL	av00, av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:
	mov	COPYZERO, outer_iters
	mov	ebp, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rsi, SRCARG		;; Reload work buf addr (FFT data)
	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
IFDEF X86_64
	lea	r9, [rdx+2*YMM_GMD]	;; Prefetch pointer for group multipliers
ENDIF
	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
av0b:	&ops
	bump	rsi, 64			;; Next cache line
	bump	rdi, 2			;; Next big/little flags
;	bump	rdx, 2*YMM_GMD		;; Next set of group multipliers
	bump	rdx, 64
	dec	rbp
	jnz	av0b
	mov	ebp, inner_iters
	dec	COPYZERO
	jnz	av0a
	jmp	exit
	ENDM


IFNDEF X86_64
ynormwpn4mac MACRO memused, memarea, ops:vararg
	ynormwpnmac memused, memarea, ops
	ENDM
ELSE
ynormwpn4mac MACRO memused, memarea, ops:vararg
	LOCAL	avxlabel
	avxlabel CATSTR <avcase>,%avx_case_num
	avx_case_num = avx_case_num + 1
	ynormwpn4mac1 avxlabel, memused, memarea, ops
	ENDM
ynormwpn4mac1 MACRO lab, memused, memarea, ops:vararg
	LOCAL	av00, av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:
	mov	COPYZERO, outer_iters
	mov	ebp, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rsi, SRCARG		;; Reload work buf addr (FFT data)
	lea	r13, [rsi+8192+256]	;; Source ptr #2
	mov	r12, norm_grp_mults	;; Addr of the group multipliers
	lea	r15, [r12+2*YMM_GMD]	;; Addr of the group multipliers #2
	lea	r9, [r15+2*YMM_GMD]	;; Prefetch pointer for group multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	lea	r14, [rdi+4096+128]	;; Addr of the big/little flags array #2
	movzx	rbx, WORD PTR [rdi]	;; Preload 4 big vs. little & fudge flags
	movzx	rcx, WORD PTR [r14]	;; Preload 4 big vs. little & fudge flags
av0b:	&ops
	bump	rsi, 64			;; Next cache line
	bump	rdi, 2			;; Next big/little flags
	bump	r13, 64			;; Next cache line
	bump	r14, 2			;; Next big/little flags
;	bump	r12, 4*YMM_GMD		;; Next set of group multipliers
	bump	r12, 64
	bump	r15, 64
	dec	rbp
	jnz	av0b
	mov	ebp, inner_iters
	dec	COPYZERO
	jnz	av0a
	jmp	exit
	ENDM
ENDIF

avx512mac MACRO	memused, memarea, rdi_incr, ops:vararg
	LOCAL	avx512label
	avx512label CATSTR <av512case>,%avx512_case_num
	avx512_case_num = avx512_case_num + 1
	avx512mac1 avx512label, memused, memarea, rdi_incr, ops
	ENDM
avx512mac1 MACRO lab, memused, memarea, rdi_incr, ops:vararg
	LOCAL	av00, av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 0			;; Offset for some avx macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf avx macros
	mov	r12, norm_grp_mults
	mov	r13, r12
	mov	r8, 128			;; d1reg for "rsc" macros
	mov	r9, 384			;; d3reg for "rsc" macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+262144+4096+64] ;; Source & dest ptr
	lea	rdx, [rsi+524288+8192+256] ;; Destination for "rsc" macros
	lea	r10, [rdx+512]		;; src4reg/dst4reg for "rsc" macros
	lea	r15, [rdx+262144+256+64] ;; s/c data for first macros
IF avx512_case_num-1 EQ IACA_TEST_CASE
IACA_START
ENDIF
av0b:	&ops
	bump	rdi, rdi_incr		;; Next sine/cosine pointer
	dec	ecx
	jnz	av0b
IF avx512_case_num-1 EQ IACA_TEST_CASE
IACA_END
ENDIF
	mov	ecx, inner_iters
	dec	eax
	jnz	av0a
	jmp	exit
	ENDM

avx512znormmac MACRO	memused, memarea, rdi_incr, ops:vararg
	LOCAL	avx512label
	avx512label CATSTR <av512case>,%avx512_case_num
	avx512_case_num = avx512_case_num + 1
	avx512znormmac1 avx512label, memused, memarea, rdi_incr, ops
	ENDM
avx512znormmac1 MACRO lab, memused, memarea, rdi_incr, ops:vararg
	LOCAL	av00, av0a, av0b
	inner_iters = memarea / (memused)
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 0			;; Offset for some avx macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf avx macros
	mov	r12, norm_grp_mults
	mov	r13, r12
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
av0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+262144+4096+64] ;; Source & dest ptr
	mov	r14, 262144+4096+64
	lea	r13, [rsi+r14*2]
	lea	r8, [rdi+4096]
	lea	r10, [r8+4096+64]
	sub	rdx, rdx
IF avx512_case_num-1 EQ IACA_TEST_CASE
IACA_START
ENDIF
av0b:	&ops
	bump	rsi, 128
	bump	r13, 128
	bump	rdi, 1
	dec	ecx
	jnz	av0b
IF avx512_case_num-1 EQ IACA_TEST_CASE
IACA_END
ENDIF
	mov	ecx, inner_iters
	dec	eax
	jnz	av0a
	jmp	exit
	ENDM

do2 MACRO ops:vararg
	&ops
	&ops
	ENDM
do5 MACRO ops:vararg
	do2 &ops
	do2 &ops
	&ops
	ENDM
do10 MACRO ops:vararg
	do2 do5 &ops
	ENDM

fiveclocks MACRO
	vaddpd	ymm0, ymm0, ymm0
	vaddpd	ymm1, ymm1, ymm1
	vaddpd	ymm2, ymm2, ymm2
	vaddpd	ymm3, ymm3, ymm3
	vaddpd	ymm4, ymm4, ymm4
	vaddpd	ymm5, ymm5, ymm5
	vaddpd	ymm6, ymm6, ymm6
	vaddpd	ymm7, ymm7, ymm7
	vaddpd	ymm8, ymm8, ymm8
	vaddpd	ymm9, ymm9, ymm9
	ENDM
tenclocks MACRO
	fiveclocks
	fiveclocks
	ENDM
hundredclocks MACRO
	do10 tenclocks
	ENDM

fiveclocks512 MACRO
	vaddpd	zmm0, zmm0, zmm0
	vaddpd	zmm1, zmm1, zmm1
	vaddpd	zmm2, zmm2, zmm2
	vaddpd	zmm3, zmm3, zmm3
	vaddpd	zmm4, zmm4, zmm4
	vaddpd	zmm5, zmm5, zmm5
	vaddpd	zmm6, zmm6, zmm6
	vaddpd	zmm7, zmm7, zmm7
	vaddpd	zmm8, zmm8, zmm8
	vaddpd	zmm9, zmm9, zmm9
	ENDM
tenclocks512 MACRO
	fiveclocks512
	fiveclocks512
	ENDM
hundredclocks512 MACRO
	do10 tenclocks512
	ENDM

_TEXT	SEGMENT

x87table: looptab case, X87_CASES
sse2table: looptab sscase, SSE2_CASES
avxtable: looptab avcase, AVX_CASES
avx512table: looptab av512case, AVX512_CASES

; gwtimeit (asm_data)
;	Time a mini benchmark
; Windows 32-bit (_gwtimeit)
; Linux 32-bit (gwtimeit)
;	Parameter asm_data = [esp+4]
; Windows 64-bit (gwtimeit)
;	Parameter asm_data = rcx
; Linux 64-bit (gwtimeit)
;	Parameter asm_data = rdi

PROCFL	gwtimeit
	ad_prolog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r11,r12,r13,r14,r15,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

	mov	rdi, PPTR [AD_BASE]	; Load work buf address
	mov	rsi, rdi

	mov	edx, DPTR [AD_BASE+8]	; Load n (which test to run)

	sub	ebx, ebx		; Clear registers
	sub	ecx, ecx
	sub	ebp, ebp

	mov	eax, X87_CASES
	cmp	edx, -1			; -1 = get num x87 cases
	je	exit
	mov	eax, SSE2_CASES
	cmp	edx, -2			; -2 = get num sse2 cases
	je	exit
	mov	eax, AVX_CASES
	cmp	edx, -3			; -3 = get num avx cases
	je	exit
	mov	eax, AVX512_CASES
	cmp	edx, -4			; -4 = get num avx512 cases
	je	exit

	cmp	edx, 1000		; Tests below 1000 are x87 code
	jl	x87

	cmp	edx, 2000		; Tests below 2000 are SSE2 code
	jl	short sse2

	cmp	edx, 3000		; Tests below 3000 are AVX code
	jl	short avx

; Init registers and jump to desired AVX test case

	vzeroall			; Clear ZMM registers
	sub	rdx, 3000
	mov	rax, OFFSET avx512table
	mov	rax, [rax+rdx*SZPTR]; Get address of test to execute
	jmp	rax

; Init registers and jump to desired AVX test case

avx:	vzeroall			; Clear YMM registers
	sub	rdx, 2000
	mov	rax, OFFSET avxtable
	mov	rax, [rax+rdx*SZPTR]; Get address of test to execute
	jmp	rax

; Init registers and jump to desired SSE2 test case

sse2:	subpd	xmm0, xmm0		; Clear XMM registers
	subpd	xmm1, xmm1
	subpd	xmm2, xmm2
	subpd	xmm3, xmm3
	subpd	xmm4, xmm4
	subpd	xmm5, xmm5
	subpd	xmm6, xmm6
	subpd	xmm7, xmm7
IFDEF X86_64
	subpd	xmm8, xmm8
	subpd	xmm9, xmm9
	subpd	xmm10, xmm10
	subpd	xmm11, xmm11
	subpd	xmm12, xmm12
	subpd	xmm13, xmm13
	subpd	xmm14, xmm14
	subpd	xmm15, xmm15
ENDIF
	sub	rdx, 1000
	mov	rax, OFFSET sse2table
	mov	rax, [rax+rdx*SZPTR]; Get address of test to execute
	jmp	rax

; Jump to desired x87 test case

x87:	mov	rax, OFFSET x87table
	mov	rax, [rax+rdx*SZPTR]; Get address of test to execute
	jmp	rax

; Time the "do-nothing" case.

IFNDEF X86_64
case0:
	jmp	exit

; Time the loop of "do-nothing" case.  1000 iterations.

case1:
	mov	ecx, 1000
c1a:	dec	ecx
	jnz	c1a
	jmp	exit

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

case2:	read1	4096, 1000	; Read 4KB
	jmp	exit

case3:	read1	96*1024, 100	; Read 96KB
	jmp	exit

case4:	read1	2048*1024, 2	; Read 2MB
	jmp	exit

case5:	x87mac	64, 4096, eight_reals_fft, 8, 16, 32
	jmp	exit

case6:	x87mac	64, 100000, eight_reals_fft, 8, 16, 32
	jmp	exit

case7:	x87mac	64, 4096, eight_reals_unfft, 8, 16, 32
	jmp	exit

case8:	x87mac	64, 100000, eight_reals_unfft, 8, 16, 32
	jmp	exit

case9:	x87mac	64, 4096, four_complex_fft, 8, 16, 32
	jmp	exit

case10:	x87mac	64, 100000, four_complex_fft, 8, 16, 32
	jmp	exit

case11:	x87mac	64, 4096, four_complex_unfft, 8, 16, 32
	jmp	exit

case12:	x87mac	64, 100000, four_complex_unfft, 8, 16, 32
	jmp	exit
ENDIF

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 32MB will test reading
; from main memory.

sscase0:
	read2	4096, 1000	; Read 4KB
	jmp	exit

sscase1:
	read2	96*1024, 100	; Read 96KB
	jmp	exit

sscase2:
	read2	32768*1024, 2	; Read 32MB
	jmp	exit

; This code writes a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 32MB will test writing
; to main memory.

sscase3:
	write1	4096, 1000	; Write 4KB
	jmp	exit

sscase4:
	write1	96*1024, 100	; Write 96KB
	jmp	exit

sscase5:
	write1	32768*1024, 2	; Write 32MB
	jmp	exit

; This code writes a block of memory non-contiguously.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 32MB will test writing
; to main memory.

sscase6:
	write2	4096, 1000	; Read 4KB
	jmp	exit

sscase7:
	write2	96*1024, 100	; Read 96KB
	jmp	exit

sscase8:
	write2	32768*1024, 2	; Read 32MB
	jmp	exit

; This code reads & writes a block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 32MB will test writing
; to main memory.

sscase9:
	readwrite1	4096, 1000	; Read 4KB
	jmp	exit

sscase10:
	readwrite1	96*1024, 100	; Read 96KB
	jmp	exit

sscase11:
	readwrite1	32768*1024, 2	; Read 32MB
	jmp	exit

; Time ~10000 iterations of the SSE2 macros in L1 and L2 caches

	sse2mac sscase12, 128, 4096, x2cl_eight_reals_fft rsi, 2*64, 64
	sse2mac sscase13, 128, 100000, x2cl_eight_reals_fft rsi, 2*64, 64
	sse2mac sscase14, 128, 4096, x2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase15, 128, 100000, x2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase16, 128, 4096, x2cl_eight_reals_fft_2 rsi, 2*64, 64
	sse2mac sscase17, 128, 100000, x2cl_eight_reals_fft_2 rsi, 2*64, 64
	sse2mac sscase18, 128, 4096, x2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase19, 128, 100000, x2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase20, 128, 4096, s2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase21, 128, 100000, s2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase22, 128, 4096, s2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase23, 128, 100000, s2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase24, 128, 4096, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2mac sscase25, 128, 100000, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2mac sscase26, 128, 4096, s2cl_eight_reals_fft_2_final rsi, 2*64, 64
	sse2mac sscase27, 128, 100000, s2cl_eight_reals_fft_2_final rsi, 2*64, 64
	sse2mac sscase28, 128, 4096, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2mac sscase29, 128, 100000, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2macbx sscase30, 128, 4096, s2cl_eight_reals_with_mult_2 rsi, 2*64, 64
	sse2macbx sscase31, 128, 100000, s2cl_eight_reals_with_mult_2 rsi, 2*64, 64
	sse2macbx sscase32, 128, 4096, s2cl_eight_reals_with_mulf_2 rsi, 2*64, 64
	sse2macbx sscase33, 128, 100000, s2cl_eight_reals_with_mulf_2 rsi, 2*64, 64
	sse2mac sscase34, 256, 4096, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase35, 256, 100000, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase36, 256, 4096, g4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	sse2mac sscase37, 256, 100000, g4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64

	sse2mac sscase38, 128, 4096, x2cl_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase39, 128, 100000, x2cl_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase40, 128, 4096, x2cl_half_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase41, 128, 100000, x2cl_half_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase42, 128, 4096, g2cl_eight_reals_unfft_2 rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase43, 128, 100000, g2cl_eight_reals_unfft_2 rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase44, 256, 4096, x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase45, 256, 100000, x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase46, 256, 4096, x4cl_eight_reals_unfft_2 rsi, 4*64, 64, 2*64
	sse2mac sscase47, 256, 100000, x4cl_eight_reals_unfft_2 rsi, 4*64, 64, 2*64
	sse2mac sscase48, 256, 4096, s4cl_eight_reals_unfft_1 rsi, 4*64, 64, 2*64
	sse2mac sscase49, 256, 100000, s4cl_eight_reals_unfft_1 rsi, 4*64, 64, 2*64

	sse2mac sscase50, 128, 4096, x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase51, 128, 100000, x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase52, 128, 4096, x2cl_two_complex_fft_in_place rsi, 2*64, 64, rdi
	sse2mac sscase53, 128, 100000, x2cl_two_complex_fft_in_place rsi, 2*64, 64, rdi

	sse2mac sscase54, 128, 4096, x2cl_two_complex_unfft rsi, 2*64, 64
	sse2mac sscase55, 128, 100000, x2cl_two_complex_unfft rsi, 2*64, 64

	sse2mac sscase56, 128, 4096, x2cl_four_complex_fft rsi, 2*64, 64
	sse2mac sscase57, 128, 100000, x2cl_four_complex_fft rsi, 2*64, 64
	sse2mac sscase58, 128, 4096, x2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase59, 128, 100000, x2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase60, 128, 4096, s2cl_four_complex_gpm_fft rsi, 2*64, 64
	sse2mac sscase61, 128, 100000, s2cl_four_complex_gpm_fft rsi, 2*64, 64
	sse2mac sscase62, 128, 4096, s2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase63, 128, 100000, s2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase64, 128, 4096, s2cl_four_complex_fft_final rsi, 2*64, 64
	sse2mac sscase65, 128, 100000, s2cl_four_complex_fft_final rsi, 2*64, 64
	sse2mac sscase66, 128, 4096, s2cl_four_complex_with_square rsi, 2*64, 64
	sse2mac sscase67, 128, 100000, s2cl_four_complex_with_square rsi, 2*64, 64
	sse2macbx sscase68, 128, 4096, s2cl_four_complex_with_mult rsi, 2*64, 64
	sse2macbx sscase69, 128, 100000, s2cl_four_complex_with_mult rsi, 2*64, 64
	sse2macbx sscase70, 128, 4096, s2cl_four_complex_with_mulf rsi, 2*64, 64
	sse2macbx sscase71, 128, 100000, s2cl_four_complex_with_mulf rsi, 2*64, 64
	sse2mac sscase72, 256, 4096, x4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase73, 256, 100000, x4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase74, 256, 4096, x4cl_four_complex_cpm_fft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase75, 256, 100000, x4cl_four_complex_cpm_fft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase76, 256, 4096, g4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	sse2mac sscase77, 256, 100000, g4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	sse2mac sscase78, 256, 4096, x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2mac sscase79, 256, 100000, x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2macbx sscase80, 256, 4096, x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase81, 256, 100000, x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase82, 256, 4096, x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	sse2macbx sscase83, 256, 100000, x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64

	sse2mac sscase84, 128, 4096, x2cl_four_complex_unfft rsi, 2*64, 64
	sse2mac sscase85, 128, 100000, x2cl_four_complex_unfft rsi, 2*64, 64
	sse2mac sscase86, 128, 4096, g2cl_four_complex_unfft rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase87, 128, 100000, g2cl_four_complex_unfft rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase88, 256, 4096, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase89, 256, 100000, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase90, 256, 4096, x4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase91, 256, 100000, x4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase92, 256, 4096, s4cl_four_complex_gpm_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase93, 256, 100000, s4cl_four_complex_gpm_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase94, 256, 4096, x4cl_four_complex_cpm_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase95, 256, 100000, x4cl_four_complex_cpm_unfft rsi, 4*64, 64, 2*64

	sse2mac sscase96, 192, 4096, x3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase97, 192, 100000, x3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase98, 192, 4096, g3cl_six_reals_first_fft rsi, 3*64, 64, rdx, 3*64, 64
	sse2mac sscase99, 192, 100000, g3cl_six_reals_first_fft rsi, 3*64, 64, rdx, 3*64, 64
	sse2mac sscase100, 192, 4096, s3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase101, 192, 100000, s3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase102, 192, 4096, x3cl_six_reals_last_unfft rsi, 3*64, 64
	sse2mac sscase103, 192, 100000, x3cl_six_reals_last_unfft rsi, 3*64, 64

	sse2mac sscase104, 192, 4096, x3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase105, 192, 100000, x3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase106, 192, 4096, s3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase107, 192, 100000, s3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase108, 192, 4096, x3cl_three_complex_last_unfft rsi, 3*64, 64
	sse2mac sscase109, 192, 100000, x3cl_three_complex_last_unfft rsi, 3*64, 64

	sse2mac sscase110, 320, 4096, x5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase111, 320, 100000, x5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase112, 320, 4096, g5cl_five_reals_first_fft rsi, 5*64, 64, rdx, 5*64, 64
	sse2mac sscase113, 320, 100000, g5cl_five_reals_first_fft rsi, 5*64, 64, rdx, 5*64, 64
	sse2mac sscase114, 320, 4096, s5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase115, 320, 100000, s5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase116, 320, 4096, x5cl_five_reals_last_unfft rsi, 5*64, 64
	sse2mac sscase117, 320, 100000, x5cl_five_reals_last_unfft rsi, 5*64, 64

	sse2mac sscase118, 448, 4096, x7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase119, 448, 100000, x7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase120, 448, 4096, g7cl_seven_reals_first_fft rsi, 7*64, 64, rdx, 7*64, 64
	sse2mac sscase121, 448, 100000, g7cl_seven_reals_first_fft rsi, 7*64, 64, rdx, 7*64, 64
	sse2mac sscase122, 448, 4096, s7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase123, 448, 100000, s7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase124, 448, 4096, x7cl_seven_reals_last_unfft rsi, 7*64, 64
	sse2mac sscase125, 448, 100000, x7cl_seven_reals_last_unfft rsi, 7*64, 64

	sse2mac sscase126, 256, 4096, x4cl_empty rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase127, 256, 100000, x4cl_empty rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase128, 256, 4096, g4cl_empty rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase129, 256, 100000, g4cl_empty rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase130, 256, 4096, g4cl_empty_nt rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase131, 256, 100000, g4cl_empty_nt rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0

	sse2mac sscase132, 128, 4096, r4_x2cl_four_complex_first_fft4_scratch rsi, 2*64, 64, rdi
	sse2mac sscase133, 128, 100000, r4_x2cl_four_complex_first_fft4_scratch rsi, 2*64, 64, rdi
	sse2mac sscase134, 256, 4096, r4_x4cl_four_complex_last_unfft4 rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase135, 256, 100000, r4_x4cl_four_complex_last_unfft4 rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase136, 256, 4096, r4_x4cl_four_complex_fft_final rsi, 4*64, 64, 2*64
	sse2mac sscase137, 256, 100000, r4_x4cl_four_complex_fft_final rsi, 4*64, 64, 2*64
	sse2mac sscase138, 256, 4096, r4_x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2mac sscase139, 256, 100000, r4_x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2macbx sscase140, 256, 4096, r4_x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase141, 256, 100000, r4_x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase142, 256, 4096, r4_x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	sse2macbx sscase143, 256, 100000, r4_x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64

	sse2mac sscase144, 128, 4096, r2_x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase145, 128, 100000, r2_x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase146, 128, 4096, r2_x2cl_two_complex_unfft rsi, 2*64, 64, rdi, 0
	sse2mac sscase147, 128, 100000, r2_x2cl_two_complex_unfft rsi, 2*64, 64, rdi, 0

	sse2mac sscase148, 128, 4096, r4_x2cl_four_complex_first_djbfft_scratch rsi, 2*64, 64, rdi, rdi
	sse2mac sscase149, 128, 100000, r4_x2cl_four_complex_first_djbfft_scratch rsi, 2*64, 64, rdi, rdi
	sse2mac sscase150, 256, 4096, r4_x4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase151, 256, 100000, r4_x4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase152, 256, 4096, r4_x4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase153, 256, 100000, r4_x4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase154, 256, 4096, r4_sg4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase155, 256, 100000, r4_sg4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2macbx sscase156, 128, 4096, r4_f2cl_four_complex_djbfft rsi, 2*64, 64, rdi
	sse2macbx sscase157, 128, 100000, r4_f2cl_four_complex_djbfft rsi, 2*64, 64, rdi

	sse2mac sscase158, 256, 4096, r4_x4cl_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase159, 256, 100000, r4_x4cl_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase160, 256, 4096, r4_x4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase161, 256, 100000, r4_x4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase162, 256, 4096, r4_x4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase163, 256, 100000, r4_x4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase164, 128, 4096, r4_sg2cl_four_complex_djbunfft rsi, 2*64, 64, rdx, 2*64, 64, rdi
	sse2mac sscase165, 128, 100000, r4_sg2cl_four_complex_djbunfft rsi, 2*64, 64, rdx, 2*64, 64, rdi

	sse2mac sscase166, 256, 4096, r4_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase167, 256, 100000, r4_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase168, 128, 4096, r4_sg2cl_four_complex_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi
	sse2mac sscase169, 128, 100000, r4_sg2cl_four_complex_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi

	sse2mac sscase170, 128, 4096, r4_x2cl_eight_reals_first_fft_scratch rsi, 2*64, 64, rdi
	sse2mac sscase171, 128, 100000, r4_x2cl_eight_reals_first_fft_scratch rsi, 2*64, 64, rdi
	sse2mac sscase172, 256, 4096, r4_x4cl_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase173, 256, 100000, r4_x4cl_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase174, 256, 4096, r4_x4cl_wpn_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase175, 256, 100000, r4_x4cl_wpn_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase176, 256, 4096, r4_sg4cl_eight_reals_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase177, 256, 100000, r4_sg4cl_eight_reals_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, rdi

	sse2mac sscase178, 256, 4096, r4_x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase179, 256, 100000, r4_x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase180, 256, 4096, r4_x4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase181, 256, 100000, r4_x4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase182, 128, 4096, r4_sg2cl_eight_reals_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi, rdi
	sse2mac sscase183, 128, 100000, r4_sg2cl_eight_reals_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi, rdi

	sse2mac sscase184, 192, 4096, r3_x3cl_six_reals_three_complex_djbfft rsi, 3*64, 64, rdi, rdi
	sse2mac sscase185, 192, 100000, r3_x3cl_six_reals_three_complex_djbfft rsi, 3*64, 64, rdi, rdi
	sse2mac sscase186, 192, 4096, r3_x3cl_six_reals_unfft rsi, 3*64, 64, rdi, 0, rdi, 0
	sse2mac sscase187, 192, 100000, r3_x3cl_six_reals_unfft rsi, 3*64, 64, rdi, 0, rdi, 0

	sse2mac sscase188, 192, 4096, r3_x3cl_three_complex_djbfft rsi, 3*64, 64, rdi
	sse2mac sscase189, 192, 100000, r3_x3cl_three_complex_djbfft rsi, 3*64, 64, rdi
	sse2mac sscase190, 192, 4096, r3_x3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, 0
	sse2mac sscase191, 192, 100000, r3_x3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, 0

	sse2mac sscase192, 320, 4096, r5_x5cl_20_reals_first_fft_scratch rsi, 5*64, 64, rdi
	sse2mac sscase193, 320, 100000, r5_x5cl_20_reals_first_fft_scratch rsi, 5*64, 64, rdi
	sse2mac sscase194, 640, 4096, r5_x10cl_20_reals_last_unfft rsi, 10*64, 64, rdi, 8*XMM_SCD9
	sse2mac sscase195, 640, 100000, r5_x10cl_20_reals_last_unfft rsi, 10*64, 64, rdi, 8*XMM_SCD9

	sse2mac sscase196, 320, 4096, r5_x5cl_five_complex_djbfft rsi, 5*64, 64, rdi
	sse2mac sscase197, 320, 100000, r5_x5cl_five_complex_djbfft rsi, 5*64, 64, rdi
	sse2mac sscase198, 320, 4096, r5_x5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, 0
	sse2mac sscase199, 320, 100000, r5_x5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, 0

	sse2mac sscase200, 448, 4096, r7_x7cl_28_reals_first_fft_scratch rsi, 7*64, 64, rdi
	sse2mac sscase201, 448, 100000, r7_x7cl_28_reals_first_fft_scratch rsi, 7*64, 64, rdi
	sse2mac sscase202, 896, 4096, r7_x14cl_28_reals_last_unfft rsi, 14*64, 64, rdi, 8*XMM_SCD13
	sse2mac sscase203, 896, 100000, r7_x14cl_28_reals_last_unfft rsi, 14*64, 64, rdi, 8*XMM_SCD13

	sse2mac sscase204, 512, 4096, r8_x8cl_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	sse2mac sscase205, 512, 100000, r8_x8cl_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	sse2mac sscase206, 512, 4096, r8_x8cl_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	sse2mac sscase207, 512, 100000, r8_x8cl_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase208, 512, 4096, r8_x8cl_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase209, 512, 100000, r8_x8cl_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase210, 512, 4096, r8_x8cl_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase211, 512, 100000, r8_x8cl_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64

	sse2mac sscase212, 512, 4096, r8_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi
	sse2mac sscase213, 512, 100000, r8_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi
	sse2mac sscase214, 256, 4096, r8_sg4cl_eight_complex_unfft8 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase215, 256, 100000, r8_sg4cl_eight_complex_unfft8 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi

	; need 16 reals cases

; Time ~10000 iterations of the AVX macros in L1 and L2 caches

INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE ymult.mac
INCLUDE yr4.mac
INCLUDE ynormal.mac

; This code reads/writes 64MB (1M cache lines) in contiguous blocks.  Timings are done
; on 4 memory sizes.  4KB will operate on the L1 cache only, 128KB will operate on the
; L2 cache only, 1MB will operate on the L3 cache and 32MB will operate on main memory.

avcase0:	read4	4096, 16384	; Read 4KB
		jmp	exit
avcase1:	read4	128*1024, 512	; Read 128KB
		jmp	exit
avcase2:	read4	1024*1024, 64	; Read 1MB
		jmp	exit
avcase3:	read4	32768*1024, 2	; Read 32MB
		jmp	exit

avcase4:	write4	4096, 16384	; Write 4KB
		jmp	exit
avcase5:	write4	128*1024, 512	; Write 128KB
		jmp	exit
avcase6:	write4	1024*1024, 64	; Write 1MB
		jmp	exit
avcase7:	write4	32768*1024, 2	; Write 32MB
		jmp	exit

avcase8:	readwrite4 4096, 16384	; Read/write 4KB
		jmp	exit
avcase9:	readwrite4 128*1024, 512 ; Read/write 128KB
		jmp	exit
avcase10:	readwrite4 1024*1024, 64 ; Read/write 1MB
		jmp	exit
avcase11:	readwrite4 32768*1024, 2 ; Read/write 32MB
		jmp	exit

	yloop_init 32			;; Dummy call to yloop_init

	avx_case_num = 12

;;12
	avxmac 256*1, 8192, 0, yr4_4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*2, 8192, 0, yr4_4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 2
	avxmac 256*1, 100000, 0, yr4_4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*2, 100000, 0, yr4_4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 2
	avxmac 256*1, 8192, 0, yr4_s4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*1, 100000, 0, yr4_s4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmacbx 256*1, 8192, 0, yr4_fs4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmacbx 256*1, 100000, 0, yr4_fs4cl_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256, 8192, 0, yr4_b4cl_csc_wpn_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256, 100000, 0, yr4_b4cl_csc_wpn_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256*1, 8192, 0, yr4_sg4cl_2sc_eight_reals_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1
	avxmac 256*1, 100000, 0, yr4_sg4cl_2sc_eight_reals_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1

;;24
	avxmac 256*1, 8192, 0, yr4_4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*2, 8192, 0, yr4_4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 2
	avxmac 256*1, 100000, 0, yr4_4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*2, 100000, 0, yr4_4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 2
	avxmac 256*1, 8192, 0, yr4_s4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*1, 100000, 0, yr4_s4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256, 8192, 0, yr4_b4cl_csc_wpn_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256, 100000, 0, yr4_b4cl_csc_wpn_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256, 8192, 0, yr4_sg4cl_2sc_eight_reals_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1
	avxmac 256, 100000, 0, yr4_sg4cl_2sc_eight_reals_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1

;;34
	avxmac 256*1, 8192, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1	;;, L1PREFETCH_ALL, 2*256*1
	avxmac 256*2, 8192, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2	;;, L1PREFETCH_ALL, 2*256*2
	avxmac 256*3, 8192, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 3	;;, L1PREFETCH_ALL, 2*256*3
	avxmac 256*4, 8192, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 4	;;, L1PREFETCH_ALL, 2*256*4
	avxmac 256*5, 8192, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 5	;;, L1PREFETCH_ALL, 2*256*5
	avxmac 256*1, 100000, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 100000, 0, yr4_4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmac 256*1, 8192, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 1
	avxmac 256*2, 8192, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 2
	avxmac 256*3, 8192, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 3
	avxmac 256*4, 8192, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 4
	avxmac 256*5, 8192, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 5
	avxmac 256*1, 100000, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 1
	avxmac 256*2, 100000, 0, yr4_b4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 2
	avxmac 256*1, 8192, 0, yr4_rb4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*2, 8192, 0, yr4_rb4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 2
	avxmac 256*1, 100000, 0, yr4_rb4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256*2, 100000, 0, yr4_rb4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 2

;;52
	avxmac 256*1, 8192, 0, yr4_s4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 8192, 0, yr4_s4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmac 256*3, 8192, 0, yr4_s4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 3
	avxmac 256*4, 8192, 0, yr4_s4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 4
	avxmac 256*1, 100000, 0, yr4_s4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 100000, 0, yr4_s4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmacbx 256*1, 8192, 0, yr4_fs4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmacbx 256*2, 8192, 0, yr4_fs4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmacbx 256*1, 100000, 0, yr4_fs4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmacbx 256*2, 100000, 0, yr4_fs4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2

;;62
	avxmac 256*1, 8192, 0, yr4_b4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 8192, 0, yr4_b4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 100000, 0, yr4_b4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 100000, 0, yr4_b4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 8192, 0, yr4_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 100000, 0, yr4_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1

;;68
	avxmac 256*1, 8192, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 8192, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmac 256*3, 8192, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 3
	avxmac 256*4, 8192, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 4
	avxmac 256*5, 8192, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 5
	avxmac 256*1, 100000, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 100000, 0, yr4_4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmac 256*1, 8192, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 1
	avxmac 256*2, 8192, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 2
	avxmac 256*3, 8192, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 3
	avxmac 256*4, 8192, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 4
	avxmac 256*5, 8192, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 5
	avxmac 256*1, 100000, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 1
	avxmac 256*2, 100000, 0, yr4_b4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2/4, 2

;; 82
	avxmac 256*1, 8192, 0, yr4_s4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 8192, 0, yr4_s4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmac 256*3, 8192, 0, yr4_s4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 3
	avxmac 256*4, 8192, 0, yr4_s4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 4
	avxmac 256*1, 100000, 0, yr4_s4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 1
	avxmac 256*2, 100000, 0, yr4_s4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD2, 2
	avxmac 256*1, 8192, 0, yr4_b4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 8192, 0, yr4_b4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 100000, 0, yr4_b4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 100000, 0, yr4_b4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 8192, 0, yr4_sg4cl_four_complex_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 100000, 0, yr4_sg4cl_four_complex_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1

;;94
	avxmac 256*1, 8192, 0, yr4_4cl_csc_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 1
	avxmac 256*2, 8192, 0, yr4_4cl_csc_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 2
	avxmac 256*1, 100000, 0, yr4_4cl_csc_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 1
	avxmac 256*2, 100000, 0, yr4_4cl_csc_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 2
	avxmac 256*1, 8192, 0, yr4_fs4cl_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 100000, 0, yr4_fs4cl_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmacbx 256*1, 8192, 0, yr4_fs4cl_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmacbx 256*1, 100000, 0, yr4_fs4cl_four_complex_first_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD4, 1

;;102
	avxmac 256*1, 8192, 0, yr4_4cl_csc_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 1
	avxmac 256*2, 8192, 0, yr4_4cl_csc_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 2
	avxmac 256*1, 100000, 0, yr4_4cl_csc_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 1
	avxmac 256*2, 100000, 0, yr4_4cl_csc_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD6, 2
	avxmac 256*1, 8192, 0, yr4_s4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 100000, 0, yr4_s4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD4, 1

;;108
	avxmac 256, 8192, 0, yr4_4cl_four_complex_fft_final rsi, 4*64, 64, 2*64
	avxmac 256, 100000, 0, yr4_4cl_four_complex_fft_final rsi, 4*64, 64, 2*64
	avxmac 256*1, 8192, 0, yr4_4cl_four_complex_with_square rsi, 4*64, 64, 2*64, 1
	avxmac 256*2, 8192, 0, yr4_4cl_four_complex_with_square rsi, 4*64, 64, 2*64, 2
	avxmac 256*1, 100000, 0, yr4_4cl_four_complex_with_square rsi, 4*64, 64, 2*64, 1
	avxmac 256*2, 100000, 0, yr4_4cl_four_complex_with_square rsi, 4*64, 64, 2*64, 2
	avxmacbx 256, 8192, 0, yr4_4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	avxmacbx 256, 100000, 0, yr4_4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	avxmacbx 256, 8192, 0, yr4_4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	avxmacbx 256, 100000, 0, yr4_4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64

;;118
	avxmac 256, 8192, 0, yr4_4cl_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256, 100000, 0, yr4_4cl_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256, 8192, 0, yr4_4cl_eight_reals_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1
	avxmac 256, 100000, 0, yr4_4cl_eight_reals_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, YMM_SCD3, 1

;;122
	avxmac 256, 8192, 0, yr4_4cl_eight_reals_four_complex_fft_final rsi, 4*64, 64, 2*64
	avxmac 256, 100000, 0, yr4_4cl_eight_reals_four_complex_fft_final rsi, 4*64, 64, 2*64
	avxmac 256, 8192, 0, yr4_4cl_eight_reals_four_complex_with_square rsi, 4*64, 64, 2*64
	avxmac 256, 100000, 0, yr4_4cl_eight_reals_four_complex_with_square rsi, 4*64, 64, 2*64
	avxmacbx 256, 8192, 0, yr4_4cl_eight_reals_four_complex_with_mult rsi, 4*64, 64, 2*64
	avxmacbx 256, 100000, 0, yr4_4cl_eight_reals_four_complex_with_mult rsi, 4*64, 64, 2*64
	avxmacbx 256, 8192, 0, yr4_4cl_eight_reals_four_complex_with_mulf rsi, 4*64, 64, 2*64
	avxmacbx 256, 100000, 0, yr4_4cl_eight_reals_four_complex_with_mulf rsi, 4*64, 64, 2*64

;;130
	avxmac 192*1, 8192, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 1
	avxmac 192*2, 8192, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 2
	avxmac 192*3, 8192, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 3
	avxmac 192*4, 8192, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 4
	avxmac 192*5, 8192, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 5
	avxmac 192*1, 100000, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 1
	avxmac 192*2, 100000, 0, yr3_3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 2
	avxmacbx 192*1, 8192, 0, yr3_f3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 1
	avxmacbx 192*2, 8192, 0, yr3_f3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 2
	avxmacbx 192*1, 100000, 0, yr3_f3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 1
	avxmacbx 192*2, 100000, 0, yr3_f3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1, 2
	avxmac 192*1, 8192, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 1
	avxmac 192*2, 8192, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 2
	avxmac 192*3, 8192, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 3
	avxmac 192*4, 8192, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 4
	avxmac 192*5, 8192, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 5
	avxmac 192*1, 100000, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 1
	avxmac 192*2, 100000, 0, yr3_b3cl_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 2

;;148
	avxmac 192*1, 8192, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 1
	avxmac 192*2, 8192, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 2
	avxmac 192*3, 8192, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 3
	avxmac 192*4, 8192, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 4
	avxmac 192*5, 8192, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 5
	avxmac 192*1, 100000, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 1
	avxmac 192*2, 100000, 0, yr3_3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1, 2
	avxmac 192*1, 8192, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 1
	avxmac 192*2, 8192, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 2
	avxmac 192*3, 8192, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 3
	avxmac 192*4, 8192, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 4
	avxmac 192*5, 8192, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 5
	avxmac 192*1, 100000, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 1
	avxmac 192*2, 100000, 0, yr3_b3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD1/4, 2

;;162
	avxmac 192*1, 8192, 0, yr3_3cl_six_reals_fft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmac 192*4, 8192, 0, yr3_3cl_six_reals_fft rsi, 3*64, 64, rdi, YMM_SCD2, 4
	avxmacbx 192*1, 8192, 0, yr3_f3cl_six_reals_fft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmacbx 192*1, 100000, 0, yr3_f3cl_six_reals_fft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmac 192*1, 8192, 0, yr3_3cl_six_reals_unfft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmac 192*4, 8192, 0, yr3_3cl_six_reals_unfft rsi, 3*64, 64, rdi, YMM_SCD2, 4

;;168
	avxmac 192, 8192, 0, yr3_3cl_six_reals_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmac 192, 100000, 0, yr3_3cl_six_reals_three_complex_djbfft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmac 192, 8192, 0, yr3_3cl_six_reals_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD2, 1
	avxmac 192, 100000, 0, yr3_3cl_six_reals_three_complex_djbunfft rsi, 3*64, 64, rdi, YMM_SCD2, 1

;;172
	avxmac 320*1, 8192, 0, yr5_5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2, 1
	avxmac 320*2, 8192, 0, yr5_5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2, 2
	avxmac 320*1, 100000, 0, yr5_5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2, 1
	avxmac 320*2, 100000, 0, yr5_5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2, 2
	avxmacbx 320*1, 8192, 0, yr5_f5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2, 1
	avxmacbx 320*1, 100000, 0, yr5_f5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2, 1
	avxmac 320*1, 8192, 0, yr5_b5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 1
	avxmac 320*2, 8192, 0, yr5_b5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 2
	avxmac 320*1, 100000, 0, yr5_b5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 1
	avxmac 320*2, 100000, 0, yr5_b5cl_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 2

;;182
	avxmac 320*1, 8192, 0, yr5_5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2, 1
	avxmac 320*2, 8192, 0, yr5_5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2, 2
	avxmac 320*1, 100000, 0, yr5_5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2, 1
	avxmac 320*2, 100000, 0, yr5_5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2, 2
	avxmac 320*1, 8192, 0, yr5_b5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 1
	avxmac 320*2, 8192, 0, yr5_b5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 2
	avxmac 320*1, 100000, 0, yr5_b5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 1
	avxmac 320*2, 100000, 0, yr5_b5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD2/4, 2

;;190
	avxmac 192, 8192, 0, yr5_5cl_ten_reals_fft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmac 192, 100000, 0, yr5_5cl_ten_reals_fft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmacbx 192, 8192, 0, yr5_f5cl_ten_reals_fft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmacbx 192, 100000, 0, yr5_f5cl_ten_reals_fft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmac 192, 8192, 0, yr5_5cl_ten_reals_unfft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmac 192, 100000, 0, yr5_5cl_ten_reals_unfft rsi, 5*64, 64, rdi, YMM_SCD4, 1

;;196
	avxmac 320, 8192, 0, yr5_5cl_ten_reals_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmac 320, 100000, 0, yr5_5cl_ten_reals_five_complex_djbfft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmac 320, 8192, 0, yr5_5cl_ten_reals_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD4, 1
	avxmac 320, 100000, 0, yr5_5cl_ten_reals_five_complex_djbunfft rsi, 5*64, 64, rdi, YMM_SCD4, 1

;;200
	avxmac 512*1, 8192, 0, yr8_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 8192, 0, yr8_sg8cl_eight_complex_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_sg8cl_eight_complex_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1

;;204
	avxmac 512, 8192, 0, yr8_8cl_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	avxmac 512, 100000, 0, yr8_8cl_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	avxmac 512, 8192, 0, yr8_8cl_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	avxmac 512, 100000, 0, yr8_8cl_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	avxmacbx 512, 8192, 0, yr8_8cl_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	avxmacbx 512, 100000, 0, yr8_8cl_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	avxmacbx 512, 8192, 0, yr8_8cl_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64
	avxmacbx 512, 100000, 0, yr8_8cl_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64

;;212
	avxmac 512*1, 8192, 0, yr8_sg8cl_2sc_sixteen_reals_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_sg8cl_2sc_sixteen_reals_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1
	avxmac 512*1, 8192, 0, yr8_sg8cl_2sc_sixteen_reals_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_sg8cl_2sc_sixteen_reals_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1

;;216
	avxmac 512*1, 8192, 0, yr8_8cl_16_reals_fft rsi, 8*64, 64, 128, 256, rdi, YMM_SCD7, 1
	avxmac 512*1, 100000, 0, yr8_8cl_16_reals_fft rsi, 8*64, 64, 128, 256, rdi, YMM_SCD7, 1
	avxmac 512*1, 8192, 0, yr8_8cl_16_reals_unfft rsi, 8*64, 64, 128, 256, rdi, YMM_SCD7, 1
	avxmac 512*1, 100000, 0, yr8_8cl_16_reals_unfft rsi, 8*64, 64, 128, 256, rdi, YMM_SCD7, 1

;;220
	avxmac 640*1, 8192, 0, yr5_10cl_20_reals_fft rsi, 10*64, 64, 128, rdi, YMM_SCD9, 1
	avxmac 640*1, 100000, 0, yr5_10cl_20_reals_fft rsi, 10*64, 64, 128, rdi, YMM_SCD9, 1
	avxmac 640*1, 8192, 0, yr5_10cl_20_reals_unfft rsi, 10*64, 64, 128, rdi, YMM_SCD9, 1
	avxmac 640*1, 100000, 0, yr5_10cl_20_reals_unfft rsi, 10*64, 64, 128, rdi, YMM_SCD9, 1

;;224
	avxmac 768*1, 8192, 0, yr6_12cl_24_reals_fft rsi, 12*64, 64, 128, rdi, YMM_SCD11, 1
	avxmac 768*1, 100000, 0, yr6_12cl_24_reals_fft rsi, 12*64, 64, 128, rdi, YMM_SCD11, 1
	avxmac 768*1, 8192, 0, yr6_12cl_24_reals_unfft rsi, 12*64, 64, 128, rdi, YMM_SCD11, 1
	avxmac 768*1, 100000, 0, yr6_12cl_24_reals_unfft rsi, 12*64, 64, 128, rdi, YMM_SCD11, 1

;;228
	avxmac 896*1, 8192, 0, yr7_14cl_28_reals_fft rsi, 14*64, 64, 128, rdi, YMM_SCD13, 1
	avxmac 896*1, 100000, 0, yr7_14cl_28_reals_fft rsi, 14*64, 64, 128, rdi, YMM_SCD13, 1
	avxmac 896*1, 8192, 0, yr7_14cl_28_reals_unfft rsi, 14*64, 64, 128, rdi, YMM_SCD13, 1
	avxmac 896*1, 100000, 0, yr7_14cl_28_reals_unfft rsi, 14*64, 64, 128, rdi, YMM_SCD13, 1

;;232
	ynormmac 64, 8192, ynorm_1d exec, exec, noexec, noexec, noexec		;; base2 ttp (the most common case)
	ynormmac 64, 8192, ynorm_1d exec, exec, noexec, exec, noexec		;; base2 ttp echk
	ynormmac 64, 8192, ynorm_1d exec, exec, noexec, noexec, exec		;; base2 ttp const
	ynormmac 64, 8192, ynorm_1d exec, exec, noexec, exec, exec		;; base2 ttp const echk
	ynormmac 64, 8192, ynorm_1d noexec, exec, noexec, noexec, noexec	;; base2 nottp
	ynormmac 64, 8192, ynorm_1d exec, noexec, noexec, noexec, noexec	;; nobase2 ttp
	ynormmac 64, 8192, ynorm_1d exec, noexec, noexec, exec, noexec		;; nobase2 ttp echk
	ynormmac 64, 8192, ynorm_1d exec, noexec, noexec, noexec, exec		;; nobase2 ttp const
	ynormmac 64, 8192, ynorm_1d exec, noexec, noexec, exec, exec		;; nobase2 ttp const echk
	ynormmac 64, 8192, ynorm_1d noexec, noexec, noexec, noexec, noexec	;; nobase2 nottp

;;242
	ynormmac 64, 8192, ynorm_1d_zpad exec, exec, noexec, noexec, noexec, exec, noexec	;; base2 ttp c1 (the most common case)
	ynormmac 64, 8192, ynorm_1d_zpad exec, exec, exec, noexec, noexec, exec, noexec		;; base2 ttp c1 echk
	ynormmac 64, 8192, ynorm_1d_zpad exec, exec, noexec, exec, noexec, noexec, noexec	;; base2 ttp const
	ynormmac 64, 8192, ynorm_1d_zpad exec, exec, exec, exec, noexec, noexec, noexec		;; base2 ttp const echk
	ynormmac 64, 8192, ynorm_1d_zpad exec, exec, noexec, noexec, noexec, exec, noexec	;; base2 ttp c1 khi
	ynormmac 64, 8192, ynorm_1d_zpad exec, noexec, noexec, noexec, noexec, exec, noexec	;; nobase2 ttp c1 (the most common case)
	ynormmac 64, 8192, ynorm_1d_zpad exec, noexec, exec, noexec, noexec, exec, noexec	;; nobase2 ttp c1 echk
	ynormmac 64, 8192, ynorm_1d_zpad exec, noexec, noexec, exec, noexec, noexec, noexec	;; nobase2 ttp const
	ynormmac 64, 8192, ynorm_1d_zpad exec, noexec, exec, exec, noexec, noexec, noexec	;; nobase2 ttp const echk
	ynormmac 64, 8192, ynorm_1d_zpad exec, noexec, noexec, noexec, noexec, exec, noexec	;; nobase2 ttp c1 khi

;;252
	ynormwpn4mac 64, 8192, ynorm_wpn exec, exec, noexec, noexec, noexec		;; base2 ttp (the most common case)
	ynormwpn4mac 64, 8192, ynorm_wpn exec, exec, noexec, exec, noexec		;; base2 ttp echk
	ynormwpn4mac 64, 8192, ynorm_wpn exec, exec, noexec, noexec, exec		;; base2 ttp const
	ynormwpn4mac 64, 8192, ynorm_wpn exec, exec, noexec, exec, exec			;; base2 ttp const echk
	ynormwpn4mac 64, 8192, ynorm_wpn noexec, exec, noexec, noexec, noexec		;; base2 nottp
	ynormwpn4mac 64, 8192, ynorm_wpn exec, noexec, noexec, noexec, noexec		;; nobase2 ttp
	ynormwpn4mac 64, 8192, ynorm_wpn exec, noexec, noexec, exec, noexec		;; nobase2 ttp echk
	ynormwpn4mac 64, 8192, ynorm_wpn exec, noexec, noexec, noexec, exec		;; nobase2 ttp const
	ynormwpn4mac 64, 8192, ynorm_wpn exec, noexec, noexec, exec, exec		;; nobase2 ttp const echk
	ynormwpn4mac 64, 8192, ynorm_wpn noexec, noexec, noexec, noexec, noexec		;; nobase2 nottp

;;262
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, exec, noexec, noexec, noexec, exec, noexec	;; base2 ttp c1 (the most common case)
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, exec, exec, noexec, noexec, exec, noexec	;; base2 ttp c1 echk
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, exec, noexec, exec, noexec, noexec, noexec	;; base2 ttp const
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, exec, exec, exec, noexec, noexec, noexec	;; base2 ttp const echk
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, exec, noexec, noexec, noexec, exec, noexec	;; base2 ttp c1 khi
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, noexec, noexec, noexec, noexec, exec, noexec ;; nobase2 ttp c1 (the most common case)
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, noexec, exec, noexec, noexec, exec, noexec	;; nobase2 ttp c1 echk
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, noexec, noexec, exec, noexec, noexec, noexec ;; nobase2 ttp const
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, noexec, exec, exec, noexec, noexec, noexec	;; nobase2 ttp const echk
	ynormwpn4mac 64, 8192, ynorm_wpn_zpad exec, noexec, noexec, noexec, noexec, exec, noexec ;; nobase2 ttp c1 khi

;;272
	avxmac 256, 8192, 0, yr4_b4cl_csc_wpn4_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256, 100000, 0, yr4_b4cl_csc_wpn4_eight_reals_fft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256, 8192, 0, yr4_b4cl_csc_wpn4_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256, 100000, 0, yr4_b4cl_csc_wpn4_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND4, 1
	avxmac 256*1, 8192, 0, yr4_b4cl_wpn4_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 8192, 0, yr4_b4cl_wpn4_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 100000, 0, yr4_b4cl_wpn4_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 100000, 0, yr4_b4cl_wpn4_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 8192, 0, yr4_b4cl_wpn4_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 8192, 0, yr4_b4cl_wpn4_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2
	avxmac 256*1, 100000, 0, yr4_b4cl_wpn4_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 1
	avxmac 256*2, 100000, 0, yr4_b4cl_wpn4_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, YMM_SCND2, 2

;;284
	avxmac 256*1, 8192, 0, yr4_rsc_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 100000, 0, yr4_rsc_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 8192, 0, yr4_rsc_sg4cl_four_complex_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 100000, 0, yr4_rsc_sg4cl_four_complex_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, 1
	avxmac 256*1, 8192, 0, yr4_rsc_sg4cl_2sc_eight_reals_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1
	avxmac 256*1, 100000, 0, yr4_rsc_sg4cl_2sc_eight_reals_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1
	avxmac 256*1, 8192, 0, yr4_rsc_sg4cl_2sc_eight_reals_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1
	avxmac 256*1, 100000, 0, yr4_rsc_sg4cl_2sc_eight_reals_unfft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, YMM_SCD4, rdi, YMM_SCD2, 1

;;292
	avxmac 512*1, 8192, 0, yr8_rsc_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_rsc_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 8192, 0, yr8_rsc_sg8cl_eight_complex_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_rsc_sg8cl_eight_complex_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, YMM_SCD8, 1
	avxmac 512*1, 8192, 0, yr8_rsc_sg8cl_2sc_sixteen_reals_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_rsc_sg8cl_2sc_sixteen_reals_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1
	avxmac 512*1, 8192, 0, yr8_rsc_sg8cl_2sc_sixteen_reals_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1
	avxmac 512*1, 100000, 0, yr8_rsc_sg8cl_2sc_sixteen_reals_unfft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi, 0, rdi, YMM_SCD8, 1



; Time ~10000 iterations of the AVX-512 macros in L1 and L2 caches

IFDEF X86_64
ARCH EQU SKX
INCLUDE zarch.mac
INCLUDE zbasics.mac
INCLUDE zmult.mac
INCLUDE zonepass.mac
INCLUDE zr4.mac
INCLUDE znormal.mac
INCLUDE znormal_zpad.mac

FMApenalty_test1 MACRO		;; Test if latency is really 4 clocks, ideal macro timing is 4 clocks
	zfmaddpd	zmm0, zmm0, zmm0, zmm0
	zfmaddpd	zmm1, zmm1, zmm1, zmm1
	zfmaddpd	zmm2, zmm2, zmm2, zmm2
	zfmaddpd	zmm3, zmm3, zmm3, zmm3
	zfmaddpd	zmm4, zmm4, zmm4, zmm4
	zfmaddpd	zmm5, zmm5, zmm5, zmm5
	zfmaddpd	zmm6, zmm6, zmm6, zmm6
	zfmaddpd	zmm7, zmm7, zmm7, zmm7
	ENDM
FMApenalty_test2 MACRO		;; Test if *average* latency is really 5 clocks, ideal macro timing is 5 clocks
	zfmaddpd	zmm0, zmm0, zmm0, zmm0
	zfmaddpd	zmm1, zmm1, zmm1, zmm1
	zfmaddpd	zmm2, zmm2, zmm2, zmm2
	zfmaddpd	zmm3, zmm3, zmm3, zmm3
	zfmaddpd	zmm4, zmm4, zmm4, zmm4
	zfmaddpd	zmm5, zmm5, zmm5, zmm5
	zfmaddpd	zmm6, zmm6, zmm6, zmm6
	zfmaddpd	zmm7, zmm7, zmm7, zmm7
	zfmaddpd	zmm8, zmm8, zmm8, zmm8
	zfmaddpd	zmm9, zmm9, zmm9, zmm9
	ENDM
FMApenalty_test3 MACRO		;; Test if *average* latency is 6 clocks, ideal macro timing is 6 clocks
	zfmaddpd	zmm0, zmm0, zmm0, zmm0
	zfmaddpd	zmm1, zmm1, zmm1, zmm1
	zfmaddpd	zmm2, zmm2, zmm2, zmm2
	zfmaddpd	zmm3, zmm3, zmm3, zmm3
	zfmaddpd	zmm4, zmm4, zmm4, zmm4
	zfmaddpd	zmm5, zmm5, zmm5, zmm5
	zfmaddpd	zmm6, zmm6, zmm6, zmm6
	zfmaddpd	zmm7, zmm7, zmm7, zmm7
	zfmaddpd	zmm8, zmm8, zmm8, zmm8
	zfmaddpd	zmm9, zmm9, zmm9, zmm9
	zfmaddpd	zmm10, zmm10, zmm10, zmm10
	zfmaddpd	zmm11, zmm11, zmm11, zmm11
	ENDM
FMApenalty_test4 MACRO		;; Test if there is a penalty switching between FMA and add/sub/mul, ideal macro timing is 6 clocks
	zfmaddpd	zmm0, zmm0, zmm0, zmm0
	zfmaddpd	zmm1, zmm1, zmm1, zmm1
	zfmaddpd	zmm2, zmm2, zmm2, zmm2
	zfmaddpd	zmm3, zmm3, zmm3, zmm3
	zfmaddpd	zmm12, zmm12, zmm12, zmm12
	zfmaddpd	zmm13, zmm13, zmm13, zmm13
	vaddpd		zmm4, zmm4, zmm4
	vaddpd		zmm5, zmm5, zmm5
	vaddpd		zmm6, zmm6, zmm6
	vaddpd		zmm7, zmm7, zmm7
	vaddpd		zmm16, zmm16, zmm16
	vaddpd		zmm17, zmm17, zmm17
	ENDM
FMApenalty_test5 MACRO		;; Test if there is a penalty using FMA vs. add/sub/mul, ideal macro timing is 6 clocks
	vaddpd		zmm0, zmm0, zmm0
	vaddpd		zmm1, zmm1, zmm1
	vaddpd		zmm2, zmm2, zmm2
	vaddpd		zmm3, zmm3, zmm3
	vaddpd		zmm12, zmm12, zmm12
	vaddpd		zmm13, zmm13, zmm13
	vaddpd		zmm4, zmm4, zmm4
	vaddpd		zmm5, zmm5, zmm5
	vaddpd		zmm6, zmm6, zmm6
	vaddpd		zmm7, zmm7, zmm7
	vaddpd		zmm16, zmm16, zmm16
	vaddpd		zmm17, zmm17, zmm17
	ENDM

do_nothing MACRO srcreg,srcinc,d1,d2,d4,screg,scinc,maxrpt,L1pt,L1pd
	vmovapd	zmm1, [srcreg+d1]	;; R2
	vmovapd	zmm2, [srcreg+d4+d1]	;; R6
	vmovapd	zmm3, [srcreg+d2+d1+64]	;; I4
	vmovapd	zmm4, [srcreg+d4+d2+d1+64] ;; I8
	vmovapd	zmm5, [srcreg+d2+d1]	;; R4
	vmovapd	zmm6, [srcreg+d4+d2+d1]	;; R8
	vmovapd	zmm7, [srcreg+d1+64]	;; I2
	vmovapd	zmm8, [srcreg+d4+d1+64]	;; I6
	vmovapd	zmm9, [srcreg]		;; R1
	vmovapd	zmm10, [srcreg+d4]	;; R5
	vmovapd	zmm11, [srcreg+d2]	;; R3
	vmovapd	zmm12, [srcreg+d4+d2]	;; R7
	vmovapd	zmm13, [srcreg+64]	;; I1
	vmovapd	zmm14, [srcreg+d4+64]	;; I5
	vmovapd	zmm15, [srcreg+d2+64]	;; I3
	vmovapd	zmm16, [srcreg+d4+d2+64]	;; I7
	zstore	[srcreg], zmm1			;; Store R1
	zstore	[srcreg+d1], zmm2		;; Store R2
	zstore	[srcreg+d2], zmm3		;; Store R3
	zstore	[srcreg+d2+d1], zmm4		;; Store R4
	zstore	[srcreg+d4], zmm5		;; Store R5
	zstore	[srcreg+d4+d1], zmm6		;; Store R6
	zstore	[srcreg+d4+d2], zmm7		;; Store R7
	zstore	[srcreg+d4+d2+d1], zmm8		;; Store R8
	zstore	[srcreg+64], zmm9		;; Store I1
	zstore	[srcreg+d1+64], zmm10		;; Store I2
	zstore	[srcreg+d2+64], zmm11		;; Store I3
	zstore	[srcreg+d2+d1+64], zmm12	;; Store I4
	zstore	[srcreg+d4+64], zmm13		;; Store I5
	zstore	[srcreg+d4+d1+64], zmm14	;; Store I6
	zstore	[srcreg+d4+d2+64], zmm15	;; Store I7
	zstore	[srcreg+d4+d2+d1+64], zmm16	;; Store I8
	bump	srcreg, srcinc
	ENDM

do_nothing_plus_sincos MACRO srcreg,srcinc,d1,d2,d4,screg,scinc,maxrpt,L1pt,L1pd
	vmovapd	zmm1, [srcreg+d1]	;; R2
	vmovapd	zmm2, [srcreg+d4+d1]	;; R6
	vmovapd	zmm3, [srcreg+d2+d1+64]	;; I4
	vmovapd	zmm4, [srcreg+d4+d2+d1+64] ;; I8
	vmovapd	zmm5, [srcreg+d2+d1]	;; R4
	vmovapd	zmm6, [srcreg+d4+d2+d1]	;; R8
	vmovapd	zmm7, [srcreg+d1+64]	;; I2
	vmovapd	zmm8, [srcreg+d4+d1+64]	;; I6
	vmovapd	zmm9, [srcreg]		;; R1
	vmovapd	zmm10, [srcreg+d4]	;; R5
	vmovapd	zmm11, [srcreg+d2]	;; R3
	vmovapd	zmm12, [srcreg+d4+d2]	;; R7
	vmovapd	zmm13, [srcreg+64]	;; I1
	vmovapd	zmm14, [srcreg+d4+64]	;; I5
	vmovapd	zmm15, [srcreg+d2+64]	;; I3
	vmovapd	zmm16, [srcreg+d4+d2+64]	;; I7
 vmovapd zmm17, [screg+1*128]		;; sine for R3/I3 and R4/I4
 vmovapd zmm17, [screg+3*128+64]	;; cosine/sine for R2/I2
 vmovapd zmm17, [screg+2*128+64]	;; cosine/sine for R7/I7 and R6/I6
 vmovapd zmm17, [screg+0*128+64]	;; cosine/sine for R5/I5 and R8/I8
 vmovapd zmm17, [screg+1*128+64]	;; cosine/sine for R3/I3 and R4/I4
 vmovapd zmm17, [screg+3*128]		;; sine for R2/I2
 vmovapd zmm17, [screg+2*128]		;; sine for R7/I7 and R6/I6
 vmovapd zmm17, [screg+0*128]		;; sine for R5/I5 and R8/I8
	bump	screg, scinc
	zstore	[srcreg], zmm1			;; Store R1
	zstore	[srcreg+d1], zmm2		;; Store R2
	zstore	[srcreg+d2], zmm3		;; Store R3
	zstore	[srcreg+d2+d1], zmm4		;; Store R4
	zstore	[srcreg+d4], zmm5		;; Store R5
	zstore	[srcreg+d4+d1], zmm6		;; Store R6
	zstore	[srcreg+d4+d2], zmm7		;; Store R7
	zstore	[srcreg+d4+d2+d1], zmm8		;; Store R8
	zstore	[srcreg+64], zmm9		;; Store I1
	zstore	[srcreg+d1+64], zmm10		;; Store I2
	zstore	[srcreg+d2+64], zmm11		;; Store I3
	zstore	[srcreg+d2+d1+64], zmm12	;; Store I4
	zstore	[srcreg+d4+64], zmm13		;; Store I5
	zstore	[srcreg+d4+d1+64], zmm14	;; Store I6
	zstore	[srcreg+d4+d2+64], zmm15	;; Store I7
	zstore	[srcreg+d4+d2+d1+64], zmm16	;; Store I8
	bump	srcreg, srcinc
	ENDM

do_nothing_plus_L1prefetch MACRO srcreg,srcinc,d1,d2,d4,screg,scinc,maxrpt,L1pt,L1pd
	L1prefetchw srcreg+L1pd, L1pt
	L1prefetchw srcreg+64+L1pd, L1pt
	L1prefetchw srcreg+d1+L1pd, L1pt
	L1prefetchw srcreg+d1+64+L1pd, L1pt
	L1prefetchw srcreg+d2+L1pd, L1pt
	L1prefetchw srcreg+d2+64+L1pd, L1pt
	L1prefetchw srcreg+d2+d1+L1pd, L1pt
	L1prefetchw srcreg+d2+d1+64+L1pd, L1pt
	L1prefetchw srcreg+d4+L1pd, L1pt
	L1prefetchw srcreg+d4+64+L1pd, L1pt
	L1prefetchw srcreg+d4+d1+L1pd, L1pt
	L1prefetchw srcreg+d4+d1+64+L1pd, L1pt
	L1prefetchw srcreg+d4+d2+L1pd, L1pt
	L1prefetchw srcreg+d4+d2+64+L1pd, L1pt
	L1prefetchw srcreg+d4+d2+d1+L1pd, L1pt
	L1prefetchw srcreg+d4+d2+d1+64+L1pd, L1pt
	vmovapd	zmm1, [srcreg+d1]	;; R2
	vmovapd	zmm2, [srcreg+d4+d1]	;; R6
	vmovapd	zmm3, [srcreg+d2+d1+64]	;; I4
	vmovapd	zmm4, [srcreg+d4+d2+d1+64] ;; I8
	vmovapd	zmm5, [srcreg+d2+d1]	;; R4
	vmovapd	zmm6, [srcreg+d4+d2+d1]	;; R8
	vmovapd	zmm7, [srcreg+d1+64]	;; I2
	vmovapd	zmm8, [srcreg+d4+d1+64]	;; I6
	vmovapd	zmm9, [srcreg]		;; R1
	vmovapd	zmm10, [srcreg+d4]	;; R5
	vmovapd	zmm11, [srcreg+d2]	;; R3
	vmovapd	zmm12, [srcreg+d4+d2]	;; R7
	vmovapd	zmm13, [srcreg+64]	;; I1
	vmovapd	zmm14, [srcreg+d4+64]	;; I5
	vmovapd	zmm15, [srcreg+d2+64]	;; I3
	vmovapd	zmm16, [srcreg+d4+d2+64]	;; I7
	zstore	[srcreg], zmm1			;; Store R1
	zstore	[srcreg+d1], zmm2		;; Store R2
	zstore	[srcreg+d2], zmm3		;; Store R3
	zstore	[srcreg+d2+d1], zmm4		;; Store R4
	zstore	[srcreg+d4], zmm5		;; Store R5
	zstore	[srcreg+d4+d1], zmm6		;; Store R6
	zstore	[srcreg+d4+d2], zmm7		;; Store R7
	zstore	[srcreg+d4+d2+d1], zmm8		;; Store R8
	zstore	[srcreg+64], zmm9		;; Store I1
	zstore	[srcreg+d1+64], zmm10		;; Store I2
	zstore	[srcreg+d2+64], zmm11		;; Store I3
	zstore	[srcreg+d2+d1+64], zmm12	;; Store I4
	zstore	[srcreg+d4+64], zmm13		;; Store I5
	zstore	[srcreg+d4+d1+64], zmm14	;; Store I6
	zstore	[srcreg+d4+d2+64], zmm15	;; Store I7
	zstore	[srcreg+d4+d2+d1+64], zmm16	;; Store I8
	bump	srcreg, srcinc
	ENDM

; this code reads/writes 64MB (1M cache lines) in contiguous blocks.  Timings are done
; on 4 memory sizes.  4KB will operate on the L1 cache only, 128KB will operate on the
; L2 cache only, 2MB will operate on the L3 cache and 32MB will operate on main memory.

av512case0:	read8	4096, 16384	; Read 4KB
		jmp	exit
av512case1:	read8	128*1024, 512	; Read 128KB
		jmp	exit
av512case2:	read8	2048*1024, 32	; Read 2MB
		jmp	exit
av512case3:	read8	32768*1024, 2	; Read 32MB
		jmp	exit

av512case4:	write8	4096, 16384	; Write 4KB
		jmp	exit
av512case5:	write8	128*1024, 512	; Write 128KB
		jmp	exit
av512case6:	write8	2048*1024, 32	; Write 2MB
		jmp	exit
av512case7:	write8	32768*1024, 2	; Write 32MB
		jmp	exit

av512case8:	readwrite8 4096, 16384	; Read/write 4KB
		jmp	exit
av512case9:	readwrite8 128*1024, 512 ; Read/write 128KB
		jmp	exit
av512case10:	readwrite8 2048*1024, 32 ; Read/write 2MB
		jmp	exit
av512case11:	readwrite8 32768*1024, 2 ; Read/write 32MB
		jmp	exit

		;; Simulate various ways to implement a 9600K FFT, pass 1 = 1280, pass 2 = 7680.
		;; Our AVX512 implementation is much, much slower in pass 1 with it's scatter reads
		;; and scatter writes vs. pass 2 with its sequential reads and sequential writes.

		;; Optimal, but not possible - sequential read and write.
av512case12:	sequential_readwrite8 1280*7680*8, 5 ; Read/write 78MB
		jmp	exit
		;; Scatter read and scatter write -- simulates pass 1 with clm = 2
av512case13:	scatter_readwrite8 1280*7680*8, 5 ; Read/write 78MB
		jmp	exit
		;; Sequential read and scatter write -- simulates rewritten pass 1 & 2
av512case14:	sequential_read_scatter_write8 1280*7680*8, 5 ; Read/write 78MB
		jmp	exit
		;; For completeness, scatter read and sequential write
av512case15:	scatter_read_sequential_write8 1280*7680*8, 5 ; Read/write 78MB
		jmp	exit


	zloop_init 32			;; Dummy call to zloop_init

	avx512_case_num = 16

;;16
	avx512mac 1, 10000, 0, tenclocks		; Macro that should take exactly 10 clocks (for calibration purposes)
	avx512mac 1, 10000, 0, hundredclocks		; Macro that should take exactly 100 clocks (for calibration purposes)
	avx512mac 1, 10000, 0, hundredclocks512		; compare 100 clocks of zmm vs ymm
	avx512mac 1, 10000, 0, do5 do5 FMApenalty_test1
	avx512mac 1, 10000, 0, do2 do10 FMApenalty_test2
	avx512mac 1, 10000, 0, do2 do10 FMApenalty_test3
	avx512mac 1, 10000, 0, do2 do10 FMApenalty_test4
	avx512mac 1, 10000, 0, do2 do10 FMApenalty_test5

;;24							; Macros that do nothing but load and store -- tests best case scenario
	avx512mac 10240, 16384, 0, do10 do_nothing rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1	;; Fastest possible L1 cache macro
	avx512mac 10240, 262144, 0, do10 do_nothing rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1	;; Fastest possible L2 cache macro
	avx512mac 10240, 4194304, 0, do10 do_nothing rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1	;; Fastest possible L3 cache macro
	avx512mac 10240, 32768*1024, 0, do10 do_nothing rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1	;; Fastest possible memory macro
	avx512mac 10240, 16384, 0, do10 do_nothing_plus_sincos rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1
	avx512mac 10240, 262144, 0, do10 do_nothing_plus_sincos rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1
	avx512mac 10240, 262144, 0, do10 do_nothing_plus_L1prefetch rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1, L1PREFETCH_ALL, 1*8*128
	avx512mac 10240, 262144, 0, do10 do_nothing_plus_L1prefetch rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1, L1PREFETCH_ALL, 2*8*128
	avx512mac 10240, 262144, 0, do10 do_nothing_plus_L1prefetch rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1, L1PREFETCH_ALL, 3*8*128
	avx512mac 10240, 262144, 0, do10 do_nothing_plus_L1prefetch rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1, L1PREFETCH_ALL, 4*8*128

;;34
	avx512mac 1024, 16384, 0, zr8_eight_complex_djbfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1
	avx512mac 1024, 262144, 0, zr8_eight_complex_djbfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1
	avx512mac 1024, 262144, 0, zr8_eight_complex_djbfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1, L1PREFETCH_ALL, 1*8*128
	avx512mac 1024, 262144, 0, zr8_eight_complex_djbfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1, L1PREFETCH_ALL, 2*8*128
	avx512mac 1024, 16384, 0, zr8_eight_complex_djbunfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1
	avx512mac 1024, 262144, 0, zr8_eight_complex_djbunfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD4, 1

;;40
	avx512mac 640, 16384, 0, zr5_five_complex_djbfft rsi, 5*128, 128, rdi, ZMM_SCD2, 1
	avx512mac 640, 262144, 0, zr5_five_complex_djbfft rsi, 5*128, 128, rdi, ZMM_SCD2, 1
	avx512mac 640, 16384, 0, zr5_five_complex_djbunfft rsi, 5*128, 128, rdi, ZMM_SCD2, 1
	avx512mac 640, 262144, 0, zr5_five_complex_djbunfft rsi, 5*128, 128, rdi, ZMM_SCD2, 1
	avx512mac 768, 16384, 0, zr6_six_complex_djbfft rsi, 6*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 768, 262144, 0, zr6_six_complex_djbfft rsi, 6*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 768, 16384, 0, zr6_six_complex_djbunfft rsi, 6*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 768, 262144, 0, zr6_six_complex_djbunfft rsi, 6*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 896, 16384, 0, zr7_seven_complex_djbfft rsi, 7*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 896, 262144, 0, zr7_seven_complex_djbfft rsi, 7*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 896, 16384, 0, zr7_seven_complex_djbunfft rsi, 7*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 896, 262144, 0, zr7_seven_complex_djbunfft rsi, 7*128, 128, rdi, ZMM_SCD3, 1
	avx512mac 1536, 16384, 0, zr12_twelve_complex_djbfft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD6, 1
	avx512mac 1536, 262144, 0, zr12_twelve_complex_djbfft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD6, 1
	avx512mac 1536, 16384, 0, zr12_twelve_complex_djbunfft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD6, 1
	avx512mac 1536, 262144, 0, zr12_twelve_complex_djbunfft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD6, 1
	avx512mac 2048, 16384, 0, zr16_sixteen_complex_djbfft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD8, 1
	avx512mac 2048, 262144, 0, zr16_sixteen_complex_djbfft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD8, 1
	avx512mac 2048, 16384, 0, zr16_sixteen_complex_djbunfft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD8, 1
	avx512mac 2048, 262144, 0, zr16_sixteen_complex_djbunfft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD8, 1

;;60
	avx512mac 640, 16384, 0, zr5_csc_ten_reals_fft rsi, 5*128, 128, rdi, ZMM_SCD5, 1
	avx512mac 640, 262144, 0, zr5_csc_ten_reals_fft rsi, 5*128, 128, rdi, ZMM_SCD5, 1
	avx512mac 640, 16384, 0, zr5_csc_ten_reals_unfft rsi, 5*128, 128, rdi, ZMM_SCD5, 1
	avx512mac 640, 262144, 0, zr5_csc_ten_reals_unfft rsi, 5*128, 128, rdi, ZMM_SCD5, 1
	avx512mac 768, 16384, 0, zr6_csc_twelve_reals_fft rsi, 6*128, 128, rdi, ZMM_SCD6, 1
	avx512mac 768, 262144, 0, zr6_csc_twelve_reals_fft rsi, 6*128, 128, rdi, ZMM_SCD6, 1
	avx512mac 768, 16384, 0, zr6_csc_twelve_reals_unfft rsi, 6*128, 128, rdi, ZMM_SCD6, 1
	avx512mac 768, 262144, 0, zr6_csc_twelve_reals_unfft rsi, 6*128, 128, rdi, ZMM_SCD6, 1
	avx512mac 896, 16384, 0, zr7_csc_fourteen_reals_fft rsi, 7*128, 128, rdi, ZMM_SCD7, 1
	avx512mac 896, 262144, 0, zr7_csc_fourteen_reals_fft rsi, 7*128, 128, rdi, ZMM_SCD7, 1
	avx512mac 896, 16384, 0, zr7_csc_fourteen_reals_unfft rsi, 7*128, 128, rdi, ZMM_SCD7, 1
	avx512mac 896, 262144, 0, zr7_csc_fourteen_reals_unfft rsi, 7*128, 128, rdi, ZMM_SCD7, 1
	avx512mac 1024, 16384, 0, zr8_csc_sixteen_reals_fft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD8, 1
	avx512mac 1024, 262144, 0, zr8_csc_sixteen_reals_fft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD8, 1
	avx512mac 1024, 16384, 0, zr8_csc_sixteen_reals_unfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD8, 1
	avx512mac 1024, 262144, 0, zr8_csc_sixteen_reals_unfft rsi, 8*128, 128, 2*128, 4*128, rdi, ZMM_SCD8, 1
	avx512mac 1536, 16384, 0, zr12_csc_twentyfour_reals_fft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD12, 1
	avx512mac 1536, 262144, 0, zr12_csc_twentyfour_reals_fft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD12, 1
	avx512mac 1536, 16384, 0, zr12_csc_twentyfour_reals_unfft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD12, 1
	avx512mac 1536, 262144, 0, zr12_csc_twentyfour_reals_unfft rsi, 12*128, 128, 2*128, 4*128, rdi, ZMM_SCD12, 1
	avx512mac 2048, 16384, 0, zr16_csc_thirtytwo_reals_fft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD16, 1
	avx512mac 2048, 262144, 0, zr16_csc_thirtytwo_reals_fft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD16, 1
	avx512mac 2048, 16384, 0, zr16_csc_thirtytwo_reals_unfft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD16, 1
	avx512mac 2048, 262144, 0, zr16_csc_thirtytwo_reals_unfft rsi, 16*128, 128, 2*128, 4*128, 8*128, rdi, ZMM_SCD16, 1

;;84
	zr64_sixtyfour_complex_fft_final_preload
	avx512mac 1024, 16384, 0, zr64_sixtyfour_complex_fft_final rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	avx512mac 1024, 262144, 0, zr64_sixtyfour_complex_fft_final rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	zr64_sixtyfour_complex_with_square_preload
	avx512mac 1024, 16384, 0, zr64_sixtyfour_complex_with_square rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	avx512mac 1024, 262144, 0, zr64_sixtyfour_complex_with_square rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	zr64_sixtyfour_complex_with_mult_preload
	avx512mac 1024, 16384, 0, zr64_sixtyfour_complex_with_mult rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	avx512mac 1024, 262144, 0, zr64_sixtyfour_complex_with_mult rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	zr64f_sixtyfour_complex_with_mulf_preload
	avx512mac 1024, 16384, 0, zr64f_sixtyfour_complex_with_mulf rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1
	avx512mac 1024, 262144, 0, zr64f_sixtyfour_complex_with_mulf rsi, 8*128, 1*128, 2*128, 4*128, rdi, 0, 1

;;92
	zr8_csc_wpn_eight_complex_first_djbfft_preload
	avx512mac 1024, 16384, 0, zr8_csc_wpn_eight_complex_first_djbfft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1024, 262144, 0, zr8_csc_wpn_eight_complex_first_djbfft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	zr8_csc_wpn_eight_complex_last_djbunfft_preload
	avx512mac 1024, 16384, 0, zr8_csc_wpn_eight_complex_last_djbunfft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1024, 262144, 0, zr8_csc_wpn_eight_complex_last_djbunfft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1

;;96
	zr8_wpn_sixteen_reals_first_fft_preload
	avx512mac 1024, 16384, 0, zr8_wpn_sixteen_reals_first_fft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1024, 262144, 0, zr8_wpn_sixteen_reals_first_fft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	zr8_wpn_sixteen_reals_last_unfft_preload
	avx512mac 1024, 16384, 0, zr8_wpn_sixteen_reals_last_unfft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1024, 262144, 0, zr8_wpn_sixteen_reals_last_unfft rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1

;;100
	zr12_csc_wpn_twelve_complex_first_djbfft_preload
	avx512mac 1536, 16384, 0, zr12_csc_wpn_twelve_complex_first_djbfft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1536, 262144, 0, zr12_csc_wpn_twelve_complex_first_djbfft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	zr12_csc_wpn_twelve_complex_last_djbunfft_preload
	avx512mac 1536, 16384, 0, zr12_csc_wpn_twelve_complex_last_djbunfft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1536, 262144, 0, zr12_csc_wpn_twelve_complex_last_djbunfft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1

;;104
	zr12_wpn_twentyfour_reals_first_fft_preload
	avx512mac 1536, 16384, 0, zr12_wpn_twentyfour_reals_first_fft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1536, 262144, 0, zr12_wpn_twentyfour_reals_first_fft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	zr12_wpn_twentyfour_reals_last_unfft_preload
	avx512mac 1536, 16384, 0, zr12_wpn_twentyfour_reals_last_unfft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1
	avx512mac 1536, 262144, 0, zr12_wpn_twentyfour_reals_last_unfft rsi, 12*128, 128, 2*128, 4*128, r15, 0, rdi, 16, rdx, 0, 1

;;108
	avx512mac 1024, 8192, 0, zr8_rsc_wpn_sgreg_eight_complex_fft8 rsi, 1024, 128, 256, 512, rdx, 1024, r8, r9, r10, rdi, ZMM_SCD8, 1
	avx512mac 1024, 262144, 0, zr8_rsc_wpn_sgreg_eight_complex_fft8 rsi, 1024, 128, 256, 512, rdx, 1024, r8, r9, r10, rdi, ZMM_SCD8, 1
	avx512mac 1024, 8192, 0, zr8_rsc_wpn_sgreg_eight_complex_unfft8 rdx, 1024, r8, r9, r10, rsi, 1024, 128, 256, 512, rdi, ZMM_SCD8, 1
	avx512mac 1024, 262144, 0, zr8_rsc_wpn_sgreg_eight_complex_unfft8 rdx, 1024, r8, r9, r10, rsi, 1024, 128, 256, 512, rdi, ZMM_SCD8, 1
	avx512mac 1024, 8192, 0, zr8_rsc_wpn_sgreg_2sc_sixteen_reals_fft8 rsi, 1024, 128, 256, 512, rdx, 1024, r8, r9, r10, rdi, 0, rdi, 0, rdi, ZMM_SCD8, 1
	avx512mac 1024, 262144, 0, zr8_rsc_wpn_sgreg_2sc_sixteen_reals_fft8 rsi, 1024, 128, 256, 512, rdx, 1024, r8, r9, r10, rdi, 0, rdi, 0, rdi, ZMM_SCD8, 1
	avx512mac 1024, 8192, 0, zr8_rsc_wpn_sgreg_2sc_sixteen_reals_unfft8 rdx, 1024, r8, r9, r10, rsi, 1024, 128, 256, 512, rdi, 0, rdi, 0, rdi, ZMM_SCD8, 1
	avx512mac 1024, 262144, 0, zr8_rsc_wpn_sgreg_2sc_sixteen_reals_unfft8 rdx, 1024, r8, r9, r10, rsi, 1024, 128, 256, 512, rdi, 0, rdi, 0, rdi, ZMM_SCD8, 1

;;116
	zsf_onepass_real_fft_wrapper_preload
	avx512mac 1024, 16384, 0, zsf_onepass_real_fft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 14, rdx, 0, 1
	avx512mac 1024, 262144, 0, zsf_onepass_real_fft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 14, rdx, 0, 1
	zs_onepass_real_unfft_wrapper_preload
	avx512mac 1024, 16384, 0, zs_onepass_real_unfft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 14, 1
	avx512mac 1024, 262144, 0, zs_onepass_real_unfft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 14, 1
	zsf_onepass_complex_fft_wrapper_preload
	avx512mac 1024, 16384, 0, zsf_onepass_complex_fft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 14, rdx, 0, 1
	avx512mac 1024, 262144, 0, zsf_onepass_complex_fft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, rdi, 14, rdx, 0, 1
	zs_onepass_complex_unfft_wrapper_preload
	avx512mac 1024, 16384, 0, zs_onepass_complex_unfft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, r15, 0, rdi, 14, 1
	avx512mac 1024, 262144, 0, zs_onepass_complex_unfft_wrapper rsi, 8*128, 128, 2*128, 4*128, r15, 0, r15, 0, rdi, 14, 1

;;124
	znorm_wpn_preload noexec, noexec, noexec, noexec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn noexec, noexec, noexec, noexec
	avx512znormmac 1536, 262144, 0, znorm_wpn noexec, noexec, noexec, noexec
	znorm_wpn_preload noexec, noexec, exec, noexec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn noexec, noexec, exec, noexec
	avx512znormmac 1536, 262144, 0, znorm_wpn noexec, noexec, exec, noexec
	znorm_wpn_preload noexec, noexec, noexec, exec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn noexec, noexec, noexec, exec
	avx512znormmac 1536, 262144, 0, znorm_wpn noexec, noexec, noexec, exec
	znorm_wpn_preload noexec, noexec, exec, exec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn noexec, noexec, exec, exec
	avx512znormmac 1536, 262144, 0, znorm_wpn noexec, noexec, exec, exec

;;132
	znorm_wpn_preload exec, noexec, noexec, noexec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn exec, noexec, noexec, noexec
	avx512znormmac 1536, 262144, 0, znorm_wpn exec, noexec, noexec, noexec
	znorm_wpn_preload exec, noexec, exec, noexec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn exec, noexec, exec, noexec
	avx512znormmac 1536, 262144, 0, znorm_wpn exec, noexec, exec, noexec
	znorm_wpn_preload exec, noexec, noexec, exec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn exec, noexec, noexec, exec
	avx512znormmac 1536, 262144, 0, znorm_wpn exec, noexec, noexec, exec
	znorm_wpn_preload exec, noexec, exec, exec	;; ttp, zero, echk, const
	avx512znormmac 1536, 16384, 0, znorm_wpn exec, noexec, exec, exec
	avx512znormmac 1536, 262144, 0, znorm_wpn exec, noexec, exec, exec

ENDIF

; Exit the timing code

exit:	ad_epilog 0,0,rbx,rbp,rsi,rdi,r8,r9,r10,r11,r12,r13,r14,r15,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

gwtimeit ENDP

_TEXT ENDS
END
