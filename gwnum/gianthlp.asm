; Copyright 1995-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This file implements helper routines for the ECM code
;

	TITLE   gianthlp

IFNDEF X86_64
	.386
	.MODEL	FLAT
ENDIF

INCLUDE unravel.mac

;
; Global variables
;

_DATA SEGMENT
TWOPOW32	DD	4294967296.0		; 2^32
BIGVAL		DD	5F000000h		; 2^63
_DATA	ENDS

_TEXT	SEGMENT

; addhlp (res, carry, val)
;	Add val to carry,res
; Windows 32-bit (_addhlp)
; Linux 32-bit (addhlp)
;	Parameter res = [esp+4]
;	Parameter carry = [esp+8]
;	Parameter val = [esp+12]
; Windows 64-bit (addhlp) - leaf routine, no unwind info necessary
;	Parameter res = rcx
;	Parameter carry = rdx
;	Parameter val = r8
; Linux 64-bit (addhlp)
;	Parameter res = rdi
;	Parameter carry = rsi
;	Parameter val = rdx

PROCL	addhlp
IFNDEF X86_64
	mov	eax, [esp+12]		; Val
	mov	edx, [esp+4]		; Load res pointer
	add	DWORD PTR [edx], eax	; Add two 32-bit integers
	mov	ecx, [esp+8]		; Load carry pointer
	adc	DWORD PTR [ecx], 0	; Add the carry
ENDIF
IFDEF WINDOWS64
	add	DWORD PTR [rcx], r8d	; Add two 32-bit integers
	adc	DWORD PTR [rdx], 0	; Add the carry
ENDIF
IFDEF LINUX64
	add	DWORD PTR [rdi], edx	; Add two 32-bit integers
	adc	DWORD PTR [rsi], 0	; Add the carry
ENDIF
	ret
addhlp	ENDP

; subhlp (res, carry, val)
;	Subtract val from carry,res
; Windows 32-bit (_subhlp)
; Linux 32-bit (subhlp)
;	Parameter res = [esp+4]
;	Parameter carry = [esp+8]
;	Parameter val = [esp+12]
; Windows 64-bit (subhlp) - leaf routine, no unwind info necessary
;	Parameter res = rcx
;	Parameter carry = rdx
;	Parameter val = r8
; Linux 64-bit (subhlp)
;	Parameter res = rdi
;	Parameter carry = rsi
;	Parameter val = rdx

PROCL	subhlp
IFNDEF X86_64
	mov	eax, [esp+12]		; Val
	mov	edx, [esp+4]		; Load res pointer
	sub	DWORD PTR [edx], eax	; Subtract two 32-bit integers
	mov	ecx, [esp+8]		; Load carry pointer
	sbb	DWORD PTR [ecx], 0	; Subtract the carry
ENDIF
IFDEF WINDOWS64
	sub	DWORD PTR [rcx], r8d	; Subtract two 32-bit integers
	sbb	DWORD PTR [rdx], 0	; Subtract the carry
ENDIF
IFDEF LINUX64
	sub	DWORD PTR [rdi], edx	; Subtract two 32-bit integers
	sbb	DWORD PTR [rsi], 0	; Subtract the carry
ENDIF
	ret
subhlp	ENDP

; muladdhlp (res, carryl, carryh, val1, val2)
;	Multiply val1 and val2 adding result to carryh,carryl,res
; Windows 32-bit (_muladdhlp)
; Linux 32-bit (muladdhlp)
;	Parameter res = [esp+4]
;	Parameter carryl = [esp+8]
;	Parameter carryh = [esp+12]
;	Parameter val1 = [esp+16]
;	Parameter val2 = [esp+20]
; Windows 64-bit (muladdhlp) - leaf routine, no unwind info necessary
;	Parameter res = rcx
;	Parameter carryl = rdx
;	Parameter carryh = r8
;	Parameter val1 = r9
;	Parameter val2 = [rsp+40]
; Linux 64-bit (muladdhlp)
;	Parameter res = rdi
;	Parameter carryl = rsi
;	Parameter carryh = rdx
;	Parameter val1 = rcx
;	Parameter val2 = r8

PROCL	muladdhlp
IFNDEF X86_64
	mov	eax, [esp+16]		; Val1
	mul	DWORD PTR [esp+20]	; Multiply by val2
	mov	ecx, [esp+4]		; Load res pointer
	add	DWORD PTR [ecx], eax	; Add result to 3 word accumulator
	mov	ecx, [esp+8]		; Load carryl pointer
	adc	DWORD PTR [ecx], edx
	mov	eax, [esp+12]		; Load carryh pointer
	adc	DWORD PTR [eax], 0	; Add the carry
ENDIF
IFDEF WINDOWS64
	mov	eax, [rsp+40]		; Load val2
	mov	r10, rdx		; Save carryl ptr
	mul	r9d			; Mul two 32-bit integers
	add	DWORD PTR [rcx], eax	; Add result to 3 word accumulator
	adc	DWORD PTR [r10], edx
	adc	DWORD PTR [r8], 0
ENDIF
IFDEF LINUX64
	mov	eax, r8d		; Load val2
	mov	r8, rdx			; Save carryh pointer
	mul	ecx			; Mul two 32-bit integers
	add	DWORD PTR [rdi], eax	; Add result to 3 word accumulator
	adc	DWORD PTR [rsi], edx
	adc	DWORD PTR [r8], 0
ENDIF
	ret
muladdhlp ENDP


; muladd2hlp (res, carryl, carryh, val1, val2)
;	Multiply val1 and val2 adding twice the result to carryh,carryl,res
; Windows 32-bit (_muladd2hlp)
; Linux 32-bit (muladd2hlp)
;	Parameter res = [esp+4]
;	Parameter carryl = [esp+8]
;	Parameter carryh = [esp+12]
;	Parameter val1 = [esp+16]
;	Parameter val2 = [esp+20]
; Windows 64-bit (muladd2hlp) - leaf routine, no unwind info necessary
;	Parameter res = rcx
;	Parameter carryl = rdx
;	Parameter carryh = r8
;	Parameter val1 = r9
;	Parameter val2 = [rsp+40]
; Linux 64-bit (muladd2hlp)
;	Parameter res = rdi
;	Parameter carryl = rsi
;	Parameter carryh = rdx
;	Parameter val1 = rcx
;	Parameter val2 = r8

PROCL	muladd2hlp
IFNDEF X86_64
	mov	eax, [esp+16]		; Val1
	mul	DWORD PTR [esp+20]	; Multiply by val2
	add	eax, eax		; Double the result
	adc	edx, edx
	mov	ecx, [esp+12]		; Load carryh pointer
	adc	DWORD PTR [ecx], 0	; Add the carry
	mov	ecx, [esp+4]		; Load res pointer
	add	DWORD PTR [ecx], eax	; Add result to 3 word accumulator
	mov	ecx, [esp+8]		; Load carryl pointer
	adc	DWORD PTR [ecx], edx
	mov	eax, [esp+12]		; Load carryh pointer
	adc	DWORD PTR [eax], 0	; Add the carry
ENDIF
IFDEF WINDOWS64
	mov	eax, [rsp+40]		; Load val2
	mov	r10, rdx		; Save carryl ptr
	mul	r9d			; Mul two 32-bit integers
	add	eax, eax		; Double the result
	adc	edx, edx
	adc	DWORD PTR [r8], 0
	add	DWORD PTR [rcx], eax	; Add result to 3 word accumulator
	adc	DWORD PTR [r10], edx
	adc	DWORD PTR [r8], 0
ENDIF
IFDEF LINUX64
	mov	eax, r8d		; Load val2
	mov	r8, rdx			; Save carryh pointer
	mul	ecx			; Mul two 32-bit integers
	add	eax, eax		; Double the result
	adc	edx, edx
	adc	DWORD PTR [r8], 0
	add	DWORD PTR [rdi], eax	; Add result to 3 word accumulator
	adc	DWORD PTR [rsi], edx
	adc	DWORD PTR [r8], 0
ENDIF
	ret
muladd2hlp ENDP


; mulsubhlp (res, carryl, carryh, val1, val2)
;	Multiply val1 and val2 subtracting result from carryh,carryl,res
; Windows 32-bit (_mulsubhlp)
; Linux 32-bit (mulsubhlp)
;	Parameter res = [esp+4]
;	Parameter carryl = [esp+8]
;	Parameter carryh = [esp+12]
;	Parameter val1 = [esp+16]
;	Parameter val2 = [esp+20]
; Windows 64-bit (mulsubhlp) - leaf routine, no unwind info necessary
;	Parameter res = rcx
;	Parameter carryl = rdx
;	Parameter carryh = r8
;	Parameter val1 = r9
;	Parameter val2 = [rsp+40]
; Linux 64-bit (mulsubhlp)
;	Parameter res = rdi
;	Parameter carryl = rsi
;	Parameter carryh = rdx
;	Parameter val1 = rcx
;	Parameter val2 = r8

PROCL	mulsubhlp
IFNDEF X86_64
	mov	eax, [esp+16]		; Val1
	mul	DWORD PTR [esp+20]	; Multiply by val2
	mov	ecx, [esp+4]		; Load res pointer
	sub	DWORD PTR [ecx], eax	; Sub result from 3 word accumulator
	mov	ecx, [esp+8]		; Load carryl pointer
	sbb	DWORD PTR [ecx], edx
	mov	eax, [esp+12]		; Load carryh pointer
	sbb	DWORD PTR [eax], 0
ENDIF
IFDEF WINDOWS64
	mov	eax, [rsp+40]		; Load val2
	mov	r10, rdx		; Save carryl ptr
	mul	r9d			; Mul two 32-bit integers
	sub	DWORD PTR [rcx], eax	; Sub result from 3 word accumulator
	sbb	DWORD PTR [r10], edx
	sbb	DWORD PTR [r8], 0
ENDIF
IFDEF LINUX64
	mov	eax, r8d		; Load val2
	mov	r8, rdx			; Save carryh pointer
	mul	ecx			; Mul two 32-bit integers
	sub	DWORD PTR [rdi], eax	; Sub result from 3 word accumulator
	sbb	DWORD PTR [rsi], edx
	sbb	DWORD PTR [r8], 0
ENDIF
	ret
mulsubhlp ENDP

; gcdhlp (ulen, udata, vlen, vdata, &struct)
;
; Routine to help in computing extended GCD quickly
;
; Do several single-precision steps for the extended GCD code
; U is larger than V.  Both are the same length or U is one larger
; than V.  Returns a structure containing A,B,C,D,ODD where A,B,C,D
; are as defined in Knuth vol 2. description of
; extended GCD for large numbers.  This was implemented in assembly
; language to better deal with 64-bit integers.

; Windows 32-bit (_gcdhlp)
; Linux 32-bit (gcdhlp)
;	Parameter ulen = [esp+4]
;	Parameter udata = [esp+8]
;	Parameter vlen = [esp+12]
;	Parameter vdata = [esp+16]
;	Parameter struct_ptr = [esp+20]
; Windows 64-bit (gcdhlp)
;	Parameter ulen = rcx
;	Parameter udata = rdx
;	Parameter vlen = r8
;	Parameter vdata = r9
;	Parameter struct_ptr = [rsp+40]
; Linux 64-bit (gcdhlp)
;	Parameter ulen = rdi
;	Parameter udata = rsi
;	Parameter vlen = rdx
;	Parameter vdata = rcx
;	Parameter struct_ptr = r8

PROCFL	gcdhlp
IFNDEF X86_64
	ah_prolog 4,16,0,rbx,rbp,rdi,rsi

; Load up to 64 bits of U and V
;	U will be in edx:eax with bits shifted in from edi
;	V will be in esi:ebp with bits shifted in from ebx

	mov	edi, [esp+push_amt+8]	; Giant U
	mov	ebx, [esp+push_amt+4]	; U->sign
	mov	ecx, [esp+push_amt+16]	; Giant V
	mov	ebp, [esp+push_amt+12]	; V->sign

	mov	edx, [edi-4][ebx*4]	; U[Ulen-1]

	xor	esi, esi
	cmp	ebx, ebp
	jne	short noload
	mov	esi, [ecx-4][ebx*4]	; V[Ulen-1]

noload:	cmp	ebx, 1			; Are there more words to shift
	jg	short multi		; bits from?
	xor	eax, eax		; No, zero out MSWs
	xor	ebp, ebp
	xchg	eax, edx
	xchg	esi, ebp
	jmp	short noshft

multi:	mov	eax, [edi-8][ebx*4]	; U[Ulen-2]
	mov	ebp, [ecx-8][ebx*4]	; V[Ulen-2]
	
	cmp	ebx, 2			; Are there more words to shift
	je	short noshft		; bits from?

	mov	edi, [edi-12][ebx*4]	; U[Ulen-3]
	mov	ebx, [ecx-12][ebx*4]	; V[Ulen-3]

	bsr	ecx, edx		; Count bits to shift U
	xor	ecx, 31			; Turn bit # into a shift count

	shld	edx, eax, cl		; Shift U
	shld	eax, edi, cl
	shld	esi, ebp, cl		; Shift V
	shld	ebp, ebx, cl

; Init extended GCD information
;	A = 1;
;	C = 0;
;	B = 0;
;	D = 1;
;	ODD = 0;

noshft:	fld1				; A
	fldz				; C
	fldz				; B
	fld1				; D

; Load the values into the FPU
;	U is in edx:eax
;	V is in esi:ebp

	mov	DWORD PTR [esp], edx
	mov	DWORD PTR [esp+4], 0
	mov	DWORD PTR [esp+8], esi
	mov	DWORD PTR [esp+12], 0
	fild	QWORD PTR [esp]		; uhi
	fild	QWORD PTR [esp+8]	; vhi, uhi
	mov	DWORD PTR [esp], eax
	mov	DWORD PTR [esp+8], ebp
	fmul	TWOPOW32
	fxch	st(1)			; uhi, vhi
	fmul	TWOPOW32
	fild	QWORD PTR [esp]		; ulo, uhi, vhi
	fild	QWORD PTR [esp+8]	; vlo, ulo, uhi, vhi
	faddp	st(3), st		; ulo, uhi, v
	faddp	st(1), st		; u, v

; Turn on truncating mode, 64-bit precision

	fstcw	WORD PTR [esp]
	or	WORD PTR [esp], 0F00h
	fldcw	WORD PTR [esp]

; Load structure pointer to return data

	mov	esi, DWORD PTR [esp+push_amt+20]

; Check if we are doing an exact GCD

	cmp	DWORD PTR [esp+push_amt+4], 2	; Test Ulen <= 2
	jle	simple

; Do as many single precision operations as we can
;	FPU contains U, V, D, B, C, A
; As Knuth suggests:
;	Compute (U-B)/(V+D), the smaller quotient
;	Compute (U+A)/(V-C), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

dloop:	fld	st(0)			; U, U, V, D, B, C, A
	fsub	st(0), st(4)		; U-B, U, V, D, B, C, A
	fld	st(2)			; V, U-B, U, V, D, B, C, A
	fadd	st(0), st(4)		; V+D, U-B, U, V, D, B, C, A
	fdivp	st(1), st(0)		; Q, U, V, D, B, C, A
	fadd	BIGVAL
	fld	st(2)			; V, Q, U, V, D, B, C, A
	fsub	st(0), st(6)		; V-C, Q, U, V, D, B, C, A
	fxch	st(1)			; Q, V-C, U, V, D, B, C, A
	fsub	BIGVAL
	fst	QWORD PTR [esp]
	fmul	st(0), st(1)		; Q*(V-C), V-C, U, V, D, B, C, A
	faddp	st(1), st(0)		; (Q+1)*(V-C), U, V, D, B, C, A
	fld	QWORD PTR [esp]		; Q, (Q+1)*(V-C), U, V, D, B, C, A
	fmul	st(0), st(4)		; Q*D, (Q+1)*(V-C), U, V, D, B, C, A
	fxch	st(1)			; (Q+1)*(V-C), Q*D, U, V, D, B, C, A
	fsub	st(0), st(2)		; (Q+1)*(V-C)-U, Q*D, U, V, D, B, C, A
	fxch	st(1)			; Q*D, (Q+1)*(V-C)-U, U, V, D, B, C, A
	fadd	st(0), st(5)		; B+Q*D,(Q+1)*(V-C)-U, U, V, D, B, C, A
	fxch	st(1)			; (Q+1)*(V-C)-U,B+Q*D, U, V, D, B, C, A
	fcomp	st(7)			; B+Q*D, U, V, D, B, C, A
	fstsw	ax			; Copy comparison results
	and	eax, 4100h		; Isolate C0 & C3 bits
	jnz	lpdone			; Break if less than or equal
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	jz	lpdone			; Break if not less than
	fxch	st(4)			; oldB, U, V, newB, newD, C, A
	fcomp	st(1)			; U, V, newB, newD, C, A
	fld	QWORD PTR [esp]		; Q, U, V, newB, newD, C, A
	fmul	st, st(5)		; Q*C, U, V, newB, newD, C, A
	fld	QWORD PTR [esp]		; Q, Q*C, U, V, newB, newD, C, A
	fmul	st, st(3)		; Q*V, Q*C, U, V, newB, newD, C, A
	fxch	st(1)			; Q*C, Q*V, U, V, newB, newD, C, A
	faddp	st(7), st		; Q*V, U, V, newB, newD, newA, newC
	fsubp	st(1), st		; newV, newU, newB, newD, newA, newC

;	Compute (U-A)/(V+C), the smaller quotient
;	Compute (U+B)/(V-D), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	fld	st(1)			; U, V, U, B, D, A, C
	fsub	st(0), st(5)		; U-A, V, U, B, D, A, C
	fld	st(1)			; V, U-A, V, U, B, D, A, C
	fadd	st(0), st(7)		; V+C, U-A, V, U, B, D, A, C
	fdivp	st(1), st(0)		; Q, V, U, B, D, A, C
	fadd	BIGVAL
	fld	st(1)			; V, Q, V, U, B, D, A, C
	fsub	st(0), st(5)		; V-D, Q, V, U, B, D, A, C
	fxch	st(1)			; Q, V-D, V, U, B, D, A, C
	fsub	BIGVAL
	fst	QWORD PTR [esp]
	fmul	st(0), st(1)		; Q*(V-D), V-D, V, U, B, D, A, C
	faddp	st(1), st(0)		; (Q+1)*(V-D), V, U, B, D, A, C
	fld	QWORD PTR [esp]		; Q, (Q+1)*(V-D), V, U, B, D, A, C
	fmul	st(0), st(5)		; Q*D, (Q+1)*(V-D), V, U, B, D, A, C
	fxch	st(1)			; (Q+1)*(V-D), Q*D, V, U, B, D, A, C
	fsub	st(0), st(3)		; (Q+1)*(V-D)-U, Q*D, V, U, B, D, A, C
	fxch	st(1)			; Q*D, (Q+1)*(V-D)-U, V, U, B, D, A, C
	fadd	st(0), st(4)		; B+Q*D,(Q+1)*(V-D)-U, V, U, B, D, A, C
	fxch	st(1)			; (Q+1)*(V-D)-U,B+Q*D, V, U, B, D, A, C
	fcomp	st(4)			; B+Q*D, V, U, B, D, A, C
	fstsw	ax			; Copy comparison results
	and	eax, 4100h		; Isolate C0 & C3 bits
	jnz	lpdone1			; Break if less than or equal
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	jz	lpdone1			; Break if not less than
	fxch	st(3)			; oldB, V, U, newD, newB, A, C
	fcomp	st(1)			; V, U, newD, newB, A, C
	fld	QWORD PTR [esp]		; Q, V, U, newD, newB, A, C
	fmul	st, st(6)		; Q*C, V, U, newD, newB, A, C
	fld	QWORD PTR [esp]		; Q, Q*C, V, U, newD, newB, A, C
	fmul	st, st(2)		; Q*V, Q*C, V, U, newD, newB, A, C
	fxch	st(1)			; Q*C, Q*V, V, U, newD, newB, A, C
	faddp	st(6), st		; Q*V, V, U, newD, newB, newC, newA
	fsubp	st(2), st		; newU, newV, newD, newB, newC, newA
	jmp	dloop

; The single precision case:
;	Compute U/V, the quotient
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

simple:	xor	ecx, ecx		; Clear EGCD_ODD flag
sloop:	fld	st(0)			; U, U, V, D, B, C, A
	fdiv	st(0), st(2)		; Q, U, V, D, B, C, A
	fadd	BIGVAL
	fsub	BIGVAL
	fst	QWORD PTR [esp]
	fmul	st(0), st(3)		; Q*D, U, V, D, B, C, A
	fadd	st(0), st(4)		; B+Q*D, U, V, D, B, C, A
	fcom	TWOPOW32		; Will new D be too large?
	fstsw	ax			; Copy comparison results
	and	eax, 0100h		; Isolate C0 bit
	jz	short lpdone2		; Break if not less than
	fxch	st(3)			; newB, U, V, newD, oldB, C, A
	fxch	st(4)			; oldB, U, V, newD, newB, C, A
	fcomp	st(1)			; U, V, newD, newB, C, A
	fld	QWORD PTR [esp]		; Q, U, V, newD, newB, C, A
	fmul	st, st(5)		; Q*C, U, V, newD, newB, C, A
	fld	QWORD PTR [esp]		; Q, Q*C, U, V, newD, newB, C, A
	fmul	st, st(3)		; Q*V, Q*C, U, V, newD, newB, C, A
	fxch	st(1)			; Q*C, Q*V, U, V, newD, newB, C, A
	faddp	st(7), st		; Q*V, U, V, newD, newB, newA, newC
	fxch	st(5)
	fxch	st(6)
	fxch	st(5)			; Q*V, U, V, newD, newB, newC, newA
	fsubp	st(1), st		; newV, newU, newD, newB, newC, newA
	inc	ecx			; Toggle EGCD_ODD flag
	ftst
	fxch	st(1)			; newU, newV, newD, newB, newC, newA
	fstsw	ax			; Copy comparison results
	and	eax, 4000h		; Isolate C3 bit
	jz	short sloop		; Loop if V is not zero

; Copy extended GCD info to globals
; FPU has U, V, D, B, C, A

	fld	st(0)			; Push value for following code to pop
lpdone2:and	ecx, 1
	mov	DWORD PTR [esi+16], ecx	; Save ODD
	jmp	short lpd1

; Copy extended GCD info to globals
; FPU has B+Q*D, U, V, D, B, C, A

lpdone:	mov	DWORD PTR [esi+16], 0	; Save ODD
lpd1:	fcompp				; Pop 3 values
	fcomp	st(3)
	fistp	QWORD PTR [esp]		; D
	fistp	QWORD PTR [esp+8]	; B
	mov	eax, DWORD PTR [esp]
	mov	edx, DWORD PTR [esp+8]
	cmp	edx, 0			; Check for special case
	je	short lpd2
	mov	DWORD PTR [esi+12], eax	; Save D
	mov	DWORD PTR [esi+4], edx	; Save B
	fistp	QWORD PTR [esp]		; C
	fistp	QWORD PTR [esp+8]	; A
	mov	eax, DWORD PTR [esp]
	mov	edx, DWORD PTR [esp+8]
	mov	DWORD PTR [esi+8], eax	; Save C
	mov	DWORD PTR [esi], edx	; Save A
	jmp	egdone

; In this case the main loop couldn't determine even a single quotient.
; Return FALSE, so caller can handle by more sophisticated means.

lpd2:	fcompp				; Pop 2 values
	xor	eax, eax
	jmp	short egdone1

; Just like lpdone, but handles the second chunk of dloop code

lpdone1:fcompp				; Pop 3 values
	fcomp	st(3)
	fistp	QWORD PTR [esp]		; B
	fistp	QWORD PTR [esp+8]	; D
	mov	eax, DWORD PTR [esp]
	mov	edx, DWORD PTR [esp+8]
	mov	DWORD PTR [esi+4], eax	; Save B
	mov	DWORD PTR [esi+12], edx	; Save D
	fistp	QWORD PTR [esp]		; A
	fistp	QWORD PTR [esp+8]	; C
	mov	eax, DWORD PTR [esp]
	mov	edx, DWORD PTR [esp+8]
	mov	DWORD PTR [esi], eax	; Save A
	mov	DWORD PTR [esi+8], edx	; Save C
	mov	DWORD PTR [esi+16], 1	; Save ODD
	jmp	short egdone

; Restore FPU control word, return TRUE

egdone:	mov	eax, 1
egdone1:fstcw	WORD PTR [esp]
	and	WORD PTR [esp], 0F3FFh
	fldcw	WORD PTR [esp]
	ah_epilog 4,16,0,rbx,rbp,rdi,rsi
ENDIF

; The x86-64 version uses integer instructions

IFDEF X86_64
	ah_prolog 4,0,0,rbx,rbp,rsi,r12,r13,r14,r15
IFDEF WINDOWS64
	mov	r11, [rsp+push_amt+40]	; Load structure pointer
ENDIF
IFDEF LINUX64
	mov	r11, r8			; Load structure pointer
	mov	r9, rcx			; Copy parameters to Windows registers
	mov	r8, rdx
	mov	rdx, rsi
	mov	rcx, rdi
ENDIF

; Load up to 64 bits of U and V
;	U will be in r14
;	V will be in r15

	mov	r14d, [rdx][rcx*4-4]	; U[Ulen-1]
	xor	r15, r15		; Zero V in case U and V not same len
	cmp	rcx, r8			; Load top V word if U and V same len
	jne	short noload
	mov	r15d, [r9][rcx*4-4]	; V[Ulen-1]

	cmp	ecx, 1			; Are there more words to shift
	je	simple			; bits from?

noload:	shl	r14, 32
	shl	r15, 32
	mov	eax, [rdx][rcx*4-8]	; U[Ulen-2]
	mov	ebx, [r9][rcx*4-8]	; V[Ulen-2]
	add	r14, rax
	add	r15, rbx

	cmp	ecx, 2			; Are there more words to shift
	je	simple			; bits from?

	mov	eax, [rdx][rcx*4-12]	; U[Ulen-3]
	mov	ebx, [r9][rcx*4-12]	; V[Ulen-3]

	bsr	rcx, r14		; Count bits to shift U
	xor	cl, 63			; Turn bit # into a shift count

	shl	rax, 32
	shl	rbx, 32
	shld	r14, rax, cl		; Shift U
	shld	r15, rbx, cl		; Shift V

; Init extended GCD information

	mov	r8, 1			; A
	xor	r9, r9			; B
	xor	r10, r10		; C
	mov	r13, 1			; D
	xor	r12, r12		; ODD

; Do as many operations as we can constrained by a maximum 32-bit result
;	FPU contains U, V, D, B, C, A
; As Knuth suggests:
;	Compute (U-B)/(V+D), the smaller quotient
;	Compute (U+A)/(V-C), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	mov	rsi, 1			; Compute 2^32
	shl	rsi, 32
	shr	r14, 1			; Lose a bit so we can use signed
	shr	r15, 1			; multiply and divide instructions
dloop:	mov	rax, r14		; U
	sub	rax, r9			; U-B
	mov	rcx, r15		; V
	add	rcx, r13		; V+D
	xor	rdx, rdx
	idiv	rcx			; Q = (U-B)/(V+D)

	mov	rcx, r15		; V
	sub	rcx, r10		; V-C
	mov	rbp, rax		; Copy Q
	imul	rbp, rcx		; Q*(V-C)
	add	rbp, rcx		; (Q+1)*(V-C)
	lea	rbx, [r14+r8]		; U+A
	cmp	rbx, rbp		; Compare U+A to (Q+1)*(V-C)
	jae	lpdone			; Break if above or equal

	mov	rbx, r13		; D
	imul	rbx, rax		; Q*D
	add	rbx, r9			; B+Q*D (newD)
	cmp	rbx, rsi		; Will new D be too large?
	jae	lpdone			; Break if newD is 32-bits or more
	mov	r9, r13			; newB = D
	mov	r13, rbx		; newD = B+Q*D

	mov	rcx, r10		; C
	imul	rcx, rax		; Q*C
	add	rcx, r8			; A+Q*C (newC)
	mov	r8, r10			; newA = C
	mov	r10, rcx		; newC = A+Q*C

	imul	rax, r15		; Q*V
	sub	r14, rax		; U-Q*V (newV)
	xchg	r14, r15		; newU = V, newV = U-Q*V

	xor	r12, 1			; Flip ODD

;	Compute (U-A)/(V+C), the smaller quotient
;	Compute (U+B)/(V-D), the larger quotient, break if not equal
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	mov	rax, r14		; U
	sub	rax, r8			; U-A
	mov	rcx, r15		; V
	add	rcx, r10		; V+C
	xor	rdx, rdx
	idiv	rcx			; Q = (U-A)/(V+C)

	mov	rcx, r15		; V
	sub	rcx, r13		; V-D
	mov	rbp, rax		; Copy Q
	imul	rbp, rcx		; Q*(V-D)
	add	rbp, rcx		; (Q+1)*(V-D)
	lea	rbx, [r14+r9]		; U+B
	cmp	rbx, rbp		; Compare U+B to (Q+1)*(V-D)
	jae	lpdone			; Break if above or equal

	mov	rbx, r13		; D
	imul	rbx, rax		; Q*D
	add	rbx, r9			; B+Q*D (newD)
	cmp	rbx, rsi		; Will new D be too large?
	jae	lpdone			; Break if newD is 32-bits or more
	mov	r9, r13			; newB = D
	mov	r13, rbx		; newD = B+Q*D

	mov	rcx, r10		; C
	imul	rcx, rax		; Q*C
	add	rcx, r8			; A+Q*C (newC)
	mov	r8, r10			; newA = C
	mov	r10, rcx		; newC = A+Q*C

	imul	rax, r15		; Q*V
	sub	r14, rax		; U-Q*V (newV)
	xchg	r14, r15		; newU = V, newV = U-Q*V

	xor	r12, 1			; Flip ODD
	jmp	dloop

; The single precision case:

; Init extended GCD information

simple:	mov	r8, 1			; A
	xor	r9, r9			; B
	xor	r10, r10		; C
	mov	r13, 1			; D
	xor	r12, r12		; ODD

;	Compute U/V, the quotient
;	Set newB = D, newD = B+Q*D, break if newD won't fit in 32 bits
;	Set newA = C, newC = A+Q*C
;	Set newU = V, newV = U-Q*V

	mov	rsi, 1			; Compute 2^32
	shl	rsi, 32
sloop:	mov	rax, r14		; U
	xor	rdx, rdx
	div	r15			; Q = U/V
	mov	rbp, rax		; Copy Q

	mul	r13			; Q*D
	add	rax, r9			; B+Q*D
	cmp	rax, rsi		; Will new D be too large?
	jae	short lpdone		; Break if D >= 2^32
	mov	r9, r13			; newB = D
	mov	r13, rax		; newD = B+Q*D

	mov	rax, r10		; C
	mul	rbp			; Q*C
	add	rax, r8			; A+Q*C (newC)
	mov	r8, r10			; newA = C
	mov	r10, rax		; newC = A+Q*C

	mov	rax, r15		; V
	mul	rbp			; Q*V
	sub	r14, rax		; U-Q*V (newV)
	xchg	r14, r15		; newU = V, newV = U-Q*V

	xor	r12, 1			; Flip ODD

	and	r15, r15		; Loop if V is not zero
	jnz	short sloop

; Copy extended GCD info to globals

lpdone:	mov	DWORD PTR [r11], r8d	; Save A
	mov	DWORD PTR [r11+4], r9d	; Save B
	mov	DWORD PTR [r11+8], r10d	; Save C
	mov	DWORD PTR [r11+12], r13d ; Save D
	mov	DWORD PTR [r11+16], r12d ; Save ODD
	and	r9, r9			; Check for special case (B = 0)
	jnz	short egdone

; In this case the main loop couldn't determine even a single quotient.
; Return FALSE, so caller can handle by more sophisticated means.

	xor	rax, rax
	jmp	short egdone1

egdone:	mov	rax, 1
egdone1:ah_epilog 4,0,0,rbx,rbp,rsi,r12,r13,r14,r15
ENDIF

gcdhlp	ENDP

_TEXT	ENDS
END
