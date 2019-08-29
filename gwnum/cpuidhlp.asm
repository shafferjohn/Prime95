; Copyright 1995-2012 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This file implements helper routines for the CPU identification code
; as well as other helper functions
;

	TITLE   cpuidhlp

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE unravel.mac

_TEXT SEGMENT

;
; Utility routine to initialize the FPU
;

; fpu_init ()
;	Initializes the x87 FPU.  Probably no longer necessary, but it was
;	necessary a long, long, time ago.
; Windows 32-bit (_fpu_init)
; Linux 32-bit (fpu_init)
;	No parameters or return value
; Windows 64-bit (fpu_init) - leaf routine, no unwind info necessary
; Linux 64-bit (fpu_init)
;	No parameters or return value

PROCL	fpu_init
IFNDEF X86_64
	fninit
ENDIF
	ret
fpu_init ENDP

;
; Utility routine to read the time stamp counter
;

; erdtsc (&hi, &lo)
;	Read the 64-bit timestamp counter and return it in 32-bit chunks
; Windows 32-bit (_erdtsc)
; Linux 32-bit (erdtsc)
;	Parameter hi = [esp+4]
;	Parameter lo = [esp+8]
; Windows 64-bit (erdtsc) - leaf routine, no unwind info necessary
;	Parameter hi = rcx
;	Parameter lo = rdx
; Linux 64-bit (erdtsc)
;	Parameter hi = rdi
;	Parameter lo = rsi

PROCL	erdtsc
IFNDEF X86_64
	rdtsc
	mov	ecx, [esp+4]		; Address to return info
	mov	[rcx], edx		; store high 32 bits
	mov	ecx, [esp+8]		; Address to return info
	mov	[rcx], eax		; store low 32 bits
ENDIF
IFDEF WINDOWS64
	mov	r8, rdx			; Copy return info address
	rdtsc
	mov	[rcx], edx		; store high 32 bits
	mov	[r8], eax		; store low 32 bits
ENDIF
IFDEF LINUX64
	rdtsc
	mov	[rdi], edx		; store high 32 bits
	mov	[rsi], eax		; store low 32 bits
ENDIF
	ret
erdtsc	ENDP


;
; Utility routine to see if CPUID is supported
; Returns non-zero if CPUID is supported
;

; ecpuidsupport ()
;	Return true if CPU supports CPUID
; Windows 32-bit (_ecpuidsupport)
; Linux 32-bit (ecpuidsupport)
;	Return value = eax
; Windows 64-bit (ecpuidsupport) - leaf routine, no unwind info necessary
; Linux 64-bit (ecpuidsupport)
;	Return value = rax

PROCL	ecpuidsupport
IFNDEF X86_64
        pushfd				; Get original EFLAGS
	pop	eax
	mov 	ecx, eax
        xor     eax, 200000h		; Flip ID bit in EFLAGS
        push    eax			; Save new EFLAGS value on stack
        popfd				; Replace current EFLAGS value
        pushfd				; Get new EFLAGS
        pop     eax			; Store new EFLAGS in EAX
        xor     eax, ecx		; Test toggled ID bit
ELSE
	mov	rax, 1			; 64-bit CPUs support CPUID
ENDIF
	ret				; Processor supports CPUID if eax != 0
ecpuidsupport ENDP

;
; Utility routine to execute CPUID
;

; ecpuid (resptr)
;	Call cpuid a fill a structure with the eax, ebx, ecx, edx values
; Windows 32-bit (_ecpuid)
; Linux 32-bit (ecpuid)
;	Parameter resptr = [esp+4]
; Windows 64-bit (ecpuid)
;	Parameter resptr = rcx
; Linux 64-bit (ecpuid)
;	Parameter res = rdi

PROCFL	ecpuid
	ah_prolog 1,0,0,rbx,rdi
IFNDEF X86_64
	mov	rdi, [esp+push_amt+4]	; Address of data struct to return info
ENDIF
IFDEF WINDOWS64
	mov	rdi, rcx		; Address of data struct to return info
ENDIF
	mov	eax, [rdi]		; EAX argument to CPUID
	mov	ecx, [rdi+8]		; ECX argument to CPUID
	cpuid				; Perform the CPUID
	mov	[rdi], eax		; Return the results
	mov	[rdi+4], ebx		; Return the results
	mov	[rdi+8], ecx		; Return the results
	mov	[rdi+12], edx		; Return the results
	ah_epilog 1,0,0,rbx,rdi
ecpuid	ENDP

;
; Utility routine to execute XGETBV
;

; exgetbv (resptr)
;	Call xgetbx a fill a structure with the eax, ebx, ecx, edx values
; Windows 32-bit (_ecpuid)
; Linux 32-bit (ecpuid)
;	Parameter resptr = [esp+4]
; Windows 64-bit (ecpuid)
;	Parameter resptr = rcx
; Linux 64-bit (ecpuid)
;	Parameter res = rdi

PROCFL	exgetbv
	ah_prolog 1,0,0,rbx,rdi
IFNDEF X86_64
	mov	rdi, [esp+push_amt+4]	; Address of data struct to return info
ENDIF
IFDEF WINDOWS64
	mov	rdi, rcx		; Address of data struct to return info
ENDIF
	mov	ecx, [rdi+8]		; ECX argument to XGETBV
	xgetbv				; Perform the XGETBV
	mov	[rdi], eax		; Return the results
	mov	[rdi+4], ebx		; Return the results
	mov	[rdi+8], ecx		; Return the results
	mov	[rdi+12], edx		; Return the results
	ah_epilog 1,0,0,rbx,rdi
exgetbv	ENDP

;
; Utility routines that hopefully take (roughly) 100,000 and 1,000,000 clock cycles
; This implementation is based on the ROR op-code which has latency of 1 or less
; on Intel and AMD architectures according to Agner Fog's documents.  The actual
; latency should be one because each ROR instruction is dependent on the previous one.
;
; The opcode chosen has a reciprocal throughput of exactly one on Intel so that when two
; hyperthreads are running on one physical core each thread will take twice as long to execute.
;

PROCFL	one_hundred_thousand_clocks_help
	ah_prolog 0,0,0,rdi,rsi
	mov	rax, 100000 / 100
biglp:	REPEAT 100
	ror	rcx, 1
	ENDM
	sub	rax, 1
	jnz	biglp
	ah_epilog 0,0,0,rdi,rsi
one_hundred_thousand_clocks_help ENDP

PROCFL	one_million_clocks_help
	ah_prolog 0,0,0,rdi,rsi
	mov	rax, 1000000 / 100
biglp2:	REPEAT 100
	ror	rcx, 1
	ENDM
	sub	rax, 1
	jnz	biglp2
	ah_epilog 0,0,0,rdi,rsi
one_million_clocks_help ENDP

_TEXT	ENDS
END
