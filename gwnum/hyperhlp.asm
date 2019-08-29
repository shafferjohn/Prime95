; Copyright 2018 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Helper routines for hyperthread prefetching
;

	TITLE   hyperhlp

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac

_TEXT SEGMENT

;; Utility routines to help with hyperthread prefetching

; prefetchL2 (addr, cachelines)
;	Prefetch into L2 cache from addr for the specified number of 64-byte cache lines
; Windows 32-bit (_prefetchL2)
; Linux 32-bit (prefetchL2)
;	Parameter addr = [esp+4]
;	Parameter cachelines = [esp+8]
; Windows 64-bit (prefetchL2) - leaf routine, no unwind info necessary
;	Parameter addr = rcx
;	Parameter cachelines = rdx
; Linux 64-bit (prefetchL2)
;	Parameter addr = rdi
;	Parameter cachelines = rsi

PROCL	prefetchL2
IFNDEF X86_64
	mov	edi, [esp+4]		; Load addr
	mov	ecx, [esp+8]		; Load count
pfloop:	prefetcht1 [edi]		; Prefetch a 64-byte cache line
;;BUG?	pause				; Pause between prefetches
	bump	edi, 64			; Next cache line
	dec	ecx			; Test count
	jnz	short pfloop
ENDIF
IFDEF WINDOWS64
pfloop:	prefetcht1 [rcx]		; Prefetch a 64-byte cache line
;;	pause				; Pause between prefetches
	bump	rcx, 64			; Next cache line
	dec	rdx			; Test count
	jnz	short pfloop
ENDIF
IFDEF LINUX64
pfloop:	prefetcht1 [rdi]		; Prefetch a 64-byte cache line
	pause				; Pause between prefetches
	bump	rdi, 64			; Next cache line
	dec	rsi			; Test count
	jnz	short pfloop
ENDIF
	ret
prefetchL2 ENDP


; pause_for_count (count)
;	Execute "count" pause ops so that hyperthread can wait for some prefetching to do
;	without impacting the compute thread.
; Windows 32-bit (_pause_for_count)
; Linux 32-bit (pause_for_count)
;	Parameter count = [esp+4]
; Windows 64-bit (pause_for_count) - leaf routine, no unwind info necessary
;	Parameter count = rcx
; Linux 64-bit (pause_for_count)
;	Parameter count = rdi

PROCL	pause_for_count
IFNDEF X86_64
	mov	ecx, [esp+4]		; Load count
ploop:	pause				; Pause a few clocks
	dec	ecx			; Test count
	jnz	short ploop
ENDIF
IFDEF WINDOWS64
ploop:	pause				; Pause a few clocks
	dec	rcx			; Test count
	jnz	short ploop
ENDIF
IFDEF LINUX64
ploop:	pause				; Pause a few clocks
	dec	rdi			; Test count
	jnz	short ploop
ENDIF
	ret
pause_for_count ENDP

_TEXT	ENDS
END
