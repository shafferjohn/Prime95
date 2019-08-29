REM
REM  Strip all debug info from the object file.  We do this because
REM  object files built by MASM 10 are rejected by MSVC 2005's linker
REM  when building an executable with debug info.
REM
c:\objconv\objconv -felf %1 dummy.o
c:\objconv\objconv -fcoff dummy.o %1
del dummy.o
