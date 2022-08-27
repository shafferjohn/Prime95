# Makefile for Linux and FreeBSD 32-bit gwnum library
#

CC = gcc
CFLAGS = -I.. -I../sqlite-amalgamation-3180000 -I/usr/local/include -std=gnu99 -Wno-unused-result -march=i486 -malign-double -O2
#CFLAGS = -g -I.. -I../sqlite-amalgamation-3180000 -I/usr/local/include -std=gnu99 -Wno-unused-result -march=i486 -malign-double -DGDEBUG -O2

CPP = g++
CPPFLAGS = -I.. -I../qd -std=c++11 -O2 -march=i486 -malign-double
#CPPFLAGS = -g -I.. -I../qd -std=c++11 -O2 -march=i486 -malign-double -DGDEBUG

AR = ar

LINUXOBJS = cpuid.o gwnum.o gwtables.o gwthread.o gwini.o gwbench.o gwutil.o gwdbldbl.o giants.o radix.o ecmstag1.o
POLYOBJS = polymult.o polymult_sse2.o polymult_avx.o polymult_fma.o

LIB = gwnum.a
POLYLIB = polymult.a

#########################################################################

all:	$(LIB) $(POLYLIB)

$(LIB): $(LINUXOBJS)
	cp linux/gwnum.a .
	$(AR) -rs $(LIB) $(LINUXOBJS)

$(POLYLIB): $(POLYOBJS)
	$(AR) -rs $(POLYLIB) $(POLYOBJS)

clean:
	rm -f $(LINUXOBJS) $(POLYOBJS)

distclean: clean
	rm -f $(LIB) $(POLYLIB)

polymult_sse2.o: polymult.c polymult.h
	$(CC) $(CFLAGS) -Dsse2 -msse2 -c -o polymult_sse2.o polymult.c

polymult_avx.o: polymult.c polymult.h
	$(CC) $(CFLAGS) -DAVX -mavx -c -o polymult_avx.o polymult.c

polymult_fma.o: polymult.c polymult.h
	$(CC) $(CFLAGS) -DFMA -mfma -c -o polymult_fma.o polymult.c

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
