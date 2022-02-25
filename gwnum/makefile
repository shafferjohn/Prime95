# Makefile for Linux and FreeBSD 32-bit gwnum library
#

CC = gcc
CFLAGS = -I.. -I../sqlite-amalgamation-3180000 -I/usr/local/include -Wno-unused-result -march=i486 -malign-double -O2

CPP = g++
CPPFLAGS = -I.. -I../qd -O2 -march=i486 -malign-double

AR = ar

OBJS = cpuid.o gwnum.o gwtables.o gwthread.o gwini.o gwbench.o gwutil.o gwdbldbl.o giants.o radix.o ecmstag1.o

LIB = gwnum.a

#########################################################################

$(LIB): $(OBJS)
	cp linux/gwnum.a .
	$(AR) -rs $(LIB) $(OBJS)

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(LIB)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
