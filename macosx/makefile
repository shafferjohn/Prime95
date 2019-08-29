# Makefile for Mac OS X mprime
#
# Hwloc can be downloaded from https://www.open-mpi.org/projects/hwloc/
# Configure hwloc with --enable-static --disable-shared --disable-pci --disable-cuda --disable-gl
#

ENVP= MACOSX_DEPLOYMENT_TARGET=10.4

CC = gcc
CFLAGS = -I.. -I../gwnum -I../linux -O2 -m32 -march=i686 -DCOMMAND_LINE_MPRIME

CPP = g++
CPPFLAGS = -I.. -I../gwnum -I../linux -O2 -m32 -march=i686 -DCOMMAND_LINE_MPRIME

LFLAGS = -Wl,-no_pie -L/usr/local/lib
LIBS   = ../gwnum/release/gwnum.a -lm -lpthread -lhwloc -lcurl -framework IOKit -framework CoreFoundation -lstdc++ -lgmp

FACTOROBJ = ../prime95/macosx/factor32.o
OBJS = prime.o menu.o
EXE = mprime

#########################################################################

$(EXE): security $(OBJS) $(FACTOROBJ)
	$(ENVP) $(CC) $(LFLAGS) -o $(EXE) $(OBJS) $(FACTOROBJ) $(LIBS)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true
	[ ! -e ../secure5.c ] && touch ../secure5.c || true

clean:
	rm -f $(EXE) $(EXE2) $(OBJS)

.c.o:
	$(ENVP) $(CC) $(CFLAGS) -c $<

.cpp.o:
	$(ENVP) $(CPP) $(CPPFLAGS) -c $<
