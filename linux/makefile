# Makefile for Linux 32-bit mprime
#
# NOTES: libcurl built by downloading latest sources and:
#	./configure --without-ssl --without-libssh2 --without-libidn --disable-ldap --disable-ldaps --without-gssapi --enable-ipv6 --without-krb4
#	make
#	make install
# After linking, check to see if too many dynamic libraries are required by doing:
#	ldd mprime
#
# Hwloc can be downloaded from https://www.open-mpi.org/projects/hwloc/
# Configure hwloc with --enable-static --disable-shared --disable-pci --disable-cuda --disable-gl --disable-opencl
#
#  Ugh, different linux variants require different makefiles.
#  The current makefile is for CentOS 6.  We prefer to link against
#  older Linux versions because linking on the latest, greatest version
#  will create an mprime executable that will not run on older
#  Linux versions because of glibc incompatibilites.
#
#  Some linux versions require some of the variations below:
#	"export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig"
#	CFLAGS = -I.. -I../gwnum -O2 -march=i486 -malign-double -Wno-unused-result
#	LFLAGS = -Wl,-M
#	LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic $(shell pkg-config --static --libs libhwloc) -Wl,-Bstatic $(shell pkg-config --static --libs libcurl) -lstdc++ -Wl,-Bdynamic -ldl

CC = gcc
CFLAGS = -I.. -I../gwnum -Wno-unused-result -march=i486 -malign-double -O2

CPP = g++
CPPFLAGS = -I.. -I../gwnum -Wno-unused-result -march=i486 -malign-double -O2

LFLAGS = -Wl,-M -Wl,-L/usr/local/lib
LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic -lhwloc -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -lz -lxml2 -ldl -lgmp
#LIBS = ../gwnum/gwnum.a ../gwnum/gwnum.ld -lm -lpthread -Wl,-Bstatic -lhwloc -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -ldl -lgmp

FACTOROBJ = factor32.o
LINUXOBJS = prime.o menu.o
EXE      = mprime

#########################################################################

$(EXE): security $(LINUXOBJS) $(FACTOROBJ)
	$(CC) $(LFLAGS) -o $(EXE) $(LINUXOBJS) $(FACTOROBJ) $(LIBS)

security:
	[ ! -e ../security.h ] && touch ../security.h || true
	[ ! -e ../security.c ] && touch ../security.c || true
	[ ! -e ../secure5.c ] && touch ../secure5.c || true

clean:
	rm -f $(EXE) $(EXE2) $(LINUXOBJS)

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
