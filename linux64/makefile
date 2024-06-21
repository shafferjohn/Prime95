# Makefile for Linux 64-bit mprime
#
# NOTES: libcurl built by downloading latest sources and:
#	./configure --without-ssl --without-libssh2 --without-libidn --disable-ldap --disable-ldaps --without-gssapi --enable-ipv6 --without-krb4
#	make
#	make install
#
# Hwloc can be downloaded from https://www.open-mpi.org/projects/hwloc/
# Configure hwloc with --enable-static --disable-shared --disable-io
#              or with --enable-static --disable-shared --disable-pci --disable-cuda --disable-gl --disable-opencl --disable-rsmi --disable-nvml
#
#  Ugh, different linux variants require different makefiles.
#  The current makefile is for CentOS 7.  We prefer to link against older Linux versions because linking on the latest, greatest version
#  creates an mprime executable that will not run on older Linux versions because of glibc incompatibilites.
#
#  Some linux versions require some of the variations below:
#	"export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig"
#	CFLAGS = -I.. -I../gwnum -DX86_64 -O2 -Wno-unused-result
#	LFLAGS = -Wl,-M
#	LIBS = ../gwnum/gwnum.a -lm -lpthread -Wl,-Bstatic $(shell pkg-config --static --libs libhwloc) -Wl,-Bstatic $(shell pkg-config --static --libs libcurl) -lstdc++ -Wl,-Bdynamic -ldl

CC = gcc
CFLAGS = -I.. -I../gwnum -I/usr/local/include -std=gnu99 -DX86_64 -Wno-unused-result -O2 -DPORT=8

CPP = g++
CPPFLAGS = -I.. -I../gwnum -I../../boost_1_76_0 -DX86_64 -std=c++11 -O2 -DPORT=8

LFLAGS = -Wl,-M -Wl,-L/usr/local/lib
#LIBS = ../gwnum/gwnum.a ../gwnum/polymult.a -lm -lpthread -Wl,-Bstatic -lhwloc -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -lz -lxml2 -ldl -lgmp
LIBS = ../gwnum/gwnum.a ../gwnum/polymult.a -lm -lpthread -Wl,-Bstatic -lhwloc -Wl,-Bstatic -lcurl -Wl,-Bdynamic -lrt -lstdc++ -ldl -lgmp

FACTOROBJ = factor64.o
LINUXOBJS = prime.o menu.o cJSON.o ecm.o exponentiate.o pair.o pm1prob.o
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

cJSON.o:
	$(CC) $(CFLAGS) -c ../cJSON.c

ecm.o:
	$(CPP) $(CPPFLAGS) -c ../ecm.cpp

exponentiate.o:
	$(CC) $(CFLAGS) -c ../exponentiate.c

pair.o:
	$(CPP) $(CPPFLAGS) -c ../pair.cpp

pm1prob.o:
	$(CC) $(CFLAGS) -c ../pm1prob.c

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CPP) $(CPPFLAGS) -c $<
