PF = -w


UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
LIBS = -lGL -lglut -lz -lX11 -ldl -lXext
LIBDIR = -L/usr/local/lib -L/usr/lib -L/usr/X11/lib -L/p/lib -L/usr/lib/x86_64-linux-gnu
INCLUDE = -I/usr/local/include -I/usr/include -I/usr/X11/include -I/p/include -I/usr/include/X11 -I../common
endif
ifeq ($(UNAME), Darwin)
LIBS = -framework OpenGL -framework GLUT
LIBDIR = -L/usr/lib -L/usr/X11/lib
INCLUDE = -I/usr/include -I../common
endif

CFLAGS = -g $(PF)

CC = g++

.cpp.o: 
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $*.o $<

ALL.O = kdtree.o

COMMON.H = ../common/bitmap.h

COMMON.O = ../common/bitmap.o

kdtree: $(ALL.O) $(COMMON.H)
	$(CC) $(CFLAGS) -o $@ $(ALL.O) $(COMMON.O) $(INCLUDE) $(LIBDIR) $(LIBS)

clean:
	rm -f *.o kdtree

