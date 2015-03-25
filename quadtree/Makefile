INCLUDE = -I/usr/local/include -I/usr/include -I/usr/X11/include -I/p/include -I/usr/include/X11
LIBDIR = -L/usr/local/lib -L/usr/lib -L/usr/X11/lib -L/p/lib -L/usr/lib/x86_64-linux-gnu

LIBS = -lGL -lglut -lz -lX11 -ldl -lXext

CFLAGS = -g

CC = g++

.cpp.o: 
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $*.o $<

ALL.O = quadtree.o

quadtree: $(ALL.O)
	$(CC) $(CFLAGS) -o $@ $(ALL.O) $(INCLUDE) $(LIBDIR) $(LIBS)

clean:
	rm -f *.o quadtree

