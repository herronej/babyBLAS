include ../Makefile.inc

all: lib

mmm.o: mmm.c
	$(CC) $(CFLAGS) -c mmm.c

lib: mmm.o
	ar -rc libpthreads_bblas.a mmm.o
	ranlib libpthreads_bblas.a

clean:
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
