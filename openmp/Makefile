include ../Makefile.inc

# OpenMP library targets

all: lib

mmm.o: mmm.c
	$(CC) $(CFLAGS) -fopenmp -c mmm.c

lib: mmm.o
	ar -rc libopenmp_bblas.a mmm.o
	ranlib libopenmp_bblas.a

clean:
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
