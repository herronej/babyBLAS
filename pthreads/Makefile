include ../Makefile.inc

all: lib

vvm.o: vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mmm.o: mmm.c
	$(CC) $(CFLAGS) -c mmm.c

mvv.o: mvv.c
	$(CC) $(CFLAGS) -c mvv.c

dot.o: dot.c
	$(CC) $(CFLAGS) -c dot.c

dls.o: dls.c
	$(CC) $(CFLAGS) -c dls.c

ils.o: ils.c
	$(CC) $(CFLAGS) -c ils.c

lib: vvm.o mmm.o mvv.o dot.o dls.o ils.o
	ar -rc libpthreads_bblas.a vvm.o mmm.o mvv.o dot.o dls.o ils.o
	ranlib libpthreads_bblas.a

clean:
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
