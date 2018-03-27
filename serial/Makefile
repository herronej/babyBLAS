include ../Makefile.inc

all: lib

vvm.o: vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mmm.o: mmm.c
	$(CC) $(CFLAGS) -C mmm.c

lib: vvm.o mmm.o
	ar -rc libserial_bblas.a vvm.o mmm.o
	ranlib libserial_bblas.a

clean:
	rm *.o
	touch *.c

pristine:
	make clean 
	rm *.a
