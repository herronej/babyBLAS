F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O2
CFLAGS = -O2
LDLIBS = -L./serial -lserial_bblas -L./lbstime -llbstime 
