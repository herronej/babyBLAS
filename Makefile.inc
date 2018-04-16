F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O2 -DLS_TEST #-DACCURACY_TEST -DLS_TEST -DITERATIVE
CFLAGS = -g -fopenmp -O2
LDLIBS = -L./pthreads -lpthreads_bblas -L./lbstime -llbstime 
SYSLIBS = -lgomp -lpthread

