F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O2 #-DLS_TEST -DITERATIVE #-DACCURACY_TEST -DLS_TEST -DITERATIVE
CFLAGS = -O2 -g -fopenmp -std=c99
LDLIBS = -L./openmp -lopenmp_bblas -L./lbstime -llbstime 
SYSLIBS = -lgomp -lpthread

