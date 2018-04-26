F90 = gfortran
CC  = gcc
FFLAGS = -cpp #-DLS_TEST -DITERATIVE #-DACCURACY_TEST -DLS_TEST -DITERATIVE
CFLAGS = -g -fopenmp -std=c99
LDLIBS = -L./openmp -lopenmp_bblas -L./lbstime -llbstime 
SYSLIBS = -lgomp -lpthread

