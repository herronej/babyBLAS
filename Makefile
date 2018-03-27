# Main makefile

include Makefile.inc

all: driver serial pthreads openmp lbstime

driver: driver.o serial lbstime pthreads openmp
	$(F90) driver.o -o driver $(MYLIBS) $(SYSLIBS)

driver.o: driver.f90
	$(F90) $(FFLAGS) driver.f90 -c

serial:
	cd serial && $(MAKE)

pthreads:
	cd pthreads && $(MAKE)

openmp:
	cd openmp && $(MAKE)

lbstime:
	cd lbstime && $(MAKE)

clean:
	cd serial && $(MAKE) clean
	cd pthreads && $(MAKE) clean
	cd openmp && $(MAKE) clean
	cd lbstime && $(MAKE) clean
	rm *.o
	touch *.f90

pristine:
	cd serial && $(MAKE) pristine
	cd pthreads && $(MAKE) pristine
	cd openmp && $(MAKE) pristine
	rm *.o
	rm rm driver
	touch *.f90

.PHONY:
	serial pthreads pthreads openmp lbstime
