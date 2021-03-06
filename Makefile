include Makefile.inc

all : driver serial openmp pthreads lbstime

driver: driver.o serial lbstime openmp pthreads
	$(F90) driver.o -o driver $(LDLIBS) $(SYSLIBS)

driver.o: driver.f90
	$(F90) $(FFLAGS) driver.f90 -c  

serial: 
	cd serial && $(MAKE)

lbstime: 
	cd lbstime && $(MAKE)

pthreads:
	cd pthreads && $(MAKE) 

openmp:
	cd openmp && $(MAKE) 

clean:
	cd serial && $(MAKE) clean
	cd lbstime && $(MAKE) clean
	cd pthreads && $(MAKE) clean
	cd openmp && $(MAKE) clean
	rm *.o
	touch *.f90

pristine:
	cd serial && $(MAKE) pristine 
	cd lbstime && $(MAKE) pristine
	cd pthreads && $(MAKE) pristine
	cd openmp && $(MAKE) pristine
	rm *.o	
	rm driver
	touch *.f90

#This next target get "made" every time
.PHONY: serial lbstime openmp pthreads
