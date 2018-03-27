program driver

integer :: NDIM

real(kind = 8) :: wall_start, wall_end
real(kind = 8) :: cpu_start, cpu_end
real(kind = 8) :: trace

integer :: startval, stopval, stepval, nthreads
real (kind=8) :: walltime
real (kind=8) :: cputime
external walltime, cputime

character(len=8) :: carg1, carg2, carg3, carg4

real (kind=8), dimension(:), allocatable :: veca, vecb
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrix

#ifdef ACCURACY_TEST

! download test files
call system("curl -s -o matrixa.dat --url
http://theochem.mercer.edu/csc435/data/matrixa.dat")
call system("curl -s -o matrixb.dat --url
http://theochem.mercer.edu/csc435/data/matrixa.dat")

NDIM = 100
nthreads = 2
allocate(matrixa(NDIM, NDIM), stat=ierr)
allocate(matrixb(NDIM, NDIM), stat=ierr)
allocate(matrixc(NDIM, NDIM), stat=ierr)
open (unit=5, file="matrixa.dat", status = "old")
do i = 1, NDIM
  do j = 1, NDIM
        read(5,*) matrix(j,i)
  enddo
enddo
close(5)
open(unit = 5, file = "matrixb.dat", status = "old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5*,) matrixa(j,i)
  enddo
enddo
close(5)

! delete matrix a and b from disk

#else

! read starting, stop, step value, number threads
call get_command_argument(1, carg1)
call get_command_argument(2, carg2)
call get_command_argument(3, carg3)
call get_command_argument(4, carg4)

! convert cargs to ints

read (carg1, '(i8)') startval
read (carg2, '(i8)') stopval
read (carg3, '(i8)') stepval
read (carg4, '(i8)') nthreads

! 

do iter = startval, stopval, stepval

NDIM = iter
allocate(veca(NDIM), stat=ierr)
allocate(vecb(NDIM), stat=ierr)
allocate(matrixa(NDIM, NDIM), stat=ierr)
allocate(matrixb(NDIM, NDIM), stat=ierr)
allocate(matrixc(NDIM, NDIM), stat=ierr)

! build veca and vecb
do i = 1, NDIM
        veca(i) = 1.0
        vecb(i) = 1.0 /sqrt(dble(NDIM))
enddo

matrixa = 0.0
matrixb = 0.0

call vvm(NDIM, veca, vecb, matrixa)
call vmm(NDIM, veca, vecb, matrixb)

#endif

wall_start = walltime()
cpu_start = cputime()

call mmm(nthreads, NDIM, matrixa, matrixb, matrixc)

cpu_end = cputime()
wall_end = walltime()

trace = 0.0

do i = 1, NDIM
        trace = trace + matrixc(i,i)
enddo

! Calculate Megaflops
mflops = 2*dble(NDIM)**3/(cpu_end-cpu_start)/1.0e6
mflops2 = 2*dble(NDIM)**3/(wall_end-wall_start)/1.0e6

print *, NDIM, trace, cpu_end-cpu_start, wall_end-wall_start, mflops, mflops2

deallocate(matrixa)
deallocate(matrixb)
deallocate(matrixc)

#ifndef ACCURACY_TEST
deallocate(veca)
deallocate(vecb)

enddo
#endif

end program driver


