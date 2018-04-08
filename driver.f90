program driver 

integer :: NDIM

real (kind=8) :: wall_start, wall_end
real (kind=8) :: cpu_start, cpu_end
real (kind=8) :: trace, rval


integer :: startval, stopval, stepval, nthreads
real (kind=8) :: walltime
real (kind=8) :: cputime 
external walltime, cputime

character (len=8) :: carg1, carg2, carg3, carg4

real (kind=8), dimension(:), allocatable :: veca, vecb, vec1, vec2, vec3, rvec,vecx
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrixc, matrix90
logical :: DIAG_DOMINANT, SPARSE_MATRIX
real (kind=8) :: residual

DIAG_DOMINANT = .false.
SPARSE_MATRIX = .false.
!This portion of code is ONLY used for verifying the accuracy of the code using
!the matrix and matrix inverse stored on the class website.

!Download the files from theochem using curl (don't store these on anvil!)
call system("curl -s -o matrixa.dat --url http://theochem.mercer.edu/csc435/data/matrixa.dat")
call system("curl -s -o matrixb.dat --url http://theochem.mercer.edu/csc435/data/matrixb.dat")

NDIM = 100  ! The test files are 100x100 double precision matrix and its inverse
nthreads = 2
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)
allocate ( matrix90(NDIM, NDIM), stat=ierr)
allocate ( vec1(NDIM), stat=ierr)
allocate ( vec2(NDIM), stat=ierr)
allocate ( vec3(3), stat=ierr)
allocate ( rvec(3), stat=ierr)
allocate ( vecx(NDIM), stat=ierr)
open (unit=5,file="matrixa.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixa(j,i)
  enddo
enddo
close(5)
open (unit=5,file="matrixb.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixb(j,i)
  enddo
enddo
close(5)

do i = 1, NDIM
        vec1(i) = 1.0*i
        vec2(i) = 1.0/vec1(i)
        !print*, vec1(i), vec2(i)
enddo

matrix90(1,1) = cos(90.0)
matrix90(1,2) = -sin(90.0)
matrix90(1,3) = 0.0
matrix90(2,1) = sin(90.0)
matrix90(2,2) = cos(90.0)
matrix90(2,3) = 0.0
matrix90(3,1) = 0.0
matrix90(3,2) = 0.0
matrix90(3,3) = 1.0

vec3(1) = 1.0
vec3(2) = 0.0
vec3(3) = 0.0

do i = 1, 3
        rvec = 0.0
enddo

! Delete the files from disk
call system("rm matrixa.dat matrixb.dat")

wall_start = walltime()
cpu_start = cputime()

print*, "calling mmm"
!call mmm(nthreads, NDIM, matrixa, matrixb, matrixc)
call dls(nthreads, NDIM, matrixa, vecb, vecx)

cpu_end = cputime()
wall_end = walltime()

trace = 0.0;

do i=1, NDIM 
     trace = trace + matrixc(i,i)
enddo

print*, "trace: ", trace

! Calculate megaflops based on CPU time and Walltime

mflops  = 2*dble(NDIM)**3/ (cpu_end-cpu_start) / 1.0e6
mflops2 = 2*dble(NDIM)**3/ (wall_end-wall_start)/ 1.0e6
 
print *, NDIM, trace, cpu_end-cpu_start, wall_end-wall_start,  mflops, mflops2


print*, "dot: "
rval = 0.0
call dot(nthreads, NDIM, vec1, vec2, rval)
print *, rval

print*, "mvv"
call mvv(3, matrix90, vec3, rvec)

do i = 1, 3 
        print*, rvec(i)
enddo

! Free the memory that was allocated based on which version of the program was
! run.

deallocate(matrixa)
deallocate(matrixb)
deallocate(matrixc)
deallocate(matrix90)
deallocate(vec1)
deallocate(vec2)
deallocate(vec3)
deallocate(rvec)
deallocate(vecx)

end program driver 
