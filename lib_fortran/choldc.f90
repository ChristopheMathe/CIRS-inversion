!**********************************************************************!
!*** INVERT MATRIX BY CHOLESKY DECOMPOSITION                        ***!
!*** -------------------------------------------------------------- ***!
!*** BASE ON : NUMERICAL RECIPES IN FORTRAN                         ***!
!***           SECOND EDITION                                       ***!
!*** -------------------------------------------------------------- ***!
!*** Given a positive-definite symmetric matrix A(1:n,1:n), this    ***!
!*** routine computes the invert matrix A by Cholesky decompostion  ***!
!*** method's : A == LL^T.                                          ***!
!*** First stage   : construct the Cholesky factor L, returned in   ***!
!***                 the lower triangle of A, except for its        ***!
!***                 diagonal elements which are returned in P(1:n).***!
!*** Second stage  : set the upper triangle of A to 0               ***!
!*** Third stage   : compute inversion for the Cholesky Factor L    ***!
!***                 (in the lower triangle of matrix A)            ***!
!*** Fourth stage  : compute inverse of matrix A by :               ***!
!***                 (A)^{-1} = (L^{-1})^T * L^{-1}                 ***!
!*** -------------------------------------------------------------- ***!
!*** Author : Christophe Mathé                                      ***!
!*** Date : 19/11/2016                                              ***!
!**********************************************************************!


SUBROUTINE choldc(a,n)

implicit none

integer, intent(in) :: n                ! Dimension of square matrix A

double precision, dimension(n,n) :: a   ! Square matrix A
double precision, dimension(n)   :: p   ! Diagonal elements of matrix A

integer :: i                            ! Loop on rows
integer :: j                            ! Loop on columns
integer :: k                            ! Loop to sum elements

double precision :: summ = 0.           ! Buffer of sum elements
double precision :: time_begin = 0.     ! Begin of the subroutine
double precision :: time_end = 0.       ! End of the subroutine
double precision :: time_elapsed = 0.   ! Time spend of the subroutine (seconds)

call cpu_time(time_begin)

!======================================================================!
!=== FIRST STAGE : CONSTRUCT CHOLESKY FACTOR L                      ===!
!======================================================================!

do i = 1, n
    do j = i, n 
        summ = a(i,j)
        do k =i-1, 1, -1
            summ = summ - a(i,k)*a(j,k)
        end do
        if (i == j) then
            if(summ.le.0) stop 'choldc failed'
            p(i) = dsqrt(summ)
        else
            a(j,i) = summ/p(i)
        end if
    end do
end do

!======================================================================!
!=== SECOND STAGE : SET UPPER TRIANGLE TO 0                         ===!
!======================================================================!

do i = 1, n
    do j = i+1, n
        a(i,j) = 0.0d0
    end do
end do

!======================================================================!
!=== THIRD STAGE : COMPUTE INVERSE OF L                             ===!
!======================================================================!

do i = 1, n
    A(i,i) = 1./P(i)
    do j = i+1, n
        summ = 0.0d0
        do k = i, j-1
            summ = summ - A(j,k)*A(k,i)
        end do
        A(j,i) = summ / P(j)
    end do
end do

!======================================================================!
!=== FOURTH STAGE : COMPUTE INVERSE OF A                            ===!
!======================================================================!

A = matmul(transpose(A),A)

call cpu_time(time_end)

time_elapsed = time_end - time_begin
write(*,*)"inverse matrix takes :",time_elapsed," s"


return
END SUBROUTINE choldc
