!~ SUBROUTINE CONV_K
!~ =====================================================================
!~     R: Raw spectrum (no convolution)
!~     H: Hamming function (Res is first zero from center)
!~     S: Sinx/x (Res is first zero from center)
!~     T: Triangular function
!~     G: Gaussian function (Res is FWHM)
!~     L: Lorentz function (Res is FWHM)
!~     C: Square function (Res is the full width)
!~ =====================================================================

subroutine conv_k(f1,pas,n_sort,rad,iconv,res,f_out,rad_out,nc)

use lib_functions

implicit none

integer :: i
integer :: j
double precision :: x = 0.
double precision :: xx = 0.
double precision :: rconv = 0.
double precision :: sum = 0.
double precision :: fp = 0.
integer, intent(in) :: n_sort
double precision, intent(in) :: f1
double precision, intent(in) :: pas
double precision, intent(in) :: res
double precision, dimension(n_sort) :: rad
double precision, dimension(n_sort) :: rad_out
double precision, dimension(n_sort) :: f_out
double precision, dimension(-100:100) :: c

character(len=1), intent(in) :: iconv
integer, intent(out) :: nc
select case(iconv)
    case ("R")
        nc=0
        sum=1.d+00
        c(0)=1.d+00
        do i=1,n_sort
            rad_out(i)=rad(i)
            f_out(i)=f1+(dfloat(i)-0.5d+00)*pas
        end do
        return
    case ("H")
        fp=res/pas
        sum=0.
        nc=idnint(1.20d+00*fp)
        c(0)=2.d+00*0.54d+00*pas/res
        do i=1,nc
            x=pas*i
            xx=2.d+00*x/res
            if(dabs(xx-1.d+00).lt.1d-15) then
                c(i)=c(0)*0.4259259259d+00
            else
                c(i)=c(0)*(dsin(pi*xx)/(pi*xx)+(0.46d+00*xx*dsin(pi*xx)/(1-xx*xx))/(0.54d+00*pi))
            endif
            c(-i)=c(i)
        end do
    case ("S")
        nc=idnint(5.40d+00*fp)
        c(0)=pas/res
        do i=1,nc
            x=pas*i
            c(i)=pas*sin(pi*x/res)/(pi*x)
            c(-i)=c(i)
        end do
    case ("T")
        nc=idint(fp)
        c(0)=pas/res
        do i=1,nc
            x=pas*i
            c(i)=c(0)*(1.d+00-x/res)
            c(-i)=c(i)
        end do
    case ("G")
        nc=idnint(1.20d+00*fp)
        c(0)=2.d+00*pas*sqrt(log(2.d+00)/pi)/res
        do i=1,nc
            x=pas*i
            c(i)=c(0)*2.d+00**(-((2.d+00*x/res)**2))
            c(-i)=c(i)
        end do
    case ("L")
        nc=idnint(5.d+00*fp)
        c(0)=2.d+00*pas/(pi*res)
        do i=1,nc
            x=pas*i
            c(i)=c(0)/(1.d+00+(2.d+00*x/res)**2)
            c(-i)=c(i)
        end do
    case ("C")
        nc=idnint(0.5d+00*(fp-1.d+00))
        c(0)=1.d+00/(2.d+00*dfloat(nc)+1.d+00)
        do i=-nc,nc,1
            c(i)=c(0)
        end do
end select

if(nc.gt.100) then
    print 201, 2*nc+1
    201  format(' WARNING: Number of points in the convolution function larger than 100:',i5)
end if

do i=-nc,nc,1
    sum=sum+c(i)
end do
!~                 Spectral convolution
do i=nc+1,n_sort-nc
    rconv=0.
    do j=-nc,nc,1
        rconv=rconv+rad(i+j)*c(j)
    end do
    rad_out(i)=rconv/sum
    f_out(i)=f1+(dfloat(i)-0.5d+00)*pas
end do

return
end subroutine conv_k
