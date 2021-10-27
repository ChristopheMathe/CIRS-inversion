!~ SUBROUTINE CONV
!~ ---------------------------------------------------------------------
!~ R: Raw spectrum (no convolution)
!~ H: Hamming function (Res is first zero from center)
!~ S: Sinx/x (Res is first zero from center)
!~ T: Triangular function
!~ G: Gaussian function (Res is FWHM)
!~ L: Lorentz function (Res is FWHM)
!~ C: Square function (Res is the full width)

subroutine conv(f1,pas,n_sort,rad,is,iconv,res,f_out,rad_out,tb_out,nc,nsec)
use lib_functions

implicit none
integer :: j = 0
integer :: i = 0
integer :: ic = 0
double precision :: x = 0.
double precision :: xx = 0.
double precision :: sum = 0.
double precision :: rconv = 0.
double precision :: fp = 0.
integer, intent(in) :: n_sort
integer, intent(in) :: nsec
integer, intent(in) :: is
integer, intent(out) :: nc
double precision, intent(in) :: f1
double precision, intent(in) :: pas
double precision, intent(in) :: res

character(len=1), intent(in) :: iconv

double precision, dimension(nsec,n_sort) :: rad
double precision, dimension(n_sort), intent(out) :: rad_out
double precision, dimension(n_sort), intent(out) :: tb_out
double precision, dimension(n_sort), intent(out) :: f_out
double precision, dimension(-100:100) :: c


select case (iconv)
    case ("R")
        nc=0
        sum=1.
        c(0)=1.
        do i=1,n_sort
            rad_out(i)=rad(is,i)
            f_out(i)=f1+(dfloat(i)-0.5)*pas
            if(rad_out(i).gt.1.d-15) then
            tb_out(i)=bright(f_out(i),rad_out(i))
            else
            tb_out(i)=0.
            endif
        end do
        go to 14
    case ("H")
        fp=res/pas
        sum=0.
        nc=idnint(1.20*fp)
        c(0)=2.*0.54*pas/res
        do i=1,nc
            x=pas*i
            xx=2.*x/res
            if(abs(xx-1.d0).lt.1d-15) then
                c(i)=c(0)*0.4259259259d+00
            else
                c(i)=c(0)*(sin(pi*xx)/(pi*xx)+(0.46d0*xx*sin(pi*xx)/(1d0-xx*xx))/(0.54d0*pi))
            endif
            c(-i)=c(i)
        end do
    case ("S")
        fp=res/pas
        sum=0.
        nc=idnint(5.40*fp)
        c(0)=pas/res
        do i=1,nc
            x=pas*i
            c(i)=pas*sin(pi*x/res)/(pi*x)
            c(-i)=c(i)
        end do
        do i=-100, -nc
            c(i) = 0.
            c(-i) = 0.
        end do
        do i= nc, 100
            c(i) = 0.
            c(-i) = 0.
        end do
        
    case ("T")
        fp=res/pas
        sum=0.
        nc=idint(fp)
        c(0)=pas/res
        do i=1,nc
            x=pas*i
            c(i)=c(0)*(1.-x/res)
            c(-i)=c(i)
        end do
    case ("G")
        fp=res/pas
        sum=0.
        nc=idnint(1.20*fp)
        c(0)=2.*pas*sqrt(alog(2.)/pi)/res
        do i=1,nc
            x=pas*i
            c(i)=c(0)*2.**(-((2.*x/res)**2))
            c(-i)=c(i)
        end do
    case ("L")
        fp=res/pas
        sum=0.
        nc=idnint(5.00*fp)
        c(0)=2.*pas/(pi*res)
        do i=1,nc
            x=pas*i
            c(i)=c(0)/(1.+(2.*x/res)**2)
            c(-i)=c(i)
        end do
    case ("C")
        nc=idnint(0.5*(fp-1.))
        c(0)=1./(2.*dfloat(nc)+1.)
        do i=-nc,nc,1
            c(i)=c(0)
        end do
end select

do i=-nc,nc,1
   sum=sum+c(i)
end do
!~                  Spectral convolution
do i=nc+1,n_sort-nc
    rconv=0.
    do j=-nc,nc,1
        rconv=rconv+rad(is,i+j)*c(j)
    end do
    rad_out(i)=rconv/sum
    f_out(i)=f1+(dfloat(i)-0.5)*pas
    if(rad_out(i) > 1.d-15) then
        tb_out(i)=bright(f_out(i),rad_out(i))
    else
        tb_out(i)=0.
    endif
end do
14   continue
if (ic.le.0)then
    write(6,200)iconv,2*nc+1,sum,(c(i),i=-nc,nc)
    200 format(/,'Convolution ',a1,':',i4,' points, Sum =',f7.4,'Function: ',/,(2x,7e11.4))
    ic=1
end if
write(6,202)is,f_out(nc+1),rad_out(nc+1),tb_out(nc+1),f_out(n_sort-nc),rad_out(n_sort-nc),tb_out(n_sort-nc)
202 format(3x,'Calculation #',I2,3X,F10.4,' cm-1: rad =',E12.4,' erg, Tb =',F7.2,' K',/,21X,F10.4,' cm-1: rad =',&
E12.4,' erg, Tb =',F7.2,' K')
return
end
