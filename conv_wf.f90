!~ SUBROUTINE CONV_WF
!~ =====================================================================
!~ 	R: Raw spectrum (no convolution)
!~ 	H: Hamming function (Res is first zero from center)
!~ 	S: Sinx/x (Res is first zero from center)
!~ 	T: Triangular function
!~ 	G: Gaussian function (Res is FWHM)
!~ 	C: Square function (Res is the full width)
!~ =====================================================================
subroutine conv_wf(f1,pas,n_sort,etau,is,iconv,res,ap,nlay,f_out,nsec,nlev)

use lib_functions

implicit none

integer :: i = 0
integer :: j = 0
integer :: k = 0
integer :: nc = 0
double precision :: x = 0.
double precision :: xx = 0.
double precision :: wf = 0.
double precision :: sum = 0.
double precision :: fp = 0.
integer, intent(in) :: nlay
integer, intent(in) :: nlev
integer, intent(in) :: n_sort
integer, intent(in) :: nsec
integer, intent(in) :: is

double precision, intent(in) :: f1
double precision, intent(in) :: pas
double precision, intent(in) :: res
character(len=1), intent(in) :: iconv

double precision, dimension(n_sort) ::  f_out
double precision, dimension(-100:100) :: c
double precision, dimension(nlev,nsec,n_sort) :: etau
double precision, dimension(nlay+1) :: ap
double precision, dimension(nlev) :: econv


select case (iconv)
    case("R")
        nc=0
        sum=1.
        c(0)=1.
        do i=1,n_sort
            f_out(i)=f1+(dfloat(i)-0.5)*pas
            write(23,204)f_out(i)
            do j=nlay+1,1,-1
                econv(j)=etau(j,is,i)
                if(j.le.nlay) then
                    wf=(econv(j+1)-econv(j))/(sum*(ap(j)-ap(j+1)))
                end if
            end do
        end do
        return
        
    case("H")
        fp=res/pas
        sum=0.
        nc=idnint(1.20*fp)
        c(0)=2.*0.54*pas/res
        do i=1,nc
            x=pas*i
            xx=2.*x/res
            c(i)=c(0)*(dsin(pi*xx)/(pi*xx)+(0.46*xx*dsin(pi*xx)/(1.-xx*xx))/(0.54*pi))
            c(-i)=c(i)
        end do
        
    case("S")
        fp=res/pas
        sum=0.
        nc=idnint(5.40*fp)
        c(0)=pas/res
        do i=1,nc
            x=pas*i
            c(i)=pas*sin(pi*x/res)/(pi*x)
            c(-i)=c(i)
        end do
        
    case("T")
        fp=res/pas
        sum=0.
        nc=idint(fp)
        c(0)=pas/res
        do i=1,nc
            x=pas*i
            c(i)=c(0)*(1.-x/res)
            c(-i)=c(i)
        end do
        
    case("G")
        fp=res/pas
        sum=0.
        nc=idnint(1.20*fp)
        c(0)=2.*pas*sqrt(alog(2.)/pi)/res
        do i=1,nc
            x=pas*i
            c(i)=c(0)*2.**(-((2.*x/res)**2))
            c(-i)=c(i)
        end do
    case("C")
        fp=res/pas
        sum=0.
        nc=idnint(0.5*fp)
        c(0)=1./(2.*dfloat(nc)+1.)
        do i=-nc,nc,1
            c(i)=c(0)
        end do
        
end select

if(nc.gt.100) then
    print 201, 2*nc+1
    201  format(' WARNING: Number of points in the convolution function larger than 100:',i5)
end if

do  i=-nc,nc,1
    sum=sum+c(i)
end do

!~                 Spectral convolution

do i=nc+1,n_sort-nc
    f_out(i)=f1+(dfloat(i)-0.5)*pas
    write(23,204)f_out(i)
    do j=nlay+1,1,-1
        econv(j)=0.
        do k=-nc,nc,1
            econv(j)=econv(j)+etau(j,is,i+k)*c(k)
        end do
        if(j.le.nlay) then
            wf=(econv(j+1)-econv(j))/(sum*(ap(j)-ap(j+1)))
        end if
    end do
end do
write(6,202)is
202 format(3X,'Calculation of weighting functions for airmass #',I2)
204 format(30X,F12.6,' cm-1')
return
end
