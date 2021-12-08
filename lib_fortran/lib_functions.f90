Module lib_functions
implicit none
double precision, parameter :: hck = 1.43868d+00
double precision, parameter :: r = 8.31441d+00
double precision, parameter :: coef = 1.190605d-05
!double precision, parameter :: spi=1d0/(dsqrt(4d0*datan(1d0)))
!double precision, parameter :: pi= 4d0*datan(1d0)
double precision, parameter :: spi = 1./dsqrt(3.1415926536d+00)
double precision, parameter :: pi= 3.1415926536
double precision, parameter :: dp=0.0025d+00
double precision, parameter :: ds=0.01d+00
double precision, parameter :: yl=1.d+00/99.d+00
double precision, parameter :: xl=399.d+00

contains

!~ =====================================================================
double precision function TINA(A,X,Y,N)

integer,intent(in)             :: N
integer                        :: J
integer                        :: JL
integer                        :: JU
integer                        :: JM
double precision               :: frac
double precision,intent(in)    :: A
double precision, dimension(N) :: X
double precision, dimension(N) :: Y

JL = 1
JU = N
do while(JU-JL > 1)
        JM = (JU+JL) / 2
        if( (X(N) > X(1)) .EQV. (A > X(JM)) ) then
            JL = JM
        else
            JU = JM
        end if
end do
J=JL
if (X(J) == X(J+1)) then
  TINA=Y(J)
else
  FRAC=(A-X(J))/(X(J+1)-X(J))
  TINA=(1.-FRAC)*Y(J)+FRAC*Y(J+1)
end if
end function TINA
!~ =====================================================================

!~ =====================================================================
double precision function bright(f,r)

implicit none

double precision, intent(in) :: f
double precision, intent(in) :: r

bright=hck*f/dlog(1.+coef*f*f*f/r)

end function bright
!~ =====================================================================

!~ =====================================================================
double precision function pla(f,t)

implicit none
double precision, intent(in) :: f
double precision, intent(in) :: t

pla=coef*f*f*f/(dexp(hck*f/t)-1.)

end function pla
!~ =====================================================================

!~ =====================================================================
double precision function dpla(f,t)

implicit none
double precision, intent(in) :: f
double precision, intent(in) :: t

dpla=coef*f*f*f*(hck*f/(t*t))*dexp(hck*f/t)/(dexp(hck*f/t)-1.)**2

end function dpla

!~ =====================================================================

!~ =====================================================================
integer function nb_line(filename) result(nbl)

integer :: ios1
integer :: ios2
integer :: id
character(*), intent(in) :: filename

ios1 = 0
open(unit=id, file=filename, status='old', action='read', iostat=ios1)
if  (ios1.ne.0) then
      print *, 'Probleme avec le fichier ', filename
      stop
end if
ios2 = 0
nbl = 0
do while(ios2.eq.0)
    read(id, *, iostat=ios2)
    nbl = nbl + 1
end do
nbl = nbl - 1
close(id)

end function nb_line
!~ =====================================================================

!~ =====================================================================
double precision function voigt(x, y, v)
  double precision, intent(in) :: x
  double precision, intent(in) :: y
  double precision, dimension(400,100), intent(in) ::  v
  integer :: ip
  integer :: is
  double precision :: x2
  double precision :: y2
  double precision :: x2y2
  double precision :: a
  double precision :: b
  double precision :: f
  double precision :: p
  double precision :: s

  if(y>5.4) then
    x2 = x*x
    y2 = y*y
    x2y2 = x2 + y2
    voigt= spi * y * (1.+(3.*x2-y2)/(2.*x2y2*x2y2)+0.75*(5.*x2*x2+y2*y2-10.*x2*y2)/(x2y2*x2y2*x2y2*x2y2))/x2y2
  else
    if(y<=yl) then
      if(abs(x) <= 16.d0) then
        b=dexp(-x*x)
      else
        b= 0.d+00
      end if
      p=1.d+00/(1.d+00+x)
      ip = int(p/dp)
      f=p/dp-ip
      if(ip>=1) then
        if(ip<400) then
          a=(1.d+00-f)*v(ip,1)+f*v(ip+1,1)
        else
          a=v(ip,1)
        endif
      else
        a=f*v(ip+1,1)
      endif
      a=a*a
      f=y/yl
      voigt=(1.d+00-f)*b+f*a
    else
      if(x>xl) then
        s=y/(1.d+00+y)
        is=idint(s/ds)
        b=v(1,is)*xl/x
        b=b*b
        a=v(1,is+1)*xl/x
        a=a*a
        f=s/ds-is
        voigt=(1.d+00-f)*b+f*a
      else
        s = y / (1d0+y)
        p = 1d0 / (1.d+00+x)
        is = int(s/ds)
        ip = int(p/dp)
        f=p/dp-ip
        if(ip<400) then
          b=(1.d+00-f)*v(ip,is)+f*v(ip+1,is)
          a=(1.d+00-f)*v(ip,is+1)+f*v(ip+1,is+1)
        else
          b=v(ip,is)
          a=v(ip,is+1)
        endif
          b=b*b
          a=a*a
          f=s/ds-is
          voigt=(1.d+00-f)*b+f*a
      endif
  endif
endif
return
end function voigt
!~ =====================================================================

!~ =====================================================================
double precision function EINT(TT)
implicit none
integer :: I
integer :: K
integer :: M
integer :: N
integer :: P
double precision :: XN
double precision :: XI
double precision :: T
double precision :: TT
double precision :: A
double precision :: E
double precision :: X
double precision :: XPR
double precision :: TM
double precision :: TEMP
M=2
T=-TT
N=M
P=N
X=T+1.0d-38
IF(X-1.0) 100,200,200
 100  A=0.0
      K=N-1
      IF(K) 5,5,1
 1    TEMP=K
      DO 2 I=1,K
      A=A+1.0/TEMP
 2    TEMP=TEMP-1.0
      TEMP=1.0
      XPR=1.0
      DO 4 I=1,K
      XPR=-X*XPR/TEMP
 4    TEMP=TEMP+1.0
      EINT=-XPR*(dLOG(X)-A+.577215664901533d+00)
      GO TO 6
 5    E=-(dLOG(X)+.577215664901533d+00)
 6    TEMP=K
      TM=0.0
      XPR=1
 7    IF(TEMP)8,10,8
 8    E=EINT+XPR/TEMP
      IF(ABS(E-EINT))9,9,10
 9    EINT=E
      RETURN
 10   TM=TM+1.0
      TEMP=TEMP-1.0
      XPR=-X*XPR/TM
      EINT=E
      GO TO 7
 200  IF(X-1.5)201,202,202
 201  K=23
      GO TO 300
 202  IF(X-2.0)203,204,204
 203  K=16
      GO TO 300
 204  IF(X-2.5)205,206,206
 205  K=12
      GO TO 300
 206  IF(X-3.0)207,208,208
 207  K=10
      GO TO 300
 208  IF(X-4.0)209,210,210
 209  K=9
      GO TO 300
 210  IF(X-5.0)211,212,212
 211  K=7
      GO TO 300
 212  IF(X-6.0)213,214,214
 213  K=6
      GO TO 300
 214  IF(X-7.0)215,216,216
 215  K=5
      GO TO 300
 216  IF(X-15.0)217,218,218
 217  K=4
      GO TO 300
 218  K=2
 300  XN=N+K
      XI=K
      EINT=XI/(X+XN)
 301  XN=XN-1.0
      XI=XI-1.0
      IF(XI)303,303,302
 302  EINT=XI/(X+XN/(EINT+1.0))
      GO TO 301
 303  EINT=dEXP(-X)/(X+XN/(1.0+EINT))
return
end function EINT
End module lib_functions
