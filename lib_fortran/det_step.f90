subroutine det_step(hw0,f,t,p,mass,ncorps,glor,elor,step,navg,pas_mult,pas_sort)

implicit none

integer, intent(out) :: navg
integer, intent(in) :: ncorps

double precision, parameter :: r=8.31441d+00
double precision, parameter :: c=2.997925d+08
double precision, parameter :: atm=1013.25
double precision, intent(in) :: pas_sort
double precision, intent(in) :: pas_mult
double precision, intent(in) :: t
double precision, intent(in) :: p
double precision, intent(inout) :: hw0
double precision, intent(in) :: f
double precision, intent(out) :: step
double precision, dimension(ncorps)  , intent(in) :: mass
double precision, dimension(ncorps,5), intent(in) :: glor
double precision, dimension(ncorps,5), intent(in) :: elor

integer :: k
double precision :: dop
double precision :: rlor
double precision :: hw

double precision :: al2 = 0.8325546112d0

do k = 1, ncorps
!~   --- Doppler halfwidths (cm-1) in layer J at frequency F1
    dop = f * al2 * dsqrt( 2.d+03*r*t / mass(k) ) / c
!~   --- Lorentz halfwidths (cm-1) in the layer J
    rlor = glor(k,1) * (p/atm) * (296./t)**elor(k,1)
    hw = dsqrt(dop**2 + rlor**2)
    if (hw < hw0) then
        hw0 = hw
    end if
end do

!~   --- HW0 is the lowest halfwidth from all absorbers
step = pas_mult * hw0
navg = 1 + int(pas_sort/step)
step = pas_sort/dfloat(navg)
return
end
