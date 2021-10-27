subroutine dtau_h(dtau, etau0h, detau0, nstep, nstep0, step, step0, et, j, is, nlev, nsec)

implicit none
integer :: i
integer :: i1
integer, intent(in) :: j
integer, intent(in) :: is
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
integer, intent(in) :: nlev
integer, intent(in) :: nsec
double precision, intent(in) :: step
double precision, intent(in) :: step0
double precision :: d
double precision :: dx
double precision, dimension(0:500) :: et
double precision, dimension(nstep), intent(out) :: dtau
double precision, dimension(nstep0), intent(out) :: detau0
double precision, dimension(nlev,nsec,nstep0),intent(inout) :: etau0h

!~  DETAU0(*)   is Transmittance of the layer (j->j+1)
!~  ETAU0H(j,is,*) is Transmittance at level j at airmass # is
do i = 1, nstep
    if(dtau(i)>10.) then
        dtau(i) = 0.
    else
        d       = 50.*dtau(i)
        i1      = idint(d)
        dx      = 0.02*(d-dfloat(i1))
        dtau(i) = et(i1)*(1.-dx*(1.-0.5*dx))
    endif
end do

call intwdl(dtau,detau0,nstep,nstep0,step,step0)
etau0h(j,is,:)=etau0h(j+1,is,:)*detau0(:)
return
end
