subroutine dtau_v(dtau,etau0,detau0,nstep,nstep0,step,step0,et,is,nsec)

implicit none
integer :: i = 0
integer :: i1 = 0
integer, intent(in) :: is
integer, intent(in) :: nsec
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
double precision, intent(in) :: step
double precision, intent(in) :: step0
double precision :: d = 0d0
double precision :: dx = 0d0
double precision, dimension(0:500), intent(in) :: et
double precision, dimension(nstep)             :: dtau
double precision, dimension(nsec,nstep0)       :: etau0
double precision, dimension(nstep0)            :: detau0

!~  DETAU0(*)   is Transmittance of the layer (j->j+1)
!~  ETAU0(is,*) is Transmittance at level (j+1) at airmass # is
do i=1,nstep
    if(dtau(i).gt.10.) then
        dtau(i)=0.
    else
        d=50.*dtau(i)
        i1=idint(d)
        dx=0.02*(d-dfloat(i1))
        dtau(i)=et(i1)*(1.-dx*(1.-0.5*dx))
    endif
end do
!~ print*,dtau(1)
call intwdl(dtau,detau0,nstep,nstep0,step,step0)      
!~ print*,etau0(is,1),detau0(1)
etau0(is,:)=etau0(is,:)*detau0(:)

return
end
