subroutine dtau_e(dtau,etau0,detau0,nstep,nstep0,step,step0,et,nsec)
implicit none

integer :: i
integer :: i1
double precision :: d
double precision :: dx
integer, intent(in) :: nsec
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
double precision, intent(in) :: step
double precision, intent(in) :: step0
double precision, dimension(0:500), intent(in) :: et
double precision, dimension(nstep0), intent(in) :: dtau
double precision, dimension(nsec,nstep0), intent(inout) ::etau0
double precision, dimension(nstep0), intent(out) :: detau0

!~ ETAU0(1,*) is Optical depth at nadir
!~ DETAU0(*) is 2E3(ETAU0(1,*))

call intwdl(dtau,detau0,nstep,nstep0,step,step0)

do i=1,nstep0
    etau0(1,i)=etau0(1,i)+detau0(i)
    if(etau0(1,i).ge.10.) then
        detau0(i)=0.
    else
        d=50.*etau0(1,i)
        i1=idint(d)
        dx=d-dfloat(i1)
        detau0(i)=et(i1)*(1.-dx)+et(i1+1)*dx
    endif
end do


return
end
