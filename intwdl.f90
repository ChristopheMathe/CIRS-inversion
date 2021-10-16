subroutine intwdl(a,b,nstep,nstep0,step,step0)
implicit none
integer :: i = 0
integer :: i1 = 0
double precision :: rap = 0d0
double precision :: f = 0d0
double precision :: d = 0d0
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
double precision, intent(in) :: step
double precision, intent(in) :: step0
double precision, dimension(nstep)  :: a
double precision, dimension(nstep0) :: b

if(nstep.eq.nstep0) THEN
    do i=1,nstep
        b(i)=a(i)
    end do
else
    rap=step0/step
    b(1)=a(1)
    b(nstep0)=a(nstep)
    do i=2,nstep0-1
        d=rap*(i-1)
        i1=1+idint(d)
        f=d+1.-i1
        b(i)=(1.-f)*a(i1)+f*a(i1+1)
!~         if(i==5248)then
!~             print*,'alba',i1,a(i1),a(i1+1),b(i)
!~         end if
    end do
end if

return
end
