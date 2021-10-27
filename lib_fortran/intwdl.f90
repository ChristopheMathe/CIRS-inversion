subroutine intwdl(a,b,nstep,nstep0,step,step0)
implicit none
integer :: i
integer :: i1
double precision :: rap
double precision :: f
double precision :: d
integer, intent(in) :: nstep
integer, intent(in) :: nstep0
double precision, intent(in) :: step
double precision, intent(in) :: step0
double precision, dimension(nstep), intent(in)  :: a
double precision, dimension(nstep0), intent(out) :: b

if(nstep==nstep0) THEN
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
    end do
end if

return
end
