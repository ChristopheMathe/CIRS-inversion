!~ SUBROUTINE AVG 

subroutine avg(etau0,etau,n_sort,nsec,nlev,nstep0,navg0,j,is)

implicit none

integer :: i
integer :: ik
integer :: i1

integer, intent(in) :: n_sort
integer, intent(in) :: nsec
integer, intent(in) :: nlev
integer, intent(in) :: nstep0
integer, intent(in) :: navg0
integer, intent(in) :: j
integer, intent(in) :: is

double precision, dimension(nstep0)          , intent(in)  :: etau0
double precision, dimension(nlev,nsec,n_sort), intent(out) :: etau

do ik=1,n_sort
    i1=(ik-1)*navg0+1
    etau(j,is,ik)=0.5*(etau0(i1)+etau0(i1+navg0))
    if(navg0 > 1) then
        do i=1,navg0-1
            etau(j,is,ik)=etau(j,is,ik)+etau0(i1+i)
        end do
        etau(j,is,ik)=etau(j,is,ik)/dfloat(navg0)
    endif
end do
return
end
