!~ SUBROUTINE AVG_K
subroutine avg_k(etau0,etau,n_sort,nlay,nstep0,navg0,j1,j2)


implicit none

integer :: i
integer :: ik
integer :: i1

integer, intent(in) :: n_sort
integer, intent(in) :: nlay
integer, intent(in) :: nstep0
integer, intent(in) :: navg0
integer, intent(in) :: j1
integer, intent(in) :: j2

double precision, dimension(nstep0)          , intent(in)  :: etau0
double precision, dimension(nlay,nlay,n_sort), intent(out) :: etau

do ik=1,n_sort ! n_sort = 105
    i1=(ik-1)*navg0+1
    etau(j1,j2,ik)=0.5*(etau0(i1)+etau0(i1+navg0))
    if(navg0.gt.1) THEN
        do i=1,navg0-1 ! navg0 = 378
            etau(j1,j2,ik)=etau(j1,j2,ik)+etau0(i1+i)
        end do
        etau(j1,j2,ik)=etau(j1,j2,ik)/dfloat(navg0)
    end if

end do

return
end
