subroutine planck(f1,pas_sort,n_sort,nlay,tl,t_ground,pl,pl_ground)

use lib_functions

implicit none

integer :: i
integer :: j
integer, intent(in) ::  n_sort
integer, intent(in) ::  nlay
double precision, intent(in) :: f1
double precision, intent(in) :: pas_sort
double precision, intent(in) :: t_ground

double precision, dimension(nlay), intent(in) :: tl
double precision, dimension(n_sort), intent(out) :: pl_ground
double precision, dimension(nlay,n_sort), intent(out) :: pl

do i = 1, n_sort
    pl_ground(i) = pla(f1+pas_sort*(dfloat(i)-0.5), t_ground)
end do

do j = 1, nlay
    do i = 1, n_sort
        pl(j,i) = pla(f1+pas_sort*(dfloat(i)-0.5), tl(j))
    end do
end do

return
end
