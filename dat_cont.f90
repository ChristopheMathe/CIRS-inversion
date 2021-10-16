subroutine dat_cont(n2n2_g,n2ch4_g,ch4ch4_g,n2h2_g,tmod,n1,n2,tl,nlay,n2n2,n2ch4,ch4ch4,n2h2,knu,nfcont,lt)

use lib_functions
implicit none
integer :: i = 0
integer :: j = 0
integer :: l = 0
integer, intent(in) :: n1
integer, intent(in) :: n2
integer, intent(in) :: nlay
integer, intent(in) :: knu
integer, intent(in) :: nfcont
integer, intent(in) :: lt

double precision, dimension(nlay),intent(in) :: tl
real, dimension(lt),intent(in) :: tmod
real, dimension(lt,knu),intent(in) :: n2n2_g
real, dimension(lt,knu),intent(in) :: n2ch4_g
real, dimension(lt,knu),intent(in) :: ch4ch4_g
real, dimension(lt,knu),intent(in) :: n2h2_g

double precision, dimension(nlay,nfcont),intent(out) :: n2n2
double precision, dimension(nlay,nfcont),intent(out) :: n2ch4
double precision, dimension(nlay,nfcont),intent(out) :: ch4ch4
double precision, dimension(nlay,nfcont),intent(out) :: n2h2

double precision, dimension(lt) :: td
double precision, dimension(lt) :: a1
double precision, dimension(lt) :: b1
double precision, dimension(lt) :: c1
double precision, dimension(lt) :: d1


do 1 j=n1,n2
    do 2 i=1,lt
        td(i)=tmod(i)
        a1(i)=n2n2_g(i,j)
        b1(i)=n2ch4_g(i,j)
        c1(i)=ch4ch4_g(i,j)
        d1(i)=n2h2_g(i,j)
    2 END DO
    do 3 l=1,nlay
        n2n2(l,j-n1+1)=tina(tl(l),td,a1,lt)
        n2ch4(l,j-n1+1)=tina(tl(l),td,b1,lt)
        ch4ch4(l,j-n1+1)=tina(tl(l),td,c1,lt)
        n2h2(l,j-n1+1)=tina(tl(l),td,d1,lt)
    3 END DO
1 END DO
return
end subroutine dat_cont
