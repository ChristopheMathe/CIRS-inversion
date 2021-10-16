subroutine read_line2(freq1,freq2,alor,smin,ttest,erot,nlor,g2lor,w_in,s_in,g_in,e_in,w_out,s_out,&
e_out,n_out,n_nl,nlines,nbl_sp_max,nbl_sp,ncorps)

implicit none
integer :: i
integer :: k = 0
integer :: nbl_sp_max
integer :: ncorps
integer :: nbl

double precision :: f1
double precision :: f2
integer,dimension(ncorps) :: nlines
integer,dimension(ncorps), intent(in) :: nlor
integer,dimension(ncorps), intent(in) :: nbl_sp
double precision :: sm
double precision :: st
double precision :: a
double precision, intent(in) :: freq1
double precision, intent(in) :: freq2

double precision ,dimension(ncorps), intent(in) :: smin
double precision ,dimension(ncorps), intent(in) :: ttest
double precision ,dimension(ncorps), intent(in) :: erot
double precision ,dimension(ncorps), intent(in) :: alor

double precision, dimension(ncorps,5) :: g2lor
double precision, dimension(ncorps,nbl_sp_max), intent(in) :: w_in
double precision, dimension(ncorps,nbl_sp_max), intent(in) :: s_in
double precision, dimension(ncorps,nbl_sp_max), intent(in) :: g_in
double precision, dimension(ncorps,nbl_sp_max), intent(in) :: e_in
double precision, dimension(ncorps,nbl_sp_max), intent(out) :: w_out
double precision, dimension(ncorps,nbl_sp_max), intent(out) :: s_out
double precision, dimension(ncorps,nbl_sp_max), intent(out) :: e_out
integer, dimension(ncorps,nbl_sp_max), intent(out) :: n_out
integer, dimension(ncorps,5),intent(out) :: n_nl

integer ::read_begin = 0
integer ::read_end = 0
integer :: cpt = 0

do k = 1, ncorps
    cpt = 0
    w_out(k,:) = 0.
    s_out(k,:) = 0.
    e_out(k,:) = 0.
    n_out(k,:) = 0
    nlines(k) = 0
    n_nl(k,:) = 0
    read_begin = 1
    read_end = 0
    f1 = freq1 - alor(k)
    f2 = freq2 + alor(k)
    sm=smin(k)/(296./ttest(k))**erot(k)
    a=1.43878d+00*(1./ttest(k)-1./296.)

    if (w_in(k,1) > f2 .or. w_in(k,nbl_sp(k)-1) < f1) then
        read_begin = 0
        read_end = 0
    else
        do i =1, nbl_sp_max
            if(w_in(k,i) >= f1)then
                read_begin = i
                exit
            else
                continue
            end if
        end do
        do i = read_begin,nbl_sp_max
            if (w_in(k,i) > f2 .or. w_in(k,i) == 0) then
                read_end = i - 1
                exit
            else if (i == nbl_sp_max) then
                read_end = nbl_sp_max
            else
                continue
            end if
        end do
    end if

    if(read_begin < read_end) then
        do 666 nbl = read_begin,read_end
            st = s_in(k,nbl)*dexp(-a*e_in(k,nbl))
            if(st >= sm) then
                cpt = cpt +1
                nlines(k) = nlines(k)+1
                w_out(k,cpt) = w_in(k,nbl)
                s_out(k,cpt) = s_in(k,nbl)
                e_out(k,cpt) = e_in(k,nbl)
                if (nlor(k) == 1) then
                    n_out(k,cpt) = 1
                    n_nl(k,1) = n_nl(k,1) + 1
                    continue
                else
                    do i=1,nlor(k)-1
                        if(g_in(k,nbl) <= g2lor(k,i)) then
                            n_out(k,cpt) = i
                            n_nl(k,i) = n_nl(k,i) + 1
                            go to 666
                        endif
                    end do
                    n_out(k,cpt) = nlor(k)
                    n_nl(k,nlor(k)) = cpt - sum(n_nl(k,1:4))
                end if
            end if
        666 end do
    endif
end do
return
end subroutine read_line2
