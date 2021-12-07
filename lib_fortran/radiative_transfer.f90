!***        Abundance inversion program      ***
!*      Radiative transfer code for Titan
!*      N2-N2, N2-CH4, N2-H2 and CH4-CH4 continuum (coefficients read from unit 12: 4 spectroscopic files)
!*      Multi cloud non-grey model (only absorption: w0=0)
!*      FREQ1: First frequency (cm-1)
!*      FREQ2: Last frequency  (cm-1)
!***                                         ***
subroutine transfer_abundance(nlev, nlay, qn2, qh2, qar, pl, tl, ml, gl, qch4, ql, p, ap, t, nview, sec, lay_min, &
                              nb_mol ,mass, erot, et, freq1, freq2, pas_mult, pas_sort, file_trans, file_out, iconv, &
                              res, itrans, w_in, s_in, g_in, e_in, icloud, taucl, smin, ttest, nlor, alor, glor, elor, &
                              g2lor, nvib, vib, ndeg, n_k, icorps_k, matrix_k, rad, ikh, nfcont, nond, nond_max, sig, &
                              qex, nbl_sp_max, nbl_sp, f1_cont, df_cont, iter, nfreq_inv, qref, n2n2, n2ch4, ch4ch4, &
                              n2h2, v, corps, path_input, limbe, mode_inversion, matrix_t)
!$ use OMP_LIB
use lib_functions, only: voigt, tina, spi, hck, r
implicit none

! PARAMETERS
integer         , parameter :: nvoigt_max = 10000
double precision, parameter :: c = 2.997925d+08
double precision, parameter :: atm = 1013.25d+00
double precision, parameter :: t0 = 273.15d+00

! INTENT(IN)

!----- INTEGERS
integer, intent(in) :: nlev
integer, intent(in) :: nlay
integer, intent(in) :: nview
integer, intent(in) :: nb_mol
integer, intent(in) :: n_k
integer, intent(in) :: iter
integer, intent(in) :: nond_max
integer, intent(in) :: nfcont
integer, intent(in) :: icloud
integer, intent(in) :: itrans
integer, intent(in) :: nbl_sp_max
integer, intent(in) :: nfreq_inv
integer, intent(in) :: ikh
integer, intent(in), dimension(icloud)            :: nond
integer, intent(in), dimension(nview)             :: lay_min
integer, intent(in), dimension(nb_mol)            :: nvib
integer, intent(in), dimension(nb_mol)            :: nlor
integer, intent(in), dimension(nb_mol,5)          :: ndeg
integer, intent(in), dimension(n_k)               :: icorps_k
integer, intent(in), dimension(nb_mol)            :: nbl_sp
!----- DOUBLE PRECISION
double precision, intent(in) :: freq1
double precision, intent(in) :: freq2
double precision, intent(in) :: qn2
double precision, intent(in) :: qh2
double precision, intent(in) :: qar
double precision, intent(in) :: res
double precision, intent(in) :: pas_mult
double precision, intent(in) :: pas_sort
double precision, intent(in) :: df_cont
double precision, intent(in) :: f1_cont
double precision, intent(in), dimension(nlay)              :: ml
double precision, intent(in), dimension(nb_mol)            :: mass
double precision, intent(in), dimension(nlay,nfcont)       :: n2n2
double precision, intent(in), dimension(nlay,nfcont)       :: n2ch4
double precision, intent(in), dimension(nlay,nfcont)       :: ch4ch4
double precision, intent(in), dimension(nlay,nfcont)       :: n2h2
double precision, intent(in), dimension(nb_mol)            :: erot
double precision, intent(in), dimension(nb_mol)            :: alor
double precision, intent(in), dimension(nb_mol,5)          :: glor
double precision, intent(in), dimension(nb_mol,5)          :: elor
double precision, intent(in), dimension(nb_mol,5)          :: g2lor
double precision, intent(in), dimension(nb_mol,5)          :: vib
double precision, intent(in), dimension(nb_mol)            :: smin
double precision, intent(in), dimension(nb_mol)            :: ttest
double precision, intent(in), dimension(400,100)           :: v
double precision, intent(in), dimension(nlev)              :: p
double precision, intent(in), dimension(nlev)              :: ap
double precision, intent(in), dimension(nlev)              :: t
double precision, intent(in), dimension(nlay)              :: pl
double precision, intent(in), dimension(nlay)              :: tl
double precision, intent(in), dimension(nlay)              :: gl
double precision, intent(in), dimension(nlay)              :: qch4
double precision, intent(in), dimension(icloud,nlay)       :: taucl
double precision, intent(in), dimension(nlay,0:nb_mol)     :: ql
double precision, intent(in), dimension(0:500)             :: et
double precision, intent(in), dimension(icloud,nond_max)   :: sig
double precision, intent(in), dimension(icloud,nond_max)   :: qex
double precision, intent(in), dimension(nb_mol,nbl_sp_max) :: w_in
double precision, intent(in), dimension(nb_mol,nbl_sp_max) :: s_in
double precision, intent(in), dimension(nb_mol,nbl_sp_max) :: g_in
double precision, intent(in), dimension(nb_mol,nbl_sp_max) :: e_in
double precision, intent(in), dimension(icloud)            :: qref
double precision, intent(inout), dimension(nlay,nview)     :: sec

!----- CHARACTERS
character(len=1)  , intent(in) :: iconv
character(len=256), intent(in) :: file_out
character(len=256), intent(in) :: file_trans
character(*)      , intent(in) :: path_input
character(len=256), intent(in), dimension(nb_mol) :: corps

!* INTENT (OUT)

!----- DOUBLE PRECISION
double precision, intent(out), dimension(nview,nfreq_inv)          :: rad
double precision, intent(out), dimension(n_k,nlay,nview,nfreq_inv) :: matrix_k
double precision, intent(out), dimension(nlay,nview,nfreq_inv)     :: matrix_t

!* LOCAL VARIABLES

!----- INTEGERS
integer :: i = 0
integer :: i1 = 0
integer :: ic = 0
integer :: icalc = 0
integer :: ik = 0
integer :: is = 0
integer :: isort = 0
integer :: j = 0
integer :: j1 = 0
integer :: j2 = 0
integer :: k = 0
integer :: l = 0
integer :: nstep = 0
integer :: nstep0 = 0
integer :: nvoigt = 0
integer :: nvoigt_i = 0
integer :: nc = 0
integer :: n_sort1 = 0
integer :: navg0 = 0
integer :: navgj = 0
integer :: n_sort
integer :: lay_min_is
integer, dimension(nlay)              :: nst
integer, dimension(:,:), allocatable  :: n_out
integer, dimension(nb_mol)            :: nlines
integer, dimension(nb_mol,5)          :: n_nl

!----- DOUBLE PRECISION
double precision :: x = 0d0
double precision :: y = 0d0
double precision :: step0 = 0d0
double precision :: stepj = 0d0
double precision :: qa = 0d0
double precision :: qn = 0d0
double precision :: qh = 0d0
double precision :: v0 = 0d0
double precision :: rlor = 0d0
double precision :: time_begin = 0d0
double precision :: time_end = 0d0
double precision :: time_elapsed = 0d0
double precision :: f_xi = 0d0
double precision :: f_mult = 0d0
double precision :: f_cent = 0d0
double precision :: f = 0d0
double precision :: f1 = 0d0
double precision :: f2 = 0d0
double precision :: fdop = 0d0
double precision :: fc = 0d0
double precision :: fac_cont = 0d0
double precision :: f_xivib = 0d0
double precision :: dt = 0d0
double precision :: cmam = 0d0
double precision :: hw0 = 0d0
double precision :: hckt = 0d0
double precision :: hck296 = 0d0
double precision :: h0 = 0d0
double precision :: pj = 0d0
double precision :: tj = 0d0
double precision :: anc = 0d0
double precision :: ann = 0d0
double precision :: anh = 0d0
double precision :: acc = 0d0
double precision :: cloudj = 0d0
double precision :: avge = 0d0
double precision, dimension(:,:)  , allocatable :: w_out
double precision, dimension(:,:)  , allocatable :: s_out
double precision, dimension(:,:)  , allocatable :: e_out
!double precision :: s_over
double precision, dimension(:), allocatable :: s_over
double precision, dimension(:)    , allocatable :: dtau                 !(nfreq)
double precision, dimension(:)    , allocatable :: etau00               !(nfreq)
double precision, dimension(:)    , allocatable :: dtauk                !(nfreq)
double precision, dimension(:)    , allocatable :: detau0               !(nfreq)
double precision, dimension(:)    , allocatable :: rad_out              !(n_sort)
double precision, dimension(:)    , allocatable :: rad_outb              !(n_sort)
double precision, dimension(:)    , allocatable :: tb_out               !(n_sort)
double precision, dimension(:)    , allocatable :: f_out                !(n_sort)
double precision, dimension(:)    , allocatable :: rad_k                !(n_sort)
double precision, dimension(:)    , allocatable :: rad_kb               !(n_sort)
double precision, dimension(:)    , allocatable :: pl_ground            !(n_sort)
double precision, dimension(:)    , allocatable :: dtauk0               !(nfreq)
double precision, dimension(:,:)  , allocatable :: vgt
double precision, dimension(:,:)  , allocatable :: etau0                !(nview,nfreq)
double precision, dimension(:,:)  , allocatable :: tb                   !(n_sort)
double precision, dimension(:,:)  , allocatable :: pl_lay               !(nlay,n_sort)
double precision, dimension(:,:,:), allocatable :: etau                 !(nlev,nview,n_sort)
double precision, dimension(:,:,:), allocatable :: dtau_k               !(3,nlay,nfreq)
double precision, dimension(:,:,:), allocatable :: etau0_k              !(nlev,nview,nfreq)
double precision, dimension(:,:,:,:,:), allocatable :: etau_k               !(nlay,nlay,n_sort)
double precision, dimension(:,:,:,:,:), allocatable :: etau_k2              !(nlay,nlay,n_sort)
double precision, dimension(:,:,:), allocatable :: dtau_k0              !(nlay,nlay,n_sort)
double precision, dimension(:,:,:), allocatable :: etau_total !(nlay, nview, nfreq)
double precision, dimension(:,:)  , allocatable :: dpl_lay              !(nlay,n_sort)
double precision, dimension(:)    , allocatable :: dpl_ground           !(n_sort)
integer, intent(in) :: mode_inversion
logical, intent(in) :: limbe
integer :: nthread, ithread, nfreq
double precision :: z_debut, z_fin, z_tranche

integer :: nfreq_cut_off, i_min, i_mid, i_max, i_line_min, i_line_max, i_line_center, nr
double precision :: cut_off, wr, sr, er
!~ =====================================================================
!~ *** BEGIN TRANSFER ***
!~ =====================================================================

call cpu_time(time_begin)
write(*,*)"Begin transfert"

!$OMP PARALLEL DEFAULT(NONE) &

!$OMP SHARED(sec, pas_mult, pas_sort, erot, w_in, s_in, g_in, e_in, nlor, nvib, glor, elor, ndeg, vib, &
!$OMP        freq1, freq2, nthread, nbl_sp_max, nbl_sp, nb_mol, tl, pl, mass, limbe, iter, smin, ttest, nlay, nview, &
!$OMP        matrix_k, rad, et, icorps_k,  n2h2, ch4ch4, p, t, gl, ml, qref, lay_min, ql, ikh,&
!$OMP        n2ch4, n2n2, df_cont, f1_cont, nlev, corps, hck296, nond, taucl, qex, sig, qn2, qar, qh2, &
!$OMP        qch4, n_k, v, icloud, alor, g2lor, nfreq_inv, nfreq, etau_total, mode_inversion, matrix_t) &

!$OMP PRIVATE(etau_k, etau_k2, pl_lay, pl_ground, f1, f2, icalc, ithread, avge, s_over, x, v0, y,&
!$OMP         rlor, fdop, f_xivib, f_xi, cloudj, f_cent, anh, acc, anc, ann, i1, fc, dtau, qa, vgt, &
!$OMP         qh, qn, nst, nstep, nstep0, step0, navg0, hw0, stepj, navgj, cmam, h0, pj, dt, hckt, tj, detau0, &
!$OMP         dtauk, dtauk0, n_sort1, lay_min_is, etau, etau00, dpl_lay, dpl_ground, &
!$OMP         dtau_k, etau0_k, etau0, f, n_sort, z_tranche, z_fin, z_debut, f_mult,&
!$OMP         isort, n_nl, nlines, nvoigt, nvoigt_i, w_out, s_out, e_out, n_out, fac_cont, dtau_k0, &
!$OMP         nfreq_cut_off, i_min, i_mid, i_max, i_line_min, i_line_max, i_line_center, cut_off, wr, sr, er, nr)

!$ nthread = OMP_GET_NUM_THREADS()
!$ ithread = OMP_GET_THREAD_NUM()
n_nl(:,:) = 0
nlines(:) = 0
f2 = 0d0
hck296 = hck/296d0

!$OMP SINGLE
f1 = freq1
if(limbe .and. iter == 0) then
    do j = 1, nlay
        do is = 2, nview
            sec(j,is) = sec(j,is) / sec(j,1)
        end do
    end do
end if

nfreq = int((freq2 - freq1)/pas_sort) + 1
write(*,fmt='(a,i5,i5,a)')'There are ',nfreq, nfreq_inv,' points'
allocate(etau_total(nlay, nview, nfreq))
!$OMP END SINGLE

!	Determination of the calculation step (cm-1)
n_sort1 = (ithread*nfreq/nthread)
z_fin  = ((ithread+1)*nfreq/nthread) - 1
n_sort = (int(z_fin) - n_sort1) + 1
write(*,fmt='(i2,a,i4,a,i4,a,i5)')ithread,', start point=',n_sort1, 'end point', int(z_fin),', number point=', n_sort

z_debut   = freq1 + (n_sort1)*pas_sort
z_fin     = freq1 + (z_fin)*pas_sort
z_tranche = (z_fin - z_debut)
write(*,*)ithread, ' Calculation from', z_debut, ' to', z_fin, ' with ', z_tranche, ' cm-1'

allocate(pl_lay(nlay,n_sort), pl_ground(n_sort), dpl_lay(nlay,n_sort), dpl_ground(n_sort))
f1 = z_debut
f2 = z_fin
f = 0.5 * (f1 + f2)
hw0 = 1.d+03
do j = nlay, 1, -1
    call det_step(hw0, f1, tl(j), pl(j), mass, nb_mol, glor, elor, step0, navg0, pas_mult, pas_sort)
end do
nstep0 = navg0*(idint((f2-f1)/pas_sort)+1) + 1
write(*,fmt='(10x,16(1H*),a,i3,x,16(1H*))')' Spectral Interval #', ithread+1
write(*,fmt='(a,f10.4,a,f10.4,a,e11.4,a)')' Calculation from', f1, ' to', f2, ' with a step of', step0, ' cm-1'

!** Input: Molecular line parameters
write(*,*)'INPUT: MOLECULAR LINE PARAMETERS'
allocate(w_out(nb_mol,nbl_sp_max), s_out(nb_mol,nbl_sp_max), e_out(nb_mol,nbl_sp_max), n_out(nb_mol,nbl_sp_max))
w_out(:,:) = 0d0
s_out(:,:) = 0d0
e_out(:,:) = 0d0
n_out(:,:) = 0
call read_line2(f1, f2, alor, smin, ttest, erot, nlor, g2lor, w_in, s_in, g_in, e_in, w_out, s_out, e_out, n_out, &
                n_nl, nlines, nbl_sp_max, nbl_sp, nb_mol)
do k=1,nb_mol
    write(*,fmt='(a4,a,i5,a,f5.1,a,e10.3,a,f4.1)')corps(k), ': ', nlines(k), ' lines with S(', ttest(k), ' K) >', &
                                                  smin(k), ' n_rot =', erot(k)
    if (nlines(k) > 0)then
        do i=1,nlor(k)
            write(*,fmt='(4x,i7,a,f6.3,a,f5.3)')n_nl(k,i), ' lines with Lorentz halfwidth ⁼', glor(k,i), &
                                               ' cm-1, n = -', elor(k,i)
        end do
        if (nvib(k) > 0) then
            write(*,fmt='(9x,i2,a,6(f8.2,a,i1,a,1x))')nvib(k),' vibrational modes: ', &
                                                      (vib(k,i), ' (', ndeg(k,i), ')', i=1, nvib(k))
        end if
    end if
end do

!**		Calculation of the optical depths in the layers
allocate(etau0(nview,nstep0), etau0_k(nlev,nview,nstep0), dtau_k(n_k,nlay,nstep0), dtauk0(nstep0), etau00(nstep0), &
         dtauk(nstep0), detau0(nstep0))
etau0(:,:)     = 0d0
etau0_k(:,:,:) = 0d0
dtau_k(:,:,:)  = 0d0
dtauk0(:) = 0d0
etau00(:) = 0d0
dtauk(:)  = 0d0

allocate(etau(nlev,nview,n_sort))
do is = 1, nview
    if(nview /= 0) then
        etau0_k(nlev,is,:) = 1d0
        etau0(is,:)        = 1d0
        etau(nlev,is,:)    = 1d0
    end if
end do

do j = nlay, 1, -1
    detau0(1:nstep0) = 0d0
    tj = tl(j)
    hckt = hck/tj
    pj = pl(j)
    h0 = 1d+07 * r * t0 / (ml(j) * gl(j))
    cmam = 2.6867d+19 * spi * h0 * sec(j,1) * (p(j)-p(j+1))/atm
    hw0 = 1d+03
    call det_step(hw0, f, tj, pj, mass, nb_mol, glor, elor, stepj, navgj, pas_mult, pas_sort)
    nstep = navgj * n_sort + 1
    nst(j) = nstep
    if (nstep > nstep0) then
        nstep = nstep0
        nst(j)= nstep
        navgj = (nstep - 1)/n_sort
        stepj = pas_sort/dfloat(navgj)
    end if

    allocate(dtau(nstep))
    dtau(1:nstep) = 0d0

    ! Opacity due to N2-N2, CH4-CH4, N2-CH4, N2-H2
    fac_cont = (t0/tj) * h0 * pj * (p(j)-p(j+1)) / (atm*atm)
    qn = ( 1d0 - qch4(j) ) / ( 1d0 + (qh2+qar)/qn2 )
    qh = ( 1d0 - qch4(j) ) / ( 1d0 + (qn2+qar)/qh2 )
    qa = 0d0
    if (qar /= 0) then
        qa = ( 1d0 - qch4(j) ) / ( 1d0 + (qn2+qh2)/qar )
    end if

    ! Assume N2-N2 and N2-Ar opacities are the same
    qn = qn + qa
    do i = 1, nstep
        fc = (f1 + stepj*dfloat(i-1) - f1_cont) / df_cont
        i1 = idint(fc) + 1
        fc = fc - i1  + 1d0
        ann = (1d0-fc) * n2n2(j,i1)   + fc*n2n2(j,i1+1)
        anc = (1d0-fc) * n2ch4(j,i1)  + fc*n2ch4(j,i1+1)
        acc = (1d0-fc) * ch4ch4(j,i1) + fc*ch4ch4(j,i1+1)
        anh = (1d0-fc) * n2h2(j,i1)   + fc*n2h2(j,i1+1)
        dtau(i) = sec(j,1) * fac_cont * ( qn * (qn*ann + qch4(j)*anc + qh*anh) + qch4(j)*qch4(j)*acc )
    end do

    ! Cloud opacity
    if(icloud /= 0) then
        do ik = 1, n_sort
            f_cent = f1 + pas_sort * (dfloat(ik) - 0.5d0)
            cloudj = 0d0
            do ic = 1, icloud
                cloudj = cloudj + sec(j,1) * taucl(ic,j) * tina(f_cent, sig(ic,:), qex(ic,:), nond(ic)) / qref(ic)
            end do
            i1 = (ik-1) * navgj
            do i = i1+1, i1+navgj
                dtau(i) = dtau(i) + cloudj
            end do
            if(ikh > 0) then
                do i = i1+1, i1+navgj
                    dtau_k(ikh,j,i) = cloudj
                end do
            end if
        end do
        dtau(nstep) = dtau(nstep) + cloudj
        if(ikh > 0) then
            dtau_k(ikh,j,nstep) = cloudj
        end if
    end if

    !*	Iteration over the nb_mol absorbers
    do k = 1, nb_mol
        if(nlines(k) > 0) then
            f_xi = (296d0 / tj)**erot(k)
            f_xivib = 1d0
            if(nvib(k) > 0) then
                do i = 1, nvib(k)
                    f_xivib = f_xivib * ( (1.d0-dexp(-vib(k,i)*hckt)) / (1.d0-dexp(-vib(k,i)*hck296)) )**ndeg(k,i)
                end do
            end if
            nvoigt = min0(nvoigt_max, 1 + nint(alor(k)/stepj))
            allocate(vgt(5,nvoigt))
            vgt(1:5,1:nvoigt) = 0d0
            dtauk(1:nstep0) = 0d0

            !*  --- 1/Doppler parameter (cm-1)
            fdop = 1. / (f * dsqrt(2.d+03*r*tj/mass(k)) / c)

            !*  --- Lorentz halfwidth (cm-1)
            nvoigt_i = nvoigt
            do ik = nlor(k), 1, -1
                rlor = glor(k,ik) * (pj/atm) * (296d0/tj)**elor(k,ik)
                y = fdop * rlor
                v0 = voigt(0d0, y, v)
                do i = 1, nvoigt
                !do i = 1, nvoigt_i
                    x = fdop * stepj * (i-1)
                    vgt(ik,i) = voigt(x, y, v)
                    if(vgt(ik,i) < 0.5d-09*v0 .and. ik == nlor(k)) then
                        !nvoigt_i = i
                        nvoigt = i
                        exit
                    end if
                end do
            end do
            allocate(s_over(nlines(k)))
            s_over(:) = 0d0
            dt = hckt - hck296
            do l = 1, nlines(k)
                s_over(l) = s_out(k,l) * dexp(-e_out(k,l)*dt) * (1.d0-dexp(-w_out(k,l)*hckt)) / &
                            (1.d0 - dexp(-w_out(k,l)*hck296))
            end do

            if(nstep < nvoigt) then
                    !* (f2-f1) < cutoff
                call tau_lines1(w_out(k,:), s_over, n_out(k,:), nlines(k), dtauk, f1, stepj, nstep, nstep0, vgt, nvoigt)
            else
                if(nstep < 2*nvoigt) then
                    !*	cutoff <= (f2-f1) < 2*cutoff
                   call tau_lines2(w_out(k,:), s_over, n_out(k,:), nlines(k), dtauk, f1, stepj, nstep, nstep0, vgt, &
                                   nvoigt)
                else
                    !* (f2-f1) >= 2*cutoff
                    call tau_lines3(w_out(k,:), s_over, n_out(k,:), nlines(k), dtauk, f1, stepj, nstep, nstep0, vgt, &
                                    nvoigt)
                end if
            end if
           ! OTHER METHOD
            ! Iteration over the NLINES(K) lines of Absorber K
!              cut_off = fdop * stepj * (nvoigt_i-1)
!              dt = hckt - hck296
!              do l = 1, nlines(k)
!                wr = w_out(k,l)
!                nr = n_out(k,l)
!                er = e_out(k,l)
!                sr = s_out(k,l)
!
!                ! Exclude lines too far from the wavenumber interval
!                if(wr > f2 + cut_off) then
!                  cycle
!                else if(wr < f1 - cut_off) then
!                  cycle
!                end if
!
!                s_over = sr * dexp(-er*dt) * (1.d0-dexp(-wr*hckt)) / &
!                         (1.d0 - dexp(-wr*hck296))
!
!!                ! Find the index of the wavenumber interval where the line center is located
!                i_line_center = nint((wr - f1) / stepj) + 1
!
!                ! Find the index of the borders of the line profile
!                i_line_min = i_line_center - nvoigt_i + 1
!                i_line_max = i_line_center + nvoigt_i - 1
!
!                ! Find the absorption cross section calculation indices
!                ! 0 <= i_mid <= n + 1 (and not 1 <= i_mid <= n) because we need to calculate at i = 1 and i = n
!                i_min = max(i_line_min, 1)
!                i_mid = min(max(i_line_center, 0), nstep0 + 1)
!                i_max = min(i_line_max, nstep)
!
!                ! Left side of the line profile (center of the line not included)
!                if(i_line_center > 1) then  ! line_wavenumber > wavenumber_min
!                  do i = i_min, i_mid - 1
!                    dtauk(i) = dtauk(i) + s_over * vgt(nr, i_line_center - i + 1)
!                  end do
!                end if
!
!                ! Right side of the line profile (center of the line not included)
!                if(i_line_center < nstep) then  ! line_wavenumber < wavenumber_max
!                  do i = i_mid + 1, i_max
!                    dtauk(i) = dtauk(i) + s_over * vgt(nr, i - i_line_center + 1)
!                  end do
!                end if
!
!                ! Center of the line
!                if(i_line_center >= 1 .and. i_line_center <= nstep) then
!                  dtauk(i_mid) = dtauk(i_mid) + s_over * vgt(nr, 1)
!                end if
!            end do

            f_mult = fdop * f_xi * f_xivib * cmam
            f_mult = f_mult * ql(j,k)

            ! Opacity due to molecule k
            dtau = dtau + f_mult * dtauk(:nstep)
            do ik = 1, n_k
                if(k == icorps_k(ik)) then
                    dtau_k(ik,j,:) = f_mult * dtauk(:)
                end if
            end do
            deallocate(vgt, s_over)
        end if
    end do
    do ik = 1, n_k
        dtauk(:) = dtau_k(ik,j,:)
        call intwdl(dtauk,dtauk0,nstep,nstep0,stepj,step0)
        dtau_k(ik,j,:) = dtauk0(:)
    end do

    if(nview == 0) then
        call dtau_e(dtau,etau0,detau0,nstep,nstep0,stepj,step0,et,nview)
        call avg(detau0,etau,n_sort,nview,nlev,nstep0,navg0,j,1)
    else
        do is = 1, nview
            lay_min_is = lay_min(is)
            if(.not.limbe .or. (limbe.and.lay_min_is == 0)) then
                if(is == 1) then
                    dtauk(:nstep) = dtau
                else
                    dtauk(:nstep) = dtau *sec(j,is)
                endif
                call dtau_v(dtauk,etau0,detau0,nstep,nstep0,stepj,step0,et,is,nview)
                dtauk = etau0(is,:)
                etau0_k(j,is,:)=etau0(is,:)
                call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,j,is)
            else
                if(is == 1) then
                    dtauk(:nstep)=dtau(:)
                    call dtau_h(dtauk,etau0_k,detau0,nstep,nstep0,stepj,step0,et,j,is,nlev,nview)
                else
                    if(lay_min_is <= j) then
                        dtauk(:nstep) = dtau(:) * sec(j,is)
                        call dtau_h(dtauk,etau0_k,detau0,nstep,nstep0,stepj,step0,et,j,is,nlev,nview)
                    endif
                endif
            end if
        end do
    end if
    deallocate(dtau)
end do ! 18

!  For horizontal viewing, compute exp(-tau(j)) - exp(-2*tau(1)) / exp(-tau(j))
if (limbe) then
    do is = 1, nview
        lay_min_is = lay_min(is)
        if(lay_min_is == 0) then
            cycle
        end if
        etau00(1:nstep0) = etau0_k(lay_min_is,is,1:nstep0)**2
        do j = nlay, lay_min_is, -1
            do i = 1, nstep0
                if(etau00(i) > 0d0) then
                    dtauk(i) = etau0_k(j,is,i) - etau00(i)/etau0_k(j,is,i)
                else
                    dtauk(i) = etau0_k(j,is,i)
                endif
            end do
            call avg(dtauk, etau, n_sort, nview, nlev, nstep0, navg0, j, is)
        end do
        dtauk(1:nstep0) = 1. - etau00(1:nstep0)
        call avg(dtauk, etau, n_sort, nview, nlev, nstep0, navg0, nlev, is)
    end do
end if

! Calculation of the Planck function (erg s-1 cm-2 sr-1/cm-1)
write(*,*)'Calculation of the Planck function'
call planck(f1, pas_sort, n_sort, nlay, tl, t(1), pl_lay, pl_ground)

! Calculation of the outgoing radiances
write(*,*)'Calculation of the outgoing radiances'
do is = 1, nview
    lay_min_is = lay_min(is)
    do i = 1, n_sort
        isort = n_sort1 + i
        rad(is,isort) = 0d0
        if(limbe .AND. lay_min_is > 0) then
            do j = lay_min_is, nlay
                rad(is,isort) = rad(is,isort) + pl_lay(j,i)*(etau(j+1,is,i)-etau(j,is,i))
            end do
        else
            do j=1,nlay
                rad(is,isort) = rad(is,isort) + pl_lay(j,i)*(etau(j+1,is,i)-etau(j,is,i))
            end do
            rad(is,isort) = rad(is,isort) + pl_ground(i)*etau(1,is,i)
        end if
    end do
end do

! mode_inversion = 2: calculation of the matrix_t (for temperature)
if (mode_inversion == 2) then
    call dplanck(f1,pas_sort,n_sort,nlay,tl,t(1),dpl_lay,dpl_ground)
    write(*,*)'Calculation of the temperature matrix'
    do is = 1, nview
        lay_min_is = lay_min(is)
        do i = 1, n_sort
            isort = n_sort1+i
            if(limbe .and. lay_min_is > 0) then
                do j = lay_min_is, nlay
                    matrix_t(j,is,isort) = (etau(j+1,is,i)-etau(j,is,i)) * dpl_lay(j,i)
                end do
                if(lay_min(is) > 1) then
                    do j = 1, lay_min_is-1
                        matrix_t(j,is,isort) = 0d0
                    end do
                end if
            else
                do j=1,nlay
                    matrix_t(j,is,isort) = (etau(j+1,is,i) - etau(j,is,i)) * dpl_lay(j,i)
                end do
                matrix_t(1,is,isort) = matrix_t(1,is,isort) + etau(1,is,i) * dpl_ground(i)
            end if
        end do
    end do
end if

! mode_inversion = 1 and 2: calculation of the matrix_k (for molecules)
if (mode_inversion /= 0) then
    write(*,*)'Calculation of the K body matrice'

    ! calcul de etau_k
    write(*,*)'Calcul of etau_k'
    allocate(etau_k(n_k,nview,nlay,nlay,n_sort), etau_k2(n_k,nview,nlay,nlay,n_sort))
    !allocate(dtau_k0(n_k,nview,nlay,nlay,nstep0)) ! impossible de faire dtau_k0(n_k,nview,nlay,nlay,nstep0) car trop de
                                                   ! memoire demandee
    allocate(dtau_k0(nstep0,nlay,nlay))
    etau_k(:,:,:,:,:)  = 0d0
    etau_k2(:,:,:,:,:) = 0d0
    dtau_k0(:,:,:)     = 0d0

    ! is = 1
    lay_min_is = lay_min(1)
    do ik = 1, n_k
        do i = 1, nstep0
            do j1 = 1, nlay
                ! layer between 1 and j1
                do j2 = 1, j1
                    dtau_k0(i,j1,j2) = dtau_k(ik,j1,i) * etau0_k(j2,1,i)
                end do
                if((.not. limbe) .OR. (limbe .AND. lay_min_is == 0)) then
                    exit
                end if
                ! layer between j1+1 and nlay
                do j2 = j1+1, nlay
                    dtau_k0(i,j1,j2) = dtau_k(ik,j1,i)*etau0_k(j2,1,i)
                end do
            end do
        end do
        do isort = 1, n_sort -1
            i1 = (isort-1)*navg0 + 1
            do j1 = 1, nlay
                do j2 = 1, nlay
                   etau_k(ik,1,j1,j2,isort) = 0.5 * (dtau_k0(i1,j1,j2)+dtau_k0(i1+navg0,j1,j2))
                   if(navg0 > 1) then
                       do i = 1, navg0-1
                           etau_k(ik,1,j1,j2,isort) = etau_k(ik,1,j1,j2,isort) + dtau_k0(i1+i,j1,j2)
                       end do
                       etau_k(ik,1,j1,j2,isort) = etau_k(ik,1,j1,j2,isort)/dfloat(navg0)
                   end if
                end do
            end do
       end do
    end do
    ! is > 2
    do is = 2, nview
        lay_min_is = lay_min(is)
        do ik = 1, n_k
            do i = 1, nstep0
                do j1 = 1, nlay
                    ! layer between 1 and j1
                    do j2 = 1, j1
                        dtau_k0(i,j1,j2) = dtau_k(ik,j1,i)*sec(j1,is)*etau0_k(j2,is,i)
                    end do

                    if((.not. limbe) .OR. (limbe .AND. lay_min_is == 0)) then
                        exit
                    end if

                    ! layer between j1+1 and nlay
                    do j2 = j1+1, nlay
                        dtau_k0(i,j1,j2) = dtau_k(ik,j1,i)*sec(j1,is)*etau0_k(j2,is,i)
                    end do
                end do
            end do
            do isort = 1, n_sort -1
                i1 = (isort-1)*navg0 + 1
                do j1 = 1, nlay
                    do j2 = 1, nlay
                       etau_k(ik,is,j1,j2,isort) = 0.5 * (dtau_k0(i1,j1,j2)+dtau_k0(i1+navg0,j1,j2))
                       if(navg0 > 1) then
                           do i = 1, navg0-1
                               etau_k(ik,is,j1,j2,isort) = etau_k(ik,is,j1,j2,isort) + dtau_k0(i1+i,j1,j2)
                           end do
                           etau_k(ik,is,j1,j2,isort) = etau_k(ik,is,j1,j2,isort)/dfloat(navg0)
                       end if
                    end do
                end do
           end do
        end do
    end do

    ! calcul de etau_k2
    write(*,*)'Calcul of etau_k2'
    do ik = 1, n_k
        lay_min_is = lay_min(1)
        dtau_k0(:,:,:) = 0d0
        do i = 1, nstep0
            do j1 = 1, nlay
                do j2 = lay_min_is, nlay
                    if (etau0_k(lay_min_is,1,i) > 0d0) then
                        dtau_k0(i,j1,j2) = dtau_k(ik,j1,i) * etau0_k(lay_min_is,1,i)**2/etau0_k(j2,1,i)
                    end if
                end do
            end do
        end do
        do isort = 1, n_sort -1
            i1 = (isort-1)*navg0 + 1
            do j1 = 1, nlay
                do j2 = 1, nlay
                   etau_k2(ik,1,j1,j2,isort) = 0.5d0 * (dtau_k0(i1,j1,j2)+dtau_k0(i1+navg0,j1,j2))
                   if(navg0 > 1) then
                       do i = 1, navg0-1
                           etau_k2(ik,1,j1,j2,isort) = etau_k2(ik,1,j1,j2,isort) + dtau_k0(i1+i,j1,j2)
                       end do
                       etau_k2(ik,1,j1,j2,isort) = etau_k2(ik,1,j1,j2,isort)/dfloat(navg0)
                   end if
                end do
            end do
        end do

        do is = 2, nview
            lay_min_is = lay_min(is)
            dtau_k0(:,:,:) = 0d0
            do i = 1, nstep0
                do j1 = 1, nlay
                    do j2 = lay_min_is, nlay
                        if (etau0_k(lay_min_is,is,i) > 0d0) then
                            dtau_k0(i,j1,j2) = dtau_k(ik,j1,i) * sec(j1,is) * etau0_k(lay_min_is,is,i)**2/etau0_k(j2,is,i)
                        endif
                    end do
                end do
            end do
            do isort = 1, n_sort -1
                i1 = (isort-1)*navg0 + 1
                do j1 = 1, nlay
                    do j2 = 1, nlay
                       etau_k2(ik,is,j1,j2,isort) = 0.5d0 * (dtau_k0(i1,j1,j2)+dtau_k0(i1+navg0,j1,j2))
                       if(navg0 > 1) then
                           do i = 1, navg0-1
                               etau_k2(ik,is,j1,j2,isort) = etau_k2(ik,is,j1,j2,isort) + dtau_k0(i1+i,j1,j2)
                           end do
                           etau_k2(ik,is,j1,j2,isort) = etau_k2(ik,is,j1,j2,isort)/dfloat(navg0)
                       end if
                    end do
                end do
           end do
        end do
    end do

    ! Calcul of matrix K
    write(*,*)'Fill matrix K'
    do ik = 1, n_k
        do is = 1, nview
            lay_min_is = lay_min(is)
            do i = 1, n_sort
                isort = n_sort1 + i
                if(.not. limbe .OR. (limbe .AND. lay_min_is == 0)) then
                    do j1 = 1, nlay
                        matrix_k(ik,j1,is,isort) = etau_k(ik,is,j1,j1,i)*pl_lay(j1,i) - etau_k(ik,is,j1,1,i)*pl_ground(i)
                        if(j1 > 1) then
                            do j2 = 1, j1-1
                             matrix_k(ik,j1,is,isort) = matrix_k(ik,j1,is,isort) - &
                                                        (etau_k(ik,is,j1,j2+1,i)-etau_k(ik,is,j1,j2,i))*pl_lay(j2,i)
                            end do
                        end if
                    end do
                else
                    do j1 = lay_min_is, nlay-1 !si pas de -1 je dépasse dim2 de etau_k2
                      matrix_k(ik,j1,is,isort) = pl_lay(j1,i)*(etau_k(ik,is,j1,j1,i) - etau_k2(ik,is,j1,j1,i) + &
                                                 2d0*etau_k2(ik,is,j1,j1+1,i))
                        if(j1 > lay_min_is) then
                            do j2 = lay_min_is, j1-1
                                matrix_k(ik,j1,is,isort) = matrix_k(ik,j1,is,isort) - &
                                                           ( etau_k(ik,is,j1,j2+1,i) - etau_k(ik,is,j1,j2,i) + &
                                                             etau_k2(ik,is,j1,j2,i) - etau_k2(ik,is,j1,j2+1,i) ) &
                                                           * pl_lay(j2,i)
                            end do
                        end if
                        if(j1 < nlay) then
                            do j2 = j1+1, nlay-1
                                matrix_k(ik,j1,is,isort) = matrix_k(ik,j1,is,isort) - (etau_k2(ik,is,j1,j2,i)&
                                                           - etau_k2(ik,is,j1,j2+1,i)) * 2d0 * pl_lay(j2,i)
                            end do
                        end if
                    end do
                    if(lay_min(is) > 1) then
                        do j1 = 1, lay_min_is-1
                            matrix_k(ik,j1,is,isort) = 0d0
                        end do
                    end if
                end if
            end do
        end do
    end do
    deallocate(etau_k, etau_k2, dtau_k0)
end if

!	Writing the averaged exp(-tau) for Calculation # 1
write(*,fmt='(10x,a,3x,a,9x,a,4x,a,3x,a,2x,a)')'P (mb)','Transmittance','P (mb)','T (K)','Airmass','N_step'
if(.not. limbe) then
    write(*,fmt='(i4,e14.3,f12.7,5x,e14.3,f8.2,f8.4,i8)')nlev, p(nlev), 1d0
else
    avge =0d0
    do i = 1, n_sort
        avge = avge + etau(nlev,1,i)
    end do
    write(*,fmt='(i4,e14.3,f12.7,5x,e14.3,f8.2,f8.4,i8)')nlev, p(nlev), avge/(n_sort*1d0)
end if

do j = nlay, 1, -1
    avge = 0d0
    do i = 1, n_sort
        avge = avge + etau(j,1,i)
    end do
    avge = avge/(n_sort*1d0)
    write(*,fmt='(i4,e14.3,f12.7,5x,e14.3,f8.2,f8.4,i8)')j, p(j), avge, pl(j), tl(j), sec(j,1), nst(j)
end do

do j = 1, nlay
    do is = 1, nview
        do i = 1, n_sort
            isort = n_sort1 + i
            etau_total(j,is,isort) = etau(j,is,i)
        end do
    end do
end do

deallocate(etau00, etau, etau0, etau0_k, detau0, dtauk0, dtauk, dtau_k, pl_lay, pl_ground, &
           w_out, s_out, e_out, n_out)
!$OMP END PARALLEL

! Convolution of the spectrum by the apparatus function
write(*,*)'Convolution of the spectrum'
open(unit=22, file=file_out, status='unknown', form='formatted')
allocate(tb(nview,nfreq), rad_out(nfreq), tb_out(nfreq), f_out(nfreq))
tb(:,:)    = 0d0
rad_out(:) = 0d0
tb_out(:)  = 0d0
f_out(:)   = 0d0
do is = 1, nview
    call conv(freq1, pas_sort, nfreq, rad, is, iconv, res, f_out, rad_out, tb_out, nc, nview)
    do i = 1+nc, nfreq-nc
        rad(is,i) = rad_out(i)
        tb(is,i)  = tb_out(i)
    end do
end do
do i = 1+nc, nfreq-nc
    write(22,fmt='(f13.6,*(e13.5,f9.3))')f_out(i), (rad(is,i), tb(is,i), is=1, nview)
end do
close(22)
deallocate(tb)

if (mode_inversion /= 0) then
    !* Convolution of the K-matrices
    write(*,*)'Convolution of the K-matrices'
    allocate(rad_k(nfreq))
    rad_k(:) = 0d0
    rad_out(:) = 0d0
    f_out(:)   = 0d0
    do ik = 1, n_k
        do j= 1, nlay
            do is = 1, nview
                do i= 1, nfreq
                    rad_k(i) = matrix_k(ik,j,is,i)
                end do
                call conv_k(freq1, pas_sort, nfreq, rad_k, iconv, res, f_out, rad_out, nc)
                do i = 1+nc, nfreq-nc
                    matrix_k(ik,j,is,i) = rad_out(i)
                end do
            end do
        end do
    end do
    deallocate(rad_k, rad_out)
end if

! mode_inversion = 2: convolution of the matrice_t
if (mode_inversion == 2) then
    write(*,*)'Convolution of the temperature matrices'
    allocate(rad_kb(nfreq),rad_outb(nfreq))
    rad_kb(:) = 0.
    rad_outb(:) = 0.
    do j=1,nlay
        do is=1,nview
            do i=1,nfreq
                rad_kb(i) = matrix_t(j,is,i)
            enddo
            call conv_k(freq1,pas_sort,nfreq,rad_kb,iconv,res,f_out,rad_outb,nc)
            do i=1+nc,nfreq-nc
                matrix_t(j,is,i)=rad_outb(i)
            enddo
        enddo
    enddo
    deallocate(rad_kb,rad_outb)
end if

! Calculation of the weighting functions for airmass # itrans
if(itrans > 0) then
    open(unit=23, file=trim(path_input)//file_trans, status='unknown')
    call conv_wf(freq1, pas_sort, nfreq, etau_total, itrans, iconv, res, ap, nlay, f_out, nview, nlev)
    close(23)
end if
deallocate(f_out)

call cpu_time(time_end)
time_elapsed = time_end - time_begin
write(*,*)'=================================='
write(*,*)' Transfert takes :',time_elapsed,' s'

return
end subroutine transfer_abundance
