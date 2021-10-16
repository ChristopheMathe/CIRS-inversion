!***		Abundance inversion program (n_k = 3 absorbers at most)      ***
!*		Radiative transfer code for Titan
!*	* N2-N2, N2-CH4, N2-H2 and CH4-CH4 continuum
!*             (coefficients read from unit 12: 4 spectroscopic files)
!*	* Multi cloud non-grey model (only absorption: w0=0)
!*      	FREQ1: First frequency (cm-1)
!*		FREQ2: Last frequency  (cm-1)
!***									     ***
subroutine transfer_transac(nlev,nlay,qn2,qh2,qar,pl,tl,ml,gl,qch4,ql,p,ap,t,nview,sec,lay_min,&
ncorps,mass,erot,et,freq1,freq2,pas_mult,pas_sort,file_trans,file_out,iconv,res,itrans,w_in,&
s_in,g_in,e_in,icloud,taucl,smin,ttest,nlor,alor,glor,elor,g2lor,nvib,vib,ndeg,n_k,icorps_k,&
matrix_k,matrix_t,rad,ikh,nfcont,nond,nond_max,sig,qex,nbl_sp_max,nbl_sp,f1_cont,df_cont,iter,nfreq_inv,qref,&
n2n2,n2ch4,ch4ch4,n2h2,v,corps,path_input,limbe)

use lib_functions
implicit none
!~ PARAMETERS
integer, parameter :: nvoigt_max = 10000
double precision, parameter :: c = 2.997925d+08
double precision, parameter :: atm = 1013.25d0
double precision, parameter :: t0 = 273.15d0
character(*) :: path_input

!~ INTENT(IN)
integer :: itrans
integer :: iter
integer :: nbl_sp_max
integer :: nlev
integer :: nlay
integer :: n_k
integer :: nond_max
integer :: nview
integer :: ncorps
integer :: nfcont
integer :: icloud
integer :: nfreq_inv
integer :: ikh

double precision :: df_cont
double precision :: f1_cont
double precision :: qar
double precision :: qn2
double precision :: qh2
double precision :: res
double precision :: pas_sort
double precision :: pas_mult
double precision :: freq2
double precision :: freq1

character(len=1)  :: iconv
character(len=80) :: file_trans
character(len=80) :: file_out

integer, dimension(ncorps,nbl_sp_max) :: n_out
integer, dimension(ncorps)            :: nlines
integer, dimension(ncorps)            :: nbl_sp
integer, dimension(nview)             :: lay_min
integer, dimension(icloud)            :: nond
integer, dimension(ncorps)            :: nvib
integer, dimension(ncorps)            :: nlor
integer, dimension(ncorps,5)          :: ndeg
integer, dimension(n_k)               :: icorps_k

double precision, dimension(0:500)             :: et
double precision, dimension(icloud,nond_max)   :: sig
double precision, dimension(icloud,nond_max)   :: qex
double precision, dimension(icloud)            :: qref
double precision, dimension(ncorps)        :: mass
double precision, dimension(nlay,nfcont)       :: n2n2
double precision, dimension(nlay,nfcont)       :: n2ch4
double precision, dimension(nlay,nfcont)       :: ch4ch4
double precision, dimension(nlay,nfcont)       :: n2h2
double precision, dimension(ncorps)        :: erot
double precision, dimension(ncorps)            :: alor
double precision, dimension(ncorps,5)          :: glor
double precision, dimension(ncorps,5)          :: elor
double precision, dimension(ncorps,5)          :: g2lor
double precision, dimension(ncorps,5)          :: vib
integer, dimension(ncorps,5)          :: n_nl
double precision, dimension(ncorps)            :: smin
double precision, dimension(ncorps)            :: ttest
double precision, dimension(400,100)           :: v
double precision, dimension(ncorps,nbl_sp_max) :: w_in
double precision, dimension(ncorps,nbl_sp_max) :: s_in
double precision, dimension(ncorps,nbl_sp_max) :: g_in
double precision, dimension(ncorps,nbl_sp_max) :: e_in
double precision, dimension(nlev)              :: p
double precision, dimension(nlev)              :: ap
double precision, dimension(nlev)              :: t
double precision, dimension(nlay)              :: pl
double precision, dimension(nlay)              :: tl
double precision, dimension(nlay)              :: gl
double precision, dimension(nlay)              :: ml
double precision, dimension(nlay)              :: qch4
double precision, dimension(icloud,nlay)       :: taucl
double precision, dimension(nlay,nview)         :: sec
double precision, dimension(nlay,0:ncorps)     :: ql

character(len=256), dimension(ncorps)       :: corps

logical :: limbe
!~ INTENT(INOUT)
double precision, dimension(n_k,nlay,nview,nfreq_inv) :: matrix_k     !(3,nlay,nview,nfreq_inv)
double precision, dimension(nlay,nview,nfreq_inv)        :: matrix_t     !(nlay,nview,nfreq_inv)
double precision, dimension(nview,nfreq_inv)             :: rad          !(nview,nfreq_inv) 

!~ LOCAL
integer :: i = 0
integer :: i1 = 0
integer :: ic = 0
integer :: ik = 0
integer :: is = 0
integer :: isort = 0
integer :: icalc = 0
integer :: j = 0
integer :: j1 = 0
integer :: j2 = 0
integer :: l = 0
integer :: k = 0
integer :: navgj = 0
integer :: navg0 = 0
integer :: n_sort1 = 0
integer :: nstep = 0
integer :: nstep0 = 0
integer :: nc = 0
integer :: nvoigt = 0
integer  :: n_sort = 0

double precision :: avge = 0d0
double precision :: anc = 0d0
double precision :: ann = 0d0
double precision :: anh = 0d0
double precision :: acc = 0d0
double precision :: cloudj = 0d0
double precision :: cmam = 0d0
double precision :: dt = 0d0
double precision :: f = 0d0
double precision :: f1 = 0d0
double precision :: f2 = 0d0
double precision :: f_mult = 0d0
double precision :: f_cent = 0d0
double precision :: fdop = 0d0
double precision :: fc = 0d0
double precision :: fac_cont = 0d0
double precision :: h0 = 0d0
double precision :: hw0 = 0d0
double precision :: hck296 = hck/296d0
double precision :: hckt = 0.d0
double precision :: qn = 0d0
double precision :: qh = 0d0
double precision :: qa = 0d0
double precision :: step0 = 0d0
double precision :: stepj = 0d0
double precision :: rlor = 0d0
double precision :: v0 = 0d0
double precision :: pj = 0d0
double precision :: tj = 0d0
double precision :: x = 0d0
double precision :: y = 0d0
double precision :: f_xivib = 0d0
double precision :: f_xi = 0d0
double precision :: time_begin = 0d0                                     ! Begin of the subroutine
double precision :: time_end = 0d0                                       ! End of the subroutine
double precision :: time_elapsed = 0d0                                   ! Time spend of the subroutine (seconds)

integer, dimension(nlay)                        :: nst

double precision, dimension(ncorps,nbl_sp_max)  :: w_out
double precision, dimension(ncorps,nbl_sp_max)  :: s_out
double precision, dimension(ncorps,nbl_sp_max)  :: e_out

double precision, dimension(:)    , allocatable :: dtau
double precision, dimension(:)    , allocatable :: etau00
double precision, dimension(:)    , allocatable :: dtauk
double precision, dimension(:)    , allocatable :: dtauk0
double precision, dimension(:,:)  , allocatable :: etau0
double precision, dimension(:)    , allocatable :: detau0
double precision, dimension(:,:,:), allocatable :: etau_k               !(nlay,nlay,n_sort)
double precision, dimension(:,:,:), allocatable :: etau_k2              !(nlay,nlay,n_sort)
double precision, dimension(:)    , allocatable :: rad_k                !(n_sort)
double precision, dimension(:)    , allocatable :: rad_kb               !(n_sort)
double precision, dimension(:,:,:), allocatable :: etau                 !(nlev,nview,n_sort) 
double precision, dimension(:,:)  , allocatable :: tb                   !(nview,n_sort)
double precision, dimension(:)    , allocatable :: rad_out              !(n_sort)
double precision, dimension(:)    , allocatable :: rad_outb             !(n_sort) 
double precision, dimension(:)    , allocatable :: tb_out               !(n_sort)
double precision, dimension(:)    , allocatable :: f_out                !(n_sort)
double precision, dimension(:,:)  , allocatable :: pl_lay               !(nlay,n_sort) 
double precision, dimension(:)    , allocatable :: pl_ground            !(n_sort)
double precision, dimension(:,:)  , allocatable :: dpl_lay              !(nlay,n_sort)
double precision, dimension(:)    , allocatable :: dpl_ground           !(n_sort)
double precision, dimension(:,:,:), allocatable :: dtau_k               !(3,nlay,nx mais pas sur)
double precision, dimension(:,:,:), allocatable :: etau0_k              
double precision, dimension(:,:)  , allocatable :: vgt                  !(5,nvoigt)
double precision, dimension(:)    , allocatable :: s_over               

!~ =====================================================================
!~ =====================================================================
!~ *** BEGIN TRANSFER ***
!~ =====================================================================
!~ =====================================================================

call cpu_time(time_begin)

f2=0d0
f1=freq1
icalc=0
n_sort1=0

if(limbe .AND. iter == 0) then
    do j=1,nlay
        do is=2,nview
            sec(j,is)=sec(j,is)/sec(j,1)
        end do
    end do
ENDIF

allocate(etau(1,1,1))
DO WHILE (f2 < freq2)
    deallocate(etau)
!** Determination of the calculation step (cm-1)
    icalc=icalc+1
    hw0=1.d+03
    do j=1,nlay
        call det_step(hw0,f1,tl(j),pl(j),mass,ncorps,glor,elor,step0,navg0,pas_mult,pas_sort)
    end do
    nstep0=min0(navg0*(idint((freq2-f1)/pas_sort)+1),navg0*((40000-1)/navg0))+1 !40000 vient de l'ancien nfreq_max = 40 000
    n_sort=(nstep0-1)/navg0
    f2=f1+(nstep0-1)*step0
    write(6,212)icalc,f1,f2,step0
    212 format(//,10X,'Spectral interval #',I3,'   ',//,'Calculation from',F10.4,'to',F10.4,' with a step of',E11.4,' cm-1')
    allocate(etau_k(nlay,nlay,n_sort),etau_k2(nlay,nlay,n_sort),pl_lay(nlay,n_sort),pl_ground(n_sort))
    allocate(dpl_lay(nlay,n_sort),dpl_ground(n_sort))

!** Input: Molecular line parameters
    write(6,*)'INPUT: MOLECULAR LINE PARAMETERS'
    call read_line2(f1,f2,alor,smin,ttest,erot,nlor,g2lor,w_in,s_in,g_in,e_in,&
    w_out,s_out,e_out,n_out,n_nl,nlines,nbl_sp_max,nbl_sp,ncorps)

    do k=1,ncorps
        write(6,238)corps(k),nlines(k),ttest(k),smin(k),erot(k)
        if (nlines(k) > 0)then
            do i=1,nlor(k)
                write(6,240),n_nl(k,i),glor(k,i),elor(k,i)
            end do
            if (nvib(k) > 0) then
                write(6,239),nvib(k),(vib(k,i),ndeg(k,i),i=1,nvib(k))
            end if
        end if
    end do
    238 format(A4,': ',I5,' lines with S(',F5.1,' K) >',E10.3,' n_rot =',F4.1)
    239 format(9X,I2,' vibrational modes: ',6(F8.2,' (',i1,')',1X))
    240 format(4X,I7,' lines with Lorentz halfwidth ⁼',F6.3,' cm-1, ','n = -',F5.3)
    f=0.5*(f1+f2)

!** Calculation of the optical depths in the layers
    allocate(etau(nlev,nview,n_sort),etau0(nview,nstep0),etau0_k(nlev,nview,nstep0))
    allocate(dtauk(nstep0),dtau_k(n_k,nlay,nstep0))

    do 29 is=1,nview
        if(nview /= 0) THEN
                etau(nlev,is,:)=1.d0
                etau0_k(nlev,is,:)=1.d0
                etau0(is,:)=1.d0
        ELSE
                etau(nlev,is,:)=1.d0
                etau0(is,:)=0d0
        ENDIF
    29 END DO

    allocate(dtauk0(nstep0),detau0(nstep0),etau00(nstep0))
    do 18 j=nlay,1,-1
        tj=tl(j)
        hckt=hck/tj
        dt=hckt-hck296
        pj=pl(j)
        h0=1.d+07*r*t0/(ml(j)*gl(j))
        cmam=2.6867d+19*spi*h0*sec(j,1)*(p(j)-p(j+1))/atm
        hw0=1.d+03
        call det_step(hw0,f,tj,pj,mass,ncorps,glor,elor,stepj,navgj,pas_mult,pas_sort)
        nstep=navgj*n_sort+1
        nst(j)=nstep
        !*		Opacity due to N2-N2, CH4-CH4, N2-CH4, N2-H2
        fac_cont=(t0/tj)*h0*pj*(p(j)-p(j+1))/(atm*atm)
        qn=(1.d0-qch4(j))/(1.d0+(qh2+qar)/qn2)
        qh=(1.d0-qch4(j))/(1.d0+(qn2+qar)/qh2)
        if (qar /= 0) then
            qa=(1.d0-qch4(j))/(1.d0+(qn2+qh2)/qar)
        else
            qa = 0d0
        end if
    !*		Assume N2-N2 and N2-Ar opacities are the same
        qn=qn+qa
        allocate(dtau(nstep))
        do 19 i=1,nstep
            fc=(f1+stepj*dfloat(i-1)-f1_cont)/df_cont
            i1=idint(fc)+1
            fc=fc-i1+1.
            ann=(1.d0-fc)*n2n2(j,i1)+fc*n2n2(j,i1+1)
            anc=(1.d0-fc)*n2ch4(j,i1)+fc*n2ch4(j,i1+1)
            acc=(1.d0-fc)*ch4ch4(j,i1)+fc*ch4ch4(j,i1+1)
            anh=(1.d0-fc)*n2h2(j,i1)+fc*n2h2(j,i1+1)
            dtau(i)=sec(j,1)*fac_cont*(qn*(qn*ann+qch4(j)*anc+qh*anh)+qch4(j)*qch4(j)*acc)
        19 END DO
    !*		Cloud opacity
        if(icloud == 0) goto 508
        do  ik=1,n_sort
            f_cent=f1+pas_sort*(dfloat(ik)-0.5)
            cloudj=0.
            do ic=1,icloud
                cloudj=cloudj+sec(j,1)*taucl(ic,j)*tina(f_cent,sig(ic,:),qex(ic,:),nond(ic))/qref(ic)
            end do
            i1=(ik-1)*navgj
            do i=i1+1,i1+navgj
                dtau(i)=dtau(i)+cloudj
            enddo
            if(ikh > 0) then
                do i=i1+1,i1+navgj
                    dtau_k(ikh,j,i)=cloudj
                enddo
            endif
        end do
        dtau(nstep)=dtau(nstep)+cloudj

        if(ikh > 0) then
            dtau_k(ikh,j,nstep)=cloudj
        end if

    !*		Iteration over the NCORPS absorbers
        508 do 20 k=1,ncorps
            if(nlines(k) <= 0) goto 20
            f_xi=(296.d0/tj)**erot(k)
            f_xivib=1.d0
            if(nvib(k) > 0) THEN
                do i=1,nvib(k)
                    f_xivib=f_xivib*((1.d0-dexp(-vib(k,i)*hckt))/(1.d0-dexp(-vib(k,i)*hck296)))**ndeg(k,i)
                end do
            ENDIF
            nvoigt = min0(nvoigt_max,1+idnint(alor(k)/stepj))
            allocate(vgt(5,nvoigt))
            dtauk(:) = 0.d0
        !*  --- 1/Doppler parameter (cm-1)
            fdop=1./(f*dsqrt(2.d+03*r*tj/mass(k))/c)
            do 21 ik=nlor(k),1,-1
            !*  --- Lorentz halfwidth (cm-1)
                rlor=glor(k,ik)*(pj/atm)*(296.d0/tj)**elor(k,ik)
                y=fdop*rlor
                v0=voigt(0d0,y,v)
                do 22 i=1,nvoigt
                    x=fdop*stepj*(i-1)
                    vgt(ik,i)=voigt(x,y,v)
                    if(vgt(ik,i) < 0.5d-09*v0 .AND. ik == nlor(k)) then
                        nvoigt=i
                        goto 21
                    endif
                22 END DO
            21 END DO
        !*		Iteration over the NLINES(K) lines of Absorber K
            allocate(s_over(nlines(k)))
            do 23 l=1,nlines(k)
                s_over(l)=s_out(k,l)*dexp(-e_out(k,l)*dt)*(1.d0-dexp(-w_out(k,l)*hckt))/(1.d0-dexp(-w_out(k,l)*hck296))
            23 END DO
            if(nstep < nvoigt) THEN
            !*			(f2-f1) < cutoff
                call tau_lines1(w_out(k,:),s_over,n_out(k,:),nlines(k),dtauk,f1,stepj,nstep,nstep0,vgt,nvoigt)
            ELSE
                if(nstep < 2*nvoigt) then
                !*		  cutoff <= (f2-f1) < 2*cutoff
                    call tau_lines2(w_out(k,:),s_over,n_out(k,:),nlines(k),dtauk,f1,stepj,nstep,nstep0,vgt,nvoigt)
                else
                !*		      (f2-f1) >= 2*cutoff
                    call tau_lines3(w_out(k,:),s_over,n_out(k,:),nlines(k),dtauk,f1,stepj,nstep,nstep0,vgt,nvoigt)
                endif
            ENDIF

            deallocate(vgt,s_over)
            f_mult=fdop*f_xi*f_xivib*cmam
            f_mult=f_mult*ql(j,k)
            dtau= dtau + f_mult*dtauk(:nstep)
            do 60 ik=1,n_k
                if(k == icorps_k(ik)) THEN
                    dtau_k(ik,j,:) = f_mult*dtauk(:)
                ENDIF
            60 END DO
        20 END DO

        do 61 ik=1,n_k
            dtauk(:) = dtau_k(ik,j,:)
            call intwdl(dtauk,dtauk0,nstep,nstep0,stepj,step0)
            dtau_k(ik,j,:)=dtauk0(:)
        61 END DO

        if(nview == 0) THEN
            call dtau_e(dtau,etau0,detau0,nstep,nstep0,stepj,step0,et,nview)
            call avg(detau0,etau,n_sort,nview,nlev,nstep0,navg0,j,1)
        ELSE
            do 31 is=1,nview
                if(.not. limbe .OR. (limbe .AND. lay_min(is) == 0)) then
                    if(is == 1) then
                        dtauk(:nstep) = dtau
                    else
                        dtauk(:nstep) = dtau * sec(j,is)
                    endif
                    call dtau_v(dtauk,etau0,detau0,nstep,nstep0,stepj,step0,et,is,nview)
                    dtauk = etau0(is,:)
                    etau0_k(j,is,:) = etau0(is,:)
                    call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,j,is)
                ELSE
                    if(is == 1) then
                        dtauk(:nstep) = dtau
                        call dtau_h(dtauk,etau0_k,detau0,nstep,nstep0,stepj,step0,et,j,is,nlev,nview)
                    else
                        if(lay_min(is) <= j) then
                            dtauk(:nstep) = dtau * sec(j,is)
                            call dtau_h(dtauk,etau0_k,detau0,nstep,nstep0,stepj,step0,et,j,is,nlev,nview)
                         endif
                    endif
                ENDIF
            31 END DO
        ENDIF
        deallocate(dtau)
    18 END DO
    deallocate(dtauk0,detau0)
    print *, 'test transfer 2'

!*            For horizontal viewing, compute
!*      DEXP(-tau(j)) - DEXP(-2*tau(1))/DEXP(-tau(j))
    if(limbe) then
        do 37 is=1,nview
            if(lay_min(is) == 0) goto 37
            etau00=etau0_k(lay_min(is),is,:)**2
            do j=nlay,lay_min(is),-1
                do i=1,nstep0
                    if(etau00(i) > 0.) then
                        dtauk(i)=etau0_k(j,is,i)-etau00(i)/etau0_k(j,is,i)
                    else
                        dtauk(i)=etau0_k(j,is,i)
                    endif
                end do
                call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,j,is)
            end do
            dtauk = 1.-etau00
            call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,nlev,is)
        37 end do
    end if
    deallocate(etau00)
    
!**	  Calculation of the Planck function (erg s-1 cm-2 sr-1/cm-1)
    call planck(f1,pas_sort,n_sort,nlay,tl,t(1),pl_lay,pl_ground)
    call dplanck(f1,pas_sort,n_sort,nlay,tl,t(1),dpl_lay,dpl_ground)

    print *, 'test transfer 3'

!**		    Calculation of the outgoing radiances
    do 46 is=1,nview
        do 47 i=1,n_sort
            isort=n_sort1+i
            rad(is,isort)=0.
            if(limbe .AND. lay_min(is) > 0) then
                do 480 j=lay_min(is),nlay
                    rad(is,isort)=rad(is,isort)+pl_lay(j,i)*(etau(j+1,is,i)-etau(j,is,i))
                480 END DO
            ELSE
                do 48 j=1,nlay
                    rad(is,isort)=rad(is,isort)+pl_lay(j,i)*(etau(j+1,is,i)-etau(j,is,i))
                48 END DO
                rad(is,isort)=rad(is,isort)+pl_ground(i)*etau(1,is,i)
            ENDIF
        47 END DO
    46 END DO

    print *, 'test transfer 4'
!            --------------------------------------------
!**		Calculation of the matrix_t (for temperature)
    do is=1,nview
        do i=1,n_sort
            isort=n_sort1+i
            if(limbe .AND. lay_min(is) > 0) then
                do j=lay_min(is),nlay
                    matrix_t(j,is,isort)=(etau(j+1,is,i)-etau(j,is,i))*dpl_lay(j,i)
                enddo
                if(lay_min(is) > 1) then
                    do j=1,lay_min(is)-1
                        matrix_t(j,is,isort)=0d0
                    enddo
                endif
            ELSE
                do j=1,nlay
                    matrix_t(j,is,isort)=(etau(j+1,is,i)-etau(j,is,i))*dpl_lay(j,i)
                enddo
                matrix_t(1,is,isort)=matrix_t(1,is,isort)+etau(1,is,i)*dpl_ground(i)
            ENDIF
        enddo
    enddo

    print *, 'test transfer 5'
    
!**		Calculation of the K matrices (for absorbers)
    do ik=1,n_k
        do is=1,nview
            do 71 j1=1,nlay
                do j2=1,j1
                    if(is == 1) then
                        dtauk = dtau_k(ik,j1,:)*etau0_k(j2,is,:)
                    else
                        dtauk = dtau_k(ik,j1,:)*sec(j1,is)*etau0_k(j2,is,:)
                    end if
                    call avg_k(dtauk,etau_k,n_sort,nlay,nstep0,navg0,j1,j2)
                end do
                
                if((.not. limbe) .OR. (limbe .AND. lay_min(is) == 0)) goto 71
                
                do j2=j1+1,nlay
                    if(is == 1) then
                        dtauk = dtau_k(ik,j1,:)*etau0_k(j2,is,:)
                    else
                        dtauk = dtau_k(ik,j1,:)*sec(j1,is)*etau0_k(j2,is,:)
                    end if
                    call avg_k(dtauk,etau_k,n_sort,nlay,nstep0,navg0,j1,j2)
                end do
                
                ! calcul de etau_k2
                do j2=lay_min(is),nlay
                    if(is == 1) then
                        do i=1,nstep0
                            if (etau0_k(lay_min(is),is,i) > 0d0) then
                                dtauk(i) = dtau_k(ik,j1,i)*etau0_k(lay_min(is),is,i)**2/etau0_k(j2,is,i)
                            else
                                dtauk(i) = 0d0
                            endif
                        end do
                    else
                        do i=1,nstep0
                            if (etau0_k(lay_min(is),is,i) > 0d0) then
                                dtauk(i) = dtau_k(ik,j1,i)*sec(j1,is)*etau0_k(lay_min(is),is,i)**2/etau0_k(j2,is,i)
                            else
                                dtauk(i) = 0d0
                            endif
                        end do
                        end if
                    call avg_k(dtauk,etau_k2,n_sort,nlay,nstep0,navg0,j1,j2)
                end do
            71 end do
            ! Fill matrix K
            do i=1,n_sort
                isort=n_sort1+i
                if(.not. limbe .OR. (limbe .AND. lay_min(is) == 0)) then
                    do j1=1,nlay
                        matrix_k(ik,j1,is,isort)=etau_k(j1,j1,i)*pl_lay(j1,i)-etau_k(j1,1,i)*pl_ground(i)
                        if (j1 > 1) then
                            do j2=1,j1-1
                                matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)-(etau_k(j1,j2+1,i)-etau_k(j1,j2,i))*pl_lay(j2,i)
                            end do
                        end if
                    end do
                else
                    do j1=lay_min(is),nlay-1!si pas de -1 je dépasse dim2 de etau_k2
                        matrix_k(ik,j1,is,isort) = pl_lay(j1,i)*(etau_k(j1,j1,i) - etau_k2(j1,j1,i) + 2.*etau_k2(j1,j1+1,i))
                        if(j1 > lay_min(is)) then
                            do j2=lay_min(is),j1-1
                                matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)- &
                                & (etau_k(j1,j2+1,i)-etau_k(j1,j2,i)+etau_k2(j1,j2,i)-etau_k2(j1,j2+1,i))*pl_lay(j2,i)
                            end do
                        end if
                        if(j1 < nlay) then
                            do j2=j1+1,nlay-1
                                matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)-(etau_k2(j1,j2,i)-&
                                & etau_k2(j1,j2+1,i))*2.d0*pl_lay(j2,i)
                            end do
                        ENDIF
                    end do
                    if(lay_min(is) > 1) then
                        do j1=1,lay_min(is)-1
                            matrix_k(ik,j1,is,isort)=0d0
                        end do
                    end if
                end if
            end do
        end do
    end do

    print *, 'test transfer 6'
    !**		Writing the averaged exp(-tau) for Calculation # 1
    write(6,222)
    if(.not. limbe) then
        write(6,221),nlev,p(nlev),1.d0
    else
        avge=0d0
        do i=1,n_sort
            avge = avge + etau(nlev,1,i)
        end do
        write(6,221),nlev,p(nlev),avge/(n_sort*1.d0)
    end if
    do j=nlay,1,-1
        avge=0d0
        do i=1,n_sort
            avge = avge + etau(j,1,i)
        end do
        avge = avge/(n_sort*1.d0)
        write(6,223),j,p(j),avge,pl(j),tl(j),sec(j,1),nst(j)
    end do
    221 format(I4,E14.3,F12.7)
    222 format(/,10X,'P (mb)',3X,'Transmittance',9X,'P (mb)',4X,'T (K)',3X,'Airmass',2X,'N_step')
    223 format(I4,E14.3,F12.7,5X,E14.3,F8.2,F8.4,I8)
    f1=f2
    n_sort1=n_sort1+n_sort
    deallocate(etau_k)
    deallocate(etau_k2)
    deallocate(pl_lay)
    deallocate(pl_ground)
    deallocate(dpl_lay)
    deallocate(dpl_ground)
    deallocate(dtauk)
    deallocate(dtau_k)
    deallocate(etau0_k)
    deallocate(etau0)
ENDDO

print *, 'test transfer 7'
!**	   Convolution of the spectrum by the apparatus function
open(unit=22,file=file_out,status='unknown',form='formatted')
allocate(tb(nview,n_sort1),rad_out(n_sort1),tb_out(n_sort1),f_out(n_sort1))
do 54 is=1,nview
    call conv(freq1,pas_sort,n_sort1,rad,is,iconv,res,f_out,rad_out,tb_out,nc,nview)
    do 55 i=1+nc,n_sort1-nc
        rad(is,i)=rad_out(i)
        tb(is,i)=tb_out(i)
    55 END DO
54 END DO
deallocate(tb_out)
do 56 i=1+nc,n_sort1-nc
    write(22,220) f_out(i),(rad(is,i),tb(is,i),is=1,nview)
56 END DO
220 format(f13.6,*(e13.5,f9.3))
close(22)
deallocate(tb)

print *, 'test transfer 8'
!**			Convolution of the K-matrices
allocate(rad_k(n_sort1))
do 85 ik=1,n_k
    do 86 j=1,nlay
        do 87 is=1,nview
            do 88 i=1,n_sort1
                rad_k(i) = matrix_k(ik,j,is,i)
            88 END DO
            call conv_k(freq1,pas_sort,n_sort1,rad_k,iconv,res,f_out,rad_out,nc)
            do 89 i=1+nc,n_sort1-nc
                matrix_k(ik,j,is,i)=rad_out(i)
            89 END DO
        87  END DO
    86  END DO
85  END DO
deallocate(rad_k,rad_out)
print *, 'test transfer 9'
!**			Convolution of the matrice_t
allocate(rad_kb(n_sort1),rad_outb(n_sort1))
do j=1,nlay
    do is=1,nview
        do i=1,n_sort1
            rad_kb(i) = matrix_t(j,is,i)
        enddo
        call conv_k(freq1,pas_sort,n_sort1,rad_kb,iconv,res,f_out,rad_outb,nc)
        do i=1+nc,n_sort1-nc
            matrix_t(j,is,i)=rad_outb(i)
        enddo
    enddo
enddo
deallocate(rad_kb,rad_outb)
print *, 'test transfer 10'
!**	   Calculation of the weighting functions for airmass # itrans
if(itrans > 0) THEN
    open(unit=23,file=trim(path_input)//file_trans,status='unknown')
    call conv_wf(freq1,pas_sort,n_sort,etau,itrans,iconv,res,ap,nlay,f_out,nview,nlev)
    close(23)
ENDIF
deallocate(f_out,etau)
print *, 'test transfer 11'


call cpu_time(time_end)
time_elapsed = time_end - time_begin
print *,'=================================='
print *,' Transfert takes :',time_elapsed,' s'

return
end subroutine transfer_transac
