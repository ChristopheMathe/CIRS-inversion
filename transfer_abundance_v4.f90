!***        Abundance inversion program      ***
!*      Radiative transfer code for Titan
!*      N2-N2, N2-CH4, N2-H2 and CH4-CH4 continuum
!*      (coefficients read from unit 12: 4 spectroscopic files)
!*      Multi cloud non-grey model (only absorption: w0=0)
!*      FREQ1: First frequency (cm-1)
!*      FREQ2: Last frequency  (cm-1)
!***                                         ***
subroutine transfer_abundance(nlev,nlay,qn2,qh2,qar,pl,tl,ml,gl,qch4,ql,&
p,ap,t,nview,sec,lay_min,nb_mol,mass,erot,et,freq1,freq2,pas_mult,pas_sort,&
file_trans,file_out,iconv,res,itrans,w_in,s_in,g_in,e_in,icloud,taucl,smin,&
ttest,nlor,alor,glor,elor,g2lor,nvib,vib,ndeg,n_k,icorps_k,matrix_k,rad,ikh,&
nfcont,nond,nond_max,sig,qex,nbl_sp_max,nbl_sp,f1_cont,df_cont,iter,nfreq_inv,&
qref,n2n2,n2ch4,ch4ch4,n2h2,v,corps,path_input,limbe)

use lib_functions
implicit none

! PARAMETERS
integer         , parameter :: nvoigt_max = 10000
double precision, parameter :: c = 2.997925d+08
double precision, parameter :: atm = 1013.25d+00
double precision, parameter :: t0 = 273.15d+00

!* INTENT(IN)

!----- INTEGERS
integer :: nlev
integer :: nlay
integer :: nview
integer :: nb_mol
integer :: n_k
integer :: iter
integer :: nond_max
integer :: nfcont
integer :: n_sort
integer :: icloud
integer :: itrans
integer :: nbl_sp_max
integer :: nfreq_inv
integer :: ikh
integer, dimension(nlay)              :: nst
integer, dimension(icloud)            :: nond
integer, dimension(nb_mol,nbl_sp_max) :: n_out
integer, dimension(nb_mol)            :: nlines
integer, dimension(nview)              :: lay_min
integer, dimension(nb_mol)            :: nvib
integer, dimension(nb_mol)            :: nlor
integer, dimension(nb_mol,5)          :: ndeg
integer, dimension(n_k)               :: icorps_k
integer, dimension(nb_mol,5)          :: n_nl
integer, dimension(nb_mol)            :: nbl_sp
!----- DOUBLE PRECISION
double precision :: freq1
double precision :: freq2
double precision :: qn2
double precision :: qh2
double precision :: qar
double precision :: anc
double precision :: ann
double precision :: anh
double precision :: acc
double precision :: cloudj
double precision :: avge
double precision :: res
double precision :: pas_mult
double precision :: pas_sort
double precision :: df_cont
double precision :: f1_cont
double precision :: f_xi
double precision :: f_mult
double precision :: f_cent
double precision :: f
double precision :: f1
double precision :: f2
double precision :: fdop
double precision :: fc
double precision :: fac_cont
double precision :: f_xivib
double precision :: dt
double precision :: cmam
double precision :: hw0
double precision :: hckt
double precision :: hck296
double precision :: h0
double precision :: pj
double precision :: tj
double precision, dimension(nlay)              :: ml
double precision, dimension(nb_mol)        :: mass
double precision, dimension(nlay,nfcont)       :: n2n2
double precision, dimension(nlay,nfcont)       :: n2ch4
double precision, dimension(nlay,nfcont)       :: ch4ch4
double precision, dimension(nlay,nfcont)       :: n2h2
double precision, dimension(nb_mol)        :: erot
double precision, dimension(nb_mol)            :: alor
double precision, dimension(nb_mol,5)          :: glor
double precision, dimension(nb_mol,5)          :: elor
double precision, dimension(nb_mol,5)          :: g2lor
double precision, dimension(nb_mol,5)          :: vib
double precision, dimension(nb_mol)            :: smin
double precision, dimension(nb_mol)            :: ttest
double precision, dimension(400,100)           :: v
double precision, dimension(nlev)              :: p
double precision, dimension(nlev)              :: ap
double precision, dimension(nlev)              :: t
double precision, dimension(nlay)              :: pl
double precision, dimension(nlay)              :: tl
double precision, dimension(nlay)              :: gl
double precision, dimension(nlay)              :: qch4
double precision, dimension(icloud,nlay)       :: taucl
double precision, dimension(nlay,nview)         :: sec
double precision, dimension(nlay,0:nb_mol)     :: ql
double precision, dimension(0:500)             :: et
double precision, dimension(icloud,nond_max)   :: sig
double precision, dimension(icloud,nond_max)   :: qex
double precision, dimension(nb_mol,nbl_sp_max) :: w_in
double precision, dimension(nb_mol,nbl_sp_max) :: s_in
double precision, dimension(nb_mol,nbl_sp_max) :: g_in
double precision, dimension(nb_mol,nbl_sp_max) :: e_in
double precision, dimension(icloud)             :: qref

!----- CHARACTERS
character(len=1)  :: iconv
character(len=256) :: file_out
character(len=256) :: file_trans
character(*)      :: path_input
character(len=256), dimension(nb_mol) :: corps

!* INTENT (OUT)

!----- DOUBLE PRECISION
double precision, dimension(nview,nfreq_inv)             :: rad
double precision, dimension(n_k,nlay,nview,nfreq_inv) :: matrix_k 

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
integer :: nc = 0
integer :: n_sort1 = 0
integer :: navg0 = 0
integer :: navgj = 0

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
double precision, dimension(nb_mol,nbl_sp_max)  :: w_out
double precision, dimension(nb_mol,nbl_sp_max)  :: s_out
double precision, dimension(nb_mol,nbl_sp_max)  :: e_out
double precision, dimension(:)    , allocatable :: s_over
double precision, dimension(:)    , allocatable :: dtau                 !(nfreq)
double precision, dimension(:)    , allocatable :: etau00               !(nfreq)
double precision, dimension(:)    , allocatable :: dtauk                !(nfreq)
double precision, dimension(:)    , allocatable :: detau0               !(nfreq)
double precision, dimension(:)    , allocatable :: rad_out              !(n_sort)
double precision, dimension(:)    , allocatable :: tb_out               !(n_sort)
double precision, dimension(:)    , allocatable :: f_out                !(n_sort)
double precision, dimension(:)    , allocatable :: rad_k                !(n_sort)
double precision, dimension(:)    , allocatable :: pl_ground            !(n_sort)
double precision, dimension(:)    , allocatable :: dtauk0               !(nfreq)
double precision, dimension(:,:)  , allocatable :: vgt
double precision, dimension(:,:)  , allocatable :: etau0                !(nview,nfreq)
double precision, dimension(:,:)  , allocatable :: tb                   !(n_sort)
double precision, dimension(:,:)  , allocatable :: pl_lay               !(nlay,n_sort)
double precision, dimension(:,:,:), allocatable :: etau                 !(nlev,nview,n_sort)
double precision, dimension(:,:,:), allocatable :: dtau_k               !(3,nlay,nfreq)
double precision, dimension(:,:,:), allocatable :: etau0_k              !(nlev,nview,nfreq)
double precision, dimension(:,:,:), allocatable :: etau_k               !(nlay,nlay,n_sort)
double precision, dimension(:,:,:), allocatable :: etau_k2              !(nlay,nlay,n_sort)

logical, intent(in) :: limbe
!~ =====================================================================
!~ *** BEGIN TRANSFER ***
!~ =====================================================================

call cpu_time(time_begin)
write(*,*)"Begin transfert"

f2=0.d+00
f1=freq1
n_sort1=0
icalc=0
!~ if(nsec > 1 .AND. iter == 0) THEN
if(limbe .AND. iter == 0) then
    do 27 j=1,nlay
        do 27 is=2,nview
            sec(j,is)=sec(j,is)/sec(j,1)
    27 end do
end if
allocate(etau(1,1,1))
DO WHILE (F2 < FREQ2)
    deallocate(etau)
!**		Determination of the calculation step (cm-1)
    icalc=icalc+1
    hw0=1.d+03
    do 17 j=1,nlay
        call det_step(hw0,f1,tl(j),pl(j),mass,nb_mol,glor,elor,step0,navg0,pas_mult,pas_sort)
    17 end do
    nstep0=min0(navg0*(idint((freq2-f1)/pas_sort)+1),navg0*((40000-1)/navg0))+1
    n_sort=(nstep0-1)/navg0

    f2=f1+(nstep0-1)*step0
    
    write(6,212) icalc,f1,f2,step0
    212 format(//,10x,16(1H*),'   Spectral Interval #',i3,'   ',16(1H*),&
    //,' Calculation from',f10.4,' to',f10.4,' with a step of',e11.4,' cm-1')
    allocate(etau_k(nlay,nlay,n_sort),etau_k2(nlay,nlay,n_sort),pl_lay(nlay,n_sort),pl_ground(n_sort))

!**             Input: Molecular line parameters
    write(6,*)'INPUT: MOLECULAR LINE PARAMETERS'
    call read_line2(f1,f2,alor,smin,ttest,erot,nlor,g2lor,w_in,s_in,g_in,e_in,&
    w_out,s_out,e_out,n_out,n_nl,nlines,nbl_sp_max,nbl_sp,nb_mol)

    do k=1,nb_mol
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

    hck296=hck/296.d0
!**		Calculation of the optical depths in the layers
    allocate(etau(nlev,nview,n_sort),etau0(nview,nstep0),etau0_k(nlev,nview,nstep0),dtau_k(n_k,nlay,nstep0))
    do 29 is=1,nview
        if(nview /= 0) then
            etau(nlev,is,:)=1.d+00
            etau0_k(nlev,is,:)=1.d+00
            etau0(is,:)=1.d+00
        else
            etau(nlev,is,:)=1.d+00
            etau0(is,:)=0.d+00
        end if
    29 end do
    allocate(dtauk0(nstep0),etau00(nstep0),dtauk(nstep0))

    do 18 j=nlay,1,-1 !nlay
        allocate(detau0(nstep0))
        detau0=0d0
        tj=tl(j)
        hckt=hck/tj
        dt=hckt-hck296
        pj=pl(j)
        h0=1.d+07*r*t0/(ml(j)*gl(j))
        cmam=2.6867d+19*spi*h0*sec(j,1)*(p(j)-p(j+1))/atm
        hw0=1.d+03
        call det_step(hw0,f,tj,pj,mass,nb_mol,glor,elor,stepj,navgj,pas_mult,pas_sort)
        nstep=navgj*n_sort+1
        nst(j)=nstep
    !*		Opacity due to N2-N2, CH4-CH4, N2-CH4, N2-H2
        fac_cont=(t0/tj)*h0*pj*(p(j)-p(j+1))/(atm*atm)
        qn=(1.d0-qch4(j))/(1.d0+(qh2+qar)/qn2)
        qh=(1.d0-qch4(j))/(1.d0+(qn2+qar)/qh2)
        if (qar /= 0) then
            qa=(1.d+00-qch4(j))/(1.d+00+(qn2+qh2)/qar)
        else
            qa = 0.d+00
        end if
    !*		Assume N2-N2 and N2-Ar opacities are the same
        qn=qn+qa
        allocate(dtau(nstep))
        do 19 i=1,nstep
            fc=(f1+stepj*dfloat(i-1)-f1_cont)/df_cont
            i1=idint(fc)+1
            fc=fc-i1+1.d+00
            ann=(1.d+00-fc)*n2n2(j,i1)+fc*n2n2(j,i1+1)
            anc=(1.d+00-fc)*n2ch4(j,i1)+fc*n2ch4(j,i1+1)
            acc=(1.d+00-fc)*ch4ch4(j,i1)+fc*ch4ch4(j,i1+1)
            anh=(1.d+00-fc)*n2h2(j,i1)+fc*n2h2(j,i1+1)
            dtau(i)=sec(j,1)*fac_cont*(qn*(qn*ann+qch4(j)*anc+qh*anh)+qch4(j)*qch4(j)*acc)
        19 END DO

    !*		Cloud opacity
        if(icloud == 0) goto 508
        do 505 ik=1,n_sort
            f_cent=f1+pas_sort*(dfloat(ik)-0.5d+00)
            cloudj=0d0
            do 506 ic=1,icloud
                cloudj=cloudj+sec(j,1)*taucl(ic,j)*tina(f_cent,sig(ic,:),qex(ic,:),nond(ic))/qref(ic)
            506 END DO
            i1=(ik-1)*navgj
            do i=i1+1,i1+navgj
                dtau(i)=dtau(i)+cloudj
            enddo
            if(ikh > 0) then
                do i=i1+1,i1+navgj
                    dtau_k(ikh,j,i)=cloudj
                enddo
            endif
        505 END DO
        dtau(nstep)=dtau(nstep)+cloudj
        if(ikh > 0) dtau_k(ikh,j,nstep)=cloudj

    !*		Iteration over the nb_mol absorbers
        508 do 20 k=1,nb_mol 
                if(nlines(k) <= 0) goto 20
                    f_xi=(296.d+00/tj)**erot(k)
                    f_xivib=1.
                if(nvib(k) > 0) then
                    do 26 i=1,nvib(k)
                        f_xivib=f_xivib*((1.d0-dexp(-vib(k,i)*hckt))/(1.d0-dexp(-vib(k,i)*hck296)))**ndeg(k,i)
                    26 end do
                end if
                nvoigt=min0(nvoigt_max,1+idnint(alor(k)/stepj))
                allocate(vgt(5,nvoigt))
                do 24 i=1,nstep0
                    dtauk(i)=0.d+00
                24 end do
            !*  --- 1/Doppler parameter (cm-1)
                fdop=1./(f*dsqrt(2.d+03*r*tj/mass(k))/c)
                do 21 ik=nlor(k),1,-1
                !*  --- Lorentz halfwidth (cm-1)
                    rlor=glor(k,ik)*(pj/atm)*(296.d+00/tj)**elor(k,ik)
                    y=fdop*rlor
                    v0=voigt(0.d+00,y,v)
                    do 22 i=1,nvoigt
                        x=fdop*stepj*(i-1)
                        vgt(ik,i)=voigt(x,y,v)
                        if(vgt(ik,i) < 0.5d-09*v0 .AND. ik == nlor(k)) then
                            nvoigt=i
                            goto 21
                        endif
                    22 end do
                21 end do
            !*		Iteration over the NLINES(K) lines of Absorber K
                allocate(s_over(nlines(k)))
                do 23 l=1,nlines(k)
                    s_over(l)=s_out(k,l)*dexp(-e_out(k,l)*dt)*(1.d0-dexp(-w_out(k,l)*hckt))/(1.d0-dexp(-w_out(k,l)*hck296))
                23 end do
                if(nstep < nvoigt) then
                !*			(f2-f1) < cutoff
                    call tau_lines1(w_out(k,:),s_over,n_out(k,:),nlines(k),dtauk,f1,stepj,nstep,nstep0,vgt,nvoigt)
                else
                    if(nstep < 2*nvoigt) then
                    !*		  cutoff <= (f2-f1) < 2*cutoff
                        call tau_lines2(w_out(k,:),s_over,n_out(k,:),nlines(k),dtauk,f1,stepj,nstep,nstep0,vgt,nvoigt)
                    else
                    !*		      (f2-f1) >= 2*cutoff
!~                         print*,'k',k,'dtauk',dtauk(4677)
                        call tau_lines3(w_out(k,:),s_over,n_out(k,:),nlines(k),dtauk,f1,stepj,nstep,nstep0,vgt,nvoigt)
!~                         if(j==39)then
!~                             print*,'j=39',dtauk(4677)
!~                         stop
!~                         end if
                    endif
                end if
                
                deallocate(vgt,s_over)
                
                f_mult=fdop*f_xi*f_xivib*cmam
                f_mult=f_mult*ql(j,k)
                dtau=dtau+f_mult*dtauk(:nstep)
                do 60 ik=1,n_k
                    if(k == icorps_k(ik)) then
                    dtau_k(ik,j,:) = f_mult*dtauk(:)
                    end if
                60 end do
        20 end do
    !**
        do 61 ik=1,n_k
            dtauk(:)=dtau_k(ik,j,:)
            call intwdl(dtauk,dtauk0,nstep,nstep0,stepj,step0)
            dtau_k(ik,j,:)=dtauk0(:)
        61 END DO
        if(nview == 0) then
            call dtau_e(dtau,etau0,detau0,nstep,nstep0,stepj,step0,et,nview)
            call avg(detau0,etau,n_sort,nview,nlev,nstep0,navg0,j,1)
        else
            do 31 is=1,nview
                if(.not. limbe .OR. (limbe .AND. lay_min(is) == 0)) then
                    if(is == 1) then
                        dtauk(:nstep)=dtau
                    else
                        dtauk(:nstep)=dtau *sec(j,is)
                    endif
                    call dtau_v(dtauk,etau0,detau0,nstep,nstep0,stepj,step0,et,is,nview)
                        dtauk = etau0(is,:)
                        etau0_k(j,is,:)=etau0(is,:)
                    call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,j,is)
                else
                    if(is == 1) then
!~                         print*,dtau(4676:4677),nstep
                        dtauk(:nstep)=dtau(:)
!~                         print*,dtau(5248),j
                        call dtau_h(dtauk,etau0_k,detau0,nstep,nstep0,stepj,step0,et,j,is,nlev,nview)
!~                         print*,etau0_k(1,1,5248)
                    else
                        if(lay_min(is) <= j) then
                            dtauk(:nstep)=dtau(:)*sec(j,is)
                            call dtau_h(dtauk,etau0_k,detau0,nstep,nstep0,stepj,step0,et,j,is,nlev,nview)
                        endif
                    endif
                end if
            31 end do
        end if
        deallocate(dtau,detau0)
    18 end do
    deallocate(dtauk0)
    !  For horizontal viewing, compute
    !  DEXP(-tau(j)) - DEXP(-2*tau(1))/DEXP(-tau(j))
    if(limbe) then
        do 37 is=1,nview
            if(lay_min(is) == 0) goto 37
            etau00(:)=etau0_k(lay_min(is),is,:)**2
!~             print*,etau00(5248)
!~             stop
            do 44 j=nlay,lay_min(is),-1
                do 39 i=1,nstep0
                    if(etau00(i) > 0d0) then
                        dtauk(i)=etau0_k(j,is,i)-etau00(i)/etau0_k(j,is,i)
                    else
                        dtauk(i)=etau0_k(j,is,i)
                    endif
!~                     if (j==26 .AND. is==1 .AND. i==5248) then
!~                         print*,i,etau0_k(j,is,i),etau00(i)
!~                     end if
                39 end do
!~                 print*,'titi',j,is,dtauk(5248)
                call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,j,is)
            44 end do
            dtauk(:)= 1. - etau00(:)
            call avg(dtauk,etau,n_sort,nview,nlev,nstep0,navg0,nlev,is)
        37 end do
    end if
!~     print*,etau(1,1,9)
!~     stop
    deallocate(etau00)
!**	  Calculation of the Planck function (erg s-1 cm-2 sr-1/cm-1)
    print*,'Calculation of the Planck function'
    call planck(f1,pas_sort,n_sort,nlay,tl,t(1),pl_lay,pl_ground)
!**		    Calculation of the outgoing radiances
    print*,'Calculation of the outgoing radiances'
    do 46 is=1,nview
        do 47 i=1,n_sort
            isort=n_sort1+i
            rad(is,isort)=0.d0
            if(limbe .AND. lay_min(is) > 0) then
                do 480 j=lay_min(is),nlay
                    rad(is,isort)=rad(is,isort)+pl_lay(j,i)*(etau(j+1,is,i)-etau(j,is,i))
                480 end do
            else
                do 48 j=1,nlay
                    rad(is,isort)=rad(is,isort)+pl_lay(j,i)*(etau(j+1,is,i)-etau(j,is,i))
                48 end do
                rad(is,isort)=rad(is,isort)+pl_ground(i)*etau(1,is,i)
            end if
        47 end do
    46 end do
!**		Calculation of the K matrices
    print*,'Calculation of the K matrices'
    do 69 ik=1,n_k
        print*,'ik: ',ik
        do 70 is=1,nview
            do 71 j1=1,nlay
                do 72 j2=1,j1
                    if(is == 1) then
                        dtauk = dtau_k(ik,j1,:)*etau0_k(j2,is,:)
                    else
                        dtauk = dtau_k(ik,j1,:)*sec(j1,is)*etau0_k(j2,is,:)
                    end if
                    call avg_k(dtauk,etau_k,n_sort,nlay,nstep0,navg0,j1,j2)
                72 end do
                
                if((.not. limbe) .OR. (limbe .AND. lay_min(is) == 0)) goto 71
                
                do 700 j2=j1+1,nlay
                    if(is == 1) then
                        dtauk = dtau_k(ik,j1,:)*etau0_k(j2,is,:)
                    else
                        dtauk = dtau_k(ik,j1,:)*sec(j1,is)*etau0_k(j2,is,:)
                    end if
                    call avg_k(dtauk,etau_k,n_sort,nlay,nstep0,navg0,j1,j2)
                700 END DO
                do 720 j2=lay_min(is),nlay
                    if(is == 1) THEN
                        do 730 i=1,nstep0
                            if (etau0_k(lay_min(is),is,i) > 0d0) then
                                dtauk(i)= dtau_k(ik,j1,i) * etau0_k(lay_min(is),is,i)**2/etau0_k(j2,is,i)
                            else
                                dtauk(i)=0.d+00
                            endif
                        730 END DO
                    ELSE
                        do 740 i=1,nstep0
                            if (etau0_k(lay_min(is),is,i) > 0d0) then
                                dtauk(i) = dtau_k(ik,j1,i)* sec(j1,is)* etau0_k(lay_min(is),is,i)**2/etau0_k(j2,is,i)
                            else
                                dtauk(i)=0.d+00
                            endif
                        740 end do
                    end if
                    call avg_k(dtauk,etau_k2,n_sort,nlay,nstep0,navg0,j1,j2)
                720 end do
            71 end do
!           Fill matrix K
            do 76 i=1,n_sort
                isort=n_sort1+i
                if(.not. limbe .OR. (limbe .AND. lay_min(is) == 0)) then
                    do 75 j1=1,nlay
                        matrix_k(ik,j1,is,isort)=etau_k(j1,j1,i)*pl_lay(j1,i)-etau_k(j1,1,i)*pl_ground(i)
                        if(j1 > 1) then
                            do 77 j2=1,j1-1
                                matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)-(etau_k(j1,j2+1,i)-etau_k(j1,j2,i))*pl_lay(j2,i)
                            77 end do
                        end if
                    75 end do
                ELSE
                    do 750 j1=lay_min(is),nlay -1!si pas de -1 je dépasse dim2 de etau_k2
!~                         if (j1 .lt. nlay) then
                        matrix_k(ik,j1,is,isort)=pl_lay(j1,i)*(etau_k(j1,j1,i)-etau_k2(j1,j1,i)+2.d+00*etau_k2(j1,j1+1,i))
!~                         else
!~                             matrix_k(ik,j1,is,isort)=pl_lay(j1,i)*(etau_k(j1,j1,i)-etau_k2(j1,j1,i))
!~                         end if
                        if(j1 > lay_min(is)) then
                            do 770 j2=lay_min(is),j1-1
                                matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)- &
                                (etau_k(j1,j2+1,i)-etau_k(j1,j2,i)+etau_k2(j1,j2,i)-etau_k2(j1,j2+1,i))*pl_lay(j2,i)
                            770 end do
                        end if
                        if(j1 < nlay) then
                            do 771 j2=j1+1,nlay-1
                                matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)-(etau_k2(j1,j2,i)&
                                -etau_k2(j1,j2+1,i))*2.d+00*pl_lay(j2,i)
                            771 end do
!~                         else
!~                             do 772 j2=j1+1,nlay
!~                                 matrix_k(ik,j1,is,isort)=matrix_k(ik,j1,is,isort)-etau_k2(j1,j2,i)*2.d+00*pl_lay(j2,i)
!~                             772 end do
                        end if
                    750 end do
                    if(lay_min(is) > 1) then
                        do 751 j1=1,lay_min(is)-1
                            matrix_k(ik,j1,is,isort)=0.d+00
                        751 end do
                    end if
                end if
            76 end do
        70 end do
    69 end do
    !    ici le writing exp(-tau)
    !**		Writing the averaged exp(-tau) for Calculation # 1
    write(6,222)
    if(.not. limbe) then
        write(6,221) nlev,p(nlev),1.
    ELSE
        avge=0.
        do 530 i=1,n_sort
            avge=avge+etau(nlev,1,i)
        530 END DO
        write(6,221) nlev,p(nlev),avge/(n_sort*1.d+00)
    ENDIF
    do 52 j=nlay,1,-1
        avge=0.
        do 53 i=1,n_sort
            avge=avge+etau(j,1,i)
        53 END DO
        avge = avge/(n_sort*1.0)
        write(6,221) j,p(j),avge,pl(j),tl(j),sec(j,1),nst(j)
    52 END DO
    221 format(i4,e14.3,f12.7,5x,e14.3,f8.2,f8.4,i8)
    222 format(/,10x,'P (mb)',3x,'Transmittance',9x,'P (mb)',4x,'T (K)', &
    & 3x,'Airmass',2x,'N_step')

    f1=f2
    n_sort1=n_sort1+n_sort
    deallocate(etau_k)
    deallocate(etau_k2)
    deallocate(pl_lay)
    deallocate(pl_ground)
    deallocate(dtauk)
    deallocate(dtau_k)
    deallocate(etau0_k)
    deallocate(etau0)
end do

!**	   Convolution of the spectrum by the apparatus function
!      open(unit=22,file=file_out,status='unknown',recl=231)
print*,rad(1,5),rad(1,7),rad(1,9)
open(unit=22,file=file_out,status='unknown',form='formatted')
allocate(tb(nview,n_sort1),rad_out(n_sort1),tb_out(n_sort1),f_out(n_sort1))
do 54 is=1,nview
    call conv(freq1,pas_sort,n_sort1,rad,is,iconv,res,f_out,rad_out,tb_out,nc,nview)
    do 55 i=1+nc,n_sort1-nc
        rad(is,i)=rad_out(i)
        tb(is,i)=tb_out(i)
    55 end do
54 end do
do 56 i=1+nc,n_sort1-nc
    write(22,220) f_out(i),(rad(is,i),tb(is,i),is=1,nview)
56 end do
220 format(f13.6,*(e13.5,f9.3))
close(22)
deallocate(tb)
!~ print*,rad(1,5),rad(1,7),rad(1,9)
!~ stop
!**			Convolution of the K-matrices
allocate(rad_k(n_sort1))
do 85 ik=1,n_k
    do 86 j=1,nlay
        do 87 is=1,nview
            do 88 i=1,n_sort1
                rad_k(i)=matrix_k(ik,j,is,i)
            88 end do
            call conv_k(freq1,pas_sort,n_sort1,rad_k,iconv,res,f_out,rad_out,nc)
            do 89 i=1+nc,n_sort1-nc
                matrix_k(ik,j,is,i)=rad_out(i)
            89 end do
        87 end do
    86 end do
85 end do
deallocate(rad_k,rad_out)

!**	   Calculation of the weighting functions for airmass # itrans
if(itrans > 0) then
    open(unit=23,file=trim(path_input)//file_trans,status='unknown')
    call conv_wf(freq1,pas_sort,n_sort,etau,itrans,iconv,res,ap,nlay,f_out,nview,nlev)
    close(23)
end if
deallocate(f_out,etau)

call cpu_time(time_end)
time_elapsed = time_end - time_begin
print *,'=================================='
print *,' Transfert takes :',time_elapsed,' s'

return
end subroutine transfer_abundance
