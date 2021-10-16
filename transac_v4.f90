!~ =====================================================================
!~ ===== ABUNDANCE PROFILE INVERSION - TITAN (TIT)                  ====
!~ =====================================================================
!*  nb_mol = Number of molecular absorbers
!*  NLEV   = Number of atmospheric levels
!*     Z(NLEV) : Altitude levels (numbered from 1 - lower boundary -
!*               to NLEV - Upper boundary -). Unit: km
!*               (z=0 at the 1-bar pressure level)
!*     P(NLEV) : Pressure levels (numbered from 1 - lower boundary -
!*               to NLEV - Upper boundary -). Unit: mbar.
!*     AP(NLEV): Log(e) of pressure of the levels
!*     T(NLEV) : Temperature of the levels
!*  NLAY   = NLEV-1 Number of atmospheric layers
!*                  (Layer i is located between levels i and i+1)
!*     PL(NLAY) : Pressure of the layers (mbar)
!*     APL(NLAY): Log(e) of pressure of the layers
!*     TL(NLAY) : Temperature of the layers (K)
!*     SEC(NLAY,nb_obs): Airmasses in the layers
!*     ML(NLAY) : Molecular weight (g).
!*     GL(NLAY) : Gravity (cm sec-2)
!*     QCH4(NLAY): CH4 mole fraction
!*     QL(NLAY,nb_mol) : Absorber Mole fractions
!*  NINPUT = Number of input points for temperature profile
!*     ppi(NINPUT) : Pressure of the input points (mbar)
!*     TI(NINPUT) : Temperature of the input points (K)
!*     API(NINPUT): Log(e) of pressure of the input points
!*     ATI(NINPUT): Log(e) of temperature of the input points
!*  NIN = Number of input points for concentration profiles
!*     PQI(NIN) : Pressure at the input points
!*     APQI(NIN): Log(e) of pressure of the input points
!*     QI(NIN)  : Mole fraction at the input points
!*     AQI(NIN) : Log(e) of mole fraction at the input points
!***                                                                    
Program TRANSAC
use lib_functions

implicit none
!----------------------------------------------------------------------!
! PARAMETERS
!----------------------------------------------------------------------!
character(*), parameter :: &
    path_input = "/home/cmathe/Documents/codes/inversion_v4/data_input/"

logical, parameter :: &
    focal_plane = .False.

!----------------------------------------------------------------------!
! INTEGERS
!----------------------------------------------------------------------!
integer :: &
    i          = 0, &
    ik         = 0, &
    ikh        = 0, &
    k          = 0, &
    iter       = 0, &
    niter      = 0, &
    n1         = 0, &
    n2         = 0, &
    j          = 0, &
    i1         = 0, &
    ic         = 0, &
    ifreq      = 0, &
    ifk        = 0, &
    nb_cloud   = 0, &
    m          = 0, &
    lt         = 0, &
    knu        = 0, &
    nbl        = 0, &
    ja         = 0, &
    ix         = 0, &
    is         = 0, &
    iprint     = 0, &
    ninput_cl  = 0, &
    nlay       = 0, &
    ninput     = 0, &
    nfcont     = 0, &
    nbl_sp_max = 0, &
    nbl_noise  = 0, &
    nx         = 0, &
    np         = 0, &
    nond_max   = 0, &
    nsp        = 0, &
    nfreq_inv  = 0, &
    itrans     = 0, &
    cpt        = 0, &
    tmp_nlev   = 0, &
    tmp_nlay   = 0
integer :: nlor_ch4
integer :: nvib_ch4
integer, parameter:: dbl_type = SELECTED_REAL_KIND(15)
real(kind=dbl_type) :: var_expo
real(kind=dbl_type) :: corr
real(kind=dbl_type) :: corr_b
integer           :: NB_MOL = 0
integer           :: NB_MOL_INV = 0
integer           :: NB_OBS = 0
integer           :: NB_DF = 0
integer           :: NB_KERNEL = 0
integer           :: NB_WEIGHTFUNC = 0
integer           :: nlev = 0
integer, dimension(1,5) :: ndeg_ch4

integer, dimension(:)  , allocatable :: &
&    icorps_k, &     !
&    n_x     , &     !
&    nbl_sp  , &     !
&    lay_min , &     !(nobs)
&    nlor    , &     !(nb_mol)
&    nvib    , &     !(nb_mol)
&    nond            !

integer, dimension(:,:), allocatable :: ndeg !(nb_mol,5)
!----------------------------------------------------------------------!
! REALS
!----------------------------------------------------------------------!
real, dimension(:), allocatable :: &
&    tmod, &
&    gnu

real, dimension(:,:), allocatable :: &
&    n2n2_g  , &
&    n2ch4_g , &
&    ch4ch4_g, &
&    n2h2_g

!----------------------------------------------------------------------!
! DOUBLE PRECISION
!----------------------------------------------------------------------!
double precision :: &
&    lat          = 0.0d0, &
&    rms          = 0.0d0, &
&    time_begin   = 0.0d0, &                                   ! Begin of the subroutine
&   time_end     = 0.0d0, &                                   ! End of the subroutine
&    time_elapsed = 0.0d0, &                                   ! Time spend of the subroutine (seconds)
&    a_1bar       = 0.0d0, &
&    air          = 0.0d0, &
&    apl1         = 0.0d0, &
&    ep           = 0.0d0, &
&    h            = 0.0d0, &
&    g_1bar       = 0.0d0, &
&    gg           = 0.0d0, &
&    f1_cont      = 0.0d0, &
&    dteta        = 0.0d0, &
&    dlnq         = 0.0d0, &
&    df_cont      = 0.0d0, &
&    d_1bar       = 0.0d0, &
&    t_1bar       = 0.0d0, &
&    tra_b          = 0.0d0,&
&    dt             = 0.0d0,&
&    corr2        = 0.0d0, &
&    corrT2          = 0.0d0,&
&    r2           = 0.0d0, &
&    r1           = 0.0d0, &
&    r_1bar       = 0.0d0, &
&    qunsat       = 0.0d0, &
&    qsat         = 0.0d0, &
&    plo          = 0.0d0, &
&    teta1        = 0.0d0, &
&    r3           = 0.0d0, &
&    x            = 0.0d0, &
&    tra          = 0.0d0, &
&    planet_radii = 0.0d0, &
&    planet_gravity = 0.0d0, &
&    planet_flatness = 0.0d0, & !aplatissement de la planete
&    planet_grav_field_rotation = 0.0d0, &   !champ de gravite (rotation de la planete)
&    planet_grav_field_momentum = 0.0d0, & !champ de gravite (moment gravite)
&    P_SURF = 0.0d0, &
&    P_TOP = 0.0d0, &
&    F_START = 0.0d0, &
&    F_END = 0.0d0, &
&    F_PAS = 0.0d0, &
&    F_pas_mult = 0.0d0, &
&    RAYLEIGH = 0.0d0, &
&    VMR_TROPO_N2 = 0.0d0, &
&    VMR_TROPO_H2 = 0.0d0, &
&    VMR_TROPO_AR = 0.0d0, &
&    VMR_TROPO_CH4 = 0.0d0
double precision, dimension(0:500) :: et
double precision,dimension(5)      :: F_CLOUD = 0.0d0
double precision,dimension(20)     :: &
&    LOS = 0.0d0, &
&    F_KERNEL =0.0d0, &
&    CORR_LEN = 0.0d0, &
&    WEI = 0.0d0
double precision,dimension(50,3) :: DF = 0.0d0
double precision, dimension(400,100)    :: v
double precision, dimension(1) :: &
&    a_ch4                , &
&    b_ch4                , &
&    mass_ch4             , &
&    expo_rot_ch4         , &
&    cross_section_min_ch4, &
&    temperature_test_ch4
double precision, dimension(1,5) :: &
&    alor_ch4 , &
&    vib_ch4  , &
&    glor_ch4 , &
&    elor_ch4 , &
&    g2lor_ch4
double precision, dimension(:), allocatable :: &
&    mass             , & !
&    a                , & !
&    b                , & !
&   expo_rot         , & !
&    wobs             , & !
&    robs             , & !
&    wnoise           , & !
&    rnoise           , & !
&    z                , & !(nlev_max)
&    qref             , & !
&   p                , & !(nlev_max)
&    ap               , & !(nlev_max)
&    t                , & !(nlev_max)
&    pl               , & !(nlay_max)
&    apl              , & !(nlay_max)
&    tl               , & !(nlay_max)
&    gl               , & !(nlay_max)
&   qch4             , & !(nlay_max)
&    tab_i            , & !(nlay_max)
&    aql              , & !(nlay_max)
&    tmp_p            , & !(nlay_max)
&    tmp_t            , & !(nlay_max)
&    tmp_ap           , & !(tmp_nlev)
&    tmp_z            , & !(tmp_nlev)
&    tmp_pl           , & !(tmp_nlay)
&    tmp_tl           , & !(tmp_nlay)
&    tmp_apl          , & !(tmp_nlay)
&    tmp_qch4         , & !(tmp_nlay)
&    tmp_ml           , & !(tmp_nlay)
&    tmp_gl           , & !(tmp_nlay)
&    tmp_aql          , & !(tmp_nlay)
&    atmp_pl          , & !
&    atmp_aql         , & !
&    ppi               , & !(ninput_max)
&    ti               , & !(ninput_max)
&    api              , & !(ninput_max)
&    ati              , & !(ninput_max)
&    at               , & !(nlev_max)
&    taucl_tot        , & !(nlay_max)
&    alor             , & !(nb_mol_max)
&    cross_section_min, & !(nb_mol_max)
&    temperature_test , & !(nb_mol_max)
&    freq_x           , & !(nf_max)
&    rad_noise        , & !(nf_max)
&    errt             , & !(nlay_max)
&    errtb            , & !(nlay_max)
&    c0               , & !(nx)
&    ml                   !(nlay_max)

double precision, dimension(:,:)     , allocatable :: &
&    tmp_ql    , & !(tmp_nlay,0:nb_mol)
&    sec       , & !(nlay_max,nb_obs_max)
&    ql        , & !(nlay_max,0:nb_mol_max)
&    pcli      , & !(ninput_max)
&    taucli    , & !(ninput_max)
&    apcli     , & !(ninput_max)
&    ataucli   , & !(ninput_max)
&    taucl     , & !(4,nlay_max)
&   w_input   , & !
&    s_input   , & !
&    g_input   , & !
&    e_input   , & !
&    glor      , & !(nb_mol_max,5)
&    elor      , & !(nb_mol_max,5)
&   g2lor     , & !(nb_mol_max,5)
&    vib       , & !(nb_mol_max,5)
&    rad_obs   , & !(nb_obs_max,nf_max)
&    rad_diff  , & !(nb_obs_max,nf_max)
&    sm        , & !(nlay_max,nx_max)
&    s_b       , & !(nlay_max,nlay_max)
&    ssk_b     , & !(nx_max,nx_max)
&    sig       , & !(nlay_max,nlay_max)
&    qex       , & !(nlay_max,nlay_max)
&    cf        , & !(nlay_max,nlay_max)
&    cf_b      , & !(nlay_max,nlay_max)
&    dlnerrqb  , & !(nlay_max,nb_mol_max)
&    errqb_sup , & !(nlay_max,nb_mol_max)
&    errqb_inf , & !(nlay_max,nb_mol_max)
&    dlnerrq   , & !(nlay_max,nb_mol_max)
&    errq_sup  , & !(nlay_max,nb_mol_max)
&    errq_inf  , & !(nlay_max,nb_mol_max)
&    c         , & !(nx,nx)
&    matrix_r  , & !(nlay_max,nx_max)
&    matrix_r2 , & !(nlay_max,nx_max)
&    rad       , & !(nb_obs_max,nfreq_inv)
&    tab_n2n2  , & !
&    tab_n2ch4 , & !
&    tab_ch4ch4, & !
&    tab_n2h2  !
double precision, dimension(:,:,:), allocatable :: s, & !(nlay_max,nlay_max,nb_mol_max)&
    matrix_t ,& !(nlay,nb_obs,nfreq_inv)
    sk        , & !(3,nlay_max,nx_max)
    ssk        !(3,nx_max,nx_max)

double precision, dimension(:,:,:,:), allocatable :: matrix_k !(3,nlay_max,nb_obs_max,nfreq_inv)

!----------------------------------------------------------------------!
! CHARACTERS
!----------------------------------------------------------------------!
character(len=256) :: &
&    filename_input = '', &
&    file_noise = ''

character(len=256) :: file_spectro_ch4
character(len=256), dimension(5) :: file_scat =''
character(len=256), dimension(:), allocatable :: file_spectro
character(len=1)    :: CONVFUNC = ''
character(len=256)  :: NAME_PLA = ''
character(len=256)  :: PROF_T = ''
character(len=256)  :: FILE_N2N2 = ''
character(len=256)  :: FILE_N2CH4 = ''
character(len=256)  :: FILE_CH4CH4 = ''
character(len=256)  :: FILE_N2H2 = ''
character(len=256)  :: FILE_WEIGHTFUNC = ''
character(len=256)  :: FILE_SPEC_SYNTH = ''
character(len=256)  :: FILE_PROFQ = ''
character(len=256)  :: FILE_KERNEL = ''
character(len=256),dimension(5)  :: prof_haze = ''
character(len=4),dimension(20)    :: NAME_MOL_INV = ''
character(len=256),dimension(20,3)  :: NAME_MOL = ''
character(len=256),dimension(20)    :: file_obs = ''

!----------------------------------------------------------------------!
! LOGICALS
!----------------------------------------------------------------------!
logical :: method_invert
logical :: LIMBE

integer :: nb_kernel_t = 0
double precision ::CORR_LEN_T,WEI_T
double precision,dimension(20) :: F_KERNEL_T
character(len=256) :: FILE_PROF_T,FILE_KERNEL_T
!======================================================================!
!===                      BEGIN PROGRAM                             ===!
!======================================================================!
call CPU_TIME(time_begin)
write(*,*)"============================================================"
write(*,*)"BEGIN PROGRAM TRANSAC"
write(*,*)"============================================================"
call GET_COMMAND_ARGUMENT(1,filename_input)
!** read file parameters
write(*,*)"Read file parameters ..."
call read_parameters(filename_input,NAME_PLA,LAT,NLEV,P_SURF,&
P_TOP,PROF_T,NB_CLOUD,F_CLOUD,PROF_HAZE,file_scat,VMR_TROPO_N2,&
VMR_TROPO_H2,VMR_TROPO_AR,VMR_TROPO_CH4,FILE_N2N2,FILE_N2CH4,FILE_CH4CH4,&
FILE_N2H2,NB_MOL,NAME_MOL,NB_OBS,LOS,FILE_OBS,FILE_NOISE,F_START,F_END,F_PAS,F_pas_mult,&
CONVFUNC,RAYLEIGH,NB_WEIGHTFUNC,FILE_WEIGHTFUNC,NB_DF,DF,NB_MOL_INV,NAME_MOL_INV,&
CORR_LEN,WEI,NITER,IPRINT,NB_KERNEL,F_KERNEL,FILE_SPEC_SYNTH,FILE_PROFQ,FILE_KERNEL,&
LIMBE,NSP,NB_KERNEL_T,F_KERNEL_T,CORR_LEN_T,WEI_T,FILE_PROF_T,FILE_KERNEL_T)
write(*,*)"Read file parameters done."
nlay = abs(nlev) - 1
allocate(file_spectro(nb_mol),cross_section_min(nb_mol),temperature_test(nb_mol))
allocate(nlor(nb_mol),alor(nb_mol))
allocate(glor(nb_mol,5),elor(nb_mol,5),g2lor(nb_mol,5))
allocate(nvib(nb_mol))
allocate(vib(nb_mol,5),ndeg(nb_mol,5))
allocate(p(abs(nlev)),ap(abs(nlev)))
allocate(pl(nlay),apl(nlay))
allocate(mass(nb_mol),a(nb_mol),b(nb_mol),expo_rot(nb_mol))
write(*,*)"--------------------------"
write(*,*)"Searching molecules parameters in the database_mol.f90"
write(*,*)"Case: ",(trim(name_mol(ik,1))//' ',ik=1,nb_mol)
call database_mol(nb_mol,name_mol,focal_plane,nlor,nvib,ndeg,a,b,mass,&
expo_rot,cross_section_min,temperature_test,alor,vib,glor,elor,g2lor,&
file_spectro)
write(*,*)"Done."
write(*,*)"--------------------------"
write(*,*)"Searching planet parameters in the database_planet.f90"
write(*,*)"Case: ",trim(name_pla)
call database_planet(name_pla,lat,g_1bar,r_1bar)
write(*,fmt='(TR5,A,F6.2,A)')'G(surface): ',g_1bar,' cm.sec-2'
write(*,fmt='(TR5,A,TR2,F7.2,A)')'R(1bar): ',r_1bar,' km'
write(*,fmt='(TR5,A,TR7,F6.2,A)')'Lat: ',lat,' deg'
write(*,*)"Done."
!======================================================================!
!=== PART 1: CALCULS ATMOSPHERIC PARAMETERS                         ===!
!======================================================================!
write(*,*)"--------------------------"
write(*,*)"Calculs atmospheric parameters"
write(*,*)"--------------------------"
!**             Input: N2, H2, Ar and CH4 tropospheric mixing ratios
write(*,*)"Input tropospheric mixing ratios"
write(*,fmt='(A,TR2,F5.4)')'N2: ',VMR_TROPO_N2
write(*,fmt='(A,TR2,F5.4)')'H2: ',VMR_TROPO_H2
write(*,fmt='(A,TR2,F8.7)')'Ar: ',VMR_TROPO_AR
write(*,fmt='(A,TR1,F5.4)')'CH4: ',VMR_TROPO_CH4
write(*,*)"Done."

!**             Input: Pressure level grid (mbar)
write(*,*)"--------------------------"
write(*,*)"Input pressure level grid (mbar)"
if(nlev > 0) then ! attention pas pris en compte cette partie a cause du p(j) lu en entree
    method_invert = .FALSE.
    do 8 j=1,nlev
        read(5,*) p(j)
        ap(j)=dlog(p(j))
    8 end do
else
    method_invert = .TRUE.
    nlev=-nlev
    p(1)= p_surf
    ap(1)=dlog(p_surf)
    do 80 j=2,nlev
        ep=dfloat(j-1)/dfloat(nlev-1)
        p(j)=p_surf*(p_top/p_surf)**ep
        ap(j)=dlog(p(j))
    80 end do
end if
write(*,*)"Done."

!**       Calculation of the pressure in the layers (mbar)
do j=1,nlay
    pl(j)=dsqrt(p(j)*p(j+1))
    apl(j)=dlog(pl(j))
end do

write(*,*)"--------------------------"
write(*,*)"Calculation of the pressure in the layers (mbar)"
do j=1,nlay
    pl(j)=dsqrt(p(j)*p(j+1))
    apl(j)=dlog(pl(j))
end do
write(*,*)"Done."

write(*,*)"--------------------------"
write(*,*)"Input temperature profile: "
write(*,fmt='(TR5,A,A)')'Target file: ',prof_t
ninput = nb_line(prof_t)
allocate(ppi(ninput),ti(ninput),api(ninput),ati(ninput))
open(unit=10,file=prof_t,status='old')
do j=1,ninput
    read(10,*) ppi(j),ti(j)
end do
close(10)
if (ppi(1) < ppi(ninput)) then !si p(1) = sol et p(ninput) = ciel
    allocate(tmp_p(ninput),tmp_t(ninput))
    tmp_p = ppi
    tmp_t = ti
    do i= 1, ninput
        ppi(i) = tmp_p(ninput+1-i)
        ti(i) = tmp_t(ninput+1-i)
    end do
    deallocate(tmp_p,tmp_t)
end if
api=dlog(ppi)
ati=dlog(ti)
write(*,*)"Done."

write(*,*)"--------------------------"
write(*,*)"Calculation of the temperature in the layers and the levels (K)"
allocate(t(nlev),tl(nlay))
do 14 j=1,nlay
    t(j)=dexp(tina(ap(j),api,ati,ninput))
    tl(j)=dexp(tina(apl(j),api,ati,ninput))
14 end do
t(nlev)=dexp(tina(ap(nlay),api,ati,ninput))
write(*,*)"Done."

write(*,*)"--------------------------"
write(*,*)"Calculation of the CH4 mixing ratio in the layers"
allocate(qch4(nlay))
call database_mol(1,'CH4',focal_plane,nlor_ch4,nvib_ch4,ndeg_ch4,a_ch4,b_ch4,mass_ch4,&
expo_rot_ch4,cross_section_min_ch4,temperature_test_ch4,alor_ch4,vib_ch4,glor_ch4,elor_ch4,g2lor_ch4,&
file_spectro_ch4)
qsat=1.01325d+03*10**(a_ch4(1)-b_ch4(1)/tl(1))/pl(1)
qch4(1)=dmin1(VMR_TROPO_CH4,qsat)
do 15 j=2,nlay
    qsat=1.01325d+03*10**(a_ch4(1)-b_ch4(1)/tl(j))/pl(j)
    qch4(j)=dmin1(qch4(j-1),qsat)
15 end do
write(*,*)"Done."

write(*,*)"--------------------------"
write(*,*)"Calculation of the molecular weight in the layers"
allocate(ml(nlay))
do 16 j=1,nlay
    ml(j)=16.043d0*qch4(j)+(28.0134d0*VMR_TROPO_N2+2.0158d0*VMR_TROPO_H2+35.9675d0*VMR_TROPO_AR)/(1.d0-VMR_TROPO_CH4+qch4(j))
16 end do
write(*,*)"Done."

write(*,*)"--------------------------"
write(*,*)"Calculation of the altitude levels (km) and of the gravity in the layers (cm sec-2)"

allocate(tab_i(nlay))

do 17 j=1,nlay
    tab_i(j)=j
17 end do

a_1bar=dlog(p(1))
d_1bar=tina(a_1bar,apl,tab_i,nlay)
t_1bar=tina(a_1bar,apl,tl,nlay)

if(d_1bar < 0.4999d+00) print 205
205 format(' WARNING: the 1-bar level is not included in the grid')
i1=idnint(d_1bar+0.0001)

allocate(z(nlev),gl(nlay))

if(i1 > 1) then
    h=1.d+02*r*0.5*(t_1bar+tl(i1))/(ml(i1)*g_1bar)
    z(i1)=h*(a_1bar-ap(i1))
    do 19 j=i1,2,-1
        gg=g_1bar*(r_1bar/(r_1bar+z(j)))**2
        h=1.d+02*r*tl(j-1)/(ml(j-1)*gg)
        z(j-1)=z(j)+h*(ap(j)-ap(j-1))
        gl(j-1)=g_1bar*(r_1bar/(r_1bar+0.5*(z(j)+z(j-1))))**2
    19 end do
else
    i1=1
    h=1.d+02*r*0.5*(t_1bar+tl(1))/(ml(1)*g_1bar)
    z(1)=h*(a_1bar-ap(1))
endif

do 18 j=i1,nlay
    gg=g_1bar*(r_1bar/(r_1bar+z(j)))**2
    h=1.d+02*r*tl(j)/(ml(j)*gg)
    z(j+1)=z(j)+h*(ap(j)-ap(j+1))
    gl(j)=g_1bar*(r_1bar/(r_1bar+0.5*(z(j)+z(j+1))))**2
18 END DO

write(*,*)"--------------------------"
write(*,*)"Printing the parameters for the atmospheric levels and layers"
write(6,206)
206 format(/,5x,'Altitude',5x,'Pressure',3x,'Temperature',2x,'Mol. weight',4x,&
'q(ch4)',5x,'Gravity',/,8x,'(km)',8x,'(mbar)',8x,'(K)',10x,'(g)',18x,'(cm.sec-2)',/)
do j=1,nlay
    write(6,207)j,z(j),p(j),pl(j),tl(j),ml(j),qch4(j),gl(j)
end do
207 format(i4,f9.2,e14.4,/,15x,e14.4,f9.2,f13.3,e15.4,f10.2)
write(6,fmt="(I4,F9.2,E14.4)")nlev,z(nlev),p(nlev)
write(*,*)"Done."
write(*,*)"--------------------------"
write(*,*)"Input: Vertical profiles  of the molecular absorbers"
!**        idist=1 : vertical profile given in NIN points
!**        idist=2 : saturated distribution above condensation level
!**                  QUNSAT is the mole fraction below condensation level
!**        idist=3 : saturated distribution below condensation level
!**                  QUNSAT is the mole fraction above condensation level
allocate(ql(nlay,0:nb_mol))
do 3 ik=1,nb_mol
    print*,trim(NAME_MOL(ik,1))
    select case (NAME_MOL(ik,2))
        case('1')
            nbl = nb_line(NAME_MOL(ik,3))
            allocate(tmp_pl(nbl),tmp_aql(nbl),atmp_pl(nbl),atmp_aql(nbl))
            open(unit=10,file=NAME_MOL(ik,3),status='old')
            do j=1,nbl
                read(10,*) tmp_pl(j),tmp_aql(j)
                atmp_pl(j) = dlog(tmp_pl(j)) 
                atmp_aql(j) = dlog(tmp_aql(j)) 
            end do
            do 5 j=1,nlay
                ql(j,ik)=dexp(tina(apl(j),atmp_pl,atmp_aql,nbl))
            5 end do
            deallocate(tmp_pl,tmp_aql,atmp_pl,atmp_aql)
            close(10)
        case('2')
            read(NAME_MOL(ik,3),*)qunsat ! convert string NAME_MOL(ik,3) to float
            qsat=1.01325d+03*10**(a(ik)-b(ik)/tl(1))/pl(1)
            ql(1,ik)=dmin1(qsat,qunsat)
            do 6 j=2,nlay
                qsat=1.01325d+03*10**(a(ik)-b(ik)/tl(j))/pl(j)
                ql(j,ik)=dmin1(ql(j-1,ik),qsat)
            6 END DO
        case('3')
            read(NAME_MOL(ik,3),*)qunsat
            qsat=1.01325d+03*10**(a(ik)-b(ik)/tl(nlay))/pl(nlay)
            ql(nlay,ik)=dmin1(qsat,qunsat)
            do 7 j=nlay-1,1,-1
                qsat=1.01325d+03*10**(a(ik)-b(ik)/tl(j))/pl(j)
                ql(j,ik)=dmin1(ql(j+1,ik),qsat)
            7 END DO
        case default
            print *, "Type '1': input on file"
            print *, "Type '2': saturated distribution above condensation level"
            print *, "Type '3': saturated distribution below condensation level"
            stop
    end select
    if(name_mol(ik,1) == 'CH4') then
        write(*,*)"--------------------------"
        write(*,*)"Changing qch4 to ql(j,numch4) if CH4 is an absorber"
        do j=1,nlay
            qch4(j)=ql(j,ik)
        enddo
    end if
    write(6,211)ik,name_mol(ik,1),mass(ik)
    if(name_mol(ik,2) == '1') then
        write(6,*) trim(name_mol(ik,3))
    else
        if (name_mol(ik,2) == '2') then
            write(6,212)qunsat
        else
            write(6,213)qunsat
        endif
    endif
3 end do
211 format(/,10x,'Absorber #',i3,': ',a4,' Mol.weight =',f6.2,' g')
212 format(/,'Saturated profile for P<P_cond: unsaturated mixing ratio =',e11.4)
213 format(/,'Saturated profile for P>P_cond: unsaturated mixing ratio =',e11.4)
write(*,*)"Done."
write(*,*)"--------------------------"
write(*,*)"Printing the mole fractions in the atmospheric layers"
write(6,208)(name_mol(i,1),i=1,nb_mol)
do j=1,nlay
    write(6,209)j,pl(j),tl(j),(ql(j,i),i=1,nb_mol)
end do
208 format(10x,'P',10x,'T',20(6x,a4,1x))
209 format(i4,e12.4,f8.2,20(e11.3))
write(*,*)"Done."

!**             Input: Geometry of the observations
!**     nb_obs>=1 : vertical viewing (nb_obs points; input airmass pertains
!**                to the 1-bar level) (starting airmass and step in
!**                airmass are given). Transmittances are first calculated
!**                for an airmass of los.
!**     nb_obs=0  : radiances are calculated with the 2nd order integral
!**     nb_obs<0  : horizontal viewing (-nb_obs points)
!**                (starting altitude -km- and step in altitude are given)
!**                If los>0, the pressure grid is redefined. Transmittances
!**                are first calculated for altitude los.

!**       Calculation of the airmasses in the layers (nview>0)
!**       Calculation of the airmasses in the layers (nb_obs>0)
if(.not. limbe) then
    write(6,214),nb_obs,los(1),los(nb_obs)
    allocate(sec(nlay,nb_obs),lay_min(nb_obs))
    lay_min = 0
    do 23 i=1,nb_obs
        air=los(i)
        teta1=dacos(1./air)
        do 24 j=1,nlay
            r2=r_1bar+z(j)
            r3=r_1bar+z(j+1)
            dteta=dasin(dsin(teta1)*r_1bar/r3)-dasin(dsin(teta1)*r_1bar/r2)
            sec(j,i)=dsqrt((r3*r3+r2*r2-2.d+00*r2*r3*dcos(dteta)))/(r3-r2)
        24 end do
    23 end do
else
    if(nb_obs == 0) then
        write(6,215)
        nb_obs=1
        allocate(sec(nlay,1))
        do 25 j=1,nlay
            sec(j,1)=1.
        25 end do
    else
        write(6,216)nb_obs,los(1),los(nb_obs)
        if(los(1) <= z(1)) goto 33 ! == nadir
    !**       Redefinition of the pressure grid (nb_obs<0)
        allocate(tmp_p(nlev),tmp_ap(nlev),tmp_t(nlev),tmp_z(nlev))
        allocate(tmp_pl(nlay),tmp_tl(nlay),tmp_apl(nlay),tmp_qch4(nlay))
        allocate(tmp_ml(nlay),tmp_gl(nlay),tmp_ql(nlay,0:nb_mol))
        allocate(tmp_aql(nlay))
        tmp_p = p
        tmp_ap = ap
        tmp_t = t
        tmp_z = z
        tmp_pl = pl
        tmp_tl = tl
        tmp_apl = apl
        tmp_qch4 = qch4
        tmp_ml = ml
        tmp_gl = gl
        tmp_ql = ql
        deallocate(p,ap,t,z,pl,tl,apl,qch4,ml,ql,gl)
        do 26 j=2,nlev
            if(los(1) >= tmp_z(j)) goto 26
            ja=j-1
            goto 27
        26 END DO
        print 218,los(1)
        218 format('Viewing altitude',f8.2,' is above the atmospheric grid')
        stop
    27  tmp_nlev = nlev
        tmp_nlay = nlay
        nlev=nlev-ja+1
        nlay=nlev-1
        allocate(p(nlev),ap(nlev),z(nlev),t(nlev))
        allocate(tl(nlay),pl(nlay),apl(nlay),qch4(nlay),ml(nlay))
        allocate(gl(nlay),ql(nlay,0:nb_mol),aql(tmp_nlay))
        p(1)=dexp(tina(los(1),tmp_z,tmp_ap,tmp_nlev))
        t(1)=dexp(tina(dlog(p(1)),api,ati,ninput))
        z(1)=los(1)
        write(6,219)los(1),p(1)
        219 format('Atmospheric grid is redefined from altitude =',f8.2,' km',/,36x,'Pressure level =',e11.4,' mbar')
        ap(1)=dlog(p(1))
        pl(1)=dsqrt(p(1)*tmp_p(ja+1))
        apl1=dlog(pl(1))
        tl(1)=dexp(tina(apl1,api,ati,ninput))
        ml(1)=tina(apl1,tmp_apl,tmp_ml,tmp_nlay)
        gl(1)=tina(apl1,tmp_apl,tmp_gl,tmp_nlay)
        do j=1,tmp_nlay
            aql(j)=dlog(tmp_qch4(j))
        end do
        qch4(1)=dexp(tina(apl1,tmp_apl,aql,tmp_nlay))
        apl(1)=apl1
        tmp_apl(1)=apl1
        do 28 ik=1,nb_mol
            do 29 j=1,tmp_nlay
                aql(j)=dlog(tmp_ql(j,ik))
            29 END DO
            ql(1,ik)=dexp(tina(apl1,tmp_apl,aql,tmp_nlay))
        28 END DO
        do 31 j=2,nlev
            z(j)=tmp_z(j+ja-1)
            p(j)=tmp_p(j+ja-1)
            t(j)=tmp_t(j+ja-1)
            ap(j)=tmp_ap(j+ja-1)
        31 end do
        do 310 j=2,nlay
            pl(j)=tmp_pl(j+ja-1)
            apl(j)=tmp_apl(j+ja-1)
            tl(j)=tmp_tl(j+ja-1)
            ml(j)=tmp_ml(j+ja-1)
            gl(j)=tmp_gl(j+ja-1)
            qch4(j)=tmp_qch4(j+ja-1)
            do 32 ik=1,nb_mol
                ql(j,ik)=tmp_ql(j+ja-1,ik)
            32 end do
        310 end do
        ql(:,0) = 0.
        deallocate(tmp_p,tmp_ap,tmp_t,tmp_z,tmp_pl,tmp_tl,tmp_apl,tmp_qch4)
        deallocate(tmp_ml,tmp_gl,tmp_ql,tmp_aql)

    !**       Calculation of the airmasses in the layers (nb_obs<0)
    !**  lay_min(i) is the layer where altitude AIR is reached
    33  allocate(lay_min(nb_obs),sec(nlay,nb_obs))
        sec = 0.
        do 34 i=1,nb_obs
            lay_min(i)=nlev
            air=los(i)
            r1=r_1bar+air
            if(air < z(1)) then
                lay_min(i)=0
                write(6,222)air
                222 format('Altitude',f10.2,' km reached below first layer')
            endif
            do 35 j=1,nlay
                r2=r_1bar+z(j)
                r3=r_1bar+z(j+1)
                if(air < z(j)) then
                    sec(j,i)=(dsqrt(r3*r3-r1*r1)-dsqrt(r2*r2-r1*r1))/(r3-r2)
                else
                    if(air < z(j+1)) then
                        sec(j,i)=dsqrt(r3*r3-r1*r1)/(r3-r2)
                        lay_min(i)=j
                        write(6,220)air,lay_min(i)
                        220 format('Altitude',f10.2,' km reached in layer #',i3)
                    end if
                end if
            35 end do
            if(lay_min(i) > nlay) then
                nb_obs=i-1
                print 221, air,nb_obs
                221 format(' WARNING: Altitude',f10.2,' km is above the atmospheric grid;'/, &
                & '   only',i3,' altitudes are hereafter considered')
                goto 36
            endif
        34 end do
        if(lay_min(nb_obs) == 0) then
            nb_obs=-nb_obs
        end if
    end if
end if

214 format(/,'Spectral calculations for',i3,' airmasses',/,' from',f8.4,' to',f8.4)
215 format(/,'Spectral calculations of the disk-averaged radiance (plan parallel geometry)')
216 format(/,'Spectral calculations for',i3,' altitudes',/,' from',f10.2,' to', f10.2,' km')
!**
36 if(nb_obs == 0) then
    do 50 i=0,500
        x=0.02*dfloat(i)
        et(i)=dexp(-x)-x*eint(-x)
    50 end do
else
    do 51 i=0,500
        et(i)=dexp(-0.02*dfloat(i))
    51 end do
end if

allocate(at(nlev))
do 52 j=1,nlev
    at(j)=dlog(t(j))
52 end do

!**             Input: Spectral range (cm-1)
nfreq_inv = nint((F_end-F_start)/F_pas)+1

write(6,1207)F_start,F_end,F_pas_mult,F_pas
1207 format(/,'Spectrum if computed from',F9.3,' to',F9.3,' cm-1',/,'Calculation step is ',F6.1,&
&' times the smallest line halfwidth',/,'Output step is ',F9.5,' cm-1')

!**		Input: Output file, convolution parameters, weighting functions
write(6,1208)file_spec_synth,convfunc,rayleigh
1208 format('Output file is: ',A80,/,'Convolution type is "',A1,'"',&
'Resolution is',F11.5,' cm-1')

!**             Input: Parameters for cloud opacity (nb_cloud clouds)
!*                      TAUCLT: Total optical depth (nadir) at wavenumber WREF
!*                              (cm-1)
!*                      PBOT, PTOP: Boundaries of the cloud (mbar)
!*                      HPHG: Particle-to-gas scale height
allocate(taucl_tot(nlay),taucl(nb_cloud,nlay))
if(nb_cloud > 0) then
    do j=1, nlay
        taucl_tot(j) = 0.
    enddo
    do ic=1,nb_cloud
        do 1009 j=1,nlay
            taucl(ic,j)=0.
        1009 end do
        ninput_cl = nb_line(PROF_HAZE(ic))
        allocate(pcli(ic,ninput_cl),taucli(ic,ninput_cl),apcli(ic,ninput_cl),ataucli(ic,ninput_cl))
        open(unit=60,file=PROF_HAZE(ic),status='old')
        do j=1,ninput_cl
            read(60,*) pcli(ic,j),taucli(ic,j)
            apcli(ic,j)=dlog(pcli(ic,j))
            ataucli(ic,j)=dlog(taucli(ic,j))
        enddo
        close(60)
        print*,'ninput_cl',ninput_cl
        do j=1,ninput_cl
            print *,'j',j,pcli(ic,j),'-',taucli(ic,j)
        end do
        do j=1, nlay
            taucl(ic,j) = dexp(tina(apl(j),apcli,ataucli,ninput_cl))
            print *,'j',j,p(j),'-',p(j+1),taucl(ic,j)
        enddo
        do j=1,nlay
            taucl_tot(j) = taucl_tot(j) + taucl(ic,j)
        enddo
    enddo
endif
!**                     Input: Spectroscopic parameters
!*        NLOR: Number of different Lorentz halfwidths for gas K
!*        ALOR: Lorentz profile is computed up to ALOR cm-1 from line center
!*        GLOR: Lorentz halfwidths (1 --> NLOR)
!*        ELOR: T-exponent (1 --> NLOR)
!*        G2LOR: GLOR is used for lines having G2LOR(i-1)< G < G2LOR(i) (296 K)
!*              Selection is made in subroutine READ_LINE
!*        NVIB: Number of vibrational modes to be considered for calculation
!*              of the partition function
!*        VIB, NDEG: Frequency and degeneracy of the vibrational modes
!*              (1 --> NVIB)
!**

do 190 k=1,nb_mol
    if(nlor(k) > 5) then
        print 2007, k,name_mol(k,1),nlor(k)
        2007 format(' Molecule #',i3,' (',a4,') must not have more than 5 line halfwidths (',i2,' asked)')
        stop 16
    end if
    if(nvib(k) > 5) then
        print 2011, k,name_mol(k,1)
        2011 format(' Molecule #',i3,' (',a4,') must not have more than 5 vibration modes')
        stop 16
    end if
190 end do

!**             Input: Absorbers, frequencies (cm-1) and airmasses
!**		       for which K matrices are calculated
! DF(1,X) = start, DF(2,X) = end, DF(3,X) = pas, NB_DF = nb de plage
! DF : ligne 1 = start, ligne 2 = end, ligne 3 = pas
allocate(n_x(nb_df))
do i=1,nb_df
    n_x(i) = nint((DF(i,2)-DF(i,1))/DF(i,3))+1
    nx = nx + n_x(i)
end do
print*,'nx =',nx

allocate(freq_x(nx))
cpt=1
do j = 1, nb_df
    do i=1,n_x(j)
        freq_x(cpt)=DF(j,1)+DF(j,3)*(i-1)
        cpt=cpt+1
    end do
end do
write(*,*)"Done."

write(*,*)"--------------------------"
write(*,*)"Open files observed"
print 1120,(file_obs(is),is=1,nb_obs)
1120 format(/,'Files of observed spectra:',(/,3x,A256))

allocate(rad_obs(nb_obs,nx))
do 810 is=1,nb_obs
    open(unit=10,file=file_obs(is),status='old')
    read(10,*)
    read(10,*) np
    read(10,*)
    allocate(wobs(np),robs(np))
    do 811 i=1,np
        read(10,*) wobs(i),robs(i)
    811 end do
    do 812 i=1,nx
        rad_obs(is,i)=1.e+07*tina(freq_x(i),wobs,robs,np)
    812 end do
    deallocate(wobs,robs)
810 end do

write(*,*)"--------------------------"
write(*,*)"Read file noise"
write(*,fmt='(TR5,A,A)')'File of NESRs: ',trim(path_input)//trim(file_noise)
write(*,fmt='(TR5,A,I3)')'Number of spectrum: ',nsp
nbl_noise = nb_line(trim(path_input)//file_noise) - 1
open(unit=10,file=trim(path_input)//file_noise,status='old')
allocate(wnoise(nbl_noise),rnoise(nbl_noise))
read(10,*) !premiere ligne est le titre
do i=1,nbl_noise
    read(10,*)wnoise(i),rnoise(i)
end do
close(10)
allocate(rad_noise(nx))
do i=1,nx
    rad_noise(i)=1.e+07*tina(freq_x(i),wnoise,rnoise,nbl_noise)/dsqrt(dfloat(nsp))
END DO
write(*,fmt='(A,TR2,F8.3,A,TR2,E12.3)')'Freq_x(1): ',freq_x(1),', rad_noise(1): ',rad_noise(1)
write(*,fmt='(A,TR1,F8.3,A,TR1,E12.3)')'Freq_x(nx): ',freq_x(nx),', rad_noise(nx): ',rad_noise(nx)
write(*,*)"Done."


print 1122,(name_mol_inv(i),corr_len(i),wei(i),i=1,nb_mol_inv)
1122 format(/,(1x,a4,': Correlations length:', f7.2,' scale height',5x,'weight =',e10.3))
print 1124, NITER,FILE_PROFQ
1124 format(/,(i2,' iterations, abundance profile written on: ',a80))
print 1126, niter,file_prof_t
1126 format(/,i2,' iterations, T profile written on: ',A80)

if(nb_kernel > 0) then
    print 1123, FILE_KERNEL,NB_KERNEL,(F_KERNEL(i),i=1,NB_KERNEL)
    1123 format(/,'Kernels are written on file: ',A80,/,I6,'frequencies: ',8F10.3)
endif
if(nb_kernel_t > 0) then
    print 1125, file_kernel_t,nb_kernel_t,(f_kernel_t(i),i=1,nb_kernel_t)
    1125 format(/,'T Kernels are written on file: ',A80,/,I6,' frequencies: ',8F10.3)
endif

!**             Spectral range, output file, convolution
!**             Spectroscopic file H2-H2, H2-He, and H2-CH4
!*	  N2N2(j,n), N2CH4(j,n), CH4CH4(j,n), N2H2(j,n) are the N2-N2,
!*        N2-CH4, CH4-CH4 and N2-H2 absorption coefficients at
!*        temperature TL(j) for the nth frequency
!*	  (f1_cont is the first frequency, df_cont is the step)
!**		Read N2-H2-CH4 spectroscopic file (unit 12)
!**		    LT temperatures and KNU frequencies
write(*,*)"--------------------------"
write(*,*)"Collisions"
!~ N2N2 collisions
open(unit=12,file=trim(path_input)//FILE_N2N2,status='old',form='unformatted',access='stream')
read(12,pos=5)lt,knu
close(12)
allocate(n2n2_g(lt,knu))
allocate(tmod(lt),gnu(knu))
open(unit=12,file=trim(path_input)//FILE_N2N2,status='old',form='unformatted')
read(12)lt,knu,(tmod(i),i=1,lt),(gnu(i),i=1,knu),((n2n2_g(i,j),i=1,lt),j=1,knu)
close(12)
write(6,800)FILE_N2N2,gnu(1),gnu(knu),lt,knu
800 format('N2-N2 spectroscopic file :',A20,/,'from',F8.2,'to',F8.2,'cm-1: ',I3,'temperatures',I4,'frequencies')
!~ N2CH4 collisions
allocate(n2ch4_g(lt,knu))
open(unit=12,file=trim(path_input)//FILE_N2CH4,status='old',form='unformatted')
read(12)lt,knu,(tmod(i),i=1,lt),(gnu(i),i=1,knu),((n2ch4_g(i,j),i=1,lt),j=1,knu)
close(12)
write(6,801)FILE_N2CH4,gnu(1),gnu(knu),lt,knu
801 format('N2-CH4 spectroscopic file :',A20,/,'from',F8.2,'to',F8.2,'cm-1: ',I3,'temperatures',I4,'frequencies')

!~ CH4CH4 collisions
allocate(ch4ch4_g(lt,knu))
open(unit=12,file=trim(path_input)//FILE_CH4CH4,status='old',form='unformatted')
read(12)lt,knu,(tmod(i),i=1,lt),(gnu(i),i=1,knu),((ch4ch4_g(i,j),i=1,lt),j=1,knu)
close(12)
write(6,802)FILE_CH4CH4,gnu(1),gnu(knu),lt,knu
802 format('CH4-CH4 spectroscopic file :',A20,/,'from',F8.2,'to',F8.2,'cm-1: ',I3,'temperatures',I4,'frequencies')

!~ N2H2 collisions
allocate(n2h2_g(lt,knu))
open(unit=12,file=trim(path_input)//FILE_N2H2,status='old',form='unformatted')
read(12)lt,knu,(tmod(i),i=1,lt),(gnu(i),i=1,knu),((n2h2_g(i,j),i=1,lt),j=1,knu)
close(12)
write(6,803)FILE_N2H2,gnu(1),gnu(knu),lt,knu
803 format('N2-H2 spectroscopic file :',A20,/,'from',F8.2,'to',F8.2,'cm-1: ',I3,'temperatures',I4,'frequencies')

df_cont=gnu(2)-gnu(1)
n1=max0(idint((F_start-gnu(1))/df_cont)+1,1)      ! F_start dans le spectre qu'on calcule
n2=min0(idint((F_end-gnu(1))/df_cont)+2,knu)    ! pareil pour F_end
if(n1 >= knu) then
    n1=knu-1
end if
if(n2 <= 1) then
    n2=2
end if
nfcont=n2-n1+1
f1_cont = gnu(n1)
write(6,806)f1_cont,gnu(n2)
806 format('Data from',F8.2,'to',F8.2,' cm-1 are used')
write(*,*)"Done."

allocate(tab_n2n2(nlay,nfcont))
allocate(tab_n2ch4(nlay,nfcont))
allocate(tab_ch4ch4(nlay,nfcont))
allocate(tab_n2h2(nlay,nfcont))

call dat_cont(n2n2_g,n2ch4_g,ch4ch4_g,n2h2_g,tmod,n1,n2,tl,nlay,tab_n2n2,tab_n2ch4,tab_ch4ch4,tab_n2h2,knu,nfcont,lt)

deallocate(tmod,gnu)
deallocate(n2n2_g,n2ch4_g,ch4ch4_g,n2h2_g)
!**             Input: File for calculation of the Voigt profile
open(unit=10,file=trim(path_input)//'VOIGT.DAT',status='old',form='unformatted')
read(10)v
close(10)
!**             Spectral extinction of cloud IC particles
allocate(qref(nb_cloud))
allocate(nond(nb_cloud))
if(nb_cloud > 0) then
    do ic=1,nb_cloud
        open(unit=9,file=trim(path_input)//file_scat(ic),status='old',form='formatted')
        read(9,*) nond(ic)
        close(9)
    end do
    nond_max = maxval(nond)
    allocate(sig(nb_cloud,nond_max),qex(nb_cloud,nond_max))
    do ic = 1,nb_cloud
        open(unit=9,file=trim(path_input)//file_scat(ic),status='old',form='formatted')
        read(9,*)
        do i=1,nond(ic)
            read(9,*) sig(ic,i),qex(ic,i)
        enddo
        close(9)
        qref(ic)=tina(f_cloud(ic),sig(ic,:),qex(ic,:),nond(ic))
    enddo
end if

!**		Absorbers, frequencies, airmasses for which K matrices are calculated
ikh=0
allocate(icorps_k(nb_mol_inv))
do 12 ik=1,nb_mol_inv
    if(name_mol_inv(ik) == 'HAZE') then
        icorps_k(ik)=0
        ikh=ik
        goto 12
    end if
    do 120 i=1,nb_mol
        if(name_mol_inv(ik) /= name_mol(i,1)) goto 120
        icorps_k(ik)=i
        goto 12
    120 end do
12 end do

allocate(nbl_sp(nb_mol))

do k = 1, nb_mol
    nbl_sp(k) = nb_line(trim(path_input)//file_spectro(k))
        if (nbl_sp(k) >= nbl_sp_max) then
        nbl_sp_max = nbl_sp(k)
    end if
end do
allocate(w_input(nb_mol,nbl_sp_max),s_input(nb_mol,nbl_sp_max),g_input(nb_mol,nbl_sp_max),e_input(nb_mol,nbl_sp_max))

do k = 1, nb_mol
    open(unit=24,file=trim(path_input)//file_spectro(k),status='old')
    do j = 1, nbl_sp(k)
        read(24,fmt="(f10.6,e10.3,f5.3,f10.3)")w_input(k,j),s_input(k,j),g_input(k,j),e_input(k,j)
    end do
    close(24)
end do

write(*,*)"--------------------------"
write(*,*)"Check inversion parameters"
nb_obs = abs(nb_obs)
write(*,fmt='(TR5,A,I2,A,20(A,X))')&
    'K matrices are calculated for ',nb_mol_inv,' absorbers: ',(trim(name_mol_inv(i)),i=1,nb_mol_inv)
write(*,fmt='(TR5,A,F8.3,A,F8.3,A,I4,A)')&
    'Inverse calculation between ',freq_x(1),' and ',freq_x(nx), ' cm-1 (',nx,' points)'
write(*,fmt='(TR5,A,I3)')'Numbers of airmasses: ',nb_obs
!~ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**			Iterations
!~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
m=nb_obs*nx
allocate(matrix_k(nb_mol_inv,nlay,nb_obs,nfreq_inv),matrix_t(nlay,nb_obs,nfreq_inv),rad(nb_obs,nfreq_inv),rad_diff(nb_obs,nx))
allocate(matrix_r(nlay,m))
allocate(s(nlay,nlay,nb_mol_inv),sk(nb_mol_inv,nlay,m),ssk(nb_mol_inv,m,m),sm(nlay,m),s_b(nlay,nlay),ssk_b(m,m))
allocate(matrix_r2(nlay,m),errt(nlay),errtb(nlay),cf(nlay,nlay),cf_b(nlay,nlay))
allocate(dlnerrqb(nlay,nb_mol),errqb_sup(nlay,nb_mol),errqb_inf(nlay,nb_mol),dlnerrq(nlay,nb_mol),errq_sup(nlay,nb_mol))
allocate(errq_inf(nlay,nb_mol))
rad_diff = 0d0
do 500 iter=0,niter
    print 307, iter
    307 format(/,30x,' ******* ITER =',i2,' *******')
     call transfer_transac(nlev,nlay,vmr_tropo_n2,vmr_tropo_h2,vmr_tropo_ar,pl,tl,ml,gl,qch4,ql,p,ap,t,nb_obs,&
    sec,lay_min,nb_mol,mass,expo_rot,et,f_start,f_end,f_pas_mult,f_pas,file_weightfunc,&
    file_spec_synth,convfunc,rayleigh,itrans,w_input,s_input,g_input,e_input,nb_cloud,taucl,cross_section_min,&
    temperature_test,nlor,alor,glor,elor,g2lor,nvib,vib,ndeg,nb_mol_inv,icorps_k,matrix_k,matrix_t,rad,ikh,nfcont,&
    nond,nond_max,sig,qex,nbl_sp_max,nbl_sp,f1_cont,df_cont,iter,nfreq_inv,qref,tab_n2n2,tab_n2ch4,tab_ch4ch4,tab_n2h2,v,&
    name_mol(:,1),path_input,limbe)

    print *, 'test transac 1'
    print 355
    355 format(/,3x,'freq_x',5x,'noise',10x,'obs',8x,'synth',3x,'obs-synth',/)
    do 311 ix=1,nx
        print*,'ix =',ix,'freq_x(ix) =',freq_x(ix),'rad_noise(ix) =',rad_noise(ix)
        i=nint((freq_x(ix)-F_start)/F_pas+0.5)
        do 317 is=1,nb_obs
            rad_diff(is,ix)=rad_obs(is,ix)-rad(is,i)
        317 END DO
        print 354,freq_x(ix),rad_noise(ix),rad_obs(1,ix),rad(1,i),rad_diff(1,ix)
        do is =2,nb_obs
            print 357, rad_obs(is,ix),rad(is,i),rad_diff(is,ix),is
        end do
    311 END DO
    354 format(F10.3,E12.3,2X,3E12.4,I5)
    357 format(24X,3E12.4,I5)
! -------------------------------------------------------------------
    print *, 'NESR**2 * nb_wavenumber =',rad_noise(nx)**2 * nx
    print *, '*** rms ***'

    do is=1,nb_obs
        rms=0.
        do 314 ix=1,nx
            rms=rms+rad_diff(is,ix)**2
        314 END DO
        print 356, is, dsqrt(rms/dfloat(nx))
        356 format(/,' is =',i2,': rms =',e10.3)
    END DO

    print *, 'test transac 3'

!**		Calculate S * Kt matrices
    do 401 ik=1,nb_mol_inv
        corr2=corr_len(ik)**2
        do i=1,nlay
            do j=1,m
                sk(ik,i,j)=0.
            enddo
        enddo
        do i=1,nlay
            do k=1,nlay
                plo=dlog(pl(i)/pl(k))
                var_expo=-0.5d+00*plo*plo/corr2
                if(var_expo < -200) then
                    corr= 0d0
                else
                    corr=wei(ik)*exp(var_expo)
                end if
                s(i,k,ik)=corr
                do j=1,m
                    is=1+(j-1)/nx
                    ix=j-(is-1)*nx
                    ifreq=nint((freq_x(ix)-F_start)/F_pas+0.5)
                    sk(ik,i,j)=sk(ik,i,j)+corr*matrix_k(ik,k,is,ifreq)
                enddo
            enddo
        enddo
    401 END DO

    print *, 'test transac 4'

!**			Calculate K * S *Kt matrices
    do 400 ik=1,nb_mol_inv
        do i=1,m
            is=1+(i-1)/nx
            ix=i-(is-1)*nx
            ifreq=nint((freq_x(ix)-F_start)/F_pas+0.5)
            do j=1,m
                ssk(ik,i,j)=0.
                do k=1,nlay
                    ssk(ik,i,j)=ssk(ik,i,j)+matrix_k(ik,k,is,ifreq)*sk(ik,k,j)
                enddo
            enddo
        enddo
        tra=0.
        do j=1,m
            tra=tra+ssk(ik,j,j)
        enddo
        print 899, ik,tra
    400 END DO
    899 format(/,'Trace of matrix a K * S * Kt  for absorber',I2,':',F10.5)

    print *, 'test transac 5'

!***    Calculate S*Mt matrice (avec Mt = transposee kernel pour inv T)
    corrT2=corr_len_T**2
    do i=1,nlay
        do j=1,m
            sm(i,j)=0.
        enddo
    enddo
    do i=1,nlay
        do k=1,nlay
            plo=dlog(pl(i)/pl(k))
            var_expo=-0.5d+00*plo*plo/corrT2
            if (var_expo < -200) then
                corr_b = 0d0
            else
                corr_b = wei_t*exp(var_expo)
            end if
            s_b(i,k)=corr_b
            do j=1,m
                is=1+(j-1)/nx
                ix=j-(is-1)*nx
                ifreq=nint((freq_x(ix)-F_start)/F_pas+0.5)
                sm(i,j)=sm(i,j)+corr_b*matrix_t(k,is,ifreq)
            enddo
        enddo
    enddo
    print *, 'test transac 6'


!**			Calculate M * S * Mt matrices
    do i=1,m
        is=1+(i-1)/nx
        ix=i-(is-1)*nx
        ifreq=nint((freq_x(ix)-F_start)/F_pas+0.5)
        do j=1,m
            ssk_b(i,j)=0.
            do k=1,nlay
                ssk_b(i,j)=ssk_b(i,j)+matrix_t(k,is,ifreq)*sm(k,j)
            enddo
        enddo
    enddo
    tra_b=0.
    do j=1,m
        tra_b=tra_b+ssk_b(j,j)
    enddo
    print 890, tra_b
    890 format(/,'Trace of matrix a M * S * Mt (for temperature):',f10.5)
    
    print *, 'test transac 7'
    
!~  Calculate Matrix C
    allocate(c(m,m))
    allocate(c0(m))
    
    do 403 i=1,nx
        do is=1,nb_obs
            j=i+(is-1)*nx
            c(j,j)=rad_noise(i)**2
            c0(j)=c(j,j)
        enddo
    403 END DO

    do 402 i=1,m
        do 402 j=1,m
            if(i /= j) c(i,j)=0.
            do ik=1,nb_mol_inv
                c(i,j)=c(i,j)+ssk(ik,i,j)+ssk_b(i,j)
            enddo
    402 END DO
     

    print *, 'test transac 8'

!**			Invert Matrix C

    if (method_invert) then
        call choldc(c,m)
    else
        call matrix_inv(c,m)
    end if

    print *, 'test transac 9'

!**	       Calculate final matrice and variations of T
    do i=1,nlay
        do j=1,m
            matrix_r2(i,j)=0.
            do k=1,m
                matrix_r2(i,j)=matrix_r2(i,j)+sm(i,k)*c(k,j)
            enddo
        enddo
    enddo
          
    print *, 'test transac 10'

    do i=nlay,1,-1
        dt=0.
        do 407 j=1,m
            is=1+(j-1)/nx
            ifreq=j-(is-1)*nx
            dt=dt+matrix_r2(i,j)*rad_diff(is,ifreq)
        407 END DO
        tl(i)=tl(i)+dt
        print 1111, i,pl(i),dt,tl(i)
    end do
    t(1)=t(1)+dt
    1111 format(I3,E14.4,F10.3,5X,F10.3)
!**		Calculate covariance matrix
    do 421 i=1,nlay
        errtb(i)=0.
        do 423 j=1,m
            errtb(i)=errtb(i)+matrix_r2(i,j)**2*c0(j)
        423 END DO
        do 421 k=1,nlay
            cf_b(i,k)=s_b(i,k)
            do  j=1,m
                cf_b(i,k)=cf_b(i,k)-matrix_r2(i,j)*sm(k,j)
            enddo
    421 END DO
       
          
    print *, 'test transac 11'
       

    do i=1,nlay
        errt(i)=dsqrt(cf_b(i,i))
        errtb(i)=dsqrt(errtb(i))
    enddo

    print *, 'test transac 12'

!          Writing T vertical profile
    open(unit=10,file=file_prof_t,status='unknown')
    do i=1,ninput
        if(ppi(i) > p(1)) THEN
            write(10,1113) ppi(i),ti(i)
        ELSE
            goto 510
        ENDIF
    enddo
    510 write(10,1113) p(1),t(1)
    write(10,1113) (pl(i),tl(i),errt(i),errtb(i),i=1,nlay)
    print *,('errt =',errt(i),i=1,nlay)
    close(10)
    1113 format(e12.4,f10.3,2f12.2)
    print *, 'test transac 13'
     
! -------------------------------------
!**	       Calculate final matrices and variations of absorbers
    do 410 ik=1,nb_mol_inv
        print *, name_mol_inv(ik)
        do 404 i=1,nlay
            do 405 j=1,m
                matrix_r(i,j)=0.
                do k=1,m
                    matrix_r(i,j)=matrix_r(i,j)+sk(ik,i,k)*c(k,j)
                enddo
            405 END DO
        404 END DO
        do 406 i=nlay,1,-1
            dlnq=0.
            do 413 j=1,m
                is=1+(j-1)/nx
                ifreq=j-(is-1)*nx
                dlnq=dlnq+matrix_r(i,j)*rad_diff(is,ifreq)
            413 END DO
            if(icorps_k(ik) > 0) then
                ql(i,icorps_k(ik))=ql(i,icorps_k(ik))*dexp(dlnq)
            !** Remove any supersaturation
                qsat=1.01325d+03*10**(a(ik)-b(ik))/tl(i)/pl(i)
                ql(i,icorps_k(ik))=dmin1(qsat,ql(i,icorps_k(ik)))
                print 1100,i,pl(i),dlnq,ql(i,icorps_k(ik))
            else
                taucl_tot(i)=taucl_tot(i)*dexp(dlnq)
                print 1100,i,pl(i),dlnq,taucl_tot(i)
                if(nb_cloud > 0) then
                    do ic=1,nb_cloud
                        taucl(ic,i)=taucl(ic,i)*dexp(dlnq)
                    enddo
                endif
            endif
        406 END DO
        1100 format(I3,E14.4,F12.4,E15.4)
         
        print *, 'test transac 14'

        if(name_mol_inv(ik) == 'HAZE') then
            do j=1,nlay
                ql(j,icorps_k(ik))=taucl_tot(j)
            enddo
        endif

        if(iter == niter) then
        !       Calcul des barres d'erreur avec la premiere methode
        !       (methode de Barney)
            do i=1,nlay
                dlnerrqb(i,ik)=0.
                do j=1,m
                    dlnerrqb(i,ik)=dlnerrqb(i,ik)+matrix_r(i,j)**2*c0(j)
                enddo
            enddo
            do i=1,nlay
                dlnerrqb(i,ik)=dsqrt(dlnerrqb(i,ik))
                errqb_sup(i,ik)=ql(i,icorps_k(ik))*dexp(dlnerrqb(i,ik))
                errqb_inf(i,ik)=ql(i,icorps_k(ik))/dexp(dlnerrqb(i,ik))
            enddo
        
            print *, 'test transac 15'


        !       Calcul des barres d'erreur avec la seconde methode
            do i=1,nlay
                do k=1,nlay
                    cf(i,k)=s(i,k,ik)
                    do j=1,m
                        cf(i,k)=cf(i,k)-matrix_r(i,j)*sk(ik,k,j)
                    enddo
                enddo
            enddo
            do i=1,nlay
                dlnerrq(i,ik)=dsqrt(cf(i,i))
                errq_sup(i,ik)=ql(i,icorps_k(ik))*dexp(dlnerrq(i,ik))
                errq_inf(i,ik)=ql(i,icorps_k(ik))/dexp(dlnerrq(i,ik))
            enddo
        !       fin calcul barres d'erreur
        endif
    410 end do

    deallocate(c0)
    deallocate(c)
500 end do
  
print *, 'test transac 16'

!-----------------------------------------------------------
!       Writing abundance vertical profiles
open(unit=10,file=file_profq,status='unknown')
write(10,1112) (name_mol_inv(ik),ik=1,nb_mol_inv), 'profT'
1112 format(28x,3(4x,a4,4x,48x), a5)
do 411 j=nlay,1,-1
    write(10,1101) 0.5*(z(j)+z(j+1)),pl(j),tl(j), &
    (ql(j,icorps_k(ik)),errq_inf(j,ik),errq_sup(j,ik), &
    errqb_inf(j,ik),errqb_sup(j,ik),ik=1,nb_mol_inv), tl(j),errt(j),errtb(j)
411 end do
close(10)

1101    format(f8.2,e11.4,f9.2,15e12.4)

print *, 'test transac 18'
!**             Calculate kernels
open(unit=11,file=file_kernel,status='unknown')
do 408 ifk=1,nb_kernel
    i=nint((f_kernel(ifk)-F_start)/F_pas+0.5)
    write(11,1103) F_start+F_pas*(dfloat(i)-0.5d+00)
    do 412 ik=1,nb_mol_inv
        write(11,1112) name_mol_inv(ik)
        do j=nlay,1,-1
            write(11,1104) 0.5*(z(j)+z(j+1)),pl(j),tl(j),(matrix_k(ik,j,is,i)/(ap(j)-ap(j+1))/rad(is,i),is=1,nb_obs)
        enddo
    412 end do
408 END DO
close(11)
1103 format(/,20x,f10.3,' cm-1')
!**

print *, 'test transac 19'
!**		Calculate kernels for T
open(unit=12,file=file_kernel_t,status='unknown')
do ik=1,nb_kernel_t
    i=nint((f_kernel_t(ik)-F_start)/F_pas+0.5)
    write(12,1102) F_start+F_pas*(dfloat(i)-0.5d+00)
    do j=nlay,1,-1
        write(12,1104) 0.5*(z(j)+z(j+1)),pl(j),tl(j),(matrix_t(j,is,i)/(ap(j)-ap(j+1))/rad(is,i),is=1,nb_obs)
    enddo
enddo
close(12)
1102    format(/,20x,f10.3,' cm-1')
1104 format(f8.2,e13.5,f8.2,2x,*(e13.5))


!~ =====================================================================
deallocate(z,p,ap,t,pl,apl,tl,gl,qch4,tab_i,lay_min,sec,ql,ppi,ti,ati,&
pcli,taucli,apcli,ataucli,taucl,at,taucl_tot,nlor,alor,glor,elor,g2lor,nvib,vib,ndeg,cross_section_min,temperature_test,&
freq_x,rad_obs,rad_diff,rad_noise,sk,sm,ssk,s,s_b,ssk_b,cf,cf_b,dlnerrqb,errqb_sup,errqb_inf,dlnerrq,errq_sup,&
errq_inf,errt,errtb,ml,matrix_k,matrix_t,matrix_r,matrix_r2,rad,tab_n2n2,tab_n2ch4,tab_ch4ch4,tab_n2h2,qref,n_x)

!**
call cpu_time(time_end)
time_elapsed = time_end - time_begin
print *,'=================================================='
print *,'Program takes :',time_elapsed,' s'


end program TRANSAC
