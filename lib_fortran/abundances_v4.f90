!***             Abundance profile inversion      Titan (TIT)
!*  nb_mol = Number of molecular absorbers
!*  NLEV   = Number of atmospheric levels
!*     Z(NLEV) : Altitude levels (numbered from 1 - lower boundary - to NLEV - Upper boundary -). Unit: km.
!*               (z=0 at the 1-bar pressure level)
!*     P(NLEV) : Pressure levels (numbered from 1 - lower boundary - to NLEV - Upper boundary -). Unit: mbar.
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
!*     PI(NINPUT) : Pressure of the input points (mbar)
!*     TI(NINPUT) : Temperature of the input points (K)
!*     API(NINPUT): Log(e) of pressure of the input points
!*     ATI(NINPUT): Log(e) of temperature of the input points
!*  NIN = Number of input points for concentration profiles
!*     PQI(NIN) : Pressure at the input points
!*     APQI(NIN): Log(e) of pressure of the input points
!*     QI(NIN)  : Mole fraction at the input points
!*     AQI(NIN) : Log(e) of mole fraction at the input points
!======================================================================!
!===                      TABLE OF CONTENTS                         ===!
!===----------------------------------------------------------------===!
!===    PART 1: READ PARAMETERS & EXTRACT DATA                      ===!
!===    PART 2: CALCULATION OF ATMOSPHERIC CONDITIONS               ===!
!===    PART 3: ITERATIVE METHODE TO RETRIEVE ABUNDANCES PROFILES   ===!
!===    PART 4: RADIATIVE TRANSFER                                  ===!
!===    PART 5: INVERSE METHODE                                     ===!
!===    PART 6: WRINTING OUTPUT                                     ===!
!======================================================================!
Program ABUNDANCE
  use OMP_LIB
  use lib_functions

  implicit none
  !----------------------------------------------------------------------!
  !--- PARAMETERS                                                     ---!
  !----------------------------------------------------------------------!
  integer, parameter:: dbl_type = SELECTED_REAL_KIND(15) !double precision

  character(*), parameter :: &
      !TODO: change to a dynamic variable
      path_input = "/home/mathe/Documents/ProgramProjects/CIRS-Titan-master/code_inversion/data/"

  logical, parameter :: &
      focal_plane = .True. !if True == fp3, False == fp4 (case CIRS data only)


  !----------------------------------------------------------------------!
  ! INTEGERS
  !----------------------------------------------------------------------!
  integer :: &
      nthread       = 0, & !number of thread use in parallel mode
      ithread       = 0, & !thread ID
      unit          = 0, & !unit file
      cpt           = 0, & !loop
      i             = 0, & !loop
      ic            = 0, & !loop on cloud
      ifreq         = 0, & !index on dimension nfreq_inv
      ifk           = 0, & !loop on kernel
      ik            = 0, & !loop on absorbers
      ikh           = 0, & !index where 'HAZE' is located on name_mol
      iprint        = 0, & !printing value
      iter          = 0, & !loop on iteration
      is            = 0, & !loop on secante
      ix            = 0, & !loop on spectral point
      i1            = 0, & !
      j             = 0, & !loop on layer
      ja            = 0, & !minimum layer reach (limb case)
      k             = 0, & !loop on molecule
      knu           = 0, & !number of frequencies on collisions
      lt            = 0, & !number of temperatures on collisions
      m             = 0, & !size of matrix C(m,m) and c0(m)
      nfcont        = 0, & !size of collision matrix
      nfreq_inv     = 0, & !size of convolution spectral range
      niter         = 0, & !number of iteration
      ninput        = 0, & !number of line on input temperature profile
      ninput_cl     = 0, & !number of line on input haze profile
      nbl           = 0, & !number of line
      nbl_noise     = 0, & !number of line on file noise
      nbl_sp_max    = 0, & !number max on spectroscopic file
      nb_cloud      = 0, & !number of cloud
      nb_df         = 0, & !number of invert spectral range
      nb_kernel     = 0, & !number of kernel on absorber
      nb_kernel_t   = 0, & !number of kernel on temperature
      nb_mol        = 0, & !number of molecule
      nb_mol_inv    = 0, & !number of molecule to retrieve (==absorber)
      nb_obs        = 0, & !number of secante
      nb_weightfunc = 0, & !number of weight function
      nlev          = 0, & !number of level on the atmospheric grid
      nlay          = 0, & !number of layer on the atmospheric grid
      nlor_ch4      = 0, & !number of different Lorentz halfwidths for CH4
      nond_max      = 0, & !number of line max on scattering file
      nvib_ch4      = 0, & !number of vibration mode on CH4
      np            = 0, & !number of line on file observations
      nsp           = 0, & !number of co-added spectra on each secante
      nx            = 0, & !size of invert spectral range
      n1            = 0, & !start index to extracting collisions data
      n2            = 0, & !end index to extracting collisions data
      time_begin    = 0.0d0, & !begin of the program (seconds)
      time_end      = 0.0d0, & !end of the program (seconds)
      time_precision= 0.0d0, & !time precision of the CPU (seconds)
      tmp_nlay      = 0, & !temporary nlay when atm grid is redefined
      tmp_nlev      = 0, & !temporary nlev when atm grid is redefined
      mode_inversion = 0   ! 0 = calcul direct, 1 = inversion molecule, 2 = inversion molecule + temperature

  integer, dimension(1,5) :: &
      ndeg_ch4 = 0 !number of degerescence on CH4

  integer, dimension(:)  , allocatable :: &
      icorps_k, & !array of absorbers (nb_mol_inv)
      lay_min , & !array where each secante reach minimal layer(nb_obs)
      nbl_sp  , & !array of number of line for each molecule(nb_mol)
      nlor    , & !number of different Lorentz halfwidths for each molecule(nb_mol)
      nond    , & !array of number of line for each cloud (nb_cloud)
      nvib    , & !array of number of vibration mode for each molecule(nb_mol)
      n_x         !array of number of invert spectral range

integer, dimension(:,:), allocatable :: &
    ndeg !array of number of degenerescence for each molecule(nb_mol,5)

!----------------------------------------------------------------------!
!--- REALS                                                          ---!
!----------------------------------------------------------------------!
real, dimension(:), allocatable :: &
    gnu, & !spectral array for collisions(knu)
    tmod   !temperature array for collisions(lt)

real, dimension(:,:), allocatable :: &
    ch4ch4_g, & !CH4-CH4 cross section collisions array(lt,knu)
    n2ch4_g , & !N2-CH4 cross section collisions array(lt,knu)
    n2h2_g  , & !N2-H2 cross section collisions array(lt,knu)
    n2n2_g      !N2-N2 cross section collisions array(lt,knu)

real(kind=dbl_type) :: &
    var_expo, & !variable for the exponent
    corr, &     !correlation lenght
    corr_b      !correlation lenght
!----------------------------------------------------------------------!
!--- DOUBLE PRECISION                                               ---!
!----------------------------------------------------------------------!
double precision :: &
    a_1bar             = 0.0d0, & !log of P_surface
    apl1               = 0.0d0, & !log of Pl_surface
    corr_len_t         = 0.0d0, & !correlation lengh of temperature
    corrT2             = 0.0d0, & !correlation lengh of temperature
    dT                 = 0.0d0, & !correlation lengh of temperature
    plo                = 0.0d0, & !correlation lengh of temperature
    tra_b              = 0.0d0, & !correlation lengh of temperature
    df_cont            = 0.0d0, & !spectral step for collisions
    dlnq               = 0.0d0, & !column density
    dteta              = 0.0d0, & !angle of secante
    d_1bar             = 0.0d0, & !
    f1_cont            = 0.0d0, & !first frequency for collisions
    f_end              = 0.0d0, & !last frequency for convolution
    f_pas              = 0.0d0, & !step frequency for convolution
    f_pas_mult         = 0.0d0, & !factor for line halfwidth
    f_start            = 0.0d0, & !first frequency for convolution
    gg                 = 0.0d0, & !
    g_1bar             = 0.0d0, & !Gravity at the surface
    h                  = 0.0d0, & !
    lat                = 0.0d0, & !latitude
    p_surf             = 0.0d0, & !Pressure at surface
    p_top              = 0.0d0, & !Pressure to the top
    qunsat             = 0.0d0, & !condensation distribution
    qsat               = 0.0d0, & !saturate distribution
    rayleigh           = 0.0d0, & !rayleigh
    rms                = 0.0d0, & !root mean square
    r1                 = 0.0d0, & !radii 1
    r2                 = 0.0d0, & !radii 2
    r3                 = 0.0d0, & !radii 3
    r_1bar             = 0.0d0, & !radii at 1 bar level
    time_elapsed       = 0.0d0, & !time spend of the program (seconds)
    cpu_time_begin     = 0.0d0, & !begin of the program (seconds CPU)
    cpu_time_elapsed   = 0.0d0, & !time spend of the program (seconds CPU)
    cpu_time_end       = 0.0d0, & !end of the program (seconds CPU)
    tra                = 0.0d0, & !trace of absorbers
    t_1bar             = 0.0d0, & !temperature at 1 bar level
    vmr_tropo_n2       = 0.0d0, & !N2 volume mixing ratio in troposphere
    vmr_tropo_h2       = 0.0d0, & !H2 volume mixing ratio in troposphere
    vmr_tropo_ar       = 0.0d0, & !Ar volume mixing ratio in troposphere
    vmr_tropo_ch4      = 0.0d0, & !CH4 volume mixing ratio in troposphere
    wei_t              = 0.0d0, & !weight for temperature
    gl_iter            = 0.d0, & ! use for gravity in layer
    x                  = 0.0d0    !

double precision, dimension(1) :: &
    a_ch4                , & ! ln(P) = A + B/T
    b_ch4                , & ! ln(P) = A + B/T
    cross_section_min_ch4, & !
    expo_rot_ch4         , & !
    mass_ch4             , & !
    temperature_test_ch4     !

double precision,dimension(5)      :: &
    f_cloud = 0.0d0 !frequency of the cloud

double precision,dimension(20)     :: &
    corr_len   = 0.0d0, & !correlation lenght for absorbers
    f_kernel   = 0.0d0, & !frequencies kernel for absorbers
    f_kernel_t = 0.0d0, & !frequencies kernel for temperature
    los        = 0.0d0, & !line of sight for each secante
    wei        = 0.0d0    !weight for each absorber

double precision, dimension(0:500) :: &
    et = 0.0d0 !

double precision, dimension(1,5) :: &
    alor_ch4 , & !freqencies from the line center of Lorentz profile
    elor_ch4 , & !T-exponent
    glor_ch4 , & !Lorentz halfwidth
    g2lor_ch4, & !
    vib_ch4      !

double precision,dimension(50,3) :: &
    df = 0.0d0 !array for invert spectral range, first frequency (col1),
               !last frequency (col2), step (col3)

double precision, dimension(400,100) :: &
    v = 0.0d0 !

double precision, dimension(:), allocatable :: &
    a                , & !ln(P) = A + B/T(nb_mol)
    alor             , & !(nb_mol)
    ap               , & !(nlev)
    api              , & !(ninput)
    apl              , & !(nlay)
    aql              , & !(nlay)
    ati              , & !(ninput)
    atmp_aql         , & !
    atmp_pl          , & !
    b                , & !ln(P) = A + B/T(nb_mol)
    cross_section_min, & !(nb_mol)
    c0               , & !(nx)
    errt             , & !(nlay)
    errtb            , & !(nlay)
    expo_rot         , & !(nb_mol)
    freq_x           , & !(nf)
    gl               , & !(nlay)
    qch4             , & !(nlay)
    qref             , & !
    mass             , & !mass for each molecule (nb_mol)
    ml               , & !(nlay)
    p                , & !(nlev)
    pii              , & !(ninput)
    pl               , & !(nlay)
    rad_noise        , & !(nf)
    robs             , & !
    rnoise           , & !
    t                , & !(nlev)
    tab_i            , & !(nlay)
    taucl_tot        , & !(nlay)
    temperature_test , & !(nb_mol)
    ti               , & !(ninput)
    tl               , & !(nlay)
    tmp_ap           , & !(tmp_nlev)
    tmp_apl          , & !(tmp_nlay)
    tmp_aql          , & !(tmp_nlay)
    tmp_gl           , & !(tmp_nlay)
    tmp_ml           , & !(tmp_nlay)
    tmp_p            , & !(tmp_nlay)
    tmp_pl           , & !(tmp_nlay)
    tmp_qch4         , & !(tmp_nlay)
    tmp_t            , & !(tmp_nlay)
    tmp_tl           , & !(tmp_nlay)
    tmp_z            , & !(tmp_nlev)
    wobs             , & !
    wnoise           , & !
    z                    !altitude array(nlev)

double precision, dimension(:,:), allocatable :: &
    apcli     , & !(ninput)
    ataucli   , & !(ninput)
    c         , & !(nx,nx)
    cf        , & !(nlay,nlay)
    cf_b      , & !(nlay,nlay)
    dlnerrq   , & !(nlay,nb_mol)
    dlnerrqb  , & !(nlay,nb_mol)
    elor      , & !(nb_mol,5)
    errqb_sup , & !(nlay,nb_mol)
    errqb_inf , & !(nlay,nb_mol)
    errq_sup  , & !(nlay,nb_mol)
    errq_inf  , & !(nlay,nb_mol)
    e_input   , & !
    glor      , & !(nb_mol,5)
    g2lor     , & !(nb_mol,5)
    g_input   , & !
    matrix_r  , & !(nlay,nx)
    matrix_r2 , & !(nlay,nx)
    pcli      , & !(ninput)
    qex       , & !(nlay,nlay)
    ql        , & !(nlay,0:nb_mol)
    rad       , & !synthetic radiance(nb_obs,nfreq_inv)
    rad_diff  , & !rad_obs -rad(nb_obs,nf)
    rad_obs   , & !observed radiance(nb_obs,nf)
    sec       , & !(nlay,nb_obs)
    sig       , & !(nlay,nlay)
    sm        , & !(nlay,nx)
    ssk_b     , & !(nx,nx)
    s_input   , & !
    s_b       , & !(nlay,nlay)
    tab_ch4ch4, & !CH4-CH4 collisions array used(nlay,nfcont)
    tab_n2ch4 , & !N2-CH4 collisions array used(nlay,nfcont)
    tab_n2h2  , & !N2-H2 collisions array used(nlay,nfcont)
    tab_n2n2  , & !N2-N2 collisions array used(nlay,nfcont)
    taucl     , & !(4,nlay)
    taucli    , & !(ninput)
    tmp_ql    , & !(tmp_nlay,0:nb_mol)
    vib       , & !frequencies of vibrational mode for each molecule(nb_mol,5)
    w_input       !frequencies array for spectroscopic data (nbl_mol,nbl_sp_max)

double precision, dimension(:,:,:), allocatable :: &
    s  , & !matrix s(nlay,nlay,nb_mol)
    sk , & !matrix sk(nb_mol_inv,nlay,nx)
    ssk    !matrix ssk(nb_mol_inv,nx,nx)

double precision, dimension(:,:,:,:), allocatable :: &
    matrix_k !matrix for absorbers (nb_mol_inv,nlay,nb_obs,nfreq_inv)

  double precision, dimension(:,:,:), allocatable :: &
    matrix_t    !matrix for temperature (nlay,nb_obs,nfreq_inv)

!----------------------------------------------------------------------!
! CHARACTERS
!----------------------------------------------------------------------!
character(len=1)    :: &
    convfunc = '' !Convolution function used

character(len=256) :: &
    filename_input   = '', & !input parameter file
    file_ch4ch4      = '', & !CH4-CH4 collision file
    file_kernel      = '', & !output abundance kernel file (output)
    file_kernel_t    = '', & !output temperature kernel file (output)
    file_n2ch4       = '', & !N2-CH4 collision file
    file_n2h2        = '', & !N2-H2 collision file
    file_n2n2        = '', & !N2-N2 collision file
    file_noise       = '', & !noise file
    file_profq       = '', & !output abundance profile (output)
    file_prof_t      = '', & !output temperature profile (output)
    file_spectro_ch4 = '', & !spectroscopic file for CH4
    file_spec_synth  = '', & !synthetic spectrum file (output)
    file_weightfunc  = '', & !weight function file
    name_pla         = '', & !name of the planet
    prof_t           = ''    !input thermal profile

character(len=256), dimension(5) :: &
    file_scat = '', & !input scattering profile
    prof_haze = ''    !input haze profile

character(len=4),dimension(20) :: &
    name_mol_inv = '' !name for each absorbers

character(len=256),dimension(20) :: &
    file_obs = '' !observation files

character(len=256),dimension(20,3) :: &
    name_mol = '' !name (col 1), type profile (col 2), input profile (col 3) for each molecule

character(len=256), dimension(:), allocatable :: &
    file_spectro !spectroscopic files for each molecule (nb_mol)
!----------------------------------------------------------------------!
! LOGICALS
!----------------------------------------------------------------------!
  logical :: &
    limbe     , & !if at least one limb == True
    method_invert !inverse matrix C, Cholesky or Gauss method
  !====================================================================================================================!
  !=== BEGIN PROGRAM                                                                                                ===!
  !====================================================================================================================!
  call system_clock(count=time_begin, count_rate=time_precision)
  call CPU_TIME(cpu_time_begin)
  write(*,*)"============================================================"
  write(*,*)"BEGIN PROGRAM ABUNDANCE"
  write(*,*)"============================================================"
  !$OMP PARALLEL PRIVATE(i)
    nthread = OMP_GET_NUM_THREADS()
    write(*,*)'The work is done with, ', nthread, 'threads.'
  !$OMP END PARALLEL
  !====================================================================================================================!
  !=== PART 1: READ PARAMETERS & EXTRACT DATA                                                                       ===!
  !====================================================================================================================!
  !--- Section 1.1: Read file parameter                                                                             ---!
  !--------------------------------------------------------------------------------------------------------------------!
  call GET_COMMAND_ARGUMENT(1, filename_input)
  write(*,*)"Read file parameters ..."
  call read_parameters(filename_input, name_pla, lat, nlev, p_surf, p_top, prof_t, nb_cloud, f_cloud, prof_haze, &
                       file_scat, vmr_tropo_n2, vmr_tropo_h2, vmr_tropo_ar, vmr_tropo_ch4, file_n2n2, file_n2ch4, &
                       file_ch4ch4, file_n2h2, nb_mol, name_mol, nb_obs, los, file_obs, file_noise, f_start, f_end, &
                       f_pas, f_pas_mult, convfunc, rayleigh, nb_weightfunc, file_weightfunc, nb_df, df, nb_mol_inv, &
                       name_mol_inv, corr_len, wei, niter, iprint, nb_kernel, f_kernel, file_spec_synth, file_profq, &
                       file_kernel, limbe, nsp, nb_kernel_t, f_kernel_t, corr_len_t, wei_t, file_prof_t, &
                       file_kernel_t, mode_inversion)
  if (mode_inversion == 0) then
    niter = 0 ! calcul direct does not need iteration
  end if
  write(*,*)"Read file parameters done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.2: Allocate variables                                                                              ---!
  !--------------------------------------------------------------------------------------------------------------------!
  nlay = abs(nlev) - 1
  allocate(file_spectro(nb_mol), cross_section_min(nb_mol), temperature_test(nb_mol), nlor(nb_mol), alor(nb_mol), &
           glor(nb_mol,5), elor(nb_mol,5), g2lor(nb_mol,5), nvib(nb_mol), vib(nb_mol,5), ndeg(nb_mol,5), p(abs(nlev)), &
           ap(abs(nlev)), pl(nlay), apl(nlay), mass(nb_mol), a(nb_mol), b(nb_mol), expo_rot(nb_mol))
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.8: Extract planet data                                                                             ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Searching planet parameters in the database_planet.f90"
  write(*,*)"Case: ",trim(name_pla)
  call database_planet(name_pla,lat,g_1bar,r_1bar)
  write(*,fmt='(TR5,A,F6.2,A)')    'G(surface): ', g_1bar, ' cm.sec-2'
  write(*,fmt='(TR5,A,TR2,F7.2,A)')'R(1bar): ', r_1bar, ' km'
  write(*,fmt='(TR5,A,TR7,F6.2,A)')'Lat: ', lat, ' deg'
  write(*,*)"Input tropospheric mixing ratios"
  write(*,fmt='(A,TR2,F5.4)')'N2: ', vmr_tropo_n2
  write(*,fmt='(A,TR2,F5.4)')'H2: ', vmr_tropo_h2
  write(*,fmt='(A,TR2,F8.7)')'Ar: ', vmr_tropo_ar
  write(*,fmt='(A,TR1,F5.4)')'CH4: ', vmr_tropo_ch4
  write(*,*)"Done."
!----------------------------------------------------------------------------------------------------------------------!
!--- Section 1.9: Extract molecule data                                                                             ---!
!---              Selection is made in subroutine READ_LINE                                                         ---!
!---                                                                                                                ---!
!---    NLOR: Number of different Lorentz halfwidths for gas K                                                      ---!
!---    ALOR: Lorentz profile is computed up to ALOR cm-1 from line center                                          ---!
!---    GLOR: Lorentz halfwidths (1 --> NLOR)                                                                       ---!
!---    ELOR: T-exponent (1 --> NLOR)                                                                               ---!
!---    G2LOR: GLOR is used for lines having G2LOR(i-1)< G < G2LOR(i) (296 K)                                       ---!
!---    NVIB: Number of vibrational modes to be considered for calculation of the partition function                ---!
!---    VIB, NDEG: Frequency and degeneracy of the vibrational modes (1 --> NVIB)                                   ---!
!----------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Searching molecules parameters in the database_mol.f90"
  do ik = 1, nb_mol
    write(*,*)trim(name_mol(ik,1)), " found."
  end do
  call database_mol(nb_mol, name_mol, focal_plane, nlor, nvib, ndeg, a, b, mass, expo_rot, cross_section_min, &
                    temperature_test, alor, vib, glor, elor, g2lor, file_spectro)
  write(*,*)"Done."
  do k = 1, nb_mol
    if(nlor(k) > 5) then
        write(*,fmt='(a,i3,a,a4,a)')' Molecule #', k,' (', name_mol(k,1), ') must not have more than 5 line halfwidths'
        exit
    end if
    if(nvib(k) > 5) then
        write(*,fmt='(a,i3,a,a4,a)')' Molecule #', k, ' (', name_mol(k,1), ') must not have more than 5 vibration modes'
      exit
    end if
  end do
  !====================================================================================================================!
  !=== PART 2: CALCULS ATMOSPHERIC PARAMETERS                                                                       ===!
  !====================================================================================================================!
  write(*,*)"=========================="
  write(*,*)"PART 2: Calculs atmospheric parameters"
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.1: Calcul of pressure grid level                                                                   ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Input pressure grid level (mbar)"
  if(nlev > 0) then ! attention pas pris en compte cette partie a cause du p(j) lu en entree
    method_invert = .FALSE.
    do j=1,nlev
        read(5,*) p(j)
        ap(j) = dlog(p(j))
    end do
  else
    method_invert = .TRUE.
    nlev = -nlev
    p(1) = p_surf
    ap(1) = dlog(p_surf)
    do j = 2, nlev
      p(j) = p_surf * (p_top/p_surf)**(dfloat(j-1)/dfloat(nlev-1))
      ap(j) = dlog(p(j))
    end do
  end if
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.2: Calcul of the pressure in the layers                                                            ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Calculation of the pressure in the layers (mbar)"
  do j = 1, nlay
    pl(j) = dsqrt(p(j)*p(j+1))
    apl(j) = dlog(pl(j))
  end do
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.3: Extract temperature profil                                                                      ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Input temperature profile: "
  write(*,fmt='(tr5,a,a)')'target file: ', prof_t
  ninput = nb_line(prof_t)
  allocate(pii(ninput), ti(ninput), api(ninput), ati(ninput))
  open(unit=10, file=prof_t, status='old')
  do j = 1, ninput
    read(10,*) pii(j), ti(j)
  end do
  close(10)
  if (pii(1) < pii(ninput)) then !si p(1) = sol et p(ninput) = ciel
    allocate(tmp_p(ninput), tmp_t(ninput))
    tmp_p = pii
    tmp_t = ti
    do i= 1, ninput
        pii(i) = tmp_p(ninput+1-i)
        ti(i) = tmp_t(ninput+1-i)
    end do
    deallocate(tmp_p, tmp_t)
  end if
  api = dlog(pii)
  ati = dlog(ti)
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.4: Calcul the temperature in the layers                                                            ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Calculation of the temperature in the layers (K)"
  allocate(t(nlev), tl(nlay))
  do j = 1, nlay
    t(j) = dexp(tina(ap(j), api, ati, ninput))
    tl(j) = dexp(tina(apl(j), api, ati, ninput))
  end do
  t(nlev) = dexp(tina(ap(nlev), api, ati, ninput))
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.5: Calcul the CH4 mixing ratio in the layers                                                       ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Calculation of the CH4 mixing ratio in the layers"
  allocate(qch4(nlay))
  call database_mol(1, 'CH4', focal_plane, nlor_ch4, nvib_ch4, ndeg_ch4, a_ch4, b_ch4, mass_ch4, expo_rot_ch4, &
                    cross_section_min_ch4, temperature_test_ch4, alor_ch4, vib_ch4, glor_ch4, elor_ch4, g2lor_ch4, &
                    file_spectro_ch4)
  qsat = 1.01325d+03 * 10**( a_ch4(1) - b_ch4(1)/tl(1) )/pl(1)
  qch4(1) = dmin1(vmr_tropo_ch4, qsat)
  do j = 2, nlay
    qsat = 1.01325d+03 * 10**( a_ch4(1)-b_ch4(1)/tl(j) )/pl(j)
    qch4(j) = dmin1(qch4(j-1), qsat)
  end do
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.6: Calcul the molecular weight in the layers                                                       ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Calculation of the molecular weight in the layers"
  allocate(ml(nlay))
  do j=1,nlay
    ml(j) = 16.043*qch4(j) + (28.0134*vmr_tropo_n2 + 2.0158*vmr_tropo_h2 + 35.9675*vmr_tropo_ar) * &
            (1.-qch4(j))/(1.-vmr_tropo_ch4)
  end do
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.7: Calcul the altitude levels and gravity in the layers                                            ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Calculation of the altitude levels (km) and the gravity in the layers (cm sec-2)"
  allocate(tab_i(nlay))
  do j = 1, nlay
    tab_i(j) = j
  end do
  a_1bar = dlog(p(1))
  d_1bar = tina(a_1bar, apl, tab_i, nlay)
  t_1bar = tina(a_1bar, apl, tl, nlay)
  if(d_1bar < 0.4999d+00) then
    write(*,*)' WARNING: the 1-bar level is not included in the grid'
  end if
  i1 = idnint(d_1bar+0.0001)
  allocate(z(nlev), gl(nlay))
  if(i1 > 1) then
    h = 1.d+02 * r * 0.5 * (t_1bar+tl(i1))/(ml(i1)*g_1bar)
    z(i1) = h * (a_1bar-ap(i1))
    do j= i1, 2, -1
        gg = g_1bar * (r_1bar/(r_1bar+z(j)))**2
        h = 1.d+02 * r * tl(j-1)/(ml(j-1)*gg)
        z(j-1) = z(j) + h * (ap(j)-ap(j-1))
        gl(j-1) = g_1bar * (r_1bar/(r_1bar+0.5*(z(j)+z(j-1))))**2
    end do
  else
    i1 = 1
    h = 1.d+02 * r * 0.5 * (t_1bar+tl(1))/(ml(1)*g_1bar)
    z(1) = h * (a_1bar-ap(1))
  end if
  ! methode iterative pour calculer g au milieu de couche et estimer au mieux les altitudes des niveaux.
  do j = i1, nlay
    gg = g_1bar * (r_1bar/r_1bar+z(j))**2
    h = 1.d+02 * r * tl(j)/(ml(j)*gg)
    z(j+1) = z(j) + h * (ap(j)-ap(j+1))
    gl(j) = g_1bar * ( r_1bar / (r_1bar + 0.5 * (z(j)+z(j+1)) ) )**2
    gl_iter = 0d0
    do while(abs(gl(j)-gl_iter) >= 1e-12)
      gl_iter = gl(j)
      h = 1.d+02 * r * tl(j)/(ml(j)*gl(j))
      z(j+1) = z(j) + h * (ap(j)-ap(j+1))
      gl(j) = g_1bar * ( r_1bar / (r_1bar + 0.5 * (z(j)+z(j+1)) ) )**2
    end do
  end do
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.8: Printing the parameters for the atmospheric levels and layers                                   ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Printing the parameters for the atmospheric levels and layers"
  write(*,*)""
  write(*,fmt='(5x,a,5x,a,3x,a,2x,a,4x,a,5x,a)')'Altitude', 'Pressure', 'Temperature', 'Mol. weight', 'q(ch4)', &
                                                'Gravity'
  write(*,fmt='(8x,a,8x,a,8x,a,10x,a,18x,a)')'(km)', '(mbar)', '(K)', '(g)', '(cm.sec-2)'
  do j= 1, nlay
    write(*,fmt='(i4,f9.2,e14.4)')j, z(j), p(j)
    write(*,fmt='(15x,e14.4,f9.2,f13.3,e15.4,f10.2)')pl(j), tl(j), ml(j), qch4(j), gl(j)
  end do
  write(*,fmt="(i4,f9.2,e14.4)")nlev, z(nlev), p(nlev)
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.9: Interpolation of molecular abundances profils                                                   ---!
  !---    name_mol(ik,2) = 1 : vertical profile given in NIN points                                                 ---!
  !---    name_mol(ik,2) = 2 : saturated distribution above condensation level QUNSAT is the mole fraction below    ---!
  !---                         condensation level                                                                   ---!
  !---    name_mol(ik,2) = 3 : saturated distribution below condensation level QUNSAT is the mole fraction above    ---!
  !---                         condensation level                                                                   ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Input: Vertical profiles  of the molecular absorbers"
  write(*,*)""
  allocate(ql(nlay,0:nb_mol))
  !$OMP PARALLEL PRIVATE(tmp_pl, tmp_aql, atmp_pl, atmp_aql, unit, ithread,  nbl)
  !$OMP DO SCHEDULE(RUNTIME)
  do ik = 1, nb_mol
    !$ ithread = OMP_GET_THREAD_NUM()
    !$ unit = 10 * ithread + ik
    write(*,*)'----------------------------------'
    write(*,*)trim(NAME_MOL(ik,1)), ': ', trim(file_spectro(ik))
    select case (NAME_MOL(ik,2))
        case('1')
            nbl = nb_line(trim(NAME_MOL(ik,3)))
            allocate(tmp_pl(nlay), tmp_aql(nlay), atmp_pl(nlay), atmp_aql(nlay))
            open(unit=unit, file=trim(NAME_MOL(ik,3)), status='old')

            do j = 1, nbl
                read(unit,*) tmp_pl(j), tmp_aql(j)
                atmp_pl(j) = dlog(tmp_pl(j)) 
                atmp_aql(j) = dlog(tmp_aql(j)) 
            end do
            if (nbl <= nlay) then
              do j = nbl+1, nlay
                atmp_pl(j) = atmp_pl(j-1)
                atmp_aql(j) = atmp_aql(nbl)
              end do
            end if
            do j = 1, nlay
              ql(j,ik) = dexp(tina(apl(j), atmp_pl, atmp_aql, nlay))
            end do
            deallocate(tmp_pl,tmp_aql,atmp_pl,atmp_aql)
            close(unit)
        case('2')
            read(NAME_MOL(ik,3),*)qunsat ! convert string NAME_MOL(ik,3) to float
            qsat = 1.01325d+03 * 10**( a(ik)-b(ik)/tl(1) ) / pl(1)
            ql(1,ik)=dmin1(qsat,qunsat)
            do j = 2, nlay
                qsat = 1.01325d+03 * 10**( a(ik)-b(ik)/tl(j) ) / pl(j)
                ql(j,ik) = dmin1(ql(j-1,ik), qsat)
            end do
        case('3')
            read(NAME_MOL(ik,3),*)qunsat
            qsat = 1.01325d+03 * 10**( a(ik)-b(ik)/tl(nlay) ) / pl(nlay)
            ql(nlay,ik) = dmin1(qsat, qunsat)
            do j = nlay-1, 1, -1
                qsat = 1.01325d+03 * 10**( a(ik)-b(ik)/tl(j) ) / pl(j)
                ql(j,ik) = dmin1(ql(j+1,ik), qsat)
            end do
        case default
            write(*,*)"Type '1': input on file"
            write(*,*)"Type '2': saturated distribution above condensation level"
            write(*,*)"Type '3': saturated distribution below condensation level"
            stop
    end select
    if(name_mol(ik,1) == 'CH4') then
        write(*,*)"--------------------------"
        write(*,*)"Changing qch4 to ql(j,numch4) if CH4 is an absorber"
        do j = 1, nlay
            qch4(j) = ql(j,ik)
        enddo
    end if
    write(*,fmt='(10x,a,i3,a,a4,a,f6.2,a)')'Absorber #', ik, ': ', name_mol(ik,1), ' Mol.weight =', mass(ik), ' g'
    if(name_mol(ik,2) == '1') then
        write(*,fmt='(10x,a)')trim(name_mol(ik,3))
    else
        if (name_mol(ik,2) == '2') then
            write(*,fmt='(10x, a, e11.4)')'Saturated profile for P<P_cond: unsaturated mixing ratio =', qunsat
        else
            write(*,fmt='(10x, a, e11.4)')'Saturated profile for P>P_cond: unsaturated mixing ratio =', qunsat
        endif
    endif
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  write(*,*)"Done."
  write(*,*)"--------------------------"
  write(*,*)"Printing the mole fractions in the atmospheric layers"
  write(*,fmt='(10x,a,10x,a,20(6x,a4,1x))')'P', 'T', (trim(name_mol(i,1)), i=1, nb_mol)
  do j=1,nlay
    write(*,fmt='(i4,e12.4,f8.2,20(e11.3))')j, pl(j), tl(j), (ql(j,i), i=1, nb_mol)
  end do
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.10: Geometry of the observations                                                                   ---!
  !---    Calculation of the airmasses in the layers (nb_obs>0)                                                     ---!
  !---    nb_obs>=1 : vertical viewing (nb_obs points; input airmass pertains to the 1-bar level) (starting airmass ---!
  !---                and step in airmass are given). Transmittances are first calculated for an airmass of los.    ---!
  !---   nb_obs=0  : radiances are calculated with the 2nd order integral                                           ---!
  !---   nb_obs<0  : horizontal viewing (-nb_obs points) (starting altitude -km- and step in altitude are given) If ---!
  !---               los>0, the pressure grid is redefined. Transmittances are first calculated for altitude in the ---!
  !---               Line of Sight (los)                                                                            ---!
  !--------------------------------------------------------------------------------------------------------------------!
  if (.not. limbe) then
    write(*,fmt='(a,i3,a,f8.4,a,f8.4)')'Spectral calculations for', nb_obs, ' airmasses, from', los(1),' to', los(nb_obs)
    allocate(sec(nlay,nb_obs), lay_min(nb_obs))
    lay_min = 0
    do i = 1, nb_obs
      do j = 1, nlay
        r2 = r_1bar + z(j)
        r3 = r_1bar + z(j+1)
        dteta = dasin( dsin( dacos(1./los(i)) )*r_1bar/r3 ) - dasin( dsin( dacos(1./los(i)) )*r_1bar/r2 )
        sec(j,i) = dsqrt( (r3**2 + r2**2 - 2d0*r2*r3*dcos(dteta)) ) / (r3-r2)
      end do
    end do
  else
    if(nb_obs == 0) then
      write(*,*)'Spectral calculations of the disk-averaged radiance (plan parallel geometry)'
      nb_obs = 1
      allocate(sec(nlay,1))
      do j = 1, nlay
        sec(j,1) = 1d0
      end do
    else
      write(*,fmt='(a,i3,a,f10.2,a,f10.2,a)')'Spectral calculations for', nb_obs, ' altitudes, from', los(1), 'to', &
                                             los(nb_obs), 'km'
      if(los(1) >= z(1)) then
        !** Redefinition of the pressure grid (nb_obs<0)
        allocate(tmp_p(nlev), tmp_ap(nlev), tmp_t(nlev), tmp_z(nlev), tmp_pl(nlay), tmp_tl(nlay), tmp_apl(nlay), &
                 tmp_qch4(nlay), tmp_ml(nlay), tmp_gl(nlay), tmp_ql(nlay,0:nb_mol), tmp_aql(nlay))
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
        j = 2
        do while (los(1) >= tmp_z(j))
          j = j + 1
        end do
        ja = j - 1
        if (j == nlev) then
          write(218,fmt='(a,f8.2,a)')'Viewing altitude', los(1),' is above the atmospheric grid'
          stop
        end if
        tmp_nlev = nlev ! ici le marker 27
        tmp_nlay = nlay
        nlev = nlev - ja + 1
        nlay = nlev - 1
        allocate(p(nlev), ap(nlev), z(nlev), t(nlev), tl(nlay), pl(nlay), apl(nlay), qch4(nlay), ml(nlay), gl(nlay), &
                 ql(nlay,0:nb_mol), aql(tmp_nlay))
        p(1) = dexp(tina(los(1), tmp_z, tmp_ap, tmp_nlev))
        t(1) = dexp(tina(dlog(p(1)), api, ati, ninput))
        z(1) = los(1)
        write(*,fmt='(a,f8.2,a,e11.4,a)')'Atmospheric grid is redefined from altitude =', los(1), 'km: P =', p(1), &
                                          ' mbar'
        ap(1) = dlog(p(1))
        pl(1) = dsqrt( p(1)*tmp_p(ja+1) )
        apl1  = dlog(pl(1))
        tl(1) = dexp( tina(apl1, api, ati, ninput) )
        ml(1) = tina(apl1, tmp_apl, tmp_ml, tmp_nlay)
        gl(1) = tina(apl1, tmp_apl, tmp_gl, tmp_nlay)
        do j = 1, tmp_nlay
            aql(j) = dlog( tmp_qch4(j) )
        end do
        qch4(1) = dexp( tina(apl1, tmp_apl, aql, tmp_nlay) )
        apl(1) = apl1
        tmp_apl(1) = apl1
        do ik = 1, nb_mol
            do j = 1, tmp_nlay
                aql(j) = dlog(tmp_ql(j,ik))
            end do
            ql(1,ik) = dexp( tina(apl1, tmp_apl, aql, tmp_nlay) )
        end do
        do j = 2, nlev
            z(j)  = tmp_z(j+ja-1)
            p(j)  = tmp_p(j+ja-1)
            t(j)  = tmp_t(j+ja-1)
            ap(j) = tmp_ap(j+ja-1)
        end do
        do j = 2, nlay
          pl(j)   = tmp_pl(j+ja-1)
          apl(j)  = tmp_apl(j+ja-1)
          tl(j)   = tmp_tl(j+ja-1)
          ml(j)   = tmp_ml(j+ja-1)
          gl(j)   = tmp_gl(j+ja-1)
          qch4(j) = tmp_qch4(j+ja-1)
          do ik = 1, nb_mol
            ql(j,ik) = tmp_ql(j+ja-1,ik)
          end do
        end do
        ql(:,0) = 0.
        deallocate(tmp_p, tmp_ap, tmp_t, tmp_z, tmp_pl, tmp_tl, tmp_apl, tmp_qch4, tmp_ml, tmp_gl, tmp_ql, tmp_aql)
      end if
      ! Calculation of the airmasses in the layers (nb_obs<0) lay_min(i) is the layer where altitude LOS is reached
      allocate(lay_min(nb_obs),sec(nlay,nb_obs))
        sec = 0.
        do i = 1, nb_obs
          lay_min(i) = nlev
          r1 = r_1bar + los(i)
          if(los(i) < z(1)) then
                lay_min(i) = 0
                write(*,fmt='(a,f10.2,a)')'Altitude', los(i), ' km reached below first layer'
          endif
          do j = 1, nlay
            r2 = r_1bar + z(j)
            r3 = r_1bar + z(j+1)
            if(los(i) < z(j)) then
              sec(j,i) = ( dsqrt(r3*r3 - r1*r1) - dsqrt(r2*r2 - r1*r1) ) / (r3-r2)
            else
              if(los(i) < z(j+1)) then
                sec(j,i) = dsqrt(r3*r3 - r1*r1) / (r3-r2)
                lay_min(i) = j
                write(*,fmt='(a,f10.2,a,i3,a,e11.4,a)')'Altitude', los(i), ' km reached in layer #', lay_min(i), &
                                                       ', at pressure layer ',pl(j), ' mbar'
              end if
            end if
          end do
          if(lay_min(i) > nlay) then
            nb_obs = i-1
            write(*,fmt='(a,f10.2,a,i3,a)')' WARNING: Altitude', los(i), ' km is above the atmospheric grid;'
            write(*,fmt='(5x,a,i3,a)')'only', nb_obs, ' altitudes are hereafter considered'
            exit !goto 36
          endif
      end do
      if(lay_min(nb_obs) == 0) then
        nb_obs = -nb_obs
      end if
    end if
  end if
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.12: Calculation of et, nfreq_inv                                                                   ---!
  !--------------------------------------------------------------------------------------------------------------------!
  if(nb_obs == 0) then
    do i = 0, 500
        x = 0.02 * dfloat(i)
        et(i) = dexp(-x) - x*eint(-x)
    end do
  else
    do i = 0, 500
        et(i) = dexp(-0.02*dfloat(i))
    end do
  end if
  !--- Spectral range (cm-1)
  nfreq_inv = nint( (F_end-F_start) / F_pas ) + 1
  write(*,fmt='(a,f9.3,a,f9.3,a)')'Spectrum if computed from', F_start, ' to', F_end, ' cm-1'
  write(*,fmt='(5x,a,f3.1,a)')'Calculation step is ', F_pas_mult,' times the smallest line halfwidth'
  write(*,fmt='(5x,a,f5.3,a)')'Output step is ', F_pas,' cm-1'
  !--- Output file, convolution parameters, weighting functions
  write(*,fmt='(a,a)')'Output file is: ', trim(file_spec_synth)
  write(*,fmt='(a,a1,a,f6.4,a)')'Convolution type is "', convfunc, '", resolution is ', rayleigh, ' cm-1'
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.13: Parameters for cloud opacity (nb_cloud clouds)                                                 ---!
  !---    TAUCLT: Total optical depth (LIMBE) at wavenumber f_cloud (cm-1)                                          ---!
  !---    PBOT, PTOP: Boundaries of the cloud (mbar)                                                                ---!
  !---    HPHG: Particle-to-gas scale height                                                                        ---!
  !--------------------------------------------------------------------------------------------------------------------!
  allocate(taucl_tot(nlay),taucl(nb_cloud,nlay))
  if(nb_cloud > 0) then
    do j = 1, nlay
      taucl_tot(j) = 0.
    end do
    do ic = 1, nb_cloud
      do j = 1, nlay
        taucl(ic,j) = 0.
      end do
      ninput_cl = nb_line(PROF_HAZE(ic)) ! here file unit = 400
      allocate(pcli(ic,ninput_cl), taucli(ic,ninput_cl), apcli(ic,ninput_cl), ataucli(ic,ninput_cl))
      open(unit=60,file=PROF_HAZE(ic),status='old')
      do j = 1, ninput_cl
        read(60,*)pcli(ic,j), taucli(ic,j)
        apcli(ic,j) = dlog(pcli(ic,j))
        ataucli(ic,j) = dlog(taucli(ic,j))
      end do
      close(60)
      write(*,*)'ninput_cl = ', ninput_cl
      do j = 1, ninput_cl
        write(*,fmt='(a,i3,5x,e11.4,x,a,x,e11.4)')'j = ', j, pcli(ic,j), '-', taucli(ic,j)
      end do
      do j = 1, nlay
        taucl(ic,j) = dexp( tina(apl(j), apcli, ataucli, ninput_cl) )
        taucl_tot(j) = taucl_tot(j) + taucl(ic,j)
        write(*,fmt='(a,i3,5x,e11.4,x,a,x,e11.4,5x,e11.4)')'j = ', j, p(j), '-', p(j+1), taucl(ic,j)
      end do
    end do
  end if
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.4: Spectral range for which K matrices are calculated                                              ---!
  !--------------------------------------------------------------------------------------------------------------------!
  allocate(n_x(nb_df))
  do i=1,nb_df
    n_x(i) = nint( (df(i,2)-df(i,1)) / df(i,3) ) + 1
    nx = nx + n_x(i)
  end do
  write(*,*)'nx = ', nx
  allocate(freq_x(nx))
  cpt = 1
  do j = 1, nb_df
    do i = 1, n_x(j)
        freq_x(cpt) = df(j,1)  + df(j,3)*(i-1)
        cpt = cpt + 1
    end do
  end do
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.9: Extract observation data                                                                        ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Open files observed"
  write(*,*)'Files of observed spectra:'
  write(*,fmt='(10x,a)')(trim(file_obs(is)),is=1,nb_obs)
  allocate(rad_obs(nb_obs,nx))
  do is = 1, nb_obs
    open(unit=10, file=file_obs(is), status='old')
    read(10,*)
    read(10,*)np
    read(10,*)
    allocate(wobs(np), robs(np))
    do i = 1, np
      read(10,*)wobs(i), robs(i)
    end do
    do i = 1, nx
      rad_obs(is,i) = 1.e+07 * tina(freq_x(i), wobs, robs, np)
    end do
    deallocate(wobs,robs)
  end do
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.10: Extract noise data                                                                             ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Read file noise"
  write(*,fmt='(tr5,a,a)')'File of NESRs: ', trim(path_input)//trim(file_noise)
  write(*,fmt='(tr5,a,i3)')'Number of spectrum: ', nsp
  nbl_noise = nb_line(trim(path_input)//file_noise) - 1 ! Here file unit is 400
  open(unit=10, file=trim(path_input)//file_noise, status='old')
  allocate(wnoise(nbl_noise), rnoise(nbl_noise))
  read(10,*) ! first line is the title
  do i = 1, nbl_noise
    read(10,*)wnoise(i), rnoise(i)
  end do
  close(10)
  allocate(rad_noise(nx))
  do i = 1, nx
    rad_noise(i) = 1.e+07 * tina(freq_x(i), wnoise, rnoise, nbl_noise) / dsqrt(dfloat(nsp))
  end do
  write(*,fmt='(a,tr2,f8.3,a,tr2,e12.3)')'Freq_x(1): ',freq_x(1),', rad_noise(1): ',rad_noise(1)
  write(*,fmt='(a,tr1,f8.3,a,tr1,e12.3)')'Freq_x(nx): ',freq_x(nx),', rad_noise(nx): ',rad_noise(nx)
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.6: Extract spectroscopic file H2-H2, H2-He and H2-CH4                                              ---!
  !---              f1_cont is the first frequency                                                                  ---!
  !---              df_cont is the step                                                                             ---!
  !---              lt temperatures                                                                                 ---!
  !---              knu frequencies                                                                                 ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Collisions"
  ! N2N2 collisions
  open(unit=12, file=trim(path_input)//FILE_N2N2, status='old', form='unformatted', access='stream')
  read(12,pos=5)lt, knu
  close(12)
  allocate(n2n2_g(lt,knu), tmod(lt), gnu(knu))
  open(unit=12,file=trim(path_input)//FILE_N2N2, status='old', form='unformatted')
  read(12)lt, knu, (tmod(i),i=1,lt), (gnu(i),i=1,knu), ((n2n2_g(i,j),i=1,lt),j=1,knu)
  close(12)
  write(*,fmt='(a,a)')'N2-N2 spectroscopic file :', trim(file_n2n2)
  write(*,fmt='(5x,a, f8.2,a,f8.2,a,i3,a,i4,a)')'from', gnu(1), ' to', gnu(knu),' cm-1: ', lt, ' temperatures', knu, &
                                                ' frequencies'
  ! N2CH4 collisions
  allocate(n2ch4_g(lt,knu))
  open(unit=12, file=trim(path_input)//FILE_N2CH4, status='old', form='unformatted')
  read(12)lt, knu, (tmod(i),i=1,lt), (gnu(i),i=1,knu), ((n2ch4_g(i,j),i=1,lt),j=1,knu)
  close(12)
  write(*,fmt='(a,a)')'N2-CH4 spectroscopic file :', trim(file_n2ch4)
  write(*,fmt='(5x,a, f8.2,a,f8.2,a,i3,a,i4,a)')'from', gnu(1), ' to', gnu(knu),' cm-1: ', lt, ' temperatures', knu, &
                                                ' frequencies'
  ! CH4CH4 collisions
  allocate(ch4ch4_g(lt,knu))
  open(unit=12, file=trim(path_input)//FILE_CH4CH4, status='old', form='unformatted')
  read(12)lt, knu, (tmod(i),i=1,lt), (gnu(i),i=1,knu), ((ch4ch4_g(i,j),i=1,lt),j=1,knu)
  close(12)
  write(*,fmt='(a,a)')'CH4-CH4 spectroscopic file :', trim(file_ch4ch4)
  write(*,fmt='(5x,a, f8.2,a,f8.2,a,i3,a,i4,a)')'from', gnu(1), ' to', gnu(knu),' cm-1: ', lt, ' temperatures', knu, &
                                                ' frequencies'
  !~ N2H2 collisions
  allocate(n2h2_g(lt,knu))
  open(unit=12, file=trim(path_input)//FILE_N2H2, status='old', form='unformatted')
  read(12)lt, knu, (tmod(i),i=1,lt), (gnu(i),i=1,knu), ((n2h2_g(i,j),i=1,lt),j=1,knu)
  close(12)
  write(*,fmt='(a,a)')'N2-H2 spectroscopic file :', trim(file_n2h2)
  write(*,fmt='(5x,a, f8.2,a,f8.2,a,i3,a,i4,a)')'from', gnu(1), ' to', gnu(knu),' cm-1: ', lt, ' temperatures', knu, &
                                                ' frequencies'
  df_cont = gnu(2 )- gnu(1)
  n1 = max0( idint( (F_start-gnu(1)) / df_cont ) + 1, 1)    ! F_start dans le spectre qu'on calcule
  n2 = min0( idint( (F_end-gnu(1)) / df_cont ) + 2, knu)    ! pareil pour F_end
  if(n1 >= knu) then
    n1 = knu-1
  end if
  if(n2 <= 1) then
    n2 = 2
  end if
  nfcont = n2 - n1 + 1
  f1_cont = gnu(n1)
  write(*,fmt='(a,f8.2,a,f8.2,a)')'Data from', f1_cont, ' to', gnu(n2), ' cm-1 are used'
  write(*,*)"Done."
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.11: DAT_CONT                                                                                       ---!
  !--------------------------------------------------------------------------------------------------------------------!
  allocate(tab_n2n2(nlay,nfcont), tab_n2ch4(nlay,nfcont), tab_ch4ch4(nlay,nfcont), tab_n2h2(nlay,nfcont))
  !$OMP PARALLEL IF (n2-n1 > 20)
  call dat_cont(n2n2_g, n2ch4_g, ch4ch4_g, n2h2_g, tmod, n1, n2, tl, nlay, tab_n2n2, tab_n2ch4, tab_ch4ch4, tab_n2h2, &
                knu, nfcont, lt)
  !$OMP END PARALLEL
  deallocate(tmod, gnu, n2n2_g, n2ch4_g, ch4ch4_g, n2h2_g)
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 1.7: Extract file for calculation of the Voigt profile                                               ---!
  !--------------------------------------------------------------------------------------------------------------------!
  open(unit=10, file=trim(path_input)//'VOIGT.DAT', status='old', form='unformatted')
  read(10)v
  close(10)
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.14: Spectral extinction of cloud IC particles                                                      ---!
  !--------------------------------------------------------------------------------------------------------------------!
  allocate(qref(nb_cloud), nond(nb_cloud))
  if(nb_cloud > 0) then
    do ic = 1, nb_cloud
      open(unit=9, file=trim(path_input)//file_scat(ic), status='old', form='formatted')
      read(9,*)nond(ic)
      close(9)
    end do
    nond_max = maxval(nond)
    allocate(sig(nb_cloud,nond_max), qex(nb_cloud,nond_max))
    do ic = 1, nb_cloud
      open(unit=9, file=trim(path_input)//file_scat(ic), status='old', form='formatted')
      read(9,*)
      do i = 1, nond(ic)
        read(9,*)sig(ic,i), qex(ic,i)
      end do
      close(9)
      qref(ic) = tina(f_cloud(ic), sig(ic,:), qex(ic,:), nond(ic))
    end do
  end if
!----------------------------------------------------------------------------------------------------------------------!
!--- Section 1.3: Select absorbers for which K matrices are calculated                                              ---!
!----------------------------------------------------------------------------------------------------------------------!
  allocate(icorps_k(nb_mol_inv))
  do ik = 1, nb_mol_inv
    if(name_mol_inv(ik) == 'HAZE') then
      icorps_k(ik) = 0
      ikh = ik
    end if
    do i = 1, nb_mol
      if(name_mol_inv(ik) == name_mol(i,1)) then
        icorps_k(ik) = i
      end if
    end do
  end do

  ! Each file_spectro have different size, => need to find the longer one
  allocate(nbl_sp(nb_mol))
  do k = 1, nb_mol
    nbl_sp(k) = nb_line(trim(path_input)//file_spectro(k)) ! here file unit id = 400
    write(*,*)k, nbl_sp(k)
    if (nbl_sp(k) >= nbl_sp_max) then
      nbl_sp_max = nbl_sp(k)
    end if
  end do
  write(*,*)nbl_sp_max
  ! Extract data
  allocate(w_input(nb_mol,nbl_sp_max), s_input(nb_mol,nbl_sp_max), g_input(nb_mol,nbl_sp_max), &
           e_input(nb_mol,nbl_sp_max))
  do k = 1, nb_mol
    open(newunit=unit, file=trim(path_input)//file_spectro(k), status='old')
    do j = 1, nbl_sp(k)
      read(unit, fmt="(f10.6,e10.3,f5.3,f10.3)")w_input(k,j), s_input(k,j), g_input(k,j), e_input(k,j)
    end do
    close(unit)
  end do
  !--------------------------------------------------------------------------------------------------------------------!
  !--- Section 2.15: Check inversion parameters                                                                     ---!
  !--------------------------------------------------------------------------------------------------------------------!
  write(*,*)"--------------------------"
  write(*,*)"Check inversion parameters"
  nb_obs = abs(nb_obs)
  write(*,fmt='(TR5,A,I2,A,20(A,X))')'K matrices are calculated for ', nb_mol_inv, ' absorbers: ', &
                                     (trim(name_mol_inv(i)),i=1,nb_mol_inv)
  write(*,fmt='(TR5,A,F7.3,A,F7.3,A,I4,A)')'Inverse calculation between ', freq_x(1), ' and ', freq_x(nx), &
                                           ' cm-1 (', nx, ' points)'
  write(*,fmt='(TR5,A,I3)')'Numbers of airmasses: ', nb_obs
  write(*,fmt='(1x,a4,a,f7.2,a,5x,a,e10.3)')(name_mol_inv(i), ': correlations length:', corr_len(i), ' scale height,', &
                                            'weight =', wei(i), i=1, nb_mol_inv)
  write(*,fmt='(a,i2,a)')' iterations, abundance profile written on: ', niter, trim(file_profq)
  if(nb_kernel > 0) then
    write(*,fmt='(a,a,i6,a,8f10.3)')'Kernels are written on file: ', trim(file_kernel), nb_kernel, &
                                    ' frequencies: ', (F_KERNEL(i), i=1, NB_KERNEL)
  end if

  !====================================================================================================================!
  !=== PART 3: ITERATIVE INVERSE METHODE                                                                            ===!
  !====================================================================================================================!
  !--- Section 3.1: Allocate variables                                                                              ---!
  !--------------------------------------------------------------------------------------------------------------------!
  m = nb_obs * nx
  allocate(matrix_k(nb_mol_inv,nlay,nb_obs,nfreq_inv), rad(nb_obs,nfreq_inv), rad_diff(nb_obs,nx), matrix_r(nlay,m), &
           s(nlay,nlay,nb_mol_inv), sk(nb_mol_inv,nlay,m), ssk(nb_mol_inv,m,m), sm(nlay,m), s_b(nlay,nlay), &
           ssk_b(m,m), matrix_r2(nlay,m), errt(nlay), errtb(nlay), cf(nlay,nlay), cf_b(nlay,nlay), &
           dlnerrqb(nlay,nb_mol), errqb_sup(nlay,nb_mol), errqb_inf(nlay,nb_mol), dlnerrq(nlay,nb_mol), &
           errq_sup(nlay,nb_mol), errq_inf(nlay,nb_mol), matrix_t(nlay,nb_obs, nfreq_inv))
  rad_diff(:,:)     = 0d0
  rad(:,:)          = 0d0
  matrix_k(:,:,:,:) = 0d0
  matrix_r(:,:)     = 0d0
  matrix_t(:,:,:)   = 0d0

  do iter = 0, niter
    write(*,fmt='(/,30X,A,I2,A)')' ******* ITER = ', iter, ' *******'
    !==================================================================================================================!
    !=== PART 4: RADIATIVE TRANSFER                                                                                 ===!
    !==================================================================================================================!
    call transfer_abundance(nlev, nlay, vmr_tropo_n2, vmr_tropo_h2, vmr_tropo_ar, pl, tl, ml, gl, qch4, ql, p, ap, t, &
                            nb_obs, sec, lay_min, nb_mol, mass, expo_rot, et, f_start, f_end, f_pas_mult, f_pas, &
                            file_weightfunc, file_spec_synth, convfunc, rayleigh, nb_weightfunc, w_input, s_input, &
                              g_input, e_input, nb_cloud, taucl, cross_section_min, temperature_test, nlor, alor, glor, &
                            elor, g2lor, nvib, vib, ndeg, nb_mol_inv, icorps_k, matrix_k, rad, ikh, nfcont, nond, &
                            nond_max, sig, qex, nbl_sp_max, nbl_sp, f1_cont, df_cont, iter, nfreq_inv, qref, tab_n2n2, &
                            tab_n2ch4, tab_ch4ch4, tab_n2h2, v, name_mol(:,1), path_input, limbe, mode_inversion, &
                            matrix_t)
    !------------------------------------------------------------------------------------------------------------------!
    !--- Section 5.7: calculate residu (observed - synthetic)                                                       ---!
    !------------------------------------------------------------------------------------------------------------------!
    write(*,fmt='(3x,a,5x,a,10x,a,8x,a,3x,a)')'freq_x', 'noise', 'obs', 'synth', 'obs - synth'
    do ix = 1, nx
        i = nint( (freq_x(ix) - F_start) / F_pas + 0.5 )
        do is = 1, nb_obs
            rad_diff(is,ix) = rad_obs(is,ix) - rad(is,i)
        end do
        write(*,fmt='(f10.3,e12.3,2x,3e12.4,i5)')freq_x(ix), rad_noise(ix), rad_obs(1,ix), rad(1,i), rad_diff(1,ix), 1
        if(nb_obs > 1) then
            do is = 2, nb_obs
                write(*,fmt='(24x,3e12.4,i5)')rad_obs(is,ix), rad(is,i), rad_diff(is,ix), is
            enddo
        endif
    end do
    !------------------------------------------------------------------------------------------------------------------!
    !--- Section 5.8: Calculate root mean square (RMS)                                                              ---!
    !------------------------------------------------------------------------------------------------------------------!
    write(*,*)'NESR**2 * nb_wavenumber * n_obs =', rad_noise(nx)**2 * m
    write(*,*)'*** rms ***'
    do  is = 1, nb_obs
        rms = 0d0
        do ix = 1, nx
            i = nint( (freq_x(ix) - F_start) / F_pas + 0.5d0)
            rms = rms + rad_diff(is,ix)**2
        end do
        write(*,fmt='(a,i2,a,e10.3)')' is = ', is, ': rms = ', dsqrt(rms/dfloat(nx))
    end do

    if (mode_inversion /= 0) then
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 4.1: print K-matrices if iprint > 1                                                              ---!
      !----------------------------------------------------------------------------------------------------------------!
      if(iprint > 1) then
          do ix = 1, nx
              i = nint( (freq_x(ix)-F_start) / F_pas_mult + 0.5)
              write(*,fmt='(a,f8.2,a)')' K-matrices at wavenumber =', freq_x(ix), ' cm-1:'
              do is = 1, nb_obs
                  write(*,fmt='(5x,a,i2)')'Secant #', is
                  do ik = 1, nb_mol_inv
                      write(*,fmt='(10x,a5)')name_mol_inv(ik)
                      if(icorps_k(ik) > 0) then
                          write(*,fmt='(12x,e14.4,f9.2,e15.4,5x,e15.4)')(pl(j), tl(j), ql(j,icorps_k(ik)), &
                                                                        matrix_k(ik,j,is,i), j = nlay, 1, -1)
                      else
                          write(*,fmt='(12x,e14.4,f9.2,e15.4,5x,e15.4)')(pl(j), tl(j), taucl_tot(j), &
                                                                        matrix_k(ik,j,is,i), j = nlay, 1, -1)
                      end if
                  end do
              end do
          end do
      end if
      !================================================================================================================!
      !=== PART 5: INVERSE METHODE                                                                                  ===!
      !================================================================================================================!
      !--- Section 5.1: Calculate S * Kt matrice                                                                    ---!
      !----------------------------------------------------------------------------------------------------------------!
      sk(:,:,:) = 0d0
      do ik = 1, nb_mol_inv
          do i = 1, nlay
              do k = 1, nlay
                  var_expo = (-0.5d0 * (dlog(pl(i)/pl(k))**2)) / (corr_len(ik)**2)
                  if (var_expo < -200) then
                      corr = 0d0
                  else
                      corr = wei(ik) * dexp(var_expo)
                  end if
                  s(i,k,ik) = corr
                  do j = 1, m
                      is = 1 + (j-1)/nx
                      ix = j - (is-1)*nx
                      ifreq = nint( (freq_x(ix) - F_start)/F_pas + 0.5 )
                      sk(ik,i,j) = sk(ik,i,j) + corr * matrix_k(ik,k,is,ifreq)
                  end do
              end do
          end do
      end do
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 5.2: Calculate K * S *Kt matrices                                                                ---!
      !----------------------------------------------------------------------------------------------------------------!
      ssk(:,:,:) = 0d0
      do ik = 1, nb_mol_inv
          do i = 1, m
              is = 1 + (i-1)/nx
              ix = i - (is-1)*nx
              ifreq = nint( (freq_x(ix) - F_start)/F_pas + 0.5 )
              do j = 1, m
                  do k = 1, nlay
                      ssk(ik,i,j) = ssk(ik,i,j) + matrix_k(ik,k,is,ifreq) * sk(ik,k,j)
                  end do
              end do
          end do
          tra = 0d0
          do j = 1, m
              tra = tra + ssk(ik,j,j)
          end do
          write(*,fmt='(a,i2,a,f10.5)')' Trace of matrix a K * S * Kt for absorber', ik, ':', tra
      end do

      ssk_b(:,:) = 0d0
      if (mode_inversion == 2) then
        !--------------------------------------------------------------------------------------------------------------!
        !--- Section 5.3: Calculate S * Mt matrice                                                                  ---!
        !--------------------------------------------------------------------------------------------------------------!
        corrT2 = corr_len_T**2
        sm(:,:) = 0d0
        do i = 1, nlay
            do k = 1, nlay
                plo = dlog(pl(i)/pl(k))
                var_expo=-0.5d0 * plo**2/corrT2
                if (var_expo < -200) then
                    corr_b = 0d0
                else
                    corr_b = wei_t * exp(var_expo)
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
        !--------------------------------------------------------------------------------------------------------------!
        !	Calculate M * S * Mt matrices
        !--------------------------------------------------------------------------------------------------------------!
        do i = 1, m
          is = 1 + (i-1)/nx
          ix = i - (is-1)*nx
          ifreq = nint( (freq_x(ix) - F_start) / F_pas + 0.5)
          do j = 1 , m
            do k = 1 , nlay
              ssk_b(i,j) = ssk_b(i,j) + matrix_t(k,is,ifreq) * sm(k,j)
            end do
          end do
        end do
        tra_b = 0d0
        do j = 1 , m
            tra_b = tra_b + ssk_b(j,j)
        end do
        write(*,fmt='(a,f10.5)')'Trace of matrix a M * S * Mt (for temperature):', tra_b
      end if
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 5.3: Calculate Matrix C                                                                          ---!
      !----------------------------------------------------------------------------------------------------------------!
      allocate(c(m,m), c0(m))
      do i = 1, nx
        do is = 1, nb_obs
          j = i + (is-1) * nx
          c(j,j) = rad_noise(i)**2
          c0(j) = c(j,j)
        end do
      end do
      do i = 1, m
        do j = 1, m
          if(i /= j) then
            c(i,j) = 0d0
          end if
          do ik = 1, nb_mol_inv
            c(i,j) = c(i,j) + ssk(ik,i,j) + ssk_b(i,j)
          end do
        end do
      end do
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 5.4: Invert Matrix C                                                                             ---!
      !----------------------------------------------------------------------------------------------------------------!
      if (method_invert) then
          call choldc(c, m)
      else
          call matrix_inv(c, m, m)
      end if
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 5.5: Calculate final matrice and variations of temperatures                                      ---!
      !----------------------------------------------------------------------------------------------------------------!
      if (mode_inversion == 2) then
        matrix_r2(:,:) = 0d0
        do i = 1, nlay
            do j = 1, m
                do k = 1, m
                    matrix_r2(i,j) = matrix_r2(i,j) + sm(i,k) * c(k,j)
                end do
            end do
        end do

        do i = nlay, 1, -1
            dt = 0d0
            do j = 1, m
                is = 1 + (j-1)/nx
                ifreq = j - (is-1)*nx
                dt = dt + matrix_r2(i,j) * rad_diff(is,ifreq)
            end do
            tl(i) = tl(i) + dt
            write(*,fmt='(i3,e14.4,f10.3,5x,f10.3)')i,pl(i),dt,tl(i)
        end do
        t(1) = t(1) + dt

        ! Calculate covariance matrix
        do i = 1, nlay
          errtb(i) = 0d0
          do j = 1, m
            errtb(i) = errtb(i) + c0(j) * matrix_r2(i,j)**2
          end do
          do k = 1, nlay
            cf_b(i,k) = s_b(i,k)
            do  j = 1, m
              cf_b(i,k) = cf_b(i,k) - matrix_r2(i,j) * sm(k,j)
            end do
          end do
        end do

        do i = 1, nlay
            errt(i) = dsqrt(cf_b(i,i))
            errtb(i) = dsqrt(errtb(i))
        end do

        ! Writing T vertical profile
        open(unit=10, file=file_prof_t, status='unknown')
        do i = 1, ninput
            if(pii(i) < p(1)) then
                write(10,fmt='(e12.4,f10.3,2f12.2)')pii(i), ti(i)
            end if
        end do
        write(10,fmt='(e12.4,f10.3,2f12.2)')p(1), t(1)
        write(10,fmt='(e12.4,f10.3,2f12.2)')(pl(i), tl(i), errt(i), errtb(i),i=1,nlay)
        write(*,*)('errt = ', errt(i), i = 1, nlay)
        close(10)
      end if
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 5.6: Calculate final matrices and variations of absorbers                                        ---!
      !----------------------------------------------------------------------------------------------------------------!
      matrix_r(:,:) = 0d0
      do ik = 1, nb_mol_inv
        write(*,fmt='(10x,a5)')name_mol_inv(ik)
        do i = 1, nlay
          do j = 1, m
            do k = 1, m
              matrix_r(i,j) = matrix_r(i,j) + sk(ik,i,k) * c(k,j)
            end do
          end do
        end do

        do i = nlay, 1, -1
          dlnq=0d0
          do j = 1, m
            is = 1 + (j-1)/nx
            ifreq = j - (is-1)*nx
            dlnq = dlnq + matrix_r(i,j) * rad_diff(is,ifreq)
          end do
          if(icorps_k(ik) > 0) then
            ql(i,icorps_k(ik)) = ql(i,icorps_k(ik)) * dexp(dlnq)
            ! Remove any supersaturation
            qsat = 1.01325d+03 * 10**( a(icorps_k(ik)) - b(icorps_k(ik))/tl(i) ) / pl(i)
            ql(i,icorps_k(ik)) = dmin1(qsat, ql(i,icorps_k(ik)))
            write(*,fmt='(i3,e14.4,f12.4,e15.4)')i, pl(i), dlnq, ql(i,icorps_k(ik))
          else
            taucl_tot(i) = taucl_tot(i) * dexp(dlnq)
            write(*,fmt='(i3,e14.4,f12.4,e15.4)')i, pl(i), dlnq, taucl_tot(i)
            if(nb_cloud > 0) then
              do ic = 1, nb_cloud
                taucl(ic,i) = taucl(ic,i) * dexp(dlnq)
              end do
            end if
          end if
        end do

        if(name_mol_inv(ik) == 'HAZE') then
          if(nb_cloud > 0) then
            do ic = 1, nb_cloud
              write(*,*)'nb_cloud= ', ic
              do i = nlay, 1, -1
                write(*,fmt='(i3,e14.4,f12.4,e15.4)')i, pl(i), taucl(ic,i)
              end do
            end do
          end if
          do j = 1, nlay
            ql(j,icorps_k(ik)) = taucl_tot(j)
          end do
        end if
      end do
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 5.6: Errorbar calculation by Barney's Method                                                     ---!
      !----------------------------------------------------------------------------------------------------------------!
      do ik = 1, nb_mol_inv
        do i=1,nlay
          dlnerrqb(i,ik) = 0d0
          do j = 1, m
            dlnerrqb(i,ik) = dlnerrqb(i,ik) + c0(j) * matrix_r(i,j)**2
          end do
        end do

        do i = 1, nlay
          dlnerrqb(i,ik) = dsqrt(dlnerrqb(i,ik))
          errqb_sup(i,ik) = ql(i,icorps_k(ik)) * dexp(dlnerrqb(i,ik))
          errqb_inf(i,ik) = ql(i,icorps_k(ik)) / dexp(dlnerrqb(i,ik))
        end do
        !--------------------------------------------------------------------------------------------------------------!
        !--- Section 5.7: Errorbar calcultion by a second method                                                    ---!
        !--------------------------------------------------------------------------------------------------------------!
        do i = 1, nlay
          do k= 1, nlay
            cf(i,k) = s(i,k,ik)
            do j = 1, m
              cf(i,k) = cf(i,k) - matrix_r(i,j) * sk(ik,k,j)
            end do
          end do
        end do
        do i = 1, nlay
          dlnerrq(i,ik) = dsqrt(cf(i,i))
          errq_sup(i,ik) = ql(i,icorps_k(ik)) * dexp(dlnerrq(i,ik))
          errq_inf(i,ik) = ql(i,icorps_k(ik)) / dexp(dlnerrq(i,ik))
        end do
      end do
      !================================================================================================================!
      !=== PART 6: WRINTING OUTPUT                                                                                  ===!
      !================================================================================================================!
      !--- Section 6.1: Writing abundance vertical profiles                                                         ---!
      !----------------------------------------------------------------------------------------------------------------!
      open(unit=10, file=file_profq, status='unknown')
      write(10,fmt='(28x,5(4x,a4,4x,48x))')(name_mol_inv(ik), ik=1, nb_mol_inv)
      do j = nlay, 1, -1
        write(10,fmt='(f8.2,e11.4,f9.2,5(5e12.4))')0.5*(z(j)+z(j+1)), pl(j), tl(j), (ql(j,icorps_k(ik)), &
                                                   errq_inf(j,ik), errq_sup(j,ik), errqb_inf(j,ik), errqb_sup(j,ik), &
                                                   ik=1, nb_mol_inv)
      end do
      close(10)
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 6.2: Calculate kernels for molecules                                                             ---!
      !----------------------------------------------------------------------------------------------------------------!
      open(unit=11, file=file_kernel, status='unknown')
      do ifk = 1, nb_kernel
        i = nint( (f_kernel(ifk) - f_start)/f_pas + 0.5 )
        write(11,fmt='(20x,f10.3,a)')f_start+f_pas*(dfloat(i)-0.5d0), ' cm-1'
        do ik = 1, nb_mol_inv
          write(11,fmt='(28x,5(4x,a4,4x,48x))')name_mol_inv(ik)
          do j = nlay, 1, -1
            write(11,fmt='(f8.2,e13.5,f8.2,2x,*(e13.5))')0.5*(z(j)+z(j+1)), pl(j), tl(j), &
                                                         (matrix_k(ik,j,is,i)/(ap(j)-ap(j+1))/rad(is,i), is=1, nb_obs)
          end do
        end do
      end do
      close(11)
      !----------------------------------------------------------------------------------------------------------------!
      !--- Section 6.3: Calculate kernels for temperature                                                           ---!
      !----------------------------------------------------------------------------------------------------------------!
      if (mode_inversion == 2) then
        open(unit=12,file=file_kernel_t,status='unknown')
        do ik = 1, nb_kernel_t
            i = nint( (f_kernel_t(ik) - F_start) / F_pas + 0.5)
            write(12,fmt='(20x,f10.3,a)')F_start + F_pas * (dfloat(i)-0.5d0), ' cm-1'
            do j = nlay, 1, -1
                write(12,fmt='(f8.2,e13.5,f8.2,2x,*(e13.5))')0.5*(z(j)+z(j+1)), pl(j), tl(j), &
                                                             (matrix_t(j,is,i)/(ap(j)-ap(j+1))/rad(is,i), is=1, nb_obs)
            end do
        end do
        close(12)
      end if
      deallocate(c, c0)
    end if
  end do

  deallocate(z, p, ap, t, pl, apl, tl, gl, qch4, tab_i, lay_min, sec, ql, pii, ti, ati, pcli, taucli, apcli, ataucli, &
             taucl, taucl_tot, nlor, alor, glor, elor, g2lor, nvib, vib, ndeg, cross_section_min, temperature_test, &
             freq_x, rad_obs, rad_diff, rad_noise, sk, sm, ssk, s, s_b, ssk_b, cf, cf_b, dlnerrqb, errqb_sup, &
             errqb_inf, dlnerrq, errq_sup, errq_inf, errt, errtb, ml, matrix_k, matrix_r, matrix_r2, rad, tab_n2n2, &
             tab_n2ch4, tab_ch4ch4, tab_n2h2, qref, n_x, matrix_t)
!======================================================================================================================!
!=== END PROGRAM                                                                                                    ===!
!======================================================================================================================!
call system_clock(count=time_end, count_rate=time_precision)
call cpu_time(cpu_time_end)
time_elapsed = real(time_end - time_begin,kind=8)/real(time_precision,kind=8)
cpu_time_elapsed = cpu_time_end - cpu_time_begin
write(*,*)'=================================================='
write(*,*)'Program takes :', time_elapsed, ' s time.'
write(*,*)'Program takes :', cpu_time_elapsed, ' s cpu time.'

end program ABUNDANCE
