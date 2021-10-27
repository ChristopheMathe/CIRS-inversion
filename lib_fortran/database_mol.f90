!**	HCN based on Atreya (near 128 K) ; HC3N from P. Romani (157-251 K)
!**     modif du 21/08/2014 : HCN based on Fray and Schmitt 2009 (120-140 K)
!**     modif du 25/06/2015 : HC3N based on Fray and Schmitt 2009 (120-140 K)
!**     modif du 02/06/2016 : C2H2 based on Fray and Schmitt 2009
!**     modif du 02/06/2016 : C2H6 based on Fray and Schmitt 2009
!**     modif du 02/06/2016 : C2H4 based on Fray and Schmitt 2009

SUBROUTINE database_mol(input_mol_nb,input_mol,focal_plane,nlor,nvib,nb_deg,a,b,mass,&
expo_rot,cross_section_min,temperature_test,alor,vib,glor,elor,g2lor,&
file_spectro)

implicit none
! LOCAL
integer, parameter :: nvibmax = 5 ! maximal number of vibration for a molecule
integer :: i = 0

! INTENT (IN)
integer, intent(in) :: input_mol_nb
character(*), dimension(input_mol_nb), intent(in) :: input_mol
logical, intent(in) :: focal_plane ! if true == fp3, false == fp4

! INTENT(OUT)
integer, dimension(input_mol_nb), intent(out) :: nlor
integer, dimension(input_mol_nb), intent(out) :: nvib
integer, dimension(input_mol_nb,nvibmax), intent(out) :: nb_deg

double precision, dimension(input_mol_nb), intent(out) :: a
double precision, dimension(input_mol_nb), intent(out) :: b
double precision, dimension(input_mol_nb), intent(out) :: mass
double precision, dimension(input_mol_nb), intent(out) :: expo_rot
double precision, dimension(input_mol_nb), intent(out) :: cross_section_min
double precision, dimension(input_mol_nb), intent(out) :: temperature_test
double precision, dimension(input_mol_nb), intent(out) :: alor
double precision, dimension(input_mol_nb,nvibmax), intent(out) :: vib
double precision, dimension(input_mol_nb,nvibmax), intent(out) :: glor
double precision, dimension(input_mol_nb,nvibmax), intent(out) :: elor
double precision, dimension(input_mol_nb,nvibmax), intent(out) :: g2lor

character(*), dimension(input_mol_nb), intent(out) :: file_spectro

!======================================================================!
!======================================================================!
do i = 1, input_mol_nb
    select case(input_mol(i))
        case('3CH3')
            mass(i) = 18.04d0
            a(i) = 4.585d0
            b(i) = 502.1d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = '13CH3D_Conor_corrige.dat'
            cross_section_min(i) = 3.000E-29 
            temperature_test(i) = 140.0d0
            nlor(i)   = 5
            alor(i)   = 35.0d0
            glor(i,:) = (/0.0450d0, 0.0500d0, 0.0550d0, 0.0600d0, 0.0650d0/)
            elor(i,:) = (/0.7500d0, 0.7500d0, 0.7500d0, 0.7500d0, 0.7500d0/)
            g2lor(i,:)= (/0.0475d0, 0.0525d0, 0.0575d0, 0.0625d0, 0.0000d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('3CH4')
            mass(i) = 17.04d0
            a(i) = 4.585d0
            b(i) = 502.1d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = '13CH4_Conor_corrige.dat'
            cross_section_min(i) = 2.000E-29
            temperature_test(i) = 160.0d0
            nlor(i)   = 5
            alor(i)   = 35.d0
            glor(i,:) = (/0.0450d0, 0.0500d0, 0.0550d0, 0.0600d0, 0.0650d0/)
            elor(i,:) = (/0.6400d0, 0.6500d0, 0.6600d0, 0.7100d0, 0.7500d0/)
            g2lor(i,:)= (/0.0475d0, 0.0525d0, 0.0575d0, 0.0625d0, 0.0000d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('CH3')
            mass(i) = 15.00d0
            a(i) = 7.00d0
            b(i) = 2000.00d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = ''
            cross_section_min(i) = 0.0d0
            temperature_test(i) = 0.0d0
            nlor(i)   = 0
            alor(i)   = 0.0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('CH3D')
            mass(i) = 17.04d0
            a(i) = 4.585d0
            b(i) = 502.1d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = '12CH3D_Conor_corrige.dat'
            cross_section_min(i) = 3.000E-29
            temperature_test(i) = 140.0d0
            nlor(i)   = 5
            alor(i)   = 35.0d0
            glor(i,:) = (/0.0450d0, 0.0500d0, 0.0550d0, 0.0600d0, 0.0650d0/)
            elor(i,:) = (/0.7500d0, 0.7500d0, 0.7500d0, 0.7500d0, 0.7500d0/)
            g2lor(i,:)= (/0.0475d0, 0.0525d0, 0.0575d0, 0.0625d0, 0.0000d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('CH4')
            mass(i) = 16.04d0
            a(i) = 4.585d0
            b(i) = 502.1d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = '12CH4_Conor_1000_1500.dat'
            cross_section_min(i) = 2.000E-29
            temperature_test(i) = 160.0d0
            nlor(i)   = 5
            alor(i)   = 35.0d0
            glor(i,:) = (/0.0450d0, 0.0500d0, 0.0550d0, 0.0600d0, 0.0650d0/)
            elor(i,:) = (/0.6400d0, 0.6500d0, 0.6600d0, 0.7100d0, 0.7500d0/)
            g2lor(i,:)= (/0.0475d0, 0.0525d0, 0.0575d0, 0.0625d0, 0.0000d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('CO')
            mass(i) = 28.01d0
            a(i) = 5.426d0
            b(i) = 424.9d0
            expo_rot(i) = 1.d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('CO2')
            mass(i) = 44.01d0
            a(i) = 7.501d0
            b(i) = 1441.5d0
            expo_rot(i) = 1.d0
            file_spectro(i) = 'co2_15mu.dat'
            cross_section_min(i) = 2.000E-28
            temperature_test(i) =160.d0
            nlor(i)   = 1
            alor(i)   = 5.d0
            glor(i,:) = (/0.080d0, 0d0, 0d0, 0d0, 0d0/)
            elor(i,:) = (/0.750d0, 0d0, 0d0, 0d0, 0d0/)
            g2lor(i,:)= (/0.000d0, 0d0, 0d0, 0d0, 0d0/)
            nvib(i)   = 1
            vib(i,:)  = (/667.d0, 0d0, 0d0, 0d0, 0d0/)
            nb_deg(i,:) = (/2     , 0  , 0  , 0  , 0  /)

        case('C2HD')
            mass(i) = 27.04d0
            a(i) = 6.979d0
            b(i) = 1284.3d0
            expo_rot(i) = 1.d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('C2H2')
            mass(i) = 26.04d0
            a(i) = 5.8138d0
            b(i) = 1101.4d0
            expo_rot(i) = 1.d0
            if (focal_plane) then
                  file_spectro(i) = 'fich_C2H2all.txt'
                  nlor(i)   = 5
                  alor(i)   = 10.d0
                  glor(i,:) = (/0.060d0, 0.070d0, 0.080d0, 0.090d0, 0.105d0/)
                  elor(i,:) = (/0.620d0, 0.710d0, 0.800d0, 0.800d0, 0.800d0/)
                  g2lor(i,:)= (/0.065d0, 0.075d0, 0.085d0, 0.097d0, 0.000d0/)
            else
                  nlor(i)   = 1
                  alor(i)   = 10.d0
                  glor(i,:) = (/0.100d0,0d0,0d0,0d0,0d0/)
                  elor(i,:) = (/0.75d0,0d0,0d0,0d0,0d0/)
                  g2lor(i,:)= (/0d0,0d0,0d0,0d0,0d0/)
                  file_spectro(i) = 'c2h2_7mu.dat'
            end if
            cross_section_min(i) = 5.000E-26
            temperature_test(i) = 160.0d0
            nvib(i)   = 2
            vib(i,:)  = (/612.0d0, 729.0d0, 0d0, 0d0, 0d0/)
            nb_deg(i,:) = (/2      , 2      , 0  , 0  , 0  /)

        case('C2H4')
            mass(i) = 28.05d0
            a(i) = 6.9540d0
            b(i) = 1026.63d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = 'c2h4_10mu.dat'
            cross_section_min(i) = 2.000E-25
            temperature_test(i) = 160.0d0
            nlor(i)   = 1
            alor(i)   = 35.d0
            glor(i,:) = (/0.087d0, 0d0, 0d0, 0d0, 0d0/)
            elor(i,:) = (/0.750d0, 0d0, 0d0, 0d0, 0d0/)
            g2lor(i,:)= (/0.000d0, 0d0, 0d0, 0d0, 0d0/)
            nvib(i)   = 3
            vib(i,:)  = (/826.d0, 943.d0, 949.d0, 0d0, 0d0/)
            nb_deg(i,:) = (/1     , 1     , 1     , 0  , 0  /)

        case('C2H6')
            mass(i) = 30.07d0
            a(i) = 7.0408d0
            b(i) = 1082.74d0
            expo_rot(i) = 1.5d0
            if (focal_plane) then
                  file_spectro(i) = 'c2h6_12mu.dat'
                  nlor(i)   = 1
                  alor(i)   = 10.0d0
                  glor(i,:) = (/0.090d0, 0.d0, 0.d0, 0.d0, 0.d0/)
                  elor(i,:) = (/0.900d0, 0.d0, 0.d0, 0.d0, 0.d0/)
                  g2lor(i,:)= (/0.000d0, 0.d0, 0.d0, 0.d0, 0.d0/)
            else
                  nlor(i)   = 2
                  alor(i)   = 10.0d0
                  glor(i,:) = (/0.080d0,0.090d0,0d0,0d0,0d0/)
                  elor(i,:) = (/0.9d0,0.9d0,0d0,0d0,0d0/)
                  g2lor(i,:)= (/0.085d0,0d0,0d0,0d0,0d0/)
                  file_spectro(i) = 'c2h6_7mu_linda.dat'
            end if
            cross_section_min(i) = 2.000E-26
            temperature_test(i) = 160.0d0
            nvib(i)   = 2
            vib(i,:)  = (/289.0d0, 821.0d0, 0.d0, 0.d0, 0.d0/)
            nb_deg(i,:) = (/1      , 2      , 0   , 0   , 0   /)

        case('C2N2')
            mass(i) = 52.035d0
            a(i) = 7.454d0
            b(i) = 1832.8d0
            expo_rot(i) = 1.d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('C3h4') !allene
            mass(i) = 40.06d0
            a(i) = 5.9025d0
            b(i) = 1339.9d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('C3H4') !methylacetylene
            mass(i) = 40.06d0
            a(i) = 6.029d0
            b(i) = 1444.7d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = 'c3h4_15mu.dat'
            cross_section_min(i) = 1.000d-23
            temperature_test(i) = 160.0d0
            nlor(i)   = 1
            alor(i)   = 35.d0
            glor(i,:) = (/0.100d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.75d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 2
            vib(i,:)  = (/331.0d0, 633.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/2    , 2    , 0    , 0    , 0    /)

        case('C3H8')
            mass(i) = 44.10d0
            a(i) = 5.191d0
            b(i) = 1164.3d0
            expo_rot(i) = 1.5d0
            if (focal_plane) then
!~                   file_spectro(i) = 'c3h8_13mu_nov08.dat'
                  file_spectro(i) = 'C3H8_lbl_Sung_FP3.txt'
            else
                  !file_spectro(i) = 'c3h8_7mu_nov08.dat'
                  file_spectro(i) = 'C3H8_Sung_Flaut_FP4.txt'
            end if
            cross_section_min(i) = 3.000E-25
            temperature_test(i) = 160.0d0
            nlor(i)   = 1
            alor(i)   = 5.d0
            glor(i,:) = (/0.12d0, 0d0, 0d0, 0d0, 0d0/)
            elor(i,:) = (/0.50d0, 0d0, 0d0, 0d0, 0d0/)
            g2lor(i,:)= (/0.00d0, 0d0, 0d0, 0d0, 0d0/)
            nvib(i)   = 5
            vib(i,:)  = (/216.d0, 268.d0, 369.d0, 748.d0, 869.d0/)
            nb_deg(i,:) = (/1     , 1     , 1     , 1     , 1     /)

        case('C4H2')
            mass(i) = 50.10d0
            a(i) = 7.070d0
            b(i) = 1897.1d0
            expo_rot(i) = 1.d0
            if (focal_plane) then
                file_spectro(i) = 'c4h2_16mu.datall'
            else
                file_spectro(i) = 'C4H2_nu6_nu8.txt'
            end if
            cross_section_min(i) = 1.000E-23
            temperature_test(i) = 160.0d0
            nlor(i)   = 1
            alor(i)   = 35.d0
            glor(i,:) = (/0.100d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.750d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.000d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 2
            vib(i,:)  = (/220.d0, 628.d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/2     , 2     , 0    , 0    , 0    /)

        case('C6H6')
            mass(i) = 78.11d0
            a(i) = 7.500d0
            b(i) = 2454.1d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = 'c6h6_15mu.dat'
            cross_section_min(i) = 3.000E-22
            temperature_test(i) = 160.0d0
            nlor(i)   = 1
            alor(i)   = 5.0d0
            glor(i,:) = (/0.110d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.750d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.000d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 5
            vib(i,:)  = (/398.0d0, 608.0d0, 674.0d0, 707.0d0, 846.0d0/)
            nb_deg(i,:) = (/2      , 2      , 1      , 1      , 2      /)

        case('HCN')
            mass(i) = 27.03d0
            a(i) = 7.466d0
            b(i) = 2061.2d0
            expo_rot(i) = 1.d0
            file_spectro(i) = 'HCN_G11_seul.txt'
!~             file_spectro(i) = 'hcn_seul_hitr.dat'
            cross_section_min(i) = 2.000E-25
            temperature_test(i) = 160.0d0
            nlor(i)   = 4
            alor(i)   = 35.d0
            glor(i,:) = (/0.10d0, 0.12d0, 0.14d0, 0.16d0, 0.d0/)
            elor(i,:) = (/0.80d0, 0.75d0, 0.80d0, 0.85d0, 0.d0/)
            g2lor(i,:)= (/0.11d0, 0.13d0, 0.15d0, 0.d0  , 0.d0/)
            nvib(i)   = 1
            vib(i,:)  = (/713.d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/2     , 0    , 0    , 0    , 0    /)

        case('hcn3')
            mass(i) = 28.0d0
            a(i) = 7.466d0
            b(i) = 2061.2d0
            expo_rot(i) = 1.d0
            file_spectro(i) = 'h13cn_G11_scaled.txt'
!~             file_spectro(i) = 'h13cn_hitr_scaled.dat'
            cross_section_min(i) = 2.000E-25
            temperature_test(i) = 160.0d0
            nlor(i)   = 4
            alor(i)   = 35.d0
            glor(i,:) = (/0.10d0, 0.12d0, 0.14d0, 0.16d0, 0.0d0/)
            elor(i,:) = (/0.80d0, 0.75d0, 0.80d0, 0.85d0, 0.0d0/)
            g2lor(i,:)= (/0.11d0, 0.13d0, 0.15d0, 0.00d0, 0.0d0/)
            nvib(i)   = 1
            vib(i,:)  = (/706.d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/2     , 0    , 0    , 0    , 0    /)
 
        case('hcn5')
            mass(i) = 28.0
            a(i) = 7.466
            b(i) = 2061.2
            expo_rot(i) = 1.
            file_spectro(i) = 'hc15n_G11_scaled.txt'
!~             file_spectro(i) = 'hc15n_hitr_scaled.dat'
            cross_section_min(i) = 2.000E-25
            temperature_test(i) = 160.0d0
            nlor(i)   = 4
            alor(i)   = 35.d0
            glor(i,:) = (/0.10d0, 0.12d0, 0.14d0, 0.16d0, 0.0d0/)
            elor(i,:) = (/0.80d0, 0.75d0, 0.80d0, 0.85d0, 0.0d0/)
            g2lor(i,:)= (/0.11d0, 0.13d0, 0.15d0, 0.00d0, 0.0d0/)
            nvib(i)   = 1
            vib(i,:)  = (/711.d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/2     , 0    , 0    , 0    , 0    /)
 
        case('HC3N')
            mass(i) = 51.00d0
            a(i) = 5.6441d0 ! valeur en atmo
            b(i) = 1922.1d0
            expo_rot(i) = 1.d0
            file_spectro(i) = 'hc3n_LISA.dat'
            cross_section_min(i) = 1.000E-23
            temperature_test(i) = 160.0d0
            nlor(i)   = 1
            alor(i)   = 35.d0
            glor(i,:) = (/0.10d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.75d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.00d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 3
            vib(i,:)  = (/223.0d0, 499.0d0, 663.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/2      , 2      , 2      , 0    , 0    /)

        case('HD')
            mass(i) = 3.01d0
            a(i) = 3.7447d0
            b(i) = 75.10d0
            expo_rot(i) = 0.89d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('H2')
            mass(i) = 2.016d0
            a(i) = 3.0441d0
            b(i) = 58.36d0
            expo_rot(i) = 0.95d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0 
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case('H2CO')
            mass(i) = 30.026d0
            a(i) = 9.5445d0
            b(i) = 2394.2d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)
 
        case('H20')
            mass(i) = 18.015d0
            a(i) = 7.5673d0
            b(i) = 2673.3d0
            expo_rot(i) = 1.5d0
            file_spectro(i) = ''
            cross_section_min(i) = 0d0
            temperature_test(i) = 0d0
            nlor(i)   = 0
            alor(i)   = 0d0
            glor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            elor(i,:) = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            g2lor(i,:)= (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nvib(i)   = 0
            vib(i,:)  = (/0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)
            nb_deg(i,:) = (/0    , 0    , 0    , 0    , 0    /)

        case default
            print *, "Unknown molecule ",input_mol(i)," in database_mol.f90"
            stop
        
    end select
end do
END SUBROUTINE database_mol
