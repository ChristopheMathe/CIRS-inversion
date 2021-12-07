subroutine read_parameters(filename_input,name_pla,lat,nlev,p_surf,&
p_top,prof_t,nb_cloud,f_cloud,prof_haze,file_scat,vmr_tropo_n2,vmr_tropo_h2,&
vmr_tropo_ar,vmr_tropo_ch4,file_n2n2,file_n2ch4,file_ch4ch4,file_n2h2,nb_mol,&
name_mol,nb_obs,los,file_obs,file_noise,f_start,f_end,f_pas,f_pas_mult,convfunc,rayleigh,&
nb_weightfunc,file_weightfunc,nb_df,df,nb_mol_inv,name_mol_inv,corr_len,&
wei,niter,iprint,nb_kernel,f_kernel,file_spec_synth,file_profq,file_kernel,limbe,&
nsp,nb_kernel_t,f_kernel_t,corr_len_t,wei_t,file_prof_t,file_kernel_t,mode_inversion)

implicit none
integer           :: nb_cloud,nb_mol,nb_mol_inv,nb_obs,nb_df,niter,iprint,nb_kernel,nb_kernel_t,nb_weightfunc,nlev,nsp
integer           :: mode_inversion
double precision  :: lat,p_surf,p_top,f_start,f_end,f_pas,f_pas_mult,rayleigh
double precision  :: vmr_tropo_n2,vmr_tropo_h2,vmr_tropo_ar,vmr_tropo_ch4
double precision   :: corr_len_t,wei_t
double precision, dimension(5) :: f_cloud
double precision, dimension(20)  :: los,f_kernel,f_kernel_t
double precision, dimension(20)   :: corr_len,wei
double precision, dimension(50,3) :: df

character(len=1)    :: convfunc
character(len=256)  :: filename_input
character(len=256)  :: name_pla
character(len=256)  :: prof_t,file_n2n2,file_n2ch4,file_ch4ch4,file_n2h2,file_weightfunc
character(len=256)  :: file_spec_synth,file_profq,file_kernel,file_prof_t,file_kernel_t
character(len=256)  :: file_noise
character(len=4)  ,dimension(20)    :: name_mol_inv
character(len=256),dimension(20)    :: file_obs
character(len=256),dimension(5)     :: prof_haze,file_scat
character(len=256),dimension(20,3)  :: name_mol

logical :: limbe

namelist /input_planet/       name_pla, lat, nlev, p_surf, p_top, prof_t, nb_cloud, f_cloud, prof_haze, file_scat,&
                              vmr_tropo_n2,vmr_tropo_h2,vmr_tropo_ar,vmr_tropo_ch4
namelist /input_collision/    file_n2n2, file_n2ch4, file_ch4ch4, file_n2h2
namelist /input_molecule/     nb_mol, name_mol
namelist /input_observations/ limbe, nb_obs, los, file_obs
namelist /input_convolution/  f_start, f_end, f_pas, f_pas_mult, convfunc, rayleigh, nb_weightfunc,file_weightfunc,file_noise,nsp
namelist /input_inversion/    mode_inversion, nb_df, df, nb_mol_inv, name_mol_inv, corr_len, wei, niter, iprint,&
                              nb_kernel, f_kernel,nb_kernel_t,f_kernel_t,wei_t,corr_len_t
namelist /output/             file_spec_synth, file_profq, file_kernel,file_prof_t,file_kernel_t

open(8, file=trim(filename_input), delim='apostrophe')
read(8, nml=input_planet)
read(8, nml=input_collision)
read(8, nml=input_molecule)
read(8, nml=input_observations)
read(8, nml=input_convolution)
read(8, nml=input_inversion)
read(8, nml=output)
close(8)
!write(*,nml=input_planet)
!write(*,nml=input_collision)
!write(*,nml=input_molecule)
!write(*,nml=input_observations)
!write(*,nml=input_convolution)
!write(*,nml=input_inversion)
!write(*,nml=output)
return
end subroutine read_parameters
