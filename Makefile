#On défini le compilateur à utiliser :
F90=/usr/bin/gfortran
#Les options à utiliser lors de la création des fichiers objets :
CFLAGS=-ffixed-line-length-none -mcmodel=large -Wall -free -fbounds-check -g -fbacktrace -ffpe-trap=all -ffpe-summary=all -pg -fcheck=all -fPIC

#!all : Compile tous les exécutables.
all : abundance transac calcul_direct abundance_debug transac_debug

debug : abundance_debug transac_debug calcul_direct_debug

transac : transac_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o \
	database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o \
	tau_lines2.o tau_lines3.o transfer_transac_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o
	$(F90) -o ../transac_v4 transac_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o \
	dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o \
	read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_transac_v4.o conv_k.o avg_k.o choldc.o \
	matrix_inv.o

abundance : abundances_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o \
	database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o \
	tau_lines2.o tau_lines3.o transfer_abundance_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o
	$(F90) -o ../abundances_v4 abundances_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o \
	dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o \
	read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_abundance_v4.o conv_k.o avg_k.o choldc.o \
	matrix_inv.o

abundance_para: abundances_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o \
	det_step.o dtau_e.o dtau_h.o \
	dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o \
	tau_lines1.o tau_lines2.o tau_lines3.o transfer_abundance_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o
	$(F90) -o ../abundances_para abundances_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o \
	det_step.o dtau_e.o dtau_h.o \
	dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o \
	tau_lines1.o tau_lines2.o tau_lines3.o transfer_abundance_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o

calcul_direct : calcul_direct_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_calcul_direct_v4.o conv_k.o avg_k.o
	$(F90) -o ../calcul_direct_v4 calcul_direct_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_calcul_direct_v4.o conv_k.o avg_k.o

transac_debug : transac_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_transac_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o
	$(F90) $(CFLAGS) -o ../transac_v4_debug transac_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_transac_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o

abundance_debug : abundances_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_abundance_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o
	$(F90) $(CFLAGS) -o ../abundances_v4_debug abundances_v4.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_abundance_v4.o conv_k.o avg_k.o choldc.o matrix_inv.o

calcul_direct_debug : calcul_direct.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_calcul_direct_v4.o conv_k.o avg_k.o
	$(F90) $(CFLAGS) -o ../calcul_direct calcul_direct.o lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o planck.o dplanck.o read_line2.o read_parameters.o tau_lines1.o tau_lines2.o tau_lines3.o transfer_calcul_direct_v4.o conv_k.o avg_k.o


abundances_mono.o: abundances_v4.f90
	$(F90) -c abundances_v4.f90

abundances_para.o: abundances_v4.f90
	$(F90) -fopenmp -c abundances_v4.f90

avg.o: avg.f90
	$(F90) -c avg.f90

calcul_direct_v4.o: calcul_direct_v4.f90
	$(F90) -c calcul_direct_v4.f90

conv.o: conv.f90
	$(F90) -c conv.f90

conv_wf.o: conv_wf.f90
	$(F90) -c conv_wf.f90

dat_cont.o: dat_cont.f90
	$(F90) -c dat_cont.f90 lib_functions.f90

det_step.o: det_step.f90
	$(F90) -c det_step.f90

dtau_e.o: dtau_e.f90
	$(F90) -c dtau_e.f90

dtau_h.o: dtau_h.f90
	$(F90) -c dtau_h.f90

dtau_v.o: dtau_v.f90
	$(F90) -c dtau_v.f90

dplanck.o: dplanck.f90
	$(F90) -c dplanck.f90

database_mol.o: database_mol.f90
	$(F90) -c database_mol.f90

database_planet.o: database_planet.f90
	$(F90) -c database_planet.f90

intwdl.o: intwdl.f90
	$(F90) -c intwdl.f90

lib_functions.o: lib_functions.f90
	$(F90) -c lib_functions.f90

planck.o: planck.f90
	$(F90) -c planck.f90 lib_functions.f90

read_line2.o: read_line2.f90
	$(F90) -c read_line2.f90

read_parameters.o: read_parameters.f90
	$(F90) -c read_parameters.f90

tau_lines1.o: tau_lines1.f90
	$(F90) -c tau_lines1.f90

tau_lines2.o: tau_lines2.f90
	$(F90) -c tau_lines2.f90

tau_lines3.o: tau_lines3.f90
	$(F90) -c tau_lines3.f90

transac_v4.o: transac_v4.f90
	$(F90) -c transac_v4.f90

transfer_transac_v4.o: transfer_transac_v4.f90
	$(F90) -c transfer_transac_v4.f90

transfer_abundance_v4.o: transfer_abundance_v4.f90
	$(F90) -c transfer_abundance_v4.f90

transfer_calcul_direct_v4.o: transfer_calcul_direct_v4.f90
	$(F90) -c transfer_calcul_direct_v4.f90

conv_k.o: conv_k.f90
	$(F90) -c conv_k.f90

avg_k.o: avg_k.f90
	$(F90) -c avg_k.f90

choldc.o: choldc.f90
	$(F90) -c choldc.f90

matrix_inv.o: matrix_inv.f90
	$(F90) -c matrix_inv.f90

#!help : Affiche cette aide.
#!
help:
	@grep -E "^#!" Makefile | sed -e 's:#!::g'

#!clean : efface les fichiers temporaires *.o
#!
clean:
	rm -f *.o
