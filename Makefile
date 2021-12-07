#On défini le compilateur à utiliser :
F90 = /usr/bin/gfortran

#Les options à utiliser lors de la création des fichiers objets :
DEBUG = -ffixed-line-length-none -mcmodel=large -Wall -free -fbounds-check -g -fbacktrace  \
	    -pg -fcheck=all -fPIC
#-ffpe-trap=all # not working with openMP

PROD = -O3

FLAG =$(PROD)
libf = lib_fortran

lib = lib_functions.o avg.o conv.o conv_wf.o dat_cont.o det_step.o dtau_e.o dtau_h.o dtau_v.o database_mol.o database_planet.o intwdl.o \
	planck.o dplanck.o read_line2.o read_parameters.o conv_k.o avg_k.o choldc.o matrix_inv.o tau_lines1.o tau_lines2.o \
	tau_lines3.o

lib_o = tmp/lib_functions.o tmp/avg.o tmp/conv.o tmp/conv_wf.o tmp/dat_cont.o tmp/det_step.o \
		tmp/dtau_e.o tmp/dtau_h.o tmp/dtau_v.o tmp/database_mol.o tmp/database_planet.o tmp/intwdl.o tmp/planck.o \
		tmp/dplanck.o tmp/read_line2.o tmp/read_parameters.o tmp/conv_k.o tmp/avg_k.o tmp/choldc.o tmp/matrix_inv.o \
		tmp/tau_lines1.o tmp/tau_lines2.o tmp/tau_lines3.o
#----------------------------------------------------------------------------------------------------------------------#
#---- All: compile all programs
#----------------------------------------------------------------------------------------------------------------------#
all : abundance_para
#----------------------------------------------------------------------------------------------------------------------#
#---- Programs
#----------------------------------------------------------------------------------------------------------------------#
abundance_para: $(lib) radiative_transfer.o
	$(F90) $(FLAG) -c -fopenmp  $(libf)/lib_functions.f90 $(libf)/abundances_v4.f90
	mv abundances_v4.o  lib_functions.o lib_functions.mod tmp/
	$(F90) $(FLAG) -fopenmp -o bin/abundances_para tmp/abundances_v4.o tmp/radiative_transfer.o $(lib_o)
#----------------------------------------------------------------------------------------------------------------------#
#---- File
#----------------------------------------------------------------------------------------------------------------------#
avg.o: $(libf)/avg.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/avg.f90
	mv avg.o tmp/

conv.o: $(libf)/conv.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/lib_functions.f90 $(libf)/conv.f90
	mv conv.o lib_functions.o lib_functions.mod  tmp/

conv_wf.o: $(libf)/conv_wf.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/lib_functions.f90 $(libf)/conv_wf.f90
	mv conv_wf.o lib_functions.o lib_functions.mod tmp/

dat_cont.o: $(libf)/dat_cont.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/lib_functions.f90 $(libf)/dat_cont.f90
	mv dat_cont.o lib_functions.mod tmp/

det_step.o: $(libf)/det_step.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/det_step.f90
	mv det_step.o tmp/

dtau_e.o: $(libf)/dtau_e.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/dtau_e.f90
	mv dtau_e.o tmp/

dtau_h.o: $(libf)/dtau_h.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/dtau_h.f90
	mv dtau_h.o tmp/

dtau_v.o: $(libf)/dtau_v.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/dtau_v.f90
	mv dtau_v.o tmp/

dplanck.o: $(libf)/dplanck.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/lib_functions.f90 $(libf)/dplanck.f90
	mv dplanck.o  lib_functions.mod tmp/

database_mol.o: $(libf)/database_mol.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/database_mol.f90
	mv database_mol.o tmp/

database_planet.o: $(libf)/database_planet.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/database_planet.f90
	mv database_planet.o tmp/

intwdl.o: $(libf)/intwdl.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/intwdl.f90
	mv intwdl.o tmp/

lib_functions.mod: lib_functions.o $(libf)/lib_functions.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/lib_functions.f90
	mv lib_functions.mod tmp/

lib_functions.o: $(libf)/lib_functions.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/lib_functions.f90
	mv lib_functions.o tmp/

planck.o: $(libf)/planck.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/lib_functions.f90 $(libf)/planck.f90
	mv planck.o lib_functions.o tmp/

read_line2.o: $(libf)/read_line2.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/read_line2.f90
	mv read_line2.o tmp/

read_parameters.o: $(libf)/read_parameters.f90
	$(F90) $(FLAG) -c  -fopenmp $(libf)/read_parameters.f90
	mv read_parameters.o tmp/

radiative_transfer.o: $(libf)/radiative_transfer.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/lib_functions.f90 $(libf)/radiative_transfer.f90
	mv radiative_transfer.o lib_functions.o lib_functions.mod tmp/

tau_lines1.o: $(libf)/tau_lines1.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/tau_lines1.f90
	mv tau_lines1.o tmp/

tau_lines2.o: $(libf)/tau_lines2.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/tau_lines2.f90
	mv tau_lines2.o tmp/

tau_lines3.o: $(libf)/tau_lines3.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/tau_lines3.f90
	mv tau_lines3.o tmp/

conv_k.o: $(libf)/conv_k.f90
	$(F90) $(FLAG) -c -fopenmp $(libf)/lib_functions.f90 $(libf)/conv_k.f90
	mv conv_k.o lib_functions.mod tmp/

avg_k.o: $(libf)/avg_k.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/avg_k.f90
	mv avg_k.o tmp/

choldc.o: $(libf)/choldc.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/choldc.f90
	mv choldc.o tmp/

matrix_inv.o: $(libf)/matrix_inv.f90
	$(F90) $(FLAG) -c -fopenmp  $(libf)/matrix_inv.f90
	mv matrix_inv.o tmp/

#!help : Affiche cette aide.
#!
help:
	@grep -E "^#!" Makefile | sed -e 's:#!::g'

#!clean : efface les fichiers temporaires *.o
clean:
	rm -f tmp/*.o tmp/*.mod
