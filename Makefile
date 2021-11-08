#On défini le compilateur à utiliser :
F90 = /usr/bin/gfortran

#Les options à utiliser lors de la création des fichiers objets :
DEBUG = -ffixed-line-length-none -mcmodel=large -Wall -free -fbounds-check -g -fbacktrace  \
	    -pg -fcheck=all -fPIC
#-ffpe-trap=all # not working with openMP

PROD = -O3

libf = lib_fortran

lib = lib_functions avg conv conv_wf dat_cont det_step dtau_e dtau_h dtau_v database_mol database_planet intwdl \
	planck dplanck read_line2 read_parameters conv_k avg_k choldc matrix_inv

lib_o = tmp/lib_functions.o tmp/avg.o tmp/conv.o tmp/conv_wf.o tmp/dat_cont.o tmp/det_step.o \
		tmp/dtau_e.o tmp/dtau_h.o tmp/dtau_v.o tmp/database_mol.o tmp/database_planet.o tmp/intwdl.o tmp/planck.o \
		tmp/dplanck.o tmp/read_line2.o tmp/read_parameters.o tmp/conv_k.o tmp/avg_k.o tmp/choldc.o tmp/matrix_inv.o
#----------------------------------------------------------------------------------------------------------------------#
#---- All: compile all programs
#----------------------------------------------------------------------------------------------------------------------#
all : abundance abundance_para transac calcul_direct abundance_debug transac_debug calcul_direct_debug
#----------------------------------------------------------------------------------------------------------------------#
#---- Programs
#----------------------------------------------------------------------------------------------------------------------#
abundance_para: radiative_transfer $(lib)
	$(F90) $(PROD) -c -fopenmp $(libf)/abundances_v4.f90
	mv abundances_v4.o tmp/
	$(F90) $(PROD) -fopenmp -o bin/abundances_para tmp/abundances_v4.o tmp/radiative_transfer.o $(lib_o)
#----------------------------------------------------------------------------------------------------------------------#
#---- File
#----------------------------------------------------------------------------------------------------------------------#
avg:
	$(F90) $(PROD) -c -fopenmp $(libf)/avg.f90
	mv avg.o tmp/

conv:
	$(F90) $(PROD) -c -fopenmp $(libf)/conv.f90
	mv conv.o tmp/

conv_wf:
	$(F90) $(PROD) -c  -fopenmp $(libf)/conv_wf.f90
	mv conv_wf.o tmp/

dat_cont:
	$(F90) $(PROD) -c  -fopenmp $(libf)/dat_cont.f90 $(libf)/lib_functions.f90
	mv dat_cont.o lib_functions.mod tmp/

det_step:
	$(F90) $(PROD) -c  -fopenmp $(libf)/det_step.f90
	mv det_step.o tmp/

dtau_e:
	$(F90) $(PROD) -c  -fopenmp $(libf)/dtau_e.f90
	mv dtau_e.o tmp/

dtau_h:
	$(F90) $(PROD) -c -fopenmp  $(libf)/dtau_h.f90
	mv dtau_h.o tmp/

dtau_v:
	$(F90) $(PROD) -c  -fopenmp $(libf)/dtau_v.f90
	mv dtau_v.o tmp/

dplanck:
	$(F90) $(PROD) -c  -fopenmp $(libf)/dplanck.f90
	mv dplanck.o tmp/

database_mol:
	$(F90) $(PROD) -c -fopenmp  $(libf)/database_mol.f90
	mv database_mol.o tmp/

database_planet:
	$(F90) $(PROD) -c -fopenmp  $(libf)/database_planet.f90
	mv database_planet.o tmp/

intwdl:
	$(F90) $(PROD) -c -fopenmp  $(libf)/intwdl.f90
	mv intwdl.o tmp/

lib_functions:
	$(F90) $(PROD) -c  -fopenmp $(libf)/lib_functions.f90
	mv lib_functions.o lib_functions.mod tmp/

planck:
	$(F90) $(PROD) -c  -fopenmp $(libf)/planck.f90 $(libf)/lib_functions.f90
	mv planck.o lib_functions.o lib_functions.mod tmp/

read_line2:
	$(F90) $(PROD) -c -fopenmp  $(libf)/read_line2.f90
	mv read_line2.o tmp/

read_parameters:
	$(F90) $(PROD) -c  -fopenmp $(libf)/read_parameters.f90
	mv read_parameters.o tmp/

radiative_transfer:
	$(F90) $(PROD) -c -fopenmp  $(libf)/radiative_transfer.f90
	mv radiative_transfer.o tmp/

conv_k:
	$(F90) $(PROD) -c -fopenmp  $(libf)/conv_k.f90
	mv conv_k.o tmp/

avg_k:
	$(F90) $(PROD) -c -fopenmp  $(libf)/avg_k.f90
	mv avg_k.o tmp/

choldc:
	$(F90) $(PROD) -c -fopenmp  $(libf)/choldc.f90
	mv choldc.o tmp/

matrix_inv:
	$(F90) $(PROD) -c -fopenmp  $(libf)/matrix_inv.f90
	mv matrix_inv.o tmp/

#!help : Affiche cette aide.
#!
help:
	@grep -E "^#!" Makefile | sed -e 's:#!::g'

#!clean : efface les fichiers temporaires *.o
clean:
	rm -f tmp/*.o
