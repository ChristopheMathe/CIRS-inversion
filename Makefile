#On défini le compilateur à utiliser :
F90 = /usr/bin/gfortran
#Les options à utiliser lors de la création des fichiers objets :
CFLAGS = -ffixed-line-length-none -mcmodel=large -Wall -free -fbounds-check -g -fbacktrace -ffpe-trap=all \
	-pg -fcheck=all -fPIC

lib = lib_functions avg conv conv_wf dat_cont det_step dtau_e dtau_h dtau_v database_mol database_planet intwdl \
	planck dplanck read_line2 read_parameters tau_lines1 tau_lines2 tau_lines3 conv_k avg_k choldc matrix_inv

lib_o = tmp/lib_functions.o tmp/avg.o tmp/conv.o tmp/conv_wf.o tmp/dat_cont.o tmp/det_step.o \
		tmp/dtau_e.o tmp/dtau_h.o tmp/dtau_v.o tmp/database_mol.o tmp/database_planet.o tmp/intwdl.o tmp/planck.o \
		tmp/dplanck.o tmp/read_line2.o tmp/read_parameters.o tmp/tau_lines1.o tmp/tau_lines2.o tmp/tau_lines3.o \
		tmp/conv_k.o tmp/avg_k.o tmp/choldc.o tmp/matrix_inv.o

#----------------------------------------------------------------------------------------------------------------------#
#---- All: compile all programs
#----------------------------------------------------------------------------------------------------------------------#
all : abundance transac calcul_direct abundance_debug transac_debug calcul_direct_debug
#----------------------------------------------------------------------------------------------------------------------#
#---- Programs
#----------------------------------------------------------------------------------------------------------------------#
transac : transfer_transac $(lib)
	$(F90) -c transac_v4.f90
	mv transac_v4.o tmp/
	$(F90) -o ../transac_v4 tmp/transac_v4.o  tmp/transfer_transac_v4.o  $(lib_o)

abundance : transfer_abundance $(lib)
	$(F90) -c abundances_v4.f90
	mv abundances_v4.o tmp/
	$(F90) -o ../abundances_v4 tmp/abundances_v4.o tmp/transfer_abundance_v4.o $(lib_o)

abundance_para: transfer_abundance $(lib)
	$(F90) -fopenmp -c abundances_v4.f90
	mv abundances_v4.o tmp/
	$(F90) -o ../abundances_para tmp/abundances_v4.o tmp/transfer_abundance_v4.o $(lib_o)

calcul_direct : transfer_calcul_direct $(lib)
	$(F90) -c calcul_direct_v4.f90
	mv calcul_direct_v4.o tmp/
	$(F90) -o ../calcul_direct_v4 tmp/calcul_direct_v4.o tmp/transfer_calcul_direct_v4.o $(lib_o)

transac_debug : transfer_transac $(lib)
	$(F90) -c transac_v4.f90
	mv transac_v4.o tmp/
	$(F90) $(CFLAGS) -o ../transac_v4_debug tmp/transac_v4.o tmp/transfer_transac_v4.o $(lib_o)

abundance_debug : transfer_abundance $(lib)
	$(F90) -c abundances_v4.f90
	mv abundances_v4.o tmp/
	$(F90) $(CFLAGS) -o ../abundances_v4_debug tmp/abundances_v4.o tmp/transfer_abundance_v4.o $(lib_o)

calcul_direct_debug : transfer_calcul_direct $(lib)
	$(F90) -c calcul_direct_v4.f90
	mv calcul_direct_v4.o tmp/
	$(F90) $(CFLAGS) -o ../calcul_direct tmp/calcul_direct_v4.o tmp/transfer_calcul_direct_v4.o $(lib_o)
#----------------------------------------------------------------------------------------------------------------------#
#---- File
#----------------------------------------------------------------------------------------------------------------------#
avg:
	$(F90) -c avg.f90
	mv avg.o tmp/

conv:
	$(F90) -c conv.f90
	mv conv.o tmp/

conv_wf:
	$(F90) -c conv_wf.f90
	mv conv_wf.o tmp/

dat_cont:
	$(F90) -c dat_cont.f90 lib_functions.f90
	mv dat_cont.o lib_functions.o tmp/

det_step:
	$(F90) -c det_step.f90
	mv det_step.o tmp/

dtau_e:
	$(F90) -c dtau_e.f90
	mv dtau_e.o tmp/

dtau_h:
	$(F90) -c dtau_h.f90
	mv dtau_h.o tmp/

dtau_v:
	$(F90) -c dtau_v.f90
	mv dtau_v.o tmp/

dplanck:
	$(F90) -c dplanck.f90
	mv dplanck.o tmp/

database_mol:
	$(F90) -c database_mol.f90
	mv database_mol.o tmp/

database_planet:
	$(F90) -c database_planet.f90
	mv database_planet.o tmp/

intwdl:
	$(F90) -c intwdl.f90
	mv intwdl.o tmp/

lib_functions:
	$(F90) -c lib_functions.f90
	mv lib_functions.o tmp/

planck:
	$(F90) -c planck.f90 lib_functions.f90
	mv planck.o lib_functions.o tmp/

read_line2:
	$(F90) -c read_line2.f90
	mv read_line2.o tmp/

read_parameters:
	$(F90) -c read_parameters.f90
	mv read_parameters.o tmp/

tau_lines1:
	$(F90) -c tau_lines1.f90
	mv tau_lines1.o tmp/

tau_lines2:
	$(F90) -c tau_lines2.f90
	mv tau_lines2.o tmp/

tau_lines3:
	$(F90) -c tau_lines3.f90
	mv tau_lines3.o tmp/

transfer_transac:
	$(F90) -c transfer_transac_v4.f90
	mv transfer_transac_v4.o tmp/

transfer_abundance:
	$(F90) -c transfer_abundance_v4.f90
	mv transfer_abundance_v4.o tmp/

transfer_calcul_direct:
	$(F90) -c transfer_calcul_direct_v4.f90
	mv transfer_calcul_direct_v4.o tmp/

conv_k:
	$(F90) -c conv_k.f90
	mv conv_k.o tmp/

avg_k:
	$(F90) -c avg_k.f90
	mv avg_k.o tmp/

choldc:
	$(F90) -c choldc.f90
	mv choldc.o tmp/

matrix_inv:
	$(F90) -c matrix_inv.f90
	mv matrix_inv.o tmp/

#!help : Affiche cette aide.
#!
help:
	@grep -E "^#!" Makefile | sed -e 's:#!::g'

#!clean : efface les fichiers temporaires *.o
clean:
	rm -f *.o
