#
# Makefile for TreeMaker
#
# flags for f90 compilation on crozes
#

# Possible options: 
# -----------------
# - BIG_RUN : for a normal simulation (as opposed to resimulation)
# - SIMPL   : Simplification of the merger tree (one son per (sub)halo
# 
# Tree simplification options: (SIMPL option needed)
# - SDAD_FOR_SUB  : recover progs for subhalos which doesn't have any
# - R_DADLESS_SUB : remove subhalos with no progenitor (will be considered as fragments by default)
# - COL_TREE      : remove all subhalos except when it solves a fly-by issue

#MACHINE = "meso"
COMPILO = "ifort"

OPTIONS = 

SIMPL_OPTIONS = -DSIMPL -DS_DAD_FOR_SUB #-DR_DADLESS_SUB 

ifeq ($(COMPILO), "pgf90") 
F90C    = pgf90
DEBUG   = -Minform=warn -g -C -Mnofprelaxed
OPTIM   = -O2 -Mnofprelaxed
FFLAGS  = -tp k8-64 -mcmodel=medium -Mlarge_arrays -Mpreprocess -Mextend ${OPTIM}
endif	

ifeq ($(COMPILO), "ifort") 
F90C    = ifort
DEBUG   = -g -CB -fpe0 -traceback -warn all -warn unused #-pg -fp-stack-check
OPTIM   = -O3 -g -traceback -xAVX -fp-model precise 
FFLAGS  =  -cpp ${OPTIM}
endif	

ifeq ($(COMPILO),"gfortran")
F90C    = gfortran -ffree-line-length-none ## -ffree-form -fopenmp
DEBUG   = -g -ffpe-trap=invalid -fbacktrace -fbounds-check -frange-check -Wall
OPTIM   = -O2
FFLAGS  = -cpp ${OPTIM}
endif

.SUFFIXES: .o .f90 
# Rules:
# %.o: %.f90
# 	$(F90C) $(FFLAGS) $(OPTIONS) ${SIMPL_OPTIONS} -c $*.f90
%.o: %.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $*.f90

normal:         TreeMaker clean TreeMaker_BR clean

#
# erase useless files:
#

clean: 	
	rm -f *.o
	rm -f *.mod
	rm -f *~

clean_data:
	rm -f tree.dat
	rm -f halos_results.???
	rm -f tree_file_???.???
	rm -f file_master_log

TreeMaker:  tree_defs.o input_output.o compute_merger_tree.o  TreeMaker.o
	${F90C} ${FFLAGS} ${OPTIONS} tree_defs.o input_output.o compute_merger_tree.o TreeMaker.o -o TreeMaker

TreeMaker_BR:  tree_defs.o input_output.o compute_merger_tree.o  TreeMaker.o
	${F90C} ${FFLAGS} ${OPTIONS} -DBIG_RUN tree_defs.o input_output.o compute_merger_tree.o TreeMaker.o -o TreeMaker_BR


TreeMaker_simpl:  tree_defs.o input_output.o compute_merger_tree.o simpl_merger_tree.o jeje.o TreeMaker.o
	${F90C} ${FFLAGS} ${OPTIONS} ${SIMPL_OPTIONS} tree_defs.o input_output.o compute_merger_tree.o simpl_merger_tree.o jeje.o TreeMaker.o -o TreeMaker_simpl


TreeMaker_bigrun_simpl:  tree_defs.o input_output.o compute_merger_tree.o simpl_merger_tree.o jeje.o TreeMaker.o
	${F90C} ${FFLAGS} -DBIG_RUN ${OPTIONS} ${SIMPL_OPTIONS} tree_defs.o input_output.o compute_merger_tree.o simpl_merger_tree.o jeje.o TreeMaker.o -o TreeMaker_bigrun

input_output.o : tree_defs.f90 

compute_merger_tree.o : tree_defs.f90 input_output.f90

simpl_merger_tree.o : tree_defs.f90 

TreeMaker.o : tree_defs.f90 input_output.f90 compute_merger_tree.f90 simpl_merger_tree.f90

blah:	tree.o
	${F90C} ${FFLAGS} ${OPTIONS} ${SIMPL_OPTIONS} tree.o -o blah
 

