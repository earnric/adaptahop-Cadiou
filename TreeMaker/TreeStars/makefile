#
# Makefile for TreeMaker
#
# flags for f90 compilation on crozes
#

# Possible options: 
# -----------------
# - BIG_RUN : for a normal simulation (as opposed to resimulation)

#MACHINE = "meso"
MACHINE = "ifort"

OPTIONS =  #-DBIG_RUN

ifeq ($(MACHINE), "meso") 
F90C    = pgf90
DEBUG   = -Minform=warn -g -C -Mnofprelaxed
OPTIM   = -O2 -Mnofprelaxed
FFLAGS  = -tp k8-64 -mcmodel=medium -Mlarge_arrays -Mpreprocess -Mextend ${OPTIM}
endif	

ifeq ($(MACHINE), "demo") 
F90C    = ifort
DEBUG   = -g -CB -fp-stack-check -fpe0 -traceback -warn #-pg
OPTIM   = -O2 -fp-model precise 
FFLAGS  =  -cpp ${OPTIM}
endif	

ifeq ($(MACHINE), "ifort") 
F90C    = ifort
DEBUG   = -g -CB -fpe0 -traceback -warn all -warn unused #-pg -fp-stack-check
OPTIM   = -O3 -xHost -fp-model precise 
MIXED   = -O3 -xHost -g -traceback -fp-model precise 
FFLAGS  =  -cpp ${MIXED}
endif	


.SUFFIXES: .o .f90 
# Rules:
%.o: %.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $*.f90

default:        TreeMaker clean


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

TreeMaker:  tree_defs.o input_output.o compute_merger_tree.o TreeMaker.o
	${F90C} ${FFLAGS} ${OPTIONS} tree_defs.o input_output.o compute_merger_tree.o TreeMaker.o -o TreeMaker

TreeMaker_bigrun:  tree_defs.o input_output.o compute_merger_tree.o TreeMaker.o
	${F90C} ${FFLAGS} ${OPTIONS} -DBIG_RUN tree_defs.o input_output.o compute_merger_tree.o TreeMaker.o -o TreeMaker_bigrun

input_output.o : tree_defs.f90

compute_merger_tree.o : tree_defs.f90 input_output.f90

TreeMaker.o : tree_defs.f90 input_output.f90 compute_merger_tree.f90

