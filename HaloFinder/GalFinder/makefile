#
# Makefile for GalaxyMaker
#
######
# 1/ RAMSES related OPTIONS:
#    ---------------------
# Use the -DMETALS flag if RAMSES has computed metallicity from stellar feedback and you are reading version 3.0 i/o format
######
# 2/ GADGET related OPTION:
#    ---------------------
# Use the -DTest_Gd flag for testing Gadget reading routines only, especially a resimulation or hydro. That's a debugging option use to use if GalaxyMaker cannot read a gadget snapshot.
######
# 4/ Other simulation related OPTIONS:
#    --------------------------------
# --- Use the -DBIG_RUN option for a normal (as opposed to resimulation) simulation whatever the i/o format. SHOULD NOT BE USED FOR STARS TODO REMOVE THIS IN THE CODE 
# --- Use the -DRENORM option if you are analyzing the DM part of a hydro simulation and want to renormalize 
# the mass of your dark matter particle so as to include the baryons (by pretending they are distributed like the DM):
# this might be useful when one wants to compare a pure DM simulation to its DM + hydro counterpart
# --- Use the -DANG_MOM_OF_R option if you want the code to output the angular momentum of concentric shells for your DM halos
######

COMPILO = "ifort"
OPTIONS = -DMETALS -DBIG_RUN #-DRENORM -DANG_MOM_OF_R
EXEC = GalFinder


ifeq ($(COMPILO),"gfortran")
F90C    = gfortran
DEBUG   = -g -C
FFLAGS  = ${OPTIONS} -fbounds-check -fbacktrace -O2 -cpp
endif

ifeq ($(COMPILO),"pgf90")
# options when compiling with pgf90 (meso, horizon1(k8-64), queeg)
F90C    = pgf90
DEBUG   = -g -C
OPTIM   = -O3 -mp=align -fast -fastsse -Mnontemporal -Minfo #-Mmovnt -Minfo #-Mvect=fuse,nosizelimit,prefetch,sse 
FFLAGS  = ${OPTIONS} ${OPTIM} -Mpreprocess -tp k8-64 -mcmodel=medium -Mlarge_arrays -Mextend 
endif

ifeq ($(COMPILO),"ifort")
# option when compiling with ifort (crozes, orac)
 F90C    = ifort
 DEBUG   = -g -fpe0 -O0 -check bounds -traceback #-debug -check -warn all -stand f90 
 # OPTIM   = -O2 -traceback -foptimize-sibling-calls -fopenmp #-fast
 OPTIM   = -O3 -xAVX -cpp -free -fopenmp -foptimize-sibling-calls
 FFLAGS  = ${OPTIONS} ${OPTIM} -free -m64 -cpp 
 #FFLAGS  = ${OPTIONS} ${DEBUG} -cpp 
endif

# need these libs to read enzo HDF4 format
LIBS    = # -L/sw/lib -ljpeg -lz -lsz -lmfhdf -ldf 


SPACE := $(EMPTY) $(EMPTY)
EXECtmp = $(EXEC) $(CEXEC)
EXECNAME := $(subst $(SPACE),$(EMPTY),$(EXECtmp))


.SUFFIXES: .o .f90 

# Rules:
%.o: %.f90
	$(F90C) $(FFLAGS) -c $*.f90

default:	GalaxyMaker snap2mass clean


#
# erase useless files:
#

clean:	
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f fort.*
	rm -f nohup.out

GalaxyMaker:  halo_defs.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o  GalaxyMaker.o
	${F90C} ${FFLAGS} ${LIBS}  halo_defs.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o GalaxyMaker.o -o ${EXECNAME}

input_output.o: halo_defs.f90

fof_mod.o: halo_defs.f90

compute_neiKDtree_mod.o: halo_defs.f90

compute_halo_props.o: halo_defs.f90 input_output.f90  fof_mod.f90 compute_neiKDtree_mod.f90

GalaxyMaker.o: compute_halo_props.f90

snap2mass: halo_defs.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o snap2mass.o
	${F90C} ${FFLAGS} halo_defs.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o snap2mass.o -o snap2mass ${clean}
