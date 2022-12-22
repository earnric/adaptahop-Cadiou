#
# Makefile for HaloMaker
#
######
# 4/ OPTIONS:
#    --------------------------------
# --- Use the -DBIG_RUN option for a normal (as opposed to resimulation) simulation whatever the i/o format. SHOULD NOT BE USED FOR STARS TODO REMOVE THIS IN THE CODE 
# --- Use the -DRENORM option if you are analyzing the DM part of a hydro simulation and want to renormalize 
# the mass of your dark matter particle so as to include the baryons (by pretending they are distributed like the DM):
# this might be useful when one wants to compare a pure DM simulation to its DM + hydro counterpart
# --- Use the -DANG_MOM_OF_R option if you want the code to output the angular momentum of concentric shells for your DM halos
# --- BN98 : if set, use a redshift-dependent value for rho_threshold to identify halos in adaptaHOP, 
#            this value is proportional to the critical overdensity given by the fit of Bryan&Norman(1998)  
#            if set, the value of rho_threshold given in the input file parameters is not used
# --- SIMPLE_OUTPUTS : output simple stuff for debugging and checks.
# --- CONTAM : compute and output contamination flag for each (sub-)halo -> tells whether a halo contains a low-res particle. This is for DM zooms. 
# --- TYPE: 0 for all, 1 for DM, 2 for stars
#
######

COMPILO = "gfortran"
OPTIONS = -DBN98 -DANG_MOM_OF_R #-DSIMPLE_OUTPUTS #-DIOGROUPSIZEREP=1200
TYPE    = 1
BIGRUN  = 1
CONTAM  = 0
LONGINT = 0
OLD     = 0

ifeq ($(TYPE),0)
 OPTIONS += -DSTARS -DALLPARTS
 EXEC = GroupFinder
else
 ifeq ($(TYPE),2)
  OPTIONS += -DSTARS
  EXEC = GalaxyFinder
 else
  EXEC = HaloFinder
 endif
endif

ifeq ($(CONTAM),1)
 OPTIONS += -DCONTAM
 CEXEC +=_contam
endif

ifeq ($(OLD),1)
 OPTIONS += -DOLDRAMSES
 CEXEC +=_nofamilies
endif

ifeq ($(BIGRUN),1)
 OPTIONS += -DBIG_RUN
 CEXEC +=_BR
endif

ifeq ($(LONGINT),1)
 OPTIONS += -DLONGINT
endif

SPACE := $(EMPTY) $(EMPTY)
EXECtmp = $(EXEC) $(CEXEC)
EXECNAME := $(subst $(SPACE),$(EMPTY),$(EXECtmp))

ifeq ($(COMPILO),"gfortran")
F90C    = gfortran -fallow-argument-mismatch
DEBUG   = -g -C
FFLAGS  = ${OPTIONS} -fbounds-check -fbacktrace -cpp ${DEBUG}
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
 F90C    = ifort -qopenmp
 #DEBUG   = -O0 -traceback -check bounds -cpp -debug -check
 DEBUG   = -O2 -g -traceback -check bounds -cpp -free -warn
 #OPTIM   = -O3 -qopt-mem-layout-trans=4 -qopt-prefetch -g -traceback -xCORE-AVX2 -cpp -free -ip -ipo -heap-arrays 0 #-qopt-report #-guide
 OPTIM   = -O3 -g -traceback -xAVX -cpp -free -heap-arrays #-qopt-report #-guide
 #OPTIM   = -fast
 FFLAGS  = ${OPTIONS} ${OPTIM}
 #FFLAGS  = ${OPTIONS} ${DEBUG} # RS DEBUG
endif


.SUFFIXES: .s .o .f90 

#
# Rules:
#

#.f90.o:
#	$(F90C) $(FFLAGS) -c $*.f90

%.o: %.f90
	$(F90C) $(FFLAGS) -c $*.f90

default:	HaloFinder snap2mass clean 


#
# erase useless files:
#

clean:	
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f fort.*
	rm -f nohup.out

input_output.o: halo_defs.f90

fof_mod.o: halo_defs.f90

compute_neiKDtree_mod.o: halo_defs.f90 stack.f90

# utils.o: halos_defs.f90

compute_halo_props.o: halo_defs.f90 utils.f90 input_output.f90 fof_mod.f90 compute_neiKDtree_mod.f90

HaloMaker.o: compute_halo_props.f90

HaloFinder:  halo_defs.o utils.o num_rec.o stack.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o  HaloMaker.o
	${F90C} ${FFLAGS} ${LIBS}  halo_defs.o utils.o num_rec.o stack.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o HaloMaker.o -o ${EXECNAME}

snap2mass: halo_defs.o utils.o subbox.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o snap2mass.o
	${F90C} ${FFLAGS} halo_defs.o utils.o subbox.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o snap2mass.o -o snap2mass	${clean}

# HaloFinder_BR:  halo_defs.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o  HaloMaker.o
# 	${F90C} ${FFLAGS} -DBIG_RUN ${LIBS} halo_defs.o num_rec.o input_output.o fof_mod.o compute_neiKDtree_mod.o compute_halo_props.o HaloMaker.o -o HaloFinder_BR 
