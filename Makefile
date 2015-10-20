FILES= consts lun_management math Nose_Hoover NormalModes MPItimer pot_mod system_mod Langevin main_stuff force_calc InputOutput nasa_mod find_neigh nasa ewald pot_spc pot_ttm potential smear PIMD

OBJS=$(addsuffix .o, $(FILES))

FC=mpif90 

WORK_DIR=$(shell pwd)
SIESTA_DIR=/home/dan/Dropbox/RESEARCH/SIESTA_myMonomerCorrected

SRC_DIR=$(SIESTA_DIR)/Src
OBJDIR=Obj
OBJ_DIR=$(SIESTA_DIR)/$(OBJDIR)
#
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH=$(SRC_DIR)
#
default: PIMD.x
#
ARCH_MAKE=$(OBJ_DIR)/arch.make
include $(ARCH_MAKE)
# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)


#
# Siesta libraries used by MPI version
#
SIESTA_LIB=libSiestaForces.a
FDF=libfdf.a
XMLPARSER=libxmlparser.a
XC=libSiestaXC.a
FoX_LIBS=`$(OBJ_DIR)/FoX/FoX-config --libs --wcml`
ALL_LIBS= $(SIESTA_LIB) $(FDF) $(WXML) $(XMLPARSER) $(XC) \
	$(MPI_INTERFACE) $(COMP_LIBS) $(FoX_LIBS) $(LIBS)
#
libs_collected:
	(cd $(OBJ_DIR) ; \
	make libSiestaForces.a ; \
	cp -f *.a *.mod $(WORK_DIR) )
	touch libs_collected

#FFLAGS=-fpp  -O3 -C -debug -traceback 
FFLAGS = -O3  -cpp -Dsiesta -Dparallel -DFC_HAVE_FLUSH -DFC_HAVE_ABORT 

all: PIMD.x

#compilation for debugging
debug: FFLAGS += --debug --backtrace -fbounds-check 
debug: PIMD.x

#compilation for profiling
profile: FFLAGS += -p
profile: PIMD.x 

#serial compilation commands
serial: FC=ifort
serial: FFLAGS = -O3  -cpp -Dsiesta -DFC_HAVE_FLUSH -DFC_HAVE_ABORT 
serial: PIMD_serial.x 


PIMD.x: libs_collected $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS) -o  ./PIMD.x  $(ALL_LIBS)

PIMD_serial.x: MPI.o $(OBJS)
	$(FC) $(FFLAGS) MPI.o  $(OBJS) -o  ./PIMD_serial.x

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
	
%.o : %.F90
	$(FC) $(FFLAGS) -c $< -o $@
	
MPI.o:  
	$(FC) $(FFLAGS) -c MPI.f90 -o MPI.o

#
clean: 
	@echo "==> Cleaning intermediate files ONLY"
	rm -f *.a *.o *.mod libs_collected
#
pristine: 
	@echo "==> Cleaning all intermediate and executable files"
	rm -f protoNEB
	rm -f *.a *.o *.mod libs_collected
