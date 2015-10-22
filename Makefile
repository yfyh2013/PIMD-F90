FILES=pxf fsiesta_pipes dans_timer consts lun_management math Nose_Hoover NormalModes pot_mod system_mod Langevin main_stuff force_calc InputOutput nasa_mod find_neigh nasa ewald pot_spc pot_ttm potential smear PIMD

OBJS=$(addsuffix .o, $(FILES))

#VPATH=/home/dan/Dropbox/RESEARCH/SIESTA_myMonomerCorrected/Src

FC=mpif90 

##Note : -DFC_HAVE_FLUSH -DFC_HAVE_ABORT  flags are for SIESTA pipes support  

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
serial: FC=gfortran #ifort not working!
serial: FFLAGS = -O3  -cpp -Dsiesta -DFC_HAVE_FLUSH -DFC_HAVE_ABORT --debug --backtrace -fbounds-check 
serial: PIMD_serial.x 


PIMD_serial.x: MPI.o $(OBJS)
	$(FC) $(FFLAGS) MPI.o  $(OBJS) -o  ./PIMD_serial.x

PIMD.x: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o  ./PIMD.x

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@
	
%.o : %.F90
	$(FC) $(FFLAGS) -c $< -o $@
	
MPI.o:  
	$(FC) $(FFLAGS) -c MPI.f90 -o MPI.o

clean: 
	rm -rf *.o *.mod 
	
pristine: 
	rm -rf *.o *.mod INPUT_TMP* *.log *.xml *.x 