FILES=consts lun_management math Nose_Hoover NormalModes MPItimer pot_mod system_mod Langevin main_stuff force_calc InputOutput nasa_mod find_neigh nasa ewald pot_spc pot_ttm potential smear PIMD

OBJS=$(addsuffix .o,   $(FILES))

FC=mpif90 

#FFLAGS=-fpp  -O3 -C -debug -traceback -ftz 
FFLAGS = -O3

all: PIMD.x

debug: FFLAGS += --debug --backtrace -fbounds-check
debug: PIMD.x

profile: FFLAGS += -p
profile: PIMD.x 

#serial compilation commands
serial: FC=gfortran 
serial: PIMD_serial.x 

PIMD_serial.x: MPI.o $(OBJS)
	$(FC) $(FFLAGS) MPI.o  $(OBJS) -o  ./PIMD_serial.x

PIMD.x: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o  ./PIMD.x

%.o:%.f90
	$(FC) $(FFLAGS) -c $< -o $@

MPI.o:  
	$(FC) $(FFLAGS) -c MPI.f90 -o MPI.o

clean: 
	rm -rf *.o *.mod 
