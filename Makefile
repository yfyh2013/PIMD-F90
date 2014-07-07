FILES=consts math Nose_Hoover NormalModes MPItimer pot_mod system_mod Langevin main_stuff force_calc InputOutput nasa_mod find_neigh nasa ewald pot_spc pot_ttm potential smear main_p

OBJS=$(addsuffix .o,   $(FILES))

FC=mpif90 

#FFLAGS=-fpp  -O3 -C -debug -traceback -ftz 
FFLAGS = -O3

all: main.x

debug: FFLAGS += --debug --backtrace -fbounds-check
debug: main.x

profile: FFLAGS += -p
profile: main.x 

#serial compilation commands
serial: FC=gfortran 
serial: serial.x 

serial.x: MPI.o $(OBJS)
	$(FC) $(FFLAGS) MPI.o  $(OBJS) -o  ./main.x

main.x: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o  ./main.x

%.o:%.f90
	$(FC) $(FFLAGS) -c $< -o $@

MPI.o:  
	$(FC) $(FFLAGS) -c MPI.f90 -o MPI.o

clean: 
	rm -rf *.o *.mod 

