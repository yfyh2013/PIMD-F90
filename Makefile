MODULES=$(addprefix ./, consts.o  Nose_Hoover.o NormalModes.o MPItimer.o pot_mod.o system_mod.o main_stuff.o nasa_mod.o calcRij.o find_neigh.o math.o nasa.o ewald.o pot_spc.o pot_ttm.o potential.o smear.o minimize.o main_p.o  )

SRC = $(MODULES)

OBJ=$(SRC:.f90=.o)

#FC=ifort
#FC=gfortran
FC=mpif90 

#FFLAGS=-fpp  -O3 -C -debug -traceback -ftz 
#FFLAGS=  -g -O3 -ip  #intel compiler only 
FFLAGS = -O3 #--debug --backtrace

all: main.x

$(OBJ) : Makefile



main.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o  ./main.x

%.o:%.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean: 
	rm -rf $(OBJ) *.mod main.x

