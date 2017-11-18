# Fast and accurate *ab-initio* Path Integral Molecular Dynamics in Fortran 90

2015-2017 Daniel C. Elton

PIMD-F90 implements path integral molecular dynamics (PIMD) to capture the nuclear quantum effects in water.

By changing the adiabaticity parameter (ficticious mass rescaling) the code can also implement Centroid Molecular Dynamics (CMD) and Ring Polymer Molecular Dynamics (RPMD).

See the example input file input128.inp for program options

## Supported potentials
### TTM2-F, TTM3-F, q-SPC/Fw
The TTM / SPC force calculation routines are modification of codes written by C. J. Burnham, G. S. Fanourgakis, and S. S. Xantheas
(potential.f90, consts.f90, find_neigh.f90, init_ttm21f.f90, pot_mod.f90, pot_ttm.f90, smear.f90 system_mod.f90 math.f90 ewald.f90 nasa.f90 nasa_mod.f90)
These are publically available on [Google code](https://code.google.com/p/nqcbabel/source/browse/trunk/extra/?r=58#extra%2Fttm-ice).
TTM-2F ref: J. Phys. Chem. A. vol 110, page 4100-4106, (2006)

### SIESTA
The PIMD code can obtain forces & energies from the density functional theory package [SIESTA](http://departments.icmab.es/leem/siesta/).
It uses unix pipes to communicate with SIESTA. A working flush() routine is required for this (see pflush.f90).


## Features
* Parallelization over the beads with MPI
* Options for a global Nose-Hoover thermostat or Anderson barostat
* virial energy and pressure estimators
* many output options
* Choice of Langevin or Nose-Hoover bead thermostats
* Option to not thermostat the centroid
* Option to perform the PIMD on the monomer only (a subcase of the ring polymer contraction scheme ref: [JCP 129, 024105 (2008)](http://dx.doi.org/10.1063/1.2953308) )
* implementation of custom PES for using the "monomer PIMD" method with DFT

## Input file specification
See [manual/PIMD_manual.pdf](manual/PIMD_manual.pdf)
