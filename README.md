#Path Integral Molecular Dynamics in Fortran 90

2015 Daniel C. Elton 

PIMD-F90 implements path integral molecular dynamics in Fortran-90.

By changing the adiabaticity parameter (facticious mass rescaling) the code can also implement Centroid Molecular Dynamics (CMD) and Ring Polymer Molecular Dynamics (RPMD)

See the example input file input256.inp for program options

The current repo includes routines for implementation of the TTM2-F, TTM3-F, and spc-like potentials.

It also can use the SIESTA code for DFT force calculation, using pipes to communicate with the siesta process. 

Other features: 
-- Parallelization over the beads with MPI
-- Options for a global Nose-Hoover thermostat or Anderson barostat 
-- virial energy and pressure estimators
-- many output options 
-- Choice of Langevin or Nose-Hoover bead thermostats 
-- Option to not thermostat the centroid 
-- Option to perform the PIMD on the monomer only (a subcase of the ring polymer contraction scheme, ref: JCP 129, 024105 (2008))


The TTM / SPC force calculation routines are modification of codes written by C. J. Burnham, G. S. Fanourgakis, and S. S. Xantheas
(potential.f90, consts.f90, find_neigh.f90, init_ttm21f.f90, pot_mod.f90, pot_ttm.f90, smear.f90 system_mod.f90 math.f90 ewald.f90 nasa.f90 nasa_mod.f90)
see https://code.google.com/p/nqcbabel/source/browse/trunk/extra/?r=58#extra%2Fttm-ice
TTM-2F ref: J. Phys. Chem. A. vol 110, page 4100-4106, (2006) 
