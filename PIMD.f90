program PIMD
!----------------------------------------------------------------------------------
!----------------- PIMD / RPMD program -------------------------------------------
!----------------------------------------------------------------------------------
!Changelog: 
!TTM3F code edited and expanded by Dan Elton 2/13
!Added PIMD parallelization 12/13
!Switched from velocity to momentum variables 3/14
!Addition of SIESTA 10/15
!
!To run the program use the following : 
! serial   : ./main.x inputfilename
! parallel :  mpirun -np _ ./main.x inputfilename
!----------------------------------------------------------------------------------
 use system_mod
 use consts
 use main_stuff
 use InputOutput
 use mpi
 use MD
 use dans_timer

 implicit none
!main variables can be found in main_stuff.f90

!---------------------- Read in input file ----------------------------------------
 call read_input_file

!---------------------- Start MPI and find number of nodes and pid --------------
 call MPI_Init(ierr)
 if (.not. (ierr .eq. 0)) write(*,*) "WARNING: MPI did not intialize correctly."
 call MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)

!---------------- Initialize / allocate some variables -------------------------
 call read_and_initialize_all_nodes
 call init_potential_all_nodes

!------------------Master node stuff --------------------------------------------
 if (pid .eq. 0) then
	call start_timer("Total time") 
	call open_files
	call master_node_init
	call read_coords_and_init  ! Read in coord data to RRc 
	call calc_uk 		       ! Initalize kinetic energy
 endif!(pid .eq. 0)

 !----------- MD Part of the program ---------------------------------------------
 call start_timer("MD")
 
 do t = 1, run_timesteps + eq_timesteps
	
	if (CONTRACTION) then
		call contracted_MD
	else if (MONOMERPIMD) then 
		call monomer_PIMD
	else
		call PIMD_VV
	endif
	
 enddo!t = 1, run_timesteps + eq_timesteps
 
 call stop_timer("MD")

 call shutdown 

end program PIMD