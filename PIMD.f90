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
use NormalModes
use mpi
use force_calc
use dans_timer

implicit none
!main variables can be found in main_stuff.f90

!---------------------- Read in input file ----------------------------------------
call read_input_file

!---------------------- Start MPI and find number of nodes and pid --------------
call MPI_Init(ierr)
if (.not. (ierr .eq. 0))	write(*,*) "WARNING: MPI did not intialize correctly."
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

!---------------------------------------------------------------------------------
!----------- MD Part of the program ---------------------------------------------
!---------------------------------------------------------------------------------
call start_timer("MD")
do t = 1, num_timesteps + eq_timesteps

	!----- Velocity-Verlet integration -------
	if (pid .eq. 0) then
		!Propagate NH chains 
		if (BEADTHERMOSTAT) call bead_thermostat
		
		if (THERMOSTAT)  then 
			call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
			PPt = PPt*s
		endif

		!update momenta a half step w/ old forces
		PPt = PPt - MASSCON*dRRt*delt2
		
		!Normal modes stuff
	    if (.not. (CONTRACTION)) then
	
			if (Nbeads .gt. 1) then
				do i = 1, Nwaters
					!Evolve ring a full step 
					Call EvolveRing(RRt(:,3*i-2,:), PPt(:,3*i-2,:), Nbeads, massO)
					Call EvolveRing(RRt(:,3*i-1,:), PPt(:,3*i-1,:), Nbeads, massH)
					Call EvolveRing(RRt(:,3*i-0,:), PPt(:,3*i-0,:), Nbeads, massH)
				enddo
			else 
				do i = 1,Nwaters
					do k = 1,Nbeads
						RRt(:,3*i-2,k) = RRt(:,3*i-2,k) + imassO*PPt(:,3*i-2,k)*delt
						RRt(:,3*i-1,k) = RRt(:,3*i-1,k) + imassH*PPt(:,3*i-1,k)*delt
						RRt(:,3*i-0,k) = RRt(:,3*i-0,k) + imassH*PPt(:,3*i-0,k)*delt
					enddo
				enddo			
			endif

			!calculate centroid positions
			RRc = sum(RRt,3)/Nbeads

			!calculate centroid momenta
			PPc = sum(PPt,3)/Nbeads 

			!check PBCs
			call PBCs(RRt, RRc)

		endif!.not. (CONTRACTION)
	
	endif!if (pid .eq. 0) 
	!------ call force routine ------------------------------------ 
	if (CONTRACTION .eqv. .false.) then	
		call full_bead_forces
	else 	
		call contracted_forces
	endif 
	!------ end call force routine ------------------------------------ 
	if (pid .eq. 0) then

		!update kinetic energy 
		call calc_uk 

		!write stuff out if necessary 
		call start_timer("WritingOut")
		call write_out
		call stop_timer("WritingOut")
		
		!update momenta a half step w/ new forces
		PPt = PPt - MASSCON*dRRt*delt2

		!Propagate NH chains 
		if (BEADTHERMOSTAT) call bead_thermostat
		
		call calc_uk 

		if (THERMOSTAT)     then 
			call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
			PPt = PPt*s
		endif
		
		if (BAROSTAT) call Pcouple


	endif!(pid .eq. 0) 

enddo!t = 1, num_timesteps + eq_timesteps
call stop_timer("MD")


call shutdown 


end program PIMD
