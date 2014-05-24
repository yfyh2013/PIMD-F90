program ttm3p
!----------------------------------------------------------------------------------
!----------------- PIMD / RPMD program 
!Changelog: 
!TTM3F code edited and expanded by Dan Elton 2/13
!Added PIMD parallelization 12/13
!Switched from velocity to momentum variables 3/14
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
implicit none
!Variables for main can be found in main_stuff.f90

!---------------------- Read in input file ----------------------------------------
call read_input_file 

!---------------------- Start MPI and find number of nodes and pid --------------
call MPI_Init(ierr)
if (.not. (ierr .eq. 0))	write(*,*) "WARNING: MPI did not intialize correctly"
call MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)

!---------------- Initialize / allocate some variables -------------------------
call initialize_all_node_variables
call init_pot

!------------------Master node stuff --------------------------------------------
if (pid .eq. 0) then
	call master_node_init
	call read_coords           ! Read in coord data to RRc 
	call initialize_beads      ! Initialize RRt
	call initialize_velocities ! Initialize PPt
	call calc_uk 		    ! Initalize kinetic energy
	call open_files
	call MPItimer(1,'start',seconds)
endif!(pid .eq. 0)

!---------------------------------------------------------------------------------
!----------- MD Part of the program ---------------------------------------------
!---------------------------------------------------------------------------------
do t = 1, num_timesteps + eq_timesteps
	!----- Verlet integration -------
	if (pid .eq. 0) then

		!Propagate NH chains 
		if (BEADTHERMOSTAT) call bead_thermostat
		
		if (THERMOSTAT)  then 
			call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms, temp)
			PPt = PPt*s
		endif

		!update momenta a half step w/ old forces
		PPt = PPt - MASSCON*dRRt*delt2
		
		!Normal modes stuff
		call MPItimer(2,'start',secondsNM)

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

		call MPItimer(2, 'stop ', secondsNM)

		!calculate centroid positions
		RRc = sum(RRt,3)/Nbeads

		!calculate centroid momenta
		PPc = sum(PPt,3)/Nbeads 

		!check PBCs
		call PBCs(RRt, RRc)
	
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
		call MPItimer(3, 'start', secondsIO)
		call write_out
		call MPItimer(3,'stop ',secondsIO)

		!update momenta a half step w/ new forces
		PPt = PPt - MASSCON*dRRt*delt2

		!Propagate NH chains 
		if (BEADTHERMOSTAT) call bead_thermostat

		call calc_uk_centroid

		if (THERMOSTAT)     then 
			call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms, temp)
			PPt = PPt*s
		endif
		
		if (BAROSTAT) call Pcouple


	endif!(pid .eq. 0) 

enddo!t = 1, num_timesteps + eq_timesteps

if (pid .eq. 0)  then
	call MPItimer(1,'write',seconds) 
	call MPItimer(2,'write',secondsNM) 
	call MPItimer(3,'write',secondsIO) 
	call print_run
endif

!if (pid .eq. 0) then
!	deallocate(RRc)
!	deallocate(VVc)
!	deallocate(RRt)
!	deallocate(VVt)
!	deallocate(RRc)
!	deallocate(RRdev)
!endif 

Call MPI_Barrier(MPI_COMM_WORLD, ierr)
Call MPI_Finalize(ierr)

end program 
