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

	
	endif!if (pid .eq. 0) 
				
	!---------------------------------- Start MPI Force calculation ----------------------------
	!we want to send a (3 x Natoms) array to each processor 
	do bat = 0, Nbatches - 1 !batch index

	if (pid .eq. 0) then

		do i = 1, Nnodes - 1 !node index
			k = bat*Nnodes + i  !bead index
			Call MPI_Send(RRt(:,:,k), counti, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, ierr)   
		enddo
		!masternode force calculation
		k = bat*Nnodes + Nnodes  

		call potential(RRt(:,:,k), Upott(k), dRRt(:,:,k), virt, dip_momIt(:,:,k), dip_momEt(:,:,k), chg, t, BAROSTAT)

		!recieve stuff from nodes
		do i = 1, Nnodes - 1
			k = bat*Nnodes + i  
			!masternode receive derivatives
			call MPI_Recv(dRRt(:,:,k), counti, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
			!masternode recieve energies
 			call MPI_Recv(Upott(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)	
			!masternode recieve dipole moments		
			if (dip_out .or. TD_out) call MPI_Recv(dip_momIt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
			if (Edip_out) call MPI_Recv(dip_momEt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		enddo
	else
		!slavenode recieve coords from master node
		call MPI_Recv(RR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
		!slavenode force calculation
		call potential(RR, Upot, dRR, virt, dip_momI, dip_momE, chg, t, BAROSTAT)
		!slavenode send back derivatives
		call MPI_Send(dRR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
		!slavenode send back energies
		call MPI_Send(Upot, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
		!slavenode send back dipole moments  
		if (dip_out .or. TD_out) call MPI_Send(dip_momI, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		if (Edip_out) call MPI_Send(dip_momE, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
	endif
	Call MPI_Barrier(MPI_COMM_WORLD, ierr)


	enddo! j = 0, Nbatches - 1 !batch index
	!---------------------------------- End MPI force calculation ----------------------------
	if (pid .eq. 0) then

		!calculate centroid positions
		RRc = sum(RRt,3)/Nbeads

		!check PBCs
		call PBCs(RRt, RRc)

		!update momenta a half step w/ new forces
		PPt = PPt - MASSCON*dRRt*delt2

		!calculate centroid momenta
		PPc = sum(PPt,3)/Nbeads 

		!update kinetic energy 
		uk = 0
		do i = 1,Nwaters
			uk = uk + imassO*sum( PPc(:,3*i-2)**2 )
			uk = uk + imassH*sum( PPc(:,3*i-1)**2 )
			uk = uk + imassH*sum( PPc(:,3*i-0)**2 )
		enddo	
		uk = .5d0*uk 

		!Propagate NH chains 
		if (BEADTHERMOSTAT) call bead_thermostat
	
		if (THERMOSTAT)     then 
			call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms, temp)
			PPt = PPt*s
		endif

		call MPItimer(3, 'start', secondsIO)
		
		!write stuff out if necessary 
		call write_out

		call MPItimer(3,'stop ',secondsIO)

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
