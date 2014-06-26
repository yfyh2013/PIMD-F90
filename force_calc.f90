module force_calc

contains
!--------------------------------------------------------------------------------------------
!---------------------------------- Default MPI force calculation --------------------------
!--------------------------------------------------------------------------------------------
subroutine full_bead_forces
 use main_stuff 
 Implicit None 

 do bat = 0, Nbatches - 1 !batch index

    !we want to send a (3 x Natoms) array to each processor 
    if (pid .eq. 0) then

 	do i = 1, Nnodes - 1 !node index
		k = bat*Nnodes + i  !bead index
		!masternode send image coordinates
		Call MPI_Send(RRt(:,:,k), counti, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, ierr)   
		!masternode send centroid coordinates to calculate centroid virial
		Call MPI_Send(RRc, 3*Natoms, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, ierr)   
	enddo
	!masternode force & virial calculation
	k = bat*Nnodes + Nnodes  

	call potential(RRt(:,:,k), RRc, Upott(k), dRRt(:,:,k), virt, virialct(k), dip_momIt(:,:,k), dip_momEt(:,:,k), chg, t, BAROSTAT)

	virialt(k) = virt(1,1) + virt(2,2) + virt(3,3)	

	!recieve stuff from nodes
	do i = 1, Nnodes - 1
		k = bat*Nnodes + i  
		!masternode receive derivatives
		call MPI_Recv(dRRt(:,:,k), counti, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		!masternode recieve energies
		call MPI_Recv(Upott(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)	
		!masternode receive virials
		call MPI_Recv(virialt(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)	
		!masternode receive centroid virials
		call MPI_Recv(virialct(k), 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)	
		!masternode recieve dipole moments		
		if (dip_out .or. TD_out) call MPI_Recv(dip_momIt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		if (Edip_out) call MPI_Recv(dip_momEt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
	enddo
    else
		!slavenode recieve coords from master node
		call MPI_Recv(RR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
		!slavenode recieve centroid coords from master node
		call MPI_Recv(RRc, 3*Natoms, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
		!slavenode force & virial calculation
		call potential(RR, RRc, Upot, dRR, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT)
		virial = virt(1,1) + virt(2,2) + virt(3,3)	
		!slavenode send back derivative
		call MPI_Send(dRR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
		!slavenode send back energy
		call MPI_Send(Upot, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
		!slavenode send back virial
		call MPI_Send(virial, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
		!slavenode send back centroid virial
		call MPI_Send(virialc, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
		!slavenode send back dipole moments 
		if (dip_out .or. TD_out) call MPI_Send(dip_momI, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
		if (Edip_out) call MPI_Send(dip_momE, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
    endif
    Call MPI_Barrier(MPI_COMM_WORLD, ierr)

 enddo! j = 0, Nbatches - 1 !batch index

if (pid .eq. 0) then
 Upot    = sum(Upott)    	!potential energy for the ENTIRE system
 virial  = sum(virialt)  	!virial for the ENTIRE system (kcal/mol)
 virialc = sum(virialct)/Nbeads !centroid virial for the ENTIRE system (kcal/mol)
endif 

end subroutine full_bead_forces


!--------------------------------------------------------------------------------------------
!----------------- Contracted force calculation with intermolecular forces on monomer-------
!--------------------------------------------------------------------------------------------
!-- This subroutine performs evaluates only the intramoleculari (fast) forces on the beads and 
!-- evaluates the intermolecular (slow) forces on the centroid. It also uses a multiple timestep method
!-- (See Tuckerman, pg 118). The momentum and positions will be updated with the intramolecular forces
!-- intra_timesteps times. For instance if the 'outer timestep' (normal timestep) is .5 ps and intra_timesteps = 5
!-- then the inner timestep is .1 ps

subroutine contracted_forces
 use main_stuff 
 use NormalModes
 Implicit None 
 double precision :: e1, virialmon, virialcmon, Umonomers
 double precision, dimension(3,Natoms,Nbeads)  ::  dRRfast
 double precision, dimension(3,Natoms) :: dRRc
 double precision, dimension(3,3)      :: dr1, r1
 double precision, dimension(3)        :: roh1, roh2, rh1m, rh2m, rM
 integer :: tintra

 if (pid .eq. 0) then

   Umonomers = 0 
   virialmon = 0 
   virialcmon = 0 
	   
   !---  intramolecular (fast) forces --------------------------------------------------- 
   do tintra = 1, intra_timesteps

	!update momenta with fast forces
        PPt = PPt - MASSCON*dRRfast*delt2fast

		call MPItimer(2,'start',secondsNM)
	!update positions with fast forces
	do i = 1, Nwaters
        	Call EvolveRing(RRt(:,3*i-2,:), PPt(:,3*i-2,:), Nbeads, massO)
        	Call EvolveRing(RRt(:,3*i-1,:), PPt(:,3*i-1,:), Nbeads, massH)
        	Call EvolveRing(RRt(:,3*i-0,:), PPt(:,3*i-0,:), Nbeads, massH)
	enddo
		call MPItimer(2, 'stop ', secondsNM)
	!update fast forces (intramolecular forces)
	!masternode calcuates the intramolecular forces, puts them in dRRfast

	Umonomers = 0 
	virialmon = 0 
	virialcmon = 0 
	do j = 1, Nbeads
		do iw = 1, Nwaters
			iO=3*iw-2; iH1 = 3*iw-1; iH2=3*iw-0

	   		r1(1:3, 1:3) = RRt(1:3, (/iO, ih1, ih2/), j)

		  	call pot_nasa(r1, dr1, e1, box, boxi)  

			dRRfast(1:3, (/iO, iH1, iH2/), j) = dr1

			!if last timestep in loop update energies 
			if (tintra .eq. intra_timesteps) then
				Umonomers = Umonomers + e1
			endif
		enddo
	enddo

        !update momenta with fast forces
        PPt = PPt - MASSCON*dRRfast*delt2fast

   enddo !tintra  = 1.. 

 !calculate centroid positions
 RRc = sum(RRt,3)/Nbeads

 !calculate centroid momenta
 PPc = sum(PPt,3)/Nbeads

 !check PBCs
 call PBCs(RRt, RRc)

 !intermolecular force calculation
  !if (Nnodes .gt. 1) then 
	!send the centroid coords out to be processed 
!	Call MPI_Send(RRc(:,:), counti, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, ierr)   
	!masternode recieve forces on centroid		
!	call MPI_Recv(dRRc, counti, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, status2, ierr)
	!masternode receive intermolecular energy 
!	call MPI_Recv(Upot, 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, status2, ierr)	
	!masternode receive virials
!	call MPI_Recv(virialt(k), 1, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, status2, ierr)		
	!masternode recieve centroid based dipole moments		
	!if (dip_out .or. TD_out) call MPI_Recv(dip_momIt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
	!if (Edip_out) call MPI_Recv(dip_momEt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
 !  else 
	!if only one node then masternode calculates forces on centroid
	call potential(RRc, RRc, Upot, dRRc, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT)
  ! endif

   !update dRRt
   do j = 1, Nbeads
	dRRt(:,:,j) = dRRc
   enddo

  !calculate dipole moments and virial
  do j = 1, Nbeads
	do iw = 1, Nwaters
		iO=3*iw-2; iH1 = 3*iw-1; iH2=3*iw-0

		!do centroid virial first
		roh1 = RRc(1:3, ih1) - RRc(1:3, iO)
		roh1 = roh1 - box*anint(roh1*boxi) !PBC
		roh2 = RRc(1:3, ih2) - RRc(1:3, iO)
		roh2 = roh2 - box*anint(roh2*boxi) !PBC

		virialcmon = virialcmon + dot_product(roh1, dr1(:,2)) 
		virialcmon = virialcmon + dot_product(roh2, dr1(:,3)) 


		!do virial
		roh1 = RRt(1:3, ih1, j) - RRt(1:3, iO, j)
       		roh1 = roh1 - box*anint(roh1*boxi) !PBC
		roh2 = RRt(1:3, ih2, j) - RRt(1:3, iO, j)
       		roh2 = roh2 - box*anint(roh2*boxi) !PBC

		virialmon = virialmon + dot_product(roh1, dr1(:,2))
		virialmon = virialmon + dot_product(roh2, dr1(:,3)) 

		!do dipole moments
		rM(:) = 0.5d0*0.46d0*( roh1(:) + roh2(:) ) + RRt(1:3, iO, j)

		rh1m = RRt(1:3, ih1, j) - rM(:)
  		rh2m = RRt(1:3, ih2, j) - rM(:)

		dip_momIt(:,iw,j) = chg(iH1)*rh1m + chg(iH2)*rh2m
	enddo
	!add polarization dipoles
   	dip_momIt(:,:,j) = dip_momIt(:,:,j) + dip_momE
   	dip_momEt(:,:,j) = dip_momE
   enddo
 
   !update Upot, virial and virialc
   Upot = Upot*Nbeads + Umonomers !potential energy for the ENTIRE system (all images)
   virial = virialmon + virt(1,1) + virt(2,2) + virt(3,3) !virial for the ENTIRE system (all images)
   virialc = virialcmon/Nbeads + virialc

 endif
 




!slave node stuff
 if (pid .eq. 1) then
	!slavenode recieve centroid coords from master node
	call MPI_Recv(RRc, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
	!slavenode force calculation
	call potential(RRc, RRc, Upot, dRRc, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT)
	!slavenode send back centroid derivatives
	call MPI_Send(dRR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
	!slavenode send back centroid energy
	call MPI_Send(Upot, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
	!slavenode send back virial
	call MPI_Send(virial, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
	!slavenode send back dipole moments  
	!if (dip_out .or. TD_out) call MPI_Send(dip_momI, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
	!if (Edip_out) call MPI_Send(dip_momE, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
 endif


end subroutine contracted_forces

end module force_calc
