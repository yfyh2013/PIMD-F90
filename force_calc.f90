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

end subroutine full_bead_forces


!--------------------------------------------------------------------------------------------
!----------------- Contracted force calculation with intermolecular forces on monomer-------
!--------------------------------------------------------------------------------------------
subroutine contracted_forces
 use main_stuff 
! use neigh_mod 
 Implicit None 
 double precision :: e1 !monomer energy 
 double precision, dimension(3,Natoms) ::  dRRc
 double precision, dimension(3,3)      :: temp1, temp2
 !type(t_neigh) :: neigh

 if (pid .eq. 0) then

 	!if we are using more than once processor, send the centroid coords out to be processed 
 	if (Nnodes .gt. 1) Call MPI_Send(RRc(:,:), counti, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, ierr)   

	!masternode calcuates the intramolecular forces, puts them in dRRt
	do j = 1, Nbeads
		do iw = 1, Nwaters
			 iO  = 3*iw-2
			 iH1 = 3*iw-1
			 iH2 = 3*iw-0
			temp1 = RRt(1:3, (/iO, iH1, iH2/), j)  
		  	
			call pot_nasa(temp1, temp2, e1, box, boxi)  

			dRRt(1:3, (/iO, iH1, iH2/), j) = temp2
		enddo
	enddo




	if (Nnodes .gt. 1) then 
		!masternode recieve forces on centroid		
		call MPI_Recv(dRRc, counti, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, status2, ierr)
		!masternode receive intermolecular energy 
		call MPI_Recv(Upot, 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)	
		!masternode recieve centroid based dipole moments		
		!if (dip_out .or. TD_out) call MPI_Recv(dip_momIt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
		!if (Edip_out) call MPI_Recv(dip_momEt(:,:,k), 3*Nwaters, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, status2, ierr)
	else 
		!if only one node then masternode calculates forces on centroid
		call potential(RRc, Upot, dRRc, virt, dip_momI, dip_momE, chg, t, BAROSTAT)
	endif
 
	!add centroid forces to each bead 
	do j = 1, Nbeads
		dRRt(:,:,j) = dRRt(:,:,j) + dRRc
	enddo


 !slave node stuff
 else
	!slavenode recieve centroid coords from master node
	call MPI_Recv(RR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status2, ierr)
	!slavenode force calculation
	call potential(RR, Upot, dRR, virt, dip_momI, dip_momE, chg, t, BAROSTAT)
	!slavenode send back centroid derivatives
	call MPI_Send(dRR, counti, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
	!slavenode send back centroid energy
	call MPI_Send(Upot, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr) 
	!slavenode send back dipole moments  
	!if (dip_out .or. TD_out) call MPI_Send(dip_momI, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
	!if (Edip_out) call MPI_Send(dip_momE, 3*Nwaters, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
 endif


end subroutine contracted_forces

end module force_calc
