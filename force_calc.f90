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
#ifdef parallel
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
#endif
    else
#ifdef parallel
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
#endif
	endif
	Call MPI_Barrier(MPI_COMM_WORLD, ierr)
 
 enddo! j = 0, Nbatches - 1 !batch index

if (pid .eq. 0) then
 Upot    = sum(Upott)    	!potential energy for the ENTIRE system
 virial  = sum(virialt)  	!virial for the ENTIRE system (kcal/mol)
 virialc = sum(virialct)/Nbeads !centroid virial for the ENTIRE system (kcal/mol)
 
 !computation of dipole moments for the SIESTA case 
 if (pot_model .eq. 6) then 
	call calc_dip_moments(dip_momIt, RRt, Nbeads, Nwaters,box)
 endif
 
 
 
endif 

end subroutine full_bead_forces


!--------------------------------------------------------------------------------------------
!------------- Contracted force calculation with intermolecular forces on monomer ----------
!-- This subroutine performs evaluates only the intramolecular (fast) forces on the beads and 
!-- evaluates the intermolecular (slow) forces on the centroid. It also uses a multiple timestep method
!-- (see Tuckerman, pg 118). The momentum and positions will be updated with the intramolecular forces
!-- every intra_timesteps times. For instance if the 'outer timestep' (normal timestep) is .5 ps and intra_timesteps = 5
!-- then the inner timestep is .1 ps
!--------------------------------------------------------------------------------------------
subroutine contracted_forces
 use main_stuff
 use NormalModes
 use pot_mod 
 use dans_timer
 Implicit None 
 double precision :: e1, virialmon, virialcmon, Umonomers, chgH1, chgH2, tmp
 double precision, dimension(3,Natoms,Nbeads)  ::  dRRfast
 double precision, dimension(3,Natoms) :: dRRc
 double precision, dimension(3,3)      :: dr1, r1
 double precision, dimension(3,3,3) :: dq3
 double precision, dimension(3)        :: roh1, roh2, rh1m, rh2m, rM, q3
 integer :: tintra, iM

 tmp = 0.5d0*gammaM/(1.d0-gammaM)
 
 if (t .eq. 1) dRRfast =0

if (pid .eq. 0) then
  call start_timer("MonomerPIMD")

  Umonomers = 0 
  virialmon = 0 
  virialcmon = 0 
	   
   !---  intramolecular (fast) forces -------------------------------------------------
  do tintra = 1, intra_timesteps

	!update momenta with fast forces
	PPt = PPt - MASSCON*dRRfast*delt2fast

	!update positions with fast forces
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
				RRt(:,3*i-2,k) = RRt(:,3*i-2,k) + imassO*PPt(:,3*i-2,k)*delt2fast
				RRt(:,3*i-1,k) = RRt(:,3*i-1,k) + imassH*PPt(:,3*i-1,k)*delt2fast
				RRt(:,3*i-0,k) = RRt(:,3*i-0,k) + imassH*PPt(:,3*i-0,k)*delt2fast
			enddo
		enddo			
	endif	  
	

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

			!if last timestep in loop update monomer energy
			!and calculate dipole moments using dip. mom. surface 
			if (tintra .eq. intra_timesteps) then
				Umonomers = Umonomers + e1

				!get charges
				call dms_nasa(r1, q3, dq3,box,boxi)

				q3 = q3*CHARGECON

				chgH1 = q3(2) + tmp*(q3(2)+q3(3))
				chgH2 = q3(3) + tmp*(q3(2)+q3(3))

				roh1 = RRt(:, iH1, j) - RRt(:, iO, j)
   				roh1 = roh1 - box*anint(roh1*boxi)!PBC
				roh2 = RRt(:, iH2, j) - RRt(:, iO, j)
  				roh2 = roh2 - box*anint(roh2*boxi)!PBC
  				
				!get m site position
				if ((pot_model .eq. 3) .or. (pot_model .eq. 4)) then 
					Rm = 0.5d0*gammaM*( roh1(:) + roh2(:) ) + RRt(:, iO, j)
	
					rh1m = RRt(:, iH1, j) - Rm
					rh1m = rh1m - box*anint(rh1m*boxi)!PBC
					rh2m = RRt(:, iH2, j) - Rm
					rh2m = rh2m - box*anint(rh2m*boxi)!PBC
  	
					dip_momIt(:,iw,j) = chgH1*rh1m + chgH2*rh2m
				else
					dip_momIt(:,iw,j) = chgH1*roh1 + chgH2*roh2
				endif
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

 call stop_timer("MonomerPIMD")
 
 !intermolecular force calculation
 call potential(RRc, RRc, Upot, dRRc, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT)

 !update dRRt
 do j = 1, Nbeads
	dRRt(:,:,j) = dRRc
 enddo

 !calculate virial
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
	enddo
	!add polarization dipoles
   	dip_momIt(:,:,j) = dip_momIt(:,:,j) + dip_momE
   	dip_momEt(:,:,j) = dip_momE
   enddo
 
   !update Upot, virial and virialc
   Upot    = Upot*Nbeads + Umonomers !potential energy for the ENTIRE system (all images)
   virial  = virialmon + virt(1,1) + virt(2,2) + virt(3,3) !virial for the ENTIRE system (all images)
   virialc = virialcmon/Nbeads + virialc

endif !(pid .eq. 0) then
 

end subroutine contracted_forces


!---------------------------------------------------------------------
!-----------------Call correct potential ----------------------------
!---------------------------------------------------------------------
subroutine potential(RR, RRc, Upot, dRR, virt, virialc, dip_momI, Edip_mom, chg, t, BAROSTAT)
use consts
use system_mod
use pot_mod
use fsiesta
use dans_timer
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR, RRc
double precision, intent(out) :: Upot, virialc
double precision, dimension(3, Natoms), intent(out) :: dRR
double precision, dimension(3, 3), intent(out) :: virt
double precision, dimension(Natoms), intent(out)  :: chg
double precision, dimension(3, NWaters), intent(out)  ::  dip_momI, Edip_mom
double precision, dimension(3) :: dip_mom
integer, intent(in) :: t
logical, intent(in) :: BAROSTAT
double precision, dimension(3,3) :: siesta_box
   
siesta_box = 0.0
siesta_box(1,1) = box(1)
siesta_box(2,2) = box(2)
siesta_box(3,3) = box(3) 

!All the stuff that depends on volume needs to be rescaled. 
 p4V = FOURPI/volume

!scale VdW long range correction due to box size change and add correction
!multiplying by a correction factor is slightly more efficient than recalculating the entire Uvdw_lrc term each timestep
 Uvdw_lrc = Uvdw_lrc0*(volume_init/volume)

!If the volume is changing than the Ewald k-vectors have to be reset every timestep
!if (BAROSTAT) call ewald_set(.false.)

if (pot_model==2 .or. pot_model==3) then
    call pot_ttm(RR, RRc, Upot, dRR, virt, virialc, dip_momI, Edip_mom, chg,t)
else if (pot_model==4 .or. pot_model==5) then
    call pot_spc(RR, Upot, dRR, virt, dip_momI, chg)
else if (pot_model==6) then 
    call start_timer("SIESTA")
  !  write(*,*) "RR=", RR
  ! write(*,*) "box=", siesta_box
    call siesta_forces( trim(sys_label), Natoms, RR, cell=siesta_box, energy=Upot, fa=dRR)
    call stop_timer("SIESTA")
    Upot = Upot*EVTOKCALPERMOLE
    dRR = -1d0*dRR*EVTOKCALPERMOLE    
endif

end subroutine potential


!---------------------------------------------------------------------
!- Calculate dipole moments using the TIP4P/2005 charges and m-site 
!- for the coordinates obtained from a SIESTA calculation 
!---------------------------------------------------------------------
subroutine calc_dip_moments(dip_momIt, RRt, Nbeads, Nwaters, box)
 use consts
 implicit none
 integer, intent(in) :: Nbeads, Nwaters
 double precision, dimension(3), intent(in) :: box
 double precision, dimension(3, Nwaters, Nbeads), intent(out) :: dip_momIt
 double precision, dimension(3, 3*Nwaters, Nbeads), intent(in) :: RRt
 real, dimension(3) :: r1, r2, r3, summ
 integer :: i, j, io, ih1, ih2
 double precision, parameter :: rOM = .1546
 double precision, parameter :: qH = .5564
 double precision, parameter :: qO = -1.1128

 do i = 1, Nbeads
	do j = 1, Nwaters
		io  = 3*j 
		ih1 = 3*j - 1 
		ih2 = 3*j - 2
		r1 = RRt(:,ih1,i) - RRt(:,io,i)
		r1 = r1 - box*anint(r1/box)!PBC
		r2 = RRt(:,ih2,i) - RRt(:,io,i)
		r2 = r2 - box*anint(r2/box)!PBC
		summ = r1 + r2!find vector to M-site
		r3=(summ/sqrt(dot_product(summ,summ)))*rOM
		dip_momIt(:,j,i) = (qH*r1 + qH*r2 + qO*r3)*a2m*e2coul/3.33564e-30!conv. to Debye
	enddo
 enddo
  
end subroutine calc_dip_moments



end module force_calc
