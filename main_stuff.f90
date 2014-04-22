module main_stuff
use consts
use system_mod !more global variables
use mpi
Implicit none
!----------------------------------------------------------------------------------!
!-------------- "global" variables used by main.f90 and the subroutines below ---
!----------------------------------------------------------------------------------!
!-old variables:---------------------------------------
double precision, dimension(:,:), allocatable :: RR, dRR
double precision, dimension(:), allocatable :: chg
double precision, dimension(3,3) :: virt
double precision, dimension(3) :: dip_mom
double precision, dimension(:,:), allocatable :: dip_momI, dip_momE
double precision :: Upot, sys_press
double precision :: tolg, tolx, stpmx
integer :: i, iw, iat,  iO, ih1, ih2, narg, ia, read_method
integer :: ix, iy, iz, nx, ny, nz
integer, external :: iargc
 character(len=2) :: ch2
 character(len=125) :: finp,fconfig,fvel,fsave 
!-new variables:--------------------------------------------
double precision, dimension(:,:), allocatable :: VV, dRRold, dRRnew
double precision, dimension(3) :: summom, sumvel
double precision :: delt, delt2,  uk,  imassO, imassH
double precision :: temp, sum_temp, sum_press,sys_temp, avg_vel, init_energy, sum_RMSenergy
double precision :: sum_energy, sum_energy2, specific_heat, avg_temp, init_temp
double precision :: avg_box, avg_box2, sum_box, sum_box2, isotherm_compress
integer, dimension(:), allocatable :: seed
 character(len=125) :: dip_file
 character(len=11)  :: bead_thermostat_type
integer :: num_timesteps, t, t_freq, tp_freq, td_freq, m, clock, eq_timesteps, TPoutStream, tt, tr 
logical :: dip_out, coord_out, TD_out, vel_out, TP_out, Edip_out
logical :: BAROSTAT, PEQUIL, BOXSIZEOUT, THERMOSTAT, GENVEL, INPVEL,PRINTFINALIMAGE
logical ::  BEADTHERMOSTAT, CENTROIDTHERMOSTAT, CALC_RADIUS_GYRATION

!N-H variables
double precision, save :: tau, tau_centroid, s, sbead
integer, save           :: global_chain_length, bead_chain_length 

!Nose-Hoover chain velocities for all beads are stored here
double precision, dimension(:,:,:,:), allocatable :: vxi_beads

!Nose-Hoover global chain velocities are stored here
double precision, dimension(:), allocatable     :: vxi_global

!-Berendsen thermostat variables
double precision :: tau_P, ref_P, press, CompFac, scale_factor

!- Variables for the paralleziation / PIMD
double precision, dimension(:,:,:), allocatable :: RRt, PPt, dip_momIt, dip_momEt, dRRt
double precision, dimension(:,:), allocatable :: RRc, PPc
double precision ::  omegan, kTN, iNbeads, setNMfreq
double precision :: radiusH, radiusO
integer :: Nnodes, Nbeads, pid, j, k, ierr, Nbatches, counti, bat
integer :: status2(MPI_STATUS_SIZE)

! timing variables
double precision :: seconds, secondsNM, secondsIO

contains



!----------------------------------------------------------------------------------!
!---------------- Initialize some variables ---------------------------------------
!----------------------------------------------------------------------------------!
subroutine initialize_variables

!--- Slave-node allocations ---- 
if (pid .ne. 0) then
	allocate(RR(3, Natoms))
	allocate(VV(3, Natoms))
	allocate(dRR(3, Natoms))
endif

!--- All-node allocations ------ 
allocate(dip_momI(3, Nwaters))
allocate(dip_momE(3, Nwaters))
allocate(chg (Natoms))
allocate(tx_dip(3,4*Nwaters, 4))

!These parameters are used later on --------
Rc2 = Rc * Rc
nx = (1.d0*Rc) / box(1) + 1 !These are integers!
ny = (1.d0*Rc) / box(2) + 1 !it will round to 1 assuming Rc = boxlength/2
nz = (1.d0*Rc) / box(3) + 1 !they are used in find_neigh.f90
repl_nx = nx 
repl_ny = ny
repl_nz = nz
nidols = nx*ny*nz !not sure what this is for. it is always =1 

Nwaters = Natoms/3
volume = box(1)*box(2)*box(3)
volume_init = volume
delt = delt/1000d0 !***CONVERT fs - > ps ***
delt2 = delt/2d0
boxi = 1.d0 / box
!inverse masses
imassO = DBLE(1/massO) 
imassH = DBLE(1/massH)
iNbeads = 1d0/DBLE(Nbeads)
CompFac = (.4477d-5*delt)/(tau_P) !Barostat var. (contains compressibility of H2O)
sum_temp = 0 
sum_press = 0 
sum_energy = 0
sum_energy2 = 0
sum_RMSenergy = 0
tt = 0 
tr = 0
counti = 3*Natoms
omegan = KB_amuA2ps2perK*temp*Nbeads/hbar
kTN = KB_amuA2ps2perK*temp*Nbeads
s = 1
sbead = 1

end subroutine initialize_variables


!----------------------------------------------------------------------------------!
!---------- Master node allocations -----------------------------------------------
!----------------------------------------------------------------------------------!
subroutine master_node_init
	use Langevin 
	use NormalModes

	!initialize random number generator
 	CALL RANDOM_SEED(size = m) !get size of seed for the system
 	ALLOCATE(seed(m))
	call system_clock(count=clock) 
	seed = clock + 357 * (/ (i - 1, i = 1, m) /)
 	call random_seed(put = seed)  !put in the seed
 	
	if (Nnodes .lt. Nbeads) then 
		if (.not. (mod(Nbeads,Nnodes) .eq. 0)) then
	           write(*,*) "ERROR: the number of beads must be a multiple of the number of nodes."
	           write(*,'(a,i4,a,i4,a)') "To run on ", Nnodes, " nodes I suggest using ", Nbeads - mod(Nbeads,Nnodes), " beads"
		stop
		endif
	else
		write(*,*) "WARNING : The number of processors is greater &
		than the number of beads!! \n Setting the number of beads to the number of processors (", Nnodes, ") "
		
		Nbeads = Nnodes
	endif
	
	write(*,'(a,i4,a,i4,a)') "Running with ", Nbeads, " beads on ", Nnodes, " nodes"

	!Master node allocations
	!only the master node (pid = 0) stores a fully copy of the
	! coords / vel for all beads and the centroid
	allocate(RRt(3, Natoms,Nbeads))
	allocate(PPt(3, Natoms,Nbeads))
	allocate(dRRt(3, Natoms,Nbeads))
	allocate(dip_momIt(3, Nwaters,Nbeads))
	allocate(dip_momEt(3, Nwaters,Nbeads))
	allocate(RRc(3, Natoms))
	allocate(PPc(3, Natoms))
	dRRt = 0 

	call InitNormalModes(Nbeads, omegan, delt, setNMfreq)

	if (THERMOSTAT)  then
		allocate(vxi_global(global_chain_length))
		vxi_global = 1 !set chain velocities to zero initially
	endif 
	if (BEADTHERMOSTAT)  then
		allocate(vxi_beads(bead_chain_length,natoms,Nbeads,3))
		vxi_beads = 0 !set chain velocities to zero initially
	endif
	if (bead_thermostat_type .eq. 'Langevin') call Init_Langevin_NM(delt2, CENTROIDTHERMOSTAT, tau_centroid, Nbeads, temp)


end subroutine master_node_init


!----------------------------------------------------------------------------------!
!---------- Periodic Boundary conditions for the beads ----------------------------
!----------------------------------------------------------------------------------!
!The PBCs follow the bead centroid - if the bead centroid crosses the edge of the box,
!then all the beads move with it. This means that at any time, some beads may lie 
!outside the box. The potential(RR, ...) subroutine must be able to handle situations
!where beads are outside the box
subroutine PBCs(RRt, RRc)
 Implicit none 
 double precision, dimension(3, Natoms,Nbeads), intent(inout) :: RRt
 double precision, dimension(3, Natoms), intent(inout) :: RRc
 double precision, dimension(3, Natoms) :: shifts 
 integer :: i

	!store the shifts here 
	shifts = box(1)*anint(RRc*boxi(1))

 	!Correct the centroids first
	RRc = RRc - shifts

	!move the beads 
	do i = 1,Nbeads	
		RRt(:,:,i) = RRt(:,:,i) - shifts
	enddo

end subroutine PBCs


!----------------------------------------------------------------------------------!
!-------------- calling thermostat for the beads ---------------------------------- 
!----------------------------------------------------------------------------------!
subroutine bead_thermostat
 use NormalModes
 use Langevin 
 Implicit None
 double precision, dimension(3,Nbeads) :: PPtr 
 double precision :: uk_bead, tau, imass
 Integer :: i, j, k 
! Calling Langevin thermostat 
 if (bead_thermostat_type .eq. 'Langevin') then
	call Langevin_NM(PPt, Nbeads)
 endif
! Nose-Hoover coupling in normal mode space ---------------------------------------
 if (bead_thermostat_type .eq. 'Nose-Hoover') then
   do i = 1,Natoms
	if (mod(i+2,3) .eq. 0) then
		imass = imassO 
	else 
		imass = imassH
	endif
	PPtr = real2NM(PPt(:,i,:),Nbeads) 
	if (CENTROIDTHERMOSTAT) then
		tau = tau_centroid
		do k = 1, 3
			uk_bead =  .5d0*imass*PPtr(k,1)**2
			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,1,k), tau, delt2, 1, Nbeads*temp)
			PPtr(k,1) = PPtr(k,1)*sbead
		enddo
	endif
	do j = 2, Nbeads 
		tau = 1d0/omegalist(j)
		do k = 1, 3
			uk_bead =  .5d0*imass*PPtr(k,j)**2
			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,j,k), tau, delt2, 1, Nbeads*temp)
			PPtr(k,j) = PPtr(k,j)*sbead
		enddo
	enddo 
	PPt(:,i,:) = NM2real(PPtr,Nbeads) 
   enddo
 endif


!----------------- Old Nose-Hoover scheme coupling in real space-------------------
! do i = 1, Natoms	
!	do j = 1, Nbeads
!		do k = 1, 3
!			if (mod(i+2,3) .eq. 0) then
!				uk_bead = ( 1d0/ (2d0*massO)  )*PPt(k,i,j)**2
!			else 
!				uk_bead = ( 1d0/ (2d0*massH)  )*PPt(k,i,j)**2
!			endif
!			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,j,k), tau_centroid, delt2, 1, Nbeads*temp)
!			PPt(k,i,j) = PPt(k,i,j)*sbead
!		enddo
!	enddo 
!enddo
end subroutine bead_thermostat



!----------------------------------------------------------------------------------!
!-------------- calculate average radius of gyration (for testing/debugging) ------
!----------------------------------------------------------------------------------!
subroutine calc_radius_of_gyration(RRt, RRc) 
 Implicit None
 double precision, dimension(3,Natoms,Nbeads),intent(in)  :: RRt
 double precision, dimension(3,Natoms),intent(in)         :: RRc
 integer    	    :: i, j, iH1, iH2, iO

 radiusH = 0d0 
 radiusO = 0d0

 do i = 1, Nwaters 
	iO = 3*i-2;  iH1 = 3*i-1;  iH2=3*i
	do j = 1, Nbeads
	radiusO = radiusO + sqrt( sum( (RRt(:,iO ,j) - RRc(:,iO ))**2, 1)  )	
	!write(*,*) i, RRt(1,iO,j), RRc(1,iO)
	radiusH = radiusH + sqrt( sum( (RRt(:,iH1,j) - RRc(:,iH1))**2, 1)  )
	radiusH = radiusH + sqrt( sum( (RRt(:,iH2,j) - RRc(:,iH2))**2, 1)  )
	enddo
 enddo
 
 radiusO = radiusO/dble(Nbeads*Nwaters)
 radiusH = radiusH/dble(2*Nbeads*Nwaters)

end subroutine calc_radius_of_gyration





!----------------------------------------------------------------------------------!
!--------------  Berendson Pressure Coupling (J. Chem. Phys. 81, 3684, 1984)-------
!----------------------------------------------------------------------------------!
subroutine Pcouple 
	scale_factor = ( 1 - CompFac*(press - sys_press)  )**.333333
	RR = RR*scale_factor
	box = box*scale_factor	
	volume = box(1)*box(2)*box(3)
	boxi = 1d0/box
end subroutine Pcouple 






!----------------------------------------------------------------------------------!
!-------------- Initialize bead positions -----------------------------------------
!----------------------------------------------------------------------------------!
subroutine initialize_beads
use NormalModes
double precision, dimension(:,:), allocatable :: RRtemp
double precision :: avgrO, avgrH
allocate(RRtemp(3,Nbeads))

!predict average radius of the ring polymer
if (Nbeads .gt. 1) then
	avgrO  = 2*PI*hbar/(PI*Sqrt(24*massO*KB_amuA2ps2perK*temp))
	avgrH  = 2*PI*hbar/(PI*Sqrt(24*1d0*KB_amuA2ps2perK*temp)) 
	write(*,'(a,f10.5,a4)') "Predicted average radius of (converged) Oxygen ring polymer is   ", avgrO, " Ang"
	write(*,'(a,f10.5,a4)') "Predicted average radius of (converged) Hydrogen ring polymer is ", avgrH, " Ang"
endif 

do i=1, Nwaters
		Call gen_rand_ring(RRtemp,massO,temp,Nbeads)	
		do j = 1,3
			RRt(j,3*i-2,:)  =  RRc(j,3*i-2) + RRtemp(j,:)
		enddo
		Call gen_rand_ring(RRtemp,massH,temp,Nbeads)	
		do j = 1,3
			RRt(j,3*i-1,:)  =  RRc(j,3*i-1) + RRtemp(j,:)
		enddo
		Call gen_rand_ring(RRtemp,massH,temp,Nbeads)	
		do j = 1,3
			RRt(j,3*i-0,:)  =  RRc(j,3*i-0) + RRtemp(j,:)
		enddo
enddo

call calc_radius_of_gyration(RRt,RRc) 

deallocate(RRtemp)

end subroutine initialize_beads


!---------------------------------------------------------------------------------
!---------------- Generate initial bead *momentum* ------------------------------
!---------------------------------------------------------------------------------
subroutine initialize_velocities
use math
summom = 0
if (INPVEL) then
 write(*,*) "ERROR: Parallel version does not support inputting velocities at this time"
 write(*,*) "       Please select a different option!"
 stop 
endif
if (GENVEL) then

!The program will generate Maxwell-Boltzmann velocities . With a small number of molecules,
!the initial temperature will never be exactly what is inputted
!(variances of +/- 10 K for 128 molecules) For this reason, the program regenerates
!the velocities until the temperature is within 1 K of what was specified. 
   do   
	summom = 0 
	do i=1, Nwaters
		do j = 1, Nbeads
			do k = 1, 3
				PPt(k,3*i-2,j) = rand_norm(Sqrt(KB_amuA2ps2perK*massO*Nbeads*temp))
				PPt(k,3*i-1,j) = rand_norm(Sqrt(KB_amuA2ps2perK*massH*Nbeads*temp))
				PPt(k,3*i,j)   = rand_norm(Sqrt(KB_amuA2ps2perK*massH*Nbeads*temp))
			enddo
			summom = summom + PPt(:,3*i-2,j) + PPt(:,3*i-1,j) + PPt(:,3*i-0,j)
		enddo
	enddo
	!remove center of momentum from system
	do i = 1, Natoms
		do j = 1, Nbeads
			PPt(:,i,j) = PPt(:,i,j) - summom/(Natoms*Nbeads)
	 	enddo
	enddo

	!calculate centroid velocities
	PPc = sum(PPt, 3)/Nbeads !centroid momenta

	!update kinetic energy 
	uk = 0
	do i = 1,Nwaters
		uk = uk + imassO*sum( PPc(:,3*i-2)**2 )
		uk = uk + imassH*sum( PPc(:,3*i-1)**2 )
		uk = uk + imassH*sum( PPc(:,3*i-0)**2 )
	enddo	
	uk = .5d0*uk 
	
	sys_temp = TEMPFACTOR*uk/Natoms

	!write(*,*) sys_temp

	if ( abs(sys_temp - temp) .lt. 1 ) then
		exit
	endif
 enddo

	
else if ( (.not. GENVEL) .and. (.not. INPVEL)) then 
  PPt = 0 
  PPc = 0 
  write(*,*) "Initial velocities set to zero." 
endif 
end subroutine initialize_velocities





end module main_stuff
