module main_stuff
use consts
use system_mod
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
integer :: num_timesteps, t, t_freq, tp_freq, td_freq, m, clock, eq_timesteps, TPoutStream, tt, tr 
logical :: dip_out, coord_out, TD_out, vel_out, TP_out, Edip_out
logical :: BAROSTAT, PEQUIL, BOXSIZEOUT, THERMOSTAT, GENVEL, INPVEL, PRINTFINVEL, BEADTHERMOSTAT, CALC_RADIUS_GYRATION

!N-H variables
double precision, save :: tau, tau_bead, s, sbead
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
double precision ::  omegan, iNbeads, setNMfreq
double precision :: radiusH, radiusO
integer :: Nnodes, pid, Nbeads, j, k, ierr, Nbatches, counti, bat
integer :: status2(MPI_STATUS_SIZE)
! timing variables
double precision :: seconds, secondsNM, secondsIO

contains

!----------------------------------------------------------------------------------!
!---------------- Read input file -------------------------------------------------
!----------------------------------------------------------------------------------!
subroutine read_input_file 

call getarg(1, finp)
narg = command_argument_count()
if (narg .eq. 0) then
	write(*,*) "ERROR: no input file argument."
	stop
endif 
open(11, file=finp, status='old')
read(11,*)read_method
read(11,*)fconfig
read(11,*)fsave
read(11,*)coord_out
read(11,*)vel_out
read(11,*)dip_out
read(11,*)Edip_out
read(11,*)TD_out
read(11,*)TP_out
read(11,*)CALC_RADIUS_GYRATION
read(11,*)pot_model 
read(11,*)Rc, rc1, eps_ewald
read(11,*)polar_maxiter, polar_sor, polar_eps, guess_initdip, print_dipiters
read(11,*)GENVEL
read(11,*)INPVEL
read(11,*)fvel
read(11,*)PRINTFINVEL
read(11,*)THERMOSTAT
read(11,*)BEADTHERMOSTAT
read(11,*)tau
read(11,*)tau_bead
read(11,*)global_chain_length
read(11,*)bead_chain_length
read(11,*)temp
read(11,*)BAROSTAT
read(11,*)tau_P
read(11,*)press 
read(11,*)BOXSIZEOUT		
read(11,*)PEQUIL
read(11,*)eq_timesteps
read(11,*)num_timesteps
read(11,*)delt
read(11,*)td_freq
read(11,*)tp_freq
read(11,*)t_freq
read(11,*)Nbeads
read(11,*)setNMfreq

close(11)  
end subroutine read_input_file 

!----------------------------------------------------------------------------------!
!---------------- Read in coordinate data to RRc ----------------------------------
!----------------------------------------------------------------------------------!
subroutine read_coords

if (BEADTHERMOSTAT .and. (Nbeads .eq. 1)) then
	write(*,*) "ERROR: Bead thermostat does not make much sense with 1 bead."
	write(*,*) "The dynamics will be unphysical since every atomic DOF will be thermostated. "
	stop
endif

if (BEADTHERMOSTAT .and. .not. (THERMOSTAT)) then
	write(*,*) "WARNING: running bead thermostating without a global thermostat is not recommended."
	write(*,*) "You may observe abnormally large temperature fluctuations in the system."
endif

if (.not. ( (box(1).eq.box(2)) .and. (box(2).eq.box(3)) ) ) then
	write(*,*) 'ERROR: program can only handle square boxes'
	stop 
endif   
if ( (INPVEL) .and. (GENVEL) ) then
	write(*,*) 'ERROR: You selected both to input velocities and generate velocities. Please choose one or the other'
	stop 
endif   
if ( Rc .gt. box(1)/2 ) then
	write(*,*) 'WARNING: cutoff radius is larger than half box size'
endif  
if (rc1 .lt. 0) then
 	write(*,*) "ERROR: start of shifted cutoff cannot be less than zero!!"
	stop
endif
if (Rc .lt. 0) then
 	write(*,*) "ERROR: Coloumb cutoff cannot be less than zero!!"
	stop
endif
if (rc1 .gt. Rc) then
 	write(*,*) "ERROR: start of shifted cutoff cannot be greater than Coloumb cuttoff!!"
	stop
endif

if (read_method==0) then
   do i=1, Nwaters
      iO = 3*i-2 
      read(10,*)ch2, RRc(1:3, iO)
   enddo
   do i=1, Nwaters
      ih1 = 3*i-1; ih2=3*i
      read(10,*)ch2, RRc(1:3, ih1)
      read(10,*)ch2, RRc(1:3, ih2)
   enddo
else if (read_method==1) then
   do i=1, Natoms
      read(10,*)ch2, RRc(1:3, i)
   enddo
endif

end subroutine read_coords



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
s = 1
sbead = 1
end subroutine initialize_variables


!----------------------------------------------------------------------------------!
!---------- Master node allocations -----------------------------------------------
!----------------------------------------------------------------------------------!
subroutine master_node_allocations
	if (Nnodes .lt. Nbeads + 1) then 
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
	allocate(massScaleFactor(Nbeads))
	dRRt = 0 
end subroutine master_node_allocations



!-----------------------------------------------------------------------------------------
!--------------Open write out files -----------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine open_files
if (coord_out) then
	open(20, file='out_'//TRIM(fsave)//'_coord.xyz', status='unknown')
endif
if (vel_out) then 	
	open(21, file='out_'//TRIM(fsave)//'_mom.dat', status='unknown')
endif
if  (dip_out) then
	open(22, file='out_'//TRIM(fsave)//'_dip.dat', status='unknown')
endif
if  (Edip_out) then
	open(26, file='out_'//TRIM(fsave)//'_Edip.dat', status='unknown')
endif
if (TD_out) then 
	open(23, file='out_'//TRIM(fsave)//'_tot_dip.dat', status='unknown')
endif
if (BOXSIZEOUT) then 
	open(25, file='out_'//TRIM(fsave)//'_box.dat', status='unknown')
endif
if (TP_out) then 
	TPoutStream = 24
	open(TPoutStream, file='out_'//TRIM(fsave)//'_TempPress.dat', status='unknown')
	write(*,*) "Supressing further terminal output to file"
else 
	 TPoutStream = 6
endif 

!write header
if (CALC_RADIUS_GYRATION) then 
	write(TPoutStream,'(a15,a10,a13,a11,a11,a23,a12,a7,a7)') "time (ps) ","temp (K) ", & 
	"press.(bar) ", "avg temp ","avg press ", "Pot. energy (kcal/mol)", &
	" Tot. energy", " Oxy r ", " Hyd r "
else
	write(TPoutStream,'(a15,a10,a13,a11,a11,a23,a12)') "time (ps) ","temp (K) ", & 
	"press.(bar) ", "avg temp ","avg press ", "Pot. energy (kcal/mol)"," Tot. energy"
endif

end subroutine open_files


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
!-------------- calling Nose-Hoover for the beads, coupling in real space (no longer used) 
!----------------------------------------------------------------------------------!
!subroutine bead_NH
! Implicit None
! double precision :: uk_bead
! Integer :: i, j, k 
! do i = 1, Natoms	
!	do j = 1, Nbeads
!		do k = 1, 3
!			if (mod(i+2,3) .eq. 0) then
!				uk_bead = ( 1d0/ (2d0*massO)  )*PPt(k,i,j)**2
!			else 
!				uk_bead = ( 1d0/ (2d0*massH)  )*PPt(k,i,j)**2
!			endif
!			call Nose_Hoover(sbead, uk_bead, bead_chain_length, vxi_beads(:,i,j,k), tau_bead, delt2, 1, Nbeads*temp)
!			PPt(k,i,j) = PPt(k,i,j)*sbead
!		enddo
!	enddo 
!enddo
!end subroutine bead_NH



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
!-------------- Write out info to file(s) -----------------------------------------
!----------------------------------------------------------------------------------!
subroutine write_out 
!for accuracy, pressure is computed at every timestep, since P fluctuations are large
	sys_press = (1/(3*volume))*( 2*uk - MASSCON*( virt(1,1)+virt(2,2)+virt(3,3) )   )!sign is important
	sys_press = PRESSCON*sys_press !convert to bar
	sum_press = sum_press + sys_press 

!print out temperature, pressure, average press & energies
if (mod(t,tp_freq) == 0) then
	sys_temp = TEMPFACTOR*uk/(Natoms)

	tt = tt + 1

	sum_temp      = sum_temp + sys_temp
	sum_energy    = sum_energy + Upot + uk*MASSCONi
	sum_energy2   = sum_energy2 + (Upot + uk*MASSCONi)**2
	sum_RMSenergy = sum_RMSenergy + (Upot + uk*MASSCONi - sum_energy/tt)**2

	!dip_mom(:) = sum(dip_momI(:, 1:Nwaters), dim=2) *DEBYE/CHARGECON
 	!write(*,'(3f12.4)') dip_mom

	if (CALC_RADIUS_GYRATION) then
		call calc_radius_of_gyration(RRt,RRc) 
		write(TPoutStream,'(1f10.4,6f12.2,2x,f6.4,2x,f6.4)') tr*delt, sys_temp ,sys_press, & 
		sum_temp/tt, sum_press/tr, Upot, Upot+uk*MASSCONi, radiusO, radiusH

	else	
		write(TPoutStream,'(1f10.4,6f12.2)') tr*delt, sys_temp ,sys_press, & 
			sum_temp/tt, sum_press/tr, Upot, Upot+uk*MASSCONi
	endif 

endif 
	tr = tr + 1

!reset averaging after equilbration ends
if  (t .eq. eq_timesteps) then
 	write(TPoutStream,*) "#----end of equilibration---restarting averages----------------------------------"
	!store average temp during equil and final energy after equil
	init_energy = Upot + uk*MASSCONi
	init_temp = sum_temp/tt

	!start new averaging of temp & energy
	tt     = 0
	tr     = 1
	sum_temp = 0 
	sum_press = 0 
	sum_energy = 0
	sum_energy2 = 0
	sum_RMSenergy = 0
endif 

!write out data during run 
if  (t .gt. eq_timesteps) then
	if (mod(t,t_freq)  == 0 ) then 
		if (coord_out) then
		     call save_XYZ(20, RRc, Upot, read_method, t, delt) 
	  	endif
	  	if (vel_out) then
	   	     call save_XYZ(21, PPc, Upot, read_method, t, delt) 
	   	endif
	   	if (dip_out) then
		     do iw=1,Nwaters
			 !calculate average dipole moment over all beads
			 do j = 1, 3
			 	dip_momI(j,iw) = sum(dip_momIt(j,iw,:))/Nbeads
			 enddo
			 !dsqrt(dot_product(dip_momI(:,iw), dip_momI(:, iw)))*DEBYE/CHARGECON
			 write(22,'(3(1x,f12.4))') dip_momI(:,iw)*DEBYE/CHARGECON 
 		     enddo
	   	endif
		if (Edip_out) then
		     do iw=1,Nwaters
		     	!calculate average dipole moment over all beads
 			do j = 1, 3
				dip_momE(j,iw) = sum(dip_momEt(j,iw,:))/Nbeads
			enddo
			write(26,'(3(1x,f12.4))') dip_momE(:,iw)*DEBYE/CHARGECON 
 		     enddo
	   	endif 
 
	endif
	if (mod(t,td_freq)  == 0 ) then 
   		if (TD_out) then
			!calculate average dipole moment over all beads
			do iw=1,Nwaters
				do j = 1, 3
					dip_momI(j,iw) = sum(dip_momIt(j,iw,:))/Nbeads
				enddo
			enddo
    	   		dip_mom(1:3) = sum(dip_momI(1:3, 1:Nwaters), dim=2) *DEBYE/CHARGECON
			write(23,'(3f12.4)') dip_mom
		endif
  	endif
endif

!box size output stuff 
if (BAROSTAT .and. BOXSIZEOUT) then 
        !for accuracy, box size computed at every timestep
        sum_box  = sum_box + box(1)
        sum_box2 = sum_box2 + box(1)**2
        if (mod(t,t_freq) == 0) write(25,*) sum_box/t
endif 

end subroutine write_out


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
!----------Print information about the run ----------------------------------------
!----------------------------------------------------------------------------------!
subroutine print_run

write(TPoutStream,'(a50, f10.3,a3)') "timestep = ", delt*1000, " fs"
if (THERMOSTAT) write(TPoutStream,'(a50, f10.3,a3)') "Nose-Hoover tau = ", tau, " ps"
if (.not. THERMOSTAT) write(TPoutStream,'(a50, a3)') "Nose-Hoover tau = ", "n/a"
if (BEADTHERMOSTAT) write(TPoutStream,'(a50, f10.3,a3)')  "Nose-Hoover tau for beads= ", tau_bead, " ps"
if (.not. BEADTHERMOSTAT) write(TPoutStream,'(a50, a3)') "Nose-Hoover tau for beads= ", "n/a"
if (BAROSTAT) write(TPoutStream,'(a50, f10.3,a3)') "Barostat tau = ", tau_p, " ps"
if (.not. BAROSTAT) write(TPoutStream,'(a50, a3)') "Barostat tau = ", "n/a"


write(TPoutStream,*) "#-----------  timing report ----------------------------"
write(TPoutStream,'(a50,i4,a,i4,a)') "Ran with ", Nbeads, " beads on ", Nnodes, " nodes"
write(TPoutStream,'(a50, i5, a7, i3, a9, i3, a8)') "Total elapsed time = ", int(real(seconds)/3600), " hours ",  & 
			int(mod(seconds,3600d0)/60), " minutes ", int(mod(seconds,60d0)), " seconds" 
write(TPoutStream,'(a50, i5, a7, i3, a9, i3, a9, i3, a3)') "Time spent on normal modes = ", int(real(secondsNM)/3600), " hours ", & 
			int(mod(secondsNM,3600d0)/60), " minutes ", int(mod(secondsNM,60d0)), " seconds ", &
			int(1000d0*mod(secondsNM,1d0)), " ms"
write(TPoutStream,'(a50, i5, a7, i3, a9, i3, a9, i3, a3)') "Time spent writing stuff out = ", int(real(secondsIO)/3600)," hours ",& 
			int(mod(secondsIO,3600d0)/60), " minutes ", int(mod(secondsIO,60d0)), " seconds ", int(1000d0*mod(secondsNM,1d0)), " ms"

write(TPoutStream,'(a50, f10.2)') "ps/hour = ", (  (num_timesteps + eq_timesteps)*delt/seconds  )*3600
write(TPoutStream,'(a50, f10.2)') "ps/day = ",  (  (num_timesteps + eq_timesteps)*delt/seconds  )*3600*24



write(TPoutStream,*) "#------------------------------------------------------"
avg_temp =  sum_temp/tt

write(TPoutStream,'(a50, 3f10.2)') "Average temperature during run (K) = ", avg_temp
write(TPoutStream,'(a50, 3f10.2)') "Average pressure during run (bar) = ", sum_press/tr
write(TPoutStream,'(a50, 3f10.2)') "Average total energy during run (kcal/mol) = ", sum_energy/tt
write(TPoutStream,'(a50, 3f10.2)') "estimated energy drift (kcal/mol/ps) = ",(Upot+uk*MASSCONi - init_energy)/(num_timesteps*delt)
write(TPoutStream,'(a50, 3f10.2)') "Temp drift (K/ps) = ", (avg_temp - init_temp)/(num_timesteps*delt)
write(TPoutStream,'(a50, 3f10.2)') "RMS energy fluctuation  (kcal/mol) = ", dsqrt( sum_RMSenergy/tt )

specific_heat = sqrt(   ( sum_energy2/tt - (sum_energy/tt)**2 )/( avg_temp * Kb )  )

write(TPoutStream,'(a50, f10.2)') "Specific heat C_V (only valid in NVT) (cal/g) = ", 1000*specific_heat/18

if (BAROSTAT .and. THERMOSTAT) then
        avg_box2 = sum_box2/(floor(real(tr/t_freq)))
        avg_box  = avg_box /(floor(real(tr/t_freq)))

        isotherm_compress = (avg_box2**3 - avg_box**3)*(10d-7)/(1.38*avg_temp*avg_box)
        write(TPoutStream,'(a50, f10.2)') "Isothermal compressibility (only valid in NPT)", isotherm_compress
else 
        write(TPoutStream,'(a50, a4)') "Isothermal compressibility (only valid in NPT)", " n/a"
endif

end subroutine print_run

!----------------------------------------------------------------------------------!
!----------Print information about the potential-----------------------------------
!----------------------------------------------------------------------------------!
subroutine print_pot(RR, Upot, dRR, virt, dip_momI, chg, TPoutStream) 
use system_mod
use consts
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR
double precision, intent(in) :: Upot
double precision, dimension(3, Natoms), intent(in) :: dRR
double precision, dimension(3,3), intent(in) :: virt
double precision, dimension(3, Nwaters), intent(in) :: dip_momI
double precision, dimension(3 ) :: dip_mom
double precision, dimension( Natoms ), intent(in) :: chg
integer, intent(in) :: TPoutStream
integer :: iw

dip_mom(1:3) = sum(dip_momI(1:3, 1:Nwaters), dim=2)

write( TPoutStream,'(a50,f14.6)')"Potential Energy (kcal/mol) = ", Upot
write( TPoutStream,'(a50,f14.6)')"Estimated Enthalpy of Vaporization &
(kJ/mol) = ", Upot*4.184 + 1000*8.3144 
write( TPoutStream,'(a50,f14.6)')"monomer Energy = ", Umon
write( TPoutStream,'(a50,f14.6)')"vdw Energy = ", Uvdw
write( TPoutStream,'(a50,f14.6)')"Long range vdw Energy = ", Uvdw_lrc
write( TPoutStream,'(a50,f14.6)')"electrostatic  Energy = ", Uelec
write( TPoutStream,'(a50,f14.6)')"induced Energy = ", Uind
write( TPoutStream,'(a50,3(1x,f12.3))')"Dipole moment [Debye] : ",dip_mom*DEBYE/CHARGECON
write( TPoutStream,'(a/3(10x,f14.5,1x))')"Virial tensor",virt(1:3,1:3)
!write(*,'(/a/)')"DERIVATIVES"
!do iw=1, Nwaters
!   write(*,'(a2,3x,3(f12.6,2x))')"O ",dRR(1:3, 3*iw-2)
!   write(*,'(a2,3x,3(f12.6,2x))')"H ",dRR(1:3, 3*iw-1)
!   write(*,'(a2,3x,3(f12.6,2x))')"H ",dRR(1:3, 3*iw-0)
!enddo
!write(*,'(/,"Dipole moments -X -Y -Z -MAG"/)')
!do iw=1, Nwaters
!   write(*,'(i4,2x,3(f12.6,1x),3x,f12.6)')iw, dip_momI(1:3,iw) / dble(nidols)  *DEBYE/CHARGECON, &
!    dsqrt(dot_product(dip_momI(1:3, iw), dip_momI(1:3, iw))) / dble(nidols)*DEBYE/CHARGECON
!enddo
end subroutine print_pot


!-----------------------------------------------------------------------------------------
!--------------Write out coords-----------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine save_XYZ(iun, RR, Upot, read_method, t, delt) 
use system_mod
use consts
implicit none
double precision, dimension(3, Natoms), intent(in) :: RR
double precision :: Upot,delt
integer ::  i, iO, ih1, ih2, t
integer :: iun, read_method

write(iun,'(i10)') Natoms !, angle
write(iun,'(f12.6,2x,f12.6,3(1x,f12.6))') t*delt, Upot,  box
if (read_method==0) then
   do i=1, Nwaters
      iO = 3*i-2
      write(iun,'(a2,3(1x,f12.6))')'O ',RR(1:3, iO)
   enddo
   do i=1, Nwaters
      ih1 = 3*i-1
      ih2 = 3*i
      write(iun,'(a2,3(1x,f12.6))')'H ',RR(1:3, ih1)
      write(iun,'(a2,3(1x,f12.6))')'H ',RR(1:3, ih2)
   enddo
else if (read_method==1) then
   do i=1, Nwaters
      iO = 3*i-2
      ih1 = 3*i-1
      ih2 = 3*i
      write(iun,'(a2,3(1x,f12.6))')'O ',RR(1:3, iO)
      write(iun,'(a2,3(1x,f12.6))')'H ',RR(1:3, ih1)
      write(iun,'(a2,3(1x,f12.6))')'H ',RR(1:3, ih2)
   enddo
endif
end subroutine save_XYZ




!----------------------------------------------------------------------------------!
!-------------- Initialize bead positions -----------------------------------------
!----------------------------------------------------------------------------------!
subroutine initialize_beads
use NormalModes
double precision, dimension(:,:), allocatable :: RRtemp
double precision :: avgrO, avgrH
allocate(RRtemp(3,Nbeads))

 CALL RANDOM_SEED(size = m) !get size of seed for the system
 ALLOCATE(seed(m))
 call system_clock(count=clock) 
	seed = clock + 357 * (/ (i - 1, i = 1, m) /)
 call random_seed(put = seed)  !put in the seed

!predict average radius of the ring polymer
avgrO  = 2*PI*hbar/(PI*Sqrt(24*massO*KB_amuA2ps2perK*temp))
avgrH  = 2*PI*hbar/(PI*Sqrt(24*1d0*KB_amuA2ps2perK*temp)) 
write(*,'(a58,f10.5,a4)') "Predicted average radius of the Oxygen ring polymer is   ", avgrO, " Ang"
write(*,'(a58,f10.5,a4)') "Predicted average radius of the Hydrogen ring polymer is ", avgrH, " Ang"

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
summom = 0
if (INPVEL) then
 write(*,*) "ERROR: Parallel version does not support inputting velocities at this time"
 write(*,*) "       Please select a different option!"
 stop 
endif
if (GENVEL) then

 if (allocated(seed) .eqv. .false.) then
 	CALL RANDOM_SEED(size = m) !get size of seed for the system
 	ALLOCATE(seed(m))
	call system_clock(count=clock) 
	seed = clock + 357 * (/ (i - 1, i = 1, m) /)
 	call random_seed(put = seed)  !put in the seed
 endif

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



!---------------------------------------------------------------------
!------------ Generate random number from Gaussian distribution -----
!------------ using Box_Muller sampling -----------------------------
!----- http://en.literateprograms.org/Box-Muller_transform_%28C%29 --
!---------------------------------------------------------------------
function rand_norm(std_dev) 
 Implicit None 
 double precision, intent(in) :: std_dev
 double precision  :: rand_norm, rand1, rand2, r

 r = 0
 do while ((r .eq. 0).or.(r .gt. 1)) 
 	call random_number(rand1)
 	call random_number(rand2)
 	rand1 = 2d0*rand1 - 1
 	rand2 = 2d0*rand2 - 1
	r = rand1*rand1 + rand2*rand2
 enddo 

 rand_norm = std_dev*rand1*Sqrt(-2d0*Log(r)/r)

end function rand_norm




end module main_stuff
