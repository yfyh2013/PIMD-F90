!-----------------------------------------------
! Split line input file reader 
! Under preliminary construction
! Not yet implemented!! 
!------------------------------------------------------
module split_line_reader
 Implicit None
 MAX_LINE_LENGTH = 1000
 
 type KeyValue
	character(len=100) :: key
	double precision   :: value
	character(len=100) :: strvalue
 end type KeyValue

character(len=*) :: line

logical :: is_digit
integer :: i 
character(1) :: test_char

 contains 

subroutine PIMD_Input_File_read()


 coord_out   = .f.     ! centroid coordinates output 
 momenta_out = .f.     ! centroid velocities output  
 dip_out     = .f.     ! all dipoles output 
 Edip_out    = .f.     ! electronic (polarization) dipole output
 TD_out      = .f.     ! total dipole ouput 
 OUTPUTIMAGES = .f.    ! coordinates for all images output
 IMAGEDIPOLESOUT = .f. ! output dipoles for all PIMD images 
 BOXSIZEOUT  = .f.     ! box size running average output 
 TP_out      = .f.     ! Temp/Press (.t. = to file, .f. = to terminal)
 ENERGYOUT     = .t.   ! output total energy to seperate file
 HISTOUT       = .t.   ! output OH histograms to seperate file 
 
 CALC_RADIUS_GYRATION .t. ! output avg radius of gyration (as column in Temp/Press file)
 DIELECTRICOUT = .t.   ! output dielectric constant from running averages <M^2> - <M>^2 (as column in Temp/Press file)
 CHARGESOUT    = .f.   ! write out charges on atoms in seperate file?
 WRITECHECKPOINTS = .t.! configurations of all beads to seperate file 
 td_freq = 1		   ! Total dipole output frequency 
 tp_freq = 10 		   ! Temp/Press output frequency
 ti_freq = 2000		   ! all images output frequency 
 t_freq  = 10          ! Output frequency for everything else

 
 siesta_name = 'siesta' !name of siesta executable 
 sys_label   =     !SIESTA label of .fdf file
 num_SIESTA_nodes = -1 ! 
 
 Rc, rc1, eps_ewald
 polar_maxiter, polar_sor, polar_eps, guess_initdip, print_dipiters
 GENVEL = .true.       !generate velocities 
 THERMOSTAT = .f.      ! Global Nose-Hoover thermostat? 
 BEADTHERMOSTAT = .f.  ! Bead thermostat? 
 CENTROIDTHERMOSTAT = .f. !thermostat bead centroid? 
 bead_thermostat_type = 'none' !bead thermostat type 
 tau = .1               !
 tau_centroid = .1      ! 
 global_chain_length = 4  ! global Nose-Hoover chain length
 bead_chain_length = 4    ! bead Nose-Hoover chain length
 temp = 300 
 BAROSTAT = .f.       ! Berendson barostat ?
 tau_P = .2           ! tau for barostat (ps)
 press = 1            ! reference pressure (bar) 
 PEQUIL = .f.		  ! pressure coupling during equilibration only?

 setNMfreq = 0        ! frequency (inverse cm) to scale normal modes to (0 for none/RPMD)
 PIMD_type = 'full'   ! type of PIMD run (full", "contracted", "monomerPIMD") 
 intra_timesteps = 10 ! ratio of slow timestep / fast timestep for contraction scheme
 massO = 15.994		  ! mass of Oxygen (au)
 massH = 1.008		  ! mass of Hyrdrogen
 RESTART = .false.    ! restart? (append to previous files) 

 CALCGEOMETRY  = .true. !computes averge geometry of h2o molecules and outputs at end
 CALCDIFFUSION = .true. !computes diffusion constant of oxygen atoms
 read_method   = 1 !read_method(=0,1) (0 for OOOO....HHHHH and 1 for OHHOHHOHH...)
 CALCIRSPECTRA = .false. !store dipole moments and calculate IR spectra at end of run
 CALCDOS       = .true. 

 SIESTA_MON_CALC = .false. !SIESTA monomer calculation 

 


MD OPTIONS 
 'siesta','h2o', 1 			    ! sys_label for SIESTA runs, # SIESTA processors 
7d0 5d0 1.d-6           ! realspace cutoff (Ang), start of VdW force switch (Ang) (TTM3F only), eps for aewald
15  0.7  1.d-2 .f.  .f. ! params for dipole iterative procedure: Max iterations, pol factor, accuracy, guess initial?, debug output
.t.					 	! generate initial velocities? 
.f.    			 		! Global thermostat? t/f
.f.      				! Thermostat non-centroid modes ? 
.f.			 			! Thermostat centroid mode ?
'Nose-Hoover'     			! type of bead thermostat - 'Nose-Hoover', 'Langevin' or 'none'
.01			 			! tau for global thermostat (ps)
.1			 			! tau for centroid thermostat (ps)
 
300	     ! temperature (Kelvin)
 

end subroutine PIMD_Input_File_read 
 
 
 
 

!-----------------------------------------------------------------
!--  Process an imput file and place into database
!-----------------------------------------------------------------
Subroutine split_line_process_file(filename) 
 Implicit None
 character(len=200), intent(in) :: filename 
 character(MAX_LINE_LENGTH) :: aline
 integer :: lun=20
 
 !!One line function declarations 
 is_digit(i) = (i .ge. 48) .and. (i .le. 57)
 is_upper(i) = (i .ge. 65) .and. (i .le. 90)
 is_lower(i) = (i .ge. 97) .and. (i .le. 122)
 is_alpha(i) = is_upper(i) .or. is_lower(i)
 is_alnum(i) = is_digit(i) .or. is_alpha(i)

 is_comment(i) = (i .eq. 33) .or. (i .eq. 35) .or. (i .eq. 59) !Comment characters:  !  #  ; 
 is_delstr(i)  = (i .eq. 34) .or. (i .eq. 39) .or. (i .eq. 96) !String delimiters: "  '  `
 
 open(lun, file=filename, status='old')

 !main loop
 do 
	  aline  !read line

	aline = trim(aline) 
		
	if is_blank(aline) goto 1000
	if is_comment(aline(1)) goto 1000
	
	for i = 1, len(aline) 
		if 
	
	
 
1000  continue

 enddo
 

End Subroutine split_line_process_file




!-----------------------------------------------------------------
!--  Make a string lower case
!-----------------------------------------------------------------
Function make_lower(astr) result(astr)
 Implicit None
 character(len=*), intent(inout) :: astr
 integer i,j

 is_upper(i) = (i .ge. 65) .and. (i .le. 90)

 do i = 1, len(astr)
	j = ichar(astr(i)) 
	if is_upper(j) then
		astr(i) = achar(j + 32)
	endif
 enddo

End Function make_lower

      
end program test_fun
