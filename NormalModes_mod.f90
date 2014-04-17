!---------------------------------------------------------------------!
!-----------------Normal mode calculation module ---------------------
!-----------------This module is self contained ----------------------
!---------------------------------------------------------------------!
module NormalModes
 use consts 
 Implicit none
 double precision, dimension(:), allocatable, save   :: A11, A12, A21, A22
 double precision, dimension(:), allocatable, save   :: freqlist
 double precision, dimension(:,:), allocatable, save :: C

contains
!---------------------------------------------------------------------
!----------- Evolve one ring a time period delta --------------------
!---------------------------------------------------------------------
subroutine EvolveRing(RR, PP, Nbeads, mass)
 Implicit none
 double precision, dimension(3,Nbeads), intent(inout)  :: RR, PP
 double precision, dimension(3,Nbeads)  :: RRtr, PPtr
 double precision, dimension(3) :: temp2
 double precision, intent(in)   :: mass
 integer, intent(in) :: Nbeads
 integer              :: i, j, k

!Transform into normal mode representation 
call real2NM(RR,RRtr,Nbeads)
call real2NM(PP,PPtr,Nbeads) 

!Propagate Normal Modes
!j = the normal mode index. The number of normal modes always equals the number of beads
do j = 1, Nbeads
	!The normal mode evolution is simply the evolution of a harmonic oscillator
	!The matrix elements of the propagation matrix are stored in Aij(j) (minus the mass terms)
	!this may not be the best way to do this , but ... I wanted to avoid matrix multiplication	
	temp2     =  A11(j)*RRtr(:,j)      + A12(j)*PPtr(:,j)/mass
	PPtr(:,j) =  A21(j)*mass*RRtr(:,j) + A22(j)*PPtr(:,j)  
 	RRtr(:,j) = temp2
enddo

!Transform normal modes back into real space
call NM2real(RR,RRtr,Nbeads)
call NM2real(PP,PPtr,Nbeads)

PPtr = PPtr

end subroutine EvolveRing

!--------------------------------------------------------------------------
!------------ Initialize Normal Mode matrices ----------------------------
!--------------------------------------------------------------------------
subroutine InitNormalModes(Nbeads,omegan,delta,setNMfreq)
 use consts
 Implicit None
 double precision, intent(in) :: omegan, delta, setNMfreq
 double precision              :: omega
 integer, intent(in) 	 :: Nbeads 
 integer 			 :: i, j, k, l

 allocate(A11(Nbeads)) 
 allocate(A12(Nbeads)) 
 allocate(A21(Nbeads)) 
 allocate(A22(Nbeads)) 
 allocate(freqlist(Nbeads))
 allocate(C(Nbeads,Nbeads))

 if (mod(Nbeads,2) .eq. 0) then
 	l = (Nbeads-2)/2  !even Nbeads case
 else
 	l = (Nbeads-1)/2  !odd Nbeads case
 endif


!---------------- construct matrix -------------------------------------------- 
!The modes are stored in the following order: 
!--negative modes 
!--positive modes
!--zero mode
!--extra mode for when Nbeads is even

 do i = 1, Nbeads  !i = column index
 	do k = 1, l
		C(k  ,i) = Sqrt(2d0/Nbeads)*Cos(k*i*2d0*PI/Nbeads) !-k modes 
		C(k+l,i) = Sqrt(2d0/Nbeads)*Sin(k*i*2d0*PI/Nbeads) !+k modes
	enddo

		C(2*l+1,i) = Sqrt(1d0/Nbeads) !zero mode
	
	!extra mode for even Nbeads case only
	if (mod(Nbeads,2) .eq. 0) then
		C(2*l+2,i) = Sqrt(1d0/Nbeads)*((-1)**i)
	endif 
 enddo


!-------------- RPMD frequencies (setNMfreq = 0) --------------------------------------
!-------------- This case is always done for reference --------------------------------
!negative k modes
 do k = 1,l
	omega = 2*omegan*Sin(PI*( abs(k)/real(Nbeads)) )  
 	!write(*,'(a5,i3,a3,f8.3,a8,f10.2,a6)') "freq ", -k, " = ", omega/(2*PI), " 1/ps = ", 33.33333d0*omega/(2*PI), " cm^-1"	
    	freqlist(k) = omega
	A11(k) = Cos(omega*delta)
	A12(k) = (1/omega)*Sin(omega*delta)
	A21(k) = -omega*Sin(omega*delta)
	A22(k) = Cos(omega*delta)	
 enddo

!positive k modes
 do k = 1,l
	omega = 2*omegan*Sin(PI*( abs(k)/real(Nbeads)) )  
  	freqlist(k+l) = omega
	A11(k+l) = Cos(omega*delta)
	A12(k+l) = (1/omega)*Sin(omega*delta)
	A21(k+l) = -omega*Sin(omega*delta)
	A22(k+l) = Cos(omega*delta)	
 enddo
!correct omega=0 mode for the omega-> 0 limit
	freqlist(2*l+1) = 0d0
	A11(2*l+1) = 1d0
	A12(2*l+1) = delta
	A21(2*l+1) = 0d0
	A22(2*l+1) = 1d0	

!extra mode for even Nbeads case only
 if (mod(Nbeads,2) .eq. 0) then 
	k = l+1
	omega = 2*omegan*Sin(PI*( abs(k)/real(Nbeads)) )  
 	!write(*,'(a5,i3,a3,f8.3,a8,f10.2,a6)') "freq ", k, " = ", omega/(2*PI), " 1/ps = ", 33.33333d0*omega/(2*PI), " cm^-1"	
    	freqlist(2*l + 2) = omega
	A11(2*l+2) = Cos(omega*delta)
	A12(2*l+2) = (1/omega)*Sin(omega*delta)
	A21(2*l+2) = -omega*Sin(omega*delta)
	A22(2*l+2) = Cos(omega*delta)	
 endif 

if (setNMfreq .eq. 0) then

	massScaleFactor = 1

 	write(*,*) "Running usuing RPMD. All beads have physical mass."
 	write(*,*) "Normal mode frequencies: (cm^-1)"
 	do k = 1,l
		omega= freqlist(k) 
		write(*,'(f10.2)')  33.33333d0*omega/(2*PI) 
	enddo
 	if (mod(Nbeads,2) .eq. 0) then 
	    	omega = freqlist(2*l + 2)
	 	write(*,'(f10.2)')  33.33333d0*omega/(2*PI) 
 	endif

endif 


!------------------------- ACMD / CMD / PIMD case (adiabaticity > 1) ------------------------------------
if (.not. (setNMfreq .eq. 0)) then

!figure out masses 
	omega = (2*PI)*setNMfreq/33.33333d0 !conv. cm-1 -> 1/ps
	
	do i = 1,Nbeads
		if (freqlist(i) == 0) then 
			massScaleFactor(i) = 1d0
		else
			massScaleFactor(i) = (freqlist(i)/omega)**2
		endif
		write(*,*) "mass scale factor", i, " = ", massScaleFactor(i)
	enddo

	write(*,*) "Adiabaticity = ", omega/omegan
 	write(*,'(a,f8.3,a8,f10.2,a6)') "All normal modes scaled to ", omega/(2*PI), " 1/ps = ", 33.33333d0*omega/(2*PI), " cm^-1"
	write(*,'(a,f8.2,a)') "Timestep should probably not be larger than ", 	(((2*PI)/omega)/4d0)*1000, " fs"
	
  	freqlist(1:2*l) = omega
	A11(1:2*l) = Cos(omega*delta)
	A12(1:2*l) = (1/omega)*Sin(omega*delta)
	A21(1:2*l) = -omega*Sin(omega*delta)
	A22(1:2*l) = Cos(omega*delta)

	!omega=0 mode for the omega-> 0 limit
	freqlist(2*l+1) = 0d0
	A11(2*l+1) = 1d0
	A12(2*l+1) = delta
	A21(2*l+1) = 0d0
	A22(2*l+1) = 1d0

 !extra mode for even Nbeads case only
 if (mod(Nbeads,2) .eq. 0) then 
    	freqlist(2*l + 2) = omega
	A11(2*l+2) = Cos(omega*delta)
	A12(2*l+2) = (1/omega)*Sin(omega*delta)
	A21(2*l+2) = -omega*Sin(omega*delta)
	A22(2*l+2) = Cos(omega*delta)	
 endif

endif


end subroutine InitNormalModes

!-----------------------------------------------------------------------------------------
!---------------- generage one ring polymer sampled from the-----------------------------
!---------------- free ring distribution at temperature T -------------------------------
!-----------------------------------------------------------------------------------------
subroutine gen_rand_ring(RR,mass,temp,Nbeads)
 Implicit None
 double precision, parameter :: hbar=6.35077993041d0 !hbar in amu*ang^2/ps
 double precision, parameter :: KB_amuA2ps2perK = .831447148d0
 double precision, dimension(3,Nbeads), intent(out) :: RR
 double precision, dimension(3,Nbeads)	         :: RRtr
 double precision, intent(in) :: mass, temp
 double precision 		 :: std_dev, omega, omegan
 integer 			 		 :: j, k, l
 integer, intent(in) 		 :: Nbeads

 omegan = KB_amuA2ps2perK*temp*real(Nbeads)/hbar

 do j = 1,Nbeads
		omega = freqlist(j)
		if (omega .eq. 0d0) then 
			RRtr(:,j) = 0d0
		else 
  			std_dev = Sqrt((omegan*hbar)/(mass*omega**2)) !removed massscalefactor
			RRtr(:,j) = (/ rand_norm(std_dev), rand_norm(std_dev), rand_norm(std_dev) /)	
		endif	
 enddo

 RR = 0
 do j = 1,Nbeads 
	do k = 1, Nbeads
		RR(1,j) = RR(1,j) + C(k,j)*RRtr(1,k)
		RR(2,j) = RR(2,j) + C(k,j)*RRtr(2,k)
		RR(3,j) = RR(3,j) + C(k,j)*RRtr(3,k)
	enddo
 enddo

end subroutine gen_rand_ring


!---------------------------------------------------------------------
!------------ generate initial coordinates for a pure normal mode ---
!------------ (extra code for testing / visualization) --------------
!---------------------------------------------------------------------
!
!subroutine MakeNormalMode(RR, i, Nbeads) 
! Implicit none
! double precision, dimension(3,Nbeads), intent(out) :: RR
! double precision, dimension(3,Nbeads)		 :: RRtr
! integer, intent(in)	::  i, Nbeads
! integer 		::  j,k
!
! RRtr = 0
! RRtr(:,i) = 6d0
!
! do j = 1,Nbeads
!	do k = 1, Nbeads
!		RR(1,j) = RR(1,j) + C(k,j)*RRtr(1,k)
!		RR(2,j) = RR(2,j) + C(k,j)*RRtr(2,k)
!		RR(3,j) = RR(3,j) + C(k,j)*RRtr(3,k)
!	enddo
! enddo
!
!end subroutine MakeNormalMode


function real2NM(AA, AAtr, Nbeads) 
 Implicit none
 double precision, dimension(3,Nbeads), intent(in)  ::AA
 double precision, dimension(3,Nbeads), intent(out)  :: AAtr
 integer, intent(in) :: Nbeads
 integer :: j, k 
	AAtr = 0 
	do j = 1,Nbeads
		do k = 1, Nbeads
			AAtr(1,j) = AAtr(1,j) + C(j,k)*AA(1,k)
			AAtr(2,j) = AAtr(2,j) + C(j,k)*AA(2,k)
			AAtr(3,j) = AAtr(3,j) + C(j,k)*AA(3,k)
		enddo
	enddo
return AAtr
end subroutine real2NM

subroutine NM2real(AA, AAtr, Nbeads) 
 Implicit none
 double precision, dimension(3,Nbeads), intent(in)  :: AAtr
 double precision, dimension(3,Nbeads), intent(out)  :: AA
 integer, intent(in) :: Nbeads
 integer :: j, k 
do j = 1,Nbeads
	AA = 0 
	do k = 1, Nbeads
		AA(1,j) = AA(1,j) + C(k,j)*AAtr(1,k)
		AA(2,j) = AA(2,j) + C(k,j)*AAtr(2,k)
		AA(3,j) = AA(3,j) + C(k,j)*AAtr(3,k)
	enddo
enddo
return AA
end subroutine NM2real


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




end module
