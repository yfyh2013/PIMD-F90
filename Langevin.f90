!---------------------------------------------------------------------!
!----------------- Langevin Thermostat for the beads -----------------
!----------------- This module is self contained ---------------------
!----------------------- D. Elton 4/22/14 ----------------------------
!------------------ Ref: JCP 133 124104 (2010) -----------------------
!---------------------------------------------------------------------!
module Langevin 
 Implicit none
 double precision, dimension(:), allocatable, save :: c1, c2
 
contains

!--------------- Langevin subroutine ---------------------------------
subroutine Langevin_NM(PPt, Nbeads)
 use NormalModes
 use math
 use consts
 use system_mod
 Implicit none
 double precision, dimension(3,Natoms,Nbeads), intent(inout) :: PPt 
 double precision, dimension(3,Nbeads)  :: PPtr
 double precision :: sqrtmass
 integer, intent(in) :: Nbeads
 integer :: i2, j, k
 
 do i2 = 1, Natoms

	!convert into normal mode space
	PPtr = real2NM(PPt(:,i2,:),Nbeads) 

	do j = 1, Nbeads 
		if (mod(i2+2,3) .eq. 0) then
			sqrtmass = sqrt(massO*MassScaleFactor(j)) 
		else 
			sqrtmass = sqrt(massH*MassScaleFactor(j))
		endif
		do k = 1, 3
			PPtr(k, j) = c1(j)*PPtr(k, j) + c2(j)*rand_norm(sqrtmass)
 		enddo 
	enddo
	
	!convert back into real space
	PPt(:,i2,:) = NM2real(PPtr,Nbeads) 

 enddo



end subroutine Langevin_NM


!--------------- Initialization --------------------------------------
subroutine Init_Langevin_NM(delt2, CENTROIDTHERMOSTAT, tau_bead, Nbeads, temp)
 use NormalModes
 use consts
 Implicit none
 double precision, intent(in) :: delt2, tau_bead, temp
 double precision :: gammak
 integer, intent(in) :: Nbeads
 integer :: i 
 logical, intent(in) :: CENTROIDTHERMOSTAT

 allocate(c1(Nbeads))
 allocate(c2(Nbeads))

 do i = 1, Nbeads
	if (i .eq. 1) then 
		if (CENTROIDTHERMOSTAT) then
			gammak = 1/tau_bead !zero mode case
 		else
			gammak = 0 !no centroid thermostat
		endif
 	else 
		gammak = 2*omegalist(i)! other modes
	endif

	c1(i) = dexp(-delt2*gammak)
	c2(i) = dsqrt(1 - c1(i)**2) 
	c2(i) = dsqrt(KB_amuA2ps2perK*temp*Nbeads)*c2(i) !include sqrt(kTN_b) factor here

 enddo
	
end subroutine Init_Langevin_NM

end module
