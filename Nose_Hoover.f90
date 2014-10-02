!----------------------------------------------------------------------------------!
!-------------- General Nose-Hoover subroutine to propagate one chain  ------------
!-------------- 2014 Dan Elton. This subroutine is self-contained.     ------------
!-------------- It propagates the chain a half time step, so it should ------------
!-------------- be called twice per timestep.                           -----------
!-------------- It could be modified to support chains of length 1      -----------
!-------------- but having only one chain is not recommended and thus prohibited --
!-------------- Ref.: Martyna, et. al. Mol. Phys. Vol. 87 5 1117 (1996) -----------
!----------------------------------------------------------------------------------!
subroutine Nose_Hoover(s, uk, chain_length, vxi, tau, delt2, L, temp)
 Implicit none 
 double precision, intent(out) 			:: s   	      !variable to couple to
 double precision, intent(inout) 		        :: uk 	      !kinetic energy
 integer, intent(in) 			        :: chain_length   !chain length
 double precision, dimension(chain_length), intent(inout) :: vxi !chain velocities
 double precision, intent(in) 		        :: tau 	      !kinetic energy
 double precision, intent(in) 		        :: delt2      !half timestep
 integer, intent(in) 			        :: L	      !# DOF being coupled to 
 double precision, intent(in) 		        :: temp	      !target temperature !DOUBLE PRECISION
 double precision 		 :: kT, delt4, delt8, M1, M2
 double precision, parameter  :: KB_amuA2ps2perK = .831447148d0
 double precision, parameter  :: hbar=6.35077993041d0 !hbar in amu*ang^2/ps
 integer 			  :: i

 if (chain_length .lt. 2) then 
	 write(*,*) "ERROR: Nose Hoover chain length cannot be less than two."
	 stop 
 endif 

 delt4 = delt2/2d0
 delt8 = delt2/4d0

 kT  = KB_amuA2ps2perK*temp 
 M1  = L*kT*tau**2 !tau = 1/omega
! M2  = /(kT*(32/hbar)**2) ! = kt/omegan**2 massive thermostating scheme recommended by Tuckerman book assuming 32 beads here. (pg 476)
 M2 = M1/L !Traditional thermostating 

 if (chain_length .gt. 2) then
 	i = chain_length
 	vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4

 	do i = chain_length - 1, 3
		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
		vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4
		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
 	enddo
 	vxi(2)  = vxi(2)*dExp(-vxi(3)*delt8)
 endif

 vxi(2) = vxi(2) + (M1*vxi(1)**2 - kT)*delt4/M2
 if (chain_length .gt. 2) vxi(2) = vxi(2)*dExp(-vxi(3)*delt8)
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)
 vxi(1) = vxi(1) + (2*uk - L*kT)*delt4/M1
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)

 s = dExp(-vxi(1)*delt2)

 !VV(i) = VV(i)*s
 uk    = uk*s*s
 
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)
 vxi(1) = vxi(1) + (2*uk - L*kT)*delt4/M1
 vxi(1) = vxi(1)*dExp(-vxi(2)*delt8)

 if (chain_length .gt. 2) vxi(2) = vxi(2)*dExp(-vxi(3)*delt8)
 vxi(2) =  vxi(2) + (M1*vxi(1)**2 - kT)*delt4/M2

 if (chain_length .gt. 2) then
	vxi(2)  = vxi(2)*dExp(-vxi(3)*delt8)

 	do i = 3, chain_length - 1
		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
		vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4
		vxi(i) = vxi(i)*dExp(-vxi(i+1)*delt8)
	enddo
	 i = chain_length
 	vxi(i) = vxi(i) + (vxi(i-1)**2 - kT/M2)*delt4
 endif 



end subroutine Nose_Hoover
