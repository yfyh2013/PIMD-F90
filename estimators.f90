!----------------------------------------------------------------------------------
!- PIMD Estimators 
!- Daniel C. Elton
!----------------------------------------------------------------------------------
module estimators 
 use consts
 use main_stuff

 contains
!----------------------------------------------------------------------------------
!- Quantum virial estimator for the energy and pressure (ref: Tuckerman, "Statistical Mechanics.." 2008 pg 485)
!- Inputs (virial, virialc, sys_temp, Upot are for the ENTIRE system) 
!----------------------------------------------------------------------------------
subroutine quantum_virial_estimators(RRt, virial, virialc, qEnergy, qPress, sys_temp, Upot) 
 Implicit none
 double precision, dimension(3,Natoms,Nbeads),intent(in) :: RRt 
 double precision, intent(in)      ::  sys_temp, Upot, virial, virialc
 double precision, intent(out)     :: qPress, qEnergy 
 double precision 		   :: qVirial, qVirial2, KE
 integer :: i, j, k 

 !Note on units : it is assumed that dRRt is in (in kcal/(mol*Ang))
 !and that Upot is in kcal/mol and that sys_temp is in Kelvin 
 !This pressure estimator assumes that the potential energy does not have volume dependence
 
 KE = 1.5*Natoms*kb*sys_temp*Nbeads!kinetic energy in kcal/mol

 !convert to kcal/(mole of mol H2O) by dividing by Nwaters
 qEnergy = ( KE  +  .5*(virial -  virialc ) +  Upot )/(Nwaters*Nbeads) 

 qPress  =  PRESSCON2*(1/(3*volume))*( 2*KE - virialc )/(Nwaters*Nbeads) !factor of 1/3 not in Tuckerman's book. (book is wrong!!)

end subroutine quantum_virial_estimators


!----------------------------------------------------------------------------------
!- Simple quantum estimators for energy & pressure 
!- these simple estimators take a long time to converge with 4+ beads
!----------------------------------------------------------------------------------
subroutine simple_quantum_estimators(RRt, virial, qEnergy, qPress, sys_temp, Upot) 
 use NormalModes !need MassScaleFactor
 implicit none
 double precision, dimension(3,Natoms,Nbeads),intent(in)  :: RRt  !coords in Ang
 double precision, intent(in)      :: sys_temp       !temp in Kelvin
 double precision, intent(in)      :: Upot 	!potential energy in kcal/mol
 double precision, intent(in)      :: virial 	!virial
 double precision, intent(out)     :: qEnergy  	!energy out in kcal/mol
 double precision, intent(out)     :: qPress  	!pressure out in bar
 double precision 		   :: KE, K0, mass
 integer :: i, j, k 

 !Note on units : it is assumed that dRRt is in
 !and that Upot is in kcal/mol and that sys_temp is in Kelvin 
 !The pressure estimator assumes that the potential energy does not have volume dependence

 K0 = 0
 do k = 1, Nbeads
	do j = 1, Natoms
		if (mod(j+2,3) .eq. 0) then
			mass = massO
		else 
			mass = massH
		endif
		do i = 1, 3
			if (k .eq. 1) then
				K0 = K0 + mass*( RRt(i,j,k) - RRt(i,j,Nbeads) )**2
			else
				K0 = K0 + mass*( RRt(i,j,k) - RRt(i,j,k-1) )**2
			endif
		enddo
	enddo
 enddo

 KE = 1.5*Natoms*Nbeads*kb*sys_temp*Nbeads !kinetic energy of the beads in kcal/mol
 K0 = .5*K0*(omegan**2)*MASSCONi   !quantum correction to kinetic energy - convert from Ang,amu,ps to kcal/mol

 qEnergy = (KE - K0 + Upot )/(Nwaters*Nbeads)   		    !kcal/(mole of mol)
 qPress  = PRESSCON2*(1/(3*volume))*(  2*(KE - K0) - virial )/Nwaters  !subtract virial to convert derivative to force

end subroutine simple_quantum_estimators



end module estimators 