!***
!***  allocate necessery arrays needed for the calculation of the potenial
!***
subroutine init_ttm21f(Nw)
use ttm21f_mod
implicit none
integer, intent(in) :: Nw
logical :: alloc


if (.not. allocated(RM)) then
   alloc = .true.
else if (Nw/=Nw_old) then
   alloc = .true.
   deallocate( RM )    ! temporary array keeping the coordinates of the Msite
   deallocate( dRM )   ! temporary array keeping the derivatives of the Msite
   deallocate( DDT )   ! dipole tensor
   deallocate( dip )   ! array containg the induced dipoles
   deallocate( pr_dip )  ! 
   deallocate( phi )    ! electrostatic potential
   deallocate( Efq )     ! electric field from charges
   deallocate( Efd )     ! electric field from dipoles
   deallocate( charge )  ! charges 
   deallocate( grdq )    ! dericatives of charge wrt monomer geometry
else
   alloc = .false.
endif

if (alloc) then
   Natsq = 3*Nw       ! # of atoms with charge   (including M-sites)
   Natsd = 3*Nw       ! # of atoms with dipole   (including M-sites)
   Nats  = 4*Nw       ! # of total atoms         (including M-sites)
   Nw_old = Nw
   fO = 1             ! index on the first oxygen
   lO = Nw            ! index on the last oxygen        
   fH = Nw+1          ! index on the first hydrogen
   lH = 3*Nw          ! index on the last hydrogen
   fM = 3*Nw+1        ! index on the first M-site
   lM = 4*Nw          ! index on the lst M-site
   fO3 = 3*fO-2  ! index on the x-comp. (out of the x,y,z) of the first oxygen
   lO3 = 3*lO    ! index on the x-comp. (out of the x,y,z) of the last oxygen
   fH3 = 3*fH-2  ! index on the x-comp. (out of the x,y,z) of the first hydrogen
   lH3 = 3*lH    ! index on the x-comp. (out of the x,y,z) of the last hydrogen
   fM3 = 3*fM-2  ! index on the x-comp. (out of the x,y,z) of the first M-site
   lM3 = 3*lM    ! index on the x-comp. (out of the x,y,z) of the last M-site

   allocate(charge(fO:lM))          
   allocate(RM (3, fM:lM))          
   allocate(dRM(3, fM:lM))         
   allocate(dip(3*Natsd))
   allocate(pr_dip(3*Natsd))
   allocate(phi(fH:lM))
   allocate(Efq(fO3:lM3))
   allocate(Efd(fO3:lH3))
   allocate( DDT(fO3:lH3, fO3:lH3) )
   allocate( grdq(Nw, 3, 3, 3) )
endif

end subroutine init_ttm21f
