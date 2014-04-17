!edited by Dan Elton to improve speed 2/2013
module neigh_mod
type t_neigh
   integer :: N, N0
   integer, dimension(2000) :: j
   double precision, dimension(3, 2000) :: Rij
   double precision, dimension(2000) :: R2
end type t_neigh

contains
   subroutine find_neigh(iO, RR, neigh)
   use system_mod
   implicit none
   double precision, dimension(3, Natoms) :: RR
   type(t_neigh) :: neigh
   integer :: i, iO,  j, N, ix, iy, iz, nx
   double precision, dimension(3) :: Ri, Rij
   double precision :: R2

   N = 0
   Ri = RR(1:3, iO)

!old inefficient method for PBC has been commented out by D. Elton
!The old code first finds the neighbors which are in the cell and not 
!in periodic images and stores the number of them in N0. It does not appear that this variable is used at all. 
!The new method uses the ANINT() function to handle PBCs. With 128 molecules the code is about 10% faster with this change 
!To use this older scheme the following variables would need to be reintroduced: 
!nx = (1.d0*Rc) / box(1) + 1 !These are integers!
!ny = (1.d0*Rc) / box(2) + 1 !it will round to 1 assuming Rc = boxlength/2
!nz = (1.d0*Rc) / box(3) + 1 !they are used in find_neigh.f90
!repl_nx = nx 
!repl_ny = ny
!repl_nz = nz
!nidols = nx*ny*nz !not sure what this is for. it is always =1 

!      do j=iO+1, Nwaters
!        Rij = Ri - RR(1:3, j)
!	! Rij = Rij - box(1)*anint(Rij/box(1)) !PBC
!         R2 = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
!         if (R2<Rc2) then
!            N = N+1
!            neigh % Rij(1:3, N)  = Rij
!            neigh % R2(N)  = R2
!            neigh % j(N) = j
!         endif
!      enddo
!      neigh % N0 = N
!         
!      do ix=-repl_nx, repl_nx
!         do iy=-repl_ny, repl_ny
!            do iz=-repl_nz, repl_nz
!               if (ix/=0 .or. iy/=0 .or. iz/=0) then
!                  do j=iO, Nwaters
!                     Rij = Ri - RR(1:3, j) + box(1:3)*dble( (/ix, iy, iz/) )
!                     R2 = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)
!                     if (R2<Rc2) then
!                        N = N+1
!                        neigh % Rij(1:3, N)  = Rij
!                        neigh % R2(N)  = R2
!                       neigh % j(N) = j
!                     endif
!                  enddo ! do j=iO, Nwaters
!               endif
!            enddo
!         enddo
!      enddo
!----------------------------------          
          


!note that the loop only goes over the atoms with j>i0. This is because the neighbor routine is called within a loop in the force/potential calculation. 
	do j=iO+1, Nwaters
        	Rij = Ri - RR(1:3, j)
		Rij = Rij - box(1)*anint(Rij*boxi(1)) !PBC
		
	           !----------------------------------
		   !Alternate method of implementing PBCs. 
		   !With some compilers this may be faster than the ANINT() method 
		   !if the program is making slow ANINT() library calls
	           !It is left here for that reason. 
	           !Rij(1) = Ri(1) - RR(1, j)
                   !   if (Rij(1) .lt. 0.0d0) then 
                   !      nx = boxi(1)*Rij(1) - 0.5
                   !   else
                   !      nx = boxi(1)*Rij(1) + 0.5
                   !   endif
                   !  Rij(1) = Rij(1) - nx*box(1)
                   !----------------------------------
                   !  Rij(2) = Ri(2) - RR(2, j)
                   !   if (Rij(2) .lt. 0.0d0) then 
                   !      nx = boxi(2)*Rij(2) - 0.5
                   !   else
                   !      nx = boxi(2)*Rij(2) + 0.5
                   !   endif
                   !  Rij(2) = Rij(2) - nx*box(2)
                   !----------------------------------
                   !      Rij(3) = Ri(3) - RR(3, j)
                   !   if (Rij(2) .lt. 0.0d0) then 
                   !      nx = boxi(3)*Rij(3) - 0.5
                   !   else
                   !      nx = boxi(3)*Rij(3) + 0.5
                   !   endif
                   !  Rij(3) = Rij(3) - nx*box(3)
		   !----------------------------------

                 R2 = Rij(1)*Rij(1) + Rij(2)*Rij(2) + Rij(3)*Rij(3)

                 if (R2<Rc2) then
	                 N = N+1
                         neigh % Rij(1:3, N) = Rij
                         neigh % R2(N)  = R2
                         neigh % j(N) = j
                 endif
	enddo !do j=iO+1, Nwaters
   neigh % N = N

end subroutine find_neigh
end module neigh_mod
