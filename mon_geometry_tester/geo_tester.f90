!
!   geometry tester
!
program geo_tester
 use fsiesta
 implicit none
 integer, parameter :: num_SIESTA_nodes = 24
 character(100), parameter :: mon_siesta_name = "siesta3.2"
 character(100), parameter :: file_tag = "PBE_DZP"
 double precision, parameter :: minr=.2, rstep = .12 
 integer, parameter          :: nr = 18
 integer, parameter          :: mint=4,  maxt = 180, tstep = 4
 double precision, dimension(3,3) :: RR, dr1
 integer, parameter :: lunc = 22, lune = 23
 double precision :: r1, r2, r2x, r2y, e1
 double precision, parameter :: Deg2Rad = 3.141592654d0/180d0 
 integer :: ri1, ri2, theta
 
 !launch SIESTA
 if (num_SIESTA_nodes .eq. 1) then  
	call siesta_launch( trim(mon_siesta_name), "monomer") !launch serial SIESTA process
 elseif (num_SIESTA_nodes .gt. 1) then  
 	call siesta_launch( trim(mon_siesta_name), "monomer", nnodes=num_SIESTA_nodes ) !launch parallel SIESTA process  
 endif

 !Oxygen coord
 RR(:,1) = (/ 0.0, 0.0, 0.0 /) 

 !open output files
 open(lunc, file="out_"//trim(file_tag)//".ANI", status='unknown') 
 open(lune, file="out_"//trim(file_tag)//".energies", status='unknown') 
 
 write(*,*) "performing", int((nr**2)*(floor((maxt-mint)/real(tstep)))), " monomer calculations"
 
 !main loop 
 do ri1 = 0, nr-1
	r1 = minr + ri1*rstep
	do ri2 = 0, nr-1 
		r2 = minr + ri2*rstep
		do theta = mint, maxt, tstep

			!set up coords
			RR(:,2) = (/ dble(r1), 0.0d0, 0.0d0 /) !1st Hydrogen on x axis
			r2x = r2*dcos(dble(theta)*Deg2Rad)
			r2y = r2*dsin(dble(theta)*Deg2Rad)
			RR(:,3) = (/ dble(r2x), dble(r2y), 0.0d0 /)
 
			!write(*,*) r1, r2, theta
			
		        call siesta_monomer(RR, dr1, e1)

			!write out
			write(lunc,'(i10)') 3  
			write(lunc,*) " "
			write(lunc,'(a2,3(1x,f12.6))')'O ', RR(:, 1)
			write(lunc,'(a2,3(1x,f12.6))')'H ', RR(:, 2)
			write(lunc,'(a2,3(1x,f12.6))')'H ', RR(:, 3)
			
			write(lune,'(f15.8)')  e1
		enddo
	enddo
 enddo

 call siesta_quit('all')
 
end program geo_tester


!---------------------------------------------------------------------
!-----------------SIESTA monomer force calculation ------------------
!---------------------------------------------------------------------
subroutine siesta_monomer(r1, dr1, e1)
 use fsiesta
 implicit none 
 double precision, dimension(3, 3), intent(in) :: r1
 double precision, dimension(3, 3), intent(out) :: dr1
 double precision, dimension(3, 3) :: siesta_box
 double precision, intent(out) :: e1
 double precision, parameter :: EVTOKCALPERMOLE = 23.0600d0 !ev to kcal/mol 

 siesta_box = 0.0
 siesta_box(1,1) = 10.0
 siesta_box(2,2) = 10.0
 siesta_box(3,3) = 10.0
  
 call siesta_forces("monomer", 3, r1, cell=siesta_box, energy=e1, fa=dr1)

 e1 = e1*EVTOKCALPERMOLE
 dr1 = -1d0*dr1*EVTOKCALPERMOLE    
  
end subroutine siesta_monomer

