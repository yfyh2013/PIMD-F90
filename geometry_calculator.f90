!---------------------- calc geometry module --------------------------------------
! this is a self contained module for calculating the avg. geometry during a run
! it is designed to be callable by any MD program 
! 
! Copyright 2015 Daniel C. Elton
!-----------------------------------------------------------------------------

module geometry_calculator
implicit none
real, save :: sum_HH=0,   sum_OH=0,  sum_HOH=0 !distances
real, save :: sum_HH2=0, sum_OH2=0, sum_HOH2=0 !distances squared

real, save :: bsum_HH=0,   bsum_OH=0,  bsum_HOH=0 !bead-bead distances
real, save :: bsum_HH2=0, bsum_OH2=0, bsum_HOH2=0 !bead-bead distances2

integer, save :: num_total=0

 contains 
 
!-------------------------------------------------------------
!-- compute geometry ----------------------------------------
!-------------------------------------------------------------
subroutine calc_geometry(RR, RRt, Nbeads)
 use system_mod !source of Natoms, Nbeads, and boxi
 implicit none
 integer, intent(in) :: Nbeads
 double precision, dimension(3,Natoms,Nbeads), intent(in) :: RRt ! 3xNatomsxNbeads
 double precision, dimension(3,Natoms), intent(in) :: RR ! 3xNatoms
 double precision, dimension(3) :: ROH1, ROH2, RHH
 double precision :: d_HH, d_OH1, d_OH2, d_HOH
 integer ::  ih1, ih2, io, i, j


 do i = 1, Nwaters
	ih1 = 3*i
	ih2 = 3*i-1
	io = 3*i-2
	
	ROH1 = RR(:,io)  - RR(:,ih1)
	ROH1 = ROH1 - box*anint(ROH1*boxi)!PBC 

	ROH2 = RR(:,io)  - RR(:,ih2)
	ROH2 = ROH2 - box*anint(ROH2*boxi)!PBC

	RHH  = RR(:,ih2) - RR(:,ih1)
	RHH  = RHH - box*anint(RHH*boxi) !PBC
	
	d_HH  = dsqrt( dot_product(RHH,RHH) ) 
	d_OH1 = dsqrt( sum(ROH1**2) ) 
	d_OH2 = dsqrt( sum(ROH2**2) )  
	d_HOH = dacos( (d_OH1**2 + d_OH2**2 - d_HH**2)/(2d0*d_OH1*d_OH2) )
	
	sum_HH  = sum_HH  + d_HH    	
	sum_OH  = sum_OH  + d_OH1 + d_OH2
	sum_HOH = sum_HOH + d_HOH
	
	sum_HH2  = sum_HH2  + d_HH**2   	
	sum_OH2  = sum_OH2  + d_OH1**2 + d_OH2**2
	sum_HOH2 = sum_HOH2 + d_HOH**2

 enddo 
 
  do j = 1, Nbeads 
	do i = 1, Nwaters
		ih1 = 3*i
		ih2 = 3*i-1
		io  = 3*i-2
	
		ROH1 = RRt(:,io,j)  - RRt(:,ih1,j)
		ROH1 = ROH1 - box*anint(ROH1*boxi)!PBC 

		ROH2 = RRt(:,io,j)  - RRt(:,ih2,j)
		ROH2 = ROH2 - box*anint(ROH2*boxi)!PBC

		RHH  = RRt(:,ih2,j) - RRt(:,ih1,j)
		RHH  = RHH - box*anint(RHH*boxi) !PBC
		
		d_HH  = dsqrt( sum(RHH**2) ) 
		d_OH1 = dsqrt( sum(ROH1**2) ) 
		d_OH2 = dsqrt( sum(ROH2**2) )  
		d_HOH = dacos( (d_OH1**2 + d_OH2**2 - d_HH**2)/(2.0*d_OH1*d_OH2) )
	
		bsum_HH  = bsum_HH  + d_HH    	
		bsum_OH  = bsum_OH  + d_OH1 + d_OH2
		
		if(abs(d_HOH) .gt. 0) bsum_HOH = bsum_HOH + d_HOH
	
		bsum_HH2  = bsum_HH2  + d_HH**2  	
		bsum_OH2  = bsum_OH2  + d_OH1**2 + d_OH2**2
		bsum_HOH2 = bsum_HOH2 + d_HOH**2

	enddo 
 enddo
 
 num_total = num_total + Nwaters

endsubroutine calc_geometry



!-------------------------------------------------------------
!-- write out geometry --------------------------------------
!-------------------------------------------------------------
subroutine write_out_geometry(iun, Nbeads)
 implicit none 
 integer, intent(in) :: iun, Nbeads
 double precision :: bnum_total
 double precision :: avg_OH, avg_HH, avg_HOH
 double precision :: avg_OH2, avg_HH2, avg_HOH2, RMS_OH, RMS_HH, RMS_HOH
 double precision :: bavg_OH, bavg_HH, bavg_HOH
 double precision :: bavg_OH2,bavg_HH2,bavg_HOH2, bRMS_OH, bRMS_HH, bRMS_HOH

 
 bnum_total = num_total*Nbeads
 
 avg_HH   = sum_HH/num_total 
 avg_OH   = sum_OH/(2.0*num_total)
 avg_HOH  = (sum_HOH/num_total)*180d0/3.141592d0
 avg_HH2  = sum_HH2/num_total
 avg_OH2  = sum_OH2/(2.0*num_total)
 avg_HOH2 = (sum_HOH2/num_total)*(180d0/3.141592d0)**2

 RMS_HH  = dsqrt(avg_HH2 - avg_HH**2) 
 RMS_OH  = dsqrt(avg_OH2 - avg_OH**2) 
 RMS_HOH = dsqrt(avg_HOH2 - avg_HOH**2) 
 
 
 bavg_HH   = bsum_HH/bnum_total 
 bavg_OH   = bsum_OH/(2.0*bnum_total)
 bavg_HOH  = (bsum_HOH/bnum_total)*180d0/3.141592d0
 bavg_HH2  = bsum_HH2/bnum_total
 bavg_OH2  = bsum_OH2/(2.0*bnum_total)
 bavg_HOH2 = (bsum_HOH2/bnum_total)*(180d0/3.141592d0)**2

 bRMS_HH  = dsqrt(bavg_HH2  - bavg_HH**2) 
 bRMS_OH  = dsqrt(bavg_OH2  - bavg_OH**2) 
 bRMS_HOH = dsqrt(bavg_HOH2 - bavg_HOH**2) 


 
 write(iun,'(a)')         "#---------- h2o geometry report -------------------------"
 write(iun,'(a)')  "                                 centroid-centroid                | bead - bead  (Ang,deg) " 
 write(iun,'(a40,f12.3,a4,2f12.3,a4,f12.3)') "average HH distance : ", avg_HH, " +/-", RMS_HH, bavg_HH, " +/-", bRMS_HH  
 write(iun,'(a40,f12.3,a4,2f12.3,a4,f12.3)') "average OH distance : ", avg_OH, " +/-", RMS_OH, bavg_OH, " +/-", bRMS_OH  
 write(iun,'(a40,f12.3,a4,2f12.3,a4,f12.3)') "  average HOH angle : ", avg_HOH, " +/-", RMS_HOH, bavg_HOH, " +/-", bRMS_HOH  
#ifdef FC_HAVE_FLUSH
 call flush(iun) 
#endif

endsubroutine write_out_geometry

 
endmodule geometry_calculator