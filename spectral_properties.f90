!---------------------- calc infrared spectrum --------------------------------
! this is a module for calculating the infrared spectrum 
! and the velocity-velocity autocorrelation function 
! at the end of a run
!
! Copyright 2015 Daniel C. Elton
!-----------------------------------------------------------------------------
module spectral_properties
 use consts
 implicit none
 Integer 			 		 :: NumPointsOut = 250 !maximum number of points to use 
 double precision, parameter :: MaxFreqOut = 4500  !maximum frequency to plot (cm^-1)
 
 contains 
 
!-------------------------------------------------------------
!-- compute infrared spectrum and write out -----------------
!-------------------------------------------------------------
subroutine calc_infrared_spectrum(dip_moms,box,timestep,fsave,temp)
 use lun_management 
 implicit none
 real, dimension(:,:), intent(in) :: dip_moms 
 double precision, dimension(3), intent(in) :: box ! box size in Ang
 double precision, intent(in) :: timestep          ! times step in PS
 double precision, intent(in) :: temp	           ! temp in Kelvin
 character(len=*), intent(in) :: fsave             !file label 

 integer :: trun, tcor, ix, i, t, n, count, tread, lun
 double precision :: omega, IR, magn, avgMag, MinFreqOut, vol, PointsAvailable
 complex, dimension(:), allocatable :: aux1,vcross,ACF
 !double precision :: pi, A0,A1,A2,A3, window
 integer :: NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 tread = size(dip_moms,2)

 !find closest power of 2 less than tread
 trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  	 )

 allocate(aux1(0:2*trun-1)) !Auxiliary array, for fourier transform
 allocate(Vcross(0:2*trun-1)) ! Stores the cross correlation in reciprocal space
 allocate(ACF(0:2*trun-1))  ! stores  the autocorrelation Fucntion (real space)

 ACF=0 

 !find correlation function 
 do ix=1,3
	! Call to the direct FFT
    aux1(0:trun-1) = cmplx(dip_moms(ix,1:trun))
    n=size(aux1)
    call four1(aux1,n,-1)
    Vcross=aux1*conjg(aux1)          
    call four1(Vcross,n,1)
    ACF = ACF + Vcross
 enddo

 !norrmalization so ACF(1) = 1
 !ACF=ACF/real(ACF(1))

 ! Save the correlation function to file
 call io_open(lun, "out_"//trim(fsave)//"_dip_corr_function.dat")
 do i = 1, tread
 	write(lun,*) i*timestep, real(ACF(i))/real(ACF(1))
 enddo	
 call io_close(lun)

 !! Fourier transform the ACF
 aux1=cmplx(ACF)
 call four1(aux1,n,-1)
 n=size(aux1) !need this 
  
 !Save the IR spectrum in the file
 call io_open(lun, "out_"//trim(fsave)//"_IR_spectra.dat")

 MinFreqOut =  1d0/(timestep*ps2s*Cspeed*n)
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)
 
 if (PointsAvailable .lt. NumPointsOut) then 
	NumPointsOut = PointsAvailable
	NumAvgOver = 1
 else 
	NumAvgOver = floor(PointsAvailable/NumPointsOut)
 endif 
	
 !write(*,*) "#-------------infrared spectrum ------------------------"
 !write(*,*) "trun = ", trun
 !write(*,*) "n, tread = ",  n, tread
 !write(*,*) "MinFreqOut = ",  MinFreqOut
 !write(*,*) "MaxFreqOut = ",  MaxFreqOut
 !write(*,*) "points available = ", PointsAvailable
 !write(*,*) "Averaging over", NumAvgOver

 Do t = 0, NumPointsOut-1

  avgMag = 0 
  
  !block averaging
  avgMag = 0 
  do i = 1, NumAvgOver
	omega=( t*NumAvgOver+i )/(timestep*ps2s*n) !get freq in 1/s (Hz
	avgMag = avgMag + omega*tanh(hbarSI*omega*Cspeed/(Kb*2.0d0*Temp))*sqrt(real(aux1(t*NumAvgOver+i) )**2 &
					+ aimag( aux1(t*NumAvgOver+i) )**2)
  enddo
  avgMag = avgMag/real(NumAvgOver)

  ! Use the prefactor with harmonic quantum correction (See Ramirez paper)
  IR = (2d0*3.14159d0*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8) 

  IR = IR/100 !convert from 1/m to 1/cm

  IR = IR/2   !fudge factor (this is probably due to the fact that we have real data in, therefore a 2x redundancy when taking the FT. The FT part was not written by me and whoever wrote it clearly wasn't very careful about normalization. -D. Elton ) 

  omega = floor((t+.5)*numAvgOver)/(timestep*ps2s*n) !get central freq in 1/s (Hz
  omega = omega/Cspeed  	! convert frequency to cm-1
  
  if (.not. (IR .eq. IR)) IR = 0 !check for NaNs
  
  write(lun,*)  omega, IR
 EndDo

 call io_close(lun)

EndSubroutine calc_infrared_spectrum





!-------------------------------------------------------------------------
!- vel-vel ACF spectrum (aka "density of states") and write out 
!-------------------------------------------------------------------------
subroutine calc_DOS(Hvelocities, box, timestep, fsave, temp)
 use lun_management 
 implicit none
 real, dimension(:,:,:), intent(in) :: Hvelocities    ! Hvelocities stored in (3,Ntimesteps,Nhyd) array
 double precision, dimension(3), intent(in) :: box ! box size in Ang
 double precision, intent(in) :: timestep          ! times step in PS
 double precision, intent(in) :: temp	           ! temp in Kelvin
 character(len=*), intent(in) :: fsave             !file label 

 integer :: trun, tcor, ix, i, iH, t, n, count, tread, lun, Nhyd	
 double precision :: omega, IR, magn, avgMag, MinFreqOut, vol, PointsAvailable
 complex, dimension(:), allocatable :: aux1,vcross,ACF
 !double precision :: pi, A0,A1,A2,A3, window
 integer :: NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 tread = size(Hvelocities,2)
 Nhyd  = size(Hvelocities,3)

 !find closest power of 2 less than tread
 trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  	 )

 allocate(aux1(0:2*trun-1))  ! Auxiliary array, for fourier transform
 allocate(Vcross(0:2*trun-1))! Stores the cross correlation in reciprocal space
 allocate(ACF(0:2*trun-1))   ! stores the autocorrelation Function (real space)

 ACF=0 
 aux1 = 0 

 !find correlation function 
 do iH = 1, Nhyd
	do ix = 1, 3
		! Call to the direct FFT
		aux1(0:trun-1) = cmplx(Hvelocities(ix, 1:trun, iH))
		n=size(aux1)
		call four1(aux1,n,-1)
		Vcross=aux1*conjg(aux1)          
		call four1(Vcross,n,1)
		ACF = ACF + Vcross
	enddo
 enddo

 ! Save the correlation function to file
 call io_open(lun, "out_"//trim(fsave)//"_vel_vel_corr_function.dat")
 do i = 1, tread
 	write(lun,*) i*timestep, real(ACF(i))/real(ACF(1))
 enddo	
 call io_close(lun)
 !norrmalization so ACF(1) = 1
 !ACF=ACF/real(ACF(1))

 !! Fourier transform the ACF
 aux1=cmplx(ACF)
 call four1(aux1,n,-1)
 n=size(aux1) !need this 
  
 !Save the IR spectrum in the file
 call io_open(lun, "out_"//trim(fsave)//"_DOS.dat")

 MinFreqOut =  1d0/(timestep*ps2s*Cspeed*n)
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)
 
 if (PointsAvailable .lt. NumPointsOut) then 
	NumPointsOut = PointsAvailable
	NumAvgOver = 1
 else 
	NumAvgOver = floor(PointsAvailable/NumPointsOut)
 endif 
	

 Do t = 0, NumPointsOut-1

 avgMag = 0 
  
  !block averaging
  avgMag = 0 
  do i = 1, NumAvgOver
	omega=( t*NumAvgOver+i )/(timestep*ps2s*n) !get freq in 1/s (Hz
	avgMag = avgMag + real(aux1( t*NumAvgOver + i ))**2 + aimag(aux1( t*NumAvgOver + i ))**2
  enddo
  avgMag = avgMag/real(NumAvgOver)

  !prefactor
  IR = avgMag/Nhyd 

  IR = IR/2   !fudge factor (this is probably due to the fact that we have real data in, therefore a 2x redundancy when taking the FT. The FT part was not written by me and whoever wrote it clearly wasn't very careful about normalization. -D. Elton ) 

  omega = floor((t+.5)*numAvgOver)/(timestep*ps2s*n) !get central freq in 1/s (Hz
  omega = omega/Cspeed  	! convert frequency to cm-1
  
  if (.not. (IR .eq. IR)) IR = 0 !check for NaNs
  
  write(lun,*)  omega, IR
 EndDo

 call io_close(lun)

EndSubroutine calc_DOS
 
EndModule spectral_properties