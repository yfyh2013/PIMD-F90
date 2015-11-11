!---------------------- calc infrared spectrum --------------------------------
! this is a module for calculating the infrared spectrum 
! at the end of a run
!
! 
! Copyright 2015 Daniel C. Elton
!-----------------------------------------------------------------------------
module infrared
implicit none

 contains 
 
!-------------------------------------------------------------
!-- compute infrared spectrum and write out -----------------
!-------------------------------------------------------------
subroutine calc_infrared_spectrum(dip_moms,box,timestep,fsave,temp)
 use lun_management 
 implicit none
 double precision, dimension(:,:), intent(in) :: dip_moms 
 double precision, dimension(3), intent(in) :: box ! box size in Ang
 double precision, intent(in) :: timestep          ! times step in PS
 double precision, intent(in) :: temp	           ! temp in Kelvin
 character(len=*), intent(in) :: fsave             !file label 

 integer :: trun, tcor, ix, i, t, n, count, tread, lun_IR
 double precision :: omega, IR, magn, avgMag, MinFreqOut, MaxFreqOut, vol, PointsAvailable
 complex, dimension(:), allocatable :: aux1,vcross,ACF
 double precision, parameter :: Cspeed=3.00d10 ! cm/s
 double precision, parameter :: Kb=1.38d-23 !J/K
 double precision, parameter :: hbar=6.626d-34 ! J*s
 double precision, parameter :: ps2s=1d-12 ! 1fs in s
 double precision, parameter :: Debye2SI = 3.33564d-30
 !double precision :: pi, A0,A1,A2,A3, window
 integer :: NumAvgOver, NumPointsOut

 NumPointsOut = 200
 MaxFreqOut = 4500 !maximum frequency to plot (cm^-1)

 tread = size(dip_moms,2)

 !find closest power of 2 less than tread
 trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  )

 allocate(aux1(0:2*trun-1)) !Auxiliary array, for fourier transform
 allocate(Vcross(0:2*trun-1)) ! Stores the cross correlation in reciprocal space
 allocate(ACF(0:2*trun-1))  ! stores  the autocorrelation Fucntion (real space)

 ACF=0.0d0
 
 do ix=1,3
	! Call to the direct FFT
    aux1(0:tread-1) = cmplx(dip_moms(ix,:))
    n=size(aux1)
    call four1(aux1,n,-1)
    Vcross=aux1*conjg(aux1)          
    call four1(Vcross,n,1)
    ACF = ACF + Vcross
 enddo

 !norrmalization so ACF(1) = 1
 !ACF=ACF/real(ACF(1))

 ! Save the correlation function to file
 call io_assign(lun_IR)
 open(unit=lun_IR, file="out_"//trim(fsave)//"_dip_corr_function.dat", form="formatted", status="unknown")
 do i = 1, tread
 	write(lun_IR,*) i*timestep, real(ACF(i))/real(ACF(1))
 enddo	
 call io_close(lun_IR)

 !! Fourier transform the ACF
 aux1=cmplx(ACF)
 call four1(aux1,n,-1)

 !Save the IR spectrum in the file
 call io_assign(lun_IR)
 open(unit=lun_IR, file="out_"//trim(fsave)//"_IR_spectra.dat", form="formatted", status="unknown")

 MinFreqOut =  1d0/(timestep*ps2s*Cspeed*n)
 PointsAvailable = floor(MaxFreqOut/MinFreqOut)
 
 if (PointsAvailable .lt. NumPointsOut) then 
	NumPointsOut = PointsAvailable
	NumAvgOver = 1
 else 
	NumAvgOver = floor(PointsAvailable/NumPointsOut)
 endif 
	
 write(*,*) "-----------infrared spectrum ------------------------"
 write(*,*) "trun = ", trun
 write(*,*) "n, tread = ",  n, tread
 write(*,*) "MinFreqOut = ",  MinFreqOut
 write(*,*) "MaxFreqOut = ",  MaxFreqOut
 write(*,*) "points available = ", PointsAvailable
 write(*,*) "Averaging over", NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 Do t = 0, NumPointsOut-1

  avgMag = 0 
  
  !block averaging
  do i = 1, NumAvgOver
	omega=( t*NumAvgOver+i )/(timestep*n*ps2s) !get freq in 1/s (Hz
	avgMag = avgMag + omega*tanh(hbar*omega*Cspeed/(Kb*2.0d0*Temp))*sqrt(real(aux1(mod(t*NumAvgOver+i,n)) )**2 &
						+ aimag( aux1(mod(t*NumAvgOver+i,n)) )**2)
  enddo
  avgMag = avgMag/real(NumAvgOver)

  ! Use the prefactor with harmonic quantum correction (See Ramirez paper)
  IR = (2d0*3.14159d0*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8) 

  IR = IR/100 !convert from 1/m to 1/cm

  IR = IR/2   !fudge factor (this is probably due to the fact that we have real data in, therefore a 2x redundancy when taking the FT. The FT part was not written by me and whoever wrote it clearly wasn't very careful about normalization. -D. Elton ) 

  omega=( floor((t+.5)*numAvgOver) )/(timestep*n*ps2s) !get central freq of block in 1/s (Hz)
  omega = omega/Cspeed  ! convert Hz to cm-1

  write(lun_IR,*)  omega, IR
 EndDo

 call io_close(lun_IR)

EndSubroutine calc_infrared_spectrum
 
EndModule infrared