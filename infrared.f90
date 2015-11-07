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
subroutine calc_infrared_spectrum(dip_moms,box,timestep,fsave)
 use lun_management 
 implicit none
 integer, parameter :: dp=kind(1.0d0)
 real(dp), dimension(:,:), intent(in) :: dip_moms 
 real(dp), dimension(3), intent(in) :: box ! box size in Ang
 real(dp), intent(in) :: timestep ! times step in fs
 character(len=*) :: fsave  
 integer :: trun, tcor, ix, i, t, n, count, tread, lun_IR
 real(dp) :: omega, IR, magn, avgMag, MaxFreqOut, vol, PointsAvailable
 complex, dimension(:), allocatable :: aux1,aux2,vcross,ACF
 real(dp), parameter :: Cspeed=3.00d10 ! cm/s
 real(dp), parameter :: Kb=1.38d-23 !J/K
 real(dp), parameter :: Temp=300 !K
 real(dp), parameter :: hbar=6.626d-34 ! J*s
 real(dp), parameter :: fs2s=1d-15 ! 1fs in s
 real(dp), parameter :: Debye2SI = 3.33564d-30
 !real(dp) :: pi, A0,A1,A2,A3, window
 integer :: NumAvgOver, NumPointsOut

 NumPointsOut = 150
 MaxFreqOut = 4000 !maximum frequency to plot (cm^-1)

 tread = size(dip_moms,2)

 !find closest power of 2 less than tread
 trun = 2**(    floor( dlog(  dble(tread) )/dlog(2d0)  )  )

 allocate(aux1(0:2*trun-1)) !Auxiliary array, for fourier transform
 allocate(aux2(0:2*trun-1)) ! second Aux arrays
 allocate(Vcross(0:2*trun-1)) ! Stores the cross correlation in reciprocal space
 allocate(ACF(0:2*trun-1))  ! stores  the autocorrelation Fucntion (real space)

 ACF=0.0d0
 do ix=1,3
	! Call to the direct FFT
    aux1(0:tread-1) = cmplx(dip_moms(ix,:))
    aux2(0:tread-1) =cmplx(dip_moms(ix,:))
    n=size(aux1)
    call four1(aux1,n,-1)
    call four1(aux2,n,-1)
    Vcross=aux1*conjg(aux2)          
    call four1(Vcross,n,1)
    Vcross=Vcross/real(size(Vcross))
    ACF = ACF + Vcross
 enddo

 !norrmalization so ACF(1) = 1
 ACF=ACF/real(ACF(1))

 !! Save the correlation function to file
 !open(unit=10, file="corr_function.dat", form="formatted", status="unknown")
 !do i = 1, tread
 !	write(10,*) i*timestep/1000, real(ACF(i))
 !enddo	
 ! close(10)

 !! Fourier transform the ACF
 aux1=cmplx(ACF)
 call four1(aux1,n,-1)

 !Save the IR spectrum in the file
 call io_assign(lun_IR)
 open(unit=lun_IR, file="out_"//fsave//"_IR_spectra.dat", form="formatted", status="unknown")

 PointsAvailable = MaxFreqOut/( (1/(timestep*fs2s))/Cspeed )*n
 if (PointsAvailable .lt. NumPointsOut) NumPointsOut = PointsAvailable
 NumAvgOver =PointsAvailable/NumPointsOut 

!write(*,*) "-----------infrared spectrum ------------------------"
!write(*,*) "n = ",  n
!write(*,*) "MaxFreqOut = ",  MaxFreqOut
!write(*,*) "points available = ", PointsAvailable
!write(*,*) "Averaging over", NumAvgOver

 vol = (box(1)*box(2)*box(3))*1d-30

 Do t = 0, NumPointsOut-1

  avgMag = 0 
  do i = 1, NumAvgOver
	omega=( t*NumAvgOver+i )/(timestep*n*fs2s) !get freq in 1/s (Hz
	avgMag = avgMag + omega*tanh(hbar*omega*Cspeed/(Kb*2.0d0*Temp))*sqrt(real(aux1(mod(t*NumAvgOver+i,n)) )**2 &
						+ aimag( aux1(mod(t*NumAvgOver+i,n)) )**2)
  enddo
  avgMag = avgMag/real(NumAvgOver)

  ! Use the prefactor with harmonic quantum correction (See Ramirez paper)
  IR = (2d0*3.14159d0*(Debye2SI**2)*avgMag)/(3d0*vol*2.99d8) 

  IR = IR/100 !convert from 1/m to 1/cm

  IR = IR/2   !fudge factor (this is probably due to the fact that we have real data in, therefore a 2x redundancy when taking the FT. The FT part was not written by me and whoever wrote it clearly wasn't very careful about normalization. -D. Elton ) 

  omega=( floor((t+.5)*numAvgOver) )/(timestep*n*fs2s) !get central freq in 1/s (Hz
  omega = omega/Cspeed  ! convert frequency to cm-1

  write(lun_IR,*)  omega, IR
 EndDo

 call io_close(lun_IR)

EndSubroutine calc_infrared_spectrum
 
EndModule infrared