!----------------------------------------------------------------------------------
!A simple MPI_timing routine using the MPI wall timer
!It can suppport up to 10 timers, designated by an integer 1,2, ... 10
!The following actions can be called : 
! -- 'start' -- start the timer
! -- 'stop ' -- stop the timer and return elapsed time in seconds
! -- 'write' -- stop timer and return the elapsed time in seconds regardless of timer state
!Copyright 12/2013 D. Elton 
!----------------------------------------------------------------------------------
subroutine MPItimer(timerNum,action,seconds)
use mpi 
Implicit None
integer, intent(in) :: timerNum
character (len=5) :: action
double precision, save :: temp(10)
double precision, save :: totaltimes(10)
real(8), intent(out)  :: seconds
logical, dimension(:), allocatable, save :: STATE

!Trick to make sure that all timers start "off" 
if (.not. allocated(STATE) ) then
	allocate(STATE(10))
	STATE = .false.
	totaltimes = 0 
endif

if (action == 'start') then
	if (STATE(timernum) .eqv. .true.) then
		write(*,*) "ERROR: There is a mistake in the timing scheme. Timer ", timernum, "is already on." 
		stop
	else 
		STATE(timernum) = .true.
		temp(timernum) = MPI_Wtime()
	endif 
endif

if (action == 'stop ') then
	if (STATE(timernum) .eqv. .true.) then
		totaltimes(timernum) = totaltimes(timernum) + MPI_Wtime() - temp(timernum)
		STATE(timernum) = .false.
		seconds = totaltimes(timernum)
	else
		write(*,*) "ERROR: There is a mistake in the timing scheme. Cannot stop timer ", timernum, "because it is not on"
		stop
	endif 
endif 

if (action == 'write') then 
	if (STATE(timernum) .eqv. .true.) then 	
		totaltimes(timernum) = totaltimes(timernum) + MPI_Wtime() - temp(timernum)
		STATE(timernum) = .false.
	endif 
	seconds = totaltimes(timernum)
endif 


end subroutine

