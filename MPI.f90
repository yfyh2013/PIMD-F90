!----------------------------------------------------------------------------------------
!This is a module of stub functions and variables which 
!allow MPI calls to be executed even if the code is compiled with a serial (traditional) compiler
!
! MPI_Send and MPI_Recv may be called with various objects (arrays, scalars, etc) in the 'input' argument
! The code is written assuming that the 'input' argument is an array, but this may not be convient 
! The use of the 'unlimited' polymorphism feature of F2003 does not appear to solve this problem 
!
! Obvioulsy this way of allowing serial code compilation may be viewed as inelegant, but it does allow 
! For consistent use of MPI calls (ie, for timing, getting the number of nodes (=1), etc) during serial code execution
!
! Alternatively, this module may be removed and preprocessor declarations (ifdef mpi .. endif) may be placed around 
! around all MPI calls instead. 
!
! Copyright 2014-2015 Daniel C. Elton 
!----------------------------------------------------------------------------------------
module mpi
Integer, parameter :: MPI_STATUS_SIZE=1
Integer, parameter :: MPI_COMM_WORLD=1
Integer, parameter :: MPI_DOUBLE_PRECISION=8
 
 contains


subroutine MPI_Init(ierr)

end subroutine MPI_Init



subroutine MPI_Finalize(ierr)

end subroutine MPI_Finalize




subroutine MPI_Send(input, counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr) 
 integer ::  counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr
 double precision, intent(inout), dimension(counti) :: input
 !class(*) :: input
end subroutine MPI_Send


subroutine MPI_Recv(input, counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, status2, ierr)
 integer  ::  status2(MPI_STATUS_SIZE)
 integer   ::  counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr
 double precision, intent(inout), dimension(counti) :: input
 !class(*)  :: input
end subroutine MPI_Recv



subroutine MPI_Comm_size(MPI_COMM_WORLD, Nnodes, ierr)
	Nnodes = 1
end subroutine MPI_Comm_size




subroutine MPI_Comm_rank(MPI_COMM_WORLD, pid, ierr)
	integer, intent(out) :: pid
	pid = 0 
end subroutine MPI_Comm_rank




subroutine MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine



function MPI_Wtime()
	double precision ::  MPI_Wtime
        call cpu_time(MPI_Wtime)
end function


end module
