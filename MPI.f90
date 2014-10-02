!This is a module of stub functions/variables for MPI for running on serial machines
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

end subroutine MPI_Send




 subroutine MPI_Recv(input, counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, status2, ierr)
 integer  ::  status2(MPI_STATUS_SIZE)
 integer   ::  counti, MPI_DOUBLE_PRECISION, i, MPI_TAG, MPI_COMM_WORLD, ierr
 double precision, intent(inout), dimension(counti) :: input

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
