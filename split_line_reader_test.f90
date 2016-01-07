!----------- Split Line Reader --------------------------------
!
! A Fortran-90 scheme for reading a variable format input file 
! 
! 2015 Daniel Elton 
!
!--------------------------------------------------------------
program split_line_reader_test
	use split_line_mod
Implicit None
integer :: lun=20
integer :: variable1
real 	:: variable2
character(len=100) :: variable3

call split_line_process_file("example.inp") 





end program split_line_reader_test



