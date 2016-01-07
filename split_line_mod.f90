!-----------------------------------------------
! Split line input file reader 
! Under preliminary construction
! Not yet implemented!! 
!------------------------------------------------------
module split_line_reader
 Implicit None
 MAX_LINE_LENGTH = 1000
 
 type KeyValue
	character(len=100) :: key
	double precision   :: value
	character(len=100) :: strvalue
 end type KeyValue



character(len=*) :: line

logical :: is_digit
integer :: i 
character(1) :: test_char

contains 


!-----------------------------------------------------------------
!--  Process an imput file and place into database
!-----------------------------------------------------------------
Subroutine split_line_process_file(filename) 
 Implicit None
 character(len=200), intent(in) :: filename 
 character(MAX_LINE_LENGTH) :: aline
 integer :: lun=20
 
 !!One line function declarations 
 is_digit(i) = (i .ge. 48) .and. (i .le. 57)
 is_upper(i) = (i .ge. 65) .and. (i .le. 90)
 is_lower(i) = (i .ge. 97) .and. (i .le. 122)
 is_alpha(i) = is_upper(i) .or. is_lower(i)
 is_alnum(i) = is_digit(i) .or. is_alpha(i)

 is_comment(i) = (i .eq. 33) .or. (i .eq. 35) .or. (i .eq. 59) !Comment characters:  !  #  ; 
 is_delstr(i)  = (i .eq. 34) .or. (i .eq. 39) .or. (i .eq. 96) !String delimiters: "  '  `
 
 open(lun, file=filename, status='old')

 !main loop
 do 
	read(lun,*)  aline  !read line

	aline = trim(aline) 
		
	if is_blank(aline) goto 1000
	if is_comment(aline(1)) goto 1000
	
	for i = 1, len(aline) 
		if 
	
	
 
1000  continue

 enddo
 

End Subroutine split_line_process_file




!-----------------------------------------------------------------
!--  Make a string lower case
!-----------------------------------------------------------------
Function make_lower(astr) result(astr)
 Implicit None
 character(len=*), intent(inout) :: astr
 integer i,j

 is_upper(i) = (i .ge. 65) .and. (i .le. 90)

 do i = 1, len(astr)
	j = ichar(astr(i)) 
	if is_upper(j) then
		astr(i) = achar(j + 32)
	endif
 enddo

End Function make_lower

      
end program test_fun
