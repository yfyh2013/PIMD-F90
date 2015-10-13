!-----------------------------------------------
! Split line reader 
!
!------------------------------------------------------

module split_line_reader
implicit none 

character(len=*) :: line

logical :: is_digit
integer :: i 
character(1) :: test_char

!!One line function declarations 
is_digit(i) = (i .ge. 48) .and. (i .le. 57)
is_upper(i) = (i .ge. 65) .and. (i .le. 90)
is_lower(i) = (i .ge. 97) .and. (i .le. 122)
is_alpha(i) = is_upper(i) .or. is_lower(i)
is_alnum(i) = is_digit(i) .or. is_alpha(i)

is_comment(i) = (i .eq. 33) .or. (i .eq. 35) .or. (i .eq. 59) !Comment characters:  !  #  ; 
is_delstr(i)  = (i .eq. 34) .or. (i .eq. 39) .or. (i .eq. 96) !String delimiters: "  '  `




Function make_lower(astr)
Implicit None
character(len=*), intent(inout) :: astr
integer i,j

 is_upper(i) = (i .ge. 65) .and. (i .le. 90)
 is_lower(i) = (i .ge. 97) .and. (i .le. 122)

 do i = 1, len(astr)
	j = ichar(astr(i)) 
	if is_upper(j) then
		astr(i) = achar(j + 32)
	endif
 enddo

End Function make_lower




      
end program test_fun