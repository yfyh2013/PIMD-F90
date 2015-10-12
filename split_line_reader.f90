program test_fun
implicit none 

logical :: is_digit
integer :: i 
character(1) :: test_char

       is_digit(i) = (i .ge. 48) .and. (i .le. 57)

test_char = " "

write(*,*) test_char
!write(*,*) real(test_char)
write(*,*) ichar(test_char)
write(*,*) is_digit(ichar(test_char))




is_upper(i) = (i .ge. 65) .and. (i .le. 90)
is_lower(i) = (i .ge. 97) .and. (i .le. 122)
is_alpha(i) = is_upper(i) .or. is_lower(i)
is_alnum(i) = is_digit(i) .or. is_alpha(i)

!Comments are signaled by:  !  #  ; 
is_comment(i) = (i .eq. 33) .or. (i .eq. 35) .or. (i .eq. 59)

!String delimiters: "  '  `
      is_delstr(i)  = (i .eq. 34) .or. (i .eq. 39) .or. (i .eq. 96)

!     List delimiters: { }
      is_dellist(i)  = (i .eq. 123) .or. (i .eq. 125)

!     Special characters which are tokens by themselves: <
      is_special(i) = (i .eq. 60)
      
      
end program test_fun