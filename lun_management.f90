!---------------------------------------------------------------
! Fortan-90 module for logical unit management 
! By using logical unit management, we ensure logical units for reading/writing 
! do not conflict with those used by other code that is interfaced with or called by this code
!
! Based on code by Alberto Garcia and Richard Maine 
!
! Copyright 2015 D. Elton 
!---------------------------------------------------------------
Module lun_management
 implicit none
 integer, save 	    :: stdout=6, stderr=0
 integer   	    :: iostat 
 integer, parameter :: min_lun=10, max_lun=99
 integer, parameter :: nunits = max_lun - min_lun+1 
 logical, save      :: lun_is_free(min_lun:max_lun)
 logical 	    :: used
 
 !initially we assume all units between 10 and 99 are free
 data lun_is_free /nunits*.true./

 contains 
 
subroutine io_get_stdout(lun)
  integer, intent(out) :: lun
  lun = stdout
end subroutine

!---------------------------------------------------
!Looks for a free unit and assigns it to lun
!---------------------------------------------------
subroutine io_assign(lun)
   integer, intent(inout) :: lun

   do lun= min_lun, max_lun
      if (lun_is_free(lun)) then
	  inquire(unit=lun, opened=used, iostat=iostat)
          if (iostat .ne. 0) used = .true.
          if (used .eqv. .true.) then 
	      lun_is_free(lun) = .false.
	  else 
	       return  
	  endif
      endif
    enddo
end subroutine io_assign    
 
!---------------------------------------------------
!Close a unit 
!---------------------------------------------------
subroutine io_close(lun)
      integer, intent(in) :: lun
      close(lun)
      if ((lun .ge. min_lun) .and. (lun .le. max_lun)) lun_is_free(lun) = .true.
end subroutine io_close

!---------------------------------------------------
!Reserve a unit (useful for compatibility with legacy code)
!eg. call io_reserve(15) before open(15) appears
!---------------------------------------------------
subroutine io_reserve(lun)
      integer, intent(in) :: lun
      inquire(unit=lun, opened=used, iostat=iostat)
      if (iostat .ne. 0) used = .true.
      if (used) then
	 write(*,*) "Error in logical unit scheme. Cannot reserve unit ", lun, " since it is in use."
      else
      	 if (lun .ge. min_lun .and. lun .le. max_lun) lun_is_free(lun) = .false.
      endif
end subroutine io_reserve 

end module 






