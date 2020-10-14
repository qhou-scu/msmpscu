      subroutine extract_substr1( String, substr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract the first substring, seperated by space or comma
!  Input:   STRING a string
!  Ouput:   SUBSTR the substring founded
!
! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer  i
      character *(*) String
      character *(*) substr

      i = 1
      substr = ""
      do while(.true.)
         if(String(i:i) .eq. ', ' .or. String(i:i) .eq. ' ' ) then
           exit
         else
           substr(i:i) = String(i:i)
         end if

         i=i+1
         if(i .gt. len_trim(string) .or. i .gt. len(substr)  ) then
            exit
         end if
      enddo

      return
      end subroutine extract_substr1
