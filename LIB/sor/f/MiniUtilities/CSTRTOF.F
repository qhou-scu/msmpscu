      subroutine CStrtoF(str)
      implicit none
      character*(*)::str
      integer::i

      do i=1, len(str)
         if(iachar(str(i:i)) .eq. 0) exit
      end do

      str(i:len(str)) = " "
      return
      end subroutine CStrtoF
