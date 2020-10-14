      subroutine FStrtoc(str)
      implicit none
      character*(*)::str
      integer::i

      i = len(str)
       str(i:i) = achar(0)

      return
      end subroutine FStrtoc
