      subroutine AvailableIOUnit(hFile)
!***  PORPOSE: to get a IO unit that is not been used
      implicit none
      integer,intent(out)::hFile
      logical::opened

        !--- check the status of the file
         hFile = 10
         do while(.true.)
           !DO hFile = 10, 99
            INQUIRE(UNIT=hFile, OPENED=opened)
            if(.not.opened) exit
            hFile = hFile + 1
         end do
         return
      end subroutine AvailableIOUnit
