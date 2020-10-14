      subroutine getFname(fullname, name)
!*** to get the file name from a full name with path
!    AUTHOR: Hou Qing, updated 2014-06-03
      implicit none
      character*(*)::fullname
      character*(*)::name
      integer::I, flag

            flag = 0
            do I=1, len_trim(fullname)
               !iachar('\') == 92
               !iachar('/') == 47
               if(iachar(fullname(I:I)) .eq. 92 .or.
     -            iachar(fullname(I:I)) .eq. 47) flag = I
            enddo
            name = adjustl(fullname(flag+1:len_trim(fullname)))
      return
      end subroutine getFname
