!************************************************
      subroutine GetPath(fullname, path)
!*** to get the path from a full name with path
      implicit none
      character*(*)::fullname, path
      integer::I, flag

            !I = scan(fullname,"\",.true.)-1
           !fullname = fullname(1:I)
            flag = 0
            do I=1, len_trim(fullname)
               !iachar('\') == 92
               !iachar('/') == 47
               if(iachar(fullname(I:I)) .eq. 92 .or.
     -            iachar(fullname(I:I)) .eq. 47) flag = I
            enddo
            if(flag .gt. 0) then
               path = fullname(1:flag)
            else
               path = ""
            end if

      return
      end subroutine GetPath
!************************************************

!************************************************
      subroutine GetParentPath(fullname, path)
!*** to get the path from a full name with path
      implicit none
      character*(*)::fullname, path
      integer::I, flag

            flag = 0
            do I=1, len_trim(fullname)
               !iachar('\') == 92
               !iachar('/') == 47
               if(iachar(fullname(I:I)) .eq. 92 .or.
     -            iachar(fullname(I:I)) .eq. 47) flag = I
            enddo
            if(flag .gt. 0) then
               path = fullname(1:flag-1)
            else
               path = ""
               return
            end if
            !--- do it above process again
            flag = 0
            do I=1, len_trim(path)
               !iachar('\') == 92
               !iachar('/') == 47
               if(iachar(path(I:I)) .eq. 92 .or.
     -            iachar(path(I:I)) .eq. 47) flag = I
            enddo
            if(flag .gt. 0) then
               path = fullname(1:flag)
            else
               path = ""
               return
            end if

      return
      end subroutine GetParentPath
