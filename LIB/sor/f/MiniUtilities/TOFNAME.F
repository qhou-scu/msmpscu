      subroutine toFname(str)
      implicit none
      character*(*)::str
      integer::I, flag

        !I = scan(str,"\",.true.)
        !str(1:len(str)-I) = str(I+1:len(str) )
        !I = scan(str,".")
        !str(I:len(str)) = " "
            do I=1, len_trim(str)
               !iachar('\') == 92
               !iachar('/') == 47
               if(iachar(str(I:I)) .eq. 92 .or.
     c            iachar(str(I:I)) .eq. 47) flag = I
            enddo
            str = str(flag+1:len_trim(str))

           ! filterout the extention
            do I=1, len_trim(str)
               if(iachar(str(I:I)) .eq. iachar('.')) then
                  flag = I
                  exit
               end if
            enddo
            str = str(1:flag-1)
      return
      end subroutine toFname
