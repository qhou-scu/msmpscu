      real function rstr(Text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This function is set to convert a text into a real number !
!  with double precision                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision tem
      character *(*) Text

       read(text,*)tem
       rstr = tem

       return
      end function
