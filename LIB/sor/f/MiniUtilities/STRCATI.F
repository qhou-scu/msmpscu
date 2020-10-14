       subroutine STRCATI(str, str0, sep, N, L)
!*** To catenate a character number to str
!
!       Author: HOU Qing, 2010-11
!
!    INPUT:  str0, the string be catenated a number
!            sep,  the seperation between str0 and the number
!            N,    the number
!            L.    the max lenght  of the number charscter
!    OUTPUT: str, the new string

!--- DUMMY VARIABLES:
      implicit none
      character *(*), intent(in)::str0, sep
      character *(*), intent(out)::str
      integer, intent(in)::N
      integer, intent(in)::L
      character*32::NSTR, ZSTR
      integer::I, M

        write(NSTR,*)N
        NSTR = ADJUSTL(NSTR)
        M = L-len_trim(NSTR)
        ZSTR = ""
        DO I=1, M
           ZSTR = ZSTR(1:len_trim(ZSTR))//"0"
        END DO
        str = str0(1:len_trim(STR0))//SEP(1:len_trim(SEP))
        str = str(1:len_trim(STR))//ZSTR(1:len_trim(ZSTR))
        str = str(1:len_trim(STR))//NSTR(1:len_trim(NSTR))
      return
      end subroutine STRCATI
