!****************************************************************
      SUBROUTINE GETINPUTSTRLINE(hFile,STR, LINE, FILTER, *)
!*** to read a non-comment line fropm a I/O unit
!*** INPUT: hFile    the I/O unit
!    OUTPUT: STR, the string on of the line
!            LINE, the accumulated number of line
!    HISTORY: to be compatiable to DOS and UNIX text file
!             changed hase been made to extract
!             only visible characters
!    AUTHOR: Hou Qing, updated 2014-06-03
      IMPLICIT NONE
      integer, intent(in)::hFile
      integer, intent(inout)::LINE
      character*(*)::STR
      character*(*), intent(in)::FILTER
      character(len=1024)::TSTR
      character(len=1)::FC(5)=''
      integer::I, IP, IC, LS

          FC= ''
          do I=1, min(size(FC), len_trim(FILTER))
             FC(I) = FILTER(I:I)
          end do

          LS = len(STR)
          do while(.true.) !.NOT.EOF(hFile))
             line = line + 1
             read(hFile,fmt="(A1024)", end=100 )TSTR
             IP=0
             STR = ""
             do I=1, len_trim(TSTR)
                IC = iachar(TSTR(I:I))
                if(IC .GE. 20 .and. IC .LE. 126) then
                   IP = IP + 1
                   STR(IP:IP) = TSTR(I:I)
                   if(IP.GE.LS) return
                end if
             end do
             STR = adjustL(STR)
             if(len_trim(STR).le.0) cycle
             if(any(FC .eq. STR(1:1))) cycle
             !if(STR(1:1) .ne. '!' .AND. LEN_TRIM(STR) .GT. 0) then
              return
             !end if
          end do
  100     return 1
       END SUBROUTINE GETINPUTSTRLINE
!****************************************************************
