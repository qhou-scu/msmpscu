module My_MiniUtilities_QMJ
!***  DESCRIPTION: This module contains several mini-utilities used in my programs.
!
!
!*** AUTHOR:
!    Mingjie Qiu, Oct., 2016
!    Thank you for the guidance of my teacher Qing Hou!
!
use MiniUtilities
implicit none
contains
!**************************************************************************
subroutine Get_FileExtent (InputFile, ext)
! Purpose :  to get the InputFile extension
!
!
  implicit none
  !--- dummy variable
  character*(*) :: InputFile
  character*(*) ,  intent(out):: ext
  !--- local variable
  character*256 :: str
  integer :: I , IC , L
  InputFile = adjustl(InputFile)
  str = InputFile(1:len_trim(InputFile))
  L   = len_trim(str)
  do I = 1, L
    if (iachar(str(I:I)) .EQ. iachar(".")) then
      IC = I
    end if
  end do
  ext = str(IC+1:L)
end subroutine Get_FileExtent
!*********************************************************************************
!*********************************************************************************
subroutine EXSTR (STRING, mxN, N, NSTR)
  ! Purpose:   to extract substrings between two ' ' .
  !   Input:   STRING a string
  !   Ouput:   N  the number of substrings found in the string
  !            NSTR , the array which is used to put substrings in
  !  Auther:   Mingjie Qiu, student of Qing Hou, Inst. of Nucl. Sci.and Tech., Sichuan Union University
  !            Thank you for the guidance of my teacher Qing Hou!
implicit none
!--- dummy variables
  character*(*), intent(in) :: string
  integer, intent(in) :: mxN
  integer, intent(out) :: N
  character*(*), dimension(mxN),intent(out) :: NSTR
!--- local variables
  integer :: I                 !*** loop index
  character*256::tstr
!********
  tstr = adjustl(string(1:len_trim(string)))
  N = 0
  I = 0
  do while(len_trim(tstr) .gt. 0.)
     I = I + 1
     if(tstr(I:I) .eq. ' ') then
        N = N + 1
        NSTR(N) = tstr(1:I-1)
        tstr = adjustl( tstr(I:len_trim(tstr)))
        I = 0
      else
        if(I .ge. len_trim(tstr)) then
           N = N + 1
           NSTR(N) = tstr(1:I)
           exit
        end if
      end if

   end do
end subroutine EXSTR
!*********************************************************************************
!*********************************************************************************
subroutine Count_NL_NC ( InputFile, NL, NC )   !*** This subroutine is not good. 2016.11.19 by QMJ
  ! Purpose:   to count the number of line and column in InputFile.
  !   Input:   InputFile
  !   Ouput:   NL , the number of data line in tne file
  !            NC , the number of data column in the file
  !  Auther:   Mingjie Qiu, student of Qing Hou, Inst. of Nucl. Sci.and Tech., Sichuan Union University
implicit none
  !--- dummy variables
  character*(*) , intent(in)  :: InputFile     !*** the input file which must have the same data lines in each column
  integer       , intent(out) :: NL            !*** the number of data lines
  integer       , intent(out) :: NC            !*** the number of data columns
  !--- local variables
  integer :: hFile, LINE , flag, N
  character*256 :: STR, STRING
  character*32  :: NSTR(30)
  !***
  call AvailableIOUnit(hFile)
    open (UNIT=hFile,FILE= InputFile, STATUS='old')
      call GetInputStrLine(hFile,STR, LINE, "!", *100)
      NL = 0
      NC = 0
      flag = 0
      do while(.true.)
        read (hFile,fmt= '(A)',end = 100) STRING       !*** get the total number of files in each column
        if (flag == 0) then
          call EXSTR (STRING, 30, N, NSTR)
          NC = N
          flag = 1
        end if
       if (len_trim(STRING).gt.0)  NL = NL + 1
      end do
  100 rewind(hFile)
      close (hFile)
  end subroutine Count_NL_NC
!*********************************************************************************
end module My_MiniUtilities_QMJ
