  module MD_TYPEDEF_RecordStamp
  !***  DESCRIPTION:
  !     This module provides a typedef that collects the information about
  !     time, and generator type of simulation boxes
  !
  !     SEE ALSO:
  !
  !
  !    Written by HOU Qing, Dec, 2016
  !
   use MD_CONSTANTS
   implicit none

      type MDRecordStamp
           character*8::    AppType     = ""
           integer::        IProcess    = 0
           integer::        ITest       = -1
           integer::        NBox        = 0
           integer::        IBox(2)     = -1
           integer::        ITime       = -1
           integer::        ISect       = -1
           real(KINDDF)::   Time        = -1
           real(KINDDF)::   ScalTime    = -1
           real(KINDDF)::   InstantTemp = -1
           integer::        ICfg(2)     = -1
           integer::        IRec(2)     = -1
           character*12::   ClockDate   =""
           character*12::   ClockTime   =""

           !logical::        NeedCheck   = .true.   !flag indicating if consistent needed to be chcek on load in data
      end type MDRecordStamp

  contains

  !****************************************************************************************
  subroutine Default_RecordStamp(Stamp)
  !***  DESCRIPTION: to set default values for a RecordStamp
  !
  implicit none
     !--- DUMMY variables
      type(MDRecordStamp)::Stamp
     !--- Local
           Stamp%ITest      = -1
           Stamp%IProcess   =  0
           Stamp%NBox       =  0
           Stamp%IBox       = -1
           Stamp%ITime      = -1
           Stamp%ISect      = -1
           Stamp%Time       = -1
           Stamp%ScalTime   = -1
           Stamp%InstantTemp= -1
           Stamp%ICfg       = -1
           Stamp%IRec       = -1

           Stamp%ClockDate  = ""
           Stamp%ClockTime  = ""

           return
   end subroutine Default_RecordStamp
  !****************************************************************************************


  !****************************************************************************************
  subroutine Copy_RecordStamp(From, To)
  !***  DESCRIPTION: to set default values for a RecordStamp
  !
  implicit none
     !--- DUMMY variables
      type(MDRecordStamp)::From, To
     !--- Local

           To%AppType     = From%AppType
           To%ITest       = From%ITest
           To%NBox        = From%NBox
           To%IBox        = From%IBox
           To%ITime       = From%ITime
           To%ISect       = From%ISect
           To%Time        = From%Time
           To%ScalTime    = From%ScalTime
           To%InstantTemp = From%InstantTemp
           To%ICfg        = From%ICfg
           To%IRec        = From%IRec

           return
   end subroutine Copy_RecordStamp
  !****************************************************************************************

  !****************************************************************************************
  subroutine Putout_RecordStamp(hFile, Stamp)
  !***  PORPOSE: to putout the current record stamp
  !     INPUT: hFile, the IO unit
  !            Stamp, the record stamp
  implicit none
      !--- DUMMY variables
      integer, intent(in):: hFile
      type(MDRecordStamp):: Stamp
      !--- Local

         call date_and_time(Date=Stamp%ClockDate, Time=Stamp%ClockTime)
         write(hFile,FMT="(A)")              '&APPTYPE        "'//Stamp%AppTYpe(1:len_trim(Stamp%AppTYpe))//'" '
         write(hFile,FMT="(A, A)")           "&DATE            ",  Stamp%ClockDate(1:4)//"-"//Stamp%ClockDate(5:6)//"-"//Stamp%ClockDate(7:8)//","//&
                                                              Stamp%ClockTime(1:2)//"h"//Stamp%ClockTime(3:4)//"m"//Stamp%ClockTime(5:6)//"s"

         write(hFile,FMT="(A, I12)")         "&TESTID          ", Stamp%ITest
         write(hFile,FMT="(A, I12,A,I6)")    "&BOXID           ", Stamp%IBox(1),", ", Stamp%IBox(2)
         write(hFile,FMT="(A, I12,A,I6)")    "&CFGID           ", Stamp%ICfg(1),", ", Stamp%ICfg(2)
         write(hFile,FMT="(A, I12,A,I6)")    "&RECID           ", Stamp%IRec(1),", ", Stamp%IRec(2)
         write(hFile,FMT="(A, I12)")         "&TIMESTEPS       ", Stamp%ITime
         write(hFile,FMT="(A, I12)")         "&TIMESECT #      ", Stamp%ISect
         write(hFile,FMT="(A, 1PE12.5, A)")  "&TIME (ps)       ", Stamp%Time
         write(hFile,FMT="(A, 1PE12.5, A)")  "&SCALTIME (ps)   ", Stamp%ScalTime
         write(hFile,FMT="(A, 1PE12.5, A)")  "&TEMPSET (K)     ", Stamp%InstantTemp
         write(hFile,FMT="(A, 1PE12.5, A)")  " "

  end subroutine Putout_RecordStamp
  !****************************************************************************************

  !****************************************************************************
  subroutine Putin_RecordStamp(hFile, Stamp)
  !***  PORPOSE: to putin the current record stamp
  !     INPUT: hFile, the IO unit
  !            Stamp, the record stamp
   use MiniUtilities
   implicit none
      !--- DUMMY variables
      integer,            intent(in) ::hFile
      type(MDRecordStamp),intent(out)::Stamp

      !--- local variables
      integer::LINE, N
      character*256::str
      character*32::KEYWORD, STRNUMB(3)
      !--- local variables for colume

         LINE = 0
            !--- check the consistent of data
            !--- and get stamp information
            do while(.true.)
                call GetInputStrLine(hFile,STR, LINE, "!",*100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                if(len_trim(KEYWORD) .le. 0) exit
                !--- deal witht the keyword
                call UpCase(KEYWORD)
                select case(KEYWORD(1:LEN_TRIM(KEYWORD)) )
                       case("&APPTYPE")
                            call Extract_Substr(STR,1,n,STRNUMB)
                            Stamp%AppType = STRNUMB(1)(1:len_trim(STRNUMB(1)))

                       case("&TESTID")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%ITest = ISTR(STRNUMB(1))

                       case("&BOXID")
                            call Extract_Numb(STR,2,n,STRNUMB)
                            Stamp%IBOX(1) = ISTR(STRNUMB(1))
                            Stamp%IBOX(2) = ISTR(STRNUMB(2))

                       case("&CFGID")
                            call Extract_Numb(STR,2,n,STRNUMB)
                            Stamp%ICfg(1) = ISTR(STRNUMB(1))
                            Stamp%ICfg(2) = ISTR(STRNUMB(2))

                       case("&TIMESTEPS")
                            call Extract_Numb(STR,3,n,STRNUMB)
                            Stamp%ITime = ISTR(STRNUMB(1))

                       case("&TIME")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%Time = DRSTR(STRNUMB(1))

                       case("&TIMESECT")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%ISect = ISTR(STRNUMB(1))

                       case("&SCALTIME")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%ScalTime = DRSTR(STRNUMB(1))

                       case("&INSTANTTEMP", "&TEMPSET")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%InstantTemp = DRSTR(STRNUMB(1))
                end select
            end do
            return

 100   close(hFile)
       write(*,*)"MDPSCU Error: error in putin time stamp"
       write(*,*)"The process to be stopped."
       stop            
  end subroutine Putin_RecordStamp
  !****************************************************************************************

  !****************************************************************************************
  subroutine Archive_RecordStamp(hFile, Stamp)
  !***  PORPOSE:
  !     INPUT: hFile, the IO unit
  !            Stamp, the record stamp
  implicit none
      !--- DUMMY variables
      integer, intent(in)            :: hFile
      type(MDRecordStamp), intent(in):: Stamp
      !--- Local

         write(hFile) Stamp%AppTYpe(1:len(Stamp%AppTYpe))
         write(hFile) Stamp%ITest, Stamp%NBox, Stamp%IBox, Stamp%ICfg, Stamp%IRec, Stamp%ITime, Stamp%ISect,Stamp%Time,Stamp%ScalTime, Stamp%InstantTemp


  end subroutine Archive_RecordStamp
  !****************************************************************************************

  !****************************************************************************************
  subroutine Restore_RecordStamp(hFile, Stamp)
  !***  PORPOSE:
  !     INPUT: hFile, the IO unit
  !            Stamp, the record stamp
  implicit none
      !--- DUMMY variables
      integer, intent(in):: hFile
      type(MDRecordStamp):: Stamp
      !--- Local

         read(hFile) Stamp%AppTYpe(1:len(Stamp%AppTYpe))
         read(hFile) Stamp%ITest, Stamp%NBox, Stamp%IBox, Stamp%ICfg, Stamp%IRec, Stamp%ITime, Stamp%ISect,Stamp%Time,Stamp%ScalTime, Stamp%InstantTemp

  end subroutine Restore_RecordStamp
  !****************************************************************************************



  !****************************************************************************************
  end module MD_TYPEDEF_RecordStamp
  !****************************************************************************************
