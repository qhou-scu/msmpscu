  module MSM_TYPEDEF_InputPaser
  !***  DESCRIPTION:
  !     This module provides a typedef that gather the strings load from a file.
  !
  !
  !
  !    Written by HOU Qing, Nov., 2015
  !
   use MSM_CONSTANTS
   use MiniUtilities
   implicit none

  !--- interface to the external routine -------------------
      type::Statementlist
           integer                     ::line = 0
           character(len=256)          ::this = ""
           type(StatementList), pointer::next=>null()
      end type StatementList

      type InputStatements
           character(len=256)          ::filename = ""    !--- the filename, which provides where the inputs are loaded from
           character(len=64)           ::stag=""          !    the Tag, which used to distingush the modules that will use the input
           type(StatementList), pointer::stat=>null()
      end type InputStatements

      private  :: GetbyKWD_StatementList, &
                  FindByKWD_StatementList
      interface   Get_StatementList
         module procedure GetbyKWD_StatementList
         module procedure GetbyIND0_StatementList
         module procedure GetbyIND1_StatementList
         module procedure FindByKWD_StatementList
         module procedure FindByInd_StatementList
      end interface

  contains

 !****************************************************************************
  recursive subroutine Add_StatementList(Statements, Str, Line)
  !***  PURPOSE:   to add a statement to a list
  !
  !     INPUT:    Statements, the statement list
  !               Str,        the staement to be added
  !               Line,       the integer flag of the statment
  !
  !     OUTPUT    Statements
  !
  implicit none
     !--- dummy varioables
     type(StatementList),       pointer  ::Statements
     character*(*), intent(in)           ::Str
     integer,       intent(in), optional ::LINE

     !----
          if(.not.associated(Statements)) then
             allocate(Statements)
             Statements%this = Str(1:min(len(Str), len(Statements%this)))
             if(present(Line)) then
                Statements%Line = Line
             else
                Statements%Line = 0
             end if
             Statements%Next => null()
          else
             if(present(Line)) then
                call Add_StatementList(Statements%Next, Str, Line)
             else
                call Add_StatementList(Statements%Next, Str)
             end if
          end if
          return
  end subroutine Add_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine Copy_StatementList(From, To,  EndKWD)
  !***  PURPOSE:   to copy a statementList
  !
  !     INPUT:     From,  the statement list copy from
  !                EndKWD, optional, the end keyword, copy the list until the keyword found
  !
  !     OUTPUT     To,    the statmentlist copy into
  implicit none
     !--- dummy varioables
     type(StatementList), intent(in), pointer ::From
     type(StatementList),             pointer ::To
     character*(*), intent(in),       optional::EndKWD
     !----
     character*256::tSTR
     character*64 ::tWord

     !----
              call Release_StatementList(To)
              if(.not.associated(From)) then
                 return
              end if
              allocate(To)
              To%line =  From%line
              To%this =  From%this
              To%next => null()
              if(present(EndKWD)) then
                 tSTR = adjustl(To%this)
                 call GetKeyWord("&", tSTR, tWord)
                 call UpCase(tWord)
                 if(tWord(1:len_trim(tWord)) .eq. EndKWD(1:len_trim(EndKWD)) ) then
                    return
                 end if
                 call Copy_StatementList(From%next, To%next, EndKWD)
              else
                 call Copy_StatementList(From%next, To%next)
              end if

              return
  end subroutine Copy_StatementList
!****************************************************************************

!****************************************************************************
  recursive subroutine Load_StatementList(hFile, Statements, STR, LINE)
  !***  PURPOSE:   to load statements from a I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     Statements
  !                STR, LINE    working space
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(StatementList),pointer::Statements
     integer::LINE
     character*(*)::STR
     !----

              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              allocate(Statements)
              Statements%line = LINE
              Statements%this = STR
              Statements%next=>null()
              call Load_StatementList(hFile, Statements%next, STR, LINE)
              return

  100         return
 !----------------------------------------------------------------------------
  end subroutine Load_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine Print_StatementList(hFile, Statements)
  !***  PURPOSE:   to write out statements to a I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     Statements
  !
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(StatementList),pointer::Statements

     !----
              if(.not.associated(Statements)) then
                 return
              end if

              write(hFile, fmt="(I8,1X,A)") Statements%line, Statements%this(1:len_trim(Statements%this))
              call Print_StatementList(hFile, Statements%next)
              return
  end subroutine Print_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine Write_StatementList(hFile, Statements)
  !***  PURPOSE:   to write out statements to a I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     Statements
  !
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(StatementList),pointer::Statements

     !----
              if(.not.associated(Statements)) then
                 return
              end if

              write(hFile, fmt="(A)") Statements%this(1:len_trim(Statements%this))
              call Write_StatementList(hFile, Statements%next)
              return
  end subroutine Write_StatementList
 !****************************************************************************

 !****************************************************************************
  subroutine Archive_StatementList(hFile, Statements)
  !***  PURPOSE:   to write out statements to a I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     Statements
  !
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(StatementList),pointer::Statements

     !---- Local
     type(StatementList),pointer::Next
     integer::N

              N = Number_StatementList(Statements)
              write(hFile) N
              Next => Statements
              do while(associated(Next))
               write(hFile) Next%LINE
                 write(hFile) Next%THIS
                 Next => Next%NEXT
              end do
              return
  end subroutine Archive_StatementList
 !****************************************************************************

 !****************************************************************************
  subroutine Restore_StatementList(hFile, Statements)
  !***  PURPOSE:   to write out statements to a I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     Statements
  !
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(StatementList),pointer::Statements

     !---- Local
     integer::I, N, LINE 
     character(len=256)          ::STR
              read(hFile) N 
              do I=1, N
                 read(hFile) LINE 
                 read(hFile) STR
                 call Add_StatementList(Statements, STR, LINE)
              end do   
              return
  end subroutine Restore_StatementList
 !****************************************************************************  

 !****************************************************************************
  recursive subroutine Release_StatementList(Statements)
  !***  PURPOSE:   to load statements from a I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     Statements
  !                STR, LINE    working space
  implicit none
     !--- dummy varioables
     type(StatementList),pointer::Statements
     !----

              if(.not.associated(Statements)) then
                 return
              end if

              call Release_StatementList(Statements%next)
              deallocate(Statements)
              Statements => null()
              return
  end subroutine Release_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetbyKWD_StatementList(keyword, Statements, STR, LINE, SPtr)
  !***  PURPOSE:   to extract the statement in headed by keyword
  !
  !     INPUT:     keyword,  the keyword
  !                Statements, the StatementList
  !
  !     OUTPUT     STR, LINE, the statment and line information
  !                SPtr, optional, the pointer to the Statement containg the keyword
  implicit none
     !--- dummy varioables
     character*(*)::keyword
     type(StatementList),pointer::Statements
     integer::LINE
     character*(*)::STR
     type(StatementList),pointer, optional::SPtr
     !----
     character*64::tWord
     character*256::tSTR

              if(.not. associated(Statements)) then
                 LINE = 0
                 STR  =""
                 if(present(SPtr)) SPtr=>null()
                 return
              end if

              tSTR = adjustl(Statements%this)
              call GetKeyWord("&", tSTR, tWord)
              call UpCase(tWord)
              if(tWord(1:len_trim(tWord)) .eq. keyword(1:len_trim(keyword)) ) then
                 STR  = adjustl(tSTR)
                 LINE = Statements%line
                 if(present(SPtr)) SPtr=>Statements
                 return
              else
                 if(.not.present(SPtr))  then
                    call GetbyKWD_StatementList(keyword, Statements%next, STR, LINE)
                 else
                    call GetbyKWD_StatementList(keyword, Statements%next, STR, LINE, SPtr)
                 end if
              end if

              return

  end subroutine GetbyKWD_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetbyIND0_StatementList(Ind, Statements, STR)
  !***  PURPOSE:   to extract the Indth statement
  !
  !     INPUT:     Ind,        the index of the statements
  !                Statements, the StatementList
  !
  !     OUTPUT     STR,
  !
  implicit none
     !--- dummy varioables
     integer::Ind
     type(StatementList),pointer::Statements
     character*(*)::STR
     !----
     integer::ITH

              if(.not. associated(Statements)) then
                 STR = ""
                 return
              end if
              ITH = IND - 1
              if(ITH .eq. 0) then
                 STR  = adjustl(Statements%this)
                 return
              else
                 call GetbyIND0_StatementList(ITH, Statements%next, STR)
              end if

              return

  end subroutine GetbyIND0_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetbyIND1_StatementList(Ind, Statements, STR, LINE, SPtr)
  !***  PURPOSE:   to extract the statement by index
  !
  !     INPUT:     Ind,  the index
  !                Statements, the StatementList
  !
  !     OUTPUT     STR, LINE, the statment and line information
  !                SPtr, optional, the pointer to the Statement containg the keyword
  implicit none
     !--- dummy varioables
     integer::Ind
     type(StatementList),pointer::Statements
     integer::LINE
     character*(*)::STR
     type(StatementList),pointer, optional::SPtr
     !----
     integer::ITH

              if(.not. associated(Statements)) then
                 LINE = 0
                 STR  = ""
                 if(present(SPtr)) SPtr=>null()
                 return
              end if
              ITH = IND - 1
              if(ITH .eq. 0) then
                 STR  = adjustl(Statements%this)
                 LINE = Statements%line
                 if(present(SPtr)) SPtr=>Statements
                 return
              else
                 if(present(SPtr))  then
                    call GetbyIND1_StatementList(ITH, Statements%next, STR, LINE, SPtr)
                 else
                    call GetbyIND1_StatementList(ITH, Statements%next, STR, LINE)
                 end if
              end if

              return

  end subroutine GetbyIND1_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine FindByKWD_StatementList(keyword, Statements, SPtr)
  !***  PURPOSE:   to find the statement in headed by keyword
  !
  !     INPUT:     keyword,    the keyword
  !                Statements, the StatementList
  !
  !     OUTPUT     SPtr,       the pointer to the Statement containg the keyword
  implicit none
     !--- dummy varioables
     character*(*)::keyword
     type(StatementList),pointer::Statements
     type(StatementList),pointer::SPtr
     !----
     character*64::tWord
     character*256::tSTR

              SPtr=>null()
              if(.not. associated(Statements)) then
                 return
              end if

              tSTR = adjustl(Statements%this)
              call GetKeyWord("&", tSTR, tWord)
              call UpCase(tWord)
              if(tWord(1:len_trim(tWord)) .eq. keyword(1:len_trim(keyword)) ) then
                 SPtr=>Statements
                 return
              else
                 call FindByKWD_StatementList(keyword, Statements%next, SPtr)
              end if

              return

  end subroutine FindByKWD_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive subroutine FindByInd_StatementList(Ind, Statements, SPtr)
  !***  PURPOSE:   to find the statement in headed by keyword
  !
  !     INPUT:     Ind,        the index
  !                Statements, the StatementList
  !
  !     OUTPUT     SPtr,       the pointer to the Indth Statement
  implicit none
     !--- dummy varioables
     integer::Ind
     type(StatementList),pointer::Statements
     type(StatementList),pointer::SPtr
     !----
     integer::ITH

              SPtr=>null()
              if(.not. associated(Statements)) then
                 return
              end if

              ITH = IND-1
              if(ITH .eq.0 ) then
                 SPtr=>Statements
                 return
              else
                 call FindByInd_StatementList(Ind, Statements%next, SPtr)
              end if

              return

  end subroutine FindByInd_StatementList
 !****************************************************************************

 !****************************************************************************
  recursive logical function HasKeyword_StatementList(keyword, Statements) result(res)
  !***  PURPOSE:   to check if the keyword has exist
  !
  !     INPUT:     keyword,  the keyword
  !                Statements, the StatementList
  !
  !     OUTPUT     logical, if the keyword exist
  implicit none
     !--- dummy varioables
     character*(*)::keyword
     type(StatementList),pointer::Statements
     !----
     character*64::tWord
     character*256::tSTR

              if(.not. associated(Statements)) then
                 res = .false.
                 return
              end if

              tSTR = adjustl(Statements%this)
              call GetKeyWord("&", tSTR, tWord)
              call UpCase(tWord)
              if(tWord(1:len_trim(tWord)) .eq. keyword(1:len_trim(keyword)) ) then
                 res = .true.
                 return
              else
                 res =  HasKeyword_StatementList(keyword, Statements%next)
              end if

              return

  end function HasKeyword_StatementList
!****************************************************************************
!****************************************************************************
  recursive integer function Number_StatementList(Statements) result(Num)
  !***  PURPOSE:   to get number of statements on the list.
  !
  !     INPUT:    Statements,    the statement list
  !
  !     OUTPUT    Num
  !
  implicit none
     !--- dummy varioables
     type(StatementList),pointer::Statements
     !----
     integer::N

          Num = 0
          if(.not.associated(Statements)) then
             return
          else
             Num = Number_StatementList(Statements%Next) + 1
          end if
          return
  end function Number_StatementList
 !****************************************************************************

!****************************************************************************
!****************************************************************************
  subroutine Load_InputStatements(fname, Inputs, stag)
  !***  PURPOSE:   to load the strings from a file
  !
  !     INPUT:     fname,  the file name
  !                stag,   optional, a string tag indicating which module this input to be used for
  !
  !     OUTPUT     Inputs
  implicit none
     !--- dummy varioables
     character*(*)::fname
     character*(*), optional::stag
     type(InputStatements)::Inputs
     !---- local variables
     integer::LINE, hFile
     character*256::STR


              Inputs%filename = fname
              if(present(stag)) then
                 Inputs%stag = stag(1:min(len(stag), len(Inputs%stag)))
              else
                 Inputs%stag = ""
              end if

              if(associated(Inputs%stat) ) then
                 call Release_StatementList(Inputs%stat)
              end if

              call AvailableIOUnit(hFile)
              open(UNIT=hFile, file = fname, status='old', err=200)
              LINE = 0
              call  Load_StatementList(hFile, Inputs%stat, STR, LINE)
              close(hFile)
              return

    200     write(*,fmt="(A)") "MDPSCU Error: fail to open file "//fname(1:len_trim(fname))
            write(*,fmt="(A)") "Process to be stopped"
            stop
            return
  end subroutine Load_InputStatements
!****************************************************************************

!****************************************************************************
  subroutine Copy_InputStatements(From, To, StartKWD, EndKWD)
  !***  PURPOSE:   to copy the statement  From to To
  !
  !     INPUT:     From,
  !                StartKWD, optional, the start keyword, copy the list will start from
  !                EndKWD,   optional, the end keyword, copy the list until the keyword found
  !     OUTPUT     To
  implicit none
     !--- dummy varioables
     type(InputStatements), intent(in)          ::From
     type(InputStatements)                      ::To
     character*(*),         intent(in), optional::StartKWD
     character*(*),         intent(in), optional::EndKWD
     !---- local variables
     type(StatementList),  pointer ::SPtr

              call Release_InputStatements(To)
              To%filename = From%filename
              To%stag     = From%stag
              SPtr        => null()
              if(present(StartKWD)) then
                 call FindByKWD_StatementList(StartKWD, From%stat, SPtr)
              else
                 SPtr => From%stat
              end if
              if(present(EndKWD)) then
                 call Copy_StatementList(Sptr, To%stat, EndKWD)
              else
                 call Copy_StatementList(Sptr, To%stat)
              end if
            return
  end subroutine Copy_InputStatements
!****************************************************************************

!****************************************************************************
  subroutine Release_InputStatements(Inputs)
  !***  PURPOSE:   to release the allocted memory
  !
  !     INPUT:     Inputs,
  !     OUTPUT     Inputs
  implicit none
     !--- dummy varioables
     type(InputStatements)::Inputs
     !---- local variables

              Inputs%filename = ""
              Inputs%stag = ""

              if(associated(Inputs%stat) ) then
                 call Release_StatementList(Inputs%stat)
              end if

              return
  end subroutine Release_InputStatements
!****************************************************************************

!****************************************************************************
  subroutine Get_InputStatements(keyword, Inputs, STR, LINE, SPtr)
  !***  PURPOSE:   to extract string  from statementlist for give keyword
  !
  !     INPUT:     keyword,  the keyword to find tghe statement list
  !                Inputs,   the input statement
  !
  !     OUTPUT     STR, LINE, the statment and line information
  !                SPtr, optional, the pointer to the Statement containg the keyword
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)           ::keyword
     type(InputStatements)               ::Inputs
     character*(*)                       ::STR
     integer                             ::LINE
     type(StatementList),pointer,optional::SPtr
     !---- local variables

              STR = ""
              LINE = 0
              if(present(SPtr)) then
                 call GetbyKWD_StatementList(keyword, Inputs%stat, STR, LINE, SPtr)
              else
                 call GetbyKWD_StatementList(keyword, Inputs%stat, STR, LINE)
              end if
            return
  end subroutine Get_InputStatements
  !****************************************************************************

  !****************************************************************************
  logical function HasKeyword_InputStatements(keyword, Inputs) result(res)
  !***  PURPOSE:   to check if the keyword has exist
  !
  !     INPUT:     keyword,  the keyword
  !                Inputs,   the input statements
  !
  !     OUTPUT     logical, if the keyword exist
  implicit none
     !--- dummy varioables
     character*(*), intent(in)  ::keyword
     type(InputStatements)      ::Inputs

     !----
     character*64::tWord
     character*256::tSTR

              res = HasKeyword_StatementList(keyword, Inputs%stat)
              return

  end function HasKeyword_InputStatements
!****************************************************************************

  end module MSM_TYPEDEF_InputPaser
  !****************************************************************************************
