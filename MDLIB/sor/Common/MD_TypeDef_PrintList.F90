  module MD_TYPEDEF_PrintList
  !***  DESCRIPTION:
  !     This module provides a typedef that collects the pointers of external proceddures
  !     that are  for printing out simulation control parameters
  !
  !     SEE ALSO:
  !             MD_SimBoxArray_ToolShell.f90
  !             MD_SimBoxArray_ToolShell_GPU.f90
  !
  !
  !    Written by HOU Qing, May, 2018
  !                  2019-09- (HOU Qing),
  !                          Add archive-restore process list
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl

   implicit none

  !--- interface to the external routine -------------------
   private::EXPRINTPROC
   abstract interface
       SUBROUTINE EXPRINTPROC(hFile, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       integer          :: hFile
       type(SimMDBox)   :: SimBox
       type(SimMDCtrl)  :: CtrlParam
       END SUBROUTINE EXPRINTPROC
   end interface

   private::EXARCHIVEPROC
   abstract interface
       SUBROUTINE EXARCHIVEPROC(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE EXARCHIVEPROC
   end interface   
  !--- END INTERFACE --------------------------------

   type PrintProcedures
        procedure(EXPRINTPROC), pointer, nopass::pExPrintProc=>null()
        type(PrintProcedures),  pointer::Next
   end type PrintProcedures

   type ArchiveProcedures
        character*64 :: Tag
        procedure(EXARCHIVEPROC), pointer, nopass::pExArchiveProc=>null()
        type(ArchiveProcedures),  pointer::Next
   end type ArchiveProcedures

  !--- END TYPEDEFINE --------------------------------

   type(PrintProcedures),   pointer, private::m_pParamPrinter=>null()
   type(ArchiveProcedures), pointer, private::m_pArchiver=>null()
   type(ArchiveProcedures), pointer, private::m_pRestore=>null()

  !--- interface to the external routine -------------------
  !---------------------------------------------------------
      private:: Add_ArchiveProcess0,   &
                Add_ArchiveProcess1,   &
                Add_ArchiveProcess2
      public :: Add_ArchiveProcess          
      interface Add_ArchiveProcess
          module procedure Add_ArchiveProcess0
          module procedure Add_ArchiveProcess1
          module procedure Add_ArchiveProcess2
      end interface Add_ArchiveProcess

      private:: Add_RestoreProcess0,   &
                Add_RestoreProcess1
      public :: Add_RestoreProcess                
      interface Add_RestoreProcess
          module procedure Add_RestoreProcess0
          module procedure Add_RestoreProcess1
          module procedure Add_ArchiveProcess2
      end interface Add_RestoreProcess


  !---------------------------------------------------------
      private:: Add_PrintProcess0,   &
                Add_PrintProcess1
      public :: Add_PrintProcess                 
      interface Add_PrintProcess
          module procedure Add_PrintProcess0
          module procedure Add_PrintProcess1
      end interface Add_PrintProcess

  !---------------------------------------------------------
      private:: Clear_PrintProcess0,   &
                Clear_PrintProcess1
      interface Clear_PrintProcess
          module procedure Clear_PrintProcess0
          module procedure Clear_PrintProcess1
      end interface Clear_PrintProcess
      
  !---------------------------------------------------------
      private:: Do_PrintProcess0,   &
                Do_PrintProcess1
      public :: Do_PrintProcess          
      interface Do_PrintProcess
          module procedure Do_PrintProcess0
          module procedure Do_PrintProcess1
      end interface Do_PrintProcess

  !---------------------------------------------------------
      private:: Do_ArchiveProcess0
      public :: Do_ArchiveProcess        

  !---------------------------------------------------------
      private:: Do_RestoreProcess0
      public::  Do_RestoreProcess
      
  !---------------------------------------------------------
      private:: Ind_ArchiveProcess0, &
                Ind_ArchiveProcess1
      public::  Ind_ArchiveProcess
      interface Ind_ArchiveProcess
          module procedure Ind_ArchiveProcess0
          module procedure Ind_ArchiveProcess1
      end interface Ind_ArchiveProcess

      public::  Ind_RestoreProcess
      interface Ind_RestoreProcess
          module procedure Ind_ArchiveProcess0
          module procedure Ind_RestoreProcess1
      end interface Ind_RestoreProcess

  !---------------------------------------------------------
      private:: ReplaceByInd_ArchiveProcess0,    &
                ReplaceByInd_ArchiveProcess1,    &
                ReplaceByTag_ArchiveProcess0,    &
                ReplaceByTag_ArchiveProcess1
      public::  Replace_ArchiveProcess
      interface Replace_ArchiveProcess
          module procedure ReplaceByInd_ArchiveProcess0
          module procedure ReplaceByInd_ArchiveProcess1
          module procedure ReplaceByTag_ArchiveProcess0
          module procedure ReplaceByTag_ArchiveProcess1
      end interface Replace_ArchiveProcess

      public::  Replace_RestoreProcess
      interface Replace_RestoreProcess
          module procedure ReplaceByInd_ArchiveProcess0
          module procedure ReplaceByTag_ArchiveProcess0
          module procedure ReplaceByInd_RestoreProcess1
          module procedure ReplaceByTag_RestoreProcess1
      end interface Replace_RestoreProcess

  contains

  !****************************************************************************************
  recursive subroutine Add_PrintProcess0(LIST, ExPrinter)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
       SUBROUTINE ExPrinter(hFile, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       integer         :: hFile
       type(SimMDBox)  :: SimBox
       type(SimMDCtrl) :: CtrlParam
       END SUBROUTINE ExPrinter
   end interface
  !--- END INTERFACE --------------------------------
   type(PrintProcedures),   pointer ::LIST
   !--- Local

           if(.not.associated(LIST)) then
               allocate(LIST)
               LIST%pExPrintProc =>ExPrinter
               LIST%Next => null()
           else
               call Add_PrintProcess(LIST%Next, ExPrinter)
           end if
           return
   end subroutine Add_PrintProcess0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Add_PrintProcess1(ExPrinter)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
       SUBROUTINE ExPrinter(hFile, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       integer         :: hFile
       type(SimMDBox)  :: SimBox
       type(SimMDCtrl) :: CtrlParam
       END SUBROUTINE ExPrinter
   end interface
  !--- END INTERFACE --------------------------------
   !--- Local

           call Add_PrintProcess0(m_pParamPrinter, ExPrinter)
           return
   end subroutine Add_PrintProcess1
  !****************************************************************************************

  !****************************************************************************************
  subroutine Ind_ArchiveProcess0(LIST, Tag, Ind) 
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    type(ArchiveProcedures), pointer    ::LIST
    character*(*),           intent(in) ::Tag
    integer,                 intent(out)::Ind
    !--- Local
     type(ArchiveProcedures), pointer::pTemp
     integer::AT

           AT  = 0
           pTemp => List 
           do while(associated(pTemp))
              AT = AT + 1
              if(LIST%Tag .eq. Tag) then
                  Ind = AT 
                  return
              end if
              pTemp => pTemp%Next
           end do
           Ind = 0
           return
   end subroutine Ind_ArchiveProcess0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Ind_ArchiveProcess1(Tag, Ind) 
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    character*(*),           intent(in) ::Tag
    integer,                 intent(out)::Ind
    !--- Local

           call Ind_ArchiveProcess0(m_pArchiver, Tag, Ind)
           return
   end subroutine Ind_ArchiveProcess1
  !****************************************************************************************   

  !****************************************************************************************
  subroutine Ind_RestoreProcess1(Tag, Ind) 
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    character*(*),           intent(in) ::Tag
    integer,                 intent(out)::Ind
    !--- Local

           call Ind_ArchiveProcess0(m_pRestore, Tag, Ind)
           return
   end subroutine Ind_RestoreProcess1
  !****************************************************************************************

  !****************************************************************************************
  subroutine ReplaceByInd_ArchiveProcess0(LIST, Ind, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    type(ArchiveProcedures), pointer ::LIST
    integer::Ind
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local
    type(ArchiveProcedures), pointer::pTemp
    integer::AT

          AT  = 0
          pTemp => List 
          do while(associated(pTemp))
             AT = AT + 1
             if(AT .eq. Ind) then
                pTemp%pExArchiveProc =>ExArchive
                return
             end if
             pTemp => pTemp%Next
          end do
          return
   end subroutine ReplaceByInd_ArchiveProcess0
  !****************************************************************************************

  !****************************************************************************************
  subroutine ReplaceByInd_ArchiveProcess1(Ind, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    type(ArchiveProcedures), pointer ::LIST
    integer::Ind
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local

          call ReplaceByInd_ArchiveProcess0(m_pArchiver, Ind, ExArchive)       
          return
   end subroutine ReplaceByInd_ArchiveProcess1
  !****************************************************************************************  

  !****************************************************************************************
  subroutine ReplaceByInd_RestoreProcess1(Ind, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    type(ArchiveProcedures), pointer ::LIST
    integer::Ind
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local

          call ReplaceByInd_ArchiveProcess0(m_pRestore, Ind, ExArchive)       
          return
   end subroutine ReplaceByInd_RestoreProcess1
  !****************************************************************************************
   
  !****************************************************************************************
  subroutine ReplaceByTag_ArchiveProcess0(LIST, Tag, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    type(ArchiveProcedures), pointer ::LIST
    character*(*)::Tag
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local
    integer::AT

          AT  = 0
          call Ind_ArchiveProcess0(LIST, Tag, AT)  
          if(AT .gt. 0) then
              call ReplaceByInd_ArchiveProcess0(LIST, AT, ExArchive)
          end if
          return
   end subroutine ReplaceByTag_ArchiveProcess0
  !****************************************************************************************   

  !****************************************************************************************
  subroutine ReplaceByTag_ArchiveProcess1(Tag, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    character*(*)::Tag
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local

          call ReplaceByTag_ArchiveProcess0(m_pArchiver, Tag, ExArchive)  
          return
   end subroutine ReplaceByTag_ArchiveProcess1
  !****************************************************************************************   

  !****************************************************************************************
  subroutine ReplaceByTag_RestoreProcess1(Tag, ExRestore)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    character*(*)::Tag
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExRestore(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExRestore
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local

          call ReplaceByTag_ArchiveProcess0(m_pRestore, Tag, ExRestore)  
          return
   end subroutine ReplaceByTag_RestoreProcess1
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine Add_ArchiveProcess0(LIST, Tag, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    type(ArchiveProcedures), pointer ::LIST
    character*(*)::Tag
    !--- interface to the external routine -------------------
    interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
    end interface
    !--- END INTERFACE --------------------------------
    !--- Local
    integer::IND

           call Ind_ArchiveProcess0(LIST, Tag, IND)  
           if(IND .gt. 0) then
              call ReplaceByInd_ArchiveProcess0(LIST, IND, ExArchive)
           else   
              if(.not.associated(LIST)) then
                 allocate(LIST)
                 LIST%Tag = Tag
                 LIST%pExArchiveProc =>ExArchive
                 LIST%Next => null()
              else
                 call Add_ArchiveProcess0(LIST%Next, Tag, ExArchive)
              end if
           end if    
           return
   end subroutine Add_ArchiveProcess0
  !****************************************************************************************
   
  !****************************************************************************************
  subroutine Add_ArchiveProcess1(Tag, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
     character*(*)::Tag
     !--- interface to the external routine -------------------
     interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
     end interface
     !--- END INTERFACE --------------------------------
     !--- Local
           call Add_ArchiveProcess0(m_pArchiver, Tag, ExArchive)
           return
   end subroutine Add_ArchiveProcess1
  !****************************************************************************************   

  !****************************************************************************************
  subroutine Add_ArchiveProcess2(Tag, ExArchive, ExRestore)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
     character*(*)::Tag
     !--- interface to the external routine -------------------
     interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
     end interface
     interface
       SUBROUTINE ExRestore(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExRestore
     end interface
     !--- END INTERFACE --------------------------------
     !--- Local
     integer::IA, IR

           call Add_ArchiveProcess1(Tag, ExArchive)
           call Add_RestoreProcess1(Tag, ExRestore)
           !--- to check if the Archive and Restore are in the same sequence
           call Ind_ArchiveProcess1(Tag, IA) 
           call Ind_RestoreProcess1(Tag, IR) 
           if(IA .ne. IR) then
              write(*,*) "MCPSCU Error: archive and restore priocess is not consistent "
              stop
           end if 
           return
   end subroutine Add_ArchiveProcess2
  !****************************************************************************************

  !****************************************************************************************
  subroutine Add_RestoreProcess0(LIST, Tag, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    character*(*)::Tag
    !--- interface to the external routine -------------------
     interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
     end interface
    !--- END INTERFACE --------------------------------
     type(ArchiveProcedures), pointer ::LIST
    !--- Local
           call Add_ArchiveProcess0(LIST, Tag, ExArchive)
           return
   end subroutine Add_RestoreProcess0
  !****************************************************************************************   


  !****************************************************************************************
  subroutine Add_RestoreProcess1(Tag, ExArchive)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
     character*(*)::Tag
     !--- interface to the external routine -------------------
     interface
       SUBROUTINE ExArchive(hFile)
       implicit none
       integer             :: hFile
       END SUBROUTINE ExArchive
     end interface
     !--- END INTERFACE --------------------------------
     !--- Local
           call Add_RestoreProcess0(m_pRestore, Tag, ExArchive)
           return
   end subroutine Add_RestoreProcess1
  !****************************************************************************************   

  !****************************************************************************************
  recursive subroutine Clear_PrintProcess0(LIST)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
   type(PrintProcedures),   pointer ::LIST
   !--- Local

           if(.not.associated(LIST)) return

           call Clear_PrintProcess(LIST%Next)
           deallocate(LIST)
           LIST =>null()
           return
   end subroutine Clear_PrintProcess0
  !****************************************************************************************

  !*************************************************************************************
  subroutine Clear_PrintProcess1()
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
   !--- Local
           call Clear_PrintProcess0(m_pParamPrinter)
           return
   end subroutine Clear_PrintProcess1
  !**************************************************************************************** 
   
  !****************************************************************************************
  recursive subroutine Clear_ArchiveProcess0(LIST)
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
   type(ArchiveProcedures),   pointer ::LIST
    !--- Local
 
            if(.not.associated(LIST)) return
 
            call Clear_ArchiveProcess0(LIST%Next)
            deallocate(LIST)
            LIST =>null()
            return
  end subroutine Clear_ArchiveProcess0
  !****************************************************************************************
 
  !*************************************************************************************
  subroutine Clear_ArchiveProcess1()
  !***  DESCRIPTION: to add an external print procedure to the list
  !
  !     INPUT:  pExPrinter,  the subroutine provided by user for printing control parameters
  !     OUTPUT: LIST
  !
  implicit none
    !--- Local
            call Clear_ArchiveProcess0(m_pArchiver)
            call Clear_ArchiveProcess0(m_pRestore)
            return
  end subroutine Clear_ArchiveProcess1
  !****************************************************************************************   

  !****************************************************************************************
  recursive subroutine Do_PrintProcess0(hFile, LIST, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  integer                        ::hFile
  type(PrintProcedures), pointer ::LIST
  type(SimMDBox)                 ::SimBox
  type(SimMDCtrl)                ::CtrlParam


      if(.not.associated(LIST)) return

         call LIST%pExPrintProc(hFile, SimBox, CtrlParam)
         if(associated(LIST%next)) then
             call Do_PrintProcess(hFile, LIST%next, SimBox, CtrlParam)
         end if
      return
  end subroutine Do_PrintProcess0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_PrintProcess1(hFile, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  integer         ::hFile
  type(SimMDBox)  ::SimBox
  type(SimMDCtrl) ::CtrlParam

      call Do_PrintProcess0(hFile, m_pParamPrinter, SimBox, CtrlParam)
      return
  end subroutine Do_PrintProcess1
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine Do_ArchiveProcess0(hFile, LIST)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  integer                        ::hFile
  type(ArchiveProcedures), pointer ::LIST

      if(.not.associated(LIST)) return
         write(hFile) LIST%Tag 
         call LIST%pExArchiveProc(hFile)
         if(associated(LIST%next)) then
             call Do_ArchiveProcess0(hFile, LIST%next)
         end if
      return
  end subroutine Do_ArchiveProcess0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_ArchiveProcess(hFile)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  integer                        ::hFile
      call Do_ArchiveProcess0(hFile, m_pArchiver) 
      return
  end subroutine Do_ArchiveProcess
  !****************************************************************************************
 
  !****************************************************************************************
  recursive subroutine Do_RestoreProcess0(hFile, LIST)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
     integer                          ::hFile
     type(ArchiveProcedures), pointer ::LIST
     !-- local
     character*64::Tag
     integer::AtWarning = CP_WARNING_PAUSE
     save AtWarning 
   
      if(.not.associated(LIST)) return
         read(hFile, end=100) Tag 
         if(Tag(1:len_trim(Tag)) .ne. LIST%Tag(1:len_trim(LIST%Tag)) ) then
            write(*,*) "MDPSCU Error: Tag is no consistent in restoring the data from file"
            stop
         end if   
         call LIST%pExArchiveProc(hFile)
         if(associated(LIST%next)) then
             call Do_RestoreProcess0(hFile, LIST%next)
         end if
         return

    100  continue
         write(*,*) "MDPSCU Warning: fail to resore data for ", LIST%Tag(1:len_trim(LIST%Tag))
         call ONWARNING(AtWarning)
      return
  end subroutine Do_RestoreProcess0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_RestoreProcess(hFile)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  integer                        ::hFile
  character*256::EXENAME

      call Do_RestoreProcess0(hFile, m_pRestore) 
      return
  end subroutine Do_RestoreProcess
  !****************************************************************************************

 end module MD_TYPEDEF_PrintList
 !****************************************************************************************
