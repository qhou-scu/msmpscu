  module MD_TYPEDEF_RecordList
  !***  DESCRIPTION:
  !     This module provides a typedef that collects the pointers of external proceddures
  !     that be usd in simulation or analysis
  !
  !     SEE ALSO:
  !             MD_SimBoxArray_ToolShell.f90
  !             MD_EM_TB_ForceTable_SHELL_CPU.f90
  !             MD_SimBoxArray_ToolShell_GPU.f90
  !             MD_EM_TB_ForceTable_SHELL_GPU.f90
  !
  !
  !    Written by HOU Qing, May, 2014
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_TYPEDEF_RecordStamp

   implicit none

  !--- interface to the external routine -------------------
   abstract interface
       subroutine AFTERONESTEP(ITime, SimBox, CtrlParam, Statu)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       use MD_TYPEDEF_RecordStamp
       implicit none
       integer                     :: ITime
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       integer                     :: Statu
       end subroutine AFTERONESTEP
   end interface

   abstract interface
       subroutine AFTERONETEST(SimBox0, Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)              :: SimBox0
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       end subroutine AFTERONETEST
   end interface

   abstract interface
       subroutine AFTERRECORD(SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)   :: SimBox
       type(SimMDCtrl)  :: CtrlParam
       end subroutine AFTERRECORD
   end interface

   abstract interface
       subroutine BEFOREONETEST(JBOX, SimBox0, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       integer                     ::JBOX
       type(SimMDBox)              :: SimBox0
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       end subroutine BEFOREONETEST
   end interface

   abstract interface
       subroutine PRERECORD(SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)  :: SimBox
       type(SimMDCtrl) :: CtrlParam
       end subroutine PRERECORD
   end interface

   abstract interface
       subroutine RECORDPROC(Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       use MD_TYPEDEF_RecordStamp
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       end subroutine RECORDPROC
   end interface

   abstract interface
       subroutine RECORDSTATU(STATU)
       implicit none
       integer::STATU
       end subroutine RECORDSTATU
   end interface

   abstract interface
       subroutine RESTARTREC(Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       use MD_TYPEDEF_RecordStamp
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       end subroutine RESTARTREC
   end interface

   abstract interface
       subroutine SAVERECSTATU(Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       use MD_TYPEDEF_RecordStamp
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       end subroutine SAVERECSTATU
   end interface


  !--- end INTERFACE --------------------------------
      private::RecordProcedures
      type   ::RecordProcedures
           !--- NOTE: if the name of a procedure pointer is too long,
           !          there may be compiling errors on linux.
           !          Note by Hou Qing 2017-11-25

           procedure(AFTERONESTEP),          pointer, nopass::pAftonestep=>null()
           procedure(AFTERONETEST),          pointer, nopass::pAftonetest=>null()
           procedure(AFTERRECORD),           pointer, nopass::pAftrecord=>null()
           procedure(BEFOREONETEST),         pointer, nopass::pBefonetest=>null()
           procedure(PRERECORD),             pointer, nopass::pPrerecord=>null()
           procedure(RECORDPROC),            pointer, nopass::pRecordproc=>null()
           procedure(RECORDSTATU),           pointer, nopass::pRecordStatu=>null()
           procedure(RESTARTREC),            pointer, nopass::pRestartRec=>null()
           procedure(SAVERECSTATU),          pointer, nopass::pSaveRecStat=>null()
      end type RecordProcedures

      type::RecordProcedureList
           type(RecordProcedures),    pointer::this=>null()
           type(RecordProcedureList), pointer::next=>null()
      end type RecordProcedureList
      !**** interface to external calling stack
      !----------------------------------------------------
      public:: Add_RecordProcedures
      !----------------------------------------------------
      public:: DoAfterRecord_List
      !----------------------------------------------------
      public:: DoAfterOneStep_List
      !----------------------------------------------------
      public:: DoAfterOneTest_List
      !----------------------------------------------------
      public:: DoBeforeOneTest_List
      !----------------------------------------------------
      public:: DoPrerecord_List
      !----------------------------------------------------
      public:: DoRecord_List
      !----------------------------------------------------
      public:: DoRestartRecord_List
      !----------------------------------------------------
      public:: DoSaveRecStatu_List
      !----------------------------------------------------
      public:: RecordStatu_List
      !----------------------------------------------------

  contains

  !****************************************************************************************
  recursive subroutine Add_RecordProcedures(LIST, pAFTONESTEP, pAFTONETEST, pAFTRECORD, pBEFONETEST, pPRERECORD, pRECORDPROC, pRECSTATU, pRESTARTREC, pSAVERECSTATU)
  !***  DESCRIPTION: to set prerecord procedure of RecordProcedures P.
  !                  the prerecord procedure to be run before the simulation
  !                  started.
  !     INPUT: PRERECORD,  the subroutine provided by user for pre-recording
  !
  implicit none
  !--- END INTERFACE --------------------------------
   type(RecordProcedureList)::LIST
   procedure(AFTERONESTEP),      pointer, intent(in)::pAftonestep
   procedure(AFTERONETEST),      pointer, intent(in)::pAftonetest
   procedure(AFTERRECORD),       pointer, intent(in)::pAftrecord
   procedure(BEFOREONETEST),     pointer, intent(in)::pBefonetest
   procedure(PRERECORD),         pointer, intent(in)::pPrerecord
   procedure(RECORDPROC),        pointer, intent(in)::pRecordproc
   procedure(RECORDSTATU),       pointer, intent(in)::pRecStatu
   procedure(RESTARTREC),        pointer, intent(in)::pRestartRec
   procedure(SAVERECSTATU),      pointer, intent(in)::pSaveRecStatu
   !--- Local

           if(.not.associated(LIST%this)) then
               allocate(LIST%this)
               LIST%this%pAftonestep   =>pAftonestep
               LIST%this%pAftonetest   =>pAftonetest
               LIST%this%pAftrecord    =>pAftrecord
               LIST%this%pBefonetest   =>pBefonetest
               LIST%this%pPrerecord    =>pPrerecord
               LIST%this%pRecordproc   =>pRecordproc
               LIST%this%pRecordStatu  =>pRecStatu
               LIST%this%pRestartRec   =>pRestartrec
               LIST%this%pSaveRecStat  =>pSaveRecStatu
               LIST%next => null()
           else
               if(.not.associated(LIST%next)) then
                 allocate(LIST%next)
                 LIST%next%this => null()
               end if
               call Add_RecordProcedures(LIST%next,                     &
                                         pAFTONESTEP   = pAftonestep,   &
                                         pAFTONETEST   = pAftonetest,   &
                                         pAFTRECORD    = pAftrecord,    &
                                         pBEFONETEST   = pBefonetest,   &
                                         pPRERECORD    = pPrerecord,    &
                                         pRECORDPROC   = pRecordproc,   &
                                         pRESTARTREC   = pRestartrec,   &
                                         pRECSTATU     = pRecStatu,     &
                                         pSAVERECSTATU = pSaveRecStatu   )


           end if
           return
   end subroutine Add_RecordProcedures
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine DoPrerecord_List(LIST, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  type(RecordProcedureList)::LIST
  type(SimMDBox)           ::SimBox
  type(SimMDCtrl)          ::CtrlParam

      if(.not.associated(LIST%this)) return

      if(associated(LIST%this%pPrerecord)) call LIST%this%pPrerecord(SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoPrerecord_List(LIST%next, SimBox, CtrlParam)
      end if
      return
  end subroutine DoPrerecord_List
 !****************************************************************************************

 !****************************************************************************************
  recursive subroutine DoAfterRecord_List(LIST, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  type(RecordProcedureList)::LIST
  type(SimMDBox)           ::SimBox
  type(SimMDCtrl)          ::CtrlParam

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pAftRecord)) call List%this%pAftRecord(SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoAfterRecord_List(LIST%next, SimBox, CtrlParam)
      end if
      return
  end subroutine DoAfterRecord_List
 !****************************************************************************************

 !****************************************************************************************
  recursive subroutine DoBeforeOneTest_List(LIST, JBOX, SimBox0, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  type(RecordProcedureList)   ::LIST
  type(SimMDBox)              ::SimBox0
  type(SimMDBox), dimension(:)::SimBox
  type(SimMDCtrl)             ::CtrlParam
  integer                     ::JBOX

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pBefonetest)) call LIST%this%pBefonetest(JBOX, SimBox0, SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoBeforeOneTest_List(LIST%next, JBOX, SimBox0, SimBox, CtrlParam)
      end if
      return
  end subroutine DoBeforeOneTest_List
 !****************************************************************************************

 !****************************************************************************************
  recursive subroutine DoAfterOneTest_List(LIST, SimBox0, Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
   type(RecordProcedureList)    ::LIST
   type(SimMDBox)               ::SimBox0
   type(MDRecordStamp)          ::Stamp
   type(SimMDBox), dimension(:) ::SimBox
   type(SimMDCtrl)              ::CtrlParam

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pAftonetest)) call LIST%this%pAftonetest(SimBox0, Stamp, SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoAfterOneTest_List(LIST%next, SimBox0, Stamp, SimBox, CtrlParam)
      end if

      return
  end subroutine DoAfterOneTest_List
 !****************************************************************************************

  !****************************************************************************************
  recursive subroutine DoAfterOneStep_List(LIST, ITime, SimBox, CtrlParam, Statu)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
   type(RecordProcedureList)     :: LIST
   integer                       :: ITime
   type(SimMDBox), dimension(:)  :: SimBox
   type(SimMDCtrl)               :: CtrlParam
   integer                       :: Statu

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pAftonestep)) call LIST%this%pAftonestep(ITime, SimBox, CtrlParam, Statu)

      if(associated(LIST%next)) then
         call DoAfterOneStep_List(LIST%next, ITime, SimBox, CtrlParam, Statu)
      end if
      return
  end subroutine DoAfterOneStep_List
 !****************************************************************************************


 !****************************************************************************************
  recursive subroutine DoRecord_List(LIST, Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
   type(RecordProcedureList)     :: LIST
   type(MDRecordStamp)           :: Stamp
   type(SimMDBox), dimension(:)  :: SimBox
   type(SimMDCtrl)               :: CtrlParam

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pRecordproc)) call LIST%this%pRecordproc(Stamp, SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoRecord_List(LIST%next,Stamp, SimBox, CtrlParam)
      end if
      return
  end subroutine DoRecord_List
 !****************************************************************************************

 !****************************************************************************************
  recursive subroutine DoRestartRecord_List(LIST, Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
   type(RecordProcedureList),   intent(in):: LIST
   type(MDRecordStamp),         intent(in):: Stamp
   type(SimMDBox), dimension(:)           :: SimBox
   type(SimMDCtrl)                        :: CtrlParam

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pRestartRec)) call LIST%this%pRestartRec(Stamp, SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoRestartRecord_List(LIST%next,Stamp, SimBox, CtrlParam)
      end if
      return
  end subroutine DoRestartRecord_List
 !****************************************************************************************

 !****************************************************************************************
  recursive subroutine DoSaveRecStatu_List(LIST, Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
   type(RecordProcedureList),  intent(in):: LIST
   type(MDRecordStamp),        intent(in):: Stamp
   type(SimMDBox),dimension(:)           :: SimBox
   type(SimMDCtrl)                       :: CtrlParam

      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pSaveRecStat)) call LIST%this%pSaveRecStat(Stamp, SimBox, CtrlParam)

      if(associated(LIST%next)) then
         call DoSaveRecStatu_List(LIST%next,Stamp, SimBox, CtrlParam)
      end if
      return
  end subroutine DoSaveRecStatu_List
 !****************************************************************************************

 !****************************************************************************************
  recursive subroutine RecordStatu_List(LIST, STATU)
  !***  DESCRIPTION:
  !     INPUT:
  implicit none
  type(RecordProcedureList)::LIST
  integer::STATU, TS


      STATU = 0
      if(.not.associated(LIST%this)) return
      if(associated(LIST%this%pRecordStatu)) then
         call LIST%this%pRecordStatu(STATU)
      end if

      if(associated(LIST%next)) then
         call RecordStatu_List(LIST%next,TS)
         STATU = IOR(STATU, TS)
      end if
      return
  end subroutine RecordStatu_List
  !****************************************************************************************


  end module MD_TYPEDEF_RecordList
  !****************************************************************************************
