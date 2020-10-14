  module MD_LocalTempMethod_GPU
  !**** DESCRIPTION: the is module encupsal a lnumber of methods of control
  !                  local temperature.
  !
  !****
  !     DEPENDENCE:  MD_EP_Coupling_GPU.F90
  !                  MD_STPLib_Register.F90
  !
  !****
  !     HISTORY:     2013-03-15(HOU Qing):
  !                          Create the first version
  !

  use CUDAFOR
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl

  implicit none
  contains

  !**************************************************************************
  subroutine Initialize_LocalTempMethod_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to initialze this module by loading calculation parameters
  !
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  use MD_EP_Coupling_GPU, only:Initialize_EPCMOD_DEV
  use MD_ST_Coupling_GPU, only:Initialize_STMOD_DEV
  use MD_STPLib_Register, only:m_STPTable, Register_STPLib, Release_STPTable
  implicit none
      !----   DUMMY Variables
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in), target::CtrlParam
      !--- local variable
       type(SimMDCtrl), pointer::nextp, next
       integer::I, flag

         !--- initialize the EPC model
         call Clear_LocalTempMethod_DEV()
         call Initialize_EPCMOD_DEV(SimBox, CtrlParam)

         !--- initialize the Stopping model
         call Release_STPTable(m_STPTable)

         !$$--- to check if we need stopping
          nextp=>CtrlParam
          do while(.TRUE.)
             do I = 1, SimBox%NGROUP
                flag = iand(nextp%LT_CTRL(I)%METH, CP_TICTRL_METH_ST)
                if(flag .eq. CP_TICTRL_METH_ST) exit
             end do
             if(flag .eq. CP_TICTRL_METH_ST) exit

             call GetNext_SimMDCtrl(nextp, 1, next)
             if(.not. associated(next)) exit
             nextp=>next
          end do

          if(flag .eq. CP_TICTRL_METH_ST) then
             call Register_STPLib(SimBox, CtrlParam, m_STPTable)
             call Initialize_STMOD_DEV(SimBox, CtrlParam, m_STPTable)
          end if

        return
   end subroutine Initialize_LocalTempMethod_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Clear_LocalTempMethod_DEV(SimBox)
  !***  PURPOSE: to clear the memories allocated on calling initialize
  !
  !     INPUT:
  !     OUTPUT:
  !
  use MD_EP_Coupling_GPU, only:Clear_EPCMOD_DEV
  use MD_ST_Coupling_GPU, only:Clear_STMOD_DEV
  implicit none
      !----   DUMMY Variables
      type(SimMDBox), dimension(:), optional::SimBox

         call Clear_EPCMOD_DEV()
         if(.not. present(SimBox)) then
            call Clear_STMOD_DEV()
         else
             call Clear_STMOD_DEV(SimBox)
         end if
        return
   end subroutine Clear_LocalTempMethod_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Do_ResetParam_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to initialze this module by loading calculation parameters
  !
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  use MD_EP_Coupling_GPU, only:Reset_EPCMOD_DEV
  use MD_ST_Coupling_GPU, only:Reset_STMOD_DEV
  implicit none
      !----   DUMMY Variables
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam
      !--- local variable

             call Reset_EPCMOD_DEV(SimBox, CtrlParam)
             call Reset_STMOD_DEV(SimBox, CtrlParam)

          return
  end subroutine Do_ResetParam_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Do_EPCForce_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to initialze this module by loading calculation parameters
  !
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  use MD_EP_Coupling_GPU, only:Do_EPCMOD_DEV
  use MD_ST_Coupling_GPU, only:Do_STMOD_DEV
  implicit none
      !----   DUMMY Variables
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
      !--- local variable
             call Do_EPCMOD_DEV(SimBox, CtrlParam)
             call Do_STMOD_DEV(SimBox, CtrlParam)

          return
  end subroutine Do_EPCForce_DEV
  !**************************************************************************

 !**************************************************************************
 end module MD_LocalTempMethod_GPU
