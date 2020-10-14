  module try
  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none

  contains

  subroutine MyMask(SimBox, CtrlParam, MASK)
  use MD_CONSTANTS
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
      implicit none
           !--- dummy vaiables
           type(SimMDBox),                 intent(in)::SimBox
           type(SimMDCtrl),                intent(in)::CtrlParam
           integer,dimension(:),allocatable          ::MASK
           integer::I

           do I=1, SimBox%NPRT
              if( SimBox%ITYP(I) .eq. 1) then
                  MASK(I) = 1
              else
                  MASK(I) = 0
              end if
          end do
   end subroutine MyMask
  end module try


 !****************************************************************************************
  Program SUBSTRATE_2014_Main
  use MD_SimBoxArray_AppShell_16_GPU
  use MD_Method_ParRep_GPU, only:SetCompareRoutine
  !---- If user-define potential to be used, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
  !     use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
  !     use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table
  use try
  implicit none
  integer::numprocs=1, procid=0

       !call SetCompareRoutine(COMPAREMASK=MyMask)
       call APPSHELL_Main(numprocs,procid)


       STOP "--------- 2014_GPU"
       stop
  end  program SUBSTRATE_2014_Main
