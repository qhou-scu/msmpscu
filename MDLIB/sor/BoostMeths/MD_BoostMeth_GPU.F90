  module MD_BoostMeth_GPU
  !**** DESCRIPTION: This module provides routines to spring-boost force
  !                  potential surface along a given direction. This module coulbe be used
  !                  in steepest descent or CG of search minima on energy surface         
  !
  !                  ______________________________________________________
  !**** HISTORY:     2019-01-16(HOU Qing),
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_Globle_Variables_GPU
  implicit none
  !**********************
 

  contains
  !****************************************************************************************
  subroutine Reset_BoostMeth()
  !***  PORPOSE:   to clear the memory allocated in  SPBState_List
  !
  !     INPUT:     theSPBState
  !     OUTPUT:    theSPBState
  !
  use MD_SpringBoost_GPU, only:Reset_SpringForce=>Reset_BoostForce
  implicit none
  !----  DUMMY Variables
  !----  Local variables
            
            call Reset_SpringForce()
      return 
  end subroutine Reset_BoostMeth
  !**************************************************************************************** 

  !****************************************************************************************
  subroutine New_BoostForce(SimBox, CtrlParam, NP, PId0, XP0)
   !***  PORPOSE:   to to add a boost state to the boost list m_BoostStates
   !
   !     INPUT:     SimBox,    
   !                CtrlParam, 
   !                NP,         the number of particles to be boosted  
   !                PId0,       the IDs of particle to be boosted
   !
   !     OUTPUT:    to add a boost state to the boost list:   the boost state list
   ! 
   use MD_SpringBoost_GPU, only:Add_BoostForce_Spring=>Add_BoostForce
   implicit none
   !----   DUMMY Variables
      type(SimMDBox),  intent(in), dimension(:)  ::SimBox
      type(SimMDCtrl), intent(in)                ::CtrlParam
      integer,         intent(in)                ::NP   
      integer,         intent(in), dimension(:)  ::PID0
      real(KINDDF),    intent(in), dimension(:,:)::XP0
   !----  Local variables
      real(KINDDF), dimension(:,:), allocatable::POS0

            allocate(POS0(NP,3))
            !IP = 0

            !do IP=1, NP
            !   POS0(IP,:) 
            !call Add_BoostForce_Spring(size(SimBox), CtrlParam%Boost_RepStr, CtrlParam%Boost_RepDist, m_XP)
           return
   end subroutine New_BoostForce
  !****************************************************************************************

  end module MD_BoostMeth_GPU
