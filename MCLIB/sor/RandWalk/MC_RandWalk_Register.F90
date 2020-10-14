
 module MC_RandWalk_Register
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The module is to associate walking engine with type of walkers  
 !
 !**** HISTORY:
 !                  version 1st 2020-07 (Hou Qing, Sichuan university)
 !
 !**********************************************************************************   
  use MC_RandWalk_GV
  implicit none


 contains
 !****************************************************************************
  subroutine Register_WalkEngine(CtrlParam, Walker)
  !***  PURPOSE:   to check the consistent of inputs
  !
  !     INPUT:     CtrlParam, the controlling parameter
  !                Walker     the walker description 
  !
  !     OUTPUT     
  ! 
    use MC_RandWalk_Exp_GPU, only:Stepto_RecTime_Exp_1D=>Stepto_RecTime_1D
    use MC_RandWalk_Pow_GPU, only:Stepto_RecTime_Pow_1D=>Stepto_RecTime_1D
    implicit none

      !--- dummy varioables
      type(RandWalkCtrl)  ::CtrlParam
      type(WalkerBase)    ::Walker(:)
      !--- local variables
      integer::WID, I, J 
      !----  
         if(dgm_WalkerType(1)%IDEV < 0 ) then
            call Initializt_RandWalk_GV(CtrlParam, Walker)
         end if
         
         do WID=1, size(Walker)
            select case(Walker(WID)%WalkType) 
                   !--- 
                   case(CP_WALKTYPE_DIM_1D)
                        !--- WTD
                        select case(Walker(WID)%JumpType(1)%WTD_TYPE) 
                               case(CP_WALKTYPE_WTD_EXP) 
                                    call SetWalkEngine_RandWalk_GV(WID, Stepto_RecTime_Exp_1D)
                             
                               case(CP_WALKTYPE_WTD_POW) 
                                    call SetWalkEngine_RandWalk_GV(WID, Stepto_RecTime_Pow_1D)
                        end select              
                   !--- 
                   case(CP_WALKTYPE_DIM_2D)
                   !--- 
                   case(CP_WALKTYPE_DIM_3D)

            end select  
         end do   
      return
  end subroutine Register_WalkEngine
  !****************************************************************************

 end module  MC_RandWalk_Register 
 !********************************************************************************

