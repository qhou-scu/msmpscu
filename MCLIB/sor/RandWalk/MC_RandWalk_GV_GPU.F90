
 module MC_RandWalk_GV
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The module contains globle varaibles for GPU 
 !
 !**** HISTORY:
 !                  version 1st 2020-07 (Hou Qing, Sichuan university)
 !
 !**********************************************************************************   
   use MC_TypeDef_RandWalkCtrl
   use MC_TypeDef_RandWalker_Base
   use MC_TypeDef_RandWalkerState
   implicit none

     real(KINDDF),constant :: dcm_BoxL(3) 
     real(KINDDF),constant :: dcm_BoxU(3) 
     integer,     constant :: dcm_IFPB(3)
     !---- descriptor of the interface to the routine handling rand walking
     private::SteptoEngine
     abstract interface
              subroutine SteptoEngine(RecTime, theWID, theWalkStat)
              !   INPUT:     RecTime,    recording time 
              !              theWID,      type id of particles that to be invoke in the kernel                 
              !              theWalkStat, time point at which the last jump occour
              !
              !   OUTPUT     PreTime,    updated in theWalkStat
              !              NextTime,   updated in theWalkStat 
              !              NJump,      updated in theWalkStat 
              !              Stat,       updated in theWalkStat   
              !              Pos,        updated in theWalkStat 
              use MC_TypeDef_RandWalkerState
              implicit none
              !----   DUMMY Variables
              real(KINDDF)            ::RecTime
              integer                 ::theWID
              type(RandWalkerStat_DEV)::theWalkStat
              end subroutine SteptoEngine  
     end interface

     !---- basic descriptor of  RandWalker
     type WalkingEngine
          procedure(SteptoEngine),  pointer, nopass::pStepto =>null()
     end type WalkingEngine

     type WalkerType_Dev
          integer::IDEV=-1                                              !
          integer::NW  = 0                                              !---  the number of the walker      
          type(WalkerBase),   device, allocatable::Walker(:)            !---  the description of walkers on device
          type(WalkingEngine),        allocatable::WKEngine(:)
     end type WalkerTYpe_Dev
     !--- the globle 
     type(WalkerType_Dev)::dgm_WalkerType(m_MXDEVICE) 

 contains
 !****************************************************************************
 !****************************************************************************
  subroutine Initializt_RandWalk_GV(CtrlParam, Walker)
  !***  PURPOSE:   to check the consistent of inputs
  !
  !     INPUT:     CtrlParam, the controlling parameter
  !                Walker     the walker description 
  !
  !     OUTPUT     
  ! 
    implicit none

      !--- dummy varioables
      type(RandWalkCtrl)  ::CtrlParam
      type(WalkerBase)    ::Walker(:)
      !--- local variables
      integer::IFPD(3), NW, CURDEV, ERR, I, J 
      !type(BoxSetup)::theBox

           call Release_RandWalk_GV() 
           !--- prepare the periodic condition
           select case(CtrlParam%Boxtype)
               case (CP_MC_BOX_INF)
                    IFPD = 0
               case (CP_MC_BOX_1D_SURF, CP_MC_BOX_1D_DSURF)
                    IFPD(1) = 0;  IFPD(2) = 0; IFPD(3) = 0
               case (CP_MC_BOX_2D_SURF, CP_MC_BOX_2D_DSURF)
                    IFPD(1) = 1;  IFPD(2) = 1; IFPD(3) = 0
               case (CP_MC_BOX_3D)
                    IFPD(1) = 1;  IFPD(2) = 1; IFPD(3) = 1
           end select   
           
           !---
           NW = size(Walker)
           ERR = cudaGetDevice(CURDEV)
           do I=1, m_NDEVICE
              ERR = cudaSetDevice(m_DEVICES(I)) 
              dcm_BoxL(1:3) =  CtrlParam%Boxlow(1:3)
              dcm_BoxU(1:3) =  CtrlParam%Boxup(1:3)
              dcm_IFPB(1:3) =  IFPD(1:3)
              
              !--- allocate the walkertype
              dgm_WalkerType(I)%IDEV = m_DEVICES(I)
              dgm_WalkerType(I)%NW = NW
              allocate(dgm_WalkerType(I)%Walker(NW), dgm_WalkerType(I)%WKEngine(NW))
              do J=1, NW
                 dgm_WalkerType(I)%Walker(J)           = Walker(J)
                 dgm_WalkerType(I)%WKEngine(J)%pStepto =>null()
              end do   
           end do
           ERR = cudaSetDevice(CURDEV) 

      return
  end subroutine Initializt_RandWalk_GV
  !****************************************************************************
    
  !****************************************************************************
  subroutine Release_RandWalk_GV()
  !***  PURPOSE:   to check the consistent of inputs
  !
  !     INPUT:     CtrlParam, the controlling parameter
  !                Walker     the walker description 
  !
  !     OUTPUT     
  ! 
    implicit none
 
      !--- dummy varioables
      !--- local variables
      integer::CURDEV, ERR, I
 
           ERR = cudaGetDevice(CURDEV)
           do I=1, size(dgm_WalkerType)
               if(dgm_WalkerType(I)%IDEV .ge. 0) then
                  ERR = cudaSetDevice(dgm_WalkerType(I)%IDEV) 
                  if(allocated(dgm_WalkerType(I)%Walker) ) then
                    deallocate(dgm_WalkerType(I)%Walker)
                  end if  
                  if(allocated(dgm_WalkerType(I)%Walker) ) then
                    deallocate(dgm_WalkerType(I)%WKEngine)
                  end if  
             end if   
           end do
           ERR = cudaSetDevice(CURDEV)
           return
   end subroutine Release_RandWalk_GV
  !****************************************************************************

  !****************************************************************************
  subroutine SetWalkEngine_RandWalk_GV(WID, theEngine)
  !***  PURPOSE:   to associate the walking engine with a walker
  !
  !     INPUT:     WID,     the index of walker
  !                Engine,  the engine
  !
  !     OUTPUT     
  ! 
  implicit none
 
    !--- dummy varioables
      integer::WID
      external theEngine
      interface
         subroutine theEngine(RecTime, theWID, theWalkStat)
           use MC_TypeDef_RandWalkerState
           implicit none
           !----   DUMMY Variables
            real(KINDDF)            ::RecTime
            integer                 ::theWalker
            type(RandWalkerStat_DEV)::theWalkStat
         end subroutine theEngine  
      end interface
     
     !--- local variables
      integer::I

      do I=1, m_NDEVICE
         dgm_WalkerType(I)%WKEngine(WID)%pStepto =>theEngine
      end do   
      return  
  end subroutine SetWalkEngine_RandWalk_GV
  !***************************************************************************

  !****************************************************************************
  subroutine ConductWalk_RandWalk_GV(RecTime, theWID, theWalkStat)
  !***  PURPOSE:   to conduct the walking with the registered procedure
  !
  !   INPUT:     RecTime,     recording time 
  !              theWID,      type id of particles that to be invoke in the kernel                 
  !              theWalkStat, time point at which the last jump occour
  !
  !   OUTPUT     PreTime,    updated in theWalkStat
  !              NextTime,   updated in theWalkStat 
  !              NJump,      updated in theWalkStat 
  !              Stat,       updated in theWalkStat   
  !              Pos,        updated in theWalkStat 
     implicit none
      !----   DUMMY Variables
      real(KINDDF)            ::RecTime
      integer                 ::theWID
      type(RandWalkerStat_DEV)::theWalkStat
     
     !--- local variables
      integer::IDEV
      procedure(SteptoEngine),  pointer::pWalkEngine

         IDEV        =  theWalkStat%DEVN 
         pWalkEngine => dgm_WalkerType(IDEV)%WKEngine(theWID)%pStepto
         if(.not.associated(pWalkEngine)) return

         call pWalkEngine(RecTime, theWID, theWalkStat)

      return  
  end subroutine ConductWalk_RandWalk_GV
!********************************************************************************  


 end module  MC_RandWalk_GV 
 !********************************************************************************

