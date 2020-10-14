
 module MC_TypeDef_RandWalkerState
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The module is used to define the data type describing the state of a walker
 !
 !**** HISTORY:
 !                  version 1st 2020-05 (Hou Qing, Sichuan university)
 !
 !**********************************************************************************   
   use MSM_MultiGPU_Basic
   use MC_Randwalk_Const
   implicit none

     type RandWalkerStat
          integer                 ::mxNP          = 0                       !max number of partilces allowed in the simulation
          integer                 ::NBox          = 1                       !number of parallele boxes
          integer                 ::MxNPPerBox    = 0                       !max number of parallele in each box mxNP/NBOX
          integer                 ::NWalkStyID      = 0                     !number of walk styles
          integer,     allocatable::FirstPIDP(:)                            !the first particle id in a box                                                     
          integer,     allocatable::CurNP(:)                                !current numnber of particles in the boxes
          integer     ,allocatable::WalkStyID(:)                            !ID of the walker style for the walker follwing
          integer     ,allocatable::Stat(:)                                 !state of the walkers
          real(KINDDF),allocatable::PreTime(:)                              !time of the last jump occurring
          real(KINDDF),allocatable::NextTime(:)                             !time of the next jump to occur
          integer,     allocatable::NextJPath(:)                            !jump path on the next time
          integer,     allocatable::NextJump(:)                             !jump vector id on the next time, refer to the definition of WalkerBase 
          real(KINDDF),allocatable::NJump(:)                                !number of jumps
          real(KINDDF),allocatable::XP(:,:)                                 !current position of particle
     end type RandWalkerStat

     type RandWalkerStat_DEV
          integer           ::DEVN          = -1
          integer           ::mxNP          =  0                            !max number of partilces allowed in the simulation
          integer           ::NBox          =  1                            !number of parallele boxes
          integer           ::mxNPPerBox    =  0                            !max number of parallele in each box mxNP/NBOX
          type(DevVec_I)    ::dFirstPIDP          
          type(DevVec_I)    ::dCurNP              
          type(DevVec_I)    ::dWalkStyID
          type(DevVec_I)    ::dStat
          type(DevVec_DF)   ::dPreTime
          type(DevVec_DF)   ::dNextTime
          type(DevVec_I)    ::dNextJPath
          type(DevVec_I)    ::dNextJump
          type(DevVec_DF)   ::dNJump
          type(DevMat_DF)   ::dXP
     end type RandWalkerStat_DEV

     !------------------------------------------------------
     private::  Archive_Config_0,         &
                Restore_Config_0
     interface  Archive_RandWalkerStat
       module procedure Archive_Config_0 
     end interface  Archive_RandWalkerStat

     interface  Restore_RandWalkerStat
       module procedure Restore_Config_0 
     end interface  Restore_RandWalkerStat

     !------------------------------------------------------
     private::  Copyin_RandWalkerStat0,   &
                Copyin_RandWalkerStat1
     interface  Copyin_RandWalkerStat
       module procedure Copyin_RandWalkerStat0 
          module procedure Copyin_RandWalkerStat1
     end interface  Copyin_RandWalkerStat

     !------------------------------------------------------
     private::  CopyOut_RandWalkerStat0
     interface  CopyOut_RandWalkerStat
       module procedure CopyOut_RandWalkerStat0 
     end interface  CopyOut_RandWalkerStat

     !------------------------------------------------------
     public::   DevAllocate_RandWalkerStat
     
     !------------------------------------------------------
     private::  Initialize_RandWalkerStat0
     interface  Initialize_RandWalkerStat
       module procedure Initialize_RandWalkerStat0 
     end interface  Initialize_RandWalkerStat

     !------------------------------------------------------
     private::  InBoxParticles0,       &
                InBoxParticles1
     interface  InBoxParticles_RandWalkerStat
          module procedure InBoxParticles0
            module procedure InBoxParticles1
     end interface  InBoxParticles_RandWalkerStat

     !------------------------------------------------------
     !   routines concerning insert atoms
     private::  Insert_Uniform_KERNEL,    &
                InsertAtoms_Uniform0,     &
                InsertAtoms_Uniform1,     &
                InsertAtoms_Uniform2
     interface  InsertAtoms_Uniform
        module procedure InsertAtoms_Uniform0     
           module procedure InsertAtoms_Uniform1
            module procedure InsertAtoms_Uniform2
     end interface  InsertAtoms_Uniform 

     private::  Insert_ByHist_KERNEL0,    &
                InsertAtoms_ByHist0,      &
                InsertAtoms_ByHist0_a,    &
                InsertAtoms_ByHist0_b,    &
                InsertAtoms_ByHist0_c
     interface  InsertAtoms_ByHis
        module procedure InsertAtoms_ByHist0     
           module procedure InsertAtoms_ByHist0_a
            module procedure InsertAtoms_ByHist0_b
              module procedure InsertAtoms_ByHist0_c
     end interface  InsertAtoms_ByHis 

     !------------------------------------------------------
     private::  Release_RandWalkerStat0, &
                Release_RandWalkerStat_DEV 
     interface  Release_RandWalkerStat
       module procedure Release_RandWalkerStat0 
          module procedure Release_RandWalkerStat_DEV
     end interface  Release_RandWalkerStat

     !------------------------------------------------------
     private::  Sweep0,                 &
                Sweep1,                 &
                Sweep2,                 &
                Sweep2D,                &
                Sweep3D,                &
                Sweep3,                 &
                Sweep4
                
     interface  Sweep_RandWalkerStat
       module procedure Sweep0 
         module procedure Sweep1
            module procedure Sweep2
               module procedure Sweep3
                  module procedure Sweep4
     end interface  Sweep_RandWalkerStat
     !------------------------------------------------------

     !---- INTERFACE LIST:
     public::   Copyin_RandWalkerStat
     public::   CopyOut_RandWalkerStat
     public::   Initialize_RandWalkerStat
     public::   InsertAtoms_Uniform
     public::   Putout_CurConfig_RandWalkerStat
     public::   Release_RandWalkerStat
     public::   Sweep_RandWalkerStat

 contains
 
 !****************************************************************************
  subroutine Release_RandWalkerStat0(WalkerState)
  !***  PURPOSE:   to release the allocated memory of RandWalkerStat
  !
  !     INPUT:     WalkerState
  !
  !     OUTPUT     WalkerState, the object with memery deallocated
   implicit none
     !--- dummy varioables
     class(RandWalkerStat) ::WalkerState
     !--- local variables
      
     !----
           WalkerState%mxNP        = 0
           WalkerState%NBOX        = 0
           WalkerState%mxNPPerBox  = 0
           if(allocated(WalkerState%FirstPIDP)) deallocate(WalkerState%FirstPIDP) 
           if(allocated(WalkerState%CurNP))     deallocate(WalkerState%CurNP) 
           if(allocated(WalkerState%WalkStyID)) deallocate(WalkerState%WalkStyID)
           if(allocated(WalkerState%Stat))      deallocate(WalkerState%Stat)
           
           if(allocated(WalkerState%PreTime))   deallocate(WalkerState%PreTime)
           if(allocated(WalkerState%NextTime))  deallocate(WalkerState%NextTime)
           if(allocated(WalkerState%NextJPath)) deallocate(WalkerState%NextJPath)
           if(allocated(WalkerState%NextJump))  deallocate(WalkerState%NextJump)
           if(allocated(WalkerState%NJump))     deallocate(WalkerState%NJump)
           if(allocated(WalkerState%XP))        deallocate(WalkerState%XP)
 
           return
  end subroutine Release_RandWalkerStat0
  !****************************************************************************

 !****************************************************************************
  subroutine Release_RandWalkerStat_DEV(WalkerState)
  !***  PURPOSE:   to release the memory of RandWalkerStat
  !
  !     INPUT:     WalkerState
  !
  !     OUTPUT     WalkerState, the object with memery deallocated
   implicit none
     !--- dummy varioables
     class(RandWalkerStat_DEV) ::WalkerState
     !--- local variables
      
            

           call DevDeAllocate(WalkerState%dFirstPIDP)
           call DevDeAllocate(WalkerState%dcurNP)
           call DevDeAllocate(WalkerState%dWalkStyID)
           call DevDeAllocate(WalkerState%dStat)
           
           call DevDeAllocate(WalkerState%dPreTime)
           call DevDeAllocate(WalkerState%dNextTime)
           call DevDeAllocate(WalkerState%dNextJPath)
           call DevDeAllocate(WalkerState%dNextJump)
           call DevDeAllocate(WalkerState%dNJump)
           call DevDeAllocate(WalkerState%dXP)

           WalkerState%mxNP       = 0
           WalkerState%NBox       = 0
           WalkerState%mxNPPerBox = 0
           WalkerState%DEVN       = -1


           return
  end subroutine Release_RandWalkerStat_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_RandWalkerStat0(WalkerState, CtrlParam, Walker)
  !***  PURPOSE:   to initialize state of the walke
  !
  !     INPUT:     NBOX,      number of box in the workspace
  !                MXNPABOX,  max number of particles in each box   
  !                MXWID,     max number of particle types    
  !
  !     OUTPUT     WalkerState, with memoery allocated
    use MC_TypeDef_RandWalkCtrl
    use MC_TypeDef_RandWalker_Base
    implicit none
     !--- dummy varioables
     class(RandWalkerStat) :: WalkerState
     type(RandWalkCtrl)    :: CtrlParam
     type(WalkerBase)      :: Walker(:)

     !--- local variables
      integer::I, NBOX, MXNPABOX, MXNP, ERR
     !----
           call Release_RandWalkerStat(WalkerState)
           call GetMXWalkerNum_RandWalkCtrl(CtrlParam, MXNPABOX, NBOX) 

           allocate(WalkerState%FirstPIDP(NBOX),  WalkerState%CurNP(NBOX))

           MXNP = NBOX*MXNPABOX
           allocate(WalkerState%WalkStyID(MXNP),   &
                    WalkerState%Stat(MXNP),        &
                    WalkerState%PreTime(MXNP),     &
                    WalkerState%NextTime(MXNP),    &
                    WalkerState%NextJPath(MXNP),   &
                    WalkerState%NextJump(MXNP),    &
                    WalkerState%NJump(MXNP),       &
                    WalkerState%XP(MXNP,3),        &
                    STAT=ERR)
           if(ERR .gt. 0) then
              write(*,*) "MCLIB Error:: fail allocating memory for host worksapce of WalkerStat"
              write(*,*) "              MXNP =", MXNP, "MXNPABOX=", MXNPABOX
              stop   
           end if 

           WalkerState%NBOX        = NBOX
           WalkerState%mxNPPerBox  = MXNPABOX
           WalkerState%mxNP        = MXNP
           WalkerState%NWalkStyID  = size(Walker)
           do I=1, NBOX
              WalkerState%FirstPIDP(I) = (I-1)*MXNPABOX + 1
           end do 
           WalkerState%CurNP       = 0
           WalkerState%WalkStyID   = CP_WALKTYPE_NOWALK
           WalkerState%Stat        = CP_WALKSTAT_NOTACTIVE
           WalkerState%PreTime     = 0.D0
           WalkerState%NextTime    = 0.D0
           WalkerState%NextJPath   = 0
           WalkerState%NextJump    = 0
           WalkerState%NJump       = 0.D0
           WalkerState%XP          = 0.D0


           return
  end subroutine Initialize_RandWalkerStat0
  !****************************************************************************

  !****************************************************************************
  subroutine DevAllocate_RandWalkerStat(DEVN, NBOX, MXNP, StatGPU)
  !***  PURPOSE:   to allocate memory for a RandWalkerStat_DEV structure
  !
  !     INPUT:     IDEV,    device ID
  !
  !     OUTPUT     StatGPU
   implicit none
     !--- dummy varioables
     integer                    ::DEVN, NBOX, MXNP 
     class(RandWalkerStat_DEV)  ::StatGPU
     !--- local variables
      
     !----
               call Release_RandWalkerStat_DEV(StatGPU)
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dFirstPIDP,    NBOX,     "dFirstPIDP"  )
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dCurNP,        NBOX,     "dCurNP"      ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dWalkStyID,    mxNP,     "dWalkStyID"  ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dStat,         mxNP,     "dStat"       )

               call DevAllocate(m_DEVICES(DEVN), StatGPU%dPreTime,      mxNP,     "dPreTime"    ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dNextTime,     mxNP,     "dNextTime"   ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dNextJPath,    mxNP,     "dNextJPath"  ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dNextJump,     mxNP,     "dNextJump"   ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dNJump,        mxNP,     "dNJump"      ) 
               call DevAllocate(m_DEVICES(DEVN), StatGPU%dXP,         (/mxNP,3/), "dXP"         ) 
               StatGPU%DEVN = DEVN
           return
  end subroutine DevAllocate_RandWalkerStat
  !****************************************************************************
  !****************************************************************************
  subroutine Copyin_RandWalkerStat0(DEVN, StatCPU, StatGPU)
  !***  PURPOSE:   to copyin the stat
  !
  !     INPUT:     IDEV,    device ID
  !                StatCPU,  
  !
  !     OUTPUT     StatGPU
   implicit none
     !--- dummy varioables
     integer                    ::DEVN 
     class(RandWalkerStat)      ::StatCPU
     class(RandWalkerStat_DEV)  ::StatGPU
     !--- local variables
      
     !----
           if(StatGPU%DEVN .ne. DEVN ) then
               call DevAllocate_RandWalkerStat(DEVN, StatCPU%NBOX, StatCPU%mxNP, StatGPU)
           end if    
           call DevCopyIn(StatCPU%FirstPIDP,  StatGPU%dFirstPIDP, StatCPU%NBOX)
           call DevCopyIn(StatCPU%CurNP,      StatGPU%dCurNP,     StatCPU%NBOX)
           call DevCopyIn(StatCPU%WalkStyID,  StatGPU%dWalkStyID, StatCPU%mxNP)
           call DevCopyIn(StatCPU%Stat,       StatGPU%dStat,      StatCPU%mxNP)
           call DevCopyIn(StatCPU%PreTime,    StatGPU%dPreTime,   StatCPU%mxNP)
           call DevCopyIn(StatCPU%NextTime,   StatGPU%dNextTime,  StatCPU%mxNP)
           call DevCopyIn(StatCPU%NextJPath,  StatGPU%dNextJPath, StatCPU%mxNP)
           call DevCopyIn(StatCPU%NextJump,   StatGPU%dNextJump,  StatCPU%mxNP)
           call DevCopyIn(StatCPU%Njump,      StatGPU%dNJump,     StatCPU%mxNP)
           call DevCopyIn(StatCPU%XP,         StatGPU%dXP,        StatCPU%mxNP*3)

           StatGPU%mxNP = StatCPU%mxNP 
           StatGPU%NBox = StatCPU%NBox 
           StatGPU%mxNPPerBox = StatCPU%mxNPPerBox 
 
           return
  end subroutine Copyin_RandWalkerStat0
  !****************************************************************************
  
  !****************************************************************************
  subroutine Copyin_RandWalkerStat1(StatCPU, StatGPU)
  !***  PURPOSE:   to copyin the stat
  !
  !     INPUT:     IDEV,    device ID
  !                StatCPU,  
  !
  !     OUTPUT     StatGPU
    implicit none
     !--- dummy varioables
     class(RandWalkerStat)     ::StatCPU
     class(RandWalkerStat_DEV) ::StatGPU
     !--- local variables
      
     !----
           if(StatGPU%DEVN .lt. 0) then
               call Copyin_RandWalkerStat0(1, StatCPU, StatGPU)
           else
               call Copyin_RandWalkerStat0(StatGPU%DEVN, StatCPU, StatGPU)    
           end if
 
           return
  end subroutine Copyin_RandWalkerStat1
  !****************************************************************************  

  !****************************************************************************
  subroutine CopyOut_RandWalkerStat0(StatCPU, StatGPU)
  !***  PURPOSE:   to initialize state of the walke
  !
  !     INPUT:     STR,    a string line
  !
  !     OUTPUT     Walker
   implicit none
     !--- dummy varioables
     class(RandWalkerStat)     ::StatCPU
     class(RandWalkerStat_DEV) ::StatGPU
     !--- local variables
      
     !----
           call DevCopyOut(StatCPU%FirstPIDP,  StatGPU%dFirstPIDP, StatCPU%NBOX)
           call DevCopyOut(StatCPU%curNP,      StatGPU%dcurNP,     StatCPU%NBOX)
           call DevCopyOut(StatCPU%WalkStyID,  StatGPU%dWalkStyID, StatCPU%mxNP)
           call DevCopyOut(StatCPU%Stat,       StatGPU%dStat,      StatCPU%mxNP)
           call DevCopyOut(StatCPU%PreTime,    StatGPU%dPreTime,   StatCPU%mxNP)
           call DevCopyOut(StatCPU%NextTime,   StatGPU%dNextTime,  StatCPU%mxNP)
           call DevCopyOut(StatCPU%NextJPath,  StatGPU%dNextJPath, StatCPU%mxNP)
           call DevCopyOut(StatCPU%NextJump,   StatGPU%dNextJump,  StatCPU%mxNP)
           call DevCopyOut(StatCPU%Njump,      StatGPU%dNJump,     StatCPU%mxNP)
           call DevCopyOut(StatCPU%XP,         StatGPU%dXP,        StatCPU%mxNP*3)

           return
  end subroutine CopyOut_RandWalkerStat0
  !****************************************************************************  
  !****************************************************************************

  !*********************************************************************************
  !**********************************************************************************
  attributes(global) subroutine Insert_Uniform_KERNEL(MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, ID, X0, X1, Y0, Y1, Z0, Z1, &
                                                      CurNP, WalkID, Stat, Pos, Pretime, NextTime, NextJPath, NextJump, NJump, RandState)
  !***  PURPOSE:   to insert partitcles in the box uniform in depth from Z0-Z1
  !                NOTE: the current number of particle,  CurNP is not changed in this kernel
  !
  !$$   INPUT:     MXNP,       the size of Pos  
  !$$              NBOX,       the number of box
  !$$              MXNPABOX,   the permitted number of particles in a box  
  !$$              When,       the time of inserting the particles
  !$$              Deltat.     the time interval of activating the partilces                                                     
  !$$              NewNP,      the number of particles to be added to each box 
  !$$              WID,        the type ID of the added particles   
  !$$              X0, X1, Y0, Y1, Z0, Z1,     the volume to insert the particles
  !$$              CurNP,      the current number of particles in boxes
  !$$              
  !$$              RandState,  random number state
  !$$
  !$$   OUTPUT:    WalkID,     the type ID of the added particles
  !$$              Pos,        position of the atoms
  !$$              Stat,       stat of the atoms, initially being active
  !$$              Pretime,    time at which the particles are inserted
  !$$              Nexttime,   time at which ot triger the next jump, here Nexttime = Pretime
  !$$              NextJPath,   next jump path ID,  initially to be 0 
  !$$              NextJump,   next jump vector ID, initially to be 0
  !$$              NJump,      number of jumps
  !$$              RandState,  updated random number state
  !
  implicit none
  !----   DUMMY Variables
          integer,value                  :: MXNP, NBOX, MXNPABOX, NewNP, ID
          real(KINDDF),value             :: X0, X1, Y0, Y1, Z0, Z1, When, Deltat 
          integer,     device            :: CurNP(*)   
          integer,     device            :: WalkID(*)
          integer,     device            :: Stat(*)
          real(KINDDF),device            :: Pos(MXNP,*)
          real(KINDDF),device            :: Pretime(*)
          real(KINDDF),device            :: Nexttime(*)
          integer,     device            :: NextJPath(*)
          integer,     device            :: NextJump(*)
          real(KINDDF),device            :: NJump(*)
          type(curandStateXORWOW), device:: RandState(*)
          

  !----   Local variables
          integer::IT, IB, IC, BS, GS, IC0, NL, L, INEW, IBOX
          real(KINDDF)::Time
  !
  !***
  !     
              !--- size the block and grid
              GS  = griddim%x*griddim%y
              BS  = blockdim%x*blockdim%y
              !--- number of loops needed for cover all atoms
              NL  = (NBOX*NewNP - 1)/(GS*BS) + 1

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x

              !IC -- the id of the random number sequence 
              IC0  = (IB-1)*BS + IT

              do L = 1, NL
                 !INEW -- the id of the new particle
                 INEW = (L-1)*(BS*GS) + IC0 

                 !IBOX -- the box id that the new particle to be added to
                 IBOX  = (INEW-1)/NewNP  + 1

                 !IC -- index of added particles in the box
                 if(IBOX .le. NBOX) then
                    IC = (INEW - (IBOX-1)*NewNP) + CurNP(IBOX) + (IBOX-1)*MXNPABOX
                    if(IC .le. MXNP) then
                      WalkID(IC) = ID  
                      Stat(IC)   = CP_WALKSTAT_WAITING
                      !position
                       Pos(IC,1) = X0 + (X1 - X0)*curand_uniform(RandState(IC0))
                       Pos(IC,2) = Y0 + (Y1 - Y0)*curand_uniform(RandState(IC0))
                       Pos(IC,3) = Z0 + (Z1 - Z0)*curand_uniform(RandState(IC0))
                       !Time and number of Jump
                       TIME         = When + Deltat*curand_uniform(RandState(IC0))
                       Pretime(IC)  = TIME 
                       Nexttime(IC) = TIME 
                       NJump(IC)    = 0.D0
                       NextJPath(IC)= 0
                       NextJump(IC) = 0
                    end if   
                 end if      
               end do

         return
  end subroutine Insert_Uniform_KERNEL
  !********************************************************************************   

  !********************************************************************************
  subroutine InsertAtoms_Uniform0(MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, WID, X0, X1, Y0, Y1, Z0, Z1, CurNP, WalkID, Stat, Pos, &
                                  PreTime, NextTime, NextJPath, NextJump, NJump, RandState)
  !***  PURPOSE:   to insert partitcles in the box uniform in depth from Z0-Z1
  !
  !     INPUT:     MXNP,       the size of Pos  
  !                NBOX,       the number of box
  !                MXNPABOX,   the permitted number of particles in a box  
  !                NewNP,      the number of particles to be added to each box 
  !                WID,        the type ID of the added particles   
  !                X), X1, Y0, Y1, Z0, Z1,     the volume to insert the particles
  !                CurNP,      the current number of particles in boxes
  !                When,       the time of inserting the particles
  !                Deltat,     the time interval of activating the partilces      
  !              
  !                RandState,  random number state
  !
  !     OUTPUT     WalkID,     type ID of the added particles
  !                Stat,       stat of the atoms, initially being active  
  !                Pos,        position of the atoms
  !                PreTime,    time at which the particles are inserted and activated
  !                NextTime,   time at which the particles are inserted and activated
  !                NextJPath,   next jump path ID,  initially to be 0 
  !                NextJump,   next jump vector ID, initially to be 0
  !                NJump,      initial number of jumps (0)   
  !                RandState,  updated random number state
  !
   implicit none
   !----   DUMMY Variables
     integer             :: MXNP, NBOX, MXNPABOX, NEWNP, WID
     real(KINDDF)        :: When, Deltat, X0, X1, Y0, Y1, Z0, Z1    
     type(DevVec_I)      :: CurNP
     type(DevVec_I)      :: WalkID
     type(DevMat_DF)     :: Pos 
     type(DevVec_I)      :: Stat
     type(DevVec_DF)     :: PreTime 
     type(DevVec_DF)     :: NextTime 
     type(DevVec_I)      :: NextJPath
     type(DevVec_I)      :: NextJump
     type(DevVec_DF)     :: NJump 
     type(DevRandState)  :: RandState
   !-- local varibales  
     integer::ERR, CURDEV 
   !--- Device variables and variables to be used in GPU
     type(dim3) :: blocks
     type(dim3) :: threads
           
   !---
       ERR = cudaGetDevice(CURDEV)
       ERR = cudaSetDevice(RandState%IDEV)
       blocks  = dim3(RandState%GRIDSIZE, 1, 1)
       threads = dim3(RandState%BLOCKSIZE,  1, 1)
       
       call Insert_Uniform_KERNEL<<<blocks, threads>>>(MXNP, NBOX, MXNPABOX, When, Deltat, NEWNP, WID, X0, X1, Y0, Y1, Z0, Z1,   &
                                                      CurNP%Data, WalkID%Data, Stat%Data, Pos%Data, PreTime%Data,NextTime%Data,  &
                                                      NextJPath%Data, NextJump%Data, NJump%Data, RandState%RandState)
       !--- because curNP is in not updated kernal, we should update it 
       !    on host and then copyin to device
       !  we have replace DevAdd by updat CURN on host, and then copyin CURN from Host
       !  See: InsertNewPart_WalkEvolve                                               
       !call DevAdd(CurNP, NEWNP)
       ERR = cudaSetDevice(CURDEV)
       return 
  end subroutine InsertAtoms_Uniform0
  !******************************************************************************** 

  !********************************************************************************
  subroutine InsertAtoms_Uniform1(DEVN, MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, WID, Vol, CurNP, WalkID, Stat, Pos, &
                                  PreTime, NextTime, NextJPath, NextJump, NJump)
  !***  PURPOSE:   to insert partitcles in the box uniform in depth from Z0-Z1
  !
  !   INPUT:     MXNP,       the size of Pos  
  !              NBOX,       the number of box
  !              MXNPABOX,   the permitted number of particles in a box  
  !              When,       the time of inserting the particles
  !              Deltat,     the time interval of activating the partilces      
  !              NewNP,      the number of particles to be added to each box 
  !              WID,        the type ID of the added particles
  !              Vol,        the volume to insert the particles
  !              CurNP,      the current number of particles in boxes
  !
  !   OUTPUT     WalkID,     walker id of the atoms
  !              Pos,        position of the atoms,   
  !              Stat,       stat of the atoms, initially being active  
  !              PreTime,    time at which the particles are inserted and activated
  !              NextTime,   time at which the particles are inserted and activated
  !              NJump,      initial number of jumps (0)   
  !
   implicit none
   !----   DUMMY Variables
     integer             :: DEVN, MXNP, NBOX, MXNPABOX, NewNP, WID
     real(KINDDF)        :: When, Deltat
     real(KINDDF)        :: Vol(6)
     type(DevVec_I)      :: CurNP
     type(DevVec_I)      :: WalkID
     type(DevVec_I)      :: Stat
     type(DevMat_DF)     :: Pos 
     type(DevVec_DF)     :: PreTime 
     type(DevVec_DF)     :: NextTime
     type(DevVec_I)      :: NextJPath
     type(DevVec_I)      :: NextJump
     type(DevVec_DF)     :: NJump 
   !-- local varibales  
       
       call InsertAtoms_Uniform0(MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, WID, Vol(1), Vol(2), Vol(3), Vol(4), Vol(5), Vol(6), &
                                 CurNP, WalkID, Stat, Pos, PreTime, NextTime, NextJPath, NextJump, NJump, gm_DevRandState(DevN))
       return 
  end subroutine InsertAtoms_Uniform1
  !********************************************************************************

  !********************************************************************************
  subroutine InsertAtoms_Uniform2(When, Deltat, NewNP, WID, Vol, WalkStat)
  !***  PURPOSE:   to insert partitcles in the box uniform in volume Vol
  !
  !   INPUT:     When,       the time at which the particles inserted 
  !              Deltat,     the time interval of activating the partilces      
  !              NewNP,      the number of pariticles to be inserted
  !              WID,        the walker ID to be inserted   
  !              Vol,        the volume where the particle to be inserted
  !
  !   OUTPUT     WalkStat,     
  !
  !
   implicit none
   !----   DUMMY Variables
     integer                   :: MXNP, CURNP,NEWNP, WID
     real(KINDDF)              :: When, Deltat, Vol(:)
     class(RandWalkerStat_DEV) :: WalkStat
   !-- local varibales  
       
       call InsertAtoms_Uniform1(WalkStat%DEVN, WalkStat%mxNP,     WalkStat%NBOX,       WalkStat%mxNPPerBox, When, Deltat, NewNP, WID, Vol, &
                                                WalkStat%dCurNP,   WalkStat%dWalkStyID, WalkStat%dStat, WalkStat%dXP,                       &
                                                WalkStat%dPreTime, WalkStat%dNextTime,  WalkStat%dNextJPath, WalkStat%dNextJump, WalkStat%dNJump)
       return 
  end subroutine InsertAtoms_Uniform2
  !********************************************************************************

  !**********************************************************************************
  attributes(global) subroutine Insert_ByHist_KERNEL0(MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, ID, X0, X1, Y0, Y1, NDepth, Depth, DDis,  &
                                                     CurNP, WalkID, Stat, Pos, Pretime, NextTime, NextJPath, NextJump, NJump,RandState)
  !***  PURPOSE:   to insert partitcles in the box by depth histogram 
  !                NOTE: For the histogram, thw width of the bins are non-constant, by the particle number in all bins
  !                      are the same: NewNP/(number of bins)
  !
  !$$   INPUT:     MXNP,       the size of Pos  
  !$$              NBOX,       the number of box
  !$$              MXNPABOX,   the permitted number of particles in a box  
  !$$              When,       the time of inserting the particles
  !$$              Deltat.     the time interval of activating the partilces                                                     
  !$$              NewNP,      the number of particles to be added to each box 
  !$$              ID,         the type ID of the added particles   
  !$$              X0, X1, Y0, Y1   the area to insert the particles                                                     
  !$$              NDepth,     the number of grids of depths
  !$$              Depth,      the depth grids
  !$$              DDis,       the integral depth distribution of Walker ID                                                     
  !$$              CurNP,      the current number of particles in boxes
  !$$              
  !$$              RandState,  random number state
  !$$
  !$$   OUTPUT     WalkID,     the type ID of the added particles
  !$$              Stat,       stat of the atoms, initially being active
  !$$              Pos,        position of the atoms
  !$$              Pretime,    time at which the particles are inserted
  !$$              Nexttime,   time at which ot triger the next jump, here Nexttime = Pretime
  !$$              NextJPath,   next jump path ID,  initially to be 0 
  !$$              NextJump,   next jump vector ID, initially to be 0
  !$$              NJump,      number of jumps
  !$$              RandState,  updated random number state
  !
  implicit none
  !----   DUMMY Variables
          integer,value                  :: MXNP, NBOX, MXNPABOX, NewNP, ID, NDepth
          real(KINDDF),value             :: When, Deltat, X0, X1, Y0, Y1 
          real(KINDDF),device            :: Depth(*)
          real(KINDDF),device            :: DDis(*)
          integer,     device            :: CurNP(*)   
          integer,     device            :: WalkID(*)
          integer,     device            :: Stat(*)
          real(KINDDF),device            :: Pos(MXNP,*)
          real(KINDDF),device            :: Pretime(*)
          real(KINDDF),device            :: Nexttime(*)
          integer,     device            :: NextJPath(*)
          integer,     device            :: NextJump(*)
          real(KINDDF),device            :: NJump(*)
          type(curandStateXORWOW), device:: RandState(*)
          

  !----   Local variables
          integer::IT, IB, IC, BS, GS, IC0, NL, L, INEW, IBOX, IR, IR0, IR1
          real(KINDDF)::R, Z0, Z1, ZM, TIME
  !
  !***
  !     
              !--- size the block and grid
              GS  = griddim%x*griddim%y
              BS  = blockdim%x*blockdim%y
              !--- number of loops needed for cover all atoms
              NL  = (NBOX*NewNP - 1)/(GS*BS) + 1

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x

              !IC -- the id of the random number sequence 
              IC0  = (IB-1)*BS + IT

              do L = 1, NL
                 !INEW -- the id of the new particle
                 INEW = (L-1)*(BS*GS) + IC0 

                 !IBOX -- the box id that the new particle to be added to
                 IBOX  = (INEW-1)/NewNP  + 1

                 !IC -- index of added particles in the box
                 if(IBOX .le. NBOX) then
                    IC = (INEW - (IBOX-1)*NewNP) + CurNP(IBOX) + (IBOX-1)*MXNPABOX
                    if(IC .le. MXNP) then
                      WalkID(IC) = ID  
                      Stat(IC)   = CP_WALKSTAT_ACTIVE
                       !position
                       Pos(IC,1) = X0 + (X1 - X0)*curand_uniform(RandState(IC0))
                       Pos(IC,2) = Y0 + (Y1 - Y0)*curand_uniform(RandState(IC0))

                       !--- randomly select a bin
                       R         = curand_uniform(RandState(IC0))
                       IR0       = 1
                       IR1       = NDepth
                       IR        = (IR0 + IR1)/2
                       do while(IR1 - IR0 .gt. 1)  
                           if(R .gt. DDis(IR)) then
                              IR0 = IR 
                           else if(R .lt. DDis(IR)) then
                              IR1 = IR 
                           end if   
                           IR = (IR0 + IR1)/2      
                       end do     
                       !IR        = int(R*NDepth) + 1

                       !--- randomly set postion Z in the bin 
                       Z0        = Depth(IR)
                       Z1        = Depth(IR+1)
                       Pos(IC,3) = Z0 + (Z1 - Z0)*curand_uniform(RandState(IC0))
                       !Time and number of Jump
                       TIME         = When + Deltat*curand_uniform(RandState(IC0))
                       Pretime(IC)  = TIME 
                       Nexttime(IC) = TIME 
                       NJump(IC)    = 0.D0
                       NextJPath(IC)= 0
                       NextJump(IC) = 0
                       
                    end if   
                 end if      
               end do

         return
  end subroutine Insert_ByHist_KERNEL0
  !********************************************************************************   

  !********************************************************************************
  subroutine InsertAtoms_ByHist0(MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, WID, X0, X1, Y0, Y1, NDepth, Depth, DDis, CurNP, WalkID, Stat, Pos, &
                                 PreTime, NextTime, NextJPath, NextJump, NJump, RandState)
  !***  PURPOSE:   to insert partitcles in the box according to given depth distribution
  !
  !     INPUT:     MXNP,       the size of Pos  
  !                NBOX,       the number of box
  !                MXNPABOX,   the permitted number of particles in a box  
  !                When,       the time of inserting the particles
  !                Deltat,     the time interval of activating the partilces      
  !                NewNP,      the number of particles to be added to each box 
  !                WID,        the type ID of the added particles   
  !                X0, X1, Y0, Y1   the area to insert the particles
  !                NDepth,     the number of grids of depths
  !                Depth,      the depth grids
  !                DDis,       the integral depth distribution of Walker ID       
  !                CurNP,      the current number of particles in boxes
  !              
  !                RandState,  random number state
  !
  !     OUTPUT     WalkID,     type ID of the added particles
  !                Stat,       stat of the atoms, initially being active  
  !                Pos,        position of the atoms
  !                PreTime,    time at which the particles are inserted and activated
  !                NextTime,   time at which the particles are inserted and activated
  !                NextJPath,   next jump path ID,  initially to be 0 
  !                NextJump,   next jump vector ID, initially to be 0
  !                NJump,      initial number of jumps (0)   
  !                RandState,  updated random number state
  !
   implicit none
   !----   DUMMY Variables
     integer             :: MXNP, NBOX, MXNPABOX, NEWNP, WID, NDepth
     real(KINDDF)        :: When, Deltat, X0, X1, Y0, Y1
     type(DevVec_DF)     :: Depth
     type(DevVec_DF)     :: DDis
     type(DevVec_I)      :: CurNP
     type(DevVec_I)      :: WalkID
     type(DevMat_DF)     :: Pos 
     type(DevVec_I)      :: Stat
     type(DevVec_DF)     :: PreTime 
     type(DevVec_DF)     :: NextTime 
     type(DevVec_I)      :: NextJPath
     type(DevVec_I)      :: NextJump
     type(DevVec_DF)     :: NJump 
     type(DevRandState)  :: RandState
   !-- local varibales  
     integer::ERR, CURDEV 
   !--- Device variables and variables to be used in GPU
     type(dim3) :: blocks
     type(dim3) :: threads
           
   !---
       ERR = cudaGetDevice(CURDEV)
       ERR = cudaSetDevice(RandState%IDEV)
       blocks  = dim3(RandState%GRIDSIZE, 1, 1)
       threads = dim3(RandState%BLOCKSIZE,  1, 1)
       
       call Insert_ByHist_KERNEL0<<<blocks, threads>>>(MXNP, NBOX, MXNPABOX, When, Deltat, NEWNP, WID, X0, X1, Y0, Y1, NDepth, Depth%Data, DDis%Data, &
                                                      CurNP%Data, WalkID%Data, Stat%Data, Pos%Data, PreTime%Data,NextTime%Data,                       &
                                                      NextJPath%Data, NextJump%Data, NJump%Data, RandState%RandState)
       !--- because curNP is in not updated kernal, we should update it 
       !    on host and then copyin to device
       !  we have replace DevAdd by updat CURN on host, and then copyin CURN from Host
       !  See: InsertNewPart_WalkEvolve                                               
       !call DevAdd(CurNP, NEWNP)
       ERR = cudaSetDevice(CURDEV)
       return 
  end subroutine InsertAtoms_ByHist0
  !******************************************************************************** 

  !********************************************************************************
  subroutine InsertAtoms_ByHist0_a(DEVN, MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, WID, Area, NDepth, Depth, DDis, CurNP, WalkID, Stat, Pos, &
                                   PreTime, NextTime, NextJPath, NextJump, NJump)
  !   PURPOSE:   to insert partitcles in the box according to given depth distribution
  !  
  !   INPUT:     MXNP,       the size of Pos  
  !              NBOX,       the number of box
  !              MXNPABOX,   the permitted number of particles in a box  
  !              When,       the time of inserting the particles
  !              Deltat,     the time interval of activating the partilces      
  !              NewNP,      the number of particles to be added to each box 
  !              WID,        the type ID of the added particles
  !              AREA,       the area to insert the particles
  !              NDepth,     the number of grids of depths
  !              Depth,      the depth grids
  !              DDis,       the integral depth distribution of Walker ID       
  !              CurNP,      the current number of particles in boxes
  !
  !   OUTPUT     WalkID,     walker id of the atoms
  !              Stat,       stat of the atoms, initially being active  
  !              Pos,        position of the atoms
  !              PreTime,    time at which the particles are inserted and activated
  !              NextTime,   time at which the particles are inserted and activated
  !              NextJPath,   next jump path ID,  initially to be 0 
  !              NextJump,   next jump vector ID, initially to be 0
  !              NJump,      initial number of jumps (0)   
  !
   implicit none
   !----   DUMMY Variables
     integer             :: DEVN, MXNP, NBOX, MXNPABOX, NewNP, WID, NDepth
     real(KINDDF)        :: When, Deltat
     real(KINDDF)        :: Area(4)
     type(DevVec_DF)     :: Depth
     type(DevVec_DF)     :: DDis
     type(DevVec_I)      :: CurNP
     type(DevVec_I)      :: WalkID
     type(DevVec_I)      :: Stat
     type(DevMat_DF)     :: Pos 
     type(DevVec_DF)     :: PreTime 
     type(DevVec_DF)     :: NextTime 
     type(DevVec_I)      :: NextJPath
     type(DevVec_I)      :: NextJump
     type(DevVec_DF)     :: NJump 
   !-- local varibales  
       
       call InsertAtoms_ByHist0(MXNP, NBOX, MXNPABOX, When, Deltat, NewNP, WID, Area(1), Area(2), Area(3), Area(4), NDepth, Depth, DDis, &
                                CurNP, WalkID, Stat, Pos, PreTime, NextTime, NextJPath, NextJump, NJump, gm_DevRandState(DevN))
       return 
  end subroutine InsertAtoms_ByHist0_a
  !********************************************************************************

  !********************************************************************************
  subroutine InsertAtoms_ByHist0_b(When, Deltat, NewNP, WID, Area, NDepth, Depth, DDis, WalkStat)
  !***  PURPOSE:   to insert partitcles in the box according to given depth distribution
  !
  !   INPUT:     When,       the time of inserting the particles
  !              Deltat,     the time interval of activating the partilces      
  !              NewNP,      the number of pariticles to be inserted
  !              WID,        the walker ID to be inserted   
  !              AREA,       the area to insert the particles
  !              NDepth,     the number of grids of depths
  !              Depth,      the depth grids
  !              DDIs,       the integral depth distribution
  !
  !   OUTPUT     WalkStat,     
  !
  !
   implicit none
   !----   DUMMY Variables
     integer                   :: NEWNP, WID, NDepth
     real(KINDDF)              :: When, Deltat, Area(:)
     type(DevVec_DF)           :: Depth
     type(DevVec_DF)           :: DDis
     class(RandWalkerStat_DEV) :: WalkStat
   !-- local varibales  
       
       call InsertAtoms_ByHist0_a(WalkStat%DEVN, WalkStat%mxNP,     WalkStat%NBOX,       WalkStat%mxNPPerBox, When, Deltat, NewNP, WID, Area, NDepth, Depth, DDis, &
                                                 WalkStat%dCurNP,   WalkStat%dWalkStyID, WalkStat%dStat,      WalkStat%dXP, &
                                                 WalkStat%dPreTime, WalkStat%dNextTime,  WalkStat%dNextJPath, WalkStat%dNextJump, WalkStat%dNJump)
       return 
  end subroutine InsertAtoms_ByHist0_b
  !********************************************************************************

  !********************************************************************************
  subroutine InsertAtoms_ByHist0_c(When, Deltat, NewNP, WID, Area, NDepth, hDepth, hDDis, WalkStat)
  !***  PURPOSE:   to insert partitcles in the box according to given depth distribution
  !
  !   INPUT:     When,       the time at which the particles inserted 
  !              Deltat,     the time interval of activating the partilces      
  !              NewNP,      the number of pariticles to be inserted
  !              WID,        the walker ID to be inserted   
  !              AREA,       the area to insert the particles
  !              NDepth,     the number of grids of depths
  !              HDepth,     the depth grids on host
  !              HDDIs,       the integral depth distribution on host
  !
  !   OUTPUT     WalkStat,     
  !
  !
   implicit none
   !----   DUMMY Variables
     integer                   :: NEWNP, WID, NDepth
     real(KINDDF)              :: When, Deltat, Area(:)
     real(KINDDF)              :: hDepth(:)
     real(KINDDF)              :: hDDis(:)
     class(RandWalkerStat_DEV) :: WalkStat
   !-- local varibales  
     type(DevVec_DF) :: dDepth, dDDis
     integer::IDEV
       
       IDEV = m_DEVICES(WalkStat%DEVN)
       call DevAllocate(IDEV,  dDepth, size(hDepth))
       call DevAllocate(IDEV,  dDDis,  size(hDDis))
       call DevCopyIn(hDepth,  dDepth, size(hDepth))
       call DevCopyIn(hDDis,   dDDis,  size(hDDis))

       call InsertAtoms_ByHist0_b(When, Deltat, NewNP, WID, Area, NDepth, dDepth, dDDis, WalkStat)                                                
       call DevDeallocate(dDepth)
       call DevDeallocate(dDDis)
                                         
       return 
  end subroutine InsertAtoms_ByHist0_c
  !********************************************************************************

  !****************************************************************************
  subroutine Sweep0(WalkerState)
  !***  PORPOSE: to clear the particles that have been in reflected or transmittent
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !  
  !     OUTPUT:  WalkerState,    with the particle in relfected or transmittent are remove from the 
  !                              the arrays

    !
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat):: WalkerState
        !--- local variables
        integer::IB, I,J, WID, IP, CURNP, CURNP0

        do IB=1,  WalkerState%NBox
           I      = WalkerState%FirstPIDP(IB)
           IP     = I-1
           CURNP  = 0
           CURNP0 = WalkerState%CurNP(IB)
           do J=1, CURNP0 
               select case(WalkerState%Stat(I))
               case(CP_WALKSTAT_NOTACTIVE)
                   !do nothing
               case(CP_WALKSTAT_REFLECT, CP_WALKSTAT_TRANS)
                    WalkerState%WalkStyID(I)  = CP_WALKTYPE_NOWALK
                    WalkerState%Stat(I)       = CP_WALKSTAT_NOTACTIVE
                    ! we cannot reset the Pretime, NextTime, and NJump and XP
                    ! they will initialize at Insert new particle
               case (CP_WALKSTAT_ACTIVE, CP_WALKSTAT_TRAPPED)
                    IP    = IP + 1
                    CURNP = CURNP + 1
                    WalkerState%Stat(IP)      = WalkerState%Stat(I) 
                    WalkerState%WalkStyID(IP) = WalkerState%WalkStyID(I)
                    WalkerState%PreTime(IP)   = WalkerState%PreTime(I)
                    WalkerState%NextTime(IP)  = WalkerState%NextTime(I)
                    WalkerState%NextJPath(IP) = WalkerState%NextJPath(I)
                    WalkerState%NextJump(IP)  = WalkerState%NextJump(I)
                    WalkerState%NJump(IP)     = WalkerState%NJump(I)  
                    WalkerState%XP(IP,:)      = WalkerState%XP(I,:)
               end select      
              I = I + 1
           end do

           WalkerState%CurNP(IB) = CURNP
       end do   
      return 
  end subroutine Sweep0
  !****************************************************************************

  !****************************************************************************
  subroutine Sweep1(WalkerState, Percent, Discard)
  !***  PORPOSE: to delete some particles with a percentage 
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !              Percent,
  !  
  !     OUTPUT:  WalkerState,    with some of the particle discarded 
  !              Discard,        the number of particles discarded

    !
    use RAND32_MODULE, only:Drand32
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat):: WalkerState
        real(KINDDF)         :: Percent 
        real(KINDDF)         ::Discard(:)
        !--- local variables
        integer::IB, I, J, CURNP0

        do IB=1,  WalkerState%NBox
           I      = WalkerState%FirstPIDP(IB)
           CURNP0 = WalkerState%CurNP(IB)
           do J=I, I + CURNP0 -1
              if( WalkerState%Stat(J) .eq. CP_WALKSTAT_ACTIVE .or. &
                  WalkerState%Stat(J) .eq. CP_WALKSTAT_TRAPPED  ) then
                  if(Drand32() <=  Percent) then
                     WalkerState%Stat(J) = CP_WALKSTAT_NOTACTIVE
                     Discard(WalkerState%WalkStyID(J)) = Discard(WalkerState%WalkStyID(J)) + 1
                  end if   
               end if      
           end do
        end do
        call Sweep0(WalkerState)
     return 
  end subroutine Sweep1
  !****************************************************************************  

  !****************************************************************************
  subroutine Sweep2D(WalkerState, Area, Discard)
  !***  PORPOSE: to delete some particles out of ith their X-Y position out of given Area
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !              Area,           the Area 
  !  
  !     OUTPUT:  WalkerState,    with the particle outised the area are discarded
  !              Discard,        the number of particles discarded
  !
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat):: WalkerState
        real(KINDDF)         :: Area(:) 
        real(KINDDF)         :: Discard(:)
        !--- local variables
        integer::IB, I, J, CURNP0


        do IB=1,  WalkerState%NBox
           I      = WalkerState%FirstPIDP(IB)
           CURNP0 = WalkerState%CurNP(IB)
           do J=I, I + CURNP0 -1
              if( WalkerState%Stat(J) .eq. CP_WALKSTAT_ACTIVE  .or.  &
                  WalkerState%Stat(J) .eq. CP_WALKSTAT_TRAPPED .or. &
                  WalkerState%Stat(J) .eq. CP_WALKSTAT_WAITING    ) then
                  if(WalkerState%XP(J,1) < Area(1)  .or.  &
                     WalkerState%XP(J,1) > Area(2)  .or.  &
                     WalkerState%XP(J,2) < Area(3)  .or.  &
                     WalkerState%XP(J,2) > Area(4)        &
                    ) then
                     WalkerState%Stat(J) = CP_WALKSTAT_NOTACTIVE
                     Discard(WalkerState%WalkStyID(J)) = Discard(WalkerState%WalkStyID(J)) + 1
                  end if   
              end if      
           end do
        end do
        call Sweep0(WalkerState)
     return 
  end subroutine Sweep2D
  !****************************************************************************  

  !****************************************************************************
  subroutine Sweep3D(WalkerState, Vol, Discard)
  !***  PORPOSE: to delete some particles out of ith their X-Y-Z position out of given vol
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !              Vol,            the volume 
  !  
  !     OUTPUT:  WalkerState,    with the particles outside the vol are discarded
  !              Discard,        the number of particles discarded
  !
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat):: WalkerState
        real(KINDDF)         :: Vol(:) 
        real(KINDDF)         ::Discard(:)
        !--- local variables
        integer::IB, I, J, CURNP0

        do IB=1,  WalkerState%NBox
           I      = WalkerState%FirstPIDP(IB)
           CURNP0 = WalkerState%CurNP(IB)
           do J=I, I + CURNP0 -1
              if( WalkerState%Stat(J) .eq. CP_WALKSTAT_ACTIVE  .or. &
                  WalkerState%Stat(J) .eq. CP_WALKSTAT_TRAPPED .or. &
                  WalkerState%Stat(J) .eq. CP_WALKSTAT_WAITING      ) then
                  if(WalkerState%XP(J,1) < Vol(1)  .or.  &
                     WalkerState%XP(J,1) > Vol(2)  .or.  &
                     WalkerState%XP(J,2) < Vol(3)  .or.  &
                     WalkerState%XP(J,2) > Vol(4)  .or.  &
                     WalkerState%XP(J,3) < Vol(5)  .or.  &
                     WalkerState%XP(J,3) > Vol(6)        &
                    ) then
                     WalkerState%Stat(J) = CP_WALKSTAT_NOTACTIVE
                     Discard(WalkerState%WalkStyID(J)) = Discard(WalkerState%WalkStyID(J)) + 1
                  end if   
              end if      
           end do
        end do
        call Sweep0(WalkerState)
     return 
  end subroutine Sweep3D
  !****************************************************************************  

  !****************************************************************************
  subroutine Sweep2(WalkerState, Vol, Discard)
  !***  PORPOSE: to delete some particles out of ith their X-Y-Z position out of given vol
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !              Area,           the Area 
  !  
  !     OUTPUT:  WalkerState,    with the particle in relfected or transmittent are remove from the 
  !                              the arrays
  !
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat):: WalkerState
        real(KINDDF)         :: Vol(:) 
        real(KINDDF)         ::Discard(:)
        !--- local variables
        integer::IB, I, IP, CURNP, CURNP0

        if(size(Vol) .eq. 6 )  then
           call Sweep3D(WalkerState, Vol, Discard)
        else if(size(Vol) .eq. 4) then
           call Sweep2D(WalkerState, Vol, Discard) 
        end if  
     return 
  end subroutine Sweep2
  !****************************************************************************  

  !****************************************************************************
  subroutine Sweep3(WalkerState, CtrlParam)
  !***  PORPOSE: to decrease the number of boxes
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !  
  !     OUTPUT:  WalkerState,    with the number of boxes delete, and the number of particles permitted in 
  !                              in a box is doubled
  !              CtrlParam,      paramters (Nbox, MXNPABOX ) is changed
  !
    use MC_TypeDef_RandWalkCtrl
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat) :: WalkerState
        type(RandWalkCtrl)    :: CtrlParam
        !--- local variables
        integer::HNbox, IB, I, J, IP, CURNP0
        integer, allocatable::FirstPIDP(:),  CurNP(:)   

             if(WalkerState%NBox .le. 1) return
             !--- for Nbox > 1
             !  
             HNbox =  WalkerState%NBox/2

             !--- to set the boxes of even index to initial stat
             allocate(FirstPIDP(HNBOX),  CurNP(HNBOX))
             do IB = 1,  HNbox
                I     =  WalkerState%FirstPIDP(IB*2)
                IP    = I-1
                CURNP0 = WalkerState%CurNP(IB*2)
                do J=1, CURNP0 
                   IP    = IP + 1
                   WalkerState%WalkStyID(IP)   = CP_WALKTYPE_NOWALK
                   WalkerState%Stat(IP)        = CP_WALKSTAT_NOTACTIVE
                   WalkerState%PreTime(IP)     = 0.D0
                   WalkerState%NextTime(IP)    = 0.D0
                   WalkerState%NJump(IP)       = 0.D0
                   WalkerState%XP(IP,:)        = 0.D0                   
                end do
                FirstPIDP(IB) = WalkerState%FirstPIDP(IB*2-1)
                CurNP(IB)     = WalkerState%CurNP(IB*2-1)
             end do   
             deallocate(WalkerState%FirstPIDP, WalkerState%CurNP)
             allocate(WalkerState%FirstPIDP(HNBOX), WalkerState%CurNP(HNBOX))

             !--- to set the boxes of even index to initial stat
             WalkerState%NBox         = HNbox
             WalkerState%MxNPPerBox   = WalkerState%MxNPPerBox*2
             WalkerState%FirstPIDP    = FirstPIDP
             WalkerState%CurNP        = CurNP
             deallocate(FirstPIDP,  CurNP) 
             !--- donot forget update the control parameter
             CtrlParam%MXNPABOX = WalkerState%MxNPPerBox 
             CtrlParam%NBox     = WalkerState%NBox
     return 
  end subroutine Sweep3
  !****************************************************************************  

  !****************************************************************************
  subroutine Sweep4(hWalkerState, dWalkerState, CtrlParam)
  !***  PORPOSE: to decrease the number of boxes
  !        
  !     INPUT:   hWalkerState,    the current state of the boxes
  !  
  !     OUTPUT:  hWalkerState,    with the number of boxes delete, and the number of particles permitted in 
  !                               in a box is doubled
  !              dWalkerState     the device version
  !              CtrlParam,      paramters (Nbox, MXNPABOX ) is changed
  !
    use MC_TypeDef_RandWalkCtrl
    implicit none
        !--- DUMMY variables
       class(RandWalkerStat)     :: hWalkerState
       class(RandWalkerStat_Dev) :: dWalkerState
       type(RandWalkCtrl)        :: CtrlParam
        !--- local variables
            
             if(hWalkerState%NBox .le. 1) return
             !--- for Nbox > 1
             !  
             call Sweep3(hWalkerState, CtrlParam)
             call Copyin_RandWalkerStat1(hWalkerState, dWalkerState)
     return 
  end subroutine Sweep4
  !****************************************************************************  

  !****************************************************************************
  subroutine InBoxParticles0(WalkerState, theNumb)
  !***  PORPOSE: to extract the number of particles remaining in box
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !  
  !     OUTPUT:  theNumb,        the number of particles remaining in the boxes

    !
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat)  :: WalkerState
        integer                :: theNumb
        !--- local variables
        integer::IB, I,J, CURNP

        theNumb = 0
        do IB=1,  WalkerState%NBox
           I     = WalkerState%FirstPIDP(IB)
           CURNP = WalkerState%CurNP(IB)
           do J=1, CURNP
               if(WalkerState%Stat(I) == CP_WALKSTAT_ACTIVE .or. &
                  WalkerState%Stat(I) == CP_WALKSTAT_TRAPPED     ) then
                    theNumb = theNumb + 1 
                end if      
                I = I + 1
           end do
       end do   
      return 
  end subroutine InBoxParticles0
  !****************************************************************************  

  !****************************************************************************
  subroutine InBoxParticles1(WalkerState, theNumb)
  !***  PORPOSE: to extract the number of particles remaining in box
  !        
  !     INPUT:   WalkerState,    the current state of the boxes
  !  
  !     OUTPUT:  theNumb,        the number of particles remaining in the boxes

    !
    implicit none
        !--- DUMMY variables
        class(RandWalkerStat)  :: WalkerState
        integer                :: theNumb(:)
        !--- local variables
        integer::IB, I,J, CURNP

        theNumb = 0
        do IB=1,  WalkerState%NBox
           I     = WalkerState%FirstPIDP(IB)
           CURNP = WalkerState%CurNP(IB)
           do J=1, CURNP
               if(WalkerState%Stat(I) == CP_WALKSTAT_ACTIVE .or. &
                  WalkerState%Stat(I) == CP_WALKSTAT_TRAPPED     ) then
                    theNumb(WalkerState%WalkStyID(I)) = theNumb(WalkerState%WalkStyID(I)) + 1 
                end if      
                I = I + 1
           end do
       end do   
      return 
  end subroutine InBoxParticles1
  !****************************************************************************  

  !********   routines concerning output
  !****************************************************************************
  subroutine Putout_CurConfig_RandWalkerStat(Fhead, WalkerState, Stamp, ExStats)
  !***  PORPOSE: to putout the instant configuaration of a simulation box
  !     INPUT:   Fhead,       the head of the file name storing the configuration
  !              WalkerState, the current stat of walkers
  !              Stamp,  the recording stamp
  !              ID, optional, the ID of atoms to be put out
  !
  use MSM_TYPEDEF_RecordStamp  
  use MSM_TYPEDEF_InputPaser
  use MSM_TYPEDEF_DataPad
  implicit none
      !--- DUMMY variables
      character*(*)                          :: Fhead
      class(RandWalkerStat)                  :: WalkerState
      type(MSMRecordStamp)                   :: Stamp
      type(StatementList),pointer, optional  :: ExStats
      !--- local variables
      integer, parameter::P_TL = 18, P_TAGL = 12
      integer::K, hFile, I,J, NDAT, ID, NCOL, IB, IS, CURNP, TNA
      character*256 ::FNAME
      character*1024::title, fmt, PSTR0, PSTR

      type(DataPad), dimension(:),   allocatable::pDat
      integer,       dimension(:,:), allocatable::DIM
      character*64                              ::TSN, TS, TSFMT
      type(DataPad), pointer                    ::tsDat
      type(DataPad)                             ::TagDat

      !$$--- to calculate the extension of the file
       if(Stamp%IRec(1) .ge. 0) then
          call STRCATI(fname, Fhead, "Cfg.", Stamp%IRec(1), 4)
       else 
         fname = Fhead
       end if   

       !$$--- to open the file
       call AvailableIOUnit(hFile)
       open(hFile, file = fname, status='unknown')

       !*** write out head and box informations
       write(hFile, fmt="(A)") PKW_OUTSTAT_FORMAT20
       call Putout_RecordStamp(hFile, Stamp)
       call InBoxParticles0(WalkerState, TNA)
       write(hFile,FMT="(A, 3(I8,1x))")      "&NATOM    total number of atoms:       ", TNA

       !$$*** to prepare the format of output
        IS    = 0
        NCOL  = 1
        TSFMT = "(1x I4, 1x, I3, 1x, A)"
        title = ""
        fmt   = "("

        title = title(1:IS)//   "&TYPE         "
        IS    = IS        + len("&TYPE         ")
        fmt   = fmt(1:len_trim(fmt))//"I8,3X"
        TS    = "&TYPECOL"
        write(TSN,fmt=TSFMT)  NCOL, 1, '"I"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        title = title(1:IS)//   "POS(LU)  (x)             (y)              (z)         "
        IS    = IS        + len("POS(LU)  (x)             (y)              (z)         ")
        fmt   = fmt(1:len_trim(fmt))//",3(1PE17.8,1X)"
        TS    = "&XYZCOL"
        write(TSN,fmt=TSFMT)  NCOL, 3, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 3

        title = title(1:IS)//   "STATU   "
        IS    = IS        + len("STATU   ")
        fmt   = fmt(1:len_trim(fmt))//",I6,2X"
        TS    = "&STATCOL"
        write(TSN,fmt=TSFMT) NCOL, 1, '"I"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        title = title(1:IS)//   "NJump             "
        IS    = IS        + len("NJump             ")
        fmt   = fmt(1:len_trim(fmt))//",1(1PE17.8,1X)"
        TS    = "&NJUMPCOL"
        write(TSN,fmt=TSFMT) NCOL, 1, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        title = title(1:IS)//   "BOXID     "
        IS    = IS        + len("BOXID     ")
        fmt   = fmt(1:len_trim(fmt))//",I8,2X"
        TS    = "&BOXIDCOL"
        write(TSN,fmt=TSFMT) NCOL, 1, '"I"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1
        !---

        fmt   = fmt(1:len_trim(fmt))//",A)"

        !$$--- to write out colume information
        write(hFile,FMT="(A, 3(1PE12.5,1x))") "!--- Table colume information:"
        do K =1, NumberofData_DataPad(tsDat)
           call GetData_DataPad(K, tsDat, TagDat)
           write(hFile,FMT="(A))")TagDat%Tag(1:len_trim(TagDat%Tag))
        end do

       !$$--- write out box information

       !$$--- write out supplement information
       if(present(ExStats) ) then
          if(associated(ExStats)) then
            write(hFile,FMT="(A, 3(1PE12.5,1x))") "!--- Supplemnetation information:" 
            call Write_StatementList(hFile, ExStats)
          end if
       end if

       !$$--- write out the configures
       write(hFile,FMT="(A, 3(1PE12.5,1x))") "!--- Configure:"
       write(hFile,FMT="(20A))")title(1:len_trim(title))
       do IB=1,  WalkerState%NBox
          I = WalkerState%FirstPIDP(IB)
          J = 1
          CURNP = WalkerState%CurNP(IB)
          do J =1, CURNP
             IS = WalkerState%Stat(I)
             if(IS == CP_WALKSTAT_ACTIVE .or. &
                IS == CP_WALKSTAT_TRAPPED     ) then

                write(hFile,fmt=fmt)WalkerState%WalkStyID(I),    & 
                                    WalkerState%XP(I,1:3),       &
                                    WalkerState%Stat(I),         &
                                    WalkerState%NJump(I),        &
                                   IB

             end if                      
             I = I + 1
          end do
       end do   
       close(hFile)
       if(allocated(pDat)) deallocate(pDat)
       if(allocated(DIM))  deallocate(DIM)
       call Release_DataPad(tsDat)

       return
  end subroutine Putout_CurConfig_RandWalkerStat
  !**************************************************************************** 
  
  !****************************************************************************
  subroutine Archive_Config_0(hFile, WalkerState)
  !***  PORPOSE: to archive the ucrrent stat for restart
  !     INPUT:   hFile,  the I/O unit
  !              SimBox, the simulation box
  !              Stamp,  the recording stamp
  !
  implicit none
      !--- DUMMY variables
      integer                                :: hFile
      class(RandWalkerStat)                  :: WalkerState
      !--- local variables
      integer::IB


       write(hFile)  WalkerState%mxNP, WalkerState%NBox, WalkerState%MxNPPerBox, WalkerState%NWalkStyID
       write(hFile)  WalkerState%FirstPIDP
       write(hFile)  WalkerState%CurNP
       write(hFile)  WalkerState%WalkStyID
       write(hFile)  WalkerState%Stat
       write(hFile)  WalkerState%PreTime
       write(hFile)  WalkerState%NextTime
       write(hFile)  WalkerState%NJump
       write(hFile)  WalkerState%XP
       write(hFile)  WalkerState%NextJPath
       write(hFile)  WalkerState%NextJump

       return
  end subroutine Archive_Config_0
  !****************************************************************************  

  !****************************************************************************
  subroutine Restore_Config_0(hFile, WalkerState)
  !***  PORPOSE: to archive the ucrrent stat for restart
  !     INPUT:   hFile,  the I/O unit
  !              SimBox, the simulation box
  !              Stamp,  the recording stamp
  !
  implicit none
      !--- DUMMY variables
      integer                                :: hFile
      class(RandWalkerStat)                  :: WalkerState
      !--- local variables
      integer::MXNP, NBOX, MXNPABOX, NWalkStyID, ERR
       
       read(hFile)  MXNP, NBOX, MXNPABOX, NWALKSTYID
       if(MXNP .ne. WalkerState%mxNP .or. NWALKSTYID .ne. WalkerState%NWalkStyID .or. &
          NBOX .ne. WalkerState%NBox .or. MXNPABOX   .ne. WalkerState%mxNPPerBox  ) then
         call Release_RandWalkerStat0(WalkerState) 
       end if  
 
       ERR = 0
       if(.not. allocated(WalkerState%FirstPIDP))  allocate(WalkerState%FirstPIDP(NBOX),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%CurNP))      allocate(WalkerState%CurNP(NBOX),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%WalkStyID))  allocate(WalkerState%WalkStyID(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%Stat) )      allocate(WalkerState%Stat(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%PreTime))    allocate(WalkerState%PreTime(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%NextTime))   allocate(WalkerState%NextTime(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%NextJPath))  allocate(WalkerState%NextJPath(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%NextJump))   allocate(WalkerState%NextJump(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%NJump))      allocate(WalkerState%NJump(MXNP),STAT=ERR)
       if(ERR .gt. 0) goto 100
       if(.not. allocated(WalkerState%XP))         allocate(WalkerState%XP(MXNP,3),STAT=ERR)
       if(ERR .gt. 0) goto 100

       WalkerState%NBOX        = NBOX
       WalkerState%mxNPPerBox  = MXNPABOX
       WalkerState%mxNP        = MXNP
       WalkerState%NWalkStyID  = NWALKSTYID
       read(hFile)  WalkerState%FirstPIDP
       read(hFile)  WalkerState%CurNP
       read(hFile)  WalkerState%WalkStyID
       read(hFile)  WalkerState%Stat
       read(hFile)  WalkerState%PreTime
       read(hFile)  WalkerState%NextTime
       read(hFile)  WalkerState%NJump
       read(hFile)  WalkerState%XP
       read(hFile)  WalkerState%NextJPath
       read(hFile)  WalkerState%NextJump

       return

  100  continue 
       write(*,*) "MCLIB Error:: fail allocating memory for host worksapce of WalkerStat"
       write(*,*) "              MXNP =", MXNP
       stop   

  end subroutine Restore_Config_0
  !****************************************************************************  

 end module  MC_TypeDef_RandWalkerState 
 !*****************************************************************************

