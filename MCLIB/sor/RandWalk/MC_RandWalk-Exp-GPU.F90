
 module MC_RandWalk_Exp_GPU
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The module is used to simulate random walk process. The waiting time ditribution
 !                  follows the exponential distribution  EXP(-t/tau)/tau 
 !
  !**** HISTORY:
 !                  version 1st 2020-05 (Hou Qing, Sichuan university)
 !
 !**********************************************************************************   
   use MC_RandWalk_GV
   implicit none

      !-------------------------------------
      private:: Stepto_RecTime_1D_Test_KERNEL,   & 
                Stepto_RecTime_1D_Test_a,        &
                Stepto_RecTime_1D_Test_b
      public::  Stepto_RecTime_Exp_1D_Test
      interface Stepto_RecTime_Exp_1D_Test
         module procedure Stepto_RecTime_1D_Test_a
         module procedure Stepto_RecTime_1D_Test_b
      end interface Stepto_RecTime_Exp_1D_Test        

      private:: Stepto_RecTime_1D_KERNEL0,       & 
                Stepto_RecTime_1D_0,             &
                Stepto_RecTime_1D_0_a               
      public::  Stepto_RecTime_1D                !-- we have public Stepto_RecTime_1D to make it accessable by functional pointer           
      public::  Stepto_RecTime_Exp_1D
      interface Stepto_RecTime_Exp_1D
         module procedure Stepto_RecTime_1D_0
         module procedure Stepto_RecTime_1D_0_a
      end interface Stepto_RecTime_Exp_1D        
      !-------------------------------------
      private:: Stepto_RecTime_WithLJ_1D_Test_KERNEL,   & 
                Stepto_RecTime_WithLJ_1D_Test_a,        &
                Stepto_RecTime_WithLJ_1D_Test_b
      public::  Stepto_RecTime_Exp_1D_WithLJ_Test
      interface Stepto_RecTime_Exp_1D_WithLJ_Test
         module procedure Stepto_RecTime_WithLJ_1D_Test_a
         module procedure Stepto_RecTime_WithLJ_1D_Test_b
      end interface Stepto_RecTime_Exp_1D_WithLJ_Test
  contains
  !**********************************************************************************
  !**********************************************************************************
  attributes(global) subroutine Stepto_RecTime_1D_Test_KERNEL(NPart, Tau, JVect, RecTime, PreTime, NextTime, NJump, Pos, RandState)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !$$   INPUT:     
  !$$              NPart,      actuall number of atoms on the device
  !$$              Tau,        Tau
  !$$              JVect,      jumping diaplacement in a jump
  !$$              RecTime,    recording time 
  !$$              PreTime,    time point at which the last jump occour
  !$$              NextTime,   time point at which the next jump to occour
  !$$              NJump,      total number of jumps have occur
  !$$              Pos,        position of the atoms
  !$$              RandState,  random number state
  !$$
  !$$   OUTPUT     PreTime,    updated 
  !$$              NextTime,   updated 
  !$$              NJump,      updated 
  !$$              Pos,        updated 
  !$$              RandState,  updated
  !
  implicit none
  !----   DUMMY Variables
          integer,      value::NPart
          real(KINDDF), value::Tau, JVect,RecTime
          real(KINDDF), device::NJump(*)
          real(KINDDF), device::PreTime(*), NextTime(*), Pos(*)
          type(curandStateXORWOW), device::RandState(*)
          

  !----   Local variables
          integer::IT, IB, IC, BS, GS, IC0, NL, L
          real(KINDDF)::R, ATT, CURSTEP, DT, DSTEP, T1, T2, NJ, XP
  !
  !***
  !     
              !--- size the block and grid
              GS  = griddim%x*griddim%y
              BS  = blockdim%x*blockdim%y
              !--- number of loops needed for cover all atoms
              NL  = (NPART-1)/(GS*BS) + 1

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x

              !IC -- the id of the random number sequence 
              IC0  = (IB-1)*BS + IT

              !IS -- the rand generator ID for this thread on this GPU
              DSTEP = C_TWO*JVect
              
              do L = 1, NL
                 !IC -- shift id of the atom for BS*GS when the total number thread cannot cover all atoms
                 IC = (L-1)*(BS*GS) + IC0
                 !---  
                 if(IC .le. NPART) then
                   T1 = PreTime(IC)
                   T2 = NextTime(IC)
                   NJ = NJUMP(IC)
                   XP = Pos(IC)
                   if(RecTime .gt. T2 ) then
                       if(NJ .gt. 0.D0) then
                          R   = curand_uniform(RandState(IC0))
                          NJ  = NJ + 1
                          T1  = T2
                          XP  = XP + (R-0.5D0)*DSTEP
                       end if 
                    
                       ATT     = 0.D0
                       CURSTEP =  RecTime - T1
                       do while(ATT .le. CURSTEP)
                          R      = curand_uniform(RandState(IC0))
                          DT     = -Tau*dlog(R)      
                          ATT    = ATT + DT
                          if(ATT .gt. CURSTEP) then
                             T2 = T1 + DT
                             exit
                          end if
                          !--- for displacement
                          R    = curand_uniform(RandState(IC0))
                          XP   = XP + (R-0.5D0)*DSTEP
                          NJ   = NJ + 1

                          !--- update the residule time
                          T1      = T1 + DT
                          ATT     = 0.D0
                          CURSTEP = RecTime - T1
                       end do
                       PreTime(IC)  = T1
                       NextTime(IC) = T2
                       NJUMP(IC)    = NJ
                       Pos(IC)      = XP
                   end if   
                 end if      
               end do

         return
  end subroutine Stepto_RecTime_1D_Test_KERNEL
  !********************************************************************************

  !********************************************************************************
  subroutine Stepto_RecTime_1D_Test_a(NPart, Tau, JVect, RecTime, PreTime, NextTime, NJump, Pos, RandState)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !   INPUT:     
  !              NPart,      actuall number of atoms on the device
  !              Tau,        Tau
  !              JVect,      jumping diaplacement in a jump
  !              RecTime,    recording time 
  !              PreTime,    time point at which the last jump occour
  !              NextTime,   time point at which the next jump to occour
  !              NJump,      total number of jumps have occur
  !              Pos,        position of the atoms
  !              RandState,  random number state
  !
  !   OUTPUT     PreTime,    updated 
  !              NextTime,   updated 
  !              NJump,      updated 
  !              Pos,        updated 
  !              RandState,  updated
   implicit none
   !----   DUMMY Variables
     integer           ::NPart
     real(KINDDF)      ::Tau, RecTime
     real(KINDDF)      ::JVect
     type(DevVec_DF)   ::PreTime
     type(DevVec_DF)   ::NextTime
     type(DevVec_DF)   ::NJump
     type(DevVec_DF)   ::Pos 
     type(DevRandState)::RandState
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
       
       call Stepto_RecTime_1D_Test_KERNEL<<<blocks, threads>>>(NPart, Tau, JVect, RecTime, PreTime%Data, NextTime%Data, NJump%Data, Pos%Data, RandState%RandState)
       ERR = cudaSetDevice(CURDEV)
       return 
  end subroutine Stepto_RecTime_1D_Test_a
  !********************************************************************************

  !********************************************************************************
  subroutine Stepto_RecTime_1D_Test_b(IDEV, NPart, Tau, JVect, RecTime, PreTime, NextTime, NJump, Pos)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !   INPUT:     IDEV,       device ID 
  !              NPart,      actuall number of atoms on the device
  !              Tau,        Tau
  !              JVect,      jumping diaplacement in a jump
  !              RecTime,    recording time 
  !              PreTime,    time point at which the last jump occour
  !              NextTime,   time point at which the next jump to occour
  !              NJump,      total number of jumps have occur
  !              Pos,        position of the atoms
  !              RandState,  random number state
  !
  !   OUTPUT     PreTime,    updated 
  !              NextTime,   updated 
  !              NJump,      updated 
  !              Pos,        updated 
  !
    implicit none
    !----   DUMMY Variables
      integer           ::IDev,NPart
      real(KINDDF)      ::Tau, RecTime
      real(KINDDF)      ::JVect
      type(DevVec_DF)   ::PreTime
      type(DevVec_DF)   ::NextTime
      type(DevVec_DF)   ::NJump
      type(DevVec_DF)   ::Pos 
    !---
        call Stepto_RecTime_1D_Test_a(NPart, Tau, JVect, RecTime, PreTime, NextTime, NJump, Pos, gm_DevRandState(IDev))
        return 
   end subroutine Stepto_RecTime_1D_Test_b  
  !********************************************************************************

  !**********************************************************************************
   attributes(global) subroutine Stepto_RecTime_WithLJ_1D_Test_KERNEL(NPart, Tau, Dav, Beta, RecTime, PreTime, NextTime, NJump, Pos, RandState)
   !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
   !
     !                EXP(-t/tau)/tau 
   !
   !                the distance in one-step would be:
   !                Beta*Dav^Beta/(Dav + d)^(1+Beta)
   !                in random direction
   !
   !$$   INPUT:     
   !$$              NPart,      actuall number of atoms on the device
   !$$              Tau,        Tau
   !$$              JVect,      jumping diaplacement in a jump
   !$$              RecTime,    recording time 
   !$$              PreTime,    time point at which the last jump occour
   !$$              NextTime,   time point at which the next jump to occour
   !$$              NJump,      total number of jumps have occur
   !$$              Pos,        position of the atoms
   !$$              RandState,  random number state
   !$$
   !$$   OUTPUT     PreTime,    updated 
   !$$              NextTime,   updated 
   !$$              NJump,      updated 
   !$$              Pos,        updated 
   !$$              RandState,  updated
   !
   implicit none
   !----   DUMMY Variables
           integer,      value::NPart
           real(KINDDF), value::Tau, Dav, Beta, RecTime
           real(KINDDF), device::NJump(*)
           real(KINDDF), device::PreTime(*), NextTime(*), Pos(*)
           type(curandStateXORWOW), device::RandState(*)
           
   !----   Local variables
           integer::IT, IB, IC, BS, GS, IC0, NL, L
           real(KINDDF)::R, ATT, CURSTEP, DT, DSTEP, T1, T2, NJ, XP
   !--- NOTE: when we calculate the pow R**P, if R or P is double precision, the kernel cannot work correctly        
           real(KINDSF)::SR, SBETAINV
   !       
   !***
   !     
               !--- size the block and grid
               GS  = griddim%x*griddim%y
               BS  = blockdim%x*blockdim%y
               !--- number of loops needed for cover all atoms
               NL  = (NPART-1)/(GS*BS) + 1
 
               IT  = (threadidx%y-1)*blockdim%x + threadidx%x
               IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x
 
               !IC -- the id of the random number sequence 
               IC0  = (IB-1)*BS + IT
 
               !IS -- the rand generator ID for this thread on this GPU
               SBETAINV = -1.D0/Beta
 
               do L = 1, NL
                  !IC -- shift id of the atom for BS*GS when the total number thread cannot cover all atoms
                  IC = (L-1)*(BS*GS) + IC0
                  !---  
                  if(IC .le. NPART) then
                    T1 = PreTime(IC)
                    T2 = NextTime(IC)
                    NJ = NJUMP(IC)
                    XP = Pos(IC)
                    if(RecTime .gt. T2 ) then
                        if(NJ .gt. 0.D0) then
                           R   = curand_uniform(RandState(IC0)) - 0.5D0
                           NJ  = NJ + 1
                           T1  = T2
                           !--- step for forward a distance sampleing from a pow-law distribution
                           SR    = curand_uniform(RandState(IC0))
                           SR    = SR**SBETAINV 
                           DSTEP = Dav*(SR-1.D0)
                           if(R .ge. 0.D0) then
                              XP = XP + DSTEP
                           else
                              XP = XP - DSTEP
                           end if      
                        end if 
                     
                        ATT     = 0.D0
                        CURSTEP =  RecTime - T1
                        do while(ATT .le. CURSTEP)
                           !--- sampling for exp WTD
                           !    if R or P is double precision, the kernel cannot work correctly
                           R      = curand_uniform(RandState(IC0))
                           DT     = -Tau*dlog(R)      
 
                           !--- accumulate the time
                           ATT   = ATT + DT
                           if(ATT .gt. CURSTEP) then
                              T2 = T1 + DT
                              exit
                           end if
                           !--- for displacement
                           !--- step for forward a distance sampleing from a pow-law distribution
                           R     = curand_uniform(RandState(IC0)) - 0.5D0
                           SR    = curand_uniform(RandState(IC0))
                           SR    = SR**SBETAINV 
                           DSTEP = Dav*(SR-1.D0)
                           if(R .ge. 0.D0) then
                              XP = XP + DSTEP
                           else
                              XP = XP - DSTEP
                           end if      
                           NJ   = NJ + 1
 
                           !--- update the residule time
                           T1      = T1 + DT
                           ATT     = 0.D0
                           CURSTEP = RecTime - T1
                        end do
                        PreTime(IC)  = T1
                        NextTime(IC) = T2
                        NJUMP(IC)    = NJ
                        Pos(IC)      = XP
                    end if   
                  end if      
                end do
 
          return
   end subroutine Stepto_RecTime_WithLJ_1D_Test_KERNEL
   !********************************************************************************   
 
   !********************************************************************************
   subroutine Stepto_RecTime_WithLJ_1D_Test_a(NPart, Tau, Dav, Beta, RecTime, PreTime, NextTime, NJump, Pos, RandState)
   !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
   !
  !                EXP(-t/tau)/tau 
   !
   !     INPUT:     
   !              NPart,      actuall number of atoms on the device
   !              Tau,        Tau
   !              JVect,      jumping diaplacement in a jump
   !              RecTime,    recording time 
   !              PreTime,    time point at which the last jump occour
   !              NextTime,   time point at which the next jump to occour
   !              NJump,      total number of jumps have occur
   !              Pos,        position of the atoms
   !              RandState,  random number state
   !
   !     OUTPUT   PreTime,    updated 
   !              NextTime,   updated 
   !              NJump,      updated 
   !              Pos,        updated 
   !              RandState,  updated
    implicit none
    !----   DUMMY Variables
      integer           ::NPart
      real(KINDDF)      ::Tau, Dav, Beta, RecTime
      real(KINDDF)      ::JVect
      type(DevVec_DF)   ::PreTime
      type(DevVec_DF)   ::NextTime
      type(DevVec_DF)   ::NJump
      type(DevVec_DF)   ::Pos 
      type(DevRandState)::RandState
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
 
        call Stepto_RecTime_WithLJ_1D_Test_KERNEL<<<blocks, threads>>>(NPart, Tau, Dav, Beta, RecTime, PreTime%Data, NextTime%Data, NJump%Data, Pos%Data, RandState%RandState)
        ERR = cudaSetDevice(CURDEV)
        return 
   end subroutine Stepto_RecTime_WithLJ_1D_Test_a
   !********************************************************************************
 
   !********************************************************************************
   subroutine Stepto_RecTime_WithLJ_1D_Test_b(IDEV, NPart, Tau, Dav, Beta, RecTime, PreTime, NextTime, NJump, Pos)
   !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
   !
   !                EXP(-t/tau)/tau 
   !
   !     INPUT:   IDEV,       ID id number of device  
   !              NPart,      actuall number of atoms on the device
   !              Tau,        Tau
   !              JVect,      jumping diaplacement in a jump
   !              RecTime,    recording time 
   !              PreTime,    time point at which the last jump occour
   !              NextTime,   time point at which the next jump to occour
   !              NJump,      total number of jumps have occur
   !              Pos,        position of the atoms
   !
   !     OUTPUT   PreTime,    updated 
   !              NextTime,   updated 
   !              NJump,      updated 
   !              Pos,        updated 
     implicit none
     !----   DUMMY Variables
       integer           ::IDev,NPart
       real(KINDDF)      ::Tau, Dav, Beta, RecTime
       real(KINDDF)      ::JVect
       type(DevVec_DF)   ::PreTime
       type(DevVec_DF)   ::NextTime
       type(DevVec_DF)   ::NJump
       type(DevVec_DF)   ::Pos 
     !---
         call Stepto_RecTime_WithLJ_1D_Test_a(NPart, Tau, Dav, Beta, RecTime, PreTime, NextTime, NJump, Pos, gm_DevRandState(IDev))
         return 
    end subroutine Stepto_RecTime_WithLJ_1D_Test_b  
   !********************************************************************************

  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  attributes(global) subroutine Stepto_RecTime_1D_KERNEL0(NPart, WID, RecTime, Walker, PreTime, NextTime,  &
                                          WalkID, NJump, Pos, Stat, RandState)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !$$   INPUT:     
  !$$              NPart,      actuall number of atoms on the device
  !$$              WID,        type id of particles that to be invoke in the kernel 
  !$$              Walker,     the walker descriptor
  !$$              RecTime,    recording time 
  !$$              PreTime,    time point at which the last jump occour
  !$$              NextTime,   time point at which the next jump to occour
  !$$              WalkID,     type id of the particles 
  !$$              NJump,      total number of jumps have occur
  !$$              Pos,        position of the atoms
  !$$              RandState,  random number state
  !$$
  !$$   OUTPUT     PreTime,    updated 
  !$$              NextTime,   updated 
  !$$              NJump,      updated 
  !$$              Stat,       updated                                                         
  !$$              Pos,        updated 
  !$$              RandState,  updated                                                         
  !
  implicit none
  !----   DUMMY Variables
          integer,                 value::NPart, WID
          real(KINDDF),            value::RecTime
          type(WalkerBase),        device::Walker
          real(KINDDF),            device::PreTime(*)
          real(KINDDF),            device::NextTime(*)
          integer,                 device::WalkID(*)
          real(KINDDF),            device::NJump(*)
          real(KINDDF),            device::Pos(NPart,*)
          integer,                 device::Stat(*)
          type(curandStateXORWOW), device::RandState(*)


  !----   Local variables
          integer::IT, IB, IC, BS, GS, IC0, NL, L, CSTAT, OUT
          real(KINDDF)::R, ATT, CURSTEP, DT, T1, T2, NJ, X, Y, Z
          real(KINDDF),shared::TAU, SX, SY, SZ, JX, JY, JZ, LBX, LBY, LBZ, HBX, HBY, HBZ
  !
  !***
             
              LBX = dcm_BoxL(1)
              HBX = dcm_BoxU(1)
              LBY = dcm_BoxL(2)   
              HBY = dcm_BoxU(2)
              LBZ = dcm_BoxL(3)  
              HBZ = dcm_BoxU(3)     
              SX  = HBX - LBX
              SY  = HBY - LBY
              SZ  = HBZ - LBZ
          !---
              Tau = Walker%JumpType(1)%WTD_TAU    
              JX  = Walker%JumpType(1)%JP_VEC(1, 1)
              JY  = Walker%JumpType(1)%JP_VEC(2, 1)
              JZ  = Walker%JumpType(1)%JP_VEC(3, 1)
              
          !--- size the block and grid
              GS  = griddim%x*griddim%y
              BS  = blockdim%x*blockdim%y
              !--- number of loops needed for cover all atoms
              NL  = (NPart-1)/(GS*BS) + 1

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x

              !IC -- the id of the random number sequence 
              IC0  = (IB-1)*BS + IT

              do L = 1, NL
                 !IC -- shift id of the atom for BS*GS when the total number thread cannot cover all atoms
                 IC = (L-1)*(BS*GS) + IC0
                 !---  
                 if(IC .le. NPART) then
                   CSTAT =  Stat(IC)
                   if( WalkID(IC) .eq.  WID) then
                   if(CSTAT .eq. CP_WALKSTAT_ACTIVE .or. CSTAT .eq. CP_WALKSTAT_WAITING) then
                     OUT = 0
                     T1  = PreTime(IC)
                     T2  = NextTime(IC)
                     if(RecTime .gt. T2 ) then
                        NJ    = NJUMP(IC)
                        X     = Pos(IC,1)
                        Y     = Pos(IC,2)
                        Z     = Pos(IC,3)
                        CSTAT = CP_WALKSTAT_ACTIVE

                         if(NJ .gt. 0.D0) then
                            R   = curand_uniform(RandState(IC0)) - 0.5D0
                            if(R .ge. 0.D0) then
                               X  = X + JX
                               Y  = Y + JY
                               Z  = Z + JZ
                            else   
                              X  = X - JX
                              Y  = Y - JY
                              Z  = Z - JZ
                            end if 
                            if(X .lt. LBX) then
                               X = X + SX 
                            else if(X .gt. HBX) then
                               X = X - SX
                            end if
                            if(Y .lt. LBY) then
                              Y = Y + SY 
                            else if(Y .gt. HBY) then
                              Y = Y - SY
                            end if       
                            if(Z .lt. LBZ) then
                               CSTAT = CP_WALKSTAT_REFLECT
                               OUT   = 1
                            else if(Z .gt. HBZ) then
                                CSTAT = CP_WALKSTAT_TRANS
                                OUT   = 1
                            end if       
                            NJ  = NJ + 1
 
                            T1  = T2
                         end if !endif(NJ .gt. 0.D0)

                         ATT     = 0.D0
                         CURSTEP =  RecTime - T1
      
                         do while(ATT .le. CURSTEP .and. OUT.eq.0)
                            R      = curand_uniform(RandState(IC0))
                            DT     = -Tau*dlog(R)      
                             !--- accumulate the time
                             ATT    = ATT + DT
                            if(ATT .gt. CURSTEP) then
                               T2 = T1 + DT
                               exit
                            end if
                            !--- for displacement
                            R    = curand_uniform(RandState(IC0)) - 0.5D0
                            if(R .ge. 0.D0) then
                               X  = X + JX
                               Y  = Y + JY
                               Z  = Z + JZ
                            else   
                               X  = X - JX
                               Y  = Y - JY
                               Z  = Z - JZ
                            end if 

                            if(X .lt. LBX) then
                               X = X + SX 
                            else if(X .gt. HBX) then
                               X = X - SX
                            end if
                            if(Y .lt. LBY) then
                               Y = Y + SY 
                            else if(Y .gt. HBY) then
                               Y = Y - SY
                            end if   

                            if(Z .lt. LBZ) then
                               CSTAT = CP_WALKSTAT_REFLECT
                               OUT   = 1
                            else if(Z .gt. HBZ) then
                               CSTAT = CP_WALKSTAT_TRANS
                               OUT   = 1
                            end if       
                            NJ   = NJ + 1
                            !--- update the residule time
                            T1      = T1 + DT
                            ATT     = 0.D0
                            CURSTEP = RecTime - T1
                            if(OUT .gt. 0) then
                               exit
                            end if

                         end do !enddo while(ATT .le. CURSTEP)

                         PreTime(IC)  = T1
                         NextTime(IC) = T2
                         NJUMP(IC)    = NJ
                         Pos(IC,1)    = X
                         Pos(IC,2)    = Y
                         Pos(IC,3)    = Z
                         Stat(IC)     = CSTAT   
                     end if ! endif(RecTime .gt. T2 )  
                   end if  ! endif(Stat(IC) .eq. CP_WALKSTAT_ACTIVE .or. Stat(IC) .eq. CP_WALKSTAT_WAITING)  
                   end if  ! endif( WalkID(IC) .eq.  WID)
                 end if  ! endif(IC .le. NPART)     
              end do

         return
  end subroutine Stepto_RecTime_1D_KERNEL0
  !********************************************************************************

  !********************************************************************************
  subroutine Stepto_RecTime_1D_0(NPart, WID, WalkerType, RecTime, PreTime, NextTime, WalkID, NJump, Pos, Stat, RandState)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !   INPUT:     
  !              NPart,      actuall number of atoms on the device
  !              WID,        type id of particles that to be invoke in the kernel  
  !              BoxSetup,   the box descriptor   
  !              WalkerType, the walker descriptor
  !              RecTime,    recording time 
  !              PreTime,    time point at which the last jump occour
  !              NextTime,   time point at which the next jump to occour
  !              WalkID,     type id of the particles 
  !              NJump,      total number of jumps have occur
  !              Pos,        position of the atoms
  !              RandState,  random number state
  !
  !   OUTPUT     PreTime,    updated 
  !              NextTime,   updated 
  !              NJump,      updated 
  !              Stat,       updated   
  !              Pos,        updated 
  !
   implicit none
   !----   DUMMY Variables
     integer                      ::NPart, WID
     type(WalkerType_Dev)         ::WalkerType
     real(KINDDF)                 ::RecTime
     type(DevVec_DF)              ::PreTime
     type(DevVec_DF)              ::NextTime
     type(DevVec_I)               ::WalkID
     type(DevVec_DF)              ::NJump
     type(DevVec_I)               ::Stat
     type(DevMat_DF)              ::Pos 
     type(DevRandState)           ::RandState
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

       call Stepto_RecTime_1D_KERNEL0<<<blocks, threads>>>(NPart, WID, RecTime, WalkerType%Walker(WID), PreTime%Data, NextTime%Data, &
                               WalkID%Data, NJump%Data, Pos%Data, Stat%Data, RandState%RandState)    
       
       ERR = cudaSetDevice(CURDEV)
       return 
  end subroutine Stepto_RecTime_1D_0 
  !********************************************************************************
 
  !********************************************************************************
  subroutine Stepto_RecTime_1D_0_a(IDEV, NPart, WID, RecTime, PreTime, NextTime, WalkID, NJump, Pos, Stat)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !   INPUT:     IDEV,       device ID
  !              WID,        type id of particles that to be invoke in the kernel  
  !              RecTime,    recording time 
  !              PreTime,    time point at which the last jump occour
  !              NextTime,   time point at which the next jump to occour
  !              WalkID,     type id of the particles 
  !              NJump,      total number of jumps have occur
  !              Pos,        position of the atoms
  !              RandState,  random number state
  !
  !   OUTPUT     PreTime,    updated 
  !              NextTime,   updated 
  !              NJump,      updated 
  !              Stat,       updated   
  !              Pos,        updated 
  !
  !
    implicit none
    !----   DUMMY Variables
      integer           ::IDEV, NPart, WID
      real(KINDDF)      ::RecTime
      type(DevVec_DF)   ::PreTime
      type(DevVec_DF)   ::NextTime
      type(DevVec_I)    ::WalkID
      type(DevVec_DF)   ::NJump
      type(DevVec_I)    ::Stat
      type(DevMat_DF)   ::Pos 
   !---
        call Stepto_RecTime_1D_0(NPart, WID, dgm_WalkerType(Idev), &
                                RecTime, PreTime, NextTime, WalkID, NJump, Pos,Stat, gm_DevRandState(IDev))

        return 
   end subroutine Stepto_RecTime_1D_0_a
  !******************************************************************************** 

  !********************************************************************************
  subroutine Stepto_RecTime_1D(RecTime, WID, theWalkStat)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
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
  !
  !
    implicit none
    !----   DUMMY Variables
      real(KINDDF)            ::RecTime
      integer                 ::WID
      type(RandWalkerStat_DEV)::theWalkStat
   !---
        call Stepto_RecTime_1D_0_a( theWalkStat%DEVN, theWalkStat%mxNP,  WID,    &
                                    RecTime,                                     &
                                    theWalkStat%dPreTime,                        &
                                    theWalkStat%dNextTime,                       & 
                                    theWalkStat%dWalkStyID,                      &
                                    theWalkStat%dNJump,                          &
                                    theWalkStat%dXP,                             &
                                    theWalkStat%dStat)
        return 
   end subroutine Stepto_RecTime_1D 
  !********************************************************************************
   
  !**********************************************************************************
  attributes(global) subroutine Stepto_RecTime_Generic3D_KERNEL0(NPart, WID, RecTime, Walker, PreTime, NextTime, NextJPath, NextJump, &
                                                          WalkID, NJump, Pos, Stat, RandState)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !$$   INPUT:     
  !$$              NPart,      actuall number of atoms on the device
  !$$              WID,        type id of particles that to be invoke in the kernel 
  !$$              Walker,     the walker descriptor
  !$$              RecTime,    recording time 
  !$$              PreTime,    time point at which the last jump occour
  !$$              NextTime,   time point at which the next jump to occour
  !$$              NextJPath,  jump path of the next time
  !$$              NextJump,   identifier of jump vector in the next jump                                                        
  !$$              WalkID,     type id of the particles 
  !$$              NJump,      total number of jumps have occur
  !$$              Pos,        position of the atoms
  !$$              RandState,  random number state
  !$$
  !$$   OUTPUT     PreTime,    updated 
  !$$              NextTime,   updated 
  !$$              NextJPath,  updated 
  !$$              NextJump,   updated 
  !$$              NJump,      updated 
  !$$              Stat,       updated                                                         
  !$$              Pos,        updated 
  !$$              RandState,  updated                                                         
  !
  implicit none
  !----   DUMMY Variables
          integer,                 value ::NPart, WID
          real(KINDDF),            value ::RecTime
          type(WalkerBase),        device::Walker
          real(KINDDF),            device::PreTime(*)
          real(KINDDF),            device::NextTime(*)
          integer,                 device::NextJPath(:)
          integer,                 device::NextJump(:)
          integer,                 device::WalkID(*)
          real(KINDDF),            device::NJump(*)
          real(KINDDF),            device::Pos(NPart,*)
          integer,                 device::Stat(*)
          type(curandStateXORWOW), device::RandState(*)


  !----   Local variables
          integer::IT, IB, IC, BS, GS, IC0, NL, L, CSTAT, OUT, I, JPP, JP
          real(KINDDF)::R, ATT, CURSTEP, DT, T1, T2, NJ, X, Y, Z
          integer,     shared::NPATH
          real(KINDDF),shared::LBX, LBY, LBZ, HBX, HBY, HBZ, SX, SY, SZ

          real(KINDDF),shared::PATHPROB(CP_WALKTYPE_MAXJPPATH)
          integer,     shared::NJVEC(CP_WALKTYPE_MAXJPPATH)
          real(KINDDF),shared::TAU(CP_WALKTYPE_MAXJPPATH)
          real(KINDDF),shared::JX(CP_WALKTYPE_MAXJPPATH,CP_WALKTYPE_MAXJVECT)
          real(KINDDF),shared::JY(CP_WALKTYPE_MAXJPPATH,CP_WALKTYPE_MAXJVECT)
          real(KINDDF),shared::JZ(CP_WALKTYPE_MAXJPPATH,CP_WALKTYPE_MAXJVECT)
  !
  !***
          IT  = (threadidx%y-1)*blockdim%x + threadidx%x
          IB  =  (blockidx%y-1)*griddim%x  + blockidx%x
         
              LBX = dcm_BoxL(1)
              HBX = dcm_BoxU(1)
              LBY = dcm_BoxL(2)   
              HBY = dcm_BoxU(2)
              LBZ = dcm_BoxL(3)  
              HBZ = dcm_BoxU(3)     
              SX  = HBX - LBX
              SY  = HBY - LBY
              SZ  = HBZ - LBZ
          !---
              NPATH = Walker%NJump
              if(IT .le. CP_WALKTYPE_MAXJPPATH) then
                 PATHPROB(IT) = Walker%JumpProbI(IT)
                 Tau(IT)      = Walker%JumpType(IT)%WTD_TAU
                 NJVEC(IT)    = Walker%JumpType(IT)%JP_NVEC
                 do I=1, NJVEC(IT) 
                    JX(IT,I)  = Walker%JumpType(IT)%JP_VEC(1, I)
                    JY(IT,I)  = Walker%JumpType(IT)%JP_VEC(2, I)
                    JZ(IT,I)  = Walker%JumpType(IT)%JP_VEC(3, I)
                 end do          
              end if
              call syncthreads()
              
          !--- size of block and grid
              GS  = griddim%x*griddim%y
              BS  = blockdim%x*blockdim%y
              !--- number of loops needed for cover all atoms
              NL  = (NPart-1)/(GS*BS) + 1

              !IC -- the id of the random number sequence 
              IC0  = (IB-1)*BS + IT

              do L = 1, NL
                 !IC -- shift id of the atom for BS*GS when the total number thread cannot cover all atoms
                 IC = (L-1)*(BS*GS) + IC0
                 !---  
                 if(IC .le. NPART) then
                  CSTAT =  Stat(IC)
                  if( WalkID(IC) .eq.  WID) then
                  if(CSTAT .eq. CP_WALKSTAT_ACTIVE .or. CSTAT .eq. CP_WALKSTAT_WAITING) then
                    OUT = 0
                    T1  = PreTime(IC)
                    T2  = NextTime(IC)
                    if(RecTime .gt. T2 ) then
                       NJ    = NJUMP(IC)
                       X     = Pos(IC,1)
                       Y     = Pos(IC,2)
                       Z     = Pos(IC,3)
                       CSTAT = CP_WALKSTAT_ACTIVE

                        if(NJ .gt. 0.D0) then
                           X = X + JX(NextJPath(IC), NextJump(IC)) 
                           Y = Y + JY(NextJPath(IC), NextJump(IC)) 
                           Z = Z + JZ(NextJPath(IC), NextJump(IC)) 
                           if(X .lt. LBX) then
                              X = X + SX 
                           else if(X .gt. HBX) then
                              X = X - SX
                           end if
                           if(Y .lt. LBY) then
                             Y = Y + SY 
                           else if(Y .gt. HBY) then
                             Y = Y - SY
                           end if       
                           if(Z .lt. LBZ) then
                              CSTAT = CP_WALKSTAT_REFLECT
                              OUT   = 1
                           else if(Z .gt. HBZ) then
                               CSTAT = CP_WALKSTAT_TRANS
                               OUT   = 1
                           end if       
                           NJ  = NJ + 1

                           T1  = T2
                        end if !endif(NJ .gt. 0.D0)

                        ATT     = 0.D0
                        CURSTEP =  RecTime - T1
     
                        do while(ATT .le. CURSTEP .and. OUT.eq.0)
                           !--- for displacement
                           !--- to determine the jump path
                           R   = curand_uniform(RandState(IC0))
                           do I=1, NPATH
                              if(R .le. PATHPROB(I)) then
                                 JPP = I
                              else
                                 exit
                              end if
                           end do 
                           !--- select a jump vector    
                           JP = curand_uniform(RandState(IC0)) *NJVEC(JPP) + 1

                           !---- determine time for the jump
                           R      = curand_uniform(RandState(IC0))
                           DT     = -Tau(JPP)*dlog(R)      
                           !--- accumulate the time
                            ATT    = ATT + DT
                           if(ATT .gt. CURSTEP) then
                              T2 = T1 + DT
                              NextJPath(IC) = JPP
                              NextJump(IC)  = JP
                              exit
                           end if

                           X  = X + JX(JPP, JP) 
                           Y  = Y + JY(JPP, JP) 
                           Z  = Z + JZ(JPP, JP) 

                           if(X .lt. LBX) then
                              X = X + SX 
                           else if(X .gt. HBX) then
                              X = X - SX
                           end if
                           if(Y .lt. LBY) then
                              Y = Y + SY 
                           else if(Y .gt. HBY) then
                              Y = Y - SY
                           end if   

                           if(Z .lt. LBZ) then
                              CSTAT = CP_WALKSTAT_REFLECT
                              OUT   = 1
                           else if(Z .gt. HBZ) then
                              CSTAT = CP_WALKSTAT_TRANS
                              OUT   = 1
                           end if       
                           NJ   = NJ + 1
                           !--- update the residule time
                           T1      = T1 + DT
                           ATT     = 0.D0
                           CURSTEP = RecTime - T1
                           if(OUT .gt. 0) then
                              exit
                           end if

                        end do !enddo while(ATT .le. CURSTEP)

                        PreTime(IC)  = T1
                        NextTime(IC) = T2
                        NJUMP(IC)    = NJ
                        Pos(IC,1)    = X
                        Pos(IC,2)    = Y
                        Pos(IC,3)    = Z
                        Stat(IC)     = CSTAT   
                    end if ! endif(RecTime .gt. T2 )  
                  end if  ! endif(Stat(IC) .eq. CP_WALKSTAT_ACTIVE .or. Stat(IC) .eq. CP_WALKSTAT_WAITING)  
                  end if  ! endif( WalkID(IC) .eq.  WID)
                end if  ! endif(IC .le. NPART)     
              end do

         return
  end subroutine Stepto_RecTime_Generic3D_KERNEL0
  !********************************************************************************
  
  !********************************************************************************
  subroutine Stepto_RecTime_Generic3D_0(NPart, WID, WalkerType, RecTime, PreTime, NextTime, NextJPath, NextJump, &
                                        WalkID, NJump, Pos, Stat, RandState)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !   INPUT:     
  !              NPart,      actuall number of atoms on the device
  !              WID,        type id of particles that to be invoke in the kernel  
  !              BoxSetup,   the box descriptor   
  !              WalkerType, the walker descriptor
  !              RecTime,    recording time 
  !              PreTime,    time point at which the last jump occour
  !              NextTime,   time point at which the next jump to occour
  !              NextJPath,  jump path of the next time
  !              NextJump,   identifier of jump vector in the next jump                                                        
  !              WalkID,     type id of the particles 
  !              NJump,      total number of jumps have occur
  !              Pos,        position of the atoms
  !              RandState,  random number state
  !
  !   OUTPUT     PreTime,    updated 
  !              NextTime,   updated 
  !              NextJPath,  updated
  !              NextJump,   updated
  !              NJump,      updated 
  !              Stat,       updated   
  !              Pos,        updated 
  !
   implicit none
   !----   DUMMY Variables
     integer                      ::NPart, WID
     type(WalkerType_Dev)         ::WalkerType
     real(KINDDF)                 ::RecTime
     type(DevVec_DF)              ::PreTime
     type(DevVec_DF)              ::NextTime
     type(DevVec_I)               ::NextJPath
     type(DevVec_I)               ::NextJump
     type(DevVec_I)               ::WalkID
     type(DevVec_DF)              ::NJump
     type(DevVec_I)               ::Stat
     type(DevMat_DF)              ::Pos 
     type(DevRandState)           ::RandState
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

       call Stepto_RecTime_Generic3D_KERNEL0<<<blocks, threads>>>(NPart, WID, RecTime, WalkerType%Walker(WID), PreTime%Data, NextTime%Data, &
                                             NextJPath%Data, NextJump%Data, WalkID%Data, NJump%Data, Pos%Data, Stat%Data, RandState%RandState)    
       
       ERR = cudaSetDevice(CURDEV)
       return 
  end subroutine Stepto_RecTime_Generic3D_0 
  !********************************************************************************
 
  !********************************************************************************
  subroutine Stepto_RecTime_Generic3D_0_a(IDEV, NPart, WID, RecTime, PreTime, NextTime, NextJPath, NextJump, WalkID, NJump, Pos, Stat)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
  !
  !   INPUT:     IDEV,       device ID
  !              WID,        type id of particles that to be invoke in the kernel  
  !              RecTime,    recording time 
  !              PreTime,    time point at which the last jump occour
  !              NextTime,   time point at which the next jump to occour
  !              NextJPath,  jump path of the next time
  !              NextJump,   identifier of jump vector in the next jump                                                        
  !              WalkID,     type id of the particles 
  !              NJump,      total number of jumps have occur
  !              Pos,        position of the atoms
  !              RandState,  random number state
  !
  !   OUTPUT     PreTime,    updated 
  !              NextTime,   updated 
  !              NextJPath,  updated
  !              NextJump,   updated
  !              NJump,      updated 
  !              Stat,       updated   
  !              Pos,        updated 
  !
  !
    implicit none
    !----   DUMMY Variables
      integer           ::IDEV, NPart, WID
      real(KINDDF)      ::RecTime
      type(DevVec_DF)   ::PreTime
      type(DevVec_DF)   ::NextTime
      type(DevVec_I)    ::NextJPath
      type(DevVec_I)    ::NextJump
      type(DevVec_I)    ::WalkID
      type(DevVec_DF)   ::NJump
      type(DevVec_I)    ::Stat
      type(DevMat_DF)   ::Pos 
   !---
        call Stepto_RecTime_Generic3D_0(NPart, WID, dgm_WalkerType(Idev), RecTime, PreTime, NextTime, NextJPath, NextJump, &
                                WalkID, NJump, Pos,Stat, gm_DevRandState(IDev))

        return 
   end subroutine Stepto_RecTime_Generic3D_0_a
  !******************************************************************************** 

  !********************************************************************************
  subroutine Stepto_RecTime_Generic3D(RecTime, WID, theWalkStat)
  !***  PURPOSE:   to step the atoms forward to recordding time, the WTD of the a jump is:
  !
  !                EXP(-t/tau)/tau 
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
  !
  !
    implicit none
    !----   DUMMY Variables
      real(KINDDF)            ::RecTime
      integer                 ::WID
      type(RandWalkerStat_DEV)::theWalkStat
   !---
        call Stepto_RecTime_Generic3D_0_a( theWalkStat%DEVN, theWalkStat%mxNP,  WID,    &
                                           RecTime,                                     &
                                           theWalkStat%dPreTime,                        &
                                           theWalkStat%dNextTime,                       & 
                                           theWalkStat%dNextJPath,                      & 
                                           theWalkStat%dNextJump,                       & 
                                           theWalkStat%dWalkStyID,                      &
                                           theWalkStat%dNJump,                          &
                                           theWalkStat%dXP,                             &
                                           theWalkStat%dStat)
        return 
   end subroutine Stepto_RecTime_Generic3D
  !****************************************************

  !********************************************************************************   
  end module  MC_RandWalk_Exp_GPU 