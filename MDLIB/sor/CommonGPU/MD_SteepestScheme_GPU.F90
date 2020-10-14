  module MD_SteepestScheme_GPU
  !**** DESCRIPTION: This module provides routines to search the minima position on
  !                  potential surface along a given direction. This module coulbe be used
  !                  in steepest descent or CG of search minima on energy surface         
  !
  !                  ______________________________________________________
  !**** HISTORY:     2018-12-05(HOU Qing),
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_Globle_Variables_GPU
  implicit none
  !**********************
      private::  Do_Steepest0_Forsteps_DEV,  &
                 Do_Steepest1_Forsteps_DEV 
      public::   Do_Steepest_Forsteps_DEV    
      
  contains
  !****************************************************************************************
  subroutine Do_Steepest0_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
  !***  PORPOSE:   to use steepest method to relax the sytem to local minima
  !
  !     INPUT:     SimBox,     the simulation box
  !                CtrlParam,  the control parameters
  !                ForceClass, the force calculation engine
  !                MXNUMSTEPS, the MXNUMSTEPS for searching the miniam
  !
  !     OUTPUT:    dm_XP:      the configurations on devices, 
  !
  use MD_Forceclass_Register_GPU
  implicit none
  !----   DUMMY Variables
         type(SimMDBox),dimension(:)       :: SimBox
         type(SimMDCtrl),       intent(in) :: CtrlParam
         type(MDForceClassGPU), intent(in) :: ForceClass
         integer,               intent(in) :: MXNUMSTEPS
         !----  Local variables
          type(DevMat_DF), dimension(m_MXDEVICE)::dPreFP, dDFP, dDXP
          type(DevVec_DF), dimension(m_MXDEVICE)::dEPOT0, dEPOT1

          integer::I, IFLAG
          real(KINDDF)::MAXDIS, MAXMOVE, MINDIS, ALPHA, DOTDF, DOTDXDF, DELEPOT, MINEPOT
          real(KINDDF)::LB(3), HB(3)

     
            !---  prepair the parameters
             MAXDIS  = CtrlParam%STEEPEST_MxStep*SimBox(1)%RR 
             MINDIS  = CtrlParam%STEEPEST_MiStep*SimBox(1)%RR
             MINEPOT = CtrlParam%STEEPEST_MiDelE*CP_EV2ERG
             LB      = SimBox(1)%BOXLOW
             HB      = SimBox(1)%BOXUP
             do I=1, 3
                if(CtrlParam%IFPD(I) .eq. 0) then
                  LB(I) = -1.D108
                  HB(I) =  1.D108
                end if  
             end do  
             ALPHA   = CtrlParam%STEEPEST_Alpha

            !--- allocate and initialze the working space
             call DevAllocate(dPreFP, (/m_NAPDEV,3/))
             call DevAllocate(dDFP,   (/m_NAPDEV,3/))
             call DevAllocate(dDXP,   (/m_NAPDEV,3/))
             call DevAllocate(dEPOT0,   m_NAPDEV    )
             call DevAllocate(dEPOT1,   m_NAPDEV    )

             !---  for the initial step
             call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
             call DevMakeCopy_noshift(dm_WorkSpace%FP, dPreFP)
             call DevMultiply_noshift(ALPHA, dPreFP, dDXP)
             call DevMaxAbsval_noshift(dDXP, MAXMOVE)
             if(MAXMOVE .gt. MAXDIS ) then
                call DevMultiply_noshift(MAXDIS/MAXMOVE, dDXP)
             else if(MAXMOVE .le. MINDIS) then
                 !--- convergence reached at the first step
                 !--- clear the allocated working space
                 call DevDeallocate(dPreFP)
                 call DevDeallocate(dDFP)
                 call DevDeallocate(dDXP)
                 call DevDeAllocate(dEPOT0)
                 call DevDeAllocate(dEPOT1)
                 write(*,fmt="(A, I)") " MDPSCU Message: steepest converge at the first steps "
                return    
             end if   
             call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
             call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)
          
             !--- begin next steps
             IFLAG = 0
             call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)
             call Synchroniz_XP_on_Devices()
               
             do I = 1, MXNUMSTEPS
                call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)

                !--- determine the alpha 
                call DevMinus_noshift(dPreFP, dm_WorkSpace%FP, dDFP)
                !call DevDot_noshift(dDXP, dDXP, DOTDX)
                call DevDot_noshift(dDXP, dDFP, DOTDXDF)
                call DevDot_noshift(dDFP, dDFP, DOTDF)
                ALPHA = DOTDXDF/DOTDF
                !ALPHA = DOTDX/DOTDXDF
                if(ALPHA .lt. 0.D0) then
                   ALPHA   = CtrlParam%STEEPEST_Alpha
                end if   
                !--- make displacement of the atoms
                call DevMultiply_noshift(ALPHA, dm_WorkSpace%FP, dDXP)  
                call DevMaxAbsval_noshift(dDXP, MAXMOVE)  
                if(MAXMOVE .gt. MAXDIS ) then
                  call DevMultiply_noshift(MAXDIS/MAXMOVE, dDXP)
                end if 
                call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)  
                call Synchroniz_XP_on_Devices() 
                
                !--- to check distance criteria
                !call DevMaxAbsval_noshift(dDXP, MAXMOVE) 
                if(MAXMOVE .le. MINDIS) then
                  IFLAG = I 
                  exit
                end if     
                !--- to check energy-difference criteria
                call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
                call DevMinus_noshift(dm_WorkSpace%EPOT, dEPOT0, dEPOT1)
                call DevMaxAbsval_noshift(dEPOT1, DELEPOT)   
                if(DELEPOT .le. MINEPOT) then
                  IFLAG = I
                  exit
                 end if
                 call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)

                !--- save the current force
                call DevMakeCopy_noshift(dm_WorkSpace%FP,dPreFP)
             end do  

             !--- clear the allocated working space
             call DevDeallocate(dPreFP)
             call DevDeallocate(dDFP)
             call DevDeallocate(dDXP)
             call DevDeallocate(dEPOT0)
             call DevDeallocate(dEPOT1)
      
             !--- if not converged , issue a wrning message
             if(IFLAG .eq. 0) then
                write(*,fmt="(A, I)")       " MDPSCU WARNING: steepest finished after max steps: ", I-1 
                write(*,fmt="(A, 1PE13.4)") "                 with the max movement (LU) of atoms:", MAXMOVE/SimBox(1)%RR
               call ONWARNING(gm_OnWarning)
             else 
               write(*,fmt="(A, I)")       " MDPSCU Message: steepest finished after steps: ", IFLAG 
               write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev): ", DELEPOT*CP_ERG2EV
             end if  
          return
  end subroutine Do_Steepest0_Forsteps_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Steepest1_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
  !***  PORPOSE:   to use steepest & linesearch method to relax the sytem to local minima
  !
  !     INPUT:     SimBox,     the simulation box
  !                CtrlParam,  the control parameters
  !                ForceClass, the force calculation engine
  !                MXNUMSTEPS, the MXNUMSTEPS for searching the miniam
  !
  !     OUTPUT:    dm_XP:      the configurations on devices, 
  !
   use MD_Forceclass_Register_GPU
  implicit none
  !----   DUMMY Variables
         type(SimMDBox),dimension(:)       :: SimBox
         type(SimMDCtrl),       intent(in) :: CtrlParam
         type(MDForceClassGPU), intent(in) :: ForceClass
         integer,               intent(in) :: MXNUMSTEPS
         !----  Local variables
          type(DevMat_DF), dimension(m_MXDEVICE)::dDXP, dDirect
          type(DevVec_DF), dimension(m_MXDEVICE)::dEPOT0, dEPOT1
          integer::I, ITER
          real(KINDDF)::MAXDIS, MINDIS, PF0, PF1, STEPSIZE, MINEPOT, DELEPOT, A2
          real(KINDDF)::LB(3), HB(3)

            !--- allocate and initialze the working space
             call DevAllocate(dDirect,  (/m_NAPDEV,3/))
             call DevAllocate(dDXP,     (/m_NAPDEV,3/))
             call DevAllocate(dEPOT0,     m_NAPDEV    )
             call DevAllocate(dEPOT1,     m_NAPDEV    )
     
            !---  prepair the parameters 
             MAXDIS  = CtrlParam%STEEPEST_MxStep*SimBox(1)%RR 
             MINDIS  = CtrlParam%STEEPEST_MiStep*SimBox(1)%RR
             MINEPOT = 0.001*CP_EV2ERG
             LB      = SimBox(1)%BOXLOW
             HB      = SimBox(1)%BOXUP
             do I=1, 3
                if(CtrlParam%IFPD(I) .eq. 0) then
                  LB(I) = -1.D108
                  HB(I) =  1.D108
                end if  
             end do  
             
             !----
             ITER = 0
             do while(ITER .le. MXNUMSTEPS)
                
                call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
                call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)
                !--- determine the normalized direction
                call DevMakeCopy_noshift(dm_WorkSpace%FP, dDirect)
                call DevNormalize(dDirect)
                
                !--- set displacement for the next step
                STEPSIZE = MAXDIS
                call DevMultiply_noshift(STEPSIZE, dDirect, dDXP)

                !--- to begine the line search
                call DevDot_noshift(dm_WorkSpace%FP, dDirect, PF0)
                do while(ITER .le. MXNUMSTEPS)
                   !--- move atoms forward for one step
                   call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)
                   call Synchroniz_XP_on_Devices()
                   call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                   call DevDot_noshift(dm_WorkSpace%FP, dDirect, PF1)

                   STEPSIZE = -STEPSIZE*PF1/(PF1 - PF0)
                   if(dabs(STEPSIZE) .gt. MAXDIS) then
                     STEPSIZE = MAXDIS*dabs(STEPSIZE)/STEPSIZE
                   end if  

                   call DevMultiply_noshift(STEPSIZE, dDirect, dDXP)
                   ITER = ITER + 1
                   if(dabs(STEPSIZE) .le. MINDIS) then
                      exit
                   end if
                   PF0 = PF1 
                 end do  
                 !--- 
                 call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
                 call DevMinus_noshift(dm_WorkSpace%EPOT, dEPOT0, dEPOT1)
                 call DevMaxAbsval_noshift(dEPOT1, DELEPOT)
                 if(DELEPOT .le. MINEPOT) then
                    exit
                 end if
             end do  

             !--- clear the allocated working space
             call DevDeallocate(dDirect)
             call DevDeallocate(dDXP)
             call DevDeallocate(dEPOT0)
             call DevDeallocate(dEPOT1)
      
             !--- if not converged , issue a wrning message
             if(ITER .gt. MXNUMSTEPS) then
                write(*,fmt="(A, I)")       " MDPSCU WARNING: steepest finished after max steps: ", ITER 
                write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev): ", DELEPOT*CP_ERG2EV
               call ONWARNING(gm_OnWarning)
             else 
               write(*,fmt="(A, I)")      " MDPSCU Message: steepest with LS finished after steps:  ", ITER 
               write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev): ", DELEPOT*CP_ERG2EV
             end if  
          return
  end subroutine Do_Steepest1_Forsteps_DEV
  !****************************************************************************************  

  !****************************************************************************************
  subroutine Do_Steepest_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS, METH)
  !***  PORPOSE:   to use steepest & linesearch method to relax the sytem to local minima
  !
  !     INPUT:     SimBox,     the simulation box
  !                CtrlParam,  the control parameters
  !                ForceClass, the force calculation engine
  !                MXNUMSTEPS, the MXNUMSTEPS for searching the miniam
  !
  !     OUTPUT:    dm_XP:      the configurations on devices, 
  !
   use MD_Forceclass_Register_GPU
  implicit none
  !----   DUMMY Variables
         type(SimMDBox),dimension(:)       :: SimBox
         type(SimMDCtrl),       intent(in) :: CtrlParam
         type(MDForceClassGPU), intent(in) :: ForceClass
         integer,               intent(in) :: MXNUMSTEPS, METH
         !----  Local variables
            call Synchroniz_XP_on_Devices()
            if(iand(METH, CP_DAMPSCHEME_LSEARCH) .gt. 0) then
               call Do_Steepest1_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
            else 
               call Do_Steepest0_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)  
            end if 
            call Synchroniz_XP_on_Devices()              
          return
  end subroutine Do_Steepest_Forsteps_DEV
  !****************************************************************************************    

  end module MD_SteepestScheme_GPU
