  module MD_CGScheme_GPU
  !**** DESCRIPTION: The GPU version of conjugated gradient method to search the 
  !                  ______________________________________________________
  !**** HISTORY:     2018-11 (HOU Qing), created the first version
  !
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_Globle_Variables_GPU
  implicit none
  !**********************


  contains
  !****************************************************************************************
  subroutine Do_CG0_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
  !***  PORPOSE:   to use CG method to relax the sytem to local minima
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
         type(SimMDBox),dimension(:)           :: SimBox
         type(SimMDCtrl),            intent(in):: CtrlParam
         type(MDForceClassGPU),      intent(in):: ForceClass
         integer,                    intent(in):: MXNUMSTEPS
         !----  Local variables
          type(DevMat_DF), dimension(m_MXDEVICE)::dDir0, dF0, dDXP
          type(DevVec_DF), dimension(m_MXDEVICE)::dEPOT0, dEPOT1

          integer::ITER, IFLAG
          real(KINDDF)::MAXDIS, MAXMOVE, MINDIS, DELEPOT, MINEPOT, F0NORM, MF1NORM, GAMA, STEPSIZE, PF0, PF1
          real(KINDDF)::LB(3), HB(3)
          real(KINDDF)::EPS =1.D-64

     
            !---  prepair the parameters
             MAXDIS  = CtrlParam%STEEPEST_MxStep*SimBox(1)%RR 
             MINDIS  = CtrlParam%STEEPEST_MiStep*SimBox(1)%RR
             MINEPOT = CtrlParam%STEEPEST_MiDelE*CP_EV2ERG
             LB      = SimBox(1)%BOXLOW
             HB      = SimBox(1)%BOXUP
             do ITER=1, 3
                if(CtrlParam%IFPD(ITER) .eq. 0) then
                  LB(ITER) = -1.D108
                  HB(ITER) =  1.D108
                end if  
             end do  

            !--- allocate and initialze the working space
             call DevAllocate(dDir0,  (/m_NAPDEV,3/))
             call DevAllocate(dF0,    (/m_NAPDEV,3/))
             call DevAllocate(dDXP,   (/m_NAPDEV,3/))
             call DevAllocate(dEPOT0,   m_NAPDEV    )
             call DevAllocate(dEPOT1,   m_NAPDEV    )

             !---  for the initial two steps
             !     at the first position 
             call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
             call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
             call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)
             call DevMakeCopy_noshift(dm_WorkSpace%FP, dDir0)
             call DevMakeCopy_noshift(dm_WorkSpace%FP, dF0)
             call DevDot_noshift(dF0, dF0, F0NORM)
             call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF0)
             if(F0NORM .le. EPS) then
                !--- for alread convergence
                call DevDeallocate(dDir0)
                call DevDeallocate(dF0)
                call DevDeallocate(dDXP)
                call DevDeAllocate(dEPOT0)
                call DevDeAllocate(dEPOT1)
                write(*,fmt="(A, I)") " MDPSCU Message: CG converge at the first steps:  "
                return
             end if   
          
             !--- begin next steps
             do ITER = 1, MXNUMSTEPS
                 !--- move to the second position
                 STEPSIZE = MAXDIS/dsqrt(F0NORM)
                 call DevMultiply_noshift(STEPSIZE, dDir0, dDXP)
                 call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)  
                 call Synchroniz_XP_on_Devices() 
                 call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                 call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF1)
 
                !--- move in Dir0 with STEPSIZE
                STEPSIZE = -STEPSIZE*PF1/(PF1 - PF0)
                if(dabs(STEPSIZE) .gt. MAXDIS) then
                  STEPSIZE = MAXDIS*dabs(STEPSIZE)/STEPSIZE
                end if  
                call DevMultiply_noshift(STEPSIZE/dsqrt(F0NORM), dDir0, dDXP)
                call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)  
                call Synchroniz_XP_on_Devices()                 

                !--- to check convergence condition
                call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
                call DevMinus_noshift(dm_WorkSpace%EPOT, dEPOT0, dEPOT1)
                call DevMaxAbsval_noshift(dEPOT1, DELEPOT)
                if(DELEPOT .le. MINEPOT) then
                   exit
                end if
 
                !--- convergence not reached
                call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)
                !--- to determine GAMA for next direction
                call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                !--- using Polak-Ribiere form for gama
                call DevMinus_noshift(dm_WorkSpace%FP, dF0)
                call DevDot_noshift(dm_WorkSpace%FP, dF0, MF1NORM)
                GAMA = MF1NORM/F0NORM
                !--- create new direction:
                call DevMultiply_noshift(GAMA, dDir0) 
                call DevAdd_noshift(dm_WorkSpace%FP, dDir0)  
                call DevMakeCopy_noshift(dm_WorkSpace%FP, dF0)
                call DevDot_noshift(dF0, dF0, F0NORM)
                call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF0)
                if(F0NORM .le. EPS) then
                   exit
                end if   
             end do  

             !--- clear the allocated working space
             call DevDeallocate(dDir0)
             call DevDeallocate(dF0)
             call DevDeallocate(dDXP)
             call DevDeAllocate(dEPOT0)
             call DevDeAllocate(dEPOT1)
  
             !--- if not converged , issue a wrning message
             if(ITER .gt. MXNUMSTEPS) then
                write(*,fmt="(A, I)")       " MDPSCU WARNING: CG finished after max steps: ", ITER-1 
                write(*,fmt="(A, 1PE13.4)") "                 with max movement (LU) of atoms:", MAXMOVE/SimBox(1)%RR
               call ONWARNING(gm_OnWarning)
             else 
               write(*,fmt="(A, I)") " MDPSCU Message: CG finished after steps: ", ITER 
               write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev): ", DELEPOT*CP_ERG2EV
             end if  
          return
  end subroutine Do_CG0_Forsteps_DEV
  !****************************************************************************************
  
  !****************************************************************************************
  subroutine Do_CG1_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
  !***  PORPOSE:   to use CG method to relax the sytem to local minima
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
         type(SimMDBox),dimension(:)           :: SimBox
         type(SimMDCtrl),            intent(in):: CtrlParam
         type(MDForceClassGPU),      intent(in):: ForceClass
         integer,                    intent(in):: MXNUMSTEPS
         !----  Local variables
          type(DevMat_DF), dimension(m_MXDEVICE)::dDir0, dF0, dDXP
          type(DevVec_DF), dimension(m_MXDEVICE)::dEPOT0, dEPOT1

          integer::ITER, IFLAG
          real(KINDDF)::MAXDIS, MAXMOVE, MINDIS, DELEPOT, MINEPOT, F0NORM, MF1NORM, GAMA, STEPSIZE, PF0, PF1
          real(KINDDF)::LB(3), HB(3)
          real(KINDDF)::EPS =1.D-64

     
            !---  prepair the parameters
             MAXDIS  = CtrlParam%STEEPEST_MxStep*SimBox(1)%RR 
             MINDIS  = CtrlParam%STEEPEST_MiStep*SimBox(1)%RR
             MINEPOT = CtrlParam%STEEPEST_MiDelE*CP_EV2ERG
             LB      = SimBox(1)%BOXLOW
             HB      = SimBox(1)%BOXUP
             do ITER=1, 3
                if(CtrlParam%IFPD(ITER) .eq. 0) then
                  LB(ITER) = -1.D108
                  HB(ITER) =  1.D108
                end if  
             end do  

            !--- allocate and initialze the working space
             call DevAllocate(dDir0,  (/m_NAPDEV,3/))
             call DevAllocate(dF0,    (/m_NAPDEV,3/))
             call DevAllocate(dDXP,   (/m_NAPDEV,3/))
             call DevAllocate(dEPOT0,   m_NAPDEV    )
             call DevAllocate(dEPOT1,   m_NAPDEV    )

             !---  for the initial two steps
             !     at the first position 
             call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
             call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
             call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)
             call DevMakeCopy_noshift(dm_WorkSpace%FP, dDir0)
             call DevMakeCopy_noshift(dm_WorkSpace%FP, dF0)
             call DevDot_noshift(dF0, dF0, F0NORM)
             call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF0)
             if(F0NORM .le. EPS) then
                !--- for alread convergence
                call DevDeallocate(dDir0)
                call DevDeallocate(dF0)
                call DevDeallocate(dDXP)
                call DevDeAllocate(dEPOT0)
                call DevDeAllocate(dEPOT1)
                write(*,fmt="(A, I)") " MDPSCU Message: CG converge at the first steps:  "
                return
             end if   
          
             !--- begin next steps
             ITER = 0
             do while(ITER .le. MXNUMSTEPS)
                !--- move to the second position
                STEPSIZE = MAXDIS/dsqrt(F0NORM)
                call DevMultiply_noshift(STEPSIZE, dDir0, dDXP)
                call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)  
                call Synchroniz_XP_on_Devices() 
                call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF1)
                ITER = ITER + 1 
                !--- do linesearch in DIR0
                do while(ITER .le. MXNUMSTEPS)  
                   STEPSIZE = -STEPSIZE*PF1/(PF1 - PF0)
                   if(dabs(STEPSIZE) .gt. MAXDIS) then
                      STEPSIZE = MAXDIS*dabs(STEPSIZE)/STEPSIZE
                   end if  
                   call DevMultiply_noshift(STEPSIZE/dsqrt(F0NORM), dDir0, dDXP)
                   call DevAdd_shift(LB, HB, dDXP, dm_WorkSpace%XP)  
                   call Synchroniz_XP_on_Devices()  
                   
                   ITER = ITER + 1
                   if(dabs(STEPSIZE) .le. MINDIS) then
                      exit
                   end if
                   PF0 = PF1  
                   call CalForce_ForceClass(SimBox, CtrlParam, ForceClass) 
                   call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF1)
                end do                 

                !--- to check convergence condition
                call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
                call DevMinus_noshift(dm_WorkSpace%EPOT, dEPOT0, dEPOT1)
                call DevMaxAbsval_noshift(dEPOT1, DELEPOT)
                if(DELEPOT .le. MINEPOT) then
                   exit
                end if
 
                !--- convergence not reached
                call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)
                !--- to determine GAMA for next direction
                call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                !--- using Polak-Ribiere form for gama
                call DevMinus_noshift(dm_WorkSpace%FP, dF0)
                call DevDot_noshift(dm_WorkSpace%FP, dF0, MF1NORM)
                GAMA = MF1NORM/F0NORM
                !--- create new direction:
                call DevMultiply_noshift(GAMA, dDir0) 
                call DevAdd_noshift(dm_WorkSpace%FP, dDir0)  
                call DevMakeCopy_noshift(dm_WorkSpace%FP, dF0)
                call DevDot_noshift(dF0, dF0, F0NORM)
                call DevDot_noshift(dm_WorkSpace%FP, dDir0, PF0)
                if(F0NORM .le. EPS) then
                   exit
                end if   
             end do  

             !--- clear the allocated working space
             call DevDeallocate(dDir0)
             call DevDeallocate(dF0)
             call DevDeallocate(dDXP)
             call DevDeAllocate(dEPOT0)
             call DevDeAllocate(dEPOT1)
  
             !--- if not converged , issue a wrning message
             if(ITER .gt. MXNUMSTEPS) then
                write(*,fmt="(A, I)")       " MDPSCU WARNING: CG-LS finished after max steps:  ", ITER-1 
                write(*,fmt="(A, 1PE13.4)") "                 with max movement (LU) of atoms:", MAXMOVE/SimBox(1)%RR
               call ONWARNING(gm_OnWarning)
             else 
               write(*,fmt="(A, I)")       " MDPSCU Message: CG-LS finished after steps:  ", ITER 
               write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev): ", DELEPOT*CP_ERG2EV
             end if  
          return
  end subroutine Do_CG1_Forsteps_DEV
  !****************************************************************************************


  !****************************************************************************************
  subroutine Do_CG_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS, METH)
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
         type(SimMDBox),dimension(:)           :: SimBox
         type(SimMDCtrl),            intent(in):: CtrlParam
         type(MDForceClassGPU),      intent(in):: ForceClass
         integer,                    intent(in):: MXNUMSTEPS, METH
         !----  Local variables
            call Synchroniz_XP_on_Devices() 
            if(iand(METH, CP_DAMPSCHEME_LSEARCH) .gt. 0) then
               call Do_CG1_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
            else 
               call Do_CG0_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)  
            end if
            call Synchroniz_XP_on_Devices()               
          return
  end subroutine Do_CG_Forsteps_DEV
  !****************************************************************************************
  

end module MD_CGScheme_GPU
