module MC_RandWalk_Evolve
   !**** DESCRIPTION: _______________________________________________________________________________________
   !                  The module is used to simulate random walk process. The waiting time ditribution
   !                  follows the exponential distribution  EXP(-t/tau)/tau 
   !
    !**** HISTORY:
   !                  version 1st 2020-05 (Hou Qing, Sichuan university)
   !
   !**********************************************************************************   
   use MC_TypeDef_RandWalkCtrl
   use MC_TypeDef_RandWalker_Base
   use MC_TypeDef_RandWalkerState
   use MC_TypeDef_RandWalkerStatNumb
   use MC_RandWalk_Register
   use MSM_TYPEDEF_RecordStamp

   implicit none

         integer,parameter,         private::CP_ResFact = 2                  
         type(RandWalkerStatNumb),  private::hm_AccuSTATNUMB
         type(RandWalkerStatNumb),  private::hm_CurrSTATNUMB
         type(RandWalkerStat),      private::hm_WALKSTAT
         type(RandWalkerStat_DEV),  private::dm_WALKSTAT

         character(len=256),        private::hm_Recpath     = ""
         type(MSMRecordStamp),      private::hm_Stamp

         !---- the current stat of evoluton
         private::EvolveStat
         type::EvolveStat
            real(KINDDF)::EVOLTIME
            real(KINDDF)::TIMESTEP
            real(KINDDF)::TESTTIME
            real(KINDDF)::TESTSTEP
            real(KINDDF)::RECTIME
            real(KINDDF)::RADTIME
            real(KINDDF)::RADSTEP
            real(KINDDF)::EVENTCKTIME
            real(KINDDF)::EVENTCKSTEP
            integer     ::NEWNP

            real*4      ::PRECPUTIME
            real*4      ::CURCPUTIME
            real*4      ::EVERYCPUTIME = 3600.*4. !-- -the cup time interval in seconds,  for every of which the  configuratios to be archived
            real(KINDDF)::ARCHIVETIME             !--- the time of archive the configurations
            real(KINDDF)::ARCHIVESTEP  = 1000     !--- the time of archive the configurations, default is 1ns
         end type EvolveStat
         type(EvolveStat),  private :: m_TimeStat
         type(RandWalkCtrl),private :: m_SwapCtrl

    contains
    !****************************************************************************
    !****************************************************************************
    subroutine Clear_WalkEvolve()
    !***  PURPOSE:   to release the allocated memory
    !
      implicit none

        call Release_RandWalkerStatNumb(hm_AccuSTATNUMB)
        call Release_RandWalkerStatNumb(hm_CurrSTATNUMB)
        call Release_RandWalkerStat(hm_WALKSTAT)
        call Release_RandWalkerStat(dm_WALKSTAT)
        call Release_RandWalk_GV()
        return
    end subroutine Clear_WalkEvolve         
    !****************************************************************************

    !****************************************************************************
    subroutine CheckInput_WalkEvolve(CtrlParam, Walker, Errflag)
    !***  PURPOSE:   to check the consistent of inputs
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !                Walker     the wlak description 
    !
    !     OUTPUT     CtrlParam, Walker
    !                Errflag,   indicating if the non-consistence existi       
      implicit none

         !--- dummy varioables
         type(RandWalkCtrl) ::CtrlParam
         type(WalkerBase)   ::Walker(:)
         integer            ::Errflag          
         !--- local variables
         integer::IR, IW, MXNP, J
          
         !----
               Errflag = 0 
               !--- associate the rad condition with walkers
               do IR=1, CP_MX_RADSTYLE
                  if(CtrlParam%Radtype(IR) .eq. CP_MC_RAD_NONE) exit
                  do IW = 1, size(Walker) 
                     if(CtrlParam%RadtypeForSymb(IR) .eq. Walker(IW)%Name) then
                        CtrlParam%RadtypeForID(IR) = IW
                     end if  
                  end do
               end do   

              !--- prepare irradition volume
               do IR=1, CtrlParam%RadSorNum
                  CtrlParam%RadVol(1, IR) = -0.5D0*CtrlParam%Boxsize(1)
                  CtrlParam%RadVol(2, IR) =  0.5D0*CtrlParam%Boxsize(1)
                  CtrlParam%RadVol(3, IR) = -0.5D0*CtrlParam%Boxsize(2)
                  CtrlParam%RadVol(4, IR) =  0.5D0*CtrlParam%Boxsize(2)
                  CtrlParam%RadVol(5, IR) =  CtrlParam%BoxLow(3) + CtrlParam%RadPos0(IR,3)
                  CtrlParam%RadVol(6, IR) =  CtrlParam%BoxLow(3) + CtrlParam%RadPos1(IR,3)
               end do    
        return
    end subroutine CheckInput_WalkEvolve
    !****************************************************************************  

    !****************************************************************************
    subroutine Initialize_WalkEvolve(CtrlParam, Walker)
    !***  PURPOSE:   to initialize the module 
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !                
    !     OUTPUT     
      implicit none 
      type(RandWalkCtrl) :: CtrlParam
      type(WalkerBase)   :: Walker(:)

      !--- local variables
       integer::MXNP, NBOX      
       character*256::Fname

          !--- clear the previously allocated memory 
          call Clear_WalkEvolve()

          !---  to initialize the globle constant variables on devices
          call Initializt_RandWalk_GV(CtrlParam, Walker)
          call Register_WalkEngine(CtrlParam, Walker)

          !--- allocate working space for WALKER stat
          call Initialize_RandWalkerStat(hm_WALKSTAT, CtrlParam, Walker)
          call Copyin_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT)  
          
          !--- allocate working space for statistical numbers
          call Initialize_RandWalkerStatNumb(hm_AccuSTATNUMB, CtrlParam, Walker)
          call Initialize_RandWalkerStatNumb(hm_CurrSTATNUMB, CtrlParam, Walker)

          !--- to initialize the I/O units for recording 
          call GetRecpath_RandWalkCtrl(CtrlParam, Fname)
          hm_Recpath = Fname(1:len_trim(Fname))
          call SetRecPath_RandWalkerStatNumb(Fname, hm_AccuSTATNUMB)
          call SetRecPath_RandWalkerStatNumb(Fname, hm_CurrSTATNUMB)

          hm_Stamp%ITime = 0
          hm_Stamp%Time  = 0
          hm_Stamp%ICfg  = 0
          hm_Stamp%IRec  = 0
          hm_Stamp%NBox  = CtrlParam%NBox

        return
    end subroutine Initialize_WalkEvolve
    !****************************************************************************

    !****************************************************************************
    subroutine ConductEvolve_WalkEvolve(CtrlParam, Walker)
    !***  PURPOSE:   to conduct walking simulation
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !
    !     OUTPUT     Walker,    the 
      implicit none

         !--- dummy varioables
         type(RandWalkCtrl) ::CtrlParam
         type(WalkerBase)   ::Walker(:)
         !--- local variables
         real*4::C1,C2                   

               call CPU_TIME(C1)

               !--- to initialize the time stat
               m_TimeStat%EVOLTIME    = 0.D0  
               m_TimeStat%RECTIME     = m_TimeStat%EVOLTIME
               m_TimeStat%RADTIME     = m_TimeStat%EVOLTIME
               m_TimeStat%EVENTCKTIME = m_TimeStat%EVOLTIME    
               m_TimeStat%ARCHIVETIME = m_TimeStat%EVOLTIME 
               m_TimeStat%NEWNP       = 0             
               !--- to initialize the preexist atoms if there are
               !
               call DevCopyIn(hm_WALKSTAT%CurNP, dm_WALKSTAT%dCurNP, hm_WALKSTAT%NBOX)

               !--- make copy of the control parameters, the parameter could be changed
               !    in following evolution, for exaple the surface area
               call Copy_RandWalkCtrl(CtrlParam, m_SwapCtrl)
               if(CtrlParam%Restart > 0) then
                  call Restore_WalkEvolve(m_SwapCtrl%Restart)
                  m_SwapCtrl%RectimeEnd = max(m_SwapCtrl%RectimeEnd, CtrlParam%RectimeEnd) 
               end if   
 
               !--- 
               call ConductWalk_WalkEvolve(CtrlParam, Walker)
               call Archive_WalkEvolve()
               !---
               CALL CPU_TIME(C2)
               print *, "RUN TIME: for ConductWalk_WalkEvolve",C2-C1
               return
    end subroutine ConductEvolve_WalkEvolve
    !****************************************************************************  

    !****************************************************************************
    subroutine ConductWalk_WalkEvolve(CtrlParam, Walker)
    !***  PURPOSE:   to conduct walking simulation
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !                Walker,    the walker description
    !
    !     OUTPUT     Walker,    the 
      use MC_RandWalk_Exp_GPU
      use MC_RandWalk_Pow_GPU
      implicit none

         !--- dummy varioables
         type(RandWalkCtrl) ::CtrlParam
         type(WalkerBase)   ::Walker(:)
         !--- working space
         
         !--- local variables
         real(KINDDF), parameter::EPSFact=0.99999999D0, EPS = 0.00000001D0
         integer::IR, Flag
         !----

            !--- 
            if(m_SwapCtrl%Restart > 0) then
               if(dabs(m_TimeStat%EVOLTIME - m_TimeStat%RECTIME) .le. EPS .and. &
                       m_TimeStat%RECTIME   - hm_Stamp%Time > EPS ) then
                  call DoRecord_WalkEvolve()
               end if
               write(*,fmt="(A, 1PE14.4,A)") "MCPSCU Message: restart the process at time  ", m_TimeStat%EVOLTIME, " (ps)"
            end if      
            !--- to begin evolution()

            !
            call CPU_TIME(m_TimeStat%PRECPUTIME)
            m_TimeStat%CURCPUTIME = m_TimeStat%PRECPUTIME
            do while(.true.)
               !--- to check if we have enough working space for next inserting
               !   if the space is not large enough, we should sweep the space
               call CheckSweep_WalkEvolve(m_SwapCtrl, m_TimeStat%EVOLTIME, m_TimeStat%NEWNP)

               !--- to implement particles
               !    determin the time interval,RADSTEP, in which CtrlParam%RadNPPerstep atoms
               !    to be inserted.
               !    Note: new CtrlParam%RadNPPerstep atoms to be inserted at time
               !          EVOLTIME,  next time for inserting atoms is RADTIME
               call GetRadTimestep_RandWalkCtrl(m_SwapCtrl, m_TimeStat%EVOLTIME, m_TimeStat%RADSTEP)
               call InsertNewPart_WalkEvolve(m_SwapCtrl, Walker, m_TimeStat%EVOLTIME, m_TimeStat%RADSTEP, m_TimeStat%RADTIME, m_TimeStat%NEWNP)
                
               !-- for debug
               !call CopyOut_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 
               
               !--- do event checking
               !
               m_TimeStat%EVENTCKSTEP = 1.D64
               m_TimeStat%EVENTCKTIME = m_TimeStat%EVOLTIME  + m_TimeStat%EVENTCKSTEP

               !---  the time for archive
               m_TimeStat%ARCHIVETIME  = m_TimeStat%EVOLTIME + m_TimeStat%ARCHIVESTEP
                              
               !--- set the time and time step of evolution time to be the min of all event checking
               m_TimeStat%TESTTIME    = min(m_TimeStat%RADTIME,     m_TimeStat%EVENTCKTIME)  
               m_TimeStat%TESTTIME    = min(m_TimeStat%ARCHIVETIME, m_TimeStat%TESTTIME)  

               m_TimeStat%TESTSTEP    = min(m_TimeStat%RADSTEP,     m_TimeStat%EVENTCKSTEP)
               m_TimeStat%TESTSTEP    = min(m_TimeStat%ARCHIVESTEP, m_TimeStat%TESTSTEP)

               !--- to evolve
               do while(m_TimeStat%EVOLTIME .lt. m_TimeStat%TESTTIME    .and. &
                        dabs(m_TimeStat%EVOLTIME - m_TimeStat%TESTTIME)  .gt. EPS)
                  !--- get time for the next recording     
                  !    in case of restarting,  m_TimeStat%RECTIME could be larger than m_TimeStat%EVOLTIME
                  if(m_TimeStat%RECTIME - m_TimeStat%EVOLTIME .lt. EPS ) then      
                     call GetRectime_RandWalkCtrl(m_SwapCtrl, m_TimeStat%RECTIME, Flag)
                  end if   

                  !--- set the time step of evolution time to be the min one of TESTSTEP and the interval 
                  !    of cuurent time to recording time 
                  m_TimeStat%TIMESTEP = min(m_TimeStat%TESTSTEP, m_TimeStat%RECTIME - m_TimeStat%EVOLTIME)

                  if(m_TimeStat%EVOLTIME + m_TimeStat%TIMESTEP .gt. m_SwapCtrl%RectimeEnd) then
                     exit
                  end if

                  do while(.true.)
                     m_TimeStat%EVOLTIME = m_TimeStat%EVOLTIME + m_TimeStat%TIMESTEP
                     if(m_TimeStat%EVOLTIME .ge. m_TimeStat%RECTIME) then
                        m_TimeStat%EVOLTIME = m_TimeStat%RECTIME
                     end if   
                     !--- walk 
                     call Walkto_WalkEvolve(m_SwapCtrl, Walker, m_TimeStat%EVOLTIME)
                     !-- for debug
                     !call CopyOut_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 

                     ! for every 2 hours, we arichive the current state 
                     call CPU_TIME(m_TimeStat%CURCPUTIME)
                     if(m_TimeStat%EVOLTIME .ge. m_TimeStat%ARCHIVETIME .or. &
                        m_TimeStat%CURCPUTIME - m_TimeStat%PRECPUTIME .gt. m_TimeStat%EVERYCPUTIME) then
                        call Archive_WalkEvolve() 
                        m_TimeStat%ARCHIVESTEP = m_TimeStat%ARCHIVESTEP*m_TimeStat%EVERYCPUTIME/(m_TimeStat%CURCPUTIME - m_TimeStat%PRECPUTIME)
                        m_TimeStat%ARCHIVETIME = m_TimeStat%EVOLTIME + m_TimeStat%ARCHIVESTEP
                        m_TimeStat%PRECPUTIME = m_TimeStat%CURCPUTIME
                     end if

                     if(m_TimeStat%EVOLTIME .ge. m_TimeStat%RECTIME*EPSFact ) then
                        exit 
                     end if

                     !--- to checek and implement particles
                     !
                     call CheckSweep_WalkEvolve(m_SwapCtrl, m_TimeStat%EVOLTIME, m_TimeStat%NEWNP)
                     call GetRadTimestep_RandWalkCtrl(m_SwapCtrl, m_TimeStat%EVOLTIME, m_TimeStat%RADSTEP)
                     call InsertNewPart_WalkEvolve(m_SwapCtrl, Walker, m_TimeStat%EVOLTIME, m_TimeStat%RADSTEP, m_TimeStat%RADTIME, m_TimeStat%NEWNP)
      
                     !-- for debug
                     !call CopyOut_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 
                     
                  end do  !---enddo do for(EVOLTIME .le. RECTIME*EPSFact)
                  !--- recording the current stat
                  call DoRecord_WalkEvolve()

                  if(m_TimeStat%EVOLTIME .ge. m_TimeStat%EVENTCKTIME*EPSFact) then
                     exit
                   end if   
               end do !---enddo do while(EVOLTIME .le. TESTTIME)
               
               !---
               if(m_TimeStat%EVOLTIME .ge. m_SwapCtrl%RectimeEnd*EPSFact) then
                  exit
               end if           
            end do   
        return
    end subroutine ConductWalk_WalkEvolve
    !****************************************************************************  

    !****************************************************************************
    subroutine Walkto_WalkEvolve(CtrlParam, Walker, Evoltime)
    !***  PURPOSE:   to invoke walk until a give 
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !                Walker,    the walker description
    !                Evoltime,  the current time of evolution
    !
    !     OUTPUT     updated dm_WALKSTAT
      use MC_RandWalk_Exp_GPU
      use MC_RandWalk_Pow_GPU
      implicit none

         !--- dummy varioables
         type(RandWalkCtrl),  intent(in)    ::CtrlParam
         type(WalkerBase),    intent(in)    ::Walker(:)
         real(KINDDF),        intent(in)    ::Evoltime
         !--- working space
         
         !--- local variables
         integer::IW, J

            !--- to initialize the preexist atoms if there are
            !call Walkto_1D_WalkEvolve(CtrlParam, Walker, Evoltime)
            do IW=1, size(Walker)
               call ConductWalk_RandWalk_GV(Evoltime, IW, dm_WALKSTAT)
            end do                                 
        return
    end subroutine Walkto_WalkEvolve
    !****************************************************************************    

    !****************************************************************************
    subroutine InsertNewPart_WalkEvolve(CtrlParam, Walker, Evoltime, Stepsize, Radtime, NewNP)
    !***  PURPOSE:   to insert new particles in the box
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !                Walker,    the walker description
    !                Evoltime,  the current time of evolution
    !                Stepsize,  the stepsize for next irraidtion
    !                Radtime,   the time at whihc the particles are expected to be inserted   
    !
    !     OUTPUT     Radtime,    the next time at whihc the particles are expected to be inserted  
    !                NewNP,      the number of particles newly added  
    !                dm_WALKSTAT,
      implicit none

         !--- dummy varioables
         type(RandWalkCtrl)    ::CtrlParam
         type(WalkerBase)      ::Walker(:)
         real(KINDDF)          ::Evoltime
         real(KINDDF)          ::Stepsize
         real(KINDDF)          ::Radtime
         integer               ::NewNP
         !--- working space
         
         !--- local variables
         integer::IB, IR, NP(CP_MX_RADSTYLE)
         real(KINDDF), parameter::EPSFact=0.99999999D0
         real(KINDDF)::VOL(6)
            !---
 
            !--- to initialize the preexist atoms if there are
            if(Evoltime .ge. Radtime*EPSFact) then
               call GetRadNPart_RandWalkCtrl(CtrlParam, Evoltime, Stepsize, NP)
               NewNP = sum(NP)
                
               do IR=1, CtrlParam%RadSorNum
                  if(NP(IR) .gt. 0) then
                     VOL(1:6) = CtrlParam%Radvol(1:6,IR)
                     !--- insert the particle acoording depth distribution
                     call InsertAtoms_Uniform(Evoltime, Stepsize, NP(IR), IR, VOL, dm_WALKSTAT)

                  end if   
               end do
               Radtime = Radtime + Stepsize  

               !--- save the CurNP on host for sweep checking
               hm_WALKSTAT%CurNP = hm_WALKSTAT%CurNP + NewNP
               call DevCopyIn(hm_WALKSTAT%CurNP, dm_WALKSTAT%dCurNP, hm_WALKSTAT%NBOX)

              !--- recording the inserted fluence
              !    Note: not all the inserted particles are in wating stat, the real fluence should
              !          be inserted partciles - waiting particles  
               do IB=1, hm_CurrSTATNUMB%NBox
                  hm_CurrSTATNUMB%Inserted(1:CtrlParam%RadSorNum) = hm_CurrSTATNUMB%Inserted(1:CtrlParam%RadSorNum) +      &
                                                 dble(NP(1:CtrlParam%RadSorNum))/(CtrlParam%Boxup(1)-CtrlParam%Boxlow(1)) /&
                                                                                 (CtrlParam%Boxup(2)-CtrlParam%Boxlow(2)) /&
                                                 dble(hm_CurrSTATNUMB%NBox)                                
               end do
            end if   

        return
    end subroutine InsertNewPart_WalkEvolve
    !****************************************************************************  

    !****************************************************************************
    subroutine CheckSweep_WalkEvolve(CtrlParam, Evoltime, NewNP)
    !***  PURPOSE:   to check if we have enough working space to insert new particles
    !                if the residual space is not enough, we should sweep thw space  
    !
    !     INPUT:     Evoltime,   the evolution time  
    !                NewNP,      the number of particles expexted to be inserted
    !
    !     OUTPUT     updated hm_ResNPRT, and hm_CurNPRT if sweep done
      implicit none

         !--- dummy varioables
      type(RandWalkCtrl),intent(inout) ::CtrlParam
      real(KINDDF),      intent(in)    ::Evoltime
      integer,           intent(in)    ::NewNP
      !--- local variables
      integer::ERR=1, I, SHRINK
      real(KINDDF), parameter::FACT2D=1.D0/2.D0**0.5D0, FACT3D=1.D0/2.D0**0.333333333333D0
      save::ERR

         !--call DevCopyOut(hm_WALKSTAT%curNP, dm_WALKSTAT%dcurNP,  hm_WALKSTAT%NBOX)
         if(all(hm_WALKSTAT%mxNPPerBox - hm_WALKSTAT%curNP > NewNP)) then
            return
         end if   
         write(*,fmt="(A, 1PE14.4,A)") "MCPSCU Message: working space filled at ", Evoltime, " (ps)"
         write(*,fmt="(A, 1PE14.4,A)") "                space has to be sweeped "

         !--- need to do sweep
         call CopyOut_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 

         !-- before sweep the space, we record the numbers of particles outof boxes
         call Extract_RandWalkerStatNumb(Evoltime, CtrlParam, hm_WALKSTAT, hm_CurrSTATNUMB)

         !--- accumulate the reflected and tranmitted particle, the reflected and tramitted particle
         !    would be missed when the sweep the working space
         hm_AccuSTATNUMB%Reflected = hm_AccuSTATNUMB%Reflected + hm_CurrSTATNUMB%Reflected
         hm_AccuSTATNUMB%Transmit  = hm_AccuSTATNUMB%Transmit  + hm_CurrSTATNUMB%Transmit
 
         !--- do sweep, deleting the relfected and transmitting particles 
         call Sweep_RandWalkerStat(hm_WALKSTAT)     

         !---  to check the space again
         SHRINK = 0
         do I=1, hm_WALKSTAT%NBox
            if((hm_WALKSTAT%mxNPPerBox - hm_WALKSTAT%curNP(I) < NewNP)) then
               write(*,fmt="(A, I6, 1x, A, 1x, 1PE14.4,A)") "MCPSCU Warning: cannot added new particels to the box ",I, ", at ", Evoltime, " (ps)"
               write(*,fmt="(A, 1PE14.4,A)") "                               try to shrink the box "
               SHRINK = 1
               exit
             end if 
         end do      
  
         !--- if need shrink the box
         if(SHRINK) then
            select case(CtrlParam%Boxtype)
                   case(CP_MC_BOX_INF, CP_MC_BOX_1D_SURF, CP_MC_BOX_1D_DSURF)
                        
                   case (CP_MC_BOX_2D_SURF, CP_MC_BOX_2D_DSURF)
                        if(hm_WALKSTAT%Nbox .le. 1) then
                           CtrlParam%Boxlow(1:2)  =  CtrlParam%Boxlow(1:2)*FACT2D
                           CtrlParam%Boxup(1:2)   =  CtrlParam%Boxup(1:2)*FACT2D
                           CtrlParam%Boxsize      =  CtrlParam%Boxup - CtrlParam%Boxlow
                           call Sweep_RandWalkerStat(hm_WALKSTAT, (/CtrlParam%Boxlow(1), CtrlParam%BoxUp(1),   &
                                                                    CtrlParam%Boxlow(2), CtrlParam%BoxUp(2)/), &
                                                     hm_CurrSTATNUMB%Discarded(:)                              &
                                                    )
                            !--- because the irradiation area has been shrinked by half, the number of inserting particle per step
                            !    has to be decreased by half                                                 
                            CtrlParam%RadNPPerstep = max(1,nint(CtrlParam%RadNPPerstep/2.D0))
                        else   
                            call Sweep_RandWalkerStat(hm_WALKSTAT, CtrlParam)
                        end if

                case (CP_MC_BOX_3D)
                     if(hm_WALKSTAT%Nbox .le. 1) then
                        CtrlParam%Boxlow(1:3) =  CtrlParam%Boxlow(1:3)*FACT3D
                        CtrlParam%Boxup(1:3)  =  CtrlParam%Boxup(1:3)*FACT3D
                        CtrlParam%Boxsize     =  CtrlParam%Boxup - CtrlParam%Boxlow
                        call Sweep_RandWalkerStat(hm_WALKSTAT, (/CtrlParam%Boxlow(1), CtrlParam%BoxUp(1),  &
                                                                 CtrlParam%Boxlow(2), CtrlParam%BoxUp(2),  &
                                                                 CtrlParam%Boxlow(3), CtrlParam%BoxUp(3)/),&
                                                  hm_CurrSTATNUMB%Discarded(:)                             &
                                                 )
                        !--- because the irradiation area has been shrinked by half, the number of inserting particle per step
                        !    has to be decreased by half                                                 
                        CtrlParam%RadNPPerstep = max(1, nint(CtrlParam%RadNPPerstep/2.D0))
                     else
                        call Sweep_RandWalkerStat(hm_WALKSTAT, CtrlParam)   
                     end if                            
            end select

         end if   
   
         call CopyIn_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 
         !--- the number of boxes or the area have been changed
         hm_CurrSTATNUMB%Area  = (CtrlParam%Boxup(1) - CtrlParam%Boxlow(1))*(CtrlParam%Boxup(2) - CtrlParam%Boxlow(2)) 
         hm_CurrSTATNUMB%NBox  =  CtrlParam%NBox
         hm_AccuSTATNUMB%Area  =  hm_CurrSTATNUMB%Area
         hm_AccuSTATNUMB%NBox  =  hm_CurrSTATNUMB%NBox

         !--- for debug
         !call Extract_RandWalkerStatNumb(Evoltime, CtrlParam, hm_WALKSTAT, hm_CurrSTATNUMB)         
        return
    end subroutine CheckSweep_WalkEvolve
    !****************************************************************************      
    
    !****************************************************************************
    subroutine Archive_WalkEvolve()
    !***  PURPOSE:   to archive the current state of evolution
    !
    !     INPUT:     Evoltime,    the evolution time  
    !                theTimeStat, the number of particles expexted to be inserted
    !
    !     OUTPUT     none
      implicit none

         !--- dummy varioables
      !--- local variables
      character*256::Fname
      integer::hFile

          call CopyOut_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 

          Fname = hm_Recpath(1:len_trim(hm_Recpath))//"Archive.dat" 
          call AvailableIOUnit(hFile)
          open(hFile, file = fname, form='unformatted', status='unknown')
          write(hFile) m_TimeStat
          write(hFile) m_SwapCtrl
          call Archive_RecordStamp(hFile, hm_Stamp)
          call Archive_Rand_DEVICES(hFile)
          !--- also make a copy of Inserted particles, this information will be used 
          !    on restarting  
          hm_AccuSTATNUMB%Inserted = hm_CurrSTATNUMB%Inserted
          call Archive_RandWalkerStatNumb(hFile, hm_AccuSTATNUMB)
          call Archive_RandWalkerStat(hFile, hm_WALKSTAT)
          
          close(hFile)
        return
    end subroutine Archive_WalkEvolve
    !****************************************************************************

    !****************************************************************************
    subroutine Restore_WalkEvolve(Restored)
    !***  PURPOSE:   to archive the current state of evolution
    !
    !     INPUT:     Evoltime,    the evolution time  
    !                theTimeStat, the number of particles expexted to be inserted
    !
    !     OUTPUT     Restored,   indicate if successfully restored
      implicit none

      !--- dummy varioables
      integer, intent(out) :: Restored
      !--- local variables
      character*256::Fname
      integer::hFile
      logical::EX

          Fname = hm_Recpath(1:len_trim(hm_Recpath))//"Archive.dat" 
          inquire(FILE=Fname, EXIST=EX)
          if(.not. EX) then
            Restored = 0
            return
          end if  

          call AvailableIOUnit(hFile)
          open(hFile, file = fname, form='unformatted', status='old')

         read(hFile) m_TimeStat
         read(hFile) m_SwapCtrl
         call Restore_RecordStamp(hFile, hm_Stamp)
         call Restore_Rand_DEVICES(hFile)

         call Restore_RandWalkerStatNumb(hFile, hm_AccuSTATNUMB)
         !--- also restore the Inserted particles, this information will be used 
         !    on restarting  
         hm_CurrSTATNUMB%Inserted = hm_AccuSTATNUMB%Inserted
         call Restore_RandWalkerStat(hFile, hm_WALKSTAT)
         close(hFile)

         call CopyIn_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 
         Restored = 1
        return
    end subroutine Restore_WalkEvolve
    !****************************************************************************    

    !****************************************************************************
    subroutine DoRecord_WalkEvolve()
    !***  PURPOSE:   to conduct walking simulation
    !
    !     INPUT:     CtrlParam, the controlling parameter
    !                Walker,    the walker description
    !
    !     OUTPUT     Walker,    the 
      implicit none

         !--- dummy varioables
         !--- working space
         
         !--- local variables
         real(KINDDF), parameter::EPSFact=0.99999999D0, EPS = 0.000000001D0
         !----

                 !--- recording the current stat
                    hm_Stamp%IRec  = hm_Stamp%IRec + 1
                    hm_Stamp%ICfg  = hm_Stamp%IRec
                    hm_Stamp%Time  = m_TimeStat%EVOLTIME
                    call CopyOut_RandWalkerStat(hm_WALKSTAT, dm_WALKSTAT) 

                    !--- to get the statistical number
                    call Extract_RandWalkerStatNumb(m_TimeStat%EVOLTIME, m_SwapCtrl, hm_WALKSTAT, hm_CurrSTATNUMB)

                    !--- accumulate the reflected and tranmitted particle
                    hm_CurrSTATNUMB%Reflected = hm_AccuSTATNUMB%Reflected + hm_CurrSTATNUMB%Reflected
                    hm_CurrSTATNUMB%Transmit  = hm_AccuSTATNUMB%Transmit  + hm_CurrSTATNUMB%Transmit

                    call Putout_RandWalkerStatNumb(hm_CurrSTATNUMB, m_SwapCtrl%Restart)

                    if(m_SwapCtrl%RecDepthBin > 0) &
                       call PutoutDepthHis_RandWalkerStatNumb(hm_CurrSTATNUMB, hm_Stamp)
                       
                    if(m_SwapCtrl%RecConfig .gt. 0) &
                       call Putout_CurConfig_RandWalkerStat(hm_Recpath, hm_WALKSTAT, hm_Stamp)

                    print *, "Do recording at ", hm_Stamp%Time, "Has NP", hm_WALKSTAT%CurNP(1), hm_WALKSTAT%mxNPPerBox
                    if(dabs(sum(hm_CurrSTATNUMB%Inserted) - sum(hm_CurrSTATNUMB%Reflected+ hm_CurrSTATNUMB%Transmit + hm_CurrSTATNUMB%Walking + hm_CurrSTATNUMB%Waiting)) &
                        .gt. EPS) then
                       print *, "Particle number is not conserved"
                       print *, sum(hm_CurrSTATNUMB%Inserted), " vs ", sum(hm_CurrSTATNUMB%Reflected+ hm_CurrSTATNUMB%Transmit + hm_CurrSTATNUMB%Walking + hm_CurrSTATNUMB%Waiting)
                       pause 
                    end if 
        return
    end subroutine DoRecord_WalkEvolve
    !****************************************************************************    
    
end module MC_RandWalk_Evolve

