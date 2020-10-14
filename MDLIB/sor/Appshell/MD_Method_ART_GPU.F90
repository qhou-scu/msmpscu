  module MD_Method_ART_GPU
  !***  DESCRIPTION:
  !     This module provides application interfaces of Activation-Relaxsing-Technology (ART)
  !
  !                  _________________________________________________________________________
  !**** HISTORY:     2018-07 (HOU Qing), created the first version
  ! 
  !
  !***  The modules included ******************************************
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_RecordStamp
  use MD_TYPEDEF_RecordList
  implicit none

      type(SimMDCtrl),                         private::m_CtrlParamDamp          !parameters control damping operation
      integer,                                 private::m_mxReplicas  = 0        !max number of replicas
      integer,                                 private::m_curReplicas = 0        !current number of replicas
      type(SimMDBox),                          private::m_SimBoxIni(1)           !the initial configuration at 0K, e.g., the RECTANT
      integer,                                 private::m_ITEvent   = 0          !the time step at which the last event occur
      integer,                                 private::m_NumEvent  = 0          !the number of accepted events
                                                                                 !when two replicas have state transition
      real(KINDDF), dimension(:), allocatable, private::m_EPOT0                  !the energy of atoms at initial configuration

      !--- the interface for selector
       private::SELECT_PROC
       abstract interface
          subroutine SELECT_PROC(CtrlParam, SimBox0, NC, IC)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDCtrl)                    ::CtrlParam
          type(SimMDBox)                     ::SimBox0
          integer                            ::NC
          integer, dimension(:), allocatable ::IC
       end subroutine SELECT_PROC
       end interface
       procedure(SELECT_PROC), pointer, private::m_pSelector=>null()

  contains
  !****************************************************************************************
  !****************************************************************************************
  subroutine SetARTSeedProc(extARTSeedPro)
  !***  DESCRIPTION: to set the procedure that used to construct active region
  !
  !--- INPUT: extARTSeedPro,   the external subroutine of create the candidate atom
  !                            for activartion
  !
  implicit none

  !--- dummy variables and subroutines
   optional::extARTSeedPro
   external::extARTSeedPro

  !--- interface to the external routine -------------------
   interface
       subroutine extARTSeedPro(CtrlParam, SimBox0, NC, IC)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          type(SimMDCtrl)                    ::CtrlParam
          type(SimMDBox)                     ::SimBox0
          integer                            ::NC
          integer, dimension(:), allocatable ::IC
       end subroutine extARTSeedPro 
   end interface 
  !--- END INTERFACE --------------------------------

  !--- Local variables
          if(present(extARTSeedPro)) then
             m_pSelector=>extARTSeedPro
          else
             m_pSelector=>null()
          end if

          return
      end subroutine SetARTSeedProc
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Test(SimBox0, CtrlParam, SimBox, Recordlist, J, processid, INICONFIGPROC)
  !***  PORPOSE: to process one sample
  !     INPUT:  SimBox0,            the element sample
  !             CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !             Recordlist          the external recoding routines
  !             J,                  the test ID
  !             processid,          id of the process, to be used if MPI used
  !             INICONFIGPROC,      the subroutine privided by user for initializing the initial configure
  !
  use MD_SimboxArray_GPU
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_DiffScheme_GPU
  use MD_ForceLib_Factory_GPU
  use MD_TYPEDEF_PrintList, only:Do_RestoreProcess

  implicit none
  !--- dummy variables
   type(SimMDBox),                          intent(in)::SimBox0
   type(SimMDCtrl),target                              ::CtrlParam
   type(SimMDBox), dimension(:), allocatable           ::SimBox
   type(RecordProcedureList)                           ::Recordlist
   integer,                                intent(in)::J,processid

   !--- interface for the external process--------------
   optional::                                     INICONFIGPROC
   external::                                     INICONFIGPROC

   interface
       SUBROUTINE INICONFIGPROC(SimBox, CtrlParam, RESTART)
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer                     ::RESTART
       END SUBROUTINE INICONFIGPROC
   end interface
  !------END INTERFACE ---

  !--- Local variables
  integer::ITIME, ITIME0,TESTLOOP0, TESTLOOP, ITEST, JOB, ISECT0, ISECT, TSTATU
  logical::EXIST

  character*256::CFile0
  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2
  type(MDRecordStamp)::STAMP0
  type(SimMDCtrl), pointer::sectCtrlParam, nextCtrlParam

             !**** to determine how many inner-test loop we need
              if(CtrlParam%INDEPBOX) then
                 TESTLOOP0 = 1
                 TESTLOOP  = 1
                 JOB       = J
              else
                 TESTLOOP0 = 1
                 TESTLOOP  = CtrlParam%TOTALBOX/CtrlParam%MULTIBOX
                 JOB = 1
              end if

              !*** if in restart mode, to check if the SWAP file exist
              EXIST = .FALSE.
              if(CtrlParam%RESTART .GT. 0) then
                 if(CtrlParam%INDEPBOX) then
                    !--- for independent box
                    call STRCATI(CFILE0, CtrlParam%f_configure, "_0K_P", processid, 4)
                    call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                    inquire(FILE=CFILE0, exist=EXIST)
                 else
                    !--- for dependent box, find the last one we have calculated
                    do ITEST=TESTLOOP, 1, -1
                       call STRCATI(CFILE0, CtrlParam%f_configure, "_0K_P", processid, 4)
                       call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                       inquire(FILE=CFILE0, exist=EXIST)

                       if(EXIST) then
                          !--- to restart from this box
                          TESTLOOP0 = ITEST
                          JOB = (ITEST-1)+J
                          exit
                       end if
                    end do
                 end if
              end if

              !*** Now we load the initial configuration
              m_SimBoxIni(1) = SimBox0

              !*** to load initial configuration
              ITIME0        = 0
              ISECT0        = 0
              m_NumEvent    = 0
              m_ITEvent     = 0

              if(CtrlParam%RESTART .eq. 0 .or. .NOT.EXIST) then
                 !--- from the very begin
                 if(Simbox0%IniCfgID .eq. C_IZERO) then
                  !--- load initial file from a single file
                    write(*,fmt="(A)")      ' MDPSCU Message: Load initial configuration from single file:'
                    call Initialize_Config_SimBoxArray(m_SimBoxIni, SimBox0%IniConfig, SimBox0%IniCfgFmt)

                 else if(Simbox0%IniCfgID .gt. C_IZERO) then
                  !--- from immediate point of the last simulation
                   call STRCATI(CFile0, SimBox0%IniConfig, "P", processid, 4)
                   call STRCATI(CFile0, CFile0, "_", JOB, 4)
                   call STRCATI(CFile0, CFile0, ".",  Simbox0%IniCfgID, 4)
                   write(*,fmt="(A)")      ' MDPSCU Message: Load initial substrate from multi-files'
                   !--- NOTE: for Parallel-Replica method, all the box start from the same initial configuration
                   !          ref. MD_Method_GenericMD_GPU.F90 for the difference
                    call  Initialize_Config_SimBoxArray(m_SimBoxIni, CFile0, SimBox0%IniCfgFmt, multbox=0)
                 end if

                 !--- now initialize the configuration
                  if(present(INICONFIGPROC))then
                     call INICONFIGPROC(m_SimBoxIni, CtrlParam, RESTART=0)
                  endif
                  if(allocated(m_EPOT0)) deallocate(m_EPOT0)
                  allocate(m_EPOT0(m_SimBoxIni(1)%NPRT))
                  m_EPOT0 = 0.D0
              !*** if we restart from a previous, restore from previous stop point
              else
                 if(present(INICONFIGPROC))then
                    if(CtrlParam%INDEPBOX) then
                       call INICONFIGPROC(m_SimBoxIni, CtrlParam, RESTART=1)
                    else
                       call INICONFIGPROC(m_SimBoxIni, CtrlParam, RESTART=JOB)
                    end if
                 end if
                 if(allocated(m_EPOT0)) deallocate(m_EPOT0)
                 allocate(m_EPOT0(m_SimBoxIni(1)%NPRT))
                 m_EPOT0 = 0.D0

                 print *, "Restore rectant configure from ", CFILE0(1:len_trim(CFILE0))
                 print *, "........ "
                 call Restore_Config_SimMDBox(m_SimBoxIni(1), STAMP0, CFile0, Do_RestoreProcess)
                 m_ITEvent     = STAMP0%ITime
                 m_NumEvent    = STAMP0%ICfg(1)
                 write(*,fmt="(A, I12, A, I3)") " MDPSCU Message: current events number =        ", m_NumEvent
                 write(*,fmt="(A, I12, A, I3)") "                 last event occur at ITIME    = ", m_ITEvent

                 ITIME0        = m_ITEvent
                 ISECT0        = STAMP0%ISect
              end if

              !*** start test loop
              do ITEST =TESTLOOP0, TESTLOOP
                 JOB       = (ITEST-1)+J
                 ITIME     = ITIME0
                 ISECT     = ISECT0

                 !*** to recorde the start time
                 call CPU_TIME(C1)

                 !*** to determine the start time section firstly
                 call GetNext_SimMDCtrl(CtrlParam, ISECT, sectCtrlParam)
                 !*** start loop for time sections
                 TSTATU = CP_TIMECYCLE_CONTINUE  
                 do while(.TRUE.)
                    call For_One_TimeSect(sectCtrlParam, SimBox, Recordlist, JOB, ISECT, ITIME, processid, TSTATU)
                    if(iand(TSTATU,CP_TIMECYCLE_END) .eq. CP_TIMECYCLE_END ) exit

                    call GetNext_SimMDCtrl(sectCtrlParam, 1, nextCtrlParam)
                    if(.not. associated(nextCtrlParam)) exit
                    sectCtrlParam=>nextCtrlParam

                    !*** if active region method to be used,  initialize the module
                    if(iand(sectCtrlParam%AR_METHOD, CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                       call Initialize_ActiveRegion_DEV(SimBox, sectCtrlParam)
                       call ActivateRegion_DEV(SimBox, sectCtrlParam)
                    else
                       call Clear_ActiveRegion_DEV()
                       call Active_All_ActiveRegion_DEV(SimBox, sectCtrlParam)
                    end if
                 end do   !end loop for time section

                 !--- if the box are dependent of each other, prepair the initial
                 !    configuration of next job
                 if(.not.CtrlParam%INDEPBOX) then
                    if(present(INICONFIGPROC))&
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=0)
                    ITIME0    = 0
                    ISECT0    = 0
                 end if

                 call CPU_TIME(C2)
                 print *, "RUN TIME FOR ONE TEST: ",C2-C1
              end do    !end the loop for itest
      return
  end subroutine For_One_Test
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_TimeSect(CtrlParam, SimBox, Recordlist, Job, ISect, ITime, Processid, Tstatu)
  !***  PORPOSE: to process one sample
  !     INPUT:  CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !             Recordlist          the external recoding routines
  !             JOB,                the ID of the current ITEST
  !             ITIME,              the current times step on starting the time section
  !             processid,          id of the process, to be used if MPI used
  !
  !     OUTPUT: TSTATU,             the flag indicating the the whole time cycle should be ended
  !
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU

  implicit none
  !--- dummy variables
   type(SimMDCtrl),target                               ::CtrlParam
   type(SimMDBox), dimension(:), allocatable            ::SimBox
   type(RecordProcedureList)                            ::Recordlist
   integer,                                  intent(in) ::Job
   integer,                                  intent(in) ::ISect
   integer                                              ::ITime
   integer,                                  intent(in) ::Processid
   integer,                                  intent(out)::Tstatu
   !--- Local variables
    character*256::CFile0
    integer::IB, IT, IE, hFile, EXITATEVENT
    type(MDRecordStamp)::STAMP
  !----
              STAMP%AppTYpe = "ART"
              m_mxReplicas  = CtrlParam%MULTIBOX
              m_curReplicas = m_mxReplicas

              !*** prepare the control parameters for damping
              call Copy_SimMDCtrl(CtrlParam, m_CtrlParamDamp)
              m_CtrlParamDamp%DAMPTIME0     = 1
              m_CtrlParamDamp%DAMPTIME1     = CtrlParam%DAMPTIME1
              m_CtrlParamDamp%IVTIME        = 0
              if(m_CtrlParamDamp%DAMPTIME1 .le. 1) then
                  m_CtrlParamDamp%DAMPTIME1 = 1000
              end if

              !*** to prepare the initial stable state
              if(m_ITEvent .le. 0) then
                 write(*,fmt="(A)") ' MDPSCU Message: do initial damping'
                 call Do_Damp(m_SimBoxIni,  m_CtrlParamDamp, REINIT=1)
                 if(.not.allocated(m_EPOT0)) then
                    allocate(m_EPOT0(m_SimBoxIni(1)%NPRT))
                 end if
                 m_EPOT0 = m_SimBoxIni(1)%EPOT

                 !*** save the reactant for restart
                 STAMP%ITest    = JOB
                 STAMP%IBox     = JOB
                 STAMP%ITime    = ITime
                 STAMP%Time     = 0
                 STAMP%ISect    = ISect
                 STAMP%ICfg(1)  = m_ITEvent
                 STAMP%ICfg(2)  = m_ITEvent
                 STAMP%IRec     = STAMP%ICfg
                 call STRCATI(CFILE0, CtrlParam%f_geometry, "_0K_P", processid, 4)
                 call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                 call Putout_Instance_Config_SimBoxArray(CFILE0, m_SimBoxIni, STAMP)
                 !--- save also the initial energy
                 call STRCATI(CFILE0, CtrlParam%f_geometry, "_0K_P", processid, 4)
                 call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                 CFILE0 = CFILE0(1:len_trim(CFILE0))//".EPOT0"
                 call AvailableIOUnit(hFile)
                 open(hFile, file = CFILE0, form='unformatted', status='unknown')
                 write(hFile) JOB, size(m_EPOT0)
                 write(hFile) m_EPOT0
                 close(hFile)
              end if 
              !*** 
              call Create_SimBoxArray(SimBox, m_mxReplicas, m_SimBoxIni(1))
              !*** start the event searching              
              IT = max(ITime, CtrlParam%IT0)
              if(CtrlParam%ExitAtEvent .le. 0) then
                  EXITATEVENT = CtrlParam%IT1
              else
                  EXITATEVENT = CtrlParam%ExitAtEvent
              end if

              do while(IT .le. CtrlParam%IT1)
                 call For_One_Trail_Dyn(CtrlParam, SimBox, Recordlist, Job, IT, processid)
                 if(m_NumEvent .ge. EXITATEVENT) then
                    TSTATU = ior(TSTATU,CP_TIMECYCLE_END)
                    exit
                 end if
                 ITime = IT
              end do

      return
  end subroutine For_One_TimeSect
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Trail_Dyn(CtrlParam, SimBox, Recordlist, ITest, ITime, processid)
  !***  PORPOSE: to process for until one event found
  !     INPUT:  
  !             CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !             Recordlist          the external recoding routines
  !             ITest               the test ID
  !             ITime,              the time step
  !             processid,          id of the process, to be used if MPI used
  !
  use MD_SimboxArray
  use MD_SimboxArray_GPU
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_DiffScheme_GPU
  use MD_ForceLib_Factory_GPU
  use RAND32_MODULE
  implicit none
  !--- dummy variables
   type(SimMDCtrl),target                  ::CtrlParam
   type(SimMDBox), dimension(:)            ::SimBox
   type(RecordProcedureList)               ::Recordlist
   integer                                 ::ITest, ITime
   integer,                      intent(in)::processid


  !--- Local variables
       integer::TSTATU, ITIMEP, IB, NPRT, IA, I, IFPD(3), SAME, CHECKT
       character*256::GFile, CFile0
       type(MDRecordStamp)::STAMP
       real(KINDDF)::BOX(3), HBOX(3), LBOX(3), UBOX(3), TH, H2S2, MXD2
       real(KINDDF)::ALPHA = 1.20D0, DSTEP=0.01*CP_A2CM
       real(KINDDF),        dimension(:), allocatable::TANGENT0, ETAB
       integer,             dimension(:), allocatable::IBT
       type(SimMDBox),      dimension(:), allocatable::SwapSimBox

       !------
                 NPRT  = m_SimBoxIni(1)%NPRT
                 BOX   = m_SimBoxIni(1)%ZL
                 HBOX  = 0.5D0*BOX
                 LBOX  = m_SimBoxIni(1)%BOXLOW
                 UBOX  = m_SimBoxIni(1)%BOXUP
                 IFPD  =  CtrlParam%IFPD
                 ALPHA = CtrlParam%ART_Alpha
                 DSTEP = CtrlParam%ART_StepLen*m_SimBoxIni(1)%RR
                 MXD2  = (DSTEP*CtrlParam%ART_StepTol)*(DSTEP*CtrlParam%ART_StepTol)

                 !*** to activate the system
                 !    NOTE: depending the activation method, the neighbor-list
                 !          is porbably calculated. The system on device
                 !          is created by m_SimBoxIni(1) 
                 call Do_Activation(CtrlParam, m_SimBoxIni(1), SimBox, DSTEP, m_curReplicas)
                   
                 !*** the actually created replicas could be different from MULTIBOX
                 call SetMultiBox_SimMDCtrl(CtrlParam, m_curReplicas)
                 call Create_SimBoxArray(SwapSimBox, m_curReplicas, m_SimBoxIni(1))

                 !*** to because the number of replicas could be changed, we need 
                 !    to reinitialize the GPU variable
                 call Initialize_Globle_Variables_DEV(SimBox(1:m_curReplicas), CtrlParam)

                 !*** to initialize force table on device
                 !    Note: m_ForceTable is defined in module MD_TYPEDEF_FORCELIB_GPU
                 call Init_Forcetable_DEV(SimBox, CtrlParam, gm_ForceClass)

                 !***  to give a initial neighbore list
                 call Initialize_NeighboreList_DEV(SimBox(1:m_curReplicas), CtrlParam)
                 call Cal_NeighBoreList_DEV(SimBox(1:m_curReplicas), CtrlParam)

                 !*** if active region method to be used,  initialize the modeul
                 if(iand(CtrlParam%AR_METHOD, CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                    call Initialize_ActiveRegion_DEV(SimBox(1:m_curReplicas), CtrlParam)
                    call ActivateRegion_DEV(SimBox(1:m_curReplicas), CtrlParam)
                 end if

                 !*** to reset the boxes on devices
                 allocate(IBT(m_curReplicas), TANGENT0(m_curReplicas),ETAB(m_curReplicas))
                 
                 !*** to calculate the force and potential
                 call CalForce_ForceClass(SimBox, CtrlParam, gm_ForceClass)
                 call CopyFPFrom_Devices_to_Host()

                !*** to calculate the projection of force along Dis
                IA = 0
                do IB=1, m_curReplicas
                   TANGENT0(IB) = sum(hm_DIS(IA+1:IA+NPRT,1:3)*hm_FP(IA+1:IA+NPRT,1:3)) !/ &
                   hm_FP(IA+1:IA+NPRT,1:3) = hm_FP(IA+1:IA+NPRT,1:3)                     &
                                           - ALPHA*TANGENT0(IB)*hm_DIS(IA+1:IA+NPRT,1:3) &
                                           /(sum(hm_DIS(IA+1:IA+NPRT,1:3) *hm_DIS(IA+1:IA+NPRT,1:3) ))
                   hm_XP1(IA+1:IA+NPRT,1:3) = 0.D0
                   IA = IA + NPRT
                end do
                call CopyFPFrom_Host_to_Devices()
                call CopyXP1From_Host_to_Devices()

                !*** start loop for time sections
                ITIMEP = ITIME-CtrlParam%IT0+1
                CtrlParam%DAMPTIME0 = CtrlParam%IT0
                CtrlParam%DAMPTIME1 = CtrlParam%IT1-CtrlParam%IT0 
                do while(ITIME .le. CtrlParam%IT1)
                    !--- Give a prediction
                    call Predictor_DEV(ITIME, SimBox(1), CtrlParam)

                    !--- update the neighbore list
                    if(mod(ITIMEP,CtrlParam%NB_UPTAB) .eq. C_IZERO) then
                       call Cal_NeighBoreList_DEV(SimBox(1:m_curReplicas), CtrlParam)
                    end if

                    !--- update the the active region
                    if(mod(ITIMEP, CtrlParam%AR_UPTAB) .eq. C_IZERO) then
                       call ActivateRegion_DEV(SimBox(1:m_curReplicas), CtrlParam)
                    end if

                    !--- calculate the current force
                    call CalForce_ForceClass(SimBox, CtrlParam, gm_ForceClass)
                    !*** to calculate the projection of force along Dis
                    call CopyFPFrom_Devices_to_Host()
                    call CopyDISFrom_Devices_to_Host()
                    IA      = 0
                    IBT     = 0
                    CHECKT  = 0
                    do IB=1, m_curReplicas
                       TANGENT0(IB) = sum(hm_DIS(IA+1:IA+NPRT,1:3)*hm_FP(IA+1:IA+NPRT,1:3)) !/ &
                       if(TANGENT0(IB) .ge. 0) then
                          IBT(IB) = IB
                          CHECKT  = 1
                       end if
                       IA = IA + NPRT
                    end do

                    if(CHECKT .gt. 0) then !--- possible transition found
                       call CopyOut_SimBox_DEV(SimBox(1:m_curReplicas))
                       do IB = 1, m_curReplicas
                          ETAB(IB) = sum(SimBox(IB)%EPOT-m_EPOT0)
                          SwapSimBox(IB) = SimBox(IB)
                       end do

                       call Do_Damp(SwapSimBox(1:m_curReplicas), m_CtrlParamDamp, REINIT=0)
                       call CopyOut_SimBox_DEV(SwapSimBox(1:m_curReplicas))
                       do IB = 1, m_curReplicas               
                          if(IBT(IB) .le. 0) cycle

                          call Do_Compare(m_SimBoxIni(1), SwapSimBox(IB), CtrlParam, SAME)    
                          if(SAME .le. 0) then
                             call Copy_SimMDBox(SwapSimBox(IB), m_SimBoxIni(1))
                             m_ITEvent      = ITime
                             STAMP%ITest    = ITest
                             STAMP%IBox     = ITest
                             STAMP%ITime    = m_ITEvent
                             STAMP%Time     = 0
                             !--- output intermdediate state
                             m_NumEvent     = m_NumEvent + 1
                             STAMP%ICfg     = m_NumEvent
                             STAMP%IRec     = STAMP%ICfg
                             call STRCATI(GFILE, CtrlParam%f_geometry, "_0K_P", processid, 4)
                             call STRCATI(GFILE, GFILE, "_", ITest, 4)
                             call Putout_Instance_Config_SimBoxArray(GFILE, SimBox(IB:IB), STAMP)
                             !--- output target state
                             m_NumEvent     = m_NumEvent + 1
                             STAMP%ICfg     = m_NumEvent
                             STAMP%IRec     = STAMP%ICfg
                             call STRCATI(GFILE, CtrlParam%f_geometry, "_0K_P", processid, 4)
                             call STRCATI(GFILE, GFILE, "_", ITest, 4)
                             call Putout_Instance_Config_SimBoxArray(GFILE, m_SimBoxIni, STAMP)
                             write(*,fmt="(A, I, A, I5, A)")' MDPSCU Message: Event #', (m_NumEvent-1)/2+1, ' found'
                             write(*,fmt="(A, I, A, 1PE14.6)")'                 at ITIME, ', ITIME, ', EPOT(eV) ', ETAB(IB)*CP_ERGEV 
                             exit
                          else
                             call Do_Activation(CtrlParam, m_SimBoxIni(1), SimBox(IB:IB), DSTEP)
                          end if
                       end do
                       !---- transition found, exit for next searching
                       if(SAME .eq. 0) then 
                          ITIME = ITIME + 1
                          ITIMEP = ITIMEP + 1
                          exit
                       end if
                       !--- no transition found, we continue climbing.
                       !    because the forces on device have been changed
                       !    on doing damping, we need recalculated them
                       call Initialize_Globle_Variables_DEV(SimBox(1:m_curReplicas), CtrlParam)
                       call Init_Forcetable_DEV(SimBox(1:m_curReplicas), CtrlParam, gm_ForceClass)
                       call Initialize_NeighboreList_DEV(SimBox(1:m_curReplicas), CtrlParam)
                       call Cal_NeighBoreList_DEV(SimBox(1:m_curReplicas), CtrlParam)
                       if(mod(ITIMEP, CtrlParam%AR_UPTAB) .eq. C_IZERO) then
                          call ActivateRegion_DEV(SimBox(1:m_curReplicas), CtrlParam)
                       end if
                       call CalForce_ForceClass(SimBox, CtrlParam, gm_ForceClass)
                      !--- to calculate the projection of force along Dis
                       call CopyFPFrom_Devices_to_Host()
                       call CopyDISFrom_Devices_to_Host()
                    end if

                    !---- continue for no transition occur
                    !---- do force projection
                    IA = 0
                    do IB=1, m_curReplicas
                       hm_FP(IA+1:IA+NPRT,1:3) = hm_FP(IA+1:IA+NPRT,1:3)                        &
                                               - ALPHA*TANGENT0(IB)* hm_DIS(IA+1:IA+NPRT,1:3) &
                                                 /(sum(hm_DIS(IA+1:IA+NPRT,1:3)*hm_DIS(IA+1:IA+NPRT,1:3)))
                       IA = IA + NPRT
                    end do
                    call CopyFPFrom_Host_to_Devices()

                    !--- check if the images have changes in one time step
                    TH     = CtrlParam%H
                    H2S2   = TH*TH*C_HALF
                    call CheckTimestep_DEV(ITIME,SimBox(1),CtrlParam, TH, H2S2, MXD2, CHECKT)
                    if(CHECKT .le. 0) then !--- the images are damped, to dispalce th images
                       write(*,fmt="(A, I5, A, I)")' MDPSCU Message: box images are displaced at', ITIME
                       call UpdateEPOT_ForceClass(SimBox, CtrlParam, gm_ForceClass)
                       call CopyEPOTFrom_Devices_to_Host()
                       call CopyXPFrom_Devices_to_Host()
                       call CopyDISFrom_Devices_to_Host()
                       !--- nomrlize the displacement vector to DSTEP
                       IA = 0
                       do IB = 1, m_curReplicas
                          !--- record the current potential
                           ETAB(IB) = sum(hm_EPOT(IA+1:IA+NPRT) -m_EPOT0(1:NPRT))
                           write(*,fmt="(A, I5, 1PE14.6)")'                 box images energy(eV)', IB, ETAB(IB)*CP_ERGEV
                           
                          !--- move the image to new position
                          do I=IA+1, IA+NPRT
                             hm_XP(I,1:3) = hm_XP(I,1:3) + hm_DIS(I,1:3)
                             if(IFPD(1) .gt. 0) then
                                if(hm_XP(I,1) .gt. UBOX(1)) then
                                   hm_XP(I,1) = hm_XP(I,1) - BOX(1)
                                else if(hm_XP(I,1) .lt. LBOX(1)) then
                                   hm_XP(I,1) = hm_XP(I,1) + BOX(1)
                                end if
                             end if 
                             if(IFPD(2) .gt. 0) then
                                if(hm_XP(I,2) .gt. UBOX(2)) then
                                   hm_XP(I,2) = hm_XP(I,2) - BOX(2)
                                else if(hm_XP(I,2) .lt. LBOX(2)) then
                                   hm_XP(I,2) = hm_XP(I,2) + BOX(2)
                                end if
                             end if 
                             if(IFPD(3) .gt. 0) then
                                if(hm_XP(I,3) .gt. UBOX(3)) then
                                   hm_XP(I,3) = hm_XP(I,3) - BOX(3)
                                else if(hm_XP(I,3) .lt. LBOX(3)) then 
                                   hm_XP(I,3) = hm_XP(I,3) + BOX(3)
                             end if
                             end if 
                          end do
                          hm_XP1(IA+1:IA+NPRT,1:3) = 0.D0
                          hm_DIS(IA+1:IA+NPRT,1:3) = DSTEP*hm_DIS(IA+1:IA+NPRT,1:3)/ &
                                      dsqrt(sum(hm_DIS(IA+1:IA+NPRT,1:3)*hm_DIS(IA+1:IA+NPRT,1:3))) 

                          IA = IA + NPRT
                       end do
                       
                       !--- the images are displaced, we need to update
                       !    the position and force on devices
                       call CopyXPFrom_Host_to_Devices()
                       call CopyDISFrom_Host_to_Devices()
                       call CopyXP1From_Host_to_Devices()
                       call Cal_NeighBoreList_DEV(SimBox(1:m_curReplicas), CtrlParam)
                       call CalForce_ForceClass(SimBox, CtrlParam, gm_ForceClass)
                       !--- to calculate the projection of force along Dis
                       call CopyFPFrom_Devices_to_Host()
                       call CopyDISFrom_Devices_to_Host()
                       IA = 0
                       do IB=1, m_curReplicas
                          TANGENT0(IB) = sum(hm_DIS(IA+1:IA+NPRT,1:3)*hm_FP(IA+1:IA+NPRT,1:3)) !/ &
                          IA = IA + NPRT
                       end do
                    end if  !--- end if images displaced

                    call Correction_DEV(ITIME,SimBox(1), CtrlParam)
                    ITIME = ITIME + 1
                    ITIMEP = ITIMEP + 1
                end do   !end loop for time step

                !--- deallocate the memory allocated in this routine
                call Release_SimBoxArray(SwapSimBox)
                deallocate(IBT, TANGENT0, ETAB, SwapSimBox)
                !--- restore the original number of replicas in control parameter
                call SetMultiBox_SimMDCtrl(CtrlParam, m_mxReplicas)
                

      return
  end subroutine For_One_Trail_Dyn
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Damp(SimBox, CtrlParam, Reinit)
  !***  PORPOSE: to damping the configuration to zero temperature.
  !
  !    INPUT:  CtrlParam,  the control parameters for simulation
  !            REINIT,     indicating if re-initialization of device memory is needed
  !    OUTPUT:  SimBox,    the box that has been damped
  !
  !
   use MD_NeighborsList_GPU
   use MD_ActiveRegion_GPU
   use MD_SimboxArray_GPU
 
   use MD_ForceLib_Factory_GPU
   use MD_LBFGSScheme_GPU,    only:Do_LBFGSB_Forsteps_DEV
   use MD_SteepestScheme_GPU, only:Do_Steepest_Forsteps_DEV
   use MD_CGScheme_GPU,       only:Do_CG_Forsteps_DEV
   use MD_DiffScheme_GPU,     only:Do_DynDamp_Forsteps_DEV

   implicit none
   !---dummy vaiables
       type(SimMDBox), dimension(:)           ::SimBox
       type(SimMDCtrl)                        ::CtrlParam
       integer,                     intent(in)::Reinit
      !Local variables
       integer::I, MULTIBOX, NB, IFLAG

          NB       = size(SimBox)
          MULTIBOX = CtrlParam%MULTIBOX
          CtrlParam%MULTIBOX = NB
          if(Reinit .gt. 0) then

             !*** to initialize the GPU variable
             call Initialize_Globle_Variables_DEV(SimBox(1:NB), CtrlParam)

             !*** to initialize force table on device
             call Init_Forcetable_Dev(SimBox(1:NB), CtrlParam, gm_ForceClass)

             !*** to initialize the neigborlist module
             call Initialize_NeighboreList_DEV(SimBox(1:NB), CtrlParam)
             call Cal_NeighBoreList_DEV(SimBox(1:NB), CtrlParam)

             if(iand(CtrlParam%AR_METHOD,CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                call Initialize_ActiveRegion_DEV(SimBox(1:NB), CtrlParam)
                call ActivateRegion_DEV(SimBox(1:NB), CtrlParam)
              end if
              !*** to calculate the force
              !    The force to be calculated in call in DO_LBFGS_FORSTEPS_DEV
              !    It in not necessary to calculate force here
              call  CalForce_ForceClass(SimBox, CtrlParam, gm_ForceClass)
           end if
           

           select case(iand(CtrlParam%Quench_Meth, CP_LWORD))
                  case( CP_DAMPSCHEME_LBFGS)
                     call DO_LBFGSB_FORSTEPS_DEV  (SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps, IFLAG)
                  case(CP_DAMPSCHEME_DYN)
                     call Do_DynDamp_Forsteps_DEV (SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps)
                  case(CP_DAMPSCHEME_ST )   
                     call Do_Steepest_Forsteps_DEV(SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
                  case(CP_DAMPSCHEME_CG )   
                     call Do_CG_Forsteps_DEV      (SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
            end select

            !--- copy back the dampped configuration for analysis
            call CalEpot_ForceClass(SimBox,  CtrlParam, gm_ForceClass)
            call CopyOut_SimBox_DEV(SimBox)
 
            !--- restore the originla control parameters
            CtrlParam%MULTIBOX  = MULTIBOX
 
  end subroutine Do_Damp
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Compare(SimBoxIni, SimBox, CtrlParam, Same, Mask)
  !***  PORPOSE: to check if the SimBox is different from SimBoxIni
  !
  !    INPUT:  SimBoxIni,    the 'standard' box to which the replicas to be compared
  !            SimBox,       the replicas
  !            CtrlParam,    the control parameters for simulation
  !            MASK,         indicating which atoms to be included in the comparing
  !
  !    OUTPUT: Same,         indicating if the  two boxes are the same
  !
  !
  implicit NONE
   !---dummy vaiables
       type(SimMDBox),                intent(in) ::SimBoxIni
       type(SimMDBox),                intent(in) ::SimBox
       type(SimMDCtrl),               intent(in) ::CtrlParam
       integer                                   ::Same
       integer, dimension(:),optional,intent(in) ::Mask
       !Local variables
       integer::NPRT, I, IFPD(3)
       real(KINDDF)::RC2, SEP(3), BOX(3), HBOX(3)

               RC2  =  CtrlParam%STRCUT_DRTol*CtrlParam%STRCUT_DRTol
               NPRT = SimBoxIni%NPRT
               BOX  = SimBoxIni%ZL
               HBOX = SimBoxIni%ZL*C_HALF
               IFPD = CtrlParam%IFPD
               Same = 1
               if(present(Mask)) then
                  do I=1, NPRT
                     if(Mask(I) .le. 0 ) cycle

                     SEP(1:3)  =  SimBoxIni%XP(I, 1:3) - SimBox%XP(I,1:3)
                     if( (IFPD(1).gt.0) .and. (dabs(SEP(1)) .gt. HBOX(1))) then
                         SEP(1) = SEP(1) - dsign(BOX(1),SEP(1))
                     end if

                     if( (IFPD(2).gt.0) .and. (dabs(SEP(2)) .gt. HBOX(2))) then
                         SEP(2) = SEP(2) - dsign(BOX(2),SEP(2))
                     end if

                     if( (IFPD(3).gt.0) .AND. (dabs(SEP(3)) .gt. HBOX(3))) then
                         SEP(3) = SEP(3) - dsign(BOX(3),SEP(3))
                     end if

                     if(sum(SEP*SEP) .gt. RC2) then
                        Same = 0
                        exit
                     end if
                  end do
               else
                  do I=1, NPRT
                     SEP(1:3)  =  SimBoxIni%XP(I, 1:3) - SimBox%XP(I,1:3)
                     if( (IFPD(1).gt.0) .and. (dabs(SEP(1)) .gt. HBOX(1))) then
                         SEP(1) = SEP(1) - dsign(BOX(1),SEP(1))
                     end if

                     if( (IFPD(2).gt.0) .and. (dabs(SEP(2)) .gt. HBOX(2))) then
                         SEP(2) = SEP(2) - dsign(BOX(2),SEP(2))
                     end if

                     if( (IFPD(3).gt.0) .AND. (dabs(SEP(3)) .gt. HBOX(3))) then
                         SEP(3) = SEP(3) - dsign(BOX(3),SEP(3))
                     end if

                     if(sum(SEP*SEP) .gt. RC2) then
                        Same = 0
                        exit
                     end if
                  end do
               end if
                                
               return
  end subroutine Do_Compare
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Activation(CtrlParam, SimBox0, SimBox, Dstep, Replicas)
  !***  PORPOSE: to activate replicas
  !     INPUT:  
  !             CtrlParam,          the control parameters
  !             SimBox0,            the initial configuration
  !             SimBox,             the simulation box array
  !
  use RAND32_MODULE
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                      ::CtrlParam
   type(SimMDBox)                       ::SimBox0
   type(SimMDBox), dimension(:)         ::SimBox
   real(KINDDF)                         ::Dstep
   integer,                     optional::Replicas
  !--- Local variables
       integer::NB, IB, NPRT, I, IA, IFPD(3), NC
       integer, dimension(:), allocatable::IC
       real(KINDDF)::BOX(3), HBOX(3), LBOX(3), UBOX(3), VECT(3), DSTEPP
       !------
                 NPRT = SimBox0%NPRT
                 BOX  = SimBox0%ZL
                 HBOX = 0.5D0*BOX
                 LBOX = SimBox0%BOXLOW
                 UBOX = SimBox0%BOXUP
                 IFPD = CtrlParam%IFPD
                 NB   = size(SimBox)

                 if(CtrlParam%ART_ActMeth .eq. CP_USERDEFACT_ART) then
                    if(.not.associated(m_pSelector)) then
                       write(*,fmt="(A)")  " MDPSCU Error: in ART method. user-defined procedure for"
                       write(*,fmt="(A)")  "                selecting activation candidates is null"
                       write(*,fmt="(A)")  " Process to be stopped"
                       stop
                    end if
                    call m_pSelector(CtrlParam, SimBox0, NC, IC)
                 else
                    call SelectCandidates(CtrlParam, SimBox0, NC, IC)
                 end if 
                 DSTEPP= Dstep/dble(NC)

                 do IB=1, NB
                    call Copy_SimMDBox(SimBox0, SimBox(IB))

                    SimBox(IB)%DIS = 0.D0
                    SimBox(IB)%XP1 = 0.D0
                    do I=1, NC
                       IA = IC(I)
                       if(IA .le. 0) exit 

                        VECT(1) = (DRAND32()-0.5D0)
                        VECT(2) = (DRAND32()-0.5D0)
                        VECT(3) = (DRAND32()-0.5D0)
                        SimBox(IB)%DIS(IA,1:3) = VECT(1:3)*DSTEPP/dsqrt(sum(VECT*VECT))
                        SimBox(IB)%XP(IA,1:3) =  SimBox(IB)%XP(IA,1:3) + SimBox(IB)%DIS(IA,1:3)
                        if(IFPD(1) .gt. 0) then
                           if(SimBox(IB)%XP(IA,1) .gt. UBOX(1)) then
                              SimBox(IB)%XP(IA,1) = SimBox(IB)%XP(IA,1) - BOX(1)
                           else if(SimBox(IB)%XP(IA,1) .lt. LBOX(1)) then
                              SimBox(IB)%XP(IA,1) = SimBox(IB)%XP(IA,1) + BOX(1)
                           end if
                        end if 
                        if(IFPD(2) .gt. 0) then
                           if(SimBox(IB)%XP(IA,2) .gt. UBOX(2)) then
                              SimBox(IB)%XP(IA,2) = SimBox(IB)%XP(IA,2) - BOX(2)
                           else if(SimBox(IB)%XP(IA,2) .lt. LBOX(2)) then
                              SimBox(IB)%XP(IA,2) = SimBox(IB)%XP(IA,2) + BOX(2)
                           end if
                        end if 
                        if(IFPD(3) .gt. 0) then
                           if(SimBox(IB)%XP(IA,3) .gt. UBOX(3)) then
                              SimBox(IB)%XP(IA,3) = SimBox(IB)%XP(IA,3) - BOX(3)
                           else if(SimBox(IB)%XP(IA,3) .lt. LBOX(3)) then
                              SimBox(IB)%XP(IA,3) = SimBox(IB)%XP(IA,3) + BOX(3)
                           end if
                        end if 
                    end do
                 end do
                 
                 if(present(Replicas)) then
                    Replicas = NB
                 end if
                 if(allocated(IC)) deallocate(IC)
                 return
   end subroutine Do_Activation
  !****************************************************************************************

  !****************************************************************************************
  subroutine SelectCandidates(CtrlParam, SimBox0, NC, IC)
  !***  PORPOSE: to select the candidates to be set initial displacement
  !     INPUT:  
  !             CtrlParam,          the control parameters
  !             SimBox0,            the initial configuration
  !     OUTPUT: NC,                 the number of cancdiate
  !             IC,                 the index of the candidate atoms 
  !
  use RAND32_MODULE
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                     ::CtrlParam
   type(SimMDBox)                      ::SimBox0
   integer                             ::NC
   integer, dimension(:), allocatable  ::IC
  !--- Local variables
   integer::IA, IT, NA0, NA
   integer, dimension(:), allocatable::ICFLAG
   real(kinddf)::P0, P

                 if(allocated(IC)) deallocate(IC)

                 if(iand(CtrlParam%ART_ActMeth, CP_CENTPARTACT_ART) .eq. CP_CENTPARTACT_ART) then
                    if(dabs(CtrlParam%ART_ActExt) .le. 1.D-10) then
                       call SelectCandidates0(CtrlParam, SimBox0, ICFLAG) 
                    else
                       call SelectCandidates1(CtrlParam, SimBox0, ICFLAG) 
                    end if
                 else
                 end if 

                !--- the number of candidates to be slected
                NA0 = count(ICFLAG .gt. 0)
                 if(iand(CtrlParam%ART_ActMeth, CP_MULTPLEACT_ART) .eq. CP_MULTPLEACT_ART) then   
                     NA = nint(CtrlParam%ART_ActPecent*NA0)
                 else
                     NA = 1
                 end if   

                 !--- the probability for an atom to be selected
                 P0 = dble(NA)/dble(NA0)
                 allocate(IC(NA))
                 IC = 0

                 NC = 0
                 do while(NC .lt. NA)
                    do IA=1, SimBox0%NPRT
                       if( ICFLAG(IA) .le. 0) cycle

                       P = DRAND32()
                       if(P .le. P0) then
                          NC         = NC + 1
                          IC(NC)     = IA
                          ICFLAG(IA) = 0
                          if(NC .ge. NA) exit
                       end if
                    end do
                 end do

                 if(NC .le. 0) then
                    write(*,fmt="(A)")  " MDPSCU Error: no any atom to be activated"
                    write(*,fmt="(A)")  "               check control parameters for ART"
                    write(*,fmt="(A)")  " Process to be stopped"
                    stop
                 else
                    write(*,fmt="(A, I6, A)")  " MDPSCU Message: ", NC, " atom(s) are activated"
                 end if

                 deallocate(ICFLAG)
                 return
   end subroutine SelectCandidates
  !****************************************************************************************

  !****************************************************************************************
  subroutine SelectCandidates0(CtrlParam, SimBox0, ICFLAG)
  !***  PORPOSE: to select the candidate to be set initial displacement
  !              the candicates are those selected by given type
  ! 
  !     INPUT:  CtrlParam,          the control parameters
  !             SimBox0,            the initial configuration
  !     OUTPUT: ICFLAG,     the flag indicating the candidates
  !
  use RAND32_MODULE
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                     ::CtrlParam
   type(SimMDBox)                      ::SimBox0
   integer, dimension(:), allocatable  ::ICFLAG
  !--- Local variables
   integer::IA
                    
                    !--- the total number of atoms involved in the selection
                    allocate(ICFLAG(SimBox0%NPRT))
                    ICFLAG = 0
                    do IA = 1, SimBox0%NPRT
                       if( CtrlParam%ART_ActSeed(SimBox0%ITYP(IA)) .le. 0) cycle
                       if(iand(SimBox0%STATU(IA), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) cycle
                       ICFLAG(IA) = 1
                    end do

                 return
   end subroutine SelectCandidates0
  !****************************************************************************************

  !****************************************************************************************
  subroutine SelectCandidates1(CtrlParam, SimBox0, ICFLAG)
  ! *** PORPOSE: to select the candidate to be set initial displacement
  !              the candicates are those selected by given type
  !
  !     INPUT:  CtrlParam,          the control parameters
  !             SimBox0,            the initial configuration
  !     OUTPUT: ICFLAG,     the flag indicating the candidates
  !
  use MD_SimboxArray_GPU
  use MD_NeighborsList_GPU
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                     ::CtrlParam
   type(SimMDBox)                      ::SimBox0
   integer, dimension(:), allocatable  ::ICFLAG
  !--- Local variables
   integer::IA, IN
   type(NEIGHBOR_LIST), allocatable::tList
   type(SimMDCtrl),  allocatable::tCtrlParam

                 !*** to create a neighbore list   
                 allocate(tList, tCtrlParam)                 
                 tCtrlParam       = CtrlParam
                 tCtrlParam%NB_RM = dabs(CtrlParam%ART_ActExt)*SimBox0%RR
                 call Initialize_Globle_Variables_DEV(SimBox0, tCtrlParam)
                 call Initialize_NeighboreList_DEV(SimBox0, tCtrlParam)
                 call Cal_NeighBoreList_DEV(SimBox0, tCtrlParam)
                 call Copyout_NeighboreList_DEV(tList, ORDER=1)
                    !--- the total number of atoms involved in the selection
                    allocate(ICFLAG(SimBox0%NPRT))
                    ICFLAG = 0
                    if(CtrlParam%ART_ActExt .gt. 0) then
                       do IA = 1, SimBox0%NPRT
                          if( CtrlParam%ART_ActSeed(SimBox0%ITYP(IA)) .le. 0) cycle
                          if(iand(SimBox0%STATU(IA), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) cycle
                          ICFLAG(IA) = 1
                          do IN = 1, tList%KVOIS(IA)
                             ICFLAG(tList%INDI(IA, IN)) = 1 
                          end do
                       end do
                     else
                     	! --- to active the neighbores of the seeds, by the seeds are deactivated 
                       do IA = 1, SimBox0%NPRT
                          if( CtrlParam%ART_ActSeed(SimBox0%ITYP(IA)) .le. 0) cycle
                          if(iand(SimBox0%STATU(IA), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) cycle
                          ICFLAG(IA) = 0
                          do IN = 1, tList%KVOIS(IA)
                          	! --- the neighboresof the active seed cannot be active
                           	 if( CtrlParam%ART_ActSeed(SimBox0%ITYP(tList%INDI(IA, IN))) .gt. 0) cycle
                             ICFLAG(tList%INDI(IA, IN)) = 1 
                          end do
                       end do
                     end if
                 call Clear_NeighboreList_DEV()
                 call Clear_NeighboreList(tList)
                 deallocate(tList, tCtrlParam)
                 return
   end subroutine SelectCandidates1
  !****************************************************************************************


  end module MD_Method_ART_GPU
