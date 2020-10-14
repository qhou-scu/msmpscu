  module MD_Method_GenericMD_GPU
  !***  DESCRIPTION:
  !     This module provides routine for generic MD evolution of a system.
  !
  !     HISTORY: Modified from MD_SimboxArray_AppShell_14_GPU.F90  by HOU Qing, august, 2014

  !***  The modules included ******************************************
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_RecordList

  implicit none
  contains

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
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU
  use MD_DiffScheme_GPU,        only:CalEKin_DEV
  use MD_LocalTempMethod_GPU,   only:Do_ResetParam_DEV
  use MD_TYPEDEF_PrintList,     only:Do_RestoreProcess

  implicit none
  !--- dummy variables
   type(SimMDBox),                          intent(in)::SimBox0
   type(SimMDCtrl),target                             ::CtrlParam
   type(SimMDBox), dimension(:), allocatable          ::SimBox
   type(RecordProcedureList)                          ::Recordlist
   integer,                                 intent(in)::J,processid

   !--- interface for the external process--------------
   OPTIONAL::                                     INICONFIGPROC
   external::                                     INICONFIGPROC

   interface
       SUBROUTINE INICONFIGPROC(SimBox, CtrlParam, RESTART)
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer::RESTART
       END SUBROUTINE INICONFIGPROC
   end interface
  !------END INTERFACE ---

  !--- Local variables
  integer::ITIME0, ITIME, TESTLOOP0, TESTLOOP, ITEST, JOB, ISECT0, ISECT, TSTATU
  logical::EXIST, HASCOPYOUT
  real(KINDDF)::TIME, TIME0

  character*256::SFILE, CFile, CFile0
  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2
  type(SimMDCtrl), pointer::sectCtrlParam, nextCtrlParam
  type(MDRecordStamp)::STAMP, STAMP0

  !----
             STAMP%AppTYpe = "GMD"

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
                    call STRCATI(CFILE, CtrlParam%f_configure, "P", processid, 4)
                    call STRCATI(CFILE, CFILE, "_", JOB, 4)
                    inquire(FILE=CFILE, exist=EXIST)
                 else
                    !--- for dependent box, find the last one we have calculated
                    do ITEST=TESTLOOP, 1, -1
                       call STRCATI(CFILE, CtrlParam%f_configure, "P", processid, 4)
                       call STRCATI(CFILE, CFILE, "_", ITEST, 4)
                       inquire(FILE=CFILE, exist=EXIST)
                       if(EXIST) then
                          !--- to restart from this box
                          TESTLOOP0 = ITEST
                          JOB = (ITEST-1)+J
                          exit
                       end if
                    end do
                 end if
              end if

              !*** Now we can create the ARRAY
              call Create_SimBoxArray(SimBox, CtrlParam%MULTIBOX,SimBox0)

              !*** to load initial configuration
              ITIME0 = 0
              TIME0  = C_ZERO

              if(CtrlParam%RESTART .EQ. 0 .or. .NOT.EXIST) then
                 !--- from the very begin
                 if(Simbox0%IniCfgID .eq. C_IZERO) then
                  !--- load initial file from a single file
                    write(*,fmt="(A)")      ' MDPSCU Message: Load initial substrate from single file:'
                    call Initialize_Config_SimBoxArray(SimBox, fname=SimBox0%IniConfig, fmt=SimBox0%IniCfgFmt)

                 else if(Simbox0%IniCfgID .lt. C_IZERO) then
                    !--- from the last simulations
                    call STRCATI(CFILE0, SimBox0%IniConfig, "P", processid, 4)
                    call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                    write(*,fmt="(A)")      ' MDPSCU Message: Load initial substrate from unformated file:'
                    write(*,fmt="(A, A)")   '                 ', CFILE0(1:len_trim(CFILE0))
                    call Restore_Config_SimBoxArray(SimBox, STAMP0, CFile0, Do_RestoreProcess)
                 else if(Simbox0%IniCfgID .gt. C_IZERO) then
                  !--- from immediate point of the last simulation
                   call STRCATI(CFile0, SimBox0%IniConfig, "P", processid, 4)
                   if(ibits(SimBox0%IniCfgFmt, CP_INPUT_SBOXBIT, 1).eq. 0) then
                      call STRCATI(CFile0, CFile0, "_", JOB, 4)
                      call STRCATI(CFile0, CFile0, ".",  Simbox0%IniCfgID, 4)
                      write(*,fmt="(A)") ' MDPSCU Message: Load initial configuration from multi-file:'
                      call Initialize_Config_SimBoxArray(SimBox,fname=CFile0, fmt=SimBox0%IniCfgFmt, multbox=1)
                   else    
                      call Initialize_Config_SimBoxArray(SimBox,fname=CFile0, JobID=JOB,CfgID=Simbox0%IniCfgID, fmt=SimBox0%IniCfgFmt)
                   end if   
                 end if
                 !--- initialize the configuration by user suppied routine
                 if(present(INICONFIGPROC))then
                    call INICONFIGPROC(SimBox, CtrlParam, RESTART=0)
                    !*** in case user have change the simulation box,we have to check the completeness of the box
                    call CheckForceTableAssigned_SimMDBox(SimBox)

                 endif

              !*** if we restart from a previous, restore from previous stop point
              else
                 if(present(INICONFIGPROC))then
                    if(CtrlParam%INDEPBOX) then
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=1)
                    else
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=JOB)
                    end if
                    !*** in case user have change the simulation box,we have to check the completeness of the box
                    call CheckForceTableAssigned_SimMDBox(SimBox)
                 end if

                 print *, "Restore configure from ", CFILE(1:len_trim(CFILE))
                 print *, "........ "
                 call Restore_Config_SimBoxArray(SimBox, STAMP0, CFile, Do_RestoreProcess)
                 ITIME0 = STAMP0%ITime
                 TIME0  = STAMP0%Time
                 ISECT0 = STAMP0%ISect
                 !*** to check if the restore data in consistent with present input
                 call GetSectID_SimMDCtrl(CtrlParam, ITIME0, ISECT)
                 if(ISECT0 .ne. ISECT) then
                    write(*,*) "MDPSCU Error: the present control parameter is not consistent with restored data "
                    write(*,fmt="('               check timesteps in sections', I3, ' vs ', I3)") ISECT, ISECT0
                    stop 'The process stop'
                 end if
                 write(*,fmt="(A, I8, A, I3)") "MDPSCU Message: to restart from ITIME = ", ITIME0, " in section ", ISECT0

                 call DoRestartRecord_List(Recordlist, STAMP0, SimBox, CtrlParam)
              end if

              !*** start test loop
              do ITEST =TESTLOOP0, TESTLOOP

                 JOB  = (ITEST-1)+J
                 ITIME= ITIME0
                 TIME = TIME0

                 !*** to determine the start time section firstly
                 !call GetSectID_SimMDCtrl(CtrlParam, ITIME, ISECT)
                 call GetNextByITime_SimMDCtrl(CtrlParam, ITIME, sectCtrlParam)

                !*** to initialize the GPU variable
                call Initialize_Globle_Variables_DEV(SimBox, sectCtrlParam)

                !*** to initialize force table on device
                !    Note: m_ForceTable is defined in module MD_TYPEDEF_FORCELIB_GPU
                call Init_Forcetable_Dev(SimBox, sectCtrlParam, gm_ForceClass)

                !***  to give a initial neighbore list
                call Initialize_NeighboreList_DEV(SimBox, sectCtrlParam)
                call Cal_NeighBoreList_DEV(SimBox, sectCtrlParam)

                !*** if active region method to be used,  initialize the module
                if(iand(sectCtrlParam%AR_METHOD, CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                    call Initialize_ActiveRegion_DEV(SimBox, sectCtrlParam)
                    call ActivateRegion_DEV(SimBox, sectCtrlParam)
                end if

                !*** if reset EPC coupling control parameters
                call Do_ResetParam_DEV(SimBox(1), sectCtrlParam)

                !*** to calculate the force and potential
                 call CalForce_ForceClass(SimBox, sectCtrlParam, gm_ForceClass)
                 call CalEpot_ForceClass(SimBox,  sectCtrlParam, gm_ForceClass)

                !--- prepare output of the initial thermal statu
                call CalEkin_DEV(SimBox(1), sectCtrlParam)
                call CopyOut_SimBox_DEV( SimBox)
                call Cal_thermal_quantities_SimBoxArray(SimBox)

                !--- to prepare the filename for thermal quantities
                call STRCATI(SFILE, CtrlParam%f_quantity, "P", processid, 4)
                call STRCATI(SFILE, SFILE, "_", JOB, 4)

                if(ITIME0 .gt. 0) then
                 !  call Putout_Instance_Config_SimBoxArray(ITIME, ISECT, ICFG, TIME0, GFILE, SimBox)
                 !  call Putout_ThermalQ_SimBoxArray(ITIME, TIME, ISECT, SFILE, SimBox,1)
                else
                   STAMP%ITest    = JOB
                   STAMP%NBox     = size(SimBox)
                   STAMP%IBox(1)  = (JOB-1)*size(SimBox)+1
                   STAMP%IBox(2)  = STAMP%IBox(1)+size(SimBox)-1
                   STAMP%ITime    = ITIME0
                   STAMP%Time     = TIME0
                   STAMP%ScalTime = TIME0
                   STAMP%ISect    = 0
                   STAMP%ICfg     = 0
                   STAMP%IRec     = 0
                   call Putout_Instance_Config_SimBoxArray(sectCtrlParam, SimBox, STAMP)
                   call Putout_ThermalQ_SimBoxArray(ITIME0, TIME0, 0, SFILE, SimBox,1)
                   !--- do external processing, firstly using ITIME = -1 could be useful
                   !    for external processing to do some initialization
                   call DoRecord_List(Recordlist, STAMP, SimBox, sectCtrlParam)
                end if

                !*** to recorde the start time
                 call CPU_TIME(C1)
                 TSTATU = CP_TIMECYCLE_CONTINUE
                 !*** start loop for time sections
                 do while(.TRUE.)
                    call For_One_TimeSect(CtrlParam, SimBox, Recordlist, JOB, ITIME, TIME, processid, TSTATU)
                    if(IAND(TSTATU,CP_TIMECYCLE_END) .eq. CP_TIMECYCLE_END ) exit

                    call GetNext_SimMDCtrl(sectCtrlParam, 1, nextCtrlParam)
                    if(.not. associated(nextCtrlParam)) exit
                    sectCtrlParam=>nextCtrlParam

                    !*** if active region method to be used,  initialize the module
                    if(iand(sectCtrlParam%AR_METHOD,CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
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
                    if(present(INICONFIGPROC)) then
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=0)
                       !*** in case user have change the simulation box,we have to check the completeness of the box
                       call CheckForceTableAssigned_SimMDBox(SimBox)
                    end if
                    ITIME0 = 0
                    TIME0  = C_ZERO
                 end if

                 call CPU_TIME(C2)
                 print *, "RUN TIME FOR ONE TEST: ",C2-C1
                 !--- close the output file
                 call Putout_ThermalQ_SimBoxArray(-1, 0.D0, 0, SFILE, SimBox)
              end do    !end the loop for itest
      return
  end subroutine For_One_Test
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_TimeSect(CtrlParam, SimBox, Recordlist, ITEST, ITIME, TIME, processid, TSTATU)
  !***  PORPOSE: to process one sample
  !     INPUT:  CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !             Recordlist          the external recoding routines
  !             ITEST,                the ID of the current ITEST
  !             TIME,               the current times step on starting the time section
  !             processid,          id of the process, to be used if MPI used
  !     OUTPUT: TSTATU,             the flag indicating the the whole time cycle should be ended
  !
  use MD_SimboxArray_GPU
  use MD_DiffScheme_GPU
  use MD_ForceLib_Factory_GPU
  use MD_TYPEDEF_PrintList,    only:Do_ArchiveProcess
  implicit none
  !--- dummy variables
   type(SimMDCtrl),target                              ::CtrlParam
   type(SimMDBox),dimension(:),allocatable             ::SimBox
   type(RecordProcedureList)                           ::Recordlist
   integer,                                intent(in)  ::ITEST
   integer                                             ::ITIME
   real(KINDDF)                                        ::TIME
   integer,                                intent(in)  ::processid
   integer,                                 intent(out)::TSTATU

  !--- Local variables
  integer::ITIMEP, ISECT, ICFG, IREC
  real(KINDDF)::INSTTEMP
  logical::HASCOPYOUT

  character*256::SFILE, CFile
  type(SimMDCtrl), pointer::sectCtrlParam
  type(MDRecordStamp)::STAMP

  !----
                 !*** to determine the start time section firstly
                 call GetSectID_SimMDCtrl(CtrlParam, ITIME+1, ISECT)
                 call GetNextByITime_SimMDCtrl(CtrlParam, ITIME+1, sectCtrlParam)
                 if(.not. associated(sectCtrlParam) ) then
                    TSTATU =  ior(TSTATU, CP_TIMECYCLE_END)
                    return
                 end if

                 STAMP%AppTYpe  = "GMD"
                 STAMP%ITest    = ITEST
                 STAMP%NBox     = size(SimBox)
                 STAMP%IBox(1)  = (ITEST-1)*size(SimBox)+1
                 STAMP%IBox(2)  = STAMP%IBox(1)+size(SimBox)-1
                 STAMP%ISect    = ISECT

                 !--- to prepare the filename for thermal quantities
                 call STRCATI(CFILE, CtrlParam%f_configure, "P", processid, 4)
                 call STRCATI(CFILE, CFILE, "_", ITEST, 4)
                 call STRCATI(SFILE, CtrlParam%f_quantity, "P", processid, 4)
                 call STRCATI(SFILE, SFILE, "_", ITEST, 4)

                 !*** to determine the start id of configure files
                 call GetCfgID_SimMDCtrl(CtrlParam, ITIME, ICFG)
                 call GetRecID_SimMDCtrl(CtrlParam, ITIME, IREC)

                   !*** start loop for time steps
                   do while(.TRUE.)
                      !--- use variable time step in scheme (I)
                      if(sectCtrlParam%IHDUP .GT. 0) then
                         sectCtrlParam%H  =  sectCtrlParam%HMI*(int((ITIME-sectCtrlParam%IT0+1)/sectCtrlParam%IHDUP)+1)
                         if( sectCtrlParam%H .gt.  sectCtrlParam%HMX)  sectCtrlParam%H =  sectCtrlParam%HMX
                      end if

                      !--- use the variable nieghbor-list update frequence
                      sectCtrlParam%NB_UPTAB = sectCtrlParam%NB_UPTABMI*(int((ITIME-sectCtrlParam%IT0+1)/sectCtrlParam%NB_DBITAB)+1)
                      if( sectCtrlParam%NB_UPTAB .gt. sectCtrlParam%NB_UPTABMX)  sectCtrlParam%NB_UPTAB = sectCtrlParam%NB_UPTABMX

                      !--- move system for one time step
                      !    NOTE: if CtrlParam%IHDUP < 0, the time step could be
                      !         changed in  this step
                      call For_One_Step(ITIME, INSTTEMP, sectCtrlParam, SimBox)

                      HASCOPYOUT = .false.
                      !--- output required information
                      ITIME  = ITIME  + 1
                      ITIMEP = ITIME  - sectCtrlParam%IT0 + 1 !ITIMEP + 1
                      TIME   = TIME   + sectCtrlParam%H*CP_S2PS
                      !--- Do something do not need copy out from GPU
                      call DoAfterOneStep_List(Recordlist, ITIME, SimBox, CtrlParam, TSTATU)

                      !--- output thermal quantities
                      if(sectCtrlParam%TIMESTPQ .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPQ).EQ.0) then
                            if(.not. HASCOPYOUT) then
                               call CalEpot_ForceClass(SimBox, sectCtrlParam, gm_ForceClass)
                               call CalEkin_DEV(SimBox(1), sectCtrlParam)
                               call CopyOut_SimBox_DEV( SimBox)
                               call Cal_thermal_quantities_SimBoxArray(SimBox)
                               HASCOPYOUT = .true.
                            end if

                            if(processid .eq. 0) then
                               write(*,*) 'TEST-TIME STEP = ',ITEST,ITIME, "in TIMESECT",  ISECT
                               call Putout_ThermalQ_SimBoxArray(ITIME,TIME, ISECT, SFILE, SimBox,1)
                            else
                               call Putout_ThermalQ_SimBoxArray(ITIME, TIME, ISECT, SFILE, SimBox)
                            end if
                         end if
                      end if

                      !--- prepair for output configurations
                      if(sectCtrlParam%TIMESTPG .gt. 0) then
                         if(mod(ITIMEP,sectCtrlParam%TIMESTPG).EQ.0) then
                            if(.not. HASCOPYOUT) then
                               call CalEpot_ForceClass(SimBox, sectCtrlParam, gm_ForceClass)
                               call CalEkin_DEV(SimBox(1), sectCtrlParam)
                               call CopyOut_SimBox_DEV(SimBox)
                               call Cal_thermal_quantities_SimBoxArray(SimBox)
                               HASCOPYOUT = .true.
                               !--- maybe we can use other process to save te data
                               !if(gm_NMPI .GE. 2) call MPI_SEND_Config_SimBoxArray(processid+1, m_Simbox, ITEST, ITIME0, TIME0)
                            end if
                            ICFG = ICFG + 1
                         end if
                      end if

                      !--- if external recording routine is provided, perform the routine
                      if(sectCtrlParam%TIMESTPR .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPR).EQ.0) then
                            if(.not. HASCOPYOUT) then
                               call CalEpot_ForceClass(SimBox, sectCtrlParam, gm_ForceClass)
                               call CalEkin_DEV(SimBox(1), sectCtrlParam)
                               call CopyOut_SimBox_DEV(SimBox)
                               call Cal_thermal_quantities_SimBoxArray(SimBox)
                               HASCOPYOUT = .true.
                               !--- maybe we can use other process to save te data
                               !if(gm_NMPI .GE. 2) call MPI_SEND_Config_SimBoxArray(processid+1, m_Simbox, JOB, ITIME0, TIME0)
                            end if
                            IREC = IREC + 1
                         end if
                      end if

                      !--- save current status
                      if(sectCtrlParam%TIMESTPSave .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPSave).EQ.0) then
                            if(.not. HASCOPYOUT) then
                               call CalEpot_ForceClass(SimBox, sectCtrlParam, gm_ForceClass)
                               call CalEkin_DEV(SimBox(1), sectCtrlParam)
                               call CopyOut_SimBox_DEV(SimBox)
                               call Cal_thermal_quantities_SimBoxArray(SimBox)
                               HASCOPYOUT = .true.
                               !--- maybe we can use other process to save te data
                               !if(gm_NMPI .GE. 2) call MPI_SEND_Config_SimBoxArray(processid+1, m_Simbox, ITEST, ITIME0, TIME0)
                            end if
                         end if
                      end if

                      !--- if external recording routine is provided, perform the routine
                      if(sectCtrlParam%TIMESTPR .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPR).EQ.0) then
                            STAMP%ITime       = ITIME
                            STAMP%Time        = TIME
                            STAMP%ScalTime    = TIME
                            STAMP%InstantTemp = INSTTEMP
                            STAMP%IRec        = IREC
                            STAMP%ICfg        = ICFG
                            call DoRecord_List(Recordlist, STAMP, SimBox, sectCtrlParam)
                            call RecordStatu_List(Recordlist, TSTATU)
                         end if
                      end if

                      !--- output configurations after do recording in case
                      !    if extended datapad is changes in recording prcocesses
                      if(sectCtrlParam%TIMESTPG .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPG).EQ.0) then
                               STAMP%ITime       = ITIME
                               STAMP%Time        = TIME
                               STAMP%ScalTime    = TIME
                               STAMP%ICfg        = ICFG
                               STAMP%IRec        = IREC
                               STAMP%InstantTemp = INSTTEMP
                               call Putout_Instance_Config_SimBoxArray(sectCtrlParam, SimBox, STAMP)
                         end if
                      end if

                      !--- save current status
                      if(sectCtrlParam%TIMESTPSave .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPSave).EQ.0) then
                            STAMP%ITime       = ITIME
                            STAMP%Time        = TIME
                            STAMP%ScalTime    = TIME
                            STAMP%InstantTemp = INSTTEMP
                            STAMP%ICfg        = ICFG
                            STAMP%IRec        = IREC
                            call Archive_Config_SimBoxArray(Simbox, STAMP, CFILE, Do_ArchiveProcess)
                            call DoSaveRecStatu_List(Recordlist, STAMP, SimBox, CtrlParam)
                         end if
                      end if

                      !--- check if the time point is reached
                      if(ITIME .GE. sectCtrlParam%IT1 .OR. TIME.GE.sectCtrlParam%TEMTM .OR. &
                         IAND(TSTATU,CP_TIMECYCLE_END) .EQ. CP_TIMECYCLE_END) then
                         exit
                      end if
                   end do   !end the loop for time step in current time section

      return
  end subroutine For_One_TimeSect
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Step(ITime, Temperature, CtrlParam, SimBox)
  !***  PORPOSE: to go on for only one step
  !     INPUT:   ITIME, the Ith time step
  !     NOTE:    if ITIME is in damping section, and damping scheme is CP_DAMPSCHEME_LBFGSL
  !              ITIME to be changed.
  !
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_LocalTempMethod_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU
  use MD_DiffScheme_GPU
  use MD_LBFGSScheme_GPU,    only:DO_LBFGSB_FORSTEPS_DEV
  use MD_SteepestScheme_GPU, only:Do_Steepest_Forsteps_DEV
  use MD_CGScheme_GPU,       only:Do_CG_Forsteps_DEV
  use RAND32_MODULE

  implicit none
  integer                    :: ITime
  real(KINDDF)               :: Temperature
  type(SimMDCtrl)            :: CtrlParam
  type(SimMDBox),dimension(:):: SimBox

  !Local variables
  integer, save::IVNUM = 0
  integer::IV, IP, ITERNUM, LBFGSFLAG = 0, ITIMEP, METH
  real(KINDDF), save::DT, CURT0, CURT, VEL(3)



     !--- if damping is required using non-dynamics scheme:
          if(ITIME+1 .ge. CtrlParam%DAMPTIME0 .and.                         &
             ITIME+1 .le. CtrlParam%DAMPTIME0 + CtrlParam%DAMPTIME1-1) then
             ITERNUM = CtrlParam%DAMPTIME1 
             METH    = CtrlParam%DAMPSCHEME 
             select case(iand(CtrlParam%DAMPSCHEME, CP_LWORD)) 
                    case(CP_DAMPSCHEME_LBFGS )   
                         LBFGSFLAG   = 0
                         call DO_LBFGSB_FORSTEPS_DEV(SimBox,CtrlParam, gm_ForceClass, ITERNUM, LBFGSFLAG)
                         ITIME = ITIME + ITERNUM - 1
                         ITIME = min( CtrlParam%IT1-1, ITIME)
                         call ResetXP1()
                         Temperature = 0.D0
                         return
                    case(CP_DAMPSCHEME_ST )   
                         call Do_Steepest_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, ITERNUM, METH)
                         ITIME = ITIME + ITERNUM - 1
                         ITIME = min( CtrlParam%IT1-1, ITIME)
                         Temperature = 0.D0
                         call ResetXP1()
                         return
                    case(CP_DAMPSCHEME_CG )   
                         call Do_CG_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, ITERNUM, METH)
                         ITIME = ITIME + ITERNUM - 1
                         ITIME = min( CtrlParam%IT1-1, ITIME)
                         call ResetXP1()
                         Temperature = 0.D0
                         return
                    !case(CP_DAMPSCHEME_DYN)     
                    !     ITERNUM = CtrlParam%DAMPTIME1
                    !     call Do_DynDamp_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, ITERNUM)
                    !     ITIME = (CtrlParam%DAMPTIME0 + CtrlParam%DAMPTIME1-1) - 1
                    !     Temperature = 0.D0
                    !     return
             end select
          end if

     !--- to thermalize if requird
          Temperature = CtrlParam%TI
          if(CtrlParam%IVTIME .gt. 0 ) then
             if(ITIME - CtrlParam%IVTIME0+1 .le. 0) then
                IVNUM = 0
               !*** get cuurent temperature
                call Cal_GlobalT_DEV(SimBox(1),CtrlParam, CURT0)
                DT   = (CtrlParam%TI - CURT0)/dble( CtrlParam%IT1 - ITIME + 1)
                DT   = DT*dble(CtrlParam%IVPAS)
             end if
             IV = ITIME - CtrlParam%IVTIME0 + 1
             if(MOD(IV,CtrlParam%IVPAS) .EQ. C_IZERO) then
                select case(CtrlParam%IVSCHEME)
                       case(CP_THERMALSCHEME_MC, CP_THERMALSCHEME_VSCAL, CP_THERMALSCHEME_PSCAL)
                            if(IVNUM .lt. CtrlParam%IVTIME) then
                              !call Thermalize_SimBoxArray(SimBox, CtrlParam%TI,CtrlParam%IVSCHEME)
                              !call CopyIn_SimBox_Vel_DEV(SimBox)
                              call Thermalizing_MC_DEV(SimBox, CtrlParam, CtrlParam%TI)
                              IVNUM =  IVNUM + 1
                            end if
                       case(CP_THERMALSCHEME_ASCAL)
                            CURT0 = CURT0 + DT
                            call VelScaling_DEV(SimBox(1),CtrlParam, CURT0)
                            Temperature = CURT0
                            !$$*** check kinetic energiese >= 0 involved
                            !$$call CalEkin_DEV(m_SimBox(1), CtrlParam)
                            !$$CURT = C_TWO*sum(m_EKIN, mask=(m_EKIN .ge. 0.D0))/dble(count(m_EKIN .ge. 0.D0))/(C_THR*CP_KB)
                            !$$print *, ITIME, DT, CURT0, CURT
                end select
             end if
          end if

     !--- Give a prediction
         call Predictor_DEV(ITIME, SimBox(1), CtrlParam)

     !--- update the neighbore list
         if(mod(ITIME-CtrlParam%IT0,CtrlParam%NB_UPTAB) .eq. C_IZERO) then
            call Cal_NeighBoreList_DEV(SimBox, CtrlParam)
         end if

     !--- update the the active region
          if(mod(ITIME-CtrlParam%IT0, CtrlParam%AR_UPTAB) .eq. C_IZERO) then
              call ActivateRegion_DEV(SimBox, CtrlParam)
          end if

     !--- calculate the current force
         if(CtrlParam%IFBPC .eq. 0) then
            call CalForce_ForceClass(SimBox, CtrlParam, gm_ForceClass)
         else
            !--- if pressure is needed, use the following subroutine
            call CalPTensor_ForceClass(SimBox, CtrlParam, gm_ForceClass)
         end if

     !---  the modification when there is E-P coupling in mode of global
         if(CtrlParam%TI_CTRL .eq. CP_TICTRL_GLOBAL) then
            if(ITIME .ge. CtrlParam%EPCTIME0 .and.   &
               ITIME .le. CtrlParam%EPCTIME0 + CtrlParam%EPCTIME1-1) then
               call Do_EPCForce_DEV(SimBox, CtrlParam)
            end if
         else
            call Do_EPCForce_DEV(SimBox, CtrlParam)
         end if

     !--- to do the correction
         call Correction_DEV(ITIME, SimBox(1), CtrlParam)

         !--- to do velocity correction if there are atoms having "free" property
         if(any(iand(SimBox(1)%Prop, CP_STATU_FREEPART).eq. CP_STATU_FREEPART )) then
            call CopyXP1From_Devices_to_Host()
            call CopyStatuFrom_Devices_to_Host()
            do IP =1, dm_NPRT
               if(iand(hm_STATU(IP), CP_STATU_FREEPART)  .eq. CP_STATU_FREEPART .and. &
                  iand(hm_STATU(IP), CP_STATU_PASSBOUND) .eq. CP_STATU_PASSBOUND) then
                  VEL(1) = DRAND32() - 0.5D0
                  VEL(2) = DRAND32() - 0.5D0
                  VEL(3) = DRAND32() - 0.5D0
                  VEL    = VEL/dsqrt(sum(VEL*VEL))
                  hm_XP1(IP,1:3) =  dsqrt(sum(hm_XP1(IP,1:3)*hm_XP1(IP,1:3)))*VEL
                  hm_STATU(IP)   =  iand(CP_STATU_LWORD, hm_STATU(IP))
               end if
            end do
            call CopyStatuFrom_Host_to_Devices()
            call CopyXP1From_Host_to_Devices()
         end if

     !--- at present the DEVICE version not support the pressure calculation
     !    for multiple-box, thus we do not use BPC coupling
     !--- the pressure coupling if constant pressure required
     !    if(CtrlParam%IFBPC .ne. 0) then
     !      call BPC_MODIFICATION_DEV(m_SimBox(1), CtrlParam)
     !    end if

     !--- to modify the statu of groups according their RSTATU
         !call Keep_Box_Zero_Translation_SimMDBox_DEV(m_Simbox)

     return
  end subroutine For_One_Step

  !****************************************************************************************

  end module MD_Method_GenericMD_GPU
