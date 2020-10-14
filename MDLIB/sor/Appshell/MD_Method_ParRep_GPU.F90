  module MD_Method_ParRep_GPU
  !***  DESCRIPTION:
  !     This module provides routines for Parallel-Replica MD evolution of a system.
  !
  !     HISTORY:   Modified from MD_SimboxArray_AppShell_14_GPU.F90  by HOU Qing, aug., 2016
  !     SEE ALSO:  MD_Method_GenericMD_GPU.F90,
  !                MD_ForceLib_Factory_GPU.F90

  !***  The modules included ******************************************
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_RecordStamp
  use MD_TYPEDEF_RecordList
  use MD_ForceClass_Register_GPU
  implicit none

      type(SimMDCtrl),                          private:: m_CtrlParamDamp          !parameters control damping operation
      type(SimMDCtrl),                          private:: m_CtrlParamDephase       !parameters control dephase
      integer,                                  private:: m_mxReplicas  = 0        !max number of replicas
      integer,                                  private:: m_curReplicas = 0        !current number of replicas
      integer,                                  private:: m_InterMXNum  = 0        !permiited number of times of stroing intermediate boxes
      integer,                                  private:: m_InterNum    = 0        !number of times of stroing intermediate boxes
      integer,       dimension(:), allocatable, private:: m_InterITime
      real(KINDDF),  dimension(:), allocatable, private:: m_InterTime
      type(SimMDBox),dimension(:), allocatable, private:: m_InterSimBox            !the intermediate boxes before change occur
      type(SimMDBox),dimension(1),              private:: m_SimBoxIni              !the initial configuration at 0K, e.g., the RECTANT
      integer,                                  private:: m_ITimeTgt
      real(KINDDF),                             private:: m_TimeTgt
      real(KINDDF),                             private:: m_ScalTimeTgt
      integer,                                  private:: m_NumTgt   = 0
      integer,                                  private:: m_NumEvent = 0           ! the number of events, it could be larger
                                                                                         ! than m_NumTgt when the changed is too fast



      !*** In different applictions, the definition of struture change could
      !    be different. User can provid user-defined routine to compare two
      !    configuration
      private:: COMPARE_PROC
      abstract interface
        subroutine COMPARE_PROC(SimBoxIni, SimBox, CtrlParam, MASK, FLAG)
        !    INPUT:   SimBoxIni,    the 'standard' box to which the replicas to be compared
        !             SimBox,       the replicas,
        !             CtrlParam,    the control parameters for simulation
        !             MASK,         integer, indicating which atom to be included.
        !                                   the dimsneion of MASK = NPRT
        !                                   MASK(i) > 0,  ithe atom to be included, otherwise not included
        !
        !    OUTPUT:  FLAG,        integer, indicating which atom position has been changes.
        !                                   the dimsneion of FLAG = NPRT*size(SimBox)
        !                                   FLAG(i) > 0,  if the ithe atom position has been changed
        !                                   FLAG(i) <= 0, if the ithe atom position has not been changed
        !
           use MD_TYPEDEF_SimMDBox
           use MD_TYPEDEF_SimMDCtrl
           implicit none
           !--- dummy vaiables
           type(SimMDBox),               intent(in)::SimBoxIni
           type(SimMDBox), dimension(:), intent(in)::SimBox
           type(SimMDCtrl),              intent(in)::CtrlParam
           integer,        dimension(:), intent(in)::MASK
           integer,        dimension(:)            ::FLAG
        end subroutine COMPARE_PROC
      end interface

      private::COMPAREMASK
      abstract interface
        subroutine COMPAREMASK(SimBox, CtrlParam, MASK)
        !    PURPOSE: to create mask which indicate which atoms to be included in
        !             structual comparison
        !
        !    INPUT:   SimBox,       the replicas,
        !             CtrlParam,    the control parameters for simulation
        !
        !    OUTPUT: MASK,         integer, indicating which atom to be included.
        !                                   the dimsneion of MASK = NPRT
        !                                   MASK(i) > 0,  ithe atom to be included, otherwise not included
        !
           use MD_TYPEDEF_SimMDBox
           use MD_TYPEDEF_SimMDCtrl
           implicit none
           !--- dummy vaiables
           type(SimMDBox),                     intent(in) ::SimBox
           type(SimMDCtrl),                    intent(in) ::CtrlParam
           integer, dimension(:), allocatable             ::MASK
        end subroutine COMPAREMASK
      end interface


      procedure(COMPARE_PROC), private, pointer::m_pCfgComp     => null()
      procedure(COMPAREMASK),  private, pointer::m_pCfgCompMask => null()
      type(MDForceClassGPU),   private         ::m_ForceClassB
      private::CopyOut_from_GPU

  contains

  !****************************************************************************************
  subroutine SetCompareRoutine(COMPARE, COMPAREMASK)
  !***  PORPOSE: to set the user provided routine for compare structures
  !
  implicit NONE
  !--- interface to the external routine -------------------
  
  
      interface
        subroutine COMPARE(SimBoxIni, SimBox, CtrlParam, MASK, FLAG)
           use MD_TYPEDEF_SimMDBox
           use MD_TYPEDEF_SimMDCtrl
           implicit none
           !--- dummy vaiables
           type(SimMDBox),               intent(in)::SimBoxIni
           type(SimMDBox), dimension(:), intent(in)::SimBox
           type(SimMDCtrl),              intent(in)::CtrlParam
           integer,        dimension(:), intent(in)::MASK
           integer,        dimension(:)            ::FLAG
        end subroutine COMPARE
      end interface

      interface
        subroutine COMPAREMASK(SimBox, CtrlParam, MASK)
           use MD_TYPEDEF_SimMDBox
           use MD_TYPEDEF_SimMDCtrl
           implicit none
           !--- dummy vaiables
           type(SimMDBox),                  intent(in) :: SimBox
           type(SimMDCtrl),                 intent(in) :: CtrlParam
           integer,dimension(:),allocatable            :: MASK
        end subroutine COMPAREMASK
      end interface

      optional::COMPARE, COMPAREMASK
      external::COMPARE, COMPAREMASK
 !--- END INTERFACE --------------------------------

            if(present(COMPARE) ) then
               m_pCfgComp =>  COMPARE
            else
               m_pCfgComp =>  null()
            end if

            if(present(COMPAREMASK) ) then
               m_pCfgCompMask =>  COMPAREMASK
            else
               m_pCfgCompMask =>  null()
            end if
            return
  end subroutine SetCompareRoutine
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
  use MD_SimboxArray
  use MD_SimboxArray_GPU
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_DiffScheme_GPU
  use MD_ForceLib_Factory_GPU
  use MD_TYPEDEF_PrintList,    only:Do_RestoreProcess
  implicit none
  !--- dummy variables
   type(SimMDBox),                            intent(in)::SimBox0
   type(SimMDCtrl),target                               ::CtrlParam
   type(SimMDBox), dimension(:), allocatable            ::SimBox
   type(RecordProcedureList)                            ::Recordlist
   integer,                                   intent(in)::J,processid

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
       integer                     ::RESTART
       END SUBROUTINE INICONFIGPROC
   end interface
  !------END INTERFACE ---

  !--- Local variables
  integer::ITIME, ITIME0,TESTLOOP0, TESTLOOP, ITEST, JOB, ISECT0, ISECT, TSTATU
  logical::EXIST, HASCOPYOUT
  real(KINDDF)::TIME, TIME0

  character*256::SFILE, GFile, CFile, CFile0
  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2
  type(SimMDCtrl), pointer::sectCtrlParam, nextCtrlParam
  type(MDRecordStamp)::STAMP, STAMP0

            !***
              STAMP%AppTYpe = "PARREP"
              m_mxReplicas  = CtrlParam%MULTIBOX
              m_curReplicas = m_mxReplicas

            !*** prepare the control parameters for damping
              call Copy_SimMDCtrl(CtrlParam, m_CtrlParamDamp)
              m_CtrlParamDamp%IVTIME        = 0
              if(m_CtrlParamDamp%Quench_Steps .le. 1) then
                  write(*,fmt="(A, I7, A)")   ' MDPSCU Warning: the number of timesteps for quenching is too small ', m_CtrlParamDamp%DAMPTIME1
                  write(*,fmt="(A, I7, A)")   '                 the default value 1000 to be set '
                  m_CtrlParamDamp%Quench_Steps = 1000
                  call ONWARNING(gm_OnWarning)
              end if

             !*** to prepare the parameters for thermalization
              call Copy_SimMDCtrl(CtrlParam, m_CtrlParamDephase)
              m_CtrlParamDephase%DAMPTIME0 = 1
              m_CtrlParamDephase%DAMPTIME1 = 0

              if(m_CtrlParamDephase%IVTIME .le. 0) then
                 write(*,fmt="(A, I7, A)")   ' MDPSCU Warning: the number of thermalization is not supplied'
                 write(*,fmt="(A, I7, A)")   '                 the default value 50 to be set '
                 call ONWARNING(gm_OnWarning)
                 m_CtrlParamDephase%IVTIME  = 50
                 m_CtrlParamDephase%IVTIME0 = 1
              end if

              if(m_CtrlParamDephase%IVPAS .le. 0) then
                 write(*,fmt="(A, I7, A)")   ' MDPSCU Warning: the timesteps for thermalization is not supplied'
                 write(*,fmt="(A, I7, A)")   '                 the default value 100 to be set '
                 call ONWARNING(gm_OnWarning)
                 m_CtrlParamDephase%IVPAS = 100
              end if

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

              !*** Now wee can create the ARRAY
              call Create_SimBoxArray(SimBox, m_mxReplicas, SimBox0)

              !*** to load initial configuration
              ITIME0        = 0
              TIME0         = 0.D0
              m_NumTgt      = 0
              m_NumEvent    = 0
              m_ITIMETgt    = 0
              m_TimeTgt     = 0.D0
              m_ScalTimeTgt = 0.D0

              if(CtrlParam%RESTART .EQ. 0 .or. .NOT.EXIST) then
                 !--- from the very begin
                 if(Simbox0%IniCfgID .eq. C_IZERO) then
                  !--- load initial file from a single file
                    write(*,fmt="(A)")      ' MDPSCU Message: Load initial substrate from single file:'
                    call Initialize_Config_SimBoxArray(SimBox, SimBox0%IniConfig, SimBox0%IniCfgFmt)

                 else if(Simbox0%IniCfgID .gt. C_IZERO) then
                  !--- from immediate point of the last simulation
                   call STRCATI(CFile0, SimBox0%IniConfig, "P", processid, 4)
                   call STRCATI(CFile0, CFile0, "_", JOB, 4)
                   call STRCATI(CFile0, CFile0, ".",  Simbox0%IniCfgID, 4)
                   write(*,fmt="(A)")      ' MDPSCU Message: Load initial substrate from multi-files'
                   !--- NOTE: for Parallel-Replica method, all the box start from the same initial configuration
                   !          ref. MD_Method_GenericMD_GPU.F90 for the difference
                    call  Initialize_Config_SimBoxArray(SimBox, CFile0, SimBox0%IniCfgFmt)
                 end if
                 !--- now initialize the configuration
                  if(present(INICONFIGPROC))then
                     call INICONFIGPROC(SimBox, CtrlParam, RESTART=0)
                     !*** in case user have change the simulation box,we have to check the completeness of the box
                     call CheckForceTableAssigned_SimMDBox(SimBox)
                  endif
                  call Copy_SimMDBox(SimBox(1), m_SimBoxIni(1))

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

                 print *, "Restore rectant configure from ", CFILE0(1:len_trim(CFILE0))
                 print *, "........ "
                 call Copy_SimMDBox(SimBox(1), m_SimBoxIni(1))
                 call Restore_Config_SimMDBox(m_SimBoxIni(1), STAMP0, CFile0, Do_RestoreProcess)
                 m_ITIMETgt    = STAMP0%ITime
                 m_TimeTgt     = STAMP0%Time
                 m_ScalTimeTgt = STAMP0%ScalTime
                 m_NumTgt      = STAMP0%ICfg(1)
                 m_NumEvent    = STAMP0%ICfg(2)
                 write(*,fmt="(A, I12, A, I3)") " MDPSCU Message: current events number =        ", m_NumTgt
                 write(*,fmt="(A, I12, A, I3)") "                 last event occur at ITIME    = ", m_ITIMETgt
                 write(*,fmt="(A, 1PE12.5)")    "                                     TIME     = ", m_TimeTgt
                 write(*,fmt="(A, 1PE12.5)")    "                                     SCALTIME = ", m_ScalTimeTgt

                 ITIME0        = m_ITIMETgt
                 TIME0         = m_TimeTgt
              end if

              !*** start test loop
              do ITEST =TESTLOOP0, TESTLOOP
                 JOB       = (ITEST-1)+J
                 ITIME     = ITIME0
                 TIME      = TIME0

                 !*** The initial configuration loaded in from files could be not in
                 !     equilibraium at zero temperature.
                 !     We need to quench the system to generate the starting confiuration
                 !     to be used checking occurance of a event.
                 if(ITIME .EQ. 0 ) then
                    write(*,fmt="(A)") ' MDPSCU Message: do initial damping'
                    call Do_Damp(m_SimBoxIni, m_CtrlParamDamp, gm_ForceClass, REINIT=1)

                   !*** save the reactant for restart
                   STAMP%ITest    = JOB
                   STAMP%IBox     = JOB
                   STAMP%ITime    = m_ITIMETgt
                   STAMP%Time     = m_TIMETgt
                   STAMP%ScalTime = m_SCALTIMETgt
                   STAMP%ISect    = 0
                   STAMP%ICfg(1)  = m_NumTgt
                   STAMP%ICfg(2)  = m_NumEvent
                   STAMP%IRec     = STAMP%ICfg
                   call STRCATI(CFILE0, CtrlParam%f_geometry, "_0K_P", processid, 4)
                   call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                   call Putout_Instance_Config_SimBoxArray(CFILE0, m_SimBoxIni, STAMP)
                   write(*,fmt="(A)") ' MDPSCU Message: to do initial dephasing'
                   call Do_DePhase(m_SimBoxIni(1), SimBox, gm_ForceClass)
                !*** for restart, we also need to do dephase
                 else if(CtrlParam%RESTART .gt. 0 .and. EXIST) then !
                   write(*,fmt="(A)") ' MDPSCU Message: to do dephasing'
                   call Do_DePhase(m_SimBoxIni(1), SimBox, gm_ForceClass)
                 end if

                 !*** the actually created replicas  could be different from MULTIBOX
                 call SetMultiBox_SimMDCtrl(CtrlParam, m_curReplicas)

                 !*** to determine the start time section firstly
                 call GetNextByITime_SimMDCtrl(CtrlParam, ITIME, sectCtrlParam)

                 !*** to initialize the GPU variable
                 call Initialize_Globle_Variables_DEV(SimBox(1:m_curReplicas), sectCtrlParam)

                 !*** to initialize force table on device
                 !    Note: m_ForceTable is defined in module MD_TYPEDEF_FORCELIB_GPU
                 !          and created on register the forcetable, here we
                 !          copy the force table on devices
                 call Init_Forcetable_Dev(SimBox, sectCtrlParam, gm_ForceClass)

                 !***  to give a initial neighbore list
                 call Initialize_NeighboreList_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                 call Cal_NeighBoreList_DEV(SimBox(1:m_curReplicas), sectCtrlParam)

                !*** if active region method to be used,  initialize the modeul
                if(iand(sectCtrlParam%AR_METHOD,CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                    call Initialize_ActiveRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                    call ActivateRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                end if

                 !*** to calculate the force and potential
                 call CalForce_ForceClass(SimBox, sectCtrlParam, m_ForceClassB)
                 call CalEpot_ForceClass(SimBox,  sectCtrlParam, m_ForceClassB)

                 !--- prepare output of the initial thermal statu
                 call CalEkin_DEV(SimBox(1), sectCtrlParam)
                 call CopyOut_SimBox_DEV(SimBox(1:m_curReplicas))
                 call Cal_thermal_quantities_SimBoxArray(SimBox(1:m_curReplicas))

                 !--- do external processing, firstly using ITIME = -1 could be useful
                 !    for external processing to do some initialization
                 STAMP%ITest    = JOB
                 STAMP%ITime    = -1
                 STAMP%Time     = -1
                 call DoRecord_List(Recordlist, STAMP, SimBox(1:m_curReplicas), sectCtrlParam)
                 STAMP%ITime    = ITIME0
                 STAMP%Time     = TIME0
                 call DoRecord_List(Recordlist, STAMP, SimBox(1:m_curReplicas), sectCtrlParam)

                 !--- to prepare the filename for thermal quantities
                 call STRCATI(SFILE, CtrlParam%f_quantity, "P", processid, 4)
                 call STRCATI(SFILE, SFILE, "_", JOB, 4)
                 call STRCATI(GFILE, CtrlParam%f_geometry, "P", processid, 4)
                 call STRCATI(GFILE, GFILE, "_", JOB, 4)
                 if(ITIME0 .le. 0) then
                    STAMP%ITest    = JOB
                    STAMP%ITime    = ITIME0
                    STAMP%Time     = TIME0
                    STAMP%ScalTime = m_SCALTIMETgt
                    STAMP%ISect    = 0
                    STAMP%ICfg     = 0
                    STAMP%IRec     = STAMP%ICfg
                    call Putout_Instance_Config_SimBoxArray(GFILE, SimBox(1:m_curReplicas), STAMP)
                    call Putout_ThermalQ_SimBoxArray(ITIME0, TIME0, 0, SFILE, SimBox(1:m_curReplicas),1)
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
                    if(iand(sectCtrlParam%AR_METHOD, CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                       call Initialize_ActiveRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                       call ActivateRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                    else
                       call Clear_ActiveRegion_DEV()
                       call Active_All_ActiveRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                    end if
                 end do   !end loop for time section

                 !*** restore the original value of  MULTIBOX
                 call SetMultiBox_SimMDCtrl(CtrlParam, m_mxReplicas)

                 !--- if the box are dependent of each other, prepair the initial
                 !    configuration of next job
                 if(.not.CtrlParam%INDEPBOX) then
                    if(present(INICONFIGPROC))&
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=0)
                    ITIME0    = 0
                    TIME0     = 0.D0
                 end if

                 call CPU_TIME(C2)
                 print *, "RUN TIME FOR ONE TEST: ",C2-C1
                 !--- close the output file
                 call Putout_ThermalQ_SimBoxArray(-1, 0.D0, 0, SFILE, SimBox(1:m_mxReplicas))
              end do    !end the loop for itest
      return
  end subroutine For_One_Test
  !****************************************************************************************

  !****************************************************************************************
  subroutine CopyOut_from_GPU(SimBox,CtrlParam)
  !***  PORPOSE: to copy the system statu from GPU
  !     INPUT:
  !
  use MD_Globle_Variables
  use MD_SimboxArray_GPU
  use MD_DiffScheme_GPU
  use MD_ForceLib_Factory_GPU

  implicit none
  !--- dummy variables
   type(SimMDCtrl),target      :: CtrlParam
   type(SimMDBox), dimension(:):: SimBox
   !--- Local variables

            call CalEpot_ForceClass(SimBox,  CtrlParam, m_ForceClassB)
            call CalEkin_DEV(SimBox(1), CtrlParam)
            call CopyOut_SimBox_DEV(SimBox)
            call Cal_thermal_quantities_SimBoxArray(SimBox)
            return
 end subroutine CopyOut_from_GPU
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_TimeSect(CtrlParam, SimBox, Recordlist, Job, ITime, Time, Processid, Tstatu)
  !***  PORPOSE: to process one sample
  !     INPUT:  CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !             Recordlist          the external recoding routines
  !             JOB,                the ID of the current ITEST
  !             ITIME,              the current times step on starting the time section
  !             TIME,               the current time
  !             processid,          id of the process, to be used if MPI used
  !
  !     OUTPUT: TSTATU,             the flag indicating the the whole time cycle should be ended
  !
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_TYPEDEF_PrintList, only:Do_ArchiveProcess
  implicit none
  !--- dummy variables
   type(SimMDCtrl),target                  ::CtrlParam
   type(SimMDBox), dimension(:)            ::SimBox
   type(RecordProcedureList)               ::Recordlist
   integer,                     intent(in) ::Job
   integer                                 ::ITime
   real(KINDDF)                            ::Time
   integer,                     intent(in) ::Processid
   integer,                     intent(out)::Tstatu
   !--- Local variables
   integer::ITIMEP, ISECT, II, NCB, NCBT, IB, IBT , ITT, IT, NPT, EXITATEVENT
   logical::HASCOPYOUT, RESETINTER

   character*256::SFILE, GFile, CFile, CFile0, GFile0
   type(SimMDCtrl), pointer::sectCtrlParam
   type(SimMDBox)          ::SimBox0K, SimBoxTK
   type(SimMDCtrl)         ::oldCtrlParam
   real(KINDDF)            ::TIME0, DT, SCALTIME
   type(MDRecordStamp)     ::STAMP

  !----
                 !*** to determine the start time section firstly
                 call GetSectID_SimMDCtrl(CtrlParam, ITIME+1, ISECT)
                 call GetNextByITime_SimMDCtrl(CtrlParam, ITIME+1, sectCtrlParam)
                 if(.not. associated(sectCtrlParam) ) then
                    TSTATU =  ior(TSTATU, CP_TIMECYCLE_END)
                    return
                 end if

                 STAMP%AppType   = gm_AppType
                 STAMP%ITest     = Job
                 STAMP%IBox      = Job
                 STAMP%ISect     = ISECT

                 !--- to prepare the filename for thermal quantities
                 call STRCATI(CFILE,  CtrlParam%f_configure, "P", processid, 4)
                 call STRCATI(CFILE,  CFILE, "_", JOB, 4)
                 call STRCATI(CFILE0, CtrlParam%f_configure, "_0K_P", processid, 4)
                 call STRCATI(CFILE0, CFILE0, "_", JOB, 4)
                 call STRCATI(SFILE,  CtrlParam%f_quantity, "P", processid, 4)
                 call STRCATI(SFILE,  SFILE, "_", JOB, 4)
                 call STRCATI(GFILE,  CtrlParam%f_geometry, "P", processid, 4)
                 call STRCATI(GFILE,  GFILE, "_", JOB, 4)
                 call STRCATI(GFILE0,  CtrlParam%f_geometry, "_0K_P", processid, 4)
                 call STRCATI(GFILE0,  GFILE0, "_", JOB, 4)

                   !***   NOTE: in parallel-replica method thermalization anymore
                   !         after dephasing
                   call Copy_SimMDCtrl(sectCtrlParam, oldCtrlParam)
                   sectCtrlParam%IVTIME    = 0
                   sectCtrlParam%DAMPTIME1 = 0

                   !***  On start a time-section, no boost force   
                   call ClearExtForce_ForceClass(m_ForceClassB)
                   call CopyPointer_ForceClass(gm_ForceClass, m_ForceClassB)
  
                   !**** to prepare space for the intermediate configuration
                   !     to be used to identify the exact transiption poinst
                   m_InterMXNum = (sectCtrlParam%PARREP_Timestep-1)/max(sectCtrlParam%PARREP_FineR,1) + 1
                   allocate(m_InterITime(m_InterMXNum), m_InterTIME(m_InterMXNum), m_InterSimBox(m_InterMXNum*m_mxReplicas) )
                   m_InterITime(1) = ITime
                   m_InterTIME(1)  = Time
                   do II=1, m_curReplicas
                      call Copy_SimMDBox(SimBox(II), m_InterSimBox(II))
                   end do
                   m_InterNum = 1
                   NPT        = 1

                   if(sectCtrlParam%ExitAtEvent .le. 0) then
                      EXITATEVENT = sectCtrlParam%IT1
                   else
                      EXITATEVENT = sectCtrlParam%ExitAtEvent
                   end if

                   TIME0      = Time
                   SCALTIME   = m_SCALTIMEtgt
                   !**** to begine the time evolution
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
                      call For_One_Step(ITIME, sectCtrlParam, SimBox(1:m_curReplicas))

                      HASCOPYOUT = .false.
                      RESETINTER = .false.

                      !--- output required information
                      ITime    = ITime  + 1
                      ITIMEP   = ITime  - sectCtrlParam%IT0 + 1
                      DT       = sectCtrlParam%H*CP_S2PS
                      SCALTIME = SCALTIME + DT*m_curReplicas
                      Time     = Time + DT

                      !--- output thermal quantities
                      if(MOD(ITIMEP,sectCtrlParam%TIMESTPQ).eq.0) then
                         if(.not. HASCOPYOUT) then
                            call CopyOut_from_GPU(SimBox(1:m_curReplicas), sectCtrlParam)
                            HASCOPYOUT = .true.
                         end if

                         if(processid .eq. 0) then
                            write(*,*) 'JOB-TIME STEP = ',Job,ITime, "in TIMESECT",  ISECT
                            call Putout_ThermalQ_SimBoxArray(ITime, Time, ISECT, SFILE, SimBox(1:m_curReplicas),1)
                         else
                            call Putout_ThermalQ_SimBoxArray(ITime, Time, ISECT, SFILE, SimBox(1:m_curReplicas))
                         end if
                      end if

                      !--- if external recording routine is provided, perform the routine
                      if(sectCtrlParam%TIMESTPR .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%TIMESTPR).EQ.0) then
                            if(.not. HASCOPYOUT) then
                               call CopyOut_from_GPU(SimBox(1:m_curReplicas), sectCtrlParam)
                               HASCOPYOUT = .true.
                            end if
                            !--- if external process provided, we do the process
                            STAMP%ITime    = ITime
                            STAMP%Time     = Time
                            STAMP%ScalTime = SCALTIME
                            call DoRecord_List(Recordlist, STAMP, SimBox(1:m_curReplicas), sectCtrlParam)
                            call RecordStatu_List(Recordlist, Tstatu)
                         end if
                      end if

                      !--- output and check transition event configurations
                      if(MOD(ITIMEP,sectCtrlParam%PARREP_Timestep).eq.0) then
                         if(.not. HASCOPYOUT) then
                             call CopyOut_from_GPU(SimBox(1:m_curReplicas), sectCtrlParam)
                             HASCOPYOUT = .true.
                         end if

                            !--- to detect if change(s) occur
                            write(*,fmt="(A, I, A)")   ' MDPSCU Message: to detect the transition of states at ITIME = ',ITime
                            call Do_ChangeDetect(m_SimBoxIni(1), SimBox(1:m_curReplicas), sectCtrlParam, IB, SimBox0K, NCB, REINITGPU=0)
                            write(*,fmt="(A, I, A)")    ' MDPSCU Message: structural change detected in', NCB, ' boxes'

                            if(NCB .le.0) then
                               write(*,fmt="(A, I, A)")   ' MDPSCU Message: no transition is detected '
                               m_InterNum      = 1
                               m_InterITime(1) = ITime
                               m_InterTIME(1)  = Time
                               do II=1, m_curReplicas
                                  call Copy_SimMDBox(SimBox(II), m_InterSimBox(II))
                               end do

                            else !--- if NCB > 0,indicating change occur in at least for one box
                               m_ITIMEtgt     = ITime
                               m_TIMEtgt      = Time
                               ITT            = m_InterNum
                               call Copy_SimMDBox(SimBox(IB), SimBoxTK)
                               !--- to refine the time at which the change occur
                               write(*,fmt="(A, I4, A)")   ' MDPSCU Message: to refine time point of  strutural change in ', m_InterNum, ' inter.configs.'
                               call RefineTranStat(m_SimBoxIni(1), sectCtrlParam, ITT, IBT, SimBox0K, NCBT)

                               if(NCBT .gt. 0) then
                                  m_ITIMEtgt     = m_InterITIME(ITT)
                                  m_TIMEtgt      = m_InterTIME(ITT)
                                  NCB            = NCBT
                                  IB             = IBT - int((IBT-1)/m_mxReplicas)*m_mxReplicas
                               end if
                               if(NCB .gt. C_IUN) then
                                  write(*,fmt="(A, I7, A)")   ' MDPSCU Warning: more than one boxes are found structurally changed '
                                  write(*,fmt="(A, I7, A)")   '                 considering to use less timesteps to for each change detecting '
                                  call ONWARNING(gm_OnWarning)
                               end if
                               !--- to reconstruct the transition trajectory
                               m_NUMTgt   = m_NUMTgt + 1
                               m_NUMEvent = m_NUMEvent + NCB

                               !--- Now we save the transition transition state
                               write(*,fmt="(A, I, A)") ' MDPSCU Message: ',m_NUMTgt, 'th jump event detected '
                               STAMP%NBox      = 1
                               STAMP%ITime     = m_ITIMETgt
                               STAMP%Time      = m_TIMETgt
                               STAMP%ScalTime  = m_SCALTIMEtgt + (m_TIMEtgt-TIME0)*m_curReplicas
                               STAMP%ICfg(1)   = m_NUMTgt
                               STAMP%ICfg(2)   = m_NUMEvent
                               STAMP%IRec      = STAMP%ICfg
                               STAMP%InstantTemp = sectCtrlParam%TI
                               call Putout_Instance_Config_SimMDBox(GFILE0, SimBox0K, STAMP)
                               !--- archive the REACTANT for restart purpose
                               !write(*,fmt="(A, I, A)") '                 to archive the configuration'
                               call Archive_Config_SimMDBox(SimBox0K, STAMP, CFILE0, Do_ArchiveProcess)

                               !--- Now we save the transition path
                               do IT = 1, ITT
                                  STAMP%ITime     = m_InterITIME(IT)
                                  STAMP%Time      = m_InterTIME(IT)
                                  STAMP%ScalTime  = m_SCALTIMETgt + (m_InterTIME(IT)-TIME0)*m_curReplicas
                                  STAMP%ICfg      = NPT
                                  STAMP%IRec      = STAMP%ICfg
                                  call Putout_Instance_Config_SimMDBox( GFILE, m_InterSimBox((IT-1)*m_mxReplicas+IB), STAMP)
                                  NPT = NPT + 1
                               end do
                               if(NCB .le. 0) then
                                  STAMP%ITime     = m_ITIMETgt
                                  STAMP%Time      = m_TIMETgt
                                  STAMP%ScalTime  = m_SCALTIMEtgt + (m_TIMEtgt-TIME0)*m_curReplicas
                                  STAMP%ICfg      = NPT
                                  STAMP%IRec      = STAMP%ICfg
                                  call Putout_Instance_Config_SimMDBox(GFILE, SimBoxTK, STAMP)
                               end if

                               !--- scal the time by current number of replicas
                               m_SCALTIMEtgt = m_SCALTIMEtgt + (m_TIMEtgt-TIME0)*m_curReplicas

                               !--- Now the detected PRODUCT configuration should be REACTANT of following evolution
                               call Copy_SimMDBox(SimBox0K, m_SimBoxIni(1))

                               if(m_NUMTgt .lt. EXITATEVENT ) then
                                  !--- in case of continuous evolution
                                  call Do_DePhase(m_SimBoxIni(1), SimBox(1:m_mxReplicas), gm_ForceClass)
                                  ITime      = m_ITIMETgt
                                  Time       = m_TIMETgt
                                  TIME0      = Time
                                  ITIMEP     = ITime - sectCtrlParam%IT0 + 1
                                  NCB        = 0
                                  NCBT       = 0
                                  !$$--- NOTE: in DEPhase, the number of particles on devices may be changed
                                  !$$          the number of replicas may also be changed
                                  !$$
                                  m_InterNum = 0
                                  call SetMultiBox_SimMDCtrl(CtrlParam, m_curReplicas)
                                  call Initialize_Globle_Variables_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                                  call Initialize_NeighboreList_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                                  call Cal_NeighBoreList_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                                  if(iand(sectCtrlParam%AR_METHOD,CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                                      call Initialize_ActiveRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                                      call ActivateRegion_DEV(SimBox(1:m_curReplicas), sectCtrlParam)
                                  end if

                                  call Init_Forcetable_Dev(SimBox(1:m_curReplicas), sectCtrlParam, gm_ForceClass)
                                  call CopyPointer_ForceClass(gm_ForceClass, m_ForceClassB)
                                  call ClearExtForce_ForceClass(m_ForceClassB)
                                  call CalForce_ForceClass(SimBox(1:m_curReplicas), sectCtrlParam, m_ForceClassB)
                               else
                                  Tstatu = IOR(NCB,CP_TIMECYCLE_END)
                               end if
                            end if
                          RESETINTER = .true.
                      end if

                      !--- recording the interdediate state
                      if(sectCtrlParam%PARREP_FineR .gt. 0) then
                         if(MOD(ITIMEP,sectCtrlParam%PARREP_FineR).eq.0) then
                            if(.not. HASCOPYOUT) then
                               call CopyOut_from_GPU(SimBox(1:m_curReplicas), sectCtrlParam)
                               HASCOPYOUT = .true.
                            end if

                            !--- save the intermediate configuration
                            if(.not. RESETINTER) then
                               m_InterNum = m_InterNum + 1
                               m_InterITime(m_InterNum) = ITime
                               m_InterTIME(m_InterNum)  = Time
                               do II=1, m_curReplicas
                                  call Copy_SimMDBox(SimBox(II), m_InterSimBox((m_InterNum-1)*m_mxReplicas+II))
                               end do
                            end if
                         end if
                      end if

                      !--- check if the time point is reached
                      if(ITime .GE. sectCtrlParam%IT1 .or. Time.GE.sectCtrlParam%TEMTM .or. &
                         IAND(TSTATU,CP_TIMECYCLE_END) .eq. CP_TIMECYCLE_END) then
                         exit
                      end if
                   end do   !end the loop for time step in current time section

                   !--- release memories in this time section
                   if(allocated(m_InterSimBox)) then
                      call Release_SimBoxArray(m_InterSimBox)
                      deallocate(m_InterSimBox, m_InterITime, m_InterTIME)
                   end if
                   call Release_SimMDBox(SimBox0K)
                   call Release_SimMDBox(SimBoxTK)

                   !--- to restore the original parameter
                   call Copy_SimMDCtrl(oldCtrlParam, sectCtrlParam)
      return
  end subroutine For_One_TimeSect
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Step(ITime, CtrlParam, SimBox)
  !***  PORPOSE: to go on for only one step
  !     INPUT:   I, the Ith time step
  !
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_DiffScheme_GPU
  use MD_LocalTempMethod_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU

  implicit none
  integer,                    intent(in)::ITime
  type(SimMDCtrl)                       ::CtrlParam
  type(SimMDBox), dimension(:)          ::SimBox

  !Local variables
  integer::ITIMEP

         ITIMEP = ITIME-CtrlParam%IT0+1
     !--- Give a prediction
         call Predictor_DEV(ITIME, SimBox(1), CtrlParam)

     !--- update the neighbore list
         if(mod(ITIMEP,CtrlParam%NB_UPTAB) .eq. C_IZERO) then
            call Cal_NeighBoreList_DEV(SimBox, CtrlParam)
         end if

     !--- update the the active region
          if(mod(ITIMEP, CtrlParam%AR_UPTAB) .eq. C_IZERO) then
              call ActivateRegion_DEV(SimBox, CtrlParam)
          end if

     !--- calculate the current force
         if(CtrlParam%IFBPC .EQ. 0) then
            call CalForce_ForceClass(SimBox, CtrlParam, m_ForceClassB)
         else
           !--- if pressure is needed, use the following subroutine
           call CalPTensor_ForceClass(SimBox, CtrlParam, m_ForceClassB)
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

     !--- to modify the statu of groups according their RSTATU
         !call Keep_Box_Zero_Translation_SimMDBox_DEV(m_Simbox)

     return
  end subroutine For_One_Step

  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Damp(SimBox, CtrlParam, ForceClass, Reinit)
   !***  PORPOSE: to damping the configuration to zero temperature.
   !
   !    INPUT:  CtrlParam,  the control parameters for simulation
   !            REINIT,     indicating if re-initialization of device memory is needed
   !    OUTPUT: SimBox,    the box that has been damped
   !    NOTE:   In do_damp, we use the force-class without boost-force. The same is true 
   !            for Dephase
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
        type(SimMDBox),dimension(:)            ::SimBox
        type(SimMDCtrl)                        ::CtrlParam
        type(MDForceClassGPU)                  ::ForceClass
        integer,                    intent(in) ::Reinit
       !Local variables
        integer::I, MULTIBOX, MXITER, NB, IFLAG
 
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
                 call CopyOut_SimBox_DEV(SimBox(1:NB))
               end if
               !*** to calculate the force
               !    The force to be calculated in call in DO_LBFGS_FORSTEPS_DEV
               !    It in not necessary to calculate force here
               call  CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
            end if
           
            select case(iand(CtrlParam%Quench_Meth, CP_LWORD))
                   case( CP_DAMPSCHEME_LBFGS)
                      call DO_LBFGSB_FORSTEPS_DEV  (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, IFLAG)
                   case(CP_DAMPSCHEME_DYN)
                      call Do_DynDamp_Forsteps_DEV (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps)
                   case(CP_DAMPSCHEME_ST )   
                      call Do_Steepest_Forsteps_DEV(SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
                   case(CP_DAMPSCHEME_CG )   
                      call Do_CG_Forsteps_DEV      (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
             end select
 
             !--- copy back the dampped configuration for analysis
             call CalEpot_ForceClass(SimBox,  CtrlParam, ForceClass)
             call CopyOut_SimBox_DEV(SimBox)
  
             !--- restore the originla control parameters
             CtrlParam%MULTIBOX  = MULTIBOX
  
   end subroutine Do_Damp
  !****************************************************************************************
 
  !****************************************************************************************
  subroutine Do_DePhase(SimBoxIni, SimBox, ForceClass)
  !***  PORPOSE: to dephase the configuration to zero temperature.
  !
  !     INPUT:  SimBoxIni, the reference configuration at zero temperature
  !
  !     OUTPUT: SimBox,    the box that has been thermalized
  !     NOTE:   On doing depahse, we use the force-class without boost-force. 
  !
  !
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU

  use MD_DiffScheme_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU

  IMPLICIT NONE
   !---dummy vaiables
       type(SimMDBox),              intent(in)  ::SimBoxIni
       type(SimMDBox), dimension(:)             ::SimBox
       type(MDForceClassGPU)                    ::ForceClass

      !Local variables
      integer::IB,  IBF0, ITIME, IVNUM, NPRT
      integer, dimension(:), allocatable::FLAG, MASK
      type(SimMDBox), dimension(:), allocatable::swapSimBox


           !*** to prepare the swapboxes
           NPRT     = SimBoxIni%NPRT
           allocate(swapSimBox(m_mxReplicas), Flag(NPRT*m_mxReplicas), Mask(NPRT))

           !*** to determine which atoms need to be included in
           !    checking structual change
           Mask = 1
           if(associated(m_pCfgCompMask) ) then
              call m_pCfgCompMask(SimBoxIni, m_CtrlParamDephase, Mask)
           end if

           !*** initialize the boxes
           do IB = 1, m_mxReplicas
               call Copy_SimMDBox(SimBoxIni, SimBox(IB))
           end do

           !*** to initialize the GPU variable
           call Initialize_Globle_Variables_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase)

           !*** to initialize force table on device
           call Init_Forcetable_Dev(SimBox(1:m_mxReplicas), m_CtrlParamDephase, gm_ForceClass)

           !*** to initialize the neigborlist module
           call Initialize_NeighboreList_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase)
           call Cal_NeighBoreList_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase)

           !*** if active region method to be used,  initialize the modeul
           if(iand(m_CtrlParamDephase%AR_METHOD,CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
               call Initialize_ActiveRegion_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase)
               call ActivateRegion_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase)
           end if

           !*** to calculate the force
           call CalForce_ForceClass(SimBox, m_CtrlParamDephase, ForceClass)

          !*** to thermalizing the system
          write(*,fmt="(A, I7, A)")  ' MDPSCU Message: to thermalize the replicas'

              IVNUM = 0
              do ITIME=1, (m_CtrlParamDephase%IVTIME+1)*m_CtrlParamDephase%IVPAS
                 if(MOD(ITIME-1, m_CtrlParamDephase%IVPAS) .eq. C_IZERO) then
                    call CalEpot_ForceClass(SimBox,  m_CtrlParamDephase, ForceClass)
                    call CalEkin_DEV(SimBox(1), m_CtrlParamDephase)
                    call CopyOut_SimBox_DEV(SimBox(1:m_mxReplicas))
                    call Cal_thermal_quantities_SimBoxArray(SimBox(1:m_mxReplicas))
                    write(*,fmt="(A, F12.3, A, I, A)") ' MDPSCU Message: temp. ', SimBox(1)%TEMPERATURE, ' K at', ITIME, ' timsteps in DEPHAS'

                    if(IVNUM .lt. m_CtrlParamDephase%IVTIME) then
                       call Thermalizing_MC_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase, m_CtrlParamDephase%TI)
                       IVNUM =  IVNUM + 1
                    end if
                 end if

                 !--- Give a prediction
                 call Predictor_DEV(ITIME, SimBox(1), m_CtrlParamDephase)

                !--- update the neighbore list
                 if(mod(ITIME, m_CtrlParamDephase%NB_UPTAB) .eq. C_IZERO) then
                    call Cal_NeighBoreList_DEV(SimBox(1:m_mxReplicas), m_CtrlParamDephase)
                 end if

                !--- update the the active region
                if(mod(ITIME, m_CtrlParamDephase%AR_UPTAB) .eq. C_IZERO) then
                   call ActivateRegion_DEV(SimBox, m_CtrlParamDephase)
                end if

                 !--- calculate the current force
                 call CalForce_ForceClass(SimBox, m_CtrlParamDephase, ForceClass)

                 !--- to do the correction
                 call Correction_DEV(ITIME, SimBox(1), m_CtrlParamDephase)
              end do

              !*** to check if all the boxes are in the same configuration
              write(*,fmt="(A, I7, A)")  ' MDPSCU Message: to check if the replicas are paralell'
              do IB = 1, m_mxReplicas
                 call Copy_SimMDBox(SimBox(IB), SwapSimBox(IB))
              end do

              call Do_Damp(SwapSimBox(1:m_mxReplicas), m_CtrlParamDamp, gm_ForceClass, REINIT=0)

              if(associated(m_pCfgComp) ) then
                 call m_pCfgComp(SimBoxIni, SwapSimBox(1:m_mxReplicas), m_CtrlParamDephase, Mask,Flag)
               else !--- do the default comparing
                 call Do_Compare(SimBoxIni, SwapSimBox(1:m_mxReplicas), m_CtrlParamDephase, Mask, Flag)
              end if

              IBF0 = 1
              do IB =1, m_mxReplicas
                 if(all(Flag((IB-1)*NPRT+1:IB*NPRT) .eq. 0) ) then
                    if(IB .gt. IBF0) call Copy_SimMDBox(SimBox(IB), SimBox(IBF0))
                    IBF0 = IBF0 + 1
                 else
                 end if
              end do

              m_curReplicas = IBF0-1
              write(*,fmt="(A,I,A)")  ' MDPSCU Message: ', m_curReplicas, ' replicas are succefully created'

           call Release_SimBoxArray(SwapSimBox)
           deallocate(SwapSimBox)
           deallocate(Flag, Mask)
     return
  end subroutine Do_DePhase
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_ChangeDetect(SimBoxIni, SimBox, CtrlParam, IBT, SimBox0K, NCB, REINITGPU)
  !***  PORPOSE: to detect the chang of configuration related to initial box SimBoxIni
  !
  !    INPUT:  SimBoxIni,    the 'standard' box to which the replicas to be compared
  !            SimBox,       the replicas
  !            CtrlParam,    the control parameters for simulation
  !            REINITGPU,    indicating if GPUs need to be re-initializing
  !
  !    OUTPUT: IBT,          the index of box that transition occur
  !            NBC,          the number of box(s) in which canges are detected
  !            SimBox0K,     the boxes queching to 0K. if only one box has been found changed. this
  !                          could be used for PRODUCT boxe
  !
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_DiffScheme_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU

  IMPLICIT NONE
   !---dummy vaiables
       type(SimMDBox),               intent(in)  :: SimBoxIni
       type(SimMDBox), dimension(:)              :: SimBox
       type(SimMDCtrl)                           :: CtrlParam
       integer                                   :: IBT
       type(SimMDBox)                            :: SimBox0K
       integer                                   :: NCB, REINITGPU

       !Local variables
       integer::IB, NB, NPRT
       integer, dimension(:), allocatable::Flag, Mask
       type(SimMDBox), dimension(:), allocatable::SwapBox

           !*** to prepare the swapboxes
           NB    = size(SimBox)
           NPRT  = SimBoxIni%NPRT

           allocate(SwapBox(NB),Flag(NPRT*NB), Mask(NPRT))
           do IB = 1, NB
               call Copy_SimMDBox(SimBox(IB), SwapBox(IB))
           end do

           call Do_Damp(SwapBox, m_CtrlParamDamp, gm_ForceClass, REINIT=REINITGPU)
           !---
           Mask = 1
           if(associated(m_pCfgCompMask) ) then
              call m_pCfgCompMask(SimBoxIni, CtrlParam, Mask)
           end if

           if(associated(m_pCfgComp) ) then
              call m_pCfgComp(SimBoxIni, SwapBox, CtrlParam, Mask, Flag)
           else !--- do the default comparing
              call Do_Compare(SimBoxIni, SwapBox, CtrlParam, Mask, Flag)
           end if

           NCB = 0
           IBT = 0
           do IB =1, NB
              if(any(Flag((IB-1)*NPRT+1:IB*NPRT) .gt. 0) ) then
                 IBT = IB
                 NCB = NCB + 1
                 call Copy_SimMDBox(SwapBox(IB), SimBox0K)
              end if
           end do

           !--- to restore the device statu  before damping
           call  CopyIn_SimBox_DEV(SimBox)
           call  Cal_NeighBoreList_DEV(SimBox, CtrlParam)

           call  Release_SimBoxArray(SwapBox)
           deallocate(SwapBox,Flag, Mask)
           return
     return
  end subroutine Do_ChangeDetect
  !****************************************************************************************

  !****************************************************************************************
  subroutine RefineTranStat(SimBoxIni, CtrlParam, ITT, IB, SimBox0K, NCB)
  !***  PORPOSE: to refine the time point at which the transition occur
  !
  !    INPUT:  SimBoxIni,    the 'standard' box to which the replicas to be compared
  !            SimBox,       the replicas
  !            CtrlParam,    the control parameters for simulation
  !
  !    OUTPUT: IBT,          the index of box that transition occur
  !            NBC,          the number of box(s) in which changes are detected
  !            SimBox0K,     the boxes queching to 0K. if only one box has been found changed. this
  !                          could be used for PRODUCT boxes
  !
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_DiffScheme_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU

  IMPLICIT NONE
   !---dummy vaiables
       type(SimMDBox), intent(in):: SimBoxIni
       type(SimMDCtrl)           :: CtrlParam
       integer                   :: ITT
       integer                   :: IB
       type(SimMDBox)            :: SimBox0K
       integer                   :: NCB

      !Local variables
      type(SimMDBox)::SwapBox0K
      integer::IT, IT0, IT1, IBT, NCBT

           !*** to prepare the swapboxes
           if(m_InterNum .le. 0 ) then
             NCB = 0
             return
           end if
           NCB   = 0
           IT0   = 1
           IT1   = m_InterNum

           do while(.true.)
              IT    = (IT0 + IT1)/2
              write(*,fmt="(A, I7, 4I5)")  ' MDPSCU Message: to detect the structural change at ITIME=',  m_InterITIME(IT)
              call Do_ChangeDetect(SimBoxIni, m_InterSimBox((IT-1)*m_mxReplicas+1:(IT-1)*m_mxReplicas+m_curReplicas), CtrlParam, IBT, SwapBox0K, NCBT, REINITGPU=1)
              write(*,fmt="(A, I7, A)")    ' MDPSCU Message: structural change detected in ', NCBT, ' boxes'

              !--- save the index of for the at which changed is detected
              if(NCBT .gt. 0) then
                 NCB = NCBT
                 ITT = IT
                 IB  = (IT-1)*m_mxReplicas + IBT
                 call Copy_SimMDBox(SwapBox0K, SimBox0K)
              end if

              if(IT .eq. IT0) exit

              if(NCBT .gt. 0) then
                 IT1 = IT
              else
                 IT0 = IT
              end if
           end do

           call Release_SimMDBox(SwapBox0K)
           return
     return
  end subroutine RefineTranStat
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Compare(SimBoxIni, SimBox, CtrlParam, MASK, Flag)
  !***  PORPOSE: to check if the SimBox is different from SimBoxIni
  !
  !    INPUT:  SimBoxIni,    the 'standard' box to which the replicas to be compared
  !            SimBox,       the replicas
  !            CtrlParam,    the control parameters for simulation
  !            MASK,         indicating which atoms to be included in the comparing
  !
  !    OUTPUT: FLAG,         indicating which atoms have been changed in position
  !
  !
  implicit NONE
   !---dummy vaiables
       type(SimMDBox),               intent(in)::SimBoxIni
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       integer,        dimension(:), intent(in)::MASK
       integer,        dimension(:)            ::Flag
       !Local variables
       integer::IB, NPRT, I, IP, IFPD(3)
       real(KINDDF)::RC2, SEP(3), BOX(3), HBOX(3)

               RC2  =  CtrlParam%STRCUT_DRTol*CtrlParam%STRCUT_DRTol
               NPRT = SimBoxIni%NPRT
               BOX  = SimBoxIni%ZL
               HBOX = SimBoxIni%ZL*C_HALF
               IFPD = CtrlParam%IFPD

               IP   =  0
               do IB=1, size(SimBox)
                  do I=1, NPRT
                     IP = IP + 1
                     Flag(IP) = 0
                     if(MASK(I) .le. 0 ) cycle

                     SEP(1:3)  =  SimBoxIni%XP(I, 1:3) - SimBox(IB)%XP(I,1:3)

                     if( (IFPD(1).GT.0) .AND. (DABS(SEP(1)) .GT. HBOX(1))) then
                        SEP(1) = SEP(1) - DSIGN(BOX(1),SEP(1))
                     end if

                     if( (IFPD(2).GT.0) .AND. (DABS(SEP(2)) .GT. HBOX(2))) then
                        SEP(2) = SEP(2) - DSIGN(BOX(2),SEP(2))
                     end if

                     if( (IFPD(3).GT.0) .AND. (DABS(SEP(3)) .GT. HBOX(3))) then
                        SEP(3) = SEP(3) - DSIGN(BOX(3),SEP(3))
                     end if

                     if(sum(SEP*SEP) .gt. RC2) then
                        Flag(IP)= 1
                      end if

                  end do
               end do

  end subroutine Do_Compare
  !****************************************************************************************


  end module MD_Method_ParRep_GPU
