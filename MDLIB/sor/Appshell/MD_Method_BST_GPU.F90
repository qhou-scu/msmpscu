  module MD_Method_BST_GPU
  !***  DESCRIPTION:
  !     This module provides application interfaces of BST - boost static calculation. 
  !     In BST, the initial configuration are queched to zero temperature, then boost-potentials
  !     are added repeatly, and new states are searched for. 
  !
  !                  _________________________________________________________________________
  !**** HISTORY:     2019-02 (HOU Qing), created the first version
  ! 
  !
  !***  The modules included ******************************************
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_RecordStamp
  use MD_TYPEDEF_RecordList
  implicit none

      type(SimMDCtrl),                           private:: m_CtrlParamDamp          !parameters control damping operation
      integer,                                   private:: m_mxReplicas  = 0        !max number of replicas
      integer,                                   private:: m_curReplicas = 0        !current number of replicas
      type(SimMDBox), dimension(:), allocatable, private:: m_SimBoxIni              !the initial configuration at 0K, e.g., the RECTANT
      integer,        dimension(:), allocatable, private:: m_ITEvent                !the time step at which the last event occur
      integer,        dimension(:), allocatable, private:: m_NumEvent               !the number of accepted events

      !--- the interface for selector
       private::SELECT_PROC
       abstract interface
          subroutine SELECT_PROC(CtrlParam, SimBox0, NC, IC, XP0)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDCtrl)                           ::CtrlParam
          type(SimMDBox),dimension(:)               ::SimBox0
          integer                                   ::NC
          integer,       dimension(:),  allocatable ::IC
          real(KINDDF),  dimension(:,:),allocatable::XP0
       end subroutine SELECT_PROC
       end interface
       procedure(SELECT_PROC), pointer, private::m_pSelector=>null()

  contains
  !****************************************************************************************
  !****************************************************************************************
  subroutine SetBoostSeedProc(extBoostSeedPro)
  !***  DESCRIPTION: to set the procedure that used to construct active region
  !
  !--- INPUT: extBoostSeedPro,   the external subroutine of create the candidate atom
  !                            for boost
  !
  implicit none

  !--- dummy variables and subroutines
   optional::extBoostSeedPro
   external::extBoostSeedPro

  !--- interface to the external routine -------------------
   interface
       subroutine extBoostSeedPro(CtrlParam, SimBox0, NC, IC, XP0)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl   
          type(SimMDCtrl)                           ::CtrlParam
          type(SimMDBox),dimension(:)               ::SimBox0
          integer                                   ::NC
          integer,       dimension(:),  allocatable ::IC
          real(KINDDF),  dimension(:,:),allocatable::XP0
       end subroutine extBoostSeedPro 
   end interface 
  !--- END INTERFACE --------------------------------

  !--- Local variables
          if(present(extBoostSeedPro)) then
             m_pSelector=>extBoostSeedPro
          else
             m_pSelector=>null()
          end if

          return
      end subroutine SetBoostSeedProc
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
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU
  use MD_DiffScheme_GPU,        only:CalEKin_DEV
  use MD_LocalTempMethod_GPU,   only:Do_ResetParam_DEV
  use MD_TYPEDEF_PrintList,     only:Do_RestoreProcess
  implicit none
  !--- dummy variables
   type(SimMDBox),                           intent(in)::SimBox0
   type(SimMDCtrl),target                              ::CtrlParam
   type(SimMDBox),dimension(:),allocatable             ::SimBox
   type(RecordProcedureList)                           ::Recordlist
   integer,                                  intent(in)::J,processid

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
  integer::ITIME0, ITIME, TESTLOOP0, TESTLOOP, ITEST, JOB, ISECT0, ISECT, TSTATU, IB
  logical::EXIST

  character*256::CFile, CFile0
  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2
  type(SimMDCtrl), pointer::sectCtrlParam, nextCtrlParam
  type(MDRecordStamp)::STAMP, STAMP0

  !----
             STAMP%AppTYpe = "BST"

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
              call Create_SimBoxArray(SimBox,      CtrlParam%MULTIBOX,SimBox0)
              call Create_SimBoxArray(m_SimBoxIni, CtrlParam%MULTIBOX,SimBox0)
              allocate(m_ITEvent(CtrlParam%MULTIBOX)) 
              m_ITEvent = 0

              !*** to load initial configuration
              ITIME0 = 0
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
                      call Initialize_Config_SimBoxArray(SimBox,CFile0,JOB,Simbox0%IniCfgID,SimBox0%IniCfgFmt)
                   end if   
                 end if
                 !--- initialize the configuration by user suppied routine
                 if(present(INICONFIGPROC))then
                    call INICONFIGPROC(SimBox, CtrlParam, RESTART=0)
                 endif

              !*** if we restart from a previous, restore from previous stop point
              else
                 if(present(INICONFIGPROC))then
                    if(CtrlParam%INDEPBOX) then
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=1)
                    else
                       call INICONFIGPROC(SimBox, CtrlParam, RESTART=JOB)
                    end if
                 end if
                 print *, "Restore configure from ", CFILE(1:len_trim(CFILE))
                 print *, "........ "
                 call Restore_Config_SimBoxArray(SimBox, STAMP0, CFile, Do_RestoreProcess)
                 ITIME0 = STAMP0%ITime
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
                 JOB       = (ITEST-1)+J
                 ITIME     = ITIME0
                 ISECT     = ISECT0

                 !*** to recorde the start time
                 call CPU_TIME(C1)

                 !*** to determine the start time section firstly
                 call GetNext_SimMDCtrl(CtrlParam, ISECT, sectCtrlParam)

                 !*** prepare the control parameters for damping
                 call Copy_SimMDCtrl(sectCtrlParam, m_CtrlParamDamp)
                 m_CtrlParamDamp%DAMPTIME0     = 1
                 m_CtrlParamDamp%DAMPTIME1     = sectCtrlParam%DAMPTIME1
                 m_CtrlParamDamp%IVTIME        = 0
                 if(m_CtrlParamDamp%DAMPTIME1 .le. 1) then
                     m_CtrlParamDamp%DAMPTIME1 = 1000
                 end if

                 !*** to save the inital stable states 
                 if(ITIME .le. 0) then
                    call Do_Damp(SimBox, CtrlParam, gm_ForceClass, Reinit=1)
                    call CopyOut_SimBox_DEV(SimBox)
                    m_SimBoxIni     = SimBox
                    STAMP%ITest     = JOB
                    STAMP%IBox(1)   = (JOB-1)*CtrlParam%MULTIBOX + 1
                    STAMP%IBox(2)   = JOB*CtrlParam%MULTIBOX 
                    STAMP%ITime     = ITime
                    STAMP%Time      = ITime
                    STAMP%ISect     = ISect
                    STAMP%ICfg      = 0
                    STAMP%IRec      = 0
                    call STRCATI(CFILE0, CtrlParam%f_geometry, "_0K_P", processid, 4)
                    call Putout_Instance_Config_SimBoxArray(CFILE0, 0, m_SimBoxIni, STAMP)
                 end if 
   
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
  use MD_ForceLib_Factory_GPU
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
    type(MDForceClassGPU)::BoostForceClass
  !----
              STAMP%AppTYpe = "BST"

              !*** start the event searching              
              IT = max(ITime, CtrlParam%IT0)
              if(CtrlParam%ExitAtEvent .le. 0) then
                  EXITATEVENT = CtrlParam%IT1
              else
                  EXITATEVENT = CtrlParam%ExitAtEvent
              end if
              call CopyPointer_ForceClass(gm_ForceClass, BoostForceClass)
              do while(IT .le. CtrlParam%IT1)
                 call Do_AddBoostPot(CtrlParam,  SimBox, BoostForceClass)
                 call Do_Damp(SimBox, CtrlParam, BoostForceClass, Reinit=0)
                 !call Do_EventCheck(m_SimBoxIni, SimBox)
                 IT = IT + 1
                 ITime = IT
              end do
      return
  end subroutine For_One_TimeSect
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Damp(SimBox, CtrlParam, ForceClass, Reinit)
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
   
   implicit none
   !---dummy vaiables
       type(SimMDBox),dimension(:)            ::SimBox
       type(SimMDCtrl)                        ::CtrlParam
       type(MDForceClassGPU),       intent(in)::ForceClass
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
             call Init_Forcetable_Dev(SimBox(1:NB), CtrlParam, ForceClass)

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
              call  CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
           end if

           select case(iand(CtrlParam%Quench_Meth, CP_LWORD))
                  case( CP_DAMPSCHEME_LBFGS)
                     call DO_LBFGSB_FORSTEPS_DEV  (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, IFLAG)
                  !case(CP_DAMPSCHEME_DYN)
                  !   call Do_DynDamp_Forsteps_DEV (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps)
                  case(CP_DAMPSCHEME_ST )   
                     call Do_Steepest_Forsteps_DEV(SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
                  case(CP_DAMPSCHEME_CG )   
                     call Do_CG_Forsteps_DEV      (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
            end select

            !--- copy back the dampped configuration for analysis
            call CalEpot_ForceClass(SimBox,  CtrlParam, ForceClass)
 
            !--- restore the originla control parameters
            CtrlParam%MULTIBOX  = MULTIBOX
 
  end subroutine Do_Damp
  !****************************************************************************************

  !****************************************************************************************
  subroutine SelectCandidates(CtrlParam, SimBox, NC, IC, XP0)
  !***  PORPOSE: to select the candidates to be set initial displacement
  !     INPUT:  
  !             CtrlParam,          the control parameters
  !             SimBox,             the current configurations
  !     OUTPUT: NC,                 the number of cancdiate
  !             IC,                 the index of the candidate atoms 
  !
  use RAND32_MODULE
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                            ::CtrlParam
   type(SimMDBox),dimension(:)                ::SimBox
   integer                                    ::NC
   integer,       dimension(:),   allocatable ::IC
   real(KINDDF),  dimension(:,:), allocatable ::XP0 
  !--- Local variables
   integer::IA, IP, NA0, NA, IB, NPRT0, I1, NM
   integer, dimension(:), allocatable::ICFLAG

   real(kinddf)::P0, P

                 if(allocated(IC))  deallocate(IC)
                 if(allocated(XP0)) deallocate(XP0)

                 if(iand(CtrlParam%ART_ActMeth, CP_CENTPARTACT_ART) .eq. CP_CENTPARTACT_ART) then
                    if(dabs(CtrlParam%ART_ActExt) .le. 1.D-10) then
                       call SelectCandidates0(CtrlParam, SimBox, ICFLAG) 
                    else
                       call SelectCandidates1(CtrlParam, SimBox, ICFLAG) 
                    end if
                 else
                 end if 

                !--- the number of candidates to be slected
                 NM = count(ICFLAG .gt. 0)
                 allocate(IC(NM), XP0(NM,3))
                 IC = 0
                 I1 = 0
                 NC = 0 
                 do IB = 1, size(SimBox) 
                    NPRT0 = SimBox(IB)%NPRT
                    NA0 = count(ICFLAG(I1+1:I1+NPRT0) .gt. 0)
                    if(iand(CtrlParam%ART_ActMeth, CP_MULTPLEACT_ART) .eq. CP_MULTPLEACT_ART) then   
                       NA = nint(CtrlParam%ART_ActPecent*NA0)
                    else
                       NA = 1
                    end if   

                    !--- the probability for an atom to be selected
                    P0 = dble(NA)/dble(NA0)
                    do while(NC .lt. NA)
                       do IA =1, NPRT0
                          IP = I1 + IA 
                          if( ICFLAG(IP) .le. 0) cycle
                          P = DRAND32()
                          if(P .le. P0) then
                             NC          = NC + 1
                             IC(NC)      = IA
                             XP0(NC,1:3) = SimBox(IB)%XP(IA,1:3)
                             ICFLAG(IP)  = 0
                             if(NC .ge. NA) exit
                           end if
                       end do
                    end do 
                    I1 = I1 + NPRT0  
                 end do

                 if(NC .le. 0) then
                    write(*,fmt="(A)")  " MDPSCU Error: no any atom to be activated"
                    write(*,fmt="(A)")  "               check control parameters for ART"
                    write(*,fmt="(A)")  " Process to be stopped"
                    stop
                 else
                    write(*,fmt="(A, I6, A)")  " MDPSCU Message: ", NC, " atom(s) are activated"
                 end if

                 if(allocated(ICFLAG)) deallocate(ICFLAG)
                 return
   end subroutine SelectCandidates
  !****************************************************************************************

  !****************************************************************************************
  subroutine SelectCandidates0(CtrlParam, SimBox, ICFLAG)
  !***  PORPOSE: to select the candidate to be set initial displacement
  !              the candicates are those selected by given type
  !     INPUT:  
  !             CtrlParam,          the control parameters
  !             SimBox0,            the initial configuration
  !     OUTPUT: ICFLAG,             the index of the candidate atoms 
  !
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                           ::CtrlParam
   type(SimMDBox), dimension(:)              ::SimBox
   integer,        dimension(:), allocatable ::ICFLAG
  !--- Local variables
   integer::IP, IA, IB
                    
                    !--- the total number of atoms involved in the selection
                    allocate(ICFLAG(SimBox(1)%NPRT*size(SimBox)))
                    ICFLAG = 0
                    IP     = 0
                    do IB = 1, size(SimBox)
                       do  IA= 1, SimBox(IB)%NPRT
                           IP= IP + 1
                           if( CtrlParam%ART_ActSeed(SimBox(IB)%ITYP(IA)) .le. 0) cycle
                           if(iand(SimBox(IB)%STATU(IA), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) cycle
                           ICFLAG(IP) = 1
                       end do    
                    end do

                 return
   end subroutine SelectCandidates0
  !****************************************************************************************

  !****************************************************************************************
  subroutine SelectCandidates1(CtrlParam, SimBox, ICFLAG)
  ! *** PORPOSE: to select the candidate to be set initial displacement
  !              the candicates are those selected by given type
  !     INPUT:  
  !             CtrlParam,  the control parameters
  !             SimBox,     the initial configuration
  !     OUTPUT: ICFLAG,             the index of the candidate atoms 
  !
  use MD_NeighborsList_GPU
  implicit none
  !--- dummy variables
   type(SimMDCtrl)                          ::CtrlParam
   type(SimMDBox),dimension(:)              ::SimBox
   integer,       dimension(:), allocatable ::ICFLAG
  !--- Local variables
   integer::IP, IB, IA, IN
   type(NEIGHBOR_LIST), allocatable::tList
   type(SimMDCtrl),  allocatable::tCtrlParam

                 !*** to create a neighbore list   
                 allocate(tList, tCtrlParam)                 
                 tCtrlParam       = CtrlParam
                 tCtrlParam%NB_RM = dabs(CtrlParam%ART_ActExt)*SimBox(1)%RR
                 call Initialize_Globle_Variables_DEV(SimBox, tCtrlParam)
                 call Initialize_NeighboreList_DEV(SimBox, tCtrlParam)
                 call Cal_NeighBoreList_DEV(SimBox, tCtrlParam)
                 call Copyout_NeighboreList_DEV(tList, ORDER=1)

                    !--- the total number of atoms involved in the selection
                    allocate(ICFLAG(SimBox(1)%NPRT*size(SimBox)))
                    ICFLAG = 0
                    IP     = 0
                    if(CtrlParam%ART_ActExt .gt. 0) then
                       do IB =1, size(SimBox) 
                          do IA = 1, SimBox(IB)%NPRT
                             IP = IP + 1
                             if( CtrlParam%ART_ActSeed(SimBox(IB)%ITYP(IA)) .le. 0) cycle
                             if(iand(SimBox(IB)%STATU(IA), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) cycle
                             ICFLAG(IP) = 1
                             do IN = 1, tList%KVOIS(IP)
                                ICFLAG(tList%INDI(IP, IN)) = 1 
                             end do
                          end do
                        end do
                     else
                        ! --- to active the neighbores of the seeds, but the seeds are deactivated 
                       do IB =1, size(SimBox) 
                          do IA = 1, SimBox(IB)%NPRT
                             IP = IP + 1
                             if( CtrlParam%ART_ActSeed(SimBox(IB)%ITYP(IA)) .le. 0) cycle
                             if(iand(SimBox(IB)%STATU(IA), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) cycle
                             ICFLAG(IP) = 0
                             do IN = 1, tList%KVOIS(IP)
                          	    ! --- the neighbores of the active seed cannot be active
                           	 if( CtrlParam%ART_ActSeed(SimBox(IB)%ITYP(tList%INDI(IP, IN))) .gt. 0) cycle
                               ICFLAG(tList%INDI(IP, IN)) = 1 
                             end do
                          end do   
                       end do
                     end if
                 call Clear_NeighboreList_DEV()
                 call Clear_NeighboreList(tList)
                 deallocate(tList, tCtrlParam)
                 return
   end subroutine SelectCandidates1
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_AddBoostPot(CtrlParam, SimBox, BoostForceClass)
  !***  PORPOSE: to add a boost potential to BoostForceClass
  !     INPUT:  
  !             CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !
  use MD_ForceLib_Factory_GPU
  use MD_BoostMeth_GPU
  implicit none
  !--- dummy variables
   type(SimMDCtrl)               ::CtrlParam
   type(SimMDBox),dimension(:)   ::SimBox
   type(MDForceClassGPU)         ::BoostForceClass
  !--- Local variables
       integer::NC
       integer,      dimension(:),   allocatable::IC
       real(KINDDF), dimension(:,:), allocatable::XP0
       !------

                 if(CtrlParam%ART_ActMeth .eq. CP_USERDEFACT_ART) then
                    if(.not.associated(m_pSelector)) then
                       write(*,fmt="(A)")  " MDPSCU Error: in ART method. user-defined procedure for"
                       write(*,fmt="(A)")  "                selecting activation candidates is null"
                       write(*,fmt="(A)")  " Process to be stopped"
                       stop
                    end if
                    call m_pSelector(CtrlParam, SimBox, NC, IC, XP0)
                 else
                    call SelectCandidates(CtrlParam, SimBox, NC, IC, XP0)
                 end if 
                 call New_BoostForce(SimBox, CtrlParam, NC, IC, XP0)

                 if(allocated(IC)) deallocate(IC)
                 return
   end subroutine Do_AddBoostPot
  !****************************************************************************************
  end module MD_Method_BST_GPU
