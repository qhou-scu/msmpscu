  module MD_SimboxArray_AppShell_16_GPU
  !***  DESCRIPTION:
  !     This program is a GPU version(2014) providing a common shell for
  !     construting an generic MD apllications. To use GPU efficiently,
  !     the calculations are performed for multiple simulation boxs which
  !     are independent from each others but share the box same size and other
  !     simulation conditions. The number of boxes in one simulation is
  !     controlled by m_CtrlParam%MULTIBOX
  !
  !
  !     SEE ALSO: MD_EM_TB_ForceTable_SHELL_12_GPU.f90
  !
  !     HISTORY: This version is adapted from  MD_EM_TB_ForceTable_SHELL_12_GPU.f90
  !
  !     by HOU Qing, June, 2014

  !***  The modules included ******************************************
  use MD_CONSTANTS
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_RecordList
  use MD_MethodClass_Register_GPU
  use MiniUtilities
  implicit none

  !___the global data _____________________________________________________________________
  !
      type(SimMDBox), dimension(:), allocatable::m_SimBox
      type(SimMDCtrl),  target                 ::m_CtrlParam
      type(MDMethodClassGPU)                   ::m_Method
      type(RecordProcedureList)                ::m_Recordlist

  !_______________________________________________________________________________________
  !
  contains
  !****************************************************************************************
  subroutine APPSHELL_Main(NMPI,processid, FORCETABLE, POTTYPE, INICONFIGPROC)
  !***  DESCRIPTION:
  !     The main process for in MDPSCU
  !
  !--- INPUT: NMPI, the number of MPI process
  !           processid,    the process ID
  !           FORCETABLE,   optional, the subroutine provided by user for force register
  !           POTTYPE,      optional, the identifier of the type of potential
  !           INICONFIGPROC,the subroutine privided by user for initializing the initial configure
  !
  use RAND32SEEDLIB_MODULE
  use RAND32_MODULE
  use MD_Globle_Variables
  use MD_NeighborsList_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU
  use MD_ActiveRegion_GPU,    only:Clear_ActiveRegion_DEV
  use MD_TYPEDEF_PrintList,   only:Do_PrintProcess, Add_ArchiveProcess, Add_RestoreProcess
  implicit none

  !--- dummy variables and subroutines
  integer,intent(in)::NMPI,processid
  character*(*), optional::POTTYPE
  !------------------------------------------
  !--- interface to the external routine -------------------
  !--- interface to the external routine -------------------
   optional::FORCETABLE, INICONFIGPROC
   external::FORCETABLE, INICONFIGPROC
   interface
     subroutine FORCETABLE(SimBox, CtrlParam, FTable, printout)
     use MD_CONSTANTS
     use MD_TYPEDEF_SimMDBox
     use MD_TYPEDEF_SimMDCtrl
     use MD_TYPEDEF_ForceTable
     implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in)   ::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(inout)::FTable
      integer,optional::printout
     end subroutine FORCETABLE
   end interface

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
  !--- END INTERFACE --------------------------------

  !------------------------------------------

  !--- Local variables
  integer::I, J,K,IDEV=0, NP=1, NDEV=1, TESTLOOP, NEXT, numdevices, IERR

  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2
  type(SimMDBox)      ::SimBox0
  type(MDRecordStamp) ::STAMP
  integer::ISEED0, ISEED(2)
  character*32::ARGV, tstr


      !*** to initialize the DEVICE
          J = COMMAND_ARGUMENT_COUNT()
          if(J.GE.2) then
             call GET_COMMAND_ARGUMENT(2,ARGV)
             read(ARGV, *) IDEV
          end if

          if(J.GE.3) then
             call GET_COMMAND_ARGUMENT(3,ARGV)
             read(ARGV, *) NDEV
          end if

          ierr = cudaGetDeviceCount(numdevices)
          if(IDEV .gt. numdevices-1) then
             write(*,*) "MDPSCU Error: the ID of the first GPU is larger than the available value."
             write(*,fmt="(A,1x,I3)") "              the number of GPUs on this machine is", numdevices
             stop
          end if

          if(NDEV .gt. numdevices-IDEV) then
             write(*,*) "MDPSCU Error: the number of GPUs to be used is larger than the available value."
             write(*,fmt="(A,1x,I3)") "              the number of GPUs on this machine", numdevices
             write(*,fmt="(A,1x,I3)") "              the ID of the first GPU is", IDEV
             write(*,fmt="(A,1x,I3)") "              the number of GPUs to be used cannot be larger than", numdevices-IDEV
             stop
          end if

          if(NDEV .le. 0) then
             write(*,*) "MDPSCU Error: the number of GPUs to be used is zero."
             stop
          end if
          call Initialize_DEVICES( IDEV, NDEV)

      !***
          gm_NMPI = NMPI
          NP = NMPI/2
          IF(NP .LT. 1) NP = 1


      !*** to recorde the start time
          !CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
          !C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000+DATE_TIME2(4)*3600*1000*24
          call CPU_TIME(C1)
          call OpenLogFile_Globle_Variables( )

      !*** to read in control parameters from IniFile
   1000   call Release_SimMDCtrl(m_CtrlParam)
          call Initialize_Globle_Variables(SimBox0, m_CtrlParam, NEXT=next)

          if(next .LE. 0) then
             write(tstr, *) iabs(next)
             tstr = adjustl(tstr)
             write(*,*) "Simulations stop after finishing job "//tstr(1:len_trim(tstr))

             !CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
             !C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000+DATE_TIME2(4)*3600*1000*24
             CALL CPU_TIME(C2)
             print *, "RUN TIME: ",C2-C1

             call End_DEVICES()
             return
          end if

      !*** to initialize the force-potential table
          if(present(FORCETABLE) .or. present(POTTYPE)) then
              if(.not. present(POTTYPE)) then
                 write(*,fmt="(' MDPSCU Error: type of user-supplied potential is not given')")
                 write(*,fmt="('               supported optential include: ')")
                 do I=1, size(PKW_POTTYPELIST)
                    write(*,fmt="('               ', A)") '"'//PKW_POTTYPELIST(I)(1:len_trim(PKW_POTTYPELIST(I)))//'"'
                 end do
                 write(*,fmt="('               check the code calling APPSHELL_Main')")
                 write(*,fmt="('               Process to be stopped')")
                 stop
              end if
              SimBox0%POTTYPE = POTTYPE(1:len_trim(POTTYPE))

              if(len_trim(SimBox0%potlibname) .gt. 0) then
                 write(*,fmt="(A,A,A,A)") 'MDPSCU Warning: used-supplied potential register of type ',     &
                                           POTTYPE(1:len_trim(POTTYPE)),                                   &
                                          ' to be used in stead of the potential lib given in box file: ', &
                                           SimBox0%potlibname(1:len_trim(SimBox0%potlibname))
                 call ONWARNING(gm_OnWarning)
              end if
              SimBox0%potlibname = "USER-SUPPLIED"
              if(present(FORCETABLE)) then
                 call RegUser_ForceLib(SimBox0, m_CtrlParam, gm_ForceClass, FORCETABLE)
              else   
                 call RegUser_ForceLib(SimBox0, m_CtrlParam, gm_ForceClass)
              end if   
          else
              if(len_trim(SimBox0%POTTYPE).le.0 .or. len_trim(SimBox0%potlibname).le.0) then
                 write(*,fmt="(A)") ' MDPSCU Error: potential type or potential lib name is missed.'
                 write(*,fmt="(A)") '               Check box file for potential type, '
                 write(*,fmt="(A)") '               Process to be stopped'
                 stop
              end if
              call Load_ForceLib(SimBox=SimBox0, CtrlParam=m_CtrlParam, Forceclass=gm_ForceClass)
          end if

      !*** to initialize the configuration
          call Initialize_SimMDBox(SimBox0,m_CtrlParam%DiffOrder)

     !*** to initialize the method
          call Register_MethodClass(gm_AppType, m_Method, SimBox0, m_CtrlParam)

      !*** to intialize the random number seed to call the random number generator once
           ISEED0 = m_CtrlParam%SEED(1)
           call GetSeed_RAND32SEEDLIB(ISEED0, ISEED(1), ISEED(2))
           ISEED0=ISEED0 + processid - 1
           call GetSeed_RAND32SEEDLIB(ISEED0, ISEED(1), ISEED(2))
           call DRAND32_PUTSEED(ISEED)
           call Initialize_Rand_DEVICES()
           call Add_ArchiveProcess("cuRanStat", Archive_Rand_DEVICES)
           call Add_RestoreProcess("cuRanStat", Restore_Rand_DEVICES)
           call Print1_Globle_Variables(6, SimBox0, m_CtrlParam)

      !**** if user provide preprocessor, we do it
            call DoPrerecord_List(m_Recordlist, SimBox0, m_CtrlParam)

      !**** if user provide preprocessor, we print then out it
           call Do_PrintProcess(6, SimBox0, m_CtrlParam)
           if(gm_hFILELOG.gt.0) call Do_PrintProcess(gm_hFILELOG, SimBox0, m_CtrlParam)


      !**** TO BEGINE LOOP ON SAMPLES
           if(m_CtrlParam%INDEPBOX) then
              TESTLOOP = m_CtrlParam%TOTALBOX/m_CtrlParam%MULTIBOX/NP
           else
              TESTLOOP = 1
           end if

           do J=1, TESTLOOP
              !*** do some preprocess before one test
              call DoBeforeOneTest_List(m_Recordlist, J, SimBox0, m_SimBox, m_CtrlParam)

              !*** in case user have change the simulation box,we have to check the completeness of the box
              call CheckForceTableAssigned_SimMDBox(SimBox0)
              if(allocated(m_SimBox) )call CheckForceTableAssigned_SimMDBox(m_SimBox)

              !*** to start the sumulation method
              if(present(INICONFIGPROC)) then
                  call m_Method%pForOnetest(SimBox0, m_CtrlParam, m_SimBox, m_Recordlist, J, processid, INICONFIGPROC=INICONFIGPROC)
              else
                  call m_Method%pForOnetest(SimBox0, m_CtrlParam, m_SimBox, m_Recordlist, J, processid)
              end if
              !--- to process the data for this test
              STAMP%AppType = gm_AppType
              STAMP%ITest   = J
              call DoAfterOneTest_List(m_Recordlist, SimBox0, STAMP, m_SimBox, m_CtrlParam)

           end do  !end the loop for samples

      !**** if external after-process provided, we do it
       call DoAfterRecord_List(m_Recordlist, SimBox0, m_CtrlParam)

      !---clear the momery allocated on device
       call Clear_NeighboreList_DEV()
       call Clear_ActiveRegion_DEV()
       call Clear_Globle_Variables_DEV()
       call Clear_ForceLib(gm_ForceClass)
       call Clear_MethodClass(m_Method)
       if(allocated(m_SimBox)) then
          call Release_SimBoxArray(m_SimBox)
          deallocate(m_SimBox)
       end if
       call Release_SimMDCtrl(m_CtrlParam)
       !**** to performat next simulation
       goto 1000


    return
  end subroutine APPSHELL_Main
  !****************************************************************************************

  !****************************************************************************************
  subroutine APPSHELL_AddRecord(AFTONESTEP, AFTONETEST,AFTRECORD, BEFONETEST, PRERECORD, RECORDPROC, RECSTATU, RESTARTREC, SAVERECSTATU)
  !***  DESCRIPTION: to add the the recording process to the shell
  !
  !--- INPUT: PRERECORD,    the subroutine provided by user to pre-process before the main-loop for tests
  !           RECORDPROC,   the subroutine provided by user to process data at time when the configure to be output
  !           BEFOREONETEST,the subroutine to provided by user to process data before one test
  !
  !           AFTERONETEST, the subroutine provided by user to process data after one test
  !           AFTRECORD,  the subroutine provided by user after main-loop for test
  !
  implicit none

  !--- dummy variables and subroutines
   optional::AFTONESTEP,AFTONETEST,AFTRECORD, BEFONETEST, PRERECORD, RECORDPROC, RECSTATU, RESTARTREC, SAVERECSTATU
   external::AFTONESTEP,AFTONETEST,AFTRECORD, BEFONETEST, PRERECORD, RECORDPROC, RECSTATU, RESTARTREC, SAVERECSTATU

  !--- interface to the external routine -------------------
   interface
       subroutine AFTONESTEP(ITime, SimBox, CtrlParam, Statu)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       use MD_TYPEDEF_RecordStamp
       implicit none
       integer                      :: ITime
       type(SimMDBox), dimension(:) :: SimBox
       type(SimMDCtrl)              :: CtrlParam
       integer                      :: Statu
       end subroutine AFTONESTEP
   end interface

   interface
       SUBROUTINE AFTONETEST(SimBox0, Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)                :: SimBox0
       type(MDRecordStamp)           :: Stamp
       type(SimMDBox), dimension(:)  :: SimBox
       type(SimMDCtrl)               :: CtrlParam
       END SUBROUTINE AFTONETEST
   end interface

   interface
       SUBROUTINE AFTRECORD(SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)   :: SimBox
       type(SimMDCtrl)  :: CtrlParam
       END SUBROUTINE AFTRECORD
   end interface

   interface
       SUBROUTINE BEFONETEST(JBOX, SimBox0, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)              :: SimBox0
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       integer::JBOX
       END SUBROUTINE BEFONETEST
   end interface

   interface
       SUBROUTINE PRERECORD(SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)    ::SimBox
       type(SimMDCtrl)   ::CtrlParam
       END SUBROUTINE PRERECORD
   end interface

   interface
       SUBROUTINE RECORDPROC(Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       END SUBROUTINE RECORDPROC
   end interface

   interface
       SUBROUTINE RECSTATU(STATU)
       implicit none
       integer::STATU
       END SUBROUTINE RECSTATU
   end interface

   interface
       SUBROUTINE RESTARTREC(Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       END SUBROUTINE RESTARTREC
   end interface

   interface
       SUBROUTINE SAVERECSTATU(Stamp, SimBox, CtrlParam)
       use MD_TYPEDEF_SimMDBox
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       END SUBROUTINE SAVERECSTATU
   end interface

  !--- END INTERFACE --------------------------------

  !--- Local variables
   procedure(AFTERONESTEP),       pointer::pAftonestep
   procedure(AFTERONETEST),       pointer::pAftonetest
   procedure(AFTERRECORD),        pointer::pAftrecord
   procedure(BEFOREONETEST),      pointer::pBefonetest
   procedure(PRERECORD),          pointer::pPrerecord
   procedure(RECORDPROC),         pointer::pRecordproc
   procedure(RECORDSTATU),        pointer::pRecordStatu
   procedure(RESTARTREC),         pointer::pRestartRec
   procedure(SAVERECSTATU),       pointer::pSaveRecStatu


          if(present(AFTONESTEP)) then
             pAftonestep=>AFTONESTEP
          else
             pAftonestep=>null()
          end if

          if(present(AFTONETEST)) then
             pAftonetest=>AFTONETEST
          else
             pAftonetest=>null()
          end if

          if(present(AFTRECORD)) then
             pAftrecord=>AFTRECORD
          else
             pAftrecord=>null()
          end if

          if(present(BEFONETEST)) then
             pBefonetest=>BEFONETEST
          else
             pBefonetest=>null()
          end if

          if(present(PRERECORD)) then
             pPrerecord=>PRERECORD
          else
             pPrerecord=>null()
          end if

          if(present(RECORDPROC)) then
             pRecordproc=>RECORDPROC
          else
             pRecordproc=>null()
          end if

          if(present(RECSTATU)) then
             pRecordStatu=>RECSTATU
          else
             pRecordStatu=>null()
          end if

          if(present(RESTARTREC)) then
             pRestartRec=>RESTARTREC
          else
             pRestartRec=>null()
          end if

          if(present(SAVERECSTATU)) then
             pSaveRecStatu=>SAVERECSTATU
          else
             pSaveRecStatu=>null()
          end if

          call Add_RecordProcedures(m_Recordlist, pAftonestep, pAftonetest, pAftrecord, pBefonetest, pPrerecord, pRecordproc, pRecordStatu, pRestartRec, pSaveRecStatu)
          return
      end subroutine APPSHELL_AddRecord
  !****************************************************************************************
  !****************************************************************************************
  end module MD_SimboxArray_AppShell_16_GPU
