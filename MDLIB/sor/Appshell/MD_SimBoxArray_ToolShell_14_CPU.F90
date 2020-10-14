  module MD_SimBoxArray_ToolShell_14_CPU
  !***  DESCRIPTION:
  !     This module provide a common interface to extract and analysis the
  !     box configure data.
  !
  !
  !    Written by HOU Qing, Mar, 2011
  use MD_Globle_Variables
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_RecordList
  use MiniUtilities

  implicit none

      type(SimMDBox), dimension(:), allocatable,private::m_SimBox
      type(RecordProcedureList),                private::m_Recordlist

      private::For_One_Test_GMD,     &
               For_One_Box_NEB,      &
               For_One_Test_ParRep
  contains

  !****************************************************************************************
  subroutine Main_ANALYSIS(processid, FORCETABLE)
  !***  DESCRIPTION: the main process running in "batch" mode
  !     INPUT: processid, the id of the process, usable only when MPI used,
  !            STARTP, ENDP, STEP, optional integers to determine the time interval in which the data to be processed
  !    NOTE: the three dummy subroutine must strictly follow the interface given below
  use MD_TYPEDEF_PrintList, only:Do_PrintProcess
  implicit none

  !--- dummy variables and subroutines
  integer,intent(in)::processid

  !--- interface to the external routine -------------------
   OPTIONAL::FORCETABLE
   external::FORCETABLE
   interface
     subroutine FORCETABLE(SimBox, CtrlParam, FTable, printout)
     use MD_CONSTANTS
     use MD_TYPEDEF_SimMDBox
     use MD_TYPEDEF_SimMDCtrl
     use MD_TYPEDEF_ForceTable
     implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(inout)::FTable
      integer,           optional     ::printout
     end subroutine FORCETABLE
   end interface
  !--- END INTERFACE --------------------------------

  !--- Local variable
      integer::J,K, K1, JB, TSTEP, next, ITEST, JLOOP

      type(SimMDBox)      ::SimBox0
      type(SimMDCtrl)     ::CtrlParam, CtrlParam0
      type(MDRecordStamp) ::Stamp
      character*32        ::ARGV


      !****
           call OpenLogFile_Globle_Variables( )

      !*** to read in control parameters from IniFile
   1000   call Release_SimMDCtrl(CtrlParam)
          call Initialize_Globle_Variables(SimBox0, CtrlParam, NEXT=next)
          if(next .LE. 0) then
             write(ARGV, *) iabs(next)
             ARGV = adjustl(ARGV)
             write(*,*) "Process stop after finishing job "//ARGV(1:len_trim(ARGV))
             return
          end if
          call Copy_SimMDCtrl(CtrlParam ,CtrlParam0)

      !*** to initialize the configuration
          call Initialize_SimMDBox(SimBox0,CtrlParam%DiffOrder)

      !*** if configure files contains more atoms than in initial file,
      !    atoms needed to be added to the box. The number of atoms to
      !    be added is give in ANALYTIC SUB-CONTROL in control file
          call AddAtoms_Globle_Variables(SimBox0, CtrlParam)

      !**** if user provide preprocessor, we do it
          call DoPrerecord_List(m_Recordlist, SimBox0, CtrlParam)
          call CommandlineBoxsel(CtrlParam)

      !**** check the consistent of input for analysis
           if(CtrlParam0%TOTALBOX/CtrlParam0%MULTIBOX .lt. CtrlParam%JOBID0) then
             write(*,fmt="(A)") ' MDPSCU Warning: the number of TESTs is smaller than  JOBI0'
             write(*,fmt="(A)") '                 no analysis to be performed'
             write(*,fmt="(A)") '                 check &JOBSEL input in analysis control file'
             write(*,fmt="(A)") '                 or parameter for option -JOB on command line.'
             call ONWARNING(gm_OnWarning)
          end if

          if(CtrlParam0%TOTALBOX/CtrlParam0%MULTIBOX .lt. CtrlParam%JOBID1) then
             write(*,fmt="(A)") ' MDPSCU Warning: the number of TESTs is smaller than  JOBI1'
             write(*,fmt="(A)") '                 some analysis to be missed'
             write(*,fmt="(A)") '                 check &JOBSEL input in analysis control file.'
             write(*,fmt="(A)") '                 or parameter for option -JOB on command line.'
             call ONWARNING(gm_OnWarning)
          end if

          if(CtrlParam0%MULTIBOX .lt. CtrlParam%STARTBOX) then
             write(*,fmt="(A)") ' MDPSCU Warning: the number of boxes in one TESTs is smaller than STARTBOX in analysis'
             write(*,fmt="(A)") '                 check &BOXSEL input in analysis control file'
             write(*,fmt="(A)") '                 or parameter for option -BOX on command line.'
             call ONWARNING(gm_OnWarning)
          end if

          if(CtrlParam0%MULTIBOX .lt. CtrlParam%ENDBOX) then
             write(*,fmt="(A)") ' MDPSCU Warning: the number of boxes in one TESTs is smaller than ENDBOX in analysis'
             write(*,fmt="(A)") '                 check &BOXSEL input in analysis control file'
             write(*,fmt="(A)") '                 or parameter for option -BOX on command line.'
             call ONWARNING(gm_OnWarning)
          end if

      !****
          call Print_AnalyPara_SimMDCtrl(6, SimBox0, CtrlParam)
          call Do_PrintProcess(6, SimBox0, CtrlParam)
          if(gm_hFILELOG.gt.0) then
             call Print_AnalyPara_SimMDCtrl(gm_hFILELOG, SimBox0, CtrlParam)
             call Do_PrintProcess(gm_hFILELOG, SimBox0, CtrlParam)
           end if

      !**** TO BEGINE LOOP ON SAMPLES
           select case(gm_AppType(1:len_trim(gm_AppType)))
                  case default
                      if(CtrlParam%TIMELOOPOUT .gt. 0) then
                         JLOOP = 1
                         ITEST = 1
                      else
                        JLOOP = CtrlParam%TOTALBOX/CtrlParam%MULTIBOX
                        ITEST = CtrlParam%JOBID0
                      end if
                  case ("GMD")
                      if(CtrlParam%TIMELOOPOUT .gt. 0) then
                         JLOOP = 1
                         ITEST = 1
                      else
                        JLOOP = CtrlParam%TOTALBOX/CtrlParam%MULTIBOX
                        ITEST = CtrlParam%JOBID0
                      end if
                  case ("PARREP")
                       !*** NOTE:   For PARREP, one test results only one trajecory,
                       !            and the recording time may be different for diffrent
                       !            trajectories.
                       CtrlParam%TIMELOOPOUT = 0
                       JLOOP = CtrlParam%TOTALBOX/CtrlParam%MULTIBOX
                       ITEST = CtrlParam%JOBID0
                  case ("NEB")
                       !*** NOTE:   For NEB, diffrent box may have different number of images,
                       !            we can only to analysis box-by-box.
                       CtrlParam%TIMELOOPOUT = 0
                       JLOOP = CtrlParam%TOTALBOX
                       ITEST = CtrlParam%JOBID0
           end select

           do J=1, JLOOP
              !--- process only required jobs
              if(J .NE. ITEST) cycle

              !--- do some preprocess before one test
              call DoBeforeOneTest_List(m_Recordlist, J, SimBox0, m_SimBox, CtrlParam)

              !---
              select case(gm_AppType(1:len_trim(gm_AppType)) )
                  case default
                       call For_One_Test_GMD(SimBox0, CtrlParam, J,processid)
                  case ("GMD")
                       call For_One_Test_GMD(SimBox0, CtrlParam, J,processid)
                  case ("PARREP")
                        call For_One_Test_PARREP(SimBox0, CtrlParam, J,processid)
                  case ("NEB")
                        call For_One_Box_NEB(SimBox0, CtrlParam, J,processid)
              end select


              !--- to postprocess the data for this test
              call Default_RecordStamp(Stamp)
              Stamp%ITest   =  J
              call DoAfterOneTest_List(m_Recordlist, SimBox0, Stamp, m_SimBox, CtrlParam)

              !--- update ITETS to be analysis
              ITEST = ITEST + CtrlParam%JOBIDSTEP
              if(CtrlParam%JOBID1 .gt. 0) then
                 if(ITEST .gt. CtrlParam%JOBID1) exit
              end if
           end do  !end the loop for samples
           call Release_SimBoxArray(m_SimBox)

           !**** if external after-process provided, we do it
           call DoAfterRecord_List(m_Recordlist, SimBox0, CtrlParam)

          !**** to performat next simulation
          goto 1000

    return
  end subroutine Main_ANALYSIS
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Test_GMD(SimBox0,CtrlParam,J,processid)
  !***  DESCRIPTION: to load and process of data.
  !     INPUT: SimBox0,     the element box
  !            CtrlParam,   the control parameters for simulation
  !            J,           the Jth sample
  !            processid,   the id of the proces
  !    NOTE: the  dummy subroutine RECORDPROC must strictly follow the interface given below


  implicit none
   integer,              intent(in)::J,processid
   type(SimMDBox)                  ::SimBox0
   type(SimMDCtrl), target         ::CtrlParam
   type(MDRecordStamp)             ::Stamp

  !--- Local variables
   integer::ITIME, ISTART, IEND, ISTEP, ICFG, ISECT, ERRMSG, NCFGS, MULTIBOX0, ITEST, JJ, J1, NJJ
  !----

              !*** Now we can create the box ARRAY
              MULTIBOX0 = CtrlParam%MULTIBOX
              if(CtrlParam%TIMELOOPOUT .gt. 0) then
                 if(CtrlParam%JOBID1 .ge. CtrlParam%JOBID0) then
                    CtrlParam%MULTIBOX = (CtrlParam%JOBID1 - CtrlParam%JOBID0 + 1)*MULTIBOX0
                 else
                    CtrlParam%MULTIBOX = (CtrlParam%TOTALBOX/MULTIBOX0 - CtrlParam%JOBID0 + 1)*MULTIBOX0
                 end if
                 J1 = CtrlParam%JOBID0
              else
                 CtrlParam%MULTIBOX = MULTIBOX0
                 J1 = J
              end if

              NJJ = CtrlParam%MULTIBOX/MULTIBOX0
              call Create_SimBoxArray(m_SimBox, CtrlParam%MULTIBOX,SimBox0)

              !*** determine the data we need to process
              ISTART = CtrlParam%STARTCFG
              ISTEP = CtrlParam%CFGSTEP
              IEND = CtrlParam%ENDCFG
              if(IEND .lt. ISTART .and. ISTEP .gt. 0) ISTEP = -ISTEP
              if(IEND .lt. 0) IEND = 0

              !*** check if we have these configures
              call NumberCfg_SimMDCtrl(CtrlParam, NCFGS)
              if(NCFGS .lt. ISTART) then
                 write(*,fmt="(A)") ' MDPSCU Warning: the number of configures is smaller than  STARTCFG'
                 write(*,fmt="(A)") '                 no analysis to be performed'
                 write(*,fmt="(A)") '                 check the number total number of times steps in the control file'
                 call ONWARNING(gm_OnWarning)
                 ISTART = NCFGS
              end if

              !*** to initialize the processor by set ITIME = -1
              STAMP%ITest = J
              STAMP%ITime = -1
              STAMP%Time  = -1
              call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)

              !*** start loop for time sections
              STAMP%ITIME    = -1
              STAMP%TIME     = -1
              STAMP%SCALTIME = -1
              do ICFG = ISTART, IEND, ISTEP
                 call GetCfgITime_SimMDCtrl(CtrlParam, ICFG, ITIME)
                 call GetSectID_SimMDCtrl(CtrlParam, ITIME, ISECT)
                 ITEST = 1
                 do JJ=1, NJJ
                    STAMP%ITest   = J1+JJ-1
                    STAMP%IBox(1) = (STAMP%ITest-1)* MULTIBOX0 + 1
                    STAMP%IBox(2) =  STAMP%IBox(1) + MULTIBOX0 - 1
                    STAMP%ITime   = ITIME
                    STAMP%ISect   = ISECT
                    STAMP%ICfg    = ICFG
                    STAMP%IRec    = STAMP%ICfg
                    call Putin_Instance_Config_SimBoxArray(CtrlParam, m_SimBox(ITEST:ITEST+MULTIBOX0-1), STAMP, ERR=ERRMSG)

                    if(ERRMSG .EQ. CP_FILE_NOTEXIST) then
                       exit
                    end if
                    ITEST = ITEST + MULTIBOX0
                 end do
                 if(ERRMSG .EQ. CP_FILE_NOTEXIST)exit

                 !--- to processs the information of box
                 STAMP%ITest = J
                 call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)
              end do   !end the loop for time step

           !call Release_SimBoxArray(SimBox)
    return
  end subroutine For_One_Test_GMD
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Test_ParRep(SimBox0,CtrlParam,J,processid)
  !***  DESCRIPTION: to load and process of data.
  !     INPUT: SimBox0,     the element box
  !            CtrlParam,   the control parameters for simulation
  !            J,           the Jth sample
  !            processid,   the id of the proces
  !    NOTE: the  dummy subroutine RECORDPROC must strictly follow the interface given below
  implicit none
   integer,                intent(in)::J,processid
   type(SimMDBox)                    ::SimBox0
   type(SimMDCtrl), target           ::CtrlParam
   type(MDRecordStamp)               ::Stamp

  !--- Local variables
   integer::ISTART, IEND, ISTEP, ICFG, ERRMSG, MULTIBOX0
   character*256::GFile
  !----

              !*** NOTE:   For ParRep, we cannot have TIMELOOPOUT, because the
              !            the time steps for the same cfg-ID could be different.
              !            for each time step, we will load only the 0K configurations:
              !
              MULTIBOX0 = CtrlParam%MULTIBOX
              CtrlParam%MULTIBOX = 1
              call Create_SimBoxArray(m_SimBox, CtrlParam%MULTIBOX,SimBox0)

              !*** determine the data we need to process
              ISTART = CtrlParam%STARTCFG
              ISTEP  = CtrlParam%CFGSTEP
              IEND   = CtrlParam%ENDCFG
              if(IEND .lt. ISTART .and. ISTEP .gt. 0) ISTEP = -ISTEP
              if(IEND .lt. 0) IEND = 0

              !*** to initialize the processor by set ITIME = -1
              STAMP%ITest =  J
              STAMP%ITime = -1
              STAMP%Time  = -1
              call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)

              !*** start loop for time sections
              STAMP%ITIME    = -1
              STAMP%TIME     = -1
              STAMP%SCALTIME = -1
              do ICFG = ISTART, IEND, ISTEP
                 STAMP%ITest   = J
                 STAMP%IBox(1) = STAMP%ITest
                 STAMP%IBox(2) = STAMP%IBox(1)
                 STAMP%ITime   = -1
                 STAMP%ISect   = -1
                 STAMP%ICfg    = ICFG
                 STAMP%IRec    = STAMP%ICfg
                 call STRCATI(GFILE, CtrlParam%f_geometry, "_0K_P", processid, 4)
                 call STRCATI(GFILE, GFILE, "_",  J, 4)
                 call STRCATI(GFile, GFile, ".", ICFG, 4)
                 call Putin_Instance_Config_SimBoxArray(GFILE, m_SimBox(1:1), Stamp, ERR=ERRMSG)

                 if(ERRMSG .EQ. CP_FILE_NOTEXIST) then
                    exit
                 end if

                 !--- to processs the information of box
                 call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)
              end do   !end the loop for time step
              CtrlParam%MULTIBOX = MULTIBOX0
    return
  end subroutine For_One_Test_ParRep
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Box_NEB(SimBox0,CtrlParam,J,processid)
  !***  DESCRIPTION: to load and process of data.
  !     INPUT: SimBox0,     the element box
  !            CtrlParam,   the control parameters for simulation
  !            J,           the Jth sample
  !            processid,   the id of the proces
  !    NOTE: the  dummy subroutine RECORDPROC must strictly follow the interface given below
  implicit none
   integer,                 intent(in) ::J,processid
   type(SimMDBox),          intent(in) ::SimBox0
   type(SimMDCtrl),target              ::CtrlParam

  !--- Local variables
   integer::ISTART, IEND, ISTEP, ICFG, ERRMSG, MULTIBOX0
   character*256::GFILE0, GFILE
   type(MDRecordStamp)::STAMP
  !----

              !*** NOTE:   For NEB, diffrent box may have different number of images,
              !            we can only to analysis box-by-box.
              MULTIBOX0 = CtrlParam%MULTIBOX
              CtrlParam%MULTIBOX = 1
              call Create_SimBoxArray(m_SimBox, CtrlParam%MULTIBOX,SimBox0)

              !*** determine the data we need to process
              ISTART = CtrlParam%STARTCFG
              ISTEP  = CtrlParam%CFGSTEP
              IEND   = CtrlParam%ENDCFG
              if(IEND .lt. ISTART .and. ISTEP .gt. 0) ISTEP = -ISTEP
              if(IEND .lt. 0) IEND = 0
              !*** to initialize the processor by set ITIME = -1
              STAMP%ITest =  J
              STAMP%ITime = -1
              STAMP%Time  = -1
              call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)

              !*** to prepare file name of NEBImg
              call GetPath(CtrlParam%f_geometry, GFILE0)
              GFILE0 = GFILE0(1:len_trim(GFILE0))//"NEBImg"

              !*** start loop for time sections
              STAMP%ITIME    = -1
              STAMP%TIME     = -1
              STAMP%SCALTIME = -1
              do ICFG = ISTART, IEND, ISTEP
                 !*** to prepare checking stamp
                 STAMP%ITest   = J
                 STAMP%IBox(1) = STAMP%ITest
                 STAMP%IBox(2) = STAMP%IBox(1)
                 STAMP%ITime   = -1
                 STAMP%ISect   = -1
                 STAMP%ICfg    = ICFG
                 STAMP%IRec    = STAMP%ICfg

                !*** to prepare the filename for configuration
                 call STRCATI(GFILE, GFILE0, "_P",  0,  4)
                 call STRCATI(GFILE, GFILE,  "_",   STAMP%IBox(1), 4)
                 call STRCATI(GFile, GFILE,  ".",   ICFG, 4)
                 call Putin_Instance_Config_SimBoxArray(GFILE, m_SimBox(1:1), STAMP, ERR=ERRMSG)

                 if(ERRMSG .EQ. CP_FILE_NOTEXIST)exit

                 !--- to processs the information of box
                 call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)
              end do !*** end the loop for time section
              CtrlParam%MULTIBOX = MULTIBOX0

    return
  end subroutine For_One_Box_NEB
  !****************************************************************************************

  !****************************************************************************************
  subroutine APPSHELL_AddRecord(AFTONETEST,AFTRECORD, BEFONETEST, PRERECORD, RECORDPROC, RECSTATU, RESTARTREC, SAVERECSTATU)
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
   optional::AFTONETEST,AFTRECORD, BEFONETEST, PRERECORD, RECORDPROC, RECSTATU, RESTARTREC, SAVERECSTATU
   external::AFTONETEST,AFTRECORD, BEFONETEST, PRERECORD, RECORDPROC, RECSTATU, RESTARTREC, SAVERECSTATU

  !--- interface to the external routine -------------------
   interface
       SUBROUTINE AFTONETEST(SimBox0, Stamp, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)              :: SimBox0
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       END SUBROUTINE AFTONETEST
   end interface

   interface
       SUBROUTINE AFTRECORD(SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)   :: SimBox
       type(SimMDCtrl)  :: CtrlParam
       END SUBROUTINE AFTRECORD
   end interface

   interface
       SUBROUTINE BEFONETEST(JBOX, SimBox0, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_SimboxArray
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
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox)   ::SimBox
       type(SimMDCtrl)  ::CtrlParam
       END SUBROUTINE PRERECORD
   end interface

   interface
       SUBROUTINE RECORDPROC(Stamp, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_SimboxArray
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
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       END SUBROUTINE RESTARTREC
   end interface

   interface
       SUBROUTINE SAVERECSTATU(Stamp, SimBox, CtrlParam)
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(MDRecordStamp)         :: Stamp
       type(SimMDBox), dimension(:):: SimBox
       type(SimMDCtrl)             :: CtrlParam
       END SUBROUTINE SAVERECSTATU
   end interface

  !--- END INTERFACE --------------------------------

  !--- Local variables
   procedure(AFTERONESTEP),       pointer::pAftonestep=>null()
   procedure(AFTERONETEST),       pointer::pAftonetest
   procedure(AFTERRECORD),        pointer::pAftrecord
   procedure(BEFOREONETEST),      pointer::pBefonetest
   procedure(PRERECORD),          pointer::pPrerecord
   procedure(RECORDPROC),         pointer::pRecordproc
   procedure(RECORDSTATU),        pointer::pRecordStatu
   procedure(RESTARTREC),         pointer::pRestartRec
   procedure(SAVERECSTATU),       pointer::pSaveRecStatu


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
  subroutine Single_ANALYSIS(processid,extCmdlinePara)
  !***  DESCRIPTION: the main process running in "single file" mode
  !     INPUT: processid, the id of the process, usable only when MPI used,
  implicit none

  !--- dummy variables and subroutines
  integer,intent(in)::processid
  optional::extCmdlinePara
  external::extCmdlinePara
  !--- interface to the external routine -------------------
   interface
       SUBROUTINE extCmdlinePara()
       END SUBROUTINE extCmdlinePara
   end interface

  !--- Local variable
      character*256::cmdline=""
      character*32::tstr(5)
      integer::I, NB
      type(SimMDCtrl)::CtrlParam
      type(MDRecordStamp)::STAMP


           call ExtractCmdlineFileNames_Globle_Variables()
           if(present(extCmdlinePara)) then
            call extCmdlinePara()
           end if   

           NB = count(len_trim(gm_cfgFileName) .gt. 0)
           allocate(m_SimBox(NB))
           m_SimBox(1:NB)%proAutoLoad = 1
           m_SimBox(1:NB)%NPRT        = 0
           !$$--- to initialize the controal parameters
            CtrlParam%JOBID0    = 1
            CtrlParam%JOBID1    = 1
            CtrlParam%JOBIDSTEP = 1
            CtrlParam%STARTBOX  = 1
            CtrlParam%BOXSTEP   = 1
            CtrlParam%ENDBOX    = 1
            CtrlParam%STARTCFG  = 0
            CtrlParam%ENDCFG    = 0
            CtrlParam%CFGSTEP   = 1
            CtrlParam%MULTIBOX  = 1

           !$$--- to load the configure
            call DoPrerecord_List(m_Recordlist, m_SimBox(1), CtrlParam)
            do I=1, NB
               write(*,*) "Load configuration from ", gm_cfgFileName(I)(1:len_trim(gm_cfgFileName(I)) )
               call Load_Config_SimMDBox(gm_cfgFileName(I), m_SimBox(I))
            end do
            if(any(CtrlParam%RU.lt.0)) then
                CtrlParam%RU  = m_SimBox(1)%RR
            else
                CtrlParam%RU = CtrlParam%RU*m_SimBox(1)%RR
            end if
            CtrlParam%NB_RM  = CtrlParam%RU*CtrlParam%NB_RM

            !--- to processs the information of box
            STAMP%ITest =  0
            STAMP%ITime =  1
            STAMP%Time  =  0.D0
            call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)

            !$$--- clear the allocated memory
            if(allocated(m_SimBox)) then
               call Release_SimBoxArray(m_SimBox)
               deallocate(m_SimBox)
             end if
             call Release_SimMDCtrl(CtrlParam)
    return
  end subroutine Single_ANALYSIS
  !****************************************************************************************

  !****************************************************************************************
  subroutine Multi_ANALYSIS(processid, ExitNofile,extCmdlinePara)
  !***  DESCRIPTION: the main process running in "multi file" mode. In contrast to "single file"
  !                  mode, the boxes data are loaded by given TEST IDs, and CFG IDs. The IDs are
  !                  are given in command line. In contrast to using Main_ANALYSIS, no 'setup file"
  !                  are needed.
  !     INPUT: processid, the id of the process, usable only when MPI used,
  implicit none

  !--- dummy variables and subroutines
  integer, intent(in)          ::processid
  integer, intent(in), optional::ExitNofile
  optional::extCmdlinePara
  external::extCmdlinePara
  !--- interface to the external routine -------------------
   interface
       SUBROUTINE extCmdlinePara()
       END SUBROUTINE extCmdlinePara
   end interface  
  !--- Local variable
      character*256      ::cmdline="", GFILE
      character*32       ::tstr(5)
      integer            ::I, NB, ICFG, ITEST, EXONNOF
      logical            ::EX
      type(SimMDCtrl) ::CtrlParam
      type(MDRecordStamp)::STAMP

           call ExtractCmdlineFileNames_Globle_Variables()
           if(present(extCmdlinePara)) then
            call extCmdlinePara()
           end if   

           call CommandlineBoxsel(CtrlParam)
           !$$--- to allocate memeory for simulation box
            NB = count(len_trim(gm_cfgFileName) .gt. 0)
            allocate(m_SimBox(NB))
            m_SimBox(1:NB)%proAutoLoad = 1

            CtrlParam%MULTIBOX  = 1
            call DoPrerecord_List(m_Recordlist, m_SimBox(1), CtrlParam)

            EXONNOF = 1 
            if(present(ExitNofile)) EXONNOF = ExitNofile

            if(CtrlParam%TIMELOOPOUT) then
               do ICFG = CtrlParam%STARTCFG, CtrlParam%ENDCFG, CtrlParam%CFGSTEP
                  do ITEST = CtrlParam%JOBID0, CtrlParam%JOBID1, CtrlParam%JOBIDSTEP
                     STAMP%ITest =  ITest
                     STAMP%IBox  =  ITest
                     STAMP%ICfg  =  ICFG
                     STAMP%IRec  =  ICFG
                     STAMP%ITime =  1
                     STAMP%Time  =  0.D0
                     m_SimBox(1:NB)%NPRT = 0
                     do I=1, NB
                        call STRCATI(GFILE, gm_cfgFileName(I), "P", processid, 4)
                        call STRCATI(GFILE, GFILE, "_", ITEST, 4)
                        call STRCATI(GFILE, GFILE, ".", ICFG, 4)
                        write(*,*) "Load configuration from ", GFILE(1:len_trim(GFILE) )
                        inquire(FILE=GFILE, EXIST=EX)
                        if(.not. EX) then
                           if(EXONNOF .gt. 0) then
                              write(*,*)  "MDPSCU Error: "//GFILE(1:len_trim(GFILE))//" not exsit"
                              write(*,*)  '               process to be stopped '
                              stop
                           else 
                              write(*,*)  'MDPSCU Warning: '//GFILE(1:len_trim(GFILE))//' not exsit'
                              exit 
                           end if
                        end if      
                        call Load_Config_SimMDBox(GFILE, m_SimBox(I),STAMP)
                     end do
                     if(.not.EX) exit
                     !--- to processs the information of box
                     if(any(CtrlParam%RU.lt.0)) then
                        CtrlParam%RU = m_SimBox(1)%RR
                     else
                        CtrlParam%RU = CtrlParam%RU*m_SimBox(1)%RR
                     end if
                     CtrlParam%NB_RM = CtrlParam%RU*CtrlParam%NB_RM
                     call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)
                  end do
               end do
            else
               do ITEST = CtrlParam%JOBID0, CtrlParam%JOBID1, CtrlParam%JOBIDSTEP
                  do ICFG = CtrlParam%STARTCFG, CtrlParam%ENDCFG, CtrlParam%CFGSTEP
                     !--- to processs the information of box
                     STAMP%ITest =  ITest
                     STAMP%IBox  =  ITest
                     STAMP%ICfg  =  ICFG
                     STAMP%IRec  =  ICFG
                     STAMP%ITime =  1
                     STAMP%Time  =  0.D0
                     m_SimBox(1:NB)%NPRT = 0
                     do I=1, NB
                        call STRCATI(GFILE, gm_cfgFileName(I), "P", processid, 4)
                        call STRCATI(GFILE, GFILE, "_", ITEST, 4)
                        call STRCATI(GFILE, GFILE, ".", ICFG, 4)
                        write(*,*) "Load configuration from ", GFILE(1:len_trim(GFILE) )
                        inquire(FILE=GFILE, EXIST=EX)
                        if(.not. EX) then
                           if(EXONNOF .gt. 0) then
                              write(*,*)  "MDPSCU Error: "//GFILE(1:len_trim(GFILE))//" not exsit"
                              write(*,*)  '               process to be stopped '
                              stop
                           else 
                              write(*,*)  'MDPSCU Warning: '//GFILE(1:len_trim(GFILE))//' not exsit'
                              exit 
                           end if
                        end if      
                        call Load_Config_SimMDBox(GFILE, m_SimBox(I), STAMP)
                     end do
                     if(.not.EX) exit

                     if(any(CtrlParam%RU.lt.0)) then
                        CtrlParam%RU = m_SimBox(1)%RR
                     else
                        CtrlParam%RU = CtrlParam%RU*m_SimBox(1)%RR
                     end if
                     CtrlParam%NB_RM = CtrlParam%RU*CtrlParam%NB_RM
                     call DoRecord_List(m_Recordlist, STAMP, m_SimBox, CtrlParam)
                  end do
               end do
            end if
            call DoAfterRecord_List(m_Recordlist, m_SimBox(1), CtrlParam)

            !$$--- clear the allocated memory
            if(allocated(m_SimBox)) then
               call Release_SimBoxArray(m_SimBox)
               deallocate(m_SimBox)
             end if
             call Release_SimMDCtrl(CtrlParam)
 
    return
  end subroutine Multi_ANALYSIS
  !****************************************************************************************



  end module MD_SimBoxArray_ToolShell_14_CPU
  !****************************************************************************************
