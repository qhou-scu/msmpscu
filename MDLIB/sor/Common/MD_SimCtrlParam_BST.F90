  module MD_SimMDCtrl_BST
  !***  DESCRIPTION: this module is to load control parameters for Activation-Relaxation
  !                  Technique (ART)
  !                  ______________________________________________________
  !     HISORY:      2018-09-03, first version, Hou Qing
  !
  !
  use MD_CONSTANTS
  use MiniUtilities
  use MD_TYPEDEF_SimMDCtrl
  implicit none

     private::Load_BSTCtl_SimMDCtrl, &
              Load_SectionParameter
  contains
  !*********************************************************************

  !****************************************************************************
  subroutine Load_BSTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(SimMDCtrl)    ::CtrlParam
     character*(*)      ::STR
     integer            ::LINE
     !--- local variables
      character*32::STRNUMB(4),KEYWORD
      integer::I, N, TR
     !----

           !*** start the time step control controlling parametgers
           TR = 0
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)" MDPSCU Warning: unknown keyword in &ARTSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case("&TERMINATE")
                          !$$*** To get total integral steps in this tiemsection
                          call Extract_Numb(STR,2,n,STRNUMB)
                          CtrlParam%ITE = ISTR(STRNUMB(1))
                          if( CtrlParam%ITE.lt.0) then
                              write(*,fmt="(' MDPSCU Error: number of time step smaller than 0, terminal time (in ps) is needed')")
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage:  &TERMINATE  max number of timestep=  ,max time in ps= ')")
                              write(*,fmt="(' Process to be stopped')")
                              stop
                          end if
                          CtrlParam%IT1 = (CtrlParam%IT0-1)+CtrlParam%ITE

                     case("&STEPSIZE")
                          !$$*** To get the length of time step (in fs)
                          call Extract_Numb(STR,4,n,STRNUMB)
                          if(N .LT. 4) then
                             write(*,fmt="(' MDPSCU Error: 4 control parameters for time step size are expected')")
                             write(*,fmt="('               check control file at line:', BZI6)") LINE
                             write(*,fmt="(' Usage:  &STEPSIZE  steps, hmi, hmx, dmx')")
                             write(*,fmt="(' Process to be stopped')")
                             stop
                          end if
                          CtrlParam%IHDUP = ISTR(STRNUMB(1))
                          CtrlParam%HMI = DRSTR(STRNUMB(2))
                          CtrlParam%HMX = DRSTR(STRNUMB(3))
                          CtrlParam%DMX = DRSTR(STRNUMB(4))
                          if(CtrlParam%DMX .LT. 1.D-6) CtrlParam%DMX = 1.D-6
                          !$$The time step start from its minimum value
                          CtrlParam%H =  CtrlParam%HMI

                   case ("&EXTRECORD")
                         !$$--- To get for how many time steps to save current statu
                         call Extract_Numb(STR,1,n,STRNUMB)
                         if(N .LT. 1) then
                            write(*,fmt="(' MDPSCU Warning: value controling external recording subprocess is expected')")
                            write(*,fmt="('                 check control file at line:', BZI6)") LINE
                            write(*,fmt="(' Usage:  &EXTRECORD number of timesteps')")
                            call ONWARNING(gm_OnWarning)
                         else
                            TR = ISTR(STRNUMB(1))
                         end if

                     case ("&ACTMETH")
                          !$$*** To get method for initial activating
                           call Extract_Substr(STR,2,N,STRNUMB)
                           if(N .ge. 1) then
                              if(len_trim(STRNUMB(1)) .gt. 0) then
                                 STRNUMB(1) = adjustl(STRNUMB(1))
                                 call UpCase(STRNUMB(1))
                                 select case(STRNUMB(1))
                                        case("USERPROC")
                                             CtrlParam%ART_ActMeth = CP_USERDEFACT_ART
                                        case("BYARREGION")
                                             CtrlParam%ART_ActMeth = ior(CP_AR_REGION_ART, CP_MULTPLEACT_ART)
                                        case("BYSINGLESEED")
                                             CtrlParam%ART_ActMeth = CP_CENTPARTACT_ART
                                        case("BYMULTSEED")
                                             CtrlParam%ART_ActMeth = ior(CP_CENTPARTACT_ART, CP_MULTPLEACT_ART)
                                        case default
                                            write(*,fmt="(A)")      " MDPSCU Warning:'"//STRNUMB(1)(1:len_trim(STRNUMB(1)))//"' is unknown method for construct active region"
                                            write(*,fmt="(A)")      '                available methods include "USERPROC", "BYSINGLESEED", "BYMULTSEED" '
                                            write(*,fmt="(A, BZI6)")'                check control file at line:',LINE
                                            call ONWARNING(gm_OnWarning)
                                end select
                              end if
                           end if
                           !$$*** To get type of atoms to be used as seed
                           if(iand(CtrlParam%ART_ActMeth, CP_CENTPARTACT_ART) .eq. CP_CENTPARTACT_ART) then
                              call Extract_Numb(STR,min(mp_MXGROUP, size(STRNUMB)),N,STRNUMB)
                              if(N .gt. 0) then
                                 CtrlParam%ART_ActSeed = 0
                                 do I=1, N
                                    CtrlParam%ART_ActSeed(ISTR(STRNUMB(I))) = 1
                                  end do
                              end if
                           end if 

                     case ("&ACTSEED")
                          !$$*** To get type of atoms to be used as seed
                           call Extract_Numb(STR,min(mp_MXGROUP, size(STRNUMB)),N,STRNUMB)
                           CtrlParam%ART_ActSeed = 0
                           do I=1, N
                              CtrlParam%ART_ActSeed(ISTR(STRNUMB(I))) = 1
                           end do

                     case ("&ACTPERCENT")
                          !$$*** To percentage of atoms to be activated
                           call Extract_Numb(STR,1, N,STRNUMB)
                           CtrlParam%ART_ActPecent = DRSTR(STRNUMB(1))

                     case ("&ALPHA")
                          !$$*** To percentage of atoms to be activated
                           call Extract_Numb(STR,1, N,STRNUMB)
                           if(N .ge. 1) then
                              CtrlParam%ART_Alpha = DRSTR(STRNUMB(1))
                              if(CtrlParam%ART_Alpha .le. 0.D0) then
                                  write(*,fmt="(A)")      " MDPSCU Error: the ALPHA parameter is equal or less than zero"
                                 write(*,fmt="(A, BZI6)") '               check control file at line:',LINE
                                 write(*,fmt="(A)")       "               the process to be stopped"
                                 stop
                              end if
                           end if

                     case ("&DISPSTEP")
                          !$$*** To percentage of atoms to be activated
                           call Extract_Numb(STR,1, N,STRNUMB)
                           if(N .ge. 1) then
                              CtrlParam%ART_StepLen = DRSTR(STRNUMB(1))
                              if(CtrlParam%ART_StepLen .le. 0.D0) then
                                 write(*,fmt="(A)")      " MDPSCU Error: the STEPLEN parameter is equal or less than zero"
                                 write(*,fmt="(A, BZI6)") '               check control file at line:',LINE
                                 write(*,fmt="(A)")       "               the process to be stopped"
                                 stop
                              end if
                           end if 

                     case ("&DISPCON")
                          !$$*** To percentage of atoms to be activated
                           call Extract_Numb(STR,1, N,STRNUMB)
                           if(N .ge. 1) then
                              CtrlParam%ART_StepTol = DRSTR(STRNUMB(1))
                              if(CtrlParam%ART_StepTol .le. 0.D0) then
                                  write(*,fmt="(A)")      " MDPSCU Error: the STEPTOL parameter is equal or less than zero"
                                 write(*,fmt="(A, BZI6)") '               check control file at line:',LINE
                                 write(*,fmt="(A)")       "               the process to be stopped"
                                 stop
                              end if
                           end if 

                     case ("&EXTEND")
                          !$$*** To get extend for constructing activation region
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%ART_ActExt = DRSTR(STRNUMB(1))

                     case ("&EVENTDETECTSUBCTL")
                          call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

              end select
           end do

           !*** to check the consistent of the inputs
           if(iand(CtrlParam%ART_ActMeth, CP_CENTPARTACT_ART) .eq. CP_CENTPARTACT_ART .and. &
              all(CtrlParam%ART_ActSeed .eq. 0) ) then
                 write(*,fmt="(A)") " MDPSCU Error: the activation methd of ART is set as by seed(s), "
                 write(*,fmt="(A)") "               but no seeds are given"
                 write(*,fmt="(A)") " Process to be stopped"
                 stop
           end if 

           if(iand(CtrlParam%ART_ActMeth, CP_CENTPARTACT_ART) .eq. CP_CENTPARTACT_ART .or. &
              iand(CtrlParam%ART_ActMeth, CP_AR_REGION_ART)   .eq. CP_AR_REGION_ART) then
              if(CtrlParam%ART_ActPecent .le. 0) then
                 write(*,fmt="(A)") " MDPSCU Error: the pecentage of activating atoms is zero,"
                 write(*,fmt="(A)") "                no activation to be performed"
                 write(*,fmt="(A)") " Process to be stopped"
                 stop
              end if
           end if 

           return
  !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading time step control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
 end subroutine Load_BSTCtl_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  subroutine Load_SectionParameter(hFile, CtrlParam, SimBox, STR, LINE, ISECT)
  !***  PURPOSE:   to load the control parameters for a time section
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer, intent(in)::       hFile

     type(SimMDCtrl)    :: CtrlParam
     type(SimMDBox)     :: SimBox
     character*(*)      :: STR
     integer            :: LINE, ISECT
     !--- local variables
      character*32::KEYWORD
      integer::N
     !----

          !**** to start load the controal parameters
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case ("&ENDSUBCTL")
                          exit

                     case default
                          write(*,*)"MDPSCU warning: unknown keyword in &SECTSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                       case ("&BOUNDSUBCTL")
                             call Load_BoundCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&NEIGHBSUBCTL")
                             call Load_NeighbCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, SimBox%NGROUP, LINE)

                       case ("&ACTIVEREGSUBCTL")
                             call Load_ActiveRegionCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case("&LBFGSSUBCTL")
                             call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case("&STEEPESTSUBCTL", "&CGSUBCTL")
                            call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&BSTSUBCTL")
                             call Load_BSTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&BOOSTSUBCTL")
                             call Load_BoostCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&EVENTDETECTSUBCTL")
                             call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                end select
              end do

              !$$***
               if(CtrlParam%AR_UPTABMI .le. 0) then
                  CtrlParam%AR_UPTABMI = CtrlParam%NB_UPTABMI
                  CtrlParam%AR_UPTABMX = CtrlParam%NB_UPTABMI
                  CtrlParam%AR_DBITAB  = CtrlParam%NB_DBITAB
                  CtrlParam%AR_UPTAB   = CtrlParam%AR_UPTABMI
               end if

  200    return
  !---------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading section parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDCTLF exist to end the data section"
         stop
  end subroutine Load_SectionParameter
 !****************************************************************************

 !****************************************************************************
  subroutine Load_Parameter_SimMDCtrl(hFile, CtrlParam, SimBox, LINE0)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer,        intent(in)  :: hFile
     type(SimMDCtrl),target      :: CtrlParam
     type(SimMDBox)              :: SimBox
     integer,optional,intent(in) :: LINE0
     !--- local variables
      character*256::STR
      character*32::KEYWORD
      integer::LINE, ISECT
      type(SimMDCtrl), pointer::curCtrlParam, next=>null(), nextp=>null()
     !----

          if(present(LINE0)) then
             LINE = LINE0
          else
             LINE = 0
          end if

          call GetInputStrLine(hFile,STR, LINE, "!", *100)
          STR = adjustl(STR)
          call GetKeyWord("&", STR, KEYWORD)
          call UpCase(KEYWORD)
          !**** to start load the controal parameters
             curCtrlParam => CtrlParam
             ISECT = 0
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case ("&ENDCTLF")
                          exit

                     case default
                          write(*,*)"MDPSCU warning: unknown keyword ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&COMMSUBCTL")
                            call Load_CommonParameter_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                            !*** to synchronize the common data in thelist
                             nextp=>CtrlParam
                             do while(.TRUE.)
                                call GetNext_SimMDCtrl(nextp, 1, next)
                                if(.not. associated(next)) exit
                                next%MULTIBOX  = CtrlParam%MULTIBOX
                                next%TOTALBOX  = CtrlParam%TOTALBOX
                                next%INDEPBOX  = CtrlParam%INDEPBOX
                                next%DIFFORDER = CtrlParam%DIFFORDER
                                nextp=>next
                             end do

                       case ("&POTENSUBCTL", "&POTSUBCTL")
                             call Load_ForceCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, SimBox%NGROUP, LINE)
                            !*** to synchronize the potential data in the list
                             nextp=>CtrlParam
                             do while(.TRUE.)
                                call GetNext_SimMDCtrl(nextp, 1, next)
                                if(.not. associated(next)) exit
                                next%RU         = CtrlParam%RU
                                next%NUMFTABR   = CtrlParam%NUMFTABR
                                next%NUMFTABE   = CtrlParam%NUMFTABE
                                next%RHOSCAL    = CtrlParam%RHOSCAL
                                next%OUTPUTFTAB = CtrlParam%OUTPUTFTAB
                                nextp=>next
                             end do

                      case("&LBFGSSUBCTL")
                             call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                             !*** to synchronize the LBFGS data in thelist
                             nextp=>CtrlParam
                             do while(.TRUE.)
                                call GetNext_SimMDCtrl(nextp, 1, next)
                                if(.not. associated(next)) exit

                                next%LBFGS_PGtol     = CtrlParam%LBFGS_PGtol
                                next%LBFGS_Factr     = CtrlParam%LBFGS_Factr
                                next%LBFGS_Msave     = CtrlParam%LBFGS_Msave
                                next%NEB_LBFGS_PGtol = CtrlParam%NEB_LBFGS_PGtol
                                next%NEB_LBFGS_Factr = CtrlParam%NEB_LBFGS_Factr
                                next%NEB_LBFGS_Msave = CtrlParam%NEB_LBFGS_Msave

                                nextp=>next
                             end do

                      case("&STEEPESTSUBCTL", "&CGSUBCTL")
                             call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                             !*** to synchronize the STEEPEST data in the list
                             nextp=>CtrlParam
                             do while(.TRUE.)
                                call GetNext_SimMDCtrl(nextp, 1, next)
                                if(.not. associated(next)) exit
 
                                next%STEEPEST_Alpha  = CtrlParam%STEEPEST_Alpha
                                next%STEEPEST_MiStep = CtrlParam%STEEPEST_MiStep
                                next%STEEPEST_MxStep = CtrlParam%STEEPEST_MxStep
                                next%STEEPEST_MiDelE = CtrlParam%STEEPEST_MiDelE
                                nextp=>next
                             end do

                       case ("&ANALYSUBCTL")
                             call Load_AnalyCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                            !*** to synchronize the common data in the list
                             nextp=>CtrlParam
                             do while(.TRUE.)
                                call GetNext_SimMDCtrl(nextp, 1, next)
                                if(.not. associated(next)) exit
                                next%JOBID0    = CtrlParam%JOBID0
                                next%JOBID1    = CtrlParam%JOBID1
                                next%JOBIDSTEP = CtrlParam%JOBIDSTEP
                                next%STARTCFG  = CtrlParam%STARTCFG
                                next%ENDCFG    = CtrlParam%ENDCFG
                                next%CFGSTEP   = CtrlParam%CFGSTEP
                                next%STARTBOX  = CtrlParam%STARTBOX
                                next%ENDBOX    = CtrlParam%ENDBOX
                                next%BOXSTEP   = CtrlParam%BOXSTEP

                                next%NEEDPOT       = CtrlParam%NEEDPOT
                                next%NEEDDAMP      = CtrlParam%NEEDDAMP
                                next%NEEDDAMPTYPE  = CtrlParam%NEEDDAMPTYPE
                                next%NEEDNEIGHB    = CtrlParam%NEEDNEIGHB
                                next%NEEDNEWATOMS  = CtrlParam%NEEDNEWATOMS

                                nextp=>next
                             end do

                       case ("&SECTSUBCTL")
                             if(ISECT .gt. 0) then
                                allocate(next)
                                next%parent =>null()
                                next%next   =>null()
                                call Copy_SimMDCtrl(curCtrlParam, next)
                                next%IT0 =  curCtrlParam%IT1+1
                                next%IT1 =  next%IT0+next%ITE-1
                                next%DAMPTIME0 = next%IT0
                                next%IVTIME0   = next%IT0
                                next%EPCTIME0  = next%IT0
                                !$$--- by default, we set none thermalization for time section
                                !$$    this is different from the default set for PARREP-MD
                                !$$    NOTE: 2018-03-06
                                next%IVTIME     = 0
                                next%DAMPTIME1  = 0
                                next%EPCTIME1   = 0

                                curCtrlParam=>next
                                call AddNext_SimMDCtrl(CtrlParam, curCtrlParam)
                             end if
                             ISECT = ISECT+1
                             call Load_SectionParameter(hFile, curCtrlParam, SimBox, STR,  LINE, ISECT)

                end select
              end do
              !--- make stopping power model to be be the same for all section
        return
  !---------------------------------------------------------------
  100    write(*,fmt="(A, I5)")"MDPSCU Error in reading control parameters at line", LINE
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDCTLF exist to end the data section"
         stop
  end subroutine Load_Parameter_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Print_Parameter_SimMDCtrl(hFile, CtrlParam, SimBox)
  !***  PURPOSE:   to print out control parameters for dnamics applications
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer, intent(in)     ::hFile
     type(SimMDCtrl),target ::CtrlParam
     type(SimMDBox)         ::SimBox

     !--- local variables
     integer::I, J, nsect, IS, CENTPART(mp_MXGROUP)
     type(SimMDCtrl), pointer::sectCtrlParam, next
     !----

        !$$**** HEADER
          write(hFile,fmt="(' !************ SIMULATION CONTROL PARAMETERS *******************')")

              if(CtrlParam%RESTART.le.0) then
                 write(hFile,FMT="(' !   A new simulation to be performed  ')")
              else
                 write(hFile,FMT="(' !   The simulation to be restarted from last stop point', 3(1PE11.4,1x))")
              end if
              write(hFile,fmt="(' !')")

              !*************************
              write(hFile,FMT="(' !    Total number of boxes to be considered........: ', 3(BZI7,1x))") CtrlParam%TotalBox
              write(hFile,FMT="(' !    Number of boxes in one test ..................: ', 3(BZI7,1x))") CtrlParam%MULTIBOX
              write(hFile,FMT="(' !    If boxes are independent .....................: ', 3(BZI7,1x))") CtrlParam%INDEPBOX
              write(hFile,fmt="(' !')")


              !*************************
              write(hFile,FMT="(' !    For force table calculation...................: ', 20(BNI8,1x))")
              write(hFile,FMT="(' !       Number of points in force table............: ', 20(BNI8,1x))") CtrlParam%NUMFTABR
              write(hFile,FMT="(' !       Number of points in embed-fun table........: ', 20(BNI8,1x))") CtrlParam%NUMFTABE
              !write(hFile,FMT="(' !    RHO-scal in embedment function table..........: ',  (1PE11.4,1x))")CtrlParam%RHOSCAL

              write(hFile,FMT="(' !       Cutoff range of interactions (in LU).......: ', 20(F8.2,1x))")
              do I=1, SimBox%NGROUP
                 do J=I, SimBox%NGROUP
                    write(hFile,FMT="(' !           between group',I2,' and group',I2,' ...........:', 20(F8.2,1x))")I, J, CtrlParam%RU(I,J)
                 end do
              end do
              write(hFile,fmt="(' !')")

              !*************************
              if(CtrlParam%SEED(1) .gt. 0) then
                 write(hFile,FMT="(' !    Random number seed ...........................: ',4(I8,1x))") CtrlParam%SEED(1)
              else
                 write(hFile,FMT="(' !    Random number seed determined automatically  ',4(I8,1x))")
              end if
              write(hFile,fmt="(' !')")

              !************************
              call Number_SimMDCtrl(CtrlParam, nsect)
              write(hFile,FMT="(' !    Number of time sections.......................: ', 4(I4,1x))") nsect
              write(hFile,fmt="(' !')")

              !*************************
              sectCtrlParam=>CtrlParam
              do IS=1, nsect
                 write(hFile,FMT="(' !    For time section # ', 4(I4,1x))") IS
                 write(hFile,FMT="(' !      The section will start at time steps........: ', 3(BZI8,1x))")   sectCtrlParam%IT0
                 write(hFile,FMT="(' !      and end at time steps.......................: ', 3(BZI8,1x))")   sectCtrlParam%IT1
                 write(hFile,FMT="(' !      Time step(fs) is ............................: ', 3(F10.2,1x))") sectCtrlParam%HMI
                 write(hFile,fmt="(' !')")

                 write(hFile,FMT="(' !      If periodic conditions used (X,Y,Z).........: ', 3(I3,1x))")   sectCtrlParam%IFPD(1:3)
                 write(hFile,fmt="(' !')")

                 !--- the input for neighbor-list calculation
                 write(hFile,FMT="(' !      Max permitted neighbors.....................: ', 3(BNI4,1x))") sectCtrlParam%NB_MXNBS
                 write(hFile,FMT="(' !      Min step number to update neighbor list.....: ', 3(BNI4,1x))") sectCtrlParam%NB_UPTABMI
                 write(hFile,FMT="(' !      Max step number to update neighbor list.....: ', 3(BNI4,1x))") sectCtrlParam%NB_UPTABMX
                 write(hFile,FMT="(' !      Step number to change updating intervale....: ', 3(BNI6,1x))") sectCtrlParam%NB_DBITAB
                 write(hFile,FMT="(' !      Cutoff range of neighbor region(in RCUT)....: ', 20(F8.2,1x))")
                 do I=1, SimBox%NGROUP
                    do J=I, SimBox%NGROUP
                        write(hFile,FMT="(' !          between group',I2,' and group',I2,'.............:', 20(F8.2,1x))")I, J, sectCtrlParam%NB_RM(I,J)
                    end do
                 end do
                 write(hFile,fmt="(' !')")

                 !--- the input for activation region control
                 CENTPART = 0
                 I = 0
                 do J=1, mp_MXGROUP
                    if(sectCtrlParam%AR_CENTPART(J) .gt. 0) then
                       I = I + 1
                       CENTPART(I) = J
                     end if
                 end do

                 if(iand(sectCtrlParam%AR_METHOD,CP_ENABLE_AR) .ne. CP_ENABLE_AR) then
                    write(hFile,FMT="(' !      Activation region calculation ..............: disable', 3(BNI4,1x))")
                 else
                    write(hFile,FMT="(' !      Activation region calculation ..............: enable', 3(BNI4,1x))")
                    write(hFile,FMT="(' !      Min step number to update active region.....: ', 3(BNI4,1x))") sectCtrlParam%AR_UPTABMI
                    write(hFile,FMT="(' !      Max step number to update active region.....: ', 3(BNI4,1x))") sectCtrlParam%AR_UPTABMX
                    write(hFile,FMT="(' !      Step number to change updating intervale....: ', 3(BNI6,1x))") sectCtrlParam%AR_DBITAB

                    if(iand(sectCtrlParam%AR_METHOD, CP_USERDEF_AR) .eq. CP_USERDEF_AR ) then
                       write(hFile,FMT="(' !      User supplied procedure to be used for AR calculation', 3(BNI6,1x))")
                    else
                       if(iand(sectCtrlParam%AR_METHOD, CP_CENTPART_AR) .eq. CP_CENTPART_AR .and.  &
                          iand(sectCtrlParam%AR_METHOD, CP_EKIN_AR) .eq. CP_EKIN_AR) then
                          write(hFile,FMT="(' !      AR build by seeds of atom type..............:', 20(BNI6,1x))")      &
                                           (CENTPART(J), J=1, count(CENTPART.gt.0))
                          write(hFile,FMT="(' !           and by kinetic energy..................:', F8.3,1x, '(eV)'))") &
                                            sectCtrlParam%AR_EKIN
                       else if( iand(sectCtrlParam%AR_METHOD, CP_CENTPART_AR).eq.CP_CENTPART_AR) then
                          write(hFile,FMT="(' !      AR build by seeds of atom type..............:', 20(BNI6,1x))")      &
                                           (CENTPART(J), J=1, count(CENTPART.gt.0))
                       else
                          write(hFile,FMT="(' !      AR build by seeds of atoms of energy above..:', F8.3,1x, '(eV)'))") &
                                            sectCtrlParam%AR_EKIN
                       end if

                       if(ibits(sectCtrlParam%AR_METHOD, CP_BYNBSETBIT_AR,1) .gt. 0) then
                          write(hFile,FMT="(' !      with an extend..............................: ', 1(BNI6,1x), &
                                            ' by neighborhood')") sectCtrlParam%AR_EXTEND
                       else
                           write(hFile,FMT="(' !      with an extend..............................: ', 1(BNI6,1x), &
                                             ' by linked-cells')") sectCtrlParam%AR_EXTEND
                       end if

                       if(iand(sectCtrlParam%AR_METHOD,CP_KEEP_AR) .eq. CP_KEEP_AR) then
                           write(hFile,FMT="(' !      keeping activation statu of atoms after disabling AR')")
                       else 
                           write(hFile,FMT="(' !      not keeping statu of atoms after disabling AR')")
                       end if    
                    end if
                 end if
                 write(hFile,fmt="(' !')")

                 !--- the input for activation region control
                 write(hFile,FMT="(' !      Control parameter for ART ..................: ', 3(BNI4,1x))")
                 write(hFile,FMT="(' !          Alpha ..................................:', F8.3)") sectCtrlParam%ART_Alpha
                 write(hFile,FMT="(' !          Step length (LU)........................:', F8.3)") sectCtrlParam%ART_StepLen
                 write(hFile,FMT="(' !          Step tolerance .........................:', F8.3)") sectCtrlParam%ART_StepTol
                 write(hFile,FMT="(' !          Displ threshold (LU) for event check....: ', 20(F8.2,1x))") sectCtrlParam%STRCUT_DRTol
                 write(hFile,FMT="(' !          Initial activation method...............:',  20(BNI6,1x))") 
                 CENTPART = 0
                 I = 0
                 do J=1, mp_MXGROUP
                    if(sectCtrlParam%ART_ActSeed(J) .gt. 0) then
                       I = I + 1
                       CENTPART(I) = J
                     end if
                 end do

                 if(iand(sectCtrlParam%ART_ActMeth, CP_CENTPARTACT_ART) .eq. CP_CENTPARTACT_ART) then
                    if(dabs(sectCtrlParam%ART_ActExt) .le. 1.D-8) then
                       write(hFile,FMT="(' !          By atoms of types.......................:',  20(BNI6,1x))") &
                                                (CENTPART(J), J=1, count(CENTPART.gt.0))
                    else if(sectCtrlParam%ART_ActExt .gt. 0) then
                        write(hFile,FMT="(' !          By atoms of types.......................:',  20(BNI6,1x))") &
                                                  (CENTPART(J), J=1, count(CENTPART.gt.0))
                        write(hFile,FMT="(' !             and their neighbores in extent(LU)...:',  F8.3)") &
                                                  sectCtrlParam%ART_ActExt
                    else
                        write(hFile,FMT="(' !          By neighbors of atoms of types..........:',  20(BNI6,1x))") &
                                                  (CENTPART(J), J=1, count(CENTPART.gt.0))
                    end if
                    if(iand(CtrlParam%ART_ActMeth, CP_MULTPLEACT_ART) .eq. CP_MULTPLEACT_ART) then
                        write(hFile,FMT="(' !          Percentage of atoms to be activated.....:', F8.3)") &
                                                  sectCtrlParam%ART_ActPecent
                    else
                        write(hFile,FMT="(' !          One of the atoms to be activated........:')")
                    end if

                 else if(sectCtrlParam%ART_ActMeth .eq. CP_AR_REGION_ART) then
                    write(hFile,FMT="(' !          By atoms in activation region...........:',  20(BNI6,1x))")
                    if(iand(CtrlParam%ART_ActMeth, CP_MULTPLEACT_ART) .eq. CP_MULTPLEACT_ART) then
                        write(hFile,FMT="(' !          Percentage of atoms to be activated.....:', F8.3)") &
                                                  sectCtrlParam%ART_ActPecent
                    else
                        write(hFile,FMT="(' !          One of the atoms to be activated........:')")
                    end if

                 else if(sectCtrlParam%ART_ActMeth .eq. CP_USERDEFACT_ART) then
                    write(hFile,FMT="(' !          By atoms in activation region...........:',  20(BNI6,1x))")
                                                                                 
                 end if
                 write(hFile,fmt="(' !')")
                 !--- to to next section
                 call GetNext_SimMDCtrl(sectCtrlParam, 1, next)
                 sectCtrlParam=>next
           end do ! end do for sections

         return
  end subroutine Print_Parameter_SimMDCtrl
  !****************************************************************************

  end module MD_SimMDCtrl_BST
