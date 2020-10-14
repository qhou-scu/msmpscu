  module MD_SimMDCtrl_PARREP
  !***  DESCRIPTION: this module is to load control parameters for Parallel-Replica MD
  !
  !                  ______________________________________________________
  !     HISORY:      2016-11-16, seperated the from MD_TypeDef_SimMDCtrl.F90
  !
  use MD_CONSTANTS
  use MiniUtilities
  use MD_TYPEDEF_SimMDCtrl
  implicit none

     private::Load_TemperaturCtl,      &
              Load_TimestepCtl,        &
              Load_SectionParameter,   &
              Load_PARREPCtl
  contains
  !*********************************************************************

  !****************************************************************************
  subroutine Load_TemperaturCtl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,        intent(in)::hFile
     type(SimMDCtrl)           ::CtrlParam
     character*(*)             ::STR
     integer                   ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

             !*** start the temprature controlling parametgers
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case( "&ENDSUBCTL")
                                exit
                         case default
                              write(*,*)" MDPSCU Warning: unknown keyword in &TEMPSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)

                         case("&TEMPERATURE")
                              !*** the temperature
                               call Extract_Numb(STR,1,n,STRNUMB)
                               CtrlParam%TI = DRSTR(STRNUMB(1))

                         case("&THERMALIZATION","&THERMAL")
                               !$$*** To get how many times of thermalization will be carried out
                                 call Extract_Numb(STR,2,n,STRNUMB)
                                 if(n .lt. 2) then
                                    write(*,fmt="(' MDPSCU Error: the time steps for thermalization should be set')")
                                    write(*,fmt="('               check control file at line:', BZI6)") LINE
                                    write(*,fmt="(' Usage:  &THERMALIZATION numbe of thermalizing=  , interval for thermalizing =, method=')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if
                                 CtrlParam%IVTIME  = ISTR(STRNUMB(1))
                                 CtrlParam%IVPAS   = ISTR(STRNUMB(2))
                                 if(CtrlParam%IVTIME .gt. 0 .and.  CtrlParam%IVPAS .eq. 0) then
                                    write(*,fmt="(' MDPSCU Error: the interval for thermalizing cannot be zero')")
                                    write(*,fmt="('               check control file at line:', BZI6)") LINE
                                    write(*,fmt="(' Usage:  &THERMALIZATION numbe of thermalizing=  , interval for thermalizing =, method=')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if

                                 call Extract_Substr(STR,1,N,STRTMP)
                                 if(n .ge. 1) then
                                    call UpCase(STRTMP(1))
                                    if(STRTMP(1)(1:len_trim(STRTMP(1))) .ne. "MC") then
                                        write(*,fmt="(A)")      " MDPSCU Warning: for PARREP MD, only MC method is used for depahsing"
                                        write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                        call ONWARNING(gm_OnWarning)
                                    end if
                                 end if

                         case("&EPC", "&EPCSUBCTRL")
                                !$$*** To get if electron-phono coupling is required
                                 write(*,fmt="(' MDPSCU Warning: for PARREP MD, E-P coupling is not supported')") LINE
                                 write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                 call ONWARNING(gm_OnWarning)

                         case("&QUICKDUMP", "&QUICKDAMP", "&QUENCH")
                              write(*,fmt="(A, A, A)")  " MDPSCU Warning: unknown keyword ", KEYWORD(1:LEN_TRIM(KEYWORD)),  " in &TIMESUBCTL section"
                              write(*,fmt="(A, A)")  "                 this keyword could appear in &PARREPSUBCTL subsection"
                              write(*,fmt="(A,BZI6)")"                 check control file at line:", LINE
                              call ONWARNING(gm_OnWarning)

                  end select
              end do
             !*** end of tmeprature control section
             return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading temperature control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_TemperaturCtl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_TimestepCtl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer,          intent(in)::hFile
     type(SimMDCtrl)             ::CtrlParam
     character*(*)               ::STR
     integer                     ::LINE
     !--- local variables
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----
           !*** start the time step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)" MDPSCU Warning: unknown keyword in &TIMESUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case("&TERMINATE")
                          !$$*** To get total integral steps in this tiemsection
                          call Extract_Numb(STR,2,n,STRNUMB)
                          CtrlParam%ITE = ISTR(STRNUMB(1))
                          if(n.lt. 2) then
                            if( CtrlParam%ITE.lt.0) then
                                write(*,fmt="(' MDPSCU Error: number of time step smaller than 0, terminal time (in ps) is needed')")
                                write(*,fmt="('               check control file at line:', BZI6)") LINE
                                write(*,fmt="(' Usage:  &TERMINATE  max number of timestep=  ,max time in ps= ')")
                                write(*,fmt="(' Process to be stopped')")
                                stop
                            end if
                            CtrlParam%TEMTM = 1.D60  !to give an large value
                          else
                            CtrlParam%TEMTM = DRSTR(STRNUMB(2))
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

                    case ("&OUTPUT_T")
                         !$$--- To get for how many time steps to print out temperature etc
                          call Extract_Numb(STR,1,n,STRNUMB)
                          if(N .LT. 1) then
                             write(*,fmt="(' MDPSCU Warning: 1 control parameter controling output of thermal quantities is expected')")
                             write(*,fmt="('               check control file at line:', BZI6)") LINE
                             write(*,fmt="(' Usage:  &OUTPUT_T  number of timesteps')")
                             call ONWARNING(gm_OnWarning)
                           else
                             CtrlParam%TimestpQ = ISTR(STRNUMB(1))
                           end if

                   case ("&EXTRECORD")
                         !$$--- To get for how many time steps to save current statu
                         call Extract_Numb(STR,2,n,STRNUMB)
                         if(N .LT. 1) then
                            write(*,fmt="(' MDPSCU Warning: value controling external recording subprocess is expected')")
                            write(*,fmt="('                 check control file at line:', BZI6)") LINE
                            write(*,fmt="(' Usage:  &EXTRECORD number of timesteps')")
                            call ONWARNING(gm_OnWarning)
                         else
                            CtrlParam%TimestpR = ISTR(STRNUMB(1))
                         end if

                   case ("&OUTPUT_C")
                        write(*,fmt="(A, A, A)")  " MDPSCU Warning: unknown keyword ", KEYWORD(1:LEN_TRIM(KEYWORD)), " in &TIMESUBCTL section"
                        write(*,fmt="(A, A)")     "                 this keyword is replaced by &EVENTCHECKSTEPS in &PARREPSUBCTL subsection"
                        write(*,fmt="(A,BZI6)")   "                 check control file at line:", LINE
                        call ONWARNING(gm_OnWarning)

              end select
           end do
           !*** end the time step control controlling parametgers
           return
  !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading time step control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
 end subroutine Load_TimestepCtl
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
     integer,intent(in):: hFile
     type(SimMDCtrl)   :: CtrlParam
     type(SimMDBox)    :: SimBox
     character*(*)     :: STR
     integer           :: LINE, ISECT
     !--- local variables
      character*32::KEYWORD
      integer::I, N
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

                       case ("&TEMPSUBCTL")
                            !*** start the temprature controlling parametgers
                            call Load_TemperaturCtl(hFile, CtrlParam, STR, LINE)

                       case ("&PRESSSUBCTL")
                            !*** start the press step control controlling parametgers
                            call Load_PressureCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&TIMESUBCTL")
                             call Load_TimestepCtl(hFile, CtrlParam, STR, LINE)

                       case ("&BOUNDSUBCTL")
                             call Load_BoundCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&NEIGHBSUBCTL")
                             call Load_NeighbCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, SimBox%NGROUP, LINE)

                       case ("&ACTIVEREGSUBCTL")
                             call Load_ActiveRegionCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&PARREPSUBCTL")
                              call Load_PARREPCtl(hFile, CtrlParam, STR, LINE)
                              
                       case ("&EVENTDETECTSUBCTL")
                             call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                end select
              end do

              !$$*** to check the consistent of the control parameters
              if(CtrlParam%PARREP_FineR .le. 0) CtrlParam%PARREP_FineR = max(CtrlParam%TimestpG, 100)/100
              if(CtrlParam%TimestpSave  .le. 0) CtrlParam%TimestpSave = CtrlParam%ITE

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
  subroutine Load_PARREPCtl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for events detection
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

           !*** start the time step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)

              !---
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)" MDPSCU Warning: unknown keyword in &PARREPSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&FINECHECKSTEPS")
                          !$$--- fine time steps of saving intermediate configuration
                          call Extract_Numb(STR,2,n,STRNUMB)
                          if(N .LT. 1) then
                             write(*,fmt="(' MDPSCU Warning: parameter controling refining event checking is missed')")
                             write(*,fmt="('                 check control file at line:', BZI6)") LINE
                             write(*,fmt="(' Usage:  &FINEREC number of timesteps')")
                             call ONWARNING(gm_OnWarning)
                          else
                             CtrlParam%PARREP_FineR = ISTR(STRNUMB(1))
                          end if

                     case ("&EVENTCHECKSTEPS")
                          !$$--- To get for how many time steps to check the structure changes
                          call Extract_Numb(STR,1,n,STRNUMB)
                          if(N .LT. 1) then
                             write(*,fmt="(' MDPSCU Warning: parameter controling event checking is missed')")
                             write(*,fmt="('                 check control file at line:', BZI6)") LINE
                             write(*,fmt="(' Usage:  &EVENTCHECKSTEPS number of timesteps')")
                             call ONWARNING(gm_OnWarning)
                           else
                             CtrlParam%PARREP_Timestep = ISTR(STRNUMB(1))
                          end if

                     case ("&TERMINATEVENT")
                         !$$--- To get if the time-section will exist at an event detected
                         call Extract_Numb(STR,1,n,STRNUMB)
                         if(N .LT. 1) then
                            write(*,fmt="(' MDPSCU Warning: value for controling parameter is missed')")
                            write(*,fmt="('                 check control file at line:', BZI6)") LINE
                            write(*,fmt="(' Usage:  &EXITATEVENT value, the value couble be 0 or 1')")
                            call ONWARNING(gm_OnWarning)
                         else
                            CtrlParam%ExitAtEvent = ISTR(STRNUMB(1))
                         end if

                     case ("&EVENTDETECTSUBCTL")
                         call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

              end select
           end do
           !*** end the time step control controlling parametgers
           return
  !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading parallel-replica control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
 end subroutine Load_PARREPCtl
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
     integer,                intent(in):: hFile
     type(SimMDCtrl),target            :: CtrlParam
     type(SimMDBox)                    :: SimBox
     integer,optional,       intent(in):: LINE0
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
                            !*** to synchronize the common data in the list
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

                       case ("&POTENSUBCTL","&POTSUBCTL")
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
                             !*** to synchronize the common data in thelist
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

                                next%NEEDPOT   = CtrlParam%NEEDPOT
                                next%NEEDDAMP  = CtrlParam%NEEDDAMP
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

                                curCtrlParam=>next
                                call AddNext_SimMDCtrl(CtrlParam, curCtrlParam)
                             end if
                             ISECT = ISECT+1
                             call Load_SectionParameter(hFile, curCtrlParam, SimBox, STR,  LINE, ISECT)

                end select
              end do

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
     integer,                intent(in)::hFile
     type(SimMDCtrl), target           ::CtrlParam
     type(SimMDBox)                    ::SimBox

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
                 write(hFile,FMT="(' !      or  terminated at time(ps)..................: ', 3(1PE11.4,1x))")sectCtrlParam%TEMTM
                 write(hFile,fmt="(' !')")

                 if(sectCtrlParam%DAMPTIME1  .LE. 0  .and. &
                    sectCtrlParam%IVTIME     .LE. 0  .and. &
                    sectCtrlParam%EPCTIME1   .LE. 0  .and. &
                    sectCtrlParam%IFNOSE     .LE. 0         ) then
                    write(hFile,FMT="(' !      Temperature control.........................:  NONE ', 3(BZF8.2,1x))")
                 else
                    write(hFile,FMT="(' !      Temperature(K)..............................: ', 3(BZF8.2,1x))")   sectCtrlParam%TI
                    if(sectCtrlParam%DAMPTIME1 .GT. 0) then
                       write(hFile,FMT="(' !       QuickDumping to be performed in section....: ', BZI5)") IS
                       write(hFile,FMT="(' !         start from the timestep .................: ', BZI8)") sectCtrlParam%DAMPTIME0
                       write(hFile,FMT="(' !         with the number of time steps............: ', BZI8)") sectCtrlParam%DAMPTIME1
                       select case( sectCtrlParam%DAMPSCHEME)
                              case (CP_DAMPSCHEME_LBFGS)
                                   write(hFile,FMT="(' !         damping scheme...........................: LNFGS method ', BZI8)")
                                   write(hFile,FMT="(' !                 the convergence criteria.........: PGTOL ', 3(1PE11.4,1x))") sectCtrlParam%LBFGS_PGtol
                                   write(hFile,FMT="(' !                                                    FACTR ', 3(1PE11.4,1x))") sectCtrlParam%LBFGS_Factr
                                   write(hFile,FMT="(' !                                                    MSAVE ', BZI8)")          sectCtrlParam%LBFGS_MSave

                              case (CP_DAMPSCHEME_DYN)
                                   write(hFile,FMT="(' !         damping scheme...........................: DYNAMICS method ', BZI8)")
                       end select
                    end if

                    if(sectCtrlParam%IVTIME .GT. 0) then
                       write(hFile,FMT="(' !       Thermalization to be performed in section..: ', BZI5)") IS
                       write(hFile,FMT="(' !          start from the time step ...............: ', BZI8)") sectCtrlParam%IVTIME0
                       write(hFile,FMT="(' !          num of thermalizing in the section......: ', BZI8)") sectCtrlParam%IVTIME
                       write(hFile,FMT="(' !          num of time steps between thermalizing..: ', BZI8)") sectCtrlParam%IVPAS
                       select case(sectCtrlParam%IVSCHEME)
                              case (CP_THERMALSCHEME_MC)
                                write(hFile,FMT="(' !          thermalizating scheme...................: MC method ', BZI8)")
                              case (CP_THERMALSCHEME_VSCAL)
                                write(hFile,FMT="(' !          thermalizating scheme...................: V-scaling method ', BZI8)")
                              case (CP_THERMALSCHEME_PSCAL)
                                write(hFile,FMT="(' !          thermalizating scheme...................: P-scaling method ', BZI8)")
                              case (CP_THERMALSCHEME_ASCAL)
                                write(hFile,FMT="(' !          thermalizating scheme...................: A-scaling method ', BZI8)")
                       end select

                    end if

                    if(sectCtrlParam%EPCTIME1 .GT. 0) then
                       write(hFile,FMT="(' !       E-P couipling to be performed in section..: ', BZI5)") IS
                       write(hFile,FMT="(' !         start from the timestep .................: ', BZI8)") sectCtrlParam%EPCTIME0
                       write(hFile,FMT="(' !         with the number of time steps............: ', BZI8)") sectCtrlParam%EPCTIME1
                       !write(hFile,FMT="(' !         with physical model defined in ..........: ', A)") sectCtrlParam%EPCPHYS(1:len_trim(CtrlParam%EPCPHYS))
                    end if
                 end if
                 write(hFile,fmt="(' !')")

                 !*************************
                 if(sectCtrlParam%IFBPC .le. 0 .and. &
                    sectCtrlParam%IFJONSHON .le. 0) then
                       write(hFile,FMT="(' !      Pressure control ...........................:  NONE ', 3(BZF8.2,1x))")
                 else
                       write(hFile,FMT="(' !      Pressure(kbar)..............................: ', 3(F8.2,1x))")   sectCtrlParam%PEX
                 end if
                  write(hFile,fmt="(' !')")

                 if(sectCtrlParam%IHDUP .eq. 0) then
                    write(hFile,FMT="(' !      Fixed time step size to be used..............: ', 3(F10.2,1x))")
                    write(hFile,FMT="(' !      Time step(fs) is ............................: ', 3(F10.2,1x))") sectCtrlParam%HMI
                 else
                    if(sectCtrlParam%IHDUP .gt. 0) then
                       write(hFile,FMT="(' !      Varibale time step (I) to be used ..........: ', 3(I7,1x))")
                       write(hFile,FMT="(' !      Number of steps of increasing time step.....: ', 3(I7,1x))")    IABS(sectCtrlParam%IHDUP)
                       write(hFile,FMT="(' !      Min time step size(fs)......................: ', 3(F10.2,1x))") sectCtrlParam%HMI
                       write(hFile,FMT="(' !      Max time step size(fs)......................: ', 3(F10.2,1x))") sectCtrlParam%HMX
                    else
                       write(hFile,FMT="(' !      Varibale time step (II) to be used .........: ', 3(I7,1x))")
                       write(hFile,FMT="(' !      Number of steps of adjusting time step......: ', 3(I7,1x))")    IABS(sectCtrlParam%IHDUP)
                       write(hFile,FMT="(' !      Max time step size(fs)......................: ', 3(F10.2,1x))") sectCtrlParam%HMX
                       write(hFile,FMT="(' !      Max displacement of atoms(A)................: ', 3(F10.2,1x))") sectCtrlParam%DMX
                    endif
                 end if
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
                 I        = 0
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
                          write(hFile,FMT="(' !      AR build by kinetic energy..................:', F8.3,1x, '(eV)'))") &
                                            sectCtrlParam%AR_EKIN
                       end if

                       if(ibits(sectCtrlParam%AR_METHOD, CP_BYNBSETBIT_AR,1) .gt. 0) then
                          write(hFile,FMT="(' !      with an extend..............................: ', 1(BNI6,1x),  &
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

                 !--- the input for output control
                 write(hFile,FMT="(' !      Number of timesteps of output of thermal quantities: ', 20(BNI8,1x))") sectCtrlParam%TimestpQ
                 write(hFile,FMT="(' !      Number of timesteps of invoking external recording.: ', 20(BNI8,1x))") sectCtrlParam%TimestpR
                 write(hFile,fmt="(' !')")

                 if(sectCtrlParam%ExitAtEvent .gt. 0) then
                 write(hFile,FMT="(' !      Exit proecess at event detected....................:  YES', 20(BNI8,1x))")
                 else
                 write(hFile,FMT="(' !      Exit proecess at event detected....................:  NO', 20(BNI8,1x))")
                 end if
                 write(hFile,FMT="(' !      Number of timesteps for detecting event............: ', 20(BNI8,1x))") sectCtrlParam%TimestpG
                 write(hFile,FMT="(' !      Number of timesteps of refining event detection....: ', 20(BNI8,1x))") sectCtrlParam%PARREP_FineR
                 write(hFile,FMT="(' !      Displacement threshold (LU) for strcuture change....: ', 20(F8.2,1x))")sectCtrlParam%STRCUT_DRTol
                 write(hFile,fmt="(' !')")


                 call GetNext_SimMDCtrl(sectCtrlParam, 1, next)
                 sectCtrlParam=>next
           end do ! end do for sections

         return
  end subroutine Print_Parameter_SimMDCtrl
  !****************************************************************************

  end module MD_SimMDCtrl_PARREP
