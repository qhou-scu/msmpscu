  module MD_SimMDCtrl_GMD
  !***  DESCRIPTION: this module is to load control parameters for GMD
  !                  ______________________________________________________
  !     HISORY:      2016-11-16, seperated the from MD_TypeDef_SimMDCtrl.F90, Hou Qing
  !
  !                 2018-03-08, add the control mod for local temperature control, HOU Qing
  !
  use MD_CONSTANTS
  use MiniUtilities
  use MD_TYPEDEF_SimMDCtrl
  implicit none

     private::Load_TemperaturCtl_Global,    &
              Load_TemperaturCtl_Local,     &
              Load_TimestepCtl,             &
              Load_SectionParameter
  contains
  !*********************************************************************

  !****************************************************************************
  subroutine Load_TemperaturCtl_Global(hFile, CtrlParam, STR, LINE)
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
      character*32::STRNUMB(10),KEYWORD
      integer::I, J, N
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
                               CtrlParam%TI         = DRSTR(STRNUMB(1))
                               CtrlParam%LT_CTRL%TI = CtrlParam%TI

                         case("&QUICKDUMP", "&QUICKDAMP", "&QUENCH")
                               !*** To get if damping is required
                                call Extract_Numb(STR,1,n,STRNUMB)
                                 if(n .lt. 1) then
                                    write(*,fmt="(' MDPSCU Error: the time steps for quenching should be set')")
                                    write(*,fmt="('               check control file at line:', BZI6)") LINE
                                    write(*,fmt="(' Usage: &QUICKDAMP (&QUENCH) time steps need for damping =')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if
                                 CtrlParam%DAMPTIME1 = ISTR(STRNUMB(1))

                                 call Extract_Substr(STR,1,N,STRNUMB)
                                 if(n .ge. 1) then
                                    call UpCase(STRNUMB(1))
                                    if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "LBFGS") then
                                       CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_LBFGS
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG") then
                                       CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_CG
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG-LS") then
                                       CtrlParam%DAMPSCHEME = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                      
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST") then
                                       CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_ST
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST-LS") then
                                       CtrlParam%DAMPSCHEME = ior(CP_DAMPSCHEME_ST, CP_DAMPSCHEME_LSEARCH)                                     
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYN" .or. &
                                            STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYNAMICS" ) then
                                       CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_DYN
                                    else
                                        write(*,fmt="(A)")      " MDPSCU Error: the damping scheme "//STRNUMB(1)(1:len_trim(STRNUMB(1)))//" is unknown"
                                        write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                        write(*,fmt="(A)")      ' Process to be stopped'
                                        stop
                                    end if
                                 end if
                          case("&STEEPESTSUBCTL", "&CGSUBCTL")
                                call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
 
                          case("&THERMALIZATION", "&THERMAL")
                               !$$*** To get how many times of thermalization will be carried out
                                 call Extract_Numb(STR,2,n,STRNUMB)
                                 if(n .lt. 1) then
                                    write(*,fmt="(' MDPSCU Error: the number of thermalization circle should be set')")
                                    write(*,fmt="('               check control file at line:', BZI6)") LINE
                                    write(*,fmt="(' Usage:  &THERMALIZATION numbe of thermalizing=  , interval for thermalizing =, method=')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if
                                 CtrlParam%IVTIME  = ISTR(STRNUMB(1))

                                 if(CtrlParam%IVTIME .gt. 0 .and.  n .lt. 2) then
                                    write(*,fmt="(' MDPSCU Error: the timestep interval for thermalizing should be set')")
                                    write(*,fmt="('               check control file at line:', BZI6)") LINE
                                    write(*,fmt="(' Usage:  &THERMALIZATION numbe of thermalizing=  , interval for thermalizing =, method=')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if
                                 if(n .ge. 2) CtrlParam%IVPAS   = ISTR(STRNUMB(2))

                                 if(CtrlParam%IVTIME .gt. 0 .and.  CtrlParam%IVPAS .eq. 0) then
                                    write(*,fmt="(' MDPSCU Error: the interval for thermalizing cannot be zero')")
                                    write(*,fmt="('               check control file at line:', BZI6)") LINE
                                    write(*,fmt="(' Usage:  &THERMALIZATION numbe of thermalizing=  , interval for thermalizing =, method=')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if

                                 call Extract_Substr(STR,1,N,STRNUMB)
                                 if(n .ge. 1) then
                                    call UpCase(STRNUMB(1))
                                    if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "MC") then
                                       CtrlParam%IVSCHEME = CP_THERMALSCHEME_MC
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "VSCAL") then
                                       CtrlParam%IVSCHEME = CP_THERMALSCHEME_VSCAL
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "PSCAL") then
                                       CtrlParam%IVSCHEME = CP_THERMALSCHEME_PSCAL
                                    else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ASCAL") then
                                       CtrlParam%IVSCHEME = CP_THERMALSCHEME_ASCAL
                                    else
                                        write(*,fmt="(A)")      " MDPSCU Error: the thermalizing scheme "//STRNUMB(1)(1:len_trim(STRNUMB(1)))//" is unknown"
                                        write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                        write(*,fmt="(A)")      '               supported thermalizing schemes are: "MC", "VSCAL",  "PSCAL" and "ASCAL"'
                                        write(*,fmt="(A)")      ' Process to be stopped'
                                        stop
                                    end if
                                 end if

                          case("&EPC", "&EPCSUBCTRL")
                                !$$*** To get if electron-phono coupling is required
                                 call Extract_Numb(STR, 1,n,STRNUMB)
                                 if(n .lt. 1) then
                                    write(*,fmt="(' MDPSCU Error:  time steps should be given for EPC at line', BZI6)") LINE
                                    write(*,fmt="(' Usage:  &EPC  timesteps, or &EPCSUBCTRL timesteps')")
                                    write(*,fmt="(' where:  timesteps:  timesteps in which E-P coupling calculation considered')")
                                    write(*,fmt="('         filename:   file  name contains the parameters for the E-P calculation at the time')")
                                    write(*,fmt="(' Process to be stopped')")
                                    stop
                                 end if

                                 if(ISTR(STRNUMB(1)) .gt. 0) then
                                    CtrlParam%EPCTIME1 = ISTR(STRNUMB(1))
                                    if(CtrlParam%EPCTIME1 .gt. 0) then
                                       !call GetInputStrLine(hFile,STR, LINE, "!", *100)
                                       !STR = adjustl(STR)
                                       !call GetKeyWord("&", STR, KEYWORD)
                                       !call UpCase(KEYWORD)
                                       !if(KEYWORD(1:LEN_TRIM(KEYWORD)) .ne. "&EPCSUBCTL") then
                                       !   write(*,fmt="(A)")       'MDPSCU Error:  EPC is eneabled,  EPC control data is missed'
                                       !   write(*,fmt="(A, BZI6)") '               add  &EPCSUBCTL should be given at line', LINE
                                       !   write(*,fmt="(A)")       ' Process to be stopped'
                                       !   stop
                                       !else
                                          CtrlParam%LT_CTRL(1)%TI = CtrlParam%TI
                                          CtrlParam%LT_CTRL(1)%METH  = CP_TICTRL_METH_EPC
                                          call Load_EPC_TiCtrlParam(hFile, CtrlParam%LT_CTRL(1), STR, LINE)
                                       !end if
                                       CtrlParam%LT_CTRL = CtrlParam%LT_CTRL(1)
                                    end if
                                    !call Extract_Substr(STR,1,n,STRTMP)
                                    !if(n .lt. 1) then
                                    !   write(*,fmt="(' MDPSCU Error: filename of physical data should be given for E-P at line', BZI6)") LINE
                                    !   write(*,fmt="(' Process to be stopped')")
                                    !   stop
                                    ! else
                                    ! CtrlParam%EPCPHYS = STRTMP(1)
                                    !end if
                                 end if

                          case("&NOSE_PARA")
                               !$$*** To get Nose parameter if Nose thermo-bath used
                                call Extract_Numb(STR,1,n,STRNUMB)
                                CtrlParam%NoseQ = DRSTR(STRNUMB(1))
                  end select
              end do
             !*** end of tmeprature control section
             return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading temperature control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_TemperaturCtl_Global
  !****************************************************************************

  !****************************************************************************
  subroutine Load_TemperaturCtl_Local(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,          intent(in)::hFile
     type(SimMDCtrl)             ::CtrlParam
     character*(*)               ::STR
     integer                     ::LINE
     !--- local variables
      character*32::STRNUMB(mp_mxGROUP),KEYWORD
      integer::I, J, N, IGS(mp_mxGROUP), IGS0(mp_mxGROUP)
      type(TiCtrlParam)::TCtrl

     !----
             !*** start the temprature controlling parametgers
             IGS0 = 0
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case( "&ENDSUBCTL")
                                exit
                         case default
                              write(*,*)" MDPSCU Warning: unknown keyword in &TEMPSUBCTL for local temperature control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('              check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)

                         case("&GROUPSUBCTL","&GRSUBCTL","&GRPSUBCTL","&GSUBCTL")
                              !*** the group IDS
                               IGS = 0
                               call Extract_Numb(STR, mp_mxGROUP,n,STRNUMB)
                               if(N .le. 0) then
                                  write(*,fmt="(' MDPSCU Error: atomic types is missed at line', BZI6)") LINE
                                  write(*,fmt="('               Usage:  &GROUPSUBCTL  t1, t2, ...,meth')")
                                  write(*,fmt="('                       tx, the atom types the temperature control applied to')")
                                  write(*,fmt="(' Process to be stopped')")
                                  stop
                               else
                                  do I=1, N
                                     J = ISTR(STRNUMB(I))
                                     if(J.le. 0 .or. J .gt. mp_mxGROUP) then
                                        write(*,fmt="(' MDPSCU Error: wrong atomic types at line', BZI6)") LINE
                                        write(*,fmt="(' Process to be stopped')")
                                        stop
                                     end if
                                     if(IGS0(J) .gt. 0) then
                                        write(*,fmt="(A, BZI6, A)") 'MDPSCU Error: atomic type ', J, ' already assigned a contral model'
                                        write(*,fmt="(A, BZI6)")    '              check control file at line:', LINE
                                        write(*,fmt="(' Process to be stopped')")
                                        stop
                                     end if
                                     IGS(I)  = J
                                     IGS0(J) = 1
                                  end do
                               end if

                               call Default_TiCtrlParam(TCtrl)
                               !--- to extract the method for local temperature control
                               call Extract_Substr(STR, size(STRNUMB), N, STRNUMB)
                               if(N .le. 0) then !EPC of all atoms is disabled
                                  write(*,fmt="(A, BZI6)") 'MDPSCU Error: the method variable is missed at line', LINE
                                  write(*,fmt="(A)")       '              Usage:  &GROUPSUBCTL  t1, t2, ...,meth'
                                  write(*,fmt="(A)")       '                      tx, the atom types the temperature control applied to'
                                  write(*,fmt="(A)")       '                      meth should be one of "EPC", "ST", "EPC-ST"'
                                  write(*,fmt="(' Process to be stopped')")
                                  stop
                               else
                                   call UpCase(STRNUMB(1))
                                   select case(STRNUMB(1))
                                          case("EPC")
                                              TCtrl%METH = CP_TICTRL_METH_EPC

                                          case("ST")
                                              TCtrl%METH = CP_TICTRL_METH_ST

                                          case("EPC-ST","EPCST","ST-EPC","STEPC")
                                               TCtrl%METH = CP_TICTRL_METH_EPCST
                                   end select
                               end if
                               call Load_TiCtrlParam(hFile, TCtrl, STR, LINE)
                               do I=1, size(IGS)
                                  if(IGS(I) .gt. 0) then
                                     CtrlParam%LT_CTRL(IGS(I)) = TCtrl
                                  end if
                               end do

                  end select
              end do
             !*** end of tmeprature control section
             return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading temperature control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_TemperaturCtl_Local
  !****************************************************************************

  !****************************************************************************
  subroutine Load_TimestepCtl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,          intent(in)::hFile
     type(SimMDCtrl)             ::CtrlParam
     character*(*)               ::STR
     integer                     ::LINE
     !--- local variables
      character*32::STRNUMB(10),KEYWORD
      integer::I, N, TR, TS
     !----

           !*** start the time step control controlling parametgers
           TR = 0
           TS = 0
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
                            CtrlParam%TEMTM = 1.D60  !to time in give an large value
                          else
                            CtrlParam%TEMTM = DRSTR(STRNUMB(2)) !the time in ps for terminating the section
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

                    case ("&OUTPUT_C")
                         !$$--- To get for how many time steps to print out configuration
                          call Extract_Numb(STR,2,N,STRNUMB)
                          if(N .LT. 1) then
                             write(*,fmt="(' MDPSCU Warning: 1 control parameter controling output of configures is expected')")
                             write(*,fmt="('                 check control file at line:', BZI6)") LINE
                             write(*,fmt="(' Usage:  &OUTPUT_C number of timesteps')")
                             call ONWARNING(gm_OnWarning)
                           else
                             CtrlParam%TimestpG = ISTR(STRNUMB(1))
                          end if

                          call Extract_Substr(STR,1, N,STRNUMB)
                          if(N .ge. 1) then
                             call UpCase(STRNUMB(1)) 
                             if(STRNUMB(1) .eq. "M" .or. STRNUMB(1) .eq. "MULT" .or. STRNUMB(1) .eq. "YES") then
                                CtrlParam%MultOutputG = 1     
                             else if(STRNUMB(1) .eq. "S" .or. STRNUMB(1) .eq. "SINGLE" .or. STRNUMB(1) .eq. "NO") then
                                CtrlParam%MultOutputG = 0
                             end if    
                          end if
                          
                   case ("&SAVE")
                         !$$--- To get for how many time steps to save current statu
                         call Extract_Numb(STR,1,n,STRNUMB)
                         if(N .LT. 1) then
                            write(*,fmt="(' MDPSCU Warning: 1 control parameter controling save of current statu is expected')")
                            write(*,fmt="('                 check control file at line:', BZI6)") LINE
                            write(*,fmt="(' Usage:  &SAVE number of timesteps')")
                            call ONWARNING(gm_OnWarning)
                          else
                             TS = ISTR(STRNUMB(1))
                             if(TS .eq. 0) TS = CtrlParam%ITE
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
                            TR = ISTR(STRNUMB(1))
                         end if

              end select
           end do
           !*** end the time step control controlling parametgers
           if(TR .lt. 0) then
              CtrlParam%TimestpR = CtrlParam%TimestpG
           else
              CtrlParam%TimestpR = TR
           end if

           if(TS .lt. 0) then
              CtrlParam%TimestpSave = CtrlParam%ITE
           else
              CtrlParam%TimestpSave = TS
           end if
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
     integer, intent(in):: hFile
     type(SimMDCtrl)    :: CtrlParam
     type(SimMDBox)     :: SimBox
     character*(*)      :: STR
     integer            :: LINE, ISECT
     !--- local variables
      character*32::KEYWORD, SUBKWD(2)
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
                             call Extract_Substr(STR,1,N,SUBKWD)
                             if(N .gt. 0) then
                                call UpCase(SUBKWD(1))
                                select case(SUBKWD(1))
                                       case default
                                            write(*,fmt="(A,BZI6)") 'MDPSCU warning: wrong temperature control type "' &
                                                       //SUBKWD(1)(1:LEN_TRIM(SUBKWD(1)))//'" at line ', LINE
                                            write(*,fmt="(A,BZI6)") '                use "GLOBAL" or "GL" or "G" for global control '
                                            write(*,fmt="(A,BZI6)") '                use "LOCAL" or "LOC" or "GROUP" or "GR" for local control '
                                            write(*,fmt="(A,BZI6)") '                use "DISABLE" or "NONE" or "NULL" to disable temp. control '
                                            call ONWARNING(gm_OnWarning)

                                       case("GLOBAL", "GL", "G")
                                             CtrlParam%TI_CTRL = CP_TICTRL_GLOBAL
                                       case("LOCAL", "LOC", "L", "GROUP", "GR")
                                             CtrlParam%TI_CTRL = CP_TICTRL_LOCAL
                                       case("DISABLE", "NONE", "NULL")
                                             CtrlParam%TI_CTRL = CP_TICTRL_NONE
                                end select
                            else
                               !$$--- add error message
                                CtrlParam%TI_CTRL = CP_TICTRL_GLOBAL
                            end if

                            !*** start the temprature controlling parametgers
                            select case(CtrlParam%TI_CTRL)
                                   case(CP_TICTRL_GLOBAL)
                                        call Load_TemperaturCtl_Global(hFile, CtrlParam, STR, LINE)
                                   case(CP_TICTRL_LOCAL)
                                        call Load_TemperaturCtl_Local(hFile, CtrlParam, STR, LINE)
                            end select

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

                       case("&LBFGSSUBCTL")
                             call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case("&STEEPESTSUBCTL", "&CGSUBCTL")
                             call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&EVENTDETECTSUBCTL")
                             call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                             
                end select
              end do

              !$$*** to check the consistent of the control parameters
              if(CtrlParam%TI_CTRL .eq. CP_TICTRL_GLOBAL) then
                 if(CtrlParam%DAMPTIME1 .gt. 0 .and. CtrlParam%IVTIME .gt. 0) then
                    write(*,fmt="('MDPSCU warning: quick-damping and thermalizing to be performed in the same timesection #', I3)") ISECT
                    call ONWARNING(gm_OnWarning)
                 end if

                 if(CtrlParam%EPCTIME1 .gt. 0 .and. CtrlParam%IVTIME .gt. 0) then
                    write(*,fmt="('MDPSCU warning: E-P coupling and thermalizing to be performed in the same timesection #', I3)") ISECT
                    call ONWARNING(gm_OnWarning)
                 end if

                 if(CtrlParam%DAMPTIME1 .gt. 0 .and. CtrlParam%EPCTIME1 .gt. 0) then
                    write(*,fmt="('MDPSCU warning: E-P coupling and quick-damping to be performed in the same timesection #', I3)") ISECT
                    call ONWARNING(gm_OnWarning)
                 end if
              else if(CtrlParam%TI_CTRL .eq. CP_TICTRL_LOCAL) then
                 if(all(CtrlParam%LT_CTRL%METH .eq. CP_TICTRL_NONE) ) then
                    write(*,*)" MDPSCU Warning: no method set to control temperature of any atomic type"
                    write(*,*)"                 check temperature control for atom groups"
                    call ONWARNING(gm_OnWarning)
                 end if
              end if
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
     integer,               intent(in):: hFile
     type(SimMDCtrl),target           :: CtrlParam
     type(SimMDBox)                   :: SimBox
     integer, optional,     intent(in):: LINE0
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

                       case ("&STSUBCTL", "&STOPSUBCTL")
                             call Load_STCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                            !*** to synchronize the stopping data in the list
                             nextp=>CtrlParam
                             do while(.TRUE.)
                                call GetNext_SimMDCtrl(nextp, 1, next)
                                if(.not. associated(next)) exit
                                next%ST_CTRL    = CtrlParam%ST_CTRL
                                nextp=>next
                             end do

                       case("&LBFGSSUBCTL")
                             call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                             !*** to synchronize the LBFGS data in the list
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

                       case ("&EVENTDETECTSUBCTL")
                             call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)


                end select
              end do
              !--- make stopping power model to be be the same for all section


              !---
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
     integer,         intent(in)::hFile
     type(SimMDCtrl), target    ::CtrlParam
     type(SimMDBox)             ::SimBox

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

                 !*************************
                 call Print_TiCtrl_SimMDCtrl(hFile, sectCtrlParam, SimBox)

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
                    write(hFile,FMT="(' !      Min step number to update activation region.....: ', 3(BNI4,1x))") sectCtrlParam%AR_UPTABMI
                    write(hFile,FMT="(' !      Max step number to update activation region.....: ', 3(BNI4,1x))") sectCtrlParam%AR_UPTABMX
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

                 !--- the input for output control
                 write(hFile,FMT="(' !      Number of timesteps of output of thermal quantities: ', 20(BNI8,1x))") sectCtrlParam%TimestpQ
                 write(hFile,FMT="(' !      Number of timesteps of output of instant configures: ', 20(BNI8,1x))") sectCtrlParam%TimestpG
                 write(hFile,FMT="(' !      Number of timesteps of save systm statu ...........: ', 20(BNI8,1x))") sectCtrlParam%TimestpSave
                 write(hFile,FMT="(' !      Number of timesteps of invoking external recording.: ', 20(BNI8,1x))") sectCtrlParam%TimestpR
                 write(hFile,fmt="(' !')")

                 call GetNext_SimMDCtrl(sectCtrlParam, 1, next)
                 sectCtrlParam=>next
           end do ! end do for sections

         return
  end subroutine Print_Parameter_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Print_TiCtrl_SimMDCtrl(hFile, CtrlParam, SimBox)
  !***  PURPOSE:   to print out control parameters for dnamics applications
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer,               intent(in)::hFile
     type(SimMDCtrl),target           ::CtrlParam
     type(SimMDBox)                   ::SimBox

     !--- local variables
     integer::I, J
     !----
              if(CtrlParam%TI_CTRL .eq. CP_TICTRL_GLOBAL) then
                 if(CtrlParam%DAMPTIME1  .LE. 0  .and. &
                    CtrlParam%IVTIME     .LE. 0  .and. &
                    CtrlParam%EPCTIME1   .LE. 0  .and. &
                    CtrlParam%IFNOSE     .LE. 0         ) then
                    write(hFile,FMT="(' !      Temperature control.........................:  NONE ', 3(BZF8.2,1x))")
                 else
                    write(hFile,FMT="(' !      Temperature control.........................:  GLOBAL ', 3(BZF8.2,1x))")
                    write(hFile,FMT="(' !      Temperature(K)..............................: ', 3(BZF8.2,1x))")   CtrlParam%TI
                    if(CtrlParam%DAMPTIME1 .GT. 0) then
                       write(hFile,FMT="(' !       QuickDamping to be performed in section....: ', BZI5)")
                       write(hFile,FMT="(' !         start from the timestep .................: ', BZI8)") CtrlParam%DAMPTIME0
                       write(hFile,FMT="(' !         with the number of time steps............: ', BZI8)") CtrlParam%DAMPTIME1
                       select case( iand(CtrlParam%DAMPSCHEME, CP_LWORD) )
                              case (CP_DAMPSCHEME_LBFGS)
                                   write(hFile,FMT="(' !         damping scheme...........................: LBFGS method ', BZI8)")
                                   write(hFile,FMT="(' !                 the convergence criteria.........: PGTOL ', 3(1PE11.4,1x))") CtrlParam%LBFGS_PGtol
                                   write(hFile,FMT="(' !                                                    FACTR ', 3(1PE11.4,1x))") CtrlParam%LBFGS_Factr
                                   write(hFile,FMT="(' !                                                    MSAVE ', BZI8)")          CtrlParam%LBFGS_MSave

                              case (CP_DAMPSCHEME_ST)
                                    if(iand(CtrlParam%DAMPSCHEME, CP_DAMPSCHEME_LSEARCH) .gt. 0) then
                                       write(hFile,FMT="(' !         damping scheme...........................: STEEPEST with linesearch ', BZI8)")
                                    else   
                                       write(hFile,FMT="(' !         damping scheme...........................: STEEPEST method ', BZI8)")
                                    end if
                                    write(hFile,FMT="(' !                 inital alpha ....................: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_Alpha
                                    write(hFile,FMT="(' !                 max move in a step MXSTEP (LU)...: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_MxStep
                                    if(CtrlParam%STEEPEST_MiStep .gt. 0) &
                                       write(hFile,FMT="(' !                 convergence criteria MISTEP (LU).: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_MiStep
                                    if(CtrlParam%STEEPEST_MiDelE .gt. 0) &
                                       write(hFile,FMT="(' !                 convergence criteria DELTAE (eV).: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_MiDelE
                              case (CP_DAMPSCHEME_CG)
                                    if(iand(CtrlParam%DAMPSCHEME, CP_DAMPSCHEME_LSEARCH) .gt. 0) then
                                       write(hFile,FMT="(' !         damping scheme...........................: CG with linesearch ', BZI8)")
                                    else   
                                       write(hFile,FMT="(' !         damping scheme...........................: CG method ', BZI8)")
                                    end if
                                    write(hFile,FMT="(' !                 max move in a step MXSTEP (LU)...: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_MxStep
                                    if(CtrlParam%STEEPEST_MiStep .gt. 0) &
                                       write(hFile,FMT="(' !                 convergence criteria MISTEP (LU).: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_MiStep
                                    if(CtrlParam%STEEPEST_MiDelE .gt. 0) &
                                       write(hFile,FMT="(' !                 convergence criteria DELTAE (eV).: ', 3(1PE11.4,1x))") CtrlParam%STEEPEST_MiDelE

                              case (CP_DAMPSCHEME_DYN)
                                    write(hFile,FMT="(' !         damping scheme...........................: DYNAMICS method ', BZI8)")
                         end select
                    end if

                    if(CtrlParam%IVTIME .GT. 0) then
                       write(hFile,FMT="(' !       Thermalization to be performed in section..: ', BZI5)")
                       write(hFile,FMT="(' !          start from the time step ...............: ', BZI8)") CtrlParam%IVTIME0
                       write(hFile,FMT="(' !          num of thermalizing in the section......: ', BZI8)") CtrlParam%IVTIME
                       write(hFile,FMT="(' !          num of time steps between thermalizing..: ', BZI8)") CtrlParam%IVPAS
                       select case(CtrlParam%IVSCHEME)
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

                    if(CtrlParam%EPCTIME1 .GT. 0) then
                       write(hFile,FMT="(' !       E-P couipling to be performed in section..: ', BZI5)")
                       write(hFile,FMT="(' !         start from the timestep .................: ', BZI8)") CtrlParam%EPCTIME0
                       write(hFile,FMT="(' !         with the number of time steps............: ', BZI8)") CtrlParam%EPCTIME1
                       call Print_EPC_TiCtrlParam(hFile, CtrlParam%LT_CTRL(1))
                       !write(hFile,FMT="(' !         with physical model defined in ..........: ', A)")    CtrlParam%EPCPHYS(1:len_trim(CtrlParam%EPCPHYS))
                    end if
                 end if
                 write(hFile,fmt="(' !')")
             else if(CtrlParam%TI_CTRL .eq. CP_TICTRL_LOCAL) then
                 write(hFile,FMT="(' !      Temperature control.........................: LOCAL', 3(BZF8.2,1x))")
                 do I=1, size(CtrlParam%LT_CTRL)
                    if(iand(CtrlParam%LT_CTRL(I)%METH, CP_LWORD) .eq. CP_TICTRL_NONE) then
                       !--- do nothing
                    else
                       write(hFile,FMT="(' !        For atomic type .........................#: ', BZI3)"), I
                       call Print_TiCtrlParam(hFile, CtrlParam%LT_CTRL(I))
                       if(iand(CtrlParam%LT_CTRL(I)%METH,CP_TICTRL_METH_ST) .eq. CP_TICTRL_METH_ST) then
                          call Print_STCtrlParam(hFile, CtrlParam%ST_CTRL)
                       end if
                    end if
                 end do
             else
                 write(hFile,FMT="(' !      Temperature control.........................:  NONE ', 3(BZF8.2,1x))")
             end if

         return
  end subroutine Print_TiCtrl_SimMDCtrl
  !****************************************************************************

  end module MD_SimMDCtrl_GMD
