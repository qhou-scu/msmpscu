  module MD_TYPEDEF_EPCCtrl
  !***  DESCRIPTION: this module is to define the data type for temperature control
  !                  ______________________________________________________
  !
  !     HISORY:      2018-03-07, created by HOU Qing
  use MD_CONSTANTS
  use MiniUtilities
  use MD_TYPEDEF_InputPaser

  implicit none
   !--- The data type for controlling local temperature
   !--- The data type for controlling local temperature
      type::STCtrlParam
           integer      ::MOD    = CP_STCTRL_MOD_LS_B     ! the model for ST
           real(KINDDF) ::LE     = 10.D0*CP_EVERG         ! the energy threshold higher than which stopping model to be is applied
           real(KINDDF) ::HE     = 10000.D0*CP_EVERG      ! the energy upper limit for the stopping tab;e
           integer      ::NTAB   = 10000                  ! the stoping table size
           real(KINDDF) ::MDEN   = 0.D0                   ! the density of media in cm^(-3).
                                                          !     <0 , the density to be calcualted locally
                                                          !     =0,  the density calculated by the number of atom divided by box size
           character*256::ExtFile                         ! the filename for load in stopping power when external stopping used

           integer      ::PrintTab  = 1                   ! the flag indicate if output stopping-energy table
           integer      ::SaveELoss = 1                   ! the flag indicate if output inelastic energy loss of atoms
      end type STCtrlParam

      type::TiCtrlParam
           !-- the parameters for local temperature control
           real(KINDDF)::TI        = 0                      ! the local temperature for atomic types
           integer     ::METH      = CP_TICTRL_NONE         ! the method of control local temperature

           integer     ::EPC_MOD   = CP_EPCCTRL_MOD_MAN     ! the model for EPC
           real(KINDDF)::EPC_HE    = 100.D0*CP_EVERG        ! the energy threshold lower than which EPC model to be applied
           real(KINDDF)::EPC_ALPHA = 1.D0*CP_PS2S           ! the epc coupling time in s
           real(KINDDF)::EPC_CUT   = 0.1D0                  ! a parameter determine cutoff temperture when the instance temperature of a atom too low
           integer     ::EPC_SAVEWK= 0                      ! the flag indicating if energy gain-loss of atoms to be output
      end type TiCtrlParam

      interface assignment (=)
          module procedure CopyFrom_TiCtrlParam
      end interface

      interface assignment (=)
          module procedure CopyFrom_STCtrlParam
      end interface

  contains
    !****************************************************************************
     elemental subroutine CopyFrom_TiCtrlParam(To, From)
      !***  PURPOSE:   to copy a control parameter to another
      !
      !     INPUT:     From, the source simulation controal
      !     OUTPUT     To, the copy of FROM
      implicit none

      !--- DUMMY variables
      type(TiCtrlParam)::To, From
      !--- Local variables
           To%TI        = From%TI
           To%METH      = From%METH
           To%EPC_MOD   = From%EPC_MOD
           To%EPC_HE    = From%EPC_HE
           To%EPC_ALPHA = From%EPC_ALPHA
           To%EPC_CUT   = From%EPC_CUT
           To%EPC_SAVEWK= From%EPC_SAVEWK
           return
      end subroutine CopyFrom_TiCtrlParam
    !****************************************************************************

    !****************************************************************************
     elemental subroutine CopyFrom_STCtrlParam(To, From)
      !***  PURPOSE:   to copy a control parameter to another
      !
      !     INPUT:     From, the source simulation controal
      !     OUTPUT     To, the copy of FROM
      implicit none

      !--- DUMMY variables
      type(STCtrlParam)::To, From
      !--- Local variables
           To%MOD     = From%MOD
           To%LE      = From%LE
           To%HE      = From%HE
           To%NTAB    = From%NTAB
           To%MDEN    = From%MDEN
           To%ExtFile = From%ExtFile

           To%PrintTab = To%PrintTab
           To%SaveELoss= To%SaveELoss

           return
      end subroutine CopyFrom_STCtrlParam
    !****************************************************************************

    !****************************************************************************
     elemental subroutine Default_TiCtrlParam(TCtrl)
      !***  PURPOSE:   to copy a control parameter to another
      !
      !     INPUT:     From, the source simulation controal
      !     OUTPUT     To, the copy of FROM
      implicit none

      !--- DUMMY variables
      type(TiCtrlParam)::TCtrl
      !--- Local variables
           TCtrl%TI        = 0.d0
           TCtrl%METH      = CP_TICTRL_NONE
           TCtrl%EPC_MOD   = CP_EPCCTRL_MOD_MAN
           TCtrl%EPC_HE    = 100.D0*CP_EVERG
           TCtrl%EPC_ALPHA = 1.D0*CP_PS2S
           TCtrl%EPC_CUT   = 0.1D0
           TCtrl%EPC_SAVEWK= 0
           return
      end subroutine Default_TiCtrlParam
  !****************************************************************************

      !****************************************************************************
     elemental subroutine Default_STCtrlParam(STCtrl)
      !***  PURPOSE:   to copy a control parameter to another
      !
      !     INPUT:     From, the source simulation controal
      !     OUTPUT     To, the copy of FROM
      implicit none

      !--- DUMMY variables
      type(STCtrlParam)::STCtrl
      !--- Local variables

           STCtrl%MOD    = CP_STCTRL_MOD_LS_B
           STCtrl%LE     = 10.D0*CP_EVERG
           STCtrl%HE     = 10000.D0*CP_EVERG
           STCtrl%MDEN   = 0.D0
           STCtrl%ExtFile = ""
           STCtrl%PrintTab = 1
           STCtrl%SaveELoss= 1
           return
      end subroutine Default_STCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Load_TiCtrlParam(hFile, TCtrl, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !                LINE,   current line  number
  !     OUTPUT     TCtrl
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(TiCtrlParam)::TCtrl
     character*(*)::STR
     integer::LINE
     !--- local variables
     character*32::STRNUMB(1),KEYWORD
           !----

               select case (TCtrl%METH)
                      case(CP_TICTRL_METH_EPC)
                           call Load_EPC_TiCtrlParam(hFile, TCtrl, STR, LINE)

                      case(CP_TICTRL_METH_ST)
                          do while(.TRUE.)
                             call GetInputStrLine(hFile,STR, LINE, "!", *100)
                             STR = adjustl(STR)
                             call GetKeyWord("&", STR, KEYWORD)
                             call UpCase(KEYWORD)
                             select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                                   case( "&ENDSUBCTL")
                                        exit
                                   case default
                                         write(*,fmt="(A)")' MDPSCU Warning: unknown keyword in the SUBCTL section '//KEYWORD(1:LEN_TRIM(KEYWORD))
                                         write(*,fmt="('               check control file at line:', BZI6)") LINE
                                         call ONWARNING(gm_OnWarning)
                              end select
                         end do
                      case(CP_TICTRL_METH_EPCST)
                           call Load_EPC_TiCtrlParam(hFile, TCtrl, STR, LINE)
               end select
             return
  !---------------------------------------------------------------
  100    write(*,fmt="(A, I5)")"MDPSCU Error in reading control parameters at line", LINE
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDCTLF exist to end the data section"
         stop

  end subroutine Load_TiCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Print_TiCtrlParam(hFile, TCtrl)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                TCtrl
  implicit none
     !--- dummy varioables
     integer,           intent(in)::hFile
     type(TiCtrlParam), intent(in)::TCtrl
     !----

              select case(TCtrl%METH)
                     case(CP_TICTRL_METH_EPC)
                          write(hFile,FMT="(' !         Control method...........................: EPC', BZI5)")
                          call Print_EPC_TiCtrlParam(hFile, TCtrl)

                     case(CP_TICTRL_METH_ST)
                          write(hFile,FMT="(' !         Control method...........................: ST', BZI5)")
                          !call Print_ST_TiCtrlParam(hFile, TCtrl)

                     case(CP_TICTRL_METH_EPCST)
                          write(hFile,FMT="(' !         Control method...........................: EPC+ST', BZI5)")
                          call Print_EPC_TiCtrlParam(hFile, TCtrl)
              end select
             return

   end subroutine Print_TiCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Load_EPC_TiCtrlParam(hFile, TCtrl, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !                LINE,   current line  number
  !     OUTPUT     TCtrl
  implicit none
     !--- dummy varioables
     integer,intent(in)::hFile
     type(TiCtrlParam) ::TCtrl
     character*(*)::STR
     integer::LINE
     !--- local variables
      character*32::STRNUMB(1),KEYWORD
      integer::I, J, N
     !---

               do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case( "&ENDSUBCTL")
                                exit
                         case default
                              write(*,fmt="(A)")' MDPSCU Warning: unknown keyword in &EPCSUBCTL section '//KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)

                         case("&TEMPERATURE", "&TEMP", "&TE")
                              !*** the temperature
                               call Extract_Numb(STR,1,n,STRNUMB)
                               TCtrl%TI = DRSTR(STRNUMB(1))

                         case("&SAVEWORK")
                              !*** the temperature
                               call Extract_Substr(STR,1,N,STRNUMB)
                               if(N.gt.0) then
                                 call UpCase(STRNUMB(1))
                                 select case(STRNUMB(1))
                                        case("YES")
                                            TCtrl%EPC_SAVEWK = 1

                                        case("NO")
                                            TCtrl%EPC_SAVEWK = 0

                                        case default
                                            TCtrl%EPC_SAVEWK = 0
                                 end select
                                 TCtrl%TI = DRSTR(STRNUMB(1))
                               end if

                         case("&MODSUBCTL")
                              call Extract_Substr(STR, size(STRNUMB), N, STRNUMB)
                              if(N.gt.0) then
                                 call UpCase(STRNUMB(1))
                                 select case(STRNUMB(1))
                                        !--- the model by Sommerfiled theory
                                        case("SMFT")
                                             TCtrl%EPC_MOD = CP_EPCCTRL_MOD_SMFT
                                             write(*,fmt="(A, BZI6)") 'MDPSCU Error: EPC by SMFT is not supported yet'
                                             write(*,fmt="(' Process to be stopped')")
                                             stop

                                         case("MAN")
                                               TCtrl%EPC_MOD = CP_EPCCTRL_MOD_MAN
                                               call Load_EPC_MAN_TiCtrlParam(hFile, TCtrl, STR, LINE)

                                         case default
                                               write(*,fmt="(A, BZI6)") 'MDPSCU Error: unknown EPC model at', LINE
                                               write(*,fmt="(A)")       '              Usage: &EPCSUBCTL model'
                                               write(*,fmt="(A)")       '                     model should be one of "SMFT", "MAN"'
                                               write(*,fmt="(' Process to be stopped')")
                                               stop
                                 end select
                              else
                                 write(*,fmt="(A, BZI6)") 'MDPSCU Error: unknown EPC model at', LINE
                                 write(*,fmt="(A)")       '              Usage: &MODSUBCTL model'
                                 write(*,fmt="(A)")       '                     model should be one of "SMFT", "MAN"'
                                 write(*,fmt="(' Process to be stopped')")
                                 stop
                              end if
                  end select
               end do
             return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading temperature control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop

   end subroutine Load_EPC_TiCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Print_EPC_TiCtrlParam(hFile, TCtrl)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                TCtrl
  implicit none
     !--- dummy varioables
     integer,           intent(in)::hFile
     type(TiCtrlParam), intent(in)::TCtrl
     !----

              select case(TCtrl%EPC_MOD)
                     case(CP_EPCCTRL_MOD_MAN)
                          call Print_EPC_MAN_TiCtrlParam(hFile, TCtrl)
              end select
             return
   end subroutine Print_EPC_TiCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Load_EPC_MAN_TiCtrlParam(hFile, TCtrl, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !                LINE,   current line  number
  !     OUTPUT     TCtrl
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(TiCtrlParam)::TCtrl
     character*(*)::STR
     integer::LINE
     !--- local variables
      character*32::STRNUMB(2),KEYWORD
      integer::N
     !----

             !*** start the temprature controlling parametgers
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case("&ALPHA")
                              !*** the alpha
                               call Extract_Numb(STR,1,n,STRNUMB)
                               if(N .le. 0) then
                               else
                                  TCtrl%EPC_ALPHA = DRSTR(STRNUMB(1))*CP_PS2S
                               end if

                          case("&EMAX")
                                 call Extract_Numb(STR, 1, N, STRNUMB)
                                 if(N .le. 0) then !EPC of all atoms is disabled
                                 else
                                   TCtrl%EPC_HE = DRSTR(STRNUMB(1))*CP_EVERG
                                end if

                          case("&LTCUT")
                                 call Extract_Numb(STR, 1, N, STRNUMB)
                                 if(N .le. 0) then
                                 else
                                   TCtrl%EPC_CUT = DRSTR(STRNUMB(1))
                                end if

                         case( "&ENDSUBCTL")
                                exit
                         case default
                              write(*,*)" MDPSCU Warning: unknown keyword in EPC control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)

                  end select
              end do
             !*** end of EPC control section
             return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading temperature control parameters."
         write(*,*)"The process to be stopped."
         stop
  end subroutine Load_EPC_MAN_TiCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Print_EPC_MAN_TiCtrlParam(hFile, TCtrl)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                TCtrl
  implicit none
     !--- dummy varioables
     integer,           intent(in)::hFile
     type(TiCtrlParam), intent(in)::TCtrl

     !----
             write(hFile,FMT="(' !         EPC-model................................: manully set coupling time')")
             write(hFile,FMT="(' !             Electron Temperature(K)..............: ', 1PE10.3)") TCtrl%Ti
             write(hFile,FMT="(' !             with alpha set (ps)..................: ', 1PE10.3)") TCtrl%EPC_ALPHA*CP_S2PS
             write(hFile,FMT="(' !             below energy (eV)....................: ', 1PE10.3)") TCtrl%EPC_HE*CP_ERGEV
             write(hFile,FMT="(' !             low temperature cut ratio............: ', 1PE10.3)") TCtrl%EPC_CUT
             if(TCtrl%EPC_SAVEWK .gt. 0) then
               write(hFile,FMT="(' !             work done by EPC to be output........: ', 1PE10.3)")
             else
               write(hFile,FMT="(' !             work done by EPC not to be output....: ', 1PE10.3)")
             end if
             return
   end subroutine Print_EPC_MAN_TiCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Load_STCtrlParam(hFile, STCtrl, STR, LINE)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !                LINE,   current line  number
  !     OUTPUT     TCtrl
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(STCtrlParam)::STCtrl
     character*(*)::STR
     integer::LINE
     !--- local variables
      character*32:: STRNUMB(2),KEYWORD
      character*256::FNAME(2)
      equivalence (FNAME(1),STRNUMB(1))
      integer::N
     !----
             !*** start the temprature controlling parametgers
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                          case("&EMIN")
                                 call Extract_Numb(STR, 1, N, STRNUMB)
                                 if(N .le. 0) then
                                 else
                                   STCtrl%LE = DRSTR(STRNUMB(1))*CP_EVERG
                                end if

                          case("&EMAX")
                                 call Extract_Numb(STR, 1, N, STRNUMB)
                                 if(N .le. 0) then
                                 else
                                   STCtrl%HE = DRSTR(STRNUMB(1))*CP_EVERG
                                end if

                          case("&MDEN")
                                 call Extract_Numb(STR, 1, N, STRNUMB)
                                 if(N .le. 0) then
                                   STCtrl%MDEN = 0.D0
                                 else
                                   STCtrl%MDEN = DRSTR(STRNUMB(1))
                                end if

                          case("&MODEL", "&MOD")
                                 call Extract_Substr(STR, 2, N, STRNUMB)
                                 if(N .gt. 0) then
                                    call UpCase(STRNUMB(1))
                                    select case(STRNUMB(1))
                                           case("LS_Z85")
                                                STCtrl%MOD = CP_STCTRL_MOD_LS_Z85

                                           case("LS_Z95")
                                                STCtrl%MOD = CP_STCTRL_MOD_LS_Z95

                                           case("LS_B")
                                                STCtrl%MOD = CP_STCTRL_MOD_LS_B


                                           case("OR")
                                                STCtrl%MOD = CP_STCTRL_MOD_OR

                                          case ("SRIM","TRIM")
                                                if(N .gt.1)then
                                                    STCtrl%MOD = CP_STCTRL_MOD_SRIM
                                                    call Extract_Substr(STR, 2, N, FNAME)
                                                    STCtrl%ExtFile = FNAME(2)
                                                else
                                                    write(*,fmt="(A, BZI6)") 'MDPSCU Error: filename for input SRIM stopping is missed at', LINE
                                                    write(*,fmt="(A)")       '              Usage: &MODSUBCTL "SRIM" filename'
                                                    write(*,fmt="(' Process to be stopped')")
                                                    stop
                                                end if
                                           case default
                                                STCtrl%MOD = CP_STCTRL_MOD_USER
                                                call Extract_Substr(STR, 1, N, FNAME)
                                                STCtrl%ExtFile = FNAME(1)

                                     end select
                                 else
                                   write(*,fmt="(A, BZI6)") 'MDPSCU Error: unknown ST model at', LINE
                                   write(*,fmt="(A)")       '              Usage: &MODSUBCTL model'
                                   write(*,fmt="(A)")       '                     model should be one of "LS", "OR"'
                                   write(*,fmt="(' Process to be stopped')")
                                   stop
                                end if

                          case("&PRINTTAB")
                                 call Extract_Substr(STR, 1, N, STRNUMB)
                                 if(N .gt. 0) then
                                    call UpCase(STRNUMB(1))
                                    select case(STRNUMB(1))
                                           case("YES")
                                                STCtrl%PrintTab = 1

                                           case("NO")
                                                STCtrl%PrintTab = 0

                                           case default
                                                STCtrl%PrintTab = 0
                                     end select
                                 else
                                   STCtrl%PrintTab = 0
                                 end if

                          case("&SAVEELOSS")
                                 call Extract_Substr(STR, 1, N, STRNUMB)
                                 if(N .gt. 0) then
                                    call UpCase(STRNUMB(1))
                                    select case(STRNUMB(1))
                                           case("YES")
                                                STCtrl%SaveELoss = 1

                                           case("NO")
                                                STCtrl%SaveELoss = 0

                                           case default
                                                STCtrl%SaveELoss = 0
                                     end select
                                 else
                                   STCtrl%SaveELoss = 0
                                 end if

                         case( "&ENDSUBCTL")
                                exit
                         case  default
                               write(*,*)" MDPSCU Warning: unknown keyword in STOP control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                               write(*,fmt="('               check control file at line:', BZI6)") LINE
                               call ONWARNING(gm_OnWarning)

                  end select
              end do
             !*** end of tmeprature control section
             return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading temperature control parameters."
         write(*,*)"The process to be stopped."
         stop
  end subroutine Load_STCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Print_STCtrlParam(hFile, STCtrl)
  !***  PURPOSE:   to load the control parameter for temperature control of groups
  !
  !     INPUT:     hFile,  I/O unit number
  !                TCtrl
  implicit none
     !--- dummy varioables
     integer,           intent(in)::hFile
     type(STCtrlParam), intent(in)::STCtrl

         select case(STCtrl%MOD)
                case(CP_STCTRL_MOD_LS_Z95)
                     write(hFile,FMT="(' !         ST-model.................................: LS_Z95')")
                case(CP_STCTRL_MOD_LS_Z85)
                     write(hFile,FMT="(' !         ST-model.................................: LS_Z85')")
                case(CP_STCTRL_MOD_LS_B)
                     write(hFile,FMT="(' !         ST-model.................................: LS_B')")
                case(CP_STCTRL_MOD_OR)
                     write(hFile,FMT="(' !         ST-model.................................: OR')")

                case(CP_STCTRL_MOD_USER)
                     write(hFile,FMT="(' !         ST-model.................................: Load from ', A)") &
                                         '"'//STCtrl%ExtFile(1:len_trim(STCtrl%ExtFile))//'"'

                case(CP_STCTRL_MOD_SRIM)
                     write(hFile,FMT="(' !         ST-model.................................: SRIM stopping from ', A)") &
                                         '"'//STCtrl%ExtFile(1:len_trim(STCtrl%ExtFile))//'"'

         end select
         write(hFile,FMT="(' !             above energy (eV)....................: ', 1PE10.3)") STCtrl%LE*CP_ERGEV
         write(hFile,FMT="(' !             below energy (eV)....................: ', 1PE10.3)") STCtrl%HE*CP_ERGEV

         return
  end subroutine Print_STCtrlParam
  !****************************************************************************

  end module MD_TYPEDEF_EPCCtrl
