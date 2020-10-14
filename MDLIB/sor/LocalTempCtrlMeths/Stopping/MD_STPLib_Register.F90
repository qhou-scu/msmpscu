  module MD_STPLib_Register
  !***  DESCRIPTION:
  !     This module provides a interface to register to diffrent stopping power
  !     calculation that to be usd in MD simulation.
  !
  !     DEPENDENCE: MD_TYPEDEF_STOPRANGE_TABLE
  !
  !     SEE ALSO:
  !
  !****
  !     HISTORY:     2012-03(HOU Qing): created:
  !                          The codes of calculating stopping power are adapted from the
  !                          codes for ion-tranpsport calculation in the folder:
  !                          Transport\Panda\PANDAPHY.  The Zigler-85 (-Z85) and Biersack
  !                          (-B) methods are implemented.
  !
  !                          The implementation has been validated by comparing with our
  !                          old stopping calulation codes that can be found in the folder
  !                          Transport/Demo/.
  !
  !
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_TYPEDEF_StopRange_Table
   implicit none

       type(MDSTPTable)::m_STPTable

      !---------------------------------------------------------
       private:: Register_Internal_STPLib
       private:: Register_External_STPLib
      !---------------------------------------------------------

  contains

  !****************************************************************************************
  subroutine Register_Internal_STPLib(SimBox, CtrlParam, STTable)
  !***  DESCRIPTION: to load a the ST lib according to control parameters
  !     INPUT:       SimBox:     the simulation box
  !                  CtrlParam:  the control parameters
  !
  !     OUPUT:       STTable:    the ST table for atoms
  !

  use STOP_B_MODULE,   only:ESTOP_B=>ESTOP
  use STOP_Z85_MODULE, only:ESTOP_Z85=>ESTOP
  use STOP_Z95_MODULE, only:ESTOP_Z95=>ESTOP
  !--- NOTE:  the units used for ESTOP: keV for Enengy, keV*cm^2 for
  !           stopping, atomic mass for mass of atoms.
  !           the unit in MDSTPTableis CGS: erg for energy,  and
  !           erg*cm^2 for stopping
  !
  !           To convert the stopping to force, the number of density
  !           should be multiplied to the stopping.
  !
  implicit none
    type(SimMDBox)   ::SimBox
    type(SimMDCtrl)  ::CtrlParam
    type(MDSTPTable) ::STTable

  !$$--- local varibales
   integer::I, J
   real(KINDDF)::M1, M2, Z1, Z2

          call Release_STPTable(STTable)
          do I=1, SimBox%NGROUP
             M1 = SimBox%CM(I)/CP_AU2G
             Z1 = SimBox%CZ(I)

             do J=1, SimBox%NGROUP
                M2 = SimBox%CM(J)/CP_AU2G
                Z2 = SimBox%CZ(J)

                select case(CtrlParam%ST_CTRL%MOD)
                        case(CP_STCTRL_MOD_LS_B)
                             call Register_STPTableProc(M1, M2, Z1, Z2, STTable, ESTOP=ESTOP_B, NOTE="LS_B stopping:")

                        case(CP_STCTRL_MOD_LS_Z85)
                             call Register_STPTableProc(M1, M2, Z1, Z2, STTable, ESTOP=ESTOP_Z85, NOTE="LS_Z85 stopping:")

                        case(CP_STCTRL_MOD_LS_Z95)
                             call Register_STPTableProc(M1, M2, Z1, Z2, STTable, ESTOP=ESTOP_Z95, NOTE="LS_Z95 stopping:")

                        case(CP_STCTRL_MOD_OR)
                             write(*,*) "Support for OR stopping is still in procceding"
                             stop
                 end select
             end do
          end do

          if(CtrlParam%ST_CTRL%PrintTab .gt. 0) then
             call Create_STPTable(SimBox, CtrlParam, STTable, "Stopping_table")
          else
             call Create_STPTable(SimBox, CtrlParam, STTable)
          end if
       return
   end subroutine Register_Internal_STPLib
  !****************************************************************************************

  !****************************************************************************************
  subroutine Register_External_STPLib(SimBox, CtrlParam, STTable)
  !***  DESCRIPTION: to load a the ST lib according to control parameters
  !     INPUT:       SimBox:     the simulation box
  !                  CtrlParam:  the control parameters
  !
  !     OUPUT:       STTable:    the ST table for atoms
  !

  !--- NOTE:  the units used for ESTOP: keV for Enengy, keV*cm^2 for
  !           stopping, atomic mass for mass of atoms.
  !           the unit in MDSTPTableis CGS: erg for energy,  and
  !           erg*cm^2 for stopping
  !
  !           To convert the stopping to force, the number of density
  !           should be multiplied to the stopping.
  !
  use STOP_SRIM_MODULE, only:Creat_SRIM_STPTable

  implicit none
    type(SimMDBox)   ::SimBox
    type(SimMDCtrl)  ::CtrlParam
    type(MDSTPTable) ::STTable

  !$$--- local varibales
    logical:: EX
    integer:: hFile
    character*256::fname,fname1


    select case(CtrlParam%ST_CTRL%MOD)
        case(CP_STCTRL_MOD_USER)
             fname = CtrlParam%ST_CTRL%ExtFile(1:len_trim(CtrlParam%ST_CTRL%ExtFile))//".stp"
             inquire(FILE=fname, EXIST=EX)
             if(.not. EX) then
                write(*,fmt="(A)") ' MDPSCU Error: cannot find required external stop-table file: '
                write(*,fmt="(A)") '               '//fname(1:len_trim(fname))
                write(*,fmt="(A)") '               process to be stop'
                stop
             end if
             call Release_STPTable(STTable)
             call Import_STPTable(SimBox, CtrlParam, STTable, fname)
             if(CtrlParam%ST_CTRL%PrintTab .gt. 0) then
                call Export_STPTable(SimBox, CtrlParam, STTable, "Stopping_table_" &
                   //CtrlParam%ST_CTRL%ExtFile(1:len_trim(CtrlParam%ST_CTRL%ExtFile)))
             endif

         case(CP_STCTRL_MOD_SRIM)
             fname = CtrlParam%ST_CTRL%ExtFile(1:len_trim(CtrlParam%ST_CTRL%ExtFile))//".srim"
             inquire(FILE=fname, EXIST=EX)
             if(.not. EX) then
                write(*,fmt="(A)") ' MDPSCU Error: cannot find required SRIM input file: '
                write(*,fmt="(A)") '               '//fname(1:len_trim(fname))
                write(*,fmt="(A)") '               process to be stop'
                stop
             end if
             call Release_STPTable(STTable)
             call Creat_SRIM_STPTable(fname,fname1)
             fname1 = fname1(1:len_trim(fname1))//".stp"
             call Import_STPTable(SimBox, CtrlParam, STTable, fname1)

             call AvailableIOUnit(hFile)
             open (UNIT = hFile , FILE = fname1)
             close(hFile, STATUS= 'delete')

             if(CtrlParam%ST_CTRL%PrintTab .gt. 0) then
                call Export_STPTable(SimBox, CtrlParam, STTable, "Stopping_table_" &
                   //CtrlParam%ST_CTRL%ExtFile(1:len_trim(CtrlParam%ST_CTRL%ExtFile)))
             endif
         end select
      return
   end subroutine Register_External_STPLib
  !****************************************************************************************

  !****************************************************************************************
  subroutine Register_STPLib(SimBox, CtrlParam, STTable)
  !***  DESCRIPTION: to load a the ST lib according to control parameters
  !     INPUT:       SimBox:     the simulation box
  !                  CtrlParam:  the control parameters
  !
  !     OUPUT:       STTable:    the ST table for atoms
  !

  implicit none
    type(SimMDBox)   ::SimBox
    type(SimMDCtrl)  ::CtrlParam
    type(MDSTPTable) ::STTable

  !$$--- local varibales

          if(CtrlParam%ST_CTRL%MOD .eq. CP_STCTRL_MOD_USER .or. &
             CtrlParam%ST_CTRL%MOD.eq.CP_STCTRL_MOD_SRIM)       then
             call Register_External_STPLib(SimBox, CtrlParam, STTable)
          else
             call Register_Internal_STPLib(SimBox, CtrlParam, STTable)
          end if
       return
   end subroutine Register_STPLib
  !****************************************************************************************


  end module MD_STPLib_Register
  !****************************************************************************************
