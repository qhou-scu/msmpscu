 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to generate and save the potential tables. The potential could
 !                  be the EAM or FS potentials internally implemented in MDPSCU. Also, the potential could
 !                  be user-defined EAM or FS potential. If the table of user-define potential will be
 !                  generated, some minor-modification in the program FORCETABLE_GENERATOR_Main
 !                  should be made (see the exmaples inlcude the program).
 !
 !                  The generated tables of pairewise functions, including pairwise potentials and the
 !                  electron densities of atoms, will be save in a file with extension ".pair".
 !                  The embedment function F(RHO) will save in a file with ".embd".
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  Commandline command:  ForceTableGenerator.exe  -LIB pname [OPTIONS]
 !                  -LIB,        indicates the name (pname) of potential. The tables will be
 !                               saved in pname.pair and pname.embd
 !                  [OPTIONE]
 !
 !                   -GR psname, the subname of potential
 !                   -RC rcut,   the cutoff range of the potential (in A)
 !                   -SZ nts,    the number of points in the potential table

 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  by HOU Qing, Dec., 2014
 !
  module MD_FORCETABLE_GENERATOR
  !***  DESCRIPTION:
  !     This program is to generate and save the potential tables implemented in
  !     MDPSCU
  !
  !     DEPENDENCE: MD_ForceLib_Factory_CPU
  !
  !     by HOU Qing, Dec., 2014

  !***  The modules included ******************************************
  use MD_CONSTANTS
  use MiniUtilities
  implicit none

  contains
  !****************************************************************************************
  subroutine Generator_Main(FORCETABLE, POTTYPE)
  !***  DESCRIPTION:
  !
  !--- INPUT:
  !           FORCETABLE,   optional, the subroutine provided by user for force register
  !           POTTYPE,      optional, the identifier of the type of potential
  !
  use MiniUtilities
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_ForceTable
  use MD_ForceLib_Register, only:Export_Internal_ForceLib, Export_External_ForceLib
  use MD_Globle_Variables, only:ExtractExeName_Globle_Variables

  implicit none

  !--- dummy variables and subroutines
  character*(*), optional::POTTYPE
  !------------------------------------------
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
      integer,optional::printout
     end subroutine FORCETABLE
   end interface
  !--- END INTERFACE --------------------------------


  !--- Local variables
  integer::I, J,K, NP=1, IERR, ARGN, count

  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2
  type(SimMDBox) ::SimBox0
  type(SimMDCtrl)::CtrlParam

  character*256::ARGV, tstr(16)=""
  character*512::cmdline


           call ExtractExeName_Globle_Variables()
           call get_command(cmdline)
           I = index(cmdline, gm_ExeName(1:len_trim(gm_ExeName)))
           I = I+len_trim(gm_ExeName)
           cmdline = cmdline(I:len_trim(cmdline))

           call extract_optstr1(cmdline, "-", "LIB", tstr(1))
           call extract_optstr(cmdline, "-", "GR",  1,count,tstr(2))
           call extract_optstr(cmdline, "-", "RC",  1,count,tstr(3))
           call extract_optstr(cmdline, "-", "SZ",  1,count,tstr(4))

           if(len_trim(tstr(1)) .eq. 0) then
             write(*,fmt="(A)") "MDPSCU error: the potential name is missed"
             write(*,fmt="(A, A, A)") 'Usage: ', gm_ExeName(1:len_trim(gm_ExeName)), ' -lib name [options]'
             write(*,fmt="(A, A, A)") '-lib name,  indicates the name of the potential for which the tables'
             write(*,fmt="(A, A, A)") '            to be generated, with the table to be saved in name.pair'
             write(*,fmt="(A, A, A)") '            and pname.embd'
             write(*,fmt="(A, A, A)") '[options] include:'
             write(*,fmt="(A, A, A)") '-gr psname, the subname of potential'
             write(*,fmt="(A, A, A)") '-rc rcut,   the cutoff range of the potential (in A)'
             write(*,fmt="(A, A, A)") '-sz nts,    the number of points in the potential table'
             stop
           else
              SimBox0%Potlibname = tstr(1)(1:len_trim(tstr(1)) )
           end if
           SimBox0%PotSublibname = tstr(2)(1:len_trim(tstr(2)) )

           if(len_trim(tstr(3)) .gt. 0) then
              CtrlParam%RU = DRSTR(tstr(3))*CP_A2CM
           else
               CtrlParam%RU = 10.D0*CP_A2CM
           end if

           if(len_trim(tstr(4)) .gt. 0) then
              CtrlParam%NUMFTABR = ISTR(tstr(4))
           else
              CtrlParam%NUMFTABR = 10000
           end if
           CtrlParam%NUMFTABE = CtrlParam%NUMFTABR

      !*** to initialize the force-potential table
            if(present(FORCETABLE)) then
              if(.not. present(POTTYPE)) then
                 write(*,fmt="(' MDPSCU Error: type of used-supplied potential is not given')")
                 write(*,fmt="('               supported optential include: ')")
                 do I=1, size(PKW_POTTYPELIST)
                    write(*,fmt="('               ', A)") '"'//PKW_POTTYPELIST(I)(1:len_trim(PKW_POTTYPELIST(I)))//'"'
                 end do
                 write(*,fmt="('               check the code calling Generator_Main')")
                 write(*,fmt="('               Process to be stopped')")
                 stop
              end if
              SimBox0%POTTYPE = POTTYPE(1:len_trim(POTTYPE))
              call Export_External_ForceLib(SimBox0, CtrlParam, FORCETABLE)
            else
              call Export_Internal_ForceLib(SimBox=SimBox0, CtrlParam=CtrlParam)
            end if

      return
  end subroutine Generator_Main
 !****************************************************************************************
 end module MD_FORCETABLE_GENERATOR

!****************************************************************************************
 Program FORCETABLE_GENERATOR_Main
 use MD_FORCETABLE_GENERATOR
  !---- If user-defined potential to be generated, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
       use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
  !---- Or:
       use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table
  implicit none
  integer::numprocs=1, procid=0

       !---- If user-define potential to be generated, modify the code as following examplee
       !     call Generator_Main(FORCETABLE=Reg1, POTTYPE="EAM_TYPE")
       !---- Or:
       !     call Generator_Main(FORCETABLE=Reg2, POTTYPE="FS_TYPE")

       !--- If internal potential to be generated, call Generator_Main as following
            call Generator_Main()


       STOP "--- FORCETABLE_GENERATOR"
       stop
  end  program FORCETABLE_GENERATOR_Main
!****************************************************************************************
