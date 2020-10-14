  module MD_ForceLib_Register
  !***  DESCRIPTION:
  !     This module provides a interface to register to diffrent force lib
  !     that to be usd in simulation or analysis. This module will be used
  !     by MD_ForceLib_Factory_CPU.F90 and MD_ForceLib_Factory_GPU.F90
  !
  !     DEPENDENCE: MD_TYPEDEF_ForceTable
  !
  !     SEE ALSO:
  !             MD_ForceLib_Factory_CPU.F90
  !             MD_ForceLib_Factory_GPU.F90
  !
  !    Written by HOU Qing, May, 2014
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none

  !--- interface to the register force routine -------------------
   private REGISTER_TABLE
   abstract interface
     subroutine REGISTER_TABLE(SimBox, CtrlParam, FTable, printout)
     use MD_TYPEDEF_SimMDBox
     use MD_TYPEDEF_SimMDCtrl
     use MD_TYPEDEF_ForceTable
     implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in)   ::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(inout)::FTable
      integer,            optional     ::printout
      
     end subroutine REGISTER_TABLE
   end interface

 !--- END INTERFACE --------------------------------

 !--- private moudle menbers:
     procedure(REGISTER_TABLE),  pointer, private::m_pRegForcetable=>null()

    !-------------------------------------------
     private:: Register_ForceLib_FS
    !-------------------------------------------
     private:: Register_ForceLib_EAM
    !-------------------------------------------
     public::  Register_External_ForceLib
    !-------------------------------------------
     public::  Register_Internal_ForceLib
    !-------------------------------------------
  contains

  !****************************************************************************************
  subroutine Register_ForceLib_FS(SimBox, CtrlParam, ForceTable, ERR)
  !***  DESCRIPTION: to register the FS force lib
  !     INPUT: libname,   the libname we have
  !            SimBox,    the simulation box
  !            CtrlParam, the control parameters
  !     OUTPUT:
  !            FTable, the force tabled to be registered
  !            ERR,    the statu of regiter
  !
  !
  use  MD_TYPEDEF_ForceTable, only:MDForceTable
  !$$--- the potential for Zr-Zr of Acland
  use  EM_TB_ForceTable_AWB_ZR, only:Reg_AWB_ZR=>Register_Interaction_Table

  !$$--- the potential for W-W and He-W of Wang Jun
  use  EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg_Wang_HeW=>Register_Interaction_Table

  !$$--- the potential of Cleri and Rosato, Phys.Rev.B48, (1993) 22
  use  EM_TB_ForceTable_CR_2014, only:Reg_CR=>Register_Interaction_Table

  !$$--- the potential of G.J.Ackland and V.Vitek, Phys.Rev.B41, (1990) 10324
  use  EM_TB_ForceTable_Ackland_Vitek_PRB41_10324, only:Reg_AV_PRB10324=>Register_Interaction_Table

  !$$--- the potential for Ti-Ti and He-Ti of Wang Jun
  use  EM_TB_ForcetTable_WangJun_Ti_He_2007, only:Reg_Wang_HeTi=>Register_Interaction_Table

  !
  implicit none
    type(SimMDBox),      intent(in)   ::SimBox
    type(SimMDCtrl),     intent(in)   ::CtrlParam
    type(MDForceTable),  intent(inout)::ForceTable
    integer,             intent(out)  ::ERR
    !--- local varibales
    character*16::tstr
    character*2::symb
    character*128::libname
    integer::I

       ForceTable%PotLibname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       if(len_trim(SimBox%potsublibname) .gt. 0) then
          ForceTable%PotSubLibname = SimBox%potsublibname(1:len_trim(SimBox%potsublibname))
       else
          ForceTable%PotSubLibname = ""
       end if

       libname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       call UpCase(libname)
       m_pRegForcetable=>null()
       ERR = 0
       select case(libname(1:len_trim(libname)))
              !--------------------------------------------------------------------
              case("EM_TB_AWB_PHIL_MAG_A71_1995_553")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "ZR") then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is for Zr-Zr interactions,')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_AWB_ZR

              !--------------------------------------------------------------------
              case("EM_TB_WANGJUN_W-HE_2010")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "HE" .and. &
                         symb(1:len_trim(symb)) .ne. "W") then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is for He-He, W-W, and He-W interactions,')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_Wang_HeW

              !--------------------------------------------------------------------
              case("EM_TB_CLERI_ROSATO_PRB48_22")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "AL" .and. &
                         symb(1:len_trim(symb)) .ne. "AU" .and. &
                         symb(1:len_trim(symb)) .ne. "TI"  ) then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is available only for Al-Al,Au-Au, Ti-Ti...')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_CR
              !--------------------------------------------------------------------

              case("EM_TB_ACKLAND_VITEK_PRB41_10324")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "CU" .and. &
                         symb(1:len_trim(symb)) .ne. "AU" ) then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is available only for Au-Au, Cu-Cu, Au-Cu...')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_AV_PRB10324

              case("EM_TB_WANGJUN_TI-HE_2007")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "TI" .and. &
                         symb(1:len_trim(symb)) .ne. "HE" ) then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is available only for Ti-Ti, He-He')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_Wang_HeTi

              !--------------------------------------------------------------------
              case default
                   if(len_trim(libname).le.0) then
                      write(*,fmt="(' MDPSCU Error: empty potential libname')")
                      write(*,fmt="('               check input boxfile in  &POTSUBCTL subsection')")
                      write(*,fmt="('               Process to be stopped')")
                      stop
                    else
                      ERR = 1
                      return
                    end if
       end select

       call m_pRegForcetable(SimBox, CtrlParam, ForceTable)
       return
   end subroutine Register_ForceLib_FS
  !****************************************************************************************

  !****************************************************************************************
  subroutine Register_ForceLib_EAM(SimBox, CtrlParam, ForceTable, ERR)
  !***  DESCRIPTION: to register the internal EAM force lib
  !     INPUT: libname,   the libname we have
  !            SimBox,    the simulation box
  !            CtrlParam, the control parameters
  !     OUTPUT:
  !            FTable, the force tabled to be registered
  !            ERR,    the statu of regiter
  !
  !
  use MD_TYPEDEF_ForceTable, only:MDForceTable
  !$$--- the potential for W-W of Marinica et al, JPCM25(2013)395502
  use EAM_WW_ForceTable_Marinica_JPCM25_2013, only:Reg_Marinica_WW=>Register_Interaction_Table
  !$$--- the potential for WHeH of Bonny, Grigorev and Terentyev, JPCM26(2014)485001
   use EAM_WHeH_ForceTable_Bonny_JPCM26_2014,  only:Reg_Bonny_WHeH=>Register_Interaction_Table
  !$$--- the potential for ZrZr of Mendelev and Acland, Phil. Mag. Lett. 87(2007)349-359
   use EAM_ZR_ForceTable_Mendelev_Phil_Mag_L87_2007,  only:Reg_Memdelev_ZrZr=>Register_Interaction_Table

  implicit none
    type(SimMDBox),      intent(in)::SimBox
    type(SimMDCtrl),     intent(in)::CtrlParam
    type(MDForceTable),  intent(inout)::ForceTable
    integer,             intent(out)::ERR
    !--- local varibales
    character*16::tstr
    character*2::symb
    character*128::libname
    integer I


       ForceTable%PotLibname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       if(len_trim(SimBox%potsublibname) .gt. 0) then
          ForceTable%PotSubLibname = SimBox%potsublibname(1:len_trim(SimBox%potsublibname))
       else
          ForceTable%PotSubLibname = ""
       end if

       libname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       call UpCase(libname)
       m_pRegForcetable=>null()
       ERR = 0
       select case(libname(1:len_trim(libname)))
              case("EM_TB_AWB_PHIL_MAG_A71_1995_553",&
                   "EM_TB_WANGJUN_W-HE_2010",        &
                   "EM_TB_CLERI_ROSATO_PRB48_22",    &
                   "EM_TB_ACKLAND_VITEK_PRB41_10324")
                   call Register_ForceLib_FS(SimBox, CtrlParam, ForceTable, ERR)
                   return

              case("EAM_WW_MARINICA_JPCM25_2013")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "W" ) then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is available only for W-W')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_Marinica_WW;

              case("EAM_WHEH_BONNY_JPCM26_2014")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "W" .and. &
                         symb(1:len_trim(symb)) .ne. "HE".and. &
                         symb(1:len_trim(symb)) .ne. "H") then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is available only for W-W, W-He, W-H. He-He, He-H, H-H')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_Bonny_WHeH;

              case("EAM_ZR_MENDELEV_PHIL_MAG_L87_2007")
                   do I=1, SimBox%NGROUP
                      tstr = SimBox%SYMB(I)(1:len_trim(SimBox%SYMB(I)))
                      call UpCase(tstr)
                      tstr = adjustl(tstr)
                      call Extract_SubSymb(tstr, symb)
                      if(symb(1:len_trim(symb)) .ne. "ZR") then
                         write(*,fmt="(' MDPSCU Warning: potential lib ', A)")libname(1:len_trim(libname))
                         write(*,fmt="('                 is available only for Zr-Zr')")
                         if(len_trim(tstr) .gt. 0) then
                            write(*,fmt="('                 while atom ',A, ' exist in the present box.')") tstr(1:len_trim(tstr))
                         else
                            write(*,fmt="('                 while unkown atom exist in the present box.')")
                         end if
                         write(*,fmt="('                 please check the input box file.')")
                         write(*,fmt="('                 ')")
                         call ONWARNING(gm_OnWarning)
                      end if
                   end do
                   m_pRegForcetable => Reg_Memdelev_ZrZr;

              case default
                   if(len_trim(libname).le.0) then
                      write(*,fmt="(' MDPSCU Error: empty potential libname')")
                      write(*,fmt="('               check input boxfile in  &POTSUBCTL subsection')")
                      write(*,fmt="('               Process to be stopped')")
                      stop
                    else
                      ERR = 1
                      return
                    end if
       end select

       call m_pRegForcetable(SimBox, CtrlParam, ForceTable)
       return
   end subroutine Register_ForceLib_EAM
  !****************************************************************************************

  !****************************************************************************************
  subroutine Register_Internal_ForceLib(SimBox, CtrlParam, FTable)
  !***  DESCRIPTION: to load a the force lib. A force class to be resitered
  !                  with given type. Then the force table to be created
  !     INPUT:       SimBox:     the simulation box
  !                  CtrlParam:  the control parameters
  !
  !     OUPUT:       Forceclass: the force class, with its routine pointers of
  !                              force calculation assigned and force table created
  !
  use MD_TYPEDEF_ForceTable, only: MDForceTable, Register_Imported_ForceTable, Release_ForceTable, Export_ForceTable
  use NIST_ForceTable,       only:NIST_Register_Interaction_Table
  implicit none
  type(SimMDBox),    intent(in)::SimBox
  type(SimMDCtrl)              ::CtrlParam
  type(MDForceTable)           ::FTable

  !$$--- local varibales
   character*12::pottype
   integer::ERR
   logical::EX1, EX2

       !$$--- register the force class first
       pottype = SimBox%pottype(1:len_trim(SimBox%pottype))
       call UpCase(pottype)

       !$$--- register the force lib and create the force table
       FTable%PotType = pottype(1:len_trim(pottype))
       ERR = 0
       select case(pottype(1:len_trim(pottype)))
              case("FS_TYPE")
                   call Register_ForceLib_FS(SimBox, CtrlParam, FTable, ERR)

              case("EAM_TYPE")
                   call Register_ForceLib_EAM(SimBox, CtrlParam, FTable, ERR)

              case default
                   ERR = 1
       end select

       if(ERR .gt. 0) then
          !$$--- to find out if the force table could be imprted from files
          inquire(FILE=SimBox%potlibname(1:len_trim(SimBox%potlibname))//'.pair', exist=EX1)
          inquire(FILE=SimBox%potlibname(1:len_trim(SimBox%potlibname))//'.embd', exist=EX2)
          if(.not.EX1 .or. .not.EX2) then

            !$$--- we further try if NIST the potential is NIST forcetable
            FTable%PotType = "EAM_TYPE"
            write(*,fmt="(A)") ' MDPSCU Message: to transfer NIST potential table '//SimBox%potlibname(1:len_trim(SimBox%potlibname))
            call NIST_Register_Interaction_Table(SimBox, CtrlParam, FTable)
            call Release_ForceTable(FTable)
         end if

         call Register_Imported_ForceTable(SimBox, CtrlParam, FTable)
         !$$--- it is better not to overwrite the table in file
         CtrlParam%OUTPUTFTAB = 0
       end if

       if(CtrlParam%OUTPUTFTAB .gt. 0) then
          call Export_ForceTable(SimBox%potlibname(1:len_trim(SimBox%potlibname)),FTable)
       end if
       return
   end subroutine Register_Internal_ForceLib
  !****************************************************************************************

  !****************************************************************************************
  subroutine Register_External_ForceLib(SimBox, CtrlParam, FTable, Register)
  !***  DESCRIPTION: to register an external force table
  !     INPUT:
  !
  !
  use MD_TYPEDEF_ForceTable, only:MDForceTable, Export_ForceTable
  implicit none
   !--- dummy variable
   external::Register
   !--- interface for external register
   interface
     subroutine Register(SimBox, CtrlParam, FTable, printout)
     use MD_TYPEDEF_SimMDBox
     use MD_TYPEDEF_SimMDCtrl
     use MD_TYPEDEF_ForceTable
     implicit none
      !--- dummy vaiables
      type(SimMDBox),      intent(in)   ::SimBox
      type(SimMDCtrl),     intent(in)   ::CtrlParam
      type(MDForceTable),  intent(inout)::FTable
      integer,             optional     ::printout
     end subroutine Register
   end interface
   !--- end interface
   type(SimMDBox),     intent(in)   ::SimBox
   type(SimMDCtrl),    intent(in)   ::CtrlParam
   type(MDForceTable), intent(inout)::FTable

  !$$--- local varibales
   character*12::pottype

       !$$--- register the force class first
       pottype = SimBox%pottype(1:len_trim(SimBox%pottype))
       call UpCase(pottype)

       !$$--- register the force lib and create the force table
       FTable%PotType = pottype(1:len_trim(pottype))

       m_pRegForcetable =>Register
       call m_pRegForcetable(SimBox, CtrlParam, FTable)
       if(CtrlParam%OUTPUTFTAB .gt. 0) then
          call Export_ForceTable(SimBox%potlibname(1:len_trim(SimBox%potlibname)),FTable)
       end if
       return
   end subroutine Register_External_ForceLib
  !****************************************************************************************

  !****************************************************************************************
  subroutine Export_Internal_ForceLib(SimBox, CtrlParam)
  !***  DESCRIPTION: to export the internal force lib
  !     INPUT: libname,   the libname we have
  !            SimBox,    the simulation box
  !            CtrlParam, the control parameters
  !     OUTPUT:
  !            FTable, the force tabled to be registered
  !
  !
   use MD_TYPEDEF_ForceTable, only:MDForceTable, Generate_Lib_ForceTable, Release_ForceTable
   implicit none
   !--- dummy variable
    type(SimMDBox)  ::SimBox
    type(SimMDCtrl) ::CtrlParam
    !--- local varibales
    character*128::libname
    type(MDForceTable)::FTable


       FTable%PotLibname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       if(len_trim(SimBox%potsublibname) .gt. 0) then
          FTable%PotSubLibname = SimBox%potsublibname(1:len_trim(SimBox%potsublibname))
       else
          FTable%PotSubLibname = ""
       end if

       libname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       call UpCase(libname)
       CtrlParam%OUTPUTFTAB = 0
       FTable%Rmax = maxval(CtrlParam%RU)
       select case(libname(1:len_trim(libname)))
              case("EM_TB_AWB_PHIL_MAG_A71_1995_553", &
                   "EM_TB_WANGJUN_W-HE_2010",         &
                   "EM_TB_CLERI_ROSATO_PRB48_22",     &
                   "EM_TB_ACKLAND_VITEK_PRB41_10324" )
                   SimBox%pottype = "FS_TYPE"
                   call Register_Internal_ForceLib(SimBox, CtrlParam, FTable)
                   call Generate_Lib_ForceTable(FTable, SimBox%potlibname)

              case("EAM_WW_MARINICA_JPCM25_2013", &
                   "EAM_WHEH_BONNY_JPCM26_2014",  &
                   "EAM_ZR_MENDELEV_PHIL_MAG_L87_2007"   )
                   SimBox%pottype = "EAM_TYPE"
                   call Register_Internal_ForceLib(SimBox, CtrlParam, FTable)
                   call Generate_Lib_ForceTable(FTable, SimBox%potlibname)

             case default
                  SimBox%pottype = "EAM_TYPE"
                  call Register_Internal_ForceLib(SimBox, CtrlParam, FTable)
       end select
       call Release_ForceTable(FTable, keepreglist=0)

       return
   end subroutine Export_Internal_ForceLib
  !****************************************************************************************
  !****************************************************************************************
  subroutine Export_External_ForceLib(SimBox, CtrlParam, Register)
  !***  DESCRIPTION: to export the external force lib
  !     INPUT: SimBox,    the simulation box
  !            CtrlParam, the control parameters
  !     OUTPUT:
  !            FTable, the force tabled to be registered
  !
  !
   use MD_TYPEDEF_ForceTable, only:MDForceTable, Generate_Lib_ForceTable, Release_ForceTable
  implicit none
   !--- dummy variable
   !--- interface for external register
   interface
     subroutine Register(tSimBox, tCtrlParam, tFTable, printout)
     use MD_TYPEDEF_SimMDBox
     use MD_TYPEDEF_SimMDCtrl
     use MD_TYPEDEF_ForceTable
     implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in)  ::tSimBox
      type(SimMDCtrl),    intent(in)  ::tCtrlParam
      type(MDForceTable), intent(inout)::tFTable
      integer,optional::printout
     end subroutine Register
   end interface
   !--- end interface
   type(SimMDBox)    ::SimBox
   type(SimMDCtrl)   ::CtrlParam
    !--- local varibales
    character*128::libname
    type(MDForceTable)::FTable


       FTable%PotLibname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       if(len_trim(SimBox%potsublibname) .gt. 0) then
          FTable%PotSubLibname = SimBox%potsublibname(1:len_trim(SimBox%potsublibname))
       else
          FTable%PotSubLibname = ""
       end if

       libname = SimBox%potlibname(1:len_trim(SimBox%potlibname))
       call UpCase(libname)

       CtrlParam%OUTPUTFTAB = 0
       call Register_External_ForceLib(SimBox, CtrlParam, FTable, Register)
       call Generate_Lib_ForceTable(FTable, SimBox%potlibname)
       call Release_ForceTable(FTable, keepreglist=0)
       return
   end subroutine Export_External_ForceLib
  !****************************************************************************************

  end module MD_ForceLib_Register
  !****************************************************************************************
