  module MD_ForceLib_Factory_GPU
  !***  DESCRIPTION:
  !     This module provides a interface to register to diffrent force lib
  !     that to be usd in simulation or analysis
  !
  !     DEPENDENCE: MD_Forceclass_Register_GPU
  !                 MD_ForceLib_Register
  !
  !     SEE ALSO:
  !
  !    Written by HOU Qing, May, 2014
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_TYPEDEF_ForceTable
   use MD_Forceclass_Register_GPU
   implicit none

  contains

  !****************************************************************************************
  subroutine Load_ForceLib(SimBox, CtrlParam, Forceclass)
  !***  DESCRIPTION: to load a the force lib. A force class to be resitered
  !                  with given type. Then the force table to be created
  !     INPUT:       SimBox:     the simulation box
  !                  CtrlParam:  the control parameters
  !
  !     OUPUT:       Forceclass: the force class, with its routine pointers of
  !                              force calculation assigned and force table created
  !

  use MD_Forcelib_Register, only:Register_Internal_ForceLib
  implicit none
  type(SimMDBox)                       ::SimBox
  type(SimMDCtrl),        intent(in)   ::CtrlParam
  type(MDForceClassGPU),  intent(inout)::Forceclass

  !$$--- local varibales
   character*12::pottype

       !$$--- reset the force class registered before
        call Release_ForceTable(Forceclass%ForceTable)

       !$$--- register the force class first
       pottype = SimBox%pottype(1:len_trim(SimBox%pottype))
       call UpCase(pottype)
       call Register_ForceClass(pottype, Forceclass)

       !$$--- register the table
       call Register_Internal_ForceLib(SimBox, CtrlParam, Forceclass%ForceTable)
       if(ForceClass%ForceTable%PotType .ne. pottype) then
          pottype = ForceClass%ForceTable%PotType
          SimBox%pottype = pottype
          call Register_ForceClass(pottype, Forceclass)
        end if

       return
   end subroutine Load_ForceLib
  !****************************************************************************************

  !****************************************************************************************
  subroutine RegUser_ForceLib(SimBox, CtrlParam, Forceclass, REGISTER)
  !***  DESCRIPTION: to load a the force lib
  !     INPUT: pottype,   the type of the potential
  !                       such as: FS_TYPE, EAM_TYPE etc
  !            libname,   the libname we have
  !
  !
  use MD_Forcelib_Register, only:Register_External_ForceLib

  implicit none
   !--- dummy variable
   optional::REGISTER
   external::REGISTER
   interface
     subroutine REGISTER(SimBox, CtrlParam, FTable, printout)
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
     end subroutine REGISTER
   end interface
   type(SimMDBox),        intent(in)   ::SimBox
   type(SimMDCtrl),       intent(in)   ::CtrlParam
   type(MDForceClassGPU), intent(inout)::Forceclass


   !$$--- local varibales
   character*12::pottype

       !$$--- reset the force class registered before
        call Release_ForceTable(Forceclass%ForceTable)

       !$$--- register the force class first
       pottype = SimBox%pottype(1:len_trim(SimBox%pottype))
       call UpCase(pottype)
       call Register_ForceClass(pottype, Forceclass)

       !$$--- register the table
       if(present(REGISTER)) then
          call Register_External_ForceLib(SimBox, CtrlParam, Forceclass%ForceTable, REGISTER)
       end if   

       return
   end subroutine RegUser_ForceLib
  !****************************************************************************************

 !****************************************************************************************
  subroutine Clear_ForceLib(Forceclass)
  !***  DESCRIPTION: to releae the memory allocated on load a the force lib.
  !     INPUT:
  !
  !     OUPUT:       Forceclass: the force class, with its  force tables deallolcated
  !

  implicit none
  type(MDForceClassGPU), intent(inout) ::Forceclass

  !$$--- local varibales
   character*12::pottype

        call Forceclass%pClrForcetable()
       return
   end subroutine Clear_ForceLib
  !****************************************************************************************
  end module MD_ForceLib_Factory_GPU
  !****************************************************************************************
