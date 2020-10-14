  module MD_ForceClass_Register_GPU
  !***  DESCRIPTION:
  !     This module provides interfaces to force calculations using different TYPE of
  !     potentials. In this module, no concrete potential is implmented, but the type
  !     potential is registered
  !
  !**** HISTORY:     2014-05 (HOU Qing), created the first version
  !
  !                  2020-07-05(HOU Qing),
  !                          In routine Register_ForceClass, added an class type "EXT_TYPE" in addition to
  !                          "FS_TYPE, and "EAM_TYPE".
  !                           If a EXT_TYPE potential is used,  the pointers in MDForceClassGPU are nullified.
  !                           The force and epot calculation will be conducted in MDExtForceList, in which
  !                           the external force and epot calculation routines are added by calling
  !                           AddExtForce_ForceClass.
  !     
  !
  !________________________________________________________________________________________________________
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_TYPEDEF_ForceTable

   implicit none

   !------------------------------------------------------------------------
   !--- interface to the initalization of force table
   private::INITIALIZE_FORCETABLE
   abstract interface
     subroutine INITIALIZE_FORCETABLE(SimBox, CtrlParam, FTable, RelaseTable, MULTIBOX)
     !***  PURPOSE:   to intialize the parameters to be used in force calculation
     !                can copy them in device
     !     INPUT:     SimBox: the simulation box
     !                CtrlParam: the control parameters
     !                FTable, the force tabled to be used
     !                RelaseTable, indicator to indicate if the force table will be deallocate after
     !                             copyin devices
     !                MULTIBOX, indicating if use multiple box
     !
     !     OUTPUT:   the working spaces allocated
     !
     use MD_TYPEDEF_SimMDBox
     use MD_TYPEDEF_SimMDCtrl
     use MD_TYPEDEF_ForceTable
     implicit none
       !--- dummy vaiables
      type(SimMDBox),     intent(inout):: SimBox
      type(SimMDCtrl),    intent(in)   :: CtrlParam
      type(MDForceTable), intent(in)   :: FTable
      integer,            optional     :: RelaseTable, MULTIBOX

     end subroutine INITIALIZE_FORCETABLE
   end interface
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  !--- interface to clearfing of force table
  private::CLEAR_FORCETABLE
   abstract interface
     subroutine CLEAR_FORCETABLE()
     !***  PURPOSE:   to clear memories allocated in force calculation
     !     INPUT:     none
     !
     implicit none

     end subroutine CLEAR_FORCETABLE
   end interface
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  !--- interface to force calculationm
  private::CALFORCE
  abstract interface
   subroutine CALFORCE(SimBox, CtrlParam)
   !***  PURPOSE:   to calculate the forces
   !
   !     INPUT:     SimBox,    the simulation box
   !                CtrlParam: the control parameters
   !
   !     OUTPUT     SimBox,    the simulation box with force updated
   !
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam
    end subroutine CALFORCE
  end interface
 !------------------------------------------------------------------------

 !------------------------------------------------------------------------
 !--- interface to external force calculationm
  private::CALFORCE_EXT
  abstract interface
   subroutine CALFORCE_EXT(SimBox, CtrlParam)
   !***  PURPOSE:   to calculate the forces
   !
   !     INPUT:     SimBox,    the simulation box
   !                CtrlParam: the control parameters
   !
   !     OUTPUT     SimBox,    the simulation box with force updated
   !
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout)::SimBox
      type(SimMDCtrl),              intent(in)   ::CtrlParam
    end subroutine CALFORCE_EXT
 end interface
 !------------------------------------------------------------------------

 !------------------------------------------------------------------------
 !--- interface to potential calculation
 private::CALEPOT
 abstract interface
   subroutine CALEPOT(SimBox, CtrlParam)
   !***  PURPOSE:   to calculate the potential for each atom, and copy the data into dm_EPOT
   !
   !     INPUT:     SimBox,    the simulation box
   !                CtrlParam: the control parameters
   !
   !     OUTPUT:   m_EPOT
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam
   end subroutine CALEPOT
 end interface
 !------------------------------------------------------------------------

 !------------------------------------------------------------------------
 !--- interface to external potential calculation
 private::CALEPOT_EXT
 abstract interface 
   subroutine CALEPOT_EXT(SimBox, CtrlParam)
   !***  PURPOSE:   to calculate the potential for each atom, and copy the data into dm_EPOT
   !
   !     INPUT:     SimBox,    the simulation box
   !                CtrlParam: the control parameters
   !
   !     OUTPUT:   m_EPOT
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout)::SimBox
      type(SimMDCtrl),              intent(in)   ::CtrlParam
   end subroutine CALEPOT_EXT
 end interface
 !------------------------------------------------------------------------

 !------------------------------------------------------------------------
 !--- interface to electronic density calculation
 private::CALEDEN
 abstract interface
   subroutine CALEDEN(SimBox, CtrlParam)
   !***  PURPOSE:   to calculate the electron densities on atoms only.
   !                Could be used for other applications, for exmaples
   !                in local stress caculations
   !
   !     INPUT:     SimBox,    the simulation box
   !                CtrlParam: the control parameters
   !     OUTPUT     dm_DEM,    the simulation box with force updated
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam
   end subroutine CALEDEN
 end interface
 !------------------------------------------------------------------------


 !------------------------------------------------------------------------
 !--- interface to force and pressure calculation
 private::CALPTENSOR
 abstract interface
   subroutine CALPTENSOR(SimBox, CtrlParam)
   !***  PURPOSE:   to calculate the forces.Virial tensor are also calculated
   !
   !     INPUT:     SimBox,    the simulation box
   !                CtrlParam: the control parameters
   !
   !     OUTPUT     SimBox,    the simulation box with force updated
   !
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl

   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam
   end subroutine CALPTENSOR
 end interface
 !------------------------------------------------------------------------

 !------------------------------------------------------------------------
 !--- interface to atomic virial stress tensor calculation
 private::CALATOMICSTRESS
 abstract interface
   subroutine CALATOMICSTRESS(IDEV, dAVP)
   !***  PURPOSE:   to calculate the stomic virial stress tensor only.
   !                Could be used for other applications, for exmaples
   !                in local stress caculations
   !
   !     INPUT:     IDEV,  the ID of device
   !
   !     OUTPUT     dAVP:  the array of atomic stress, should be allocated before allocated the routione
   !
   use MD_CONSTANTS
    implicit none
      !--- dummy vaiables
      integer,                              intent(in) ::IDEV
      real(KINDDF), device, dimension(:,:), intent(out)::dAVP
   end subroutine CALATOMICSTRESS
 end interface
 !------------------------------------------------------------------------

 !--- END INTERFACES --------------------------------

 !--- type definition of external force list
      private::MDExtForce
      type MDExtForce
           procedure(CALEPOT),              pointer, nopass::pEpot    =>null()
           procedure(CALFORCE),             pointer, nopass::pForce   =>null()
           procedure(CALEPOT_EXT),          pointer, nopass::pEpotEx  =>null()
           procedure(CALFORCE_EXT),         pointer, nopass::pForceEx =>null()
      end type MDExtForce

      private::MDExtForceList
      type MDExtForceList
           character(len=16)                       :: Tag  = ""
           type(MDExtForce),                pointer:: this    =>null()
           type(MDExtForceList),            pointer:: next    =>null()
      end type MDExtForceList

 !--- type definition of FORCE CLASS
      type::MDForceClassGPU
           procedure(INITIALIZE_FORCETABLE), pointer, nopass::pIniForcetable =>null()
           procedure(CLEAR_FORCETABLE),      pointer, nopass::pClrForcetable =>null()
           procedure(CALFORCE),              pointer, nopass::pCalForce      =>null()
           procedure(CALEPOT),               pointer, nopass::pCalEpot0      =>null()
           procedure(CALEPOT),               pointer, nopass::pCalEpot       =>null()
           procedure(CALEDEN),               pointer, nopass::pCalEDen       =>null()
           procedure(CALPTENSOR),            pointer, nopass::pCalPTensor    =>null()
           procedure(CALATOMICSTRESS),       pointer, nopass::pCalAVStress   =>null()
           type(MDForceTable)::                               ForceTable
           type(MDExtForceList)::                             ExtForceList
      end type MDForceClassGPU
      !--- public accessible members
      type(MDForceClassGPU)::gm_ForceClass
  !**********************************************************************
  !--- the Interfaces to external calls
  !    
  !---------------------------------------------------------
      private::  Add_ExtForce,             &
                 Clear_ExtForce,           &
                 AddExtForce_ForceClass0,  &
                 AddExtForce_ForceClass1,  & 
                 ClearExtForce_ForceClass0,&
                 ClearExtForce_ForceClass1

      public::   AddExtForce_ForceClass
      interface  AddExtForce_ForceClass
                 module procedure AddExtForce_ForceClass0
                 module procedure AddExtForce_ForceClass1
      end interface  AddExtForce_ForceClass

      public::   ClearExtForce_ForceClass
      interface  ClearExtForce_ForceClass
                 module procedure ClearExtForce_ForceClass0
                 module procedure ClearExtForce_ForceClass1
      end interface  ClearExtForce_ForceClass

  !---------------------------------------------------------
      private::  CalEpot_ForceClass0,   &
                 CalEpot_ForceClass1
      public::   CalEpot_ForceClass
      interface  CalEpot_ForceClass
                 module procedure CalEpot_ForceClass0
                 module procedure CalEpot_ForceClass1
      end interface CalEpot_ForceClass

  !---------------------------------------------------------
      private::  CalForce_ForceClass0,   &
                 CalForce_ForceClass1
      public::   CalForce_ForceClass
      interface  CalForce_ForceClass
                 module procedure CalForce_ForceClass0
                 module procedure CalForce_ForceClass1
      end interface CalForce_ForceClass

  !---------------------------------------------------------
      private::  CalPTensor_ForceClass0,   &
                 CalPTensor_ForceClass1
      public::   CalPTensor_ForceClass
      interface  CalPTensor_ForceClass
                 module procedure CalPTensor_ForceClass0
                 module procedure CalPTensor_ForceClass1
      end interface CalPTensor_ForceClass

  !---------------------------------------------------------
      public::   CopyPointer_ForceClass
  !---------------------------------------------------------
      private::  Cumulate_ExtForce0,   &
                 Cumulate_ExtForce1
      public::   Cumulate_ExtForce
      interface  Cumulate_ExtForce
                 module procedure Cumulate_ExtForce0
                 module procedure Cumulate_ExtForce1
      end interface Cumulate_ExtForce

  !---------------------------------------------------------
      private::  Cumulate_ExtEpot0,   &
                 Cumulate_ExtEpot1
      public::   Cumulate_ExtEpot
      interface  Cumulate_ExtEpot
                 module procedure Cumulate_ExtEpot0
                 module procedure Cumulate_ExtEpot1
      end interface Cumulate_ExtEpot

  !---------------------------------------------------------
      public::   IfInit_ForceClass
      private::  Init_Forcetable0,   &
                 Init_Forcetable1,   &
                 Init_Forcetable0_a, &
                 Init_Forcetable1_a
      public::   Init_Forcetable_Dev
      interface  Init_Forcetable_Dev
                 module procedure Init_Forcetable0
                 module procedure Init_Forcetable1
                 module procedure Init_Forcetable0_a
                 module procedure Init_Forcetable1_a
       
      end interface Init_Forcetable_Dev

  !---------------------------------------------------------
      public::   Register_ForceClass

  !---------------------------------------------------------
      private::  UpdateEpot_ForceClass0,   &
                 UpdateEpot_ForceClass1
      public::   UpdateEpot_ForceClass
      interface  UpdateEpot_ForceClass
                 module procedure UpdateEpot_ForceClass0
                 module procedure UpdateEpot_ForceClass1
      end interface UpdateEpot_ForceClass
      
  !---------------------------------------------------------
  contains

  !****************************************************************************************
  subroutine Register_ForceClass(classname, ForceClass)
  !***  DESCRIPTION:   to regist a force class
  !     INPUT:         classname,   the name of the force class
  !
  !
 !$$-------FS moudle-----------------------
  use  MD_FS_Force_Table_GPU, only:FS_INIT    =>INITIALIZE_FS_Force_Table_DEV,     &
                                   FS_CALF    =>CALFORCE_FS_Force_Table2A_DEV,     &
                                   FS_CALE0   =>UpdateEPOT_FS_Force_Table2A_DEV,   &
                                   FS_CALE    =>CALEPOT_FS_Force_Table2A_DEV,      &
                                   FS_CALP    =>CALPTENSOR_FS_Force_Table2A_DEV,   &
                                   FS_CLR     =>Clear_FS_Force_Table_DEV,          &
                                   FS_CALD    =>CALDEN_FS_Force_Table2A_DEV,       &
                                   FS_CALAVS  =>Cal_FS_AtomicStressTensor_DEV

 !$$-------EAM moudle-----------------------
  use  MD_EAM_Force_Table_GPU, only:EAM_INIT  =>INITIALIZE_EAM_Force_Table_DEV,   &
                                    EAM_CALF  =>CALFORCE_EAM_Force_Table2A_DEV,   &
                                    EAM_CALE0 =>UpdateEPOT_EAM_Force_Table2A_DEV, &
                                    EAM_CALE  =>CALEPOT_EAM_Force_Table2A_DEV,    &
                                    EAM_CALP  =>CALPTENSOR_EAM_Force_Table2A_DEV, &
                                    EAM_CLR   =>Clear_EAM_Force_Table_DEV,        &
                                    EAM_CALD  =>CALDEN_EAM_Force_Table2A_DEV,     &
                                    EAM_CALAVS=>Cal_EAM_AtomicStressTensor_DEV


  !$$------------------------------
  implicit none
  character*(*)::classname
  type(MDForceClassGPU)::ForceClass

  !$$--- local varibales

       select case(classname(1:len_trim(classname)))
              case("FS_TYPE")
                   ForceClass%ForceTable%PotType = "FS_TYPE"
                   ForceClass%pIniForcetable     =>FS_INIT
                   ForceClass%pCalForce          =>FS_CALF
                   ForceClass%pCalEpot0          =>FS_CALE0
                   ForceClass%pCalEpot           =>FS_CALE
                   ForceClass%pCalPTensor        =>FS_CALP
                   ForceClass%pClrForcetable     =>FS_CLR
                   ForceClass%pCalEDen           =>FS_CALD
                   ForceClass%pCalAVStress       =>FS_CALAVS

              case("EAM_TYPE")
                   ForceClass%ForceTable%PotType = "EAM_TYPE"
                   ForceClass%pIniForcetable     =>EAM_INIT
                   ForceClass%pCalForce          =>EAM_CALF
                   ForceClass%pCalEpot0          =>EAM_CALE0
                   ForceClass%pCalEpot           =>EAM_CALE
                   ForceClass%pCalPTensor        =>EAM_CALP
                   ForceClass%pClrForcetable     =>EAM_CLR
                   ForceClass%pCalEDen           =>EAM_CALD
                   ForceClass%pCalAVStress       =>EAM_CALAVS

              case("EXT_TYPE")
                   ForceClass%ForceTable%PotType = "EXT_TYPE"
                   ForceClass%pIniForcetable     =>NULL()
                   ForceClass%pCalForce          =>NULL()
                   ForceClass%pCalEpot0          =>NULL()
                   ForceClass%pCalEpot           =>NULL()
                   ForceClass%pCalPTensor        =>NULL()
                   ForceClass%pClrForcetable     =>NULL()
                   ForceClass%pCalEDen           =>NULL()
                   ForceClass%pCalAVStress       =>NULL()

              case default
                   if(len_trim(classname).le.0) then
                      write(*,fmt="(' MDPSCU Error: empty type of potential', A16)")
                    else
                      write(*,fmt="(' MDPSCU Error: unsupported type of potential ', A16)") classname
                    end if
                    write(*,fmt="('               check input box file in  &POTSUBCTL subsection', A16)")
                    write(*,fmt="('   Process to be stopped', A16)")
                   stop
       end select

       return
   end subroutine Register_ForceClass
  !*********************************************************************************

  !*********************************************************************************
  subroutine CopyPointer_ForceClass(From, To)
  !***  DESCRIPTION:   to pass the function pointers in ForceClass From
  !                    to ForceClass To. The ForceClass From and To share
  !                    the same force table and core functions for force 
  !                    and energy calculation. But the ExtForceList
  !                    could be diffrent. This designe is for the conveniencec
  !                    of switching force calculations when different 
  !                    external forces, for example boost force, are used.     
  !                       
  !     INPUT:         From,  the source force class
  !
  !     OUTPUT:        To,    the target force class
  !

   implicit none
   !--- dummy arguments
   type(MDForceClassGPU), intent(in)  ::From 
   type(MDForceClassGPU), intent(out) ::To
   !---

       To%pIniForcetable => From%pIniForcetable
       To%pClrForcetable => From%pClrForcetable
       To%pCalForce      => From%pCalForce
       To%pCalEpot0      => From%pCalEpot0
       To%pCalEpot       => From%pCalEpot
       To%pCalEDen       => From%pCalEDen
       To%pCalPTensor    => From%pCalPTensor
       To%pCalAVStress   => From%pCalAVStress
       !--- we to not copy core Force table which has been
       !    copyin to device on registering forceclass From.
       !    The target force class may also has different external
       !    force list, we do not copy ExtForcelist 
       return
  end subroutine CopyPointer_ForceClass
  !*********************************************************************************
  
  !*********************************************************************************
  subroutine Init_Forcetable0(SimBox, CtrlParam, ForceClass, RelaseTable)
  implicit none
   !--- dummy vaiables
        type(SimMDBox),       intent(inout) :: SimBox
        type(SimMDCtrl),      intent(in)    :: CtrlParam
        type(MDForceClassGPU),intent(in)    :: ForceClass
        integer, optional                   :: RelaseTable
    !--- local
           
           if(ForceClass%ForceTable%PotType == "EXT_TYPE") return;
           !--- for EAM, FS type
           if(.not.associated(ForceClass%pIniForcetable)  .or.   &
              .not.associated(ForceClass%pCalForce)       .or.   &
              len_trim(ForceClass%ForceTable%PotType) .le. 0)   then
              write(*,fmt="(' MDPSCU Error: the force class is not registered on calling Init_Forcetable_Dev')")
              stop
           end if  

            if(present(RelaseTable) ) then
               call ForceClass%pIniForcetable(SimBox, CtrlParam, ForceClass%ForceTable, RelaseTable, MULTIBOX=0)
            else  
               call ForceClass%pIniForcetable(SimBox, CtrlParam, ForceClass%ForceTable, MULTIBOX=0)  
            end if   
 
          return
  end subroutine Init_Forcetable0
  !*********************************************************************************

  !*********************************************************************************
  subroutine Init_Forcetable1(SimBox, CtrlParam, ForceClass, RelaseTable)
   implicit none
    !--- dummy vaiables
         type(SimMDBox),dimension(:), intent(inout):: SimBox
         type(SimMDCtrl),             intent(in)   :: CtrlParam
         type(MDForceClassGPU),       intent(in)   :: ForceClass
         integer,                     optional     :: RelaseTable
    !--- local

         if(.not.associated(ForceClass%pIniForcetable)  .or.   &
            .not.associated(ForceClass%pCalForce)       .or.   &
             len_trim(ForceClass%ForceTable%PotType) .le. 0)   then
             write(*,fmt="(' MDPSCU Error: the force class is not registered on calling Init_Forcetable_Dev')")
             stop
           end if  

           if(present(RelaseTable)) then
              call ForceClass%pIniForcetable(SimBox(1), CtrlParam, ForceClass%ForceTable, RelaseTable, MULTIBOX=1)
           else  
              call ForceClass%pIniForcetable(SimBox(1), CtrlParam, ForceClass%ForceTable, MULTIBOX=1)  
           end if   
  
           return
  end subroutine Init_Forcetable1
  !*********************************************************************************
 
  !*********************************************************************************
  subroutine Init_Forcetable0_a(SimBox, CtrlParam)
  implicit none
    !--- dummy vaiables
         type(SimMDBox),  intent(inout) :: SimBox
         type(SimMDCtrl), intent(in)    :: CtrlParam
    !--- local 
              call Init_Forcetable0(SimBox, CtrlParam, gm_ForceClass)
           return
  end subroutine Init_Forcetable0_a
  !*********************************************************************************
  
  !*********************************************************************************
  subroutine Init_Forcetable1_a(SimBox, CtrlParam)
  implicit none
    !--- dummy vaiables
         type(SimMDBox),dimension(:), intent(inout):: SimBox
         type(SimMDCtrl),             intent(in)   :: CtrlParam
    !--- local 
              call Init_Forcetable1(SimBox, CtrlParam, gm_ForceClass)
           return
  end subroutine Init_Forcetable1_a
  !*********************************************************************************

  !*********************************************************************************
  logical function IfInit_ForceClass(ForceClass) result(YES)
  use  MD_FS_Force_Table_GPU,  only:IfInit_FS_Force_Table
  use  MD_EAM_Force_Table_GPU, only:IfInit_EAM_Force_Table
  implicit none
   type(MDForceClassGPU)::ForceClass

          select case(ForceClass%ForceTable%PotType(1:len_trim(ForceClass%ForceTable%PotType)))
                case("FS_TYPE")
                   YES = IfInit_FS_Force_Table()

                case("EAM_TYPE")
                   YES = IfInit_EAM_Force_Table()
                case("EXT_TYPE")   
                   YES = .false.
          end select

          return
  end function IfInit_ForceClass
  !*********************************************************************************

  !****************************************************************************************
  recursive subroutine Cumulate_ExtForce0(SimBox, CtrlParam, LIST)
  !***  DESCRIPTION:   to cumulate the external force
  !     INPUT:
  !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
  !                    the external routines should calculate the forces of atoms that have been
  !                    partitioned.
  !                    see also MD_NeighborsList_GPU.F90
   implicit none
      !--- dummy vaiables
      type(SimMDBox),       intent(inout)::SimBox
      type(SimMDCtrl),      intent(in)   ::CtrlParam
      type(MDExtForceList), intent(in)   ::LIST

            if(.not.associated(LIST%this)) then
               return
            else
              if( associated(LIST%this%pForce)) call LIST%this%pForce(SimBox, CtrlParam)
              if( associated(LIST%next) ) then
                  call Cumulate_ExtForce0(SimBox, CtrlParam, LIST%next)
              end if
            end if
            return
   end subroutine Cumulate_ExtForce0
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine Cumulate_ExtForce1(SimBox, CtrlParam, LIST)
  !***  DESCRIPTION:   to cumulate the external force
  !     INPUT:
  !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
  !                    the external routines should calculate the forces of atoms that have been
  !                    partitioned.
  !                    see also MD_NeighborsList_GPU.F90
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  dimension(:), intent(inout)::SimBox
      type(SimMDCtrl),               intent(in)   ::CtrlParam
      type(MDExtForceList),          intent(in)   ::LIST

            if(.not.associated(LIST%this)) then
               return
            else
              if( associated(LIST%this%pForceEx)) call LIST%this%pForceEx(SimBox, CtrlParam)
              if( associated(LIST%next) ) then
                  call Cumulate_ExtForce1(SimBox, CtrlParam, LIST%next)
              end if
            end if
            return
   end subroutine Cumulate_ExtForce1
   !****************************************************************************************

   !****************************************************************************************
   subroutine Calforce_ForceClass0(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to calculate forces of atoms
   !     INPUT:
   !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
   !                    the external routines should calculate the forces of atoms that have been
   !                    partitioned.
   !                    see also MD_NeighborsList_GPU.F90
   !
   implicit none
      !--- dummy vaiables
      type(SimMDBox),        intent(inout) :: SimBox
      type(SimMDCtrl),       intent(in)    :: CtrlParam
      type(MDForceClassGPU), intent(in)    :: ForceClass

            if(associated(ForceClass%pCalForce)) then
               call ForceClass%pCalForce(SimBox, CtrlParam)
            end if
            call Cumulate_ExtForce0(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine Calforce_ForceClass0
   !****************************************************************************************

   !****************************************************************************************
   subroutine Calforce_ForceClass1(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to calculate forces of atoms
   !     INPUT:
   !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
   !                    the external routines should calculate the forces of atoms that have been
   !                    partitioned.
   !                    see also MD_NeighborsList_GPU.F90
   !
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout) :: SimBox
      type(SimMDCtrl),              intent(in)    :: CtrlParam
      type(MDForceClassGPU),        intent(in)    :: ForceClass

            if(associated(ForceClass%pCalForce)) then
               call ForceClass%pCalForce(SimBox(1), CtrlParam)
            end if
            call Cumulate_ExtForce1(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine Calforce_ForceClass1
   !****************************************************************************************

   !****************************************************************************************
   subroutine CalPTensor_ForceClass0(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to calculate forces of atoms
   !     INPUT:
   !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
   !                    the external routines should calculate the forces of atoms that have been
   !                    partitioned.
   !                    see also MD_NeighborsList_GPU.F90
   !
   implicit none
      !--- dummy vaiables
      type(SimMDBox),        intent(inout) :: SimBox
      type(SimMDCtrl),       intent(in)    :: CtrlParam
      type(MDForceClassGPU), intent(in)    :: ForceClass

            if(associated(ForceClass%pCalPTensor)) then
               call ForceClass%pCalPTensor(SimBox, CtrlParam)
            end if
            call Cumulate_ExtForce0(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine CalPTensor_ForceClass0
   !****************************************************************************************

   !****************************************************************************************
   subroutine CalPTensor_ForceClass1(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to calculate forces of atoms
   !     INPUT:
   !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
   !                    the external routines should calculate the forces of atoms that have been
   !                    partitioned.
   !                    see also MD_NeighborsList_GPU.F90
   !
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout) :: SimBox
      type(SimMDCtrl),              intent(in)    :: CtrlParam
      type(MDForceClassGPU),        intent(in)    :: ForceClass

            if(associated(ForceClass%pCalPTensor)) then
               call ForceClass%pCalPTensor(SimBox(1), CtrlParam)
            end if
            call Cumulate_ExtForce1(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine CalPTensor_ForceClass1
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine Cumulate_ExtEpot0(SimBox, CtrlParam, LIST)
  !***  DESCRIPTION:   to cumulate the potentials on atoms
  !     INPUT:
  !
   implicit none
      !--- dummy vaiables
      type(SimMDBox),       intent(inout) :: SimBox
      type(SimMDCtrl),      intent(in)    :: CtrlParam
      type(MDExtForceList), intent(in)    :: LIST

            if(.not.associated(LIST%this)) then
               return
            else
              if( associated(LIST%this%pEpot)) call LIST%this%pEpot(SimBox, CtrlParam)
              if( associated(LIST%next) ) then
                  call Cumulate_ExtEpot0(SimBox, CtrlParam, LIST%next)
              end if
            end if
            return
   end subroutine Cumulate_ExtEpot0
   !****************************************************************************************

  !****************************************************************************************
  recursive subroutine Cumulate_ExtEpot1(SimBox, CtrlParam, LIST)
  !***  DESCRIPTION:   to cumulate the potentials on atoms
  !     INPUT:
  !
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout) :: SimBox
      type(SimMDCtrl),              intent(in)    :: CtrlParam
      type(MDExtForceList),         intent(in)    :: LIST

            if(.not.associated(LIST%this)) then
               return
            else
              if( associated(LIST%this%pEpotEx)) call LIST%this%pEpotEx(SimBox, CtrlParam)
              if( associated(LIST%next) ) then
                  call Cumulate_ExtEpot1(SimBox, CtrlParam, LIST%next)
              end if
            end if
            return
   end subroutine Cumulate_ExtEpot1
   !****************************************************************************************

   !****************************************************************************************
   subroutine CalEpot_ForceClass0(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to cumulate the potential
   !     INPUT:         SimBox, CtrlParam, ForceClass
   !
   !     NOTE:          the potentials on atoms BEFORE and AFTER partitioning are availbale.
   !
   !
   !                    see also: UpdateEpot_ForceClass
   implicit none
      !--- dummy vaiables
      type(SimMDBox),        intent(inout) :: SimBox
      type(SimMDCtrl),       intent(in)    :: CtrlParam
      type(MDForceClassGPU), intent(in)    :: ForceClass

            if(associated(ForceClass%pCalEpot)) then
               call ForceClass%pCalEpot(SimBox, CtrlParam)
            end if
            call Cumulate_ExtEpot0(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine CalEpot_ForceClass0
   !****************************************************************************************

   !****************************************************************************************
   subroutine CalEpot_ForceClass1(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to cumulate the potential
   !     INPUT:         SimBox, CtrlParam, ForceClass
   !
   !     NOTE:          the potentials on atoms BEFORE and AFTER partitioning are availbale.
   !
   !
   !                    see also: UpdateEpot_ForceClass
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout) :: SimBox
      type(SimMDCtrl),              intent(in)    :: CtrlParam
      type(MDForceClassGPU),        intent(in)    :: ForceClass

            if(associated(ForceClass%pCalEpot)) then
               call ForceClass%pCalEpot(SimBox(1), CtrlParam)
            end if
            call Cumulate_ExtEpot1(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine CalEpot_ForceClass1
   !****************************************************************************************

   !****************************************************************************************
   subroutine UpdateEpot_ForceClass0(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to update the potentials
   !     INPUT:         classname,   the name of the force class
   !
   !     NOTE:          the potentials are the potentials on atoms that HAVE BEEN partitioned on devices.
   !                    the external routines should calculate the potnetial of atoms that have been
   !                    partitioned.
   !                    see also: MD_NeighborsList_GPU.F90
   !
   !                     see also: CalEpot_ForceClass

   implicit none
      !--- dummy vaiables
      type(SimMDBox),       intent(inout) :: SimBox
      type(SimMDCtrl),      intent(in)    :: CtrlParam
      type(MDForceClassGPU),intent(in)    :: ForceClass

            if(associated(ForceClass%pCalEpot0)) then
               call ForceClass%pCalEpot0(SimBox, CtrlParam)
            end if
            call Cumulate_ExtEpot0(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine UpdateEpot_ForceClass0
  !****************************************************************************************

  !****************************************************************************************
   subroutine UpdateEpot_ForceClass1(SimBox, CtrlParam,  ForceClass)
   !***  DESCRIPTION:   to update the potentials
   !     INPUT:         classname,   the name of the force class
   !
   !     NOTE:          the potentials are the potentials on atoms that HAVE BEEN partitioned on devices.
   !                    the external routines should calculate the potnetial of atoms that have been
   !                    partitioned.
   !                    see also: MD_NeighborsList_GPU.F90
   !
   !                     see also: CalEpot_ForceClass

   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout) :: SimBox
      type(SimMDCtrl),              intent(in)    :: CtrlParam
      type(MDForceClassGPU),        intent(in)    :: ForceClass

            if(associated(ForceClass%pCalEpot0)) then
               call ForceClass%pCalEpot0(SimBox(1), CtrlParam)
            end if
            call Cumulate_ExtEpot1(SimBox, CtrlParam, ForceClass%ExtForceList)

            return
    end subroutine UpdateEpot_ForceClass1
  !****************************************************************************************

  !*********************************************************************************
  recursive subroutine Add_ExtForce(List, Tag, pForceArray, pEpotArray, pForce, pEpot)
  !***  DESCRIPTION:   to add a externale force to a external force-list
  !     INPUT:         LIST,   the externale force list class
  !                    Tag,    the module name of the routines
  !                    pForce, the external routine to calcualte external force
  !                    pEpot,  the external routine to calcualte external potential
  !
  !
  implicit none
      type(MDExtForceList)::List
      character(*),             intent(in)::Tag
      procedure(CALFORCE_EXT),  pointer::pForceArray
      procedure(CALEPOT_EXT),   pointer::pEpotArray
      procedure(CALFORCE),      pointer::pForce
      procedure(CALEPOT),       pointer::pEpot

     !$$--- local varibales

            if(.not.associated(List%this)) then
               allocate(List%this)
               List%Tag            =  Tag(1:min(len_trim(Tag),len_trim(LIST%Tag)))
               LIST%this%pForce    => pForce
               LIST%this%pEpot     => pEpot
               List%this%pForceEx  => pForceArray
               List%this%pEpotEx   => pEpotArray
            else
              !--- replace the preexisting procedures
              if(List%Tag .eq. Tag) then
                 LIST%this%pForce   => pForce
                 LIST%this%pEpot    => pEpot
                 List%this%pForceEx => pForceArray
                 List%this%pEpotEx  => pEpotArray
              else
                 if(.not. associated(List%next) ) then
                     allocate(List%next)
                     List%next%this => null() 
                 end if
                 call Add_ExtForce(List%next, Tag, pForceArray, pEpotArray, pForce, pEpot)
              end if
            end if
            return
   end subroutine Add_ExtForce
  !*********************************************************************************

  !*********************************************************************************
  recursive subroutine Del_ExtForce(List, Tag)
  !***  DESCRIPTION:   to delete a externale force to a external force-list
  !     INPUT:         LIST,   the externale force list class
  !                    Tag,    the module name of the routines
  !
  implicit none
      type(MDExtForceList)     ::List
      character(*),  intent(in)::Tag
     !$$--- local varibales
      type(MDExtForce),     pointer::THIS
      type(MDExtForceList), pointer::NEXT

              if(LIST%Tag .eq. Tag) then
                 THIS  => List%this
                 NEXT  => List%next
                 if(associated(NEXT)) then
                    List%Tag  =  NEXT%Tag
                    List%this => NEXT%This
                    List%next => NEXT%Next
                 end if
                 if(associated(THIS)) deallocate(THIS)
                 if(associated(NEXT)) deallocate(NEXT)
              else
                 if(associated(List%next)) then
                    call Del_ExtForce(List%next, Tag)
                 end if
              end if
            return
   end subroutine Del_ExtForce
  !****************************************************************************************

  !*********************************************************************************
  recursive subroutine Clear_ExtForce(List)
  !***  DESCRIPTION:   to clear a external force class list
  !     INPUT:         List,   the externale force list class
  !     OUTPUT:        List,   the externale force list class
  !
  !
  implicit none
      type(MDExtForceList)::List

     !$$--- local varibales

            if(associated(List%Next)) then
               call Clear_ExtForce(List%Next)
               deallocate(List%Next)
            end if    

            if(associated(List%this)) then
               deallocate(List%this)
            end if   
            List%Tag = "" 
            return
  end subroutine Clear_ExtForce
  !*********************************************************************************

  !*********************************************************************************
  subroutine AddExtForce_ForceClass0(ForceClass, Tag, ForceArray, EpotArray, Force, Epot)
  !***  DESCRIPTION:   to add a externale force to a external force-list
  !     INPUT:         ForceClass,  the MDForceClassGPU force class
  !                    Tag,         the module name of the routines
  !                    CALFORCE,    the external routine to calcualte external force
  !                    CALPOTE,     the external routine to calcualte external potential
  !
  implicit none
      !--- dummy variables
      type(MDForceClassGPU)  :: ForceClass
      character(*),intent(in)::Tag
      optional::ForceArray, EpotArray, Force, Epot
      external::ForceArray, EpotArray, Force, Epot

     !--- interface to the external routine -------------------
      interface
        subroutine Force(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),   intent(inout)::SimBox
          type(SimMDCtrl),  intent(in)   ::CtrlParam
        end subroutine Force
      end interface

      interface
        subroutine ForceArray(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),dimension(:), intent(inout)::SimBox
          type(SimMDCtrl),             intent(in)   ::CtrlParam
        end subroutine ForceArray
      end interface

      !--- interface to potential calculation
      interface
        subroutine Epot(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),  intent(inout)::SimBox
          type(SimMDCtrl), intent(in)   ::CtrlParam
        end subroutine Epot
      end interface

      interface
        subroutine EpotArray(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),dimension(:), intent(inout)::SimBox
          type(SimMDCtrl),             intent(in)   ::CtrlParam
        end subroutine EpotArray
      end interface

      !--- local varibales
      procedure(CALFORCE_EXT), pointer::pF
      procedure(CALEPOT_EXT),  pointer::pE
      procedure(CALFORCE),     pointer::pF0
      procedure(CALEPOT),      pointer::pE0

            if(present(Force)) then
               pF0 => Force
            else
               pF0 =>null()
            end if

            if(present(Epot)) then
               pE0 => Epot
            else
               pE0 =>null()
            end if

            if(present(ForceArray)) then
               pF => ForceArray
            else
               pF =>null()
            end if

            if(present(EpotArray)) then
               pE => EpotArray
            else
               pE =>null()
            end if
            
            call Add_ExtForce(LIST=ForceClass%ExtForceList, &
                              Tag = Tag,                    &
                              pForceArray  = pF,            &
                              pEpotArray   = pE,            &
                              pForce       = pF0,           &
                              pEpot        = pE0)
            return
   end subroutine AddExtForce_ForceClass0
  !****************************************************************************************

  !*********************************************************************************
  subroutine AddExtForce_ForceClass1(Tag, ForceArray, EpotArray, Force, Epot)
  !***  DESCRIPTION:   to add a externale force to a external force-list
  !     INPUT:         ForceClass,  the MDForceClassGPU force class
  !                    Tag,         the module name of the routines
  !                    CALFORCE,    the external routine to calcualte external force
  !                    CALPOTE,     the external routine to calcualte external potential
  !
  implicit none
      !--- dummy variables
      type(MDForceClassGPU)  :: ForceClass
      character(*),intent(in)::Tag
      optional::ForceArray, EpotArray, Force, Epot
      external::ForceArray, EpotArray, Force, Epot

     !--- interface to the external routine -------------------
      interface
        subroutine Force(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),   intent(inout)::SimBox
          type(SimMDCtrl),  intent(in)   ::CtrlParam
        end subroutine Force
      end interface

      interface
        subroutine ForceArray(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),dimension(:), intent(inout)::SimBox
          type(SimMDCtrl),             intent(in)   ::CtrlParam
        end subroutine ForceArray
      end interface

      !--- interface to potential calculation
      interface
        subroutine Epot(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),  intent(inout)::SimBox
          type(SimMDCtrl), intent(in)   ::CtrlParam
        end subroutine Epot
      end interface

      interface
        subroutine EpotArray(SimBox, CtrlParam)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          !--- dummy vaiables
          type(SimMDBox),dimension(:), intent(inout)::SimBox
          type(SimMDCtrl),             intent(in)   ::CtrlParam
        end subroutine EpotArray
      end interface

      !--- local varibales
      procedure(CALFORCE_EXT), pointer::pF
      procedure(CALEPOT_EXT),  pointer::pE
      procedure(CALFORCE),     pointer::pF0
      procedure(CALEPOT),      pointer::pE0

            if(present(Force)) then
               pF0 => Force
            else
               pF0 =>null()
            end if

            if(present(Epot)) then
               pE0 => Epot
            else
               pE0 =>null()
            end if

            if(present(ForceArray)) then
               pF => ForceArray
            else
               pF =>null()
            end if

            if(present(EpotArray)) then
               pE => EpotArray
            else
               pE =>null()
            end if
            
            call Add_ExtForce(LIST=gm_ForceClass%ExtForceList, &
                              Tag = Tag,                       &
                              pForceArray  = pF,               &
                              pEpotArray   = pE,               &
                              pForce       = pF0,              &
                              pEpot        = pE0)
            return
   end subroutine AddExtForce_ForceClass1
  !****************************************************************************************   

  !*********************************************************************************
  subroutine ClearExtForce_ForceClass0(ForceClass)
  !***  DESCRIPTION:   to clearthe externale forces in  ForceClass
  !
  implicit none
      !--- dummy variables
      type(MDForceClassGPU)  :: ForceClass
            
            call Clear_ExtForce(ForceClass%ExtForceList)
            return
   end subroutine ClearExtForce_ForceClass0
  !****************************************************************************************

  !*********************************************************************************
  subroutine ClearExtForce_ForceClass1()
  !***  DESCRIPTION:   to clear the external force-list in gm_ForceClass
  !
  implicit none
      !--- dummy variables
      type(MDForceClassGPU)  :: ForceClass
            
            call Clear_ExtForce(gm_ForceClass%ExtForceList)
            return
   end subroutine ClearExtForce_ForceClass1
  !****************************************************************************************  

  end module MD_ForceClass_Register_GPU
  !****************************************************************************************
