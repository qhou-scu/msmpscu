  module MD_ForceClass_Register
  !***  DESCRIPTION:
  !     This module provides interfaces to force calculations using different TYPE of
  !     potentials. In this module, no concrete potential is implmented, but the type
  !     potential is registered
  !
  !    Written by HOU Qing, May, 2014
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_TYPEDEF_ForceTable

   implicit none

  !--- interface to the initalization of force table
   private::INIT_FORCETABLE
   abstract interface
    subroutine INIT_FORCETABLE(SimBox, CtrlParam, FTable, RelaseTable, MULTIBOX)
    !***  PURPOSE:   to intialize the parameters to be used in force calculation
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
      type(SimMDBox),     intent(inout)::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(in)   ::FTable
      integer,            optional     ::RelaseTable, MULTIBOX

    end subroutine INIT_FORCETABLE
  end interface

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

  !--- interface to force calculationm
   private::CALFORCE
   abstract interface
    subroutine CALFORCE(SimBox, List)
    !***  PURPOSE:   to begine calculate the force, Newton's third law taken into account
    !
    !     INPUT:     SimBox,    the simulation box
    !                List,      the neighbore list
    !     OUTPUT     SimBox,    the simulation box with force updated
    !
      use MD_TYPEDEF_SimMDBox
      use MD_NeighborsList
      implicit none
      !--- dummy vaiables
      type(SimMDBox),             intent(inout)::SimBox
      type(NEIGHBOR_LIST),target, intent(in)   ::List
    end subroutine CALFORCE
  end interface

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

 !--- END INTERFACE --------------------------------

 !--- type definition of external force list
      type,private::MDExtForce
           procedure(CALEPOT),              pointer, nopass::pEpot  =>null()
           procedure(CALFORCE),             pointer, nopass::pForce  =>null()
      end type MDExtForce

      type, private::MDExtForceList
           character(len=16)                       :: Tag  = ""
           type(MDExtForce),                pointer:: this =>null()
           type(MDExtForceList),            pointer:: next =>null()
      end type MDExtForceList

 !--- type definition of FORCE CLASS
      type MDForceClassCPU
           procedure(INI_FORCETABLE),        pointer, nopass::pIniForcetable=>null()
           procedure(CLEAR_FORCETABLE),      pointer, nopass::pClrForcetable=>null()
           procedure(CALFORCE),              pointer, nopass::pCalForce=>null()
           type(MDForceTable)                               ::ForceTable
           type(MDExtForceList)                             ::ExtForceList
      end type MDForceClassCPU

      private::Add_ExtForce,  Cumulate_ExtForce

 !--- public accessible members
      type(MDForceClassCPU)::gm_ForceClass
  contains

  !****************************************************************************************
  subroutine Register_ForceClass(classname, ForceClass)
  !***  DESCRIPTION:   to regist a force class
  !     INPUT:         classname,   the name of the force class
  !
  !
  !$$-------FS moudle-----------------------
  use  MD_FS_Force_Table, only:FS_INIT=>INITIALIZE_FS_Force_Table, &
                                    FS_CALF=>CALFORCE_FS_Force_Table,   &
                                    FS_CLR=>RELEASE_FS_Force_Table

  !use  MD_FS_Force_Table, only:EAM_INIT=>INITIALIZE_EAM_Force_Table_DEV,            &
  !                                                EAM_CALF=>CALFORCE_EAM_Force_Table2A_DEV,   &
  !                                                EAM_CALE=>CALEPOT_EAM_Force_Table2A_DEV,    &
  !                                                EAM_CALP=>CALPTENSOR_EAM_Force_Table2A_DEV, &
  !                                                EAM_CLR=>Clear_EAM_Force_Table_DEV,         &
  !                                                EAM_CALD=>CALDEN_EAM_Force_Table2A_DEV,     &
  !                                                EAM_CALAVS=>Cal_EAM_AtomicStressTensor_DEV


  !$$------------------------------
  implicit none
  character*(*)::classname
  type(MDForceClassCPU)::ForceClass

  !$$--- local varibales

       select case(classname(1:len_trim(classname)))
              case("FS_TYPE")
                   ForceClass%ForceTable%PotType = "FS_TYPE"
                   ForceClass%pIniForcetable=>FS_INIT
                   ForceClass%pCalForce=>FS_CALF
                   ForceClass%pClrForcetable=>FS_CLR

              case("EAM_TYPE")
                  ! ForceClass%ForceTable%PotType = "EAM_TYPE"
                  ! ForceClass%pIniForcetable=>EAM_INIT
                  ! ForceClass%pCalForce=>EAM_CALF
                  ! ForceClass%pClrForcetable=>EAM_CLR
                  write(*,fmt="(' MDPSCU Error: EAM type not supportted in CPU version', A16)")
              case("EXT_TYPE")
                  ForceClass%ForceTable%PotType = "EXT_TYPE"
                  ForceClass%pIniForcetable     =>NULL()
                  ForceClass%pCalForce          =>NULL()
                  ForceClass%pClrForcetable     =>NULL()
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
  !****************************************************************************************

  !*********************************************************************************
  logical function IfInit_ForceClass(ForceClass) result(YES)
  use  MD_FS_Force_Table, only:IfInit_FS_Force_Table
  !use  MD_EAM_Force_Table_2012_GPU, only:IfInit_EAM_Force_Table
  implicit none
   type(MDForceClassCPU)::ForceClass

          select case(ForceClass%ForceTable%PotType(1:len_trim(ForceClass%ForceTable%PotType)))
                case("FS_TYPE")
                   YES = IfInit_FS_Force_Table()

                case("EAM_TYPE")
                    YES = .false.
                   !YES = IfInit_EAM_Force_Table()

                case("EXT_TYPE")
                    YES = .false.
          end select

          return
  end function IfInit_ForceClass
  !*********************************************************************************

   !****************************************************************************************
   subroutine Calforce_ForceClass(SimBox, NLIST,  ForceClass)
   !***  DESCRIPTION:   to calculate forces of atoms
   !     INPUT:
   !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
   !                    the external routines should calculate the forces of atoms that have been
   !                    partitioned.
   !                    see also MD_NeighborsList_GPU.F90
   !
   use MD_NeighborsList, only:NEIGHBOR_LIST
   implicit none
      !--- dummy vaiables
      type(SimMDBox),        intent(inout):: SimBox
      type(NEIGHBOR_LIST),   intent(in)   :: NLIST
      type(MDForceClassCPU), intent(in)   :: ForceClass

            if(associated(ForceClass%pCalForce)) then
               call ForceClass%pCalForce(SimBox, NLIST)
            end if
            call Cumulate_ExtForce(SimBox, NLIST, ForceClass%ExtForceList)

            return
    end subroutine Calforce_ForceClass
   !****************************************************************************************

  !*********************************************************************************
  recursive subroutine Add_ExtForce(LIST, Tag, pForce, pEpot)
  !***  DESCRIPTION:   to add a externale force to a external force-list
  !     INPUT:         LIST,     the externale force list class
  !                    Tag,      the module name of the routines
  !                    pForce,   the external routine to calcualte external force
  !                    pEpot,    the external routine to calcualte external potential
  !
  !
  implicit none

      type(MDExtForceList)::LIST
      character(*),                 intent(in)::Tag
      procedure(CALFORCE), pointer, intent(in)::pForce
      procedure(CALEPOT),  pointer, intent(in)::pEpot

     !$$--- local varibales

            if(.not.associated(LIST%this)) then
               allocate(LIST%this)
               LIST%Tag          =  Tag(1:min(len_trim(Tag),len_trim(LIST%Tag)))
               LIST%this%pForce  => pForce
               LIST%this%pEpot   => pEpot
            else
              !--- replace the preexisting procedures
              if(LIST%Tag .eq. Tag) then
                 LIST%this%pForce  => pForce
                 LIST%this%pEpot   => pEpot
              else
                 if(.not. associated(LIST%next) ) then
                     allocate(LIST%next)
                     LIST%Tag         =  ""
                     LIST%next%this   => null()
                 end if
                 call Add_ExtForce(LIST%next, Tag, pForce, pEpot)
              end if
            end if
            return
   end subroutine Add_ExtForce
  !****************************************************************************************

  !*********************************************************************************
  subroutine AddExtForce_ForceClass(ForceClass, Tag, pForce, pEpot)
  !***  DESCRIPTION:   to add a externale force to a external force-list
  !     INPUT:         ForceClass,   the MDForceClassGPU force class
  !                    Tag,          the module name of the routines
  !                    pForce,       the external routine to calcualte external force
  !                    pEpot,        the external routine to calcualte external potential
  !
  implicit none
      !--- dummy variables
      type(MDForceClassCPU)                             ::ForceClass
      character(*),                           intent(in)::Tag
      procedure(CALFORCE), pointer, optional, intent(in)::pForce
      procedure(CALEPOT),  pointer, optional, intent(in)::pEpot

      !--- local varibales
      procedure(CALFORCE), pointer::pF
      procedure(CALEPOT),  pointer::pE

            if(present(pForce)) then
               pF => pForce
            else
               pF =>null()
            end if

            if(present(pEpot)) then
               pE => pEpot
            else
               pE =>null()
            end if

            call Add_ExtForce(ForceClass%ExtForceList, Tag, pF, pE)
            return
   end subroutine AddExtForce_ForceClass
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine Cumulate_ExtForce(SimBox, NLIST, FLIST)
  !***  DESCRIPTION:   to cumulate the external force
  !     INPUT:         SimBox, the simulation box
  !                    NLIST,  the neighbor list
  !                    FLIST,  the external force list
  !     NOTE:          the forces are the forces on atoms that have been partitioned on devices.
  !                    the external routines should calculate the forces of atoms that have been
  !                    partitioned.
  !                    see also MD_NeighborsList_GPU.F90
   use MD_NeighborsList, only:NEIGHBOR_LIST
   implicit none
      !--- dummy vaiables
      type(SimMDBox),       intent(inout):: SimBox
      type(NEIGHBOR_LIST),  intent(in)   :: NLIST
      type(MDExtForceList), intent(in)   :: FLIST

            if(.not.associated(FLIST%this)) then
               return
            else
              if( associated(FLIST%this%pForce)) call FLIST%this%pForce(SimBox, NLIST)
              if( associated(FLIST%next) ) then
                  call Cumulate_ExtForce(SimBox, NLIST, FLIST%next)
              end if
            end if
            return
   end subroutine Cumulate_ExtForce
   !****************************************************************************************

  end module MD_ForceClass_Register
  !****************************************************************************************
