  module MD_MethodClass_Register_GPU
  !***  DESCRIPTION:
  !     This module provides interfaces to force calculations using different TYPE of
  !     potentials. In this module, no concrete potential is implmented, but the type
  !     potential is registered
  !
  !    Written by HOU Qing, May, 2014
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl

   implicit none

  !--- interface to the initalization of force table
   private::FOR_ONE_TEST
   abstract interface
        subroutine FOR_ONE_TEST(SimBox0, CtrlParam, SimBox, Recordlist, J, processid, INICONFIGPROC)
         !***  PORPOSE: to process one sample
         !     INPUT:  SimBox0,            the element sample
         !             CtrlParam,          the control parameters
         !             Recordlist          the external recoding routines
         !             J,                  the test ID
         !             processid,          id of the process, to be used if MPI used
         !             INICONFIGPROC,      the subroutine privided by user for initializing the initial configure
         !
         use MD_TYPEDEF_SimMDBox
         use MD_TYPEDEF_SimMDCtrl
         use MD_TYPEDEF_RecordList

         implicit none
         !--- dummy variables
         type(SimMDBox),                            intent(in)::SimBox0
         type(SimMDCtrl),target                               ::CtrlParam
         type(SimMDBox), dimension(:), allocatable            ::SimBox
         type(RecordProcedureList)                            ::Recordlist
         integer,                                   intent(in)::J,processid

         !--- interface for the external process--------------
         OPTIONAL::                                     INICONFIGPROC
         external::                                     INICONFIGPROC

         interface
             subroutine INICONFIGPROC(SimBox, CtrlParam, RESTART)
             use MD_TYPEDEF_SimMDBox
             use MD_TYPEDEF_SimMDCtrl
             implicit none
             type(SimMDBox), dimension(:)::SimBox
             type(SimMDCtrl)             ::CtrlParam
             integer                     ::RESTART
             end subroutine INICONFIGPROC
         end interface
      end subroutine FOR_ONE_TEST
   end interface
  !------END INTERFACE ---

  !--- type definition of APP CLASS
      type MDMethodClassGPU
           character(len=8)::name
           type(SimMDBox),          pointer         ::pSimBox=>null()
           type(SimMDCtrl),         pointer         ::pCtrlParam=>null()
           procedure(FOR_ONE_TEST), pointer, nopass ::pForOnetest=>null()
      end type MDMethodClassGPU

  contains

  !****************************************************************************************
  subroutine Register_MethodClass(classname, MethodClass, SimBox, CtrlParam)
  !***  DESCRIPTION:   to regist a app class
  !     INPUT:         classname,   the name of the simulation method class
  !
  !
  use MD_LocalTempMethod_GPU,  only:Initialize_LocalTempMethod_DEV
  use MD_DiffScheme_GPU,      only:Initialize_Diff_Scheme_DEV

  !--- the available methods
  use MD_Method_GenericMD_GPU, only:For_One_Test_GMD      =>For_One_Test
  use MD_Method_ParRep_GPU,    only:For_One_Test_ParRep   =>For_One_Test
  use MD_Method_TAD_GPU,       only:For_One_Test_TAD      =>For_One_Test
  use MD_Method_NEB_GPU,       only:For_One_Test_NEB      =>For_One_Test
  use MD_Method_ART_GPU,       only:For_One_Test_ART      =>For_One_Test
  use MD_Method_BST_GPU,       only:For_One_Test_BST      =>For_One_Test

  implicit none
       character*(*)           :: classname
       type(MDMethodClassGPU)  :: MethodClass
       type(SimMDBox),target   :: SimBox
       type(SimMDCtrl),target  :: CtrlParam

  !$$--- local varibales
        MethodClass%pSimBox    => SimBox
        MethodClass%pCtrlParam => CtrlParam
        select case(classname(1:len_trim(classname)))
               case("GMD")
                   MethodClass%name          =  "GMD"
                   MethodClass%pForOnetest   =>  For_One_Test_GMD

               case("PARREP")
                    MethodClass%name          =  "PARREP"
                    MethodClass%pForOnetest   =>  For_One_Test_ParRep

               case("TAD")
                    MethodClass%name          =  "TAD"
                    MethodClass%pForOnetest   =>  For_One_Test_TAD

               case("NEB")
                    MethodClass%name          =  "NEB"
                    MethodClass%pForOnetest   =>  For_One_Test_NEB

               case("ART")
                    MethodClass%name          =  "ART"
                    MethodClass%pForOnetest   =>  For_One_Test_ART

               case("BST")
                    MethodClass%name          =  "BST"
                    MethodClass%pForOnetest   =>  For_One_Test_BST

               case default
                   if(len_trim(classname).le.0) then
                      write(*,fmt="(' MDPSCU Error: empty method', A16)")
                    else
                      write(*,fmt="(' MDPSCU Error: unsupported method ', A16)") classname
                    end if
                   stop
        end select

        select case(classname(1:len_trim(classname)))
               case("GMD","PARREP","TAD")
                   !*** to initialize the electron phonon coupling module
                    call Initialize_LocalTempMethod_DEV(MethodClass%pSimBox, MethodClass%pCtrlParam)

                   !*** to initialize the differential scheme
                    call Initialize_Diff_Scheme_DEV(MethodClass%pSimBox, MethodClass%pCtrlParam)
               case("NEB")
                   !*** to initialize the differential scheme
                   !  call Initialize_Diff_Scheme_DEV(MethodClass%pSimBox, MethodClass%pCtrlParam)

               case("ART")
                   !*** to initialize the differential scheme
                    call Initialize_Diff_Scheme_DEV(MethodClass%pSimBox, MethodClass%pCtrlParam)

               case("BST")
                    !*** to initialize the differential scheme
                    !  call Initialize_Diff_Scheme_DEV(MethodClass%pSimBox, MethodClass%pCtrlParam)
         end select

       return
   end subroutine Register_MethodClass
  !****************************************************************************************


  !****************************************************************************************
  subroutine Clear_MethodClass(MethodClass)
  !***  DESCRIPTION:   to clear the momeory allocated when register the method
  !     INPUT:         MethodClass,   the method class
  !
  !
  use MD_DiffScheme_GPU, only:Clear_Diff_Scheme_DEV
  !--- the available methods
  implicit none
       type(MDMethodClassGPU)::MethodClass

  !$$--- local varibales
        MethodClass%pSimBox       => null()
        MethodClass%pCtrlParam    => null()
        MethodClass%pForOnetest   => null()
        select case(MethodClass%name(1:len_trim(MethodClass%name)))
               case("GMD","PARREP", "TAD")
                   !*** to initialize the differential scheme
                    call Clear_Diff_Scheme_DEV()
               case("NEB")
                   !*** to initialize the differential scheme
                    call Clear_Diff_Scheme_DEV()
               case("ART")
                   !*** to initialize the differential scheme
                    call Clear_Diff_Scheme_DEV()
       end select

       return
   end subroutine Clear_MethodClass
  !****************************************************************************************


  end module MD_MethodClass_Register_GPU
  !****************************************************************************************
