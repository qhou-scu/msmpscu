  module MD_SpringBoost_GPU
  !**** DESCRIPTION: This module provides routines to spring-boost force
  !                  on potential surface along a given direction. This module could 
  !                  be be used in steepest descent or CG of search minima on energy surface         
  !
  !                  ______________________________________________________
  !**** HISTORY:     2019-01-16(HOU Qing),
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_Globle_Variables_GPU
  implicit none
  !**********************
        !--- the data type define the sping-boost-force
        private::SPBPotent
        type::SPBPotent
            type(DevVec_DF), dimension(m_MXDEVICE)::K0 
            type(DevVec_DF), dimension(m_MXDEVICE)::D0 
            type(DevVec_I),  dimension(m_MXDEVICE)::PId     !--- id of patiticles to be boosted, the ID is in riginal box (before partition)
            type(DevMat_DF), dimension(m_MXDEVICE)::Pos     !--- position of particles to be bbosted
        end type SPBPotent
        
        #define   LIISTNAME SPBPotent_List
        #define   DATATYPE  type(SPBPotent) 
        #define   DEALLOCTE_DATA  Clear_SPBPotent
        private   SPBPotent_List
        #include "MSM_MultiGPU_DatList_Header.F90"

        !--- the workspace 
        type::SPBWKSpace
            integer::NPOT = 0
            type(MDDEVWorkSpace), pointer::MDWKSpace=>null()
            type(SPBPotent),      pointer::LastPot =>null()
            type(SPBPotent_List)         ::BoostPotList
        end type SPBWKSpace
        type(SPBWKSpace), private::m_SPBWKSpace            !--- default workspace
        !--- the interafce
        private::Clear_SPBPotent
        
        !---------------------------------
        private:: Add_DefaultSPBPotent, &
                  Add_SPBPotent0
        public::  Add_BoostForce
        interface Add_BoostForce
           module procedure Add_DefaultSPBPotent
           module procedure Add_SPBPotent0
        end interface Add_BoostForce 
        !---------------------------------
        private:: Attach0,  &
                  Attach1
        public::  Reset_BoostForce
        interface Reset_BoostForce
           module procedure Attach0
           module procedure Attach1
        end interface Reset_BoostForce
        !---------------------------------
        private:: Clear_DefaultSPBWKSpace,  &
                  Clear_SPBWKSpace
        public::  Clear_BoostWKSpace
        interface Clear_BoostWKSpace
           module procedure Clear_DefaultSPBWKSpace
           module procedure Clear_SPBWKSpace
        end interface Clear_BoostWKSpace

        

  contains
  !*************************************************************************************
  #include "MSM_MultiGPU_DatList_Body.F90"
  !****************************************************************************************
  subroutine Clear_SPBPotent(theSPBPotent)
  !***  PORPOSE:   to clear the memory allocated in  SPBPotent_List
  !
  !     INPUT:     theSPBPotent
  !     OUTPUT:    theSPBPotent
  !
  implicit none
  !----   DUMMY Variables
          type(SPBPotent) :: theSPBPotent
  !----  Local variables
            call DevDeallocate(theSPBPotent%Pos)

  end subroutine Clear_SPBPotent
  !****************************************************************************************
  
  !****************************************************************************************
  subroutine Clear_SPBWKSpace(theSPBWKSpace)
   !***  PORPOSE:   to clear the memory allocated in SPBWKSpace
   !
   !     INPUT:     SPBWKSpace
   !     OUTPUT:    SPBWKSpace
   !
   implicit none
   !----  DUMMY Variables
          type(SPBWKSpace) :: theSPBWKSpace
   !----  Local variables
             
           call ClearDataFromList(theSPBWKSpace%BoostPotList)
           theSPBWKSpace%LastPot   => null()
           theSPBWKSpace%MDWKSpace => null()
           theSPBWKSpace%NPOT    = 0
   end subroutine Clear_SPBWKSpace
  !****************************************************************************************

  !****************************************************************************************
  subroutine Clear_DefaultSPBWKSpace()
  !***  PORPOSE:   to clear the memory allocated in the default workspace
  !
  implicit none
  !----   DUMMY Variables
  !----  Local variables
            
          call Clear_SPBWKSpace(m_SPBWKSpace)
  end subroutine Clear_DefaultSPBWKSpace
  !****************************************************************************************

  !****************************************************************************************
  subroutine Attach0(MDWKS, theSPBWKSpace)
   !***  PORPOSE:   to connect the SPB workspace to MD workspace
   !
   !     INPUT:     
   !     OUTPUT:  
   !
   implicit none
   !----   DUMMY Variables
           type(MDDEVWorkSpace), target::MDWKS
           type(SPBWKSpace)           ::theSPBWKSpace
   !----  Local variables
             call Clear_SPBWKSpace(theSPBWKSpace)
             theSPBWKSpace%MDWKSpace => MDWKS
           return
   end subroutine Attach0
   !****************************************************************************************  

   !****************************************************************************************
   subroutine Attach1()
   implicit none
   !----   DUMMY Variables
   !----  Local variables
              call Attach0(dm_WorkSpace, m_SPBWKSpace)
           return
   end subroutine Attach1
  !****************************************************************************************   

  !****************************************************************************************
  subroutine Add_SPBPotent0(NB, K0, D0, Pid0, Pos0, theSPBWKSpace)
  !***  PORPOSE:   to add a boost state to the boost list
  !
  !     INPUT:     K0,        the "spring" constant
  !                D0,        the length of the spring force
  !                Pos0,      the configuration of minimum configuration of a minia on  potenmtial  surface
  !
  !     OUTPUT:    theList:   the boost state list
  !
  implicit none
  !----   DUMMY Variables
          integer                                  ::NB
          real(KINDDF),  intent(in)                ::K0, D0
          integer,       intent(in),dimension(:)   ::PId0
          real(KINDDF),  intent(in),dimension(:,:) ::Pos0
          type(SPBWKSpace)                         ::theSPBWKSpace
  !----  Local variables
          type(SPBPotent), pointer::pNewState
          integer::NPRT
            
            allocate(pNewState)
            NPRT = size(PId0) 
            call DevAllocate(pNewState%POS, (/NPRT,3/), "BOOST_STATE") 
            call DevAllocate(pNewState%K0,   NB,   "BOOST_STRENGTH")
            call DevAllocate(pNewState%D0,   NB,   "BOOST_LENGTH")
            call AddDataToList(theSPBWKSpace%BoostPotList, pNewState)
            theSPBWKSpace%LastPot => pNewState
            theSPBWKSpace%NPOT    =  theSPBWKSpace%NPOT+1

            call DevCopyIn(PId0, pNewState%PId,  NPRT)
            call DevCopyIn(POS0, pNewState%POS,3*NPRT)
            call DevSet(pNewState%K0, K0)
            call DevSet(pNewState%D0, D0)
          return
  end subroutine Add_SPBPotent0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Add_DefaultSPBPotent(NB, K0, D0, PId0, Pos0)
  !***  PORPOSE:   to add a boost state to the boost list
  !
  !     INPUT:     K0,        the "spring" constant
  !                D0,        the length of the spring force
  !                Pos0,      the configuration of minimum configuration of a minia on  potenmtial  surface
  !
  !     OUTPUT:    theList:   the boost state list
  !
  implicit none
  !----   DUMMY Variables
          integer                                  ::NB
          real(KINDDF),  intent(in)                ::K0, D0
          integer,       intent(in),dimension(:)   ::PId0
          real(KINDDF),  intent(in),dimension(:,:) ::Pos0
  !----  Local variables
            
            call Add_SPBPotent0(NB, K0, D0, PId0, Pos0, m_SPBWKSpace)
          return
  end subroutine Add_DefaultSPBPotent
  !****************************************************************************************



  end module MD_SpringBoost_GPU
