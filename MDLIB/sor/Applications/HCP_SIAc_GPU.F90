 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to simulate the diffusion of SIA clusters in HCP crystal.
 !                  A HCP template box need to be load as a reference lattice to identify the localtion of SIA (LSIA).
 !                  If a LSIA is identified, a new atom is added to one of the neighbore lattice of the LSIA.
 !                  The module RefVoronoiVacancy_GPU to be used to identify LSIA.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       RefVoronoiVA_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  by HOU Qing, Oct., 2015
 !
  module HCP_SIAc_GPU_
  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_NeighborsList

  implicit none
         integer, private, parameter::m_MXNEARSITE  = 14
         integer, private::m_NEARSITE    = 14
         integer, parameter, private::mp_MXI  = 10    ! the permitted number of interstital in a box
         integer, parameter, private::mp_MXAI = 10    ! the permitted number of atoms occupy one site

         type(NEIGHBOR_LIST)::m_SiteNeighbor

  contains

  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate memories, and create the
   !               neighbor-list of the reference lattice
   !
   use RefVoronoiVacancy_GPU, only:Initialize_ReferenceVacancy_DEV, hm_RefSimBox
   use MD_NeighborsList_GPU, only:Reorder_NeighBoreList_Nearest_DEV, COPYOUT_NeighboreList_DEV
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables

         !$$--- loading control parameters
         call Initialize_ReferenceVacancy_DEV(SimBox, CtrlParam)

         !$$--- We keep the neighbor list of the reference lattice
         call Reorder_NeighBoreList_Nearest_DEV(hm_RefSimBox, CtrlParam, m_NEARSITE, m_SiteNeighbor)
         return
   end subroutine MyInitialize
 !**********************************************************************************

  !**********************************************************************************
  subroutine Insert_NewAtom(SimBox, CtrlParam)
  !***  DESCRIPTION: to insert a new atom to the box
  !
  !    INPUT:  SimBox, the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use MD_SimboxArray_GPU
  use MD_NeighborsList_GPU, only:Initialize_NeighboreList_DEV, Cal_NeighBoreList_DEV
  use RefVoronoiVacancy_GPU, only:Cal_Occupied_Voronoi_DEV, hm_RefSimBox
  use RAND32_MODULE

  implicit none
       !--- dummy variables
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl), intent(in) ::CtrlParam

       !--- local variables
        type(SimMDBox)::trackSimBox
        integer::I, J, IB, NBOX, IP, IA, NA, NPRT0, NPRT, SHIFT, NI, IS, IS0
        integer, dimension(:),allocatable::hSITE, hSTAT, hSITEi
        real(KINDDF)::hSitePOS(2,3), APOS(1,3), VECT(3), LATT, BS(3), HBS(3),BL(3),BU(3), LAM
        character*256::GFILE1

  !-------------
          call Initialize_Globle_Variables_DEV(SimBox, CtrlParam)
          call Initialize_NeighboreList_DEV(SimBox, CtrlParam)
          call Cal_NeighBoreList_DEV(SimBox, CtrlParam)

           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)
           NPRT  = SimBox(1)%NPRT
           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           BL    = hm_RefSimBox%BOXLOW
           BU    = hm_RefSimBox%BOXUP
           HBS   = C_HALF*BS

         !$$--- calculate the occupation of lattice sites
           allocate(hSITE(NPRT*NBOX), hSTAT(NPRT0*NBOX), hSITEi(NPRT0))
           call Cal_Occupied_Voronoi_DEV(SimBox(1), CtrlParam, hSITE)

          !$$--- to calculate occupation state
           hSTAT  = 0
           IP     = 0
           do IB=1, NBOX
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT
                 IP = IP + 1
                 if(hSITE(IP)>0) then
                    hSTAT(hSITE(IP)+SHIFT) = hSTAT(hSITE(IP)+SHIFT) + 1
                 end if
              end do
           end do

           IP = 0
           do IB=1, NBOX
              !$$--- to extract the number of position of interstitial
              NI = 0
              NA = 0
              hSITEi  = 0
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT0
                 IP = IP + 1
                 if(hSTAT(IP) .ge. 2) then
                    !$$--- check if the lattice site already on the list
                     if(any(hSITEi(1:NI) .eq. I)) cycle

                    !$$--- this site containing interstitial atoms
                    NI = NI + 1
                    hSITEi(NI) = I

                 end if
              end do !--- end loop for search LSIA

              !$$--- to place an atom in one of the randomly slected LSIA on the
               if(NI .gt. 0) then
                  IS0 = int(DRAND32()*dble(NI)) + 1
                  IS0 = hSITEi(IS0)

               else
                  IS0 = DRAND32()*NPRT0 + 1
               end if
               IS  = DRAND32()*m_SiteNeighbor%KVOIS(IS0) + 1
               IS  = m_SiteNeighbor%INDI(IS0, IS)
               hSitePOS(1,1:3) = hm_RefSimBox%XP(IS0, 1:3)
               hSitePOS(2,1:3) = hm_RefSimBox%XP(IS, 1:3)
               VECT(1:3) = hSitePOS(2,1:3)  - hSitePOS(1,1:3)
               if(dabs(VECT(1)) .gt. HBS(1)) VECT(1) = VECT(1) - DSIGN(BS(1),VECT(1))
               if(dabs(VECT(2)) .gt. HBS(2)) VECT(2) = VECT(2) - DSIGN(BS(2),VECT(2))
               if(dabs(VECT(3)) .gt. HBS(3)) VECT(3) = VECT(3) - DSIGN(BS(3),VECT(3))
               LAM = DRAND32()
               APOS(1,1) = hSitePOS(1,1) + VECT(1)*LAM
               APOS(1,2) = hSitePOS(1,2) + VECT(2)*LAM
               APOS(1,3) = hSitePOS(1,3) + VECT(3)*LAM

               if(APOS(1,1) .lt. BL(1) ) APOS(1,1) = APOS(1,1) + BS(1)
               if(APOS(1,1) .gt. BU(1) ) APOS(1,1) = APOS(1,1) - BS(1)
               if(APOS(1,2) .lt. BL(2) ) APOS(1,2) = APOS(1,2) + BS(2)
               if(APOS(1,2) .gt. BU(2) ) APOS(1,2) = APOS(1,2) - BS(2)
               if(APOS(1,3) .lt. BL(3) ) APOS(1,3) = APOS(1,3) + BS(3)
               if(APOS(1,3) .gt. BU(3) ) APOS(1,3) = APOS(1,3) - BS(3)


               call AddAtoms_SimMDBox(SimBox(IB), N=1, ITYPE=SimBox(IB)%NGROUP, RXP=APOS, TYPEORDER=1)

           end do !--- end loop for  boxes

           deallocate(hSITE, hSTAT, hSITEi)
          return
  end subroutine Insert_NewAtom
  !**********************************************************************************

  !**********************************************************************************
  subroutine Insert_Diffusor(SimBox, CtrlParam, RESTART)
  !***  PORPOSE: to initialize the configure of the box by emdedding an atom
  !              in a box
  !     INPUT:  SimBox0, the original substarte
  !             SimBox,  the box array to be created
  !             CtrlParam, the control parameter
  !             RESTART, indictor to indicating if restarting a new session
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   use RAND32_MODULE

   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer                     ::RESTART

  !--- Local variables
       integer::I

           !---
           if(RESTART.eq. 0) then
              call Insert_NewAtom(SimBox, CtrlParam)
           else
              do I=1, RESTART
                 call AddOneAtom_SimBoxArray(SimBox,ITYPE=SimBox(1)%NGROUP)
              end do
           end if
       return
  end subroutine Insert_Diffusor
  !****************************************************************************************

 !****************************************************************************************
 end module HCP_SIAc_GPU_


 !****************************************************************************************
  Program Atom_Interface_GPU_Main
  use MD_SimBoxArray_AppShell_16_GPU
  !---- If user-define potential to be used, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
  !     use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
  !     use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table

  use HCP_SIAc_GPU_
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord(PRERECORD=MyInitialize)
       !---- If user-define potential to be generated, modify the code as following examplee
            ! call APPSHELL_Main(numprocs,procid, FORCETABLE=Reg1, POTTYPE="EAM_TYPE", INICONFIGPROC=IniConfig_EMBEDMENT)
       !---  else use internal potentials
       call APPSHELL_Main(numprocs,procid, INICONFIGPROC=Insert_Diffusor)

       stop
  end  program Atom_Interface_GPU_Main
