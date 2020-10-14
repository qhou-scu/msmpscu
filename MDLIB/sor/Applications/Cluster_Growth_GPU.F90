 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to simulate cluster growth.
 !                  The growth of a cluster is simulated by place a new atom in the neigbor region of
 !                  the cluster. The Wigner-Seitz cells occupied by the cluster are firstly calcualted,
 !                  Then, the new atoms will be added to the WS sites that are a few WS site away from the
 !                  cluster. The module RefVoronoiVacancy_GPU to be used to identify the occupied WS cells.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       Cfg_Track_Cluster_GPU,F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  by HOU Qing, Oct., 2015
 !
  module GROWTH_CLUSTER_GPU_
  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use Cfg_Track_Cluster_GPU, only:&
                                  mp_SITE_VAC,       &
                                  mp_SITE_HOSTATOM,  &
                                  mp_SITE_IMPURITY,  &
                                  mp_SITE_COMBCLUS,  &
                                  mp_SITE_SIA,       &
                                  m_MXAI,            &
                                  m_SiteNeighbor,    &
                                  hm_RefSimBox,      &

                                  BefOneTest_Cfg_Track_Cluster, &
                                  AftOneTest_Cfg_Track_Cluster, &
                                  Clear_Cfg_Track_Cluster,      &
                                  Record_Cfg_Track_Cluster
  implicit none
         integer, private::m_DEFCTTYP   = 0
         integer, private::m_INSERTTYP  = 0
         integer, private::m_EXTSHELLS  = 2
         integer, private::m_INSSHELLS  = 0
  contains

  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate memories, and create the
   !               neighbor-list of the reference lattice
   !
   use Cfg_Track_Cluster_GPU, only:Intilize_Cfg_Track_Cluster, Get_InputTag
   use MD_TYPEDEF_InputPaser, only:Get_InputStatements
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

   !----   Local variable
    integer::LINE, N, I, ITYP, IERR
    character*256::STR
    character*64 ::INTAG
    character*32::STRNUMB(mp_MXGROUP)
    type(InputStatements)::INPUTS

         call Intilize_Cfg_Track_Cluster(SimBox, CtrlParam)
         call Get_InputTag(INTAG)
         call GetExInput_SimMDCtrl(CtrlParam, INTAG, INPUTS, IERR)

         call Get_InputStatements("&EXTSHELL", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
               m_EXTSHELLS = ISTR(STRNUMB(1))
         end if

         call Get_InputStatements("&INSSHELL", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
               m_INSSHELLS = ISTR(STRNUMB(1))
         end if

         call Get_InputStatements("&INSTYP", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
               m_INSERTTYP  = ISTR(STRNUMB(1))
         end if

         m_DEFCTTYP = 0
         call Get_InputStatements("&DEFECTTYP", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,mp_MXGROUP,N,STRNUMB)
            do I=1, N
               ITYP   = ISTR(STRNUMB(I))
               m_DEFCTTYP = ior(m_DEFCTTYP,2**(ITYP-1))
            end do
            if(m_DEFCTTYP .gt. 0) then
               m_DEFCTTYP = ior(m_DEFCTTYP,2**(SimBox%NGROUP + mp_SITE_COMBCLUS))
            end if

            call Extract_Substr(STR,mp_MXGROUP,N,STRNUMB)
            do I=1, N
               call UPCASE(STRNUMB(I))
               select case ( STRNUMB(I)(1:len_trim(STRNUMB(I)) ) )
                     case ("VAC")
                           ITYP =  SimBox%NGROUP + 1 + mp_SITE_VAC
                     case ("SIA")
                           ITYP =  SimBox%NGROUP + 1 + mp_SITE_SIA
                     case ("COMBCLUS")
                           ITYP =  SimBox%NGROUP + 1 + mp_SITE_COMBCLUS
               end select
               m_DEFCTTYP = ior(m_DEFCTTYP,2**(ITYP-1))
            end do
         end if
         !--- check the completeness of input
          if(m_INSERTTYP .eq. 0) then
             write(*,*) "MDPSCU Warning: the type of insert diffor is missed"
             write(*,*) "                add in your input file: &INSTYP t "
             write(*,*) "                t > 0 for insert an atom of type x"
             write(*,*) "                t < 0 for insert a vacancy"
             call ONWARNING(gm_OnWarning)
          end if
          call Release_InputStatements(INPUTS)
         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine CreateExtendedCluster(SimBox, CtrlParam, ExtShell, FillShell, Mask)
  !***  DESCRIPTION: to create extended cluster. An extended cluster is the cluster with
  !                  margin added to a cluster define by m_DEFCTTYP
  !
  !    INPUT:
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulatio
  !            ExtShell,    the number of shell to be exteneded
  !            ExtShell,    the number of shell to used to insert new atoms
  !
  !    OUTPUT:
  !            Mask,        the marker which indicating the statu of WS cells
  !                         Mask > 0, the site belong to a defect cluster
  !                         Mask = 0, the site dose not belog to any cluster
  !
  use MD_Globle_Variables_GPU,   only:Initialize_Globle_Variables_DEV, Clear_Globle_Variables_DEV, hm_ITYP, hm_XP
  use MD_NeighborsList_GPU,      only:Initialize_NeighboreList_DEV,Cal_NeighBoreList_DEV, Clear_NeighboreList_DEV
  use RefVoronoiVacancy_GPU,     only:Cal_Site_State_DEV
  implicit none
       !--- dummy variables
       type(SimMDBox), dimension(:), intent(in) :: SimBox
       type(SimMDCtrl),              intent(in) :: CtrlParam
       integer,                      intent(in) :: ExtShell, FillShell
       integer,dimension(:),         intent(out):: Mask

       !--- local variables
        integer:: NBOX, NPRTR, NG0, IB, IP, NC, NI, I, J, ITYP, JJ, L, K, IND, INS
        integer, dimension(:), allocatable::SITESTAT, ATOMID
       !$$--- debug box
        !$$type(SimMDBox)::swapSimBox
        !$$type(MDRecordStamp):: LSTAMP
       !-------------

           NPRTR  = hm_RefSimBox%NPRT
           NBOX   = size(SimBox)
           NG0    = SimBox(1)%NGROUP
          !--- allocate memeory to store the type of sites and the atoms on the sites
           allocate(SITESTAT(NPRTR*NBOX), ATOMID(NPRTR*NBOX*m_MXAI))

           !*** to initialize the GPU variable
           !    NOTE:
           call Initialize_Globle_Variables_DEV(SimBox, CtrlParam)
           call Initialize_NeighboreList_DEV(SimBox, CtrlParam)
           call Cal_NeighBoreList_DEV(SimBox, CtrlParam)
           call Cal_Site_State_DEV(SimBox(1), CtrlParam, SITESTAT, ATOMID)

           Mask = 0
           !--- identiy defect sites defined by m_DEFCTTYP
           do J=1, NPRTR*NBOX
              if(SITESTAT(J) .lt. 1) then ! for vacancy defect
                 if(iand(m_DEFCTTYP, 2**(NG0 + mp_SITE_VAC)) .gt. 0) &
                    Mask(J) = 1
              else if(SITESTAT(J) .eq. 1) then ! for impurity
                    ITYP = hm_ITYP( AtomID((J-1)*m_MXAI + 1) )
                    if(iand(2**(ITYP-1),m_DEFCTTYP) .gt. 0) Mask(J) = 1
              else if(SITESTAT(J) .ge. 2) then
                      NC = 0
                      NI = 0
                      do I=(J-1)*m_MXAI+1, J*m_MXAI
                         if(AtomID(I) .eq. 0) exit
                         ITYP = hm_ITYP( AtomID(I) )

                         if(iand(2**(ITYP-1), m_DEFCTTYP) .gt. 0 ) then
                            NC = NC + 1     ! number impurity atoms on the site
                         else
                            NI = NI + 1     ! number of SIA on the site
                         end if
                      end do
                      if( NI .gt. 0 .and. NC .gt.0 ) then ! for combined site
                          if(iand(m_DEFCTTYP, 2**(NG0 + mp_SITE_COMBCLUS)) .gt. 0) Mask(J) = 1
                      else
                         if(NC .eq. 0) then  ! for pure SIA cluster
                            if(iand(m_DEFCTTYP, 2**(NG0 + mp_SITE_SIA)) .gt. 0) Mask(J) = 1
                         else                ! for pure impurity cluster
                            Mask(J) = 1
                         end if
                      end if
              end if
           end do

          !--- to extend the site clusters
           do IB=1, NBOX
              do L=1, ExtShell
                 do J=(IB-1)*NPRTR+1, IB*NPRTR
                    if(Mask(J) .eq. L) then
                       JJ = J - NPRTR*(IB-1)
                       do K=1, m_SiteNeighbor%KVOIS(JJ)
                          IND = m_SiteNeighbor%INDI(JJ, K) + NPRTR*(IB-1)
                          if(Mask(IND) .le. 0) then
                             Mask(IND) = Mask(J)+1
                          endif
                       end do
                    end if
                 end do
              end do
           end do

           !--- to make the filling region
           if(FillShell .gt. 0) then
              do IB=1, NBOX
                 do L=1, FillShell
                    do J=(IB-1)*NPRTR+1, IB*NPRTR
                       if(Mask(J) .gt. 0) then
                          JJ = J - NPRTR*(IB-1)
                          do K=1, m_SiteNeighbor%KVOIS(JJ)
                             IND = m_SiteNeighbor%INDI(JJ, K) + NPRTR*(IB-1)
                             if(Mask(IND) .le. 0) then
                                Mask(IND) = -1
                             endif
                          end do
                       end if
                    end do
                 end do
              end do

           else if(FillShell .eq. 0) then !--- the defect sites to be set to accepting new atoms
              do IB=1, NBOX
                 do J=(IB-1)*NPRTR+1, IB*NPRTR
                    if(Mask(J) .gt. 0) then
                       Mask(J) = -1
                    end if
                 end do
              end do

           else if(FillShell .lt. 0) then !--- all site outof extended clusters to be set to region accepting new atoms
              do J=1, NPRTR*NBOX
                 if(Mask(J) .eq. 0) Mask(J) = -1
              end do
           end if

           deallocate(SITESTAT, ATOMID)
           call Clear_NeighboreList_DEV()
           call Clear_Globle_Variables_DEV()

           !$$--- debug
           !$$do IB=1, NBOX
           !$$   call Copy_SimMDBox(hm_RefSimBox, swapSimBox)
           !$$   swapSimBox%ITYP(1:NPRTR) = Mask((IB-1)*NPRTR+1:IB*NPRTR) + 2
           !$$   LSTAMP%ITest = 1
           !$$   LSTAMP%IBox =  1
           !$$   LSTAMP%ICfg  = 1
           !$$  call Putout_Instance_Config_SimMDBox("DebugCfg", swapSimBox, LSTAMP)
          !$$end do

          return
  end subroutine CreateExtendedCluster
  !**********************************************************************************

  !**********************************************************************************
  subroutine Insert_NewAtom(SimBox, CtrlParam)
  !***  DESCRIPTION: to insert a new atom to the box
  !
  !    INPUT:  SimBox, the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use RAND32_MODULE
  use RefVoronoiVacancy_GPU,   only:hm_RefSimBox
  use Cfg_Track_Cluster_GPU,   only:m_SiteNeighbor

  implicit none
       !--- dummy variables
       type(SimMDBox), dimension(:) ::SimBox
       type(SimMDCtrl), intent(in)  ::CtrlParam

       !--- local variables
        integer::J, J0, J1, IB, NBOX, IP, IA, NPRT0, NPRT, IS0, NS, INS
        integer, dimension(:),allocatable::MASK
        real(KINDDF)::hSitePOS(2,3), APOS(1,3), VECT(3), LATT, BS(3), HBS(3),BL(3),BU(3), LAM

       !--------
           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)
           NPRT  = SimBox(1)%NPRT

           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           BL    = hm_RefSimBox%BOXLOW
           BU    = hm_RefSimBox%BOXUP
           HBS   = C_HALF*BS

           INS = m_INSERTTYP
           allocate(MASK(NPRT0*NBOX))
           call CreateExtendedCluster(SimBox, CtrlParam, m_EXTSHELLS, m_INSSHELLS, MASK)

           !--- have defect clusters
           do IB=1, NBOX
              NS  = count(Mask((IB-1)*NPRT0+1:IB*NPRT0) .lt. 0)
              if(NS .gt. 0) then   !--- there are preexisting defect cluster
                 IS0 = DRAND32()*NS + 1
                 IP  = 0
                 do J=(IB-1)*NPRT0+1, IB*NPRT0
                    if(MASK(J) .lt. 0) then
                       IP = IP + 1
                       if(IP .eq. IS0) then
                          J0 =  J - (IB-1)*NPRT0
                          J1  = DRAND32()*m_SiteNeighbor%KVOIS(J0) + 1
                          J1  = m_SiteNeighbor%INDI(J0, J1)
                          hSitePOS(1,1:3) = hm_RefSimBox%XP(J0,  1:3)
                          hSitePOS(2,1:3) = hm_RefSimBox%XP(J1, 1:3)
                          VECT(1:3)       = hSitePOS(2,1:3)  - hSitePOS(1,1:3)
                          if(dabs(VECT(1)) .gt. HBS(1)) VECT(1) = VECT(1) - DSIGN(BS(1),VECT(1))
                          if(dabs(VECT(2)) .gt. HBS(2)) VECT(2) = VECT(2) - DSIGN(BS(2),VECT(2))
                          if(dabs(VECT(3)) .gt. HBS(3)) VECT(3) = VECT(3) - DSIGN(BS(3),VECT(3))
                          LAM = DRAND32()*C_HALF
                          APOS(1,1) = hSitePOS(1,1) + VECT(1)*LAM
                          APOS(1,2) = hSitePOS(1,2) + VECT(2)*LAM
                          APOS(1,3) = hSitePOS(1,3) + VECT(3)*LAM

                          if(APOS(1,1) .lt. BL(1) ) APOS(1,1) = APOS(1,1) + BS(1)
                          if(APOS(1,1) .gt. BU(1) ) APOS(1,1) = APOS(1,1) - BS(1)
                          if(APOS(1,2) .lt. BL(2) ) APOS(1,2) = APOS(1,2) + BS(2)
                          if(APOS(1,2) .gt. BU(2) ) APOS(1,2) = APOS(1,2) - BS(2)
                          if(APOS(1,3) .lt. BL(3) ) APOS(1,3) = APOS(1,3) + BS(3)
                          if(APOS(1,3) .gt. BU(3) ) APOS(1,3) = APOS(1,3) - BS(3)
                          call AddAtoms_SimMDBox(SimBox(IB), N=1, ITYPE=INS, RXP=APOS, TYPEORDER=1)
                       end if
                    end if
                 end do
              else !--- there is no preexisting cluster, place the new atom around origin
                 APOS(1,1) = (DRAND32() - 0.5D0)*LATT
                 APOS(1,2) = (DRAND32() - 0.5D0)*LATT
                 APOS(1,3) = (DRAND32() - 0.5D0)*LATT
                 call AddAtoms_SimMDBox(SimBox(IB), N=1, ITYPE=INS, RXP=APOS, TYPEORDER=1)
              end if
           end do !--- end loop for  boxes
           deallocate(MASK)
          return
  end subroutine Insert_NewAtom
  !**********************************************************************************

  !**********************************************************************************
  subroutine Insert_NewVac(SimBox, CtrlParam)
  !***  DESCRIPTION: to insert a new vacancy to the box
  !
  !    INPUT:  SimBox, the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use RAND32_MODULE
  use RefVoronoiVacancy_GPU,   only:hm_RefSimBox
  use Cfg_Track_Cluster_GPU,   only:m_SiteNeighbor

  implicit none
       !--- dummy variables
       type(SimMDBox),  dimension(:)::SimBox
       type(SimMDCtrl), intent(in)  ::CtrlParam

       !--- local variables
        integer::J, J0, IB, NBOX, IP,  NPRT0, NPRT, IS0, NS, INS, IA, K
        integer, dimension(:),allocatable::MASK
        real(KINDDF)::hSitePOS(3), LATT, BS(3), HBS(3),BL(3),BU(3), R2, R2MX, VECT(3)

       !--------
           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)
           NPRT  = SimBox(1)%NPRT

           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           BL    = hm_RefSimBox%BOXLOW
           BU    = hm_RefSimBox%BOXUP
           HBS   = C_HALF*BS

           INS = iabs(m_INSERTTYP)
           allocate(MASK(NPRT0*NBOX))
           call CreateExtendedCluster(SimBox, CtrlParam, m_EXTSHELLS, m_INSSHELLS, MASK)

           !--- have defect clusters
           do IB=1, NBOX
              NS  = count(Mask((IB-1)*NPRT0+1:IB*NPRT0) .lt. 0)
              if(NS .gt. 0) then   !--- there are preexisting defect cluster
                 IS0 = DRAND32()*NS + 1
                 IP  = 0
                 do J=(IB-1)*NPRT0+1, IB*NPRT0
                    if(MASK(J) .lt. 0) then
                       IP = IP + 1
                       if(IP .eq. IS0) then
                          J0 =  J - (IB-1)*NPRT0
                          hSitePOS(1:3) = hm_RefSimBox%XP(J0, 1:3)
                       end if
                    end if
                 end do
              else !--- there is no preexisting cluster, place the vacancy the site close est to the origin
                 hSitePOS(1:3) = 0.D0
              end if

              !--- to find out the atoms closest to hSitePOS
              R2MX  = 1.D64
              IA    = 0
              do J=1, NPRT0
                 if(SimBox(IB)%ITYP(J) .eq. INS) then
                    VECT(1:3) = SimBox(IB)%XP(J,1:3) - hSitePOS(1:3)
                    do K=1, 3
                       if(dabs(VECT(K)) .GT. HBS(K)) then
                          VECT(K) = VECT(K) - DSIGN(BS(K),VECT(K))
                       end if
                    end do
                    if(sum(VECT*VECT) .lt. R2MX) then
                       R2MX = sum(VECT*VECT)
                       IA = J
                    end if
                 end if
              end do
              if(IA .gt. 0) then
                 call DelAtom_SimMDBox(SimBox(IB), IA)
              end if
           end do !--- end loop for  boxes

           deallocate(MASK)
          return
  end subroutine Insert_NewVac
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
              if(m_INSERTTYP .gt. 0) then
                 call Insert_NewAtom(SimBox, CtrlParam)
              else if(m_INSERTTYP .lt. 0) then
                 call Insert_NewVac(SimBox, CtrlParam)
              end if
           else
              if(m_INSERTTYP .gt. 0) then
                 do I=1, RESTART
                    call AddOneAtom_SimBoxArray(SimBox,ITYPE=SimBox(1)%NGROUP)
                 end do
              else if(m_INSERTTYP .lt. 0) then
                 do I=1, RESTART
                    call DeleteOneAtom_SimBoxArray(SimBox, iabs(m_INSERTTYP))
                 end do
             end if
           end if
       return
  end subroutine Insert_Diffusor
  !****************************************************************************************

 !****************************************************************************************
 end module GROWTH_CLUSTER_GPU_


 !****************************************************************************************
  Program GROWTH_CLUSTER_GPU_Main
  use MD_SimBoxArray_AppShell_16_GPU
  use GROWTH_CLUSTER_GPU_
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord( PRERECORD =MyInitialize,                 &
                                BEFONETEST=BefOneTest_Cfg_Track_Cluster, &
                                AFTONETEST=AftOneTest_Cfg_Track_Cluster, &
                                AFTRECORD= Clear_Cfg_Track_Cluster,      &
                                RECORDPROC=Record_Cfg_Track_Cluster)
       call APPSHELL_Main(numprocs,procid, INICONFIGPROC=Insert_Diffusor)

       stop
  end  program GROWTH_CLUSTER_GPU_Main
