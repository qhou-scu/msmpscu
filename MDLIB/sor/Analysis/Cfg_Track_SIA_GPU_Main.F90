 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract the track of self-interstitials.
 !                  A reference lattice will be loaded in, the localtion of interstitials are thus
 !                  identified by the Wigner_Seize cells of the reference lattice in which more than
 !                  one atoms exist. In addition to the interstitial atoms, the neighbores of the
 !                  interstitial atoms, and the neighbour lattces are to be output.
 !
 !                  This program can be used to depict the trajectories of vacancies and interstitials
 !                  as well as the projectile produced by cascade collision of energtic projectil in materials.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       RefVoronoiVA_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one line needs to be add:
 !
 !                    &AUXF_RVVIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_RVVOUT filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                  &REFCFG        indication of the filename of the reference configurqation.
 !                                  if the reference configuration is not assigned explicitly,
 !                                  the configuration at zero time step will be the reference
 !                                  configuration.
 !
 !                                  usage:   &REFCFG filename
 !
 !                                  example: &REFCFG "Myreference.crystal"
 !
 !                  If the filename of the reference configuration is given in the control file described
 !                  above, the third file storing the reference configuration needs to be prepaired. The
 !                  the file could be in the format of "&CFGXYZ", "&BOXCFG14" or "t-x-y-z" file whithout
 !                  header.
 !
 !                  With the input file(s) are ready, the program can be run on command line:
 !
 !                  Cfg_Track_SIA_GPU.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-11 (Hou Qing, Sichuan university)
 !
 !

 module Cfg_Track_SIA_GPU
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use MD_NeighborsList

 implicit none
         integer, private::m_NREC                                        ! number of actually record
         character(len=256),  private::m_OUTPATH =""                     ! the path output data
         character(len=256),  private::m_TRACKPATH =""                   ! the path output tracking data
         character(len=256),  private::m_CLUSTPATH =""                   ! the path output clusteringing data
         integer,             private::m_OUPUTREF  = 1                   ! >0 output reference lattice
                                                                         ! otherwise, not output reference lattice

         integer, parameter, private::mp_MXAI =4                         ! the permitted number of atoms occupy one site

         integer, parameter, private::mp_SITE_SIA = 1
         integer, parameter, private::mp_SITE_NEIB_SIA = 2
         integer, parameter, private::mp_SITE_NEIB_NORMAL = 3
         integer, parameter, private::mp_SITE_SIA_ATOM = 4
         integer, parameter, private::mp_SITE_NEIHB_SITEATOM = 5
         integer, parameter, private::mp_SITE_VAC = 6
         integer, parameter, private::mp_SITE_NEIB_VAC = 7
         integer, parameter, private::mp_PROJ_ATOM = 8
         integer, parameter, private::mp_HOST_LATT = 9

         type(SimMDBox), dimension(:), allocatable::m_trackDEF
         type(SimMDBox), dimension(:), allocatable::m_trackPROJ    ! the tracks of projectiles
         integer, private::m_PROJTYP      = 0                           ! the type projectile atoms
         integer, private::m_MXNEARSITE   = 70
         type(NEIGHBOR_LIST)::m_SiteNeighbor

  contains
  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to Initialize the calculation by allocate memories, and create the
   !               neighbor-list of the reference lattice
   !
   use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables
   use RefVoronoiVacancy_GPU, only:Initialize_ReferenceVacancy_DEV, Get_filename_Output, Get_InputTag, hm_RefSimBox
   use MD_NeighborsList_GPU, only:Reorder_NeighBoreList_Nearest_DEV, COPYOUT_NeighboreList_DEV
   use MD_TYPEDEF_InputPaser
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
  integer::N, I, LINE, IERR
  character*256::STR
  character*64::INTAG
  character*32::STRNUMB(5)
  type(MDRecordStamp)::STAMP
  type(InputStatements)::INPUTS

         !$$--- loading control parameters
         call Initialize_ReferenceVacancy_DEV(SimBox, CtrlParam)
         call Get_InputTag(INTAG)
         call GetExInput_SimMDCtrl(CtrlParam, INTAG, INPUTS, IERR)

         call Get_InputStatements("&NEARSITE", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
            m_MXNEARSITE = ISTR(STRNUMB(1))
            if(m_MXNEARSITE .le. 0) then
               write(*,fmt="(A,I4,A)")  ' MDPSCU Error: max number of nearest sites is zero in '//gm_ExeName(1:len_trim(gm_ExeName))
               write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
               write(*,fmt="(A)")       '               Process to be stopped'
               stop
            end if
         end if

         call Get_InputStatements("&PROJTYP", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
            m_PROJTYP = ISTR(STRNUMB(1))
         end if

         call Get_InputStatements("&OUTPUTLATT", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
            m_OUPUTREF = ISTR(STRNUMB(1))
         end if

         call Get_filename_Output(m_OUTPATH)
         call GetPath(m_OUTPATH, m_OUTPATH)
         m_TRACKPATH = m_OUTPATH(1:len_trim(m_OUTPATH))//"TRACK/"
         m_CLUSTPATH = m_OUTPATH(1:len_trim(m_OUTPATH))//"CLUST/"
         call CreateDataFolder_Globle_Variables(m_TRACKPATH)
         call CreateDataFolder_Globle_Variables(m_CLUSTPATH)

         write(*, fmt="(A, 1x, I4)")     " !    The max number of nearest sites is: ", m_MXNEARSITE
         if(m_PROJTYP .gt. 0) then
            write(*, fmt="(A, 1x, I4)")  " !    The type ID of projectile is: ", m_PROJTYP
         end if

         !$$--- We keep the neighbor list of the reference lattice
         call Reorder_NeighBoreList_Nearest_DEV(m_MXNEARSITE)
         call COPYOUT_NeighboreList_DEV(m_SiteNeighbor, ORDER=0)

         !$$--- put out the reference box for the convience of post analysis
         !$$    NOTE: the particles in the reference box have been partitioned
         !$$           on calling Initialize_ReferenceVacancy_DEV
          STAMP%ITime =0
          STAMP%ISect =1
          STAMP%ICfg  =0
          STAMP%IRec  =0
          STAMP%Time  =0
          call Putout_Instance_Config_SimMDBox(m_CLUSTPATH(1:len_trim(m_CLUSTPATH))//"RefCfg", hm_RefSimBox, STAMP)
          call Putout_NeighboreList(m_CLUSTPATH(1:len_trim(m_CLUSTPATH))//"RefNNList", m_SiteNeighbor)

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine Defect_Track(Stamp, SimBox, CtrlParam, hSTAT, hOCCUa)
  !***  DESCRIPTION: to identify and putout the SIA track. It is assumed the
  !                  the Cal_Site_State_DEV has be called.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !            hSTAT,       the occupation state of the reference sites
  !            hOCCUa,      the ID of atoms occupying the sites
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT,hm_XP
  use RefVoronoiVacancy_GPU, only:hm_RefSimBox
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       integer, dimension(:),        intent(in)::hSTAT, hOCCUa

       !--- local variables
        type(SimMDBox)::TRACKBOX
        integer, dimension(:), allocatable::hSITEi, hOCCUai
        integer::I, II, J, K, IB, NBOX, IP, IA, NA, NPRT0, NPRT, SHIFT, NI, IS, ISECT, ICFG, IFDEF
        real(KINDDF)::hSitePOS(1,3), hSiteVECT(1,3), APOS(1,3), AVEC(1,3), LATT, BS(3), HBS(3)
        character*256::GFILE1

  !-------------
           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(hSTAT)/NPRT0
           NPRT  = dm_NPRT/NBOX
           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           HBS   = C_HALF*BS

          !--- calculate the occupation of lattice sites
           allocate(hSITEi(NPRT0), hOCCUai(NPRT))

          !--- to extract the position of interstitial
           IP = 0
           do IB=1, NBOX
              call Release_SimMDBox(TRACKBOX)
              call CopyInformation_SimMDBox(hm_RefSimBox, TRACKBOX)
              TRACKBOX%NA = 0
              TRACKBOX%NGROUP = size(TRACKBOX%NA)
              NI = 0
              NA = 0
              hSITEi  = 0
              hOCCUai = 0
              SHIFT = (IB-1)*NPRT0

              do I=1, NPRT0
                 IP = IP + 1
                 if(hSTAT(IP) .ne. 1) then
                    !--- this site containing interstitial atoms
                    NI = NI + 1
                    hSITEi(NI) = I
                    hSitePOS(1,1:3) = hm_RefSimBox%XP(I, 1:3)
                    if(hSTAT(IP) .eq. 0) then
                       call AddAtoms_SimMDBox(TRACKBOX, 1,  ITYPE=mp_SITE_VAC, RXP=hSitePOS(1:1,1:3), TYPEORDER=0)
                    else
                       call AddAtoms_SimMDBox(TRACKBOX, 1,  ITYPE=mp_SITE_SIA, RXP=hSitePOS(1:1,1:3), TYPEORDER=0)
                    end if
                    !--- the atoms on the site
                    do J=1, mp_MXAI
                       IA = hOCCUa((IP-1)*mp_MXAI+J)
                       if(IA .gt. 0 .and. .not.any(hOCCUai(1:NA) .eq. IA)) then
                          NA = NA + 1
                          hOCCUai(NA) = IA
                          APOS(1,1:3) = hm_XP(IA, 1:3)
                          AVEC(1,1:3) = hSitePOS(1,1:3) - hm_XP(IA, 1:3)
                          if(dabs(AVEC(1,1)) .GT. HBS(1))  then
                              AVEC(1,1) = AVEC(1,1) - DSIGN(BS(1),AVEC(1,1))
                              APOS(1,1) = hSitePOS(1,1) + AVEC(1,1)
                          end if
                          if(dabs(AVEC(1,2)) .GT. HBS(2))  then
                             AVEC(1,2) = AVEC(1,2) - DSIGN(BS(2),AVEC(1,2))
                             APOS(1,2) = hSitePOS(1,2) + AVEC(1,2)
                          end if
                          if(dabs(AVEC(1,3)) .GT. HBS(3))  then
                             AVEC(1,3) = AVEC(1,3) - DSIGN(BS(3),AVEC(1,3))
                             APOS(1,3) = hSitePOS(1,3) + AVEC(1,3)
                          end if
                          AVEC = AVEC/LATT
                          call AddAtoms_SimMDBox(TRACKBOX, 1, ITYPE=mp_SITE_SIA_ATOM, &
                                                      RXP=APOS(1:1,1:3), RXP1=AVEC(1:1,1:3), TYPEORDER=0)
                       end if
                    end do
                 end if
              end do !--- end looking for LDEFECT

              !--- looking for the neigboring sites of LDEFECT
              do II=1, NI
                 I = hSITEi(II)
                    !--- add the neigbor lattices
                    do K=1, m_SiteNeighbor%KVOIS(I)
                       IS = m_SiteNeighbor%INDI(I, K)
                       hSitePOS(1,1:3) = hm_RefSimBox%XP(IS, 1:3)
                       if(any(hSITEi(1:NI) .eq. IS)) then
                          hSiteVECT(1,1:3) = hm_RefSimBox%XP(I, 1:3) -  hSitePOS(1,1:3)
                          if(dabs(hSiteVECT(1,1)) .GT. HBS(1)) then
                             hSiteVECT(1,1) = hSiteVECT(1,1) - DSIGN(BS(1),hSiteVECT(1,1))
                          end if

                          if(dabs(hSiteVECT(1,2)) .GT. HBS(2)) then
                             hSiteVECT(1,2) = hSiteVECT(1,2) - DSIGN(BS(2),hSiteVECT(1,2))
                          end if

                          if(dabs(hSiteVECT(1,3)) .GT. HBS(3)) then
                             hSiteVECT(1,3) = hSiteVECT(1,3) - DSIGN(BS(3),hSiteVECT(1,3))
                          end if
                          hSiteVECT = hSiteVECT/LATT
                          if(hSTAT(SHIFT + IS) .eq. 0) then
                             call AddAtoms_SimMDBox(TRACKBOX, 1, ITYPE=mp_SITE_NEIB_VAC, RXP=hSitePOS(1:1,1:3), &
                                                         RXP1=hSiteVECT(1:1,1:3), TYPEORDER=0)
                          else
                             call AddAtoms_SimMDBox(TRACKBOX, 1, ITYPE=mp_SITE_NEIB_SIA, RXP=hSitePOS(1:1,1:3), &
                                                         RXP1=hSiteVECT(1:1,1:3), TYPEORDER=0)
                          end if
                       else
                          call AddAtoms_SimMDBox(TRACKBOX, 1, ITYPE=mp_SITE_NEIB_NORMAL, RXP=hSitePOS(1:1,1:3), &
                                                      TYPEORDER=0)
                       end if

                       do J=1, mp_MXAI
                          IA = hOCCUa((SHIFT+IS-1)*mp_MXAI+J)
                          if(IA .gt. 0 .and. .not.any(hOCCUai(1:NA) .eq. IA)) then
                             NA = NA + 1
                             hOCCUai(NA) = IA
                             APOS(1,1:3) = hm_XP(IA, 1:3)
                             AVEC(1,1:3) = hSitePOS(1,1:3) - hm_XP(IA, 1:3)
                             if(dabs(AVEC(1,1)) .GT. HBS(1))  then
                                AVEC(1,1) = AVEC(1,1) - DSIGN(BS(1),AVEC(1,1))
                                APOS(1,1) = hSitePOS(1,1) + AVEC(1,1)
                             end if
                             if(dabs(AVEC(1,2)) .GT. HBS(2))  then
                                AVEC(1,2) = AVEC(1,2) - DSIGN(BS(2),AVEC(1,2))
                                APOS(1,2) = hSitePOS(1,2) + AVEC(1,2)
                             end if
                             if(dabs(AVEC(1,3)) .GT. HBS(3))  then
                                AVEC(1,3) = AVEC(1,3) - DSIGN(BS(3),AVEC(1,3))
                                APOS(1,3) = hSitePOS(1,3) + AVEC(1,3)
                             end if
                             AVEC = AVEC/LATT
                             call AddAtoms_SimMDBox(TRACKBOX, 1, ITYPE=mp_SITE_NEIHB_SITEATOM, &
                                                     RXP=APOS(1:1,1:3), RXP1=AVEC(1:1,1:3), TYPEORDER=0)

                          end if
                       end do
                    end do
              end do  !--- end loop for site reference sites
              !--- merge to tracking box
              !--- to put out the configure
               !call GetSectID_SimMDCtrl(CtrlParam, ITIME, ISECT)
               !call GetCfgID_SimMDCtrl(CtrlParam, ITIME, ICFG)
               call Merge_SimMDBox(TRACKBOX, m_trackDEF(IB))

               if(m_PROJTYP .gt. 0) then
                  call Merge_SimMDBox(m_trackPROJ(IB), TRACKBOX)
               end if

               call STRCATI(GFILE1, m_TRACKPATH(1:len_trim(m_TRACKPATH)), "P", 0, 4)
               call STRCATI(GFILE1, GFILE1, "_", (Stamp%ITest-1)*NBOX+IB, 4)
               call Putout_Instance_Config_SimMDBox(GFILE1, TRACKBOX, Stamp)
               call STRCATI(GFILE1, GFILE1, ".", Stamp%IRec(1), 4)
               write(*,fmt="(A, I6, A)") " box ",(Stamp%ITest-1)*NBOX+IB, "-> "//GFILE1(1:len_trim(GFILE1))

           end do !--- end loop for  boxes
           call Release_SimMDBox(TRACKBOX)
           deallocate(hSITEi, hOCCUai)
          return
  end subroutine Defect_Track
  !**********************************************************************************

  !**********************************************************************************
  subroutine Defect_Site_Clustering(hSTAT, hCNUM, hSQNC)
  !***  DESCRIPTION: to identify the connected cluster shape. It is assumed the
  !                  the Cal_Site_State_DEV has be called.
  !
  !    INPUT:  hSTAT,       the occupation state of the reference sites
  !
  !    OUTPUT: hCNUM        the total number of clusters
  !            hSQNC,       the squence of cluster, the integer of each element in the arry
  !                         is the site ID of the reference lattice. A cluster is ended with zero,
  !                         the clusters in a box is ended by -1
  !
  !            NOTE:       the arraies hSQNC and hSTYP are allocted in this routine
  !
  use RefVoronoiVacancy_GPU, only:hm_RefSimBox
  implicit none
       !--- dummy variables
       integer, dimension(:), intent(in)::hSTAT
       integer::hCNUM
       integer, dimension(:)::hSQNC, hOCCUP

       !--- local variables
        integer, dimension(:), allocatable::VISITED, SQNC
        integer::I, II, J, K, IB, NBOX, NN, NPRT0, SHIFT, NI, IS, STEP, IP, NI0

  !-------------
           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(hSTAT)/NPRT0

          !--- calculate the occupation of lattice sites
           allocate( VISITED(NPRT0), SQNC(NPRT0))

          !--- to extract the position of interstitial
           hCNUM = 0
           NI0 = 0
           do IB=1, NBOX
              NI = 0
              VISITED = 0
              SQNC    = -1
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT0
                 if(VISITED(I) .gt. 0) cycle

                 if(hSTAT(SHIFT+I) .ne. 1) then
                    !$$--- this site containing interstitial atoms or vacancy
                    NI = NI + 1
                    SQNC(NI)   = I
                    VISITED(I)  = 1
                    IP = NI
                    STEP = 0
                    do while(IP .le. NPRT0)

                       II = SQNC(IP)
                       if(II .lt. 0) exit

                       if(hSTAT(SHIFT + II) .ne. 1) then
                          do K=1, m_SiteNeighbor%KVOIS(II)
                             IS = m_SiteNeighbor%INDI(II, K)
                             if(VISITED(IS) .gt. 0) cycle

                             STEP  = STEP + 1
                             SQNC(NI + STEP)  = IS
                             VISITED(IS)  = 1
                          end do
                       end if
                       IP = IP + 1
                    end do
                    NI = NI + STEP + 1
                    SQNC(NI) = 0
                    hCNUM = hCNUM + 1
                 end if
              end do !--- end looking for LDEFECT

              hSQNC(NI0+1:NI0+NI) = SQNC(1:NI)
              NI0 = NI0+NI+1
              hSQNC(NI0) = -1

           end do !--- end loop for  boxes
           deallocate(VISITED, SQNC)
          return
  end subroutine Defect_Site_Clustering
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to
  !                  interfaced to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the
  !                  the neighbor list routine has be called.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use MD_SimBoxArray
  use RefVoronoiVacancy_GPU, only:Cal_Site_State_DEV, hm_RefSimBox
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
        integer:: NBOX, NPRT0, I, IP, hCNUM
        real(KINDDF)::POS(1,3), VECT(1,3), BS(3), HBS(3)
        integer, dimension(:),allocatable::hSTAT, hOCCUa, hSQNC

       !-------------
           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)

           if(Stamp%ITime .lt. 0) then
              m_NREC = 0
              !--- write out the header for iste clustering
              call Record_Site_Clustering(Stamp, SimBox, CtrlParam)
              return
           end if

           if(.not.allocated(m_trackDEF) ) then
              allocate(m_trackDEF(NBOX))
           end if

           if(.not.allocated(m_trackPROJ) .and. m_PROJTYP .gt. 0) then
              allocate(m_trackPROJ(NBOX))
           end if

           if(m_NREC .eq. 0) then
              call Release_SimBoxArray(m_trackDEF)
              if(m_PROJTYP .gt. 0) then
                 call Release_SimBoxArray(m_trackPROJ)
                 do I=1, NBOX
                    call CopyInformation_SimMDBox(SimBox(I), m_trackPROJ(I))
                    m_trackPROJ(I)%NA = 0
                    m_trackPROJ(I)%NGROUP = size(m_trackPROJ(I)%NA)
                 end do
              end if
           end if

          !---  prepair the trajectoies of the projectiles
           if(m_PROJTYP .gt. 0) then
              do I=1, NBOX
                 BS    = hm_RefSimBox%ZL
                 HBS   = C_HALF*BS
                 do IP=1, SimBox(I)%NPRT
                    if(SimBox(I)%ITYP(IP) .eq. m_PROJTYP) then

                      if(m_NREC .eq. 0) then
                         VECT(1,1:3) = 0
                         POS(1,1:3) = SimBox(I)%XP(IP,1:3)
                      else
                         POS(1,1:3) = SimBox(I)%XP(IP,1:3)
                         VECT(1,1:3) = POS(1,1:3) - m_trackPROJ(I)%XP(m_trackPROJ(I)%NPRT, 1:3)
                         if(dabs(VECT(1,1)) .GT. HBS(1))  then
                            VECT(1,1) = VECT(1,1) - DSIGN(BS(1),VECT(1,1))
                         end if
                         if(dabs(VECT(1,2)) .GT. HBS(2))  then
                            VECT(1,2) = VECT(1,2) - DSIGN(BS(2),VECT(1,2))
                         end if
                         if(dabs(VECT(1,3)) .GT. HBS(3))  then
                            VECT(1,3) = VECT(1,3) - DSIGN(BS(3),VECT(1,3))
                         end if
                         m_trackPROJ(I)%XP1(m_trackPROJ(I)%NPRT, 1:3) = VECT(1,1:3)/SimBox(I)%RR
                      end if
                      call AddAtoms_SimMDBox(m_trackPROJ(I), 1, ITYPE=mp_PROJ_ATOM, RXP=POS(1:1,1:3), RXP1=VECT(1:1,1:3), &
                                                         TYPEORDER=0)
                    end if
                 end do
              end do
           end if

          !--- allocate memeory to store the STAT of Sites and the atoms on the sites
           allocate(hSTAT(NPRT0*NBOX), hOCCUa(NPRT0*NBOX*mp_MXAI))

          !--- to calculate occupation state
           write(*,fmt="(A,A)") ' MDPSCU Message: the occupation state to be calculated'
           call  Cal_Site_State_DEV(SimBox(1), CtrlParam, hSTAT, hOCCUa)
           call Defect_Track(Stamp, SimBox, CtrlParam, hSTAT, hOCCUa)

           allocate(hSQNC(NPRT0*NBOX))
           call Defect_Site_Clustering(hSTAT, hCNUM, hSQNC)
           call Record_Site_Clustering(Stamp, SimBox, CtrlParam, hSTAT, hSQNC)

           deallocate(hSTAT, hOCCUa)
           if(allocated(hSQNC)) deallocate(hSQNC)

           m_NREC = m_NREC + 1
          return
  end subroutine MyRECORD
  !**********************************************************************************

  !**********************************************************************************
  subroutine Record_Site_Clustering(Stamp, SimBox, CtrlParam, hSTAT, hSQNC)
  !***  DESCRIPTION: to putout the clustering results to output file.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_SimBoxArray
  use RefVoronoiVacancy_GPU, only:Cal_Site_State_DEV, Get_filename_RefCfg, hm_RefSimBox
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,            intent(in):: Stamp
       type(SimMDBox), dimension(:),    intent(in):: SimBox
       type(SimMDCtrl),                 intent(in):: CtrlParam
       integer, dimension(:), optional, intent(in):: hSTAT, hSQNC
       !--- local variables
        integer:: NBOX, NPRT0, IB, IP, ITYP, CID, hFile, NSIA, SHIFT,JOB,ITIME, ICFG
        real(KINDDF)::POS(3), TIME
        character*256::GFILE1


  !-------------
           NPRT0  = hm_RefSimBox%NPRT
           NBOX   = size(SimBox)
           JOB    = Stamp%ITest
           ITIME  = Stamp%ITime
           TIME   = Stamp%Time
           ICFG   = Stamp%IRec(1)

           if(ITIME .lt. 0 .or. .not.present(hSQNC) .or. .not.present(hSTAT)) then
              do IB=1, NBOX
               call STRCATI(GFILE1, m_CLUSTPATH(1:len_trim(m_CLUSTPATH)), "P", 0, 4)
               call STRCATI(GFILE1, GFILE1, "_", (JOB-1)*NBOX+IB, 4)
               call AvailableIOUnit(hFile)
               open(UNIT=hFile, file = GFILE1)
               write(*,fmt="(A, I6, A)") " clustering result to be output for box ",(JOB-1)*NBOX+IB, "-> "//GFILE1(1:len_trim(GFILE1))

               write(hFile, fmt="(A, I8)")    '!--- THE CLUSTERING SEQUENCE VS TIME CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
               write(hFile, fmt="(A)")        '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
               write(hFile, fmt="(A)")        '!    AUTHOR: HOU Qing'
               write(hFile, fmt="(A)")        '!    '
               write(hFile, fmt="(A, I8)")    '&LATT_CLUST_SEQ'

               write(hFile, fmt="(A,1X,I8)")            "!--- REFERENCE LATTICE: "
               call Get_filename_RefCfg(GFILE1)
               write(hFile, fmt="(A,1X,I8)")            "&REFCFG     '"//m_CLUSTPATH(1:len_trim(m_CLUSTPATH))//"RefCfg"
               write(hFile, fmt="(A,1X,I8)")            "&NATOM      ", hm_RefSimBox%NPRT
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE    ", hm_RefSimBox%ZL/hm_RefSimBox%RR
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW     ", hm_RefSimBox%BOXLOW/hm_RefSimBox%RR
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT (A)   ", hm_RefSimBox%RR*CP_CM2A

               write(hFile, fmt="(A,1X,I8)")  "!--- JOB AND BOX ID: "
               write(hFile, fmt="(A, I8)")    '&JOB ID:      ', JOB
               write(hFile, fmt="(A, I8)")    '&BOX ID:      ', IB
               close(hFile)

              end do
              return
           end if

           IP = 0
           do IB=1, NBOX
               SHIFT = (IB-1)*NPRT0
               call STRCATI(GFILE1, m_CLUSTPATH(1:len_trim(m_CLUSTPATH)), "P", 0, 4)
               call STRCATI(GFILE1, GFILE1, "_", (JOB-1)*NBOX+IB, 4)
               call AvailableIOUnit(hFile)
               open(UNIT=hFile, file = GFILE1, POSITION='APPEND')
               write(hFile, fmt="(A, I8)")
               write(hFile, fmt="(A, I8)")      '&ICFG:       ',  ICFG
               write(hFile, fmt="(A, I8)")      '&ITIME:      ',  ITIME
               write(hFile, fmt="(A, 1PE12.6)") '&TIME (ps):  ',  TIME
               write(hFile, fmt="(A11, 1x, A9, 1x, 3(A12,1x), A10, 1x, A10)") '&STARTSQN' , 'SITE-TYPE', 'POSX(LATT)', 'POSY', 'POSZ', 'LATTCE-ID', 'SIAs'
               do while(.true.)
                  IP = IP + 1

                  if(hSQNC(IP) .gt. 0) then
                    POS(1:3) = hm_RefSimBox%XP(hSQNC(IP), 1:3)
                    NSIA     = hSTAT(hSQNC(IP)+SHIFT)
                    if(NSIA .gt. 1) then
                       ITYP = mp_SITE_SIA
                    else if(NSIA .le. 0) then
                       ITYP = mp_SITE_VAC
                    else
                       ITYP = mp_SITE_NEIB_NORMAL
                    end if
                    write(hFile, fmt="(11x, 1x, I9, 1x, 3(1PE12.4,1x), I10, 1x, I10)") ITYP, POS(1:3)/hm_RefSimBox%RR, hSQNC(IP), NSIA

                  else if(hSQNC(IP) .eq. 0) then
                    write(hFile, fmt="(11x, 1x, I9, 1x, 3(1PE12.4,1x), I10)")   0

                  else if(hSQNC(IP) .eq. -1) then
                     exit
                  end if
               end do
               write(hFile, fmt="(A, I8)")  '&ENDSQN'
               close(hFile)
               write(*,fmt="(A, I6, A)") " clustering result output to "//GFILE1(1:len_trim(GFILE1))
           end do

          return
  end subroutine Record_Site_Clustering
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRecord_FROM_TRACK(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to
  !                  interfaced to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the
  !                  the neighbor list routine has be called.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    NOTE:  The difference between this and MyRECORD is that the configurations
  !           are loaded from  the track file generated by DEFECT_TRACK.
  !
  !           This routine can only be used for SIA indentification.
  !
  !
  use MD_SimBoxArray
  use MD_Globle_Variables_GPU, only:dm_NPRT, hm_ITYP
  use RefVoronoiVacancy_GPU, only:Cal_Occupied_Voronoi_DEV, hm_RefSimBox
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        integer:: NBOX, NPRT0, IB, I, IP, SHIFT, NPRT, hCNUM
        real(KINDDF)::POS(1,3), VECT(1,3), BS(3), HBS(3)
        integer, dimension(:),allocatable::hSTAT, hSQNC
        integer, dimension(:),allocatable::hSITE

  !-------------

           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)
           if(NBOX .gt. 1) then
               write(*,fmt="(A)") "MDPSCU Error: RECORD_FROM_TRACK in Cfg_Track_SIA_GPU supports only one-TEST-one-BOX"
               write(*,fmt="(A)") "              Process to be stopped"
               stop
           end if

           if(Stamp%ITime .lt. 0) then
              !--- write out the header
              call Record_Site_Clustering(Stamp, SimBox, CtrlParam)
              m_NREC = 0
              return
           end if

           NPRT  = dm_NPRT/NBOX
           if(m_NREC .eq. 0) then
           end if

         !--- calculate the occupation of lattice sites
           allocate(hSTAT(NPRT0*NBOX))

          !--- to calculate occupation state
           write(*,fmt="(A,A)") ' MDPSCU Message: the occupation state to be calculated'
           allocate(hSITE(dm_NPRT))
           call Cal_Occupied_Voronoi_DEV(SimBox(1), CtrlParam, hSITE)

           hSTAT  = 1  !--- This lead to more occupations on SIA sites
           IP     = 0
           do IB=1, NBOX
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT
                 IP = IP + 1
                 select case(hm_ITYP(IP))
                        case (mp_SITE_SIA_ATOM)
                              hSTAT(hSITE(IP)+SHIFT) = hSTAT(hSITE(IP)+SHIFT) + 1
                              !print *, 'SIA', hSITE(IP)+SHIFT
                        case (mp_SITE_VAC)
                              hSTAT(hSITE(IP)+SHIFT) =0
                              !print *, 'VAC', hSITE(IP)+SHIFT
                 end select
              end do
           end do

           !--- Note: SIA sites have been added more than once
           where(hSTAT .gt. 1)
                 hSTAT = hSTAT - 1
           endwhere
           deallocate(hSITE)

          !print *, count(hSTAT .eq. 1), count(hSTAT .eq. 2), count(hSTAT .eq. 0)
          !----
           allocate(hSQNC(NPRT0*NBOX))
           call Defect_Site_Clustering(hSTAT, hCNUM, hSQNC)
           call Record_Site_Clustering(Stamp, SimBox, CtrlParam, hSTAT, hSQNC)

          !--- to put out the cluster configuration

           deallocate(hSTAT)
           if(allocated(hSQNC)) deallocate(hSQNC)

           m_NREC = m_NREC + 1
          return
  end subroutine MyRecord_FROM_TRACK
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyPutoutSummary(ITEST)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use RefVoronoiVacancy_GPU, only:hm_RefSimBox
   implicit none
       !--- dummy variables
       integer::ITEST
       !--- local variables
       character*256::GFILE1
       integer::IB, NBOX
       type(MDRecordStamp)::STAMP
       type(SimMDBox)::REFLATT, TRACKBOX

              !--- to put out the configure
               NBOX = size(m_trackDEF)
               if(m_OUPUTREF .gt. 0) then
                    call Copy_SimMDBox(hm_RefSimBox, REFLATT)
                    REFLATT%ITYP = mp_HOST_LATT
               end if

               do IB=1, NBOX
                  call STRCATI(GFILE1, m_TRACKPATH(1:len_trim(m_TRACKPATH))//"merged_", "P", 0, 4)
                  call STRCATI(GFILE1, GFILE1, "_", (ITEST-1)*NBOX+IB, 4)
                  write(*,fmt="(A, I6, A)") " Merged track configuration saved to  "//GFILE1(1:len_trim(GFILE1))//".0000"
                  STAMP%ITest = ITEST
                  STAMP%ITime = 0
                  STAMP%ISect = 1
                  STAMP%ICfg  = 0
                  STAMP%IRec  = 0
                  STAMP%Time  = 0
                  call Copy_SimMDBox(m_trackDEF(IB), TRACKBOX)
                  if(m_OUPUTREF .gt. 0) call Merge_SimMDBox(REFLATT, TRACKBOX)
                  call Putout_Instance_Config_SimMDBox(GFILE1, TRACKBOX, STAMP)
              end do
              call Release_SimMDBox(TRACKBOX)
              call Release_SimMDBox(REFLATT)
          return
  end subroutine MyPutoutSummary
  !**********************************************************************************


  !**********************************************************************************
  subroutine MyAfterOneTest(SimBox0, Stamp, SimBox, CtrlParam)
  !***  PURPOSE:  conduct after one test
  !    INPUT:
  !               SimBox0,    the simulation box, useless here
  !               Stamp,      the record stamp
  !               SimBox,     the simulation box array, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
    implicit none
    type(SimMDBox),              intent(in) ::SimBox0
    type(MDRecordStamp),         intent(in) ::Stamp
    type(SimMDBox), dimension(:)            ::SimBox
    type(SimMDCtrl)                         ::CtrlParam
    integer                                 ::ITEST

          call MyPutoutSummary(Stamp%ITest)
          m_NREC = 0

  end subroutine MyAfterOneTest
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
  use MD_SimBoxArray
  use RefVoronoiVacancy_GPU, only:Clear_WorkingArray_DEV
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam


          if(allocated(m_trackDEF)) then
             call Release_SimBoxArray(m_trackDEF)
             deallocate(m_trackDEF)
          end if

          if(allocated(m_trackPROJ)) then
             call Release_SimBoxArray(m_trackPROJ)
             deallocate(m_trackPROJ)
          end if

          call Clear_NeighboreList(m_SiteNeighbor)
          call Clear_WorkingArray_DEV(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module Cfg_Track_SIA_GPU


 !**********************************************************************************
 Program Cfg_Track_SIA_GPU_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use Cfg_Track_SIA_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD=MyInitialize,            &
                                AFTONETEST=MyAfterOneTest,         &
                                AFTRECORD= MyCleaner,              &
                                RECORDPROC=MyRECORD)
                                !RECORDPROC=MyRecord_FROM_TRACK)
       call Main_ANALYSIS(0)

       stop
 End program Cfg_Track_SIA_GPU_Main
