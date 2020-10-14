 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract the track of clusters.
 !                  A reference lattice will be loaded in, the sites that contain interstitials and
 !                  projectiles are allocated. A few kind of sites are defined:
 !
 !                  1. A SIA sites are identified by the Wigner_Seize cells containing more than one than
 !                     one substrate atoms.
 !                  2. A impuritie site is defined as a site containing one or more impurity atoms.
 !                  3. A combined site is defined as a site containing at least one impurity atom and
 !                     one substrate atom.
 !                  4. A vacancy site is defined as a site containing no atom.
 !
 !                  The connectivity between the sites are checked. A site cluster is defined as connected
 !                  if they are cloestly neighboring. Two sites are neighboring when they are Wigner-Seitz
 !                  cells share the same faces.
 !
 !                  This program can be used to depict the trajectories of clusters defined above.
 !
 !                  as well as the projectile produced by cascade collision of energtic projectil in materials.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       RefVoronoiVA_GPU.F90
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2017-04 (Hou Qing, Sichuan university)
 !
 !

 module Cfg_Track_Cluster_GPU
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use MD_NeighborsList,      only:NEIGHBOR_LIST
 use RefVoronoiVacancy_GPU, only:hm_RefSimBox

 implicit none
         integer,             private::m_NREC                            ! number of actually record
         character(len=256),  private::m_OUTPATH =""                     ! the path output data
         character(len=256),  private::m_TRACKPATH =""                   ! the path output tracking data
         character(len=256),  private::m_CLUSTPATH =""                   ! the path output clusteringing data


         integer::m_OUTMERGEREFB = 0
         integer::m_OUTVORVOL    = 0
         integer::m_MXAI =12                                             ! the permitted number of atoms occupy one site
         integer::m_MXCLUS = 5                                           ! the permitted number of seperated clusters in a box

         integer, parameter::mp_SITE_VAC      = 0
         integer, parameter::mp_SITE_HOSTATOM = 1
         integer, parameter::mp_SITE_IMPURITY = 2
         integer, parameter::mp_SITE_COMBCLUS = 3
         integer, parameter::mp_SITE_SIA      = 4

         integer, parameter::mp_Cluster_EndFlag = 0
         integer, parameter::mp_Box_EndFlag     = -1

         integer, private     ::m_PROJTYP  = 0                            ! the type of atoms to be included in clustering
         integer, private     ::m_DEFTYP   = 0                            ! the type of defect to be included in clustering

         !--- clustering method
         integer, parameter, private ::mp_DELAUNAY = 1
         integer, parameter, private ::mp_BOND     = 2
         integer           , private ::m_ClsMethod = mp_DELAUNAY          ! indcating the method used for clustering
                                                                          ! mp_DELAUNAY, use Voronoi tessellation, default
                                                                          ! mp_BOND, use bond lenght given by user
         real(KINDDF),       private ::m_BondLen   = 0.65                 ! bond length only when m_ClsMethod = mp_BOND


         type(NEIGHBOR_LIST)  ::m_SiteNeighbor

         integer, dimension(:),allocatable, private::m_hFiles             ! the IO unit for files to store the evolution of  number of defects
         integer, dimension(:),allocatable, private::m_SQNC               ! the cluster current sequences

         !--------------------------------------------------------------
         private::  PrintControlParameters

  contains

  !**********************************************************************************
   integer function GetPROJTYP()
     implicit none
            GetPROJTYP =  m_PROJTYP
            return
   end function GetPROJTYP
  !**********************************************************************************

  !****************************************************************************
  subroutine Get_InputTag(Tag)
  !***  PURPOSE:   to extract the input tag of output
  !
  !     INPUT:
  use RefVoronoiVacancy_GPU,   only:InputTag=>Get_InputTag
      implicit none
      !--- dummy variables
      character*(*)::Tag
           call InputTag(Tag)
           return
  end subroutine Get_InputTag
 !****************************************************************************

 !****************************************************************************
   subroutine Intilize_Cfg_Track_Cluster(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate memories, and create the
   !               neighbor-list of the reference lattice
   !
   use MD_Globle_Variables,     only:CreateDataFolder_Globle_Variables
   use MD_Globle_Variables_GPU
   use MD_NeighborsList,        only:Putout_NeighboreList
   use MD_NeighborsList_GPU
   use RefVoronoiVacancy_GPU,   only:Initialize_ReferenceVacancy_DEV, Get_filename_Output
   use VoronoiTessellation_GPU, only:Cal_Delaunay_NeighborList, Clear_VoronoiTessellation_DEV
   use MD_TYPEDEF_PrintList,    only:Add_PrintProcess

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
  integer::N, I, LINE, ITYP, IERR
  character*256::STR
  character*64 ::INTAG
  character*32::STRNUMB(mp_MXGROUP)
  type(MDRecordStamp)::STAMP
  type(InputStatements)::INPUTS
  type(SimMDCtrl)::tCtrlParam

         !$$--- loading control parameters and create reference box and the neighbore list
         call Initialize_ReferenceVacancy_DEV(SimBox, CtrlParam)
         call Get_InputTag(INTAG)
         call GetExInput_SimMDCtrl(CtrlParam, INTAG, INPUTS, IERR)

         m_PROJTYP = 0
         m_DEFTYP  = 0
         call Get_InputStatements("&PROJTYP", INPUTS, STR, LINE)
         if(LINE .le. 0) then
            call Get_InputStatements("&PROP_TYPE", INPUTS, STR, LINE)
         end if
         if(LINE .gt. 0) then
            call Extract_Numb(STR,mp_MXGROUP,N,STRNUMB)
            do I=1, N
               ITYP   = ISTR(STRNUMB(I))
               m_PROJTYP = ior(m_PROJTYP,2**(ITYP-1))
            end do
            !--- if atom type to be indentifed are given, impurity should be defect to be indentified
            if(m_PROJTYP .gt. 0) then
               m_DEFTYP = ior(m_DEFTYP,2**mp_SITE_IMPURITY)
               m_DEFTYP = ior(m_DEFTYP,2**mp_SITE_COMBCLUS)
            end if

            call Extract_Substr(STR,mp_MXGROUP,N,STRNUMB)
            do I=1, N
               call UpCase(STRNUMB(I))
               select case ( STRNUMB(I)(1:len_trim(STRNUMB(I)) ) )
                     case ("VAC")
                           ITYP =  SimBox%NGROUP + 1 + mp_SITE_VAC
                     case ("SIA")
                           ITYP =  SimBox%NGROUP + 1 + mp_SITE_SIA
                     case ("COMBCLUS")
                           ITYP =  SimBox%NGROUP + 1 + mp_SITE_COMBCLUS
               end select
               m_PROJTYP = ior(m_PROJTYP,2**(ITYP-1))
               m_DEFTYP  = ior(m_DEFTYP, 2**(ITYP - SimBox%NGROUP - 1))
            end do
         end if

         call Get_InputStatements("&MXATOMATSITE", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
            if(N.ge. 1) m_MXAI = ISTR(STRNUMB(1))
         else
            call Get_InputStatements("&MXNAPERSITE", INPUTS, STR, LINE)
            call Extract_Numb(STR,1,N,STRNUMB)
            if(N.ge. 1) m_MXAI = ISTR(STRNUMB(1))
         end if

         call Get_InputStatements("&MXCLUSTER", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Numb(STR,1,N,STRNUMB)
            if(N.ge. 1) m_MXCLUS = ISTR(STRNUMB(1))
         end if

         call Get_InputStatements("&CLUSMETH", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Substr(STR,2,n,STRNUMB)
            m_ClsMethod = 0
            if(N.ge.1) then
              call UpCase(STRNUMB(1))
              if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "USEDT") then
                 m_ClsMethod = mp_DELAUNAY
              else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "USEBD") then
                m_ClsMethod = mp_BOND
                call Extract_Numb(STR,1,n,STRNUMB)
                if(n .ge. 1) m_BONDLEN = DRSTR(STRNUMB(1))
              end if
            end if

            if(m_ClsMethod .le. 0) then
               write(*,fmt="(A,I4,A)")  ' MDPSCU Error: method to determine connectity is not given'
               write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
               write(*,fmt="(A, BZI6)") '        Usage: &CLUSMETH  "USEDT", '
               write(*,fmt="(A, BZI6)") '               to use Delaunay tessellation to determine connectity of defects'
               write(*,fmt="(A, BZI6)") '        Or:    &CLUSMETH  "USEBD" bondlen"'
               write(*,fmt="(A, BZI6)") '               to use bond length to determine connectity of defects'
               write(*,fmt="(A)")       '               Process to be stopped'
               stop
            end if
         end if

         call Get_InputStatements("&SAVE_VOL", INPUTS, STR, LINE)
         if(LINE .le. 0) call Get_InputStatements("&SAVE_VORVOL", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Substr(STR,1,N,STRNUMB)
            if(N.ge.1) then
              call UpCase(STRNUMB(1))
              if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "YES") then
                 m_OUTVORVOL = 1
              else
                m_OUTVORVOL = 0
              end if
            else
               call Extract_Numb(STR,1,N,STRNUMB)
               if(N.ge. 1) m_OUTVORVOL = ISTR(STRNUMB(1))
            end if

         end if

         call Get_InputStatements("&SAVE_REFBOX", INPUTS, STR, LINE)
         if(LINE .le. 0) call Get_InputStatements("&SAVE_RB", INPUTS, STR, LINE)
         if(LINE .gt. 0) then
            call Extract_Substr(STR,1,N,STRNUMB)
            if(N.ge.1) then
              call UpCase(STRNUMB(1))
              if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "YES") then
                 m_OUTMERGEREFB = 1
              else
                 m_OUTMERGEREFB = 0
              end if
            else
               call Extract_Numb(STR,1,N,STRNUMB)
               if(N.ge. 1) m_OUTMERGEREFB = ISTR(STRNUMB(1))
            end if
         end if

         call Release_InputStatements(INPUTS)
         !--- after get control parameters

         call Get_filename_Output(m_OUTPATH)
         call GetPath(m_OUTPATH, m_OUTPATH)
         m_TRACKPATH = m_OUTPATH(1:len_trim(m_OUTPATH))//"TRACK/"
         m_CLUSTPATH = m_OUTPATH(1:len_trim(m_OUTPATH))//"CLUST/"
         call CreateDataFolder_Globle_Variables(m_TRACKPATH)
         call CreateDataFolder_Globle_Variables(m_CLUSTPATH)
         call CreateDataFolder_Globle_Variables(m_OUTPATH(1:len_trim(m_OUTPATH))//"VORONOI/")

         !$$--- create the neighbore list for clustering
         call Copy_SimMDCtrl(CtrlParam, tCtrlParam)
         if(m_ClsMethod .eq. mp_DELAUNAY) then
            !$$--- to create the Delaunay neighbor-list
            call Cal_Delaunay_NeighborList(hm_RefSimBox, tCtrlParam, m_SiteNeighbor)
            call Clear_VoronoiTessellation_DEV(hm_RefSimBox, tCtrlParam)
         else
            !$$--- to create the by distance
            tCtrlParam%NB_RM = m_BondLen*hm_RefSimBox%RR
            call Initialize_Globle_Variables_DEV(hm_RefSimBox, tCtrlParam)
            call Initialize_NeighboreList_DEV(hm_RefSimBox, tCtrlParam)
            call Cal_NeighBoreList_DEV(hm_RefSimBox, tCtrlParam)
            call Copyout_NeighboreList_DEV(m_SiteNeighbor, ORDER=1)
            call Clear_NeighboreList_DEV()
         end if
         call Clear_Globle_Variables_DEV()

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

          allocate(m_SQNC(hm_RefSimBox%NPRT*CtrlParam%MULTIBOX), m_hFiles(CtrlParam%MULTIBOX))
          m_NREC = 0
          m_hFiles = 0

          call Add_PrintProcess(PrintControlParameters)
         return
   end subroutine Intilize_Cfg_Track_Cluster
  !**********************************************************************************

  !**********************************************************************************
   subroutine PrintControlParameters(hFile, SimBox, CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
   implicit none
  !----   DUMMY Variables
   integer,         intent(in)::hFile
   type(SimMDBox),  intent(in)::SimBox
   type(SimMDCtrl), intent(in)::CtrlParam

  !----   Local variables
    integer::I
    character*256::str

         !$$--- print out the control parameters
          write(hFile,*) "!************ Cfg_Track_Cluster_GPU module to be used **********"
          write(hFile,*) "!    With the control parameters:"

          if(m_PROJTYP .gt. 0) then
             write(hFile, fmt="(A, 1x, 10(I4,1x))")  " !    Site type contained in clustering are....: "
             do I=1,SimBox%NGROUP
                if(iand(2**(I-1),m_PROJTYP) .gt. 0) then
                   write(hFile, fmt="(A, 1x, 10(I4,1x))") " !    atom type................................: " , I
                end if
             end do
             str = ""
             if(iand(2**(SimBox%NGROUP+mp_SITE_VAC),m_PROJTYP)      .ne. 0 ) str = str(1:len_trim(str))//" VACANCY, "
             if(iand(2**(SimBox%NGROUP+mp_SITE_SIA),m_PROJTYP)      .ne. 0 ) str = str(1:len_trim(str))//" SIA, "
             if(iand(2**(SimBox%NGROUP+mp_SITE_COMBCLUS),m_PROJTYP) .ne. 0 ) str = str(1:len_trim(str))//" COMBSITE "
             if(len_trim(str) .gt. 0) then
                str = " !    and......................................: "//str(1:len_trim(str))
                 write(hFile, fmt="(A, 1x, 10(I4,1x))")  str(1:len_trim(str))
             end if
          end if

          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !    Max number of atoms in each site.........: ", m_MXAI
          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !    Max number of clusters in each box.......: ", m_MXCLUS

          if(m_ClsMethod .eq. mp_DELAUNAY) then
          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !    Clutering method.........................: DT"
          else if(m_ClsMethod .eq. mp_BOND) then
          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !    Clutering method.........................: BOND"
          write(hFile, fmt="(A, 1x, F8.2,1x))")   " !         with bond-length (LU) ..............: ",m_BONDLEN
          end if

          if(m_OUTVORVOL .gt. 0) then
          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !    Save cluter envelopes....................: YES"
          else
          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !    Save cluter envelopes....................: NO"
          end if
          write(hFile, fmt="(A, 1x, 10(I6,1x))")  " !  "

         return
   end subroutine PrintControlParameters
  !**********************************************************************************

  !**********************************************************************************
  subroutine IdentifySites(SiteType, AtomID, Siteflag)
  !***  DESCRIPTION: to identify the type of sites.
  !
  !    INPUT:
  !            SiteType,    the stat of sites, indicating how many atoms and the types the sites contain
  !                         this array will be overwrited by new stat value which indicating the site
  !                         to be vacancy, or hostatom, or SIA or impurity
  !             AtomID,     the atom id in a site
  !
  !    OUTPUT: SiteType,    the type of sites
  !            Siteflag,    the flag indicating if the site to be included in clustering
  !
  use MD_Globle_Variables_GPU, only:hm_ITYP
  implicit none
       !--- dummy variables
       integer                          :: NSite
       integer, dimension(:)            :: SiteType
       integer, dimension(:), intent(in):: AtomID
       integer, dimension(:)            :: Siteflag
       !--- local variables
        integer:: NPRT0, J, I, ITYP, NC, NI

       !-------------

              NPRT0 = hm_RefSimBox%NPRT
              do J=1, NPRT0
                 if(SiteType(J) .lt. 1) then
                    SiteType(J) = mp_SITE_VAC

                 else if(SiteType(J) .eq. 1) then
                      ITYP = hm_ITYP( AtomID((J-1)*m_MXAI + 1) )
                      if(iand(2**(ITYP-1),m_PROJTYP) .ne. 0) then
                         SiteType(J) = mp_SITE_IMPURITY
                      else
                         SiteType(J) = mp_SITE_HOSTATOM
                      end if

                 else if(SiteType(J) .ge. 2) then
                      NC = 0
                      NI = 0
                      do I=(J-1)*m_MXAI+1, J*m_MXAI
                         if(AtomID(I) .eq. 0) exit

                         ITYP = hm_ITYP( AtomID(I) )

                         if(iand(2**(ITYP-1), m_PROJTYP) .ne. 0 ) then
                            NC = NC + 1
                         else
                            NI = NI + 1
                         end if

                      end do
                      if( NI .gt. 0 .and. NC .gt.0 ) then
                          SiteType(J) = mp_SITE_COMBCLUS
                      else
                        if(NC .eq. 0) then
                           SiteType(J) = mp_SITE_SIA
                        else
                           SiteType(J) = mp_SITE_IMPURITY
                        end if
                      end if
                 end if

                 if(iand(2**SiteType(J),m_DEFTYP) .gt. 0) then
                    Siteflag(J) = 1
                 else
                    Siteflag(J) = 0
                 end if
              end do
              !print *, "m_DEFTYP", m_DEFTYP, m_PROJTYP
              !print *, "SiteType0", count(SiteType(1:NPRT0) .eq. 0)
              !print *, "SiteType1", count(SiteType(1:NPRT0) .eq. 1)
              !print *, "SiteType2", count(SiteType(1:NPRT0) .eq. 2)
              !print *, "SiteType3", count(SiteType(1:NPRT0) .eq. 3)
              !print *, "SiteType4", count(SiteType(1:NPRT0) .eq. 4)
              !pause
          return
  end subroutine IdentifySites
  !**********************************************************************************

  !**********************************************************************************
  subroutine Defect_Site_Clustering(SiteType, Siteflag, NumC, Sqnc)
  !***  DESCRIPTION: to identify the connected cluster shape. It is assumed the
  !                  the Cal_Site_State_DEV has be called.
  !
  !    INPUT:  SiteType,     the type of sites
  !            Siteflag,     the flag indicating if the Site to be included in cluster
  !
  !    OUTPUT: hNumC         the total number of clusters
  !            hSqnc,        the squence of cluster, the integer of each element in the arry
  !                          is the site ID of the reference lattice. A cluster is ended with zero,
  !                          the clusters in a box is ended by -1
  !
  !
  implicit none
       !--- dummy variables
       integer, dimension(:), intent(in)::SiteType, Siteflag
       integer              ::NumC
       integer, dimension(:)::Sqnc

       !--- local variables
        integer, dimension(:), allocatable::VISITED
        integer::I, II, K, NPRT0, NI, IS, STEP, IP

       !-------------

           NPRT0 = hm_RefSimBox%NPRT
          !--- calculate the occupation of lattice sites
           allocate( VISITED(NPRT0))

          !--- to extract the position of interstitial
              NumC = 0
              NI = 0
              VISITED = 0
              Sqnc(1:NPRT0)   = mp_Box_EndFlag

              do I=1, NPRT0
                 if(VISITED(I) .gt. 0) cycle

                 !if(SiteType(I) .ne. mp_SITE_HOSTATOM) then
                 if(Siteflag(I) .gt. 0) then
                    !$$--- this site containing interstitial atoms or vacancy or impurity atom
                    NI        = NI + 1
                    Sqnc(NI)  = I

                    VISITED(I) = 1
                    IP         = NI
                    STEP       = 0
                    do while(IP .le. NPRT0)
                       II = Sqnc(IP)
                       if(II .lt. 0) exit

                       !if(SiteType(II) .ne. mp_SITE_HOSTATOM) then
                       if(Siteflag(II) .gt. 0) then
                          do K=1, m_SiteNeighbor%KVOIS(II)
                             IS = m_SiteNeighbor%INDI(II, K)
                             if(VISITED(IS) .gt. 0) cycle

                             !if(SiteType(IS) .ne. mp_SITE_HOSTATOM) then
                             if(Siteflag(IS) .gt. 0) then
                                STEP              = STEP + 1
                                VISITED(IS)       = 1
                                Sqnc(NI + STEP)  = IS
                             end if
                          end do
                       end if
                       IP = IP + 1
                    end do
                    NI       = NI + STEP + 1
                    NumC     = NumC + 1
                    Sqnc(NI) = mp_Cluster_EndFlag
                 end if
              end do

           deallocate(VISITED)
          return
  end subroutine Defect_Site_Clustering
  !**********************************************************************************

  !**********************************************************************************
  subroutine SeperationOfTwoDefectClusters(Nprt, List, Sqnc1, Sqnc2, Sep)
  !***  DESCRIPTION: to calculate the seperation between two clusters
  !                  the seperation is defined by the number of neighbor sites needed
  !                  to connect these two clsuter
  !
  !    INPUT:  Nprt,        the number of sites in the system
  !            List,        the neighbor list of the sites
  !            Sqnc1,       the site sequence of cluster #1
  !            Sqnc2,       the site sequence of cluster #2
  !
  !    OUTPUT:
  !            Sep          the seperation between the two cluster
  !
  !
  implicit none
       !--- dummy variables
       integer,              intent(in)::Nprt
       type(NEIGHBOR_LIST),  intent(in)::List
       integer, dimension(:),intent(in)::Sqnc1, Sqnc2
       integer                         ::Sep

       !--- local variables
        integer::I, J, ID1, ITY1, NN, IP, FLAG, LL, ISQNF, ISQNT
        integer, dimension(:), allocatable::SQN1, SFLAG


              allocate(SQN1(Nprt), SFLAG(Nprt) )

              SEP =  0
              FLAG = 0
              SFLAG = 0
              IP = 0

              !--- get size of cluster 1
              IP    = 1
              do while(Sqnc1(IP) .ne. mp_Cluster_EndFlag)
                 SQN1(IP)         = Sqnc1(IP)
                 SFLAG(Sqnc1(IP)) = 1
                 IP = IP + 1
              end do
              ISQNT = IP - 1
              !--- begining check the overlap of the two cluster
              ISQNF = 1
              LL = 0
              do while(.true.)
                 I = 1
                 do while(Sqnc2(I) .ne. mp_Cluster_EndFlag )
                    if(any(SQN1(ISQNF:ISQNT) .eq. Sqnc2(I)) ) then
                      FLAG = 1
                      exit
                   end if
                   I = I + 1
                 end do
                 if(FLAG .gt. 0) exit

                 LL = LL + 1
                 IP = ISQNT
                 do I = ISQNF, ISQNT
                    NN =  m_SiteNeighbor%KVOIS(SQN1(I))
                    do J=1, NN
                       ID1 = m_SiteNeighbor%INDI(SQN1(I),J)
                       if(SFLAG(ID1) .gt. 0) cycle
                       IP = IP + 1
                       SQN1(IP) = ID1
                       SFLAG(ID1) = 1
                    end do
                 end do

                ISQNF = ISQNT + 1
                ISQNT = IP
                if(ISQNF .gt. ISQNT) exit
              end do
              SEP = LL - 1
              deallocate(SQN1, SFLAG)

             return
  end subroutine SeperationOfTwoDefectClusters
  !**********************************************************************************

  !**********************************************************************************
  subroutine Putout_Site_Clustering(Stamp, SimBox, CtrlParam, NP, NC, ITYP, XP, VECT, SITESZ, LID, LTY, ICLS)
  !***  DESCRIPTION: to putout the clustering results to output file.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !            NP,          the number of points
  !            NC,          the number of clusters
  !            ITYP,        the type of points
  !            XP,          the position of the points
  !            VECT,        the vectors to locations
  !            SITESIZ,     the number of atoms contained in a site
  !            LID,         the ID of Wigner-Seitz conataing the points
  !            LTY,         the type reference lattice
  !            ICLS,        the cluster the point belong to
  !
  use MD_SimboxArray
  use RefVoronoiVacancy_GPU, only:Cal_Site_State_DEV, Get_filename_RefCfg
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,        intent(in):: Stamp
       type(SimMDBox),              intent(in):: SimBox
       type(SimMDCtrl),             intent(in):: CtrlParam
       integer,                     intent(in):: NP, NC
       integer, dimension(:),       intent(in):: SITESZ, ITYP, ICLS, LID, LTY
       real(KINDDF), dimension(:,:),intent(in):: XP, VECT
       !--- local variables
        integer:: NBOX, NPRT0, IB, I, hFile, JOB, NG
        real(KINDDF)::POS(3), TIME
        character*256::GFILE1


  !-------------
           JOB    = Stamp%ITest
           NG     = SimBox%NGROUP
           NBOX   = CtrlParam%MULTIBOX

           do IB  = Stamp%IBox(1), Stamp%IBox(2)
               call STRCATI(GFILE1, m_CLUSTPATH(1:len_trim(m_CLUSTPATH)), "P", 0, 4)
               call STRCATI(GFILE1, GFILE1, "_", (JOB-1)*NBOX+IB, 4)
               call STRCATI(GFILE1, GFILE1, ".", Stamp%IRec(1), 4)
               call AvailableIOUnit(hFile)
               open(UNIT=hFile, file = GFILE1)
               write(*,fmt="(A, I6, A)") " clustering result to be output for box ",(JOB-1)*NBOX+IB, "-> "//GFILE1(1:len_trim(GFILE1))

               write(hFile, fmt="(A, I8)")    '!--- THE  CLUSTER CONFIGURATION CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
               write(hFile, fmt="(A)")        '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
               write(hFile, fmt="(A)")        '!    AUTHOR: HOU Qing'
               write(hFile, fmt="(A)")        '!    '
               write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
               call Putout_RecordStamp(hFile, Stamp)
               write(hFile,fmt="(A)")           "!---  SYMBOL EXPLAINATION:"
               write(hFile,fmt="(A)")           "!     TYPECOL:   are the types of location"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", 1, " to type ", NG, " are the type of atoms"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", (NG+1)+mp_SITE_VAC,       " is location of vacancy"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", (NG+1)+mp_SITE_HOSTATOM,  " is location of single host atom"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", (NG+1)+mp_SITE_IMPURITY,  " is location of impurity atom(s)"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", (NG+1)+mp_SITE_COMBCLUS,  " is location of combinded cluster"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", (NG+1)+mp_SITE_SIA,       " is location of SIA"
               write(hFile,fmt="(A,I3,A,I3,A)") "!                type ", (NG+1)+mp_SITE_SIA + 1,   " is background lattice"
               write(hFile,fmt="(A)")           "!     XYZCOL:    are the coordinates for atoms and W-S cell locations"
               write(hFile,fmt="(A)")           "!     VECTCOL:   are vectors of atoms to the W-S cell location that the atoms in"
               write(hFile,fmt="(A)")           "!     SITESZ:    are number of atoms contained in the defect locations"
               write(hFile,fmt="(A)")           "!     LID:       are lattice ID that atoms and W-S cell locations belong to"
               write(hFile,fmt="(A)")           "!     LTY:       are lattice type that atoms and W-S cell locations belong to"
               write(hFile,fmt="(A)")           "!     CLSCOL:    are cluster ID that atoms and W-S cell locations belong to"
               write(hFile,fmt="(A)")           "!     "


               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE    ", hm_RefSimBox%ZL/hm_RefSimBox%RR
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW     ", hm_RefSimBox%BOXLOW/hm_RefSimBox%RR
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT       ", hm_RefSimBox%RR*CP_CM2A
               write(hFile, fmt="(A,1X,(I4,1X))")       "&TYPECOL    ", 1
               write(hFile, fmt="(A,1X,2(I4,1X),A)")    "&XYZCOL     ", 2
               write(hFile, fmt="(A,1X,2(I4,1X),A)")    "&VECTCOL    ", 5,  3, ', "D"'
               write(hFile, fmt="(A,1X,2(I4,1X),A)")    "&SITESZCOL  ", 8,  1, ', "I"'
               write(hFile, fmt="(A,1X,2(I4,1X),A)")    "&LIDCOL     ", 9,  1, ', "I"'
               write(hFile, fmt="(A,1X,2(I4,1X),A)")    "&LTYCOL     ", 10, 1, ', "I"'
               write(hFile, fmt="(A,1X,2(I4,1X),A)")    "&CLSCOL     ", 11, 1, ', "I"'
               write(hFile, fmt="(A,1X,(I4,1X))")       ""

               if(m_OUTMERGEREFB .gt. 0) then
                 write(hFile, fmt="(A,1X,I8, A)")       "&NATOM           ", NP + hm_RefSimBox%NPRT
               else
                 write(hFile, fmt="(A,1X,I8, A)")       "&NATOM           ", NP
               end if
               write(hFile, fmt="(A,1X,I8, A)")         "&NLATTICE        ", hm_RefSimBox%NPRT
               write(hFile, fmt="(A,1X,I8, A)")         "&NCLUS           ", NC
               write(hFile, fmt="(A,1X,I8, I8)")        "&ATOMTYP         ", 1, NG
               write(hFile, fmt="(A,1X,I8, A)")         "&SITETYP-VAC     ", (NG+1)+mp_SITE_VAC
               write(hFile, fmt="(A,1X,I8, A)")         "&SITETYP-HOSTATOM", (NG+1)+mp_SITE_HOSTATOM
               write(hFile, fmt="(A,1X,I8, A)")         "&SITETYP-IMPURITY", (NG+1)+mp_SITE_IMPURITY
               write(hFile, fmt="(A,1X,I8, A)")         "&SITETYP-COMBCLUS", (NG+1)+mp_SITE_COMBCLUS
               write(hFile, fmt="(A,1X,I8, A)")         "&SITETYP-SIA     ", (NG+1)+mp_SITE_SIA

               write(hFile, fmt="(A9, 1x, 6(A15,1x), A8, 1x, A8, 1x, A8, 1x, A8)") &
                                                                           '!TYPE' , 'POSX(LATT)', 'POSY', 'POSZ', &
                                                                           'VECX(LATT)', 'VECY', 'VECZ', 'SITESZ', &
                                                                           'LID', 'LTY','CLUSID'
               do I=1, NP
                  write(hFile, fmt="(I9, 1x, 6(1PE15.4,1x), I8,1x, I8,1x,I8,1x,I8))") &
                                                                          ITYP(I), XP(I,1:3)/hm_RefSimBox%RR, &
                                                                          VECT(I,1:3)/hm_RefSimBox%RR,        &
                                                                          SITESZ(I), LID(I), LTY(I), ICLS(I)
               end do

               if(m_OUTMERGEREFB .gt. 0) then
                  do I=1, hm_RefSimBox%NPRT
                     write(hFile, fmt="(I9, 1x, 6(1PE15.4,1x), I8,1x,I8,1x,I8,1x,I8))") &
                                                              (NG+1)+mp_SITE_SIA + 1,  hm_RefSimBox%XP(I,1:3)/hm_RefSimBox%RR, &
                                                              0.D0, 0.D0, 0.D0, 0, I, hm_RefSimBox%ITYP(I), 0
                 end do
               end if
               close(hFile)
           end do

          return
  end subroutine Putout_Site_Clustering
  !**********************************************************************************

  !**********************************************************************************
  subroutine Record_Cfg_Track_Cluster(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to
  !                  interfaced to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the
  !                  the neighbor list routine has be called.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use MD_Globle_Variables_GPU, only:hm_ITYP, hm_XP
  use RefVoronoiVacancy_GPU,   only:Cal_Site_State_DEV
  use VoronoiTessellation_GPU, only:Record_Voronoi_ClusterEnvelope_Type
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
        integer:: NBOX, NPRTR, NPRT0, NG0, NG1, IB, IP, NC, NC1, NP, NPS, IS, SHIFT, I, J, ERR
        real(KINDDF)::BS(3), HBS(3), SPOS(3)
        integer,      dimension(:),   allocatable::SITETYP, SITESIZ, ATOMID, SITEFLAG, SQNC, ITYP, ICLS, LID, LTY, CFLAG, CS
        integer,      dimension(:,:), allocatable::SEP
        real(KINDDF), dimension(:,:), allocatable::XP, VECT
        type(MDRecordStamp):: TSTAMP
        character*256::STRSWP
        character*16::NUMSTR

       !-------------

           if(Stamp%ITime .lt. 0) return

           NPRTR  = hm_RefSimBox%NPRT
           NBOX   = size(SimBox)
           NPRT0  =  SimBox(1)%NPRT
           NG0    =  SimBox(1)%NGROUP
           NG1    =  NG0 + 1
           call Copy_RecordStamp(Stamp, TSTAMP)
          !--- allocate memeory to store the type of sites and the atoms on the sites
           allocate(SITETYP(NPRTR*NBOX), ATOMID(NPRTR*NBOX*m_MXAI), SITEFLAG(NPRTR), SQNC(NPRTR), SITESIZ(NPRTR+NPRT0), &
                    ITYP(NPRTR+NPRT0), ICLS(NPRTR+NPRT0), XP(NPRTR+NPRT0,3), VECT(NPRTR+NPRT0,3), LID(NPRTR+NPRT0),     &
                    LTY(NPRTR+NPRT0), CFLAG(m_MXCLUS), SEP(m_MXCLUS, m_MXCLUS), CS(m_MXCLUS), stat=ERR )

           if(ERR .gt. 0) then
              write(*,fmt="(A)") "MDPSCU Error: fail to allocate memory in Record_Cfg_Track_Cluster"
              write(*,fmt="(A)") "              Process to be stopped"
              stop
           end if

           call Cal_Site_State_DEV(SimBox(1), CtrlParam, SITETYP, ATOMID)

           do IB=1, NBOX
              BS    = hm_RefSimBox%ZL
              HBS   = C_HALF*BS
              call IdentifySites(SITETYP((IB-1)*NPRTR+1:), ATOMID((IB-1)*NPRTR*m_MXAI+1:), SITEFLAG)
              call Defect_Site_Clustering(SITETYP((IB-1)*NPRTR+1:), SITEFLAG, NC, SQNC)
              if(NC .gt. m_MXCLUS) then
                 write(*,fmt="(A, 1x, I6)")        "MDPSCU Error: number of clusters larger than permitted value, in box #", (Stamp%ITest-1)*NBOX+IB
                 write(*,fmt="(A, 1x, I6, A, I6)") "              NCLUSE :", NC, " vs MXCLUSTER", m_MXCLUS
                 write(*,fmt="(A)")                "              Consider increase the value for keyword &MXCLUSTER."
                 write(*,fmt="(A)") "              Process to be stopped"
                 stop
              end if

              !--- calculate the number of sites and atoms to be output
              SHIFT   = (IB-1)*NPRTR
              NC      = 0
              NC1     = NC + 1
              NP      = 0
              CFLAG   = 1
              SITESIZ = 0
              do J=1, NPRTR
                 IS  = SQNC(J)
                 if(IS .eq. mp_Box_EndFlag) then
                    exit
                 end if

                 if(IS .eq. mp_Cluster_EndFlag) then
                    NC          = NC + 1
                    NC1         = NC + 1
                    CFLAG(NC+1) = J  + 1
                    cycle
                 end if

                 !---  add the site
                 NP         = NP + 1
                 ITYP(NP)   = SITETYP(SHIFT+IS) + NG1
                 ICLS(NP)   = NC1
                 LID(NP)    = IS
                 LTY(NP)    = hm_RefSimBox%ITYP(IS)
                 NPS        = NP

                 !--- NOTE: the lattice postion in reference box have been reordered
                 !          after call Neighborelist calculation
                 SPOS(1:3)   =  hm_RefSimBox%XP(IS,1:3)
                 XP(NP,1:3)  =  SPOS(1:3)
                 VECT(NP,1:3)= 0.D0

                 !--- add the atom in the site
                 do I=(SHIFT+IS-1)*m_MXAI+1, (SHIFT+IS)*m_MXAI
                    if(ATOMID(I) .eq. 0) exit

                    NP           = NP + 1
                    ICLS(NP)     = NC1
                    ITYP(NP)     = hm_ITYP( ATOMID(I) )
                    LID(NP)      = IS
                    LTY(NP)      = hm_RefSimBox%ITYP(IS)
                    XP(NP,1:3)   = hm_XP( ATOMID(I), 1:3)
                    VECT(NP,1:3) = SPOS(1:3) - XP(NP,1:3)
                    if(dabs(VECT(NP,1)) .GT. HBS(1))  then
                       VECT(NP,1) = VECT(NP,1) - DSIGN(BS(1),VECT(NP,1))
                    end if
                    if(dabs(VECT(NP,2)) .GT. HBS(2))  then
                        VECT(NP,2) = VECT(NP,2) - DSIGN(BS(2),VECT(NP,2))
                    end if
                    if(dabs(VECT(NP,3)) .GT. HBS(3))  then
                       VECT(NP,3) = VECT(NP,3) - DSIGN(BS(3),VECT(NP,3))
                    end if
                    SITESIZ(NPS) = SITESIZ(NPS) + 1
                 end do

              end do
              !--- to output the configures
              TSTAMP%NBox    = 1
              TSTAMP%IBox    = IB
              call Putout_Site_Clustering(TSTAMP, SimBox(IB), CtrlParam, NP, NC, ITYP, XP, VECT, SITESIZ, LID, LTY, ICLS)

              !--- put out the number of clusters
              !--- get the cluster size
              CS  = 0
              SEP = 0
              do I=1, NC
                 CS(I) = CFLAG(I+1) - CFLAG(I) - 1
              end do
              do I=1, NC-1
                 do J=I+1, NC
                    call SeperationOfTwoDefectClusters(NPRTR, m_SiteNeighbor, SQNC(CFLAG(I):), SQNC(CFLAG(J):), SEP(I,J))
                    SEP(J,I) = SEP(I,J)
                 end do
               end do

               write(NUMSTR,*)m_MXCLUS
               NUMSTR= adjustl(NUMSTR)
               STRSWP = "(2I12, 1x, 1PE12.5, 1x, A, I6, 1x, A," &
                         //NUMSTR(1:len_trim(NUMSTR))//"(I6,1x),A,"
               STRSWP = STRSWP(1:len_trim(STRSWP))//NUMSTR(1:len_trim(NUMSTR))//&
                        "("//NUMSTR(1:len_trim(NUMSTR))//"(I6,1x),A))"
               write(m_hFiles(IB), fmt=STRSWP) TSTAMP%IRec(1), TSTAMP%ITIME, TSTAMP%TIME, &
                                               ",",   NC, ",", (CS(I), I=1, m_MXCLUS), &
                                               "," , ((SEP(I,J), I=1, m_MXCLUS), ",", J=1, m_MXCLUS)

           end do  !--- end loop for box
           deallocate(SITETYP, SITESIZ, ATOMID, SQNC, SITEFLAG, ITYP, ICLS, LID, LTY, XP, VECT, CFLAG, SEP, CS)

           !--- record Voronoi volumes
           if(m_OUTVORVOL) then
              write(*,fmt="(A,A)") ' MDPSCU Message: to calculate the cluster evelope...'
              call Record_Voronoi_ClusterEnvelope_Type(m_OUTPATH(1:len_trim(m_OUTPATH))//"VORONOI/", &
                                                    TSTAMP, SimBox(1), CtrlParam, m_PROJTYP)
           end if
          return
  end subroutine Record_Cfg_Track_Cluster
  !**********************************************************************************

  !**********************************************************************************
  subroutine BefOneTest_Cfg_Track_Cluster(ITest, SimBox0, SimBox, CtrlParam)
  !***  PURPOSE:  conduct before one test
  !    INPUT:     ITest,      the test ID
  !               SimBox0,    the simulation box, useless here
  !               SimBox,     the simulation box array, useless here
  !               CtrlParam,  the control parameter, useless here
  !   OUTPUT:
  !
    use RefVoronoiVacancy_GPU, only:hm_RefSimBox
    implicit none
    integer                     ::Itest
    type(SimMDBox), intent(in)  ::SimBox0
    type(SimMDBox), dimension(:)::SimBox
    type(SimMDCtrl)            ::CtrlParam
    !---
    integer::I, IFLAG, J
    character*256::GFILE, STRSWP
    character*32::NUMSTR

             if(.not. allocated(m_SQNC))    allocate(m_SQNC(hm_RefSimBox%NPRT*CtrlParam%MULTIBOX))
             if(.not. allocated(m_hFiles))  then
                 allocate(m_hFiles(CtrlParam%MULTIBOX))
                 m_hFiles = 0
             else
                 do I=1, size(m_hFiles)
                    if(m_hFiles(I) .gt. 0) close(m_hFiles(I))
                 end do
             end if

             do I=1, size(m_hFiles)
                call STRCATI(GFILE, m_TRACKPATH(1:len_trim(m_TRACKPATH)), "Defects_P", 0, 4)
                call STRCATI(GFILE, GFILE, "_", (Itest-1)*size(m_hFiles)+I, 4)
                call AvailableIOUnit(m_hFiles(I))
                if(CtrlParam%RESTART .le. 0) then
                   open(UNIT=m_hFiles(I), file = GFILE)
                   write(m_hFiles(I), fmt="(A, I8)")    '!--- THE  CLUSTER STATISTIC CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
                   write(m_hFiles(I), fmt="(A)")        '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
                   write(m_hFiles(I), fmt="(A)")        '!    '
                   write(m_hFiles(I), fmt="(A)")        '!--- Data interpretation:'
                   write(m_hFiles(I), fmt="(A)")        '!    COL #1: ID of recording'
                   write(m_hFiles(I), fmt="(A)")        '!    COL #2: time steps at the recording'
                   write(m_hFiles(I), fmt="(A)")        '!    COL #3: time (ps) at the recording'
                   write(m_hFiles(I), fmt="(A)")        '!    COL #4: number of clusters'
                   IFLAG = 4

                   IFLAG = IFLAG + 1
                   STRSWP = "!    COL #"
                   write(NUMSTR, *) IFLAG
                   NUMSTR = adjustl(NUMSTR)
                   STRSWP = STRSWP(1:len_trim(STRSWP))//NUMSTR(1:len_trim(NUMSTR))
                   IFLAG = IFLAG + m_MXCLUS-1
                   write(NUMSTR, *) IFLAG
                   NUMSTR = adjustl(NUMSTR)
                   STRSWP = STRSWP(1:len_trim(STRSWP))//" to #"//NUMSTR(1:len_trim(NUMSTR))//": cluster sizes for"
                   write(NUMSTR, *) 1
                   NUMSTR = adjustl(NUMSTR)
                   STRSWP = STRSWP(1:len_trim(STRSWP))//" C"//NUMSTR(1:len_trim(NUMSTR))//"-C"
                   write(NUMSTR, *) m_MXCLUS
                   NUMSTR = adjustl(NUMSTR)
                   STRSWP = STRSWP(1:len_trim(STRSWP))//NUMSTR(1:len_trim(NUMSTR))
                   write(m_hFiles(I), fmt="(A)")  STRSWP(1:len_trim(STRSWP))

                   do J=1, m_MXCLUS
                      IFLAG = IFLAG + 1
                      STRSWP = "!    COL #"
                      write(NUMSTR, *) IFLAG
                      NUMSTR = adjustl(NUMSTR)
                      STRSWP = STRSWP(1:len_trim(STRSWP))//NUMSTR(1:len_trim(NUMSTR))
                      IFLAG = IFLAG + m_MXCLUS-1
                      write(NUMSTR, *) IFLAG
                      NUMSTR = adjustl(NUMSTR)
                      STRSWP = STRSWP(1:len_trim(STRSWP))//" to #"//NUMSTR(1:len_trim(NUMSTR))//": seperation between "
                      write(NUMSTR, *) J
                      NUMSTR = adjustl(NUMSTR)
                      STRSWP = STRSWP(1:len_trim(STRSWP))//" C"//NUMSTR(1:len_trim(NUMSTR))// " and C"//NUMSTR(1:len_trim(NUMSTR))//"-"
                      write(NUMSTR, *) m_MXCLUS
                      NUMSTR = adjustl(NUMSTR)
                      STRSWP = STRSWP(1:len_trim(STRSWP))//NUMSTR(1:len_trim(NUMSTR))
                      write(m_hFiles(I), fmt="(A)")  STRSWP(1:len_trim(STRSWP))
                   end do
                   write(m_hFiles(I), fmt="(A)")        '!----    '
                   write(m_hFiles(I), fmt="(A)")
                else
                   open(UNIT=m_hFiles(I), file = GFILE,  POSITION='APPEND')
                end if
             end do

             return
  end subroutine BefOneTest_Cfg_Track_Cluster
  !**********************************************************************************

  !**********************************************************************************
  subroutine AftOneTest_Cfg_Track_Cluster(SimBox0, Stamp, SimBox, CtrlParam)
  !***  PURPOSE:  conduct after one test
  !    INPUT:
  !               SimBox0,    the simulation box, useless here
  !               Stamp,      the record stamp
  !               SimBox,     the simulation box array, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
    implicit none
    type(SimMDBox),      intent(in) ::SimBox0
    type(MDRecordStamp), intent(in) ::Stamp
    type(SimMDBox), dimension(:)    ::SimBox
    type(SimMDCtrl)                 ::CtrlParam
    integer::I
             if(allocated(m_SQNC)) deallocate(m_SQNC)
             if(allocated(m_hFiles) ) then
                do I=1, size(m_hFiles)
                   if(m_hFiles(I) .gt. 0) then
                      close(m_hFiles(I))
                   end if
                end do
                deallocate(m_hFiles)
             end if
             return
  end subroutine AftOneTest_Cfg_Track_Cluster
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_Cfg_Track_Cluster(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
  use MD_NeighborsList,        only:Clear_NeighboreList
  use RefVoronoiVacancy_GPU,   only:Clear_WorkingArray_DEV
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_NeighboreList(m_SiteNeighbor)
          call Clear_WorkingArray_DEV(SimBox, CtrlParam)
          return
  end subroutine Clear_Cfg_Track_Cluster
  !**********************************************************************************
 end module Cfg_Track_Cluster_GPU
