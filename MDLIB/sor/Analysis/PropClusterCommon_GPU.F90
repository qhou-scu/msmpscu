  module PropClusterCommon_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to common framework for cluster atoms according to their type and assinged
  !                  property. The type and assinged property are loaded from files created by MD simulation
  !                  or other analysis tool, for example, RefVoronoiVA. The atoms that are in neioghbor and
  !                  share the same type or property will be placeed in a cluster. The neighbor relation is
  !                  define by Delaunay tessellation. The clutsreing process may be performed on CPU. However,
  !                  the GPU version for Delanuay tessellation to be used, cuda fortran support is required.
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       VoronoiTessellationM_GPU.F90
  !                       MD_Globle_Variables_GPU.F90
  !
  !                  SEE ALSO____________________________________________________________________________
  !                       VoronoiTessellationM_GPU.F90
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2014-10 (Hou Qing)
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
  !-----------------------------------------------
    implicit none

      integer::m_processid=0
      !$$--- the id of  filename get from SimMDCtrl for I/O.
      !$$    Refering the definition of f_others in SimMDCtrl.
      character(len=12),parameter, private::mp_FTAGI="&AUXF_PCLSIN"
      character(len=13),parameter, private::mp_FTAGO="&AUXF_PCLSOUT"
      character(len=256)::m_INFILE =""                       ! filename of input control data
      character(len=256)::m_OUTFILE =""                      ! filename of output data
      character(len=256), private::m_CLUSTPATH =""           ! the path output clusteringing data
      character(len=256), private::m_VOLPATH =""             ! the path output cluster volum
      character(len=256), private::m_DTPATH =""              ! the path output DT data

      integer::m_OUTPUTDT           = 0                      ! flag indicteing if intermediate Delaunay vertice to be saved
      integer::m_OUTPUTVOL          = 0                      ! flag indicteing if output envelope volume

      integer, parameter, private       ::mp_Cluster_EndFlag = 0
      integer, parameter, private       ::mp_Box_EndFlag     = -1
      integer                           ::m_NCS              ! number of clusters
      integer, dimension(:), allocatable::m_SCS              ! arrary store size of clusters
      integer, dimension(:), allocatable::m_ASQCS            ! array  store the atom ID sequence of clusters
                                                             ! NOTE:  the atom ID is the ID of atoms of partitioned system
      integer, dimension(:), allocatable::m_STATCS           ! array  store the statu of the clusters, used only when
                                                             ! m_Method = mp_DELAUNAY

      integer            ::m_HISBIN    = 1                   ! number of bins for histogram of cluster size distribution
      integer, parameter ::mp_DELAUNAY = 1
      integer, parameter ::mp_BOND     = 2
      integer            ::m_Method    = mp_DELAUNAY          ! indcating the method used for clustering
                                                             ! 1, use Voronoi tessellation, default
                                                             ! 2, use bond lenght given by user
      !$$--- bond length only when m_Method = mp_BOND
      real(KINDDF)::m_BONDLEN      = -1


      !$$--- the type atoms to be clustered
      integer                               ::m_ATYP(mp_MXGROUP) = 0
      real(KINDDF),dimension(:), allocatable::m_EPOT
      real(KINDDF),dimension(:), allocatable::m_EKIN

      !$$--- the interface for Flag procedure.
      private::FLAGPROC
      abstract interface
        subroutine FLAGPROC(SimBox, CtrlParam, FLAG)
        !***  PURPOSE:   to create the FLAG witch identify the atoms
        !                that will be included in clustering
        !     INPUT:     SimBox:
        !                CtrlParam:
        !
        !     OUTPUT:    FLAG,
        !
          use MD_CONSTANTS
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)             ::SimBox
          type(SimMDCtrl)                          ::CtrlParam
          integer,        dimension(:), allocatable::FLAG

       end subroutine FLAGPROC
      end interface
      procedure(FLAGPROC), pointer, private::m_pFlagProc=>null()

      !--- calculation control parameters
      integer, private::m_INITED = 0

      private::Clear_WorkingArray

  contains

  !*****************************************************************************
  subroutine SetFlagProc(EXFLAGPROC)
  !***  DESCRIPTION: to set the user defrined FlagProc, which will be used to
  !                  create the flags of atoms that will be considered in clustering.
  !
  !     INPUT: PRERECORD,  the subroutine provided by user for pre-recording
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
        subroutine EXFLAGPROC(SimBox, CtrlParam, FLAG)
        !***  PURPOSE:   to create the FLAG which identify the atoms
        !                that will be included in clustering
        !     INPUT:     SimBox:
        !                CtrlParam:
        !
        !     OUTPUT:    FLAG,
        !
          use MD_CONSTANTS
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)             ::SimBox
          type(SimMDCtrl)                          ::CtrlParam
          integer,        dimension(:), allocatable::FLAG
       end subroutine EXFLAGPROC
  end interface
  !--- END INTERFACE --------------------------------
               m_pFlagProc  =>EXFLAGPROC
           return
   end subroutine SetFlagProc
  !*********************************************************************************

  !*****************************************************************************
  integer function IfInit_PropCluster() result(init)
    implicit none
          init = m_INITED

          return
  end function IfInit_PropCluster
  !*****************************************************************************

  !*****************************************************************************
  subroutine Allocate_WorkingArray()
  !***  DESCRIPTION: to allocate working memory
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT
  implicit none
  integer::IERR

        allocate(m_SCS(dm_NPRT), m_ASQCS(dm_NPRT), m_STATCS(dm_NPRT), STAT=IERR)
        if(IERR) then
           write(*,fmt="(A)") "MDPSCU Error: fail to allocate working memory in RECORD_PropCluster"
           write(*,fmt="(A)") "              Process to be stopped"
           stop
         end if
         m_INITED = IOR(m_INITED,1)

        return
  end subroutine Allocate_WorkingArray
  !*****************************************************************************

  !*****************************************************************************
  subroutine Clear_WorkingArray()
  !***  DESCRIPTION: to allocate working memory
  !
  implicit none
          if(allocated(m_SCS))    deallocate(m_SCS)
          if(allocated(m_ASQCS))  deallocate(m_ASQCS)
          if(allocated(m_STATCS)) deallocate(m_STATCS)

          if(allocated(m_EPOT))   deallocate(m_EPOT)
          if(allocated(m_EKIN))   deallocate(m_EKIN)

          m_NCS = 0
          m_INITED = 0
          return
  end subroutine Clear_WorkingArray
  !*****************************************************************************

  !****************************************************************************************
  subroutine Clear(SimBox,CtrlParam)
  !***  DESCRIPTION: to allocate mempry for diffusion calculations
  !
  !
  !***  The modules included ******************************************
  implicit none
     type(SimMDBox) ::SimBox
     type(SimMDCtrl)::CtrlParam

       call Clear_WorkingArray()
       return
  end subroutine Clear
  !*****************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, SimBox, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    SimBox,   the simulation box, with the member propColTitle
  !                          could be changed.
  !                CtrlParam, the control parameters, with the member NEEDDAMP
  !                           could be changed.
  !
  !
  !
  use MD_Globle_Variables,only:Load_ExAnalyCtl_SimMDCtrl
  implicit none
  !----   DUMMY Variables
   character*(*)::fname
   type(SimMDBox)::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
   integer::N, LINE, I, NN, IS
   character*256::STR
   character*32::STRNUMB(mp_mxGROUP), KEYWORD
   type(StatementList), pointer::StatList

  !--------
           call Load_ExAnalyCtl_SimMDCtrl(mp_FTAGI, Fname, SimBox, CtrlParam, &
                                             OutTag=mp_FTAGO, OutFile=m_OUTFILE, TheInput=StatList)

            do IS=1, Number_StatementList(StatList)
               call Get_StatementList(IS, StatList, STR, Line=LINE)
               call GetKeyWord("&", STR, KEYWORD)
               call UpCase(KEYWORD)

                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case("&HISBINS")
                            !$$*** To get number of bins to output size distribution
                            call Extract_Numb(STR,1,N,STRNUMB)
                            m_HISBIN = ISTR(STRNUMB(1))
                            if(m_HISBIN .le. 0) m_HISBIN = 1

                        case("&METHOD")
                            !$$*** To get the method that determines the connectivity of atoms
                            call Extract_Substr(STR,2,N,STRNUMB)
                            m_Method = 0
                            if(N.ge.1) then
                               if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "USEDT") then
                                  m_Method = mp_DELAUNAY
                                  if(N.ge.2 .and. STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "SAVE_DV") then
                                       m_OUTPUTDT = 1
                                  else
                                       m_OUTPUTDT = 0
                                  endif
                                else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "USEBD") then
                                   m_Method = mp_BOND
                                   call Extract_Numb(STR,1,N,STRNUMB)
                                   if(n .ge. 1) m_BONDLEN = DRSTR(STRNUMB(1))
                                   if(m_BONDLEN .le. 0) then
                                      write(*,fmt="(A,I4,A)")  ' MDPSCU Error: bond length use to determine the connectity '
                                      write(*,fmt="(A,I4,A)")  '               bond length should be larger than zero'
                                      write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                      write(*,fmt="(A)")       '               Process to be stopped'
                                      stop
                                   end if
                                end if
                           end if

                           if(m_Method .le. 0) then
                                write(*,fmt="(A,I4,A)")  ' MDPSCU Error: method to determine connectity is not given'
                                write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                write(*,fmt="(A, BZI6)") '        Usage: &METHOD  "USEDT", or &METHOD "USEDT" "SAVE_DV"'
                                write(*,fmt="(A, BZI6)") '               to use Delaunay tessellation to determine connectity of atoms'
                                write(*,fmt="(A, BZI6)") '               if having "SVAE_DV", Delaunay vertice will be saved'
                                write(*,fmt="(A, BZI6)") '        Or:    &METHOD  "USEBD" bondlen"'
                                write(*,fmt="(A, BZI6)") '               to use bond length to determine connectity of atoms'
                                write(*,fmt="(A, BZI6)") '               if bondlen is not given, force cutoff distance will be used as bond length'
                                write(*,fmt="(A)")       '               Process to be stopped'
                                stop
                           end if

                        case( "&SAVE_DV")
                             !*** get if Delaunay vertice to be save
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_OUTPUTDT = 0
                             else
                                m_OUTPUTDT = ISTR(STRNUMB(1))
                             end if

                        case( "&SAVE_VOL")
                             !*** get if envelpe volume to be save
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_OUTPUTVOL = 0
                             else
                                m_OUTPUTVOL = ISTR(STRNUMB(1))
                             end if

                         case( mp_FTAGO)
                              call Extract_Substr(STR,1,n,m_OUTFILE)

                         !*** to get the properties to be included in clustering
                         case("&PROP_TYPE")
                              call Extract_Numb(STR,SimBox%nGroup,N,STRNUMB)
                              do I=1, N
                                 NN = ISTR(STRNUMB(I))
                                 if(NN .gt.0) m_ATYP(NN) = 1
                              end do

                         case("&PROP_EPOT")
                              call Extract_Numb(STR,1,n,STRNUMB)
                              if(N.ge.1) then
                                 NN = 2*ISTR(STRNUMB(1))
                                 if(.not.allocated(m_EPOT))allocate(m_EPOT(NN))
                                 call Extract_Numb(STR,NN+1,N,STRNUMB)
                                 if(N.LE.NN+1) then
                                    write(*,fmt="(A,I4,A)")  ' MDPSCU Error: there are ',NN/2, ' pairs of EPOT range are required'
                                    write(*,fmt="(A,I4,A)")  '               but only ', (N-1)/2, ' pairs are availbel'
                                    write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                    write(*,fmt="(A, BZI6)") '        Usage: n, ( E1, E2 ), ( E3, E4 ),...,(En-1, En)'
                                    write(*,fmt="(A)")       '               Process to be stopped'
                                    stop
                                 end if
                                 do I=1, NN
                                    m_EPOT(I*2-1) = DRSTR(STRNUMB(I*2))
                                    m_EPOT(I*2)   = DRSTR(STRNUMB(I*2+1))
                                 end do
                              end if

                         case("&PROP_EKIN")
                              call Extract_Numb(STR,1,n,STRNUMB)
                              if(N.ge.1) then
                                 NN = 2*ISTR(STRNUMB(1))
                                 if(.not.allocated(m_EKIN))allocate(m_EKIN(NN))
                                 call Extract_Numb(STR,NN+1,N,STRNUMB)
                                 if(N.LE.NN+1) then
                                    write(*,fmt="(A,I4,A)")  ' MDPSCU Error: there are ',NN/2, ' pairs of EKIN range are required'
                                    write(*,fmt="(A,I4,A)")  '               but only ', (N-1)/2, ' pairs are availbel'
                                    write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                    write(*,fmt="(A, BZI6)") '        Usage: n, ( E1, E2 ), ( E3, E4 ),...,(En-1, En)'
                                    write(*,fmt="(A)")       '               Process to be stopped'
                                    stop
                                 end if
                                 do I=1, NN
                                    m_EKIN(I*2-1) = DRSTR(STRNUMB(I*2))
                                    m_EKIN(I*2)   = DRSTR(STRNUMB(I*2+1))
                                 end do
                              end if
                  end select
              end do

            !$$--- check input consistent
             if(m_BONDLEN .gt. 0) then
                m_BONDLEN = m_BONDLEN*SimBox%RR
             else
                m_BONDLEN = -minval(CtrlParam%RU)
             end if

             if(m_Method .ne. mp_DELAUNAY) then
                m_OUTPUTDT = 0
             end if
            return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in PropClusterCommon_14_GPU module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop

            return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine PrintControlParameters(hFile, SimBox, CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
  !

   implicit none
  !----   DUMMY Variables
   integer,        intent(in)::hFile
   type(SimMDBox), intent(in)::SimBox
   type(SimMDCtrl),intent(in)::CtrlParam

  !----   Local variables
   integer::I, J
   character*32::KEYWORD

        !$$--- print out the control parameters
         write(hFile,fmt="(A)")     " !************ Control parameters for PropClustering module **********"
         write(hFile,fmt="(A)")     " !     Properties of atoms included for clustering: "
         do I=1, NumberofData_DataPad(SimBox%proKWDList)
            call Tag_DataPad(I, SimBox%proKWDList, KEYWORD)
            select case(KEYWORD)
                   case ("&TYPECOL")
                         do J=1, size(m_ATYP)
                            if(m_ATYP(J) .gt. 0) then
                               write(hFile,fmt="(A,10(I2,1X))") " !     TYPE of atoms for type: ",J
                            end if
                         end do

                   case ("&EPOTCOL")
                         write(hFile,fmt="(A,10('( 'E12.4,',',E12.4' )',1X))") &
                                                          " !     EPOT of atoms in ranges: ", m_EPOT(1:size(m_EPOT))
                   case ("&EKINCOL")
                         write(hFile,fmt="(A,10('( 'E12.4,',',E12.4' )',1X))") &
                                                          " !     EKIN of atoms in ranges: ", m_EKIN(1:size(m_EKIN))
            end select
         end do
         write(hFile,fmt="(A)")        " !     "

         write(hFile,fmt="(A)")        " !    Connectity of atoms to be determined by method: "
         select case(m_METHOD)
                case(mp_DELAUNAY)
                     write(hFile,fmt="(A)")     " !      Delaunay tesselleation (default) with: "
                     if(m_OUTPUTDT .gt. 0) then
                        write(hFile,fmt="(A)")  " !      intermediate Delaunay vertices to be saved"
                     else
                        write(hFile,fmt="(A)")  " !      intermediate Delaunay vertices not to be saved"
                     end if
                case(mp_BOND)
                     if(m_BONDLEN .gt. 0) then
                        write(hFile,fmt="(A, 1X, F7.3, A)") " !      Bond length: ", m_BONDLEN/SimBox%RR, " (LU)"
                     else
                        write(hFile,fmt="(A, 1X, F7.3, A")  " !      Bond length is force cutoff distance", -m_BONDLEN/SimBox%RR, " (LU)"
                        m_BONDLEN = dabs(m_BONDLEN)
                     end if
         end select
         write(hFile,fmt="(A)")  " !     "

         return
   end subroutine PrintControlParameters
  !**********************************************************************************

  !**********************************************************************************
   subroutine Initialize_PropClusteringCommon(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the module by allocate memries
   use MD_Globle_Variables_GPU, only:Check_DEVICES
   use MD_Globle_Variables,     only:gm_hFILELOG, CreateDataFolder_Globle_Variables
   use MD_TYPEDEF_PrintList,    only:Add_PrintProcess
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, IFILE

         !$$--- to check device version is used
          call Check_DEVICES()

         !$$--- to clear the memory allocated before if there is
         if(m_INITED .gt. 0) then
           call Clear_WorkingArray()
         end if

        !$$--- to findout the I/O unit
         IFILE = 0
         do I=1, SIZE(CtrlParam%f_tag)
            if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                m_OUTFILE = CtrlParam%f_others(I)
            end if
            if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
               IFILE = I
            end if
         end do
         !$$--- check the input files
         if(IFILE .le. 0) then
            write(*,fmt="(' MDPSCU Error: the control file is missed in PropClustering module')")
            write(*,*) "                 add the keyword in SETUP file:", mp_FTAGI
            write(*,fmt="(' Process to be stopped')")
            stop
         end if

         m_INFILE = CtrlParam%f_others(IFILE)
         write(*,fmt="(A)") " !**** Loading control data from: "//m_INFILE(1:len_trim(m_INFILE))
         if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                         m_INFILE(1:len_trim(m_INFILE))
         call LoadControlParameters(m_INFILE, SimBox, CtrlParam)

         call GetPath(m_OUTFILE,  m_CLUSTPATH)
         call GetPath(m_OUTFILE,  m_VOLPATH)
         call GetPath(m_OUTFILE,  m_DTPATH)
         call GetFname(m_OUTFILE, m_OUTFILE)

         m_CLUSTPATH = m_CLUSTPATH(1:len_trim(m_CLUSTPATH))//"CLUSTER/"//m_OUTFILE(1:len_trim(m_OUTFILE))
         m_VOLPATH   = m_VOLPATH(1:len_trim(m_VOLPATH))//"CLUSTVOL/"//m_OUTFILE(1:len_trim(m_OUTFILE))
         m_DTPATH    = m_DTPATH(1:len_trim(m_DTPATH))//"VORONOI/"//m_OUTFILE(1:len_trim(m_OUTFILE))
         call CreateDataFolder_Globle_Variables(m_CLUSTPATH)
         if(m_OUTPUTVOL .gt. 0) &
            call CreateDataFolder_Globle_Variables(m_VOLPATH)
         if(m_OUTPUTDT .gt. 0) &
            call CreateDataFolder_Globle_Variables(m_DTPATH)

        !$$--- print out the control parameters
         call Add_PrintProcess(PrintControlParameters)
         return
   end subroutine Initialize_PropClusteringCommon
  !**********************************************************************************

  !*********************************************************************************
  subroutine Do_Delaunay_Clustering(SimBox, CtrlParam, NC, CS, SQNC, STATCS)
  !***  PURPOSE:  to performing the clustering with connectity of atoms
  !               determined by Delaunay tessellation
  !
  !     INPUT:   SimBox,    the simulation boxs
  !              CtrlParam, the control parameters for simulation
  !     OUTPUT:  NC,        the number of seperated clusters
  !              CS,        the size of the clusters
  !              SQNC,      the id sequence of atoms that form the clusters
  !              STATCS,    the flag indicating if the cluster is opend
  !
  !     NOTE:    the indice of atoms are the system that have been ordered by call to
  !              NeighboreList calculations
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList,  only:NEIGHBOR_LIST, Clear_NeighboreList
   use MD_NeighborsList_GPU
   use VoronoiTessellation_GPU, only: Cal_Delaunay_Vertice,GetMaxNF,mp_VSTAT_OPENED
   implicit none
   !---dummy vaiables
       type(SimMDBox),dimension(:), intent(in)::SimBox
       type(SimMDCtrl),             intent(in)::CtrlParam
       integer                                ::NC
       integer,       dimension(:), intent(out)::SQNC, CS
       integer,       dimension(:), intent(out)::STATCS
   !--- local variables
       integer, dimension(:),   allocatable::STAT
       integer, dimension(:),   allocatable::DNF            ! number of  Delanuay vertice of atoms
       integer, dimension(:,:), allocatable::DVTS           ! ID of atoms constructing the Delanuay vertice
       integer, dimension(:),   allocatable::DVTST          ! stat of Delanuay volume

       type(NEIGHBOR_LIST)::List
       integer::I, J, K, IP0, IP1, STEP, IA, IS, NN, IERR, OFLAG

           !$$--- allocate working
           allocate(STAT(dm_NPRT), DNF(dm_NPRT), DVTS(dm_NPRT,GetMaxNF()), DVTST(dm_NPRT), STAT=IERR )

          if(IERR) then
             write(*,fmt="(A)") "MDPSCU Error: fail to allocate working memory in Do_Delaunay_Clustering"
             write(*,fmt="(A)") "              Process to be stopped"
             stop
          end if

          !$$--- get the nieghbore list created by MD_NeighborsList2C_12_GPU module
          !$$    the neighbore list is require to create the atom ID of Delaunay vertice
          call Copyout_NeighboreList_DEV(List, ORDER=0)
          call Cal_Delaunay_Vertice(SimBox(1), CtrlParam, List, DNF, DVTS, hSTAT=DVTST, hFile=m_OUTPUTDT)
          call Clear_NeighboreList(List)

          !$$--- when a substrate is open, the number of DNF could be
          !$$    larger than m_MXNF, the permitted number of Delaunay faces,
          !$$    for example, in free cluster calculations. The algorithm
          !$$    can not handle this case.
          if(maxval(DNF) .gt. size(DVTS, dim=2)) then
             write(*,fmt="(A, I)")         ' MDPSCU Error: atoms have Delaunay neighbores more than permiited value'
             write(*,fmt="(A, I5, A, I5)") '              ', maxval(DNF), ' vs ', size(DVTS, dim=2)
             write(*,fmt="(A, I5, A, I5)") '               &METHOD = "USEBD" is recommanded.'
             write(*,fmt="(A)")            '               Process to be stopped'
             stop
          end if

          !$$--- to filter the atoms that will be included in clustering
          if( .not. associated(m_pFlagProc)) then
              write(*,fmt="(' MDPSCU Error: property identifying routine is not provided in PropClustering module')")
              write(*,fmt="('               Process to be stopped')")
              stop
          end if
          call m_pFlagProc(SimBox, CtrlParam, STAT)

          SQNC  = mp_Cluster_EndFlag
          IP0   = 0
          CS    = 0
          NC    = 0
          OFLAG = 0
          do I=1, dm_NPRT
             if(STAT(I) .le. 0) cycle

             IP0       = IP0 + 1
             SQNC(IP0) = I
             IS        = IP0
             NN        = 1
             if(DVTST(I) .eq. mp_VSTAT_OPENED) OFLAG = 1
             !$$--- atom I has been visited
             STAT(I) = 0
             do while(IP0 .le. dm_NPRT)     !Loop=1, dm_NPRT
                if(SQNC(IP0) .eq. mp_Cluster_EndFlag) exit   !The end of the cluaster
                K    = SQNC(IP0)
                STEP = 0
                do J=1, DNF(K)
                   IA = DVTS(K, J)
                   if(STAT(IA).le.0 ) then
                     cycle
                   end if
                   STEP            = STEP + 1
                   SQNC(IS+STEP)   = IA
                   STAT(IA)        = 0
                   if(DVTST(IA) .eq. mp_VSTAT_OPENED) OFLAG = 1
                end do
                IP0 = IP0 + 1
                IS  = IS  + STEP
                NN  = NN  + STEP
             end do
             IP0        = IS
             NC         = NC+1
             CS(NC)     = NN
             STATCS(NC) = OFLAG
             OFLAG      = 0
          end do
          if(allocated(STAT))  deallocate(STAT)
          if(allocated(DNF))   deallocate(DNF)
          if(allocated(DVTS))  deallocate(DVTS)
          if(allocated(DVTST)) deallocate(DVTST)

      return
  end subroutine Do_Delaunay_Clustering
  !****************************************************************************************

 !*********************************************************************************
  subroutine Do_Bond_Clustering(SimBox, CtrlParam, NC, CS, SQNC)
  !***  PURPOSE:  to performing the clustering with connectity of atoms
  !               determined by Delaunay tessellation
  !
  !     INPUT:   SimBox,    the simulation boxs
  !              CtrlParam, the control parameters for simulation
  !              IFLAG,     the numeric value that indentify the cluster of atoms belong to
  !     OUTPUT:  NC,        the number of seperated clusters
  !              CS,        the size of the clusters
  !              SQNC,      the id sequence of atoms that form the clusters
  !
  !     NOTE:    the indice of atoms are the system that have been ordered by call to
  !              NeighboreList calculations
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList, only:NEIGHBOR_LIST, Clear_NeighboreList
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       type(SimMDBox),dimension(:), intent(in) ::SimBox
       type(SimMDCtrl),             intent(in) ::CtrlParam
       integer                                 ::NC
       integer,       dimension(:), intent(out)::SQNC, CS
   !--- local variables
       integer, dimension(:), allocatable::STAT
       type(NEIGHBOR_LIST)::List
       integer::I, J, K, IP0, IP1, STEP, IA, IS, NN, IM, IERR, IFPDX, IFPDY, IFPDZ
       real(KINDDF)::XP0(3), XP(3), bond2, BS(3), HBS(3)

           !$$--- allocate working
           IFPDX = CtrlParam%IFPD(1)
           IFPDY = CtrlParam%IFPD(2)
           IFPDZ = CtrlParam%IFPD(3)
           BS    = SimBox(1)%ZL
           HBS   = BS*C_HALF
           bond2 = m_BONDLEN*m_BONDLEN
           allocate(STAT(dm_NPRT), STAT=IERR)
           if(IERR) then
             write(*,fmt="(A)") "MDPSCU Error: fail to allocate working memory in Do_Bond_Clustering"
             write(*,fmt="(A)") "              Process to be stopped"
             stop
           end if

          !$$--- get the nieghbore list created by MD_NeighborsList2C_12_GPU module
          call Copyout_NeighboreList_DEV(List, ORDER=0)

          !$$--- to filter the atoms that will be included in clustering
          if( .not. associated(m_pFlagProc)) then
              write(*,fmt="(' MDPSCU Error: property identifying routine is not provided in PropClustering module')")
              write(*,fmt="('               Process to be stopped')")
              stop
          end if
          call m_pFlagProc(SimBox, CtrlParam, STAT)
          SQNC = mp_Cluster_EndFlag
          IP0  = 0
          CS   = 0
          NC   = 0
          do I=1, dm_NPRT
             if(STAT(I) .le. 0) cycle
             IP0 = IP0 + 1
             SQNC(IP0) = I
             IS = IP0
             NN = 1
             !--- atom I has been visited
             STAT(I) = 0
             do while(IP0 .le. dm_NPRT) !Loop=1, dm_NPRT
                if(SQNC(IP0) .eq. mp_Cluster_EndFlag) exit
                K   = SQNC(IP0)
                XP0(1:3) = hm_XP(K,1:3)
                IM  = List%KVOIS(K)
                STEP = 0
                do J=1, IM
                   IA = List%INDI(K, J)
                   if(STAT(IA).le.0 ) then
                     cycle
                   end if

                   XP(1:3) = hm_XP(IA,1:3) - XP0(1:3)
                   if(IFPDX.GT.0 .AND.(DABS(XP(1)) .GT. HBS(1))) then
                      XP(1) = XP(1) - DSIGN(BS(1),XP(1))
                   end if

                   if(IFPDY.GT.0 .AND.(DABS(XP(2)) .GT. HBS(2))) then
                      XP(2) = XP(2) - DSIGN(BS(2),XP(2))
                   end if

                   if(IFPDZ.GT.0 .AND.(DABS(XP(3)) .GT. HBS(3))) then
                      XP(3) = XP(3) - DSIGN(BS(3),XP(3))
                   end if

                   if(sum(XP*XP) .gt. bond2) cycle
                   STEP          = STEP + 1
                   SQNC(IS+STEP) = IA
                   STAT(IA)      = 0
               end do
               IP0 = IP0 + 1
               IS  = IS+STEP
               NN  = NN+STEP
             end do
             IP0    = IS
             NC     = NC+1
             CS(NC) = NN
          end do

          if(allocated(STAT))  deallocate(STAT)
          call Clear_NeighboreList(List)

      return
  end subroutine Do_Bond_Clustering
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_PropCluster(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to do the propcluster.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_Globle_Variables_GPU
  use VoronoiTessellation_GPU, only: OUTPUT_Delaunay_Vertices_Header, mp_VSTAT_OPENED
  implicit none
       !--- dummy varaibles
       type(MDRecordStamp),          intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        character*256::GFILE
        integer::IERR, JOB, ICFG



           write(*,*) "MDPSCU Message: start clustering..."
           if(iand(m_INITED, 1) .eq. 0) then
              call Allocate_WorkingArray()
           end if

           JOB  = Stamp%ITest
           ICFG = Stamp%IRec(1)
           select case(m_METHOD)
                  case(mp_DELAUNAY)
                      !*** if intermeidate Delaunay vertice is reuired to output
                      !    we need output the header here
                      if(m_OUTPUTDT>0) then
                        !$$--- output for Delaunay Tesslation
                        call STRCATI(GFILE, m_DTPATH, "P", m_processid, 4)
                        call STRCATI(GFILE, GFILE, "_", JOB, 4)
                        call STRCATI(GFILE, GFILE, ".", ICFG, 4)
                        call AvailableIOUnit(m_OUTPUTDT)
                        open(UNIT=m_OUTPUTDT, file = GFILE, status='unknown')
                        write(*,*) "Delaunay vertices to be saved in "//GFILE(1:len_trim(GFILE))
                        !$$--- write out the hearder of the file
                        call OUTPUT_Delaunay_Vertices_Header(m_OUTPUTDT)
                      end if
                      call Do_Delaunay_Clustering(SimBox, CtrlParam, m_NCS, m_SCS, m_ASQCS, m_STATCS)
                      !$$---
                      if(m_OUTPUTDT>0) close(m_OUTPUTDT)
                  case(mp_BOND)
                     call Do_Bond_Clustering(SimBox, CtrlParam, m_NCS, m_SCS, m_ASQCS)

                  case default
                      write(*,fmt="(' MDPSCU Error: no method was selected for clustering in PropClustering module')")
                      write(*,*) "                  check the control file "//m_INFILE(1:len_trim(m_INFILE))
                      write(*,fmt="('               Process to be stopped')")
                      stop

           end select

          return
  end subroutine Do_PropCluster
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_PropCluster(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    SEE ALSO:
  !            MD_SimBoxArray_AppShell_xx_GPU.F90
  !            MD_SimBoxArray_ToolShell_xx_GPU.F90
  !
  use MD_Globle_Variables_GPU
  use VoronoiTessellation_GPU, only: mp_VSTAT_OPENED, record_voronoi_clusterenvelope_sqn
  implicit none
       !--- dummy varaibles
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       integer, dimension(:), allocatable::SD, SB
       character*256::GFILE
       integer::hFile, I, IC, IP, SMI, SMX, NBIN, SH, IS, NS, JOB,ICFG


           JOB   = Stamp%ITest
           ICFG  = Stamp%IRec(1)
           call Do_PropCluster(Stamp, SimBox, CtrlParam)

           !$$--- to output the clusters
           call STRCATI(GFILE, m_CLUSTPATH, "P", m_processid, 4)
           call STRCATI(GFILE, GFILE, "_", JOB, 4)
           call STRCATI(GFILE, GFILE, ".", ICFG, 4)

           call AvailableIOUnit(hFile)
           open(UNIT=hFile, file = GFILE, status='unknown')
           write(*,*) "Instant clustering results  to be saved in "//GFILE(1:len_trim(GFILE))

           !$$--- get size distribution
           SMI = minval(m_SCS, m_SCS.gt.0)
           SMX = maxval(m_SCS, m_SCS.gt.0)
           NBIN = (SMX - SMI + 1)/m_HISBIN + 1
           allocate(SD(NBIN), SB(NBIN+1))
           SD    = 0
           SB(1) = SMI
           do I=1, NBIN
              SB(I+1) = SB(I) + m_HISBIN
              SD(I) = count(m_SCS(1:m_NCS) .ge. SB(I) .and. m_SCS(1:m_NCS).lt.SB(I+1))
            end do

           write(hFile, fmt="(A, I8)")           '!--- THE CLUSTERING RESULTS CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
           write(hFile, fmt="(A)")               '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
           write(hFile, fmt="(A)")               '!    AUTHOR: HOU Qing'
           write(hFile, fmt="(A)")               '!    '

           write(hFile, fmt="(A, I8)")           '!--- Number of atoms in clusters:', sum(m_SCS(1:m_NCS))
           write(hFile, fmt="(A, I8)")           '!--- Number of connective clusters:', m_NCS
           write(hFile, fmt="(A, I8, A, I8)")    '!---     with min size: ', SMI, ', max size: ', SMX
           if(m_METHOD .eq. mp_DELAUNAY) then
              write(hFile, fmt="(A, I8, A, I8)") '!---     with number of opened clusters: ', &
                                                 count(m_STATCS(1:m_NCS) .eq. mp_VSTAT_OPENED)
              write(hFile, fmt="(A, I8, A, I8)") '!---     with number of free atoms: ',      &
                                                  count(m_STATCS(1:m_NCS) .eq. mp_VSTAT_OPENED .and. m_SCS(1:m_NCS).eq.1)
           end if

           write(hFile, fmt="(A, I4, 1x, A, I4)")'!--- Histogram of size distribution with bin width:', m_HISBIN
           write(hFile, fmt="(A, 4x, A9, 1x, A9, 1x A9)")   '!--- S-Type', 'size from', 'to', 'count'
           NS = 0
           do I=1, NBIN
              if(SD(I) .gt. 0) then
                 NS = NS + 1
                 write(hFile, fmt="(A, I3, 4x, I9, 4x, I9,1x,I9))") '!    ', NS, SB(I), SB(I+1), SD(I)
                 SD(NS) = SB(I)
              end if
           end do

           !$$--- write out the XYZ format
           write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
           write(hFile, fmt="(A, I12)")             "&TESTID   ", Stamp%ITest
           write(hFile, fmt="(A, I12,A,I6)")        "&BOXID    ", Stamp%IBox(1),", ", Stamp%IBox(2)
           write(hFile, fmt="(A, I12,A,I6)")        ""
           write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", sum(m_SCS(1:m_NCS))
           write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox(1)%ZL/SimBox(1)%RR
           write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox(1)%BOXLOW/SimBox(1)%RR
           write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox(1)%RR*CP_CM2A
           write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 1, 2, 3
           write(hFile, fmt="(A,1X,I4, A)")         "&ID0COL   ", 4, " !--- ID of atoms in  original box"
           write(hFile, fmt="(A,1X,I4, A)")         "&IDPCOL   ", 5, " !--- ID of atoms in  partitioned box"
           write(hFile, fmt="(A,1X,I4, A)")         "&CLUSCOL  ", 6, " !--- Cluster # of the atoms belong to"
           write(hFile, fmt="(A,1X,I4, A)")         "&TYPECOL  ", 7, " !--- Atomic type"
           write(hFile, fmt="(A,1X,I4, A)")         "&STYPECOL ", 8, " !--- Size of the cluster"

           if(m_METHOD .eq. mp_DELAUNAY) then
           write(hFile, fmt="(A,1X,I4,A))")         "&CLSSTATCOL ",  9, " !--- State of the cluster"
           write(hFile, fmt="(12A)")                "!---",     "   X(latt.)  ",&
                                                                "   Y(latt.)  ",&
                                                                "   Z(latt.)  ",&
                                                                "     ID0     ",&
                                                                "     IDP     ",&
                                                                "     CLUS    ",&
                                                                "    A-Type   ",&
                                                                "    S-Type   ",&
                                                                "    CLSSTAT  "
           else
           write(hFile, fmt="(12A)")                "!---",    "   X(latt.)  ", &
                                                                "   Y(latt.)  ",&
                                                                "   Z(latt.)  ",&
                                                                "     ID0     ",&
                                                                "     IDP     ",&
                                                                "     CLUS    ",&
                                                                "    A-Type   ",&
                                                                "    S-Type   "
           end if

           IP = 0
           do IC=1, m_NCS
              !$$--- get S-Type of CS
              do IS=1, NS-1
                 if(m_SCS(IC) .GE. SD(IS) .and. m_SCS(IC) .lt. SD(IS+1)) then
                    exit
                 end if
              end do

              do I=1, m_SCS(IC)
                 IP = IP + 1
                 if(m_METHOD .eq. mp_DELAUNAY) then
                    write(hFile, fmt="(2x, 3(1PE13.4), 5X, I8, 2X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8)")   &
                                                  hm_XP(m_ASQCS(IP),1:3)/SimBox(1)%RR,                             &
                                                  hm_GID(m_ASQCS(IP)), m_ASQCS(IP), IC, hm_ITYP(m_ASQCS(IP)), IS,  &
                                                  m_STATCS(IC)
                 else
                    write(hFile, fmt="(2x, 3(1PE13.4), 5X, I8, 2X, I8, 4X, I8, 4X, I8, 4X, I8, 4X, I8)")   &
                                                  hm_XP(m_ASQCS(IP),1:3)/SimBox(1)%RR,                     &
                                                  hm_GID(m_ASQCS(IP)), m_ASQCS(IP), IC, hm_ITYP(m_ASQCS(IP)), IS

                 end if
              end do
           end do
           close(hFile)

           if(m_OUTPUTVOL .gt. 0) then
              call Record_Voronoi_ClusterEnvelope_Sqn(m_VOLPATH, Stamp, SimBox(1), CtrlParam, m_ASQCS)
           end if

           deallocate(SD, SB)
          return
  end subroutine RECORD_PropCluster
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_PropClusterTool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !                  comapred to RECORD_PropCluster.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    SEE ALSO:
  !            MD_SimBoxArray_ToolShell_xx_GPU.F90
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables


            if(Stamp%ITime .LT. 0) return
            call RECORD_PropCluster(Stamp, SimBox, CtrlParam)
            return
  end subroutine RECORD_PropClusterTool
  !****************************************************************************************

  !****************************************************************************************
   subroutine SelectAtoms(SimBox, CtrlParam, FLAG)
  !***  DESCRIPTION: to determin which atoms to be used in clustering
  !                  this is a default routine of selecting atoms.    ,
  !                  one may provide a specific routine if specific rule
  !                  of selecting atoms is implemented, and then call
  !                  SetFlagProc in its initialization
  !                  routine.
  !
  !    INPUT:  SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !    OUTPUT: FLAG,      =0, or >0, maker of selected atoms
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT, hm_ITYP, hm_GID, gm_EPOT=>hm_EPOT, gm_EKIN=>hm_EKIN
   implicit none
   !----   DUMMY Variables
      type(SimMDBox), dimension(:)             ::SimBox
      type(SimMDCtrl)                          ::CtrlParam
      integer,        dimension(:), allocatable::FLAG
   !---- Local variables
      integer::I, GID, J, NP, NK, NPROP, IP, IB, NPRT, TFLAG

            !$$---- NOTE, the atoms ID have been grouped by cell in creating neighbore-list
            !$$           ref.to MD_Globle_Variables_GPU,
            !$$           and    MD_NeighborsList2C_12_GPU.F90
            FLAG = 0
            do I=1, dm_NPRT
               if(m_ATYP(hm_ITYP(I)) .gt. 0) then
                  FLAG(I) = 1
               end if
            end do

           !$$--- check the energy critical
            if(allocated(m_EPOT)) then
               NP = size(m_EPOT)/2
               do I=1, dm_NPRT
                  if(FLAG(I) .le. 0) cycle !already excluded
                  GID = hm_GID(I)
                  TFLAG = 0
                  do J=1, NP
                     if(gm_EPOT(GID) .le. m_EPOT(2*J-1) .and. gm_EPOT(GID) .le. m_EPOT(2*J)) then
                       TFLAG = 1
                        exit
                     end if
                  end do
                  FLAG(I) = TFLAG
               end do
            end if

            if(allocated(m_EKIN)) then
               NP = size(m_EKIN)/2
               do I=1, dm_NPRT
                  if(FLAG(I) .le. 0) cycle !already excluded
                  GID = hm_GID(I)
                  TFLAG = 0
                  do J=1, NP
                     if(gm_EPOT(GID) .ge. m_EKIN(2*J-1) .and. gm_EPOT(GID) .le. m_EKIN(2*J)) then
                        TFLAG = 1
                        exit
                     end if
                  end do
                  FLAG(I) = TFLAG
               end do
            end if

            !NPROP = count(len_trim(SimBox(1)%proColTitle) .gt. 0)
            !if(NPROP .gt. 0) then
             !$$--- select atom by other properties
             !$$--- NOTE, proTable is not grouped
            !  NPRT = SimBox(1)%NPRT
            !  do I=1, dm_NPRT
            !     if(FLAG(I) .le. 0) cycle !already excluded
            !     GID = hm_GID(I)
            !     IB  = (GID-1)/NPRT
            !     IP  = GID - IB*NPRT
            !     IB = IB  + 1
            !     TFLAG = 0
            !     do J=1, NPROP
            !        if(SimBox(IB)%proTable(IP,J) .ge. m_PROTHD(2*J-1) .and.  &
            !           SimBox(IB)%proTable(IP,J) .le. m_PROTHD(2*J) ) then
            !           TFLAG = TFLAG + 1
            !        end if
            !     end do
            !     if(TFLAG .ne. NPROP) FLAG(I) = 0
            !  end do
            !end if
            return
  end subroutine SelectAtoms
  !****************************************************************************************
  end module PropClusterCommon_GPU
