 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This is a MC simulation program used to caculate the visiting lattices for an defect
 !                  in HCP lattice.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                  MC_VisitLattice_HPC.exe -s "filename"
 !                  where:
 !                        filename  - the name of the file storing the coordinates of HCP lattices.
 !                                    The file could be created by, for example MakeBox tool, or
 !                                    a configuration file created by MDPSCU
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-06 (Hou Qing, Sichuan university)
 !
 !

 module MC_VisitLattice_HPC_module
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none

      !--  the specific parameters
      integer, parameter, private::mp_NPT_A1D=2
      integer, parameter, private::mp_NPT_A2D=4
      integer, parameter, private::mp_NPT_A  =mp_NPT_A1D + mp_NPT_A2D
      integer, parameter, private::mp_NPT_C=3*2
      integer, parameter, private::mp_NPT=mp_NPT_A + mp_NPT_C

      integer, parameter, private::mp_MXNEVENTS = mp_NPT
      integer, parameter, private::mp_EVENT_FOR      = 1
      integer, parameter, private::mp_EVENT_BACK     = 2
      integer, parameter, private::mp_EVENT_OFFLINE  = 3
      integer, parameter, private::mp_EVENT_C        = 4

      !--  the axsis in crystal
      integer, private:: m_AXSISX(3) = (/1, 0, 0/)
      integer, private:: m_AXSISY(3) = (/0, 1, 0/)
      integer, private:: m_AXSISZ(3) = (/0, 0, 1/)

      integer,      private::m_NEVENTS = 4                     !the actual number of events involved
      real(KINDDF), private::m_JUMPT = 1.D0                    !the transition time of events in ps
      real(KINDDF), private::m_JUMPF(mp_MXNEVENTS)= 1.D0       !the relative occurrence frequency of each kind of events
      real(KINDDF), private::m_cumJUMPR(0:mp_MXNEVENTS)= 0.D0  !the accumulative transition rate of events

      integer, private::m_NTEST = 100                          !the number of simulation box
      integer, private::m_NJUMP = 10000                        !the number of simulation jump
      integer, private::m_EXNBOX(3) = (/50,50,50/)             !the extension of the simulation box

      real(KINDDF), dimension(:,:), allocatable::m_XPL         !the positions of template lattices
      integer,      dimension(:,:), allocatable::m_INDL        !the neighbore lattice id of lattices

      !--- statistical quantities
      integer, private::m_cumEVENT(mp_MXNEVENTS)= 0            !the cumlalted number of events
      integer, private::m_cumVISIT = 0                         !the cumlatted number of visisted lattice
      character(len=256),  private::m_OUTFILE =""              !filename of output data
      real(KINDDF), dimension(:,:), private,allocatable::m_SMD !standard mean displacement

      integer,      private::mp_SaveTrack = 1
      integer,      private::m_SaveFlag = 0

      real(KINDDF), private::m_EPSLON = 1.D-3
      integer,      private::m_hFileSum = 0

  contains

  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to load control parameters and allocate memories needed
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_Globle_Variables,only:CreateDataFolder_Globle_Variables
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, ERR, hFile, LINE, N
   character*256::STR
   character*32::STRNUMB(32), KEYWORD


          CtrlParam%RU      = SimBox%RR*1.8
          CtrlParam%NB_RM   = CtrlParam%RU
          CtrlParam%IFPD    = 1
         !--- load control parameters
          if(len_trim(gm_ctrlFileName(1)) .le. 0) then
             write(*,fmt="(A)")          " MDPSCU Warning: no control file for running "//gm_ExeName(1:len_trim(gm_ExeName))//" is given."
             write(*,fmt="(A,A,A)")      "                add -c ctrfilename in command line"
             call ONWARNING(gm_OnWarning)
             write(*,fmt="(A,A,A)")      " Default control parameter to be used"
          else
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = gm_ctrlFileName(1), status='old', err=200)

            LINE = 0
            do while(.TRUE.)
               call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
               STR = adjustl(STR)
               call GetKeyword("&", STR, KEYWORD)
               call UpCase(KEYWORD)
               select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                      case( "&NUMTESTS")
                           !*** get the number of test
                             call Extract_Numb(STR,1, N, STRNUMB)
                             m_NTEST = ISTR(STRNUMB(1))

                      case( "&NUMEVENTS")
                           !*** get the number of test
                             call Extract_Numb(STR,1, N, STRNUMB)
                             m_NJUMP = ISTR(STRNUMB(1))

                      case( "&EVENTTIME")
                           !*** get the average time for a event occurring
                             call Extract_Numb(STR,1, N, STRNUMB)
                             m_JUMPT = DRSTR(STRNUMB(1))

                      case( "&EVENTFRQ")
                           !*** get the patition frequency of events
                             call Extract_Numb(STR,mp_MXNEVENTS, N, STRNUMB)
                             do I=1, N
                                m_JUMPF(I) = DRSTR(STRNUMB(I))
                             end do

                      case( "&BOXEXPAND")
                           !*** get the patition frequency of events
                             call Extract_Numb(STR,3, N, STRNUMB)
                             if(N.eq.1) then
                                 m_EXNBOX(I) = ISTR(STRNUMB(1))
                             else if(N.eq.2) then
                                  m_EXNBOX(1) = ISTR(STRNUMB(1))
                                  m_EXNBOX(2) = ISTR(STRNUMB(2))
                                  m_EXNBOX(3) = m_EXNBOX(2)
                             else if(N.eq.3) then
                                do I=1, N
                                   m_EXNBOX(I) = ISTR(STRNUMB(I))
                                end do
                             end if
                         case( "&OUTFILE")
                              call Extract_Substr(STR,1,N,m_OUTFILE)

                         case( "&SAVE_TRACK")
                              call Extract_Substr(STR,1,N,STRNUMB)
                              STRNUMB(1) = adjustl(STRNUMB(1))
                              call UpCase(STRNUMB(1))
                              if(STRNUMB(1)(1:len_trim("YES")).eq. "YES") &
                                 m_SaveFlag = IOR(m_SaveFlag,mp_SaveTrack)
               end select
            end do
            !--- assigne the transition rate
 100        close(hFile)
          end if

          m_NEVENTS = 4
          m_cumJUMPR(0) = 0.D0
          do I=1, m_NEVENTS
             m_cumJUMPR(I) = sum(m_JUMPF(1:I))/sum(m_JUMPF(1:m_NEVENTS))
          end do

          m_EXNBOX = int(m_EXNBOX/2)*2 + 1
          call PrintControlParameters(6)

           !--- allocate memory and calculate the Neighbor sites
            if(allocated(m_XPL)) deallocate(m_XPL)
            if(allocated(m_INDL)) deallocate(m_INDL)
            allocate(m_XPL(SimBox%NPRT,3), m_INDL(SimBox%NPRT,mp_NPT), STAT=ERR)
            if(ERR) then
               write(*,fmt='(A)') "MDPSCU Error: Fail to allocate memory for neighbor list in MC_VisitLattice_HPC_module"
               write(*,fmt='(A)') "              You may want use less box-cells"
               write(*,fmt='(A)') "              Process to be stoped"
               stop
           end if
          call CreateNeighborSites(SimBox, CtrlParam, m_INDL, m_XPL)
          print *, "Neighbore sites created"

          !--- allocate memory to store the MSD
          if(allocated(m_SMD)) deallocate(m_SMD)
          allocate(m_SMD(m_NTEST, 6))
          m_SMD = 0.D0

         !--- to open out put file
         if(len_trim(m_OUTFILE) .le. 0) then
            m_OUTFILE =  gm_ExeName(1:len_trim(gm_ExeName))
         end if
         call CreateDataFolder_Globle_Variables(m_OUTFILE)
         call AvailableIOUnit(m_hFileSum)
         open(UNIT=m_hFileSum, file = m_OUTFILE(1:len_trim(m_OUTFILE))//".sum", status='unknown')
         call PrintControlParameters(m_hFileSum)
         write(m_hFileSum, fmt="(A18,7(1x, A10), 18(A14))")&
                "Count list for BOX", "Forward", "Backward", "Off-line", "Off-plane", "Tot.Event","NUM.Visit", "Cum.Visit", &
                 "N.V./N.E.", "SMDxx", "SMDxy", "SMDxz", "SMDyy", "SMDyz", "SMDzz"
!                 "N.V./N.E.", "SMDxx", "DSMDxx", "SMDxy", "DSMDxy", "SMDxz", "DSMDxz", &
!                 "SMDyy", "DSMDyy", "SMDyz", "DSMDyz", "SMDzz", "DSMDzz"
         return

    200  write(*,fmt="(' MDPSCU Error: fail to load control file in MC_VisitLattice_HPC_module')")
         write(*,fmt="('               check the existence of file: ', A)") gm_ctrlFileName(1)(1:len_trim(gm_ctrlFileName(1)))
         write(*,fmt="(' Process to be stopped')")
         stop

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine PrintControlParameters(hFile)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile,
  !     OUTPUT:
  !
  implicit none

     !--- dummy vaiables
     integer::hFile
     !--- local varaibles

           !$$--- print out the control parameters
            write(hFile,FMT="(A)")             "!************ "//gm_ExeName(1:len_trim(gm_ExeName))//" to be performed **********"
            write(hFile,FMT="(A)")             "!    With the control parameters:"
            write(hFile,FMT="(A ,I5)")         "!    Number of tests: ", m_NTEST
            write(hFile,FMT="(A ,I5)")         "!    Number of events for each test: ", m_NJUMP
            write(hFile,FMT="(A ,3(I5,1x))")   "!    Expansion of the template box: ",  m_EXNBOX(1:3)
            write(hFile,FMT="(A ,F9.4)")       "!    Average time(ps) of event: ",      m_JUMPT
            write(hFile,FMT="(A ,99(F9.4,1x))")"!    Relative occurence frequence: ",   m_JUMPF(1:m_NEVENTS)
            write(hFile,FMT="(A ,99(F9.4,1x))")"!    "

  end subroutine PrintControlParameters
  !**********************************************************************************

  !**********************************************************************************
  subroutine CreateNeighborSites(SimBox, CtrlParam, IND, XP)
  !***  DESCRIPTION: to create the TetrahdralSites
  !
  !    INPUT:  SimBox, the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    OUTPUT: IND, the indices of the first 12 neighbores of lattices
  !            XP,  the position of the lattice in lattice units
  !
  use MD_NeighborsList
  implicit none
       type(SimMDBox),   intent(in)::SimBox
       type(SimMDCtrl),  intent(in)::CtrlParam
       integer,      dimension(:,:)::IND
       real(KINDDF), dimension(:,:)::XP

!-----  Local varaible
        type(NEIGHBOR_LIST)::List
        integer::I, IP, J, K, NPRT, NNA, IAP, IA(mp_NPT), NN
        real(KINDDF)::RR, DIST, SEP(3), XPP(3), BL(3), BU(3), BS(3), HBS(3), PTX(mp_NPT), PTY(mp_NPT), PTZ(mp_NPT), RC(mp_NPT)


            !$$--- create neighbor-list
            call Cal_NeighboreList2C(SimBox, CtrlParam, List)
            NPRT = SimBox%NPRT

            !$$--- get the secondary neighbors for each atoms
            IND = 0
            RR  = SimBox%RR
            BS  = SimBox%ZL/RR
            HBS = 0.5D0*BS
            BL  = SimBox%BOXLOW/RR
            BU  = SimBox%BOXUP/RR
            XP(1:NPRT, 1:3) = SimBox%XP(1:NPRT, 1:3)/RR

            !$$--- calculate the neighbor sites in basal-plane
            do I=1, NPRT

               NNA = List%KVOIS(I)
               NN  = 0
               RC  = 1.D32
               IA  = 0
               do J=1, NNA
                  IAP = List%INDI(I, J)
                  XPP(1:3) = XP(IAP,1:3) - XP(I,1:3)
                  if(dabs(XPP(1)) .gt. HBS(1)) &
                     XPP(1) = XPP(1) - DSIGN(BS(1),XPP(1))
                  if(dabs(XPP(2)) .gt. HBS(2)) &
                     XPP(2) = XPP(2) - DSIGN(BS(2),XPP(2))
                  if(dabs(XPP(3)) .gt. HBS(3)) &
                     XPP(3) = XPP(3) - DSIGN(BS(3),XPP(3))

                  if(dabs(XPP(3)) .le. m_EPSLON) then !--- for in basal sites
                     do  IP=1, mp_NPT_A
                         DIST = XPP(1)*XPP(1) + XPP(2)*XPP(2)
                         if(DIST .le. RC(IP)) then
                            IA(IP+1:mp_NPT_A) = IA(IP:mp_NPT_A -1)
                            RC(IP+1:mp_NPT_A) = RC(IP:mp_NPT_A -1)
                            IA(IP) = J
                            RC(IP) = DIST
                            exit
                         end if
                     end do
                  else !--- for off basal sites
                     do  IP=mp_NPT_A+1, mp_NPT
                         DIST = XPP(1)*XPP(1) + XPP(2)*XPP(2) + XPP(3)*XPP(3)
                         if(DIST .le. RC(IP)) then
                            IA(IP+1:mp_NPT) = IA(IP:mp_NPT -1)
                            RC(IP+1:mp_NPT) = RC(IP:mp_NPT -1)
                            IA(IP) = J
                            RC(IP) = DIST
                            exit
                         end if
                     end do

                  end if
               end do

               do IP=1, mp_NPT
                  IND(I,IP) =  List%INDI(I, IA(IP))
               end do
            end do

            call Clear_NeighboreList(List)
          return
  end subroutine CreateNeighborSites
  !**********************************************************************************

  !**********************************************************************************
  subroutine NPE_For_One_MC_Test(SimBox, CtrlParam, IB, IL0)
  !***  DESCRIPTION: to simulate the random walk of particle on hcp-lattice for a given
  !                  templeate lattice box. To obtaine the number of lattices visited
  !                  per event
  !
  !    INPUT:  SimBox,      the template box
  !            CtrlParam,   the control parameters for simulation, needed if neighbore list needed
  !                         to be upadates
  !            IB,          the ID of the test
  !            IL0,         optional, the initial lattice for starting random walk
  !                         if without given, start point will be randomly selected
  !                         from the template box
  !
  use RAND32_MODULE
  implicit none
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam
       integer,         intent(in)::IB
       integer,         optional  ::IL0
       !--- local variables
        character*256::GFILE
!-----  Local varaible
        integer::IJ, I, IS0, IS1, IS2, IS, IC1(3), IC2(3), ICXYZ1, ICXYZ2, NPRT0
        integer, dimension(:), allocatable::VISIT
        real(KINDDF)::RC, POS01(3), POS02(3), POS0(3), POS1(3), POS2(3), ORIG(3), VECT(3), V(3), R1, R2, COSDIR, RR, SMD(6), DSMD(6)
        real(KINDDF)::BL(3), BU(3), BH(3), BCS(3), HBCS(3)
        integer::IE, preEVENT, NEXTSITE, J, JJ, NCELL, ERR, hFile, NEWNV

        integer::dbgIC(3), dbgICXYX, dbgIS
        real(KINDDF)::dbgPos(3)

            !call AvailableIOUnit(hFile)
            !GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))
            !call STRCATI(GFILE,  GFILE, "_track.", IB, 4)
            !open(UNIT=hFile, file =GFILE(1:len_trim(GFILE)), status='unknown')
            !  write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            !  write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", NV
            !  write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", BU-BL
            !  write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", BL
            !  write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox%RR*CP_CM2A
            !  write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
            !  write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
            !  write(hFile, fmt="(A,1X,3(I4,1X))")      "&NVISITCOL ",5


           !--- determine the tempelate box size
            BCS    =  SimBox%ZL/SimBox%RR
            HBCS   =  0.5D0*BCS
            NCELL =  m_EXNBOX(1)*m_EXNBOX(2)*m_EXNBOX(3)
            BH    =  BCS*m_EXNBOX*C_HALF
            BL    = -BH
            BU    =  BH
            NPRT0 =  SimBox%NPRT
           !--- check the consistent
            if(allocated(m_XPL)) then
               if(size(m_XPL, dim=1) .ne. NPRT0) then
                  write(*,fmt='(A)') "MDPSCU Error: The size of template lattice box is not consistent"
                  write(*,fmt='(A, I, A, I)') "            ", NPRT0, "vs", size(m_XPL, dim=1)
                  write(*,fmt='(A)') "              Process to be stoped"
                  stop
               end if
            end if

           !--- allocate memory and calculate the Neighbor sites
           if(.not. allocated(m_XPL)) then
              allocate(m_XPL(NPRT0,3), m_INDL(NPRT0,mp_NPT), STAT=ERR)
              if(ERR) then
                write(*,fmt='(A)') "MDPSCU Error: Fail to allocate memory for neighbor list in MC_VisitLattice_HPC_module"
                write(*,fmt='(A)') "              You may want use less box-cells"
                write(*,fmt='(A)') "              Process to be stoped"
                stop
               end if
               call CreateNeighborSites(SimBox, CtrlParam, m_INDL, m_XPL)
               print *, "Neighbore sites created"
            end if

           !--- allocate memory for SMD
            if(.not. allocated(m_SMD)) then
               allocate(m_SMD(m_NTEST, 6))
               m_SMD = 0.D0
            end if

           !--- allocate memory for storing visited sites
            allocate(VISIT(NPRT0*NCELL),STAT=ERR)
            if(ERR .or. size(VISIT) .le. 0) then
              write(*,fmt='(A)') "MDPSCU Error: Fail to allocate memory for VISIT in MC_VisitLattice_HPC_module"
              write(*,fmt='(A)') "              You may want use less box-cells"
              write(*,fmt='(A)') "              Process to be stoped"
              stop
            end if

            !--- determine the initial site in template box
            if(present(IL0)) then
                IS0 = IL0
            else
                IS0 = DRAND32()*NPRT0
            end if

            !--- place the initial site closest to origin
             VISIT = 0
             NEWNV = 0
                 !--- determine the box-cell
                 POS1(1:3) = m_XPL(IS0,1:3)
                 IC1(1:3)  = int( (POS1(1:3) - BL(1:3))/BCS(1:3) - 0.000001)
                 ICXYZ1    = (IC1(3)*m_EXNBOX(2)+IC1(2))*m_EXNBOX(1) + IC1(1)
                 IS        = ICXYZ1*NPRT0 + IS0
!                VISIT(IS) = VISIT(IS) + 1
                 ORIG(1:3) = IC1(1:3)*BCS(1:3) + HBCS(1:3) + BL(1:3)
                 POS1(1:3) = POS1(1:3) + ORIG(1:3)
                 POS0(1:3) = POS1(1:3)

            !--- start the first moving
             R1 = DRAND32()
             do I=1, m_NEVENTS
                if(R1 .gt. m_cumJUMPR(I-1) .and. R1 .le. m_cumJUMPR(I)) then
                   IE = I
                   exit
                 end if
             end do

             select case(IE)
                    case( mp_EVENT_FOR,  mp_EVENT_BACK, mp_EVENT_OFFLINE)
                         I = DRAND32()*mp_NPT_A + 1.000001
                    case( mp_EVENT_C)
                        I = DRAND32()*mp_NPT_C + 1.000001 + mp_NPT_A
             end select
             IS1  = m_INDL(IS0, I)
             IS2  = IS1
             POS01(1:3) = m_XPL(IS0,1:3)
             POS02(1:3) = m_XPL(IS2,1:3)
             POS2(1:3)  = POS1(1:3) + (POS02(1:3) - POS01(1:3))
             IC2(1:3)  = int( (POS2(1:3) - BL(1:3))/BCS(1:3) - 0.000001)
             ICXYZ2    = (IC2(3)*m_EXNBOX(2)+IC2(2))*m_EXNBOX(1) + IC2(1)
             IS        = ICXYZ2*NPRT0 + IS2
             VISIT(IS) = VISIT(IS) + 1
             NEWNV     = NEWNV + 1
             m_cumEVENT(IE) = m_cumEVENT(IE) + 1
             preEVENT   = IE
            !--- start the next jumps
             do IJ = 1, m_NJUMP-1
                !--- the moving direction in previous step
                 VECT(1:3) = POS02(1:3) - POS01(1:3)
                 if(dabs(VECT(1)) .gt. HBCS(1)) VECT(1) = VECT(1) - DSIGN(BCS(1),VECT(1))
                 if(dabs(VECT(2)) .gt. HBCS(2)) VECT(2) = VECT(2) - DSIGN(BCS(2),VECT(2))
                 if(dabs(VECT(3)) .gt. HBCS(3)) VECT(3) = VECT(3) - DSIGN(BCS(3),VECT(3))
                 RR        = dsqrt(sum(VECT(1:3)*VECT(1:3)) )
                 VECT(1:3) = VECT(1:3)/RR

                 POS01(1:3) = POS02(1:3)
                 POS1(1:3)  = POS2(1:3)
               !--- deteremine the event type
                 R1 = DRAND32()
                 do I=1, m_NEVENTS
                    if(R1 .gt. m_cumJUMPR(I-1) .and. R1 .le. m_cumJUMPR(I)) then
                       IE = I
                       exit
                    end if
                 end do

               !--- if previous jump is in Z direction, direction lose
               !    we randomly select a nieghbor site in plane
                 !print *, "preEVENT", preEVENT, IS1, IS2
                 if(preEVENT .eq. mp_EVENT_C) then
                    J  =  DRAND32()*mp_NPT_A + 1.000001
                    IS2 = m_INDL(IS1, J)
                 else
                  !--- if the previous jump in in-plane
                    select case(IE)
                         case( mp_EVENT_FOR)
                              do J=1, mp_NPT_A
                                 V(1:3) = m_XPL(m_INDL(IS1, J), 1:3) - POS01(1:3)
                                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                                 COSDIR = sum(V(1:3)*VECT(1:3))/dsqrt(sum(V(1:3)*V(1:3)))
                                 if(COSDIR .gt. 0.99) then
                                    IS2 = m_INDL(IS1, J)
                                    exit
                                 end if

                              end do

                         case( mp_EVENT_BACK)
                           !--- determine the backward site
                              do J=1, mp_NPT_A
                                 V(1:3) = m_XPL(m_INDL(IS1, J), 1:3) - POS01(1:3)
                                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                                 COSDIR = sum(V(1:3)*VECT(1:3))/dsqrt(sum(V(1:3)*V(1:3)))
                                 if(COSDIR .lt. -0.99) then
                                    IS2 = m_INDL(IS1, J)
                                    exit
                                 end if
                              end do

                          case( mp_EVENT_OFFLINE)
                           !--- determine the offline site
                              JJ =  int(DRAND32()*mp_NPT_A2D) + 1
                              IS2 = 0
                              do J=1, mp_NPT_A
                                 V(1:3) = m_XPL(m_INDL(IS1, J), 1:3) - POS01(1:3)
                                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                                 COSDIR = sum(V(1:3)*VECT(1:3))/dsqrt(sum(V(1:3)*V(1:3)))
                                 if(dabs(COSDIR) .gt. 0.99 ) cycle

                                 IS2 = IS2 + 1
                                 if(IS2 .eq. JJ) then
                                    IS2 = m_INDL(IS1, J)
                                    exit
                                 end if
                              end do

                         case( mp_EVENT_C)
                           !--- determine the C jump site
                             J = int(DRAND32()*mp_NPT_C) + 1 + mp_NPT_A
                             IS2 = m_INDL(IS1, J)
                             !print *, IS1, J, IND(IS1, J)
                             !do J=1, mp_NPT
                             !   write(*,fmt="(I,3(1x,1pE13.6))") J, XP(IND(IS1, J), 1:3)-XP(IS1,1:3)
                             !end do
                             !pause
                    end select
                 end if
                 m_cumEVENT(IE) = m_cumEVENT(IE) + 1
                 preEVENT   = IE
                !--- the displacement
                 POS02(1:3) =  m_XPL(IS2, 1:3)
                 V(1:3)     =  POS02(1:3) - POS01(1:3)
                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                 POS2(1:3)  =  POS1(1:3) + V(1:3)

                !--- the global index of the site
                 if(POS2(1) .lt. BL(1) .or. POS2(1) .ge. BU(1) .or.  &
                    POS2(2) .lt. BL(2) .or. POS2(2) .ge. BU(2) .or.  &
                    POS2(3) .lt. BL(3) .or. POS2(3) .ge. BU(3)  ) then
                     write(*,fmt="(A, I5, A, 1x, I, A)") &
                              " MDPSCU Warning: the particle move out of box in BOX",IB, "after ",IJ, "jumps"
                     write(*,fmt="(A, 1PE14.5, 1x, 1PE12.5, 1x, A, 1x, 1PE12.5)") &
                              "                 Box range X: ",BL(1), BU(1), " vs X ", POS2(1)
                     write(*,fmt="(A, 1PE14.5, 1x, 1PE12.5, 1x, A, 1x, 1PE12.5)") &
                              "                 Box range Y: ",BL(2), BU(2), " vs Y ", POS2(2)
                     write(*,fmt="(A, 1PE14.5, 1x, 1PE12.5, 1x, A, 1x, 1PE12.5)") &
                              "                 Box range Z: ",BL(3), BU(3), " vs Z ", POS2(3)
                       call ONWARNING(gm_OnWarning)
                       exit
                 end if
                 IC2(1:3)   = int( (POS2(1:3) - BL(1:3))/BCS(1:3) - 0.000001)
                 ICXYZ2     = (IC2(3)*m_EXNBOX(2)+IC2(2))*m_EXNBOX(1) + IC2(1)
                 IS     = ICXYZ2*NPRT0 + IS2
                 if(VISIT(IS) .le. 0) NEWNV = NEWNV + 1
                 VISIT(IS)  = VISIT(IS) + 1
                 IS1 = IS2
!--- for debug
!                 write(hFile, fmt="(I5, 1x, 3(1PE13.5, 1x),I7)") 1, POS2(1:3), VISIT(IS)
!                 dbgICXYX = (IS-1)/NPRT0
!                 dbgIS     = IS - dbgICXYX*NPRT0
!                 dbgIC(3) = dbgICXYX/(m_EXNBOX(1)*m_EXNBOX(2))
!                 dbgIC(2) = (dbgICXYX- dbgIC(3)*m_EXNBOX(1)*m_EXNBOX(2) )/m_EXNBOX(1)
!                 dbgIC(1) = dbgICXYX- (dbgIC(3)*m_EXNBOX(2)+dbgIC(2))*m_EXNBOX(1)
!                 ORIG(1:3) = dbgIC(1:3)*BCS(1:3) + HBCS(1:3) + BL(1:3)
!                 dbgPOS(1:3)  = m_XPL(dbgIS,1:3) + ORIG(1:3)
!                 if(any(dabs(dbgPOS - POS2) .gt. 0.001)) then
!                    print *, "dbg    ",  IB,IJ, IS, dbgICXYX, dbgIS, IS2
!                    print *, "IC2    ",  IC2
!                    print *, "dbgIC  ",  dbgIC
!                    print *, "POS2   ", POS2(1:3)
!                    print *, "dbgPOS ", dbgPOS(1:3)
!                    print *, "BL     ", BL
!                    print *, "HBCS     ", HBCS
!                    print *, "m_XPL     ", m_XPL(dbgIS,1:3)
!                    pause
!                 end if
!                 write(hFile, fmt="(I5, 1x, 3(1PE13.5, 1x),I7)") 1, dbgPOS(1:3), VISIT(IS)

             end do !--- end loop for total events

            !--- extract statistical quantities
             m_cumVISIT = m_cumVISIT + NEWNV
            !--- accumulate the SMD
             m_SMD(IB, 1) = (POS2(1) - POS0(1))*(POS2(1) - POS0(1))
             m_SMD(IB, 2) = (POS2(1) - POS0(1))*(POS2(2) - POS0(2))
             m_SMD(IB, 3) = (POS2(1) - POS0(1))*(POS2(3) - POS0(3))
             m_SMD(IB, 4) = (POS2(2) - POS0(2))*(POS2(2) - POS0(2))
             m_SMD(IB, 5) = (POS2(2) - POS0(2))*(POS2(3) - POS0(3))
             m_SMD(IB, 6) = (POS2(3) - POS0(3))*(POS2(3) - POS0(3))
             SMD(1) = sum(m_SMD(1:IB, 1))/dble(IB)
             SMD(2) = sum(m_SMD(1:IB, 2))/dble(IB)
             SMD(3) = sum(m_SMD(1:IB, 3))/dble(IB)
             SMD(4) = sum(m_SMD(1:IB, 4))/dble(IB)
             SMD(5) = sum(m_SMD(1:IB, 5))/dble(IB)
             SMD(6) = sum(m_SMD(1:IB, 6))/dble(IB)
             DSMD(1)= sum( (m_SMD(1:IB, 1) - SMD(1))*(m_SMD(1:IB, 1) - SMD(1)))/dble(IB)
             DSMD(2)= sum( (m_SMD(1:IB, 2) - SMD(2))*(m_SMD(1:IB, 2) - SMD(2)))/dble(IB)
             DSMD(3)= sum( (m_SMD(1:IB, 3) - SMD(3))*(m_SMD(1:IB, 3) - SMD(3)))/dble(IB)
             DSMD(4)= sum( (m_SMD(1:IB, 4) - SMD(4))*(m_SMD(1:IB, 4) - SMD(4)))/dble(IB)
             DSMD(5)= sum( (m_SMD(1:IB, 5) - SMD(5))*(m_SMD(1:IB, 5) - SMD(5)))/dble(IB)
             DSMD(6)= sum( (m_SMD(1:IB, 6) - SMD(6))*(m_SMD(1:IB, 6) - SMD(6)))/dble(IB)
             DSMD   = dsqrt(DSMD)

             write(m_hFileSum, fmt="(12x,I6,7(3x,I8),3x,18(1PE13.4,1x))") IB, m_cumEVENT(1:m_NEVENTS), sum(m_cumEVENT(1:m_NEVENTS)),NEWNV, &
                                                    m_cumVISIT, dble(m_cumVISIT)/sum(m_cumEVENT(1:m_NEVENTS)), (SMD(J), J=1,6)


            !--- to extract the trajectory of migration if required
            if( IAND(m_SaveFlag,mp_SaveTrack) .eq. mp_SaveTrack) then
               call AvailableIOUnit(hFile)
               GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))
               call STRCATI(GFILE,  GFILE, "_track.", IB, 4)

               open(UNIT=hFile, file =GFILE(1:len_trim(GFILE)), status='unknown')
               write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
               write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", NEWNV
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", BU-BL
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", BL
               write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox%RR*CP_CM2A
               write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
               write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
               write(hFile, fmt="(A,1X,3(I4,1X))")      "&NVISITCOL ",5
               do I=1, NPRT0*NCELL
                  if(VISIT(I) .le. 0) cycle
                  !--- determine the box-cell
                  ICXYZ1 = (I-1)/NPRT0
                  IS     = I - ICXYZ1*NPRT0
                  IC1(3) = ICXYZ1/(m_EXNBOX(1)*m_EXNBOX(2))
                  IC1(2) = (ICXYZ1- IC1(3)*m_EXNBOX(1)*m_EXNBOX(2) )/m_EXNBOX(1)
                  IC1(1) = ICXYZ1- (IC1(3)*m_EXNBOX(2)+IC1(2))*m_EXNBOX(1)
                  ORIG(1:3) = IC1(1:3)*BCS(1:3) + HBCS(1:3) + BL(1:3)
                  POS1(1:3)  = m_XPL(IS,1:3) + ORIG(1:3)
                  write(hFile, fmt="(I5, 1x, 3(1PE13.5, 1x),I7)") 1, POS1(1:3), VISIT(I)
               end do
               close(hFile)
            end if

            deallocate(VISIT)
          return
  end subroutine NPE_For_One_MC_Test
  !**********************************************************************************

  !**********************************************************************************
  subroutine SEGLEN_For_One_MC_Test(SimBox, CtrlParam, IL0)
  !***  DESCRIPTION: to simulate the random walk of particle on hcp-lattice for a given
  !                  templeate lattice box. To obtaine the distribution of  number of lattices
  !                  on the line-segments of particle visiting
  !
  !    INPUT:  SimBox,      the template box
  !            CtrlParam,   the control parameters for simulation, needed if neighbore list needed
  !                         to be upadates
  !            IB,          the ID of the test
  !            IL0,         optional, the initial lattice for starting random walk
  !                         if without given, start point will be randomly selected
  !                         from the template box
  !
  use RAND32_MODULE
  implicit none
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam
       integer,         optional  ::IL0
       !--- local variables
        character*256::GFILE
!-----  Local varaible
        integer::IJ, I, J, IS0, IS1, IS2, IS, IC1(3), IC2(3), ICXYZ1, ICXYZ2, NPRT0
        integer, dimension(:),   allocatable::VISIT, HISNV, HISNE, HISNVPE
        integer, dimension(:,:), allocatable::HISMAP
        real(KINDDF)::RC, POS01(3), POS02(3), POS0(3), POS1(3), POS2(3), ORIG(3), VECT(3), V(3), R1, R2, COSDIR, RR
        real(KINDDF)::BL(3), BU(3), BH(3), BCS(3), HBCS(3), NVMI, NVMX, NEMI, NEMX, PNVMI, PNVMX, PNV, DPNV
        integer::IE,  preEVENT, NCELL, NEWNV, NEWEV, HNVMI, HNEMI, FLAG, IBIN, NREC
        integer::ERR, hFile1, hFile2, hFile3, hFile4, hFileSwap

            call AvailableIOUnit(hFileSwap)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//"_segments.dat"
            open(UNIT=hFileSwap, file =GFILE(1:len_trim(GFILE)), status='unknown')
            call AvailableIOUnit(hFile1)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//"_His_NV.dat"
            open(UNIT=hFile1, file =GFILE(1:len_trim(GFILE)), status='unknown')
            call AvailableIOUnit(hFile2)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//"_His_NE.dat"
            open(UNIT=hFile2, file =GFILE(1:len_trim(GFILE)), status='unknown')
            call AvailableIOUnit(hFile3)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//"_His_PNV.dat"
            open(UNIT=hFile3, file =GFILE(1:len_trim(GFILE)), status='unknown')
            call AvailableIOUnit(hFile4)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//"_His_MAP.dat"
            open(UNIT=hFile4, file =GFILE(1:len_trim(GFILE)), status='unknown')
            !call PrintControlParameters(hFileSwap)
            call PrintControlParameters(hFile1)
            call PrintControlParameters(hFile2)
            call PrintControlParameters(hFile3)
            call PrintControlParameters(hFile4)

           !--- determine the tempelate box size
            BCS    =  SimBox%ZL/SimBox%RR
            HBCS   =  0.5D0*BCS
            NCELL =  m_EXNBOX(1)*m_EXNBOX(2)*m_EXNBOX(3)
            BH    =  BCS*m_EXNBOX*C_HALF
            BL    = -BH
            BU    =  BH
            NPRT0 =  SimBox%NPRT
           !--- check the consistent
            if(allocated(m_XPL)) then
               if(size(m_XPL, dim=1) .ne. NPRT0) then
                  write(*,fmt='(A)') "MDPSCU Error: The size of template lattice box is not consistent"
                  write(*,fmt='(A, I, A, I)') "            ", NPRT0, "vs", size(m_XPL, dim=1)
                  write(*,fmt='(A)') "              Process to be stoped"
                  stop
               end if
            end if

           !--- allocate memory and calculate the Neighbor sites
           if(.not. allocated(m_XPL)) then
              allocate(m_XPL(NPRT0,3), m_INDL(NPRT0,mp_NPT), STAT=ERR)
              if(ERR) then
                write(*,fmt='(A)') "MDPSCU Error: Fail to allocate memory for neighbor list in MC_VisitLattice_HPC_module"
                write(*,fmt='(A)') "              You may want use less box-cells"
                write(*,fmt='(A)') "              Process to be stoped"
                stop
               end if
               call CreateNeighborSites(SimBox, CtrlParam, m_INDL, m_XPL)
               print *, "Neighbore sites created"
            end if

           !--- allocate memory for storing visited sites
            allocate(VISIT(NPRT0*NCELL),STAT=ERR)
            if(ERR .or. size(VISIT) .le. 0) then
              write(*,fmt='(A)') "MDPSCU Error: Fail to allocate memory for VISIT in MC_VisitLattice_HPC_module"
              write(*,fmt='(A)') "              You may want use less box-cells"
              write(*,fmt='(A)') "              Process to be stoped"
              stop
            end if

            !--- determine the initial site in template box
            if(present(IL0)) then
                IS0 = IL0
            else
                IS0 = DRAND32()*NPRT0
            end if

            !--- place the initial site closest to origin
            NVMI = 1.D64
            NVMX = 0
            NEMI = 1.D64
            NEMX = 0
            NREC = 0
            do IJ=1, m_NTEST
               VISIT  = 0
               NEWNV  = 0
               NEWEV  = 0
               !--- determine the box-cell
               POS1(1:3) = m_XPL(IS0,1:3)
               IC1(1:3)  = int( (POS1(1:3) - BL(1:3))/BCS(1:3) - 0.000001)
               ICXYZ1    = (IC1(3)*m_EXNBOX(2)+IC1(2))*m_EXNBOX(1) + IC1(1)
               IS        = ICXYZ1*NPRT0 + IS0
               ORIG(1:3) = IC1(1:3)*BCS(1:3) + HBCS(1:3) + BL(1:3)
               POS1(1:3) = POS1(1:3) + ORIG(1:3)
               POS0(1:3) = POS1(1:3)
               NEWNV     = NEWNV + 1
               VISIT(IS) = VISIT(IS) + 1

               !--- start the first moving
               R1 = DRAND32()
               do I=1, m_NEVENTS
                  if(R1 .gt. m_cumJUMPR(I-1) .and. R1 .le. m_cumJUMPR(I)) then
                     IE = I
                     exit
                  end if
               end do

               select case(IE)
                      case( mp_EVENT_FOR,  mp_EVENT_BACK, mp_EVENT_OFFLINE)
                         I = DRAND32()*mp_NPT_A + 1.000001
                      case( mp_EVENT_C)
                         exit
                         I = DRAND32()*mp_NPT_C + 1.000001 + mp_NPT_A
               end select
               IS1  = m_INDL(IS0, I)
               IS2  = IS1
               POS01(1:3) = m_XPL(IS0,1:3)
               POS02(1:3) = m_XPL(IS2,1:3)
               POS2(1:3)  = POS1(1:3) + (POS02(1:3) - POS01(1:3))
               IC2(1:3)  = int( (POS2(1:3) - BL(1:3))/BCS(1:3) - 0.000001)
               ICXYZ2    = (IC2(3)*m_EXNBOX(2)+IC2(2))*m_EXNBOX(1) + IC2(1)
               IS        = ICXYZ2*NPRT0 + IS2
               VISIT(IS) = VISIT(IS) + 1
               NEWNV     = NEWNV + 1
               NEWEV     = NEWEV + 1
               preEVENT  = IE
               FLAG      = 0
               !--- start the next jumps
               do while(.true.)
                  !--- the moving direction in previous step
                  VECT(1:3) = POS02(1:3) - POS01(1:3)
                  if(dabs(VECT(1)) .gt. HBCS(1)) VECT(1) = VECT(1) - DSIGN(BCS(1),VECT(1))
                  if(dabs(VECT(2)) .gt. HBCS(2)) VECT(2) = VECT(2) - DSIGN(BCS(2),VECT(2))
                  if(dabs(VECT(3)) .gt. HBCS(3)) VECT(3) = VECT(3) - DSIGN(BCS(3),VECT(3))
                  RR        = dsqrt(sum(VECT(1:3)*VECT(1:3)) )
                  VECT(1:3) = VECT(1:3)/RR

                  POS01(1:3) = POS02(1:3)
                  POS1(1:3)  = POS2(1:3)
                  !--- deteremine the event type
                  R1 = DRAND32()
                  do I=1, m_NEVENTS
                     if(R1 .gt. m_cumJUMPR(I-1) .and. R1 .le. m_cumJUMPR(I)) then
                        IE = I
                        exit
                     end if
                  end do

                 !--- if the previous jump in in-plane
                 select case(IE)
                        case( mp_EVENT_FOR)
                              do J=1, mp_NPT_A
                                 V(1:3) = m_XPL(m_INDL(IS1, J), 1:3) - POS01(1:3)
                                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                                 COSDIR = sum(V(1:3)*VECT(1:3))/dsqrt(sum(V(1:3)*V(1:3)))
                                 if(COSDIR .gt. 0.99) then
                                    IS2 = m_INDL(IS1, J)
                                    exit
                                 end if

                              end do

                        case( mp_EVENT_BACK)
                           !--- determine the backward site
                              do J=1, mp_NPT_A
                                 V(1:3) = m_XPL(m_INDL(IS1, J), 1:3) - POS01(1:3)
                                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                                 COSDIR = sum(V(1:3)*VECT(1:3))/dsqrt(sum(V(1:3)*V(1:3)))
                                 if(COSDIR .lt. -0.99) then
                                    IS2 = m_INDL(IS1, J)
                                    exit
                                 end if
                              end do

                        case( mp_EVENT_OFFLINE)
                              exit

                        case( mp_EVENT_C)
                             exit

                            !--- determine the C jump site
                             J = int(DRAND32()*mp_NPT_C) + 1 + mp_NPT_A
                             IS2 = m_INDL(IS1, J)
                 end select

                !--- the displacement
                 POS02(1:3) =  m_XPL(IS2, 1:3)
                 V(1:3)     =  POS02(1:3) - POS01(1:3)
                 if(dabs(V(1)) .gt. HBCS(1)) V(1) = V(1) - DSIGN(BCS(1),V(1))
                 if(dabs(V(2)) .gt. HBCS(2)) V(2) = V(2) - DSIGN(BCS(2),V(2))
                 if(dabs(V(3)) .gt. HBCS(3)) V(3) = V(3) - DSIGN(BCS(3),V(3))
                 POS2(1:3)  =  POS1(1:3) + V(1:3)

                !--- the global index of the site
                 if(POS2(1) .lt. BL(1) .or. POS2(1) .ge. BU(1) .or.  &
                    POS2(2) .lt. BL(2) .or. POS2(2) .ge. BU(2) .or.  &
                    POS2(3) .lt. BL(3) .or. POS2(3) .ge. BU(3)  ) then
                     write(*,fmt="(A, I5, A, 1x, I, A)") &
                              " MDPSCU Warning: the particle move out of box in BOX",IJ, "after ",NEWEV, "jumps"
                     write(*,fmt="(A, 1PE14.5, 1x, 1PE12.5, 1x, A, 1x, 1PE12.5)") &
                              "                 Box range X: ",BL(1), BU(1), " vs X ", POS2(1)
                     write(*,fmt="(A, 1PE14.5, 1x, 1PE12.5, 1x, A, 1x, 1PE12.5)") &
                              "                 Box range Y: ",BL(2), BU(2), " vs Y ", POS2(2)
                     write(*,fmt="(A, 1PE14.5, 1x, 1PE12.5, 1x, A, 1x, 1PE12.5)") &
                              "                 Box range Z: ",BL(3), BU(3), " vs Z ", POS2(3)
                       call ONWARNING(gm_OnWarning)
                       FLAG = 1
                       exit
                 end if
                 IC2(1:3)   = int( (POS2(1:3) - BL(1:3))/BCS(1:3) - 0.000001)
                 ICXYZ2     = (IC2(3)*m_EXNBOX(2)+IC2(2))*m_EXNBOX(1) + IC2(1)
                 IS     = ICXYZ2*NPRT0 + IS2
                 if(VISIT(IS) .le. 0) NEWNV = NEWNV + 1
                 VISIT(IS)  = VISIT(IS) + 1
                 IS1 = IS2

                 NEWEV = NEWEV + 1
                 preEVENT   = IE
               end do !--- end loop for total events

               !--- save the length to temp file
               if(FLAG .eq. 0) then
                  if(NVMI .gt. NEWNV) NVMI = NEWNV
                  if(NVMX .lt. NEWNV) NVMX = NEWNV
                  if(NEMI .gt. NEWEV) NEMI = NEWEV
                  if(NEMX .lt. NEWEV) NEMX = NEWEV
                  NREC = NREC + 1
                  write(hFileSwap,fmt="(100I)")  NEWNV, NEWEV
                  print *, "TEST", IJ, NEWNV, NEWEV
               end if
            end do ! end do for test
            deallocate(VISIT)

            allocate(HISNV(int(NVMI+0.0001)-1:int(NVMX+0.0001)+1),HISNE(int(NEMI+0.0001)-1:int(NEMX+0.0001)+1),STAT=ERR)
            allocate(HISMAP(int(NVMI+0.0001)-1:int(NVMX+0.0001)+1, int(NEMI+0.0001)-1:int(NEMX+0.0001)+1) )
            HISNV  = 0
            HISNE  = 0
            HISMAP = 0
            PNVMI = 1.D64
            PNVMX = 0
            rewind(hFileSwap)
            do I=1, NREC
               read(hFileSwap, *, ERR=100)  NEWNV, NEWEV
               HISNV(NEWNV) = HISNV(NEWNV) + 1
               HISNE(NEWEV) = HISNE(NEWEV) + 1
               HISMAP(NEWNV, NEWEV) = HISMAP(NEWNV, NEWEV) + 1
               PNV = dble(NEWNV)/dble(NEWEV)
               if(PNVMI .gt. PNV) PNVMI = PNV
               if(PNVMX .le. PNV) PNVMX = PNV
            end do

     100    rewind(hFileSwap)
            allocate(HISNVPE(100))
            DPNV = (PNVMX - PNVMI)/dble(size(HISNVPE))
            HISNVPE = 0
            do I=1, NREC
               read(hFileSwap, *, ERR=200)  NEWNV, NEWEV
               PNV = dble(NEWNV)/dble(NEWEV)
               IBIN = (PNV-PNVMI)/DPNV + 1
               HISNVPE(IBIN) = HISNVPE(IBIN) + 1
            end do

     200    do I= LBOUND(HISNV,dim=1), UBOUND(HISNV,dim=1)
               write(hFile1, fmt="(2(I8,1x))") I, HISNV(I)
            end do
            close(hFile1)

            do I= LBOUND(HISNE, dim=1), UBOUND(HISNE, dim=1)
               write(hFile2, fmt="(2(I8,1x))") I, HISNE(I)
            end do
            close(hFile2)

            do I= LBOUND(HISNVPE, dim=1), UBOUND(HISNVPE, dim=1)
               write(hFile3, fmt="((1PE13.5,1x, I8))") PNVMI+(I-0.5)*DPNV, HISNVPE(I)
            end do
            close(hFile3)

            do I= LBOUND(HISMAP, dim=1), UBOUND(HISMAP, dim=1)
               write(hFile4, fmt="(2000(I6,1x))") HISMAP(I, LBOUND(HISMAP, dim=2):UBOUND(HISMAP, dim=2))
            end do
            close(hFile4)


            close(hFileSwap)

          return
  end subroutine SEGLEN_For_One_MC_Test
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use RAND32_MODULE
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        character*256::GFILE
!-----  Local varaible
        integer::I, IB, ERR, NPRT0, IS0
        real(KINDDF)::RC


            NPRT0 =  SimBox(1)%NPRT
           !--- allocate memory and calculate the Neighbor sites
           if(.not. allocated(m_XPL)) then
              allocate(m_XPL(NPRT0,3), m_INDL(NPRT0,mp_NPT), STAT=ERR)
              if(ERR) then
                write(*,fmt='(A)') "MDPSCU Error: Fail to allocate memory for neighbor list in MC_VisitLattice_HPC_module"
                write(*,fmt='(A)') "              You may want use less box-cells"
                write(*,fmt='(A)') "              Process to be stoped"
                stop
               end if
               call CreateNeighborSites(SimBox(1), CtrlParam, m_INDL, m_XPL)
               print *, "Neighbore sites created"
            end if

            !--- determine the initial site in template box
            RC = 1.D64
            do I=1, size(m_XPL,dim=1)
               if(sum(m_XPL(I,1:3)*m_XPL(I,1:3)) .lt. RC) then
                  RC  = sum(m_XPL(I,1:3)*m_XPL(I,1:3))
                  IS0 = I
                end if
            end do

            call SEGLEN_For_One_MC_Test(SimBox(1), CtrlParam, IS0)
            return


            m_cumVISIT   = 0
            m_cumEVENT = 0
            do IB=1, m_NTEST
               print *, "BOX #", IB
               call NPE_For_One_MC_Test(SimBox(1), CtrlParam, IB, IS0)
            end do !--- end loop for box

          return
  end subroutine MyRECORD
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
       integer::I

          if(m_hFileSum .gt. 0) close(m_hFileSum)

          if(allocated(m_XPL))  deallocate(m_XPL)
          if(allocated(m_INDL)) deallocate(m_INDL)
          if(allocated(m_SMD))  deallocate(m_SMD)

          return
  end subroutine MyCleaner
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyAfterRECORD(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
          !call MyOutput(SimBox, CtrlParam)
          call MyCleaner(SimBox, CtrlParam)
          return
  end subroutine MyAfterRECORD
  !**********************************************************************************
 end module MC_VisitLattice_HPC_module


 !**********************************************************************************
 Program MC_VisitLattice_HPC_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use MC_VisitLattice_HPC_module

 implicit none

       call APPSHELL_AddRecord( PRERECORD= MyInitialize, &
                                RECORDPROC=MyRECORD,     &
                                AFTRECORD= MyAfterRECORD)

       call Single_ANALYSIS(0)

       stop
 End program MC_VisitLattice_HPC_Main
