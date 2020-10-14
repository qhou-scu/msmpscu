 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to create  tetrahdral lattice for a BCC crystal. The BCC lattice
 !                  is to be loaded from a configuration file.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                  CreateTetraLattice_BCC.exe -s "filename"
 !                  where:
 !                        filename  - the name of the file storing the coordinates of BCC lattices.
 !                                    The file could be created by, for example MakeBox tool, or
 !                                    a configuration file created by MDPSCU
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-06 (Hou Qing, Sichuan university)
 !
 !

 module CreateTetraLattice_BCC
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none

      !--  the axsis in crystal
      integer, private:: m_AXSISX(3) = (/1, 0, 0/)
      integer, private:: m_AXSISY(3) = (/0, 1, 0/)
      integer, private:: m_AXSISZ(3) = (/0, 0, 1/)
      integer, parameter, private::mp_NPT=6
      integer, parameter, private::mp_NPT2=mp_NPT*2
      integer, parameter, private::mp_NPT4=mp_NPT*4
      real(KINDDF),parameter, private::mp_PVECTX(mp_NPT)=(/1.0D0, -1.0D0, 0.0D0,  0.0D0,  0.0D0,  0.0D0/)
      real(KINDDF),parameter, private::mp_PVECTY(mp_NPT)=(/0.0D0,  0.0D0, 1.0D0, -1.0D0,  0.0D0,  0.0D0/)
      real(KINDDF),parameter, private::mp_PVECTZ(mp_NPT)=(/0.0D0,  0.0D0, 0.0D0,  0.0D0,  1.0D0, -1.0D0/)

      integer,parameter, private::mp_PIAR(mp_NPT4) =(/1, 3,  1, 4,  1, 5,  1, 6,  3, 5,  3, 6, &
                                                      2, 3,  2, 4,  2, 5,  2, 6,  4, 5,  4, 6/)

      real(KINDDF), private::m_RPVECTX(mp_NPT)
      real(KINDDF), private::m_RPVECTY(mp_NPT)
      real(KINDDF), private::m_RPVECTZ(mp_NPT)
      real(KINDDF), private::m_EPSLON = 1.D-3
  contains

  !**********************************************************************************
   subroutine Initialize_TetraSite(SimBox, CtrlParam)
   !***  PURPOSE:  to load control parameters and allocate memories needed
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I
   real(KINDDF)::RMAT(3,3)

          CtrlParam%RU      = SimBox%RR*1.2
          CtrlParam%NB_RM   = CtrlParam%RU
          CtrlParam%IFPD    = 1

          !$$--- to create the crystal coordinates
           RMAT(1,1:3) = m_AXSISX(1:3)/dsqrt(sum( dble(m_AXSISX(1:3)) * dble(m_AXSISX(1:3))))
           RMAT(2,1:3) = m_AXSISY(1:3)/dsqrt(sum( dble(m_AXSISY(1:3)) * dble(m_AXSISY(1:3))))
           RMAT(3,1:3) = m_AXSISZ(1:3)/dsqrt(sum( dble(m_AXSISZ(1:3)) * dble(m_AXSISZ(1:3))))
           do I=1, mp_NPT
              m_RPVECTX(I) = mp_PVECTX(I)*RMAT(1,1) + mp_PVECTY(I)*RMAT(1,2) + mp_PVECTZ(I)*RMAT(1,3)
              m_RPVECTY(I) = mp_PVECTX(I)*RMAT(2,1) + mp_PVECTY(I)*RMAT(2,2) + mp_PVECTZ(I)*RMAT(2,3)
              m_RPVECTZ(I) = mp_PVECTX(I)*RMAT(3,1) + mp_PVECTY(I)*RMAT(3,2) + mp_PVECTZ(I)*RMAT(3,3)
           end do
         return
   end subroutine Initialize_TetraSite
  !**********************************************************************************

  !**********************************************************************************
  subroutine CreateTetraSites(SimBox, CtrlParam, TN, TXP)
  !***  DESCRIPTION: to create theTetrahdralSites
  !
  !    INPUT:  SimBox, the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    OUTPUT: TN,     the number of  tetrahdral sites
  !            TXP,    the coordinates of the tetrahdral sites
  !
  use MD_NeighborsList
  implicit none
       type(SimMDBox),             intent(in)::SimBox
       type(SimMDCtrl),            intent(in)::CtrlParam
       integer                               ::TN
       real(KINDDF), dimension(:,:)          ::TXP

!-----  Local varaible
        type(NEIGHBOR_LIST)::List
        integer::I, IP, J, K, NPRT, NNA, IAP, IA(mp_NPT), J1, J2, TN0
        integer,      dimension(:),   allocatable::ISP
        integer,      dimension(:,:), allocatable::IND
        real(KINDDF), dimension(:,:), allocatable::XP
        real(KINDDF)::RR, DIST, SEP(3), XPP(3), BL(3), BU(3), BS(3), HBS(3), PTX(mp_NPT), PTY(mp_NPT), PTZ(mp_NPT), RC(mp_NPT)


            !$$--- create neighbor-list
            call Cal_NeighboreList2C(SimBox, CtrlParam, List)
            NPRT = SimBox%NPRT
            allocate(IND(NPRT, mp_NPT), XP(NPRT,3), ISP(NPRT*mp_NPT4))

            !$$--- get the secondary neighbors for each atoms
            IND = 0
            RR  = SimBox%RR
            BS  = SimBox%ZL/RR
            HBS = 0.5D0*BS
            BL  = SimBox%BOXLOW/RR
            BU  = SimBox%BOXUP/RR
            XP(1:NPRT, 1:3) = SimBox%XP(1:NPRT, 1:3)/RR
            do I=1, NPRT
               !$$--- calculate the static vertice points of the Tetraheral
               do IP=1, mp_NPT
                  PTX(IP) = XP(I,1) + m_RPVECTX(IP)
                  PTY(IP) = XP(I,2) + m_RPVECTY(IP)
                  PTZ(IP) = XP(I,3) + m_RPVECTZ(IP)
                  RC(IP)  = 1.D32
               end do

               !$$--- select out the lattice positions closest to the vertice points
               NNA = List%KVOIS(I)
               do J=1, NNA
                  IAP = List%INDI(I, J)
                  XPP(1:3) = XP(IAP,1:3)
                  if(dabs(XPP(1) - XP(I,1)) .gt. HBS(1)) &
                     XPP(1) = XPP(1) - DSIGN(BS(1),XPP(1) - XP(I,1))
                  if(dabs(XPP(2) - XP(I,2)) .gt. HBS(2)) &
                     XPP(2) = XPP(2) - DSIGN(BS(2),XPP(2) - XP(I,2))
                  if(dabs(XPP(3) - XP(I,3)) .gt. HBS(3)) &
                     XPP(3) = XPP(3) - DSIGN(BS(3),XPP(3) - XP(I,3))

                  do IP=1, mp_NPT
                      SEP(1) = XPP(1) - PTX(IP)
                      SEP(2) = XPP(2) - PTY(IP)
                      SEP(3) = XPP(3) - PTZ(IP)
                      DIST = sum(SEP(1:3)*SEP(1:3) )
                      if(DIST .lt. RC(IP)) then
                        IA(IP) = J
                        RC(IP) = DIST
                      end if
                  end do
               end do

               do IP=1, mp_NPT
                  IND(I,IP) =  List%INDI(I, IA(IP))
               end do
            end do

            !$$--- start construct the tetrahdral sites
            TN0 = 0
            do I=1, NPRT
               !$$--- calculate the instant vertice points of the Tetraheral
               do IP=1, mp_NPT
                  J = IND(I, IP)
                  IA(IP) = J

                  PTX(IP) = XP(J, 1) - XP(I,1)
                  if(dabs(PTX(IP)) .gt. HBS(1)) PTX(IP) = PTX(IP) - DSIGN(BS(1), PTX(IP))

                  PTY(IP) = XP(J, 2) - XP(I,2)
                  if(dabs(PTY(IP)) .gt. HBS(2)) PTY(IP) = PTY(IP) - DSIGN(BS(2), PTY(IP))

                  PTZ(IP) = XP(J, 3) - XP(I,3)
                  if(dabs(PTZ(IP)) .gt. HBS(3)) PTZ(IP) = PTZ(IP) - DSIGN(BS(3), PTZ(IP))
               end do

               !$$--- the sites formed by PTs
                do K=1, mp_NPT2
                   J1  = mp_PIAR(2*K-1)
                   J2  = mp_PIAR(2*K)
                   TN0 = TN0 + 1
                   TXP(TN0, 1) = 0.5D0*PTX(J1) + 0.25D0*PTX(J2) + XP(I,1)
                   TXP(TN0, 2) = 0.5D0*PTY(J1) + 0.25D0*PTY(J2) + XP(I,2)
                   TXP(TN0, 3) = 0.5D0*PTZ(J1) + 0.25D0*PTZ(J2) + XP(I,3)

                   TN0 = TN0 + 1
                   TXP(TN0, 1) = 0.25D0*PTX(J1) + 0.5D0*PTX(J2) + XP(I,1)
                   TXP(TN0, 2) = 0.25D0*PTY(J1) + 0.5D0*PTY(J2) + XP(I,2)
                   TXP(TN0, 3) = 0.25D0*PTZ(J1) + 0.5D0*PTZ(J2) + XP(I,3)
                end do
            end do

            !$$--- to delete the same sites
            ISP = 1
            do I=1, TN0
               if(ISP(I) .le. 0) cycle
               do J=I+1, TN0
                  if(ISP(J) .le. 0) cycle

                  SEP(1:3) = TXP(J,1:3) - TXP(I,1:3)
                  if(dabs(SEP(1)) .gt. HBS(1)) SEP(1) = SEP(1) - DSIGN(BS(1), SEP(1))
                  if(dabs(SEP(2)) .gt. HBS(2)) SEP(2) = SEP(2) - DSIGN(BS(2), SEP(2))
                  if(dabs(SEP(3)) .gt. HBS(3)) SEP(3) = SEP(3) - DSIGN(BS(3), SEP(3))
                  if(sum(SEP*SEP) .lt. m_EPSLON ) then
                     ISP(J) = 0
                  end if
               end do
            end do

            !$$--- to remove the sites out side the box into box
            TN = 0
            do I=1, TN0
               if(ISP(I) .gt. 0) then
                  if(TXP(I, 1) .lt. BL(1)) TXP(I, 1) = TXP(I, 1) + BS(1)
                  if(TXP(I, 1) .ge. BU(1)) TXP(I, 1) = TXP(I, 1) - BS(1)

                  if(TXP(I, 2) .lt. BL(2)) TXP(I, 2) = TXP(I, 2) + BS(2)
                  if(TXP(I, 2) .ge. BU(2)) TXP(I, 2) = TXP(I, 2) - BS(2)

                  if(TXP(I, 3) .lt. BL(3)) TXP(I, 3) = TXP(I, 3) + BS(3)
                  if(TXP(I, 3) .ge. BU(3)) TXP(I, 3) = TXP(I, 3) - BS(3)

                  TN = TN + 1
                  TXP(TN,1:3) = TXP(I,1:3)*RR
               end if
            end do
            deallocate(IND, XP, ISP)
          return
  end subroutine CreateTetraSites
  !**********************************************************************************

  !**********************************************************************************
  subroutine CreateTetraSiteBox(SimBox, CtrlParam, TetraBox)
  !***  DESCRIPTION: to create the TetrahdralSites box
  !
  !    INPUT:  SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    OUTPUT: TetraBox,     the simulation box with atom position at the tetrahedral sites
  !                         of box SimBox
  !
  use MD_NeighborsList
  implicit none
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam
       type(SimMDBox)             ::TetraBox

!-----  Local varaible
       integer::TN
        real(KINDDF), dimension(:,:), allocatable::TXP

            allocate(TXP(SimBox%NPRT*mp_NPT4,3))
            call CreateTetraSites(SimBox, CtrlParam, TN, TXP)
            call AddAtoms_SimMDBox(TetraBox, TN, ITYPE=SimBox%NGROUP+1, RXP=TXP, TYPEORDER=0)
            deallocate(TXP)
            return
  end subroutine CreateTetraSiteBox
  !**********************************************************************************

  !**********************************************************************************
  subroutine Record_TetraSite(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_NeighborsList
  implicit none
       type(MDRecordStamp),          intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        character*256::GFILE
!-----  Local varaible
        type(NEIGHBOR_LIST)::List
        integer::I, TN, NBOX
        type(SimMDBox)::SWAPBOX
        type(MDRecordStamp)::SWAPSTAMP

            call Copy_SimMDBox(SimBox(1), SWAPBOX)
            call CreateTetraSiteBox(SimBox(1), CtrlParam, SWAPBOX)
            SWAPSTAMP%ITime = 0
            SWAPSTAMP%ITEST = 1
            SWAPSTAMP%ISect = 1
            SWAPSTAMP%ICfg  = 0
            SWAPSTAMP%IRec  = 0
            SWAPSTAMP%Time  = 0
            write(*,fmt="(A, I6, A)") " Tetrahedral latiice to be output to "//&
                                       gm_outFileName(1)(1:len_trim(gm_outFileName(1)))//".0000"
            call Putout_Instance_Config_SimMDBox(gm_outFileName(1), SWAPBOX,SWAPSTAMP)
            call Release_SimMDBox(SWAPBOX)
          return
  end subroutine Record_TetraSite
  !**********************************************************************************
 end module CreateTetraLattice_BCC
