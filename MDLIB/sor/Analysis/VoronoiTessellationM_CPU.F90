  module VoronoiTessellation_CPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  To do the Voronoi tessellation of atoms base on the algorithm
  !                  of Medvedev (J. Comput. Phys. 67(1986)223-229)
  !
  !                  The difference between this code and VoronoiTessellationM.F90 is
  !                  that this is a CPU version, while VoronoiTessellationM.cuf is a GPU
  !                  version.
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.cuf
  !                       MD_Globle_Variables_GPU.cuf
  !
  !                  REFERENCE____________________________________________________________________________
  !                       Medvedev, J. Comput. Phys. 67(1986)223-229
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2013-08 (Hou Qing)
  !
  !________________________________________________________________________________________________________



  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_ForceTable
  !-----------------------------------------------
    implicit none

      integer::m_processid=0
      !$$--- the id of  filename get from SimMDCtrl for I/O.
      !$$    Refering the definition of f_others in SimMDCtrl.
      !$$    By dafualt, m_VT_OFILE = 1. If the keyword "&AUXF_VTO"
      !$$    is found in member array f_tag of SimMDCtrl,
      !$$    m_VT_OFILE will be reassigned corresponding value.
      integer::m_VT_IFILE =-1
      character(len=9),parameter, private::mp_FTAGI="&AUXF_VTI"
      integer::m_VT_OFILE = 1
      character(len=9),parameter, private::mp_FTAGO="&AUXF_VTO"

      !--- control parameters determine what will be calculated
      !    and ouput
      integer, private::m_OUTPUTVOL = 1                      ! if >0, atomic volum will be calculated and output
      integer, private::m_OUTPUTDT  = 0                      ! if >0, DT will be ouput
      integer, private::m_OUTPUTVT  = 0                      ! if >0, VT will be ouput
      integer, dimension(:), allocatable, private::m_CONTOUR ! if m_CONTOUR(I)>0, the 3D contour formed by atom type I will be calculated

      !--- calculation control parameters
      integer, private::m_INITED = 0

      !--- box information
      real(KINDDF), private::m_BOXSHAPE(3,3)                 ! siumulation box shape
      real(KINDDF), private::m_BOXSIZE(3)                    ! simulation box size
      integer, private::m_IFPD(3)                            ! the periodic condition

      !--- Voronoi tessellation vertices of atoms  on host
      integer, parameter, private::mp_INITED_VOLS = 1
      integer, parameter, private::mp_INITED_VVTS = 2
      integer, private::m_INITED_V = 0
      integer, parameter, private::mp_MXNA    = 128*128      !Max number of atoms included in a call of tesellation calculation
      integer, parameter, private::mp_MXNF    = 40           !Max permitted number of faces of a Voronio volume
      integer, parameter, private::mp_MXNV    = mp_MXNF      !Max permitted number of vertice of a face
      integer, parameter, private::mp_MXKVOIS = 800

      integer, private::m_MXNF    = mp_MXNF                  !Actual max number  of faces of Voronio volumes
      integer, private::m_MXNV    = mp_MXNV                  !Actual max number  of vertice of faces of Voronio volumes
      integer, private::m_MXKVOIS = mp_MXKVOIS


      !--- Delaunay tessellation vertices
      integer,    dimension(:),  allocatable,private::m_DVTS
      !--- Working spaces
      integer,    dimension(:),  allocatable,private::m_WVFS
      integer(1), dimension(:),  allocatable, private::m_WFLG

      !--- Voronoi volum and vertices
      real(KINDSF), dimension(:),  allocatable::m_VOLS
      real(KINDSF), dimension(:,:),  allocatable::m_VVTS

      !--- states of Voronoi volumes
      integer, parameter::mp_VSTAT_ENCLOSED = 0
      integer, parameter::mp_VSTAT_OPENED   = 1
      integer(1), dimension(:),  allocatable::m_VSTAT
      !----
  contains

  !**********************************************************************************
  subroutine Clear_VoronoiTessellation(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the allocated memories
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
      !--- Local variables

      !***
      if(allocated(m_DVTS)) deallocate(m_DVTS)
      if(allocated(m_WVFS)) deallocate(m_WVFS)
      if(allocated(m_WFLG)) deallocate(m_WFLG)
      if(allocated(m_VOLS)) deallocate(m_VOLS)
      if(allocated(m_VVTS)) deallocate(m_VVTS)
      if(allocated(m_VSTAT)) deallocate(m_VSTAT)

      RETURN
  end subroutine Clear_VoronoiTessellation
  !*********************************************************************************

   !*********************************************************************************
    SUBROUTINE LoadControlParameters(fname)
    !***  PURPOSE:   to readin control parameters from a file
    !     INPUT:     fname: the file name
    !
    !     OUTPUT:
    !
    use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables
    implicit none
    character*(*)::fname
    integer::hFile, N, I, LINE
    character*256::STR
    character*32::STRNUMB(10)

            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old')

            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,1, N, STRNUMB)
            m_OUTPUTVOL = ISTR(STRNUMB(1))

            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,1, N, STRNUMB)
            m_OUTPUTVT  = ISTR(STRNUMB(1))

            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,1, N, STRNUMB)
            m_OUTPUTDT  = ISTR(STRNUMB(1))

            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,size(m_CONTOUR), N, STRNUMB)
            DO I=1, N
               m_CONTOUR(I) = ISTR(STRNUMB(I))
            END DO
            close(hFile)

    100     write(*,*) "!************ Information for Voronoi tessellation  **********"
            write(*,FMT="(' !    output atomic Voronoi volume:           ', 3(I3,1x))") m_OUTPUTVOL
            write(*,FMT="(' !    output intermediate Voronoi vertice:    ', 3(I3,1x))") m_OUTPUTVT
            write(*,FMT="(' !    output intermediate Delaunay vertices:  ', 3(I3,1x))") m_OUTPUTDT
            !DO I=1, size(m_CONTOUR)
            !   write(*,FMT="(' !    output contour of clusters: ', 3(I3,1x))") I, m_CONTOUR(I)
            !END DO

            return
    END SUBROUTINE LoadControlParameters
  !*********************************************************************************


  !*********************************************************************************
    SUBROUTINE Initialize_VoronoiVolume(SimBox, CtrlParam, NONE)
    !***  PURPOSE:   to allocate the primaritive memory for Voronoi Tesselation
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !
    !     OUTPUT:   the working spaces allocated
    !
    implicit none
      !--- dummy vaiables
      type(SimMDBox), intent(in)::SimBox
      type(SimMDCtrl),intent(in)::CtrlParam
      OPTIONAL::NONE
      external::NONE

      !--- Local variables
      integer::NPRT, I, ERR

           if(m_INITED .gt. 0) then
             call Clear_VoronoiTessellation(SimBox, CtrlParam)
           end if
           allocate(m_CONTOUR(SimBox%nGroup))
           m_CONTOUR = 0

           !$$--- to findout the I/O unit
           DO I=1, SIZE(CtrlParam%f_tag)
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                 m_VT_OFILE = I
              end if
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                 m_VT_IFILE = I
                 call LoadControlParameters(CtrlParam%f_others(m_VT_IFILE))
              end if
           END DO

           !if(CtrlParam%mKVOIS .GT. mp_MXNN) then
           !   write(*,*) "MDPSCU Warning: the permitted max number of neighbors larger than mp_MXNN in Voronoi tessellation."
           !   write(*,*) "                the m_MXNN is set to the smaller value.", min(CtrlParam%mKVOIS,mp_MXNN)
           !      call ONWARNING(gm_OnWarning)
           !end if
           m_MXNF    = min(CtrlParam%NB_MXNBS,mp_MXNF)
           m_MXKVOIS = min(CtrlParam%NB_MXNBS, mp_MXKVOIS)

           NPRT    = SimBox%NPRT
           allocate(m_DVTS(mp_MXNA*m_MXNF*m_MXNV),  m_WVFS(m_MXNF), m_WFLG(m_MXKVOIS), STAT=err)
           if(ERR) goto 100

           allocate(m_VVTS(mp_MXNA*m_MXNF*m_MXNV,3), m_VOLS(NPRT), m_VSTAT(NPRT), STAT=err)
           if(ERR) goto 100

           m_INITED = .TRUE.
           m_INITED_V = IOR(m_INITED_V, mp_INITED_VOLS)
           m_INITED_V = IOR(m_INITED_V, mp_INITED_VVTS)
           return

     100   write(*,*) "MDPSCU Error: fail to allocate memories for Voronio tessellation"
           stop

           return
    END SUBROUTINE Initialize_VoronoiVolume
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Cal_Delaunay_Vertice(IM, XP, IA0, NPRT, MXNF, KVOIS,INDI, MXNV, VERT, STAT, FACE, FLAG)
  !***  PURPOSE:   to create the primary faces and vertices of Delaunay tesstellation
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              XP:     the positions of atoms
  !$$              IA0:    the index (in the whole box)of the fisrt atom on in this call
  !$$              NPRT:   the actual number of atoms on the in this call
  !$$              MXNF:   the max permitted number of neighbores for each atom
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$
  !$$  OUTPUT:     VERT:   the neighbore ID of atoms on Delaunay vertices
  !$$              STAT:   the statu indicating if the tessellation is end closed
  !$$              FACE, FLAG: working space
  !$$
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, intent(in)::IM, IA0, NPRT, MXNF, MXNV
      real(KINDDF), dimension(:,:), intent(in)::XP
      integer,      dimension(:),   intent(in)::KVOIS
      integer,      dimension(:,:), intent(in)::INDI
      integer,      dimension(:)              ::VERT, FACE
      integer(1),   dimension(:)              ::FLAG, STAT

      !Local variables
      real(KINDSF)::RR1, RR2, RR3, R2T, R2MI, XC, YC, ZC, XCP, YCP, ZCP, NMX, NMY, NMZ, NMD
      real(KINDSF)::XY12, YZ12, ZX12, RX12, RY12, RZ12, DET2
      real(KINDSF)::SEPX, SEPY, SEPZ
      real(KINDSF)::X0, Y0, Z0
      real(KINDSF)::X1, Y1, Z1
      real(KINDSF)::X2, Y2, Z2
      real(KINDSF)::X3, Y3, Z3
      real(KINDSF)::X4, Y4, Z4
      real(KINDSF)::SEP(mp_MXKVOIS,3)

      real(KINDSF)::HBX, HBY, HBZ, BX, BY, BZ
      integer::IFPDX, IFPDY, IFPDZ

      ! Local variables
      real(KINDSF), parameter::scal=1.E8, underflow=1.E-8, rescal=0.5E-8
      integer::IA, STATU, MXNV2
      integer::J, IW, IW1, IIW, IFACE0, IVERT0, IFLG0, IVERT, NVERT,I1, I2, I3, I4, LF

             !$--- scal the coordination to prevent underflow
             BX = m_BOXSIZE(1)*scal
             BY = m_BOXSIZE(2)*scal
             BZ = m_BOXSIZE(3)*scal

             HBX = BX*C_HALF
             HBY = BY*C_HALF
             HBZ = BZ*C_HALF

             IFPDX = m_IFPD(1)
             IFPDY = m_IFPD(2)
             IFPDZ = m_IFPD(3)
             MXNV2 = MXNV-2
             VERT = 0

             DO IA=1, NPRT
                !NOTE: in CPU version, size of VERT is define
                !VERT(IVERT0+1:IVERT0+MXNN*MXNV) = 0
                !$$--- get the position of the atom IC
                !$$    to prevent underflow, the coordination are scaled
                IIW = KVOIS(IA+IA0)
                if(IIW .LE. 0) then
                   STAT(IA+IA0) =  mp_VSTAT_OPENED
                   cycle
                end if

                IVERT0 = (IA-1)*MXNF*MXNV
                X0 = XP(IA+IA0,1)*scal
                Y0 = XP(IA+IA0,2)*scal
                Z0 = XP(IA+IA0,3)*scal

                DO IW =1, min(IIW, mp_MXKVOIS)
                   J=INDI(IA+IA0,IW)
                   SEP(IW,1) = XP(J,1)*scal - X0
                   SEP(IW,2) = XP(J,2)*scal - Y0
                   SEP(IW,3) = XP(J,3)*scal - Z0
                   if(IFPDX.GT.0) then
                      if(ABS(SEP(IW,1)) .GT. HBX) SEP(IW,1) = SEP(IW,1) - SIGN(BX,SEP(IW,1))
                   end if

                   if(IFPDY.GT.0) then
                      if(ABS(SEP(IW,2)) .GT. HBY) SEP(IW,2) = SEP(IW,2) - SIGN(BY,SEP(IW,2))
                   end if

                   if(IFPDZ.GT.0) then
                      if(ABS(SEP(IW,3)) .GT. HBZ) SEP(IW,3) = SEP(IW,3) - SIGN(BZ,SEP(IW,3))
                   end if
                 END DO

                 !$$--- STEP1: scan the neighbores to findout the nearest neighbore of atom IA
                 !$$           and thus the first face of VT
                  R2MI = 1.D32
                  DO IW=1, min(IIW, mp_MXKVOIS)
                     !$$--- To calculate the seperation between particle IA
                     !$$    and its IWth neighbore
                     SEPX = SEP(IW,1)
                     SEPY = SEP(IW,2)
                     SEPZ = SEP(IW,3)
                     RR1 = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                     !$$--- NOTE:DO we need to accont only the neighbore in potential cutoff reqion?
                     if(RR1 .LE. R2MI) then
                        R2MI = RR1
                        I1   = IW
                        X1 = SEPX
                        Y1 = SEPY
                        Z1 = SEPZ
                     end if
                  ENDDO
                 !$$--- store the the associated atom ID on the first face
                  VERT(IVERT0+1)  = I1
                  RR1 = R2MI

                 !$$--- STEP2: scan the neighbores to findout the starting point
                 !$$
                  R2MI = 1.D32
                  DO IW=1, min(IIW, mp_MXKVOIS)
                     if(IW .EQ. I1) cycle

                     SEPX = SEP(IW,1)
                     SEPY = SEP(IW,2)
                     SEPZ = SEP(IW,3)
                     RR2 = (SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ)

                     !$$--- the foot from half*(X1, Y1, Z1) to the intersecting line
                     !$$    of face I1 and I2
                     XY12 = X1*SEPY - Y1*SEPX
                     YZ12 = Y1*SEPZ - Z1*SEPY
                     ZX12 = Z1*SEPX - X1*SEPZ
                     RX12 = RR1*SEPX - RR2*X1
                     RY12 = RR1*SEPY - RR2*Y1
                     RZ12 = RR1*SEPZ - RR2*Z1
                     DET2 = XY12*XY12 + YZ12*YZ12 + ZX12*ZX12
                     if(ABS(DET2) .LE. underflow) cycle

                     XC   = (RY12*XY12 - RZ12*ZX12)/DET2 - X1
                     YC   = (RZ12*YZ12 - RX12*XY12)/DET2 - Y1
                     ZC   = (RX12*ZX12 - RY12*YZ12)/DET2 - Z1

                     R2T  = XC*XC + YC*YC + ZC*ZC
                     if(R2T .LT. R2MI) then
                        R2MI  = R2T
                        I2    = IW
                        X2 = SEPX
                        Y2 = SEPY
                        Z2 = SEPZ
                     end if
                  ENDDO
                  VERT(IVERT0+2)  = I2

                  !$$--- get the real project foot:(X1, Y1, Z1)
                  !$$    on face I2
                  !$$    NOTE: the factor C_HALF is not included
                  RR2  = (X2*X2 + Y2*Y2 + Z2*Z2)
                  XY12 = X1*Y2 - X2*Y1
                  YZ12 = Y1*Z2 - Y2*Z1
                  ZX12 = Z1*X2 - Z2*X1
                  RX12 = RR1*X2 - RR2*X1
                  RY12 = RR1*Y2 - RR2*Y1
                  RZ12 = RR1*Z2 - RR2*Z1
                  DET2 = XY12*XY12 + YZ12*YZ12 + ZX12*ZX12

                  XC   = (RY12*XY12 - RZ12*ZX12)/DET2
                  YC   = (RZ12*YZ12 - RX12*XY12)/DET2
                  ZC   = (RX12*ZX12 - RY12*YZ12)/DET2

                  !$$--- STEP3; from the foot to findout the first vertex
                  !$$           on the intersecting lin of face I1 and I2
                   R2MI = 1.D32
                   DO IW=1, min(IIW, mp_MXKVOIS)
                      if(IW .EQ. I1 .OR. IW.EQ.I2) cycle
                      SEPX = SEP(IW,1)
                      SEPY = SEP(IW,2)
                      SEPZ = SEP(IW,3)
                      RR3 = (SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ)

                      DET2 = SEPX*YZ12 + SEPY*ZX12 + SEPZ*XY12
                      XCP = (RR3*YZ12 + RY12*SEPZ - RZ12*SEPY)/DET2 - XC
                      YCP = (RR3*ZX12 + RZ12*SEPX - RX12*SEPZ)/DET2 - YC
                      ZCP = (RR3*XY12 + RX12*SEPY - RY12*SEPX)/DET2 - ZC
                      R2T =  XCP*XCP + YCP*YCP + ZCP*ZCP
                      if(R2T .LT. R2MI) then
                         R2MI  = R2T
                         I3 = IW
                      end if
                   END DO
                   VERT(IVERT0+3) = I3

                   !$$--- I2 is also the second face ID and the first associated atom
                   !$$    of the second face
                   !$$    I1 , I3, are also the assocaied atoms on the second face.
                   I4 = IVERT0+MXNV
                   VERT(I4+1) = I2
                   VERT(I4+2) = I3
                   VERT(I4+3) = I1

                   !$$--- I3 is also the second face ID and the first associated atom
                   !$$    of the third face
                   !$$    I2 , I3, are also the associated atoms on the third face.
                   I4 = I4 + MXNV
                   VERT(I4+1) = I3
                   VERT(I4+2) = I1
                   VERT(I4+3) = I2

                   !$$ reset working space
                   FACE = 0
                   FACE(1:3) = VERT(IVERT0+1:IVERT0+3)
                   STATU =  mp_VSTAT_ENCLOSED
                   !$$**********************************************************
                   !$$---  Now we have the first three faces and starting vertex
                   !$$     scan the faces to construt the faces
                    DO LF =1, min(IIW, MXNF)
                       !$$--- check if we have further face needed to be constructed
                       if(FACE(LF) .LE. 0) exit

                       !$$--- get the position of the first vertex (XC, YC, ZC) on face LF
                       IVERT = IVERT0 + (LF-1)*MXNV
                       I1 = VERT(IVERT + 1)
                       I2 = VERT(IVERT + 2)
                       I3 = VERT(IVERT + 3)

                       !$$--- reset the working space for visited neighbors
                       FLAG = 0
                       FLAG(I1) = 1
                       FLAG(I2) = 1
                       FLAG(I3) = 1

                       !$$--- get the position (Xc, YC, ZC0 of the first vertex
                       X1 = SEP(I1,1)
                       Y1 = SEP(I1,2)
                       Z1 = SEP(I1,3)
                       X2 = SEP(I2,1)
                       Y2 = SEP(I2,2)
                       Z2 = SEP(I2,3)
                       X3 = SEP(I3,1)
                       Y3 = SEP(I3,2)
                       Z3 = SEP(I3,3)
                       RR1  = (X1*X1 + Y1*Y1 + Z1*Z1)
                       RR2  = (X2*X2 + Y2*Y2 + Z2*Z2)
                       RR3  = (X3*X3 + Y3*Y3 + Z3*Z3)
                       DET2 =  Y2*(X1*Z3 - X3*Z1)  + Y3*(X2*Z1 - X1*Z2)  + Y1*(X3*Z2 - X2*Z3)
                       XC   = (RR1*(Y2*Z3 - Y3*Z2) + RR2*(Y3*Z1 - Y1*Z3) + RR3*(Y1*Z2 - Y2*Z1))/DET2
                       YC   = (RR1*(Z2*X3 - Z3*X2) + RR2*(Z3*X1 - Z1*X3) + RR3*(Z1*X2 - Z2*X1))/DET2
                       ZC   = (RR1*(X2*Y3 - X3*Y2) + RR2*(X3*Y1 - X1*Y3) + RR3*(X1*Y2 - X2*Y1))/DET2

                      !$$--- prepare to find the next vertex
                      !
                      !$$--- the normal of the face I2
                       NMD = SQRT(RR2)
                       NMX = X2/NMD
                       NMY = Y2/NMD
                       NMZ = Z2/NMD

                      !$$--- the determinants to be used to find the next vertex
                       X2  = X3
                       Y2  = Y3
                       Z2  = Z3
                       RR2 = RR3
                       XY12 = X1*Y2 - X2*Y1
                       YZ12 = Y1*Z2 - Y2*Z1
                       ZX12 = Z1*X2 - Z2*X1
                       RX12 = RR1*X2 - RR2*X1
                       RY12 = RR1*Y2 - RR2*Y1
                       RZ12 = RR1*Z2 - RR2*Z1

                       !$$--- STEP4: scan the neighbores to construct other vertex
                        IVERT   = IVERT + 3
                        NVERT   = 1
                        DO IW1=1, MXNV
                           R2MI = 1.E32
                           I4   = -1
                           DO IW=1, min(IIW, mp_MXKVOIS)
                              !$$--- check if the neighbor has been on the vertex list
                               if(NVERT .LT. 2) then
                                  if(FLAG(IW) .GT. 0) cycle
                               else
                                  if(FLAG(IW) .GT. 0 .AND. IW.NE.I2) cycle
                               endif

                              !$$--- IW is not still on the vertex list
                               SEPX = SEP(IW,1)
                               SEPY = SEP(IW,2)
                               SEPZ = SEP(IW,3)

                              !$$---
                               RR3  = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                               DET2 = SEPX*YZ12 + SEPY*ZX12 + SEPZ*XY12
                               XCP  = (RR3*YZ12 + RY12*SEPZ - RZ12*SEPY)/DET2
                               YCP  = (RR3*ZX12 + RZ12*SEPX - RX12*SEPZ)/DET2
                               ZCP  = (RR3*XY12 + RX12*SEPY - RY12*SEPX)/DET2
                               !$$--- check if the vertex and I0 are on the same of the
                               !$$    of current moving face
                               if(NMX*XCP + NMY*YCP + NMZ*ZCP .GT. NMD) cycle

                               XCP  = XCP - XC
                               YCP  = YCP - YC
                               ZCP  = ZCP - ZC
                               R2T  = XCP*XCP + YCP*YCP + ZCP*ZCP
                               if(R2T .LT. R2MI) then
                                  R2MI  = R2T
                                  I4 = IW
                                  X4 = SEPX
                                  Y4 = SEPY
                                  Z4 = SEPZ
                                  X3 = XCP + XC
                                  Y3 = YCP + YC
                                  Z3 = ZCP + ZC
                               end if
                           END DO !END LOOP FOR IW

                           !$$--- check if the number of vertice to be larger then permitted number
                           NVERT = NVERT + 1
                           if(NVERT .GT. MXNV2) then
                              STATU  = mp_VSTAT_OPENED
                              exit
                           end if

                           !$$--- add the new vertex to the vertex list of face LF
                           IVERT   = IVERT + 1
                           VERT(IVERT) = I4

                           !$$--- check if the face has been closed
                           if(I4 .LT. 0 ) exit
                           if(I4 .EQ. I2) exit

                           !$$--- add I4 to face list if I4 is still not on the list
                            FLAG(I4) = 1
                            DO J=1, MXNF
                               if(FACE(J).EQ. I4) then
                                  exit
                               else
                                if(FACE(J).LE.0) then
                                   FACE(J) = I4
                                   !$$--- set the first vertex on the face I4
                                   IW = IVERT0+(J-1)*MXNV
                                   VERT(IW+1) = I4
                                   VERT(IW+2) = I3
                                   VERT(IW+3) = I1
                                  exit
                                end if
                               end if
                            END DO

                            !$$--- prepare for the next vertex
                            !$$--- store the positon of current vertex
                            NMD = SQRT(X2*X2+Y2*Y2+Z2*Z2)
                            NMX = X2/NMD
                            NMY = Y2/NMD
                            NMZ = Z2/NMD
                            I3  = I4
                            XC  = X3
                            YC  = Y3
                            ZC  = Z3
                            X2  = X4
                            Y2  = Y4
                            Z2  = Z4
                            RR2 = X2*X2 + Y2*Y2 + Z2*Z2
                            XY12 = X1*Y2 - X2*Y1
                            YZ12 = Y1*Z2 - Y2*Z1
                            ZX12 = Z1*X2 - Z2*X1
                            RX12 = RR1*X2 - RR2*X1
                            RY12 = RR1*Y2 - RR2*Y1
                            RZ12 = RR1*Z2 - RR2*Z1
                        END DO !END LOOP FOR IW1
                        if(STATU .EQ. mp_VSTAT_OPENED) exit
                     END DO !END LOOP FOR LF
                     STAT(IA+IA0) = STATU

            END DO  !END LOOP FOR IA
        RETURN
  END SUBROUTINE Cal_Delaunay_Vertice
  !*********************************************************************************

  !*********************************************************************************
   SUBROUTINE Cal_Voronoi_Vertice(IM, XP, IA0, NPRT, MXNF, KVOIS,INDI, MXNV, DVERT, STAT, VOL, VVERT)
  !***  PURPOSE:   to calculate the Voronoi volume of atoms. Called after calling
  !$$              Cal_Delaunay_Vertice
  !$$
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              XP:     the positions of atoms
  !$$              IA0:    the index (in the whole box)of the fisrt atom on in this call
  !$$              NPRT:   the actual number of atoms on the in this call
  !$$              MXNF:   the max permitted number of neighbores for each atom
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$              DVERT:  the neighbore ID of atoms on Delaunay vertices
  !$$              STAT:   the stat indicating if the Voronoi volume are enclosed
  !$$
  !$$   OUTPUT:    VOL,    Voronoi volume of atoms
  !$$              VVERT,  Voronoi vertice of atoms
  !$$
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, intent(in)::IM, IA0, NPRT, MXNF, MXNV
      real(KINDDF), dimension(:,:), intent(in)::XP
      integer,      dimension(:),   intent(in)::KVOIS
      integer,      dimension(:,:), intent(in)::INDI
      integer,      dimension(:)              ::DVERT, FACE
      integer(1),   dimension(:)              ::STAT

      real(KINDSF), dimension(:,:)::VVERT
      real(KINDSF), dimension(:)::VOL
      !Local variables
      real(KINDSF)::RR1, RR2, RR3, R2T, R2MI, VP1X, VP1Y, VP1Z, VP2X, VP2Y, VP2Z, VP3X, VP3Y, VP3Z
      real(KINDSF)::XY12, YZ12, ZX12, RX12, RY12, RZ12, DET2
      real(KINDSF)::X0, Y0, Z0
      real(KINDSF)::X1, Y1, Z1
      real(KINDSF)::X2, Y2, Z2
      real(KINDSF)::X3, Y3, Z3
      real(KINDSF)::SEP(mp_MXKVOIS,3)
      real(KINDSF)::VVV

      real(KINDSF)::HBX, HBY, HBZ, BX, BY, BZ
      integer::IFPDX, IFPDY, IFPDZ
      integer::NT, NB

      ! Local variables
      integer::IA, MXNV4
      integer::J, IW, IW1, IIW, IFACE0, IVERT0, IVERT, NVERT,I1, I2, I3, I4, LF
      real(KINDSF), parameter::scal=1.E8, underflow=1.E-8, rescal=0.5E-8
      !real(KINDSF), parameter::p_VFACT = 1.D0/48.D0    ! used without scal
       real(KINDSF), parameter::p_VFACT = 1.D-24/48.D0   ! used with scal

            BX = m_BOXSIZE(1)*scal
            BY = m_BOXSIZE(2)*scal
            BZ = m_BOXSIZE(3)*scal

            HBX = BX*C_HALF
            HBY = BY*C_HALF
            HBZ = BZ*C_HALF

            IFPDX = m_IFPD(1)
            IFPDY = m_IFPD(2)
            IFPDZ = m_IFPD(3)
            MXNV4 = MXNV-4

            VVERT = 1.E32
            DO IA=1, NPRT
               VVV = 0.D0
               IVERT0 = (IA-1)*MXNF*MXNV

               if(STAT(IA+IA0) .EQ. mp_VSTAT_ENCLOSED) then
                !$$--- get the position of the atom IA
                 X0 = XP(IA+IA0,1)*scal
                 Y0 = XP(IA+IA0,2)*scal
                 Z0 = XP(IA+IA0,3)*scal
                 IIW = KVOIS(IA+IA0)
                 DO IW =1, min(IIW, mp_MXKVOIS)
                    J=INDI(IA+IA0,IW)
                    SEP(IW,1) = XP(J,1)*scal - X0
                    SEP(IW,2) = XP(J,2)*scal - Y0
                    SEP(IW,3) = XP(J,3)*scal - Z0
                    if(IFPDX.GT.0) then
                      if(ABS(SEP(IW,1)) .GT. HBX) SEP(IW,1) = SEP(IW,1) - SIGN(BX,SEP(IW,1))
                    end if

                    if(IFPDY.GT.0) then
                      if(ABS(SEP(IW,2)) .GT. HBY) SEP(IW,2) = SEP(IW,2) - SIGN(BY,SEP(IW,2))
                    end if

                   if(IFPDZ.GT.0) then
                      if(ABS(SEP(IW,3)) .GT. HBZ) SEP(IW,3) = SEP(IW,3) - SIGN(BZ,SEP(IW,3))
                   end if
                   END DO

                  DO LF =1, MXNF
                     !$$--- get the position of the first vertex on face LF
                      IVERT = IVERT0 + (LF-1)*MXNV
                      I1 = DVERT(IVERT + 1)
                      if(I1 .LE. 0) exit

                      I2 = DVERT(IVERT + 2)
                      I3 = DVERT(IVERT + 3)

                      X1 = SEP(I1,1)
                      Y1 = SEP(I1,2)
                      Z1 = SEP(I1,3)
                      X2 = SEP(I2,1)
                      Y2 = SEP(I2,2)
                      Z2 = SEP(I2,3)
                      X3 = SEP(I3,1)
                      Y3 = SEP(I3,2)
                      Z3 = SEP(I3,3)

                      RR1  = (X1*X1 + Y1*Y1 + Z1*Z1)
                      RR2  = (X2*X2 + Y2*Y2 + Z2*Z2)
                      XY12 = X1*Y2 - X2*Y1
                      YZ12 = Y1*Z2 - Y2*Z1
                      ZX12 = Z1*X2 - Z2*X1
                      RX12 = RR1*X2 - RR2*X1
                      RY12 = RR1*Y2 - RR2*Y1
                      RZ12 = RR1*Z2 - RR2*Z1

                      RR3  = (X3*X3 + Y3*Y3 + Z3*Z3)
                      DET2 = X3*YZ12 + Y3*ZX12 + Z3*XY12
                      VP1X  = (RR3*YZ12 + RY12*Z3 - RZ12*Y3)/DET2
                      VP1Y  = (RR3*ZX12 + RZ12*X3 - RX12*Z3)/DET2
                      VP1Z  = (RR3*XY12 + RX12*Y3 - RY12*X3)/DET2

                     !$$--- save the position
                      VVERT(IVERT+1,1) = VP1X*rescal
                      VVERT(IVERT+1,2) = VP1Y*rescal
                      VVERT(IVERT+1,3) = VP1Z*rescal
                     !$$--- get the position of the second vertex on face LF
                      X2  = X3
                      Y2  = Y3
                      Z2  = Z3
                      RR2 = RR3
                      XY12 = X1*Y2 - X2*Y1
                      YZ12 = Y1*Z2 - Y2*Z1
                      ZX12 = Z1*X2 - Z2*X1
                      RX12 = RR1*X2 - RR2*X1
                      RY12 = RR1*Y2 - RR2*Y1
                      RZ12 = RR1*Z2 - RR2*Z1

                      IVERT= IVERT+4
                      I3 = DVERT(IVERT)
                      X3 = SEP(I3,1)
                      Y3 = SEP(I3,2)
                      Z3 = SEP(I3,3)

                      RR3  = X3*X3 + Y3*Y3 + Z3*Z3
                      DET2 = X3*YZ12 + Y3*ZX12 + Z3*XY12
                      VP2X  = (RR3*YZ12 + RY12*Z3 - RZ12*Y3)/DET2
                      VP2Y  = (RR3*ZX12 + RZ12*X3 - RX12*Z3)/DET2
                      VP2Z  = (RR3*XY12 + RX12*Y3 - RY12*X3)/DET2

                     !$$--- save the position
                      VVERT(IVERT-2,1) = VP2X*rescal
                      VVERT(IVERT-2,2) = VP2Y*rescal
                      VVERT(IVERT-2,3) = VP2Z*rescal

                     !$$--- prepare for next vertex
                      X2  = X3
                      Y2  = Y3
                      Z2  = Z3
                      RR2 = RR3
                      XY12 = X1*Y2 - X2*Y1
                      YZ12 = Y1*Z2 - Y2*Z1
                      ZX12 = Z1*X2 - Z2*X1
                      RX12 = RR1*X2 - RR2*X1
                      RY12 = RR1*Y2 - RR2*Y1
                      RZ12 = RR1*Z2 - RR2*Z1

                     !$$---  scan the next neighbores to construct other vertex
                       DO IW=1, MXNV4
                          IVERT= IVERT+1
                          I3 = DVERT(IVERT)
                          if(I3 .LE. 0) then
                             VVERT(IVERT-2,1) = 2.E31
                             VVERT(IVERT-2,2) = 2.E31
                             VVERT(IVERT-2,3) = 2.E31
                             exit
                          end if
                          X3 = SEP(I3,1)
                          Y3 = SEP(I3,2)
                          Z3 = SEP(I3,3)

                          RR3  = X3*X3 + Y3*Y3 + Z3*Z3
                          DET2 = X3*YZ12 + Y3*ZX12 + Z3*XY12
                          VP3X  = (RR3*YZ12 + RY12*Z3 - RZ12*Y3)/DET2
                          VP3Y  = (RR3*ZX12 + RZ12*X3 - RX12*Z3)/DET2
                          VP3Z  = (RR3*XY12 + RX12*Y3 - RY12*X3)/DET2

                         !$$--- accumulating volume of tetrahedron
                          VVV = VVV + ABS( (VP1X*(VP2Y*VP3Z - VP2Z*VP3Y) + &
                                            VP1Y*(VP2Z*VP3X - VP2X*VP3Z) + &
                                            VP1Z*(VP2X*VP3Y - VP2Y*VP3X) ) )

                         !$$--- save the position
                          VVERT(IVERT-2,1) = VP3X*rescal
                          VVERT(IVERT-2,2) = VP3Y*rescal
                          VVERT(IVERT-2,3) = VP3Z*rescal

                         !$$--- prepare for next
                          VP2X = VP3X
                          VP2Y = VP3Y
                          VP2Z = VP3Z
                          X2   = X3
                          Y2   = Y3
                          Z2   = Z3
                          RR2  = RR3
                          XY12 = X1*Y2 - X2*Y1
                          YZ12 = Y1*Z2 - Y2*Z1
                          ZX12 = Z1*X2 - Z2*X1
                          RX12 = RR1*X2 - RR2*X1
                          RY12 = RR1*Y2 - RR2*Y1
                          RZ12 = RR1*Z2 - RR2*Z1
                      END DO
                  END DO !END LOOP FOR LF
                  VOL(IA+IA0)  = VVV*p_VFACT
                else
                  VOL(IA+IA0)  = -1.0
                end if
            END DO
        RETURN
  END SUBROUTINE Cal_Voronoi_Vertice
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Cal_VoronoiVolume(SimBox, CtrlParam, List, hDVTS, hVVTS, hVOLS)
  !***  PURPOSE:  to calculate Voronoi tesslation
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT
  !
   use MD_NeighborsList
   implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(inout)::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(NEIGHBOR_LIST),intent(in)::List
      integer,      dimension(:,:,:),   optional::hDVTS
      real(KINDSF), dimension(:,:,:,:), optional::hVVTS
      real(KINDSF), dimension(:),       optional::hVOLS

      !--- Local variables
      integer::IA0, NPRT, I, J, K, IP0, IN, NN

          if(m_INITED .eq. 0) then
             call Initialize_VoronoiVolume(SimBox, CtrlParam)
          end if

          !$$--- If pressure effect acounted for, the box size could be changed, we
          !$$    need to recopy the boxsize
             m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
             m_BOXSIZE(1:3)      = SimBox%ZL(1:3)
             m_IFPD = CtrlParam%IFPD
             IA0 = 0
             NPRT = mp_MXNA
             DO WHILE(.TRUE.)
                if(IA0+NPRT.GE.SimBox%NPRT) then
                   NPRT = SimBox%NPRT - IA0
                end if
                call Cal_Delaunay_Vertice(SimBox%NPRT, SimBox%XP, IA0, NPRT, m_MXNF, List%KVOIS, List%INDI, m_MXNV, m_DVTS, m_VSTAT, m_WVFS, m_WFLG)
                call Cal_Voronoi_Vertice(SimBox%NPRT, SimBox%XP, IA0, NPRT, m_MXNF, List%KVOIS, List%INDI, m_MXNV, m_DVTS, m_VSTAT, m_VOLS, m_VVTS)

                !$$--- ouput current tessellation if I/O unit are provided
                if(m_OUTPUTDT  > 0) call Ouput_DelaunayVertice(IA0, NPRT, LIST, SimBox%RR, m_DVTS)
                if(m_OUTPUTVT  > 0) call Ouput_VoronoiVertice(IA0, NPRT, LIST, SimBox%RR, m_DVTS, m_VVTS)

                if(present(hDVTS)) then
                   IP0= 1
                   DO I=1, NPRT
                      DO J=1, m_MXNF
                        DO K=1, m_MXNV
                           IN = m_DVTS(IP0)
                           if(IN .GT. 0) then
                              NN = List%INDI(I+IA0, IN)
                              hDVTS(I+IA0,J,K) = NN
                           else
                             hDVTS(I+IA0,J,K) =0
                           end if
                           IP0 = IP0 + 1
                        END DO
                      END DO
                   END DO
                end if

                if(present(hVVTS)) then
                   IP0= 1
                   DO I=1, NPRT
                      DO J=1, m_MXNF
                        DO K=1, m_MXNV
                           hVVTS(I+IA0,J,K,1:3) = m_VVTS(IP0,1:3)
                           IP0 = IP0 + 1
                        END DO
                      END DO
                   END DO
                end if
                IA0 = IA0 + NPRT
                if(IA0.GE.SimBox%NPRT) exit
             END DO

          return
  end SUBROUTINE Cal_VoronoiVolume
  !*********************************************************************************

   !**********************************************************************************
  subroutine Ouput_DelaunayVertice(IA0, NPRT, LIST, LATT, hDVTS)
  !***  PURPOSE:  to output Delaunay vertices get in Cal_Delaunay_Vertice
  !               for atoms from IA0+1 to IA0+NPRT to a file uint. T
  !
  !     INPUT:    IA0:  the starting atom
  !               NPRT: the number of atoms
  !               LIST: the neighbor-list
  !               LATT: length unit
  !               hDVTS: Delaunay vertices
  !     OUTPUT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
       !--- dymmy varibales
       integer,               intent(in)::IA0, NPRT
       type(NEIGHBOR_LIST),   intent(in)::List
       real(KINDDF),          intent(in)::LATT
       integer, dimension(:), intent(in)::hDVTS
       !--- local variables
       integer, dimension(:), allocatable::NVERTS
       integer::NFACES, I, J, K, IA, IP0, IN, IP

       !---
            allocate(NVERTS(m_MXNF))
            DO I=1, NPRT
               NFACES = 0
               NVERTS = 0
               IP0 = (I-1)*m_MXNF*m_MXNV
               !get the number of faces and vertices
               DO J=1, m_MXNF
                  IP = IP0+(J-1)*m_MXNV+1
                  IN = hDVTS(IP)
                  if(IN.LE.0) exit
                  NFACES = NFACES + 1
                  DO K=1, m_MXNV
                     NVERTS(J) = NVERTS(J) + 1
                     IP = IP + 1
                     IN = hDVTS(IP)
                     if(IN.LE.0) exit
                  END DO
               END DO
               write(m_OUTPUTDT, *) "****************************************************************"
               write(m_OUTPUTDT, fmt="(3(I8, 1x), 5(1PE15.6,2x))") I+IA0, I+IA0, m_ITYP(I+IA0),m_XP(I+IA0,1:3)/LATT
               write(m_OUTPUTDT, fmt="(I8, 1x, 5(1PE15.6,2x))") NFACES
               DO J=1, NFACES
                  write(m_OUTPUTDT, fmt="(I8, 1x, 5(1PE15.6,2x))") NVERTS(J)
                  IP = IP0+(J-1)*m_MXNV
                  DO K=1,NVERTS(J)
                     IP  = IP + 1
                     IN  = hDVTS(IP)
                     IA  = List%INDI(I+IA0,IN)
                     write(m_OUTPUTDT, fmt="(4(I8, 1x), 5(1PE15.6,2x))") &
                              K, IA, IA, m_ITYP(IA), (m_XP(IA,1:3)- m_XP(I+IA0,1:3))/LATT
                  END DO
               END DO

            END DO
           deallocate(NVERTS)
           return
  end subroutine Ouput_DelaunayVertice
  !***********************************************************************************

  !**********************************************************************************
  subroutine Ouput_VoronoiVertice(IA0, NPRT, LIST, LATT, hDVTS, hVVTS)
  !***  PURPOSE:  to output Voronoi vertices get in Cal_Voronoi_Volume
  !               for atoms from IA0+1 to IA0+NPRT to a file uint.
  !
  !
  !     INPUT:    IA0:  the starting atom
  !               NPRT: the number of atoms
  !               LIST: the neighbor-list
  !               LATT: length unit
  !               hVVTS: Voronoi vertices
  !
  !     OUTPUT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
       !--- dymmy varibales
       integer, intent(in)::IA0, NPRT
       type(NEIGHBOR_LIST),          intent(in)::List
       real(KINDDF),                 intent(in)::LATT
       real(KINDSF), dimension(:,:), intent(in)::hVVTS
       integer,      dimension(:),   intent(in)::hDVTS

       !--- local variables
       integer, dimension(:), allocatable::NVERTS
       integer::NFACES, I, J, K, INDEX, IP0, IP, IA, IN

       !---
           allocate(NVERTS(m_MXNF))
            DO I=1, 1!NPRT
               NFACES = 0
               NVERTS = 0
               IP0 = (I-1)*m_MXNF*m_MXNV
               !get the number of faces and vertices
               DO J=1, m_MXNF
                  IP = IP0+(J-1)*m_MXNV+1
                  if(ANY(hVVTS(IP,1:3) .GE. 1.E31)) then
                     exit
                  end if
                  NFACES = NFACES + 1
                  DO K=1, m_MXNV
                     NVERTS(J) = NVERTS(J) + 1
                     IP = IP + 1
                     if(ANY(hVVTS(IP,1:3) .GE. 1.E31)) exit
                  END DO
               END DO

               write(m_OUTPUTVT, *) "****************************************************************"
               write(m_OUTPUTVT, fmt="(3(I8, 1x), 5(1PE15.6,2x))") IA0+I, IA0+I, m_ITYP(IA0+I), m_XP(IA0+I,1:3)/LATT
               write(m_OUTPUTVT, fmt="(I8, 1x, 5(1PE15.6,2x))") NFACES
               DO J=1, NFACES
                  ! we also output the normal of the face
                  IP  = IP0+(J-1)*m_MXNV
                  IN  = hDVTS(IP+1)
                  IA  = List%INDI(IA0+I,IN)
                  write(m_OUTPUTVT, fmt="(2(I8, 1x), 5(1PE15.6,2x))") NVERTS(J),  m_ITYP(IA), (m_XP(IA,1:3)-m_XP(IA0+I,1:3))/LATT
                  DO K=1,NVERTS(J)
                     IP  = IP + 1
                     write(m_OUTPUTVT, fmt="(I5, 1x, I5, 1x, 5(1PE15.6,2x))") J, K, hVVTS(IP,1:3)/LATT
                  END DO
               END DO
            END DO

           deallocate(NVERTS)
           return
  end subroutine Ouput_VoronoiVertice
  !****************************************************************************************

 !**********************************************************************************
  subroutine Output_VoronoiVol(hFile, List, LATT, VOLS, STATU)
  !***  DESCRIPTION: to putout the Delaunay tesselation results to output file.
  !
  !     INPUT:  hFile,  unit of the output file
  !
  !
  use MD_Globle_Variables_GPU
  use MD_NeighborsList
   implicit none
       integer, intent(in)::hFile
       type(NEIGHBOR_LIST),intent(in)::List
       real(KINDDF),intent(in)::LATT
       real(KINDSF), dimension(:)::VOLS
       integer(1), dimension(:), allocatable::STATU

       !--- local variables
       integer::I
       real(KINDSF)::TOTALV, UNIT

       !---
            TOTALV = 0.0
            UNIT   = LATT**3.D0
            DO I=1, dm_NPRT
               if(STATU(I) .EQ. mp_VSTAT_ENCLOSED) then
                  TOTALV = TOTALV+VOLS(I)/UNIT
               end if
            END DO
            print *, "TOTALV", TOTALV

            DO I=1, dm_NPRT
               if(STATU(I) .EQ. mp_VSTAT_ENCLOSED) then
                  write(hFile, fmt="(I8, 1x, 6(1PE15.6,2x))") m_ITYP(I), m_XP(I,1:3)/LATT, VOLS(I)/UNIT, UNIT/VOLS(I)
               else if(STATU(I) .EQ. mp_VSTAT_OPENED) then
                  write(hFile, fmt="(I8, 1x, 6(1PE15.6,2x))") m_ITYP(I)+dm_NGROUP, m_XP(I,1:3)/LATT, VOLS(I)/UNIT, UNIT/VOLS(I)
               end if
            END DO


           return
  end subroutine Output_VoronoiVol
  !**********************************************************************************

  !**********************************************************************************
  subroutine RECORD_VoronoiTessellation(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results of Voronoi tesselation to output file.
  !                  This routine is to interfaced to and MD_EM_TB_ForceTable_SHELL_12_GPU.
  !                  It is assumed the the neighbor-ist routine has be called
  !
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
   use MD_NeighborsList
   use MD_NeighborsList_GPU
   implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE, FMT
       equivalence (GFILE, FMT)
       type(NEIGHBOR_LIST)::List
       integer::OUTPUTDT, OUTPUTVT, hFILE

       character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
       integer::DATE_TIME1(8),DATE_TIME2(8)
       real*4::C1,C2,TIME1, TIME2
       integer::LOOP=1, I, JOB, ICFG
       !----
          if(m_INITED .eq. 0) then
             call Initialize_VoronoiVolume(SimBox(1), CtrlParam)
          end if
          JOB   = Stamp%ITest
          ICFG  = Stamp%IRec(1)

           !$$--- prepare for output
           OUTPUTDT   = m_OUTPUTDT
           OUTPUTVT   = m_OUTPUTVT
           if(m_VT_OFILE > 0 .and. (OUTPUTDT>0 .OR. OUTPUTVT>0)) then
              if(OUTPUTDT>0) then
                 !$$--- output for Delaunay Tesslation
                  call STRCATI(GFILE, CtrlParam%f_others(m_VT_OFILE), "_D_P", m_processid, 4)
                  call STRCATI(GFILE, GFILE, "_", JOB, 4)
                  call STRCATI(GFILE, GFILE, ".", ICFG, 4)
                  call AvailableIOUnit(m_OUTPUTDT)
                  open(UNIT=m_OUTPUTDT, file = GFILE, status='unknown')
                  write(*,*) "Delaunay vertices to be saved in "//GFILE(1:len_trim(GFILE))

               end if

              if(OUTPUTVT>0) then
               !$$--- output Voronoi Tesslation
                 call STRCATI(GFILE, CtrlParam%f_others(m_VT_OFILE), "_V_P", m_processid, 4)
                 call STRCATI(GFILE, GFILE, "_", JOB, 4)
                 call STRCATI(GFILE, GFILE, ".", ICFG, 4)
                 call AvailableIOUnit(m_OUTPUTVT)
                 open(UNIT=m_OUTPUTVT, file = GFILE, status='unknown')
                 write(*,*) "Voronoi vertices to be saved in "//GFILE(1:len_trim(GFILE))
              end if
           end if

           !$$--- in cpu version this neighbore-list is also calculated on device
           !$$     need to copy it out
           call Copyout_NeighboreList_DEV(List, ORDER=1)

           CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
           C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000
           DO I=1, LOOP
              call Cal_VoronoiVolume(SimBox(1), CtrlParam, List)
           END DO
           CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
           C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
           TIME2 = C2-C1
           write(*,*) "The elapsed time for Voronoi Tessellation(ms):", TIME2

           !restore the output control
            if(m_OUTPUTDT>0 .OR. m_OUTPUTVT>0) then
              call Clear_NeighboreList(List)
              if(m_OUTPUTDT .GT. 0)  close(m_OUTPUTDT)
              if(m_OUTPUTVT .GT. 0)  close(m_OUTPUTVT)
            end if
            m_OUTPUTDT  = OUTPUTDT
            m_OUTPUTVT  = OUTPUTVT


           if(m_OUTPUTVOL >0) then
             !$$--- output atomic volume
              call STRCATI(GFILE, CtrlParam%f_others(m_VT_OFILE), "_VOL_P", m_processid, 4)
              call STRCATI(GFILE, GFILE, "_", JOB, 4)
              call STRCATI(GFILE, GFILE, ".", ICFG, 4)

              call AvailableIOUnit(hFile)
              write(*,*) "Save Atomic volume  to "//GFILE(1:len_trim(GFILE))
              open(UNIT=hFile, file = GFILE, status='unknown')
              call Output_VoronoiVol(hFile, List, SimBox(1)%RR, m_VOLS, m_VSTAT)
              close(hFile)
           end if

           return
  end subroutine RECORD_VoronoiTessellation
  !****************************************************************************************

 !****************************************************************************************
  subroutine RECORD_VoronoiTessellation_TOOL(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calculate and output the  the results of Voronoi tesselation to output file.
  !                  This routine is to interfaced MD_SimBoxArray_ToolShell_12_GPU.
  !                  It is assumed the the neighbor-list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
   use MD_NeighborsList_GPU
   implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables

           if(Stamp%ITime .LT. 0) then
              call Initialize_NeighboreList_DEV(SimBox, CtrlParam)
              return
            end if
            call Cal_NeighBoreList_DEV(SimBox, CtrlParam)
            call RECORD_VoronoiTessellation(Stamp, SimBox, CtrlParam)
          return
  end subroutine RECORD_VoronoiTessellation_TOOL
  !****************************************************************************************


  end module VoronoiTessellation_CPU





