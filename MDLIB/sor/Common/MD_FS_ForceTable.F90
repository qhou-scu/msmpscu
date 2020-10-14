  module MD_FS_Force_Table
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                 To calculate the forces on atoms, the force is assumed in tight-binding form of Finnis:
  !                    Ec=(RO**(1/2)+SUMij(Vij)
  !                 with RO the total electron density on atoms; Vij, the pairwise portential on atom i by atom j.
  !                 The force calculations are performed by using FORCETABLE previuosly calculated.
  !                 A FORCETABLE should be register before calling this module.
  !
  !                 REF___________________________________________________________________________________
  !                    (1) M.W.Finnis and J.E.Sinclair Phil.Mag. A50 (1984) 45
  !                    (2) G.J.Ackland and V.Vitek, Phys.Rev.B41, (1990) 10324
  !
  !                 SEE ALSO______________________________________________________________________________
  !                      MDTypeDef_ForceTable.f90
  !                      Cal_TB_EMPIRICAL_ForceTable.f90( version 2000)
  !                 ______________________________________________________________________________________
  !                 HOU Qing, Mar, 2010
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_ForceTable
    use MD_NeighborsList
  !-----------------------------------------------
    implicit none

      logical,private::m_INITED = .false.

      integer,private::m_IFPD(3)                                        ! if minima image invention periodic condition
      real(KINDDF), private::m_HBOX(3)                                  ! half size of box
      real(KINDDF), private::m_BOX(3)                                   ! size of simulation box, could in centimeter, lattice unit and UNIT size
      real(KINDDF), private::m_BOXSHAPE(3,3)                            !the box shape, would be used in Rahmann-Parrilino Scheme

      real(KINDDF), private::m_Rmax                                     ! the maxma value of R in the table
      real(KINDDF), dimension(mp_MXGROUP,mp_MXGROUP),private::m_RCUT2   ! the square of RCUT
      real(KINDDF), private::m_RmaxSQRT                                 ! the square root of Rmax
      real(KINDDF)::m_CSI                                               ! the inverse of m_RmaxSQRT, e.g. NTAB/m_RmaxSQRT
      real(KINDDF)::m_CSIV                                              ! the inverse of CSI         1/CSI
      integer, dimension(:,:), allocatable::m_KPAIR                     ! the index of table for atom1 and atom2
      real(KINDDF),dimension(:,:),allocatable::m_POTR                   ! the potential table for Nucl.-Nucl.  interaction with size of  NKIND*NTAB
      real(KINDDF),dimension(:,:),allocatable::m_FPOTR                  ! the differential of potential table POTR
      real(KINDDF),dimension(:,:),allocatable::m_POTB                   ! the potential table
      real(KINDDF),dimension(:,:),allocatable::m_FPOTB                  ! the differential on r of electron density
      type(NEIGHBOR_LIST), pointer::m_NLISTP                            ! the pointer pointing the neighbore list target,
                                                                        ! could be used in analysis routines

      real(KINDDF), dimension(:), allocatable, private::m_ER            ! the total pairwise potential on atom SUMj(Vij)
      real(KINDDF), dimension(:), allocatable, private::m_DEN           ! the total electron density on atom SUMj(ROij)
      real(KINDDF), dimension(:), allocatable, private::m_DDEN          ! the derivetive of F(RHO), for Finnis-Sinlar potential, it is 1/DSQRT(m_DEN)
      !****
      private::preCALFOR, CALFOR, preCALFOR2, CALFOR2, GET_DIFF_FFUN_FS
  contains

  !*********************************************************************************
    SUBROUTINE INITIALIZE_FS_Force_Table(SimBox, CtrlParam, FTable, RelaseTable)
    !***  PURPOSE:   to intialize the parameters to be used in force calculation
    !    INPUT:     SimBox: the simulation box
    !               CtrlParam: the control parameters
    !               FTable, the force tabled to be used
    !
    !    OUTPUT
    implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(inout)::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(in)   ::FTable
      integer,           optional     ::RelaseTable
      !--- Local variables
      integer::I,J

           !---
              call RELEASE_FS_Force_Table()

           !$$--- to copy simulation box information
               m_BOX =  SimBox%ZL
               m_HBOX = SimBox%ZL*C_HALF
               m_BOXSHAPE = SimBox%BOXSHAPE
               m_IFPD = CtrlParam%IFPD
               allocate(m_ER(SimBox%NPRT),m_DEN(SimBox%NPRT),m_DDEN(SimBox%NPRT))

           !$$--- to copy force table information
               m_RMAX = FTable%RMAX
               m_RMAXSQRT = FTable%RMAXSQRT
               m_CSI = FTable%CSI
               m_CSIV = FTable%CSIV
               if(SIZE(CtrlParam%RU) .NE. SIZE(m_RCUT2)) then
                 print *, SIZE(CtrlParam%RU)
                 print *, SIZE(m_RCUT2)
                 print *, "Size of RU matrix is inconsistent in force calculation. Check the program"
                 stop
               end if
               m_RCUT2 = CtrlParam%RU*CtrlParam%RU

               allocate(m_KPAIR(SIZE(FTable%KPAIR,dim=1),  SIZE(FTable%KPAIR,dim=2)), &
                        m_POTR(SIZE(FTable%POTR,dim=1),    SIZE(FTable%POTR,dim=2)),  &
                        m_POTB(SIZE(FTable%POTB,dim=1),    SIZE(FTable%POTB,dim=2)),  &
                        m_FPOTR(SIZE(FTable%FPOTR,dim=1),  SIZE(FTable%FPOTR,dim=2)), &
                        m_FPOTB(SIZE(FTable%FPOTB,dim=1),  SIZE(FTable%FPOTB,dim=2))  )

               m_KPAIR = FTable%KPAIR
               m_POTR = FTable%POTR
               m_POTB = FTable%POTB
               m_FPOTR = FTable%FPOTR
               m_FPOTB = FTable%FPOTB

               if(present(RelaseTable)) then
                 if(RelaseTable) call Release_ForceTable(FTable, keepreglist=1)
               end if

               m_INITED = .TRUE.
           return
    END SUBROUTINE INITIALIZE_FS_Force_Table

  !*********************************************************************************
    SUBROUTINE RELEASE_FS_Force_Table()
    !***  PURPOSE:   to release memory allocated in this module
    !     INPUT:
    !     OUTPUT
    implicit none
       !--- dummy vaiables


               if(allocated(m_ER)) then
                  deallocate(m_ER)
               end if

               if(allocated(m_DEN)) then
                   deallocate(m_DEN)
               end if

               if(allocated(m_DDEN)) then
                   deallocate(m_DDEN)
               end if

               if(allocated(m_KPAIR)) deallocate(m_KPAIR)
               if(allocated(m_POTR)) deallocate(m_POTR)
               if(allocated(m_POTB)) deallocate(m_POTB)
               if(allocated(m_FPOTR)) deallocate(m_FPOTR)
               if(allocated(m_FPOTB)) deallocate(m_FPOTB)

               m_INITED = .FALSE.
           return
    END SUBROUTINE RELEASE_FS_Force_Table

  !*********************************************************************************
  SUBROUTINE CALFORCE_FS_Force_Table(SimBox, List)
  !***  PURPOSE:   to begine calculate the force, Newton's third law taken into account
  !
  !     INPUT:     SimBox,    the simulation box
  !                List,      the neighbore list
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
   use MD_NeighborsList
   implicit none
      !--- dummy vaiables
      type(SimMDBox),            intent(inout)::SimBox
      type(NEIGHBOR_LIST),target,intent(in)   ::List


      !--- Local variables
               !$$-- keep the pointer of the current neioghbore list
                m_NLISTP =>List

               !$$--- The boxsize could be changed due to pressure
               !$$    we should update it at each time
                m_BOX =  SimBox%ZL
                m_HBOX = SimBox%ZL*C_HALF
                m_BOXSHAPE = SimBox%BOXSHAPE

                !$$--- To calculate the electron densitis
                call PRECALFOR(SimBox%NPRT, List%KVOIS, List%INDI, SimBox%ITYP, SimBox%XP)

                !$$--- To get the force
                call CALFOR(SimBox%NPRT, List%KVOIS,List%INDI, SimBox%ITYP, SimBox%XP, SimBox%FP, SimBox%VTENSOR)

                !$$--- at the same time we update potential energy of the atoms
                 SimBox%EPOT = m_ER - DSQRT(m_DEN)
                 return

  end SUBROUTINE CALFORCE_FS_Force_Table
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE preCALFOR(IM,KVOIS,INDI,ITYP,XP)
  !***  PURPOSE:   to begine calculate the pairwise potential m_ER and the electron density m_DEN
  !     INPUT:     IM:     the number of particles concerned
  !                KVOIS:  the number of for each atom in style1
  !                INDI:   the index for neighbores
  !                ITYPE:  the type of atom corresponding to INDI
  !     OUTPUT     m_ER
  !                m_DEN

  !
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, intent(in)::IM
      integer, dimension(:), intent(in)::KVOIS
      integer, dimension(:,:), intent(in)::INDI
      integer, dimension(:), intent(in)::ITYP
      real(KINDDF), dimension(:,:), intent(in)::XP

      !---Local variables
      integer::I,J,K,KK, IW, IIW, KTAB,KTAB1, K1, ITYPI, ITYPJ
      real(KINDDF)::SK,DK,EXP1,EXP2
      real(KINDDF)::SEP(3), DXYZ(3), R2, R


           m_ER = C_ZERO
           m_DEN = C_ZERO

           DO I=1, IM
              IIW = KVOIS(I)
              ITYPI  = ITYP(I)
              DO IW=1, IIW
                 J=INDI(I,IW)

                 !$$--- To calculate the seperation between particle I and its IWth neighbore
                 DO K=1, 3
                    SEP(K) = XP(I,K)-XP(J,K)
                    IF((m_IFPD(K).GT.0) .AND. (DABS(SEP(K)) .GT. m_HBOX(K))) THEN
                       SEP(K) = SEP(K) - DSIGN(m_BOX(K),SEP(K))
                    END IF
                 END DO

                 !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                 !$$          only when the shape of the box is cubic,
                 !$$          DXYZ = SEP
                 !
                 DO K=1, 3
                    DXYZ(K) =  sum(m_BOXSHAPE(K,1:3)*SEP(1:3))
                 END DO
                 R2  = sum(DXYZ*DXYZ)

                 !$$--- To calculate electron density on atom I
                 ITYPJ=ITYP(J)
                 if(R2 .LE. m_RCUT2(ITYPI,ITYPJ)) then
                    KTAB = m_KPAIR(ITYPI,ITYPJ)
                    KTAB1= m_KPAIR(ITYPJ,ITYPI)

                    !$$ SK=R2*m_CSI
                    !$$ SK = (C_TWO*m_Rmax*DSQRT(R2)-R2)*m_CSI
                    !$$--- Change made 2010-12-05
                    R = DSQRT(R2)
                    SK= DSQRT(R)*m_CSI
                    KK=SK
                    K1=KK+1
                    DK=SK-DBLE(KK)
                    !$$--- NOTE by HOU Qing: Dec 4,2014
                    !$$          The force table POTR, POTB are V(r), RHO(r)*0.5, in OLD VERSION
                    !$$          In the new version. we define POTR, POTB as V(r)*r and RHO(r)
                    !EXP1 = m_POTR(KTAB,KK)+DK*(m_POTR(KTAB,K1)-m_POTR(KTAB,KK))
                    EXP1 = (m_POTR(KTAB,KK)+DK*(m_POTR(KTAB,K1)-m_POTR(KTAB,KK)) )/R
                    m_ER(I) = m_ER(I)+EXP1
                    m_ER(J) = m_ER(J)+EXP1

                    m_DEN(I)= m_DEN(I) + (m_POTB(KTAB,KK) + DK*(m_POTB(KTAB,K1)-m_POTB(KTAB,KK))   )
                    m_DEN(J)= m_DEN(J) + (m_POTB(KTAB1,KK)+ DK*(m_POTB(KTAB1,K1)-m_POTB(KTAB1,KK)) )

                end if
            ENDDO
         END DO

         !$$*** to get the derivetive of DSQRT(m_DEN)
          DO I=1, IM
             call GET_DIFF_FFUN_FS(m_DEN(I), m_DDEN(I))
         END DO
        RETURN
  END SUBROUTINE preCALFOR
  !*************************************************************************************

  !*************************************************************************************
  SUBROUTINE CALFOR(IM, KVOIS,INDI,ITYP,XP, FP,VTENSOR)
  !***  PURPOSE:   to begine calculate the force, Newton's third law taken into account
  !     INPUT:     IM:     the number of particles
  !                KVOIS:  the number of neighbores in stle0
  !                INDI:   the index for neighbores
  !                ITYP:   the type of atom corresponding to INDI
  !                XP:     the coordinate of atoms
  !     OUTPUT     FP:     the force on atoms
  !                VTENSOR: the virial tensor
  !                EPOT:   the potential of a atom
  !
  use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, intent(in)::IM
      integer, dimension(:),intent(in)::KVOIS
      integer, dimension(:,:), intent(in)::INDI
      integer, dimension(:), intent(in)::ITYP
      real(KINDDF), dimension(:,:),  intent(in)::XP
      real(KINDDF), dimension(:,:), intent(out)::FP
      real(KINDDF), intent(out)::VTENSOR(3,3)

         !---Local variables
      integer::I,J,K,K1,N,IW,IIW,KTAB,KTAB1,KK, KK1,ITYPI, ITYPJ
      real(KINDDF)::SK,DK, EXP5,EXP6,EXP7, DENKI, DENKJ, FORTOT, FOR(3)
      real(KINDDF)::SEP(3), DXYZ(3), R2, R


        !$$--- now, we have the electron densities on atoms, to calculate the force
          FP = C_ZERO
          VTENSOR = C_ZERO
          DO I=1, IM
             IIW = KVOIS(I)
             ITYPI  = ITYP(I)
             !$$IF(m_DEN(I) .GT. 0) then
             !$$   DENKI = C_UN/DSQRT(m_DEN(I))
             !$$ELSE
             !$$  DENKI = C_ZERO
             !$$END IF
             DENKI = m_DDEN(I)

             !$$--- to begine scanning the neighbores
             DO IW=1, IIW
                !$$--- To calculate the seperation between particle I and its IWth neighbore
                J=INDI(I,IW)
                DO K=1, 3
                   SEP(K) = XP(I,K)-XP(J,K)
                   IF(m_IFPD(K) .AND. DABS(SEP(K)) .GT. m_HBOX(K)) THEN
                      SEP(K) = SEP(K) - DSIGN(m_BOX(K),SEP(K))
                   END IF
                END DO

                !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                !$$          only when the shape of the box is cubic,
                !$$          DXYZ = SEP
                DO K=1, 3
                   DXYZ(K) =  sum(m_BOXSHAPE(K,1:3)*SEP(1:3))
                END DO
                R2  = sum(DXYZ*DXYZ)

                !$$--- IF R2 is smaller than cutoff range, to calculate the force
                ITYPJ=ITYP(J)
                if(R2 .LE. m_RCUT2(ITYPI,ITYPJ)) then
                   KTAB  = m_KPAIR(ITYPI,ITYPJ)
                   KTAB1 = m_KPAIR(ITYPJ,ITYPI)
                   !$$SK = R2*m_CSI
                   !$$SK = (C_TWO*m_Rmax*DSQRT(R2)-R2)*m_CSI
                   !$$--- Change made 2010-12-05
                   R = DSQRT(R2)
                   SK= DSQRT(R)*m_CSI
                   KK = SK
                   KK1 = KK+1
                   DK=SK-dble(KK)
                   !$$-- interpolation of force table
                   !$$-- for  pairwise section
                   !EXP5 = m_FPOTR(KTAB,KK)  + DK*(m_FPOTR(KTAB,KK1) - m_FPOTR(KTAB,KK))

                   !$$-- for many-body section
                   !$$    Note; the electron density of atom I on atom J could not be the same as
                   !$$    from atom J on atom J, if atom I and J are not the same type.
                   !
                   !$$-- from the type of atom I to the type of atom J
                   !EXP6 = m_FPOTB(KTAB,KK)  + DK*(m_FPOTB(KTAB,KK1) - m_FPOTB(KTAB,KK))

                   !$$-- from the type of atom J to the type of atom I
                   !EXP7 = m_FPOTB(KTAB1,KK) + DK*(m_FPOTB(KTAB1,KK1) - m_FPOTB(KTAB1,KK))

                   !$$FORTOT=EXP5-EXP6*(DENKI+DENKJ)

                   !$$--- NOTE by Hou Qing: Dec 4,2014
                   !$$          The force table FPOTR, FPOTB are dV(r)/dr/r, and dRHO/dr/r, in OLD VERSION
                   !$$          In the new version. we define FPOTR, FPOTB as (dV(r)/dr)*r and (dRHO/dr)*r
                   FORTOT= (m_FPOTR(KTAB,KK)    + DK*(m_FPOTR(KTAB,KK1)  - m_FPOTR(KTAB,KK)))/R2        +   &
                           ( (m_FPOTB(KTAB1,KK) + DK*(m_FPOTB(KTAB1,KK1) - m_FPOTB(KTAB1,KK)))*m_DDEN(J)+   &
                             (m_FPOTB(KTAB,KK)  + DK*(m_FPOTB(KTAB,KK1)  - m_FPOTB(KTAB,KK)))*DENKI )/R

                   do K=1, 3
                      FOR(K) = FORTOT*SEP(K)
                      FP(I,K) = FP(I,K)+FOR(K)
                      FP(J,K) = FP(J,K)-FOR(K)
                   end do

                   !$$- at the same time, we can calculate the pressure tensor.
                   !$$   NOTE: here we have consider periodic condition, thus
                   !$$         we cannot laterly use PRESS=SUM(XP*FP) to get the pressure.
                   !$$         It would be convenience and efficient to have the pressure
                   !$$         calculated here.
                   DO K=1, 3
                      DO K1=1, 3
                         VTENSOR(K,K1) = VTENSOR(K,K1)+DXYZ(K1)*DXYZ(K)*FORTOT
                      END DO
                   END DO
                   !-------------------------------------------------
                 end if
            ENDDO
        END DO

        RETURN
  END SUBROUTINE CALFOR
  !*************************************************************************************


  !*********************************************************************************
  SUBROUTINE CALFORCE_FS_Force_Table2(SimBox, List)
  !***  PURPOSE:   to begine calculate the force
  !                compared to CALFORCE_FS_Force_Table, WIHTOUT Newton's third law consiered here
  !     INPUT:     SimBox,    the simulation box
  !                List,      the neighbore list
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
   use MD_NeighborsList
   implicit none
      !--- dummy vaiables
      type(SimMDBox),            intent(inout)::SimBox
      type(NEIGHBOR_LIST),target,intent(in)   ::List


      !--- Local variables
               !$$-- keep the pointer of the current neioghbore list
                 m_NLISTP =>List

               !$$-- The boxsize could be changed due to pressure
               !$$    we should update it at each time
                m_BOX =  SimBox%ZL
                m_HBOX = SimBox%ZL*C_HALF
                m_BOXSHAPE = SimBox%BOXSHAPE

                !$$-- To calculate the electron densitis
                call PRECALFOR2(SimBox%NPRT, List%KVOIS, List%INDI, SimBox%ITYP, SimBox%XP)

                !$$-- To calculate the electron densitis
                call CALFOR2(SimBox%NPRT, List%KVOIS,List%INDI, SimBox%ITYP, SimBox%XP, SimBox%FP, SimBox%VTENSOR)

                !$$-- at the same time we update potential energy of the atoms
                 SimBox%EPOT = m_ER - DSQRT(m_DEN)

  end SUBROUTINE CALFORCE_FS_Force_Table2
  !*********************************************************************************


  !*********************************************************************************
  SUBROUTINE preCALFOR2(IM,KVOIS,INDI,ITYP,XP)
  !***  PURPOSE:   to begine calculate the pairwise potential m_ER and the electron density m_DEN
  !                compared to preCALFOR, WIHTOUT Newton's third law consiered
  !     INPUT:     IM:     the number of particles concerned
  !                KVOIS:  the number of for each atom in style1
  !                INDI:   the index for neighbores
  !                ITYPE:  the type of atom corresponding to INDI
  !
  !     OUTPUT     m_ER
  !                m_DEN

  !
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, intent(in)::IM
      integer, dimension(:), intent(in)::KVOIS
      integer, dimension(:,:), intent(in)::INDI
      integer, dimension(:), intent(in)::ITYP
      real(KINDDF), dimension(:,:), intent(in)::XP

      !--Local variables
      integer::I,J,K,KK, IW, IIW, KTAB, K1, ITYPI, ITYPJ
      real(KINDDF)::SK,DK,EXP1,EXP2
      real(KINDDF)::SEP(3), DXYZ(3), R2, R


           m_ER = C_ZERO
           m_DEN = C_ZERO

           DO I=1, IM
              IIW = KVOIS(I)
              ITYPI  = ITYP(I)
              DO IW=1, IIW
                 J=INDI(I,IW)

                 !$$-- To calculate the seperation between particle I and its IWth neighbore
                 DO K=1, 3
                    SEP(K) = XP(I,K)-XP(J,K)
                    IF((m_IFPD(K).GT.0) .AND. (DABS(SEP(K)) .GT. m_HBOX(K))) THEN
                       SEP(K) = SEP(K) - DSIGN(m_BOX(K),SEP(K))
                    END IF
                 END DO

                 !$$-- NOTE: To converte the length-unit into absolute unit (cm)
                 !$$          only when the shape of the box is cubic,
                 !$$          DXYZ = SEP
                 DO K=1, 3
                    DXYZ(K) =  sum(m_BOXSHAPE(K,1:3)*SEP(1:3))
                 END DO
                 R2  = sum(DXYZ*DXYZ)
                 !$$-- To calculate electron density on atom I
                 ITYPJ=ITYP(J)
                 if(R2 .LE. m_RCUT2(ITYPI,ITYPJ)) then
                    KTAB = m_KPAIR(ITYPI,ITYPJ)
                    !$$SK=R2*m_CSI
                    !$$SK = (C_TWO*m_Rmax*DSQRT(R2)-R2)*m_CSI
                    !$$--- Change made 2010-12-05
                    R = DSQRT(R2)
                    SK= DSQRT(R)*m_CSI

                    KK=SK
                    K1=KK+1
                    DK=SK-DBLE(KK)
                    !EXP1 = m_POTR(KTAB,KK)+DK*(m_POTR(KTAB,K1)-m_POTR(KTAB,KK))
                    !EXP2 = m_POTB(KTAB,KK)+DK*(m_POTB(KTAB,K1)-m_POTB(KTAB,KK))
                    !m_ER(I) = m_ER(I)+EXP1
                    !m_DEN(I)= m_DEN(I)+EXP2

                    !$$--- NOTE by HOU Qing: Dec 4,2014
                    !$$          The force table POTR, POTB are V(r), RHO(r)*0.5, in OLD VERSION
                    !$$          In the new version. we define POTR, POTB as V(r)*r and RHO(r)
                    m_ER(I) = m_ER(I)  + (m_POTR(KTAB,KK)+DK*(m_POTR(KTAB,K1)-m_POTR(KTAB,KK)))/R
                    m_DEN(I)= m_DEN(I) + (m_POTB(KTAB,KK)+DK*(m_POTB(KTAB,K1)-m_POTB(KTAB,KK)))

                end if
            ENDDO
         END DO

         !$$*** to get the derivetive of DSQRT(m_DEN)
          DO I=1, IM
             call GET_DIFF_FFUN_FS(m_DEN(I), m_DDEN(I))
         END DO

        RETURN
  END SUBROUTINE preCALFOR2
  !*************************************************************************************

  !*************************************************************************************
  SUBROUTINE CALFOR2(IM, KVOIS,INDI,ITYP,XP, FP,VTENSOR)
  !***  PURPOSE:   to begine calculate the force, Newton's third law taken into account
  !                compared to CALFOR, WIHTOUT Newton's third law consiered
  !     INPUT:     IM:     the number of particles
  !                KVOIS:  the number of neighbores in stle0
  !                INDI:   the index for neighbores
  !                ITYP:   the type of atom corresponding to INDI
  !                XP:     the coordinate of atoms
  !     OUTPUT     FP:     the force on atoms
  !                VTENSOR: the virial tensor
  !                EPOT:   the potential of a atom
  !
  use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, intent(in)::IM
      integer, dimension(:),intent(in)::KVOIS
      integer, dimension(:,:), intent(in)::INDI
      integer, dimension(:), intent(in)::ITYP
      real(KINDDF), dimension(:,:),  intent(in)::XP
      real(KINDDF), dimension(:,:), intent(out)::FP
      real(KINDDF), intent(out)::VTENSOR(3,3)

         !---Local variables
      integer::I,J,K,K1,N,IW,IIW,KTAB,KTAB1,KK, KK1,ITYPI, ITYPJ
      real(KINDDF)::SK,DK, EXP5,EXP6,EXP7, DENKI, DENKJ, FORTOT, FOR(3)
      real(KINDDF)::SEP(3), DXYZ(3), R2, R


        !$$--- now, we have the electron densities on atoms, to calculate the force
          FP = C_ZERO
          VTENSOR = C_ZERO
          DO I=1, IM
             IIW = KVOIS(I)
             ITYPI  = ITYP(I)
             DENKI = m_DDEN(I)

             !$$-- to begine scanning the neighbores
             DO IW=1, IIW
                !$$-- To calculate the seperation between particle I and its IWth neighbore
                J=INDI(I,IW)
                DO K=1, 3
                   SEP(K) = XP(I,K)-XP(J,K)
                   IF(m_IFPD(K) .AND. DABS(SEP(K)) .GT. m_HBOX(K)) THEN
                      SEP(K) = SEP(K) - DSIGN(m_BOX(K),SEP(K))
                   END IF
                END DO

                 !$$-- NOTE: To converte the length-unit into absolute unit (cm)
                 !$$          only when the shape of the box is cubic,
                 !$$          DXYZ = SEP
                DO K=1, 3
                   DXYZ(K) =  sum(m_BOXSHAPE(K,1:3)*SEP(1:3))
                END DO
                R2  = sum(DXYZ*DXYZ)

                !$$-- IF R2 is smaller than cutoff range, to calculate the force
                ITYPJ=ITYP(J)
                if(R2 .LE. m_RCUT2(ITYPI,ITYPJ)) then
                   KTAB  = m_KPAIR(ITYPI,ITYPJ)
                   KTAB1 = m_KPAIR(ITYPJ,ITYPI)
                   !$$SK = R2*m_CSI
                   !$$SK = (C_TWO*m_Rmax*DSQRT(R2)-R2)*m_CSI
                   !$$--- Change made 2010-12-05
                   R = DSQRT(R2)
                   SK= DSQRT(R)*m_CSI
                   KK = SK
                   KK1 = KK+1
                   DK=SK-dble(KK)
                   !$$-- interpolation of force table
                   !$$-- for  pairwise section
                   !EXP5 = m_FPOTR(KTAB,KK)  + DK*(m_FPOTR(KTAB,KK1) - m_FPOTR(KTAB,KK))

                   !$$-- for many-body section
                   !$$    Note; the electron density of atom I on atom J could not be the same as
                   !$$    from atom J on atom J, if atom I and J are not the same type.
                   !
                   !$$--- from the type of atom I to the type of atom J
                   !EXP6 = m_FPOTB(KTAB,KK)  + DK*(m_FPOTB(KTAB,KK1) - m_FPOTB(KTAB,KK))

                   !$$--- from the type of atom J to the type of atom I
                   !EXP7 = m_FPOTB(KTAB1,KK) + DK*(m_FPOTB(KTAB1,KK1) - m_FPOTB(KTAB1,KK))

                   !FORTOT= EXP5 + EXP7*DENKI + EXP6* m_DDEN(J)

                   !$$--- NOTE by Hou Qing: Dec 4,2014
                   !$$          The force table FPOTR, FPOTB are dV(r)/dr/r, and dRHO/dr/r, in OLD VERSION
                   !$$          In the new version. we define FPOTR, FPOTB as (dV(r)/dr)*r and (dRHO/dr)*r
                   FORTOT= (m_FPOTR(KTAB,KK)    + DK*(m_FPOTR(KTAB,KK1)  - m_FPOTR(KTAB,KK)))/R2       +   &
                           ( (m_FPOTB(KTAB1,KK) + DK*(m_FPOTB(KTAB1,KK1) - m_FPOTB(KTAB1,KK)))*DENKI   +   &
                             (m_FPOTB(KTAB,KK)  + DK*(m_FPOTB(KTAB,KK1)  - m_FPOTB(KTAB,KK)))*m_DDEN(J) )/R

                   DO K=1, 3
                      FOR(K) = FORTOT*SEP(K)
                      FP(I,K) = FP(I,K)+FOR(K)
                   END DO

                   !$$- at the same time, we can calculate the pressure tensor.
                   !$$   NOTE: here we have consider periodic condition, thus
                   !$$         we cannot laterly use PRESS=SUM(XP*FP) to get the pressure.
                   !$$         It would be convenience and efficient to have the pressure
                   !$$         calculated here.
                   DO K=1, 3
                      DO K1=1, 3
                         VTENSOR(K,K1) = VTENSOR(K,K1)+DXYZ(K1)*DXYZ(K)*FORTOT*C_HALF
                      END DO
                   END DO
                   !-------------------------------------------------
                 end if
            ENDDO
        END DO

        RETURN
  END SUBROUTINE CALFOR2
  !*************************************************************************************
  !  The following section supplies the defferential of F(RHO) in
  !  EAM or modified FS tight-binding
  !
  !*************************************************************************************
  SUBROUTINE GET_DIFF_FFUN_FS(DEN, DDEN)
  !***  PURPOSE:
  !
  use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      real(KINDDF), intent(in)::DEN
      real(KINDDF), intent(out)::DDEN

         !$$*** to get the derivetive of DSQRT(m_DEN)
         !$$--- NOTE by Hou Qing, Dec 4,2014
         !$$       in old version  the force term FPOTB has been multiplied by 0.5
         !$$       in new version, we donnot multiply 0.5 on FPOTB
         !$$       thus we should chnage
         !$$         m_DDEN(I) = -C_UN/DSQRT(m_DEN(I))
         !$$       to
         !$$         m_DDEN(I) = -C_HLAF/DSQRT(m_DEN(I))
             IF(DEN .GT. 0) then
                !DDEN = -C_UN/DSQRT(DEN)
                DDEN = -C_HALF/DSQRT(DEN)
             ELSE
                DDEN = C_ZERO
             END IF
  END SUBROUTINE GET_DIFF_FFUN_FS
  !*************************************************************************************

  !*************************************************************************************
  logical function IfInit_FS_Force_Table() result(YES)
  implicit none
          YES = m_INITED
          return
  end function IfInit_FS_Force_Table
   !*************************************************************************************

  end module MD_FS_Force_Table





