  module MD_NeighborsList
  !**** DESCRIPTION: to calculate the neighbor list. instead of using
  !                  matrix INDI and KVOIS that are used in Cal_Neighbors_List module
  !                  (CalNeighbores.f90), we use user defined data type LIST whick
  !                   encapsual these two index. The size for INDI could be various for
  !                   each atoms
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010

  !****
  use MD_CONSTANTS
  use MiniUtilities
  implicit none


      type::NEIGHBOR_LIST
           integer, dimension(:,:),allocatable::INDI    ! the index of neighbores
           integer, dimension(:),  allocatable::KVOIS     ! the number of neighbores at current time step
           integer::mxKVOIS=400                         ! maxmum permitted number of neighbores
      end type

  contains

  !************************************************************************************************
  subroutine Clear_NeighboreList(List)
  !***  PURPOSE:  to deallocate the memory allocated in List
  !
  !     INPUT:     List, the neighbore list
  !     OUTPUT     List, the neighbore list with memory deallocated

  use MD_CONSTANTS
  implicit none
  !----   DUMMY Variables
          type(NEIGHBOR_LIST)::List
  !----   Local variables
  !***
             if(allocated(List%KVOIS)) then
                deallocate(List%KVOIS)
             end if

             if(allocated(List%INDI)) then
                deallocate(List%INDI)
             end if


        RETURN
  end subroutine Clear_NeighboreList
  !************************************************************************************************


  !************************************************************************************************
  subroutine Cal_NeighboreList(SimBox, CtrlParam, List)
  !***  PURPOSE:  to update the neighbore list of atoms
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                List, the neighbore list
  !     OUTPUT     SimBox,    the simulation box with new position and velocity

  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)     ::SimBox
          type(SimMDCtrl)    ::CtrlParam
          type(NEIGHBOR_LIST)::List

  !----   Local variables
         integer I, J, K, IW, NPART, MXKVOIS, ITYP, JTYP
         real(KINDDF)::CXYZ(3), BOXSIZE(3), HBOXSIZE(3)
         real(KINDDF)::RCUT2(10,10)
         integer, dimension(:), allocatable::INDI
         logical::PD(3)

  !***
         NPART    = SimBox%NPRT
         BOXSIZE =  SimBox%ZL
         HBOXSIZE =  BOXSIZE*C_HALF

         PD = CtrlParam%IFPD
         DO I=1, SimBox%NGROUP
            DO J=1, SimBox%NGROUP
               RCUT2(I,J) = CtrlParam%NB_RM(I,J)*CtrlParam%NB_RM(I,J)
            END DO
         END DO

         !$$--- to check the consistent of allocated memory
         if(List%mxKVOIS .NE. CtrlParam%NB_MXNBS) then
             if(allocated(List%KVOIS)) then
                deallocate(List%KVOIS)
             end if

             if(allocated(List%INDI)) then
                deallocate(List%INDI)
             end if
             List%mxKVOIS = CtrlParam%NB_MXNBS
         end if

         if(allocated(List%KVOIS)) then
             if(size(List%KVOIS) .ne. NPART) then
                deallocate(List%KVOIS)
                deallocate(List%INDI)
             end if
         end if

         !$$--- check if the memory has been allocated
         if(.not.allocated(List%KVOIS)) then
             allocate(List%KVOIS(NPART))
             List%KVOIS = 0
         end if

         if(.not.allocated(List%INDI)) then
             allocate(List%INDI(NPART,List%mxKVOIS))
         end if

         !---
         allocate(INDI(NPART))
         MXKVOIS = List%mxKVOIS
         DO I=1,NPART
            IW=C_IZERO
            DO J=I+1,NPART
               DO K=1, 3
                  CXYZ(K) = SimBox%XP(I,K) - SimBox%XP(J,K)
                  IF(DABS(CXYZ(K)).GT.HBOXSIZE(K) .and. PD(K)) CXYZ(K)=CXYZ(K)-DSIGN(BOXSIZE(K),CXYZ(K))
               END DO

               ITYP = SimBox%ITYP(I)
               JTYP = SimBox%ITYP(J)
               IF(SUM(CXYZ*CXYZ).LT.RCUT2(I,J)) THEN
                  IW=IW+1
                  IF(IW .GT. MXKVOIS) THEN
                     WRITE(*,*) "MDPSCU Error: Number of neighbores larger than permitted value:",MXKVOIS
                     WRITE(*,*) "Process to be stoped"
                     stop
                  END IF
                  INDI(IW)=J
               END IF
           ENDDO
           !---
           List%KVOIS(I)=IW
           List%INDI(I,1:IW) = INDI(1:IW)
        END DO
        deallocate(INDI)


        RETURN
  end subroutine Cal_NeighboreList
  !************************************************************************************************

  !************************************************************************************************
  subroutine Cal_NeighboreListC(SimBox, CtrlParam, List)
  !***  PURPOSE:  to update the neighbore list of atoms with head-link
  !               technique used.
  !               NOTE: Newton's 3rd law is considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should be applied

  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                List, the neighbore list
  !     OUTPUT     SimBox,    the simulation box with new position and velocity

  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)     ::SimBox
          type(SimMDCtrl)    ::CtrlParam
          type(NEIGHBOR_LIST)::List

  !----   Local variables
         integer I, J, K, KK, N, ID, NPART, PD(3),NCELL(3),NC, IXYZ(3), IX, IY, IZ, IC0, IC, JXYZ(3), N0, NN, MXKVOIS, IERR,ITP, JTP, MISSED
         real(KINDDF)::BOXSIZE(3), LBOX(3),BOXSHAPE(3,3), RCUT2(10,10), CXYZ(3), CXYZ1(3)
         integer, dimension(:), allocatable::HEAD
         integer, dimension(:), allocatable::LINK    ! link list
         integer, dimension(:), allocatable::IDENT
         real(KINDDF), dimension(:,:), allocatable::XPT
         integer, dimension(:), allocatable::INDI
         integer, dimension(:), allocatable::ITYP

         integer,parameter::NIX(14)=(/0,-1,-1,-1, 0, 0, -1, 1,-1, 0, 1,-1, 0, 1/)
         integer,parameter::NIY(14)=(/0, 0,-1, 1, 1, 0,  0, 0,-1,-1,-1, 1, 1, 1/)
         integer,parameter::NIZ(14)=(/0, 0, 0, 0, 0, 1,  1, 1, 1, 1, 1, 1, 1, 1/)
         real(KINDDF), parameter::eps=1.d-04
  !***
         NPART    = SimBox%NPRT
         BOXSIZE =  SimBox%ZL
         LBOX = SimBox%BOXLOW
         BOXSHAPE = SimBox%BOXSHAPE
         PD = CtrlParam%IFPD
         DO I=1, SimBox%NGROUP
            DO J=1, SimBox%NGROUP
               RCUT2(I,J) = CtrlParam%NB_RM(I,J)*CtrlParam%NB_RM(I,J)
            END DO
         END DO

         !$$--- to check the consistent of allocated memory
         if(List%mxKVOIS .NE. CtrlParam%NB_MXNBS) then
             if(allocated(List%KVOIS)) then
                deallocate(List%KVOIS)
             end if

             if(allocated(List%INDI)) then
                deallocate(List%INDI)
             end if
             List%mxKVOIS = CtrlParam%NB_MXNBS
         end if

         if(allocated(List%KVOIS)) then
             if(size(List%KVOIS) .ne. NPART) then
                deallocate(List%KVOIS)
                deallocate(List%INDI)
             end if
         end if

         !$$--- check if the memory has been allocated
         if(.not.allocated(List%KVOIS)) then
             allocate(List%KVOIS(NPART))
             List%KVOIS = 0
         end if

         if(.not.allocated(List%INDI)) then
             allocate(List%INDI(NPART,List%mxKVOIS))
         end if

         !$$--- to determine how many cells we need
         DO K=1,3
            !if(PD(K)) then
               NCELL(K) = int(BOXSIZE(K)/(1.0D0*maxval(CtrlParam%NB_RM)) - EPS)
               !$$if(NCELL(K) .lt. CtrlParam%NC(K)) NCELL(K)=CtrlParam%NC(K)
               if(NCELL(K) .lt. C_ITHR) NCELL(K) = C_ITHR
            !else
            !   NCELL(K) = CtrlParam%NC(K)
            !   if(NCELL(K) .lt. C_ITHR) NCELL(K) = C_ITHR
            !end if
         END DO

         NC = NCELL(1)*NCELL(2)*NCELL(3)
         !$$--- allocate the temperory memory
         allocate(HEAD(NC),LINK(NPART),IDENT(NPART),XPT(NPART,3),INDI(NPART),ITYP(NPART),STAT=IERR)
         if(IERR) then
            write(*,*) "Fail to allocate memory in neighbore calculation."
            write(*,*) "Process to be stoped"
            stop
         end if

         !------------------------------------------------------------------------------------
         !$$--- to create the LINKED-CELL
         LINK = 0
         HEAD = 0
         MXKVOIS = List%mxKVOIS
         MISSED = 0
         do I=1, NPART
            if(IAND(SimBox%STATU(I),CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX) then
               do K=1,3
                  IXYZ(K) = int( (SimBox%XP(I,K)-LBOX(K))/BOXSIZE(K)*dble(NCELL(K)) -eps)
               end do

               !$$--- if a particle is marked as active and outof box
               if((IXYZ(1) .lt. 0 .or. IXYZ(1) .gt. NCELL(1)) .or. &
                  (IXYZ(2) .lt. 0 .or. IXYZ(2) .gt. NCELL(2)) .or. &
                  (IXYZ(3) .lt. 0 .or. IXYZ(3) .gt. NCELL(3))) then
                     write(*,fmt="(A, I8, A)") "MDPSCU Error: particle ",I," marked as active but out of box"
                     write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "              box range x ", &
                                                                        LBOX(1), LBOX(1)+BOXSIZE(1), ", particle x ", SimBox%XP(I, 1)
                     write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "              box range y ", &
                                                                        LBOX(2), LBOX(2)+BOXSIZE(2), ", particle y ", SimBox%XP(I, 2)
                     write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "              box range z ", &
                                                                        LBOX(3), LBOX(3)+BOXSIZE(3), ", particle z ", SimBox%XP(I, 3)
                     write(*,fmt="(A)")"Process to be stopped"
                     stop
               end if

               IC = C_IUN+(IXYZ(1)+NCELL(1)*(IXYZ(2)+NCELL(2)*IXYZ(3)) )
               J = HEAD(IC)
               HEAD(IC) = I
               LINK(I)   = J
           else
             if(IAND(SimBox%STATU(I),CP_STATU_OUTOFBOX) .eq. CP_STATU_OUTOFBOX) then
                List%KVOIS(I) = 0
                MISSED = MISSED + 1
             end if
           end if
         end do

          if(MISSED .GT. 0) then
            write(*,fmt="(A, I8, A)") " MDPSCU Warning: there are ",MISSED," atoms out of box found in MD_NeighborsList."
            call ONWARNING(gm_OnWarning)
          end if

         !------------------------------------------------------------------------------------
         !$$--- NOW SCAN ON ALL CELLS
         DO IZ=1, NCELL(3)
            DO IY=1, NCELL(2)
               DO IX=1, NCELL(1)
                  !$$--- for the first, we get the particles in cell(IX, IY, IZ)
                    IC0 = IX+NCELL(1)*((IY-1)+NCELL(2)*(IZ-1))
                    ID = HEAD(IC0)
                    IF(ID.EQ. C_IZERO) then
                       cycle !$$ no particle in the cell
                    end if

                    N = 0
                    DO WHILE(.TRUE.)
                       N = N + 1
                       IDENT(N) = ID
                       XPT(N,1:3) = SimBox%XP(ID,1:3)
                       ITYP(N) = SimBox%ITYP(ID)
                       ID = LINK(ID)
                       IF(ID.LE.C_IZERO) EXIT
                    END DO
                    !$$--- record the number of particles in cell(IX, IY, IZ)
                    N0 = N
                    !$$-- Now we perform a loop over the neighbouring cells
                    DO IC = 2,14
                       JXYZ(3) = IZ+NIZ(IC)
                       JXYZ(2) = IY+NIY(IC)
                       JXYZ(1) = IX+NIX(IC)
                       !$$---  check PBC's : minimum image convention
                       CXYZ =C_ZERO
                        DO KK=1,3
                           If(PD(KK)) Then
                              IF( JXYZ(KK).GT.NCELL(KK) )THEN
                                  JXYZ(KK) = C_IUN
                                  CXYZ(KK) = BOXSIZE(KK)
                             ELSE IF (JXYZ(KK).LT.C_IUN) THEN
                                  JXYZ(KK) = NCELL(KK)
                                  CXYZ(KK) = -BOXSIZE(KK)
                             ENDIF
                           End If
                       END DO

                       IF( ANY(JXYZ .GT. NCELL) ) cycle
                       IF( ANY(JXYZ .LT. C_IUN )) cycle

                      !$$--- NO WE HAVE A NEIGHBOR CELL, ADD ITS PARTICLES TO THE SUBBOX
                       ID  = HEAD(JXYZ(1)+NCELL(1)*((JXYZ(2)-1)+NCELL(2)*(JXYZ(3)-1)))
                       IF(ID.GT.C_IZERO) THEN
                          DO WHILE(.TRUE.)
                             N = N+1
                             IDENT(N) = ID
                             XPT(N,1:3) = SimBox%XP(ID,1:3)+CXYZ(1:3)
                             ITYP(N) = SimBox%ITYP(ID)
                             ID = LINK(ID)
                             IF(ID.LE.C_IZERO) EXIT
                           END DO
                       END IF
                    END DO !$$ End the loop of IC

                    !__________________________________________________________________
                    !$$   Now we have all the particles in the cells, to  start the calculation
                    !$$   of neighbore list
                    !$$   NOTE: we calculate the neighbores in present cell only
                    DO I=1,N0
                       ID = IDENT(I)
                       ITP = ITYP(I)
                       NN = 0
                       DO J=I+1,N
                          CXYZ1(1:3) = XPT(I,1:3) - XPT(J,1:3)
                          CXYZ(1)    = sum(BOXSHAPE(1,1:3)*CXYZ1(1:3))
                          CXYZ(2)    = sum(BOXSHAPE(2,1:3)*CXYZ1(1:3))
                          CXYZ(3)    = sum(BOXSHAPE(3,1:3)*CXYZ1(1:3))
                          JTP = ITYP(J)
                          IF(SUM(CXYZ*CXYZ).LT.RCUT2(ITP,JTP)) THEN
                             NN = NN+1
                             IF(NN .GT. MXKVOIS) THEN
                                WRITE(*,*) "MDPSCU Error: Number of neighbores larger than permitted value:",NN, MXKVOIS
                                WRITE(*,*) "Process to be stoped"
                                stop
                             END IF
                             INDI(NN) = IDENT(J)
                          END IF
                       END DO !End loop for J
                      !---
                      List%KVOIS(ID) = NN
                      List%INDI(ID,1:NN) = INDI(1:NN)
                   END DO !$$End loop for I
               END DO !$$End loop for IX cell
            END DO !$$End loop for IY cell
         END DO  !$$End loop for IZ cell

        deallocate(HEAD,LINK,IDENT,XPT,INDI,ITYP,STAT=IERR)
        if(IERR) then
            write(*,*) "MDPSCU Warning:Fail to deallocate memory in neighbore calculation."
            call ONWARNING(gm_OnWarning)
        end if

        RETURN
  end subroutine Cal_NeighboreListC
  !************************************************************************************************


  !************************************************************************************************
  subroutine Cal_NeighboreList2C(SimBox, CtrlParam, List)
  !***  PURPOSE:  to update the neighbore list of atoms with head-link
  !               technique used.
  !               NOTE: Newton's 3rd law NOT is considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should be applied

  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                List, the neighbore list
  !     OUTPUT     SimBox,    the simulation box with new position and velocity

  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)     ::SimBox
          type(SimMDCtrl)    ::CtrlParam
          type(NEIGHBOR_LIST)::List

  !----   Local variables
         integer I, J, K, KK, N, ID, NPART, PD(3),NCELL(3),NC, IXYZ(3), IX, IY, IZ, IC0, IC, JXYZ(3), N0, NN, MXKVOIS,IERR,ITP, JTP, MISSED
         real(KINDDF)::BOXSIZE(3), LBOX(3), BOXSHAPE(3,3), RCUT2(10,10), CXYZ(3), CXYZ1(3)
         integer, dimension(:), allocatable::HEAD
         integer, dimension(:), allocatable::LINK    ! link list
         integer, dimension(:), allocatable::IDENT
         real(KINDDF), dimension(:,:), allocatable::XPT
         integer, dimension(:), allocatable::INDI
         integer, dimension(:), allocatable::ITYP

         integer,parameter::NIX(27)=(/0,-1,-1,-1, 0, 0, -1, 1,-1, 0, 1,-1, 0, 1, 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/)
         integer,parameter::NIY(27)=(/0, 0,-1, 1, 1, 0,  0, 0,-1,-1,-1, 1, 1, 1, 0, 1,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/)
         integer,parameter::NIZ(27)=(/0, 0, 0, 0, 0, 1,  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
         real(KINDDF), parameter::eps=0.0001
  !***
         NPART    = SimBox%NPRT
         BOXSIZE =  SimBox%ZL
         !HBOXSIZE =  BOXSIZE*C_HALF
         LBOX = SimBox%BOXLOW
         BOXSHAPE = SimBox%BOXSHAPE
         PD = CtrlParam%IFPD
         DO I=1, SimBox%NGROUP
            DO J=1, SimBox%NGROUP
               RCUT2(I,J) = CtrlParam%NB_RM(I,J)*CtrlParam%NB_RM(I,J)
            END DO
         END DO

         !$$--- to check the consistent of allocated memory
         if(List%mxKVOIS .NE. CtrlParam%NB_MXNBS) then
             if(allocated(List%KVOIS)) then
                deallocate(List%KVOIS)
             end if

             if(allocated(List%INDI)) then
                deallocate(List%INDI)
             end if
             List%mxKVOIS = CtrlParam%NB_MXNBS
         end if

         if(allocated(List%KVOIS)) then
             if(size(List%KVOIS) .ne. NPART) then
                deallocate(List%KVOIS)
                deallocate(List%INDI)
             end if
         end if

         !$$--- check if the memory has been allocated
         if(.not.allocated(List%KVOIS)) then
             allocate(List%KVOIS(NPART))
             List%KVOIS = 0
         end if

         if(.not.allocated(List%INDI)) then
             allocate(List%INDI(NPART,List%mxKVOIS))
         end if

         !$$--- to determine how many cells we need
         DO K=1,3
            !if(PD(K)) then
               NCELL(K) = int(BOXSIZE(K)/(1.0D0*maxval(CtrlParam%NB_RM)) - EPS)
               !$$if(NCELL(K) .lt. CtrlParam%NC(K)) NCELL(K)=CtrlParam%NC(K)
               if(NCELL(K) .lt. C_ITHR) NCELL(K) = C_ITHR
            !else
            !   NCELL(K) = CtrlParam%NC(K)
            !   if(NCELL(K) .lt. C_ITHR) NCELL(K) = C_ITHR
            !end if
         END DO

         NC = NCELL(1)*NCELL(2)*NCELL(3)
         !$$--- allocate the temperory memory
         allocate(HEAD(NC),LINK(NPART),IDENT(NPART),XPT(NPART,3),INDI(NPART),ITYP(NPART))

         !------------------------------------------------------------------------------------
         !$$--- to create the LINKED-CELL
         LINK = 0
         HEAD = 0
         MXKVOIS = List%mxKVOIS
         MISSED = 0
         do I=1, NPART
            if(IAND(SimBox%STATU(I),CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX) then
               do K=1,3
                  IXYZ(K) = int( (SimBox%XP(I,K)-LBOX(K))/BOXSIZE(K)*dble(NCELL(K)) -eps)
               end do

               !$$--- if a particle is marked as active and outof box
               if((IXYZ(1) .lt. 0 .or. IXYZ(1) .gt. NCELL(1)) .or. &
                  (IXYZ(2) .lt. 0 .or. IXYZ(2) .gt. NCELL(2)) .or. &
                  (IXYZ(3) .lt. 0 .or. IXYZ(3) .gt. NCELL(3))) then
                     write(*,fmt="(A, I8, A)") "MDPSCU Error: particle ",I," marked as active but out of box"
                     write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "              box range x ", &
                                                                        LBOX(1), LBOX(1)+BOXSIZE(1), ", particle x ", SimBox%XP(I, 1)
                     write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "              box range y ", &
                                                                        LBOX(2), LBOX(2)+BOXSIZE(2), ", particle y ", SimBox%XP(I, 2)
                     write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "              box range z ", &
                                                                        LBOX(3), LBOX(3)+BOXSIZE(3), ", particle z ", SimBox%XP(I, 3)
                     write(*,fmt="(A)")"Process to be stopped"
                     stop
               end if

               IC = C_IUN+(IXYZ(1)+NCELL(1)*(IXYZ(2)+NCELL(2)*IXYZ(3)) )
               J = HEAD(IC)
               HEAD(IC) = I
               LINK(I)   = J
           else
             if(IAND(SimBox%STATU(I),CP_STATU_OUTOFBOX) .eq. CP_STATU_OUTOFBOX) then
                List%KVOIS(I) = 0
                MISSED = MISSED + 1
             end if
           end if
         end do

         if(MISSED .GT. 0) then
            write(*,fmt="(A, I8, A)") " MDPSCU Warning: there are ",MISSED," atoms out of box found in MD_NeighborsList."
            call ONWARNING(gm_OnWarning)
          end if

         !------------------------------------------------------------------------------------
         !$$--- NOW SCAN ON ALL CELLS
         DO IZ=1, NCELL(3)
            DO IY=1, NCELL(2)
               DO IX=1, NCELL(1)
                  !$$--- for the first, we get the particles in cell(IX, IY, IZ)
                    IC0 = IX+NCELL(1)*((IY-1)+NCELL(2)*(IZ-1))
                    ID = HEAD(IC0)
                    IF(ID.EQ. C_IZERO) cycle !$$ no particle in the cell

                    N = 0
                    DO WHILE(.TRUE.)
                       N = N + 1
                       IDENT(N) = ID
                       XPT(N,1:3) = SimBox%XP(ID,1:3)
                       ITYP(N) = SimBox%ITYP(ID)
                       ID = LINK(ID)
                       IF(ID.LE.C_IZERO) EXIT
                    END DO
                    !$$--- record the number of particles in tcell(IX, IY, IZ)
                    N0 = N
                    !$$-- Now we perform a loop over the neighbouring cells
                    DO IC = 2,27
                       JXYZ(3) = IZ+NIZ(IC)
                       JXYZ(2) = IY+NIY(IC)
                       JXYZ(1) = IX+NIX(IC)
                       !$$---  check PBC's : minimum image convention
                       CXYZ =C_ZERO
                        DO KK=1,3
                           If(PD(KK)) Then
                              IF( JXYZ(KK).GT.NCELL(KK) )THEN
                                  JXYZ(KK) = C_IUN
                                  CXYZ(KK) = BOXSIZE(KK)
                             ELSE IF (JXYZ(KK).LT.C_IUN) THEN
                                  JXYZ(KK) = NCELL(KK)
                                  CXYZ(KK) = -BOXSIZE(KK)
                             ENDIF
                           End If
                       END DO

                       IF( ANY(JXYZ .GT. NCELL) ) cycle
                       IF( ANY(JXYZ .LT. C_IUN )) cycle

                      !$$--- NO WE HAVE A NEIGHBOR CELL, ADD ITS PARTICLES TO THE SUBBOX
                       ID  = HEAD(JXYZ(1)+NCELL(1)*((JXYZ(2)-1)+NCELL(2)*(JXYZ(3)-1)))
                       IF(ID.GT.C_IZERO) THEN
                          DO WHILE(.TRUE.)
                             N = N+1
                             IDENT(N) = ID
                             XPT(N,1:3) = SimBox%XP(ID,1:3)+CXYZ(1:3)
                             ITYP(N) = SimBox%ITYP(ID)
                             ID = LINK(ID)
                             IF(ID.LE.C_IZERO) EXIT
                           END DO
                       END IF
                    END DO !$$ End the loop of IC

                    !__________________________________________________________________
                    !$$  Now we have all the particles in the cells, to  start the calculation
                    !$$   of neighbore list
                    !$$   NOTE: we calculate the neighbores in present cell only
                    DO I=1,N0
                       ID = IDENT(I)
                       ITP = ITYP(I)
                       NN = 0
                       DO J=1,N
                          if(J.NE.I) then
                             CXYZ1(1:3) = XPT(I,1:3) - XPT(J,1:3)
                             CXYZ(1)    = sum(BOXSHAPE(1,1:3)*CXYZ1(1:3))
                             CXYZ(2)    = sum(BOXSHAPE(2,1:3)*CXYZ1(1:3))
                             CXYZ(3)    = sum(BOXSHAPE(3,1:3)*CXYZ1(1:3))
                             JTP = ITYP(J)
                             IF(SUM(CXYZ*CXYZ).LT.RCUT2(ITP,JTP)) THEN
                                NN = NN+1
                                !---
                                IF(NN .GT. MXKVOIS) THEN
                                   WRITE(*,*) "MDPSCU Error: Number of neighbores larger than permitted value:",MXKVOIS
                                   WRITE(*,*) "Process to be stoped"
                                   stop
                                END IF
                                INDI(NN) = IDENT(J)
                             END IF
                          end if
                       END DO !End loop for J

                      !---
                      List%KVOIS(ID) = NN
                      List%INDI(ID,1:NN) = INDI(1:NN)
                   END DO !End loop for I
               END DO !End loop for IX cell
            END DO !End loop for IY cell
         END DO  !End loop for IZ cell

        deallocate(HEAD,LINK,IDENT,XPT,INDI,ITYP,STAT=IERR)
        if(IERR) then
            write(*,*) "MDPSCU Warning: Fail to deallocate memory in neighbore calculation."
            call ONWARNING(gm_OnWarning)
        end if

        RETURN
  end subroutine Cal_NeighboreList2C
  !************************************************************************************************

  !************************************************************************************************
  subroutine Putout_NeighboreList(fname, List)
  !***  PURPOSE:  to output the neighbor list to file
  !
  !     INPUT:     fname,    the file name
  !                List,     the neighbore list
  !     OUTPUT

  use MD_CONSTANTS
  implicit none
  !----   DUMMY Variables
          character*(*)::fname
          type(NEIGHBOR_LIST)::List

  !----   Local variables
         integer::hFile, I, NPART, MXN, IW
  !***

            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname)

               write(hFile, fmt="(A)")    '!--- THE NEIGHBORE-LIST CREATED BY MDPSCU'
               write(hFile, fmt="(A)")    '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
               write(hFile, fmt="(A)")    '!    AUTHOR: HOU Qing'
               write(hFile, fmt="(A)")    '!    '
               write(hFile, fmt="(A)")    '&NEIGHBOR_LIST'
               if(allocated(List%KVOIS) ) then
                  NPART = size(List%KVOIS)
                  MXN   = maxval(List%KVOIS)
               else
                  NPART = 0
               end if

               write(hFile, fmt="(A, I8)")        '&NATOM ', NPART
               write(hFile, fmt="(A, I8, I8)")    '&MAXNB ', List%mxKVOIS, MXN
               do I=1, NPART
                  IW = List%KVOIS(I)
                  write(hFile, fmt="(I8, 2x, I8)")   I, IW
                  write(hFile,*) List%INDI(I, 1:IW)
               end do
             close(hFile)

        return
  end subroutine Putout_NeighboreList
  !************************************************************************************************

  !************************************************************************************************
  subroutine Load_NeighboreList(fname, List)
  !***  PURPOSE:  to load a neighbor list from file
  !
  !     INPUT:     fname,    the file name
  !                List,     the neighbore list
  !     OUTPUT
  implicit none
  !----   DUMMY Variables
          character*(*)::fname
          type(NEIGHBOR_LIST)::List

  !----   Local variables
         integer::hFile, I, IW, J, NPART, MXKVOIS, MXN, LINE, N
         character*256::STR
         character*64::KEYWORD, STRNUMB(2)
  !***
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old')

            LINE = 0
            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            STR = adjustl(STR)
            call GetKeyWord("&", STR, KEYWORD)
            call UpCase(KEYWORD)
            if( KEYWORD(1:LEN_TRIM(KEYWORD)) .ne. '&NEIGHBOR_LIST') goto 100

            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,1,n,STRNUMB)
            NPART = ISTR(STRNUMB(1))

            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,2,n,STRNUMB)
            MXKVOIS = ISTR(STRNUMB(1))
            MXN     = ISTR(STRNUMB(2))

            if(.not. allocated(List%KVOIS) ) then
               allocate(List%KVOIS(NPART), List%INDI(NPART, MXN))
            else
               if(List%MXKVOIS .lt. MXN) then
                  write(*,fmt="(A)")     "MDPSCU Error: The max permitted number of neihbore in List is maller than that of loade List"
                  write(*,fmt="(A, I, A, I)")     "              ", List%MXKVOIS, " vs ", MXN
                  write(*,fmt="(A)")     "              Process to be stopped"
                  stop
               end if
             end if

             List%KVOIS = 0
             do I=1, NPART
                  read(hFile, fmt="(I8, 2x, I8)")   J, IW
                  List%KVOIS(I) = IW
                  if(IW .gt. 0) read(hFile,*) List%INDI(I, 1:IW)
             end do
             close(hFile)
             return

100    continue
       write(*,fmt="(A)")     "MDPSCU Error: Wrong file format for the neighbor-lsit"
       write(*,fmt="(A)")     "              Keyword "//'&NEIGHBOR_LIST'//"not found"
       write(*,fmt="(A, I7)") "              Check file "//fname(1:len_trim(fname))//"at line", LINE
       write(*,fmt="(A)")     "              Process to be stopped"
       stop
  end subroutine Load_NeighboreList
  !************************************************************************************************




  end module MD_NeighborsList
