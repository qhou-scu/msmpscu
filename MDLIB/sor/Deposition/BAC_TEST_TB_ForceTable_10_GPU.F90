  module BAC_TEST_TB_ForceTable_GPU
  !**** DESCRIPTION: To calculate the force used for test BAC for energetic particle transport
  !                  in a substarte.
  !                  The force calculation is the same as in module MD_FINNIS_EM_TB_Force_Table_2010_GPU.
  !                  The force on energetic atom is calculated for only the neareast target atoms.
  !                  The number of the  neareast target atoms is defined by a parameter mp_NEARIST.

  !
  !                  SEE ALSO______________________________________________________________________________
  !                       MDTypeDef_ForceTable.f90
  !                       MD_FINNIS_EM_TB_Force_Table_2010_GPU.f90
  !                  ______________________________________________________________________________________
  !                  HOU Qing, Oct, 2011
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS_2010
    use MD_TYPEDEF_SIMULAIONBOX_2010
    use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
    use MD_TYPEDEF_FORCETABLE_2010
    use CUDAFOR
  !-----------------------------------------------
    implicit none

      integer, private::m_INITED = .false.
      !--- C2050
      !integer, parameter, private::mp_MXSHARESIZE = 3*3*224
      !--- C1060
      integer, parameter, private::mp_MXSHARESIZE = 3*3*64

      real(KINDDF),constant::dm_BOXSHAPE(3,3)                   ! siumulation box boundary, in unit of atomic radiua
      real(KINDDF), constant::dm_BOX(3)                         ! simulation box size
      real(KINDDF), constant::dm_HBOX(3)                        ! half box size

      integer, constant::dm_IFPD(3)                             !The periodic condition indicator


      real(KINDDF), private::m_Rmax                              ! the maxma value of R in the table
      real(KINDDF), private::m_CSI                               ! the inverse of RmaxSQRT e.g. NTAB/RmaxSQRT
      real(KINDDF), private::m_CSIV                              ! the inverse of CSI         1/CSI
      integer, parameter, private::mp_MXGROUP = 4
      integer, constant::dm_KPAIR(mp_MXGROUP,mp_MXGROUP)         ! the index of table for atom1 and atom2
      real(KINDDF), constant::dm_RU2(mp_MXGROUP,mp_MXGROUP)      ! the cutoff range
      real(KINDDF), device, dimension(:,:),allocatable::dm_POTR  ! the potential table for Nucl.-Nucl.  interaction with size of  NKIND*NTAB
      real(KINDDF), device, dimension(:,:),allocatable::dm_FPOTR ! the differential of potential table POTR
      real(KINDDF), device, dimension(:,:),allocatable::dm_POTB  ! the potential table
      real(KINDDF), device, dimension(:,:),allocatable::dm_FPOTB ! the differe


      real(KINDDF), device, dimension(:), allocatable::dm_DEN    ! the total electron density on atom SUMj(ROij)

      !**** to change the parameter mp_NEARIST to test how many neighbore we need
      integer, parameter, private::mp_NEARIST = 64

      !****
      private::PRECALFOR_KERNAL, CALFORCE_KERNAL, CALEPOT_KERNALA
  contains

  !*********************************************************************************
    SUBROUTINE INITIALIZE_FINNIS_EM_TB_Force_Table_DEV(SimBox, CtrlParam, FTable,RelaseTable, MULTIBOX)
    !***  PURPOSE:   to intialize the parameters to be used in force calculation
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !                FTable, the force tabled to be used
    !
    !     OUTPUT
     use MD_Globle_Variables_2010_GPU
    implicit none
       !--- dummy vaiables
      type(SimMDBox),            intent(inout)::SimBox
      type(SimulationCtrlParam), intent(in)   ::CtrlParam
      type(MDForceTable),        intent(in)   ::FTable
      integer,                   optional     ::RelaseTable, MULTIBOX
      !--- Local variables
      integer::ERR, NPRT

               call CLEAR_FINNIS_EM_TB_Force_Table_DEV()

           !--- allocate memory for electron density
               if(present(MULTIBOX)) then
                  NPRT = SimBox%NPRT*CtrlParam%MULTIBOX
               else
                  NPRT = SimBox%NPRT
               end if

               if(dm_NPRT .LT. NPRT) then
                  write(*,*) "Error:number of particles in device is smaller than Simulation box", dm_NPRT,NPRT
                  stop
               end if

               allocate(dm_DEN(NPRT), STAT=ERR)
               if(ERR) then
                  write(*,*) "ERROR: fail to allocate device memory in FINIS_EM_TB, 1"
               end if
           !--- copy to device constant memory
              dm_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
              dm_BOX(1:3)      = SimBox%ZL(1:3)
              dm_HBOX(1:3)     = SimBox%ZL(1:3)*C_HALF
              dm_IFPD(1:3)     = CtrlParam%IFPD(1:3)
              dm_RU2           = CtrlParam%RU*CtrlParam%RU

           !--- to copy force table information
               m_RMAX = FTable%RMAX
               m_CSI = FTable%CSI
               m_CSIV = FTable%CSIV

               allocate(dm_POTR(SIZE(FTable%POTR,dim=1),    SIZE(FTable%POTR,dim=2)),  &
                        dm_POTB(SIZE(FTable%POTB,dim=1),    SIZE(FTable%POTB,dim=2)),  &
                        dm_FPOTR(SIZE(FTable%FPOTR,dim=1),  SIZE(FTable%FPOTR,dim=2)), &
                        dm_FPOTB(SIZE(FTable%FPOTB,dim=1),  SIZE(FTable%FPOTB,dim=2)), STAT=ERR  )

               if(ERR) then
                  write(*,*) "ERROR: fail to allocate device memory in FINIS_EM_TB, 2"
               end if

               dm_KPAIR(1:SIZE(FTable%KPAIR,dim=1),  1:SIZE(FTable%KPAIR,dim=2)) =  &
                        FTable%KPAIR(1:SIZE(FTable%KPAIR,dim=1),  1:SIZE(FTable%KPAIR,dim=2))
               dm_POTR  = FTable%POTR
               dm_POTB  = FTable%POTB
               dm_FPOTR = FTable%FPOTR
               dm_FPOTB = FTable%FPOTB

               if(present(RelaseTable)) then
                 if(RelaseTable) call Release_ForceTable(FTable)
               end if


           return
    END SUBROUTINE INITIALIZE_FINNIS_EM_TB_Force_Table_DEV
    !*********************************************************************************

    !*********************************************************************************
    SUBROUTINE CLEAR_FINNIS_EM_TB_Force_Table_DEV()
    !***  PURPOSE:   to deallocate device memory allocated before
    !     INPUT:
    !
    !     OUTPUT
    implicit none
    integer::ERR=0


         if(allocated(dm_DEN)) deallocate(dm_DEN, STAT=ERR)
         if(ERR) goto 100

  !       if(allocated(dm_KPAIR)) deallocate(dm_KPAIR, STAT=ERR)
  !       if(ERR) goto 100

         if(allocated(dm_POTR)) deallocate(dm_POTR, STAT=ERR)
         if(ERR) goto 100

         if(allocated(dm_POTB)) deallocate(dm_POTB, STAT=ERR)
         if(ERR) goto 100

         if(allocated(dm_FPOTR)) deallocate(dm_FPOTR, STAT=ERR)
         if(ERR) goto 100

         if(allocated(dm_FPOTB)) deallocate(dm_FPOTB, STAT=ERR)
         if(ERR) goto 100

         return

  100    write(*,*) "Error in deallocate device memory in INNIS_EM_TB_Force_Table"
         return
    end SUBROUTINE CLEAR_FINNIS_EM_TB_Force_Table_DEV
  !*********************************************************************************



  !*********************************************************************************
   attributes(global) SUBROUTINE PRECALFOR_KERNAL(IM,KVOIS,INDI,ITYP,XP, KT, RMAX, CSI,POTR, POTB, DEN)
  !***  PURPOSE:   to begine calculate the inverse of  the electron density DEN,
  !                 the difference between this and PRECALFOR_KERNAL A is that
  !                 here we assuem the simulation box is cubic, e.g, without consider the
  !                 deformation of box
  !     INPUT:     IM:     the number of particles concerned
  !                KVOIS:  the number of for each atom in style1
  !                INDI:   the index for neighbores
  !                ITYPE:  the type of atom corresponding to INDI
  !                XP:     the coordinate of atoms
  !                KT,     the first dimnesion of lookup table POTR and POTB
  !                RMAX,   the max range of the potential table
  !                CSI,    the interval for interpolation
  !                POTR,   the look-up table for pairwise potential
  !                POTB,   the look-up table for tight-binding potential
  !                DEN,    electron density on atom given by PREFOR_KERNALA

  !     OUTPUT     DEN,   inverse of electron density on atoms

  !
   use MD_CONSTANTS_2010
      implicit none
      !--- DUMMY VARIABLES
      real(KINDDF),value, intent(in)::RMAX, CSI
      integer, value, intent(in)::IM, KT
      integer, device, intent(in)::KVOIS(IM)
      integer, device, intent(in)::INDI(IM,*)
      integer, device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      real(KINDDF), device, intent(in)::POTR(KT,*)
      real(KINDDF), device, intent(in)::POTB(KT,*)
      real(KINDDF), device::DEN(*)

      !Local variables
      integer::J,K,KK, K1,IW, IIW, KTAB, ITYPI, ITYPJ, IFPDX, IFPDY, IFPDZ
      real(KINDDF)::SK,DK,R2,R, DEN0
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ
      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ
      !Local variables: ID for neast neigbors
      real(KINDDF)::NBR2(mp_NEARIST)
      integer::NBID(mp_NEARIST)

      ! Local variables
      integer::IC, IT, IB


              IFPDX = dm_IFPD(1)
              IFPDY = dm_IFPD(2)
              IFPDZ = dm_IFPD(3)

              HBX = dm_HBOX(1)
              HBY = dm_HBOX(2)
              HBZ = dm_HBOX(3)

              BX = dm_BOX(1)
              BY = dm_BOX(2)
              BZ = dm_BOX(3)

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              !POS-- position of the atom
              POSX = XP(IC,1)
              POSY = XP(IC,2)
              POSZ = XP(IC,3)

              if(IC.LT.IM) then
                 !**** for target atoms only
                 !-- start calculation of electron density
                 IIW = KVOIS(IC)
                 ITYPI  = ITYP(IC)
                 DEN0 = C_ZERO
                 DO IW=1, IIW
                    J=INDI(IC,IW)

                    !--- To calculate the seperation between particle IC and its IWth neighbore
                    SEPX = POSX - XP(J,1)
                    SEPY = POSY - XP(J,2)
                    SEPZ = POSZ - XP(J,3)
                    if(IFPDX.GT.0 .AND.(DABS(SEPX) .GT. HBX)) then
                       SEPX = SEPX - DSIGN(BX,SEPX)
                    end if

                    if(IFPDY.GT.0 .AND.(DABS(SEPY) .GT. HBY)) then
                       SEPY = SEPY - DSIGN(BY,SEPY)
                    end if

                    if(IFPDZ.GT.0 .AND.(DABS(SEPZ) .GT. HBZ)) then
                       SEPZ = SEPZ - DSIGN(BZ,SEPZ)
                    end if

                    R2 = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                    !--- To calculate electron density on atom I
                    ITYPJ=ITYP(J)
                    if(R2 .LE. dm_RU2(ITYPI,ITYPJ)) then
                       KTAB = dm_KPAIR(ITYPI,ITYPJ)
                       !SK=R2*CSI
                       !SK = (C_TWO*RMAX*DSQRT(R2)-R2)*CSI
                       !--- Change made 2010-12-05
                       R = DSQRT(R2)
                       SK= DSQRT(R)*CSI

                       KK=SK
                       DK=SK-KK
                       DEN0= DEN0 + POTB(KTAB,KK)+DK*(POTB(KTAB,KK+1)-POTB(KTAB,KK))
                    end if
                 ENDDO
                 DEN0 = DSQRT(DEN0)
                 if(DEN0 .GT. C_ZERO) then
                    DEN(IC) = C_UN/DEN0
                 else
                    DEN(IC) = C_ZERO
                 end if

               !***** for the projectile
              else if(IC.EQ.IM) then

                 !-- start calculation of electron density
                 IIW = KVOIS(IC)
                 ITYPI  = ITYP(IC)
                 DEN0 = C_ZERO

                 NBR2(1:mp_NEARIST) = 1.D64
                 NBID(1:mp_NEARIST) = 0
                 DO IW=1, IIW
                    J=INDI(IC,IW)

                    !--- To calculate the seperation between particle IC and its IWth neighbore
                    SEPX = POSX - XP(J,1)
                    SEPY = POSY - XP(J,2)
                    SEPZ = POSZ - XP(J,3)
                    if(IFPDX.GT.0 .AND.(DABS(SEPX) .GT. HBX)) then
                       SEPX = SEPX - DSIGN(BX,SEPX)
                    end if

                    if(IFPDY.GT.0 .AND.(DABS(SEPY) .GT. HBY)) then
                       SEPY = SEPY - DSIGN(BY,SEPY)
                    end if

                    if(IFPDZ.GT.0 .AND.(DABS(SEPZ) .GT. HBZ)) then
                       SEPZ = SEPZ - DSIGN(BZ,SEPZ)
                    end if

                    R2 = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                    !--- to sort the neighbos
                     DO K=1, mp_NEARIST
                        if(R2.LT.NBR2(K)) then
                           exit
                        end if
                     END DO

                     if(K.LE.mp_NEARIST) then
                        DO K1=mp_NEARIST, K+1,-1
                           NBR2(K1) = NBR2(K1-1)
                           NBID(K1) = NBID(K1-1)
                        END DO
                        NBR2(K) = R2
                        NBID(K) = J
                     end if
                 END DO

                 !--- to calculate e density contributed from nearest neighbors
                 DO IW=1, mp_NEARIST
                    if(NBID(IW) .LT. 1) exit
                    J = NBID(IW)
                    R2 = NBR2(IW)

                    !--- To calculate electron density on atom I
                    ITYPJ=ITYP(J)
                    if(R2 .LE. dm_RU2(ITYPI,ITYPJ)) then
                       KTAB = dm_KPAIR(ITYPI,ITYPJ)
                       R = DSQRT(R2)
                       SK= DSQRT(R)*CSI
                       KK=SK
                       DK=SK-KK
                       DEN0= DEN0 + POTB(KTAB,KK)+DK*(POTB(KTAB,KK+1)-POTB(KTAB,KK))
                    end if
                 ENDDO
                 DEN0 = DSQRT(DEN0)
                 if(DEN0 .GT. C_ZERO) then
                    DEN(IC) = C_UN/DEN0
                 else
                    DEN(IC) = C_ZERO
                 end if

              end if

            !EPOT(IC) = ER0-DEN0
        RETURN
  END SUBROUTINE PRECALFOR_KERNAL
  !*************************************************************************************


  !*************************************************************************************
   attributes(global) SUBROUTINE CALFORCE_KERNAL(IM,KVOIS,INDI,ITYP,XP, KT, RMAX, CSI, FPOTR, FPOTB, DEN,FP)
  !***  PURPOSE:   to calculate the force.
  !                 the difference between this and CALFORCE_KERNALA is that
  !                 here we assuem the simulation box is cubic, e.g, without consider the
  !                 deformation of box
  !
  !     INPUT:     IM:     the number of particles
  !                KVOIS:  the number of neighbores in stle0
  !                INDI:   the index for neighbores
  !                ITYP:   the type of atom corresponding to INDI
  !                XP:     the coordinate of atoms
  !                KT,     the first dimnesion of lookup table FPOTR and FPOTB
  !                RMAX,   the max range of the potential table
  !                CSI,    the interval for interpolation
  !                FPOTR,  the look-up table for pairwise interaction
  !                FPOTB,  the look-up table for tight-binding interaction
  !                DEN,    electron density on atom given by PREFOR_KERNALA
  !
  !     OUTPUT     FP:     the force on atoms
  !
  use MD_CONSTANTS_2010
      implicit none
      !--- DUMMY VARIABLES
      real(KINDDF),value, intent(in)::RMAX, CSI
      integer, value, intent(in)::IM, KT
      integer, device, intent(in)::KVOIS(IM)
      integer, device, intent(in)::INDI(IM,*)
      integer, device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)
      real(KINDDF), device, intent(in)::FPOTR(KT,*)
      real(KINDDF), device, intent(in)::FPOTB(KT,*)
      real(KINDDF), device::DEN(IM)
      real(KINDDF), device::FP(IM,*)

      !Local variables:
      integer::J,K,KK, IW, IIW, KTAB, K1, ITYPI, ITYPJ, IFPDX, IFPDY, IFPDZ
      real(KINDDF)::DENKI,DENKJ,SK,DK, R2,R, FORTOT, SEPX, SEPY, SEPZ, POSX, POSY,POSZ,FPX, FPY, FPZ
      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ
      !Local variables: ID for neast neigbors
      real(KINDDF)::TSEPX(mp_NEARIST),TSEPY(mp_NEARIST),TSEPZ(mp_NEARIST),NBR2(mp_NEARIST)
      integer::NBID(mp_NEARIST)

      ! Local variables
      integer::IC, IT, IB


              IFPDX = dm_IFPD(1)
              IFPDY = dm_IFPD(2)
              IFPDZ = dm_IFPD(3)

              HBX = dm_HBOX(1)
              HBY = dm_HBOX(2)
              HBZ = dm_HBOX(3)

              BX = dm_BOX(1)
              BY = dm_BOX(2)
              BZ = dm_BOX(3)

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              !POS-- position of the atom
              POSX = XP(IC,1)
              POSY = XP(IC,2)
              POSZ = XP(IC,3)
              FPX = C_ZERO
              FPY = C_ZERO
              FPZ = C_ZERO

              if(IC.LT.IM) then
              !*** to calculate the force on target atoms
              !    it is assumed that the last atom is the projectile

               !--- now, we have the electron densities on atoms, to calculate the force
               IIW = KVOIS(IC)
               ITYPI  = ITYP(IC)
               DENKI = DEN(IC)
               !--- to begine scanning the neighbores
               DO IW=1, IIW
                  !--- To calculate the seperation between particle I and its IWth neighbore
                  J=INDI(IC,IW)

                  !--- interaction from projectile accounted later
                  if(J.GE.IM) cycle


                  SEPX = POSX - XP(J,1)
                  SEPY = POSY - XP(J,2)
                  SEPZ = POSZ - XP(J,3)
                  if(IFPDX.GT.0 .AND.(DABS(SEPX) .GT. HBX)) then
                     SEPX = SEPX - DSIGN(BX,SEPX)
                  end if

                  if(IFPDY.GT.0 .AND.(DABS(SEPY) .GT. HBY)) then
                     SEPY = SEPY - DSIGN(BY,SEPY)
                  end if

                  if(IFPDZ.GT.0 .AND.(DABS(SEPZ) .GT. HBZ)) then
                     SEPZ = SEPZ - DSIGN(BZ,SEPZ)
                  end if
                  R2  = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ

                 !--- IF R2 is smaller than cutoff range, to calculate the force
                 ITYPJ=ITYP(J)
                 if(R2 .LE. dm_RU2(ITYPI,ITYPJ)) then
                    KTAB  = dm_KPAIR(ITYPI,ITYPJ)
                    !SK = R2*CSI
                    !SK = (C_TWO*RMAX*DSQRT(R2)-R2)*CSI
                    !--- Change made 2010-12-05
                    R = DSQRT(R2)
                    SK= DSQRT(R)*CSI

                    KK = SK
                    DK=SK-KK

                    DENKJ = DEN(J)
                    FORTOT=(FPOTR(KTAB,KK) + DK*(FPOTR(KTAB,KK+1) - FPOTR(KTAB,KK))) - &
                           (FPOTB(KTAB,KK) + DK*(FPOTB(KTAB,KK+1) - FPOTB(KTAB,KK)))*(DENKI+DENKJ)
                    FPX = FPX + FORTOT*SEPX
                    FPY = FPY + FORTOT*SEPY
                    FPZ = FPZ + FORTOT*SEPZ

                   !-------------------------------------------------
                 end if
               ENDDO
               FP(IC,1) = FPX
               FP(IC,2) = FPY
               FP(IC,3) = FPZ
          end if

          !call syncthreads()
        RETURN
  END SUBROUTINE CALFORCE_KERNAL
  !*************************************************************************************


  !*************************************************************************************
   attributes(global) SUBROUTINE CALFORCE_PROJECTILE_KERNAL(IM,KVOIS,INDI,ITYP,XP, KT, RMAX, CSI, FPOTR, FPOTB, DEN,FP)
  !***  PURPOSE:   to calculate the force between projectile and target atoms.
  !
  !     INPUT:     IM:     the number of particles
  !                KVOIS:  the number of neighbores
  !                INDI:   the index for neighbores
  !                ITYP:   the type of atom corresponding to INDI
  !                XP:     the coordinate of atoms
  !                KT,     the first dimnesion of lookup table FPOTR and FPOTB
  !                RMAX,   the max range of the potential table
  !                CSI,    the interval for interpolation
  !                FPOTR,  the look-up table for pairwise interaction
  !                FPOTB,  the look-up table for tight-binding interaction
  !                DEN,    electron density on atom given by PREFOR_KERNALA
  !
  !     OUTPUT     FP:     the force on atoms
  !
  use MD_CONSTANTS_2010
      implicit none
      !--- DUMMY VARIABLES
      real(KINDDF),value, intent(in)::RMAX, CSI
      integer, value, intent(in)::IM, KT
      integer, device, intent(in)::KVOIS(IM)
      integer, device, intent(in)::INDI(IM,*)
      integer, device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)
      real(KINDDF), device, intent(in)::FPOTR(KT,*)
      real(KINDDF), device, intent(in)::FPOTB(KT,*)
      real(KINDDF), device::DEN(IM)
      real(KINDDF), device::FP(IM,*)

      !Local variables:
      integer::J,K,KK, IW, IIW, KTAB, K1, ITYPI, ITYPJ, IFPDX, IFPDY, IFPDZ
      real(KINDDF)::DENKI,DENKJ,SK,DK, R2,R, FORTOT, SEPX, SEPY, SEPZ, POSX, POSY,POSZ,FPX, FPY, FPZ
      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ
      !Local variables: ID for neast neigbors
      real(KINDDF)::TSEPX(mp_NEARIST),TSEPY(mp_NEARIST),TSEPZ(mp_NEARIST),NBR2(mp_NEARIST)
      integer::NBID(mp_NEARIST)

      ! Local variables
      integer::IC, IT, IB


              IFPDX = dm_IFPD(1)
              IFPDY = dm_IFPD(2)
              IFPDZ = dm_IFPD(3)

              HBX = dm_HBOX(1)
              HBY = dm_HBOX(2)
              HBZ = dm_HBOX(3)

              BX = dm_BOX(1)
              BY = dm_BOX(2)
              BZ = dm_BOX(3)

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !IC -- the id of the atom
              !NOTE: it is assumed that theprojectile atom is thelast one
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC .NE. IM) return

              !POS-- position of the atom
              POSX = XP(IC,1)
              POSY = XP(IC,2)
              POSZ = XP(IC,3)
              FPX = C_ZERO
              FPY = C_ZERO
              FPZ = C_ZERO

              !*** to calculate the force on the projectile atom
              !--- now, we have the electron densities on atoms, to calculate the force
               IIW = KVOIS(IC)
               ITYPI  = ITYP(IC)
               DENKI = DEN(IC)
               NBR2(1:mp_NEARIST) = 1.D64
               NBID(1:mp_NEARIST) = 0
               !--- to begine scanning the neighbores
               DO IW=1, IIW
                  !--- To calculate the seperation between particle I and its IWth neighbore
                  J=INDI(IC,IW)
                  SEPX = POSX - XP(J,1)
                  SEPY = POSY - XP(J,2)
                  SEPZ = POSZ - XP(J,3)
                  if(IFPDX.GT.0 .AND.(DABS(SEPX) .GT. HBX)) then
                     SEPX = SEPX - DSIGN(BX,SEPX)
                  end if

                  if(IFPDY.GT.0 .AND.(DABS(SEPY) .GT. HBY)) then
                     SEPY = SEPY - DSIGN(BY,SEPY)
                  end if

                  if(IFPDZ.GT.0 .AND.(DABS(SEPZ) .GT. HBZ)) then
                     SEPZ = SEPZ - DSIGN(BZ,SEPZ)
                  end if
                  R2  = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ

                  !--- to sort the neighbos
                  DO K=1, mp_NEARIST
                     if(R2.LT.NBR2(K)) then
                        exit
                     end if
                  END DO
                  if(K.LE.mp_NEARIST) then
                     DO K1=mp_NEARIST, K+1, -1
                        NBR2(K1) = NBR2(K1-1)
                        NBID(K1) = NBID(K1-1)
                        TSEPX(K1) = TSEPX(K1-1)
                        TSEPY(K1) = TSEPY(K1-1)
                        TSEPZ(K1) = TSEPZ(K1-1)
                     END DO
                     NBR2(K) = R2
                     TSEPX(K) = SEPX
                     TSEPY(K) = SEPY
                     TSEPZ(K) = SEPZ
                     NBID(K) = J
                  end if
                END DO

                !--- to calculate force from nearest neighbors
                DO IW=1, mp_NEARIST
                   if(NBID(IW) .LT. 1) exit
                   J = NBID(IW)
                   R2 = NBR2(IW)
                   SEPX = TSEPX(IW)
                   SEPY = TSEPY(IW)
                   SEPZ = TSEPZ(IW)

                   ITYPJ=ITYP(J)
                   if(R2 .LE. dm_RU2(ITYPI,ITYPJ)) then
                      KTAB  = dm_KPAIR(ITYPI,ITYPJ)
                      R = DSQRT(R2)
                      SK= DSQRT(R)*CSI
                      KK = SK
                      DK=SK-KK
                      DENKJ = DEN(J)
                      FORTOT=(FPOTR(KTAB,KK) + DK*(FPOTR(KTAB,KK+1) - FPOTR(KTAB,KK))) - &
                             (FPOTB(KTAB,KK) + DK*(FPOTB(KTAB,KK+1) - FPOTB(KTAB,KK)))*(DENKI+DENKJ)
                      FPX = FPX + FORTOT*SEPX
                      FPY = FPY + FORTOT*SEPY
                      FPZ = FPZ + FORTOT*SEPZ

                      FP(J,1) = FP(J,1) - FORTOT*SEPX
                      FP(J,2) = FP(J,2) - FORTOT*SEPY
                      FP(J,3) = FP(J,3) - FORTOT*SEPZ
                  end if
               ENDDO
               FP(IC,1) = FPX
               FP(IC,2) = FPY
               FP(IC,3) = FPZ
            !*********************************************

          !call syncthreads()
        RETURN
  END SUBROUTINE CALFORCE_PROJECTILE_KERNAL
  !*************************************************************************************

  SUBROUTINE CALFORCE_FINNIS_EM_TB_Force_Table2A_DEV(SimBox, CtrlParam, CopyIn, CopyOut, MULTIBOX)
  !***  PURPOSE:   to begine calculate the force.
  !                this routine do not consider the deformation of box
  !                REF:
  !                CALFTENSOR_FINNIS_EM_TB_Force_Table2A_DEV
  !                CALPTENSOR_FINNIS_EM_TB_Force_Table2A_DEV
  !
  !     INPUT:     SimBox,    the simulation box
  !                List,      the neighbore list
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
  !     NOTE:      The neigbore list created by Cal_NeighboreList_DEVKERNAL2C is used
  !                The atom index in NL cretated by Cal_NeighboreList_DEVKERNAL2C is
  !                global index in simulation box. Ref to Cal_NeighboreList_DEVKERNAL2CP
  !                where the atom index in NL is index in a group of neighboring cells
  !
  !     SEE ALSO:  CALFTENSOR_FINNIS_EM_TB_Force_Table2_DEV
  !                Cal_NeighboreList_DEVKERNAL2C
  !
  !
   use MD_Globle_Variables_2010_GPU
   use MD_NeighborsList2C_2010_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),            intent(inout)::SimBox
      type(SimulationCtrlParam), intent(in)::CtrlParam
      integer,                   optional::CopyIn
      integer,                   optional::CopyOut
      integer,                   optional::MULTIBOX

      !--- Local variables
         integer::I, K, K1, NPART, BX, BY, NB, NBX, NBY, KT, ERR
         !--- C2050:
         !integer,parameter::BLOCKSIZE = 512, BLOCKDIMX=32
         !-- C1060:
         integer,parameter::BLOCKSIZE = 64, BLOCKDIMX=64
      !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

         !---
         if(present(CopyIn)) then
             select case(CopyIn)
                case(1)
                   call CopyIn_SimMDBox_XP0_DEV(SimBox)
                case(2,3,4,5,6)
                    call CopyIn_SimMDBox_DEV(SimBox)
                case default
                   call CopyIn_SimMDBox_XP0_DEV(SimBox)
            end select
            if(present(MULTIBOX)) then
               write(*,*) "WARNING: Only one SimMDBox is copied in DEVICE, although MULTIBOX is enabled"
            end if
         end if

        !--- to check the consistent of particle number on host and device
        if(present(MULTIBOX)) then
           NPART = SimBox%NPRT*CtrlParam%MULTIBOX
        else
           NPART = SimBox%NPRT
        end if

        if(dm_NPRT .LT. NPART) then
           write(*,*) "Error:number of particles in device is smaller than Simulation box"
           stop
         end if

        !--- Get the size of the regitered force table
        KT = size(dm_POTR, dim=1)

        !--- If pressure effect acounted for, the box size could be changed, we
        !    need to recopy the boxsize
           dm_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
           dm_BOX(1:3)      = SimBox%ZL(1:3)
           dm_HBOX(1:3)     = SimBox%ZL(1:3)*C_HALF

        !--- to determine size of a block (the number of threads in a block)
            BX = min(BLOCKSIZE, mp_MXSHARESIZE)
            BY = 1
        !-- to determine the dimension of blocks( the number of blocks in a grid)
            NB  = (NPART-1)/(BX*BY)+1
            NBX = min(NB,  BLOCKDIMX)
            NBY = (NB-1)/NBX+1
            NB  = NBX*NBY

            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)

            call PRECALFOR_KERNAL<<<blocks, threads>>>(NPART,dm_KVOIS,dm_INDI, dm_ITYP, dm_XP, KT,m_RMAX, m_CSI,dm_POTR, dm_POTB, dm_DEN)
            err = cudaThreadSynchronize();

            call CALFORCE_KERNAL<<<blocks, threads>>>(NPART,dm_KVOIS,dm_INDI, dm_ITYP, dm_XP, KT, m_RMAX, m_CSI,dm_FPOTR, dm_FPOTB, dm_DEN, dm_FP)
            err = cudaThreadSynchronize();

        !-- to calculate the force between a projectile and targets atoms
            call CALFORCE_PROJECTILE_KERNAL<<<blocks, threads>>>(NPART,dm_KVOIS,dm_INDI, dm_ITYP, dm_XP, KT, m_RMAX, m_CSI,dm_FPOTR, dm_FPOTB, dm_DEN, dm_FP)
            err = cudaThreadSynchronize();


            !--- top copy out the force  if required
            !    NOTE: only the forces in one box is copy out even if MULTIBOX is enabled
            if(present(CopyOut)) then
               if(CopyOut .GT. C_IZERO) then
                  call CopyOut_SimMDBox_Force_DEV(SimBox)

                  if(present(MULTIBOX)) then
                     write(*,*) "WARNING: Only the force of one SimMDBox is copied out DEVICE, although MULTIBOX is enabled"
                  end if
               end if
           end if


          return
  end SUBROUTINE CALFORCE_FINNIS_EM_TB_Force_Table2A_DEV
  !*********************************************************************************


  !*********************************************************************************
   attributes(global) SUBROUTINE CALEPOT_KERNALA(IM,KVOIS,INDI,ITYP,XP, KT, RMAX, CSI,POTR, POTB, EPOT)
  !***  PURPOSE:   to calculate cohesive ennergy
  !     INPUT:     IM:     the number of particles concerned
  !                KVOIS:  the number of for each atom in style1
  !                INDI:   the index for neighbores
  !                ITYPE:  the type of atom corresponding to INDI
  !                XP:     the coordinate of atoms
  !                KT,     the first dimnesion of lookup table POTR and POTB
  !                RMAX,   the max range of the potential table
  !                CSI,    the interval for interpolation
  !                POTR,   the look-up table for pairwise potential
  !                POTB,   the look-up table for tight-binding potential
  !                DEN,    electron density on atom given by PREFOR_KERNALA

  !     OUTPUT     EPOT,  potential on atoms
  !                DEN,   inverse of electron density on atoms

  !
   use MD_CONSTANTS_2010
      implicit none
      !--- DUMMY VARIABLES
      real(KINDDF),value, intent(in)::RMAX, CSI
      integer, value, intent(in)::IM, KT
      integer, device, intent(in)::KVOIS(IM)
      integer, device, intent(in)::INDI(IM,*)
      integer, device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)
      real(KINDDF), device, intent(in)::POTR(KT,*)
      real(KINDDF), device, intent(in)::POTB(KT,*)
      real(KINDDF), device::EPOT(*)

      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, ITYPI, ITYPJ, IFPDX, IFPDY, IFPDZ
      real(KINDDF)::SK,DK,R2,R,ER0, DEN0
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z
      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ, BS11, BS12, BS13, BS21, BS22, BS23, BS31, BS32, BS33

      ! Local variables
      integer::IC, IT, IB


              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC.GT.IM) return

              !POS-- position of the atom
              POSX = XP(IC,1)
              POSY = XP(IC,2)
              POSZ = XP(IC,3)

              IFPDX = dm_IFPD(1)
              IFPDY = dm_IFPD(2)
              IFPDZ = dm_IFPD(3)

              HBX = dm_HBOX(1)
              HBY = dm_HBOX(2)
              HBZ = dm_HBOX(3)

              BX = dm_BOX(1)
              BY = dm_BOX(2)
              BZ = dm_BOX(3)

              BS11 =  dm_BOXSHAPE(1,1)
              BS12 =  dm_BOXSHAPE(1,2)
              BS13 =  dm_BOXSHAPE(1,3)
              BS21 =  dm_BOXSHAPE(2,1)
              BS22 =  dm_BOXSHAPE(2,2)
              BS23 =  dm_BOXSHAPE(2,3)
              BS31 =  dm_BOXSHAPE(3,1)
              BS32 =  dm_BOXSHAPE(3,2)
              BS33 =  dm_BOXSHAPE(3,3)

              !-- start calculation of electron density
              IIW = KVOIS(IC)
              ITYPI  = ITYP(IC)
              ER0 = C_ZERO
              DEN0 = C_ZERO
              DO IW=1, IIW
                 J=INDI(IC,IW)

                 !--- To calculate the seperation between particle IC and its IWth neighbore
                 SEPX = POSX - XP(J,1)
                 SEPY = POSY - XP(J,2)
                 SEPZ = POSZ - XP(J,3)
                 if(IFPDX.GT.0 .AND.(DABS(SEPX) .GT. HBX)) then
                    SEPX = SEPX - DSIGN(BX,SEPX)
                 end if

                 if(IFPDY.GT.0 .AND.(DABS(SEPY) .GT. HBY)) then
                    SEPY = SEPY - DSIGN(BY,SEPY)
                 end if

                 if(IFPDZ.GT.0 .AND.(DABS(SEPZ) .GT. HBZ)) then
                    SEPZ = SEPZ - DSIGN(BZ,SEPZ)
                 end if

                 !--- NOTE: To converte the length-unit into absolute unit (cm)
                 !          only when the shape of the box is cubic,
                 !          DXYZ = SEP
                 DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                 DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                 DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                 R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                 !--- To calculate electron density on atom I
                 ITYPJ=ITYP(J)
                 if(R2 .LE. dm_RU2(ITYPI,ITYPJ)) then
                    KTAB = dm_KPAIR(ITYPI,ITYPJ)
                    !SK=R2*CSI
                    !SK = (C_TWO*RMAX*DSQRT(R2)-R2)*CSI
                    !--- Change made 2010-12-05
                    R = DSQRT(R2)
                    SK= DSQRT(R)*CSI

                    KK=SK
                    DK=SK-KK
                    ER0 = ER0  + POTR(KTAB,KK)+DK*(POTR(KTAB,KK+1)-POTR(KTAB,KK))
                    DEN0= DEN0 + POTB(KTAB,KK)+DK*(POTB(KTAB,KK+1)-POTB(KTAB,KK))
                end if
            ENDDO

            EPOT(IC) = ER0-DSQRT(DEN0)

        RETURN
  END SUBROUTINE CALEPOT_KERNALA
  !*************************************************************************************



  !*********************************************************************************
  SUBROUTINE CALEPOT_FINNIS_EM_TB_Force_Table2A_DEV(SimBox, CtrlParam, CopyIn, CopyOut, MULTIBOX)
  !***  PURPOSE:   to begine calculate the force. the box shape could be uncubic.
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                List,      the neighbore list
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
  !     NOTE:      The neigbore list created by Cal_NeighboreList_DEVKERNAL2C is used
  !                The atom index in NL cretated by Cal_NeighboreList_DEVKERNAL2C is
  !                global index in simulation box. Ref to Cal_NeighboreList_DEVKERNAL2CP
  !                where the atom index in NL is index in a group of neighboring cells
  !
  !     SEE ALSO:  CALFORCE_FINNIS_EM_TB_Force_Table2A_DEV
  !                CALFTENSOR_FINNIS_EM_TB_Force_Table2A_DEV
  !                Cal_NeighboreList_DEVKERNAL2C
  !
   use MD_Globle_Variables_2010_GPU
   use MD_NeighborsList2C_2010_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),            intent(inout)::SimBox
      type(SimulationCtrlParam), intent(in)   ::CtrlParam
      integer,                   optional     ::CopyIn
      integer,                   optional     ::CopyOut
      integer,                   optional     ::MULTIBOX

      !--- Local variables
         integer I, K, K1, NPART, BX, BY, NB, NBX, NBY, KT, ERR
         !--- C2050
         !integer,parameter::BLOCKSIZE = 512, BLOCKDIMX=32
         !--- C1060
         integer,parameter::BLOCKSIZE = 64, BLOCKDIMX=64
      !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

         !---
         if(present(CopyIn)) then
             select case(CopyIn)
                case(1)
                   call CopyIn_SimMDBox_XP0_DEV(SimBox)
                case(2,3,4,5,6)
                    call CopyIn_SimMDBox_DEV(SimBox)
                case default
                   call CopyIn_SimMDBox_XP0_DEV(SimBox)
            end select

            if(present(MULTIBOX)) then
               write(*,*) "WARNING: Only one SimMDBox is copied in DEVICE, although MULTIBOX is enabled"
            end if
         end if


        !--- to check the consistent of particle number on host and device
        if(present(MULTIBOX)) then
           NPART = SimBox%NPRT*CtrlParam%MULTIBOX
        else
           NPART = SimBox%NPRT
        end if

        if(dm_NPRT .LT. NPART) then
           write(*,*) "Error:number of particles in device is smaller than Simulation box"
           stop
         end if

        !--- Get the size of the regitered force table
        KT = size(dm_POTR, dim=1)

        !--- If pressure effect acounted for, the box size could be changed, we
        !    need to recopy the boxsize
           dm_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
           dm_BOX(1:3)      = SimBox%ZL(1:3)
           dm_HBOX(1:3)     = SimBox%ZL(1:3)*C_HALF

        !--- to determine size of a block (the number of threads in a block)
            BX = min(BLOCKSIZE, mp_MXSHARESIZE)
            BY = 1
        !-- to determine the dimension of blocks( the number of blocks in a grid)
            NB  = (NPART-1)/(BX*BY)+1
            NBX = min(NB,  BLOCKDIMX)
            NBY = (NB-1)/NBX+1
            NB  = NBX*NBY

            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)

            call CALEPOT_KERNALA<<<blocks, threads>>>(NPART,dm_KVOIS,dm_INDI, dm_ITYP, dm_XP, KT, m_RMAX, m_CSI,dm_POTR, dm_POTB, dm_EPOT)
            err = cudaThreadSynchronize();

            !--- top copy out the EPOT  if required
            if(present(CopyOut)) then
               if(CopyOut .GT. C_IZERO) then
                  call CopyOut_SimMDBox_EPOT_DEV(SimBox)

                  if(present(MULTIBOX)) then
                     write(*,*) "WARNING: Only the EPOT of one SimMDBox is copied out DEVICE, although MULTIBOX is enabled"
                  end if
               end if
           end if


          return
  end SUBROUTINE CALEPOT_FINNIS_EM_TB_Force_Table2A_DEV
  !*********************************************************************************


  end module BAC_TEST_TB_ForceTable_GPU





