  module MD_EAM_Force_Table_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  To calculate the forces on atoms USING GPU paralle comuputing, the force is assumed in
  !                  EAM: Ec=-F(RO)+SUMij(Vij) with RO the total electron density on atoms; Vij,the pairwise
  !                  portential on atom i by atom j.
  !
  !                  The force calculations are performed by table-lookup methods with the FORCETABLE previuosly
  !                  calculated. A FORCETABLE should be registered before calling this module.
  !                  Also should be noted is that in paralle calculations, we cannot apply Newton's
  !                  third law. Be careful to use correct neighbore list(ref to MD_NeighborsList2C_12_GPU.cuf)
  !
  !
  !                  SEE ALSO______________________________________________________________________________
  !                       MDTypeDef_ForceTable.F90
  !                       MD_FS_Force_Table_GPU.F90
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.f90
  !                       MD_Globle_Variables_GPU
  !
  !                  ______________________________________________________________________________________
  !                  HOU Qing, April, 2012
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:     The version is addapted from  MD_FS_Force_Table_GPU
  !                  2012-04 (Hou Qing)
  !
  !
  !                  2017-10-19,  statu checking of atoms are added. Only the force of activative
  !                               atoms are calculated (Hou Qing).
  !
  !                  2018-09-17,  in COPYOUT_VIRIALTENSOR, the virial tensor copied out is the
  !                               total virial tensor divided by the number of boxes (Hou Qing)
  ! 
  !                  2018-12-31,  remove the constant memory dcm_BOXSHAPE, and dcm_BOXSIZE, and thus
  !                               thus the copyin operations of the constant memory in PRECALFOR_template .
  !                               The boxsize and shape are passed to kernel by BX0,BY0,....
  !                               Based on DDT-MAP, This modification maybe profit for acceleratting the 
  !                               force calculations.
  !     
  !                  2018-01-02,  copyout of DEN on device that was performed in PRECALFOR_template are 
  !                               replaced by Synchroniz_DEN_on_Devices, called outside PRECALFOR_template.
  !     
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_ForceTable
    use MD_MultiGPU_Basic
  !-----------------------------------------------
    implicit none

      logical, private::m_INITED = .false.
      !--- C2050
      !integer, parameter, private::mp_MXSHARESIZE = 3*3*224
      !integer,parameter::mp_BLOCKSIZE = 512, mP_BLOCKDIMX=32
      !--- C1060
      integer, parameter, private::mp_MXSHARESIZE = 3*3*64
      integer, parameter, private::mp_BLOCKSIZE = 256, mp_BLOCKDIMX=128

      real(KINDDF)  ::m_CSI                                            ! the inverse of RmaxSQRT e.g. NTAB/RmaxSQRT
      real(KINDDF)  ::m_CSIV                                           ! the inverse of CSI,        1/CSI
      integer       ::m_KPAIRDIM(2)                                    ! the dimension of KPAIR
      integer       ::m_TABDIM(2)                                      ! the dimension of the force table
      integer       ::m_KPAIR(mp_MXGROUP,mp_MXGROUP)                   ! the index of table for atom-atom interaction
      real(KINDDF)  ::m_RU2(mp_MXGROUP,mp_MXGROUP)                     ! the squqre of cutoff disances
      real(KINDDF)  ::m_mxRU2                                          ! the largest value of m_RU2

      real(KINDDF)  ::m_RHOD                                           ! the interval of RO for embedment function table
      integer       ::m_KEMBDDIM                                       ! the dimension of KEMBD
      integer       ::m_TABEMDDIM(2)                                   ! the dimension of the embedment function table
      integer       ::m_KEMBD(mp_MXGROUP)                              ! the index of table for embedment function of atom1

      !--- box information
      real(KINDDF), private::m_BOXSHAPE(3,3)                           ! siumulation box shape
      real(KINDDF), private::m_BOXSIZE(3)                              ! simulation box size
      integer,      private::m_IFPD(3)                                 ! the periodic condition

      !$$--- box and cuoff information on devices
      integer,      constant, private::dcm_IFPD(3)
      integer,      constant, private::dcm_KPAIR(mp_MXGROUP,mp_MXGROUP) ! the index of table for atom-atom interaction

      !$$real(KINDDF),constant, private::dcm_RU2(mp_MXGROUP,mp_MXGROUP)! the cutoff range
      !$$ because the limitation of constant memory,
      !$$ dcm_RU2 is replaced by m_mxRU2 by HOU Qing Dec 22, 2017
      integer,     constant, private::dcm_KEMBD(mp_MXGROUP)            ! the index of table for embedment function of atom1

      real(KINDDF), dimension(:), allocatable, private::hm_DEN ! the embedded term (e-density),the subscript is the index in the whole box

      private::FTableDEVWS
      type::FTableDEVWS
            real(KINDDF), device, dimension(:,:),  allocatable::POTR       ! the potential table for Nucl.-Nucl.  interaction with size of  NKIND*NTAB
            real(KINDDF), device, dimension(:,:),  allocatable::FPOTR      ! the differential of potential table POTR
            real(KINDDF), device, dimension(:,:),  allocatable::POTB       ! the potential table
            real(KINDDF), device, dimension(:,:),  allocatable::FPOTB      ! the table of differential of RHO on r
            real(KINDDF), device, dimension(:,:),  allocatable::FEMBD      ! the table of emebedment function
            real(KINDDF), device, dimension(:,:),  allocatable::DFEMBD     ! the table of differential of emebedment function
            real(KINDDF), device, dimension(:),    allocatable::DEN        ! the electron density when calculate EPOT,
            real(KINDDF), device, dimension(:,:,:),allocatable::VTENSOR    !the virial tensor
      end type FTableDEVWS
      type(FTableDEVWS), dimension(m_MXDEVICE), private::dm_FTableWS

      !****
      private::&
               Allocate_ForceTable_template,  &
               CALEPOT_template,              &
               CALFORCE_template,             &
               CALPTENSOR_template,           &
               CLEAR_Force_Table_template,    &
               PRECALFOR_template,            &
               Synchroniz_DEN_on_Devices

      private::&
               CALEPOT_KERNEL,    &
               CALFORCE_KERNEL,   &
               CALPTENSOR_KERNEL, &
               PRECALFOR_KERNEL
               
  contains

  !*********************************************************************************
    subroutine Allocate_ForceTable_template (IDEV, FTable, POTR,FPOTR, POTB, FPOTB, FEMBD, DFEMBD, DEN, VTENSOR)
    !***  PURPOSE:   to allocate device memory for force calculation
    !     INPUT:     IDEV,   the ID of device
    !                FTable, the force tabled to be used
    !
    !     OUTPUT:    POTR,FPOTR, POTB, FPOTB, DEN, VTENSOR
    !                also initialize the constant memories: dcm_KAIR, dcm_IFPD
    !
    use MD_Globle_Variables_GPU
    implicit none
       !--- dummy vaiables
       integer::IDEV
       type(MDForceTable),intent(in)::FTable
       real(KINDDF), device, dimension(:,:),allocatable::POTR
       real(KINDDF), device, dimension(:,:),allocatable::FPOTR
       real(KINDDF), device, dimension(:,:),allocatable::POTB
       real(KINDDF), device, dimension(:,:),allocatable::FPOTB
       real(KINDDF), device, dimension(:,:),allocatable::FEMBD
       real(KINDDF), device, dimension(:,:),allocatable::DFEMBD
       real(KINDDF), device, dimension(:),  allocatable::DEN
       real(KINDDF), device, dimension(:,:,:), allocatable::VTENSOR

      !--- Local variables
      integer::ERR, CURDEV, BX, BY, NBX, NBY, NB

               ERR = cudaGetDevice(CURDEV)
               ERR = cudaSetDevice(IDEV)

           !$$--- allocate memory for electron density and the virial terms
           !$$    NOTE: the subsrcipt of DEN is the index of partciels in the whole box
           !$$
               allocate( DEN(dm_NPRT),VTENSOR(3,3,mp_BLOCKDIMX), STAT=ERR)
               if(ERR) then
                  write(*,fmt="(A, I3)") "MDPSCU Error: fail to allocate memory for DEN in MD_EAM_ForceTable on device", IDEV
                  write(*,fmt="(A, I)") "              number of atoms is ", dm_NPRT
                  stop "1"
               end if
               DEN = C_ZERO
               VTENSOR = C_ZERO

               dcm_KPAIR = m_KPAIR
               dcm_KEMBD = m_KEMBD
               dcm_IFPD  = m_IFPD

               allocate(POTR(m_TABDIM(1), m_TABDIM(2)), POTB(m_TABDIM(1), m_TABDIM(2)),              &
                        FPOTR(m_TABDIM(1),m_TABDIM(2)), FPOTB(m_TABDIM(1),m_TABDIM(2)),              &
                        FEMBD(m_TABEMDDIM(1),m_TABEMDDIM(2)), DFEMBD(m_TABEMDDIM(1),m_TABEMDDIM(2)), &
                        STAT=ERR  )

               if(ERR) then
                  write(*,*) "MDPSCU Error: fail to allocate device memory for force tables in MD_EAM_ForceTable", IDEV
                  stop "2"
               end if

               POTR   = FTable%POTR
               POTB   = FTable%POTB
               FPOTR  = FTable%FPOTR
               FPOTB  = FTable%FPOTB
               FEMBD  = FTable%FEMBD
               DFEMBD = FTable%DFEMBD

               ERR = cudaSetDevice(CURDEV)

           return
    end subroutine Allocate_ForceTable_template
  !*********************************************************************************

  !*********************************************************************************
    subroutine INITIALIZE_EAM_Force_Table_DEV(SimBox, CtrlParam, FTable,RelaseTable, MULTIBOX)
    !***  PURPOSE:   to intialize the parameters to be used in force calculation
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !                FTable, the force tabled to be used
    !                RelaseTable, indicator to indicate if the force table will be deallocate after
    !                             copyin devices
    !                MULTIBOX, indicating if use multiple box
    !
    !     OUTPUT:   the working spaces allocated
    !
     use MD_Globle_Variables_GPU
     use MiniUtilities, only:ONWARNING
     implicit none
       !--- dummy vaiables
      type(SimMDBox),     intent(inout)::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(in)   ::FTable
      integer,            optional     ::RelaseTable, MULTIBOX
      !--- Local variables
      integer::ERR, NPRT, I
      logical::pinnedFlag

               call Clear_EAM_Force_Table_DEV()

           !$$--- allocate memory for electron density
               if(present(MULTIBOX)) then
                  if(MULTIBOX.gt. 0) then
                     NPRT = SimBox%NPRT*CtrlParam%MULTIBOX
                  else
                     NPRT = SimBox%NPRT
                  end if
               else
                  NPRT = SimBox%NPRT
               end if

               if(dm_NPRT .LT. NPRT) then
                  write(*,*) "MDPSCU Error:number of particles in device is smaller than Simulation box", dm_NPRT,NPRT
                  stop
               end if

           !$$--- to copy force table information
               m_CSI         = FTable%CSI
               m_CSIV        = FTable%CSIV
               m_KPAIRDIM(1) = SIZE(FTable%KPAIR,dim=1)
               m_KPAIRDIM(2) = SIZE(FTable%KPAIR,dim=2)
               m_TABDIM(1)   = SIZE(FTable%POTR,dim=1)
               m_TABDIM(2)   = SIZE(FTable%POTR,dim=2)

               m_RHOD        = FTable%RHOD
               m_KEMBDDIM    = SIZE(FTable%KEMBD)
               m_TABEMDDIM(1)= SIZE(FTable%FEMBD, dim=1)
               m_TABEMDDIM(2)= SIZE(FTable%FEMBD, dim=2)

               allocate(hm_DEN(NPRT), STAT=ERR) !, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of DEN failed in INITIALIZE_EAM_Force_Table'
                  write(*,fmt="(A)") '               Process to be stopped'
                  stop
               !else
               !   if (.not. pinnedFlag) then
               !       write(*,fmt="(A)") ' MDPSCU Warnning; Pinned allocation failed of DEN'
               !       call ONWARNING(gm_OnWarning)
               !   end if
               end if

               m_KPAIR(1:m_KPAIRDIM(1),1:m_KPAIRDIM(2)) =  FTable%KPAIR(1:m_KPAIRDIM(1), 1:m_KPAIRDIM(2))
               m_KEMBD(1:m_KEMBDDIM) = FTable%KEMBD(1:m_KEMBDDIM)

               m_RU2   = CtrlParam%RU*CtrlParam%RU
               m_mxRU2 = maxval(m_RU2)
               m_IFPD  = CtrlParam%IFPD
           !$$---  create the table on devices
               do I=1,  m_NDEVICE
                 call Allocate_ForceTable_template(m_DEVICES(I), FTable,                                                &
                                                   dm_FTableWS(I)%POTR,  dm_FTableWS(I)%FPOTR,   dm_FTableWS(I)%POTB,   &
                                                   dm_FTableWS(I)%FPOTB, dm_FTableWS(I)%FEMBD,   dm_FTableWS(I)%DFEMBD, &
                                                   dm_FTableWS(I)%DEN,   dm_FTableWS(I)%VTENSOR)
               end do
               if(present(RelaseTable)) then
                 if(RelaseTable) call Release_ForceTable(FTable, keepreglist=1)
               end if

               m_INITED = .TRUE.

           return
    end subroutine INITIALIZE_EAM_Force_Table_DEV
    !*********************************************************************************

    !*********************************************************************************
    subroutine Clear_Force_Table_template(IDEV, POTR,FPOTR, POTB, FPOTB, FEMBD, DFEMBD, DEN,VTENSOR)
    !***  PURPOSE:   to deallocate memories allocated on a device
    !     INPUT:     IDEV, the ID of the device
    !
    !     OUTPUT:    POTR,FPOTR, POTB, FPOTB, DEN,VTENSOR
    implicit none
       !--- dummy vaiables
       integer::IDEV

       real(KINDDF), device, dimension(:,:),   allocatable::POTR
       real(KINDDF), device, dimension(:,:),   allocatable::FPOTR
       real(KINDDF), device, dimension(:,:),   allocatable::POTB
       real(KINDDF), device, dimension(:,:),   allocatable::FPOTB
       real(KINDDF), device, dimension(:,:),   allocatable::FEMBD
       real(KINDDF), device, dimension(:,:),   allocatable::DFEMBD
       real(KINDDF), device, dimension(:),     allocatable::DEN
       real(KINDDF), device, dimension(:,:,:), allocatable::VTENSOR

       !--- local vairbales
      integer::ERR, CURDEV

              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

             if(allocated(DEN)) deallocate(DEN, STAT=ERR)
             if(ERR) goto 100

             if(allocated(POTR)) deallocate(POTR, STAT=ERR)
             if(ERR) goto 100

             if(allocated(POTB)) deallocate(POTB, STAT=ERR)
             if(ERR) goto 100

             if(allocated(FPOTR)) deallocate(FPOTR, STAT=ERR)
             if(ERR) goto 100

             if(allocated(FPOTB)) deallocate(FPOTB, STAT=ERR)
             if(ERR) goto 100

             if(allocated(FEMBD)) deallocate(FEMBD, STAT=ERR)
             if(ERR) goto 100

             if(allocated(DFEMBD)) deallocate(DFEMBD, STAT=ERR)
             if(ERR) goto 100

             if(allocated(VTENSOR)) deallocate(VTENSOR, STAT=ERR)
             if(ERR) goto 100


             ERR = cudaSetDevice(CURDEV)

             m_INITED = .FALSE.
         return

  100    write(*,*) "MDPSCU Error in deallocate device memory in FINNIS_EM_TB_Force_Table"
         return
    end subroutine Clear_Force_Table_template
   !*********************************************************************************


   !*********************************************************************************
    subroutine Clear_EAM_Force_Table_DEV()
    !***  PURPOSE:   to deallocate memories allocated before
    !     INPUT:
    !
    !     OUTPUT
    implicit none
    integer::I

           !---
              if(allocated(hm_DEN)) then
                 deallocate(hm_DEN)
              end if
           !$$---  clear the table on devices
              do I=1, m_NDEVICE 
                 call Clear_Force_Table_template(m_DEVICES(I), dm_FTableWS(I)%POTR,  dm_FTableWS(I)%FPOTR, dm_FTableWS(I)%POTB,   &
                                                               dm_FTableWS(I)%FPOTB, dm_FTableWS(I)%FEMBD, dm_FTableWS(I)%DFEMBD, &
                                                               dm_FTableWS(I)%DEN,   dm_FTableWS(I)%VTENSOR)
              end do
         return

  100    write(*,*) "MDPSCU Error in deallocate device memory in FINNIS_EM_TB_Force_Table"
         return
    end subroutine Clear_EAM_Force_Table_DEV
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine PRECALFOR_KERNEL(IM, ITYP,XP, NAPDEV, STATU, IA0, NPRT, KVOIS,INDI, CSI, RU2, KT, POTB,       &
                                                  RHOD, KT1, DFEMBD, DEN,                                                      &                                        
                                                  BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033 )
  !***  PURPOSE:   to calculate the inverse of  the electron density DEN
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              ITYP :  the type of atoms corresponding to INDI
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              STATU:  the statu of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              NPRT:   the actual number of atoms on the device
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              CSI,    the interval for interpolation
  !$$              RU2,    the square of cutoff distance
  !$$              KT,     the first dimension of lookup table POTR and POTB
  !$$              POTB,   the look-up table for electron density
  !$$              RHOD,   the RHO interval in emebedment function table
  !$$              KT1,    the first dimension of FEMBD and DFEMBD
  !$$              DFEMBD, the look-up table for differential of embedment function
  !
  !$$   OUTPUT:    DEN,   dF(RHO)/d(RHO)
  !

      implicit none
      !--- DUMMY VARIABLES
      integer, value, intent(in)::IM
      integer, device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer, value,  intent(in)::NAPDEV, IA0, NPRT
      integer, device, intent(in)::STATU(NAPDEV)
      integer, device, intent(in)::KVOIS(NAPDEV)
      integer, device, intent(in)::INDI(NAPDEV,*)

      real(KINDDF), value,  intent(in)::CSI, RU2
      integer,      value,  intent(in)::KT
      real(KINDDF), device, intent(in)::POTB(KT,*)

      real(KINDDF), value, intent(in)::RHOD
      integer,      value, intent(in)::KT1
      real(KINDDF), device, intent(in)::DFEMBD(KT1,*)

      real(KINDDF), device::DEN(*)

      real(KINDDF), value::BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033
      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, ITYPI, ITYPJ
      real(KINDDF)::DEN0, SK, R2,R
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z

      real(KINDDF),shared::HBX, HBY, HBZ, BX, BY, BZ, BS11, BS12, BS13, BS21, BS22, BS23, BS31, BS32, BS33
      integer,     shared::IFPDX, IFPDY, IFPDZ
      integer,     shared::NT, NB

      ! Local variables
      integer::IC, IT, IB, LL, NAL,NL

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              !if(IT .EQ. 1) then
                 BX = BX0 !dcm_BOXSIZE(1)
                 BY = BY0 !dcm_BOXSIZE(2)
                 BZ = BZ0 !dcm_BOXSIZE(3)

                 HBX = BX*C_HALF
                 HBY = BY*C_HALF
                 HBZ = BZ*C_HALF

                 BS11 =  BS011 !dcm_BOXSHAPE(1,1)
                 BS12 =  BS012 !dcm_BOXSHAPE(1,2)
                 BS13 =  BS013 !dcm_BOXSHAPE(1,3)
                 BS21 =  BS021 !dcm_BOXSHAPE(2,1)
                 BS22 =  BS022 !dcm_BOXSHAPE(2,2)
                 BS23 =  BS023 !dcm_BOXSHAPE(2,3)
                 BS31 =  BS031 !dcm_BOXSHAPE(3,1)
                 BS32 =  BS032 !dcm_BOXSHAPE(3,2)
                 BS33 =  BS033 !dcm_BOXSHAPE(3,3)

                 IFPDX = dcm_IFPD(1)
                 IFPDY = dcm_IFPD(2)
                 IFPDZ = dcm_IFPD(3)

                 !$$--- size of Block
                 NT = blockdim%x*blockdim%y

                 !$$--- size of grid
                 NB = griddim%x*griddim%y
            !end if
            !call syncthreads()

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1

            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
            DO LL=1, NL

              !$$IC -- the id of the atom on the device
              !$$      NOTE: IC is not the id of the same atom in the whole box
              IC= (IB-1)*NT+IT +(LL-1)*NAL

              if(IC.LE.NPRT) then
                 !$$--- NOTE by HOU Qing: Oct 10,2017
                 !$$         The activation examination added. only the activative atoms
                 !$$         need to calculate force
                 DEN0 = C_ZERO
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then

                    !$$POS-- position of the atom
                    !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                    POSX = XP(IC+IA0,1)
                    POSY = XP(IC+IA0,2)
                    POSZ = XP(IC+IA0,3)
                    ITYPI  = ITYP(IC+IA0)

                    !$$-- start calculation of electron density
                    IIW = KVOIS(IC)
                    do IW=1, IIW
                       !$$--- NOTE: the particles index of neighbore-list is
                       !$$          the index of particle in the whole box
                       J=INDI(IC,IW)

                       !$$--- To calculate the seperation between particle IC
                       !$$    and its IWth neighbore
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

                       !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                       !$$          only when the shape of the box is cubic,
                       !$$          DXYZ = SEP
                       DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                       DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                       DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                       R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                       !$$--- To calculate electron density on atom I
                       ITYPJ=ITYP(J)
                       !$$if(R2 .LE. dcm_RU2(ITYPI,ITYPJ)) then
                       if(R2 .LE. RU2) then
                          KTAB =  dcm_KPAIR(ITYPI,ITYPJ)
                          R = DSQRT(R2)
                          SK= DSQRT(R)*CSI
                          KK=SK
                          !$$--- NOTE by HOU Qing: Dec 4,2014
                          !$$          The force table POTR, POTB are V(r), RHO(r)*0.5, in OLD VERSION
                          !$$          In the new version. we define POTR, POTB as V(r)*r and RHO(r)
                          DEN0= DEN0 + (POTB(KTAB,KK)+(SK-KK)*(POTB(KTAB,KK+1)-POTB(KTAB,KK)))
                       end if
                    end do

                    !$$**** obtaine dF(RHO)/d(RHO)  by interpolation
                    if(DEN0 .GT. C_ZERO) then
                       KTAB =  dcm_KEMBD(ITYPI)
                       !the range of KK is from 1 to KTABEMBDDIM(2)
                       SK   =  DEN0/RHOD + 1.D0
                       KK   =  SK + 0.000001D0
                       DEN0 =  DFEMBD(KTAB, KK) + (SK-KK)*(DFEMBD(KTAB, KK+1) -  DFEMBD(KTAB, KK) )
                    end if
                 end if
                 !$$--- NOTE: the subscript of DEN is the indice of atoms in the whole box
                 DEN(IC+IA0) = DEN0
              end if
            END DO !End loop for LL
        RETURN
  end subroutine PRECALFOR_KERNEL
  !*********************************************************************************

  !*********************************************************************************
  subroutine PRECALFOR_template(IDEV, STARTCELL, ENDCELL, KVOIS, INDI, ITYP, XP, STATU, POTB,DFEMBD,DEN)
  !***  PURPOSE:   to calculate Finnis-Sinclar term.
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL, the ID of the first cell on the device
  !                ENDCELL,   the ID of the last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                STATU:  the statu of atoms
  !                POTB,   the look-up table for tight-binding potential
  !                DFEMBD, the look-up table for differential of embedment function
  !
  !     OUTPUT:   DEN,     the Finnis-Sinlar term
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      integer,      device, dimension(:)  ::KVOIS
      integer,      device, dimension(:,:)::INDI
      integer,      device, dimension(:)  ::ITYP
      real(KINDDF), device, dimension(:,:)::XP
      integer,      device, dimension(:)  ::STATU
      real(KINDDF), device, dimension(:,:)::POTB
      real(KINDDF), device, dimension(:,:)::DFEMBD
      real(KINDDF), device, dimension(:)::DEN

      !--- Local variables
         integer::ERR, STARTA, IA0, ENDA, NPRT

      !$$--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

             ERR = cudaSetDevice(IDEV)
             STARTA   = hm_IA1th(STARTCELL)
             ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
             NPRT = ENDA - STARTA + 1

             blocks  = dim3(mp_BLOCKDIMX, 1, 1)
             threads = dim3(mp_BLOCKSIZE, 1, 1)
             IA0     = STARTA-1

             call PRECALFOR_KERNEL<<<blocks, threads>>>(dm_NPRT,ITYP,XP, m_NAPDEV, STATU, IA0, NPRT, KVOIS,INDI,           &
                                                        m_CSI, m_MXRU2, m_TABDIM(1), POTB, m_RHOD, m_TABEMDDIM(1), DFEMBD, &
                                                        DEN,                                                               &
                                                        m_BOXSIZE(1),    m_BOXSIZE(2),    m_BOXSIZE(3),                    &
                                                        m_BOXSHAPE(1,1), m_BOXSHAPE(1,2), m_BOXSHAPE(1,3),                 &
                                                        m_BOXSHAPE(2,1), m_BOXSHAPE(2,2), m_BOXSHAPE(2,3),                 &
                                                        m_BOXSHAPE(3,1), m_BOXSHAPE(3,2), m_BOXSHAPE(3,3)    )

            !-- but we need to copyin DEN, because the subscript of
            !   DEN is the index pf particles in the whole box, just like XP
            !   we need copy the DEN on device into host
            !   ERR = cudaMemcpyAsync(hm_DEN(STARTA),DEN(STARTA),NPRT)
            !   This operatione is performed by call Synchroniz_DEN_on_Devices(), /2019-01-02 HOU Qing                                                        
          return
  end subroutine PRECALFOR_template
  !*********************************************************************************

  !****************************************************************************
  subroutine Synchroniz_DEN_on_Devices()
   !***  PURPOSE:  because different devices handled different segment of XP.
   !                this routine is to make XP on all devices to be same
   !
       use MD_Globle_Variables_GPU
       implicit none
       !--- dummy variables
       !----   Local variables
       integer::I, ERR, CURDEV

            if(m_NDEVICE .gt. 1) then
               ERR = cudaGetDevice(CURDEV)
               do I=1, m_NDEVICE
                  ERR = cudaSetDevice(m_DEVICES(I))
                  ERR = cudaMemcpyAsync(hm_DEN(m_STARTA(I)),dm_FTableWS(I)%DEN(m_STARTA(I)), m_NPA(I))               
               enddo   
               call SynchronizeDevices()
               do I=1, m_NDEVICE
                  ERR = cudaSetDevice(m_DEVICES(I))
                  ERR = cudaMemcpyAsync(dm_FtableWS(I)%DEN, hm_DEN, dm_NPRT)
               end do
               ERR = cudaSetDevice(CURDEV)
            endif   

          return
    end subroutine Synchroniz_DEN_on_Devices
  !****************************************************************************  
  
  !*********************************************************************************
   attributes(global) subroutine CALFORCE_KERNEL(IM, ITYP,XP, NAPDEV, STATU, IA0, NPRT,KVOIS,INDI, &
                                                 CSI, RU2, KT, FPOTR, FPOTB, DEN,FP,               &
                                                 BX0, BY0, BZ0)
  !***  PURPOSE:   to  calculate the forces on atoms
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              ITYP :  the type of atoms corresponding to INDI
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              STATU:  the statu of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              NPRT:   the actual number of atoms on the device
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              CSI,    the interval for interpolation
  !$$              RU2,    the square of cutoff distance
  !$$              KT,     the first dimension of lookup table POTR and POTB
  !$$              FPOTR,  the look-up table for pairwise interaction
  !$$              FPOTB,  the look-up table for tight-binding interaction
  !$$              DEN,    the Sinclar-Finnis term
  !
  !$$   OUTPUT:    FP,     the force on atoms
  !

      implicit none
      !--- DUMMY VARIABLES
      integer,      value, intent(in)::IM
      integer,      device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer,      value,  intent(in)::NAPDEV, IA0,NPRT
      integer,      device, intent(in)::STATU(NAPDEV)
      integer,      device, intent(in)::KVOIS(NAPDEV)
      integer,      device, intent(in)::INDI(NAPDEV,*)

      real(KINDDF), value,  intent(in)::CSI, RU2
      integer,      value,  intent(in)::KT
      real(KINDDF), device, intent(in)::FPOTR(KT,*)
      real(KINDDF), device, intent(in)::FPOTB(KT,*)
      real(KINDDF), device, intent(in)::DEN(IM)
      real(KINDDF), value             ::BX0, BY0, BZ0
      !Output:
      real(KINDDF), device, intent(out)::FP(NAPDEV,*)

      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, KTAB1, ITYPI, ITYPJ
      real(KINDDF)::FORTOT,FPX,FPY,FPZ, SK, DK, R2, R ,EXP5, EXP6, EXP7, DENKI, DENKJ
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ

      real(KINDDF), shared::HBX, HBY, HBZ, BX, BY, BZ
      integer,      shared::IFPDX, IFPDY, IFPDZ
      integer,      shared::NT, NB
      ! Local variables
      integer::IC, IT, IB, LL, NAL,NL

            IT  = (threadidx%y-1)*blockdim%x + threadidx%x

            !if(IT .EQ. 1) then
               BX = BX0  !dcm_BOXSIZE(1)
               BY = BY0  !dcm_BOXSIZE(2)
               BZ = BZ0  !dcm_BOXSIZE(3)

               HBX = BX*C_HALF
               HBY = BY*C_HALF
               HBZ = BZ*C_HALF

               IFPDX = dcm_IFPD(1)
               IFPDY = dcm_IFPD(2)
               IFPDZ = dcm_IFPD(3)

               !$$--- size of Block
               NT = blockdim%x*blockdim%y

               !$$--- size of grid
               NB = griddim%x*griddim%y
            !end if
            !call syncthreads()

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1

            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
            DO LL=1, NL
              !$$IC -- the id of the atom on the device
              !$$      NOTE: IC is not the id of the same atom in the whole box
              IC= (IB-1)*NT+IT+(LL-1)*NAL
              if(IC.LE.NPRT) then
                 !$$--- NOTE by HOU Qing: Oct 10,2017
                 !$$         The activation examination added. only the activative atoms
                 !$$         need to calculate force
                 FPX = C_ZERO
                 FPY = C_ZERO
                 FPZ = C_ZERO
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE ) then
                    !$$POS-- position of the atom
                    !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                    POSX = XP(IC+IA0,1)
                    POSY = XP(IC+IA0,2)
                    POSZ = XP(IC+IA0,3)
                    ITYPI  = ITYP(IC+IA0)
                    DENKI = DEN(IC+IA0)

                    !$$-- start calculation of electron density
                    IIW = KVOIS(IC)
                    do IW=1, IIW
                       !$$--- NOTE: the particles index of neighbore-list is
                       !$$          the index of particle in the whole box
                       J=INDI(IC,IW)

                       !$$--- To calculate the seperation between particle IC
                       !$$    and its IWth neighbore
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

                       R2 = SEPX*SEPX+SEPY*SEPY+SEPZ*SEPZ
                       !$$--- To calculate electron density on atom I
                       ITYPJ=ITYP(J)
                       !$$if(R2 .LE. dcm_RU2(ITYPI,ITYPJ)) then
                       if(R2 .LE. RU2) then
                          KTAB  = dcm_KPAIR(ITYPI,ITYPJ)
                          KTAB1 = dcm_KPAIR(ITYPJ,ITYPI)
                          R     = DSQRT(R2)
                          SK    = DSQRT(R)*CSI
                          KK    = SK
                          DK    = SK-KK

                          DENKJ = DEN(J)
                          !$$--- interpolation of force table
                          !$$--- for  pairwise section
                          !EXP5 = FPOTR(KTAB,KK) + DK*(FPOTR(KTAB,KK+1) - FPOTR(KTAB,KK))

                          !$$--- for many-body section
                          !$$    Note; the electron density of atom I on atom J could not be the same as
                          !$$    from atom J on atom J, if atom I and J are not the same type.
                          !$$
                          !$$--- from the type of atom J to the type of atom I
                          !EXP6 = FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK+1) - FPOTB(KTAB,KK))

                          !$$--- from the type of atom I to the type of atom J
                          !EXP7 = FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK+1) - FPOTB(KTAB1,KK))

                          !FORTOT= (EXP5 + EXP7*DENKJ + EXP6*DENKI)

                          !$$--- NOTE by Hou Qing: Dec 4,2014
                          !$$          The force table FPOTR, FPOTB are dV(r)/dr/r, and dRHO/dr/r, in OLD VERSION
                          !$$          In the new version. we define FPOTR, FPOTB as (dV(r)/dr)*r and (dRHO/dr)
                          !$$--- NOTE by Hou Qing: Nov,7, 2014
                          !$$         A bug exist for cases when the e-density contribution is not symmetry for I->J and J->I.
                          !$$         This  is corrected by switching the position of  DENKI and DENKJ

                          FORTOT= (FPOTR(KTAB,KK)  + DK*(FPOTR(KTAB,KK+1)  - FPOTR(KTAB,KK)))/R2        + &
                                  ( (FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK+1)  - FPOTB(KTAB,KK))) *DENKI  + &
                                    (FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK+1) - FPOTB(KTAB1,KK)))*DENKJ  )/R


                          FPX = FPX + FORTOT*SEPX
                          FPY = FPY + FORTOT*SEPY
                          FPZ = FPZ + FORTOT*SEPZ
                       end if
                    end do
                 end if
                 FP(IC,1) = FPX
                 FP(IC,2) = FPY
                 FP(IC,3) = FPZ
              end if
            END DO !$$end loop for LL

        RETURN
  end subroutine CALFORCE_KERNEL
  !*************************************************************************************

  !*********************************************************************************
  subroutine CALFORCE_template(IDEV, STARTCELL, ENDCELL, KVOIS,INDI, ITYP, XP, STATU, FPOTR, FPOTB, DEN, FP)
  !***  PURPOSE:   to begine calculate the force.
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL, the ID of the first on the device
  !                ENDCELL,   the ID of last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                STATU:  the statu of atoms
  !                FPOTR,  the look-up table for pairwise force
  !                FPOTB,  the look-up table for tight-binding force
  !                DEN,    the Finnis-Sinclar term

  !     OUTPUT:    FP,      the forces on atoms
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      integer,      device, dimension(:)  ::KVOIS
      integer,      device, dimension(:,:)::INDI
      integer,      device, dimension(:)  ::ITYP
      real(KINDDF), device, dimension(:,:)::XP
      integer,      device, dimension(:)  ::STATU
      real(KINDDF), device, dimension(:,:)::FPOTR
      real(KINDDF), device, dimension(:,:)::FPOTB
      real(KINDDF), device, dimension(:)  ::DEN
      !--- the output
      real(KINDDF), device, dimension(:,:)::FP

      !--- Local variables
         integer::ERR, STARTA, ENDA, NPRT
      !$$--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

             ERR = cudaSetDevice(IDEV)

             STARTA  = hm_IA1th(STARTCELL)
             ENDA    = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
             NPRT    = ENDA - STARTA + 1

             blocks  = dim3(mp_BLOCKDIMX, 1, 1)
             threads = dim3(mp_BLOCKSIZE, 1, 1)
             STARTA  = STARTA-1
             call CALFORCE_KERNEL<<<blocks, threads>>>(dm_NPRT,ITYP,XP, m_NAPDEV, STATU, STARTA, NPRT, KVOIS,INDI, &
                                                       m_CSI, m_MXRU2, m_TABDIM(1), FPOTR, FPOTB, DEN, FP,         &
                                                       m_BOXSIZE(1), m_BOXSIZE(2), m_BOXSIZE(3))

          return
  end subroutine CALFORCE_template
  !*********************************************************************************

  !*********************************************************************************
  subroutine CALDEN_EAM_Force_Table2A_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to calculate the differential of emebedment function on atoms:
  !
  !                -dF(RHOO)/d(RHO)
  !
  !                Could be used for other applications, for exmaples
  !                in local stress caculations
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !     OUTPUT     dm_DEM,    the simulation box with force updated
  !
  !     SEE ALSO:  CALFORCE_EAM_Force_Table2A_DEV
  !                CALFTENSOR_EAM_Force_Table2A_DEV
  !                Cal_NeighboreList_DEVKERNEL2C
  !
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam

      !--- Local variables
         integer::ERR, CURDEV, I

             ERR = cudaGetDevice(CURDEV)

        !$$--- If pressure effect acounted for, the box size could be changed, we
        !$$    need to recopy the boxsize
             m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
             m_BOXSIZE(1:3)      = SimBox%ZL(1:3)

        !$$--- caluclate the -dFi(RHO)/d(RHO) on atoms
             do I=1, m_NDEVICE
                 call PRECALFOR_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                       &
                                         dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,                            &
                                         dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data,  dm_WorkSpace%STATU(I)%Data, &
                                         dm_FtableWS(I)%POTB,   dm_FtableWS(I)%DFEMBD, dm_FtableWS(I)%DEN)
             end do
                        
            !$$-- Because the subscript of DEN on host and on device are the index in the whole box,
            !$$   we need wait for all the DENs are ready, and copy the m_DEN on host into device
            !call SynchronizeDevices()
            !do I=1, m_NDEVICE
            !   ERR = cudaSetDevice(m_DEVICES(I))
            !   ERR = cudaMemcpyAsync(dm_FtableWS(I)%DEN, hm_DEN,dm_NPRT)
            !end do
             call Synchroniz_DEN_on_Devices()   

           ERR = cudaSetDevice(CURDEV)

          return
  end subroutine CALDEN_EAM_Force_Table2A_DEV
  !*********************************************************************************


  !*********************************************************************************
  subroutine CALFORCE_EAM_Force_Table2A_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to calculate the forces. For speedup, no virial tensor
  !                are calculted (ref;CALPTENSOR_EAM_Force_Table2A_DEV)
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
  !     NOTE:      The neigbore list created by Cal_NeighboreList_DEVKERNEL2C is used
  !                The atom index in NL cretated by Cal_NeighboreList_DEVKERNEL2C is
  !                the index in the whole simulation box.
  !
  !     SEE ALSO:  CALFTENSOR_EAM_Force_Table2A_DEV
  !                Cal_NeighboreList_DEVKERNEL2C
  !
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam

      !--- Local variables
         integer::ERR, CURDEV, I

             ERR = cudaGetDevice(CURDEV)

        !$$--- If pressure effect acounted for, the box size could be changed, we
        !$$    need to recopy the boxsize
             m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
             m_BOXSIZE(1:3)      = SimBox%ZL(1:3)

        !$$--- caluclate the -dFi(RHO)/d(RHO) on atoms
             do I=1, m_NDEVICE
                 call PRECALFOR_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                       &
                                         dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,                            &
                                         dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data,  dm_WorkSpace%STATU(I)%Data, &
                                         dm_FTableWS(I)%POTB,   dm_FTableWS(I)%DFEMBD, dm_FTableWS(I)%DEN)
             end do

            !$$-- Because the subscript of DEN are the indexed in the whole box, just like XP
            !$$   we need wait for all the DENs are ready.
            !$$   The DEN on host, hm_DEN, will be copy into devices in CALFORCE_template
            !  call SynchronizeDevices()
              call Synchroniz_DEN_on_Devices()  

        !$$--- caluclate the forces on atoms
             do I=1, m_NDEVICE
                 call CALFORCE_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                       &
                                        dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,                            &
                                        dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data, dm_WorkSpace%STATU(I)%Data,  &
                                        dm_FTableWS(I)%FPOTR,  dm_FTableWS(I)%FPOTB, dm_FTableWS(I)%DEN,                  &
                                        dm_WorkSpace%FP(I)%Data )
             end do

             ERR = cudaSetDevice(CURDEV)

          return
  end subroutine CALFORCE_EAM_Force_Table2A_DEV
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine CALPTENSOR_KERNEL(IM, ITYP, XP, NAPDEV, STATU, IA0, NPRT,KVOIS,INDI, CSI, RU2, KT, FPOTR, FPOTB, DEN,FP, VTENSOR, &
                                                   BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033)
  !***  PURPOSE:   to begine calculate the forces on atoms, and also the virial tensor
  !                Compared to CALFORCE_KERNEL, the box could be non-cubic.
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              ITYP :  the type of atoms corresponding to INDI
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              STATU:  the statu of atoms
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              NPRT:   the actual number of atoms on the device
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              CSI,    the interval for interpolation
  !$$              RU2,    the square of cutoff distance
  !$$              KT,     the first dimension of lookup table POTR and POTB
  !$$              FPOTR,  the look-up table for pairwise interaction
  !$$              FPOTB,  the look-up table for tight-binding interaction
  !$$              DEN,    the Finnis-Sinclar term
  !$$
  !$$   OUTPUT:   FP,      the force on on atoms
  !$$             VTENSOR,  the virial tensor obatined in one block
  !

      implicit none
      !--- DUMMY VARIABLES
      integer,      value,  intent(in)::IM
      integer,      device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer,      value,  intent(in)::NAPDEV, IA0,NPRT
      integer,      device, intent(in)::STATU(NAPDEV)
      integer,      device, intent(in)::KVOIS(NAPDEV)
      integer,      device, intent(in)::INDI(NAPDEV,*)

      real(KINDDF), value,  intent(in)::CSI, RU2
      integer,      value,  intent(in)::KT
      real(KINDDF), device, intent(in)::FPOTR(KT,*)
      real(KINDDF), device, intent(in)::FPOTB(KT,*)
      real(KINDDF), device, intent(in)::DEN(IM)
      real(KINDDF), value::BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033

      !Output:
      real(KINDDF), device, intent(out)::FP(NAPDEV,*)
      real(KINDDF), device, intent(out)::VTENSOR(3,3,*)

      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, KTAB1, ITYPI, ITYPJ
      real(KINDDF)::FORTOT,FPX,FPY,FPZ, SK, DK, R2, R ,EXP5, EXP6, EXP7, DENKI, DENKJ
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z
      real(KINDDF)::P11, P12, P13, P21, P22, P23, P31, P32, P33


      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ, BS11, BS12, BS13, BS21, BS22, BS23, BS31, BS32, BS33
      integer::IFPDX, IFPDY, IFPDZ
      real(KINDDF), shared::SP11(mp_MXSHARESIZE/9),SP12(mp_MXSHARESIZE/9),SP13(mp_MXSHARESIZE/9)
      real(KINDDF), shared::SP21(mp_MXSHARESIZE/9),SP22(mp_MXSHARESIZE/9),SP23(mp_MXSHARESIZE/9)
      real(KINDDF), shared::SP31(mp_MXSHARESIZE/9),SP32(mp_MXSHARESIZE/9),SP33(mp_MXSHARESIZE/9)
      ! Local variables
      integer::IC, IT, IB, LL, NT, NB, NL, NAL

            IT  = (threadidx%y-1)*blockdim%x + threadidx%x
            
               BX = BX0 !dcm_BOXSIZE(1)
               BY = BY0 !dcm_BOXSIZE(2)
               BZ = BZ0 !dcm_BOXSIZE(3)

               HBX = BX*C_HALF
               HBY = BY*C_HALF
               HBZ = BZ*C_HALF

               BS11 =  BS011 !dcm_BOXSHAPE(1,1)
               BS12 =  BS012 !dcm_BOXSHAPE(1,2)
               BS13 =  BS013 !dcm_BOXSHAPE(1,3)
               BS21 =  BS021 !dcm_BOXSHAPE(2,1)
               BS22 =  BS022 !dcm_BOXSHAPE(2,2)
               BS23 =  BS023 !dcm_BOXSHAPE(2,3)
               BS31 =  BS031 !dcm_BOXSHAPE(3,1)
               BS32 =  BS032 !dcm_BOXSHAPE(3,2)
               BS33 =  BS033 !dcm_BOXSHAPE(3,3)

               IFPDX = dcm_IFPD(1)
               IFPDY = dcm_IFPD(2)
               IFPDZ = dcm_IFPD(3)

               !$$--- size of Block
               NT = blockdim%x*blockdim%y

               !$$--- size of grid
               NB = griddim%x*griddim%y

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1

            !$$--- initialize the pressure
            P11= C_ZERO
            P12= C_ZERO
            P13= C_ZERO
            P21= C_ZERO
            P22= C_ZERO
            P23= C_ZERO
            P31= C_ZERO
            P32= C_ZERO
            P33= C_ZERO

            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
            DO LL=1, NL
              !$$IC -- the id of the atom on the device
              !$$      NOTE: IC is not the id of the same atom in the whole box
              IC= (IB-1)*NT+IT + (LL-1)*NAL
              if(IC.LE.NPRT) then
                 !$$--- NOTE by HOU Qing: Oct 10,2017
                 !$$         The activation examination added. only the activative atoms
                 !$$         need to calculate force
                 FPX = C_ZERO
                 FPY = C_ZERO
                 FPZ = C_ZERO
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE ) then
                    !$$POS-- position of the atom
                    !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                    POSX  = XP(IC+IA0,1)
                    POSY  = XP(IC+IA0,2)
                    POSZ  = XP(IC+IA0,3)
                    ITYPI = ITYP(IC+IA0)
                    DENKI = DEN(IC+IA0)

                    !$$-- start calculation of electron density
                    IIW = KVOIS(IC)
                    do IW=1, IIW
                       !$$--- NOTE: the particles index of neighbore-list is
                       !$$          the index of particle in the whole box
                       J=INDI(IC,IW)

                       !$$--- To calculate the seperation between particle IC
                       !$$    and its IWth neighbore
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

                       !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                       !$$          only when the shape of the box is cubic,
                       !$$          DXYZ = SEP
                       DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                       DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                       DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                       R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                       !$$--- To calculate electron density on atom I
                       ITYPJ=ITYP(J)
                       !$$if(R2 .LE. dcm_RU2(ITYPI,ITYPJ)) then
                       if(R2 .LE. RU2) then
                          KTAB  = dcm_KPAIR(ITYPI,ITYPJ)
                          KTAB1 = dcm_KPAIR(ITYPJ,ITYPI)
                          R     = DSQRT(R2)
                          SK    = DSQRT(R)*CSI
                          KK    = SK
                          DK    = SK-KK

                          DENKJ = DEN(J)
                          !$$--- interpolation of force table
                          !$$--- for  pairwise section
                          !EXP5 = FPOTR(KTAB,KK) + DK*(FPOTR(KTAB,KK+1) - FPOTR(KTAB,KK))

                          !$$--- for many-body section
                          !$$    Note; the electron density of atom I on atom J could not be the same as
                          !$$    from atom J on atom J, if atom I and J are not the same type.
                          !$$
                          !$$--- from the type of atom I to the type of atom J
                          !EXP6 = FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK+1) - FPOTB(KTAB,KK))

                          !$$--- from the type of atom J to the type of atom I
                          !EXP7 = FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK+1) - FPOTB(KTAB1,KK))

                          !FORTOT= EXP5 + EXP7*DENKJ + EXP6*DENKI

                          !$$--- NOTE by Hou Qing: Dec 4,2014
                          !$$          The force table FPOTR, FPOTB are dV(r)/dr/r, and dRHO/dr/r, in OLD VERSION
                          !$$          In the new version. we define FPOTR, FPOTB as (dV(r)/dr)*r and (dRHO/dr)
                          !$$--- NOTE by Hou Qing: Nov,7, 2014
                          !$$         A bug exist for cases when the e-density contribution is not symmetry for I->J and J->I.
                          !$$         This  is corrected by switching the position of  DENKI and DENKJ
                          FORTOT= (FPOTR(KTAB,KK)  + DK*(FPOTR(KTAB,KK+1)  - FPOTR(KTAB,KK))) /R2      + &
                                  ( (FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK+1)  - FPOTB(KTAB,KK))) *DENKI + &
                                    (FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK+1) - FPOTB(KTAB1,KK)))*DENKJ  )/R

                          FPX = FPX + FORTOT*SEPX
                          FPY = FPY + FORTOT*SEPY
                          FPZ = FPZ + FORTOT*SEPZ

                          !$$-- at the same time, we can calculate the pressure tensor.
                          !$$   NOTE: here we have consider periodic condition, thus
                          !$$         we cannot laterly use PRESS=SUM(XP*FP) to get the pressure.
                          !$$         It would be convenience and efficient to have the pressure
                          !$$         calculated here.
                          FORTOT = FORTOT*C_HALF

                          P11 = P11+DXYZ_X*DXYZ_X*FORTOT
                          P12 = P12+DXYZ_X*DXYZ_Y*FORTOT
                          P13 = P13+DXYZ_X*DXYZ_Z*FORTOT
                          P21 = P21+DXYZ_Y*DXYZ_X*FORTOT
                          P22 = P22+DXYZ_Y*DXYZ_Y*FORTOT
                          P23 = P23+DXYZ_Y*DXYZ_Z*FORTOT
                          P31 = P31+DXYZ_Z*DXYZ_X*FORTOT
                          P32 = P32+DXYZ_Z*DXYZ_Y*FORTOT
                          P33 = P33+DXYZ_Z*DXYZ_Z*FORTOT
                       end if
                    end do
                 end if
                 FP(IC,1) = FPX
                 FP(IC,2) = FPY
                 FP(IC,3) = FPZ
              end if
            END DO !$$end loop for LL

            !$$--- COPY VIRIAL TO SHARED MEMOMRY
            SP11(IT) = P11
            SP12(IT) = P12
            SP13(IT) = P13
            SP21(IT) = P21
            SP22(IT) = P22
            SP23(IT) = P23
            SP31(IT) = P31
            SP32(IT) = P32
            SP33(IT) = P33
            call syncthreads()

           !$$--- Summming the tensor
            if(IT.eq. 1) then
              !$$VTENSOR(1:3,1:3,IB) = C_ZERO
              P11 =  C_ZERO
              P12 =  C_ZERO
              P13 =  C_ZERO
              P21 =  C_ZERO
              P22 =  C_ZERO
              P23 =  C_ZERO
              P31 =  C_ZERO
              P32 =  C_ZERO
              P33 =  C_ZERO
              DO IW=1, blockdim%x*blockdim%y
                 P11 = P11 + SP11(IW)
                 P12 = P12 + SP12(IW)
                 P13 = P13 + SP13(IW)
                 P21 = P21 + SP21(IW)
                 P22 = P22 + SP22(IW)
                 P23 = P23 + SP23(IW)
                 P31 = P31 + SP31(IW)
                 P32 = P32 + SP32(IW)
                 P33 = P33 + SP33(IW)
              END DO
              VTENSOR(1,1,IB) =  P11
              VTENSOR(1,2,IB) =  P12
              VTENSOR(1,3,IB) =  P13
              VTENSOR(2,1,IB) =  P21
              VTENSOR(2,2,IB) =  P22
              VTENSOR(2,3,IB) =  P23
              VTENSOR(3,1,IB) =  P31
              VTENSOR(3,2,IB) =  P32
              VTENSOR(3,3,IB) =  P33
            end if

        return
  end subroutine CALPTENSOR_KERNEL
  !*************************************************************************************


  !*********************************************************************************
  subroutine CALPTENSOR_template(IDEV, STARTCELL, ENDCELL, KVOIS,INDI, ITYP, XP, STATU, FPOTR, FPOTB, DEN, FP, VTENSOR)

  !***  PURPOSE:   to calculate the force. At the same time, the virial tensor
  !                is also calculated and the box shape could be uncubic.
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL, the ID of the first on the device
  !                ENDCELL,   the ID of last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                STATU:  the statu of atoms
  !                FPOTR,  the look-up table for pairwise force
  !                FPOTB,  the look-up table for tight-binding force
  !                DEN,    the embeded function

  !     OUTPUT:    FORCE,  the forces on atoms
  !                VTENSOR, the virial tensor
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      integer,      device, dimension(:)::KVOIS
      integer,      device, dimension(:,:)::INDI
      integer,      device, dimension(:)::ITYP
      real(KINDDF), device, dimension(:,:)::XP
      integer,      device, dimension(:)  ::STATU
      real(KINDDF), device, dimension(:,:)::FPOTR
      real(KINDDF), device, dimension(:,:)::FPOTB
      real(KINDDF), device, dimension(:)  ::DEN
      !--- the output
      real(KINDDF), device, dimension(:,:)  ::FP
      real(KINDDF), device, dimension(:,:,:)::VTENSOR

      !--- Local variables
         integer::BX, BY, NBX, NBY, ERR, STARTA, ENDA, NPRT, I
      !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

            ERR = cudaSetDevice(IDEV)

            STARTA = hm_IA1th(STARTCELL)
            ENDA   = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT   = ENDA - STARTA + 1

            !$$--- to determine size of a block (the number of threads in a block)
            !$$    Note the difference between this and that calculating force
            BX = mp_MXSHARESIZE/9
            BY  = 1
            NBX = mp_BLOCKDIMX
            NBY = 1
            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call CALPTENSOR_KERNEL<<<blocks, threads>>>(dm_NPRT,ITYP,XP, m_NAPDEV, STATU, STARTA, NPRT, KVOIS,INDI,    &
                                                        m_CSI, m_MXRU2, m_TABDIM(1), FPOTR, FPOTB, DEN, FP, VTENSOR,   &
                                                        m_BOXSIZE(1),    m_BOXSIZE(2),    m_BOXSIZE(3),                &
                                                        m_BOXSHAPE(1,1), m_BOXSHAPE(1,2), m_BOXSHAPE(1,3),             &
                                                        m_BOXSHAPE(2,1), m_BOXSHAPE(2,2), m_BOXSHAPE(2,3),             &
                                                        m_BOXSHAPE(3,1), m_BOXSHAPE(3,2), m_BOXSHAPE(3,3)    )

          return
  end subroutine CALPTENSOR_template
  !*********************************************************************************


  !*********************************************************************************
  subroutine CALPTENSOR_EAM_Force_Table2A_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to calculate the forces.Virial tensor are also calculated
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
  !     NOTE:      The neigbore list created by Cal_NeighboreList_DEVKERNEL2C is used
  !                The atom index in NL cretated by Cal_NeighboreList_DEVKERNEL2C is
  !                global index in the whole simulation box.
  !
  !     SEE ALSO:  CALFORCE_EAM_Force_Table2A_DEV
  !                Cal_NeighboreList_DEVKERNEL2C
  !
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam

      !--- Local variables
         integer::ERR, CURDEV, I

             ERR = cudaGetDevice(CURDEV)
        !$$--- If pressure effect acounted for, the box size could be changed, we
        !$$    need to recopy the boxsize
             m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
             m_BOXSIZE(1:3)      = SimBox%ZL(1:3)

        !$$--- caluclate the -dFi(RHO)/d(RHO) on atoms
             do I=1, m_NDEVICE
                 call PRECALFOR_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                        &
                                         dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,                             &
                                         dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data,   dm_WorkSpace%STATU(I)%Data, &
                                         dm_FTableWS(I)%POTB,        dm_FTableWS(I)%DFEMBD,     dm_FTableWS(I)%DEN)
             end do


            !$$-- Because the subscript of DEN are the indexed in the whole box, just like XP
            !$$   we need wait for all the DENs are ready.
            !$$   The DEN on host, hm_DEN, will be copy into devices in CALPTENSOR_template
             ! call SynchronizeDevices()
             call Synchroniz_DEN_on_Devices()  

            !$$--- caluclate the forces on atoms, and viriral press
             do I=1, m_NDEVICE
                 call CALPTENSOR_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                        &
                                          dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,                             &
                                          dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data,  dm_WorkSpace%STATU(I)%Data,  &
                                          dm_FTableWS(I)%FPOTR,       dm_FTableWS(I)%FPOTB,     dm_FTableWS(I)%DEN,          &
                                          dm_WorkSpace%FP(I)%Data,    dm_FTableWS(I)%VTENSOR )
             end do


             !$$--- following summing the tensor obtained by all devices
             call SynchronizeDevices()
             call COPYOUT_VIRIALTENSOR(SimBox%VTENSOR)
             ERR = cudaSetDevice(CURDEV)

          return
  end subroutine CALPTENSOR_EAM_Force_Table2A_DEV
  !*********************************************************************************

  !*********************************************************************************
  subroutine COPYOUT_VIRIALTENSOR(hVTENSOR)
  !***  PURPOSE:   to copy the virial tensor on a device to host array
  !
  !
  !     INPUT:
  !     OUTPUT     hVTENSOR
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      real(KINDDF), dimension(:,:)::hVTENSOR


      !--- Local variables
         integer     ::ERR, CURDEV, I
         real(KINDDF)::VTENSOR(3,3,mp_BLOCKDIMX)


             ERR = cudaGetDevice(CURDEV)
             hVTENSOR = C_ZERO
             do I=1, m_NDEVICE
                ERR = cudaSetDevice(m_DEVICES(I))
                VTENSOR = dm_FTableWS(I)%VTENSOR
                hVTENSOR = hVTENSOR+sum(VTENSOR,dim=3)
             end do
             ERR = cudaSetDevice(CURDEV)

             !-- NOTE: the virial tensor have been devided by the number of boxs
             hVTENSOR = hVTENSOR/dble(dm_NPRT/dm_NAPB)

          return
  end subroutine COPYOUT_VIRIALTENSOR
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine CALEPOT_KERNEL(IM, ITYP,XP, NAPDEV, STATU, IA0, NPRT, KVOIS,INDI, CSI, RU2, KT, POTR, POTB, RHOD, KT1, FEMBD, EPOT,&
                                                BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033)
  !***  PURPOSE:   to calculate cohesive ennergy
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              ITYP :  the type of atoms corresponding to INDI
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              STATU:  the statu of atoms
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              NPRT:   the actual number of atoms on the device
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              CSI,    the interval for interpolation
  !$$              RU2,    the square of cutoff distance
  !$$              KT,     the first dimension of lookup table POTR and POTB
  !$$              POTR,   the look-up table for pairwise potential
  !$$              POTB,   the look-up table for tight-binding potential
  !
  !$$   OUTPUT:   EPOT,   the potential of atoms
  !$$   SEE ALSO: PRECALFOR_KERNEL
  !

      implicit none
      !--- DUMMY VARIABLES
      integer,      value,  intent(in)::IM
      integer,      device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer,      value,  intent(in)::NAPDEV, IA0, NPRT
      integer,      device, intent(in)::STATU(NAPDEV)
      integer,      device, intent(in)::KVOIS(NAPDEV)
      integer,      device, intent(in)::INDI(NAPDEV,*)

      real(KINDDF), value,  intent(in)::CSI, RU2
      integer,      value,  intent(in)::KT
      real(KINDDF), device, intent(in)::POTR(KT,*)
      real(KINDDF), device, intent(in)::POTB(KT,*)

      real(KINDDF), value,  intent(in)::RHOD
      integer,      value,  intent(in)::KT1
      real(KINDDF), device, intent(in)::FEMBD(KT1,*)
      real(KINDDF), value::BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033
      !--- output
      real(KINDDF), device::EPOT(*)

      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, ITYPI, ITYPJ
      real(KINDDF)::ER0, DEN0, SK,DK,R2,R
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z

      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ, BS11, BS12, BS13, BS21, BS22, BS23, BS31, BS32, BS33
      integer::IFPDX, IFPDY, IFPDZ
      integer::NT, NB

      ! Local variables
      integer::IC, IT, IB !, LL, NAL,NL

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              BX = BX0 !dcm_BOXSIZE(1)
              BY = BY0 !dcm_BOXSIZE(2)
              BZ = BZ0 !dcm_BOXSIZE(3)

              HBX = BX*C_HALF
              HBY = BY*C_HALF
              HBZ = BZ*C_HALF

              BS11 =  BS011 !dcm_BOXSHAPE(1,1)
              BS12 =  BS012 !dcm_BOXSHAPE(1,2)
              BS13 =  BS013 !dcm_BOXSHAPE(1,3)
              BS21 =  BS021 !dcm_BOXSHAPE(2,1)
              BS22 =  BS022 !dcm_BOXSHAPE(2,2)
              BS23 =  BS023 !dcm_BOXSHAPE(2,3)
              BS31 =  BS031 !dcm_BOXSHAPE(3,1)
              BS32 =  BS032 !dcm_BOXSHAPE(3,2)
              BS33 =  BS033 !dcm_BOXSHAPE(3,3)

              IFPDX = dcm_IFPD(1)
              IFPDY = dcm_IFPD(2)
              IFPDZ = dcm_IFPD(3)

              !$$--- size of Block
              NT = blockdim%x*blockdim%y

              !$$--- size of grid
              NB = griddim%x*griddim%y

              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x

              !$$IC -- the id of the atom on the device
              !$$      NOTE: IC is not the id of the same atom in the whole box
              IC= (IB-1)*NT+IT
              if(IC.LE.NPRT) then
                 !$$--- NOTE by HOU Qing: Oct 10,2017
                 !$$         The activation examination added. only the activative atoms
                 !$$         need to calculate force
                 ER0  = C_ZERO
                 DEN0 = C_ZERO
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE ) then

                    !$$POS-- position of the atom
                    !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                    POSX = XP(IC+IA0,1)
                    POSY = XP(IC+IA0,2)
                    POSZ = XP(IC+IA0,3)
                    ITYPI  = ITYP(IC+IA0)

                    !$$-- start calculation of electron density
                    IIW = KVOIS(IC)
                    do IW=1, IIW
                       !$$--- NOTE: the particles index of neighbore-list is
                       !$$          the index of particle in the whole box
                       J=INDI(IC,IW)

                       !$$--- To calculate the seperation between particle IC
                       !$$    and its IWth neighbore
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

                       !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                       !$$          only when the shape of the box is cubic,
                       !$$          DXYZ = SEP
                       DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                       DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                       DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                       R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                       !$$--- To calculate electron density on atom I
                       ITYPJ=ITYP(J)
                       !$$if(R2 .LE. dcm_RU2(ITYPI,ITYPJ)) then
                       if(R2 .LE. RU2) then
                          KTAB =  dcm_KPAIR(ITYPI,ITYPJ)
                          R    = DSQRT(R2)
                          SK   = DSQRT(R)*CSI
                          KK   = SK
                          DK   = SK-KK

                          !$$--- NOTE by HOU Qing: Dec 4,2014
                          !$$          The force table POTR, POTB are V(r), RHO(r)*0.5, in OLD VERSION
                          !$$          In the new version. we define POTR, POTB as V(r)*r and RHO(r)
                          ER0 = ER0  + (POTR(KTAB,KK)+DK*(POTR(KTAB,KK+1)-POTR(KTAB,KK)) )/R
                          DEN0= DEN0 + (POTB(KTAB,KK)+DK*(POTB(KTAB,KK+1)-POTB(KTAB,KK)) )
                       end if
                    end do

                    !$$**** NOTE:the subscript of EPOT is the index of atoms on the device
                    KTAB =  dcm_KEMBD(ITYPI)
                    SK   =  DEN0/RHOD + 1.D0
                    KK   =  SK + 0.000001D0
                    DEN0 =  FEMBD(KTAB, KK) + (SK-KK)*(FEMBD(KTAB, KK+1) -  FEMBD(KTAB, KK) )
                 end if
                 EPOT(IC) = ER0+DEN0
              end if
        RETURN
  end subroutine CALEPOT_KERNEL
  !*********************************************************************************

  !*********************************************************************************
  subroutine CALEPOT_template(IDEV, STARTCELL, ENDCELL, KVOIS, INDI, ITYP, XP, STATU, POTR, POTB, FEMBD, dEPOT)
  !***  PURPOSE:   to calculate potential of atom.
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL, the ID of the first on the device
  !                ENDCELL,   the ID of last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                STATU:  the statu ofd atoms
  !                POTR,   the look-up table for pairwise potential
  !                POTB,   the look-up table for tight-binding potential
  !                FEMBD,  the look-up table for embedment function

  !     OUTPUT:   dEPOT,     the potentials of atoms
  !
  !     SEE ALSO: PRECALFOR_template
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      integer,      device, dimension(:)  ::KVOIS
      integer,      device, dimension(:,:)::INDI
      integer,      device, dimension(:)  ::ITYP
      real(KINDDF), device, dimension(:,:)::XP
      integer,      device, dimension(:)  ::STATU
      real(KINDDF), device, dimension(:,:)::POTR
      real(KINDDF), device, dimension(:,:)::POTB
      real(KINDDF), device, dimension(:,:)::FEMBD
      real(KINDDF), device, dimension(:)  ::dEPOT

      !--- Local variables
         integer::BX, BY, NB, NBX, NBY, ERR, STARTA, ENDA, NPRT, I
         !$$--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer,parameter::BLOCKSIZE = 64, BLOCKDIMX=64


            ERR = cudaSetDevice(IDEV)
            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1

        !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            BX = min(BLOCKSIZE, mp_MXSHARESIZE)
            BY = 1

            NB  = (NPRT-1)/(BX*BY)+1
            NBX = min(NB,  BLOCKDIMX)
            NBY = (NB-1)/NBX+1
            NB  = NBX*NBY

            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call CALEPOT_KERNEL<<<blocks, threads>>>(dm_NPRT, ITYP, XP, m_NAPDEV, STATU, STARTA, NPRT, KVOIS,INDI,                   &
                                                      m_CSI, m_MXRU2, m_TABDIM(1), POTR, POTB, m_RHOD, m_TABEMDDIM(1), FEMBD, dEPOT, &
                                                      m_BOXSIZE(1),    m_BOXSIZE(2),    m_BOXSIZE(3),                                &
                                                      m_BOXSHAPE(1,1), m_BOXSHAPE(1,2), m_BOXSHAPE(1,3),                             &
                                                      m_BOXSHAPE(2,1), m_BOXSHAPE(2,2), m_BOXSHAPE(2,3),                             &
                                                      m_BOXSHAPE(3,1), m_BOXSHAPE(3,2), m_BOXSHAPE(3,3)   )
          return
  end subroutine CALEPOT_template
  !*********************************************************************************

  !*********************************************************************************
  subroutine UpdateEPOT_EAM_Force_Table2A_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to calculate the potential for each atom, 
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT:    dxm_EPOT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam
      !--- Local variables
       integer::ERR,  CURDEV, I


        !$$--- If pressure effect acounted for, the box size could be changed, we
        !$$    need to recopy the boxsize
             m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
             m_BOXSIZE(1:3)      = SimBox%ZL(1:3)

             ERR = cudaGetDevice(CURDEV)
        !$$--- caluclate the EPOT of atoms on devices
             do I=1, m_NDEVICE
                call CALEPOT_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                           &
                                         dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,                             &
                                         dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data,   dm_WorkSpace%STATU(I)%Data, &
                                         dm_FTableWS(I)%POTR,        dm_FTableWS(I)%POTB,       dm_FTableWS(I)%FEMBD,       &
                                         dm_WorkSpace%EPOT(I)%Data)
             end do

             ERR = cudaSetDevice(CURDEV)
          return
  end subroutine UpdateEPOT_EAM_Force_Table2A_DEV
  !*********************************************************************************

  !*********************************************************************************
  subroutine CALEPOT_EAM_Force_Table2A_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to calculate the potential for each atom, and copy the data into m_EPOT
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT:   m_EPOT, hm_EPOT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam
      !--- Local variables


            call UpdateEPOT_EAM_Force_Table2A_DEV(SimBox, CtrlParam)
            call CopyEPOTFrom_Devices_to_Host()

          return
  end subroutine CALEPOT_EAM_Force_Table2A_DEV
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine CAL_EAM_AtomicStress_KERNEL(IM, ITYP, XP, NAPDEV, STATU, IA0, NPRT, KVOIS, INDI, CSI, RU2, &
                                                             KT, FPOTR, FPOTB, DEN, AP,                                     &
                                                             BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033)
  !***  PURPOSE:   to begine calculate the the atomic virial tensor.
  !
  !   INPUT:     IM:     the number of particles concerned
  !              ITYP :  the type of atoms corresponding to INDI
  !              XP:     the positions of atoms
  !              NAPDEV: the max number of atoms on a device
  !              STATU:  the statu of atoms
  !              IA0:    the index (in the whole box)of the fisrt atom on the device
  !              NPRT:   the actual number of atoms on the device
  !              KVOIS:  the number of neighbors for atoms
  !              INDI:   the index for the neighbores
  !              CSI,    the interval for interpolation
  !              RU2,    the square of cutoff distance
  !              KT,     the first dimension of lookup table POTR and POTB
  !              FPOTR,  the look-up table for pairwise interaction
  !              FPOTB,  the look-up table for tight-binding interaction
  !              DEN,    the embeded function
  !
  !   OUTPUT:    AP,     the atomic pressure on on atoms
  !
  !

      implicit none
      !--- DUMMY VARIABLES
      integer,      value,  intent(in)::IM
      integer,      device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer,      value,  intent(in)::NAPDEV, IA0,NPRT
      integer,      device, intent(in)::STATU(NAPDEV)
      integer,      device, intent(in)::KVOIS(NAPDEV)
      integer,      device, intent(in)::INDI(NAPDEV,*)

      real(KINDDF), value,  intent(in)::CSI, RU2
      integer,      value,  intent(in)::KT
      real(KINDDF), device, intent(in)::FPOTR(KT,*)
      real(KINDDF), device, intent(in)::FPOTB(KT,*)
      real(KINDDF), device, intent(in)::DEN(IM)
      real(KINDDF), value,  intent(in)::BX0, BY0, BZ0, BS011, BS012, BS013, BS021, BS022, BS023, BS031, BS032, BS033
      !Output:
      real(KINDDF), device, intent(out)::AP(NAPDEV,*)

      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, KTAB1, ITYPI, ITYPJ
      real(KINDDF)::FORTOT, SK, DK, R2, R ,EXP5, EXP6, EXP7, DENKI, DENKJ
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z
      real(KINDDF)::P11, P12, P13, P21, P22, P23, P31, P32, P33

      real(KINDDF), shared::HBX, HBY, HBZ, BX, BY, BZ, BS11, BS12, BS13, BS21, BS22, BS23, BS31, BS32, BS33
      integer,      shared::IFPDX, IFPDY, IFPDZ
      integer,      shared::NT, NB

      ! Local variables
      integer::IC, IT, IB, LL, NL, NAL


            !if(IT .EQ. 1) then
               BX = BX0 !dcm_BOXSIZE(1)
               BY = BY0 !dcm_BOXSIZE(2)
               BZ = BZ0 !dcm_BOXSIZE(3)

               HBX = BX*C_HALF
               HBY = BY*C_HALF
               HBZ = BZ*C_HALF

               BS11 =  BS011 !dcm_BOXSHAPE(1,1)
               BS12 =  BS012 !dcm_BOXSHAPE(1,2)
               BS13 =  BS013 !dcm_BOXSHAPE(1,3)
               BS21 =  BS021 !dcm_BOXSHAPE(2,1)
               BS22 =  BS022 !dcm_BOXSHAPE(2,2)
               BS23 =  BS023 !dcm_BOXSHAPE(2,3)
               BS31 =  BS031 !dcm_BOXSHAPE(3,1)
               BS32 =  BS032 !dcm_BOXSHAPE(3,2)
               BS33 =  BS033 !dcm_BOXSHAPE(3,3)

               IFPDX = dcm_IFPD(1)
               IFPDY = dcm_IFPD(2)
               IFPDZ = dcm_IFPD(3)

               !$$--- size of Block
               NT = blockdim%x*blockdim%y

               !$$--- size of grid
               NB = griddim%x*griddim%y
            !end if
            !call syncthreads()

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             IB  =  (blockidx%y-1)*griddim%x +  blockidx%x
             DO LL=1, NL
              !$$IC -- the id of the atom on the device
              !$$      NOTE: IC is not the id of the same atom in the whole box
               IC= (IB-1)*NT+IT + (LL-1)*NAL
               if(IC.LE. NPRT) then
                 !$$--- NOTE by HOU Qing: Oct 10,2017
                 !$$         The activation examination added. only the activative atoms
                 !$$         need to calculate force
                 !$$--- initialize the pressure
                 P11= C_ZERO
                 P12= C_ZERO
                 P13= C_ZERO
                 P21= C_ZERO
                 P22= C_ZERO
                 P23= C_ZERO
                 P31= C_ZERO
                 P32= C_ZERO
                 P33= C_ZERO
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE ) then
                    !$$POS-- position of the atom
                    !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                    POSX = XP(IC+IA0,1)
                    POSY = XP(IC+IA0,2)
                    POSZ = XP(IC+IA0,3)
                    ITYPI  = ITYP(IC+IA0)
                    DENKI = DEN(IC+IA0)

                    !$$-- start calculation of electron density
                    IIW = KVOIS(IC)
                    do IW=1, IIW
                       !$$--- NOTE: the particles index of neighbore-list is
                       !$$          the index of particle in the whole box
                       J=INDI(IC,IW)

                       !$$--- To calculate the seperation between particle IC
                       !$$    and its IWth neighbore
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

                       !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                       !$$          only when the shape of the box is cubic,
                       !$$          DXYZ = SEP
                       DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                       DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                       DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                       R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                       !$$--- To calculate electron density on atom I
                       ITYPJ=ITYP(J)
                       !$$if(R2 .LE. dcm_RU2(ITYPI,ITYPJ)) then
                       if(R2 .LE. RU2) then
                          KTAB  = dcm_KPAIR(ITYPI,ITYPJ)
                          KTAB1 = dcm_KPAIR(ITYPJ,ITYPI)
                          R     = DSQRT(R2)
                          SK    = DSQRT(R)*CSI
                          KK    = SK
                          DK    = SK-KK

                          DENKJ = DEN(J)
                          !$$--- interpolation of force table
                          !$$--- for  pairwise section
                          !EXP5 = FPOTR(KTAB,KK) + DK*(FPOTR(KTAB,KK+1) - FPOTR(KTAB,KK))

                          !$$--- for many-body section
                          !$$    Note; the electron density of atom I on atom J could not be the same as
                          !$$    from atom J on atom J, if atom I and J are not the same type.
                          !$$
                          !$$--- from the type of atom I to the type of atom J
                          !EXP6 = FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK+1) - FPOTB(KTAB,KK))

                          !$$--- from the type of atom J to the type of atom I
                          !EXP7 = FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK+1) - FPOTB(KTAB1,KK))

                          !FORTOT= EXP5 + EXP7*DENKJI + EXP6*DENKI

                          !$$--- NOTE by Hou Qing: Dec 4,2014
                          !$$          The force table FPOTR, FPOTB are dV(r)/dr/r, and dRHO/dr/r, in OLD VERSION
                          !$$          In the new version. we define FPOTR, FPOTB as (dV(r)/dr)*r and (dRHO/dr)
                          !$$--- NOTE by Hou Qing: Nov,7, 2014
                          !$$         A bug exist for cases when the e-density contribution is not symmetry for I->J and J->I.
                          !$$         This  is corrected by switching the position of  DENKI and DENKJ
                          FORTOT= ( (FPOTR(KTAB,KK)  + DK*(FPOTR(KTAB,KK+1)  - FPOTR(KTAB,KK)))/R2      + &
                                  ( (FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK+1)  - FPOTB(KTAB,KK))) *DENKI  + &
                                    (FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK+1) - FPOTB(KTAB1,KK)))*DENKJ  )/R)

                          P11 = P11+DXYZ_X*DXYZ_X*FORTOT
                          P12 = P12+DXYZ_X*DXYZ_Y*FORTOT
                          P13 = P13+DXYZ_X*DXYZ_Z*FORTOT
                          P21 = P21+DXYZ_Y*DXYZ_X*FORTOT
                          P22 = P22+DXYZ_Y*DXYZ_Y*FORTOT
                          P23 = P23+DXYZ_Y*DXYZ_Z*FORTOT
                          P31 = P31+DXYZ_Z*DXYZ_X*FORTOT
                          P32 = P32+DXYZ_Z*DXYZ_Y*FORTOT
                          P33 = P33+DXYZ_Z*DXYZ_Z*FORTOT
                       end if
                    end do
                 end if
                 AP(IC,1) = P11
                 AP(IC,2) = P12
                 AP(IC,3) = P13
                 AP(IC,4) = P21
                 AP(IC,5) = P22
                 AP(IC,6) = P23
                 AP(IC,7) = P31
                 AP(IC,8) = P32
                 AP(IC,9) = P33
              end if
            END DO !$$end loop for LL

        RETURN
  end subroutine CAL_EAM_AtomicStress_KERNEL
  !*************************************************************************************

  !*********************************************************************************
  subroutine CAL_EAM_AtomicStress_template(IDEV, STARTCELL, ENDCELL, KVOIS, INDI, ITYP, XP, STATU, FPOTR, FPOTB, DEN, AVP)

  !***  PURPOSE:   to calculate the force. At the same time, the virial tensor
  !                is also calculated and the box shape could be uncubic.
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL, the ID of the first on the device
  !                ENDCELL,   the ID of last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                STATU:  the statu of atoms
  !                FPOTR,  the look-up table for pairwise force
  !                FPOTB,  the look-up table for tight-binding force
  !                DEN,    the Finnis-Sinlar term

  !     OUTPUT:    AVP,    the atomic pressure tensor
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      integer,      device, dimension(:)  ::KVOIS
      integer,      device, dimension(:,:)::INDI
      integer,      device, dimension(:)  ::ITYP
      real(KINDDF), device, dimension(:,:)::XP
      integer,      device, dimension(:)  ::STATU
      real(KINDDF), device, dimension(:,:)::FPOTR
      real(KINDDF), device, dimension(:,:)::FPOTB
      real(KINDDF), device, dimension(:)  ::DEN
      !--- the output
      real(KINDDF), device, dimension(:,:)::AVP

      !--- Local variables
         integer::BX, BY, NB, NBX, NBY, ERR, STARTA, ENDA, NPRT, I
      !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

             ERR = cudaSetDevice(IDEV)
            !$$--- The boxsize, boxshape have been copyin device
            !$$    in calculating DEN, we do not need to copy again
            !$$    BOXSHAPE = m_BOXSHAPE
            !$$    BOX      = m_BOXSIZE
            !$$    IFPD     = m_IFPD
            !$$    RU2      = m_RU2
            !$$    KPAIR    = m_KPAIR

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1

            !$$--- to determine size of a block (the number of threads in a block)
            !$$    Note the difference between this and that calculating force
            BX      = mp_BLOCKSIZE
            BY      = 1
            NBX     = mp_BLOCKDIMX
            NBY     = 1
            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call CAL_EAM_AtomicStress_KERNEL<<<blocks, threads>>>(dm_NPRT, ITYP, XP, m_NAPDEV, STATU, STARTA, NPRT, KVOIS,INDI, &
                                                        m_CSI, m_MXRU2, m_TABDIM(1), FPOTR, FPOTB, DEN, AVP,                    &
                                                        m_BOXSIZE(1),    m_BOXSIZE(2),    m_BOXSIZE(3),                         &
                                                        m_BOXSHAPE(1,1), m_BOXSHAPE(1,2), m_BOXSHAPE(1,3),                      &
                                                        m_BOXSHAPE(2,1), m_BOXSHAPE(2,2), m_BOXSHAPE(2,3),                      &
                                                        m_BOXSHAPE(3,1), m_BOXSHAPE(3,2), m_BOXSHAPE(3,3)   )
          return
  end subroutine CAL_EAM_AtomicStress_template
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_EAM_AtomicStressTensor_DEV(IDEV, dAVP)
  !***  PURPOSE:   to calculate the atomic virial tensor
  !
  !     INPUT:     IDEV:      the ID of the device
  !     OUTPUT     dAVP:      the atomic virital tensor on device IDEV
  !                hAVP:      optioanl, the virial tensor on host
  !
  !     NOTE:      It has been assumed that atomic density have been calculated
  !                by calling subroutine CALDEN_EAM_Force_Table2A_DEV
  !                or calling subroutine CALFORCE_EAM_Force_Table2A_DEV
  !     SEE ALSO:  CALDEN_EAM_Force_Table2A_DEV
  !                CALFORCE_EAM_Force_Table2A_DEV
  !
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU,  only: dm_Neighbors
   implicit none
      !--- dummy vaiables
      integer, intent(in)::IDEV
      real(KINDDF), device, dimension(:,:), intent(out)::dAVP
      !--- Local variables

         !$$--- Now, we start caluclate the atomic virial
         !$$    NOTE: it is assumed that the dxm_DEN are complete
         !$$          and available. Refer to force calculation
         !$$

          call CAL_EAM_AtomicStress_template(m_DEVICES(IDEV), m_STARTCELL(IDEV), m_ENDCELL(IDEV),                                        &
                                             dm_Neighbors%KVOIS(IDEV)%Data, dm_Neighbors%INDI(IDEV)%Data,                                &
                                             dm_WorkSpace%ITYP(IDEV)%Data,  dm_WorkSpace%XP(IDEV)%Data,   dm_WorkSpace%STATU(IDEV)%Data, &
                                             dm_FTableWS(IDEV)%FPOTR,  dm_FTableWS(IDEV)%FPOTB, dm_FTableWS(IDEV)%DEN,    &
                                             dAVP)
          return
  end subroutine Cal_EAM_AtomicStressTensor_DEV
  !*********************************************************************************

  !*********************************************************************************
  logical function IfInit_EAM_Force_Table() result(YES)
  implicit none
          YES = m_INITED
          return
  end function IfInit_EAM_Force_Table
  !*********************************************************************************


  end module MD_EAM_Force_Table_GPU





