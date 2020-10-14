  module MD_TYPEDEF_SimBox_DEV
  !***  DESCRIPTION: this module is to define the simulation box for GPU.
  !                  This version is adapted from MD_Globle_Variables_GPU.F90
  !                  for the convenice of structuring the darta.
  !                  ______________________________________________________
  !                  HOU Qing, April, 2012
  ! **** HOSTORY:
  !       * 2018-09-25 (HOU Qing):  First version
  !

  use CUDAFOR
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none

          !--- The list of devices to be used
          integer, parameter::m_MXDEVICE=6                                 !The number limitation of devices to be used.
          integer           ::m_NDEVICE=0                                  !The actual number of devices to be used
          integer           ::m_DEVICES(0:m_MXDEVICE)=(/0,0,1,2,3,4,5/)    !The index of devices to be used

           !--- The GPU random number generator
           !--- uniform generator
           integer(kind=int_ptr_kind())::dm_URanGen(m_MXDEVICE) = 0
           !--- gaussian generator
           integer(kind=int_ptr_kind())::dm_GRanGen(m_MXDEVICE) = 0

          !--- the data type used for storing temp. data ---------------
          type::DeviceSwapData
                integer                                          ::IDEV
                real(KINDDF), device, dimension(:),  allocatable::Swap2D
          end   type DeviceSwapData

          !--- the data type used for storing configuration data on a device
           type::SimBoxData_DEV
                 !$$---  The sorted configuration on device 1
                 integer                                          ::IDEV
                 integer                                          ::NPRT 
                 integer,      device, dimension(:),   allocatable::ITYP
                 real(KINDDF), device, dimension(:,:), allocatable::XP
                 real(KINDDF), device, dimension(:,:), allocatable::XP1
                 real(KINDDF), device, dimension(:,:), allocatable::XP2
                 real(KINDDF), device, dimension(:,:), allocatable::XP3
                 real(KINDDF), device, dimension(:,:), allocatable::XP4
                 real(KINDDF), device, dimension(:,:), allocatable::XP5
                 real(KINDDF), device, dimension(:,:), allocatable::DIS
                 real(KINDDF), device, dimension(:,:), allocatable::FP
                 real(KINDDF), device, dimension(:),   allocatable::EPOT
                 real(KINDDF), device, dimension(:),   allocatable::EKIN
                 integer,      device, dimension(:),   allocatable::STATU

                 integer,      device, dimension(:),   allocatable::NAC
                 integer,      device, dimension(:),   allocatable::IA1th
                 integer,      device,  dimension(:),  allocatable::GID
          end    type SimBoxData_DEV


          !--- the data encapsual the GPU configurations
          integer, parameter::mp_NOTPART = 0
          integer, parameter::mp_HASPART = 1
                                     
          type::MDDevConfig
                 !--- The list of devices to be used
                 integer::NDEVICE=0                                  !The actual number of devices to be used
                 integer::DEVICES(0:m_MXDEVICE)=(/0,0,1,2,3,4,5/)    !The index of devices to be used

                 !--- The host variables of a simulation box
                 integer           ::HasPart = mp_NOTPART            !The flag indicating if the originalconfiguration
                                                                     !has been partioned.
                                                                     !if hm_HasPart = mp_NOTPART, the neighbore-list
                                                                     !routine will used m_XP to do partition,
                                                                     !if hm_HasPart >= mp_HASPART, the neighbore-list
                                                                     !will copyout the cuurent configure on devices, and
                                                                     !then performing partiion
                                                                     !
                                                                     !Any time the orginal
                                                                     !simulations boxes has been changed externally
                                                                     !hm_HasPart will be set to mp_NOTPART.
                                                                     !See also MD_NeighborsList_GPU.F90
                                                                     !

                 integer::dm_NPRT     = 0                            !The total number of particles in whole simulation system
                 integer::dm_NPPB     = 0                            !The number of particles in one box
                 integer::dm_NGROUP   = 0                            !The number of types of particles
                 integer,      dimension(:),   allocatable::m_ITYP   !The type IDs of the particles
                 real(KINDDF), dimension(:,:), allocatable::m_XP     !The positions of the particles
                 real(KINDDF), dimension(:,:), allocatable::m_XP1    !The velocities of the particles
                 real(KINDDF), dimension(:,:), allocatable::m_XP2    !The second order of XP, if Gear scheme to be used
                 real(KINDDF), dimension(:,:), allocatable::m_XP3    !The third order  of XP, if Gear scheme to be used
                 real(KINDDF), dimension(:,:), allocatable::m_XP4    !The fourth order of XP, if Gear scheme to be used
                 real(KINDDF), dimension(:,:), allocatable::m_XP5    !The fifith order of XP, if Gear scheme to be used
                 real(KINDDF), dimension(:,:), allocatable::m_DIS    !The displacements of particles
                 integer,      dimension(:),   allocatable::m_STATU  !The status of particles

           !--- The device variables needed for initial partition the system
                 real(KINDDF), device, dimension(:,:), allocatable:: dm_XP
                 integer,      device, dimension(:),   allocatable:: dm_STATU

           !--- The cell IDS the devices will handle
           !    The cell IDS to be generated by calling Neighbor-list routines
                integer::m_STARTCELL(m_MXDEVICE)=(/0,1,2,3,4,5/)            !The cell ID of the first cell on a device
                integer::m_ENDCELL(m_MXDEVICE)  =(/0,1,2,3,4,5/)            !The cell ID of the last cell on a device
                integer::m_NAPDEV= 0                                        !The max number of particles for a device handling

                real(KINDDF)::m_RNAPDEV =1.2                                !The redundance of particle on a device
                                                                            ! m_NAPDEV = (number of particle/number of devices)*m_RNAPDEV

           !---  The working spaceing on host
           !
           !---  The copies of configurations of the box, however, the atoms
           !     are sorted by cells. Sorting are perforemed in the neighbore
           !     list calculation.
           !     SEE: MD_NeighborsList_GPU.F90
           !
           !     NOTE: We only need the array size of GID, ITYP and XP to
           !           cover all particles in the whole simulation box, that is
           !           size of GID, ITYP and XP is dm_NPRT;
           !
           !           For XP1-XP5, DIS, FP, EPOT  etc cover only the paticles
           !           on devices, the sizes of these arraies are thus m_NAPDEV
           !
           !           On initializing this module, hm_EPOT and dx_EPOT are also
           !           allocated. But the  values of hm_EPOT and  dx_EPOT are only
           !           available only after when UpdateEPOT_xxx, or CalEPOT_xxx
           !           is called (see  MD_FS_ForceTable_GPU.F90 adn MD_EAM_ForceTable_GPU.F90)
           !
           !            The number of particles in cells hm_NAC (dx_NAC) and
           !            the index of starting particles in each cell hm_IA1th (dx_IA1th) are
           !            allocated on calling Set_BoxCell_Number in Initialize_NeighboreList_DEV.
           !            Their values are assigned only after calling Cal_NeighBoreList_DEV.
           !
           !    SEE:   Allocate_Working_Variables in this module
                 integer, dimension(:), allocatable::hm_GID                 !Index transformation for the partitioned box to original box.
                                                                            !That is, the Ith particle in partitioned box is the hm_GIG(I)th
                                                                            !particle in the original box
                                                                            !hm_GID is created by neighbor list calculation.

                 integer, dimension(:), allocatable::hm_GIDINV              !The inverse index transformation:
                                                                            !The Ith particle in the orginal box is the hm_GIDINVth pariticle
                                                                            !in the partitioned box

                                                                            !grouped by cells. Created by neighbor list calculation.
                 integer,      pinned,dimension(:),    allocatable::hm_ITYP !The type ID of the particles that are grouped by cells
                 real(KINDDF), pinned, dimension(:,:), allocatable::hm_XP   !The positions  of the particles that are grouped by cells
                 real(KINDDF),         dimension(:,:), allocatable::hm_XP1  !The velocities of the particles that are grouped by cells
                 real(KINDDF),         dimension(:,:), allocatable::hm_XP2  !The second order of hm_XP
                 real(KINDDF),         dimension(:,:), allocatable::hm_XP3  !The third order of hm_XP
                 real(KINDDF),         dimension(:,:), allocatable::hm_XP4  !The fourth order of hm_XP
                 real(KINDDF),         dimension(:,:), allocatable::hm_XP5  !The fifith order of hm_XP
                 real(KINDDF),         dimension(:,:), allocatable::hm_DIS  !The displacements of paricles that are grouped by cells
                 real(KINDDF),         dimension(:,:), allocatable::hm_FP   !The forces on paricles that are grouped by cells
                 real(KINDDF),         dimension(:),   allocatable::hm_EPOT !The potential of paricles that are grouped by cells
                 real(KINDDF),         dimension(:),   allocatable::hm_EKIN !The kinetic energy of paricles that are grouped by cells
                 integer,      pinned, dimension(:),   allocatable::hm_STATU !The status of paricles that are grouped by cells

                 integer, dimension(:), allocatable::hm_NAC                  !The number of particles in a cell, created in neighborlist calculations
                 integer, dimension(:), allocatable::hm_IA1th                !the subscrition index of the first particles in a cell

                 type(SimBoxData_DEV),  dimension(:), allocatable::dm_SimBox !The configuration data on devices
      
          end type MDDevConfig   
          

  !--- interface of routines in this module -------------------
  !---------------------------------------------------------
          private:: Allocate_Working_Variables
  !---------------------------------------------------------
          private:: Allocate_Working_Variables_DevData
  !---------------------------------------------------------
          private:: Allocate_Cell_Variables_template

  !---------------------------------------------------------
          public::  Clear_DevConfig
  !---------------------------------------------------------
          public::  Clear_DeviceSwapData
  !---------------------------------------------------------
          private:: Clear_Working_Variables
  !---------------------------------------------------------
          private:: Clear_Working_Variables_DevData

  !---------------------------------------------------------
  !       to initialize MDDevConfig
          private:: Initialize_A_DevConfig, &
                    Initialize_B_DevConfig
          interface Initialize_DevConfig
             module procedure Initialize_A_DevConfig
             module procedure Initialize_B_DevConfig
          end interface


 contains

  !*********************************************************************
          subroutine Clear_DeviceSwapData(SwapData)
          implicit none
              type(DeviceSwapData)::SwapData
              !--- local varaiables
              integer::CURDEV, ERR

                 if( SwapData%IDEV .lt. 0) return

                 ERR = cudaGetDevice(CURDEV)
                 ERR = cudaSetDevice(SwapData%IDEV)
                 SwapData%IDEV  = -1 
                 if(allocated(SwapData%Swap2D))  deallocate( SwapData%Swap2D)
                 ERR = cudaSetDevice(CURDEV)
                 return
          end subroutine Clear_DeviceSwapData
  !*********************************************************************

  !*********************************************************************
          subroutine Clear_SimBoxData_DEV(SimBox)
          implicit none
              type(SimBoxData_DEV)::SimBox
              !--- local varaiables
              integer::CURDEV, ERR

                 if( SimBox%IDEV .lt. 0) return

                 ERR = cudaGetDevice(CURDEV)
                 ERR = cudaSetDevice(SimBox%IDEV)
                 SimBox%IDEV  = -1 
                 SimBox%NPRT  = 0
                 if(allocated(SimBox%ITYP ))   deallocate( SimBox%ITYP )
                 if(allocated(SimBox%XP   ))   deallocate( SimBox%XP   )
                 if(allocated(SimBox%XP1  ))   deallocate( SimBox%XP1  )
                 if(allocated(SimBox%XP2  ))   deallocate( SimBox%XP2  )
                 if(allocated(SimBox%XP3  ))   deallocate( SimBox%XP3  )
                 if(allocated(SimBox%XP4  ))   deallocate( SimBox%XP4  )
                 if(allocated(SimBox%XP5  ))   deallocate( SimBox%XP5  )
                 if(allocated(SimBox%DIS  ))   deallocate( SimBox%DIS  )
                 if(allocated(SimBox%FP   ))   deallocate( SimBox%FP   )
                 if(allocated(SimBox%EPOT ))   deallocate( SimBox%EPOT )
                 if(allocated(SimBox%EKIN ))   deallocate( SimBox%EKIN )
                 if(allocated(SimBox%STATU))   deallocate( SimBox%STATU)
                 if(allocated(SimBox%NAC))     deallocate( SimBox%NAC)
                 if(allocated(SimBox%IA1th))   deallocate( SimBox%IA1th)
                 ERR = cudaSetDevice(CURDEV)
                 return
          end subroutine Clear_SimBoxData_DEV
  !*********************************************************************

  !****************************************************************************
  subroutine Initialize_B_DevConfig( SimBox, CtrlParam, SimBoxDev)
  !***  PURPOSE:  to initialize this module by allocating global device memories
  !               and working memories to be used
  !     INPUT:    SimBox,    the Simulation box
  !               CtrlParam, the control parameters
  !
  !     OUTPUT:
      implicit none
      !--- dummy variables
      type(SimMDBox)    :: SimBox
      type(SimMDCtrl)   :: CtrlParam
      type(MDDevConfig) :: SimBoxDev
      !--- Local vairables
      integer::ERR, CURDEV, NPRT, I

         !---
              ERR = cudaGetDevice(CURDEV)
         !---
              call Clear_MDDevConfig(SimBoxDev)
         !$$--- copy the simulation box to device memory
              NPRT                = SimBox%NPRT
              SimBoxDev%dm_NPRT   = NPRT
              SimBoxDev%dm_NPPB   = NPRT
              SimBoxDev%dm_NGROUP = SimBox%NGROUP

              !$$--- allocated host memory
              allocate(SimBoxDev%m_ITYP(NPRT), SimBoxDev%m_XP(NPRT,3), SimBoxDev%m_XP1(NPRT,3),SimBoxDev%m_DIS(NPRT,3), &
                       SimBoxDev%m_STATU(NPRT), STAT=ERR)
              if(ERR) goto 100


              !$$--- allocate memory if Nordic scheme to be used
              if(allocated(SimBox%XP2)) then
                 allocate(SimBoxDev%m_XP2(NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox%XP3)) then
                 allocate(SimBoxDev%m_XP3(NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox%XP4)) then
                 allocate(SimBoxDev%m_XP4(NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox%XP5)) then
                 allocate(SimBoxDev%m_XP5(NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              !$$--- allocated device memory
              call Allocate_Working_Variables(SimBoxDev)

             !$$--- copy the simulation box to device memory
              SimBoxDev%m_ITYP (1:NPRT)      = SimBox%ITYP( 1:NPRT)
              SimBoxDev%m_XP   (1:NPRT,1:3)  = SimBox%XP   (1:NPRT,1:3)
              SimBoxDev%m_XP1  (1:NPRT,1:3)  = SimBox%XP1  (1:NPRT,1:3)
              SimBoxDev%m_DIS  (1:NPRT,1:3)  = SimBox%DIS  (1:NPRT,1:3)
              SimBoxDev%m_STATU(1:NPRT)      = SimBox%STATU(1:NPRT)

              !$$---  also copy in position and statu into device
              !$$NOTE: it seems that the segment copy between host  and
              !$$      device could not work correctly when the two arraies
              !$$      are two-dimensional and have different size.
              !$$      To be safe, we change the assigment statement to
              !$$      cudaMemcp
              ERR = cudaSetDevice(SimBoxDev%DEVICES(0))
              allocate(SimBoxDev%dm_XP(NPRT,3), SimBoxDev%dm_STATU(NPRT), STAT=ERR)
              ERR = cudaMemcpyAsync(SimBoxDev%dm_XP(1,1),  SimBoxDev%m_XP(1,1), 3*NPRT)
              ERR = cudaMemcpyAsync(SimBoxDev%dm_STATU(1), SimBoxDev%m_STATU(1),  NPRT)
              SimBoxDev%HasPart = mp_NOTPART

              ERR = cudaSetDevice(CURDEV)

              !$$--- allocate the simulation box data type on for devices
              SimBoxDev%NDEVICE = m_NDEVICE                                 !The actual number of devices to be us
              SimBoxDev%DEVICES = m_DEVICES
              allocate (SimBoxDev%dm_SimBox(SimBoxDev%NDEVICE))
              do I = 1, SimBoxDev%NDEVICE
                 SimBoxDev%dm_SimBox(I)%IDEV = SimBoxDev%DEVICES(I)
              end do

              return

    100       write(*,*) "MDPSCU Error in allocating global device memory"
              stop
  end subroutine Initialize_B_DevConfig
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_A_DevConfig( SimBox, CtrlParam, SimBoxDev)
  !***  PURPOSE:  to initialize this module by allocating global device memories
  !               and working memories to be used
  !     INPUT:    SimBox,    the Simulation box
  !               CtrlParam, the control parameters
  !
  !     OUTPUT:
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:):: SimBox
      type(SimMDCtrl)             :: CtrlParam
      type(MDDevConfig)           :: SimBoxDev
      !--- Local vairables
      integer::ERR, CURDEV, NB, NPRT, I, IS
         !---
              ERR = cudaGetDevice(CURDEV)
         !---
              call Clear_MDDevConfig(SimBoxDev)

              NPRT                = SimBox(1)%NPRT
              NB                  = size(SimBox)
              SimBoxDev%dm_NPRT   = NB*NPRT
              SimBoxDev%dm_NPPB   = NPRT
              SimBoxDev%dm_NGROUP = SimBox(1)%NGROUP
              !$$--- allocated host memory
              allocate(SimBoxDev%m_ITYP(NB*NPRT), SimBoxDev%m_XP(NB*NPRT,3),SimBoxDev%m_XP1(NB*NPRT,3),SimBoxDev%m_DIS(NB*NPRT,3),  &
                       SimBoxDev%m_STATU(NB*NPRT), STAT=ERR)
              if(ERR) goto 100


              !$$--- allocate memory if Nordic scheme to be used
              if(allocated(SimBox(1)%XP2)) then
                 allocate(SimBoxDev%m_XP2(NB*NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox(1)%XP3)) then
                 allocate(SimBoxDev%m_XP3(NB*NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox(1)%XP4)) then
                 allocate(SimBoxDev%m_XP4(NB*NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox(1)%XP5)) then
                 allocate(SimBoxDev%m_XP5(NB*NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              !$$--- allocated device memory
              call Allocate_Working_Variables(SimBoxDev)


             !$$--- copy the simulation box to device memory
              IS = 0
              do I=1, NB
                 SimBoxDev%m_ITYP(1+IS:NPRT+IS)    = SimBox(I)%ITYP(1:NPRT)
                 SimBoxDev%m_XP(1+IS:NPRT+IS,1:3)  = SimBox(I)%XP(1:NPRT,1:3)
                 SimBoxDev%m_XP1(1+IS:NPRT+IS,1:3) = SimBox(I)%XP1(1:NPRT,1:3)
                 SimBoxDev%m_DIS(1+IS:NPRT+IS,1:3) = SimBox(I)%DIS(1:NPRT,1:3)
                 SimBoxDev%m_STATU(1+IS:NPRT+IS)   = SimBox(I)%STATU(1:NPRT)
                 IS = IS+NPRT
              end do

              !$$---  also copy in position and statu into device
              !$$NOTE: it seems that the segment copy between host  and
              !$$      device could not work correctly when the two arraies
              !$$      are two-dimensional and have different size.
              !$$      To be safe, we change the assigment statement to
              !$$      cudaMemcp
              ERR = cudaSetDevice(m_DEVICES(0))
              allocate(SimBoxDev%dm_XP(NB*NPRT,3), SimBoxDev%dm_STATU(NB*NPRT), STAT=ERR)
              ERR = cudaMemcpyAsync(SimBoxDev%dm_XP(1,1), SimBoxDev%m_XP(1,1), 3*NB*NPRT)
              ERR = cudaMemcpyAsync(SimBoxDev%dm_STATU(1),SimBoxDev%m_STATU(1),  NB*NPRT)
              SimBoxDev%HasPart = mp_NOTPART
              ERR = cudaSetDevice(CURDEV)

              !$$--- allocate the simulation box data type on for devices
              SimBoxDev%NDEVICE = m_NDEVICE                                 !The actual number of devices to be us
              SimBoxDev%DEVICES = m_DEVICES
              allocate (SimBoxDev%dm_SimBox(SimBoxDev%NDEVICE))
              do I = 1, SimBoxDev%NDEVICE
                 SimBoxDev%dm_SimBox(I)%IDEV = SimBoxDev%DEVICES(I)
              end do

              return

    100       write(*,*) "MDPSCU Error in allocating global device memory"
              stop
  end subroutine Initialize_A_DevConfig
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_DevConfig( SimBoxDev )
  !***  PURPOSE:   to clear the allocated global device memories and working memories
  !     INPUT
  !
  !
      implicit none
      !--- dummy variables
      type(MDDevConfig)      :: SimBoxDev
      !--- Local vairables
      integer::ERR, CURDEV

       !---
              ERR = cudaGetDevice(CURDEV)

       !$$---   clear the device memory on the major device
              SimBoxDev%dm_NPRT  = 0
              SimBoxDev%dm_NPPB  = 0
              SimBoxDev%dm_nGROUP = 0

              if(allocated(SimBoxDev%dm_XP)) then
                 deallocate(SimBoxDev%dm_XP)
              end if

              if(allocated(SimBoxDev%dm_STATU)) then
                 deallocate(SimBoxDev%dm_STATU)
              end if

        !$$--- dealoocate memory on host
              if(allocated(SimBoxDev%m_ITYP)) then
                 deallocate(SimBoxDev%m_ITYP)
              end if

              if(allocated(SimBoxDev%m_XP)) then
                 deallocate(SimBoxDev%m_XP)
              end if

              if(allocated(SimBoxDev%m_XP1)) then
                 deallocate(SimBoxDev%m_XP1)
              end if

              if(allocated(SimBoxDev%m_XP2)) then
                 deallocate(SimBoxDev%m_XP2)
              end if

              if(allocated(SimBoxDev%m_XP3)) then
                 deallocate(SimBoxDev%m_XP3)
              end if

              if(allocated(SimBoxDev%m_XP4)) then
                 deallocate(SimBoxDev%m_XP4)
              end if

              if(allocated(SimBoxDev%m_XP5)) then
                 deallocate(SimBoxDev%m_XP5)
              end if

              if(allocated(SimBoxDev%m_DIS)) then
                 deallocate(SimBoxDev%m_DIS)
              end if

              if(allocated(SimBoxDev%m_STATU)) then
                 deallocate(SimBoxDev%m_STATU)
              end if

              call Clear_Working_Variables(SimBoxDev)

              return
  end subroutine Clear_DevConfig
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Variables_DevData(theData, NPRT, NAPDEV)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     theData:  the working data on a device
  !
  !     OUTPUT:    theData,  the working data with allocated ITYP, XP , XP1, XP2, XP3, XP4, XP5, DIS, FP, STATU
  !
      implicit none
      !--- dummy variables
           type(SimBoxData_DEV)      :: theData
           integer,        intent(in):: NPRT, NAPDEV
           integer, device, dimension(:),        allocatable:: ITYP, STATU, GID
           real(KINDDF), device, dimension(:,:), allocatable:: XP , XP1, XP2, XP3, XP4, XP5, DIS, FP
           real(KINDDF), device, dimension(:),   allocatable:: EPOT, EKIN

      !--- Local vairables
           integer::CURDEV,ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(theData%IDEV)
             !$$---  Note the differencde between the size of GID, ITYP XP and the size of XP1... on GPU
              allocate(theData%ITYP(NPRT),   theData%XP(NPRT,3),  theData%XP1(NAPDEV,3), theData%DIS(NAPDEV,3), theData%FP(NAPDEV,3), &
                       theData%EPOT(NAPDEV), theData%EKIN(NAPDEV),theData%STATU(NAPDEV), theData%GID(NAPDEV),                         &
                       STAT=ERR)

              !if(allocated(hm_XP2)) then
              !   allocate(XP2(NAPDEV,3),XP3(NAPDEV,3),XP4(NAPDEV,3),XP5(NAPDEV,3),STAT=ERR)
              !end if

              !$$--- return back to origial device
              ERR = cudaSetDevice(CURDEV)
              if(ERR) then
    100          write(*,*) "MDPSCU Error in allocating working device memory on device",theData%IDEV
                 stop
              end if

       return
  end subroutine Allocate_Working_Variables_DevData
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Variables(SimBoxDev)
  !***  PURPOSE:   to allocate working space
  !
      implicit none
      !--- dummy variables
           type(MDDevConfig)      :: SimBoxDev
      !--- Local vairables
           integer::CURDEV, ERR, NPRT, NAPDEV, I
           logical::pinnedFlag


           call Clear_Working_Variables(SimBoxDev)

         !$$--- determine the maxma number of particles on each device
              NPRT = SimBoxDev%dm_NPRT

              if(SimBoxDev%NDEVICE .EQ. 1) then
                 NAPDEV = NPRT
              else
                 NAPDEV = (dble(NPRT)/dble(SimBoxDev%NDEVICE))*SimBoxDev%m_RNAPDEV
              end if
              SimBoxDev%m_NAPDEV = NAPDEV

         !$$--- allocate working space on host
              allocate(SimBoxDev%hm_GID(NPRT), SimBoxDev%hm_GIDINV(NPRT), STAT = ERR)
              if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_GID failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
              end if

               allocate(SimBoxDev%hm_ITYP(NPRT), STAT=ERR, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_ITYP failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               else
                  if (.not. pinnedFlag) then
                      write(*,fmt="(A)") 'MDPSCU Warnning; Pinned allocation failed of hm_ITYP'
                      call ONWARNING(gm_OnWarning)
                  end if
               end if

               allocate(SimBoxDev%hm_XP(NPRT,3), STAT=ERR, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_XP failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               else
                  if (.not. pinnedFlag) then
                      write(*,fmt="(A)") 'MDPSCU Warnning; Pinned allocation failed of hm_XP'
                      call ONWARNING(gm_OnWarning)
                  end if
               end if

               allocate(SimBoxDev%hm_STATU(NPRT), STAT=ERR, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_STATU failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               else
                  if (.not. pinnedFlag) then
                      write(*,fmt="(A)") 'MDPSCU Warnning; Pinned allocation failed of hm_STATU'
                      call ONWARNING(gm_OnWarning)
                  end if
               end if


               allocate(SimBoxDev%hm_XP1(NPRT,3), SimBoxDev%hm_DIS(NPRT,3), SimBoxDev%hm_FP(NPRT,3), SimBoxDev%hm_EPOT(NPRT), SimBoxDev%hm_EKIN(NPRT), STAT=ERR)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_XP1 failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               end if

              !---
              if(allocated(SimBoxDev%m_XP2)) then
                  allocate(SimBoxDev%hm_XP2(NPRT,3), SimBoxDev%hm_XP3(NPRT,3), SimBoxDev%hm_XP4(NPRT,3), SimBoxDev%hm_XP5(NPRT,3), STAT=ERR)
                   if (ERR) then
                       write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_XP2 failed. '
                       write(*,fmt="(A)")  '               Process to be stopped'
                       stop
                   end if
              end if

         !$$--- allocate working space on devices
              do I = 1, size(SimBoxDev%dm_SimBox)
                 call Allocate_Working_Variables_DevData(SimBoxDev%dm_SimBox(I), NPRT, NAPDEV) 
              end do
              return

    100       write(*,*) "MDPSCU Error in allocating working device memory on device",ERR
              stop

       return
  end subroutine Allocate_Working_Variables
  !****************************************************************************


  !****************************************************************************
  subroutine Clear_Working_Variables_DevData(theData)
  !***  PURPOSE:   to clear the working memories allocated before on a device
  !     INPUT;     theData
  !
  !     OUPUT:     theData, with allocated memories deallocated
  !
      implicit none
      !--- dummy variables
           type(SimBoxData_DEV)      :: theData

      !--- Local vairables
      integer::CURDEV,IERR

              IERR = cudaGetDevice(CURDEV)
       !---   clear the working device memory on the first device
              IERR = cudaSetDevice(theData%IDEV)

              if(allocated(theData%ITYP)) then
                 deallocate(theData%ITYP)
              end if

              if(allocated(theData%XP)) then
                 deallocate(theData%XP)
              end if

              if(allocated(theData%XP1)) then
                 deallocate(theData%XP1)
              end if

              if(allocated(theData%XP2)) then
                 deallocate(theData%XP2)
              end if

              if(allocated(theData%XP3)) then
                 deallocate(theData%XP3)
              end if

              if(allocated(theData%XP4)) then
                 deallocate(theData%XP4)
              end if

              if(allocated(theData%XP5)) then
                 deallocate(theData%XP5)
              end if

              if(allocated(theData%DIS)) then
                 deallocate(theData%DIS)
              end if

              if(allocated(theData%FP)) then
                 deallocate(theData%FP)
              end if

              if(allocated(theData%EPOT)) then
                 deallocate(theData%EPOT)
              end if

              if(allocated(theData%EKIN)) then
                 deallocate(theData%EKIN)
              end if

              if(allocated(theData%STATU)) then
                 deallocate(theData%STATU)
              end if

              if(allocated(theData%NAC)) then
                 deallocate(theData%NAC)
              end if

              if(allocated(theData%IA1th)) then
                 deallocate(theData%IA1th)
              end if
              IERR = cudaSetDevice(CURDEV)


              return
  end subroutine Clear_Working_Variables_DevData
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_Working_Variables( SimBoxDev )
  !***  PURPOSE:   to clear the working memories allocated before on all devices
  !     INPUT
  !
      implicit none
      !--- dummy variables
           type(MDDevConfig)      :: SimBoxDev
      !--- Local vairables
      integer::I

       !$$---   clear the working  memory on the host
              if(allocated(SimBoxDev%hm_GID)) then
                 deallocate(SimBoxDev%hm_GID)
              end if

              if(allocated(SimBoxDev%hm_GIDINV)) then
                 deallocate(SimBoxDev%hm_GIDINV)
              end if

              if(allocated(SimBoxDev%hm_ITYP)) then
                 deallocate(SimBoxDev%hm_ITYP)
              end if

              if(allocated(SimBoxDev%hm_XP)) then
                 deallocate(SimBoxDev%hm_XP)
              end if

              if(allocated(SimBoxDev%hm_XP1)) then
                 deallocate(SimBoxDev%hm_XP1)
              end if

              if(allocated(SimBoxDev%hm_XP2)) then
                 deallocate(SimBoxDev%hm_XP2)
              end if

              if(allocated(SimBoxDev%hm_XP3)) then
                 deallocate(SimBoxDev%hm_XP3)
              end if

              if(allocated(SimBoxDev%hm_XP4)) then
                 deallocate(SimBoxDev%hm_XP4)
              end if

              if(allocated(SimBoxDev%hm_XP5)) then
                 deallocate(SimBoxDev%hm_XP5)
              end if

              if(allocated(SimBoxDev%hm_DIS)) then
                 deallocate(SimBoxDev%hm_DIS)
              end if

              if(allocated(SimBoxDev%hm_FP)) then
                 deallocate(SimBoxDev%hm_FP)
              end if

              if(allocated(SimBoxDev%hm_EPOT)) then
                 deallocate(SimBoxDev%hm_EPOT)
              end if

              if(allocated(SimBoxDev%hm_EKIN)) then
                 deallocate(SimBoxDev%hm_EKIN)
              end if

              if(allocated(SimBoxDev%hm_STATU)) then
                 deallocate(SimBoxDev%hm_STATU)
              end if


       !$$---   clear the working device memory on the first device
         !$$--- allocate working space on devices
              do I = 1, size(SimBoxDev%dm_SimBox)
                 call Clear_Working_Variables_DevData(SimBoxDev%dm_SimBox(I)) 
              end do
              SimBoxDev%HasPart = mp_NOTPART

      return
  end subroutine Clear_Working_Variables
  !****************************************************************************


  end module MD_TYPEDEF_SimBox_DEV

