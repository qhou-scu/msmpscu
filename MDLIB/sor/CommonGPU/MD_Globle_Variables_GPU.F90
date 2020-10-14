  module MD_Globle_Variables_GPU
  !***  DESCRIPTION: this module is to define the device global variables to
  !                  be used in GPU calculations of MD. This version is adapted
  !                  from MD_Globle_Variables_10_GPU.cuf (HOU Qing, Mar, 2010),
  !                  with the major modification of adding device arrays for
  !                  multple-device calculations.
  !                  ______________________________________________________
  !                  HOU Qing, April, 2012
  ! **** HOSTORY:
  !       * April, 2012(HOU Qing):    Original version
  !
  !       * Oct 19  2016(HOU Qing): Add the array hm_GIDINV, for the convenice
  !                                  for calculating spring-force between MEB images
  !
  !       * June    2018(HOU Qing): Discard the instant variables m_FP, m_EPOT, m_EKIN,
  !                                  and some device_to_Dev0 routines
  !
  !       * Sep.    2018(HOU Qing): Add the varaible hm_NPPB, dxm_GID, 
  !                                 parameters mp_ArrayOp_2Power, mp_ArrayOp_Blocksize, mp_ArrayOp_Gridsize
  !                                 temperory device array dxm_Array2DSwap. 
  !                                 and GPU kernels for performing dot_product on GPU
  !
  !       * Sep.27, 2018(HOU Qing): Some codes were move to MD_MultiGPU_Basic.F90
  !                                 Add a new data type MDDEVWorkSpace, the dxm_XP..., in old version
  !                                 are replaced by the MDDEVWorkSpace type array dm_WorkSpace.
  !                                 Correspondingly, the dxm_XXX arraies in other modules using
  !                                 this module are replaced by arraies of extend data types
  !
  !       * Dec.,   2018(HOU Qing): Some codes were move to MD_MultiGPU_Basic.F90
  !                                 Appying the new data type defined in MD_MultiGPU_Basic.F90,
  !                                 DevMat_DF, DevVec_I, for the dxm_XP..., in old version
  ! 
  !       * Jan.,   2019(HOU Qing): add host members: STARTCELL, ENDCELL, ENDA, NPA in MDDEVWorkSpace
  !                                 for the conveniece of using features appled by MD_MultiGPU_Basic.F90
  !                                 We will use MDDEVWorkSpace%STARTCELL etc in place of m_STARTA in
  !                                 therefater coding. 
  !       * Jan.29, 2019(HOU Qing): 
  !                                 m_NAAP, the number of ACTIVATIVE atom in linked-cells, is added in 
  !                                 the working space. This feature could be improve the efficient of
  !                                 calculating neighbor-list.
  !       * Sep.23, 2019(HOU Qing): 
  !                                 restore the variable m_FP. Althoug m_FP is not an accumulation variables of time,
  !                                 we need restore force of the main process after the main process is interrupted by 
  !                                 ,for example, a damping process that may be issued in event-check. 
  !                                 Add: CopyFPFrom_Host1_to_Device, CopyFPFrom_Host2_to_Device
  !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_MultiGPU_Basic
   implicit none


      !--- The original configuration variables of simulation boxes
      !    before parttioning. They are related to the config hm_XXX
      !    after partition by index list hm_GID and hm_INVGOD to be
      !    defined later 
      !    
           integer, parameter::mp_NOTPART = 0
           integer, parameter::mp_HASPART = 1
           integer           ::hm_HasPart = mp_NOTPART                 !The flag indicating if the originalconfiguration
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
           integer::dm_NPRT     = 0                                    !The total number of particles in whole simulation system
           integer::dm_NAPB     = 0                                    !The number of particles in one box
           integer::dm_NGROUP   = 0                                    !The number of types of particles
           real(KINDDF), dimension(mp_mxGROUP), private::m_CM = 0.D0   !
           integer,      dimension(:),   allocatable::m_ITYP           !The type IDs of the particles
           real(KINDDF), dimension(:,:), allocatable::m_XP             !The positions of the particles
           real(KINDDF), dimension(:,:), allocatable::m_XP1            !The velocities of the particles
           real(KINDDF), dimension(:,:), allocatable::m_XP2            !The second order of XP, if Gear scheme to be used
           real(KINDDF), dimension(:,:), allocatable::m_XP3            !The third order  of XP, if Gear scheme to be used
           real(KINDDF), dimension(:,:), allocatable::m_XP4            !The fourth order of XP, if Gear scheme to be used
           real(KINDDF), dimension(:,:), allocatable::m_XP5            !The fifith order of XP, if Gear scheme to be used
           real(KINDDF), dimension(:,:), allocatable::m_DIS            !The displacements of particles
           real(KINDDF), dimension(:,:), allocatable::m_FP             !The forces on paricles that are grouped by cells
           integer,      dimension(:),   allocatable::m_STATU          !The status of particles

     !--- The cell IDS the devices will handle.
     !    They are be generated by calling Neighbor-list routines
           integer::m_STARTCELL(m_MXDEVICE)=(/0,1,2,3,4,5/)            !The cell ID of the first cell on a device
           integer::m_ENDCELL(m_MXDEVICE)  =(/0,1,2,3,4,5/)            !The last ID of the last cell on a device
           integer::m_STARTA(m_MXDEVICE)   = 0                         !The atom ID of the first atom handled by a device
           integer::m_ENDA(m_MXDEVICE)     = 0                         !The atom ID of the last  atom handled by a device
           integer::m_NPA(m_MXDEVICE)      = 0                         !The number of atoms acutally handled by a device 
                                                                       ! m_NPA = m_ENDA - m_STARTA + 1
           integer::m_NAPDEV               = 0                         !The max number of particles for a device handling

           real(KINDDF),private::m_RNAPDEV =1.1                        !The redundance of particle on a device
                                                                       ! m_NAPDEV = (number of particle/number of devices)*m_RNAPDEV

           integer, dimension(:), allocatable::hm_NAC                  !The number of particles in a cell, created in neighborlist calculations
           integer, dimension(:), allocatable::hm_NAAC                 !The number of ACTIVE particles in a cell, created in neighborlist calculations
                                                                       !Usally hm_NAAC=hm_NAC, whhen activation region is neabled hm_NAAC is probably not equal to hm_NAC.  
           integer, dimension(:), allocatable::hm_IA1th                !The subscrition index of the first particles in a cell

     !---  The working spaceing on host
     !
     !---  The copies of configurations of the boxes, however, the atoms
     !     are sorted by cells. Sorting are perforemed in the neighbore
     !     list calculation.
     !     SEE: MD_NeighborsList_GPU.F90
     !
     !     NOTE: We only need the size of GID, ITYP and XP to
     !           cover all particles in the whole simulation box, the
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
     !           The number of particles in cells hm_NAC (dx_NAC) and
     !           the index of starting particles in each cell hm_IA1th (dx_IA1th) are
     !           allocated on calling Set_BoxCell_Number in Initialize_NeighboreList_DEV.
     !           Their values are assigned only after calling Cal_NeighBoreList_DEV.
     !
     !    SEE:   Allocate_Working_Variables in this module
          integer,      dimension(:), allocatable::hm_GID              !Index transformation for the partitioned box to original box.
                                                                       !That is, the Ith particle in partitioned box is the hm_GIG(I)th
                                                                       !particle in the original box
                                                                       !hm_GID is created by neighbor list calculation.

          integer,      dimension(:), allocatable::hm_GIDINV           !The inverse index transformation:
                                                                       !The Ith particle in the orginal box is the hm_GIDINVth pariticle
                                                                       !in the partitioned box
                                                                       !grouped by cells. Created by neighbor list calculation.
          integer,      dimension(:),   allocatable::hm_ITYP           !The type ID of the particles that are grouped by cells
          real(KINDDF), dimension(:,:), allocatable::hm_XP             !The positions  of the particles that are grouped by cells
          integer,      dimension(:),   allocatable::hm_STATU          !The status of paricles that are grouped by cells
                                                                       !NOTE: remove the pinned property for  hm_ITYP, hm_XP and hm_STATU, (2019-01-05) 
           
          real(KINDDF), dimension(:,:), allocatable::hm_XP1            !The velocities of the particles that are grouped by cells
          real(KINDDF), dimension(:,:), allocatable::hm_XP2            !The second order of hm_XP
          real(KINDDF), dimension(:,:), allocatable::hm_XP3            !The third order of hm_XP
          real(KINDDF), dimension(:,:), allocatable::hm_XP4            !The fourth order of hm_XP
          real(KINDDF), dimension(:,:), allocatable::hm_XP5            !The fifith order of hm_XP
          real(KINDDF), dimension(:,:), allocatable::hm_DIS            !The displacements of paricles that are grouped by cells
          real(KINDDF), dimension(:,:), allocatable::hm_FP             !The forces on paricles that are grouped by cells
          real(KINDDF), dimension(:),   allocatable::hm_EPOT           !The potential of paricles that are grouped by cells
          real(KINDDF), dimension(:),   allocatable::hm_EKIN           !The kinetic energy of paricles that are grouped by cells
           


     !---  The type definition of configure on device 
          type::MDDEVWorkSpace
                integer                                :: NPRT     = 0 !The total number of particles
                integer                                :: NBOX     = 0 !The number of boxes
                integer                                :: NAPB     = 0 !The number of particles in a box
                integer                                :: NAPDEV   = 0 !The max number of particles deal wtih on a device
                integer                                :: NC       = 0 !The total number of linked-cells in all boxes
                integer,         dimension(3)          :: NCELL    = 0 !The number of linked-cells in X, Y, Z direction in one box
                integer,         dimension(m_MXDEVICE) :: STARTCELL= 0 !The cell ID of the first cell on a device
                integer,         dimension(m_MXDEVICE) :: ENDCELL  = 0 !The last ID of the last cell on a device
                integer,         dimension(m_MXDEVICE) :: STARTA   = 0 !The atom ID of the first ACTIVE atom on a device
                integer,         dimension(m_MXDEVICE) :: ENDA     = 0 !The atom ID of the last ACTIVE atom on a device
                integer,         dimension(m_MXDEVICE) :: NPA      = 0 !The number of ACTIVE atoms on a device

                !--- The device variables on each devices
                type(DevVec_DF), dimension(m_MXDEVICE) :: CM
                type(DevVec_I),  dimension(m_MXDEVICE) :: ITYP
                type(DevMat_DF), dimension(m_MXDEVICE) :: XP
                !type(DevMat_DF), dimension(m_MXDEVICE) :: XPSTP       !The XP increment in one time step 
                type(DevMat_DF), dimension(m_MXDEVICE) :: XP1
                type(DevMat_DF), dimension(m_MXDEVICE) :: XP2
                type(DevMat_DF), dimension(m_MXDEVICE) :: XP3
                type(DevMat_DF), dimension(m_MXDEVICE) :: XP4
                type(DevMat_DF), dimension(m_MXDEVICE) :: XP5
                type(DevMat_DF), dimension(m_MXDEVICE) :: DIS
                type(DevMat_DF), dimension(m_MXDEVICE) :: FP
                type(DevVec_DF), dimension(m_MXDEVICE) :: EPOT
                type(DevVec_DF), dimension(m_MXDEVICE) :: EKIN
                type(DevVec_I),  dimension(m_MXDEVICE) :: STATU
                type(DevVec_I),  dimension(m_MXDEVICE) :: GID
                type(DevVec_I),  dimension(m_MXDEVICE) :: GIDINV

                type(DevVec_I),  dimension(m_MXDEVICE) :: NAC          ! Number of atoms in a cell, including active and non-active atoms
                type(DevVec_I),  dimension(m_MXDEVICE) :: NAAC         ! Number of atoms in a cell, including active atoms only
                type(DevVec_I),  dimension(m_MXDEVICE) :: IA1th        ! ID of the first atom in a call
                type(DevVec_I),  dimension(m_MXDEVICE) :: IC           ! Cell ID for atoms in
          end type MDDEVWorkSpace
          type(MDDEVWorkSpace), target::dm_WorkSpace

  !--- interface of routines in this module -------------------

  !---------------------------------------------------------
          private:: Allocate_Working_Variables

  !---------------------------------------------------------
          private::Allocate_Working_Variables_template

  !---------------------------------------------------------
  !       to check if device have been initialized
          public:: Check_DEVICES

   !---------------------------------------------------------
          private:: Clear_Working_Variables

  !---------------------------------------------------------
  !       to clear the allocated global device memories and working memories
          public:: Clear_Globle_Variables_DEV

  !---------------------------------------------------------
          public:: CopyAllFrom_Devices_to_Host

  !---------------------------------------------------------
          public:: CopyAllFrom_Host_to_Devices

  !---------------------------------------------------------
  !       Note: the NOSHFIT version could be called before
  !             the system partioning
          private  :: COPY_X1_IN_NOSHIFT_template,   &
                      COPY_XN_IN_NOSHIFT_template,   &
                      COPY_I4D1_IN_NOSHIFT_template, &
                      COPY_I2D1_IN_NOSHIFT_template

          interface   COPY_IN_NOSHIFT_template
             module procedure COPY_X1_IN_NOSHIFT_template
             module procedure COPY_XN_IN_NOSHIFT_template
             module procedure COPY_I4D1_IN_NOSHIFT_template
             module procedure COPY_I2D1_IN_NOSHIFT_template
          end interface

  !---------------------------------------------------------
  !       Note: the SHFIT version must be called after
  !             the system partioning
          private  :: COPY_X1_IN_SHIFT_template,   &
                      COPY_XN_IN_SHIFT_template,   &
                      COPY_I4D1_IN_SHIFT_template, &
                      COPY_I2D1_IN_SHIFT_template
          interface   COPY_IN_SHIFT_template
             module procedure COPY_X1_IN_SHIFT_template
             module procedure COPY_XN_IN_SHIFT_template
             module procedure COPY_I4D1_IN_SHIFT_template
             module procedure COPY_I2D1_IN_SHIFT_template
          end interface

  !---------------------------------------------------------
  !       to copy the memory on each host to a device, hm_X, m_X
  !---------------------------------------------------------
          private  :: COPY_X1_OUT_NOSHIFT_template,   &
                      COPY_XN_OUT_NOSHIFT_template
          private  :: COPY_I4D1_OUT_NOSHIFT_template, &
                      COPY_I2D1_OUT_NOSHIFT_template
          private  :: COPY_SFX1_OUT_NOSHIFT_template,   &
                      COPY_SFXN_OUT_NOSHIFT_template
          interface   COPY_OUT_NOSHIFT_template
             module procedure COPY_X1_OUT_NOSHIFT_template
             module procedure COPY_XN_OUT_NOSHIFT_template
             module procedure COPY_SFX1_OUT_NOSHIFT_template
             module procedure COPY_SFXN_OUT_NOSHIFT_template
             module procedure COPY_I4D1_OUT_NOSHIFT_template
             module procedure COPY_I2D1_OUT_NOSHIFT_template
          end interface

  !---------------------------------------------------------
          private  :: COPY_X1_OUT_SHIFT_template,   &
                      COPY_XN_OUT_SHIFT_template
          private  :: COPY_I4D1_OUT_SHIFT_template, &
                      COPY_I2D1_OUT_SHIFT_template
          private  :: COPY_I4DN_OUT_SHIFT_template, &
                      COPY_I2DN_OUT_SHIFT_template
          private  :: COPY_SFX1_OUT_SHIFT_template, &
                      COPY_SFXN_OUT_SHIFT_template
          interface   COPY_OUT_SHIFT_template
             module procedure COPY_X1_OUT_SHIFT_template
             module procedure COPY_XN_OUT_SHIFT_template
             module procedure COPY_SFX1_OUT_SHIFT_template
             module procedure COPY_SFXN_OUT_SHIFT_template
             module procedure COPY_I4D1_OUT_SHIFT_template
             module procedure COPY_I2D1_OUT_SHIFT_template
             module procedure COPY_I4DN_OUT_SHIFT_template
             module procedure COPY_I2DN_OUT_SHIFT_template
          end interface COPY_OUT_SHIFT_template

  !---------------------------------------------------------
  !       the following routine devoted to copyin/out needed 
  !       physics quantities 
          private:: CopyDISFrom_Devices_to_Host0, &
                    CopyDISFrom_Devices_to_Host1
          public::  CopyDISFrom_Devices_to_Host
          interface CopyDISFrom_Devices_to_Host
             module procedure CopyDISFrom_Devices_to_Host0
             module procedure CopyDISFrom_Devices_to_Host1
          end interface
          public:: CopyDISFrom_Host_to_Devices

    !---------------------------------------------------------
          private:: CopyEPOTFrom_Devices_to_Host0, &
                    CopyEPOTFrom_Devices_to_Host1
          interface CopyEPOTFrom_Devices_to_Host
             module procedure CopyEPOTFrom_Devices_to_Host0
             module procedure CopyEPOTFrom_Devices_to_Host1
          end interface

          public::  CopyEPOTFrom_Host_to_Devices

  !---------------------------------------------------------
          private:: CopyFPFrom_Devices_to_Host0, &
                    CopyFPFrom_Devices_to_Host1, &
                    CopyFPFrom_Devices_to_Host2
          public::  CopyFPFrom_Devices_to_Host
          interface CopyFPFrom_Devices_to_Host
             module procedure CopyFPFrom_Devices_to_Host0
             module procedure CopyFPFrom_Devices_to_Host1
             module procedure CopyFPFrom_Devices_to_Host2
          end interface

          private:: CopyFPFrom_Host0_to_Devices, &
                    CopyFPFrom_Host1_to_Devices, &
                    CopyFPFrom_Host2_to_Devices
          public::  CopyFPFrom_Host_to_Devices
          interface CopyFPFrom_Host_to_Devices
             module procedure CopyFPFrom_Host0_to_Devices
             module procedure CopyFPFrom_Host1_to_Devices
             module procedure CopyFPFrom_Host2_to_Devices
          end interface

  !---------------------------------------------------------
          public:: CopyStatuFrom_Devices_to_Host
          public:: CopyStatuFrom_Host_to_Devices

  !---------------------------------------------------------
          private:: CopyXPFrom_Devices_to_Host0, &
                    CopyXPFrom_Devices_to_Host1, &
                    CopyXPFrom_Devices_to_Host2
          public::  CopyXPFrom_Devices_to_Host
          interface CopyXPFrom_Devices_to_Host
             module procedure CopyXPFrom_Devices_to_Host0
             module procedure CopyXPFrom_Devices_to_Host1
             module procedure CopyXPFrom_Devices_to_Host2
          end interface

          private:: CopyXPFrom_Host0_to_Devices, &
                    CopyXPFrom_Host1_to_Devices, &
                    CopyXPFrom_Host2_to_Devices
          public::  CopyXPFrom_Host_to_Devices
          interface CopyXPFrom_Host_to_Devices
             module procedure CopyXPFrom_Host0_to_Devices
             module procedure CopyXPFrom_Host1_to_Devices
             module procedure CopyXPFrom_Host2_to_Devices
          end interface

          public:: Synchroniz_XP_on_Devices

  !---------------------------------------------------------
          private:: CopyXP1From_Devices_to_Host0, &
                    CopyXP1From_Devices_to_Host1
          public::  CopyXP1From_Devices_to_Host
          interface CopyXP1From_Devices_to_Host
             module procedure CopyXP1From_Devices_to_Host0
             module procedure CopyXP1From_Devices_to_Host1
          end interface

          private:: CopyXP1From_Host0_to_Devices, &
                    CopyXP1From_Host1_to_Devices
          public:: CopyXP1From_Host_to_Devices
          interface CopyXP1From_Host_to_Devices
             module procedure CopyXP1From_Host0_to_Devices
             module procedure CopyXP1From_Host1_to_Devices
          end interface CopyXP1From_Host_to_Devices

          private:: ResetXP1_0,  &
                    ResetXP1_1
          public    ResetXP1          
          interface ResetXP1 
             module procedure ResetXP1_0
             module procedure ResetXP1_1    
          end interface ResetXP1                  

  !---------------------------------------------------------
  !----   extented operations for device data
  !
         !---------------------------------------------------
         !   to add two vector on devices
          private:: &
                    Add2d_noshift_template0_c,   &
                    Add2d_noshift_template1_c
          public::  DevAdd_noshift
          interface DevAdd_noshift
                    module procedure Add2d_noshift_template0_c
                    module procedure Add2d_noshift_template1_c 
          end interface DevAdd_noshift

          private:: &
                    AddBD2d_shift_template0_b, &
                    AddBD2d_shift_template1_b
          public::  DevAdd_shift
          interface DevAdd_shift
                    module procedure AddBD2d_shift_template0_b
                    module procedure AddBD2d_shift_template1_b
          end interface DevAdd_shift
          !---------------------------------------------------
          !   to calculate the dot-product of two vector on devices
          private:: &
                    Dot_product2d_noshift_template0_c
          public::  DevDot_noshift
          interface DevDot_noshift
                    module procedure Dot_product2d_noshift_template0_c
          end interface DevDot_noshift
         !---------------------------------------------------------
         !  to Minus two vector on devices
          private:: &
                    Minus2d_noshift_template0_c,  &
                    Minus2d_noshift_template1_c,  &
                    MinusVec_noshift_template0_c, &
                    MinusVec_noshift_template1_c
          public::  DevMinus_noshift
          interface DevMinus_noshift
                    module procedure Minus2d_noshift_template0_c
                    module procedure Minus2d_noshift_template1_c 

                    module procedure MinusVec_noshift_template0_c
                    module procedure MinusVec_noshift_template1_c 
          end interface DevMinus_noshift
         !---------------------------------------------------------
         !       to make a copy of vector on devices
          private:: &
                    MakeCopy2d_noshift_template0_c,  &
                    MakeCopyVec_noshift_template0_c
          public::  DevMakeCopy_noshift
          interface DevMakeCopy_noshift
                    module procedure MakeCopy2d_noshift_template0_c
                    module procedure MakeCopyVec_noshift_template0_c
          end interface DevMakeCopy_noshift
          !---------------------------------------------------
          !   to calculate the dot-product of two vector on devices
          private:: &
                    Normalize2d_noshift_template0_c
          public::  DevNormalize
          interface DevNormalize
                    module procedure Normalize2d_noshift_template0_c
          end interface DevNormalize
         !---------------------------------------------------------
         !       to calculate the scalar-product of vectors on devices
          private:: &
                   Scalar_product2d_noshift_template0_c, &
                   Scalar_product2d_noshift_template1_c
          public::  DevMultiply_noshift
          interface DevMultiply_noshift
                    module procedure Scalar_product2d_noshift_template0_c
                    module procedure Scalar_product2d_noshift_template1_c
          end interface DevMultiply_noshift
         !---------------------------------------------------------
         !  to extract the maxval of arries on devices
          private:: &
                    MaxAbsVal2d_noshift_template0_c,  &
                    MaxAbsValVec_noshift_template0_c
          public::  DevMaxAbsVal_noshift
          interface DevMaxAbsVal_noshift
                    module procedure MaxAbsVal2d_noshift_template0_c
                    module procedure MaxAbsValVec_noshift_template0_c
          end interface DevMaxAbsVal_noshift

          private:: &
                    MaxVal2d_noshift_template0_c
          public::  DevMaxVal_noshift
          interface DevMaxVal_noshift
                    module procedure MaxVal2d_noshift_template0_c
          end interface DevMaxVal_noshift          
          
  !---------------------------------------------------------
          private:: IndexToOriginalBox0_X1, &
                    IndexToOriginalBox0_XN
          private:: IndexToOriginalBox0_I1, &
                    IndexToOriginalBox0_IN
          private:: IndexToOriginalBox0_SFX1, &
                    IndexToOriginalBox0_SFXN
          public::  IndexToOriginalBox
          interface IndexToOriginalBox
             module procedure IndexToOriginalBox0_X1
             module procedure IndexToOriginalBox0_XN
             module procedure IndexToOriginalBox0_I1
             module procedure IndexToOriginalBox0_IN
             module procedure IndexToOriginalBox0_SFX1
             module procedure IndexToOriginalBox0_SFXN
          end interface

  !---------------------------------------------------------
  !       to initialize this module by allocating global device memories
          private:: Initialize_GB_A_DEV, &
                    Initialize_GB_B_DEV
          interface Initialize_Globle_Variables_DEV
             module procedure Initialize_GB_A_DEV
             module procedure Initialize_GB_B_DEV
          end interface

  
  !---------------------------------------------------------
  !       to allocate memories on all devices for storing linked-cell informations.
  !        called by Neighborelist calculation.
          private:: Set_BoxCell_Number0, &
                    Set_BoxCell_Number1
          public::  Set_BoxCell_Number
          interface Set_BoxCell_Number
             module procedure Set_BoxCell_Number0
             module procedure Set_BoxCell_Number1 
          end interface Set_BoxCell_Number

  contains

  !****************************************************************************
  subroutine Initialize_GB_B_DEV( SimBox, CtrlParam)
  !***  PURPOSE:  to initialize this module by allocating global device memories
  !               and working memories to be used
  !     INPUT:    SimBox,    the Simulation box
  !               CtrlParam, the control parameters
  !
  !     OUTPUT:
      implicit none
      !--- dummy variables
      type(SimMDBox)  ::SimBox
      type(SimMDCtrl) ::CtrlParam
      !--- Local vairables
      integer::ERR

         !---
              
              call Clear_Globle_Variables_DEV()
         !$$--- copy the simulation box to device memory
              dm_NPRT   = SimBox%NPRT
              dm_NAPB   = SimBox%NPRT
              dm_NGROUP = SimBox%NGROUP
               m_CM     = SimBox%CM

              !$$--- allocated host memory
              allocate(m_ITYP(dm_NPRT), m_XP(dm_NPRT,3), m_XP1(dm_NPRT,3),m_DIS(dm_NPRT,3), m_FP(dm_NPRT,3), &
                       m_STATU(dm_NPRT), STAT=ERR)
              if(ERR) goto 100

              !$$--- allocate memory if Nordic scheme to be used
              if(allocated(SimBox%XP2)) then
                 allocate(m_XP2(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox%XP3)) then
                 allocate(m_XP3(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox%XP4)) then
                 allocate(m_XP4(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox%XP5)) then
                 allocate(m_XP5(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              !$$--- allocated device memory
              call Allocate_Working_Variables()

             !$$--- prepair thw working space
              m_ITYP(1:dm_NPRT)     = SimBox%ITYP(1:dm_NPRT)
              m_XP(1:dm_NPRT,1:3)   = SimBox%XP(1:dm_NPRT,1:3)
              m_XP1(1:dm_NPRT,1:3)  = SimBox%XP1(1:dm_NPRT,1:3)
              m_DIS(1:dm_NPRT,1:3)  = SimBox%DIS(1:dm_NPRT,1:3)
              m_FP(1:dm_NPRT,1:3)   = SimBox%FP(1:dm_NPRT,1:3)
              m_STATU(1:dm_NPRT)    = SimBox%STATU(1:dm_NPRT)

              hm_HasPart = mp_NOTPART

              return

    100       write(*,*) "MDPSCU Error in allocating global device memory"
              stop
  end subroutine Initialize_GB_B_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_GB_A_DEV( SimBox, CtrlParam)
  !***  PURPOSE:  to initialize this module by allocating global device memories
  !               and working memories to be used
  !     INPUT:    SimBox,    the Simulation box
  !               CtrlParam, the control parameters
  !
  !     OUTPUT:
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      type(SimMDCtrl)             ::CtrlParam
      !--- Local vairables
      integer::ERR,  NB, NPRT, I, IS
         !---
              call Clear_Globle_Variables_DEV()

              NPRT      = SimBox(1)%NPRT
              NB        = size(SimBox)
              dm_NPRT   = NB*NPRT
              dm_NAPB   = NPRT
              dm_NGROUP = SimBox(1)%NGROUP
              m_CM     = SimBox(1)%CM
              !$$--- allocated host memory
              allocate(m_ITYP(dm_NPRT), m_XP(dm_NPRT,3),m_XP1(dm_NPRT,3),m_DIS(dm_NPRT,3),  m_FP(dm_NPRT,3), &
                       m_STATU(dm_NPRT), STAT=ERR)
              if(ERR) goto 100


              !$$--- allocate memory if Nordic scheme to be used
              if(allocated(SimBox(1)%XP2)) then
                 allocate(m_XP2(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox(1)%XP3)) then
                 allocate(m_XP3(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox(1)%XP4)) then
                 allocate(m_XP4(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              if(allocated(SimBox(1)%XP5)) then
                 allocate(m_XP5(dm_NPRT,3), STAT=ERR)
                    if(ERR) goto 100
              end if

              !$$--- allocated device memory
              call Allocate_Working_Variables()

              !$$--- copy the simulation box to working memory
              IS = 0
              do I=1, NB
                 m_ITYP(1+IS:NPRT+IS)    = SimBox(I)%ITYP(1:NPRT)
                 m_XP(1+IS:NPRT+IS,1:3)  = SimBox(I)%XP(1:NPRT,1:3)
                 m_XP1(1+IS:NPRT+IS,1:3) = SimBox(I)%XP1(1:NPRT,1:3)
                 m_DIS(1+IS:NPRT+IS,1:3) = SimBox(I)%DIS(1:NPRT,1:3)
                 m_FP(1+IS:NPRT+IS,1:3)  = SimBox(I)%FP(1:NPRT,1:3)
                 m_STATU(1+IS:NPRT+IS)   = SimBox(I)%STATU(1:NPRT)
                 IS = IS+NPRT
              end do

              hm_HasPart = mp_NOTPART
              return

    100       write(*,*) "MDPSCU Error in allocating global device memory"
              stop
  end subroutine Initialize_GB_A_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_Globle_Variables_DEV( )
  !***  PURPOSE:   to clear the allocated global device memories and working memories
  !     INPUT
  !
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
       !---
       !$$---   clear the device memory on the major device
              dm_NPRT  = 0
              dm_nGROUP = 0

        !$$--- dealoocate memory on host
              if(allocated(m_ITYP)) then
                 deallocate(m_ITYP)
              end if

              if(allocated(m_XP)) then
                 deallocate(m_XP)
              end if

              if(allocated(m_XP1)) then
                 deallocate(m_XP1)
              end if

              if(allocated(m_XP2)) then
                 deallocate(m_XP2)
              end if

              if(allocated(m_XP3)) then
                 deallocate(m_XP3)
              end if

              if(allocated(m_XP4)) then
                 deallocate(m_XP4)
              end if

              if(allocated(m_XP5)) then
                 deallocate(m_XP5)
              end if

              if(allocated(m_DIS)) then
                 deallocate(m_DIS)
              end if

              if(allocated(m_FP)) then
               deallocate(m_FP)
              end if

              if(allocated(m_STATU)) then
                 deallocate(m_STATU)
              end if
              call Clear_Working_Variables()

              return
  end subroutine Clear_Globle_Variables_DEV
  !****************************************************************************


  !****************************************************************************
  subroutine Allocate_Working_Variables_template(WorkSpace, NPRT, NAPDEV, NAPB)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     WorkSpace:  the working data on a device
  !
  !     OUTPUT:    WorkSpace,  the working data with allocated ITYP, XP , XP1, XP2, XP3, XP4, XP5, DIS, FP, STATU
  !
      implicit none
      !--- dummy variables
           type(MDDEVWorkSpace)::WorkSpace
           integer,        intent(in):: NPRT, NAPDEV, NAPB

      !--- Local vairables
           integer:: I

           !$$---  Note the differencde between the size of GID, ITYP XP and the size of XP1... on GPU
              WorkSpace%NPRT   = NPRT
              WorkSpace%NAPDEV = NAPDEV
              WorkSpace%NAPB   = NAPB
              WorkSpace%NBOX   = NPRT/NAPB
              call DevAllocate(WorkSpace%CM,     mp_mxGROUP, "DM_WS_CM")
              call DevAllocate(WorkSpace%ITYP,   NPRT,       "DM_WS_ITYP")
              !call DevAllocate(WorkSpace%XPSTP,(/NAPDEV,3/), "DM_WS_XPSTP")
              call DevAllocate(WorkSpace%XP,   (/NPRT,  3/), "DM_WS_XP")
              call DevAllocate(WorkSpace%XP1,  (/NAPDEV,3/), "DM_WS_XP1")
              call DevAllocate(WorkSpace%DIS,  (/NAPDEV,3/), "DM_WS_DIS")
              call DevAllocate(WorkSpace%FP,   (/NAPDEV,3/), "DM_WS_FP")
              call DevAllocate(WorkSpace%EPOT,   NAPDEV,     "DM_WS_EPOT")
              call DevAllocate(WorkSpace%EKIN,   NAPDEV,     "DM_WS_EKIN")
              call DevAllocate(WorkSpace%STATU,  NAPDEV,     "DM_WS_STATU")
              call DevAllocate(WorkSpace%GID,    NPRT,       "DM_WS_GID")
              call DevAllocate(WorkSpace%GIDINV, NPRT,       "DM_WS_GIDINV")
              call DevAllocate(WorkSpace%IC,     NAPDEV,     "DM_WS_IC")

              if(allocated(hm_XP2)) then
                 call DevAllocate(WorkSpace%XP2,  (/NAPDEV,3/), "DM_WS_XP2") 
                 call DevAllocate(WorkSpace%XP3,  (/NAPDEV,3/), "DM_WS_XP3") 
                 call DevAllocate(WorkSpace%XP4,  (/NAPDEV,3/), "DM_WS_XP4") 
                 call DevAllocate(WorkSpace%XP5,  (/NAPDEV,3/), "DM_WS_XP5") 
              end if
              do I=1, m_NDEVICE
                  call DevCopyIn(m_CM, 1, WorkSpace%CM(I), 1, mp_mxGROUP)
              end do    
       return
  end subroutine Allocate_Working_Variables_template
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Variables()
  !***  PURPOSE:   to allocate working space
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::CURDEV, ERR, NPRT, NAPDEV, NAPB
      logical::pinnedFlag


           call Clear_Working_Variables()

         !$$--- determine the maxma number of particles
         !$$    each device
              NPRT = dm_NPRT
              NAPB = dm_NAPB

              if(m_NDEVICE .EQ. 1) then
                 NAPDEV = NPRT
              else
                 NAPDEV = (dble(NPRT)/dble(m_NDEVICE))*m_RNAPDEV
              end if
              m_NAPDEV = NAPDEV

         !$$--- allocate working space on host
              allocate(hm_GID(NPRT), hm_GIDINV(NPRT), STAT = ERR)
              if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_GID failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
              end if

               allocate(hm_ITYP(NPRT), STAT=ERR) !, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_ITYP failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               !else
               !   if (.not. pinnedFlag) then
               !       write(*,fmt="(A)") 'MDPSCU Warnning; Pinned allocation failed of hm_ITYP'
               !       call ONWARNING(gm_OnWarning)
               !   end if
               end if

               allocate(hm_XP(NPRT,3), STAT=ERR) !, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_XP failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               !else
               !   if (.not. pinnedFlag) then
               !       write(*,fmt="(A)") 'MDPSCU Warnning; Pinned allocation failed of hm_XP'
               !       call ONWARNING(gm_OnWarning)
               !   end if
               end if

               allocate(hm_STATU(NPRT), STAT=ERR) !, PINNED=pinnedFlag)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_STATU failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               !else
               !   if (.not. pinnedFlag) then
               !       write(*,fmt="(A)") 'MDPSCU Warnning; Pinned allocation failed of hm_STATU'
               !       call ONWARNING(gm_OnWarning)
               !   end if
               end if


               allocate(hm_XP1(NPRT,3), hm_DIS(NPRT,3), hm_FP(NPRT,3), hm_EPOT(NPRT), hm_EKIN(NPRT), STAT=ERR)
               if (ERR) then
                  write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_XP1 failed. '
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               end if

              !---
              if(allocated(m_XP2)) then
                  allocate(hm_XP2(NPRT,3), hm_XP3(NPRT,3), hm_XP4(NPRT,3), hm_XP5(NPRT,3), STAT=ERR)
                   if (ERR) then
                       write(*,fmt="(A)")  ' MDPSCU Error: allocation of hm_XP2 failed. '
                       write(*,fmt="(A)")  '               Process to be stopped'
                       stop
                   end if
              end if

            !$$--- allocate working space on devices
              call Allocate_Working_Variables_template(dm_WorkSpace, NPRT, NAPDEV, NAPB)

       return
  end subroutine Allocate_Working_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_Working_Variables( )
  !***  PURPOSE:   to clear the working memories allocated before on all devices
  !     INPUT
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::I
      
       !$$---   clear the working  memory on the host
              if(allocated(hm_GID)) then
                 deallocate(hm_GID)
              end if

              if(allocated(hm_GIDINV)) then
                 deallocate(hm_GIDINV)
              end if

              if(allocated(hm_ITYP)) then
                 deallocate(hm_ITYP)
              end if

              if(allocated(hm_XP)) then
                 deallocate(hm_XP)
              end if

              if(allocated(hm_XP1)) then
                 deallocate(hm_XP1)
              end if

              if(allocated(hm_XP2)) then
                 deallocate(hm_XP2)
              end if

              if(allocated(hm_XP3)) then
                 deallocate(hm_XP3)
              end if

              if(allocated(hm_XP4)) then
                 deallocate(hm_XP4)
              end if

              if(allocated(hm_XP5)) then
                 deallocate(hm_XP5)
              end if

              if(allocated(hm_DIS)) then
                 deallocate(hm_DIS)
              end if

              if(allocated(hm_FP)) then
                 deallocate(hm_FP)
              end if

              if(allocated(hm_EPOT)) then
                 deallocate(hm_EPOT)
              end if

              if(allocated(hm_EKIN)) then
                 deallocate(hm_EKIN)
              end if

              if(allocated(hm_STATU)) then
                 deallocate(hm_STATU)
              end if

       !$$---   clear the working device memory on the devices
              call DevDeallocate(dm_WorkSpace%CM )
              call DevDeallocate(dm_WorkSpace%ITYP )
              !call DevDeallocate(dm_WorkSpace%XPSTP)
              call DevDeallocate(dm_WorkSpace%XP   )
              call DevDeallocate(dm_WorkSpace%XP1  )
              call DevDeallocate(dm_WorkSpace%XP2  )
              call DevDeallocate(dm_WorkSpace%XP3  )
              call DevDeallocate(dm_WorkSpace%XP4  )
              call DevDeallocate(dm_WorkSpace%XP5  )
              call DevDeallocate(dm_WorkSpace%DIS  )
              call DevDeallocate(dm_WorkSpace%FP   )
              call DevDeallocate(dm_WorkSpace%EPOT )
              call DevDeallocate(dm_WorkSpace%EKIN )
              call DevDeallocate(dm_WorkSpace%STATU)
              call DevDeallocate(dm_WorkSpace%GID)
              call DevDeallocate(dm_WorkSpace%GIDINV)
              call DevDeallocate(dm_WorkSpace%NAC)
              call DevDeallocate(dm_WorkSpace%NAAC)
              call DevDeallocate(dm_WorkSpace%IA1th)
              call DevDeallocate(dm_WorkSpace%IC)
              dm_WorkSpace%NPRT    = 0
              dm_WorkSpace%NAPDEV  = 0

       !$$---   Mark the system to be non-partitioned      
             hm_HasPart = mp_NOTPART

      return
  end subroutine Clear_Working_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Set_BoxCell_Number0(WorkSpace, NC, NCELL)
  !***  PURPOSE:   to allocate memories on all devices for storing cell
  !                information to be used for neighbore calculations.The
  !                range of cells to be deal with by the devices also
  !                are calculated.
  !
  !
      implicit none
      !--- dummy variables
        integer, intent(in)  ::NC
        integer, intent(in)  ::NCELL(:)
        type(MDDEVWorkSpace) ::WorkSpace
      !--- Local vairables

          !$$--- set the number cell on host on devices
              WorkSpace%NC        = NC
              WorkSpace%NCELL(1:3) = NCELL(1:3)
              call DevAllocate(WorkSpace%NAC,    NC,  "DM_WS_NAC")
              call DevAllocate(WorkSpace%NAAC,   NC,  "DM_WS_NAAC")
              call DevAllocate(WorkSpace%IA1th,  NC,  "DM_WS_IA1TH")
              call DevSet(WorkSpace%NAC,   0)
              call DevSet(WorkSpace%NAAC,  0)
              call DevSet(WorkSpace%IA1th, 0)
       return
  end subroutine Set_BoxCell_Number0
  !****************************************************************************

  !****************************************************************************
  subroutine Set_BoxCell_Number1(NC, NCELL)
   !***  PURPOSE:   to allocate memories on all devices for storing cell
   !                information to be used for neighbore calculations.The
   !                range of cells to be deal with by the devices also
   !                are calculated.
   !
   !
       implicit none
       !--- dummy variables
            integer, intent(in)  ::NC, NCELL(:)
       !--- Local vairables
            integer::IERR
 
            call Set_BoxCell_Number0(dm_WorkSpace, NC, NCELL)
            !$$-- set the number cell on host
               if(allocated(hm_NAC))    deallocate(hm_NAC)
               if(allocated(hm_NAAC))   deallocate(hm_NAAC)
               if(allocated(hm_IA1th))  deallocate(hm_IA1th)
               allocate(hm_NAC(NC),hm_NAAC(NC),hm_IA1th(NC),STAT=IERR)
               if (IERR) then
                   write(*,*) 'Allocation of hm_NAC and hm_IA1TH failed'
                   stop
               end if
               hm_NAC      = 0
               hm_NAAC     = 0 
               hm_IA1th    = 0

               m_STARTCELL = 0
               m_ENDCELL   =0
               m_ENDA      = 0 
               m_NPA       = 0 
        return
   end subroutine Set_BoxCell_Number1
   !****************************************************************************

  !****************************************************************************
  subroutine CopyAllFrom_Devices_to_Host()
  !***  PURPOSE:   To copy the configuration from devices to host.
  !
  !     NOTE:      Since the positions and other vairables on devices were grouped
  !                by cells. The configuration variables on devices will be converted
  !                to the original variables inwhole box
  !
  !     INPUT:     NONE
  !     OUTPUT:    hm_ITYP, hm_STATU, XP, XP1, XP2, XP3, XP4, XP5, DIS, FP
      implicit none
      !--- dummy variables


      !----   Local variables
          integer::I,  STARTA, ENDA, NA

           !$$---
            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy from device to host before partition the system"
               stop
            end if

           !$$*** copy the array on devices to host intermediate arraies , hm_XP...
            do I=1, m_NDEVICE
               !$$--- the first atom on the device
               STARTA = hm_IA1th(m_STARTCELL(I))

               !$$--- the last atom on the device
               ENDA = hm_IA1th(m_ENDCELL(I)) + hm_NAC(m_ENDCELL(I))-1

               !$$--- the number of atoms on the device
               NA = ENDA - STARTA + 1

               call DevCopyOut(hm_XP,     STARTA, dm_WorkSpace%XP(I),    STARTA, (/NA,3/)) 
               call DevCopyOut(hm_XP1,    STARTA, dm_WorkSpace%XP1(I),   1,      (/NA,3/)) 
               call DevCopyOut(hm_DIS,    STARTA, dm_WorkSpace%DIS(I),   1,      (/NA,3/)) 
               call DevCopyOut(hm_FP,     STARTA, dm_WorkSpace%FP(I),    1,      (/NA,3/)) 
               call DevCopyOut(hm_STATU,  STARTA, dm_WorkSpace%STATU(I), 1,        NA    ) 

               if(allocated(hm_XP2)) then
                  call DevCopyOut(hm_XP2, STARTA, dm_WorkSpace%XP2(I),   1,      (/NA,3/)) 
                  call DevCopyOut(hm_XP3, STARTA, dm_WorkSpace%XP3(I),   1,      (/NA,3/)) 
                  call DevCopyOut(hm_XP4, STARTA, dm_WorkSpace%XP4(I),   1,      (/NA,3/)) 
                  call DevCopyOut(hm_XP5, STARTA, dm_WorkSpace%XP5(I),   1,      (/NA,3/)) 
               end if
               !call COPYOUTALL_template(dm_WorkSpace(I)) 
            end do 

            return
   end subroutine CopyAllFrom_Devices_to_Host
  !****************************************************************************

  !****************************************************************************
  subroutine CopyAllFrom_Host_to_Devices()
  !***  PURPOSE:   To copy the configuration from host to devices.
  !
  !     NOTE:      Since the positions and other vairables on devices were grouped
  !                by cells. The configuration variables on devices will be converted
  !                to the original variables in whole box.
  !
  !                This routine should be called after m_GID generated by calling Neighborelist
  !
  !     INPUT:     NONE
  !     OUTPUT:    dITYP, dSTATU, dXP, dXP1, dXP2, dXP3, dXP4, dXP5, dDIS, dFP: the configures on devices
      implicit none
      !--- dummy variables


      !----   Local variables
            integer::I,  STARTA, ENDA, NA


           !$$*** copy the array on host to devices
            do I=1, m_NDEVICE
               !$$--- the first atom on the device
               STARTA = hm_IA1th(m_STARTCELL(I))

               !$$--- the last atom on the device
               ENDA = hm_IA1th(m_ENDCELL(I)) + hm_NAC(m_ENDCELL(I))-1

               !$$--- the number of atoms on the device
               NA = ENDA - STARTA + 1

               call DevCopyIn(hm_XP,     STARTA, dm_WorkSpace%XP(I),    STARTA, (/NA,3/)) 
               call DevCopyIn(hm_XP1,    STARTA, dm_WorkSpace%XP1(I),   1,      (/NA,3/)) 
               call DevCopyIn(hm_DIS,    STARTA, dm_WorkSpace%DIS(I),   1,      (/NA,3/)) 
               call DevCopyIn(hm_FP,     STARTA, dm_WorkSpace%FP(I),    1,      (/NA,3/)) 
               call DevCopyIn(hm_STATU,  STARTA, dm_WorkSpace%STATU(I), 1,        NA    ) 

               if(allocated(hm_XP2)) then
                  call DevCopyIn(hm_XP2, STARTA, dm_WorkSpace%XP2(I),   1,      (/NA,3/)) 
                  call DevCopyIn(hm_XP3, STARTA, dm_WorkSpace%XP3(I),   1,      (/NA,3/)) 
                  call DevCopyIn(hm_XP4, STARTA, dm_WorkSpace%XP4(I),   1,      (/NA,3/)) 
                  call DevCopyIn(hm_XP5, STARTA, dm_WorkSpace%XP5(I),   1,      (/NA,3/)) 
               end if
               !call COPYINALL_template(dm_WorkSpace(I)) 
            end do 

            return
   end subroutine CopyAllFrom_Host_to_Devices
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I4D1_IN_SHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy one-dimension integer array from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, device, dimension(:)::dI
           integer,         dimension(:)::hI
           integer, optional, intent(in)::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(dI(1), hI(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(dI(1), hI(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I4D1_IN_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I4D1_IN_NOSHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy one-dimension integer array from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hI,        the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dI,  the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           integer, device, dimension(:)::dI
           integer, dimension(:)::hI

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- the number of atoms on the device
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(dI(1), hI(1), dm_NPRT, stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(dI(1), hI(1), dm_NPRT)
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I4D1_IN_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I2D1_IN_SHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy one-dimension integer array from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           integer(2), device, dimension(:)::dI
           integer(2), dimension(:)::hI

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(dI(1), hI(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(dI(1), hI(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if
                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I2D1_IN_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I2D1_IN_NOSHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy one-dimension integer array from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer(2), device, dimension(:)::dI
           integer(2),         dimension(:)::hI
           integer, optional,  intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV,  NA, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  !--- the number of atoms on the device
                  NA = dm_NPRT

                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(dI(1), hI(1), NA, stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(dI(1), hI(1), NA)
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I2D1_IN_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_X1_IN_NOSHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy one-dimension arrayt from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:)::dXP
           real(KINDDF),         dimension(:)::hXP
           integer, optional,    intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, NA, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  !--- the TOTAL number of atoms on the device
                  NA = dm_NPRT

                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(dXP(1), hXP(1), NA, stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(dXP(1), hXP(1), NA)
                   end if
                  ERR = cudaSetDevice(CURDEV)
                 return

  end subroutine COPY_X1_IN_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_X1_IN_SHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy one-dimension arrayt from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:)::dXP
           real(KINDDF),         dimension(:)::hXP
           integer,    optional, intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(dXP(1), hXP(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(dXP(1), hXP(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_X1_IN_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_XN_IN_NOSHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a 3-dimension array from host to device.
   !               the atoms have been clustered by cells.
 !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:,:):: dXP
           real(KINDDF),         dimension(:,:):: hXP
           integer, optional,    intent(in)    :: CPYSTREAM

      !--- Local vairables
           integer::CURDEV, NA, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- the TOTAL number of atoms on the device
                  NA = dm_NPRT
                  N  = size(hXP, dim=2)

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     do I=1, N
                        ERR = cudaMemcpyAsync(dXP(1,I), hXP(1,I), NA, stream=CPYSTREAM)
                     end do
                   else
                     do I=1, N
                        ERR = cudaMemcpyAsync(dXP(1,I), hXP(1,I), NA)
                     end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_XN_IN_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_XN_IN_SHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a 3-dimension array from host to device.
  !               the atoms have been clustered by cells.
  !
  !     INPUT:     IDEV,      the ID of device
  !                hXP,       the variables on host
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     dXP, the copy on device
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:,:)::dXP
           real(KINDDF),         dimension(:,:)::hXP
           integer, optional,    intent(in)    ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- the number of atoms on the device
                  N  = size(hXP, dim=2)
                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     do I=1, N
                        ERR = cudaMemcpyAsync(dXP(1,I), hXP(m_NPA(IDEV),I), m_NPA(IDEV), stream=CPYSTREAM)
                     end do
                   else
                     do I=1, N
                        ERR = cudaMemcpyAsync(dXP(1,I), hXP(m_NPA(IDEV),I), m_NPA(IDEV))
                     end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_XN_IN_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I4D1_OUT_NOSHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy a one-dimension integer array from each device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dI,        the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hI,        the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, device,   dimension(:)::dI
           integer,           dimension(:)::hI
           integer, optional, intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I4D1_OUT_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I4D1_OUT_SHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy a one-dimension integer array from each device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dI,        the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hI,        the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, device,   dimension(:)::dI
           integer,           dimension(:)::hI
           integer, optional, intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(1), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(1), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I4D1_OUT_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I4DN_OUT_SHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension integer array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer,             intent(in)  ::IDEV
           integer*4, device, dimension(:,:)::dI
           integer*4,         dimension(:,:)::hI
           integer, optional, intent(in)    ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  N  = size(hI, dim=2)

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                      do I=1, N
                         ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV),I), dI(1, I), m_NPA(IDEV), stream=CPYSTREAM)
                      end do
                   else
                      do I=1, N
                         ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV),I), dI(1, I), m_NPA(IDEV))
                      end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_I4DN_OUT_SHIFT_template
  !****************************************************************************


  !*****************************************************************************
  subroutine COPY_I2D1_OUT_NOSHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension integer array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer,             intent(in)::IDEV
           integer*2, device, dimension(:)::dI
           integer*2,         dimension(:)::hI
           integer, optional, intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

           return
  end subroutine COPY_I2D1_OUT_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I2D1_OUT_SHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension integer array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer,             intent(in)::IDEV
           integer*2, device, dimension(:)::dI
           integer*2,         dimension(:)::hI
           integer, optional, intent(in)  ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(1), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV)), dI(1), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)
           return
   end subroutine COPY_I2D1_OUT_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_I2DN_OUT_SHIFT_template(IDEV, dI, hI, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension integer array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer,           intent(in)    ::IDEV
           integer*2, device, dimension(:,:)::dI
           integer*2,         dimension(:,:)::hI
           integer, optional, intent(in)    ::CPYSTREAM

      !--- Local vairables
           integer::CURDEV, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  N  = size(hI, dim=2)
                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                      do I=1, N
                         ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV),I), dI(1, I), m_NPA(IDEV), stream=CPYSTREAM)
                      end do
                   else
                      do I=1, N
                         ERR = cudaMemcpyAsync(hI(m_STARTA(IDEV),I), dI(1, I), m_NPA(IDEV))
                      end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

          return
   end subroutine COPY_I2DN_OUT_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_X1_OUT_NOSHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           real(KINDDF), device, dimension(:)::dXP
           real(KINDDF), dimension(:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

          return
  end subroutine COPY_X1_OUT_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_X1_OUT_SHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer,      intent(in)          ::IDEV
           integer,      optional, intent(in)::CPYSTREAM
           real(KINDDF), device, dimension(:)::dXP
           real(KINDDF),         dimension(:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(1), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(1), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

           return
  end subroutine COPY_X1_OUT_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_SFX1_OUT_NOSHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           real(KINDSF), device, dimension(:)::dXP
           real(KINDSF), dimension(:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(m_STARTA(IDEV)), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(m_STARTA(IDEV)), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_SFX1_OUT_NOSHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_SFX1_OUT_SHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a one dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           real(KINDSF), device, dimension(:)::dXP
           real(KINDSF), dimension(:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(1), m_NPA(IDEV), stream=CPYSTREAM)
                   else
                     ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV)), dXP(1), m_NPA(IDEV))
                   end if

                  ERR = cudaSetDevice(CURDEV)

                 return

  end subroutine COPY_SFX1_OUT_SHIFT_template
  !****************************************************************************

  !*****************************************************************************
  subroutine COPY_XN_OUT_NOSHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a three dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           real(KINDDF), device, dimension(:,:)::dXP
           real(KINDDF), dimension(:,:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  N  = size(hXP, dim=2)

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(m_STARTA(IDEV),I), m_NPA(IDEV), stream=CPYSTREAM)
                      end do
                   else
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(m_STARTA(IDEV),I), m_NPA(IDEV))
                      end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

           return
  end subroutine COPY_XN_OUT_NOSHIFT_template
  !*********************************************************************************

  !*****************************************************************************
  subroutine COPY_XN_OUT_SHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a three dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           real(KINDDF), device, dimension(:,:)::dXP
           real(KINDDF), dimension(:,:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  N  = size(hXP, dim=2)

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(1,I), m_NPA(IDEV), stream=CPYSTREAM)
                      end do
                   else
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(1,I), m_NPA(IDEV))
                      end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

           return
  end subroutine COPY_XN_OUT_SHIFT_template
  !*********************************************************************************

  !*****************************************************************************
  subroutine COPY_SFXN_OUT_NOSHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a three dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer,              intent(in)    ::IDEV
           integer, optional,    intent(in)    ::CPYSTREAM
           real(KINDSF), device, dimension(:,:)::dXP
           real(KINDSF),         dimension(:,:)::hXP

      !--- Local vairables
           integer::CURDEV, STARTA, ENDA, NA, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  N  = size(hXP, dim=2)

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(m_STARTA(IDEV),I), m_NPA(IDEV), stream=CPYSTREAM)
                      end do
                   else
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(m_STARTA(IDEV),I), m_NPA(IDEV))
                      end do
                   end if

                  ERR = cudaSetDevice(CURDEV)
         return
  end subroutine COPY_SFXN_OUT_NOSHIFT_template
  !*********************************************************************************

  !*****************************************************************************
  subroutine COPY_SFXN_OUT_SHIFT_template(IDEV, dXP, hXP, CPYSTREAM)
  !***  PURPOSE:  to copy a three dimension real array from a device to host.
  !               the atoms have been clustered by cells.
  !
  !
  !     INPUT:     IDEV,      the ID of device
  !                dXP,       the variables on device
  !                CPYSTREAM, optional, the device stream
  !
  !     OUTPUT     hXP,  the copy on host
  !
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, optional, intent(in)::CPYSTREAM
           real(KINDSF), device, dimension(:,:)::dXP
           real(KINDSF), dimension(:,:)::hXP

      !--- Local vairables
           integer::CURDEV, ERR, N, I

                  ERR = cudaGetDevice(CURDEV)
                  ERR = cudaSetDevice(m_DEVICES(IDEV))
                  N  = size(hXP, dim=2)

                  !--- NOTE the size dXP1 etc, are different from that of hm_XP1
                   if(present(CPYSTREAM)) then
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(1,I), m_NPA(IDEV), stream=CPYSTREAM)
                      end do
                   else
                      do I=1, N
                         ERR = cudaMemcpyAsync(hXP(m_STARTA(IDEV),I), dXP(1,I), m_NPA(IDEV))
                      end do
                   end if

                  ERR = cudaSetDevice(CURDEV)

           return
  end subroutine COPY_SFXN_OUT_SHIFT_template
  !*********************************************************************************

  !****************************************************************************
  subroutine Synchroniz_XP_on_Devices()
   !***  PURPOSE:  because different devices handled different segment of XP.
   !                this routine is to make XP on all devices to be same
   !
       implicit none
       !--- dummy variables
       !----   Local variables

            if(m_NDEVICE .gt. 1) then
                call CopyXPFrom_Devices_to_Host0()
                call SynchronizeDevices()
                call CopyXPFrom_Host0_to_Devices()
            endif      
          return
    end subroutine Synchroniz_XP_on_Devices
  !****************************************************************************  

  !****************************************************************************
  subroutine CopyXPFrom_Host0_to_Devices()
  !***  PURPOSE:  copy position on device0 to the arries on devices.
  !               NOTE: It is actually to copy the host version hm_XP to devices
  !                     One should note the difference between m_XP and hm_XP
  !
      implicit none
      !--- dummy variables


      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               !call COPY_XN_IN_NOSHIFT_template(I, dm_WorkSpace%XP(I)%Data, hm_XP)
               call DevCopyIn(hm_XP, dm_WorkSpace%XP(I), 3*dm_NPRT)
            end do

            return
   end subroutine CopyXPFrom_Host0_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXPFrom_Host1_to_Devices(hXP)
  !***  PURPOSE:  copy position on host to devices 
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hXP 
      !----   Local variables

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP from host to devices before partition the system"
               stop
            end if
            hm_XP(1:dm_NPRT,1:3)   = hXP(hm_GID(1:dm_NPRT),1:3)
            call CopyXPFrom_Host0_to_Devices()

           return
   end subroutine CopyXPFrom_Host1_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXPFrom_Host2_to_Devices(hXP)
  !***  PURPOSE:  copy position in one dimension on host to devices
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:)::hXP
      !----   Local variables
      integer::I, J

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP from host to devices before partition the system"
               stop
            end if

            do I=1, dm_NPRT
               J = (hm_GID(I)-1)*3
               hm_XP(I,1:3) = hXP( J+1:J+3)
            end do
            call CopyXPFrom_Host0_to_Devices()
           return
   end subroutine CopyXPFrom_Host2_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXPFrom_Devices_to_Host0()
  !***  PURPOSE:  copy position on devices to the arries on device0.
  !
      implicit none
      !--- dummy variables
      !--- Local variables
      integer::I

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP from device to host before partition the system"
               stop
            end if

           !*** copy the array on device to host
            do I=1, m_NDEVICE
               call DevCopyOut(hm_XP, m_STARTA(I), dm_WorkSpace%XP(I), m_STARTA(I), (/m_NPA(I),3/) )
            end do
            return
   end subroutine CopyXPFrom_Devices_to_Host0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXPFrom_Devices_to_Host1(hXP)
  !***  PURPOSE:  copy position on devices to the arries on host, with the atom order restored.
  !               this subroutine could be used when thermalization is performed on host.
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hXP
      !----   Local variables

            call CopyXPFrom_Devices_to_Host0()
            call SynchronizeDevices()

            hXP(1:dm_NPRT,1:3)   = hm_XP(hm_GIDINV(1:dm_NPRT),1:3)
           return
   end subroutine CopyXPFrom_Devices_to_Host1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXPFrom_Devices_to_Host2(hXP)
  !***  PURPOSE:  copy position on devices to the arries on host, with the atom order restored.
  !               this subroutine could be used when thermalization is performed on host.
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:)::hXP
      !----   Local variables
          integer:: I, J

            call CopyXPFrom_Devices_to_Host0()
            call SynchronizeDevices()

            do I=1, dm_NPRT
               J = (I-1)*3
               hXP(J+1:J+3) = hm_XP(hm_GIDINV(I),1:3)
            end do
           return
   end subroutine CopyXPFrom_Devices_to_Host2
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXP1From_Host0_to_Devices()
  !***  PURPOSE:  copy velocity on device0 (the host version hm_XP1 to the arries on devices.
  !               this subroutine could be used when thermalization is performed on host.
  !
  !
      implicit none
      !--- dummy variables

      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyIn(hm_XP1, m_STARTA(I), dm_WorkSpace%XP1(I), 1, (/m_NPA(I),3/) )
            end do

            return
   end subroutine CopyXP1From_Host0_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXP1From_Host1_to_Devices(hXP1)
  !***  PURPOSE:  copy velocity on device0 (the host version hm_XP1 to the arries on devices.
  !               this subroutine could be used when thermalization is performed on host.
  !
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hXP1 
      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyIn(hm_XP1, m_STARTA(I), dm_WorkSpace%XP1(I), 1, (/m_NPA(I),3/) )
            end do

      !----   Local variables

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP1 from host to devices before partition the system"
               stop
            end if
            hm_XP1(1:dm_NPRT,1:3)   = hXP1(hm_GID(1:dm_NPRT),1:3)
            call CopyXP1From_Host0_to_Devices()
            return
   end subroutine CopyXP1From_Host1_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine ResetXP1_0(WorkSpace)
  !***  PURPOSE:  to reset a new val to the velocity
  !
   implicit none
      !--- dummy variables
      type(MDDEVWorkSpace) ::WorkSpace
      !----   Local variables
           call DevSet(WorkSpace%XP1, 0.D0)
           return
  end subroutine ResetXP1_0
  !****************************************************************************

  !****************************************************************************
  subroutine ResetXP1_1()
  !***  PURPOSE:  to reset a new val to the velocity
  !
   implicit none
      !--- dummy variables
      !----   Local variables
           call ResetXP1_0(dm_WorkSpace)
           return
  end subroutine ResetXP1_1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXP1From_Devices_to_Host0()
  !***  PURPOSE:  copy velocity on device to the arries on device0.
  !               this subroutine could be used when thermalization is performed on host.
  !
      implicit none
      !--- dummy variables
      !--- Local variables
          integer::I, STARTA, ENDA, NA

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP1 from device to host before partition the system"
               stop
            end if

           !*** copy the array on devices to host
            do I=1, m_NDEVICE
               call DevCopyOut(hm_XP1, m_STARTA(I), dm_WorkSpace%XP1(I), 1, (/m_NPA(I),3/) )
            end do

            return
   end subroutine CopyXP1From_Devices_to_Host0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyXP1From_Devices_to_Host1(hXP)
  !***  PURPOSE:  copy position on devices to the arries on host, with the atom order restored.
  !               this subroutine could be used when thermalization is performed on host.
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hXP 
      !----   Local variables

            call CopyXP1From_Devices_to_Host0()
            call SynchronizeDevices()

            hXP(1:dm_NPRT,1:3)   = hm_XP1(hm_GIDINV(1:dm_NPRT),1:3)
           return
   end subroutine CopyXP1From_Devices_to_Host1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyFPFrom_Host0_to_Devices()
  !***  PURPOSE:  copy force on device0 to the arries on devices.
  !
      implicit none
      !--- dummy variables
      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyIn(hm_FP, m_STARTA(I), dm_WorkSpace%FP(I), 1, (/m_NPA(I),3/) )
            end do

            return
   end subroutine CopyFPFrom_Host0_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyFPFrom_Host1_to_Devices(hFP)
  !***  PURPOSE:  copy force on host to devices 
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hFP 
      !----   Local variables

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy FP from host to devices before partition the system"
               stop
            end if
            hm_FP(1:dm_NPRT,1:3)   = hFP(hm_GID(1:dm_NPRT),1:3)
            call CopyFPFrom_Host0_to_Devices()

           return
   end subroutine CopyFPFrom_Host1_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyFPFrom_Host2_to_Devices(hFP)
  !***  PURPOSE:  copy force in one dimension on host to devices
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:)::hFP
      !----   Local variables
      integer::I, J

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP from host to devices before partition the system"
               stop
            end if

            do I=1, dm_NPRT
               J = (hm_GID(I)-1)*3
               hm_FP(I,1:3) = hFP( J+1:J+3)
            end do
            call CopyFPFrom_Host0_to_Devices()
           return
   end subroutine CopyFPFrom_Host2_to_Devices
  !****************************************************************************   

  !****************************************************************************
  subroutine CopyFPFrom_Devices_to_Host0()
  !***  PURPOSE:  copy force from devices to devices.
  !
      implicit none
      !--- dummy variables
      !----   Local variables
       integer::I

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy FP from device to host before partition the system"
               stop
            end if

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyOut(hm_FP, m_STARTA(I), dm_WorkSpace%FP(I), 1, (/m_NPA(I),3/) )
            end do

           return
   end subroutine CopyFPFrom_Devices_to_Host0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyFPFrom_Devices_to_Host1(hFP)
  !***  PURPOSE:  copy force from devices to devices.
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hFP

      !----   Local variables
            call CopyFPFrom_Devices_to_Host0()
            call SynchronizeDevices()
            hFP(1:dm_NPRT,1:3) = hm_FP(hm_GIDINV(1:dm_NPRT),1:3)
            return
   end subroutine CopyFPFrom_Devices_to_Host1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyFPFrom_Devices_to_Host2(hFP)
  !***  PURPOSE:  copy force from devices to devices.
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:)::hFP

      !----   Local variables
      integer::I, J

            call CopyFPFrom_Devices_to_Host0()
            call SynchronizeDevices()

            do I=1, dm_NPRT
               J = (I-1)*3
               hFP(J+1:J+3) = hm_FP(hm_GIDINV(I),1:3)
            end do
            return
   end subroutine CopyFPFrom_Devices_to_Host2
  !****************************************************************************

  !****************************************************************************
  subroutine CopyStatuFrom_Host_to_Devices()
  !***  PURPOSE:  copy statu of atoms on device0 to the arries on devices.
  !               this subroutine could be used to modify statu of atom acoording
  !               to their current statu
  !     SEE:      the usage of this subroutine in
  !               MD_SimboxArray_AppShell_12_GPU
  !               MD_SimboxArray_AppShell_14_GPU
  !
  !
      implicit none
      !--- dummy variables


      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyIn(hm_STATU, m_STARTA(I), dm_WorkSpace%STATU(I), 1, m_NPA(I))
            end do
            return
   end subroutine CopyStatuFrom_Host_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyStatuFrom_Devices_to_Host()
  !***  PURPOSE:  copy statu on devices to the arries on host.
  !
  !
      implicit none
      !--- dummy variables
      !----   Local variables
       integer::I

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy STATU from device to host before partition the system"
               stop
            end if

           !*** copy the array on device to host
            do I=1, m_NDEVICE
               call DevCopyOut(hm_STATU, m_STARTA(I), dm_WorkSpace%STATU(I), 1, m_NPA(I))
            end do

            return
   end subroutine CopyStatuFrom_Devices_to_Host
  !****************************************************************************

  !****************************************************************************
  subroutine CopyEPOTFrom_Host_to_Devices()
  !***  PURPOSE:  copy potential on device0 to the arries on devices.
  !               this subroutine could be used in some analysis tools
  !               where, the forces are loaded from files
  !
      implicit none
      !--- dummy variables


      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyIn(hm_EPOT, m_STARTA(I), dm_WorkSpace%EPOT(I), 1,  m_NPA(I))
            end do

          return
   end subroutine CopyEPOTFrom_Host_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyEPOTFrom_Devices_to_Host0()
  !***  PURPOSE:  copy EPOT on devices to the arries on host.
  !
      implicit none
      !--- dummy variables
      !----   Local variables
      integer::I

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy EPOT from device to host before partition the system"
               stop
            end if

           !*** copy the array on devices to host
            do I=1, m_NDEVICE
               call DevCopyOut(hm_EPOT, m_STARTA(I), dm_WorkSpace%EPOT(I), 1,  m_NPA(I))
            end do

            return
   end subroutine CopyEPOTFrom_Devices_to_Host0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyEPOTFrom_Devices_to_Host1(hEPOT)
  !***  PURPOSE:  copy EPOT on devices to the arries on host.
  !
      implicit none
      !--- dummy variables
      real(KINDDF),dimension(:)::hEPOT
      !----   Local variables

            call CopyEPOTFrom_Devices_to_Host0()
            call SynchronizeDevices()
            hEPOT(1:dm_NPRT) = hm_EPOT(hm_GIDINV(1:dm_NPRT))

            return
   end subroutine CopyEPOTFrom_Devices_to_Host1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyDISFrom_Host_to_Devices()
  !***  PURPOSE:  copy velocity on device0 (the host version hm_XP1 to the arries on devices.
  !               this subroutine could be used when thermalization is performed on host.
  !
  !
      implicit none
      !--- dummy variables

      !----   Local variables
          integer::I

           !*** copy the array on host to device
            do I=1, m_NDEVICE
               call DevCopyIn(hm_DIS, m_STARTA(I), dm_WorkSpace%DIS(I), 1, (/m_NPA(I),3/))
            end do

           return
   end subroutine CopyDISFrom_Host_to_Devices
  !****************************************************************************

  !****************************************************************************
  subroutine CopyDISFrom_Devices_to_Host0()
  !***  PURPOSE:  copy velocity on device to the arries on device0.
  !               this subroutine could be used when thermalization is performed on host.
  !
      implicit none
      !--- dummy variables
      !--- Local variables
      integer::I

            if(hm_HasPart .lt. mp_HASPART) then
               write(*,*) "MDPSCU Error: copy XP1 from device to host before partition the system"
               stop
            end if

           !*** copy the array on devices to host
            do I=1, m_NDEVICE
               call DevCopyOut(hm_DIS, m_STARTA(I), dm_WorkSpace%DIS(I), 1, (/m_NPA(I),3/))
            end do
            return
   end subroutine CopyDISFrom_Devices_to_Host0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyDISFrom_Devices_to_Host1(hDIS)
  !***  PURPOSE:  copy position on devices to the arries on host, with the atom order restored.
  !               this subroutine could be used when thermalization is performed on host.
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::hDIS 
      !----   Local variables

            call CopyDISFrom_Devices_to_Host0()
            call SynchronizeDevices()

            hDIS(1:dm_NPRT,1:3)   = hm_DIS(hm_GIDINV(1:dm_NPRT),1:3)
           return
   end subroutine CopyDISFrom_Devices_to_Host1
  !****************************************************************************
  
  !****************************************************************************
  subroutine IndexToOriginalBox0_I1(From, To)
  !***  PURPOSE:  to order the array indices in partitioned system to the indice 
  !               original system
  !               
  !
      implicit none
      !--- dummy variables
      integer, dimension(:)::From 
      integer, dimension(:), allocatable::To
      !----   Local variables
      integer::NPRT

            NPRT = size(FROM)
            if(NPRT .ne. dm_NPRT) then
               write(*,fmt="(A)")           ' MDPSCU Error: size of array different from number of atoms in partitioned system'
               write(*,fmt="(A, I, A, I)")  '               ', NPRT, ' vs ', dm_NPRT
               write(*,fmt="(A)")           '               Process to be stopped'
            end if 
            
            if(.not. allocated(To)) then
               allocate(To(NPRT))
            end if
            To(1:dm_NPRT) = From(hm_GIDINV(1:dm_NPRT))
           return
   end subroutine IndexToOriginalBox0_I1
  !****************************************************************************

  !****************************************************************************
  subroutine IndexToOriginalBox0_IN(From, To)
  !***  PURPOSE:  to order the array indices in partitioned system to the indice 
  !               original system
  !               
  !
      implicit none
      !--- dummy variables
      integer, dimension(:,:)::From 
      integer, dimension(:,:), allocatable::To

      !----   Local variables
      integer::NPRT, DIM2

            NPRT = size(FROM, dim=1)
            DIM2 = size(FROM, dim=2)
            if(NPRT .ne. dm_NPRT) then
               write(*,fmt="(A)")           ' MDPSCU Error: size of array different from number of atoms in partitioned system'
               write(*,fmt="(A, I, A, I)")  '               ', NPRT, ' vs ', dm_NPRT
               write(*,fmt="(A)")           '               Process to be stopped'
            end if 
            
            if(.not. allocated(To)) then
               allocate(To(NPRT, DIM2))
            end if
            To(1:dm_NPRT,1:DIM2) = From(hm_GIDINV(1:dm_NPRT), 1:DIM2)
           return
   end subroutine IndexToOriginalBox0_IN
  !****************************************************************************

  !****************************************************************************
  subroutine IndexToOriginalBox0_X1(From, To)
  !***  PURPOSE:  to order the array indices in partitioned system to the indice 
  !               original system
  !               
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:)::From 
      real(KINDDF), dimension(:), allocatable::To
      !----   Local variables
      integer::NPRT

            NPRT = size(FROM)
            if(NPRT .ne. dm_NPRT) then
               write(*,fmt="(A)")           ' MDPSCU Error: size of array different from number of atoms in partitioned system'
               write(*,fmt="(A, I, A, I)")  '               ', NPRT, ' vs ', dm_NPRT
               write(*,fmt="(A)")           '               Process to be stopped'
            end if 
            
            if(.not. allocated(To)) then
               allocate(To(NPRT))
            end if
            To(1:dm_NPRT) = From(hm_GIDINV(1:dm_NPRT))
           return
   end subroutine IndexToOriginalBox0_X1
  !****************************************************************************

  !****************************************************************************
  subroutine IndexToOriginalBox0_XN(From, To)
  !***  PURPOSE:  to order the array indices in partitioned system to the indice 
  !               original system
  !               
  !
      implicit none
      !--- dummy variables
      real(KINDDF), dimension(:,:)::From 
      real(KINDDF), dimension(:,:), allocatable::To
      !----   Local variables
      integer::NPRT, DIM2

            NPRT = size(FROM, dim=1)
            DIM2 = size(FROM, dim=2)
            if(NPRT .ne. dm_NPRT) then
               write(*,fmt="(A)")           ' MDPSCU Error: size of array different from number of atoms in partitioned system'
               write(*,fmt="(A, I, A, I)")  '               ', NPRT, ' vs ', dm_NPRT
               write(*,fmt="(A)")           '               Process to be stopped'
            end if 
            
            if(.not. allocated(To)) then
               allocate(To(NPRT, DIM2))
            end if
            To(1:dm_NPRT,1:DIM2) = From(hm_GIDINV(1:dm_NPRT), 1:DIM2)
           return
   end subroutine IndexToOriginalBox0_XN
  !****************************************************************************

  !****************************************************************************
  subroutine IndexToOriginalBox0_SFX1(From, To)
  !***  PURPOSE:  to order the array indices in partitioned system to the indice 
  !               original system
  !               
  !
      implicit none
      !--- dummy variables
      real(KINDSF), dimension(:)::From 
      real(KINDSF), dimension(:), allocatable::To
      !----   Local variables
      integer::NPRT

            NPRT = size(FROM)
            if(NPRT .ne. dm_NPRT) then
               write(*,fmt="(A)")           ' MDPSCU Error: size of array different from number of atoms in partitioned system'
               write(*,fmt="(A, I, A, I)")  '               ', NPRT, ' vs ', dm_NPRT
               write(*,fmt="(A)")           '               Process to be stopped'
            end if 
            
            if(.not. allocated(To)) then
               allocate(To(NPRT))
            end if
            To(1:dm_NPRT) = From(hm_GIDINV(1:dm_NPRT))
           return
   end subroutine IndexToOriginalBox0_SFX1
  !****************************************************************************

  !****************************************************************************
  subroutine IndexToOriginalBox0_SFXN(From, To)
  !***  PURPOSE:  to order the array indices in partitioned system to the indice 
  !               original system
  !               
  !
      implicit none
      !--- dummy variables
      real(KINDSF), dimension(:,:)::From 
      real(KINDSF), dimension(:,:), allocatable::To
      !----   Local variables
      integer::NPRT, DIM2

            NPRT = size(FROM, dim=1)
            DIM2 = size(FROM, dim=2)
            if(NPRT .ne. dm_NPRT) then
               write(*,fmt="(A)")           ' MDPSCU Error: size of array different from number of atoms in partitioned system'
               write(*,fmt="(A, I, A, I)")  '               ', NPRT, ' vs ', dm_NPRT
               write(*,fmt="(A)")           '               Process to be stopped'
            end if 
            
            if(.not. allocated(To)) then
               allocate(To(NPRT, DIM2))
            end if
            To(1:dm_NPRT,1:DIM2) = From(hm_GIDINV(1:dm_NPRT), 1:DIM2)
           return
   end subroutine IndexToOriginalBox0_SFXN
  !****************************************************************************

  !**************************************************************************
   subroutine Dot_product2d_noshift_template0_c( V1, V2, PRDCT)
   !***  PURPOSE:  to calculate the dot between two array:: sum(V1(:,1:3)*V2(:,1:3))
   !               on GPU  
   !     INPUT:     V1, V2,    the two arraies for them dot-pruduct to be performed
   !     OUTPUT     PRDCT,     the dot_product of arraies handled by the block
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: V1, V2
          real(KINDDF)                  :: PRDCT
          !--- Local vairables
                call DevDot_product(m_NPA, V1, V2, PRDCT) 
           return
       end subroutine Dot_product2d_noshift_template0_c
  !***************************************************************************  
       
  !**************************************************************************
   subroutine Normalize2d_noshift_template0_c( V )
   !***  PURPOSE:  to calculate the dot between two array:: sum(V1(:,1:3)*V2(:,1:3))
   !               on GPU  
   !     INPUT:     V1, V2,    the two arraies for them dot-pruduct to be performed
   !     OUTPUT     PRDCT,     the dot_product of arraies handled by the block
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: V
          real(KINDDF)                  :: PRDCT
          !--- Local vairables
          real(KINDDF)::NORM
                call Dot_product2d_noshift_template0_c( V, V, NORM) 
                call Scalar_product2d_noshift_template0_c(1.D0/dsqrt(NORM), V)
           return
       end subroutine Normalize2d_noshift_template0_c
  !***************************************************************************       
       
  !**************************************************************************
   subroutine Scalar_product2d_noshift_template0_c(Scalar, V)
   !***  PURPOSE:  to scaled a vactor on GPU  
   !
   !     INPUT:     Scalar,    the scalar
   !                V,         the vecor
   !     OUTPUT     V,         the modified vector
      implicit none
          !----   DUMMY Variables
          real(KINDDF)                     :: Scalar
          type(DevMat_DF), dimension(:)    :: V
          !--- Local vairables
                call DevScalar_product(m_NPA, scalar, V) 
           return
   end subroutine Scalar_product2d_noshift_template0_c
  !***************************************************************************
       
  !**************************************************************************
   subroutine Scalar_product2d_noshift_template1_c(Scalar, V, Res)
   !***  PURPOSE:  to scaled a vactor on GPU  
   !
   !     INPUT:     Scalar,    the scalar
   !                V,         the vecor
   !     OUTPUT     V,         the modified vector
      implicit none
          !----   DUMMY Variables
          real(KINDDF)                     :: Scalar
          type(DevMat_DF), dimension(:)    :: V, Res
          !--- Local vairables
                call DevScalar_product(m_NPA, scalar, V, Res) 
           return
   end subroutine Scalar_product2d_noshift_template1_c
  !***************************************************************************   
     
  !**************************************************************************
   subroutine MakeCopy2d_noshift_template0_c(Sor, Tgt)
   !***  PURPOSE:  to make a copy of Sor on GPU  
   !
   !     INPUT:     Sor,       the sorce array
   !     OUTPUT     Tgt,       the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor, Tgt
          !--- Local vairables
               call DevMakeCopy(m_NPA, Sor, Tgt)
           return
   end subroutine MakeCopy2d_noshift_template0_c
  !***************************************************************************

  !**************************************************************************
   subroutine MakeCopyVec_noshift_template0_c(Sor, Tgt)
   !***  PURPOSE:  to make a copy of Sor on GPU  
   !
   !     INPUT:     Sor,       the sorce array
   !     OUTPUT     Tgt,       the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevVec_DF), dimension(:) :: Sor, Tgt
          !--- Local vairables
               call DevMakeCopy(m_NPA, Sor, Tgt)
           return
   end subroutine MakeCopyVec_noshift_template0_c
  !***************************************************************************   

  !**************************************************************************
   subroutine MaxAbsVal2d_noshift_template0_c(Sor, Res)
   !***  PURPOSE:  to extract max value of an array 
   !
   !     INPUT:     Sor,       the array
   !     OUTPUT     Res,       the max value
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor
          real(KINDDF)                  :: Res
          !--- Local vairables
               call DevMaxAbsVal(m_NPA, Sor, Res)
          return
   end subroutine MaxAbsVal2d_noshift_template0_c
  !*************************************************************************** 
   
  !**************************************************************************
   subroutine MaxAbsValVec_noshift_template0_c(Sor, Res)
      !***  PURPOSE:  to extract max value of an array 
      !
      !     INPUT:     Sor,       the array
      !     OUTPUT     Res,       the max value
      !
         implicit none
             !----   DUMMY Variables
             type(DevVec_DF), dimension(:) :: Sor
             real(KINDDF)                  :: Res
             !--- Local vairables
                  call DevMaxAbsVal(m_NPA, Sor, Res)
             return
      end subroutine MaxAbsValVec_noshift_template0_c
     !***************************************************************************   

  !**************************************************************************
   subroutine MaxVal2d_noshift_template0_c(Sor, Res)
   !***  PURPOSE:  to extract max value of an array 
   !
   !     INPUT:     Sor,       the array
   !     OUTPUT     Res,       the max value
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor
          real(KINDDF)                  :: Res
          !--- Local vairables
               call DevMaxVal(m_NPA, Sor, Res)
          return
   end subroutine MaxVal2d_noshift_template0_c
  !***************************************************************************

  !**************************************************************************
   subroutine Add2d_noshift_template0_c(Sor1, Sor2, Res)
   !***  PURPOSE:  to add matrix Sor1 and Sor
   !
   !     INPUT:     Sor1, Sor2   the sorce arraies
   !     OUTPUT     Tgt,         the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor1, Sor2, Res

          !--- Local vairables
               call DevAdd(m_NPA, Sor1, Sor2, Res)
          return
   end subroutine Add2d_noshift_template0_c
  !***************************************************************************
   
  !**************************************************************************
   subroutine Add2d_noshift_template1_c(Sor, Tgt)
   !***  PURPOSE:  to make a copy of Sor on GPU  
   !
   !     INPUT:     Sor,       the sorce array
   !     OUTPUT     Tgt,       the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor, Tgt

          !--- Local vairables
               call DevAdd(m_NPA, Sor, Tgt)
          return
   end subroutine Add2d_noshift_template1_c
  !***************************************************************************

  !**************************************************************************
   subroutine AddBD2d_shift_template0_b(LB, HB, Sor1, Sor2, Res)
   !***  PURPOSE:  to add matrix Sor1 and Sor
   !
   !     INPUT:     Sor1, Sor2   the sorce arraies
   !     OUTPUT     Tgt,         the tgt array
   !
      implicit none
          !----   DUMMY Variables
          real(KINDDF),    intent(in), dimension(:) :: LB, HB
          type(DevMat_DF),             dimension(:) :: Sor1, Sor2, Res

          !--- Local vairables
               call DevAdd(m_STARTA, m_NPA, LB, HB, Sor1, Sor2, Res)
          return
   end subroutine AddBD2d_shift_template0_b
  !***************************************************************************
   
  !**************************************************************************
   subroutine AddBD2d_shift_template1_b(LB, HB, Sor, Tgt)
   !***  PURPOSE:  to make a copy of Sor on GPU  
   !
   !     INPUT:     Sor,       the sorce array
   !     OUTPUT     Tgt,       the tgt array
   !
      implicit none
          !----   DUMMY Variables
          real(KINDDF),    intent(in), dimension(:) :: LB, HB
          type(DevMat_DF),             dimension(:) :: Sor, Tgt

          !--- Local vairables
               call DevAdd(m_STARTA, m_NPA, LB, HB, Sor, Tgt)
          return
   end subroutine AddBD2d_shift_template1_b
  !***************************************************************************   
   
  !**************************************************************************
   subroutine Minus2d_noshift_template0_c(Sor1, Sor2, Res)
   !***  PURPOSE:  to Minus matrix Sor1 and Sor
   !
   !     INPUT:     Sor1, Sor2   the sorce arraies
   !     OUTPUT     Tgt,         the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor1, Sor2, Res

          !--- Local vairables
               call DevMinus(m_NPA, Sor1, Sor2, Res)
          return
   end subroutine Minus2d_noshift_template0_c
  !***************************************************************************
   
  !**************************************************************************
   subroutine Minus2d_noshift_template1_c(Sor, Tgt)
   !***  PURPOSE:  to make a copy of Sor on GPU  
   !
   !     INPUT:     Sor,       the sorce array
   !     OUTPUT     Tgt,       the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevMat_DF), dimension(:) :: Sor, Tgt

          !--- Local vairables
               call DevMinus(m_NPA, Sor, Tgt)
          return
   end subroutine Minus2d_noshift_template1_c
  !***************************************************************************

  !**************************************************************************
   subroutine MinusVec_noshift_template0_c(Sor1, Sor2, Res)
   !***  PURPOSE:  to Minus vector Sor1 and Sor
   !
   !     INPUT:     Sor1, Sor2   the sorce arraies
   !     OUTPUT     Tgt,         the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevVec_DF), dimension(:) :: Sor1, Sor2, Res

          !--- Local vairables
               call DevMinus(m_NPA, Sor1, Sor2, Res)
          return
   end subroutine MinusVec_noshift_template0_c
  !***************************************************************************
   
  !**************************************************************************
   subroutine MinusVec_noshift_template1_c(Sor, Tgt)
   !***  PURPOSE:  to make a copy of Sor on GPU  
   !
   !     INPUT:     Sor,       the sorce array
   !     OUTPUT     Tgt,       the tgt array
   !
      implicit none
          !----   DUMMY Variables
          type(DevVec_DF), dimension(:) :: Sor, Tgt

          !--- Local vairables
               call DevMinus(m_NPA, Sor, Tgt)
          return
   end subroutine MinusVec_noshift_template1_c
  !***************************************************************************   
   
      
  !**************************************************************************
  end module MD_Globle_Variables_GPU

