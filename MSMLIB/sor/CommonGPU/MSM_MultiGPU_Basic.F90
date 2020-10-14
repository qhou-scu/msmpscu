  module MSM_MultiGPU_Basic
  !***  DESCRIPTION: this module is to define a number of functions for the use
  !                  of multiple GPU. This version is eaxtracted from old version
  !                  of MD_Globle_Variables_GPU.F90. 
  !                  ______________________________________________________
  !
  ! **** HOSTORY:
  !       * Sep.    2018(HOU Qing): Seperated from MD_Globle_Variables_GPU.F90, 
  !                                 Keep some generic functions 
  !
  !       * Dec.,   2018(HOU Qing): Add new data types for Device arrayes and their operations:
  !                                 DevVec_I, DevVec_DF, DevMat_I, DevMat_DF
  ! 
  !
  use cudafor
  #ifndef NODEVRAN
  use curand_device
  #endif
  use MSM_CONSTANTS
  implicit none

   !--- The list of parameters for GPU operations
          integer,parameter::m_MXDEVICE=6                             !The number limitation of devices to be used.
                                                                      !without the major device counted

          integer::m_NDEVICE=0                                        !The actual number of devices to be used
          integer::m_DEVICES(0:m_MXDEVICE)=(/0,0,1,2,3,4,5/)          !The index of devices to be used

          !---  the parameter for runing kernel using share memory
          integer, parameter, private::mp_ArrayOp_2Power      = 6
          integer, parameter, private::mp_ArrayOp_Blocksize   = 2**mp_ArrayOp_2Power
          integer, parameter, private::mp_ArrayOp_Gridsize    = 512
          !---  the parameter for runing kernel without share memory, for example Scalar_product
          !     or Add_
          integer, parameter, private::mp_BlockSize            = 256   !--- The GPU random number generator


          !---  the parameter for runing random number in kernels 
          !     Note:mp_RandOp_Blocksize*mp_RandOp_Gridsize cannot be larger than 64K.
          !          otherwise, there would be server error. 
          !          
          integer, parameter, private::mp_RandOp_Blocksize   = mp_ArrayOp_Blocksize
          integer, parameter, private::mp_RandOp_Gridsize    = mp_ArrayOp_Gridsize

        #ifdef NODEVRAN
          !--- uniform generator
          integer(kind=int_ptr_kind())::dm_URanGen(m_MXDEVICE) = 0
          !--- gaussian generator
          integer(kind=int_ptr_kind())::dm_GRanGen(m_MXDEVICE) = 0
        #else
          type::DevRandState
               integer::IDEV      = -1 
               integer::BLOCKSIZE = 0
               integer::GRIDSIZE  = 0
               integer::SEQ0      = 0  
               type(curandStateXORWOW),device, allocatable,dimension(:)::RandState
          end type DevRandState
          type(DevRandState), dimension(m_MXDEVICE)::gm_DevRandState
        #endif  

  !**********************************************************************************
  !   interface of routines for initialize the GPU
  !---------------------------------------------------------
  !       to check if device have been initialized
          public:: Check_DEVICES
  !---------------------------------------------------------
  !       to exit devices to release all resource on device
          public:: End_DEVICES
  !---------------------------------------------------------
  !       to set the number and index of devices to be used
          public:: Initialize_DEVICES
  !---------------------------------------------------------
  !       to synchronize all devices
          public:: SynchronizeDevices

  !   interface to random numbers
          public:: Archive_Rand_DEVICES
          public:: Initialize_Rand_DEVICES
          public:: Restore_Rand_DEVICES


  !***********************************************************************************
  !       some parameter and generic array operations
  !       
          type::DevVec_I
                character(len=32)                               ::Tag  =""
                integer                                         ::IDEV = -1
                integer,      device, dimension(:),  allocatable::Data
          end type DevVec_I

          type::DevVec_DF
                character(len=32)                               ::Tag  =""
                integer                                         ::IDEV = -1
                real(KINDDF), device, dimension(:),  allocatable::Data
          end type DevVec_DF

          type::DevMat_I
                character(len=32)                                 ::Tag  =""
                integer                                           ::IDEV = -1
                integer,      device, dimension(:,:),  allocatable::Data
          end type DevMat_I

          type::DevMat_DF
                character(len=32)                                 ::Tag  =""
                integer                                           ::IDEV = -1
                real(KINDDF), device, dimension(:,:),  allocatable::Data
          end type DevMat_DF

          !--- the data type used for storing temp. data ---------------
          private DeviceSwapData
          type, extends(DevVec_DF)::DeviceSwapData
                real(KINDDF), pinned ::RESDF(mp_ArrayOp_Gridsize)
          end type DeviceSwapData
          type(DeviceSwapData), dimension(m_MXDEVICE), private::dm_Array2DSwap
  !---------------------------------------------------------
          private:: & 
                    Allocate_DevMat_DF_0, &
                    Allocate_DevMat_I_0,  &
                    Allocate_DevVec_DF_0, &
                    Allocate_DevVec_I_0,  &

                    Allocate_DevMat_DF,   &
                    Allocate_DevMat_I,    &
                    Allocate_DevVec_DF,   &
                    Allocate_DevVec_I,    &

                    Allocate_DevMat_DF_b, &
                    Allocate_DevMat_I_b,  &
                    Allocate_DevVec_DF_b, &
                    Allocate_DevVec_I_b
          public::  DevAllocate
          interface DevAllocate
                    module procedure Allocate_DevMat_DF_0
                    module procedure Allocate_DevMat_I_0
                    module procedure Allocate_DevVec_DF_0
                    module procedure Allocate_DevVec_I_0

                    module procedure Allocate_DevMat_DF
                    module procedure Allocate_DevMat_I
                    module procedure Allocate_DevVec_DF
                    module procedure Allocate_DevVec_I

                    module procedure Allocate_DevMat_DF_b
                    module procedure Allocate_DevMat_I_b 
                    module procedure Allocate_DevVec_DF_b
                    module procedure Allocate_DevVec_I_b  
          end interface DevAllocate 
   !---------------------------------------------------------
          private:: &
                    Copyin_DevMat_DF_template0_a, &
                    Copyin_DevMat_DF_template0_b, &
                    Copyin_DevMat_DF_template0_c, &
                    Copyin_DevMat_DF_template0_d, &

                    Copyin_DevMat_I_template0_a,  &
                    Copyin_DevMat_I_template0_b,  &
                    Copyin_DevMat_I_template0_c,  &

                    Copyin_DevVec_DF_template0_a, &
                    Copyin_DevVec_DF_template0_c, &
                    Copyin_DevVec_DF_template0_d, &

                    Copyin_DevVec_I_template0_a,  &
                    Copyin_DevVec_I_template0_c,  &
                    Copyin_DevVec_I_template0_d,  &
                    Copyin_DevVec_I_template0_e
          public::  DevCopyIn
          interface DevCopyIn
                    module procedure CopyIn_DevMat_DF_template0_a
                    module procedure CopyIn_DevMat_DF_template0_b
                    module procedure CopyIn_DevMat_DF_template0_c
                    module procedure Copyin_DevMat_DF_template0_d

                    module procedure CopyIn_DevMat_I_template0_a
                    module procedure CopyIn_DevMat_I_template0_b
                    module procedure CopyIn_DevMat_I_template0_c

                    module procedure CopyIn_DevVec_DF_template0_a
                    module procedure CopyIn_DevVec_DF_template0_c
                    module procedure CopyIn_DevVec_DF_template0_d

                    module procedure CopyIn_DevVec_I_template0_a 
                    module procedure CopyIn_DevVec_I_template0_c  
                    module procedure CopyIn_DevVec_I_template0_d 
                    module procedure CopyIn_DevVec_I_template0_e
          end interface          
   !---------------------------------------------------------
          private:: &
                    Copyout_DevMat_DF_template0_a, &
                    Copyout_DevMat_DF_template0_b, &
                    Copyout_DevMat_DF_template0_c, &
                    Copyout_DevMat_DF_template0_d, &

                    Copyout_DevMat_I_template0_a,  &
                    Copyout_DevMat_I_template0_b,  &
                    Copyout_DevMat_I_template0_c,  &

                    Copyout_DevVec_DF_template0_a, &
                    Copyout_DevVec_DF_template0_c, &

                    Copyout_DevVec_I_template0_a,  &
                    Copyout_DevVec_I_template0_c,  &
                    Copyout_DevVec_I_template0_d
          public::  DevCopyOut
          interface DevCopyOut
                    module procedure Copyout_DevMat_DF_template0_a
                    module procedure Copyout_DevMat_DF_template0_b
                    module procedure Copyout_DevMat_DF_template0_c
                    module procedure Copyout_DevMat_DF_template0_d 

                    module procedure Copyout_DevMat_I_template0_a
                    module procedure Copyout_DevMat_I_template0_b
                    module procedure Copyout_DevMat_I_template0_c

                    module procedure Copyout_DevVec_DF_template0_a
                    module procedure Copyout_DevVec_DF_template0_c

                    module procedure Copyout_DevVec_I_template0_a 
                    module procedure Copyout_DevVec_I_template0_c
                    module procedure Copyout_DevVec_I_template0_d
          end interface          
   !---------------------------------------------------------
          private:: Deallocate_DevMat_DF,   &
                    Deallocate_DevMat_I,    &
                    Deallocate_DevVec_DF,   &
                    Deallocate_DevVec_I,    &

                    Deallocate_DevMat_DF_b, &
                    Deallocate_DevMat_I_b,  &
                    Deallocate_DevVec_DF_b, &
                    Deallocate_DevVec_I_b                  
          public::  DevDeallocate
          interface DevDeallocate
                    module procedure Deallocate_DevMat_DF
                    module procedure Deallocate_DevMat_I
                    module procedure Deallocate_DevVec_DF
                    module procedure Deallocate_DevVec_I

                    module procedure Deallocate_DevMat_DF_b
                    module procedure Deallocate_DevMat_I_b
                    module procedure Deallocate_DevVec_DF_b
                    module procedure Deallocate_DevVec_I_b
          end interface DevDeallocate                    
  !---------------------------------------------------------
          private:: AddScalar_DevVec_I_KERNEL0,        &
                    AddScalar_DevVec_I_template0,      &
                    AddScalar_DevVec_I_template0_a

          private:: Add_DevVec_DF_KERNEL0
          private:: Add_DevMat_DF_template0,   &
                    Add_DevMat_DF_template0_a, &
                    Add_DevMat_DF_template0_b, &
                    Add_DevMat_DF_template0_c, &

                    Add_DevMat_DF_template1_a, &
                    Add_DevMat_DF_template1_b, &
                    Add_DevMat_DF_template1_c

          private:: AddBD_DevVec_DF_KERNEL0
          private:: AddBD_DevMat_DF_template0,   &
                    AddBD_DevMat_DF_template0_a, &
                    AddBD_DevMat_DF_template0_b, &

                    AddBD_DevMat_DF_template1_a, &
                    AddBD_DevMat_DF_template1_b

          public::  DevAdd 
          interface DevAdd 
                    module procedure AddScalar_DevVec_I_template0                 
                    module procedure AddScalar_DevVec_I_template0_a                 
 
                    module procedure Add_DevMat_DF_template0
                    module procedure Add_DevMat_DF_template0_a
                    module procedure Add_DevMat_DF_template0_b
                    module procedure Add_DevMat_DF_template0_c 

                    module procedure Add_DevMat_DF_template1_a
                    module procedure Add_DevMat_DF_template1_b
                    module procedure Add_DevMat_DF_template1_c 

                    module procedure AddBD_DevMat_DF_template0
                    module procedure AddBD_DevMat_DF_template0_a
                    module procedure AddBD_DevMat_DF_template0_b

                    module procedure AddBD_DevMat_DF_template1_a
                    module procedure AddBD_DevMat_DF_template1_b
              end interface DevAdd     
              
  !---------------------------------------------------------
          private:: Minus_DevVec_DF_KERNEL0
          private:: Minus_DevMat_DF_template0,   &
                    Minus_DevMat_DF_template0_a, &
                    Minus_DevMat_DF_template0_b, &
                    Minus_DevMat_DF_template0_c, &

                    Minus_DevMat_DF_template1_a, &
                    Minus_DevMat_DF_template1_b, &
                    Minus_DevMat_DF_template1_c, &

                    Minus_DevVec_DF_template0,   &
                    Minus_DevVec_DF_template0_a, &
                    Minus_DevVec_DF_template0_b, &
                    Minus_DevVec_DF_template0_c, &

                    Minus_DevVec_DF_template1_a, &
                    Minus_DevVec_DF_template1_b, &
                    Minus_DevVec_DF_template1_c                    
          public::  DevMinus 
          interface DevMinus
                    module procedure Minus_DevMat_DF_template0
                    module procedure Minus_DevMat_DF_template0_a
                    module procedure Minus_DevMat_DF_template0_b
                    module procedure Minus_DevMat_DF_template0_c 

                    module procedure Minus_DevMat_DF_template1_a
                    module procedure Minus_DevMat_DF_template1_b
                    module procedure Minus_DevMat_DF_template1_c 

                    module procedure Minus_DevVec_DF_template0
                    module procedure Minus_DevVec_DF_template0_a
                    module procedure Minus_DevVec_DF_template0_b
                    module procedure Minus_DevVec_DF_template0_c 

                    module procedure Minus_DevVec_DF_template1_a
                    module procedure Minus_DevVec_DF_template1_b
                    module procedure Minus_DevVec_DF_template1_c 

              end interface DevMinus                            

  !---------------------------------------------------------

          private:: &
                    Dot_DevVec_DF_KERNEL0,   &
                    Dot_DevVec_DF_KERNEL1
          private:: &
                    Dot_DevMat_DF_template0,  &
                    Dot_DevMat_DF_template1,  &
                    Dot_DevMat_DF_template0_a,& 
                    Dot_DevMat_DF_template1_a,& 
                    Dot_DevMat_DF_template0_b,& 
                    Dot_DevMat_DF_template1_b,&
                    Dot_DevMat_DF_template0_c,& 
                    Dot_DevMat_DF_template1_c 
          public::  DevDot_product
          interface DevDot_product
                    module procedure Dot_DevMat_DF_template0
                    module procedure Dot_DevMat_DF_template1
                    module procedure Dot_DevMat_DF_template0_a
                    module procedure Dot_DevMat_DF_template1_a
                    module procedure Dot_DevMat_DF_template0_b
                    module procedure Dot_DevMat_DF_template1_b
                    module procedure Dot_DevMat_DF_template0_c
                    module procedure Dot_DevMat_DF_template1_c
              end interface DevDot_product
  !---------------------------------------------------------
          private:: Collect_SwapData    
          public::  Get_SumData,        &
                    Get_MaxData,        &
                    Get_MinData
  !---------------------------------------------------------
          private:: &
                    MakeCopy_DevVec_DF_KERNEL0,     &
                    MakeCopy_DevVec_DF_KERNEL1,     &
                    MakeCopy_DevVec_I_KERNEL0,      &
                    MakeCopy_DevVec_I_KERNEL1                    
          private:: &
                    MakeCopy_DevMat_DF_template0,   &
                    MakeCopy_DevMat_DF_template0_a, &
                    MakeCopy_DevMat_DF_template0_b, &
                    MakeCopy_DevMat_DF_template0_c, &
                    
                    MakeCopy_DevMat_DF_template1,   &
                    MakeCopy_DevMat_DF_template1_a, &
                    MakeCopy_DevMat_DF_template1_b, &
                    MakeCopy_DevMat_DF_template1_c, &

                    MakeCopy_DevVec_DF_template0,   &
                    MakeCopy_DevVec_DF_template0_a, &
                    MakeCopy_DevVec_DF_template0_b, &
                    MakeCopy_DevVec_DF_template0_c, &
                    
                    MakeCopy_DevVec_DF_template1,   &
                    MakeCopy_DevVec_DF_template1_a, &
                    MakeCopy_DevVec_DF_template1_b, &
                    MakeCopy_DevVec_DF_template1_c, &                    

                    MakeCopy_DevMat_I_template0,    &
                    MakeCopy_DevMat_I_template0_a,  &
                    MakeCopy_DevMat_I_template0_b,  &
                    MakeCopy_DevMat_I_template0_c,  &
                    
                    MakeCopy_DevMat_I_template1,    &
                    MakeCopy_DevMat_I_template1_a,  &
                    MakeCopy_DevMat_I_template1_b,  &
                    MakeCopy_DevMat_I_template1_c,  &

                    MakeCopy_DevVec_I_template0,    &
                    MakeCopy_DevVec_I_template0_a,  &
                    MakeCopy_DevVec_I_template0_b,  &
                    MakeCopy_DevVec_I_template0_c                  
          public::  DevMakeCopy
          interface DevMakeCopy
                    module procedure MakeCopy_DevMat_DF_template0
                    module procedure MakeCopy_DevMat_DF_template0_a
                    module procedure MakeCopy_DevMat_DF_template0_b
                    module procedure MakeCopy_DevMat_DF_template0_c
                    
                    module procedure MakeCopy_DevMat_DF_template1
                    module procedure MakeCopy_DevMat_DF_template1_a
                    module procedure MakeCopy_DevMat_DF_template1_b
                    module procedure MakeCopy_DevMat_DF_template1_c  

                    module procedure MakeCopy_DevVec_DF_template0
                    module procedure MakeCopy_DevVec_DF_template0_a
                    module procedure MakeCopy_DevVec_DF_template0_b
                    module procedure MakeCopy_DevVec_DF_template0_c
                    
                    module procedure MakeCopy_DevVec_DF_template1
                    module procedure MakeCopy_DevVec_DF_template1_a
                    module procedure MakeCopy_DevVec_DF_template1_b
                    module procedure MakeCopy_DevVec_DF_template1_c  
                     
                    module procedure MakeCopy_DevMat_I_template0
                    module procedure MakeCopy_DevMat_I_template0_a
                    module procedure MakeCopy_DevMat_I_template0_b
                    module procedure MakeCopy_DevMat_I_template0_c
                     
                    module procedure MakeCopy_DevMat_I_template1
                    module procedure MakeCopy_DevMat_I_template1_a
                    module procedure MakeCopy_DevMat_I_template1_b
                    module procedure MakeCopy_DevMat_I_template1_c  
                         
                    module procedure MakeCopy_DevVec_I_template0
                    module procedure MakeCopy_DevVec_I_template0_a
                    module procedure MakeCopy_DevVec_I_template0_b
                    module procedure MakeCopy_DevVec_I_template0_c

          end interface DevMakeCopy

  !---------------------------------------------------------
          private:: Max_DevMat_DF_KERNEL0
          private:: Max_DevMat_DF_template0,      &
                    Max_DevMat_DF_template0_a,    &
                    Max_DevMat_DF_template0_a13,  &
                    Max_DevMat_DF_template0_b,    &
                    Max_DevMat_DF_template0_b13,  &
                    Max_DevMat_DF_template0_c13
          public::  DevMaxVal 
          interface DevMaxVal
                    module procedure Max_DevMat_DF_template0_a 
                    module procedure Max_DevMat_DF_template0_a13
                    module procedure Max_DevMat_DF_template0_b  
                    module procedure Max_DevMat_DF_template0_b13
                    module procedure Max_DevMat_DF_template0_c13  
          end interface DevMaxVal
             
  !---------------------------------------------------------
          private:: MaxAbs_DevMat_DF_KERNEL0,        &
                    MaxAbs_DevVec_DF_KERNEL0
          private:: MaxAbs_DevMat_DF_template0,      &
                    MaxAbs_DevMat_DF_template0_a,    &
                    MaxAbs_DevMat_DF_template0_a13,  &
                    MaxAbs_DevMat_DF_template0_b,    &
                    MaxAbs_DevMat_DF_template0_b13,  &
                    MaxAbs_DevMat_DF_template0_c13,  &

                    MaxAbs_DevVec_DF_template0,      &
                    MaxAbs_DevVec_DF_template0_a,    &
                    MaxAbs_DevVec_DF_template0_b,    &
                    MaxAbs_DevVec_DF_template0_c                  
          public::  DevMaxAbsVal 
          interface DevMaxAbsVal
                    module procedure MaxAbs_DevMat_DF_template0_a 
                    module procedure MaxAbs_DevMat_DF_template0_a13
                    module procedure MaxAbs_DevMat_DF_template0_b  
                    module procedure MaxAbs_DevMat_DF_template0_b13
                    module procedure MaxAbs_DevMat_DF_template0_c13  

                    module procedure MaxAbs_DevVec_DF_template0_a 
                    module procedure MaxAbs_DevVec_DF_template0_b  
                    module procedure MaxAbs_DevVec_DF_template0_c
          end interface DevMaxAbsVal

  !---------------------------------------------------------
          private:: Normal2d_template0,   &
                    Normal2d_template0_a, &
                    Normal2d_template0_b, &
                    Normal2d_template0_c 
          public::  DevNormal
          interface DevNormal
                    module procedure Normal2d_template0
                    module procedure Normal2d_template0_a
                    module procedure Normal2d_template0_b
                    module procedure Normal2d_template0_c
          end interface DevNormal
          
  !---------------------------------------------------------
          private:: Scalar_DevVec_DF_KERNEL0
          private:: Scalar_DevMat_DF_template0,   &
                    Scalar_DevMat_DF_template0_a, &
                    Scalar_DevMat_DF_template0_b, &
                    Scalar_DevMat_DF_template0_c, &

                    Scalar_DevMat_DF_template1_a, &
                    Scalar_DevMat_DF_template1_b, &
                    Scalar_DevMat_DF_template1_c                    
          public::  DevScalar_product
          interface DevScalar_product
                    module procedure Scalar_DevMat_DF_template0
                    module procedure Scalar_DevMat_DF_template0_a
                    module procedure Scalar_DevMat_DF_template0_b
                    module procedure Scalar_DevMat_DF_template0_c

                    module procedure Scalar_DevMat_DF_template1_a
                    module procedure Scalar_DevMat_DF_template1_b
                    module procedure Scalar_DevMat_DF_template1_c
          end interface DevScalar_product             

  !---------------------------------------------------------
          private:: &
                    Sep_DevMat_DF_KERNEL0, &
                    Sep_DevMat_DF_KERNEL1
          private:: &
                    Sep_DevMat_DF_template0,   &
                    Sep_DevMat_DF_template0_1, &
                    Sep_DevMat_DF_template0_a, &
                    Sep_DevMat_DF_template0_b, &
                    Sep_DevMat_DF_template1,   &
                    Sep_DevMat_DF_template1_1, &
                    Sep_DevMat_DF_template1_a, &
                    Sep_DevMat_DF_template1_b
          public::  DevSep_product
          interface DevSep_product
                    module procedure Sep_DevMat_DF_template0
                    module procedure Sep_DevMat_DF_template0_1
                    module procedure Sep_DevMat_DF_template0_a
                    module procedure Sep_DevMat_DF_template0_b
                    module procedure Sep_DevMat_DF_template1
                    module procedure Sep_DevMat_DF_template1_1
                    module procedure Sep_DevMat_DF_template1_a
                    module procedure Sep_DevMat_DF_template1_b
              end interface DevSep_product
  !---------------------------------------------------------
  !---------------------------------------------------------
          private:: &
                    Set_DevVec_I_template0,   &
                    Set_DevVec_DF_template0,  &
                    Set_DevMat_I_template0,   &
                    Set_DevMat_DF_template0,  &

                    Set_DevVec_I_template0_b, &
                    Set_DevVec_DF_template0_b,&
                    Set_DevMat_I_template0_b, &
                    Set_DevMat_DF_template0_b
          public::  DevSet 
          interface DevSet
                    module procedure Set_DevVec_I_template0
                    module procedure Set_DevVec_DF_template0
                    module procedure Set_DevMat_I_template0
                    module procedure Set_DevMat_DF_template0

                    module procedure Set_DevVec_I_template0_b
                    module procedure Set_DevVec_DF_template0_b
                    module procedure Set_DevMat_I_template0_b
                    module procedure Set_DevMat_DF_template0_b 
          end interface DevSet                
  !---------------------------------------------------------

  contains
  !****************************************************************************
  subroutine Initialize_DEVICES( FIRSTDEV, NDEV)
  !***  PURPOSE:  to set the number and index of devices to be used
  !     INPUT     FIRSTDEV,  the first device ID used by this host process
  !               NDEV,      the number of devices involved in the host process
  !
      implicit none
      !--- dummy variables
      integer::FIRSTDEV,NDEV
      !--- Local vairables
      integer::ERR, K

             if(NDEV .GT. m_MXDEVICE) then
                write(*,*) "MDPSCU Error: the number of devices larger than permitted value ", m_MXDEVICE
                stop
              end if

              m_DEVICES(0) = FIRSTDEV
              m_NDEVICE    = NDEV
              m_DEVICES(1) = m_DEVICES(0)
              do K=2, m_NDEVICE
                 m_DEVICES(K) = m_DEVICES(K-1)+1
              end do

              !--- allocate memery for temp arrays
              do K=1, m_NDEVICE 
                  call DevAllocate(m_DEVICES(K), dm_Array2DSwap(K), mp_ArrayOp_Gridsize)
              end do
            return
  end subroutine Initialize_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine Check_DEVICES( )
  !***  PURPOSE:  to check if device have been initialized
  !     INPUT
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::ERR, K


               if(m_NDEVICE .le. 0) then
                  write(*,fmt="(A)")  ' MDPSCU Error: GPU version of MDPSCU is to be used. But devices are not initialized.'
                  write(*,fmt="(A)")  '               Process to be stopped'
                  stop
               end if
            return
  end subroutine Check_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine End_DEVICES( )
  !***  PURPOSE:  to exit devices to release all resource on device
  !     INPUT
  !
      use CudaRandomC2F_M
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::ERR, K

              !--- to destroy a sequence random number generator on devices
      #ifdef NODEVRAN
              do K=1, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(K) )
                 if(dm_URanGen(K) .gt. 0) ERR = curandDestroyGenerator(dm_URanGen(K))
                 if(dm_GRanGen(K) .gt. 0) ERR = curandDestroyGenerator(dm_GRanGen(K))
              end do
              dm_URanGen = 0
              dm_GRanGen = 0
      #else       
              do K=1, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(K) )
                 gm_DevRandState(K)%IDEV      = -1
                 gm_DevRandState(K)%BLOCKSIZE = 0
                 deallocate(gm_DevRandState(K)%RandState)
              end do
      #endif
              do K=0, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(K) )
                 ERR = cudaDeviceReset()
              end do

            return
  end subroutine End_DEVICES
  !****************************************************************************

  #ifndef NODEVRAN
  !****************************************************************************
  attributes(global) subroutine Initialize_Rand_Kernal(Seed, Seq0, DevRandState)
  !***  PURPOSE:  to intialize the uniform random number generator on devices
  !     INPUT     seed,         the first device ID used by this host process
  !               NS,           the number of serials
  !     OUTPUT:   DevRandState  the current state of random number genrator

   implicit none
    !---Dummy Vars---
    integer(8), value::Seed
    integer,    value::Seq0
    type(curandStateXORWOW), device, dimension(:)::DevRandState
     !---Local Vars---
     integer::IT,IB,IC, SEQ, OFFSET
     !---Body---
  
        IT     = (threadidx%y - 1)*blockdim%x + threadidx%x
        IB     = (blockidx%y  - 1)*griddim%x + blockidx%x
        IC     = (blockdim%x*blockdim%y)*(IB - 1) + IT
        SEQ    =  IC + SEQ0-1
        OFFSET = 0
        call curandInitXORWOW(Seed,SEQ, OFFSET, DevRandState(IC))
      return
  end subroutine Initialize_Rand_Kernal
  !****************************************************************************
   #endif
  !****************************************************************************
  subroutine Initialize_Rand_DEVICES()
  !***  PURPOSE:  to intialize the random numbers on devices
  !     INPUT     FIRSTDEV,  the first device ID used by this host process
  !               NDEV,      the number of devices involved in the host process
  !
      use RAND32SEEDLIB_MODULE
      use RAND32_MODULE
      use CudaRandomC2F_M
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::I, ERR, CURDEV
      integer(kind=8)::SEED
      integer(kind=4)::ISEED(2)
      equivalence(SEED, ISEED(1))
      type(dim3) :: blocks
      type(dim3) :: threads
!--- for rand number test
      !integer::N = 32*10000
      !real(kind=8), device, dimension(:),allocatable :: dRandomArray1
      !real(kind=8), device, dimension(:),allocatable :: dRandomArray2
      !real(kind=8), dimension(:),allocatable :: hRandomArray1
      !real(kind=8), dimension(:),allocatable :: hRandomArray2


              ERR = cudaGetDevice(CURDEV )
              ! --- to create a sequence random number generator on devices
             #ifdef NODEVRAN 
              do I=1, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(I) )
                 if(dm_URanGen(I) .gt. 0) ERR = curandDestroyGenerator(dm_URanGen(I))
                 ERR = curandCreateGenerator(dm_URanGen(I),CURAND_RNG_PSEUDO_XORWOW)
                 ERR = DRand32()*RAND32SEEDLIB_SIZE+1
                 call GetSeed_RAND32SEEDLIB(ERR, ISEED(1), ISEED(2))
                 ERR = curandSetPseudoRandomGeneratorSeed(dm_URanGen(I),SEED)

                 ! --- to crerate gaussian generator
                 if(dm_GRanGen(I) .gt. 0) ERR = curandDestroyGenerator(dm_GRanGen(I))
                 ERR = curandCreateGenerator(dm_GRanGen(I),CURAND_RNG_QUASI_SCRAMBLED_SOBOL64)
                 ERR = DRand32()*RAND32SEEDLIB_SIZE+1
                 call GetSeed_RAND32SEEDLIB(ERR, ISEED(1), ISEED(2))
                 ERR = curandSetPseudoRandomGeneratorSeed(dm_GRanGen(I),SEED)

              end do
              
              !---
            #else 
              ERR = DRand32()*RAND32SEEDLIB_SIZE+1
              call GetSeed_RAND32SEEDLIB(ERR, ISEED(1), ISEED(2))
              do I=1, m_NDEVICE
                 ERR = cudaSetDevice(m_DEVICES(I) )
                 gm_DevRandState(I)%IDEV      = m_DEVICES(I) 
                 gm_DevRandState(I)%BLOCKSIZE = mp_RandOp_Blocksize
                 gm_DevRandState(I)%GRIDSIZE  = mp_RandOp_Gridsize
                 gm_DevRandState(I)%SEQ0      = (I-1)*mp_RandOp_Gridsize*mp_RandOp_Blocksize + 1
                 allocate(gm_DevRandState(I)%RandState(mp_RandOp_Gridsize*mp_RandOp_Blocksize))
                 blocks  = dim3(mp_RandOp_Gridsize,  1, 1)
                 threads = dim3(mp_RandOp_Blocksize, 1, 1)
                 call Initialize_Rand_Kernal<<<blocks,threads>>>(SEED, gm_DevRandState(I)%SEQ0, gm_DevRandState(I)%RandState)
              end do
           #endif
              ERR = cudaSetDevice(CURDEV)
            return
  end subroutine Initialize_Rand_DEVICES
  !****************************************************************************

  !****************************************************************************
  subroutine Archive_Rand_DEVICES(hFile)
  !***  PURPOSE:  to save the random nuber status, that could be used for resoring calculation
  !     INPUT     hFile,  I/O unid
  !
      implicit none
      !--- dummy variables
      integer::hFile
      !--- Local vairables
      integer::I, ERR, CURDEV
      #ifdef NODEVRAN 
      #else 
         type(curandStateXORWOW), allocatable, dimension(:)::hRandState
      #endif
  !--- for rand number test

          ERR = cudaGetDevice(CURDEV )
              ! --- 
          #ifdef NODEVRAN 
                             
              !---
          #else 
              allocate(hRandState(mp_RandOp_Gridsize*mp_RandOp_Blocksize))
              write(hFile) m_NDEVICE
              do I=1, m_NDEVICE
                 write(hFile) gm_DevRandState(I)%IDEV
                 write(hFile) gm_DevRandState(I)%BLOCKSIZE
                 write(hFile) gm_DevRandState(I)%GRIDSIZE
                 write(hFile) gm_DevRandState(I)%SEQ0
                 
                 ERR = cudaSetDevice(gm_DevRandState(I)%IDEV)
                 hRandState(1:size(hRandState)) = gm_DevRandState(I)%RandState(1:size(hRandState))
                 write(hFile) hRandState
              end do
              deallocate(hRandState)
          #endif
          ERR = cudaSetDevice(CURDEV)
          return
  end subroutine Archive_Rand_DEVICES
  !****************************************************************************  

  !****************************************************************************
  subroutine Restore_Rand_DEVICES(hFile)
  !***  PURPOSE:  to save the random nuber status, that could be used for resoring calculation
  !     INPUT     hFile,  I/O unid
  !
      implicit none
      !--- dummy variables
      integer::hFile
      !--- Local vairables
      integer::I, ERR, CURDEV, NDEVICE, IDEV, BLOCKSIZE, GRIDSIZE, SEQ0
      #ifdef NODEVRAN 
      #else 
         type(curandStateXORWOW), allocatable, dimension(:)::hRandState
      #endif
  !--- for rand number test

          ERR = cudaGetDevice(CURDEV )
              ! --- 
          #ifdef NODEVRAN 
                             
              !---
          #else 
              allocate(hRandState(mp_RandOp_Gridsize*mp_RandOp_Blocksize))
              read(hFile) NDEVICE
              do I=1, NDEVICE
                 read(hFile) IDEV
                 read(hFile) BLOCKSIZE
                 read(hFile) GRIDSIZE
                 read(hFile) SEQ0
                 read(hFile) hRandState
                
                 if(gm_DevRandState(I)%IDEV .gt. 0) then
                    ERR = cudaSetDevice(gm_DevRandState(I)%IDEV)
                    gm_DevRandState(I)%RandState(1:size(hRandState)) = hRandState(1:size(hRandState))
                 end if   
              end do
              deallocate(hRandState)
          #endif
          ERR = cudaSetDevice(CURDEV)
          return
  end subroutine Restore_Rand_DEVICES
  !****************************************************************************  

  !****************************************************************************
  subroutine SynchronizeDevices()
  ! ***  PURPOSE:   to synchronize all devices
  ! 
      implicit none
      ! --- dummy variables
           integer::CURDEV,I,ERR

           ERR = cudaGetDevice(CURDEV)
           do I=0, m_NDEVICE
              err = cudaSetDevice(m_DEVICES(I))
              err = cudaThreadSynchronize();
           end do
           ERR = cudaSetDevice(CURDEV)

  end subroutine SynchronizeDevices
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevMat_DF_0(IDEV, Mat, TheDim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       integer                ::IDEV
       class(DevMat_DF)       ::Mat
       integer                ::TheDim(2)
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV
   
             call DevDeallocate(Mat)
             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(IDEV)   
             allocate(Mat%Data(TheDim(1),TheDim(2)), STAT=ERR)
             if(ERR .gt. 0) then
               if(present(Name)) then
                  write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU for "//Name(1:len_trim(Name))
               else
                  write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU"
               end if
               stop   
             end if 
                 
             Mat%IDEV = IDEV  
             if(present(Name)) Mat%Tag = Name
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Allocate_DevMat_DF_0
   !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevMat_DF(Mat, TheDim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_DF)       ::Mat
       integer                ::TheDim(2)
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV
   
             ERR = cudaGetDevice(CURDEV)
             if(present(Name)) then
                call Allocate_DevMat_DF_0(CURDEV, Mat, TheDim, Name)   
             else
                call Allocate_DevMat_DF_0(CURDEV, Mat, TheDim) 
             end if        
             return
  end subroutine Allocate_DevMat_DF
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevMat_DF_b(Mat, TheDim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_DF), dimension(:) :: Mat
       integer                        :: TheDim(2)
       character*(*), optional        ::Name
       !--- Local vairables
       integer::I
   
             do I=1, m_NDEVICE
                if(present(Name)) then
                   call Allocate_DevMat_DF_0(m_DEVICES(I), Mat(I), TheDim, Name)   
                else
                   call Allocate_DevMat_DF_0(m_DEVICES(I), Mat(I), TheDim) 
                end if 
            end do       
             return
  end subroutine Allocate_DevMat_DF_b
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevMat_I_0(IDEV, Mat, TheDim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       integer                ::IDEV
       class(DevMat_I)        ::Mat
       integer                ::TheDim(2)
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV
   
             call DevDeallocate(Mat)
             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(IDEV)   
             allocate(Mat%Data(TheDim(1),TheDim(2)), STAT=ERR)
             if(ERR .gt. 0) then
               if(present(Name)) then
                  write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU for "//Name(1:len_trim(Name))
               else
                  write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU"
               end if
               stop   
             end if 
                 
             Mat%IDEV = IDEV  
             if(present(Name)) Mat%Tag = Name
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Allocate_DevMat_I_0
   !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevMat_I(Mat, TheDim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_I)        ::Mat
       integer                ::TheDim(2)
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV
   
             ERR = cudaGetDevice(CURDEV)
             if(present(Name)) then
                call Allocate_DevMat_I_0(CURDEV, Mat, TheDim, Name)   
             else
                call Allocate_DevMat_I_0(CURDEV, Mat, TheDim) 
             end if        
             return
  end subroutine Allocate_DevMat_I
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevMat_I_b(Mat, TheDim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_I), dimension(:) :: Mat
       integer                       :: TheDim(2)
       character*(*), optional       ::Name
       !--- Local vairables
       integer::I
   
             do I=1, m_NDEVICE
                if(present(Name)) then
                   call Allocate_DevMat_I_0(m_DEVICES(I), Mat(I), TheDim, Name)   
                else
                   call Allocate_DevMat_I_0(m_DEVICES(I), Mat(I), TheDim) 
                end if 
            end do       
             return
  end subroutine Allocate_DevMat_I_b
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevVec_DF_0(IDEV, Vect, Dim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       integer                ::IDEV
       class(DevVec_DF)       ::Vect
       integer                ::Dim
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV
             
             call DevDeallocate(Vect)
             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(IDEV)
             
             allocate(Vect%Data(Dim), STAT=ERR)
             if(ERR .gt. 0) then
               if(present(Name)) then
                 write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU for "//Name(1:len_trim(Name))
               else
                  write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU"
               end if
               stop   
             end if 
             Vect%IDEV = IDEV  
             if(present(Name)) Vect%Tag = Name
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Allocate_DevVec_DF_0
  !****************************************************************************

  !****************************************************************************
  subroutine Allocate_DevVec_DF(Vect, Dim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevVec_DF)       ::Vect
       integer                ::Dim
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV

             ERR = cudaGetDevice(CURDEV)
             if(present(Name)) then
                call Allocate_DevVec_DF_0(CURDEV, Vect, Dim, Name)
             else 
               call Allocate_DevVec_DF_0(CURDEV, Vect, Dim)
             end if    
             
             return
   end subroutine Allocate_DevVec_DF
   !****************************************************************************

   !****************************************************************************
   subroutine Allocate_DevVec_DF_b(Vect, Dim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevVec_DF), dimension(:) :: Vect
       integer                        :: Dim
       character*(*), optional        :: Name
       !--- Local vairables
       integer::I
   
             do I=1, m_NDEVICE
                if(present(Name)) then
                   call Allocate_DevVec_DF_0(m_DEVICES(I), Vect(I), Dim, Name)   
                else
                   call Allocate_DevVec_DF_0(m_DEVICES(I), Vect(I), Dim) 
                end if 
            end do       
             return
  end subroutine Allocate_DevVec_DF_b
  !****************************************************************************

   !****************************************************************************
   subroutine Allocate_DevVec_I_0(Idev, Vect, Dim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       integer                ::Idev
       class(DevVec_I)        ::Vect
       integer                ::Dim
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV
             
             call DevDeallocate(Vect)
             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(Idev)
             allocate(Vect%Data(Dim), STAT=ERR)
             if(ERR .gt. 0) then
               if(present(Name)) then
                 write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU for "//Name(1:len_trim(Name))
               else
                  write(*,*) "Error:: fail allocating memory in MSM_MutltiGPU"
               end if
               stop   
             end if 
             Vect%IDEV = IDEV  
             if(present(Name)) Vect%Tag = Name
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Allocate_DevVec_I_0
   !****************************************************************************

   !****************************************************************************
   subroutine Allocate_DevVec_I(Vect, Dim, Name)
   !***  PURPOSE:  to allocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevVec_I)        ::Vect
       integer                ::Dim
       character*(*), optional::Name
       !--- Local vairables
       integer::ERR, CURDEV

             ERR = cudaGetDevice(CURDEV)
             if(present(Name)) then
                call Allocate_DevVec_I_0(CURDEV, Vect, Dim, Name)
             else
               call Allocate_DevVec_I_0(CURDEV, Vect, Dim)    
             end if  
             return
   end subroutine Allocate_DevVec_I
   !****************************************************************************

   !****************************************************************************
   subroutine Allocate_DevVec_I_b(Vect, Dim, Name)
      !***  PURPOSE:  to allocate memory for device data
      !     INPUT     
      !
          implicit none
          !--- dummy variables
          class(DevVec_I), dimension(:) :: Vect
          integer                        :: Dim
          character*(*), optional        :: Name
          !--- Local vairables
          integer::I
      
                do I=1, m_NDEVICE
                   if(present(Name)) then
                      call Allocate_DevVec_I_0(m_DEVICES(I), Vect(I), Dim, Name)   
                   else
                      call Allocate_DevVec_I_0(m_DEVICES(I), Vect(I), Dim) 
                   end if 
               end do       
                return
     end subroutine Allocate_DevVec_I_b
   !****************************************************************************

   !****************************************************************************
   subroutine Deallocate_DevMat_DF(Mat)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_DF)        ::Mat
       !--- Local vairables
       integer::ERR, CURDEV
             
             if(Mat%IDEV .lt. 0) return

             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(Mat%IDEV)
             if(allocated(Mat%Data)) then
                deallocate(Mat%Data, STAT=ERR)
                if(ERR .gt. 0) then
                   write(*,*) "Error:: fail deallocating memory in MSM_MutltiGPU for "//Mat%Tag(1:len_trim(Mat%Tag))
                   stop   
                end if   
             end if 
             Mat%IDEV = -1  
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Deallocate_DevMat_DF
   !****************************************************************************

   !****************************************************************************
   subroutine Deallocate_DevMat_DF_b(Mat)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_DF),dimension(:)::Mat
       !--- Local vairables
       integer::I
             
             do I=1, size(Mat)
                call Deallocate_DevMat_DF(Mat(I))
             end do
             return
   end subroutine Deallocate_DevMat_DF_b
   !****************************************************************************

   !****************************************************************************
   subroutine Deallocate_DevMat_I(Mat)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_I)        ::Mat
       !--- Local vairables
       integer::ERR, CURDEV
             
             if(Mat%IDEV .lt. 0) return

             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(Mat%IDEV)
             if(allocated(Mat%Data)) then
                deallocate(Mat%Data, STAT=ERR)
                if(ERR .gt. 0) then
                   write(*,*) "Error:: fail deallocating memory in MSM_MutltiGPU for "//Mat%Tag(1:len_trim(Mat%Tag))
                   stop   
                end if   
             end if 
             Mat%IDEV = -1  
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Deallocate_DevMat_I
   !****************************************************************************

   !****************************************************************************
   subroutine Deallocate_DevMat_I_b(Mat)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevMat_I),dimension(:)::Mat
       !--- Local vairables
       integer::I
             
             do I=1, size(Mat)
                call Deallocate_DevMat_I(Mat(I))
             end do
             return
   end subroutine Deallocate_DevMat_I_b
   !****************************************************************************
   
   !****************************************************************************
   subroutine Deallocate_DevVec_DF(Vect)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevVec_DF)        ::Vect
       !--- Local vairables
       integer::ERR, CURDEV
             
             if(Vect%IDEV .lt. 0) return

             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(Vect%IDEV)
             if(allocated(Vect%Data)) then
                deallocate(Vect%Data, STAT=ERR)
                if(ERR .gt. 0) then
                   write(*,*) "Error:: fail deallocating memory in MSM_MutltiGPU for "//Vect%Tag(1:len_trim(Vect%Tag))
                   stop   
                end if   
             end if 
             Vect%IDEV = -1 
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Deallocate_DevVec_DF
   !****************************************************************************

   !****************************************************************************
   subroutine Deallocate_DevVec_DF_b(Vect)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevVec_DF),dimension(:)::Vect
       !--- Local vairables
       integer::I
             
             do I=1, size(Vect)
                call Deallocate_DevVec_DF(Vect(I))
             end do
             return
   end subroutine Deallocate_DevVec_DF_b
   !****************************************************************************
  
   !****************************************************************************
   subroutine Deallocate_DevVec_I(Vect)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
       implicit none
       !--- dummy variables
       class(DevVec_I)        ::Vect
       !--- Local vairables
       integer::ERR, CURDEV
             
             if(Vect%IDEV .lt. 0) return

             ERR = cudaGetDevice(CURDEV)
             ERR = cudaSetDevice(Vect%IDEV)
             if(allocated(Vect%Data)) then
                deallocate(Vect%Data, STAT=ERR)
                if(ERR .gt. 0) then
                   write(*,*) "Error:: fail deallocating memory in MSM_MutltiGPU for "//Vect%Tag(1:len_trim(Vect%Tag))
                   stop   
                end if   
             end if 
             Vect%IDEV = -1 
             ERR = cudaSetDevice(CURDEV)
             return
   end subroutine Deallocate_DevVec_I
   !****************************************************************************

   !****************************************************************************
   subroutine Deallocate_DevVec_I_b(Vect)
   !***  PURPOSE:  to deallocate memory for device data
   !     INPUT     
   !
          implicit none
          !--- dummy variables
          class(DevVec_I), dimension(:) :: Vect
          !--- Local vairables
          integer::I
                
                do I=1, size(Vect)
                   call Deallocate_DevVec_I(Vect(I))
                end do   
                return
   end subroutine Deallocate_DevVec_I_b
   !****************************************************************************
   
   !****************************************************************************
   subroutine Collect_SwapData()
   !***  PURPOSE:  to extract the data produced in the last 
   !               calling of Product    
   !
       implicit none
       !--- dummy variables
       real(KINDDF)::Res
       !--- Local vairables
       integer::I, ERR, CURDEV
             
             ERR  = cudaGetDevice(CURDEV) 
             do I = 1, m_NDEVICE
                ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)
             end do
             
             call SynchronizeDevices() 
             ERR  = cudaSetDevice(CURDEV)             
             return
   end subroutine Collect_SwapData
   !****************************************************************************

   !****************************************************************************
   subroutine Get_SumData(Res)
      !***  PURPOSE:  to extract the data produced in the last 
      !               calling of Product    
      !
          implicit none
          !--- dummy variables
          real(KINDDF)::Res
          !--- Local vairables
          integer::I
                
                call Collect_SwapData() 
                Res = 0.D0
                do I = 1, m_NDEVICE
                    Res = Res + sum(dm_Array2DSwap(I)%RESDF) 
                end do     
                return
      end subroutine Get_SumData
   !****************************************************************************
   
   !****************************************************************************
   subroutine Get_MaxData(Res)
   !***  PURPOSE:  to extract the data produced in the last 
   !               calling of Product    
   !
       implicit none
       !--- dummy variables
       real(KINDDF)::Res
       !--- Local vairables
       integer::I
             
             call Collect_SwapData()
             Res = -1.D108
             do I = 1, m_NDEVICE
                 Res = max(Res, maxval(dm_Array2DSwap(I)%RESDF) )
             end do     
             return
   end subroutine Get_MaxData
   !****************************************************************************

   !****************************************************************************
   subroutine Get_MinData(Res)
   !***  PURPOSE:  to extract the data produced in the last 
   !               calling of Product    
   !
       implicit none
       !--- dummy variables
       real(KINDDF)::Res
       !--- Local vairables
       integer::I
             
             call Collect_SwapData()             
             Res = 1.D108
             do I = 1, m_NDEVICE
                 Res = min(Res, minval(dm_Array2DSwap(I)%RESDF) )
             end do     
             return
   end subroutine Get_MinData
   !****************************************************************************   

   !****************************************************************************
   attributes(global) subroutine Dot_DevVec_DF_KERNEL0(X1, X2, V1, V2, NPB, BRES)
   !***  PURPOSE:   KERNEL    to calculate the dot between two array: sum(V1(X1:X2)*V2(X1:X2))
   !
   !$$   INPUT:     X1, X2,    the bounds of the array
   !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
   !$$              NPB,       the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,       the dot_product of arraies handled by the block 
   !
   implicit none
   !----   DUMMY Variables
          integer,     value                :: NPB, X1, X2
          real(KINDDF),device, dimension(:) :: V1, V2
          real(KINDDF),device, dimension(:) :: BRES
  
  !----   Local variables
          integer        :: IT, IB, IM, I, J, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM + IT + X1 -1
                 if(J.le.X2) then
                    S(IT) = S(IT) + V1(J)*V2(J) 
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = BRES(IB) + S(IT)
          return
  end subroutine Dot_DevVec_DF_KERNEL0
  !****************************************************************************

  !****************************************************************************
  attributes(global) subroutine Dot_DevVec_DF_KERNEL1(X1, X2, V1, V2, GID, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2)*V2(GID(X1:X2))
  !
  !$$   INPUT:     X1, X2,    the bounds of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !$$              GID,       the index of V2 
  !$$              NPB,       the number of cross-block operations needed
  !$$
  !$$   OUTPUT     RES,       the dot_product of arraies handled by the block 
  !
  implicit none
  !----   DUMMY Variables
          integer,     value                :: NPB, X1, X2
          integer,     device, dimension(:) :: GID
          real(KINDDF),device, dimension(:) :: V1, V2
          real(KINDDF),device, dimension(:) :: BRES
  
  !----   Local variables
          integer             :: IT, IB, I, J, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + X1- 1
                 if(J .le. X2) then
                    S(IT) = S(IT) + V1(J)*V2(GID(J))
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = BRES(IB) + S(IT)
        return
  end subroutine Dot_DevVec_DF_KERNEL1
  !**************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template0(IDEV, X1, X2, Y1, Y2, V1, V2, PRDCT)
  !***  PURPOSE:  to calculate the dot between two array:: sum(V1(X1:X2,Y1:Y2)*V2(X1:X2,Y1:Y2))
  !                    
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,     the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,      intent(in)            :: IDEV, X1, X2, Y1, Y2
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I, Y

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (X2 - X1 + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  !--- don't forget initalize dm_Array2DSwap
                  ERR  = cudaMemsetAsync(dm_Array2DSwap(I)%Data, 0.D0, mp_ArrayOp_Gridsize, 0)

                  do Y = Y1, Y2
                     call Dot_DevVec_DF_KERNEL0<<<blocks, threads>>>(X1, X2, V1(:,Y), V2(:,Y), NPB, dm_Array2DSwap(I)%Data)
                  end do   

                  !-- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)
                  if(present(PRDCT)) then
                     PRDCT= sum(dm_Array2DSwap(I)%Data)
                  end if   
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             return
   end subroutine Dot_DevMat_DF_template0
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template1(IDEV, X1, X2, Y1, Y2, V1, V2, GID, PRDCT)
  !***  PURPOSE:   KERNEL    to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,Y1:Y2)*V2(GID(X1:X2),Y1:Y2)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,      intent(in)            :: IDEV, X1, X2, Y1, Y2
       real(KINDDF), device, dimension(:,:):: V1, V2
       integer,      device, dimension(:)  :: GID
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I, Y 

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (X2 - X1 + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  !--- don't forget initalize dm_Array2DSwap
                  ERR  = cudaMemsetAsync(dm_Array2DSwap(I)%Data, 0.D0, mp_ArrayOp_Gridsize, 0)

                  do Y = Y1, Y2
                     call Dot_DevVec_DF_KERNEL1<<<blocks, threads>>>(X1, X2, V1(:,Y), V2(:,Y), GID(:), NPB, dm_Array2DSwap(I)%Data)
                  end do   
                  !--- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)
                  if(present(PRDCT)) then
                     PRDCT= sum(dm_Array2DSwap(I)%Data)
                  end if
                  exit
                end if
             end do
             ERR   = cudaSetDevice(CURDEV)
             return
   end subroutine Dot_DevMat_DF_template1
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template0_a(X1, X2, V1, V2, PRDCT)
  !***  PURPOSE:   to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,:)*V2(GID(X1:X2),:)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,      intent(in)  :: X1, X2
       class(DevMat_DF)          :: V1, V2
       real(KINDDF), optional    :: PRDCT
       !--- Device variables and variables to be used in GPU
                   
             if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0) return
             if(V1%IDEV .ne. V2%IDEV) then
                write(*,*) "Error:: dot-product of two vector on different devices in MSM_MutltiGPU"
                stop
             end if  

             if(present(PRDCT)) then
                call Dot_DevMat_DF_template0(V1%IDEV, X1, X2, 1, 3, V1%Data, V2%Data, PRDCT)
             else
                call Dot_DevMat_DF_template0(V1%IDEV, X1, X2, 1, 3, V1%Data, V2%Data)  
             end if  
             return
   end subroutine Dot_DevMat_DF_template0_a
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template1_a(X1, X2, V1, V2, GID, PRDCT)
  !***  PURPOSE:   to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,:)*V2(GID(X1:X2),:)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,      intent(in)  :: X1, X2
       class(DevMat_DF)          :: V1, V2
       class(DevVec_I)           :: GID
       real(KINDDF), optional    :: PRDCT
       !--- Device variables and variables to be used in GPU

             if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0 .or. GID%IDEV .lt. 0) return
             if(V1%IDEV .ne. V2%IDEV  .or. V1%IDEV .ne. GID%IDEV) then
               write(*,*) "Error:: dot-product of two vector on different devices in MSM_MutltiGPU"
               stop
             end if  
             
             if(present(PRDCT)) then
                call Dot_DevMat_DF_template1(V1%IDEV, X1, X2, 1, 3, V1%Data, V2%Data, GID%Data, PRDCT)
             else
                call Dot_DevMat_DF_template1(V1%IDEV, X1, X2, 1, 3, V1%Data, V2%Data, GID%Data)  
             end if  
             return
   end subroutine Dot_DevMat_DF_template1_a
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template0_b(X1, X2, V1, V2, PRDCT)
  !***  PURPOSE:   to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,:)*V2(GID(X1:X2),:)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,          dimension(:),intent(in) :: X1, X2
       class(DevMat_DF), dimension(:)            :: V1, V2
       real(KINDDF)                              :: PRDCT
       !--- Device variables and variables to be used in GPU
       integer::I
             
             do I = 1, m_NDEVICE
                call Dot_DevMat_DF_template0_a(X1(I), X2(I), V1(I), V2(I))  
             end do
             call Get_SumData(PRDCT)  
             return
   end subroutine Dot_DevMat_DF_template0_b
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template0_c(X2, V1, V2, PRDCT)
  !***  PURPOSE:   to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,:)*V2(GID(X1:X2),:)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,          dimension(:),intent(in) :: X2
       class(DevMat_DF), dimension(:)            :: V1, V2
       real(KINDDF)                              :: PRDCT
       !--- Device variables and variables to be used in GPU
       integer::I
             
             do I = 1, m_NDEVICE
                call Dot_DevMat_DF_template0_a(1, X2(I), V1(I), V2(I))  
             end do
             call Get_SumData(PRDCT)  
             return
   end subroutine Dot_DevMat_DF_template0_c
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template1_b(X1, X2, V1, V2, GID, PRDCT)
  !***  PURPOSE:   to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,:)*V2(GID(X1:X2),:)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,          dimension(:),intent(in) :: X1, X2
       class(DevMat_DF), dimension(:)            :: V1, V2
       class(DevVec_I),  dimension(:)            :: GID
       real(KINDDF)                              :: PRDCT
       !--- Device variables and variables to be used in GPU
       integer::I
             
             do I = 1, m_NDEVICE
                call Dot_DevMat_DF_template1_a(X1(I), X2(I), V1(I), V2(I), GID(I))  
             end do
             call Get_SumData(PRDCT)  
             return
   end subroutine Dot_DevMat_DF_template1_b
  !***************************************************************************

  !**************************************************************************
  subroutine Dot_DevMat_DF_template1_c(X2, V1, V2, GID, PRDCT)
  !***  PURPOSE:   to calculate the dot between two array wiht mismatched index:
  !                            sum(V1(X1:X2,:)*V2(GID(X1:X2),:)
  !
  !$$   INPUT:     X1, X2     the bound of the array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     PRDCT,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer,          dimension(:),intent(in) :: X2
       class(DevMat_DF), dimension(:)            :: V1, V2
       class(DevVec_I),  dimension(:)            :: GID
       real(KINDDF)                              :: PRDCT
       !--- Device variables and variables to be used in GPU
       integer::I
             
             do I = 1, m_NDEVICE
                call Dot_DevMat_DF_template1_a(1, X2(I), V1(I), V2(I), GID(I))  
             end do
             call Get_SumData(PRDCT)  
             return
   end subroutine Dot_DevMat_DF_template1_c
  !***************************************************************************

  !**************************************************************************
  subroutine Normal2d_template0(IDEV, X1, X2, V, NORM)
  !***  PURPOSE:   to calculate the normal of segement of an array
  !
  !$$   INPUT:     X1, X2,   the segement size of array
  !$$              V,        the arrary 
  !
  !$$   OUTPUT     NORM,     the NORM of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV, X1, X2
       real(KINDDF), device, dimension(:,:):: V
       real(KINDDF), optional              :: NORM

          if(present(NORM)) then
             call Dot_DevMat_DF_template0(IDEV, X1, X2, 1, 3, V, V, NORM)
          else
             call Dot_DevMat_DF_template0(IDEV, X1, X2, 1, 3, V, V)  
          end if   
          return
  end subroutine Normal2d_template0
  !**************************************************************************

  !**************************************************************************
  subroutine Normal2d_template0_a(X1, X2, V, NORM)
   !***  PURPOSE:   to calculate the normal of an array
   !
   !$$   INPUT:     X1, X2,   the segement size of array
   !$$              V,        the two arraies for them dot-pruduct to be performed
   !
   !$$   OUTPUT     NORM,     the Normal
   implicit none
       !----   DUMMY Variables
        integer                 :: X1, X2
        class(DevMat_DF)        :: V
        real(KINDDF), optional  :: NORM
 
           if(present(NORM)) then
              call Dot_DevMat_DF_template0_a(X1, X2, V, V, NORM)
           else
              call Dot_DevMat_DF_template0_a(X1, X2, V, V)  
           end if   
           return
   end subroutine Normal2d_template0_a
  !**************************************************************************

  !**************************************************************************
  subroutine Normal2d_template0_b(X1, X2, V, NORM)
  !***  PURPOSE:   to calculate the normal of an array
  !
   !$$   INPUT:     X1, X2,   the segement size of array
   !$$              V,        the two arraies for them dot-pruduct to be performed
   !
   !$$   OUTPUT     NORM,     the Normal
   implicit none
      !----   DUMMY Variables
       integer         , dimension(:) :: X1, X2
       class(DevMat_DF), dimension(:) :: V
       real(KINDDF)                   :: NORM

          call Dot_DevMat_DF_template0_b(X1, X2, V, V, NORM)
          return
  end subroutine Normal2d_template0_b
  !**************************************************************************

    !**************************************************************************
  subroutine Normal2d_template0_c(X2, V, NORM)
  !***  PURPOSE:   to calculate the normal of an array from 1:X2
  !
  !$$   INPUT:     X2,   the size of array
  !$$              V,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     NORM,  the normal
  implicit none
      !----   DUMMY Variables
       integer         , dimension(:) :: X2
       class(DevMat_DF), dimension(:) :: V
       real(KINDDF)                   :: NORM

          call Dot_DevMat_DF_template0_c( X2, V, V, NORM)
          return
  end subroutine Normal2d_template0_c
  !**************************************************************************
 
  !****************************************************************************
  attributes(global) subroutine Sep_DevMat_DF_KERNEL0(X1, X2, V1, V2 , BSX, BSY, BSZ, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot of difference of two arraies
  !                          sum((V1(X1:X2,:)-V2(X1:X2,:)*(V1(X1:X2,:)-V2(X1:X2,:))  
  !$$   INPUT:     X1, X2        the bounds of array
  !$$              V1, V2,       the two arraies for them dot-pruduct to be performed
  !$$              BSX, BSY, BSZ the bouindary of box
  !$$   OUTPUT     RES,          the dot_product of arraies handled by the block 
  !
  implicit none
  !----   DUMMY Variables
          integer,      value                 :: NPB, X1, X2
          real(KINDDF), value                 :: BSX, BSY, BSZ
          real(KINDDF), device, dimension(:,:):: V1, V2
          real(KINDDF), device, dimension(:)  :: BRES
  
  !----   Local variables
          integer             :: IT, IB, I, J, J1, J2, IM, OFFSET    
          real(KINDDF)        :: SEPX, SEPY, SEPZ, HBSX, HBSY, HBSZ
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB

              HBSX = 0.5D0*BSX
              HBSY = 0.5D0*BSY
              HBSZ = 0.5D0*BSZ
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + X1 -1
                 if(J.le.X2) then
                    SEPX  = V1(J,1) - V2(J,1)
                    SEPY  = V1(J,2) - V2(J,2)
                    SEPZ  = V1(J,3) - V2(J,3)
                    if(dabs(SEPX) .gt. HBSX) SEPX = SEPX - dsign(BSX, SEPX)
                    if(dabs(SEPY) .gt. HBSY) SEPY = SEPY - dsign(BSY, SEPY)
                    if(dabs(SEPZ) .gt. HBSZ) SEPZ = SEPZ - dsign(BSZ, SEPY)
                    S(IT) = S(IT) + SEPX*SEPX + SEPY*SEPZ + SEPY*SEPZ
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
        return
  end subroutine Sep_DevMat_DF_KERNEL0
  !****************************************************************************

  !****************************************************************************
  attributes(global) subroutine Sep_DevMat_DF_KERNEL1(X1, X2, V1, V2, GID, BSX, BSY, BSZ, NPB, BRES)
  !***  PURPOSE:   KERNEL    to calculate the dot of difference of two arraies
  !                          sum((V1(X1:X2,:)-V2(GID(X1:X2),:)*(V1(X1:X2,:)-V2(GID(X1:X2),:))  
  !
  !$$   INPUT:     X1, X2        the bounds of array
  !$$              V1, V2,       the two arraies for them dot-pruduct to be performed
  !$$              GID,          the index reference between V1 and V2
  !$$              BSX, BSY, BSZ the bouindary of box
  !
  !$$   OUTPUT     BRES,       the dot_product of arraies handled by the block 
  !
  implicit none
  !----   DUMMY Variables
          integer,      value                 :: NPB, X1, X2
          real(KINDDF), value                 :: BSX, BSY, BSZ
          real(KINDDF), device, dimension(:,:):: V1, V2
          integer,      device, dimension(:)  :: GID
          real(KINDDF), device, dimension(:)  :: BRES

  
  !----   Local variables
          integer             :: IT, IB, I, J, IM, OFFSET    
          real(KINDDF)        :: SEPX, SEPY, SEPZ, HBSX, HBSY, HBSZ
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = 0.D0
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + X1 - 1
                 if(J.le.X2) then
                    SEPX  = V1(J,1) - V2(GID(J),1)
                    SEPY  = V1(J,2) - V2(GID(J),2)
                    SEPZ  = V1(J,3) - V2(GID(J),3)
                    if(dabs(SEPX) .gt. HBSX) SEPX = SEPX - dsign(BSX, SEPX)
                    if(dabs(SEPY) .gt. HBSY) SEPY = SEPY - dsign(BSY, SEPY)
                    if(dabs(SEPZ) .gt. HBSZ) SEPZ = SEPZ - dsign(BSZ, SEPY)
                    S(IT) = S(IT) + SEPX*SEPX + SEPY*SEPZ + SEPY*SEPZ
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = S(IT) + S(IT+IM)
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
        return
  end subroutine Sep_DevMat_DF_KERNEL1
  !**************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template0(IDEV, X1, X2, V1, V2, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:  to calculate the dot of difference of two arraies
  !               sum((V1(X1:X2,:)-V2(X1:X2,:)*(V1(X1:X2,:)-V2(X1:X2,:))  
  !
  !     INPUT:   IDEV,       the ID of the device  
  !              X1, X2,    the bound of the subset
  !              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !     OUTPUT   PRDCT,     the results
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV, X1, X2
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: N, NPB, CURDEV, ERR, I 

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (X2-X1+1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Sep_DevMat_DF_KERNEL0<<<blocks, threads>>>(X1, X2, V1, V2, BSX, BSY, BSZ, NPB, dm_Array2DSwap(I)%Data)
                  !--- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)
                  if(present(PRDCT)) then
                     PRDCT = sum(dm_Array2DSwap(I)%Data)
                  end if   
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
         return
   end subroutine Sep_DevMat_DF_template0
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template1(IDEV, X1, X2, V1, V2, GID, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   KERNEL    to calculate the dot of difference of two arraies
  !                          sum((V1(X1:X2,:)-V2(GID(X1:X2),:)*(V1(X1:X2,:)-V2(GID(X1:X2),:))  
  !
  !$$   INPUT:     N,         the size of array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     RES,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV, X1, X2
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       integer,      device, dimension(:)  :: GID
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (X2 - X1 +1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Sep_DevMat_DF_KERNEL1<<<blocks, threads>>>(X1, X2, V1, V2, GID, BSX, BSY, BSZ, NPB, dm_Array2DSwap(I)%Data)
                  !--- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)
                  if(present(PRDCT)) then
                     PRDCT = sum(dm_Array2DSwap(I)%Data)
                  end if   
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             return
   end subroutine Sep_DevMat_DF_template1
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template0_1(IDEV, V1, V2, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   to calculate the dot of difference of two arraies
  !                   sum((V1(:,:)-V2(:,:)*(V1(:,:)-V2(:,:))  
  !
  !$$   INPUT:     N,         the size of array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     RES,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU
             if(present(PRDCT)) then
                call Sep_DevMat_DF_template0(IDEV, 1, size(V1,dim=1), V1, V2, BSX, BSY, BSZ, PRDCT)
             else  
               call Sep_DevMat_DF_template0(IDEV, 1, size(V1,dim=1), V1, V2, BSX, BSY, BSZ)  
             end if  
             return
   end subroutine Sep_DevMat_DF_template0_1
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template1_1(IDEV, V1, V2, GID, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:   to calculate the dot of difference of two arraies
  !                sum((V1(:,:)-V2(GID(:),:)*(V1(:,:)-V2(GID(:),:))  
  !
  !$$   INPUT:     N,         the size of array
  !$$              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !$$   OUTPUT     RES,       the dot_product of arraies handled by the block 
  implicit none
      !----   DUMMY Variables
       integer                             :: IDEV
       real(KINDDF), device, dimension(:,:):: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       integer,      device, dimension(:)  :: GID
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU
             
             if(present(PRDCT)) then
                call Sep_DevMat_DF_template1(IDEV, 1, size(V1,dim=1), V1, V2, GID, BSX, BSY, BSZ, PRDCT)
             else   
               call Sep_DevMat_DF_template1(IDEV, 1, size(V1,dim=1), V1, V2, GID, BSX, BSY, BSZ)
             end if  
             return
   end subroutine Sep_DevMat_DF_template1_1
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template0_a(X1, X2, V1, V2, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:  to calculate the dot of difference of two arraies
  !               sum((V1(X1:X2,:)-V2(X1:X2,:)*(V1(X1:X2,:)-V2(X1:X2,:))  
  !
  !     INPUT:   
  !              X1, X2,    the bound of the subset
  !              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !     OUTPUT   PRDCT,     the results
  implicit none
      !----   DUMMY Variables
       integer                             :: X1, X2
       class(DevMat_DF)                    :: V1, V2
       real(KINDDF)                        :: BSX, BSY, BSZ
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU

             if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0) return

             if(V1%IDEV .ne. V2%IDEV) then
               write(*,*) "Error:: dot-product of two vector on different devices in MSM_MutltiGPU"
               stop
             end if  

             if(present(PRDCT)) then
                call Sep_DevMat_DF_template0(V1%IDEV, X1, X2, V1%Data, V2%Data, BSX, BSY, BSZ, PRDCT) 
             else
                call Sep_DevMat_DF_template0(V1%IDEV, X1, X2, V1%Data, V2%Data, BSX, BSY, BSZ)   
             end if    
         return
   end subroutine Sep_DevMat_DF_template0_a
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template1_a(X1, X2, V1, V2, GID, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:  to calculate the dot of difference of two arraies
  !               sum((V1(X1:X2,:)-V2(X1:X2,:)*(V1(X1:X2,:)-V2(X1:X2,:))  
  !
  !     INPUT:   
  !              X1, X2,    the bound of the subset
  !              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !     OUTPUT   PRDCT,     the results
  implicit none
      !----   DUMMY Variables
       integer                             :: X1, X2
       class(DevMat_DF)                    :: V1, V2
       class(DevVec_I)                     :: GID
       real(KINDDF)                        :: BSX, BSY, BSZ
       real(KINDDF), optional              :: PRDCT
       !--- Device variables and variables to be used in GPU

             if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0 .or. GID%IDEV .lt. 0) return
             if(V1%IDEV .ne. V2%IDEV  .or. V1%IDEV .ne. GID%IDEV) then
               write(*,*) "Error:: dot-product of two vector on different devices in MSM_MutltiGPU"
               stop
             end if  

             if(present(PRDCT)) then
                call Sep_DevMat_DF_template1(V1%IDEV, X1, X2, V1%Data, V2%Data, GID%Data, BSX, BSY, BSZ, PRDCT) 
             else
                call Sep_DevMat_DF_template1(V1%IDEV, X1, X2, V1%Data, V2%Data, GID%Data, BSX, BSY, BSZ)   
             end if    
         return
   end subroutine Sep_DevMat_DF_template1_a
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template0_b(X1, X2, V1, V2, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:  to calculate the dot of difference of two arraies
  !               sum((V1(X1:X2,:)-V2(X1:X2,:)*(V1(X1:X2,:)-V2(X1:X2,:))  
  !
  !     INPUT:   
  !              X1, X2,    the bound of the subset
  !              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !     OUTPUT   PRDCT,     the results
  implicit none
      !----   DUMMY Variables
       integer,         dimension(:)  :: X1, X2
       class(DevMat_DF),dimension(:)  :: V1, V2
       real(KINDDF)                   :: BSX, BSY, BSZ
       real(KINDDF)                   :: PRDCT
       !--- Device variables and variables to be used in GPU
       integer::I 

             do I = 1, m_NDEVICE
                call Sep_DevMat_DF_template0_a(X1(I), X2(I), V1(I), V2(I), BSX, BSY, BSZ)
             end do
             call Get_SumData(PRDCT)     
         return
   end subroutine Sep_DevMat_DF_template0_b
  !***************************************************************************

  !**************************************************************************
  subroutine Sep_DevMat_DF_template1_b(X1, X2, V1, V2, GID, BSX, BSY, BSZ, PRDCT)
  !***  PURPOSE:  to calculate the dot of difference of two arraies
  !               sum((V1(X1:X2,:)-V2(X1:X2,:)*(V1(X1:X2,:)-V2(X1:X2,:))  
  !
  !     INPUT:   
  !              X1, X2,    the bound of the subset
  !              V1, V2,    the two arraies for them dot-pruduct to be performed
  !
  !     OUTPUT   PRDCT,     the results
  implicit none
      !----   DUMMY Variables
       integer,         dimension(:)  :: X1, X2
       class(DevMat_DF),dimension(:)  :: V1, V2
       class(DevVec_I), dimension(:)  :: GID
       real(KINDDF)                   :: BSX, BSY, BSZ
       real(KINDDF)                   :: PRDCT
       !--- Device variables and variables to be used in GPU
       integer::I 

             do I = 1, m_NDEVICE
                call Sep_DevMat_DF_template1_a(X1(I), X2(I), V1(I), V2(I), GID(I), BSX, BSY, BSZ)
             end do
             call Get_SumData(PRDCT)     
         return
   end subroutine Sep_DevMat_DF_template1_b
  !***************************************************************************

  !***************************************************************************
   attributes(global) subroutine Scalar_DevVec_DF_KERNEL0(N, Scalar, V, Res)
   !***  PURPOSE:   KERNEL    to mulply scalar to a vector:
   !                          Res = Scalar*V(1:N)
   !
   !$$   INPUT:     X1, X2    the bounds of the array
   !$$              Scalar,   the scalar 
   !$$              NPB,      the number of cross-block operations needed
   !$$              V,        the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !
   implicit none
   !----   DUMMY Variables
           integer,     value                 :: N
           real(KINDDF),value                 :: Scalar           
           real(KINDDF),device, dimension(:)  :: V, Res
   
   !----   Local variables
           integer :: IT, IB, J
   
               IT   = (threadidx%y-1)*blockdim%x + threadidx%x
               IB   = (blockidx%y-1) * griddim%x +  blockidx%x
               J     = (IB-1)*mp_BlockSize + IT
               if(J .le. N) then
                  Res(J) = Scalar*V(J)
               end if
         return
   end subroutine Scalar_DevVec_DF_KERNEL0
   !**************************************************************************

   !**************************************************************************
   subroutine Scalar_DevMat_DF_template0(IDEV, X1, X2, Y1, Y2, Scalar, V, Res)
   !***  PURPOSE:   to mulply scalar to a vector: Res = Scalar*V(X1:X2,Y1:Y2)
   !
   !$$   INPUT:     X1, X2, Y1,Y2,  the bounds of the array
   !$$              Scalar,         the scalar 
   !$$              V,              the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,      intent(in)            :: IDEV, X1, X2, Y1, Y2
        real(KINDDF), intent(in)            :: Scalar
        real(KINDDF), device, dimension(:,:):: V, Res
        !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, Y, N
 
              ERR     = cudaGetDevice(CURDEV)
              ERR     = cudaSetDevice(IDEV)
              N       = X2 - X1 + 1
              blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
              threads = dim3(mp_BlockSize, 1, 1)
              
              do Y = Y1, Y2
                   call Scalar_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, Scalar, V(X1:,Y), Res(X1:,Y))
              end do     
              ERR   = cudaSetDevice(CURDEV)
 
              return
    end subroutine Scalar_DevMat_DF_template0
   !*************************************************************************** 
    
   !**************************************************************************
    subroutine Scalar_DevMat_DF_template0_a(Dimx, Dimy, Scalar, V, Res)
      !***  PURPOSE:   to mulply scalar to a vector: Res = Scalar*V(X1:X2,Y1:Y2)
      !
      !$$   INPUT:     Dimx, Dimy,  the bounds of the array
      !$$              Scalar,      the scalar 
      !$$              V,           the array
      !$$
      !$$   OUTPUT     Res,         the scaled vector
      !   
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: Dimx(2), Dimy(2)
           real(KINDDF), intent(in)  :: Scalar
           class(DevMat_DF)          :: V, Res
           !--- Device variables and variables to be used in GPU
    
                 if(V%IDEV .lt. 0) return

                 if(Res%IDEV .lt. 0) then
                    call Allocate_DevMat_DF_0 (V%IDEV, Res, (/size(V%Data,dim=1), size(V%Data,dim=2)/))
                 end if  

                 if(V%IDEV .ne. Res%IDEV) then
                    write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                    stop
                 end if  
   
                 call Scalar_DevMat_DF_template0(V%IDEV, Dimx(1), Dimx(2), Dimy(1), Dimy(2), Scalar, V%Data, Res%Data)
                 return
       end subroutine Scalar_DevMat_DF_template0_a
   !***************************************************************************

   !**************************************************************************
   subroutine Scalar_DevMat_DF_template0_b(X1, X2, Scalar, V, Res)
   !***  PURPOSE:   to mulply scalar to a vector: Res = Scalar*V(X1:X2,1:3)
   !
   !$$   INPUT:     X1, X2,    the bounds of the array
   !$$              Scalar,    the scalar 
   !$$              V,         the array
   !$$
   !$$   OUTPUT     V,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        real(KINDDF),                   intent(in) :: Scalar
        class(DevMat_DF), dimension(:)             :: V, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Scalar_DevMat_DF_template0_a((/X1(I), X2(I)/), (/1,3/), Scalar, V(I), Res(I))
              end do   
              return
    end subroutine Scalar_DevMat_DF_template0_b
   !***************************************************************************

   !**************************************************************************
   subroutine Scalar_DevMat_DF_template0_c(X2, Scalar, V, Res)
   !***  PURPOSE:   to mulply scalar to a vector: Res = Scalar*V(1:X2,1:3)
   !
   !$$   INPUT:     X2,        the bounds of the array
   !$$              Scalar,    the scalar 
   !$$              V,         the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        real(KINDDF),                   intent(in) :: Scalar
        class(DevMat_DF), dimension(:)             :: V, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Scalar_DevMat_DF_template0_a((/1,X2(I)/), (/1,3/), Scalar, V(I), Res(I))
              end do   
              return
    end subroutine Scalar_DevMat_DF_template0_c
   !***************************************************************************    
 
   !**************************************************************************
    subroutine Scalar_DevMat_DF_template1_a(Dimx, Dimy, Scalar, V)
         implicit none
         !----   DUMMY Variables
          integer,      intent(in)  :: Dimx(2), Dimy(2)
          real(KINDDF), intent(in)  :: Scalar
          class(DevMat_DF)          :: V
   
                call Scalar_DevMat_DF_template0_a(Dimx, Dimy, Scalar, V, V)
                return
    end subroutine Scalar_DevMat_DF_template1_a
   !***************************************************************************

   !**************************************************************************
   subroutine Scalar_DevMat_DF_template1_b(X1, X2, Scalar, V)
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        real(KINDDF),                   intent(in) :: Scalar
        class(DevMat_DF), dimension(:)             :: V

             call Scalar_DevMat_DF_template0_b(X1, X2, Scalar, V, V)
             return
    end subroutine Scalar_DevMat_DF_template1_b
   !***************************************************************************

   !**************************************************************************
   subroutine Scalar_DevMat_DF_template1_c(X2, Scalar, V)
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        real(KINDDF),                   intent(in) :: Scalar
        class(DevMat_DF), dimension(:)             :: V
        
              call Scalar_DevMat_DF_template0_c(X2, Scalar, V, V)
              return
    end subroutine Scalar_DevMat_DF_template1_c
   !***************************************************************************    

   
   !***************************************************************************
    attributes(global) subroutine MakeCopy_DevVec_DF_KERNEL0(N, Sor, Tgt)
    !***  PURPOSE:   KERNEL    to copy array SOR to TGT
    !                          TGT(1:N) = SOR(1:N)
    !
    !$$   INPUT:     N,         the bounds of the array
    !$$              NPB,       the number of cross-block operations needed
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !
    implicit none
    !----   DUMMY Variables
            integer,     value                 :: N
            real(KINDDF),device, dimension(:)  :: Sor, Tgt
    
    !----   Local variables
            integer  :: IT, IB, J  

                IT   = (threadidx%y-1)*blockdim%x + threadidx%x
                IB   = (blockidx%y-1) * griddim%x +  blockidx%x
                J     = (IB-1)*mp_BlockSize + IT
                if(J .le. N) then
                   Tgt(J) = Sor(J)
                end if                
          return
    end subroutine MakeCopy_DevVec_DF_KERNEL0
    !**************************************************************************
 
   !***************************************************************************
    attributes(global) subroutine MakeCopy_DevVec_DF_KERNEL1(N, Sor, Gid, Tgt)
    !***  PURPOSE:   KERNEL    to copy array SOR to TGT
    !                          TGT(X1:X2) = SOR(GID(X1:X2),:)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              NPB,       the number of cross-block operations needed
    !$$              SOR,       the source array
    !$$              GID,       the index of SOR 
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !
    implicit none
    !----   DUMMY Variables
            integer,     value                 :: N
            integer,     device, dimension(:)  :: GID
            real(KINDDF),device, dimension(:)  :: SOR, TGT
    
    !----   Local variables
            integer             :: IT, IB, I, J, IM, OFFSET    
    
                IT   = (threadidx%y-1)*blockdim%x + threadidx%x
                IB   = (blockidx%y-1) * griddim%x +  blockidx%x
                J     = (IB-1)*mp_BlockSize + IT
                if(J .le. N) then
                   Tgt(J) = Sor(Gid(J))
                end if                
          return
    end subroutine MakeCopy_DevVec_DF_KERNEL1
   !**************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template0(IDEV, X1, X2, Y1, Y2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(X1:X2,Y1:Y2)
    !
    !$$   INPUT:     X1:X2,Y1:Y2   the bounds of the array
    !$$              SOR,          the source array
    !$$              
    !$$   OUTPUT     TGT,          the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer,      intent(in)            :: IDEV, X1, X2, Y1,Y2
         real(KINDDF), device, dimension(:,:):: Sor, Tgt
           !--- Device variables and variables to be used in GPU
             type(dim3) :: blocks
             type(dim3) :: threads
             integer    :: CURDEV, ERR, N, Y
    
                ERR     = cudaGetDevice(CURDEV)
                ERR     = cudaSetDevice(IDEV)
                N       = X2 - X1 + 1
                blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
                threads = dim3(mp_BlockSize, 1, 1)
  
                do Y = Y1, Y2
                    call   MakeCopy_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, Sor(X1:,Y), Tgt(X1:,Y))
                end do    
                ERR   = cudaSetDevice(CURDEV)
         return
    end subroutine MakeCopy_DevMat_DF_template0
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template1(IDEV, X1, X2, Y1, Y2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2, Y1:Y2) = SOR(X1:X2,Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              NPB,       the number of cross-block operations needed
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)            :: IDEV, X1, X2, Y1, Y2
           integer,      device, dimension(:)  :: Gid
           real(KINDDF), device, dimension(:,:):: Sor, Tgt
           !--- Device variables and variables to be used in GPU
             type(dim3) :: blocks
             type(dim3) :: threads
             integer    :: CURDEV, ERR, N, Y
    
                ERR     = cudaGetDevice(CURDEV)
                ERR     = cudaSetDevice(IDEV)
                N       = X2 - X1 + 1
                blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
                threads = dim3(mp_BlockSize, 1, 1)
  
                do Y = Y1, Y2
                    call   MakeCopy_DevVec_DF_KERNEL1<<<blocks, threads>>>(N, Sor(X1:,Y), Gid(X1:), Tgt(X1:,Y))
                end do    
                ERR   = cudaSetDevice(CURDEV)    
          return
    end subroutine MakeCopy_DevMat_DF_template1
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template0_a(X1, X2, Y1, Y2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2, Y1:Y2) = SOR(X1:X2,Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer, intent(in)           :: X1, X2, Y1, Y2
         class(DevMat_DF)              :: Sor, Tgt
           !--- Device variables and variables to be used in GPU

                if(Sor%IDEV .lt. 0) return 
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevMat_DF_0(Sor%IDEV, Tgt, (/size(Sor%Data,dim=1),size(Sor%Data,dim=2)/) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if    

                call MakeCopy_DevMat_DF_template0(Sor%IDEV, X1, X2, Y1,Y2, Sor%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevMat_DF_template0_a
   !***************************************************************************    

   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template0_b(X1, X2, Y1, Y2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2,:)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,         dimension(:), intent(in) :: X1, X2, Y1, Y2
           class(DevMat_DF),dimension(:)             :: Sor, Tgt
           !--- Device variables and variables to be used in GPU
           integer::I 
                
                do I =1, m_NDEVICE
                   call MakeCopy_DevMat_DF_template0_a(X1(I), X2(I), Y1(I), Y2(I), Sor(I), Tgt(I))
                end do   
         return
    end subroutine MakeCopy_DevMat_DF_template0_b
   !*************************************************************************    

   !*************************************************************************
    subroutine MakeCopy_DevMat_DF_template0_c(N, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2,:)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
          integer,          dimension(:), intent(in) :: N
          class(DevMat_DF), dimension(:)             :: Sor, Tgt
          !--- Device variables and variables to be used in GPU
          integer, dimension(m_MXDEVICE)::X10, X20, Y10, Y20 
                  
                X10 = 1
                X20 = X10 + N -1
                Y10 = 1
                Y20 = 3
                call MakeCopy_DevMat_DF_template0_b(X10, X20, Y10, Y20, Sor, Tgt)
                     
           return
      end subroutine MakeCopy_DevMat_DF_template0_c
   !*************************************************************************** 
  
   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template1_a(X1, X2, Y1, Y2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(GID(X1:X2),Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,           intent(in) :: X1, X2, Y1, Y2
           class(DevMat_DF)              :: Sor, Tgt
           class(DevVec_I)               :: GID
           !--- Device variables and variables to be used in GPU

                if(Sor%IDEV .lt. 0) return 
                if(GID%IDEV .lt. 0) return
                if(Sor%IDEV .ne. GID%IDEV) then
                  write(*,*) "Error:: Sor and Gid vector allocated on different devices in MSM_MutltiGPU"
                  stop
                end if  
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevMat_DF_0(Sor%IDEV, Tgt, (/size(Sor%Data,dim=1),size(Sor%Data,dim=2)/) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if      
                  
                call MakeCopy_DevMat_DF_template1(Sor%IDEV, X1, X2, Y1, Y2, Sor%Data, Gid%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevMat_DF_template1_a
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template1_b(X1, X2, Y1, Y2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(GID(X1:X2),Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,         dimension(:), intent(in) :: X1, X2, Y1, Y2
           class(DevMat_DF),dimension(:)             :: Sor, Tgt
           class(DevVec_I), dimension(:)             :: Gid
           !--- local variables 
           integer::I 
                
                 do I =1, m_NDEVICE
                     call MakeCopy_DevMat_DF_template1_a(X1(I), X2(I), Y1(I), Y2(I), Sor(I), Gid(I), Tgt(I))
                 end do   
         return
    end subroutine MakeCopy_DevMat_DF_template1_b
   !***************************************************************************   
    
   !**************************************************************************
    subroutine MakeCopy_DevMat_DF_template1_c(N, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(1:N, 1:3) = SOR(GID(1:N),1:3)
    !
    !$$   INPUT:     N,         the size to be copy
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,          dimension(:), intent(in) :: N
           class(DevMat_DF),  dimension(:)  :: Sor, Tgt
           class(DevVec_I),   dimension(:)  :: Gid
           !--- local variables
            integer, dimension(m_MXDEVICE)::X10, X20, Y10, Y20 
                   
                 X10 = 1
                 X20 = X10 + N -1
                 Y10 = 1
                 Y20 = 3           
                
                 call MakeCopy_DevMat_DF_template1_b(X10, X20, Y10, Y20, Sor, Gid, Tgt)
         return
    end subroutine MakeCopy_DevMat_DF_template1_c
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template0(IDEV, X1, X2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1:X2,        the bounds of the array
    !$$              SOR,          the source array
    !$$              
    !$$   OUTPUT     TGT,          the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer,      intent(in)          :: IDEV, X1, X2
         real(KINDDF), device, dimension(:):: Sor, Tgt
           !--- Device variables and variables to be used in GPU
             type(dim3) :: blocks
             type(dim3) :: threads
             integer    :: CURDEV, ERR, N
    
                ERR     = cudaGetDevice(CURDEV)
                ERR     = cudaSetDevice(IDEV)
                N       = X2 - X1 + 1
                blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
                threads = dim3(mp_BlockSize, 1, 1)
  
                  call MakeCopy_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, Sor(X1:), Tgt(X1:))
                ERR   = cudaSetDevice(CURDEV)
         return
    end subroutine MakeCopy_DevVec_DF_template0
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template1(IDEV, X1, X2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)            :: IDEV, X1, X2
           integer,      device, dimension(:)  :: Gid
           real(KINDDF), device, dimension(:)  :: Sor, Tgt
           !--- local variables
             type(dim3) :: blocks
             type(dim3) :: threads
             integer    :: CURDEV, ERR, N
    
                ERR     = cudaGetDevice(CURDEV)
                ERR     = cudaSetDevice(IDEV)
                N       = X2 - X1 + 1
                blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
                threads = dim3(mp_BlockSize, 1, 1)
  
                  call   MakeCopy_DevVec_DF_KERNEL1<<<blocks, threads>>>(N, Sor(X1:), Gid(X1:), Tgt(X1:))
                ERR   = cudaSetDevice(CURDEV)    
          return
    end subroutine MakeCopy_DevVec_DF_template1
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template0_a(X1, X2,Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer, intent(in)           :: X1, X2
         class(DevVec_DF)              :: Sor, Tgt
           !--- local varibales

                if(Sor%IDEV .lt. 0) return 
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevVec_DF_0(Sor%IDEV, Tgt, size(Sor%Data) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if    

                call MakeCopy_DevVec_DF_template0(Sor%IDEV, X1, X2, Sor%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevVec_DF_template0_a
   !***************************************************************************    

   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template0_b(X1, X2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,         dimension(:), intent(in) :: X1, X2
           class(DevVec_DF),dimension(:)             :: Sor, Tgt
           !--- Device variables and variables to be used in GPU
           integer::I 
                
                do I =1, m_NDEVICE
                   call MakeCopy_DevVec_DF_template0_a(X1(I), X2(I), Sor(I), Tgt(I))
                end do   
         return
    end subroutine MakeCopy_DevVec_DF_template0_b
   !*************************************************************************    

   !*************************************************************************
    subroutine MakeCopy_DevVec_DF_template0_c(N, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
          integer,          dimension(:), intent(in) :: N
          class(DevVec_DF), dimension(:)             :: Sor, Tgt
          !--- local varaibles
          integer, dimension(m_MXDEVICE)::X10, X20
                  
                X10 = 1
                X20 = X10 + N -1
                call MakeCopy_DevVec_DF_template0_b(X10, X20, Sor, Tgt)
                     
           return
      end subroutine MakeCopy_DevVec_DF_template0_c
   !*************************************************************************** 
  
   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template1_a(X1, X2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(GID(X1:X2))
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,           intent(in) :: X1, X2
           class(DevVec_DF)              :: Sor, Tgt
           class(DevVec_I)               :: GID
           !--- local varaibles

                if(Sor%IDEV .lt. 0) return 
                if(GID%IDEV .lt. 0) return
                if(Sor%IDEV .ne. GID%IDEV) then
                  write(*,*) "Error:: Sor and Gid vector allocated on different devices in MSM_MutltiGPU"
                  stop
                end if  
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevVec_DF_0(Sor%IDEV, Tgt, size(Sor%Data) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if      
                  
                call MakeCopy_DevVec_DF_template1(Sor%IDEV, X1, X2, Sor%Data, Gid%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevVec_DF_template1_a
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template1_b(X1, X2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(GID(X1:X2),Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,         dimension(:), intent(in) :: X1, X2
           class(DevVec_DF),dimension(:)             :: Sor, Tgt
           class(DevVec_I), dimension(:)             :: Gid
           !--- local variables 
           integer::I 
                
                 do I =1, m_NDEVICE
                     call MakeCopy_DevVec_DF_template1_a(X1(I), X2(I), Sor(I), Gid(I), Tgt(I))
                 end do   
         return
    end subroutine MakeCopy_DevVec_DF_template1_b
   !***************************************************************************   
    
   !**************************************************************************
    subroutine MakeCopy_DevVec_DF_template1_c(N, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(1:N) = SOR(GID(1:N))
    !
    !$$   INPUT:     N,         the size to be copy
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,           dimension(:), intent(in) :: N
           class(DevVec_DF),  dimension(:)  :: Sor, Tgt
           class(DevVec_I),   dimension(:)  :: Gid
           !--- local variables
            integer, dimension(m_MXDEVICE)::X10, X20
                   
                 X10 = 1
                 X20 = X10 + N -1
                
                 call MakeCopy_DevVec_DF_template1_b(X10, X20,  Sor, Gid, Tgt)
         return
    end subroutine MakeCopy_DevVec_DF_template1_c
   !***************************************************************************    

   !***************************************************************************
    attributes(global) subroutine MakeCopy_DevVec_I_KERNEL0(N, SOR, TGT)
    !***  PURPOSE:   KERNEL    to copy array SOR to TGT
    !                          TGT(1:N) = SOR(1:N)
    !
    !$$   INPUT:     N,        the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !
    implicit none
    !----   DUMMY Variables
            integer, value                 :: N
            integer ,device, dimension(:)  :: Sor, Tgt
    
    !----   Local variables
            integer  :: IT, IB, J  

                IT   = (threadidx%y-1)*blockdim%x + threadidx%x
                IB   = (blockidx%y-1) * griddim%x +  blockidx%x
                J     = (IB-1)*mp_BlockSize + IT
                if(J .le. N) then
                   Tgt(J) = Sor(J)
                end if                

          return
    end subroutine MakeCopy_DevVec_I_KERNEL0
    !**************************************************************************
 
   !***************************************************************************
    attributes(global) subroutine MakeCopy_DevVec_I_KERNEL1(N, SOR, GID, TGT)
    !***  PURPOSE:   KERNEL    to copy array SOR to TGT
    !                          TGT(1:N) = SOR(GID(1:N))
    !
    !$$   INPUT:     N,         the bounds of the array
    !$$              SOR,       the source array
    !$$              GID,       the index of SOR 
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !
    implicit none
    !----   DUMMY Variables
            integer, value                 :: N
            integer ,device, dimension(:)  :: Sor, Gid, Tgt
    
    !----   Local variables
            integer  :: IT, IB, J  

                IT   = (threadidx%y-1)*blockdim%x + threadidx%x
                IB   = (blockidx%y-1) * griddim%x +  blockidx%x
                J     = (IB-1)*mp_BlockSize + IT
                if(J .le. N) then
                   Tgt(J) = Sor(Gid(J))
                end if                
          return
    end subroutine MakeCopy_DevVec_I_KERNEL1
   !**************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template0(IDEV, X1, X2, Y1, Y2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(X1:X2, Y1:Y2)
    !
    !     INPUT:     X1, X2,    the bounds of the array
    !                SOR,       the source array
    !     OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer, intent(in)            :: IDEV, X1, X2, Y1, Y2
         integer, device, dimension(:,:):: Sor, Tgt

          !--- local vraibles
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, N, Y
    
             ERR     = cudaGetDevice(CURDEV)
             ERR     = cudaSetDevice(IDEV)
             N       = X2 - X1 + 1
             blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
             threads = dim3(mp_BlockSize, 1, 1)
  
             do Y = Y1, Y2
                call   MakeCopy_DevVec_I_KERNEL0<<<blocks, threads>>>(N, Sor(X1:,Y), Tgt(X1:,Y))
             end do    
             ERR   = cudaSetDevice(CURDEV)                 
    
         return
    end subroutine MakeCopy_DevMat_I_template0
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template0_a(X1, X2, Y1, Y2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(X1:X2, Y1:Y2)
    !
    !$$   INPUT:     X1, X2, Y1, Y2, the bounds of the array
    !$$              SOR,            the source array
    !$$   OUTPUT     TGT,           the target array
    !   
         implicit none
          !----   DUMMY Variables
          integer, intent(in)   :: X1, X2, Y1, Y2
          class(DevMat_I)       :: Sor, Tgt
           !--- local variables

                if(Sor%IDEV .lt. 0) return 
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevMat_I_0(Sor%IDEV, Tgt, (/size(Sor%Data,dim=1),size(Sor%Data,dim=2)/) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if    
                call MakeCopy_DevMat_I_template0(Sor%IDEV, X1, X2, Y1, Y2, Sor%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevMat_I_template0_a
   !***************************************************************************    

   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template0_b(X1, X2, Y1, Y2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(X1:X2,Y1:Y2)
    !
    !$$   INPUT:     X1, X2, Y1, Y2   the bounds of the array
    !$$              SOR,             the source array
    !$$              
    !$$   OUTPUT     TGT,             the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer,         dimension(:), intent(in) :: X1, X2, Y1, Y2
         class(DevMat_I), dimension(:)             :: Sor, Tgt
           !--- local variables
           integer::I 
                
               do I =1, m_NDEVICE
                  call MakeCopy_DevMat_I_template0_a(X1(I), X2(I), Y1(I), Y2(I), Sor(I), Tgt(I))
               end do   
         return
    end subroutine MakeCopy_DevMat_I_template0_b
   !***************************************************************************    

   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template0_c(NX, NY, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(1:NX,1:NY) = SOR(1:X2,1:NY)
    !
    !$$   INPUT:     N, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
          integer,          intent(in)      :: NX, NY
          class(DevMat_I),  dimension(:)    :: Sor, Tgt

          !--- Device variables and variables to be used in GPU
          integer, dimension(m_MXDEVICE)::X10, X20, Y10, Y20 
                  
               X10 = 1
               X20 = X10 + NX - 1
               Y10 = 1
               Y20 = Y10 + NY - 1
               call MakeCopy_DevMat_I_template0_b(X10, X20, Y10, Y20, Sor, Tgt)
           return
    end subroutine MakeCopy_DevMat_I_template0_c
   !*************************************************************************** 

   !**************************************************************************
    subroutine MakeCopy_DevVec_I_template0(IDEV, X1, X2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !     INPUT:     X1, X2,    the bounds of the array
    !                SOR,       the source array
    !     OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer, intent(in)          :: IDEV, X1, X2
         integer, device, dimension(:):: Sor, Tgt

          !--- local vraibles
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, N, Y
    
             ERR     = cudaGetDevice(CURDEV)
             ERR     = cudaSetDevice(IDEV)
             N       = X2 - X1 + 1
             blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
             threads = dim3(mp_BlockSize, 1, 1)
  
             call   MakeCopy_DevVec_I_KERNEL0<<<blocks, threads>>>(N, Sor(X1:), Tgt(X1:))
             ERR   = cudaSetDevice(CURDEV)                 
    
         return
    end subroutine MakeCopy_DevVec_I_template0
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevVec_I_template0_a(X1, X2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1, X2, the bounds of the array
    !$$              SOR,    the source array
    !$$   OUTPUT     TGT,    the target array
    !   
         implicit none
          !----   DUMMY Variables
          integer, intent(in)   :: X1, X2
          class(DevVec_I)       :: Sor, Tgt
           !--- local variables

                if(Sor%IDEV .lt. 0) return 
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevVec_I_0(Sor%IDEV, Tgt, size(Sor%Data) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if    
                call MakeCopy_DevVec_I_template0(Sor%IDEV, X1, X2, Sor%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevVec_I_template0_a
   !***************************************************************************    

   !**************************************************************************
    subroutine MakeCopy_DevVec_I_template0_b(X1, X2, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2) = SOR(X1:X2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
         integer,         dimension(:), intent(in) :: X1, X2
         class(DevVec_I), dimension(:)             :: Sor, Tgt
           !--- local variables
           integer::I 
                
               do I =1, m_NDEVICE
                  call MakeCopy_DevVec_I_template0_a(X1(I), X2(I), Sor(I), Tgt(I))
               end do   
         return
    end subroutine MakeCopy_DevVec_I_template0_b
   !***************************************************************************    

   !**************************************************************************
    subroutine MakeCopy_DevVec_I_template0_c(N, Sor, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(1:N) = SOR(1:N)
    !
    !$$   INPUT:     N,        the sizeof the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
          integer,          intent(in)      :: N
          class(DevVec_I),  dimension(:)    :: Sor, Tgt

          !--- Device variables and variables to be used in GPU
          integer, dimension(m_MXDEVICE)::X10, X20
                  
               X10 = 1
               X20 = X10 + N - 1
               call MakeCopy_DevVec_I_template0_b(X10, X20, Sor, Tgt)
           return
    end subroutine MakeCopy_DevVec_I_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template1(IDEV, X1, X2, Y1, Y2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(X1:X2, Y1:Y2)
    !
    !     INPUT:     X1, X2, Y1, Y2, the bounds of the array
    !                SOR,            the source array
    !     OUTPUT     TGT,            the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,  intent(in)            :: IDEV, X1, X2, Y1, Y2
           integer,  device, dimension(:)  :: Gid
           integer,  device, dimension(:,:):: Sor, Tgt
           !--- local vraibles
             type(dim3) :: blocks
             type(dim3) :: threads
             integer    :: CURDEV, ERR, N, Y
    
                ERR     = cudaGetDevice(CURDEV)
                ERR     = cudaSetDevice(IDEV)
                N       = X2 - X1 + 1
                blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
                threads = dim3(mp_BlockSize, 1, 1)
     
                do Y = Y1, Y2
                   call   MakeCopy_DevVec_I_KERNEL1<<<blocks, threads>>>(N, Sor(X1:,Y), Gid(X1:), Tgt(X1:,Y))
                end do    
                ERR   = cudaSetDevice(CURDEV)                                
    
          return
    end subroutine MakeCopy_DevMat_I_template1
   !***************************************************************************
  
   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template1_a(X1, X2, Y1, Y2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2, Y1:Y2) = SOR(GID(X1:X2),Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,           intent(in) :: X1, X2, Y1, Y2
           class(DevMat_I)               :: Sor, Tgt
           class(DevVec_I)               :: GID
           !--- Device variables and variables to be used in GPU

                if(Sor%IDEV .lt. 0) return 
                if(GID%IDEV .lt. 0) return
                if(Sor%IDEV .ne. GID%IDEV) then
                  write(*,*) "Error:: Sor and Gid vector allocated on different devices in MSM_MutltiGPU"
                  stop
                end if  
                   
                
                if(Tgt%IDEV .lt. 0) then
                   call Allocate_DevMat_I_0(Sor%IDEV, Tgt, (/size(Sor%Data,dim=1),size(Sor%Data,dim=2)/) )
                else  
                  if(Sor%IDEV .ne. Tgt%IDEV) then
                     write(*,*) "Error:: the sor and tgt data allocated on different devices in MSM_MutltiGPU"
                     stop
                   end if  
                end if      
                call MakeCopy_DevMat_I_template1(Sor%IDEV, X1, X2, Y1, Y2, Sor%Data, Gid%Data, Tgt%Data)
         return
    end subroutine MakeCopy_DevMat_I_template1_a
   !***************************************************************************

   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template1_b(X1, X2, Y1, Y2, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(X1:X2,Y1:Y2) = SOR(GID(X1:X2),Y1:Y2)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,         dimension(:), intent(in) :: X1, X2, Y1, Y2
           class(DevMat_I), dimension(:)             :: Sor, Tgt
           class(DevVec_I), dimension(:)             :: Gid
           !--- Device variables and variables to be used in GPU
           integer::I 
                
                do I =1, m_NDEVICE
                   call MakeCopy_DevMat_I_template1_a(X1(I), X2(I), Y1(I), Y2(I), Sor(I), Gid(I), Tgt(I))
                end do   
         return
    end subroutine MakeCopy_DevMat_I_template1_b
   !***************************************************************************   
    
   !**************************************************************************
    subroutine MakeCopy_DevMat_I_template1_c(NX, NY, Sor, Gid, Tgt)
    !***  PURPOSE:   to copy array SOR to TGT
    !                TGT(1:NX,1:NY) = SOR(GID(1:NX),1:NY)
    !
    !$$   INPUT:     X1, X2,    the bounds of the array
    !$$              SOR,       the source array
    !$$              
    !$$   OUTPUT     TGT,       the target array
    !   
         implicit none
          !----   DUMMY Variables
           integer,           intent(in)    :: NX, NY
           class(DevMat_I),   dimension(:)  :: Sor, Tgt
           class(DevVec_I),   dimension(:)  :: Gid
           !--- local variables
           integer, dimension(m_MXDEVICE)::X10, X20, Y10, Y20 
                  
               X10 = 1
               X20 = X10 + NX - 1
               Y10 = 1
               Y20 = Y10 + NY - 1
               call MakeCopy_DevMat_I_template1_b(X10, X20, Y10, Y20, Sor, Gid, Tgt)

         return
    end subroutine MakeCopy_DevMat_I_template1_c
   !***************************************************************************
   
   !**************************************************************************
    subroutine Set_DevMat_DF_template0( Mat, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             real(KINDDF)     :: Val
             class(DevMat_DF) :: Mat
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR

                if(Mat%IDEV .lt. 0)  return

                ERR   = cudaGetDevice(CURDEV)
                ERR   = cudaSetDevice(Mat%IDEV)
                !Mat%Data = Val
                ERR   = cudaMemsetAsync(Mat%Data(1,1), Val, size(Mat%Data), 0) 
                ERR   = cudaSetDevice(CURDEV)
           return
     end subroutine Set_DevMat_DF_template0
   !***************************************************************************

   !**************************************************************************
    subroutine Set_DevMat_I_template0( Mat, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             integer         :: Val
             class(DevMat_I) :: Mat
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR

                if(Mat%IDEV .lt. 0)  return

                ERR   = cudaGetDevice(CURDEV)
                ERR   = cudaSetDevice(Mat%IDEV)
                !Mat%Data = Val
                ERR   = cudaMemsetAsync(Mat%Data(1,1), Val, size(Mat%Data), 0) 
                ERR   = cudaSetDevice(CURDEV)
           return
     end subroutine Set_DevMat_I_template0
   !***************************************************************************

   !**************************************************************************
    subroutine Set_DevVec_DF_template0( Vect, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             real(KINDDF)     :: Val
             class(DevVec_DF) :: Vect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR

                if(Vect%IDEV .lt. 0)  return

                ERR   = cudaGetDevice(CURDEV)
                ERR   = cudaSetDevice(Vect%IDEV)
                !Vect%Data = Val
                ERR   = cudaMemsetAsync(Vect%Data(1), Val, size(Vect%Data), 0) 
                ERR   = cudaSetDevice(CURDEV)
           return
     end subroutine Set_DevVec_DF_template0
   !***************************************************************************

   !**************************************************************************
    subroutine Set_DevVec_I_template0( Vect, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             integer         :: Val
             class(DevVec_I) :: Vect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR

                if(Vect%IDEV .lt. 0)  return

                ERR   = cudaGetDevice(CURDEV)
                ERR   = cudaSetDevice(Vect%IDEV)
                !Vect%Data = Val
                ERR   = cudaMemsetAsync(Vect%Data(1), Val, size(Vect%Data), 0) 
                ERR   = cudaSetDevice(CURDEV)
           return
     end subroutine Set_DevVec_I_template0
   !***************************************************************************     

   !**************************************************************************
    subroutine Set_DevMat_DF_template0_b( Mat, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             real(KINDDF)                  :: Val
             class(DevMat_DF),dimension(:) :: Mat
             !--- Device variables and variables to be used in GPU
             integer :: I
 
                do I=1, m_NDEVICE
                   call Set_DevMat_DF_template0( Mat(I), Val)
                end do
           return
     end subroutine Set_DevMat_DF_template0_b
   !***************************************************************************

   !**************************************************************************
    subroutine Set_DevMat_I_template0_b( Mat, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: Val
             class(DevMat_I),dimension(:) :: Mat
             !--- Device variables and variables to be used in GPU
             integer :: I
 
                do I=1, m_NDEVICE
                   call Set_DevMat_I_template0( Mat(I), Val)
                end do
           return
     end subroutine Set_DevMat_I_template0_b
   !***************************************************************************

   !**************************************************************************
    subroutine Set_DevVec_DF_template0_b( Vect, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             real(KINDDF)                  :: Val
             class(DevVec_DF),dimension(:) :: Vect
             !--- Device variables and variables to be used in GPU
             integer :: I
 
                do I=1, m_NDEVICE
                   call Set_DevVec_DF_template0( Vect(I), Val)
                end do
           return
     end subroutine Set_DevVec_DF_template0_b
   !***************************************************************************

   !**************************************************************************
    subroutine Set_DevVec_I_template0_b( Vect, Val)
    !***  PURPOSE:   to set value for the Array
    !                Mat = Value
    !
    !$$   INPUT:     Value, the value to set
    !$$              
    !$$   OUTPUT     Mat,   the Mat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: Val
             class(DevVec_I),dimension(:) :: Vect
             !--- Device variables and variables to be used in GPU
             integer :: I
 
                do I=1, m_NDEVICE
                   call Set_DevVec_I_template0( Vect(I), Val)
                end do
           return
     end subroutine Set_DevVec_I_template0_b
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevMat_DF_template0_a(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: hStart(2), dStart(2), TheDim(2)
             real(KINDDF), dimension(:,:) :: hMat
             class(DevMat_DF)             :: dMat
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dMat%IDEV .lt. 0)  return

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dMat%IDEV)
                do I=0, TheDim(2)-1
                   ERR  =  cudaMemcpyAsync(hMat(hStart(1),hStart(2)+I), dMat%Data(dStart(1),dStart(2)+I), TheDim(1))
                end do   
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyout_DevMat_DF_template0_a
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevMat_DF_template0_b(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: hStart, dStart, TheDim(2)
             real(KINDDF), dimension(:,:) :: hMat
             class(DevMat_DF)             :: dMat
             !--- Device variables and variables to be used in GPU
             call Copyout_DevMat_DF_template0_a(hMat, (/hStart,1/), dMat, (/dStart,1/), TheDim)
           return
     end subroutine Copyout_DevMat_DF_template0_b
   !***************************************************************************
     
   !**************************************************************************
    subroutine Copyout_DevMat_DF_template0_c(hMat, dMat, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: NC
             real(KINDDF), dimension(:,:) :: hMat
             class(DevMat_DF)             :: dMat
             !--- Device variables and variables to be used in GPU

                call Copyout_DevMat_DF_template0_a(hMat, (/1,1/), dMat, (/1,1/), (/NC, 1/) ) 
           return
     end subroutine Copyout_DevMat_DF_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevMat_DF_template0_d(hMat, hStart, dMat, dStart, NC, Dim2)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer,          dimension(:)   :: hStart, dStart, NC
             integer,          optional       :: Dim2
             real(KINDDF),     dimension(:,:) :: hMat
             class(DevMat_DF), dimension(:)   :: dMat
             !--- Device variables and variables to be used in GPU
             integer::I 
                if(present(Dim2)) then
                   do I=1, m_NDEVICE
                      call Copyout_DevMat_DF_template0_b(hMat, hStart(I), dMat(I), dStart(I), (/NC(I),Dim2/))
                   end do 
                else
                  do I=1, m_NDEVICE
                     call Copyout_DevMat_DF_template0_b(hMat, hStart(I), dMat(I), dStart(I), (/NC(I),3/))
                  end do 
               end if 
           return
     end subroutine Copyout_DevMat_DF_template0_d
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevMat_I_template0_a(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                 :: hStart(2), dStart(2), TheDim(2)
             integer, dimension(:,:) :: hMat
             class(DevMat_I)         :: dMat
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dMat%IDEV .lt. 0)  return

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dMat%IDEV)
                do I=0, TheDim(2)-1
                   ERR  =  cudaMemcpyAsync(hMat(hStart(1),hStart(2)+I), dMat%Data(dStart(1),dStart(2)+I), TheDim(1))
                end do   
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyout_DevMat_I_template0_a
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevMat_I_template0_b(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                 :: hStart, dStart, TheDim(2)
             integer, dimension(:,:) :: hMat
             class(DevMat_I)         :: dMat
             !--- Device variables and variables to be used in GPU
             call Copyout_DevMat_I_template0_a(hMat, (/hStart,1/), dMat, (/dStart,1/), TheDim)
           return
     end subroutine Copyout_DevMat_I_template0_b
   !***************************************************************************
     
   !**************************************************************************
    subroutine Copyout_DevMat_I_template0_c(hMat, dMat, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                 :: NC
             integer, dimension(:,:) :: hMat
             class(DevMat_I)        :: dMat
             !--- Device variables and variables to be used in GPU

                call Copyout_DevMat_I_template0_a(hMat, (/1,1/), dMat, (/1,1/), (/NC, 1/) ) 
           return
     end subroutine Copyout_DevMat_I_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevVec_DF_template0_a( hVect, hStart, dVect, dStart, W)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                    :: hStart, dStart, W
             real(KINDDF), dimension(:) :: hVect
             class(DevVec_DF)           :: dVect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dVect%IDEV .lt. 0)  return

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dVect%IDEV)
                ERR  = cudaMemcpyAsync(hVect(hStart), dVect%Data(dStart), W)
                ERR  = cudaSetDevice(CURDEV)
           return
    end subroutine Copyout_DevVec_DF_template0_a
   !***************************************************************************
     
   !**************************************************************************
    subroutine Copyout_DevVec_DF_template0_c( hVect, dVect, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                    :: NC
             real(KINDDF), dimension(:) :: hVect
             class(DevVec_DF)           :: dVect
             !--- Device variables and variables to be used in GPU

                call Copyout_DevVec_DF_template0_a( hVect, 1, dVect, 1, NC) 
           return
     end subroutine Copyout_DevVec_DF_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevVec_I_template0_a( hVect, hStart, dVect, dStart, W)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                        :: hStart, dStart, W
             integer,          dimension(:) :: hVect
             class(DevVec_I)                :: dVect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dVect%IDEV .lt. 0)  return

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dVect%IDEV)
                ERR  = cudaMemcpyAsync(hVect(hStart), dVect%Data(dStart), W)
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyout_DevVec_I_template0_a
   !***************************************************************************

   !**************************************************************************
    subroutine Copyout_DevVec_I_template0_c( hVect, dVect, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                        :: NC
             integer,          dimension(:) :: hVect
             class(DevVec_I)                :: dVect
             !--- 
               call Copyout_DevVec_I_template0_a( hVect, 1, dVect, 1, NC)

           return
     end subroutine Copyout_DevVec_I_template0_c
   !***************************************************************************   
     
   !**************************************************************************
    subroutine Copyout_DevVec_I_template0_d( hVect, hStart, dVect, dStart, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
           !----   DUMMY Variables
           integer,          dimension(:) :: hStart, dStart, NC
           integer,          dimension(:) :: hVect
           class(DevVec_I),  dimension(:) :: dVect

           !--- local variable
           integer::I
            
               do I=1, m_NDEVICE
                  call Copyout_DevVec_I_template0_a( hVect, hStart(I), dVect(I), dStart(I), NC(I))
               end do 
           return
    end subroutine Copyout_DevVec_I_template0_d
   !***************************************************************************
     
   !**************************************************************************
    subroutine Copyin_DevMat_DF_template0_a(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: hStart(2), dStart(2), TheDim(2)
             real(KINDDF), dimension(:,:) :: hMat
             class(DevMat_DF)             :: dMat
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dMat%IDEV .lt. 0)  then
                   call DevAllocate(dMat, (/size(hMat,dim=1),size(hMat,dim=2)/) )
                end if   

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dMat%IDEV)
                do  I= 0, TheDim(2)-1
                   ERR  =  cudaMemcpyAsync(dMat%Data(dStart(1),dStart(2)+I), hMat(hStart(1),hStart(2)+I), TheDim(1))
                end do   
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyin_DevMat_DF_template0_a
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevMat_DF_template0_b(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: hStart, dStart, TheDim(2)
             real(KINDDF), dimension(:,:) :: hMat
             class(DevMat_DF)             :: dMat
             !--- Device variables and variables to be used in GPU

                call Copyin_DevMat_DF_template0_a(hMat, (/hStart,1/), dMat, (/dStart,1/), TheDim)
           return
     end subroutine Copyin_DevMat_DF_template0_b
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevMat_DF_template0_c(hMat, dMat, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                      :: NC
             real(KINDDF), dimension(:,:) :: hMat
             class(DevMat_DF)             :: dMat
             !--- Device variables and variables to be used in GPU

                call Copyin_DevMat_DF_template0_a(hMat, (/1,1/), dMat, (/1,1/), (/NC,1/) )
           return
     end subroutine Copyin_DevMat_DF_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevMat_DF_template0_d(hMat, dMat, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                         :: NC
             real(KINDDF),     dimension(:,:):: hMat
             class(DevMat_DF), dimension(:)  :: dMat
             !--- loacal variables
             integer::I   
             
                do I=1, m_NDEVICE 
                   call Copyin_DevMat_DF_template0_c(hMat, dMat(I), NC )
                end do   
           return
     end subroutine Copyin_DevMat_DF_template0_d
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevMat_I_template0_a(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                  :: hStart(2), dStart(2), TheDim(2)
             integer,  dimension(:,:) :: hMat
             class(DevMat_I)          :: dMat
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dMat%IDEV .lt. 0)  then
                   call DevAllocate(dMat, (/size(hMat,dim=1),size(hMat,dim=2)/) )
                end if   

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dMat%IDEV)
                do I=0, TheDim(2)-1
                   ERR  =  cudaMemcpyAsync(dMat%Data(dStart(1),dStart(2)+I), hMat(hStart(1),hStart(2)+I), TheDim(1))
                end do   
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyin_DevMat_I_template0_a
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevMat_I_template0_b(hMat, hStart, dMat, dStart, TheDim)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                 :: hStart, dStart, TheDim(2)
             integer, dimension(:,:) :: hMat
             class(DevMat_I)         :: dMat
             !--- Device variables and variables to be used in GPU

                call Copyin_DevMat_I_template0_a(hMat, (/hStart,1/), dMat, (/dStart,1/), TheDim)
           return
     end subroutine Copyin_DevMat_I_template0_b
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevMat_I_template0_c(hMat, dMat, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                 :: NC
             integer, dimension(:,:) :: hMat
             class(DevMat_I)        :: dMat
             !--- Device variables and variables to be used in GPU

                call Copyin_DevMat_I_template0_a(hMat, (/1,1/), dMat, (/1,1/), (/NC,1/) )
           return
     end subroutine Copyin_DevMat_I_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevVec_DF_template0_a( hVect, hStart, dVect, dStart, W)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                    :: hStart, dStart, W
             real(KINDDF), dimension(:) :: hVect
             class(DevVec_DF)           :: dVect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dVect%IDEV .lt. 0)  then
                  call DevAllocate(dVect, size(hVect))
                end if

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dVect%IDEV)
                ERR  = cudaMemcpyAsync(dVect%Data(dStart), hVect(hStart), W)
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyin_DevVec_DF_template0_a
   !***************************************************************************
     
   !**************************************************************************
    subroutine Copyin_DevVec_DF_template0_c( hVect, dVect,  NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                    :: NC
             real(KINDDF), dimension(:) :: hVect
             class(DevVec_DF)           :: dVect
             !--- Device variables and variables to be used in GPU
                 
                call Copyin_DevVec_DF_template0_a( hVect, 1, dVect, 1, NC) 
           return
     end subroutine Copyin_DevVec_DF_template0_c
   !***************************************************************************   

   !**************************************************************************
    subroutine Copyin_DevVec_DF_template0_d(hVec, dVec, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                         :: NC
             real(KINDDF),     dimension(:)  :: hVec
             class(DevVec_DF), dimension(:)  :: dVec
             !--- loacal variables
             integer::I   
             
                do I=1, m_NDEVICE 
                   call Copyin_DevVec_DF_template0_c(hVec, dVec(I), NC )
                end do   
           return
     end subroutine Copyin_DevVec_DF_template0_d
   !***************************************************************************     

   !**************************************************************************
    subroutine Copyin_DevVec_I_template0_a( hVect, hStart, dVect, dStart, W)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                        :: hStart, dStart, W
             integer,          dimension(:) :: hVect
             class(DevVec_I)                :: dVect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                if(dVect%IDEV .lt. 0)  then
                   call DevAllocate(dVect, size(hVect))
                end if

                ERR  = cudaGetDevice(CURDEV)
                ERR  = cudaSetDevice(dVect%IDEV)
                ERR  = cudaMemcpyAsync(dVect%Data(dStart), hVect(hStart), W)
                ERR  = cudaSetDevice(CURDEV)
           return
     end subroutine Copyin_DevVec_I_template0_a
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevVec_I_template0_c( hVect, dVect, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                        :: NC
             integer,          dimension(:) :: hVect
             class(DevVec_I)                :: dVect
             !--- Device variables and variables to be used in GPU
             integer :: CURDEV, ERR, I

                call Copyin_DevVec_I_template0_a( hVect, 1, dVect, 1, NC) 
           return
     end subroutine Copyin_DevVec_I_template0_c
   !***************************************************************************

   !**************************************************************************
    subroutine Copyin_DevVec_I_template0_d( hVec, hStart, dVec, dStart, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer,          dimension(:) :: hStart, dStart, NC
             integer,          dimension(:) :: hVec
             class(DevVec_I),  dimension(:) :: dVec
             !--- Device variables and variables to be used in GPU
             integer :: I

             do I=1, m_NDEVICE 
               call Copyin_DevVec_I_template0_a(hVec, hStart(I), dVec(I), dStart(I), NC(I) )
             end do   
       return
     end subroutine Copyin_DevVec_I_template0_d
   !***************************************************************************
     
   !**************************************************************************
    subroutine Copyin_DevVec_I_template0_e(hVec, dVec, NC)
    !***  PURPOSE:   to copy out from device to host
    !
    !$$   INPUT:     hMat, dMat
    !$$              
    !$$   OUTPUT     hMat
    !   
           implicit none
            !----   DUMMY Variables
             integer                       :: NC
             integer,         dimension(:) :: hVec
             class(DevVec_I), dimension(:) :: dVec
             !--- loacal variables
             integer::I   
             
                do I=1, m_NDEVICE 
                   call Copyin_DevVec_I_template0_c(hVec, dVec(I), NC )
                end do   
           return
     end subroutine Copyin_DevVec_I_template0_e
   !***************************************************************************     


     
   !****************************************************************************
   attributes(global) subroutine Max_DevMat_DF_KERNEL0(DimX1, DimX2, DimY1, DimY2, Mat, NPB, BRES)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     DimX, DimY,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     BRES,        the max of arraies handled by the block 
   !
   implicit none
   !----   DUMMY Variables
          integer,     value                 :: NPB, DimX1,DimX2,DimY1,DimY2
          real(KINDDF),device, dimension(:,:):: Mat
          real(KINDDF),device, dimension(:)  :: BRES
  
  !----   Local variables
          integer        :: IT, IB, I, J, K, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = -1.D108
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + DimX1 -1
                 if(J .le. DimX2) then
                     do K=DimY1, DimY2
                       if(S(IT) .lt. Mat(J,K) ) then
                          S(IT) = Mat(J,K)
                       end if
                     end do      
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = max(S(IT), S(IT+IM))
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
          return
  end subroutine Max_DevMat_DF_KERNEL0
  !****************************************************************************

  !**************************************************************************
  subroutine Max_DevMat_DF_template0(IDEV, Fromx, Tox, Fromy, Toy, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Fromy, Toy
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
   implicit none
      !----   DUMMY Variables
       integer,      intent(in)            :: IDEV, Fromx, Tox, Fromy, Toy
       real(KINDDF), device, dimension(:,:):: Mat
       real(KINDDF), optional              :: Res
       !--- local variables        
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (Tox - Fromx + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call Max_DevMat_DF_KERNEL0<<<blocks, threads>>>(Fromx, Tox, Fromy, Toy, Mat, NPB, dm_Array2DSwap(I)%Data)
                  !--- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)

                  if(present(Res)) then
                     Res= Maxval(dm_Array2DSwap(I)%Data)
                  end if   
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             return
   end subroutine Max_DevMat_DF_template0
  !***************************************************************************

  !***************************************************************************
   subroutine Max_DevMat_DF_template0_a(Fromx, Tox, Fromy, Toy, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Fromy, Toy
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,        intent(in)  :: Fromx, Tox, Fromy, Toy
       class(DevMat_DF)            :: Mat
       real(KINDDF),   optional    :: Res
       !--- local variables

             if(Mat%IDEV .lt. 0) return

             if(present(Res)) then
                call Max_DevMat_DF_template0(Mat%IDEV, Fromx, Tox, Fromy, Toy, Mat%Data, Res)
             else
                call Max_DevMat_DF_template0(Mat%IDEV, Fromx, Tox, Fromy, Toy, Mat%Data) 
             end if  
             return
   end subroutine Max_DevMat_DF_template0_a
  !***************************************************************************   

  !***************************************************************************
   subroutine Max_DevMat_DF_template0_a13(Fromx, Tox, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     DimX, DimY,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,        intent(in)  :: Fromx, Tox
       class(DevMat_DF)            :: Mat
       real(KINDDF),   optional    :: Res
       !--- local variables

             if(present(Res)) then
                call Max_DevMat_DF_template0_a(Fromx, Tox, 1, 3, Mat, Res)
                
             else
               call Max_DevMat_DF_template0_a(Fromx, Tox, 1, 3, Mat)
             end if  
             return
   end subroutine Max_DevMat_DF_template0_a13
  !***************************************************************************   

  !***************************************************************************
   subroutine Max_DevMat_DF_template0_b(Fromx, ToX, Fromy, Toy, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  Fromy, Toy, the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
      integer,     intent(in) , dimension(:) :: Fromx, Tox
      integer,     intent(in)                :: Fromy, Toy
      class(DevMat_DF),         dimension(:) :: Mat
       real(KINDDF)                          :: Res
       !--- local variables
       integer::I 
                   
             do I = 1, m_NDEVICE
                call Max_DevMat_DF_template0_a(Fromx(I), Tox(I), Fromy, Toy, Mat(I))
             end do
             call Get_MaxData(Res)
             return
   end subroutine Max_DevMat_DF_template0_b
  !***************************************************************************  
  
  !***************************************************************************
   subroutine Max_DevMat_DF_template0_b13(Fromx, Tox, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
      integer,     intent(in) , dimension(:) :: Fromx, Tox
      class(DevMat_DF),         dimension(:) :: Mat
       real(KINDDF)                          :: Res
       !--- local variables
                   
             call Max_DevMat_DF_template0_b(Fromx, Tox, 1, 3, Mat, Res)
             return
   end subroutine Max_DevMat_DF_template0_b13
  !***************************************************************************   

  
  !***************************************************************************
   subroutine Max_DevMat_DF_template0_c13(N, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,     intent(in) , dimension(:) :: N
       class(DevMat_DF),         dimension(:) :: Mat
       real(KINDDF)                          :: Res
       !--- local variables
       integer::I 
                   
             do I = 1, m_NDEVICE
                call Max_DevMat_DF_template0_a(1, N(I), 1, 3, Mat(I))
             end do
             call Get_MaxData(Res)
             return
   end subroutine Max_DevMat_DF_template0_c13
  !***************************************************************************
   
   !****************************************************************************
   attributes(global) subroutine MaxAbs_DevMat_DF_KERNEL0(DimX1, DimX2, DimY1, DimY2, Mat, NPB, BRES)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     DimX, DimY,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     BRES,        the max of arraies handled by the block 
   !
   implicit none
   !----   DUMMY Variables
          integer,     value                 :: NPB, DimX1,DimX2,DimY1,DimY2
          real(KINDDF),device, dimension(:,:):: Mat
          real(KINDDF),device, dimension(:)  :: BRES
  
  !----   Local variables
          integer        :: IT, IB, I, J, K, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = -1.D108
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + DimX1 -1
                 if(J .le. DimX2) then
                     do K=DimY1, DimY2
                       if(S(IT) .lt. dabs(Mat(J,K)) ) then
                          S(IT) = dabs(Mat(J,K))
                       end if
                     end do      
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = max(S(IT), S(IT+IM))
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
          return
  end subroutine MaxAbs_DevMat_DF_KERNEL0
  !****************************************************************************

  !**************************************************************************
  subroutine MaxAbs_DevMat_DF_template0(IDEV, Fromx, Tox, Fromy, Toy, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Fromy, Toy
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
   implicit none
      !----   DUMMY Variables
       integer,      intent(in)            :: IDEV, Fromx, Tox, Fromy, Toy
       real(KINDDF), device, dimension(:,:):: Mat
       real(KINDDF), optional              :: Res
       !--- local variables        
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (Tox - Fromx + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call MaxAbs_DevMat_DF_KERNEL0<<<blocks, threads>>>(Fromx, Tox, Fromy, Toy, Mat, NPB, dm_Array2DSwap(I)%Data)
                  !--- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)

                  if(present(Res)) then
                     Res= Maxval(dm_Array2DSwap(I)%Data)
                  end if   
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             return
   end subroutine MaxAbs_DevMat_DF_template0
  !***************************************************************************

  !***************************************************************************
   subroutine MaxAbs_DevMat_DF_template0_a(Fromx, Tox, Fromy, Toy, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Fromy, Toy
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,        intent(in)  :: Fromx, Tox, Fromy, Toy
       class(DevMat_DF)            :: Mat
       real(KINDDF),   optional    :: Res
       !--- local variables

             if(Mat%IDEV .lt. 0) return

             if(present(Res)) then
                call MaxAbs_DevMat_DF_template0(Mat%IDEV, Fromx, Tox, Fromy, Toy, Mat%Data, Res)
             else
                call MaxAbs_DevMat_DF_template0(Mat%IDEV, Fromx, Tox, Fromy, Toy, Mat%Data) 
             end if  
             return
   end subroutine MaxAbs_DevMat_DF_template0_a
  !***************************************************************************   

  !***************************************************************************
   subroutine MaxAbs_DevMat_DF_template0_a13(Fromx, Tox, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     DimX, DimY,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,        intent(in)  :: Fromx, Tox
       class(DevMat_DF)            :: Mat
       real(KINDDF),   optional    :: Res
       !--- local variables

             if(present(Res)) then
                call MaxAbs_DevMat_DF_template0_a(Fromx, Tox, 1, 3, Mat, Res)
                
             else
               call MaxAbs_DevMat_DF_template0_a(Fromx, Tox, 1, 3, Mat)
             end if  
             return
   end subroutine MaxAbs_DevMat_DF_template0_a13
  !***************************************************************************   

  !***************************************************************************
   subroutine MaxAbs_DevMat_DF_template0_b(Fromx, ToX, Fromy, Toy, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  Fromy, Toy, the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
      integer,     intent(in) , dimension(:) :: Fromx, Tox
      integer,     intent(in)                :: Fromy, Toy
      class(DevMat_DF),         dimension(:) :: Mat
       real(KINDDF)                          :: Res
       !--- local variables
       integer::I 
                   
             do I = 1, m_NDEVICE
                call MaxAbs_DevMat_DF_template0_a(Fromx(I), Tox(I), Fromy, Toy, Mat(I))
             end do
             call Get_MaxData(Res)
             return
   end subroutine MaxAbs_DevMat_DF_template0_b
  !***************************************************************************  
  
  !***************************************************************************
   subroutine MaxAbs_DevMat_DF_template0_b13(Fromx, Tox, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
      integer,     intent(in) , dimension(:) :: Fromx, Tox
      class(DevMat_DF),         dimension(:) :: Mat
       real(KINDDF)                          :: Res
       !--- local variables
                   
             call MaxAbs_DevMat_DF_template0_b(Fromx, Tox, 1, 3, Mat, Res)
             return
   end subroutine MaxAbs_DevMat_DF_template0_b13
  !***************************************************************************   

  
  !***************************************************************************
   subroutine MaxAbs_DevMat_DF_template0_c13(N, Mat, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,     intent(in) , dimension(:) :: N
       class(DevMat_DF),         dimension(:) :: Mat
       real(KINDDF)                           :: Res
       !--- local variables
       integer::I 
                   
             do I = 1, m_NDEVICE
                call MaxAbs_DevMat_DF_template0_a(1, N(I), 1, 3, Mat(I))
             end do
             call Get_MaxData(Res)
             return
   end subroutine MaxAbs_DevMat_DF_template0_c13
  !***************************************************************************
   
   !****************************************************************************
   attributes(global) subroutine MaxAbs_DevVec_DF_KERNEL0(DimX1, DimX2, V, NPB, BRES)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     DimX,        the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     BRES,        the max of arraies handled by the block 
   !
   implicit none
   !----   DUMMY Variables
          integer,     value                :: NPB, DimX1,DimX2
          real(KINDDF),device, dimension(:) :: V
          real(KINDDF),device, dimension(:) :: BRES
  
  !----   Local variables
          integer        :: IT, IB, I, J, IM, OFFSET    
          real(KINDDF), shared:: S(mp_ArrayOp_Blocksize)
  
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1) * griddim%x +  blockidx%x
              IM  =  mp_ArrayOp_Blocksize
              OFFSET = (IB-1)*IM*NPB
              
              S(IT) = -1.D108
              do I = 1, NPB
                 J = OFFSET + (I-1)*IM+IT + DimX1 -1
                 if(J .le. DimX2) then
                    if(S(IT) .lt. dabs(V(J)) ) then
                       S(IT) = dabs(V(J))
                    end if
                 end if
              end do
              call syncthreads()

              do I=1, mp_ArrayOp_2Power
                 IM = IM/2
                 call syncthreads()
                 if(IT .le. IM) then
                    S(IT) = max(S(IT), S(IT+IM))
                 end if
              end do
              call syncthreads()
              if(IT .eq. 1) BRES(IB) = S(IT)
          return
  end subroutine MaxAbs_DevVec_DF_KERNEL0
  !****************************************************************************

  !**************************************************************************
  subroutine MaxAbs_DevVec_DF_template0(IDEV, Fromx, Tox, V, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Fromy, Toy
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
   implicit none
      !----   DUMMY Variables
       integer,      intent(in)           :: IDEV, Fromx, Tox
       real(KINDDF), device, dimension(:) :: V
       real(KINDDF), optional             :: Res
       !--- local variables        
         type(dim3) :: blocks
         type(dim3) :: threads
         integer    :: NPB, CURDEV, ERR, I

             ERR     = cudaGetDevice(CURDEV)

             blocks  = dim3(mp_ArrayOp_Gridsize, 1, 1)
             threads = dim3(mp_ArrayOp_Blocksize, 1, 1)
             NPB     =  (Tox - Fromx + 1) /(mp_ArrayOp_Gridsize*mp_ArrayOp_Blocksize) + 1

             do I = 1, m_NDEVICE
                if(IDEV .eq. dm_Array2DSwap(I)%IDEV) then 
                  ERR  = cudaSetDevice(dm_Array2DSwap(I)%IDEV)
                  call MaxAbs_DevVec_DF_KERNEL0<<<blocks, threads>>>(Fromx, Tox, V, NPB, dm_Array2DSwap(I)%Data)
                  !--- ERR  = cudaMemcpyAsync(dm_Array2DSwap(I)%RESDF, dm_Array2DSwap(I)%Data, mp_ArrayOp_Gridsize)

                  if(present(Res)) then
                     Res= Maxval(dm_Array2DSwap(I)%Data)
                  end if   
                  exit
                end if
             end do

             ERR   = cudaSetDevice(CURDEV)
             return
   end subroutine MaxAbs_DevVec_DF_template0
  !***************************************************************************

  !***************************************************************************
   subroutine MaxAbs_DevVec_DF_template0_a(Fromx, Tox, V, Res)
   !***  PURPOSE:   KERNEL    to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$              NPB,         the number of cross-block operations needed
   !$$
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
       integer,        intent(in)  :: Fromx, Tox
       class(DevVec_DF)            :: V
       real(KINDDF),   optional    :: Res
       !--- local variables

             if(V%IDEV .lt. 0) return

             if(present(Res)) then
                call MaxAbs_DevVec_DF_template0(V%IDEV, Fromx, Tox, V%Data, Res)
             else
                call MaxAbs_DevVec_DF_template0(V%IDEV, Fromx, Tox, V%Data) 
             end if  
             return
   end subroutine MaxAbs_DevVec_DF_template0_a
  !***************************************************************************

  !***************************************************************************
   subroutine MaxAbs_DevVec_DF_template0_b(Fromx, ToX, V, Res)
   !***  PURPOSE:   to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     Fromx, Tox,  the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
      integer,     intent(in) , dimension(:) :: Fromx, Tox
      class(DevVec_DF),         dimension(:) :: V
       real(KINDDF)                          :: Res
       !--- local variables
       integer::I 
                   
             do I = 1, m_NDEVICE
                call MaxAbs_DevVec_DF_template0_a(Fromx(I), Tox(I), V(I))
             end do
             call Get_MaxData(Res)
             return
   end subroutine MaxAbs_DevVec_DF_template0_b
  !***************************************************************************

  !***************************************************************************
   subroutine MaxAbs_DevVec_DF_template0_c(ToX, V, Res)
   !***  PURPOSE:   to extract the maximum value of in the array segement 
   !
   !$$   INPUT:     1, Tox,      the segemnt dimension of the array
   !$$              Mat,         the array from which the max value value to be extracted
   !$$   OUTPUT     RES,         the max value obtained by device IDEV
   !$$                           this is optional, could be obtained later by call Get_MaxData 
      implicit none
      !----   DUMMY Variables
      integer,     intent(in) , dimension(:) :: Tox
      class(DevVec_DF),         dimension(:) :: V
       real(KINDDF)                          :: Res
       !--- local variables
       integer::I 
                   
             do I = 1, m_NDEVICE
                call MaxAbs_DevVec_DF_template0_a(1, Tox(I), V(I))
             end do
             call Get_MaxData(Res)
             return
   end subroutine MaxAbs_DevVec_DF_template0_c
  !***************************************************************************

  !***************************************************************************
   attributes(global) subroutine AddScalar_DevVec_I_KERNEL0(N, V, Scalar, Res)
   !***  PURPOSE:   KERNEL    to add a scalar to all the elements of a vector
   !                          Res(1:N) = V1(1:N) + V2(1:N)
   !
   !$$   INPUT:     N         
   !$$              V1, V2    
   !$$   OUTPUT     Res,      
   !
   implicit none
   !----   DUMMY Variables
           integer,        value       :: N, Scalar
           integer,device, dimension(:):: V, Res
   
   !----   Local variables
           integer :: IT, IB, J, K
   
               IT   = (threadidx%y-1)*blockdim%x + threadidx%x
               IB   = (blockidx%y-1) * griddim%x +  blockidx%x
               J     = (IB-1)*mp_BlockSize + IT
               if(J .le. N) then
                  Res(J) = V(J) + Scalar
               end if
         return
   end subroutine AddScalar_DevVec_I_KERNEL0
   !**************************************************************************

   !**************************************************************************
   subroutine AddScalar_DevVec_I_template0(IDEV, N, V, Val, Res)
   !***  PURPOSE:   to add a value to a vector
   !                Res(1:N) = V(1:N) + Val
   !
   !    INPUT:     IDEV,     device id
   !               N,        size of the vector
   !               Val,      the valuse to be added 
   !               V,        the vector   
   !    OUTPUT     Res,      the resulting  vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,         intent(in)   :: IDEV
        integer,         intent(in)   :: N
        integer,         intent(in)   :: Val
        integer, device, dimension(:) :: V, Res
        !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR
 
            
              ERR     = cudaGetDevice(CURDEV)
              blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
              threads = dim3(mp_BlockSize, 1, 1)

              ERR  = cudaSetDevice(IDEV)
              call AddScalar_DevVec_I_KERNEL0<<<blocks, threads>>>(N, V, Val, Res)
              ERR  = cudaSetDevice(CURDEV)
              return
    end subroutine AddScalar_DevVec_I_template0
   !***************************************************************************   

   !**************************************************************************
   subroutine AddScalar_DevVec_I_template0_a(V, Val)
   !***  PURPOSE:   to add a value to a vector
   !                V(1:N) = V(1:N) + Val
   !
   !    INPUT:     V,        the typedefined vector
   !               Val,      the valuse to be added 
   !    OUTPUT     V,      the resulting  vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,         intent(in) :: Val
        class(DevVec_I)             :: V
      !--- Device variables and variables to be used in GPU
      
      
              if(V%IDEV .lt. 0 ) return
 
              call AddScalar_DevVec_I_template0(V%IDEV, size(V%Data), V%Data, Val, V%Data)
              return
    end subroutine AddScalar_DevVec_I_template0_a
   !***************************************************************************    

  !***************************************************************************
   attributes(global) subroutine Add_DevVec_DF_KERNEL0(N, V1, V2, Res)
   !***  PURPOSE:   KERNEL    to add two mats
   !                          Res(1:N) = V1(1:N) + V2(1:N)
   !
   !$$   INPUT:     N         
   !$$              V1, V2    
   !$$   OUTPUT     Res,      
   !
   implicit none
   !----   DUMMY Variables
           integer,     value               :: N
           real(KINDDF),device, dimension(:):: V1, V2, Res
   
   !----   Local variables
           integer :: IT, IB, J, K
   
               IT   = (threadidx%y-1)*blockdim%x + threadidx%x
               IB   = (blockidx%y-1) * griddim%x +  blockidx%x
               J     = (IB-1)*mp_BlockSize + IT
               if(J .le. N) then
                  Res(J) = V1(J) + V2(J)
               end if
         return
   end subroutine Add_DevVec_DF_KERNEL0
   !**************************************************************************

   !**************************************************************************
   subroutine Add_DevMat_DF_template0(IDEV, X1, X2, Y1, Y2, V1, V2, Res)
   !***  PURPOSE:   to add two mats to  scalar to a vector:
   !                Res(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
   !
   !$$   INPUT:     X1, X2, Y1,Y2,  the bounds of the array
   !$$              V1, V2          the array
   !$$
   !$$   OUTPUT     Res,            the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,      intent(in)            :: IDEV, X1, X2, Y1, Y2
        real(KINDDF), device, dimension(:,:):: V1, V2, Res
        !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, Y, N
 
              ERR     = cudaGetDevice(CURDEV)
              N       = X2 - X1 + 1
              blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
              threads = dim3(mp_BlockSize, 1, 1)

              ERR     = cudaSetDevice(IDEV)
              do Y = Y1, Y2
                  call   Add_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, V1(X1:,Y), V2(X1:,Y), Res(X1:,Y))
              end do    
              ERR   = cudaSetDevice(CURDEV)
              return
    end subroutine Add_DevMat_DF_template0
   !*************************************************************************** 
    
   !**************************************************************************
   subroutine Add_DevMat_DF_template0_a(X1, X2, Y1, Y2, V1, V2, Res)
   !***  PURPOSE:  to add two mats to  scalar to a vector:
   !               Res(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
   !
   !$$   INPUT:     X1, X2, Y1, Y2,  the bounds of the array
   !$$              V1, V2,          the arraies
   !$$   OUTPUT     Res,             
   !   
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: X1, X2, Y1, Y2
           class(DevMat_DF)          :: V1, V2, Res
           !--- 
    
                 if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0) return
                 if(V1%IDEV .ne. V2%IDEV) then
                  write(*,*) "Error:: the two arries allocated on different devices ot be add in MSM_MutltiGPU"
                  stop
                 end if  

                 if(Res%IDEV .lt. 0) then
                    call Allocate_DevMat_DF_0 (V1%IDEV, Res, (/size(V1%Data,dim=1), size(V1%Data,dim=2)/))
                 end if  
   
                 call Add_DevMat_DF_template0(Res%IDEV, X1, X2, Y1, Y2, V1%Data, V2%Data, Res%Data)
                 return
       end subroutine Add_DevMat_DF_template0_a
   !***************************************************************************

   !**************************************************************************
   subroutine Add_DevMat_DF_template0_b(X1, X2, V1, V2, Res)
   !***  PURPOSE:   to add two mats to  scalar to a vector:
   !                Res(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
   !
   !$$   INPUT:     X1, X2,    the bounds of the array
   !$$              V1, V2     the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        class(DevMat_DF), dimension(:)             :: V1, V2, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Add_DevMat_DF_template0_a(X1(I), X2(I), 1,3, V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine Add_DevMat_DF_template0_b
   !***************************************************************************

   !**************************************************************************
   subroutine Add_DevMat_DF_template0_c(X2, V1, V2, Res)
   !***  PURPOSE:   to add two mats to  scalar to a vector:
   !                Res(1:X2,Y1:Y2) = V1(1:X2,Y1:Y2) + V2(1:X2,Y1:Y2)
   !
   !$$   INPUT:     X2,        the bounds of the array
   !$$              V1, V2     the array
   !$$
   !$$   OUTPUT     Res,       the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        class(DevMat_DF), dimension(:)             :: V1, V2, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Add_DevMat_DF_template0_a(1,X2(I), 1,3, V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine Add_DevMat_DF_template0_c
   !***************************************************************************

   !**************************************************************************
   subroutine Add_DevMat_DF_template1_a(X1, X2, Y1, Y2, V1, V2)
   !***  PURPOSE:  to add two mats to  scalar to a vector:
   !               V2(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
   !
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: X1, X2, Y1, Y2
           class(DevMat_DF)          :: V1, V2
           !--- 
                 call Add_DevMat_DF_template0_a(X1, X2, Y1, Y2, V1, V2, V2)
                 return
       end subroutine Add_DevMat_DF_template1_a
   !***************************************************************************

   !**************************************************************************
   subroutine Add_DevMat_DF_template1_b(X1, X2, V1, V2)
   !***  PURPOSE:  to add two mats to  scalar to a vector:
   !                V2(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        class(DevMat_DF), dimension(:)             :: V1, V2
        !--- 
        
              call Add_DevMat_DF_template0_b(X1, X2, V1, V2, V2)
              return
    end subroutine Add_DevMat_DF_template1_b
   !***************************************************************************

   !**************************************************************************
   subroutine Add_DevMat_DF_template1_c(X2, V1, V2)
   !***  PURPOSE:   to add two mats to  scalar to a vector:
   !                V2(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
   !
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        class(DevMat_DF), dimension(:)             :: V1, V2
        
              call Add_DevMat_DF_template0_c(X2, V1, V2, V2)
              return
    end subroutine Add_DevMat_DF_template1_c
  !***************************************************************************

  !***************************************************************************
   attributes(global) subroutine Minus_DevVec_DF_KERNEL0(N, V1, V2, Res)
   !***  PURPOSE:   KERNEL    to minus two mats:
   !                          Res(1:N) = V1(1:N) - V2(1:N)
   !
   !$$   INPUT:     N             the bounds of the array
   !$$              V1, V2        the array
   !$$
   !$$   OUTPUT     Res,          the scaled vector
   !
   implicit none
   !----   DUMMY Variables
           integer,     value               :: N
           real(KINDDF),device, dimension(:):: V1, V2, Res
   
   !----   Local variables
           integer :: IT, IB, J
   
               IT   = (threadidx%y-1)*blockdim%x + threadidx%x
               IB   = (blockidx%y-1) * griddim%x +  blockidx%x
               J     = (IB-1)*mp_BlockSize + IT
               if(J .le. N) then
                  Res(J) = V1(J) - V2(J)
               end if

         return
   end subroutine Minus_DevVec_DF_KERNEL0
   !**************************************************************************

   !**************************************************************************
   subroutine Minus_DevMat_DF_template0(IDEV, X1, X2, Y1, Y2, V1, V2, Res)
   !***  PURPOSE:   to add two mats to  scalar to a vector:
   !                Res(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) - V2(X1:X2,Y1:Y2)
   !
   !$$   INPUT:     X1, X2, Y1,Y2,  the bounds of the array
   !$$              V1, V2          the array
   !$$
   !$$   OUTPUT     Res,            the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,      intent(in)            :: IDEV, X1, X2, Y1, Y2
        real(KINDDF), device, dimension(:,:):: V1, V2, Res
        !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, N, Y
 
              ERR     = cudaGetDevice(CURDEV)
              ERR     = cudaSetDevice(IDEV)

              N       = X2 - X1 + 1
              blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
              threads = dim3(mp_BlockSize, 1, 1)
              do Y=Y1, Y2              
                 call   Minus_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, V1(X1:,Y), V2(X1:,Y), Res(X1:,Y))
              end do   
              ERR    = cudaSetDevice(CURDEV)
 
              return
    end subroutine Minus_DevMat_DF_template0
   !*************************************************************************** 
    
   !**************************************************************************
   subroutine Minus_DevMat_DF_template0_a(X1, X2, Y1, Y2, V1, V2, Res)
   !***  PURPOSE:  to calculate:
   !               Res(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) - V2(X1:X2,Y1:Y2)
   !
   !$$   INPUT:     Dimx, Dimy,  the bounds of the array
   !$$              V1, V2       the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: X1, X2, Y1, Y2
           class(DevMat_DF)          :: V1, V2, Res
           !--- 
    
                 if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0) return
                 if(V1%IDEV .ne. V2%IDEV) then
                  write(*,*) "Error:: the two arries allocated on different devices ot be add in MSM_MutltiGPU"
                  stop
                 end if  

                 if(Res%IDEV .lt. 0) then
                    call Allocate_DevMat_DF_0 (V1%IDEV, Res, (/size(V1%Data,dim=1), size(V1%Data,dim=2)/))
                 end if  
   
                 call Minus_DevMat_DF_template0(Res%IDEV, X1, X2, Y1, Y2, V1%Data, V2%Data, Res%Data)
                 return
       end subroutine Minus_DevMat_DF_template0_a
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevMat_DF_template0_b(X1, X2, V1, V2, Res)
   !***  PURPOSE:   to caluate :
   !                Res(X1:X2,Y1:Y2) = V1(X1:X2,1:3) - V2(X1:X2,1:3)
   !
   !$$   INPUT:     X1, X2,    the bounds of the array
   !$$              V1, V2     the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        class(DevMat_DF), dimension(:)             :: V1, V2, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Minus_DevMat_DF_template0_a(X1(I), X2(I), 1,3, V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine Minus_DevMat_DF_template0_b
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevMat_DF_template0_c(X2, V1, V2, Res)
   !***  PURPOSE:   to calculate:
   !                Res(1:X2,1:3) = V1(1:X2,1:3) - V2(1:X2,1:3)
   !
   !$$   INPUT:     X2,        the bounds of the array
   !$$              V1, V2     the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        class(DevMat_DF), dimension(:)             :: V1, V2, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Minus_DevMat_DF_template0_a(1,X2(I), 1, 3, V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine Minus_DevMat_DF_template0_c
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevMat_DF_template1_a(X1, X2, Y1, Y2, V1, V2)
   !***  PURPOSE:  to calculate:
   !               V2(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) - V2(X1:X2,Y1:Y2)
   !
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: X1, X2, Y1, Y2
           class(DevMat_DF)          :: V1, V2
           !--- 
                 call Minus_DevMat_DF_template0_a(X1, X2, Y1, Y2, V1, V2, V2)
                 return
       end subroutine Minus_DevMat_DF_template1_a
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevMat_DF_template1_b(X1, X2, V1, V2)
   !***  PURPOSE:  to calculate:
   !                V2(X1:X2,1:3) = V1(X1:X2,1:3) - V2(X1:X2,1:3)
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        class(DevMat_DF), dimension(:)             :: V1, V2
        !--- 
        
              call Minus_DevMat_DF_template0_b(X1, X2, V1, V2, V2)
              return
    end subroutine Minus_DevMat_DF_template1_b
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevMat_DF_template1_c(X2, V1, V2)
   !***  PURPOSE:   to calculate:
   !                V2(1:X2,1:3) = V1(1:X2,1:3) + V2(1:X2,1:3)
   !
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        class(DevMat_DF), dimension(:)             :: V1, V2
        
              call Minus_DevMat_DF_template0_c(X2, V1, V2, V2)
              return
    end subroutine Minus_DevMat_DF_template1_c
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevVec_DF_template0(IDEV, X1, X2, V1, V2, Res)
   !***  PURPOSE:   to calculate:
   !                Res(X1:X2) = V1(X1:X2) - V2(X1:X2)
   !
   !$$   INPUT:     X1, X2,   the bounds of the array
   !$$              V1, V2    the array
   !$$   OUTPUT     Res,      the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,      intent(in)          :: IDEV, X1, X2
        real(KINDDF), device, dimension(:):: V1, V2, Res
        !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, N
 
              ERR     = cudaGetDevice(CURDEV)
              ERR     = cudaSetDevice(IDEV)

              N       = X2 - X1 + 1
              blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
              threads = dim3(mp_BlockSize, 1, 1)
                call   Minus_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, V1(X1:), V2(X1:), Res(X1:))
              ERR    = cudaSetDevice(CURDEV)
 
              return
    end subroutine Minus_DevVec_DF_template0
   !*************************************************************************** 
    
   !**************************************************************************
   subroutine Minus_DevVec_DF_template0_a(X1, X2, V1, V2, Res)
   !***  PURPOSE:  to calculate:
   !               Res(X1:X2) = V1(X1:X2) - V2(X1:X2)
   !
   !$$   INPUT:     Dimx, Dimy,  the bounds of the array
   !$$              V1, V2       the array
   !$$   OUTPUT     Res,         the scaled vector
   !   
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: X1, X2
           class(DevVec_DF)          :: V1, V2, Res
           !--- 
    
                 if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0) return
                 if(V1%IDEV .ne. V2%IDEV) then
                  write(*,*) "Error:: the two arries allocated on different devices ot be add in MSM_MutltiGPU"
                  stop
                 end if  

                 if(Res%IDEV .lt. 0) then
                    call Allocate_DevVec_DF_0 (V1%IDEV, Res, size(V1%Data))
                 end if  
   
                 call Minus_DevVec_DF_template0(Res%IDEV, X1, X2, V1%Data, V2%Data, Res%Data)
                 return
       end subroutine Minus_DevVec_DF_template0_a
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevVec_DF_template0_b(X1, X2, V1, V2, Res)
   !***  PURPOSE:   to caluate :
   !                Res(X1:X2) = V1(X1:X2) - V2(X1:X2)
   !
   !$$   INPUT:     X1, X2,    the bounds of the array
   !$$              V1, V2     the array
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        class(DevVec_DF), dimension(:)             :: V1, V2, Res
        !---local variable
        integer::I

              do I = 1, m_NDEVICE
                 call Minus_DevVec_DF_template0_a(X1(I), X2(I), V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine Minus_DevVec_DF_template0_b
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevVec_DF_template0_c(X2, V1, V2, Res)
   !***  PURPOSE:   to calculate:
   !                Res(1:X2) = V1(1:X2) - V2(1:X2)
   !
   !$$   INPUT:     X2,        the bounds of the array
   !$$              V1, V2     the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        class(DevVec_DF), dimension(:)             :: V1, V2, Res
        !--- Device variables and variables to be used in GPU
        integer::I

              do I = 1, m_NDEVICE
                 call Minus_DevVec_DF_template0_a(1,X2(I), V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine Minus_DevVec_DF_template0_c
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevVec_DF_template1_a(X1, X2, V1, V2)
   !***  PURPOSE:  to calculate:
   !               V2(X1:X2) = V1(X1:X2) - V2(X1:X2)
   !
         implicit none
          !----   DUMMY Variables
           integer,      intent(in)  :: X1, X2
           class(DevVec_DF)          :: V1, V2
           !--- 
                 call Minus_DevVec_DF_template0_a(X1, X2,  V1, V2, V2)
                 return
       end subroutine Minus_DevVec_DF_template1_a
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevVec_DF_template1_b(X1, X2, V1, V2)
   !***  PURPOSE:  to calculate:
   !               V2(X1:X2) = V1(X1:X2) - V2(X1:X2)
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X1, X2
        class(DevVec_DF), dimension(:)             :: V1, V2
        !--- 
        
              call Minus_DevVec_DF_template0_b(X1, X2, V1, V2, V2)
              return
    end subroutine Minus_DevVec_DF_template1_b
   !***************************************************************************

   !**************************************************************************
   subroutine Minus_DevVec_DF_template1_c(X2, V1, V2)
   !***  PURPOSE:   to add two mats to  scalar to a vector:
   !                V2(X1:X2,Y1:Y2) = V1(X1:X2,Y1:Y2) + V2(X1:X2,Y1:Y2)
   !
      implicit none
       !----   DUMMY Variables
        integer,          dimension(:), intent(in) :: X2
        class(DevVec_DF), dimension(:)             :: V1, V2
        
              call Minus_DevVec_DF_template0_c(X2, V1, V2, V2)
              return
    end subroutine Minus_DevVec_DF_template1_c
   !***************************************************************************

   !***************************************************************************
   attributes(global) subroutine AddBD_DevVec_DF_KERNEL0(N, LB, HB, V1, V2, Res)
   !***  PURPOSE:   KERNEL    to add two mats with boundary:
   !                          Res(1:N) = V1(1:N) + V2(1:N)
   !
   !$$   INPUT:     X1, X2,   the bounds of the array
   !$$              V,            the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !
   implicit none
   !----   DUMMY Variables
           integer,      value                 :: N
           real(KINDDF), value                 :: LB, HB
           real(KINDDF), device, dimension(:)  :: V1, V2, Res
   
   !----   Local variables
           integer :: IT, IB, J
           real(KINDDF) :: RT
   
               IT   = (threadidx%y-1)*blockdim%x + threadidx%x
               IB   = (blockidx%y-1) * griddim%x +  blockidx%x
               J    = (IB-1)*mp_BlockSize + IT
               
               if(J .le. N) then
                  RT = V1(J) + V2(J)
                  if(RT .gt. HB) then
                     RT = RT - (HB - LB)
                  else if(RT .lt. LB) then
                     RT = RT + (HB - LB)   
                  end if 
                  Res(J) = RT
               end if
         return
   end subroutine AddBD_DevVec_DF_KERNEL0
   !**************************************************************************

   !**************************************************************************
   subroutine AddBD_DevMat_DF_template0(IDEV, Shift, N, Y1, Y2, LB, HB, V1, V2, Res)
   !***  PURPOSE:   to add two mats with boundary
   !                Res(1+Shift:N+Shfit,Y1:Y2) = V1(1+Shift:N+Shfit,Y1:Y2) + V2(1+Shift:N+Shfit,Y1:Y2)
   !
   !$$   INPUT:     N,              the bounds of the array
   !$$              V1, V2          the array
   !$$
   !$$   OUTPUT     Res,            the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,      intent(in)               :: IDEV, Shift, N, Y1, Y2
        real(KINDDF), intent(in),dimension(:)  :: LB, HB 
        real(KINDDF), device,    dimension(:,:):: V1, V2, Res
        !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer    :: CURDEV, ERR, Y
 
              ERR     = cudaGetDevice(CURDEV)
 
              blocks  = dim3((N-1)/mp_BlockSize+1, 1, 1)
              threads = dim3(mp_BlockSize, 1, 1)

              ERR     = cudaSetDevice(IDEV)
              do Y = Y1, Y2
                 call   AddBD_DevVec_DF_KERNEL0<<<blocks, threads>>>(N, LB(Y), HB(Y), V1(1:,Y), V2(Shift:,Y), Res(Shift:,Y))
              end do   
              ERR   = cudaSetDevice(CURDEV)
 
              return
   end subroutine AddBD_DevMat_DF_template0
   !***************************************************************************

   !**************************************************************************
   subroutine AddBD_DevMat_DF_template0_a(Shift, N, Y1, Y2, LB, HB, V1, V2, Res)
   !***  PURPOSE:  to add two mats to  scalar to a vector:
   !               Res(1+Shift:N+Shfit,Y1:Y2) = V1(1+Shift:N+Shfit,Y1:Y2) + V2(1+Shift:N+Shfit,Y1:Y2)
   !
   !$$   INPUT:     Dimx, Dimy,  the bounds of the array
   !$$              LB, HB,      the value boundary
   !$$              V1, V2,      the arraies
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
      !----   DUMMY Variables
         integer,      intent(in)               :: Y1, Y2, Shift, N
         real(KINDDF), intent(in), dimension(:) :: LB, HB
         class(DevMat_DF)                       :: V1, V2, Res
         !--- 
    
                if(V1%IDEV .lt. 0 .or. V2%IDEV .lt. 0) return
                if(V1%IDEV .ne. V2%IDEV) then
                  write(*,*) "Error:: the two arries allocated on different devices ot be add in MSM_MutltiGPU"
                  stop
                end if 
                
                if(size(LB) .lt. Y2 .or. size(HB) .lt. Y2) then
                  write(*,*) "Error:: the vector for boundary condition is smaller than required size in MSM_MutltiGPU"
                  stop
                endif 

                if(Res%IDEV .lt. 0) then
                   call Allocate_DevMat_DF_0 (V2%IDEV, Res, (/size(V2%Data,dim=1), size(V2%Data,dim=2)/))
                end if  
   
                call AddBD_DevMat_DF_template0(Res%IDEV, Shift, N, Y1, Y2, LB, HB, V1%Data, V2%Data, Res%Data)
                return
   end subroutine AddBD_DevMat_DF_template0_a
   !***************************************************************************

   !**************************************************************************
   subroutine AddBD_DevMat_DF_template0_b(Shift, N, LB, HB, V1, V2, Res)
   !***  PURPOSE:   to add two mats with boundary condition:
   !                Res(1+Shift:N+Shfit,1:3) = V1(1+Shift:N+Shfit,1:3) + V2(1+Shift:N+Shfit,1:3)
   !
   !$$   INPUT:     X1, X2,    the bounds of the array
   !$$              V1, V2     the array
   !$$
   !$$   OUTPUT     Res,         the scaled vector
   !   
      implicit none
       !----   DUMMY Variables
        integer,          intent(in), dimension(:) :: Shift, N
        real(KINDDF),     intent(in), dimension(:) :: LB, HB
        class(DevMat_DF),             dimension(:) :: V1, V2, Res
        !--- Device variables and variables to be used in GPU
        integer::I
        
              do I = 1, m_NDEVICE
                 call AddBD_DevMat_DF_template0_a(Shift(I), N(I), 1,3, LB, HB, V1(I), V2(I), Res(I))
              end do   
              return
    end subroutine AddBD_DevMat_DF_template0_b
   !***************************************************************************

   !**************************************************************************
   subroutine AddBD_DevMat_DF_template1_a(Shift, N, Y1, Y2, LB, HB, V1, V2)
   !***  PURPOSE:  to add two mats to  scalar to a vector:
   !               V2(1+Shift:N+Shfit,Y1:Y2) = V1(1+Shift:N+Shfit,Y1:Y2) + V2(1+Shift:N+Shfit,Y1:Y2)
   !
         implicit none
          !----   DUMMY Variables
           integer,         intent(in)               :: Shift, N, Y1, Y2
           real(KINDDF),    intent(in), dimension(:) :: LB, HB
           class(DevMat_DF)                          :: V1, V2
           !--- 
                 call AddBD_DevMat_DF_template0_a(Shift, N, Y1, Y2, LB, HB, V1, V2, V2)
                 return
       end subroutine AddBD_DevMat_DF_template1_a
   !***************************************************************************

   !**************************************************************************
   subroutine AddBD_DevMat_DF_template1_b(Shift, N, LB, HB, V1, V2)
   !***  PURPOSE:  to add two mats to  scalar to a vector:
   !               V2(1+Shift:N+Shfit,1:3) = V1(1+Shift:N+Shfit,1:3) + V2(1+Shift:N+Shfit,1:3)
      implicit none
       !----   DUMMY Variables
        integer,          intent(in), dimension(:) :: Shift, N
        real(KINDDF),     intent(in), dimension(:) :: LB, HB
        class(DevMat_DF),             dimension(:) :: V1, V2
        !--- 
        
              call AddBD_DevMat_DF_template0_b(Shift, N, LB, HB, V1, V2, V2)
              return
    end subroutine AddBD_DevMat_DF_template1_b
   !***************************************************************************


  end module MSM_MultiGPU_Basic

