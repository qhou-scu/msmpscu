module MSM_MemAllocator_GPU
  ! ***  DESCRIPTION: this file provides a number of template routines for allocating memeory on CPU, 
  !                   which are parts of file MCLIB_Utilities.F90 originally written by Zhai Lei.
  !    
  !                   see also: MSM_MemAllocate_Temp.F90
  !                   
  !                  ______________________________________________________
  !
  ! **** HOSTORY:
  !       * Sep.    2018(Zhai Lei): original version of MCLIB_Utilities.F90 
  !                                 
  !
  !       * Nov.    2018(HOU Qing): modify abd merge to MSMPLIB
  !            
  !     
  use MSM_CONSTANTS 
  implicit none

  private:: ResizeArrayi_OneDim,  &
            ResizeArrayr_OneDim,  &
            ResizeArrayd_OneDim,  &
            ResizeArrayi_TwoDim,  &
            ResizeArrayr_TwoDim,  &
            ResizeArrayd_TwoDim
  interface ResizeArray_GPU
    module procedure ResizeArrayi_OneDim
    module procedure ResizeArrayr_OneDim
    module procedure ResizeArrayd_OneDim

    module procedure ResizeArrayi_TwoDim
    module procedure ResizeArrayr_TwoDim
    module procedure ResizeArrayd_TwoDim

  end interface ResizeArray_GPU

  !-----------------------
  private:: AllocateOneDimi,            &
            AllocateOneDimr,            &
            AllocateOneDimd,            &
            AllocateTwoDimi,            &
            AllocateTwoDimr,            &
            AllocateTwoDimd,            &
            AllocateThreeDimi,          &
            AllocateThreeDimr,          &
            AllocateThreeDimd
  interface AllocateArray_GPU
    module procedure AllocateOneDimi
    module procedure AllocateOneDimr
    module procedure AllocateOneDimd

    module procedure AllocateTwoDimi
    module procedure AllocateTwoDimr
    module procedure AllocateTwoDimd

    module procedure AllocateThreeDimi
    module procedure AllocateThreeDimr
    module procedure AllocateThreeDimd

  end interface AllocateArray_GPU

  !------------------
  private:: DeAllocateOneDimi,             &
            DeAllocateOneDimr,             &
            DeAllocateOneDimd,             &
                                                
            DeAllocateTwoDimi,             &
            DeAllocateTwoDimr,             &
            DeAllocateTwoDimd,             &
                                           
            DeAllocateThreeDimi,           &
            DeAllocateThreeDimr,           &
            DeAllocateThreeDimd
  interface DeAllocateArray_GPU
    module procedure DeAllocateOneDimi
    module procedure DeAllocateOneDimr
    module procedure DeAllocateOneDimd

    module procedure DeAllocateTwoDimi
    module procedure DeAllocateTwoDimr
    module procedure DeAllocateTwoDimd

    module procedure DeAllocateThreeDimi
    module procedure DeAllocateThreeDimr
    module procedure DeAllocateThreeDimd

  end interface DeAllocateArray_GPU

  private:: DuplicateArrayi_OneDim,        &
            DuplicateArrayr_OneDim,        &
            DuplicateArrayd_OneDim,        &
            
            DuplicateArrayi_TwoDim,        &
            DuplicateArrayr_TwoDim,        &
            DuplicateArrayd_TwoDim
  interface DuplicateArray_Host
    module procedure  DuplicateArrayi_OneDim
    module procedure  DuplicateArrayr_OneDim
    module procedure  DuplicateArrayd_OneDim
    module procedure  DuplicateArrayi_TwoDim
    module procedure  DuplicateArrayr_TwoDim
    module procedure  DuplicateArrayd_TwoDim  
  end interface DuplicateArray_Host

  contains
  #define  ARRAYDEV  device,
  #include "MSM_MemAllocate_Temp.F90" 

end module MSM_MemAllocator_GPU
