  module MSM_MultiGPU_DatList_IVec
  !***  DESCRIPTION: this module is to define a number of operations on
  !                  the device data list, the extension data types,
  !                  defined in MSM_MultiGPU_Basic.F90. 
  !                  ______________________________________________________
  !
  ! **** HOSTORY:
  !       * Jan.    2019(HOU Qing): first version
  !                                 
  ! 
  !
  use MSM_MultiGPU_Basic
  implicit none
  !---------------------------------------------------------
          #define LIISTNAME DevVec_I_List
          #define DATATYPE  type(DevVec_I),dimension(:) 
          #define DEALLOCTE_DATA  DevDeallocate
          #include "MSM_MultiGPU_DatList_Header.F90"
 !****************************************************************************
  contains
          #include "MSM_MultiGPU_DatList_Body.F90"          
  end module MSM_MultiGPU_DatList_IVec

