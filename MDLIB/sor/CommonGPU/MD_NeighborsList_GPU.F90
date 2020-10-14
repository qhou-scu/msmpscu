  module MD_NeighborsList_GPU
  !**** DESCRIPTION: to calculate the neighbor list using GPU parallel calculation.
  !                   It is possible that multiple devices can be used.
  !     NOTE:         This is adapted form MD_NeighborsList_GPU,
  !                   but the algorithm is different, the meaning of the niegbore
  !                   list has also been changed. The force calculation should be
  !                   changed correspondingly.
  !
  !
  !                  SEE ALSO:  MD_NeighborsList_GPU
  !                             modification notes in MD_Globle_Variables_GPU

  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2012
  !
  !**** HISTORY:
  !       1. Mar.,   2012 (HOU Qing): adapted from  MD_NeighborsList_GPU, 2012-04 (Hou Qing)
  !       2. Oct 19, 2016 (HOU Qing): Add the calculation of hm_GIDINV, which is the inverse transformation
  !                                   of hm_GID.
  !       3. June 17, 2018 (HOU Qing): Add interfaces Cal_NeighBoreList_A_DEV and Cal_NeighBoreList_B_DEV.
  !                                    _A version indicating the dummy input is simulation boxes array
  !                                    _B version indicating the dummy input is single simulation box
  !
  !       * Sep.27, 2018(HOU Qing): 
  !                                 Add a new data type MDDEVNeighbor_list, the dxm_KVOIS, in old version
  !                                 are replaced by the MDDEVNeighbor_list type array dm_Neighbors.
  !                                 Correspondingly, the dxm_XXX arraies in other modules using
  !                                 this module are replaced by arraies of extend data types.
  !                                 See also the modification made in MD_Globle_Variables_GPU.F90
  !       * Jan.29, 2019(HOU Qing): 
  !                                 Using m_NAAP, the number of ACTIVATIVE atom in cells, in neighbore-list
  !                                 calculation. If m_NAAP of a cell is 0, the neighors of the atoms in the 
  !                                 cell are not neccearyly calculated and thus may improve calculation 
  !                                 efficiency. 
  !       * Aug.29, 2019(HOU Qing): 
  !                                 Add new interfaces for Reorder_NeighBoreList_Nearest_Dev.
  !                                 The old Reorder_NeighBoreList_Nearest_Dev was replaced by Reorder_NeighBoreList_Nearest_0,
  !                                 a routine Reorder_NeighBoreList_Nearest_1 was added
  !                                 An interface Reorder_NeighBoreList_Nearest_Dev was added that is interfaced to these
  !                                 routines .
  !                                   
                                    
 
  !****
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_NeighborsList
  use MD_Globle_Variables_GPU
  
  implicit none

           !--- C2050
           integer, parameter, private::mp_BLOCKSIZE = 128
           !--- C1060
           !integer, parameter, private::mp_BLOCKSIZE = 256

           !---
           integer, parameter, private::mp_NNC=27
           integer, parameter, private::mp_MXNEAREST = 512                  ! the max permitted number of nearest neighbors
                                                                            ! used in kernal of extracting nearest neighbors
           real(KINDDF), private::m_LBOX(3)
           real(KINDDF), private::m_BOXSIZE(3)
           real(KINDDF), private::m_BOXSHAPE(3,3)
           real(KINDDF), private::m_RM2(mp_MXGROUP,mp_MXGROUP)
           integer,      private::m_PD(3) = 1
           integer,      private::m_mxKVOIS= 256                            ! the max number of neighbore permitted

           real(KINDSF), constant,private::dcm_BOXSHAPE(3,3)
           real(KINDSF), constant,private::dcm_RM2(mp_MXGROUP,mp_MXGROUP)
           integer,      constant,private::dcm_NIX(mp_NNC)
           integer,      constant,private::dcm_NIY(mp_NNC)
           integer,      constant,private::dcm_NIZ(mp_NNC)
           integer,               private::m_CPYSTRDEVS(m_MXDEVICE)=0
           !--- the structure storing the list on device:
           type:: MDDEVNeighbor_list
                  integer::NPRT    = 0                                   !The total number of atoms, dm_NPRT
                  integer::NPRTB   = 0                                   !The number of  atoms in each box
                  integer::MXKVOIS = 256
                  integer::STARTA(m_MXDEVICE)   = 0                      !The atom ID of the first atom handled by a device
                  integer::ENDA(m_MXDEVICE)     = 0                      !The atom ID of the last  atom handled by a device
                  integer::NPA(m_MXDEVICE)      = 0                      !The number of atoms acutally handled by a device 
                  type(DevVec_I), dimension(m_MXDEVICE)::KVOIS
                  type(DevMat_I), dimension(m_MXDEVICE)::INDI
           end    type MDDEVNeighbor_list
           type(MDDEVNeighbor_list)::dm_Neighbors

           !-- varables used when head-linked cell technique to be used
           !
           integer, private::m_NCELL(3)                                 !the number of cells in three dimension
           integer, private::m_NC0                                      !the number of cells in a simulation box
           integer, private::m_NC                                       !the number of cells in multiple boxes, if no multiple box is used, m_NC=m_NC0
           integer, private::m_MXNAC0                                   !the actual max number of atoms in a cell  (to be used by force calculation, set to public)

           integer, dimension(:), allocatable, private::m_HEAD
           integer, dimension(:), allocatable, private::m_LINK          ! link list
           integer, dimension(:), allocatable, private::m_INC           ! cell id for atoms
           integer, dimension(:), allocatable, private::hm_INC          ! cell id for atoms
           integer, private::m_NUMOUT =0                                ! the number of atoms out of box

           !--- temp memoery used for copyin/out the list on GPU
           type(NEIGHBOR_LIST), dimension(m_MXDEVICE), private::m_SwapList
           !integer, dimension(:),  allocatable, private::m_SwapKVOIS
           !integer, dimension(:,:),allocatable, private::m_SwapINDI

 !--------------------------------------------------------
 !    The GPU-kernels
           private:: Cal_NearestNeighbor_Kernel
           private:: Cal_NeighboreList_Kernel2C
           private:: Copy_From_Host_To_Devices_template
           private:: NeighboreList_IC_KERNEL
           private:: NeighboreList_IC_template0
           private:: NeighboreList_IC_template1
           private:: StartOnDevice_Reoder_template
           private:: StartOnDevice_template

 !----- Interfaces: 
 !---------------------------------------------------------
          private  :: Cal_NeighBoreList_A_DEV,   &
                      Cal_NeighBoreList_B_DEV
          public   :: Cal_NeighBoreList_DEV
          interface   Cal_NeighBoreList_DEV
             module procedure Cal_NeighBoreList_A_DEV
             module procedure Cal_NeighBoreList_B_DEV
          end interface

 !---------------------------------------------------------
          private  :: Cal_NeighBoreList2C_A_DEV,   &
                      Cal_NeighBoreList2C_B_DEV,   &
                      Cal_NeighBoreList2C_0_DEV
          interface   Cal_NeighBoreList2C_DEV
             module procedure Cal_NeighBoreList2C_A_DEV
             module procedure Cal_NeighBoreList2C_B_DEV
             module procedure Cal_NeighBoreList2C_0_DEV
          end interface
  !---------------------------------------------------------
          private  :: Clear_NeighboreList_0,  &
                      Clear_NeighboreList_1
          public   :: Clear_NeighboreList_DEV
          interface   Clear_NeighboreList_DEV
             module procedure Clear_NeighboreList_0
             module procedure Clear_NeighboreList_1  
          end interface

  !---------------------------------------------------------
          private:: &
                    Copyout_NeighboreList_A0,     &
                    Copyout_NeighboreList_A1,     &
                    Copyout_NeighboreList_A2,     &
                    Copyout_NeighboreList_B0,     &  
                    Copyout_NeighboreList_B1,     &
                    Copyout_NeighboreList_B2,     &   
                    Copyout_NeighboreList_B3
          public::  Copyout_NeighboreList_DEV
          interface Copyout_NeighboreList_DEV
             module procedure Copyout_NeighboreList_A0
             module procedure Copyout_NeighboreList_A1
             module procedure Copyout_NeighboreList_A2
             module procedure Copyout_NeighboreList_B0
             module procedure Copyout_NeighboreList_B1
             module procedure Copyout_NeighboreList_B2
             module procedure Copyout_NeighboreList_B3
          end interface        

  !---------------------------------------------------------
          public::  Get_Number_of_AtomPerCell

 !---------------------------------------------------------
          public::  Get_Number_of_Cells

  !---------------------------------------------------------
          public::  Get_Number_of_PermittedNeighbor

  !---------------------------------------------------------
          private::  Initial_Dev_Constants
          private::  Initialize_NeighboreList_A_DEV, &
                     Initialize_NeighboreList_B_DEV
          public ::  Initialize_NeighboreList_DEV
          interface  Initialize_NeighboreList_DEV
             module procedure Initialize_NeighboreList_A_DEV
             module procedure Initialize_NeighboreList_B_DEV
          end interface
          
  !---------------------------------------------------------
          private::  Reorder_NeighBoreList_Nearest_0,   &
                     Reorder_NeighBoreList_Nearest_0_A, &
                     Reorder_NeighBoreList_Nearest_1_A, &                    
                     Reorder_NeighBoreList_Nearest_0_B, &
                     Reorder_NeighBoreList_Nearest_1_B, &
                     Reorder_NeighBoreList_Nearest_0_C, &                     
                     Reorder_NeighBoreList_Nearest_A_C, &                     
                     Reorder_NeighBoreList_Nearest_B_C                     
          public ::  Reorder_NeighBoreList_Nearest_Dev
          interface  Reorder_NeighBoreList_Nearest_Dev
             module procedure Reorder_NeighBoreList_Nearest_0
             module procedure Reorder_NeighBoreList_Nearest_0_A
             module procedure Reorder_NeighBoreList_Nearest_1_A
             module procedure Reorder_NeighBoreList_Nearest_0_B
             module procedure Reorder_NeighBoreList_Nearest_1_B
             module procedure Reorder_NeighBoreList_Nearest_0_C
             module procedure Reorder_NeighBoreList_Nearest_A_C
             module procedure Reorder_NeighBoreList_Nearest_B_C
          end interface

  !---------------------------------------------------------
contains

  !****************************************************************************
  subroutine Initial_Dev_Constants()
  !***  PURPOSE:   to set constant parameters to be used
  !
  !     INPUT:     
  !     OUPUT:     dcm_NIX,dcm_NIY,dcm_NIZ
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
          integer::ERR, CURDEV, I
          integer::NIX(mp_NNC)=(/0,-1,-1,-1, 0, 0, -1, 1,-1, 0, 1,-1, 0, 1, 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/)
          integer::NIY(mp_NNC)=(/0, 0,-1, 1, 1, 0,  0, 0,-1,-1,-1, 1, 1, 1, 0, 1,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/)
          integer::NIZ(mp_NNC)=(/0, 0, 0, 0, 0, 1,  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
 

              ERR = cudaGetDevice(CURDEV)
              do I=1, m_NDEVICE 
                 ERR = cudaSetDevice(m_DEVICES(I)) 
                 dcm_NIX = NIX
                 dcm_NIY = NIY
                 dcm_NIZ = NIZ
                 !if(m_CPYSTRDEVS(I) .ne. 0) then
                 !   ERR = cudaStreamDestroy(m_CPYSTRDEVS(I))
                 !end if   
                 !ERR  = cudaStreamCreate(m_CPYSTRDEVS(I))
              end do 
                 
              ERR = cudaSetDevice(CURDEV)
          return
  end subroutine Initial_Dev_Constants
  !**********************************************************************************

  !**********************************************************************************
  subroutine Initialize_NeighboreList_A_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to allocate and initialize the device memories for device calucaltion
  !               of Neighbore Table
  !     INPUT:    SimBox,    the simulation box
  !               CtrlParam, the control parameters
  !               MULTIBOX, optional, indicating if multiple-box to be calculated conccurently
  !     OUTPUT:
  !
  !     NOTE:    Inialization of Neighborlikst calculation should be called after
  !              calling Initialize_Globle_Variables_DEV (see MD_Globle_Variables_GPU.F90)
  !
   implicit none
      !
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      type(SimMDCtrl)             ::CtrlParam

      !--- Local variables
      integer::NA0,ERR,K 
      real(KINDDF), parameter::EPS=0.0001


          !$$*** to allocate device memory
            !$$--- Check if the number of particle is
            !$$    consistent with that global device memeroy
            NA0 = SimBox(1)%NPRT*size(SimBox)
            if(dm_NPRT .LT. NA0) then
               write(*,*) "MDPSCU Error:number of particles in device is smaller than Simulation box"
               stop
            end if

           !$$--- clear the memory allocated before
           call Clear_NeighboreList_DEV()

           !$$--- allocate the device memory on device 0
           call Initial_Dev_Constants()

           !$$--- allocate memory for neighbore list on devices
           m_mxKVOIS            = CtrlParam%NB_MXNBS
           dm_Neighbors%NPRT    = NA0
           dm_Neighbors%NPRTB   = SimBox(1)%NPRT
           dm_Neighbors%MXKVOIS = m_mxKVOIS
           call DevAllocate(dm_Neighbors%KVOIS,   m_NAPDEV,              "NB_KOVIS")
           call DevAllocate(dm_Neighbors%INDI,  (/m_NAPDEV, m_mxKVOIS/), "NB_INDI")
           call DevSet(dm_Neighbors%KVOIS,  0)
           call DevSet(dm_Neighbors%INDI,   0)

           !$$---   to determine how many cells we need
           do K=1,3
              !$$if(CtrlParam%IFPD(K)) then
                 m_NCELL(K) = int(SimBox(1)%ZL(K)/(1.0D0*maxval(CtrlParam%NB_RM)) - EPS)
                 if(m_NCELL(K) .lt. C_ITHR) m_NCELL(K) = C_ITHR
              !$$--- we can discard these statements
              !$$else
              !$$   m_NCELL(K) = CtrlParam%NC(K)
              !$$   if(m_NCELL(K) .lt. C_ITHR) m_NCELL(K) = C_ITHR
              !$$end if
           end do

           m_NC0 = m_NCELL(1)*m_NCELL(2)*m_NCELL(3)
           m_NC  = m_NC0*size(SimBox)

          !$$--- to set the number of cells, memories defined in
          !$$    MD_Globle_Variables_GPU to be allocated
          call Set_BoxCell_Number(m_NC, m_NCELL)

          !$$--- allocate arrays on host for HEAD-LIST
          allocate(m_HEAD(m_NC),m_LINK(NA0),m_INC(NA0), hm_INC(NA0), STAT=ERR)
          if(ERR) then
             write(*,*) "MDPSCU Error: fail to allocate memery for linked-cell in Initialize_NeighboreList"
             stop
          endif

          hm_HasPart = mp_NOTPART

      return
  end subroutine Initialize_NeighboreList_A_DEV
  !**********************************************************************************

  !**********************************************************************************
  subroutine Initialize_NeighboreList_B_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to allocate and initialize the device memories for device calucaltion
  !               of Neighbore Table
  !     INPUT:    SimBox,    the simulation box
  !               CtrlParam, the control parameters
  !               MULTIBOX, optional, indicating if multiple-box to be calculated conccurently
  !     OUTPUT:
  !
  !     NOTE:    Inialization of Neighborlikst calculation should be called after
  !              calling Initialize_Globle_Variables_DEV (see MD_Globle_Variables_GPU.F90)
  !
   implicit none
      !
      !--- dummy variables
      type(SimMDBox)  ::SimBox
      type(SimMDCtrl) ::CtrlParam

      !--- Local variables
      integer::NA0,ERR, K 
      real(KINDDF), parameter::EPS=0.0001


          !$$*** to allocate device memory
            !$$--- Check if the number of particle is
            !$$    consistent with that global device memeroy
            NA0 = SimBox%NPRT
            if(dm_NPRT .LT. NA0) then
               write(*,*) "MDPSCU Error:number of particles in device is smaller than Simulation box"
               stop
            end if

           !$$--- clear the memory allocated before
           call Clear_NeighboreList_DEV()
           !$$--- allocate the device memory on device 0
           call Initial_Dev_Constants()

           !$$--- allocate memory for neighbore list on devices
           m_mxKVOIS            = CtrlParam%NB_MXNBS
           dm_Neighbors%NPRT    = NA0
           dm_Neighbors%NPRTB   = SimBox%NPRT
           dm_Neighbors%MXKVOIS = m_mxKVOIS
           call DevAllocate(dm_Neighbors%KVOIS,   m_NAPDEV,              "NB_KOVIS")
           call DevAllocate(dm_Neighbors%INDI,  (/m_NAPDEV, m_mxKVOIS/), "NB_INDI")
           call DevSet(dm_Neighbors%KVOIS,  0)
           call DevSet(dm_Neighbors%INDI,   0)

           !$$---   to determine how many cells we need
           do K=1,3
              !$$if(CtrlParam%IFPD(K)) then
                 m_NCELL(K) = int(SimBox%ZL(K)/(1.0D0*maxval(CtrlParam%NB_RM)) - EPS)
                 if(m_NCELL(K) .lt. C_ITHR) m_NCELL(K) = C_ITHR
              !$$--- we can discard these statements
              !$$else
              !$$   m_NCELL(K) = CtrlParam%NC(K)
              !$$   if(m_NCELL(K) .lt. C_ITHR) m_NCELL(K) = C_ITHR
              !$$end if
           end do

           m_NC0 = m_NCELL(1)*m_NCELL(2)*m_NCELL(3)
           m_NC = m_NC0

          !$$--- to set the number of cells, memories defined in
          !$$    MD_Globle_Variables_GPU to be allocated
          call Set_BoxCell_Number(m_NC, m_NCELL)

          !$$--- allocate arrays on host for HEAD-LIST
          allocate(m_HEAD(m_NC),m_LINK(NA0),m_INC(NA0), hm_INC(NA0), STAT=ERR)
          if(ERR) then
            write(*,*) "MDPSCU Error: fail to allocate memery for linked-cell in Initialize_NeighboreList"
            stop
          endif

          hm_HasPart = mp_NOTPART

      return
  end subroutine Initialize_NeighboreList_B_DEV
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_NeighboreList_0(List)
   !***  PURPOSE:  to deallocate the device memories allocated for
   !               device calucaltion of Neighbore Table
   !    INPUT:
   !
   !   OUTPUT:
   !
    implicit none
      type(MDDEVNeighbor_list)::List 
       !--- Local variables

       List%NPRT = 0 
       List%NPRTB = 0
       List%STARTA = 0
       List%ENDA = 0
       List%NPA = 0
       call DevDeallocate(List%KVOIS)  
       call DevDeallocate(List%INDI) 
       return
  end  subroutine Clear_NeighboreList_0
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_NeighboreList_1()
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calucaltion of Neighbore Table
  !    INPUT:
  !
  !   OUTPUT:
  !
   use MD_NeighborsList, only:Clear_NeighboreList
   implicit none

      !--- Local variables
      integer::ERR, CURDEV, I
      !---------------------------------------
           if(allocated(m_HEAD)) then
              deallocate(m_HEAD,m_LINK,m_INC, hm_INC, STAT=ERR)
              if(ERR) then
                 write(*,*) "MDPSCU Warning: fail to deallocate memory for HEAD-LIST"
                 call ONWARNING(gm_OnWarning)
              end if
           end if

           !ERR = cudaGetDevice(CURDEV)
           !do I=1, m_NDEVICE
           !   ERR = cudaSetDevice(m_DEVICES(I)) 
           !   if(m_CPYSTRDEVS(I) .ne. 0) then
           !      ERR = cudaStreamDestroy(m_CPYSTRDEVS(I))
           !    end if 
           !    m_CPYSTRDEVS(I) = 0   
           !end do
           !ERR = cudaSetDevice(CURDEV)

           call Clear_NeighboreList_0(dm_Neighbors)
           do I=1, size(m_SwapList)
              call Clear_NeighboreList(m_SwapList(I)) 
           end do   
           
           hm_HasPart = mp_NOTPART
      return
  end subroutine Clear_NeighboreList_1
  !**********************************************************************************

  !**********************************************************************************
  subroutine Copyout_NeighboreList_A0(hList, dList, GID)
  !***  PURPOSE:  to copy the neighbore list in device to a host neighbore list
  !
  !     INPUT:     dList,     the list on devices
  !                              
  !     OUTPUT     hList,     the list on host
  !
   implicit none
   !----   DUMMY Variables
           type(NEIGHBOR_LIST), dimension(:)::hList
           type(MDDEVNeighbor_list)         ::dList
           integer, dimension(:), optional  ::GID
   !----   Local variables
           integer::I, IB, J, J0, IA, IA0, IW, NP, ERR, STARTA, NPA, KVOIS
 
   !***
           NP = dList%NPRTB
           if(size(hList) .lt. dList%NPRT/dList%NPRTB) goto 200

         !$$--- determine the array size of the neighbore list
             do IB=1, size(hList)
                if(allocated(hList(IB)%KVOIS)) then
                   if(size(hList(IB)%KVOIS) .lt. NP)  deallocate(hList(IB)%KVOIS)
                 end if
                 if(allocated(hList(IB)%INDI)) then
                    if(size(hList(IB)%INDI) .lt. NP*dList%MXKVOIS) deallocate(hList(IB)%INDI)
                 end if

                 hList(IB)%MXKVOIS = dList%MXKVOIS
                 if(.not. allocated(hList(IB)%KVOIS) ) then
                     allocate(hList(IB)%KVOIS(NP) ,STAT=ERR)
                     if(ERR) goto 100
                 end if
                 if(.not. allocated(hList(IB)%INDI) ) then
                     allocate(hList(IB)%INDI(NP, dList%MXKVOIS), STAT=ERR)
                     if(ERR) goto 100
                 end if                 
                 !$$--- set List%KVOIS =0, in order there are atoms out of box
                 hList(IB)%KVOIS = 0
             end do

             do I=1, m_NDEVICE
                if(.not.allocated(m_SwapList(I)%KVOIS) ) then
                   allocate(m_SwapList(I)%KVOIS(m_NAPDEV), STAT=ERR)
                   if(ERR) goto 100
                 end if    
                 if(.not.allocated(m_SwapList(I)%INDI) ) then
                   allocate(m_SwapList(I)%INDI(m_NAPDEV,m_mxKVOIS), STAT=ERR)
                   if(ERR) goto 100
                 end if
             end do

             do I=1, m_NDEVICE
                call DevCopyOut(m_SwapList(I)%KVOIS, dList%KVOIS(I), m_NAPDEV) 
                call DevCopyOut(m_SwapList(I)%INDI,  dList%INDI(I),  m_NAPDEV*dList%MXKVOIS )
             end do
             call SynchronizeDevices()

             do I=1, m_NDEVICE    
                STARTA = dList%STARTA(I)
                NPA    = dList%NPA(I)
                 if(present(GID)) then
                    do J=1, NPA
                       IA    =  GID(J+STARTA-1)
                       IB    = (IA-1)/NP + 1
                       IA0   = IA - (IB-1)*NP
                       KVOIS = m_SwapList(I)%KVOIS(J)
                       hList(IB)%KVOIS(IA0) = KVOIS
                       do IW = 1, KVOIS
                          hList(IB)%INDI(IA0,IW) = GID(m_SwapList(I)%INDI(J,IW)) - (IB-1)*NP
                       end do
                    end do
                  else 
                     do J=1, NPA
                        IA =  J+STARTA-1
                        IB  = (IA-1)/NP + 1
                        IA0 = IA - (IB-1)*NP
                        KVOIS = m_SwapList(I)%KVOIS(J)
                        hList(IB)%KVOIS(IA0) = KVOIS
                        do IW = 1, KVOIS
                           hList(IB)%INDI(IA0,IW) = m_SwapList(I)%INDI(J,IW) - (IB-1)*NP
                        end do
                     end do
                  end if    
             end do
             return
 
   100       write(*,fmt="(A,  I7, A)")  " MDPSCU Error: fail to allocate memory for copyout neighbor list"
             write(*,fmt="(A)")          "               Process to be stopped"
             stop

   200       write(*,fmt="(A,  I7, A)")  " MDPSCU Error: the size of neighbo-list array is smaller than neighbor-list on devices"
             write(*,fmt="(A)")          "               Process to be stopped"
             stop
         return
  end subroutine Copyout_NeighboreList_A0
  !**********************************************************************************  

  !**********************************************************************************
  subroutine Copyout_NeighboreList_A1(hList,dList)
  !***  PURPOSE:  to copy the neighbore list in device to a host neighbore list
  !
  !     INPUT:     dList,     the list on devices
  !                              
  !     OUTPUT     hList,     the list on host
  !
   implicit none
   !----   DUMMY Variables
           type(NEIGHBOR_LIST), dimension(:)::hList
           type(MDDEVNeighbor_list)         ::dList
   !----   Local variables
 
   !***
             call Copyout_NeighboreList_A0(hList,dList, hm_GID)
             return
  end subroutine Copyout_NeighboreList_A1
  !**********************************************************************************  

  !**********************************************************************************
  subroutine Copyout_NeighboreList_A2(List)
  !***  PURPOSE:  to copy the neighbore list in device to a host neighbore list
  !               in old style. This routine can be used to test
  !               the neighbor calculation, it should be noted,
  !               in this GPU implement, the neighbor  list is
  !               the list with atoms are sorted by cells.
  !
  !     OUTPUT     List,      the neighbore list in host
  !
   implicit none
   !----   DUMMY Variables
           type(NEIGHBOR_LIST), dimension(:)::List
 
   !----   Local variables
             call Copyout_NeighboreList_A1(List, dm_Neighbors) 
         return
   end subroutine Copyout_NeighboreList_A2
  !**********************************************************************************

  !**********************************************************************************
  subroutine Copyout_NeighboreList_B0(hList, dList, GID)
  !***  PURPOSE:  to copy the neighbore list in device to a host neighbore list
  !               in old style. This routine can be used to test
  !               the neighbor calculation, it should be noted,
  !               in this GPU implement, the neighbor  list is
  !               the list with atoms are sorted by cells.
  !
  !
  !     INPUT:     dList,     the list on devices
  !                              
  !     OUTPUT     hList,     the list on host
  !
  implicit none
  !----   DUMMY Variables
          type(NEIGHBOR_LIST)              ::hList
          type(MDDEVNeighbor_list)         ::dList
          integer, dimension(:), optional  ::GID

  !----   Local variables
         integer::ERR, I, J, IA, IW, NP, NPA, STARTA, KVOIS

  !***
         NP = dList%NPRT
         !$$--- determine the array size of the neighbore list
         if(allocated(hList%KVOIS)) then
            if(size(hList%KVOIS) .lt. NP)  deallocate(hList%KVOIS)
         end if
         if(allocated(hList%INDI)) then
            if(size(hList%INDI) .lt. NP*dList%MXKVOIS) deallocate(hList%INDI)
         end if

         hList%MXKVOIS = dList%MXKVOIS
         if(.not. allocated(hList%KVOIS) ) then
            allocate(hList%KVOIS(NP) ,STAT=ERR)
            if(ERR) goto 100
         end if
         if(.not. allocated(hList%INDI) ) then
            allocate(hList%INDI(NP,m_mxKVOIS), STAT=ERR)
            if(ERR) goto 100
         end if                 

         !$$--- if swap array is nor allocated , allocate them
         do I=1, m_NDEVICE
            if(.not.allocated(m_SwapList(I)%KVOIS) ) then
                allocate(m_SwapList(I)%KVOIS(m_NAPDEV), STAT=ERR)
                if(ERR) goto 100
             end if    
             if(.not.allocated(m_SwapList(I)%INDI) ) then
                allocate(m_SwapList(I)%INDI(m_NAPDEV,m_mxKVOIS), STAT=ERR)
                if(ERR) goto 100
             end if
         end do

         do I=1, m_NDEVICE
            call DevCopyOut(m_SwapList(I)%KVOIS, dList%KVOIS(I), m_NAPDEV) 
            call DevCopyOut(m_SwapList(I)%INDI,  dList%INDI(I),  m_NAPDEV*dList%MXKVOIS )
         end do
         call SynchronizeDevices()

         do I=1, m_NDEVICE  
            STARTA = dList%STARTA(I)
            NPA    = dList%NPA(I)
            if(present(GID)) then
               do J=1, NPA
                  IA    = GID(J+STARTA-1)
                  KVOIS = m_SwapList(I)%KVOIS(J)
                  hList%KVOIS(IA) = KVOIS
                  do IW = 1, KVOIS
                     hList%INDI(IA,IW) = GID(m_SwapList(I)%INDI(J,IW))
                  end do
               end do
            else 
               do J=1, NPA
                  IA = J+STARTA-1
                  KVOIS = m_SwapList(I)%KVOIS(J)
                  hList%KVOIS(IA) = KVOIS
                  do IW = 1, KVOIS
                     hList%INDI(IA,IW) = m_SwapList(I)%INDI(J,IW)
                  end do
               end do   
             end if    
         end do
         return
 

   100  write(*,fmt="(A,  I7, A)")  " MDPSCU Error: fail to allocate memory for copyout neighbor list"
        write(*,fmt="(A)")          "               Process to be stopped"
        stop
        return
  end subroutine Copyout_NeighboreList_B0
  !**********************************************************************************  

  !**********************************************************************************
  subroutine Copyout_NeighboreList_B1(hList, dList)
  !***  PURPOSE:  to copy the neighbore list in device to a host neighbore list;
  !
  !     INPUT:     dList,     the list on devices
  !                              
  !     OUTPUT     hList,     the list on host
  !
  implicit none
  !----   DUMMY Variables
          type(NEIGHBOR_LIST)              ::hList
          type(MDDEVNeighbor_list)         ::dList

  !----   Local variables

  !***
          call  Copyout_NeighboreList_B0(hList, dList, hm_GID) 
          return
  end subroutine Copyout_NeighboreList_B1
  !**********************************************************************************  

  !**********************************************************************************
  subroutine Copyout_NeighboreList_B2(hList, Order)
  !***  PURPOSE:  to dm_Neighbors to a host neighbore list;
  !
  !     INPUT:     order,     the integer indicating if need to reorder the list
  !                              
  !     OUTPUT     hList,     the list on host
  !
   implicit none
   !----   DUMMY Variables
           type(NEIGHBOR_LIST)              ::hList
           integer                          ::Order 
   !----   Local variables
 
   !***
          if(Order) then
            call  Copyout_NeighboreList_B0(hList, dm_Neighbors, hm_GID) 
          else   
            call  Copyout_NeighboreList_B0(hList, dm_Neighbors) 
          end if   
   end subroutine Copyout_NeighboreList_B2
  !**********************************************************************************  

  !**********************************************************************************
  subroutine Copyout_NeighboreList_B3(hList)
  !***  PURPOSE:  to copy the dm_Neighbors to a host neighbore list;
  !
  !     INPUT:     
  !     OUTPUT     hList,     the list on host
  !
   implicit none
   !----   DUMMY Variables
           type(NEIGHBOR_LIST)              ::hList
   !----   Local variables
 
   !***
       call Copyout_NeighboreList_B2(hList, 0) 
       return 
   end subroutine Copyout_NeighboreList_B3
  !**********************************************************************************   
 
  !**********************************************************************************
   attributes(global) subroutine NeighboreList_IC_KERNEL(IM, IA0, XP, NPART, NAPB, GID, STATU, NCX, NCY, NCZ, NC, &
                                                         LBX, LBY, LBZ, BSX, BSY, BSZ, INC)
  !***  PURPOSE:  to determine the cell ID that the atoms inD
  !$$   INPUT:    IM,       the total number of particles (size of XP)
  !$$             IA0,      the first atom ID handled by the current device                                                        
  !$$             XP,       the position of the particles
  !$$             NPART,    the number of particles handled by the current device 
  !$$             NAPB,     the number of particles in each box
  !$$             GID,      the original ID of atoms of atoms after parition                                                         
  !$$             STATU,    the current statu of the particles
  !$$             NCX, NCY, NCZ, the cell number in each directions of a box
  !$$             NC = NCX*NCY*NCZ, the total number of cells in a box
  !$$             LBX, LBY, LBZ, the low bouandary of the box
  !$$             BSX, BSY, BSZ ,the size of the box
  !$$
  !$$   OUTPUT:   INC,      the index of cells that the particles in
  !
  implicit none
      !
      !--- DUMMY VARIABLES
      integer,      value::IM, NPART, NAPB, NCX, NCY, NCZ, NC, IA0
      real(KINDDF), value::LBX, LBY, LBZ, BSX, BSY, BSZ

      real(KINDDF), device :: XP(IM,*)
      integer,      device :: GID(*) 
      integer,      device :: STATU(*)
      integer,      device :: INC(*)

      !--- Local variables
      integer::IB, IT, IC, IC1, IX, IY, IZ, INB 
      real(KINDDF), parameter::eps=0.0001

      !---

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$ IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              !$$ IB -- now IB is the (id-1) of box the atom in
              IC1 = IC + IA0
              IB  = (GID(IC1) - 1)/NAPB
              if(IC .le. NPART)then
                 if(iand(STATU(IC),CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX) then
                    IX = int( (XP(IC1,1)-LBX)/BSX*dble(NCX) -eps)
                    IY = int( (XP(IC1,2)-LBY)/BSY*dble(NCY) -eps)
                    IZ = int( (XP(IC1,3)-LBZ)/BSZ*dble(NCZ) -eps)
                    if((IX .lt. 0 .or. IX .ge. NCX) .or. &
                       (IY .lt. 0 .or. IY .ge. NCY) .or. &
                       (IZ .lt. 0 .or. IZ .ge. NCZ)  ) then
                        INC(IC) = -2
                    else
                        INC(IC) = C_IUN + (IX+NCX*(IY+NCY*IZ)) + IB*NC
                    end if
                 else
                    INC(IC) = -1
                 end if
              end if

      return
  end subroutine NeighboreList_IC_KERNEL
  !**********************************************************************************

  !**********************************************************************************
  subroutine NeighboreList_IC_template0(IDEV, TNA, NAPB, NAPDEV, STARTA, XP, GID, STATU, IC)
  !***  PURPOSE:  to calculate the cell ID that the atoms in
  !
  !     INPUT:     IDEV,    the device ID
  !                TNA,     the total number of atoms 
  !                NAPB,    the number of atoms in each box
  !                NAPDEV,  the number of atoms handled by the device
  !                STARTA,  the index of the first atom
  !                XP,      the position of atoms
  !                GID,      the original ID of atoms of atoms after parition                                                            
  !                STATU,   the statu of atoms
  !
  !     OUTPUT     IC,       the cell ID that the atoms in
  !

  implicit none
  !----   DUMMY Variables
        integer,      intent(in)             ::IDEV, TNA, NAPB, NAPDEV, STARTA
        real(KINDDF), device, dimension(:,:) ::XP
        integer,      device, dimension(:)   ::GID
        integer,      device, dimension(:)   ::STATU
        integer,      device, dimension(:)   ::IC

  !----   Local variables

  !--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads
         integer::NB, BX, BY, ERR, IA0
  !***


         !$$--- first we create the cell id for atoms in device, this is performed
         !$$    on the major device 0.
             ERR = cudaSetDevice(IDEV)

             !$$-- to determine the dimension of blocks( the number of blocks in a grid)
             BX = mp_BLOCKSIZE
             BY = 1
             NB  = (NAPDEV-1)/(BX*BY)+1
             blocks  = dim3(NB, 1, 1)
             threads = dim3(BX, BY, 1)
             IA0     = STARTA - 1
             call NeighboreList_IC_KERNEL<<<blocks, threads>>>(TNA, IA0, XP, NAPDEV, NAPB, GID, STATU,            &
                                                               m_NCELL(1),   m_NCELL(2),   m_NCELL(3), m_NC0,     &
                                                               m_LBOX(1),    m_LBOX(2),    m_LBOX(3),             &   
                                                               m_BOXSIZE(1), m_BOXSIZE(2), m_BOXSIZE(3),          &
                                                               IC)

             return
  end subroutine NeighboreList_IC_template0
  !**********************************************************************************

  !**********************************************************************************
  subroutine NeighboreList_IC_template1(TNA, NAPB, NAPDEV, STARTA, XP, GID, STATU, IC)
  !***  PURPOSE:  to calculate the cell ID that the atoms in
  !
  !     INPUT:     XP,         the position of atoms
  !                STATU,      the statu of atoms
  !
  !     OUTPUT     IC,        the cell ID that the atoms in
  !

  implicit none
  !----   DUMMY Variables
        integer, intent(in)                ::TNA, NAPB
        integer, intent(in), dimension(:)  ::NAPDEV, STARTA
        type(DevMat_DF),     dimension(:)  ::XP
        type(DevVec_I),      dimension(:)  ::GID
        type(DevVec_I),      dimension(:)  ::STATU
        type(DevVec_I),      dimension(:)  ::IC
 
  !----   Local variables
         integer::I, IDEV

          do I=1, m_NDEVICE
             call NeighboreList_IC_template0(IC(I)%IDEV, TNA, NAPB, NAPDEV(I), STARTA(I),  &
                                          XP(I)%Data, GID(I)%Data, STATU(I)%Data, IC(I)%Data)     
          end do  
       return
  end subroutine NeighboreList_IC_template1
  !**********************************************************************************  

  !**********************************************************************************
   attributes(global) subroutine Cal_NeighboreList_Kernel2C(NBPC, NPART, XP, ITYP, NC, NAC, NAAC, IA1th0, IA1th, NCX, NCY, NCZ,  &
                                                            BSX, BSY, BSZ, PDX, PDY, PDZ, CFROM, CTO, mxNAPDEV,                  &
                                                            mxKVOIS, KVOIS, INDI)
  !***  PURPOSE:  to update the neighbore list of particles in an cell.
  !               Newton's third law is NOT considered . Be careful to use
  !               corresponding force calculation where Newton's 3th should NOT be applied
  !
  !$$   INPUT:    NBPC,     the number of blocks needed for covering all particles in a cell
  !$$             NPART,    the total number of particles in the whole simulation box
  !$$             XP,       the position of the particles
  !$$             ITYP,     the type of atoms, used when the cutoff distance between different type of particles are different
  !$$             NC,       the number of cells in the whole box
  !$$             NAC,      the number of particles in the cells on the device
  !$$             NAAC,     the number of ACTIVE particles in the cells                                                         
  !$$             IA1th0,   the index of the first particle on the device
  !$$             IA1th,    the index of the first particle in a cell
  !$$             NCX, NCY, NCZ, the cell number in X, Y, Z direction
  !$$             BSX, BSY, BSZ, the size of the box in X, Y, Z, dierection
  !$$             PDX, PDY,PDZ,  the perodic condition in X, Y, Z, dierection
  !$$             CFROM,     the start cell on the device
  !$$             CTO,       the end cell on the device
  !$$             mxNAPDEV,  the max number of particles on a device, e.g, the dimension of KVOIS
  !$$             mxKVOIS,   the max permitted number of neighbors, e.g, the second dimension of INDI
  !$$
  !$$     OUTPUT: KVOIS,    the number of neighbores of particles
  !$$             INDI,     the index of neighbore particles
  !
  implicit none
      !
      !--- DUMMY VARIABLES
      integer,      value::NBPC,IP0, NC,NCX,NCY,NCZ, PDX, PDY, PDZ, NPART, CFROM, CTO, mxNAPDEV, mxKVOIS, IA1th0
      real(KINDDF), value::BSX, BSY,BSZ
      real(KINDDF), device::XP(NPART,*)
      integer,      device::ITYP(NPART), NAC(NC), NAAC(NC), IA1th(NC)

      integer,      device::KVOIS(mxNAPDEV)
      integer,      device::INDI(mxNAPDEV,*)

      !--- Local variables
         !nonshared by threads
         real(KINDSF)::POS1, POS2, POS3,  DXYZ1, DXYZ2, DXYZ3, SEP1, SEP2, SEP3
         integer::IB, IB0, IT, IA, IA0,IA00, JA, NN, I, J, K, ITY

         !variables share by all thread
         integer, shared::NB, IC, IS0, STARTA, NCXY, NCXYZ,IX0, IY0, IZ0,IC0, NACC0
         integer, shared::NS,NACC, IAC, IACE, FROM, TO
         integer, shared, dimension(mp_NNC)::CID, IX,IY, IZ, OUT

         real(KINDSF), shared, dimension(mp_NNC)::CXYZ1, CXYZ2, CXYZ3
         real(KINDSF), shared, dimension(mp_BLOCKSIZE)::SPOS1, SPOS2, SPOS3
         integer,      shared, dimension(mp_BLOCKSIZE)::JTY


      !$$--- start process
              !$$IB -- the index of block
              !$$
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x-1

              !$$NB -- the size of a block
              NB = blockdim%x*blockdim%y

              !$$IB0 -- the index of the cells, if the blocksize is smaller than the
              !$$       number of atoms in a cell (NBPC=1), each block cover all atoms
              !$$       in a cell, otherwise, use NBPC blocks to cover the atoms.
              IB0  =  IB/NBPC

              !$$--- beofre move to the globla index of the cell
              !$$    we get the starting atom in the block
              IP0 = (IB-IB0*NBPC)*NB

              !$$--- move IB0 to the globla index of cell
              IB0  =  IB0 + CFROM-1
              if(IB0 .ge. CTO) return

              !$$--- if no active paritcl in the cell, we return
              if(NAAC(IB0+1) .le. 0) return

              !$$IT -- the thread ID
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              !$$-- calculate the index of the cell in original box, performed by the first thread
              if(IT .eq. 1) then
                 NCXY = NCX*NCY
                 NCXYZ = NCXY*NCZ
                 !$$--- IS0, the ID of a box in a multiple-boxs style
                 IS0 = IB0/NCXYZ
                 IC  = IB0-IS0*NCXYZ

                 !$$-- 3D index of the original cell in original box
                 IZ0 = IC/NCXY
                 IY0 = (IC-IZ0*NCXY)/NCX
                 IX0 = IC-IZ0*NCXY-IY0*NCX
                 IZ0 = IZ0 + 1
                 IY0 = IY0 + 1
                 IX0 = IX0 + 1

                 !$$IC-- now IC is the index of the cell in the whole box
                 IC  = IB0 + 1

                 !$$STARTA -- the globle index of atom in cell CFROM
                 !$$          it is also the first atom on the device
                 STARTA =  IA1th0 !IA1th(CFROM)

                 !$$NACC0 -- the number of atoms in the first cell
                 NACC0 = NAC(IC)

              end if
              !$$--- to insure all thread has readed here
              call syncthreads()

              !$$--- if the cell is empty, return
              if(NACC0 .le. 0) return

              !$$--- to get the position of neighbore cells
              !$$      of the original cell
               if(IT .le. mp_NNC) then
                 OUT(IT)= 0
                 IZ(IT) = IZ0 + dcm_NIZ(IT)
                 IY(IT) = IY0 + dcm_NIY(IT)
                 IX(IT) = IX0 + dcm_NIX(IT)

                 CXYZ1(IT) = C_ZERO
                 CXYZ2(IT) = C_ZERO
                 CXYZ3(IT) = C_ZERO
                 if(PDX .and. IT.gt.1) then
                     if( IX(IT).gt.NCX )then
                         IX(IT) = C_IUN
                         CXYZ1(IT) = BSX
                     else if (IX(IT).LT.C_IUN) then
                         IX(IT) = NCX
                         CXYZ1(IT) = -BSX
                     end if
                  end if

                  if(PDY .and. IT.gt.1) then
                     if( IY(IT).gt.NCY )then
                         IY(IT) = C_IUN
                         CXYZ2(IT) = BSY
                     else if (IY(IT).LT.C_IUN) then
                         IY(IT) = NCY
                         CXYZ2(IT) = -BSY
                     end if
                  end if

                  if(PDZ .and. IT.gt.1) then
                     if( IZ(IT).gt.NCZ )then
                         IZ(IT) = C_IUN
                         CXYZ3(IT) = BSZ
                     else if (IZ(IT).LT.C_IUN) then
                         IZ(IT) = NCZ
                         CXYZ3(IT) = -BSZ
                     end if
                  end if

                  if( IX(IT) .gt. NCX .OR. IX(IT) .LT. C_IUN) OUT(IT) = 1
                  if( IY(IT) .gt. NCY .OR. IY(IT) .LT. C_IUN) OUT(IT) = 1
                  if( IZ(IT) .gt. NCZ .OR. IZ(IT) .LT. C_IUN) OUT(IT) = 1
                  if(OUT(IT) .EQ. 0) then
                     CID(IT) = NCXY*(IZ(IT)-1)+NCX*(IY(IT)-1)+IX(IT)+IS0*NCXYZ
                  else
                     CID(IT) = 0
                  end if
               end if
               call syncthreads()
               !$$--- we start scanning for atom to calculate their neighbors
               !$$--- IA:   the global index of the atom in the original cell.
               !$$--- IA0:  the index of the atom on the DEVICE
               !$$--- IA00: the index of the atom in its cell
               !$$
               IA    = (IT-1) + IA1th(IC) + IP0
               IA0   = IA - STARTA   + 1
               IA00  = IA - IA1th(IC)+ 1
               ITY   = ITYP(IA)
               POS1  = XP(IA, 1)
               POS2  = XP(IA, 2)
               POS3  = XP(IA, 3)
               NN    = 0

               !$$--- serach in cell the atom in
               K=1
                  !$$NACC -- the number of atoms in neighboring cell K
                  NACC = NAC(CID(K))
                  !$$NS-- the number of segment of do loop when scanning cell K
                  NS = (NACC-1)/NB+1
                  !$$IAC-- the index of start atom in cell K
                  IAC = IA1th(CID(K))

                  !$$IACE-- the index of end atom in cell K
                  IACE = IAC + NACC -1

                  do J=1, NS
                     FROM      = min((J-1)*NB+IAC,IACE)
                     TO        = min(FROM+NB-1, IACE)
                     SPOS1(IT) = XP(IT+FROM-1, 1) + CXYZ1(K)
                     SPOS2(IT) = XP(IT+FROM-1, 2) + CXYZ2(K)
                     SPOS3(IT) = XP(IT+FROM-1, 3) + CXYZ3(K)
                     JTY(IT)   = ITYP(IT+FROM-1)
                     call syncthreads()
                     !$$--- In each block we calculate the neigbores of atoms in an original cell
                     if(IA00.le.NACC0) then
                        do I=FROM, TO
                              JA = I-FROM+1
                              SEP1 = POS1 - SPOS1(JA)
                              SEP2 = POS2 - SPOS2(JA)
                              SEP3 = POS3 - SPOS3(JA)

                              !DXYZ(1) = sum(dcm_BOXSHAPE(1,1:3)*(POS(1:3) - SPOS(JA,1:3)))
                              !DXYZ(2) = sum(dcm_BOXSHAPE(2,1:3)*(POS(1:3) - SPOS(JA,1:3)))
                              !DXYZ(3) = sum(dcm_BOXSHAPE(3,1:3)*(POS(1:3) - SPOS(JA,1:3)))

                              DXYZ1 = dcm_BOXSHAPE(1,1)*SEP1 + dcm_BOXSHAPE(1,2)*SEP2 + dcm_BOXSHAPE(1,3)*SEP3
                              DXYZ2 = dcm_BOXSHAPE(2,1)*SEP1 + dcm_BOXSHAPE(2,2)*SEP2 + dcm_BOXSHAPE(2,3)*SEP3
                              DXYZ3 = dcm_BOXSHAPE(3,1)*SEP1 + dcm_BOXSHAPE(3,2)*SEP2 + dcm_BOXSHAPE(3,3)*SEP3

                             !$$-- NOTE: periodic condition has been considered when we create cells,
                             !$$         here we do not need to check it.
                              if( DXYZ1*DXYZ1+DXYZ2*DXYZ2+DXYZ3*DXYZ3 .le. dcm_RM2(ITY,JTY(JA)) ) then
                                  if(I.ne.IA) then
                                     NN  = NN + 1
                                     if(NN.le.mxKVOIS) then
                                        !$$---NOTE:  IA0 is the index of atoms on the current device
                                        !$$          I is the globle index of atom in the whole box
                                        INDI(IA0, NN) = I
                                     end if
                                 end if
                             end if
                        end do
                     end if
                     call syncthreads()
                  end do

               !$$--- search in neighbor cells
               do K=2, mp_NNC
                  if(OUT(K)) cycle
                  !$$NACC -- the number of atoms in neighboring cell K
                  NACC = NAC(CID(K))

                  !$$NS-- the number of segment of do loop when scanning cell K
                  NS = min((NACC-1)/NB+1, NACC)

                  !$$IAC-- the index of start atom in cell K
                  IAC = IA1th(CID(K))

                  !$$IACE-- the index of end atom in cell K
                  IACE = IAC + NACC -1
                  call syncthreads()

                  do J=1, NS
                     FROM       = min((J-1)*NB+IAC,IACE)
                     TO         = min(FROM+NB-1, IACE)
                     SPOS1(IT) = XP(IT+FROM-1, 1) + CXYZ1(K)
                     SPOS2(IT) = XP(IT+FROM-1, 2) + CXYZ2(K)
                     SPOS3(IT) = XP(IT+FROM-1, 3) + CXYZ3(K)
                     JTY(IT)    = ITYP(IT+FROM-1)
                     call syncthreads()
                     !$$--- In each block we calculate the neigbores of atoms in an original cell
                     if(IA00.le.NACC0) then
                        do I=FROM, TO
                           JA = I-FROM+1
                           SEP1 = POS1 - SPOS1(JA)
                           SEP2 = POS2 - SPOS2(JA)
                           SEP3 = POS3 - SPOS3(JA)

                           !DXYZ(1) = sum(dcm_BOXSHAPE(1,1:3)*(POS(1:3) - SPOS(JA,1:3)))
                           !DXYZ(2) = sum(dcm_BOXSHAPE(2,1:3)*(POS(1:3) - SPOS(JA,1:3)))
                           !DXYZ(3) = sum(dcm_BOXSHAPE(3,1:3)*(POS(1:3) - SPOS(JA,1:3)))
                           DXYZ1 = dcm_BOXSHAPE(1,1)*SEP1 + dcm_BOXSHAPE(1,2)*SEP2 + dcm_BOXSHAPE(1,3)*SEP3
                           DXYZ2 = dcm_BOXSHAPE(2,1)*SEP1 + dcm_BOXSHAPE(2,2)*SEP2 + dcm_BOXSHAPE(2,3)*SEP3
                           DXYZ3 = dcm_BOXSHAPE(3,1)*SEP1 + dcm_BOXSHAPE(3,2)*SEP2 + dcm_BOXSHAPE(3,3)*SEP3

                          !$$-- NOTE: periodic condition has been considered when we create cells,
                          !$$         here we do not need to check it.
                           if( DXYZ1*DXYZ1+DXYZ2*DXYZ2+DXYZ3*DXYZ3 .le. dcm_RM2(ITY,JTY(JA)) ) then
                                NN  = NN + 1
                                if(NN.le.mxKVOIS) then
                                   !$$---NOTE:  IA0 is the index of atoms on the current device
                                   !$$          I is the globle index of atom in the whole box
                                    INDI(IA0, NN) = I
                                end if
                           end if
                        end do
                     end if
                     call syncthreads()
                  end do

                end do
             !--------------
             if(IA00 .le. NACC0) then
                KVOIS(IA0) = min(NN, mxKVOIS)
             end if

      return
  end subroutine Cal_NeighboreList_Kernel2C
  !**********************************************************************************

  !**********************************************************************************
  subroutine StartOnDevice_template(IDEV, CFROM, CTO, XP, ITYP, NAC, NAAC, IA1th, KVOIS, INDI)
  !***  PURPOSE:  to start neighboring calculation on a device
  !
  !               Newton's 3rd law is NOT considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should NOT be applied
  !
  !    INPUT(CPU):  IDEV,      the index of the device
  !                 CFROM,     the start cell on the device
  !                 CTO,       the end cell on the device
  !    INPUT(GUP):
  !                XP,        the position of the atoms
  !                ITYP,      the type of atoms, used when the cutoff distance between different kinds of atoms are different
  !                NAC,       the number of atoms in the cells
  !                NAAC,      the number of ACTIVE particles in the cells 
  !                IA1th,     the index of the first atom in a cell
  !
  !   OUTPUT(GPU): KVOIS,     the number of neighbores of atoms
  !                INDI,      the index of neighbore atoms
  !
  implicit none
  !----   DUMMY Variables
          integer::IDEV,CFROM, CTO

          integer,      device, dimension(:)   ::NAC,NAAC,IA1th
          real(KINDDF), device, dimension(:,:) ::XP
          integer,      device, dimension(:)   ::ITYP
          integer,      device, dimension(:)   ::KVOIS
          integer,      device, dimension(:,:) ::INDI

  !----   Local variables
         type(dim3) :: blocks
         type(dim3) :: threads
         integer::NS, BX, BY, IP0, ERR, NADEV, CURDEV, C0, C1, NCINB, RESC, STARTA
  !--- start

               !$$--- copy the cell informations into device memeory
               ERR = cudaSetDevice(IDEV)

               STARTA = hm_IA1th(CFROM)
               !$$--- to determine of a block (the number of threads in a block)
               !$$     mp_BLOCKSIZE must be larger than 27
               BX = min(mp_BLOCKSIZE, m_MXNAC0)
               BX = max(BX, mp_NNC)
               BY = 1
               threads = dim3(BX, BY, 1)
               NS = (m_MXNAC0-1)/(BX*BY)+1

               !$$-- to determine the dimension of blocks( the number of blocks in a grid)
               !$$   note the max gridsize is 65535*65535*65535, in C2050
               RESC  = CTO - CFROM + 1
               NCINB = min(RESC, 65535)
               C0    = CFROM
               C1    = C0 + NCINB - 1
               do while(.true.)
                  if(C0 .le. 0) exit
                  blocks  = dim3(NCINB, NS, 1)
                  call Cal_NeighboreList_Kernel2C<<<blocks,threads>>>(NS, dm_NPRT,XP, ITYP, m_NC, NAC,NAAC, STARTA, IA1th, &
                                                                      m_NCELL(1), m_NCELL(2), m_NCELL(3),                  &
                                                                      m_BOXSIZE(1), m_BOXSIZE(2), m_BOXSIZE(3),            &
                                                                      m_PD(1), m_PD(2), m_PD(3), C0, C1,                   &
                                                                      m_NAPDEV, m_mxKVOIS, KVOIS, INDI)

                  RESC  = RESC - NCINB
                  NCINB = min(RESC, 65535)
                  if(NCINB .le. 0) exit
                  C0 = C1 + 1
                  C1 = C0 + NCINB - 1
                end do
              
       return
   end subroutine StartOnDevice_template
  !****************************************************************************

  !**********************************************************************************
  subroutine Copy_From_Host_To_Devices_template(IDEV, CPYSTREAM, dXP1, dXP2, dXP3, dXP4, dXP5, dDIS, dSTATU)
  !***  PURPOSE:  to copy the configuration of particles, that have been clustered by cells,
  !               from host to a device.
  !
  !               NOTE: since ITYP and XP have been copied in, here we do not need to copyin
  !                     these arraies. FP is an instant quantity, here we do not copyin
  !                     either. However, dXP1, dXP2, dXP3, dXP4, dXP5, dDIS, dSTATU are
  !                     accumulation quantities, thus we need to copy them in.
  !
  !     INPUT:     IDEV,      the ID of device
  !                CPYSTREAM, the stream ID for the copy operation
  !
  !     OUTPUT     dXP1, dXP2, dXP3, dXP4, dXP5, dDIS, dFP, dSTATU:
  !                the corresponding vairbales on devices
  !
  !
      implicit none
      !--- dummy variables
           integer,      intent(in)            ::IDEV,CPYSTREAM
           integer,      device, dimension(:)  ::dSTATU
           real(KINDDF), device, dimension(:,:)::dXP1, dXP2, dXP3, dXP4, dXP5, dDIS

      !--- Local vairables
           integer::STARTA, ENDA, NA, ERR

              !$$--- get the current device

                  !$$--- the first atom on the device
                  STARTA = hm_IA1th(m_STARTCELL(IDEV))

                  !$$--- the last atom on the device
                  ENDA = hm_IA1th(m_ENDCELL(IDEV))+hm_NAC(m_ENDCELL(IDEV))-1

                  !$$--- the number of atoms on the device
                  NA = ENDA - STARTA + 1
                  !$$--- NOTE the size dXP1 etc, are different from that of hm_XP1...
                   ERR = cudaMemcpyAsync(dXP1(1,1), hm_XP1(STARTA,1), NA) !, stream=cpystream)
                   ERR = cudaMemcpyAsync(dXP1(1,2), hm_XP1(STARTA,2), NA) !, stream=cpystream)
                   ERR = cudaMemcpyAsync(dXP1(1,3), hm_XP1(STARTA,3), NA) !, stream=cpystream)
                   ERR = cudaMemcpyAsync(dDIS(1,1), hm_DIS(STARTA,1), NA) !, stream=cpystream)
                   ERR = cudaMemcpyAsync(dDIS(1,2), hm_DIS(STARTA,2), NA) !, stream=cpystream)
                   ERR = cudaMemcpyAsync(dDIS(1,3), hm_DIS(STARTA,3), NA) !, stream=cpystream)

                   ERR = cudaMemcpyAsync(dSTATU(1), hm_STATU(STARTA), NA) !, stream=cpystream)
                  if(allocated(hm_XP2)) then
                     ERR = cudaMemcpyAsync(dXP2(1,1), hm_XP2(STARTA,1), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP2(1,2), hm_XP2(STARTA,2), NA) ! , stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP2(1,3), hm_XP2(STARTA,3), NA) !, stream=cpystream)

                     ERR = cudaMemcpyAsync(dXP3(1,1), hm_XP3(STARTA,1), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP3(1,2), hm_XP3(STARTA,2), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP3(1,3), hm_XP3(STARTA,3), NA) !, stream=cpystream)

                     ERR = cudaMemcpyAsync(dXP4(1,1), hm_XP4(STARTA,1), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP4(1,2), hm_XP4(STARTA,2), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP4(1,3), hm_XP4(STARTA,3), NA) !, stream=cpystream)

                     ERR = cudaMemcpyAsync(dXP5(1,1), hm_XP5(STARTA,1), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP5(1,2), hm_XP5(STARTA,2), NA) !, stream=cpystream)
                     ERR = cudaMemcpyAsync(dXP5(1,3), hm_XP5(STARTA,3), NA) !, stream=cpystream)

                  end if

                 return

  end subroutine Copy_From_Host_To_Devices_template
  !****************************************************************************

  !**********************************************************************************
  subroutine Cal_NeighBoreList2C_0_DEV(SimBox, CtrlParam, MM)
  !***  PURPOSE:  to update the neighbore list of particles with head-link
  !               technique used.
  !
  !               Newton's 3rd law is NOT considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should NOT be applied
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                MM,        number of boxes
  !
  !     OUTPUT     simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          type(SimMDBox)    ::SimBox
          type(SimMDCtrl)   ::CtrlParam
          integer           ::MM

  !----   Local variables
         integer I, J, K, M, ID, IP, NPART, NAPB, IX, IY, IZ, IC, NN, BCFLG

  !--- Device variables and variables to be used in GPU
         integer::ERR, CURDEV,  NAPDEV
         integer, parameter::UN(m_MXDEVICE)=1
         real(KINDSF)::BOXSHAPE(3,3) = 0.0
         real(KINDSF)::RM2(mp_MXGROUP,mp_MXGROUP) = 0.0

         character*1::tinput=""
         save::tinput, BOXSHAPE, RM2
         real(KINDDF), parameter::EPS=(0.01*CP_A2CM)*(0.01*CP_A2CM)
  !***
         !$$--- set the constant memories on devices
         m_LBOX     = SimBox%BOXLOW
         m_BOXSIZE  = SimBox%ZL
         m_PD       = CtrlParam%IFPD
         m_RM2      = CtrlParam%NB_RM*CtrlParam%NB_RM
         m_BOXSHAPE = SimBox%BOXSHAPE
         BCFLG      = 0
         if(any(dabs(m_BOXSHAPE - BOXSHAPE) .gt. eps) .or. &
            any(dabs(m_RM2 - RM2)           .gt. eps)     ) then
            BCFLG    = 1
            BOXSHAPE = m_BOXSHAPE
            RM2      = m_RM2
         end if      

         NAPB  = SimBox%NPRT
         NPART = NAPB*MM
         if(NPART .ne. dm_NPRT) then
            write(*,*) "MDPSCU Error: the particle number in neighbor calculation is different from that in device0", NPART, dm_NPRT
            write(*,*) "Process to be stopped"
            stop
         end if

         !$$--- keep the current device and create streeam to be used for host-device
         !$$    data transfer
         ERR = cudaGetDevice(CURDEV)
         do I=1, m_NDEVICE
            ERR = cudaSetDevice(m_DEVICES(I))
            !ERR  = cudaStreamCreate(CPYSTRDEVS(I))
         !--- if the box shape have been changed, we reset the shape
            !    and cuttof radius
            if(BCFLG .gt. 0) then
               dcm_BOXSHAPE = BOXSHAPE
               dcm_RM2      = RM2
            end if   
         end do
         
         !$$--- If is not the first time of partitioning the system
         !$$    the configuration should be updated
         if(hm_HasPart .ge. mp_HASPART) then
            !$$--- NOTE: 2019-07-03, when there are atom out of box, these atoms are placed at the 
            !$$          the end of paritioned XP, thus GID are needed to identify the box ID that 
            !$$          the atoms belong to. If the box ID was determined by the index of XP, erro would occur  
            call NeighboreList_IC_template1(NPART, NAPB, m_NPA, m_STARTA, dm_WorkSpace%XP, dm_WorkSpace%GID,  &
                                            dm_WorkSpace%STATU, dm_WorkSpace%IC)     
            call DevCopyOut(hm_INC,    m_STARTA, dm_WorkSpace%IC,    UN,       m_NPA)
            call DevCopyOut(hm_XP,     m_STARTA, dm_WorkSpace%XP,    m_STARTA, m_NPA)
            call DevCopyOut(hm_STATU,  m_STARTA, dm_WorkSpace%STATU, UN,       m_NPA)
            !$$--- to initial the LINKED-CELL
            hm_NAC   = 0
            hm_NAAC  = 0
            m_HEAD   = 0
            
            !--- restore the orignal order 
            call SynchronizeDevices()
            m_XP   ( hm_GID(1:NPART), :)  =  hm_XP    (1:NPART, :)
            m_STATU( hm_GID(1:NPART))     =  hm_STATU (1:NPART)
            m_INC  ( hm_GID(1:NPART))     =  hm_INC(1:NPART)
         !$$--- first we create the cell id for atoms in device, this is performed
         !$$    on the major device 0.
         else if(hm_HasPart .eq. mp_NOTPART) then
                !$$--- for the first time of partition, we assign each GPU with approximate equal
                !$$    atoms to handle
                NAPDEV      = (dble(NPART+m_NDEVICE)/dble(m_NDEVICE) + 0.00001)
                m_STARTA(1) = 1
                m_NPA(1)    =  min(NAPDEV, NPART)
                NN = m_STARTA(1) + m_NPA(1) - 1
                do I=2, m_NDEVICE
                   m_STARTA(I) = NN + 1
                   m_NPA(I)    = min(NAPDEV, NPART-NN)
                   NN          = m_STARTA(I) + m_NPA(I) - 1
                end do 
                
                do I=1, NPART
                   hm_GID(I) = I 
                end do   
                do I=1, m_NDEVICE
                   call DevCopyIn(m_XP,                  dm_WorkSpace%XP(I),    3*NPART)
                   call DevCopyIn(m_STATU(m_STARTA(I):), dm_WorkSpace%STATU(I), m_NPA(I))
                   call DevCopyIn(hm_GID,                dm_WorkSpace%GID(I),     NPART)
                end do 
                call NeighboreList_IC_template1(NPART, NAPB, m_NPA, m_STARTA, dm_WorkSpace%XP, dm_WorkSpace%GID,  &
                                               dm_WorkSpace%STATU, dm_WorkSpace%IC)     

                call DevCopyOut(m_INC, m_STARTA, dm_WorkSpace%IC, UN,  m_NPA)
                !$$--- to initial the LINKED-CELL
                hm_NAC   = 0
                hm_NAAC  = 0
                m_HEAD   = 0
                call SynchronizeDevices()
         end if

             !$$--- During partioning the system on host, we
             !$$    issue the asynchronized copy of arraies beyond XP and STATU
             if(hm_HasPart .ge. mp_HASPART) then
                call  DevCopyOut(hm_XP1, m_STARTA, dm_WorkSpace%XP1, UN, m_NPA)
                call  DevCopyOut(hm_DIS, m_STARTA, dm_WorkSpace%DIS, UN, m_NPA)
                if(allocated(hm_XP2)) then
                  call  DevCopyOut(hm_XP2, m_STARTA, dm_WorkSpace%XP2, UN, m_NPA)
                  call  DevCopyOut(hm_XP3, m_STARTA, dm_WorkSpace%XP3, UN, m_NPA)
                  call  DevCopyOut(hm_XP4, m_STARTA, dm_WorkSpace%XP4, UN, m_NPA)
                  call  DevCopyOut(hm_XP5, m_STARTA, dm_WorkSpace%XP5, UN, m_NPA)
                end if
             end if   

             !$$--- to do pratition on the host
             !$$--- to create the LINKED-CELL
             m_NUMOUT = 0
             do I=1, NPART
                if(m_INC(I) .gt. 0) then
                   IC = m_INC(I)
                   J  = m_HEAD(IC)
                   m_HEAD(IC) = I
                   m_LINK(I)  = J
                   hm_NAC(IC) = hm_NAC(IC)+1
                   if(iand(m_STATU(I),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                     !--- NOTE: hm_NAC contain atoms, acitve or inactive,  in a cell 
                     !          hm_NAAC contain only ACTIVE atoms 
                     hm_NAAC(IC) = hm_NAAC(IC) + 1
                   end if    
                else
                   !$$--- if a particle is marked as active and outof box
                   !$$    m_INC = -2, and an error message should issued
                   if(m_INC(I) .lt. -1) then
                      if(tinput(1:1) .ne. 'C') then
                         write(*,fmt="(A, I8, A)") " MDPSCU Error: particle ",I," marked as active but out of box"
                         write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "               box range x ", &
                                                                        m_LBOX(1), m_LBOX(1)+m_BOXSIZE(1), ", particle x ", m_XP(I, 1)
                         write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "               box range y ", &
                                                                        m_LBOX(2), m_LBOX(2)+m_BOXSIZE(2), ", particle y ", m_XP(I, 2)
                         write(*,fmt="(A, 1PE12.5,1PE12.5, A, 1PE12.5)") "               box range z ", &
                                                                        m_LBOX(3), m_LBOX(3)+m_BOXSIZE(3), ", particle z ", m_XP(I, 3)

                         write(*,fmt="(A)")"You can input 'C' to marked all out-of-box atoms as inactive and continue the process."
                         write(*,fmt="(A)")"Or input Enter to stop the process:"
                         read(*, fmt="(A1)") tinput
                         if(tinput(1:1) .eq. 'c') tinput(1:1) = 'C'
                         if(tinput(1:1) .ne. 'C') then
                            write(*,fmt="(A)")"Process to be stopped"
                            stop
                         end if
                      end if
                      m_STATU(I) = CP_STATU_OUTOFBOX
                   end if
                   m_NUMOUT = m_NUMOUT + 1
                end if
             end do
             if(m_NUMOUT .GT. 0) then
                write(*,fmt="(A, I8, A)") " MDPSCU Warning: there are ",m_NUMOUT," atoms out of box found in MD_NeighborsList."
                call ONWARNING(gm_OnWarning)
             end if
             !---- 
             call DevCopyIn(hm_NAC,  dm_WorkSpace%NAC,  m_NC)
             call DevCopyIn(hm_NAAC, dm_WorkSpace%NAAC, m_NC)
             m_MXNAC0 = maxval(hm_NAC)

           !*******************************************************************************
           !$$--- to determine the global index of atoms in cells
           !$$    We need ITYP and XP intermediately, thus we should have them ready
            IP = 0
            hm_IA1th = 1
            do M=1, MM
               do IZ=1, m_NCELL(3)
                  do IY=1, m_NCELL(2)
                     do IX=1, m_NCELL(1)
                        !$$--- for the first, we get the ID of the original cell
                        IC = IX+m_NCELL(1)*((IY-1)+m_NCELL(2)*(IZ-1))+(M-1)*m_NC0
                        if(IC.GT.C_IUN) then
                           hm_IA1th(IC) = hm_IA1th(IC-1)+hm_NAC(IC-1)
                        end if

                        ID = m_HEAD(IC)
                        if(ID.EQ. C_IZERO) then
                           cycle !$$ no particle in the cell
                        end if

                        do while(.TRUE.)
                           IP = IP + 1
                           hm_GID(IP)    = ID
                           hm_ITYP(IP)   = m_ITYP(ID)
                           hm_XP(IP,1:3) = m_XP(ID,1:3)
                           ID = m_LINK(ID)
                           if(ID.LE.C_IZERO) exit
                        end do

                     end do
                  end do
               end do
            end do
            !---
            call DevCopyIn(hm_ITYP,  dm_WorkSpace%ITYP,  NPART)
            call DevCopyIn(hm_XP,    dm_WorkSpace%XP,    NPART*3)
            call DevCopyIn(hm_IA1th, dm_WorkSpace%IA1th, m_NC)

           !$$--- to assign the cells that the devices will handle.
           !$$    we assigne each device with same atoms as possible 
            NAPDEV = (dble(NPART)/dble(m_NDEVICE))
            m_STARTCELL(1)  = 1
            K = 1
            do I=1, m_NDEVICE-1
               NN = 0
               do J=K, m_NC
                  NN = NN + hm_NAC(J)
                  if(NN.GE.NAPDEV .and. NN.LE.m_NAPDEV) then
                     m_ENDCELL(I) = J
                     K = min(m_ENDCELL(I)+1,m_NC)
                     m_STARTCELL(I+1) = K
                     exit
                  end if
               end do
            end do
            m_ENDCELL(m_NDEVICE) = m_NC

           !$$--- check if the number of atoms on a device larger than
           !$$    permitted value
            do I=1, m_NDEVICE
               m_STARTA(I) = hm_IA1th(m_STARTCELL(I))
               m_ENDA(I)   = hm_IA1th(m_ENDCELL(I)) + hm_NAC(m_ENDCELL(I))-1
               m_NPA(I)    = m_ENDA(I) - m_STARTA(I) + 1  
               if(m_NPA(I) .GT. m_NAPDEV) then
                  write(*,*) "MDPSCU Error: the number of particles on each device is larger than evaluated value", m_NAPDEV
                  stop
               end if
            end do

            !$$--- to make sure the copyout of arries beyond XP and STATUS
            !      have finished
            call SynchronizeDevices()
            !*******************************************************************************
            !$$--- to start neighbor calculation on devices,
            !$$    NOTE: XP, ITYP, NAC, IA1th on devices are to be copyin duing 
            !$$          the calculation
            do I=1, m_NDEVICE
               call StartOnDevice_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),  &
                                           dm_WorkSpace%XP(I)%Data,                     &
                                           dm_WorkSpace%ITYP(I)%Data,                   &
                                           dm_WorkSpace%NAC(I)%Data,                    &
                                           dm_WorkSpace%NAAC(I)%Data,                   &
                                           dm_WorkSpace%IA1th(I)%Data,                  &
                                           dm_Neighbors%KVOIS(I)%Data,                  &
                                           dm_Neighbors%INDI(I)%Data   )
            end do 
             !$$--- if there are atoms are out-of-box, we need sort these atoms
             !$$    backward
             !$$    NOTE: hm_XP for atom-box have been ordered
             if(m_NUMOUT .gt. 0) then
                IP = NPART+1
                do I=1, NPART
                   if(m_INC(I) .lt. 0) then
                      IP          = IP -1
                      hm_GID(IP)  = I
                      hm_ITYP(IP) = m_ITYP(I)
                      hm_XP(IP,:) = m_XP(I,:)
                   end if
                end do
             end if
             !$$---  we copyin the cell id for the atoms located in,
             !$$     that could be used in acivation region calculation 
             !else
             hm_INC(1:NPART) = m_INC(hm_GID(1:NPART))
             call DevCopyIn(hm_INC, m_STARTA, dm_WorkSpace%IC, UN, m_NPA)

             !$$--- During calculating neighbolist, we need update the arries 
             !$$    beyond XP and STATU on host
             !$$    NOTE: at this moment, hm_GID has been changed, 
             !$$          but hm_GIDINV is still no changed
             !$$    NOTE: m_STATU has been update bfore calculating INC. But at this moment
             !$$          hm_STATU s still not updated
             !$$ 
             if(hm_HasPart .ge. mp_HASPART) then
                m_XP1(1:NPART, :) = hm_XP1(hm_GIDINV(1:NPART), :)
                m_DIS(1:NPART, :) = hm_DIS(hm_GIDINV(1:NPART), :)
                if(allocated(hm_XP2) .and. allocated(m_XP2)) then
                   m_XP2(1:NPART, :) =  hm_XP2(hm_GIDINV(1:NPART), :)
                   m_XP3(1:NPART, :) =  hm_XP3(hm_GIDINV(1:NPART), :)
                   m_XP4(1:NPART, :) =  hm_XP4(hm_GIDINV(1:NPART), :)
                   m_XP5(1:NPART, :) =  hm_XP5(hm_GIDINV(1:NPART), :)
                end if
             end if   

             !$$--- update other after-ordered arries
             hm_XP1(1:NPART, :)  = m_XP1(hm_GID(1:NPART),:)
             hm_DIS(1:NPART, :)  = m_DIS(hm_GID(1:NPART),:)
             hm_STATU(1:NPART)   = m_STATU(hm_GID(1:NPART))

             if(allocated(hm_XP2)) then
                hm_XP2(1:NPART, :)  = m_XP2(hm_GID(1:NPART), :)
                hm_XP3(1:NPART, :)  = m_XP3(hm_GID(1:NPART), :)
                hm_XP4(1:NPART, :)  = m_XP4(hm_GID(1:NPART), :)
                hm_XP5(1:NPART, :)  = m_XP5(hm_GID(1:NPART), :)
             end if

           !$$--- copy the ordered arries to devices
            do I=1, m_NDEVICE
               ERR = cudaSetDevice(m_DEVICES(I))
               !ERR = cudaStreamCreate(CPYSTRDEVS(I))
               call Copy_From_Host_To_Devices_template(I,m_CPYSTRDEVS(I), dm_WorkSpace%XP1(I)%Data,     &
                                                                          dm_WorkSpace%XP2(I)%Data,     &
                                                                          dm_WorkSpace%XP3(I)%Data,     &
                                                                          dm_WorkSpace%XP4(I)%Data,     &
                                                                          dm_WorkSpace%XP5(I)%Data,     &
                                                                          dm_WorkSpace%DIS(I)%Data,     &
                                                                          dm_WorkSpace%STATU(I)%Data )
            end do

            !$$---Now, we update hm_GIDINV , copy hm_GIDINV  and hm_GID into devices
            do I=1, NPART
               hm_GIDINV(hm_GID(I)) = I
            end do
            do I=1, m_NDEVICE
               ERR = cudaSetDevice(m_DEVICES(I))
               ERR = cudaMemcpyAsync(dm_WorkSpace%GID(I)%Data(1),   hm_GID(1),    dm_NPRT) !, stream=m_CPYSTRDEVS(I))
               ERR = cudaMemcpyAsync(dm_WorkSpace%GIDINV(I)%Data(1),hm_GIDINV(1), dm_NPRT) !, stream=m_CPYSTRDEVS(I))
            end do 

            !$$--- reserve cell information
            dm_WorkSpace%STARTCELL= m_STARTCELL
            dm_WorkSpace%ENDCELL  = m_ENDCELL
            dm_WorkSpace%STARTA   = m_STARTA
            dm_WorkSpace%ENDA     = m_ENDA
            dm_WorkSpace%NPA      = m_NPA
            dm_Neighbors%STARTA   = m_STARTA
            dm_Neighbors%ENDA     = m_ENDA
            dm_Neighbors%NPA      = m_NPA

            !--- for the first time of doing partition, we check if the 
            !    size of neighbor-list is large enough, hm_INC is temperarely used 
            if(hm_HasPart .eq. mp_NOTPART) then
               call DevCopyOut(hm_INC, m_STARTA, dm_Neighbors%KVOIS, UN, m_NPA)
               call SynchronizeDevices()
               if(any(hm_INC .ge. m_mxKVOIS)) then
                  write(*,fmt="(A,I5)") " MDPSCU Warning: the number of neighbors reaches the max permitted value: ", m_mxKVOIS
                  write(*,fmt="(A,I5)") "                 possible miss of neighbores."
                  write(*,fmt="(A,I5)") "                 consider reducing the cutoff radius or increasing the permitted number of neighbors"
                  call ONWARNING(gm_OnWarning)
               end if   
            end if  

            !--- accumulate the number of performing partition 
            hm_HasPart = hm_HasPart + 1

            !NOTE-2019-09-23: since the main process may interapted by, damping, in EventSearch
            !                 or, parall-replicaa MD, we need forcees at the point before the interaption.
            !                 thus we restor the force.
            if(hm_HasPart .eq. C_IUN) then
               call CopyFPFrom_Host_to_Devices(m_FP)
            end if   
                
            ERR = cudaSetDevice(CURDEV) 
        return
  end subroutine Cal_NeighBoreList2C_0_DEV
  !********************************************************************************************

  !**********************************************************************************
  subroutine Cal_NeighBoreList2C_A_DEV(SimBox, CtrlParam, List)
  !***  PURPOSE:  to update the neighbore list of particles with head-link
  !               technique used.
  !
  !               Newton's 3rd law is NOT considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should NOT be applied
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          type(SimMDBox),      dimension(:)          ::SimBox
          type(SimMDCtrl)                            ::CtrlParam
          type(NEIGHBOR_LIST), dimension(:), optional::List

  !----   Local variables

  !--- Device variables and variables to be used in GPU
            call Cal_NeighBoreList2C_0_DEV(SimBox(1), CtrlParam, size(SimBox))
            if(present(List)) then
               call Copyout_NeighboreList_A2(List)
            end if

       return
   end subroutine Cal_NeighBoreList2C_A_DEV
  !*********************************************************************************

  !**********************************************************************************
  subroutine Cal_NeighBoreList2C_B_DEV(SimBox, CtrlParam, List)
  !***  PURPOSE:  to update the neighbore list of particles with head-link
  !               technique used.
  !
  !               Newton's 3rd law is NOT considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should NOT be applied
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          type(SimMDBox)                ::SimBox
          type(SimMDCtrl)               ::CtrlParam
          type(NEIGHBOR_LIST), optional ::List

  !----   Local variables

  !--- Device variables and variables to be used in GPU
            call Cal_NeighBoreList2C_0_DEV(SimBox, CtrlParam, 1)
            if(present(List)) then
               call Copyout_NeighboreList_B2(List, Order=1)
            end if
       return
   end subroutine Cal_NeighBoreList2C_B_DEV
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine  Cal_NearestNeighbor_Kernel(IM, XP, NAPDEV,IA0, NPRT, BSX, BSY, BSZ, PDX, PDY,PDZ, NEAREST, KVOIS,INDI, &
                                                             TKVOIS, TINDI)
  !***  PURPOSE:   to extract the first NEAR nearest neighbors from the neighbore list created on call
  !               Cal_NeighboreList_Kernel2C
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              NPRT:   the actual number of atoms on the device
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !
  !$$   OUTPUT:    TKVOIS,  the number of nerest neighbors
  !                TINDI:   the index for the neighbores
  !
      implicit none
      !--- DUMMY VARIABLES
      integer, value, intent(in)::IM
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer, value, intent(in)::NAPDEV, IA0, NPRT, NEAREST
      real(KINDDF), value::BSX, BSY, BSZ
      integer, value::PDX, PDY,PDZ
      integer, device, intent(in)   ::KVOIS(NAPDEV)
      integer, device, intent(in)   ::INDI(NAPDEV,*)
      integer, device, intent(inout)::TKVOIS(NAPDEV)
      integer, device, intent(inout)::TINDI(NAPDEV,*)

      !Local variables
      integer::IW, IIW
      real(KINDDF)::R2, SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z

      real(KINDDF)::HBX, HBY, HBZ, BX, BY, BZ
      integer::IFPDX, IFPDY, IFPDZ
      integer::NT, NB

      ! Local variables
      integer::IC, IT, IB, LL, NAL, NL, NN, I, J, K, N1
      real(KINDDF)::R2SWAP(mp_MXNEAREST)
      integer::NEARSWAP(mp_MXNEAREST)

            IT  = (threadidx%y-1)*blockdim%x + threadidx%x

            BX = BSX
            BY = BSY
            BZ = BSZ

            HBX = BX*C_HALF
            HBY = BY*C_HALF
            HBZ = BZ*C_HALF

            IFPDX = PDX
            IFPDY = PDY
            IFPDZ = PDZ

           !$$--- size of Block
            NT = blockdim%x*blockdim%y
           !$$--- size of grid
            NB = griddim%x*griddim%y

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1

            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
            do LL=1, NL

              !$$IC -- the id of the atom on the device
              !$$      NOTE: IC is not the id of the same atom in the whole box
              IC= (IB-1)*NT+IT +(LL-1)*NAL

              if(IC.LE.NPRT) then
                 !$$POS-- position of the atom
                 !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                 POSX = XP(IC+IA0,1)
                 POSY = XP(IC+IA0,2)
                 POSZ = XP(IC+IA0,3)

                 !$$-- start calculation of electron density
                 IIW    = KVOIS(IC)
                 NN     = 0
                 R2SWAP(1:mp_MXNEAREST) = 1.D32
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
                    !DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                    !DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                    !DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                    !R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                     R2 = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                     N1 = min(NN+1, NEAREST)
                     if(R2 .lt. R2SWAP(N1)) then
                        do I=1, N1
                           if(R2 .lt. R2SWAP(I) ) then
                              exit
                           end if
                        end do
                        !if(I .le. NEAREST) then
                        !do K=NEAREST, I+1, -1
                         do K=N1, I+1, -1
                            R2SWAP(K)   = R2SWAP(K-1)
                            NEARSWAP(K) = NEARSWAP(K-1)
                        end do
                        R2SWAP(I)   = R2
                        NEARSWAP(I) = J
                        NN          = NN + 1
                        if(NN  .gt. NEAREST) NN = NEAREST
                        !end if
                     end if   
                  end do
                 !$$--- Update neighbor list
                  TKVOIS(IC) = NN
                  TINDI(IC,1:NN) = NEARSWAP(1:NN)
              end if
            end do
        return
  end subroutine Cal_NearestNeighbor_Kernel
  !*********************************************************************************   

  !**********************************************************************************
  subroutine StartOnDevice_Reoder_template(IDEV, CFROM, CTO, XP, NAC, IA1th, NEAREST, KVOIS, INDI, TKVOIS, TINDI)
  !***  PURPOSE:  to exttract the first NEAR nearest neighbors from the neighbore list created on call
  !               Cal_NeighboreList_Kernel2C
  !
  !               NOTE: atom J is on the neighbore list of atom I atom,  atom I is NOT necessarily on
  !                     the neighbore list of atom J
  !
  !    INPUT(CPU):  IDEV,      the index of the device
  !                 CFROM,     the start cell on the device
  !                 CTO,       the end cell on the device
  !                 NEAREST,   the number of nearest neighbors
  !   INPUT(GUP):
  !                XP,        the position of the atoms
  !                NAC,       the number of atoms in the cells
  !                IA1th,     the index of the first atom in a cell
  !
  !   OUTPUT(GPU): TKVOIS,     the number of neighbores of atoms
  !                TINDI,      the index of neighbore atoms
  !
  implicit none
  !----   DUMMY Variables
          integer, intent(in)::IDEV,CFROM, CTO, NEAREST

          integer,      device, dimension(:)   ::NAC,IA1th
          real(KINDDF), device, dimension(:,:) ::XP
          integer,      device, dimension(:)   ::ITYP
          integer,      device, dimension(:)   ::KVOIS
          integer,      device, dimension(:,:) ::INDI
          integer,      device, dimension(:)   ::TKVOIS
          integer,      device, dimension(:,:) ::TINDI

  !----   Local variables
         type(dim3) :: blocks
         type(dim3) :: threads
         integer::STARTA, ENDA, NPRT, CURDEV, BX, BY, NBX, NBY, ERR
         integer, parameter::mp_BLOCKDIMX=256

  !--- start

         !$$--- copy the cell informations into device memeory
            ERR     = cudaGetDevice(CURDEV)
            ERR     = cudaSetDevice(IDEV)
            
            STARTA  = hm_IA1th(CFROM)
            ENDA    = (hm_IA1th(CTO)-1) + hm_NAC(CTO)
            NPRT    = ENDA - STARTA + 1

            BX      = mp_BLOCKSIZE
            BY      = 1
            NBX     = mp_BLOCKDIMX
            NBY     = 1
            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            
            call Cal_NearestNeighbor_Kernel<<<blocks, threads>>>(dm_NPRT, XP, m_NAPDEV,STARTA, NPRT,       &
                                                                 m_BOXSIZE(1), m_BOXSIZE(2), m_BOXSIZE(3), &
                                                                 m_PD(1), m_PD(2), m_PD(3),                &
                                                                 NEAREST, KVOIS, INDI, TKVOIS, TINDI)
            ERR = cudaSetDevice(CURDEV)

         return
   end subroutine StartOnDevice_Reoder_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_0(Nearest)
  !***  PURPOSE:  to sort the neighbore sorted by distances,
  !               only the first NEAREST neighbors are extracted.
  !
  !
  !     INPUT:     Nearest,   the number of neighbors expected
  !
  !     OUTPUT     Simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          integer::Nearest

  !----   Local variables
          integer::ERR, CURDEV, I

            ERR = cudaGetDevice(CURDEV)

            if(Nearest .gt.  mp_MXNEAREST) then
               write(*,*) "MDPSCU Error: the number of required neighbors (NEAREST) in Reorder_NeighBoreList_Nearest_DEV is larger than the permitted value"
               write(*,*) "              NEAREST = ", Nearest, ' vs MXNEAREST = ', mp_MXNEAREST
               write(*,*) "Process to be stopped"
               stop
            end if

            if(Nearest .gt.  m_mxKVOIS) then
               write(*,*) "MDPSCU Error: the number of required neighbors (NEAREST) in Reorder_NeighBoreList_Nearest_DEV is larger than the permitted number of neighbors"
               write(*,*) "              NEAREST = ", Nearest, " vs mxKVOIS = ", m_mxKVOIS
               call ONWARNING(gm_OnWarning)
            end if


           !$$--- to start neighbor calculation on devices
            do I=1, m_NDEVICE
               call StartOnDevice_Reoder_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), &
                                                  dm_WorkSpace%XP(I)%Data,                    &
                                                  dm_WorkSpace%NAC(I)%Data,                   &
                                                  dm_WorkSpace%IA1th(I)%Data,                 &
                                                  Nearest,                                    &
                                                  dm_Neighbors%KVOIS(I)%Data,                 &
                                                  dm_Neighbors%INDI(I)%Data,                  &
                                                  dm_Neighbors%KVOIS(I)%Data,                 &
                                                  dm_Neighbors%INDI(I)%Data)
            end do 

            ERR = cudaSetDevice(CURDEV)
         return
  end subroutine Reorder_NeighBoreList_Nearest_0
  !**********************************************************************************

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_0_A(SimBox, Nearest, List)
  !***  PURPOSE:  to sort the neighbore sorted by distances,
  !               only the first NEAREST neighbors are extracted.
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                Nearest,   the number of neighbors expected
  !
  !     OUTPUT     Simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          type(SimMDBox),      dimension(:)::SimBox
          integer                          ::Nearest
          type(NEIGHBOR_LIST), dimension(:)::List
           
          if(dm_NPRT .ne. SimBox(1)%NPRT*size(Simbox) .or. size(SimBox) .ne. size(List)) then
            write(*,*) "MDPSCU Error: the particle number in Reorder_NeighBoreList_Nearest_DEV is not consistent with that in device0"
            write(*,*) "              the particle number in device0: ", dm_NPRT
            write(*,*) "              the particle number in a box:   ", SimBox(1)%NPRT
            write(*,*) "              the number of SimBox:            ", size(SimBox)
            write(*,*) "              the number of List:              ", size(List)
            write(*,*) "Process to be stopped"
            stop
          end if

          call Reorder_NeighBoreList_Nearest_0(Nearest)
          call Copyout_NeighboreList_A2(List)

          return
  end subroutine Reorder_NeighBoreList_Nearest_0_A        
  !**********************************************************************************

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_0_B(SimBox, Nearest, List)
  !***  PURPOSE:  to sort the neighbore sorted by distances,
  !               only the first NEAREST neighbors are extracted.
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                Nearest,   the number of neighbors expected
  !
  !     OUTPUT     Simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          type(SimMDBox)     ::SimBox
          integer            ::Nearest
          type(NEIGHBOR_LIST)::List

          if(dm_NPRT .ne. SimBox%NPRT) then
            write(*,*) "MDPSCU Error: the particle number in Reorder_NeighBoreList_Nearest_DEV is not consistent with that in device0"
            write(*,*) "              the particle number in device0: ", dm_NPRT
            write(*,*) "              the particle number in a box:   ", SimBox%NPRT
            write(*,*) "Process to be stopped"
            stop
          end if
          call Reorder_NeighBoreList_Nearest_0(Nearest)
          call Copyout_NeighboreList_B2(List, order=1)

          return
  end subroutine Reorder_NeighBoreList_Nearest_0_B        
  !**********************************************************************************  

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_1_A(SimBox, CtrlParam, List)
  !***  PURPOSE:  to sort the neighbore sorted by distances,
  !               only the first NEAREST neighbors are extracted.
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                Nearest,   the number of neighbors expected
  !
  !     OUTPUT     Simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !

  implicit none
  !----   DUMMY Variables
          type(SimMDBox),      dimension(:) ::SimBox
          type(SimMDCtrl)                   ::CtrlParam
          type(NEIGHBOR_LIST), dimension(:) ::List

          integer::NEAREST
 
          if(CtrlParam%NB_NEAREST .gt. 0) then
            NEAREST = CtrlParam%NB_NEAREST
          else
            NEAREST = mp_MXNEAREST
          end if  
          
          call Reorder_NeighBoreList_Nearest_0_A(SimBox, NEAREST, List)
          return
  end subroutine Reorder_NeighBoreList_Nearest_1_A        
  !**********************************************************************************

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_1_B(SimBox, CtrlParam, List)
  !***  PURPOSE:  to sort the neighbore sorted by distances,
  !               only the first NEAREST neighbors are extracted.
  !
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                Nearest,   the number of neighbors expected
  !
  !     OUTPUT     List
  !
 
   implicit none
   !----   DUMMY Variables
           type(SimMDBox)     ::SimBox
           type(SimMDCtrl)    ::CtrlParam
           type(NEIGHBOR_LIST)::List
 
   !----   Local variables
           integer::NEAREST
 
             if(CtrlParam%NB_NEAREST .gt. 0) then
               NEAREST = CtrlParam%NB_NEAREST
             else
               NEAREST = mp_MXNEAREST
             end if  
             
             call Reorder_NeighBoreList_Nearest_0_B(SimBox, NEAREST, List)
          return
   end subroutine Reorder_NeighBoreList_Nearest_1_B
  !***********************************************************************************

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_0_C(Nearest, dList)
  !***  PURPOSE:  to extract the nearest the neighbores sorted by distances.
  !               the output list also on GPU 
  !
  !     INPUT:     Nearest,   the number of nearest neighbors to be extracted
  !
  !     OUTPUT     dList,      the extracted nearest neoighbors
  !

  implicit none
  !----   DUMMY Variables
          integer::Nearest
          type(MDDEVNeighbor_list)::dList

  !----   Local variableds 
          integer::ERR, CURDEV, I

           
          if(Nearest .gt.  mp_MXNEAREST) then
            write(*,*) "MDPSCU Error: the number of required neighbors (NEAREST) in Reorder_NeighBoreList_Nearest_DEV is larger than the permitted value"
            write(*,*) "              NEAREST = ", Nearest, ' vs MXNEAREST = ', mp_MXNEAREST
            write(*,*) "Process to be stopped"
            stop
          end if

          ERR = cudaGetDevice(CURDEV)
          if(dList%NPRT .ne. dm_Neighbors%NPRT) then
             call DevDeallocate(dList%KVOIS)  
             call DevDeallocate(dList%INDI) 
             dList%NPRT    = dm_Neighbors%NPRT
             dList%NPRTB   = dm_Neighbors%NPRTB
             dList%STARTA  = dm_Neighbors%STARTA
             dList%ENDA    = dm_Neighbors%ENDA
             dList%NPA     = dm_Neighbors%NPA
             dList%MXKVOIS = Nearest
             call DevAllocate(dList%KVOIS,   m_NAPDEV,            "NB_KOVIS")
             call DevAllocate(dList%INDI,  (/m_NAPDEV, Nearest/), "NB_INDI")
          end if   

          !$$--- to start neighbor calculation on devices
          do I=1, m_NDEVICE
             call StartOnDevice_Reoder_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), &
                                                dm_WorkSpace%XP(I)%Data,                    &
                                                dm_WorkSpace%NAC(I)%Data,                   &
                                                dm_WorkSpace%IA1th(I)%Data,                 &
                                                Nearest,                                    &
                                                dm_Neighbors%KVOIS(I)%Data,                 &
                                                dm_Neighbors%INDI(I)%Data,                  &
                                                dList%KVOIS(I)%Data,                        &
                                                dList%INDI(I)%Data)
           end do 
           ERR = cudaSetDevice(CURDEV)
          return
  end subroutine Reorder_NeighBoreList_Nearest_0_C        
  !**********************************************************************************
  
  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_A_C(Nearest, dList, hList)
  !***  PURPOSE:  to extract the nearest the neighbores sorted by distances.
  !               and copyout to host
  !
  !     INPUT:     Nearest,   the number of nearest neighbors to be extracted
  !
  !     OUTPUT     dList,      the extracted nearest neoighbors
  !

  implicit none
  !----   DUMMY Variables
          integer                            ::Nearest
          type(MDDEVNeighbor_list)           ::dList
          type(Neighbor_list), dimension(:)  ::hList
          
  !----   Local variableds 

          call Reorder_NeighBoreList_Nearest_0_C(Nearest, dList) 
          call Copyout_NeighboreList_A1(hList,dList)
          return
  end subroutine Reorder_NeighBoreList_Nearest_A_C        
  !**********************************************************************************

  !**********************************************************************************
  subroutine Reorder_NeighBoreList_Nearest_B_C(Nearest, dList, hList)
  !***  PURPOSE:  to extract the nearest the neighbores sorted by distances.
  !               and copyout to host
  !
  !     INPUT:     Nearest,   the number of nearest neighbors to be extracted
  !
  !     OUTPUT     dList,      the extracted nearest neoighbors
  !

  implicit none
  !----   DUMMY Variables
          integer                            ::Nearest
          type(MDDEVNeighbor_list)           ::dList
          type(Neighbor_list)                ::hList
          
  !----   Local variableds 

          call Reorder_NeighBoreList_Nearest_0_C(Nearest, dList) 
          call Copyout_NeighboreList_B1(hList,dList)
          return
  end subroutine Reorder_NeighBoreList_Nearest_B_C        
  !**********************************************************************************


  !**********************************************************************************
  subroutine Cal_NeighBoreList_A_DEV(SimBox, CtrlParam, List)
  !***  PURPOSE:  to calculate the neighbore-list
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !
  implicit none
  !----   DUMMY Variables
          type(SimMDBox),      dimension(:)          ::SimBox
          type(SimMDCtrl)                            ::CtrlParam
          type(NEIGHBOR_LIST), dimension(:), optional::List

  !----   Local variables
             call Cal_NeighBoreList2C_A_DEV(SimBox, CtrlParam)
             if(CtrlParam%NB_NEAREST .gt. 0) then
                if(present(List)) then
                   call Reorder_NeighBoreList_Nearest_0_A(SimBox, CtrlParam%NB_NEAREST, List)
                else 
                   call Reorder_NeighBoreList_Nearest_0(CtrlParam%NB_NEAREST)   
                end if   
             end if
         return
  end subroutine Cal_NeighBoreList_A_DEV
  !********************************************************************************************

  !**********************************************************************************
  subroutine Cal_NeighBoreList_B_DEV(SimBox, CtrlParam, List)
  !***  PURPOSE:  to calculate the neighbore-list
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     simulation box with neighbor-list updated
  !                List,      optional, the neighbore-list on CPU,
  !                           could be used for debug purpuse.
  !
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)              ::SimBox
          type(SimMDCtrl)             ::CtrlParam
          type(NEIGHBOR_LIST),optional::List

  !----   Local variables
             call Cal_NeighBoreList2C_B_DEV(SimBox, CtrlParam)
             if(CtrlParam%NB_NEAREST .gt. 0) then
                if(present(List)) then
                  call Reorder_NeighBoreList_Nearest_0_B(SimBox, CtrlParam%NB_NEAREST, List)
                else
                  call Reorder_NeighBoreList_Nearest_0(CtrlParam%NB_NEAREST)
                end if
             end if
         return
  end subroutine Cal_NeighBoreList_B_DEV
  !*********************************************************************************
  
  !********************************************************************************************
  integer function Get_Number_of_Cells()
  implicit none
           Get_Number_of_Cells = m_NC
           return
  end function Get_Number_of_Cells

  !********************************************************************************************
  integer function Get_Number_of_AtomPerCell()
  implicit none
           Get_Number_of_AtomPerCell = sum(hm_NAC)/m_NC
           return
  end function Get_Number_of_AtomPerCell
  !********************************************************************************************

  !********************************************************************************************
  subroutine GetCellInform(TNC, NCELL, MXNAC)
  implicit none
  integer::TNC, NCELL(3), MXNAC

            TNC = m_NC
            NCELL = m_NCELL
            MXNAC = m_MXNAC0
            return
 end subroutine GetCellInform
 !********************************************************************************************

 !********************************************************************************************
  integer function Get_Number_of_PermittedNeighbor()
  implicit none
           Get_Number_of_PermittedNeighbor = m_mxKVOIS
           return
  end function Get_Number_of_PermittedNeighbor

!********************************************************************************************
end module MD_NeighborsList_GPU


