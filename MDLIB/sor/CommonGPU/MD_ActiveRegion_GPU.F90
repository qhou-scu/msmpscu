  module MD_ActiveRegion_GPU
  !**** DESCRIPTION: to make particles active or deactive
  !
  !                  ______________________________________________________
  !                  HOU Qing, Oct, 2017

  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_Globle_Variables_GPU
  implicit none
  !**********************
         integer,        parameter, private::mp_BLOCKSIZE=256
         integer,        parameter, private::mp_BLOCKDIMX=128
         !--- working array to mark the activation seeds from which the activative
         !    region to be calculated
         integer,        dimension(:),   allocatable, private::hm_ActiveSeed
         integer,        dimension(:,:), allocatable, private::hm_SwapSeed

         !---  The type definition of work space for activation region on device 
         type(DevVec_I), dimension(m_MXDEVICE),  private::dm_ActiveSeed
         type(DevVec_I), dimension(m_MXDEVICE),  private::dm_ActiveMark
         !---
         integer,                private:: m_KEEPFLAG = 0
         integer,                private:: m_CENTPART(mp_MXGROUP)
         real(KINDDF),           private:: m_CM(mp_MXGROUP)
         integer,      constant, private:: dcm_CENTPART(mp_MXGROUP)
         real(KINDDF), constant, private:: dcm_CM(mp_MXGROUP)

        !__ interface for procedure of setting activative region
        private::CREATEACTIVESEEDS
        abstract interface
           subroutine CREATEACTIVESEEDS(SimBox, CtrlParam, Seeds)
           use MD_CONSTANTS
           use MD_SimboxArray
           use MD_TYPEDEF_SimMDCtrl
           implicit none
             type(SimMDBox), dimension(:) ::SimBox
             type(SimMDCtrl)              ::CtrlParam
             integer, dimension(:)        ::Seeds
           end subroutine CREATEACTIVESEEDS
        end interface
        procedure(CREATEACTIVESEEDS), pointer, private::m_pCreatAS=>null()
        !--- End interface

      !--- Interfaces for external calls  
      !-------------------------------------------------
        private:: Allocate_Working_Variables
      !-------------------------------------------------
        private:: Active_All_Kernel,            &
                  Active_All_ActiveRegion0,     &
                  Active_All_ActiveRegion1
        public::  Active_All_ActiveRegion_DEV
        interface Active_All_ActiveRegion_DEV
           module procedure Active_All_ActiveRegion0  
           module procedure Active_All_ActiveRegion1    
        end interface Active_All_ActiveRegion_DEV
      !-------------------------------------------------
        private:: MarkSeedNeighbore_Kernel,     &
                  MarkSeedNeighbors_template,   &
                  Active_Marked_Kernel,         &
                  Active_Marked_template,       &
                  ActiveByNeigbors0,            &
                  ActiveByNeigbors1
        public::  Activate_RegionByNeigbors_DEV
        interface Activate_RegionByNeigbors_DEV
           module procedure ActiveByNeigbors0
           module procedure ActiveByNeigbors1 
        end interface Activate_RegionByNeigbors_DEV
      !-------------------------------------------------
        private:: MarkSeedCell_Kernel0,         &
                  MarkSeedCell_template0,       &
                  Active_InCells_Kernel,        &               
                  Active_InCells_template0,     &
                  ActiveByCells0,               &
                  ActiveByCells1
        public::  Activate_RegionByCells_DEV
        interface Activate_RegionByCells_DEV
           module procedure ActiveByCells0
           module procedure ActiveByCells1
        end interface Activate_RegionByCells_DEV
      !-------------------------------------------------
        public::  ActivateRegion_DEV 

      !-------------------------------------------------
        public::  Clear_ActiveRegion_DEV
      !-------------------------------------------------
        private:: CreateSeedByType_Kernel,   &
                  CreateSeedByType0,         & 
                  CreateSeedByType1,         &
                  CreateSeedByType_H
        public::  CreateActiveSeedByType_Dev
        interface CreateActiveSeedByType_Dev
           module procedure CreateSeedByType0 
           module procedure  CreateSeedByType_H
        end interface CreateActiveSeedByType_Dev  
      !-------------------------------------------------
        private:: CreateSeedByEkin_Kernel,   &
                  CreateSeedByEKin0,         & 
                  CreateSeedByEKin1,         &
                  CreateSeedByEKin_H
        public::  CreateActiveSeedByEKin_Dev
        interface CreateActiveSeedByEKin_Dev
           module procedure CreateSeedByEKin0 
           module procedure CreateSeedByEKin_H 
        end interface CreateActiveSeedByEKin_Dev  
      !-------------------------------------------------
        private:: DeActive_All_Kernel,        &
                  DeActive_All_ActiveRegion0, &
                  DeActive_All_ActiveRegion1 
        public::  DeActive_All_ActiveRegion_DEV
        interface DeActive_All_ActiveRegion_DEV
           module procedure DeActive_All_ActiveRegion0
           module procedure DeActive_All_ActiveRegion1
        end interface DeActive_All_ActiveRegion_DEV
      !------------------------------------
        public:: Initialize_ActiveRegion_DEV
      !------------------------------------
        public:: SetCreatSeedProc_ActiveRegion_DEV
  contains

  !****************************************************************************************
  !****************************************************************************************
  subroutine SetCreatSeedProc_ActiveRegion_DEV(ExtCreatSeedPro)
  !***  DESCRIPTION: to set the procedure that used to construct active region
  !
  !--- INPUT: ACTIVEREGIONPRO,    the external subroutine of updating active region
  !
  implicit none

  !--- dummy variables and subroutines
   optional::ExtCreatSeedPro
   external::ExtCreatSeedPro

  !--- interface to the external routine -------------------
   interface
       subroutine ExtCreatSeedPro(SimBox, CtrlParam, Seeds)
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
           type(SimMDBox), dimension(:)::SimBox
           type(SimMDCtrl)             ::CtrlParam
           integer,        dimension(:)::Seeds
       end subroutine ExtCreatSeedPro
   end interface
  !--- END INTERFACE --------------------------------

  !--- Local variables
          if(present(ExtCreatSeedPro)) then
             m_pCreatAS=>ExtCreatSeedPro
          else
             m_pCreatAS=>null()
          end if

          return
      end subroutine SetCreatSeedProc_ActiveRegion_DEV
  !****************************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Variables(WorkSpace)
  !***  PURPOSE:   to allocate working space
  !
      implicit none
      !--- dummy variables
      type(MDDEVWorkSpace) ::WorkSpace
      !--- Local vairables
      integer::ERR

             !$$--- to check if particle number have been initialized
             !$$    each device
              if(WorkSpace%NPRT .le. 0) then
                 write(*,fmt="(A)")  ' MDPSCU Error: dm_NPRT is found zero in ActiveRegion_module. '
                 write(*,fmt="(A)")  '               Process to be stopped'
                 stop
              end if

              allocate(hm_ActiveSeed(WorkSpace%NPRT), hm_SwapSeed(WorkSpace%NPRT,m_NDEVICE), STAT=ERR)
              if(ERR) then
                  write(*,fmt="(A)") "MDPSCU Error: fail to allocate memory for ActiveSeed in ActiveRegion module"
                  write(*,fmt="(A)") '              Process to be stopped'
                  stop "1"
              end if
              hm_ActiveSeed = 1

             !$$--- allocate working space on device 1
              call DevAllocate(dm_ActiveSeed, WorkSpace%NPRT, "dm_ActiveSeed")
              call DevAllocate(dm_ActiveMark, WorkSpace%NPRT, "dm_ActiveMark")

       return
  end subroutine Allocate_Working_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_ActiveRegion_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to intialize the parameters to be used in active region calculation
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT:   the working spaces allocated
  !
  !     NOTE:    Inialization of ActiveRegion calculation should be called after
  !              calling Initialize_NeighboreList_DEV (see: MD_NeighborsList_GPU.F90)
  !
    implicit none
       !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(in) ::SimBox
      type(SimMDCtrl),              intent(in) ::CtrlParam
      !--- Local variables

               call Clear_ActiveRegion_DEV()
               call Allocate_Working_Variables(dm_WorkSpace)

               m_CM       = SimBox(1)%CM
               m_KEEPFLAG = 0
               return
  end subroutine Initialize_ActiveRegion_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_ActiveRegion_DEV()
    !***  PURPOSE:   to intialize the parameters to be used in active region calculation
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !
    !     OUTPUT:   the working spaces allocated
    !
    implicit none
       !--- dummy vaiables
       !--- Local variables
       integer::I

              if(allocated(hm_ActiveSeed)) deallocate(hm_ActiveSeed)
              if(allocated(hm_SwapSeed))   deallocate(hm_SwapSeed)

         !$$--- allocate working space on device 1
              call DevDeallocate(dm_ActiveSeed)
              call DevDeallocate(dm_ActiveMark)
              return
  end subroutine Clear_ActiveRegion_DEV
  !****************************************************************************

  !****************************************************************************
  attributes(global) subroutine DeActive_All_Kernel(NAPDEV, NPRT, STATU)
  !***  PURPOSE:   Kernel to deactive all atoms by modifying the STATU of atoms
  !
  !$$   INPUT:     NAPDEV,    permitted number of atoms on the device
  !$$              NPRT,      actuall number of atoms on the device
  !
  !$$   OUTPUT     STATU,    statu with bit 1 set to zero
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in) :: NAPDEV, NPRT
          integer, device, intent(out):: STATU(*)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms thisKernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT+IT +(LL-1)*NAL

                if(IC .le. NPRT) then
                   STATU(IC)  = ibclr(STATU(IC), CP_STATU_ACTIVE_BITPOS)
                end if
            end do
        return
  end subroutine DeActive_All_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine DeActive_All_ActiveRegion0(WorkSpace)
  !***  PORPOSE: to interface to host process
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  implicit none
  !----   DUMMY Variables
          type(MDDEVWorkSpace)              ::WorkSpace
  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
   !--- local variables       
          integer::I, ERR, CURDEV

          blocks  = dim3(mp_BLOCKDIMX, 1, 1)
          threads = dim3(mp_BLOCKSIZE, 1, 1)

          ERR = cudaGetDevice(CURDEV)
          do I = 1, m_NDEVICE 
             ERR = cudaSetDevice(m_DEVICES(I))
             call DeActive_All_Kernel<<<blocks, threads>>>(WorkSpace%NAPDEV, WorkSpace%NPA(I), WorkSpace%STATU(I)%Data)
          end do 
          ERR = cudaSetDevice(CURDEV)
     return
  end subroutine DeActive_All_ActiveRegion0
  !****************************************************************************************

  !****************************************************************************************
  subroutine DeActive_All_ActiveRegion1(SimBox,CtrlParam)
  !***  PORPOSE: to interface to host process
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:) ::SimBox
          type(SimMDCtrl)              ::CtrlParam
  !--- 

          call DeActive_All_ActiveRegion0(dm_WorkSpace)
     return
  end subroutine DeActive_All_ActiveRegion1
  !****************************************************************************************  

  !*******************************************************************************************
  attributes(global) subroutine Active_All_Kernel(NAPDEV, NPRT, STATU)
  !***  PURPOSE:   Kernel to active all atoms IN box
  !
  !$$   INPUT:     NAPDEV,    permitted number of atoms on the device
  !$$              NPRT,      actuall number of atoms on the device
  !
  !$$   OUTPUT     STATU,    statu with bit 1 set to 1
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in) :: NAPDEV, NPRT
          integer, device, intent(out):: STATU(*)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms this Kernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT + IT + (LL-1)*NAL

                if(IC .le. NPRT) then
                   if(iand(STATU(IC), CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX ) then
                      STATU(IC)  = ibset(STATU(IC), CP_STATU_ACTIVE_BITPOS)
                   end if
                end if
            end do
        return
  end subroutine Active_All_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine Active_All_ActiveRegion0(WorkSpace)
  !***  PORPOSE: to active all in-box atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !
  implicit none
  !----   DUMMY Variables
          type(MDDEVWorkSpace)              ::WorkSpace
  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
  !--- Local variables
          integer::I, ERR, CURDEV

          blocks  = dim3(mp_BLOCKDIMX, 1, 1)
          threads = dim3(mp_BLOCKSIZE, 1, 1)

          ERR = cudaGetDevice(CURDEV)
          do I = 1, m_NDEVICE 
             ERR = cudaSetDevice(m_DEVICES(I))
             call Active_All_Kernel<<<blocks, threads>>>(WorkSpace%NAPDEV, WorkSpace%NPA(I), WorkSpace%STATU(I)%Data)
          end do 
          ERR = cudaSetDevice(CURDEV)

      return
  end subroutine Active_All_ActiveRegion0
  !****************************************************************************************

  !****************************************************************************************
  subroutine Active_All_ActiveRegion1(SimBox,CtrlParam)
  !***  PORPOSE: to active all in-box atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:) ::SimBox
          type(SimMDCtrl)              ::CtrlParam
  !--- Device variables and variables to be used in GPU
  !--- Local variables

          call Active_All_ActiveRegion0(dm_WorkSpace)
       return
  end subroutine Active_All_ActiveRegion1
  !****************************************************************************************  

  !****************************************************************************************
  attributes(global) subroutine Active_Marked_Kernel(IM, NPRT, IA0, MARK, NAPDEV, STATU)
  !***  PURPOSE:   Kernel to active all atoms IN box
  !
  !$$   INPUT:     IM:        the number of atom in whole box (size of MARK)
  !$$              NPRT,      actuall number of atoms on the device
  !$$              IA0:       the index (in the whole box)of the fisrt atom on the device
  !$$              MARK,      mark of atoms to be activated
  !$$              NAPDEV,    permitted number of atoms on the device
  !
  !$$   OUTPUT     STATU,    statu with bit 1 set to zero
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in) :: IM, NPRT, IA0, NAPDEV
          integer, device, intent(in) :: MARK(IM)
          integer, device, intent(out):: STATU(NAPDEV)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms this Kernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT+IT +(LL-1)*NAL

                if(IC .le. NPRT) then
                   if(MARK(IC+IA0) .gt. 0 ) then
                      STATU(IC)  = ibset(STATU(IC), CP_STATU_ACTIVE_BITPOS)
                   end if
                end if
            end do
        return
  end subroutine Active_Marked_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine Active_Marked_template(MARK, WorkSpace)
  !***  PORPOSE:  template for devices invoking the Kernel
  !
  !     INPUT:     IDEV,      ID of the device
  !                MARK,      mark of atoms to be activated
  !     OUTPUT:    STATU,     statu of atoms
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          type(DevVec_I),     dimension(:)::MARK
          type(MDDEVWorkSpace)            ::WorkSpace

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::I, ERR, CURDEV, IA0

            ERR = cudaGetDevice(CURDEV)
            blocks  = dim3(mp_BLOCKDIMX, 1, 1)
            threads = dim3(mp_BLOCKSIZE, 1, 1)
            do I=1, m_NDEVICE
               ERR = cudaSetDevice(m_DEVICES(I)) 
               IA0 = WorkSpace%STARTA(I) - 1
               call Active_Marked_Kernel<<<blocks, threads>>>(WorkSpace%NPRT, WorkSpace%NPA(I), IA0, &
                                                              MARK(I)%Data,   WorkSpace%NAPDEV,      &
                                                              WorkSpace%STATU(I)%Data)
            end do    
            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine Active_Marked_template
  !****************************************************************************************

  !*******************************************************************************************
  attributes(global) subroutine CreateSeedByType_Kernel(IM, NPRT, IA0, NAPDEV, ITYPE, STATU, SEEDS)
  !***  PURPOSE:   Kernel to create the seeds by type
  !
  !$$   INPUT:     IM:        the number of atom in whole box (size of MARK)
  !$$              NPRT,      actuall number of atoms on the device
  !$$              IA0:       the index (in the whole box)of the fisrt atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              ITYPE,     type of atoms
  !$$              STATU,     statu of atoms
  !
  !$$   OUTPUT:    SEEDS,    mark of atoms to be activated
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in)  :: IM, NPRT, IA0, NAPDEV
          integer, device, intent(in)  :: ITYPE(IM)
          integer, device, intent(in)  :: STATU(NAPDEV)
          integer, device, intent(out) :: SEEDS(IM)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL,ITYP, IA
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms thisKernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT+IT +(LL-1)*NAL
                
                if(IC .le. NPRT) then
                   IA=  IC + IA0
                   if( (iand(STATU(IC), CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX) .and. &
                       (dcm_CENTPART(ITYPE(IA)) .gt. 0) ) then
                        SEEDS(IA) = SEEDS(IA) + 1
                   end if
                end if
            end do
        return
  end subroutine CreateSeedByType_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine CreateSeedByType0(ATYP, WorkSpace, Seeds)
  !***  PORPOSE: to create active seed by atoms type
  !
  !     INPUT:   ATYP,  the flag indicating which type of atoms to be active,
  !
  implicit none
  !----   DUMMY Variables
          integer, dimension(:), intent(in) :: ATYP
          type(MDDEVWorkSpace)              :: WorkSpace
          type(DevVec_I), dimension(:)      :: Seeds
  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
  !--- local variables 
          integer::I, IA0, ERR, CURDEV

            m_CENTPART = 0
            m_CENTPART(1:min(size(m_CENTPART), size(ATYP)) ) = ATYP(1:min(size(m_CENTPART), size(ATYP)) )

            blocks  = dim3(mp_BLOCKDIMX, 1, 1)
            threads = dim3(mp_BLOCKSIZE, 1, 1)
            ERR = cudaGetDevice(CURDEV)
            !$$--- mark the seed
            do I=1, m_NDEVICE
               ERR = cudaSetDevice(m_DEVICES(I))
               dcm_CENTPART = m_CENTPART            
               IA0 = WorkSpace%STARTA(I) - 1                                  
               call CreateSeedByType_Kernel<<<blocks, threads>>>(WorkSpace%NPRT, WorkSpace%NPA(I), IA0, WorkSpace%NAPDEV, &
                                                                 WorkSpace%ITYP(I)%Data, WorkSpace%STATU(I)%Data,         &
                                                                 Seeds(I)%Data)                                                 
            end do
            ERR = cudaSetDevice(CURDEV)
      return
  end subroutine CreateSeedByType0
  !****************************************************************************************

  !****************************************************************************************
  subroutine CreateSeedByType1(ATYP)
  !***  PORPOSE: to create active seed by atoms type
  !
  !     INPUT:     ATYP,  the flag indicating which type of atoms to be active,
  !
  implicit none
  !----   DUMMY Variables
          integer, dimension(:), intent(in) :: ATYP
  !--- Device variables and variables to be used in GPU
  !--- local variables 

          call CreateSeedByType0(ATYP, dm_WorkSpace, dm_ActiveSeed)

      return
  end subroutine CreateSeedByType1
  !****************************************************************************************

  !****************************************************************************************
  subroutine CreateSeedByType_H(ATYP, hSeeds)
  !***  PORPOSE: to create active seed by atoms type
  !
  !     INPUT:     ATYP,  the flag indicating which type of atoms to be active,
  !
  implicit none
  !----   DUMMY Variables
          integer, dimension(:), intent(in) :: ATYP
          integer, dimension(:)             :: hSeeds
  !--- Device variables and variables to be used in GPU
  !--- local variables 
          integer::I

          call CreateSeedByType1(ATYP)
          do I=1, m_NDEVICE
            call DevCopyOut(hm_SwapSeed(:,I), dm_ActiveSeed(I), dm_NPRT)
         end do   
         hSeeds = 0
         do I=1, m_NDEVICE
            hSeeds(:) = hSeeds(:) + hm_SwapSeed(:,I)
         end do  

      return
  end subroutine CreateSeedByType_H
  !****************************************************************************************

  !*******************************************************************************************
  attributes(global) subroutine CreateSeedByEkin_Kernel(IM, NPRT, IA0, NAPDEV, ITYP, XP1, STATU, E0, SEEDS)
  !***  PURPOSE:   Kernel to create the seeds by type
  !
  !$$   INPUT:     IM:        the number of atom in whole box (size of MARK)
  !$$              NPRT,      actuall number of atoms on the device
  !$$              IA0:       the index (in the whole box)of the fisrt atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              XP1,       velocity of atoms
  !$$              STATU,     statu of atoms
  !$$              E0,        energy threthhold
  !$$
  !$$   OUTPUT:    SEEDS,    mark of atoms to be activated
  !
  implicit none
  !----   DUMMY Variables
          integer,      value,  intent(in)  :: IM, NPRT, IA0, NAPDEV
          integer,      device, intent(in)  :: ITYP(IM)
          real(KINDDF), device, intent(in)  :: XP1(NAPDEV,*)
          integer,      device, intent(in)  :: STATU(NAPDEV)
          real(KINDDF), value,  intent(in)  :: E0
          integer,      device, intent(out) :: SEEDS(IM)

  !----   Local variables
          integer      ::IT, IB, IC, NT, NB, NAL, NL, LL,ITY, IA
          real(KINDDF) ::EK, CM0

  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms thisKernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT+IT +(LL-1)*NAL

                if(IC .le. NPRT) then
                   IA = IC + IA0
                   if( (iand(STATU(IC), CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX)) then
                        EK = C_HALF*dcm_CM(ITYP(IA))*(XP1(IC,1)*XP1(IC,1)+XP1(IC,2)*XP1(IC,2)+XP1(IC,3)*XP1(IC,3))
                        if(EK .ge. E0) then
                           SEEDS(IA)  = SEEDS(IA) + 1
                        end if
                   end if
                end if
            end do
        return
  end subroutine CreateSeedByEkin_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine CreateSeedByEKin0(E0, WorkSpace, Seeds)
  !***  PORPOSE: to create active seed by atoms type
  !
  !     INPUT:     ATYP,  the flag indicating which type of atoms to be active
  !
  use MD_Globle_Variables_GPU, ONLY:m_STARTCELL, m_ENDCELL, dm_WorkSpace, COPY_IN_NOSHIFT_template
  implicit none
  !----   DUMMY Variables
          real(KINDDF), intent(in)      :: E0
          type(MDDEVWorkSpace)          :: WorkSpace
          type(DevVec_I), dimension(:)  :: Seeds
  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::I, IA0, ERR, CURDEV

           blocks  = dim3(mp_BLOCKDIMX, 1, 1)
           threads = dim3(mp_BLOCKSIZE, 1, 1)
           ERR = cudaGetDevice(CURDEV)
           !$$--- mark the seed
            do I=1, m_NDEVICE 
               ERR    = cudaSetDevice(m_DEVICES(I))
               dcm_CM = m_CM
               IA0    = WorkSpace%STARTA(I)-1
               call CreateSeedByEkin_Kernel<<<blocks, threads>>>(WorkSpace%NPRT, WorkSpace%NPA(I), IA0, WorkSpace%NAPDEV,                     &
                                                                WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%STATU(I)%Data, &
                                                                E0, Seeds(I)%Data)
            end do
            ERR = cudaSetDevice(CURDEV)
      return
  end subroutine CreateSeedByEKin0
  !****************************************************************************************

  !****************************************************************************************
  subroutine CreateSeedByEKin1(E0)
  !***  PORPOSE: to create active seed by atoms type
  !
  !     INPUT:     ATYP,  the flag indicating which type of atoms to be active
  !
  implicit none
   !----   DUMMY Variables
           real(KINDDF), intent(in) :: E0
   !---    Device variables and variables to be used in GPU
 
             call CreateSeedByEKin0(E0, dm_WorkSpace, dm_ActiveSeed)
       return
   end subroutine CreateSeedByEKin1
  !****************************************************************************************

  !****************************************************************************************
  subroutine CreateSeedByEKin_H(E0, hSeeds)
  !***  PORPOSE: to create active seed by atoms type
  !
  !     INPUT:     ATYP,  the flag indicating which type of atoms to be active,
  !
  implicit none
  !----   DUMMY Variables
          real(KINDDF), intent(in) :: E0
          integer, dimension(:)    :: hSeeds
  !--- Device variables and variables to be used in GPU
  !--- local variables 
          integer::I

          call CreateSeedByEKin1(E0)
          do I=1, m_NDEVICE
            call DevCopyOut(hm_SwapSeed(:,I), dm_ActiveSeed(I), dm_NPRT)
         end do   
         hSeeds = 0
         do I=1, m_NDEVICE
            hSeeds(:) = hSeeds(:) + hm_SwapSeed(:,I)
         end do  

      return
  end subroutine CreateSeedByEKin_H
  !****************************************************************************************
 
  !*******************************************************************************************
  attributes(global) subroutine MarkSeedNeighbore_Kernel(IM, NPRT, IA0, SEED, NAPDEV, KVOIS, INDI, MARK)
  !***  PURPOSE:   to mark the seed and their neighbore as atoms to be activated
  !
  !$$   INPUT:     IM:        the number of atom in whole box (size of MARK)
  !$$              NPRT,      actuall number of atoms on the device
  !$$              IA0:       the index (in the whole box)of the fisrt atom on the device
  !$$              GID,       global ID of particles on devices
  !$$              SEED,      id of atoms taken as a seed for construct active region
  !$$              NAPDEV:     the max number of atoms on a device
  !$$              KVOIS:     the number of neighbors for atoms
  !$$              INDI:      the index for the neighbores
  !$$
  !
  !$$   OUTPUT     MARK,      statu with bit 1 set to zero
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in) :: IM, NPRT, IA0, NAPDEV
          integer, device, intent(in) :: SEED(IM)
          integer, device, intent(in) :: KVOIS(NAPDEV)
          integer, device, intent(in) :: INDI(NAPDEV,*)

          integer, device, intent(out) :: MARK(IM)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL, IW, IIW, J
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms thiskernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT+IT +(LL-1)*NAL

                if(IC .le. NPRT) then
                   if(SEED(IC+IA0) .gt. 0) then
                      !$$--- mark the seed
                      MARK(IC+IA0) = 1
                      !$$--- mark their neighbores as to be activated
                      IIW = KVOIS(IC)
                      do IW=1, IIW
                         !$$--- NOTE: the particles index of neighbore-list is
                         !$$          the index of particle in the whole box
                         J = INDI(IC,IW)
                         MARK(J) = MARK(J) + 1
                      end do
                   end if
                end if
            end do
        return
  end subroutine MarkSeedNeighbore_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine MarkSeedNeighbors_template(WorkSpace, Neighbors, SEED, MARK)
  !***  PORPOSE:  template for devices invoking the Kernel
  !
  !     INPUT:     IDEV,      ID of the device
  !                KVOIS:     the number of neighbors for atoms
  !                INDI:      the index for the neighbores
  !$$              SEED,      id of atoms taken as a seed for construct active region
  !     OUTPUT:    MAKR,      mark of atoms to be actived
  !
  use MD_NeighborsList_GPU
  implicit none
  !----   DUMMY Variables
          type(MDDEVWorkSpace)           ::WorkSpace
          type(MDDEVNeighbor_list)       ::Neighbors
          type(DevVec_I),  dimension(:)  ::SEED
          type(DevVec_I),  dimension(:)  ::MARK

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::I, ERR, CURDEV, IA0

            ERR = cudaGetDevice(CURDEV)
            blocks  = dim3(mp_BLOCKDIMX, 1, 1)
            threads = dim3(mp_BLOCKSIZE,  1, 1)
            do I=1, m_NDEVICE
               ERR = cudaSetDevice(m_DEVICES(I)) 
               IA0 = WorkSpace%STARTA(I) - 1
               call MarkSeedNeighbore_Kernel<<<blocks, threads>>>(WorkSpace%NPRT, WorkSpace%NPA(I), IA0,     &
                                                    SEED(I)%Data, WorkSpace%NAPDEV, Neighbors%KVOIS(I)%Data, &
                                                    Neighbors%INDI(I)%Data, MARK(I)%Data)
            end do
            ERR = cudaSetDevice(CURDEV)

       return
  end subroutine MarkSeedNeighbors_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine ActiveByNeigbors0(SimBox, CtrlParam, WorkSpace, Neighbors)
  !***  PORPOSE: to active all in-box atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !
  use MD_NeighborsList_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:) ::SimBox
          type(SimMDCtrl)              ::CtrlParam
          type(MDDEVWorkSpace)         ::WorkSpace
          type(MDDEVNeighbor_list)     ::Neighbors
  !--- Device variables and variables to be used in GPU
          integer::I, L
          !integer, dimension(:),allocatable::Idebug

          !allocate(Idebug(WorkSpace%NPRT))
          !$$--- create the inital seeds
          if(iand(CtrlParam%AR_METHOD, CP_USERDEF_AR) .eq. CP_USERDEF_AR) then
             if(.not.associated(m_pCreatAS) ) then
                write(*,fmt="(' MDPSCU Error: user defined procedure is needed for creating active region', A16)")
                write(*,fmt="('               but the procedure is not given', A16)")
                write(*,fmt="('               check input of &ACTIVEREGSUBCTL subsection', A16)")
                write(*,fmt="('   Process to be stopped', A16)")
                stop
             end if
             call m_pCreatAS(SimBox, CtrlParam, hm_ActiveSeed)
             do I=1, m_NDEVICE
               call DevCopyIn(hm_ActiveSeed, dm_ActiveSeed(I), WorkSpace%NPRT)
             end do   

          else if(iand(CtrlParam%AR_METHOD, CP_LWORD) .eq. 0) then
               !$$--- if no intrinic method of creating seed, we return 
               return
          else  
             call DevSet(dm_ActiveSeed, 0)
             if(iand(CtrlParam%AR_METHOD, CP_CENTPART_AR) .eq. CP_CENTPART_AR)  then
                call CreateSeedByType0(CtrlParam%AR_CENTPART, WorkSpace, dm_ActiveSeed)
             end if
             if(iand(CtrlParam%AR_METHOD, CP_EKIN_AR) .eq. CP_EKIN_AR) then
                !--- NOTE: the input value of AR_EKIN is in eV
                call CreateSeedByEkin0(CtrlParam%AR_EKIN*CP_EVERG, WorkSpace, dm_ActiveSeed)
             end if
          end if
          call DevMakeCopy(WorkSpace%NPRT, dm_ActiveSeed, dm_ActiveMark)

          
          !$$-- marking the active particles
          do L=1, CtrlParam%AR_Extend
             !$$--- mark the seed and their neighbores
             call MarkSeedNeighbors_template(WorkSpace, Neighbors, dm_ActiveSeed, dm_ActiveMark)
             if(m_NDEVICE .gt. 1) then
               !$$--- because the neighbores could be on different devices when
               !$$    multiple devices are used. we should synchronize the dm_ActiveSeed
               !$$    on different devices
               do I=1, m_NDEVICE
                   call DevCopyOut(hm_SwapSeed(:,I), dm_ActiveMark(I), WorkSpace%NPRT)
                end do   
                hm_ActiveSeed = 0
                do I=1, m_NDEVICE
                  hm_ActiveSeed(:) = hm_ActiveSeed(:) + hm_SwapSeed(:,I)
                end do  
                do I=1, m_NDEVICE
                  call DevCopyIn(hm_ActiveSeed, dm_ActiveSeed(I), WorkSpace%NPRT)
                end do   
             else
               call DevMakeCopy(WorkSpace%NPRT, dm_ActiveMark, dm_ActiveSeed)
             end if
          end do

          !$$-- actvate the marked seeds for active region
          !$$   if KEEPFLG > 0, the statu of atom will be keep active
          call Active_Marked_template(dm_ActiveMark, WorkSpace)
          !deallocate(Idebug)
     return
  end subroutine ActiveByNeigbors0
  !****************************************************************************************

  !****************************************************************************************
  subroutine ActiveByNeigbors1(SimBox,CtrlParam)
  !***  PORPOSE: to active all in-box atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !
  use MD_NeighborsList_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:) ::SimBox
          type(SimMDCtrl)              ::CtrlParam
  !--- Device variables and variables to be used in GPU
          integer::I, L

          if(iand(CtrlParam%AR_METHOD, CP_KEEP_AR) .ne. CP_KEEP_AR) then
             m_KEEPFLAG = 0
          end if

          !$$--- NOTE: initially , m_KEEPFLAG is set as 0
          if(m_KEEPFLAG .le. 0) then
             call DeActive_All_ActiveRegion0(dm_WorkSpace)
             if(iand(CtrlParam%AR_METHOD, CP_KEEP_AR) .eq. CP_KEEP_AR) then
                m_KEEPFLAG = 1
             end if
          end if

          call ActiveByNeigbors0(SimBox, CtrlParam, dm_WorkSpace, dm_Neighbors)
     return
  end subroutine ActiveByNeigbors1
  !****************************************************************************************

  !*******************************************************************************************
  attributes(global) subroutine MarkSeedCell_Kernel0(IM, NPRT, IA0, SEEDS, NAPDEV, INC, MC)
  !***  PURPOSE:   KERNEL to mark the cells that the activation seeds in
  !
  !$$   INPUT:     IM,        total number of atoms, size of SEEDS
  !$$              NPRT,      actuall number of atoms on the device
  !$$              SEEDS,     ID of seed particles
  !$$              NAPDEV,    permitted number of atoms on the device, size for INC
  !$$              INC,       cell ID of particles, calculated in partitioning      
  !
  !$$   OUTPUT     MC,        cell ID that the seeds in
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in) :: IM, NAPDEV, NPRT, IA0
          integer, device, intent(in) :: SEEDS(IM)
          integer, device, intent(in) :: INC(NAPDEV)
          integer, device, intent(out):: MC(*)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms this Kernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT + IT + (LL-1)*NAL

                if(IC .le. NPRT ) then
                   if(SEEDS(IC+IA0) .gt. 0) then
                     MC(INC(IC)) = 1
                   end if
                end if
            end do
        return
  end subroutine MarkSeedCell_Kernel0
  !****************************************************************************************

  !****************************************************************************************
  subroutine MarkSeedCell_template0(WorkSpace, SEED, MC)
   !***  PORPOSE:  template for devices invoking the Kernel
   !
   !     INPUT:     WorkSpace,  workspace on devices
   !                SEED,     id of atoms taken as a seed for construct active region
   !     OUTPUT:    MC,     marked cells to be included in activation
   !
   implicit none
   !----   DUMMY Variables
           type(MDDEVWorkSpace)           ::WorkSpace
           type(DevVec_I),  dimension(:)  ::SEED
           type(DevVec_I),  dimension(:)  ::MC
 
   !--- Device variables and variables to be used in GPU
           type(dim3) :: blocks
           type(dim3) :: threads
 
   !--- Local variables
           integer::I, ERR, CURDEV, IA0
 
             ERR = cudaGetDevice(CURDEV)
             blocks  = dim3(mp_BLOCKDIMX, 1, 1)
             threads = dim3(mp_BLOCKSIZE,  1, 1)
             do I=1, m_NDEVICE
                ERR = cudaSetDevice(m_DEVICES(I)) 
                IA0 = WorkSpace%STARTA(I) - 1
                call MarkSeedCell_Kernel0<<<blocks, threads>>>(WorkSpace%NPRT, WorkSpace%NPA(I), IA0, SEED(I)%Data,  &
                                                               WorkSpace%NAPDEV, WorkSpace%IC(I)%Data, MC(I)%Data)
             end do
             ERR = cudaSetDevice(CURDEV)
 
        return
   end subroutine MarkSeedCell_template0
  !****************************************************************************************

  !*******************************************************************************************
  attributes(global) subroutine Active_InCells_Kernel(NPRT, MC, INC, STATU)
  !***  PURPOSE:   Kernel to make atoms in marked cells activative
  !
  !$$   INPUT:     NPRT,      actuall number of atoms on the device
  !$$              MC,        ID of seed CELLS
  !$$              INC,       cell ID of particles, calculated in partitioning      
  !
  !$$   OUTPUT     STATU,     statu of atoms
  !
  implicit none
  !----   DUMMY Variables
          integer, value,  intent(in) :: NPRT
          integer, device, intent(in) :: MC(*)
          integer, device, intent(in) :: INC(*)
          integer, device, intent(out):: STATU(*)

  !----   Local variables
          integer   ::IT, IB, IC, NT, NB, NAL, NL, LL
  !
             IT  = (threadidx%y-1)*blockdim%x + threadidx%x
             !$$--- size of Block
             NT = blockdim%x*blockdim%y
             !$$--- size of grid
             NB = griddim%x*griddim%y

             !$$--- max number of atoms this Kernel can handle
             NAL = NT*NB
             NL = (NPRT-1)/NAL+1

             IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
             do LL=1, NL
                !$$IC -- the id of the atom on the device
                !$$      NOTE: IC is not the id of the same atom in the whole box
                IC= (IB-1)*NT + IT + (LL-1)*NAL

                if(IC .le. NPRT ) then
                   if(MC(INC(IC)) .gt. 0) then
                     STATU(IC)  = ibset(STATU(IC), CP_STATU_ACTIVE_BITPOS)
                   end if
                end if
            end do
        return
  end subroutine Active_InCells_Kernel
  !****************************************************************************************

  !****************************************************************************************
  subroutine Active_InCells_template0(MC, WorkSpace)
  !***  PORPOSE:  template for devices invoking the Kernel
  !
  !     INPUT:     MC,          marked cells in which the atom to be activated
  !
  !     OUTPUT:    WorkSpace,   workspace with STATUS of atoms in cells activated
  !
   implicit none
   !----   DUMMY Variables
           type(DevVec_I),  dimension(:)  ::MC
           type(MDDEVWorkSpace)           ::WorkSpace
 
   !--- Device variables and variables to be used in GPU
           type(dim3) :: blocks
           type(dim3) :: threads
 
   !--- Local variables
           integer::I, ERR, CURDEV
 
             ERR = cudaGetDevice(CURDEV)
             blocks  = dim3(mp_BLOCKDIMX, 1, 1)
             threads = dim3(mp_BLOCKSIZE,  1, 1)
             do I=1, m_NDEVICE
                ERR = cudaSetDevice(m_DEVICES(I)) 
                call Active_InCells_Kernel<<<blocks, threads>>>(WorkSpace%NPA(I),  MC(I)%Data, WorkSpace%IC(I)%Data, &
                                                                WorkSpace%STATU(I)%Data)
             end do
             ERR = cudaSetDevice(CURDEV)
 
        return
   end subroutine Active_InCells_template0
  !****************************************************************************************

  !****************************************************************************************
  subroutine ActiveByCells0(SimBox, CtrlParam, WorkSpace)
  !***  PORPOSE: to active all in-box atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !
   implicit none
   !----   DUMMY Variables
           type(SimMDBox), dimension(:) ::SimBox
           type(SimMDCtrl)              ::CtrlParam
           type(MDDEVWorkSpace)         ::WorkSpace
   !--- Device variables and variables to be used in GPU
           integer, parameter::mp_NNC=27
           integer::NIX(mp_NNC)=(/0,-1,-1,-1, 0, 0, -1, 1,-1, 0, 1,-1, 0, 1, 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/)
           integer::NIY(mp_NNC)=(/0, 0,-1, 1, 1, 0,  0, 0,-1,-1,-1, 1, 1, 1, 0, 1,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/)
           integer::NIZ(mp_NNC)=(/0, 0, 0, 0, 0, 1,  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
             
           integer::I, L, IX, IY, IZ, IX1, IY1, IZ1, NCX, NCY, NCZ, NCXY, NCXYZ, NC, IC0, IC, IB, IN
           integer::PD(3)

           !integer, dimension(:), allocatable::IDebug
           !allocate(IDebug(WorkSpace%NPRT))
   !-----
           !$$--- create the inital seeds
           if(iand(CtrlParam%AR_METHOD, CP_USERDEF_AR) .eq. CP_USERDEF_AR) then
              if(.not.associated(m_pCreatAS) ) then
                 write(*,fmt="(' MDPSCU Error: user defined procedure is needed for creating active region', A16)")
                 write(*,fmt="('               but the procedure is not given', A16)")
                 write(*,fmt="('               check input of &ACTIVEREGSUBCTL subsection', A16)")
                 write(*,fmt="('   Process to be stopped', A16)")
                 stop
              end if
              call m_pCreatAS(SimBox, CtrlParam, hm_ActiveSeed)
              do I=1, m_NDEVICE
                call DevCopyIn(hm_ActiveSeed, dm_ActiveSeed(I), dm_NPRT)
              end do   
     
           else if(iand(CtrlParam%AR_METHOD, CP_LWORD) .eq. 0) then
                !$$--- if no intrinic method of creating seed, we return 
                return
           else  
              call DevSet(dm_ActiveSeed, 0)
              if(iand(CtrlParam%AR_METHOD, CP_CENTPART_AR) .eq. CP_CENTPART_AR)  then
                 call CreateSeedByType1(CtrlParam%AR_CENTPART)
              end if
              if(iand(CtrlParam%AR_METHOD, CP_EKIN_AR) .eq. CP_EKIN_AR) then
                 !--- NOTE: the input value of AR_EKIN is in eV
                 call CreateSeedByEkin1(CtrlParam%AR_EKIN*CP_EVERG)
              end if
           end if

           !$$--- mark the cells that the seeds located in
           call DevSet(dm_ActiveMark, 0)
           call MarkSeedCell_template0(WorkSpace, dm_ActiveSeed, dm_ActiveMark)

           !$$--- mark the neighboring cells. Performing this on CPU is more efficient?
           NC = WorkSpace%NC
           do I=1, m_NDEVICE    
              call DevCopyOut(hm_SwapSeed(:,I), dm_ActiveMark(I), NC)
           end do 
           hm_ActiveSeed(1:NC) = 0
           do I=1, m_NDEVICE
              hm_ActiveSeed(1:NC) = hm_ActiveSeed(1:NC) + hm_SwapSeed(1:NC, I)  
           end do 
           
           NCX   = WorkSpace%NCELL(1)
           NCY   = WorkSpace%NCELL(2)
           NCZ   = WorkSpace%NCELL(3)
           NCXY  = NCY*NCX
           NCXYZ = NCZ*NCXY
           PD    = CtrlParam%IFPD
           hm_SwapSeed(1:NC,1) = 0
           do L=1, CtrlParam%AR_Extend
              IC0 = 0
              do IB=1, WorkSpace%NBOX
                 IC = IC0  
                 do IZ=1, NCZ
                    do IY=1, NCY
                       do IX=1, NCX
                          IC = IC + 1
                          if(hm_ActiveSeed(IC) .gt. 0) then
                             do IN = 1, mp_NNC
                                IX1 = IX + NIX(IN)
                                if(IX1 .gt. NCX )then
                                   if(PD(1) .gt. 0) IX1 = 1
                                else if (IX1 .lt. 1) then
                                   if(PD(1) .gt. 0) IX1 = NCX
                                end if
                                
                                IY1 = IY + NIY(IN)
                                if(IY1 .gt. NCY )then
                                   if(PD(2) .gt. 0) IY1 = 1
                                else if (IY1 .lt. 1) then
                                   if(PD(2) .gt. 0) IY1 = NCY
                                end if

                                IZ1 = IZ + NIZ(IN)
                                if(IZ1 .gt. NCZ )then
                                   if(PD(3) .gt. 0) IZ1 = 1
                                else if (IZ1 .lt. 1) then
                                   if(PD(3) .gt. 0) IZ1 = NCZ
                                end if
                                if(IX1 .lt. 1 .or. IX1 .gt. NCX .or. &
                                   IY1 .lt. 1 .or. IY1 .gt. NCY .or. &
                                   IZ1 .lt. 1 .or. IZ1 .gt. NCZ       ) then
                                    cycle
                                end if    
                                hm_SwapSeed(NCXY*(IZ1-1)+NCX*(IY1-1)+IX1+IC0,1) = 1
                             end do
                          end if 
                       end do  
                    end do 
                 end do  
                 IC0 = IC0 + NCXYZ !--- end loop for boxes    
              end do
              hm_ActiveSeed(1:NC) = hm_SwapSeed(1:NC,1)
               
           end do  

           !--- to active the atoms in marked cells
           do I=1, m_NDEVICE 
               call DevCopyIn(hm_ActiveSeed, dm_ActiveSeed(I), NC)  
           end do     
           call Active_InCells_template0(dm_ActiveSeed, WorkSpace)
           !deallocate(idebug)
       return
   end subroutine ActiveByCells0
  !****************************************************************************************

  !****************************************************************************************
  subroutine ActiveByCells1(SimBox,CtrlParam)
  !***  PORPOSE: to active all in-box atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:) ::SimBox
          type(SimMDCtrl)              ::CtrlParam
  !--- Device variables and variables to be used in GPU
          integer::I, L


          if(iand(CtrlParam%AR_METHOD, CP_KEEP_AR) .ne. CP_KEEP_AR) then
             m_KEEPFLAG = 0
          end if

          !$$--- NOTE: initially , m_KEEPFLAG is set as 0
          if(m_KEEPFLAG .le. 0) then
             call DeActive_All_ActiveRegion0(dm_WorkSpace)
             if(iand(CtrlParam%AR_METHOD, CP_KEEP_AR) .eq. CP_KEEP_AR) then
                m_KEEPFLAG = 1
             end if
          end if

          call ActiveByCells0(SimBox, CtrlParam, dm_WorkSpace)
       return
  end subroutine ActiveByCells1
  !****************************************************************************************   

  !****************************************************************************************
  subroutine ActivateRegion_DEV(SimBox,CtrlParam)
   !***  PORPOSE: to active all in-box atoms
   !
   !     INPUT:     SimBox,    simulation box, not used
   !                CtrlParam, control parameters, not used
   !
   !
   use MD_Globle_Variables_GPU, only:dm_WorkSpace, dm_NPRT, m_STARTA, m_NPA
   use MD_NeighborsList_GPU,    only:dm_Neighbors
   implicit none
   !----   DUMMY Variables
           type(SimMDBox), dimension(:) ::SimBox
           type(SimMDCtrl)              ::CtrlParam
   
           if(iand(CtrlParam%AR_METHOD,CP_ENABLE_AR) .ne. CP_ENABLE_AR) return

           if(ibits(CtrlParam%AR_METHOD, CP_BYNBSETBIT_AR,1) .gt. 0) then
              call ActiveByNeigbors1(SimBox,CtrlParam) 
           else    
              call ActiveByCells1(SimBox,CtrlParam)
           end if  
      return
   end subroutine ActivateRegion_DEV
   !****************************************************************************************  
  end module MD_ActiveRegion_GPU
