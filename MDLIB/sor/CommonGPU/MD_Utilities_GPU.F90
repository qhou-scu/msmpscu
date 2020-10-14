  module MD_Utilities_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                 Provide a number routines that use the GPU globale vairables to calculate
  !                 commonly used quantities.
  !    DEPENDENCE:  _______________________________________________________________________________________
  !                 MD_NeighborsList_GPU.F90
  !                 MD_Globle_Variables_GPU.F90
  !                 MD_ForceClass_Register_GPU
  !
  !    REFERENCE:   _______________________________________________________________________________________
  !                 HOU Qing, Nov, 2013
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl

  !-----------------------------------------------
    implicit none
  !---------------------------------------------------------
      private::                               & 
                Cal_EKinTensor_KERNEL,        &
                Cal_EKinTensor_template,      &
                Cal_Atomic_KineticTensor_OLD, &
                Cal_Atomic_KineticTensor_0
      public::  Cal_Atomic_KineticTensor_DEV
      interface Cal_Atomic_KineticTensor_DEV
          module procedure Cal_Atomic_KineticTensor_OLD
          module procedure Cal_Atomic_KineticTensor_0
      end interface Cal_Atomic_KineticTensor_DEV

  contains
  !*********************************************************************************
  attributes(global) subroutine Cal_EKinTensor_KERNEL(IM,IA0,ITYP, NAPDEV,NPRT, XP1, STATU, CM, AKP)
  !***  PURPOSE:   to calculate the kinetic energies tensor of atoms
  !
  !     INPUT:     IM,         total number of atoms in the whole box
  !                IA0,        starting atom on the device
  !                ITYP,       indicator of type of atoms
  !                NAPDEV,     permitted number of atoms on the device
  !                NPRT,       actuall number of atoms on the device
  !                XP1,        velocity of atoms
  !                STATU,      statu ot atoms
  !
  !     OUTPUT     AKP,        kinetic tensor of atoms
  !

  !
  implicit none
  !----   DUMMY Variables
          integer, value, intent(in)::IM, NAPDEV,NPRT, IA0
          integer,device, intent(in)::ITYP(IM)
          integer,device, intent(in)::STATU(NAPDEV)
          real(KINDDF), device, intent(in)::XP1(NAPDEV,*),CM(*)
          real(KINDDF), device, intent(out)::AKP(NAPDEV,*)

  !----   Local variables
          integer::NT, NB, NAL, NL, LL, IT, IB, IC,ITY
          real(KINDDF)::XP1x,XP1y, XP1z, CM0
          real(KINDDF)::KET11, KET12, KET13
          real(KINDDF)::KET21, KET22, KET23
          real(KINDDF)::KET31, KET32, KET33
  !***
  !
             !$$--- size of Block
              NT = blockdim%x*blockdim%y

             !$$--- size of grid
              NB = griddim%x*griddim%y

             !$$--- how many loop needed
              NAL = NT*NB
              NL = (NPRT-1)/NAL+1

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              do LL=1, NL
                 !$$IC -- the id of the atom on the device
                 IC= (IB-1)*NT+IT + (LL-1)*NAL
                 if(IC .LE. NPRT) then
                    KET11 = C_ZERO
                    KET12 = C_ZERO
                    KET13 = C_ZERO
                    KET21 = C_ZERO
                    KET22 = C_ZERO
                    KET23 = C_ZERO
                    KET31 = C_ZERO
                    KET32 = C_ZERO
                    KET33 = C_ZERO
                    if(IAND(STATU(IC),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then

                       ITY =  ITYP(IC+IA0)
                       CM0 = CM(ITY)      !--- NOTE: donot multiply C_HALF
                       XP1x = XP1(IC,1)
                       XP1y = XP1(IC,2)
                       XP1z = XP1(IC,3)
                       KET11 = XP1x*XP1x*CM0
                       KET12 = XP1x*XP1y*CM0
                       KET13 = XP1x*XP1z*CM0
                       KET21 = XP1y*XP1x*CM0
                       KET22 = XP1y*XP1y*CM0
                       KET23 = XP1y*XP1z*CM0
                       KET31 = XP1z*XP1x*CM0
                       KET32 = XP1z*XP1y*CM0
                       KET33 = XP1z*XP1z*CM0
                    end if
                    AKP(IC,1) = KET11
                    AKP(IC,2) = KET12
                    AKP(IC,3) = KET13
                    AKP(IC,4) = KET21
                    AKP(IC,5) = KET22
                    AKP(IC,6) = KET23
                    AKP(IC,7) = KET31
                    AKP(IC,8) = KET32
                    AKP(IC,9) = KET33
                 end if
               end do
         return
  end subroutine Cal_EKinTensor_KERNEL
  !******************************************************************************************

  !****************************************************************************************
  SUBROUTINE Cal_EKinTensor_template(IDEV, NPRT, NAPDEV, STARTA, NPA, CM, ITYP, XP1, STATU, AKP)
  !***  PORPOSE:   template to invoke CALEKIN_KERNEL on devices
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                CM,         mass of atoms
  !                ITYP,       type id of atoms
  !                XP1,        velocities of atoms on device IDEV
  !                STATU,      statu of atoms on  device IDEV

  !     OUTPUT     EKIN,       kinetic energies of atoms on device IDEV
  !
  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
           integer,              intent(in)    ::IDEV, NPRT, NAPDEV, STARTA, NPA
           integer,      device, dimension(:)  ::ITYP, STATU
           real(KINDDF), device, dimension(:,:)::XP1
           real(KINDDF), device, dimension(:)  ::CM
           real(KINDDF), device, dimension(:,:)::AKP

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads
          integer, parameter::mp_BLOCKSIZE=256

  !--- Local variables
           integer::IA0,  BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            !$$--- to determine size of a block (the number of threads in a block)
            BX      = mp_BLOCKSIZE
            BY      = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            IA0     = STARTA-1
            call CAL_EKINTENSOR_KERNEL<<<blocks, threads>>>(NPRT, IA0, ITYP, NAPDEV, NPA, XP1, STATU, CM, AKP)
            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine Cal_EKinTensor_template
  !********************************************************************************

  !*********************************************************************************
  subroutine Cal_Atomic_KineticTensor_Old(IDEV, dCM, dAKP)
  !***  PURPOSE:   to calculate the atomic kinetic tensor ondevice IDEV
  !
  !     INPUT:     IDEV:      the id of device
  !                dCM:       the mas of particles
  !     OUTPUT     dxm_AKP:   the atomic kinetic energy
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV
      real(KINDDF), device, dimension(:)   ::dCM
      real(KINDDF), device, dimension(:,:) ::dAKP

      !--- Local variables
               call Cal_EKinTensor_template(m_DEVICES(IDEV), dm_NPRT, m_NAPDEV, m_STARTA(IDEV), m_NPA(IDEV), dCM, &
                                            dm_WorkSpace%ITYP(IDEV)%Data, dm_WorkSpace%XP1(IDEV)%Data,  &
                                            dm_WorkSpace%STATU(IDEV)%Data,                              & 
                                            dAKP)
              return
  end subroutine Cal_Atomic_KineticTensor_Old
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_Atomic_KineticTensor_0(WorkSpace, dAKP)
  !***  PURPOSE:   to calculate the atomic kinetic tensor ondevice IDEV
  !
  !     INPUT:     IDEV:      the id of device
  !                dCM:       the mas of particles
  !     OUTPUT     dxm_AKP:   the atomic kinetic energy
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      type(MDDEVWorkSpace),  intent(in)  ::WorkSpace
      type(DevMat_DF),       intent(out) ::dAKP(:)

      !--- Local variables
        integer::I 

               do I=1, m_NAPDEV
                  call Cal_EKinTensor_template(m_DEVICES(I),                                       &
                                               WorkSpace%NPRT,          WorkSpace%NAPDEV,          &
                                               WorkSpace%STARTA(I),     WorkSpace%NPA(I),          &
                                               WorkSpace%CM(I)%Data,    WorkSpace%ITYP(I)%Data,    &
                                               WorkSpace%XP1(I)%Data,   WorkSpace%STATU(I)%Data,   & 
                                               dAKP(I)%Data)
               end do                                 
              return
  end subroutine Cal_Atomic_KineticTensor_0
  !*********************************************************************************


  end module MD_Utilities_GPU





