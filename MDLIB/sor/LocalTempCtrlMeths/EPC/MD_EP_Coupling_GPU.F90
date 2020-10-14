  module MD_EP_Coupling_GPU
  !**** DESCRIPTION: to calculate EPC model by algorithm of:
  !                  M.W.Finnis etc, Phys. Rev., B44(1991)567
  !
  !****
  !     HISTORY:     2012-03(HOU Qing):
  !                          Updated from  MD_EP_Coupling_GPU.cuf with multpleGPU supported
  !                  2018-03-12(HOU Qing):
  !                           delete the Initialize_EPC_MODIFICATION in the old version
  !                           the PEC physical parameters is passed from CtrlParam
  !                           which defined the parameters for atom groups.
  !
  !                           change the control variables from CONSTANT to DEVICE
  !
  !                  2018-06-12(HOU Qing):
  !                           add calculation of works on atom caused by EPC
  !

  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_MultiGPU_Basic


  implicit none
  !**********************
         integer, parameter,private::mp_BLOCKSIZE=256
         integer, parameter,private::mp_BLOCKDIMX=64

         logical,      private::hm_NEEDDO   = .false.
         logical,      private::hm_SAVEWORK = .false.
         integer,      private::hm_ENABLE(mp_mxGROUP)
         real(KINDDF), private::hm_TE(mp_mxGROUP)
         real(KINDDF), private::hm_EPA(mp_mxGROUP)
         real(KINDDF), private::hm_V2TI(mp_mxGROUP)
         real(KINDDF), private::hm_EPACUT(mp_mxGROUP)
         real(KINDDF), private::hm_EPUPPER(mp_mxGROUP)

         private::EPCoupingDEVWS
         type::EPCoupingDEVWS
               integer,      device, dimension(:), allocatable::ENABLE
               real(KINDDF), device, dimension(:), allocatable::TE
               real(KINDDF), device, dimension(:), allocatable::EPA
               real(KINDDF), device, dimension(:), allocatable::V2TI
               real(KINDDF), device, dimension(:), allocatable::EPACUT
               real(KINDDF), device, dimension(:), allocatable::EPUPPER
               real(KINDDF), device, dimension(:), allocatable::WORK
         end type EPCoupingDEVWS
         type(EPCoupingDEVWS), dimension(m_MXDEVICE), private::dm_EPCoupingDEVWS

         real(KINDDF),         dimension(:), allocatable, private::hm_WORK
         integer,                                         private::m_Init

         private::Allocate_Working_Constants_template,  &
                  Clear_Working_Constants_template,     &
                  Allocate_Working_Constants,           &
                  Clear_Working_Constants,              &
                  Allocate_Working_Varibales_template,  &
                  Clear_Working_Varibales_template,     &
                  Allocate_Working_Varibales,           &
                  Clear_Working_Varibales,              &
                  CopyEPCWORKFrom_Devices_to_Host,      &
                  Copyin_EPC_MOD_template,              &
                  EPC_Mod_KERNEL,                       &
                  EPC_Mod_template,                     &
                  EPC_Mod_Work_KERNEL,                  &
                  EPC_Mod_Work_template
  contains

  !****************************************************************************************
  subroutine Allocate_Working_Constants_template(IDEV, ENABLE, TE, EPA, V2TI, ELR, EUP)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory ENABLE, EPA, V2TI, ELR, EUP
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer,      device, dimension(:), allocatable::ENABLE
           real(KINDDF), device, dimension(:), allocatable::TE
           real(KINDDF), device, dimension(:), allocatable::EPA
           real(KINDDF), device, dimension(:), allocatable::V2TI
           real(KINDDF), device, dimension(:), allocatable::ELR
           real(KINDDF), device, dimension(:), allocatable::EUP

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              allocate(ENABLE(mp_mxGROUP), TE(mp_mxGROUP), EPA(mp_mxGROUP), V2TI(mp_mxGROUP), ELR(mp_mxGROUP), EUP(mp_mxGROUP), STAT=ERR)
              ERR = cudaSetDevice(CURDEV)
              if(ERR) then
    100         write(*,*) "MDPSCU Error in allocating working device memory in EP_Coupling module",IDEV
                stop
              end if
       return
  end subroutine Allocate_Working_Constants_template
  !**************************************************************************

  !****************************************************************************************
  subroutine Clear_Working_Constants_template(IDEV, ENABLE, TE, EPA, V2TI, ELR, EUP)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory ENABLE, EPA, V2TI, ELR, EUP
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer,      device, dimension(:), allocatable::ENABLE
           real(KINDDF), device, dimension(:), allocatable::TE
           real(KINDDF), device, dimension(:), allocatable::EPA
           real(KINDDF), device, dimension(:), allocatable::V2TI
           real(KINDDF), device, dimension(:), allocatable::ELR
           real(KINDDF), device, dimension(:), allocatable::EUP

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              if(allocated(ENABLE)) deallocate(ENABLE, STAT=ERR)
              if(ERR) goto 100
              if(allocated(TE))     deallocate(TE, STAT=ERR)
              if(ERR) goto 100
              if(allocated(EPA))    deallocate(EPA, STAT=ERR)
              if(ERR) goto 100
              if(allocated(V2TI))   deallocate(V2TI, STAT=ERR)
              if(ERR) goto 100
              if(allocated(ELR))    deallocate(ELR, STAT=ERR)
              if(ERR) goto 100
              if(allocated(EUP))    deallocate(EUP, STAT=ERR)
              if(ERR) goto 100
              ERR = cudaSetDevice(CURDEV)
              return

    100       write(*,*) "MDPSCU Error in deallocating working device memory in EP_Coupling module",IDEV
              stop
       return
  end subroutine Clear_Working_Constants_template
  !**************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Constants()
  !***  PURPOSE:   to allocate working space
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::I

               
             !$$--- allocate working space on devices
              do I = 1, m_NDEVICE
                 call Allocate_Working_Constants_template(m_DEVICES(I), dm_EPCoupingDEVWS(I)%ENABLE, dm_EPCoupingDEVWS(I)%TE,  &
                                                                        dm_EPCoupingDEVWS(I)%EPA,    dm_EPCoupingDEVWS(I)%V2TI,& 
                                                                        dm_EPCoupingDEVWS(I)%EPACUT, dm_EPCoupingDEVWS(I)%EPUPPER)
              end do

       return
  end subroutine Allocate_Working_Constants
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_Working_Constants()
  !***  PURPOSE:   to allocate working space
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer::I


             !$$--- deallocate working space on devices
              do I = 1, m_NDEVICE
                 call Clear_Working_Constants_template(m_DEVICES(I), dm_EPCoupingDEVWS(I)%ENABLE, dm_EPCoupingDEVWS(I)%TE,  &
                                                                     dm_EPCoupingDEVWS(I)%EPA,    dm_EPCoupingDEVWS(I)%V2TI,& 
                                                                     dm_EPCoupingDEVWS(I)%EPACUT, dm_EPCoupingDEVWS(I)%EPUPPER)
              end do

              hm_NEEDDO   = .false.
              hm_SAVEWORK = .false.
       return
  end subroutine Clear_Working_Constants
  !****************************************************************************

  !**************************************************************************
  subroutine Allocate_Working_Varibales_template(Idev, Nprt, W)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory ELOSS
  !
      implicit none
      !--- dummy variables
           integer,                             intent(in)::Idev
           integer,                             intent(in)::Nprt
           real(KINDDF), device, dimension(:), allocatable::W

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              allocate(W(Nprt), STAT=ERR)
              ERR = cudaSetDevice(CURDEV)
              if(ERR) then
    100         write(*,*) "MDPSCU Error in allocating working device memory in ST_Coupling module for constant",IDEV
                stop
              end if
       return
  end subroutine Allocate_Working_Varibales_template
  !**************************************************************************

  !**************************************************************************
  subroutine Clear_Working_Varibales_template(IDev, W)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory Eloss
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:),   allocatable::W

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              ERR = 0
              if(allocated(W)) deallocate(W, STAT=ERR)
              if(ERR) then
                 write(*,*) "MDPSCU Error in deallocating working device memory in ST_Coupling module",IDEV
                 stop
              end if
              ERR = cudaSetDevice(CURDEV)
              return
  end subroutine Clear_Working_Varibales_template
  !**************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Varibales()
  !***  PURPOSE:   to allocate working space
  !
    use MD_Globle_Variables_GPU, only:m_NAPDEV, dm_NPRT
      implicit none
      !--- dummy variables
      !--- Local vairables
      integer:: I

             !$$--- allocate working space on device 1
              do I = 1, m_NDEVICE
                 call Allocate_Working_Varibales_template(m_DEVICES(I), m_NAPDEV, dm_EPCoupingDEVWS(I)%WORK )
              end do
              allocate(hm_WORK(dm_NPRT))
              m_INIT = 1
       return
  end subroutine Allocate_Working_Varibales
 !****************************************************************************

 !****************************************************************************
  subroutine Clear_Working_Varibales()
  !***  PURPOSE:   to allocate working space
  !
      implicit none
      !--- dummy variables
      !--- Local vairables
            integer:: I

             !$$--- allocate working space on device 1
              do I = 1, m_NDEVICE
                 call Clear_Working_Varibales_template(m_DEVICES(I), dm_EPCoupingDEVWS(I)%WORK )
              end do
              if(allocated(hm_WORK))deallocate(hm_WORK)

              m_INIT = 0
       return
  end subroutine Clear_Working_Varibales
  !****************************************************************************

  !**************************************************************************
  subroutine Initialize_EPCMOD_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to initialze this module by loading calculation parameters
  !
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  implicit none
      !----   DUMMY Variables
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam

         call Clear_EPCMOD_DEV()
         call Allocate_Working_Constants()
         call Reset_EPCMOD_DEV(SimBox, CtrlParam)
        return
   end subroutine Initialize_EPCMOD_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Clear_EPCMOD_DEV()
  !***  PURPOSE: to clear the memories allocated on calling initialize
  !
  !     INPUT:
  !     OUTPUT:
  !
  implicit none
      !----   DUMMY Variables
         call Clear_Working_Constants()
         call Clear_Working_Varibales()
        return
   end subroutine Clear_EPCMOD_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Copyin_EPC_MOD_template(IDEV, ENABLE, TE, EPA, V2TI, EPACUT, EUP)
  !***  PURPOSE: to reset EPC controal parameters
  !
  !     INPUT:
  !
  !     OUTPUT:
  !
  implicit none
      !--- dummy variables
           integer,      intent(in)::IDEV
           integer,      device, dimension(:)::ENABLE
           real(KINDDF), device, dimension(:)::TE
           real(KINDDF), device, dimension(:)::EPA
           real(KINDDF), device, dimension(:)::V2TI
           real(KINDDF), device, dimension(:)::EPACUT
           real(KINDDF), device, dimension(:)::EUP

      !--- Local vairables
          integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              ERR = cudaMemcpyAsync(ENABLE, hm_ENABLE, mp_mxGROUP)
              ERR = cudaMemcpyAsync(TE,     hm_TE,     mp_mxGROUP)
              ERR = cudaMemcpyAsync(EPA,    hm_EPA,    mp_mxGROUP)
              ERR = cudaMemcpyAsync(V2TI,   hm_V2TI,   mp_mxGROUP)
              ERR = cudaMemcpyAsync(EPACUT, hm_EPACUT, mp_mxGROUP)
              ERR = cudaMemcpyAsync(EUP,    hm_EPUPPER,mp_mxGROUP)

              ERR = cudaSetDevice(CURDEV)
              return
   end subroutine Copyin_EPC_MOD_template
  !**************************************************************************

  !**************************************************************************
  subroutine Reset_EPCMOD_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to reset EPC controal parameters
  !
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  use MD_Globle_Variables_GPU, only:m_NDEVICE, m_DEVICES
  implicit none
      !----   DUMMY Variables
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam
       !--- local
       integer::I


          do I=1, SimBox%NGROUP
              hm_ENABLE(I) = iand(CtrlParam%LT_CTRL(I)%METH, CP_TICTRL_METH_EPC)
              hm_TE(I)     = CtrlParam%LT_CTRL(I)%TI
              hm_V2TI(I)   = SimBox%CM(I)*C_UTH/CP_KB
              hm_EPA(I)    = SimBox%CM(I)/CtrlParam%LT_CTRL(I)%EPC_Alpha
              hm_EPACUT(I) = hm_TE(I)*CtrlParam%LT_CTRL(I)%EPC_CUT
              hm_EPUPPER(I)= 2.D0*CtrlParam%LT_CTRL(I)%EPC_HE/SimBox%CM(I)
          end do

           !$$--- copyin the control parameters
           do I=1, m_NDEVICE 
              call Copyin_EPC_MOD_template(m_DEVICES(I), dm_EPCoupingDEVWS(I)%ENABLE, dm_EPCoupingDEVWS(I)%TE,   &
                                                         dm_EPCoupingDEVWS(I)%EPA,    dm_EPCoupingDEVWS(I)%V2TI, &
                                                         dm_EPCoupingDEVWS(I)%EPACUT, dm_EPCoupingDEVWS(I)%EPUPPER)
           end do

           if(all(hm_ENABLE .eq. 0)) then
              hm_NEEDDO = .false.
              call Clear_Working_Varibales()
           else
              hm_NEEDDO  = .true.
              if(any(CtrlParam%LT_CTRL(1:SimBox%NGROUP)%EPC_SAVEWK .gt. 0) ) then
                 hm_SAVEWORK = .true.
              else
                 hm_SAVEWORK = .false.
                 call Clear_Working_Varibales()
              endif
           end if

        return
   end subroutine Reset_EPCMOD_DEV
  !**************************************************************************

  !**************************************************************************
  attributes(global) subroutine EPC_MOD_KERNEL(IM, IA0, NAPDEV, NPART, ENABLE, TE, &
                                               EPA, V2TI, TCUT, EUP, ITYP, STATU, XP1, FP)
  !***  PURPOSE:   KERNEL to modify the force by E-P coupling
  !
  !$$   INPUT:     IM,        the total number of atoms in the whole box
  !$$              IA0,       the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$
  !$$              ENABLE,    variable indicating if EPC calculation is enabled for types of atom types
  !$$              TE,        electron temperature
  !$$              EPA,       coupling coefficient
  !$$              ITYP,      indicator of type of atoms
  !$$              STATU,     statu of atoms
  !$$              XP1,       velocity of atoms
  !$$              FP,        force on atoms
  !
  !$$   OUTPUT     FP,        updated force
  !
  implicit none
  !----   DUMMY Variables
          integer, value::IM, IA0, NAPDEV, NPART

          integer,     device::ENABLE(*)
          real(KINDDF),device::TE(*), EPA(*), V2TI(*), TCUT(*), EUP(*)
          integer,     device::ITYP(IM)
          integer,     device::STATU(NAPDEV)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)

  !----   Local variables
          integer::KK, IT, IB, IC
          real(KINDDF)::XP1SQ, TM,  MU, VX,VY,VZ
          real(KINDDF), shared::SEUP(mp_MXGROUP), SEPA(mp_MXGROUP), SV2TI(mp_MXGROUP), STE(mp_MXGROUP), STCUT(mp_MXGROUP)
          integer,      shared::FLAG(mp_mxGROUP)
  !

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT
              if(IT .le. mp_MXGROUP) then
                 SEUP(IT)   = EUP(IT)
                 SEPA(IT)   = EPA(IT)
                 SV2TI(IT)  = V2TI(IT)
                 STE(IT)    = TE(IT)
                 STCUT(IT)  = TCUT(IT)
                 FLAG(IT)   = ENABLE(IT)
              end if
              call syncthreads()

              if(IC.LE.NPART) then
                 KK = ITYP(IC+IA0)
                 if(FLAG(KK) .gt. 0) then
                   if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE ) then
                     !TECUT = dcm_EPACUT(KK)*TE
                      VX = XP1(IC,1)
                      VY = XP1(IC,2)
                      VZ = XP1(IC,3)
                      XP1SQ = VX*VX+VY*VY+VZ*VZ
                      if(XP1SQ .le. SEUP(KK) ) then
                        !$$-- ion temperature
                         TM = XP1SQ*SV2TI(KK)  !dcm_CM(KK)*C_UTH/CP_KB
                         MU = SEPA(KK)*(TM-STE(KK))/max(TM, STCUT(KK))

                         FP(IC,1) = FP(IC,1) - MU*VX
                         FP(IC,2) = FP(IC,2) - MU*VY
                         FP(IC,3) = FP(IC,3) - MU*VZ
                      end if
                   end if
                 end if
             end if
        return
  end subroutine EPC_MOD_KERNEL
  !****************************************************************************************

  !****************************************************************************************
  subroutine EPC_MOD_template(IDEV, STARTCELL, ENDCELL, ENABLE, TE, &
                                               EPA, V2TI, EPACUT, EUP, ITYP, STATU, XP1, FP)
  !***  PORPOSE:  template for devices invoking the KERNEL
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                ENABLE,    variable indicating if EPC calculation is enabled for types of atom types
  !                TE,        assumed electron temperature
  !                EPA,       coupling coefficient
  !                ITYP,      type of atoms
  !                XP1,       velocity of atoms
  !                FP,        force on the atoms
  !
  !     OUTPUT     FP,        updated force
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)::IDEV, STARTCELL, ENDCELL
          integer,      device, dimension(:)::ENABLE
          real(KINDDF), device, dimension(:)::TE, EPA, V2TI, EPACUT, EUP

          integer,      device, dimension(:)::ITYP, STATU
          real(KINDDF), device, dimension(:,:)::XP1, FP

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT     = ENDA - STARTA + 1

           !$$--- to determine size of a block (the number of threads in a block)
            BX      = mp_BLOCKSIZE
            BY      = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call EPC_MOD_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, m_NAPDEV, NPRT,          &
                                           ENABLE, TE, EPA, V2TI, EPACUT, EUP, ITYP, STATU, XP1,FP)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine EPC_MOD_template
  !****************************************************************************************

  !**************************************************************************
  attributes(global) subroutine EPC_MOD_WORK_KERNEL(IM, IA0, NAPDEV, NPART, ENABLE, TE, &
                                               EPA, V2TI, TCUT, EUP, ITYP, STATU, XP1, FP, DT, WORK)
  !***  PURPOSE:   KERNEL to modify the force by E-P coupling
  !
  !$$   INPUT:     IM,        the total number of atoms in the whole box
  !$$              IA0,       the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$
  !$$              ENABLE,    variable indicating if EPC calculation is enabled for types of atom types
  !$$              TE,        electron temperature
  !$$              EPA,       coupling coefficient
  !$$              ITYP,      indicator of type of atoms
  !$$              STATU,     statu of atoms
  !$$              XP1,       velocity of atoms
  !$$              FP,        force on atoms
  !$$              DT,        the size of time step
  !
  !$$   OUTPUT     FP,        updated force
  !$$              WORK,      works done by EPC
  !
  implicit none
  !----   DUMMY Variables
          integer, value::IM, IA0, NAPDEV, NPART
          real(KINDDF),value::DT

          integer,     device::ENABLE(*)
          real(KINDDF),device::TE(*), EPA(*), V2TI(*), TCUT(*), EUP(*)
          integer,     device::ITYP(IM)
          integer,     device::STATU(NAPDEV)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)
          real(KINDDF),device::WORK(NAPDEV)

  !----   Local variables
          integer::KK, IT, IB, IC
          real(KINDDF)::XP1SQ, TM,  MU, WW, VX,VY,VZ
          real(KINDDF), shared::SEUP(mp_MXGROUP), SEPA(mp_MXGROUP), SV2TI(mp_MXGROUP), STE(mp_MXGROUP), STCUT(mp_MXGROUP)
          integer,      shared::FLAG(mp_mxGROUP)
  !

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT
              if(IT .le. mp_MXGROUP) then
                 SEUP(IT)   = EUP(IT)
                 SEPA(IT)   = EPA(IT)
                 SV2TI(IT)  = V2TI(IT)
                 STE(IT)    = TE(IT)
                 STCUT(IT)  = TCUT(IT)
                 FLAG(IT)   = ENABLE(IT)
              end if
              call syncthreads()


              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC.LE.NPART) then
                 KK = ITYP(IC+IA0)
                 WW = 0.D0
                 if(FLAG(KK) .gt. 0) then
                    if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                     !TECUT = dcm_EPACUT(KK)*TE
                      VX = XP1(IC,1)
                      VY = XP1(IC,2)
                      VZ = XP1(IC,3)
                      XP1SQ = VX*VX+VY*VY+VZ*VZ
                      if(XP1SQ .le. SEUP(KK) ) then
                        !$$-- ion temperature
                         TM = XP1SQ*SV2TI(KK)  !dcm_CM(KK)*C_UTH/CP_KB
                         MU = SEPA(KK)*(TM-STE(KK))/max(TM, STCUT(KK))

                         FP(IC,1) = FP(IC,1) - MU*VX
                         FP(IC,2) = FP(IC,2) - MU*VY
                         FP(IC,3) = FP(IC,3) - MU*VZ
                         WW       = MU*dsqrt(XP1SQ)*DT
                      end if
                    end if
                 end if
                 WORK(IC) = WW
              end if
        return
  end subroutine EPC_MOD_WORK_KERNEL
  !****************************************************************************************

  !****************************************************************************************
  subroutine EPC_MOD_WORK_template(IDEV, STARTCELL, ENDCELL, ENABLE, TE, &
                                               EPA, V2TI, EPACUT, EUP, ITYP, STATU, XP1, FP, DT, WORK)
  !***  PORPOSE:  template for devices invoking the KERNEL
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                ENABLE,    variable indicating if EPC calculation is enabled for types of atom types
  !                TE,        assumed electron temperature
  !                EPA,       coupling coefficient
  !                ITYP,      type of atoms
  !                XP1,       velocity of atoms
  !                FP,        force on the atoms
  !                DT,        the size of time step
  !
  !     OUTPUT     FP,        updated force
  !                WORK,      works done by EPC
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)::IDEV, STARTCELL, ENDCELL
          integer,      device, dimension(:)::ENABLE
          real(KINDDF), device, dimension(:)::TE, EPA, V2TI, EPACUT, EUP

          integer,      device, dimension(:)  ::ITYP, STATU
          real(KINDDF), device, dimension(:,:)::XP1, FP
          real(KINDDF), intent(in)            ::DT
          real(KINDDF), device, dimension(:)  ::WORK

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT     = ENDA - STARTA + 1

           !$$--- to determine size of a block (the number of threads in a block)
            BX      = mp_BLOCKSIZE
            BY      = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call EPC_MOD_WORK_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, m_NAPDEV, NPRT,   &
                                           ENABLE, TE, EPA, V2TI, EPACUT, EUP, ITYP, STATU, XP1, FP, DT, WORK)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine EPC_MOD_WORK_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_EPCMOD0_DEV(SimBox,CtrlParam)
  !***  PORPOSE: to calculate the friction force of PEC
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox),dimension(:)::SimBox
          type(SimMDCtrl)            ::CtrlParam
  !--- Device variables and variables to be used in GPU
          integer::I

          if(.not.hm_NEEDDO) return

          if(m_Init .le. 0) then
             call Allocate_Working_Varibales()
          end if

          do I=1, m_NDEVICE
              call EPC_MOD_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                  &
                                    dm_EPCoupingDEVWS(I)%ENABLE,  dm_EPCoupingDEVWS(I)%TE,       &
                                    dm_EPCoupingDEVWS(I)%EPA,     dm_EPCoupingDEVWS(I)%V2TI,     &
                                    dm_EPCoupingDEVWS(I)%EPACUT,  dm_EPCoupingDEVWS(I)%EPUPPER,  &
                                    dm_WorkSpace%ITYP(I)%Data,    dm_WorkSpace%STATU(I)%Data,    &
                                    dm_WorkSpace%XP1(I)%Data,     dm_WorkSpace%FP(I)%Data)
          end do

     return
  end subroutine Do_EPCMOD0_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_EPCMOD1_DEV(SimBox,CtrlParam)
  !***  PORPOSE: to calculate the friction force of PEC and also the word done
  !              by EPC
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox),dimension(:)::SimBox
          type(SimMDCtrl)            ::CtrlParam
  !--- Device variables and variables to be used in GPU
          real(KINDDF)::DT
          integer::I

          if(.not.hm_NEEDDO) return

          !--- NOTE: the time step CtrlParam%H could be varing
          !          ref: MD_DiffScheme_GPU.F90
          DT = CtrlParam%H
          do I=1, m_NDEVICE
              call EPC_MOD_WORK_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),            &
                                    dm_EPCoupingDEVWS(I)%ENABLE,  dm_EPCoupingDEVWS(I)%TE,      &
                                    dm_EPCoupingDEVWS(I)%EPA,     dm_EPCoupingDEVWS(I)%V2TI,    &
                                    dm_EPCoupingDEVWS(I)%EPACUT,  dm_EPCoupingDEVWS(I)%EPUPPER, &
                                    dm_WorkSpace%ITYP(I)%Data,    dm_WorkSpace%STATU(I)%Data,   &
                                    dm_WorkSpace%XP1(I)%Data,     dm_WorkSpace%FP(I)%Data,      &
                                    DT, dm_EPCoupingDEVWS(I)%WORK                                )
          end do

     return
  end subroutine Do_EPCMOD1_DEV
  !****************************************************************************************
  !****************************************************************************************
  subroutine Do_EPCMOD_DEV(SimBox,CtrlParam)
  !***  PORPOSE: to interface to host process
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  use MD_Globle_Variables_GPU, only:hm_GIDINV
  implicit none
  !----   DUMMY Variables
          type(SimMDBox),dimension(:)::SimBox
          type(SimMDCtrl)            ::CtrlParam
  !--- Device variables and variables to be used in GPU
          integer::I, IA1, IA2

          if(.not.hm_NEEDDO) return

          if(hm_SAVEWORK) then
             call Do_EPCMOD1_DEV(SimBox, CtrlParam)
             call CopyEPCWORKFrom_Devices_to_Host()
             IA1 = 1
             IA2 = SimBox(1)%NPRT
             do I=1, size(SimBox)
                !$$--- NOTE: we use the whole array hm_ELOSS, but not its subset(IA1:IA2)
                call AccumDatPad_SimMDBox(SimBox(I), "EPC-WORK (ev)", hm_WORK*CP_ERGEV, &
                                               hm_GIDINV(IA1:IA2))
                IA1 = IA2 + 1
                IA2 = IA2 + SimBox(I)%NPRT
             end do
          else
             call Do_EPCMOD0_DEV(SimBox, CtrlParam)
          end if

       return
  end subroutine Do_EPCMOD_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine CopyEPCWORKFrom_Devices_to_Host()
  !***  PURPOSE:  copy Eloss on devices to the arries on host.
  !
    use MD_Globle_Variables_GPU, only:m_NDEVICE, COPY_OUT_SHIFT_template
    implicit none
      !--- dummy variables
      !----   Local variables
       integer::I

            do I=1, m_NDEVICE
               call COPY_OUT_SHIFT_template(I, dm_EPCoupingDEVWS(I)%WORK, hm_WORK)
            end do

            return
   end subroutine CopyEPCWORKFrom_Devices_to_Host
  !****************************************************************************


  end module MD_EP_Coupling_GPU
