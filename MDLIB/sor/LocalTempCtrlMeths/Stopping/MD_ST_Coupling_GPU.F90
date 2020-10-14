  module MD_ST_Coupling_GPU
  !**** DESCRIPTION: to calculate stopping force and inelastic energy loss.
  !                  The stopping force are treated as a kind of friction
  !                  force without changing the moving direction of atom,
  !                  and are added to the force calculated by potential in
  !                  generic force calculations.
  !
  !****
  !     HISTORY:     first version 2018-03(HOU Qing):
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_MultiGPU_Basic

  implicit none
  !**********************
         integer, parameter,private::mp_BLOCKSIZE=256

         integer,parameter,private::mp_STMOD_G   = 1
         integer,parameter,private::mp_STMOD_L   = 2
         logical,          private::hm_NEEDDO    = .false.
         logical,          private::hm_SAVEELOSS = .false.
         integer,          private::hm_STMOD     = mp_STMOD_L


         integer,      private::hm_ENABLE(mp_mxGROUP)
         integer,      private::hm_DENMETH                          ! the style of calculating media density
         real(KINDDF), private::hm_MDEN(mp_mxGROUP)                 ! the media density

         integer,      private::m_NE, m_NK

         private::STCoupingDEVWS
         type::STCoupingDEVWS
               integer,      device, dimension(:),   allocatable::ENABLE
               integer,      device, dimension(:,:), allocatable::KPAIR
               real(KINDDF), device, dimension(:),   allocatable::CM2
               real(KINDDF), device, dimension(:),   allocatable::MDEN
               real(KINDDF), device, dimension(:,:), allocatable::LVOL
               real(KINDDF), device, dimension(:),   allocatable::ETAB
               real(KINDDF), device, dimension(:,:), allocatable::STAB
               real(KINDDF), device, dimension(:),   allocatable::ELOSS
         end type STCoupingDEVWS
         type(STCoupingDEVWS), dimension(m_MXDEVICE), private::dm_STCoupingDEVWS

         real(KINDDF),         dimension(:),   allocatable, private::hm_ELOSS
         integer,                                           private::m_INIT = 0

         private::Allocate_Working_Constants_template,  &
                  Clear_Working_Constants_template,     &
                  Allocate_Working_Constants,           &
                  Clear_Working_Constants,              &
                  Allocate_Working_Varibales_template,  &
                  Clear_Working_Varibales_template,     &
                  Allocate_Working_Varibales,           &
                  Clear_Working_Varibales,              &
                  CopyElossFrom_Devices_to_Host,        &
                  ST_MOD_GDEN_KERNEL,                   &
                  ST_MOD_GDEN_template,                 &
                  ST_MOD_LDEN_KERNEL,                   &
                  ST_MOD_LDEN_template,                 &
                  ST_MOD_ELOSS_GDEN_KERNEL,             &
                  ST_MOD_ELOSS_GDEN_template,           &
                  ST_MOD_ELOSS_LDEN_KERNEL,             &
                  ST_MOD_ELOSS_LDEN_template
  contains

  !****************************************************************************************
  subroutine Allocate_Working_Constants_template(IDEV, STPTable, ENABLE, KPAIR, MDEN, CM2, &
                                                 LVOL, ETAB, STAB)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory ENABLE, EPA, V2TI, ELR, EUP
  !
      use MD_TYPEDEF_StopRange_Table
      implicit none
      !--- dummy variables
           integer,          intent(in)::IDEV
           type(MDSTPTable), intent(in)::STPTable

           integer,      device, dimension(:),   allocatable::ENABLE
           integer,      device, dimension(:,:), allocatable::KPAIR
           real(KINDDF), device, dimension(:),   allocatable::MDEN
           real(KINDDF), device, dimension(:),   allocatable::CM2
           real(KINDDF), device, dimension(:,:), allocatable::LVOL
           real(KINDDF), device, dimension(:),   allocatable::ETAB
           real(KINDDF), device, dimension(:,:), allocatable::STAB

      !--- Local vairables
      integer::CURDEV, ERR, NG

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              NG    = size(STPTable%KPAIR, dim=1)
              m_NE  = size(STPTable%STPWR, dim=1)
              m_NK  = size(STPTable%STPWR, dim=2)
              allocate(ENABLE(mp_mxGROUP), CM2(mp_mxGROUP), MDEN(mp_mxGROUP), LVOL(mp_mxGROUP,mp_mxGROUP), &
                       KPAIR(NG,NG), ETAB(m_NE), STAB(m_NE,m_NK), STAT=ERR)

              KPAIR  = STPTable%KPAIR
              ETAB   = STPTable%E
              STAB   = STPTable%STPWR
              ERR = cudaSetDevice(CURDEV)
              if(ERR) then
    100         write(*,*) "MDPSCU Error in allocating working device memory in ST_Coupling module for constant",IDEV
                stop
              end if
       return
  end subroutine Allocate_Working_Constants_template
  !********************************************************************************

  !********************************************************************************
  subroutine Clear_Working_Constants_template(IDEV, ENABLE, KPAIR, MDEN, CM2, LVOL, ETAB, STAB)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory ENABLE, EPA, V2TI, ELR, EUP
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer,      device, dimension(:),   allocatable::ENABLE
           integer,      device, dimension(:,:), allocatable::KPAIR
           real(KINDDF), device, dimension(:),   allocatable::MDEN
           real(KINDDF), device, dimension(:,:), allocatable::CM2
           real(KINDDF), device, dimension(:,:), allocatable::LVOL
           real(KINDDF), device, dimension(:),   allocatable::ETAB
           real(KINDDF), device, dimension(:,:), allocatable::STAB

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              if(allocated(ENABLE))    deallocate(ENABLE, STAT=ERR)
              if(ERR) goto 100
              if(allocated(KPAIR))     deallocate(KPAIR, STAT=ERR)
              if(ERR) goto 100
              if(allocated(MDEN))      deallocate(MDEN, STAT=ERR)
              if(ERR) goto 100
              if(allocated(CM2))       deallocate(CM2, STAT=ERR)
              if(ERR) goto 100
              if(allocated(LVOL))     deallocate(LVOL, STAT=ERR)
              if(ERR) goto 100
              if(allocated(ETAB))      deallocate(ETAB, STAT=ERR)
              if(ERR) goto 100
              if(allocated(STAB))      deallocate(STAB, STAT=ERR)
              if(ERR) goto 100

              ERR = cudaSetDevice(CURDEV)
              return

    100       write(*,*) "MDPSCU Error in deallocating working device memory in ST_Coupling module",IDEV
              stop
       return
  end subroutine Clear_Working_Constants_template
  !**************************************************************************

  !****************************************************************************
  subroutine Allocate_Working_Constants(STPTable)
  !***  PURPOSE:   to allocate working space
  !
    use MD_TYPEDEF_StopRange_Table
    use MD_Globle_Variables_GPU, only:m_NDEVICE, m_DEVICES
      implicit none
      !--- dummy variables
       type(MDSTPTable),    intent(in)::STPTable
      !--- Local vairables
       integer::I

             !$$--- allocate working space on device 1
              do I=1, m_NDEVICE
                 call Allocate_Working_Constants_template(m_DEVICES(I), STPTable,                                   &
                                                          dm_STCoupingDEVWS(I)%ENABLE, dm_STCoupingDEVWS(I)%KPAIR,  &
                                                          dm_STCoupingDEVWS(I)%MDEN,   dm_STCoupingDEVWS(I)%CM2,    &
                                                          dm_STCoupingDEVWS(I)%LVOL,   dm_STCoupingDEVWS(I)%ETAB,   &
                                                          dm_STCoupingDEVWS(I)%STAB)
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
              do I=1, m_NDEVICE
                 call Clear_Working_Constants_template(m_DEVICES(I), dm_STCoupingDEVWS(I)%ENABLE, dm_STCoupingDEVWS(I)%KPAIR, &
                                                                     dm_STCoupingDEVWS(I)%MDEN,   dm_STCoupingDEVWS(I)%CM2,   &
                                                                     dm_STCoupingDEVWS(I)%LVOL,   dm_STCoupingDEVWS(I)%ETAB,  &
                                                                     dm_STCoupingDEVWS(I)%STAB)
              end do

              hm_NEEDDO    = .false.
              hm_SAVEELOSS = .false.
       return
  end subroutine Clear_Working_Constants
  !**************************************************************************

  !**************************************************************************
  subroutine Allocate_Working_Varibales_template(Idev, Nprt, Eloss)
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
           real(KINDDF), device, dimension(:), allocatable::Eloss

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              allocate(Eloss(Nprt), STAT=ERR)
              ERR = cudaSetDevice(CURDEV)
              if(ERR) then
    100         write(*,*) "MDPSCU Error in allocating working device memory in ST_Coupling module for constant",IDEV
                stop
              end if
       return
  end subroutine Allocate_Working_Varibales_template
  !**************************************************************************

  !**************************************************************************
  subroutine Clear_Working_Varibales_template(IDev, Eloss)
  !***  PURPOSE:   to allocate working space on a device
  !
  !     INPUT:     IDEV     the index of a device
  !
  !     OUTPUT:    allocated memory Eloss
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:),   allocatable::Eloss

      !--- Local vairables
      integer::CURDEV, ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)

              ERR = 0
              if(allocated(Eloss)) deallocate(Eloss, STAT=ERR)
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
      integer::I

             !$$--- allocate working space on devices
              do I=1, m_NDEVICE
                 call Allocate_Working_Varibales_template(m_DEVICES(I), m_NAPDEV, dm_STCoupingDEVWS(I)%ELOSS)
              end do

              allocate(hm_ELOSS(dm_NPRT))
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
      integer::I

             !$$--- deallocate working space on devices
              do I=1, m_NDEVICE
                 call Clear_Working_Varibales_template(m_DEVICES(I), dm_STCoupingDEVWS(I)%ELOSS)
              end do
              if(allocated(hm_ELOSS))deallocate(hm_ELOSS)

              m_INIT = 0
       return
  end subroutine Clear_Working_Varibales
  !****************************************************************************

  !**************************************************************************
  subroutine Clear_STMOD_DEV(SimBox)
  !***  PURPOSE: to clear the memories allocated on calling initialize
  !
  !     INPUT:
  !     OUTPUT:
  !
  implicit none
      !----   DUMMY Variables
      type(SimMDBox), dimension(:), optional::SimBox

         call Clear_Working_Constants()
         call Clear_Working_Varibales()
        return
   end subroutine Clear_STMOD_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Initialize_STMOD_DEV(SimBox, CtrlParam, STPTable)
  !***  PURPOSE: to initialze this module by loading calculation parameters
  !
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  use MD_TYPEDEF_StopRange_Table
  implicit none
      !----   DUMMY Variables
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam
       type(MDSTPTable),intent(in)::STPTable
      !---
         real(KINDDF), dimension(:,:), pointer::ELOSS

         call Clear_STMOD_DEV()

         call Allocate_Working_Constants(STPTable)
         call Reset_STMOD_DEV(SimBox, CtrlParam)

        return
   end subroutine Initialize_STMOD_DEV
  !**************************************************************************

  !**************************************************************************
  subroutine Reset_STMOD_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to reset EPC controal parameters
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
       !--- local
       integer::I, ERR, CURDEV

          hm_ENABLE = iand(CtrlParam%LT_CTRL%METH, CP_TICTRL_METH_ST)
          !$$--- to determine the number density of atoms
          hm_MDEN   = 0.D0
          if(CtrlParam%ST_CTRL%MDEN .gt. 0.000001D0) then
              hm_MDEN  = CtrlParam%ST_CTRL%MDEN*dble(SimBox%NA)/dble(SimBox%NPRT)
              hm_STMOD = mp_STMOD_G
          else
            !$$--- if MDEN = 0, the density is calculated by number of atoms divided by box size
             if(dabs(CtrlParam%ST_CTRL%MDEN) .le. 0.00001D0) then
                if(all(CtrlParam%IFPD .gt. 0) ) then
                   hm_MDEN  = dble(SimBox%NPRT)/(Simbox%ZL(1)*Simbox%ZL(2)*Simbox%ZL(3))
                   hm_MDEN  = hm_MDEN*dble(SimBox%NA)/dble(SimBox%NPRT)
                   hm_STMOD = mp_STMOD_G
                else
                   if(hm_STMOD .ne. mp_STMOD_L) then
                      write(*,fmt="(' MDPSCU Warning: all periodic boundaries condition is not applied to the box', BZI6)")
                      write(*,fmt="('                 the atom number density in the box may be incorrect for STOPPING force calculation ')")
                      call ONWARNING(gm_OnWarning)
                   end if
                end if

             else
                hm_STMOD = mp_STMOD_L
             end if
          end if

           !$$--- copyin the control parameters
          ERR = cudaGetDevice(CURDEV)
          do I=1, m_NDEVICE
             ERR = cudaSetDevice(m_DEVICES(I))
             dm_STCoupingDEVWS(I)%ENABLE = hm_ENABLE
             dm_STCoupingDEVWS(I)%MDEN   = hm_MDEN
             dm_STCoupingDEVWS(I)%CM2    = 0.5D0*SimBox%CM
             dm_STCoupingDEVWS(I)%LVOL   = CP_4PI3*(CtrlParam%NB_RM**3.D0)
          end do
          ERR = cudaSetDevice(CURDEV)

          call Clear_Working_Varibales()
          if(all(hm_ENABLE .eq. 0)) then
             hm_NEEDDO = .false.
          else
             hm_NEEDDO  = .true.
             if(CtrlParam%ST_CTRL%SaveEloss .gt. 0) then
                hm_SAVEELOSS = .true.
             else
                hm_SAVEELOSS = .false.
             endif

          end if

        return
   end subroutine Reset_STMOD_DEV
  !**************************************************************************

  !**************************************************************************
  attributes(global) subroutine ST_MOD_GDEN_KERNEL(IM, IA0, NAPDEV, NPART, NG,   &
                                   ENABLE, MDEN, CM2, KPAIR, NE, NK, ETAB, STAB, &
                                   ITYP, STATU, XP1, FP)
  !***  PURPOSE:   KERNEL to modify the force by add the stopping power to forces.
  !$$                     The density of atoms are previously calculated by divding
  !$$                     the number of in-box atoms/box size
  !$$
  !$$   INPUT:     IM,        the total number of atoms in the whole box
  !$$              IA0,       the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$
  !$$              NG,        the number of atom types
  !$$              ENABLE,    variable indicating if stopping calculation is enabled for types of atom types
  !$$              MDEN,      the density of each type of atoms
  !$$              CM2,       0.5*mass of atoms
  !$$              KPAIR
  !$$              NE,        the number of energy points in the stopping table
  !$$              NK,        the number of E-ST tables
  !$$              ETAB,      the energy grids
  !$$              STAB,      the ST tables
  !$$              ITYP,      indicator of type of atoms
  !$$              STATU,     statu of atoms
  !$$              XP1,       velocity of atoms
  !$$              FP,        force on atoms
  !
  !$$   OUTPUT     FP,        updated force
  !
  !$$              NOTE:  the stopping to be added to force calculated by force
  !$$                     calculation. Thus the present routine should be called
  !$$                     force routine are called.
  implicit none
  !----   DUMMY Variables
          integer, value::IM, IA0, NAPDEV, NPART, NG, NE, NK

          integer,     device::ENABLE(NG)
          real(KINDDF),device::MDEN(NG)
          real(KINDDF),device::CM2(NG)
          integer,     device::KPAIR(NG,*)
          real(KINDDF),device::ETAB(NE)
          real(KINDDF),device::STAB(NE,*)
          integer,     device::ITYP(IM)
          integer,     device::STATU(NAPDEV)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)

  !----   Local variables
          integer::KK, KK1, IT, IB, IC, IG, IG1, IG2, IK, IK1, KP
          real(KINDDF)::EK, VX, VY, VZ, VV, FF, EMIN, EMAX
          integer,      shared::KS(mp_mxGROUP*mp_mxGROUP)
          real(KINDDF), shared::SK(mp_mxGROUP), DEINV(1)
          !integer,      shared::KS(NG*NG)
          !real(KINDDF), shared::SDEN(NG)
  !---
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IT .eq. 1) then
                 DEINV(1) = 1.D0/(ETAB(2) - ETAB(1))
              end if

              !$$--- the linear interpolation coefficients
              if(IT .le. NG) then
                 SK(IT) =  MDEN(IT)/(ETAB(2) - ETAB(1))
              end if

              if(IT .le. NG*NG) then
                 IG1 = (IT-1)/NG
                 IG2 = IT - IG1*NG
                 KS(IT) =  KPAIR(IG1+1, IG2)
              end if
              call syncthreads()

              EMIN  = ETAB(1)
              EMAX  = ETAB(NE)
              if(IC.LE.NPART) then
                 KK = ITYP(IC+IA0)
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE .and. &
                    ENABLE(KK) .gt. 0 ) then
                     VX  = XP1(IC,1)
                     VY  = XP1(IC,2)
                     VZ  = XP1(IC,3)
                     VV  = VX*VX+VY*VY+VZ*VZ
                     EK  = CM2(KK)*VV

                     if(EK .ge. EMIN .and. EK .le. EMAX) then
                        IK  = int((EK-EMIN)*DEINV(1)) + 1
                        IK1 = IK + 1
                        FF  = 0.D0
                        KK1 = (KK-1)*NG
                        do IG = 1, NG
                            KP  = KS(KK1+IG)
                            FF  = FF + SK(IG)* &
                                 ( ((EK-ETAB(IK))*STAB(IK1,KP) + (ETAB(IK1)-EK)*STAB(IK,KP)) )
                        end do
                        VV  = dsqrt(VV)

                        FP(IC,1) = FP(IC,1) - FF*VX/VV
                        FP(IC,2) = FP(IC,2) - FF*VY/VV
                        FP(IC,3) = FP(IC,3) - FF*VZ/VV
                     end if
                 end if
             end if
        return
  end subroutine ST_MOD_GDEN_KERNEL
  !****************************************************************************************

  !****************************************************************************************
  subroutine ST_MOD_GDEN_template(IDEV, STARTCELL, ENDCELL, ENABLE, MDEN, CM2, KPAIR, &
                                  ETAB, STAB, ITYP, STATU, XP1, FP)
  !***  PORPOSE:  template for devices invoking the KERNEL
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                TE,        assumed electron temperature
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
          integer,      device, dimension(:)  ::ENABLE
          real(KINDDF), device, dimension(:)  ::MDEN
          real(KINDDF), device, dimension(:)  ::CM2
          integer,      device, dimension(:,:)::KPAIR
          real(KINDDF), device, dimension(:)  ::ETAB
          real(KINDDF), device, dimension(:,:)::STAB
          integer,      device, dimension(:)  ::ITYP
          integer,      device, dimension(:)  ::STATU
          real(KINDDF), device, dimension(:,:)::XP1,FP

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
            BX = mp_BLOCKSIZE
            BY = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
           !$$-- NOTE: the upper limit of grid size is 65535
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call ST_MOD_GDEN_KERNEL<<<blocks, threads>>>(dm_NPRT,STARTA, m_NAPDEV, NPRT, dm_NGROUP,   &
                                                 ENABLE, MDEN, CM2, KPAIR, m_NE, m_NK, ETAB, STAB,    &
                                                 ITYP, STATU, XP1, FP)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine ST_MOD_GDEN_template
  !****************************************************************************************

  !****************************************************************************************
  attributes(global) subroutine ST_MOD_LDEN_KERNEL(IM, IA0, NAPDEV, NPART, KVOIS, INDI,&
                                   NG, ENABLE, LV, CM2, KPAIR, NE, NK, ETAB, STAB,    &
                                   ITYP, STATU, XP1, FP)
  !***  PURPOSE:   KERNEL to modify the force by add the stopping power to forces.
  !$$                     The densities of atoms are calculated locally:
  !$$                     the number of neighboring atoms/sphereical neighboring region.
  !$$                     The neighbor-list of atoms is required for stopping calculation.
  !$$
  !$$   INPUT:     IM,        the total number of atoms in the whole box
  !$$              IA0,       the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$              KVOIS:     the number of neighbors for atoms
  !$$              INDI:      the index for the neighbores
  !$$
  !$$              NG,        the number of atom types
  !$$              ENABLE,    variable indicating if stopping calculation is enabled for types of atom types
  !$$              LV,        the local neighboring volume
  !$$              CM2,       0.5*mass of atoms
  !$$              KPAIR
  !$$              NE,        the number of energy points in the stopping table
  !$$              NK,        the number of E-ST tables
  !$$              ETAB,      the energy grids
  !$$              STAB,      the ST tables
  !$$              ITYP,      indicator of type of atoms
  !$$              STATU,     statu of atoms
  !$$              XP1,       velocity of atoms
  !$$              FP,        force on atoms
  !
  !$$   OUTPUT     FP,        updated force
  !
  !$$              NOTE:  the stopping to be added to force calculated by force
  !$$                     calculation. Thus the present routine should be called
  !$$                     force routine are called.
  implicit none
  !----   DUMMY Variables
          integer, value::IM, IA0, NAPDEV, NPART, NG, NE, NK

          integer,     device::ENABLE(NG)
          real(KINDDF),device::LV(mp_mxGROUP,*)
          real(KINDDF),device::CM2(NG)
          integer,     device::KPAIR(NG,*)
          real(KINDDF),device::ETAB(NE)
          real(KINDDF),device::STAB(NE,*)

          integer,     device, intent(in)::KVOIS(NAPDEV)
          integer,     device, intent(in)::INDI(NAPDEV,*)
          integer,     device::ITYP(IM)
          integer,     device::STATU(NAPDEV)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)

  !----   Local variables
          integer     ::ITYPI, ITYPJ, KK1, IT, IB, IC, IG, IG1, IG2, IK, IK1, KP, IW, IIW
          real(KINDDF)::EK, VX, VY, VZ, VV, FF, EMIN, EMAX, DEN(mp_mxGROUP)
          integer,      shared::KS(mp_mxGROUP*mp_mxGROUP)
          real(KINDDF), shared::SK(mp_mxGROUP),DEINV(1),ILV(mp_mxGROUP*mp_mxGROUP)
  !---
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IT .eq. 1) then
                 DEINV(1) = 1.D0/(ETAB(2) - ETAB(1))
              end if

              if(IT .le. NG*NG) then
                 IG1     = (IT-1)/NG
                 IG2     = IT - IG1*NG
                 KS(IT)  =  KPAIR(IG1+1, IG2)
                 ILV(IT) =  1.D0/LV(IG1+1, IG2)/(ETAB(2) - ETAB(1))
              end if
              call syncthreads()

              EMIN  = ETAB(1)
              EMAX  = ETAB(NE)
              if(IC.LE.NPART) then
                 ITYPI = ITYP(IC+IA0)
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE .and. &
                    ENABLE(ITYPI) .gt. 0 ) then
                     VX  = XP1(IC,1)
                     VY  = XP1(IC,2)
                     VZ  = XP1(IC,3)
                     VV  = VX*VX+VY*VY+VZ*VZ
                     EK  = CM2(ITYPI)*VV
                     IIW = KVOIS(IC)

                     if(EK .ge. EMIN .and. EK .le. EMAX .and. IIW.gt.0) then
                        !$$--- first, we need to calculate the local density of atoms
                        DEN(1:NG) = 0.D0
                        do IW=1, IIW
                            IG      = ITYP(INDI(IC,IW))
                            DEN(IG) = DEN(IG) + 1.D0
                        end do

                        IK  = int((EK-EMIN)*DEINV(1)) + 1
                        IK1 = IK + 1
                        FF  = 0.D0
                        KK1 = (ITYPI-1)*NG
                        do IG = 1, NG
                            KP  = KS(KK1+IG)
                            FF  = FF + DEN(IG)*ILV(KK1+IG)* &
                                 ( ((EK-ETAB(IK))*STAB(IK1,KP) + (ETAB(IK1)-EK)*STAB(IK,KP)) )
                        end do
                        VV  = dsqrt(VV)

                        FP(IC,1) = FP(IC,1) - FF*VX/VV
                        FP(IC,2) = FP(IC,2) - FF*VY/VV
                        FP(IC,3) = FP(IC,3) - FF*VZ/VV
                     end if
                 end if
             end if
        return
  end subroutine ST_MOD_LDEN_KERNEL
  !****************************************************************************************

  !****************************************************************************************
  subroutine ST_MOD_LDEN_template(IDEV, STARTCELL, ENDCELL, ENABLE, KVOIS, INDI, LV, CM2, KPAIR, &
                                  ETAB, STAB, ITYP, STATU, XP1, FP)
  !***  PORPOSE:  template for devices invoking the KERNEL
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                TE,        assumed electron temperature
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
          integer,      device, dimension(:)  ::ENABLE
          integer,      device, dimension(:)  ::KVOIS
          integer,      device, dimension(:,:)::INDI
          real(KINDDF), device, dimension(:,:)::LV
          real(KINDDF), device, dimension(:)  ::CM2
          integer,      device, dimension(:,:)::KPAIR
          real(KINDDF), device, dimension(:)  ::ETAB
          real(KINDDF), device, dimension(:,:)::STAB
          integer,      device, dimension(:)  ::ITYP
          integer,      device, dimension(:)  ::STATU
          real(KINDDF), device, dimension(:,:)::XP1,FP

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
            BX = mp_BLOCKSIZE
            BY = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
           !$$-- NOTE: the upper limit of grid size is 65535
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1 
            call ST_MOD_LDEN_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, m_NAPDEV, NPRT, KVOIS, INDI,  &
                                                 dm_NGROUP, ENABLE, LV, CM2, KPAIR, m_NE, m_NK,         &
                                                 ETAB, STAB, ITYP, STATU, XP1, FP)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine ST_MOD_LDEN_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_STMOD_Force_DEV(SimBox, CtrlParam)
  !***  PORPOSE: to calculate the friction force due to ELOSS and add them to the total force
  !              of atoms
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  use MD_Globle_Variables_GPU, only: dm_WorkSpace, m_STARTCELL, m_ENDCELL
  use MD_NeighborsList_GPU,    only: dm_Neighbors
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)               ::CtrlParam
  !--- Device variables and variables to be used in GPU
          integer::I

          select case(hm_STMOD)
                 case(mp_STMOD_G)
                     do I=1, m_NDEVICE
                        call ST_MOD_GDEN_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),       &
                                     dm_STCoupingDEVWS(I)%ENABLE, dm_STCoupingDEVWS(I)%MDEN,        &
                                     dm_STCoupingDEVWS(I)%CM2,    dm_STCoupingDEVWS(I)%KPAIR,       &
                                     dm_STCoupingDEVWS(I)%ETAB,   dm_STCoupingDEVWS(I)%STAB,        &
                                     dm_WorkSpace%ITYP(I)%Data,   dm_WorkSpace%STATU(I)%Data,       &
                                     dm_WorkSpace%XP1(I)%Data,    dm_WorkSpace%FP(I)%Data           )
                     end do

                 case(mp_STMOD_L)
                     do I=1, m_NDEVICE
                        call ST_MOD_LDEN_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),       &
                                     dm_STCoupingDEVWS(I)%ENABLE,                                   &
                                     dm_Neighbors%KVOIS(I)%Data,  dm_Neighbors%INDI(I)%Data,        &
                                     dm_STCoupingDEVWS(I)%LVOL,   dm_STCoupingDEVWS(I)%CM2,         &
                                     dm_STCoupingDEVWS(I)%KPAIR,  dm_STCoupingDEVWS(I)%ETAB,        &
                                     dm_STCoupingDEVWS(I)%STAB,                                     &
                                     dm_WorkSpace%ITYP(I)%Data,   dm_WorkSpace%STATU(I)%Data,       &
                                     dm_WorkSpace%XP1(I)%Data,    dm_WorkSpace%FP(I)%Data           )

                     end do
          end select

     return
  end subroutine Do_STMOD_Force_DEV
  !****************************************************************************************

  !**************************************************************************
  attributes(global) subroutine ST_MOD_ELOSS_GDEN_KERNEL(IM, IA0, NAPDEV, NPART, NG,   &
                                   ENABLE, MDEN, CM2, KPAIR, NE, NK, ETAB, STAB,       &
                                   ITYP, STATU, XP1, FP, DT, ELOSS)
  !$$--   PURPOSE: KERNEL  to modify the force by add the stopping power to forces,
  !$$                       and with the energy loss in time interval DT saved
  !$$
  !$$                       The density of atoms are previously calculated by divding
  !$$                       the number of in-box atoms/box size
  !$$
  !$$   INPUT:     IM,        the total number of atoms in the whole box
  !$$              IA0,       the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$
  !$$              NG,        the number of atom types
  !$$              ENABLE,    variable indicating if stopping calculation is enabled for types of atom types
  !$$              MDEN,      the density of each type of atoms
  !$$              CM2,       0.5*mass of atoms
  !$$              KPAIR
  !$$              NE,        the number of energy points in the stopping table
  !$$              NK,        the number of E-ST tables
  !$$              ETAB,      the energy grids
  !$$              STAB,      the ST tables
  !$$              ITYP,      indicator of type of atoms
  !$$              STATU,     statu of atoms
  !$$              XP1,       velocity of atoms
  !$$              DT,        the time interval
  !$$              FP,        force on atoms
  !
  !$$   OUTPUT     FP,        updated force
  !$$              ELOSS,     the energy loss
  !$$              NOTE:      the stopping to be added to force calculated by force
  !$$                         calculation. Thus the present routine should be called
  !$$                         after force routine are called.
  implicit none
  !----   DUMMY Variables
          integer,     value::IM, IA0, NAPDEV, NPART, NG, NE, NK
          real(KINDDF),value::DT

          integer,     device::ENABLE(NG)
          real(KINDDF),device::MDEN(NG)
          real(KINDDF),device::CM2(NG)
          integer,     device::KPAIR(NG,*)
          real(KINDDF),device::ETAB(NE)
          real(KINDDF),device::STAB(NE,*)
          integer,     device::ITYP(IM)
          integer,     device::STATU(NAPDEV)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)
          real(KINDDF),device::ELOSS(*)

  !----   Local variables
          integer::KK, KK1, IT, IB, IC, IG, IG1, IG2, IK, IK1, KP
          real(KINDDF)::EK, VX, VY, VZ, VV, FF, EMIN, EMAX
          integer,      shared::KS(mp_mxGROUP*mp_mxGROUP)
          real(KINDDF), shared::SK(mp_mxGROUP), DEINV(1)
  !---
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IT .eq. 1) then
                 DEINV(1) = 1.D0/(ETAB(2) - ETAB(1))
              end if

              !$$--- the linear interpolation coefficients
              if(IT .le. NG) then
                 SK(IT) =  MDEN(IT)/(ETAB(2) - ETAB(1))
              end if

              if(IT .le. NG*NG) then
                 IG1 = (IT-1)/NG
                 IG2 = IT - IG1*NG
                 KS(IT) =  KPAIR(IG1+1, IG2)
              end if
              call syncthreads()

              EMIN  = ETAB(1)
              EMAX  = ETAB(NE)
              if(IC.LE.NPART) then
                 KK = ITYP(IC+IA0)
                 ELOSS(IC) = 0.D0

                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE .and. &
                    ENABLE(KK) .gt. 0 ) then
                     VX  = XP1(IC,1)
                     VY  = XP1(IC,2)
                     VZ  = XP1(IC,3)
                     VV  = VX*VX+VY*VY+VZ*VZ
                     EK  = CM2(KK)*VV
                     if(EK .ge. EMIN .and. EK .le. EMAX) then
                        IK  = int((EK-EMIN)*DEINV(1)) + 1
                        IK1 = IK + 1
                        KK1 = (KK-1)*NG
                        FF  = 0.D0
                        do IG = 1, NG
                            KP  = KS(KK1+IG)
                            FF  = FF + SK(IG)* &
                                 ( ((EK-ETAB(IK))*STAB(IK1,KP) + (ETAB(IK1)-EK)*STAB(IK,KP)) )
                        end do
                        VV  = dsqrt(VV)

                        FP(IC,1)  = FP(IC,1) - FF*VX/VV
                        FP(IC,2)  = FP(IC,2) - FF*VY/VV
                        FP(IC,3)  = FP(IC,3) - FF*VZ/VV
                        ELOSS(IC) = FF*VV*DT
                     end if
                 end if
             end if
        return
  end subroutine ST_MOD_ELOSS_GDEN_KERNEL
  !****************************************************************************************

  !****************************************************************************************
  subroutine ST_MOD_ELOSS_GDEN_template(IDEV, STARTCELL, ENDCELL, ENABLE, MDEN, CM2, KPAIR, &
                                  ETAB, STAB, ITYP, STATU, XP1, FP, DT, ELOSS)
  !***  PORPOSE:   template for devices invoking the KERNEL
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                TE,        assumed electron temperature
  !                ITYP,      type of atoms
  !                XP1,       velocity of atoms
  !                FP,        force on the atoms
  !
  !     OUTPUT     FP,        updated force
  !                ELOSS,     electronic energy loss
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)::IDEV, STARTCELL, ENDCELL
          integer,      device, dimension(:)  ::ENABLE
          real(KINDDF), device, dimension(:)  ::MDEN
          real(KINDDF), device, dimension(:)  ::CM2
          integer,      device, dimension(:,:)::KPAIR
          real(KINDDF), device, dimension(:)  ::ETAB
          real(KINDDF), device, dimension(:,:)::STAB
          integer,      device, dimension(:)  ::ITYP
          integer,      device, dimension(:)  ::STATU
          real(KINDDF), device, dimension(:,:)::XP1,FP
          real(KINDDF), intent(in)            ::DT
          real(KINDDF), device, dimension(:)  ::ELOSS

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
            BX = mp_BLOCKSIZE
            BY = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
           !$$-- NOTE: the upper limit of grid size is 65535
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call ST_MOD_ELOSS_GDEN_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, m_NAPDEV, NPRT, dm_NGROUP, &
                                                 ENABLE, MDEN, CM2, KPAIR, m_NE, m_NK, ETAB, STAB,         &
                                                 ITYP, STATU, XP1, FP, DT, ELOSS)

            ERR = cudaSetDevice(CURDEV)
      return
  end subroutine ST_MOD_ELOSS_GDEN_template
  !****************************************************************************************

  !****************************************************************************************
  attributes(global) subroutine ST_MOD_ELOSS_LDEN_KERNEL(IM, IA0, NAPDEV, NPART, KVOIS, INDI,&
                                   NG, ENABLE, LV, CM2, KPAIR, NE, NK, ETAB, STAB,    &
                                   ITYP, STATU, XP1, FP, DT, ELOSS)
  !$$***PURPOSE:   KERNEL to modify the force by add the stopping power to forces.
  !$$                     and with the energy loss in time interval DT saved
  !$$                     the number of neighboring atoms/sphereical neighboring region.
  !$$                     The neighbor-list of atoms is required for stopping calculation.
  !$$
  !$$   INPUT:     IM,        the total number of atoms in the whole box
  !$$              IA0,       the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$              KVOIS:     the number of neighbors for atoms
  !$$              INDI:      the index for the neighbores
  !$$
  !$$              NG,        the number of atom types
  !$$              ENABLE,    variable indicating if stopping calculation is enabled for types of atom types
  !$$              LV,        the local neighboring volume
  !$$              CM2,       0.5*mass of atoms
  !$$              KPAIR
  !$$              NE,        the number of energy points in the stopping table
  !$$              NK,        the number of E-ST tables
  !$$              ETAB,      the energy grids
  !$$              STAB,      the ST tables
  !$$              ITYP,      indicator of type of atoms
  !$$              STATU,     statu of atoms
  !$$              XP1,       velocity of atoms
  !$$              FP,        force on atoms
  !
  !$$   OUTPUT     FP,        updated force
  !$$              ELOSS,     the energy loss
  !
  !$$              NOTE:  the stopping to be added to force calculated by force
  !$$                     calculation. Thus the present routine should be called
  !$$                     force routine are called.
  implicit none
  !----   DUMMY Variables
          integer,     value::IM, IA0, NAPDEV, NPART, NG, NE, NK
          real(KINDDF),value::DT
          integer,     device::ENABLE(NG)
          real(KINDDF),device::LV(mp_mxGROUP,*)
          real(KINDDF),device::CM2(NG)
          integer,     device::KPAIR(NG,*)
          real(KINDDF),device::ETAB(NE)
          real(KINDDF),device::STAB(NE,*)

          integer,     device, intent(in)::KVOIS(NAPDEV)
          integer,     device, intent(in)::INDI(NAPDEV,*)
          integer,     device::ITYP(IM)
          integer,     device::STATU(NAPDEV)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)
          real(KINDDF),device::ELOSS(*)

  !----   Local variables
          integer     ::ITYPI, ITYPJ, KK1, IT, IB, IC, IG, IG1, IG2, IK, IK1, KP, IW, IIW
          real(KINDDF)::EK, VX, VY, VZ, VV, FF, EMIN, EMAX, DEN(mp_mxGROUP)
          integer,      shared::KS(mp_mxGROUP*mp_mxGROUP)
          real(KINDDF), shared::SK(mp_mxGROUP),DEINV(1),ILV(mp_mxGROUP*mp_mxGROUP)
  !---
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IT .eq. 1) then
                 DEINV(1) = 1.D0/(ETAB(2) - ETAB(1))
              end if

              if(IT .le. NG*NG) then
                 IG1     = (IT-1)/NG
                 IG2     = IT - IG1*NG
                 KS(IT)  =  KPAIR(IG1+1, IG2)
                 ILV(IT) =  1.D0/LV(IG1+1, IG2)/(ETAB(2) - ETAB(1))
              end if
              call syncthreads()

              EMIN  = ETAB(1)
              EMAX  = ETAB(NE)
              if(IC.LE.NPART) then
                 ITYPI = ITYP(IC+IA0)
                 ELOSS(IC)= 0.D0
                 if(iand(STATU(IC), CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE .and. &
                    ENABLE(ITYPI) .gt. 0 ) then
                     VX  = XP1(IC,1)
                     VY  = XP1(IC,2)
                     VZ  = XP1(IC,3)
                     VV  = VX*VX+VY*VY+VZ*VZ
                     EK  = CM2(ITYPI)*VV
                     IIW = KVOIS(IC)

                     if(EK .ge. EMIN .and. EK .le. EMAX .and. IIW.gt.0) then
                        !$$--- first, we need to calculate the local density of atoms
                        DEN(1:NG) = 0.D0
                        do IW=1, IIW
                            IG      = ITYP(INDI(IC,IW))
                            DEN(IG) = DEN(IG) + 1.D0
                        end do

                        IK  = int((EK-EMIN)*DEINV(1)) + 1
                        IK1 = IK + 1
                        FF  = 0.D0
                        KK1 = (ITYPI-1)*NG
                        do IG = 1, NG
                            KP  = KS(KK1+IG)
                            FF  = FF + DEN(IG)*ILV(KK1+IG)* &
                                 ( ((EK-ETAB(IK))*STAB(IK1,KP) + (ETAB(IK1)-EK)*STAB(IK,KP)) )
                        end do
                        VV  = dsqrt(VV)

                        FP(IC,1) = FP(IC,1) - FF*VX/VV
                        FP(IC,2) = FP(IC,2) - FF*VY/VV
                        FP(IC,3) = FP(IC,3) - FF*VZ/VV
                        ELOSS(IC)= FF*VV*DT
                     end if
                 end if
             end if
        return
  end subroutine ST_MOD_ELOSS_LDEN_KERNEL
  !****************************************************************************************

  !****************************************************************************************
  subroutine ST_MOD_ELOSS_LDEN_template(IDEV, STARTCELL, ENDCELL, ENABLE, KVOIS, INDI, LV, CM2, KPAIR, &
                                  ETAB, STAB, ITYP, STATU, XP1, FP, DT, ELOSS)
  !***  PORPOSE:  template for devices invoking the KERNEL
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                TE,        assumed electron temperature
  !                ITYP,      type of atoms
  !                XP1,       velocity of atoms
  !                FP,        force on the atoms
  !                DT,        the size of time step
  !
  !     OUTPUT     FP,        updated force
  !                ELOSS,     electronic energy loss
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)::IDEV, STARTCELL, ENDCELL
          integer,      device, dimension(:)  ::ENABLE
          integer,      device, dimension(:)  ::KVOIS
          integer,      device, dimension(:,:)::INDI
          real(KINDDF), device, dimension(:,:)::LV
          real(KINDDF), device, dimension(:)  ::CM2
          integer,      device, dimension(:,:)::KPAIR
          real(KINDDF), device, dimension(:)  ::ETAB
          real(KINDDF), device, dimension(:,:)::STAB
          integer,      device, dimension(:)  ::ITYP
          integer,      device, dimension(:)  ::STATU
          real(KINDDF), device, dimension(:,:)::XP1,FP
          real(KINDDF), intent(in)            ::DT
          real(KINDDF), device, dimension(:)  ::ELOSS

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
            BX = mp_BLOCKSIZE
            BY = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
           !$$-- NOTE: the upper limit of grid size is 65535
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1 
            call ST_MOD_ELOSS_LDEN_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, m_NAPDEV, NPRT, KVOIS, INDI, &
                                                 dm_NGROUP, ENABLE, LV, CM2, KPAIR, m_NE, m_NK,              &
                                                 ETAB, STAB, ITYP, STATU, XP1, FP, DT, ELOSS)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine ST_MOD_ELOSS_LDEN_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_STMOD_Force_Eloss_DEV(SimBox, CtrlParam)
  !***  PORPOSE: to calculate the force caused by stopping and the also
  !              the energy loss
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  !     OUTPUT     dxm_FP,    forces of atoms on all devices
  !
  use MD_Globle_Variables_GPU
  use MD_NeighborsList_GPU,   only:dm_Neighbors
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
  !--- Device variables and variables to be used in GPU
          real(KINDDF)::DT
          integer::I 

          if(m_Init .le. 0) then
             call Allocate_Working_Varibales()
          end if
          !--- NOTE: the time step CtrlParam%H could be varing
          !          ref: MD_DiffScheme_GPU.F90
          DT = CtrlParam%H
          select case(hm_STMOD)
                 !$$--- in case the atom density is calculated by NPRT/BOXSIZE
                 case(mp_STMOD_G)
                     do I=1, m_NDEVICE
                        call ST_MOD_ELOSS_GDEN_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), &
                                     dm_STCoupingDEVWS(I)%ENABLE, dm_STCoupingDEVWS(I)%MDEN,        &
                                     dm_STCoupingDEVWS(I)%CM2,    dm_STCoupingDEVWS(I)%KPAIR,       &
                                     dm_STCoupingDEVWS(I)%ETAB,   dm_STCoupingDEVWS(I)%STAB,        &
                                     dm_WorkSpace%ITYP(I)%Data,   dm_WorkSpace%STATU(I)%Data,       &
                                     dm_WorkSpace%XP1(I)%Data,    dm_WorkSpace%FP(I)%Data,          &
                                     DT, dm_STCoupingDEVWS(I)%ELOSS                                )
                     end do

                 !$$--- in case the local atom density is used
                 case(mp_STMOD_L)
                     do I=1, m_NDEVICE
                        call ST_MOD_ELOSS_LDEN_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), &
                                     dm_STCoupingDEVWS(I)%ENABLE,                                   &
                                     dm_Neighbors%KVOIS(I)%Data,  dm_Neighbors%INDI(I)%Data,        &
                                     dm_STCoupingDEVWS(I)%LVOL,   dm_STCoupingDEVWS(I)%CM2,         &
                                     dm_STCoupingDEVWS(I)%KPAIR,  dm_STCoupingDEVWS(I)%ETAB,        &
                                     dm_STCoupingDEVWS(I)%STAB,                                     &
                                     dm_WorkSpace%ITYP(I)%Data,   dm_WorkSpace%STATU(I)%Data,       &
                                     dm_WorkSpace%XP1(I)%Data,    dm_WorkSpace%FP(I)%Data,          &
                                     DT, dm_STCoupingDEVWS(I)%ELOSS                                 )

                     end do
          end select

     return
  end subroutine Do_STMOD_Force_Eloss_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_STMOD_DEV(SimBox, CtrlParam)
  !***  PORPOSE: to interface to host process
  !
  !     INPUT:     SimBox,    simulation box, not used
  !                CtrlParam, control parameters, not used
  !
  use MD_Globle_Variables_GPU, only:hm_GIDINV
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
  !---  local
          integer::I, IA1, IA2

          if(.not.hm_NEEDDO) return

          if(hm_SAVEELOSS) then
             call Do_STMOD_Force_Eloss_DEV(SimBox, CtrlParam)
             call CopyElossFrom_Devices_to_Host()
             IA1 = 1
             IA2 = SimBox(1)%NPRT
             do I=1, size(SimBox)
                !$$--- NOTE: we use the whole array hm_ELOSS, but not its subset(IA1:IA2)
                call AccumDatPad_SimMDBox(SimBox(I), "Eloss (ev)", hm_ELOSS*CP_ERGEV, &
                                               hm_GIDINV(IA1:IA2))
                IA1 = IA2 + 1
                IA2 = IA2 + SimBox(I)%NPRT
             end do
          else
             call Do_STMOD_Force_DEV(SimBox, CtrlParam)
          end if
          return
  end subroutine Do_STMOD_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine CopyElossFrom_Devices_to_Host()
  !***  PURPOSE:  copy Eloss on devices to the arries on host.
  !
    use MD_Globle_Variables_GPU, only:m_NDEVICE, COPY_OUT_SHIFT_template
    implicit none
      !--- dummy variables
      !----   Local variables
       integer::I

            do I=1, m_NDEVICE
               call COPY_OUT_SHIFT_template(I, dm_STCoupingDEVWS(I)%ELOSS, hm_ELOSS)
            end do 

            return
   end subroutine CopyElossFrom_Devices_to_Host
  !****************************************************************************
  end module MD_ST_Coupling_GPU
