  module MD_DiffScheme_GPU
  !**** DESCRIPTION: Use the differetial scheme proposed by Swope et al,J.Chem.Phys.76(1982)637-649
  !                  Berendson pressure coupling is realized in this module, HOWEVER,
  !                  the coordinate system is assumed to be Cartesian.
  !
  !                  To consider non-cubic box, Gear Scheme may be needed
  !
  !                  This version was updated from MD_SWOPE_Scheme_GPU
  !                  with mult-device suported
  !                  ______________________________________________________
  !**** HISTORY:     2012-03 (HOU Qing), created the first version
  !
  !                  2018-05-31(HOU Qing),
  !                          remove the FIXTEMP operation in Correction subroutine.
  !                          Insteadly, the local temperature for atoms is handled
  !                          in  MD_LocalTempMethod module.
  !
  !                  2018-12-31,  remove the constant memory dcm_LBOX, dcm_HBOX and dcm_BOXSIZE, and thus
  !                               thus the copyin operations of the constant memory in PRECALFOR_template .
  !                               The boxsize and shape are passed to kernel by BX0,BY0,....
  !                               Based on DDT-MAP, This modification maybe profit for acceleratting the 
  !                               force calculations.
  !     

  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_MultiGPU_Basic
  implicit none
  !**********************
         integer, parameter,private::mp_BLOCKSIZE=256
         integer, parameter,private::mp_BLOCKDIMX=64

         integer,      private::m_NGROUP
         real(KINDDF), private::m_TH                             !time step size (in second)
         real(KINDDF), private::m_H2S2
         real(KINDDF), private::m_HS2
         real(KINDDF), private::m_MXD2

         real(KINDDF), private::m_BOXSIZE(3), m_LBOX(3), m_UBOX(3)
         integer,      private::m_IFPD(3)
         integer,      private, constant::dcm_IFPD(3)
         real(KINDDF), private::m_CM(mp_mxGroup)

         private:: Diff_DEVWS
         type:: Diff_DEVWS
             logical,      device, dimension(:), allocatable::IFlag
         end type Diff_DEVWS
         type(Diff_DEVWS), dimension(m_MXDEVICE), private::dm_DiffWS

        private::COPYOUTEKIN_template
        private::CALEKIN_KERNEL
        private::CALEKIN_template
        private::CheckTimestep_KERNEL
        private::CheckTimestep_template
        private::Correction_KERNEL
        private::Correction_template
        private::DAMPing_KERNEL
        private::DAMPing_template
        private::Predictor_KERNEL0
        private::Predictor_KERNEL1
        private::Predictor_template
        private::Thermalizing_MC_KERNEL
        private::Thermalizing_MC_template
        private::VelScaling_KERNEL
        private::VelScaling_template
  contains

  !****************************************************************************************
  subroutine Initialize_Diff_Scheme_DEV(SimBox,CtrlParam)
  !***  PORPOSE:   to initialize this module by allocating neccessary memories
  !                and initialize parameters
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT:
  !
  implicit none
  !----   DUMMY Variables
          type(SimMDBox) ::SimBox
          type(SimMDCtrl)::CtrlParam
  !----  Local variables
          integer::ERR, CURDEV, NG, I

               m_NGROUP = SimBox%NGROUP
               m_CM(1:m_NGROUP) = SimBox%CM(1:m_NGROUP)
               m_IFPD  = CtrlParam%IFPD

               ERR = cudaGetDevice(CURDEV)
               do I=1, m_NDEVICE
                  ERR = cudaSetDevice(m_DEVICES(I))
                  if(allocated(dm_DiffWS(I)%IFlag)) deallocate(dm_DiffWS(I)%IFlag)
                  allocate(dm_DiffWS(I)%IFlag(1))
                  dcm_IFPD(1:3) = m_IFPD(1:3)
               end do
               ERR = cudaSetDevice(CURDEV)

     return
  end subroutine Initialize_Diff_Scheme_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine Clear_Diff_Scheme_DEV()
  !***  PORPOSE:   to release the resource allocated by the module
  !     INPUT:
  !     OUTPUT
  !
  implicit none
  !----   DUMMY Variables
  !----  Local variables
      integer::ERR, CURDEV, I

               ERR = cudaGetDevice(CURDEV)

               do I=1, m_NDEVICE 
                  ERR = cudaSetDevice(m_DEVICES(I))
                  if(allocated(dm_DiffWS(I)%IFlag)) deallocate(dm_DiffWS(I)%IFlag)
               end do
               ERR = cudaSetDevice(CURDEV)
               return
  end subroutine Clear_Diff_Scheme_DEV
  !*****************************************************************************************

  !*****************************************************************************************
  attributes(global) subroutine DAMPING_KERNEL(NAPDEV, NPART,XP1,FP,STATU)
  !***  PURPOSE:   to damp the movements of atoms
  !
  !     INPUT:     NAPDEV,   the permitted number of atoms on a device
  !                NPART,    the actual number of atoms on the device
  !                XP1,       velocity of atoms
  !                FP,        force on atoms
  !
  !     OUTPUT     XP1,   updated velocity

  !
  implicit none
  !----   DUMMY Variables
          integer,     value ::NAPDEV, NPART
          real(KINDDF),device::XP1(NAPDEV,3),FP(NAPDEV,3)
           integer,    device::STATU(*)

  !----   Local variables
          integer::K, IT, IB, IC, ITY, STAT
          real(KINDDF)::XP1x,XP1y, XP1z, FPx, FPy, FPz
  !

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC.LE.NPART) then
                 STAT = STATU(IC)

                 XP1x  = XP1(IC,1)
                 XP1y  = XP1(IC,2)
                 XP1z  = XP1(IC,3)

                 FPx  = FP(IC,1)
                 FPy  = FP(IC,2)
                 FPz  = FP(IC,3)

                 !---
                 if(XP1x*FPx  < C_ZERO) then
                    XP1x = C_ZERO
                 end if
                 if(XP1y*FPy  < C_ZERO) then
                    XP1y = C_ZERO
                 end if
                 if(XP1z*FPz  < C_ZERO) then
                    XP1z = C_ZERO
                 end if

                 if(IAND(STAT,CP_STATU_FIXPOSX) .eq. CP_STATU_FIXPOSX ) XP1x = C_ZERO
                 if(IAND(STAT,CP_STATU_FIXPOSY) .eq. CP_STATU_FIXPOSY ) XP1y = C_ZERO
                 if(IAND(STAT,CP_STATU_FIXPOSZ) .eq. CP_STATU_FIXPOSZ ) XP1z = C_ZERO

                 XP1(IC, 1) = XP1x
                 XP1(IC, 2) = XP1y
                 XP1(IC, 3) = XP1z
             end if


        return
  end subroutine DAMPING_KERNEL
  !************************************************************************************************

  !****************************************************************************************
  subroutine DAMPING_template(IDEV, STARTCELL, ENDCELL, XP1, FP, STATU)
  !***  PORPOSE:   the template to invoke the duming KERNEL for devices
  !
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                XP1,       velocity of atoms
  !                FP,        force on the atoms
  !     OUTPUT     XP1,       updated velocity
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,       intent(in)           ::IDEV, STARTCELL, ENDCELL
          real(KINDDF), device, dimension(:,:)::XP1, FP
          integer,      device, dimension(:)  ::STATU


  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1


           !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

           !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB  = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)

            call DAMPing_KERNEL<<<blocks, threads>>>(m_NAPDEV, NPRT,XP1,FP, STATU)

            ERR = cudaSetDevice(CURDEV)


         return
  end subroutine DAMPING_template
  !****************************************************************************************

  !************************************************************************************************
  attributes(global) subroutine Predictor_KERNEL0(IM,IA0,ITYP,XP,NAPDEV,NPART,XP1,FP,DIS,STATU,&
                                                 CM, TH,H2S2, HS2,                             &
                                                 LBX,LBY,LBZ,UBX,UBY,UBZ,BSX,BSY,BSZ)
  !***  PURPOSE:   to predict the postion and velocity of the atomss in the first half time step
  !
  !$$   INPUT:     IM,         total number of atoms in the whole box
  !$$              IA0,        starting atom on the device
  !$$              ITYP,       indicator of type of atoms
  !$$              XP,         position of atoms
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              FP,         force on atoms
  !$$              DIS,        the accumulated dispacement of atoms
  !$$              STATU,      statu ot atoms
  !$$              CM,         mass of atoms
  !$$              TH,         half time step
  !$$              H2S2
  !$$              HS2
  !
  !$$   OUTPUT     XP, XP1,     updated position and velocity
  !$$              STATU,,      statu of the atoms
  !
  implicit none
  !----   DUMMY Variables
          integer,      value ::IM, NAPDEV,NPART, KP, IA0, CHECKSTATU
          real(KINDDF), value ::TH,H2S2, HS2
          real(KINDDF), value ::LBX,LBY,LBZ,UBX,UBY,UBZ,BSX,BSY,BSZ
          integer,      device::ITYP(*)
          real(KINDDF), device::XP(IM,*),XP1(NAPDEV,*),FP(NAPDEV,*),DIS(NAPDEV,*)
          real(KINDDF), device::CM(*)
          integer,      device::STATU(*)
          

  !----   Local variables
          integer, shared::IFPDX,IFPDY,IFPDZ
          

          integer::K, IT, IB, IC, IC1, ITY, STAT
          real(KINDDF)::CM0,XPx, XPy, XPz,XP1x,XP1y, XP1z, FPx, FPy, FPz,DISx, DISy, DISz,DIS1, DIS2, DIS3

              IT    = (threadidx%y-1)*blockdim%x + threadidx%x
              IFPDX = dcm_IFPD(1)
              IFPDY = dcm_IFPD(2)
              IFPDZ = dcm_IFPD(3)

              IB    =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC    = (IB-1)*(blockdim%x*blockdim%y)+IT
              IC1   = IC + IA0

              if(IC.LE.NPART) then
                 STAT = STATU(IC)
                 if(IAND(STAT,CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                    XPx   = XP(IC1,1)
                    XPy   = XP(IC1,2)
                    XPz   = XP(IC1,3)
                    XP1x  = XP1(IC,1)
                    XP1y  = XP1(IC,2)
                    XP1z  = XP1(IC,3)
                    DIS1  = DIS(IC,1)
                    DIS2  = DIS(IC,2)
                    DIS3  = DIS(IC,3)

                    ITY =  ITYP(IC1)
                    CM0 =  CM(ITY)

                    FPx  = FP(IC,1)/CM0
                    FPy  = FP(IC,2)/CM0
                    FPz  = FP(IC,3)/CM0

                    !$$--- check periodic boundary condition
                    DISx = TH*XP1x + H2S2*FPx
                    if(IAND(STAT,CP_STATU_FIXPOSX) .eq. CP_STATU_FIXPOSX ) DISx = 0.D0
                    XPx = XPx + DISx
                    if(IFPDX) then
                       if(XPx .GT. UBX) then
                          XPx = XPx - BSX
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       else if (XPx .LT. LBX) then
                          XPx = XPx + BSX
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       end if
                    end if

                    DISy = TH*XP1y + H2S2*FPy
                    if(IAND(STAT,CP_STATU_FIXPOSY) .eq. CP_STATU_FIXPOSY ) DISy = 0.D0
                    XPy = XPy + DISy
                    if(IFPDY) then
                       if(XPy .GT. UBY) then
                          XPy = XPy - BSY
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       else if (XPy .LT. LBY) then
                          XPy = XPy + BSY
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       end if
                    end if

                    DISz = TH*XP1z + H2S2*FPz
                    if(IAND(STAT,CP_STATU_FIXPOSZ) .eq. CP_STATU_FIXPOSZ ) DISz = 0.D0
                    XPz = XPz + DISz
                    if(IFPDZ) then
                       if(XPz .GT. UBZ) then
                          XPz = XPz - BSZ
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       else if (XPz .LT. LBZ) then
                          XPz = XPz + BSZ
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       end if
                    end if

                    if(IAND(STAT,CP_STATU_FIXVELX) .eq. C_IZERO .and. &
                       IAND(STAT,CP_STATU_FIXPOSX) .eq. C_IZERO )     &
                       XP1x = XP1x + HS2*FPx
                    if(IAND(STAT,CP_STATU_FIXVELY) .eq. C_IZERO .and. &
                       IAND(STAT,CP_STATU_FIXPOSY) .eq. C_IZERO)      &
                       XP1y = XP1y + HS2*FPy
                    if(IAND(STAT,CP_STATU_FIXVELZ) .eq. C_IZERO .and. &
                       IAND(STAT,CP_STATU_FIXPOSZ) .eq. C_IZERO)      &
                       XP1z = XP1z + HS2*FPz

                    !$$--- copy back to global memory
                    XP(IC1, 1)  = XPx
                    XP(IC1, 2)  = XPy
                    XP(IC1, 3)  = XPz
                    XP1(IC, 1)     = XP1x
                    XP1(IC, 2)     = XP1y
                    XP1(IC, 3)     = XP1z

                    DIS(IC, 1)     = DIS1 + DISx
                    DIS(IC, 2)     = DIS2 + DISy
                    DIS(IC, 3)     = DIS3 + DISz

                    if(XPx .GT. UBX .OR. XPy .GT. UBY .OR.  XPz .GT. UBZ ) then
                       STATU(IC) = IOR(CP_STATU_OUTOFBOX,CP_STATU_REFLECT)
                    else if(XPx .LT. LBX  .OR. XPy .LT. LBY .OR. XPz .LT. LBZ ) then
                       STATU(IC) = IOR(CP_STATU_OUTOFBOX, CP_STATU_TRANSMIT)
                    end if
                 end if
              end if

        return
  end subroutine Predictor_KERNEL0
  !****************************************************************************************

  !************************************************************************************************
  attributes(global) subroutine Predictor_KERNEL1(IM,IA0,ITYP,XP,NAPDEV,NPART,XPSTP,XP1,FP,DIS,STATU,&
                                                 CM, TH,H2S2, HS2,                                  &
                                                 LBX,LBY,LBZ,UBX,UBY,UBZ,BSX,BSY,BSZ)
  !***  PURPOSE:   to predict the postion and velocity of the atomss in the first half time step
  !
  !$$   INPUT:     IM,         total number of atoms in the whole box
  !$$              IA0,        starting atom on the device
  !$$              ITYP,       indicator of type of atoms
  !$$              XP,         position of atoms
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              FP,         force on atoms
  !$$              DIS,        the accumulated dispacement of atoms
  !$$              STATU,      statu ot atoms
  !$$              CM,         mass of atoms
  !$$              TH,         half time step
  !$$              H2S2
  !$$              HS2
  !
  !$$   OUTPUT     XP, XP1,     updated position and velocity
  !$$              STATU,,      statu of the atoms
  !
  implicit none
  !----   DUMMY Variables
          integer,      value ::IM, NAPDEV,NPART, KP, IA0, CHECKSTATU
          real(KINDDF), value ::TH,H2S2, HS2
          real(KINDDF), value ::LBX,LBY,LBZ,UBX,UBY,UBZ,BSX,BSY,BSZ
          integer,      device::ITYP(*)
          real(KINDDF), device::XP(IM,*),XPSTP(NAPDEV,*), XP1(NAPDEV,*),FP(NAPDEV,*),DIS(NAPDEV,*)
          real(KINDDF), device::CM(*)
          integer,      device::STATU(*)
          

  !----   Local variables
          integer, shared::IFPDX,IFPDY,IFPDZ
          

          integer::K, IT, IB, IC, IC1, ITY, STAT
          real(KINDDF)::CM0,XPx, XPy, XPz,XP1x,XP1y, XP1z, FPx, FPy, FPz,DISx, DISy, DISz,DIS1, DIS2, DIS3

              IT    = (threadidx%y-1)*blockdim%x + threadidx%x
              IFPDX = dcm_IFPD(1)
              IFPDY = dcm_IFPD(2)
              IFPDZ = dcm_IFPD(3)

              IB    =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC    = (IB-1)*(blockdim%x*blockdim%y)+IT
              IC1   = IC + IA0

              if(IC.LE.NPART) then
                 STAT = STATU(IC)
                 if(IAND(STAT,CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                    XPx   = XP(IC1,1)
                    XPy   = XP(IC1,2)
                    XPz   = XP(IC1,3)
                    XP1x  = XP1(IC,1)
                    XP1y  = XP1(IC,2)
                    XP1z  = XP1(IC,3)
                    DIS1  = DIS(IC,1)
                    DIS2  = DIS(IC,2)
                    DIS3  = DIS(IC,3)

                    ITY =  ITYP(IC1)
                    CM0 =  CM(ITY)

                    FPx  = FP(IC,1)/CM0
                    FPy  = FP(IC,2)/CM0
                    FPz  = FP(IC,3)/CM0

                    !$$--- check periodic boundary condition
                    DISx = TH*XP1x + H2S2*FPx
                    if(IAND(STAT,CP_STATU_FIXPOSX) .eq. CP_STATU_FIXPOSX ) DISx = 0.D0
                    XPx = XPx + DISx
                    if(IFPDX) then
                       if(XPx .GT. UBX) then
                          XPx = XPx - BSX
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       else if (XPx .LT. LBX) then
                          XPx = XPx + BSX
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       end if
                    end if

                    DISy = TH*XP1y + H2S2*FPy
                    if(IAND(STAT,CP_STATU_FIXPOSY) .eq. CP_STATU_FIXPOSY ) DISy = 0.D0
                    XPy = XPy + DISy
                    if(IFPDY) then
                       if(XPy .GT. UBY) then
                          XPy = XPy - BSY
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       else if (XPy .LT. LBY) then
                          XPy = XPy + BSY
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       end if
                    end if

                    DISz = TH*XP1z + H2S2*FPz
                    if(IAND(STAT,CP_STATU_FIXPOSZ) .eq. CP_STATU_FIXPOSZ ) DISz = 0.D0
                    XPz = XPz + DISz
                    if(IFPDZ) then
                       if(XPz .GT. UBZ) then
                          XPz = XPz - BSZ
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       else if (XPz .LT. LBZ) then
                          XPz = XPz + BSZ
                          STAT = IOR(STAT, CP_STATU_PASSBOUND)
                       end if
                    end if

                    if(IAND(STAT,CP_STATU_FIXVELX) .eq. C_IZERO .and. &
                       IAND(STAT,CP_STATU_FIXPOSX) .eq. C_IZERO )     &
                       XP1x = XP1x + HS2*FPx
                    if(IAND(STAT,CP_STATU_FIXVELY) .eq. C_IZERO .and. &
                       IAND(STAT,CP_STATU_FIXPOSY) .eq. C_IZERO)      &
                       XP1y = XP1y + HS2*FPy
                    if(IAND(STAT,CP_STATU_FIXVELZ) .eq. C_IZERO .and. &
                       IAND(STAT,CP_STATU_FIXPOSZ) .eq. C_IZERO)      &
                       XP1z = XP1z + HS2*FPz

                    !$$--- copy back to global memory
                    XP(IC1, 1)  = XPx
                    XP(IC1, 2)  = XPy
                    XP(IC1, 3)  = XPz
                    XP1(IC, 1)     = XP1x
                    XP1(IC, 2)     = XP1y
                    XP1(IC, 3)     = XP1z

                    XPSTP(IC, 1)   = DISx
                    XPSTP(IC, 2)   = DISy
                    XPSTP(IC, 3)   = DISz

                    DIS(IC, 1)     = DIS1 + DISx
                    DIS(IC, 2)     = DIS2 + DISy
                    DIS(IC, 3)     = DIS3 + DISz

                    if(XPx .GT. UBX .OR. XPy .GT. UBY .OR.  XPz .GT. UBZ ) then
                       STATU(IC) = IOR(CP_STATU_OUTOFBOX,CP_STATU_REFLECT)
                    else if(XPx .LT. LBX  .OR. XPy .LT. LBY .OR. XPz .LT. LBZ ) then
                       STATU(IC) = IOR(CP_STATU_OUTOFBOX, CP_STATU_TRANSMIT)
                    end if
                 end if
              end if

        return
  end subroutine Predictor_KERNEL1
  !****************************************************************************************

  !****************************************************************************************
  subroutine PREDICTOR_template(IDEV, STARTCELL, ENDCELL, CM, ITYP,XP, XP1, FP, DIS, STATU, XPSTP)
  !***  PORPOSE:   template to invoke the predictor KERNEL
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                CM,         mass of atoms
  !                ITYP,       type id of atoms
  !
  !     OUTPUT     XP,XP1,DIS,STATU
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,intent(in)                  ::IDEV, STARTCELL, ENDCELL
          real(KINDDF), device, dimension(:)  ::CM
          integer,      device, dimension(:)  ::ITYP, STATU
          real(KINDDF), device, dimension(:,:)::XP, XP1, FP, DIS
          real(KINDDF), device, dimension(:,:),optional::XPSTP


  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV, IA0

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            !dcm_LBOX(1:3)    = m_LBOX(1:3)
            !dcm_UBOX(1:3)    = m_UBOX(1:3)
            !dcm_BOXSIZE(1:3) = m_BOXSIZE(1:3)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1

            !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            IA0     = STARTA-1
            if(present(XPSTP)) then
                call Predictor_KERNEL1<<<blocks, threads>>>(dm_NPRT, IA0, ITYP, XP, m_NAPDEV, NPRT, XPSTP, XP1, FP, DIS, STATU, &
                                                       CM, m_TH, m_H2S2, m_HS2,                                                 &
                                                       m_LBOX(1), m_LBOX(2), m_LBOX(3), m_UBOX(1), m_UBOX(2), m_UBOX(3),        &
                                                       m_BOXSIZE(1), m_BOXSIZE(2), m_BOXSIZE(3))
            else 
               call Predictor_KERNEL0<<<blocks, threads>>>(dm_NPRT, IA0, ITYP, XP, m_NAPDEV, NPRT, XP1, FP, DIS, STATU,         &
                                                       CM, m_TH, m_H2S2, m_HS2,                                                 &
                                                       m_LBOX(1), m_LBOX(2), m_LBOX(3), m_UBOX(1), m_UBOX(2), m_UBOX(3),        &
                                                       m_BOXSIZE(1), m_BOXSIZE(2), m_BOXSIZE(3))
            end if          
            ERR    = cudaSetDevice(CURDEV)
     return
  end subroutine PREDICTOR_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine Predictor_DEV(ITIME,SimBox,CtrlParam)
  !***  PORPOSE:  to forward one half of a time step
  !
  !     INPUT:     ITIME,     the time step
  !                SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !     OUTPUT     dxm_XP, dxm_XP1, the positions and velocities of atoms on devices
  !                hm_XP,           the positions of atoms on host
  !
  use MD_Globle_Variables_GPU, only:dm_WorkSpace, m_STARTCELL, m_ENDCELL, Synchroniz_XP_on_Devices

  implicit none
  !----   DUMMY Variables
          integer,intent(in)::ITIME
          type(SimMDBox)    ::SimBox
          type(SimMDCtrl)   ::CtrlParam
  !---   Local variables
          integer::IT, IVSECT, IFLAG, I


         !$$--- if damping required, do it
         if(ITIME .ge. CtrlParam%DAMPTIME0 .and.   &
            ITIME .le. CtrlParam%DAMPTIME0 + CtrlParam%DAMPTIME1-1) then
              do I=1, m_NDEVICE
                 call DAMPING_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), &
                                       dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%FP(I)%Data, dm_WorkSpace%STATU(I)%Data)
              end do 
         end if

       !$$--- if variable time step (scheme II) is used, need to check and adjust the time step
           if(CtrlParam%IHDUP .LT. 0) then
              !if(MOD(ITIME-CtrlParam%IT0+1,IABS(CtrlParam%IHDUP)).EQ. C_IUN .OR. &
              !   IABS(CtrlParam%IHDUP).EQ. C_IUN) then
              if(MOD(ITIME-CtrlParam%IT0+1,IABS(CtrlParam%IHDUP)).EQ. C_ZERO) then
                m_MXD2 = CtrlParam%DMX*CtrlParam%DMX
                m_TH =  CtrlParam%HMX
                m_HS2  = m_TH*C_HALF
                m_H2S2 = m_TH*m_TH*C_HALF
                do while(.TRUE.)
                   call CheckTimestep_DEV(ITIME,SimBox,CtrlParam, m_TH, m_H2S2, m_MXD2, IFLAG)

                   if(IFLAG .gt. 0) then
                      m_TH   = m_TH*C_HALF
                      m_HS2  = m_TH*C_HALF
                      m_H2S2 = m_TH*m_TH*C_HALF
                   else
                      !--- Save the current time step
                      CtrlParam%H = m_TH
                      exit
                   end if

                end do
              end if
           end if

       !$$--- update the time step
           m_TH   = CtrlParam%H
           m_HS2  = m_TH*C_HALF
           m_H2S2 = m_TH*m_TH*C_HALF

       !$$--- update boxsize
           m_LBOX(1:3) = SimBox%BOXLOW(1:3)
           m_UBOX(1:3) = SimBox%BOXUP(1:3)
           m_BOXSIZE(1:3) = SimBox%ZL(1:3)

       !$$--- predicte the position and velocity
            do I=1, m_NDEVICE
               call PREDICTOR_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,         &
                                      dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP(I)%Data, dm_WorkSpace%XP1(I)%Data, &
                                      dm_WorkSpace%FP(I)%Data, dm_WorkSpace%DIS(I)%Data,  dm_WorkSpace%STATU(I)%Data)
                                      !dm_WorkSpace%XPSTP(I)%Data)                                      
            end do

       !$$--- hm_XP has bee modified in PREDICTOR_template, we should sychonize
       !$$    XP on devices for the next force calculations
            call Synchroniz_XP_on_Devices()

     return
  end subroutine Predictor_DEV
  !****************************************************************************************

  !************************************************************************************************
  attributes(global) subroutine Correction_KERNEL(IM,IA0,ITYP,NAPDEV,NPART,XP1,FP,STATU, CM, TH,H2S2, HS2)
  !***  PURPOSE:   to to correct the velocity of the atoms in the second half time step
  !
  !$$   INPUT:     IM,         total number of atoms in the whole box
  !$$              IA0,        starting atom on the device
  !$$              ITYP,       indicator of type of atoms
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              FP,         force on atoms
  !$$              DIS,        accumulated dispacement of atoms
  !$$              STATU,      statu ot atoms
  !$$              CM,         mass of atoms
  !$$              TH,         half time step
  !$$              H2S2
  !$$              HS2
  !
  !$$   OUTPUT     XP1,      updated velocity
  !
  implicit none
  !----   DUMMY Variables
          integer,     value::IM, NAPDEV,NPART, IA0
          integer,     device::ITYP(*), STATU(*)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*),CM(*)
          real(KINDDF), value::TH, H2S2, HS2


  !----   Local variables
          integer::IT, IB, IC, ITY, STAT
          real(KINDDF)::XP1x,XP1y, XP1z, FPx, FPy, FPz, CM0
  !
  !***
  !
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC .LE. NPART) then
                 STAT = STATU(IC)

                if(IAND(STAT,CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then

                 XP1x  = XP1(IC,1)
                 XP1y  = XP1(IC,2)
                 XP1z  = XP1(IC,3)
                 ITY   = ITYP(IC+IA0)
                 CM0   = CM(ITY)

                 FPx  = FP(IC,1)/CM0
                 FPy  = FP(IC,2)/CM0
                 FPz  = FP(IC,3)/CM0

                 if(IAND(STAT,CP_STATU_FIXVELX) .eq. C_IZERO .and. &
                    IAND(STAT,CP_STATU_FIXPOSX) .eq. C_IZERO )     &
                    XP1x = XP1x + HS2*FPx

                 if(IAND(STAT,CP_STATU_FIXVELY) .eq. C_IZERO .and. &
                    IAND(STAT,CP_STATU_FIXPOSY) .eq. C_IZERO )     &
                    XP1y = XP1y + HS2*FPy

                 if(IAND(STAT,CP_STATU_FIXVELZ) .eq. C_IZERO .and. &
                    IAND(STAT,CP_STATU_FIXPOSZ) .eq. C_IZERO )     &
                    XP1z = XP1z + HS2*FPz

                 XP1(IC,1) = XP1x
                 XP1(IC,2) = XP1y
                 XP1(IC,3) = XP1z

                end if
              end if
         return
  end subroutine Correction_KERNEL
  !******************************************************************************************


  !****************************************************************************************
  subroutine CORRECTION_template(IDEV, STARTCELL, ENDCELL, CM, ITYP, XP1, FP, STATU)
  !***  PORPOSE:   template to invoke the correction KERNEL on devices
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                CM,         mass of atoms
  !                ITYP,       type id of atoms
  !                XP1,        velocities of atoms on device IDEV
  !                FP,         force on atoms on  device IDEV
  !                STATU,      statu of atoms on  device IDEV
  !
  !     OUTPUT     XP1,       updated velocities
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,intent(in)::IDEV, STARTCELL, ENDCELL
          integer, device, dimension(:)::ITYP, STATU
          real(KINDDF), device, dimension(:,:)::XP1, FP
          real(KINDDF), device, dimension(:):: CM


  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1


            !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call Correction_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, ITYP, m_NAPDEV, NPRT, XP1, FP, STATU,&
                                                        CM, m_TH, m_H2S2,m_HS2)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine CORRECTION_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine Correction_DEV(ITIME,SimBox,CtrlParam)
  !***  PORPOSE: to forward one time step
  !
  !     INPUT:     ITIME,     the time step
  !                SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !     OUTPUT     SimBox,    the simulation box with new velocity
  !
  use MD_Globle_Variables_GPU, only:dm_WorkSpace, m_STARTCELL, m_ENDCELL

  implicit none
  !----   DUMMY Variables
          integer,intent(in)::ITIME
          type(SimMDBox)    ::SimBox
          type(SimMDCtrl)   ::CtrlParam
  !---    Local varibales
          integer::I


            !$$--- to do the correction
            do I=1, m_NDEVICE
               call CORRECTION_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,         &
                                       dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%FP(I)%Data, &
                                       dm_WorkSpace%STATU(I)%Data)
            end do
     return
  end subroutine Correction_DEV
  !****************************************************************************************

  !************************************************************************************************
  attributes(global) subroutine CALEKIN_KERNEL(IM,IA0,ITYP,NAPDEV,NPART,XP1,STATU, CM, EKIN)
  !***  PURPOSE:   to calculate the kinetic energies of atoms
  !
  !     INPUT:     IM,         total number of atoms in the whole box
  !                IA0,        starting atom on the device
  !                ITYP,       indicator of type of atoms
  !                NAPDEV,     permitted number of atoms on the device
  !                NPART,      actuall number of atoms on the device
  !                XP1,        velocity of atoms
  !                STATU,      statu ot atoms
  !                CM,         mass of atoms
  !
  !     OUTPUT     EKIN,       kinetic energies of atoms
  !

  !
  implicit none
  !----   DUMMY Variables
          integer, value::IM, NAPDEV,NPART, IA0
          integer,device::ITYP(*), STATU(*)
          real(KINDDF),device::XP1(NAPDEV,*),CM(*), EKIN(*)


  !----   Local variables
          integer::IT, IB, IC,ITY
          real(KINDDF)::XP1x,XP1y, XP1z, CM0, EK
  !
  !***
  !
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC .LE. NPART) then
                EK = -1.D32
                if(IAND(STATU(IC),CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE .and. &
                   IAND(STATU(IC),CP_STATU_FIXPOS)  .eq. C_ZERO ) then

                 ITY =  ITYP(IC+IA0)
                 CM0 = CM(ITY)

                 XP1x = XP1(IC,1)
                 XP1y = XP1(IC,2)
                 XP1z = XP1(IC,3)
                 EK   = C_HALF*CM0*(XP1x*XP1x+XP1y*XP1y+XP1z*XP1z)
               end if
               EKIN(IC) = EK

             end if
         return
  end subroutine CALEKIN_KERNEL
  !******************************************************************************************

  !****************************************************************************************
  subroutine CALEKIN_template(IDEV, STARTCELL, ENDCELL, CM, ITYP, XP1, STATU, EKIN)
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
           integer,intent(in)                  ::IDEV, STARTCELL, ENDCELL
           integer,      device, dimension(:)  ::ITYP, STATU
           real(KINDDF), device, dimension(:,:)::XP1
           real(KINDDF), device, dimension(:)  :: CM
           real(KINDDF), device, dimension(:)  :: EKIN


  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1


            !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call CALEKIN_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, ITYP, m_NAPDEV, NPRT, XP1, STATU, CM, EKIN)
            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine CALEKIN_template
  !****************************************************************************************

  !*********************************************************************************
  subroutine COPYOUTEKIN_template(IDEV, STARTCELL, ENDCELL, dEKIN, hEKIN)
  !***  PURPOSE:   to copy the EKIN array on a device to host array
  !
  !
  !     INPUT:     IDEV,                ID of device
  !                STARTCELL, ENDCELL,  ID of the first and last cell on device IDEV
  !                dEKIN,               kinetic energy array on device IDEV
  !
  !     OUTPUT     hEKIN                kinetic energy array on host
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      real(KINDDF), device, dimension(:), allocatable::dEKIN
      real(KINDDF),         dimension(:), allocatable::hEKIN


      !--- Local variables
         integer::ERR, STARTA, ENDA, NPRT, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT     = ENDA - STARTA + 1
            !$$--- NOTE: the subscript of dEKIN is the index of atom on device
            !$$hEKIN(STARTA:ENDA) = dEKIN(1:NPRT)

            ERR = cudaMemcpyAsync(hEKIN(STARTA),dEKIN(1), NPRT)
            ERR = cudaSetDevice(CURDEV)
          return
  end subroutine COPYOUTEKIN_template
  !*********************************************************************************

  !*********************************************************************************
  subroutine CALEKIN_DEV(SimBox,CtrlParam, dbgEKIN)
  !***  PORPOSE:   to calculate the kinetic energy
  !
  !     INPUT:     SimBox,    the simulation box,   not used
  !                CtrlParam, the control parameters, not used
  !     OUTPUT     dm_EKIN,   kinetic eenrgy on device0
  !                dbgEKINE,  kinetic energy array on host, optional, for debug purpose
  !
  use MD_Globle_Variables_GPU, only:dm_WorkSpace, m_STARTCELL, m_ENDCELL, hm_EKIN, hm_GIDINV

  implicit none
  !----   DUMMY Variables
          type(SimMDBox)                      ::SimBox
          type(SimMDCtrl)                     ::CtrlParam
          real(KINDDF), dimension(:), optional::dbgEKIN

  !--- Local variables
          integer::I
  
            do I=1, m_NDEVICE       
               call CALEKIN_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,             &
                                     dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%STATU(I)%Data, &
                                     dm_WorkSpace%EKIN(I)%Data)
            end do

            hm_EKIN = -1.D32
            do I=1, m_NDEVICE
               call COPYOUTEKIN_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),  &
                                         dm_WorkSpace%EKIN(I)%Data, hm_EKIN)
            end do
            call SynchronizeDevices()

             !$$--- copy the EKIN in the global array of EKIN
             if(present(dbgEKIN)) then
                  dbgEKIN = hm_EKIN(hm_GIDINV)
             end if

     return
  end subroutine CALEKIN_DEV
  !****************************************************************************************

  !*********************************************************************************
  subroutine Cal_GlobalT_DEV(SimBox,CtrlParam, CURT)
  !***  PORPOSE:   to calculate the global temperature
  !                NOTE:        the temperature is the average tempertaure of the all boxes
  !
  !     INPUT:     SimBox,     the simulation box,   not used
  !                CtrlParam,  the control parameters, not used
  !     OUTPUT     CURT,       current temperature
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          type(SimMDBox)  :: SimBox
          type(SimMDCtrl) :: CtrlParam
          real(KINDDF)    :: CURT

  !--- Local variables
  !--- Device variables and variables to be used in GPU
             call CALEKIN_DEV(SimBox,CtrlParam)

             CURT = C_TWO*sum(hm_EKIN, mask=(hm_EKIN .ge. 0))/dble(count(hm_EKIN .ge. 0.D0))/(C_THR*CP_KB)
     return
  end subroutine Cal_GlobalT_DEV
  !****************************************************************************************

  !************************************************************************************************
  attributes(global) subroutine CheckTimestep_KERNEL(IM,IA0,ITYP,NAPDEV,NPART,XP1,FP,STATU,CM, &
                                                 TH,H2S2,MXD2,IFLAG)
  !***  PURPOSE:   to check the if any displacement to be larger than given value
  !
  !$$   INPUT:     IM,         total number of atoms in the whole box
  !$$              IA0,        starting atom on the device
  !$$              ITYP,       indicator of type of atoms
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              FP,         force on atoms
  !$$              STATU,      statu ot atoms
  !$$              CM,         mass of atoms
  !$$              TH,         half time step
  !$$              H2S2
  !$$              MXD2,       max permitted displacement
  !
  !$$   OUTPUT     IFLAG,     indicator to indicating if need adjust step size
  !
  implicit none
  !----   DUMMY Variables
          integer, value::IM, NAPDEV,NPART, KP, IA0
          logical,device::IFLAG(*)
          integer,device::ITYP(*)
          real(KINDDF),device::XP1(NAPDEV,*),FP(NAPDEV,*)
          real(KINDDF),device::CM(*)
          integer,device::STATU(*)
          real(KINDDF), value::TH, H2S2, MXD2

  !----   Local variables
          integer::K, IT, IB, IC, ITY
          real(KINDDF)::CM0,XP1x,XP1y, XP1z, FPx, FPy, FPz,DISx, DISy, DISz

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              !--- already konw that the time step is too large
              if(IFLAG(1)) return

               DISx = 0.D0
               DISy = 0.D0
               DISz = 0.D0

              !--- check if the displacemnet is larger than give value.
              if(IC.LE.NPART) then
              if(IAND(STATU(IC),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE ) then
                 XP1x  = XP1(IC,1)
                 XP1y  = XP1(IC,2)
                 XP1z  = XP1(IC,3)

                 ITY =  ITYP(IC+IA0)
                 CM0 =  CM(ITY)

                 FPx  = FP(IC,1)/CM0
                 FPy  = FP(IC,2)/CM0
                 FPz  = FP(IC,3)/CM0

                 if(IAND(STATU(IC),CP_STATU_FIXPOSX) .eq. C_IZERO ) &
                    DISx = TH*XP1x + H2S2*FPx
                 if(IAND(STATU(IC),CP_STATU_FIXPOSY) .eq. C_IZERO ) &
                    DISy = TH*XP1y + H2S2*FPy
                 if(IAND(STATU(IC),CP_STATU_FIXPOSZ) .eq. C_IZERO ) &
                    DISz = TH*XP1z + H2S2*FPz

                 if(DISx*DISx+DISy*DISy+DISz*DISz .gt. MXD2) IFLAG(1) = .true.
              end if
              end if

        return
  end subroutine CheckTimestep_KERNEL
  !************************************************************************************************

  !****************************************************************************************
  subroutine CheckTimestep_template(IDEV, STARTCELL, ENDCELL, CM, ITYP, XP1, FP, STATU, &
                                    TH, H2S2, DMX2,dIFlag)
  !***  PORPOSE:   template to invoke the CheckTimestep_KERNEL
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                CM,         mass of atoms
  !                ITYP,       type id of atoms
  !                XP1,        velocity of atoms
  !                FP,         force of atoms
  !                STATU,      state of atoms
  !                TH,         current time step
  !                H2S2,       
  !                HS2,        
  !                DMX2,       square of the max displacement in a time step  
  !
  !     OUTPUT     dIFLAG,    flag indicating if there is any atom wiht estimated displacement larger
  !                           than permitted value
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
           integer,intent(in)                  ::IDEV, STARTCELL, ENDCELL
           real(KINDDF), intent(in)            ::TH, H2S2, DMX2 
           integer,      device, dimension(:)  ::ITYP, STATU
           real(KINDDF), device, dimension(:,:)::XP1, FP
           real(KINDDF), device, dimension(:)  :: CM
           logical,      device, dimension(:)  ::dIFLAG

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            !dcm_LBOX(1:3)    = m_LBOX(1:3)
            !dcm_UBOX(1:3)    = m_UBOX(1:3)
            !dcm_BOXSIZE(1:3) = m_BOXSIZE(1:3)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1


            !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB        = (NPRT-1)/(BX*BY)+1
            blocks    = dim3(NB, 1, 1)
            threads   = dim3(BX, BY, 1)
            dIFLAG(1) = .false.
            STARTA    = STARTA-1 
            call CheckTimestep_KERNEL<<<blocks, threads>>>(dm_NPRT, STARTA, ITYP, m_NAPDEV, NPRT, XP1, FP, STATU,&
                                                           CM, TH, H2S2, DMX2, dIFLAG)

            ERR = cudaSetDevice(CURDEV)


     return
  end subroutine CheckTimestep_template
  !****************************************************************************************

  !****************************************************************************************
  subroutine CheckTimestep_DEV(ITIME,SimBox,CtrlParam, TH, H2S2, DMX2, Iflag)
  !***  PORPOSE:  to check if there is atom with its displacement larger than permitted value
  !
  !     INPUT:     ITIME,     the time step
  !                SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !                DM2,       the square of the permiited displacemdent
  !
  !     OUTPUT     IFlah,     the flag indicating if any atom with its displacement
  !                           larger than the permitted value
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,         intent(in)::ITIME
          type(SimMDBox),  intent(in)::SimBox
          type(SimMDCtrl)            ::CtrlParam
          real(KINDDF),    intent(in)::TH, H2S2, DMX2
          integer                    ::IFlag 
  !---   Local variables
          integer::ERR, CURDEV, I
          logical::hIFLAG(m_MXDEVICE)=.false.

             ERR = cudaGetDevice(CURDEV)

                 do I=1,  m_NDEVICE
                    call CheckTimestep_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,                          &
                                    dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%FP(I)%Data, dm_WorkSpace%STATU(I)%Data, &
                                    TH, H2S2, DMX2, dm_DiffWS(I)%IFlag)
                 end do
                 hIFLAG = .false.
                 IFlag =  0
                 call SynchronizeDevices()

                 do I=1,  m_NDEVICE
                    ERR = cudaSetDevice(m_DEVICES(I))
                    hIFLAG(I:I) = dm_DiffWS(I)%IFlag(1:1)
                    if(hIFLAG(I)) IFlag = 1
                 end do

            ERR = cudaSetDevice(CURDEV)
          return 
     return
  end subroutine CheckTimestep_DEV
  !****************************************************************************************

  !************************************************************************************************
  attributes(global) subroutine VelScaling_KERNEL(IA0, NAPDEV, NPART, XP1, STATU, SCAL, NPRTPB)
  !***  PURPOSE:   to correct the velocity of the atoms by multiply a facter to velocities
  !                of atoms
  !
  !$$   INPUT:     
  !$$              IA0,        starting atom on the device
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              STATU,      statu ot atoms
  !$$              SCAL,       the scaling factor
  !$$
  !$$   OUTPUT     XP1,      updated velocity
  !
  implicit none
  !----   DUMMY Variables
          integer,     value  ::IA0, NAPDEV,NPART, NPRTPB
          real(KINDDF),device ::XP1(NAPDEV,*)
          integer,     device ::STATU(*)
          real(KINDDF),device ::SCAL(*)

  !----   Local variables
          integer::IT, IB, IC, IB0
          real(KINDDF)::XP1x,XP1y, XP1z, TSCAL

  !
  !***
  !
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC .LE. NPART) then
                if(IAND(STATU(IC),CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE ) then
                   XP1x = 0.D0
                   XP1y = 0.D0
                   XP1z = 0.D0
                   IB0 = (IA0 + IC-1)/NPRTPB + 1
                   TSCAL = SCAL(IB0)

                   if(IAND(STATU(IC),CP_STATU_FIXPOSX) .eq. C_IZERO .and. &
                      IAND(STATU(IC),CP_STATU_FIXVELX) .eq. C_IZERO) then
                      XP1x  = XP1(IC,1)*TSCAL
                   end if

                   if(IAND(STATU(IC),CP_STATU_FIXPOSY) .eq. C_IZERO .and. &
                      IAND(STATU(IC),CP_STATU_FIXVELY) .eq. C_IZERO) then
                       XP1y  = XP1(IC,2)*TSCAL
                   end if

                   if(IAND(STATU(IC),CP_STATU_FIXPOSZ) .eq. C_IZERO .and. &
                      IAND(STATU(IC),CP_STATU_FIXVELZ) .eq. C_IZERO) then
                      XP1z  = XP1(IC,3)*TSCAL
                   end if

                   XP1(IC,1) = XP1x
                   XP1(IC,2) = XP1y
                   XP1(IC,3) = XP1z
                end if
             end if
         return
  end subroutine VelScaling_KERNEL
  !******************************************************************************************

  !****************************************************************************************
  subroutine VelScaling_template(IDEV, STARTCELL, ENDCELL, XP1, STATU, hSCAL, dSCAL, NPRTPB)
  !***  PORPOSE:   template to invoke VelScaling_KERNEL on devices
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                XP1,        velocities of atoms on device IDEV
  !                STATU,      statu of atoms on  device IDEV
  !                hSCAL,      scal factor on devices
  !                NPRTPB,     number of atoms in each box
  !     TEMP:
  !                dSCAL,      scal factor on devices
  !
  !     OUTPUT     XP1,        velocities updated
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)              ::IDEV, STARTCELL, ENDCELL
          integer,      device, dimension(:)    ::ITYP, STATU
          real(KINDDF), device, dimension(:,:)  ::XP1
          real(KINDDF), dimension(:), intent(in)::hSCAL
          real(KINDDF), device, dimension(:)    ::dSCAL
          integer::NPRTPB

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV, MULTBOX

            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1


            !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB      = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)

            MULTBOX = size(hSCAL)
            STARTA  = STARTA-1
            ERR     =  cudaMemcpyAsync(dSCAL, hSCAL, MULTBOX)
            call VelScaling_KERNEL<<<blocks, threads>>>(STARTA, m_NAPDEV, NPRT,XP1,STATU, dSCAL, NPRTPB)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine VelScaling_template
 !*********************************************************************************

 !*********************************************************************************
  subroutine VelScaling_DEV(SimBox,CtrlParam, DT)
  !***  PORPOSE:   to scaling the temperature of the system to
  !                given value
  !
  !     INPUT:     SimBox,    the simulation box,   not used
  !                CtrlParam, the control parameters, not used
  !                DT,        the termperature (in K) to scaled to
  !     OUTPUT     dm_XP1,    velocities of atoms on devices
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          type(SimMDBox) :: SimBox
          type(SimMDCtrl):: CtrlParam
          real(KINDDF)   :: DT
  !--- Local variables
          integer I, J, NPRT, IPA, ERR, CURDEV, MULTIBOX

  !--- Device variables and variables to be used in GPU
            real(KINDDF), dimension(:),allocatable::cEKIN

            ERR = cudaGetDevice(CURDEV)

            call CALEKIN_DEV(SimBox,CtrlParam)
            NPRT     = SimBox%NPRT
            MULTIBOX = dm_NPRT/NPRT
            allocate(cEKIN(MULTIBOX))

            !$$--- to get the current temperature (in erg)
             IPA = 0
             do I=1, MULTIBOX
                cEKIN(I) = sum(hm_EKIN(hm_GIDINV(IPA+1:IPA+NPRT)), mask=(hm_EKIN(hm_GIDINV(IPA+1:IPA+NPRT)) .ge. 0.D0))/ &
                           dble(count(hm_EKIN(hm_GIDINV(IPA+1:IPA+NPRT)) .ge. 0.D0))

               IPA = IPA + NPRT
             end do
             !cEKIN = sum(hEKIN(1:dm_NPRT), mask=(hEKIN(1:dm_NPRT) .ge. 0.D0))/dble(count(hEKIN .ge. 0.D0))
             if(any(cEKIN .le. 0.D0)) then
                write(*,fmt="(A)")          " MDPSCU Error: current temperature of the system is zero in scaling velocity"
                write(*,fmt="(A, 1pE12.4)") "               current temperature(erg) = ", cEKIN
                write(*,fmt="(A)")          "               Process to be stopped"
                stop
             end if
             cEKIN = dsqrt(DT*C_THR*CP_KB*C_HALF/cEKIN)

             do I=1, m_NDEVICE
                call VelScaling_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%STATU(I)%Data, &
                                         cEKIN, dm_WorkSpace%EKIN(I)%Data, NPRT)
             end do
              !--- NOTE: here, we use dm_WorkSpace(I)%EKIN only for temperoray working spwce 
              !            dm_WorkSpace(I)%EKIN has been changed when call VelScaling_template.
              !          Each time we need EKIN, we should call CALEKIN_DEV 
             deallocate(cEKIN)

     return
  end subroutine VelScaling_DEV
  !********************************************************************************

  #ifdef NODEVRAN
  !************************************************************************************************
  attributes(global) subroutine Thermalizing_MC_KERNEL(IA0,ITYP,NAPDEV,NPART,XP1,STATU,CM, TI,Z1,Z2)
  !***  PURPOSE:   to scal the velocitie of Gaussin ditribution of unit half width to
  !                Maxswell distribution of given temperature
  !
  !$$   INPUT:     
  !$$              IA0,        starting atom on the device
  !$$              ITYP,       indicator of type of atoms
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              STATU,      statu ot atoms
  !$$              CM,         mass of atoms
  !$$              TI,         tmeprature to thermalize to
  !
  !$$   OUTPUT     XP1,       updated velocity
  !
  implicit none
  !----   DUMMY Variables
          integer,      value::NAPDEV,NPART, IA0
          integer,      device::ITYP(*), STATU(*)
          real(KINDDF), device::XP1(NAPDEV,*),CM(*)
          real(KINDDF), value::TI
          real(KINDDF), device::Z1(NPART,*),Z2(NPART,*)

  !----   Local variables
          integer::IT, IB, IC, STAT
          real(KINDDF)::V0
  !
  !***
  !
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC .LE. NPART) then
                STAT = STATU(IC)

                if(IAND(STAT,CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                  !$$--- NOTE: the stderr of Gaussian distribution ise set to 1.d0
                  V0 = DSQRT( C_TWO*TI*CP_KB/CM(ITYP(IC+IA0)))
                  !V0 = DSQRT( TI*CP_KB/CM(ITYP(IC+IA0)))

                   if(IAND(STAT,CP_STATU_FIXVELX) .eq. C_IZERO .and. &
                      IAND(STAT,CP_STATU_FIXPOSX) .eq. C_IZERO )     &
                      XP1(IC,1) = V0*DSQRT(-DLOG(Z1(IC,1)))*DCOS(CP_TWOPI*Z2(IC,1))

                   if(IAND(STAT,CP_STATU_FIXVELY) .eq. C_IZERO .and. &
                      IAND(STAT,CP_STATU_FIXPOSY) .eq. C_IZERO )     &
                      XP1(IC,2) = V0*DSQRT(-DLOG(Z1(IC,2)))*DCOS(CP_TWOPI*Z2(IC,2))

                   if(IAND(STAT,CP_STATU_FIXVELZ) .eq. C_IZERO .and. &
                      IAND(STAT,CP_STATU_FIXPOSZ) .eq. C_IZERO )     &
                      XP1(IC,3) = V0*DSQRT(-DLOG(Z1(IC,3)))*DCOS(CP_TWOPI*Z2(IC,3))
                else
                   XP1(IC,1) = 0.D0
                   XP1(IC,2) = 0.D0
                   XP1(IC,3) = 0.D0
                end if
             end if
         return
  end subroutine Thermalizing_MC_KERNEL
  !********************************************************************************

  !********************************************************************************
  subroutine Thermalizing_MC_template(IDEV, STARTCELL, ENDCELL, RGEN, CM, ITYP, STATU, XP1, TI)
  !***  PORPOSE:   template to generating and Gaussian ditribution for XP1
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                RGEN,       the handle of randnumber generator
  !                CM,         mass of atoms
  !                ITYP,       type id of atoms
  !                STATU,      statu of atoms on  device IDEV
  !                TI,         temprature to be thermalized to
  !
  !     OUTPUT     XP1,        velocities updated
  !
  use RAND32_MODULE, only:RANF=>DRAND32
  use CudaRandomC2F_M
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)              ::IDEV, STARTCELL, ENDCELL
          integer(kind=int_ptr_kind())          ::RGEN          
          real(KINDDF), device, dimension(:)    ::CM
          integer,      device, dimension(:)    ::ITYP, STATU
          real(KINDDF), device, dimension(:,:)  ::XP1
          real(KINDDF)                          ::TI

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, BX, BY, NB, ERR, CURDEV, I
           real(c_double), parameter::mean=0.D0, stddev=1.D0
           real(KINDDF), device, dimension(:,:),allocatable::Z1, Z2
             
            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            !--- first, we generate uniform random number

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT     = ENDA - STARTA + 1

            allocate(Z1(NPRT,3), Z2(NPRT,3))
            ERR = curandGenerateUniformDouble(RGEN, Z1(1:NPRT,1), NPRT)
            ERR = curandGenerateUniformDouble(RGEN, Z1(1:NPRT,2), NPRT)
            ERR = curandGenerateUniformDouble(RGEN, Z1(1:NPRT,3), NPRT)
            ERR = curandGenerateUniformDouble(RGEN, Z2(1:NPRT,1), NPRT)
            ERR = curandGenerateUniformDouble(RGEN, Z2(1:NPRT,2), NPRT)
            ERR = curandGenerateUniformDouble(RGEN, Z2(1:NPRT,3), NPRT)
            !call curandGenerateNormalDouble(RGEN, XP1(1:NPRT,1), NPRT, mean, stddev)
            !call curandGenerateNormalDouble(RGEN, XP1(1:NPRT,2), NPRT, mean, stddev)
            !call curandGenerateNormalDouble(RGEN, XP1(1:NPRT,3), NPRT, mean, stddev)
            !do I=1, NPRT
            !   Z1 = RANF()
            !   Z2 = RANF()
            !   hm_XP1(I,1)=SQRT(-LOG(Z1))*COS(CP_TWOPI*Z2)
            !end do
            !do I=1, NPRT
            !   Z1 = RANF()
            !   Z2 = RANF()
            !   hm_XP1(I,2)=SQRT(-LOG(Z1))*COS(CP_TWOPI*Z2)
            !end do
            !do I=1, NPRT
            !   Z1 = RANF()
            !   Z2 = RANF()
            !   hm_XP1(I,3)=SQRT(-LOG(Z1))*COS(CP_TWOPI*Z2)
            !end do
            !call CopyXP1From_Host_to_Devices()
            
            !$$--- to determine size of a block (the number of threads in a block)
            BX = mp_BLOCKSIZE
            BY = 1

            !$$-- to determine the dimension of blocks( the number of blocks in a grid)
            NB     = (NPRT-1)/(BX*BY)+1
            blocks  = dim3(NB, 1, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call Thermalizing_MC_KERNEL<<<blocks, threads>>>(STARTA,ITYP,m_NAPDEV, NPRT,XP1,STATU, CM, TI,Z1,Z2)

            deallocate(Z1,Z2)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine Thermalizing_MC_template
  !*********************************************************************************
  #else  
  !**********************************************************************************
  attributes(global) subroutine Thermalizing_MC_KERNEL(IA0,ITYP,NAPDEV,NPART,XP1,STATU, CM, TI, RandState)
  !***  PURPOSE:   to scal the velocitie of Gaussin ditribution of unit half width to
  !                Maxswell distribution of given temperature
  !
  !$$   INPUT:     
  !$$              IA0,        starting atom on the device
  !$$              ITYP,       indicator of type of atoms
  !$$              NAPDEV,     permitted number of atoms on the device
  !$$              NPART,      actuall number of atoms on the device
  !$$              XP1,        velocity of atoms
  !$$              STATU,      statu ot atoms
  !$$              CM,         mass of atoms
  !$$              TI,         tmeprature to thermalize to
  !
  !$$   OUTPUT     XP1,       updated velocity
  !
  implicit none
  !----   DUMMY Variables
          integer,      value::NAPDEV, NPART, IA0
          integer,      device::ITYP(*), STATU(*)
          real(KINDDF), device::XP1(NAPDEV,*),CM(*)
          real(KINDDF), value::TI
          type(curandStateXORWOW), device::RandState(*)
          

  !----   Local variables
          integer::IT, IB, IC, STAT, BS, GS, IC0, NL, L
          real(KINDDF)::V0, Z1, Z2
  !
  !***
  !     
              !--- size the block and grid
              GS  = griddim%x*griddim%y
              BS  = blockdim%x*blockdim%y
              !--- number of loops needed for cover all atoms
              NL  = (NPART-1)/(GS*BS) + 1

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  =  (blockidx%y-1)*griddim%x  +  blockidx%x

              !IC -- the id of the random number sequence 
              IC0  = (IB-1)*BS + IT

              !IS -- the rand generator ID for this thread on this GPU
              do L = 1, NL
                 !IC -- shift id of the atom for BS*GS when the total number thread cannot cover all atoms
                 IC = (L-1)*(BS*GS) + IC0
                 if(IC .le. NPART) then
                    STAT = STATU(IC)
                    if(IAND(STAT,CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                       !$$--- NOTE: the stderr of Gaussian distribution ise set to 1.d0
                       V0 = DSQRT( C_TWO*TI*CP_KB/CM(ITYP(IC+IA0)))
                       !V0 = DSQRT( TI*CP_KB/CM(ITYP(IC+IA0)))

                       if(IAND(STAT,CP_STATU_FIXVELX) .eq. C_IZERO .and. &
                          IAND(STAT,CP_STATU_FIXPOSX) .eq. C_IZERO )   then
                           Z1 = curand_uniform(RandState(IC0))
                           Z2 = curand_uniform(RandState(IC0))
                           XP1(IC,1) = V0*DSQRT(-DLOG(Z1))*DCOS(CP_TWOPI*Z2)
                       end if   

                       if(IAND(STAT,CP_STATU_FIXVELY) .eq. C_IZERO .and. &
                          IAND(STAT,CP_STATU_FIXPOSY) .eq. C_IZERO )   then
                           Z1 = curand_uniform(RandState(IC0))
                           Z2 = curand_uniform(RandState(IC0))
                           XP1(IC,2) = V0*DSQRT(-DLOG(Z1))*DCOS(CP_TWOPI*Z2)
                       end if    

                       if(IAND(STAT,CP_STATU_FIXVELZ) .eq. C_IZERO .and. &
                          IAND(STAT,CP_STATU_FIXPOSZ) .eq. C_IZERO )   then
                           Z1 = curand_uniform(RandState(IC0))
                           Z2 = curand_uniform(RandState(IC0))
                           XP1(IC,3) = V0*DSQRT(-DLOG(Z1))*DCOS(CP_TWOPI*Z2)
                       end if    
                    else
                      XP1(IC,1) = 0.D0
                      XP1(IC,2) = 0.D0
                      XP1(IC,3) = 0.D0
                    end if
                 end if
              end do     
         return
  end subroutine Thermalizing_MC_KERNEL
  !********************************************************************************

  !********************************************************************************
  subroutine Thermalizing_MC_template(IDEV, STARTCELL, ENDCELL, CM, ITYP, STATU, XP1, TI, RandState)
  !***  PORPOSE:   template to generating and Gaussian ditribution for XP1
  !
  !     INPUT:     IDEV,       device ID
  !                STARTCELL,  starting cell on the device
  !                ENDCELL,    last cell on the device
  !                RGEN,       the handle of randnumber generator
  !                CM,         mass of atoms
  !                ITYP,       type id of atoms
  !                STATU,      statu of atoms on  device IDEV
  !                TI,         temprature to be thermalized to
  !
  !     OUTPUT     XP1,        velocities updated
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)              ::IDEV, STARTCELL, ENDCELL
          real(KINDDF), device, dimension(:)    ::CM
          integer,      device, dimension(:)    ::ITYP, STATU
          real(KINDDF), device, dimension(:,:)  ::XP1
          real(KINDDF)                          ::TI
          type(DevRandState)                    ::RandState
  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::STARTA, ENDA, NPRT, ERR, CURDEV
             
            ERR = cudaGetDevice(CURDEV)
            ERR = cudaSetDevice(IDEV)

            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT     = ENDA - STARTA + 1

            blocks  = dim3(RandState%GRIDSIZE, 1, 1)
            threads = dim3(RandState%BLOCKSIZE,  1, 1)
            STARTA  = STARTA-1
            call Thermalizing_MC_KERNEL<<<blocks, threads>>>(STARTA,ITYP,m_NAPDEV, NPRT,XP1,STATU,CM, TI, RandState%RandState)

            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine Thermalizing_MC_template

  #endif
 !*********************************************************************************

 !*********************************************************************************
  subroutine Thermalizing_MC_DEV(SimBox, CtrlParam, TI)
  !***  PORPOSE:   to do global thermalizing by MC
  !
  !     INPUT:     SimBox,    the simulation box,   not used
  !                CtrlParam, the control parameters, not used
  !                TI,        the termperature (in K) to thermalize to
  !     OUTPUT     dm_XP1,    velocities of atoms on devices
  !
  !$$ use MD_SimBoxArray,     only:Thermalize_SimBoxArray
  !$$ use MD_SimBoxArray_GPU, only:CopyIn_SimBox_Vel_DEV, CopyOut_SimBox_Vel_DEV
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
          real(KINDDF)                ::TI
  !--- Local variables
          integer::I, NPRT, IB, NB, IA
          real(KINDDF)::WT, VT(3)

          !real(KINDDF)::DBINX(-100:100), DBINY(-100:100), DBINZ(-100:100), ddd, xmi, xmx
  !--- 

             do I=1,  m_NDEVICE 
               #ifdef NODEVRAN 
                call Thermalizing_MC_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_URanGen(I), dm_WorkSpace%CM(I)%Data, &
                                              dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%STATU(I)%Data, dm_WorkSpace%XP1(I)%Data, TI)
               #else 
                call Thermalizing_MC_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,              &
                                              dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%STATU(I)%Data, dm_WorkSpace%XP1(I)%Data,  &
                                              TI, gm_DevRandState(I))
               #endif                               
             end do

             !$$-- to keep velocity of the center of mass to be zero
             call CopyXP1From_Devices_to_Host(m_XP1)
             NPRT = SimBox(1)%NPRT
             NB   = size(SimBox) 
             do IB = 1, NB
                 WT = 0.d0
                 VT(1:3) = 0.D0
                 IA = (IB-1)*NPRT
                 do I=1, NPRT
                    WT = WT + m_CM(m_ITYP(IA+I))
                 end do
                 do I=1, NPRT
                     VT(1:3) = VT(1:3) +  m_CM(m_ITYP(IA+I))*m_XP1(IA+I,1:3)/WT 
                 end do
                 m_XP1(IA+1:IA+NPRT,1)  = m_XP1(IA+1:IA+NPRT,1) - VT(1)
                 m_XP1(IA+1:IA+NPRT,2)  = m_XP1(IA+1:IA+NPRT,2) - VT(2)
                 m_XP1(IA+1:IA+NPRT,3)  = m_XP1(IA+1:IA+NPRT,3) - VT(3)
             end do
             hm_XP1(1:dm_NPRT,1) = m_XP1(hm_GID(1:dm_NPRT),1)
             hm_XP1(1:dm_NPRT,2) = m_XP1(hm_GID(1:dm_NPRT),2)
             hm_XP1(1:dm_NPRT,3) = m_XP1(hm_GID(1:dm_NPRT),3)
             call CopyXP1From_Host_to_Devices()

          return
  end subroutine Thermalizing_MC_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine Do_DynDamp_Forsteps_DEV(SimBox, CtrlParam, ForceClass, MXNUMSTEPS)
   !***  PORPOSE:   to use steepest & linesearch method to relax the sytem to local minima
   !
   !     INPUT:     SimBox,     the simulation box
   !                CtrlParam,  the control parameters
   !                ForceClass, the force calculation engine
   !                MXNUMSTEPS, the MXNUMSTEPS for searching the miniam
   !
   !     OUTPUT:    dm_XP:      the configurations on devices, 
   !
    use MD_Forceclass_Register_GPU
    use MD_Globle_Variables_GPU
   implicit none
   !----   DUMMY Variables
          type(SimMDBox),dimension(:)           :: SimBox
          type(SimMDCtrl),            intent(in):: CtrlParam
          type(MDForceClassGPU),      intent(in):: ForceClass
          integer,                    intent(in):: MXNUMSTEPS
          !----  Local variables
           type(DevVec_DF), dimension(m_MXDEVICE)::dEPOT0, dEPOT1
           integer::ITER, I
           real(KINDDF)::DELEPOT, MINEPOT

             !---
             MINEPOT = CtrlParam%STEEPEST_MiDelE*CP_EV2ERG
             call DevAllocate(dEPOT0,   m_NAPDEV    )
             call DevAllocate(dEPOT1,   m_NAPDEV    )
             call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
             call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0) 

             do ITER=1, MXNUMSTEPS
                do I=1, m_NDEVICE
                   call DAMPING_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), &
                                        dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%FP(I)%Data, dm_WorkSpace%STATU(I)%Data)

                   call PREDICTOR_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,            &
                                          dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP(I)%Data, dm_WorkSpace%XP1(I)%Data,    &
                                          dm_WorkSpace%FP(I)%Data,   dm_WorkSpace%DIS(I)%Data,dm_WorkSpace%STATU(I)%Data)
                end do
      
                !--- hm_XP has bee modified in PREDICTOR_template, we should sychonize
                !    XP on devices for the next force calculations
                call Synchroniz_XP_on_Devices()
        
                !--- calculate the current force
                call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                call UpdateEPOT_ForceClass(SimBox,  CtrlParam, ForceClass)
                call DevMinus_noshift(dm_WorkSpace%EPOT, dEPOT0, dEPOT1) 
                call DevMaxAbsval_noshift(dEPOT1, DELEPOT)
                if(DELEPOT .le. MINEPOT) then
                   exit
                end if
                call DevMakeCopy_noshift(dm_WorkSpace%EPOT, dEPOT0)

                !--- to do the correction
                do I=1, m_NDEVICE
                   call CORRECTION_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), dm_WorkSpace%CM(I)%Data,         &
                                           dm_WorkSpace%ITYP(I)%Data, dm_WorkSpace%XP1(I)%Data, dm_WorkSpace%FP(I)%Data, &
                                           dm_WorkSpace%STATU(I)%Data)
                end do
             end do
             !--- clear the alloc ated memeory 
             call DevDeallocate(dEPOT0)
             call DevDeallocate(dEPOT1)

             if(ITER .gt. MXNUMSTEPS) then
               write(*,fmt="(A, I)")       " MDPSCU WARNING: DnyDamp finished after max steps:  ", ITER-1 
               write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev):", DELEPOT*CP_ERG2EV
              call ONWARNING(gm_OnWarning)
             else 
              write(*,fmt="(A, I)")       " MDPSCU Message: DnyDamp finished after steps:  ", ITER 
              write(*,fmt="(A, 1PE13.4)") "                 with max energy uncertainty(ev): ", DELEPOT*CP_ERG2EV
             end if  
          return
   end subroutine Do_DynDamp_Forsteps_DEV
   !****************************************************************************************

  end module MD_DiffScheme_GPU
