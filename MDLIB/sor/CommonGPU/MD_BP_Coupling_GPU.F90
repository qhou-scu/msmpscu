  module MD_BP_Coupling_GPU
  !**** DESCRIPTION: to calculate presure model by algorithm of:
  !                  H.J.C. Berendsen etc, J. Chem. Phys., 81(1984)3684
  !****
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_ForceTable
  use MD_MultiGPU_Basic
  implicit none

  !***********************************************************
         integer, parameter,private::mp_BLOCKSIZE=256
         integer, parameter,private::mp_BLOCKDIMX=64

         real(KINDDF), private::m_Comp          ! The compressibility
         real(KINDDF), private::m_CPTIME0       ! The charcateristic time input (in ps)
         real(KINDDF), private::m_CPINVT        ! time step H/m_CPTIME0
         real(KINDDF), private::m_MINCORR=0.001 ! min value of correction to mu factor

         private::BPC_MODIFICATION_KERNEL
         private::BPC_MODIFICATION_template

  contains

  !**************************************************************************
  subroutine Initialize_BPC_MODIFICATION_DEV(SimBox, CtrlParam)
  !***  PURPOSE: to calculate the pressure coupling with an external path
  !**** DESCRIPTION: to calulate the pressure coupling based on the algorithm by
  !                  H.J.C. Berendsen etc, J. Chem. Phys., 81(1984)3684
  !****
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box
  !
  use MD_BP_Coupling
  implicit none
      !----   DUMMY Variables
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
      !--- local variables
          character*256::STR
       character*32::STRNUMB(30)
       integer::hFile,I,N, LINE
       logical::opened

         !$$**** USE THE CPU VERSIONE HAS NO OBVIOUS EFFECT ON EFFICIENCY
         !$$     THUS WE USE CPU VERSION
         call Initialize_BPC_MODIFICATION(SimBox, CtrlParam)
         return
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if(.not.CtrlParam%IFBPC) return

            if(len_trim(CtrlParam%BPCDATA) .le. 0  ) then
               write(*,*) "MDPSCU Error: No input file for Berenson pressure coupling"
               stop
            end if
             LINE = 0
            !$$--- to load the parameters
             DO hFile = 10, 99
               INQUIRE(UNIT=hFile, OPENED=opened)
               if(.not.opened) exit
            END DO
            open(hFile, file = CtrlParam%BPCDATA, status='old')

               !$$--- get compressibility
               call GetInputStrLine(hFile,STR, LINE, "!", *100)
               call Extract_Numb(STR,1,n,STRNUMB)
               m_Comp = DRSTR(STRNUMB(1))

               !$$--- get chracteristic time
               call GetInputStrLine(hFile,STR, LINE, "!", *100)
               call Extract_Numb(STR,1,n,STRNUMB)
               m_CPTIME0 = DRSTR(STRNUMB(1))*CP_PS2S

               !$$--- get min corrector
               call GetInputStrLine(hFile,STR, LINE, "!", *100)
               call Extract_Numb(STR,1,n,STRNUMB)
               m_MINCORR = DRSTR(STRNUMB(1))

            close(hFile)

            !----

            if(m_Comp .lt. C_ZERO) then
               !$$--- NOTE: CtrlParam%H has been assumed in second
               m_CPINVT = CtrlParam%H/m_CPTIME0
            else
              m_CPINVT = m_Comp*CtrlParam%H/m_CPTIME0
            end if
         return

  100    STOP "MDPSCU Error in read parameters for B_P coupling"

   end subroutine Initialize_BPC_MODIFICATION_DEV
  !***************************************************************************

  !*******************************************************************************
  attributes(global) subroutine BPC_MODIFICATION_KERNEL(IM,IA0, NAPDEV, NPART,SX,SY,SZ,XP)
  !***  PURPOSE:   KERNEL to rescal the system by Berendson algorithm
  !
  !$$   INPUT:     IM,       the total number of atoms in the whole box
  !$$              IA0,      the starting atom on the device
  !$$              NAPDEV,    permitted number of atoms on the device
  !$$              NPART,     actuall number of atoms on the device
  !$$              SX,        the rescal factor in x
  !$$              SY,        the rescal factor in y
  !$$              SZ,        the rescal factor in z
  !$$              XP,        position of atoms
  !
  !$$   OUTPUT     XP,        updated position and velocity
  !
  implicit none
  !----   DUMMY Variables
          integer,      value, intent(in)::IM
          integer,      value, intent(in)::IA0
          integer,      value, intent(in)::NAPDEV
          integer,      value, intent(in)::NPART
          real(KINDDF), value, intent(in)::SX,SY,SZ
          real(KINDDF), device::XP(IM,*)

  !----   Local variables
          integer:: IT, IB, IC
  !

              IT  = (threadidx%y-1)*blockdim%x  + threadidx%x
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
              !$$IC -- the id of the atom
              IC= (IB-1)*(blockdim%x*blockdim%y)+IT

              if(IC.LE.NPART) then
                 !$$--- TO scale the box
                 !$$    NOTE: Boxshape is unchanged here. We assumed all the length scal
                 !$$          is in absolute coordination (cm)
                 !$$          Be careful to ref the calculations of the neighbor-list and
                 !$$          and forces
               XP(IC+IA0,1) = XP(IC,1)*SX
               XP(IC+IA0,2) = XP(IC,2)*SY
               XP(IC+IA0,3) = XP(IC,3)*SZ
             end if
        return
  end subroutine BPC_MODIFICATION_KERNEL
  !************************************************************************************************


  !****************************************************************************************
  SUBROUTINE BPC_MODIFICATION_template(IDEV, STARTCELL, ENDCELL, SX,SY, SZ, XP)
  !***  PORPOSE:   to forward one time step
  !     INPUT:     IDEV,      ID of the device
  !                STARTCELL, starting cell on the device
  !                ENDCELL,   last cell on the device
  !                SX,        the rescal factor in x
  !                SY,        the rescal factor in y
  !                SZ,        the rescal factor in z
  !                XP,        position of atoms
  !
  !     OUTPUT     XP,        updated positions
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer,      intent(in)::IDEV, STARTCELL, ENDCELL
          real(KINDDF), intent(in)::SX,SY, SZ
          real(KINDDF), device, dimension(:,:)::XP

  !--- Device variables and variables to be used in GPU
          type(dim3) :: blocks
          type(dim3) :: threads

  !--- Local variables
           integer::ERR, CURDEV, STARTA, ENDA, NPRT, BX, BY, NB, IA0

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
            IA0 = STARTA-1
            call EPC_MODIFICATION_KERNEL<<<blocks, threads>>>(dm_NPRT, IA0, m_NAPDEV,NPRT, SX,SY,SZ,XP)

            ERR = cudaMemcpyAsync(hm_XP(STARTA,1),XP(STARTA,1), NPRT)
            ERR = cudaMemcpyAsync(hm_XP(STARTA,2),XP(STARTA,2), NPRT)
            ERR = cudaMemcpyAsync(hm_XP(STARTA,3),XP(STARTA,3), NPRT)


            ERR = cudaSetDevice(CURDEV)

     return
  end subroutine BPC_MODIFICATION_template
  !****************************************************************************************

  !****************************************************************************************
  SUBROUTINE BPC_MODIFICATION_DEV(SimBox,CtrlParam)
  !***  PORPOSE: to forward one time step
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !     OUTPUT     SimBox,    the simulation box with new position and velocity
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          type(SimMDBox) ::SimBox
          type(SimMDCtrl)::CtrlParam
  !--- Device variables and variables to be used in GPU
          real(KINDDF)::MUFACTOR(3), TG
          integer::I


         !$$--- To calculate the scaling vector
           do I=1, 3
              TG = (CtrlParam%PEXTENSOR(I,I)-SimBox%PTENSOR(I,I))*m_CPINVT*C_UTH
              if(CtrlParam%PEXTENSOR(I,I) .eq. C_ZERO) then
                 TG =  DSIGN(m_MINCORR,TG)
              else
                 TG = TG/DABS(CtrlParam%PEXTENSOR(I,I))
                 if(DABS(TG) .gt. m_MINCORR ) TG = DSIGN(m_MINCORR,TG)
              end if
              MUFACTOR(I) = C_UN - TG
            end do

           do I=1, m_NDEVICE
              call BPC_MODIFICATION_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I), MUFACTOR(1),MUFACTOR(2),MUFACTOR(3),dm_WorkSpace%XP(I)%Data)
           end do
     !$$--- hm_XP has bee modified by needBPC_MODIFICATION_template
     !$$    we need to update XP on devices
           call CopyXPFrom_Host_to_Devices()
           call SynchronizeDevices()

     !$$--- need also update the boxsize
           SimBox%ZL = MUFACTOR*SimBox%ZL
           SimBox%BOXLOW = MUFACTOR*SimBox%BOXLOW
           SimBox%BOXUP =  MUFACTOR*SimBox%BOXUP

     return
  end subroutine BPC_MODIFICATION_DEV
  !****************************************************************************************




  end module MD_BP_Coupling_GPU
