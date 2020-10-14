  module MD_SimboxArray_GPU
  !***  DESCRIPTION: this module provides a number of interfaces data transfer between
  !                  Simulation boxes and the module MD_Globle_Variables_GPU
  !
  !                  NOTE: this modeul must be called after calling Initialize_Globle_Variables_DEV
  !                        which is defined in MD_Globle_Variables_GPU.F90
  !
  !                  SEE ALSO: MD_Globle_Variables_GPU.F90.cuf
  !
  ! **** HOSTORY:
  !       1. May, 2012(HOU Qing):    Original version
  !
  !       2. June, 2018(HOU Qing):   the data transfer are made to directly copyin/out
  !                                  between SimMDBox and the partitioned system
  !

  use MD_Globle_Variables_GPU
  implicit none

  !---------------------------------------------------------
          private:: CopyIn_SimBoxB
          private:: CopyIn_SimBoxA
          public :: CopyIn_SimBox_DEV
          interface CopyIn_SimBox_DEV
             module procedure CopyIn_SimBoxB
             module procedure CopyIn_SimBoxA
          end interface CopyIn_SimBox_DEV

  !---------------------------------------------------------
          private:: CopyIn_SimBoxB_XP1
          private:: CopyIn_SimBoxA_XP1
          public :: CopyIn_SimBox_Vel_DEV
          interface CopyIn_SimBox_Vel_DEV
             module procedure CopyIn_SimBoxB_XP1
             module procedure CopyIn_SimBoxA_XP1
          end interface CopyIn_SimBox_Vel_DEV

  !---------------------------------------------------------
          private:: CopyOut_SimBoxB
          private:: CopyOut_SimBoxA
          private:: CopyOut_SimBoxA1
          interface CopyOut_SimBox_DEV
             module procedure CopyOut_SimBoxB
             module procedure CopyOut_SimBoxA
             module procedure CopyOut_SimBoxA1
          end interface CopyOut_SimBox_DEV

  !---------------------------------------------------------
          private:: CopyOut_SimBoxB_DIS
          private:: CopyOut_SimBoxA_DIS
          interface CopyOut_SimBox_Dis_DEV
             module procedure CopyOut_SimBoxB_DIS
             module procedure CopyOut_SimBoxA_DIS
          end interface CopyOut_SimBox_Dis_DEV

  !---------------------------------------------------------
          private:: CopyOut_SimBoxB_EPOT
          private:: CopyOut_SimBoxA_EPOT
          interface CopyOut_SimBox_EPOT_DEV
             module procedure CopyOut_SimBoxB_EPOT
             module procedure CopyOut_SimBoxA_EPOT
          end interface CopyOut_SimBox_EPOT_DEV

  !---------------------------------------------------------
          private:: CopyOut_SimBoxB_FP
          private:: CopyOut_SimBoxA_FP
          interface CopyOut_SimBox_For_DEV
             module procedure CopyOut_SimBoxB_FP
             module procedure CopyOut_SimBoxA_FP
          end interface CopyOut_SimBox_For_DEV

  !---------------------------------------------------------
          private:: CopyOut_SimBoxB_XP0
          private:: CopyOut_SimBoxA_XP0
          interface CopyOut_SimBox_Pos_DEV
             module procedure CopyOut_SimBoxB_XP0
             module procedure CopyOut_SimBoxA_XP0
          end interface CopyOut_SimBox_Pos_DEV
  !---------------------------------------------------------
          private:: CopyOut_SimBoxB_XP1
          private:: CopyOut_SimBoxA_XP1
          interface CopyOut_SimBox_Vel_DEV
             module procedure CopyOut_SimBoxB_XP1
             module procedure CopyOut_SimBoxA_XP1
          end interface CopyOut_SimBox_Vel_DEV

  !---------------------------------------------------------
  contains

  !****************************************************************************
  subroutine CopyIn_SimBoxA( SimBox)
  !***  PURPOSE:   to copyin the host simulation boxes to the device memories on
  !                the major device
  !
  !     INPUT     SimBox,    the array of the simulation boxes
  !
      implicit none
      !--- dummy variables
      type(SimMDBox),dimension(:)::SimBox
      !--- Local vairables
      integer::ERR, NPRT, TNP,I, IS, NSB !, CURDEV

         !$$--- prepare copy the simulation box to device memory
             NPRT = SimBox(1)%NPRT
             TNP = NPRT*size(SimBox)
             if(dm_NPRT .ne. TNP) then
                write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyIn_SimBoxA", dm_NPRT, TNP
                stop
             end if

         !$$--- first copy the boxes to a temperaroty arraies
              NSB = size(SimBox)
              IS = 0
              do I=1, NSB
                 m_ITYP(1+IS:NPRT+IS)    = SimBox(I)%ITYP(1:NPRT)
                 IS = IS+NPRT
              end do

              IS = 0
              do I=1, NSB
                 m_XP(1+IS:NPRT+IS,1:3)  = SimBox(I)%XP(1:NPRT,1:3)
                 IS = IS+NPRT
              end do

              IS = 0
              do I=1, NSB
                 m_XP1(1+IS:NPRT+IS,1:3)  = SimBox(I)%XP1(1:NPRT,1:3)
                 m_DIS(1+IS:NPRT+IS,1:3)  = SimBox(I)%DIS(1:NPRT,1:3)
                 m_FP(1+IS:NPRT+IS,1:3)   = SimBox(I)%FP(1:NPRT,1:3)
                 IS = IS+NPRT
              end do

              IS = 0
              do I=1, NSB
                 m_STATU(1+IS:NPRT+IS)  = SimBox(I)%STATU(1:NPRT)
                 IS = IS+NPRT
              end do

              if(allocated(m_XP2)) then
                 IS = 0
                 do I=1, NSB
                    m_XP2(1+IS:NPRT+IS,1:3) = SimBox(I)%XP2(1:NPRT,1:3)
                    m_XP3(1+IS:NPRT+IS,1:3) = SimBox(I)%XP3(1:NPRT,1:3)
                    m_XP4(1+IS:NPRT+IS,1:3) = SimBox(I)%XP4(1:NPRT,1:3)
                    m_XP5(1+IS:NPRT+IS,1:3) = SimBox(I)%XP5(1:NPRT,1:3)
                    IS = IS+NPRT
                 end do
              end if

              !$$--- NOTE: because orignal positions have been changed.
              !$$          the system needed to be partiioned
              hm_HasPart = mp_NOTPART

              return
  end subroutine CopyIn_SimBoxA
  !****************************************************************************

  !****************************************************************************
  subroutine CopyIn_SimBoxB( SimBox)
  !***  PURPOSE:   to copyin the host simulation boxes to the device memories on
  !                the major device
  !
  !     INPUT     SimBox,    the array of the simulation boxes
  !
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::ERR, NPRT, I !, CURDEV


         !$$--- prepare copy the simulation box to device memory
             NPRT    = SimBox%NPRT
             if(dm_NPRT .ne. NPRT) then
                write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyIn_SimBoxB", dm_NPRT, NPRT
                stop
             end if
         !$$--- first copy the boxes to a temperaroty arraies
             m_ITYP(1:NPRT)     = SimBox%ITYP(1:NPRT)
             m_XP  (1:NPRT,1:3) = SimBox%XP(1:NPRT,1:3)
             m_XP1(1:NPRT,1:3)  = SimBox%XP1(1:NPRT,1:3)
             m_DIS(1:NPRT,1:3)  = SimBox%DIS(1:NPRT,1:3)
             m_FP(1:NPRT,1:3)   = SimBox%FP(1:NPRT,1:3)
             m_STATU(1:NPRT)    = SimBox%STATU(1:NPRT)

             if(allocated(m_XP2)) then
                m_XP2(1:NPRT,1:3) = SimBox%XP2(1:NPRT,1:3)
                m_XP3(1:NPRT,1:3) = SimBox%XP3(1:NPRT,1:3)
                m_XP4(1:NPRT,1:3) = SimBox%XP4(1:NPRT,1:3)
                m_XP5(1:NPRT,1:3) = SimBox%XP5(1:NPRT,1:3)
              end if

              !$$--- NOTE: because orignal positions have been changed.
              !$$          the system needed to be partiioned
              hm_HasPart = mp_NOTPART

              return
  end subroutine CopyIn_SimBoxB
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxA(SimBox)
  !***  PURPOSE:   to copyout the global device memories and
  !                reset the simulation box on host
  !
  !     INPUT     SimBox,    the Simulation box
  !
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_DEV", dm_NPRT, NPRT*NB
                 stop
              end if

              call CopyAllFrom_Devices_to_Host()
              call SynchronizeDevices()

              IS = 0
              do I=1, NB
                 SimBox(I)%ITYP(1:NPRT)    = hm_ITYP(hm_GIDINV(IS+1:IS+NPRT))
                 SimBox(I)%XP(1:NPRT,1:3)  = hm_XP(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%XP1(1:NPRT,1:3) = hm_XP1(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%DIS(1:NPRT,1:3) = hm_DIS(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%STATU(1:NPRT)   = hm_STATU(hm_GIDINV(IS+1:IS+NPRT))
                 SimBox(I)%FP(1:NPRT,1:3)  = hm_FP(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%EPOT(1:NPRT)    = hm_EPOT(hm_GIDINV(IS+1:IS+NPRT))
                 SimBox(I)%EKIN(1:NPRT)    = hm_EKIN(hm_GIDINV(IS+1:IS+NPRT))
                 IS = IS + NPRT
              end do
            return
  end subroutine CopyOut_SimBoxA
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxA1(SimBox, IB0)
  !***  PURPOSE:   to copyout the global device memories and
  !                reset the simulation box on host
  !
  !     INPUT     SimBox,    the Simulation box
  !
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)           ::SimBox
      integer,                    intent(in) ::IB0                
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if((IB0-1+NB)*NPRT .gt. dm_NPRT) then
                write(*,*) "MDPSCU Error: number of atoms on device smaller than that to be copyout in CopyOut_SimBoxA1", dm_NPRT, (IB0-1+NB)*NPRT
                stop
              end if
              call CopyAllFrom_Devices_to_Host()
              call SynchronizeDevices()

              IS = (IB0 - 1)*NPRT
              do I=1, NB
                 SimBox(I)%ITYP(1:NPRT)    = hm_ITYP(hm_GIDINV(IS+1:IS+NPRT))
                 SimBox(I)%XP(1:NPRT,1:3)  = hm_XP(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%XP1(1:NPRT,1:3) = hm_XP1(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%DIS(1:NPRT,1:3) = hm_DIS(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%STATU(1:NPRT)   = hm_STATU(hm_GIDINV(IS+1:IS+NPRT))
                 SimBox(I)%FP(1:NPRT,1:3)  = hm_FP(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 SimBox(I)%EPOT(1:NPRT)    = hm_EPOT(hm_GIDINV(IS+1:IS+NPRT))
                 SimBox(I)%EKIN(1:NPRT)    = hm_EKIN(hm_GIDINV(IS+1:IS+NPRT))
                 IS = IS + NPRT
              end do
            return
  end subroutine CopyOut_SimBoxA1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxB(SimBox)
  !***  PURPOSE:   to copyout the global device memories and
  !                reset the simulation box on host
  !
  !     INPUT     SimBox,    the Simulation box
  !
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::err, NPRT

         !$$--- copy the simulation box to device memory
              NPRT = SimBox%NPRT
              if(dm_NPRT .ne. NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_DEV", dm_NPRT, NPRT
                 stop
              end if

              call CopyAllFrom_Devices_to_Host()
              call SynchronizeDevices()

              SimBox%ITYP(1:NPRT)    = hm_ITYP(hm_GIDINV(1:NPRT))
              SimBox%XP(1:NPRT,1:3)  = hm_XP(hm_GIDINV(1:NPRT),1:3)
              SimBox%XP1(1:NPRT,1:3) = hm_XP1(hm_GIDINV(1:NPRT),1:3)
              SimBox%DIS(1:NPRT,1:3) = hm_DIS(hm_GIDINV(1:NPRT),1:3)
              SimBox%STATU(1:NPRT)   = hm_STATU(hm_GIDINV(1:NPRT))
              SimBox%FP(1:NPRT,1:3)  = hm_FP(hm_GIDINV(1:NPRT),1:3)
              SimBox%EPOT(1:NPRT)    = hm_EPOT(hm_GIDINV(1:NPRT))
              SimBox%EKIN(1:NPRT)    = hm_EKIN(hm_GIDINV(1:NPRT))

            return
  end subroutine CopyOut_SimBoxB
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxA_XP0( SimBox)
  !***  PURPOSE: to copyout the position array in device memory to
  !              to the position of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_X0_DEV", dm_NPRT, NPRT*NB
                 stop
              end if

              call CopyXPFrom_Devices_to_Host()
              call SynchronizeDevices()
              IS = 0
              do I=1, NB
                 SimBox(I)%XP(1:NPRT,1:3)  = hm_XP(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 IS = IS + NPRT
              end do
              return
  end subroutine CopyOut_SimBoxA_XP0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxB_XP0( SimBox)
  !***  PURPOSE: to copyout the position array in device memory to
  !              to the position of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::err, NPRT


              NPRT = SimBox%NPRT
              if(dm_NPRT .ne. NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_X0_DEV", dm_NPRT, NPRT
                 stop
              end if
              call CopyXPFrom_Devices_to_Host()
              call SynchronizeDevices()
              SimBox%XP(1:NPRT,1:3) = hm_XP(hm_GIDINV(1:NPRT),1:3)
              return
  end subroutine CopyOut_SimBoxB_XP0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyIn_SimBoxA_XP1( SimBox)
  !***  PURPOSE:   to copyin the velocity arraies of simulation boxes on host to
  !                to the device memory storing velocity.
  !                this subroutine could be used when thermalization is performed on host.
  !                the usage of this subroutine can be found in MD_SimboxArray_AppShell_12_GPU
  !                and MD_SimboxArray_AppShell_14_GPU

  !     INPUT     SimBox,    the simulation box
  !     OUTPUT    hm_XP1,    the velocities of partioned system on host
  !               dxm_XP1,   the velocity on devices
  !
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::NPRT, NB, I, IS

         !$$---
              if(hm_HasPart .lt. mp_HASPART) then
                 write(*,*) "MDPSCU Error: system has not been partitioned. Can not copy XP1 from original boxes to device"
                 stop
              end if

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyIn_SimBox_XP1", dm_NPRT, NPRT*NB
                 stop
              end if

              IS = 0
              do I=1, NB
                 !m_XP1(1+IS:NPRT+IS,1:3)  = SimBox(I)%XP1(1:NPRT,1:3)
                 hm_XP1(hm_GIDINV(IS+1:IS+NPRT), 1:3) = SimBox(I)%XP1(1:NPRT,1:3)
                 IS = IS+NPRT
              end do
              call CopyXP1From_Host_to_Devices()
              return
  end subroutine CopyIn_SimBoxA_XP1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyIn_SimBoxB_XP1( SimBox)
  !***  PURPOSE:   to copyin the velocity arraies of simulation boxes on host to
  !                to the device memory storing velocity.
  !                this subroutine could be used when thermalization is performed on host.
  !                the usage of this subroutine can be found in MD_SimboxArray_AppShell_12_GPU
  !                and MD_SimboxArray_AppShell_14_GPU

  !     INPUT     SimBox,    the simulation box
  !     OUTPUT    hm_XP1,    the velocities of partioned system on host
  !               dxm_XP1,   the velocity on devices
  !
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::NPRT

         !$$---
              if(hm_HasPart .lt. mp_HASPART) then
                 write(*,*) "MDPSCU Error: system has not been partitioned. Can not copy XP1 from original boxes to device"
                 stop
              end if

         !$$--- copy the simulation box to device memory
              NPRT = SimBox%NPRT
              if(NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyIn_SimBox_XP1", dm_NPRT, NPRT
                 stop
              end if

              !m_XP1(1+IS:NPRT+IS,1:3)  = SimBox(I)%XP1(1:NPRT,1:3)
              hm_XP1(hm_GIDINV(1:NPRT), 1:3) = SimBox%XP1(1:NPRT,1:3)
              call CopyXP1From_Host_to_Devices()
              return
  end subroutine CopyIn_SimBoxB_XP1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxA_XP1( SimBox)
  !***  PURPOSE: to copyout the velocity array in device memory to
  !              to the velocity of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_Vel_DEV", dm_NPRT, NPRT*NB
                 stop
              end if

              call CopyXP1From_Devices_to_Host()
              call SynchronizeDevices()
              IS = 0
              do I=1, NB
                 SimBox(I)%XP1(1:NPRT,1:3)  = hm_XP1(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 IS = IS + NPRT
              end do
              return
  end subroutine CopyOut_SimBoxA_XP1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxB_XP1( SimBox)
  !***  PURPOSE: to copyout the velocity array in device memory to
  !              to the velocity of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::err, NPRT


              NPRT = SimBox%NPRT
              if(dm_NPRT .ne. NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_Vel_DEV", dm_NPRT, NPRT
                 stop
              end if

              call CopyXP1From_Devices_to_Host()
              call SynchronizeDevices()
              SimBox%XP1(1:NPRT,1:3) = hm_XP1(hm_GIDINV(1:NPRT),1:3)
              return
  end subroutine CopyOut_SimBoxB_XP1
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxA_FP( SimBox)
  !***  PURPOSE: to copyout the force array in device memory to
  !              to the force of simulation boxes on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_For_DEV", dm_NPRT, NPRT*NB
                 stop
              end if

              call CopyFPFrom_Devices_to_Host()
              call SynchronizeDevices()
              IS = 0
              do I=1, NB
                 SimBox(I)%FP(1:NPRT,1:3)  = hm_FP(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 IS = IS + NPRT
              end do
              return
  end subroutine CopyOut_SimBoxA_FP
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxB_FP( SimBox)
  !***  PURPOSE: to copyout the force array in device memory to
  !              to the force of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::err, NPRT


              NPRT = SimBox%NPRT
              if(dm_NPRT .ne. NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_For_DEV", dm_NPRT, NPRT
                 stop
              end if

              call CopyFPFrom_Devices_to_Host()
              call SynchronizeDevices()
              SimBox%FP(1:NPRT,1:3) = hm_FP(hm_GIDINV(1:NPRT),1:3)
              return
  end subroutine CopyOut_SimBoxB_FP
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxA_EPOT( SimBox)
  !***  PURPOSE: to copyout the EPOT array in device memory to
  !              to the EPOT of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_Epot_DEV", dm_NPRT, NPRT*NB
                 stop
              end if

              call CopyEPOTFrom_Devices_to_Host()
              call SynchronizeDevices()
              IS = 0
              do I=1, NB
                 SimBox(I)%EPOT(1:NPRT)  = hm_EPOT(hm_GIDINV(IS+1:IS+NPRT))
                 IS = IS + NPRT
              end do
              return
  end subroutine CopyOut_SimBoxA_EPOT
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxB_EPOT( SimBox)
  !***  PURPOSE: to copyout the epot array in device memory to
  !              to the EPOT of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::err, NPRT


              NPRT = SimBox%NPRT
              if(dm_NPRT .ne. NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_For_DEV", dm_NPRT, NPRT
                 stop
              end if

              call CopyEPOTFrom_Devices_to_Host()
              call SynchronizeDevices()
              SimBox%EPOT(1:NPRT) = hm_EPOT(hm_GIDINV(1:NPRT))
              return
  end subroutine CopyOut_SimBoxB_EPOT
  !****************************************************************************


  !****************************************************************************
  subroutine CopyOut_SimBoxA_DIS( SimBox)
  !***  PURPOSE: to copyout the displacement array in device memory to
  !              to the displacement arraies simulation boxes on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox), dimension(:)::SimBox
      !--- Local vairables
      integer::err, NPRT, NB, I, IS

         !$$--- copy the simulation box to device memory
              NPRT = SimBox(1)%NPRT
              NB   = size(SimBox)
              if(NB*NPRT .ne. dm_NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_Vel_DEV", dm_NPRT, NPRT*NB
                 stop
              end if

              call CopyDISFrom_Devices_to_Host()
              call SynchronizeDevices()
              IS = 0
              do I=1, NB
                 SimBox(I)%DIS(1:NPRT,1:3)  = hm_DIS(hm_GIDINV(IS+1:IS+NPRT),1:3)
                 IS = IS + NPRT
              end do
              return
  end subroutine CopyOut_SimBoxA_DIS
  !****************************************************************************

  !****************************************************************************
  subroutine CopyOut_SimBoxB_DIS( SimBox)
  !***  PURPOSE: to copyout the displacement array in device memory to
  !              to the displacement of a simulation box on the  host
  !
  !     INPUT
  !     OUTPUT   SimBox, the Simulation box with position array has
  !                      been updated
      implicit none
      !--- dummy variables
      type(SimMDBox)::SimBox
      !--- Local vairables
      integer::err, NPRT


              NPRT = SimBox%NPRT
              if(dm_NPRT .ne. NPRT) then
                 write(*,*) "MDPSCU Error: simulation box in host is difference from that in device in CopyOut_SimBox_Vel_DEV", dm_NPRT, NPRT
                 stop
              end if

              call CopyDISFrom_Devices_to_Host()
              call SynchronizeDevices()
              SimBox%DIS(1:NPRT,1:3) = hm_DIS(hm_GIDINV(1:NPRT),1:3)
              return
  end subroutine CopyOut_SimBoxB_DIS
  !****************************************************************************

  end module MD_SimboxArray_GPU

