  module EMBEDMENT_Sandwich
  !***  DESCRIPTION:
  !     This program is a tool to make a sandwich box. The program will read in a box(s)
  !     that has been relaxed using SUBSTRATE, MEBEDMENT.
  !     Then two replicas of the box of the box are  made, but with embedments atoms
  !     are removed. The two repica will be place on the top and bottom of the in input
  !     box, and the combined box will output. This program will place in the framwork
  !     of MD_SimBoxArray_ToolShell_14_CPU.F90
  !
  !    Authored by HOU Qing
  !
  !
  !    SEE ALSO:      MD_SimBoxArray_ToolShell_14_CPU.F90
  use MD_CONSTANTS
  use MD_Globle_Variables
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use EMBEDMENT_TypeDef_CtrlParam_2010
  implicit none

  !_______________________________________________________________________________________
  !
  !---
       type(SimMDBox), dimension(:), allocatable,private ::m_SimBoxA
       type(SimMDBox), dimension(:), allocatable,private ::m_SimBoxB
       integer::m_RMType = 2
       real(KINDDF)::m_EXRatio = 1.8
  !_______________________________________________________________________________________
       character(len=11),parameter, private::mp_FTAGI="&AUXF_EMBED"
       type(EMBEDCtrlParam), private::m_EMBEDCtrlParam

  contains

  !****************************************************************************************
  subroutine Initialize(SimBox,CtrlParam)
  !***  DESCRIPTION: to allocate mempry for diffusion calculations
  !
  !
  implicit none
      type(SimMDBox)            ::SimBox
      type(SimMDCtrl),intent(in)::CtrlParam
      RETURN
  END subroutine Initialize
  !****************************************************************************************

  !****************************************************************************************
  subroutine Clear(SimBox,CtrlParam)
  !***  DESCRIPTION: to allocate mempry for diffusion calculations
  !
  !
  !***  The modules included ******************************************
  implicit none
      type(SimMDBox) ::SimBox
      type(SimMDCtrl)::CtrlParam
      !--- Local variables

      call Release_SimBoxArray(m_SimBoxA)
      call Release_SimBoxArray(m_SimBoxB)

      RETURN
  END subroutine Clear
  !****************************************************************************************

  !****************************************************************************************
  subroutine Creat_Sandwich_Box(JBOX,ITIME, TIME, SimBox, CtrlParam )
  !***  PORPOSE: to get the migration distance  of clusters
  !     INPUT:  SimBox, the  samples
  !     OUTPUT: RAV, RMI, RMA
  implicit none
   !--- input
   type(SimMDBox),dimension(:)           ::SimBox
   type(SimMDCtrl),            intent(in)::CtrlParam
   integer                               ::JBOX,ITIME
   real(KINDDF)                          ::TIME

   !
   !--- Local variables
   integer::NB, ID, hFile, I, J, K, TNP0, TNP
   real(KINDDF)::THICK, VC0(3)
   character*256::GFILE, ext
  !----

         if(ITIME.LE.0) return

          NB = size(SimBox)
          if(.not. allocated(m_SimBoxA)) then
              allocate(m_SimBoxA(NB))
          end if

          if(.not. allocated(m_SimBoxB)) then
             allocate(m_SimBoxB(NB))
          end if
          !--- release allocated memory if ITIME LE 0
          do I=1, NB
             call Copy_SimMDBox(SimBox(I), m_SimBoxA(I))
             call Copy_SimMDBox(SimBox(I), m_SimBoxB(I))
             m_SimBoxA(I)%XP(:,3) = m_SimBoxA(I)%XP(:,3) + m_SimBoxA(I)%ZL(3)
             m_SimBoxB(I)%XP(:,3) = m_SimBoxB(I)%XP(:,3) - m_SimBoxA(I)%ZL(3)
             !--- to make translation velocity of the combined box to be zere
             VC0 = 0.D0
             TNP = 0
             do K=1, SimBox(I)%NPRT
                if(m_SimBoxA(I)%ITYP(K) .NE.m_RMType .AND. &
                     DABS(m_SimBoxA(I)%XP(K,3)) .LE. 0.5D0*THICK ) then
                     VC0(1:3) = VC0(1:3) + m_SimBoxA(I)%XP1(K,1:3)
                     TNP = TNP + 1
                 end if

                 if(DABS(SimBox(I)%XP(K,3)) .LE. 0.5D0*THICK ) then
                     VC0(1:3) = VC0(1:3) + SimBox(I)%XP1(K,1:3)
                     TNP = TNP + 1
                 end if

                if(m_SimBoxB(I)%ITYP(K) .NE.m_RMType .AND. &
                     DABS(m_SimBoxB(I)%XP(K,3)) .LE. 0.5D0*THICK ) then
                     VC0(1:3) = VC0(1:3) + m_SimBoxA(I)%XP1(K,1:3)
                     TNP = TNP + 1
                 end if
             end do
             VC0 = VC0/dble(TNP)

             do K=1, SimBox(I)%NPRT
                m_SimBoxA(I)%XP1(K,1:3) = m_SimBoxA(I)%XP1(K,1:3)-VC0(1:3)
                SimBox(I)%XP1(K,1:3)    = SimBox(I)%XP1(K,1:3)   -VC0(1:3)
                m_SimBoxB(I)%XP1(K,1:3) = m_SimBoxB(I)%XP1(K,1:3)-VC0(1:3)
             end do

          end do


          ext =  CtrlParam%f_geometry(1:len_trim( CtrlParam%f_geometry))//"_CB"
          call STRCATI(GFILE, ext, "P", 0, 4)
          call STRCATI(GFILE, GFILE, "_", JBOX, 4)
          call GetCfgID_SimMDCtrl(CtrlParam, ITIME, ID)
          call STRCATI(GFILE, GFILE, ".", ID, 4)

          !$$--- to open the file
          call AvailableIOUnit(hFile)
          open(hFile, file = GFILE, status='unknown')
          print *, "Output combined box to ", GFILE(1:len_trim(GFILE))

          TNP0 = 0
          do I=1, NB
             THICK = m_SimBoxA(I)%ZL(3)*m_EXRatio
             TNP  = 0
             do J=1, SimBox(I)%NGROUP
                do K=1, SimBox(I)%NPRT

                  if(m_SimBoxA(I)%ITYP(K) .EQ. J .AND. J.NE.m_RMType .AND. &
                     DABS(m_SimBoxA(I)%XP(K,3)) .LE. 0.5D0*THICK ) then
                     write(hFile,100)m_SimBoxA(I)%ITYP(K), m_SimBoxA(I)%XP(K,1:3)/m_SimBoxA(I)%RR, m_SimBoxA(I)%XP1(K,1:3), m_SimBoxA(I)%STATU(K), &
                                   m_SimBoxA(I)%fP(K,1:3), -m_SimBoxA(I)%EPOT(K)*CP_ERGEV, m_SimBoxA(I)%EKIN(K)*CP_ERGEV,                    &
                                   m_SimBoxA(I)%DIS(K,1:3)/m_SimBoxA(I)%RR
                   if(I.EQ.1) TNP0 = TNP0 + 1
                   if(I.GT.1) TNP  = TNP + 1
                  end if

                  if(SimBox(I)%ITYP(K) .EQ. J .AND.  &
                     DABS(SimBox(I)%XP(K,3)) .LE. 0.5D0*THICK ) then
                     write(hFile,100)SimBox(I)%ITYP(K), SimBox(I)%XP(K,1:3)/SimBox(I)%RR, SimBox(I)%XP1(K,1:3), SimBox(I)%STATU(K), &
                                   SimBox(I)%fP(K,1:3), -SimBox(I)%EPOT(K)*CP_ERGEV, SimBox(I)%EKIN(K)*CP_ERGEV,                  &
                                   SimBox(I)%DIS(K,1:3)/SimBox(I)%RR
                      if(I.EQ.1) TNP0 = TNP0 + 1
                      if(I.GT.1) TNP  = TNP + 1
                  end if

                  if(m_SimBoxB(I)%ITYP(K) .EQ. J .AND. J.NE.m_RMType .AND. &
                     DABS(m_SimBoxB(I)%XP(K,3)) .LE. 0.5D0*THICK ) then
                     write(hFile,100)m_SimBoxB(I)%ITYP(K), m_SimBoxB(I)%XP(K,1:3)/m_SimBoxB(I)%RR, m_SimBoxB(I)%XP1(K,1:3), m_SimBoxB(I)%STATU(K), &
                                   m_SimBoxB(I)%fP(K,1:3), -m_SimBoxB(I)%EPOT(K)*CP_ERGEV, m_SimBoxB(I)%EKIN(K)*CP_ERGEV,                    &
                                   m_SimBoxB(I)%DIS(K,1:3)/m_SimBoxB(I)%RR
                     if(I.EQ.1) TNP0 = TNP0 + 1
                     if(I.GT.1) TNP  = TNP + 1
                  end if

                end do
             end do
             if(I.GT.1 .AND. TNP0 .NE. TNP) then
                write(*,*) "MDPSCU Warning: not consistent atom number", I, TNP0, TNP
                call ONWARNING(gm_OnWarning)
             end if
          end do


       close(hFile)

  100  format(I8,2X,6(1PE13.4,1X),I6,2X,15(1PE13.4,1X))
  200  format(I8,2X,I4,2X, 6(1PE13.4,1X),I6,2X,15(1PE13.4,1X))


    return
  end subroutine Creat_Sandwich_Box
  !****************************************************************************************
  end module EMBEDMENT_Sandwich

