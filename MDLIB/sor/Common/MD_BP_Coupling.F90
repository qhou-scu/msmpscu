  module MD_BP_Coupling
  !**** DESCRIPTION: to calculate presure model by algorithm of:
  !                  H.J.C. Berendsen etc, J. Chem. Phys., 81(1984)3684
  !****
  use MD_CONSTANTS
  implicit none

  !***********************************************************
         real(KINDDF), private::m_Comp          !-- The compressibility
         real(KINDDF), private::m_CPTIME0       !-- The charcateristic time input (in ps)
         real(KINDDF), private::m_CPINVT        !-- time step H/m_CPTIME0
         real(KINDDF), private::m_MINCORR=0.001 !-- min value of correction to mu factor

  contains

  !**************************************************************************
  subroutine Initialize_BPC_MODIFICATION(SimBox, CtrlParam)
  !***  PURPOSE: to calculate the pressure coupling with an external path
  !**** DESCRIPTION: to calulate the pressure coupling based on the algorithm by
  !                  H.J.C. Berendsen etc, J. Chem. Phys., 81(1984)3684
  !****
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MiniUtilities, only:GetInputStrLine, ISTR, DRSTR
  implicit none
       !---   DUMMY Variables
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

       !--- local variables
       character*256::STR,STRTMP(1)
       character*32::STRNUMB(30)
       integer::hFile,I,N, LINE
       logical::opened

         if(.not.CtrlParam%IFBPC) return

            if(len_trim(CtrlParam%BPCDATA) .le. 0  ) then
               write(*,*) "MDPSCU Error: No input file for Berenson pressure coupling"
               stop
            end if

            !$$--- to load the parameters
             DO hFile = 10, 99
               INQUIRE(UNIT=hFile, OPENED=opened)
               if(.not.opened) exit
            END DO
            open(hFile, file = CtrlParam%BPCDATA, status='old')

               LINE = 0
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
               !--- NOTE: CtrlParam%H has been assumed in second
               m_CPINVT = CtrlParam%H/m_CPTIME0
            else
              m_CPINVT = m_Comp*CtrlParam%H/m_CPTIME0
            end if
         return

  100    STOP "MDPSCU Error in read parameters for B_P coupling"

   end subroutine Initialize_BPC_MODIFICATION
   !****************************************************************************


  !**************************************************************************
  subroutine BPC_MODIFICATION(SimBox, CtrlParam)
  !***  PURPOSE: to calculate the electron-phonon coupling
  !**** DESCRIPTION: to calulate the electron-phonon coupling
  !                  based on the algorithm  by M.W.Finnis etc,
  !                  Phys. Rev., B44(1991)567
  !****
  !     INPUT: SimBox      , the simulation box
  !            CtrlParamP,   the control parameters
  !
  !     OUTPUT:SimBox,      the simulation box with force been updated
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none
  !----   DUMMY Variables
          type(SimMDBox) ::SimBox
          type(SimMDCtrl)::CtrlParam

  !----   Local variables
          integer::I
          real(KINDDF)::MUFACTOR(3), TG
  !
  !***
  !
          !$$--- To get current pressure
          call Cal_thermal_quantities_SimMDBox(SimBox)

          !$$--- To calculate the scaling vector
          DO I=1, 3
             TG = (CtrlParam%PEXTENSOR(I,I)-SimBox%PTENSOR(I,I))*m_CPINVT*C_UTH
             if(CtrlParam%PEXTENSOR(I,I) .eq. C_ZERO) then
               TG =  DSIGN(m_MINCORR,TG)
             else
               TG = TG/DABS(CtrlParam%PEXTENSOR(I,I))
               if(DABS(TG) .gt. m_MINCORR ) TG = DSIGN(m_MINCORR,TG)
             end if
             MUFACTOR(I) = C_UN - TG
          END DO

          !$$--- TO scale the box
          !$$    NOTE: Boxshape is unchanged here. We assumed all the length scal
          !$$          is in obsolute coordination (cm)
          !$$          Be careful to ref the calculations of the neighbor-list and
          !$$          and forces
          SimBox%ZL = MUFACTOR*SimBox%ZL
          SimBox%BOXLOW = MUFACTOR*SimBox%BOXLOW
          SimBox%BOXUP =  MUFACTOR*SimBox%BOXUP

          DO I=1, SimBox%NPRT
             SimBox%XP(I,1:3) = SimBox%XP(I,1:3)*MUFACTOR(1:3)
          END DO
         return

  end subroutine BPC_MODIFICATION
  !**************************************************************************

  end module MD_BP_Coupling
