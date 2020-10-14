  module MD_EP_Coupling
  !**** DESCRIPTION: to calculate EPC model by algorithm of:
  !                  M.W.Finnis etc, Phys. Rev., B44(1991)567
  !****
  !     HISTORY:    2018-03-12(HOU Qing):
  !                           delete the Initialize_EPC_MODIFICATION in the old version
  !                           the PEC physical parameters is passed from CtrlParam
  !                           which defined the parameters for atom groups

  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none

  !***********************************************************
  contains

  !**************************************************************************
  subroutine EPC_MODIFICATION(SimBox, CtrlParam)
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
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam

      !local variables
      integer::I,K,N
      real(KINDDF)::C, XP1SQ, TM, MU, TE, TECUT,ALPHAM, EMAX
      integer::first = 0


        !$$TECUT = 0.1D0 !0.25*TE    !The cutoff temperature
        N = 0
        do K=1,SimBox%NGROUP
           TE     = CtrlParam%LT_CTRL(K)%TI
           C      = SimBox%CM(K)*C_UTH/CP_KB
           !ALPHAM = SimBox%EPA(K)
           !TECUT  =  SimBox%EPACUT(K)*TE
           !SimBox%EPA = EMM * 10.D0**CVEMM/ALPHA
           ALPHAM =  SimBox%CM(K)/CtrlParam%LT_CTRL(K)%EPC_Alpha
           TECUT  =  TE*CtrlParam%LT_CTRL(K)%EPC_CUT
           EMAX   =  2.D0*CtrlParam%LT_CTRL(K)%EPC_HE
           do I=1, SimBox%NA(K)
              N = N + 1
              XP1SQ = sum(SimBox%XP1(N,1:3)*SimBox%XP1(N,1:3))
              if(XP1SQ .le. EMAX) then
                 TM = XP1SQ*C
                 if(TM.lt.TECUT) then
                    MU = ALPHAM*(TM-TE)/TECUT
                 else
                    MU = ALPHAM*(TM-TE)/TM
                 end if
                 !MU = MU*1.D3   ! convert to CGS
                 SimBox%FP(N,1:3) = SimBox%FP(N,1:3) - MU*SimBox%XP1(N,1:3)
              end if
           end do
        end do

        return
  end subroutine EPC_MODIFICATION
  !**************************************************************************

  end module MD_EP_Coupling
