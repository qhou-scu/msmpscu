  module STOP_B_MODULE
  !***  this module is to calculate the stoping power of energitic atom in a media.
  !     this module is adapted from one of my old version code used in the bipartitional
  !     model for atom transport in materials.
  !
  !     REFERNCE: Wolfgang Eckstein, Computer simulation of  Ion-Solid intercation, p67-68,
  !                                  Springer-Verlag, Berlin 1991
  !
  !
  !     HISTORY:     adopted by Hou Qing, Feb., 2018
  !
  use MD_CONSTANTS
  implicit none

    contains
!**************************************************************
     subroutine NSTOP(X,A1,A2,Z1,Z2, SN)
!*** PURPOSE: To calculate the nuclear stopping for ion of atomic charge Z1
!             and atomic mass A1 of energy X incident on matter of atomic
!             number Z2 and atomic mass A2.The formula for stopping is
!             from Biersack
!*** INPUT:   X  energy (keV)
!             A1 atomic mass of ion
!             A2 atomic mass of target atom
!             Z1 atomic number of ion
!             Z2 atomic number of target atom
!*** OUPUT:   SN the nuclear stopping in unit of (keV*cm^2)
      implicit none
      !--- dummy variables
        real(KINDDF), intent(in) ::X, A1, A2, Z1, Z2
        real(KINDDF), intent(out)::SN
      !--- local varibales
        real(KINDDF)::ESP, W0, W1, RU

        ESP=32.53*A2*X
        W0=(Z1**(2./3)+Z2**(2./3))**0.5*(A1+A2)
        ESP=ESP/(Z1*Z2*W0)
        W1=(8.462*Z1*Z2*A1)/W0*0.5*DLOG(1.+ESP)
        RU=W1/(ESP+0.10718*ESP**0.37544)

        SN=RU*1.D-18
      return
      end subroutine NSTOP

!----------------------------------------------------------------
      subroutine ESTOP(X,A1,A2,Z1,Z2,SE)
!*** PURPOSE: To calculate the electronic stopping for ion of atomic charge Z1
!             and atomic mass A1 of energy X incident on matter of atomic
!             number Z2 and atomic mass A2. The stopping formular is from
!             Bierscak
!*** INPUT:   X  energy (in keV)
!             A1 atomic mass of ion
!             A2 atomic mass of target atom
!             Z1 atomic number of ion
!             Z2 atomic number of target atom
!*** OUPUT:   SE the electronic  stopping in unit of (kev*cm^2)
      implicit none
      !--- dummy variables
        real(KINDDF), intent(in) ::X, A1, A2, Z1, Z2
        real(KINDDF), intent(out)::SE
      !--- local varibales
      !** The correction factor for hydrogen
      real(KINDDF),dimension(92), parameter::CK=(/              &
      0.93D0, 0.67D0, 0.67D0, 0.97D0, 1.01D0, 1.03D0, 1.11D0,   &
      0.97D0, 0.75D0, 0.68D0, 0.88D0, 1.29D0, 1.40D0, 1.38D0,   &
      1.06D0, 1.12D0, 1.63D0, 1.84D0, 1.64D0, 1.75D0, 1.64D0,   &
      1.52D0, 1.40D0, 1.24D0, 1.07D0, 1.08D0, 0.96D0, 1.09D0,   &
      1.13D0, 1.28D0, 1.52D0, 1.68D0, 1.60D0, 1.76D0, 1.68D0,   &
      1.92D0, 1.70D0, 1.89D0, 1.90D0, 1.99D0, 2.04D0, 1.90D0,   &
      2.00D0, 1.80D0, 1.74D0, 1.54D0, 1.65D0, 1.70D0, 1.82D0,   &
      1.87D0, 2.18D0, 2.03D0, 2.24D0, 2.39D0, 2.11D0, 2.28D0,   &
      2.32D0, 2.16D0, 2.10D0, 2.04D0, 1.99D0, 1.93D0, 1.88D0,   &
      1.93D0, 1.78D0, 1.58D0, 1.49D0, 1.45D0, 1.41D0, 1.37D0,   &
      1.40D0, 1.43D0, 1.35D0, 1.30D0, 1.48D0, 1.44D0, 1.40D0,   &
      1.27D0, 1.38D0, 1.22D0, 1.34D0, 1.51D0, 1.68D0, 1.74D0,   &
      1.75D0, 1.75D0, 1.96D0, 2.12D0, 2.16D0, 2.17D0, 2.09D0,   &
      2.05D0/)
      !----
        integer::IZ, IZ2
        real(KINDDF)::C, AI0, EB, SL, SB

        IZ=Z1+0.0000001
        IZ2=Z2+0.000001

        if(IZ.LT.3) then
           C=100.*Z1/Z2
        else
           C=5.
        end if

        if(IZ.LT.13) then
           AI0=12.+7./Z2
        else
           AI0=9.76+58.5*Z2**(-1.19)
        end if
        EB=2.19432*X/AI0/Z2
        SB=8*3.142592*Z1*Z1*1.439968**2/AI0/EB*DLOG(EB+1.+C/EB)*1.E-17
        !---Lindhard-Scharff stoping
        SL=3.8455E-18*Z1**(7./6)*Z2/(Z1**(2./3)+Z2**(2./3))**(3./2)*(X/A1)**0.5

        !--- the correction for hydrogen
        !if(IZ .eq. 1) then
           SL=SL*CK(IZ2)
        !end if

        SE=SL*SB/(SL+SB)
      return
      end subroutine ESTOP

      end module STOP_B_MODULE
