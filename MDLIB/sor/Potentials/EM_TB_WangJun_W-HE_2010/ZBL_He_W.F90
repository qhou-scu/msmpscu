  module Pot_ZBL_He_W
  !***  this module is to calculate the potential and force between Au-Au atoms
  !REFERENCE:

  use MD_CONSTANTS
        implicit none


     !***  the parameter block for nuclear-nuclear part in EAM ***************
      real(KINDDF), parameter, private::IGCLONG   = 8.854187817D-12
      real(KINDDF), parameter, private::ELECHARGE  = 1.60219D-19
      real(KINDDF), parameter, private::z1=2.0D0
      real(KINDDF), parameter, private::z2=74.0D0
      real(KINDDF), parameter, private::a0=5.291772108D-11


  contains

  subroutine ZBL_NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by Born-Mayer model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR

  implicit none
      real(KINDDF)::r  !,intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      real(KINDDF)::XX,a, X, faix, dfaix


           a=0.88534D0*a0/(Z1**0.23D0+Z2**0.23D0)
           X= r*1.0D-2/a

           faix=0.18175D0*DEXP(-3.1998D0*X)+0.50986D0*DEXP(-0.94229D0*X)+0.28022D0*DEXP(-0.4029D0*X)+0.028171D0*DEXP(-0.20162D0*X)
           POTR=ELECHARGE*Z1*Z2*faix/(4.0D0*CP_PI*IGCLONG*r*1.0D-2)
           POTR=0.5D0*POTR*CP_EVERG

           dfaix=0.18175D0*DEXP(-3.1998D0*X)*(-3.1998D0)+0.50986D0*DEXP(-0.94229D0*X)*(-0.94229D0)+0.28022D0*DEXP(-0.4029D0*X)*(-0.4029D0)+0.028171D0*DEXP(-0.20162D0*X)*(-0.20162D0)
           dfaix=dfaix/a

           FPOTR=ELECHARGE*Z1*Z2/(4.0D0*CP_PI*IGCLONG)*(dfaix*r-faix*1.0D2)/r/r
           !modified accordding to NOTE-2014-12-04 in forcetable:
           !FPOTR=-CP_EVERG*FPOTR/R
           FPOTR = -CP_EVERG*FPOTR

      return
  end subroutine ZBL_NN_Interaction


  subroutine ZBL_NE_Interaction(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the attractive potential by EAM
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB



      POTB  =  0
      FPOTB  = 0

      return
  end subroutine ZBL_NE_Interaction

  end module Pot_ZBL_He_W
