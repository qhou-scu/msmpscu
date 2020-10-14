  module Pot_LJ_HeHe
  !***  this module is to calculate the potential and force between He-He atoms
  !REFERENCE: Wang Jun,

  use MD_CONSTANTS
        implicit none


        !***  the parameter block for pairwise interaction of Lennard-Jones type
        real(KINDDF),parameter,private::pha0  = 0.000876*CP_EVERG
        real(KINDDF),parameter,private::r0    = 2.559D-8

  contains

  subroutine NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by LJ model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !******
      real(KINDDF)::FX,FX1,redx

       redx=r0/r

       POTR=0.5*pha0*(redx**12-2.D0*redx**6)
      !modified accordding to NOTE-2014-12-04 in forcetable
      ! FPOTR = 0.5*C_TWO*pha0*(12.D0*redx**12/r-12.D0*redx**6/r)/r
       FPOTR = 0.5*C_TWO*pha0*(12.D0*redx**12/r-12.D0*redx**6/r)

      return
  end subroutine NN_Interaction


  subroutine NE_Interaction(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the attractive potential by EAM
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB


      POTB   =  0.D0
      FPOTB  =  0.D0

      return
  end subroutine NE_Interaction

  end module Pot_LJ_HeHe
