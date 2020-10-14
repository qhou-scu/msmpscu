module Pot_TB_CR_TiTi
!***  this module is to calculate the potential and force between Ti-Ti atoms
!REFERENCE: F.Cleri and V.Rosato, Phys.Rev.B48, (1993) 22

use MD_CONSTANTS
      implicit none


      !***  the parameter block for nuclear-nuclear part in EAM ***************
      real(KINDDF),parameter,private::A= 0.1519*CP_EVERG
      real(KINDDF),parameter,private::ETA= 1.8112*CP_EVERG
      real(KINDDF),parameter,private::P= 8.620
      real(KINDDF),parameter,private::Q= 2.390

      real(KINDDF),parameter,private::BETA = 1.5874
      real(KINDDF),parameter,private::RIJ0 = 2.951D-8*(1./3.+BETA*BETA/4.)**0.5

contains

subroutine TB_NN_Interaction(r,POTR, FPOTR)
!***  PURPOSE: to calculate the potential and force by Born-Mayer model
!     INPUT:   r, the distance between two particle
!     OUTPUT:  POTR, FPOTR

implicit none
    real(KINDDF),intent(in)::r
    real(KINDDF),intent(out)::POTR, FPOTR
    !******
    real(KINDDF), parameter::PSR0 = P/RIJ0

        POTR  = A * DEXP(-P*((R/RIJ0)-C_UN))
        !modified accordding to NOTE-2014-12-04 in forcetable
        !FPOTR = TWO*POTR*PSR0/R
        FPOTR = C_TWO*POTR*PSR0

    return
end subroutine TB_NN_Interaction


subroutine TB_NE_Interaction(r,POTB, FPOTB)
!***  PURPOSE: to calculate the attractive potential by EAM
!
!     INPUT:   r, the distance between two particle
!     OUTPUT:
implicit none
    real(KINDDF),intent(in)::r
    real(KINDDF),intent(out)::POTB, FPOTB
    real(KINDDF), parameter::QSR0 = Q/RIJ0
    real(KINDDF), parameter::ETA2 = ETA*ETA
    real(KINDDF), parameter::TWOQ = C_TWO*Q


    POTB  =  ETA2 * DEXP(-TWOQ*((R/RIJ0)-C_UN))
   !modified accordding to NOTE-2014-12-04 in forcetable
   ! FPOTB  = POTB * QSR0/R
    FPOTB  = POTB * QSR0*C_TWO

    return
end subroutine TB_NE_Interaction

end module Pot_TB_CR_TiTi
