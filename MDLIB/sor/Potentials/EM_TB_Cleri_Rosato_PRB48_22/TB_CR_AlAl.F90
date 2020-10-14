 module Pot_TB_CR_AlAl
 !***  this module is to calculate the potential and force between Al-Al atoms
 !REFERENCE: F.Cleri and V.Rosato, Phys.Rev.B48, (1993) 22
  use MD_CONSTANTS
  implicit none


       !***  the parameter block for nuclear-nuclear part in EAM ***************
       real(KINDDF),parameter,private::A= 0.1221D0*CP_EVERG
       real(KINDDF),parameter,private::ETA= 1.316D0*CP_EVERG
       real(KINDDF),parameter,private::P= 8.612D0
       real(KINDDF),parameter,private::Q= 2.516D0
       real(KINDDF),parameter,private::RIJ0= 4.050*CP_HALFRSQ2*CP_A2CM

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
         !FPOTR = C_TWO*POTR*PSR0/R
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
     !FPOTB  = POTB * QSR0/R
     FPOTB  = POTB *QSR0*C_TWO

     return
 end subroutine TB_NE_Interaction

 end module Pot_TB_CR_AlAl
