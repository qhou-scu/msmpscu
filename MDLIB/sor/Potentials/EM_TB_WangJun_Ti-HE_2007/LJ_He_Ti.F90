module Pot_LJ_He_Ti
!***  this module is to calculate the potential and force between Au-Au atoms
!REFERENCE: F.Cleri and V.Rosato, Phys.Rev.B48, (1993) 22

use MD_CONSTANTS
      implicit none


      !***  the parameter block for pairwise interaction of Lennard-Jones type
      real(KINDDF),parameter,private::pha0   = 0.02617D0*CP_EVERG
      real(KINDDF),parameter,private::r0     = 4.72879D-8
      !real(KINDDF),parameter,private::cons  = 0.07199*EVERG

contains

!******************************************************************

!******************************************************************
subroutine NN_HeTi_Interaction(r,POTR, FPOTR)
!***  PURPOSE: to calculate the potential and force by Born-Mayer model
!     INPUT:   r, the distance between two particle
!     OUTPUT:  POTR, FPOTR

implicit none
    real(KINDDF),intent(in)::r
    real(KINDDF),intent(out)::POTR, FPOTR
    !******

    real(KINDDF)::redx

     redx=r0/r

     POTR=0.5*pha0*(redx**4.35-3.211829D0*redx**1.4)  !+cons
    !modified accordding to NOTE-2014-12-04 in forcetable:
     !FPOTR =0.5*C_TWO*pha0*(4.35D0*redx**4.35/r-3.211829D0*1.4D0*redx**1.4/r)/r
     FPOTR =0.5*C_TWO*pha0*(4.35D0*redx**4.35/r-3.211829D0*1.4D0*redx**1.4/r)

    return
end subroutine NN_HeTi_Interaction
!******************************************************************

!******************************************************************
subroutine NE_HeTi_Interaction(r,POTB, FPOTB)
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
end subroutine NE_HeTi_Interaction
!**********************************************************************************

!**********************************************************************************
  subroutine ZBL_HeHe_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by ZBL model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_CorePotent
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !******
       real(KINDDF), parameter::HeZ = 2.D0

       call ZBL_Pot(HeZ, HeZ, r,POTR, FPOTR)

      return
  end subroutine ZBL_HeHe_Interaction
!**********************************************************************************

!**********************************************************************************
  subroutine EXP6_HeHe_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by EXP6 model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !******
      real(KINDDF)::FX,FX1,redx
      real(KINDDF),parameter::CE = 0.000930644D0*CP_EVERG, &
                              R0 = 0.29673D-7,             &
                              alpha=13.1D0

       redx=R0/r

       POTR=6.0D0/(alpha-6.0D0)*DEXP(alpha*(1.0D0-1.0D0/redx))-(alpha/(alpha-6.0D0))*redx**6.0D0
       POTR=0.5*CE*POTR

       FPOTR =6.0D0/(alpha-6.0D0)*DEXP(alpha*(1.0D0-1.0D0/redx))*alpha/R0-    &
              (alpha/(alpha-6.0D0))*6.0D0*redx**5.0D0*R0/r/r
       !modified accordding to NOTE-2014-12-04 in forcetable:
       !FPOTR =CE*FPOTR/r
       FPOTR = CE*FPOTR

      return
  end subroutine EXP6_HeHe_Interaction
  !**********************************************************************************

  !**********************************************************************************
  subroutine NN_HeHe_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by EXP6 model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !******
      real(KINDDF),parameter::m_RZBL=1.D-8, m_REXP6=1.3D-8
      real(KINDDF),parameter::m_SC(4) = (/3.5443059983603055D-011, &
                                       -4.8619762395294636D-003, &
                                   200609.6243721997,            &
                           -2069784026620.333/)


       if(r < m_RZBL) then
          call ZBL_HeHe_Interaction(R,POTR, FPOTR)
       else if (r < m_REXP6) then
            POTR =   (m_SC(1) + m_SC(2)*R + m_SC(3)*R*R + m_SC(4)*R*R*R)*C_HALF
           !modified accordding to NOTE-2014-12-04 in forcetable:
           ! FPOTR = -(m_SC(2) + 2.D0*m_SC(3)*R + 3.D0*m_SC(4)*R*R)/R
            FPOTR = -(m_SC(2) + 2.D0*m_SC(3)*R + 3.D0*m_SC(4)*R*R)
       else
          call EXP6_HeHe_Interaction(R,POTR, FPOTR)
       end if

      return
  end subroutine NN_HeHe_Interaction
  !**********************************************************************************

  !**********************************************************************************
  subroutine NE_HeHe_Interaction(r,POTB, FPOTB)
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
  end subroutine NE_HeHe_Interaction
  !**********************************************************************************


end module Pot_LJ_He_Ti
