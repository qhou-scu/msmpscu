  module Pot_ZBL_EXP6_HeHe
  !***  this module is to calculate the potential and force between He-He atoms
  !     in EXP6 form
  !REFERENCE: Wang Jun,

  use MD_CONSTANTS
  implicit none
  real(KINDDF),private::m_RZBL=1.D-8, m_REXP6=1.3D-8
  real(KINDDF),private::m_SC(4) = (/3.5443059983603055D-011, &
                                   -4.8619762395294636D-003, &
                                    200609.6243721997,       &
                                   -2069784026620.333/)

  contains

 !**********************************************************************************
  subroutine ZBL_Interaction(r,POTR, FPOTR)
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
  end subroutine ZBL_Interaction
  !**********************************************************************************

  !**********************************************************************************
  subroutine EXP6_Interaction(r,POTR, FPOTR)
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
  end subroutine EXP6_Interaction
  !**********************************************************************************

  !**********************************************************************************
  subroutine NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by EXP6 model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !******

       if(r < m_RZBL) then
          call ZBL_Interaction(R,POTR, FPOTR)
       else if (r < m_REXP6) then
            POTR =   (m_SC(1) + m_SC(2)*R + m_SC(3)*R*R + m_SC(4)*R*R*R)*C_HALF
           !modified accordding to NOTE-2014-12-04 in forcetable:
           ! FPOTR = -(m_SC(2) + 2.D0*m_SC(3)*R + 3.D0*m_SC(4)*R*R)/R
            FPOTR = -(m_SC(2) + 2.D0*m_SC(3)*R + 3.D0*m_SC(4)*R*R)
       else
          call EXP6_Interaction(R,POTR, FPOTR)
       end if

      return
  end subroutine NN_Interaction
  !**********************************************************************************

  !**********************************************************************************
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
  !**********************************************************************************

  !**********************************************************************************
  subroutine Smooth_Connect_ZBL_TOOL()
  !***  PURPOSE:   to smoothly connect exp6 to ZBL potential
  !     INPUT:     FORCE1,  the force function in short distance
  !                R1,      the connection point (in cm)
  !                FORCE2,  the force function in large distance
  !                R2,      the second connection point (in cm)
  !     OUTPUT:    C,       the fitting parameters
  !
  !*********************************************************************************
  implicit none
      real(KINDDF)::R1 = 1.0D-8
      real(KINDDF)::R2 = 1.3D-8
      real(KINDDF)::C(4)

        real(KINDDF)::POTR1, POTR2, FPOTR1, FPOTR2

  !****
          !--- to do the fitting
          call ZBL_Interaction(R1, POTR1, FPOTR1)
          call EXP6_Interaction(R2, POTR2, FPOTR2)

          POTR1 = POTR1*2.D0
          POTR2 = POTR2*2.D0
          FPOTR1= FPOTR1*R1
          FPOTR2= FPOTR2*R2

          C(4)=(-(FPOTR1+FPOTR2)-2.D0*(POTR1-POTR2)/(R1-R2)) /(R1-R2)/(R1-R2)
          C(3)=-(FPOTR1-FPOTR2)/2.0D0/(R1-R2)-C(4)*1.5D0*(R1+R2)
          C(2)= -FPOTR1-2.0D0*C(3)*R1-3.0D0*C(4)*R1*R1
          C(1)= POTR1-R1*(C(2)+C(3)*R1+C(4)*R1*R1)


          !--- to test
          m_RZBL  = R1
          m_REXP6 = R2
          m_SC    = C

          print *, "R1= ", R1
          print *, "R2= ", R2
          print *, "C=  ", C(1:4)

          print *, "ZBL  at R1= ", POTR1*C_HALF, FPOTR1/R1
          print *, "EXP6 at R2= ", POTR2*C_HALF, FPOTR2/R2
          POTR1  =   (m_SC(1) + m_SC(2)*R1 + m_SC(3)*R1*R1 + m_SC(4)*R1*R1*R1)*C_HALF
          FPOTR1 =  -(m_SC(2) + 2.D0*m_SC(3)*R1 + 3.D0*m_SC(4)*R1*R1)/R1
          POTR2  =   (m_SC(1) + m_SC(2)*R2 + m_SC(3)*R2*R2 + m_SC(4)*R2*R2*R2)*C_HALF
          FPOTR2 =  -(m_SC(2) + 2.D0*m_SC(3)*R2 + 3.D0*m_SC(4)*R2*R2)/R2
          print *, "Connection at R1", POTR1, FPOTR1
          print *, "Connection at R2", POTR2, FPOTR2

          return
  end subroutine Smooth_Connect_ZBL_TOOL
  !*********************************************************************************
 end module Pot_ZBL_EXP6_HeHe
