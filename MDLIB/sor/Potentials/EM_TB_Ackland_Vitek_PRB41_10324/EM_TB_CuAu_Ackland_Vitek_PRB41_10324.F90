  module EM_TB_CuAu_Ackland_Vitek_PRB41_10324
  !***  this module is to calculate the potential and force between Cu-Au atoms
  !REFERENCE: G.J.Ackland and V.Vitek, Phys.Rev.B41, (1990) 10324

   use MD_CONSTANTS
   use EM_TB_CuCu_Ackland_Vitek_PRB41_10324, only: RCT1CU=>RCT1, RCT2CU=>RCT2,     &
                                                   AT1CU=>AT1, AT2CU=>AT2
   use EM_TB_AuAu_Ackland_Vitek_PRB41_10324, only: RCT1AU=>RCT1, RCT2AU=>RCT2,     &
                                                   AT1AU=>AT1, AT2AU=>AT2
   implicit none

   !***  the parameter block for nuclear-nuclear part in the potential ***************
        real(KINDDF),parameter, private::AR1 = -0.0855455166D24*CP_EVERG*C_HALF
        real(KINDDF),parameter, private::AR2 =  0.1928358800D24*CP_EVERG*C_HALF
        real(KINDDF),parameter, private::AR3 =  0.7593228600D24*CP_EVERG*C_HALF
        real(KINDDF),parameter, private::AR4 =  0.D0
        real(KINDDF),parameter, private::AR5 =  0.D0
        real(KINDDF),parameter, private::AR6 =  0.D0
        real(KINDDF),parameter, private::RCR1= 4.3098160D-08
        real(KINDDF),parameter, private::RCR2= 4.0474790D-08
        real(KINDDF),parameter, private::RCR3= 3.2979464D-08
        real(KINDDF),parameter, private::RCR4= 0.D0
        real(KINDDF),parameter, private::RCR5= 0.D0
        real(KINDDF),parameter, private::RCR6= 0.D0

  !******* Parameter block for potential
        real(KINDDF),parameter, private::R11 = 1.85D-8
        real(KINDDF),parameter, private::R21 = 2.05D-8

  contains

  subroutine Moliere_Pot(r,POTR,FPOTR)
  !***  PURPOSE: to calculate the potential and force by Moliere model for r <r11
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR,FPOTR
      !******* parameters: *****************************************
      real(KINDDF), parameter::alpha(3)=(/0.35D0, 0.55D0, 0.1D0/)
      real(KINDDF), parameter::as = 0.0745D0*1.D-8
      real(KINDDF), parameter::zk = 9.0d9
      real(KINDDF), parameter::e  = 1.6029d-19
      real(KINDDF), parameter::z1  = 29.D0
      real(KINDDF), parameter::z2  = 79.D0
      real(KINDDF), parameter::zze2 = zk*z1*z2*e*e*1.D9
      !**************************************************************

      real(KINDDF)::be1, be2, be3, TPEXP

          be1=DEXP(-0.3D0*R/as)
          be2=(be1)**4.D0
          be3=(be2)**5.D0
          TPEXP=zze2*(alpha(1)*be1+alpha(2)*be2+alpha(3)*be3)/R
          POTR = TPEXP/2.D0
          !modified accordding to NOTE-2014-12-04 in forcetable
          !FPOTR = TPEXP/(R*R) + zze2*(alpha(1)*0.3*be1+            &
          !        alpha(2)*1.2D0*be2+alpha(3)*6.D0*be3)/(R*R*as)
          FPOTR = TPEXP/R + zze2*(alpha(1)*0.3*be1+            &
                  alpha(2)*1.2D0*be2+alpha(3)*6.D0*be3)/(R*as)

      return
  end subroutine Moliere_Pot

  subroutine First_Spline_Pot(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential for r11<r<r12 by spline
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !***   the parameter block *********************************
      real(KINDDF), parameter::ca0=0.172724d-08, ca1=-0.253573d0, ca2=0.124675d+08, ca3=-0.204896d+15
      !***********************************************************
      !  Local variable
       real(KINDDF)::TPEXP

          TPEXP = ca0 + ca1*r + ca2*r**2 + ca3*r**3
          POTR  = TPEXP/2.D0
          !modified accordding to NOTE-2014-12-04 in forcetable
          !FPOTR  = -(ca1 + 2.D0*ca2*r + 3.D0*ca3*r**2 ) / r
          FPOTR  = -(ca1 + 2.D0*ca2*r + 3.D0*ca3*r**2 )

      return
  end subroutine First_Spline_Pot

  subroutine TB_EM_NN_Pot(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the repulsive potential for r>r12 by EAM
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR

          POTR  = 0.D0
          FPOTR = 0.D0

          IF(R.LT.RCR6)THEN
             POTR=AR6*((RCR6-R)**3.D0)+                       &
                    AR5*((RCR5-R)**3.D0)+                     &
                    AR4*((RCR4-R)**3.D0)+                     &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)
            !modified accordding to NOTE-2014-12-04 in forcetable
             !FPOTR=(C_SIX/R)*AR6*(RCR6-R)**2.D0+                &
             !        (C_SIX/R)*AR5*(RCR5-R)**2.D0+              &
             !        (C_SIX/R)*AR4*(RCR4-R)**2.D0+              &
             !        (C_SIX/R)*AR3*(RCR3-R)**2.D0+              &
             !        (C_SIX/R)*AR2*(RCR2-R)**2.D0+              &
             !        (C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR6*(RCR6-R)**2.D0+                &
                    C_SIX*AR5*(RCR5-R)**2.D0+                &
                    C_SIX*AR4*(RCR4-R)**2.D0+                &
                    C_SIX*AR3*(RCR3-R)**2.D0+                &
                    C_SIX*AR2*(RCR2-R)**2.D0+                &
                    C_SIX*AR1*(RCR1-R)**2.D0

          ENDIF

          IF(R.GE.RCR6.AND.R.LT.RCR5)THEN
             POTR=AR5*((RCR5-R)**3.D0)+                       &
                    AR4*((RCR4-R)**3.D0)+                     &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)

            !modified accordding to NOTE-2014-12-04 in forcetable
            ! FPOTR=(C_SIX/R)*AR5*(RCR5-R)**2.D0+                &
            !         (C_SIX/R)*AR4*(RCR4-R)**2.D0+              &
            !         (C_SIX/R)*AR3*(RCR3-R)**2.D0+              &
            !         (C_SIX/R)*AR2*(RCR2-R)**2.D0+              &
            !         (C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR5*(RCR5-R)**2.D0+                &
                    C_SIX*AR4*(RCR4-R)**2.D0+                &
                    C_SIX*AR3*(RCR3-R)**2.D0+                &
                    C_SIX*AR2*(RCR2-R)**2.D0+                &
                    C_SIX*AR1*(RCR1-R)**2.D0

          ENDIF

          IF(R.GE.RCR5.AND.R.LT.RCR4)THEN
             POTR=AR4*((RCR4-R)**3.D0)+                       &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)

            !modified accordding to NOTE-2014-12-04 in forcetable
            ! FPOTR=(C_SIX/R)*AR4*(RCR4-R)**2.D0+                &
            !         (C_SIX/R)*AR3*(RCR3-R)**2.D0+              &
            !         (C_SIX/R)*AR2*(RCR2-R)**2.D0+              &
            !         (C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR4*(RCR4-R)**2.D0+                &
                    C_SIX*AR3*(RCR3-R)**2.D0+                &
                    C_SIX*AR2*(RCR2-R)**2.D0+                &
                    C_SIX*AR1*(RCR1-R)**2.D0

          ENDIF

          IF(R.GE.RCR4.AND.R.LT.RCR3)THEN
             POTR=AR3*((RCR3-R)**3.D0)+                       &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)

            !modified accordding to NOTE-2014-12-04 in forcetable
            ! FPOTR=(C_SIX/R)*AR3*(RCR3-R)**2.D0+                &
            !         (C_SIX/R)*AR2*(RCR2-R)**2.D0+              &
            !         (C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR3*(RCR3-R)**2.D0+                &
                    C_SIX*AR2*(RCR2-R)**2.D0+                &
                    C_SIX*AR1*(RCR1-R)**2.D0

          ENDIF

          IF(R.GE.RCR3.AND.R.LT.RCR2)THEN
             POTR=AR2*((RCR2-R)**3.D0)+                       &
                    AR1*((RCR1-R)**3.D0)

            !modified accordding to NOTE-2014-12-04 in forcetable
             !FPOTR=(C_SIX/R)*AR2*(RCR2-R)**2.D0+                &
             !        (C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR2*(RCR2-R)**2.D0+                &
                    C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR2.AND.R.LT.RCR1)THEN
             POTR=AR1*((RCR1-R)**3.D0)
            !modified accordding to NOTE-2014-12-04 in forcetable
            ! FPOTR=(C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF
      return
  end subroutine TB_EM_NN_Pot

  subroutine TB_EM_NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the nuclear-nuclear interaction
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR

        if ( r .lt. r11 ) then
             call Moliere_Pot(r,POTR, FPOTR)
        else if( r .lt. r21) then
             call First_Spline_Pot(r,POTR, FPOTR)
        else
             call TB_EM_NN_Pot(r,POTR, FPOTR)
        endif
        return
  end subroutine TB_EM_NN_Interaction


  subroutine TB_EM_NE_Interaction(r,POTB,FPOTB)
  !***  PURPOSE: to calculate the attractive potential by EAM
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB
      !Local varaibles
      real(KINDDF)::X1,Y1,X2,Y2

          X1 = (RCT1cu-R)
          Y1 = (RCT2cu-R)
          X2 = (RCT1au-R)
          Y2 = (RCT2au-R)
          IF(R.GT.C_ZERO.AND.R.LT.RCT2au)THEN
             POTB = AT1au*AT1cu*(X2*X1)**3.D0+                     &
                    AT1au*AT2cu*(X2*Y1)**3.D0+                     &
                    AT2au*AT1cu*(Y2*X1)**3.D0+                     &
                    AT2au*AT2cu*(Y2*Y1)**3.D0
             POTB = DSQRT(POTB)
             FPOTB = 1.5D0*AT1au*AT1cu*(X2+X1)*(X2*X1)**2.D0+      &
                     1.5D0*AT1au*AT2cu*(X2+Y1)*(X2*Y1)**2.D0+      &
                     1.5D0*AT2au*AT1cu*(Y2+X1)*(Y2*X1)**2.D0+      &
                     1.5D0*AT2au*AT2cu*(Y2+Y1)*(Y2*Y1)**2.D0
            !modified accordding to NOTE-2014-12-04 in forcetable
             !FPOTB = (FPOTB/POTB)*C_HALF/R
          ENDIF

          IF(R.GE.RCT2au.AND.R.LT.RCT2cu)THEN
             POTB = AT1au*AT1cu*(X2*X1)**3.D0+                     &
                    AT1au*AT2cu*(X2*Y1)**3.D0
             POTB  = DSQRT(POTB)
             FPOTB = 1.5D0*AT1au*AT1cu*(X2+X1)*(X2*X1)**2.D0+      &
                     1.5D0*AT1au*AT2cu*(X2+Y1)*(X2*Y1)**2.D0
            !modified accordding to NOTE-2014-12-04 in forcetable
             !FPOTB = (FPOTB/POTB)*C_HALF/R
          ENDIF

          IF(R.GE.RCT2cu.AND.R.LT.RCT1cu)THEN
             POTB = AT1au*AT1cu*(X2*X1)**3.D0
             POTB = DSQRT(POTB)
             FPOTB =1.5D0*AT1au*AT1cu*(X2+X1)*(X2*X1)**2.D0
            !modified accordding to NOTE-2014-12-04 in forcetable
             !FPOTB =(FPOTB/POTB)*C_HALF/R
          ENDIF

          return
  end subroutine TB_EM_NE_Interaction

  end module EM_TB_CuAu_Ackland_Vitek_PRB41_10324
