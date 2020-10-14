  module EM_TB_ZrZr_AWB_PHIL_MAG_A71_1995_533
  !***  this module is to calculate the potentialand force between Zr-Zr atoms
  !
  !     REFERENCE: G.J.Ackland, S. J. Wooding and D. J. Bacon, Philo. Mag. A71 (1995)553-565
  !                U. Pinsook and G. J. Ackland, Phys. Rev. B58(1998)11252
  !
  !     Author:   HOU Qing, Oct, 2012

  use MD_CONSTANTS
        implicit none


        real(KINDDF), private, parameter::RRFIT= 3.249D-08
        real(KINDDF), private, parameter::ZrZ   = 40.00D0                        !atomic number of Zr

        !***  the parameter block for nuclear-nuclear part  ***************
        real(KINDDF),parameter, private::AR1=-58.480671*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter, private::AR2= 83.692945*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter, private::AR3=-20.958026*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter, private::AR4= -1.309533*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter, private::AR5= 65.726356*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter, private::AR6=139.411255*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        !In original paper of AWB, AR7 is 700, the potential is not continue at R21
        !thus we change the parameter to 433
        !real(KINDDF),parameter, private::AR7=700.000000*CP_EVERG/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter, private::AR7=443.442749*CP_EVERG/((RRFIT*CP_RSQ2)**3)

        real(KINDDF),parameter, private::RCR1=1.22*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter, private::RCR2=1.20*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter, private::RCR3=1.14*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter, private::RCR4=0.95*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter, private::RCR5=0.80*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter, private::RCR6=0.70*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter, private::RCR7=0.63*(RRFIT*CP_RSQ2)

      !****** Parameter block for nuclear-electron part *******************
        real(KINDDF),parameter::RCT1 = 1.22*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter::RCT2 = 1.05*(RRFIT*CP_RSQ2)
        real(KINDDF),parameter::AT1  = 48.288836*(CP_EVERG**2)/((RRFIT*CP_RSQ2)**3)
        real(KINDDF),parameter::AT2  = -0.850556*(CP_EVERG**2)/((RRFIT*CP_RSQ2)**3)

       !******* Parameter block for small distance (ZBL+spline) potential
        real(KINDDF),parameter, private::R11 = 1.20D-8
        real(KINDDF),parameter, private::R21 = 2.55D-8

        private::First_Spline_Pot

  contains
  !*********************************************************************


  !*********************************************************************
  subroutine First_Spline_Pot(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential for r11<r<r12 by spline
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !***   the parameter block *********************************
      real(KINDDF), parameter::B1=13.4037, B2=-12.2190, B3=5.0261, B4=-0.8435
      !***********************************************************
      !  Local variable
       real(KINDDF)::RA, TPEXP

          RA=R*CP_CM2A
          TPEXP = CP_EVERG*DEXP(B1 + B2*RA + B3*RA**2 + B4*RA**3)
          POTR  = TPEXP/2.D0
          !modified accordding to NOTE-2014-12-04 in forcetable
          !FPOTR  = -(B2 + 2.D0*B3*RA + 3.D0*B4*RA**2 )*CP_CM2A*TPEXP/R
          FPOTR  = -(B2 + 2.D0*B3*RA + 3.D0*B4*RA**2 )*CP_CM2A*TPEXP

      return
  end subroutine First_Spline_Pot
  !*********************************************************************

  !*********************************************************************
  subroutine TB_EM_NN_Pot(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the repulsive potential for r>r12 by EAM
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
  !--- local variables
      real(KINDDF)::DR

         !modified accordding to NOTE-2014-12-04 in forcetable
         ! change DR = 3.D0/R to DR = 3.D0
         !DR = 3.D0/R
          DR = 3.D0
          IF(R.LT.RCR7)THEN
           POTR=    AR7*((RCR7-R)**3.D0)+                     &
                    AR6*((RCR6-R)**3.D0)+                     &
                    AR5*((RCR5-R)**3.D0)+                     &
                    AR4*((RCR4-R)**3.D0)+                     &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)
           FPOTR=   DR*AR7*(RCR7-R)**2.D0+                    &
                    DR*AR6*(RCR6-R)**2.D0+                    &
                    DR*AR5*(RCR5-R)**2.D0+                    &
                    DR*AR4*(RCR4-R)**2.D0+                    &
                    DR*AR3*(RCR3-R)**2.D0+                    &
                    DR*AR2*(RCR2-R)**2.D0+                    &
                    DR*AR1*(RCR1-R)**2.D0
          ENDIF


          IF(R.GE.RCR7.AND.R.LT.RCR6)THEN
           POTR=    AR6*((RCR6-R)**3.D0)+                     &
                    AR5*((RCR5-R)**3.D0)+                     &
                    AR4*((RCR4-R)**3.D0)+                     &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)
           FPOTR=   DR*AR6*(RCR6-R)**2.D0+                    &
                    DR*AR5*(RCR5-R)**2.D0+                    &
                    DR*AR4*(RCR4-R)**2.D0+                    &
                    DR*AR3*(RCR3-R)**2.D0+                    &
                    DR*AR2*(RCR2-R)**2.D0+                    &
                    DR*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR6.AND.R.LT.RCR5)THEN
           POTR=AR5*((RCR5-R)**3.D0)+                         &
                    AR4*((RCR4-R)**3.D0)+                     &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)
           FPOTR=   DR*AR5*(RCR5-R)**2.D0+                    &
                    DR*AR4*(RCR4-R)**2.D0+                    &
                    DR*AR3*(RCR3-R)**2.D0+                    &
                    DR*AR2*(RCR2-R)**2.D0+                    &
                    DR*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR5.AND.R.LT.RCR4)THEN
           POTR=AR4*((RCR4-R)**3.D0)+                         &
                    AR3*((RCR3-R)**3.D0)+                     &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)
           FPOTR=DR*AR4*(RCR4-R)**2.D0+                       &
                     DR*AR3*(RCR3-R)**2.D0+                   &
                     DR*AR2*(RCR2-R)**2.D0+                   &
                     DR*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR4.AND.R.LT.RCR3)THEN
           POTR=AR3*((RCR3-R)**3.D0)+                         &
                    AR2*((RCR2-R)**3.D0)+                     &
                    AR1*((RCR1-R)**3.D0)
           FPOTR=DR*AR3*(RCR3-R)**2.D0+                       &
                     DR*AR2*(RCR2-R)**2.D0+                   &
                     DR*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR3.AND.R.LT.RCR2)THEN
           POTR=AR2*((RCR2-R)**3.D0)+                         &
                    AR1*((RCR1-R)**3.D0)
           FPOTR=DR*AR2*(RCR2-R)**2.D0+                       &
                     DR*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR2.AND.R.LT.RCR1)THEN
             POTR=AR1*((RCR1-R)**3.D0)
             FPOTR=DR*AR1*(RCR1-R)**2.D0
          ENDIF

          POTR= C_HALF*POTR
      return
  end subroutine TB_EM_NN_Pot
  !*********************************************************************

  !*********************************************************************
  subroutine Second_Spline_Pot(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential for r21<r<RCR7 by spline
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !***   the parameter block *********************************
      !***********************************************************
      !  Local variable
       real(KINDDF)::POTR1, POTR2, FPOTR1, FPOTR2
       real(KINDDF)::C3=-19388949878932.86D0, C2=1706444.034874689D0, C1=-5.0454065151782766D-002, &
                     C0=5.0256487041006336D-010

          POTR =   (C0 + C1*R + C2*R**2.D0 + C3*R**3.D0 )*C_HALF
          !modified accordding to NOTE-2014-12-04 in forcetable
          !FPOTR = -(C1 + 2.D0*C2*R + 3.D0*C3*R**2.D0)/R
          FPOTR = -(C1 + 2.D0*C2*R + 3.D0*C3*R**2.D0)

      return
  end subroutine Second_Spline_Pot
  !*********************************************************************


  !*********************************************************************
  subroutine TB_EM_NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the nuclear-nuclear interaction
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  !
  use MD_CorePotent
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR

        if ( r .lt. r11 ) then
            call ZBL_Pot(ZrZ, ZrZ, r,POTR, FPOTR)
        else if( r .le. r21) then
                 call First_Spline_Pot(r,POTR, FPOTR)
        else if(r  .le. RCR7) then
                 call Second_Spline_Pot(r,POTR, FPOTR)
         else
                 call TB_EM_NN_Pot(r,POTR, FPOTR)
        endif

        !*** to test the Second_Spline_Pot, deleted
        !    after testing
        !call First_Spline_Pot(r21,POTR, FPOTR)
        !print *,  POTR, FPOTR
        !call Second_Spline_Pot(r21,POTR, FPOTR)
        !print *,  POTR, FPOTR

        !call TB_EM_NN_Pot(RCR7,POTR, FPOTR)
        !print *,  POTR, FPOTR
        !call Second_Spline_Pot(RCR7,POTR, FPOTR)
        !print *,  POTR, FPOTR
        !pause

        return
  end subroutine TB_EM_NN_Interaction
  !*********************************************************************

  !*********************************************************************
  subroutine TB_EM_NE_Interaction(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the attractive potential by EAM
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB

          IF(R.GT.C_ZERO.AND.R.LT.RCT2)THEN
            POTB=AT2*(RCT2-R)**3.D0 + AT1*(RCT1-R)**3.D0
            !modified accordding to NOTE-2014-12-04 in forcetable
            !FPOTB=(C_THR*C_HALF/R)*AT2*(RCT2-R)**2.D0 + (C_THR*C_HALF/R)*AT1*(RCT1-R)**2.D0
            FPOTB= C_THR*AT2*(RCT2-R)**2.D0 + C_THR*AT1*(RCT1-R)**2.D0
          ENDIF

          IF(R.GE.RCT2.AND.R.LT.RCT1)THEN
            POTB=AT1*(RCT1-R)**3.D0
            !modified accordding to NOTE-2014-12-04 in forcetable
            !FPOTB=(C_THR*C_HALF/R)*AT1*(RCT1-R)**2.D0
            FPOTB= C_THR*AT1*(RCT1-R)**2.D0
          ENDIF

          return
  end subroutine TB_EM_NE_Interaction
  !*********************************************************************


  end module EM_TB_ZrZr_AWB_PHIL_MAG_A71_1995_533
