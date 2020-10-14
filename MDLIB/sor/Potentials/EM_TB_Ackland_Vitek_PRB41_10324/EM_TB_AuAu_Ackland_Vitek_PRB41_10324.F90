  module EM_TB_AuAu_Ackland_Vitek_PRB41_10324
  !***  this module is to calculate the potential and force between Au-Au atoms
  !REFERENCE: G.J.Ackland and V.Vitek, Phys.Rev.B41, (1990) 10324

  use MD_CONSTANTS
        implicit none


        real(KINDDF),parameter::RRFIT= 4.078D-08
        !***  the parameter block for nuclear-nuclear part in the potential ***********
        real(KINDDF),parameter, private::AR1= 29.059066*CP_EVERG*C_HALF/(RRFIT**3)
        real(KINDDF),parameter, private::AR2=-153.14779*CP_EVERG*C_HALF/(RRFIT**3)
        real(KINDDF),parameter, private::AR3= 148.17881*CP_EVERG*C_HALF/(RRFIT**3)
        real(KINDDF),parameter, private::AR4= -22.20508*CP_EVERG*C_HALF/(RRFIT**3)
        real(KINDDF),parameter, private::AR5= 72.71465*CP_EVERG*C_HALF/(RRFIT**3)
        real(KINDDF),parameter, private::AR6= 199.26269*CP_EVERG*C_HALF/(RRFIT**3)
        real(KINDDF),parameter, private::RCR1=1.2247449*RRFIT
        real(KINDDF),parameter, private::RCR2=1.1547054*RRFIT
        real(KINDDF),parameter, private::RCR3=1.1180065*RRFIT
        real(KINDDF),parameter, private::RCR4=1.0000000*RRFIT
        real(KINDDF),parameter, private::RCR5=0.8660254*RRFIT
        real(KINDDF),parameter, private::RCR6=0.7071068*RRFIT

      !****** Parameter block for nuclear-electron part the potential *****************
        real(KINDDF),parameter::RCT1 = 1.1180065*RRFIT
        real(KINDDF),parameter::RCT2 = 0.8660254*RRFIT
        real(KINDDF),parameter::AT1 = 21.930125*(CP_EVERG**2)/(RRFIT**3)
        real(KINDDF),parameter::AT2 = 284.99631*(CP_EVERG**2)/(RRFIT**3)

       !******* Parameter block for potential
        real(KINDDF),parameter, private::R11 = 1.85D-8
        real(KINDDF),parameter, private::R21 = 2.10D-8


  contains
  !*********************************************************************

  !*********************************************************************
  subroutine Moliere_Pot(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by Moliere model for r <r11
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !******* parameters: *****************************************
      real(KINDDF), parameter::alpha(3)=(/0.35D0, 0.55D0, 0.1D0/)
      real(KINDDF), parameter::as = 0.0752D0*1.D-8
      real(KINDDF), parameter::zk = 9.0d9
      real(KINDDF), parameter::e  = 1.6029d-19
      real(KINDDF), parameter::z  = 79.D0
      real(KINDDF), parameter::zze2 = zk*z*z*e*e*1.D9
      !**************************************************************

      real(KINDDF)::be1, be2, be3, TPEXP

          be1=DEXP(-0.3D0*R/as)
          be2=(be1)**4.D0
          be3=(be2)**5.D0
          TPEXP=zze2*(alpha(1)*be1+alpha(2)*be2+alpha(3)*be3)/R
          POTR = TPEXP/2.D0
          !modified accordding to NOTE-2014-12-04 in forcetable
          !FPOTR = TPEXP/(R*R) + zze2*(alpha(1)*0.3*be1+       &
          !        alpha(2)*1.2D0*be2+alpha(3)*6.D0*be3)/(R*R*as)
          FPOTR = TPEXP/R + zze2*(alpha(1)*0.3*be1+       &
                  alpha(2)*1.2D0*be2+alpha(3)*6.D0*be3)/(R*as)

      return
  end subroutine Moliere_Pot
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
      real(KINDDF), parameter::ca0=0.333482d-08, ca1=-0.476248d0, ca2=0.228187d+08, ca3=-0.365939d+15
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

          POTR  = 0.D0
          FPOTR = 0.D0

          IF(R.LT.RCR6)THEN
           POTR=AR6*((RCR6-R)**3.D0)+                              &
                    AR5*((RCR5-R)**3.D0)+                          &
                    AR4*((RCR4-R)**3.D0)+                          &
                    AR3*((RCR3-R)**3.D0)+                          &
                    AR2*((RCR2-R)**3.D0)+                          &
                    AR1*((RCR1-R)**3.D0)

          !modified accordding to NOTE-2014-12-04 in forcetable
           !FPOTR=(C_SIX/R)*AR6*(RCR6-R)**2.D0+                          &
           !          (C_SIX/R)*AR5*(RCR5-R)**2.D0+                      &
           !          (C_SIX/R)*AR4*(RCR4-R)**2.D0+                      &
           !          (C_SIX/R)*AR3*(RCR3-R)**2.D0+                      &
           !          (C_SIX/R)*AR2*(RCR2-R)**2.D0+                      &
           !          (C_SIX/R)*AR1*(RCR1-R)**2.D0
           FPOTR = C_SIX*AR6*(RCR6-R)**2.D0+                      &
                   C_SIX*AR5*(RCR5-R)**2.D0+                      &
                   C_SIX*AR4*(RCR4-R)**2.D0+                      &
                   C_SIX*AR3*(RCR3-R)**2.D0+                      &
                   C_SIX*AR2*(RCR2-R)**2.D0+                      &
                   C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR6.AND.R.LT.RCR5)THEN
           POTR=AR5*((RCR5-R)**3.D0)+                             &
                    AR4*((RCR4-R)**3.D0)+                         &
                    AR3*((RCR3-R)**3.D0)+                         &
                    AR2*((RCR2-R)**3.D0)+                         &
                    AR1*((RCR1-R)**3.D0)
          !modified accordding to NOTE-2014-12-04 in forcetable
          ! FPOTR=(C_SIX/R)*AR5*(RCR5-R)**2.D0+                          &
          !           (C_SIX/R)*AR4*(RCR4-R)**2.D0+                      &
          !           (C_SIX/R)*AR3*(RCR3-R)**2.D0+                      &
          !           (C_SIX/R)*AR2*(RCR2-R)**2.D0+                      &
          !           (C_SIX/R)*AR1*(RCR1-R)**2.D0
           FPOTR = C_SIX*AR5*(RCR5-R)**2.D0+                      &
                   C_SIX*AR4*(RCR4-R)**2.D0+                      &
                   C_SIX*AR3*(RCR3-R)**2.D0+                      &
                   C_SIX*AR2*(RCR2-R)**2.D0+                      &
                   C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR5.AND.R.LT.RCR4)THEN
           POTR=AR4*((RCR4-R)**3.D0)+                             &
                    AR3*((RCR3-R)**3.D0)+                         &
                    AR2*((RCR2-R)**3.D0)+                         &
                    AR1*((RCR1-R)**3.D0)
          !modified accordding to NOTE-2014-12-04 in forcetable
           !FPOTR=(C_SIX/R)*AR4*(RCR4-R)**2.D0+                          &
           !          (C_SIX/R)*AR3*(RCR3-R)**2.D0+                      &
           !          (C_SIX/R)*AR2*(RCR2-R)**2.D0+                      &
           !          (C_SIX/R)*AR1*(RCR1-R)**2.D0
           FPOTR= C_SIX*AR4*(RCR4-R)**2.D0+                      &
                  C_SIX*AR3*(RCR3-R)**2.D0+                      &
                  C_SIX*AR2*(RCR2-R)**2.D0+                      &
                  C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR4.AND.R.LT.RCR3)THEN
           POTR=AR3*((RCR3-R)**3.D0)+                            &
                    AR2*((RCR2-R)**3.D0)+                        &
                    AR1*((RCR1-R)**3.D0)
          !modified accordding to NOTE-2014-12-04 in forcetable
           !FPOTR=(C_SIX/R)*AR3*(RCR3-R)**2.D0+                          &
           !          (C_SIX/R)*AR2*(RCR2-R)**2.D0+                      &
           !          (C_SIX/R)*AR1*(RCR1-R)**2.D0
          FPOTR= C_SIX*AR3*(RCR3-R)**2.D0+                       &
                 C_SIX*AR2*(RCR2-R)**2.D0+                       &
                 C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR3.AND.R.LT.RCR2)THEN
           POTR=AR2*((RCR2-R)**3.D0)+                            &
                    AR1*((RCR1-R)**3.D0)
          !modified accordding to NOTE-2014-12-04 in forcetable
           !FPOTR=(C_SIX/R)*AR2*(RCR2-R)**2.D0+                          &
           !          (C_SIX/R)*AR1*(RCR1-R)**2.D0
           FPOTR= C_SIX*AR2*(RCR2-R)**2.D0+                      &
                  C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF

          IF(R.GE.RCR2.AND.R.LT.RCR1)THEN
             POTR=AR1*((RCR1-R)**3.D0)
          !modified accordding to NOTE-2014-12-04 in forcetable
             !FPOTR=(C_SIX/R)*AR1*(RCR1-R)**2.D0
             FPOTR= C_SIX*AR1*(RCR1-R)**2.D0
          ENDIF
      return
  end subroutine TB_EM_NN_Pot
  !*********************************************************************

  !*********************************************************************
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

          POTB  = 0.D0
          FPOTB = 0.D0
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
            FPOTB = C_THR*AT1*(RCT1-R)**2.D0
          ENDIF

          return
  end subroutine TB_EM_NE_Interaction
  !*********************************************************************


  end module EM_TB_AuAu_Ackland_Vitek_PRB41_10324
