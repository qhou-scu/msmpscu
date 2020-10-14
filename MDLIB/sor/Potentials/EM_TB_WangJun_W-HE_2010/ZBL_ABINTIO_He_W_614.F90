 module Pot_ZBL_ABINTIO_He_W_614
 !***  this module is to calculate the potential and force between W-W W-He
 !     atoms. Supplied by Wang Jun (2011-06-14)
 !REFERENCE:

 use MD_CONSTANTS
      implicit none
   !real(KINDDF)::potentialv,a,r,faix,X,dfaix,dpotentialv,ZBLPOTEN,CUTRR=5.8D0,CUTORIN=1.5D0

     real(KINDDF), parameter,private::IGCLONG   = 8.854187817D-12
     real(KINDDF), parameter,private::ELECHARGE  = 1.60219D-19
     real(KINDDF), parameter,private::z1=2.0D0
     real(KINDDF), parameter,private::z2=74.0D0
     real(KINDDF), parameter,private::a0=5.291772108D-11

     real(KINDDF), parameter,private:: ZBLA=0.D0
     real(KINDDF), parameter,private:: ZBLB=0.D0
     real(KINDDF), parameter,private:: ZBLC=7599886.97245717D0
     real(KINDDF), parameter,private:: ZBLD=8.10683667659760D0
     real(KINDDF), parameter,private:: ZBLE=0.325999200344086D0

     real(KINDDF), parameter,private:: MA=122.160513080050D0
     real(KINDDF), parameter,private:: MB=-169.940064472736D0
     real(KINDDF), parameter,private:: MC=82.4270155974909D0
     real(KINDDF), parameter,private:: MD=-13.8730267156615D0

     real(KINDDF), parameter,private::MMA=0.656521958414464D0
     real(KINDDF), parameter,private::MMB=0.305144430204926D0
     real(KINDDF), parameter,private::MMC=-0.324284189377027D0
     real(KINDDF), parameter,private::MMD=5.243602313441372D-002


    real(KINDDF), parameter,private::EA=101.222166575487D0
    real(KINDDF), parameter,private::EB=2.46948345005512D0
    real(KINDDF), parameter,private::EC=2.009000388561760D-2


    real(KINDDF), parameter,private::EXP6A=0.0084600D0
    real(KINDDF), parameter,private::EXP6B=9.6004204435D0
    real(KINDDF), parameter,private::EXP6C=4.135D0

    real(KINDDF), parameter,private::CUTRR=5.8D0
    real(KINDDF), parameter,private::CUTORIN=1.5D0

 contains

 subroutine ZBL_ABINTIO_NN_Interaction(r,POTR, FPOTR)
 !***  PURPOSE: to calculate the nuclear- nclear interaction
 !     INPUT:   r, the distance between two particle
 !     OUTPUT:  POTR, FPOTR

 implicit none
     real(KINDDF)::r,fc,dfc !,intent(in)::r
     real(KINDDF),intent(out)::POTR, FPOTR
     real(KINDDF)::XX,RR,POTRMID
     real(KINDDF)::a,faix,dfaix,X,ZBLPOTEN

      rr=r*1.0D8

      if((r.ge.0.0D0).and.(r.lt.1.65D0*1.0D-8)) then
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
          FPOTR=-CP_EVERG*FPOTR

      ELSEIF((r.ge.1.65D0*1.D-8).and.(r.lt.1.978D0*1.D-8)) then
          POTR  =MA+MB*rr+MC*rr**2.0D0+MD*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR
          FPOTR =(MB+2.0D0*MC*rr+3.0D0*MD*rr**2.0D0)*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR/R
          FPOTR =-CP_EVERG*FPOTR
     ELSEIF((r.ge.1.978D0*1.D-8).and.(r.lt.2.20227904915810D-8)) then
          ZBLPOTEN=ZBLA/rr*dexp(-rr*ZBLB)
          POTR  =ZBLPOTEN+ZBLC*dexp(-rr*ZBLD)+ZBLE
          POTR  =0.5D0*CP_EVERG*POTR
          FPOTR =-ZBLA*dexp(-rr*ZBLB)*(1.0d0/rr/rr+ZBLB/RR)
          FPOTR=FPOTR-ZBLC*dexp(-rr*ZBLD)*ZBLD
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-FPOTR*CP_EVERG*1.0D8/R
            FPOTR =-FPOTR*CP_EVERG*1.0D8

     ELSEIF((r.ge.2.20227904915810D-8).and.(r.lt.2.73403934657574D-8))  then
          POTR  =EA*dexp(-rr*EB)+Ec
          POTR  =0.5D0*CP_EVERG*POTR
          FPOTR =-Ea*Eb*dexp(-rr*Eb)*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR/r
          FPOTR =-CP_EVERG*FPOTR
     ELSEIF((r.ge.2.73403934657574D-8).and.(r.lt.3.5D-8))  then
         ! --------correct
          POTR  =MMA+MMB*rr+MMC*rr**2.0D0+MMD*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR
          FPOTR =(MMB+2.0D0*MMC*rr+3.0D0*MMD*rr**2.0D0)
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8

     ELSEIF(r.ge.3.5D-8)  then

        if(RR.le.(CUTRR-CUTORIN)) then
           fc=1.D0
           dfc=0.D0
        else if((CUTRR+CUTORIN).le.RR)     then
           fc=0.D0
           dfc=0.D0
        else
           fc=3.141592654D0*0.5D0*(RR-CUTRR)/CUTORIN
           dfc=-0.5D0*dcos(fc)
           fc=0.5D0*(1.0D0-DSIN(fc))
        endif

         ZBLPOTEN=6.0D0/(EXP6B-6.0D0)*DEXP(EXP6B*(1.0D0-RR/EXP6C))
         POTRMID=EXP6A*(ZBLPOTEN-EXP6B/(EXP6B-6.0D0)*(EXP6C/RR)**6.0D0)
         POTR=POTRMID*fc
         POTR  =CP_EVERG*0.5D0*POTR

         FPOTR=-ZBLPOTEN*EXP6B/EXP6C+EXP6B/(EXP6B-6.0D0)*6.0D0*(EXP6C/RR)**5.0D0*EXP6C/RR/RR
         FPOTR=EXP6A*FPOTR
         FPOTR=FPOTR*FC+POTRMID*DFC
         !modified accordding to NOTE-2014-12-04 in forcetable:
         !FPOTR=-CP_EVERG*FPOTR*1.0D8/R
         FPOTR=-CP_EVERG*FPOTR*1.0D8
     END IF



     return
 end subroutine ZBL_ABINTIO_NN_Interaction


 subroutine ZBL_ABINTIO_NE_Interaction(r,POTB, FPOTB)
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
 end subroutine ZBL_ABINTIO_NE_Interaction

 end module Pot_ZBL_ABINTIO_He_W_614
