 module Pot_ZBL_ABINTIO_He_W_930
 !***  this module is to calculate the potential and force between W-W W-He
 !     atoms. Supplied by Wang Jun (2011-09-30)

 use MD_CONSTANTS
     implicit none

     !***  the parameter block for nuclear-nuclear part in EAM ***************
     real(KINDDF), parameter,private::IGCLONG   = 8.854187817D-12
     real(KINDDF), parameter,private::ELECHARGE  = 1.60219D-19
     real(KINDDF), parameter,private::z1=2.0D0
     real(KINDDF), parameter,private::z2=74.0D0
     real(KINDDF), parameter,private::a0=5.291772108D-11

     real(KINDDF), parameter,private::CUTRR=5.8D0
     real(KINDDF), parameter,private::CUTORIN=1.5D0


    real(KINDDF), parameter,private::EXP6A=0.0084600D0
    real(KINDDF), parameter,private::EXP6B=9.6004204435D0
    real(KINDDF), parameter,private::EXP6C=4.135D0

    real(KINDDF), parameter,private::MIDPOINT03=1.57783601880074D0
    real(KINDDF), parameter,private::MIDPOINT02=1.94810507446527D0
    real(KINDDF), parameter,private::MIDPOINT01=2.10799896717072D0
    real(KINDDF), parameter,private::MIDPOINT=2.22277348041534D0
    real(KINDDF), parameter,private::MIDPOINTadd=2.32289959788322D0
    real(KINDDF), parameter,private::MIDPOINT2=3.05974050760269D0

    real(KINDDF), parameter,private::MMA0=-54.6920134504158D0
    real(KINDDF), parameter,private::MMB0=147.764553002780D0
    real(KINDDF), parameter,private::MMC0=-106.604073601512D0
    real(KINDDF), parameter,private::MMD0=23.3457856042901D0

    real(KINDDF), parameter,private::MMA=8.13894099222489D0
    real(KINDDF), parameter,private::MMB=-6.77821270976276D0
    real(KINDDF), parameter,private::MMC= 2.38825082778931D0
    real(KINDDF), parameter,private::MMD=-0.378954410552979D0

    real(KINDDF), parameter,private::MA2=11.8316959842284D0
    real(KINDDF), parameter,private::MB2=-9.92672851472125D0
    real(KINDDF), parameter,private::MC2=2.88240909576416D0
    real(KINDDF), parameter,private::MD2=-0.299054384231567D0

    real(KINDDF), parameter,private::MA1=18.5956897803305D0
    real(KINDDF), parameter,private::MB1= -18.0913941109865D0
    real(KINDDF), parameter,private::MC1=6.12169742584229D0
    real(KINDDF), parameter,private::MD1=-0.719758272171021D0

    real(KINDDF), parameter,private::MA0=11.5965650160971D0
    real(KINDDF), parameter,private::MB0=-9.38901671011024D0
    real(KINDDF), parameter,private::MC0=2.52039670944214D0
    real(KINDDF), parameter,private::MD0=-0.223789215087891D0

    real(KINDDF), parameter,private::MA=11.0095425643702D0
    real(KINDDF), parameter,private::MB=-9.03994004997097D0
    real(KINDDF), parameter,private::MC=2.48033046722412D0
    real(KINDDF), parameter,private::MD=-0.227488279342651D0

 contains

 subroutine ZBL_ABINTIO_NN_Interaction(r,POTR, FPOTR)
 !***  PURPOSE: to calculate the potential and force by Born-Mayer model
 !     INPUT:   r, the distance between two particle
 !     OUTPUT:  POTR, FPOTR

 implicit none
     real(KINDDF),intent(in)::r
     real(KINDDF),intent(out)::POTR, FPOTR
     real(KINDDF)::rr,POTRMID,fc,dfc, a,x, faix,dfaix,ZBLPOTEN


      rr=r*1.0D8

      if((r.ge.0.0D0).and.(r.lt.MIDPOINT03*1.0D-8)) then

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

      else if((r.ge.MIDPOINT03*1.0D-8).and.(r.lt.MIDPOINT02*1.0D-8)) then
          POTR  =MMA0+MMB0*rr+MMC0*rr**2.0D0+MMD0*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR

          FPOTR =(MMB0+2.0D0*MMC0*rr+3.0D0*MMD0*rr**2.0D0)
         !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8

      else if((r.ge.MIDPOINT02*1.0D-8).and.(r.lt.MIDPOINT01*1.0D-8)) then
          POTR  =MMA+MMB*rr+MMC*rr**2.0D0+MMD*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR

          FPOTR =(MMB+2.0D0*MMC*rr+3.0D0*MMD*rr**2.0D0)
         !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8

      else if((r.ge.MIDPOINT01*1.0D-8).and.(r.lt.MIDPOINT*1.0D-8)) then
          POTR  =MA2+MB2*rr+MC2*rr**2.0D0+MD2*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR

          FPOTR =(MB2+2.0D0*MC2*rr+3.0D0*MD2*rr**2.0D0)
         !modified accordding to NOTE-2014-12-04 in forcetable:
         ! FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8

      else if((r.ge.MIDPOINT*1.0D-8).and.(r.lt.MIDPOINTadd*1.0D-8))  then
          POTR  =MA1+MB1*rr+MC1*rr**2.0D0+MD1*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR

          FPOTR =(MB1+2.0D0*MC1*rr+3.0D0*MD1*rr**2.0D0)
         !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8

      else if((r.ge.MIDPOINTadd*1.0D-8).and.(r.lt.MIDPOINT2*1.0D-8))  then
          POTR  =MA0+MB0*rr+MC0*rr**2.0D0+MD0*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR

          FPOTR =(MB0+2.0D0*MC0*rr+3.0D0*MD0*rr**2.0D0)
         !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8

      else if((r.ge.MIDPOINT2*1.0D-8).and.(r.lt.3.5D-8))  then
          POTR  =MA+MB*rr+MC*rr**2.0D0+MD*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR

          FPOTR =(MB+2.0D0*MC*rr+3.0D0*MD*rr**2.0D0)
         !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8
      else  if(r.ge.3.5D-8)  then

         ! --------correct
         if(RR.le.(CUTRR-CUTORIN)) then
               fc=1.D0
               dfc=0.D0
        else if((CUTRR+CUTORIN).le.RR) then
               fc=0.D0
               dfc=0.D0
        else
               fc=3.141592654D0*0.5D0*(RR-CUTRR)/CUTORIN
               dfc=-0.5D0*dcos(fc)
               fc=0.5D0*(1.0D0-DSIN(fc))
        endif
        !------CORRECT
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
      end if


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

 end module Pot_ZBL_ABINTIO_He_W_930
