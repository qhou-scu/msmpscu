 module Pot_ZBL_ABINTIO_He_W_513
 !***  this module is to calculate the potential and force between W-W W-He
 !     atoms. Supplied by Wang Jun (2011-05-13)
 !REFERENCE:

 use MD_CONSTANTS
 !use MD90_Globle_Variables,only:ZBLA,ZBLB,ZBLC,ZBLD,ZBLe,ma,mb,mc,md,MIDPOINT,MIDPOINT2,ea,eb,ec
      implicit none

       !***  the parameter block for nuclear-nuclear part in EAM ***************
      !real(KINDDF)::potentialv,a,r,faix,X,dfaix,dpotentialv,ZBLPOTEN
      real(KINDDF), parameter,private::CUTRR=5.8D0
      real(KINDDF), parameter,private::CUTORIN=0.7D0

      real(KINDDF), parameter,private:: ZBLA=3127199.4113922D0
      real(KINDDF), parameter,private:: ZBLB=14.51015D0
      real(KINDDF), parameter,private:: ZBLC=7592776.41773224D0
      real(KINDDF), parameter,private:: ZBLD=8.10747385D0
      real(KINDDF), parameter,private:: ZBLE=0.3208359D0


      real(KINDDF), parameter,private::MA=1.09951027587204D0   !7.92041989584602D0 6.43823829937402
      real(KINDDF), parameter,private::MB=-0.180173550401268D0    !-5.46574252078757D0
      real(KINDDF), parameter,private::MC=-0.151558013341505D0    !1.56132350481559D0
      real(KINDDF), parameter,private::MD=3.237684631862348D-2    !-0.150061294708069D0


      real(KINDDF), parameter,private::EA=106.690696859600D0         !6.717046865181104D89
      real(KINDDF), parameter,private::EB=2.49999690055847D0         !8.569512844085693D-3
      real(KINDDF), parameter,private::EC=2.098244335294092D-2        !56.8629783391953D0


      real(KINDDF), parameter,private::LA=0.007D0
      real(KINDDF), parameter,private::LB=4.914989D0
      real(KINDDF), parameter,private::LC=8.261011D0
      real(KINDDF), parameter,private::LD=1.670001D0
      real(KINDDF), parameter,private::LE=6.738495D0
      real(KINDDF), parameter,private::LF=1.482181596472D0

 contains

 subroutine ZBL_ABINTIO_NN_Interaction(r,POTR, FPOTR)
 !***  PURPOSE: to calculate the potential and force by Born-Mayer model
 !     INPUT:   r, the distance between two particle
 !     OUTPUT:  POTR, FPOTR

 implicit none
     real(KINDDF), intent(in)::r
     real(KINDDF),intent(out)::POTR, FPOTR
     real(KINDDF)::RR,POTRMID,ZBLPOTEN,fc,dfc

      RR=R*1.0D8
      if((r.ge.0.0D0).and.(r.lt.2.2022866D-8)) then

          ZBLPOTEN=ZBLA/rr*dexp(-rr*ZBLB)
          POTR  =ZBLPOTEN+ZBLC*dexp(-rr*ZBLD)+ZBLE
          POTR  =0.5D0*CP_EVERG*POTR

          FPOTR =-ZBLA*dexp(-rr*ZBLB)*(1.0d0/rr/rr+ZBLB/RR)
          FPOTR=FPOTR-ZBLC*dexp(-rr*ZBLD)*ZBLD
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-FPOTR*CP_EVERG*1.0D8/R
           FPOTR =-FPOTR*CP_EVERG*1.0D8

     else if((r.ge.2.2022866D-8).and.(r.lt.2.7423008D-8))  then
          POTR  =EA*dexp(-rr*EB)+Ec
          POTR  =0.5D0*CP_EVERG*POTR
          FPOTR =-Ea*Eb*dexp(-rr*Eb)*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR/r
          FPOTR =-CP_EVERG*FPOTR

      else if((r.ge.2.7423008D-8).and.(r.lt.3.5D-8))  then
 ! --------correct
          POTR  =MA+MB*rr+MC*rr**2.0D0+MD*rr**3.0D0
          POTR  =CP_EVERG*0.5D0*POTR
          FPOTR =(MB+2.0D0*MC*rr+3.0D0*MD*rr**2.0D0)*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR/r
          FPOTR =-CP_EVERG*FPOTR
      else  if(r.ge.3.5D-8)  then
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
          POTRMID  =LA*((LB/rr)**LC-LD*(LB/rr)**LE)  !+LF
          POTR  =POTRMID*fc

          POTR  =CP_EVERG*0.5D0*POTR
          FPOTR =-LA*(LC*LB*(LB/rr)**(LC-1.0D0)-LD*LB*LE*(LB/rr)**(LE-1.0D0))/rr/rr
          FPOTR = FPOTR*FC+POTRMID*DFC
          !modified accordding to NOTE-2014-12-04 in forcetable:
          !FPOTR =-CP_EVERG*FPOTR*1.0D8/R
          FPOTR =-CP_EVERG*FPOTR*1.0D8
      end if

     return
 end subroutine ZBL_ABINTIO_NN_Interaction


 subroutine ZBL_ABINTIO_NE_Interaction(r,POTB, FPOTB)
 !***  PURPOSE: to calculate the attractive potential by Embbedment energy
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

 end module Pot_ZBL_ABINTIO_He_W_513
