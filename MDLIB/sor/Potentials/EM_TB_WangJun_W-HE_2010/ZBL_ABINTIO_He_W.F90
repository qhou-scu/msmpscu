  module Pot_ZBL_ABINTIO_He_W
  !***  this module is to calculate the potential and force between W-He
  !     atoms. Supplied by Wang Jun
  !REFERENCE:

  use MD_CONSTANTS
        implicit none


      !***  the parameter block for nuclear-nuclear part in EAM ***************
      real(KINDDF), parameter,private::IGCLONG   = 8.854187817D-12
      real(KINDDF), parameter,private::ELECHARGE  = 1.60219D-19
      real(KINDDF), parameter,private::z1=2.0D0
      real(KINDDF), parameter,private::z2=74.0D0
      real(KINDDF), parameter,private::a0=5.291772108D-11
      !INTEGER::I
      !real(KINDDF)::potentialv,a,r,faix,X,dfaix,dpotentialv



      real(KINDDF), parameter,private::MA=49.7192120535695D0   !7.92041989584602D0 6.43823829937402
      real(KINDDF), parameter,private::MB=-52.5892870966472D0    !-5.46574252078757D0
      real(KINDDF), parameter,private::MC=18.5897462728089D0    !1.56132350481559D0
      real(KINDDF), parameter,private::MD=-2.19112033728284D0     !-0.150061294708069D0


      real(KINDDF), parameter,private::EA=45.7916129377437D0         !6.717046865181104D89
      real(KINDDF), parameter,private::EB=0.4625661D0         !8.569512844085693D-3
      real(KINDDF), parameter,private::EC=-2.322303987550913D-2         !56.8629783391953D0
  !    real(KINDDF), parameter::ED=         !0.4467921D0
  !    real(KINDDF), parameter::EE=1.45839595913096977D0  ! -2.378563734103023E-002

      real(KINDDF), parameter,private::LA=0.007D0
      real(KINDDF), parameter,private::LB=4.914989D0
      real(KINDDF), parameter,private::LC=8.261011D0
      real(KINDDF), parameter,private::LD=1.670001D0
      real(KINDDF), parameter,private::LE=6.738495D0
      real(KINDDF), parameter,private::LF=1.482181596472D0

  contains

  subroutine ZBL_ABINTIO_NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by Nuclear-Nuclear interaction
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR

  implicit none
      real(KINDDF)::r  !,intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      real(KINDDF)::RR, X, A, faix, dfaix

       RR=R*1.0D8
       if((r.ge.0.0).and.(r.lt.1.8D-8)) then

           A=0.88534D0*a0/(Z1**0.23D0+Z2**0.23D0)
           X= R*1.0D-2/A

           faix=0.18175D0*DEXP(-3.1998D0*X)+0.50986D0*DEXP(-0.94229D0*X)+0.28022D0*DEXP(-0.4029D0*X)+0.028171D0*DEXP(-0.20162D0*X)
           POTR=ELECHARGE*Z1*Z2*faix/(4.0D0*CP_PI*IGCLONG*R*1.0D-2)
           POTR=0.5D0*POTR*CP_EVERG


           dfaix=0.18175D0*DEXP(-3.1998D0*X)*(-3.1998D0)+0.50986D0*DEXP(-0.94229D0*X)*(-0.94229D0)+0.28022D0*DEXP(-0.4029D0*X)*(-0.4029D0)+0.028171D0*DEXP(-0.20162D0*X)*(-0.20162D0)
           dfaix=dfaix/A

           FPOTR=ELECHARGE*Z1*Z2/(4.0D0*CP_PI*IGCLONG)*(dfaix*r-faix*1.0D2)/R/R
          !modified accordding to NOTE-2014-12-04 in forcetable:
           !FPOTR=-CP_EVERG*FPOTR/R
           FPOTR=-CP_EVERG*FPOTR

       else if((r.ge.1.8D-8).and.(r.lt.2.6D-8))  then
           POTR  =MA+MB*RR+MC*RR**2.0D0+MD*RR**3.0D0
           POTR  =CP_EVERG*0.5D0*POTR
           FPOTR =(MB+2.0D0*MC*RR+3.0D0*MD*rr**2.0D0)*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
           !FPOTR =-CP_EVERG*FPOTR/R
           FPOTR =-CP_EVERG*FPOTR
       else if((r.ge.2.6D-8).and.(r.lt.3.5D-8))    then
           POTR  =EA*dexp(-RR/EB)+Ec
           POTR  =0.5D0*CP_EVERG*POTR
           FPOTR =-Ea/Eb*dexp(-RR/Eb)*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
          ! FPOTR =-CP_EVERG*FPOTR/R
           FPOTR =-CP_EVERG*FPOTR

       else  if(r.ge.3.5D-8)  then

           POTR  =LA*((LB/RR)**LC-LD*(LB/RR)**LE)  !+LF
           POTR  =CP_EVERG*0.5D0*POTR
           FPOTR =-LA*(LC*LB*(LB/RR)**(LC-1.0D0)-LD*LB*LE*(LB/rr)**(LE-1.0D0))/RR/RR*1.0D8
          !modified accordding to NOTE-2014-12-04 in forcetable:
          ! FPOTR =-CP_EVERG*FPOTR/R
           FPOTR =-CP_EVERG*FPOTR
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

  end module Pot_ZBL_ABINTIO_He_W
