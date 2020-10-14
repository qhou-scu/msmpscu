  module Pot_FS_ACK_WW
  !***  this module is to calculate the potential and force between W-W atoms
  !REFERENCE: G. J. Ackland and R. Thetford, Philosophical Magazine A Vol.56(1987)15-30.
  !
  !           Originally adopted by Wang Jun, and modified by Hou Qing

  use MD_CONSTANTS
        implicit none


        real(KINDDF),parameter,private::A= 1.896373D8*CP_EVERG
        real(KINDDF),parameter,private::lattidis=3.1652D-8
        real(KINDDF),parameter,private::c=3.25D-8
        real(KINDDF),parameter,private::c0=47.1346499D16*CP_EVERG
        real(KINDDF),parameter,private::c1=-33.7665655D24*CP_EVERG
        real(KINDDF),parameter,private::c2=6.2541999D32*CP_EVERG
        real(KINDDF),parameter,private::d=4.400224D-8
        real(KINDDF),parameter,private::ACKB=90.3D24*CP_EVERG
        real(KINDDF),parameter,private::ACKA=1.2D8
        real(KINDDF),parameter,private::ACKB0=2.7411D-8
  contains

  subroutine TB_NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by Born-Mayer model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      real(KINDDF)::ACKFS,DACKFS
      !******

       IF(R.LT.ACKB0) THEN
            ACKFS=ACKB*(ACKB0-R)**3.0D0*DEXP(-ACKA*R)
            DACKFS=-ACKB*(3.0D0*(ACKB0-R)**2.0D0*DEXP(-ACKA*R)+ACKA*(ACKB0-R)**3.0D0*DEXP(-ACKA*R))
       ELSE
           ACKFS=0.0D0
           DACKFS=0.0D0
       END IF


      if (R.LE.C) then
          POTR=(R-C)**2.0D0*(C0+C1*R+C2*R**2.0D0)+ACKFS
          POTR=0.5D0*POTR

          FPOTR= 2.0D0*(R-C)*(C0+C1*R+C2*R**2.0D0)+(R-C)**2.0D0*(C1+2.0D0*C2*R)+DACKFS
          !modified accordding to NOTE-2014-12-04 in forcetable
          !FPOTR=-FPOTR/R
          FPOTR=-FPOTR
      else
          POTR=0.0D0
          FPOTR=0.0D0
      end if


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

      if (R.LE.D) then
          POTB=A*A*(R-D)**2.0D0

          FPOTB=A*A*2.0D0*(R-D)

          !--- NOTE: the force term is multiplied by 0.5
          !          this is from the differential of Finis-Sinclar form
          !          SQRT(RHO).
          !          Be careful in the later calculation of forces.
          !--- old version
          !FPOTB=-0.5D0*FPOTB/R

          !modified accordding to NOTE-2014-12-04 in forcetable
          FPOTB= -FPOTB
      else
          POTB=0.0D0
          FPOTB=0.0D0
      end if


      return
  end subroutine TB_NE_Interaction

  end module  Pot_FS_ACK_WW
