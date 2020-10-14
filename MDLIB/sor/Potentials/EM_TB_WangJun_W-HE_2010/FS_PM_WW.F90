  module Pot_FS_PM_WW
  !***  this module is to calculate the potential and force between W-W atoms
  !REFERENCE: P.M.Derlet, D.Nguyen-Manh and L.Dudarev, Phys.Rev.B76(2007) 054107
  !
  !           Originally adopted by Wang Jun, and modified by Hou Qing

  use MD_CONSTANTS
        implicit none


        !***  the parameter block for nuclear-nuclear part in EAM ***************
        real(KINDDF),parameter,private::A= 10.84238200368439*CP_EVERG

        real(KINDDF),parameter,private::rf1=3.3216D-8,                &
                                        rf2=3.165333333333333D-8,     &
                                        rf3=3.009066666666667D-8,     &
                                        rf4=2.8528D-8,                &
                                        rf5=2.7411D-8,                &
                                        rf6=2.604045D-8,              &
                                        rf7=2.46699D-8

        real(KINDDF),parameter,private::f1=1.677334871575606D24,      &
                                        f2=-1.673137509507131D24,     &
                                        f3=-4.737051994974169D24,     &
                                        f4=7.870789420780674D24,      &
                                        f5=-3.7924230525563378D20,    &
                                        f6=6.190816158916646D24,      &
                                        f7=-0.9565713891199151D24

        real(KINDDF),parameter,private::rv1=4.2689D-8,                &
                                        rv2=3.98568D-8,               &
                                        rv3=3.70246D-8,               &
                                        rv4=3.41924D-8,               &
                                        rv5=3.13602D-8,               &
                                        rv6=2.8528D-8,                &
                                        rv7=2.74110D-8,               &
                                        rv8=2.604045D-8,              &
                                        rv9=2.46699D-8

        real(KINDDF),parameter,private::v1=-0.1036435865158945D24*CP_EVERG, &
                                        v2=-0.2912948318493851D24*CP_EVERG, &
                                        v3=-2.096765499656263D24*CP_EVERG,  &
                                        v4=19.1604545270101D24*CP_EVERG,    &
                                        v5=-41.01619862085917D24*CP_EVERG,  &
                                        v6=46.05205617244703D24*CP_EVERG,   &
                                        v7=26.42203930654883D24*CP_EVERG,   &
                                        v8=15.35211507804088D24*CP_EVERG,   &
                                        v9=14.12806259323987D24*CP_EVERG


  contains

  subroutine TB_NN_Interaction(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by Born-Mayer model
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR

  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      real(KINDDF)::theda1,theda2,theda3,theda4,theda5,theda6,theda7,theda8,theda9
      !******
      if (r.GT.rv1) then
          theda1=0.0D0
      else
          theda1=1.0D0
      end if

      if (r.GT.rv2) then
          theda2=0.0D0
      else
          theda2=1.0D0
      end if

      if (r.GT.rv3) then
          theda3=0.0D0
      else
          theda3=1.0D0
      end if

      if (r.GT.rv4) then
          theda4=0.0D0
      else
          theda4=1.0D0
      end if

      if (r.GT.rv5) then
          theda5=0.0D0
      else
          theda5=1.0D0
      end if

      if (r.GT.rv6) then
          theda6=0.0D0
      else
          theda6=1.0D0
      end if

      if (r.GT.rv7) then
          theda7=0.0D0
      else
          theda7=1.0D0
      end if

      if (r.GT.rv8) then
          theda8=0.0D0
      else
          theda8=1.0D0
      end if

      if (r.GT.rv9) then
          theda9=0.0D0
      else
          theda9=1.0D0
      end if

      POTR= v1*(rv1-r)*(rv1-r)*(rv1-r)*theda1+v2*(rv2-r)*(rv2-r)*(rv2-r)*theda2+v3*(rv3-r)*(rv3-r)*(rv3-r)*theda3+v4*(rv4-r)*(rv4-r)*(rv4-r)*theda4 +&
            v5*(rv5-r)*(rv5-r)*(rv5-r)*theda5+v6*(rv6-r)*(rv6-r)*(rv6-r)*theda6+v7*(rv7-r)*(rv7-r)*(rv7-r)*theda7+v8*(rv8-r)*(rv8-r)*(rv8-r)*theda8  +&
            v9*(rv9-r)*(rv9-r)*(rv9-r)*theda9

      POTR=0.5D0*POTR

      FPOTR= 3.0D0*V1*(rv1-r)*(rv1-r)*theda1+3.0D0*V2*(rv2-r)*(rv2-r)*theda2+3.0D0*V3*(rv3-r)*(rv3-r)*theda3+3.0D0*V4*(rv4-r)*(rv4-r)*theda4+  &
             3.0D0*V5*(rv5-r)*(rv5-r)*theda5+3.0D0*V6*(rv6-r)*(rv6-r)*theda6+3.0D0*V7*(rv7-r)*(rv7-r)*theda7+3.0D0*V8*(rv8-r)*(rv8-r)*theda8 +&
             3.0D0*V9*(rv9-r)*(rv9-r)*theda9

     !modified accordding to NOTE-2014-12-04 in forcetable
     !FPOTR=FPOTR/R

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
      real(KINDDF)::theda1,theda2,theda3,theda4,theda5,theda6,theda7


      if (r.GT.rF1) then
          theda1=0.0D0
      else
          theda1=1.0D0
      end if

      if (r.GT.rF2) then
          theda2=0.0D0
      else
          theda2=1.0D0
      end if

      if (r.GT.rF3) then
          theda3=0.0D0
      else
          theda3=1.0D0
      end if

      if (r.GT.rF4) then
          theda4=0.0D0
      else
          theda4=1.0D0
      end if

      if (r.GT.rF5) then
          theda5=0.0D0
      else
          theda5=1.0D0
      end if

      if (r.GT.rF6) then
          theda6=0.0D0
      else
          theda6=1.0D0
      end if

      if (r.GT.rF7) then
          theda7=0.0D0
      else
          theda7=1.0D0
      end if


      POTB=f1*(rf1-r)*(rf1-r)*(rf1-r)*theda1+f2*(rf2-r)*(rf2-r)*(rf2-r)*theda2+f3*(rf3-r)*(rf3-r)*(rf3-r)*theda3+f4*(rf4-r)*(rf4-r)*(rf4-r)*theda4 +&
           f5*(rf5-r)*(rf5-r)*(rf5-r)*theda5+f6*(rf6-r)*(rf6-r)*(rf6-r)*theda6+f7*(rf7-r)*(rf7-r)*(rf7-r)*theda7

      POTB=POTB*A*A

      FPOTB=3.0D0*f1*(rf1-r)*(rf1-r)*theda1+3.0D0*f2*(rf2-r)*(rf2-r)*theda2+3.0D0*f3*(rf3-r)*(rf3-r)*theda3+3.0D0*f4*(rf4-r)*(rf4-r)*theda4+  &
            3.0D0*f5*(rf5-r)*(rf5-r)*theda5+3.0D0*f6*(rf6-r)*(rf6-r)*theda6+3.0D0*f7*(rf7-r)*(rf7-r)*theda7

     !--- NOTE: the force term is multiplied by 0.5
     !          this is from the differential of Finis-Sinclar form
     !          SQRT(RHO).
     !          Be careful in the later calculation of forces.
     ! FPOTB=0.5D0*FPOTB*A*A/R

     !modified accordding to NOTE-2014-12-04 in forcetable
      FPOTB = FPOTB*A*A

      return
  end subroutine TB_NE_Interaction

  end module  Pot_FS_PM_WW
