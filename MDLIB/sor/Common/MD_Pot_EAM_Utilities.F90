  module MD_Pot_EAM_Utilities
  !***  this module contained a number of routines that are often used in creating
  !     EAM and FS potentials and froces.
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:     2014-12 (HOU Qing), created the first version
  !
  !                  2018-07-05(HOU Qing),
  !                          add NN_to_ZBL to connect the pairwise potential to ZBL potential
  !                              RHO_to_CONST1 to connect the RHO to a zero-grgient point
  !
  !
  use MD_CONSTANTS
  implicit none

        !**** the variables block ****************************************************
   private::H
   public::NN_FuncPoly3
   public::NN_FunExp1
   public::NN_FunExp2
   public::NN_Johnson
   public::NN_to_ZBL

   public::RHOZero_Func
   public::RHO_FuncPoly3
   public::RHO_FunExp1
   public::RHO_to_CONST1

   public::EMBEDZero_Func
   public::EMBED_FS_PLOY_Func
   
  contains
  !*****************************************************************
  ! the step function
  function H(X) result(RES)
  implicit none
  real(KINDDF), intent(in)::x
  real(KINDDF)::RES
     if (x.GE.0) then
         RES = 1.D0
     else
         RES = 0.D0
     end if
  end function H
  !*****************************************************************

  !*****************************************************************
  subroutine NN_FuncPoly3(RT,POTR, FPOTR, A, RI)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !              the V(r) is assumed in the form:
  !
  !              V(r) = A1*(R1-r)**3*H(R1-r) + A2*(R2-r)**3*H(R2-r)+...
  !
  !              where H(x) is step function.
  !
  !     INPUT:   RT, the distance between two particles
  !     OUTPUT:  POTR, FPOTR
  !
  !     NOTE:    the potetntial in literatures usually in eV and r is in A,
  !              we transfer the unit into EGS unit
  !
  !     MODIFY: dummy variable R changed to RT, RT in the same unit of RI, 1016/10/13/

  implicit none
      real(KINDDF),intent(in)::RT
      real(KINDDF),intent(out)::POTR, FPOTR
      real(KINDDF), dimension(:), intent(in)::A, RI
      !******
      real(KINDDF)::STP
      integer::I

        POTR  = 0.D0
        FPOTR = 0.D0
        !RT    = r*CP_CM2A
        do I=1, size(A)
           STP   = H(RI(I) - RT)
           POTR  = POTR  + A(I)*(RI(I) - RT)*(RI(I) - RT)*(RI(I) - RT)*STP
           FPOTR = FPOTR + A(I)*(RI(I) - RT)*(RI(I) - RT)*STP
        end do
        POTR   = 0.5D0*POTR*CP_EVERG
        !--- donnot forget multiply 3, the factor due to differential
        FPOTR  = 3.0D0*FPOTR*CP_EVERG/CP_A2CM

      return
  end subroutine NN_FuncPoly3
  !*****************************************************************

  !*****************************************************************
  subroutine NN_Johnson(RT, POTR, FPOTR, NNa, RHOa, NNb, RHOb)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !              for A-B atoms, where A, B denote different kind
  !              of elements.
  !              See: Johnson, PRB39(1989)12554:
  !
  !              Vab(r) = 0.5*[RHOb(r)*Vaa(r)/RHOa(r) + RHOa(r)*Vbb(r)/RHOb(r)]
  !
  !     INPUT:   RT,    the distance between two particles
  !              NNa,   the external function of pairwise potential for A-A interaction
  !              RHOa,  the external function of electron denstiy for atom of type A
  !              NNb,   the external function of pairwise potential for B-B interaction
  !              RHOb,  the external function of electron denstiy for atom of type B
  !     OUTPUT:  POTR, FPOTR
  !
  !     NOTE:    the potetntial in literatures usually in eV and r is in A,
  !              we transfer the unit into EGS unit
  !
  !     MODIFY: dummy variable R changed to RT, RT in the same unit of RI, 1016/10/13/

  implicit none
     real(KINDDF),intent(in)::RT
     real(KINDDF),intent(out)::POTR, FPOTR
     !---- the intereface
     interface
       subroutine NNa(R,POTR, FPOTR)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
         real(KINDDF),intent(in)::R
         real(KINDDF),intent(out)::POTR, FPOTR

        end subroutine NNa
     end interface

     interface
       subroutine RHOa(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
         real(KINDDF),intent(in)::R
         real(KINDDF),intent(out)::POTB, FPOTB

        end subroutine RHOa
     end interface

     interface
       subroutine NNb(R,POTR, FPOTR)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
         real(KINDDF),intent(in)::R
         real(KINDDF),intent(out)::POTR, FPOTR

        end subroutine NNb
     end interface

     interface
       subroutine RHOb(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
         real(KINDDF),intent(in)::R
         real(KINDDF),intent(out)::POTB, FPOTB

        end subroutine RHOb
     end interface
     !---- the end of intereface
      real(KINDDF)::A, B, DA, DB, WBSA, WASB, DWBSA, DWASB
      integer::I

        POTR  = 0.D0
        FPOTR = 0.D0

        call RHOa(RT,A,DA)
        call RHOb(RT,B,DB)
        WASB  = A/B
        WBSA  = B/A
        DWASB = (DA*B-A*DB)/B/B
        DWBSA = (DB*A-B*DA)/A/A
        call NNa(RT,A,DA)
        call NNb(RT,B,DB)
        POTR  = C_HALF*(WASB*B+WBSA*A)
        FPOTR = C_HALF*(WASB*DB+WBSA*DA)+(DWASB*B+DWBSA*A)
      return
  end subroutine NN_Johnson
  !*****************************************************************

  !*****************************************************************
  subroutine NN_FunExp1(RR,POTR, FPOTR, A, ALPHA, Re, Kapa, M) 
  !***  PURPOSE: to calculate the potential using exponential form
  !  
  !              V(r) = A*exp(-alpha*(r/re-1))/(1+(r/re-kapa)^m)
  !                   
  !     INPUT:   RR, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  !
  !     NOTE:    the potetntial in literatures usually in eV and r is in A,
  !              we transfer the unit into EGS unit
  !

     implicit none
     real(KINDDF),intent(in) ::RR, A, ALPHA, Re, Kapa, M
     real(KINDDF),intent(out)::POTR, FPOTR
       !****
       real(KINDDF)::R 
       
       R = RR !*CP_CM2A
       POTR  =A*Dexp(-Alpha*(R/Re-1.D0))/(1.D0+(R/Re-Kapa)**M)
       FPOTR =A*Dexp(-ALPHA*(R/Re-1.D0))*(-Alpha*(1.D0+(R/Re-Kapa)**M)-M*(R/Re-Kapa)**(M-1.D0))/ &
             (Re*(1.D0+(R/Re-Kapa)**M)**2.D0)
              
       POTR  =   C_HALF*POTR*CP_EVERG
       FPOTR  = -FPOTR*CP_EVERG/CP_A2CM
      
    return
  end subroutine NN_FunExp1  
  !*****************************************************************

  !******************************************************************
  subroutine NN_FunExp2(RR, POTR, FPOTR, A1, ALPHA1, Re1, Kapa1, M1, A2, ALPHA2, Re2, Kapa2, M2)
  !***  PURPOSE: to calculate the potential using exponential form
  !  
  !              V(r) = A1*exp(-alpha1*(r/re1-1))/(1+(r/re1-kapa1)^m1) 
  !                    +A2*exp(-alpha2*(r/re2-1))/(1+(r/re2-kapa2)^m2)          
  !                   
  !     INPUT:   RR, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  !
  !     NOTE:    the potetntial in literatures usually in eV and r is in A,
  !              we transfer the unit into EGS unit
  !

    implicit none
     real(KINDDF),    intent(in) ::RR, A1, ALPHA1, Re1, Kapa1, M1, A2, ALPHA2, Re2, Kapa2, M2
     real(KIND=8),intent(out)::POTR, FPOTR
     real(KINDDF)::POTR0, FPOTR0
     
        call  NN_FunExp1(RR,POTR0, FPOTR0, A1, ALPHA1, Re1, Kapa1, M1) 
        POTR =  POTR0
        FPOTR=  FPOTR0

        call  NN_FunExp1(RR,POTR0, FPOTR0, A2, ALPHA2, Re2, Kapa2, M2) 
        POTR = POTR -  POTR0
        FPOTR= FPOTR - FPOTR0
 
    return
  end subroutine NN_FunExp2
  !*****************************************************************
  
  !*****************************************************************
  subroutine RHOZero_Func(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the electron density RHO(r) and  -dRHO(r)/dr
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB
      !******
      real(KINDDF)::STP, RT
      integer::I

        POTB  = 0.D0
        FPOTB = 0.D0
      return
  end subroutine RHOZero_Func
  !*****************************************************************

  !*****************************************************************
  subroutine RHO_FuncPoly3(RT,POTB, FPOTB, A, RI)
  !***  PURPOSE: to calculate the electron density RHO(r) and  -dRHO(r)/dr
  !              the RHO(r) is assumed in the form:
  !
  !              RHO(r) = A1*(R1-r)**3*H(R1-r) + A2*(R2-R)**3*H(R2-r)+...
  !
  !              where H(x) is step function.
  !
  !     INPUT:   RT, the distance between two particle
  !     OUTPUT:
  !     NOTE:    the r in literatures usually is in A,
  !              we transfer the unit into EGS unit
  !     MODIFY: dummy variable R changed to RT, RT in the same unit of RI, 1016/10/13/
  implicit none
      real(KINDDF),intent(in)::RT
      real(KINDDF),intent(out)::POTB, FPOTB
      real(KINDDF), dimension(:), intent(in)::A, RI
      !******
      real(KINDDF)::STP
      integer::I

        POTB  = 0.D0
        FPOTB = 0.D0
        !RT    = r*CP_CM2A
        do I=1, size(A)
           STP   = H(RI(I) - RT)
           POTB  = POTB  + A(I)*(RI(I) - RT)*(RI(I) - RT)*(RI(I) - RT)*STP
           FPOTB = FPOTB + A(I)*(RI(I) - RT)*(RI(I) - RT)*STP
        end do
        FPOTB  = 3.D0*FPOTB/CP_A2CM
      return
  end subroutine RHO_FuncPoly3
  !*****************************************************************

  !*****************************************************************
  subroutine RHO_FunExp1(R,POTB, FPOTB, A, ALPHA, Re, Kapa, M) 
  !***  PURPOSE: to calculate the potential using exponential form
  !  
  !              RHO(r) = A*exp(-alpha*(r/re-1))/(1+(r/re-kapa)^m)
  !                   
  !     INPUT:   RR, the distance between two particle
  !     OUTPUT:  POTB, FPOTB
     implicit none
     real(KINDDF),intent(in) ::R, A, ALPHA, Re, Kapa, M
     real(KINDDF),intent(out)::POTB, FPOTB
        !****        
         POTB  =A*Dexp(-Alpha*(R/Re-1.D0))/(1.D0+(R/Re-Kapa)**M)
         FPOTB =A*Dexp(-ALPHA*(R/Re-1.D0))*(-Alpha*(1.D0+(R/Re-Kapa)**M)-M*(R/Re-Kapa)**(M-1.D0))/ &
               (Re*(1.D0+(R/Re-Kapa)**M)**2.D0)
      
    return
  end subroutine RHO_FunExp1  
  !*****************************************************************


  !*****************************************************************
  subroutine EMBEDZero_Func(rho,FRHO, DFRHO)
  !***  PURPOSE: to calculate the embedment function F(rho) and  dF(rho)/drho
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
  implicit none
      real(KINDDF),intent(in)::rho
      real(KINDDF),intent(out)::FRHO, DFRHO

          FRHO  = 0.D0
          DFRHO = 0.D0

      return
  end subroutine EMBEDZero_Func
  !*****************************************************************

  !*****************************************************************
  subroutine EMBED_FS_PLOY_Func(rho,FRHO, DFRHO,A)
  !***  PURPOSE: to calculate the embedment function F(rho) and  dF(rho)/drho
  !              the embedment function is assumed
  !
  !              F(rho) = A1*sqrt(rho) + A2*rho^2 + A3*rho^3 + ...
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
  !
  !     NOTE:    the potetntial in literatures usually in eV, we
  !              transfer the unit into EGS unit

  implicit none
      real(KINDDF),intent(in)::rho
      real(KINDDF),intent(out)::FRHO, DFRHO
      real(KINDDF), dimension(:), intent(in)::A
      integer::I
      real(KINDDF)::TRHO

      !******
        if(rho .le. 0.D0) then
           FRHO  = 0.D0
           DFRHO = 0.D0
        else
           TRHO = rho
           FRHO  = A(1)*dsqrt(TRHO)
           DFRHO = C_HALF*A(1)/dsqrt(TRHO)
           do I=2, size(A)
              FRHO  = FRHO  + A(I)*rho*TRHO
              DFRHO = DFRHO + dble(I)*A(I)*TRHO
              TRHO  = rho*TRHO
           end do
        end if
        FRHO  = FRHO*CP_EVERG
        DFRHO = DFRHO*CP_EVERG
      return
  end subroutine EMBED_FS_PLOY_Func
  !*****************************************************************

  !*****************************************************************
  subroutine NN_to_ZBL(RT, POTR, FPOTR, Z1, Z2, NN0, R0, R1)
  !***  PURPOSE: to conect the pairwizw potential to ZBL potential
  !
  !              V(r) = lamda*Vzbl(r) + (1-lamda)*NN0(r)
  !              where NN0(r) is the orignal potential to be connected to ZBL
  !              lamda is a function of r:
  !              lamda(r) = C1+C2*r+C3*r^2+C4*r^3
  !
  !     INPUT:   RT,     the distance between two particles
  !              NN0,    the external function of pairwise potential
  !              Z1, Z2, the atomic number of atoms to be used for ZBL potential
  !              R0, R1, the endpoints at which the function is connected
  !     OUTPUT:  POTR, FPOTR
  !
  use MD_CorePotent
  implicit none
     real(KINDDF),intent(in)::RT, Z1, Z2, R0, R1
     real(KINDDF),intent(out)::POTR, FPOTR
     !---- the intereface
     interface
       subroutine NN0(R,POTR, FPOTR)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
         real(KINDDF),intent(in)::R
         real(KINDDF),intent(out)::POTR, FPOTR

        end subroutine NN0
     end interface

     !---- the end of intereface
      real(KINDDF)::VZBL, VNN0, FZBL, FNN0, LAM, DLAM, CC(4)

        call ZBL_Pot(Z1, Z2, RT, VZBL, FZBL)
        call NN0(RT, VNN0, FNN0)

        CC(1) =  R0*R0*(R0 - 3.D0*R1)
        CC(2) =  6.D0*R0*R1
        CC(3) = -3.D0*(R0 + R1)
        CC(4) =  2.D0
        CC    =  CC/(R0 -  R1)**3.D0
        LAM   = CC(1) + CC(2)*RT + CC(3)*RT*RT    + CC(4)*RT*RT*RT
        DLAM  =         CC(2)    + CC(3)*RT*2.D0  + CC(4)*RT*RT*3.D0

        if(RT .lt. R0) then
           LAM  = 0.D0
           DLAM = 0.D0
        else if(RT .gt. R1) then
           LAM =  1.D0   
           DLAM = 0.D0
        end if   
        POTR  = LAM*VNN0 + (1-LAM)*VZBL
        FPOTR = DLAM*VNN0 + LAM*FNN0 - DLAM*VZBL + (1-LAM)*FZBL
      return
  end subroutine NN_to_ZBL
  !*****************************************************************

  !*****************************************************************
  subroutine RHO_to_CONST1(RT, POTB, FPOTB, RHO0, R0, R1)
  !***  PURPOSE: to conect the RHO smoothly to a constant
  !
  !              RHO(r) = F(r)*RHO0(r)
  !              where F(r) is a function satisfying:
  !              RHO(R1)      = RHO0(R1)
  !              dRHO(R1)/dR1 = dRHO0(R1)/dR1
  !              dRHO(R0)/dR0 = 0
  !
  !
  !     INPUT:   RT,     the distance between two particles
  !              RHO0,   the external function of RHO
  !
  !     OUTPUT:  POTB, FPOTB
  !
  use MD_CorePotent
  implicit none
     real(KINDDF),intent(in) ::RT, R0, R1
     real(KINDDF),intent(out)::POTB, FPOTB
     !---- the intereface
     interface
       subroutine RHO0(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
         real(KINDDF),intent(in)::R
         real(KINDDF),intent(out)::POTB, FPOTB

        end subroutine RHO0
     end interface

     !---- the end of intereface
      real(KINDDF)::POTB0, FPOTB0, A, RR

       !--- to calculate alpha by end point 
        call RHO0(R0, POTB0, FPOTB0)
        A  = - FPOTB0/(2.D0*POTB0/(R0-R1) + FPOTB0)
        RR =  (R0 - R1)*(R0 - R1)

        !--- to calculate the rho and frho at RT
        if(RT .lt. R0) then
          POTB  = (A + 1.D0)*POTB0
          FPOTB = 0.D0
        else 
          call RHO0(RT, POTB0, FPOTB0)  
          if(RT .lt. R1) then 
             POTB  = (A*(RT - R1)*(RT - R1)/RR + 1.D0)*POTB0
             FPOTB = (A*(RT - R1)*2.D0/RR)*POTB0 + (A*(RT - R1)*(RT - R1)/RR + 1.D0)*FPOTB0
          else    
            POTB  = POTB0
            FPOTB = FPOTB0
          end if 
        end if      
      return
  end subroutine RHO_to_CONST1
  !*****************************************************************


  end module  MD_Pot_EAM_Utilities
