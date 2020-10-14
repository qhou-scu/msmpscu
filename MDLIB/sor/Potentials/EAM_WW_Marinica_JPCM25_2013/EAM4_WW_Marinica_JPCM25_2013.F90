  module EAM4_WW_Marinica_JPCM25_2013
  !***  this module is to calculate the potential and force between W-W atoms,
  !     using EAM2 provided in reference
  !
  !     REFERENCE: M-C Marinica, Lisa Ventelon, M R Gilbert, L Proville, S L Dudarev,
  !                J Marian, G Bencteux and F Willaime, J. Phys. Condens. Metter 25(2013)395502
  !
  !     HISTORY:     adopted by Hou Qing, Dec 9, 2014
  !                  table created by Li Min, Dec 12, 2014
  !
  use MD_CONSTANTS
  implicit none

  contains
  !*****************************************************************
  subroutine NN_Func(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !***  the parameter block for nuclear-nuclear part in EAM ***************
      real(KINDDF),parameter::AFAI(15) = (/0.954071477542914D2, &
                                          -0.181161004448916D3, &
                                           0.930215233132627D2, &
                                          -0.108428907515182D2, &
                                           0.112027468539573D2, &
                                          -0.312459176640227D1, &
                                           0.123028140617302D1, &
                                           0.154767467307454D1, &
                                          -0.128861387780439D1, &
                                          -0.843327493551467D0, &
                                           0.214009882965042D1, &
                                          -0.102898314532388D1, &
                                           0.138163259601912D1, &
                                          -0.360872433001551D1, &
                                           0.217655968740690D1/)

        real(KINDDF),parameter::AFAIr(15) = (/2.564897500000000,  &
                                              2.629795000000000,  &
                                              2.694692500000000,  &
                                              2.866317500000000,  &
                                              2.973045000000000,  &
                                              3.079772500000000,  &
                                              3.516472500000000,  &
                                              3.846445000000000,  &
                                              4.176417500000000,  &
                                              4.700845000000000,  &
                                              4.895300000000000,  &
                                              5.089755000000000,  &
                                              5.342952500000000,  &
                                              5.401695000000000,  &
                                              5.460437500000000/)
      !*****
         call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, AFAI, AFAIr)
      return
  end subroutine NN_Func
  !*****************************************************************

  !*****************************************************************
  subroutine RHO_Func(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the electron density RHO(r) and  -dRHO(r)/dr
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  use MD_Pot_EAM_Utilities, only:RHO_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB
      !***  the parameter block concerning RHO in EAM ***************
      real(KINDDF),parameter::ARHO(4) = (/-0.420429107805055D-1,  &
                                           0.518217702261442D0,   &
                                           0.562720834534370D-1,  &
                                           0.344164178842340D-1/)

      real(KINDDF),parameter::ARHOr(4) =(/2.500000000000000D0,    &
                                          3.100000000000000D0,    &
                                          3.500000000000000D0,    &
                                          4.900000000000000D0/)

      real(KINDDF),parameter::RC = 2.002970124727D0*CP_A2CM
      !****

         if(r .le. RC) then
            call RHO_FuncPoly3(RC*CP_CM2A,POTB, FPOTB, ARHO, ARHOr)
            FPOTB = 0.D0
          else
            call RHO_FuncPoly3(r*CP_CM2A,POTB, FPOTB, ARHO, ARHOr)
          end if

      return
  end subroutine RHO_Func
  !*****************************************************************

  !*****************************************************************
  subroutine EMBED_Func(rho,FRHO, DFRHO)
  !***  PURPOSE: to calculate the embedment function F(rho) and  dF(rho)/drho
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
  use MD_Pot_EAM_Utilities, only:EMBED_FS_PLOY_Func
  implicit none
      real(KINDDF),intent(in)::rho
      real(KINDDF),intent(out)::FRHO, DFRHO
      !***  the parameter block concerning Embedment function in EAM ***************
        real(KINDDF),parameter::AF(2) = (/-5.553986589859130D0, -0.045691157657292D0/)

      !******
        call EMBED_FS_PLOY_Func(rho,FRHO, DFRHO, AF)

      return
  end subroutine EMBED_Func
  !*****************************************************************

  end module  EAM4_WW_Marinica_JPCM25_2013
