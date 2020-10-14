  module EAM2_WW_Marinica_JPCM25_2013
  !***  this module is to calculate the potential and force between W-W atoms,
  !     using EAM2 provided in reference
  !
  !     REFERENCE: M-C Marinica, Lisa Ventelon, M R Gilbert, L Proville, S L Dudarev,
  !                J Marian, G Bencteux and F Willaime, J. Phys. Condens. Metter 25(2013)395502
  !
  !     HISTORY:     adopted by Hou Qing, Dec 9, 2014
  !
  use MD_CONSTANTS
  implicit none


  contains

  !*****************************************************************
  subroutine NN_Func(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle in CM
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !***  the parameter block for nuclear-nuclear part in EAM ***************
       real(KINDDF),parameter::AFAI(15) = (/0.960851701343041D2, &
                                           -0.184410923895214D3, &
                                            0.935784079613550D2, &
                                           -0.798358265041677D1, &
                                            0.747034092936229D1, &
                                           -0.152756043708453D1, &
                                            0.125205932634393D1, &
                                            0.163082162159425D1, &
                                           -0.141854775352260D1, &
                                           -0.819936046256149D0, &
                                            0.198013514305908D1, &
                                           -0.696430179520267D0, &
                                            0.304546909722160D-1,&
                                           -0.163131143161660D1, &
                                            0.138409896486177D1/)

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
     !---
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
      !--- local varibales
      !***  the parameter block concerning RHO in EAM ***************
      real(KINDDF),parameter::ARHO(4) = (/-0.420429107805055D1,  &
                                           0.518217702261442D0,  &
                                           0.562720834534370D-1, &
                                           0.344164178842340D-1/)

      real(KINDDF),parameter::ARHOr(4) =(/2.500000000000000,     &
                                          3.100000000000000,     &
                                          3.500000000000000,     &
                                          4.900000000000000/)
      real(KINDDF),parameter::RC = 2.002970124727D0*CP_A2CM
      !-----

         if(r .le. RC) then
            call RHO_FuncPoly3(RC*CP_CM2A, POTB, FPOTB, ARHO, ARHOr)
            FPOTB = 0.D0
          else
            call RHO_FuncPoly3(r*CP_CM2A, POTB, FPOTB, ARHO, ARHOr)
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
      !---  the parameter block concerning Embedment function in EAM ***************
      real(KINDDF),parameter::AF(2) = (/-5.946454472402710D0,  &
                                        -0.049477376935239D0/)


      !******
        call EMBED_FS_PLOY_Func(rho,FRHO, DFRHO, AF)

      return
  end subroutine EMBED_Func
  !*****************************************************************

  end module  EAM2_WW_Marinica_JPCM25_2013
