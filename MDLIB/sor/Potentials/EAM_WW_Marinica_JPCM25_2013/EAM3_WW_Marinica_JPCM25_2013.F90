  module EAM3_WW_Marinica_JPCM25_2013
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
       real(KINDDF),parameter::AFAI(14) = (/0.335180090207171D3,  &
                                           -0.630187843239795D3,  &
                                            0.351436490059390D3,  &
                                           -0.187104554936536D-1, &
                                            0.599698179163370D-1, &
                                           -0.210635074070024D-1, &
                                            0.258575483268791D-1, &
                                            0.389896994857220D-1, &
                                           -0.139737073858817D-1, &
                                            0.184785876328131D-1, &
                                           -0.505971175928021D0,  &
                                           -0.890738292956047D1,  &
                                            0.495992901621487D-1, &
                                            0.377639491131043D-1/)

        real(KINDDF),parameter::AFAIr(14) = (/2.610000000000000,  &
                                              2.680000000000000,  &
                                              2.718160253915930,  &
                                              2.758522910939240,  &
                                              2.892427681390960,  &
                                              3.401517108884100,  &
                                              3.675283167162630,  &
                                              4.248392613157640,  &
                                              4.631029170086540,  &
                                              4.822645095086770,  &
                                              5.021214835934140,  &
                                              5.282885848627100,  &
                                              5.283029247305360,  &
                                              5.290925630219230/)

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
      real(KINDDF),parameter::ARHO(4) = (/-0.420429107805055D-1, &
                                           0.518217702261442D0,  &
                                           0.562720834534370D-1, &
                                           0.344164178842340D-1/)

      real(KINDDF),parameter::ARHOr(4) =(/2.500000000000000,     &
                                          3.100000000000000,     &
                                          3.500000000000000,     &
                                          4.900000000000000/)

      real(KINDDF),parameter::RC = 2.002970124727D0*CP_A2CM
      !***
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
      real(KINDDF),parameter::AF(2) = (/-5.073237905708500, -0.039263610952181/)
      !******

        call EMBED_FS_PLOY_Func(rho,FRHO, DFRHO, AF)

      return
  end subroutine EMBED_Func
  !*****************************************************************

  end module  EAM3_WW_Marinica_JPCM25_2013
