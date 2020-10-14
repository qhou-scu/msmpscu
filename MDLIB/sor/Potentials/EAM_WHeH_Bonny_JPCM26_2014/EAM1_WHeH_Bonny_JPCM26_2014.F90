  module EAM1_WHeH_Bonny_JPCM26_2014
  !***  this module is to calculate the potential and force in WHeH system
  !
  !     REFERENCE:   G Bonny, P Grigorev and D Terentyev, JPCM26(2014)485001
  !
  !     HISTORY:     adopted by Hou Qing, Dec 9, 2014
  !
  use MD_CONSTANTS
        implicit none

        !**** the variables block ****************************************************
  contains

  !*****************************************************************
  subroutine NN_WW_Func(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use EAM2_WW_Marinica_JPCM25_2013, only:NN_Marinica_Func=>NN_Func, &
                                         RHO_Marinica_Func=>RHO_Func
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF)::RHO, FRHO
      real(KINDDF), parameter::CP_C = 1.848055990D0*CP_EVERG

         call NN_Marinica_Func(r,POTR, FPOTR)

        !$$--- do gauge transformation
         call RHO_Marinica_Func(r, RHO, FRHO)
        !note: POTR has been multiplied C_HAlF in NN_Marinica_Func
        !      see also the definition of POTR in MD_TYPEDEF_ForceTable.F90
         POTR  = POTR - CP_C*RHO
         FPOTR = FPOTR - 2.D0*CP_C*FRHO

      return
  end subroutine NN_WW_Func
  !*****************************************************************

  !*****************************************************************
  subroutine RHO_WW_Func(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the electron density RHO(r) and  -dRHO(r)/dr
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  use EAM2_WW_Marinica_JPCM25_2013, only:RHO_Marinica_Func=>RHO_Func
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB
      !******
      real(KINDDF), parameter::CP_S = 2.232322602D-1

        call RHO_Marinica_Func(r,POTB, FPOTB)

        !$$--- do gauge transformation
        POTB  = POTB*CP_S
        FPOTB = FPOTB*CP_S
      return
  end subroutine RHO_WW_Func
  !*****************************************************************

  !*****************************************************************
  subroutine EMBED_WW_Func(rho,FRHO, DFRHO)
  !***  PURPOSE: to calculate the embedment function F(rho) and  dF(rho)/drho
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
  use EAM2_WW_Marinica_JPCM25_2013, only:EMBED_Marinica_Func=>EMBED_Func
  implicit none
      real(KINDDF),intent(in)::rho
      real(KINDDF),intent(out)::FRHO, DFRHO

      !******
      real(KINDDF), parameter::CP_S    =  2.232322602D-01 !*CP_EVERG
      real(KINDDF), parameter::CP_IS   =  1.D+00/CP_S

      real(KINDDF), parameter::CP_C    =  1.848055990D+00*CP_EVERG
      real(KINDDF), parameter::CP_COS  =  CP_C*CP_IS
      real(KINDDF), parameter::CP_RHOI =  1.359141225D0 !*CP_EVERG
      real(KINDDF), parameter::CP_A0   = -5.524855802D+00, &
                               CP_A1   =  2.317313103D-01, &
                               CP_A2   = -3.665345949D-02, &
                               CP_A3   =  8.989367404D-03


      real(KINDDF)::TRHO, TFRHO, TDFRHO
      integer::I

       if(rho .le. CP_RHOI) then
          TRHO = rho*CP_IS
          call EMBED_Marinica_Func(TRHO,TFRHO, TDFRHO)
          !$$--- do gauge transformation
           FRHO  = TFRHO  + CP_C*TRHO
           DFRHO = TDFRHO*CP_IS + CP_COS
       else
           TRHO  = rho !/CP_EVERG
           FRHO  = (CP_A0 + TRHO*(CP_A1 + TRHO*(CP_A2 + TRHO*CP_A3)) )*CP_EVERG
           DFRHO = (CP_A1 + TRHO*(2.D0*CP_A2 + 3.D0*CP_A3*TRHO) )*CP_EVERG
       end if

      return
  end subroutine EMBED_WW_Func
  !*****************************************************************

  !*****************************************************************
  subroutine NN_HeHe_Func1(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(2), parameter::RI = (/2.D+00, 3.D+00/)
      real(KINDDF), dimension(2), parameter::A =  (/2.106615791D+00, -2.217639348D-01/)

        call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_HeHe_Func1
  !*****************************************************************

  !*****************************************************************
  subroutine NN_HeHe_Func2(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(4), parameter::RI = (/1.800000000D+00, &
                                                    2.000000000D+00, &
                                                    2.100000000D+00, &
                                                    3.000000000D+00/)
      real(KINDDF), dimension(4), parameter::A =  (/2.000000000D+01, &
                                                    1.051327582D+00, &
                                                   -4.000000000D+00, &
                                                    3.203322671D-02/)

        call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_HeHe_Func2
  !*****************************************************************

  !*****************************************************************
  subroutine NN_HH_Func1(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(2), parameter::RI = (/2.D+00, 3.D+00/)
      real(KINDDF), dimension(2), parameter::A =  (/4.862785907D-01, 1.018797872D-01/)

          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_HH_Func1
  !*****************************************************************

  !*****************************************************************
  subroutine NN_HH_Func2(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(4), parameter::RI = (/1.600000000D+00, &
                                                    2.000000000D+00, &
                                                    2.200000000D+00, &
                                                    3.000000000D+00/)

      real(KINDDF), dimension(4), parameter::A =  (/4.000000000D+01, &
                                                    2.315124670D-01, &
                                                   -2.000000000D-01, &
                                                    5.180584543D-02/)

          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_HH_Func2
  !*****************************************************************

  !*****************************************************************
  subroutine EMBED_HH_Func2(rho,FRHO, DFRHO)
  !***  PURPOSE: to calculate the embedment function F(rho) and  dF(rho)/drho
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
  use  MD_Pot_EAM_Utilities, only:EMBED_FS_PLOY_Func
  implicit none
      real(KINDDF),intent(in)::rho
      real(KINDDF),intent(out)::FRHO, DFRHO

      !******
      real(KINDDF), parameter::A(2) = (/-2.610066441D+01, 4.688963869D-01/)
      real(KINDDF)::TRHO

       TRHO = rho
       call  EMBED_FS_PLOY_Func(TRHO,FRHO, DFRHO,A)
       return

       !if(rho .gt. 0.D0) then
       !    FRHO  = (A(1)*dsqrt(TRHO) +A(2)*TRHO*TRHO )
       !    DFRHO = (C_HALF*A(1)/dsqrt(TRHO) + 2.D0*A(2)*TRHO)
       !else
       !    FRHO = 0.D0
       !    DFRHO = 0.D0
       !end if
       !FRHO  = FRHO  * CP_EVERG
       !DFRHO = DFRHO * CP_EVERG
       !FRHO  = 0
       !DFRHO = 0

      return
  end subroutine EMBED_HH_Func2
  !*****************************************************************


  !*****************************************************************
  subroutine NN_HHe_Func1(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(3), parameter::RI = (/1.8D+00, 2.0D+00, 3.0D+00/)
      real(KINDDF), dimension(3), parameter::A =  (/1.500000000D+01, 2.563700119D-01,-4.489510592D-02/)

          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_HHe_Func1
  !*****************************************************************

  !*****************************************************************
  subroutine NN_HHe_Func2(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(3), parameter::RI = (/2.000000000D+00, &
                                                    2.500000000D+00, &
                                                    3.000000000D+00/)
      real(KINDDF), dimension(3), parameter::A =  (/3.256370012D+00, &
                                                   -4.000000000D-01, &
                                                   -4.489510592D-02/)


          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_HHe_Func2
  !*****************************************************************

  !*****************************************************************
  subroutine NN_WHe_Func1(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(3), parameter::RI = (/1.9D+00, 2.2D+00, 3.5D+00/)
      real(KINDDF), dimension(3), parameter::A =  (/2.100000000D+01, &
                                                    8.565323293D-01, &
                                                    2.750099819D-01/)

          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_WHe_Func1
  !*****************************************************************

  !*****************************************************************
  subroutine NN_WHe_Func2(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(4), parameter::RI = (/1.900000000D+00, &
                                                    2.000000000D+00, &
                                                    2.200000000D+00, &
                                                    3.500000000D+00/)

      real(KINDDF), dimension(4), parameter::A =  (/0.000000000D+00, &
                                                    1.400000000D+01, &
                                                   -3.712116187D+00, &
                                                    3.105031456D-01/)


          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_WHe_Func2
  !*****************************************************************


  !*****************************************************************
  subroutine NN_WH_Func1(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(2), parameter::RI = (/2.0D+00, 3.0D+00/)
      real(KINDDF), dimension(2), parameter::A =  (/1.375733214D+01, 1.296071475D-01/)

          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_WH_Func1
  !*****************************************************************

  !*****************************************************************
  subroutine NN_WH_Func2(r,POTR, FPOTR)
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTR, FPOTR
      !--- local
      real(KINDDF), dimension(3), parameter::RI = (/2.000000000D+00, &
                                                    2.647500000D+00, &
                                                    3.295000000D+00/)

      real(KINDDF), dimension(3), parameter::A =  (/4.424459079D+01, &
                                                   -4.993477782D+00, &
                                                    1.461712984D+00/)


          call NN_FuncPoly3(r*CP_CM2A,POTR, FPOTR, A, RI)

      return
  end subroutine NN_WH_Func2
  !*****************************************************************

  !*****************************************************************
  subroutine RHO_WH_Func(r,POTB, FPOTB)
  !***  PURPOSE: to calculate the electron density RHO(r) and  -dRHO(r)/dr
  !
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:
  use EAM2_WW_Marinica_JPCM25_2013, only:RHO_Marinica_Func=>RHO_Func
  implicit none
      real(KINDDF),intent(in)::r
      real(KINDDF),intent(out)::POTB, FPOTB
      !******

         call RHO_WW_Func(r,POTB, FPOTB)
      return
  end subroutine RHO_WH_Func
  !*****************************************************************


  end module  EAM1_WHeH_Bonny_JPCM26_2014
