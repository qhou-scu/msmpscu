  module MD_CorePotent
  !***  this module is to calculate the potential and force between
  !     atoms in close-distance. Especially used for MD simulation of
  !     cascade collisions.
  !
  !     The module should be used in combination with the Force_Table_module
  !     for example  the module MD_FS_Force_Table_2010.
  !
  !
  !     REFERENCE:
  !
  !     SEE ALSO: MD_FS_Force_Table.f90
  !
  !    Author:   HOU Qing, Oct, 2012

  use MD_CONSTANTS
        implicit none


  contains
  !*********************************************************************

  !*********************************************************************
  subroutine ZBL_Pot(Z1, Z2, R, POTR, FPOTR)
  !***  PURPOSE: to calculate the potential and force by ZBL model
  !     INPUT:   Z1, Z2, the atomic numbers of the two atoms
  !              R, the distance between the two atoms, in unit of centimeter
  !
  !     OUTPUT:  POTR,  the relation between potential and the distance
  !              FPOTR, the first differential of the potential
  !

  implicit none
      real(KINDDF),intent(in)::R, Z1, Z2
      real(KINDDF),intent(out)::POTR, FPOTR
      !******* parameters: *****************************************
      real(KINDDF), parameter::ZK    = 9.0d9
      real(KINDDF), parameter::ELECHARGE  = 1.60219D-19
      real(KINDDF), parameter::KE2 = ZK*ELECHARGE*ELECHARGE*1.D7*1.D2    ! 1.D7 and 1.D2 comes from unit convert: J=>Erg, m=>cm

      real(KINDDF), parameter::A1    = 3.1998D0
      real(KINDDF), parameter::C1    = 0.18175D0
      real(KINDDF), parameter::A2    = 0.94229D0
      real(KINDDF), parameter::C2    = 0.50986D0
      real(KINDDF), parameter::A3    = 0.4029D0
      real(KINDDF), parameter::C3    = 0.28022D0
      real(KINDDF), parameter::A4    = 0.20162D0
      real(KINDDF), parameter::C4    = 0.028171D0

       real(KINDDF), parameter::AC1  = A1*C1
       real(KINDDF), parameter::AC2  = A2*C2
       real(KINDDF), parameter::AC3  = A3*C3
       real(KINDDF), parameter::AC4  = A4*C4

      !**************************************************************

      real(KINDDF)::ZZE2, as, X,faix,dfaix

           AS = 0.88534D0*CP_BOHR/(Z1**0.23D0+Z2**0.23D0)
           X= R/AS

           ZZE2=Z1*Z2*KE2

           faix=C1*DEXP(-A1*X)+C2*DEXP(-A2*X)+C3*DEXP(-A3*X)+C4*DEXP(-A4*X)
           faix=ZZE2*faix/R
           POTR=C_HALF*faix

           dfaix=AC1*DEXP(-A1*X)+AC2*DEXP(-A2*X)+AC3*DEXP(-A3*X)+AC4*DEXP(-A4*X)
           dfaix=ZZE2*dfaix/AS
           !modified accordding to NOTE-2014-12-04 in forcetable
           !FPOTR= (faix + dfaix)/(R*R)
           FPOTR= (faix + dfaix)/R

      return
  end subroutine ZBL_Pot
  !*********************************************************************


  end module MD_CorePotent
