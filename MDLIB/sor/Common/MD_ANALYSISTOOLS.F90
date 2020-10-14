  module MD_ANALYSISTOOLS
  use MD_CONSTANTS
  implicit none

  contains

  SUBROUTINE PairDistribution(N,ATOM1,ATOM2,BOX,ITYP,XP,DELR,MAXBIN,G)
  !***  PURPOSE: to calculate the pair distance distribution between ATOM1 and ATOM2
  !    INPUT:    N, the number of particle(size of xp)
  !               ATOM1,    the type of atom1
  !               ATOM2,    the type of atom2
  !               BOX,      the box size, in CGS
  !               ITYP,     the type index for each atom
  !               XP,       the coordinate of particles in CGS
  !               DELR,     the thick of shell. (for r**2-r**2+width as a shell)
  !               MAXBIN,   the maxma number of bins
  !   OUTPUT:     G,          the number of particle in 4PI*R**2*dR
  IMPLICIT NONE
  integer::N, ATOM1, ATOM2, MAXBIN
  integer, dimension(:)::ITYP
  real(KINDDF)::DELR
  real(KINDDF), dimension(:,:)::XP
  real(KINDDF), dimension(:)::BOX
  real(KINDDF), dimension(:)::G
  !---local variables
  integer I,J,IP, K
  real(KINDDF)  CX, CY, CZ, RIJ

      G = 0.D0
       DO I=1, N-1
         if(ITYP(I) .EQ. ATOM1) then
            if(ATOM1 .EQ. ATOM2) then !The same kind of atoms
               IP=I+1
            else                      !Diffferent kind of atoms
               IP = 1
            end if
            DO J=IP, N
               if(ITYP(J) .EQ. ATOM2) then
                    CX=XP(I,1)-XP(J,1)
                    CY=XP(I,2)-XP(J,2)
                    CZ=XP(I,3)-XP(J,3)
                    IF(DABS(CX).GT.0.5D0*BOX(1) ) CX=CX-DSIGN(BOX(1),CX)
                    IF(DABS(CY).GT.0.5D0*BOX(2) ) CY=CY-DSIGN(BOX(2),CY)
                    IF(DABS(CZ).GT.0.5D0*BOX(3) ) CZ=CZ-DSIGN(BOX(3),CZ)
                    RIJ=CX*CX + CY*CY + CZ*CZ
                    RIJ=DSQRT(RIJ)
                    K =INT(RIJ/DELR)+1
                    IF(K .LE.MAXBIN) THEN
                       G(K) = G(K)+1.D0
                    END IF
                 end if
             END DO
             end if
        END DO

        RETURN
  END SUBROUTINE PairDistribution

  SUBROUTINE Distribution_Around_aPoint(X0, Y0, Z0, N,ATOM1,BOX,ITYP,XP,DELR,MAXBIN,G,GI)
  !***  PURPOSE:  to calculate the particle distribution around a fixed point
  !    INPUT:      X0, Y0, the point around which we will find the distribution
  !                N, the number of particle(size of xp)
  !               ATOM1,    the type of atom1 will be concerned
  !               BOX,      the box size, in CGS
  !               ITYP,     the type index for each atom
  !               XP,       the coordinate of particles, in CGS
  !               DELR,     the thick of shell. (for r**2-r**2+width as a shell)
  !               MAXBIN,   the maxma number of bins
  !   OUTPUT:     G,          the number of particle in 4PI*R**2*dR
  !               GI,          the number of particle with sphere of radia R
  IMPLICIT NONE
  real(KINDDF)::X0, Y0, Z0
  integer::N, ATOM1, MAXBIN
  integer, dimension(:)::ITYP
  real(KINDDF)::DELR
  real(KINDDF), dimension(:,:)::XP
  real(KINDDF), dimension(:)::BOX
  real(KINDDF), dimension(:)::G
  real(KINDDF), dimension(:)::GI
  !---local variables
  integer I, K
  real(KINDDF)  CX, CY, CZ, RIJ

      G = 0.D0

       DO I=1, N
         if(ITYP(I) .EQ. ATOM1) then
            CX=XP(I,1)-X0
            CY=XP(I,2)-Y0
            CZ=XP(I,3)-Z0
            IF(DABS(CX).GT.0.5D0*BOX(1) ) CX=CX-DSIGN(BOX(1),CX)
            IF(DABS(CY).GT.0.5D0*BOX(2) ) CY=CY-DSIGN(BOX(2),CY)
            IF(DABS(CZ).GT.0.5D0*BOX(3) ) CZ=CZ-DSIGN(BOX(3),CZ)
            RIJ=CX*CX + CY*CY + CZ*CZ
            RIJ=DSQRT(RIJ)
            K =INT(RIJ/DELR)+1
            IF(K .LE.MAXBIN) THEN
               G(K) = G(K)+1.D0
            END IF
          end if
        END DO

           GI = 0.D0
        DO K=1, MAXBIN
           GI(K) = sum(G(1:K))
        END DO
        RETURN
  END SUBROUTINE Distribution_Around_aPoint

  end module MD_ANALYSISTOOLS
