
C*DECK DTOS
      SUBROUTINE DTOS(XD, XS,N)
C***BEGIN PROLOGUE  DTOS
C***PURPOSE:  To convert double precision array to single precision array
C***LIBRARY   OTHERS
C***INPUT:    XD:       the double precision array
C***OUPUT:    XS:       the single precision array
C***AUTHOR    Hou Qing, Insti. of Nucle. Sci. and Tech., Sichuan University
C
       DOUBLE PRECISION XD(*)
       REAL             XS(*)
       INTEGER          I,N

       DO I=1,N
          XS(I)=XD(I)
       ENDDO
       RETURN
       END SUBROUTINE





























