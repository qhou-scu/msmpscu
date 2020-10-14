C*DECK STOD
      SUBROUTINE STOD(XS, XD,N)
C***BEGIN PROLOGUE  STOD
C***PURPOSE:  To convert double precision array to single precision array
C***LIBRARY   OTHERS
C***INPUT:    XS:       the SINGLE precision array
C***OUPUT:    XD:       the DOUBLE precision array
C***AUTHOR    Hou Qing, Insti. of Nucle. Sci. and Tech., Sichuan University
C
       DOUBLE PRECISION XD(*)
       REAL             XS(*)
       INTEGER          I,N

       DO I=1,N
          XD(I)=XS(I)
      ENDDO
      RETURN
      END SUBROUTINE
