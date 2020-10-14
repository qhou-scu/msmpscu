      INTEGER FUNCTION NSTOKEN(STR,CHA)
C**** PURPOSE: to iobtain the index of the token in a str
C     INPUT:   STR a string
C              CHA the token
C     OUTPUT:  the index where CHA is
C     AUTHER:  Hou Qing
C     DATE:    Aug.16, 1995

      CHARACTER *(*) STR
      CHARACTER *1   CHA
      INTEGER I

      DO I=1,LEN(STR)
         IF(STR(I:I).EQ.CHA) goto 1
      ENDDO
 1    NSTOKEN=I

      RETURN
      END FUNCTION
