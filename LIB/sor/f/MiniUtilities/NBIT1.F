c      DO I=0,10
c      print * 1023409,NBIT1(1023409,I)
c      ENDDO
c      STOP
c      END

      INTEGER FUNCTION NBIT1(NUMBER,N)
C**** PURPOSE: to get the number at N bit of NUMBER
C               Example: NBIT1(10283,0)=3
C                        NBIT1(10283,1)=8
C                        NBIT1(10283,4)=1
C                        NBIT1(10283,5)=0  ...
C
C     INPUT:   NUMBER a inetegr
C              N      at Nhe position
C     OUTPUT:  NBIT1
C     AUTHER:  Hou Qing
C     DATE:    Jan.10, 1995

      INTEGER N,I,NUMBER

      I=NUMBER/10**N
      NBIT1=I-(I/10)*10

      RETURN
      END FUNCTION
