      SUBROUTINE TIMER(HOUR,MINUTS,SECONDS)
C**** PURPOSE: To get the time between this calling of the subroutine
C              and the previous calling.
C     OUTPUT:  HOUR, MINUTS, SECODS,TOTAL
C     AUTHER:  Hou Qing
C     DATE:    Aug.16, 1995
      INTEGER HOUR, MINUTS, SECONDS, TOTAL, BEGING
      DATA    BEGING/0/
      SAVE    BEGING

        CALL TIME(TOTAL)
        IF(BEGING.EQ.0) BEGING=TOTAL
        HOUR=(TOTAL-BEGING)/3600
        MINUTS=(TOTAL-BEGING)/60-HOUR*60
        SECONDS=(TOTAL-BEGING)-HOUR*3600-MINUTS*60
C       BEGING=TOTAL

C       WRITE(6,*)'THE TIMER IS ',HOUR,':',MINUTS,':',SECONDS

        RETURN
      END SUBROUTINE

