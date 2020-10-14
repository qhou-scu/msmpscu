C      CHARACTER *80 STR
C      STR='ABDSsd L A sadsa SA 12 Slkmasa 131'
C      print * str
C      call lowcase(str)
C      print * str
C      stop
C      end

      SUBROUTINE LOWCASE(STR)
C**** PURPOSE: to convert a string into lowcase
C     INPUT:   STR a string with upper or low case
C     OUTPUT:  STR the string with low case
C     AUTHER:  Hou Qing
C     DATE:    DEC.27, 1995

      CHARACTER *(*) STR
      INTEGER I,IA

      DO I=1,LEN(STR)
         IA=iachar(STR(I:I))
         IF(IA.GE.iachar('A').AND.IA.LE.iachar('Z')) THEN
            IA=IA+iachar('a')-iachar('A')
            STR(I:I)=achar(IA)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE
