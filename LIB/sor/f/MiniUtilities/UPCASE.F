      SUBROUTINE UPCASE(STR)
!**** PURPOSE: to convert a string into upper case
!     INPUT:   STR a string with upper or low case
!     OUTPUT:  STR the string with upper case
!     AUTHER:  Hou Qing
!     DATE:    DEC.27, 1995

      CHARACTER *(*) STR
      INTEGER I,IA

        DO I=1,LEN_TRIM(STR)
           IA=iachar(STR(I:I))
           IF(IA.GE.iachar('a').AND.IA.LE.iachar('z')) THEN
              IA=IA+iachar('A')-iachar('a')
              STR(I:I)=achar(IA)
           ENDIF
        ENDDO

      RETURN
      END SUBROUTINE UPCASE
