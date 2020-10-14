      SUBROUTINE Combine( outstr, s)
!***   PORPOSE: to connect two string with  a space
      CHARACTER * (*) outstr, s
      INTEGER LenO, LenS, i
      LenO = LEN_TRIM(outstr) + 2
      LenS = LEN_TRIM(s)

      DO i = 1, Lens
        outstr(LenO:LenO) = s(i:i)
        LenO = LenO + 1
      END DO
      END SUBROUTINE
