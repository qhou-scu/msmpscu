      SUBROUTINE Catenate( outstr, s)
!***  purpose: to cat tow strings
!     input:   outstr,s these strings
!     output:  outstr:  the string of outstr(in)+s
      CHARACTER * (*) outstr, s
      INTEGER LenO, LenS, i
      LenO = LEN_TRIM(outstr) + 1
      LenS = LEN_TRIM(s)

      DO i = 1, Lens
        outstr(LenO:LenO) = s(i:i)
        LenO = LenO + 1
      END DO
      END SUBROUTINE
