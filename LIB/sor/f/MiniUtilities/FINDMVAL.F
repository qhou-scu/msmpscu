      SUBROUTINE FindMinMax (dataset, numdat, minval, maxval)
      IMPLICIT NONE
      REAL dataset(*), minval, maxval
      INTEGER numdat, i

      minval = dataset(1)
      maxval = dataset(1)
      DO I = 1, numdat
        IF (dataset(I) .LT. minval)  minval = dataset(i)
        IF (dataset(I) .GT. maxval)  maxval = dataset(i)
      END DO
      END SUBROUTINE
