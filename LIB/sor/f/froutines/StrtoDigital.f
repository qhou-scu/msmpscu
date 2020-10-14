      

      subroutine StrtoInts(str, len, n, var)
      !ms$attributes dllexport :: STRTOINTS

!***  PURPOSE:   to extract number from a string
!
!     INPUT:    str,
!               len, len of the string
!     OUTPUT     n, number of digital gotten
!                var, the numbers
	implicit none
	integer::ISTR        ! an external function converting string to integer
	integer::len, n
	integer(1)::STR(*)
	integer*4::var(*)

	!local variables
	character*256::temp
	character*32::STRNUMB(20)
	integer I,J

      !allocate( STRNUMB(128))

      J = min(len,256)
	temp = ""
	DO I=1,J
	   temp(I:I) = char(str(I))
 	   if(str(I) .EQ. 0) exit
      END DO 

       call EXTRACT_NUMB(temp, 20, n,STRNUMB)

	    DO I=1, n
             var(I) = ISTR(STRNUMB(I))
          END DO

	!deallocate(STRNUMB) 
	return
	end

      subroutine StrtoFloats(str, len, n, var)
      !ms$attributes dllexport :: STRTOFLOATS

!***  PURPOSE:   to extract number from a string
!
!     INPUT:    str,
!               len, len of the string
!     OUTPUT     n, number of digital gotten
!                var, the numbers
	implicit none
	real*8::DRSTR        ! an external function converting string to integer
	integer::len, n
	integer(1)::STR(*)
	real*4::var(*)

	!local variables
	character*256::temp
	character*32::STRNUMB(20)
	integer I,J

      !allocate( STRNUMB(128))

      J = min(len,256)
	temp = ""
	DO I=1,J
	   temp(I:I) = char(str(I))
	   if(str(I) .EQ. 0) exit
      END DO 

       call EXTRACT_NUMB(temp, 20, n,STRNUMB)

	 DO I=1, n
          var(I) = DRSTR(STRNUMB(I))
       END DO

	!deallocate(STRNUMB) 
	return
	end