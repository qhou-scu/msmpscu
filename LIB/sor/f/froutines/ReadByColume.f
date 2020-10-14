      subroutine NumberofColume(fname, nCol)
      !ms$attributes dllexport :: NUMBEROFCOLUME
!***  PURPOSE:   to check how many colume the table in text file has
!
!     INPUT:   fname,
!     OUTPUT:  nCol, number of colume

	implicit none
	integer(1)::fname(*)
	integer::nCol

	!local variables
	  character*256::temp
	  character*(32*100)::str
	  character*32::STRNUMB(400)

	  integer I,N

	  I = 0
	  temp = ""
	  DO WHILE(.TRUE.)
	   I = I+1
	   if(fname(I) .EQ. 0) exit
	   if(I .GT. 256) exit
	   temp(I:I) = char(fname(I))
      END DO 

       open(unit = 1, file = temp, status='old') 
	   DO WHILE(1)
	      READ(1,1, err=100) str
	      if(len_trim(str) .le. 0) cycle
 	      str = adjustl(str)
	      if(str(1:1) .NE. '!') exit
	   END DO
1	   FORMAT(A3200)
100    call EXTRACT_NUMB(str, 100, n, STRNUMB)
	  close(1)
	nCol = n
	return
	end
      
      subroutine NumberofLine(fname, nLine)
      !ms$attributes dllexport :: NUMBEROFLINE
!***  PURPOSE:   to check how many lines the table in text file has
!
!     INPUT:   fname,
!     OUTPUT:  nLine, number of line

	implicit none
	integer(1)::fname(*)
	integer::nLine

	!local variables
	character*256::temp
	character*(32*100)::str

	integer I

	I = 0
	temp = ""
	DO WHILE(.TRUE.)
	   I = I+1
	   if(fname(I) .EQ. 0) exit
	   if(I .GT. 256) exit
	   temp(I:I) = char(fname(I))
      END DO 

      open(unit = 1, file = temp, status='old') 
	!Read the head
	DO WHILE(1)
	    READ(1,1,err=100) str
	    if(len_trim(str) .le. 0) cycle
 	    str = adjustl(str)
	    if(str(1:1) .NE. '!') exit
	END DO

100	nLine = 1
	DO WHILE(1)
	   READ(1,1,err=200) str
	   if(len_trim(str) .gt. 0) then
	      nLine = nLine+1
         else
	      cycle
         end if
1	   FORMAT(A3200)
	END DO
200	close(1)

	return
	end

      subroutine ReadByColume(fname, nLine, nStep, nCol, nVar,flag,var)
      !ms$attributes dllexport :: READBYCOLUME
!***  PURPOSE:   to read from file by colum set by flag
!
!     INPUT:   fname,
!              nCol, number of colume
!              nLine, number of the lines in the file
!              nStep, for how much lines to read once
!              flag, the identity of colume to be read
!     OUTPUT     
!              var, buffer
	implicit none
	integer::nLine, nStep,nCol, nVar
	integer(1)::fname(*)
	integer::flag(*)
	real*4::var(*)

	!local variables
	character*256::temp
	character*(32*100)::str
	real*4, dimension(:), pointer::tvar
	integer I,J,K

	I = 0
	temp = ""
	DO WHILE(.TRUE.)
	   I = I+1
	   if(fname(I) .EQ. 0) exit
	   if(I .GT. 256) exit
	   temp(I:I) = char(fname(I))
      END DO 

      open(unit = 1, file = temp, status='old') 

	!Read the head
	 allocate(tvar(nCol))
	 DO WHILE(1)
	    READ(1,1,err=100) str
	    if(len_trim(str) .le. 0) cycle
 	    str = adjustl(str)
	    if(str(1:1) .NE. '!') exit
	  END DO
1	  FORMAT(A3200)

100   BACKSPACE (1)

	    K = 0
	    DO I=1, nLine
	       READ(1,*) tvar(1:nCol)
	       DO J=1, nVar
	            K= K+1
	            var(K) = tvar(flag(J))
             ENDDO
	       DO J=1, nStep-1
	          READ(1,*) tvar(1:nCol)
             END DO 
	    END DO
	  deallocate(tvar)

	close(1)
	return
	end
