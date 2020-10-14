
      subroutine FortranRAND(var)
      !ms$attributes dllexport :: FORTRANRAND

!***  PURPOSE:   to to create random number between 0-1
	implicit none
	real*4::var
	integer::flag = 0

	 if(flag .eq. 0) then
	    call RANDOM_SEED()
	    flag = 1
        end if

	  call RANDOM_NUMBER(var);
	return
	end

