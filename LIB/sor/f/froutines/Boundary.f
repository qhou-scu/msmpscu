      subroutine IfInBoundary(eq, x,y,z, ERR)
      !ms$attributes dllexport :: IFINBOUNDARY
      
	  use Express_Module
      use Express_LOGICAL_Module
      use STRING_OP_module
      implicit none

	 integer(1) eq(*)
	 real*4::x, y, z
	 !Local variables
	  character*80::str 
	  character*32,dimension(3)::V=(/'X','Y','Z'/)
      real*8,dimension(3,1)::VV
	  integer I,ERR
	  logical res(1)

      VV(1,1) = x 
      VV(2,1) = y 
      VV(3,1) = z
	
      str = ""
      DO I=1,80
	   if(eq(I) .eq. 0) exit
	   str(I:I) = char(eq(I))
      END DO

	ERR = GetLogExpress(str,res,V,VV)
	if(res(1)) then
	   ERR = 1
	else  
	   ERR = 0
      end if
         
	end

