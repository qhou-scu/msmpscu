      subroutine GetValue(eq, value,ERR)
      !ms$attributes dllexport :: GETVALUE
	use Express_Module
      implicit none

	integer(1) eq(*)
	real*4 value
	!Local variables
	character*80::str 
	integer I,ERR
	double precision::res(1)
	
	
      str = ""
      DO I=1,80
	   if(eq(I) .eq. 0) exit
	   str(I:I) = char(eq(I))
      END DO

	ERR = GetExpress(EQ=str, Result=res)
	value = res(1)
         
	end


      subroutine Transform(eq,n,x,y,f,ERR)
      !ms$attributes dllexport :: TRANSFORM
      
	use Express_Module
      use STRING_OP_module
      implicit none

	integer(1) eq(*)
	integer::n
	real*8::x(*), y(*), f(*)
	!Local variables
	character*80::str 
	character*32,dimension(2)::V=(/'X','Y'/)
      real*8,dimension(:,:), pointer::VV
	real*8, dimension(:),pointer::res

	integer I,ERR


      str = ""
      DO I=1,80
	   if(eq(I) .eq. 0) exit
	   str(I:I) = char(eq(I))
      END DO

      
	allocate(VV(2,n), res(n))

	DO I=1,n
          VV(1,I) = x(I) 
          VV(2,I) = y(I)
      ENDDO
			 
	ERR = GetExpress(str,res,V,VV)


	DO I=1,n
          f(I) = res(I) 
      ENDDO

	deallocate(VV, res);
               
	end
