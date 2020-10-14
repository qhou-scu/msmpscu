recursive subroutine Analyse_Sentence(I0)
use STRING_OP_module
IMPLICIT NONE
integer::I0,I,IB,IN
integer*4::IFCODE
integer*2::WORD(2)
equivalence(IFCODE,WORD)

I = I0
DO while(.true.)
   I = I+1
   if(I.gt.Lines) then
      if(ifCount .gt. 0) write(*,*) 'Error: IF statement without ENDIF' 
      if(ifCount .lt. 0) write(*,*) 'Error: ENDIF statement without IF' 
      exit
   endif
   curL = I
   if(Sentence(I).eq.'ELSE') then
      return
   end if
   if(Sentence(I).eq.'ENDIF') then
      ifCount = ifCount-1
      return
   end if
   
   if(it_is_ifState(Sentence(I)) ) then
	  IB = ifCount
      ifCount = ifCount + 1
	  IN = I
      call Analyse_Sentence(IN)
	  IN = curL
	  if(ifCount .gt. IB) then
	     WORD(1)=curL
         call Analyse_Sentence(IN)
         if(ifCount .ne. IB) then
		    write(*,*) 'Error: IF without ENDIF'
         else
		    WORD(2) =curL
         end if
      end if
	  Lstatu(I) = IFCODE
	  I = curL
   else
      if(scan(Sentence(I),'=') .eq. 0) then
 	     write(*,*) 'Syntax Error at:'
		 write(*,*) '      ', Sentence(I)
      end if
      if(get_Vars_Tab(Sentence(I)) .ne. 0) then
	     write(*,*) 'Syntax Error at:'
		 write(*,*) '      ', Sentence(I)
      end if
   end if 
  
end do
return
end subroutine Analyse_Sentence

subroutine Begin_Calculation() 
use Express_Module
use STRING_OP_module
use DATUTIL_MODULE 
IMPLICIT NONE
integer::I, J, N,ERR
character*80::Str
character*32:: Tstr
character*32,  dimension(:),   allocatable::Varibles
double precision, dimension(:, :), allocatable::VValues
double precision, dimension(:), allocatable::Res


!*** to beging the calculation
allocate(Varibles(1:numb_of_Strs(Vars)) )
allocate(VValues(numb_of_Strs(Vars),dim))
allocate(Res(dim))
DO I = 1, numb_of_Strs(Vars)
   Str = getStr(I, Vars)
   Varibles(I) = ""
   call catenate(Varibles(I),Str)
END DO

DO I=1,lines
    Str = Sentence(I)
	 J = Index(Str,"=")
     Tstr = Str(1:J-1)
	 Str = Str(J+1:len(str))
	 ERR = GetExpress(Str, Res, Varibles, VValues)
	 N = 0
   	 Do J=1, numb_of_Strs(Vars)
	      if(Tstr .eq. Varibles(J)) then
		     N = J
		     exit
          end if 
     end do
	 DO J=1, dim
	    VValues(N,J)=Res(J) 
     END DO
END DO

DO I=1, numb_of_Strs(Vars)
   Res(1:dim) = VValues(I,1:dim) 
   ERR = resetData(I, VVars, Res)
END DO 

deallocate(Varibles,VValues, Res)

return
end subroutine Begin_Calculation
