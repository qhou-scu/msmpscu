module Express_Ext_Module
!--- DESCRIPTION:  This module is used to perform character calculations,
!                  an extension of  Express_Module and Express_LOGICAL_Module.
!                     
!                  The expression can be a set of statments with conditional
!                  statement.
!      
!                  The expression can be read from a file.  
!                  
!
!    AUTHOR: Qing HOU, 1998
!    LAST MODIFICATION: Jan., 2012       

use DATUTIL_MODULE 

integer(4), parameter, private::DELIMITER = 1
integer(4), parameter, private::VARIABLE  = 2
integer(4), parameter, private::NUMBER    = 3
integer(4), parameter, private::KEYWORD   = 4
character*32, private::token = ""
integer(4), private::ERR

integer(4), private::Lines = 0, curL = 0
integer(4), parameter, private::MxLine = 200
character*128,dimension(MxLine),private::Sentence=""
integer*4, dimension(MxLine),private::Lstatu=0
type(strDATA), pointer, private::Vars
type(xDATA), pointer, private::VVars

!**
character*32,  dimension(:),   allocatable::Varibles
double precision, dimension(:, :), allocatable::VValues
double precision, dimension(:), allocatable::Res
logical, dimension(:), allocatable::LRes
!real*8, dimension(:),allocatable, private::T
integer(4), private::dim = 1

public::Get_SubProc_form_File
public::Begin_Calculation
public::Putout_results_to_File
public::Analyse_Sentence
public::reset_Vars
public::get_Vars
private::get_Vars_Tab
private::isDelim
private::isAlpha
private::isDigit
private::isKeyWord
private::it_is_ifState
private::it_is_Equation
private::Begin_Calculate

contains
!**********************************************************************
subroutine End_Express_Ext_Module()
!--- PURPOSE: to release the momery allocated
! 
use DATUTIL_MODULE 
IMPLICIT NONE
integer::ERR

ERR = delDataChain(Vvars) 
ERR = delStrChain(Vars)
Lines = 0
curL = 0
Sentence=""
Lstatu=0

return
end subroutine End_Express_Ext_Module
!**********************************************************************


!**********************************************************************
subroutine Get_SubProc_form_File(Fname) 
!--- PURPOSE: to read sentence from a file
!    INPUT,   Fname,  the file name providing the statements
!
use STRING_OP_module
IMPLICIT NONE
   character*(*),intent(in)::Fname

!--- Local vairables
     integer::hFile
     character*128::Str

     call AvailableIOUnit(hFile) 

     open(hFile, File = Fname, status='old')
    
     DO WHILE(.true.)
        call GETINPUTSTRLINE(hFile,STR, *100)
        call delSpaces(STR)
        call UPCASE(STR)
        call Add_Sentence(STR,*100)  
     END DO
100  close(hFile)

   call Analyse_Sentence()
return
end subroutine Get_SubProc_form_File 
!**********************************************************************

!**********************************************************************
subroutine Add_Sentence(str,*) 
!--- PURPOSE: to add a statement
!    INPUT,   str,  the statement to add
!
IMPLICIT NONE
    character*(*)::str

    if(Lines.GE.MxLine) then
       print *, "Error in Express_Ext_Module: the number of statement excced the permitted number"
       return 1
    end if
    Lines = Lines+1
    Sentence(Lines)= STR
    
    return
end subroutine Add_Sentence
!**********************************************************************

!**********************************************************************
subroutine Analyse_Sentence()
!--- PURPOSE: to analyse the grammer of the sentenses
!    
use STRING_OP_module
IMPLICIT NONE
integer::I,J, ifC
character*128::Str
integer, dimension(:),allocatable::II
integer*2::word0(2)
integer*4::code
equivalence(code,word0)
integer*2, dimension(:),allocatable::word

   ifC = 0
   allocate(II(lines),word(2*lines))
   II = 0
   word = 0
   Lstatu = 0
   DO I = 1, lines
      if(Sentence(I).eq.'ELSE') then
         WORD(2*II(ifC)-1) = I
         cycle
      end if

      if(Sentence(I).eq.'ENDIF') then
         WORD(2*II(ifC)) = I
	     if(WORD(2*II(ifC)-1) .eq. 0) WORD(2*II(ifC)-1) = I
	     !Lstatu(II(ifC)) = code
	      ifC = ifC - 1
	     !code = 0
         cycle
      end if
   

      if(it_is_ifState(Sentence(I)) ) then
         J = len_trim(Sentence(I)) -5
         Str=Sentence(I)(4:J)
	     ifC = ifC + 1
	     II(ifC) = I
      else
         Str = Sentence(I)
      end if	     

      if(get_Vars_Tab(Str) .ne. 0) then
	     write(*,*) 'Syntax Error at:'
	     write(*,*) '      ', Str
       end if
   END DO 

   if(ifC .ne. 0) then
      stop 'Error: the IF/ENDIF mismatched.'
   endif 

   DO I=1,lines
      if(word(2*I) .ne. 0) then
         word0(1) = word(2*I-1) 
         word0(2) = word(2*I) 
	     Lstatu(I)=code
       end if
  END DO

  deallocate(II,word) 

return
end subroutine Analyse_Sentence
!**********************************************************************

!**********************************************************************
subroutine Begin_Calculation() 
use DATUTIL_MODULE 
IMPLICIT NONE
integer::I,ERR

!*** to beging the calculation
  allocate(Varibles(1:numb_of_Strs(Vars)) )
  allocate(VValues(numb_of_Strs(Vars),dim))
  allocate(Res(dim))
  allocate(LRes(dim))
  DO I = 1, numb_of_Strs(Vars)
     Varibles(I) = getStr(I, Vars)
     ERR = get_Vars(Varibles(I), Res)
	 VValues(I,1:dim) = Res(1:dim)
  END DO
  call Begin_Calculate(0,lines)

  DO I=1, numb_of_Strs(Vars)
     Res(1:dim) = VValues(I,1:dim) 
     ERR = resetData(I, VVars, Res)
  END DO 

  deallocate(Varibles, VValues,Res,LRes)
  return
end subroutine Begin_Calculation
!**********************************************************************

!**********************************************************************
recursive subroutine Begin_Calculate(I0,L0) 
use Express_Module
use Express_LOGICAL_Module,only:GetLogExpress
use STRING_OP_module
use DATUTIL_MODULE 
IMPLICIT NONE
integer::I,J, I0, L0, limit,ERR
character*128::Str
integer*2::word(2)
integer*4::code = 0
equivalence(code,word)

I = I0
DO WHILE(.true.)
   I = I+1
   if(I.gt.L0) then
      return
   endif
   curL = I
   if(Lstatu(I) .ne. 0) then
      code = Lstatu(I)
	  J = len_trim(Sentence(I))
	  Str=Sentence(I)(4:J-5)
      ERR=GetLogExpress(STR, LRes, Varibles, VValues)
	  if(LRes(1)) then
         code = Lstatu(I)
		 Limit = word(1)
         call Begin_Calculate(I,Limit) 
		 I = curL
		 Limit = word(2)
		 do while(.true.)
		    I=I+1
		    if(I.gt.Limit) exit
			curL = I
         end do
	  else
         code = Lstatu(I)
		 Limit = word(1)
		 do while(.true.)
		    I=I+1
		    if(I.eq.Limit) exit
			curL = I
         end do	     
		 Limit = word(2)
         call Begin_Calculate(I,Limit) 
		 I = curL
	  end if   
    else
	  call Sentence_Calculate(Sentence(I))
   end if

END DO

return
end subroutine Begin_Calculate
!**********************************************************************

!**********************************************************************
recursive subroutine Sentence_Calculate(Str) 
use Express_Module
use Express_LOGICAL_Module
use STRING_OP_module
use DATUTIL_MODULE 
IMPLICIT NONE
integer::J, N,ERR
character*(*)::Str
character*32::Tstr
character*128::Sexp

	  J = Index(Str,"=")
      Tstr = Str(1:J-1)
	  Sexp = Str(J+1:len(str))
	  ERR = GetExpress(Sexp, Res, Varibles, VValues)
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
    return
end subroutine Sentence_Calculate
!**********************************************************************

!**********************************************************************
subroutine Putout_results_to_File(Fname) 
use STRING_OP_module
IMPLICIT NONE
integer::I
character*(*)::Fname
character*128::Str
character*32::Var
double precision, dimension(:), allocatable::Value

allocate(Value(1:dim))
open(1, File = Fname, status='unknown')
DO I =1, numb_of_Strs(Vars)
    Str = getStr(I, Vars)
	Var =""
	call catenate(Var,Str)
    if( get_Vars(Var, Value) .eq. 0) then
        write(1,100) Var,Value
    end if
100 format((A,1x,E14.7))
END DO
close(1)
deallocate(Value)
return
end subroutine Putout_results_to_file
!**********************************************************************

!**********************************************************************
integer function reset_Vars(Var, NewValue) result(ERR)
IMPLICIT NONE
integer::I
character*(*)::Var
real*8,dimension(:)::NewValue

    ERR = 0
	if(.not.if_exist_Str(Vars,Var) ) then
       ERR = 1
       return
    else
	   I =where_the_Str(Vars,Var) 	
	   ERR = resetData(I, VVars, NewValue) 
    end if
	return
end function reset_Vars
!**********************************************************************

!**********************************************************************
integer function get_Vars(Var, Value) result(ERR)

IMPLICIT NONE
integer::I
character*(*)::Var
real*8,dimension(:)::Value

    ERR = 0
	if(.not.if_exist_Str(Vars,Var) ) then
       ERR = 1
       return
    else
	   I =where_the_Str(Vars,Var) 	
	   ERR = getData(I, VVars, Value) 
    end if
	return
end function get_Vars
!**********************************************************************

!**********************************************************************
integer function get_Vars_Tab(Exps) result(ERR)

IMPLICIT NONE
integer::I, I1
character*(*)::Exps
real*8, dimension(:),allocatable::T
   
   ERR = 0
   I = 0
   DO WHILE(.TRUE.)
      I = I+1 
      if(I .gt. len_trim(Exps)) then
	     return
	  end if

      if(isalpha(Exps(I:I)) ) then
           I1 = 0
		   token = ""
           do while(.not.isdelim(Exps,I+I1) )
		            if(I1.lt.len(token)) then
					   token(I1+1:I1+1) = Exps(I+I1:I+I1) 
                    else
					   ERR = 1
					   return
                    end if
					I1 = I1 + 1					    
           end do
		   I=I+I1-1
		   if (isKeyWord(token) ) then
           else
               if(.not.if_exist_Str(Vars,token) ) then	
                  call addSTR(VarS, token)
				  allocate(T(1:dim))
				  T= 0.D0
				  call addDATA(VVars, T)
				  deallocate(T)
               end if
           end if
       else 
	       if(isdigit(Exps(I:I)) ) then
	             I1 = 0 	          
				 token = ""
                 do while(.not.isdelim(Exps, I+I1) )
		            if(I1.lt.len(token)) then
					   token(I1+1:I1+1) = Exps(I+I1:I+I1) 
                    else
					   ERR = 1
					   return
                    end if	
					I1 = I1+1
                  end do
				  I=I+I1-1
            end if
      end if           
   end do   
return
end function get_Vars_tab
!**********************************************************************

!**********************************************************************
logical function isDelim(C, I)
implicit none
character*(*)::C
integer::I, J
   
     J = len_trim(C)  
     if(I .gt. J) then
	    isDelim = .true.
		return
     end if 		 
	 if(C(I:I) .eq. '=' .or. &
	    C(I:I) .eq. '+' .or. &
	    C(I:I) .eq. '-' .or. &
	    C(I:I) .eq. '*' .or. &
	    C(I:I) .eq. '^' .or. &
	    C(I:I) .eq. '/' .or. &
	    C(I:I) .eq. '%' .or. &
	    C(I:I) .eq. '(' .or. &
	    C(I:I) .eq. ')' .or. &
	    C(I:I) .eq. '>' .or. &
	    C(I:I) .eq. '<' .or. &
		C(I:I) .eq. '&' .or. &
		C(I:I) .eq. '|'          ) then

		isDelim = .true.
		if(C(I:I).eq.'+' .or. C(I:I).eq.'-') then
		   if(C(I-1:I-1).eq.'E') then
		      isDelim = .false.
           end if
         end if		               

     else
		isDelim = .false.
     end if
	 return
end function isDelim
!**********************************************************************

!**********************************************************************
logical function isAlpha(C)
implicit none
character*1::C
     
	 if(ichar(C(1:1)) .ge. ichar('A') .and. ichar(C(1:1)) .le. ichar('Z') ) then
		isAlpha = .true.
     else
		isAlpha = .false.
     end if
	 return
end function isAlpha
!**********************************************************************

!**********************************************************************
logical function isDigit(C)
implicit none
character*1::C
     
	 if(ichar(C(1:1)) .ge. ichar('0') .and. ichar(C(1:1)) .le. ichar('9') ) then
		isDigit = .true.
     else
		isDigit= .false.
     end if
	 return
end function isDigit
!**********************************************************************

!**********************************************************************
logical function isKeyWord(Str)
implicit none
character*(*)::Str
     
	 if(Str(1:len_trim(str)) .eq. 'SIN'      .or.  &
	    Str(1:len_trim(str)) .eq. 'COS'      .or.  &
	    Str(1:len_trim(str)) .eq. 'ASIN'     .or.  &
	    Str(1:len_trim(str)) .eq. 'ACOS'     .or.  &
	    Str(1:len_trim(str)) .eq. 'TAN'      .or.  &
	    Str(1:len_trim(str)) .eq. 'ATAN'     .or.  &
	    Str(1:len_trim(str)) .eq. 'EXP'      .or.  &
	    Str(1:len_trim(str)) .eq. 'LOG'      .or.  &
	    Str(1:len_trim(str)) .eq. 'LOG10'        ) then
	    
		isKeyWord = .true.
     else
		isKeyWord = .false.
     end if
	 return
end function isKeyWord
!*********************************************************************

!*********************************************************************
logical function it_is_ifState(Str) result(ERR)
implicit none
character*(*)::Str
character*4::S
integer::I,J
   
     J =0
	 DO I=1,len_trim(STR)
	    if(Str(I:I) .eq. '(') J=J+1
        if(STR(I:I) .eq. ')' ) J=J-1
        if(J .lt. 0) then
		   Stop 'Error: mismatch of bracket in IF statement' 		        
        end if
     END DO

     if(J .gt. 0) then
		Stop 'Error: mismatch of bracket in IF statement' 		        
     end if

	 J = scan(STR,'IF(') 
	 if(J .eq. 1) then
	    ERR = .true.
        J = len_trim(STR)
	    S=str(J-3:J)
	    if(S.ne.'THEN') then
		   Stop 'Error: Sytax Error in IF statement' 		        
        end if
     else
	    ERR = .false.
     end if
	 return
end function it_is_ifState
!*********************************************************************

!*********************************************************************
logical function it_is_EQUATION(Str) result(ERR)
implicit none
character*(*)::Str
integer::I,J
   
     J =0
	 DO I=1,len_trim(STR)
	    if(Str(I:I) .eq. '(') J=J+1
        if(STR(I:I) .eq. ')' ) J=J-1
        if(J .lt. 0) then
		   write(*,*) 'Error: mismatch of bracket' 		        
        end if
     END DO
     if(J .gt. 0) then
		write(*,*) 'Error: mismatch of bracket' 		        
     end if
	 if(scan(STR,'=') .gt. 0) then
	    ERR = .true.
     else
	    ERR = .false.
     end if
	 return
end function it_is_EQUATION
!**********************************************************************


end module Express_Ext_Module
