module Express_Module
!--- DESCRIPTION:  This module is used to perform character calculations
!                  of a given expression. 
!                  
!
!    AUTHOR: Qing HOU, 1998
!    LAST MODIFICATION: Jan., 2012       


integer(4), parameter, private::DELIMITER = 1
integer(4), parameter, private::VARIABLE  = 2
integer(4), parameter, private::NUMBER    = 3
integer(4), parameter, private::KEYWORD   = 4
character*32, private::token = ""
integer(4), private::pToken = 0
integer(4), private::tok_type = 0
integer(4), private::nVars = 0
integer(4), private::ERR

character*1,  dimension(:), pointer, private::ExpS
character*32, dimension(:), pointer, private::Vars
double precision, dimension(:,:), pointer, private::VVars

public::GetExpress
private::get_token
private::isDelim
private::isAlpha
private::isKeyWord
private::Level2
private::Level3
private::Level4
private::Level5
private::Level6
private::Primitive
private::arith
private::unary
private::Which_var
private::Def_func


contains

!*****************************************************************************
integer  function  GetExpress(EQ, Result, Varibles, VValues) result(ERR)
use STRING_OP_module, only:catenate, upcase
implicit none
 character*(*)::EQ
 integer::I, COUNT
 character*32,  dimension(:),   optional::Varibles
 double precision, dimension(:, :), optional::VValues
 double precision,dimension(:)::Result

 ERR = 0 
 allocate(ExpS(len_trim(EQ)) )
 ExpS = ' '
 count = 0
 do I =1, len_trim(EQ)
    if(EQ(I:I) .ne. ' ') then
	   count = count + 1
	   ExpS(count) = EQ(I:I)
    end if 
	if(ichar(EQ(I:I)) .ge. ichar('a')  .and. ichar(EQ(I:I)) .le. ichar('z')) then
	   ExpS(count) = char(ichar(EQ(I:I))-ichar('a')+ichar('A'))

    end if
 end do
!*** to get the number of variable
 nVars = 0 
 if(present(Varibles)) then 
    nVars = size(Varibles) 
    count = size(result)
    allocate(Vars(nVars),VVars(nVars,count) )
    do I=1,nVars
       VarS(I)=Varibles(I)(1:len_trim(Varibles(I)))
       call upcase(VarS(I)) 
	   VVars(I,1:count) = VValues(I,1:count)
    end do 
 end if

!*** to beging the calculation
 token = ""
 pToken  = 1
 ERR = get_token()
 ERR = Level2(result)

 deallocate(ExpS)
 if(associated(Vars)) deallocate(Vars)
 if(associated(VVars)) deallocate(VVars)

 return
end function GetExpress
!*****************************************************************************

!*****************************************************************************
integer function get_token() result(ERR)
IMPLICIT NONE
integer::I, I1
   
   ERR = 0     
   I = pToken 
   if(I .gt. size(Exps)) then
      tok_type = DELIMITER
	  ERR = -1
	  return
   end if

   if(isDelim(Exps,I) ) then
      tok_type = DELIMITER
      token    = Exps(I)
	  pToken   = I+1
   else if(isalpha(Exps(I)) ) then
           I1 = 0
		   token = ""
           !--- to extract the token
           do while(.not.isdelim(Exps,I+I1) )
		            if(I1.lt.len(token)) then
					   token(I1+1:I1+1) = Exps(I+I1) 
                    else
					   ERR = 1
					   return
                    end if
					I1 = I1 + 1					    
           end do
           !--- check the type of  the token
		   if (isKeyWord(token) ) then
		       tok_type = KEYWORD
           else
		       tok_type = VARIABLE
           end if
           pToken = I+I1

        else if(isdigit(Exps(I)) ) then
	             I1 = 0 	          
				 token = ""
                 do while(.not.isdelim(Exps, I+I1) )
		            if(I1.lt.len(token)) then
					   token(I1+1:I1+1) = Exps(I+I1) 
                    else
					   ERR = 1
					   return
                    end if	
					I1 = I1+1				    
                  end do
		          tok_type = NUMBER
                  pToken    = I+I1
   end if           
      

return
end function get_token
!*****************************************************************************

!*****************************************************************************
logical function isDelim(C, I)
implicit none
character*1, dimension(:)::C
integer::I, J
   
     J = size(C)  
     if(I .gt. J) then
	    isDelim = .true.
		return
     end if 		 
	 if(C(I) .eq. '+' .or. &
	    C(I) .eq. '-' .or. &
	    C(I) .eq. '*' .or. &
	    C(I) .eq. '^' .or. &
	    C(I) .eq. '/' .or. &
	    C(I) .eq. '%' .or. &
	    C(I) .eq. '(' .or. &
	    C(I) .eq. ')'         ) then

		isDelim = .true.
		if(C(I).eq.'+' .or. C(I).eq.'-') then
		   if(C(I-1).eq.'E') then
		      isDelim = .false.
           end if
         end if		               

     else
		isDelim = .false.
     end if
	 return
end function isDelim
!*****************************************************************************

!*****************************************************************************
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
!*****************************************************************************

!*****************************************************************************
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
!*****************************************************************************


!*****************************************************************************
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
!*****************************************************************************

!*****************************************************************************
recursive integer function Level2(result) result(ERR)
implicit none
  double precision, dimension(:)::result
  double precision, dimension(:), allocatable::hold
  character*1::OP

    ERR = Level3(result)
	allocate(hold(size(result)) )
    do while(.true.)
	   OP = token(1:1)
	   if( OP.eq.'+' .or. OP.eq.'-') then
	       ERR = get_token()
	       ERR = Level3(hold)
	       call arith(OP, result, hold)
       else if( OP.eq.')') then
	        ERR = get_token()
            exit
       else
	       exit
       end if
    end do
	deallocate(hold)
	return

end function Level2
!*****************************************************************************

!*****************************************************************************
recursive integer function Level3(result) result(ERR)
implicit none
  double precision, dimension(:)::result
  double precision, dimension(:), allocatable::hold
  character*1::OP

    ERR = Level4(result)
	allocate(hold(size(result)) )
	do while (.true.)
	   OP = token(1:1)
	   if(OP.eq.'*' .or. OP.eq.'/' .or. OP.eq.'%') then
	      ERR = get_token()
	      ERR = Level4(hold)
	      call arith(OP, result, hold)
        else
		  exit
       end if
    end do
	deallocate(hold)
	return
end function Level3
!*****************************************************************************

!*****************************************************************************
recursive integer function Level4(result) result(ERR)
implicit none
  double precision, dimension(:)::result
  double precision, dimension(:), allocatable::hold
  character*1::OP

    ERR = Level5(result)
	allocate(hold(size(result)) )
	OP = token(1:1)
	if(OP .eq. '^') then
	   ERR = get_token()
	   ERR = Level4(hold)
	   call arith(OP, result, hold)
    end if 
	deallocate(hold)
	return
end function Level4
!*****************************************************************************

!*****************************************************************************
recursive integer function Level5(result) result(ERR)
implicit none
  double precision, dimension(:)::result
  character*1::OP

    OP = ''
	if(tok_type .eq. DELIMITER .and. (token(1:1).eq.'+' .or. token(1:1).eq.'-') ) then
	   OP = token(1:1)
	   ERR = get_token()
    end if
	ERR = level6(result)
	if(len_trim(OP) .gt. 0) then
	   call unary(OP, result)
    end if
	return
end function Level5
!*****************************************************************************

!*****************************************************************************
recursive integer function Level6(result) result(ERR)
implicit none
  double precision, dimension(:)::result

	if(tok_type .eq. DELIMITER .and. token(1:1).eq.'(' ) then
	   ERR = get_token()
	   ERR = Level2(result)
	   !if(token(1:1) .ne. ')') then
	   !   ERR = 60
       !end if
    else
	   ERR = primitive(result) 
    end if

	return
end function Level6
!*****************************************************************************

!*****************************************************************************
integer function Primitive(result) result(ERR)
implicit none
 double precision, dimension(:)::result
 double precision::R

   select case(tok_type)
          case(NUMBER) 
               read(token, *)R
               result = R
	           ERR = get_token()

          case(VARIABLE) 
               result = VVars(Which_var(token),1:size(result))
	           ERR = get_token()

          case(KEYWORD) 
		       ERR = Def_func(result)
               !function end with ), pToken has move advanced 
               !We do not call get_token here
	           !ERR = get_token()
		       
   end select	           
end function Primitive
!*****************************************************************************

!*****************************************************************************
subroutine arith(OP, result, hold)
implicit none
character*1::OP
double precision, dimension(:)::result, hold
   
   select case(OP)
          case('-')
		       result = result - hold
          case('+')
		       result = result + hold
          case('*')
		       result = result * hold
          case('/')
		       result = result / hold
          case('%')
		       result = mod(result,hold)
          case('^')
		       result = result ** hold
   end select
   return          

end subroutine arith
!*****************************************************************************

!*****************************************************************************
subroutine unary(OP, RESULT)
implicit none
character*1::OP
double precision, dimension(:)::result
   
      if(OP .eq. '-') then
	     result = - result
      end if
	  return
end subroutine unary	  		     
!*****************************************************************************

!*****************************************************************************
integer function Which_var(varName) result(IRES)
IMPLICIT NONE
character*(*)::varName
integer::I
character*32::tstr
   
   do I=1, nVars
      if(varName(1:len_trim(varName)) .eq. VarS(I)(1:len_trim(VarS(I))) ) then
	     exit
      end if
   end do
   if(I .gt. nVars) I = 0 
   IRES = I
   return
end function Which_var
!*****************************************************************************
	      
!*****************************************************************************
recursive integer function Def_func(result) result(ERR)
implicit none
character*32::OP
double precision, dimension(:)::result
double precision, dimension(:), allocatable::hold

   ERR = 0
   OP  = token
   select case(OP)
          case('SIN')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dsin(hold)
			   deallocate(hold)
			   return
          case('COS')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dcos(hold)
			   deallocate(hold)
			   return

          case('TAN')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dtan(hold)
			   deallocate(hold)
			   return

          case('ASIN')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dasin(hold)
			   deallocate(hold)
			   return
          case('ACOS')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dacos(hold)
			   deallocate(hold)
			   return

          case('ATAN')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = datan(hold)
			   deallocate(hold)
			   return
!*** Log and exp
          case('LOG')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dlog(hold)
			   deallocate(hold)
			   return

          case('LOG10')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dlog10(hold)
			   deallocate(hold)
			   return

          case('EXP')
		       allocate(hold(1:size(result)))
	           ERR = get_token()
			   ERR = Level6(hold)
		       result = dexp(hold)
			   deallocate(hold)
			   return


   end select
   return          

end function Def_func
!*****************************************************************************



end module Express_Module
