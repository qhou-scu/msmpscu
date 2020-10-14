module Express_LOGICAL_Module
!--- DESCRIPTION:  This module is used to perform logical calculations
!                  of a given inequations. For exmaple, to calculate  
!                  if inequations or equation "A*B < C^2" is true or not. 
!                  For more complex, if the statement "A*B < C^2 && D^3=exp(C)" 
!                  is true or not. 
!                  
!
!    AUTHOR: Qing HOU, 1998
!    LAST MODIFICATION: Jan., 2012       
!--- 
!--- Local parameter and variables
integer(4), parameter, private::DELIMITER = 1
integer(4), parameter, private::TERM      = 2
character*256, private::token = ""
integer(4), private::pToken = 1
integer(4), private::tok_type = DELIMITER
integer(4), private::nVars = 0
integer(4),private::dim = 1

character*1,  dimension(:), pointer, private::ExpS              !--- The inequation expression
character*32, dimension(:), pointer, private::Vars              !--- The variable strings in the inequation
double precision, dimension(:,:), pointer, private::VVars       !--- The values of the variables
!logical, dimension(:,:), pointer, private::result

private::get_token
private::Get_Logical_Vaulen
private::getTermValue
private::isComp_true
private::Comp_number
private::isCompSymbol
private::isDelim
private::isLogicalNo
private::arith

public::GetLogExpress

!--------
contains

!**************************************************************************************
integer  function  GetLogExpress(EQ, Result, Varibles, VValues) result(ERR)
!--- DESCRIPTION: to get the value of a logical expression
!    INPUT:       EQ, the string of the inequations (equation)
!                 Varibles, optional, the variable characters
!                 VValues , optional,the value of the variables. The size of this 
!                 two dimensional array must agree with the size of Variables*Result 
!        
!    OUTPUT:      Result, the logical result of the inequation, the size
!                 of the array must be at lest 1
!
use STRING_OP_module, only:catenate
implicit none
 character*(*)::EQ
 integer::I, COUNT
 character*32,  dimension(:),   optional::Varibles
 double precision, dimension(:, :), optional::VValues
 logical,dimension(:)::Result

 ERR = 0 
 allocate(ExpS(len_trim(EQ)) )
 ExpS = ' '
 count = 0
 !*** to delete space and change the character to upcase
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
       VarS(I)=""
	   call catenate(VarS(I),Varibles(I))
	   VVars(I,1:count) = VValues(I,1:count)
    end do 
 end if

!*** to beging the calculation
 token = ""
 pToken  = 1
 !--- get the first inequation term in the statement 
 ERR = get_token()
 !--- get the logical value of the statement
 ERR = Get_Logical_Vaulen(result)

 deallocate(ExpS)
 if(associated(Vars)) deallocate(Vars)
 if(associated(VVars)) deallocate(VVars)

 return
end function GetLogExpress
!**************************************************************************************

!**************************************************************************************
recursive integer function Get_Logical_Vaulen(result) result(ERR)
!--- DESCRIPTION: to get the value of a logical expression defined by
!                 the module variable ExpS
!    INPUT:        
!        
!    OUTPUT:      Result, the logical result of ExpS with VarS and their value VVars

implicit none
logical, dimension(:)::result
logical, dimension(:),allocatable::hold
character*2::OP

    !--- get the vaule of the first term in the statement
    ERR = getTermValue(result)
	allocate(hold(size(result)) )

    !--- to perform logical calculation with the 
    !    second (or more) terms
    do while(.true.)
	   OP = token(1:2)
	   if( OP.eq.'&&' .or. OP.eq.'||') then
	       ERR = get_token()
           ERR = getTermValue(hold)
	       call arith(OP, result, hold)
       else
	       exit
       end if
    end do
	deallocate(hold)
    return
end function Get_Logical_Vaulen
!**************************************************************************************

!**************************************************************************************
integer function Get_token() result(ERR)
implicit none
integer::I,I1
         
	ERR = 0     
    I = pToken 
    if(I .gt. size(Exps)) then
      tok_type = DELIMITER
	  ERR = -1
	  return
    end if

	if(isDelim(Exps, I) )then
	   tok_type = DELIMITER
       token    = Exps(I)
	   if(Exps(I) .ne. Exps(I+1)) then
	      write(*,*)'Syntax Error in EXPRESS_LOGICAL_MODULE'
       end if
       I=I+1
	   token(2:2)=Exps(I)
	   pToken = I+1
    else
       I1 = 0
	   token = ""
       do while(.not.isdelim(Exps,I+I1) )
		       if(I1.lt.len(token)) then
				  token(I1+1:I1+1) = Exps(I+I1) 
                else
				   write(*,*)'Syntax Error in Express_LOGICAL_Module'
				   ERR = 1
				   return
                end if
				I1 = I1 + 1					    
       end do
	   tok_type = TERM
       pToken = I+I1
     end if	    
      
	 return
end function Get_token
!**************************************************************************************


!**************************************************************************************
integer function getTermValue(Y) result(ERR)
!--- DESCRIPTION: to get the logical value of a term in the statement.
!                 the term string, token, is determined in Get_token
!                 
!    INPUT:        
!        
!    OUTPUT:      

use Express_Module
implicit none
logical,dimension(:)::Y
character*256::T
   
   if(tok_type .ne. TERM) then
      return
   end if 
   !T = ""
   !T(1:1)= token(1:1)
   !T(2:2)= token(2:2)
   !T(3:3)= token(3:3)
   !T(4:4)= token(4:4)
   !T(5:5)= token(5:5)
   !if(T .eq. '.NOT.') then
   !   T=token(6:len_trim(token))
   !   call isComp_true(token,Y)
   !   Y = .not.Y   
   !else
   !   T=token
      call isComp_true(token,Y)
   !end if

   ERR = Get_token()

   return
end function getTermValue
!**************************************************************************************

!**************************************************************************************
subroutine isComp_true(term,Y)
use Express_Module
implicit none
character*(*)::term
character*2  ::Op
integer::I,I0
real*8,dimension(:),allocatable::number1,number2
logical,dimension(:)::Y
character*256::left,right
     
	 Op=""
	 do I=1,len_trim(term)
	    if(isCompSymbol(term(I:I)) ) then
		  if(I.eq.len(term)) then
		     stop 'Syntax Error in Express_LOGICAL_Module'
          end if
		  Op(1:1) = term(I:I)
		  I0=1
	      if(isCompSymbol(term(I+1:I+1)) )then
		     Op(2:2) = term(I+1:I+1)
			 I0 = 2
          end if
		  exit
         end if
      end do 

      left = term(1:I-1)
	  right = term(I+I0:len(term))

      allocate(number1(1:dim),number2(1:dim))
      I0 = GetExpress(left, number1, Vars, VVars) 
      I0 = GetExpress(right, number2, Vars, VVars) 
      call Comp_number(number1, number2, Op, Y)
      deallocate(number1,number2)

	 return
end subroutine isComp_true
!**************************************************************************************

!**************************************************************************************
subroutine Comp_number(number1, number2, op, Y)
implicit none
character*2  ::Op
real*8,dimension(:)::number1, number2
logical,dimension(:)::Y
!--- Local variables
integer::I
     
     select case(Op)
	        case('==')
                 DO I=1, SIZE(Y)
			        if(number1(I).eq.number2(I)) then 
				       Y(I) = .true.
                    else 
				       Y(I) = .false.
                    end if
                 END DO  

			case('> ')
                 DO I=1, SIZE(Y)
			        if(number1(I).gt.number2(I)) then
				       Y(I) = .true.
                    else 
				       Y(I) = .false.
                    end if
                 END DO

	        case('< ')
                 DO I=1, SIZE(Y)
			        if(number1(I).lt.number2(I)) then
				       Y(I) = .true.
                    else 
				       Y(I) = .false.
                    end if
                 END DO

	        case('>=', '=>')
                 DO I=1, SIZE(Y)
			        if(number1(I).ge.number2(I)) then
				       Y(I) = .true.
                    else 
				       Y(I) = .false.
                    end if
                 END DO

	        case('<=', '=<')
                 DO I=1, SIZE(Y)
			        if(number1(I).le.number2(I)) then
				       Y(I) = .true.
                    else
				       Y(I) = .false.
                    end if
                 END DO
     end select
  
	 return
end subroutine
!**************************************************************************************

!**************************************************************************************
logical function isCompSymbol(C)
implicit none
character*1::C
     
	 if(ichar(C(1:1)) .eq. ichar('>') .or.  &
	    ichar(C(1:1)) .eq. ichar('<') .or.  &
	    ichar(C(1:1)) .eq. ichar('=') ) then
		isCompSymbol = .true.
     else
		isCompSymbol = .false.
     end if
	 return
end function isCompSymbol
!**************************************************************************************

!**************************************************************************************
logical function isDelim(C, I)
implicit none
character*1, dimension(:),pointer::C
integer::I, J
   
     J = size(C)  
     if(I .gt. J) then
	    isDelim = .true.
		return
     end if 		 
	 if(C(I) .eq. '&' .or. C(I) .eq. '|') then
		isDelim = .true.
     else
		isDelim = .false.
     end if
	 return
end function isDelim
!**************************************************************************************

!**************************************************************************************
logical function isLogicalNo(C)
implicit none
character*1::C
     
	 if(C .eq. '!') then
	    isLogicalNo = .true. 
     else
	    isLogicalNo = .false. 
     end if
	 return
end function isLogicalNo
!**************************************************************************************

!**************************************************************************************
subroutine arith(OP, result, hold)
implicit none
character*2::OP
logical, dimension(:)::result, hold
integer::I
   
   select case(OP)
          case('&&')
              DO I=1, SIZE(result)
		         if(result(I).and.hold(I)) then
		            result(I) = .true.
                 else
			        result(I) = .false.
                 end if
              END DO

          case('||')
              DO I=1, SIZE(result)
   		         if(result(I).or.hold(I)) then
		            result(I) = .true.
                 else
			        result(I) = .false.
                 end if
              END DO
   end select
   return          
end subroutine arith
!**************************************************************************************

end module Express_LOGICAL_Module
