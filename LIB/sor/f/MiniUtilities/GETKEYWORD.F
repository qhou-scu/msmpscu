      subroutine getKeyword(prefix, str, kwd)
!*** to get the keyword indicated by given prefix
!**** AUTHOR: Hou Qing, updated 2014-06-03
      implicit none
      character*(*), intent(in)::str
      character*(*)::kwd
      character::prefix
      integer::I, flag1, flag2

      flag1 = 0
      flag2 = 0
      do I=1, len_trim(str)
         if(flag1 .eq. 0) then
            if(str(I:I) .eq. prefix) then
               flag1 = I
            end if
         else
            if(str(I:I) .eq. ' ')then
               exit
            end if
            flag2 = I
         end if
      end do
      if(flag2 .gt. flag1) then
         if(flag2-flag1+1 .gt. len(kwd)) flag2 = flag1+len(kwd)-1
         kwd = str(flag1:flag2)
      else
         kwd = ""
      end if
      return
      end subroutine getKeyword
!**************************************************************************************************

!**************************************************************************************************
      subroutine getNthKeyword(str,nth,kwd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract the n-th substring, seperated by space or comma
!  Input:   STRING, a string
!           n,      the rank-number of substring wanted
!  Ouput:   SUBSTR, the substring founded
!

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      character*(*), intent(in)::str
      character*(*):: kwd
      integer::i,j
      integer::n,nth
      logical::key
      i   = 1
	  n   = 0
	  j   = 0
      kwd = ""
      key = 0
      do while(.true.)
         if(str(i:i) .eq. ', ' .or. str(i:i) .eq. ' ' ) then
            key = 0  
         else
            if (key.eq.0)then
                key = 1
                n = n+1
		            if(n.gt.nth)exit
            end if
		       if (n.eq.nth)then
		          j = j+1
              kwd(j:j) = str(i:i)		
		       end if  	  
         end if
 
         i = i+1
         if(i .gt. len_trim(str) .or. j.ge. len(kwd)) then
            exit
         end if
      end do

      return
      end subroutine getNthKeyword
!**************************************************************************************************

!**************************************************************************************************
      subroutine getFirstword(str, kwd)
!*** to get the first word in a string
!**** AUTHOR: Hou Qing, updated 2014-06-03
      implicit none
      character*(*), intent(in)::str
      character*(*)::kwd
      integer::I, flag1, flag2

      flag1 = 1
      flag2 = 0
      do I=1, len_trim(str)
         if(str(I:I) .eq. ' ' .or. str(I:I) .eq. ',' .or.
     -      str(I:I) .eq. ';') then
            exit
         end if
         flag2 = I
      end do
      kwd = str(flag1:flag2)
      kwd = adjustl(kwd)
      return
      end subroutine getFirstword
!**************************************************************************************************

!**************************************************************************************************
      logical function SameFirstword(str1, str2) result(yes)
!*** to get the first word in a string
!**** AUTHOR: Hou Qing, updated 2014-06-03
      implicit none
      character*(*), intent(in)::str1,str2
      integer::I, flag1(2), flag2(2), nc, ic1, ic2, i1, i2

      flag1 = 0
      flag2 = 0
      do I=1, len_trim(str1)
         if(flag1(1) .eq. 0) then
             !--- skip the space
             if(str1(I:I) .eq. ' ') cycle
             flag1(1) = I
             flag1(2) = I
         else
           if(str1(I:I) .eq. ' ' .or.
     -        str1(I:I) .eq. ',' .or.
     -        str1(I:I) .eq. ';'  ) exit
           flag1(2) = I
         end if
      end do

      do I=1, len_trim(str2)
         if(flag2(1) .eq. 0) then
             !--- skip the space
             if(str2(I:I) .eq. ' ') cycle
             flag2(1) = I
             flag2(2) = I
         else
             if(str2(I:I) .eq. ' ' .or.
     -          str2(I:I) .eq. ',' .or.
     -          str2(I:I) .eq. ';') exit
             flag2(2) = I
         end if
      end do

      if(flag2(2)-flag2(1) .ne. flag1(2)-flag1(1)) then
         yes = .false.
         return
      end if

      nc = flag2(2)-flag2(1)+1
      i1 = flag1(1)
      i2 = flag2(1)
      yes = .true.
      do I=1, nc
         ic1 = iachar(str1(i1:i1))
         ic2 = iachar(str2(i2:i2))
         !--- the make case insensitive
         if(ic1.ge.iachar('a').and.ic1.le.iachar('z')) then
              ic1 = ic1 + iachar('A') - iachar('a')
         end if

         if(ic2.ge.iachar('a').and.ic2.le.iachar('z')) then
              ic2 = ic2 + iachar('A') - iachar('a')
         end if
         if(ic1 .ne. ic2) then
            yes =.false.
            return
         end if
         i1  = i1 + 1
         i2  = i2 + 1
      end do
      return
      end function SameFirstword


