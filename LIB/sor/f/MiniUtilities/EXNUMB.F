CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c             program test
c            character *250 a
c            character *20 b
c            dimension b(20)
cc  1          read (*,*) a
c            a='this is a try 1234.6 56756.,a1,a2,a3, -98765.875'
c            call extract_numb( a, n,b)
c            do i=1,n
c              write(6,*) n,i,'  ',b(i)
c                enddo
c             if(n.eq.0) stop
C              goto 1
c              end
C
C_______________________________________________________________


      subroutine extract_numb( String, mxCount,count,real_numb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract number from a string
!  Input:   STRING a string
!  Ouput:   COUNT  the number of numbers found in the string
!           REAL_NUMB the number text founded
!  Note:    the REAL_NUMBi are character varible. These character number
!           could be transformed to real number by DRSTR or RSTRi
! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       implicit none
       integer  i, flag, index1, count, mxCount
       character *(*) String
       character *(*) real_numb
       dimension real_numb(mxCount)

       count=0
       index1=0
       i = 1
       flag   = 0
       do while(.true.)
          if(String(i:i) .eq. ' ' .or.
     -    String(i:i) .eq. ',' .or.
     -    String(i:i) .eq. '(' .or.
     -    String(i:i) .eq. ')' .or.
     -    String(i:i) .eq. '[' .or.
     -    String(i:i) .eq. ']' .or.
     -    String(i:i) .eq. '{' .or.
     -    String(i:i) .eq. '}' .or.
     -    String(i:i) .eq. '='  ) then
            if(index1 .gt. 0) then
               real_numb(count)=string(i-index1:i-1)
               index1 = 0
            end if
            flag = 0
         else
           if((iachar(String(i:i)) .ge.iachar('0') .and.
     -           iachar(String(i:i)) .le.iachar('9')) .or.
     -           iachar(String(i:i)) .eq.iachar('-') .or.
     -           iachar(String(i:i)) .eq.iachar('+') .or.
     -           iachar(String(i:i)) .eq.iachar('.') )then
                 if(flag.eq.0 ) then
                    flag = 1
                    count = count + 1
                    if(count .gt. mxCount) exit
                 end if
           else
               if(flag.ne.1) flag = 2
           end if
         end if

        if(flag.eq.1) then
           index1 = index1+1
        end if

        i=i+1
        if(i .gt. len_trim(string) ) then
           if(index1 .gt. 0) then
              real_numb(count)=string(i-index1:i-1)
           end if
           exit
        end if

      enddo
      return
      end subroutine


      subroutine extract_numb1( String, count,real_numb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract number from a string
!  Input:   STRING a string
!  Ouput:   COUNT  the number of numbers found in the string
!           REAL_NUMB the number text founded
!  Note:    the REAL_NUMBi are character varible. These character number
!           could be transformed to real number by DRSTR or RSTRi
! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer  i, flag,index1,index2, count
      character *(*) String
      character *(*) real_numb
      dimension real_numb(*)

      count=0
      flag=1
      index1=0
      index2=0
      do i=1,len( String)
         if(iachar(String(i:i)).eq.45.or.iachar(String(i:i)).eq.46.or.
     -   (iachar(String(i:i)).lt.58.and.iachar(String(i:i)).gt.47))then
          if(flag.eq.1) then
             index1=i
          !--- ASCII 32 - ' ', 44 -','
             if( iachar(String(index1-1:index1-1)).eq.32.or.
     -               iachar(String(index1-1:index1-1)).eq.44 ) then
                 flag=-1
                 count=count+1
             endif
          endif
          index2=i
         else
           if(flag.eq.-1) then
              real_numb(count)=string(index1:index2)
              flag=1
          endif
        endif
      enddo
      return
      end subroutine
