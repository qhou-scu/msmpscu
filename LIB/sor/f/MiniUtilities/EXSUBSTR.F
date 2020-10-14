!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
        subroutine extract_substr( String, mxCount,count,substr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract substring indicated ""
!  Input:   STRING a string
!  Ouput:   COUNT  the number of substring found in the string
!           SUBSTR the substring founded
!
! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       implicit none
       integer::i, index1, count, mxCount
       character *(*)::String
       character *(*)::substr
       dimension substr(mxCount)
       integer::curchar

       count=0
       index1=0
       i = 1
       substr = ""
       do while(.true.)
          if( iachar(String(i:i)) .eq. curchar .and. index1.gt.0 ) then
             count = count + 1
             substr(count)=string(index1:i-1)
             index1 = 0
          else
             if(String(i:i) .eq. '"' .or. String(i:i).eq."'"  )then
                index1 = i+1
                curchar = iachar(String(i:i))
             end if
          end if

          i=i+1
          if(i .gt. len_trim(string) .or. count .ge. mxCount) then
             exit
          end if
       enddo
       return
       end subroutine extract_substr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!

      subroutine extract_subsymb(String, symb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!
!  Purpose: to extract header symble of a string
!  Input:   STRING a string
!
!            example: String = He3, output Symb = He
!
!  Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       implicit none
       integer  i, j, n, flag
       character *(*)::String, Symb

       Symb = ""
       n    = len(symb)
       flag = 0
       j    = 0
       do i = 1, len_trim(String)
          !--- skip the space
          if(flag .le. 0) then
             if(String(I:I) .eq. ' ') cycle
             flag = 1
          end if

          if((iachar(String(i:i)) .ge. iachar('A')  .and.
     c        iachar(String(i:i)) .le. iachar('Z')) .or.
     c       (iachar(String(i:i)) .ge. iachar('a')  .and.
     c        iachar(String(i:i)) .le. iachar('z'))  )then
              j = j + 1
              Symb(j:j) = String(i:i)
              if(j .ge. n) exit
              cycle
           end if
           exit
       end do
      end subroutine extract_subsymb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine extract_subnum(String, Sbn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Purpose: to extract sub-string which is numberical characters
!  Input:   STRING a string
!
!            example: String = He3, output Sbn = 3
!
!  Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       implicit none
       integer  i, j, n
       character *(*)::String, Sbn

       Sbn = ""
       n   = len(Sbn)

       j = 0
       do i = 1, len_trim(String)
          if( (iachar(String(i:i)).ge.iachar('0'))  .and.
     c        (iachar(String(i:i)).le.iachar('9')) ) then
              j = j + 1
              Sbn(j:j) = String(i:i)
              if(j.ge.n) exit
              cycle
           end if
           if(j. gt. 0) then
              exit
           end if
       end do
      end subroutine extract_subnum


