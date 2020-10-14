        subroutine extract_optstr( String,SEP,KWD,mxCount,count,substr,
     -                             ind)
!**************************************************************************
!  Purpose: to extract string following a keyword.
!
!  Input:   STRING, a string
!           KWD ,   the keyword
!           SEP,    the seperator between keyword
!
!  Ouput:   COUNT  the number of substring found in the string
!           SUBSTR the substring founded
!
! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
!
!
       implicit none
       integer::count, mxCount
       character *(*)::String, KWD, SEP
       character *(*)::substr
       dimension substr(mxCount)
       integer, optional::ind

       !Local
       integer::I, I0, I1, FLAG
       character*256::tstr, tstr0
       character*16::tkwd

        tstr0 = String(1:len_trim(String))
        tstr  = tstr0

        call UPCASE(tstr)
        tkwd = SEP(1:1)//KWD(1:len_trim(KWD))
        call UPCASE(tkwd)

        count = 0
        flag  = 0
        do while(.true.)
           tstr = adjustl(tstr)
           !--- to find out the position of the kwd
           I0 = index(tstr, tkwd(1:len_trim(tkwd)) )
           if(I0 .eq. 0) exit
           I1 = I0 + len_trim(tkwd)
           !--- check if the keyword followed by space or ','
           if(iachar(tstr(I1:I1)) .ne. 32 .and.
     -        iachar(tstr(I1:I1)) .ne. 44) then
              !--- not the keyword
              tstr  = tstr(I1:len_trim(tstr))
              tstr0 = tstr0(I1:len_trim(tstr0))
              cycle
           end if
           flag = 1
           exit
        end do
        if(.not.flag) return

           !--- find the keyword, extract the substring
           if(present(ind)) ind = I0

           tstr  = tstr(I1+1:len_trim(tstr))
           tstr0 = tstr0(I1+1:len_trim(tstr0))

           do while(.true.)
              if(count .ge. mxCount) exit
              !--- check if it is a seperator of keywords
              tstr  = adjustl(tstr)
              tstr0 = adjustl(tstr0)
              if(tstr(1:1) .eq. SEP(1:1)) exit

              I0 = 1
              do I=1, len_trim(tstr)+1
                 I1 = I
                 if(iachar(tstr(I:I)) .eq. 32 .or.
     -              iachar(tstr(I:I)) .eq. 44) then
                    exit
                 end if
              end do
              count = count + 1
              substr(count) = tstr0(I0:min(I1,len(substr(count))) )

              if(iachar(tstr(I:I)) .eq. 44) I1 = I1+1
              tstr = tstr(I1:len_trim(tstr))
              tstr = adjustl(tstr)
              tstr0 = tstr0(I1:len_trim(tstr0))
              tstr0 = adjustl(tstr0)
              if(len_trim(tstr).le. 0) exit
           end do

       return
       end subroutine extract_optstr

       subroutine extract_optstr1( String,SEP,KWD,substr)
    !**************************************************************************
    !  Purpose: to extract string following a keyword.
    !
    !  Input:   STRING, a string
    !           KWD ,   the keyword
    !           SEP,    the seperator between keyword
    !
    !  Ouput:   COUNT  the number of substring found in the string
    !           SUBSTR the substring founded
    !
    ! Auther:   Hou Qing, Inst. of Nucl. Sci.and Tech., Sichuan Union University
    !
    !
           implicit none
           character *(*)::String, KWD, SEP
           character *(*)::substr
    
           !Local
           integer::I, I0, I1, FLAG
           character*256::tstr, tstr0
           character*16::tkwd
    
            tstr0 = String(1:len_trim(String))
            tstr  = tstr0
            call UPCASE(tstr)
            tkwd = SEP(1:1)//KWD(1:len_trim(KWD))
            call UPCASE(tkwd)
    
            flag  = 0
            do while(.true.)
               tstr = adjustl(tstr)
               !--- to find out the position of the kwd
               I0 = index(tstr, tkwd(1:len_trim(tkwd)) )
               if(I0 .eq. 0) exit
               I1 = I0 + len_trim(tkwd)
               !--- check if the keyword followed by space or ','
               if(iachar(tstr(I1:I1)) .ne. 32 .and.    
     -            iachar(tstr(I1:I1)) .ne. 44) then
                  !--- not the keyword
                  tstr  = tstr(I1:len_trim(tstr))
                  tstr0 = tstr0(I1:len_trim(tstr0))
                  cycle
               end if
               flag = 1
               exit
            end do
            if(.not.flag) return
    
            !--- extract the substring
             tstr  = tstr(I1+1:len_trim(tstr))
             tstr0 = tstr0(I1+1:len_trim(tstr0))
             do I=1, len_trim(tstr)
                if(tstr(I:I) .eq. SEP(1:1)) then
                   exit
                end if    
                substr(I:I) = tstr0(I:I)
             end do 
           return
           end subroutine extract_optstr1
    
