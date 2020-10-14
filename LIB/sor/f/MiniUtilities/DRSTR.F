c         program test
c        implicit real*8 (a-h,o-z)
c         write(6,*) drstr('sds1224323a, 12sdsadsa345.6768')
c         stop
c         end


      double precision function drstr(Text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This function is set to convert a text into a real number !
!  with double precision
!i Auther: Hou Qing, Inst. of Nucl. Sci. and Tech,
!          Sichuan Union University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision  tem
      character *(*) Text

       read(Text,*, err=100)tem
       drstr = tem
       return
  100  write(*, fmt="(A)") "Error when convert string'"//
     -                     Text(1:len_trim(Text))//
     -                     "'into real number"
       write(*, fmt="(A)") "Process to be stpped"
       stop
       end function


      double precision function drstr1(Text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This function is set to convert a text into a real number !
!  with double precision
!i Auther: Hou Qing, Inst. of Nucl. Sci. and Tech,
!          Sichuan Union University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer Text_Len, Count, pointm,numb,flag,sign
      double precision  tem
      character *(*) Text

       Text_Len=len(Text)
       count=0
       flag=1
       pointm=0
       sign=1
       tem=0.0

       do I=1,Text_Len
          Count=Count+1
          if((iachar(Text(Count:Count)).gt.47).and.
     -           (iachar(Text(Count:Count)).lt.58) ) then
            numb=iachar(Text(Count:Count))-48
            if(flag.gt.0) then
               tem=tem*10+numb
            else
               pointm=pointm+1
               tem=tem+dble(numb)/10**pointm
            endif
          else
            if(iachar(Text(Count:Count)).eq.46) flag=-1
            if(iachar(Text(1:1)).eq.45) sign=-1
          endif
       enddo
       drstr1=sign*tem
       return
       end function
