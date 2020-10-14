       integer function istr(Text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This function is set to convert a text into a integer number
!  Author: Hou Qing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer Count,numb,flag,sign
      integer tem
      character *(*) Text

        count=0
        flag=1
        sign=1
        tem=0

        do I=1,Len(Text)
           Count=Count+1
           if((iachar(Text(Count:Count)).gt.47).and.
     -        (iachar(Text(Count:Count)).lt.58) ) then
               numb=iachar(Text(Count:Count))-48
               if(flag.gt.0) tem=tem*10+numb
          else
              if(iachar(Text(Count:Count)).eq.46) flag=-1
              if(iachar(Text(1:1)).eq.45) sign=-1
         endif
       enddo
       istr=sign*tem
       return
       end function
