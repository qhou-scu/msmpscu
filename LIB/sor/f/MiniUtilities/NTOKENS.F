      integer function ntokens( String, token)
      integer i, n, token
      character *(*) String

        N= 0
        do I=1, len_trim(String)
           if(iachar(string(I:I)) .eq. token) then
              N = N+1
           end if
        end do
        ntokens = N

      return
      end function
