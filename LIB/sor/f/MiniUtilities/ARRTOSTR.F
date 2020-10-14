          subroutine ArrytoStr(A,str)
          implicit none
          character*(*)::str
          character*1, dimension(:)::A
          integer::I

         do I = 1, min(size(A),len(str))
          str(I:I) = A(I)
          end do
          return
          end subroutine ArrytoStr



