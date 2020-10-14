      subroutine ONWARNING(action)
!*** to used to response wanring message
!**** AUTHOR: Hou Qing, updated 2014-06-03
      implicit none
      integer::action
      character*1::CC

      select case(action)
         case (1)
            write(*,fmt="(A)")'Pause on warning'
            write(*,fmt="(A)")'   Input C to continue without pause '//
     -                        'on later warning message'
            write(*,fmt="(A)")'   Input S to stop '//
     -                        'on later warning message'
            write(*,fmt="(A)")'   Or press Enter to continue '//
     -                        'with pause on warning message'

            read(*,fmt="(A1)") CC
            if(CC .eq. 'C' .or. CC.eq. 'c') then
               action = 0
            end if
            if(CC .eq. 'S' .or. CC.eq. 's') then
               action = 2
               stop
            end if
        case (2)
            stop "Process stop on warning"
      end select
      return
      end subroutine ONWARNING
