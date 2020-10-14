!======================================================================
      subroutine rqsort(n,a,p)
!======================================================================
!     Return integer array p which indexes array a in increasing order.
!     Array a is not disturbed.  The Quicksort algorithm is used.
!
!     B. G. Knapp, 86/12/23
!
!     Reference: N. Wirth, Algorithms and Data Structures,
!     Prentice-Hall, 1986
!======================================================================
      implicit none

!    Input:
      integer   n
      real      a(n)

!    Output:
      integer   p(n)

!    Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
!        (LGN = log base 2 of maximum n;
!         Q = smallest subfile to use quicksort on)

!     Local:
      real      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

!     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

!     Initialize the pointer array
      do i=1,n
         p(i)=i
      end do

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

!           Use straight insertion
            do i=l+1,r
               t = p(i)
               x = a(t)
               do j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
               end do
               j=l-1
    5          p(j+1) = t
           end do
         else

!           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

!           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

!           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end subroutine rqsort
      
!======================================================================
      subroutine dqsort(n,a,p)
!======================================================================
!     Return integer array p which indexes array a in increasing order.
!     Array a is not disturbed.  The Quicksort algorithm is used.
!
!     B. G. Knapp, 86/12/23
!
!     Reference: N. Wirth, Algorithms and Data Structures,
!     Prentice-Hall, 1986
!======================================================================
      implicit none

!    Input:
      integer   n
      real*8    a(n)

!    Output:
      integer   p(n)

!    Constants
      integer   LGN, Q
      parameter (LGN=64, Q=11)
!        (LGN = log base 2 of maximum n;
!         Q = smallest subfile to use quicksort on)

!     Local:
      real*8    x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

!     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

!     Initialize the pointer array
      do i=1,n
         p(i)=i
      end do

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

!           Use straight insertion
            do i=l+1,r
               t = p(i)
               x = a(t)
               do j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
               end do
               j=l-1
    5          p(j+1) = t
           end do
         else

!           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

!           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

!           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end subroutine dqsort      