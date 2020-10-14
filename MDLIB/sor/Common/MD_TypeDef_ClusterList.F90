  module MD_TYPEDEF_ClusterList
  !***  DESCRIPTION:
  !     This module provides a typedef for the connected particles or clusters.
  !     A cluster is defined as a collection of IDs of atoms that connected.
  !     The positions of the atoms are not included in the cluster object.
  !     The position can not be get from the atoms directed by the IDs.
  !     Thus, a cluster object is usually used along with a SimulatioBox.
  !
  !    Written by HOU Qing, Nov., 2015
  !
   use MD_CONSTANTS
   implicit none

      type PartcleClusterList
           integer::NPRT = 0                             !--- the number of particles in the cluster
           integer, dimension(:), allocatable::PTYP      !--- the type of the particles (atom, or, point defect)
           integer, dimension(:), allocatable::IND       !--- the ID of the atoms
           type(PartcleClusterList), pointer::prev=>null()
           type(PartcleClusterList), pointer::next=>null()
      end type PartcleClusterList

  contains

 !****************************************************************************
  recursive subroutine Add_PartcleClusters(Clusters, NPRT, IND, PTYP)
  !***  PURPOSE:   to add a cluster to the list
  !
  !     INPUT:     NPRT,      the number of particles
  !                IND,       the index of the atoms
  !                PTYPE,     the type of the atoms (defect)
  !
  !     OUTPUT     Clsuter,   the PartcleClusterList
  !
  implicit none
     !--- dummy varioables
     type(PartcleClusterList),target::Clusters
     integer, intent(in)::NPRT
     integer, dimension(:), intent(in)::IND
     integer, dimension(:), intent(in)::PTYP
     !----
              if(NPRT .le. 0) return

              if(Clusters%NPRT .eq. 0) then
                 Clusters%NPRT = NPRT
                 if(allocated(Clusters%IND) ) deallocate(Clusters%IND)
                 if(allocated(Clusters%PTYP)) deallocate(Clusters%PTYP)
                 allocate(Clusters%IND(NPRT), Clusters%PTYP(NPRT))
                 Clusters%IND(1:NPRT)  = IND(1:NPRT)
                 Clusters%PTYP(1:NPRT) = PTYP(1:NPRT)
                 Clusters%next =>null()
                 Clusters%prev =>null()
                 return
              else
                 if(.not. associated(Clusters%next) ) then
                    allocate(Clusters%next)
                    Clusters%next%NPRT = 0
                    call Add_PartcleClusters(Clusters%next, NPRT, IND, PTYP)
                    Clusters%next%prev => Clusters
                 else
                    call Add_PartcleClusters(Clusters%next, NPRT, IND, PTYP)
                 end if
              end if

              return
  end subroutine Add_PartcleClusters
 !****************************************************************************

 !****************************************************************************
  subroutine Expand_PartcleClusters(Cluster, NPRT, IND, PTYP)
  !***  PURPOSE:   to expand a cluster by adding more particles
  !
  !     INPUT:     NPRT,      the number of new particles
  !                IND,       the index of the new particles
  !                PTYPE,     the type of the new particles
  !
  !     OUTPUT     Clsuter,   the Cluster
  !
  !     NOTE:      this routine operating on one cluster
  !
  implicit none
     !--- dummy varioables
     type(PartcleClusterList)::Cluster
     integer, intent(in)::NPRT
     integer, dimension(:), intent(in)::IND
     integer, dimension(:), intent(in)::PTYP
     !---- Local variables
     integer::newNPRT, I, NP
     integer, dimension(:), allocatable::newIND, newPTYP

             if(NPRT .le. 0) return

              if(Cluster%NPRT .eq. 0) then
                 call Add_PartcleClusters(Cluster, NPRT, IND, PTYP)
                 return
              end if

                 newNPRT = Cluster%NPRT+NPRT
                 allocate(newIND(newNPRT),  newPTYP(newNPRT))
                 newIND(1:Cluster%NPRT)  = Cluster%IND(1:Cluster%NPRT)
                 newPTYP(1:Cluster%NPRT) = Cluster%PTYP(1:Cluster%NPRT)

                 NP = Cluster%NPRT
                 do I=1, NPRT
                    if(any(newIND(1:NP) .ne. IND(I)) ) then
                      NP = NP + 1
                      newIND(NP)  = IND(I)
                      newPTYP(NP) = PTYP(I)
                    end if
                 end do

                 deallocate(Cluster%IND, Cluster%PTYP)
                 allocate(Cluster%IND(NP), Cluster%PTYP(NP))

                 Cluster%NPRT       = NP
                 Cluster%IND(1:NP)  = newIND(1:NP)
                 Cluster%PTYP(1:NP) = newPTYP(1:NP)

              deallocate(newIND, newPTYP)
              return
  end subroutine Expand_PartcleClusters
 !****************************************************************************

 !****************************************************************************
  subroutine NumberCluster_PartcleClusters(Clusters, NUM)
  !***  PURPOSE:   to get the number of clusters on the list
  !
  !     INPUT:      Clsuters,   the cluster list
  !     OUTPUT:     NUM,        the number of clusters on the list
  !
  implicit none
     !--- dummy varioables
     type(PartcleClusterList)::Clusters
     integer::NUM
     !---- local variables
      type(PartcleClusterList), pointer::tC

           if(Clusters%NPRT .eq. 0) then
              NUM = 0
              return
           end if

           NUM = 1
           tC =>Clusters%next
           do while(associated(tC))
              NUM = NUM + 1
              tC => tC%next
           end do
           return
  end subroutine NumberCluster_PartcleClusters
 !****************************************************************************

!****************************************************************************
  subroutine NumberPart_PartcleClusters(Clusters, NUM, TYP)
  !***  PURPOSE:   to get the number of particles on the list
  !
  !     INPUT:      Clsuters,   the cluster list
  !                 TYP,        optional,the type of interesting particle
  !
  !     OUTPUT:     NUM,        the number of clusters on the list
  !
  implicit none
     !--- dummy varioables
     type(PartcleClusterList)::Clusters
     integer::NUM
     integer, optional::TYP
     !---- local variables
     integer::N
      type(PartcleClusterList), pointer::tC

           if(Clusters%NPRT .eq. 0) then
              NUM = 0
              return
           end if

           if(present(TYP)) then
              NUM = count(Clusters%PTYP(1:Clusters%NPRT) .eq. TYP)
              tC =>Clusters%next
              do while(associated(tC))
                 NUM = NUM + count(tC%PTYP(1:Clusters%NPRT) .eq. TYP)
                 tC => tC%next
              end do
           else
             NUM = Clusters%NPRT
             tC =>Clusters%next
             do while(associated(tC))
                NUM = NUM + tC%NPRT
                tC => tC%next
             end do
           end if

           return
  end subroutine NumberPart_PartcleClusters
 !****************************************************************************

 !****************************************************************************
  subroutine GetNext_PartcleClusters(ClustList, ITH, Cluster)
  !***  PURPOSE:   to get the pointer of ith clusters
  !
  !     INPUT:      List,      the cluster list
  !                 ITH,       position of the clusters on the list
  !     OUTPUT:     Cluster,   the pointer of cluster
  !
  implicit none
     !--- dummy varioables
     integer::ITH
     type(PartcleClusterList), target::ClustList
     type(PartcleClusterList), pointer::Cluster
     !---- local variables
      integer::I
      type(PartcleClusterList), pointer::tC

            tC => ClustList
            do I=1, ITH
               tC => tC%next
            end do
            Cluster => tC
           return
  end subroutine GetNext_PartcleClusters
 !****************************************************************************

 !****************************************************************************
  subroutine Delete_PartcleClusters(ClustList, ITH)
  !***  PURPOSE:   to delete the ith cluster on the list
  !
  !     INPUT:      ClustList,   the cluster list
  !                 ITH,         position of the clusters on the list
  !     OUTPUT:     ClustList,
  !
  implicit none
     !--- dummy varioables
     integer::ITH
     type(PartcleClusterList), target::ClustList
     !---- local variables
      integer::I
      type(PartcleClusterList), pointer::todel, next, prev

           todel => ClustList
           if(ITH .eq. 0) then
              if(allocated(todel%IND))   deallocate(todel%IND)
              if(allocated(todel%PTYP))  deallocate(todel%PTYP)
              todel%NPRT = 0
              next => todel%next
              if(associated(next)) then
                 call Add_PartcleClusters(todel, next%NPRT, next%IND, next%PTYP)
                 todel%next => next%next

                 if(allocated(next%IND))   deallocate(next%IND)
                 if(allocated(next%PTYP))  deallocate(next%PTYP)
                 next%NPRT = 0
                 deallocate(next)
              end if
              return
           end if

           !$$--- for ITH >= 1
           do I=1, ITH
               todel => todel%next
           end do
           if(.not.associated(todel)) return

           prev => todel%prev
           next => todel%next
           prev%next => next
           if(associated(next)) then
               next%prev => prev
           end if

           if(allocated(todel%IND))   deallocate(todel%IND)
           if(allocated(todel%PTYP))  deallocate(todel%PTYP)
           todel%NPRT = 0
           deallocate(todel)
           return
  end subroutine Delete_PartcleClusters
 !****************************************************************************

 !****************************************************************************
  subroutine DeleteAll_PartcleClusters(Clusters)
  !***  PURPOSE:   to release the list
  !
  !     INPUT:      Clsuter,   the cluster to be release
  !     OUTPUT:     Clsuter,   the cluster released
  !
  implicit none
     !--- dummy varioables
     type(PartcleClusterList)::Clusters
     !----
      type(PartcleClusterList), pointer::todel, next

            next=>Clusters%next
            if(allocated(Clusters%IND))   deallocate(Clusters%IND)
            if(allocated(Clusters%PTYP))  deallocate(Clusters%PTYP)
            Clusters%NPRT = 0
            Clusters%prev => null()
            Clusters%next => null()
            todel => next


            do while(associated(todel))
               next => todel%next
               if(allocated(todel%IND))   deallocate(todel%IND)
               if(allocated(todel%PTYP))  deallocate(todel%PTYP)

               todel%NPRT = 0
               todel%prev => null()
               todel%next => null()
               deallocate(todel)
               todel =>next
            end do

           return
  end subroutine DeleteAll_PartcleClusters
 !****************************************************************************

  !*****************************************************************************
  subroutine CopyClusterList_PartcleClusters(ClusterS, ClusterT)
  !***  DESCRIPTION: to copy a cluster (S) to a a target cluster
  !
  !    INPUT:        ClusterS, the source cluster
  !                  ClusterT, the target cluster
  !
  !    OUTPUT:      ClusterT, the target cluster
  !
  implicit none
       !--- dummy variables
       type(PartcleClusterList)::ClusterS
       type(PartcleClusterList)::ClusterT
       !--- local variables
        type(PartcleClusterList), pointer::tC

             call DeleteAll_PartcleClusters(ClusterT)
             call Add_PartcleClusters(ClusterT, ClusterS%NPRT, ClusterS%IND, ClusterS%PTYP)

             tC => ClusterS%next
             do while(associated(tC) )
                call Add_PartcleClusters(ClusterT, tC%NPRT, tC%IND, tC%PTYP)
                tC => tC%next
             end do

             return
  end subroutine CopyClusterList_PartcleClusters
  !**********************************************************************************

  !**********************************************************************************
  subroutine CheckConnected_PartcleClusters(Cluster1, Cluster2, Stat, NN, NLIST)
  !***  DESCRIPTION: to check if two clusters are connected. Two cluster are defined
  !                  as connected if they share at least one same particle ID
  !                  OR
  !                  they have at least two particles are in neast neighbore of each other
  !
  !    INPUT:        Cluster1, Cluter2, the clusters to be check
  !                  NN, NLIST, optional, the neighbor-list of particles.
  !                             if NN or NLIST is not supplied, if two clusters is checked
  !                             by their shared particles.
  !                             otherwise, the the nearest neighbores are checked.
  !
  !    OUTPUT:     Stat,  =0   not connected
  !                       =1   share the same particles
  !                       =2   attached, nearest neighbor
  !
  implicit none
       !--- dummy variables
       type(PartcleClusterList), intent(in)::Cluster1
       type(PartcleClusterList), intent(in)::Cluster2
       integer, dimension(:),   optional::NN
       integer, dimension(:,:), optional::NLIST
       integer::STAT
       !--- local variables
        integer::I, J, IA2, NN2

            STAT = 0
            do I=1, Cluster2%NPRT
               if(any( Cluster1%IND(1:Cluster1%NPRT) .eq. Cluster2%IND(I) )  ) then
                  STAT = 1
                  exit
               end if
             end do

             if(STAT .gt. 0) return

             !$$--- check if attached
             if(present(NLIST) .and. present(NN)) then
                do I=1, Cluster2%NPRT
                   IA2 = Cluster2%IND(I)
                   do J=1, NN(IA2)
                      if(any(Cluster1%IND(1:Cluster1%NPRT) .eq. NLIST(IA2, J) )  ) then
                         STAT = 2
                         return
                    end if
                   end do
                end do
             end if

             return
  end subroutine CheckConnected_PartcleClusters
  !**********************************************************************************

  !*****************************************************************************
  logical function MergeTwoClusters_PartcleClusters(ClusterS, ClusterT, NN, NLIST) result(YES)
  !***  DESCRIPTION: to merge cluster (S) to Cluster T if they are connected.
  !                  Two cluster are defined as connected if they share at
  !                  least one same particle ID.
  !
  !    INPUT:        ClusterS, the source cluster
  !                  ClusterT, the target cluster
  !                  NN, NLIST, optional, the neighbor-list of particles.
  !                             if NN or NLIST is not supplied, if two clusters is checked
  !                             by their shared particles.
  !                             otherwise, the the nearest neighbores are checked.
  !
  !    OUTPUT:      CluterT,   the cluster with ClusterS mereged to.
  !                 YES,       = .true. if merged
  !                            = .false. if not merged
  !
  implicit none
       !--- dummy variables
       type(PartcleClusterList)::ClusterS
       type(PartcleClusterList)::ClusterT
       integer, dimension(:),   optional::NN
       integer, dimension(:,:), optional::NLIST
       !--- local variables
       integer::Connect


            if(present(NN) .and. present(NLIST)) then
              call CheckConnected_PartcleClusters(ClusterT, ClusterS, Connect, NN, NLIST)
            else
              call CheckConnected_PartcleClusters(ClusterT, ClusterS, Connect)
            end if

            if(Connect .eq. 0) then
               YES = .false.
               return
            end if

            call Expand_PartcleClusters(ClusterT, ClusterS%NPRT, ClusterS%IND, ClusterS%PTYP)
            YES = .true.

            return
  end function MergeTwoClusters_PartcleClusters
  !**********************************************************************************

  !**********************************************************************************
  subroutine MergeClusterList_PartcleClusters(ClusterList, NN, NLIST)
  !***  DESCRIPTION: to merge the clusters on list that are connected
  !                  Two cluster are defined as connected if they share at
  !                  least one same particle ID.
  !                  NN, NLIST, optional, the neighbor-list of particles.
  !                             if NN or NLIST is not supplied, if two clusters is checked
  !                             by their shared particles.
  !                             otherwise, the the nearest neighbores are checked.
  !
  !    INPUT:  ClusterList,   the cluster list
  !
  !
  !    OUTPUT: ClusterList,   the cluster list with the connected clusters merged
  !
  implicit none
       !--- dummy variables
       type(PartcleClusterList),target::ClusterList
       integer, dimension(:),   optional::NN
       integer, dimension(:,:), optional::NLIST
       !--- local variables
        integer::FLAG
        logical::MERGED
        type(PartcleClusterList), pointer::tC1, tC2, tC3


             tC1 => ClusterList
             do while(.true.)
                FLAG = 1
                do while(.true.)
                   call GetNext_PartcleClusters(tC1, FLAG, tC2)

                   if(.not. associated(tC2) ) then
                      exit
                   end if

                   if(present(NN) .and. present(NLIST) ) then
                      MERGED = MergeTwoClusters_PartcleClusters(tC2, tC1, NN, NLIST)
                   else
                      MERGED = MergeTwoClusters_PartcleClusters(tC2, tC1)
                   end if

                   if( MERGED ) then
                       call Delete_PartcleClusters(tC1, FLAG)
                   else
                       FLAG = FLAG + 1
                   end if
                end do

               !$$--- move to the next cluster
                tC1 => tC1%next
                if(.not.associated(tC1)) then
                   exit
                end if
             end do

             return
  end subroutine MergeClusterList_PartcleClusters
  !**********************************************************************************

  !**********************************************************************************
  subroutine Array2ClusterList_PartcleClusters(SQSIZE, SQNC, SQTYP, ClusterList)
  !***  DESCRIPTION: to transfer a cluster sequence recoded in array to ClusterList.
  !                  a cluster sequence is an array with each element
  !                  is the particle ID. A zero element indicates the end
  !                  of a cluster
  !
  !    INPUT:  SQSIZE,   the size of the sequence
  !            SQNC,     the particle sequence
  !            SQTYP,    the type of the particles
  !
  !    OUTPUT: ClusterList,   the cluster list
  !
  implicit none
       !--- dummy variables
       integer, intent(in)::SQSIZE
       integer, dimension(:), intent(in)::SQNC, SQTYP
       type(PartcleClusterList)::ClusterList
       !--- local variables
        integer::I, FLAG

             call DeleteAll_PartcleClusters(ClusterList)

             FLAG = 0
             do I=1, SQSIZE
                if(SQNC(I) .le. 0) then
                   call Add_PartcleClusters(ClusterList, I-1-FLAG, SQNC(FLAG+1:I-1), SQTYP(FLAG+1:I-1))
                   FLAG = I
                end if
             end do

             return
  end subroutine Array2ClusterList_PartcleClusters
  !**********************************************************************************

  !**********************************************************************************
  subroutine ClusterList2Array_PartcleClusters(ClusterList, SQSIZE, SQNC, SQTYP, FILTER)
  !***  DESCRIPTION: to transfer a cluster list to a cluster sequence.
  !                  A cluster sequence is an array with each element
  !                  is the particle ID. A zero element indicates the end
  !                  of a cluster
  !
  !    INPUT: ClusterList,   the cluster list
  !           FILTER,        optional, the type of particles to NOT be transfer to the array
  !
  !    OUTPUT: SQSIZE,   the size of the sequence
  !            SQNC,     the particle sequence
  !            SQTYP,    the type of the particles
  !
  implicit none
       !--- dummy variables
       type(PartcleClusterList), target, intent(in)::ClusterList
       integer::SQSIZE
       integer, dimension(:), allocatable::SQNC, SQTYP
       integer, optional::FILTER
       !--- local variables
        integer::NC, NA, FLAG, I
        type(PartcleClusterList), pointer::tC

             if(allocated(SQNC))  deallocate(SQNC)
             if(allocated(SQTYP)) deallocate(SQTYP)
             SQSIZE = 0

             call NumberCluster_PartcleClusters(ClusterList, NC)
             if(NC .le. 0) return

             call NumberPart_PartcleClusters(ClusterList, NA)
             if(NA .le. 0) return

             SQSIZE = NA + NC
             allocate(SQNC(SQSIZE), SQTYP(SQSIZE))
             SQNC  = 0
             SQTYP = 0

             tC => ClusterList
             FLAG = 1

             if(present(FILTER)) then
               do while(associated(tC))
                  do I=1, tC%NPRT
                     if(tC%PTYP(I) .ne. FILTER) then
                        SQNC(FLAG)  = tC%IND(I)
                        SQTYP(FLAG) = tC%PTYP(I)
                        FLAG = FLAG +  1
                     end if
                  end do
                  FLAG = FLAG + 1
                  tC => tC%next
                end do
             else
               do while(associated(tC))
                  SQNC(FLAG:FLAG  + tC%NPRT - 1)  = tC%IND(1:tC%NPRT)
                  SQTYP(FLAG:FLAG + tC%NPRT - 1)  = tC%PTYP(1:tC%NPRT)
                  FLAG = FLAG +  tC%NPRT + 1
                  tC => tC%next
                end do
             end if
             SQSIZE  = FLAG
             return
  end subroutine ClusterList2Array_PartcleClusters
  !**********************************************************************************


  end module MD_TYPEDEF_ClusterList
  !****************************************************************************************
