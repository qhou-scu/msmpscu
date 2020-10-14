  !****************************************************************************
  recursive subroutine Add_0(TheList, TheData)
  !***  PURPOSE:   to add a Mat to the list
  !
  !     INPUT:    TheList,    the data list
  !               TheData,    the data to be added to the list
  !
  !     OUTPUT    TheList
  !
  implicit none
     !--- dummy varioables
     type(LIISTNAME)   ::TheList
     DATATYPE, pointer ::TheData

  !----
       if(.not.associated(TheList%TheData)) then
          TheList%TheData => TheData
          TheList%Next    => null()
       else
          if(.not. associated(TheList%Next) ) then
             allocate(TheList%Next) 
             TheList%Next%TheData => null()
             TheList%Next%Next    => null()
          end if  
          call Add_0(TheList%Next, TheData)
       end if
       return
  end subroutine Add_0
  !****************************************************************************

  !****************************************************************************
  recursive subroutine Add_1(TheList, Tag, TheData)
  !***  PURPOSE:   to add a Mat to the list
  !
  !     INPUT:    TheList,    the data list
  !               TheData,    the data to be added to the list
  !
  !     OUTPUT    TheList
  !
  implicit none
     !--- dummy varioables
     type(LIISTNAME)  ::TheList
     character*(*)    ::Tag
     DATATYPE, pointer::TheData

  !----
       if(.not.associated(TheList%TheData)) then
          TheList%Tag  = Tag
          TheList%TheData => TheData
          TheList%Next    => null()
       else
          if(.not. associated(TheList%Next) ) then
             allocate(TheList%Next) 
             TheList%Next%TheData => null()
             TheList%Next%Next    => null()
          end if  
          call Add_1(TheList%Next, Tag, TheData)
       end if
       return
  end subroutine Add_1
  !****************************************************************************

  !****************************************************************************
  recursive subroutine Numberof_0(TheList, N)
  !***  PURPOSE:   to get the number of data on the list
  !     INPUT:    TheList,    the data list
  !     OUTPUT    N
  !
    implicit none
     !--- dummy varioables
     type(LIISTNAME)::TheList
     integer        ::N
     !--- Local
     integer::NN 
  !----
       NN = 0
       if(.not. associated(TheList%TheData)) then
          N = NN
          return
       end if
       
       if(associated(TheList%Next)) then
          call Numberof_0(TheList%Next, NN)
       end if   
       N = NN + 1
       return
  end subroutine Numberof_0
!****************************************************************************
!****************************************************************************
 recursive subroutine Get_0(TheList, Index, TheData)
 !***  PURPOSE:   to get the data with Tag
 !     INPUT:    TheList,    the data list
 !               Index,        
 !     OUTPUT    TheData
 !
   implicit none
    !--- dummy varioables
     type(LIISTNAME)  ::TheList
     integer          ::Index
     DATATYPE, pointer::TheData
     !----
       integer::I
       I = Index - 1 
       if(I .eq. 0)  then
          TheData => TheList%TheData
          return
       else   
          if(associated(TheList%Next)) then
             call Get_0(TheList%Next, I, TheData)
           end if   
       end if

       return
  end subroutine Get_0
  !****************************************************************************

  !****************************************************************************
  recursive subroutine Get_1(TheList, Tag, TheData)
  !***  PURPOSE:   to get the data with Tag
  !     INPUT:    TheList,    the data list
  !               Index,        
  !     OUTPUT    TheData
  !
   implicit none
    !--- dummy varioables
     type(LIISTNAME)  ::TheList
     character*(*)    ::Tag
     DATATYPE, pointer::TheData
     !----
       if(TheList%Tag .eq. Tag)  then
          TheData => TheList%TheData
          return
       else   
          if(associated(TheList%Next)) then
             call Get_1(TheList%Next, Tag, TheData)
           end if   
       end if
       return
  end subroutine Get_1
  !****************************************************************************

  !****************************************************************************
 recursive subroutine Clear_0(TheList)
 !***  PURPOSE:   to clear the data added in the list
 !     INPUT:    TheList,    the data list
 !     OUTPUT    TheData
 !
   implicit none
    !--- dummy varioables
     type(LIISTNAME)  ::TheList

     !----

       if(associated(TheList%Next)) then
         call Clear_0(TheList)
         deallocate(TheList%Next)
         TheList%Next => null()
       end if
       
       if(associated(TheList%TheData)) then
          #ifdef DEALLOCTE_DATA
          call DEALLOCTE_DATA(TheList%TheData)
          #endif
          TheList%TheData => null()
       end if   
       return
  end subroutine Clear_0
  !****************************************************************************

