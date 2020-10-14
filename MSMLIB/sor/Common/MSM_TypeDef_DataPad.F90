  module MSM_TYPEDEF_DataPad
  !***  DESCRIPTION:
  !     This module provides a typedef for data operation.
  !
  !
  !
  !    Written by HOU Qing, Nov., 2015
  !
   use MSM_CONSTANTS
   use MiniUtilities
   implicit none

      integer, parameter, private::mp_TAGLEN = 64
      type::DataPad
           character(len=mp_TAGLEN)            ::Tag
           integer,     dimension(:),   pointer::IDATA1 = null()
           real(KINDDF),dimension(:),   pointer::DDATA1 = null()
           integer,     dimension(:,:), pointer::IDATA2 = null()
           real(KINDDF),dimension(:,:), pointer::DDATA2 = null()
           type(DataPad),               pointer::prev=>null()
           type(DataPad),               pointer::next=>null()
      end type DataPad

  !--- interface to the external routine -------------------
      !--------------------------------------------
      private:: AccumulateIByTag_DataPad,    &
                AccumulateDByTag_DataPad,    &
                AccumulateI2ByTag_DataPad,   &
                AccumulateD2ByTag_DataPad
      interface Accumulate_DataPad
             module procedure AccumulateIByTag_DataPad
             module procedure AccumulateDByTag_DataPad
             module procedure AccumulateI2ByTag_DataPad
             module procedure AccumulateD2ByTag_DataPad
      end interface Accumulate_DataPad
      !--------------------------------------------
      public::  Archive_DataPad
      !--------------------------------------------
      private:: ChangeSize_DataPad0,       &
                ChangeSizeByInd_DataPad,   &
                ChangeSizeByTag_DataPad,   &
                ChangeSizeAll_DataPad
      interface ChangeSize_DataPad
           module procedure ChangeSizeByInd_DataPad
           module procedure ChangeSizeByTag_DataPad
           module procedure ChangeSizeAll_DataPad
      end interface ChangeSize_DataPad
      !--------------------------------------------
      private:: CheckI1EQ_ByTag_DataPad, &
                CheckI1GT_ByTag_DataPad, &
                CheckI1GE_ByTag_DataPad, &
                CheckI1LT_ByTag_DataPad, &
                CheckI1LE_ByTag_DataPad, &

                CheckI2EQ_ByTag_DataPad, &
                CheckI2GT_ByTag_DataPad, &
                CheckI2GE_ByTag_DataPad, &
                CheckI2LT_ByTag_DataPad, &
                CheckI2LE_ByTag_DataPad, &

                CheckD1EQ_ByTag_DataPad, &
                CheckD1GT_ByTag_DataPad, &
                CheckD1GE_ByTag_DataPad, &
                CheckD1LT_ByTag_DataPad, &
                CheckD1LE_ByTag_DataPad, &

                CheckD2EQ_ByTag_DataPad, &
                CheckD2GT_ByTag_DataPad, &
                CheckD2GE_ByTag_DataPad, &
                CheckD2LT_ByTag_DataPad, &
                CheckD2LE_ByTag_DataPad

      interface CheckEQ_DataPad
          module procedure CheckI1EQ_ByTag_DataPad
          module procedure CheckD1EQ_ByTag_DataPad
          module procedure CheckI2EQ_ByTag_DataPad
          module procedure CheckD2EQ_ByTag_DataPad
      end interface CheckEQ_DataPad
      interface CheckGT_DataPad
          module procedure CheckI1GT_ByTag_DataPad
          module procedure CheckD1GT_ByTag_DataPad
          module procedure CheckI2GT_ByTag_DataPad
          module procedure CheckD2GT_ByTag_DataPad
      end interface CheckGT_DataPad
      interface CheckGE_DataPad
          module procedure CheckI1GE_ByTag_DataPad
          module procedure CheckD1GE_ByTag_DataPad
          module procedure CheckI2GE_ByTag_DataPad
          module procedure CheckD2GE_ByTag_DataPad
      end interface CheckGE_DataPad

      interface CheckLT_DataPad
          module procedure CheckI1LT_ByTag_DataPad
          module procedure CheckD1LT_ByTag_DataPad
          module procedure CheckI2LT_ByTag_DataPad
          module procedure CheckD2LT_ByTag_DataPad
      end interface
      interface CheckLE_DataPad
          module procedure CheckI1LE_ByTag_DataPad
          module procedure CheckD1LE_ByTag_DataPad
          module procedure CheckI2LE_ByTag_DataPad
          module procedure CheckD2LE_ByTag_DataPad
      end interface CheckLE_DataPad
      !--------------------------------------------
      private:: CopyFrom0_DataPad, &
                CopyFrom1_DataPad, &
                CopyFromA_DataPad, &
                CopyFromB_DataPad
      interface assignment (=)
          module procedure CopyFromA_DataPad
      end interface
      interface CopyFrom_DataPad
          module procedure CopyFromA_DataPad
          module procedure CopyFromB_DataPad
      end interface CopyFrom_DataPad
      !--------------------------------------------
      private:: DelElement0_DataPad,     &
                DelElementByInd_DataPad, &
                DelElementByTag_DataPad, &
                DelElementAll_DataPad
      interface DelElement_DataPad
           module procedure DelElementByInd_DataPad
           module procedure DelElementByTag_DataPad
           module procedure DelElementAll_DataPad
      end interface DelElement_DataPad

      private:: &
                DelElementsByInd_DataPad, &
                DelElementsByTag_DataPad, &
                DelElementsAll_DataPad
      interface DelElements_DataPad
           module procedure DelElementsByInd_DataPad
           module procedure DelElementsByTag_DataPad
           module procedure DelElementsAll_DataPad
      end interface DelElements_DataPad
      !--------------------------------------------
      private:: GetDataByTag_DataPad,     &
                GetDataPByTag_DataPad,    &
                GetDataIByTag_DataPad,    &
                GetDataI2ByTag_DataPad,   &
                GetDataDByTag_DataPad,    &

                GetDataByInd_DataPad,     &
                GetDataPByInd_DataPad,    &
                GetDataIByInd_DataPad,    &
                GetDataI2ByInd_DataPad,   &
                GetDataDByInd_DataPad,    &
                GetDataD2ByInd_DataPad
      interface GetData_DataPad
          module procedure GetDataByInd_DataPad
          !module procedure GetDataPByInd_DataPad
          module procedure GetDataIByInd_DataPad
          module procedure GetDataI2ByInd_DataPad
          module procedure GetDataDByInd_DataPad
          module procedure GetDataD2ByInd_DataPad
          module procedure GetDataByTag_DataPad
          !module procedure GetDataPByTag_DataPad
          module procedure GetDataIByTag_DataPad
          module procedure GetDataI2ByTag_DataPad
          module procedure GetDataDByTag_DataPad
          module procedure GetDataD2ByTag_DataPad
      end interface GetData_DataPad
      !--------------------------------------------
      private:: Dim_DataPad,              &
                GetSizeByInd_DataPad,     &
                GetSizeByTag_DataPad
      interface GetSize_DataPad
           module procedure GetSizeByInd_DataPad
           module procedure GetSizeByTag_DataPad
      end interface
      !--------------------------------------------
      public::HasTag_DataPad
      !--------------------------------------------
      public::Ind_DataPad
      !--------------------------------------------
      public::IsInteger_DataPad
      !--------------------------------------------
      public::IsDouble_DataPad
      !--------------------------------------------
      private:: Merge0_DataPad,     &
                MergeByTag_DataPad, &
                MergeAll_DataPad
      interface Merge_DataPad
           module procedure MergeByTag_DataPad
           module procedure MergeAll_DataPad
      end interface Merge_DataPad
      !--------------------------------------------
      public::  NumberofData_DataPad
      !--------------------------------------------
      private:: New_DataPad_NULL, &
                New_DataPad_I,    &
                New_DataPad_D,    &
                New_DataPad_I2,   &
                New_DataPad_D2,   &
                New_DataPad_Type, &
                New_DataPad_Type2
      interface New_DataPad
          module procedure New_DataPad_NULL
          module procedure New_DataPad_I
          module procedure New_DataPad_D
          module procedure New_DataPad_I2
          module procedure New_DataPad_D2
          module procedure New_DataPad_Type
          module procedure New_DataPad_Type2
      end interface New_DataPad
      !--------------------------------------------
      public::Tag_DataPad
      !--------------------------------------------
      private:: ReadByRow_DataPad
      interface Read_DataPad
          module procedure ReadByRow_DataPad
      end interface Read_DataPad
      !--------------------------------------------
      private:: ReleaseAll_DataPad, &
                ReleaseByTag_DataPad
      interface Release_DataPad
          module procedure ReleaseAll_DataPad
          module procedure ReleaseByTag_DataPad
      end interface Release_DataPad
      !--------------------------------------------
      private:: ResetByInd_DataPad,   &
                ResetByTag_DataPad
      interface Reset_DataPad
          module procedure ResetByTag_DataPad
          module procedure ResetByInd_DataPad
      end interface Reset_DataPad
      !--------------------------------------------
      private:: SetValByInd_DataPad_1I_A,   &
                SetValByTag_DataPad_1I_A,   &
                SetValByInd_DataPad_1I_B,   &
                SetValByTag_DataPad_1I_B,   &
                SetValByInd_DataPad_1D_A,   &
                SetValByTag_DataPad_1D_A,   &
                SetValByInd_DataPad_1D_B,   &
                SetValByTag_DataPad_1D_B,   &
                SetValByTag_DataPad_2D_A,   &
                SetValByInd_DataPad_2D_A
      interface SetVal_DataPad
          module procedure SetValByInd_DataPad_1I_A
          module procedure SetValByTag_DataPad_1I_A
          module procedure SetValByInd_DataPad_1I_B
          module procedure SetValByTag_DataPad_1I_B
          module procedure SetValByInd_DataPad_1D_A
          module procedure SetValByTag_DataPad_1D_A
          module procedure SetValByInd_DataPad_1D_B
          module procedure SetValByTag_DataPad_1D_B
          module procedure SetValByInd_DataPad_2D_A
          module procedure SetValByTag_DataPad_2D_A
      end interface SetVal_DataPad

      !--------------------------------------------
      public:: Restore_DataPad
      !--------------------------------------------
      private:: WriteByRowFmt_DataPad,  &
                WriteByRowUnFmt_DataPad
      interface Write_DataPad
          module procedure WriteByRowFmt_DataPad
          module procedure WriteByRowUnFmt_DataPad
      end interface Write_DataPad
      !--------------------------------------------

  contains
  !****************************************************************************
  !****************************************************************************
  recursive subroutine Tag_DataPad(Ind, DatList, Tag)
  !***  PURPOSE:   to get the Tag of data by index
  !
  !     INPUT:    Index,      the index for the data
  !               DatList,    the data list
  !     OUTPUT    Tag,        the Tag of the data
  !
  implicit none
     !--- dummy varioables
     integer,       intent(in) ::Ind
     type(DataPad), pointer    ::DatList
     character*(*)             ::Tag
     !----
     integer::ITH

          if(.not.associated(DatList)) then
             Tag = ""
             return
          else
             ITH = IND - 1
             if(ITH .eq. 0) then
                call GetFirstWord(DatList%Tag, Tag)
                return
             else
                 call Tag_DataPad(ITH, DatList%Next, Tag)
             end if
          end if
          return
  end subroutine Tag_DataPad
  !****************************************************************************

  !****************************************************************************
  subroutine Ind_DataPad(Tag, DatList, Ind)
  !***  PURPOSE:   to get the Index of data with tag
  !
  !     INPUT:    Tag,        the Tag of the data
  !               DatList,    the data list
  !     OUTPUT    Index,      the index for the data
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in) ::Tag
     type(DataPad), pointer    ::DatList
     integer                   ::Ind
     !----
     character(len=mp_TAGLEN) ::tTag
     integer::I

          Ind = 0
          do I = 1, NumberofData_DataPad(DatList)
             call Tag_DataPad(I, DatList, tTag)
             if(SameFirstWord(Tag, tTag)) then
                Ind = I
                return
             end if
          end do
          return
  end subroutine Ind_DataPad
  !****************************************************************************


  !****************************************************************************
  recursive logical function HasTag_DataPad(Tag, DatList) result(yes)
  !***  PURPOSE:   to find if the Tag of data exist
  !
  !     INPUT:    DatList,    the data list
  !               Tag,        the Tag of the data
  !     OUTPUT
  !
  implicit none
     !--- dummy varioables
     type(DataPad), pointer  ::DatList
     character*(*)           ::Tag
     !----

          if(.not.associated(DatList)) then
             yes = .false.
             return
          else
             if(SameFirstWord(Tag, DatList%Tag) ) then
                yes = .true.
                return
             else
                yes = HasTag_DataPad(Tag, DatList%Next)
             end if
          end if
          return
  end function HasTag_DataPad
  !****************************************************************************

  !****************************************************************************
  subroutine New_DataPad_NULL(Tag, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     character*(*)                 ::Tag
     type(DataPad),         pointer::DatList

     !----
              call New_DataPad_Type(Tag, "I", 0, DatList)
       return
  end subroutine New_DataPad_NULL
  !****************************************************************************

  !****************************************************************************
  subroutine New_DataPad_Type(Tag, Ch, Dim, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     character*(*)                 ::Tag
     character*(1)                 ::Ch
     integer, intent(in)           ::Dim
     type(DataPad),         pointer::DatList

     !----
     integer,      dimension(:), pointer::PadI
     real(KINDDF), dimension(:), pointer::PadD

           if(Ch(1:1) .eq. 'I' .or. Ch(1:1) .eq. 'i') then
              call New_DataPad_I(Tag, Dim, PadI, DatList)
           else if(Ch(1:1) .eq. 'D' .or. Ch(1:1) .eq. 'd') then
              call New_DataPad_D(Tag, Dim, PadD, DatList)
           end if
       return
  end subroutine New_DataPad_Type
  !****************************************************************************
  !****************************************************************************
  subroutine New_DataPad_Type2(Tag, Ch, Dim, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     character*(*)                 ::Tag
     character*(1)                 ::Ch
     integer, intent(in)           ::Dim(2)
     type(DataPad),         pointer::DatList

     !----
     integer,      dimension(:,:), pointer::PadI
     real(KINDDF), dimension(:,:), pointer::PadD

           if(Ch(1:1) .eq. 'I' .or. Ch(1:1) .eq. 'i') then
              if(Dim(2).le.1) then
                 call New_DataPad_I(Tag, Dim(1), PadI, DatList)
              else   
                 call New_DataPad_I2(Tag, Dim, PadI, DatList)
              end if   
           else if(Ch(1:1) .eq. 'D' .or. Ch(1:1) .eq. 'd') then
              if(Dim(2).le.1) then
                 call New_DataPad_D(Tag,  Dim(1), PadD, DatList)
              else    
                 call New_DataPad_D2(Tag, Dim, PadD, DatList)
              end if   
           end if
       return
  end subroutine New_DataPad_Type2
  !****************************************************************************

  !****************************************************************************
  recursive subroutine New_DataPad_I(Tag, Dim, Pad, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     character*(*)                 ::Tag
     integer, intent(in)           ::Dim
     integer, dimension(:), pointer::Pad
     type(DataPad),         pointer::DatList

     !----

          if(.not.associated(DatList)) then
             allocate(DatList)
             DatList%Tag     = Tag(1:min(len(Tag), len(DatList%Tag )))
             if(Dim .gt. 0) then
               allocate(Pad(Dim))
               DatList%IDATA1  => Pad
             else
               DatList%IDATA1 => null()
             end if
             DatList%DDATA1  => null()
             DatList%IDATA2  => null()
             DatList%DDATA2  => null()
             DatList%Next    => null()
             DatList%Prev    => null()
          else
             if(SameFirstWord(Tag, DatList%Tag) ) then
                if(associated(DatList%IDATA1))    then
                   if(Dim .ne. size(DatList%IDATA1)) then
                      deallocate(DatList%IDATA1)
                      DatList%IDATA1 => null()
                   end if
                end if

                if(.not. associated(DatList%IDATA1)) then
                   if(Dim .gt. 0) then
                      allocate(Pad(Dim))
                      DatList%IDATA1  => Pad
                   else
                      DatList%IDATA1  => null()
                   end if
                   DatList%DDATA1  => null()
                   DatList%IDATA2  => null()
                   DatList%DDATA2  => null()
                else
                   Pad => DatList%IDATA1
                end if
             else
                call New_DataPad_I(Tag, Dim, Pad, DatList%Next)
                DatList%Next%Prev=>DatList
             end if
          end if
          return
  end subroutine New_DataPad_I
 !****************************************************************************

  !****************************************************************************
  recursive subroutine New_DataPad_I2(Tag, Dim, Pad, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     character*(*)                            ::Tag
     integer, dimension(2), intent(in)        ::Dim
     integer, dimension(:,:),          pointer::Pad
     type(DataPad),                    pointer::DatList

     !----
          if(Dim(2) .le. 1) then
            call New_DataPad_I(Tag, Dim(1), Pad, DatList)
            return  
          end if  

          if(.not.associated(DatList)) then
             allocate(DatList)
             DatList%Tag     = Tag(1:min(len(Tag), len(DatList%Tag )))
             DatList%IDATA1  => null()
             if(all(Dim .gt. 0)) then
                allocate(Pad(Dim(1), Dim(2)))
                DatList%IDATA2  => Pad
             else
                DatList%IDATA2  => null()
             end if
             DatList%DDATA1  => null()
             DatList%DDATA2  => null()
             DatList%Next    => null()
             DatList%Prev    => null()
          else
             if(SameFirstWord(Tag, DatList%Tag) ) then
                if(associated(DatList%IDATA2)) then
                   if(Dim(1) .ne. size(DatList%IDATA2,dim=1) .or. &
                      Dim(2) .ne. size(DatList%IDATA2,dim=2)) then
                      deallocate(DatList%IDATA2)
                      DatList%IDATA2 => null()
                   end if
                end if
                if(.not. associated(DatList%IDATA2)) then
                   DatList%IDATA1  => null()
                   if(all(Dim .gt. 0) ) then
                      allocate(Pad(Dim(1),Dim(2)))
                      DatList%IDATA2  => Pad
                   else
                      DatList%IDATA2  => null()
                   end if
                   DatList%DDATA1  => null()
                   DatList%DDATA2  => null()
                else
                   Pad => DatList%IDATA2
                end if
             else
                call New_DataPad_I2(Tag, Dim, Pad, DatList%Next)
                DatList%Next%Prev=>DatList
             end if
          end if
          return
  end subroutine New_DataPad_I2
 !****************************************************************************

  !****************************************************************************
  recursive subroutine New_DataPad_D(Tag, Dim, Pad, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     !--- dummy varioables
     character*(*)                     ::Tag
     integer,     intent(in)           ::Dim
     real(KINDDF),dimension(:), pointer::Pad
     type(DataPad),             pointer::DatList

     !----
          if(.not.associated(DatList)) then
             allocate(DatList)
             DatList%Tag     = Tag(1:min(len(Tag), len(DatList%Tag )))
             DatList%IDATA1  => null()
             DatList%IDATA2  => null()
             if(Dim .gt. 0) then
                allocate(Pad(Dim))
                DatList%DDATA1  => Pad
             else
                DatList%DDATA1  => null()
             end if
             DatList%DDATA2  => null()
             DatList%Next    => null()
             DatList%Prev    => null()
          else
             if(SameFirstWord(Tag, DatList%Tag) ) then
                if(associated(DatList%DDATA1)) then
                   if(Dim .ne. size(DatList%DDATA1)) then
                      deallocate(DatList%DDATA1)
                      DatList%DDATA1 => null()
                   end if
                end if

                if(.not.associated(DatList%DDATA1)) then
                   DatList%IDATA1  => null()
                   DatList%IDATA2  => null()
                   if(Dim .gt. 0) then
                      allocate(Pad(Dim))
                      DatList%DDATA1  => Pad
                   else
                      DatList%DDATA1  => null()
                   end if
                   DatList%DDATA2  => null()
                else
                   Pad => DatList%DDATA1
                end if
             else
                call New_DataPad_D(Tag, Dim, Pad, DatList%Next)
                DatList%Next%Prev=>DatList
             end if
          end if
          return
  end subroutine New_DataPad_D
 !****************************************************************************

 !****************************************************************************
  recursive subroutine New_DataPad_D2(Tag, Dim, Pad, DatList)
  !***  PURPOSE:   to create a data pad and link to a list
  !
  !     INPUT:    Tag,        the Tag for the data
  !               Dim,        dimension of Pad
  !
  !     OUTPUT    DatList
  !               Pad,
  !
  implicit none
     !--- dummy varioables
     !--- dummy varioables
     character*(*)                          ::Tag
     integer,     dimension(2),   intent(in)::Dim
     real(KINDDF),dimension(:,:), pointer   ::Pad
     type(DataPad),               pointer   ::DatList

     !----

           if(Dim(2) .le. 1) then
              call New_DataPad_D(Tag, Dim(1), Pad, DatList)
              return  
           end if  

          if(.not.associated(DatList)) then
             allocate(DatList)
             DatList%Tag     = Tag(1:min(len(Tag), len(DatList%Tag )))
             DatList%IDATA1  => null()
             DatList%IDATA2  => null()
             DatList%DDATA1  => null()
             if(all(Dim .gt. 0) ) then
                allocate(Pad(Dim(1),Dim(2)))
                DatList%DDATA2  => Pad
             else
                DatList%DDATA2  => null()
             end if
             DatList%Next    => null()
             DatList%Prev    => null()
          else
             if(SameFirstWord(Tag, DatList%Tag) ) then
                if(associated(DatList%DDATA2)) then
                   if(Dim(1) .ne. size(DatList%DDATA2,dim=1) .or. &
                      Dim(2) .ne. size(DatList%DDATA2,dim=2)) then
                      deallocate(DatList%DDATA2)
                      DatList%DDATA2 =>null()
                   endif
                end if

                if(.not. associated(DatList%DDATA2)) then
                   allocate(Pad(Dim(1),Dim(2)))
                   DatList%IDATA1  => null()
                   DatList%DDATA1  => null()
                   DatList%IDATA2  => null()
                   if(all(Dim .gt. 0) ) then
                      allocate(Pad(Dim(1),Dim(2)))
                      DatList%DDATA2  => Pad
                   else
                      DatList%DDATA2  => null()
                   end if
                else
                   Pad => DatList%DDATA2
                end if
             else
                call New_DataPad_D2(Tag, Dim, Pad, DatList%Next)
                DatList%Next%Prev=>DatList
             end if
          end if
          return
  end subroutine New_DataPad_D2
  !****************************************************************************

  !****************************************************************************
  recursive subroutine ReleaseByTag_DataPad(Tag, DatList)
  !***  PURPOSE:   to release a data pad list.
  !                Note: the allocated memeory do not deallocated, we set only
  !                      the pointers to null.
  !
  !     INPUT:    Tag,        the Tag for the data
  !
  !     OUTPUT    DatList
  !
  implicit none
     !--- dummy varioables
     character*(*)                                 ::Tag
     type(DataPad), pointer                        ::DatList
     !----

          if(.not.associated(DatList)) then
             return
          else
             if(SameFirstWord(Tag, DatList%Tag) ) then
                if(associated(DatList%Prev) ) then
                   DatList%Prev%Next =>DatList%Next
                end if

                if(associated(DatList%Next)) then
                   DatList%Next%Prev=>DatList%Prev
                end if
                !$$ deallocate memory
                if(associated(DatList%IDATA1)) then
                   deallocate(DatList%IDATA1)
                end if
                if(associated(DatList%IDATA2)) then
                   deallocate(DatList%IDATA2)
                end if
                if(associated(DatList%DDATA1)) then
                   deallocate(DatList%DDATA1)
                end if
                if(associated(DatList%DDATA2)) then
                   deallocate(DatList%DDATA2)
                end if
                deallocate(DatList)

                DatList=>null()
                return
             else
                 call ReleaseByTag_DataPad(Tag, DatList%Next)
             end if
          end if
          return
  end subroutine ReleaseByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  recursive subroutine ReleaseAll_DataPad(DatList)
  !***  PURPOSE:   to release a data pad list.
  !                Note: the allocated DatPad memeory pointerd by the DataList is deallocated
  !
  !     INPUT:    DatList
  !
  !     OUTPUT    DatList
  !
  implicit none
     !--- dummy varioables
     type(DataPad), pointer::DatList
     !----

          if(.not.associated(DatList)) then
             return
          else
             call Release_DataPad(DatList%Next)
             !$$ deallocate memory
             if(associated(DatList%IDATA1)) then
                 deallocate(DatList%IDATA1)
             end if
             if(associated(DatList%IDATA2)) then
                deallocate(DatList%IDATA2)
             end if
             if(associated(DatList%DDATA1)) then
                 deallocate(DatList%DDATA1)
             end if
             if(associated(DatList%DDATA2)) then
                deallocate(DatList%DDATA2)
             end if
             deallocate(DatList)
             DatList=>null()
          end if
          return
  end subroutine ReleaseAll_DataPad
 !****************************************************************************

  !****************************************************************************
   subroutine CopyFrom0_DataPad(To, From, At, For)
  !***  PURPOSE:   to copy part of a data pad of From to another data pad To
  !                Note: new memories are allocated in datapad To
  !
  !     INPUT:    From,
  !               At,   the starting index of array
  !               For,  the size to be copy
  !
  !     OUTPUT    To,   the data pad with the same tag of From
  !                     but with new memories allocated
  !
  implicit none
     !--- dummy varioables
     type(DataPad), intent(in), pointer::From
     type(DataPad),             pointer::To
     integer, intent(in)               ::At
     integer, optional                 ::For

     !--- local varibales
     integer,     dimension(:),   pointer::IDATA1
     real(KINDDF),dimension(:),   pointer::DDATA1
     integer,     dimension(:,:), pointer::IDATA2
     real(KINDDF),dimension(:,:), pointer::DDATA2
     integer::IP0, IP1, DIM(2)


          if(.not.associated(From)) then
              return
          end if

             call Dim_DataPad(DIM, From)
             IP0 = At
             if(present(For)) then
                IP1 = IP0 + For -1
                if(IP1 .gt. DIM(1)) IP1 = DIM(1)
             else
                IP1 = DIM(1)
             end if

             if( associated(From%IDATA1)) then
                 call New_DataPad(From%Tag, IP1-IP0+1, IDATA1, To)
                 IDATA1 = From%IDATA1(IP0:IP1)

             else if( associated(From%DDATA1) ) then
                 call New_DataPad(From%Tag, IP1-IP0+1,  DDATA1, To)
                 DDATA1 = From%DDATA1(IP0:IP1)

             else if( associated(From%IDATA2) ) then
                 call New_DataPad(From%Tag, (/IP1-IP0+1,size(From%IDATA2,dim=2)/),  IDATA2, To)
                 IDATA2 = From%IDATA2(IP0:IP1,:)

             else if( associated(From%DDATA2) ) then
                 call New_DataPad(From%Tag, (/IP1-IP0+1,size(From%DDATA2,dim=2)/),  DDATA2, To)
                 DDATA2 = From%DDATA2(IP0:IP1,:)

             else
                 call New_DataPad_NULL(From%Tag, To)
             end if

          return
  end subroutine CopyFrom0_DataPad
 !****************************************************************************

   !****************************************************************************
   subroutine CopyFrom1_DataPad(To, From)
  !***  PURPOSE:   to copy a data pad From to another data pad To
  !                Note: new memories are allocated in data pad To
  !
  !     INPUT:    From,
  !
  !     OUTPUT    To,   the data pad with the same tag of From
  !                     but with new memories allocated
  !
  implicit none
     !--- dummy varioables
     type(DataPad), intent(in), pointer::From
     type(DataPad),             pointer::To

          if(.not.associated(From)) then
              return
          end if

          call CopyFrom0_DataPad(To, From, 1)

          return
  end subroutine CopyFrom1_DataPad
 !****************************************************************************

 !****************************************************************************
   recursive subroutine CopyFromA_DataPad(To, From)
  !***  PURPOSE:   to copy a data pad From to another data pad To
  !                Note: new memories are allocated in data pad To
  !
  !     INPUT:    From,
  !
  !     OUTPUT    To,   the data pad with the same tag of From
  !                     but with new memories allocated
  !
  implicit none
     !--- dummy varioables
     type(DataPad), intent(in), pointer::From
     type(DataPad),             pointer::To

     !----
          call Release_DataPad(To)
          if(.not.associated(From)) then
              return
          end if

             call CopyFrom1_DataPad(To, From)
             call CopyFromA_DataPad(To%Next, From%Next)
          return
  end subroutine CopyFromA_DataPad
 !****************************************************************************

  !****************************************************************************
   recursive subroutine CopyFromB_DataPad(To, From, At, For)
  !***  PURPOSE:   to copy a subsection of a data pad From to another data pad To
  !                Note: new memories are allocated in data pad To
  !
  !     INPUT:    From,
  !
  !     OUTPUT    To,   the data pad with the same tag of From
  !                     but with new memories allocated
  !
  implicit none
     !--- dummy varioables
     type(DataPad), intent(in), pointer::From
     type(DataPad),             pointer::To
     integer, intent(in)               ::At
     integer, intent(in), optional     ::For

     !----
          call Release_DataPad(To)
          if(.not.associated(From)) then
              return
          end if

             if(present(For)) then
                call CopyFrom0_DataPad(To, From, At, For)
                call CopyFromB_DataPad(To%Next, From%Next, For)
             else
                call CopyFrom0_DataPad(To, From, At)
                call CopyFromB_DataPad(To%Next, From%Next, At)
             end if
          return
  end subroutine CopyFromB_DataPad
 !****************************************************************************

 !****************************************************************************
  recursive integer function NumberofData_DataPad(DatList) result(Num)
  !***  PURPOSE:   to get number of padded data on the list.
  !
  !     INPUT:    DatList,        the Tag for the data
  !
  !     OUTPUT    Num
  !
  implicit none
     !--- dummy varioables
     type(DataPad), pointer ::DatList
     !----
     integer::N

          Num = 0
          if(.not.associated(DatList)) then
             return
          else
             Num = NumberofData_DataPad(DatList%Next) + 1
          end if
          return
  end function NumberofData_DataPad
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetDataPByTag_DataPad(Tag, DatList, pDat)
  !***  PURPOSE:   to the the DataPad by tag
  !
  !     INPUT:    Tag,        the Tag for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the DataPad pointer
  !
  implicit none
     !--- dummy varioables
     character*(*)           ::Tag
     type(DataPad), pointer  ::DatList
     type(DataPad), pointer  ::pDat
     !----

          if(.not.associated(DatList)) then
             pDat => null()
             return
          else
             if(SameFirstWord(Tag, DatList%Tag)) then
                pDat => DatList
             else
                 call GetDataPByTag_DataPad(Tag, DatList%Next, pDat)
             end if
          end if
          return
  end subroutine GetDataPByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine GetDataByTag_DataPad(Tag, DatList, Dat)
  !***  PURPOSE:   to the the DataPad by tag
  !
  !     INPUT:    Tag,        the Tag for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     character*(*)                 ::Tag
     type(DataPad), pointer        ::DatList
     type(DataPad)                 ::Dat
     !----
     type(DataPad), pointer  ::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat) ) then
             Dat%Tag     = pDat%Tag
             Dat%IDATA1  => pDat%IDATA1
             Dat%IDATA2  => pDat%IDATA2
             Dat%DDATA1  => pDat%DDATA1
             Dat%DDATA2  => pDat%DDATA2
          else
             Dat%Tag     = ""
             Dat%IDATA1  => null()
             Dat%IDATA2  => null()
             Dat%DDATA1  => null()
             Dat%DDATA2  => null()
          end if
          return
  end subroutine GetDataByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine GetDataIByTag_DataPad(Tag, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by tag
  !
  !     INPUT:    Tag,        the Tag for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     character*(*)                 ::Tag
     type(DataPad), pointer        ::DatList
     integer, dimension(:), pointer::pArray
     !----
     type(DataPad),pointer ::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat) ) then
             pArray => pDat%IDATA1
          else
             pArray => null()
          end if
          return
  end subroutine GetDataIByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine GetDataI2ByTag_DataPad(Tag, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by tag
  !
  !     INPUT:    Tag,        the Tag for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     character*(*)                   ::Tag
     type(DataPad), pointer          ::DatList
     integer, dimension(:,:), pointer::pArray
     !----
     type(DataPad), pointer::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat) ) then
             pArray => pDat%IDATA2
          else
             pArray => null()
          end if
          return
  end subroutine GetDataI2ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine GetDataDByTag_DataPad(Tag, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by tag
  !
  !     INPUT:    Tag,        the Tag for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     character*(*)                      ::Tag
     type(DataPad), pointer             ::DatList
     real(KINDDF), dimension(:), pointer::pArray
     !----
     type(DataPad), pointer::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat) ) then
             pArray => pDat%DDATA1
          else
             pArray => null()
          end if

          return
  end subroutine GetDataDByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine GetDataD2ByTag_DataPad(Tag, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by tag
  !
  !     INPUT:    Tag,        the Tag for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     character*(*)                        ::Tag
     type(DataPad), pointer               ::DatList
     real(KINDDF), dimension(:,:), pointer::pArray
     !----
     type(DataPad), pointer::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat) ) then
             pArray => pDat%DDATA2
          else
             pArray => null()
          end if

          return
  end subroutine GetDataD2ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  recursive subroutine GetDataPByInd_DataPad(Ind, DatList, pDat)
  !***  PURPOSE:   to the the DataPad by index
  !
  !     INPUT:    Index,      the index for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the DataPad pointer
  !
  implicit none
     !--- dummy varioables
     integer                 ::Ind
     type(DataPad), pointer  ::DatList
     type(DataPad), pointer  ::pDat
     !----
     integer::ITH

          pDat => null()
          if(.not.associated(DatList)) then
             return
          else
             ITH = IND - 1
             if(ITH .eq. 0) then
                pDat => DatList
                return
             else
                 call GetDataPByInd_DataPad(ITH, DatList%Next, pDat)
             end if
          end if
          return
  end subroutine GetDataPByInd_DataPad
  !****************************************************************************

 !****************************************************************************
  subroutine GetDataByInd_DataPad(Ind, DatList, Dat)
  !***  PURPOSE:   to the the DataPad by index
  !
  !     INPUT:    Index,      the index for the data
  !               DatList,    the data list
  !     OUTPUT    Dat,        the DataPad copy
  !
  implicit none
     !--- dummy varioables
     integer                 ::Ind
     type(DataPad), pointer  ::DatList
     type(DataPad)           ::Dat
     !----
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) then
             Dat%Tag    =  pDat%Tag
             Dat%IDATA1 => pDat%IDATA1
             Dat%DDATA1 => pDat%DDATA1
             Dat%IDATA2 => pDat%IDATA2
             Dat%DDATA2 => pDat%DDATA2
          else
             Dat%Tag = ""
             Dat%IDATA1 => null()
             Dat%DDATA1 => null()
             Dat%IDATA2 => null()
             Dat%DDATA2 => null()
          end if
          return
  end subroutine GetDataByInd_DataPad
  !****************************************************************************

  !****************************************************************************
  subroutine GetDataIByInd_DataPad(Ind, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by index
  !
  !     INPUT:    Ind,        the index for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     integer                       ::Ind
     type(DataPad), pointer        ::DatList
     integer, dimension(:), pointer::pArray
     !----
     type(DataPad), pointer  ::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) then
             pArray => pDat%IDATA1
          else
             pArray => null()
          end if

          return
  end subroutine GetDataIByInd_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine GetDataI2ByInd_DataPad(Ind, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by index
  !
  !     INPUT:    Ind,        the index for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     integer                         ::Ind
     type(DataPad), pointer          ::DatList
     integer, dimension(:,:), pointer::pArray
     !----
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) then
             pArray => pDat%IDATA2
          else
             pArray => null()
          end if

          return
  end subroutine GetDataI2ByInd_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine GetDataDByInd_DataPad(Ind, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by index
  !
  !     INPUT:    Ind,        the index for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     integer                             ::Ind
     type(DataPad), pointer              ::DatList
     real(KINDDF),  dimension(:), pointer::pArray
     !----
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) then
             pArray => pDat%DDATA1
          else
             pArray => null()
          end if

          return
  end subroutine GetDataDByInd_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine GetDataD2ByInd_DataPad(Ind, DatList, pArray)
  !***  PURPOSE:   to the the DataPad by index
  !
  !     INPUT:    Ind,        the index for the data
  !               DatList,    the data list
  !     OUTPUT    pDat,       the array pointer
  !
  implicit none
     !--- dummy varioables
     integer                              ::Ind
     type(DataPad), pointer               ::DatList
     real(KINDDF), dimension(:,:), pointer::pArray
     !----
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) then
             pArray => pDat%DDATA2
          else
             pArray => null()
          end if

          return
  end subroutine GetDataD2ByInd_DataPad
 !****************************************************************************

  !****************************************************************************
  logical function IsInteger_DataPad(PDat) result(YES)
  !***  PURPOSE:   to require if the data type is integer
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     type(DataPad)::PDat

          if( associated(PDat%IDATA1) .or. associated(PDat%IDATA2) ) then
             YES  = .true.
          else
             YES  = .false.
          end if
          return
  end function IsInteger_DataPad
 !****************************************************************************

  !****************************************************************************
  logical function IsDouble_DataPad(PDat) result(YES)
  !***  PURPOSE:   to require if the data type is double
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     type(DataPad)::PDat

          if( associated(PDat%DDATA1) .or. associated(PDat%DDATA2) ) then
             YES  = .true.
          else
             YES  = .false.
          end if
          return
  end function IsDouble_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine Dim_DataPad(Dim, PDat)
  !***  PURPOSE:   to require dimension of the data
  !     INPUT:    DatList
  !     OUTPUT    Dim
  !
  implicit none
     !--- dummy varioables
     type(DataPad), pointer::PDat
     integer      ::Dim(2)

          if(associated(PDat%IDATA1))      then
             Dim(1)  = size(PDat%IDATA1, dim=1)
             Dim(2)  = 1
          else if(associated(PDat%IDATA2) ) then
             Dim(1)  = size(PDat%IDATA2, dim=1)
             Dim(2)  = size(PDat%IDATA2, dim=2)

          else if(associated(PDat%DDATA1) ) then
             Dim(1)  = size(PDat%DDATA1, dim=1)
             Dim(2)  = 1

          else if(associated(PDat%DDATA2) ) then
             Dim(1)  = size(PDat%DDATA2, dim=1)
             Dim(2)  = size(PDat%DDATA2, dim=2)
          end if
          return
  end subroutine Dim_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine GetSizeByTag_DataPad(Tag, DatList, Dim)
  !***  PURPOSE:   to require dimension of the data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     character*(*)          ::Tag
     integer                ::Dim(2)
     type(DataPad),  pointer::DatList

     !--- local
     type(DataPad), pointer::pDat


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat)) then
              call Dim_DataPad(Dim, pDat)
          else
              Dim = 0
          end if

          return
  end subroutine GetSizeByTag_DataPad
 !****************************************************************************

!****************************************************************************
  subroutine GetSizeByInd_DataPad(Ind, DatList, Dim)
  !***  PURPOSE:   to require dimension of the data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                ::Ind
     integer                ::Dim(2)
     type(DataPad),  pointer::DatList

     !--- local
     type(DataPad), pointer::pDat


          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) then
              call Dim_DataPad(Dim, pDat)
          else
              Dim = 0
          end if

          return
  end subroutine GetSizeByInd_DataPad
 !****************************************************************************


  !****************************************************************************
  subroutine ChangeSize_DataPad0(Dim, pDat)
  !***  PURPOSE:   to change the size of data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Dim
     type(DataPad), pointer    ::pDat
     !--- local
     integer,     dimension(:),   pointer::IDATA1 = null()
     real(KINDDF),dimension(:),   pointer::DDATA1 = null()
     integer,     dimension(:,:), pointer::IDATA2 = null()
     real(KINDDF),dimension(:,:), pointer::DDATA2 = null()


          if(associated(PDat%IDATA1))  then
             allocate(IDATA1(Dim))
             IDATA1(1:min(Dim,size(PDat%IDATA1))) = pDat%IDATA1(1:min(Dim,size(PDat%IDATA1)))
             deallocate(PDat%IDATA1)
             PDat%IDATA1 => IDATA1

          else if(associated(PDat%IDATA2) ) then
             allocate(IDATA2(Dim, size(PDat%IDATA2, dim=2)))
             IDATA2(1:min(Dim,size(PDat%IDATA2, dim=1)),:) = pDat%IDATA2(1:min(Dim,size(PDat%IDATA2, dim=1)),:)
             deallocate(PDat%IDATA2)
             PDat%IDATA2 => IDATA2

          else if(associated(PDat%DDATA1) ) then
             allocate(DDATA1(Dim))
             DDATA1(1:min(Dim,size(PDat%DDATA1))) = pDat%DDATA1(1:min(Dim,size(PDat%DDATA1)))
             deallocate(PDat%DDATA1)
             PDat%DDATA1 => DDATA1

          else if(associated(PDat%DDATA2) ) then
             allocate(DDATA2(Dim, size(PDat%DDATA2, dim=2)))
             DDATA2(1:min(Dim,size(PDat%DDATA2, dim=1)),:) = pDat%DDATA2(1:min(Dim,size(PDat%DDATA2, dim=1)),:)
             deallocate(PDat%DDATA2)
             PDat%DDATA2 => DDATA2

          end if
          return
  end subroutine ChangeSize_DataPad0
 !****************************************************************************

 !****************************************************************************
  subroutine ChangeSizeByTag_DataPad(Tag, Dim, DatList)
  !***  PURPOSE:   to ChangeSize the size of data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     character*(*)             ::Tag
     integer                   ::Dim
     type(DataPad), pointer    ::DatList
     !--- local
     type(DataPad), pointer::pDat


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat)) call ChangeSize_DataPad0(Dim, pDat)
          return
  end subroutine ChangeSizeByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine ChangeSizeByInd_DataPad(Ind, Dim, DatList)
  !***  PURPOSE:   to ChangeSize the size of data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Ind
     integer                   ::Dim
     type(DataPad), pointer    ::DatList
     !--- local
     type(DataPad), pointer::pDat


          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat))  call ChangeSize_DataPad0(Dim, pDat)
          return
  end subroutine ChangeSizeByInd_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine ChangeSizeAll_DataPad(Dim, DatList)
  !***  PURPOSE:   to ChangeSize the size of data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Dim
     type(DataPad), pointer    ::DatList
     !--- local
     integer::I

          do I=1, NumberofData_DataPad(DatList)
             call ChangeSizeByInd_DataPad(I, Dim, DatList)
          end do
          return
  end subroutine ChangeSizeAll_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine WriteByRowUnfmt_DataPad(pDat, Row, Str)
  !***  PURPOSE:  to write a data to a string
  !     INPUT:    pDat
  !               Row,  the index of the element
  !               fmt,  optional the fmt
  !     OUTPUT    Str
  !
  implicit none
     !--- dummy varioables
     type(DataPad)           ::pDat
     integer                 ::Row
     character*(*)           ::Str
     !--- local

          if(associated(pDat%IDATA1))      then
             write(Str,*) pDat%IDATA1(Row)

          else if(associated(pDat%IDATA2) ) then
             write(Str,*) pDat%IDATA2(Row,:)

          else if(associated(pDat%DDATA1) ) then
             write(Str,*) pDat%DDATA1(Row)

          else if(associated(pDat%DDATA2) ) then
             write(Str,*) pDat%DDATA2(Row,:)
          end if

        return
  end subroutine WriteByRowUnfmt_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine WriteByRowFmt_DataPad(pDat, Row, Str, Fmt)
  !***  PURPOSE:  to write a data to a string
  !     INPUT:    pDat
  !               Row,  the index of the element
  !               fmt,  the fmt
  !     OUTPUT    Str
  !
  implicit none
     !--- dummy varioables
     type(DataPad)           ::pDat
     integer                 ::Row
     character*(*)           ::Str
     character*(*)           ::Fmt
     !--- local
     character*32::SN, TFMT

          if(associated(pDat%IDATA1))      then
             TFMT = "("//Fmt(1:len_trim(Fmt))//")"
             write(Str,fmt=TFMT) pDat%IDATA1(Row)

          else if(associated(pDat%IDATA2) ) then
             write(SN,*) size(pDat%IDATA2, dim=2)
             SN = adjustl(SN)
             TFMT = "("//SN(1:len_trim(SN))//"("//Fmt(1:len_trim(Fmt))//"))"
             write(Str,fmt=TFMT) pDat%IDATA2(Row,:)

          else if(associated(pDat%DDATA1) ) then
             TFMT = "("//Fmt(1:len_trim(Fmt))//")"
             write(Str,fmt=TFMT) pDat%DDATA1(Row)

          else if(associated(pDat%DDATA2) ) then
             write(SN,*) size(pDat%DDATA2, dim=2)
             SN = adjustl(SN)
             TFMT = "("//SN(1:len_trim(SN))//"("//Fmt(1:len_trim(Fmt))//"))"
             write(Str,fmt=TFMT) pDat%DDATA2(Row,:)
          end if

        return
  end subroutine WriteByRowFmt_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine ReadByRow_DataPad(pDat, Row, Str)
  !***  PURPOSE:  to read a data to a string
  !     INPUT:    pDat
  !               Row,  the index of the element
  !               Str
  !     OUTPUT    pDat
  !
  implicit none
     !--- dummy varioables
     type(DataPad)             ::pDat
     integer                   ::Row
     character*(*),dimension(:)::Str
     !--- local
     integer::I

          if(associated(pDat%IDATA1))      then
             read(Str(1),*) pDat%IDATA1(Row)

          else if(associated(pDat%IDATA2) ) then
             do I=1,  size(pDat%IDATA2, dim=2)
                read(Str(I),*) pDat%IDATA2(Row, I)
             end do
          else if(associated(pDat%DDATA1) ) then
             read(Str(1),*) pDat%DDATA1(Row)

          else if(associated(pDat%DDATA2) ) then
             do I=1,  size(pDat%DDATA2, dim=2)
                read(Str(I),*)pDat%DDATA2(Row,I)
             end do
          end if

        return
  end subroutine ReadByRow_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine DelElement0_DataPad(AT, NE, pDat)
  !***  PURPOSE:   to delelet an element of data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::AT, NE
     type(DataPad), pointer    ::pDat
     !--- local
     integer,     dimension(:),   pointer::IDATA1 = null()
     real(KINDDF),dimension(:),   pointer::DDATA1 = null()
     integer,     dimension(:,:), pointer::IDATA2 = null()
     real(KINDDF),dimension(:,:), pointer::DDATA2 = null()


          if(associated(PDat%IDATA1))  then
             allocate(IDATA1(size(PDat%IDATA1)-NE))
             IDATA1(1:AT-1) = pDat%IDATA1(1:AT-1)
             IDATA1(AT:size(IDATA1)) = pDat%IDATA1(AT+NE:size(PDat%IDATA1))
             deallocate(PDat%IDATA1)
             PDat%IDATA1 => IDATA1

          else if(associated(PDat%IDATA2) ) then
             allocate(IDATA2(size(PDat%IDATA2, dim=1)-NE, size(PDat%IDATA2, dim=2)))
             IDATA2(1:AT-1,:) = pDat%IDATA2(1:AT-1,:)
             IDATA2(AT:size(IDATA2, dim=1),:) = pDat%IDATA2(AT+NE:size(PDat%IDATA2, dim=1),:)
             deallocate(PDat%IDATA2)
             PDat%IDATA2 => IDATA2

          else if(associated(PDat%DDATA1) ) then
             allocate(DDATA1(size(PDat%DDATA1)-NE))
             DDATA1(1:AT-1) = pDat%DDATA1(1:AT-1)
             DDATA1(AT:size(DDATA1)) = pDat%DDATA1(AT+NE:size(PDat%DDATA1))
             deallocate(PDat%DDATA1)
             PDat%DDATA1 => DDATA1

          else if(associated(PDat%DDATA2) ) then
             allocate(DDATA2(size(PDat%DDATA2, dim=1)-NE, size(PDat%DDATA2, dim=2)))
             DDATA2(1:AT-1,:) = pDat%DDATA2(1:AT-1,:)
             DDATA2(AT:size(DDATA2, dim=1),:) = pDat%DDATA2(AT+NE:size(PDat%DDATA2, dim=1),:)
             deallocate(PDat%DDATA2)
             PDat%DDATA2 => DDATA2

          end if
          return
  end subroutine DelElement0_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine DelElementByTag_DataPad(Tag, AT, DatList)
  !***  PURPOSE:   to delelet an element of data
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     character*(*)             ::Tag
     integer                   ::AT
     type(DataPad), pointer    ::DatList
     !--- local
     type(DataPad), pointer    ::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat))call DelElement0_DataPad(At,1,pDat)

          return
  end subroutine DelElementByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine DelElementsByTag_DataPad(Tag, AT, NE, DatList)
  !***  PURPOSE:   to delete an element for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     character*(*)             ::Tag
     integer                   ::AT, NE
     type(DataPad), pointer    ::DatList
     !--- local
     type(DataPad), pointer    ::pDat

          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(associated(pDat)) call DelElement0_DataPad(At,NE,pDat)


          return
  end subroutine DelElementsByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine DelElementByInd_DataPad(Ind, At, DatList)
  !***  PURPOSE:   to delete an element for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Ind
     integer                   ::At
     type(DataPad), pointer    ::DatList
     !--- local
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) call DelElement0_DataPad(At,1,pDat)

          return
  end subroutine DelElementByInd_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine DelElementsByInd_DataPad(Ind, At, Ne, DatList)
  !***  PURPOSE:   to delete elements for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Ind
     integer                   ::At, Ne
     type(DataPad), pointer    ::DatList
     !--- local
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(associated(pDat)) call DelElement0_DataPad(At,Ne,pDat)

          return
  end subroutine DelElementsByInd_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine DelElementAll_DataPad(At, DatList)
  !***  PURPOSE:   to delete an element for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Ind
     integer                   ::At
     type(DataPad), pointer    ::DatList
     !--- local
     integer::I

          do I=1, NumberofData_DataPad(DatList)
             call DelElementByInd_DataPad(I, At, DatList)
          end do
          return
  end subroutine DelElementAll_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine DelElementsAll_DataPad(At, Ne, DatList)
  !***  PURPOSE:   to delete elements for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     integer                   ::Ind
     integer                   ::At, Ne
     type(DataPad), pointer    ::DatList
     !--- local
     integer::I

          do I=1, NumberofData_DataPad(DatList)
             call DelElementsByInd_DataPad(I, At, Ne, DatList)
          end do
          return
  end subroutine DelElementsAll_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine Merge0_DataPad(pDatFrom, pDatTo)
  !***  PURPOSE:   to meger the data of From to To
  !     INPUT:    pDatFrom, pDatTo
  !     OUTPUT    pDatTo
  !
  implicit none
     !--- dummy varioables
     type(DataPad), pointer    ::pDatFrom, pDatTo
     !--- local
     integer,     dimension(:),   pointer::IDATA1 = null()
     real(KINDDF),dimension(:),   pointer::DDATA1 = null()
     integer,     dimension(:,:), pointer::IDATA2 = null()
     real(KINDDF),dimension(:,:), pointer::DDATA2 = null()
     integer::NT, NF

          if(associated(pDatFrom%IDATA1))  then
            NF = size(pDatFrom%IDATA1)
            NT = 0
             if(associated(pDatTo%IDATA1) ) then
                NT = size(pDatTo%IDATA1)
             end if
             allocate(IDATA1(NT+NF))

             if(NT .gt. 0 ) then
                IDATA1(1:NT) = pDatTo%IDATA1(1:NT)
                deallocate(pDatTo%IDATA1)
             end if
             IDATA1(NT+1:NT+NF) = pDatFrom%IDATA1(1:NF)
             PDatTo%IDATA1 => IDATA1

          else if(associated(pDatFrom%IDATA2) ) then
            NF = size(pDatFrom%IDATA2, dim=1)
            NT = 0
             if(associated(pDatTo%IDATA2) ) then
                NT = size(pDatTo%IDATA2, dim=1)
             end if
             allocate(IDATA2(NT+NF, size(pDatFrom%IDATA2, dim=2)))

             if(NT .gt. 0 ) then
                IDATA2(1:NT,:) = pDatTo%IDATA2(1:NT,:)
                deallocate(pDatTo%IDATA2)
             end if
             IDATA2(NT+1:NT+NF,:) = pDatFrom%IDATA2(1:NF,:)
             PDatTo%IDATA2 => IDATA2

          else if(associated(pDatFrom%DDATA1) ) then
            NF = size(pDatFrom%DDATA1)
            NT = 0
             if(associated(pDatTo%DDATA1) ) then
                NT = size(pDatTo%DDATA1)
             end if
             allocate(DDATA1(NT+NF))

             if(NT .gt. 0 ) then
                DDATA1(1:NT) = pDatTo%DDATA1(1:NT)
                deallocate(pDatTo%DDATA1)
             end if
             DDATA1(NT+1:NT+NF) = pDatFrom%DDATA1(1:NF)
             PDatTo%DDATA1 => DDATA1

          else if(associated(pDatFrom%DDATA2) ) then
            NF = size(pDatFrom%DDATA2, dim=1)
            NT = 0
             if(associated(pDatTo%DDATA2) ) then
                NT = size(pDatTo%DDATA2, dim=1)
             end if
             allocate(DDATA2(NT+NF, size(pDatFrom%DDATA2, dim=2)))

             if(NT .gt. 0 ) then
                DDATA2(1:NT,:) = pDatTo%DDATA2(1:NT,:)
                deallocate(pDatTo%DDATA2)
             end if
             DDATA2(NT+1:NT+NF,:) = pDatFrom%DDATA2(1:NF,:)
             PDatTo%DDATA2 => DDATA2

          end if
          return
  end subroutine Merge0_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine MergeByTag_DataPad(Tag, From, To)
  !***  PURPOSE:   to delete an element for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     character*(*)             ::Tag
     type(DataPad), pointer    ::From, To
     !--- local
     type(DataPad), pointer::pDatF, pDatT

          call GetDataPByTag_DataPad(Tag, From, pDatF)
          if(.not.associated(pDatF) ) then
             return
          end if

          call GetDataPByTag_DataPad(Tag, To, pDatT)
          if(.not. associated(pDatT)) then
             if(associated(pDatF%IDATA1)) then
                call New_DataPad(Tag, "I", 0, To)

             else if(associated(pDatF%IDATA2)) then
                call New_DataPad(Tag, "I", (/0,0/), To)

             else if(associated(pDatF%DDATA1)) then
                call New_DataPad(Tag, "D", 0, To)

             else if(associated(pDatF%DDATA2)) then
                call New_DataPad(Tag, "D", (/0,0/), To)
             end if
          end if
          call Merge0_DataPad(pDatF, pDatT)
          return
  end subroutine MergeByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine MergeAll_DataPad(From, To)
  !***  PURPOSE:   to delete an element for data array
  !     INPUT:    DatList
  !     OUTPUT    YES
  !
  implicit none
     !--- dummy varioables
     type(DataPad), pointer    ::From, To
     !--- local
     type(DataPad), pointer::pDatF
     integer I

          do I=1, NumberofData_DataPad(From)
             call GetDataPByInd_DataPad(I, From, pDatF)
             call MergeByTag_DataPad(pDatF%Tag, From, To)
          end do
          return
  end subroutine MergeAll_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI1EQ_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%IDATA1)) then
             N = size(pDat%IDATA1)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%IDATA1(I) .eq. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckI1EQ_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI1GT_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%IDATA1)) then
             N = size(pDat%IDATA1)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%IDATA1(I) .gt. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckI1GT_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI1GE_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%IDATA1)) then
             N = size(pDat%IDATA1)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%IDATA1(I) .ge. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckI1GE_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckI2EQ_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Col,     the col of 2-dimension data
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%IDATA2)) then
             N = size(pDat%IDATA2)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%IDATA2(I,Col) .eq. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckI2EQ_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI2GT_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Col,     the col of 2-dimension data
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%IDATA2)) then
             N = size(pDat%IDATA2)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%IDATA2(I,Col) .gt. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckI2GT_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI2GE_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Col,     the col of 2-dimension data
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%IDATA2)) then
             N = size(pDat%IDATA2)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%IDATA2(I,Col) .ge. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckI2GE_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckD1EQ_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%DDATA1)) then
             N = size(pDat%DDATA1)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%DDATA1(I) .eq. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckD1EQ_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckD1GT_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%DDATA1)) then
             N = size(pDat%DDATA1)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%DDATA1(I) .gt. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckD1GT_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckD1GE_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%DDATA1)) then
             N = size(pDat%DDATA1)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%DDATA1(I) .ge. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckD1GE_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckD2EQ_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Col,     the col of 2-dimension data
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%DDATA2)) then
             N = size(pDat%DDATA2)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%DDATA2(I,Col) .eq. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckD2EQ_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckD2GT_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Col,     the col of 2-dimension data
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%DDATA2)) then
             N = size(pDat%DDATA2)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%DDATA2(I,Col) .gt. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckD2GT_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckD2GE_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data larger than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Col,     the col of 2-dimension data
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local
     type(DataPad), pointer::pDat
     integer::I, N


          call GetDataPByTag_DataPad(Tag, DatList, pDat)
          if(.not. associated(pDat)) then
              return
          end if

          if(associated(pDat%DDATA2)) then
             N = size(pDat%DDATA2)
             do I=1, N
                if(Flag(I) .le. 0) cycle

                if(PDat%DDATA2(I,Col) .ge. Val) then
                   Flag(I) = 1
                else
                   Flag(I) = 0
                end if
             end do
          end if
          return
  end subroutine CheckD2GE_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI1LT_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckI1GT_ByTag_DataPad(Tag, DatList, -Val, Flag)
          return
  end subroutine CheckI1LT_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckI1LE_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckI1GE_ByTag_DataPad(Tag, DatList, -Val, Flag)
          return
  end subroutine CheckI1LE_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckI2LT_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckI2GT_ByTag_DataPad(Tag, DatList, Col, -Val, Flag)
          return
  end subroutine CheckI2LT_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckI2LE_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     integer,       intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckI2GE_ByTag_DataPad(Tag, DatList, Col, -Val, Flag)
          return
  end subroutine CheckI2LE_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckD1LT_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckD1GT_ByTag_DataPad(Tag, DatList, -Val, Flag)
          return
  end subroutine CheckD1LT_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckD1LE_ByTag_DataPad(Tag, DatList, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckD1GE_ByTag_DataPad(Tag, DatList, -Val, Flag)
          return
  end subroutine CheckD1LE_ByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine CheckD2LT_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckD2GT_ByTag_DataPad(Tag, DatList, Col, -Val, Flag)
          return
  end subroutine CheckD2LT_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine CheckD2LE_ByTag_DataPad(Tag, DatList, Col, Val, Flag)
  !***  PURPOSE:   to check if the data smaller than value
  !
  !     INPUT:    Tag,     the tag of the data to be checked
  !               DatList, the data list
  !               Val,     the value to compare to
  !               Flag,    the check results
  !
  !     OUTPUT    Flag,    the check results
  !
  !
  implicit none
     !--- dummy varioables
     character*(*), intent(in)         ::Tag
     type(DataPad), intent(in), pointer::DatList
     integer,       intent(in)         ::Col
     real(KINDDF),  intent(in)         ::Val
     integer, dimension(:)             ::Flag
     !--- local

          call CheckD2GE_ByTag_DataPad(Tag, DatList, Col, -Val, Flag)
          return
  end subroutine CheckD2LE_ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine AccumulateIByTag_DataPad(Tag, DatList, Val, Order)
  !***  PURPOSE:   to accumulated the data of TAG
  !     INPUT:    Tag,     the data tag
  !               DatList, the data list
  !               Val,     the data to be added to Tag data
  !               Order,   the optional order array
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  !               NOTE: the size of the array Val should be GE the size of
  !                     the data size. the size of Order is the same size
  !                     of the TAG data.
  !                     This routine is desinged so in the consideration that
  !                     the data generated by GPU with multi[le boxes, the
  !                     size of array on GPU is larger than the arraies in a
  !                     simulation box. And the indice of the arraies have been
  !                     reordered.
  !
  implicit none
     !--- dummy varioables
     character*(*)                  ::Tag
     type(DataPad), pointer         ::DatList
     integer, dimension(:)          ::Val
     integer, optional, dimension(:)::Order
     !--- local
     type(DataPad), pointer::pDat
     integer, dimension(:), pointer::Pad
     integer::N, I

          call GetDataPByTag_DataPad(Tag, DatList, pDat)

          if(.not.associated(pDat) ) then
             if(present(Order)) then
                call New_DataPad(Tag, size(Order), Pad, DatList)
             else
                call New_DataPad(Tag, size(Val), Pad, DatList)
             end if
             call GetDataPByTag_DataPad(Tag, DatList, pDat)
             pDat%IDATA1 = 0
          end if

          N = size(pDat%IDATA1)
          if(present(Order)) then
             pDat%IDATA1(1:N) = pDat%IDATA1(1:N) + Val(Order(1:N))
          else
             pDat%IDATA1(1:N) = pDat%IDATA1(1:N) + Val(1:N)
          end if
          return
  end subroutine AccumulateIByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine AccumulateDByTag_DataPad(Tag, DatList, Val, Order)
  !***  PURPOSE:   to accumulated the data of TAG
  !     INPUT:    Tag,     the data tag
  !               DatList, the data list
  !               Val,     the data to be added to Tag data
  !               Order,   the optional order array
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  !               NOTE: the size of the array Val should be GE the size of
  !                     the data size. the size of Order is the same size
  !                     of the TAG data.
  !                     This routine is desinged so in the consideration that
  !                     the data generated by GPU with multi[le boxes, the
  !                     size of array on GPU is larger than the arraies in a
  !                     simulation box. And the indice of the arraies have been
  !                     reordered.
  !
  implicit none
     !--- dummy varioables
     character*(*)                  ::Tag
     type(DataPad), pointer         ::DatList
     real(KINDDF), dimension(:)     ::Val
     integer, optional, dimension(:)::Order
     !--- local
     type(DataPad), pointer::pDat
     real(KINDDF),  dimension(:), pointer::Pad
     integer::N, I

          call GetDataPByTag_DataPad(Tag, DatList, pDat)

          if(.not.associated(pDat) ) then
             if(present(Order)) then
                call New_DataPad(Tag, size(Order), Pad, DatList)
             else
                call New_DataPad(Tag, size(Val), Pad, DatList)
             end if
             call GetDataPByTag_DataPad(Tag, DatList, pDat)
             pDat%DDATA1 = 0.D0
          end if

          N = size(pDat%DDATA1)
          if(present(Order)) then
             pDat%DDATA1(1:N) = pDat%DDATA1(1:N) + Val(Order(1:N))
          else
             pDat%DDATA1(1:N) = pDat%DDATA1(1:N) + Val(1:N)
          end if
          return
  end subroutine AccumulateDByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine AccumulateI2ByTag_DataPad(Tag, DatList, Val, Order)
  !***  PURPOSE:   to accumulated the data of TAG
  !     INPUT:    Tag,     the data tag
  !               DatList, the data list
  !               Val,     the data to be added to Tag data
  !               Order,   the optional order array
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  !               NOTE: the size of the array Val should be GE the size of
  !                     the data size. the size of Order is the same size
  !                     of the TAG data.
  !                     This routine is desinged so in the consideration that
  !                     the data generated by GPU with multi[le boxes, the
  !                     size of array on GPU is larger than the arraies in a
  !                     simulation box. And the indice of the arraies have been
  !                     reordered.
  !
  implicit none
     !--- dummy varioables
     character*(*)                  ::Tag
     type(DataPad), pointer         ::DatList
     integer, dimension(:,:)        ::Val
     integer, optional, dimension(:)::Order
     !--- local
     type(DataPad), pointer::pDat
     integer, dimension(:,:), pointer::Pad
     integer::N, I

          call GetDataPByTag_DataPad(Tag, DatList, pDat)

          if(.not.associated(pDat) ) then
             if(present(Order)) then
                call New_DataPad(Tag, (/size(Order), size(Val,dim=2)/), Pad, DatList)
             else
                call New_DataPad(Tag, (/size(Val,dim=1), size(Val,dim=2)/), Pad, DatList)
             end if
             call GetDataPByTag_DataPad(Tag, DatList, pDat)
             pDat%IDATA2 = 0
          end if

          N = size(pDat%IDATA2, dim=1)
          if(present(Order)) then
             pDat%IDATA2(1:N,:) = pDat%IDATA2(1:N,:) + Val(Order(1:N),:)
          else
             pDat%IDATA2(1:N,:) = pDat%IDATA2(1:N,:) + Val(1:N,:)
          end if
          return
  end subroutine AccumulateI2ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine AccumulateD2ByTag_DataPad(Tag, DatList, Val, Order)
  !***  PURPOSE:   to accumulated the data of TAG
  !     INPUT:    Tag,     the data tag
  !               DatList, the data list
  !               Val,     the data to be added to Tag data
  !               Order,   the optional order array
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  !               NOTE: the size of the array Val should be GE the size of
  !                     the data size. the size of Order is the same size
  !                     of the TAG data.
  !                     This routine is desinged so in the consideration that
  !                     the data generated by GPU with multi[le boxes, the
  !                     size of array on GPU is larger than the arraies in a
  !                     simulation box. And the indice of the arraies have been
  !                     reordered.
  !
  implicit none
     !--- dummy varioables
     character*(*)                  ::Tag
     type(DataPad), pointer         ::DatList
     real(KINDDF), dimension(:,:)   ::Val
     integer, optional, dimension(:)::Order
     !--- local
     type(DataPad), pointer::pDat
     real(KINDDF), dimension(:,:), pointer::Pad
     integer::N, I

          call GetDataPByTag_DataPad(Tag, DatList, pDat)

          if(.not.associated(pDat) ) then
             if(present(Order)) then
                call New_DataPad(Tag, (/size(Order),size(Val,dim=2)/), Pad, DatList)
             else
                call New_DataPad(Tag, (/size(Val,dim=1),size(Val,dim=2)/), Pad, DatList)
             end if
             call GetDataPByTag_DataPad(Tag, DatList, pDat)
             pDat%DDATA2 = 0.D0
          end if

          N = size(pDat%DDATA2, dim=2)
          if(present(Order)) then
             pDat%DDATA2(1:N,:) = pDat%DDATA2(1:N,:) + Val(Order(1:N),:)
          else
             pDat%DDATA2(1:N,:) = pDat%DDATA2(1:N,:) + Val(1:N,:)
          end if
          return
  end subroutine AccumulateD2ByTag_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine Archive_DataPad(hFile, DatList)
  !***  PURPOSE:   to archive a data pad
  !
  !     INPUT:    hFile,      the I/O unist
  !               Ind,        the Index to be archive
  !               DatList,
  !
  !     OUTPUT
  implicit none
     !--- dummy varioables
     integer, intent(in)           ::hFile
     type(DataPad),         pointer::DatList
     !---- local
     character(len=mp_TAGLEN)::Tag
     integer::I, N, Dim(2)
     character::dtyp
     type(DataPad), pointer::pDat

             Tag = "PADDATA"
             N   = NumberofData_DataPad(DatList)
             write(hFile)Tag
             write(hFile)N
             do I=1, N
                call GetDataPByInd_DataPad(I, DatList, pDat)
                call Dim_DataPad(Dim, pDat)
                if(IsInteger_DataPad(pDat)) then
                   dtyp = 'I'
                else
                   dtyp = 'D'
                end if
                write(hFile) pDat%Tag
                write(hFile) dtyp
                write(hFile) Dim
                if(associated(pDat%IDATA1)) then
                   write(hFile) pDat%IDATA1

                else if(associated(pDat%IDATA2)) then
                   write(hFile) pDat%IDATA2

                else if(associated(pDat%DDATA1)) then
                   write(hFile) pDat%DDATA1

                else if(associated(pDat%DDATA2)) then
                    write(hFile) pDat%DDATA2
                end if

             end do
          return
  end subroutine Archive_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine Restore_DataPad(hFile, DatList)
  !***  PURPOSE:   to restore a data pad
  !
  !     INPUT:    hFile,      the I/O unist
  !               Ind,        the Index to be archive
  !               DatList,
  !
  !     OUTPUT
  implicit none
     !--- dummy varioables
     integer, intent(in)           ::hFile
     type(DataPad),         pointer::DatList
     !---- local
     character(len=mp_TAGLEN)::Tag
     integer                 ::I, N, Dim(2)
     character               ::dtyp
     integer,      dimension(:),   pointer::PadI
     real(KINDDF), dimension(:),   pointer::PadD
     integer,      dimension(:,:), pointer::PadI2
     real(KINDDF), dimension(:,:), pointer::PadD2


            call Release_DataPad(DatList)
            Read(hFile) Tag
            Read(hFile) N
            do I=1, N
               Read(hFile) Tag
               Read(hFile) dtyp
               Read(hFile) Dim

               if(Dim(2) .gt. 1) then
                  if(dtyp .eq. 'I' .or. dtyp .eq. 'i') then
                    call New_DataPad_I2(Tag, Dim, PadI2, DatList)
                    read(hFile)PadI2

                  else if(dtyp.eq. 'D' .or. dtyp .eq. 'd') then
                    call New_DataPad_D2(Tag, Dim, PadD2, DatList)
                    read(hFile)PadD2
                  end if
              else
                  if(dtyp .eq. 'I' .or. dtyp .eq. 'i') then
                    call New_DataPad_I(Tag, Dim(1), PadI, DatList)
                    read(hFile)PadI

                  else if(dtyp .eq. 'D' .or. dtyp .eq. 'd') then
                    call New_DataPad_D(Tag, Dim(1), PadD, DatList)
                    read(hFile)PadD
                  end if
              end if
            end do
          return
  end subroutine Restore_DataPad
 !****************************************************************************
 !****************************************************************************
  subroutine ResetByTag_DataPad(Tag, DatList)
  !***  PURPOSE:   to reset the data of TAG to zero
  !     INPUT:    Tag,      the data tag
  !               DatList,  the data list
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  implicit none
     !--- dummy varioables
     character*(*)                  ::Tag
     type(DataPad), pointer         ::DatList
     !--- local
     integer::IND

          call Ind_DataPad(Tag, DatList, IND)
          call ResetByInd_DataPad(IND, DatList)
          return
  end subroutine ResetByTag_DataPad
 !****************************************************************************

  !****************************************************************************
  subroutine ResetByInd_DataPad(Ind, DatList)
  !***  PURPOSE:  to reset the data of TAG to zero
  !     INPUT:    Ind,      the index of the data 
  !               DatList,  the data list
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  implicit none
     !--- dummy varioables
     integer                  ::Ind
     type(DataPad), pointer   ::DatList
     !--- local
     type(DataPad), pointer::pDat

          call GetDataPByInd_DataPad(Ind, DatList, pDat)
          if(.not.associated(pDat)) return

          if(associated(pDat%IDATA1)) then
             pDat%IDATA1 = 0
          end if

          if(associated(pDat%IDATA2)) then
             pDat%IDATA2 = 0
          end if

          if(associated(pDat%DDATA1)) then
             pDat%DDATA1 = 0.D0
          end if
          if(associated(pDat%DDATA2)) then
             pDat%DDATA2 = 0.D0
          end if
          return
  end subroutine ResetByInd_DataPad
 !****************************************************************************

 !****************************************************************************
  subroutine SetValByTag_DataPad_1I_A(Tag, DatList, From, To, Val)
  !***  PURPOSE:   to reset the data of TAG to zero
  !     INPUT:     Tag,      the data tag
  !                DatList,  the data list
  !               from, to, the index of the segement
  !               Val,      the val to be set 
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  implicit none
     !--- dummy varioables
     character*(*)            ::Tag
     type(DataPad), pointer   ::DatList
     integer                  ::From, To
     integer                  ::Val 
    !--- local
     integer::IND

          call Ind_DataPad(Tag, DatList, IND)
          call SetValByInd_DataPad_1I_A(IND, DatList, From, To, Val)
          return
  end subroutine SetValByTag_DataPad_1I_A
  !****************************************************************************

  !****************************************************************************
  subroutine SetValByInd_DataPad_1I_A(Ind, DatList, From, To, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Ind,      the data tag
   !               DatList,  the data list
   !               from, to, the index of the segement
   !               Val,      the val to be set 
   !     OUTPUT    DatList, the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      integer                  ::Ind
      type(DataPad), pointer   ::DatList
      integer                  ::From, To
      integer                  ::Val 
      !--- local
      type(DataPad), pointer::pDat
 
           call GetDataPByInd_DataPad(Ind, DatList, pDat)
           if(.not.associated(pDat)) return
 
           if(associated(pDat%IDATA1)) then
              pDat%IDATA1(From:To) = Val
           end if
 
           if(associated(pDat%DDATA1)) then
              pDat%DDATA1(From:To) = dble(Val)
           end if
           return
   end subroutine SetValByInd_DataPad_1I_A
  !****************************************************************************

  !****************************************************************************
  subroutine SetValByTag_DataPad_1I_B(Tag, DatList, I, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Tag,      the data tag
   !               DatList,  the data list
   !               from, to, the index of the segement
   !               Val,      the val to be set 
   !     OUTPUT    DatList, the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      character*(*)            ::Tag
      type(DataPad), pointer   ::DatList
      integer                  ::I
      integer                  ::Val 
      !--- local
      integer::IND

          call Ind_DataPad(Tag, DatList, IND)
          call SetValByInd_DataPad_1I_A(Ind, DatList, I, I, Val) 
          return
   end subroutine SetValByTag_DataPad_1I_B
  !****************************************************************************   

  !****************************************************************************
  subroutine SetValByInd_DataPad_1I_B(Ind, DatList, I, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Ind,      the data tag
   !               DatList,  the data list
   !               I,         the index of the segement
   !               Val,      the val to be set 
   !     OUTPUT    DatList, the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      integer                  ::Ind
      type(DataPad), pointer   ::DatList
      integer                  ::I
      integer                  ::Val 
      !--- local
      type(DataPad), pointer::pDat
 
          call SetValByInd_DataPad_1I_A(Ind, DatList, I, I, Val) 
           return
   end subroutine SetValByInd_DataPad_1I_B
  !****************************************************************************

 !****************************************************************************
  subroutine SetValByTag_DataPad_1D_A(Tag, DatList, From, To, Val)
  !***  PURPOSE:   to reset the data of TAG to zero
  !     INPUT:     Tag,      the data tag
  !                DatList,  the data list
  !               from, to, the index of the segement
  !               Val,      the val to be set 
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  implicit none
     !--- dummy varioables
     character*(*)            ::Tag
     type(DataPad), pointer   ::DatList
     integer                  ::From, To
     real(KINDDF)             ::Val 
    !--- local
     integer::IND

          call Ind_DataPad(Tag, DatList, IND)
          call SetValByInd_DataPad_1D_A(IND, DatList, From, To, Val)
          return
  end subroutine SetValByTag_DataPad_1D_A
  !****************************************************************************

  !****************************************************************************
  subroutine SetValByInd_DataPad_1D_A(Ind, DatList, From, To, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Ind,      the data tag
   !               DatList,  the data list
   !               from, to, the index of the segement
   !               Val,      the val to be set 
   !     OUTPUT    DatList, the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      integer                  ::Ind
      type(DataPad), pointer   ::DatList
      integer                  ::From, To
      real(KINDDF)              ::Val 
      !--- local
      type(DataPad), pointer::pDat
 
           call GetDataPByInd_DataPad(Ind, DatList, pDat)
           if(.not.associated(pDat)) return
 
           if(associated(pDat%IDATA1)) then
              pDat%IDATA1(From:To) = Val
           end if
 
           if(associated(pDat%DDATA1)) then
              pDat%DDATA1(From:To) = Val
           end if
           return
   end subroutine SetValByInd_DataPad_1D_A
  !**************************************************************************** 
   
  !****************************************************************************
  subroutine SetValByTag_DataPad_1D_B(Tag, DatList, I, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Tag,      the data tag
   !               DatList,  the data list
   !               from, to, the index of the segement
   !               Val,      the val to be set 
   !     OUTPUT    DatList, the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      character*(*)            ::Tag
      type(DataPad), pointer   ::DatList
      integer                  ::I
      real(KINDDF)             ::Val 
      !--- local
      integer::IND

          call Ind_DataPad(Tag, DatList, IND)
          call SetValByInd_DataPad_1D_A(Ind, DatList, I, I, Val) 
          return
   end subroutine SetValByTag_DataPad_1D_B
  !****************************************************************************   

  !****************************************************************************
  subroutine SetValByInd_DataPad_1D_B(Ind, DatList, I, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Ind,      the data tag
   !               DatList,  the data list
   !               I,         the index of the segement
   !               Val,      the val to be set 
   !     OUTPUT    DatList, the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      integer                  ::Ind
      type(DataPad), pointer   ::DatList
      integer                  ::I
      real(KINDDF)             ::Val 
      !--- local
      type(DataPad), pointer::pDat
 
          call SetValByInd_DataPad_1D_A(Ind, DatList, I, I, Val) 
          return
   end subroutine SetValByInd_DataPad_1D_B
  !****************************************************************************  
   
 !****************************************************************************
  subroutine SetValByTag_DataPad_2D_A(Tag, DatList, At, Val)
  !***  PURPOSE:   to reset the data of TAG to zero
  !     INPUT:     Tag,      the data tag
  !                DatList,  the data list
  !                At,       the index of the element
  !                Val,      the val to be set 
  !     OUTPUT    DatList, the DadaPad with its data of TAG updated
  !
  implicit none
     !--- dummy varioables
     character*(*)             ::Tag
     type(DataPad), pointer    ::DatList
     integer                   ::At
     real(KINDDF), dimension(:)::Val 
    !--- local
     type(DataPad), pointer::pDat

             call GetDataPByTag_DataPad(Tag, DatList, pDat)
             if(.not.associated(pDat)) return
     
             if(associated(pDat%IDATA2)) then
                pDat%IDATA2(At, :) = Val(:)
             end if
     
             if(associated(pDat%DDATA2)) then
                pDat%DDATA2(At, :) = Val(:)
             end if          
          return
  end subroutine SetValByTag_DataPad_2D_A
  !****************************************************************************

  !****************************************************************************
  subroutine SetValByInd_DataPad_2D_A(Ind, DatList, At, Val)
   !***  PURPOSE:   to set the data sehgement of Ind to givenvalue
   !     INPUT:    Ind,      the data tag
   !               DatList,  the data list
   !               At,       the index of the element
   !               Val,      the val to be set 
   !     OUTPUT    DatList,  the DadaPad with its data of TAG updated
   !
   implicit none
      !--- dummy varioables
      integer                  ::Ind
      type(DataPad), pointer   ::DatList
      integer                  ::At
      real(KINDDF), dimension(:)::Val 
      !--- local
      type(DataPad), pointer::pDat
 
           call GetDataPByInd_DataPad(Ind, DatList, pDat)
           if(.not.associated(pDat)) return
 
           if(associated(pDat%IDATA2)) then
              pDat%IDATA2(At, :) = Val(:)
           end if
 
           if(associated(pDat%DDATA2)) then
              pDat%DDATA2(At, :) = Val(:)
           end if
         return
   end subroutine SetValByInd_DataPad_2D_A
  !**************************************************************************** 

  end module MSM_TYPEDEF_DataPad
 !****************************************************************************************
