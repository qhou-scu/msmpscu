  ! ***  DESCRIPTION: this file provides a number of template routines for allocating memeory, 
  !                   which are parts of file MCLIB_Utilities.F90 originally written by Zhai Lei.
  !                   this file can be included in a file with ARRAYDEV defined as empty for CPU or
  !                   device for GPU
  !                   see also: MSM_MemAllocator.F90
  !                             MSM_MemAllocator_GPU.F90
  !                   
  !                  ______________________________________________________
  !
  ! **** HOSTORY:
  !       * Sep.    2018(Zhai Lei): original version of MCLIB_Utilities.F90 
  !                                 
  !
  !       * Nov.    2018(HOU Qing): Seperate from MCLIB_Utilities.F90, 
  !                                 define macro for allocate memory on host
  !                                 and GPU
  !            
  !                                 
  
  !*************************************************************
  subroutine AllocateOneDimi(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    integer,ARRAYDEV dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars--- 
    integer::istat
    !---Body---
    call DeAllocateOneDimi(Array,Name)

    allocate(Array(Length),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateOneDimi

  !*************************************************************
  subroutine AllocateOneDimr(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF),ARRAYDEV dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateOneDimr(Array,Name)

    allocate(Array(Length),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateOneDimr

  !*************************************************************
  subroutine AllocateOneDimd(Array,Length,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF),ARRAYDEV dimension(:),allocatable::Array
    integer,intent(in)::Length
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateOneDimd(Array,Name)

    allocate(Array(Length),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateOneDimd


  !*************************************************************
  subroutine AllocateTwoDimi(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    integer,ARRAYDEV dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimi(Array,Name)

    allocate(Array(LengthX,LengthY),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateTwoDimi

  !*************************************************************
  subroutine AllocateTwoDimr(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF),ARRAYDEV dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimr(Array,Name)

    allocate(Array(LengthX,LengthY),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateTwoDimr

  !*************************************************************
  subroutine AllocateTwoDimd(Array,LengthX,LengthY,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF),ARRAYDEV dimension(:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateTwoDimd(Array,Name)

    allocate(Array(LengthX,LengthY),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateTwoDimd


  !*************************************************************
  subroutine AllocateThreeDimi(Array,LengthX,LengthY,LengthZ,Name)
    implicit none
    !---Dummy Vars---
    integer,ARRAYDEV dimension(:,:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    integer,intent(in)::LengthZ
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateThreeDimi(Array,Name)

    allocate(Array(LengthX,LengthY,LengthZ),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateThreeDimi

  !*************************************************************
  subroutine AllocateThreeDimr(Array,LengthX,LengthY,LengthZ,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF),ARRAYDEV dimension(:,:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    integer,intent(in)::LengthZ
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateThreeDimr(Array,Name)

    allocate(Array(LengthX,LengthY,LengthZ),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateThreeDimr

  !*************************************************************
  subroutine AllocateThreeDimd(Array,LengthX,LengthY,LengthZ,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF),ARRAYDEV dimension(:,:,:),allocatable::Array
    integer,intent(in)::LengthX
    integer,intent(in)::LengthY
    integer,intent(in)::LengthZ
    character(*)::Name
    !---Dummy Vars---
    integer::istat
    !---Body---
    call DeAllocateThreeDimd(Array,Name)

    allocate(Array(LengthX,LengthY,LengthZ),STAT=istat)
    if(istat /=0) then
        write(*,*) "MSMPSCU error: The Array :",Name,"allocate Failed !"
        pause
        stop
    end if

    return
  end subroutine AllocateThreeDimd


  !*************************************************************
  subroutine DeAllocateOneDimi(Array,Name)
    implicit none
    !---Dummy Vars---
    integer,ARRAYDEV dimension(:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateOneDimi


  !*************************************************************
  subroutine DeAllocateOneDimr(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF),ARRAYDEV dimension(:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateOneDimr

  !*************************************************************
  subroutine DeAllocateOneDimd(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF),ARRAYDEV dimension(:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateOneDimd

  !*************************************************************
  subroutine DeAllocateTwoDimi(Array,Name)
    implicit none
    !---Dummy Vars---
    integer,ARRAYDEV dimension(:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateTwoDimi

  !*************************************************************
  subroutine DeAllocateTwoDimr(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF),ARRAYDEV dimension(:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateTwoDimr

  !*************************************************************
  subroutine DeAllocateTwoDimd(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF),ARRAYDEV dimension(:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateTwoDimd


  !*************************************************************
  subroutine DeAllocateThreeDimi(Array,Name)
    implicit none
    !---Dummy Vars---
    integer,ARRAYDEV dimension(:,:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateThreeDimi

  !*************************************************************
  subroutine DeAllocateThreeDimr(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF),ARRAYDEV dimension(:,:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateThreeDimr

  !*************************************************************
  subroutine DeAllocateThreeDimd(Array,Name)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF),ARRAYDEV dimension(:,:,:),allocatable::Array
    character(*)::Name
    !---Local Vars---
    integer::istat
    !---Body---

    if(allocated(Array)) then
        deallocate(Array,STAT=istat)
        if(istat /=0) then
            write(*,*) "MSMPSCU error: The Array :",Name,"Deallocate Failed !"
            pause
            stop
        end if
    end if

    return
  end subroutine DeAllocateThreeDimd


  !*****************************************
  subroutine ResizeArrayi_OneDim(TheArray,NewSize)
    implicit none
    !---Dummy Vars---
    integer, ARRAYDEV dimension(:), allocatable,intent(inout)::TheArray
    integer,intent(in)::NewSize
    !---Local Vars---
    integer::OldSize(1)
    integer, ARRAYDEV dimension(:), allocatable::tempArray
    integer::istat
    !---Body----
    OldSize = shape(TheArray)

    if(OldSize(1) .ne. NewSize) then

        ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
        !       the tempArray need not need to be allocated, the compiler would
        !       allocate the array automatic that is marked as "allocatable" based on the
        !       assigned array's size
        call AllocateOneDimi(tempArray,OldSize(1),"tempArray")
             

        tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

        call AllocateOneDimi(TheArray,NewSize,"TheArray")

        TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[0])

        call DeAllocateOneDimi(tempArray,"tempArray")

    end if


    return
  end subroutine ResizeArrayi_OneDim

  !*****************************************
  subroutine ResizeArrayr_OneDim(TheArray,NewSize)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF), ARRAYDEV dimension(:), allocatable,intent(inout)::TheArray
    integer,intent(in)::NewSize
    !---Local Vars---
    integer::OldSize(1)
    real(kind=KINDSF), ARRAYDEV dimension(:), allocatable::tempArray
    integer::istat
    real(kind=KINDSF)::zero_Single
    !---Body----
    zero_Single = 0.D0
    OldSize = shape(TheArray)

    if(OldSize(1) .ne. NewSize) then

        ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
        !       the tempArray need not need to be allocated, the compiler would
        !       allocate the array automatic that is marked as "allocatable" based on the
        !       assigned array's size
        call AllocateOneDimr(tempArray,OldSize(1),"tempArray")

        tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

        call AllocateOneDimr(TheArray,NewSize,"TheArray")

        TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[zero_Single])

        call DeAllocateOneDimr(tempArray,"tempArray")

    end if


    return
  end subroutine ResizeArrayr_OneDim

  !*****************************************
  subroutine ResizeArrayd_OneDim(TheArray,NewSize)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF), ARRAYDEV dimension(:), allocatable,intent(inout)::TheArray
    integer,intent(in)::NewSize
    !---Local Vars---
    integer::OldSize(1)
    real(kind=KINDDF), ARRAYDEV dimension(:), allocatable::tempArray
    integer::istat
    real(kind=KINDDF)::zero_Double
    !---Body----
    zero_Double = 0.D0
    OldSize = shape(TheArray)

    if(OldSize(1) .ne. NewSize) then

        ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
        !       the tempArray need not need to be allocated, the compiler would
        !       allocate the array automatic that is marked as "allocatable" based on the
        !       assigned array's size
        call AllocateOneDimd(tempArray,OldSize(1),"tempArray")

        tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

        call AllocateOneDimd(TheArray,NewSize,"TheArray")

        TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[zero_Double])

        call DeAllocateOneDimd(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayd_OneDim

  !**********************************************************
  subroutine ResizeArrayi_TwoDim(TheArray,NewX,NewY)
    implicit none
    !---Dummy Vars---
    integer, ARRAYDEV dimension(:,:), allocatable::TheArray
    integer::NewX
    integer::NewY
    !---Local Vars---
    integer::oldShape(2)
    integer, ARRAYDEV dimension(:,:), allocatable::tempArray
    integer::istat
    !---Body---

    oldShape = shape(TheArray)

    if(oldShape(1) .ne. NewX .and. oldShape(2) .ne. NewY) then


        if(NewX .LE. oldShape(1)) then

            if(NewY .LE. oldShape(2)) then

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateTwoDimi(tempArray,NewX,NewY,"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(1:NewX,1:NewY)],SHAPE=[NewX,NewY])

                call AllocateTwoDimi(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY])
            else
                ! note: Here, we use the "reshape" to copy the array. The main reason for this is that we need to consider the
                !       case that compiler option "-Mallocatable=03" is used, for this compiler options, base on our test, if we
                !       try to copy the array1 to array2 by the form of array2 = array1 and the array2 is marked as "allocatable",
                !       if these two array own same shape, the array2 would be reshaped to the size saming with array1, for instant,
                !       if array1 with size of 1xs1 and array with size of 1xs2 and s1 .ne. s2, in this case, after the operation,
                !       the array2 would be changed to the size of 1xs1 !!!!!
                !       Thus, We apply the "reshape" here, reshape([source1,source2,...],[dim1,dim2,dim3...]) to copy the value and
                !       ensure that array2 would keep the size of 1xs2 after the copying operation.
                !       If s1< s2, we use the form of array2 = reshape([array1],[s2]) here,infact ,there are tow steps are
                !       implicated in this opertion:
                !           1.The reshape operation would construct return an "tempArray" with size 1xs2, and the arrays elements
                !             from "tempArray(1)" to "tempArray(s1)" would be filled by array1
                !           2, array2 = "tempArray" would be executed,
                !       So, based on this startegy, we can keep the size for array2 and copy value from array1 at the same time.
                !       (If s1> s2, we should use array2 = reshape([array1(1:s2)],[s2]) to avoid overflow)
                !
                !       Additional, if s1< s2, the pgfortran compiler would do two things:
                !            1: copy the contents with the size of s1 from to array1 to array2
                !            2: the reminded contents with the sizeof (s2 -s1 +1) in array2
                !               would be assigned to some value automaticlly(the reminded
                !               contents would be filled based on the value in array1, for
                !               instant, if the value in array1 is 1,2,3...., the remineded
                !               contents would be assigned to s1+1,s1+2,....,s2)
                !       In fact, we hope that the value in array2 would be array(1:s1),0(s1+1:s2) ,which means the remined value
                !       shoud not be assigned by some values, thus we should do the additional operation: array2(s1+1:s2) = 0

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateTwoDimi(tempArray,NewX,oldShape(2),"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(NewX,oldShape(2))],SHAPE=[NewX,oldShape(2)])


                call AllocateTwoDimi(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[0],ORDER=[1,2])

            end if

        else

            if(NewY .LE. oldShape(2)) then

                call AllocateTwoDimi(tempArray,NewY,oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[NewY,oldShape(1)],PAD=[0],ORDER=[2,1])

                call AllocateTwoDimi(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[0],ORDER=[2,1])

            else

                call AllocateTwoDimi(tempArray,oldShape(2),oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[oldShape(2),oldShape(1)],PAD=[0],ORDER=[2,1])

                call AllocateTwoDimi(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[reshape(SOURCE=[tempArray],&
                                                   SHAPE=[NewX,oldShape(2)],&
                                                   PAD=[0],&
                                                   ORDER=[2,1])],&
                                   SHAPE=[NewX,NewY],&
                                   PAD=[0],&
                                   ORDER=[1,2])

            end if

        end if

        call DeAllocateTwoDimi(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayi_TwoDim

  !**********************************************************
  subroutine ResizeArrayr_TwoDim(TheArray,NewX,NewY)
    implicit none
    !---Dummy Vars---
    real(kind=KINDSF), ARRAYDEV dimension(:,:), allocatable::TheArray
    integer::NewX
    integer::NewY
    !---Local Vars---
    integer::oldShape(2)
    real(kind=KINDSF), ARRAYDEV dimension(:,:), allocatable::tempArray
    real(kind=KINDSF)::zero_Single
    integer::istat
    !---Body---

    zero_Single = 0.D0

    oldShape = shape(TheArray)

    if(oldShape(1) .ne. NewX .and. oldShape(2) .ne. NewY) then


        if(NewX .LE. oldShape(1)) then

            if(NewY .LE. oldShape(2)) then

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateTwoDimr(tempArray,NewX,NewY,"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(1:NewX,1:NewY)],SHAPE=[NewX,NewY])

                call AllocateTwoDimr(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY])
            else
                ! note: Here, we use the "reshape" to copy the array. The main reason for this is that we need to consider the
                !       case that compiler option "-Mallocatable=03" is used, for this compiler options, base on our test, if we
                !       try to copy the array1 to array2 by the form of array2 = array1 and the array2 is marked as "allocatable",
                !       if these two array own same shape, the array2 would be reshaped to the size saming with array1, for instant,
                !       if array1 with size of 1xs1 and array with size of 1xs2 and s1 .ne. s2, in this case, after the operation,
                !       the array2 would be changed to the size of 1xs1 !!!!!
                !       Thus, We apply the "reshape" here, reshape([source1,source2,...],[dim1,dim2,dim3...]) to copy the value and
                !       ensure that array2 would keep the size of 1xs2 after the copying operation.
                !       If s1< s2, we use the form of array2 = reshape([array1],[s2]) here,infact ,there are tow steps are
                !       implicated in this opertion:
                !           1.The reshape operation would construct return an "tempArray" with size 1xs2, and the arrays elements
                !             from "tempArray(1)" to "tempArray(s1)" would be filled by array1
                !           2, array2 = "tempArray" would be executed,
                !       So, based on this startegy, we can keep the size for array2 and copy value from array1 at the same time.
                !       (If s1> s2, we should use array2 = reshape([array1(1:s2)],[s2]) to avoid overflow)
                !
                !       Additional, if s1< s2, the pgfortran compiler would do two things:
                !            1: copy the contents with the size of s1 from to array1 to array2
                !            2: the reminded contents with the sizeof (s2 -s1 +1) in array2
                !               would be assigned to some value automaticlly(the reminded
                !               contents would be filled based on the value in array1, for
                !               instant, if the value in array1 is 1,2,3...., the remineded
                !               contents would be assigned to s1+1,s1+2,....,s2)
                !       In fact, we hope that the value in array2 would be array(1:s1),0(s1+1:s2) ,which means the remined value
                !       shoud not be assigned by some values, thus we should do the additional operation: array2(s1+1:s2) = 0

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateTwoDimr(tempArray,NewX,oldShape(2),"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(NewX,oldShape(2))],SHAPE=[NewX,oldShape(2)])

                call AllocateTwoDimr(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Single],ORDER=[1,2])

            end if

        else

            if(NewY .LE. oldShape(2)) then

                call AllocateTwoDimr(tempArray,NewY,oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[NewY,oldShape(1)],PAD=[zero_Single],ORDER=[2,1])

                call AllocateTwoDimr(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Single],ORDER=[2,1])

            else
                call AllocateTwoDimr(tempArray,oldShape(2),oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[oldShape(2),oldShape(1)],PAD=[zero_Single],ORDER=[2,1])

                call AllocateTwoDimr(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[reshape(SOURCE=[tempArray],&
                                                   SHAPE=[NewX,oldShape(2)],&
                                                   PAD=[zero_Single],&
                                                   ORDER=[2,1])],&
                                   SHAPE=[NewX,NewY],&
                                   PAD=[zero_Single],&
                                   ORDER=[1,2])

            end if

        end if

        call DeAllocateTwoDimr(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayr_TwoDim


  !**********************************************************
  subroutine ResizeArrayd_TwoDim(TheArray,NewX,NewY)
    implicit none
    !---Dummy Vars---
    real(kind=KINDDF), ARRAYDEV dimension(:,:), allocatable::TheArray
    integer::NewX
    integer::NewY
    !---Local Vars---
    integer::oldShape(2)
    real(kind=KINDDF), ARRAYDEV dimension(:,:), allocatable::tempArray
    real(kind=KINDDF)::zero_Double
    integer::istat
    !---Body---

    zero_Double = 0.D0

    oldShape = shape(TheArray)

    if(oldShape(1) .ne. NewX .and. oldShape(2) .ne. NewY) then


        if(NewX .LE. oldShape(1)) then

            if(NewY .LE. oldShape(2)) then

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateTwoDimd(tempArray,NewX,NewY,"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(1:NewX,1:NewY)],SHAPE=[NewX,NewY])

                call AllocateTwoDimd(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY])
            else
                ! note: Here, we use the "reshape" to copy the array. The main reason for this is that we need to consider the
                !       case that compiler option "-Mallocatable=03" is used, for this compiler options, base on our test, if we
                !       try to copy the array1 to array2 by the form of array2 = array1 and the array2 is marked as "allocatable",
                !       if these two array own same shape, the array2 would be reshaped to the size saming with array1, for instant,
                !       if array1 with size of 1xs1 and array with size of 1xs2 and s1 .ne. s2, in this case, after the operation,
                !       the array2 would be changed to the size of 1xs1 !!!!!
                !       Thus, We apply the "reshape" here, reshape([source1,source2,...],[dim1,dim2,dim3...]) to copy the value and
                !       ensure that array2 would keep the size of 1xs2 after the copying operation.
                !       If s1< s2, we use the form of array2 = reshape([array1],[s2]) here,infact ,there are tow steps are
                !       implicated in this opertion:
                !           1.The reshape operation would construct return an "tempArray" with size 1xs2, and the arrays elements
                !             from "tempArray(1)" to "tempArray(s1)" would be filled by array1
                !           2, array2 = "tempArray" would be executed,
                !       So, based on this startegy, we can keep the size for array2 and copy value from array1 at the same time.
                !       (If s1> s2, we should use array2 = reshape([array1(1:s2)],[s2]) to avoid overflow)
                !
                !       Additional, if s1< s2, the pgfortran compiler would do two things:
                !            1: copy the contents with the size of s1 from to array1 to array2
                !            2: the reminded contents with the sizeof (s2 -s1 +1) in array2
                !               would be assigned to some value automaticlly(the reminded
                !               contents would be filled based on the value in array1, for
                !               instant, if the value in array1 is 1,2,3...., the remineded
                !               contents would be assigned to s1+1,s1+2,....,s2)
                !       In fact, we hope that the value in array2 would be array(1:s1),0(s1+1:s2) ,which means the remined value
                !       shoud not be assigned by some values, thus we should do the additional operation: array2(s1+1:s2) = 0

                ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
                !       the tempArray need not need to be allocated, the compiler would
                !       allocate the array automatic that is marked as "allocatable" based on the
                !       assigned array's size
                call AllocateTwoDimd(tempArray,NewX,oldShape(2),"tempArray")

                !---restore
                tempArray = reshape(SOURCE=[TheArray(NewX,oldShape(2))],SHAPE=[NewX,oldShape(2)])

                call AllocateTwoDimd(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Double],ORDER=[1,2])

            end if

        else

            if(NewY .LE. oldShape(2)) then

                call AllocateTwoDimd(tempArray,NewY,oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[NewY,oldShape(1)],PAD=[zero_Double],ORDER=[2,1])

                call AllocateTwoDimd(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewX,NewY],PAD=[zero_Double],ORDER=[2,1])

            else

                call AllocateTwoDimd(tempArray,oldShape(2),oldShape(1),"tempArray")

                !---restore and reverse
                tempArray = reshape(SOURCE=[TheArray],SHAPE=[oldShape(2),oldShape(1)],PAD=[zero_Double],ORDER=[2,1])

                call AllocateTwoDimd(TheArray,NewX,NewY,"TheArray")

                TheArray = reshape(SOURCE=[reshape(SOURCE=[tempArray],&
                                                   SHAPE=[NewX,oldShape(2)],&
                                                   PAD=[zero_Double],&
                                                   ORDER=[2,1])],&
                                   SHAPE=[NewX,NewY],&
                                   PAD=[zero_Double],&
                                   ORDER=[1,2])

            end if

        end if

        call DeAllocateTwoDimd(tempArray,"tempArray")

    end if

    return
  end subroutine ResizeArrayd_TwoDim


    !*******************************************
    subroutine DuplicateArrayi_OneDim(TheArray,DuplicateNum)
        implicit none
        !---Dummy Vars---
        integer, ARRAYDEV dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DuplicateNum
        !---Local Vars---
        integer::OldSize(1)
        integer, ARRAYDEV dimension(:), allocatable::tempArray
        integer::NewSize
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize = DuplicateNum*OldSize(1)

        if(DuplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateOneDimi(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateOneDimi(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[tempArray])

            call DeAllocateOneDimi(tempArray,"tempArray")

        end if

        return

    end subroutine DuplicateArrayi_OneDim

    !*******************************************
    subroutine DuplicateArrayr_OneDim(TheArray,DuplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KINDSF), ARRAYDEV dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DuplicateNum
        !---Local Vars---
        integer::OldSize(1)
        real(kind=KINDSF), ARRAYDEV dimension(:), allocatable::tempArray
        integer::NewSize
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize = DuplicateNum*OldSize(1)

        if(DuplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateOneDimr(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateOneDimr(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[tempArray])

            call DeAllocateOneDimr(tempArray,"tempArray")

        end if

        return

    end subroutine DuplicateArrayr_OneDim

    !*******************************************
    subroutine DuplicateArrayd_OneDim(TheArray,DuplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KINDDF), ARRAYDEV dimension(:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DuplicateNum
        !---Local Vars---
        integer::OldSize(1)
        real(kind=KINDDF), ARRAYDEV dimension(:), allocatable::tempArray
        integer::NewSize
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize = DuplicateNum*OldSize(1)

        if(DuplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateOneDimd(tempArray,OldSize(1),"tempArray")

            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(1)])

            call AllocateOneDimd(TheArray,NewSize,"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize],PAD=[tempArray])

            call DeAllocateOneDimd(tempArray,"tempArray")

        end if

        return

    end subroutine DuplicateArrayd_OneDim

    !*******************************************
    subroutine DuplicateArrayi_TwoDim(TheArray,DuplicateNum)
        implicit none
        !---Dummy Vars---
        integer, ARRAYDEV dimension(:,:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DuplicateNum
        !---Local Vars---
        integer::OldSize(2)
        integer, ARRAYDEV dimension(:,:), allocatable::tempArray
        integer::NewSize(2)
        integer::istat
        !---Body----

        OldSize = shape(TheArray)

        NewSize(1) = DuplicateNum*OldSize(1)
        NewSize(2) = OldSize(2)

        if(DuplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateTwoDimi(tempArray,OldSize(2),OldSize(1),"tempArray")

            !---Restore and reverse
            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(2),OldSize(1)],PAD=[0],ORDER=[2,1])

            call AllocateTwoDimi(TheArray,NewSize(1),NewSize(2),"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize(1),NewSize(2)],PAD=[tempArray],ORDER=[2,1])

            call DeAllocateTwoDimi(tempArray,"tempArray")

        end if

        return

    end subroutine DuplicateArrayi_TwoDim

    !*******************************************
    subroutine DuplicateArrayr_TwoDim(TheArray,DuplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KINDSF), ARRAYDEV dimension(:,:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DuplicateNum
        !---Local Vars---
        integer::OldSize(2)
        real(kind=KINDSF), ARRAYDEV dimension(:,:), allocatable::tempArray
        real(kind=KINDSF)::zero_Single
        integer::NewSize(2)
        integer::istat
        !---Body----

        zero_Single = 0.D0

        OldSize = shape(TheArray)

        NewSize(1) = DuplicateNum*OldSize(1)
        NewSize(2) = OldSize(2)

        if(DuplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateTwoDimr(tempArray,OldSize(2),OldSize(1),"tempArray")

            !---Restore and reverse
            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(2),OldSize(1)],PAD=[zero_Single],ORDER=[2,1])

            call AllocateTwoDimr(TheArray,NewSize(1),NewSize(2),"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize(1),NewSize(2)],PAD=[tempArray],ORDER=[2,1])

            call DeAllocateTwoDimr(tempArray,"tempArray")

        end if

        return

    end subroutine DuplicateArrayr_TwoDim


    !*******************************************
    subroutine DuplicateArrayd_TwoDim(TheArray,DuplicateNum)
        implicit none
        !---Dummy Vars---
        real(kind=KINDDF), ARRAYDEV dimension(:,:), allocatable,intent(inout)::TheArray
        integer,intent(in)::DuplicateNum
        !---Local Vars---
        integer::OldSize(2)
        real(kind=KINDDF), ARRAYDEV dimension(:,:), allocatable::tempArray
        real(kind=KINDDF)::zero_Double
        integer::NewSize(2)
        integer::istat
        !---Body----

        zero_Double = 0.D0

        OldSize = shape(TheArray)

        NewSize(1) = DuplicateNum*OldSize(1)
        NewSize(2) = OldSize(2)

        if(DuplicateNum .GT. 1) then
            ! note: in fortran2003, while the compiler option -Mallocatable=03 is used
            !       the tempArray need not need to be allocated, the compiler would
            !       allocate the array automatic that is marked as "allocatable" based on the
            !       assigned array's size
            call AllocateTwoDimd(tempArray,OldSize(2),OldSize(1),"tempArray")

            !---Restore and reverse
            tempArray = reshape(SOURCE=[TheArray],SHAPE=[OldSize(2),OldSize(1)],PAD=[zero_Double],ORDER=[2,1])

            call AllocateTwoDimd(TheArray,NewSize(1),NewSize(2),"TheArray")

            TheArray = reshape(SOURCE=[tempArray],SHAPE=[NewSize(1),NewSize(2)],PAD=[tempArray],ORDER=[2,1])

            call DeAllocateTwoDimd(tempArray,"tempArray")

        end if

        return

    end subroutine DuplicateArrayd_TwoDim
