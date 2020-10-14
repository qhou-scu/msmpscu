module Filedatas_Func_Lspt
!*** DESCRIPTION:
!    This module is used for ".spt" file to Spline tables.
!
!*** AUTHOR:
!    Mingjie Qiu, Oct., 2016
!    Thank you for the guidance of my teacher Qing Hou!
!***
!
use MATH90A_MODULE
use MD_TYPEDEF_ForceTable

implicit none
  private :: filedata_DBparameter
  type :: filedata_DBparameter
  !*** read file datas
    character*256 :: filename             !*** This type is defined to read file data and get the corresponding parameters of DBINT4 and DBVALU.
    real(KINDDF), allocatable, dimension(:) :: RI , FI
    real(KINDDF) :: RIMAX = 0.D0
    real(KINDDF) :: RIMIN = 0.D0
    real(KINDDF) :: FIMAX = 0.D0
    real(KINDDF) :: FIMIN = 0.D0
    integer :: NL                         !*** the number of data lines
  !*** parameters of DBINT4 and DBVALU
    integer :: IBCL     = 2
    integer :: IBCR     = 2
    integer :: KNTOPT   = 1
    integer :: IDERIV0  = 0
    integer :: IDERIV1  = 1
    integer :: INBV  = 1                               !*** INBV must be set to 1 the first time DBVALU is called .
    integer :: K
    integer :: N
    real(KINDDF)  :: FBCL = 0.D0
    real(KINDDF)  :: FBCR = 0.D0
    real(KINDDF), allocatable, dimension(:)  ::  BCOEF
    real(KINDDF), allocatable, dimension(:)  ::  T
    real(KINDDF), allocatable, dimension(:)  ::  W1
    real(KINDDF), allocatable, dimension(:,:)::  W2
  end type filedata_DBparameter

  private :: element_data                             !*** This type is defined to read each element information.
  type :: element_data
    character*10 :: element
    integer :: NA         !*** atomic number
    real(KINDDF) :: mass , constant
    character*10 :: crystal_type

    type(filedata_DBparameter) :: F_rho
    type(filedata_DBparameter) :: Rho_r
    type(filedata_DBparameter),allocatable,dimension(:) :: V_r      !*** one element has N kinds of V_r with all the N elements, so an array is defined.

  end type element_data

  integer, private :: m_NE
  integer, private :: m_ITAB
  integer, allocatable, dimension(:), private :: m_FRHO_ZERO         !***Flag, "1" means FRHO = 0
  type(element_data), allocatable, dimension(:),private :: m_element_datas

  private::get_filedatas_lspt
  private::Release_data_parameter, Release_element

  private::GET_DBINT4_DBVALU_PARAMETER
  private::Vr_Spline, Rhor_Spline, Frho_Spline
  private::NN_Spline, RHO_Spline,  EMBED_Spline
contains
!***************************************************************************************
!***********************************************************************************************************
subroutine Release ()
!***PURPOSE:
!***
implicit none
  integer :: I

  if (allocated(m_element_datas)) then
    do I = 1, size(m_element_datas)
      call Release_element(m_element_datas(I))
    end do
    deallocate (m_element_datas)
  end if
  if (allocated(m_FRHO_ZERO))  deallocate (m_FRHO_ZERO)
  return
end subroutine Release
!***********************************************************************************************************
!***********************************************************************************************************
subroutine Release_element(the_element)
!***
implicit none
  integer :: I

  type (element_data) :: the_element

    call Release_data_parameter (the_element%F_rho)
    call Release_data_parameter (the_element%RHO_r)
    do I = 1, size(the_element%V_r)
      call Release_data_parameter (the_element%V_r(I))
    end do
    deallocate (the_element%V_r)
    return
end subroutine Release_element
!***********************************************************************************************************
!***********************************************************************************************************
subroutine Release_data_parameter(the_data_parameter)
!***
implicit none
  type (filedata_DBparameter) :: the_data_parameter

    if (allocated(the_data_parameter%RI))      deallocate (the_data_parameter%RI)
    if (allocated(the_data_parameter%FI))      deallocate (the_data_parameter%FI)
    if (allocated(the_data_parameter%BCOEF))   deallocate (the_data_parameter%BCOEF)
    if (allocated(the_data_parameter%T))       deallocate (the_data_parameter%T)
    if (allocated(the_data_parameter%W1))      deallocate (the_data_parameter%W1)
    if (allocated(the_data_parameter%W2))      deallocate (the_data_parameter%W2)

  return
end subroutine Release_data_parameter
!***********************************************************************************************************
!***********************************************************************************************************
subroutine get_filedatas_lspt (InputFile, FTable, elements)
!*** PURPOSE: This subroutine is called in "Select_Creat_Table" to read all the ".spt" files and the DBINT4 will be called
!             to get all the paremeters. All the datas will be stored in type m_filedatas.
!
implicit none
  !--- dummy variable
    character*(*), intent(in)  :: InputFile
    type(MDForceTable),  intent(inout)::FTable
    character*(*), allocatable, dimension(:), intent(inout) :: elements
  !--- local variable
    integer       ::  hFile,hFile_Frho,hFile_RHOr,hFile_Vr , I , J, K, N, status, LINE,  NL
    integer       ::  NE            !*** the kinds of elements
    character(len=5), parameter :: STARTFLAG = "&LSPT"
    character*256 :: STR, string, filename, PATH
    character*32::STRNUMB(10), KEYWORD,SUBSTR(10), INUM
    real(KINDDF), allocatable, dimension(:) :: R ,RHO
    real(KINDDF) :: RMIN, RMAX, RHOMIN, RHOMAX

  !*** get filenames
    RHOMAX = 0.D0                      !*** set the initial value of RHOMAX and RMAX
    RMAX = 0.D0

    call GetPath(InputFile, PATH)
    call AvailableIOUnit(hFile)
    open (UNIT=hFile,FILE= InputFile, STATUS='old')
      call GetInputStrLine(hFile,STR, LINE, "!", *100)
      STR = adjustl(STR)
      call GetKeyWord("&", STR, KEYWORD)
      call UpCase(KEYWORD)
      if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. STARTFLAG) then         !*** get the start flag "&LSPT"
    !**** to start lspt file datas ****************
        call GetInputStrLine(hFile,STR, LINE, "!", *100)
        STR = adjustl(STR)
        call EXTRACT_NUMB(STR,1,N,STRNUMB)                          !*** get the number of elements
        if(N .lt. 1) then
          write(*,fmt="(' MDPSCU Error: the number of groups should be given')")
          write(*,fmt="(' Process to be stopped')")
          stop
        end if
        NE = ISTR(STRNUMB(1))
        m_NE = NE
        allocate (m_element_datas(NE),  m_FRHO_ZERO( NE),elements(NE))                              !*** allocate the number of type elemets_data
        m_FRHO_ZERO = 0
        call GetInputStrLine(hFile,STR, LINE, "!", *100)
        STR = adjustl(STR)
        call Extract_Substr(STR,10,N,SUBSTR)
        if(N /= NE) then
          write(*,fmt="(' MDPSCU Warning: the number of elements is different from the given number, please check again.')")
          write(*,fmt="(' Process to be stopped')")
          stop
        end if
        do I = 1, NE
          m_element_datas(I)%element = adjustl(SUBSTR(I)(1:len_trim(SUBSTR(I))))        !*** get element symbol
          elements(I) = adjustl(SUBSTR(I)(1:len_trim(SUBSTR(I))))
        end do

        do I = 1, NE
          !*** read F_rho data
          call GetInputStrLine(hFile,STR, LINE, "!", *100)
          STR = adjustl(STR)
          call Extract_Substr(STR,1,N,SUBSTR)
          m_element_datas(I)%F_rho%filename = adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1)))) !*** get F_rho file name
          filename = adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1))))
          call UpCase ( filename )
          if  ( filename .EQ. "NA" ) then
            write(INUM,*) I
            !write(*,*) "The file of F_rho("//INUM//") is NA"
            m_FRHO_ZERO(I) = 1
          else
            call AvailableIOUnit(hFile_Frho)                                            !*** read file data
            open( UNIT=hFile_Frho,FILE= PATH(1:len_trim(PATH))//adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1)))), STATUS='old')
              call GetInputStrLine(hFile_Frho,STR, LINE, "#", *100)
              backspace (hFile_Frho)
              NL = 0
              do while(.true.)
                read (hFile_Frho,*,end=1000) str
                NL = NL+1
              end do
         1000 rewind (hFile_Frho)
              m_element_datas(I)%F_rho%NL = NL
              allocate ( m_element_datas(I)%F_rho%RI(NL),m_element_datas(I)%F_rho%FI(NL) )  !**** allocate  m_filedatas
              call GetInputStrLine(hFile_Frho,STR, LINE, "#", *100)
              backspace (hFile_Frho)
              do K = 1, NL
                read(hFile_Frho, *) m_element_datas(I)%F_rho%RI(K),m_element_datas(I)%F_rho%FI(K)
              end do
            close(hFile_Frho)
              m_element_datas(I)%F_rho%RIMAX = maxval ( m_element_datas(I)%F_rho%RI )
              m_element_datas(I)%F_rho%RIMIN = minval ( m_element_datas(I)%F_rho%RI )
              !*** get RHOMAX
              if ( RHOMAX  .lt. m_element_datas(I)%F_rho%RIMAX ) RHOMAX = m_element_datas(I)%F_rho%RIMAX
              !*** get parameters of DBINT4 and DBVALU
              call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%F_rho )
          end if

          !*** read Rho_r data
          call GetInputStrLine(hFile,STR, LINE, "!", *100)
          STR = adjustl(STR)
          call Extract_Substr(STR,1,N,SUBSTR)
          m_element_datas(I)%Rho_r%filename = adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1))))   !*** get Rho_r filename
          filename = adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1))))
          call UpCase ( filename )
          if  ( filename .EQ. "NA" ) then
            write(INUM,*) I
            !write(*,*) "The file of R_rho("//INUM//") is NA"
          else
            call AvailableIOUnit(hFile_RHOr)
             open( UNIT=hFile_RHOr,FILE= PATH(1:len_trim(PATH))//adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1)))), STATUS='old')
              call GetInputStrLine(hFile_RHOr,STR, LINE, "#", *100)
              backspace (hFile_RHOr)
              NL = 0
              do while(.true.)
                read (hFile_RHOr,*,end=2000) str
                NL = NL+1
              end do
         2000 rewind (hFile_RHOr)
              m_element_datas(I)%Rho_r%NL = NL
              allocate ( m_element_datas(I)%Rho_r%RI(NL), m_element_datas(I)%Rho_r%FI(NL) )  !**** allocate  m_filedatas
              call GetInputStrLine(hFile_RHOr,STR, LINE, "#", *100)
              backspace(hFile_RHOr)
              do K = 1, NL
                read(hFile_RHOr, *) m_element_datas(I)%Rho_r%RI(K), m_element_datas(I)%Rho_r%FI(K)
              end do
            close(hFile_RHOr)
              m_element_datas(I)%Rho_r%RIMAX = maxval ( m_element_datas(I)%Rho_r%RI )
              m_element_datas(I)%Rho_r%RIMIN = minval ( m_element_datas(I)%Rho_r%RI )
              !*** get parameters of DBINT4 and DBVALU
              call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%Rho_r )
          end if
        end do

        !*** read V_r data
        do I = 1, NE
          call GetInputStrLine(hFile,STR, LINE, "!", *100)
          STR = adjustl(STR)
          call Extract_Substr(STR,10,N,SUBSTR)
          if (N .NE. I)  then
            write(INUM,*) I
            write(*,*) "The VR file number in line"//INUM//"is not right, please check it!"
          else
            allocate ( m_element_datas(I)%V_r(I))
            do J = 1, I
              m_element_datas(I)%V_r(J)%filename = adjustl(SUBSTR(J)(1:len_trim(SUBSTR(J))))
              filename = adjustl(SUBSTR(J)(1:len_trim(SUBSTR(J))))
              call AvailableIOUnit(hFile_Vr)
              open( UNIT=hFile_Vr,FILE= PATH(1:len_trim(PATH))//adjustl(SUBSTR(1)(1:len_trim(SUBSTR(1)))), STATUS='old')
                call GetInputStrLine(hFile_Vr,STR, LINE, "#", *100)
                backspace (hFile_Vr)
                NL = 0
                do while(.true.)
                  read (hFile_Vr,*,end=3000) str
                  NL = NL+1
                end do
           3000 rewind (hFile_Vr)
                m_element_datas(I)%V_r(J)%NL = NL
                allocate ( m_element_datas(I)%V_r(J)%RI(NL),m_element_datas(I)%V_r(J)%FI(NL))
                call GetInputStrLine(hFile_Vr,STR, LINE, "#", *100)
                backspace(hFile_Vr)
                do K = 1, NL
                  read(hFile_Vr, *)  m_element_datas(I)%V_r(J)%RI(K),m_element_datas(I)%V_r(J)%FI(K)
                end do
                m_element_datas(I)%V_r(J)%RIMAX = maxval( m_element_datas(I)%V_r(J)%RI )
                m_element_datas(I)%V_r(J)%RIMIN = minval( m_element_datas(I)%V_r(J)%RI )
              close(hFile_Vr)
                if (RMAX .lt. m_element_datas(I)%V_r(J)%RIMAX ) RMAX = m_element_datas(I)%V_r(J)%RIMAX
                !*** get parameters of DBINT4 and DBVALU
                call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%V_r(J) )
              end do
          end if
        end do
     end if

     FTable%RHOMX = RHOMAX
     FTable%RMAX = RMAX*CP_A2CM

   close(hFile)
   return
 !---------------------------------------------------------------
 100 close(hFile)
     write(*,*)"MDPSCU Error in getting file."
     write(*,*)"The process to be stopped."
     stop

end subroutine get_filedatas_lspt
!*******************************************************************************
!*******************************************************************************
subroutine GET_DBINT4_DBVALU_PARAMETER (the_data_parameter)
!**
implicit none
!--- dummy variable
type(filedata_DBparameter):: the_data_parameter
!--- local variable

  allocate ( the_data_parameter%T     ( the_data_parameter%NL + 6 ), &                                         !*** allocate space used in DBINT4
             the_data_parameter%BCOEF ( the_data_parameter%NL + 2 ), &
             the_data_parameter%W2(5, ( the_data_parameter%NL + 2 )))

  call  DBINT4 ( the_data_parameter%RI,    the_data_parameter%FI, the_data_parameter%NL, &
                 the_data_parameter%IBCL,  the_data_parameter%IBCR, &
                 the_data_parameter%FBCL,  the_data_parameter%FBCR, &
                 the_data_parameter%KNTOPT,the_data_parameter%T, &
                 the_data_parameter%BCOEF, the_data_parameter%N, &
                 the_data_parameter%K,     the_data_parameter%W2 )

 allocate ( the_data_parameter%W1( 3 * the_data_parameter%K ) )           !*** allocate space used in DBVALU

end subroutine GET_DBINT4_DBVALU_PARAMETER
!*******************************************************************************
!****************************************************************************
subroutine Notify_func(ITAB)
  !***  PURPOSE: ger an m_ITAB according IFORCE to choose the corresponding m_filedatas(I,J)
  !     INPUT:   ITAB
  !     OUTPUT:  m_ITAB, used in this entire module
implicit none
  !----dummy variables
    integer,intent(in)::ITAB
  !--- local variables
    m_ITAB = ITAB
    return
 end subroutine Notify_func
!****************************************************************************
!****************************************************************************
subroutine Vr_Spline( r, F, dF)
  !***  PURPOSE: to spline the force table.
  !     INPUT:   r, the distance between two particle
  !
  !     OUTPUT:  F, dF
implicit none
  !----dummy variables
    real(KINDDF),intent(in) :: r
    real(KINDDF),intent(out):: F, dF
  !--- local variables
    integer ::I,  J ,temp

    !**************
    I = (m_ITAB-1)/m_Ne + 1
    J = m_ITAB-(I-1)*m_Ne

      if (J.GT.I) then
        temp = I
        I = J
        J = temp
      end if

      if ( r .ge. m_element_datas(I)%V_r(J)%RIMIN .and. r .le. m_element_datas(I)%V_r(J)%RIMAX ) then

        F =  DBVALU ( m_element_datas(I)%V_r(J)%T, m_element_datas(I)%V_r(J)%BCOEF,&
                      m_element_datas(I)%V_r(J)%N, m_element_datas(I)%V_r(J)%K,&
                      m_element_datas(I)%V_r(J)%IDERIV0, r, &
                      m_element_datas(I)%V_r(J)%INBV, m_element_datas(I)%V_r(J)%W1 )

        dF = DBVALU ( m_element_datas(I)%V_r(J)%T, m_element_datas(I)%V_r(J)%BCOEF,&
                      m_element_datas(I)%V_r(J)%N, m_element_datas(I)%V_r(J)%K,&
                      m_element_datas(I)%V_r(J)%IDERIV1, r, &
                      m_element_datas(I)%V_r(J)%INBV, m_element_datas(I)%V_r(J)%W1 )
      else if ( r.lt. m_element_datas(I)%V_r(J)%RIMIN ) then

        F =  DBVALU ( m_element_datas(I)%V_r(J)%T, m_element_datas(I)%V_r(J)%BCOEF,&
                      m_element_datas(I)%V_r(J)%N, m_element_datas(I)%V_r(J)%K,&
                      m_element_datas(I)%V_r(J)%IDERIV0, m_element_datas(I)%V_r(J)%RIMIN, &
                      m_element_datas(I)%V_r(J)%INBV, m_element_datas(I)%V_r(J)%W1 )
        dF = 0.D0
      else
        F =  0.D0
        dF = 0.D0
      end if
    return
 end subroutine Vr_Spline
!*******************************************************************************************
!*******************************************************************************************
subroutine Rhor_Spline( r, F, dF)
  !***  PURPOSE: to spline the force table.
  !     INPUT:   r, the distance between two particle
  !
  !     OUTPUT:  F, dF
implicit none
  !----dummy variables
    real(KINDDF),intent(in) :: r
    real(KINDDF),intent(out):: F, dF
  !--- local variables
    integer ::I,  J
    character*256 :: filename

    !**************
    I = (m_ITAB-1)/m_Ne + 1
    J = m_ITAB-(I-1)*m_Ne


      if (.not.m_FRHO_ZERO(I)) then
        filename = m_element_datas(J)%Rho_r%filename
        call UpCase (filename)
        if ( filename .EQ. "NA" )  then
           F = 0.D0
          dF = 0.D0
        else
          if ( r .ge. m_element_datas(J)%Rho_r%RIMIN .and. r .le. m_element_datas(J)%Rho_r%RIMAX ) then

            F =  DBVALU ( m_element_datas(J)%Rho_r%T, m_element_datas(J)%Rho_r%BCOEF,&
                          m_element_datas(J)%Rho_r%N, m_element_datas(J)%Rho_r%K,&
                          m_element_datas(J)%Rho_r%IDERIV0, r, &
                          m_element_datas(J)%Rho_r%INBV, m_element_datas(J)%Rho_r%W1 )

            dF = DBVALU ( m_element_datas(J)%Rho_r%T, m_element_datas(J)%Rho_r%BCOEF,&
                          m_element_datas(J)%Rho_r%N, m_element_datas(J)%Rho_r%K,&
                          m_element_datas(J)%Rho_r%IDERIV1, r, &
                          m_element_datas(J)%Rho_r%INBV, m_element_datas(J)%Rho_r%W1 )

          else if (r .lt. m_element_datas(J)%Rho_r%RIMIN ) then

            F =  DBVALU ( m_element_datas(J)%Rho_r%T, m_element_datas(J)%Rho_r%BCOEF,&
                          m_element_datas(J)%Rho_r%N, m_element_datas(J)%Rho_r%K,&
                          m_element_datas(J)%Rho_r%IDERIV0, m_element_datas(J)%Rho_r%RIMIN, &
                          m_element_datas(J)%Rho_r%INBV, m_element_datas(J)%Rho_r%W1 )
            dF = 0.D0

          else
            F =  0.D0
            dF = 0.D0
          end if
        end if
      else
        F  = 0.D0
        dF = 0.D0
     end if
   return
 end subroutine Rhor_Spline
!***************************************************************************************
!***************************************************************************************
subroutine Frho_Spline( rho, F, dF)
  !***  PURPOSE: to spline the force table.
  !     INPUT:   rho, the electron density
  !
  !     OUTPUT:  F, dF
implicit none
  !----dummy variables
    real(KINDDF),intent(in) :: rho
    real(KINDDF),intent(out):: F, dF
  !--- local variables
    integer ::I,  J
    character*256 :: filename
    !**************
    I = (m_ITAB-1)/m_Ne + 1
    J = m_ITAB-(I-1)*m_Ne

      if (I==J) then
        filename = m_element_datas(I)%F_rho%filename
        call UpCase (filename)
        if ( filename .EQ. "NA" ) then
          F = 0.D0
          dF = 0.D0
        else
          if ( rho .ge. m_element_datas(I)%F_rho%RIMIN .and. rho .le. m_element_datas(I)%F_rho%RIMAX ) then

            F =  DBVALU ( m_element_datas(I)%F_rho%T, m_element_datas(I)%F_rho%BCOEF,&
                          m_element_datas(I)%F_rho%N, m_element_datas(I)%F_rho%K,&
                          m_element_datas(I)%F_rho%IDERIV0, rho, &
                          m_element_datas(I)%F_rho%INBV, m_element_datas(I)%F_rho%W1 )

            dF = DBVALU ( m_element_datas(I)%F_rho%T, m_element_datas(I)%F_rho%BCOEF,&
                          m_element_datas(I)%F_rho%N, m_element_datas(I)%F_rho%K,&
                          m_element_datas(I)%F_rho%IDERIV1, rho, &
                          m_element_datas(I)%F_rho%INBV, m_element_datas(I)%F_rho%W1 )

          else if ( rho .lt. m_element_datas(I)%F_rho%RIMIN  ) then
            F = DBVALU ( m_element_datas(I)%F_rho%T, m_element_datas(I)%F_rho%BCOEF,&
                         m_element_datas(I)%F_rho%N, m_element_datas(I)%F_rho%K,&
                         m_element_datas(I)%F_rho%IDERIV0, m_element_datas(I)%F_rho%RIMIN , &
                         m_element_datas(I)%F_rho%INBV, m_element_datas(I)%F_rho%W1 )
            dF = 0.D0
          else if ( rho .gt. m_element_datas(I)%F_rho%RIMAX ) then
            F = DBVALU ( m_element_datas(I)%F_rho%T, m_element_datas(I)%F_rho%BCOEF,&
                         m_element_datas(I)%F_rho%N, m_element_datas(I)%F_rho%K,&
                         m_element_datas(I)%F_rho%IDERIV0, m_element_datas(I)%F_rho%RIMAX , &
                         m_element_datas(I)%F_rho%INBV, m_element_datas(I)%F_rho%W1 )
            dF = 0.D0
          end if
        end if
      else
        F  = 0.D0
        dF = 0.D0
      end if
    return
 end subroutine Frho_Spline
!********************************************************************************************
!********************************************************************************************
subroutine NN_Spline(r,POTR, FPOTR)
  !***  PURPOSE: to spline the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
implicit none
  !----dummy variables
    real(KINDDF),intent(in)::r
    real(KINDDF),intent(out)::POTR, FPOTR
  !--- local variables

    call Vr_Spline( r*CP_CM2A,POTR, FPOTR)

    POTR = 0.5D0*POTR * CP_EVERG
    FPOTR = -1.D0*FPOTR * CP_EVERG*CP_CM2A

    return
 end subroutine NN_Spline
!****************************************************************************
!****************************************************************************
subroutine RHO_Spline(r,POTB, FPOTB)
  !***  PURPOSE: to spline the pair potential RHO(r) and  -dRHO(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTB, FPOTB
implicit none
  !----dummy variables
    real(KINDDF),intent(in)::r
    real(KINDDF),intent(out)::POTB, FPOTB
  !--- local variables

    call Rhor_Spline( r*CP_CM2A,POTB, FPOTB)

    FPOTB = (-1)*FPOTB * CP_CM2A

    return
 end subroutine RHO_Spline
!****************************************************************************
!****************************************************************************
subroutine EMBED_Spline(rho,FRHO, DFRHO)
  !***  PURPOSE: to spline the embedment function F(rho) and  dF(rho)/drho
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
implicit none
  !----dummy variables
    real(KINDDF),intent(in)::rho
    real(KINDDF),intent(out)::FRHO, DFRHO
  !--- local variables

    call Frho_Spline( rho, FRHO, DFRHO)

    FRHO  = FRHO * CP_EVERG
    DFRHO = DFRHO * CP_EVERG

    return
end subroutine EMBED_Spline
!****************************************************************************
!****************************************************************************
subroutine Register_ForceTableProc_SPT(Fname, FTable)
  !***  PURPOSE:   to regieter forcetable process
  !
  !     INPUT:     Fname: the filename for lspt
  !
  !     OUTPUT:    FTable, the table with potential processes registered
  !----------------------------------------------------------------------------------------
  implicit none
      !--- dummy vaiables
      character*(*), intent(in)::       Fname
      type(MDForceTable),intent(inout)::FTable
      !--- local vaiables
      integer     :: I, J, IT
      character*256::NOTE
      character*10, allocatable, dimension(:):: elements

             call get_filedatas_lspt (Fname, FTable,  elements )
             IT = 0
             do I=1,  m_NE
                do J=1, m_NE
                   IT = IT + 1
                   !if(I .eq. J) then
                   !   NOTE = "spt table: V(r)"//m_element_datas(I)%V_r(J)%filename(1:len_trim(m_element_datas(I)%V_r(J)%filename))
                   !   NOTE = NOTE(1:len_trim(NOTE))//", RHO(r)"//m_element_datas(I)%RHO_r%filename(1:len_trim(m_element_datas(I)%RHO_r%filename))
                   !   NOTE = NOTE(1:len_trim(NOTE))//", F(rho)"//m_element_datas(I)%F_rho%filename(1:len_trim(m_element_datas(I)%F_rho%filename))
                   !else
                   !   if(J.le.I) then
                   !      NOTE = "spt table: V(r)"//m_element_datas(I)%V_r(J)%filename(1:len_trim(m_element_datas(I)%V_r(J)%filename))
                   !   else
                   !      NOTE = "spt table: V(r)"//m_element_datas(J)%V_r(I)%filename(1:len_trim(m_element_datas(J)%V_r(I)%filename))
                   !   end if
                   !end if
                    NOTE= "lspt table: "//elements(I)(1:len_trim(elements(I)))//"<-"  &
                                        //elements(J)(1:len_trim(elements(J)))//"  ( " &
                                        //Fname(1:len_trim(Fname))//" )"


                   call Register_ForceTableProc(IT, FTable, NNFORCE=NN_Spline, NEFORCE=RHO_Spline,     &
                                                           EMBDF = EMBED_Spline, NOTIFY=Notify_func,   &
                                                           NOTE= NOTE)

                end do
             end do
        return
end subroutine Register_ForceTableProc_SPT
!****************************************************************************
end module Filedatas_Func_Lspt
