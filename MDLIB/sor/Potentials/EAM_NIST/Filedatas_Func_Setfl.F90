module Filedatas_Func_Setfl
!*** DESCRIPTION:
!    This module is used for ".setfl" file to Spline tables.
!
!*** AUTHOR:
!    Mingjie Qiu, Oct., 2016
!    Thank you for the guidance of my teacher Qing Hou!
!***
!
use MATH90A_MODULE
use MD_TYPEDEF_ForceTable
use My_MiniUtilities_QMJ
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

  private::get_filedatas_setfl
  private::Release_data_parameter, Release_element

  private::GET_DBINT4_DBVALU_PARAMETER
  private::Vr_Spline, Rhor_Spline, Frho_Spline
  private::NN_Spline, RHO_Spline,  EMBED_Spline

contains
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
subroutine get_filedatas_setfl (InputFile,  FTable, elements)
!***PURPOSE: To get file datas form Lammps file
!***INPUT:   InputFile
!***OUTPUT:  m_element_datas
!***
implicit none
!--- dummy variable
  character*(*), intent(in) :: InputFile
  type(MDForceTable),intent(inout)::FTable
  character*(*), allocatable, dimension(:), intent(inout) :: elements
!--- local variables
  integer :: hFile, Nrho, Nr, NE, N, I, J, K
  real(KINDDF) :: drho, dr, cutoff, minR
  character*256 :: str
  character*256 :: comment(3)
  character*32  :: NUMB(10), ISTR, JSTR
  real(KINDDF), allocatable, dimension(:) :: R, RHO   !*** used in DBINT4
!*****
  call AvailableIOUnit(hFile)
  open (UNIT=hFile, FILE=InputFile, STATUS='old')
    read (hFile,fmt='(A)') comment(1)                  !*** read comment in the first three lines
    read (hFile,fmt='(A)') comment(2)
    read (hFile,fmt='(A)') comment(3)
    read (hFile,fmt='(A)') str                         !*** read elements
    !call Extract_Numb(str, 10, N, NUMB)
    !NE = DRSTR( NUMB(1))                               !*** the number of element
    !m_NE = NE                                          !*** m_NE is used in other subroutine of this moudle
    call EXSTR (str, 10, N, NUMB)
    NE = DRSTR( NUMB(1))
    m_NE = NE
    allocate ( m_element_datas(NE), m_FRHO_ZERO( NE), elements(NE) )                   !*** allocate m_element_datas
    m_FRHO_ZERO = 0
    do I = 1, NE
      m_element_datas(I)%element = NUMB(I+1)
      elements(I) = NUMB(I+1)
    end do

    read (hFile,*) Nrho, drho, Nr, dr, cutoff

    FTable%RHOMX = Nrho*drho
    FTable%RMAX = cutoff*CP_A2CM

    allocate ( RHO(Nrho), R(Nr) )
    do I = 1, Nrho
      RHO(I) = (I-1) * drho
    end do
    minR = cutoff/Nr
    do I = 1, Nr
      R(I) = (I-1) * minR
    end do
    !*** read F_rho and Rho_r
    do I = 1, NE
      read (hFile,fmt='(A)') str
      call Extract_Numb(str, 10, N, NUMB)
      m_element_datas(I)%NA = DRSTR (NUMB(1))
      m_element_datas(I)%mass = DRSTR (NUMB(2))
      m_element_datas(I)%constant = DRSTR (NUMB(3))
 !    m_element_datas(I)%crystal_type = NSTR(4)

      allocate (m_element_datas(I)%F_rho%RI(Nrho), m_element_datas(I)%Rho_r%RI(Nr) ,&
                m_element_datas(I)%F_rho%FI(Nrho), m_element_datas(I)%Rho_r%FI(Nr))
      m_element_datas(I)%F_rho%NL = Nrho
      m_element_datas(I)%Rho_r%NL = Nr
      m_element_datas(I)%F_rho%RI(1:Nrho) = RHO(1:Nrho)
      m_element_datas(I)%Rho_r%RI(1:Nr)   = R(1:Nr)
      m_element_datas(I)%F_rho%RIMIN = RHO(1)
      m_element_datas(I)%F_rho%RIMAX = RHO(Nrho)
      m_element_datas(I)%Rho_r%RIMIN = R(1)
      m_element_datas(I)%Rho_r%RIMAX = R(Nr)
      !*** read F_rho
      read (hFile,*) m_element_datas(I)%F_rho%FI(1:Nrho)
      if ( maxval(m_element_datas(I)%F_rho%FI).eq.0 .and. minval(m_element_datas(I)%F_rho%FI).eq.0 )   m_FRHO_ZERO(I) = 1
      !***
      call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%F_rho )
      !*** read Rho_r
      read (hFile,*) str
      str = adjustl(str)
      str = str(1:len_trim(str))
      backspace (hFile)
      if (str .EQ."INF") then
        write(ISTR,*) I
        write(*,*) "r=0,","rho"//ISTR//"=INF"
        read (hFile,fmt='(A)') str
        call Extract_Numb(str, 10, N, NUMB)
        do K = 1, N
          m_element_datas(I)%Rho_r%FI(K+1) = DRSTR (NUMB(K))
        end do
        m_element_datas(I)%Rho_r%FI(1) = m_element_datas(I)%Rho_r%FI(2)
        read (hFile,*) m_element_datas(I)%Rho_r%FI(N+2:Nr)
        !*** get parameters of DBINT4 and DBVALU
        call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%Rho_r )
      else
        read (hFile,*) m_element_datas(I)%Rho_r%FI(1:Nr)
        !*** get parameters of DBINT4 and DBVALU
        call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%Rho_r )
      end if
    end do
      !*** read V_r data
      do I = 1, NE                                                                   !*** read V_r
        allocate(m_element_datas(I)%V_r(I))
                                                                                !*** (I,J) determine V_r data of (I,J) elements
        do J = 1, I
          allocate(m_element_datas(I)%V_r(J)%RI(Nr),m_element_datas(I)%V_r(J)%FI(Nr))
          m_element_datas(I)%V_r(J)%NL = Nr
          m_element_datas(I)%V_r(J)%RI(1:Nr) = R(1:Nr)
          m_element_datas(I)%V_r(J)%RIMIN = R(1)
          m_element_datas(I)%V_r(J)%RIMAX = R(Nr)

          read (hFile,*) str
          str = adjustl(str)
          str = str(1:len_trim(str))
          backspace (hFile)
          if (str .EQ."NAN") then
            write(ISTR,*) I
            write(JSTR,*) J
            write(*,*) "r=0,","V"//"("//ISTR//JSTR//")"//"=NAN"
            read (hFile,fmt='(A)') str
            call Extract_Numb(str, 10, N, NUMB)
            do K = 1, N
              m_element_datas(I)%V_r(J)%FI(K+1) = DRSTR (NUMB(K))
            end do
            m_element_datas(I)%V_r(J)%FI(1) = m_element_datas(I)%V_r(J)%FI(2)
            read (hFile,*) m_element_datas(I)%V_r(J)%FI(N+2:Nr)
            !*** get parameters of DBINT4 and DBVALU
            call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%V_r(J) )
          else
            read (hFile,*) m_element_datas(I)%V_r(J)%FI(1:Nr)
            !*** get parameters of DBINT4 and DBVALU
            call GET_DBINT4_DBVALU_PARAMETER( m_element_datas(I)%V_r(J) )
          end if
        end do
      end do
     close(hFile)
   return
end subroutine get_filedatas_setfl
!************************************************************************************
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
!****************************************************************************
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
!**********************************************************************************************
!***********************************************************************************************************
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

    !**************
    I = (m_ITAB-1)/m_Ne + 1
    J = m_ITAB-(I-1)*m_Ne


     if (.not.m_FRHO_ZERO(I)) then
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
      else
        F  = 0.D0
        dF = 0.D0
     end if
   return
 end subroutine Rhor_Spline
!**********************************************************************************************************
!***********************************************************************************************************
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

    !**************
    I = (m_ITAB-1)/m_Ne + 1
    J = m_ITAB-(I-1)*m_Ne

      if (I==J) then
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
        else
          F  = 0.D0
          dF = 0.D0
        end if
    return
 end subroutine Frho_Spline
!***********************************************************************************************************
!***********************************************************************************************************
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

    POTR = POTR/(r*CP_CM2A)
    FPOTR = ( FPOTR - POTR )/(r*CP_CM2A)
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

    FPOTB = -1.D0*FPOTB * CP_CM2A

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

    FRHO = FRHO * CP_EVERG
    DFRHO = DFRHO * CP_EVERG

    return
end subroutine EMBED_Spline
!****************************************************************************
!****************************************************************************
subroutine Register_ForceTableProc_Setfl(Fname, FTable)
  !***  PURPOSE:   to regieter forcetable process
  !
  !     INPUT:     Fname: the filename for Setfl
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


             call get_filedatas_setfl(Fname, FTable,  elements )
             IT = 0
             do I=1,  m_NE
                do J=1, m_NE
                   IT = IT + 1

                     !if(J.le.I) then
                     !    NOTE = "setfl table: V(r) "//m_element_datas(I)%element(1:len_trim(m_element_datas(I)%element)) // &
                     !            "-"//m_element_datas(J)%element(1:len_trim(m_element_datas(J)%element))
                     !else
                     !  NOTE = "setfl table: V(r) "//m_element_datas(J)%element(1:len_trim(m_element_datas(J)%element)) // &
                     !            "-"//m_element_datas(I)%element(1:len_trim(m_element_datas(I)%element))
                     !end if
                   NOTE= "setfl table: "//elements(I)(1:len_trim(elements(I)))//"<-"  &
                                        //elements(J)(1:len_trim(elements(J)))//" ( " &
                                        //Fname(1:len_trim(Fname))//" )"

                   call Register_ForceTableProc(IT, FTable, NNFORCE=NN_Spline, NEFORCE=RHO_Spline,       &
                                                            EMBDF = EMBED_Spline, NOTIFY=Notify_func,    &
                                                            NOTE= NOTE)


                end do
             end do
        return
end subroutine Register_ForceTableProc_Setfl
!****************************************************************************
end module Filedatas_Func_Setfl
