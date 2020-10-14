module Filedatas_Func_Moldy
!*** DESCRIPTION:
!    This module is used for ".moldy" file to Spline tables.
!
!*** AUTHOR:
!    Mingjie Qiu, Oct., 2016
!    Thank you for the guidance of my teacher Qing Hou!
!***
!
use MD_TYPEDEF_ForceTable
implicit none
  character*12          ,private :: symb
  real(KINDDF), allocatable,dimension(:) ,private:: m_ak_V , m_rk_V
  real(KINDDF), allocatable,dimension(:) ,private:: m_Ak_RHO , m_Rk_RHO
  real(KINDDF), allocatable,dimension(:) ,private:: m_const                      !*** m_const(1) must be lattice constant

  private::get_filedatas_moldy
  private::EmbedFunc
  private::PairFunc
  private::RhoFunc
contains
!*********************************************************************************
subroutine Release()
!***
implicit none
  if (allocated(m_ak_V))  deallocate (m_ak_V)
  if (allocated(m_rk_V))  deallocate (m_rk_V)
  if (allocated(m_Ak_RHO))deallocate (m_Ak_RHO)
  if (allocated(m_const)) deallocate (m_const)
  return
end subroutine Release
!*********************************************************************************
  subroutine get_filedatas_moldy (IniFile)
  !***  PURPOSE:  This subroutine is called in "Select_Creat_Table" to read all the ".moldy" files
  !               Datas in the file will be stored in the corresponding module variables.
  !
  implicit none
  !--- dummy variable
    character*(*),              intent(in)  :: IniFile
  !--- local variable
    integer       :: I , hFile, N
    character*256::string
    character*32::STRNUMB(30)
  !**************
    call AvailableIOUnit(hFile)
    open (UNIT=hFile, FILE= IniFile, STATUS='old')
      read(hFile,*)  symb
      read(hFile,fmt="(A)") string
        call Extract_Numb(string,30,N,STRNUMB)
        allocate ( m_ak_V(N) )
        do I = 1 , N
          m_ak_V(I) = DRSTR ( STRNUMB(I) )
        end do
      read(hFile,fmt="(A)") string
        call Extract_Numb(string,30,N,STRNUMB)
        allocate ( m_rk_V(N) )
        do I = 1 , N
          m_rk_V(I) = DRSTR ( STRNUMB(I) )
        end do
      read(hFile,fmt="(A)") string
        call Extract_Numb(string,30,N,STRNUMB)
        allocate ( m_Ak_RHO(N) )
        do I = 1 , N
          m_Ak_RHO(I) = DRSTR ( STRNUMB(I) )
        end do
      read(hFile,fmt="(A)") string
        call Extract_Numb(string,30,N,STRNUMB)
        allocate ( m_Rk_RHO(N) )
        do I = 1 , N
          m_Rk_RHO(I) = DRSTR ( STRNUMB(I) )
        end do
      read(hFile,fmt="(A)") string
        call Extract_Numb(string,30,N,STRNUMB)
        allocate ( m_const(N) )
        do I = 1 , N
          m_const(I) = DRSTR ( STRNUMB(I) )
        end do
    close(hFile)
      m_rk_V   = m_rk_V * m_const(1)
      m_Rk_RHO = m_Rk_RHO * m_const(1)
      m_ak_V   = m_ak_V / (m_const(1)**3.D0)
      m_Ak_RHO = m_Ak_RHO / (m_const(1)**6.D0)
  end subroutine get_filedatas_moldy
!*********************************************************************************
!*********************************************************************************
  subroutine PairFunc(  R, POTR, FPOTR )
  !***  PURPOSE: to calculate the pair potential V(r) and  -dV(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT:  POTR, FPOTR
  use MD_Pot_EAM_Utilities, only:NN_FuncPoly3
  implicit none
  !--- dummy variable
    real(KINDDF), intent(in) :: R
    real(KINDDF), intent(out):: POTR, FPOTR
    !******
      call NN_FuncPoly3( R*CP_CM2A, POTR, FPOTR, m_ak_V, m_rk_V )
      return
  end subroutine PairFunc
!*********************************************************************************
!*********************************************************************************
  subroutine RhoFunc(R, POTB, FPOTB )
  !***  PURPOSE: to calculate the electron density RHO(r) and  -dRHO(r)/dr
  !     INPUT:   r, the distance between two particle
  !     OUTPUT: POTB, FPOTB
  use MD_Pot_EAM_Utilities, only:RHO_FuncPoly3
  implicit none
  ! ---- dummy variable
    real(KINDDF), intent(in) :: R
    real(KINDDF), intent(out):: POTB, FPOTB
    !******
      call RHO_FuncPoly3( R*CP_CM2A,POTB, FPOTB, m_Ak_RHO, m_Rk_RHO )
      return
  end subroutine RhoFunc
!*********************************************************************************
!*********************************************************************************
  subroutine EmbedFunc(rho,FRHO, DFRHO)
  !***  PURPOSE: to calculate the embedment function F(rho) and  dF(rho)/drho
  !
  !     INPUT:   rho, the electron density
  !     OUTPUT:  FRHO, DFRHO
  use MD_Pot_EAM_Utilities, only:EMBED_FS_PLOY_Func
  implicit none
      real(KINDDF),intent(in)::rho
      real(KINDDF),intent(out)::FRHO, DFRHO
      real(KINDDF), parameter ::AF(1) = (/-1.d0/)
      !******
        call EMBED_FS_PLOY_Func(rho,FRHO, DFRHO, AF)
        !print*, rho,FRHO
        !pause
       ! FRHO  = - DSQRT(rho)
       ! DFRHO = - (1.D0/2.D0)*(1.D0/DSQRT(rho))
      return
  end subroutine EmbedFunc
!****************************************************************************
!****************************************************************************
subroutine Register_ForceTableProc_Moldy(Fname, FTable)
  !***  PURPOSE:   to regieter forcetable process
  !
  !     INPUT:     Fname: the filename for moldy parameters
  !
  !     OUTPUT:    FTable, the table with potential processes registered
  !----------------------------------------------------------------------------------------
  implicit none
      !--- dummy vaiables
      character*(*), intent(in)::       Fname
      type(MDForceTable),intent(inout)::FTable
      !--- local vaiables
      integer     :: I, NE
             call get_filedatas_moldy(Fname)
             call Register_ForceTableProc(1, FTable, NNFORCE=PairFunc, NEFORCE=RhoFunc, EMBDF =EmbedFunc,  &
                                               NOTE= "Moldy table:"//Fname(1:len_trim(Fname)))
        return
end subroutine Register_ForceTableProc_Moldy
!****************************************************************************
end module Filedatas_Func_Moldy
