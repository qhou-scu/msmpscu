  module BUBBLE_BREAKTIME_14
  !***  DESCRIPTION:
  !     This program is used to calculate the time for a bubble breaking
  !
  !    Authored by HOU Qing
  !
  !    DEPENDENCY:
  !
  !    HOSTORY:  2014-09
  !
  !    SEE ALSO:      MD_SimBoxArray_ToolShell_14_GPU.F90
  use MD_Globle_Variables
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use EMBEDMENT_TypeDef_CtrlParam_2010
  implicit none

  !_______________________________________________________________________________________
  !
      !---
      integer::m_NDATA = 0
      !---
      integer,private::m_BubbleID = 2
      real(KINDDF),dimension(:), allocatable, private::m_NA            !the intial number of atoms in a bubble
      real(KINDDF),dimension(:), allocatable, private::m_TIME1         !the time point for a bubble has one atom  release
      real(KINDDF),dimension(:), allocatable, private::m_TIME10        !the time point for a bubble has 10% release
      real(KINDDF),dimension(:), allocatable, private::m_TIME90        !the time point for a bubble has 90% release

  !_______________________________________________________________________________________
       character(len=11),parameter, private::mp_FTAGI="&AUXF_EMBED"

  contains

  !****************************************************************************************
  subroutine Initialize(SimBox,CtrlParam)
  !***  DESCRIPTION: to allocate mempry for diffusion calculations
  !
  !
  implicit none
      type(SimMDBox), intent(in)::SimBox
      type(SimMDCtrl),intent(in)::CtrlParam
      !--- Local variables
      integer::nbox


      !*** allocate meomery storing the displacements
           nbox = CtrlParam%TOTALBOX

             allocate(m_NA(nbox),                    &
                      m_TIME1(nbox),                 &
                      m_TIME10(nbox),                &
                      m_TIME90(nbox))

             m_NA     = 0
             m_TIME1  = 1.D64
             m_TIME10 = 1.D64
             m_TIME90 = 1.D64
             m_NDATA  = 0
             return

  end subroutine Initialize
  !****************************************************************************************

  !****************************************************************************************
  subroutine Clear(SimBox,CtrlParam)
  !***  DESCRIPTION: to allocate mempry for diffusion calculations
  !
  !
  !***  The modules included ******************************************
  implicit none
      type(SimMDBox) ::SimBox
      type(SimMDCtrl)::CtrlParam
      !--- Local variables

             if(allocated(m_NA) )    deallocate(m_NA)
             if(allocated(m_TIME1))  deallocate(m_TIME1)
             if(allocated(m_TIME10)) deallocate(m_TIME10)
             if(allocated(m_TIME90)) deallocate(m_TIME90)

      RETURN
  END subroutine Clear
  !****************************************************************************************

  !****************************************************************************************
  subroutine Record_(JBOX,ITIME, TIME, SimBox, CtrlParam )
  !***  PORPOSE: to get the migration distance  of clusters
  !     INPUT:  SimBox, the  samples
  !     OUTPUT: RAV, RMI, RMA
  implicit none
   !--- input
   type(SimMDBox), dimension(:), intent(in)::SimBox
   type(SimMDCtrl),              intent(in)::CtrlParam
   integer::JBOX,ITIME
   real(KINDDF)::TIME

   !
   integer::NP
   !--- Local variables
   integer::I0, I,J, K, NS, NR

  !----
              NS = size(SimBox)
              if(ITIME .LE. 0) then
                 !--- record the initial number of atoms in a bubble
                 do I= 1, NS
                    m_NA((JBOX-1)*NS+I) = SimBox(I)%NA(m_BubbleID)
                 end do
                 m_NDATA = m_NDATA+NS
              end if

              I0 = (JBOX-1)*NS
              do I=1, NS
                 NR = 0
                 do J=1, SimBox(I)%NPRT
                    if(SimBox(I)%ITYP(J) .eq. m_BubbleID) then
                       if(IAND(SimBox(I)%STATU(J), CP_STATU_ACTIVE) .ne. CP_STATU_ACTIVE) then
                          NR = NR + 1
                        end if
                    end if
                 end do

                 if(m_TIME1(I0+I) .ge. 1.D63) then
                    if(NR .gt. 0) then
                       m_TIME1(I0+I) = TIME
                     end if
                 end if

                 if(m_TIME10(I0+I) .ge. 1.D63) then
                    if(dble(NR)/dble(m_NA(I0+I)) .ge. 0.1D0) then
                       m_TIME10(I0+I) = TIME
                     end if
                 end if

                 if(m_TIME90(I0+I) .ge. 1.D63) then
                    if(dble(NR)/dble(m_NA(I0+I)) .ge. 0.9D0) then
                       m_TIME90(I0+I) = TIME
                     end if
                 end if
             end do

    return
  end subroutine Record_
  !****************************************************************************************

  !****************************************************************************************
  subroutine AfterRecord_(SimBox, CtrlParam )
  !***  PORPOSE: to print out diffusion parameters
  !     INPUT:  SimBox, the  samples
  implicit none
   !--- input
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

   !--- Local variables
   character*16::str(10)=""
   character*64::fmt
   character*32::strnum
   character*256::Outfile
   integer::hFile,ERR
   logical::opened, EX
   real(KINDDF)::TMAX, DT
   INTEGER::I,IS, NTYP, NBINS=10

         !--- determine the output filename
           OutFile = CtrlParam%f_geometry(1:len_trim(CtrlParam%f_geometry))//".diff"

          !--- check the status of the file
          if(m_NDATA .ne. CtrlParam%TOTALBOX) then
             print *, "MDPSCU Error: the number of accumulated data is not equal to the number of boxes"
          end if

          TMAX = 0.D0
          do I=1, m_NDATA
             if(m_TIME90(I) .lt. 1.D60) then
                if(m_TIME90(I) .gt. TMAX) TMAX = m_TIME90(I)
             end if
          end do

          DT = TMAX /dble(NBINS)

          !*** to release allocated memory
          call Clear(SimBox, CtrlParam )
    return
  end subroutine AfterRecord_
  !****************************************************************************************
  end module BUBBLE_BREAKTIME_14

