 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to track the atoms with displacement larger than
 !                  a given threshold (in A)
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !                       AtomSelectorCommon.F90
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one lines need to be add:
 !
 !                    &AUXF_SELOUT filename
 !
 !                  in the setup file. The filename is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  With the input file(s) are ready, the program can be run on the command line:
 !
 !                  AtomSelectorDisp.exe SETUP-file
 !                  or:
 !                  AtomSelectorDisp.exe SETUP-file -T(est) t1, t2, ts -C(fg) c1, c2, c3 -B(ox) b1, b2, b3
 !
 !                  where:  SETUP-file is the name of setup file used by MD simulations in MDPSCU.
 !                  OPTIONS:
 !                        -T(est) t1, t2, t3 : indenetify the test IDs to be involved
 !                        -C(fg)  c1, c2, c3 : indenetify the configure IDs to be involved
 !                        -B(ox)  b1, b2, b3 : indenetify the boxes IDs in a TEST
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2015-04 (Hou Qing Sichuan university)
 !
 !

 module AtomTrackDisp
  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MiniUtilities
    implicit none
       real(KINDDF)::m_RCUT2 = 0.D0
  contains
  !**********************************************************************************
   subroutine Initialize_(SimBox, CtrlParam)
   !***  PURPOSE:  to Initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use AtomTrackCommon, only:Initialize_AtomTrackCommon, SetTrackProc, m_INFILE
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables

         call Initialize_AtomTrackCommon(SimBox, CtrlParam)
         call LoadControlParameters(m_INFILE, SimBox, CtrlParam)
         call SetTrackProc(MySelectProc)
         return
   end subroutine Initialize_
  !**********************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, SimBox, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    SimBox,   the simulation box, with the member propColTitle
  !                          could be changed.
  !                CtrlParam, the control parameters, with the member NEEDDAMP
  !                           could be changed.
  !
  !
  !
  implicit none
  !----   DUMMY Variables
   character*(*)  ::fname
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
   integer::hFile, N, LINE
   character*256::STR
   character*32::STRNUMB(10), KEYWORD

            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old', err=200)
            m_RCUT2 = 3.D0
            !*** start loading control parametgers specific for this module
              LINE = 0
              do while(.TRUE.)
                  call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyword("&", STR, KEYWORD)
                  call UpCase(KEYWORD)

                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                        case("&DISPCUT")
                            !$$*** To get box range to be analysis
                            call Extract_Numb(STR,1,n,STRNUMB)
                            m_RCUT2 = DRSTR(STRNUMB(1))
                  end select
              end do
    100     close(hFile)
            m_RCUT2 = m_RCUT2*CP_A2CM
            m_RCUT2 = m_RCUT2*m_RCUT2
            return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in BoxSelector module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop
            return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine MySelectProc(SimBox, CtrlParam, FLAG)
  !***  DESCRIPTION: to select which box staisfying given condition
  !
  !    INPUT:  SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !    OUTPUT: FLAG,      =0, or >0, maker of selected boxs
  !
   implicit none
   !----   DUMMY Variables
      type(SimMDBox), dimension(:)::SimBox
      type(SimMDCtrl)             ::CtrlParam
      integer,        dimension(:)::FLAG
   !---- Local variables
      integer::I, J, IP

            IP = 0
            do I=1, size(SimBox)
               do J=1, SimBox(I)%NPRT
                  IP = IP + 1
                  if(sum(SimBox(I)%DIS(J,1:3)*SimBox(I)%DIS(J,1:3)) .gt. m_RCUT2) then
                     FLAG(IP) = 1
                  else
                     FLAG(IP) = 0
                  end if
               end do
            end do
       return
  end subroutine MySelectProc
  !**********************************************************************************

  !****************************************************************************************
  subroutine Record_MySelectool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use AtomTrackCommon, only:Record_AtomTrackCommon
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables

            if(Stamp%ITime .lt. 0) return
            call Record_AtomTrackCommon(Stamp, SimBox, CtrlParam)
            return
  end subroutine Record_MySelectool
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use AtomTrackCommon, only:Clear_AtomTrackCommon

   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_AtomTrackCommon(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module AtomTrackDisp

 !***********************************************************************************
 Program AtomSelectorDisp_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use AtomTrackDisp

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_,          &
                                RECORDPROC=Record_MySelectool,  &
                                AFTRECORD=MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program AtomSelectorDisp_Main
