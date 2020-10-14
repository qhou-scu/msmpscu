 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to seperate multibox configure into a set of independent boxes that
 !                  could be more convenient for visualization.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !                       BoxSelector.F90
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
 !                  MultBox2Boxes.exe SETUP-file
 !                  or:
 !                  MultBox2Boxes.exe SETUP-file -T(est) t1, t2, ts -C(fg) c1, c2, c3 -B(ox) b1, b2, b3
 !
 !                  where:  SETUP-file is the name of setup file used by MD simulations in MDPSCU.
 !                  OPTIONS:
 !                        -T(est) t1, t2, t3 : indenetify the test IDs to be involved
 !                        -C(fg)  c1, c2, c3 : indenetify the configure IDs to be involved
 !                        -B(ox)  b1, b2, b3 : indenetify the boxes IDs in a TEST
 !
 !
 !                  If you needed additional control on the output, you should add a line in
 !                  setup fileL
 !
 !                        &AUXF_SELIN filename
 !
 !                  For example, in the AUXF_SELIN file, you can give the filename of a reference configuration.
 !                  The reference configuration will be merged into the boxes. The 'atoms' in reference configuration
 !                  will be assigned a type identification NGROUP+1.
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2015-04 (Hou Qing Sichuan university)
 !
 !

 module MultBox2Boxes
 use BoxSelector
 implicit none

  contains
  !**********************************************************************************
   subroutine Initialize_(SimBox, CtrlParam)
   !***  PURPOSE:  to Initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use BoxSelector
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables

         call Initialize_BoxSelector(SimBox, CtrlParam)
         call SetSelectProc(MySelectProc)
         return
   end subroutine Initialize_
  !**********************************************************************************

  !**********************************************************************************
   subroutine MySelectProc(SimBox, CtrlParam, FLAG)
  !***  DESCRIPTION: to select which box staisfying given condition
  !
  !    INPUT:  SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !    OUTPUT: FLAG,      =0, or >0, maker of selected boxs
  !
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl

   implicit none
   !----   DUMMY Variables
      type(SimMDBox), dimension(:)::SimBox
      type(SimMDCtrl)             ::CtrlParam
      integer,        dimension(:)::FLAG
   !---- Local variables
      integer::I, IB0, IB1, IBS

            IB0 = CtrlParam%STARTBOX
            if(CtrlParam%ENDBOX.le.0) then
               IB1 = size(SimBox)
            else
               IB1 = CtrlParam%ENDBOX
            end if
            IBS = max(1,CtrlParam%BOXSTEP)

            FLAG = 0
            do I=IB0, IB1, IBS
               FLAG(I) = 1
            end do
            return
  end subroutine MySelectProc
  !****************************************************************************************

  !****************************************************************************************
  subroutine Record_MySelectool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables

            if(Stamp%ITime .lt. 0) return

            call Record_BoxSelector(Stamp, SimBox, CtrlParam)
            return
  end subroutine Record_MySelectool
  !****************************************************************************************

  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use BoxSelector, only:Clear_BoxSelector

   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_BoxSelector(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module MultBox2Boxes

 !****************************************************************************************
 Program MultBoxSeperator_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use MultBox2Boxes

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_,          &
                                RECORDPROC=Record_MySelectool,  &
                                AFTRECORD=MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program MultBoxSeperator_Main
