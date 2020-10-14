 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program select out the boxs by centro-symmetry
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       BoxSelector.F90.F90
 !                       CentroSymParam_for_BCC_GPU.F90
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
 !                  BoxSelector_CSP.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2014-11 (Hou Qing Sichuan university)
 !                  modified      2014-11 (Wang Xiaoshuan Sichuan university)
 !
 !

 module BoxSelector_CSP
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

  !****************************************************************************************
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

   !--- use which module depend on your selection rule
   use MD_Globle_Variables_GPU, only:dm_NPRT
   use CentroSymParam_for_BCC_GPU, only:Cal_CentroSymmetry_DEV
   implicit none
   !----   DUMMY Variables
      type(SimMDBox), dimension(:)            ::SimBox
      type(SimMDCtrl)                         ::CtrlParam
      integer(1),    dimension(:), allocatable::FLAG
   !---- Local variables
      real(KINDSF), dimension(:), allocatable::CP
      integer::I, IA

            allocate(CP(dm_NPRT))
            call Cal_CentroSymmetry_DEV(SimBox(1), CtrlParam, CP)

            do I=1, size(SimBox)
               IA = (I-1)*SimBox(1)%NPRT + SimBox(1)%NPRT
               if(CP(IA) .lt. 0.1) then
                  FLAG(I) = 1
               end if
            end do
            deallocate(CP)

            return
  end subroutine MySelectProc
  !****************************************************************************************

  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use BoxSelector, only:Clear_BoxSelector
   use CentroSymParam_for_BCC_GPU, only:Clear_CentroSymmetry_DEV

   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_CentroSymmetry_DEV(SimBox, CtrlParam)
          call Clear_BoxSelector(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module BoxSelector_CSP


 !****************************************************************************************
 Program BoxSelector_CSP_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use BoxSelector_CSP

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_,              &
                                RECORDPROC=Record_BoxSelectorTool,  &
                                AFTRECORD=MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program BoxSelector_CSP_Main
