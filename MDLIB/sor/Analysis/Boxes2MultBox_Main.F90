 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to combine boxes into a multibox that could be used sometimes
 !                  for the convenient of visualization.
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
 !                  MultBoxSeperator.exe SETUP-file
 !                  or:
 !                  MultBoxSeperator.exe SETUP-file -T(est) t1, t2, ts -C(fg) c1, c2, c3 -B(ox) b1, b2, b3
 !
 !                  where:  SETUP-file is the name of setup file used by MD simulations in MDPSCU.
 !                  OPTIONS:
 !                        -T(est) t1, t2, t3 : indenetify the test IDs to be involved
 !                        -C(fg)  c1, c2, c3 : indenetify the configure IDs to be involved
 !                        -B(ox)  b1, b2, b3 : indenetify the boxes IDs in a TEST
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2015-04 (Hou Qing Sichuan university)
 !
 !
 module Boxes2MultBox
 use BoxSelector
 implicit none
        integer,private::m_MULTBOX0 = 1

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
         m_MULTBOX0 = CtrlParam%MULTIBOX
         CtrlParam%TIMELOOPOUT = 1

         return
   end subroutine Initialize_
  !**********************************************************************************

  !**********************************************************************************
  subroutine Record_MySelectool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_SimBoxArray, only:Putout_Instance_Config_SimBoxArray
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE
       integer, dimension(:), allocatable::BID
       integer::I, J, NB, IB0, IB1, IBS, NM, JOB,ITIME, ICFG

            JOB   = Stamp%ITest
            ITIME = Stamp%ITime
            ICFG  = Stamp%IRec(1)

            if(ITIME .lt. 0) return
               allocate(BID(size(SimBox)) )
               BID = 0
               !$$--- to create the ID of boxes
               NM = size(SimBox)/m_MULTBOX0

               IB0 = CtrlParam%STARTBOX
               if(CtrlParam%ENDBOX.le.0) then
                  IB1 = m_MULTBOX0
               else
                  IB1 = CtrlParam%ENDBOX
               end if
               IBS = max(1,CtrlParam%BOXSTEP)

               NB = 0
               do J=1, NM
                  do I= (J-1)*m_MULTBOX0+IB0, (J-1)*m_MULTBOX0+IB1, IBS
                     NB = NB + 1
                     BID(NB) = I
                  end do
               end do

             !$$--- write out summay
             call STRCATI(GFILE,  m_OUTFILE, "P", m_processid, 4)
             call STRCATI(GFILE,  GFILE, "_", 1, 4)
             call STRCATI(GFILE, GFILE, ".", ICFG, 4)

             write(*,fmt="(A)")         "Combined boxes output to "//GFILE(1:len_trim(GFILE))
             write(m_SumUNIT,fmt="(A)") "Combined boxes output to "//GFILE(1:len_trim(GFILE))
             write(m_SumUNIT,fmt="(A)") "   with the original boxes list: "
             do J= 1, NM
                do I= IB0, IB1, IBS
                   write(m_SumUNIT,fmt="(A,I7,A,I7,A,I8)") "  TEST ID:", CtrlParam%JOBID0+J-1, " BOX ID:", I, " CFG ID:", ICFG
                end do
             end do
             call STRCATI(GFILE,  m_OUTFILE, "P", m_processid, 4)
             call STRCATI(GFILE,  GFILE, "_", 1, 4)
             call Putout_Instance_Config_SimBoxArray(GFILE, SimBox, Stamp, NB=count(BID.gt.0), IB=BID)

             deallocate(BID)

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
   use BoxSelector, only:Clear_BoxSelector

   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_BoxSelector(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module Boxes2MultBox

 !***********************************************************************************
 Program MultBoxSeperator_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use Boxes2MultBox

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_,          &
                                RECORDPROC=Record_MySelectool,  &
                                AFTRECORD=MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program MultBoxSeperator_Main
