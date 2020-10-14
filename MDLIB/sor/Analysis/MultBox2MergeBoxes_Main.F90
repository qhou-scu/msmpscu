 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to merge boxes into a large box. The difference between this program
 !                  and Boxes2MultBox is that, in Boxes2MultBox, the generated multibox contains acutally a number
 !                  of independent boxes with the number of atoms in each of the boxes is the NPRT, in a MERGEDBOX,
 !                  one boxes to be generated with the number of atoms being NPRT*number of boxes.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !                       BoxSelector.F90
 !
 !                  SEE ALSO____________________________________________________________________________
 !                       Boxes2MultBox.F90
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one line need to be add:
 !
 !                    &AUXF_SELOUT filename
 !
 !                  in the setup file. The filename is the file name for output results.
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
 !                        -T(est) t1, t2, t3 : indentify the test IDs to be involved
 !                        -C(fg)  c1, c2, c3 : indentify the configure IDs to be involved
 !                        -B(ox)  b1, b2, b3 : indentify the boxes IDs in a TEST
 !
 !                  If you needed additional control on the output, you should add a line in
 !                  setup fileL
 !
 !                        &AUXF_SELIN filename
 !
 !                  in the AUXF_SELIN file, you can specify the types of atoms that will be output.
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2015-04 (Hou Qing Sichuan university)
 !
 !
 module MultiBox2MergeBoxes
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use AtomTrackCommon

 implicit none
        integer, private::m_MULTBOX0 = 1

  contains
  !**********************************************************************************
   subroutine Initialize_(SimBox, CtrlParam)
   !***  PURPOSE:  to Initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables

         call Initialize_AtomTrackCommon(SimBox, CtrlParam)
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
  use MD_SimBoxArray, only:Putout_Instance_Config_SimBoxArray, Release_SimBoxArray
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE
       integer, dimension(:), allocatable::BID
       integer::I, J, NB, IB0, IB1, IBS, NM, ISECT, ICFG
       type(SimMDBox)::SWAPBOX(1)
        type(MDRecordStamp)::SWAPSTAMP


            if(Stamp%ITime .lt. 0) return

               ICFG = Stamp%IRec(1)
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
                     BID(NB) = 1
                  end do
               end do
             !$$--- merge the boxes
              call MergeBoxs_SimMDBox(SimBox, SWAPBOX(1), MASK=BID)

             if(.not.allocated(m_NUMINIATOMS)) then
                 call Identify_Atoms(SWAPBOX, CtrlParam)
             end if

             !$$--- write out summay
             !call GetSectID_SimMDCtrl(CtrlParam, ITIME, ISECT)
             !call GetCfgID_SimMDCtrl(CtrlParam, ITIME, ICFG)

             call STRCATI(GFILE,  m_OUTFILE, "P", m_processid, 4)
             call STRCATI(GFILE,  GFILE, "_", 1, 4)
             call STRCATI(GFILE, GFILE, ".", ICFG, 4)

             write(*,fmt="(A, 2x, I8)") "Combined boxes output to "//GFILE(1:len_trim(GFILE))//", with num of atoms", m_NUMINIATOMS(1)
             write(m_SumUNIT,fmt="(A)") "Combined boxes output to "//GFILE(1:len_trim(GFILE))
             write(m_SumUNIT,fmt="(A)") "   with the original boxes list: "
             do J= 1, NM
                do I= IB0, IB1, IBS
                   write(m_SumUNIT,fmt="(A,I7,A,I7,A,I8)") "  TEST ID:", CtrlParam%JOBID0+J-1, " BOX ID:", I, " CFG ID:", ICFG
                end do
             end do
!             call STRCATI(GFILE,  m_OUTFILE, "P", m_processid, 4)
!             call STRCATI(GFILE,  GFILE, "_", 1, 4)
!             call Putout_Instance_Config_SimBoxArray(ITIME, ISECT, ICFG, TIME, GFILE, SimBoxSwap)
              SWAPSTAMP%ITest    = 1
              SWAPSTAMP%IBox(1)  = 1
              SWAPSTAMP%IBox(2)  = NB
              SWAPSTAMP%ITime    = Stamp%ITime
              SWAPSTAMP%Time     = Stamp%Time
              SWAPSTAMP%ICfg     = Stamp%ICfg
              SWAPSTAMP%IRec     = Stamp%IRec
              call Record_AtomTrackCommon(SWAPSTAMP, SWAPBOX, CtrlParam)

             call Release_SimBoxArray(SWAPBOX)
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
   implicit none
       type(SimMDBox)::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_AtomTrackCommon(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module MultiBox2MergeBoxes

 !***********************************************************************************
 Program MultBoxSeperator_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use MultiBox2MergeBoxes

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_,          &
                                RECORDPROC=Record_MySelectool,  &
                                AFTRECORD=MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program MultBoxSeperator_Main
