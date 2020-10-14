 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to select out the boxes with bubbles burting.
 !                  We use Delauney tessellation to generate atomic clusters. If a cluster is
 !                  marked as OPENED, the cluster is burst.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       BoxSelector.F90
 !                       PropCluster_AtomType.F90
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
 !                  BoxSelector_BubbleBurst.exe SETUP-file
 !                  or:
 !                  BoxSelector_BubbleBurst.exe SETUP-file -T(est) t1, t2, ts -C(fg) c1, c2, c3 -B(ox) b1, b2, b3
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

 module BoxSelector_BubbleBurst
  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
 implicit none

  contains
  !**********************************************************************************
   subroutine Initialize_(SimBox, CtrlParam)
   !***  PURPOSE:  to Initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use BoxSelector, only:Initialize_BoxSelector, SetBoxSelectFlagProc=>SetSelectProc
   use PropClusterCommon_GPU, only:Initialize_PropClusteringCommon, SetClusteringFlagProc=>SetFlagProc, &
                                      m_ClusteringMeth=>m_METHOD, mp_DELAUNAY, m_INFILE, m_ATYP, &
                                      SelectAtomsDefault=>SelectAtoms
   use VoronoiTessellation_GPU, only:SetAtomSelector
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
  integer::I

         call Initialize_BoxSelector(SimBox, CtrlParam)
         call Initialize_PropClusteringCommon(SimBox, CtrlParam)
              if(all(m_ATYP.le.0) ) then
                 write(*,fmt="(A)")   " MDPSCU Warning: no atom type is specified"
                 write(*,fmt="(A)")   "                 all type of atoms to be used in clutering"
                 m_ATYP(1:SimBox%NGROUP) = 1
                 call ONWARNING(gm_OnWarning)
              end if

              do I=1, size(m_ATYP)
                 if(m_ATYP(I) .gt. 0 .and. I.gt.SimBox%NGROUP) then
                   write(*,fmt="(A, I2, A)")   " MDPSCU Warning: atom type", I, " larger than number of types in the box"
                   write(*,fmt="(A)")          "        check input for &PROP_TYPE in input file "//m_INFILE(1:len_trim(m_INFILE) )
                   call ONWARNING(gm_OnWarning)
                 end if
              end do
              m_ClusteringMeth = mp_DELAUNAY

         call SetClusteringFlagProc(SelectAtomsDefault)
         call SetBoxSelectFlagProc(MySelectProc)
         call SetAtomSelector(ASelector=MyMaskPRoc)
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
   use MD_Globle_Variables_GPU, only:dm_NPRT
   use PropClusterCommon_GPU, only:Do_Delaunay_Clustering
   use VoronoiTessellation_GPU, only: mp_VSTAT_OPENED
   implicit none
   !----   DUMMY Variables
      type(SimMDBox), dimension(:)::SimBox
      type(SimMDCtrl)             ::CtrlParam
      integer,        dimension(:)::FLAG
   !---- Local variables
      integer::NCS
      integer, dimension(:), allocatable::SCS
      integer, dimension(:), allocatable::ASQCS
      integer, dimension(:),allocatable::STATCS
      integer::I, NPRT, IB, IP, IERR

            allocate(SCS(dm_NPRT), ASQCS(dm_NPRT), STATCS(dm_NPRT), STAT=IERR)
            if(IERR) then
               write(*,fmt="(A)") "MDPSCU Error: fail to allocate working memory in BoxSelector_BubbleBurst"
               write(*,fmt="(A)") "              Process to be stopped"
               stop
            end if
           call Do_Delaunay_Clustering(SimBox, CtrlParam, NCS, SCS, ASQCS, STATCS)

            NPRT = SimBox(1)%NPRT
            FLAG = 0
            IP   = 1
            do I=1, NCS
               if(STATCS(I) .eq. mp_VSTAT_OPENED) then
                  IB = ASQCS(IP)/NPRT + 1
                  FLAG(IB) = 1
               end if
               IP = IP + SCS(I)
            end do
            deallocate(SCS, ASQCS, STATCS, STAT=IERR)

            return
  end subroutine MySelectProc
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyMaskPRoc(SimBox, CtrlParam, MASK)
  !***  PURPOSE:   to create the FLAG witch identify the atoms
  !                that will be included in Voronoi calculation
  !     INPUT:     SimBox:
  !                CtrlParam:
  !
  !     OUTPUT:    FLAG,
  !
     use MD_Globle_Variables_GPU, only:dm_NPRT, hm_XP
     implicit none
     type(SimMDBox)                    ::SimBox
     type(SimMDCtrl)                   ::CtrlParam
     integer, dimension(:), allocatable::MASK
     integer I

         return
         do I=1, dm_NPRT
            if(hm_XP(I,3) .ge.  SimBox%RR*15.0) then
                MASK(I) = 0
            else
                MASK(I) = 1
            end if
         end do
      return
  end subroutine MyMaskPRoc
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
  use BoxSelector, only:Record_BoxSelector
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables

            if(Stamp%ITime .lt. 0) return
            call Record_BoxSelector(Stamp, SimBox, CtrlParam)
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
   use PropClusterCommon_GPU, only:Clear_Clustering=>Clear

   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_BoxSelector(SimBox, CtrlParam)
          call Clear_Clustering(SimBox, CtrlParam)
          return
  end subroutine MyCleaner
  !**********************************************************************************
 end module BoxSelector_BubbleBurst

 !***********************************************************************************
 Program BoxSelector_BubbleBurst_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use BoxSelector_BubbleBurst

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_,          &
                                RECORDPROC=Record_MySelectool,  &
                                AFTRECORD=MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program BoxSelector_BubbleBurst_Main
