  module PropCluster_CoordNum
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to cluster atoms by their coordination number using PropClusterCommon_14_GPU module.
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       PropClusterCommon_GPU.F90
  !                       MD_Globle_Variables_2012_GPU.F90
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2014-11 (Hou Qing)
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
  !-----------------------------------------------
    implicit none

  contains

  !**********************************************************************************
   subroutine Initialize_PropClustering(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC
   use PropClusterCommon_GPU, only:Initialize_PropClusteringCommon, SetFlagProc, &
                                      m_INFILE, m_ATYP, &
                                      SelectAtomsDefault=>SelectAtoms
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::hFile, N, LINE, ISTR, I
   character*256::STR
   character*32::STRNUMB(10), KEYWORD
   external ISTR

            call SetFlagProc(SelectAtomsDefault)

            !$$--- loading control parameters
            call Initialize_PropClusteringCommon(SimBox, CtrlParam)

              if(m_ATYP(1).le.0) then
                write(*,fmt="(A)") " MDPSCU Error: at least one atom type should be given for running "  &
                                      //gm_ExeName(1:len_trim(gm_ExeName))
                write(*,fmt="(A)") "        Usage: add &PROP_TYPE type1, type2... in input file "//m_INFILE(1:len_trim(m_INFILE) )
                write(*,fmt="(A)") "        Process to be stopped"
                stop
              end if

      return
  end subroutine Initialize_PropClustering
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_PropClusterTool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use PropClusterCommon_GPU, only:RECORD_PropCluster_Common=>RECORD_PropCluster
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables


            if(Stamp%ITime .LT. 0) return
            call RECORD_PropCluster_Common(Stamp, SimBox, CtrlParam)

            return
  end subroutine RECORD_PropClusterTool
  !****************************************************************************************



  end module PropCluster_CoordNum
