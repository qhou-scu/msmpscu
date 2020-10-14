  module PropCluster_AtomType
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to cluster atoms by their type using PropClusterCommon_14_GPU module.____
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       PropClusterCommon_GPU.F90
  !                       MD_Globle_Variables_2012_GPU.F90
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2014-10 (Hou Qing)
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
   !***  PURPOSE:  to initialize the module
   !
   use PropClusterCommon_GPU, only:Initialize_PropClusteringCommon, SetFlagProc, &
                                    m_INFILE, m_ATYP, &
                                    SelectAtomsDefault=>SelectAtoms
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I

            call SetFlagProc(SelectAtomsDefault)

            !$$--- loading control parameters
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



  end module PropCluster_AtomType
