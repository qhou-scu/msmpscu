 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to quench the instant configuration generated in MD calculations
 !                  to zero temperature.
 !                  The quenched system may be used for analysis of stable configutations.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  The usage of this prgram is similar to thta of MD simulations, with the same input
 !                  files used in MD simulation. However, the analysis section should be added in the
 !                  control file of MD simulation, that is the following section should be added
 !                  in the control file:
 !
 !                  &ANALYSUBCTL  - indicate the start of this subsection of control parameters
 !
 !                     &JOBSEL      - indicate which TESTs (ref to user guid) to be included in analysis.
 !                                    without this keyword, all TESTs will be analyzed.
 !
 !                                    usage: &JOBSEL  J0, J1, JS
 !                                           where J0 is the id of the frist TEST
 !                                                 J1 is the id of the end  TEST
 !                                                 JS is the intervale of the TESTs
 !
 !                                   example: &JOBSEL  1, 99, 2 indicating TESTs #1, #3, #5,...,#99 will be included for analyzing
 !
 !                    &CFGSEL    - indicate which configurations to be included in analysis for included TESTs
 !                                 without this keyword, all configurations will be analyzed.
 !
 !                                 usage: &CFGSEL  C0, C1, CS
 !                                        where C0 is the id of the frist configuration in a test
 !                                              C1 is the id of the end configuration in a test
 !                                              CS is the intervale of the configurations
 !
 !                                 example: &CFGSEL  5, 100, 5 indicating configuration #5, #10, ...,#100 will be included for analyzing
 !
 !                    &BOXSEL    -  indicate which boxes to be included in analysis for included TESTs
 !                                  without this keyword, all boxes in a TEST will be analyzed.
 !                                  The &BOXEL is applied when a TEST contain multiple boxes(see the parameters for keyword &BOX).
 !                                  This control parameter is not necessary for all analysis tool.
 !
 !                                  usage: &BOXSEL  B0, B1, BS
 !                                  where B0 is the id of the frist box in a test
 !                                        B1 is the id of the end box in a test
 !                                        BS is the intervale of the boxes
 !
 !                                  example: &BOXSEL  1, 50, 5 indicating boxes #1, #6, ...,#46 will be included for analyzing
 !
 !                   &QUICKDamp -  indicate if quenchng will be performed before the analysis.
 !
 !                                 usage: &QUICKDamp  steps
 !                                 where steps is the times steps needed for doing queching
 !
 !                                 example: &QUICKDamp 2000, indicating the system will be queching for 2000 timesteps beforeing
 !                                          going of analysis
 !
 !                 &ENDSUBCTL    -  indicate end of a subsection
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2013-10 (Hou Qing, Sichuan university)
 !                  update:     2014-11 (Hou Qing, Sichuan university)
 !
 module Anneal_GPU
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use MiniUtilities

 implicit none
      character(len=12),parameter, private::mp_FTAGO="&AUXF_ANNEAL"
      character(len=256), private::m_OUTFILE=""

  contains
  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate memories, and create the
   !               neighbor-list of the reference lattice
   !
   use MD_Globle_Variables,only:CreateDataFolder_Globle_Variables
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
  integer::I
  character*256::GFILE

         !$$--- loading control parameters


          if(CtrlParam%NEEDDAMP .le. 0) then
             write(*,fmt="(A,I4,A)")  ' MDPSCU Message: the number of timesteps for quenching is not supplied'
             write(*,fmt="(A, BZI6)") '                 the number of timestep will set to the default value 2000'
             call SetNeedDamp_SimMDCtrl(CtrlParam, 2000, CtrlParam%NEEDDAMPTYPE)
          end if

          m_OUTFILE = ""
           do I=1, size(CtrlParam%f_tag)
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                 m_OUTFILE = CtrlParam%f_others(I)
               end if
          end do

          if(len_trim(m_OUTFILE) .le. 0) then
            call GetPath(CtrlParam%f_geometry, GFILE)
            GFILE = GFILE(1:len_trim(GFILE)-1)//"_Quench/ "
            m_OUTFILE = GFILE(1:len_trim(GFILE))
            call GetFname(CtrlParam%f_geometry, GFILE)
            m_OUTFILE =  m_OUTFILE(1:len_trim(m_OUTFILE))//GFILE(1:len_trim(GFILE))
            write(*,fmt="(A,I4,A)")  ' MDPSCU Message: the number of timesteps for quenching is not supplied'
            write(*,fmt="(A, BZI6)") '                 the number of timestep will set to the default value 2000'
          end if
          call CreateDataFolder_Globle_Variables(m_OUTFILE)

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to
  !                  interfaced to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the
  !                  the neighbor list routine has be called.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use MD_SimBoxArray
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in) :: Stamp
       type(SimMDBox), dimension(:), intent(in) :: SimBox
       type(SimMDCtrl),              intent(in) :: CtrlParam
       !--- local variables
       character*256::GFILE

           if(Stamp%ITime .lt. 0) then
              return
           end if

           call STRCATI(GFILE,  m_OUTFILE, "P", 0, 4)
           call Putout_Instance_Config_SimBoxArray(GFILE, CtrlParam%MultOutputG, SimBox, Stamp)

          return
  end subroutine MyRECORD
  !**********************************************************************************

 end module Anneal_GPU
!************************************************************************************
 Program Anneal_GPU_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use Anneal_GPU
 implicit none

       call APPSHELL_AddRecord( PRERECORD=MyInitialize,   &
                                RECORDPROC=MyRECORD)

       call Main_ANALYSIS(PROCESSID=0) !, FORCETABLE=Register_Interaction_Table, POTTYPE="FS_TYPE")

       stop
 end  program Anneal_GPU_Main
