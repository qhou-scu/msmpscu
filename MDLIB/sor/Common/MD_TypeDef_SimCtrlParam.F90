  module MD_TYPEDEF_SimMDCtrl
  !***  DESCRIPTION: this module is to define the data type for Simulation Box.
  !
  !                  ______________________________________________________
  !
  !     HISORY:      HOU Qing, Mar, 2010,adopted from an MDGvar.f90( version 2000)
  !
  !                  2016-11-16, move the routines loading parameters for different methods
  !                              to seperated files:
  !                              MD_SimMDCtrl_ART.F90
  !                              MD_SimMDCtrl_GMD.F90
  !                              MD_SimMDCtrl_NEB.F90
  !                              MD_SimMDCtrl_PARREP.F90
  !                              MD_SimMDCtrl_TAD.F90
  !
  !                 2018-03-08, add the control mod for local temperature control, HOU Qing
  !
  use MD_CONSTANTS
  use MiniUtilities
  use MD_TYPEDEF_InputPaser
  use MD_TYPEDEF_EPCCtrl

  implicit none

   !--- The type of a simulation control parameter
      type::SimMDCtrl
          !***
           integer::RESTART=0                           !=0,   start a new session
                                                        !=-1,  restart from the end point of the previous calculation
                                                        !=xxx, restart from a xxx time step of the previouss calculation
          !*** the information about time step
          !
           integer::DIFFORDER = 2                       !the differential scheme order,
                                                        ! <=2, Swope scheme to be used
                                                        ! > 2, Gear scheme to be used
           integer::IT0  = 1                            !the start point of time steps in the whole time sections
           integer::IT1  = 1                            !the end  point of time steps  in the whole time sections
           integer::ITE  = 0                            !number of integration steps in a time section: ITE = IT1-IT0 + 1
           real(KINDDF)::TEMTM                          !the assumed terminal time (ps), the process will stop
                                                        !at first point of time step ITE and time TEMTIME
           integer::IHDUP=0                             != 0, use fixed time step, HMI
                                                        !> 0, the time step is increased for every IHDUP steps
                                                        !< 0, the time step adjusted according the max displacements of atoms
           real(KINDDF)::H                              !the current time step of integration in second
           real(KINDDF)::HMI                            !the minimum of the time step
           real(KINDDF)::HMX                            !the maximum of the time step
           real(KINDDF)::DMX                            !the maximum of the permitted displacement of an atoms
                                                        !used to adjust the varibale time step


           integer::SEED(4) = -1                        !the seed for random numner, <0, to use time as a seed


          !*** the information concerning periodic condition
          !
           integer::IFPD(3)  =1                         !if periodic condition a long x, y, z
           integer::MULTIBOX =1                         !the number of sub-box, this variable would be used in GPU computing
           integer::TOTALBOX =1                         !the total number of indentical boxes, this variable would be used in GPU computing
           integer::INDEPBOX =1                         !if all the box are independent

          !*** the information control temperature
          !
           integer::TI_CTRL = CP_TICTRL_GLOBAL            !--- variable indicating if temperature control to be applied
                                                          !    = CP_TICTRL_NONE,  the whole will freely relaxing without temperature control
                                                          !    = CP_TICTRL_GLOBAL the whole system to be  thermalized to given temperature
                                                          !    = CP_TICTRL_LOCAL, partitions of system to be thermalized
                                                          !
           !-- the parameters for global temperature control
           real(KINDDF)::TI                               !--- the external temperature in K

           integer::DAMPSCHEME = CP_DAMPSCHEME_DYN        !--- the damp method
                                                          !     =  CP_DAMPSCHEME_LBFGS, to damp the system using LBFGS method
                                                          !     =  CP_DAMPSCHEME_DYN,   to damp the system by set the velocities
                                                          !                             of atoms to when the force and the velocities
                                                          !                             are in difference direction
                                                          !     =  CP_DAMPSCHEME_CG     to damp the system using CG method
                                                          !     =  CP_DAMPSCHEME_ST     to damp the system using steepest method
           integer::DAMPTIME0  = 1                        !--- start time point for damping
           integer::DAMPTIME1  = 0                        !--- time steps of performing damping

           !--- the parameters contral thermalizartion. the thermalization
           !
           integer::IVSCHEME = CP_THERMALSCHEME_MC        !the thermalizing method
                                                          ! = cp_THERMALSCHEME_MC,  to thermalizing the system by assigne the velocity of atoms by Maxwell distribution
                                                          ! = cp_THERMALSCHEME_VSCAL, to thermalizing the system by sacling the velocity of atoms
                                                          ! = cp_THERMALSCHEME_PSCAL, to thermalizing the system by sacling the position of atoms
                                                          ! = cp_THERMALSCHEME_ASCAL, to thermalizing the system by sacling the velocity with a increment

           integer::IVTIME0  = 1                          !start time point of thermalization will be carried out
           integer::IVTIME   = 0                          !number of thermalization will be carried out in ONE segements
           integer::IVPAS    = 50                         !number of time steps between two thermalization

           integer::IFNOSE    = 0                         !if Nose module to be used to control the temperature
           real(KINDDF)::NoseQ                            !the Nose parameter


           !-- the parameters for local temperature or EPC control
           integer::EPCTIME0                              ! start time point for E-P coup[ling
           integer::EPCTIME1                              ! time steps of performing E-P coupling
           type(TiCtrlParam)::LT_CTRL(mp_MXGROUP)         ! local temperature control for atom types
           type(STCtrlParam)::ST_CTRL


          !*** the information control temperature and pressure
          !
           integer::IFJONSHON = 0                                !if use jonshon method to control pressure
           real(KINDDF)::JONSHONW                                !the Nose parameter
           real(KINDDF)::PEX                                     !the external pressure in CGS
           real(KINDDF)::PEXTENSOR(3,3)=0.D0                     !the external pressure tensor in CGS

           integer::IFBPC         = 0                            !if Berendsen pressure coupling is used
           character*256::BPCDATA = ""

          !*** the control parameters for neighborelist calculations
          !
           integer::NB_MXNBS       = 800                         !max number of neighbores
           integer::NB_NEAREST     = 0                           !=0, the neihbores is determined by Rcut
                                                                 !>0, only the nearest neighbore to be included in neighborelist
           real(KINDDF)::NB_RM(mp_MXGROUP,mp_MXGROUP)=0.D0       !The cutoff range for neighbore calculations
           integer::NB_UPTABMI     = 10                          !minimum time steps, the neighbore list will be updated
           integer::NB_UPTABMX     = 10                          !maxmum  time steps, the neighbore list will be updated
           integer::NB_DBITAB      = 100000                      !time steps, UPITAB to be doubled
           integer::NB_UPTAB       = 10                          !Each NB_UPITAB steps, the neighbore list will be updated

          !*** the control parameter for update activation region
          !
           integer     ::AR_UPTABMI     = -10                    !minimum time steps for updating activation region
           integer     ::AR_UPTABMX     = -10                    !maxmum  time steps for updating activation region
           integer     ::AR_DBITAB      = -100000                !time steps, AR_UPITAB to be doubled
           integer     ::AR_UPTAB       = -10                    !each AR_UPITAB steps, the activation region will be updated
           integer     ::AR_EXTEND      = 1                      !extention of activation region
           integer     ::AR_METHOD      = 0                      !method of construct the activation region
           integer     ::AR_CENTPART(mp_MXGROUP) = 0             !flag indicate the atom type to be used for activation region seeds
           real(KINDDF)::AR_EKIN        = C_ZERO                 !kinetic energy thresthold indicate the atom type to be used for activation region seeds

          !*** the control parameters for force calculations
          !
           real(KINDDF)::RU(mp_MXGROUP,mp_MXGROUP)=0.D0          !The cutoff range for interactions
           integer     ::NUMFTABR      = 10000                   !The number of points for distance in generating force table
           integer     ::NUMFTABE      = 10000                   !The number of pointe for RHO in generating embedment function ( used for EAM type potential)
           real(KINDDF)::RHOSCAL       = 20                      !The scal that determine the table range of RHO for embedment function
           integer     ::OUTPUTFTAB    = 0                       !Indicate if force table to be output

           !*** control the output of thermal quantites
           !
           integer::TimestpR      = -1                           !control parameter for performing external recoding routine
           integer::TimestpG      = 1000                         !control parameter for output of configuration
           integer::TimestpQ      = 1000                         !control parameter for output of temperature etc
           integer::TimestpSave   = 10000                        !control the period to save the current statu
           integer::MultOutputG   = 0                            !control configuration output when MULTIBOX > 1 is used
                                                                 ! =0, the MULTIBOX box to be output in a single file
                                                                 ! =1, one box for one file


           !*** control variables for the analysis application
           !    this variables take no effect in on-flight analysis
           integer::TIMELOOPOUT   =  0                   !indicator to indicate if time loop is the loop out of JOB loop
                                                         ! used in analysis tools. usuall, the parameter is modified in
                                                         ! initialization routine of analysis tools.
                                                         ! =0 (default), jobloop encloses of timeloop
                                                         ! =1, timeloop enclose job loop
           integer::JOBID0        = -1                   !the id of the starting job to be analyzed, <0, all jobs to be analysed
           integer::JOBID1        = -1                   !the id of the ending job to be analyzed, <0, all jobs to be analysed
           integer::JOBIDSTEP     =  1                   !the step of id to be analyzed, <0, all jobs to be analysed
           integer::STARTCFG      = -1                   !the start configure to be analysed
           integer::ENDCFG        = -1                   !the end configure to be analysed, <0 all configures untill termination to be analysis
           integer::CFGSTEP       =  1                   !the id step between two configres
           integer::STARTBOX      = -1                   !the start box to be analysed, if the multiple-box was used in MD simulations
           integer::ENDBOX        = -1                   !the end box to be analysed, <0 all boxes  in a multi-box ocnfoure to be analysis
           integer::BOXSTEP       =  1                   !the id step between two boxes

           integer::NEEDPOT       =  0                   !variable indicating if potential or force caluclations are required in analysis
                                                         !  modification of this variable left to user-supplied routine
           integer::NEEDDAMP      =  0                   !variable indicating if damping is needed in analysis
                                                         !  modification of this variable left to user-supplied routine
           integer::NEEDDAMPTYPE  =  CP_DAMPSCHEME_ST    !the damp method

           integer::NEEDNEIGHB    =  1                   !variable indicating if neighbore list is needed in analysis
                                                         !  modification of this variable left to user-supplied routine

           integer::NEEDNEWATOMS =  0                    !variable indicatin if new atoms to be to the initial boxes, in cases of some applications
                                                         !the number of atoms is changed by deposition, embedment new atoms

          !***  the output files
          character*256::f_quantity    = ''              !file name for output basic thermal quantities
          character*256::f_geometry    = ''              !file name for output configuarations
          character*256::f_configure   = ''              !file name for svae current configuration of the system, that can be used for restart
          character(len=32), dimension(10)::f_tag =''
          character(len=256),dimension(10)::f_others =''
          type(InputStatements):: ExInputStat(10)         !the  statement list containing control information for external
                                                          ! analysis applications


          type(SimMDCtrl), pointer::parent=>null()
          type(SimMDCtrl), pointer::next=>null()

          !*** control parameters for some specific applications **************
          !
           integer::      ExitAtEvent      = 0             ! parameter indicating if the processs should exit the time-section from IT0 to IT1.
                                                           ! this parameter could be used in PARALELL-REPLICA application.

           real(KINDDF):: STRCUT_DRTol     = 0.03D0        ! distance tolerance to check if two configurations are the same
                                                           ! if any displacement of atoms is > DRTolEvent (in LU), the config. are considered different
           integer     :: Quench_Steps    = 1000           ! number of iterration for LBFGS, Steepest or CG algorithm  
           integer     :: Quench_Meth     = CP_DAMPSCHEME_ST! method of quench

           !--- control parameters specific for  LBFGS
           real(KINDDF):: LBFGS_PGtol     = 0.D0           ! PGTOL for LBFGS optimization
           real(KINDDF):: LBFGS_Factr     = 0.D0           ! FACTR for LBFGS optimization
           integer     :: LBFGS_Msave     = 7              ! M parameter for LBFGS

           !--- control parameters specific for Steepest and CG
           real(KINDDF):: STEEPEST_Alpha   = 0.1D0         ! alpha parameter for Steepest optimization
           real(KINDDF):: STEEPEST_MiStep  = 1.D-5         ! min step (in LU) of an atom moving, if the largest movement of atoms
                                                           ! in a step is smaller then this value, convergence is considered reached
           real(KINDDF):: STEEPEST_MxStep  = 0.1D0         ! permitted max step (in LU) of an atom moving in a step
           real(KINDDF):: STEEPEST_MiDelE  = 0.001D0       ! potential tolerance (eV) when the defference between two step or two line-search loop
                                                           ! is smaller than this value, convergence is considered reached   

           !--- control parameters specific for  NEB
           character*4 :: NEB_Meth         ="NEB"
           integer     :: NEB_NumImg       = 21            ! number of images for each path between a pair of stable state
           real(KINDDF):: NEB_KValue       = 10.D0         ! spring constant in NEB
           integer     :: NEB_Climb        = 1             ! flag determining if climbing search is permitted
           integer     :: NEB_FindMin      = 1             ! flag determining if more minia searching on path to be performed
           integer     :: NEB_StartCfg     = -1            ! initial configuration ID, used for NEB calculationfor a serial connected events
           integer     :: NEB_CfgStep      = 1             ! step for intermediate configure, used for NEB calculationfor a serial connected events
           integer     :: NEB_EndCfg       = -1            ! end configuration ID, used for NEB calculationfor a serial connected events
           integer     :: NEB_OutIniCfg    = 0             ! flag determining if output initial configurations
           real(KINDDF):: NEB_LBFGS_PGtol  = 1.D-12        ! PGTOL for LBFGS optimization
           real(KINDDF):: NEB_LBFGS_Factr  = 10.D0         ! FACTR for LBFGS optimization
           integer     :: NEB_LBFGS_Msave  = 7             ! M parameter for LBFGS

          !*** the control parameter for boost-force
          !
           integer     :: Boost_Enable     = CP_DISABLE_BOOST
           real(KINDDF):: Boost_RepStr     = 0.D0              ! strength of repulsive boost potential
           real(KINDDF):: Boost_RepDist    = 0.05D0            ! length of "spring" or half width of gauss distribution
           integer     :: Boost_Timestep   = 100               ! timesteps of adding Boost potential

           !--- control parameters specific for PARREP and TAD
           integer     :: PARREP_FineR     = 10                ! step for recodrding intermediate configuration for PARREP and TAD
           integer     :: PARREP_Timestep  = 1000              ! time steps of checking transition event

           real(KINDDF):: TAD_TILow        = 0.D0              ! the low temperature to be explorated to
           real(KINDDF):: TAD_MUmin        = 1.D12             ! the miniua prefactor
           real(KINDDF):: TAD_Delta        = 0.001D0           !

           !--- control parameters specific for ART
           real(KINDDF):: ART_Alpha       = 1.1D0              ! the climbing facter for force caculation in ART
           real(KINDDF):: ART_StepLen     = 0.05D0             ! the step length for activating and moving an image in lattice unit
           real(KINDDF):: ART_StepTol     = 0.1D0              ! the tolerance when the displacement of the fastest atom is smaller than 
                                                               ! than this value, the image displaced 
           integer     :: ART_ActMeth     = CP_CENTPARTACT_ART ! the method of initally activating system
                                                               ! = CP_CENTPARTACT_ART, by atoms of type and their extention
                                                               ! = CP_USERDEFACT_ART,  by user supplying routine
           integer     :: ART_ActSeed(mp_MXGROUP) = 0          ! flag indicating the atom type to be included in initialy activation in ART
           real(KINDDF):: ART_ActExt     = 0.D0                ! extention of activation region
                                                               ! > 0, the neighboring atoms of ART_ActSeed also included in activation
                                                               ! = 0, the the atoms of ART_ActSeed are included in activation
                                                               ! < 0, only the neighboring atoms of ART_ActSeed to be included in activation
           real(KINDDF):: ART_ActPecent  = 1.0D0               ! pecentage of atoms to be activated when multiple activation is enabled
                                 
      end type SimMDCtrl

      interface assignment (=)
          module procedure CopyFrom_SimMDCtrl
      end interface

  contains
  !*********************************************************************
  !****************************************************************************
   subroutine CopyFrom_SimMDCtrl(To, From)
      !***  PURPOSE:   to copy a control parameter to another
      !
      !     INPUT:     From, the source simulation controal
      !     OUTPUT     To, the copy of FROM
      implicit none

      !--- DUMMY variables
      type(SimMDCtrl)::From, TO
      !--- Local variables
         call Copy_SimMDCtrl(From, To)
         return
   end subroutine CopyFrom_SimMDCtrl
  !*********************************************************************

  !****************************************************************************
   subroutine Copy_SimMDCtrl(From, To)
      !***  PURPOSE:   to copy a control parameter to another
      !
      !     INPUT:     From, the source simulation controal
      !     OUTPUT     To, the copy of FROM
      implicit none

      !--- DUMMY variables
      type(SimMDCtrl)::From, TO
      !--- Local variables
      integer::I

          !***
           To%RESTART        = From%RESTART

          !*** the information about time step
          !
           To%DIFFORDER      = From%DIFFORDER
           To%IT0            = From%IT0
           To%IT1            = From%IT1
           To%ITE            = From%ITE
           To%TEMTM          = From%TEMTM

           To%IHDUP          = From%IHDUP
           To%H              = From%H
           To%HMI            = From%HMI
           To%HMX            = From%HMX
           To%DMX            = From%DMX

           To%SEED           = From%SEED

          !*** the information concerning periodic condition
          !
           To%IFPD           = From%IFPD
           To%MULTIBOX       = From%MULTIBOX
           To%TOTALBOX       = From%TOTALBOX
           To%INDEPBOX       = From%INDEPBOX

          !*** the information control temperature
          !
           To%TI_CTRL        = From%TI_CTRL
           To%TI             = From%TI
           To%DAMPSCHEME     = From%DAMPSCHEME
           To%DAMPTIME0      = From%DAMPTIME0
           To%DAMPTIME1      = From%DAMPTIME1

           To%IVSCHEME       = From%IVSCHEME
           To%IVTIME0        = From%IVTIME0
           To%IVTIME         = From%IVTIME
           To%IVPAS          = From%IVPAS

           To%EPCTIME0       = From%EPCTIME0
           To%EPCTIME1       = From%EPCTIME1
           !To%EPCPHYS       = From%EPCPHYS

           To%IFNose         = From%IFNose
           To%NoseQ          = From%NoseQ

           To%LT_CTRL        = From%LT_CTRL
           To%ST_CTRL        = From%ST_CTRL

          !*** the information control pressure
          !
           To%PEX            = From%PEX
           To%PEXTENSOR      = From%PEXTENSOR
           To%IFJONSHON      = From%IFJONSHON
           To%JONSHONW       = From%JONSHONW
           To%IFBPC          = From%IFBPC
           To%BPCDATA        = From%BPCDATA

          !*** the information for neighbore table
          !
           To%NB_RM          = From%NB_RM
           To%NB_MXNBS       = From%NB_MXNBS
           To%NB_NEAREST     = From%NB_NEAREST
           To%NB_UPTABMI     = From%NB_UPTABMI
           To%NB_UPTABMX     = From%NB_UPTABMX
           To%NB_DBITAB      = From%NB_DBITAB
           To%NB_UPTAB       = From%NB_UPTAB

          !*** the information for activation region
          !
           To%AR_UPTABMI     = From%AR_UPTABMI
           To%AR_UPTABMX     = From%AR_UPTABMX
           To%AR_DBITAB      = From%AR_DBITAB
           To%AR_UPTAB       = From%AR_UPTAB
           To%AR_EXTEND      = From%AR_EXTEND
           To%AR_METHOD      = From%AR_METHOD
           To%AR_CENTPART    = From%AR_CENTPART
           To%AR_EKIN        = From%AR_EKIN


          !*** the information for  force calculations
          !
           To%NUMFTABR       = From%NUMFTABR
           To%NUMFTABE       = From%NUMFTABE
           To%RHOSCAL        = From%RHOSCAL
           To%RU             = From%RU
           To%OUTPUTFTAB     = From%OUTPUTFTAB


          !*** control the output of thermal quantites
           To%TimestpR        = From%TimestpR
           To%TimestpG        = From%TimestpG
           To%TimestpQ        = From%TimestpQ
           To%TimestpSave     = From%TimestpSave
           To%MultOutputG     = From%MultOutputG

          !*** control the analysis
           To%JOBID0          = From%JOBID0
           To%JOBID1          = From%JOBID1
           To%JOBIDSTEP       = From%JOBIDSTEP
           To%STARTCFG        = From%STARTCFG
           To%ENDCFG          = From%ENDCFG
           To%CFGSTEP         = From%CFGSTEP
           To%STARTBOX        = From%STARTBOX
           To%ENDBOX          = From%ENDBOX
           To%BOXSTEP         = From%BOXSTEP
           To%NEEDPOT         = From%NEEDPOT
           To%NEEDDAMP        = From%NEEDDAMP
           To%NEEDDAMPTYPE    = From%NEEDDAMPTYPE
           To%NEEDNEIGHB      = From%NEEDNEIGHB
           To%NEEDNEWATOMS    = From%NEEDNEWATOMS


          !***  the output files
           To%f_quantity      = From%f_quantity
           To%f_geometry      = From%f_geometry
           To%f_configure     = From%f_configure
           To%f_tag           = From%f_tag
           To%f_others        = From%f_others
           do I=1, size(From%ExInputStat)
              call Copy_InputStatements(From%ExInputStat(I), To%ExInputStat(I))
           end do
          !*** parameter conrolling Event detect
           To%ExitAtEvent     = From%ExitAtEvent
           To%STRCUT_DRTol    = From%STRCUT_DRTol
           To%Quench_Steps    = From%Quench_Steps
           To%Quench_Meth     = From%Quench_Meth
           !*** parameter conrolling LBFGS
           To%LBFGS_PGtol     = From%LBFGS_PGtol
           To%LBFGS_Factr     = From%LBFGS_Factr
           To%LBFGS_Msave     = From%LBFGS_Msave

          !*** parameter conrolling Steepest method
           To%STEEPEST_Alpha  = From%STEEPEST_Alpha        
           To%STEEPEST_MiStep = From%STEEPEST_MiStep      
           To%STEEPEST_MxStep = From%STEEPEST_MxStep    
           To%STEEPEST_MiDelE = From%STEEPEST_MiDelE    

          !*** parameters controlling NEB 
           To%NEB_LBFGS_PGtol = From%NEB_LBFGS_PGtol
           To%NEB_LBFGS_Factr = From%NEB_LBFGS_Factr
           To%NEB_LBFGS_Msave = From%NEB_LBFGS_Msave
           To%NEB_Meth        = From%NEB_Meth
           To%NEB_NumImg      = From%NEB_NumImg
           To%NEB_KValue      = From%NEB_KValue
           To%NEB_Climb       = From%NEB_Climb
           To%NEB_FindMin     = From%NEB_FindMin
           To%NEB_StartCfg    = From%NEB_StartCfg
           To%NEB_CfgStep     = From%NEB_CfgStep
           To%NEB_EndCfg      = From%NEB_EndCfg
           To%NEB_OutIniCfg   = From%NEB_OutIniCfg

          !*** parameters for  boost force
          !
           To%Boost_Enable    = From%Boost_Enable
           To%Boost_RepStr    = From%Boost_RepStr
           To%Boost_RepDist   = From%Boost_RepDist
           To%Boost_Timestep  = From%Boost_Timestep


          !*** parameters controlling PAREP and TAD 
           To%PARREP_FineR    = From%PARREP_FineR
           To%PARREP_Timestep = From%PARREP_Timestep

           To%TAD_TILow       = From%TAD_TILow
           To%TAD_MUmin       = From%TAD_MUmin
           To%TAD_Delta       = From%TAD_Delta

          !*** parameters controlling ART
           To%ART_Alpha       = From%ART_Alpha
           To%ART_StepLen     = From%ART_StepLen
           To%ART_StepTol     = From%ART_StepTol
           To%ART_ActMeth     = From%ART_ActMeth
           To%ART_ActSeed     = From%ART_ActSeed
           To%ART_ActExt      = From%ART_ActExt                   
           To%ART_ActPecent   = From%ART_ActPecent

          return
  end subroutine Copy_SimMDCtrl
  !****************************************************************************


  !****************************************************************************
  recursive subroutine AddNext_SimMDCtrl(parent, CtrlParam)
  !***  PURPOSE:   to add next control parameter to parent parameter
  !
  !     INPUT:     parent,    the head of the control parameter list
  !                CtrlParam, the head of the control parameter structure to be appended on the list
  !     OUTPUT     parent,
  !
     type(SimMDCtrl), target::parent
     type(SimMDCtrl), pointer::CtrlParam

          if(.not.associated(parent%next)) then
             parent%next => CtrlParam
             CtrlParam%parent => parent
          else
             call AddNext_SimMDCtrl(parent%next, CtrlParam)
          end if
          return
  end subroutine AddNext_SimMDCtrl
  !****************************************************************************


  !****************************************************************************
  recursive subroutine Number_SimMDCtrl(CtrlParam, num)
  !***  PURPOSE:   to get number control parameter structure in the list
  !
  !     INPUT:     CtrlParam,  the head of the control parameter list
  !     OUTPUT     num,        the number of parameter structure in the list
  !
  implicit none
     type(SimMDCtrl)::CtrlParam
     integer::num, n

          num = 1
          if(associated(CtrlParam%next)) then
            call Number_SimMDCtrl(CtrlParam%next, n)
            num = num + n
          end if
          return
  end subroutine Number_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  recursive subroutine GetNext_SimMDCtrl(CtrlParam, id, res)
  !***  PURPOSE:   to get the IDth control parameter structure in the list
  !
  !     INPUT:     CtrlParam,  the head of the control parameter list
  !                id,         the number of level of next
  !                            if id = 0, res is the CtrlParam
  !     OUTPUT     res,      the pointer to the next id control parameter
  !
  implicit none
     type(SimMDCtrl), target::CtrlParam
     type(SimMDCtrl), pointer::res
     integer::id, cur

          if(id .eq. 0) then
             res=>CtrlParam
             return
          end if

          if(associated(CtrlParam%next)) then
             cur = id -1
             call GetNext_SimMDCtrl(CtrlParam%next, cur, res)
          else
            res => null()
          end if
          return
  end subroutine GetNext_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  recursive subroutine GetNextByITime_SimMDCtrl(CtrlParam, ITIME, res)
  !***  PURPOSE:   to get the control parameter structure at ITIME
  !
  !     INPUT:     CtrlParam,  the head of the control parameter list
  !                ITIME,      the number time steps
  !
  !     OUTPUT     res,      the pointer to the next id control parameter
  !
  implicit none
     type(SimMDCtrl), target     ::CtrlParam
     integer,         intent(in) ::ITIME
     type(SimMDCtrl), pointer    ::res

          res=>null()
          if(ITIME .eq. 0)then
             res => CtrlParam
             return
          end if

          if(ITIME .le. CtrlParam%IT1 .and. ITIME .ge. CtrlParam%IT0 ) then
             res => CtrlParam
             return
          else
             if(associated(CtrlParam%next)) then
               call GetNextByITime_SimMDCtrl(CtrlParam%next, ITIME, res)
               return
             end if
          end if
          return
  end subroutine GetNextByITime_SimMDCtrl
  !****************************************************************************


  !****************************************************************************
  recursive subroutine NumberCfg_SimMDCtrl(CtrlParam, num)
  !***  PURPOSE:   to get the number of configure file to be created using
  !                this control parameter list
  !
  !     INPUT:     CtrlParam, the head of the list
  !
  !     OUTPUT     CtrlParam
  implicit none
     type(SimMDCtrl)::CtrlParam
     integer::num, n

          if(CtrlParam%TimestpG .gt. 0) then
              num = (CtrlParam%IT1 - CtrlParam%IT0 + 1)/CtrlParam%TimestpG
          else
              num = 0
          end if

          if(associated(CtrlParam%next)) then
            call NumberCfg_SimMDCtrl(CtrlParam%next, n)
            num = num + n
          end if
          return
  end subroutine NumberCfg_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  recursive subroutine NumberRec_SimMDCtrl(CtrlParam, num)
  !***  PURPOSE:   to get the number of records to be performed using
  !                this control parameter list
  !
  !     INPUT:     CtrlParam, the head of the list
  !
  !     OUTPUT     CtrlParam
  implicit none
     type(SimMDCtrl)::CtrlParam
     integer::num, n

          if(CtrlParam%TimestpR .gt. 0) then
              num = (CtrlParam%IT1 - CtrlParam%IT0 + 1)/CtrlParam%TimestpR
          else
              num = 0
          end if

          if(associated(CtrlParam%next)) then
            call NumberRec_SimMDCtrl(CtrlParam%next, n)
            num = num + n
          end if
          return
  end subroutine NumberRec_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetCfgITime_SimMDCtrl(CtrlParam, ICFG, ITIME)
  !***  PURPOSE:   to get the ITIME corresponding to ICFG
  !
  !     INPUT:     CtrlParam, the head of the list
  !                ICFG,      the configure ID
  !
  !     OUTPUT     ITIEM,     the time step corresponding to ICFG
  implicit none
     type(SimMDCtrl)       ::CtrlParam
     integer, intent(inout)::ICFG
     integer, intent(out)  ::ITIME
     !--- local variables
     integer::num, nextid, IT

          IT = ICFG
          if(CtrlParam%TimestpG .gt. 0) then
              num = (CtrlParam%IT1 - CtrlParam%IT0 + 1)/CtrlParam%TimestpG
          end if

          if(IT - num .le. 0) then
             ITIME = CtrlParam%TimestpG*IT + CtrlParam%IT0 - 1
             return
          else
             IT = IT - num
             if(associated(CtrlParam%next)) then
               call GetCfgITime_SimMDCtrl(CtrlParam%next, IT, ITIME)
             end if
             return
          end if

          return
  end subroutine GetCfgITime_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetRecITime_SimMDCtrl(CtrlParam, IREC, ITIME)
  !***  PURPOSE:   to get the ITIME corresponding to IREC
  !
  !     INPUT:     CtrlParam, the head of the list
  !                IREC,      the configure ID
  !
  !     OUTPUT     ITIEM,     the time step corresponding to ICFG
  implicit none
     type(SimMDCtrl)       ::CtrlParam
     integer, intent(inout)::IREC
     integer, intent(out)  ::ITIME
     !--- local variables
     integer::num, nextid, IT

          IT = IREC
          if(CtrlParam%TimestpR .gt. 0) then
              num = (CtrlParam%IT1 - CtrlParam%IT0 + 1)/CtrlParam%TimestpR
          end if

          if(IT - num .le. 0) then
             ITIME = CtrlParam%TimestpR*IT + CtrlParam%IT0 - 1
             return
          else
             IT = IT - num
             if(associated(CtrlParam%next)) then
               call GetRecITime_SimMDCtrl(CtrlParam%next, IT, ITIME)
             end if
             return
          end if

          return
  end subroutine GetRecITime_SimMDCtrl
 !****************************************************************************


 !****************************************************************************
  recursive subroutine GetCfgID_SimMDCtrl(CtrlParam, ITIME, ID)
  !***  PURPOSE:   to get the ID of a configure file at time step ITIME
  !
  !     INPUT:     CtrlParam, the head of the list
  !                ITIME,     the number of time steps
  !
  !     OUTPUT     ID,        the ID of configure
  implicit none
     type(SimMDCtrl)       ::CtrlParam
     integer, intent(in)   ::ITIME
     integer, intent(inout)::ID
     !--- local variables
     integer::num, nextid

          if(ITIME .le. CtrlParam%IT1 ) then
             if(CtrlParam%TimestpG .gt. 0) then
                num = (ITIME - CtrlParam%IT0 + 1)/CtrlParam%TimestpG
             end if
             ID = num
             return
          else
             num = (CtrlParam%IT1 - CtrlParam%IT0 + 1)/CtrlParam%TimestpG
             nextid = 0
             if(associated(CtrlParam%next)) then
               call GetCfgID_SimMDCtrl(CtrlParam%next, ITIME, nextid)
             end if
             ID = num + nextid
          end if

          return
  end subroutine GetCfgID_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetRecID_SimMDCtrl(CtrlParam, ITIME, ID)
  !***  PURPOSE:   to get the ID of a recording at time step ITIME
  !
  !     INPUT:     CtrlParam, the head of the list
  !                ITIME,     the number of time steps
  !
  !     OUTPUT     ID,        the ID of configure
  implicit none
     type(SimMDCtrl)       ::CtrlParam
     integer, intent(in)   ::ITIME
     integer, intent(inout)::ID
     !--- local variables
     integer::num, nextid

          if(ITIME .le. CtrlParam%IT1 ) then
             if(CtrlParam%TimestpR .gt. 0) then
                num = (ITIME - CtrlParam%IT0 + 1)/CtrlParam%TimestpR
             else 
                num = 0   
             end if
             ID = num
             return
          else
             if(CtrlParam%TimestpR .gt. 0) then
               num = (CtrlParam%IT1 - CtrlParam%IT0 + 1)/CtrlParam%TimestpR
             else 
               num = 0
             end if      
             nextid = 0
             if(associated(CtrlParam%next)) then
               call GetRecID_SimMDCtrl(CtrlParam%next, ITIME, nextid)
             end if
             ID = num + nextid
          end if

          return
  end subroutine GetRecID_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  recursive subroutine GetSectID_SimMDCtrl(CtrlParam, ITIME, ID)
  !***  PURPOSE:   to get the ID of a time section at time step ITIME
  !
  !     INPUT:     CtrlParam, the head of the list
  !                ITIME,     the number of time step
  !
  !     OUTPUT     ID,        the ID tof time section
  implicit none
     type(SimMDCtrl)       ::CtrlParam
     integer, intent(in)   ::ITIME
     integer, intent(inout)::ID
     !--- local variables
     integer::nextid

          ID = 0
          if(ITIME .eq. 0) return

          if(ITIME .le. CtrlParam%IT1 .and. ITIME .ge. CtrlParam%IT0 ) then
             ID = ID + 1
             return
          else
             nextid = 0
             if(associated(CtrlParam%next)) then
               call GetSectID_SimMDCtrl(CtrlParam%next, ITIME, nextid)
               ID = nextid + 1
             end if
          end if

          return
  end subroutine GetSectID_SimMDCtrl
  !****************************************************************************


 !****************************************************************************
  recursive subroutine TotalTSteps_SimMDCtrl(CtrlParam, timesteps)
  !***  PURPOSE:   to get the total number of time steps
  !
  !     INPUT:     CtrlParam, the head of the list
  !
  !     OUTPUT     CtrlParam
  implicit none
     type(SimMDCtrl)::CtrlParam
     integer::timesteps
     !--- local variables

          if(associated(CtrlParam%next)) then
            call TotalTSteps_SimMDCtrl(CtrlParam%next, timesteps)
          else
            timesteps =  CtrlParam%IT1
          end if
          return
  end subroutine TotalTSteps_SimMDCtrl
  !****************************************************************************


  !****************************************************************************
  recursive subroutine Release_SimMDCtrl(CtrlParam)
  !***  PURPOSE:   to release the memories when add the control parameters in the list
  !
  !     INPUT:     CtrlParam,  the head of the list
  !
  !     OUTPUT     CtrlParam
  implicit none
     type(SimMDCtrl)::CtrlParam
     !---
     integer::I

          if(associated(CtrlParam%next)) then
            call Release_SimMDCtrl(CtrlParam%next)
            deallocate(CtrlParam%next)
            CtrlParam%next => null()
          end if

          do I=1, size(CtrlParam%ExInputStat)
             call Release_InputStatements(CtrlParam%ExInputStat(I))
          end do

          return
  end subroutine Release_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  recursive subroutine Default_Parameter_SimMDCtrl(CtrlParam)
  !***  PURPOSE:   to set the default value of control parameters
  !
  !     INPUT:     CtrlParam,  the head of the list
  !
  !     OUTPUT     CtrlParam
  implicit none
     type(SimMDCtrl)::CtrlParam


          call Release_SimMDCtrl(CtrlParam)
          !*** the information about time step
          !
           CtrlParam%DIFFORDER = 2
           CtrlParam%IT0   = 1
           CtrlParam%IT1   = 1
           CtrlParam%ITE   = 0
           CtrlParam%TEMTM = 1.D60
           CtrlParam%IHDUP = 0
           CtrlParam%H     = 1.D0
           CtrlParam%HMI   = 1.D0
           CtrlParam%HMX   = 1.D0
           CtrlParam%DMX   = 0.1D0

          !*** the information concerning periodic condition
          !
           CtrlParam%IFPD     =1
           CtrlParam%MULTIBOX =1
           CtrlParam%TOTALBOX =1
           CtrlParam%INDEPBOX =1

          !*** the information control temperature and pressure
          !
           CtrlParam%TI_CTRL    = CP_TICTRL_GLOBAL
           CtrlParam%TI         = 0.D0
           CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_DYN
           CtrlParam%DAMPTIME0  = 1
           CtrlParam%DAMPTIME1  = 0

           !--- E-P coupling control
           CtrlParam%EPCTIME0  = 0
           CtrlParam%EPCTIME1  = 0
           !CtrlParam%EPCPHYS  = ""

           !--- For local temparature control
           call Default_TiCtrlParam(CtrlParam%LT_CTRL)
           call Default_STCtrlParam(CtrlParam%ST_CTRL)

           !*** the parameters contral thermalizartion. the thermalization
           CtrlParam%IVSCHEME=  CP_THERMALSCHEME_MC
           CtrlParam%IVTIME0 = 1
           CtrlParam%IVTIME  = 0
           CtrlParam%IVPAS   = 50

           CtrlParam%IFNOSE  = 0
           CtrlParam%NoseQ   = 1.D20

          !*** the information control temperature and pressure
          !
           CtrlParam%IFJONSHON = 0
           CtrlParam%JONSHONW  = 1.D64
           CtrlParam%PEX       = 0.D0
           CtrlParam%PEXTENSOR =0.D0

           CtrlParam%IFBPC   = 0
           CtrlParam%BPCDATA = ""

          !*** the information for neighbore table
          !
           CtrlParam%NB_MXNBS       = 800
           CtrlParam%NB_NEAREST     = 0
           CtrlParam%NB_RM          = 0.D0
           CtrlParam%NB_UPTABMI     = 10
           CtrlParam%NB_UPTABMX     = 10
           CtrlParam%NB_DBITAB      = 100000
           CtrlParam%NB_UPTAB       = 10

          !*** the information for activation region
          !
           CtrlParam%AR_UPTABMI     = -CtrlParam%NB_UPTABMI
           CtrlParam%AR_UPTABMX     = -CtrlParam%NB_UPTABMX
           CtrlParam%AR_DBITAB      = -CtrlParam%NB_DBITAB
           CtrlParam%AR_UPTAB       = -CtrlParam%NB_UPTAB
           CtrlParam%AR_EXTEND      = 1
           CtrlParam%AR_METHOD      = 0
           CtrlParam%AR_CENTPART    = 0
           CtrlParam%AR_EKIN        = 0.D0

          !*** the information for force calculations
          !
           CtrlParam%RU         = 0.D0
           CtrlParam%NUMFTABR   = 10000
           CtrlParam%NUMFTABE   = 10000
           CtrlParam%RHOSCAL    =  20
           CtrlParam%OUTPUTFTAB = 0

           !*** control the output of thermal quantites
           CtrlParam%TimestpR      = -1
           CtrlParam%TimestpG      = 1000
           CtrlParam%TimestpQ      = 1000
           CtrlParam%TimestpSave   = 1000
           CtrlParam%MultOutputG   = 1

           CtrlParam%JOBID0        = -1
           CtrlParam%JOBID1        = -1
           CtrlParam%JOBIDSTEP     =  1
           CtrlParam%STARTCFG      = -1
           CtrlParam%ENDCFG        = -1
           CtrlParam%CFGSTEP       =  1
           CtrlParam%STARTBOX      = -1
           CtrlParam%ENDBOX        = -1
           CtrlParam%BOXSTEP       =  1

           CtrlParam%NEEDPOT       = 0
           CtrlParam%NEEDDAMP      = 0
           CtrlParam%NEEDDAMPTYPE  = CP_DAMPSCHEME_LBFGS
           CtrlParam%NEEDNEIGHB    = 1
           CtrlParam%NEEDNEWATOMS  = 0

          !***  the output files
           CtrlParam%f_quantity    = ''
           CtrlParam%f_geometry    = ''
           CtrlParam%f_configure   = ''
           CtrlParam%f_tag         =''
           CtrlParam%f_others      =''

          !*** other parameters for some specific application
           CtrlParam%ExitAtEvent     = 0
           CtrlParam%STRCUT_DRTol    = 3.D-2
           CtrlParam%Quench_Steps    = 1000
           CtrlParam%Quench_Meth     = CP_DAMPSCHEME_ST
           !--- for LBFGS
           CtrlParam%LBFGS_PGtol     = 0.D0
           CtrlParam%LBFGS_Factr     = 0.D0
           CtrlParam%LBFGS_Msave     = 7

          !--- for Steepest
           CtrlParam%STEEPEST_Alpha  = 1.D-1        
           CtrlParam%STEEPEST_MiStep = 1.D-3     
           CtrlParam%STEEPEST_MxStep = 1.D-1         
           CtrlParam%STEEPEST_MiDelE = 1.D-3

           !--- for NEB
           CtrlParam%NEB_LBFGS_PGtol = 1.D-12
           CtrlParam%NEB_LBFGS_Factr = 1.D+1
           CtrlParam%NEB_LBFGS_Msave = 7

           CtrlParam%NEB_Meth        = "NEB"
           CtrlParam%NEB_NumImg      = 21
           CtrlParam%NEB_KValue      = 10.D0
           CtrlParam%NEB_Climb       = 1
           CtrlParam%NEB_FindMin     = 1
           CtrlParam%NEB_StartCfg    = -1
           CtrlParam%NEB_CfgStep     =  1
           CtrlParam%NEB_EndCfg      = -1
           CtrlParam%NEB_OutIniCfg   = 0

           CtrlParam%Boost_Enable    = CP_DISABLE_BOOST
           CtrlParam%Boost_RepStr    = 0.D0
           CtrlParam%Boost_RepDist   = 0.05D0
           CtrlParam%Boost_Timestep  = 1000

           CtrlParam%PARREP_FineR    = 10
           CtrlParam%PARREP_Timestep = 1000
           CtrlParam%TAD_TILow       = 0.D0
           CtrlParam%TAD_MUmin       = 1.D12
           CtrlParam%TAD_Delta       = 0.001D0

           CtrlParam%ART_Alpha       = 1.1D0              
           CtrlParam%ART_StepLen     = 0.05D0  
           CtrlParam%ART_StepTol     = 0.1D0      
           CtrlParam%ART_ActMeth     = CP_CENTPARTACT_ART 
           CtrlParam%ART_ActSeed(mp_MXGROUP) = 0          
           CtrlParam%ART_ActExt      = 0                   
           CtrlParam%ART_ActPecent   = 1.D0

          return
  end subroutine Default_Parameter_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
   subroutine SetNeedDamp_SimMDCtrl(CtrlParam, damp, scheme)
   !***  PURPOSE:   to set the NEEDDAMP parameter for all ctrlParam of the
   !                 the NEEDDAMP parameter could be used in analysis tools
      !
      !     INPUT:     damp,      the time steps for time
      !                scheme,    the scheme of performing damp
      !     OUTPUT     CtrlParam, the control parametert with NEEDDAMP updated
      implicit none

      !--- DUMMY variables
      type(SimMDCtrl), target::CtrlParam
      integer                ::damp, scheme

      !--- Local variables
      type(SimMDCtrl), pointer::nextp, next

           nextp=>CtrlParam
           do while(.TRUE.)
               nextp%NEEDDAMP      = damp
               nextp%NEEDDAMPTYPE  = scheme
               call GetNext_SimMDCtrl(nextp, 1, next)
               if(.not. associated(next)) exit
               nextp=>next
           end do
           return
   end subroutine SetNeedDamp_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
   subroutine SetMultiBox_SimMDCtrl(CtrlParam, Multibox)
   !***  PURPOSE:   to set the MULTBOX parameter for all CtrlParam in the list
      !
      !     INPUT:     Multbox,   the value to be set for MULTBOX
      !
      !     OUTPUT     CtrlParam, the control parametert with NEEDDAMP updated
      implicit none

      !--- DUMMY variables
      type(SimMDCtrl), target::CtrlParam
      integer                ::Multibox

      !--- Local variables
      type(SimMDCtrl), pointer::nextp, next

           nextp=>CtrlParam
           do while(.TRUE.)
               nextp%MULTIBOX     = Multibox
               call GetNext_SimMDCtrl(nextp, 1, next)
               if(.not. associated(next)) exit
               nextp=>next
           end do
           return
   end subroutine SetMultiBox_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_PressureCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for pressure control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(SimMDCtrl)    ::CtrlParam
     character*(*)      ::STR
     integer            ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Warning: unknown keyword in &PRESSSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&PRESSURE")
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%PEX = DRSTR(STRNUMB(1))
                           CtrlParam%PEXTENSOR   = C_ZERO
                           CtrlParam%PEXTENSOR(1,1) = CtrlParam%PEX
                           CtrlParam%PEXTENSOR(2,2) = CtrlParam%PEX
                           CtrlParam%PEXTENSOR(3,3) = CtrlParam%PEX
                     case ("&JONHSON_PARA")
                           !$$*** To get the Jonhson parameter for pressure
                           call Extract_Numb(STR,2,n,STRNUMB)
                           CtrlParam%IFJONSHON = ISTR(STRNUMB(1))
                           CtrlParam%JONSHONW = DRSTR(STRNUMB(2))
                     case ("&BERENDSEN_PARA")
                           !$$*** To get if Berendsen pressure couping is used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%IFBPC = ISTR(STRNUMB(1))

                           if(CtrlParam%IFBPC) then
                             !$$--- get characteristic time of coupling
                              call Extract_Substr(STR,1,n,STRTMP)
                             CtrlParam%BPCDATA = STRTMP(1)
                           end if

              end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading pressure control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_PressureCtl_SimMDCtrl
 !****************************************************************************

   !****************************************************************************
  subroutine Load_BoostCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for pressure control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,        intent(in)::hFile
     type(SimMDCtrl)           ::CtrlParam
     character*(*)             ::STR
     integer                   ::LINE
     !--- local variables
      character*32::STRNUMB(2),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           call Extract_SubStr(STR,2,N,STRNUMB)
           if(N .le. 0) then
              CtrlParam%Boost_Enable = CP_DISABLE_BOOST
           else if(N .gt. 0) then
              call UpCase(STRNUMB(1))
              select case (STRNUMB(1))
                     case ("DISABLE")  
                          CtrlParam%Boost_Enable  = CP_DISABLE_BOOST
                     case ("SPRING")
                           CtrlParam%Boost_Enable = CP_SPRING_BOOST
                     case ("NORMAL", "GAUSS")
                           CtrlParam%Boost_Enable = CP_GAUSS_BOOST
                     case ("BOND")
                           CtrlParam%Boost_Enable = CP_BOND_BOOST
              end select       
           end if

           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Warning: unknown keyword in &BOOSTSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)
                     case ("&BOOSTTYPE")
                        call Extract_SubStr(STR,2,N,STRNUMB)
                        call UpCase(STRNUMB(1))
                        select case (STRNUMB(1))
                            case ("DISABLE")  
                                 CtrlParam%Boost_Enable  = CP_DISABLE_BOOST
                            case ("SPRING")
                                  CtrlParam%Boost_Enable = CP_SPRING_BOOST
                            case ("NORMAL", "GAUSS")
                                  CtrlParam%Boost_Enable = CP_GAUSS_BOOST
                            case ("BOND")
                                  CtrlParam%Boost_Enable = CP_BOND_BOOST
                        end select       
       
                     case ("&STRENGTH","&KVAL")
                           call Extract_Numb(STR,2,n,STRNUMB)
                           CtrlParam%Boost_RepStr  = drstr(STRNUMB(1)) 

                     case ("&DISPLACEMENT", "&DISP")
                           call Extract_Numb(STR,2,n,STRNUMB)
                           CtrlParam%Boost_RepDist  = drstr(STRNUMB(1)) 
              end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading boost control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_BoostCtl_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  subroutine Load_BoundCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for box control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Error: unknown keyword in &BOUNDSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop
                     case ("&PERIDIC")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,3,n,STRNUMB)
                           if(n.lt. 3) then
                              write(*,fmt="('MDPSCU Error: peroidic conditions should be given for three direction')")
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              stop
                           end if
                           CtrlParam%IFPD(1) = ISTR(STRNUMB(1))
                           CtrlParam%IFPD(2) = ISTR(STRNUMB(2))
                           CtrlParam%IFPD(3) = ISTR(STRNUMB(3))

              end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading box control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_BoundCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_NeighbCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, NGROUP, LINE)
  !***  PURPOSE:   to load the control parameter for neighbor list calculation
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,          intent(in)::hFile
     type(SimMDCtrl)             ::CtrlParam
     character*(*)               ::STR
     integer                     ::LINE, NGROUP
     !--- local variables
      character*32::STRNUMB(mp_MXGROUP*mp_MXGROUP),KEYWORD
      integer::I, J, IJ, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Warning: unknown keyword in &NEIGHBSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&MAXNB")
                          !$$*** To get largest permitted number of neighbores
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NB_MXNBS = ISTR(STRNUMB(1))

                     case ("&NEARESTNB")
                          !$$*** To getpermitted  number of nearest neighbores
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NB_NEAREST = ISTR(STRNUMB(1))

                     case ("&UPDATEFRE", "&UPDATESTP")
                         !$$*** To get how many time steps between updates of neighboure list
                          call Extract_Numb(STR,3,n,STRNUMB)
                          if(n.ge.3) then
                             CtrlParam%NB_UPTABMI = ISTR(STRNUMB(1))
                             CtrlParam%NB_UPTABMX = ISTR(STRNUMB(2))
                             CtrlParam%NB_DBITAB  = ISTR(STRNUMB(3))
                          else if(n.eq.2) then
                             CtrlParam%NB_UPTABMI = ISTR(STRNUMB(1))
                             CtrlParam%NB_UPTABMX = ISTR(STRNUMB(2))
                             CtrlParam%NB_DBITAB  = 100000
                          else if(n.eq.1) then
                             CtrlParam%NB_UPTABMI = ISTR(STRNUMB(1))
                             CtrlParam%NB_UPTABMX = CtrlParam%NB_UPTABMI
                             CtrlParam%NB_DBITAB  = 100000
                          end if
                          CtrlParam%NB_UPTAB = CtrlParam%NB_UPTABMI

                    case ("&CUTOFF")
                         !$$*** To get cutoff range for neighbores
                          call Extract_Numb(STR,NGROUP*NGROUP,N,STRNUMB)
                          IJ=0
                          do I=1, NGROUP
                             do J=1, NGROUP
                                IJ = IJ + 1
                                if(IJ .GT. N) IJ = N
                                CtrlParam%NB_RM(I,J) = DRSTR(STRNUMB(IJ))
                             end do
                          end do
              end select
           end do

           return
  !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading control parameters for neigboring calculation."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_NeighbCutoffCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_ActiveRegionCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for activation region calculation
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,          intent(in)::hFile
     type(SimMDCtrl)             ::CtrlParam
     character*(*)               ::STR
     integer                     ::LINE
     !--- local variables
      character*32::STRNUMB(20),KEYWORD
      integer::I, J, N
     !----

           !*** start the press step control controlling parametgers
           call Extract_Substr(STR,1,N,STRNUMB)
           if(N .le. 0) then
              CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
           else
              call UpCase(STRNUMB(1))
              select case(STRNUMB(1))
                     case("ENABLE")
                          CtrlParam%AR_METHOD = ibset(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
                     case("DISABLE")
                           CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
              end select             
           end if   
 
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Warning: unknown keyword in &ARSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&METHOD")
                          !$$*** To get largest permitted number of neighbores
                           call Extract_Substr(STR,6,N,STRNUMB)
                           if(N .le. 0) then
                              write(*,fmt="(A)")      " MDPSCU Warning: no method is set for constructing activation region"
                              write(*,fmt="(A)")      '                 activation region to be disabled'
                              CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
                              call ONWARNING(gm_OnWarning)
                           else if(N .ge. 1) then
                              do I=1, N
                                 STRNUMB(I) = adjustl(STRNUMB(I))
                                 call UpCase(STRNUMB(I))
                                 select case(STRNUMB(I))
                                        case("USERPROC")
                                             CtrlParam%AR_METHOD = ior(CtrlParam%AR_METHOD, CP_USERDEF_AR)
                                        case("BYSEED")
                                             CtrlParam%AR_METHOD = ior(CtrlParam%AR_METHOD, CP_CENTPART_AR)
                                        case("BYEKIN")
                                             CtrlParam%AR_METHOD = ior(CtrlParam%AR_METHOD, CP_EKIN_AR)
                                        case("DISABLE")
                                             CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
                                        case("ENABLE")
                                             CtrlParam%AR_METHOD = ibset(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
                                        case("KEEPAR")
                                             CtrlParam%AR_METHOD = ior(CtrlParam%AR_METHOD,CP_KEEP_AR)
                                        case("BYNB")
                                             CtrlParam%AR_METHOD = ibset(CtrlParam%AR_METHOD, CP_BYNBSETBIT_AR)
                                        case("BYCELL","BYLC")
                                             CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD, CP_BYNBSETBIT_AR)
                                        case default
                                             write(*,fmt="(A)")      " MDPSCU Warning:'"//STRNUMB(I)(1:len_trim(STRNUMB(I)))//"' is unknown method for construct active region"
                                             write(*,fmt="(A)")      '                available methods include "USERPROC", "BYSEED", "BYEKIN" '
                                             write(*,fmt="(A, BZI6)")'                check control file at line:',LINE
                                             call ONWARNING(gm_OnWarning)
                                 end select
                              end do
                           end if
                           if(iand(CtrlParam%AR_METHOD, CP_USERDEF_AR) .eq. CP_USERDEF_AR  .and. &
                              iand(CtrlParam%AR_METHOD, CP_ENABLE_AR)  .eq. CP_ENABLE_AR   .and. &
                             (iand(CtrlParam%AR_METHOD,CP_CENTPART_AR) .eq. CP_CENTPART_AR .or.  &
                              iand(CtrlParam%AR_METHOD,CP_EKIN_AR).eq.CP_EKIN_AR) ) then
                              write(*,fmt="(A)")      ' MDPSCU Warning: method "USERPROC" can not be applied together with "BYSEED"and "BYEKIN" '
                              write(*,fmt="(A)")      '                 only method "USERPROC" to be applied '
                              write(*,fmt="(A, BZI6)")'                 check control file at line:',LINE
                              call ONWARNING(gm_OnWarning)
                              CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD,CP_ENABLEBIT_AR)
                           end if

                     case ("&BYSEED")
                          !$$*** To get type of atoms to be used as seed
                           call Extract_Numb(STR,min(mp_MXGROUP, size(STRNUMB)),N,STRNUMB)
                           do I=1, N
                              J = ISTR(STRNUMB(I))
                              CtrlParam%AR_CENTPART(J) = 1
                           end do

                     case ("&BYEKIN")
                          !$$*** To get type of atoms to be used as seed
                           call Extract_Numb(STR,1,N,STRNUMB)
                           if(N .ge. 1) then
                              CtrlParam%AR_EKIN = DRSTR(STRNUMB(1))
                           else
                              CtrlParam%AR_EKIN = 0
                           end if

                     case ("&EXTEND")
                          !$$*** To get extend for constructing active region
                           call Extract_Numb(STR,1,N,STRNUMB)
                           if(N.ge.1) CtrlParam%AR_EXTEND = ISTR(STRNUMB(1))
                           call Extract_Substr(STR,6,N,STRNUMB)
                           if(N .ge. 1) then
                             call UpCase(STRNUMB(1))
                             select case(STRNUMB(1))
                                    case("BYNB")
                                        CtrlParam%AR_METHOD = ibset(CtrlParam%AR_METHOD, CP_BYNBSETBIT_AR)
                                    case("BYCELL","BYLC")
                                       CtrlParam%AR_METHOD = ibclr(CtrlParam%AR_METHOD, CP_BYNBSETBIT_AR)
                                    case default
                                       write(*,fmt="(A)")      " MDPSCU Warning:'"//STRNUMB(1)(1:len_trim(STRNUMB(1)))//"' is method for construct active region"
                                       write(*,fmt="(A)")      '                available methods include "BYNB", "BYLC"'
                                       write(*,fmt="(A, BZI6)")'                check control file at line:',LINE
                                       call ONWARNING(gm_OnWarning)
                             end select
                           end if 

                     case ("&UPDATEFRE", "&UPDATESTP")
                         !$$*** To get how many time steps between updates of neighboure list
                         call Extract_Numb(STR,3,n,STRNUMB)
                         if(n.ge.3) then
                            CtrlParam%AR_UPTABMI = ISTR(STRNUMB(1))
                            CtrlParam%AR_UPTABMX = ISTR(STRNUMB(2))
                            CtrlParam%AR_DBITAB  = ISTR(STRNUMB(3))
                         else if(n.eq.2) then
                            CtrlParam%AR_UPTABMI = ISTR(STRNUMB(1))
                            CtrlParam%AR_UPTABMX = ISTR(STRNUMB(2))
                            CtrlParam%AR_DBITAB  = 100000
                         else if(n.eq.1) then
                            CtrlParam%AR_UPTABMI = ISTR(STRNUMB(1))
                            CtrlParam%AR_UPTABMX = CtrlParam%AR_UPTABMI
                            CtrlParam%AR_DBITAB  = 100000
                         end if
                         CtrlParam%AR_UPTAB = CtrlParam%AR_UPTABMI
              end select
           end do
           
           !--- check the complete and consistent of input
           if(iand(CtrlParam%AR_METHOD,CP_ENABLE_AR)  .eq. CP_ENABLE_AR .and. &
              iand(CtrlParam%AR_METHOD,CP_USERDEF_AR) .ne. CP_USERDEF_AR)  then
               if(iand(CtrlParam%AR_METHOD, CP_CENTPART_AR) .eq. CP_CENTPART_AR ) then
                  if(count(CtrlParam%AR_CENTPART .gt. 0) .le. 0) then
                     write(*,fmt="(A)")      " MDPSCU Warning: atom type for SEED of active region missed"
                     write(*,fmt="(A, BZI6)")'                 check control subsection &ACTIVEREGSUBCTL'
                    call ONWARNING(gm_OnWarning)
                  end if
               end if   
           end if

           return
  !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading control parameters for neigboring calculation."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_ActiveRegionCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_ForceCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, NGROUP, LINE)
  !***  PURPOSE:   to load the control parameter for cutoff control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE,NGROUP
     !--- local variables
      character*128::STRTMP(1)=""
      character*32::STRNUMB(mp_MXGROUP*mp_MXGROUP),KEYWORD
      integer::I, J, IJ, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU warning: unknown keyword in &POTENSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)


                    case ("&CUTOFF")
                          !$$*** To get the cutoff range for interaction
                          call Extract_Numb(STR,NGROUP*NGROUP,N,STRNUMB)
                          IJ=0
                          do I=1, NGROUP
                             do J=1, NGROUP
                                IJ = IJ + 1
                                if(IJ .GT. N) IJ = N
                                CtrlParam%RU(I,J) = DRSTR(STRNUMB(IJ))
                             end do
                          end do

                    case ("&TABLESIZE")
                         !$$*** To get the table size of potential
                          call Extract_Numb(STR,2,n,STRNUMB)
                          if(n.ge.1) CtrlParam%NUMFTABR = ISTR(STRNUMB(1))
                          if(n.ge.2) CtrlParam%NUMFTABE = ISTR(STRNUMB(2))
                          if(n.ge.3) CtrlParam%RHOSCAL  = DRSTR(STRNUMB(3))
                          call Extract_Substr(STR,1,n,STRTMP)
                          if(n.gt.0 ) then
                             STRTMP(1) = adjustl(STRTMP(1))
                             call UpCase(STRTMP(1))
                             if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "YES") then
                                CtrlParam%OUTPUTFTAB = 1
                             else
                                CtrlParam%OUTPUTFTAB = 0
                             end if
                          else
                             CtrlParam%OUTPUTFTAB = 0
                          end if

                          if(CtrlParam%NUMFTABR .lt. 1000 .or. CtrlParam%NUMFTABE .lt. 1000) then
                             write(*,fmt="A, BXI6") "MDPSCU warning: the force table size is smaller than 1000"
                             write(*,fmt="A, BZI6") "                accuracy for force calculation may loss"
                             write(*,fmt="A, BZI6") "                check control file at line: ", LINE
                             call ONWARNING(gm_OnWarning)
                          end if
              end select
           end do
           return
  !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading force cutoff control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_ForceCutoffCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_CommonParameter_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameters for all time sections
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in):: hFile
     type(SimMDCtrl)            :: CtrlParam
     character*(*)              :: STR
     integer                    :: LINE
     !--- local variables
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

          !**** to start load the controal parameters
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit

                     case default
                          write(*,*)"MDPSCU warning: unknown keyword in &COMMSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&BOXS")
                           !$$*** To get if multiple box to be used in the calculation
                           call Extract_Numb(STR,3,n,STRNUMB)
                           CtrlParam%MULTIBOX = ISTR(STRNUMB(1))
                           if( CtrlParam%MULTIBOX .LE. C_IZERO) CtrlParam%MULTIBOX = 1
                           if(N .GE. 2) CtrlParam%TOTALBOX = ISTR(STRNUMB(2))
                           if(CtrlParam%TOTALBOX .lt.CtrlParam%MULTIBOX) CtrlParam%TOTALBOX =CtrlParam%MULTIBOX
                           CtrlParam%TOTALBOX = (CtrlParam%TOTALBOX/CtrlParam%MULTIBOX)*CtrlParam%MULTIBOX
                           if(N .GE. 3) CtrlParam%INDEPBOX = ISTR(STRNUMB(3))

                     case ("&RANDSEED")
                          !$$*** To get seed for random number
                           call Extract_Numb(STR,4,n,STRNUMB)
                           do I=1, n
                              CtrlParam%SEED(I) = ISTR(STRNUMB(I))
                           end do

                end select
             end do

          return
   !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading force common control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
   end subroutine Load_CommonParameter_SimMDCtrl
 !****************************************************************************

 !****************************************************************************
  subroutine Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for LBFGS minimizer
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Error: unknown keyword in &LBFGSSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop

                     case ("&PGTOL")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%LBFGS_PGtol = 0.D0
                           if(n.ge. 1) then
                              CtrlParam%LBFGS_PGtol = DRSTR(STRNUMB(1))
                           end if

                           if(CtrlParam%LBFGS_PGtol .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &GFTOL input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for LBFGS_GFTol"
                              call ONWARNING(gm_OnWarning)
                           end if

                     case ("&FACTR")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%LBFGS_Factr = 1.D+1
                           if(n.ge. 1) then
                              CtrlParam%LBFGS_Factr = DRSTR(STRNUMB(1))
                           end if

                           if(CtrlParam%LBFGS_Factr .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &FACTR input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for LBFGS_Factr"
                              call ONWARNING(gm_OnWarning)
                           end if

                     case ("&MSAVE")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%LBFGS_Msave = 7
                           if(n.ge. 1) then
                              CtrlParam%LBFGS_Msave = ISTR(STRNUMB(1))
                           end if

                           if(CtrlParam%LBFGS_Msave .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &MSAVE input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for LBFGS_Msave"
                              call ONWARNING(gm_OnWarning)
                           end if

                     case ("&PGTOL_NEB")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_LBFGS_PGtol = 0.D0
                           if(n.ge. 1) then
                              CtrlParam%NEB_LBFGS_PGtol = DRSTR(STRNUMB(1))
                           end if

                           if(CtrlParam%NEB_LBFGS_PGtol .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &GFTOL_NEB input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for NEB_LBFGS_GFTol"
                              call ONWARNING(gm_OnWarning)
                           end if

                     case ("&FACTR_NEB")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_LBFGS_Factr = 1.D+1
                           if(n.ge. 1) then
                              CtrlParam%NEB_LBFGS_Factr = DRSTR(STRNUMB(1))
                           end if

                           if(CtrlParam%NEB_LBFGS_Factr .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &FACTR input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for NEB_LBFGS_Factr"
                              call ONWARNING(gm_OnWarning)
                           end if

                     case ("&MSAVE_NEB")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_LBFGS_Msave = 7
                           if(n.ge. 1) then
                              CtrlParam%NEB_LBFGS_Msave = ISTR(STRNUMB(1))
                           end if

                           if(CtrlParam%NEB_LBFGS_Msave .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &MSAVE input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for LBFGS_Msave"
                              call ONWARNING(gm_OnWarning)
                           end if
              end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading LBFGS control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_LBFGSCtl_SimMDCtrl
  !****************************************************************************

 !****************************************************************************
  subroutine Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for steepest-descant minimizer
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(10),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Error: unknown keyword in &STEEPSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop
               
                     case ("&ALPHA")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           if(n.ge. 1) then
                              CtrlParam%STEEPEST_Alpha = DRSTR(STRNUMB(1))
                           end if

                           if(CtrlParam%STEEPEST_Alpha .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &ALPHA input at line", LINE
                              write(*,fmt="(A, BZI6)") "                 default value to be set for STEEPEST_Alpha"
                              call ONWARNING(gm_OnWarning)
                              stop
                           end if

                     case ("&STEPCOND", "&STEPBOUND")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,2,n,STRNUMB)
                           if(n.eq. 1) then
                              CtrlParam%STEEPEST_MxStep = DRSTR(STRNUMB(1))
                           else if(n .ge. 2) then
                              CtrlParam%STEEPEST_MiStep = min(DRSTR(STRNUMB(1)), DRSTR(STRNUMB(2)))
                              CtrlParam%STEEPEST_MxStep = max(DRSTR(STRNUMB(1)), DRSTR(STRNUMB(2)))
                           end if

                           if(CtrlParam%STEEPEST_MxStep .lt. 0 .or. n.lt.1) then
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &STEPBOUND input at line", LINE
                              stop
                           end if

                     case ("&DELTAPOT", "&POTCOND")
                           !$$*** To get if periodic conditions will be used
                           call Extract_Numb(STR,1,n,STRNUMB)
                           if(n.ge. 1) then
                              CtrlParam%STEEPEST_MiDelE = DRSTR(STRNUMB(1))
                           else 
                              write(*,fmt="(A, BZ16)") " MDPSCU Warning: wrong &DELTAPOT input at line", LINE
                              stop
                           end if
              end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading STEEPEST control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_STEEPESTCtl_SimMDCtrl
  !****************************************************************************

 !****************************************************************************
  subroutine Load_CGCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for CG minimizer
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE
     !--- local variables
     !----
            call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  end subroutine Load_CGCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_STCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for box control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in)::hFile
     type(SimMDCtrl)            ::CtrlParam
     character*(*)              ::STR
     integer                    ::LINE
     !--- local variables

           call Load_STCtrlParam(hFile, CtrlParam%ST_CTRL, STR, LINE)
  end subroutine Load_STCtl_SimMDCtrl
  !****************************************************************************

 !****************************************************************************
  subroutine Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameter for steepest-descant minimizer
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,          intent(in)::hFile
     type(SimMDCtrl)             ::CtrlParam
     character*(*)               ::STR
     integer                     ::LINE
     !--- local variables
      character*256::STRTMP(1)=""
      character*32::STRNUMB(4),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MDPSCU Error: unknown keyword in &EVENTDETECTSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop
               
                     case( "&DRTOL")
                           !$$--- To get if the tolerance displacement for event detection
                           call Extract_Numb(STR,1,n,STRNUMB)
                           if(N .LT. 1) then
                              write(*,fmt="(' MDPSCU Warning: parameter controlling the structure change detetction is missed')")
                              write(*,fmt="('                 check control file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage:  &DRTOL value ')")
                              call ONWARNING(gm_OnWarning)
                           else
                              CtrlParam%STRCUT_DRTol = DRSTR(STRNUMB(1))
                           end if
                           
                     case("&LBFGSSUBCTL")
                           call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                     case("&STEEPESTSUBCTL", "&CGSUBCTL")
                          call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  
                     case( "&QUENCHSTEP")
                           !$$*** queching method and step for event checking
                            call Extract_Numb(STR,1,n,STRNUMB)
                            if(n .lt. 1) then
                               write(*,fmt="(' MDPSCU Error: the time steps for quenching should be set')")
                               write(*,fmt="('               check control file at line:', BZI6)") LINE
                               write(*,fmt="(' Usage: &QUICKDAMP (&QUENCH) steps need for damping =')")
                               write(*,fmt="(' Process to be stopped')")
                               stop
                            end if
                            CtrlParam%Quench_Steps = ISTR(STRNUMB(1))
  
                            call Extract_Substr(STR,1,N,STRTMP)
                            if(n .ge. 1) then
                               call UpCase(STRTMP(1))
                               if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "LBFGS") then
                                  CtrlParam%Quench_Meth = CP_DAMPSCHEME_LBFGS
                               else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "CG") then
                                   CtrlParam%Quench_Meth = CP_DAMPSCHEME_CG
                               else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "CG-LS") then
                                   CtrlParam%Quench_Meth = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                 
                               else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "ST") then
                                  CtrlParam%Quench_Meth = CP_DAMPSCHEME_ST
                               else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "ST-LS") then
                                  CtrlParam%Quench_Meth = ior(CP_DAMPSCHEME_ST, CP_DAMPSCHEME_LSEARCH)                                     
                               else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "DYN" .or. &
                                       STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "DYNAMICS" ) then
                                  CtrlParam%Quench_Meth = CP_DAMPSCHEME_DYN
                               else
                                  write(*,fmt="(A)")      " MDPSCU Error: the damping scheme "//STRTMP(1)(1:len_trim(STRTMP(1)))//" is unknown"
                                  write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                  write(*,fmt="(A)")      ' Process to be stopped'
                                  stop
                               end if
                            end if
                end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading EventDetect control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_EventDetectCtl_SimMDCtrl
  !****************************************************************************  

 !****************************************************************************
  subroutine Load_AnalyCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameters for all time sections
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,          intent(in) ::hFile
     type(SimMDCtrl)              ::CtrlParam
     character*(*)                ::STR
     integer                      ::LINE
     !--- local variables
      character*32::STRNUMB(10),KEYWORD
      integer::I, N, NEEDDAMP, DAMPSCHEME,NEEDNEWATOMS
      integer(2)::NEWATOMS(2)
      equivalence(NEWATOMS(1), NEEDNEWATOMS)
     !----

          !**** to start load the controal parameters
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit

                     case default
                          write(*,*)"MDPSCU warning: unknown keyword in &ANALYSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&JOBSEL")
                           !$$*** To get job range to be analysis
                           call Extract_Numb(STR,3,n,STRNUMB)
                           if(N.gt. 0) then
                              CtrlParam%JOBID0 = ISTR(STRNUMB(1))
                              if(N .GE. 2) CtrlParam%JOBID1  = ISTR(STRNUMB(2))
                              if(N .GE. 3) CtrlParam%JOBIDSTEP = ISTR(STRNUMB(3))
                           end if
                     case ("&CFGSEL")
                           !$$*** To get cfg range to be analysis
                           call Extract_Numb(STR,3,n,STRNUMB)
                           if(N.gt.0) then
                              CtrlParam%STARTCFG = ISTR(STRNUMB(1))
                              if(N .GE. 2) CtrlParam%ENDCFG  = ISTR(STRNUMB(2))
                              if(N .GE. 3) CtrlParam%CFGSTEP = ISTR(STRNUMB(3))
                           end if

                     case ("&BOXSEL")
                           !$$*** To get cfg range to be analysis
                           call Extract_Numb(STR,3,n,STRNUMB)
                           if(N.gt.0) then
                              CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                              if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                              if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))
                           end if

                     case ("&NEWATOMS")
                           !*** get the controal paraemter of output average CP values
                           call Extract_Numb(STR,2, N, STRNUMB)
                           if(N.gt.0) then
                              NEWATOMS(1) = ISTR(STRNUMB(1))
                              NEWATOMS(2) = 1
                              if(N.ge.2) NEWATOMS(2)  = ISTR(STRNUMB(2))
                              CtrlParam%NEEDNEWATOMS  = NEEDNEWATOMS
                           end if

                     case ("&QUICKDUMP", "&QUICKDAMP")
                           !*** get the controal paraemter of output average CP values
                           NEEDDAMP   = CtrlParam%NEEDDAMP
                           DAMPSCHEME = CtrlParam%NEEDDAMPTYPE
                           call Extract_Numb(STR,1, N, STRNUMB)
                           if(N.gt.0) NEEDDAMP = ISTR(STRNUMB(1))

                          call Extract_Substr(STR,1,N,STRNUMB)
                          if(N .ge. 1) then
                             call UpCase(STRNUMB(1))
                             if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "LBFGS") then
                                DAMPSCHEME = CP_DAMPSCHEME_LBFGS
                             else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYN" .or. &
                                STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYNAMICS" ) then
                                DAMPSCHEME = CP_DAMPSCHEME_DYN
                             else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST" ) then
                                 DAMPSCHEME = CP_DAMPSCHEME_ST
                             else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST-LS") then
                                 DAMPSCHEME = ior(CP_DAMPSCHEME_ST, CP_DAMPSCHEME_LSEARCH)                                     
                             else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG" ) then
                                 DAMPSCHEME = CP_DAMPSCHEME_CG 
                             else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG-LS") then
                                 DAMPSCHEME = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                                                
                             else
                                write(*,fmt="(A)")      " MDPSCU Error: the damping scheme "//STRNUMB(1)(1:len_trim(STRNUMB(1)))//" is unknown"
                                write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                write(*,fmt="(A)")      ' Process to be stopped'
                                stop
                             end if
                          end if
                          call SetNeedDamp_SimMDCtrl(CtrlParam, NEEDDAMP, DAMPSCHEME)
                end select
             end do

          return
   !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading analysis control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_AnalyCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine LoadExInput_SimMDCtrl(Tag, Fname, CtrlParam)
  !***  PURPOSE:   to load the external control parameters from a file,
  !                the external parameters to be used in analysis application
  !
  !     INPUT:
  !                Tag,   the Tag for the input
  !                Fname, the filename
  !
  !     OUTPUT     CtrlParam, with the member ExInputStat inputed
  implicit none
      !--- dummy varioables
      character*(*), intent(in) ::Fname
      character*(*), intent(in) ::Tag
      type(SimMDCtrl)           ::CtrlParam
      !--- local variables
      integer::I, L, IN

          !$$--- to check if the statments have been loaded
          IN = 0
          do I=1, size(CtrlParam%ExInputStat)
             L = len_trim(CtrlParam%ExInputStat(I)%filename)
             if(CtrlParam%ExInputStat(I)%filename(1:L) .eq. fname(1:len_trim(fname))) then
                IN = I
                exit
             end if
          end do
          !$$--- if not assigned an statement structure to the external input
          if(IN .le. 0) then
             do I=1, size(CtrlParam%ExInputStat)
                L = len_trim(CtrlParam%ExInputStat(I)%filename)
                if(L.le.0) then
                   IN = I
                   exit
                end if
             end do
          end if

          call Load_InputStatements(fname, CtrlParam%ExInputStat(IN), Tag)
          return
  end subroutine LoadExInput_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine GetExInput_SimMDCtrl(CtrlParam, Tag, Input, Statu)
  !***  PURPOSE:   to get a copy of the statements for given Tag
  !
  !     INPUT:
  !                Tag,       the Tag for the input
  !                CtrlParam, control parameters
  !
  !     OUTPUT     Input,     the statementlist for Tag
  !                Statu      =0, connot find the statemnets for Tag
  !                           =1, get statments
  implicit none
      !--- dummy varioables
      character*(*),          intent(in)::Tag
      type(SimMDCtrl),        intent(in)::CtrlParam
      type(InputStatements)             ::Input
      integer                           ::Statu
      !--- local variables
      integer::I, L

          Statu = 0
          do I=1, size(CtrlParam%ExInputStat)
             L = len_trim(CtrlParam%ExInputStat(I)%Stag)
             if(CtrlParam%ExInputStat(I)%Stag(1:L) .eq. Tag(1:len_trim(Tag))) then
                call Copy_InputStatements(CtrlParam%ExInputStat(I), Input)
                Statu = 1
                exit
             end if
          end do
          return
  end subroutine GetExInput_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine ExtractAnalyComm_SimMDCtrl(CtrlParam)
  !***  PURPOSE:   to extract the common control parameters for analysis
  !                routines
  !
  !     INPUT:     CtrlParam, control parameters
  !
  !     OUTPUT     CtrlParam, control parameter with parameters for analysis updated
  !
   implicit none
      !--- dummy varioables
      type(SimMDCtrl)::CtrlParam
      !--- local variables
      integer::I, N, LINE, NEEDDAMP, DAMPSCHEME
      type(InputStatements)::INPUTS
      character*256::STR
      character*32::STRNUMB(3)

          !***  NOTE: the control parameters of previous inputs will be overwriten
          !           by the later one
          do I=1, size(CtrlParam%ExInputStat)
             if(.not.associated (CtrlParam%ExInputStat(I)%Stat) ) cycle

                call Copy_InputStatements(CtrlParam%ExInputStat(I), INPUTS)

                !*** get the control parameters for dummping
                 call Get_InputStatements("&QUICKDAMP", INPUTS, STR, LINE)
                 if(LINE .eq. 0) call Get_InputStatements("&QUICKDAMP", INPUTS, STR, LINE)
                 if(LINE .gt. 0) then
                    call Extract_Numb(STR,1, N, STRNUMB)
                    NEEDDAMP = ISTR(STRNUMB(1))
                    call Extract_Substr(STR,1,N,STRNUMB)
                    if(N .ge. 1) then
                       call UpCase(STRNUMB(1))
                       if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "LBFGS") then
                          DAMPSCHEME = CP_DAMPSCHEME_LBFGS
                       else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYN" .or. &
                          STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYNAMICS" ) then
                          DAMPSCHEME = CP_DAMPSCHEME_DYN
                       else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG") then
                           DAMPSCHEME = CP_DAMPSCHEME_CG
                       else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG-LS") then
                           DAMPSCHEME = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                 
                       else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST") then
                           DAMPSCHEME = CP_DAMPSCHEME_ST
                       else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST-LS") then
                           DAMPSCHEME = ior(CP_DAMPSCHEME_ST, CP_DAMPSCHEME_LSEARCH)                                     
                       else
                          write(*,fmt="(A)")         " MDPSCU Error: the damping scheme "//STRNUMB(1)(1:len_trim(STRNUMB(1)))//" is unknown"
                          write(*,fmt="(A, BZI6, A)")'               check control file at line:', LINE, ' in '// &
                                                  INPUTS%filename(1:len_trim(INPUTS%filename))
                          write(*,fmt="(A)")        ' Process to be stopped'
                          stop
                       end if
                    end if
                    call SetNeedDamp_SimMDCtrl(CtrlParam, NEEDDAMP, DAMPSCHEME)
                 end if

                 !*** To get job range to be analysis
                 call Get_InputStatements("&JOBSEL", INPUTS, STR, LINE)
                 if(LINE .gt. 0) then
                    call Extract_Numb(STR,3,n,STRNUMB)
                    CtrlParam%JOBID0 = ISTR(STRNUMB(1))
                    if(N .GE. 2) CtrlParam%JOBID1  = ISTR(STRNUMB(2))
                    if(N .GE. 3) CtrlParam%JOBIDSTEP = ISTR(STRNUMB(3))
                 end if

                 !$$*** To get cfg range to be analysis
                 call Get_InputStatements("&CFGSEL", INPUTS, STR, LINE)
                 if(LINE .gt. 0) then
                    call Extract_Numb(STR,3,n,STRNUMB)
                    CtrlParam%STARTCFG = ISTR(STRNUMB(1))
                    if(N .GE. 2) CtrlParam%ENDCFG  = ISTR(STRNUMB(2))
                    if(N .GE. 3) CtrlParam%CFGSTEP = ISTR(STRNUMB(3))
                 end if

                 !$$*** To get cfg range to be analysis
                 call Get_InputStatements("&BOXSEL", INPUTS, STR, LINE)
                 if(LINE .gt. 0) then
                    call Extract_Numb(STR,3,n,STRNUMB)
                    CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                    if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                     if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))
                 end if
                 call Release_InputStatements(INPUTS)
          end do
          return
  end subroutine ExtractAnalyComm_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Parameter_SimMDCtrl_OLD(hFile, CtrlParam, SimBox)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(SimMDCtrl)    ::CtrlParam
     type(SimMDBox)     ::SimBox

     !--- local variables
      character*256::STR,STRTMP(1)=""
      character*32::STRNUMB(50)
      integer::I,J, IJ, N, LINE
     !----
         LINE = 0
         !$$*** To get the temprature
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,1,n,STRNUMB)
              CtrlParam%TI = DRSTR(STRNUMB(1))
         !$$*** To get the Nose parameter for temprature
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,1,n,STRNUMB)
              CtrlParam%NoseQ = DRSTR(STRNUMB(1))
         !$$*** To get the pressure
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,1,n,STRNUMB)
              CtrlParam%PEX = DRSTR(STRNUMB(1))
              CtrlParam%PEXTENSOR   = C_ZERO
              CtrlParam%PEXTENSOR(1,1) = CtrlParam%PEX
              CtrlParam%PEXTENSOR(2,2) = CtrlParam%PEX
              CtrlParam%PEXTENSOR(3,3) = CtrlParam%PEX

         !$$*** To get the Jonhson parameter for pressure
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,1,n,STRNUMB)
              CtrlParam%JONSHONW = DRSTR(STRNUMB(1))

        !$$*** To get total integral steps
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,2,n,STRNUMB)
             CtrlParam%ITE = ISTR(STRNUMB(1))
             if(N.lt. 2) then
               if( CtrlParam%ITE.lt.0) then
                write(*,*) " Error in reading control parameters for time"
                write(*,*) " terminal time (in ps) is needed"
                stop
               end if
                CtrlParam%TEMTM = 1.D60  !to give an large value
             else
               CtrlParam%TEMTM = DRSTR(STRNUMB(2))
             end if
             CtrlParam%IT1 = (CtrlParam%IT0-1)+CtrlParam%ITE

        !$$*** To get the length of time step (in fs)
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,4,n,STRNUMB)
             if(N .LT. 4) then
                write(*,*) " Error in reading control parameters for time step"
                write(*,*) " 4 parameters are expected"
                stop
             end if
             CtrlParam%IHDUP = ISTR(STRNUMB(1))
             CtrlParam%HMI = DRSTR(STRNUMB(2))
             CtrlParam%HMX = DRSTR(STRNUMB(3))
             CtrlParam%DMX = DRSTR(STRNUMB(4))
             if(CtrlParam%DMX .LT. 1.D-6) CtrlParam%DMX = 1.D-6
             !$$The time step start from its minimum value
             CtrlParam%H =  CtrlParam%HMI
        !$$*** To get how many times of thermalization will be carried out
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,1,n,STRNUMB)
             CtrlParam%IVTIME = ISTR(STRNUMB(1))
        !$$*** To get how many time steps between two thermalizations
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,2,n,STRNUMB)
             CtrlParam%IVPAS = ISTR(STRNUMB(1))
             if(n .ge. 2) then
                CtrlParam%IVTIME0 = ISTR(STRNUMB(2))
             end if
        !$$*** To get if damping is required
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,1,n,STRNUMB)
             CtrlParam%DAMPTIME1 = ISTR(STRNUMB(1))

        !$$*** To get if electron-phono coupling is required
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,1,n,STRNUMB)
             CtrlParam%EPCTIME1 = ISTR(STRNUMB(1))

            if(CtrlParam%EPCTIME1 .gt. 0) then
               !$$--- get the filename storing EPC coupling
               call GetInputStrLine(hFile,STR, LINE, "!", *100)
               call Extract_Substr(STR,1,n,STRTMP)
               !CtrlParam%EPCPHYS = STRTMP(1)
               CtrlParam%EPCTIME1 = 0
            end if

        !$$*** To get if Berendsen pressure couping is used
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,1,n,STRNUMB)
             CtrlParam%IFBPC = ISTR(STRNUMB(1))

             if(CtrlParam%IFBPC) then
               !$$--- get characteristic time of coupling
               call GetInputStrLine(hFile,STR, LINE, "!", *100)
               call Extract_Substr(STR,1,n,STRTMP)
               CtrlParam%BPCDATA = STRTMP(1)
             end if

        !$$*** To get if periodic conditions will be used
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,3,n,STRNUMB)
             CtrlParam%IFPD(1) = ISTR(STRNUMB(1))
             CtrlParam%IFPD(2) = ISTR(STRNUMB(2))
             CtrlParam%IFPD(3) = ISTR(STRNUMB(3))

       !$$*** To get if multiple box to be used in the calculation
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,3,n,STRNUMB)
             CtrlParam%MULTIBOX = ISTR(STRNUMB(1))
             if( CtrlParam%MULTIBOX .LE. C_IZERO) CtrlParam%MULTIBOX = 1

             if(N .GE. 2) CtrlParam%TOTALBOX = ISTR(STRNUMB(2))
             if(CtrlParam%TOTALBOX .lt.CtrlParam%MULTIBOX) CtrlParam%TOTALBOX =CtrlParam%MULTIBOX
             CtrlParam%TOTALBOX = (CtrlParam%TOTALBOX/CtrlParam%MULTIBOX)*CtrlParam%MULTIBOX

             if(N .GE. 3) CtrlParam%INDEPBOX = ISTR(STRNUMB(3))


       !$$*** To get the cell numbers along x,y,z
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,3,n,STRNUMB)
             !CtrlParam%NC(1) = ISTR(STRNUMB(1))
             !CtrlParam%NC(2) = ISTR(STRNUMB(2))
             !CtrlParam%NC(3) = ISTR(STRNUMB(3))
       !$$*** To get largest permitted number of neighbores
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,1,n,STRNUMB)
             CtrlParam%NB_MXNBS = ISTR(STRNUMB(1))
       !$$*** To get how many time steps between updates of neighboure list
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,3,n,STRNUMB)
             if(n.ge.3) then
                CtrlParam%NB_UPTABMI = ISTR(STRNUMB(1))
                CtrlParam%NB_UPTABMX = ISTR(STRNUMB(2))
                CtrlParam%NB_DBITAB  = ISTR(STRNUMB(3))
             else
                CtrlParam%NB_UPTABMI = ISTR(STRNUMB(1))
                CtrlParam%NB_UPTABMX = CtrlParam%NB_UPTABMI
                CtrlParam%NB_DBITAB  = 100000
             end if
             CtrlParam%NB_UPTAB = CtrlParam%NB_UPTABMI
       !$$*** To get cutoff range for neighbores
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,SimBox%NGROUP*SimBox%NGROUP,N,STRNUMB)
             IJ=0
             DO I=1, SimBox%NGROUP
                DO J=1, SimBox%NGROUP
                   IJ = IJ + 1
                   if(IJ .GT. N) IJ = N
                   CtrlParam%NB_RM(I,J) = DRSTR(STRNUMB(IJ))
                END DO
             END DO

       !$$*** To get the cutoff range for interaction
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,SimBox%NGROUP*SimBox%NGROUP,N,STRNUMB)

             IJ=0
             DO I=1, SimBox%NGROUP
                DO J=1, SimBox%NGROUP
                   IJ = IJ + 1
                   if(IJ .GT. N) IJ = N
                   CtrlParam%RU(I,J) = DRSTR(STRNUMB(IJ))
                END DO
             END DO

      !$$*** To get the table size of potential
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,2,n,STRNUMB)
             if(N.ge.1) CtrlParam%NUMFTABR = ISTR(STRNUMB(1))
             if(N.ge.2) CtrlParam%NUMFTABE = ISTR(STRNUMB(2))
             if(N.ge.3) CtrlParam%RHOSCAL  = ISTR(STRNUMB(3))

      !$$*** To get seed for random number
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             call Extract_Numb(STR,4,n,STRNUMB)
             DO I=1, n
                CtrlParam%SEED(I) = ISTR(STRNUMB(I))
             END DO

      !$$*** To get seed for random number
            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call Extract_Numb(STR,1,n,STRNUMB)
            CtrlParam%RESTART = ISTR(STRNUMB(1))

      !$$*** To get output contral
          !$$--- To get for how many time steps to print out temperature etc
           call GetInputStrLine(hFile,STR, LINE, "!", *100)
           call Extract_Numb(STR,1,n,STRNUMB)
           CtrlParam%TimestpQ = ISTR(STRNUMB(1))

          !$$--- To get for how many time steps to print out configuration
           call GetInputStrLine(hFile,STR, LINE, "!", *100)
           call Extract_Numb(STR,2,n,STRNUMB)
           CtrlParam%TimestpG = ISTR(STRNUMB(1))
           if(n.ge.2) then
             CtrlParam%TimestpR = ISTR(STRNUMB(2))
           else
             CtrlParam%TimestpR = 0
           end if


         !$$--- To get for how many time steps to save current statu
          call GetInputStrLine(hFile,STR, LINE, "!", *100)
          call Extract_Numb(STR,2,n,STRNUMB)
          CtrlParam%TimestpSave= ISTR(STRNUMB(1))
          if(CtrlParam%TimestpSave .eq. 0)CtrlParam%TimestpSave = CtrlParam%ITE

         !$$--- the output filenames are left to
         !      MD_Globle_Variables. ref.
         !      Initialize_Globle_Variables


  200    return
  !---------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading control parameters."
         write(*,*)"The process to be stopped."
         stop
  end subroutine Load_Parameter_SimMDCtrl_OLD
  !****************************************************************************



  end module MD_TYPEDEF_SimMDCtrl
