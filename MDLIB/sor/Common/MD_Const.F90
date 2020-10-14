  module MD_CONSTANTS
  !***  DESCRIPTION: this module is to define the constants to be used
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010
  ! 
  use MSM_CONSTANTS
  implicit none

  !*** type define
      integer, parameter::CP_TIMECYCLE_END      = CP_HBIT
      integer, parameter::CP_TIMECYCLE_CONTINUE = 0

      integer, parameter::CP_FILE_NOTEXIST      = 0
      integer, parameter::CP_FILE_EXIST         = 1

      !--- parameters controling the process
      integer, parameter::CP_WARNING_CONTIN   = 0
      integer, parameter::CP_WARNING_PAUSE    = 1
      integer, parameter::CP_WARNING_STOP     = 2
      integer, parameter::CP_DISABLE          = 0

      !--- parameters denote methods of controling temperature of box
      integer, parameter::CP_DAMPSCHEME_LBFGS   = 0
      integer, parameter::CP_DAMPSCHEME_DYN     = 1
      integer, parameter::CP_DAMPSCHEME_CG      = 2
      integer, parameter::CP_DAMPSCHEME_ST      = 3
      integer, parameter::CP_DAMPSCHEME_LSEARCH = 2**16

      integer, parameter::CP_THERMALSCHEME_MC    = 0
      integer, parameter::CP_THERMALSCHEME_VSCAL = 1
      integer, parameter::CP_THERMALSCHEME_PSCAL = 2
      integer, parameter::CP_THERMALSCHEME_ASCAL = 3

      integer, parameter::CP_TICTRL_NONE        = 0
      integer, parameter::CP_TICTRL_GLOBAL      = 1
      integer, parameter::CP_TICTRL_LOCAL       = 2
      integer, parameter::CP_TICTRL_METH_EPC    = 1
      integer, parameter::CP_TICTRL_METH_ST     = 2
      integer, parameter::CP_TICTRL_METH_EPCST  = CP_TICTRL_METH_EPC+CP_TICTRL_METH_ST
      integer, parameter::CP_TICTRL_METH_ASCAL  = 4
      integer, parameter::CP_EPCCTRL_MOD_MAN    = 1
      integer, parameter::CP_EPCCTRL_MOD_SMFT   = 2
      integer, parameter::CP_STCTRL_MOD_LS_B    = 1
      integer, parameter::CP_STCTRL_MOD_LS_Z85  = 2
      integer, parameter::CP_STCTRL_MOD_LS_Z95  = 3
      integer, parameter::CP_STCTRL_MOD_USER    = 4
      integer, parameter::CP_STCTRL_MOD_SRIM    = 5
      integer, parameter::CP_STCTRL_MOD_OR      = 2**16
      integer, parameter::CP_STCTRL_MOD_DEN_G   = 1
      integer, parameter::CP_STCTRL_MOD_DEN_L   = 2

      !--- parameters denote methods of construct activation regionm
      integer, parameter::CP_ENABLEBIT_AR        = 16
      integer, parameter::CP_ENABLE_AR           = 2**(CP_ENABLEBIT_AR)
      integer, parameter::CP_CENTPART_AR         = 1
      integer, parameter::CP_EKIN_AR             = 2
      integer, parameter::CP_BYNBSETBIT_AR       = 27
      integer, parameter::CP_BYNB_AR             = 2**(CP_BYNBSETBIT_AR)
      integer, parameter::CP_KEEP_AR             = 2**24
      integer, parameter::CP_USERDEF_AR          = 2**30

      !--- parameters denote methods of boost force
      integer, parameter::CP_DISABLE_BOOST       = 0
      integer, parameter::CP_SPRING_BOOST        = 1
      integer, parameter::CP_GAUSS_BOOST         = 2
      integer, parameter::CP_BOND_BOOST          = 4

      !--- parameters denote method for activation system in ART simulations
      integer, parameter::CP_CENTPARTACT_ART     = 1
      integer, parameter::CP_MULTPLEACT_ART      = 2
      integer, parameter::CP_AR_REGION_ART       = 2**2
      integer, parameter::CP_USERDEFACT_ART      = 2**30


      !--- parameters denote the state of simulation box ----
      !
      !--- max number of atom types
       integer, parameter::mp_MXGROUP             = 16
       integer, parameter::CP_STATU_ACTIVE_BITPOS = 0
       integer, parameter::CP_STATU_DEACTIVE = 0                             ! NOTE: active flag changed from 0 to 1 by HOU Qing, 2011-03-11
       integer, parameter::CP_STATU_ACTIVE   = 1
       integer, parameter::CP_STATU_FIXPOSX  = 2**1
       integer, parameter::CP_STATU_FIXPOSY  = 2**2
       integer, parameter::CP_STATU_FIXPOSZ  = 2**3
       integer, parameter::CP_STATU_FIXPOS   = CP_STATU_FIXPOSX+CP_STATU_FIXPOSY+CP_STATU_FIXPOSZ
       integer, parameter::CP_STATU_FIXVELX  = 2**4
       integer, parameter::CP_STATU_FIXVELY  = 2**5
       integer, parameter::CP_STATU_FIXVELZ  = 2**6
       integer, parameter::CP_STATU_FREEPART = 2**8                           !-- if an atom is maked as free particle, the moving direction of the particle
                                                                              !   will be changed randomly when it pass through the boundary of the simulation box
       integer, parameter::CP_STATU_ARSEEDPART = 2**9                         !-- if an atom is maked as CP_STATU_ARSEEDPART particle, it can be used for creating
                                                                              !   activation region, see also CP_CENTPART_AR
       integer, parameter::CP_STATU_HASEPC     = 2**10                        !-- if an atom is maked as CP_STATU_HASEPC, the EPC calculation will be performed for
                                                                              !   the atom

       integer, parameter::CP_STATU_OUTOFBOX  = 2**16
       integer, parameter::CP_STATU_REFLECT   = 2*CP_STATU_OUTOFBOX
       integer, parameter::CP_STATU_TRANSMIT  = 4*CP_STATU_OUTOFBOX
       integer, parameter::CP_STATU_PASSBOUND = 8*CP_STATU_OUTOFBOX           !indicate the particle has pass through the boudary

       integer, parameter::CP_STATU_HWORD    = CP_HWORD               !0xFFFE0000
       integer, parameter::CP_STATU_LWORD    = CP_LWORD               !0x0000FFFF

      !---KEYWORDS LIST
       character*9, parameter::PKW_OUTCFG_FORMAT14  = "&BOXCFG14"
       character*9, parameter::PKW_OUTCFG_FORMAT15  = "&BOXCFG15"
       character*9, parameter::PKW_OUTCFG_FORMAT16  = "&BOXCFG16"
       character*9, parameter::PKW_OUTCFG_FORMAT18  = "&BOXCFG18"

       character*7, parameter::PKW_OUTCFG_FORMATXYZ = "&CFGXYZ"
       character*7, parameter::PKW_INCFG_FORMATXYZ  = PKW_OUTCFG_FORMATXYZ

      !--- KEYWORDS may appear in CFGXYZ files
       character*5,  parameter::PKW_LATT = "&LATT"
       character*10, parameter::PKW_LENTHUNIT = "&LENTHUNIT"

      !--- KEYWORDs for type of potentials
       character(len=8), parameter::PKW_POTTYPELIST(2) = (/"FS_TYPE", "EAM_TYPE"/)

      !--- KEYWORDs for type of applications
       character(len=8), parameter::PKW_APPTYPELIST(4) = (/"GMD",            &  !generic MD
                                                           "PARREP",         &  !parallel replica
                                                           "NEB",            &  !neb
                                                           "ART"             &  !art
                                                          /)  

      !--- some global variables
       integer::        gm_NMPI = 1                       ! the number of MPI processes
       integer::        gm_hFILELOG = 0                   ! the IO unit of log file
       character*256::  gm_ExeName = ""                   ! the name of excutable
       integer::        gm_OnWarning = CP_WARNING_PAUSE   ! the flag indicating what action will be taken on wraning message
       character*8::    gm_AppType = "GMD"                ! the application type, one in the keword list PKW_APPTYPELIST


       character*256::  gm_cfgFileName(4) = ""            ! the filenames of configure files used when non-"setupfile" is used for
                                                          ! analysis routines. ref.MD_SimBoxArray_ToolShell_14_GPU.F90 and
                                                          ! MD_SimBoxArray_ToolShell_14_CPU.F90
       character*256::  gm_ctrlFileName(4) = ""             ! the filenames of control files used when non-"setupfile" is used for
                                                          ! analysis routines. ref.MD_SimBoxArray_ToolShell_14_GPU.F90 and
                                                          ! MD_SimBoxArray_ToolShell_14_CPU.F90
       character*256::  gm_outFileName(4) = ""              ! the filenames of input files used when non-"setupfile" is used for
                                                          ! analysis routines. ref.MD_SimBoxArray_ToolShell_14_GPU.F90 and
                                                          ! MD_SimBoxArray_ToolShell_14_CPU.F90

  end module MD_CONSTANTS

