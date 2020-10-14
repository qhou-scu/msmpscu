  module MC_Randwalk_Const
  !***  DESCRIPTION: this module is to define the constants to be used
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010
  ! 
  use MSM_CONSTANTS
  implicit none

  !*** type define
  integer,   parameter::CP_MX_RADSTYLE       =  2
  integer,   parameter::CP_MX_PREEXITSTYLE   =  5
 
  !*** type of walk style
  integer(KINDINT),   parameter::CP_WALKTYPE_NOWALK        = 0
  integer(KINDINT),   parameter::CP_WALKTYPE_WTD_EXP       = 1
  integer(KINDINT),   parameter::CP_WALKTYPE_WTD_POW       = 2
  integer(KINDINT),   parameter::CP_WALKTYPE_WTD_MASK      = 15           !1111

  integer(KINDINT),   parameter::CP_WALKTYPE_DIM_1D        = 16
  integer(KINDINT),   parameter::CP_WALKTYPE_DIM_2D        = CP_WALKTYPE_DIM_1D*2
  integer(KINDINT),   parameter::CP_WALKTYPE_DIM_3D        = CP_WALKTYPE_DIM_2D*2
  integer(KINDINT),   parameter::CP_WALKTYPE_DIM_MASK      = 255 - CP_WALKTYPE_WTD_MASK
  integer(KINDINT),   parameter::CP_WALKTYPE_MAXJPPATH     = 8           ! max discrete jumping path
  integer(KINDINT),   parameter::CP_WALKTYPE_MAXJVECT      = 4           ! max possible jumping vectors, these jumps haave the same probabbility

  integer(KINDINT),   parameter::CP_WALKTYPE_DIS_FIXSTEP   = 2**17
  integer(KINDINT),   parameter::CP_WALKTYPE_DIS_POWDIS    = CP_WALKTYPE_DIS_FIXSTEP*2  

  !*** state of walker
  integer(KINDINT),   parameter::CP_WALKSTAT_NOTACTIVE     = 0
  integer(KINDINT),   parameter::CP_WALKSTAT_ACTIVE        = 1
  integer(KINDINT),   parameter::CP_WALKSTAT_REFLECT       = 2
  integer(KINDINT),   parameter::CP_WALKSTAT_TRANS         = 4
  integer(KINDINT),   parameter::CP_WALKSTAT_TRAPPED       = 8
  integer(KINDINT),   parameter::CP_WALKSTAT_WAITING       = 16

  !**** out put format
  character*7,        parameter::PKW_OUTSTAT_FORMATXYZ = "&CFGXYZ"
  character*9,        parameter::PKW_OUTSTAT_FORMAT20  = "&MCSTAT20"


  end module MC_Randwalk_Const

