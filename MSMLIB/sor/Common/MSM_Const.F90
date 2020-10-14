  module MSM_CONSTANTS
  !***  DESCRIPTION: this module is to define the constants to be used
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010

  implicit none

      !--- type define
      integer,          parameter::KINDDF      = 8
      integer,          parameter::KINDSF      = 4
      integer,          parameter::KINDINT     = 4
      integer,          parameter::KINDINT2    = 2
      integer,          parameter::KINDINTL    = 8
      integer,          parameter::KINDINTS    = KINDINT2

      !--- Memory size  
      integer,          parameter::C_BYTE      = 8    ! (8 bits)
      integer,          parameter::C_KBYTES    = 1024*C_BYTE
      integer,          parameter::C_MBYTES    = 1024*C_KBYTES
      integer,          parameter::C_GBYTES    = 1024*C_MBYTES
  
      !--- numbers                             
      real(KINDDF),     parameter::C_ZERO      = 0.D0
      real(KINDDF),     parameter::C_UTH       = 1.D0/3.D0
      real(KINDDF),     parameter::C_HALF      = 0.5D0
      real(KINDDF),     parameter::C_UN        = 1.D0
      real(KINDDF),     parameter::C_UPF       = 1.5D0
      real(KINDDF),     parameter::C_TWO       = 2.D0
      real(KINDDF),     parameter::C_THR       = 3.D0
      real(KINDDF),     parameter::C_FOUR      = 4.D0
      real(KINDDF),     parameter::C_FIVE      = 5.D0
      real(KINDDF),     parameter::C_SIX       = 6.D0
      real(KINDDF),     parameter::C_HUIT      = 8.D0
      real(KINDDF),     parameter::C_TEN       = 10.D0
      real(KINDDF),     parameter::C_SIXTIN    = 16.D0
      real(KINDDF),     parameter::C_4OVER3    = 4.D0/3.D0
     
      integer(KINDINT), parameter::C_IZERO     = 0
      integer(KINDINT), parameter::C_IUN       = 1
      integer(KINDINT), parameter::C_ITWO      = 2
      integer(KINDINT), parameter::C_ITHR      = 3
      integer(KINDINT), parameter::C_IDIX      = 10
      integer(KINDINT), parameter::C_ICENT     = 100

      real(KINDDF),     parameter::CP_RSQ2     = 1.414213562D0
      real(KINDDF),     parameter::CP_HALFRSQ2 = 0.707106781187D0
      real(KINDDF),     parameter::CP_RSQ3     = 1.732050808D0
      real(KINDDF),     parameter::CP_HALFRSQ3 = CP_RSQ3/2.D0
      real(KINDDF),     parameter::CP_RSQ4     = 2.D0

      integer(KINDINT), parameter::CP_HWORD    = -131072             ! 0xFFFE0000   bits 11111111111111110000000000000000
      integer(KINDINT), parameter::CP_LWORD    =  65535              ! 0x0000FFFF        00000000000000001111111111111111
      integer(KINDINT), parameter::CP_HBIT     = -2147483648         ! 0x80000000        10000000000000000000000000000000

      !--- math constants
      real(KINDDF),     parameter::CP_PI       = 3.141592654D0
      real(KINDDF),     parameter::CP_4PI3     = 4.D0*CP_PI/3.D0
      real(KINDDF),     parameter::CP_4PIOVER3 = CP_4PI3
      real(KINDDF),     parameter::CP_TWOPI    = 2.D0*CP_PI
      real(KINDDF),     parameter::CP_FOURPI   = 4.D0*CP_PI

     !--- anglur
     real(KINDDF),      parameter::CP_DEG2ARC  = CP_PI/180.D0
     real(KINDDF),      parameter::CP_ARC2DEG  = 180.D0/CP_PI

     !--- length and time
      real(KINDDF),     parameter::CP_A2CM     = 1.D-8              ! Unit conversion: Anstro to Centimeter
      real(KINDDF),     parameter::CP_CM2A     = 1.D8               ! Unit conversion: Centimeter to Anstro
      real(KINDDF),     parameter::CP_CM2NM    = 1.D7               ! Unit conversion: Centimeter to nanometer
      real(KINDDF),     parameter::CP_NM2CM    = 1.D-7              ! Unit conversion: nanometer to centimeter
      real(KINDDF),     parameter::CP_SQRCM2NM  = CP_CM2NM*CP_CM2NM ! Unit conversion: square centimeter to square nanometer
      real(KINDDF),     parameter::CP_SQRTNM2CM = CP_NM2CM*CP_NM2CM ! Unit conversion: nanometer to centimeter

      real(KINDDF),     parameter::CP_AU2G     = 1.66053D-24        ! Mass unit conversion: atomic unit to CGS
      real(KINDDF),     parameter::CP_G2AU     = 1.D0/CP_AU2G       ! Mass unit conversion: CGS to atomic unit
      real(KINDDF),     parameter::CP_FS2S     = 1.D-15             ! Time unit conversion: fecosecond to second
      real(KINDDF),     parameter::CP_PS2S     = 1.D-12             ! Time unit conversion: picosecond to second
      real(KINDDF),     parameter::CP_S2FS     = 1.D+15             ! Time unit conversion: second to fecosecond
      real(KINDDF),     parameter::CP_S2PS     = 1.D+12             ! Time unit conversion: second to pecosecond
                                                                    
      !--- phys constants
      real(KINDDF),     parameter::CP_KB       = 1.38054D-16        ! Boltzmann constant (erg/K)
      real(KINDDF),     parameter::CP_EVERG    = 1.60219D-12        ! unit conversion from eV->Erg
      real(KINDDF),     parameter::CP_EV2ERG   = CP_EVERG           ! unit conversion from eV->Erg
      real(KINDDF),     parameter::CP_ERGEV    = 1.D0/CP_EVERG      ! unit conversion from Erg->eV
      real(KINDDF),     parameter::CP_ERG2EV   = CP_ERGEV           ! unit conversion from Erg->eV
      real(KINDDF),     parameter::CP_KEVERG   = CP_EVERG*1.D3      ! unit conversion from keV->Erg
      real(KINDDF),     parameter::CP_KEV2ERG  = CP_KEVERG          ! unit conversion from keV->Erg
      real(KINDDF),     parameter::CP_ERGKEV   = 1.D0/CP_KEVERG     ! unit conversion from Erg->keV
      real(KINDDF),     parameter::CP_ERG2KEV  = CP_ERGKEV          ! unit conversion from Erg->keV
      real(KINDDF),     parameter::CP_EV2KEV   = 1.D3               ! unit conversion from ev->keV
      real(KINDDF),     parameter::CP_KEV2EV   = 1.D-3              ! unit conversion from keV->ev
      real(KINDDF),     parameter::CP_KEV2MEV  = 1.D3               ! unit conversion from kev->MeV
      real(KINDDF),     parameter::CP_MEV2KEV  = 1.D-3              ! unit conversion from MeV->kev
                                               
      real(KINDDF),     parameter::CP_EV2K     = CP_EVERG/CP_KB     ! unit conversion of kinetice energy ERG->K
                                               
      real(KINDDF),     parameter::CP_DYN2EVA  = CP_ERGEV/CP_CM2A   ! unit conversion of force DYN->eV/A
      real(KINDDF),     parameter::CP_BOHR     = 5.291772108D-9     ! Bohr radius (in unit cm)
      real(KINDDF),     parameter::CP_AVOG     = 6.0221367D23       ! Avigado constants (/mol))
         
      !--- unit conversion for pressure                                         
      real(KINDDF),     parameter::CP_CGS2BAR  = 1.0D-6
      real(KINDDF),     parameter::CP_CGS2KBAR = 1.0D-9
      real(KINDDF),     parameter::CP_BAR2CGS  = 1.0D6
      real(KINDDF),     parameter::CP_KBAR2CGS = 1.0D9
      real(KINDDF),     parameter::CP_BARCGS   = 1.0D6              ! cgs dyn/cm^2
      real(KINDDF),     parameter::CP_BAR2PA   = 1.0D5              ! cgs dyn/cm^2
      real(KINDDF),     parameter::CP_BAREVA   = 1.D-18/CP_EVERG    ! ev/A**3
                       
      real(KINDDF),     parameter::CP_KE2TEMP  = C_TWO*C_UTH/CP_KB  ! conversion from kinetic energy to temperature
                       
  end module MSM_CONSTANTS

