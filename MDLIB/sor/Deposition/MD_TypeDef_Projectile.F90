  module MD_TypeDef_Projectile
  !***  DESCRIPTION: this module is to define the data type for Depostion Control Parameter.
  !
   !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  2010-11 (Hou Qing)
  !
  !                  2018-01 (Hou Qing) modified to accept keyword driven input
  !
  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none

   !--- TAG for input filename
      character(len=13),parameter, private::mp_FTAGI="&AUXF_DEPOSIT"

   !--- The type describes the data controlling the deposition
   !    simulation
      integer,parameter,private::CP_DEP_STYPE_EXT          = 0   !external depsotion
      integer,parameter,private::CP_DEP_STYPE_INT          = 1   !internal depsotion
      integer,parameter,private::CP_DEP_STYPE_PKA          = 2   !internal PKA depsotion

      integer,parameter,private::CP_EK_STYLE_MONO   = 0          ! the incident energy style
      integer,parameter,private::CP_EK_STYLE_GAUSS  = 1          ! gauss energy spectrum
      integer,parameter,private::CP_EK_STYLE_EXT    = 2          ! external energy spectrum

      integer,parameter,private::CP_ORIENT_STYLE_MONODIR   = 0   !incident orientation is in given polar angle and azimuth angle
      integer,parameter,private::CP_ORIENT_STYLE_POLAR     = 1   !incident orientation is in given polar angle but random azimuth angle
      integer,parameter,private::CP_ORIENT_STYLE_RANDPOLAR = 3   !incident orientation is in randomly selected polar angle
      integer,parameter,private::CP_ORIENT_STYLE_CENTER    = 4   !incident oritenation is in dirtection from entry point to box center
      integer,parameter,private::CP_ORIENT_STYLE_LATTDIR   = 5   !the incident orientation is given in lattice direction

      integer,parameter,private::CP_ENTRY_STYLE_BOX        = 1   !the incident position is at given place
      integer,parameter,private::CP_ENTRY_STYLE_CENTER     = 2   !the incident position is in the center of box for internale irradiation
                                                                 ! or at the top surface center for external irraidtion

      integer,parameter,private::CP_ENTRY_STYLE_CORNER     = 3   !the incident position is at the left-up corner of the box
      !--- the parameters for an incident atoms
      type ProjectileParam
           integer         ::DEPSTYLE = CP_DEP_STYPE_EXT        !deposition style =0 , external; =1 internal
           integer         ::ATYPE    = 0                       !type of incident atom
           character(len=8)::SYMB     = ""
           !--- parameters for incident energy spectrum
           integer      ::EKSTYPE    = CP_EK_STYLE_MONO          !type of energy spectrum =0 for monoenergy
           real(KINDDF) ::EK         =0.D0                       !average incident kinetic energy for each kind of atoms
           real(KINDDF) ::EKSD       =0.D0                       !with of guass distribution of incident energy
           character*256::EKSFNAME   =""                         !files for energy spectrum

           !--- parameters for incident angular distribution
           integer      ::DIRTYPE    = CP_ORIENT_STYLE_MONODIR   !type of incident dirtection
           real(KINDDF) ::THETA      =0.D0                       !polar angle of incident orientations
           real(KINDDF) ::FAI        =0.D0                       !azimuth angle of incident orientation
           real(KINDDF) ::THETAMI    =0.D0                       !range of polar angle for RANDPOLAR incident
           real(KINDDF) ::THETAMX    =0.D0                       !range of polar angle for RANDPOLAR incident

           real(KINDDF) ::MILLER(3)  =0                          ! the Miller index lattice direction that the energtic atoms directed in.
                                                                 ! with ANDType = CP_ORIENT_STYLE_GAUSS
           character*256::ANDFNAME   =""                         !files for angular distrubtion

           !--- parameters for incident position
           integer      ::ENTRYTYPE   = CP_ENTRY_STYLE_BOX       !type of entry position
           real(KINDDF) ::ENTRYBOX(6) = (/-0.5D0, 0.5D0, &       !the box entry position, the incident atom will be randomly placed
                                          -0.5D0, 0.5D0, &
                                          -0.5D0, 0.5D0/)
           !--- parameters for incident position
           real(KINDDF)  ::RATIO       = 0.D0
      end type ProjectileParam

      !--- the parameters for deposition calculations
      type(ProjectileParam), dimension(:), allocatable::m_PROJ    !
     !--------------------------------------------------

  contains
  !**********************************************************************************
  !**********************************************************************************
   subroutine Initialize_Projectile(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate memories, and create the
   !               neighbor-list of the reference lattice
   !
   use MD_TYPEDEF_InputPaser, only:Get_InputStatements
   implicit none
  !----   DUMMY Variables
   type(SimMDBox),  intent(in)::SimBox
   type(SimMDCtrl), intent(in)::CtrlParam

   !----   Local variable
    integer::I, IFILE, STATU
    type(InputStatements)::INPUTS

         !$$--- to findout the I/O unit
          IFILE = 0
          do I=1, size(CtrlParam%f_tag)
             if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                IFILE = I
                exit
             end if
          enddo

         !$$--- check the input files
          if(IFILE .le. 0) then
             write(*,fmt="(A)") ' MDPSCU Warning: the control file for deposition is missed'
             write(*,fmt="(A)") '               add the keyword in SETUP file:' //mp_FTAGI
             call ONWARNING(gm_OnWarning)
             write(*,fmt="(A)") ' Process to be continue without depositing any atom'
             return
          end if

         !$$--- load input statements
          write(*,*) "!**** Load control parameters for deposition simulation from:"
          write(*,*) "!**** ",CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
          write(*,*)

          call LoadExInput_SimMDCtrl(mp_FTAGI, CtrlParam%f_others(IFILE), CtrlParam)
          call GetExInput_SimMDCtrl(CtrlParam, mp_FTAGI, INPUTS, STATU)
          call ExtractParameter_Projectile(SimBox, CtrlParam, INPUTS, m_PROJ)

          !$$---
          write(6,*) "!************ PROJECTILE CONTROL PARAMETERS **********"
          call PrintParameter_Projectile(6, m_PROJ)

          if(gm_hFILELOG .gt. 0) call PrintParameter_Projectile(gm_hFILELOG, m_PROJ)
          call Release_InputStatements(INPUTS)
         return
   end subroutine Initialize_Projectile
  !**********************************************************************************

  !**********************************************************************************
   subroutine ExtractParameter_Projectile(SimBox, CtrlParam, Input, Proj)
   !***  PURPOSE:  to extract control parameter for projectiles
   !
   implicit none
  !----   DUMMY Variables
    type(SimMDBox),                                 intent(in)::SimBox
    type(SimMDCtrl),                                intent(in)::CtrlParam
    type(InputStatements),                          intent(in)::Input
    type(ProjectileParam),dimension(:),allocatable            ::Proj

   !----   Local variable
    integer::LINE, NP, N, I, J
    real(KINDDF)::FSWAP1, FSWAP2
    character*256::STR, SUBSTR(2)
    character*32::KWD, KWD1,STRNUMB(9)
    equivalence(STRNUMB(1), SUBSTR(2))
    type(InputStatements)::tIn


          !--- to check how many kind of incident atoms
          NP = 0
          do while (.true.)
             write(STRNUMB(1),*) NP+1
             STRNUMB(1) = adjustl(STRNUMB(1))
             KWD  = "&PROJ"//STRNUMB(1)(1:len_trim(STRNUMB(1)))
             KWD1 = "&ENDPROJ"//STRNUMB(1)(1:len_trim(STRNUMB(1)))
             call Copy_InputStatements(Input, tIn, StartKWD=KWD, EndKWD=KWD1)
             if(.not.associated(tIn%stat)) then
                exit
             end if
             !--- to check if the input structure is complete
             if(.not. HasKeyword_InputStatements(KWD1, tIn) ) then
                write(*,fmt="(A)")    ' MDPSCU Error: projectile data is not ended with keyword '//KWD1(1:len_trim(KWD1))
                call Print_StatementList(6, tIn%stat)
                write(*,fmt="(A)")    '               Check you control file for depostion calculations'
                write(*,fmt="(A)")    'Process to be stopped'
                stop
             end if
             NP = NP + 1
          end do
          if(NP .le. 0) then
             write(*,fmt="(A)")    ' MDPSCU Error: projectile data is missed'
             write(*,fmt="(A)")    '               Check you control file for depostion calculations'
             write(*,fmt="(A)")    'Process to be stopped'
             stop
          end if

          if(allocated(Proj)) then
             deallocate(Proj)
          end if
          allocate(Proj(NP))

          !--- extracte the projectile properties
          do I=1, NP
             write(STRNUMB(1),*) I
             STRNUMB(1) = adjustl(STRNUMB(1))
             KWD  = "&PROJ"//STRNUMB(1)(1:len_trim(STRNUMB(1)))
             KWD1 = "&ENDPROJ"//STRNUMB(1)(1:len_trim(STRNUMB(1)))
             call Copy_InputStatements(Input, tIn, StartKWD=KWD, EndKWD=KWD1)

             !--- extract irradiation style
             call Get_InputStatements("&INCIDENT_TYPE", tIN, STR, LINE)
             if(LINE .eq. 0) call Get_InputStatements("&INCID_TYPE", tIN, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Substr(STR,2,N,STRNUMB)
                call UpCase(STRNUMB(1))

                select case ( STRNUMB(1)(1:len_trim(STRNUMB(1)) ) )
                       case ("EXT", "EXTERNAL")
                             Proj(I)%DEPSTYLE = CP_DEP_STYPE_EXT

                       case ("INT", "INTERNAL")
                              Proj(I)%DEPSTYLE = CP_DEP_STYPE_INT

                       case ("PKA")
                              Proj(I)%DEPSTYLE = CP_DEP_STYPE_PKA

                       case default
                             write(*,fmt="(A)")       ' MDPSCU Error: wrong incident projectile type.'
                             write(*,fmt="(A)")       '               The incident type should be one of "EXT", "INT", "PKA"'
                             write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                             write(*,*) 'Process to be stopped'
                             stop
                end select
                Proj(I)%ATYPE = 0
                if(N .ge. 2) then
                   call UpCase(STRNUMB(2))
                   do J=1, SimBox%NGROUP
                      STRNUMB(1) = SimBox%SYMB(J)
                      call UpCase(STRNUMB(1))
                      if(STRNUMB(2)(1:len_trim(STRNUMB(2))) .eq. STRNUMB(1)(1:len_trim(STRNUMB(1))) ) then
                         Proj(I)%ATYPE = J
                         exit
                      end if
                   end do
                   if(Proj(I)%ATYPE .le. 0) then
                      write(*,fmt="(A)") ' MDPSCU Error: wrong atomic symbol "'//STRNUMB(2)(1:len_trim(STRNUMB(2)))//'"'
                      write(*,fmt="(A)") ' Usage: &INCID_TYPE "'//STRNUMB(1)(1:len_trim(STRNUMB(1)))//'" asymb'
                      write(*,fmt="(A)") '        where asymb is the atomic symbol of incident atom'
                      write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                      write(*,*) 'Process to be stopped'
                      stop
                   end if
                else
                   call Extract_Numb(STR,1,N,STRNUMB)
                   if(N .ge. 1) then
                      Proj(I)%ATYPE = ISTR(STRNUMB(1))
                   else
                      write(*,fmt="(A)") ' MDPSCU Error: the type of projectile is missed'
                      write(*,fmt="(A)") ' Usage: &INCID_TYPE "'//STRNUMB(1)(1:len_trim(STRNUMB(1)))//'" atype'
                      write(*,fmt="(A)") '        where atype is the type ID of incident atom'
                      write(*,*) 'Process to be stopped'
                      stop
                   end if
                end if

                if( Proj(I)%ATYPE .le. 0 .or.  Proj(I)%ATYPE .gt. SimBox%NGROUP) then
                   write(*,fmt="(A, I4, A)") ' MDPSCU Error: the type of projectile',Proj(I)%ATYPE, ' is not defined in simulation box'
                   write(*,fmt="(A, I4, A)") '               the number of groups in SIMBOX is ', SimBox%NGROUP
                   write(*,*) 'Process to be stopped'
                   stop
                end if
                Proj(I)%SYMB = SimBox%SYMB(Proj(I)%ATYPE)
             else
                write(*,fmt="(A, I4)") ' MDPSCU Error: cannot find incident type for of projectile #', I
                write(*,fmt="(A)")     ' Usage: &INCID_TYPE radtyp atype, in your deposition file'
                write(*,fmt="(A)")     '        where radtype is the type type of irradition,'
                write(*,fmt="(A)")     '        which should be one of "EXT", "INT", "PKA";'
                write(*,fmt="(A)")     '        atype is the type ID of incident atom'
                write(*,fmt="(A)")     ' Process to be stopped'
                stop
             end if

             !*** second to extract the incident energies
             Proj(I)%EKSType = CP_EK_STYLE_MONO
             Proj(I)%EK      = 0.D0
             Proj(I)%EKSD    = 0.D0
             Proj(I)%EKSFNAME= ""
             call Get_InputStatements("&EN_SPECT", tIN, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Substr(STR,1,N,STRNUMB)
                if(N .gt. 0) then
                   call UpCase(STRNUMB(1))
                   select case(STRNUMB(1)(1:len_trim(STRNUMB(1)) ) )
                          case ("MONO")
                                Proj(I)%EKSTYPE = CP_EK_STYLE_MONO
                                call Extract_Numb(STR,1,n,STRNUMB)
                                if(N.ge.1) then
                                   Proj(I)%EK = DRSTR(STRNUMB(1))
                                end if

                          case ("GAUSS")
                                Proj(I)%EKSTYPE = CP_EK_STYLE_GAUSS
                                call Extract_Numb(STR,2,N,STRNUMB)
                                if(N.ge.1) then
                                   Proj(I)%EK = DRSTR(STRNUMB(1))
                                   if(N .ge. 2) then
                                      Proj(I)%EKSD = DRSTR(STRNUMB(2))
                                   end if
                                end if
                                if(Proj(I)%EKSD .le. 0.D0) then
                                   Proj(I)%EKSTYPE = CP_EK_STYLE_MONO
                                end if

                          case ("EXT", "EXTERN", "EXTERNAL")
                                Proj(I)%EKSTYPE = CP_EK_STYLE_EXT
                                call Extract_Substr(STR,2,N,SUBSTR)
                                if(N .ge. 2) then
                                   Proj(I)%EKSFNAME = SUBSTR(2)
                                else
                                   write(*,fmt="(A)") ' MDPSCU Error: the file name for energy spectrum of incident atoms is missed'
                                   write(*,fmt="(A)") '               Usage: &EN_SPECT "EXT", fname'
                                   write(*,fmt="(A)") '                      where fname is the name of file storing energy spectrum of the incident atom'
                                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                                   write(*,fmt="(A)") 'Process to be stopped'
                                   stop
                                end if
                          case default
                             write(*,fmt="(A)")       ' MDPSCU Error: wrong energy spectrum type of incident projectile'
                             write(*,fmt="(A)")       '               The spectrum type should be one of "MONO", "GAUSS", "EXT"'
                             write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                             write(*,*) 'Process to be stopped'
                             stop
                   end select
                else
                   write(*,fmt="(A)")       ' MDPSCU Error: energy spectrum type of incident projectile is missed'
                   write(*,fmt="(A)")       '               Usage: &EN_SPECT st'
                   write(*,fmt="(A)")       '                      where st is one of "MONO", "GAUSS", "EXT"'
                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                   write(*,*) 'Process to be stopped'
                   stop
                end if
             end if

             !*** to extract entry point in normalized unit
             if( Proj(I)%DEPSTYLE .eq. CP_DEP_STYPE_EXT) then
                 Proj(I)%ENTRYTYPE = CP_ENTRY_STYLE_CENTER
             else if( Proj(I)%DEPSTYLE .eq. CP_DEP_STYPE_INT .or. &
                     Proj(I)%DEPSTYLE .eq. CP_DEP_STYPE_PKA ) then
                 Proj(I)%ENTRYTYPE = CP_ENTRY_STYLE_CORNER
             end if
             Proj(I)%ENTRYBOX = (/-0.5D0, 0.5D0, -0.5D0, 0.5D0, -0.5D0, 0.5D0/)

             call Get_InputStatements("&INCID_POS", tIN, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Substr(STR,1,N,STRNUMB)
                if(N .gt. 0) then
                   call UpCase(STRNUMB(1))
                   select case(STRNUMB(1)(1:len_trim(STRNUMB(1)) ) )
                          case ("BOX")
                                Proj(I)%ENTRYTYPE = CP_ENTRY_STYLE_BOX
                                call Extract_Numb(STR,6,N,STRNUMB)
                                if(N .eq. 1) then
                                   Proj(I)%ENTRYBOX(1) = DRSTR(STRNUMB(1))
                                else if(N .eq. 2) then
                                   Proj(I)%ENTRYBOX(1) = DRSTR(STRNUMB(1))
                                   Proj(I)%ENTRYBOX(2) = DRSTR(STRNUMB(2))
                                else if(N .eq. 3) then
                                   Proj(I)%ENTRYBOX(1) = DRSTR(STRNUMB(1))
                                   Proj(I)%ENTRYBOX(2) = DRSTR(STRNUMB(2))
                                   Proj(I)%ENTRYBOX(3) = DRSTR(STRNUMB(3))
                                else if(N .eq. 4) then
                                   Proj(I)%ENTRYBOX(1) = DRSTR(STRNUMB(1))
                                   Proj(I)%ENTRYBOX(2) = DRSTR(STRNUMB(2))
                                   Proj(I)%ENTRYBOX(3) = DRSTR(STRNUMB(3))
                                   Proj(I)%ENTRYBOX(4) = DRSTR(STRNUMB(4))
                                else if(N .eq. 5) then
                                   Proj(I)%ENTRYBOX(1) = DRSTR(STRNUMB(1))
                                   Proj(I)%ENTRYBOX(2) = DRSTR(STRNUMB(2))
                                   Proj(I)%ENTRYBOX(3) = DRSTR(STRNUMB(3))
                                   Proj(I)%ENTRYBOX(4) = DRSTR(STRNUMB(4))
                                   Proj(I)%ENTRYBOX(5) = DRSTR(STRNUMB(5))
                                else if(N .ge. 6) then
                                   Proj(I)%ENTRYBOX(1) = DRSTR(STRNUMB(1))
                                   Proj(I)%ENTRYBOX(2) = DRSTR(STRNUMB(2))
                                   Proj(I)%ENTRYBOX(3) = DRSTR(STRNUMB(3))
                                   Proj(I)%ENTRYBOX(4) = DRSTR(STRNUMB(4))
                                   Proj(I)%ENTRYBOX(5) = DRSTR(STRNUMB(5))
                                   Proj(I)%ENTRYBOX(6) = DRSTR(STRNUMB(6))
                                end if
                                Proj(I)%ENTRYBOX(1) = max(Proj(I)%ENTRYBOX(1), -0.5D0)
                                Proj(I)%ENTRYBOX(2) = min(Proj(I)%ENTRYBOX(2),  0.5D0)
                                FSWAP1 = min(Proj(I)%ENTRYBOX(1), Proj(I)%ENTRYBOX(2))
                                FSWAP2 = max(Proj(I)%ENTRYBOX(1), Proj(I)%ENTRYBOX(2))
                                Proj(I)%ENTRYBOX(1) = FSWAP1
                                Proj(I)%ENTRYBOX(2) = FSWAP2

                                Proj(I)%ENTRYBOX(3) = max(Proj(I)%ENTRYBOX(3), -0.5D0)
                                Proj(I)%ENTRYBOX(4) = min(Proj(I)%ENTRYBOX(4),  0.5D0)
                                FSWAP1 = min(Proj(I)%ENTRYBOX(3), Proj(I)%ENTRYBOX(4))
                                FSWAP2 = max(Proj(I)%ENTRYBOX(3), Proj(I)%ENTRYBOX(4))
                                Proj(I)%ENTRYBOX(3) = FSWAP1
                                Proj(I)%ENTRYBOX(4) = FSWAP2

                                Proj(I)%ENTRYBOX(5) = max(Proj(I)%ENTRYBOX(5), -0.5D0)
                                Proj(I)%ENTRYBOX(6) = min(Proj(I)%ENTRYBOX(6),  0.5D0)
                                FSWAP1 = min(Proj(I)%ENTRYBOX(5), Proj(I)%ENTRYBOX(5))
                                FSWAP2 = max(Proj(I)%ENTRYBOX(6), Proj(I)%ENTRYBOX(6))
                                Proj(I)%ENTRYBOX(5) = FSWAP1
                                Proj(I)%ENTRYBOX(6) = FSWAP2

                          case ("CENTER")
                                Proj(I)%ENTRYTYPE = CP_ENTRY_STYLE_CENTER
                          case ("CORNER")
                                Proj(I)%ENTRYTYPE = CP_ENTRY_STYLE_CORNER
                   end select
                else
                   write(*,fmt="(A)")       ' MDPSCU Error: parameters for incident position is missed'
                   write(*,fmt="(A)")       '               Usage: &INCID_POS st, x0, x1, y0, y1, z0, z1'
                   write(*,fmt="(A)")       '                      where st is one of "BOX", "CENTER", "CORNER"'
                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                   write(*,*) 'Process to be stopped'
                   stop
                end if
             end if

             !*** to extract angular distribution
             Proj(I)%DIRTYPE = CP_ORIENT_STYLE_MONODIR
             Proj(I)%THETA   = 0.D0
             Proj(I)%FAI     = 0.D0
             Proj(I)%THETAMI = 0.D0
             Proj(I)%THETAMX = 0.D0
             call Get_InputStatements("&INCID_DIR", tIN, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Substr(STR,1,N,STRNUMB)
                if(N .gt. 0) then
                   call UpCase(STRNUMB(1))
                   select case(STRNUMB(1)(1:len_trim(STRNUMB(1)) ) )
                          case ("MONO")
                                Proj(I)%DIRTYPE = CP_ORIENT_STYLE_MONODIR
                                call Extract_Numb(STR,2,N,STRNUMB)
                                if(N .ge. 2) then
                                   Proj(I)%THETA  = DRSTR(STRNUMB(1))
                                   Proj(I)%FAI    = DRSTR(STRNUMB(2))
                                else
                                   write(*,fmt="(A)") ' MDPSCU Error: incident orientation is missed'
                                   write(*,fmt="(A)") '               Usage: &INCID_DIR "MONO", theta, fai'
                                   write(*,fmt="(A)") '                      where theta is polar angle and fai the azimuth angle'
                                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                                   write(*,fmt="(A)") 'Process to be stopped'
                                   stop
                                end if

                          case ("POLAR")
                                Proj(I)%DIRTYPE = CP_ORIENT_STYLE_POLAR
                                call Extract_Numb(STR,1,N,STRNUMB)
                                if(N .ge. 1) then
                                   Proj(I)%THETA  = DRSTR(STRNUMB(1))
                                else
                                   write(*,fmt="(A)") ' MDPSCU Error: incident orientation is missed'
                                   write(*,fmt="(A)") '               Usage: &INCID_DIR "POLAR", theta'
                                   write(*,fmt="(A)") '                      where theta is polar angle'
                                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                                   write(*,fmt="(A)") 'Process to be stopped'
                                   stop
                                end if

                          case ("RANDPOLAR")
                                Proj(I)%DIRTYPE = CP_ORIENT_STYLE_RANDPOLAR
                                call Extract_Numb(STR,2,N,STRNUMB)
                                if(N .eq. 1) then
                                   Proj(I)%THETAMI = 0.D0
                                   Proj(I)%THETAMX = DRSTR(STRNUMB(1))
                                else if(N .ge. 2) then
                                   Proj(I)%THETAMI = min(DRSTR(STRNUMB(1)), DRSTR(STRNUMB(2)) )
                                   Proj(I)%THETAMX = max(DRSTR(STRNUMB(1)), DRSTR(STRNUMB(2)) )
                                else if(N .lt. 1) then
                                   write(*,fmt="(A)") ' MDPSCU Error: range of polar angle is missed'
                                   write(*,fmt="(A)") '               Usage: &INCID_DIR "RANDPOLAR", r'
                                   write(*,fmt="(A)") '                      where r is the range for polar angle'
                                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                                   write(*,fmt="(A)") 'Process to be stopped'
                                   stop
                                end if

                          case ("CENTER")
                                Proj(I)%DIRTYPE = CP_ORIENT_STYLE_CENTER

                          case default
                               write(*,fmt="(A)")       ' MDPSCU Error: wrong parameter for incident direction '
                               write(*,fmt="(A)")       '               Usage: &INCID_DIR st, var'
                               write(*,fmt="(A)")       '                      where st is one of "MONO", "POLAR", "RANDPOLAR", "CENTER"'
                               write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                               write(*,*) 'Process to be stopped'
                               stop

                   end select
                else
                   write(*,fmt="(A)")       ' MDPSCU Error: parameters for incident angle is missed'
                   write(*,fmt="(A)")       '               Usage: &INCID_DIR st, var'
                   write(*,fmt="(A)")       '                      where st is one of "MONO", "POLAR", "RANDPOLAR", "CENTER"'
                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                   write(*,*) 'Process to be stopped'
                   stop
                end if
             end if

             !*** to extract percentage of the projectile
             Proj(I)%RATIO =   0.D0
             call Get_InputStatements("&INCID_RATIO", tIN, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,1,N,STRNUMB)
                if(N .ge. 1) then
                   Proj(I)%RATIO  = DRSTR(STRNUMB(1))
                else
                   write(*,fmt="(A)") ' MDPSCU Error: ratio for the projectile is missed'
                   write(*,fmt="(A)") '               Usage: &INCID_RATIO r'
                   write(*,fmt="(A)") '                      where r is ratio of the projectile in the beam'
                   write(*,fmt="(A, BZI6)") ' Check input file at line:', LINE
                   write(*,fmt="(A)") 'Process to be stopped'
                   stop
                end if
             end if
             !*****
          end do

          !*** check input consistency
          FSWAP1 = 0.D0
          do I=1, size(Proj)
             FSWAP1 = FSWAP1 + Proj(I)%RATIO
          end do

          if(FSWAP1 .le. 0.D0) then
             write(*,fmt="(A)") ' MDPSCU Warning: percentage for all projectile are zero'
             write(*,fmt="(A)") '                 no projectile will be issued'
             call ONWARNING(gm_OnWarning)
          else
             do I=1, size(Proj)
                Proj(I)%RATIO =  Proj(I)%RATIO/FSWAP1
             end do
          end if

          call Release_InputStatements(tIN)
         return
   end subroutine ExtractParameter_Projectile
  !****************************************************************************

  !****************************************************************************
  subroutine PrintParameter_Projectile(hFile, Proj)
  !***  PURPOSE:   to print projectile parameters into I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
      type(ProjectileParam), dimension(:)::Proj

     !--- local variables
      character*256::STR,STRTMP(1)=""
      character*32::STRNUM, STRNUM1
      integer::I, J, N
      real(KINDDF)::RATIO
     !----
             write(STRNUM,*) size(Proj)
             STRNUM = adjustl(STRNUM)
             write(hFile, fmt="(A)") ' !    Number of projectile types..............: '//STRNUM(1:len_trim(STRNUM))
             do I=1, size(Proj)
                !--- irradiation style
                write(STRNUM,*) I
                STRNUM = adjustl(STRNUM)
                write(hFile, fmt="(A)") ' !      Projecile ...........................#: '//STRNUM(1:len_trim(STRNUM))&
                                        //', symbol "'//Proj(I)%SYMB(1:len_trim(Proj(I)%SYMB))//'"'
                select case(Proj(I)%DEPSTYLE)
                       case (CP_DEP_STYPE_EXT)
                             write(hFile, fmt="(A)") ' !        irradiation style...................: external'
                       case (CP_DEP_STYPE_INT, CP_DEP_STYPE_PKA)
                             write(hFile, fmt="(A)") ' !        irradiation style...................: internal'
                       !case (CP_DEP_STYPE_PKA)
                       !      write(hFile, fmt="(A)") ' !        irradiation style...................: PKA'
                end select

                !--- energy spectrum
                select case(Proj(I)%EKSTYPE)
                       case(CP_EK_STYLE_MONO)
                            write(STRNUM, fmt="(1PE12.3)") Proj(I)%EK
                            STRNUM = adjustl(STRNUM)
                            write(hFile, fmt="(A)") ' !          energy spectrum ..................: mono-energy, '//&
                                                    STRNUM(1:len_trim(STRNUM))//'(keV)'
                       case(CP_EK_STYLE_GAUSS)
                            write(STRNUM, fmt="(1PE12.3)") Proj(I)%EK
                            STRNUM = adjustl(STRNUM)
                            write(STRNUM1, fmt="(1PE12.3)") Proj(I)%EKSD
                            STRNUM1 = adjustl(STRNUM1)
                            write(hFile, fmt="(A)") ' !          energy spectrum ..................: gaussian, '//&
                                                                 STRNUM(1:len_trim(STRNUM))//'(keV)'
                            write(hFile, fmt="(A)") ' !                                   with width: '//STRNUM1(1:len_trim(STRNUM1))
                end select

                !--- entry point
                select case(Proj(I)%ENTRYTYPE)
                       case(CP_ENTRY_STYLE_CENTER)
                            write(hFile, fmt="(A)") ' !          incident at ......................: box center'
                       case(CP_ENTRY_STYLE_CORNER)
                            write(hFile, fmt="(A)") ' !          incident at ......................: box left-up corner'
                       case(CP_ENTRY_STYLE_BOX)
                            if(Proj(I)%DEPSTYLE .eq. CP_DEP_STYPE_EXT) then
                            write(hFile, fmt="(A, 3(F5.2,','), F5.2, A)") ' !          incident in region (norm-nuit)....: (', &
                                                      Proj(I)%ENTRYBOX(1:4) , ')'
                            else
                            write(hFile, fmt="(A, 5(F5.2,','), F5.2, A)") ' !          incident in box (norm-nuit).......: (', &
                                                      Proj(I)%ENTRYBOX(1:6) , ')'
                            end if
                end select

                !--- entry direction
                select case(Proj(I)%DIRTYPE)
                       case(CP_ORIENT_STYLE_MONODIR)
                            write(hFile, fmt="(A, 2(F6.1,','), A)") &
                                                   ' !          in fixed direction.................: (', &
                                                    Proj(I)%THETA, Proj(I)%FAI, ') in deg'
                       case(CP_ORIENT_STYLE_POLAR)
                            write(hFile, fmt="(A, 1(F6.1,','), A)") &
                                                   ' !          in fixed polar angle...............: ', &
                                                    Proj(I)%THETA,  'in deg'
                       case(CP_ORIENT_STYLE_RANDPOLAR)
                            write(hFile, fmt="(A, 1(F6.1,','), F6.1, A)") &
                                                   ' !          in range of polar angle...........: (', &
                                                    Proj(I)%THETAMI, Proj(I)%THETAMX,')in deg'
                       case(CP_ORIENT_STYLE_CENTER)
                            write(hFile, fmt="(A)") ' !         in direction to center of box......: '
                end select

             end do

             !--- pecentage of projectiles
             do I=1, size(Proj)
                write(STRNUM,*) I
                STRNUM = adjustl(STRNUM)
                write(hFile, fmt="(A, F7.3)") ' !    Percentage of projectile..............#'//STRNUM(1:len_trim(STRNUM))&

                                         //':', Proj(I)%RATIO*100.D0
             end do
             write(hFile,*)

         return
  end subroutine PrintParameter_Projectile
  !****************************************************************************

  !****************************************************************************
  subroutine Create_Projectile(SimBox, CtrlParam, Proj, Pos, Vel)
  !***  PORPOSE: to deposit an atom in a box
  !     INPUT:  SimBox,  the box array to be created
  !             Proj,    the projectiles
  !
   use RAND32_MODULE
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox)         ::SimBox
       type(SimMDCtrl)        ::CtrlParam
       type(ProjectileParam)  ::Proj
       real(KINDDF)           ::POS(3),VEL(3)
  !--- local varibales
       integer::I, ITYP, NA
       real(KINDDF)::CENTER(3), BOX(6), ZMAX, EN, VEL0


            !---
            ZMAX      = -1.D64
            CENTER(3) = 0.D0
            NA        = 0
            do I =1, SimBox%NPRT
               if(IAND(SimBox%STATU(I),CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX) then
                  NA = NA + 1
                  CENTER(1:3) = CENTER(1:3) + SimBox%XP(I,1:3)

                  if(ZMAX .lt. SimBox%XP(I,3)) then
                     ZMAX = SimBox%XP(I,3)
                  end if
               end if
            end do
            CENTER = CENTER/dble(NA)
            BOX(1) = 0.5D0*(SimBox%BOXUP(1) + SimBox%BOXLOW(1)) + &
                     Proj%ENTRYBOX(1)*(SimBox%BOXUP(1) - SimBox%BOXLOW(1))
            BOX(2) = 0.5D0*(SimBox%BOXUP(1) + SimBox%BOXLOW(1)) + &
                     Proj%ENTRYBOX(2)*(SimBox%BOXUP(1) - SimBox%BOXLOW(1))
            BOX(3) = 0.5D0*(SimBox%BOXUP(2) + SimBox%BOXLOW(2)) + &
                     Proj%ENTRYBOX(3)*(SimBox%BOXUP(2) - SimBox%BOXLOW(2))
            BOX(4) = 0.5D0*(SimBox%BOXUP(2) + SimBox%BOXLOW(2)) + &
                     Proj%ENTRYBOX(4)*(SimBox%BOXUP(2) - SimBox%BOXLOW(2))
            BOX(5) = 0.5D0*(SimBox%BOXUP(3) + SimBox%BOXLOW(3)) + &
                     Proj%ENTRYBOX(5)*(SimBox%BOXUP(3) - SimBox%BOXLOW(3))
            BOX(6) = 0.5D0*(SimBox%BOXUP(3) + SimBox%BOXLOW(3)) + &
                     Proj%ENTRYBOX(6)*(SimBox%BOXUP(3) - SimBox%BOXLOW(3))

            !--- determine the position of the projectile
             select case(Proj%ENTRYTYPE)
                    case(CP_ENTRY_STYLE_CENTER)
                         POS(1) = CENTER(1)
                         POS(2) = CENTER(2)
                         if(Proj%DEPSTYLE .eq. CP_DEP_STYPE_EXT) then
                            POS(3) = ZMAX + maxval(CtrlParam%RU)
                         else
                            POS(3) = CENTER(3)
                         end if
                    case(CP_ENTRY_STYLE_CORNER)
                         POS(1) = SimBox%BOXLOW(1)
                         POS(2) = SimBox%BOXLOW(2)
                         if(Proj%DEPSTYLE .eq. CP_DEP_STYPE_EXT) then
                            POS(3) = ZMAX + maxval(CtrlParam%RU)
                         else
                            POS(3) = ZMAX - maxval(CtrlParam%RU)
                         end if
                    case(CP_ENTRY_STYLE_BOX)
                         POS(1) = BOX(1) + DRAND32()*(BOX(2)-BOX(1))
                         POS(2) = BOX(3) + DRAND32()*(BOX(4)-BOX(3))
                         if(Proj%DEPSTYLE .eq. CP_DEP_STYPE_EXT) then
                            POS(3) = ZMAX + maxval(CtrlParam%RU)
                         else
                            POS(3) = BOX(5) + DRAND32()*(BOX(6)-BOX(5))
                         end if
             end select

            !--- determine the incident direction
             select case(Proj%DIRTYPE)
                    case(CP_ORIENT_STYLE_MONODIR)
                         VEL(1) = dsin(Proj%THETA*CP_DEG2ARC)*dcos(Proj%FAI*CP_DEG2ARC)
                         VEL(2) = dsin(Proj%THETA*CP_DEG2ARC)*dsin(Proj%FAI*CP_DEG2ARC)
                         VEL(3) = dcos(Proj%THETA*CP_DEG2ARC)
                         VEL    = -VEL
                       case(CP_ORIENT_STYLE_POLAR)
                         Proj%FAI   = DRAND32()*CP_TWOPI*CP_ARC2DEG
                         VEL(1) = dsin(Proj%THETA*CP_DEG2ARC)*dcos(Proj%FAI*CP_DEG2ARC)
                         VEL(2) = dsin(Proj%THETA*CP_DEG2ARC)*dsin(Proj%FAI*CP_DEG2ARC)
                         VEL(3) = dcos(Proj%THETA*CP_DEG2ARC)
                         VEL    = -VEL
                       case(CP_ORIENT_STYLE_RANDPOLAR)
                         Proj%THETA = Proj%THETAMI + DRAND32()*(Proj%THETAMX - Proj%THETAMI)
                         Proj%FAI   = DRAND32()*CP_TWOPI*CP_ARC2DEG
                         VEL(1) = dsin(Proj%THETA*CP_DEG2ARC)*dcos(Proj%FAI*CP_DEG2ARC)
                         VEL(2) = dsin(Proj%THETA*CP_DEG2ARC)*dsin(Proj%FAI*CP_DEG2ARC)
                         VEL(3) = dcos(Proj%THETA*CP_DEG2ARC)
                         VEL    = -VEL
                       case(CP_ORIENT_STYLE_CENTER)
                         VEL(1) = CENTER(1) - POS(1)
                         VEL(2) = CENTER(2) - POS(2)
                         VEL(3) = CENTER(3) - POS(3)
                         VEL    = VEL/dsqrt(sum(VEL*VEL))
                end select

            !--- determine the incident direction
                EN   = Proj%EK*CP_KEVERG
                VEL0 = DSQRT(C_TWO*EN/SimBox%CM(Proj%ATYPE))
                VEL  = VEL*VEL0
                select case(Proj%EKSTYPE)
                       case (CP_EK_STYLE_MONO)
                       case (CP_EK_STYLE_GAUSS)
                       case (CP_EK_STYLE_EXT)
                end select


           return
  end subroutine Create_Projectile
  !****************************************************************************
  !****************************************************************************
  subroutine Deposit_Projectile(SimBox, CtrlParam, RESTART)
  !***  PORPOSE: to deposit an atom in a box
  !     INPUT:  SimBox,  the box array to be created
  !             CtrlParam, the control parameter
  !             RESTART, indictor to indicating if restarting a new session
   use RAND32_MODULE
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer                     ::RESTART
  !--- local varibales
       integer::I, J, IP
       real(KINDDF)::R, SR, POS(1,3), VEL(1,3)

           do I=1, size(SimBox)
              if(RESTART .eq. 0) then
                 R  = DRAND32()
                 SR = 0.d0
                 IP = 0
                 do J=1, size(m_PROJ)
                    SR = SR + m_PROJ(J)%RATIO
                    if(SR .ge. R) then
                       IP = J
                       exit
                    end if
                 end do
                 if(IP .gt. 0) then
                    call Create_Projectile(SimBox(I), CtrlParam, m_PROJ(IP), POS, VEL)
                    select case (m_PROJ(IP)%DEPSTYLE )
                       case (CP_DEP_STYPE_EXT, CP_DEP_STYPE_INT)
                             call AddAtoms_SimMDBox(SimBox(I),1,m_PROJ(IP)%ATYPE, TYPEORDER=1, RXP=POS, RXP1=VEL)

                       case (CP_DEP_STYPE_PKA)
                              call ReplaceAtom_SimMDBox(SimBox(I), POS(1,:), m_PROJ(IP)%ATYPE, VEL(1,:))
                    end select
                 end if
              else
                 J= size(m_PROJ)
                 if(sum(m_PROJ(1:J)%RATIO) .gt. 0.d0) then
                    call AddAtoms_SimMDBox(B=SimBox(I), N=RESTART, ITYPE=SimBox(I)%NGROUP, TYPEORDER=1)
                 end if
              end if
           end do

           return
  end subroutine Deposit_Projectile
  !****************************************************************************************

  end module MD_TypeDef_Projectile
