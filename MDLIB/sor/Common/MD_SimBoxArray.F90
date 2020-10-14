  module MD_SimboxArray
  !***  DESCRIPTION: this module is to define the data type for Simulation Box.
  !                  Thie version is adopted from an MDGvar.f90( version 2000)
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none

    !---------------------------------------------------------
    private:: Copy_SimBoxArray0,         &
              CopyFrom_SimBoxArray0
    public::  Copy_SimBoxArray
    interface Copy_SimBoxArray
           module procedure Copy_SimBoxArray0
    end interface Copy_SimBoxArray 
    
    !--- interface to the external routine -------------------
    interface assignment (=)
       module procedure CopyFrom_SimBoxArray0
    end interface

    !---------------------------------------------------------
    private:: Initialize_Config_SimBoxArray0,          &
              Initialize_Config_SimBoxArray1
    public::  Initialize_Config_SimBoxArray
    interface Initialize_Config_SimBoxArray
           module procedure Initialize_Config_SimBoxArray0
           module procedure Initialize_Config_SimBoxArray1
    end interface Initialize_Config_SimBoxArray

    !---------------------------------------------------------
    private:: Putin_Instance_Config_SimBoxArray0,      &
              Putin_Instance_Config_SimBoxArray1,      &
              Putin_Instance_Config_SimBoxArray2
    public::  Putin_Instance_Config_SimBoxArray 
    interface Putin_Instance_Config_SimBoxArray
            module procedure Putin_Instance_Config_SimBoxArray0 
            module procedure Putin_Instance_Config_SimBoxArray1
            module procedure Putin_Instance_Config_SimBoxArray2
    end interface Putin_Instance_Config_SimBoxArray
    !---------------------------------------------------------

    !---------------------------------------------------------
    private:: Putout_Instance_Config_SimBoxArray0,      &
              Putout_Instance_Config_SimBoxArray1,      &
              Putout_Instance_Config_SimBoxArray2
    public::  Putout_Instance_Config_SimBoxArray 
    interface Putout_Instance_Config_SimBoxArray
            module procedure Putout_Instance_Config_SimBoxArray0 
            module procedure Putout_Instance_Config_SimBoxArray1
            module procedure Putout_Instance_Config_SimBoxArray2
    end interface Putout_Instance_Config_SimBoxArray
    !---------------------------------------------------------


  contains
  !*********************************************************************
      subroutine Release_SimBoxArray(B)
      !
      implicit none
      type(SimMDBox), dimension(:),allocatable::B

      !--- local variables
          integer::I, NB
          
          if(.not.allocated(B)) return
          NB = size(B)
          do I=1, NB
             call Release_SimMDBox(B(I))
          end do
          return
       end subroutine Release_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Create_SimBoxArray(SimBox, BS, SimBox0)
  !***  PURPOSE:   to create a simulation box array from a template box
  !
  !     INPUT:    SimBox0,   the template
  !               BS,        optional for the size of the array
  !     OUTPUT:   SimBox, the resultant array
  implicit none
     !--- dummy vaiables
       type(SimMDBox), dimension(:), allocatable::SimBox
       integer,                         optional::BS
       type(SimMDBox), intent(in),      optional::SimBox0
     !--- local variables
     integer::I

       if(present(BS)) then
          if(allocated(SimBox)) then
             if(BS .NE. SIZE(SimBox)) then
                call Release_SimBoxArray(SimBox)
                deallocate(SimBox)
             end if
          end if
          if( .not. allocated(SimBox)) then
            allocate(SimBox(BS))
          end if
       end if

       if(present(SimBox0)) then
          do I=1, size(SimBox)
             call Copy_SimMDBox(SimBox0, SimBox(I))
          end do
       end if
       return
  end subroutine Create_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Copy_SimBoxArray0(Sor, Tgt)
  !***  PURPOSE:   to copy a simulation box array from sor to tgt
  !
  !     INPUT:    Sor,    the source
  !               
  !     OUTPUT:   Tgt,    the target
  implicit none
     !--- dummy vaiables
       type(SimMDBox), dimension(:), intent(in) ::Sor
       type(SimMDBox), dimension(:), allocatable::Tgt
     !--- local variables
     integer::I

       if(allocated(Tgt)) then
          if(size(Tgt) .ne. size(sor)) then
            call Release_SimBoxArray(Tgt)
            deallocate(Tgt)
          end if
       end if   

       if(.not.allocated(Tgt)) then
          allocate(Tgt(size(Sor))) 
       end if   
       do I=1, size(Sor)
          call Copy_SimMDBox(Sor(I),Tgt(I))
       end do
       return
  end subroutine Copy_SimBoxArray0
  !****************************************************************************

  !****************************************************************************
  subroutine CopyFrom_SimBoxArray0(Tgt, Sor)
  !***  PURPOSE:   to copy a simulation box array from sor to tgt
  !
  !     INPUT:    Sor,    the source
  !               
  !     OUTPUT:   Tgt,    the target
  implicit none
     !--- dummy vaiables
       type(SimMDBox), dimension(:), intent(in) ::Sor
       type(SimMDBox), dimension(:), allocatable::Tgt
     !--- local variables

       call Copy_SimBoxArray0(Sor, Tgt) 
       return
  end subroutine CopyFrom_SimBoxArray0
  !****************************************************************************


  !*********************************************************************
   subroutine AddOneAtom_SimBoxArray(SimBox, ITYPE, DISMI0, CENTER0, INR0, OUTR0 )
      !***  PURPOSE:   to add new atoms to a simulation boxs
      !     INPUT:     SimBox, the box array
      !                itype,  the type of the atoms to be added
      !                DISMI0, the permitted min distance between atoms
      !                CENTER0, the center that the atom will be placed around
      !                INR0,    the radiu around the center, out  which the atom will be placed
      !                OUTR0,   the radiu around the center, in which the atom will be placed
      !     OUTPUT     SimBox, the box with number of atoms updated
      use RAND32_MODULE
      implicit none
      type(SimMDBox),dimension(:)      ::SimBox
      integer,               intent(in)::ITYPE
      real(KINDDF),           optional ::DISMI0, CENTER0(3), INR0, OUTR0
      !--- LOCAL variabels
      integer::I, J, L, NPRT
      integer, parameter::MXLOOP=200
      real(KINDDF)::XP(1,3)
      real(KINDDF)::DISMI, CENTER(3), INR, OUTR, DIS, RMIN

                 if(present(DISMI0) ) then
                    DISMI = DISMI0*DISMI0
                 else
                    DISMI = 1.D64 !$$---without limit
                 end if

                 if(present(CENTER0) ) then
                    CENTER = CENTER0
                 else
                    CENTER = C_ZERO
                 end if

                 if(present(INR0) ) then
                    INR = INR0*INR0
                 else
                    INR = 1.D64  !$$--- without limit
                 end if

                 if(present(OUTR0) ) then
                    OUTR = OUTR0*OUTR0
                 else
                    OUTR = C_ZERO
                 end if
  !---
                 do I=1, SIZE(SimBox)
                    NPRT = SimBox(I)%NPRT
                    RMIN = 1.D64
                    do L =1, MXLOOP
                       !$$--- to check if two atoms too close
                       XP(1,1) = DRAND32()*SimBox(I)%ZL(1)+SimBox(I)%BOXLOW(1)
                       XP(1,2) = DRAND32()*SimBox(I)%ZL(2)+SimBox(I)%BOXLOW(2)
                       XP(1,3) = DRAND32()*SimBox(I)%ZL(3)+SimBox(I)%BOXLOW(3)

                       DIS = SUM((XP(1,1:3)-CENTER(1:3))*(XP(1,1:3)-CENTER(1:3)))
                       if(DIS .LT. OUTR) cycle
                       if(DIS .GT. INR)  cycle

                       do J=1, NPRT
                          DIS = SUM((XP(1,1:3)-SimBOX(I)%XP(J,1:3))*(XP(1,1:3)-SimBOX(I)%XP(J,1:3)))
                          if(DIS .lt. RMIN) RMIN = DIS
                       end do
                       IF(RMIN .LT. DISMI .AND. DISMI.LE.1.D30) cycle
                       exit
                    end do
                    if(L .GT. MXLOOP) then
                       print *, "MDPSCU Warning: Added atom not meet limitations", DSQRT(RMIN), "vs", DSQRT(DISMI)
                      call ONWARNING(gm_OnWarning)
                    end if

                    call AddAtoms_SimMDBox(SimBox(I), 1, ITYPE, TYPEORDER=1, RXP=XP)
                 end do


  end subroutine AddOneAtom_SimBoxArray
  !****************************************************************************

  !*********************************************************************
   subroutine AppendOneAtom_SimBoxArray(SimBox, ITYPE, DISMI0, CTYPE0, INR0, OUTR0 )
      !***  PURPOSE:   to append a new atom to exist cluster
      !     INPUT:     SimBox, the box array
      !                itype,  the type of the atoms to be added
      !                DISMI0, the permitted min distance between atoms
      !                CENTER0, the center that the atom will be placed around
      !                INR0,    the radiu around the center, in which the atom will be placed
      !                OUTR0,   the radiu around the center, our of which the atom will be placed
      !     OUTPUT     SimBox, the box with number of atoms updated
      use RAND32_MODULE
      implicit none
      type(SimMDBox), dimension(:)::SimBox

      integer,      intent(in)::ITYPE
      integer,      optional  ::CTYPE0
      real(KINDDF), optional  ::DISMI0, INR0, OUTR0
      !--- LOCAL variabels
      integer::I, J, L, NPRT, IBEST, CTYPE, NA0, IA0
      integer, parameter::MXLOOP=200
      real(KINDDF)::XP(1,3)
      real(KINDDF)::DISMI,DISMI2, CENTER0(3), CENTER(3),CSIZE, INR, OUTR, DIS, DIR(3), RMIN, RMIMX, XPT(MXLOOP,3)

                 if(present(DISMI0) ) then
                    DISMI = DISMI0
                 else
                    DISMI = 1.D64 !---without limit
                 end if
                 DISMI2 = DISMI*DISMI

                 if(present(CTYPE0)) then
                    CTYPE = CTYPE0
                 else
                    CTYPE = ITYPE
                 end if

                 !$$--- determine the region in which the atom to be placed
                 if(present(OUTR0) ) then
                    OUTR = OUTR0
                 else
                    OUTR = DISMI
                 end if

                 if(present(INR0) ) then
                    INR = INR0
                 else
                    INR = DISMI
                 end if

  !---
                 do I=1, SIZE(SimBox)
                    NPRT = SimBox(I)%NPRT
                    RMIN = 1.D64
                    RMIMX=0.D0
                    !--- first to get the center of mass of the cluster
                    !--- get the initial center of mass
                    IA0 = 0
                    NA0 = 0
                    CENTER0 = 0
                    do J = 1, NPRT
                       if(SimBox(I)%ITYP(J).eq. CTYPE) then
                          CENTER0(1:3) = SimBox(I)%XP(J,1:3)
                          IA0 = J
                          NA0 = 1
                          exit
                       end if
                    end do

                    CENTER = CENTER0
                    do J = 1, NPRT
                       if(SimBox(I)%ITYP(J).eq. CTYPE .and. J.ne.IA0) then
                          XP(1,1:3) = SimBox(I)%XP(J,1:3) - CENTER0(1:3)

                          if(dabs(XP(1,1)) .GT. C_HALF*SimBox(I)%ZL(1) ) &
                             XP(1,1) = XP(1,1) - DSIGN(SimBox(I)%ZL(1), XP(1,1))

                          if(dabs(XP(1,2)) .GT. C_HALF*SimBox(I)%ZL(2) ) &
                             XP(1,2) = XP(1,2) - DSIGN(SimBox(I)%ZL(2), XP(1,2))

                          if(dabs(XP(1,3)) .GT. C_HALF*SimBox(I)%ZL(3) ) &
                             XP(1,3) = XP(1,3) - DSIGN(SimBox(I)%ZL(3), XP(1,3))

                          CENTER(1:3) = CENTER(1:3) + (XP(1,1:3) + CENTER0(1:3))
                          NA0 = NA0 + 1
                          CENTER0(1:3) = CENTER(1:3)/dble(NA0)

                        end if
                    end do

                    do J=1, 3
                       if(CENTER0(J) .LT. SimBox(I)%BOXLOW(J) ) then
                          CENTER0(J) = CENTER0(J) + SimBox(I)%ZL(J)
                       else if(CENTER0(J) .GT. SimBox(I)%BOXUP(J)) then
                          CENTER0(J) = CENTER0(J)-SimBox(I)%ZL(J)
                       end if
                    end do
                    CENTER = CENTER0


                    !$$--- then to get the radiu of the cluster
                    CSIZE = 0.D0
                    do J=1, 1, SimBox(I)%NPRT
                       if(SimBox(I)%ITYP(J).eq. CTYPE) then
                          XP(1,1:3) = SimBox(I)%XP(J,1:3) - CENTER(1:3)

                          if(dabs(XP(1,1)) .GT. C_HALF*SimBox(I)%ZL(1) ) &
                             XP(1,1) = XP(1,1) - DSIGN(SimBox(I)%ZL(1), XP(1,1))

                          if(dabs(XP(1,2)) .GT. C_HALF*SimBox(I)%ZL(2) ) &
                             XP(1,2) = XP(1,2) - DSIGN(SimBox(I)%ZL(2), XP(1,2))

                          if(dabs(XP(1,3)) .GT. C_HALF*SimBox(I)%ZL(3) ) &
                             XP(1,3) = XP(1,3) - DSIGN(SimBox(I)%ZL(3), XP(1,3))

                          DIS = SUM(XP(1,1:3)*XP(1,1:3))
                          if(CSIZE .le. DIS) CSIZE = DIS
                       end if
                    end do
                    INR = CSIZE+INR
                    OUTR = CSIZE+OUTR

                    !$$--- the distance from the center
                    DIS =  INR+DRAND32()*(OUTR-INR)
                    do L =1, MXLOOP

                       !$$--- the direction placed
                       DIR(1) = DRAND32() - C_HALF
                       DIR(2) = DRAND32() - C_HALF
                       DIR(3) = DRAND32() - C_HALF
                       DIR = DIR/DSQRT(SUM(DIR*DIR))
                       do J=1, 3
                          XPT(L,J) = CENTER(J)+DIR(J)*DIS
                          if(XPT(L,J) .LT. SimBox(I)%BOXLOW(J) ) then
                             XPT(L,J) = XPT(L,J)+SimBox(I)%ZL(J)
                          else if(XPT(L,J) .GT. SimBox(I)%BOXUP(J)) then
                             XPT(L,J) = XPT(L,J)-SimBox(I)%ZL(J)
                          end if
                       end do

                       !$$--- to check if the atom is too close to other atoms
                       do J=1, NPRT
                          DIR(1) = XPT(L,1)-SimBOX(I)%XP(J,1)
                          if(DABS(DIR(1)) .GE. C_HALF*SimBox(I)%ZL(1)) then
                             DIR(1) = DIR(1) - DSIGN(SimBox(I)%ZL(1),DIR(1))
                          end if

                          DIR(2) = XPT(L,2)-SimBOX(I)%XP(J,2)
                          if(DABS(DIR(2)) .GE. C_HALF*SimBox(I)%ZL(2)) then
                             DIR(2) = DIR(2) - DSIGN(SimBox(I)%ZL(2),DIR(2))
                          end if

                          DIR(3) = XPT(L,3)-SimBOX(I)%XP(J,3)
                          if(DABS(DIR(3)) .GE. C_HALF*SimBox(I)%ZL(3)) then
                             DIR(3) = DIR(3) - DSIGN(SimBox(I)%ZL(3),DIR(3))
                          end if

                          DIS = SUM(DIR*DIR)
                          if(DIS .lt. RMIN) RMIN = DIS
                       end do

                       if(RMIN .LT. DISMI2) then
                          if(RMIN.GT.RMIMX) then
                             IBEST = L
                             RMIMX = RMIN
                          end if
                          cycle
                       else
                          IBEST = L
                          exit
                       end if
                    end do
                    if(L .GT. MXLOOP) then
                       print *, "WARNING: Added atom not meet limitations", DSQRT(RMIMX), "vs", DSQRT(DISMI2)
                       call ONWARNING(gm_OnWarning)
                    end if

                    XP(1,1:3) = XPT(IBEST,1:3)
                    call AddAtoms_SimMDBox(SimBox(I), 1, ITYPE, TYPEORDER=1, RXP=XP)
                 end do

                 return
  end subroutine AppendOneAtom_SimBoxArray
  !****************************************************************************

  !*********************************************************************
   subroutine ReplaceOneAtom_SimBoxArray(SimBox, ITYPE, IA, OLDTY )
      !***  PURPOSE:   to change the type of an atom in simulation boxs
      !     INPUT:     SimBox, the box array
      !                itype,  the new type of the atom
      !                IA,     the ID of the atom to be changed
      !                OLDTY,  the original type of the atom
      !     OUTPUT     SimBox, the box with number of atoms updated
      use RAND32_MODULE
      implicit none
      type(SimMDBox), dimension(:)::SimBox

      integer,  intent(in)         ::ITYPE
      integer,  intent(in),optional::IA, OLDTY
      !--- LOCAL variabels
      integer::I, J, L, NPRT

                 !$$--- if IA given we changed type of IA
                 if(present(IA) ) then
                    do I=1, SIZE(SimBox)
                       call ReplaceAtom_SimMDBox(SimBox(I), IA, ITYPE)
                    end do
                    return
                 end if

                 !$$--- if old type is given, we replace an randomly choicen atom of this type
                 if(present(OLDTY) ) then
                    do I=1, SIZE(SimBox)
                       J = DRAND32()*SimBox(I)%NA(OLDTY)+1
                       call ReplaceAtom_SimMDBox(SimBox(I), J, ITYPE)
                    end do
                    return
                 end if

                 !$$---  randomly replace the type an atom in the whole box
                 do I=1, SIZE(SimBox)
                    J = DRAND32()*SimBox(I)%NPRT+1
                    call ReplaceAtom_SimMDBox(SimBox(I), J, ITYPE)
                 end do
                 return

  end subroutine ReplaceOneAtom_SimBoxArray
  !****************************************************************************

  !*********************************************************************
   subroutine DeleteOneAtom_SimBoxArray(SimBox, ITYPE)
      !***  PURPOSE:   to delete the an atom of give type to make a vacancy
      !     INPUT:     SimBox, the box array
      !                itype,  the  type of the atom
      !     OUTPUT     SimBox, the box with number of atoms updated
      use RAND32_MODULE
      implicit none
      type(SimMDBox), dimension(:)::SimBox
      integer,          intent(in)::ITYPE
      !--- LOCAL variabels
      integer::I, J, IA, ID


                 !$$---  randomly delete an atom of given type
                 do I=1, SIZE(SimBox)
                    if(SimBox(I)%NA(ITYPE) .eq. 0) cycle

                    ID = DRAND32()*SimBox(I)%NA(ITYPE)+1
                    IA = 0
                    do J=1, SimBox(I)%NPRT
                       if(SimBox(I)%ITYP(J) .eq. ITYPE) then
                          IA = IA + 1
                          if(IA .eq. ID) then
                             call DelAtom_SimMDBox(SimBox(I), J)
                             exit
                          end if
                       end if
                    end do
                 end do
                 return

  end subroutine DeleteOneAtom_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_Config_SimBoxArray0(SimBox, fname, fmt, multbox)
  !***  PURPOSE:   to initialize configuration by create a new or
  !                restore a configuration from previous calcualtion.
  !                The configuration to be read from file SimBox%IniConfig.
  !
  !     INPUT:    fname,   the file storing the configuration
  !
  !     OUTPUT    SimBox,  the simulation box
  implicit none
     !--- dummy vaiables
       type(SimMDBox), dimension(:)::SimBox

       character*(*), intent(in)::fname
       integer,       intent(in)::fmt
       integer,       optional  ::multbox

     !--- local variables
       integer::I, hFile, LINE
       character*256::TSTR
       character*32::KEYWORD0, KEYWORD1
       logical::EX

       !$$-- to check if the file existing
       inquire(FILE=fname, EXIST=EX)
       if(.not.EX) then
          write(*,*) "MDPSCU Error: cannot find initial configuration file ", fname(1:len_TRIM(fname))
          stop
       else
          write(*,fmt="(A)") ' MDPSCU Message: Load from '//fname(1:len_trim(fname))   
       end if
       call AvailableIOUnit(hFile)
       open(unit=hFile, file=fname, status='old')

        LINE = 0
        call GetInputStrLine(hFile,TSTR, LINE, "!", *100)
        TSTR = adjustl(TSTR)
        call GetKeyWord("&", TSTR, KEYWORD0)
        call UpCase(KEYWORD0)
        select case(KEYWORD0(1:LEN_TRIM(KEYWORD0)))
            !--- has the header
            case(PKW_OUTCFG_FORMAT14, PKW_OUTCFG_FORMAT15,PKW_OUTCFG_FORMAT16, PKW_OUTCFG_FORMAT18)
                 if(present(multbox)) then
                    if(multbox .gt. 0) then
                       !--- first, we check the consistent of input configure
                       LINE = 0
                       rewind(hFILE)
                       !--- read in header of file and check the number of atoms is consistent
                       !    change of number of atoms is not allowed
                       call Read_ConfigFileHeader_SimMDBox(Simbox(1), hFile, LINE, CFG=KEYWORD0)

                       LINE = 0
                       rewind(hFILE)
                       do I=1, SIZE(SimBox)
                          do while(.true.)  !skip the header
                             call GetInputStrLine(hFile,TSTR, LINE, "!", *100)
                             TSTR = adjustl(TSTR)
                             call GetKeyWord("&", TSTR, KEYWORD1)
                             call UpCase(KEYWORD1)
                            if(KEYWORD1(1:LEN_TRIM(KEYWORD1)) .EQ. "&TYPE") exit
                          end do
                          call Read_Initial_Config_SimMDBox0(Simbox(I), hFile, fmt, LINE)
                          !--- NOTE: for PKW_OUTCFG_FORMAT18, the unit of velocity and force
                          !          needed to be converted to CGS
                          if(KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq.  PKW_OUTCFG_FORMAT18) then
                             if(iand(fmt,CP_INPUT_VEL) .gt. C_IZERO) &
                                Simbox(I)%XP1 = Simbox(I)%XP1*SimBox(I)%RR*CP_S2PS
                          end if
                       end do
                    else
                       !--- first, we check the consistent of input configure
                       LINE = 0
                       rewind(hFILE)
                       !--- read in header of file and check the number of atoms is consistent
                       !    change of number of atoms is not allowed
                       call Read_ConfigFileHeader_SimMDBox(Simbox(1), hFile, LINE, CFG=KEYWORD0)
                       LINE = 0
                       rewind(hFILE)
                       do while(.true.)  !skip the header
                          call GetInputStrLine(hFile,TSTR, LINE, "!", *100)
                          TSTR = adjustl(TSTR)
                          call GetKeyWord("&", TSTR, KEYWORD1)
                          call UpCase(KEYWORD1)
                          if(KEYWORD1(1:LEN_TRIM(KEYWORD1)) .EQ. "&TYPE") exit
                       end do

                       call Read_Initial_Config_SimMDBox0(Simbox(1), hFile, fmt, LINE)
                       !--- NOTE: for PKW_OUTCFG_FORMAT18, the unit of velocity and force
                       !          needed to be converted to CGS
                       if(KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq.  PKW_OUTCFG_FORMAT18) then
                          if(iand(fmt,CP_INPUT_VEL) .gt. C_IZERO) &
                             Simbox(1)%XP1 = Simbox(1)%XP1*SimBox(1)%RR*CP_S2PS
                       end if
                       do I=2, SIZE(SimBox)
                          call Copy_SimMDBox(SimBox(1), SimBox(I))
                       end do
                    end if
                 else  !--- all the box use the same initial configure
                    !--- first, we check the consistent of input configure
                    LINE = 0
                    rewind(hFILE)
                    !--- read in header of file and check the number of atoms is consistent
                    !    change of number of atoms is not allowed
                    call Read_ConfigFileHeader_SimMDBox(Simbox(1), hFile, LINE, CFG=KEYWORD0)
                    LINE = 0
                    rewind(hFILE)
                    do while(.true.)  !skip the header
                       call GetInputStrLine(hFile,TSTR, LINE, "!", *100)
                       TSTR = adjustl(TSTR)
                       call GetKeyWord("&", TSTR, KEYWORD1)
                       call UpCase(KEYWORD1)
                       if(KEYWORD1(1:LEN_TRIM(KEYWORD1)) .EQ. "&TYPE") exit
                    end do

                    call Read_Initial_Config_SimMDBox0(Simbox(1), hFile, fmt, LINE)
                    !--- NOTE: for PKW_OUTCFG_FORMAT18, the unit of velocity and force
                    !          needed to be converted to CGS
                    if(KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq.  PKW_OUTCFG_FORMAT18) then
                       if(iand(fmt,CP_INPUT_VEL) .gt. C_IZERO) &
                          Simbox(1)%XP1 = Simbox(1)%XP1*SimBox(1)%RR*CP_S2PS
                    end if

                    do I=2, SIZE(SimBox)
                       call Copy_SimMDBox(SimBox(1), SimBox(I))
                    end do
                end if

            case(PKW_OUTCFG_FORMATXYZ)
                 rewind(hFile)
                 !--- read in initial configuration from XYZ configure,
                 !    change of number of atoms is not allowed
                 call Read_ConfigByCol18_SimMDBox(Simbox(1), hFile)
                 do I=2, SIZE(SimBox)
                     call Copy_SimMDBox(SimBox(1), SimBox(I))
                 end do

            !--- in old foramt without header
            case default
                 if(len_trim(KEYWORD0) .gt. 0) then
                    write(*,fmt="(A,A)")        ' MDPSCU Error: unknown file format: ', KEYWORD0(1:len_trim(KEYWORD0))
                    write(*,fmt="(A,A, A, I5)") '               check file ', fname(1:len_trim(fname)), ' at line:', LINE
                    write(*,fmt="(A)")          '               process to be stopped'
                    stop
                 end if

                !$$--- load confiuration without header
                !$$    skip the discription line, if there is
                 rewind(hFile)
                 do I=1, LINE -1
                    read(hFile,*) TSTR
                 end do

                 if(present(multbox)) then
                    do I=1, SIZE(SimBox)
                       call Read_Initial_Config_SimMDBox0(Simbox(I), hFile, fmt, LINE, mod=0)
                    end do
                 else  !--- all the box use the same initial configure
                    call Read_Initial_Config_SimMDBox0(Simbox(1), hFile, fmt, LINE, mod=0)
                    do I=2, SIZE(SimBox)
                       call Copy_SimMDBox(SimBox(1), SimBox(I))
                    end do
                end if
        end select
        close(hFile)
        return

  100   write(*,fmt="(A,A)")        ' MDPSCU Error: error on reading from initial configuration file '//fname(1:len_TRIM(fname))
        write(*,fmt="(A,A, A, I5)") '               check file ', fname(1:len_trim(fname)), ' at line:', LINE
        write(*,fmt="(A)")          '               process to be stopped'
        stop
      return
  end subroutine Initialize_Config_SimBoxArray0
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_Config_SimBoxArray1(SimBox, Fname, JobID, CfgID, Fmt)
  !***  PURPOSE:   to initialize configuration by create a new or
  !                restore a configuration from previous calcualtion.
  !                The configuration to be read from a serial files with
  !                name: Fname_BoxID.CfgID, where the BoxIDs are calculated
  !                from JobID
  !
  !     INPUT:    fname,   the file storing the configuration
  !
  !     OUTPUT    SimBox,  the simulation box
  implicit none
     !--- dummy vaiables
       type(SimMDBox), dimension(:)::SimBox
       character*(*)               ::fname
       integer, intent(in)         ::JobID, CfgID, Fmt

     !--- local variables
       integer::I
       character*256::CFile0

       !$$-- to check if the file existing
       do I=1, size(SimBox)
          call STRCATI(CFile0, Fname, "_",   (JobID-1)*size(SimBox)+I, 4)
          call STRCATI(CFile0, CFile0, ".",  CfgID,                    4)
          call Initialize_Config_SimBoxArray0(SimBox(I:I), CFile0, Fmt, multbox=0)
       end do          
       return
  end subroutine Initialize_Config_SimBoxArray1
  !****************************************************************************

  !****************************************************************************
  subroutine Putout_Instance_Config_SimBoxArray0(Fname, SimBox, Stamp, NB, IB)
  !***  PORPOSE: to putout the instant configuaration of a simulation box
  !
  !     INPUT:   Fname, the head of the file name storing the configuration
  !              SimBox, the simulation box
  !              Stamp,  the recording stamp
  !              NB, optional, the number of boxs to be putout
  !              IB, optional, the ID of simulation box to be put out
  !
  implicit none
      !--- DUMMY variables
      character*(*),                intent(in):: Fname
      type(SimMDBox),dimension(:),  intent(in):: SimBox
      type(MDRecordStamp)                     :: Stamp

      integer,             intent(in),optional:: NB
      integer,dimension(:),intent(in),optional:: IB
      !--- local variables
      integer, parameter::P_TL = 18, P_TAGL = 12
      integer::K, hFile, I,J, NP1, IS, NDAT, ID, NCOL
      character*1024::title, fmt, PSTR0, PSTR

      type(DataPad), dimension(:),   allocatable::pDat
      integer,       dimension(:,:), allocatable::DIM
      character*64                              ::TSN, TS, TSFMT
      type(DataPad), pointer                    ::tsDat
      type(DataPad)                             ::TagDat


       !$$--- to open the file
       call AvailableIOUnit(hFile)
       open(hFile, file = Fname, status='unknown')
       if(present(NB)) then
          NP1 = NB
       else
          NP1 = size(SimBox)
       end if

      !*** write out head and box informations
       write(hFile, fmt="(A)") PKW_OUTCFG_FORMAT18
       call Putout_RecordStamp(hFile, Stamp)

       !$$*** to prepare the format of output
        IS    = 0
        NCOL  = 1
        TSFMT = "(1x I4, 1x, I3, 1x, A)"
        title = ""
        fmt   = "("

        title = title(1:IS)//   "&TYPE         "
        IS    = IS        + len("&TYPE         ")
        fmt   = fmt(1:len_trim(fmt))//"I8,3X"
        TS    = "&TYPECOL"
        write(TSN,fmt=TSFMT)  NCOL, 1, '"I"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        title = title(1:IS)//   "POS(LU)  (x)             (y)              (z)         "
        IS    = IS        + len("POS(LU)  (x)             (y)              (z)         ")
        fmt   = fmt(1:len_trim(fmt))//",3(1PE17.8,1X)"
        TS    = "&XYZCOL"
        write(TSN,fmt=TSFMT)  NCOL, 3, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 3

        title = title(1:IS)//   "VEL(LU/ps)(vx)           (vy)             (vz)        "
        IS    = IS        + len("VEL(LU/ps)(vx)           (vy)             (vz)        ")
        fmt   = fmt(1:len_trim(fmt))//",3(1PE17.8,1X)"
        TS    = "&VELCOL"
        write(TSN,fmt=TSFMT) NCOL, 3, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 3

        title = title(1:IS)//   "STATU   "
        IS    = IS        + len("STATU   ")
        fmt   = fmt(1:len_trim(fmt))//",I6,2X"
        TS    = "&STATCOL"
        write(TSN,fmt=TSFMT) NCOL, 1, '"I"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        title = title(1:IS)//   "FOR(ev/LU)(fx)           (fy)             (fz)        "
        IS    = IS        + len("FOR(ev/LU)(fx)           (fy)             (fz)        ")
        fmt   = fmt(1:len_trim(fmt))//",3(1PE17.8,1X)"
        TS    = "&FPCOL"
        write(TSN,fmt=TSFMT) NCOL, 3, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 3

        title = title(1:IS)//   "POT(ev).          "
        IS    = IS        + len("POT(ev).          ")
        fmt   = fmt(1:len_trim(fmt))//",1(1PE17.8,1X)"
        TS    = "&EPOTCOL"
        write(TSN,fmt=TSFMT) NCOL, 1, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        title = title(1:IS)//   "K.E.(ev)          "
        IS    = IS        + len("K.E.(ev)          ")
        fmt   = fmt(1:len_trim(fmt))//",1(1PE17.8,1X)"
        TS    = "&EKINCOL"
        write(TSN,fmt=TSFMT) NCOL, 1, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 1

        !title = title(1:IS)//   "ENERGY(ev)        "
        !IS    = IS        + len("ENERGY(ev)        ")
        !fmt   = fmt(1:len_trim(fmt))//",1(1PE17.8,1X)"

        title = title(1:IS)//   "DISPLACE(dx )            (dy)               (dz)       "
        IS    = IS        + len("DISPLACE(dx )            (dy)               (dz)       ")
        fmt   = fmt(1:len_trim(fmt))//",3(1PE17.8,1X)"
        TS    = "&DISCOL"
        write(TSN,fmt=TSFMT) NCOL,3, '"D"'
        TS    = TS(1:P_TAGL)//TSN
        call New_DataPad(TS, tsDat)
        NCOL  = NCOL + 3

        NDAT = 0
        if(associated(SimBox(1)%ptrDatPad) ) then
           NDAT = NumberofData_DataPad(SimBox(1)%ptrDatPad)
           allocate(pDat(NDAT), DIM(NDAT,2))
           do ID=1, NDAT
              call GetData_DataPad(ID, SimBox(1)%ptrDatPad, pDat(ID))
              call GetSize_DataPad(ID, SimBox(1)%ptrDatPad, DIM(ID,:))
              title = title(1:IS)//pDat(ID)%Tag(1:len_trim(pDat(ID)%Tag))

              call GetFirstWord(pDat(ID)%Tag, TS)
              call UpCase(TS)
              TS    = "&"//TS(1:len_trim(TS))//"COL"
              if(IsInteger_DataPad(pDat(ID)) ) then
                 write(TSN,fmt=TSFMT) NCOL,Dim(ID,2), '"I"'
              else
                 write(TSN,fmt=TSFMT) NCOL,Dim(ID,2), '"D"'
              end if
              TS    = TS(1:P_TAGL)//TSN
              call New_DataPad(TS, tsDat)
              NCOL  = NCOL + DIM(ID,2)
              IS = IS + P_TL*DIM(ID,2)
           end do
       end if
       fmt   = fmt(1:len_trim(fmt))//",A)"

       !$$--- to write out colume information
        write(hFile,FMT="(A, 3(1PE12.5,1x))") "!--- Table colume information:"
        do K =1, NumberofData_DataPad(tsDat)
           call GetData_DataPad(K, tsDat, TagDat)
           write(hFile,FMT="(A))")TagDat%Tag(1:len_trim(TagDat%Tag))
        end do

       if(.not.present(IB)) then
          do I=1, NP1
             !$$--- write out box information
             write(hFile,FMT="(A))")
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "!--- Box information:"
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&LATT     lattice length (in A):       ",  SimBox(I)%RR*CP_CM2A
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&BOXLOW   low boundary of box (in LU): ",  SimBox(I)%BOXLOW(1:3)/SimBox(I)%RR
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&BOXSIZE  boxsize (in LU):             ",  SimBox(I)%ZL(1:3)/SimBox(I)%RR
             write(hFile,FMT="(A, 3(I8,1x))")      "&NATOM    total number of atoms:       ",  SimBox(I)%NPRT
             write(hFile,FMT="(A, 3(I8,1x))")      "&NGROUP   number of group of atoms:    ",  SimBox(I)%NGROUP
             write(hFile,FMT="(A, 100(I8,1x))")    "    &NA   number of atoms in groups:   ",  (SimBox(I)%NA(J), J=1, SimBox(I)%NGROUP)
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&TEMPCAL  instant temperature(K):      ", SimBox(I)%Temperature

             !$$--- write out the configures
             write(hFile,FMT="(20A))")title(1:len_trim(title))
             do K=1, SimBox(I)%NPRT
                PSTR = ""
                IS  = 0
                do ID=1, NDAT
                   !call GetData_DataPad(ID, SimBox(I)%ptrDatPad, pDat(ID))
                   !call GetSize_DataPad(ID, SimBox(I)%ptrDatPad, DIM)
                   if(IsInteger_DataPad(pDat(ID)) ) then
                      call Write_DataPad(pDat(ID), K, PSTR0, Fmt="5x,I8,5x")
                    else if(IsDouble_DataPad(pDat(ID)) ) then
                      call Write_DataPad(pDat(ID), K, PSTR0, Fmt="1PE17.8,1X")
                    end if
                    PSTR = PSTR(1:IS)//PSTR0(1:DIM(ID,2)*P_TL)
                    IS   = IS + DIM(ID,2)*P_TL
                end do

                 write(hFile,fmt=fmt)SimBox(I)%ITYP(K), SimBox(I)%XP(K,1:3)/SimBox(I)%RR,           &
                                     SimBox(I)%XP1(K,1:3)/SimBox(I)%RR/CP_S2PS, SimBox(I)%STATU(K), &
                                     SimBox(I)%FP(K,1:3)*CP_ERGEV*SimBox(I)%RR,                     &
                                    -SimBox(I)%EPOT(K)*CP_ERGEV, SimBox(I)%EKIN(K)*CP_ERGEV,        &
                                     SimBox(I)%DIS(K,1:3)/SimBox(I)%RR,                             &
                                     PSTR(1:len_trim(PSTR))

             end do
          end do

       else
          do I=1, NP1
             !$$--- write out boxinformation
             write(hFile,FMT="(A))")
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&LATT     lattice length (in A):       ",  SimBox(IB(I))%RR*CP_CM2A
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&BOXLOW   low boundary of box (in LU): ",  SimBox(IB(I))%BOXLOW(1:3)/SimBox(IB(I))%RR
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&BOXSIZE  boxsize (in LU):             ",  SimBox(IB(I))%ZL(1:3)/SimBox(IB(I))%RR
             write(hFile,FMT="(A, 3(I8,1x))")      "&NATOM    total number of atoms:       ", SimBox(IB(I))%NPRT
             write(hFile,FMT="(A, 3(I8,1x))")      "&NGROUP   number of group of atoms:    ", SimBox(IB(I))%NGROUP
             write(hFile,FMT="(A, 100(I8,1x))")    "  &NA     number of atoms in groups:   ", (SimBox(IB(I))%NA(J), J=1, SimBox(IB(I))%NGROUP)
             write(hFile,FMT="(A, 3(1PE12.5,1x))") "&TEMPCAL  instant temperature(K):      ", SimBox(I)%Temperature

             !$$--- write out the configures
             write(hFile,FMT="(20A))")title(1:len_trim(title))
             do K=1, SimBox(IB(I))%NPRT
                PSTR = ""
                IS  = 0
                do ID=1, NDAT
                   if(IsInteger_DataPad(pDat(ID)) ) then
                      call Write_DataPad(pDat(ID), K, PSTR0, Fmt="5x,I8,5x")
                    else if(IsDouble_DataPad(pDat(ID)) ) then
                      call Write_DataPad(pDat(ID), K, PSTR0, Fmt="1PE17.8,1X")
                    end if
                    PSTR = PSTR(1:IS)//PSTR0(1:DIM(ID,2)*P_TL)
                    IS   = IS + DIM(ID,2)*P_TL
                end do

                 write(hFile,fmt=fmt)SimBox(IB(I))%ITYP(K), SimBox(IB(I))%XP(K,1:3)/SimBox(IB(I))%RR,           &
                                     SimBox(IB(I))%XP1(K,1:3)/SimBox(IB(I))%RR/CP_S2PS, SimBox(IB(I))%STATU(K), &
                                     SimBox(IB(I))%FP(K,1:3)*CP_ERGEV*SimBox(I)%RR,                             &
                                    -SimBox(IB(I))%EPOT(K)*CP_ERGEV,  SimBox(IB(I))%EKIN(K)*CP_ERGEV,           &
                                     SimBox(IB(I))%DIS(K,1:3)/SimBox(IB(I))%RR,                                 &
                                     PSTR(1:len_trim(PSTR))

             end do
             !write(hFile,FMT="(A, 100(I8,1x))")
          end do

       end if
       close(hFile)
       if(allocated(pDat)) deallocate(pDat)
       if(allocated(DIM))  deallocate(DIM)
       call Release_DataPad(tsDat)

  100  format(I8,3X,3(1PE17.8,1X),I6,2X,16(1PE17.8,1X))
       return
  end subroutine Putout_Instance_Config_SimBoxArray0
  !****************************************************************************

  !****************************************************************************
  subroutine Putout_Instance_Config_SimBoxArray1(Fhead, Multi, SimBox, Stamp)
  !***  PORPOSE: to putout the instant configuaration of a simulation box
  !
  !     INPUT:   Fhead, the head of the file name storing the configuration
  !              Multi, indicating the consigures to be output in one combined file or a set of seperated files  
  !              SimBox, the simulation boxes
  !              Stamp,  the recording stamp
  !
  implicit none
      !--- DUMMY variables
      character*(*),                   intent(in):: Fhead
      integer,                         intent(in):: Multi
      type(SimMDBox),dimension(:),     intent(in):: SimBox
      type(MDRecordStamp)                        :: Stamp
      !--- local variables
      character*256 ::FNAME
      type(MDRecordStamp)::tSTAMP
      integer::IB, IB0

      !---
        !--- to output the configures in a single file
        if(Multi .gt. 0) then
           call STRCATI(FNAME, Fhead, "_", Stamp%ITEST, 4)
           !$$--- to calculate the extension of the file
           call STRCATI(FNAME, FNAME, ".", Stamp%ICFG(1), 4)
           call Putout_Instance_Config_SimBoxArray0(FNAME, SimBox, Stamp) 
           return
        end if     
        !--- to output the configures in seperated files
        call Copy_RecordStamp(Stamp, tSTAMP)
        do IB= Stamp%IBox(1), Stamp%IBox(2)
           tSTAMP%IBox = IB
           call STRCATI(FNAME, Fhead, "_", IB, 4)
           !$$--- to calculate the extension of the file
           call STRCATI(FNAME, FNAME, ".", Stamp%ICFG(1), 4)

           IB0 = IB - Stamp%IBox(1) + 1
           call Putout_Instance_Config_SimBoxArray0(FNAME, SimBox(IB0:IB0), tSTAMP)
        end do

       return
  end subroutine Putout_Instance_Config_SimBoxArray1
  !****************************************************************************

  !****************************************************************************
  subroutine Putout_Instance_Config_SimBoxArray2(CtrlParam, SimBox, Stamp)
  !***  PORPOSE: to putout the instant configuaration of a simulation box
  !
  !     INPUT:   CtrlParam, 
  !              SimBox,   the simulation box
  !              Stamp,    the recording stamp
  !
  implicit none
      !--- DUMMY variables
      type(SimMDCtrl),              intent(in):: CtrlParam
      type(SimMDBox),dimension(:),  intent(in):: SimBox
      type(MDRecordStamp)                     :: Stamp
      !--- local variables
      character*256 ::FNAME
      !---
        call STRCATI(FNAME, CtrlParam%f_geometry, "P", 0, 4)
        call Putout_Instance_Config_SimBoxArray1(FNAME, CtrlParam%MultOutputG, SimBox, Stamp)
       return
  end subroutine Putout_Instance_Config_SimBoxArray2
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_Instance_Config_SimBoxArray0(Fname, SimBox, Stamp, ERR)
  !***  PORPOSE: to put in the instant configuaration of a simulation box
  !     INPUT:   Fname, the filname of configuration
  !              Stamp, the time stamp created in calling routine
  !
  !     OUTPUT:  SimBox, the box with configuration loaded from the file
  !              Stamp,  the stamp load from the file
  !              ERR, message indicating the statu
  !
  implicit none
      !--- DUMMY variables
      character*(*),               intent(in)   ::Fname
      type(SimMDBox), dimension(:),intent(inout)::SimBox
      type(MDRecordStamp),         intent(inout)::Stamp
      integer::ERR
      !--- local variables
      integer::hFile, LINE, I
      logical::opened
      character*256::str
      character*64::KEYWORD0, KEYWORD1
      type(SimMDBox)::SimBoxSwap

       !--- check if the file exit
       ERR = CP_FILE_EXIST
       inquire(FILE=fname, exist=opened)
       if(.not.opened) then
          ERR = CP_FILE_NOTEXIST
          write(*,*) "MDPSCU Warning: file not found ", fname(1:len_trim(fname))
          call ONWARNING(gm_OnWarning)
          return
       end if

       !--- to open the file
       write(*,*) "Load from ", fname(1:len_trim(fname))
       call  AvailableIOUnit(hFile)
       open(hFile, file = fname, status='old')
       LINE = 0
       call GetInputStrLine(hFile,STR, LINE, "!", *100)
       STR = adjustl(STR)
       call GetKeyWord("&", STR, KEYWORD0)
       call UpCase(KEYWORD0)

       select case(KEYWORD0(1:LEN_TRIM(KEYWORD0)))
              case (PKW_OUTCFG_FORMAT14, PKW_OUTCFG_FORMAT15, PKW_OUTCFG_FORMAT16)
                    rewind(hFile)
                    !$$--- check if the required data is available
                    if(NumberofData_DataPad(SimBox(1)%proKWDList) .gt. 0) then
                       do I=1, NumberofData_DataPad(SimBox(1)%proKWDList)
                          call Tag_DataPad(I, SimBox(1)%proKWDList, KEYWORD1)
                          if(.not.IsBasicKWD(KEYWORD1)) then
                              write(*,fmt="(A,A,A)") ' MDPSCU Warning: property data keywords:'
                              write(*,fmt="(A,A,A)") '                                   '//KEYWORD1(1:len_trim(KEYWORD1))
                              write(*,fmt="(A,A,A)") ' is unvailable from the configure file.'
                              call ONWARNING(gm_OnWarning)
                          end if
                       end do

                    end if
                    select case(KEYWORD0(1:LEN_TRIM(KEYWORD0)))
                           case (PKW_OUTCFG_FORMAT14)
                                call Putin_FORMAT14_Config_SimBoxArray(hFile, Stamp%ITime, Stamp%ISect, Stamp%Time, &
                                                                                     Stamp%ScalTime, SimBox)

                           case (PKW_OUTCFG_FORMAT15)
                               call Putin_FORMAT15_Config_SimBoxArray(hFile, Stamp%ITime, Stamp%ISect, Stamp%Time, &
                                                                                Stamp%ScalTime, SimBox)

                           case (PKW_OUTCFG_FORMAT16)
                               call Putin_FORMAT16_Config_SimBoxArray(hFile, SimBox, Stamp)
                    end select

              case (PKW_OUTCFG_FORMAT18)
                    rewind(hFile)
                    call Putin_FORMAT18_Config_SimBoxArray(hFile, SimBox, Stamp)

              case(PKW_OUTCFG_FORMATXYZ)
                   rewind(hFile)
                   call CopyInformation_SimMDBox(SimBox(1), SimBoxSwap)
                   SimBoxSwap%NPRT = SimBox(1)%NPRT*size(SimBox)
                   if(allocated(SimBox(1)%XP5) ) then
                      call Initialize_SimMDBox(SimBoxSwap,scheme=5)
                   else
                      call Initialize_SimMDBox(SimBoxSwap,scheme=2)
                   end if
                   !&&--- the XYZCOL is always needed
                   call AddDataProKWD_SimMDBox(SimBoxSwap, "XYZCOL")
                   call Read_ConfigByCol18_SimMDBox(SimBoxSwap, hFile, mod=1)
                   call Seperate_SimMDBox(SimBoxSwap, SimBox)
                   call Release_SimMDBox(SimBoxSwap)

              case default
                   if(len_trim(KEYWORD0) .gt. 0) then
                      write(*,fmt="(A,A)")        ' MDPSCU Error: unknown file format: ', KEYWORD0(1:len_trim(KEYWORD0))
                      write(*,fmt="(A,A, A, I5)") '               check file ', fname(1:len_trim(fname)), ' at line:', LINE
                      write(*,fmt="(A)")          '               process to be stopped'
                      stop
                   end if
                   rewind(hFile)
                   call Putin_FORMAT18_Config_SimBoxArray(hFile, SimBox, Stamp)

       end select
       close(hFile)
       return

 100   close(hFile)
       write(*,*)"MDPSCU Error: error in putin configuration"
       write(*,*)"The process to be stopped."
       stop

       return
  end subroutine Putin_Instance_Config_SimBoxArray0
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_Instance_Config_SimBoxArray1(Fhead, Multi, SimBox, Stamp, Err)
  !***  PORPOSE: to putin the instant configuaration of a simulation box
  !
  !     INPUT:   Fhead, the head of the file name storing the configuration
  !              Multi,  indicating the consigures to be output in one combined file or a set of seperated files  
  !              SimBox, the simulation boxes
  !              Stamp,  the recording stamp
  !
  implicit none
      !--- DUMMY variables
      character*(*),               intent(in):: Fhead
      integer,                     intent(in):: Multi
      type(SimMDBox),dimension(:)            :: SimBox
      type(MDRecordStamp)                    :: Stamp
      integer                                :: Err     
      !--- local variables
      character*256 ::FNAME
      type(MDRecordStamp)::tSTAMP
      integer::IB, IB0, IBOX(2)

      !---
        !--- to output the configures in a single file
        if(Multi .gt. 0) then
           call STRCATI(FNAME, Fhead, "_", Stamp%ITEST, 4)
           !$$--- to calculate the extension of the file
           call STRCATI(FNAME, FNAME, ".", Stamp%ICFG(1), 4)
           call Putin_Instance_Config_SimBoxArray0(FNAME, SimBox, Stamp, ERR) 
           return
        end if     
        !--- to putin the configures in seperated files
        call Copy_RecordStamp(Stamp, tSTAMP)
        IBOX = Stamp%IBox
        do IB= IBOX(1), IBOX(2)
           tSTAMP%IBox = IB
           call STRCATI(FNAME, Fhead, "_", IB, 4)
           call STRCATI(FNAME, FNAME, ".", Stamp%ICFG(1), 4)
           IB0 = IB - IBOX(1) + 1
           call Putin_Instance_Config_SimBoxArray0(FNAME, SimBox(IB0:IB0), tSTAMP,Err)
           call Copy_RecordStamp(tSTAMP, Stamp)
        end do
        Stamp%IBox =IBOX

       return
  end subroutine Putin_Instance_Config_SimBoxArray1
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_Instance_Config_SimBoxArray2(CtrlParam, SimBox, Stamp, Err)
  !***  PORPOSE: to putout the instant configuaration of a simulation box
  !
  !     INPUT:   CtrlParam, 
  !              SimBox,   the simulation box
  !              Stamp,    the recording stamp
  !
  implicit none
      !--- DUMMY variables
      type(SimMDCtrl),            intent(in):: CtrlParam
      type(SimMDBox),dimension(:)           :: SimBox
      type(MDRecordStamp)                   :: Stamp
      integer                               :: Err
      !--- local variables
      character*256 ::FNAME
      !---
        call STRCATI(FNAME, CtrlParam%f_geometry, "P", 0, 4)
        call Putin_Instance_Config_SimBoxArray1(FNAME, CtrlParam%MultOutputG, SimBox, Stamp, Err)
       return
  end subroutine Putin_Instance_Config_SimBoxArray2
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_FORMAT14_Config_SimBoxArray(hFile, ITIME, ISECT, TIME, SCALTIME, SimBox)
  !***  PORPOSE: to put in the instant configuaration of a simulation box
  !     INPUT:   hFile,   the IO unit
  !              ITIME,   the ith time step
  !              ISECT,   the time section ID
  !
  !
  !     OUTPUT:  ITIME,   the ith time step
  !              ISECT,   the time section ID
  !              TIME,    the time
  !              SCATIME, the scaling time, that could generated by PARREP
  !              SimBox,  the simulation box
  !
  implicit none
      !--- DUMMY variables
      integer,                    intent(in)   ::hFile, ITIME, ISECT
      real(KINDDF),               intent(inout)::TIME, SCALTIME
      type(SimMDBox), dimension(:),intent(inout)::SimBox
      !--- local variables
      integer::K, I, LINE, ITIME0, ISECT0, N, NATOM
      character*256::str
      character*32::KEYWORD0, KEYWORD, STRNUMB(4)
      character*1::CC

       LINE = 0
       NATOM = 0
       CC    = ""
       call GetInputStrLine(hFile,STR, LINE, "!",*100)
       STR = adjustl(STR)
       call GetKeyWord("&", STR, KEYWORD0)
       call UpCase(KEYWORD0)
        if( KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq. PKW_OUTCFG_FORMAT14) then
            !--- check the consistent of data
            !--- and get time information
            do while(.true.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&TIMESTEPS") then
                   call Extract_Numb(STR,3,n,STRNUMB)
                   ITIME0 = ISTR(STRNUMB(1))
                   if(ITIME0 .NE. ITIME) then
                      write(*,fmt="(A)")" MDPSCU Warning: time step in putin configure is not consistent "
                      write(*,fmt="(A, I8, A, I8)")"                 ITIME0", ITIME0, " vs ITIME ", ITIME
                      call ONWARNING(gm_OnWarning)
                   end if
                   TIME     = DRSTR(STRNUMB(2))
                   SCALTIME = TIME
                   ISECT0   = ISTR(STRNUMB(3))
                   if(ISECT0 .NE. ISECT) then
                      write(*,fmt="(A)")" MDPSCU Warning: time section in putin configure is not consistent "
                      write(*,fmt="(A,I8,A,I8)")"                 ISECT0", ISECT0, " vs ISECT ", ISECT
                      call ONWARNING(gm_OnWarning)
                   end if
                else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&NATOM") then
                     !--- should be changed in future
                      call Extract_Numb(STR,3,N,STRNUMB)
                      N = ISTR(STRNUMB(2))
                      if(N .ne. SimBox(1)%NPRT .and. NATOM.le.0) then
                         write(*,fmt="(A,I9)") ' MDPSCU Warning: number of atoms in the configure file is: ', N
                         write(*,fmt="(A,I9)") '                 while, the number of atoms in a box is: .',  SimBox(1)%NPRT
                         call ONWARNING(gm_OnWarning)
                         CC = 'C'
                         write(*,fmt="(A, I9)")    '                 the number of atoms in boxes are reset to ',   N
                         if(CC .eq. 'c' .or. CC .eq. 'C') then
                            NATOM = N
                            if(NATOM .gt. Simbox(1)%NPRT) then !-- reallocate memory is needed
                               do I=1, SIZE(SimBox)
                                  call AddAtoms_SimMDBox(Simbox(I), NATOM-Simbox(I)%NPRT, ITYPE=1, TYPEORDER=1)
                               end do
                            else  !-- reallocate memory is not needed, just set the number of atoms as NATOM
                               do I=1, SIZE(SimBox)
                                  Simbox(I)%NPRT = NATOM
                               end do
                            end if

                         else
                            stop ' Process stopped by user'
                         end if
                      end if
                      exit
                else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&TYPE") then
                    exit
                end if
            end do
            rewind(hFile)
            LINE = 0
        else
            !--- no header found, rewind to the top
            rewind(hFile)
            !--- skip over the header if there is
            do I=1, LINE-1
               read(hFile,fmt="(A256)") STR
            end do
        end if

        do I=1, size(SimBox)
           if( KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq. PKW_OUTCFG_FORMAT14) then
             !--- skip over the header if there is
              do while(.true.)  !skip the header
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 STR = adjustl(STR)
                 call GetKeyWord("&", STR, KEYWORD)
                 call UpCase(KEYWORD)
                 if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&TYPE") exit
              end do
           end if
           do K=1, SimBox(I)%NPRT
              LINE = LINE + 1
              read(hFile,*, ERR=100) SimBox(I)%ITYP(K), SimBox(I)%XP(K,1:3), SimBox(I)%XP1(K,1:3), SimBox(I)%STATU(K), SimBox(I)%fP(K,1:3), &
                            SimBox(I)%EPOT(K), SimBox(I)%EKIN(K), SimBox(I)%DIS(K,1:3)
           end do

         !$$--- change the UNIT back to CGS
          SimBox(I)%XP(1:SimBox(I)%NPRT,1:3) = SimBox(I)%XP(1:SimBox(I)%NPRT,1:3)*SimBox(I)%RR
          SimBox(I)%DIS(1:SimBox(I)%NPRT,1:3)= SimBox(I)%DIS(1:SimBox(I)%NPRT,1:3)*SimBox(I)%RR
          SimBox(I)%EPOT(1:SimBox(I)%NPRT)   = -SimBox(I)%EPOT(1:SimBox(I)%NPRT)/CP_ERGEV
          SimBox(I)%EKIN(1:SimBox(I)%NPRT)   = SimBox(I)%EKIN(1:SimBox(I)%NPRT)/CP_ERGEV

         !$$--- PART in boxes have been change, we should check if the group information consistent
          if(CC .eq. 'c' .or. CC .eq. 'C' .and. I.eq.1) then
             call Check_Config_Consistent_SimMDBox(Simbox(1), 3)
             do k=2, SIZE(SimBox)
                call CopyInformation_SimMDBox(Simbox(1), Simbox(k))
                Simbox(k)%NA = Simbox(1)%NA
             end do
          end if
       end do
       return

  100  write(*,fmt="(A,I9)") ' MDPSCU Error: error in reading configuration from a file at line ', LINE
       write(*,fmt="(A)")    '               Process to be stopped'
       stop
  end subroutine Putin_FORMAT14_Config_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_FORMAT15_Config_SimBoxArray(hFile, ITIME, ISECT, TIME, SCALTIME, SimBox)
  !***  PORPOSE: to put in the instant configuaration of a simulation box
  !     INPUT:   hFile,   the IO unit
  !              ITIME,   the ith time step
  !              ISECT,   the time section ID
  !
  !
  !     OUTPUT:  ITIME,   the ith time step
  !              ISECT,   the time section ID
  !              TIME,    the time
  !              SCATIME, the scaling time, that could generated by PARREP
  !              SimBox,  the simulation box
  !
  implicit none
      !--- DUMMY variables
      integer,                     intent(in)   ::hFile, ITIME, ISECT
      real(KINDDF),                intent(inout)::TIME, SCALTIME
      type(SimMDBox), dimension(:),intent(inout)::SimBox
      !--- local variables
      integer::K, I, LINE, ITIME0, ISECT0, N, NATOM
      character*256::str
      character*32::KEYWORD0, KEYWORD, STRNUMB(4)
      character*1::CC
      real*8::ENERGY

       LINE = 0
       NATOM = 0
       CC    = ""
       call GetInputStrLine(hFile,STR, LINE, "!",*100)
       STR = adjustl(STR)
       call GetKeyWord("&", STR, KEYWORD0)
       call UpCase(KEYWORD0)
        if( KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq. PKW_OUTCFG_FORMAT15) then
            !--- check the consistent of data
            !--- and get time information
            do while(.true.)
                call GetInputStrLine(hFile,STR, LINE, "!",*100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&TIMESTEPS") then
                   call Extract_Numb(STR,3,n,STRNUMB)
                   ITIME0 = ISTR(STRNUMB(1))
                   if(ITIME0 .NE. ITIME) then
                      write(*,fmt="(A)")" MDPSCU Warning: time step in putin configure is not consistent "
                      write(*,fmt="(A, I8, A, I8)")"                 ITIME0", ITIME0, " vs ITIME ", ITIME
                      call ONWARNING(gm_OnWarning)
                   end if
                   TIME     = DRSTR(STRNUMB(2))
                   ISECT0   = ISTR(STRNUMB(3))
                   if(ISECT0 .NE. ISECT) then
                      write(*,fmt="(A)")" MDPSCU Warning: time section in putin configure is not consistent "
                      write(*,fmt="(A,I8,A,I8)")"                 ISECT0", ISECT0, " vs ISECT ", ISECT
                      call ONWARNING(gm_OnWarning)
                   end if

                else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&SCALTIME") then
                     call Extract_Numb(STR,3,n,STRNUMB)
                     SCALTIME     = DRSTR(STRNUMB(2))

                else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&NATOM") then
                     !--- should be changed in future
                      call Extract_Numb(STR,3,N,STRNUMB)
                      N = ISTR(STRNUMB(2))
                      if(N .ne. SimBox(1)%NPRT .and. NATOM.le.0) then
                         write(*,fmt="(A,I9)") ' MDPSCU Warning: number of atoms in the configure file is: ', N
                         write(*,fmt="(A,I9)") '                 while, the number of atoms in a box is: .',  SimBox(1)%NPRT
                         call ONWARNING(gm_OnWarning)
                         CC = 'C'
                         write(*,fmt="(A, I9)")    '                 the number of atoms in boxes are reset to ',   N
                         if(CC .eq. 'c' .or. CC .eq. 'C') then
                            NATOM = N
                            if(NATOM .gt. Simbox(1)%NPRT) then !-- reallocate memory is needed
                               do I=1, SIZE(SimBox)
                                  call AddAtoms_SimMDBox(Simbox(I), NATOM-Simbox(I)%NPRT, ITYPE=1, TYPEORDER=1)
                               end do
                            else  !-- reallocate memory is not needed, just set the number of atoms as NATOM
                               do I=1, SIZE(SimBox)
                                  Simbox(I)%NPRT = NATOM
                               end do
                            end if

                         else
                            stop ' Process stopped by user'
                         end if
                      end if
                      exit
                else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&TYPE") then
                    exit
                end if
            end do
            rewind(hFile)
            LINE = 0
        else
            !--- no header found, rewind to the top
            rewind(hFile)
            !--- skip over the header if there is
            do I=1, LINE-1
               read(hFile,fmt="(A256)") STR
            end do
        end if

        do I=1, SIZE(SimBox)
           if( KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq. PKW_OUTCFG_FORMAT15) then
             !--- skip over the header if there is
              do while(.true.)  !skip the header
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 STR = adjustl(STR)
                 call GetKeyWord("&", STR, KEYWORD)
                 call UpCase(KEYWORD)
                 if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&TYPE") exit
              end do
           end if
           do K=1, SimBox(I)%NPRT
              LINE = LINE + 1
              read(hFile,*, ERR=100) SimBox(I)%ITYP(K), SimBox(I)%XP(K,1:3), SimBox(I)%XP1(K,1:3), SimBox(I)%STATU(K), SimBox(I)%fP(K,1:3), &
                            SimBox(I)%EPOT(K), SimBox(I)%EKIN(K), ENERGY, SimBox(I)%DIS(K,1:3)
           end do

         !$$--- change the UNIT back to CGS
          SimBox(I)%XP(1:SimBox(I)%NPRT,1:3) = SimBox(I)%XP(1:SimBox(I)%NPRT,1:3)*SimBox(I)%RR
          SimBox(I)%DIS(1:SimBox(I)%NPRT,1:3)= SimBox(I)%DIS(1:SimBox(I)%NPRT,1:3)*SimBox(I)%RR
          SimBox(I)%EPOT(1:SimBox(I)%NPRT)   = -SimBox(I)%EPOT(1:SimBox(I)%NPRT)/CP_ERGEV
          SimBox(I)%EKIN(1:SimBox(I)%NPRT)   = SimBox(I)%EKIN(1:SimBox(I)%NPRT)/CP_ERGEV

         !$$--- PART in boxes have been change, we should check if the group information consistent
          if(CC .eq. 'c' .or. CC .eq. 'C' .and. I.eq.1) then
             call Check_Config_Consistent_SimMDBox(Simbox(1), 3)
             do k=2, SIZE(SimBox)
                call CopyInformation_SimMDBox(Simbox(1), Simbox(k))
                Simbox(k)%NA = Simbox(1)%NA
             end do
          end if
       end do
       return

  100  write(*,fmt="(A,I9)") ' MDPSCU Error: error in reading configuration from a file at line ', LINE
       write(*,fmt="(A)")    '               Process to be stopped'
       stop
  end subroutine Putin_FORMAT15_Config_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_FORMAT16_Config_SimBoxArray(hFile, SimBox, Stamp)
  !***  PORPOSE: to put in the instant configuaration of a simulation box
  !     INPUT:   hFile, the IO unit
  !              Stamp, the Stamp
  !
  !     OUTPUT:  SimBox, the simulation box
  !              Stamp,  the stamp updated
  !
  implicit none
      !--- DUMMY variables
      integer, intent(in)::hFile
      type(SimMDBox),dimension(:),intent(inout)::SimBox
      type(MDRecordStamp),        intent(inout)::Stamp

      !--- local variables
      integer::K, I, LINE, ITEST0, IBOX0(2), ICFG0(2), ITIME0, ISECT0, N, NATOM
      character*256::str
      character*32::KEYWORD0, KEYWORD, STRNUMB(4)
      character*1::CC
      real*8::ENERGY

       !---

       LINE = 0
       NATOM = 0
       CC    = ""
       call GetInputStrLine(hFile,STR, LINE, "!",*100)
       STR = adjustl(STR)
       call GetKeyWord("&", STR, KEYWORD0)
       call UpCase(KEYWORD0)
       if( KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq. PKW_OUTCFG_FORMAT16) then
            !--- check the consistent of data
            !--- and get stamp information
            do while(.true.)
                call GetInputStrLine(hFile,STR, LINE, "!",*100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case(KEYWORD(1:LEN_TRIM(KEYWORD)) )
                       case("&APPTYPE")
                            call Extract_Substr(STR,1,n,STRNUMB)
                            Stamp%AppType = STRNUMB(1)(1:len_trim(STRNUMB(1)))

                       case("&TESTID")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            ITEST0 = ISTR(STRNUMB(1))
                            if(Stamp%ITest .ge. C_IZERO) then
                               if(ITEST0 .ne. Stamp%ITest) then
                                   write(*,fmt="(A)")           " MDPSCU Warning: Job ID in putin configure is not consistent "
                                   write(*,fmt="(A, I8, A, I8)")"                 ITest0", ITEST0, " vs ITest ", Stamp%ITest
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ITest = ITEST0

                       case("&BOXID")
                            call Extract_Numb(STR,2,n,STRNUMB)
                            IBOX0(1) = ISTR(STRNUMB(1))
                            IBOX0(2) = ISTR(STRNUMB(2))
                            if(any(Stamp%IBox(1:2) .ge. 0)) then
                               if(any(IBOX0(1:2) .ne. Stamp%IBox(1:2))) then
                                   write(*,fmt="(A)")           " MDPSCU Warning: Box ID in putin configure is not consistent "
                                   write(*,fmt="(A, 2(I8,1x), A, 2(I8,1x))")"                 IBox0", IBOX0, " vs IBox ", Stamp%IBox
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%IBox = IBOX0

                       case("&CFGID")
                            call Extract_Numb(STR,2,n,STRNUMB)
                            ICFG0(1) = ISTR(STRNUMB(1))
                            ICFG0(2) = ISTR(STRNUMB(2))
                            if(any(Stamp%ICfg(1:2).gt.0)) then
                               if(any(ICFG0(1:2) .ne. Stamp%ICfg(1:2))) then
                                   write(*,fmt="(A)")           " MDPSCU Warning: Configure ID in putin configure is not consistent "
                                   write(*,fmt="(A, 2(I8,1x), A, 2(I8,1x))")"                 ICfg0", ICfg0, " vs ICfg ", Stamp%ICfg
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ICfg = ICFG0

                       case("&TIMESTEPS")
                            call Extract_Numb(STR,3,n,STRNUMB)
                            ITIME0 = ISTR(STRNUMB(1))
                            if(Stamp%ITime .ge. 0) then
                               if(ITIME0 .ne. Stamp%ITime) then
                                   write(*,fmt="(A)")" MDPSCU Warning: time step in putin configure is not consistent "
                                   write(*,fmt="(A, I8, A, I8)")"                 ITIME0", ITIME0, " vs ITIME ", Stamp%ITime
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ITime = ITIME0

                       case("&TIME")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%Time = DRSTR(STRNUMB(1))

                       case("&TIMESECT")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            ISECT0   = ISTR(STRNUMB(1))
                            if(Stamp%ISect .gt. 0) then
                               if(ISECT0 .ne. Stamp%ISect) then
                                  write(*,fmt="(A)")" MDPSCU Warning: time section in putin configure is not consistent "
                                  write(*,fmt="(A,I8,A,I8)")"                 ISECT0", ISECT0, " vs ISECT ", Stamp%ISect
                                  call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ISect = ISECT0

                       case("&SCALTIME")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%ScalTime = DRSTR(STRNUMB(1))

                       case("&INSTANTTEMP", "&TEMPSET")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%InstantTemp = DRSTR(STRNUMB(1))

                       case("&NATOM")
                            call Extract_Numb(STR,1,N,STRNUMB)
                            N = ISTR(STRNUMB(1))
                            if(N .ne. SimBox(1)%NPRT .and. NATOM.le.0) then
                               write(*,fmt="(A,I9)") ' MDPSCU Warning: number of atoms in the configure file is: ', N
                               write(*,fmt="(A,I9)") '                 while, the number of atoms in a box is: .',  SimBox(1)%NPRT
                               call ONWARNING(gm_OnWarning)
                               !--- to be change this default seting?
                               CC = 'C'
                               write(*,fmt="(A, I9)")    '                 the number of atoms in boxes are reset to ',   N
                               if(CC .eq. 'c' .or. CC .eq. 'C') then
                                  NATOM = N
                                  if(NATOM .gt. Simbox(1)%NPRT) then !-- reallocate memory is needed
                                     do I=1, SIZE(SimBox)
                                        call AddAtoms_SimMDBox(Simbox(I), NATOM-Simbox(I)%NPRT, ITYPE=1, TYPEORDER=1)
                                     end do
                                  else  !-- reallocate memory is not needed, just set the number of atoms as NATOM
                                     do I=1, SIZE(SimBox)
                                        Simbox(I)%NPRT = NATOM
                                      end do
                                  end if

                               else
                                  stop ' Process stopped by user'
                               end if
                            end if

                       case("&TYPE")
                            exit
                end select
            end do
            rewind(hFile)
            LINE = 0
        else
            !--- no header found, rewind to the top
            rewind(hFile)
            !--- skip over the header if there is
            do I=1, LINE-1
               read(hFile,fmt="(A256)") STR
            end do
        end if

        do I=1, SIZE(SimBox)
           if( KEYWORD0(1:LEN_TRIM(KEYWORD0)) .eq. PKW_OUTCFG_FORMAT16) then
             !--- skip over the header if there is
              do while(.true.)  !skip the header
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 STR = adjustl(STR)
                 call GetKeyWord("&", STR, KEYWORD)
                 call UpCase(KEYWORD)
                 if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&REALTEMP" .or. &
                    KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&TEMPCAL")  then
                    call Extract_Numb(STR,1,n,STRNUMB)
                    SimBox(I)%Temperature = DRSTR(STRNUMB(1))
                 else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&TYPE") then
                    exit
                 end if
              end do
           end if
           do K=1, SimBox(I)%NPRT
              LINE = LINE + 1
              read(hFile,*, ERR=100) SimBox(I)%ITYP(K), SimBox(I)%XP(K,1:3), SimBox(I)%XP1(K,1:3), SimBox(I)%STATU(K), SimBox(I)%fP(K,1:3), &
                            SimBox(I)%EPOT(K), SimBox(I)%EKIN(K), ENERGY, SimBox(I)%DIS(K,1:3)
           end do

         !$$--- change the UNIT back to CGS
          SimBox(I)%XP(1:SimBox(I)%NPRT,1:3) = SimBox(I)%XP(1:SimBox(I)%NPRT,1:3)*SimBox(I)%RR
          SimBox(I)%DIS(1:SimBox(I)%NPRT,1:3)= SimBox(I)%DIS(1:SimBox(I)%NPRT,1:3)*SimBox(I)%RR
          SimBox(I)%EPOT(1:SimBox(I)%NPRT)   = -SimBox(I)%EPOT(1:SimBox(I)%NPRT)/CP_ERGEV
          SimBox(I)%EKIN(1:SimBox(I)%NPRT)   = SimBox(I)%EKIN(1:SimBox(I)%NPRT)/CP_ERGEV

         !$$--- PART in boxes have been change, we should check if the group information consistent
          if(CC .eq. 'c' .or. CC .eq. 'C' .and. I.eq.1) then
             call Check_Config_Consistent_SimMDBox(Simbox(1), 3)
             do k=2, SIZE(SimBox)
                call CopyInformation_SimMDBox(Simbox(1), Simbox(k))
                Simbox(k)%NA = Simbox(1)%NA
             end do
          end if
       end do
       return

  100  write(*,fmt="(A,I9)") ' MDPSCU Error: error in reading configuration from a file at line ', LINE
       write(*,fmt="(A)")    '               Process to be stopped'
       stop
  end subroutine Putin_FORMAT16_Config_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Putin_FORMAT18_Config_SimBoxArray(hFile, SimBox, Stamp)
  !***  PORPOSE: to put in the instant configuaration of a simulation box
  !              The differece from FOMAT16 is that the DataPad to be read in
  !              utomatically.
  !
  !     INPUT:   hFile, the IO unit
  !              Stamp, the Stamp
  !
  !     OUTPUT:  SimBox, the simulation box
  !              Stamp,  the stamp updated
  !
  implicit none
      !--- DUMMY variables
      integer,                     intent(in)   ::hFile
      type(SimMDBox), dimension(:),intent(inout)::SimBox
      type(MDRecordStamp),         intent(inout)::Stamp

      !--- local variables
      integer::K, I, LINE, ITEST0, IBOX0(2), ICFG0(2), ITIME0, ISECT0, N, NATOM, EXT
      real(KINDDF)::TIME0
      character*256::str
      character*32::KEYWORD, STRNUMB(3)
      character*1::CC, TFMT
      type(DataPad), pointer::TSDat
      !--- local variables for colume
      integer:: NCOL, XPCOL, XP1COL, FPCOL, TYPCOL, POTCOL, KENCOL, STATCOL, DISCOL
      integer:: NDAT
      integer,  dimension(:), allocatable::PDatCol


         LINE = 0
         NATOM = 0
         CC    = ""

            !--- check the consistent of data
            !--- and get stamp information
            do while(.true.)
                call GetInputStrLine(hFile,STR, LINE, "!",*100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                if(len_trim(KEYWORD) .le. 0) exit
                !--- deal witht the keyword
                call UpCase(KEYWORD)
                select case(KEYWORD(1:LEN_TRIM(KEYWORD)) )
                       case("&APPTYPE")
                            call Extract_Substr(STR,1,n,STRNUMB)
                            Stamp%AppType = STRNUMB(1)(1:len_trim(STRNUMB(1)))

                       case("&TESTID")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            ITEST0 = ISTR(STRNUMB(1))
                            if(Stamp%ITest .ge. C_IZERO) then
                               if(ITEST0 .ne. Stamp%ITest) then
                                   write(*,fmt="(A)")           " MDPSCU Warning: Job ID in putin configure is not consistent "
                                   write(*,fmt="(A, I8, A, I8)")"                 ITest0", ITEST0, " vs ITest ", Stamp%ITest
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ITest = ITEST0

                       case("&BOXID")
                            call Extract_Numb(STR,2,n,STRNUMB)
                            IBOX0(1) = ISTR(STRNUMB(1))
                            IBOX0(2) = ISTR(STRNUMB(2))
                            if(any(Stamp%IBox(1:2) .ge. 0)) then
                               if(any(IBOX0(1:2) .ne. Stamp%IBox(1:2))) then
                                   write(*,fmt="(A)")           " MDPSCU Warning: Box ID in putin configure is not consistent "
                                   write(*,fmt="(A, 2(I8,1x), A, 2(I8,1x))")"                 IBox0", IBOX0, " vs IBox ", Stamp%IBox
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%IBox = IBOX0

                       case("&CFGID")
                            call Extract_Numb(STR,2,n,STRNUMB)
                            ICFG0(1) = ISTR(STRNUMB(1))
                            ICFG0(2) = ISTR(STRNUMB(2))
                            if(any(Stamp%ICfg(1:2).gt.0)) then
                               if(any(ICFG0(1:2) .ne. Stamp%ICfg(1:2))) then
                                   write(*,fmt="(A)")           " MDPSCU Warning: Configure ID in putin configure is not consistent "
                                   write(*,fmt="(A, 2(I8,1x), A, 2(I8,1x))")"                 ICfg0", ICfg0, " vs ICfg ", Stamp%ICfg
                                   call OnWarning(gm_OnWarning)
                               end if
                            end if
                            Stamp%ICfg = ICFG0

                       case("&TIMESTEPS")
                            call Extract_Numb(STR,3,n,STRNUMB)
                            ITIME0 = ISTR(STRNUMB(1))
                            if(Stamp%ITime .ge. 0) then
                               if(ITIME0 .ne. Stamp%ITime) then
                                   write(*,fmt="(A)")" MDPSCU Warning: time step in putin configure is not consistent "
                                   write(*,fmt="(A, I8, A, I8)")"                 ITIME0", ITIME0, " vs ITIME ", Stamp%ITime
                                   call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ITime = ITIME0

                       case("&TIME")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            TIME0 =  DRSTR(STRNUMB(1))
                            if(TIME0 .lt. Stamp%Time) then
                                  write(*,fmt="(A)")" MDPSCU Warning: time in putin configure is not consistent "
                                  write(*,fmt="(A, 1PE12.4, A, 1PE12.4)")"                  TIME0", TIME0, " smaller than TIME ", Stamp%Time
                                  call ONWARNING(gm_OnWarning)
                            end if
                            Stamp%Time = TIME0

                       case("&TIMESECT")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            ISECT0   = ISTR(STRNUMB(1))
                            if(Stamp%ISect .gt. 0) then
                               if(ISECT0 .ne. Stamp%ISect) then
                                  write(*,fmt="(A)")" MDPSCU Warning: time section in putin configure is not consistent "
                                  write(*,fmt="(A,I8,A,I8)")"                 ISECT0", ISECT0, " vs ISECT ", Stamp%ISect
                                  call ONWARNING(gm_OnWarning)
                               end if
                            end if
                            Stamp%ISect = ISECT0

                       case("&SCALTIME")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            TIME0 =  DRSTR(STRNUMB(1))
                            if(TIME0 .lt. Stamp%ScalTime) then
                                  write(*,fmt="(A)")" MDPSCU Warning: scaled time in putin configure is not consistent "
                                  write(*,fmt="(A, I8, A, I8)")"                  TIME0", TIME0, " smaller than ", Stamp%ScalTime
                                  call ONWARNING(gm_OnWarning)
                            end if
                            Stamp%ScalTime = TIME0

                       case("&INSTANTTEMP", "&TEMPSET")
                            call Extract_Numb(STR,1,n,STRNUMB)
                            Stamp%InstantTemp = DRSTR(STRNUMB(1))

                       case("&NATOM")
                            call Extract_Numb(STR,1,N,STRNUMB)
                            N = ISTR(STRNUMB(1))
                            if(N .ne. SimBox(1)%NPRT .and. NATOM.le.0) then
                               write(*,fmt="(A,I9)") ' MDPSCU Warning: number of atoms in the configure file is: ', N
                               write(*,fmt="(A,I9)") '                 while, the number of atoms in a box is: .',  SimBox(1)%NPRT
                               call ONWARNING(gm_OnWarning)
                               !--- to be change this default seting?
                               CC = 'C'
                               write(*,fmt="(A, I9)")    '                 the number of atoms in boxes are reset to ',   N
                               if(CC .eq. 'c' .or. CC .eq. 'C') then
                                  NATOM = N
                                  if(NATOM .gt. Simbox(1)%NPRT) then !-- reallocate memory is needed
                                     do I=1, SIZE(SimBox)
                                        call AddAtoms_SimMDBox(Simbox(I), NATOM-Simbox(I)%NPRT, ITYPE=1, TYPEORDER=1)
                                     end do
                                  else  !-- reallocate memory is not needed, just set the number of atoms as NATOM
                                     do I=1, SIZE(SimBox)
                                        Simbox(I)%NPRT = NATOM
                                      end do
                                  end if

                               else
                                  stop ' Process stopped by user'
                               end if
                            end if

                       case default
                            !$$--- add the colume tag to a list. We do not add the data pad here,
                            !$$    considering the possibility thta NATOM could change the size.
                            if(IsBasicKWD(KEYWORD)) cycle
                            if(SameFirstWord(KEYWORD(len_trim(KEYWORD)-2:len_trim(KEYWORD)), 'COL')) then
                               call New_DataPad(KEYWORD(2:len_trim(KEYWORD)-3), TSDat)
                            end if
                end select
            end do

            !$$-- to check if the data we needed are available
             do I=1, NumberofData_DataPad(Simbox(1)%ptrDatPad)
               call Tag_DataPad(I, Simbox(1)%ptrDatPad, KEYWORD)
               call UpCase(KEYWORD)
               if(.not.HasTag_DataPad(KEYWORD, TSDat)) then
                  write(*,fmt="(A)") ' MDPSCU Warning: cannot find  keyword &'//KEYWORD(1:len_trim(KEYWORD))// &
                                     ' in the configure file.'
                  write(*,fmt="(A)") '                  which is needed for loading data '//KEYWORD(1:len_trim(KEYWORD))
                  call OnWarning(gm_OnWarning)
               end if
             end do

             !$$-- Now, we add the available extented data to the data pad
             rewind(hFile)
             LINE = 0
             do while(.true.)
                call GetInputStrLine(hFile,STR, LINE, "!",*100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                if(len_trim(KEYWORD) .le. 0) exit
                !$$--- check if it is a data tag
                if(SameFirstWord(KEYWORD(len_trim(KEYWORD)-2:len_trim(KEYWORD)), 'COL')) then
                   !$$--- in TSDat, the tag is upcased, and HasTag_DataPad is case
                   !$$    insensitive, here we keep the original case for title of
                   !$$    the data pad
                   if(HasTag_DataPad(KEYWORD(2:len_trim(KEYWORD)-3), TSDat) ) then
                      call Extract_Numb(STR, 2, N, STRNUMB)
                      EXT =ISTR(STRNUMB(2))
                      call Extract_Substr(STR, 1, N, STRNUMB)
                      TFMT = STRNUMB(1)(1:1)
                      call GetFirstWord(STR, KEYWORD)
                      do I = 1, size(Simbox)
                         call NewDatPad_SimMDBox(Simbox(I), KEYWORD(2:len_trim(KEYWORD)-3), TFMT, EXT)
                     end do
                   end if
                end if
             end do
             call Release_DataPad(tsDat)

             !$$--- to prepair the colume information
             NDAT = NumberofData_DataPad(SimBox(1)%ptrDatPad)
             if(NDAT .gt. 0) allocate(PDatCol(NDAT*2))
             !$$---
             LINE = 0
             rewind(hFile)
             call Read_ConfigByCol18A_SimMDBox(Simbox(1), hFile, LINE, NCOL, TYPCOL, XPCOL, XP1COL,  &
                                                    STATCOL, FPCOL, POTCOL, KENCOL, DISCOL,               &
                                                    PDatCol)

            rewind(hFile)
            LINE = 0
            do I=1, size(SimBox)
                !$$--- skip the header strings
                do while(.true.)
                   call GetInputStrLine(hFile,STR, LINE, "!",*100)
                    STR = adjustl(STR)
                    call GetKeyWord("&", STR, KEYWORD)
                    if(len_trim(KEYWORD) .le. 0) then
                       backspace(hFile)
                       exit
                    end if

                  if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&REALTEMP" .or. &
                     KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&TEMPCAL")  then
                     call Extract_Numb(STR,1,n,STRNUMB)
                     SimBox(I)%Temperature = DRSTR(STRNUMB(1))
                 end if
                end do

                call Read_ConfigByCol18B_SimMDBox(Simbox(I), hFile, NCOL, TYPCOL, XPCOL, XP1COL,  &
                                                       STATCOL, FPCOL, POTCOL, KENCOL, DISCOL,         &
                                                       PDatCol, LUNIT=1.D0)

                !$$--- change the unit back to CGS in CFGFMT18
                !$$    NOTE: because LUNIT on calling Read_ConfigByCol18B_SimMDBox is 1
                !$$    we should change the lenght unit to CGS.
                !$$    See also: the note in subroutine Load_Config_SimMDBox
                SimBox(I)%XP    = SimBox(I)%XP*SimBox(I)%RR
                SimBox(I)%DIS   = SimBox(I)%DIS*SimBox(I)%RR
                SimBox(I)%XP1   = SimBox(I)%XP1*SimBox(I)%RR*CP_S2PS
                SimBox(I)%FP    = SimBox(I)%FP*CP_EVERG/SimBox(I)%RR
                SimBox(I)%EPOT  = SimBox(I)%EPOT*CP_EVERG
                SimBox(I)%EKIN  = SimBox(I)%EKIN*CP_EVERG

               !$$--- PART in boxes have been changed, we should check if the group information consistent
               if(CC .eq. 'c' .or. CC .eq. 'C' .and. I.eq.1) then
                  call Check_Config_Consistent_SimMDBox(Simbox(1), 3)
                  do K=2, size(SimBox)
                     call CopyInformation_SimMDBox(Simbox(1), Simbox(K))
                     Simbox(k)%NA = Simbox(1)%NA
                  end do
               end if
            end do
       return

  100  write(*,fmt="(A,I9)") ' MDPSCU Error: error in reading configuration from a file at line ', LINE
       write(*,fmt="(A)")    '               Process to be stopped'
       stop
  end subroutine Putin_FORMAT18_Config_SimBoxArray
  !****************************************************************************

  !****************************************************************************
   subroutine Restore_Config_SimBoxArray(Simbox, Stamp, fname, ExtRestore)
   !***  PURPOSE:   to restore configuration of a simulation box created by last simulation
   !     INPUT:     fname,  the filename storing the configuration
   !
   !     OUTPUT:    Simbox, the simulation box
   !                Stamp,  the record stamp
   !
      implicit none

      !--- DUMMY variables
      character*(*),              intent(in):: fname
      type(SimMDBox),dimension(:)           :: SimBox
      type(MDRecordStamp)                   :: Stamp
      interface 
        subroutine ExtRestore(hFile)
        implicit none
         !--- dummy vaiables
         integer, intent(in)::hFile
        end subroutine ExtRestore
      end interface 

      !--- Local variables
      integer::hFile, NB, I

         call AvailableIOUnit(hFile)
         if(Stamp%NBox .le. 0) then
            NB = size(Simbox)
         else
            NB = Stamp%NBox
         end if

         open(hFile, file = fname, form='unformatted', status='old')
             call Restore_RecordStamp(hFile, Stamp)
             do I=1, NB
                call Restore_Config_SimMDBox(hFile, Simbox(I))
                call Restore_DataPad(hFile,Simbox(I)%ptrDatPad)
             end do
             call ExtRestore(hFile)
         close(hFile)
        return
  end subroutine Restore_Config_SimBoxArray
  !*********************************************************************

  !***************************************************************************
   subroutine Archive_Config_SimBoxArray(Simbox, Stamp, fname, ExtArchive)
   !***  PURPOSE:   to save configuration of a simulation box created by last simulation
   !     INPUT:     Simbox, the simulation box
   !                Stamp,  the record stamp
   !                fname,  the filename storing the configuration
      implicit none

      !--- DUMMY variables
      type(SimMDBox),dimension(:), intent(in):: SimBox
      type(MDRecordStamp),         intent(in):: Stamp
      character*(*),               intent(in):: fname
      interface 
        subroutine ExtArchive(hFile)
        implicit none
         !--- dummy vaiables
         integer, intent(in)::hFile
        end subroutine ExtArchive
      end interface 

      !--- Local variables
      integer::hFile, NB, I

         call AvailableIOUnit(hFile)
         if(Stamp%NBox .le. 0) then
            NB = size(Simbox)
         else
            NB = Stamp%NBox
         end if

         open(hFile, file = fname, form='unformatted', status='unknown')
            call Archive_RecordStamp(hFile, Stamp)
            do I=1, NB
               call Archive_Config_SimMDBox(hFile, Simbox(I))
               call Archive_DataPad(hFile,Simbox(I)%ptrDatPad)
            end do
            call ExtArchive(hFile)
         close(hFile)

        return
  end subroutine Archive_Config_SimBoxArray
  !*********************************************************************

  !****************************************************************************
  subroutine Thermalize_SimBoxArray(SimBox, TL, METHOD, IK0)
  !***  PURPOSE:   to thermalize a set of clusters or the whole simulation box by
  !                by selected method
  !
  !     INPUT:     SimBox, the simulation box to be thermalized
  !                TL,     the temperature the box to be thermalized to
  !                METHOD, the method to thermalize the box
  !                IK0,    the logical flags to determine if a group of atoms will be thermalized,
  !                        this is optional, if IK0 is not present, the whole box will be thermalized
  !
  !     OUTPUT     SimBox, with velocities of atoms assigned
  IMPLICIT NONE
      !--- DUMMY variables
      type(SimMDBox),dimension(:)::SimBox
      real(KINDDF), intent(in)::TL
      integer, intent(in)::METHOD
      logical, dimension(:),intent(in), optional::IK0

      !---Loacal variables
      integer::I

      !******
         if(present(IK0)) then
            do I=1, SIZE(SimBox)
               call Thermalize_SimMDBox(SimBox(I), TL, METHOD, IK0)
            end do
         else
            do I=1, SIZE(SimBox)
               call Thermalize_SimMDBox(SimBox(I), TL, METHOD)
            end do
         end if

      RETURN
  end subroutine Thermalize_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  !**** to get the current thermal quantities
  subroutine Cal_thermal_quantities_SimBoxArray(SimBox)
  !***  PORPOSE: to get the thermal quantities
  !     INPUT:   SimBox
  !     OUTPUT:  SimBOX, with thermal quantities updated
  IMPLICIT NONE
      !--- DUMMY VARIABLES
      type(SimMDBox),dimension(:)::SimBox

      !--- Local variables
      integer::I

        do I=1, size(SimBox)
           call Cal_thermal_quantities_SimMDBox(SimBox(I))
        end do

       return
  end subroutine Cal_thermal_quantities_SimBoxArray
  !****************************************************************************

  !****************************************************************************
  subroutine Putout_ThermalQ_SimBoxArray(ITIME, TIME, ISECT, FNAME, SimBox,con)
  !***  PORPOSE: to putout the instant thermal quantities of a simulation box
  !     INPUT:   ITIME, the ith time step
  !              TIME,  the time
  !              FNAME, the file name storing the configuration
  !              SimBox, the simulation box
  !              CON,    optional, indicating if print out on concsole
  !
  implicit none
      !--- DUMMY variables
      integer, intent(in)::ITIME, ISECT
      real(KINDDF),INTENT(IN)::TIME
      character*(*)::FNAME
      type(SimMDBox),dimension(:), intent(in)::SimBox
      integer,optional::con
      !--- local variables
      integer, save::hFile =0 , I
      logical::opened = .false.
      character*16::str(10)=""
      character*32::strnum
      character*64::fmt
      real(KINDDF)::values(10)

       !$$--- to determine the output format for the head line
         write(strnum,*) SIZE(SimBox)
         strnum = adjustl(strnum)
         fmt ="(10x"//strnum(1:len_trim(strnum))//"(9(A16)))"
       !$$--- at the first time the routine is called
       !$$    we should open the output file
         if(ITIME .eq. 0) then
            if(hFile .gt. 0) close(hFile)
            !$$--- to open the file
            call AvailableIOUnit(hFile)
            open(hFile, file = fname, status='unknown')
            str(1) = "  TIMESECTION   "
            str(2) = "    TIME(PS)    "
            str(3) = "    TEMP.(K)    "
            str(4) = "  VOULUME(LU)   "
            str(5) = "   PRESS0(kb)   "
            str(6) = "   PRESS1(kb)   "
            str(7) = "   PRESST(kb)   "
            str(8) = "    C.E.(ev)    "
            str(9) = "  HARMILT.(cgs) "
            write(hFile, FMT=fmt) str(1:2),(str(3:9), I=1,SIZE(SimBox))
         else if(ITIME .lt. 0) then
              if(hFile .gt. 0) close(hFile)
              hFile = 0
         else
          !$$ for cases restart the simulation
              if(hFile .le. 0) then
                 call AvailableIOUnit(hFile)
                 open(hFile, file = fname, status='unknown')
                !$$--- read the first line
                 read(hFile, *, end=100) str(1)
                 goto 200
                !$$--- the file is empty
  100            str(1) = "  TIMESECTION   "
                 str(2) = "    TIME(PS)    "
                 str(3) = "    TEMP.(K)    "
                 str(4) = "  VOULUME(LU)   "
                 str(5) = "   PRESS0(kb)   "
                 str(6) = "   PRESS1(kb)   "
                 str(7) = "   PRESST(kb)   "
                 str(8) = "    C.E.(ev)    "
                 str(9) = "  HARMILT.(cgs) "
                 write(hFile, FMT=fmt) str(1:2),(str(3:9), I=1,SIZE(SimBox))
                 goto 300

                !$$--- goto the end of the file
  200           continue
                do while(.true.)
                   read(hFile, *, end=300) str(1)
                end do
  300           write(hFile,*)
              end if
         end if

         if(hFile .le. 0) return;

         !$$--- redetermine the format
         fmt ="(1x,I8, 4x, I4, 8x, 1pE14.5,2x, "//strnum(1:len_trim(strnum))//"(7(1pE14.5,2x)))"
         write(hFile, FMT =fmt) ITIME, ISECT,TIME, &
                                (SimBox(I)%TEMPERATURE, SimBox(I)%VOLUME/SimBox(I)%RR**3, &
                                 SimBox(I)%SPRESS0, SimBox(I)%SPRESS1, SimBox(I)%SPRESS, &
                                 SimBox(I)%AVEPOT, SimBox(I)%HARMIL, I=1,SIZE(SimBox))

         if(present(con)) then
            write(*, FMT=fmt) ITIME, ISECT, TIME, &
                              SimBox(1)%TEMPERATURE, &
                              SimBox(1)%VOLUME/SimBox(1)%RR**3, SimBox(1)%SPRESS, &
                              SimBox(1)%AVEPOT, SimBox(1)%HARMIL
         end if
      return
  end subroutine Putout_ThermalQ_SimBoxArray
  !****************************************************************************


  end module MD_SimboxArray
