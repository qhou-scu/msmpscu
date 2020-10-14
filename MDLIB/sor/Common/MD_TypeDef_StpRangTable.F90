  module MD_TYPEDEF_STOPRANGE_TABLE
  !**** DESCRIPTION: this module is to define the data type for stopping power and CSDA
  !                  range.
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2018
  !
  ! **** HOSTORY:
  !                  1. Feb 28, 2018 (HOU Qing):       the first version
  !                  2. Apr 24, 2018 (HOU Qing, Cui JC):
  !                     Add the routine, for loading stopping data from external file.
  !
  !
  !
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    implicit none

     !--- the interface to external routine calculating stoping power
     abstract interface
       subroutine STOPPWR(E, M1,M2, Z1, Z2, ST)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::E, M1,M2, Z1, Z2
          real(KINDDF),intent(out)::ST
        end subroutine STOPPWR
     end interface

     !--- the interface to external routine calculating F(RHO), dF(RHO)/dRHO
     abstract interface
       subroutine NOTIFY_STOPPWR(IDS)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          integer,dimension(:), intent(in)::IDS
        end subroutine NOTIFY_STOPPWR
     end interface

     !--- definition of structures for stopping calculation routines
     integer, parameter, private::mp_IDSSIZE = 2
     type, private::MDStopPowerExtProc
          character*128::STPDescript = ""                          ! the description of availble stopping power
          integer::IDS(mp_IDSSIZE)            = 0
          procedure(STOPPWR),        pointer, nopass::pNSTOP=>null()
          procedure(STOPPWR),        pointer, nopass::pESTOP=>null()
          procedure(NOTIFY_STOPPWR), pointer, nopass::pNOTIFY=>null()
     end type MDStopPowerExtProc

     type, private::MDStopPowerProcList
           type(MDStopPowerExtProc),  pointer::this=>null()
           type(MDStopPowerProcList), pointer::next=>null()
     end type MDStopPowerProcList

     !--- definition of stop-range table structure
     type::MDSTPTable
          integer     ::NTAB=10000                                ! the size of the table
          real(KINDDF)::EMIN=10.D0                                ! the energy range needed (in keV)
          real(KINDDF)::EMAX=1000.D0
          real(KINDDF)::DLTE=(1000.D0-10.D0)/10000

          integer,      dimension(:,:), allocatable::KPAIR        ! the index of table for atom-target
          real(KINDDF), dimension(:),   allocatable::E            ! the energy points
          real(KINDDF), dimension(:,:), allocatable::STPWR        ! the stoping power table

          type(MDStopPowerProcList)::REGPROCS                     ! the list of registered external proceces to generating stopping tables
     end type MDSTPTable

  !--- interface list to the external routine --------------
  !---------------------------------------------------------
     private:: Add_StopPowerProcList
  !---------------------------------------------------------
     private:: CreateKPair_STPTable
  !---------------------------------------------------------
     public::  Create_DefaultExtProcID
  !---------------------------------------------------------
     public::  Create_STPTable
  !---------------------------------------------------------
     public::  Export_STPTable
  !---------------------------------------------------------
     private:: GetExtProc_StopPowerProcList
  !---------------------------------------------------------
     private:: GetIthExtProc_StopPowerProcList
  !---------------------------------------------------------
     private:: IsExProcID
  !---------------------------------------------------------
     public::  Import_STPTable
  !---------------------------------------------------------
     private:: Load_STPTable
  !---------------------------------------------------------
     public::  NULL_STOPPING
  !---------------------------------------------------------
     private:: Release_StopPowerProcList
  !---------------------------------------------------------
     public::  Release_STPTable
  !---------------------------------------------------------
  contains

   !*********************************************************************************
   subroutine NULL_STOPPING(E, M1, M2, Z1, Z2, ST)
   implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::E, M1,M2, Z1, Z2
          real(KINDDF),intent(out)::ST
              ST = 0.D0
              return
   end subroutine NULL_STOPPING
   !*********************************************************************************

   !*********************************************************************************
   subroutine Create_DefaultExtProcID(M1, M2, Z1, Z2, IDS)
   implicit none
       !--- dummy vaiables
       real(KINDDF),intent(in)::M1,M2, Z1, Z2
       integer, dimension(:)  ::IDS

            IDS    = 0
            IDS(1) = int(Z1 + 0.00001)*1000 + int(M1+0.001)
            IDS(2) = int(Z2 + 0.00001)*1000 + int(M2+0.001)
        return
   end subroutine Create_DefaultExtProcID
  !*********************************************************************************

  !****************************************************************************************
  logical function IsExProcID(pExtProc, IDS) result(yes)
  !***  DESCRIPTION:   to get the external process has the IDS
  !     INPUT:   pExtProc,  an external process class
  !              IDS,       type IDS of StopPowerExtProc class to be found
  !
  !     OUTPUT:  pExtProc,  the pointer to a MDStopPowerExtProc class with ID = IDS
  !
  implicit none
     integer, dimension(:)            , intent(in)::IDS
     type(MDStopPowerExtProc), pointer, intent(in)::pExtProc
     !--- Local
     integer::n1, n2

           yes = .false.
           if(.not.associated(pExtProc))  then
               return
           end if

           n1   = count(pExtProc%IDS .gt. 0)
           n2   = count(IDS .gt. 0)
           if(n1 .eq. n2) then
              if(all(pExtProc%IDS(1:n1) .eq. IDS(1:n1))  )then
                 yes = .true.
              end if
           end if
           return
   end function IsExProcID
   !****************************************************************************************

   !*********************************************************************************
   recursive subroutine Release_StopPowerProcList(List)
   !***  PURPOSE:   to to relase the memory allocated in the StopPowerProcList
   !     INPUT:     List: the MDStopPowerProcList
   !     OUTPUT     List: the MDStopPowerProcList
   implicit none
       type(MDStopPowerProcList)::List

            if(associated(List%next)) then
               call Release_StopPowerProcList(List%next)
               deallocate(List%next)
           end if
           if(associated(List%this) ) deallocate(List%this)
           return
  end subroutine Release_StopPowerProcList
  !*********************************************************************************

  !*********************************************************************************
  recursive subroutine Add_StopPowerProcList(LIST, IDS, pESTOP, pNSTOP, pNOTIFY, NOTE)
  !***  DESCRIPTION: to add external sopting power procedures to MDStopPowerProcList
  !     INPUT:   LIST,      the MDStopPowerProcList
  !              ID,        type double IDs of the stopping process
  !              pESTOP,    the external procedure creating electronic ST(E)
  !              pNSTOP,    the external procedure creating nuclear ST(E)
  !              NOTE,      the description about the external routines
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDStopPowerProcList)::List
     integer, dimension(:),       intent(in)::IDS
     procedure(STOPPWR), pointer, intent(in)::pESTOP
     procedure(STOPPWR), pointer, intent(in)::pNSTOP
     procedure(NOTIFY_STOPPWR), pointer, optional, intent(in)::pNOTIFY
     character*(*), optional, intent(in)::NOTE
     !--- Local
       integer::I

           if(.not.associated(LIST%this)) then
               allocate(LIST%this)
               LIST%this%pNSTOP  =>pNSTOP
               LIST%this%pESTOP  =>pESTOP
               LIST%this%IDS = 0
               do I=1, min(size(IDS), size(LIST%this%IDS))
                  LIST%this%IDS(I) = IDS(I)
               end do

               if(present(pNOTIFY) ) then
                  LIST%this%pNOTIFY => pNOTIFY
               else
                  LIST%this%pNOTIFY => null()
               end if

               if(present(NOTE)) then
                  LIST%this%STPDescript = NOTE(1:len_trim(NOTE))
               else
                  LIST%this%STPDescript = ""
               end if

               LIST%next => null()
           else
               if(.not.associated(LIST%next)) then
                 allocate(LIST%next)
               end if
               if(present(NOTE)) then
                 if(present(pNOTIFY) )then
                    call Add_StopPowerProcList(LIST%next, IDS, pESTOP, pNSTOP, pNOTIFY=pNOTIFY, NOTE=NOTE)
                 else
                    call Add_StopPowerProcList(LIST%next, IDS, pESTOP, pNSTOP, NOTE=NOTE)
                 end if
               else

                 if(present(pNOTIFY) )then
                    call Add_StopPowerProcList(LIST%next, IDS, pESTOP, pNSTOP, pNOTIFY=pNOTIFY)
                 else
                    call Add_StopPowerProcList(LIST%next, IDS, pESTOP, pNSTOP)
                 end if
               end if
           end if
           return
   end subroutine Add_StopPowerProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetExtProc_StopPowerProcList(LIST, IDS, pExtProc)
  !***  DESCRIPTION: to get the StopPowerExtProc class of ID ISTOP
  !     INPUT:   LIST,      the StopPowerExtProcList
  !              IDS,       type IDS of StopPowerExtProc class to be found
  !
  !     OUTPUT:  pExtProc,  the pointer to a MDStopPowerExtProc class with ID = IDS
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDStopPowerProcList)::LIST
     integer, dimension(:), intent(in)::IDS
     type(MDStopPowerExtProc), pointer, intent(out)::pExtProc
     !--- Local

           pExtProc => null()
           if(.not.associated(LIST%this)) return

           if(IsExProcID(LIST%this, IDS)) then
              pExtProc => LIST%this
              return
           else
              if(associated(LIST%next)) then
                 call GetExtProc_StopPowerProcList(LIST%next, IDS, pExtProc)
              end if
           end if
           return
   end subroutine GetExtProc_StopPowerProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetIthExtProc_StopPowerProcList(LIST, ITH, pExtProc)
  !***  DESCRIPTION: to get the ith StopPowerExtProc class in LIST
  !                  Note the different betweent this GetExtProc_StopPowerProcList
  !     INPUT:   LIST,      the MDStopPowerExtProc
  !              ITH,       ID of MDStopPowerExtProc class to be get
  !
  !     OUTPUT:  pExtProc,  the pointer to a MDForceExtProc class
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDStopPowerProcList)::List
     integer, intent(in)::ITH
     type(MDStopPowerExtProc), pointer, intent(out)::pExtProc
     !--- Local
     integer::N

           pExtProc => null()
           N = ITH-1
           if(N .eq. 0) then
              pExtProc=>  LIST%this
              return
           else
              if(associated(LIST%next)) then
                 call GetIthExtProc_StopPowerProcList(LIST%next, N, pExtProc)
              end if
           end if
           return
   end subroutine GetIthExtProc_StopPowerProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetNumExtProc_StopPowerProcList(LIST, NUM)
  !***  DESCRIPTION: to get the number ForceExtProc class in LIST
  !     INPUT:   LIST,      the MDStopPowerExtProc
  !
  !     OUTPUT:  NUM,      the number of  MDStopPowerExtProc class in LIST
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDStopPowerProcList)::List
     integer, intent(out)::NUM
     !--- Local
     integer::N

           NUM = 0
           N   = 0
           if(associated(LIST%next)) then
              call GetNumExtProc_StopPowerProcList(LIST%next, N)
           end if

           if(associated(LIST%this)) then
              NUM = N + 1
           end if
           return
   end subroutine GetNumExtProc_StopPowerProcList
  !*********************************************************************************

  !*********************************************************************************
   subroutine Release_STPTable(STPTable, keepreglist)
   !***  PURPOSE:   to to relase the memory allocated in the stop-range table
   !
   !     INPUT:     keepreglist, optional, =1, without delete the registered force
   !     OUTPUT     MDSTPTable:  the table with the stoping, range arraies deallocated
   implicit none
       type(MDSTPTable)::STPTable
       integer, optional::keepreglist
       integer::ERR

       if(allocated(STPTable%KPAIR)) then
          deallocate(STPTable%KPAIR)
       end if

       if(allocated(STPTable%E)) then
          deallocate(STPTable%E)
       end if

       if(allocated(STPTable%STPWR)) then
          deallocate(STPTable%STPWR)
       end if

       if(.not. present(keepreglist)) then
          call Release_StopPowerProcList(STPTable%REGPROCS)
       else
          if(keepreglist .eq.0 ) call Release_StopPowerProcList(STPTable%REGPROCS)
       end if

       return
   end subroutine Release_STPTable
  !***************************************************************************************

  !*********************************************************************************
  subroutine Register_STPTableProc(M1, M2, Z1, Z2, STPTable, ESTOP, NSTOP, NOTIFY, NOTE)
  !***  PURPOSE:   to register the interaction processes between two kinds of atoms
  !     INPUT:     M1,M2, Z1, Z2:  atomic maxx and numbers
  !                ESTOP,    the elctronic stopping function
  !                NSTOP,    optional, the nuclear stopping function
  !                NOTIFY,   optional, the notification function
  !                NOTE,     a description for the stopping
  !     OUTPUT:    STPTable, the stopping table to be registered
  !
  !***************************************************************************************
  implicit none
      real(KINDDF),intent(in)            ::M1,M2, Z1, Z2
      type(MDSTPTable),     intent(inout)::STPTable
  !--- interface to the external routine -------------------
     interface
       subroutine NSTOP(E, M1,M2, Z1, Z2, ST)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::E, M1,M2, Z1, Z2
          real(KINDDF),intent(out)::ST
        end subroutine NSTOP
     end interface

     interface
       subroutine ESTOP(E, M1,M2, Z1, Z2, ST)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::E, M1,M2, Z1, Z2
          real(KINDDF),intent(out)::ST
        end subroutine ESTOP
     end interface

     interface
       subroutine NOTIFY(IDS)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          integer,dimension(:), intent(in)::IDS
        end subroutine NOTIFY
     end interface

      external::ESTOP, NSTOP, NOTIFY
      optional::       NSTOP, NOTIFY
  !--- end interface to the external routine -------------------
      character*(*),optional             ::NOTE

  !--- local variables
      integer::IDS(mp_IDSSIZE)
      procedure(STOPPWR),        pointer::pNSTOP
      procedure(STOPPWR),        pointer::pESTOP
      procedure(NOTIFY_STOPPWR), pointer::pNOTIFY
      type(MDStopPowerExtProc),  pointer::pExtProc


             !$$--- check if IFORCE alread have
              call Create_DefaultExtProcID(M1, M2, Z1, Z2, IDS)
              call GetExtProc_StopPowerProcList(STPTable%REGPROCS, IDS, pExtProc)
              if(associated(pExtProc)) then
                 return
              end if

              pESTOP => ESTOP

              if(present(NSTOP)) then
                 pNSTOP => NSTOP
              else
                 pNSTOP => null()
              end if

              if(present(NOTIFY)) then
                 pNOTIFY => NOTIFY
              else
                 pNOTIFY => null()
              end if

              if(present(NOTE)) then
                 call Add_StopPowerProcList(STPTable%REGPROCS, IDS, pESTOP, pNSTOP, pNOTIFY, NOTE)
              else
                 call Add_StopPowerProcList(STPTable%REGPROCS, IDS, pESTOP, pNSTOP, pNOTIFY)
              end if

            return
  end subroutine Register_STPTableProc
  !***************************************************************************************

  !***************************************************************************************
  subroutine CreateKPair_STPTable(SimBox, KPAIR, IDS)
  !***  PURPOSE:   to create the KPAIR table for atom-taget combinations
  !     INPUT:     SimBox: the simulation box
  !
  !      OUTPUT:   KPAIR,  the KPAIR table
  !                IDS,    the digitial ID of atom-target combinations
  !
  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(SimMDBox),       intent(in) ::SimBox
      integer,dimension(:,:)           ::KPAIR
      integer,dimension(:,:)           ::IDS
      !local variables
      integer::I,J,K, NK, IDS0(mp_IDSSIZE)
      !--
      real(KINDDF)::M1, M2, Z1, Z2

           !$$--- determine how many kind of interactions are concerned
           KPAIR = 0
           IDS    = 0
           NK     = 0
           do I=1, SimBox%NGROUP
              M1 = SimBox%CM(I)/CP_AU2G
              Z1 = SimBox%CZ(I)
              do J=1, SimBox%NGROUP
                  M2 = SimBox%CM(J)/CP_AU2G
                  Z2 = SimBox%CZ(J)
                  call Create_DefaultExtProcID(M1, M2, Z1, Z2, IDS0)
                  do K=1, NK
                     if(all(IDS0(1:mp_IDSSIZE) .eq. IDS(1:mp_IDSSIZE, K)) ) then
                       KPAIR(I,J) = K
                       exit
                    end if
                  end do

                  if(K.GT.NK) then
                     NK           = NK + 1
                     KPAIR(I,J)   = NK
                     IDS(1:mp_IDSSIZE, NK) = IDS0(1:mp_IDSSIZE)
                  end if
               end do
           end do

           return
  end subroutine CreateKPair_STPTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Create_STPTable(SimBox, CtrlParam, STPTable, FNAME)
  !***  PURPOSE:   to register the interaction between two kinds of atoms
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !      OUTPUT:   STBTable, the stopping table registered
  !

  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDSTPTable),  intent(inout)::STPTable
      character*(*),     optional     ::FNAME

      !local variables
      integer::I,J,K,N, NK,IERR,IDS(mp_IDSSIZE, mp_mxGROUP)
      integer::HASTAB(mp_mxGROUP*mp_mxGROUP)

      !--
      real(KINDDF)::STP, DE, M1, M2, Z1, Z2
      type(MDStopPowerExtProc),  pointer::pExtProc

           call Release_STPTable(STPTable, keepreglist=1)
           !$$--- NOTE: the unit conversion
           STPTable%NTAB = CtrlParam%ST_CTRL%NTAB
           STPTable%EMIN = CtrlParam%ST_CTRL%LE*CP_ERGKEV
           STPTable%EMAX = CtrlParam%ST_CTRL%HE*CP_ERGKEV

           !$$--- determine how many kind of interactions are concerned
           allocate(STPTable%KPAIR(SimBox%NGROUP, SimBox%NGROUP))

           !$$--- determine how many kind of interactions are concerned
           call CreateKPair_STPTable(SimBox, STPTable%KPAIR, IDS)
           NK = maxval(STPTable%KPAIR)

           !--- allocate memory for the tables
           allocate(STPTable%E(0:STPTable%NTAB), STPTable%STPWR(0:STPTable%NTAB, NK),STAT=IERR )

           !$$--- create enengy grids
           !$$    NOTE: the unit of E is keV
            DE = (STPTable%EMAX - STPTable%EMIN)/dble(STPTable%NTAB)
            do I = 0, STPTable%NTAB
               STPTable%E(I) = STPTable%EMIN + DE*dble(I)
            end do

           HASTAB = 0
           STPTable%STPWR = 0.D0
           do I=1, SimBox%NGROUP
              !$$---  NOTE: stopping prcesses use atomic mass
              !$$           we should conver the atom mass,
              !$$           which is in g in the simulation box.
              M1 = SimBox%CM(I)/CP_AU2G
              Z1 = SimBox%CZ(I)

               do J=1, SimBox%NGROUP
                  K = STPTable%KPAIR(I,J)
                  if(HASTAB(K) .gt. 0) cycle
                  HASTAB(K) = 1
                  M2 = SimBox%CM(J)/CP_AU2G
                  Z2 = SimBox%CZ(J)

                  call GetExtProc_StopPowerProcList(STPTable%REGPROCS, IDS(:,K), pExtProc)
                  if(associated(pExtProc%pNOTIFY)) then
                     call pExtProc%pNOTIFY(IDS(:,K))
                  end if
                  !$$--- create electronic stopping table

                  do N = 0, STPTable%NTAB
                     call pExtProc%pESTOP(STPTable%E(N),M1,M2,Z1,Z2,STP)
                     STPTable%STPWR(N,K) = STPTable%STPWR(N,K) + STP
                  end do

                  !$$--- create nuclear stopping table
                  if(associated(pExtProc%pNSTOP)) then
                     do N = 0, STPTable%NTAB
                        call pExtProc%pNSTOP(STPTable%E(N),M1,M2,Z1,Z2,STP)
                        STPTable%STPWR(N,K) = STPTable%STPWR(N, K) + STP
                     end do
                  end if
                !---
              end do
           end do
           !$$--- in MD simulation, we use CGS, we convert the keV-cm units
           !      to CGS
           STPTable%EMIN  = STPTable%EMIN*CP_KEVERG
           STPTable%EMAX  = STPTable%EMAX*CP_KEVERG
           STPTable%E     = STPTable%E*CP_KEVERG
           STPTable%STPWR = STPTable%STPWR*CP_KEVERG
           STPTable%DLTE  = DE*CP_KEVERG

           !---- output for debug
           if(present(FNAME)) then
              call Export_STPTable(SimBox, CtrlParam, STPTable, FNAME)
           end if
           return
  end subroutine Create_STPTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Export_STPTable(SimBox, CtrlParam, STPTable, FNAME)
  !***  PURPOSE:   to export the generated ST table
  !     INPUT:     SimBox:    the simulation box
  !                CtrlParam: the control parameters
  !                STBTable,  the stopping table registered
  !                FNAME,     the filename to export the table to
  !      OUTPUT:
  !
  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(SimMDBox),   intent(in) ::SimBox
      type(SimMDCtrl),  intent(in) ::CtrlParam
      type(MDSTPTable), intent(in) ::STPTable
      character*(*),    intent(in) ::FNAME

      !local variables
      integer::I,J, K, hFile, COL
      integer::HASTAB(mp_mxGROUP*mp_mxGROUP)
      character*256::FILE
      character*16::STAB(mp_mxGROUP*mp_mxGROUP), KEYWORD
      equivalence(FILE, STAB)
      !--

              FILE = FNAME(1:len_trim(FNAME))//".stp"
              call AvailableIOUnit(hFile)
              open(FILE=FILE, UNIT=hFile)
                !$$--- write a indentify header
                write(*,fmt="(A,A)")   ' MDPSCU Message: stopping table save to:', FILE(1:len_trim(FILE))
                write(hFile, fmt="(A)") "&MDPSCU_STPTAB.stp"
                write(hFile, fmt="(A, I7, A, 100I5)") "! In the energy-stopping relation table"
                write(hFile, fmt="(A, I7, A, 100I5)") "! energy  in the unit of keV"
                write(hFile, fmt="(A, I7, A, 100I5)") "! stoping in the unit of keV*cm^2"
                write(hFile, fmt="(A, I7, A, 100I5)") "! "
                write(hFile, fmt="(A, I7, A, 100I5)") "&NUMTABLE ", size(STPTable%STPWR,dim=2)
                write(hFile, fmt="(A, I7)")           "&NUMPOINT ", size(STPTable%STPWR,dim=1)

                HASTAB = 0
                COL    = 1
                do I=1, SimBox%NGROUP
                   do J=1, SimBox%NGROUP
                      K = STPTable%KPAIR(I,J)
                      if(HASTAB(K) .gt. 0) cycle
                      STAB(K)   = SimBox%Symb(I)(1:len_trim(SimBox%Symb(I)))//"->"//&
                                  SimBox%Symb(J)(1:len_trim(SimBox%Symb(J)))
                      HASTAB(K) = 1

                      COL = COL + 1
                      write(hFile, fmt="(A, A I4, 1X, 4(F9.2,1X))") "&"//STAB(K)(1:10), " with COL# ", COL, &
                                                              SimBox%CZ(I), SimBox%CM(I)*CP_G2AU, &
                                                              SimBox%CZ(J), SimBox%CM(J)*CP_G2AU
                   end do
               end do


               write(hFile, fmt="(A16,3x, 64(A16))")  "!--- ENERGY     ", STAB(1:size(STPTable%STPWR,dim=2)+1)
               do J=0, STPTable%NTAB
                  write(hFile, fmt="(1x, 40(1PE16.8))")STPTable%E(J)/CP_KEVERG, &
                               (STPTable%STPWR(J,I)/CP_KEVERG, I=1,size(STPTable%STPWR,dim=2))
               end do

              close(hFile)
           return
  end subroutine Export_STPTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Load_STPTable(SimBox, CtrlParam, STPTable, FNAME)
  !***  PURPOSE:   to export the generated ST table
  !     INPUT:     SimBox:    the simulation box
  !                CtrlParam: the control parameters
  !                STBTable,  the stopping table registered
  !                FNAME,     the filename to load the table from
  !      OUTPUT:
  !
  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in) ::SimBox
      type(SimMDCtrl),    intent(in) ::CtrlParam
      type(MDSTPTable)               ::STPTable
      character*(*),      intent(in) ::FNAME

      !local variables
      integer::I,J, N, NK, NCOL, K, IERR, hFile, LINE, IDS(mp_IDSSIZE, mp_mxGROUP),IDS0(mp_IDSSIZE)
      character*256::STR
      character*32 ::STRNUMB(6), KEYWORD, SYMB
      real(KINDDF) ::M1, M2, Z1, Z2
      !--
      integer,      dimension(:), allocatable::HASTAB, COL
      real(KINDDF), dimension(:), allocatable::VAL

             call  Release_STPTable(STPTable, keepreglist=0)

             !$$--- determine how many kind of interactions are concerned
             allocate(STPTable%KPAIR(SimBox%NGROUP, SimBox%NGROUP))
             call CreateKPair_STPTable(SimBox, STPTable%KPAIR, IDS)
             NK = maxval(STPTable%KPAIR)

              call AvailableIOUnit(hFile)
              open(FILE=FNAME, UNIT=hFile, STATUS='OLD')
                !$$--- read a identify header
                 LINE  = 0
                 call GetInputStrLine(hFile, STR, LINE, "!", *100)
                 STR = adjustl(STR)
                 call GetKeyWord("&", STR, KEYWORD)
                 call UpCase(KEYWORD)
                 if(KEYWORD .ne. "&MDPSCU_STPTAB.STP") then
                    write(*,fmt="(A)") ' MDPSCU Error: unknown format for external stop-table'
                    write(*,fmt="(A)") '               the header keyword should be: &MDPSCU_STPTAB.STP'
                    write(*,fmt="(A)") '               process to be stop'
                    stop
                 end if

                 do while(.true.)
                    call GetInputStrLine(hFile, STR, LINE, "!", *100)
                    STR = adjustl(STR)
                    call GetKeyWord("&", STR, KEYWORD)
                    if(len_trim(KEYWORD) .le. 0) exit

                    call UpCase(KEYWORD)
                    select case(KEYWORD(1:len_trim(KEYWORD)))
                           case( "&NUMTABLE")
                                call Extract_Numb(STR,1,N,STRNUMB(1:1))
                                if(N .le. 0 .or. Istr(STRNUMB(1)).le.0) then
                                   write(*,fmt="(A)")       ' MDPSCU Error: the number of table is missed'
                                   write(*,fmt="(A)")       '        Usage: &NUMTABLE  nt '
                                   write(*,fmt="(A, BZI6)") '        check the stopping table file at line:', LINE
                                   write(*,fmt="(A)")       '        process to be stopped'
                                  stop
                                end if
                                NCOL = Istr(STRNUMB(1))
                                if(NCOL .lt. NK) then
                                   write(*,fmt="(A, BZI6)") ' MDPSCU Error: the number of table required is ', maxval(STPTable%KPAIR)
                                   write(*,fmt="(A, BZI6)") '        but the number of taable in this file is ', NK
                                   write(*,fmt="(A)")       '        process to be stopped'
                                  stop
                                end if

                           case( "&NUMPOINT")
                                call Extract_Numb(STR,1,N,STRNUMB(1:1))
                                if(N .le. 0 .or. Istr(STRNUMB(1)).le.0) then
                                   write(*,fmt="(A)")       ' MDPSCU Error: the number of enrgy points is missed'
                                   write(*,fmt="(A)")       '        Usage: NUMPOINT  ne '
                                   write(*,fmt="(A, BZI6)") '        check the stopping table file at line:', LINE
                                   write(*,fmt="(A)")       '        process to be stopped'
                                  stop
                                end if
                                STPTable%NTAB = Istr(STRNUMB(1))-1
                    end select
                 end do

                 !$$--- to check if we have required table
                 allocate(HASTAB(NK), COL(NK))
                 HASTAB = 0
                 COL    = 0
                 do I=1, SimBox%NGROUP
                    do J=1, SimBox%NGROUP
                       K = STPTable%KPAIR(I,J)
                       if(HASTAB(K) .gt. 0) cycle

                         SYMB   = SimBox%Symb(I)(1:len_trim(SimBox%Symb(I)))//"->"//&
                                  SimBox%Symb(J)(1:len_trim(SimBox%Symb(J)))
                         call UpCase(SYMB)

                         rewind(hFile)
                         LINE = 0
                         do while(.true.)
                            call GetInputStrLine(hFile, STR, LINE, "!", *100)
                            STR = adjustl(STR)
                            call GetKeyWord("&", STR, KEYWORD)
                            if(len_trim(KEYWORD) .le. 0) exit

                            call UpCase(KEYWORD)
                            if(KEYWORD(2:len_trim(KEYWORD)) .eq. SYMB(1:len_trim(SYMB)) ) then
                               call Extract_Numb(STR,5,N,STRNUMB)
                               COL(K) = Istr(STRNUMB(1))
                               Z1     = Drstr(STRNUMB(2))
                               M1     = Drstr(STRNUMB(3))
                               Z2     = Drstr(STRNUMB(4))
                               M2     = Drstr(STRNUMB(5))
                               call Create_DefaultExtProcID(M1, M2, Z1, Z2, IDS0)
                               if(any(IDS0(:) .ne. IDS(:,K)) ) then
                                  write(*,fmt="(A)")       ' MDPSCU Error: the atomic number or mass is not consistent for '//&
                                                           SYMB(1:len_trim(SYMB))//' table'
                                  write(*,fmt="(A, BZI6)") '        check the your boxfile and stopping table file'
                                  write(*,fmt="(A)")       '        process to be stopped'
                                  stop
                               end if
                               HASTAB(K) = 1
                             end if
                         end do
                         if(HASTAB(K) .le. 0) then
                            write(*,fmt="(A)")       ' MDPSCU Error: cannot find stopping cross section for '//&
                                                       SYMB(1:len_trim(SYMB))
                            write(*,fmt="(A, BZI6)") '        check the stopping file: '//FNAME(1:len_trim(FNAME))
                            write(*,fmt="(A)")       '        process to be stopped'
                            stop
                         end if
                   end do
               end do

               !$$---  to begin loading the table
               allocate(VAL(1+NCOL))
               allocate(STPTable%E(0:STPTable%NTAB), STPTable%STPWR(0:STPTable%NTAB, NK),STAT=IERR )
               rewind(hFile)
               LINE = 0
               do while(.true.)
                   call GetInputStrLine(hFile, STR, LINE, "!", *100)
                   STR = adjustl(STR)
                   call GetKeyWord("&", STR, KEYWORD)
                   if(len_trim(KEYWORD) .le. 0) exit
               end do

               backspace(hFile)
               do  I=0, STPTable%NTAB
                    read(hFile, *) VAL(1:NCOL+1)
                    STPTable%E(I)           = VAL(1)
                    STPTable%STPWR(I, 1:NK) = VAL(COL(1:NK))
               end do

          close(hFile)
          if(allocated(HASTAB)) deallocate(HASTAB)
          if(allocated(COL))    deallocate(COL)
          if(allocated(VAL))    deallocate(VAL)
          return

  100     continue
          write(*,fmt="(A)")       ' MDPSCU Error: fail to read in the stopping cross section table '
          return
  end subroutine Load_STPTable
  !***************************************************************************************

  !***************************************************************************************
  !***************************************************************************************
  subroutine Import_STPTable(SimBox, CtrlParam, STPTable, FNAME)
  !***  PURPOSE:   to export the generated ST table
  !     INPUT:     SimBox:    the simulation box
  !                CtrlParam: the control parameters
  !                FNAME,     the filename to export the table to
  !
  !      OUTPUT:   STBTable,  the stopping table registered
  !
  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(SimMDBox),   intent(in) ::SimBox
      type(SimMDCtrl),  intent(in) ::CtrlParam
      type(MDSTPTable)             ::STPTable
      character*(*),    intent(in) ::FNAME

      !local variables
        type(MDSTPTable)::STPTable0
        integer::I,J,K, NK,IERR
        real(KINDDF), dimension(:), allocatable::WW, WA, WB, WC
        real(KINDDF)::DE, TAB(3)
        integer::IOP(2)=5
      !--
            write(*,fmt="(A,A)") ' MDPSCU Message: load stopping cross section from: ', FNAME(1:len_trim(FNAME))
            call Load_STPTable(SimBox, CtrlParam, STPTable0, FNAME)

            !$$--- determine how many kind of interactions are concerned
            call  Release_STPTable(STPTable, keepreglist=1)
            allocate(STPTable%KPAIR(size(STPTable0%KPAIR,dim=1), size(STPTable0%KPAIR,dim=2)))

            !$$--- determine how many kind of interactions are concerned
            STPTable%KPAIR = STPTable0%KPAIR
            NK = maxval(STPTable%KPAIR)
            STPTable%EMIN = CtrlParam%ST_CTRL%LE/CP_KEVERG
            STPTable%EMAX = CtrlParam%ST_CTRL%HE/CP_KEVERG
            STPTable%NTAB = CtrlParam%ST_CTRL%NTAB

            if (STPTable%EMAX .gt. maxval(STPTable0%E) )then
                write(*,*) "MDPSCU Error: The maximum energy is larger than the upper limit of input stopping data "
                write(*,*) "Check input data and &EMAX input in control file"
                write(*,*) "Process to be stopped"
                stop
            end if
            if (STPTable%EMIN .lt. minval(STPTable0%E) )then
                write(*,*) "MDPSCU Error: The minimum energy is less than the low limit of input stopping data "
                write(*,*) "Check input data and &EMIN input in control file"
                write(*,*) "Process to be stopped"
                stop
            end if

           !$$--- allocate memeory for stopping data
           allocate(STPTable%E(0:STPTable%NTAB), STPTable%STPWR(0:STPTable%NTAB, NK),STAT=IERR )

           !$$--- allocate working space for interpolation
           allocate(WW(STPTable0%NTAB+1), WA(STPTable0%NTAB+1), WB(STPTable0%NTAB+1), WC(STPTable0%NTAB+1))

           !$$--- create enengy grids
           !$$    NOTE: the unit of E is keV
            DE = (STPTable%EMAX - STPTable%EMIN)/dble(STPTable%NTAB)
            do I = 0, STPTable%NTAB
               STPTable%E(I) = STPTable%EMIN + DE*dble(I)
            end do

            do J=1,  NK
                !$$--- interpolation for STPTable
                call SPLID1(STPTable0%NTAB+1, STPTable0%E, STPTable0%STPWR(0:STPTable0%NTAB, J), WW, IOP, 1,WA,WB,WC)
                do K=0, STPTable%NTAB
                   call SPLID2(STPTable0%NTAB+1, STPTable0%E,  STPTable0%STPWR(0:STPTable0%NTAB,J), WW, 1, STPTable%E(K),TAB)
                   STPTable%STPWR(K,J) = TAB(1)
                end do
           end do

           STPTable%E     = STPTable%E*CP_KEVERG
           STPTable%STPWR = STPTable%STPWR*CP_KEVERG
           deallocate(WW, WA, WB, WC)
            !---
            call Release_STPTable(STPTable0)
  end subroutine Import_STPTable
  !***************************************************************************************


  !***************************************************************************************
  end module MD_TYPEDEF_STOPRANGE_TABLE
