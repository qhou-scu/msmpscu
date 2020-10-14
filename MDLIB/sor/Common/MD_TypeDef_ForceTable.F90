  module MD_TYPEDEF_ForceTable
  !**** DESCRIPTION: this module is to define the data type for force table.
  !                  Thie version is adopted from an MD_TB_EM_Poten_Force_Table.f90( version 2000)
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010
  !
  ! **** HOSTORY:
  !       1. Dec 4, 2014:
  !          In old version,      POTR  ->  0.5*V(r)
  !                               FPOTR -> -dV(r)/dr/r
  !                               POTB  ->  RHO(r)
  !                               FPOTB -> -(dRHO(r)/dr/r)*0.5
  !
  !          In new version,      POTR   ->  0.5*V(r)*r
  !                               FPOTR  -> -(dV(r)/dr)*r
  !                               POTB   ->  RHO(r)
  !                               FPOTB  -> -(dRHO(r)/dr)
  !      2. Dec 21, 2014 (HOU Qing):
  !         The subroutine of exporting internal potential has been
  !         added, subroutine PUTOUT is reomved.
  !
  !      3. Oct 18, 2016 (HOU Qing):
  !         In MDForceExtProc structure, add an member of notifying function. If notifying function
  !         is registered, an notfiy message will be issued when the tabled is created. This is useful
  !         when create external potentials
  !
  !      4. Apr 18, 2018 (HOU Qing):
  !         Add the symbole ID (POTSymbID) in MDForceExtProc structure. Thus, user can register force table
  !         by calling
  !         either:
  !           Register_ForceTableProc(IFORCE, FTable, NNFORCE, NEFORCE, EMBDF, NOTIFY, NOTE) ,
  !         or:
  !           Register_ForceTableProc(Symb1, Symb2, FTable, NNFORCE, NEFORCE, EMBDF, NOTIFY, NOTE)
  !         where Symb1 is the atomic symbol of source atom, and Symb2 is the atomic symbol of
  !         target atom.
  !         In the boxfile of an application, the interaction IDs between atom groups are not necessarily
  !         setup, only if the atomic symbol of atom groups are assinged
  !
  !
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    implicit none

     !--- the interface to external routine calculating V(r), -dV(r)/dr
     private::NNFORCE
     abstract interface
       subroutine NNFORCE(R,POTR, FPOTR)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTR, FPOTR
        end subroutine NNFORCE
     end interface

     !--- the interface to external routine calculating RHO(r), -dRHO(r)/dr
     private::NEFORCE
     abstract interface
       subroutine NEFORCE(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTB, FPOTB
        end subroutine NEFORCE
     end interface

     !--- the interface to external routine calculating F(RHO), dF(RHO)/dRHO
     private::EMBDF
     abstract interface
       subroutine EMBDF(RHO,FRHO, DFRHO)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::RHO
          real(KINDDF),intent(out)::FRHO, DFRHO
        end subroutine EMBDF
     end interface

     !--- the interface to external routine calculating F(RHO), dF(RHO)/dRHO
     private::NOTIFY_EXTID
     abstract interface
       subroutine NOTIFY_EXTID(IFROCE)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          integer,intent(in)::IFROCE
        end subroutine NOTIFY_EXTID
     end interface

     !--- definition of structures for force calculation routines
     type::PairInteractCtrl
           character(len=8), dimension(2) ::SYMB = ""                 ! atomic symbole of interaction atom pair
           real(KINDDF),     dimension(2) ::CZ=0.D0                   ! atomic number of atom pair
           real(KINDDF),     dimension(2) ::CM=0.D0                   ! atomic mass   of atom pair
           real(KINDDF),     dimension(2) ::SR=0.D0                   ! short-range  region used to connected to core potential such as ZBL potential
           real(KINDDF),     dimension(2) ::LR=1.D64                  ! long-range region used to smoothly connecting the potential to
                                                                      ! cut-off distance where the potential and force tend to be zero
     end type PairInteractCtrl

     type, private::MDForceExtProc
           integer                        ::POTTypeID = 0             ! the type id  of available potentials
           character(len=8)               ::POTSymbID = ""            ! the type id  (symbol) of available potentials
           character(len=128)             ::POTDescript = ""          ! the description of availble potentials
           procedure(NNFORCE),      pointer, nopass::pNNFORCE=>null()
           procedure(NEFORCE),      pointer, nopass::pNEFORCE=>null()
           procedure(EMBDF),        pointer, nopass::pEMBDF=>null()
           procedure(NOTIFY_EXTID), pointer, nopass::pNOTIFY=>null()
     end type MDForceExtProc

     type,private::MDForceExtProcList
           type(MDForceExtProc),     pointer::this=>null()
           type(MDForceExtProcList), pointer::next=>null()
     end type MDForceExtProcList

     !--- definition of force table structure
     type::MDForceTable
          character*12 ::PotType = "FS_TYPE"                      ! the name of  potential type, default is FS_TYPE
          character*128::PotLibname = ""                          ! the name of the potential library
          character*32 ::PotSubLibname = ""                       ! the sub name of the potential library

          integer                                  ::NKIND=0      ! the number of types of pair function concerned
          integer                                  ::NTAB=10000   ! the size of the table

          real(KINDDF)                             ::Rmax  = -1   ! the maxma value of R in the table
          real(KINDDF)                             ::RmaxSQRT     ! the square root of Rmax
          real(KINDDF)                             ::CSI          ! the inverse of RmaxSQRT, e.g. NTAB/RmaxSQRT
          real(KINDDF)                             ::CSIV         ! the inverse of CSI         1/CSI
          integer,      dimension(:),   allocatable::FPAIR        ! the pair of force which has been registered
          integer,      dimension(:),   allocatable::HASREG       ! the indicator to indicate if the force table has been created
          integer,      dimension(:,:), allocatable::KPAIR        ! the index of table for atom1 and atom2

          real(KINDDF), dimension(:),   allocatable::R            ! the R points, used when external table are used
          real(KINDDF), dimension(:,:), allocatable::POTR         ! the potential table for Nucl.-Nucl.  interaction with size of  NKIND*NTAB
          real(KINDDF), dimension(:,:), allocatable::FPOTR        ! the differential of potential table POTR
          real(KINDDF), dimension(:,:), allocatable::POTB         ! the table for e-density distribution
          real(KINDDF), dimension(:,:), allocatable::FPOTB        ! the table for differential of density distribution
                                                                  ! *** NOTE: POTR  = V(r) in unit erg
                                                                  !          FPOTR = dV(r)/dr in unit erg/cm
                                                                  !          POTB
                                                                  !          FPOTB = dRHO(r)/dr in unit of POTB/cm

          integer                                  ::NEMBD=10000
          integer                                  ::NKIND1=0     ! the number of typed of embedment function concerned
          integer,      dimension(:),   allocatable::FPAIR1       ! the embedment function which has been registered
          integer,      dimension(:),   allocatable::HASREG1      ! the indicator to indicate if the embedment function has been created
          integer,      dimension(:),   allocatable::KEMBD        ! the index of table for embedment function of atom1
          real(KINDDF), dimension(:),   allocatable::RHO          ! the RHO points, used when external table are used
          real(KINDDF), dimension(:,:), allocatable::FEMBD        ! the table for embemdment function
          real(KINDDF), dimension(:,:), allocatable::DFEMBD       ! the table for differential of th embedment function
          real(KINDDF)::RHOMX = -1                                ! the estimated max value of electron density
          real(KINDDF)::RHOD  = 1                                 ! the estimated RHO interval for F(RHO) table in EAM

          type(MDForceExtProcList)::REGPROCS                      ! the list of registered external proceces to generating force tables
     end type MDForceTable

     real(KINDDF), private::m_RHOMULT=20

  !--- interface list to the external routine --------------
  !---------------------------------------------------------
     private:: Add_ForceExtProcList
  !---------------------------------------------------------
     public::  Create_Interaction_ForceTable
  !---------------------------------------------------------
     public::  Create_Pairwise_ForceTable
  !---------------------------------------------------------
     public::  Create_EMBDFUNTable
  !---------------------------------------------------------
     public::  Export_ForceTable
  !---------------------------------------------------------
     public::  Generate_Lib_ForceTable
  !---------------------------------------------------------
     private:: GetExtProcByTypeID_ForceExtProcList, &
               GetExtProcBySymbID_ForceExtProcList, &
               GetExtProc_ForceExtProcList

     interface GetExtProc_ForceExtProcList
          module procedure GetExtProcByTypeID_ForceExtProcList
          module procedure GetExtProcBySymbID_ForceExtProcList
     end interface GetExtProc_ForceExtProcList
  !---------------------------------------------------------
     public::  GetNumExtProc_ForceExtProcList
  !---------------------------------------------------------
     private:: GetIthExtProc_ForceExtProcList
  !---------------------------------------------------------
     private:: Import_ForceTable
  !---------------------------------------------------------
     private:: New_ForceTable
  !---------------------------------------------------------
     public::  NULL_NNFORCE
  !---------------------------------------------------------
     public::  NULL_NEFORCE
  !---------------------------------------------------------
     public::  NULL_EMBDF
  !---------------------------------------------------------
     public:: Print_Available_Force_ForceTable
  !---------------------------------------------------------
     public:: Print_Concerned_Force_ForceTable
  !---------------------------------------------------------
     private:: RegisterByTypeID_ForceTableProc, &
               RegisterBySymbID_ForceTableProc
     interface Register_ForceTableProc
         module procedure RegisterByTypeID_ForceTableProc
         module procedure RegisterBySymbID_ForceTableProc
     end interface Register_ForceTableProc
  !---------------------------------------------------------
     private:: Release_ForceExtProcList
  !---------------------------------------------------------
     private:: SymbID2TypeID_ForceTableProc
  !---------------------------------------------------------

  !--- end interface to the external routine ---------------

  contains
   !*********************************************************************************
   !*********************************************************************************
   subroutine NULL_NNFORCE(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTB, FPOTB
              POTB  = 0.D0
              FPOTB = 0.D0
              return
   end subroutine NULL_NNFORCE

   subroutine NULL_NEFORCE(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTB, FPOTB

              POTB  = 0.D0
              FPOTB = 0.D0
              return
   end subroutine NULL_NEFORCE

   subroutine NULL_EMBDF(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTB, FPOTB
              POTB  = 0.D0
              FPOTB = 0.D0
          return
   end subroutine NULL_EMBDF

   !*********************************************************************************
   recursive subroutine Release_ForceExtProcList(LIST)
   !***  PURPOSE:   to to relase the memory allocated in the ForceExtProcList
   !     INPUT:     List: the MDForceExtProcList
   !     OUTPUT     List: the MDForceExtProcList
   implicit none
       type(MDForceExtProcList)::List

            if(associated(LIST%next)) then
               call Release_ForceExtProcList(LIST%next)
               deallocate(LIST%next)
           end if
           if(associated(LIST%this) ) deallocate(LIST%this)
           return
  end subroutine Release_ForceExtProcList
  !*********************************************************************************

  !*********************************************************************************
  recursive subroutine Add_ForceExtProcList(LIST, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY, NOTE)
  !***  DESCRIPTION: to add external force procedures to MDForceExtProcList
  !     INPUT:   LIST,      the MDForceExtProcList
  !              IFORCE,    type index of the force table to be created
  !              pNNFORCE,  the external procedure creating V(r) and -dV(r)/r
  !              pNEFORCE,  the external procedure creating RHO(r) and -dRHO(r)/r
  !              pEMBDF,    the external procedure creating F(RHO) and dF(RHO)/dRHO
  !              NOTE,      the description about the external force routines
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDForceExtProcList)::List
     integer, intent(in)::IFORCE
     procedure(NNFORCE), pointer, intent(in)::pNNFORCE
     procedure(NEFORCE), pointer, intent(in)::pNEFORCE
     procedure(EMBDF),   pointer, intent(in)::pEMBDF
     procedure(NOTIFY_EXTID), pointer, optional, intent(in)::pNOTIFY
     character*(*), optional, intent(in)::NOTE
     !--- Local

           if(.not.associated(LIST%this)) then
               allocate(LIST%this)
               LIST%this%pNNFORCE  =>pNNFORCE
               LIST%this%pNEFORCE  =>pNEFORCE
               LIST%this%pEMBDF    =>pEMBDF
               LIST%this%POTTypeID = IFORCE
               if(present(pNOTIFY) ) then
                  LIST%this%pNOTIFY => pNOTIFY
               else
                  LIST%this%pNOTIFY => null()
               end if

               if(present(NOTE)) then
                  LIST%this%POTDescript = NOTE(1:len_trim(NOTE))
               else
                  LIST%this%POTDescript = ""
               end if

               LIST%next => null()
           else
               if(.not.associated(LIST%next)) then
                 allocate(LIST%next)
               end if
               if(present(NOTE)) then
                 if(present(pNOTIFY) )then
                    call Add_ForceExtProcList(LIST%next, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY=pNOTIFY, NOTE=NOTE)
                 else
                    call Add_ForceExtProcList(LIST%next, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, NOTE=NOTE)
                 end if
               else

                 if(present(pNOTIFY) )then
                    call Add_ForceExtProcList(LIST%next, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY=pNOTIFY)
                 else
                    call Add_ForceExtProcList(LIST%next, IFORCE, pNNFORCE, pNEFORCE, pEMBDF)
                 end if
               end if
           end if
           return
   end subroutine Add_ForceExtProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetExtProcByTypeID_ForceExtProcList(LIST, IFORCE, pExtProc)
  !***  DESCRIPTION: to get the ForceExtProc class of ID IFORCE
  !     INPUT:   LIST,      the MDForceExtProcList
  !              IFORCE,    type index of ForceExtProc class to be found
  !
  !     OUTPUT:  pExtProc,  the pointer to a MDForceExtProc class with PotTypeID = IFORCE
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDForceExtProcList)::List
     integer, intent(in)::IFORCE
     type(MDForceExtProc), pointer, intent(out)::pExtProc
     !--- Local

           pExtProc => null()
           if(.not.associated(LIST%this)) return

           if(LIST%this%POTTypeID .eq. IFORCE) then
              pExtProc => LIST%this
              return
           else
              if(associated(LIST%next)) then
                 call GetExtProcByTypeID_ForceExtProcList(LIST%next, IFORCE, pExtProc)
              end if
           end if
           return
   end subroutine GetExtProcByTypeID_ForceExtProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetExtProcBySymbID_ForceExtProcList(LIST, Symb, pExtProc)
  !***  DESCRIPTION: to get the ForceExtProc class of Symb ID
  !     INPUT:   LIST,      the MDForceExtProcList
  !              Symb,      Symb of ForceExtProc class to be found
  !
  !     OUTPUT:  pExtProc,  the pointer to a MDForceExtProc class with PotTypeID = IFORCE
  !
  use MiniUtilities, only:SameFirstWord
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDForceExtProcList) ::List
     character*(*), intent(in)::Symb
     type(MDForceExtProc), pointer, intent(out)::pExtProc
     !--- Local

           pExtProc => null()
           if(.not.associated(LIST%this)) return

           if(SameFirstWord(LIST%this%POTSymbID, Symb)) then
              pExtProc => LIST%this
              return
           else
              if(associated(LIST%next)) then
                 call GetExtProcBySymbID_ForceExtProcList(LIST%next, Symb, pExtProc)
              end if
           end if
           return
   end subroutine GetExtProcBySymbID_ForceExtProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetIthExtProc_ForceExtProcList(LIST, ITH, pExtProc)
  !***  DESCRIPTION: to get the ith ForceExtProc class in LIST
  !                  Note the different betweent this GetExtProc_ForceExtProcList
  !     INPUT:   LIST,      the MDForceExtProcList
  !              ITH,       ID of ForceExtProc class to be get
  !
  !     OUTPUT:  pExtProc,  the pointer to a MDForceExtProc class
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDForceExtProcList)::List
     integer, intent(in)::ITH
     type(MDForceExtProc), pointer, intent(out)::pExtProc
     !--- Local
     integer::N

           pExtProc => null()
           N = ITH-1
           if(N .eq. 0) then
              pExtProc=>  LIST%this
              return
           else
              if(associated(LIST%next)) then
                 call GetIthExtProc_ForceExtProcList(LIST%next, N, pExtProc)
              end if
           end if
           return
   end subroutine GetIthExtProc_ForceExtProcList
  !****************************************************************************************

  !****************************************************************************************
  recursive subroutine GetNumExtProc_ForceExtProcList(LIST, NUM)
  !***  DESCRIPTION: to get the number ForceExtProc class in LIST
  !     INPUT:   LIST,      the MDForceExtProcList
  !
  !     OUTPUT:  NUM,      the number of  ForceExtProc class in LIST
  !
  implicit none
  !--- END INTERFACE --------------------------------
     type(MDForceExtProcList)::List
     integer, intent(out)::NUM
     !--- Local
     integer::N

           NUM = 0
           N   = 0
           if(associated(LIST%next)) then
              call GetNumExtProc_ForceExtProcList(LIST%next, N)
           end if

           if(associated(LIST%this)) then
              NUM = N + 1
           end if
           return
   end subroutine GetNumExtProc_ForceExtProcList
  !****************************************************************************************

  !*********************************************************************************
   subroutine Release_ForceTable(FTable, keepreglist)
   !***  PURPOSE:   to to relase the memory allocated in the force table
   !
   !     INPUT:     keepreglist, optional, =1, without delete the registered force
   !     OUTPUT     FTable: the table with the force deallocated
   implicit none
       type(MDForceTable)::FTable
       integer, optional::keepreglist
       integer::ERR

       if(allocated(FTable%FPAIR)) then
          deallocate(FTable%FPAIR)
       end if

       if(allocated(FTable%HASREG)) then
          deallocate(FTable%HASREG)
       end if

       if(allocated(FTable%KPAIR)) then
          deallocate(FTable%KPAIR)
       end if

       if(allocated(FTable%R)) then
          deallocate(FTable%R)
       end if

       if(allocated(FTable%POTR)) then
          deallocate(FTable%POTR)
       end if

       if(allocated(FTable%POTB)) then
          deallocate(FTable%POTB)
       end if

       if(allocated(FTable%FPOTR)) then
          deallocate(FTable%FPOTR)
       end if

       if(allocated(FTable%FPOTB)) then
          deallocate(FTable%FPOTB)
       end if


       if(allocated(FTable%FPAIR1)) then
          deallocate(FTable%FPAIR1)
       end if

       if(allocated(FTable%HASREG1)) then
          deallocate(FTable%HASREG1)
       end if

       if(allocated(FTable%KEMBD)) then
          deallocate(FTable%KEMBD)
       end if

       if(allocated(FTable%RHO)) then
          deallocate(FTable%RHO)
       end if

       if(allocated(FTable%FEMBD)) then
          deallocate(FTable%FEMBD)
       end if

       if(allocated(FTable%DFEMBD)) then
          deallocate(FTable%DFEMBD)
       end if


       if(.not. present(keepreglist)) then
          call Release_ForceExtProcList(FTable%REGPROCS)
       else
          if(keepreglist .eq.0 ) call Release_ForceExtProcList(FTable%REGPROCS)
       end if

       return
   end subroutine Release_ForceTable
  !***************************************************************************************

  !***************************************************************************************
   subroutine New_ForceTable(SimBox, CtrlParam, FTable)
   !***  PURPOSE:   to allocate required memory for the force table
   !
   !     INPUT:     SimBox: the simulation box
   !                CtrlParam: the control parameters
   !                FTable, the force tabled to be registered
   !
   !     OUTPUT:    FTable, the table with memory allocated
   implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in)   ::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(inout)::FTable

      !--- local vaiables
      integer::NKIND, N, NR, NE, NKIND1
      integer:: I, J, K, ierr= 0

      !***
            call Release_ForceTable(FTable, keepreglist=1)
      !$$*** to start initialize the
            !$$--- determine how many kind of interactions are concerned
            N     = SimBox%NGROUP
            allocate(FTable%FPAIR(N*N),FTable%FPAIR1(N),FTable%HASREG(N*N), FTable%HASREG1(N))

            FTable%FPAIR  =  0
            FTable%FPAIR1 =  0
            FTable%HASREG =  0
            FTable%HASREG1=  0
            !$$--- list out the type of pairwise potential and e-density we needed
            NKIND  = 0
            do I=1, N
               do J=1,N
                  IERR = SimBox%PTYPE(I,J)
                  do K=1, NKIND
                     if(IERR .EQ. FTable%FPAIR(K)) then
                        exit
                     end if
                  end do
                  if(K.GT.NKIND) then
                     NKIND = NKIND+1
                     FTable%FPAIR(NKIND) = IERR
                  end if
               end do
            end do

            !--- the table size
            NR    = CtrlParam%NUMFTABR
            allocate(FTable%KPAIR(N,N), FTable%KEMBD(N) ,STAT=IERR )
            allocate(FTable%POTR(NKIND,NR), FTable%FPOTR(NKIND,NR),        &
                     FTable%POTB(NKIND,NR), FTable%FPOTB(NKIND,NR),        &
                     STAT=ierr )

            FTable%NTAB   = NR
            FTable%NKIND  = NKIND
            FTable%KPAIR  = 0
            FTable%POTR   = 0.D0
            FTable%FPOTR  = 0.D0
            FTable%POTB   = 0.D0
            FTable%FPOTB  = 0.D0

            FTable%Rmax  = maxval(CtrlParam%RU) !maxval(CtrlParam%NB_RM)
            FTable%RmaxSQRT = DSQRT(FTable%Rmax)
            FTable%CSI   = dble(FTable%NTAB)/FTable%RmaxSQRT
            FTable%CSIV  = 1.D0/FTable%CSI

            !$$--- list out the type of embedment function table if we needed
            if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "EAM_TYPE") then
               NE    = CtrlParam%NUMFTABE
               NKIND1 = 0
               do I=1, N
                  IERR = SimBox%PTYPE(I,I)
                  do K=1, NKIND1
                     if(IERR .EQ. FTable%FPAIR1(K)) then
                        exit
                     end if
                  end do
                  if(K.GT.NKIND1) then
                     NKIND1 = NKIND1+1
                     FTable%FPAIR1(NKIND1) = IERR
                  end if
               end do

               allocate(FTable%FEMBD(NKIND1,NE), FTable%DFEMBD(NKIND1,NE),STAT=ierr )
               FTable%NEMBD  = NE
               FTable%NKIND1 = NKIND1
               FTable%KEMBD  = 0
               FTable%FEMBD  = 0.D0
               FTable%DFEMBD = 0.D0
               m_RHOMULT = CtrlParam%RHOSCAL
           end if
           return
   end subroutine New_ForceTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine RegisterByTypeID_ForceTableProc(IFORCE, FTable, NNFORCE, NEFORCE, EMBDF, NOTIFY, NOTE)
  !***  PURPOSE:   to register the interaction processes between two kinds of atoms
  !     INPUT:     IFORCE: type index of the force table to be created
  !                NNFORCE, the pairwise potential function
  !                NEFORCE, optional, the pairwise electron density
  !                EMBDF,   optional, the embedment function
  !                NOTIFY,  optional, the notification function
  !                NOTE,    a description for the force
  !     OUTPUT:    FTable, the force table to be registered
  !
  !***************************************************************************************
  implicit none
      integer,           intent(in)   ::IFORCE
      type(MDForceTable),intent(inout)::FTable
      character*(*),     optional     ::NOTE
  !--- interface to the external routine -------------------
     interface
       subroutine NNFORCE(R,POTR, FPOTR)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTR, FPOTR
        end subroutine NNFORCE
     end interface

     interface
       subroutine NEFORCE(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTB, FPOTB
        end subroutine NEFORCE
     end interface

     interface
       subroutine EMBDF(RHO,FRHO, DFRHO)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::RHO
          real(KINDDF),intent(out)::FRHO, DFRHO
        end subroutine EMBDF
     end interface

     interface
       subroutine NOTIFY(EXTID)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          integer,intent(in)::EXTID
        end subroutine NOTIFY
     end interface

      external::NNFORCE, NEFORCE, EMBDF, NOTIFY
      optional::NEFORCE, EMBDF, NOTIFY
  !--- end interface to the external routine -------------------
  !local variables
     procedure(NNFORCE),      pointer::pNNFORCE
     procedure(NEFORCE),      pointer::pNEFORCE
     procedure(EMBDF),        pointer::pEMBDF
     procedure(NOTIFY_EXTID), pointer::pNOTIFY
     type(MDForceExtProc), pointer::pExtProc
     character*32::TSTR

             !$$--- check if IFORCE alread have
              call GetExtProc_ForceExtProcList(FTable%REGPROCS, IFORCE, pExtProc)
              if(associated(pExtProc)) then
                 write(TSTR,*) IFORCE
                 TSTR = adjustl(TSTR)
                 write(*,fmt="(A)") "MDPSCU Error: the potential type #"//TSTR(1:len_trim(TSTR))//" has been registed"
                 write(*,fmt="(A)") "              check the code where RegisterBySymbID_ForceTableProc have been called."
                 write(*,fmt="(A)") "              Process to be stopped"
                 stop
              end if

              pNNFORCE => NNFORCE

              if(present(NEFORCE)) then
                 pNEFORCE => NEFORCE
              else
                 pNEFORCE => null()
              end if

              if(present(EMBDF)) then
                 pEMBDF => EMBDF
              else
                 pEMBDF => null()
              end if

              if(present(NOTIFY)) then
                 pNOTIFY => NOTIFY
              else
                 pNOTIFY => null()
              end if

              if(present(NOTE)) then
                 call Add_ForceExtProcList(FTable%REGPROCS, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY, NOTE=NOTE)
              else
                 call Add_ForceExtProcList(FTable%REGPROCS, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY)
              end if

            return
  end subroutine RegisterByTypeID_ForceTableProc
  !***************************************************************************************

  !***************************************************************************************
  subroutine RegisterBySymbID_ForceTableProc(Sor, Tgt, FTable, NNFORCE, NEFORCE, EMBDF, NOTIFY, NOTE)
  !***  PURPOSE:   to register the interaction processes between two kinds of atoms
  !     INPUT:     Sor:     symbol of source atoms
  !                Tgt:     symbol of target atoms
  !                NNFORCE, the pairwise potential function
  !                NEFORCE, optional, the pairwise electron density
  !                EMBDF,   optional, the embedment function
  !                NOTIFY,  optional, the notification function
  !                NOTE,    a description for the force
  !     OUTPUT:    FTable, the force table to be registered
  !
  !***************************************************************************************
  use MINIUTILITIES, only:UpCase
  implicit none
      character*(*),      intent(in)   ::Sor, Tgt
      type(MDForceTable), intent(inout)::FTable
      character*(*),      optional     ::NOTE
  !--- interface to the external routine -------------------
     interface
       subroutine NNFORCE(R,POTR, FPOTR)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTR, FPOTR
        end subroutine NNFORCE
     end interface

     interface
       subroutine NEFORCE(R,POTB, FPOTB)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::R
          real(KINDDF),intent(out)::POTB, FPOTB
        end subroutine NEFORCE
     end interface

     interface
       subroutine EMBDF(RHO,FRHO, DFRHO)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::RHO
          real(KINDDF),intent(out)::FRHO, DFRHO
        end subroutine EMBDF
     end interface

     interface
       subroutine NOTIFY(EXTID)
         use MD_CONSTANTS
         implicit none
         !--- dummy vaiables
          integer,intent(in)::EXTID
        end subroutine NOTIFY
     end interface

      external::NNFORCE, NEFORCE, EMBDF, NOTIFY
      optional::NEFORCE, EMBDF, NOTIFY
  !--- end interface to the external routine -------------------
  !local variables
     procedure(NNFORCE),      pointer::pNNFORCE
     procedure(NEFORCE),      pointer::pNEFORCE
     procedure(EMBDF),        pointer::pEMBDF
     procedure(NOTIFY_EXTID), pointer::pNOTIFY
     type(MDForceExtProc),    pointer::pExtProc
     character(Len=8)::SYMB
     integer::I, NUM, IFORCE

             !$$--- check if IFORCE alread have
              SYMB = Tgt(1:len_trim(Tgt))//"<-"//Sor(1:len_trim(Sor))
              call UpCase(SYMB)
              call GetExtProc_ForceExtProcList(FTable%REGPROCS, SYMB, pExtProc)
              if(associated(pExtProc)) then
                 write(*,fmt="(A)") "MDPSCU Error: the potential type by symb"//SYMB(1:len_trim(SYMB))//" has been registed"
                 write(*,fmt="(A)") "              check the code where Register_ForceTableProc have been called."
                 write(*,fmt="(A)") "              Process to be stopped"
                 stop
              end if

              !$$--- automatically assigne a type ID to the symb ID procedure
              call GetNumExtProc_ForceExtProcList(FTable%REGPROCS, NUM)
              IFORCE = 0
              do I=1, NUM
                 call GetIthExtProc_ForceExtProcList(FTable%REGPROCS, I, pExtProc)
                 IFORCE = max(IFORCE, pExtProc%POTTypeID)
              end do
              IFORCE = IFORCE + 1

              pNNFORCE => NNFORCE

              if(present(NEFORCE)) then
                 pNEFORCE => NEFORCE
              else
                 pNEFORCE => null()
              end if

              if(present(EMBDF)) then
                 pEMBDF => EMBDF
              else
                 pEMBDF => null()
              end if

              if(present(NOTIFY)) then
                 pNOTIFY => NOTIFY
              else
                 pNOTIFY => null()
              end if

              if(present(NOTE)) then
                 call Add_ForceExtProcList(FTable%REGPROCS, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY, NOTE=NOTE)
              else
                 call Add_ForceExtProcList(FTable%REGPROCS, IFORCE, pNNFORCE, pNEFORCE, pEMBDF, pNOTIFY)
              end if
              call GetExtProc_ForceExtProcList(FTable%REGPROCS, IFORCE, pExtProc)
              pExtProc%POTSymbID = SYMB

            return
  end subroutine RegisterBySymbID_ForceTableProc
  !***************************************************************************************

   !***************************************************************************************
  integer function SymbID2TypeID_ForceTableProc(Sor, Tgt, FTable) result(IFORCE)
  !***  PURPOSE:   to get the type ID of atom symb
  !     INPUT:     Sor:     symbol of source atoms
  !                Tgt:     symbol of target atoms
  !                FTable, the force table to be registered
  !     OUTPUT     IFORCE
  !***************************************************************************************
  use MINIUTILITIES, only:UpCase
  implicit none
      character*(*),      intent(in)   ::Sor, Tgt
      type(MDForceTable), intent(inout)::FTable
     !local variables
     type(MDForceExtProc),    pointer::pExtProc
     character(Len=8)::SYMB
     integer::I, NUM

             !$$--- check if IFORCE alread have
              SYMB = Tgt(1:len_trim(Tgt))//"<-"//Sor(1:len_trim(Sor))
              call UpCase(SYMB)

              !$$--- automatically assigne a type ID to the symb ID procedure
              call GetNumExtProc_ForceExtProcList(FTable%REGPROCS, NUM)
              IFORCE = 0
              do I=1, NUM
                 call GetIthExtProc_ForceExtProcList(FTable%REGPROCS, I, pExtProc)
                 if(SYMB(1:len_trim(SYMB)) .eq. pExtProc%POTSymbID(1:len_trim(pExtProc%POTSymbID))) then
                    IFORCE =  pExtProc%POTTypeID
                 end if
              end do
            return
  end function SymbID2TypeID_ForceTableProc
  !***************************************************************************************

  !***************************************************************************************
  subroutine Create_Pairwise_ForceTable(IFORCE, FTable)
  !***  PURPOSE:   to register the interaction between two kinds of atoms
  !     INPUT:     IFORCE: type index of the force table to be created
  !     OUTPUT:    FTable, the force table to be registered
  !
  implicit none
      integer, intent(in)::IFORCE
      type(MDForceTable),intent(inout)::FTable

  !local variables
      integer::I, K,N, NTAB
      real(KINDDF)::R, T, CSIV
      procedure(NNFORCE),      pointer::pNNFORCE
      procedure(NEFORCE),      pointer::pNEFORCE
      procedure(NOTIFY_EXTID), pointer::pNOTIFY
      type(MDForceExtProc), pointer::pExtProc
      character*32::TSTR

              !$$--- find out the table index for IFORCE
              K = 0
              do I=1, SIZE(FTable%FPAIR)
                 if(FTable%FPAIR(I) .EQ. IFORCE) then
                    K = I
                    exit
                 end if
              end do

              if(K .GT.0) then
                 !$$--- need this table
                 !$$    if alread registered, return
                 if(FTable%HASREG(K).GT.0) then
                    return
                 end if
              else
                 !$$--- DO not need this table
                 return
              end if

              !$$--- to findout the table generator
              call GetExtProc_ForceExtProcList(FTable%REGPROCS, IFORCE, pExtProc)
              if(.not.associated(pExtProc) ) then
                    call Print_Available_Force_ForceTable(6, FTable)
                    if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
                    write(*,*)
                    write(TSTR,*) IFORCE
                    TSTR = adjustl(TSTR)
                    write(*,fmt="(A)") "MDPSCU Error: the potential type #"//TSTR(1:len_trim(TSTR))//" is not available"
                    write(*,fmt="(A)") "              check if this potential has been registered in Register_process"
                    write(*,fmt="(A)") "              Process to be stopped"
                    stop
              end if

              pNNFORCE => pExtProc%pNNFORCE
              pNEFORCE => pExtProc%pNEFORCE
              pNOTIFY  => pExtProc%pNOTIFY
              if(associated(pNOTIFY)) then
                 call pNOTIFY(IFORCE)
              end if

              CSIV = FTable%CSIV
              NTAB = FTable%NTAB
              do I=1,NTAB
                 T = dble(i)*CSIV
                 R = T*T
                 call  pNNFORCE(R,FTable%POTR(K,i), FTable%FPOTR(K,i))
                 FTable%POTR(K,i)  = FTable%POTR(K,i)*R
                 FTable%FPOTR(K,i) = FTable%FPOTR(K,i)*R
              end do

              if(associated((pNEFORCE))) then
                 do I=1,NTAB
                    T = dble(i)*CSIV
                    R = T*T
                    call pNEFORCE(R,FTable%POTB(K,i), FTable%FPOTB(K,i))
                    !$$--- estimate a range of RHO and the interval for linear interpolation
                    if(FTable%RHOMX .lt. FTable%POTB(K,i)*m_RHOMULT) FTable%RHOMX = FTable%POTB(K,i)*m_RHOMULT

                    !NOTE: POTB and FPOTB do not multiplied by R
                    !FTable%POTR(K,i)  = FTable%POTR(K,i)*R
                    !FTable%FPOTR(K,i) = FTable%FPOTR(K,i)*R
                 end do
              else
                 FTable%POTB(K,:)  = 0.D0
                 FTable%FPOTB(K,:) = 0.D0
              end if

              FTable%RHOD  = FTable%RHOMX/dble(FTable%NEMBD)

              !$$--- indicate the force has been registerd
              FTable%HASREG(K) = 1

            return
  end subroutine Create_Pairwise_ForceTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Create_EMBDFUNTable(IFORCE, FTable)
  !***  PURPOSE:   to create the embedment function of atoms
  !     INPUT:     IFORCE: type index of the force table to be created
  !     OUTPUT:    FTable, the force table to be registered
  !
  implicit none
      integer, intent(in)::IFORCE
      type(MDForceTable),intent(inout)::FTable

     !local variables
      integer::I, K,N, NTAB
      real(KINDDF)::R
      procedure(EMBDF),        pointer::pEMBDF
      procedure(NOTIFY_EXTID), pointer::pNOTIFY
      type(MDForceExtProc),    pointer::pExtProc
      character*32::TSTR

              !$$--- the table index for IFORCE
              K = 0
              do I=1, SIZE(FTable%FPAIR1)
                 if(FTable%FPAIR1(I) .EQ. IFORCE) then
                    K = I
                    exit
                 end if
              end do

              if(K .GT.0) then
                 !$$--- need this table
                 !$$    if alread registered, return
                 if(FTable%HASREG1(K).GT.0) then
                    return
                 end if
              else
                 !$$--- DO not need this table
                 return
              end if

              !$$--- to findout the table generator
              call GetExtProc_ForceExtProcList(FTable%REGPROCS, IFORCE, pExtProc)
              if(.not.associated(pExtProc) ) then
                    call Print_Available_Force_ForceTable(6, FTable)
                    if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
                    write(*,*)
                    write(TSTR,*) IFORCE
                    TSTR = adjustl(TSTR)
                    write(*,fmt="(A)") "MDPSCU Error: the embedment function of potential type #"//TSTR(1:len_trim(TSTR))//" is not available"
                    write(*,fmt="(A)") "              check if this function has been registered by Register_ForceTableProc"
                    write(*,fmt="(A)") "              Process to be stopped"
                    stop
              end if

             pEMBDF   => pExtProc%pEMBDF
             pNOTIFY  => pExtProc%pNOTIFY
             if(associated(pNOTIFY)) then
                call pNOTIFY(IFORCE)
             end if

             if(associated(pEMBDF)) then
                NTAB = FTable%NEMBD
                do I=1,NTAB
                   R = dble(I-1)*FTable%RHOD
                   call  pEMBDF(R,FTable%FEMBD(K,I), FTable%DFEMBD(K,I))
                end do
             else
               FTable%FEMBD(K,:)  = 0.D0
               FTable%DFEMBD(K,:) = 0.D0
             end if
             FTable%HASREG1(K) = 1

            return
  end subroutine Create_EMBDFUNTable
  !***************************************************************************************

  !***************************************************************************************
   subroutine Print_Available_Force_ForceTable(hFile, FTable)
  !***  PURPOSE:   to print out the available inteaction
  !
  !     INPUT:     hFile, the IO unit for output
  !                FTable, the force table has been registered
  !
  implicit none
     integer::hFile
     type(MDForceTable),intent(inout)::FTable

     !local variables
     character*256::title
     character*32::tstr
     integer::I,J, NUM
     type(MDForceExtProc), pointer::pExtProc

            title = "!************ THE AVAILALE FORCE TABLES *** "
            write(hFile,*) title(1:len_trim(title))
            title = "!     Type of potential: "//FTable%PotType(1:len_trim(FTable%PotType))
            write(hFile,*) title(1:len_trim(title))

            title = "!     Library name of potential: "//FTable%PotLibname(1:len_trim(FTable%PotLibname))
            write(hFile,*) title(1:len_trim(title))

            if(len_trim(FTable%PotSubLibname) .gt. 0) then
               title = "!     subtype in the library: "//FTable%PotSubLibname(1:len_trim(FTable%PotSubLibname))
               write(hFile,*) title(1:len_trim(title))
            end if


            call GetNumExtProc_ForceExtProcList(FTable%REGPROCS, NUM)
            do I=1, NUM
               call GetIthExtProc_ForceExtProcList(FTable%REGPROCS, I, pExtProc)
               write(tstr,*) pExtProc%POTTypeID
               tstr = adjustl(tstr)
               if(len_trim(pExtProc%POTSymbID) .gt. 0) then
                  title = " !     Potential ID #"//tstr(1:len_trim(tstr))//"("// &
                            pExtProc%POTSymbID(1:len_trim(pExtProc%POTSymbID))//"): "
               else
                  title = " !     Potential ID #"//tstr(1:len_trim(tstr))//": "
               end if
               title = title(1:32)//pExtProc%POTDescript(1:len_trim(pExtProc%POTDescript) )
               write(hFile,fmt="(A)")title(1:len_trim(title))
            end do

            return
  end subroutine Print_Available_Force_ForceTable
  !***************************************************************************************

  !***************************************************************************************
   subroutine Print_Concerned_Force_ForceTable(hFile, FTable)
  !***  PURPOSE:   to print out the concerned inteaction in cuurent simualtion
  !
  !     INPUT:    FTable, the force table has been registered
  !
  implicit none
     integer::hFile
     type(MDForceTable),intent(in)::FTable
     !local variables
     character*256::title
     character*32::tstr
     integer::I,J
     type(MDForceExtProc), pointer::pExtProc


            title = "!************ THE CONCERNED FORCE TABLES IN THIS CALCULATION *** "
            write(hFile,*) title(1:len_trim(title))
            DO I=1, size(FTable%FPAIR)
               J = FTable%FPAIR(I)
               if( J.LE. 0) exit
               call  GetExtProc_ForceExtProcList(FTable%REGPROCS, J, pExtProc)

               write(tstr,*) J
               tstr = adjustl(tstr)
               if(len_trim(pExtProc%POTSymbID) .gt. 0) then
                  title = " !     Potential ID #"//tstr(1:len_trim(tstr))//"("// &
                            pExtProc%POTSymbID(1:len_trim(pExtProc%POTSymbID))//"): "
               else
                  title = " !     Potential ID #"//tstr(1:len_trim(tstr))//": "
               end if

               title = title(1:32)//pExtProc%POTDescript(1:len_trim(pExtProc%POTDescript) )
               write(hFile,fmt="(A)")title(1:len_trim(title))
            END DO
            write(hFile,*)""

            return
  end subroutine Print_Concerned_Force_ForceTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Create_Interaction_ForceTable(SimBox, CtrlParam, FTable, FNAME)
  !***  PURPOSE:   to register the interaction between two kinds of atoms
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !      OUTPUT:    FTable, the force table registered
  !

  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(SimMDBox)                    ::SimBox
      type(SimMDCtrl),     intent(in)   ::CtrlParam
      type(MDForceTable),  intent(inout)::FTable
      character*(*),       optional     ::FNAME

      !local variables
      integer::I,J,K, IFORCE
      character(len=2)::SYMB1, SYMB2

           !$$--- first, we check if each atom pair are given interaction ID
           do I=1, SimBox%NGROUP
              do J=1,SimBox%NGROUP
                 if(SimBox%PTYPE(I,J) .eq. 0) then
                    call Extract_SubSymb(SimBox%Symb(J), SYMB1)
                    call Extract_SubSymb(SimBox%Symb(I), SYMB2)
                    call UpCase(SYMB1)
                    call UpCase(SYMB2)
                    SimBox%PTYPE(I,J) = SymbID2TypeID_ForceTableProc(Sor=SYMB1, Tgt=SYMB2, FTable=FTable)
                 end if
                 if(SimBox%PTYPE(I,J) .eq. 0) then
                    write(*,fmt="(A, I4, 1x, A, I4)") "MDPSCU Error: cannot find the registed force table for between atom groups", &
                                                       I, " and ", J
                    write(*,fmt="(A)")                "              Process to be stopped"
                    stop
                 end if
             end do
           end do

           !$$--- then, we create the table
           call New_ForceTable(SimBox, CtrlParam, FTable)
           FTable%KPAIR = -1
           do I=1, SimBox%NGROUP
              do J=1,SimBox%NGROUP
                 IFORCE = SimBox%PTYPE(I,J)
                 call Create_Pairwise_ForceTable(IFORCE, FTable)

                 do K=1, FTable%NKIND
                    if(IFORCE .EQ. FTable%FPAIR(K)) then    ! alread resgister
                       FTable%KPAIR(I, J) = K
                       exit
                    end if
                 end do
                !---
              end do
           end do

           if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "EAM_TYPE") then
              FTable%KEMBD = -1
              do I=1, SimBox%NGROUP
                 IFORCE = SimBox%PTYPE(I,I)
                 call Create_EMBDFUNTable(IFORCE, FTable)

                 do K=1, FTable%NKIND1
                    if(IFORCE .EQ. FTable%FPAIR1(K)) then    ! alread resgister
                       FTable%KEMBD(I) = K
                      exit
                    end if
                 end do

              end do
              !$$--- if FTable%RHOD=0, indicating only pairwire potential
              !$$    hase been registered, we set an finite value to m_RHOD
              if(FTable%RHOD .le. 1.D-64) FTable%RHOD = 1.D0
           end if

           if(present(FNAME)) then
              call Export_ForceTable(FNAME, FTable)
           end if
           return
  end subroutine Create_Interaction_ForceTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Generate_Lib_ForceTable(FTable, FNAME)
  !***  PURPOSE:   to generate a force table that could be used to by other MD packaged.
  !                The difference between this routine and Create_Interaction_ForceTable
  !                is that all the registered force table will be generated, and output.
  !
  !
  !                 NOTE: This routine dose not need SimBox, the range of R and
  !                       the number of grid pointer of the table should be given before
  !                       calling this routine.
  !
  !      INPUT:     FNAME: the name of the library
  !                 Fable: the table with its forces have been registered
  !
  !      OUTPUT:    FTable, the force table store in a file
  !

  !***************************************************************************************
  implicit none
      !--- dummy vaiables
      type(MDForceTable),intent(inout)::FTable
      character*(*), optional::         FNAME

      !local variables
      integer::IFORCE, NUM, NR, NE, IERR

            call Release_ForceTable(FTable, keepreglist=1)

            if(FTable%Rmax .le. 0.D0) then
              write(*,fmt="(A)") "MDPSCU Error: the distance range is missed in generating force table"
              write(*,fmt="(A)") "              Process to be stopped"
              stop
            end if

            if(FTable%NTAB .le. 0.D0) then
              write(*,fmt="(A)") "MDPSCU Error: the number of grid is less than zero in generating force table"
              write(*,fmt="(A)") "              Process to be stopped"
              stop
            end if

            FTable%RmaxSQRT = DSQRT(FTable%Rmax)
            FTable%CSI   = dble(FTable%NTAB)/FTable%RmaxSQRT
            FTable%CSIV  = 1.D0/FTable%CSI

            call GetNumExtProc_ForceExtProcList(FTable%REGPROCS, NUM)
            NR = FTable%NTAB

            allocate(FTable%FPAIR(NUM), FTable%HASREG(NUM), STAT=IERR )
            allocate(FTable%POTR(NUM,NR), FTable%FPOTR(NUM,NR),        &
                     FTable%POTB(NUM,NR), FTable%FPOTB(NUM,NR),        &
                     STAT=IERR )
            FTable%POTR   = 0.D0
            FTable%FPOTR  = 0.D0
            FTable%POTB   = 0.D0
            FTable%FPOTB  = 0.D0
            FTable%FPAIR  = 0
            FTable%HASREG = 0
            do IFORCE=1, NUM
              FTable%FPAIR(IFORCE) = IFORCE
              call Create_Pairwise_ForceTable(IFORCE, FTable)
            end do

            NE = FTable%NEMBD
            allocate(FTable%FEMBD(NUM,NE), FTable%DFEMBD(NUM,NE),FTable%FPAIR1(NUM), FTable%HASREG1(NUM), STAT=IERR)
            FTable%FEMBD  = 0.D0
            FTable%DFEMBD = 0.D0
            FTable%HASREG1 = 0
            FTable%FPAIR1  = 0
            do IFORCE=1, NUM
               FTable%FPAIR1(IFORCE) = IFORCE
               call Create_EMBDFUNTable(IFORCE, FTable)
            end do

            if(present(FNAME) ) then
               call Export_ForceTable(FNAME, FTable)
            end if

           return
  end subroutine Generate_Lib_ForceTable
  !***************************************************************************************

  !***************************************************************************************
   subroutine Export_ForceTable(FNAME, FTable)
  !***  PURPOSE:   to write force table to a file
  !
  !     INPUT:     FNAME, the filename for the forcetable will write to
  !                FTable, the force table has been registered
  !
  implicit none
      character*(*)::FNAME
      type(MDForceTable),intent(inout)::FTable

  !local variables
     character*1024::title
     character*32::tstr
     integer::hFile,I, J,IT, NC
     real(KINDDF)::R, T, RHOUNIT
     real(KINDDF), dimension(:), allocatable::DEN, DENI
     integer, dimension(:), allocatable::TID

            !$$--- putout the pairwise functions
            title = FNAME(1:len_trim(FNAME))//".pair"
            call AvailableIOUnit(hFile)
            open(FILE=title, UNIT=hFile)
            !$$--- write a indentify header
            write(hFile, fmt="(A)") "&MDPSCU_POTTAB.Pair"
            call Print_Available_Force_ForceTable(hFile, FTable)
            call Print_Concerned_Force_ForceTable(hFile, FTable)
            write(*,fmt="(A,A)") ' MDPSCU Message: pair forcetable save to:', title(1:len_trim(title))

            NC = count(FTable%FPAIR.gt. 0)
            if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "FS_TYPE") then
               RHOUNIT = CP_ERGEV*CP_ERGEV
            else
               RHOUNIT = C_UN
            end if

           !$$--- The FPAIR is not neceesarily in an incremental order
           !$$    we reorder the table ID in an incremental order
            allocate(TID(NC))
            do I=1, NC
               J = count(FTable%FPAIR .le. FTable%FPAIR(I) .and. &
                         FTable%FPAIR .gt. 0  )
               TID(J) = I
            end do

            write(hFile, fmt="(A, A)")            "!     "
            write(hFile, fmt="(A, A)")            "!     NOTE: For FS potential, RHO is in eV^2, the potential calculated by"
            write(hFile, fmt="(A, A)")            "!           -sqrt(RHO)+V(r) "
            write(hFile, fmt="(A, A)")            "!           is thus in unit eV, and the d(RHO)/dr in unit eV^2/A."
            write(hFile, fmt="(A, A)")            "!           "
            write(hFile, fmt="(A, A)")            "!           For EAM potential, not specific unit is assigned for RHO,"
            write(hFile, fmt="(A, A)")            "!           and for d(RHO)/dr. The potential calculated by "
            write(hFile, fmt="(A, A)")            "!           F(RHO)+V(r) is in unit eV."
            write(hFile, fmt="(A, A)")            "!           "

            write(hFile, fmt="(A, A)")            "&POTTYPE ", '"'//FTable%POTTYPE(1:len_trim(FTable%POTTYPE))//'"'
            write(hFile, fmt="(A, I7, A, 100I5)") "&NUMTABLE ", NC, " table IDs: ", (FTable%FPAIR(TID(J)), J=1, NC)
            write(hFile, fmt="(A, I7)")           "&NUMPOINT ", FTable%NTAB
            title = ""
            do J=1, NC
               write(tstr,*) FTable%FPAIR(TID(J))
               tstr = adjustl(tstr)
               title = title(1:len_trim(title))//"      r*V"//tstr(1:len_trim(tstr))//"(r)[eV*A] "   &
                                               //"       -r*dV"//tstr(1:len_trim(tstr))//"/dr[eV] "  &
                                               //"           RHO"//tstr(1:len_trim(tstr))//"(r) "    &
                                               //"           -dRHO"//tstr(1:len_trim(tstr))//"/dr[/A] "
            end do
            write(hFile, fmt="(A7,A21,1x,A)") "&#     ","      Rij(A)     ", title(1:len_trim(title))

            write(tstr,*) NC
            title = "(I6, 1PE12.5,1x"//tstr(1:len_trim(tstr))//"(1PE13.5,1x))"
            do IT=1, FTable%NTAB
               T = dble(IT)*FTable%CSIV
               R = (T*T)*CP_CM2A
               !*** NOTE: POTR  = V(r)*r,     in unit erg*cm
               !          FPOTR = dV(r)/dr*r, in unit erg
               !          POTB
               !          FPOTB = dRHO(r)/dr  in unit RHO/cm
                write(hFile, fmt="(1x,I6, 1x, 256(1PE21.9))")IT,R,(FTable%POTR(TID(J),IT)*2.D0*CP_ERGEV*CP_CM2A,       &
                                                                   FTable%FPOTR(TID(J),IT)*CP_ERGEV,                   &
                                                                   FTable%POTB(TID(J),IT)*RHOUNIT,                     &
                                                                   FTable%FPOTB(TID(J),IT)*RHOUNIT*CP_A2CM,            &
                                                                   J=1, NC)

            end do
            close(hFile)

            !$$--- putout the embedment functions
            title = FNAME(1:len_trim(FNAME))//".embd"
            call AvailableIOUnit(hFile)
            open(FILE=title, UNIT=hFile)
            !$$--- write a indentify header
            write(hFile, fmt="(A)") "&MDPSCU_POTTAB.Embd"
            call Print_Available_Force_ForceTable(hFile, FTable)
            call Print_Concerned_Force_ForceTable(hFile, FTable)
            write(*,fmt="(A,A)") ' MDPSCU Message: embedment function save to:', title(1:len_trim(title))

            NC = count(FTable%FPAIR1.gt. 0)
            do I=1, NC
               J = count(FTable%FPAIR1 .le. FTable%FPAIR1(I) .and. &
                         FTable%FPAIR1 .gt. 0 )
               TID(J) = I
            end do
            if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "FS_TYPE" .and. NC .eq. 0) NC = 1

            write(hFile, fmt="(A, A)")            "&POTTYPE ", '"'//FTable%POTTYPE(1:len_trim(FTable%POTTYPE))//'"'
            write(hFile, fmt="(A, I7, A, 100I5)") "&NUMTABLE ", NC, " table IDs: ", (FTable%FPAIR1(TID(J)), J=1, NC)
            write(hFile, fmt="(A, I7)")           "&NUMPOINT ", FTable%NEMBD
            title = ""
            do J=1, NC
               write(tstr,*) FTable%FPAIR1(TID(J))
               tstr = adjustl(tstr)
               title = title(1:len_trim(title))//"         F"//tstr(1:len_trim(tstr))//"(RHO)" &
                                               //"      dF"//tstr(1:len_trim(tstr))//"/d(RHO)"
            end do
            write(hFile, fmt="(A7,A14,1x,A)") "&#     ","      RHO      ", title(1:len_trim(title))

            write(tstr,*) NC
            title = "(I6, 1PE12.5,1x"//tstr(1:len_trim(tstr))//"(1PE13.5,1x))"

            allocate (DEN(NC), DENI(NC))
            do IT=1, FTable%NEMBD
               R = dble(IT-1)*FTable%RHOD*RHOUNIT
               if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "FS_TYPE") then
                  do J=1, NC
                     if(R .gt. 0.D0) then
                        DEN(J)  = -dsqrt(R)/CP_ERGEV
                        DENI(J) = C_HALF/DEN(J)
                     else
                        DEN(J) = 0.D0
                        DENI(J) = 0.D0
                     end if
                  end do
               else
                  do J=1, NC
                     DEN(J)  = FTable%FEMBD(TID(J),IT)
                     DENI(J) = FTable%DFEMBD(TID(J),IT)
                  end do
               end if
               write(hFile, fmt="(1x,I6, 1x, 256(1PE16.8))")IT,R,(DEN(J)*CP_ERGEV,  DENI(J),  J=1, NC)
           end do
           deallocate (DEN, DENI, TID)
           close(hFile)

           return
  end subroutine Export_ForceTable
  !***************************************************************************************

  !***************************************************************************************
   subroutine Import_ForceTable(FNAME, FTable)
  !***  PURPOSE:   to imp[ort force table from a file
  !
  !     INPUT:     FNAME,  the filename for the forcetable will write to
  !     OUPUT:     FTable, the force table has been registered
  !
  use MINIUTILITIES
  implicit none
  !--- dummy variables
     character*(*)::FNAME
     type(MDForceTable),intent(inout)::FTable

  !--- local variables
     character*256::STR
     character*32::STRNUMB(mp_MXGROUP*mp_MXGROUP+1), KEYWORD
     integer::hFile,I, J,IT, N, NC, NT, LINE, IERR
     real(KINDDF)::RHOUNIT

            call Release_ForceTable(FTable)
            !$$--- putout the pairwise functions
            call AvailableIOUnit(hFile)
            open(FILE=FNAME(1:len_trim(FNAME))//".pair", UNIT=hFile, STATUS='old')
            write(*,fmt="(A,A)") ' MDPSCU Message: import pair forcetable from: ', FNAME(1:len_trim(FNAME))//".pair"

            !$$--- extract the number of types of forces
            FTable%POTType = ""
            LINE = 0
            do while(.true.)
               call GetInputStrLine(hFile,STR, LINE, "!", *100)
               STR = adjustl(STR)
               call GetKeyWord("&", STR, KEYWORD)
               call UpCase(KEYWORD)
               select case (KEYWORD(1:len_trim(KEYWORD)) )
                      case ( "&NUMTABLE ")
                           call Extract_Numb(STR,1,N,STRNUMB)
                           NC = ISTR(STRNUMB(1))
                           allocate(FTable%FPAIR(NC))
                           FTable%FPAIR = 0
                           call Extract_Numb(STR,NC+1,N,STRNUMB)
                           do I=1, N-1
                              FTable%FPAIR(I) = ISTR(STRNUMB(I+1))
                           end do

                      case ( "&NUMPOINT ")
                           call Extract_Numb(STR,1,N,STRNUMB)
                           NT = ISTR(STRNUMB(1))

                      case ( "&POTTYPE ")
                           call Extract_Substr(STR,1,N,STRNUMB)
                           FTable%POTTYPE = STRNUMB(1)

                      case ( "&#")
                           exit
               end select
            end do
            allocate(FTable%POTR(NC,NT),   FTable%FPOTR(NC,NT),        &
                     FTable%POTB(NC,NT),   FTable%FPOTB(NC,NT),        &
                     FTable%R(NT), STAT=ierr )

            FTable%NKIND = NC
            FTable%NTAB  = NT
            !*** not the unit
            if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "FS_TYPE") then
               RHOUNIT  = CP_ERGEV*CP_ERGEV
            else
               RHOUNIT  = C_UN
            end if

            do J=1, FTable%NTAB
               LINE = LINE + 1
               read(hFile, *, err=100) IT, FTable%R(J),(FTable%POTR(I,J), FTable%FPOTR(I,J), &
                                                FTable%POTB(I,J), FTable%FPOTB(I,J), I=1, NC)
            end do
            !*** NOTE: the imported table:
            !          R                   in unit A
            !          POTR  = V(r)*r,     in unit eV*A
            !          FPOTR = dV(r)/dr*r, in unit eV
            !          POTB
            !          FPOTB = dRHO(r)/dr  in unit RHO/A
            FTable%R     = FTable%R*CP_A2CM
            FTable%POTR  = FTable%POTR*C_HALF*CP_EVERG*CP_A2CM
            FTable%FPOTR = FTable%FPOTR*CP_EVERG
            FTable%POTB  = FTable%POTB/RHOUNIT
            FTable%FPOTB = FTable%FPOTB*CP_CM2A/RHOUNIT
            close(hFile)

            !$$--- import the embedment functions
            call AvailableIOUnit(hFile)
            open(FILE= FNAME(1:len_trim(FNAME))//".embd", UNIT=hFile, STATUS='OLD')
            write(*,fmt="(A,A)") ' MDPSCU Message: import embedment function from: ',  FNAME(1:len_trim(FNAME))//".embd"

            !$$--- extract the number of types of forces
            LINE = 0
            do while(.true.)
               call GetInputStrLine(hFile,STR, LINE, "!", *200)
               STR = adjustl(STR)
               call GetKeyWord("&", STR, KEYWORD)
               call UpCase(KEYWORD)
               select case (KEYWORD(1:len_trim(KEYWORD)) )
                      case ( "&NUMTABLE ")
                           call Extract_Numb(STR,1,N,STRNUMB)
                           NC = ISTR(STRNUMB(1))
                           allocate(FTable%FPAIR1(NC))
                           FTable%FPAIR1 = 0
                           call Extract_Numb(STR,NC+1,N,STRNUMB)
                           do I=1, N-1
                              FTable%FPAIR1(I) = ISTR(STRNUMB(I+1))
                           end do

                      case ( "&NUMPOINT ")
                           call Extract_Numb(STR,1,N,STRNUMB)
                           NT = ISTR(STRNUMB(1))

                      case ( "&POTTYPE ")
                           call Extract_Substr(STR,1,N,STRNUMB)
                           FTable%POTTYPE = STRNUMB(1)

                      case ( "&#")
                           exit
               end select
            end do
            allocate (FTable%FEMBD(NC, NT), FTable%DFEMBD(NC,NT),FTable%RHO(NT), STAT=ierr )

            FTable%NKIND1= NC
            FTable%NEMBD = NT
            do J=1, FTable%NEMBD
               LINE = LINE + 1
               read(hFile, *, err=200) IT,FTable%RHO(J),(FTable%FEMBD(I, J), FTable%DFEMBD(I, J),  I=1, NC)
            end do
            FTable%RHO = FTable%RHO/RHOUNIT
            FTable%FEMBD = FTable%FEMBD*CP_EVERG
            close(hFile)

           return
           !*************************************************
100        close(hFile)
           write(*,fmt="(A, I8)") " MDPSCU Error: error in loading external force table"// FNAME(1:len_trim(FNAME))&
                                   //".pair, at line ", LINE
           write(*,fmt="(A)") "               Process to be stopped"
           stop
200        close(hFile)
           write(*,fmt="(A, I8)") " MDPSCU Error: error in loading external force table"// FNAME(1:len_trim(FNAME))&
                                  //".embd, at line", LINE
           write(*,fmt="(A)") "               Process to be stopped"
           stop
  end subroutine Import_ForceTable
  !***************************************************************************************

  !***************************************************************************************
  subroutine Register_Imported_ForceTable(SimBox, CtrlParam, FTable, printout)
  !***  PURPOSE:   to import a force table and register the force table to be used
  !                in force calculations
  !
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !                FTable, the force tabled to be registered
  !
  !     OUTPUT:    FTable, the registered table
  !-----------------------------------------------------------------------------------------
  use MINIUTILITIES
  implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(inout)::FTable
      integer,           optional     ::printout

      !--- local vaiables
      type(MDForceTable)::tFTable
      real(KINDDF), dimension(:), allocatable::WW, WA, WB, WC
      real(KINDDF)::CSIV, T, TT, TAB(3), TRMAX
      integer::I, J, K, KTAB0, KTAB, IFORCE, IJ, IOP(2)=5
      character*32::tstr


      !*** register the table according to user supplied routines
            call Import_ForceTable(SimBox%PotLibname, tFTable)
            do I=1, SimBox%NGROUP
              do J=1,SimBox%NGROUP
                 if(.not. any(tFTable%FPAIR .eq. SimBox%PTYPE(I,J)) ) then
                     write(*,fmt="(A,I2,A,I2,A,I2,A)") ' MDPSCU Error: cannot find force table # ',SimBox%PTYPE(I,J), ' in lib '// &
                                                        SimBox%PotLibname(1:len_trim(SimBox%PotLibname))//'.pair.'
                     write(*, fmt="(A)")               '                check box file or the force table file'
                     write(*, fmt="(A)")               '                Process to be stopped'
                     stop
                 end if
              end do
           end do

           do I=1, SimBox%NGROUP
               if(.not. any(tFTable%FPAIR1 .eq. SimBox%PTYPE(I,I)) ) then
                  write(*,fmt="(A,I2,A,I2,A,I2,A)") ' MDPSCU Error: cannot find force table # ',SimBox%PTYPE(I,I), ' in lib '// &
                                                  SimBox%PotLibname(1:len_trim(SimBox%PotLibname))//'.embd'
                  write(*, fmt="(A)")               '                check box file or the force table file'
                  write(*, fmt="(A)")               '                Process to be stopped'
                  stop
              end if
           end do


           write(*,fmt="(A,I2,A,I2,A,I2,A)") ' MDPSCU Message: to register force table...'
           do K=1, size(tFTable%FPAIR)
              if(tFTable%FPAIR(K) .gt. 0) then
                 write(tstr,*) tFTable%FPAIR(K)
                 tstr = adjustl(tstr)
                 call Register_ForceTableProc(tFTable%FPAIR(K), FTable, NNFORCE = NULL_NNFORCE, NOTE="Table # "//tstr(1:len_trim(tstr)))
              end if
           end do
      !****  create the force table
           FTable%PotType = tFTable%PotType
           call New_ForceTable(SimBox, CtrlParam, FTable)
           FTable%KPAIR = -1
           FTable%KEMBD = -1
           CSIV = FTable%CSIV

           !$$--- to perform interpolation of pair functions
           !$$--- allocate working space for interpolation
           allocate(WW(tFTable%NTAB), WA(tFTable%NTAB), WB(tFTable%NTAB), WC(tFTable%NTAB))
           do I=1, SimBox%NGROUP
              do J=1,SimBox%NGROUP
                !$$--- find out the table index for IFORCE
                 KTAB = 0
                 IFORCE = SimBox%PTYPE(I,J)
                 do K=1, size(FTable%FPAIR)
                    if(FTable%FPAIR(K) .EQ. IFORCE) then
                       KTAB = K
                      exit
                   end if
                 end do

                if(KTAB .GT. 0 )  then
                  !$$--- need this table
                  !$$    if alread registered, return
                   if(FTable%HASREG(KTAB).GT.0) cycle
                else
                  !$$--- DO not need this table
                  cycle
                end if

                do K=1, size(tFTable%FPAIR)
                   if(tFTable%FPAIR(K) .EQ. IFORCE) then
                      KTAB0 = K
                      exit
                   end if
                end do

                !$$--- interpolation for POTR
                TRMAX = maxval(tFTable%R)
                call SPLID1(tFTable%NTAB,tFTable%R, tFTable%POTR(KTAB0, 1:tFTable%NTAB), WW, IOP, 1,WA,WB,WC)
                do K=1, FTable%NTAB
                   T  = dble(K)*CSIV
                   TT = T*T
                   if(TT .le. TRMAX) then
                      call SPLID2(tFTable%NTAB, tFTable%R, tFTable%POTR(KTAB0, 1:tFTable%NTAB), WW, 1, TT,TAB)
                      FTable%POTR(KTAB,K) = TAB(1)
                   else
                      FTable%POTR(KTAB,K) = 0.D0
                   end if
                end do

                !$$--- interpolation for FPOTR
                call SPLID1(tFTable%NTAB,tFTable%R, tFTable%FPOTR(KTAB0, 1:tFTable%NTAB), WW, IOP, 1,WA,WB,WC)
                do K=1, FTable%NTAB
                   T = dble(K)*CSIV
                   TT = T*T
                   if(TT .le. TRMAX) then
                      call SPLID2(tFTable%NTAB, tFTable%R, tFTable%FPOTR(KTAB0, 1:tFTable%NTAB), WW, 1, TT,TAB)
                      FTable%FPOTR(KTAB,K) = TAB(1)
                   else
                      FTable%FPOTR(KTAB,K) = 0.D0
                   end if
                end do

                !$$--- interpolation for POTB
                call SPLID1(tFTable%NTAB,tFTable%R, tFTable%POTB(KTAB0, 1:tFTable%NTAB), WW, IOP, 1,WA,WB,WC)
                do K=1, FTable%NTAB
                   T  = dble(K)*CSIV
                   TT = T*T
                   if(TT .le. TRMAX) then
                      call SPLID2(tFTable%NTAB, tFTable%R, tFTable%POTB(KTAB0, 1:tFTable%NTAB), WW, 1, TT,TAB)
                      FTable%POTB(KTAB,K) = TAB(1)
                   else
                      FTable%POTB(KTAB,K) = 0.D0
                   end if
                end do

                !$$--- interpolation for FPOTB
                call SPLID1(tFTable%NTAB,tFTable%R, tFTable%FPOTB(KTAB0, 1:tFTable%NTAB), WW, IOP, 1,WA,WB,WC)
                do K=1, FTable%NTAB
                   T = dble(K)*CSIV
                   TT = T*T
                   if(TT .le. TRMAX) then
                      call SPLID2(tFTable%NTAB, tFTable%R, tFTable%FPOTB(KTAB0, 1:tFTable%NTAB), WW, 1, T*T,TAB)
                      FTable%FPOTB(KTAB,K) = TAB(1)
                   else
                      FTable%FPOTB(KTAB,K) = 0.D0
                   end if
                end do
                FTable%HASREG(KTAB) = 1

              end do
           end do
           deallocate(WW, WA, WB, WC)

           do I=1, SimBox%NGROUP
              do J=1,SimBox%NGROUP
                 IFORCE = SimBox%PTYPE(I,J)
                 do K=1, FTable%NKIND
                    if(IFORCE .EQ. FTable%FPAIR(K)) then    ! alread resgister
                       FTable%KPAIR(I, J) = K
                       exit
                    end if
                 end do
             end do
           end do

           !$$--- to perform interpolation of embedment function
           !$$--- allocate working space for interpolation
           allocate(WW(tFTable%NEMBD), WA(tFTable%NEMBD), WB(tFTable%NEMBD), WC(tFTable%NEMBD))

           !$$--- for imported table, we used the maxvals of the tables provided max value
           !$$    insted the statement in old version:
           !$$    FTable%RHOMX = maxval(FTable%POTB)*m_RHOMULT
           !$$
           FTable%RHOMX = (tFTable%RHO(tFTable%NEMBD)/dble(tFTable%NEMBD - 1))*dble(tFTable%NEMBD)
           FTable%RHOD  = FTable%RHOMX/dble(FTable%NEMBD)
           if(FTable%RHOD .le. 1.D-64) FTable%RHOD = 1.D0
           do I=1, SimBox%NGROUP
              !$$--- find out the table index for IFORCE
              KTAB = 0
              IFORCE = SimBox%PTYPE(I,I)
              do K=1, size(FTable%FPAIR1)
                 if(FTable%FPAIR1(K) .EQ. IFORCE) then
                    KTAB = K
                    exit
                 end if
              end do

              if(KTAB .GT. 0 )  then
                 !$$--- need this table
                 !$$    if alread registered, return
                  if(FTable%HASREG1(KTAB).GT.0) cycle
              else
                 !$$--- DO not need this table
                 cycle
              end if

              do K=1, size(tFTable%FPAIR1)
                 if(tFTable%FPAIR1(K) .EQ. IFORCE) then
                    KTAB0 = K
                    exit
                 end if
              end do

              !$$--- interpolation for FEMBD
              call SPLID1(tFTable%NEMBD, tFTable%RHO, tFTable%FEMBD(KTAB0, 1:tFTable%NEMBD), WW, IOP, 1,WA,WB,WC)
              do K=1, FTable%NEMBD
                 T = dble(K-1)*FTable%RHOD
                 call SPLID2(tFTable%NEMBD, tFTable%RHO, tFTable%FEMBD(KTAB0, 1:tFTable%NEMBD), WW, 1, T,TAB)
                 FTable%FEMBD(KTAB, K) = TAB(1)
              end do

              call SPLID1(tFTable%NEMBD, tFTable%RHO, tFTable%DFEMBD(KTAB0, 1:tFTable%NEMBD), WW, IOP, 1,WA,WB,WC)
              do K=1, FTable%NEMBD
                 T = dble(K-1)*FTable%RHOD
                 call SPLID2(tFTable%NEMBD, tFTable%RHO, tFTable%DFEMBD(KTAB0, 1:tFTable%NEMBD), WW, 1, T,TAB)
                 FTable%DFEMBD(KTAB, K) = TAB(1)
              end do

              FTable%HASREG1(KTAB) = 1
           end do
           deallocate(WW, WA, WB, WC)

           do I=1, SimBox%NGROUP
              IFORCE = SimBox%PTYPE(I,I)
              do K=1, FTable%NKIND1
                 if(IFORCE .EQ. FTable%FPAIR1(K)) then    ! alread resgister
                    FTable%KEMBD(I) = K
                    exit
                  end if
              end do
           end do

           call Release_ForceTable(tFTable)

     !*** print out information about the avalable potential
           call Print_Available_Force_ForceTable(6, FTable);
           call Print_Concerned_Force_ForceTable(6, FTable)
           if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
           if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)
           return

  end subroutine Register_Imported_ForceTable
  !**********************************************************************************

  !**********************************************************************************
   subroutine FS_EMEDMENTFUN(RHO,FRHO, DFRHO)
  !***  PURPOSE:   the embedment function of FS potential, used for debug
  !                EAM module
  !
  !     INPUT:     RHO:       the electron density function
  !                FRHO:      the embedment function
  !                FTable,    the differential of embedment function

  use MD_CONSTANTS
  implicit none
         !--- dummy vaiables
          real(KINDDF),intent(in)::RHO
          real(KINDDF),intent(out)::FRHO, DFRHO

              if(RHO .gt. 0.D0) then
                FRHO  =  -dsqrt(RHO)
                DFRHO =  -C_HALF/dsqrt(RHO)
              else
                FRHO  =  0.D0
                DFRHO =  0.D0
              end if

          return
   end subroutine FS_EMEDMENTFUN
  !**********************************************************************************

  end module MD_TYPEDEF_ForceTable
