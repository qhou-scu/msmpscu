  module DEPOSITION_COMMON_2010
  !***  DESCRIPTION:
  !     To simulate the depostion of atoms on surface
  !
  !    Adopted by HOU Qing, Otc, 2010

  !***  The modules included ******************************************
  use MD_CONSTANTS
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use DEPOSITION_TypeDef_CtrlParam_2010

  !*** the application specific module ********************************
  implicit none

     type(DepositionCtrlParam),private::m_DepCtrlParam
     character(len=13),parameter, private::mp_FTAGI="&AUXF_DEPOSIT"

  contains
  !**************************************************************************
  subroutine Initialize_DEPOSITION(fname,CtrlParam)
  !***  PURPOSE: to load control parameters for depostion calculations
  !****
  !     INPUT: fname      , the file name storing the parameters
  !
  !     OUTPUT:CtrlParam,   the control parameters for deposition
  !
  implicit none
      !----   DUMMY Variables
       character*(*)::fname
       type(DepositionCtrlParam)::CtrlParam

       integer::hFile

         !*** to load controling parameters
          write(*,*) "!**** Load depostion parameters from:"
          write(*,*) "!**** ",fname(1:len_trim(fname))
          call AvailableIOUnit(hFile)
          open(hFile, file = fname, status='old')
          call  Load_Parameter_DepositionCtrlParam(hFile, CtrlParam)
          close(hFile)
        return
   end subroutine Initialize_DEPOSITION
  !****************************************************************************************

  !****************************************************************************************
  SUBROUTINE DepositOneAtom_DEPOSITION(SimBox0,SimBox,CtrlParam,DepCtrlParam)
  !***  PORPOSE: to deposite an atom by place the atom above the surface
  !              with a distance about the cutoff range or initialize an atom
  !              inside the box with given kinetic energy
  !
  !     INPUT:   SimBox0,     the intial substrate box
  !              CtrlParam,   the control parameters for whole box
  !              DepCtrlParam,the deposition parameters
  !     OUTPUT:  SimBox,      the box with atoms to be deposition in position
  !
  use RAND32_MODULE
  IMPLICIT NONE
   type(SimMDBox), intent(in),optional::SimBox0
   type(SimMDBox)                     ::SimBox
   type(SimMDCtrl)                    ::CtrlParam
   type(DepositionCtrlParam)          ::DepCtrlParam

  !Local variables
   integer::ITYP, IP
   real(KINDDF)::POS(1,3)=C_ZERO, VEL(1,3)=C_ZERO, DIR(3), ZMAX, ZMIN, CUTOFF, KE, VEL0

          !*** copy the initial substrate back to active box
          ! this was used in old version for
          ! DEPOSITION_COMMOM_10_GPU.f90 and DEPOSITION_COMMOM_10_CPU.f90
          if(present(SimBox0) ) then
             call Copy_SimMDBox(SimBox0, SimBox)
          end if

          if(DepCtrlParam%NPRT .LE. 0 .OR. DepCtrlParam%nType .LE. 0) return


          !*** Generate the particle
          call GenerateParticle_DepositionCtrlParam(DepCtrlParam, ITYP, KE, DIR)
          KE = KE*1000 ! from kEv to eV

          !*** then we create the velocity of the incident atom
            VEL0 = DSQRT(C_TWO*KE*CP_EVERG/SimBox%CM(ITYP))
            VEL(1,1:3) = VEL0*DIR(1:3)

          !*** to place the atom at right place according to depsotion style
             select case(DepCtrlParam%DEPSTYLE)
                    case(CP_DEP_STYPE_EXT)
                        !*** to place adatom to the surface
                        !    first we look for the surface coordinate of the substrate
                        ZMAX = maxval(SimBox%XP(1:SimBox%NPRT,3),                                      &
                                      mask=IAND(SimBox%STATU,CP_STATU_ACTIVE).EQ.CP_STATU_ACTIVE .AND. &
                                                SimBox%ITYP.EQ.1 )
                        CUTOFF = maxval(CtrlParam%RU*0.99)
                        POS(1,1) = DRAND32()*SimBox%ZL(1)+SimBox%BOXLOW(1)
                        POS(1,2) = DRAND32()*SimBox%ZL(2)+SimBox%BOXLOW(2)
                        POS(1,3) = ZMAX+CUTOFF
                        call AddAtoms_SimMDBox(SimBox,1,ITYP,TYPEORDER=1,RXP=POS, RXP1=VEL)
                        SimBox%EKIN(SimBox%NPRT) = KE*CP_EVERG
                    case(CP_DEP_STYPE_INT)
                        !*** to place an atom in the box
                        !--- randomly select a particle of type ITYP in the box
                        !    and then set up its kinetic energy
                        IP = SimBox%NA(ITYP)*DRAND32()+0.000001
                        IP = SimBox%IPA(ITYP) + IP
                        SimBox%XP1(IP,1:3) = VEL(1,1:3)
                        SimBox%EKIN(IP)    = KE*CP_EVERG
             end select

             !$$--- for externel depsotion, the position of surface
             !$$    may be change because of reflection of transimittion
             !$$    we should check and change the bozsize, oherwise
             !$$    there will be erro in calculating neighbore
             !ZMIN = minval(SimBox%XP(1:SimBox%NPRT,3))*1.1
             !ZMAX = maxval(SimBox%XP(1:SimBox%NPRT,3))*1.1
             !ZMAX = max(ZMAX, SimBox%BOXUP(3))
             !ZMIN = min(ZMIN, SimBox%BOXLOW(3))
             !SimBox%ZL(3) = ZMAX - ZMIN
             !SimBox%BOXLOW(3) = ZMIN
             !SimBox%BOXUP(3)  = ZMAX
             if(SimBox%XP(SimBox%NPRT,3) .GT. SimBox%BOXUP(3)) then
                SimBox%BOXUP(3) = SimBox%BOXUP(3) + SimBox%RR
             end if

             if(SimBox%XP(SimBox%NPRT,3) .LT. SimBox%BOXLOW(3)) then
                SimBox%BOXLOW(3) = SimBox%BOXLOW(3) - SimBox%RR
             end if
             SimBox%ZL(3)    = SimBox%BOXUP(3) - SimBox%BOXLOW(3)
          return
  END SUBROUTINE DepositOneAtom_DEPOSITION
  !****************************************************************************************

  !****************************************************************************************
  subroutine LoadControlParameter12_DEPOSITION(SimBox,CtrlParam, FORCETABLE)
  !***  DESCRIPTION:  to load controal parameters and initialize the parameters
  !     INPUT:        Simbox, the description of the simulation box
  !                   CtrlParam, the generic control parameters
   use MD_CONSTANTS
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam
       OPTIONAL        ::FORCETABLE
       external        ::FORCETABLE

   !--- local variables
       call LoadControlParameter14_DEPOSITION(SimBox,CtrlParam)
    return
  end subroutine LoadControlParameter12_DEPOSITION
  !****************************************************************************************

  !****************************************************************************************
  subroutine LoadControlParameter14_DEPOSITION(SimBox,CtrlParam)
  !***  DESCRIPTION:  to load controal parameters and initialize the parameters
  !     INPUT:        Simbox, the description of the simulation box
  !                   CtrlParam, the generic control parameters
   use MD_CONSTANTS
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam

   !--- local variables
       integer::I


      !*** to load control parameters specific for deposition calculations
           do I=1, size(CtrlParam%f_tag)
               if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                  exit
               end if
           enddo
           if(I.gt. size(CtrlParam%f_tag)) goto 100
           if(len_trim(CtrlParam%f_others(I)) .le. 0) goto 100

            call Initialize_DEPOSITION(CtrlParam%f_others(I),m_DepCtrlParam)
      !*** to load control parameters specific for embeddment calculations
           call Print_Parameter_DepositionCtrlParam(6, m_DepCtrlParam)
           return
  100      write(*,*)"MDPSCU error: control file for deposition calculations  is not given"
           write(*,*)"              check the setup file for keyword ", mp_FTAGI
           stop
    return
  end subroutine LoadControlParameter14_DEPOSITION
  !****************************************************************************************

!****************************************************************************************
  subroutine DepositAtoms_DEPOSITION(SimBox, CtrlParam, RESTART)
  !***  PORPOSE: to deposit an atom in a box
  !     INPUT:  SimBox,  the box array to be created
  !             CtrlParam, the control parameter
  !             RESTART, indictor to indicating if restarting a new session
   use MD_CONSTANTS
   use MD_Globle_Variables
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer::RESTART
  !--- local varibales
       integer::I
           DO I=1, size(SimBox)
              if(RESTART .EQ. 0) then
                 call DepositOneAtom_DEPOSITION(SimBox=SimBox(I),  CtrlParam=CtrlParam, DepCtrlParam=m_DepCtrlParam)
              else
                 call AddAtoms_SimMDBox(B=SimBox(I), N=RESTART, ITYPE=SimBox(I)%NGROUP, TYPEORDER=1)
              end if
           END DO
           return

  end subroutine DepositAtoms_DEPOSITION
  !****************************************************************************************

  end module DEPOSITION_COMMON_2010
