 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to simulate the effects of grain-boundary on the diffusion behavior
 !                  of embedment atoms. It is assumed that a boundary has been created (using, for example,
 !                  the program MakeTwin) in the initial substrate and the substrate is divided into two
 !                  parts in Z direction. Then, in the present program, two atom to be injected into the
 !                  two parts. The types of the two atoms are the last two type IDs.  For exmaple, the
 !                  substrate usually contains three types of atoms: type 1 is for the atoms on the interface,
 !                  type 2 is for the bulk atoms, type 3 is for the atom fixed on top and bottom of the
 !                  simulation box. The typoe of the embeded atoms are type 4 and type 5. Thus, in our
 !                  box-file, the NGROUP should be at least 5.
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  by HOU Qing, July., 2015
 !
  module Embedment_Atom_Interface_
  use MD_CONSTANTS
  implicit none

  contains
  !****************************************************************************************
  subroutine Initialize(SimBox,CtrlParam)
  !***  DESCRIPTION:  to load controal parameters and initialize the parameters
  !                   when the APPShell module MD_SimBoxArray_ToolShell_14_GPU to be used
  !     INPUT:        Simbox, the description of the simulation box
  !                   CtrlParam, the generic control parameters
   use MD_CONSTANTS
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
       !--- local
       integer::I
              !--- Do nothing at the present time

    return
  end subroutine Initialize
 !****************************************************************************************

 !****************************************************************************************
  subroutine Insert_Diffusor(SimBox, CtrlParam, RESTART)
  !***  PORPOSE: to initialize the configure of the box by emdedding an atom
  !              in a box
  !     INPUT:  SimBox0, the original substarte
  !             SimBox,  the box array to be created
  !             CtrlParam, the control parameter
  !             RESTART, indictor to indicating if restarting a new session
   use MD_SimboxArray
   use MD_TYPEDEF_SimMDCtrl
   use RAND32_MODULE

   implicit none
   !--- dummy variables and subroutines
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer                     ::RESTART

  !--- Local variables
  integer::I,NL, IL, NGROUP
  real(KINDDF)::RR, XP(1,3), ZMIN, ZMAX, THICK, BOUNDZ0,BOUNDZ1

           !---
           RR = SimBox(1)%RR
           NGROUP = SimBox(1)%NGROUP

          !--- change the values if other orietantion of boudary is used
           BOUNDZ0 = 0.224*3*RR
           BOUNDZ1 = 0.224*2*RR
           THICK   = 0.224*RR
           ZMIN = minval(SimBox(1)%XP(:,3), mask=((IAND(SimBox(1)%STATU,CP_STATU_FIXPOSZ) .EQ. 0) ) )
           ZMAX = maxval(SimBox(1)%XP(:,3), mask=((IAND(SimBox(1)%STATU,CP_STATU_FIXPOSZ) .EQ. 0) ) )
           NL   = (0.5D0*(ZMAX - ZMIN) - BOUNDZ0 - BOUNDZ1)/THICK

           !--- now we embeded an atom to the substrate
           IL = 0
           do I=1, SIZE(SimBox)
              ZMIN = minval(SimBox(I)%XP(:,3), mask=((IAND(SimBox(I)%STATU,CP_STATU_FIXPOSZ) .EQ. 0) ) )
              ZMAX = maxval(SimBox(I)%XP(:,3), mask=((IAND(SimBox(I)%STATU,CP_STATU_FIXPOSZ) .EQ. 0) ) )
              !--- add atom to the low part
              XP(1,1) = DRAND32()*SimBox(I)%ZL(1)+SimBox(I)%BOXLOW(1)
              XP(1,2) = DRAND32()*SimBox(I)%ZL(2)+SimBox(I)%BOXLOW(2)
              XP(1,3) = DRAND32()*(ZMAX - BOUNDZ0 - BOUNDZ1) + BOUNDZ0
              !XP(1,3) = BOUNDZ0 + IL*THICK

              call AddAtoms_SimMDBox(SimBox(I), 1, NGROUP-1, TYPEORDER=1, RXP=XP)
              !--- add atom to the upper part
              XP(1,1) = DRAND32()*SimBox(I)%ZL(1)+SimBox(I)%BOXLOW(1)
              XP(1,2) = DRAND32()*SimBox(I)%ZL(2)+SimBox(I)%BOXLOW(2)
              XP(1,3) = -( DRAND32()*(ZMAX - BOUNDZ0 - BOUNDZ1) + BOUNDZ0)
              !XP(1,1) = - XP(1,1)
              !XP(1,2) = - XP(1,2)
              !XP(1,3) = - XP(1,3)
              call AddAtoms_SimMDBox(SimBox(I), 1, NGROUP, TYPEORDER=1, RXP=XP)

              IL = IL + 1
              if(IL .gt. NL) IL = 0
           end do
       return
  end subroutine Insert_Diffusor
  !****************************************************************************************

 !****************************************************************************************
 end module Embedment_Atom_Interface_


 !****************************************************************************************
  Program Atom_Interface_GPU_Main
  use MD_SimBoxArray_AppShell_16_GPU
  !---- If user-define potential to be used, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
  !     use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
  !     use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table

  use Embedment_Atom_Interface_
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord(PRERECORD=Initialize)
       !---- If user-define potential to be generated, modify the code as following examplee
            ! call APPSHELL_Main(numprocs,procid, FORCETABLE=Reg1, POTTYPE="EAM_TYPE", INICONFIGPROC=IniConfig_EMBEDMENT)
       !---  else use internal potentials
       call APPSHELL_Main(numprocs,procid, INICONFIGPROC=Insert_Diffusor)

       stop
  End  program Atom_Interface_GPU_Main
