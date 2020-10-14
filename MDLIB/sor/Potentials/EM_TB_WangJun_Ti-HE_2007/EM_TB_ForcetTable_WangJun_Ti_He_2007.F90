 module EM_TB_ForcetTable_WangJun_Ti_He_2007
  !**** DESCRIPTION: To register the forces table of between He-Ti.
  !                  The force is assumed in tight-binding form of Finnis.
  !                  The fitting parameters for Ti-Ti used are from:
  !                  F.Cleri and V.Rosato, Phys.Rev.B48, (1993)22.
  !
  !                  Adapted from Wangjun's program by HOU Qing, Jan, 2015
  !
  !    REFERENCE:    F.Cleri and V.Rosato, Phys.Rev.B48, (1993)22.
  !
  !    SEE ALSO:      MD_TYPEDEF_ForceTable_10.f90
  !
  !****
  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_ForceTable
  !-----------------------------------------------
    implicit none


  contains

  !***************************************************************************************
  subroutine Register_Interaction_Table(SimBox, CtrlParam, FTable, printout)
  !***  PURPOSE:   to initialize the force table to be used in force calculations
  !
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !                FTable, the force tabled to be registered
  !
  !     OUTPUT:    FTable, the registered table
  !***** to include the following interaction model ***********************
  use Pot_TB_CR_TiTi, only:TiTi_TB_NN_Interaction=>TB_NN_Interaction,     &
                           TiTi_TB_NE_Interaction=>TB_NE_Interaction

  use Pot_LJ_He_Ti,   only:HeTi_TB_NN_Interaction=>NN_HeTi_Interaction,  &
                           HeTi_TB_NE_Interaction=>NE_HeTi_Interaction,  &
                           HeHe_NN_Interaction=>NN_HeHe_Interaction,     &
                           HeHe_NE_Interaction=>NE_HeHe_Interaction

  !-----------------------------------------------------------------------------------------

  implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(inout)::FTable
      integer,           optional     ::printout

      !--- local vaiables

      !*** register the table according to user suplied routines
           call  Register_ForceTableProc("Ti", "Ti", FTable, NNFORCE= TiTi_TB_NN_Interaction, NEFORCE = TiTi_TB_NE_Interaction,&
                                    NOTE= "Ti-Ti(Cleri, Rosato, Phys.Rev.B48,(1993)22)")

           call  Register_ForceTableProc("He", "He", FTable, NNFORCE=HeHe_NN_Interaction,      &
                                    NOTE= "He-He (EXP6 potential)")


           call  Register_ForceTableProc("Ti", "He", FTable, NNFORCE= HeTi_TB_NN_Interaction,  &
                                    NOTE= "Ti-He(Wang Jun et al,  J.Appl.Phys.102(2007)093510)")

           call  Register_ForceTableProc("He", "Ti", FTable, NNFORCE= HeTi_TB_NN_Interaction,  &
                                    NOTE= "He-Ti(Wang Jun et al,  J.Appl.Phys.102(2007)093510)")


      !****  create to force table
           if(present(printout))  then
              call Create_Interaction_ForceTable(SimBox,CtrlParam,FTable, FNAME="FORCETABLE.DAT")
           else
      !***  if you need the detailed information of force, you can putout the table to a file:
              call Create_Interaction_ForceTable(SimBox,CtrlParam,FTable)
           end if

      !*** print out information about the avalable potential
           call Print_Available_Force_ForceTable(6, FTable);
           call Print_Concerned_Force_ForceTable(6,FTable)
           if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
           if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)

           return
  end subroutine Register_Interaction_Table
  !**********************************************************************************

 end module EM_TB_ForcetTable_WangJun_Ti_He_2007





