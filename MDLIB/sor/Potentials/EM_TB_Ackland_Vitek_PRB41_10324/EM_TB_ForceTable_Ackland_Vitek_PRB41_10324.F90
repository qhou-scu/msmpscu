  module EM_TB_ForceTable_Ackland_Vitek_PRB41_10324
  !**** DESCRIPTION: To calculate the forces on atoms, the force is assumed in tight-binding form of FS.
  !                  Thie version is adopted from an Cal_TB_EMPIRICAL_ForceTable.f90( version 2000)
  !                  HOU Qing, Mar, 2010
  !
  !                 REF___________________________________________________________________________________
  !                    (1) M.W.Finnis and J.E.Sinclair Phil.Mag. A50 (1984) 45
  !                    (2) G.J.Ackland and V.Vitek, Phys.Rev.B41, (1990) 10324

  !****
  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_CONSTANTS
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_ForceTable
  !-----------------------------------------------
    implicit none

    !***  the list of the force avalible **********************
      !integer, parameter, public::TB_EM_CuCu = 1       !
      !integer, parameter, public::TB_EM_AuAu = 2       !
      !integer, parameter, public::TB_EM_CuAu = 3       !
      !integer, parameter, public::TB_EM_AuCu = 4       !
    !***********************************************************

  contains

  !***************************************************************************************
  subroutine Register_Interaction_Table(SimBox, CtrlParam, FTable, printout)

  !***** to include the following interaction model ***********************
  use EM_TB_CuCu_Ackland_Vitek_PRB41_10324,only: CUCU_TB_EM_NN_Interaction=>TB_EM_NN_Interaction,    &
                                                 CUCU_TB_EM_NE_Interaction=>TB_EM_NE_Interaction

  use EM_TB_AuAu_Ackland_Vitek_PRB41_10324,only: AUAU_TB_EM_NN_Interaction=>TB_EM_NN_Interaction,    &
                                                 AUAU_TB_EM_NE_Interaction=>TB_EM_NE_Interaction

  use EM_TB_CuAu_Ackland_Vitek_PRB41_10324,only: CUAU_TB_EM_NN_Interaction=>TB_EM_NN_Interaction,    &
                                                 CUAU_TB_EM_NE_Interaction=>TB_EM_NE_Interaction,    &
                                                 AUCU_TB_EM_NN_Interaction=>TB_EM_NN_Interaction,    &
                                                 AUCU_TB_EM_NE_Interaction=>TB_EM_NE_Interaction
  !REMARK:   in the algorithm of G.J.Ackland and V.Vitek, the AuCU_NE_interaction
  !          is the same as CuAu_NE_interaction. This is not true if you use
  !          Johnson's algorithm.
  !************************************************************************

  implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(inout)::FTable
      integer,           optional     ::printout

      !--- local vaiables

      !*** register the table according to user suplied routines
           call  Register_ForceTableProc("Cu", "Cu", FTable, NNFORCE = CUCU_TB_EM_NN_Interaction, NEFORCE= CUCU_TB_EM_NE_Interaction,&
                                                        NOTE= "Cu-Cu(Ackland,Vitek, Phys.Rev.B41(1990)10324)")
           call  Register_ForceTableProc("Au", "Au", FTable, NNFORCE = AUAU_TB_EM_NN_Interaction, NEFORCE=AUAU_TB_EM_NE_Interaction,&
                                                        NOTE= "Au-Au(Ackland,Vitek, Phys.Rev.B41(1990)10324)")
           call  Register_ForceTableProc("Au", "Cu", FTable, NNFORCE = CUAU_TB_EM_NN_Interaction, NEFORCE=CUAU_TB_EM_NE_Interaction,&
                                                        NOTE= "Cu-Au(Ackland,Vitek, Phys.Rev.B41(1990)10324)")
           call  Register_ForceTableProc("Cu", "Au", FTable, NNFORCE = AUCU_TB_EM_NN_Interaction, NEFORCE=AUCU_TB_EM_NE_Interaction,&
                                                        NOTE= "Au-Cu(Ackland,Vitek, Phys.Rev.B41(1990)10324)")

      !****  create the force table
          call Create_Interaction_ForceTable(SimBox,CtrlParam,FTable)

      !*** print out information about the avalable potential
           if(.not.present(printout)) then
               call Print_Available_Force_ForceTable(6, FTable);
               call Print_Concerned_Force_ForceTable(6, FTable)
               if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
               if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)
           else
               if(printout)  then
                  call Print_Available_Force_ForceTable(6, FTable);
                  call Print_Concerned_Force_ForceTable(6,FTable)
                  if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
                  if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)
                end if
           end if

           return

  end subroutine Register_Interaction_Table
  !**********************************************************************************


  end module EM_TB_ForceTable_Ackland_Vitek_PRB41_10324





