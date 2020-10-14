  module EM_TB_ForceTable_WangJun_W_HE_2010
  !**** DESCRIPTION: To calculate the forces on atoms, the force is assumed in tight-binding form of Finnis.
  !                  Thie version is adopted from an Cal_TB_EMPIRICAL_ForceTable.f90( version 2000)
  !                  Adapted from Wangjun's program by by HOU Qing, Mar, 2010
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
  !REMARK:   There are serveral version for He-W interaction supplied by WANG Jun.
  !          The choice of HE-W potential should be associated with the choice of W-W potential.
  !          Pot_ZBL_ABINTIO_He_W_614 should be associated with FS_PM_WW for W-W potential
  !          Pot_ZBL_ABINTIO_He_W_818 should be asscociated with POT_FS_ACK_WW for W-W potential
  !
  !          It has been noted that EXP6 potential is better for LJ potential for  HE-HE interaction

  !-----------------------------------------------------------------------------------------
  !-- combination 1: Pot_ZBL_ABINTIO_He_W_614, with FS_PM_WW  ------------------------
  !use POT_LJ_HeHe,only:HeHe_NN_Interaction=>NN_Interaction,     &
  !                     HeHe_NE_Interaction=>NE_Interaction
  !use Pot_EXP6_HeHe,only:HeHe_NN_Interaction=>NN_Interaction,     &
  !                     HeHe_NE_Interaction=>NE_Interaction

  !use POT_FS_PM_WW,only:WW_TB_CR_NN_Interaction=>TB_CR_NN_Interaction,     &
  !                      WW_TB_CR_NE_Interaction=>TB_CR_NE_Interaction

  !use Pot_ZBL_ABINTIO_He_W_614,only:WHe_ZBL_NN_Interaction=>ZBL_ABINTIO_NN_Interaction, &
  !                              WHe_ZBL_NE_Interaction=>ZBL_ABINTIO_NE_Interaction,     &
  !                              HeW_ZBL_NN_Interaction=>ZBL_ABINTIO_NN_Interaction,     &
  !                              HeW_ZBL_NE_Interaction=>ZBL_ABINTIO_NE_Interaction
  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  !-- combination 2: Pot_ZBL_ABINTIO_He_W_818, with FS_ACK_WW and EXP6 for He-He
  !use Pot_EXP6_HeHe,only:HeHe_NN_Interaction=>NN_Interaction,     &
  !                       HeHe_NE_Interaction=>NE_Interaction

  !use POT_FS_ACK_WW,only:WW_TB_CR_NN_Interaction=>TB_CR_NN_Interaction,     &
  !                       WW_TB_CR_NE_Interaction=>TB_CR_NE_Interaction

  !use Pot_ZBL_ABINTIO_He_W_818,only:WHe_ZBL_NN_Interaction=>ZBL_ABINTIO_NN_Interaction, &
  !                              WHe_ZBL_NE_Interaction=>ZBL_ABINTIO_NE_Interaction,     &
  !                              HeW_ZBL_NN_Interaction=>ZBL_ABINTIO_NN_Interaction,     &
  !                              HeW_ZBL_NE_Interaction=>ZBL_ABINTIO_NE_Interaction

  !-----------------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------------------
  !-- combination 3: Pot_ZBL_ABINTIO_He_W_930, with FS_ACK_WW and EXP6 for He-He
  use Pot_ZBL_EXP6_HeHe,only:HeHe_NN_Interaction=>NN_Interaction,     &
                             HeHe_NE_Interaction=>NE_Interaction

  use POT_FS_ACK_WW,only:WW_TB_NN_Interaction=>TB_NN_Interaction,     &
                         WW_TB_NE_Interaction=>TB_NE_Interaction

  use Pot_ZBL_ABINTIO_He_W_930,only:WHe_ZBL_NN_Interaction=>ZBL_ABINTIO_NN_Interaction, &
                                WHe_ZBL_NE_Interaction=>ZBL_ABINTIO_NE_Interaction,     &
                                HeW_ZBL_NN_Interaction=>ZBL_ABINTIO_NN_Interaction,     &
                                HeW_ZBL_NE_Interaction=>ZBL_ABINTIO_NE_Interaction

  use POT_FS_PM_WW,only:WW_TB_DND_NN_Interaction=>TB_NN_Interaction,     &
                        WW_TB_DND_NE_Interaction=>TB_NE_Interaction
  !
  !-----------------------------------------------------------------------------------------

  implicit none
      !--- dummy vaiables
      type(SimMDBox),      intent(in)::SimBox
      type(SimMDCtrl),     intent(in)::CtrlParam
      type(MDForceTable),  intent(inout)::FTable
      integer,             optional::printout

      !--- local vaiables


      !*** register the table according to user suplied routines
           if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "FS_TYPE") then
             call  Register_ForceTableProc(1, FTable, NNFORCE=WW_TB_NN_Interaction, NEFORCE= WW_TB_NE_Interaction,         &
                                    NOTE= "W-W (Ackland, Thetford, Phil.Mag.56(1987)15)")
             call  Register_ForceTableProc(2, FTable, NNFORCE=HeHe_NN_Interaction, NEFORCE= HeHe_NE_Interaction,           &
                                    NOTE= "He-He (EXP6 potential)")
             call  Register_ForceTableProc(3, FTable, NNFORCE=WHe_ZBL_NN_Interaction, NEFORCE= WHe_ZBL_NE_Interaction,     &
                                    NOTE= "He-W")
             call  Register_ForceTableProc(4, FTable, NNFORCE=WHe_ZBL_NN_Interaction, NEFORCE= WHe_ZBL_NE_Interaction,     &
                                    NOTE= "W-He")
             call  Register_ForceTableProc(5, FTable, NNFORCE=WW_TB_DND_NN_Interaction, NEFORCE= WW_TB_DND_NE_Interaction, &
                                    NOTE= "W-W (Derlet et al, Phys.Rev.B76(2007)054107)")

           else if(FTable%PotType(1:len_trim(FTable%PotType)) .eq. "EAM_TYPE") then
             call  Register_ForceTableProc(1, FTable, NNFORCE=WW_TB_NN_Interaction, NEFORCE= WW_TB_NE_Interaction,         &
                                    EMBDF = FS_EMEDMENTFUN,                                                         &
                                    NOTE= "W-W (Ackland, Thetford, Phil.Mag.56(1987)15)")
             call  Register_ForceTableProc(2, FTable, NNFORCE=HeHe_NN_Interaction, NEFORCE= HeHe_NE_Interaction,           &
                                    EMBDF = FS_EMEDMENTFUN,                                                         &
                                    NOTE= "He-He (EXP6 potential)")
             call  Register_ForceTableProc(3, FTable, NNFORCE=WHe_ZBL_NN_Interaction, NEFORCE= WHe_ZBL_NE_Interaction,     &
                                    EMBDF = FS_EMEDMENTFUN,                                                         &
                                    NOTE= "He-W")
             call  Register_ForceTableProc(4, FTable, NNFORCE=WHe_ZBL_NN_Interaction, NEFORCE= WHe_ZBL_NE_Interaction,     &
                                    EMBDF = FS_EMEDMENTFUN,                                                         &
                                    NOTE= "W-He")
             call  Register_ForceTableProc(5, FTable, NNFORCE=WW_TB_DND_NN_Interaction, NEFORCE= WW_TB_DND_NE_Interaction, &
                                    EMBDF = FS_EMEDMENTFUN,                                                         &
                                    NOTE= "W-W (Derlet et al, Phys.Rev.B76(2007)054107)")
           end if

      !****  create the force table
           if(present(printout))  then
              !***  if you need the detailed information of force, you can putout the table to a file:
              call Create_Interaction_ForceTable(SimBox,CtrlParam,FTable, FNAME="FORCETABLE.DAT")
           else
              call Create_Interaction_ForceTable(SimBox,CtrlParam,FTable)
           end if

      !*** print out information about the avalable potential
           call Print_Available_Force_ForceTable(6, FTable);
           call Print_Concerned_Force_ForceTable(6, FTable)
           if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
           if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)
           return

  end subroutine Register_Interaction_Table
  !**********************************************************************************


  end module EM_TB_ForceTable_WangJun_W_HE_2010





