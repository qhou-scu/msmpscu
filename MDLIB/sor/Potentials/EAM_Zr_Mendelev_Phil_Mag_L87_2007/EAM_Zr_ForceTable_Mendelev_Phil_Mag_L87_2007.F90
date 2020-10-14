  module EAM_ZR_ForceTable_Mendelev_Phil_Mag_L87_2007
  !**** DESCRIPTION: To calculate the forces on atoms, the force is assumed in tight-binding FS form
  !                  HOU Qing, Oct, 2012
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

  use EAM_ZrZr_Mendelev_Phil_Mag_L87_2007,  only:ZRZR_NN_Interaction=>NN_Interaction, &
                                                 ZRZR_NE_Interaction=>NE_Interaction, &
                                                 ZrZr_EMBED=>EMBED_Func

  !-----------------------------------------------------------------------------------------


  implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in)   ::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(inout)::FTable
      integer,            optional     ::printout

      !--- local vaiables


      !*** register the table according to user suplied routines
           call Register_ForceTableProc("Zr", "Zr", FTable, NNFORCE= ZRZR_NN_Interaction, &
                                                   NEFORCE= ZRZR_NE_Interaction, &
                                                   EMBDF =  ZRZR_EMBED,          &
                                        NOTE= "Zr-Zr (Mendelev, Ackland,Phil.Mag.Lett.87(2007)349)")

      !****  create the force table
          call Create_Interaction_ForceTable(SimBox,CtrlParam,FTable)

      !*** print out information about the avalable potential
           if(.not.present(printout)) then
               call Print_Available_Force_ForceTable(6,FTable);
               call Print_Concerned_Force_ForceTable(6,FTable)
               if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
               if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)

           else
               if(printout)  then
                  call Print_Available_Force_ForceTable(6, FTable);
                  call Print_Concerned_Force_ForceTable(6, FTable)
                  if(gm_hFILELOG .gt. 0) call Print_Available_Force_ForceTable(gm_hFILELOG, FTable)
                  if(gm_hFILELOG .gt. 0) call Print_Concerned_Force_ForceTable(gm_hFILELOG, FTable)
               end if
           end if

           return

  end subroutine Register_Interaction_Table
  !**********************************************************************************
end module EAM_ZR_ForceTable_Mendelev_Phil_Mag_L87_2007
