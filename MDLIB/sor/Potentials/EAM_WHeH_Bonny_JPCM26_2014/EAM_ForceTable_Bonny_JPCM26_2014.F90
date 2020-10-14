  module EAM_WHeH_ForceTable_Bonny_JPCM26_2014
  !**** DESCRIPTION: To generate the force table using the WHeH EAM potential
  !                  given by G Bonny, P Grigorev and D Terentyev, JPCM26(2014)485001
  !                  HOU Qing, Dec 9, 2014
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
  use EAM1_WHeH_Bonny_JPCM26_2014,only: WW_NN =>NN_WW_Func,  WW_RHO=>RHO_WW_Func, WW_EMBED=>EMBED_WW_Func,&

                                        EAM1_HeHe=>NN_HeHe_Func1,   &
                                        EAM1_HH  =>NN_HH_Func1,     &
                                        EAM1_HHe =>NN_HHe_Func1,    &
                                        EAM1_WHe =>NN_WHe_Func1,    &
                                        EAM1_WH  =>NN_WH_Func1,     &

                                        EAM2_HeHe=>NN_HeHe_Func2,   &
                                        EAM2_HH  =>NN_HH_Func2,  HH_EMBED => EMBED_HH_Func2,   &
                                        EAM2_HHe =>NN_HHe_Func2,    &
                                        EAM2_WHe =>NN_WHe_Func2,    &
                                        EAM2_WH  =>NN_WH_Func2,     &
                                        WH_RHO   =>RHO_WH_Func

  use MD_Pot_EAM_Utilities, only: RHOZero_Func, EMBEDZero_Func
  !-----------------------------------------------------------------------------------------
  use MiniUtilities, only:UpCase
  implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(inout)::FTable
      integer,           optional     ::printout

      !--- local vaiables
      character*32::subtype

      !*** register the table according to user suplied routines
           subtype = SimBox%PotSubLibname(1:len_trim(SimBox%PotSubLibname))
           call UpCase(subtype)
           select case(subtype)
                  case default
                    call  Register_ForceTableProc("W", "W", FTable, NNFORCE=WW_NN, NEFORCE=WW_RHO, EMBDF =WW_EMBED,  &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "W", FTable, NNFORCE=EAM1_WH,       &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "W", FTable, NNFORCE=EAM1_WHe,     &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")

                    call  Register_ForceTableProc("W","H", FTable, NNFORCE=EAM1_WH,        &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "H", FTable, NNFORCE=EAM1_HH,       &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He","H", FTable, NNFORCE=EAM1_HHe,      &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")

                    call  Register_ForceTableProc("W","He", FTable, NNFORCE=EAM1_WHe,      &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H","He", FTable, NNFORCE=EAM1_HHe,      &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "He", FTable, NNFORCE=EAM1_HeHe,   &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")

                    FTable%PotSubLibname = "EAM1"
                  case ("EAM1")
                    call  Register_ForceTableProc("W", "W", FTable, NNFORCE=WW_NN, NEFORCE=WW_RHO, EMBDF =WW_EMBED,  &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "W", FTable, NNFORCE=EAM1_WH,       &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "W", FTable, NNFORCE=EAM1_WHe,     &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")

                    call  Register_ForceTableProc("W","H", FTable, NNFORCE=EAM1_WH,        &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "H", FTable, NNFORCE=EAM1_HH,       &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He","H", FTable, NNFORCE=EAM1_HHe,      &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")

                    call  Register_ForceTableProc("W","He", FTable, NNFORCE=EAM1_WHe,      &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H","He", FTable, NNFORCE=EAM1_HHe,      &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "He", FTable, NNFORCE=EAM1_HeHe,   &
                                    NOTE= "(EAM1 in Bonny et al, JPCM26(2014)485001)")
                    FTable%PotSubLibname = "EAM1"

                  case ("EAM2")
                    call  Register_ForceTableProc("W", "W", FTable, NNFORCE=WW_NN, NEFORCE=WW_RHO, EMBDF=WW_EMBED,  &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "W", FTable, NNFORCE=EAM2_WH, NEFORCE=RHOZero_Func, &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "W", FTable, NNFORCE=EAM2_WHe,       &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")

                    call  Register_ForceTableProc("W", "H", FTable, NNFORCE=EAM2_WH, NEFORCE=WH_RHO, &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "H", FTable, NNFORCE=EAM2_HH, NEFORCE=RHOZero_Func, EMBDF=HH_EMBED, &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "H", FTable, NNFORCE=EAM2_HHe,  &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")

                    call  Register_ForceTableProc("W", "He", FTable, NNFORCE=EAM2_WHe,  &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("H", "He", FTable, NNFORCE=EAM2_HHe,  &
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc("He", "He", FTable, NNFORCE=EAM2_HeHe,&
                                    NOTE= "(EAM2 in Bonny et al, JPCM26(2014)485001)")
                    call  Register_ForceTableProc(10, FTable, NNFORCE=NULL_NNFORCE,     &
                                    NOTE= " NULL FORCE")
                    FTable%PotSubLibname = "EAM2"
            end select
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
  end module EAM_WHeH_ForceTable_Bonny_JPCM26_2014





