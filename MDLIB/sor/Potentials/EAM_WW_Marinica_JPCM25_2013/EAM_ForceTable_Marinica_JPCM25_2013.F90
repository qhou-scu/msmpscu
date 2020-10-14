  module EAM_WW_ForceTable_Marinica_JPCM25_2013
  !**** DESCRIPTION: To generate the force table using the EAM potential
  !                  given by Marinica et al, JPCM25(2013)395502
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
  use EAM2_WW_Marinica_JPCM25_2013,only:EAM2_WW_NN=>NN_Func,     &
                                        EAM2_WW_RHO=>RHO_Func,   &
                                        EAM2_WW_EMBED=>EMBED_Func
  use EAM3_WW_Marinica_JPCM25_2013,only:EAM3_WW_NN=>NN_Func,     &
                                        EAM3_WW_RHO=>RHO_Func,   &
                                        EAM3_WW_EMBED=>EMBED_Func
  use EAM4_WW_Marinica_JPCM25_2013,only:EAM4_WW_NN=>NN_Func,     &
                                        EAM4_WW_RHO=>RHO_Func,   &
                                        EAM4_WW_EMBED=>EMBED_Func
  !-----------------------------------------------------------------------------------------
  use MiniUtilities, only:UpCase
  implicit none
      !--- dummy vaiables
      type(SimMDBox),    intent(in)   ::SimBox
      type(SimMDCtrl),   intent(in)   ::CtrlParam
      type(MDForceTable),intent(inout)::FTable
      integer,           optional     ::printout

      !--- local vaiables

      !*** register the table according to user suplied routines
      !--- local vaiables
      character*32::subtype

      !*** register the table according to user suplied routines
           subtype = SimBox%PotSubLibname(1:len_trim(SimBox%PotSubLibname))
           call UpCase(subtype)
           select case(subtype)
                  case default
                        call  Register_ForceTableProc("W","W", FTable, NNFORCE=EAM2_WW_NN, NEFORCE=EAM2_WW_RHO, EMBDF =EAM2_WW_EMBED,  &
                                    NOTE= "EAM2 (Marinica et al, JPCM25(2013)395502)")
                  case("EAM2")
                        call  Register_ForceTableProc("W","W", FTable, NNFORCE=EAM2_WW_NN, NEFORCE=EAM2_WW_RHO, EMBDF =EAM2_WW_EMBED,  &
                                    NOTE= "EAM2 (Marinica et al, JPCM25(2013)395502)")
                  case("EAM3")
                        call  Register_ForceTableProc("W","W", FTable, NNFORCE=EAM3_WW_NN, NEFORCE=EAM3_WW_RHO, EMBDF =EAM3_WW_EMBED,  &
                                    NOTE= "EAM3 (Marinica et al, JPCM25(2013)395502)")
                  case("EAM4")
                        call  Register_ForceTableProc("W","W", FTable, NNFORCE=EAM4_WW_NN, NEFORCE=EAM4_WW_RHO, EMBDF =EAM4_WW_EMBED,  &
                                    NOTE= "EAM4 (Marinica et al, JPCM25(2013)395502)")
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
  end module EAM_WW_ForceTable_Marinica_JPCM25_2013





