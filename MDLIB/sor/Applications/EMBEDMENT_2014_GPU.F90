  Program EMBEDMENT_HE_W_2010_Main
  use MD_SimBoxArray_AppShell_16_GPU
  use EMBEDMENT_COMMON_2010
  !---- If user-define potential to be used, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
       use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
       use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord(PRERECORD=Initialize14_EMBEDMENT)

       !---- If user-define potential to be generated, modify the code as following examplee
            ! call APPSHELL_Main(numprocs,procid, FORCETABLE=Reg1, POTTYPE="EAM_TYPE", INICONFIGPROC=IniConfig_EMBEDMENT)
       !---  else use internal potentials
       call APPSHELL_Main(numprocs,procid, INICONFIGPROC=IniConfig_EMBEDMENT)

       stop
  End  program EMBEDMENT_HE_W_2010_Main
