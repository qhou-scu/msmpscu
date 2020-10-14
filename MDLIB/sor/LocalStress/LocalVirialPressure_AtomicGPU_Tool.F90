!***  DESCRIPTION:  
!     This program is used to calculate the centrosymmetry distribution
!     in a BCC system

Program LocalVirialPressure_Atomi_Tool_main
use MD_SimBoxArray_ToolShell_12_GPU
use LOCALVIRIALPRESSURE_ATOMIC_MODULE
use EM_TB_ForceTable_WangJun_W_HE_2010

implicit none
 
      call Main_ANALYSIS(PROCESSID=0, FORCETABLE=Register_Interaction_Table,        & 
	                                  PRERECORD=Init_Atomic_StressTensor_Tool,      &
                                      RECORDPROC=RECORD_Atomic_StressTensor_TOOL,   &
                                      AFTRECORD=CLEAR_LocalVirialPressure_DEV, STARTP=0, ENDP=1) 
     
      stop 
end	program LocalVirialPressure_Atomi_Tool_main