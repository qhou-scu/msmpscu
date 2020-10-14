!***  DESCRIPTION:  
!     This program is used to anaylsis the local pressure
!     using virial formular  

Program LOCAL_VIRIAL_EM_TB_Force_Table_Main
use MD_SimBoxArray_ToolShell_2010
use LOCAL_VIRIAL_EM_TB_Force_Table_2013
use EM_TB_ForceTable_WangJun_W_HE_2010
!use EM_TB_ForceTable_Ackland_Vitek_PRB41_10324

implicit none
 
    !**** to give calculation parameter
        call Main_ANALYSIS(PROCESSID=0, FORCETABLE=Register_Interaction_Table,  &
                                        PRERECORD=Init_CalLocalStress_Tool,     &
                                        RECORDPROC=RECORD_LocalStress_Tool,     &
                                        AFTERRECORD=RECORD_Average_LocalStress,   &
                                        STARTP=0, ENDP=0)
     
       stop 
end	program LOCAL_VIRIAL_EM_TB_Force_Table_Main
