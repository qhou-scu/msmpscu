!***  DESCRIPTION:  
!     This program is used to anaylsis the local stress
!     using control-volume method

Program LocalCVStress_EM_TB_ForceTable_Main
use MD_SimBoxArray_ToolShell_2010
use LOCAL_MOPSTRESS_EM_TB_Force_Table_2013
use EM_TB_ForceTable_WangJun_W_HE_2010
!use EM_TB_ForceTable_Ackland_Vitek_PRB41_10324

implicit none
 
    !**** to give calculation parameter
    
        !if(COMMAND_ARGUMENT_COUNT() .ge. 2) then
        !   call get_command_argument(2,m_OUTPUTFILE)
        !end if 

        call Main_ANALYSIS(PROCESSID=0, FORCETABLE=Register_Interaction_Table,  &
                                        PRERECORD=Init_CalLocalStress_Tool,     &
                                        RECORDPROC=RECORD_LocalStress_Tool,     &
                                        AFTERRECORD=RECORD_Average_LocalStress,   &
                                        STARTP=2, ENDP=2) !STARTP=2 )
     
       stop 
end	program LocalCVStress_EM_TB_ForceTable_Main
