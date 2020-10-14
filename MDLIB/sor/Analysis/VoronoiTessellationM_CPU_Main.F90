 !***  DESCRIPTION:
 !     This program is used to calculate the centrosymmetry distribution
 !     in a BCC system

 Program VoronoiTessllation_Tool_main
 use MD_SimBoxArray_ToolShell_12_GPU
 use VoronoiTessllation_2013

 implicit none

       !call Main_ANALYSIS(PROCESSID=0, PRERECORD=Initialize_CoordNumber_DEV, &
       !                                RECORDPROC=RECORD_CoordNumber_TOOL,   &
       !                                AFTRECORD=CLEAR_WorkingArray_CoordNum_DEV)

       call Main_ANALYSIS(PROCESSID=0,  &
                                       RECORDPROC=RECORD_VoronoiTesselation_TOOL,   &
                                       AFTRECORD=CLEAR_VoronoiTessllation, STARTP=7, ENDP=7)

       stop
 end    program VoronoiTessllation_Tool_main
