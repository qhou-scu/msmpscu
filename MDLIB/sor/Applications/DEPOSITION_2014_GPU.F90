  Program DEPOSITION_EM_TB2010_Main
  use MD_SimBoxArray_AppShell_16_GPU
  use DEPOSITION_COMMON_2010

  !********************************************************************

  !--- local vairaibles
  integer::procid=0,numprocs=1

            call APPSHELL_AddRecord(PRERECORD=LoadControlParameter14_DEPOSITION)
            call APPSHELL_Main(numprocs,procid, INICONFIGPROC=DepositAtoms_DEPOSITION)

         stop
  end    program DEPOSITION_EM_TB2010_Main
