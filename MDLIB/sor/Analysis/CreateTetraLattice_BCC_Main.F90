 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to create  tetrahdral lattice for a BCC crystal. The BCC lattice
 !                  is to be loaded from a configuration file.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                  CreateTetraLattice_BCC.exe -s "filename"
 !                  where:
 !                        filename  - the name of the file storing the coordinates of BCC lattices.
 !                                    The file could be created by, for example MakeBox tool, or
 !                                    a configuration file created by MDPSCU
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-06 (Hou Qing, Sichuan university)
 !
 !

 !**********************************************************************************
 Program CreateTetraLattice_BCC_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use CreateTetraLattice_BCC

 implicit none

       call APPSHELL_AddRecord( PRERECORD= Initialize_TetraSite, &
                                RECORDPROC= Record_TetraSite)

       call Single_ANALYSIS(0)

       stop
 End program CreateTetraLattice_BCC_Main
