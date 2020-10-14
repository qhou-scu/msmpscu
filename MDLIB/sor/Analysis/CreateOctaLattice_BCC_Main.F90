 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to create  octahedral lattice for a BCC crystal. The BCC lattice
 !                  is to be loaded from a configuration file.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                  CreateOctaLattice_Main.F90.exe -s "filename"
 !                  where:
 !                        filename  - the name of the file storing the coordinates of BCC lattices.
 !                                    The file could be created by, for example MakeBox tool, or
 !                                    a configuration file created by MDPSCU
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2018-01 (Hou Qing, Sichuan university)
 !
 !

 !**********************************************************************************
 Program CreateOctaLattice_BCC_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use CreateOctaLattice_BCC

 implicit none

       call APPSHELL_AddRecord( PRERECORD= Initialize_OctaSite, &
                                RECORDPROC=Record_OctaSite)

       call Single_ANALYSIS(0)

       stop
 End program CreateOctaLattice_BCC_Main
