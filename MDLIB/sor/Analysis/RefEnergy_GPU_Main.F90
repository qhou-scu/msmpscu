 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to create the displacement map of configuration
 !
 !**** ALGORITHM:
 !                  A configuration is initially loaded in and taken as a reference configuration (RC).
 !                  Each "atom" in the reference configuration (usually a perfect crystgal) is a site.
 !                  The MD configuration is laterlly loaded. The sites in the reference configuration
 !                  are scaned to find out the sites that are closest to the atoms and in the same type
 !                  of the atoms . The displacement vectors are define as the related positions of
 !                  the atoms to the their closest sites.
 !
 !                  The search of the closest sites of atoms is perfromed on GPUs. Thus, cuda
 !                  support is required for the computer.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       RefDisplace_GPU_Main.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one line needs to be add:
 !
 !                    &AUXF_REFEGIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_REFEGOUT filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                  &REFCFG        indication of the filename of the reference configurqation.
 !                                 if the reference configuration is not assigned explicitly,
 !                                 the configuration at zero time step will be the reference
 !                                 configuration.
 !
 !                                 usage:   &REFCFG filename
 !
 !                                 example: &REFCFG "Myreference.crystal"
 !
 !
 !                 &AUXF_DVOUT    indicating of filename of for outputing results. If the filename is
 !                                given here, the filename2 (if given in SETUP file) will be replaced.
 !
 !                  If the filename of the reference configuration is given in the control file described
 !                  above, the third file storing the reference configuration needs to be prepaired. The
 !                  the file could be in the format of "&CFGXYZ", "&BOXCFG14" or "t-x-y-z" file whithout
 !                  header.
 !
 !                  With the input file(s) are ready, the program can be run on command line:
 !
 !                  RefDisplace_GPU.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-06 (Hou Qing, Sichuan university)
 !
 !
 Program REFENERGY_TOOL_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use RefEnergy_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD  =PreRecord_RefEnergy ,   &
                                BEFONETEST =BefOneTest_RefEnergy,   &
                                RECORDPROC =Record_RefEnergy,       &
                                AFTRECORD  =AftRecord_RefEnergy)

       call Main_ANALYSIS(0)

       stop
 End program REFENERGY_TOOL_Main
