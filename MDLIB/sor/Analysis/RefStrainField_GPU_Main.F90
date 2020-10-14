 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to energy spectrum of configurations created by MD simulations
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one line needs to be add:
 !
 !                    &AUXF_REFSTRIANIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_REFSTRIANOUT filename2
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
 !                 &AUXF_REFSTRIANOUT  indicating of filename of for outputing results. If the filename is
 !                                     given here, the filename2 (if given in SETUP file) will be replaced.
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
 Program REFDISPLACE_TOOL_MAIN
 use MD_SimBoxArray_ToolShell_14_GPU
 use RefStrainField_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize,                 &
                                RECORDPROC=Record_StrainField_TOOL,   &
                                AFTRECORD=Clear_WorkingArray_DEV)

       call Main_ANALYSIS(0)

       stop
 End program REFDISPLACE_TOOL_MAIN
