 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to calculate the centrosymmetry distribution in a BCC system
 !
 !**** ALGORITHM:
 !
 !
 !                  The calculation of the centrosymmetry distribution is perfromed on GPUs. Thus,
 !                  cuda support is required for the computer.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       CentroSymParam_for_BCC_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at liest two lines need to be add:
 !
 !                    &AUXF_CSPI  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_CSPO filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                    &AXSISX    the coordinate axsis X of simulation box, expressed by
 !                               Miller crystal index. Without this input, X axsis will be (1,0,0)
 !                    &AXSISY    the coordinate axsis Y of simulation box, expressed by
 !                               Miller crystal index. Without this input, Y axsis will be (0,1,0)
 !                    &AXSISZ    the coordinate axsis Z of simulation box, expressed by
 !                               Miller crystal index. Without this input, Z axsis will be (0,0,1)
 !                    &CPLEVELS  NL, level1, levl2, ...
 !                               where NL is the number of thresholds of CSP values, which
 !                               can be used to distinguish atoms. This could be useful
 !                               for visualization. Without this input, the atoms will not be
 !                               distinguished by CSP values.
 !                    &CAL_CPAV  the time steps for output average CSP, if the parameter < 0
 !                               no average of CSPs will be output
 !
 !                    &AUXF_CSPO filename for outputing results. If the filename is given hare,
 !                               the filename2 (if given in SETUP file) will be replaced.
 !
 !                  With the input file(s) are prepaired, the program is run on the command line:
 !
 !                  CentroSymParam_for_BCC_GPU.exe SETUP-file dev0  ndev
 !                  SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                  dev0        - the ID of the first GPU to be used
 !                  ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2013-05 (Hou Qing, Sichuan university)
 !                  update:     2014-11 (Hou Qing, Sichuan university)
 !
 !


 Program CentroSymParam_for_BCC_MAIN
 use MD_SimBoxArray_ToolShell_14_GPU
 use CentroSymParam_for_BCC_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_CentroSymmetry_DEV, &
                                RECORDPROC=Record_CentroSymmetryTool,    &
                                AFTRECORD=Clear_CentroSymmetry_DEV)

       call Main_ANALYSIS(PROCESSID=0)

       stop
 end    program CentroSymParam_for_BCC_MAIN
