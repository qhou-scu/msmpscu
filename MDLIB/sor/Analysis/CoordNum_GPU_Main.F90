 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to calculate the coordination number of atoms
 !                  The calculation of the coordination number of atoms distribution is perfromed on GPUs.
 !                  Thus, cuda support is required for the computer.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       COORDNUMB_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at liest two lines need to be add:
 !
 !                    &AUXF_CNI  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_CNO filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                  KEYWORDS (common for most anaylysis tools):
 !                    &QUICKDAMP     indicate if damping to be performed before running analysis routine
 !                                   = 1 for yes, =0 for no (default)
 !                                   example: QUICKDAMP  1
 !
 !                    &JOBSEL        indicate which TESTs (ref to user guid) to be included in analysis.
 !                                   without this keyword, all TESTs will be analyzed.
 !
 !                                   usage: &JOBSEL  J0, J1, JS
 !                                          where J0 is the id of the frist TEST
 !                                          J1 is the id of the end  TEST
 !                                          JS is the intervale of the TESTs
 !
 !                                   example: &JOBSEL  1, 99, 2 indicating TESTs #1, #3, #5,...,#99 will
 !                                            be included for analyzing
 !
 !                    &CFGSEL        indicate which configurations to be included in analysis for included TESTs
 !                                   without this keyword, all configurations will be analyzed.
 !
 !                                   usage: &CFGSEL  C0, C1, CS
 !                                          where C0 is the id of the frist configuration in a test
 !                                          C1 is the id of the end configuration in a test
 !                                          CS is the intervale of the configurations
 !
 !                                   example: &CFGSEL  5, 100, 5 indicating configuration #5, #10, ...,#100
 !                                            will be included for analyzing
 !
 !                  KEYWORDS (specific):
 !                    &BONDLEN       indicate the input of bond length. One can define the diffrent bond
 !                                   length for different pair of atom types.
 !
 !                                   usage: &BOND  ty1 ,  ty2,   bond
 !                                          where ty1 and ty2 are the types of atoms.
 !                                          bond is bond length in lattice unit.
 !                                   Without this input, the bond length by atoms of ty1 and ty2
 !                                   will be their force cutoff distance.
 !
 !                    &AUXF_CNO     The file name of output data. If the filename is given hare,
 !                                  the filename2 (if given in SETUP file) will be replaced.
 !                                  usage: &AUXF_CNO  "fname"
 !
 !                  With the input file(s) are prepaired, the program is run on the command line:
 !
 !                  CoordNum_GPU.exe SETUP-file dev0  ndev
 !                  SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                  dev0        - the ID of the first GPU to be used
 !                  ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2013-05 (Hou Qing, Sichuan university)
 !                  update:     2014-11 (Hou Qing, Sichuan university)
 !
 Program CoordNum_GPU_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use COORDNUMB_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_CoordNumber_DEV,  &
                                RECORDPROC=Record_CoordNumber_TOOL,    &
                                AFTRECORD=Clear_WorkingArray_CoordNum_DEV)

       call Main_ANALYSIS(PROCESSID=0)

       stop
 end    program CoordNum_GPU_Main
