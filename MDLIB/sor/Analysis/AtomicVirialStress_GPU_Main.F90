 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to calculate the atomic virial stress of atoms. The volume of
 !                  atoms are determined b y Voronoi tessellation. The stress in cells of given size and
 !                  the stress in sphere region with given center and radius can be also calculated. The volume
 !                  of a cell or a sphere in the summation of the volume of atoms contained in the cell or
 !                  the sphere.
 !                  The calculation of the atomic stress atoms is perfromed on GPUs. Thus, cuda support is
 !                  required for the computer.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       AtomicVirialStress_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at liest two lines need to be add:
 !
 !                    &AUXF_AVSI  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_AVSO filename2
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
 !***********************************************************************
 !                  KEYWORDS (specific):
 !                    &ASTRESS       indicating if the atomic virial stress of atoms will be output.
 !
 !                                   usage: &ASTRESS out
 !                                   where out=1 for yes, out=0 for no.
 !                                   By default , the parameter is 0.
 !
 !                    &ASTRESSAV    indicating if the average atomic virial stress of atoms will be output.
 !
 !                                  usage: &ASTRESSAV  num
 !                                  where num is the number of configurations that the average of
 !                                  atomic stress will be perofrmed on. If num =0, no average stress
 !                                  will be calculated.
 !                                  By default , the parameter is 0.
 !
 !                    &CSTRESS      indicating if the atomic virial stress in cells will be output.
 !
 !                                  usage: &CSTRESS out, cx, cy, cz
 !                                  where out=1 for yes, out=0 for no. cx, cy, cz is the size of cells in LU.
 !                                  By default , our = 0, cx=cy=cz =1 LU
 !
 !                    &CSTRESSAV    indicating if the average atomic virial stress in cells will be output.
 !
 !                                  usage: &CSTRESSAV  num
 !                                  where num is the number of configurations that the average will
 !                                  be perofrmed on. If num =0, no average stress will be calculated.
 !                                  By default , the parameter is 0.
 !
 !                    &SSTRESS      indicating if the atomic virial stress in spheres in a sphere region will be output.
 !
 !                                  usage: &SSTRESS out, cx, cy, cz, r, num
 !                                  where out=1 for yes, out=0 for no. (cx, cy, cz) in LU is the center of the sphere region.
 !                                  r is the radius of the region, num is the number of shperes that the sphere region contain.
 !                                  The radia of the ith sphere will be (r/num)*i.
 !                                  By default , our = 0, cx=cy=cz=0, r=5, num = 10
 !
 !                    &SSTRESSAV    indicating if the average atomic virial stress in spheres will be output.
 !
 !                                  usage: &SSTRESSAV  num
 !                                  where num is the number of configurations that the average will
 !                                  be perofrmed on. If num =0, no average stress will be calculated.
 !                                  By default , the parameter is 0.
 !
 !                    &AUXF_AVSO   the filename of output results
 !
 !
 !                  With the input file(s) are prepaired, the program is run on the command line:
 !
 !                  AtomicVirialStress_GPU_Main.exe SETUP-file dev0  ndev
 !                  SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                  dev0        - the ID of the first GPU to be used
 !                  ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2013-10 (Hou Qing, Sichuan university)
 !                  update:     2014-11 (Hou Qing, Sichuan university)
 !

 Program AtomicVirialStress_Tool_main
 use MD_SimBoxArray_ToolShell_14_GPU
 use AtomicVirialStress_GPU
 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_AtomicVirialStress_DEV,   &
                                RECORDPROC=Record_AtomicVirialStress_TOOL,     &
                                AFTRECORD=Clear_AtomicVirialStress_DEV)

       call Main_ANALYSIS(PROCESSID=0) !, FORCETABLE=Register_Interaction_Table, POTTYPE="FS_TYPE")

       stop
 end   program AtomicVirialStress_Tool_main
