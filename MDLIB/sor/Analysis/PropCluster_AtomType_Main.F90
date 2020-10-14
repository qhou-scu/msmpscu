 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program uses PropClusterCommon_GPU and PropCluster_AtomType module to identify
 !                  clustering atoms according to their types, potentials, kinetic energies and connections
 !                  in the space.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       PropCluster_AtomType.F90
 !                       PropClusterCommon_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least two lines need to be add:
 !
 !                    &AUXF_PCLSIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_PCLSOUT filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                    &METHOD       indicate which method to be used to determine the connection
 !                                  between atoms.
 !
 !                                  usage: &METHOD  name parameters
 !                                  where name is a string "USEDT" (default) or "USEBD".
 !
 !                                  If name = "USEDT", the connection to be determined by Delaunay tessellation.
 !                                  If a string "SAVE_DV" appears after "USEDT", the intermediate result
 !                                  of Delaunay vertices will be output to a file (this could be time consuming).
 !
 !                                  If name = "USEBD", the bond length will be used to determine the connections.
 !                                  If the parameter is a real number, then the parameter will be used as the
 !                                  bond length in lattice unit.
 !                                  If without any parameter, the bond length will be the force cutoff distance.
 !
 !                                  example1: &METHOD "USEDT"
 !                                  example2: &METHOD "USEBD"  1.0 (LU)
 !
 !                    &SAVE_DV      indicate if the intermediate result of Delaunay vertices will be output to a file.
 !                                  This parameter takes effect only the METHOD is "USEDT" (default)
 !
 !                    &PROP_TYPE    indicating the type(s) of atoms will in used in clustering and the program will
 !                                  checking wether type data of atoms are available in the configuration files.
 !                                  Without this keyword, the checking will not ber performed.
 !
 !                                  usage:  &PROP_TYPE typ1, type2.
 !                                  example:&PROP_TYPE 2, 4, indicating atoms of type 2 and type 4 will be used in clustering.
 !
 !                    &PROP_EPOT    indicating the potential(s) of atoms will in used in clustering and the program will
 !                                  checking wether potential data of atoms are available in the configuration files.
 !                                  Without this keyword, the checking will not be performed.
 !
 !                                  usage: &PROP_EPOT ( E0, E1 ), ( E2, E3 )...
 !                                  where the parameters E0, E1... (in eV) define a number of energy ranges
 !                                  that could be used in clustering.
 !
 !                                  example:&PROP_EPOT ( 0.1, 0.5 ), ( 1.5, 1.8)...
 !
 !                    &PROP_EKIN    indicating the kinetic energy(s) of atoms will in used in clustering and
 !                                  the program will checking wether kinetic energy data of atoms are available
 !                                  in the configuration files.
 !                                  Without this keyword, the checking will not be performed.
 !
 !                                  usage: &PROP_EKIN ( K0, K1 ), ( K2, K3 )...
 !                                  where the parameters K0, K1... (in eV) define a number of energy range
 !                                  that could be used in clustering.
 !
 !                                  example:&PROP_EKIN ( 0.1, 0.5 ), ( 1.5, 1.8)...
 !
 !                    &PROP_XXX    indicating the property XXX of atoms will in used in clustering and
 !                                 the program will checking wether date of property XXX of atoms are
 !                                 available in configuration files.
 !
 !                    &AUXF_PCLSOUT filename for outputing results. If the filename is is given hare,
 !                                  the filename2 (if given in SETUP file) will be replaced.
 !
 !
 !                  With the input file(s) are ready, the program can be run on the command line:
 !
 !                  ProCluster_AtomType.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2014-11 (Hou Qing, Sichuan university)
 !
 !

 Program PropCluster_AtomType_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use PropCluster_AtomType

 implicit none

       call APPSHELL_AddRecord( PRERECORD=Initialize_PropClustering, &
                                RECORDPROC=Record_PropClusterTool)

       call Main_ANALYSIS(0)

       stop
 End program PropCluster_AtomType_Main
