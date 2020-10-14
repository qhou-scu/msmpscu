 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract the track of clusters.
 !                  A reference lattice will be loaded in, the sites that contain interstitials and
 !                  projectiles are allocated. Afew kind of sites are defined:
 !
 !                  1. A SIA sites are identified by the Wigner_Seize cells containing more than one than
 !                     one substrate atoms.
 !                  2. A impuritie site is defined as a site containing one or more impurity atoms.
 !                  3. A combined site is defined as a site containing at least one impurity atom and
 !                     one substrate atom.
 !                  4. A vacancy site is defined as a site containing no atom.
 !
 !                  The connectivity between the sites are checked. A site cluster is defined as connected
 !                  if they are cloestly neighboring. Two sites are neighboring when they are Wigner-Seitz
 !                  cells share the same faces.
 !
 !                  This program can be used to depict the trajectories of clusters defined above.
 !
 !                  as well as the projectile produced by cascade collision of energtic projectil in materials.
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       Cfg_Track_Cluster_GPU.F90

 !                  HISTORY____________________________________________________________________________
 !                       first version: March, 2017 by HOU Qing
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one line needs to be add:
 !
 !                    &AUXF_RVVIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_RVVOUT filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                  &REFCFG        indication of the filename of the reference configurqation.
 !                                  if the reference configuration is not assigned explicitly,
 !                                  the configuration at zero time step will be the reference
 !                                  configuration.
 !
 !                                  usage:   &REFCFG filename
 !
 !                                  example: &REFCFG "Myreference.crystal"
 !
 !                  If the filename of the reference configuration is given in the control file described
 !                  above, the third file storing the reference configuration needs to be prepaired. The
 !                  the file could be in the format of "&CFGXYZ", "&BOXCFG16" or "t-x-y-z" file whithout
 !                  header.
 !
 !                  With the input file(s) are ready, the program can be run on command line:
 !
 !                  Cfg_Track_Cluster_GPU.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-11 (Hou Qing, Sichuan university)
 !
 !

 !**********************************************************************************
 Program Cfg_Track_Cluster_GPU_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use Cfg_Track_Cluster_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD =Intilize_Cfg_Track_Cluster,   &
                                BEFONETEST=BefOneTest_Cfg_Track_Cluster, &
                                AFTONETEST=AftOneTest_Cfg_Track_Cluster, &
                                AFTRECORD = Clear_Cfg_Track_Cluster,     &
                                RECORDPROC=Record_Cfg_Track_Cluster)
       call Main_ANALYSIS(0)

       stop
 End program Cfg_Track_Cluster_GPU_Main
