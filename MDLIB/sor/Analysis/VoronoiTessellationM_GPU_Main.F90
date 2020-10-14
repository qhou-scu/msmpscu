 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program uses VoronoiTessellation_13_GPU module to do the Voronoi tessellation of
 !                  with the seeds are atoms created by MD simulations. Using this program one can
 !                  obtain the Voronoi volume and surface area of atoms. One can also optionally
 !                  output the Delaunay vertice and Voronoi vertice of atoms. The Voronoi volume can be
 !                  also used for the calculations of local stress of atoms
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       VoronoiTessellationM_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one lines needs to be add:
 !
 !                    &AUXF_VTI  filename1
 !
 !                  where filename1 is the file that provide some control parameters.Optionally,
 !                  one can also add:

 !                    &AUXF_VTO filename2

 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                    &SAVE_VOL      indicating if the Voronoi volumes of atoms will be calculated
 !                                   and out put =1 for yes, =0 for no.
 !
 !                    &SAVE_VV      indicating if the Voronoi vertice of atoms will be svaed to file(s).
 !                                  = 1 for yes, = 0 for no.
 !
 !                    &SAVE_DV      indicating if the save Delaunay vertice of atoms will be saved to file(s).
 !                                  = 1 for yes, = 0 for no
 !
 !                    &MXFACE     the parameter controling the max permitted number of faces of a Voronoi
 !                                volume. This parameter has effects on the memory needed for in calculation.
 !
 !                    &MXVERT     the parameter controling the max permitted number of vertice of a Voronoi
 !                                volume. This parameter has effects on the memory needed for in calculation.
 !                                The parameter of MXVERT should be larger than MXFACE.
 !
 !
 !                    &AUXF_VTO  filename for outputing results. If the filename is given hare,
 !                               the filename2 (if given in SETUP file) will be replaced.

 !                  With the input file(s) are ready, the program can be run  on the command line:
 !
 !                  VoronoiTessellation_GPU.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2013-08     (Hou Qing, Sichuan university)
 !                  updated on 2014-10-21   (Hou Qing, Sichuan university)
 !
 !

!------------------------- External module ------------------
   module MySelectForVoronoiTessellation
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
  !-----------------------------------------------
    implicit none
    contains

        SUBROUTINE SelectAtoms(SimBox, CtrlParam, MASK)
        !***  PURPOSE:   to create the AMASK witch identify the atoms
        !                that will be included in calculation
        !     INPUT:     SimBox:
        !                CtrlParam:
        !
        !     OUTPUT:    FLAG,
        !
         use MD_Globle_Variables_GPU,  only:dm_NPRT, hm_STATU, hm_ITYP
         use  VoronoiTessellation_GPU, only:mp_AMASK_UNSEL,mp_AMASK_SEL
          implicit none
          type(SimMDBox),            intent(in) ::SimBox
          type(SimMDCtrl),           intent(in) ::CtrlParam
          integer,    dimension(:),  intent(out)::MASK
          !--- local
          integer::I

          do I=1, dm_NPRT
             if(IAND(hm_STATU(I),CP_STATU_OUTOFBOX) .eq. CP_STATU_OUTOFBOX) then
                MASK(I) = mp_AMASK_UNSEL
             else
                MASK(I)    = mp_AMASK_SEL
                if(hm_ITYP(I) .ne. 1) then
                   MASK(I) = mp_AMASK_UNSEL
                end if
             end if
           end do

       END SUBROUTINE SelectAtoms

  end module MySelectForVoronoiTessellation
!****************************************************************

!****************************************************************
 Program VoronoiTessellation_Tool_main
 use MD_SimBoxArray_ToolShell_14_GPU
 use VoronoiTessellation_GPU
 use MySelectForVoronoiTessellation
 implicit none

       call SetAtomSelector(HSelector=SelectAtoms)

       call APPSHELL_AddRecord( PRERECORD=Initialize_VoronoiTessellation_DEV, &
                                RECORDPROC=Record_VoronoiTessellation_TOOL,   &
                                AFTRECORD=Clear_VoronoiTessellation_DEV)

       call Main_ANALYSIS(PROCESSID=0)

       stop
 End program VoronoiTessellation_Tool_main
