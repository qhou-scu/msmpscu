 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program uses PropClusterCommon_GPU to extract the rate of dissociation of free clusters.
 !                  If a cluster is sepearted into two or more parts, the cluster is identified as dissociated.
 !                  Because the seperated parts may collids and coalence again due toe periodic boundary
 !                  only the first occurence of dissolication is idenfied.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       PropClusterCommon_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least two lines need to be added:
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
 !                  FreeCluster_Dissociation_Main.exe SETUP-file dev0  ndev
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

 module FreeCluster_Dissociation
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none
         integer, private::m_NREC                                        ! number of actually record
         real(KINDDF), dimension(:), allocatable, private::m_cumTime     ! the time
         integer, dimension(:,:),    allocatable, private::m_cumDis      ! the number of occurrence of dissoication at give time
         integer, dimension(:),      allocatable, private::m_1stDis      ! the size of small cluster dissociated from the large cluster for the first time
                                                                         ! occurence of dissociation.
                                                                         ! this variable also plays flag indicating the occurence dissociation of cluster
         integer, parameter, private::mp_NUMACLS = 10                    ! the max permitted size of cluster segement diisociated from a large cluster
  contains
  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to Initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use PropClusterCommon_GPU, only:Initialize_PropClusteringCommon, SetFlagProc, &
                                      m_INFILE, m_ATYP, &
                                      SelectAtomsDefault=>SelectAtoms

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::NCFG, NCLS, NREC

         call SetFlagProc(SelectAtomsDefault)

         !$$--- loading control parameters
         call Initialize_PropClusteringCommon(SimBox, CtrlParam)

         if(all(m_ATYP.le.0) ) then
             write(*,fmt="(A)")   " MDPSCU Error: no atom type is specified"
             write(*,fmt="(A)")   "               process to be stoped"
             stop
         end if

         !$$--- allocated the memory storing the information we needed
         NCLS =  CtrlParam%TOTALBOX
         call NumberCfg_SimMDCtrl(CtrlParam, NCFG)
         call NumberRec_SimMDCtrl(CtrlParam, NREC)

         allocate(m_cumTime(max(NCFG, NREC)+1), m_1stDis(NCLS), m_cumDis(max(NCFG, NREC)+1, mp_NUMACLS))
         m_cumTime = 0.D0
         m_1stDis  = 0
         m_cumDis  = 0
         m_NREC    = 0
         return
   end subroutine MyInitialize
  !**********************************************************************************

  !****************************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_AppShell_14_GPU.F90
  !            MD_SimBoxArray_ToolShell_14_GPU.F90
  use PropClusterCommon_GPU, only:Record_PropCluster, m_SCS, m_NCS, m_ASQCS
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        integer::I, IC, IP, NPRT0, NBOX, BSHIFT, IB, IC1, IC2, IS, SNREC
        integer::JOB
        real(KINDDF)::TIME

           if(Stamp%ITime .lt. 0) then
              m_NREC = 0
              return
           end if
           JOB = Stamp%ITest
           TIME = Stamp%Time

           call Record_PropCluster(Stamp, SimBox, CtrlParam)
           NBOX  = size(SimBox)
           NPRT0 = SimBox(1)%NPRT
           m_NREC = m_NREC + 1
           SNREC  = size(m_cumDis, dim=1)
           BSHIFT = (JOB-1)*NBOX

           m_cumTime(m_NREC) = TIME
           IP = 0
           IC1 = 1
           do I=1, NBOX
              IP = (I-1)*NPRT0
              do IC = IC1, m_NCS
                 IP = IP + m_SCS(IC)
                 if((IP-1)/NPRT0  + 1 .gt. I) then
                    exit
                 end if
                 IC2 = IC
              end do

              if(m_1stDis(BSHIFT + I) .le. 0) then ! not dissocaited before
                if(IC2 .gt. IC1) then
                   IS = minval(m_SCS(IC1:IC2))
                   if(IS .gt. mp_NUMACLS) IS = mp_NUMACLS
                   m_1stDis(BSHIFT + I) = IS
                   m_cumDis(m_NREC:SNREC, IS)   = m_cumDis(m_NREC:SNREC, IS) + 1
                   if(IC2 - IC1 .gt. 1) then
                      write(*,fmt="(A)") "MDPSCU Warning: more than one dissociation events occur in a time inverval."
                      write(*,fmt="(A)") "                this may influence the accuracy of the number of occurence of dissociation"
                      write(*,fmt="(A)") "                using more intense recoding time steps is suggested"
                      write(*,fmt="(A)") "                by changing the set for &OUTPUT_C in the control file"
                      call ONWARNING(gm_OnWarning)
                   end if
                end if
              end if
              IC1 = IC2 + 1
           end do

          return
  end subroutine MyRECORD
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyOutput(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
  use PropClusterCommon_GPU, only:m_OUTFILE
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
       !--- local variables
        character*256::GFILE
        integer::hFile, I, J

            !$$--- prepare the output file
            call AvailableIOUnit(hFile)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//".sum"
            open(hFile, file = GFILE(1:len_trim(GFILE)), status='unknown')
            print *, "Summarised cluster to be output to: ",GFILE(1:len_trim(GFILE))

             !--- determine output format
            write(hFile, fmt="(A, I8)")           '!--- THE DISSOCIATION RESULTS CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)")               '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)")               '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)")               '!    '

            write(hFile, fmt="(A, I8)")           '!--- Symbol explaination:'
            write(hFile, fmt="(A, I8)")           '!--- S1, S2, ...., are the accumulated number of small clusters dissociated from a large cluster'
            write(hFile, fmt="(A, I8)")           '!--- where the number 1,2,... indicating the size of small clusters dissociated from the large clustwer'

            write(hFile,FMT="(12x,A6, 2x, A14,2x, 10(10(A8,2x)))") "I","TIME(s)",&
                                       "S1","S2","S3","S4", "S5","S6","S7", "S8","S9","S>9"

            do I=1, m_NREC
                write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" ) &
                                    I, m_cumTime(I)*CP_PS2S, (m_cumDis(I,J), J=1,mp_NUMACLS)
            end do
            close(hFile)
          return
  end subroutine MyOutput
  !**********************************************************************************


  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
  use PropClusterCommon_GPU, only:Clear_PropClusterCommon=>Clear
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
          call Clear_PropClusterCommon(SimBox, CtrlParam)
          if(allocated(m_cumTime)) deallocate(m_cumTime)
          if(allocated(m_1stDis))  deallocate(m_1stDis)
          if(allocated(m_cumDis))  deallocate(m_cumDis)
          m_NREC = 0
          return
  end subroutine MyCleaner
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyAfterRECORD(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
          call MyOutput(SimBox, CtrlParam)
          call MyCleaner(SimBox, CtrlParam)
          return
  end subroutine MyAfterRECORD
  !**********************************************************************************


 end module FreeCluster_Dissociation


 !**********************************************************************************
 Program FreeCluster_Dissociation_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use FreeCluster_Dissociation

 implicit none

       call APPSHELL_AddRecord( PRERECORD=MyInitialize, &
                                RECORDPROC=MyRECORD,    &
                                AFTRECORD=MyAfterRECORD)

       call Main_ANALYSIS(0)

       stop
 End program FreeCluster_Dissociation_Main
