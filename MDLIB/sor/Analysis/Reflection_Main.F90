 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract reflection and transmission of projectiles on a surfaces
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_CPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least two lines need to be added:
 !
 !                    &AUXF_REFIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_REFOUT filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                    &PROP_TYPE    indicating the atomic type(s) of projectiles
 !
 !                                  usage:  &PROP_TYPE typ1, type2.
 !                                  example:&PROP_TYPE 2, 4, indicating atoms of type 2 and type 4 are projectiles
 !
 !
 !                    &DEPTDIS      indicate the position of upper and low surfaces.
 !
 !                                  usage: &BOUND  z1, z2, nbin
 !                                  where z1 is the position of the upper surface in lattice unit,
 !                                        z2 is the position of the lower surface in lattice unit,
 !                                        nbin is number of bins for histogram of depth distribution
 !
 !                                  example: &DEPTDIS  z1 = 10, z2 = -10, nbin= 100
 !
 !                    &RSPECT       control parameters for reflection energy spectrum
 !
 !                                  usage: &RSPECT  emin, emax, nbin
 !                                  where emin and emax define the energy range for the energy spectrum.
 !                                         nbin is number of bins for histogram of the energy spectrum.
 !
 !                                  example: &RSPECT  0, 10, 32
 !
 !
 !                    &TSPECT       control parameters for transmission energy spectrum
 !
 !                                  usage: &TSPECT  emin, emax, nbin
 !                                  where emin and emax define the energy range for the energy spectrum.
 !                                         nbin is number of bins for histogram of the energy spectrum.
 !
 !                                  example: &TSPECT  0, 10, 32
 !
 !                    &NUMPROJ     control parameters for permitted number of projectiles in each box.
 !                                 the default value is 1.
 !
 !                                  usage: &NUMPROJ num
 !                                  where  num define the number of projectile in each box.
 !
 !                                  example: &NUMPROJ 2
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
 !                  version 1st 2015-06 (Hou Qing, Sichuan university)
 !
 !

 module Reflection
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use MiniUtilities
 implicit none

         character(len=12),parameter, private::mp_FTAGI="&AUXF_REFIN"
         character(len=13),parameter, private::mp_FTAGO="&AUXF_REFOUT"
         character(len=256)::m_INFILE =""                                ! filename of input control data
         character(len=256)::m_OUTFILE =""                               ! filename of output data

         integer, private::m_ATYP(mp_MXGROUP) = 0
         real(KINDDF), private::m_UPSURF= -10.D0
         real(KINDDF), private::m_LOSURF=  10.D0
         integer::m_REFC(mp_MXGROUP)=0, m_TRANSC(mp_MXGROUP)=0
         real(KINDDF), private::m_REFE(mp_MXGROUP)=0.D0, m_TRANSE(mp_MXGROUP)=0.D0
         integer, private::m_rHBINS = 32, m_tHBINS = 32, m_dHBINS = 100
         real(KINDDF), private::m_rEMIN = 0.D0, m_rEMAX = -100.D0
         real(KINDDF), private::m_tEMIN = 0.D0, m_tEMAX = -100.D0
         integer, dimension(:,:), allocatable, private::m_rSPECT, m_tSPECT, m_DEPTDIS
         !--- the working spaces:
         integer, private::m_nProject = 1                             ! the number of projectile in each box
         integer, dimension(:),   allocatable, private::m_hasCount    ! the flag indicating if the particle has been included in
  contains

!*********************************************************************************
  subroutine MyLoadControlParameters(fname, SimBox, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    SimBox,   the simulation box, with the member propColTitle
  !                          could be changed.
  !                CtrlParam, the control parameters, with the member NEEDDAMP
  !                           could be changed.
  !
  !
  !
  use MD_Globle_Variables,only:CreateDataFolder_Globle_Variables
  implicit none
  !----   DUMMY Variables
   character*(*)::fname
   type(SimMDBox)::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
   integer::hFile, N, IC, LINE, I, NN
   character*256::STR
   character*32::STRNUMB(10), KEYWORD
   real(KINDDF)::UPSURF, LOSURF


            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old', err=200)

            !*** start loading control parametgers specific for this module
              LINE = 0
              IC = 0
              do while(.TRUE.)
                  call GETINPUTSTRLINE(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyword("&", STR, KEYWORD)
                  call UpCase(KEYWORD)

                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case("&JOBSEL")
                             !$$*** To get job range to be analysis
                             call Extract_Numb(STR,3,n,STRNUMB)
                             CtrlParam%JOBID0 = ISTR(STRNUMB(1))
                             if(N .GE. 2) CtrlParam%JOBID1  = ISTR(STRNUMB(2))
                             if(N .GE. 3) CtrlParam%JOBIDSTEP = ISTR(STRNUMB(3))

                         case("&CFGSEL")
                            !$$*** To get cfg range to be analysis
                            call Extract_Numb(STR,3,n,STRNUMB)
                            CtrlParam%STARTCFG = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDCFG  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%CFGSTEP = ISTR(STRNUMB(3))

                         case("&BOXSEL")
                            !$$*** To get box range to be analysis
                            call Extract_Numb(STR,3,n,STRNUMB)
                            CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))

                         case( mp_FTAGO)
                              call Extract_Substr(STR,1,n,m_OUTFILE)

                         case("&PROP_TYPE")
                              call Extract_Numb(STR,SimBox%nGroup,n,STRNUMB)
                              do I=1, N
                                 NN = ISTR(STRNUMB(I))
                                 if(NN .gt.0) m_ATYP(NN) = 1
                              end do
                              call AddDataProKWD_SimMDBox(SimBox, "TYPECOL")

                         case("&RSPECT")
                            !$$*** To get number of bins to output size distribution
                            call Extract_Numb(STR,3,n,STRNUMB)
                            if(n.ge. 1) &
                               m_rEMIN = DRSTR(STRNUMB(1))

                            if(n.ge. 2) &
                               m_rEMAX = DRSTR(STRNUMB(2))

                            if(n.ge. 3) then
                              m_rHBINS = ISTR(STRNUMB(3))
                              if(m_rHBINS .le. 0) m_rHBINS = 32
                            end if

                         case("&TSPECT")
                            !$$*** To get number of bins to output size distribution
                            call Extract_Numb(STR,3,n,STRNUMB)
                            if(n.ge. 1) &
                               m_tEMIN = DRSTR(STRNUMB(1))

                            if(n.ge. 2) &
                               m_tEMAX = DRSTR(STRNUMB(2))

                            if(n.ge. 3) then
                              m_tHBINS = ISTR(STRNUMB(3))
                              if(m_tHBINS .le. 0) m_tHBINS = 32
                            end if

                         case("&DEPTDIS")
                              call Extract_Numb(STR,3,n,STRNUMB)
                              UPSURF = -1.D32
                              LOSURF =  1.D32

                              if(n .ge. 1) &
                                 LOSURF = DRSTR(STRNUMB(1))
                              if(n .ge. 2) &
                                 UPSURF = DRSTR(STRNUMB(2))
                              if(n .ge. 3) &
                                 m_dHBINS = ISTR(STRNUMB(3))

                         case("&NUMPROJ")
                              call Extract_Numb(STR,1,n,STRNUMB)
                              if(n .ge. 1) &
                                 m_nProject = ISTR(STRNUMB(1))

                  end select
              end do
    100     close(hFile)

             !--- unit transformation
             if(UPSURF .gt. -1.D30)  then
               m_UPSURF = UPSURF*SimBox%RR
             end if

             if(LOSURF .lt. 1.D30) then
                m_LOSURF = LOSURF*SimBox%RR
             end if

            !$$--- check input consistent
             if(m_UPSURF .le. m_LOSURF) then
                write(*,fmt="(A)")          " MDPSCU Error: the upper surface is lower than the bottom surface"
                write(*,fmt="(' Process to be stopped')")
                stop
             end if

             if(m_rEMIN .ge. m_rEMAX) then
                write(*,fmt="(A)")          " MDPSCU Error: the energy range for reflection energy spectrum is not defined"
                write(*,fmt="(' Process to be stopped')")
                stop
             end if

             if(m_tEMIN .ge. m_tEMAX) then
                write(*,fmt="(A)")          " MDPSCU Error: the energy range for transmission energy spectrum is not defined"
                write(*,fmt="(' Process to be stopped')")
                stop
             end if

             if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for "//gm_ExeName(1:len_trim(gm_ExeName))// " is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
             else
               call CreateDataFolder_Globle_Variables(m_OUTFILE)
             end if

            return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in Reflection module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop

            return
    end subroutine MyLoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
  !----   DUMMY Variables
   type(SimMDBox)::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, NCFG, NCLS, NREC, IFILE


        !$$--- to findout the I/O unit
         IFILE = 0
         do I=1, SIZE(CtrlParam%f_tag)
            if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                m_OUTFILE = CtrlParam%f_others(I)
            end if
            if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
               IFILE = I
            end if
         end do
         !$$--- check the input files
         if(IFILE .le. 0) then
            write(*,fmt="(' MDPSCU Error: the control file is missed in Reflection module')")
            write(*,*) "                 add the keyword in SETUP file:", mp_FTAGI
            write(*,fmt="(' Process to be stopped')")
            stop
         end if

         m_INFILE = CtrlParam%f_others(IFILE)
         write(*,fmt="(A)") " !**** Loading control data from: "//m_INFILE(1:len_trim(m_INFILE))
         call MyLoadControlParameters(m_INFILE, SimBox, CtrlParam)

         if(all(m_ATYP.le.0) ) then
             write(*,fmt="(A)")   " MDPSCU Error: no atom type is specified for projectiles"
             write(*,fmt="(' Process to be stopped')")
             stop
         end if


         !$$--- allocated the memory storing the information we needed
         CtrlParam%TIMELOOPOUT = 1

         NCLS =  CtrlParam%TOTALBOX*m_nProject
         allocate(m_hasCount(NCLS), m_rSPECT(SimBox%nGroup, m_rHBINS), m_tSPECT(SimBox%nGroup, m_tHBINS), m_DEPTDIS(SimBox%nGroup, m_dHBINS) )
         m_hasCount  = 0
         m_rSPECT    = 0
         m_tSPECT    = 0
         m_REFC      = 0
         m_TRANSC    = 0

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
        integer::I, J, IA, NPRT0, NBOX, ITYP, IBIN, NPRO
        real(KINDDF)::rDE, tDE, DD

           if(Stamp%ITime .lt. 0) then
              return
           end if

           NBOX  = size(SimBox)
           NPRT0 = SimBox(1)%NPRT
           NPRO  = size(m_hasCount)

           !$$--- determine the number of  reflection or transmission projectiles
           rDE = (m_rEMAX - m_rEMIN)/dble(m_rHBINS)
           tDE = (m_tEMAX - m_tEMIN)/dble(m_tHBINS)
           DD  = (m_UPSURF - m_LOSURF)/dble(m_dHBINS)
           IA = 0
           m_DEPTDIS = 0
           do I=1, NBOX
              do J = 1, NPRT0
                 ITYP = SimBox(I)%ITYP(J)
                 if(m_ATYP(ITYP).gt.0) then
                    IA = IA + 1

                    if(IA .gt. NPRO) then
                       write(*,fmt="(A)")   " MDPSCU Error: number of projectiles larger than the permitted value"
                       write(*,fmt="(A)")   "               check the input value for &NUMPROJ in your control file"
                       write(*,fmt="(' Process to be stopped')")
                       stop
                    end if

                    if(m_hasCount(IA) .eq. 0) then

                       if(SimBox(I)%XP(J,3) .gt. m_UPSURF .and. iand(SimBox(I)%STATU(J),CP_STATU_ACTIVE).eq.CP_STATU_ACTIVE ) then
                          m_REFC(ITYP) = m_REFC(ITYP) + 1
                          m_REFE(ITYP) = m_REFE(ITYP) + SimBox(I)%EKIN(J)
                          IBIN = int((SimBox(I)%EKIN(J) - m_rEMIN)/rDE) + 1
                          if(IBIN .gt. m_rHBINS) IBIN = m_rHBINS
                          m_rSPECT(ITYP, IBIN) = m_rSPECT(ITYP, IBIN) + 1
                          m_hasCount(IA) = 1

                       else if(SimBox(I)%XP(J,3) .lt. m_LOSURF .and. iand(SimBox(I)%STATU(J),CP_STATU_ACTIVE).eq.CP_STATU_ACTIVE) then
                          m_TRANSC(ITYP) = m_TRANSC(ITYP) + 1
                          m_TRANSE(ITYP) = m_TRANSE(ITYP) + SimBox(I)%EKIN(J)
                          IBIN = int((SimBox(I)%EKIN(J) - m_tEMIN)/tDE) + 1
                          if(IBIN .gt. m_tHBINS) IBIN = m_tHBINS
                          m_tSPECT(ITYP, IBIN) = m_tSPECT(ITYP, IBIN) + 1
                          m_hasCount(IA) = 1

                       else
                          IBIN = int(dabs(SimBox(I)%XP(J,3) - m_UPSURF)/DD) + 1
                          m_DEPTDIS(ITYP,IBIN) = m_DEPTDIS(ITYP,IBIN) + 1
                       end if

                    end if
                 end if
              end do
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
   implicit none
       type(SimMDBox)::SimBox
       type(SimMDCtrl)::CtrlParam
       !--- local variables
        character*256::GFILE
        integer::hFile, I, J, NR(mp_MXGROUP), NT(mp_MXGROUP), NIN(mp_MXGROUP), NP(mp_MXGROUP), ITYP(mp_MXGROUP), NTYP
        real(KINDDF)::rDE, tDE, DD


            !$$--- total particles of each kind
             NP = 0
             do I=1, SimBox%NGROUP
                NP(I) = m_REFC(I) + m_TRANSC(I) + sum(m_DEPTDIS(I,:))
             end do

            !$$--- prepare the output file for reflection
            call AvailableIOUnit(hFile)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//".sum"
            open(hFile, file = GFILE(1:len_trim(GFILE)), status='unknown')
            print *, "Reflection to be output to: ",GFILE(1:len_trim(GFILE))

             !--- determine output format
            write(hFile, fmt="(A, I8)")           '!--- THE REFLECTION RESULTS CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)")               '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)")               '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)")               '!    '

            write(hFile, fmt="(A, I8)")           '!--- Total number of projectiles: ', sum(NP)
            do I=1, SimBox%NGROUP
               if(m_ATYP(I) .gt. 0) then
                  write(hFile, fmt="(A, I8)")     '!---   For projectile of type:                 ', I
                  write(hFile, fmt="(A, I8)")     '!---       effective number of projectiles=    ', NP(I)
                  write(hFile, fmt="(A, I8)")     '!---       number of reflection projectiles=   ', m_REFC(I)
                  write(hFile, fmt="(A, I8)")     '!---       number of transmission projectiles= ', m_TRANSC(I)
                  write(hFile, fmt="(A, I8)")     '!---       number of retention projectiles=    ', sum(m_DEPTDIS(I,:))
                  write(hFile, fmt="(A, I8)")     '!---       refelction energy (eV)=             ', m_REFE(I)*CP_ERGEV
                  write(hFile, fmt="(A, I8)")     '!---       transmission energy (eV)=           ', m_TRANSE(I)*CP_ERGEV
               end if
            end do


            !$$--- determine the number of  reflectg or transmission projectiles
            rDE = (m_rEMAX - m_rEMIN)/dble(m_rHBINS)
            tDE = (m_tEMAX - m_tEMIN)/dble(m_tHBINS)
            DD  = (m_UPSURF - m_LOSURF)/dble(m_dHBINS)
            ITYP = 0
            J=0
            do I=1, size(m_ATYP)
               if(m_ATYP(I) .gt. 0) then
                  J = J + 1
                  ITYP(J) = I
               end if
            end do
            NTYP = J

           !$$--- reflection
           write(hFile, fmt="(A, I8)")           '!--- Refelction energy spectrum:'
            write(hFile,FMT="(12x,A6, 2x, A14,2x, 10(10(A8,2x)))") "I","Energy(eV)", "Counts"
            do I=1, m_rHBINS
                write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" ) &
                                    I, ((I-1)*rDE+m_rEMIN)*CP_ERGEV, (m_rSPECT(ITYP(J),I), J=1,NTYP)
            end do
            write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" )    &
                                    m_rHBINS+1, m_rEMAX*CP_ERGEV

           !$$--- transimission
            write(hFile, fmt="(A, I8)")           '!--- Transmission energy spectrum:'
            write(hFile,FMT="(12x,A6, 2x, A14,2x, 10(10(A8,2x)))") "I","Energy(eV)", "Counts"
            do I=1, m_tHBINS
                write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" ) &
                                    I, ((I-1)*tDE+m_tEMIN)*CP_ERGEV, (m_tSPECT(ITYP(J),I), J=1,NTYP)
            end do
            write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" )    &
                                    m_tHBINS+1, m_tEMAX*CP_ERGEV

           !$$--- depth
            write(hFile, fmt="(A, I8)")           '!--- Depth distribution:'
            write(hFile,FMT="(12x,A6, 2x, A14,2x, 10(10(A8,2x)))") "I","Depth(LU)", "Counts"
            do I=1, m_dHBINS
                write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" ) &
                                    I, dble(I-1)*DD/SimBox%RR, (m_DEPTDIS(ITYP(J),I), J=1,NTYP)
            end do
            write(hFile, FMT = "(12x, I6, 2x, 1pE14.5,2x, 20(I8,2x))" )    &
                                    m_dHBINS+1, dble(m_dHBINS)*DD/SimBox%RR


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
   implicit none
       type(SimMDBox)::SimBox
       type(SimMDCtrl)::CtrlParam

          if(allocated(m_hasCount)) deallocate(m_hasCount)
          if(allocated(m_rSPECT))   deallocate(m_rSPECT)
          if(allocated(m_tSPECT))   deallocate(m_tSPECT)
          if(allocated(m_DEPTDIS))  deallocate(m_DEPTDIS)

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
       type(SimMDBox)::SimBox
       type(SimMDCtrl)::CtrlParam
          call MyOutput(SimBox, CtrlParam)
          call MyCleaner(SimBox, CtrlParam)
          return
  end subroutine MyAfterRECORD
  !**********************************************************************************
 end module Reflection


 !**********************************************************************************
 Program Reflection_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use Reflection

 implicit none

       call APPSHELL_AddRecord( PRERECORD= MyInitialize, &
                                RECORDPROC=MyRECORD,     &
                                AFTRECORD= MyAfterRECORD)

       call Main_ANALYSIS(0)

       stop
 End program Reflection_Main
