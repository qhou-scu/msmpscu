 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract particle distribution from a boundary. By default, the
 !                  boundary is assumed at position of Z = =0, and the simulation box is symetry to the
 !                  boundary.
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
 !                    &AUXF_PDISIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_PDISOUT filename2
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
 !                    &LTHICK      indicating the thickness of bin layer for the distribution.
 !                                 Because the grainbound is assumed at middle in Z direction of the
 !                                 simulation box.  The number of bins to be calcualted by half size
 !                                 of the simulation box divided by the the layer thickness.
 !
 !                                  usage: &LTHICK  thick
 !                                  where thick the thichness in LU
 !
 !                                  example: &LTHICK  0.224
 !
 !
 !                  With the input file(s) are ready, the program can be run on the command line:
 !
 !                  PartDisFromGB.exe SETUP-file dev0  ndev
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

 module PartDisFromGB
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none

         integer, private::m_processid = 0
         character(len=12),parameter, private::mp_FTAGI="&AUXF_PDISIN"
         character(len=13),parameter, private::mp_FTAGO="&AUXF_PDISOUT"
         character(len=256)::m_INFILE =""                                ! filename of input control data
         character(len=256)::m_OUTFILE =""                               ! filename of output data

         integer, private::m_ATYP(mp_MXGROUP) = 0
         real(KINDDF), private::m_LTHICK =  0.5D0                        ! thickness of a bin for the distribution
         !--- the working spaces:
         integer, dimension(:,:), allocatable, private::m_cumPARTDIS     ! the particle distribution integrated over time
         integer, dimension(:,:), allocatable, private::m_avINCPARTDIS   ! the increment of particle in bins averaged over time
         integer, dimension(:,:), allocatable, private::m_cumINTPARTDIS  ! the increment of particle integrated over bins and over time
         integer, dimension(:,:), allocatable, private::m_PARTDIS0       ! the particle distribution at initial time
         integer, dimension(:,:), allocatable, private::m_prePARTDIS     ! the particle distribution at previous recording time
         integer, private::m_INIT = 0
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
   character*(*)  ::fname
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
   integer::hFile, N, IC, LINE, I, NN
   character*256::STR
   character*32::STRNUMB(10), KEYWORD
   real(KINDDF)::LTHICK

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

                         case("&LTHICK")
                            !$$*** To get the thickness of a bin layer for the distribution
                            call Extract_Numb(STR,3,n,STRNUMB)
                            LTHICK = -1
                            if(n.ge. 1) LTHICK = DRSTR(STRNUMB(1))

                  end select
              end do
    100     close(hFile)

            !$$--- check input consistent
             if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for "//gm_ExeName(1:len_trim(gm_ExeName))// " is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
             else
               call CreateDataFolder_Globle_Variables(m_OUTFILE)
             end if

             !--- unit transformation
             if(LTHICK .gt. 0.D0) then
                m_LTHICK = LTHICK*SimBox%RR
             else
                m_LTHICK = 0.5D0*SimBox%RR
             end if
            return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in PartDisFromGB module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop

            return
    end subroutine MyLoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to load control parameters and allocate memories needed
   !
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, NBIN, IFILE


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
            write(*,fmt="(' MDPSCU Error: the control file is missed in PartDisFromGB module')")
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

         !$$--- NOTE: it is assumed that the box is symmetry to Z=0
         !$$          kind
         NBIN =  int(C_HALF*SimBox%ZL(3)/m_LTHICK +0.00001)
         allocate(m_cumPARTDIS(SimBox%nGroup, NBIN), m_PARTDIS0(SimBox%nGroup, NBIN), m_prePARTDIS(SimBox%nGroup, NBIN), &
                  m_cumINTPARTDIS(SimBox%nGroup, NBIN), m_avINCPARTDIS(SimBox%nGroup, NBIN) )
         m_cumPARTDIS  = 0
         m_PARTDIS0    = 0
         m_prePARTDIS  = 0
         m_avINCPARTDIS= 0
         m_cumINTPARTDIS=0
         m_INIT        = 1

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
       !--- dummy varibales
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
        integer::I, J, NPRT0, NBOX, NGROUP, NBIN, IBIN, DIR=3, ICFG, hFile, ITYP(mp_MXGROUP), NTYP
        integer, dimension(:,:), allocatable::PARTDIS, INTPARTDIS
        real(KINDDF)::RR
        character*256::GFILE


           if(Stamp%ITime .lt. 0) then
              return
           end if

           ICFG   = Stamp%IRec(1)
           NBOX   = size(SimBox)
           NPRT0  = SimBox(1)%NPRT
           NGROUP = SimBox(1)%NGROUP
           NBIN   = size(m_cumPARTDIS, dim=2)
           RR     = SimBox(1)%RR
           !$$--- determine the number of  reflectg or transmission projectiles
            ITYP = 0
            J=0
            do I=1, size(m_ATYP)
               if(m_ATYP(I) .gt. 0) then
                  J = J + 1
                  ITYP(J) = I
               end if
            end do
            NTYP = J

           !$$--- determine the number of  reflection or transmission projectiles
           allocate(PARTDIS(size(m_cumPARTDIS, dim=1),NBIN), INTPARTDIS(size(m_cumPARTDIS, dim=1),NBIN) )
           PARTDIS = 0
           do I=1, NBOX
              do J = 1, NPRT0
                 if(m_ATYP( SimBox(I)%ITYP(J)).gt.0) then
                    IBIN = dabs(SimBox(I)%XP(J,DIR))/m_LTHICK + 1
                    if(IBIN .le. NBIN) then
                       PARTDIS( SimBox(I)%ITYP(J), IBIN) = PARTDIS( SimBox(I)%ITYP(J), IBIN) + 1
                    end if
                 end if

              end do
           end do

           !$$--- save the inital distribution
           if(iand(m_INIT, 2) .eq. 0) then
              m_PARTDIS0   = PARTDIS
              m_prePARTDIS = PARTDIS
              m_INIT = ior(m_INIT, 2)
            end if

           do J=1, size(INTPARTDIS, dim=1)
              do IBIN=1, NBIN
                 INTPARTDIS(J, IBIN) = SUM(PARTDIS(J, 1:IBIN)-m_PARTDIS0(J,1:IBIN))
              end do
           end do
           m_cumPARTDIS   = m_cumPARTDIS    + PARTDIS
           m_cumINTPARTDIS= m_cumINTPARTDIS + INTPARTDIS

           do J=1, size(INTPARTDIS, dim=1)
              do IBIN=1, NBIN
                 m_avINCPARTDIS(J, IBIN) = m_avINCPARTDIS(J, IBIN) + SUM(PARTDIS(J, 1:IBIN)-m_prePARTDIS(J,1:IBIN))
              end do
           end do

           !$$--- output the instant distribution
           call STRCATI(GFILE,  m_OUTFILE, "P", m_processid, 4)
           call STRCATI(GFILE,  GFILE, "_", 1, 4)
           call STRCATI(GFILE, GFILE, ".", ICFG, 4)
           write(*,fmt="(A, 2x, I8)") "Instant particle distribution output to "//GFILE(1:len_trim(GFILE))
            call AvailableIOUnit(hFile)
            open(hFile, file = GFILE(1:len_trim(GFILE)), status='unknown')

            write(hFile, fmt="(A, I8)")             '!--- THE PARTICLE DISTRIBUTION FROM GRAIN BOUNDARY CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)")                 '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)")                 '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)")                 '!    '
            write(hFile, fmt="(A, I8)")             '!---   The number of boxes:             ', NBOX
            do I=1, NGROUP
               if(m_ATYP(I) .gt. 0) then
                 write(hFile, fmt="(A, I8, A, I8)") '!---   The types of particles included: ', I, ', with number ', sum(PARTDIS(I, :))
               end if
            end do

            write(hFile, fmt="(A, I8)")           '!--- Distance distribution:'
            write(hFile,fmt="(A)")                "!---  SYMBOL EXPLAINATION:"
            write(hFile,fmt="(A)")                "!     I:            the # of bins"
            write(hFile,fmt="(A)")                "!     Distance(LU):  the position of the bins from the boundary"
            write(hFile,fmt="(A)")                "!     Distrib:       the cuurenct particle distribution"
            write(hFile,fmt="(A)")                "!     delta.Distrib: the increment of particles in bins related to the initial distribution"
            write(hFile,fmt="(A)")                "!     Int.Inc.Dist:  the integrated increment of particles in bins related to the initial distribution"
            write(hFile,fmt="(A)")                "!     Av.Inc.Dist:   the time and space integration increment of particles"
            write(hFile,fmt="(A)")                "!     Cum.Int.Dist:  the time and space integration increment related to the initial distribution"

            write(hFile,FMT="(12x,A6, 2x, A16, 2x, 100(10(A12,2x)))") "I","Distance(LU)",  &
                                                    ("Distrib", "delta.Dist", "Int.Inc.Dist", "Av.Inc.Dist", "Cum.Int.Dist", J=1,NTYP)
            do I=1, NBIN
                write(hFile, FMT = "(12x, I6, 3x, 1pE14.5,3x, 100(I12,2x))" )      &
                                    I, dble(I-1)*m_LTHICK/RR,                      &
                                   (PARTDIS(ITYP(J),I),                            &
                                    PARTDIS(ITYP(J),I)-m_PARTDIS0(ITYP(J),I),      &
                                    INTPARTDIS(ITYP(J),I),                         &
                                    m_avINCPARTDIS(ITYP(J),I),                     &
                                    m_cumINTPARTDIS(ITYP(J),I),                    &
                                    J=1,NTYP)

            end do

            close(hFile)

           !--- save the current particle distribution
           m_prePARTDIS = PARTDIS
           deallocate(PARTDIS, INTPARTDIS)

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
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
       !--- local variables
        character*256::GFILE
        integer::hFile, I, J, ITYP(mp_MXGROUP), NTYP


            ITYP = 0
            J=0
            do I=1, size(m_ATYP)
               if(m_ATYP(I) .gt. 0) then
                  J = J + 1
                  ITYP(J) = I
               end if
            end do
            NTYP = J

            !$$--- prepare the output file for reflection
            call AvailableIOUnit(hFile)
            GFILE = m_OUTFILE(1:len_trim(m_OUTFILE))//".sum"
            open(hFile, file = GFILE(1:len_trim(GFILE)), status='unknown')
            print *, "Accumlation particle distribtion ",GFILE(1:len_trim(GFILE))

            write(hFile, fmt="(A, I8)")             '!--- THE PARTICLE DISTRIBUTION FROM GRAIN BOUNDARY CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)")                 '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)")                 '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)")                 '!    '
            do I=1, SimBox%NGROUP
               if(m_ATYP(I) .gt. 0) then
                 write(hFile, fmt="(A, I8, A, I8)") '!---   The types of particles included: ', m_ATYP(I), ', with number ', sum(m_cumPARTDIS(I, :))
               end if
            end do

            write(hFile, fmt="(A, I8)")            '!--- Distance distribution:'
            write(hFile,FMT="(A1, 11x,A6, 2x, A16,2x, 100(A13,2x))") "!","I","Distance(LU)", &
                  ( "Cum.part.dist", J=1,NTYP), ( "Av.Inc.Dist", J=1,NTYP), ("Cum.Int.Dist", J=1,NTYP)
            do I=1, size(m_cumPARTDIS, dim=2)
                write(hFile, FMT = "(12x, I6, 3x, 1pE14.5,3x, 100(I13,2x))" ) &
                             I, dble(I-1)*m_LTHICK/SimBox%RR,                 &
                            (m_cumPARTDIS(ITYP(J),I), J=1,NTYP),              &
                            (m_avINCPARTDIS(ITYP(J),I), J=1,NTYP),            &
                            (m_cumINTPARTDIS(ITYP(J),I), J=1,NTYP)

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
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          if(allocated(m_cumPARTDIS))    deallocate(m_cumPARTDIS)
          if(allocated(m_avINCPARTDIS))  deallocate(m_avINCPARTDIS)
          if(allocated(m_cumINTPARTDIS)) deallocate(m_cumINTPARTDIS)
          if(allocated(m_PARTDIS0))      deallocate(m_PARTDIS0)
          if(allocated(m_prePARTDIS))    deallocate(m_prePARTDIS)
          m_INIT = 0

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
 end module PartDisFromGB


 !**********************************************************************************
 Program PartDisFromGB_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use PartDisFromGB

 implicit none

       call APPSHELL_AddRecord( PRERECORD= MyInitialize, &
                                RECORDPROC=MyRECORD,     &
                                AFTRECORD= MyAfterRECORD)

       call Main_ANALYSIS(0)

       stop
 End program PartDisFromGB_Main
