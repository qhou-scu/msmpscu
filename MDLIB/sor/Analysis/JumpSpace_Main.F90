 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract jump rate of a give type of atoms. A 'jump' of an atom is
 !                  identified when the displacement of the atom is larger than a threthhold value given
 !                  by a user. A jump is identified, the displacement is reset to zero, and the displacement
 !                  is to be calculated in later time steps.
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
 !                  JumpRate_Main.exe SETUP-file dev0  ndev
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

 module JumpSpace_module
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none

         integer, private::m_processid = 0
         character(len=12),parameter, private::mp_FTAGI="&AUXF_JUMPSI"
         character(len=12),parameter, private::mp_FTAGO="&AUXF_JUMPSO"
         character(len=256)::m_INFILE =""                                ! filename of input control data
         character(len=256)::m_OUTFILE =""                               ! filename of output data

         integer, private::m_ATYP(mp_MXGROUP) = 0
         integer, private::m_ITYP(mp_MXGROUP) = 0
         integer, private::m_NTYP             = 0

         real(KINDDF), private::m_DR          =  0.5D0                    ! the size of bins to generate the histogram
         real(KINDDF), private::m_RCUT        =  10.D0                    ! the largest displacement to generate the histogram
         integer, private::m_TimestpR                                     ! number of configure for each record

         !--- the working spaces:
         real(KINDDF), dimension(:,:), allocatable, private::m_POS0     ! the initial positionm in a jump
         real(KINDDF), dimension(:,:), allocatable, private::m_pPOS       ! the previous position for a time step
         real(KINDDF), dimension(:),   allocatable, private::m_DISPL      ! the dispalcement length values for a jump
         real(KINDDF), dimension(:),   allocatable, private::m_PATH       ! the pathe length betwwn two jump
         integer, dimension(:,:),      allocatable::m_cumPARTDIS, m_cumPARTPATH

         type(SimMDBox)::m_SimBoxSwap
         type(SimMDBox)::m_cumSimBoxSwap
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
   real(KINDDF)::DR, RCUT, JSTEP

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

                         case("&DELTADIS")
                            !$$*** To get the thickness of a bin layer for the distribution
                            call Extract_Numb(STR,1,n,STRNUMB)
                            DR = -1
                            if(n.ge. 1) DR = DRSTR(STRNUMB(1))

                         case("&MAXDIS")
                            !$$*** To get the thickness of a bin layer for the distribution
                            call Extract_Numb(STR,1,n,STRNUMB)
                            RCUT = -1
                            if(n.ge. 1) RCUT = DRSTR(STRNUMB(1))

                         case("&RECSTEPS")
                            !$$*** To get the thickness of a bin layer for the distribution
                            call Extract_Numb(STR,1,n,STRNUMB)
                            m_TimestpR = -1
                            if(n.ge. 1) m_TimestpR = DRSTR(STRNUMB(1))

                  end select
              end do
    100     close(hFile)

            !--- unit transformation
             if(DR .gt. 0.D0) then
               m_DR   = DR*SimBox%RR
             else
               m_DR   = 0.5D0*SimBox%RR
             end if

             if(RCUT .gt. 0.D0) then
               m_RCUT = RCUT*SimBox%RR
             else
              m_RCUT = 5.D0*SimBox%RR
             end if

             if(m_TimestpR .le. 0) then
               m_TimestpR = 1
             end if

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

             call AddDataProKWD_SimMDBox(SimBox, "TYPECOL")
             call AddDataProKWD_SimMDBox(SimBox, "DISPCOL")
            return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in PartDisFromGB module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop
    300     continue

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


         CtrlParam%TIMELOOPOUT = 1
        !$$--- determine the number of  reflectg or transmission projectiles
         m_ITYP = 0
         m_NTYP=0
         do I=1, size(m_ATYP)
            if(m_ATYP(I) .gt. 0) then
               m_NTYP = m_NTYP + 1
               m_ITYP(m_NTYP) = I
            end if
         end do

         m_INIT        = 1

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRecord_FOR_TIMESEG(Stamp, SimBox, CtrlParam)
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
        integer::I, J, IP, NPRT0, NBOX, NGROUP, NBIN, IBIN, ICFG, hFile, ITIME
        integer, dimension(:,:), allocatable::PARTDIS, PARTPATH
        real(KINDDF)::RR, DISP
        character*256::GFILE

           ITIME = Stamp%ITime
           ICFG  = Stamp%IRec(1)

           NBOX   = size(SimBox)
           NPRT0  = SimBox(1)%NPRT
           NGROUP = SimBox(1)%NGROUP
           NBIN =  int(m_RCUT/m_DR +0.00001)
           RR     = SimBox(1)%RR

           !$$--- determine the number of  reflection or transmission projectiles
           allocate(PARTDIS(NGROUP,NBIN), PARTPATH(NGROUP,NBIN) )
           PARTDIS = 0
           PARTPATH = 0
           IP = 0
           do I=1, NBOX
              do J = 1, NPRT0
                 if(m_ATYP( SimBox(I)%ITYP(J)).gt.0) then
                    IP = IP + 1
                    DISP = m_DISPL(IP)

                    IBIN = DISP/m_DR + 1
                    if(IBIN .gt. NBIN) IBIN = NBIN
                    PARTDIS( SimBox(I)%ITYP(J), IBIN) = PARTDIS( SimBox(I)%ITYP(J), IBIN) + 1

                    DISP = m_PATH(IP)
                    IBIN = DISP/m_DR + 1
                    if(IBIN .gt. NBIN) IBIN = NBIN
                    PARTPATH( SimBox(I)%ITYP(J), IBIN) = PARTPATH( SimBox(I)%ITYP(J), IBIN) + 1
                 end if

              end do
           end do
           m_cumPARTDIS  = m_cumPARTDIS   + PARTDIS
           m_cumPARTPATH  = m_cumPARTPATH + PARTPATH

           !$$--- output the instant displacement distribution
           call STRCATI(GFILE,  m_OUTFILE, "P", m_processid, 4)
           call STRCATI(GFILE,  GFILE, "_", 1, 4)
           call STRCATI(GFILE, GFILE, ".", ICFG, 4)
           write(*,fmt="(A, 2x, I8)") "Instant displacement distribution output to "//GFILE(1:len_trim(GFILE))
            call AvailableIOUnit(hFile)
            open(hFile, file = GFILE(1:len_trim(GFILE)), status='unknown')

            write(hFile, fmt="(A, I8)")             '!--- THE DISPLACEMENT DISTRIBUTION CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)")                 '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)")                 '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)")                 '!    '
            write(hFile, fmt="(A, I8)")             '!---   The number of boxes:             ', NBOX
            do I=1, NGROUP
               if(m_ATYP(I) .gt. 0) then
                 write(hFile, fmt="(A, I8, A, I8)") '!---   The types of particles included: ', I, ', with number ', sum(PARTDIS(I, :))
               end if
            end do

            write(hFile, fmt="(A, I8)")           '!--- Displacement distribution:'
            write(hFile,fmt="(A)")                "!---  SYMBOL EXPLAINATION:"
            write(hFile,fmt="(A)")                "!     I:            the # of bins"
            write(hFile,fmt="(A)")                "!     Displace(LU):  the displacement values of the bins"
            write(hFile,fmt="(A)")                "!     Disp.Dis:      the currenct displacement distribution"
            write(hFile,fmt="(A)")                "!     Path.Dis:      the currenct path length distribution"
            write(hFile,fmt="(A)")                "!     Disp.Dis:      the currenct displacement distribution"
            write(hFile,fmt="(A)")                "!     Disp.Cum:      the accumulated displacemednt distribution"
            write(hFile,fmt="(A)")                "!     Disp.Cum:      the accumulated path-length distribution"

            write(hFile,FMT="(A1, 11x,A6, 2x, A16, 2x, 100(10(A12,2x)))") "!", "I","Displace(LU)",             &
                                                                         ("Disp.Dis", "Path.Dis", J=1,m_NTYP), &
                                                                         ("Disp.Cum", "Path.Cum", J=1,m_NTYP)
            do I=1, NBIN
                write(hFile, FMT = "(12x, I6, 3x, 1pE14.5,3x, 100(I12,2x))" )                   &
                                    I, dble(I-1)*m_DR/RR,                                       &
                                   (PARTDIS(m_ITYP(J),I), PARTPATH(m_ITYP(J),I),J=1,m_NTYP),    &
                                   (m_cumPARTDIS(m_ITYP(J),I), m_cumPARTPATH(m_ITYP(J),I),J=1,m_NTYP)
            end do
            close(hFile)

           !--- save the current particle distribution
           deallocate(PARTDIS, PARTPATH)


          return
  end subroutine MyRecord_FOR_TIMESEG
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
        integer::NPRT0, NBOX, ICFG, ISEG,I, J, K, IP, NP, NBIN
        real(KINDDF)::DS
        integer, dimension(:),allocatable::JSTAT, IDJUMP
        save NP
        character*256::GFILE
!-----
           if(Stamp%ITime .lt. 0) then
              return
           end if

           ICFG   = Stamp%IRec(1)
           NBOX   = size(SimBox)
           NPRT0  = SimBox(1)%NPRT

           !$$--- to allocate the memories
           if(iand(m_INIT, 2) .eq. 0) then
              NP = 0
              do I=1, NBOX
                 do J = 1, NPRT0
                    if(m_ATYP( SimBox(I)%ITYP(J)).gt.0) then
                       NP = NP + 1
                    end if
                 end do
              end do
              allocate(m_POS0(NP,3), m_pPOS(NP,3), m_DISPL(NP), m_PATH(NP))

              !$$--- to prepare the displacement configure
              call CopyInformation_SimMDBox(SimBox(1), m_SimBoxSwap)
              m_SimBoxSwap%NPRT   =  NP
              m_SimBoxSwap%NGROUP  = 0
              call Initialize_SimMDBox(m_SimBoxSwap,scheme=2)
              m_SimBoxSwap%NGROUP  = SimBox(1)%NGROUP
              m_SimBoxSwap%NA = 0

              m_POS0  = 0.D0
              m_pPOS  = 0.D0
              m_DISPL = 0.D0
              m_PATH  = 0.D0
              IP = 0
              do I=1, NBOX
                 do J = 1, NPRT0
                    if(m_ATYP( SimBox(I)%ITYP(J)).gt.0) then
                       IP = IP + 1
                       m_POS0(IP, 1)  = SimBox(I)%DIS(J,1)
                       m_POS0(IP, 2)  = SimBox(I)%DIS(J,2)
                       m_POS0(IP, 3)  = SimBox(I)%DIS(J,3)
                       m_pPOS(IP,1:3) = SimBox(I)%DIS(J,1:3)
                       m_SimBoxSwap%ITYP(IP)     = SimBox(I)%ITYP(J)
                       m_SimBoxSwap%NA(SimBox(I)%ITYP(J)) = m_SimBoxSwap%NA(SimBox(I)%ITYP(J)) + 1
                       m_SimBoxSwap%XP(IP, 1:3)  = 0.D0
                       m_SimBoxSwap%DIS(IP, 1:3) = 0.D0
                    end if
                 end do
              end do

              NBIN =  int(m_RCUT/m_DR +0.00001)
              allocate(m_cumPARTDIS(SimBox(1)%NGROUP,NBIN), m_cumPARTPATH(SimBox(1)%NGROUP,NBIN) )
              m_cumPARTDIS = 0
              m_cumPARTPATH = 0
              m_INIT = ior(m_INIT, 2)
           end if

           !$$--- to get the current displacement and path length
           IP = 0
           do I=1, NBOX
              do J = 1, NPRT0
                 if(m_ATYP( SimBox(I)%ITYP(J)).gt.0) then
                    IP = IP + 1
                    m_SimBoxSwap%DIS(IP, 1:3) = SimBox(I)%DIS(J,1:3) - m_POS0(IP,1:3)
                    m_DISPL(IP) = sum(m_SimBoxSwap%DIS(IP, 1:3)*m_SimBoxSwap%DIS(IP, 1:3) )
                    m_DISPL(IP) = dsqrt(m_DISPL(IP))

                    DS = sum( (SimBox(I)%DIS(J,1:3) - m_pPOS(IP,1:3))*(SimBox(I)%DIS(J,1:3) - m_pPOS(IP,1:3))  )
                    DS = dsqrt(DS)
                    m_PATH(IP)  = m_PATH(IP) + DS

                    m_SimBoxSwap%XP(IP, 1:3)  = m_SimBoxSwap%DIS(IP, 1:3) - m_POS0(IP,1:3)
                    m_SimBoxSwap%XP1(IP, 1:3) = SimBox(I)%XP1(J,1:3)
                    m_SimBoxSwap%FP(IP, 1:3)  = SimBox(I)%FP(J,1:3)
                    m_SimBoxSwap%STATU(IP)    = SimBox(I)%STATU(J)
                    m_SimBoxSwap%EPOT(IP)     = SimBox(I)%EPOT(J)
                    m_SimBoxSwap%EKIN(IP)     = SimBox(I)%EKIN(J)
                    m_pPOS(IP,1:3) =  SimBox(I)%DIS(J,1:3)
                 end if
              end do
           end do
           call Merge_SimMDBox(m_SimBoxSwap, m_cumSimBoxSwap)
           call MyRecord_FOR_TIMESEG(Stamp, SimBox, CtrlParam)
           ISEG = ICFG - CtrlParam%STARTCFG
           if(mod(ISEG,m_TimestpR).eq.0)  then
               call STRCATI(GFILE,  m_OUTFILE, "_Cfg_P", 1, 4)
               call STRCATI(GFILE,  GFILE, "_", 1, 4)
               call Putout_Instance_Config_SimMDBox(GFILE, m_cumSimBoxSwap, Stamp)
           end if

          return
  end subroutine MyRECORD
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

          call Release_SimMDBox(m_SimBoxSwap)
          call Release_SimMDBox(m_cumSimBoxSwap)

          if(allocated(m_POS0))  deallocate(m_POS0)
          if(allocated(m_pPOS))  deallocate(m_pPOS)
          if(allocated(m_DISPL))  deallocate(m_DISPL)
          if(allocated(m_PATH))   deallocate(m_PATH)
          if(allocated(m_cumPARTDIS))   deallocate(m_cumPARTDIS)
          if(allocated(m_cumPARTPATH))  deallocate(m_cumPARTPATH)

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
          !call MyOutput(SimBox, CtrlParam)
          call MyCleaner(SimBox, CtrlParam)
          return
  end subroutine MyAfterRECORD
  !**********************************************************************************
 end module JumpSpace_module


 !**********************************************************************************
 Program RJumpRate_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 use JumpSpace_module

 implicit none

       call APPSHELL_AddRecord( PRERECORD= MyInitialize, &
                                RECORDPROC=MyRECORD,     &
                                AFTRECORD= MyAfterRECORD)

       call Main_ANALYSIS(0)

       stop
 End program RJumpRate_Main
