  module MD_Globle_Variables
  !***  DESCRIPTION: this module is to define the global variables
  !                  ______________________________________________________
  !                  HOU Qing, Mar, 2010
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MiniUtilities

  implicit none

  contains

  !****************************************************************************
  subroutine ExtractExeName_Globle_Variables( )
  !***  PURPOSE:   to get the excutable name
  !     INPUT
  !
  !     OUTPUT
  !---
      implicit none
      !--- local variables
      integer I, flag
      character(len=4)::EXT

            !$$--- to determine the log file name
             call get_command_argument(0,gm_ExeName)
             flag = 0
             do I=1, len_trim(gm_ExeName)
                !iachar('\') == 92
                !iachar('/') == 47
                if(iachar(gm_ExeName(I:I)) .eq. 92 .or. &
                   iachar(gm_ExeName(I:I)) .eq. 47) flag = I
             enddo

             gm_ExeName = gm_ExeName(flag+1:len_trim(gm_ExeName))
             gm_ExeName = adjustl(gm_ExeName)
             I = len_trim(gm_ExeName)

             EXT(1:4) = gm_ExeName(I-3:I)
             call UpCase(EXT)
             if(EXT(1:4) .eq. ".EXE") then
                gm_ExeName = gm_ExeName(1:I-4)
             end if
             return
  end subroutine ExtractExeName_Globle_Variables
  !****************************************************************************

  !****************************************************************************************
  subroutine ExtractCmdlineFileNames_Globle_Variables()
  !***  DESCRIPTION: to get the command line information for the main process running in "single file" mode
  implicit none

  !--- Local variable
      character*256::cmdline=""
      character*32::tstr(5)
      integer::I, NB

           call ExtractExeName_Globle_Variables()

          !$$--- to find out the input files
           call get_command(cmdline)
           gm_cfgFileName = ""
           call extract_optstr(cmdline, "-","I",  size(gm_cfgFileName), NB,gm_cfgFileName)
           if(NB .le. 0) call extract_optstr(cmdline, "-","IF", size(gm_cfgFileName), NB,gm_cfgFileName)
           if(NB .le. 0) then
              write(*,*) "MDPSCU Error: No configuration file available"
              write(*,*) "              Usage: "//gm_ExeName(1:len_trim(gm_ExeName))//" -I filename"
              write(*,*) "              Process to be stopped"
              stop
           end if

           !$$--- to find out the output files
           gm_outFileName = ""
           call extract_optstr(cmdline, "-","O",  size(gm_outFileName), NB,gm_outFileName)
           if(NB .le. 0) call extract_optstr(cmdline, "-","OF", size(gm_outFileName), NB, gm_outFileName)
           if(NB.le.0) then
              do I=1, size(gm_outFileName)
                 write(tstr(1),*) I
                 tstr = adjustl(tstr(1))
                 gm_outFileName(I) = gm_ExeName(1:len_trim(gm_ExeName))//"_out"//tstr(1)(1:len_trim(tstr(1)))
              end do
           else
              do I=1, NB
                 if(len_trim(gm_outFileName(I)) .gt. 0) then
                    call CreateDataFolder_Globle_Variables(gm_outFileName(I))
                 end if
              end do
           end if

           !$$--- to find out the control files
           gm_ctrlFileName = ""
           call extract_optstr(cmdline, "-","CF",  size(gm_ctrlFileName), NB, gm_ctrlFileName)
           if(NB .le. 0) call extract_optstr(cmdline, "-","CTRF", size(gm_ctrlFileName), NB, gm_ctrlFileName)

    return
  end subroutine ExtractCmdlineFileNames_Globle_Variables
  !****************************************************************************************

  !****************************************************************************************
  subroutine CommandlineBoxsel(CtrlParam)
  !***  DESCRIPTION: to extract JOBID, CFGID and BOXID from a command line
  !                  these parameter in CtrlParam will be overwrited
  !
  !
   implicit none
   !---dummy vaiables
       type(SimMDCtrl), target::CtrlParam

   !--- ;ocal variables
       character*512::cmdline
       character*8::tstr(3)=""
       integer::I, NP, LT, LC
       type(SimMDCtrl), pointer::nextp=>null(), next=>null()

           call get_command(cmdline)
           I = index(cmdline, gm_ExeName(1:len_trim(gm_ExeName)))
           I = I+len_trim(gm_ExeName)
           cmdline = cmdline(I:len_trim(cmdline))

           LT = 0
           LC = 0
           !$$--- the configuration to be determine by the setup file
           call extract_optstr(cmdline, "-","T", 3, NP, tstr, ind=LT)
           if(NP .le. 0) call extract_optstr(cmdline, "-", "TEST", 3, NP, tstr, ind=LT)
           if(NP .le. 0) call extract_optstr(cmdline, "-", "J",    3, NP, tstr, ind=LT)
           if(NP .le. 0) call extract_optstr(cmdline, "-", "Job",  3, NP, tstr, ind=LT)
           if(NP .ge. 1) CtrlParam%JOBID0    = ISTR(tstr(1))
           if(NP .ge. 2) CtrlParam%JOBID1    = ISTR(tstr(2))
           if(NP .ge. 3) CtrlParam%JOBIDSTEP = ISTR(tstr(3))
           if(CtrlParam%JOBID1 .le. 0) CtrlParam%JOBID1 = CtrlParam%JOBID0


           call extract_optstr(cmdline, "-", "C", 3, NP,tstr, ind=LC)
           if(NP .le. 0) call extract_optstr(cmdline, "-","CFG", 3, NP, tstr, ind=LC)
           if(NP .ge. 1) CtrlParam%STARTCFG = ISTR(tstr(1))
           if(NP .ge. 2) CtrlParam%ENDCFG   = ISTR(tstr(2))
           if(NP .ge. 3) CtrlParam%CFGSTEP  = ISTR(tstr(3))
           if(CtrlParam%ENDCFG .le. 0) CtrlParam%ENDCFG = CtrlParam%STARTCFG

           call extract_optstr(cmdline, "-","B",    3, NP,tstr)
           if(NP .le. 0) call extract_optstr(cmdline, "-","BOX", 3, NP,tstr)
           if(NP .ge. 1) CtrlParam%STARTBOX = ISTR(tstr(1))
           if(NP .ge. 2) CtrlParam%ENDBOX   = ISTR(tstr(2))
           if(NP .ge. 3) CtrlParam%BOXSTEP  = ISTR(tstr(3))

           if(LT.gt.0 .and. LC.gt.0) then
              if(LT .gt. LC) then
                 CtrlParam%TIMELOOPOUT = 1
              else
                 CtrlParam%TIMELOOPOUT = 0
              end if
           end if

           !**** if not such parameter input, created the default values
           if(CtrlParam%JOBID0 .lt. 0 .and. CtrlParam%JOBID1 .lt. 0) then
              CtrlParam%JOBID0 = 1
              CtrlParam%JOBIDSTEP = 1
              CtrlParam%JOBID1 = CtrlParam%TOTALBOX/CtrlParam%MULTIBOX
           end if

           if(CtrlParam%STARTBOX .lt. 0 .and. CtrlParam%ENDBOX .lt. 0) then
              CtrlParam%STARTBOX = 1
              CtrlParam%BOXSTEP  = 1
              CtrlParam%ENDBOX   = CtrlParam%MULTIBOX
           end if


           if(CtrlParam%STARTCFG .lt. 0 .and. CtrlParam%ENDCFG .lt. 0) then
              CtrlParam%STARTCFG = 0
              CtrlParam%CFGSTEP  = 1
              call NumberCfg_SimMDCtrl(CtrlParam, CtrlParam%ENDCFG)
           end if

           nextp => CtrlParam
           do while(.TRUE.)
              call GetNext_SimMDCtrl(nextp, 1, next)
              if(.not. associated(next)) exit
              next%JOBID0    = CtrlParam%JOBID0
              next%JOBID1    = CtrlParam%JOBID1
              next%JOBIDSTEP = CtrlParam%JOBIDSTEP
              next%STARTCFG  = CtrlParam%STARTCFG
              next%ENDCFG    = CtrlParam%ENDCFG
              next%CFGSTEP   = CtrlParam%CFGSTEP
              next%STARTBOX  = CtrlParam%STARTBOX
              next%ENDBOX    = CtrlParam%ENDBOX
              next%BOXSTEP   = CtrlParam%BOXSTEP
              nextp=>next
           end do
           return
     end subroutine CommandlineBoxsel
  !****************************************************************************************


  !****************************************************************************************
  subroutine OpenLogFile_Globle_Variables( )
  !***  PURPOSE:   to create the file handle of log file
  !     INPUT
  !
  !     OUTPUT
  !---
      implicit none
      !--- local variables
      integer I, flag

            if(gm_hFILELOG .gt. 0) return

             call ExtractExeName_Globle_Variables()
             call AvailableIOUnit(gm_hFILELOG)
             open(UNIT=gm_hFILELOG, FILE=gm_ExeName(1:len_trim(gm_ExeName))//".log", POSITION='APPEND')
             return
  end subroutine OpenLogFile_Globle_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine CreateDataFolder_Globle_Variables(path)
  !***  PURPOSE:   to create the file handle of log file
  !     INPUT :    path, the path name
  !                NOTE: the path is indicated by "\" in windows
  !                      "/" linux
  !
  !     OUTPUT
  !---
      implicit none
      character*(*)::path
      !--- local variables
      character*256::swapPath
      integer I, ERR
      logical EX

             swapPath = path(1:len_trim(path))
             if(path(2:2) .eq. ':') then
                inquire(FILE=path(1:2), EXIST=EX)
                if(.not.EX) then
                   write(*,fmt="(A)") "MDPSCU Error: Device "//path(1:2)//" not found"
                   write(*,fmt="(A)") "              check the filename: ",path(1:len_trim(path))
                   write(*,fmt="(A)") "              Process to be stopped"
                   stop
                end if
             end if

             do I=1, len_trim(path)
                ! iachar('\') == 92
                ! iachar('/') == 47
                if(iachar(path(I:I)) .eq. 92 .or. &
                   iachar(path(I:I)) .eq. 47) then
                   ! $$--- the system dependence function
                   #ifdef WIN_FILE
                   swapPath(I:I) = achar(92)
                   #else
                   swapPath(I:I) = achar(47)
                   #endif

                   inquire(FILE=swapPath(1:I-1), EXIST=EX)
                   if(.not. EX) then
                      write(*,fmt="(A,A)") ' MDPSCU Message: create data folder:', swapPath(1:I-1)
                      call SYSTEM("mkdir "//swapPath(1:I-1))
                      !call EXECUTE_COMMAND_LINE
                   end if
                end if
             enddo
             !path = swapPath(1:len_trim(swapPath))
             return
  end subroutine CreateDataFolder_Globle_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_Globle_Variables( SimBox, CtrlParam, next)
  !***  PURPOSE:   to read in contral parameters from a file
  !     INPUT
  !
  !     OUTPUT     SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !                next,      optional, the flag to show if still having next job
  !
  !     SEE ALSO   Initialize_Globle_Variables_old
  !
      use RAND32_MODULE
      use MD_TYPEDEF_PrintList, only:Clear_PrintProcess
  !---
      implicit none
      type(SimMDBox)    ::SimBox
      type(SimMDCtrl)   ::CtrlParam
      integer,optional  ::next

      !--- local variables
      integer::N, hFile,hFileS, hFileC, NOTHER, ERR, first=0, LINE=0, I
      character*256::STR, BFILE='', CFILE='', INFILE, STRTMP(4)
      character*8::SDATE
      character*10::STIME

      character*32::KEYWORD
      equivalence(INFILE, STRTMP(1))
      logical::opened, EX
      save first, hFile, BFILE, CFILE

          call Default_Parameter_SimMDBox(SimBox)
          call Default_Parameter_SimMDCtrl(CtrlParam)
          call Clear_PrintProcess()
          if(first.le.0) then

          if(command_argument_count().lt.1) then
          !$$ the file name is not given by the command line
          !$$ we give a default one
              INFile = "Sample.dat"
          else
              call get_command_argument(1,INFile)
          end if

            !$$--- check the status of the file
            inquire(FILE=INFile, EXIST=EX)
            if(.NOT.EX) then
               write(*,*) "MDPSCU Error:The input file "//INFile(1:len_trim(INFile))//" not found"
               stop 'The process stop'
            end if

            !$$--- find out an I/O unit
            call AvailableIOUnit(hFile)

            !$$--- to begine loading the parameters
            call DATE_AND_TIME (DATE=SDATE, TIME=STIME)
            write(*,fmt="(' !**** Date:', A, ' at time:', A)") SDATE, STIME
            write(*,*) "!**** Load Simulation parameters from:"
            write(*,*) "!**** ",INFile(1:len_trim(INFile))
            write(*,*)
            open(UNIT=hFile, FILE=INFile, STATUS='OLD')

            !$$--- print out log file
            if(gm_hFileLog .gt. 0) then
              write(gm_hFileLog,*)
              write(gm_hFileLog,*)
              write(gm_hFileLog,fmt="(' !*****************************************************')")
              write(gm_hFileLog,fmt="(' !**** Date:', A, ' at time:', A)") SDATE, STIME
              write(gm_hFileLog,fmt="(' !**** Load Simulation parameters from:')")
              write(gm_hFileLog,fmt="(' !**** ', A)")INFile(1:len_trim(INFile))
              write(gm_hFileLog,*)
            end if

          end if

          !************************************************
          !$$*** to find out start flag
          if(first .le. 0) then
             call GetInputStrLine(hFile,STR, LINE, "!", *200)
          else
             call GetInputStrLine(hFile,STR, LINE, "!", *300)
          end if
          STR = adjustl(STR)
          call GetKeyWord("&", STR, KEYWORD)
          call UpCase(KEYWORD)
          select case(KEYWORD(1:LEN_TRIM(KEYWORD)))
                 case default
                    if(present(NEXT)) then
                       NEXT = -FIRST
                    end if
                    write(*,fmt="('MDPSCU Warning: unknown keyword at line in setup file:', BZI6)") LINE
                    write(*,fmt="('                the keyword should be one of: ', A, A, A, A)") "&START         ", "&RESTART         "
                    write(*,fmt="('                                              ', A, A, A, A)") "&START_ART     ", "&RESTART_ART     "
                    write(*,fmt="('                                              ', A, A, A, A)") "&START_BST     ", "&RESTART_BST     "
                    write(*,fmt="('                                              ', A, A, A, A)") "&START_GMD     ", "&RESTART_GMD     "
                    write(*,fmt="('                                              ', A, A, A, A)") "&START_NEB     ", "&RESTART_NEB     "
                    write(*,fmt="('                                              ', A, A, A, A)") "&START_PARREP  ", "&RESTART_PARREP  "

                    write(*,fmt="('                process to be stopped')")
                    close(hFile)
                    return
                 case("&START")
                      CtrlParam%RESTART = 0

                 case("&RESTART")
                      CtrlParam%RESTART = 1

                 case(&
                      "&START_ART",       &
                      "&START_BST",       &
                      "&START_GMD",       &
                      "&START_NEB",       &
                      "&START_PARREP",    &
                      "&START_TAD",       &

                      "&START-ART",       &
                      "&START-BST",       &
                      "&START-GMD",       &
                      "&START-NEB",       &
                      "&START-PARREP",    &
                      "&START-TAD",       &

                      "&START:ART",       &                      
                      "&START:BST",       &                      
                      "&START:GMD",       &
                      "&START:NEB",       &
                      "&START:PARREP",    &
                      "&START:TAD"        &
                       )

                      CtrlParam%RESTART = 0
                      gm_AppType = KEYWORD(LEN_TRIM("&START_")+1:LEN_TRIM(KEYWORD))

                 case(&
                      "&RESTART_ART",     &
                      "&RESTART_BST",     &
                      "&RESTART_GMD",     &
                      "&RESTART_NEB",     &
                      "&RESTART_PARREP",  &
                      "&RESTART_TAD",     &

                      "&RESTART-ART",     &
                      "&RESTART-BST",     &
                      "&RESTART-GMD",     &
                      "&RESTART-NEB",     &
                      "&RESTART-PARREP",  &
                      "&RESTART-TAD",     &

                      "&RESTART:ART",     &                      
                      "&RESTART:BST",     &                      
                      "&RESTART:GMD",     &
                      "&RESTART:NEB",     &
                      "&RESTART:PARREP",  &
                      "&RESTART:TAD"      &
                       )
                      CtrlParam%RESTART = 1
                      gm_AppType = KEYWORD(LEN_TRIM("&RESTART_")+1:LEN_TRIM(KEYWORD))
          end select


          !$$*** to find out process to be pasue or stop on waring messge
          call Extract_Substr(STR,1,n,STRTMP)
          if(n.gt.0) then
             call  UpCase(STRTMP(1))
             STRTMP(1) = adjustl(STRTMP(1))
             if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "PAUSE") then
                gm_OnWarning = CP_WARNING_PAUSE
             else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "STOP") then
                gm_OnWarning = CP_WARNING_STOP
             else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "CONTIN") then
                gm_OnWarning = CP_WARNING_CONTIN
             else
                write(*,*)"MDPSCU Warning: unknown parameter in START line ", STRTMP(1)(1:LEN_TRIM(STRTMP(1)))
                write(*,fmt="('               check control file at line:', BZI6)") LINE
                write(*,fmt="(' Usgae: ', A14, ' or no parameter - pause when warning message issued')") "&START 'PAUSE'"
                write(*,fmt="('        ', A13, ' - stop  when warning message issued')") "&START 'STOP'"
                write(*,fmt="('        ', A15, ' - continuie when warning message issued')") "&START 'CONTIN'"
                call ONWARNING(CP_WARNING_PAUSE)
             end if
          else !by default
             gm_OnWarning = CP_WARNING_PAUSE
          end if

          !$$*** accumulate the number of jobs
           FIRST = FIRST + 1
           write(*,fmt="(' !********************************************************** ', A)")
           write(*,fmt="(' !  Start JOB # ', I4)")  FIRST
           write(*,fmt="(' ! ', A)")
           select case(gm_AppType(1:LEN_TRIM(gm_AppType)))
                 case("ART")
                       write(*,fmt="(' !  ', A)")           "ART based siumlation to be performed"
                       write(gm_hFileLog,fmt="(' !  ', A)") "ART based siumlation to be performed"
                 case("BST")
                       write(*,fmt="(' !  ', A)")           "BST based siumlation to be performed"
                       write(gm_hFileLog,fmt="(' !  ', A)") "BST based siumlation to be performed"                       
                 case("GMD")
                       write(*,fmt="(' !  ', A)")           "GENERIC MD based siumlation to be performed"
                       write(gm_hFileLog,fmt="(' !  ', A)") "GENERIC MD based siumlation to be performed"
                 case("NEB")
                       write(*,fmt="(' !  ', A)")           "NEB based siumlation to be performed"
                       write(gm_hFileLog,fmt="(' !  ', A)") "NEB based siumlation to be performed"
                 case("PARREP")
                       write(*,fmt="(' !  ', A)")           "PARALELL REPLICE based siumlation to be performed"
                       write(gm_hFileLog,fmt="(' !  ', A)") "PARALELL REPLICE based siumlation to be performed"
                 case("TAD")
                       write(*,fmt="(' !  ', A)")           "TAD based siumlation to be performed"
                       write(gm_hFileLog,fmt="(' !  ', A)") "TAD based siumlation to be performed"
           end select
          if(present(NEXT)) NEXT = FIRST

          !************************************************
          !$$--- to locate the key words and fill the filenames
          NOTHER = 0
          do while(.TRUE.)
             call GetInputStrLine(hFile,STR, LINE,  "!", *100)
             STR = adjustl(STR)
             call GetKeyWord("&", STR, KEYWORD)
             call UpCase(KEYWORD)

             if(KEYWORD(1:LEN_TRIM("&START"))   .EQ. "&START" .OR.  &
                KEYWORD(1:LEN_TRIM("&RESTART")) .EQ. "&RESTART" ) then
                  write(*,*) "MDPSCU Error: key word "// KEYWORD(1:LEN_TRIM(KEYWORD))//" repeated before key word &END"
                  write(*,fmt="('               check setup file at line:', BZI6)") LINE
                  stop 'The process stop'

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&BOXF") then
                  call Extract_Substr(STR,1,N,STRTMP)
                  BFILE = STRTMP(1)
                  !$$--- check if file for box set is available
                  inquire(FILE=BFILE, EXIST=EX)
                  if(.NOT.EX) then
                     write(*,*) "MDPSCU Error: file for box set: "//BFile(1:len_trim(BFile))//" not found"
                     write(*,fmt="('               check setup file at line:', BZI6)") LINE
                     stop 'The process stop'
                   end if

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&CTLF") then
                  call Extract_Substr(STR,1,N,STRTMP)
                  CFILE = STRTMP(1)
                  !$$--- check if file for control parameter is available
                  inquire(FILE=CFile, EXIST=EX)
                  if(.NOT.EX) then
                     write(*,*) "MDPSCU Error: file for control parameter "//CFile(1:len_trim(CFile))//" not found"
                     write(*,fmt="('               check setup file at line:', BZI6)") LINE
                     stop 'The process stop'
                  end if

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&END") then
                  exit

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&INIF") then
                  call Extract_Substr(STR,4,N,STRTMP)
                  SimBox%IniConfig = STRTMP(1)
                  do I=2, N
                     call UpCase(STRTMP(I))
                     if(STRTMP(I) .eq. "YES" .or. STRTMP(I) .eq. "VEL" ) then
                        SimBox%IniCfgFmt = ior(SimBox%IniCfgFmt, CP_INPUT_VEL)
                     else if(STRTMP(I) .eq. "M") then
                        SimBox%IniCfgFmt = ibclr(SimBox%IniCfgFmt, CP_INPUT_SBOXBIT)
                     else if(STRTMP(I) .eq. "S") then
                        SimBox%IniCfgFmt = ibset(SimBox%IniCfgFmt, CP_INPUT_SBOXBIT)
                     end if      
                  end do
                  !--- check if use multiple initial file
                  call Extract_Numb(STR,3,N,STRTMP)
                  if(N .GT. 0) then
                     SimBox%IniCfgID = ISTR(STRTMP(1))
                     if(gm_AppType(1:len_trim(gm_AppType)) .eq. "NEB" ) then
                        if(N.ge.2) then
                           CtrlParam%NEB_StartCfg = ISTR(STRTMP(1))
                           CtrlParam%NEB_EndCfg   = ISTR(STRTMP(2))
                        end if
                        if(N.ge.3) then
                           CtrlParam%NEB_CfgStep  = ISTR(STRTMP(3))
                        end if
                     end if
                  end if

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&TGTF") then
                  call Extract_Substr(STR,4,N,STRTMP)
                  SimBox%TgtConfig = STRTMP(1)
                  !--- check if use multiple tgt file
                  call Extract_Numb(STR,1,n,STRTMP)
                  if(N .GT. 0) then
                     SimBox%MultiTgtConfig = ISTR(STRTMP(1))
                  end if

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&TOUT") then
                   call Extract_Substr(STR,1,N,STRTMP)
                   CtrlParam%f_quantity = STRTMP(1)
                   call CreateDataFolder_Globle_Variables(CtrlParam%f_quantity)

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&COUT") then
                  call Extract_Substr(STR,1,N,STRTMP)
                  CtrlParam%f_geometry = STRTMP(1)
                  call CreateDataFolder_Globle_Variables(CtrlParam%f_geometry)

             else if(KEYWORD(1:LEN_TRIM(KEYWORD)) .EQ. "&SWAP") then
                  call Extract_Substr(STR,1,N,STRTMP)
                  CtrlParam%f_configure = STRTMP(1)
                  call CreateDataFolder_Globle_Variables(CtrlParam%f_configure)

             else if(KEYWORD(1:LEN_TRIM("&AUXF")) .EQ. "&AUXF") then
                  call Extract_Substr(STR,1,N,STRTMP)
                  NOTHER = NOTHER + 1
                  if(NOTHER .le. size(CtrlParam%f_others)) then
                     CtrlParam%f_others(NOTHER) = STRTMP(1)
                     call Extract_Substr1(STR, CtrlParam%f_tag(NOTHER))
                  else
                     write(*,fmt="('MDPSCU Warning: number of &AUXF files larger than permitted value ',I3)") size(CtrlParam%f_others)
                     write(*,fmt="('               check control file at line:', BZI6)") LINE
                     call ONWARNING(CP_WARNING_PAUSE)
                  end if
             end if
          end do

          !************************************************
          if(len_trim(BFILE) .le. 0) then
              write(*,*) "MDPSCU Error: box file file is not given"
              write(*,fmt="('               check setup file for keyword &BOXF', BZI6)")
              stop 'The process stop'
          end if

          call AvailableIOUnit(hFileS)
          write(*,*)
          write(*,*) "!**** Load Simulation box setup from:"
          write(*,*) "!**** ",BFile(1:len_trim(BFile))
          write(*,*)
          open(UNIT=hFileS, FILE=BFILE, STATUS='OLD')
          !$$--- also write to log file
          if(gm_hFileLog .gt. 0) then
             write(gm_hFileLog,*)
             write(gm_hFileLog,fmt="(' !**** Load Simulation box setup from:')")
             write(gm_hFileLog,fmt="(' !****',A)")BFile(1:len_trim(BFile))
             write(gm_hFileLog,*)
          end if
          !$$--- to load the parameters concerning simulation box
          call Load_Parameter_SimMDBox(hFileS, SimBox)

          !************************************************
          if(len_trim(CFILE) .le. 0) then
              write(*,*) "MDPSCU Error: control file is not given"
              write(*,fmt="('               check setup file for keyword &CTLF', BZI6)")
              if(gm_hFileLog .gt. 0) then
                write(*,fmt="('MDPSCU Error: control file is not given')")
                write(*,fmt="('               check setup file for keyword &CTLF', BZI6)")
              end if
              stop 'The process stop'
          end if

          call AvailableIOUnit(hFileC)
          write(*,*) "!**** Load control parameters from:"
          write(*,*) "!**** ",CFile(1:len_trim(CFile))
          write(*,*)
          open(UNIT=hFileC, FILE=CFile, STATUS='OLD')
          if(gm_hFileLog .gt. 0) then
             write(gm_hFileLog,fmt="(' !**** Load control parameters from:')")
             write(gm_hFileLog,fmt="(' !****',A)")CFile(1:len_trim(CFile))
             write(gm_hFileLog,*)
          end if

          !-- to load the simulation controlling parameters
          if(len_trim(CtrlParam%f_geometry) .LE. 0  ) then
             write(*,*) "MDPSCU Error: filename for output configuaration is missed"
             if(gm_hFileLog .gt. 0) then
                 write(gm_hFileLog,*) "MDPSCU Error: filename for output configuaration is missed"
                close(gm_hFileLog)
             end if
             stop
          end if

          if(len_trim(CtrlParam%f_quantity) .LE. 0  ) then
             call GetPath(CtrlParam%f_geometry, CtrlParam%f_quantity)
             CtrlParam%f_quantity = CtrlParam%f_quantity(1:len_trim(CtrlParam%f_quantity))//"therm"
          end if

          if(len_trim(CtrlParam%f_configure) .LE. 0  ) then
             call GetPath(CtrlParam%f_geometry, CtrlParam%f_configure)
             CtrlParam%f_configure = CtrlParam%f_configure(1:len_trim(CtrlParam%f_configure))//"swap"
          end if

          call Load_Parameter_SimMDCtrl(hFileC, CtrlParam,SimBox)
          call DRAND32_PUTSEED(CtrlParam%SEED)

          !-- to check if all file name required have been given
          if(len_trim(SimBox%IniConfig)  .le. 0  ) then
             write(*,*) "MDPSCU Error: filename for initial configuration is missed"
             write(*,*) "              add statement in the setup file: &INIF  filename"
             if(gm_hFileLog .gt. 0) then
                write(gm_hFileLog,*) "MDPSCU Error: filename for initial configuration is missed"
                close(gm_hFileLog)
             end if
             stop
          end if


          !--- for NEB type app, we need check target files
          if(gm_AppType(1:len_trim(gm_AppType)) .eq. "NEB") then
             if(CtrlParam%NEB_StartCfg .lt. 0) then !--- in case of none-serial NEB
                if(len_trim(SimBox%TgtConfig)  .le. 0 ) then
                    write(*,*) "MDPSCU Error: filename for target configuration of NEB application is missed"
                    write(*,*) "              add statement in the setup file: &TGTF  filename"
                    if(gm_hFileLog .gt. 0) then
                       write(gm_hFileLog,*)  "MDPSCU Error: filename for target configuration of NEB application is missed"
                       close(gm_hFileLog)
                    end if
                    stop
                end if
                if(SimBox%TgtConfig .eq. SimBox%IniConfig ) then
                    write(*,*) "MDPSCU Error: filename for initial and target configurations of NEB application is same"
                    write(*,*) "              Process to be stopped"
                    if(gm_hFileLog .gt. 0) then
                       write(gm_hFileLog,*)  "MDPSCU Error: filename for initial and target configurations of NEB application is same"
                       close(gm_hFileLog)
                    end if
                    stop
                end if
                !--- in case given ini- and tgt-file, we have one test
                CtrlParam%MultiBox = 1
                CtrlParam%TotalBox = 1
             else
                if(len_trim(SimBox%TgtConfig)  .le. 0 ) then
                   SimBox%TgtConfig = SimBox%IniConfig
                endif
                if(SimBox%TgtConfig .ne. SimBox%IniConfig) then
                   write(*,fmt="(A)") " MDPSCU Error: for serial NEB calculations, target file and initial file should have the same name"
                   write(*,fmt="(A)") "               check setup file for keyword &INIF and &TGTF"

                   if(gm_hFileLog .gt. 0) then
                      write(gm_hFileLog,*)  "MDPSCU Error: for serial NEB calculations, target file and initial file should have the same name"
                      close(gm_hFileLog)
                   end if
                   stop
                end if
             end if
          end if

         !--- close the input file
         close(hFileS)
         close(hFileC)

         !--- print out process information
         call Print_Globle_Variables(6, SimBox, CtrlParam)
         if(gm_hFileLog.gt.0) then
            call Print_Globle_Variables(gm_hFileLog, SimBox, CtrlParam)
            write(gm_hFileLog,*)
         end if

         !--- check if the input parameters are consistent
         call Check_Globle_Variables( SimBox, CtrlParam)
        return

  100   close(hFile)
        write(*,*) 'MDPSCU Error in reading setup file. Please check if &END statement is missed for JOB',FIRST
        if(gm_hFileLog .gt. 0) then
           write(gm_hFileLog,*) 'MDPSCU Error in reading setup file. Please check if &END statement is missed for JOB',FIRST
           close(gm_hFileLog)
        end if
        STOP

  200   close(hFile)
        write(*,*) 'MDPSCU Error in reading setup file. Please check if &START or &RESTART statement is missed'
        if(gm_hFileLog .gt. 0) then
           write(gm_hFileLog,*) 'MDPSCU Error in reading setup file. Please check if &START or &RESTART statement is missed'
           close(gm_hFileLog)
        end if
        stop

  300   if(present(NEXT)) then
            NEXT = -FIRST
        end if
        return
  end subroutine Initialize_Globle_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_Globle_Variables_old( SimBox, CtrlParam, CtrlFile, printout)
  !***  PURPOSE:   to read in contral parameters from a file
  !     INPUT      CtrlFile,  optional, the filename of input parameters
  !
  !     OUTPUT     SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !
  !---
      use RAND32_MODULE
      implicit none
      type(SimMDBox)          ::SimBox
      type(SimMDCtrl)         ::CtrlParam
      character*(*),  optional::CtrlFile
      integer,        optional::printout

      !--- local variables
      integer::I, J, N, hFile,ERR
      character*256::STR, INFILE, STRTMP
      character*32::STRNUMB(30)
      logical::opened, EX, givenname

         !$$--- to get the input file name
         givenname = .false.
         if(present(CtrlFile)) then
            if(len_trim(CtrlFile) .gt. 0) then
               givenname = .true.
            end if
         end if

         if(.not. givenname) then
            !$$---  if the file name is not given
            !$$     we check if we can get it by the command line argument
            if(command_argument_count().lt.1) then
               !$$ the file name is not given by the command line
               !$$ we give a default one
               INFile = "Sample.dat"
            else
              call get_command_argument(1,INFile)
            end if
          else
             !$$ the file name given by the program
             INFile = CtrlFile
          end if

          !$$--- check the status of the file
          inquire(FILE=INFile, EXIST=EX)
          if(.NOT.EX) then
             write(*,*) "MDPSCU Error:The input file "//INFile(1:len_trim(INFile))//" not found"
             stop 'The process stop'
          end if

         !$$--- find out an I/O unit
          do hFile = 10, 99
            inquire(UNIT=hFile, OPENED=opened)
            if(.not.opened) exit
          end do

          !$$--- to begine loading the parameters
           write(*,*) "!**** Load Control parameters from:"
           write(*,*) "!**** ",INFile(1:len_trim(INFile))
           write(*,*)
          open(UNIT=hFile, FILE=INFile, STATUS='OLD')
  10      format(A256)

          !$$--- to load the parameters concerning simulation box
          call Load_Parameter_SimMDBox_OLD(hFile, SimBox)

          !$$-- to load the simulation controlling parameters
          call Load_Parameter_SimMDCtrl_Old(hFile, CtrlParam,SimBox)
          call DRAND32_PUTSEED(CtrlParam%SEED)

  !$$*** close the input file
  200   close(hFile)
         !$$*** print out process information
         if(.not.present(printout)) then
              call Print_Globle_Variables(6, SimBox, CtrlParam)
         else
            if(printout) call Print_Globle_Variables(6, SimBox, CtrlParam)
         end if

         call Check_Globle_Variables( SimBox, CtrlParam)
        return

  100   close(hFile)
        write(*,*) 'Error in read ',CtrlFile
        stop
  end subroutine Initialize_Globle_Variables_old
  !****************************************************************************

  !****************************************************************************
  subroutine Print_Globle_Variables(hFILE, SimBox, CtrlParam)
  !***  PURPOSE:   to print out simulation parameters
  !     INPUT      SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !                hFile, the output unit
  !    OUTPUT
  !
  !
  !---
      implicit none
      integer         ::hFILE
      type(SimMDBox)  ::SimBox
      type(SimMDCtrl) ::CtrlParam
      !-- Local variables

        call Print_Parameter_SimMDBox(hFile, SimBox);
        call Print_Parameter_SimMDCtrl(hFile,  CtrlParam, SimBox)
        return
  end subroutine Print_Globle_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Print1_Globle_Variables(hFILE, SimBox, CtrlParam)
  !***  PURPOSE:   to print out extra simulation parameters
  !     INPUT      SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !                hFile, the output unit
  !     OUTPUT
  !
  !
  !---
      implicit none
      integer         ::hFILE
      type(SimMDBox)  ::SimBox
      type(SimMDCtrl) ::CtrlParam
      !--- Local variables

        write(hFILE,FMT="(' !*************************************************************')")
        if(CtrlParam%RESTART .eq. 0) then
            write(hFILE,FMT="(' !    A NEW SIMULATION TO BE PEFORMED')")
            write(hFILE,*)
       else if(CtrlParam%RESTART .eq. -1) then
            write(hFILE,FMT="(' !    SIMULATION TO BE PEFORMED FROM PREVIOUS STOP-POINT')")
            write(hFILE,*)
       else
            write(hFILE,FMT="(' !    SIMULATION TO BE PEFORMED FROM PREVIOUS STOP-POINT')")
            write(hFILE,*)
       end if

  end subroutine Print1_Globle_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Check_Globle_Variables( SimBox, CtrlParam)
  !***  PURPOSE:   to convert the unit to CGS, simulation box data and
  !                control parameter data.
  !                length: from latice unit to CM
  !                mass:   from atomic unit to G
  !                tiem:   from fecosecond to second
  !
  !     INPUT
  !     OUTPUT     SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !
      implicit none
      !--- Dummy vairiables
      type(SimMDBox)        ::SimBox
      type(SimMDCtrl),target::CtrlParam
      !--- local variables
      integer::I
      type(SimMDCtrl), pointer::sectCtrlParam, nextCtrlParam

        !--- if particle number is self-consistent
         if(sum(SimBox%NA) .NE. SimBox%NPRT) then
            write(*,*)"MDPSCU Error: the sum of particle number not equal to:",SimBox%NPRT
            write(*,*)"       The process to be stopped."
            stop
         else
           !--- we give the default itype for the atoms
           !    the ITYP will be overwrited in initialzing configure
            SimBox%IPA(1) = 1
            do I=2,SimBox%NGROUP+1
               SimBox%IPA(I) = SimBox%IPA(I-1)+SimBox%NA(I-1)
            end do
         end if

         !--- The box
         SimBox%RR   = SimBox%RR*CP_A2CM
         SimBox%ZL   = SimBox%LATT*SimBox%RR

         if(any(SimBox%BOXLOW .gt. 1.D32)) then
            SimBox%BOXLOW = -SimBox%ZL*C_HALF
            SimBox%BOXUP  =  SimBox%ZL*C_HALF
         else
            SimBox%BOXLOW =  SimBox%BOXLOW*SimBox%RR
            SimBox%BOXUP  =  SimBox%BOXLOW + SimBox%ZL
         end if

         !--- The mass of atoms
         SimBox%CM = SimBox%CM*CP_AU2G

         !--- change contrlparameter lis
         sectCtrlParam=>CtrlParam
         do while(.true.)
            !--- Convert the press unit to CGS
            sectCtrlParam%PEX = sectCtrlParam%PEX/CP_CGS2KBAR
            sectCtrlParam%PEXTENSOR = C_ZERO
            sectCtrlParam%PEXTENSOR(1,1) = sectCtrlParam%PEX
            sectCtrlParam%PEXTENSOR(2,2) = sectCtrlParam%PEX
            sectCtrlParam%PEXTENSOR(3,3) = sectCtrlParam%PEX

            !--- Convert the time step to second
            sectCtrlParam%H   = sectCtrlParam%H*CP_FS2S
            sectCtrlParam%HMI = sectCtrlParam%HMI*CP_FS2S
            sectCtrlParam%HMX = sectCtrlParam%HMX*CP_FS2S
            sectCtrlParam%DMX = sectCtrlParam%DMX*CP_A2CM

            !--- The cutoff range
            sectCtrlParam%RU    = sectCtrlParam%RU*SimBox%RR
            sectCtrlParam%NB_RM = sectCtrlParam%NB_RM*CtrlParam%RU
            !CtrlParam%RU = DMIN1(minval(SimBox%ZL*C_HALF),CtrlParam%RU)
            !CtrlParam%RM = DMIN1(minval(SimBox%ZL*C_HALF),CtrlParam%RM)

            !--- Covert NEB parameters
             sectCtrlParam%STRCUT_DRTol = sectCtrlParam%STRCUT_DRTol*SimBox%RR

             !--- go to next section
             call GetNext_SimMDCtrl(sectCtrlParam, 1, nextCtrlParam)
             if(.not. associated(nextCtrlParam)) exit
             sectCtrlParam=>nextCtrlParam

         end do
         return
  end subroutine Check_Globle_Variables
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Parameter_SimMDCtrl(hFile, CtrlParam, SimBox)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  use MD_SimMDCtrl_ART,    only:Load_Parameter_SimMDCtrl_ART=>Load_Parameter_SimMDCtrl
  use MD_SimMDCtrl_BST,    only:Load_Parameter_SimMDCtrl_BST=>Load_Parameter_SimMDCtrl
  use MD_SimMDCtrl_GMD,    only:Load_Parameter_SimMDCtrl_GMD=>Load_Parameter_SimMDCtrl
  use MD_SimMDCtrl_NEB,    only:Load_Parameter_SimMDCtrl_NEB=>Load_Parameter_SimMDCtrl
  use MD_SimMDCtrl_PARREP, only:Load_Parameter_SimMDCtrl_PARREP=>Load_Parameter_SimMDCtrl
  use MD_SimMDCtrl_TAD,    only:Load_Parameter_SimMDCtrl_TAD=>Load_Parameter_SimMDCtrl

  implicit none
     !--- dummy varioables
     integer,        intent(in)::hFile
     type(SimMDCtrl)           ::CtrlParam
     type(SimMDBox)            ::SimBox

     !--- local variables
      integer::LINE
      character*256::STR
      character*32::KEYWORD
      character(len=5), parameter::CTLSTARTFLG = "&CTLF"
     !----
           LINE = 0
           call GetInputStrLine(hFile,STR, LINE, "!", *100)
           STR = adjustl(STR)
           call GetKeyWord("&", STR, KEYWORD)
           call UpCase(KEYWORD)
           if(KEYWORD(1:LEN_TRIM(CTLSTARTFLG)) .NE. CTLSTARTFLG) then
              rewind(hFile)
              goto 200
              !call Load_Parameter_SimMDCtrl_Old(hFile, CtrlParam, SimBox)
           else
              rewind(hFile)
              select case(gm_AppType(1:len_trim(gm_AppType)) )
                     case default
                          goto 200
                          !call Load_Parameter_SimMDCtrl_GMD(hFile, CtrlParam, SimBox)
                     case ("ART")
                            if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"_ART"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"-ART"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//":ART"          ) then
                               call Load_Parameter_SimMDCtrl_ART(hFile, CtrlParam, SimBox)
                            else
                               goto 200
                            end if

                     case ("BST")
                            if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"_BST"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"-BST"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//":BST"          ) then
                               call Load_Parameter_SimMDCtrl_BST(hFile, CtrlParam, SimBox)
                            else
                               goto 200
                            end if

                     case ("GMD")
                            if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG             .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"_GMD"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"-GMD"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//":GMD"        ) then
                               call Load_Parameter_SimMDCtrl_GMD(hFile, CtrlParam, SimBox)
                            else
                              goto 200
                            end if

                     case ("NEB")
                            if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"_NEB"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"-NEB"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//":NEB"          ) then
                               call Load_Parameter_SimMDCtrl_NEB(hFile, CtrlParam, SimBox)
                            else
                               goto 200
                            end if

                     case ("PARREP")
                            if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"_PARREP"  .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"-PARREP"  .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//":PARREP"     ) then
                               call Load_Parameter_SimMDCtrl_PARREP(hFile, CtrlParam, SimBox)
                            else
                              goto 200
                            end if

                     case ("TAD")
                            if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"_TAD"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//"-TAD"     .or. &
                               KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. CTLSTARTFLG//":TAD"       ) then
                               call Load_Parameter_SimMDCtrl_TAD(hFile, CtrlParam, SimBox)
                            else
                              goto 200
                            end if

               end select
           end if
           return
  !---------------------------------------------------------------
  100    write(*,fmt="(A,I5)")" MDPSCU Error in reading control parameters at line ", LINE
         write(*,fmt="(A)")   " The process to be stopped."
         stop

  200    close(hFile)
         write(*,fmt="(A,I5)")" MDPSCU Error:  the title of the controal file is incompatible with APPTYPE"
         write(*,fmt="(A,I5)")"                the APPTYPE is "//gm_AppType(1:len_trim(gm_AppType))
         write(*,fmt="(A,I5)")"                corresponding title should be "//CTLSTARTFLG//"_"//gm_AppType(1:len_trim(gm_AppType))
         write(*,fmt="(A)")   " The process to be stopped."
         stop
  end subroutine Load_Parameter_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Print_Parameter_SimMDCtrl(hFile, CtrlParam, SimBox)
  !***  PURPOSE:   to print out  the control parameters for a application
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
  use MD_TYPEDEF_SimMDBox
  use MD_SimMDCtrl_ART,    only:Print_Parameter_SimMDCtrl_ART=>Print_Parameter_SimMDCtrl
  use MD_SimMDCtrl_BST,    only:Print_Parameter_SimMDCtrl_BST=>Print_Parameter_SimMDCtrl
  use MD_SimMDCtrl_GMD,    only:Print_Parameter_SimMDCtrl_GMD=>Print_Parameter_SimMDCtrl
  use MD_SimMDCtrl_NEB,    only:Print_Parameter_SimMDCtrl_NEB=>Print_Parameter_SimMDCtrl
  use MD_SimMDCtrl_PARREP, only:Print_Parameter_SimMDCtrl_PARREP=>Print_Parameter_SimMDCtrl
  use MD_SimMDCtrl_TAD,    only:Print_Parameter_SimMDCtrl_TAD=>Print_Parameter_SimMDCtrl
  

  implicit none
     !--- dummy varioables
     integer, intent(in):: hFile
     type(SimMDCtrl)    :: CtrlParam
     type(SimMDBox)     :: SimBox

     !--- local variables


               select case(gm_AppType(1:len_trim(gm_AppType)) )
                     case  default
                            call Print_Parameter_SimMDCtrl_GMD(hFile, CtrlParam, SimBox)
                     case ("ART")
                           call Print_Parameter_SimMDCtrl_ART(hFile, CtrlParam, SimBox)
                     case ("BST")
                           call Print_Parameter_SimMDCtrl_BST(hFile, CtrlParam, SimBox)
                     case ("GMD")
                            call Print_Parameter_SimMDCtrl_GMD(hFile, CtrlParam, SimBox)
                     case ("NEB")
                           call Print_Parameter_SimMDCtrl_NEB(hFile, CtrlParam, SimBox)
                     case ("PARREP")
                            call Print_Parameter_SimMDCtrl_PARREP(hFile, CtrlParam, SimBox)
                     case ("TAD")
                            call Print_Parameter_SimMDCtrl_TAD(hFile, CtrlParam, SimBox)
               end select

         return
  end subroutine Print_Parameter_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Print_AnalyPara_SimMDCtrl(hFile, SimBox, CtrlParam)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,              intent(in)::hFile
     type(SimMDCtrl),target          ::CtrlParam
     type(SimMDBox)                  ::SimBox

     !--- local variables
      integer::ncfg, I
     !----

        !$$**** HEADER
          write(hFile,fmt="(' !************ ANALYSIS CONTROL PARAMETERS *******************')")
              write(hFile,fmt="(' !')")

              !*************************
              write(hFile,FMT="(' !    Starting box-group ID..........................: ', 3(BZI7,1x))") CtrlParam%JOBID0
              write(hFile,FMT="(' !    Ending box-group ID............................: ', 3(BZI7,1x))") CtrlParam%JOBID1
              write(hFile,FMT="(' !    Step of box-group ID..... .....................: ', 3(BZI7,1x))") CtrlParam%JOBIDSTEP
              write(hFile,fmt="(' !')")

              !*************************
              write(hFile,FMT="(' !    Starting configure ID..........................: ', 3(BZI7,1x))") CtrlParam%STARTCFG
              write(hFile,FMT="(' !    Ending configure ID............................: ', 3(BZI7,1x))") CtrlParam%ENDCFG
              write(hFile,FMT="(' !    Step of configure ID..... .....................: ', 3(BZI7,1x))") CtrlParam%CFGSTEP
              write(hFile,fmt="(' !')")

              !*************************
              if(CtrlParam%MULTIBOX .gt. 1) then
                 write(hFile,FMT="(' !    Starting box ID..........................: ', 3(BZI7,1x))") CtrlParam%STARTBOX
                 write(hFile,FMT="(' !    Ending box ID............................: ', 3(BZI7,1x))") CtrlParam%ENDBOX
                 write(hFile,FMT="(' !    Step of box ID..... .....................: ', 3(BZI7,1x))") CtrlParam%BOXSTEP
                 write(hFile,fmt="(' !')")
              end if

         return
  end subroutine Print_AnalyPara_SimMDCtrl
  !****************************************************************************

  !*********************************************************************************
  subroutine Load_ExAnalyCtl_SimMDCtrl(Tag, Fname, SimBox, CtrlParam, OutTag, OutFile, TheInput)
  !***  PURPOSE:   to load the external control parameters for analysis routine.
  !
  !                This routine may be used when extened files containing parameters
  !                for analysis routine is provided.
  !                This routine provides an alternative approach of Load_AnalyCtl_SimMDCtrl
  !                of loading analysis control parameters.
  !                SEE also:  Load_AnalyCtl_SimMDCtrl in MD_TypeDef_SimMDCtrl.F90
  !
  !     INPUT:     Tag:   the TAG of control section for an analysis module
  !                Fname: the file name
  !
  !     OUTPUT:    SimBox,   the simulation box, with the member KWDList
  !                          could be changed.
  !                CtrlParam, the control parameters, with the member parameters
  !                           controlling analysis members updated.
  !
  !                OutFile,   optional, the filename for the analysis module to output data
  !
  implicit none
  !----   DUMMY Variables
   character*(*)                         ::Tag, Fname
   type(SimMDBox)                        ::SimBox
   type(SimMDCtrl)                       ::CtrlParam
   character*(*),       optional         ::OutTag, OutFile
   type(StatementList), optional, pointer::TheInput
  !--- local
   integer::hFile, LINE, I, N, STATU
   character*256::STR
   character*32::KEYWORD
   type(InputStatements)::INPUTS

            !$$--- load input statements
            call LoadExInput_SimMDCtrl(Tag, Fname, CtrlParam)
            call ExtractAnalyComm_SimMDCtrl(CtrlParam)

            !$$--- to extract control parameter specific for this module
            call GetExInput_SimMDCtrl(CtrlParam, Tag, INPUTS, STATU)
            if(present(TheInput)) then
               call Copy_StatementList(INPUTS%Stat, TheInput)
            end if

            !$$--- to extract the output path
            if(present(OutFile) .and. present(OutTag)) then
               call Get_InputStatements(OutTag, INPUTS, STR, LINE)
               if(LINE .gt. 0) then
                  call Extract_Substr(STR,1,N,OutFile)
               end if
            end if

            !$$--- to scan the input to find out the needed keywords
            do I=1, Number_StatementList(INPUTS%Stat)
               call Get_StatementList(I, INPUTS%Stat, STR, Line=LINE)
               call GetKeyWord("&", STR, KEYWORD)
               call UpCase(KEYWORD)
               if(KEYWORD(1:len_trim("&PROP_")) .eq. "&PROP_") then
                  !$$--- add the keword to the keyword list required for loading from the file
                  !$$    I do not know the type of the property, I cannot create the datapad.
                  !$$    The creationg of datapad left to user
                   call AddDataProKWD_SimMDBox(Simbox, KEYWORD(LEN_TRIM("&PROP_")+1:LEN_TRIM(KEYWORD))//"COL")
               end if
            end do
            call Release_InputStatements(INPUTS)
            return
    end subroutine Load_ExAnalyCtl_SimMDCtrl
  !*********************************************************************************

  !****************************************************************************
  subroutine AddAtoms_Globle_Variables( SimBox, CtrlParam)
  !***  PURPOSE:   to add atoms to simulation box. This subroutine could be called
  !                in analysis tool when MD simulation deposited or embedment atoms
  !                to the box.
  !
  !     INPUT:     SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT:    SimBox,    the Simulation box
  !                CtrlParam, the control parameters
  !
      implicit none
      !--- Dummy vairiables
      type(SimMDBox)        ::SimBox
      type(SimMDCtrl),target::CtrlParam
      !--- local variables
      integer::NEEDNEWATOMS, NA, ITYP
      integer(2)::NEWATOMS(2)
      equivalence(NEWATOMS(1), NEEDNEWATOMS)

           NEEDNEWATOMS = CtrlParam%NEEDNEWATOMS

           NA           = NEWATOMS(1)
           ITYP         = NEWATOMS(2)

           if(NA .le. 0) return
           call AddAtoms_SimMDBox(Simbox, NA, ITYP, TYPEORDER=1)

         return
  end subroutine AddAtoms_Globle_Variables
  !****************************************************************************
  end module MD_Globle_Variables

