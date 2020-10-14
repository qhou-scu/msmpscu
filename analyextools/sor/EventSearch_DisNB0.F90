 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract the time histogram from the data generagted by 
 !                  DisNeighbor0.F90
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  EventSearch_DisNB0.F90.exe -I "filename" -T(est) t1, t2, ts 
 !                  OPTIONS
 !                        -T(est) t1, t2:   - the range for the tests to be involved in analysis
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2019-09 (Hou Qing Sichuan university)
 !
 !
 !**********************************************************************************
 Program DisNeighbore0_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 implicit none
    type(SimMDCtrl)    ::CtrlParam
    character*256      ::InFname  = "", OutFname="", Fname="", cmdline 
    character*256      ::GFile0="", GFile1   = "" 
    integer            ::IB, LINE, COUNT, hFile, hFile1, CHANGED
    logical            ::EX0, EX1
    integer, parameter ::mxCount = 128*2 + 2
    integer, parameter ::mxLen   = 20
    character(len=mxCount*mxLen)::STR 
    character(len=mxLen)::STRN(mxCount)
    integer::ITIMEPRE, ITIME, IBIN, NBIN

    real(kinddf)::HisTW = 1
    integer, dimension(:), allocatable::TimeHis 

      !*************************************
      call ExtractCmdlineFileNames_Globle_Variables()
      call CommandlineBoxsel(CtrlParam)
      InFname = gm_cfgFileName(1)
      call GetPath(InFname, OutFname)
      call GetFname(InFname, Fname)
      OutFname = OutFname(1:len_trim(OutFname))//Fname(1:len_trim(Fname))

        NBIN = 0
        do IB = CtrlParam%JOBID0, CtrlParam%JOBID1
             !---- LOADING THE RAW DATA
              call STRCATI(GFILE0, InFname, "P", 0, 4)
              call STRCATI(GFILE0, GFILE0, "_", IB,  4)
              GFILE0 = GFILE0(1:len_trim(GFILE0))//".DisNBs0"
              inquire(FILE=GFILE0, EXIST=EX0)
              if(.not. EX0) then
                  write(*,*)  'MDPSCU Error: '//GFILE0(1:len_trim(GFILE0))//' not exsit'
                  exit 
               end if
               write(*,*) "Load data from ", GFILE0(1:len_trim(GFILE0) )
               call AvailableIOUnit(hFile)
               open(UNIT=hFile, FILE=GFile0)      
               !First to get the number of lines 
               if(IB .eq. CtrlParam%JOBID0) then
                  LINE  = 0     
                  do while(.true.)
                     call GetInputStrLine(hFile, STR, LINE, "!", *100)
                  end do

 100              allocate(TimeHis(LINE/HisTW + 1))  
                  TimeHis = 0  
                  rewind(hFile)
               end if     

               print *, "BOX", IB    
               ITIMEPRE = 0
               LINE  = 0  
               do while(.true.)
                  call GetInputStrLine(hFile, STR, LINE, "!", *200)
                  call Extract_Numb( STR, mxCount,count,STRN)
                  CHANGED = ISTR(STRN(COUNT))
                  if(CHANGED) then
                     ITIME    = LINE - ITIMEPRE
                     IBIN     = ITIME/HisTW+1
                     TimeHis(IBIN) = TimeHis(IBIN) + 1
                     if(IBIN .gt. NBIN) NBIN = IBIN
                     ITIMEPRE = LINE
                     print *, "Chenage at ", IB, LINE, ITIME, IBIN, NBIN
                  end if   
               end do
   200         close(hFile)    
         end do   
         !--- write out the histogram
         call AvailableIOUnit(hFile)
         GFile1 =  OutFname(1:len_trim(OutFname))//".TimeHis"
         open(UNIT=hFile, FILE=GFile1)
         do IBIN=1, NBIN
            write(hFile,fmt="(I7,1x, I8)")  IBIN, TimeHis(IBIN)
         end do   
         close(hFile)
         deallocate(TimeHis)
       stop
 End program DisNeighbore0_Main

