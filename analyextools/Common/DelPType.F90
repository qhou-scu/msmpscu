 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to delete give type of particles from the box
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  ChangePType.F90.exe -I "filename" -t(est) t1, t2, ts -c(fg) c1, c2, cs -p pt
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2019-07 (Hou Qing Sichuan university)
 !
 !

 !**********************************************************************************
 Program DelPType_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 implicit none
    type(SimMDCtrl)    ::CtrlParam
    type(SimMDBox)     ::SimBox0, SimBox
    type(MDRecordStamp)::Stamp0
    character*256      ::InFname  = "", OutFname="", Fname="" 
    character*256      ::GFile0="", GFile   = "", GFile1 

    integer            ::IB, IC, IA, NP, TY0, TYMASK(mp_MXGROUP), C0000, HASPTAB
    logical            ::EX0, EX1
    character*512::cmdline
    character*8::tstr(3)=""
    

      !*************************************
        call ExtractCmdlineFileNames_Globle_Variables()
        call CommandlineBoxsel(CtrlParam)
        call GetPath(gm_cfgFileName(1),InFname)
        call GetFname(gm_cfgFileName(1), Fname)   
        InFname = InFname(1:len_trim(InFname))//Fname(1:len_trim(Fname))
        call GetPath(InFname, OutFname)
        call GetFname(InFname, Fname)        
        OutFname = OutFname(1:len_trim(OutFname))//Fname(1:len_trim(Fname))

      !*****************************************
        call get_command(cmdline)
        IC = index(cmdline, gm_ExeName(1:len_trim(gm_ExeName)))
        IC = IC+len_trim(gm_ExeName)
        cmdline = cmdline(IC:len_trim(cmdline))
        call extract_optstr(cmdline, "-","P", 2, NP, tstr)
        if(NP .lt. 1) then
           write(*,*) "Error: particle type you want to delete are missed"
           stop
        end if    
        TY0 = ISTR(tstr(1))

        HASPTAB = 1
        call extract_optstr(cmdline, "-", "T", 3, NP, tstr)    
        if(NP .le. 0) HASPTAB = 0

        C0000 = 1
        call extract_optstr(cmdline, "-", "C", 3, NP, tstr)    
        if(NP .le. 0) call extract_optstr(cmdline, "-","CFG", 3, NP, tstr)
        if(NP .le. 0) C0000 = 0 
        
        TYMASK      = 1
        TYMASK(TY0) = 0
        do IB = CtrlParam%JOBID0, CtrlParam%JOBID1
           do IC = CtrlParam%STARTCFG, CtrlParam%ENDCFG 
              if(HASPTAB .gt. 0) then
                 call STRCATI(GFILE0, InFname, "P", 0, 4)
                 call STRCATI(GFILE0, GFILE0, "_", IB,  4)
              else
               GFILE0 = InFname   
              end if
                 
              if(C0000 .gt. 0) call STRCATI(GFILE0, GFILE0, ".", IC, 4)
              inquire(FILE=GFILE0, EXIST=EX0)
              if(.not. EX0) then
                 write(*,*)  'MDPSCU Warning: '//GFILE0(1:len_trim(GFILE0))//' not exsit'
                 exit 
               end if

               write(*,*) "Load configuration from ", GFILE0(1:len_trim(GFILE0) )
               SimBox0%proAutoload = 1
               call Load_Config_SimMDBox(GFILE0, SimBox0, Stamp0)
               call Copy_SimMDBox(SimBox0, SimBox, TYMASK)
               
               !if(OutFname .eq. InFname) then
                  if(HASPTAB .ge. 1) then
                    call STRCATI(GFILE, OutFname, "_DELP_P", 0, 4)
                    call STRCATI(GFILE, GFILE, "_", IB, 4)
                  else
                     GFILE = OutFname(1:len_trim(OutFname))//"_DELP"    
                  end if 

                  GFILE1 = ""  
                  if(C0000.gt.0) then
                     Stamp0%ICfg = IC
                     call STRCATI(GFILE1, GFILE, ".", IC, 4)
                  else 
                     Stamp0%ICfg = -1
                     GFILE1 = GFILE
                  end if  
                  write(*,*) "Changed configuration output to: ", GFILE1(1:len_trim(GFILE1) )    
                  call Putout_Instance_Config_SimMDBox(GFILE1, SimBox, Stamp0)
               !end if   
 
           end do           
        end do     
       !***********
        call Release_SimMDBox(SimBox0)
        call Release_SimMDBox(SimBox)
       stop
 End program DelPType_Main

