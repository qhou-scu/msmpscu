 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to seperate tracjectories of defects atoms by change the type 
 !                  the defect atoms. The seperated trajectories could be used for visualization.
 !      
 !
 !**** NOTE:        The configuration are the TRACK files created by using Trajectories_GPU.F90.
 !                  SEE: Trajectories_GPU.F90 
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  DisNeighbore.F90.exe -I "filename" -T(est) t1, t2
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
    type(SimMDCtrl)     ::CtrlParam
    type(SimMDBox)      ::SimBox0, SimBoxOut
    type(MDRecordStamp) ::Stamp0

    character*256     ::InFname  = "", OutFname="", Fname="",cmdline 
    character*256     ::GFIle0="", GFile=""

    integer           ::IB,  IA, IP, IC0, I, J, NN0, NN1, NS, NP, LINE
    logical           ::EX0, EX1

    !--- NOTE the variable DefTyp subjected to changes
    integer            ::DEFTYP(mp_MXGROUP) = 0, DEFFLAG(mp_MXGROUP)=0, TYPMX =0 
    character*256      ::STR 
    character*32       ::tstr(2*mp_MXGROUP)=""
    character*72       ::ofmt

      !*************************************
        call ExtractCmdlineFileNames_Globle_Variables()

        call CommandlineBoxsel(CtrlParam)
        InFname = gm_cfgFileName(1)
        call GetPath(InFname, OutFname)
        call GetFname(InFname, Fname)
        OutFname = OutFname(1:len_trim(OutFname))//Fname(1:len_trim(Fname))

        Stamp0%ICfg = -1
        SimBox0%proAutoload = 1
        do IB = CtrlParam%JOBID0, CtrlParam%JOBID1
             !---- LOADING THE RAW DATA
             call STRCATI(GFILE0, InFname, "P", 0, 4)
             call STRCATI(GFILE0, GFILE0, "_", IB,  4)
             inquire(FILE=GFILE0, EXIST=EX0)
             if(.not. EX0) then
                write(*,*)  'MDPSCU Error: '//GFILE0(1:len_trim(GFILE0))//' not exsit'
                exit 
             end if
             write(*,*) "Load configuration from ", GFILE0(1:len_trim(GFILE0) )
             call Load_Config_SimMDBox(GFILE0, SimBox0)

  
             !--- the default number of  defect atom and the nearest neighber
             NN0 = 1
             NN1 = 8
             call Get_StatementList("&NDEF", SimBox0%proKWDStatment, STR, LINE)
             DEFTYP  = 0
             if(LINE .gt. 0) then
                call Extract_Numb(STR, mp_MXGROUP+1,I,TSTR)
                if(I.gt.0) NN0 = IStr(TSTR(1))
                do J=2, I 
                   DEFTYP(J-1) = IStr(TSTR(J))
                end do   
             end if
             DEFFLAG = 0
             do J=1, NN0
               DEFFLAG(DEFTYP(J)) = 1
             end do  

             call Get_StatementList("&NEAREST", SimBox0%proKWDStatment, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR, 1,I,TSTR)
                if(I.gt.0) NN1 = IStr(TSTR(1))
             end if   

             !--- the number of recording
             NS = SimBox0%NPRT/(NN0*(NN1+1))
             
             !--- determine TYPMX
             TYPMX = maxval(SimBox0%ITYP)

             !--- to change the type of defect atoms
             call CopyInformation_SimMDBox(SimBox0, SimBoxOut)
             call Initialize_SimMDBox(NN0*NS, SimBoxOut)
             IA = 0
             do IC0 = 1, NS 
                do I=1, NN0
                   IA = IA + 1
                   call AddAtoms_SimMDBox(SimBoxOut, (/1/),  (/I + TYPMX/))
                   IP = SimBoxOut%NPRT
                   SimBoxOut%XP(IP,:)  = SimBox0%XP(IA,:)
                   SimBoxOut%XP1(IP,:) = SimBox0%XP1(IA,:)
                   SimBoxOut%FP(IP,:)  = SimBox0%FP(IA,:)
                   SimBoxOut%DIS(IP,:) = SimBox0%DIS(IA,:)
                   SimBoxOut%EKIN(IP)  = SimBox0%EKIN(IA)
                   SimBoxOut%EPOT(IP)  = SimBox0%EPOT(IA)                   
                   IA = IA + NN1  
                end do 
             end do 
             call STRCATI(GFILE, OutFname, "_SEP_P", 0, 4)
             call STRCATI(GFILE, GFILE, "_", IB, 4)
             write(*,*) "Changed configuration output to: ", GFILE(1:len_trim(GFILE) )

            call Putout_Instance_Config_SimMDBox(GFILE, SimBoxOut, Stamp0)             
            call Release_SimMDBox(SimBoxOut)
         end do  !--- end loop for IB         
       !***********
       stop
 End program DisNeighbore0_Main

