 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract the distance between atoms with the virtual site where
 !                  the point defect is.
 !                  The program different from DisNeighbore0.F90 is that: in DisNeighbore0.F90 calculates the 
 !                  distance between the point defect and its neigbors.
 !
 !**** NOTE:        The configuration are created by using Trajectories.F90, or NearestNeighborsPTyp.F90
 !**** USAGE:       _______________________________________________________________________________________
 !                  DisNeighbore.F90.exe -I "filename" -T(est) t1, t2, ts -C(fg) c1, c2, cs -P deftype
 !                  OPTIONS
 !                        -T(est) t1, t2:   - the range for the tests to be involved in analysis
 !                        -C(est) c1, c2:   - the range for the configurations to be involved
 !                        -P,               - the type of defect atom
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2019-09 (Hou Qing Sichuan university)
 !
 !
 !**********************************************************************************
 Program DisNeighbore1_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 implicit none
    type(SimMDCtrl) ::CtrlParam
    type(SimMDBox)  ::SimBox, RefSimBox
    character*256   ::InFname  = "", OutFname="", Fname="", cmdline 
    character*256   ::GFile0="", GFile1   = "", RefBoxName(1) 

    integer         ::IB, IC0, IA, IN, K, hFile, I, J, NN0, NS, NP, NNN, SITEID
    logical         ::EX0, EX1
    real(KINDDF)    ::LATT, HBS(3), BS(3), LB(3), UB(3), VECT1(3), XP(3), CENTER(3), MIDIS, DIS0
    real(KINDDF), dimension(:), allocatable::DIS
    !--- NOTE the variable DefTyp subjected to changes
    integer            ::DefTyp=2
    character*8        ::tstr(3)=""
    character*72       ::ofmt

      !****************** the collect the input parameters
        call ExtractCmdlineFileNames_Globle_Variables()
        call get_command(cmdline)
        I = index(cmdline, gm_ExeName(1:len_trim(gm_ExeName)))
        I = I+len_trim(gm_ExeName)
        cmdline = cmdline(I:len_trim(cmdline))
        call extract_optstr(cmdline, "-","P", 1, NP, tstr)
        if(NP .ge. 1) then
           DefTyp = ISTR(tstr(1))
        end if    

        call CommandlineBoxsel(CtrlParam)
        InFname = gm_cfgFileName(1)
        call GetPath(InFname, OutFname)
        call GetFname(InFname, Fname)
        OutFname = OutFname(1:len_trim(OutFname))//Fname(1:len_trim(Fname))

        !******** to load the reference box
        call extract_optstr(cmdline, "-","RB", 1, NP, RefBoxName)
        if(NP .lt. 1) call extract_optstr(cmdline, "-","REF", 1, NP, RefBoxName)
        if(NP .lt. 1) then
           write(*,*) "Error: Reference box is missed"
           stop
        end if 
        write(*,*) "Load reference configuration from ", RefBoxName(1)(1:len_trim(RefBoxName(1)) )   
        call Load_Config_SimMDBox(RefBoxName(1), RefSimBox)

        do IB = CtrlParam%JOBID0, CtrlParam%JOBID1
             !---- LOADING THE RAW DATA
              call STRCATI(GFILE0, InFname, "P", 0, 4)
              call STRCATI(GFILE0, GFILE0, "_", IB,  4)
              inquire(FILE=GFILE0, EXIST=EX0)
              if(.not. EX0) then
                  write(*,*)  'MDPSCU Error: '//GFILE0(1:len_trim(GFILE0))//' not exsit'
                  stop 
               end if
               write(*,*) "Load configuration from ", GFILE0(1:len_trim(GFILE0) )
               call Load_Config_SimMDBox(GFILE0, SimBox)
               if(DefTyp .ne. SimBox%ITYP(1)) then
                  write(*,*)  'MDPSCU Error: wrong configuration. The type of the first atom is not the defect type'
                  stop 
               end if 

              !---- prepair the output
               LATT = SimBox%RR
               LB = SimBox%BOXLOW
               UB = SimBox%BOXUP
               BS = SimBox%ZL
               HBS = 0.5D0*BS
               
               NN0 = 0
               NNN = 0
               do I=1, SimBox%NPRT
                  if(NNN .gt. 0) then
                     if(SimBox%ITYP(I) .eq. DefTyp) exit
                     NNN = NNN + 1
                  else
                     if(SimBox%ITYP(I) .eq. DefTyp) then
                        NN0 = NN0 + 1
                     else  
                        NNN = NNN + 1
                     end if    
                  end if   
               end do
               NS = SimBox%NPRT/(NN0 + NNN)

               !--- to begin out put
               allocate(DIS(NN0+NNN))
               call AvailableIOUnit(hFile)
               GFile1 =  GFILE0(1:len_trim(GFILE0))//".DisNBs1"
               open(UNIT=hFile, FILE=GFile1)
               IA = 0
               do IC0 = 1, NS 
                  !--- get the site ID closest to the defect
                  XP = SimBox%XP(IA+1,:)
                  MIDIS = 1.D64
                  do I=1, RefSimBox%NPRT
                     do K=1,3
                        VECT1(K) = RefSimBox%XP(I,K) - XP(K)
                        if(dabs(VECT1(K)) .GT. HBS(K)) then
                           VECT1(K) = VECT1(K) - DSIGN(BS(K),VECT1(K))
                        end if
                     end do
                     DIS0 = sum(VECT1*VECT1) 
                     if(DIS0 .lt. MIDIS) then
                        SITEID = I
                        MIDIS  = DIS0
                     end if    
                  end do
                  CENTER = RefSimBox%XP(SITEID,:)
               
                  !--- get the distance from the site to the atoms
                  IN = 0
                  do I=IA+1, IA+(NN0+NNN)
                     XP = SimBox%XP(I,:)
                     do K=1,3
                        VECT1(K) = XP(K) - CENTER(K)
                        if(dabs(VECT1(K)) .GT. HBS(K)) then
                           VECT1(K) = VECT1(K) - DSIGN(BS(K),VECT1(K))
                        end if
                     end do
                     IN = IN + 1
                     DIS(IN) = dsqrt(sum(VECT1*VECT1))
                  end do

                  IA = IA + NN0 + NNN
                  write(tstr(1),fmt="(I7)") 3+NN0+NNN
                  ofmt ="(I8,1x,I8,1x,"//tstr(1)(1:len_trim(tstr(1)))//"(1pE14.6,1x))"
                  write(hFile, fmt=ofmt) IC0, SITEID, CENTER(1:3)/LATT, (DIS(I)/LATT, I=1, NN0+NNN)
               end do 
           close(hFile)
           if(allocated(DIS)) deallocate(DIS)
        end do           
       !***********
       
       stop
 End program DisNeighbore1_Main

