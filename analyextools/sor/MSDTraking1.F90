 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to calculate time-averaged MSD of single trjectories
 !                  and also ensemble averaged MSD of TESTS.
 !
 !                  In this program, the total evolution time is fixed, thus the output is 
 !                  the time-averaged MSD vs lag time. This is different from MSDTracking2.F90.
 !                  MSDTracking2.F90 output a set of time-averaged MSD, for each of them the
 !                  lag time is fixed, but the evolution time is variable (SEE MSDTracking2.F90).
 !                            
 !
 !                  The trajectories to be read from the files, xxx.DisTracking, gnerated by 
 !                  DisTracking0.F90 or DisNeighbore0.F90
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  MSDTraking1.F90.exe -I "filename" -t(est) t1, t2, ts -c(fg) c1, c2, cs 
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2019-04 (Hou Qing Sichuan university)
 !
 !
 
 !**********************************************************************************
 Program MSDTraking_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 implicit none
    type(SimMDCtrl)    ::CtrlParam
    character*256      ::InFname  = "", OutFname="", tstr="" 
    character*256      ::GFile0="", GFile1   = "" 

    integer            ::IB, IB0, IC, IC0, IC1, NIC, ITIME0, hFile, INB, NS, NDATA1, ISS=1, LAGSTEP=1, ISS1, NT0, NT=100
    real(KINDDF)       ::TIME0, VEC(3), DIS0(3), EKIN, EPOT, MSD0
    real(KINDDF), allocatable::DISA(:,:), TIME(:), MSD(:,:), MSDE(:), MSDET(:), EB(:)


      !*************************************
        call ExtractCmdlineFileNames_Globle_Variables()
        call CommandlineBoxsel(CtrlParam)
        InFname = gm_cfgFileName(1)
        call GetFname(InFname, OutFname)
        print *, OutFname

        call AvailableIOUnit(hFile)
        open(UNIT=hFile, FILE=InFname)
        read(hFile, *) tstr

        !--- number of configures, that is also number of time points
        NIC = CtrlParam%ENDCFG - CtrlParam%STARTCFG + 1

        !--- number of boxes
        NT0 = (CtrlParam%JOBID1 -  CtrlParam%JOBID0 + 1)/CtrlParam%JOBIDSTEP 

        !--- number of curves to be output for MSDt vs time
        !    MSD, for MSDT of singel trajectoty for give time interval (NIC*timesteps)
        !    MSDE,  for ensembe-averaged MSD
        !    MSDET, for ensembe-averaged of MSDT   
        NT  = 100
        allocate(DISA(NIC,3), TIME(NIC), MSD(min(NT0, NT), NIC), MSDE(NIC), MSDET(NIC), EB(NIC))
        MSD    = 0.D0

        !***********************************************************************************
        !$$- to calculate the MSD for single trajectories 
        ! 
        NT = 0
        do IB0 = CtrlParam%JOBID0, CtrlParam%JOBID1, CtrlParam%JOBIDSTEP 
           do while(.true.)
              read(hFile, *, END=100) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)
              if(IB .eq. IB0) then
                exit
              end if
            end do   

           if(IB .eq. IB0) then
              NT = NT + 1
              if(NT .gt. size(MSD, dim=1)) then
                 NT = NT -1
                 exit
              end if  
              !--- if IC0 < starting cfg id, skip these line
              do while(IC0 .lt. CtrlParam%STARTCFG) 
                 read(hFile, *) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3), EKIN, EPOT
              end do   
              if(IB .ne. IB0 .or. IC0 .ne. CtrlParam%STARTCFG) then
                 write(*,*) "Error 1: ", IB0, IB, IC0, CtrlParam%STARTCFG
                 stop
              end if

              IC   = 1
              TIME(IC)     = TIME0
              DISA(IC,1:3) = DIS0(1:3)
              do while(IC0 .lt. CtrlParam%ENDCFG)
                 read(hFile, *) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)   
                 if(IB .ne. IB0) then
                    write(*,*) "Error 2: ", IB0, IB, IC0, IC
                    stop
                 end if
                 IC = IC + 1
                 DISA(IC,1:3) = DIS0(1:3)
                 TIME(IC)     = TIME0
              end do          
              print *, "Dis loaded for IB=",IB  
              !--- to calculate MSDT for the present trajectory
              !--- loop ing the time interval
              ISS  = LAGSTEP
              do while(ISS .lt. NIC)             
                 MSD0  = 0.D0
                 IC1   = 1
                 if(IC1+ISS .gt. NIC) exit
                 do while(IC1 + ISS .le. NIC)
                    IC = IC1 + ISS
                    VEC(1:3) = (DISA(IC, 1:3) - DISA(IC1, 1:3))
                    MSD0     = MSD0 + sum(VEC*VEC)    
                    IC1      = IC1 + 1
                 end do  
                 MSD0 = MSD0/dble(IC-ISS)
                 ISS1 = ISS/LAGSTEP
                 MSD(NT, ISS1) = MSD(NT,ISS1) + MSD0
                 ISS = ISS + LAGSTEP
              end do
           end if
           
        end do   
        
        !***********************************************************************************
        !$$- to calculate the ensemble average MSD of all trajectories 
        ! 
        rewind(hFile)
        read(hFile, *) tstr
        MSDE   = 0.D0
        MSDET  = 0.D0
        NDATA1 = 0

        do IB0 = CtrlParam%JOBID0, CtrlParam%JOBID1, 1
           do while(.true.)
              read(hFile, *, END=100) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)
              if(IB .eq. IB0) then
                exit
              end if
           end do   
           print *, "Dis loaded for MSD IB=",IB  

           if(IB .eq. IB0) then
              !--- if IC0 < starting cfg id, skip these line
              do while(IC0 .lt. CtrlParam%STARTCFG) 
                 read(hFile, *) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3), EKIN, EPOT
              end do   
              if(IB .ne. IB0 .or. IC0 .ne. CtrlParam%STARTCFG) then
                 write(*,*) "Error 1: ", IB0, IB, IC0, CtrlParam%STARTCFG
                 stop
              end if

              IC   = 1
              TIME(IC)     = TIME0
              DISA(IC,1:3) = DIS0(1:3)
              do while(IC0 .lt. CtrlParam%ENDCFG)
                 read(hFile, *) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)   
                 if(IB .ne. IB0) then
                    write(*,*) "Error 2: ", IB0, IB, IC0, IC
                    stop
                 end if
                 IC = IC + 1
                 DISA(IC,1:3) = DIS0(1:3)
                 TIME(IC)     = TIME0
              end do         

              !--- to calculate ensemble average of SD
              ISS = LAGSTEP
              do while(ISS .lt. NIC)
                 VEC(1:3) = (DISA(ISS, 1:3) - DISA(1, 1:3))
                 ISS1 = ISS/LAGSTEP
                 MSDE(ISS1) = MSDE(ISS1) + sum(VEC*VEC) 
                 ISS  = ISS + LAGSTEP  
              end do  

              !--- to calculate MSDET 
              !--- looping the time interval
              ISS = LAGSTEP
              do while(ISS .lt. NIC)             
                 MSD0  = 0.D0
                 IC1   = 1
                 if(IC1+ISS .gt. NIC) exit
                 do while(IC1 + ISS .le. NIC)
                    IC = IC1 + ISS
                    VEC(1:3) = (DISA(IC, 1:3) - DISA(IC1, 1:3))
                    MSD0     = MSD0 + sum(VEC*VEC)    
                    IC1      = IC1 + 1
                 end do  
                 MSD0 = MSD0/dble(IC-ISS)
                 ISS1 = ISS/LAGSTEP
                 MSDET(ISS1) = MSDET(ISS1) + MSD0
                 ISS = ISS + LAGSTEP
              end do
              NDATA1 = NDATA1 + 1
           end if
        end do        

       !--- calculate the MSDET averaged on ensembles 
        ISS = LAGSTEP
        do while(ISS .lt. NIC) 
           ISS1 = ISS/LAGSTEP
           MSDE(ISS1)  = MSDE(ISS1)/dble(NDATA1)
           MSDET(ISS1) = MSDET(ISS1)/dble(NDATA1)
           ISS = ISS + LAGSTEP   
        end do 

       !--- calculate the Ergodicity breaking parameter of MSDT 
        rewind(hFile)
        read(hFile, *) tstr
        NDATA1   = 0
        EB = 0.D0
        do IB0 = CtrlParam%JOBID0, CtrlParam%JOBID1, 1
           do while(.true.)
              read(hFile, *, END=100) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)
              if(IB .eq. IB0) then
                exit
              end if
           end do   

           if(IB .eq. IB0) then
              !--- if IC0 < starting cfg id, skip these line
              do while(IC0 .lt. CtrlParam%STARTCFG) 
                 read(hFile, *) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3), EKIN, EPOT
              end do   
              if(IB .ne. IB0 .or. IC0 .ne. CtrlParam%STARTCFG) then
                 write(*,*) "Error 1: ", IB0, IB, IC0, CtrlParam%STARTCFG
                 stop
              end if

              IC   = 1
              TIME(IC)     = TIME0
              DISA(IC,1:3) = DIS0(1:3)
              do while(IC0 .lt. CtrlParam%ENDCFG)
                 read(hFile, *) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)   
                 if(IB .ne. IB0) then
                    write(*,*) "Error 2: ", IB0, IB, IC0, IC
                    stop
                 end if
                 IC = IC + 1
                 DISA(IC,1:3) = DIS0(1:3)
                 TIME(IC)     = TIME0
              end do         

              !--- looping the time interval
              ISS = LAGSTEP
              do while(ISS .lt. NIC)             
                 MSD0  = 0.D0
                 IC1   = 1
                 if(IC1+ISS .gt. NIC) exit
                 do while(IC1 + ISS .le. NIC)
                    IC = IC1 + ISS
                    VEC(1:3) = (DISA(IC, 1:3) - DISA(IC1, 1:3))
                    MSD0     = MSD0 + sum(VEC*VEC)    
                    IC1      = IC1 + 1
                 end do  
                 MSD0 = MSD0/dble(IC-ISS)
                 ISS1 = ISS/LAGSTEP
                 EB(ISS1) = EB(ISS1) + MSD0*MSD0
                 ISS = ISS + LAGSTEP
              end do
              NDATA1 = NDATA1 + 1
           end if
        end do        
        ISS = LAGSTEP
        do while(ISS .lt. NIC) 
           ISS1 = ISS/LAGSTEP
           EB(ISS1) = EB(ISS1)/dble(NDATA1)
           if(MSDET(ISS1) .gt. 0.D0) EB(ISS1) = EB(ISS1)/(MSDET(ISS1)*MSDET(ISS1)) - 1.D0
           ISS = ISS + LAGSTEP   
        end do 

     
100     continue
        close(hFile)        
       !***********
 
       call AvailableIOUnit(hFile)
       open(UNIT=hFile, FILE=OutFname(1:len_trim(OutFname))//"_T1.MSD1")
       ISS = LAGSTEP
       do while(ISS .lt. NIC) 
          ISS1 = ISS/LAGSTEP
           write(hFile, fmt="(I8, 1X, 5000(1PE14.4,1x) )")  ISS, TIME(ISS)-TIME(1), &
                        (TIME(ISS)-TIME(1))/(TIME(NIC)-TIME(1)),                    &
                         MSDE(ISS1), MSDET(ISS1),EB(ISS1), (MSD(IC,ISS1), IC=1, NT) 
           ISS = ISS + LAGSTEP   
       end do 
       close(hFile)  
       deallocate(DISA, TIME, MSD, MSDE, MSDET, EB)     
       stop
 End program MSDTraking_Main

