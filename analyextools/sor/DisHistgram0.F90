 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract histgram of displacement from the file cretaed by 
 !                  DisTracking0.F90. Along with, the ensemble-averaged MSD to be also output
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  DisHistgram0.exe -s "filename" -T(est) t1, t2, ts -C(fg) c1, c2, cs -B(ox) b1, b2, b3
 !                  OPTIONS
 !                        -T(est) t1, t2, ts:   - the range for the tests to be involved in analysis
 !                        -C(est) c1, c2, cs:   - the range for the configurations to be involved
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2019-04 (Hou Qing Sichuan university)
 !
 !
 
 !**********************************************************************************
 Program DisHistgram0_Main
 use MD_SimBoxArray_ToolShell_14_CPU
 implicit none
    type(SimMDCtrl)    ::CtrlParam
    character*256      ::InFname  = "", OutFname="", tstr="" 
    character*256      ::GFile0="", GFile1   = "" 

    integer            ::IB, IB0, IC, IC0, IC1, NIC, ITIME0, ITIME1, hFile, INB, NS, NP, NDATA, IS, ISS=1
    real(KINDDF)       ::TIME0, TIME1, DTIME, VEC(3), DIS00(3), DIS0(3), DIS1(3), DIS, EKIN, EPOT, MSD0
    real(KINDDF), allocatable::DISA(:,:), TIME(:), MSD(:)

    real(KINDDF),parameter::BINS = 0.7*(0.5D0/2.D0**0.5)
    integer,     parameter::BINN = 2000
    integer,dimension(BINN)::His

      !*************************************
        call ExtractCmdlineFileNames_Globle_Variables()
        call CommandlineBoxsel(CtrlParam)
        InFname = gm_cfgFileName(1)
        call GetFname(InFname, OutFname)
        print *, OutFname

        call AvailableIOUnit(hFile)
        open(UNIT=hFile, FILE=InFname)
        read(hFile, *) tstr

        NIC = CtrlParam%ENDCFG - CtrlParam%STARTCFG + 1
        NS  = 0
        allocate(DISA(NIC,3), TIME(NIC), MSD(NIC))
        His = 0
        MSD = 0.D0
        IB0 = 0
        NP  = 0 
        do while(IB0 .le. CtrlParam%JOBID1) 
           read(hFile, *, END=100) IB, IC0, ITIME0, TIME0, VEC(1:3), DIS0(1:3)
           if(IB .gt. IB0) then
              IB0 = IB
              NP  = NP + 1
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
              !--- to calculate the histogram
              IC = 1
              NS = 0
              do while(IC+CtrlParam%CFGSTEP .le. size(DISA, dim=1))
                 IC1 = IC + CtrlParam%CFGSTEP 
                 VEC(1:3) = (DISA(IC1, 1:3) - DISA(IC, 1:3))
                 DIS      = dsqrt(sum(VEC*VEC))
                 INB      = int(DIS/BINS)+ 1
                 His(INB) = His(INB) + 1     
                 IC = IC + 1
                 NS = NS + 1
              end do
              DTIME = TIME(1+CtrlParam%CFGSTEP) - TIME(1) 

              !--- to calculate MSD
              IS  = ISS 
              do while(.true.) 
                 NDATA = 0
                 MSD0  = 0.D0
                 IC1   = 1
                 if(IC1+IS .gt. NIC) exit
                 do while(.true.)
                    IC = IC1 + IS
                    if(IC .gt. NIC) exit
                    VEC(1:3) = (DISA(IC, 1:3) - DISA(IC1, 1:3))
                    MSD0     = MSD0 + sum(VEC*VEC)    
                    NDATA    = NDATA + 1
                    IC1      = IC1 + 1
                    !MSD(IS)  = MSD(IS) + sum(VEC*VEC)    
                    !IC = IC + 1
                 end do  
                 MSD0 = MSD0/dble(NDATA)
                 !print *, IB, "MSD0", MSD0, IS, IS/ISS, IC, IC1
                 MSD(IS/ISS) = MSD(IS/ISS) + MSD0
                 IS = IS + ISS
              end do

           end if 
           cycle
        end do    
     
100     continue
        close(hFile)        
       !***********
       call AvailableIOUnit(hFile)
       open(UNIT=hFile, FILE=OutFname(1:len_trim(OutFname))//".DisHis")
       write(hFile, *) "! DISHis from "//InFname(1:len_trim(InFname))
       write(hFile, *) "! For CFG=", CtrlParam%STARTCFG, CtrlParam%ENDCFG, CtrlParam%CFGSTEP 
       write(hFile, *) "! Time interval=", DTIME
       do IB = 1, size(His)
         write(hFile, fmt="(I8,1x,1PE14.4,1x,1PE14.4,1x,1PE14.4)")  IB, (IB-1)*BINS, dble(His(IB))/dble(NS), dble(His(IB))/dble(NS)/CP_FOURPI/(IB*BINS)**2
       end do 
       close(hFile)

       call AvailableIOUnit(hFile)
       open(UNIT=hFile, FILE=OutFname(1:len_trim(OutFname))//"_A.MSD")
       do IC = 1, size(MSD)
         if(MSD(IC) .gt. 0) then
         write(hFile, fmt="(I8, 1X, I8,1x,1PE14.4,1x,1PE14.4)")  ISS, IC, TIME(IC*ISS+1)-TIME(1), MSD(IC)/dble(NP)
       end if
       end do 
       close(hFile)       
       stop
 End program DisHistgram0_Main

