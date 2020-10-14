 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to calculate time-averaged MSD calulated by a set of lag time with
 !                  vraiable total evolution time, also ensemble averaged MSD of TESTS will be caculated.
 !
 !                  In this program, the total evolution time is variable.  This is the difference 
 !                  between this prgraom and from MSDTracking1.F90.
 !                  SEE also MSDTracking1.F90.
 !                            
 !
 !                  The trajectories to be read from the files, xxx.DdisTracking, gnerated by 
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
  
      integer            ::IB, IB0, IC, IC0, IC1, NIC, ITIME0, hFile, NS, NDATA1, ISS=1, NLAG, LAGSTEP, ILAG, L0, L1
      real(KINDDF)       ::TIME0, DTIME, VEC(3), DIS0(3), EKIN, EPOT
      real(KINDDF), allocatable::DISA(:,:), TIME(:), TT(:,:), MSD(:,:), MSDE(:), MSD0(:)
  
  
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
          NIC     = CtrlParam%ENDCFG - CtrlParam%STARTCFG + 1
          LAGSTEP = CtrlParam%CFGSTEP
          NLAG    = min(NIC/LAGSTEP, 50)
  
          !--- number of curves to be output for MSDt vs time
          !    MSD, for MSDT of singel trajectoty for give time interval (NIC*timesteps)
          !    MSDE,  for ensembe-averaged MSD
          !    MSDET, for ensembe-averaged of MSDT   
          allocate(DISA(NIC,3), TIME(NIC), TT(NLAG,NIC), MSD(NLAG, NIC), MSDE(NIC), MSD0(0:NIC))
          MSD    = 0.D0
          
          !***********************************************************************************
          !$$- to calculate the ensemble average MSD of all trajectories 
          ! 
          rewind(hFile)
          read(hFile, *) tstr
          MSD    = 0.D0
          MSDE   = 0.D0
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
                IC1 = 1
                do while(IC1 .lt. NIC)
                   VEC(1:3) = (DISA(IC1, 1:3) - DISA(1, 1:3))
                   MSDE(IC1) = MSDE(IC1) + sum(VEC*VEC) 
                   IC1  = IC1 + 1 
                end do  
  
                !--- to calculate MSDET 
                !--- looping the time interval
                DTIME = TIME(2)-TIME(1)
                do ILAG= 1, NLAG
                   ISS  = LAGSTEP*ILAG
                   IC1  = 1
                   MSD0 = 0.D0
                   do while(IC1 + ISS .le. NIC)
                      VEC(1:3)  = (DISA(IC1 + ISS, 1:3) - DISA(IC1, 1:3))
                      MSD0(IC1) = MSD0(IC1-1) + sum(VEC*VEC)
                      IC1      = IC1 + 1
                   end do  

                   IC1  = 1
                   do while(IC1 + ISS .le. NIC)
                      MSD(ILAG,IC1) = MSD(ILAG,IC1) + MSD0(IC1)/dble(IC1)
                      TT(ILAG,IC1)  = (IC1+ISS)*DTIME
                      IC1 = IC1 + 1
                   end do   
                end do
                NDATA1 = NDATA1 + 1 
             end if
             
          end do        
       
  100     continue
          close(hFile)        
         !***********
  
         call AvailableIOUnit(hFile)
         open(UNIT=hFile, FILE=OutFname(1:len_trim(OutFname))//"_T1.MSD2")

         do IC1 =1, NIC
            L1 = 0
            do ILAG=1, NLAG
               ISS = LAGSTEP*ILAG 
               if(ISS+IC1 .gt. NIC) then
                  exit
               end if
               L1 = ILAG
            end do           
            write(hFile, fmt="(I8, 1X, 1PE14.4,1x,1PE14.4, 1x, 5000(1PE14.4,1x) )")  IC1, TIME(IC1)-TIME(1), &
                MSDE(IC1)/dble(NDATA1), ((TT(ILAG,IC1), MSD(ILAG, IC1)/dble(NDATA1)), ILAG=1, L1)
         end do 
         close(hFile) 
         deallocate(DISA, TIME, TT, MSD, MSDE, MSD0)      
         stop
   End program MSDTraking_Main