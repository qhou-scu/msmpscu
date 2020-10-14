 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to simulate the displacement distribution for an particle walk on 1D lattice.
 !                  The purpose is to test if the distribution following Gaussian distribution.
 !
 !                  The wating time distribution (WTD) of jump follows POWER-LAW distribution.
 !                   
 !                  This is a GPU version, in comparison with the CPU verision Walk-1D-MSD-Test2-Pow.F90.
 !
 !                  Note: the model parameters are programmatically set. Thus, this program is replaced by
 !                        Walk-1D-MSD-Test1.
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                   Walk-1D-MSD-Test2-Pow-GPU.exe 
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2020-05 (Hou Qing, Sichuan university)
 !
 !

 !**********************************************************************************

 !**********************************************************************************
 program Walk_1D_MSD_Test2_Main
    use MSM_MultiGPU_Basic
    use RAND32_MODULE
    use RAND32SEEDLIB_MODULE,only:GetSeed_RAND32SEEDLIB
    use MC_RandWalk_Pow_GPU
    use MiniUtilities
    implicit none
    integer, parameter::NB   = 100000      !the number of boxes in test
    integer::ISEED0=12345, ISEED(2)=(/3555, 5443/)

    real(KINDDF)::R, DT, TAU=1.D0, MSD, DSTEP=1.D0, TSTEP=1.D0, TTIME=1000.D0, CURT, ATT, STIME, MT, MN, MT2, MN2, CURSTEP, POW=3.0, &
                  DBIN, DJBIN, MD, NJMI, NJMX
    real(KINDDF),allocatable::XP(:), PreTime(:), NextTime(:),NJUMP(:), DDISH(:), NJH(:)
    character*256:: fname="Walk_1D_MSD_Test1", OutFname="", OutFname1 
    integer::hFile1,hFile2, IB, N, M, NBIN=100, L1, L2, NJBIN=1000, IBIN
    character*18::SN1, SN2, SN3, SN4, SN5, SNS
    !--- Device variable
     character*32::ARGV
     integer::IDEV
     type(DevVec_DF)   ::dPreTime
     type(DevVec_DF)   ::dNextTime
     type(DevVec_DF)   ::dNJump
     type(DevVec_DF)   ::dXP
     !---

       TSTEP  = 0.01D0 ! (ps)
       TTIME  = 1.D8
       TAU    = 0.1
       POW    = 0.5
       write(SN1,fmt="(I)") nint(TAU*100)
       write(SN2,fmt="(I)") nint(POW*100)
       write(SN5,fmt="(I)") NB
       write(SNS,fmt="(I)") ISEED0

       write(SN3,fmt="(I)") nint(TTIME)
       write(SN4,fmt="(I)") nint(TSTEP*1000)
       SN1 = adjustl(SN1)
       SN2 = adjustl(SN2)
       SN3 = adjustl(SN3)
       SN4 = adjustl(SN4)
       SN5 = adjustl(SN5)
       SNS = adjustl(SNS)
       
       OutFname = "POW_TAU"//SN1(1:len_trim(SN1))//"_POW"//SN2(1:len_trim(SN2))//"_NB"//SN5(1:len_trim(SN5))//"_TS"//SN4(1:len_trim(SN4))//"_TT"//SN3(1:len_trim(SN3))&
                //"_S"//SNS(1:len_trim(SNS))

       call AvailableIOUnit(hFile1)
       open(UNIT=hFile1, FILE=OutFname(1:len_trim(OutFname))//".FixTime")

     !*** to initialize the DEVICE
       IDEV = COMMAND_ARGUMENT_COUNT()
       if(IDEV.GE.1) then
          call GET_COMMAND_ARGUMENT(1,ARGV)
          read(ARGV, *) IDEV
       end if
       call Initialize_DEVICES( IDEV, 1) 
       call GetSeed_RAND32SEEDLIB(ISEED0,ISEED(1),ISEED(2))
       call DRAND32_PUTSEED(ISEED) 
       call Initialize_Rand_DEVICES()

       allocate(PreTime(NB), NextTime(NB), NJUMP(NB), XP(NB), DDISH(-20*NBIN:20*NBIN),NJH(0:NJBIN))
       call DevAllocate(IDEV,dPreTime,   NB) 
       call DevAllocate(IDEV,dNextTime,  NB) 
       call DevAllocate(IDEV,dNJump,     NB) 
       call DevAllocate(IDEV,dXP,        NB) 

       XP      = 0.D0
       PreTime = 0.D0
       NextTime= 0.D0
       NJUMP   = 0
       call DevCopyIn(PreTime,  dPreTime,   NB)
       call DevCopyIn(NextTime, dNextTime,  NB)
       call DevCopyIn(NJUMP,    dNJump,     NB)
       call DevCopyIn(XP,       dXP,        NB)

       N = 0
       M = 0
       CURT   = 0.D0
       do while(CURT .le. TTIME)
          N = N+1
          if(N .ge. 10) then
             N = 1
             M = M + 1
          end if   
          CURT = TSTEP*dble(N)*10.D0**dble(M)
          if(CURT .gt. TTIME) exit
          !CURT = CURT + TSTEP
          print *, "CURT ", CURT, TTIME 
          call Stepto_RecTime_Pow_1D_Test(1, NB, Tau, Pow, DSTEP, CURT, dPreTime, dNextTime, dNJump, dXP)
          call DevCopyOut(PreTime,  dPreTime,   NB)
          call DevCopyOut(NextTime, dNextTime,  NB)
          call DevCopyOut(NJUMP,    dNJump,     NB)
          call DevCopyOut(XP,       dXP,        NB)
          
          !--- write out results
          MSD  = sum(XP(1:NB)*XP(1:NB))/dble(NB)
          MN   = sum(NJUMP)/dble(NB)
          MT   = sum(PreTime)/dble(NB)
          NJMI = minval(NJUMP)
          NJMX = maxval(NJUMP)

          if(MN .gt. 0.D0) then
             MN2 = sum((NJUMP(1:NB)-MN)*(NJUMP(1:NB)-MN))/dble(NB)
             MN2 = dsqrt(MN2)/MN
          end if    
          if(MT .gt. 0.D0) then
             MT2 = sum((PreTime(1:NB)-MT)*(PreTime(1:NB)-MT))/dble(NB)
             MT2 = dsqrt(MT2)/MT
          end if   
          write(hFile1,fmt="(20(1PE16.5,1X))") CURT, MSD, MN, NJMI, NJMX, MN2, MT, MT2     

          !--- write out the distance distribution and jump number distribution
          if(N .eq. 1) then
            L1     = LBOUND(DDISH,dim=1) 
            L2     = UBOUND(DDISH,dim=1) 
            DBIN   = dsqrt(MSD)/dble(NBIN)
            MD     = sum(XP(1:NB))/dble(NB)
            DDISH  = 0.D0 

            DJBIN  = max(1.D0,NJMX/NJBIN)
            NJH    = 0.D0
            do IB=1, NB
               IBIN = nint((XP(IB) - MD)/DBIN)
               if(IBIN .ge. L1 .and. IBIN .le. L2) then
                  DDISH(IBIN) = DDISH(IBIN) + 1.D0
               end if

               IBIN      = NJUMP(IB)/DJBIN
               NJH(IBIN) = NJH(IBIN) + 1.D0
            end do 
            DDISH = DDISH/sum(DDISH) 
            NJH   = NJH/sum(NJH)

            call STRCATI(OutFname1, OutFname,  "_DIS_N", N, 4)
            call STRCATI(OutFname1, OutFname1, "_M", M, 4)
            call AvailableIOUnit(hFile2)
            open(UNIT=hFile2, FILE=OutFname1(1:len_trim(OutFname1)))
            do IBIN=L1, L2
               write(hFile2, fmt="(1x,I,10(1x,1PE13.5))") IBIN, IBIN*DBIN, DDISH(IBIN)/DBIN, DDISH(IBIN)
            end do  
            close(hFile2)

            call STRCATI(OutFname1, OutFname,  "_NJH_N", N, 4)
            call STRCATI(OutFname1, OutFname1, "_M", M, 4)
            call AvailableIOUnit(hFile2)
            open(UNIT=hFile2, FILE=OutFname1(1:len_trim(OutFname1)))
            do IBIN=0, NJBIN
               write(hFile2, fmt="(1x,I,10(1x,1PE13.5))") IBIN, IBIN*DJBIN, NJH(IBIN)/DJBIN, NJH(IBIN)
            end do  
            close(hFile2)
          end if  
       end do
       
            
       deallocate(PreTime, NextTime, NJUMP, XP)
       close(hFile1)
       call DevDeAllocate(dPreTime) 
       call DevDeAllocate(dNextTime) 
       call DevDeAllocate(dNJump) 
       call DevDeAllocate(dXP) 

       stop   
 end program Walk_1D_MSD_Test2_Main

