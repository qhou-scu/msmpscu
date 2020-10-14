 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to simulate the displacement distribution for an particle walk on 1D lattice.
 !                  The purpose is to test if the distribution following Gaussian distribution.
 !
 !                  This program replace integrate the function of Walk-1D-MSD-Test2-Exp, and Walk-1D-MSD-Test2-Pow,
 !                  with the simulation parameters are inputted from command line, or file. 
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                   cat Walk-1D-MSD-Test1.stdin | Walk-1D-MSD-Test1.exe 
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2020-05 (Hou Qing, Sichuan university)
 ! 
 !                 The subroutines of the random-walk model are placed into a static library.  
 !
 !

 !**********************************************************************************

 !**********************************************************************************
 program Walk_1D_MSD_Test1_Main
    use MSM_MultiGPU_Basic
    use RAND32_MODULE
    use RAND32SEEDLIB_MODULE,only:GetSeed_RAND32SEEDLIB
    use MC_TypeDef_RandWalkCtrl
    use MC_TypeDef_RandWalker_Base
    use MC_RandWalk_Exp_GPU
    use MC_RandWalk_Pow_GPU
    use MiniUtilities
    implicit none
    integer::NB = 100000      !the number of boxes in test
    integer::ISEED0=12345, ISEED(2)=(/3555, 5443/)

    real(KINDDF)::MSD, CURT, ATT, MT, MN, MT2, MN2, &
                  DBIN, DJBIN, MD, NJMI, NJMX
    real(KINDDF),allocatable::XP(:), PreTime(:), NextTime(:),NJUMP(:), DDISH(:), NJH(:)
    character*256:: fname="Walk_1D_MSD_Test1", OutFname="", OutFname1 
    integer::hFile1,hFile2, IB, N, M, NBIN=100, L1, L2, NJBIN=1000, IBIN, FLAG
    character*18::SN1, SN2, SN3, SN4, SN5, SN6, SN7, SNS, SS
    !--- Device variable
     character*32::ARGV, STRTMP(1)
     character*256::cmdline
     type(RandWalkCtrl)::CtrlParam
     type(WalkerBase)  ::Walker
     integer::IDEV
     type(DevVec_DF)   ::dPreTime
     type(DevVec_DF)   ::dNextTime
     type(DevVec_DF)   ::dNJump
     type(DevVec_DF)   ::dXP

     !---
     !*** to initialize the DEVICE
      IDEV = COMMAND_ARGUMENT_COUNT()
      if(IDEV.GE.1) then
        call GET_COMMAND_ARGUMENT(1,ARGV)
        read(ARGV, *) IDEV
      end if
      call Initialize_DEVICES( IDEV, 1) 

      read(5,fmt="(A256)", end=100)  cmdline
      call extract_optstr(cmdline, "-","NB", 1, N, STRTMP)
      if(N .ge. 1) then
          NB = ISTR(STRTMP(1))
      end if    

      call extract_optstr(cmdline, "-","RS", 1, N, STRTMP)
      if(N .ge. 1) then
         ISEED0 = ISTR(STRTMP(1))
      end if    
      write(SN2,fmt="(I)") NB
      SN2 = adjustl(SN2)
      write(SNS,fmt="(I)") ISEED0
      SNS = adjustl(SNS)
      call GetSeed_RAND32SEEDLIB(ISEED0,ISEED(1),ISEED(2))
      call DRAND32_PUTSEED(ISEED) 
      call Initialize_Rand_DEVICES()

       do while(.true.)
          read(5,fmt="(A256)", end=100)  cmdline
          call Load_Rectime_RandWalkCtrl(cmdline, CtrlParam)
          call Load_RandWalker(cmdline, Walker) 
    
          write(SN4,fmt="(I)") CtrlParam%RectimeStepNum
          SN4 = adjustl(SN4)
          write(SN5,fmt="(I)") nint(CtrlParam%RectimeFixStep*10)
          SN5 = adjustl(SN5)

          select case(Walker%JumpType(1)%WTD_TYPE)
          case (CP_WALKTYPE_WTD_EXP)
                write(SN1,fmt="(I)") nint(Walker%JumpType(1)%WTD_TAU*100)
                SN1 = adjustl(SN1)
                if(CtrlParam%RectimeStepType .eq. CP_WALK_RECTIMEINC_LOG10) then
                   OutFname = "EXP_TAU"//SN1(1:len_trim(SN1))//"_NB"//SN2(1:len_trim(SN2))//"_RECLOGN"//SN4(1:len_trim(SN4)) &
                            //"_S"//SNS(1:len_trim(SNS))
                else if(CtrlParam%RectimeStepType .eq. CP_WALK_RECTIMEINC_LINEAR) then
                   OutFname = "EXP_TAU"//SN1(1:len_trim(SN1))//"_NB"//SN2(1:len_trim(SN2))//"_RECLIN"//SN5(1:len_trim(SN5))  &
                           //"_S"//SNS(1:len_trim(SNS))
                end if           

          case (CP_WALKTYPE_WTD_POW)
                write(SN1,fmt="(I)") nint(Walker%JumpType(1)%WTD_TAU*100)
                write(SN3,fmt="(I)") nint(Walker%JumpType(1)%WTD_Alpha*100)
                SN1 = adjustl(SN1)
                SN3 = adjustl(SN3)
                OutFname = "POW_TAU"//SN1(1:len_trim(SN1))//"_POW"//SN3(1:len_trim(SN3))
                if(Walker%JumpType(1)%DISD_TYPE .eq. CP_WALKTYPE_DIS_POWDIS) then
                  write(SN1,fmt="(I)") nint(Walker%JumpType(1)%DISD_AV*100)
                  write(SN3,fmt="(I)") nint(Walker%JumpType(1)%DISD_Beta*100)
                  SN1 = adjustl(SN1)
                  SN3 = adjustl(SN3)
                  OutFname = OutFname(1:len_trim(OutFname))//"POWDIS_"//SN1(1:len_trim(SN1))//"_"//SN3(1:len_trim(SN3)) 
                end if   
                if(CtrlParam%RectimeStepType .eq. CP_WALK_RECTIMEINC_LOG10) then
                    OutFname = OutFname(1:len_trim(OutFname))//"_NB"//SN2(1:len_trim(SN2))//"_RECLOGN"//SN4(1:len_trim(SN4))//"_S"//SNS(1:len_trim(SNS))
                else if(CtrlParam%RectimeStepType .eq. CP_WALK_RECTIMEINC_LINEAR) then
                     OutFname = OutFname(1:len_trim(OutFname))//"_NB"//SN2(1:len_trim(SN2))//"_RECLIN"//SN5(1:len_trim(SN5)) &
                             //"_S"//SNS(1:len_trim(SNS))
                end if           
            end select 
            call AvailableIOUnit(hFile1)
            open(UNIT=hFile1, FILE=OutFname(1:len_trim(OutFname))//".FixTime")


             allocate(PreTime(NB), NextTime(NB), NJUMP(NB), XP(NB), DDISH(-20*NBIN:20*NBIN),NJH(0:NJBIN))
             call DevAllocate(IDEV, dPreTime,   NB) 
             call DevAllocate(IDEV, dNextTime,  NB) 
             call DevAllocate(IDEV, dNJump,     NB) 
             call DevAllocate(IDEV, dXP,        NB) 

             XP      = 0.D0
             PreTime = 0.D0
             NextTime= 0.D0
             NJUMP   = 0
             call DevCopyIn(PreTime,  dPreTime,   NB)
             call DevCopyIn(NextTime, dNextTime,  NB)
             call DevCopyIn(NJUMP,    dNJump,     NB)
             call DevCopyIn(XP,       dXP,        NB)

             do while(.true.)
                call GetRectime_RandWalkCtrl(CtrlParam, CURT, Flag)
                if(CURT .gt. CtrlParam%RectimeEnd) then
                  exit
                end if  
                print *, "CURT ", CURT, CtrlParam%RectimeEnd 
                select case(Walker%JumpType(1)%WTD_TYPE)
                case (CP_WALKTYPE_WTD_EXP)
                      if(Walker%JumpType(1)%DISD_TYPE .eq. CP_WALKTYPE_DIS_POWDIS) then
                         call Stepto_RecTime_Exp_1D_WithLJ_Test(1, NB, Walker%JumpType(1)%WTD_TAU, Walker%JumpType(1)%DISD_AV,  &
                                                                Walker%JumpType(1)%DISD_Beta, CURT, dPreTime, dNextTime, dNJump, dXP)
                      else
                        call Stepto_RecTime_Exp_1D_Test(1, NB, Walker%JumpType(1)%WTD_TAU, Walker%JumpType(1)%DISD_AV, CURT, dPreTime, dNextTime, dNJump, dXP)
                      end if     
                                              
                case (CP_WALKTYPE_WTD_POW)
                      if(Walker%JumpType(1)%DISD_TYPE .eq. CP_WALKTYPE_DIS_POWDIS) then
                          call Stepto_RecTime_Pow_1D_WithLJ_Test(1, NB, Walker%JumpType(1)%WTD_TAU, Walker%JumpType(1)%WTD_Alpha, Walker%JumpType(1)%DISD_AV, &
                                                                 Walker%JumpType(1)%DISD_Beta, CURT, dPreTime, dNextTime, dNJump, dXP)
                      else
                          call Stepto_RecTime_Pow_1D_Test(1, NB, Walker%JumpType(1)%WTD_TAU, Walker%JumpType(1)%WTD_Alpha, Walker%JumpType(1)%DISD_AV, &
                                                         CURT, dPreTime, dNextTime, dNJump, dXP)  
                      end if    
                end select      
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
                if(CtrlParam%RectimeStepType .eq. CP_WALK_RECTIMEINC_LOG10 .and.  Flag .eq. 1) then !&
                   !CtrlParam%RecCurSteps .eq. 1  .or. mod(CtrlParam%RecCurSteps, CtrlParam%RectimeStepNum) .eq. 0) then
                   M      = CtrlParam%RecCurSteps/CtrlParam%RectimeStepNum 
                   !N      = CtrlParam%RecCurSteps - M*CtrlParam%RectimeStepNum  
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

                   call STRCATI(OutFname1, OutFname,  "_DIS_M", M, 4)
                   !call STRCATI(OutFname1, OutFname1, "_M", M, 4)
                   call AvailableIOUnit(hFile2)
                   open(UNIT=hFile2, FILE=OutFname1(1:len_trim(OutFname1)))
                   do IBIN=L1, L2
                      write(hFile2, fmt="(1x,I,10(1x,1PE13.5))") IBIN, IBIN*DBIN, DDISH(IBIN)/DBIN, DDISH(IBIN)
                   end do  
                   close(hFile2)

                   call STRCATI(OutFname1, OutFname,  "_NJH_M", M, 4)
                   !call STRCATI(OutFname1, OutFname1, "_M", M, 4)
                   call AvailableIOUnit(hFile2)
                   open(UNIT=hFile2, FILE=OutFname1(1:len_trim(OutFname1)))
                   do IBIN=0, NJBIN
                      write(hFile2, fmt="(1x,I,10(1x,1PE13.5))") IBIN, IBIN*DJBIN, NJH(IBIN)/DJBIN, NJH(IBIN)
                   end do  
                   close(hFile2)
                end if  
             end do
                  
             deallocate(PreTime, NextTime, NJUMP, XP, DDISH, NJH)
             close(hFile1)
             call DevDeAllocate(dPreTime) 
             call DevDeAllocate(dNextTime) 
             call DevDeAllocate(dNJump) 
             call DevDeAllocate(dXP) 
             call Reset_Rectime_RandWalkCtrl(CtrlParam)
       end do
  100  stop   
 end program Walk_1D_MSD_Test1_Main

