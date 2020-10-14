 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to simulate the displacement distribution for an particle walk on 1D lattice.
 !                  The purpose is to test if the distribution following Gaussian distribution.
 !
 !                  The wating time distribution (WTD) of jump is exponential.
 !                   
 !                  This is a CPU version, in comparison with the GPU verision Walk-1D-MSD-Test2-Exp-GPU.F90.
 !
 !                  Note: the model parameters are programmatically set. Thus, this program is replaceed by
 !                        Walk-1D-MSD-Test1.
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                   Walk-1D-MSD-Test2-Exp.exe 
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2019-06 (Hou Qing, Sichuan university)
 !
 !                  
 !
 !

 !**********************************************************************************
 program Walk_1D_MSD_Test2_Main
    use MD_TYPEDEF_SimMDBox
    use RAND32_MODULE
    use RAND32SEEDLIB_MODULE,only:GetSeed_RAND32SEEDLIB
    implicit none
    integer, parameter::NB   = 10000      !the number of boxes in test
    integer::ISEED0=12345, ISEED(2)=(/3555, 5443/)

    real(KINDDF)::R, DT, TAU=1.D0, MSD, DSTEP=2., TSTEP=1.D0, TTIME=1000.D0, CURT, ATT, STIME, MT, MN, MT2, MN2, CURSTEP
    real(KINDDF),allocatable::XP(:), PreTime(:), NextTime(:)
    integer,     allocatable::NJUMP(:)
    character*256:: fname="Walk_1D_MSD_Test1", OutFname="" 
    integer::hFile1,hFile2, IB
    character*18::SN1, SN2, SN3, SN4, SN5, SNS


       TSTEP  = 10. ! (ps)
       TTIME  = 100000.D0
       TAU    = 0.1
       write(SN1,fmt="(I)") nint(TAU*100)
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
       
       OutFname = "EXP_TAU"//SN1(1:len_trim(SN1))//"_NB"//SN5(1:len_trim(SN5))//"_TS"//SN4(1:len_trim(SN4))//"_TT"//SN3(1:len_trim(SN3))&
                //"_S"//SNS(1:len_trim(SNS))

       call AvailableIOUnit(hFile1)
       open(UNIT=hFile1, FILE=OutFname(1:len_trim(OutFname))//".FixTime")

       !--- 
       allocate(PreTime(NB), NextTime(NB), NJUMP(NB), XP(NB))
       call GetSeed_RAND32SEEDLIB(ISEED0,ISEED(1),ISEED(2))
       call DRAND32_PUTSEED(ISEED) 

       XP      = 0.D0
       PreTime = 0.D0
       NextTime= 0.D0
       NJUMP   = 0

       CURT   = 0.D0
       do while(CURT .le. TTIME)
          CURT = CURT + TSTEP
          print *, "CURT ", CURT, TTIME 
          do IB=1, NB 
             if(CURT .gt. NextTime(IB)) then
               !--- jump the next step 
                if(NJUMP(IB) .gt. 0) then
                   R           = (DRAND32()-0.5)*DSTEP
                   XP(IB)      = XP(IB) + R
                   NJUMP(IB)   = NJUMP(IB) + 1
                   PreTime(IB) = NextTime(IB)
                end if  

                ATT    =  0.D0
                CURSTEP =  CURT - PreTime(IB) 
                do  while(ATT .le. CURSTEP)
                    !--  for exponetial WTD
                    R      = DRAND32()
                    DT     = -TAU*dlog(R)                
                    ATT    = ATT + DT
                    if(ATT .gt. CURSTEP) then
                       NextTime(IB) = PreTime(IB) + DT
                       exit
                    end if
                    !--- for displacement
                    R         = (DRAND32()-0.5)*DSTEP
                    XP(IB)    = XP(IB) + R
                    NJUMP(IB) = NJUMP(IB) + 1

                    !--- update the residule time
                    PreTime(IB)  = PreTime(IB) + DT
                    ATT    =  0.D0
                    CURSTEP =  CURT - PreTime(IB) 
                end do       
             end if   
          end do   
          
          !--- write out results
          MSD = 0.D0
          MN  = 0.D0
          MT  = 0.D0
          do IB =1, NB
             MSD = MSD + (XP(IB)*XP(IB))
             MN  = MN  +  NJUMP(IB)
             MT  = MT  + PreTime(IB) 
          end do
          MSD = MSD/dble(NB)
          MN  = MN/dble(NB)
          MT  = MT/dble(NB)
          MN2 = 0.D0
          MT2 = 0.D0
          do IB =1, NB
              MN2 = MN2 + (NJUMP(IB)-MN)*(NJUMP(IB)-MN)  
              MT2 = MT2 + (PreTime(IB)-MT)*(PreTime(IB)-MT)  
         enddo 
         MN2 = dsqrt(MN2/dble(NB))/MN
         MT2 = dsqrt(MT2/dble(NB))/MT
         
         write(hFile1,fmt="(20(1PE16.5,1X))") CURT, MSD, MN, MT, MN2, MT2     
       end do
            
       deallocate(PreTime, NextTime, NJUMP, XP)
       close(hFile1)
       stop   
 end program Walk_1D_MSD_Test2_Main

