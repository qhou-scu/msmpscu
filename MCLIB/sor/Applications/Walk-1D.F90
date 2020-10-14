 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to simulate the diffusion of an particle walk on 1D lattice.
 !                  The purpose is to test if the distribution following Gaussian distribution.
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
 program Walk_1D__Main
    use MSM_MultiGPU_Basic
    use RAND32_MODULE
    use RAND32SEEDLIB_MODULE,only:GetSeed_RAND32SEEDLIB
    use MC_TypeDef_RandWalkCtrl
    use MC_TypeDef_RandWalker_Base
    use MC_RandWalk_Evolve
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


    integer::MXNP
    !--- Device variable
     character*256::cmdline
     character*32::ARGV
     character*256::CtrFile, WalkFile, STRTMP(1)
     type(RandWalkCtrl)  ::CtrlParam
     type(WalkerBase)    ::Walkers(1)

     integer::IDEV, LINE, ERR
     type(DevVec_DF)   ::dPreTime
     type(DevVec_DF)   ::dNextTime
     type(DevVec_DF)   ::dNJump
     type(DevVec_DF)   ::dXP

     !---
     !*** to initialize the DEVICE
      IDEV = 0
      if(command_argument_count() .ge. 1) then
        call get_command_argument(1,ARGV)
        read(ARGV, *) IDEV
      end if
      call Initialize_DEVICES(IDEV, 1) 

      do while(.true.)
         read(5,fmt="(A256)", end=100)  cmdline
         call extract_optstr(cmdline, "-","CF", 1, N, STRTMP)
         if(N .ge. 1) then
            CtrFile =STRTMP(1)
         end if    

         call extract_optstr(cmdline, "-","WF", 1, N, STRTMP)
         if(N .ge. 1) then
            WalkFile =STRTMP(1)
         end if   
         
         !--- load control paramter
         call LoadSetup_RandWalkCtrl(CtrFile, CtrlParam)

         !--- loal walker description
         call LoadSetup_RandWalker(WalkFile, Walkers)

         !--- some parameters can be replaced be command line arguments
         !    the restart flag
         call extract_optstr(cmdline, "-","RESTART", 1, N, STRTMP)
         if(N .ge. 1) then
            STRTMP(1) = adjustl(STRTMP(1))
            call Upcase(STRTMP(1))
            if(STRTMP(1) .eq. "Y" .or. STRTMP(1) .eq. "YES") then
                CtrlParam%Restart = 1
            else   
               CtrlParam%Restart  = 0
            end if     
         end if    

         !---- the output folder
         call extract_optstr(cmdline, "-","OUTTO", 1, N, STRTMP)
         if(N .ge. 1) then
            CtrlParam%RecPath = STRTMP(1)
         end if    

         !--- check input consistent before conduct evolution
         call CheckInput_WalkEvolve(CtrlParam, Walkers, ERR)

         !--- to initialize rand number
         ISEED0 = CtrlParam%SEED(1)
         call GetSeed_RAND32SEEDLIB(ISEED0, ISEED(1), ISEED(2))
         call DRAND32_PUTSEED(ISEED)
         call Initialize_Rand_DEVICES()
 
         !--- save the control parameter
         call Transfer_RandWalkCtrl(CtrlParam)
         call Transfer_RandWalkCtrl(CtrlParam, WalkFile)

         !--- start the evolution
         call Initialize_WalkEvolve(CtrlParam, Walkers)
         call ConductEvolve_WalkEvolve(CtrlParam, Walkers)
         
       end do
  100  stop   
 end program Walk_1D__Main

