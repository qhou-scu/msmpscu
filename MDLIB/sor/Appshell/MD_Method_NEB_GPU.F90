  module MD_Method_NEB_GPU
  !***  DESCRIPTION:
  !     This module provides routines for NEB analysis of diffusion path(s) between
  !     two substable stats.
  !
  !     ALGORITHM REF:
  !               S.A.Trygubenko and D.J. Wales, J. Chem Phys 120(2004)2082
  !               D. Sheppard,et al, J. Chem. Phys. 136(2012)074103
  !
  !     HISTORY:   created by HOU Qing, sep., 2016
  !     SEE ALSO:  MD_Method_GenericMD_GPU.F90,
  !                MD_ForceLib_Factory_GPU.F90

  !***  The modules included ******************************************
  use MD_SimboxArray
  use MD_TYPEDEF_SimMDCtrl
  use MD_TYPEDEF_Recordstamp
  use MD_TYPEDEF_RecordList
  implicit none

  contains
  !****************************************************************************************
  subroutine For_One_Test(SimBox0, CtrlParam0, SimBox, Recordlist, J, processid, INICONFIGPROC)
  !***  PORPOSE: to process one sample
  !     INPUT:  SimBox0,            the element sample
  !             CtrlParam,          the control parameters
  !             SimBox,             the simulation box array
  !             Recordlist          the external recoding routines
  !             J,                  the test ID
  !             processid,          id of the process, to be used if MPI used
  !             INICONFIGPROC,      the subroutine privided by user for initializing the initial configure
  !
  use MD_Globle_Variables
  use MD_NeighborsList_GPU
  use MD_ActiveRegion_GPU
  use MD_SimboxArray_GPU
  use MD_ForceLib_Factory_GPU

  implicit none
  !--- dummy variables
   type(SimMDBox),                          intent(in)::SimBox0
   type(SimMDCtrl),target                             ::CtrlParam0
   type(SimMDBox),dimension(:),allocatable            ::SimBox
   type(RecordProcedureList)                          ::Recordlist
   integer,                                 intent(in)::J,processid

   !--- interface for the external process--------------
   optional::INICONFIGPROC
   external::INICONFIGPROC

   interface
       SUBROUTINE INICONFIGPROC(SimBox, CtrlParam, RESTART)
       use MD_CONSTANTS
       use MD_SimboxArray
       use MD_TYPEDEF_SimMDCtrl
       implicit none
       type(SimMDBox), dimension(:)::SimBox
       type(SimMDCtrl)             ::CtrlParam
       integer                     ::RESTART
       END SUBROUTINE INICONFIGPROC
   end interface
  !------END INTERFACE ---

  !--- Local variables
  integer::ITIME0, ITIME, REACTLOOP0, REACTLOOP, IREACT, JOB, ISECT0, IB, I0, IFLAG
  real(KINDDF)::TIME, TIME0
  logical::EX
  character*256::SFILE, GFile, CFile, CFile0
  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2

  integer::NPATH, NSTAT, NIMG, IM0, IM, FINDNODE, MAXITER
  type(SimMDCtrl)::CtrlParam
  type(SimMDBox), dimension(:), allocatable::SimBoxSwap
  integer,        dimension(:), allocatable::ICFG, IBOX, NUMPATH


             !**** to determine how many inner-test loop we need
             !     for NEB, all boxes are independent. In NEB we
             !     use REACT instead of TEST
             JOB        = J
             REACTLOOP0 = -1
             REACTLOOP  = -1
             if(CtrlParam0%NEB_StartCfg .ge. 0) then
                REACTLOOP0 = CtrlParam0%NEB_StartCfg
                REACTLOOP  = CtrlParam0%NEB_EndCfg - CtrlParam0%NEB_CfgStep
                if(REACTLOOP .lt. 0) REACTLOOP = 999999
             end if

             !*** save the oringal control parameters for damping
             !    We will looking the state on the path.
             if(CtrlParam0%NEB_FindMin .gt. 0) then
                MAXITER = CtrlParam0%NEB_FindMin
                if(MAXITER .le. 1)  MAXITER =2
              else
                MAXITER = 1
             end if
             call Copy_SimMDCtrl(CtrlParam0, CtrlParam)

             !*** modify the control parameter for present NEB
             CtrlParam%MULTIBOX = CtrlParam0%MULTIBOX*CtrlParam0%NEB_NUMIMG
             if( CtrlParam%DAMPTIME1 .eq. 0)  CtrlParam%DAMPTIME1 = 1000

              !*** start test loop
              call CPU_TIME(C1)
              IREACT = REACTLOOP0
              do while(IREACT .le. REACTLOOP)

                 !*** modify the control parameter for present NEB
                 call Copy_SimMDCtrl(CtrlParam0, CtrlParam)
                 CtrlParam%MULTIBOX = CtrlParam0%MULTIBOX*CtrlParam0%NEB_NUMIMG

                 !*** create the ARRAY for start and end states
                 call Create_SimBoxArray(SimBoxSwap, CtrlParam0%MULTIBOX, SimBox0)
                 call Create_SimBoxArray(SimBox,     CtrlParam%MULTIBOX,  SimBox0)

                 !*** to prepare configurations for start state
                 if(IREACT .ge. 0) then
                  !--- load initial file from a serial of single files
                    call STRCATI(CFile0, SimBox0%IniConfig, "P", processid, 4)
                    call STRCATI(CFile0, CFile0, "_",  JOB, 4)
                    call STRCATI(CFile0, CFile0, ".",  IREACT, 4)
                    inquire(FILE=CFile0, EXIST=EX)
                    if(.not.EX) then
                        write(*,fmt="(A)")      ' MDPSCU Message: initial configuration'//CFILE0(1:len_trim(CFILE0))//' not found'
                        write(*,fmt="(A)")      '                 move to next test'
                        call Release_SimBoxArray(SimBoxSwap)
                        call Release_SimBoxArray(SimBox)
                        deallocate(SimBoxSwap, SimBox)
                        exit
                    end if
                    write(*,fmt="(A)")      ' MDPSCU Message: Load initial configuration from multi-file:'
                    call  Initialize_Config_SimBoxArray(SimBoxSwap, CFile0, CP_INPUT_POSONLY, multbox=1)
                 else
                    if(Simbox0%IniCfgID .gt. C_IZERO) then
                       call STRCATI(CFile0, SimBox0%IniConfig, "P", processid, 4)
                       call STRCATI(CFile0, CFile0, "_", JOB, 4)
                       call STRCATI(CFile0, CFile0, ".",  Simbox0%IniCfgID, 4)
                       write(*,fmt="(A)")      ' MDPSCU Message: Load initial configuration from multi-file:'
                       call  Initialize_Config_SimBoxArray(SimBoxSwap, CFile0, CP_INPUT_POSONLY, multbox=1)
                    else
                       CFile0 = SimBox0%IniConfig
                       write(*,fmt="(A)")      ' MDPSCU Message: Load initial configuration from single file:'
                       call Initialize_Config_SimBoxArray(SimBoxSwap, CFile0, CP_INPUT_POSONLY)
                    end if
                 end if

                 !*** to create the starting configurations for each path
                 do IB=1, CtrlParam0%MULTIBOX
                    SimBoxSwap(IB)%STATU = IOR(SimBox(IB)%STATU, CP_STATU_FIXPOS)
                    !*** copy the swap box to the starting box of each path in a box
                    I0 = (IB-1)*CtrlParam0%NEB_NUMIMG + 1
                    call Copy_SimMDBox(SimBoxSwap(IB),  SimBox(I0))
                 end do

                 !*** to prepare configurations for end state
                 if(IREACT .ge. 0) then
                    call STRCATI(CFile0, SimBox0%TgtConfig, "P", processid, 4)
                    call STRCATI(CFile0, CFile0, "_", JOB, 4)
                    call STRCATI(CFile0, CFile0, ".", IREACT+CtrlParam0%NEB_CfgStep, 4)
                    write(*,fmt="(A)")      ' MDPSCU Message: Load end configuration from multi-files'
                    write(*,fmt="(A, A)")   '                 ', CFILE0(1:len_trim(CFILE0))
                    inquire(FILE=CFile0, EXIST=EX)
                    if(.not.EX) then
                        write(*,fmt="(A)")      ' MDPSCU Message: target configuration'//CFILE0(1:len_trim(CFILE0))//' not found'
                        write(*,fmt="(A)")      '                 move to next test'
                        call Release_SimBoxArray(SimBoxSwap)
                        call Release_SimBoxArray(SimBox)
                        deallocate(SimBoxSwap, SimBox)
                        exit
                    end if
                    call Initialize_Config_SimBoxArray(SimBoxSwap, CFile0, CP_INPUT_POSONLY, multbox=1)
                 else
                    if(Simbox0%MultiTgtConfig .gt. C_IZERO) then
                       call STRCATI(CFile0, SimBox0%TgtConfig, "P", processid, 4)
                       call STRCATI(CFile0, CFile0, "_", JOB, 4)
                       call STRCATI(CFile0, CFile0, ".",  Simbox0%MultiTgtConfig, 4)
                       write(*,fmt="(A)")      ' MDPSCU Message: Load end configuration from multi-file:'
                       call Initialize_Config_SimBoxArray(SimBoxSwap, CFile0, CP_INPUT_POSONLY, multbox=1)
                    else
                       CFile0 = SimBox0%TgtConfig
                       write(*,fmt="(A)")      ' MDPSCU Message: Load end configuration from single file:'
                      call Initialize_Config_SimBoxArray(SimBoxSwap, CFile0, CP_INPUT_POSONLY)
                    end if
                 end if

                 !*** to create the end configurations for each path
                 do IB=1, CtrlParam0%MULTIBOX
                    SimBoxSwap(IB)%STATU = IOR(SimBox(IB)%STATU, CP_STATU_FIXPOS)
                    !*** copy the swap box to the end box of each path in a box
                    I0 = IB*CtrlParam0%NEB_NUMIMG
                    call Copy_SimMDBox(SimBoxSwap(IB),  SimBox(I0))
                    call Release_SimMDBox(SimBoxSwap(IB))
                 end do
                 deallocate(SimBoxSwap)

                 !--- now initialize the configurations of intermediate  images
                 call Initialize_NEBImages(SimBox0, SimBox, CtrlParam0, processid, JOB, IREACT, OUTPUT=CtrlParam0%NEB_OutIniCfg)

                 !*** To begin the NEB calculations ********************************************
                 !*** Initially, we search for the possible substable states on the path
                 !
                 allocate(NUMPATH(CtrlParam0%MULTIBOX))
                 NUMPATH    = 1
                 NPATH      = sum(NumPath)
                 FINDNODE  = CtrlParam0%NEB_FindMin
                 do ITIME=1, MAXITER
                    if(ITIME .eq. MAXITER) then
                       if(FINDNODE .gt. 0) then
                          write(*,fmt="(A, I2)") " MDPSCU Warning: minimum nodes have been changed in the last iteration "
                          call ONWARNING(gm_OnWarning)
                       end if
                       FINDNODE = 0
                    end if

                    !*** to initialize the GPU variable
                    call Initialize_Globle_Variables_DEV(SimBox, CtrlParam)

                    !*** to initialize force table on device
                    !    Note: m_ForceTable is defined in module MD_TYPEDEF_FORCELIB_GPU
                    call Init_Forcetable_Dev(SimBox, CtrlParam, gm_ForceClass)

                    !***  to give a initial neighbore list
                    call Initialize_NeighboreList_DEV(SimBox, CtrlParam)
                    call Cal_NeighBoreList_DEV(SimBox, CtrlParam)

                    !*** if active region method to be used,  initialize the module
                    if(iand(CtrlParam%AR_METHOD,CP_ENABLE_AR) .eq. CP_ENABLE_AR) then
                       call Initialize_ActiveRegion_DEV(SimBox, CtrlParam)
                       call ActivateRegion_DEV(SimBox, CtrlParam)
                       call CopyOut_SimBox_DEV(SimBox)
                    end if

                    !*** start do NEB  for time sections
                    IFLAG = 0
                    write(*,fmt="(A,I5)") ' MDPSCU Message: do NEB path searching for iteration #', ITIME
                    !*** queching the boxes and then check if there are more stable states
                    if(FINDNODE .gt. 0 ) then

                       allocate(ICFG(size(SimBox)), IBOX(size(SimBox)) )
                       !*** to find out path
                        call Do_NEBLBFGS_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, CtrlParam%DAMPTIME1, IFLAG, HASSPRING=1, METH=CtrlParam0%NEB_Meth,CLIMB=0)

                       !*** to find out if there is intermediate minia on the path
                       write(*,fmt="(A,A)") ' MDPSCU Message: to check if intermediate substable states exist'

                       !**** quench the boxes without spring force
                       !***  before quench, do we need introducing some disturbing

                       !IFLAG = 0
                       !call Do_NEBLBFGS_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, CtrlParam%DAMPTIME1, IFLAG, HASSPRING=0, METH=CtrlParam0%NEB_Meth, CLIMB=0)
                       call Do_Damp(SimBox,CtrlParam, gm_ForceClass)

                       call CopyOut_SimBox_DEV(SimBox)
                       call CheckStableState(SimBox, CtrlParam0, NUMPATH, NSTAT, ICFG, IBOX, CtrlParam0%STRCUT_DRTol)
                       write(*,fmt="(A, I6, A, I6, A)") ' MDPSCU Message: ', NSTAT, ' stable states found in ', &
                                                                             CtrlParam0%MULTIBOX, ' boxes'

                      !$$--- recaculate the number of paths
                       do IB=1, CtrlParam0%MULTIBOX
                          NUMPATH(IB)    = count(IBOX(1:NSTAT) .eq. IB) - 1
                       end do

                       if(sum(NUMPATH) .eq. NPATH) then
                          !$$--- no new stable states found, we restore the climbing
                          FINDNODE = 0
                       end if
                       NPATH = sum(NUMPATH)
                       if(NPATH .lt. C_IUN) then
                          write(*,fmt="(A, I6, A, I6, A)") ' MDPSCU Error: no transition path '
                          write(*,fmt="(A, F8.2,1x, A)")   '               the value of &DRTOL ',CtrlParam%STRCUT_DRTol/SimBox(1)%RR, ' in the control file is probably too large '
                          write(*,fmt="(A, F8.2,1x, A)")   '               or  the end-states are not stable, check the confgure files'
                          write(*,fmt="(A, I6, A, I6, A)") ' Process to be stopped '
                          stop
                       end if

                       !$$--- save the stable states
                       call  Create_SimBoxArray(SimBoxSwap, NSTAT)
                       do IM0=1, NSTAT
                          call Copy_SimMDBox(SimBox(ICFG(IM0)), SimBoxSwap(IM0))
                          SimBoxSwap(IM0)%STATU = ior(SimBox(IM0)%STATU, CP_STATU_FIXPOS)
                       end do

                       !$$--- re-initialize the image boxes
                       CtrlParam%MULTIBOX = CtrlParam0%NEB_NUMIMG*NPATH
                       call Create_SimBoxArray(SimBox, CtrlParam%MULTIBOX, SimBox0)

                         !$$--- reset the node images
                         IM0 = 1
                         IM  = 1
                         do IB=1, CtrlParam0%MULTIBOX
                            do I0=1, NUMPATH(IB)
                               call Copy_SimMDBox(SimBoxSwap(IM0),   SimBox(IM) )
                               call Copy_SimMDBox(SimBoxSwap(IM0+1), SimBox(IM+CtrlParam0%NEB_NUMIMG-1 ) )
                               call Init_NEBImages_Linear(SimBox(IM:IM+CtrlParam0%NEB_NUMIMG-1), CtrlParam0)
                               IM0 = IM0 + 1
                               IM  = IM  + CtrlParam0%NEB_NUMIMG
                            end do
                            IM0  = IM0 + 1
                         end do

                       call Release_SimBoxArray(SimBoxSwap)
                       deallocate(SimBoxSwap)
                       deallocate(ICFG,IBOX)
                       !$$--- the number boxes may be changed
                       CtrlParam%MULTIBOX   = size(SimBox)
                       CtrlParam%NEB_NUMIMG = CtrlParam0%NEB_NUMIMG
                       !$$--- the images are re-initialized
                    else
                      call Do_NEBLBFGS_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, CtrlParam%DAMPTIME1, IFLAG, HASSPRING=1, METH=CtrlParam0%NEB_Meth, CLIMB=CtrlParam0%NEB_Climb)
                      exit
                    end if
                 end do
                !--- recalculate the static energy of images on the path
                call Cal_NeighBoreList_DEV(SimBox, CtrlParam)
                call CalEpot_ForceClass(SimBox, CtrlParam, gm_ForceClass)
                call CopyOut_SimBox_DEV(SimBox)
                if(IREACT .eq. REACTLOOP0) then
                   call Putout_NEBImages(SimBox, CtrlParam0, JOB, NUMPATH,Flag=0)
                   call Putout_NEBEnergy(SimBox, CtrlParam0, JOB, NUMPATH,Flag=0, Statu=IFLAG)
                else
                   call Putout_NEBImages(SimBox, CtrlParam0, JOB, NUMPATH,Flag=1)
                   call Putout_NEBEnergy(SimBox, CtrlParam0, JOB, NUMPATH,Flag=1, Statu=IFLAG)
                end if
                if(allocated(NUMPATH)) deallocate(NUMPATH)

                !---
                IREACT = IREACT + CtrlParam0%NEB_CfgStep
              end do    !end the loop for ireact
              !--- to clear allocated memory by set Flag = -1
              call Putout_NEBImages(SimBox, CtrlParam0, JOB, NUMPATH, Flag= -1)
              call Putout_NEBEnergy(SimBox, CtrlParam0, JOB, NUMPATH, Flag= -1, Statu=IFLAG)

            !---- record the end time
             call CPU_TIME(C2)
             print *, "RUN TIME FOR ONE TEST: ",C2-C1


      return
  end subroutine For_One_Test
  !*********************************************************************************

  !*********************************************************************************
  subroutine Init_NEBImages_Linear(SimBox, CtrlParam)
  !***  PURPOSE:   to initialize the NEB images, by linearly insert
  !                the images between the two end point
  !
  !     INPUT:     SimBox,    the simulation boxes
  !                CtrlParam: the control parameters
  !
  !     OUTPUT     SimBox,    the simulation box with positions updated
  !
  !
  !
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:),intent(inout):: SimBox
      type(SimMDCtrl),             intent(in)   :: CtrlParam

      !--- Local variables
           integer::NB, IB, NPRT, IA, IFPD(3),K
           real(KINDDF):: BS(3), HBS(3), LB(3), HB(3), DIS, DEL

             NPRT      = SimBox(1)%NPRT
             BS        = SimBox(1)%ZL
             HBS       = C_HALF*BS
             LB        = SimBox(1)%BOXLOW
             HB        = LB + BS
             IFPD      = CtrlParam%IFPD

             NB = size(SimBox)
             do IB = 2, NB-1
                DEL = dble(IB-1)/dble(NB-1)
                do IA=1, NPRT
                   do K=1, 3
                      DIS = SimBox(NB)%XP(IA,K) - SimBox(1)%XP(IA,K)
                      if((IFPD(K).gt.0) .and. (dabs(DIS) .gt. HBS(K)) ) then
                          DIS = DIS - DSIGN(BS(K),DIS)
                      end if
                      DIS = SimBox(1)%XP(IA, K) + DIS*DEL
                      if(IFPD(K) .gt. 0) then
                         if(DIS .lt. LB(K)) then
                            DIS = DIS + BS(K)
                          else if(DIS .gt. HB(K)) then
                            DIS = DIS - BS(K)
                         end if
                      end if
                      SimBox(IB)%XP(IA,K) = DIS
                   end do
                end do
             end do

      return
  end subroutine Init_NEBImages_Linear
  !*********************************************************************************

  !*********************************************************************************
  subroutine Initialize_NEBImages(SimBox0, SimBox, CtrlParam, processid, JOB, IREACT, OUTPUT)
  !***  PURPOSE:   to initialize the NEB images
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !     OUTPUT     SimBox,    the simulation box with force updated
  !
  !
  !
   implicit none
      !--- dummy vaiables
      type(SimMDBox)                            :: SimBox0
      type(SimMDBox), dimension(:),intent(inout):: SimBox
      type(SimMDCtrl),             intent(in)   :: CtrlParam
      integer,                     intent(in)   :: processid, JOB, IREACT
      integer,                     intent(in)   :: OUTPUT

      !--- Local variables
           character*256::GFile, GFile0
           integer::NIMG, IB, IMG, IB0, IB1, NS, I0, I1, IEND
           logical::EX
           type(MDRecordStamp)::STAMP

                 !*** to prepare configurations for start state
                 NIMG = CtrlParam%NEB_NumImg
                 if(IREACT .ge. 0) then
                    IEND = IREACT + CtrlParam%NEB_CfgStep
                    NS   = CtrlParam%NEB_NumImg/CtrlParam%NEB_CfgStep
                    if(NS .lt. 1) NS = 1

                    I0  =  IREACT
                    I1  =  I0 + 1
                    IB0 =  1
                    IB1 =  min(IB0 + NS, NIMG)
                    do while(I1 .lt. IEND)  
                       !--- load initial file from a serial of single files
                       call STRCATI(GFile0, SimBox0%IniConfig, "P", processid, 4)
                       call STRCATI(GFile0, GFile0, "_",  JOB, 4)
                       call STRCATI(GFile0, GFile0, ".",  I1, 4)
                       inquire(FILE=GFile0, EXIST=EX)
                       if(.not.EX) then
                          write(*,fmt="(A)")      ' MDPSCU Message: intermdiate configuration'//GFILE0(1:len_trim(GFILE0))//' not found'
                          write(*,fmt="(A)")      '                 processto be stopped'
                          stop
                       end if
                       write(*,fmt="(A)")      ' MDPSCU Message: Load intermdiate configuration from multi-files'
                       !NOTE: for serial NEB calculations, only single box are loaded
                       call Initialize_Config_SimBoxArray(SimBox(IB1:IB1), fname=GFile0, fmt=CP_INPUT_POSONLY,multbox=1)
                       call Init_NEBImages_Linear(SimBox(IB0:IB1),CtrlParam)
                       I0  = I1
                       I1  = I0 + 1
                       IB0 = IB1
                       IB1 = min(IB0 + NS, NIMG)
                    end do
                    call Init_NEBImages_Linear(SimBox(IB0:IB1),CtrlParam)
                    
                 else
                    do IB=1, CtrlParam%MULTIBOX
                       IB0 = (IB-1)*NIMG + 1
                       IB1 = IB0 + NIMG - 1
                       !--- for linearly inserted images
                      call Init_NEBImages_Linear(SimBox(IB0:IB1),CtrlParam)
                    end do
                 end if


                 if(OUTPUT .gt. 0) then
                   !--- to output initial paths
                   IB0 = (JOB-1)*CtrlParam%MULTIBOX
                   call GetPath(CtrlParam%f_geometry, GFILE0)
                   call STRCATI(GFILE0, GFILE0, "NEBIntial_P", 0, 4)

                   STAMP%AppType   = gm_AppType
                   STAMP%ITest     = JOB
                   STAMP%ITime     = 0
                   STAMP%ISect     = 1
                   STAMP%Time      = 0
                   STAMP%ScalTime  = 0
                   do IB=1, CtrlParam%MULTIBOX
                      do I0 = 1, NIMG
                         call STRCATI(GFILE, GFILE0, "_box",  IB0+IB, 4)
                         STAMP%IBox = IB0+IB
                         STAMP%ICfg = I0
                         STAMP%IRec = STAMP%ICfg
                         call Putout_Instance_Config_SimMDBox(GFILE, SimBox((IB-1)*NIMG+I0), Stamp)
                      end do
                   end do
                 end if

      return
  end subroutine Initialize_NEBImages
  !*********************************************************************************

  !*********************************************************************************
  subroutine NEBSpringForce_CPU(SimBox, CtrlParam, Ipath, Epot, Fp, Spot, Meth, Climb)
  !***  PURPOSE:   to calculate the "spring" between images forces for a path
  !                and add to the real force
  !
  !     INPUT:     SimBox,    the simulation box (images)
  !                CtrlParam: the control parameters
  !                IPATH:     the path ID
  !                Epot:      the potential
  !                Meth:
  !                CLIMB;     the flag indicating if climbing eneabled
  !
  !     OUTPUT:    FP:        the modified force
  !                SPOT:      the potential of the "spring"
  !
  !     NOTE:   Positions of atoms, m_XP in the original box should be copy out before calling this
  !             subroutine
  !
   use MD_Globle_Variables_GPU, only:m_XP
   implicit none
      !--- dummy vaiables
      type(SimMDBox),dimension(:)              ::SimBox
      type(SimMDCtrl),               intent(in)::CtrlParam
      integer,                       intent(in)::Ipath
      real(KINDDF),dimension(:)                ::Epot
      real(KINDDF),dimension(:,:)              ::Fp
      real(KINDDF)                             ::Spot
      character*(*)                            ::Meth
      integer,                       intent(in)::Climb

      !--- Local variables
           integer::NIMG, IM, NPRT0, IERR, IFPD(3), FROMA, IH, IL, &
                    I, K, IA1, IA2, IAM
           real(KINDDF):: KValue, NORM, BS(3),HBS(3), SEP1(3), SEP2(3), DISH, DISL, DVMIN, DVMAX
           real(KINDDF), dimension(:),  allocatable::IMGPOT
           real(KINDDF), dimension(:,:),allocatable::FS, TANGENT, SVECL, SVECH
           integer,      dimension(:),  allocatable::CLFLG
     !----
             NIMG   =  CtrlParam%NEB_NUMIMG
             NPRT0  =  SimBox(1)%NPRT
             BS     =  SimBox(1)%ZL
             HBS    =  C_HALF*BS
             IFPD   =  CtrlParam%IFPD
             FROMA  =  (Ipath-1)*NIMG*NPRT0 + 1
             allocate(IMGPOT(NIMG), CLFLG(NIMG), FS(NPRT0,3), TANGENT(NPRT0,3), SVECH(NPRT0,3), SVECL(NPRT0,3),stat=IERR)
             !$$--- to extract current potentials of images
              IA1    = FROMA
              do IM=1, NIMG
                 IA2        = IA1 + NPRT0 -1
                 IMGPOT(IM) = sum(EPOT(IA1:IA2) - EPOT(FROMA:FROMA+NPRT0-1))
                 IA1        = IA1 + NPRT0
              end do

             !$$--- to determine the KValue
              if( CtrlParam%NEB_KValue .gt. 0) then
                  KValue =  CtrlParam%NEB_KValue*CP_EVERG/(CP_A2CM*CP_A2CM)
              else
                  KValue =  sum(dabs(EPOT(FROMA:FROMA+NPRT0*NIMG)))/dble(NPRT0*NIMG)/(CP_A2CM*CP_A2CM)
              endif

             !$$ to determine which images to be climbing
              CLFLG = 0
              if(Climb .gt. 0) then
                 do IM=2, NIMG-1
                    if(IMGPOT(IM).gt.IMGPOT(IM-1) .and. IMGPOT(IM).gt.IMGPOT(IM+1) ) then
                       select case(Climb)
                              case(C_IUN)
                                   CLFLG(IM)   = 1
                              case(C_ITWO)
                                   CLFLG(IM-1) = 1
                                   ClFlG(IM+1) = 1
                              case(C_ITHR)
                                   CLFLG(IM-1) = 1
                                   CLFLG(IM)   = 1
                                   ClFlG(IM+1) = 1
                       end select
                    end if
                 end do
              end if

              SPOT = 0.D0
              do IM=2, NIMG-1
                 !$$--- to find out the neighbore images of higher potential
                 if(IMGPOT(IM-1) .le. IMGPOT(IM+1) ) then
                    IH = IM + 1
                    IL = IM - 1
                 else
                    IH = IM - 1
                    IL = IM + 1
                 end if
                 DVMAX = max(dabs(IMGPOT(IH)-IMGPOT(IM)), dabs(IMGPOT(IL)-IMGPOT(IM)) )
                 DVMIN = min(dabs(IMGPOT(IH)-IMGPOT(IM)), dabs(IMGPOT(IL)-IMGPOT(IM)) )

                 !$$--- construct the "spring" vector
                 IA1 = FROMA+(IH-1)*NPRT0
                 IAM = FROMA+(IM-1)*NPRT0
                 IA2 = FROMA+(IL-1)*NPRT0
                 do I=1, NPRT0
                    do K=1, 3
                       SVECH(I,K) = m_XP(IA1, K) - m_XP(IAM, K)
                       SVECL(I,K) = m_XP(IA2, K) - m_XP(IAM, K)
                       if((IFPD(K).GT.0)) then
                          if(dabs(SVECH(I,K)) .GT. HBS(K)) then
                             SVECH(I,K) = SVECH(I,K) - DSIGN(BS(K),SVECH(I,K))
                          end if
                          if(dabs(SVECL(I,K)) .GT. HBS(K)) then
                             SVECL(I,K) = SVECL(I,K) - DSIGN(BS(K),SVECL(I,K))
                          end if
                       end if
                    end do
                    IA1 = IA1 + 1
                    IAM = IAM + 1
                    IA2 = IA2 + 1
                 end do
                 DISH    = dsqrt(sum(SVECH * SVECH))
                 DISL    = dsqrt(sum(SVECL * SVECL))

                 !$$--- to accumulate the spring force
                 IAM = FROMA+(IM-1)*NPRT0
                 if(CLFLG(IM)) then
                     !--- do climbing
                     TANGENT =  SVECH*DVMAX - SVECL*DVMIN
                     TANGENT = TANGENT/dsqrt(sum(TANGENT * TANGENT))
                     FP(IAM:IAM+NPRT0-1,1:3) = FP(IAM:IAM+NPRT0-1,1:3) -  &
                                               C_TWO * sum(FP(IAM:IAM+NPRT0-1,1:3)*TANGENT(1:NPRT0,1:3)) * TANGENT(1:NPRT0,1:3)

                 else !--- no climbing
                    select case(Meth(1:len_trim(Meth)))
                           case("DNEB")
                                TANGENT =  SVECH - SVECL
                                TANGENT =  TANGENT/dsqrt(sum(TANGENT * TANGENT))

                               !$$--- to calculate the verticel force
                                FP(IAM:IAM+NPRT0-1,1:3) = FP(IAM:IAM+NPRT0-1,1:3) - &
                                                          sum(FP(IAM:IAM+NPRT0-1,1:3)*TANGENT(1:NPRT0,1:3))*TANGENT(1:NPRT0,1:3)

                                !$$--- calculate the perpendicular spring force
                                FS(1:NPRT0,1:3) = KValue*(SVECH(1:NPRT0,1:3) + SVECL(1:NPRT0,1:3))
                                FS(1:NPRT0,1:3) = FS(1:NPRT0,1:3) - sum(FS(1:NPRT0,1:3)*TANGENT(1:NPRT0,1:3)) * TANGENT(1:NPRT0,1:3)

                                !$$-- accumulating parallel spring force
                                 FP(IAM:IAM+NPRT0-1,1:3)  = FP(IAM:IAM+NPRT0-1,1:3) + FS(1:NPRT0,1:3) + KValue*(DISH - DISL)*TANGENT(1:NPRT0,1:3)
                                 SPOT                     = SPOT                    + C_HALF*KValue*(DISL*DISL+DISH*DISH)

                           case("NEB")
                               TANGENT =  SVECH
                               TANGENT =  TANGENT/dsqrt(sum(TANGENT * TANGENT))
                               !$$--- to calculate the verticel force
                                FP(IAM:IAM+NPRT0-1,1:3) = FP(IAM:IAM+NPRT0-1,1:3) - &
                                                          sum(FP(IAM:IAM+NPRT0-1,1:3)*TANGENT(1:NPRT0,1:3))*TANGENT(1:NPRT0,1:3)

                               !$$-- accumulating the parallel spring force
                                FP(IAM:IAM+NPRT0-1,1:3) = FP(IAM:IAM+NPRT0-1,1:3) + KValue*(DISH - DISL)*TANGENT(1:NPRT0,1:3)
                                SPOT                    = SPOT                       + C_HALF*KValue*(DISH-DISL)*(DISH-DISL)

                           case("EB")
                              !$$-- accumulating spring force
                                FS(1:NPRT0,1:3)         = KValue*(SVECH(1:NPRT0,1:3) + SVECL(1:NPRT0,1:3))
                                FP(IAM:IAM+NPRT0-1,1:3) = FP(IAM:IAM+NPRT0-1,1:3)    + FS(1:NPRT0,1:3)
                                SPOT                    = SPOT                       + C_HALF*KValue*(DISL*DISL+DISH*DISH)

                           case default
                               write(*,fmt="(A)") " MDPSCU Error: unknown NEB method: "//Meth(1:len_trim(Meth))
                               write(*,fmt="(A)") "               process to be stopped"
                               stop
                    end select
                 end if
              end do

             deallocate(IMGPOT, CLFLG, FS, TANGENT, SVECH, SVECL)

      return
  end subroutine NEBSpringForce_CPU
  !*********************************************************************************

  !****************************************************************************************
  subroutine Do_NEBLBFGS_Forsteps_DEV(SimBox,CtrlParam, ForceClass, MXnumsteps, Iflag, Hasspring, Meth, Climb)
  !***  PORPOSE:  to forward one half of a time step
  !
  !     INPUT:     SimBox,      the simulation box
  !                CtrlParam,   the control parameters
  !                ForceClass,  the force clase
  !                MXNUMSTEPS,  the max number of calling LBFGS
  !                IFLAG,       flag indicating if the first time calling LBFGS
  !                HASSPRING,   flag indicating if spring force to be included
  !                             if HASSPRING = 0, minima searching to be performed
  !                CLIMB,       flag indicating if climbing to be performed
  !     OUTPUT     IFSD,        the flag indicating if a image is a saddle
  !
  use MATH_LBFGSB_MODULE, only:SETULB
  use MD_SimboxArray_GPU
  use MD_NeighborsList_GPU
  use MD_Forceclass_Register_GPU

  implicit none
  !----   DUMMY Variables
          type(SimMDBox),dimension(:)             :: SimBox
          type(SimMDCtrl),              intent(in):: CtrlParam
          type(MDForceClassGPU),target, intent(in):: ForceClass
          integer,                      intent(in):: MXNUMSTEPS
          integer                                 :: IFLAG
          integer                                 :: Hasspring
          character*(*)                           :: Meth
          integer                                 :: Climb
  !---   Local variables
          integer::I, IA, IP, NPRT0, NPRT, NPRTP, IERR, IFPD(3), NIMG, NPATH, IM, IA1
          integer, dimension(:), allocatable::IDX, IDY, IDZ
          integer::NIDX, NIDY, NIDZ

          real(KINDDF)::LB(3), UB(3), BS(3), OBJF, SPOT, EPREV, EMID, ENEXT
          real(KINDDF), dimension(:),   allocatable::EPOT
          real(KINDDF), dimension(:,:), allocatable::FP      !the gradient
!---      Declare the variables needed by LBFGSB.
!
          character*60::TASK, CSAVE
          logical     :: LSAVE(4)
          integer     :: IPRINT,ISAVE(44), NMAX, IT, MSAVE = 13
          real(KINDDF)::FACTR, PGTOL, DSAVE(29)
          integer,      dimension(:), allocatable::NBD    !the boundary type
          integer,      dimension(:), allocatable::IWA    !the integer working space
          real(KINDDF), dimension(:), allocatable::X      !the variables
          real(KINDDF), dimension(:), allocatable::G      !the gradient
          real(KINDDF), dimension(:), allocatable::L      !the low boundary of X
          real(KINDDF), dimension(:), allocatable::U      !the upper boundary of X
          real(KINDDF), dimension(:), allocatable::WA     !the working space
          
 !---
                if(HASSPRING) then
                   FACTR    = CtrlParam%NEB_LBFGS_Factr !1.0d+1
                   PGTOL    = CtrlParam%NEB_LBFGS_PGtol !1.0d-12
                   MSAVE    = CtrlParam%NEB_LBFGS_MSave
                else
                   FACTR    = CtrlParam%LBFGS_Factr
                   PGTOL    = CtrlParam%LBFGS_PGtol
                   MSAVE    = CtrlParam%LBFGS_MSave
                end if
                IPRINT   = -1
                IFLAG    = 0

                NPRT0    = SimBox(1)%NPRT
                NPRT     = dm_NPRT
                IFPD     = CtrlParam%IFPD
                LB       = SimBox(1)%BOXLOW
                UB       = SimBox(1)%BOXUP
                BS       = SimBox(1)%ZL
                NIMG     = CtrlParam%NEB_NUMIMG
                NPATH    = CtrlParam%MULTIBOX/NIMG

                NMAX = NPRT*3
                allocate(IDX(NPRT), IDY(NPRT), IDZ(NPRT), NBD(NMAX), IWA(3*NMAX), X(NMAX), L(NMAX), U(NMAX), G(NMAX), &
                         WA(2*MSAVE*NMAX + 5*NMAX + 11*MSAVE*NMAX + 8*MSAVE),STAT=IERR )
                allocate(EPOT(NPRT),FP(NPRT,3),STAT=IERR )
                EPOT = 0.D0
                FP   = 0.D0

                if(IERR .gt. 0) then
                    write(*,fmt="(A)") " MDPSCU Error: fail to allocate working space in DO_LBFGSB_FORSTEP_DEV"
                    write(*,fmt="(A)") "               process to be stopped"
                    stop
                end if

                NIDX = 0
                do IA = 1, NPRT
                   if(iand(m_STATU(IA), CP_STATU_FIXPOSX) .eq. 0 .and. &
                      iand(m_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(m_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDX = NIDX + 1
                      IDX(NIDX) = IA
                   end if
                end do

                NIDY = 0
                do IA = 1, NPRT
                   if(iand(m_STATU(IA), CP_STATU_FIXPOSY) .eq. 0 .and. &
                      iand(m_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(m_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDY = NIDY + 1
                      IDY(NIDY) = IA
                   end if
                end do

                NIDZ = 0
                do IA = 1, NPRT
                   if(iand(m_STATU(IA), CP_STATU_FIXPOSZ) .eq. 0 .and. &
                      iand(m_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(m_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDZ = NIDZ + 1
                      IDZ(NIDZ) = IA
                   end if
                end do

                X(1:NIDX)                      = m_XP(IDX(1:NIDX),1)
                X(NIDX+1:      NIDX+NIDY)      = m_XP(IDY(1:NIDY),2)
                X(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = m_XP(IDZ(1:NIDZ),3)
                NPRTP = NIDX + NIDY + NIDZ

                !$$--- set the boundary
                !$$    we constrain the variables changes in curoff ranges
                !$$    NOTE: it cannot work correctly if constrain is applied
                NBD(1:NPRTP)                   = 0
                L(1:NIDX)                      = LB(1) - BS(1)*C_HALF
                U(1:NIDX)                      = UB(1) + BS(1)*C_HALF
                L(NIDX+1:      NIDX+NIDY)      = LB(2) - BS(2)*C_HALF
                U(NIDX+1:      NIDX+NIDY)      = UB(2) + BS(2)*C_HALF
                L(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = LB(3) - BS(3)*C_HALF
                U(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = UB(3) + BS(3)*C_HALF
                TASK = 'START'
                do I=1, MXNUMSTEPS
                   !call LBFGS(NPRTP,M, X, OBJF, G, DIAG, IPRINT,EPS,XTOL,WORK,IFLAG)
                   !--- call to the L-BFGS-B code
                   call SETULB(NPRTP,MSAVE,X,L,U,NBD,OBJF,G,FACTR,PGTOL,WA,IWA,TASK,IPRINT,CSAVE,LSAVE,ISAVE,DSAVE)
                   !print *, "TASK", I, TASK(1:len_trim(TASK))

                     !----
                      if(IFPD(1)) then
                         do IA=1, NIDX
                            if(X(IA) .lt. LB(1) )then
                               m_XP(IDX(IA), 1) = X(IA) + BS(1)
                            else if(X(IA) .gt. UB(1) )then
                               m_XP(IDX(IA), 1) = X(IA) - BS(1)
                            else
                               m_XP(IDX(IA), 1) = X(IA)
                            end if
                         end do
                      else
                         m_XP(IDX(1:NIDX), 1) = X(1:NIDX)
                      end if

                      if(IFPD(2)) then
                         do IA=NIDX+1, NIDX+NIDY
                            if(X(IA) .lt. LB(2) )then
                               m_XP(IDY(IA-NIDX), 2) = X(IA) + BS(2)
                            else if(X(IA) .gt. UB(2) )then
                               m_XP(IDY(IA-NIDX), 2) = X(IA) - BS(2)
                            else
                               m_XP(IDY(IA-NIDX), 2) = X(IA)
                            end if
                         end do
                      else
                         m_XP(IDY(1:NIDY), 2) = X(NIDX+1:NIDX+NIDY)
                      end if

                      if(IFPD(3)) then
                         do IA=NIDX+NIDY+1, NIDX+NIDY+NIDZ
                            if(X(IA) .lt. LB(3) )then
                               m_XP(IDZ(IA-(NIDX+NIDY)), 3) = X(IA) + BS(3)
                            else if(X(IA) .gt. UB(3) )then
                               m_XP(IDZ(IA-(NIDX+NIDY)), 3) = X(IA) - BS(3)
                            else
                               m_XP(IDZ(IA-(NIDX+NIDY)), 3) = X(IA)
                            end if
                         end do
                      else
                         m_XP(IDZ(1:NIDZ),3) = X(NIDX+NIDY+1: NIDX+NIDY+NIDZ)
                      end if
                      !m_XP(IDX(1:NIDX), 1) = X(1:           NIDX)
                      !m_XP(IDY(1:NIDY), 2) = X(NIDX+1:      NIDX+NIDY)
                      !m_XP(IDZ(1:NIDZ), 3) = X(NIDX+NIDY+1: NIDX+NIDY+NIDZ)
                      call CopyXPFrom_Host_to_Devices(m_XP)
                      if(TASK(1:2) .eq. 'FG') then
                      !$$ --- to calculate potential
                         call CalForce_ForceClass(SimBox,   CtrlParam, ForceClass)
                         call CopyFPFrom_Devices_to_Host(FP)
                         call UpdateEpot_ForceClass(SimBox, CtrlParam, ForceClass)
                         call CopyEPOTFrom_Devices_to_Host(EPOT)

                         OBJF = 0.D0
                         do  IP=1, NPATH
                             IA1 = (IP-1)*NIMG*NPRT0+1
                             do IM=1, NIMG
                                IA = IA1 + (IM-1)*NPRT0
                                OBJF = OBJF + sum(EPOT(IA:IA+NPRT0-1)-EPOT(IA1:IA1+NPRT0-1)) !sum(m_EPOT(1:NPRT))
                             end do
                             !$$--- to add spring force and spring potential
                             if(HASSPRING .gt. 0) then
                               call NEBSpringForce_CPU(SimBox, CtrlParam, IP, EPOT, FP, SPOT, Meth, Climb)
                               OBJF = OBJF + SPOT
                             end if
                         end do

                         G(1:NIDX)                      = -FP(IDX(1:NIDX),1)
                         G(NIDX+1:      NIDX+NIDY)      = -FP(IDY(1:NIDY),2)
                         G(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = -FP(IDZ(1:NIDZ),3)

                         cycle
                      end if  !--- end if TASK = "FG"

                      if(TASK(1:5) .eq. 'NEW_X') then
                          cycle
                      end if

                   write(*,fmt="(A, I6)") " MDPSCU Message: LBFG performed for steps ", I
                   write(*,fmt="(A, I6)") "                 exit with code> "// TASK(1:len_trim(TASK))
                   exit
                end do
                deallocate(IDX, IDY, IDZ, NBD, IWA, X, L, U, G, WA,STAT=IERR)
                deallocate(EPOT, FP,STAT=IERR)

                if(I .gt. MXNUMSTEPS) then
                   IFLAG = 1
                   write(*,fmt="(A, I2, A, I7, A)") " MDPSCU Warning: LBFG exit with code after ", MXNUMSTEPS, " steps"
                   write(*,fmt="(A, I2)")           "                 more steps may be needed for better convergence"
                   call ONWARNING(gm_OnWarning)
                end if

                if(TASK(1:8) .eq. 'ABNORMAL' ) then
                   IFLAG = -1
                   write(*,fmt="(A, I2)") " MDPSCU Warning: LBFG exit with code "//TASK(1:len_trim(TASK))
                   write(*,fmt="(A, I2)") "                 optimizing may be fail"
                   call ONWARNING(gm_OnWarning)
                end if

      return
  end subroutine Do_NEBLBFGS_Forsteps_DEV
  !*********************************************************************************

  !****************************************************************************************
  subroutine Do_Damp(SimBox, CtrlParam, ForceClass)
   !***  PORPOSE: to damping the configuration to zero temperature.
   !
   !    INPUT:  CtrlParam,  the control parameters for simulation
   !            REINIT,     indicating if re-initialization of device memory is needed
   !    OUTPUT:  SimBox,    the box that has been damped
   !
   !
   use MD_LBFGSScheme_GPU,    only:Do_LBFGSB_Forsteps_DEV
   use MD_SteepestScheme_GPU, only:Do_Steepest_Forsteps_DEV
   use MD_CGScheme_GPU,       only:Do_CG_Forsteps_DEV
   use MD_ForceLib_Factory_GPU
   IMPLICIT NONE
    !---dummy vaiables
       type(SimMDBox),dimension(:)            ::SimBox
       type(SimMDCtrl),            intent(in) ::CtrlParam
       type(MDForceClassGPU),      intent(in) ::ForceClass
    !--- local
        integer::IFLAG
 
            select case(iand(CtrlParam%Quench_Meth, CP_LWORD))
                   case(CP_DAMPSCHEME_LBFGS)
                        !call Do_LBFGSB_Forsteps_DEV(SimBox,CtrlParam, gm_ForceClass, MXnumsteps, IFLAG)
                       call Do_NEBLBFGS_Forsteps_DEV (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, IFLAG, HASSPRING=0, METH=CtrlParam%NEB_Meth, CLIMB=0)
                   case(CP_DAMPSCHEME_ST )   
                        call Do_Steepest_Forsteps_DEV(SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
                   case(CP_DAMPSCHEME_CG )   
                        call Do_CG_Forsteps_DEV      (SimBox, CtrlParam, ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
             end select
       
   end subroutine Do_Damp
   !****************************************************************************************

  !*********************************************************************************
  subroutine CheckStableState(SimBox, CtrlParam, CURNUMPATH, NMINCFG, ICFG, IBOX, DRTOL) !, IFSD)
  !***  PURPOSE:   to check number of diffferent configurations contained in the
  !                SimBox array array
  !
  !
  !     INPUT:   SimBox,    the simulation box
  !              CtrlParam: the control parameters
  !              CURNUMPATH: current number of pathes in boxes

  !     NOTE:   Positions of atoms, m_XP in the original box should be copy out before calling this
  !             subroutine
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),dimension(:),       intent(in) ::SimBox
      type(SimMDCtrl),                   intent(in) ::CtrlParam
      integer,dimension(:),              intent(in) ::CURNUMPATH
      integer,                           intent(out)::NMINCFG
      integer,dimension(:)                          ::ICFG, IBOX
      real(KINDDF)                                  ::DRTOL
      !--- Local variables
           integer::IB, IPATH, NPRT, IP, IP0, IP1, IM, IM0, IM1, NIMG, IA, curMinCfg,IFPD(3),K, BFLG, CFLG
           real(KINDDF):: RMD, SEP(3), BS(3), HBS(3), E0, EC, DRTOL2

             NPRT      = SimBox(1)%NPRT
             BS        = SimBox(1)%ZL
             HBS       = C_HALF*BS
             IFPD      = CtrlParam%IFPD
             NIMG      = CtrlParam%NEB_NUMIMG
             DRTOL2    = DRTOL*DRTOL

             NMINCFG = 0
             ICFG    = 0
             IBOX    = 0
             IP0     = 1

             do IB=1, CtrlParam%MULTIBOX

                IP1 = IP0 + CURNUMPATH(IB) - 1
                do IP=IP0, IP1     !--- Loop for paths
                   IM0       = (IP-1)*NIMG + 1
                   IM1       = IM0 + NIMG -1
                   curMinCfg = IM0

                   !$$--- skip the first image on each path
                   if(IP .eq. IP0) then
                      NMINCFG       = NMINCFG + 1
                      ICFG(NMINCFG) = curMinCfg
                      IBOX(NMINCFG) = IB
                   end if
                   BFLG = 0
                   CFLG = 0
                   do IM = IM0+1, IM1
                      !if(.not.IFSD(IM)) then  !--- we skip the image of maxma
                      CFLG = CFLG+1
                      if(CFLG .gt. 1) then
                         do IA=1, NPRT
                            do k=1, 3
                               SEP(K) = SimBox(IM)%XP(IA,K) - SimBox(curMinCfg)%XP(IA,K)
                               if((IFPD(K).GT.0) .and. (dabs(SEP(K)) .GT. HBS(K)) ) then
                                   SEP(K) = SEP(K) - DSIGN(BS(K),SEP(K))
                               end if
                             end do
                             if(sum(SEP*SEP) .gt. DRTOL2) then
                                curMinCfg     = IM
                                NMINCFG       = NMINCFG + 1
                                ICFG(NMINCFG) = IM
                                IBOX(NMINCFG) = IB
                                BFLG          = 1
                                CFLG          = 0
                                exit
                            end if
                         end do
                      end if

                      !$$--- to set the last node as the original configure
                      if(IM .eq. IM1) then
                         if(BFLG .eq. 0) then
                            !--- the intermediate image collapse to the first image
                            !    we must have the last one to be the node
                            NMINCFG = NMINCFG + 1
                         end if
                         ICFG(NMINCFG) = IM
                      end if
                   end do
                end do
                IP0 = IP1 + 1
              end do

      return
  end subroutine CheckStableState
  !*********************************************************************************

  !*********************************************************************************
  subroutine Putout_NEBImages(SimBox, CtrlParam, Job, Path, Flag)
  !***  PURPOSE:   to initialize the NEB images
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !                Job:       the ID of TEST
  !                Path:      the node-states of boxes
  !                Flag:      an flag controling output
  !     OUTPUT
  !
  !
  !
   implicit none
      !--- dummy vaiables
      type(SimMDBox),dimension(:),   intent(inout):: SimBox
      type(SimMDCtrl),               intent(in)   :: CtrlParam
      integer,                       intent(in)   :: Job
      integer, dimension(:),         intent(in)   :: Path
      integer,                       intent(in)   :: Flag

      !--- Local variables
           character*256::GFile, GFile0, PFILE0, PFILE
           integer::NIMG, IMG, IMG0, IB, IB0, IP, IA, K, IFPD(3)
           real(KINDDF)::BS(3), HBS(3), SEP
           type(MDRecordStamp)::STAMP
           integer,      dimension(:), allocatable::IG
           save::IG

             if(Flag .eq. 0 ) then
                allocate(IG(CtrlParam%MULTIBOX))
                IG  = 0
             else if(Flag .lt. 0 ) then
                if(allocated(IG))  deallocate(IG)
                return
             end if

             NIMG     = CtrlParam%NEB_NUMIMG
             !--- to output the paths
             call GetPath(CtrlParam%f_geometry, GFILE0)
             GFILE0 = GFILE0(1:len_trim(GFILE0))//"NEBImg"
             write(*,fmt="(A,A)") " MDPSCU Message: save NEB images"

             !--- create stamp
             STAMP%AppType   = gm_AppType
             STAMP%ITest     = Job
             STAMP%ITime     = 0
             STAMP%ISect     = 1
             STAMP%Time      = 0
             STAMP%ScalTime  = 0

             BS   = SimBox(1)%ZL
             HBS  = C_HALF*BS
             IFPD = CtrlParam%IFPD

             IB0 =  1
             do IB=1, CtrlParam%MULTIBOX
                STAMP%IBox = IB+(Job-1)*CtrlParam%MULTIBOX

                do IP=1, Path(IB)
                   call STRCATI(GFILE, GFILE0, "_P",  0,  4)
                   call STRCATI(GFILE, GFILE,  "_",   STAMP%IBox(1), 4)

                   if(IP .gt. 1) then
                      SimBox(IB0)%XP1 = SimBox(IB0-1)%XP1
                      SimBox(IB0)%DIS = SimBox(IB0-1)%DIS
                      IMG0 = 2
                      IB0  = IB0 + 1
                   else
                      IMG0 = 1
                   end if

                   do IMG = IMG0, NIMG
                      IG(IB) = IG(IB) + 1

                      STAMP%ICfg = IG(IB)
                      STAMP%IRec = STAMP%ICfg
                      !--- to calcualte the displacement to move
                      if(IMG .lt. NIMG)then
                         SimBox(IB0)%XP1 = (SimBox(IB0+1)%XP - SimBox(IB0)%XP)

                      else !--- on IMG = NIMG
                         if(IP .lt. Path(IB))then
                            SimBox(IB0)%XP1 = (SimBox(IB0+2)%XP - SimBox(IB0)%XP)
                         else
                            SimBox(IB0)%XP1 = 0.D0
                         end if
                      end if

                      do IA=1, SimBox(IB0)%NPRT
                         do K=1, 3
                            if(IFPD(K)  .and. dabs(SimBox(IB0)%XP1(IA,K)) .gt. HBS(K) ) then
                               SimBox(IB0)%XP1(IA,K) = SimBox(IB0)%XP1(IA,K) - DSIGN(BS(K),SimBox(IB0)%XP1(IA,K))
                            end if
                         end do
                      end do
                      SimBox(IB0)%XP1  = SimBox(IB0)%XP1/SimBox(IB0)%RR

                      !--- to cummulate the displacement
                      if(IMG .eq. 1) then
                         SimBox(IB0)%DIS = 0.D0
                      else
                         do IA=1, SimBox(IB0)%NPRT
                            do K=1, 3
                               SEP =  SimBox(IB0)%XP(IA,K) - SimBox(IB0-1)%XP(IA,K)
                               if(IFPD(K)  .and. dabs(SEP) .gt. HBS(K) ) then
                                  SEP = SEP - DSIGN(BS(K),SEP)
                               end if
                               SimBox(IB0)%DIS(IA,K) = SimBox(IB0-1)%DIS(IA,K) + SEP
                            end do
                         end do
                      end if

                      call Putout_Instance_Config_SimMDBox(GFILE, SimBox(IB0), STAMP)
                      IB0 = IB0 + 1
                   end do

                end do

             end do

      return
  end subroutine Putout_NEBImages
  !*********************************************************************************

  !*********************************************************************************
  subroutine Putout_NEBEnergy(SimBox, CtrlParam, Job, Path, Flag, Statu)
  !***  PURPOSE:   to initialize the NEB images
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam: the control parameters
  !                Job:       the  ID of TEST
  !                Path:      the node-states of boxes
  !                Flag
  !                Statu:     the return value of optimization routines
  !                           indicating if convergence reached
  !     OUTPUT
  !
  !
  !
   implicit none
      !--- dummy vaiables
      type(SimMDBox), dimension(:), intent(inout):: SimBox
      type(SimMDCtrl),              intent(in)   :: CtrlParam
      integer,                      intent(in)   :: Job
      integer, dimension(:),        intent(in)   :: Path
      integer,                      intent(in)   :: Flag, Statu

      !--- Local variables
           character*256::PFILE
           integer::NIMG, IMG, IMG0, IB, IB0, IP, IG, I0, K, IA, IFPD(3)
           real(KINDDF)::PATHL, SEP(3), BS(3), HBS(3), E00, EPREV, ROT
           integer,      dimension(:),   allocatable:: hFILEP, PNIMG
           real(KINDDF), dimension(:),   allocatable:: E0, PPATHL
           real(KINDDF), dimension(:,:), allocatable:: VECTP, VECT
           save::hFILEP, PNIMG, E0, PPATHL


             !--- if Flag = 0, to open the IO units
             if(Flag .eq. 0) then
               allocate(hFILEP(CtrlParam%MULTIBOX), PNIMG(CtrlParam%MULTIBOX),E0(CtrlParam%MULTIBOX), PPATHL(CtrlParam%MULTIBOX))

               do IB=1, CtrlParam%MULTIBOX
                  call GetPath(CtrlParam%f_geometry, PFILE)
                  PFILE = PFILE(1:len_trim(PFILE))//"NEBPot"
                  call STRCATI(PFILE, PFILE, "_P",  0,  4)
                  call STRCATI(PFILE, PFILE,  "_",  IB+(Job-1)*CtrlParam%MULTIBOX, 4)
                   call AvailableIOUnit(hFileP(IB))
                  open(UNIT=hFileP(IB), FILE=PFILE)
                  write(hFileP(IB), fmt="(A)") "!--- DIFFUSION PATH CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
                  write(hFileP(IB), fmt="(A, I5)") '!--- The box #:  ', IB+(Job-1)*CtrlParam%MULTIBOX
                  write(hFilep(IB), fmt="(2A8, 7(1x, A19))") "T-State#", "State#", "G-Path len(latt)",  &
                                                                                   "Path len(latt)  ",  &
                                                                                   "Potential(eV)   ",  &
                                                                                   "G-Potential(eV)  ", &
                                                                                   "L-Potential(eV)  ", &
                                                                                   "Cos(ang)         ", &
                                                                                   "Convergence      "

                  write(*,fmt="(A,A)") " MDPSCU Message: Potentials saved to: ",PFILE(1:len_trim(PFILE))

                  PNIMG(IB)  = 0
                  PPATHL(IB) = 0.D0
                  E0(IB)     = 0.D0
               end do
             else if(Flag .lt. 0) then
                  if(allocated(hFILEP)) then
                     do IB=1, CtrlParam%MULTIBOX
                        close(hFILEP(IB))
                     end do
                     deallocate(hFILEP)
                  end if
                  if(allocated(PNIMG))  deallocate(PNIMG)
                  if(allocated(E0))     deallocate(E0)
                  if(allocated(PPATHL)) deallocate(PPATHL)
                  return
             end if

             write(*,fmt="(A,A)") " MDPSCU Message: save potentials of images"
             NIMG = CtrlParam%NEB_NUMIMG
             BS   = SimBox(1)%ZL
             HBS  = C_HALF*BS
             IFPD = CtrlParam%IFPD
             allocate(VECTP(SimBox(1)%NPRT,3),VECT(SimBox(1)%NPRT,3))

             IB0  =  1
             do IB=1, CtrlParam%MULTIBOX
                IG    = 0
                I0    = IB0
                PATHL = 0.D0
                do IP=1, Path(IB)
                   if(IP .gt. 1) then
                      IMG0  = 2
                      IB0   = IB0 + 1
                      EPREV = E00
                   else
                      IMG0  = 1
                      EPREV = 0
                   end if

                   do IMG = IMG0, NIMG
                      IG = IG + 1
                      if(IMG .gt. 1) then
                         do IA=1, SimBox(IB0)%NPRT
                            do K=1,3
                               SEP(K) = SimBox(IB0)%XP(IA,K) - SimBox(IB0-1)%XP(IA,K)
                               if((IFPD(K).GT.0) .and. dabs(SEP(K)) .GT. HBS(K)) then
                                   SEP(K) = SEP(K) - DSIGN(BS(K),SEP(K))
                               end if
                            end do
                            PATHL = PATHL + dsqrt(sum(SEP*SEP))
                            VECT(IA,1:3) = SEP(1:3)
                         end do
                         VECT = VECT/dsqrt(sum(VECT*VECT))
                      end if
                      E00 = sum(SimBox(IB0)%EPOT-SimBox(I0)%EPOT)
                      if(IMG .gt. 2 .or. IP.gt.1) then
                         ROT = sum(VECT*VECTP)
                      else
                         ROT = 1
                      endif

                      write(hFILEP(IB), fmt="(2(I7,1x), 6(1x,1PE19.7),(1x,I4))")                &
                                        PNIMG(IB)+IG, IG,                                       &
                                       (PATHL+PPATHL(IB))/SimBox(IB0)%RR, PATHL/SimBox(IB0)%RR, &
                                        sum(SimBox(IB0)%EPOT)*CP_ERGEV/dble(SimBox(IB0)%NPRT),  &
                                       (E00+E0(IB))*CP_ERGEV, (E00-EPREV)*CP_ERGEV,             &
                                        ROT,                                                    &
                                        Statu

                      VECTP = VECT
                      IB0 = IB0 + 1
                   end do
                end do
                PNIMG(IB)  = PNIMG(IB)  + IG
                PPATHL(IB) = PPATHL(IB) + PATHL
                E0(IB)     = E0(IB)     + E00
             end do
             deallocate(VECTP,VECT)
      return
  end subroutine Putout_NEBEnergy
  !*********************************************************************************

  end module MD_Method_NEB_GPU
