 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  This program is used to search possible events that may thermal-dynamically occur in a
 !                  many-body system. The EventSearch_module provides a number of interfaces to
 !                  MD_SimBoxArray_AppShell_16_GPU.F90. The program should be run in GMD method.
 !
 !
 !                  SEE ALSO____________________________________________________________________________
 !                       MD_Method_GenericMD_GPU.F90
 !                       MD_SimBoxArray_ToolShell_16_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  by HOU Qing, Dec., 2016
 !
  module EventSearch_module
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl

  implicit none

         integer::m_NREC=0
         integer, dimension(:),        allocatable, private:: m_EVENTID
         integer, dimension(:),        allocatable, private:: m_CBFLAG, m_CAFLAG, m_MASK
         type(SimMDBox), dimension(:), allocatable, private:: m_SimBoxBef
         type(SimMDBox), dimension(:), allocatable, private:: m_SimBoxNow
  contains

  !**********************************************************************************
   subroutine Initialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize
   !
   implicit none
  !----   DUMMY Variables
   type(SimMDBox)  :: SimBox
   type(SimMDCtrl) :: CtrlParam

  !----   Local variables
   integer::ERR

         !--- to check if the input is in GMD format
         if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") then
            call Clear(SimBox, CtrlParam)
            return
         end if
         allocate(m_EVENTID(CtrlParam%MultiBox), m_CBFLAG(CtrlParam%MultiBox), m_SimBoxBef(CtrlParam%MultiBox), m_SimBoxNow(CtrlParam%MultiBox), STAT=ERR)
         if(ERR .gt. 0) then
            write(*,*) "MDPSCU Error: in allocate memory in Initializa of EventSearch_module"
            stop
         end if   
         m_EVENTID     = 0
         m_NREC        = 0
         return
   end subroutine Initialize
 !**********************************************************************************

 !**********************************************************************************
   subroutine Clear(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize
   !
   use MD_SimboxArray
   implicit none
  !----   DUMMY Variables
   type(SimMDBox)  :: SimBox
   type(SimMDCtrl) :: CtrlParam

  !----   Local variables

         !--- to check if the input is in GMD format
         if(allocated(m_EVENTID))   deallocate(m_EVENTID)
         if(allocated(m_CBFLAG))    deallocate(m_CBFLAG)
         if(allocated(m_CAFLAG))    deallocate(m_CAFLAG)
         if(allocated(m_MASK))      deallocate(m_MASK)

         if(allocated(m_SimBoxBef)) then
            call Release_SimBoxArray(m_SimBoxBef)
            deallocate(m_SimBoxBef)
         endif
         if(allocated(m_SimBoxNow)) then
            call Release_SimBoxArray(m_SimBoxNow)
            deallocate(m_SimBoxNow)
         endif
         m_NREC    = 0

         return
   end subroutine Clear
!**********************************************************************************

!**********************************************************************************
   subroutine BeforeOneTest(ITEST, SimBox0, SimBox, CtrlParam)
   !***  PURPOSE:  to reset the EVENT count
   !
   implicit none
  !----   DUMMY Variables
   integer                    :: ITEST
   type(SimMDBox)             :: SimBox0
   type(SimMDBox),dimension(:):: SimBox
   type(SimMDCtrl)            :: CtrlParam
  !----   Local variables

         !--- to check if the input is in GMD format
         if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") return

         m_EVENTID =  0
         m_NREC    = 0
         return
   end subroutine BeforeOneTest
 !**********************************************************************************

 !**********************************************************************************
   subroutine DoEventCheck(Stamp, SimBox, CtrlParam)
   use MD_LBFGSScheme_GPU,      only:DO_LBFGSB_FORSTEPS_DEV
   use MD_ForceLib_Factory_GPU, only:gm_ForceClass, CalEpot_ForceClass
   use MD_SimboxArray
   use MD_SimBoxArray_GPU
   use MD_NeighborsList_GPU
   implicit none
  !----   DUMMY Variables
    type(MDRecordStamp)         :: Stamp
    type(SimMDBox), dimension(:):: SimBox
    type(SimMDCtrl)             :: CtrlParam
  !----   Local variables
    integer::I, NB, hFILE, ERR
    type(MDRecordStamp):: LSTAMP
    character*256::GFile, CFile, GFile0, GFileT
   
      !---------------
         if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") return
         if(Stamp%ITime .le. 0) return

         write(*,fmt="(A,I4, A, I8)") ' MDPSCU Message: to do event check for TEST: ',Stamp%ITest, ' at ITIME ', Stamp%ITime
         write(*,fmt="(A,I4, A, I8)") '                 current record number: ',m_NREC
         NB = size(SimBox)
         call Copy_RecordStamp(Stamp, LSTAMP)
         LSTAMP%NBox = 1
         call STRCATI(GFILE0, CtrlParam%f_geometry, "_0K_P", Stamp%IProcess, 4)
         call STRCATI(GFILET, CtrlParam%f_geometry, "_TK_P", Stamp%IProcess, 4)

         !*** for the first time, we create the initial configuration
         if(m_NREC .eq. 0) then
            do I=1, NB
               call Copy_SimMDBox(SimBox(I), m_SimBoxBef(I))
               call Copy_SimMDBox(SimBox(I), m_SimBoxNow(I))
            end do
            !*** we need to quench the system
            !call DO_LBFGSB_FORSTEPS_DEV(SimBox, CtrlParam, gm_ForceClass, MXITER, IFLAG)
            call Do_Damp(m_SimBoxBef, CtrlParam) 
            
            call CalEpot_ForceClass(m_SimBoxBef,  CtrlParam, gm_ForceClass)
            call CopyOut_SimBox_DEV(m_SimBoxBef)

            !*** the statu on GPUs have been changed, we need to restore them
            call CopyIn_SimBox_DEV(SimBox)
            call Cal_NeighBoreList_DEV(SimBox, CtrlParam)
          
            !--- save the first configurations
            do I=1, NB
               LSTAMP%IBox = (LSTAMP%ITest-1)*NB+I
               LSTAMP%ICfg = m_EVENTID(I)
               call STRCATI(GFILE, GFILE0, "_", LSTAMP%IBox(1), 4)
               call Putout_Instance_Config_SimMDBox(GFILE, m_SimBoxBef(I), LSTAMP)
               !call STRCATI(GFILE, GFILET, "_", LSTAMP%IBox(1), 4)
               !call Putout_Instance_Config_SimMDBox(GFILE, SimBox(I), LSTAMP)
            end do

            !---
            ERR = 0
            if(.not.allocated(m_MASK)) then
               allocate(m_MASK(m_SimBoxBef(1)%NPRT), STAT=ERR )
               if(ERR) goto 1000
            end if   

            if(.not.allocated(m_CAFLAG)) then
               allocate(m_CAFLAG(size(m_SimBoxBef)*m_SimBoxBef(1)%NPRT), STAT=ERR )
               if(ERR) goto 1000
            end if    
            m_NREC = m_NREC + 1
            return
         end if

         !**** For already having the initial boxes.
         !---  Before to detetc event, we quench the current statu
         !---  NOTE: in Appshell, the boxes have been copyout from GPU. We do not 
         !           do copyout again
            !m_SimBoxNow =  SimBox
            call Copy_SimBoxArray(SimBox, m_SimBoxNow)
            call Do_Damp(m_SimBoxNow, CtrlParam)
            call CalEpot_ForceClass(m_SimBoxNow,  CtrlParam, gm_ForceClass)
            call CopyOut_SimBox_DEV(m_SimBoxNow)

            !*** do event detecting
            m_CBFLAG  = 0
            m_CAFLAG  = 0
            m_Mask    = 1
            call DoCompare(m_SimBoxBef, m_SimBoxNow, CtrlParam, m_CBFLAG, m_CAFLAG, m_Mask)
            do I=1, NB
               if(m_CBFLAG(I) .gt. 0) then
                  m_EVENTID(I) = m_EVENTID(I) + 1
                  LSTAMP%IBox = (LSTAMP%ITest-1)*NB+I
                  LSTAMP%ICfg  = m_EVENTID(I)

                  write(*,fmt="(A,I5, A, I7)") ' MDPSCU Message: event # ',m_EVENTID(I), ' found in box # ', LSTAMP%IBox(1)
         
                  call STRCATI(GFILE, GFILE0, "_", LSTAMP%IBox(1), 4)
                  call Putout_Instance_Config_SimMDBox(GFILE, m_SimBoxNow(I), LSTAMP)

                  !call STRCATI(GFILE, GFILET, "_", LSTAMP%IBox(1), 4)
                  !call Putout_Instance_Config_SimMDBox(GFILE, SimBox(I), LSTAMP)
                  call Copy_SimMDBox(m_SimBoxNow(I), m_SimBoxBef(I))
               end if
            end do
            m_NREC = m_NREC + 1

            !*** the statu on GPUs have been changed, we need to restore them
            !write(*,fmt="(A,I5, A, I7)") ' MDPSCU Message: Copy boxes back'
            call  CopyIn_SimBox_DEV(SimBox)
            !write(*,fmt="(A,I5, A, I7)") ' MDPSCU Message: Restore neighbor-list'
            call  Cal_NeighBoreList_DEV(SimBox, CtrlParam)
            return

   1000     write(*, *) "MDPSCU Error: in allocate memory in DoEventCheck of EventSearch_module"             
   end subroutine DoEventCheck
 !*********************************************************************************

 !**********************************************************************************
   subroutine DoArchive(Stamp, SimBox, CtrlParam)
   use MD_SimboxArray
   implicit none
  !----   DUMMY Variables
    type(MDRecordStamp)         :: Stamp
    type(SimMDBox), dimension(:):: SimBox
    type(SimMDCtrl)             :: CtrlParam
  !----   Local variables
    integer::hFILE
    character*256::CFile

         if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") return
         if(Stamp%ITime .le. 0) return

            !*** store EVENTID for each box for restart
            call AvailableIOUnit(hFILE)
            call GetPath(CtrlParam%f_configure,CFile)
            CFile = CFile(1:len_trim(CFile))//gm_ExeName(1:len_trim(gm_ExeName))
            call STRCATI(CFILE, CFile, "P", Stamp%IProcess, 4)
            call STRCATI(CFILE, CFILE, "_", Stamp%ITest, 4)
            open(hFile, file = CFILE, form='unformatted', status='unknown')
              call Archive_RecordStamp(hFile, Stamp)
              write(hFile)m_NREC, m_EVENTID
            close(hFile)

            return
   end subroutine DoArchive
 !*********************************************************************************

 !**********************************************************************************
   subroutine DoRestore(Stamp, SimBox, CtrlParam)
   use MD_SimboxArray
   implicit none
  !----   DUMMY Variables
    type(MDRecordStamp)         :: Stamp
    type(SimMDBox), dimension(:):: SimBox
    type(SimMDCtrl)             :: CtrlParam
  !----   Local variables
    integer::hFILE, I, NB, ERRMSG
    type(MDRecordStamp)::LSTAMP
    character*256::GFile, CFile, GFile0

         if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") return
         if(Stamp%ITime .le. 0) return

            !*** store EVENTID for each box for restart
            call AvailableIOUnit(hFILE)
            call GetPath(CtrlParam%f_configure,CFile)
            CFile = CFile(1:len_trim(CFile))//gm_ExeName(1:len_trim(gm_ExeName))
            call STRCATI(CFILE, CFile, "P", Stamp%IProcess, 4)
            call STRCATI(CFILE, CFILE, "_", Stamp%ITest, 4)
            open(hFile, file = CFILE, form='unformatted', status='old')
               call Restore_RecordStamp(hFile, LSTAMP)
               read(hFile)m_NREC, m_EVENTID
            close(hFile)

            !*** restore the "before" configuration
            NB = size(m_SimBoxBef)
            call STRCATI(GFILE0, CtrlParam%f_geometry, "_0K_P", Stamp%IProcess, 4)
            do I=1, NB
               call Copy_SimMDBox(SimBox(I), m_SimBoxBef(I))
               call Copy_SimMDBox(SimBox(I), m_SimBoxNow(I))
               call Default_RecordStamp(LSTAMP)
               LSTAMP%ITest = Stamp%ITest
               LSTAMP%IBox  = (LSTAMP%ITest-1)*NB+I
               LSTAMP%ICfg  = m_EVENTID(I)
               call STRCATI(GFILE, GFILE0, "_", (LSTAMP%ITest-1)*NB+I, 4)
               call STRCATI(GFILE, GFILE,  ".", m_EVENTID(I), 4)
               call Putin_Instance_Config_SimBoxArray(GFILE, m_SimBoxBef(I:I), LSTAMP, ERR=ERRMSG)
            end do
            if(.not.allocated(m_MASK)) &
               allocate(m_MASK(m_SimBoxBef(1)%NPRT) )
            if(.not.allocated(m_CAFLAG)) &
               allocate(m_CAFLAG(size(m_SimBoxBef)*m_SimBoxBef(1)%NPRT) )

            return
   end subroutine DoRestore
 !****************************************************************************************

  !****************************************************************************************
  subroutine Do_Damp(SimBox, CtrlParam)
   !***  PORPOSE: to damping the configuration to zero temperature.
   !
   !    INPUT:  CtrlParam,  the control parameters for simulation
   !            REINIT,     indicating if re-initialization of device memory is needed
   !    OUTPUT:  SimBox,    the box that has been damped
   !
   !
    use MD_ForceLib_Factory_GPU
    use MD_LBFGSScheme_GPU,    only:Do_LBFGSB_Forsteps_DEV
    use MD_SteepestScheme_GPU, only:Do_Steepest_Forsteps_DEV
    use MD_CGScheme_GPU,       only:Do_CG_Forsteps_DEV
    use MD_DiffScheme_GPU,     only:Do_DynDamp_Forsteps_DEV
 
    implicit none
    !---dummy vaiables
        type(SimMDBox), dimension(:)::SimBox
        type(SimMDCtrl)             ::CtrlParam
       !Local variables
        integer::IFLAG
     
            select case(iand(CtrlParam%Quench_Meth, CP_LWORD))
                   case( CP_DAMPSCHEME_LBFGS)
                      call DO_LBFGSB_FORSTEPS_DEV  (SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps, IFLAG)
                   case(CP_DAMPSCHEME_DYN)
                      call Do_DynDamp_Forsteps_DEV (SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps)
                   case(CP_DAMPSCHEME_ST )   
                      call Do_Steepest_Forsteps_DEV(SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
                   case(CP_DAMPSCHEME_CG )   
                      call Do_CG_Forsteps_DEV      (SimBox, CtrlParam, gm_ForceClass, CtrlParam%Quench_Steps, CtrlParam%Quench_Meth)
             end select
  
   end subroutine Do_Damp
  !****************************************************************************************   

 !****************************************************************************************
  subroutine DoCompare(SimBoxBef, SimBoxNow, CtrlParam, BFlag, AFlag, Mask)
  !***  PORPOSE: to compare two boxarray
  !
  !    INPUT:  SimBoxBef,    the boxes to compared to
  !            SimBoxNow,    the current boxes
  !
  !    OUTPUT: BFLAG,         the flag indicating which boxes have changed
  !            AFLAG,         the flag indicating which atoms have changed
  !

  implicit none
   !---dummy vaiables
       type(SimMDBox), dimension(:)             :: SimBoxBef, SimBoxNow
       type(SimMDCtrl),              intent(in) :: CtrlParam
       integer,        dimension(:), intent(out):: BFlag, AFlag
       integer,        dimension(:), intent(in) :: Mask

      !Local variables
      integer::IB, NB, NPRT, I, IFPD(3), K, IP
      real(KINDDF)::RC2, SEP(3), BOX(3), HBOX(3), DRTOL


           !*** to prepare the swapboxes
           NB    = size(SimBoxBef)
           NPRT  = SimBoxBef(1)%NPRT
           BOX   = SimBoxBef(1)%ZL
           HBOX  = SimBoxBef(1)%ZL*C_HALF
           IFPD  = CtrlParam%IFPD

           RC2  =  CtrlParam%STRCUT_DRTol*CtrlParam%STRCUT_DRTol
           IP = 0
           do IB=1, NB
              do I=1, NPRT
                 IP = IP + 1
                 AFlag(IP) = 0
                 if(MASK(I) .le. 0 ) cycle

                 do K=1, 3
                    SEP(K)  =  SimBoxNow(IB)%XP(I,K) - SimBoxBef(IB)%XP(I,K)
                    SimBoxNow(IB)%XP1(I,K) = SimBoxNow(IB)%XP(I,K) + SEP(K)
                    if( (IFPD(K).GT.0) .AND. (DABS(SEP(K)) .GT. HBOX(K))) then
                         SEP(K) = SEP(K) - DSIGN(BOX(K),SEP(K))
                    end if
                 end do

                 if(sum(SEP*SEP) .gt. RC2) then
                    AFlag(IP) = 1
                    BFlag(IB) = 1
                 end if

              end do
           end do
     return
  end subroutine DoCompare
 !****************************************************************************************

 !****************************************************************************************
 end module EventSearch_module


 !****************************************************************************************
  Program EventSearch_GPU_Main
  use MD_SimBoxArray_AppShell_16_GPU, only:APPSHELL_AddRecord, APPSHELL_Main
  !---- If user-define potential to be used, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
  !     use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
  !     use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table

  use EventSearch_module
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord(PRERECORD=Initialize,      &
                               RECORDPROC=DoEventCheck,   &
                               BEFONETEST=BeforeOneTest,  &
                               AFTRECORD=Clear,           &
                               RESTARTREC=DoRestore,      &
                               SAVERECSTATU=DoArchive)
       call APPSHELL_Main(numprocs,procid)

       stop
  end  program EventSearch_GPU_Main
