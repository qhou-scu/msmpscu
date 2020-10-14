 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to extract trajectory of atoms of certain type generated from MD 
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  Trajectory.F90.exe -I "filename" -T(est) t1, t2, ts -C(fg) c1, c2, cs -B(ox) b1, b2, b3
 !
 !                  OPTIONS
 !                        -T(est) t1, t2, ts:   - the range for the tests to be involved in analysis
 !                        -C(est) c1, c2, cs:   - the range for the configurations to be involved
 !                        -B(ox)  b1, b2, b3 : indenetify the boxes IDs in a TEST
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  1st version   2019-01 (Hou Qing Sichuan university)
 !
 !

 module Trajectory
 use MD_SimboxArray
 use MD_TYPEDEF_SimMDCtrl
 use MD_NeighborsList_GPU

 implicit none
 
      !--- NOTE the variable DefTyp subjected to changes
       integer::m_Nearest =8

       character*256                    ::m_InFname  = ""
       type(MDRecordStamp)              ::m_CurStamp
       type(MDDEVNeighbor_List)         ::m_dNNList
       integer                          ::m_TOTREC = 0
       integer                          ::m_NDEF   = 0
       integer                          ::m_DEFFLG(mp_MXGROUP) = 0
       integer                          ::m_DEFTYP(mp_MXGROUP) = 0
       integer, dimension(:),allocatable::m_DefAtomID
       type(StatementList),pointer      ::m_ExInfor=>null()

       type(SimMDBox),      dimension(:), allocatable ::m_MergedSimBox
       type(Neighbor_List), dimension(:), allocatable ::m_NNList
  contains
  
   !**********************************************************************************
   subroutine MyPreRecord(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation 
   implicit none
  !----   DUMMY Variables
   type(SimMDBox)  ::SimBox
   type(SimMDCtrl) ::CtrlParam

  !----   Local variables
  character*256 path
  
            if(len_trim(gm_cfgFileName(1)) .gt. 0) then
               m_InFname = gm_cfgFileName(1)
               m_TOTREC  = (CtrlParam%ENDCFG - CtrlParam%STARTCFG)/CtrlParam%CFGSTEP + 1
            else   
               m_InFname = CtrlParam%f_geometry
               call NumberRec_SimMDCtrl(CtrlParam, m_TOTREC)
            end if   
            
            m_CurStamp%ITime = - 1
            m_CurStamp%ITest = 0 
            
         return
   end subroutine MyPreRecord
  !********************************************************************************** 
   
  !**********************************************************************************
  subroutine RecordTrack(SimBox, CtrlParam, SimBoxTrack, NDef, DefAtomID, NNb, List, Time)
    !***  DESCRIPTION: to calcualte and putout the occupation state of lattice and interstials
    !                  identified using method#1
    !
    !    INPUT:  Stamp,       the recoding stamp
    !            SimBox,      the simulation boxs
    !            CtrlParam,   the control parameters for simulation
    !
     implicit none
           !--- dummy variables
           type(SimMDBox),        intent(in)    ::SimBox
           type(SimMDCtrl),       intent(in)    ::CtrlParam
           type(SimMDBox),        intent(inout) ::SimBoxTrack
           integer,               intent(in)    ::NDef, NNb 
           integer, dimension(:), intent(in)    ::DefAtomID
           type(NEIGHBOR_LIST),   intent(in)    ::List 
           real(KINDDF),          intent(in)    ::Time                              
           !--- local variables
           integer::I0, I, J, NN, IND, IW, K
           real(KINDDF)::VECT1(3), BS(3), HBS(3), LATT
           integer::IP, IP0, NPRT0, IDN
    
              LATT  = SimBox%RR
              BS    = SimBox%ZL
              HBS   = 0.5D0*BS
              NPRT0 = SimBoxTrack%NPRT

              call AddAtoms_SimMDBox(SimBoxTrack, (/NDef,NNb/),  (/2,1/))
              IP0 = max(0, NPRT0-(NDef+NNb))
              IP  = NPRT0
              do I0=1,  NDef   !~SimBox%NPRT
                 I = DefAtomID(I0)
                 !if(SimBox%ITYP(I) .eq. DType) then
                     IP0 = IP0 + 1
                     IP  = IP  + 1	
                     SimBoxTrack%XP(IP,:)  = SimBox%XP(I,:)
                     SimBoxTrack%ITYP(IP)  = SimBox%ITYP(I)
                     SimBoxTrack%XP1(IP,:) = SimBox%XP1(I,:)
                     SimBoxTrack%FP(IP,:)  = SimBox%FP(I,:)
                     SimBoxTrack%DIS(IP,:) = SimBox%DIS(I,:)
                     SimBoxTrack%EKIN(IP)  = SimBox%EKIN(I)
                     SimBoxTrack%EPOT(IP)  = SimBox%EPOT(I)
                     call SetVal_DataPad("IDN", SimBoxTrack%ptrDatPad, IP, I)
                     !--- to calculate the displcement vector in the time interval of record
                     do K=1,3
                        VECT1(K) = SimBoxTrack%XP(IP,K) - SimBoxTrack%XP(IP0,K)
                        if(dabs(VECT1(K)) .GT. HBS(K)) then
                           VECT1(K) = VECT1(K) - DSIGN(BS(K),VECT1(K))
                       end if
                     end do
                     VECT1(:) = VECT1(:)/LATT
                     call SetVal_DataPad("DVEC",SimBoxTrack%ptrDatPad, IP0, VECT1)
                     !call SetVal_DataPad("TIME",SimBoxTrack%ptrDatPad, IP, Time)

                     !--- to extract the neighbors
                     NN = min(List%Kvois(I),NNb)
                     do J=1, NN
                        IP0   = IP0 + 1
                        IP    = IP  + 1
                        IDN   = List%Indi(I, J)
                        SimBoxTrack%ITYP(IP)  = SimBox%ITYP(IDN)
                        SimBoxTrack%XP(IP,:)  = SimBox%XP(IDN,:)
                        SimBoxTrack%XP1(IP,:) = SimBox%XP1(IDN,:)
                        SimBoxTrack%FP(IP,:)  = SimBox%FP(IDN,:)
                        SimBoxTrack%DIS(IP,:) = SimBox%DIS(IDN,:)
                        SimBoxTrack%EKIN(IP)  = SimBox%EKIN(IDN)
                        SimBoxTrack%EPOT(IP)  = SimBox%EPOT(IDN)
                        call SetVal_DataPad("IDN",SimBoxTrack%ptrDatPad, IP, IDN)
                        
                        do K=1,3
                           VECT1(K) = SimBoxTrack%XP(IP,K) - SimBoxTrack%XP(IP0,K)
                           if(dabs(VECT1(K)) .GT. HBS(K)) then
                              VECT1(K) = VECT1(K) - DSIGN(BS(K),VECT1(K))
                          end if
                        end do
                        VECT1(:) = VECT1(:)/LATT
                        call SetVal_DataPad("DVEC",SimBoxTrack%ptrDatPad, IP0, VECT1)
                        !call SetVal_DataPad("TIME",SimBoxTrack%ptrDatPad, IP, Time)
                     end do  

                     do J=NN+1, NNb
                       IP    = IP + 1
                       SimBoxTrack%ITYP(IP)  = SimBox%NGROUP+1
                       SimBoxTrack%XP(IP,:)  = 0
                       SimBoxTrack%XP1(IP,:) = 0
                       SimBoxTrack%FP(IP,:)  = 0
                       SimBoxTrack%DIS(IP,:) = 0
                       SimBoxTrack%EKIN(IP)  = 0
                       SimBoxTrack%EPOT(IP)  = 0
                       call SetVal_DataPad("IDN",SimBoxTrack%ptrDatPad,  IP, 0)
                     end do
                 !end if	  
              end do
              
              return              
    end subroutine RecordTrack
  !***************************************************************************   

  !**********************************************************************************
  subroutine MyRecord(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calcualte and putout the occupation state of lattice and interstials
  !                  identified using method#1
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
   
   implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in)::Stamp
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables
       integer::IB, I, K, IP, ITY
       integer::RECN=0, NB=0
       save RECN, NB
       character*256 ::GFILE, STR

       !***************
          if(Stamp%ITIME .le. 0) return

          if(m_CurStamp%ITest    .ne. Stamp%ITest .or. &
             any(m_CurStamp%IBox .ne. Stamp%Ibox) ) then
             NB = size(SimBox)	
             call Release_SimBoxArray(m_MergedSimBox)	
             allocate(m_MergedSimBox(NB), m_NNList(NB))

             call Release_StatementList(m_ExInfor)
             
             !--- determine ID of defects
             m_DEFFLG(2) = 1
             m_DEFFLG(3) = 1
             m_NDEF = 0 !    = SimBox(1)%NA(m_DefTyp) !count(SimBox(1)%ITYP .eq. m_DefTyp)
             IP     = 0
             do I=1, size(m_DEFFLG) 
                if(m_DEFFLG(I).gt. 0)then
                  m_NDEF = m_NDEF + SimBox(1)%NA(I)
                  IP = IP + 1
                  m_DEFTYP(IP) = I
                end if   
             end do  
             if(allocated(m_DefAtomID) ) deallocate(m_DefAtomID)
             allocate(m_DefAtomID(m_NDEF))
             IP = 0
             do I=1, SimBox(1)%NPRT
                ITY = SimBox(1)%ITYP(I)
                if(m_DEFFLG(ITY) .gt. 0) then
                  IP              = IP + 1
                  m_DefAtomID(IP) = I
                end if
             end do   

             !--- begin recording the track
             m_CurStamp = Stamp
             RECN       = 1
             !call Reorder_NeighBoreList_Nearest_Dev(SimBox, CtrlParam%NB_MXNBS, m_NNList)
             call Reorder_NeighBoreList_Nearest_Dev(m_Nearest, m_dNNList, m_NNList)
             do IB = 1, NB
                call Initialize_SimMDBox((m_Nearest+m_NDEF)*m_TOTREC, m_MergedSimBox(IB))
                call NewDatPad_SimMDBox( m_MergedSimBox(IB), "IDN",   "I")
                call NewDatPad_SimMDBox( m_MergedSimBox(IB), "DVEC",   "D", 3)
                !call NewDatPad_SimMDBox( m_MergedSimBox(IB), "TIME",   "D")
                call CopyInformation_SimMDBox(SimBox(IB), m_MergedSimBox(IB))
                call RecordTrack(SimBox(IB), CtrlParam, m_MergedSimBox(IB), m_NDEF, m_DefAtomID, m_Nearest, m_NNList(IB), Stamp%TIME)
             end do
             
             write(STR,*) "&TOTREC  ", m_TOTREC
             call Add_StatementList(m_ExInfor, STR)
             write(STR,*) "&NDEF    ",  m_NDEF, m_DEFTYP(1:m_NDEF)
             call Add_StatementList(m_ExInfor, STR)
             write(STR,*) "&NEAREST ",  m_Nearest
             call Add_StatementList(m_ExInfor, STR)
             write(STR,*) "&EXTRECORD",  CtrlParam%TimestpR
             call Add_StatementList(m_ExInfor, STR)
             write(STR,*) "&STEPSIZE",  CtrlParam%HMI*CP_S2PS
             call Add_StatementList(m_ExInfor, STR)
             return
          end if
          !call Reorder_NeighBoreList_Nearest_Dev(SimBox, CtrlParam%NB_MXNBS, m_NNList)
          RECN = RECN + 1 
          call Reorder_NeighBoreList_Nearest_Dev(m_Nearest, m_dNNList, m_NNList)
          do IB = 1, NB
             call RecordTrack(SimBox(IB), CtrlParam, m_MergedSimBox(IB), m_NDEF, m_DefAtomID, m_Nearest, m_NNList(IB), Stamp%TIME)
          end do   
          
          if(RECN .eq.  m_TOTREC) 	then
             do IB=1, NB
                call STRCATI(GFILE, m_InFname, "_Track_P", 0, 4)
                K = (m_CurStamp%ITest - 1)*NB + IB
                call STRCATI(GFILE, GFILE, "_", K, 4)
                m_CurStamp%ICfg = -1
                call Putout_Instance_Config_SimMDBox(GFILE, m_MergedSimBox(IB), m_CurStamp, ExStats=m_ExInfor)
             end do    
          end if	
          
 end subroutine MyRecord
 !***************************************************************************
 
 !**********************************************************************************
   subroutine MyAftOneTest(SimBox0, Stamp, SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation 
   implicit none
  !----   DUMMY Variables
   type(SimMDBox)               :: SimBox0
   type(MDRecordStamp)          :: Stamp
   type(SimMDBox), dimension(:) :: SimBox
   type(SimMDCtrl)              ::CtrlParam

  !----   Local variables
   character*256 ::GFILE
   integer::IB, K, NB

              NB = size(m_MergedSimBox) 
              do IB = 1, NB
                 call STRCATI(GFILE, m_InFname, "_Track_P", 0, 4)
                 K = (m_CurStamp%ITest - 1)*NB + IB
                 call STRCATI(GFILE, GFILE, "_", K, 4)
                 m_CurStamp%ICFG = -1
                 call Putout_Instance_Config_SimMDBox(GFILE, m_MergedSimBox(IB), m_CurStamp, ExStats=m_ExInfor)
              end do
              call Release_SimBoxArray(m_MergedSimBox)
              call Release_StatementList(m_ExInfor)

              call Clear_NeighboreList_Dev(m_dNNList)
              if(allocated(m_NNList)) then
                 do IB=1, NB
                    call Clear_NeighboreList(m_NNList(IB))
                 end do
                 deallocate(m_NNList)
               end if   
              return
   end subroutine MyAftOneTest
 !********************************************************************************** *

 !**********************************************************************************
   subroutine MyAftRecord(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation 
   implicit none
  !----   DUMMY Variables
   type(SimMDBox)  ::SimBox
   type(SimMDCtrl) ::CtrlParam

  !----   Local variables

         return
   end subroutine MyAftRecord
  !********************************************************************************** *

 
 end module Trajectory
 !**********************************************************************************

 !**********************************************************************************
 Program Trajectory_Main
 use MD_SimboxArray_AppShell_16_GPU
 use Trajectory
 implicit none
 integer::numprocs=1, procid=0


       call APPSHELL_AddRecord(AFTRECORD=MyAftRecord, PRERECORD=MyPreRecord, RECORDPROC=MyRecord, AFTONETEST=MyAftOneTest)
       call APPSHELL_Main(numprocs,procid)
       !call Multi_ANALYSIS(0,ExitNofile=0)

       stop
 End program Trajectory_Main

