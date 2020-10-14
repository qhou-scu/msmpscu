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
 
       character(len=13),parameter, private::mp_FTAGI="&AUXF_TRAJECTORY"
      !--- NOTE the variable DefTyp subjected to changes
       integer,                           private::m_Nearest =8

       character*256,                     private::m_CfgFname  = ""
       type(MDRecordStamp),               private::m_CurStamp
       type(MDDEVNeighbor_List),          private::m_dNNList
       integer,                           private::m_TOTREC = 0
       integer,                           private::m_NDEF   = 0
       integer,                           private::m_DEFFLG(mp_MXGROUP) = 0
       integer,                           private::m_DEFTYP(mp_MXGROUP) = 0
       integer, dimension(:),allocatable, private::m_DefAtomID
       type(StatementList),pointer,       private::m_ExInfor=>null()
       integer,                           private::m_RECN = 0

       type(SimMDBox),      dimension(:), allocatable, private::m_MergedSimBox
       type(Neighbor_List), dimension(:), allocatable, private ::m_NNList

       private::Clear
       private::PrintParameter
       private::RecordTrack
  contains
  
   !**********************************************************************************
   subroutine Initialize_Trajectory(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation 
   use MD_TYPEDEF_PrintList, only:Add_ArchiveProcess    
   implicit none
   !----   DUMMY Variables
   type(SimMDBox)  ::SimBox
   type(SimMDCtrl) ::CtrlParam

   !----   Local variables
   integer::I, IFILE, STATU, LINE, N
   type(InputStatements)::INPUTS
   character*256::STR
   character*32::STRNUMB(mp_MXGROUP)

           IFILE = 0
           !--- used in analysis tool
            if(len_trim(gm_cfgFileName(1)) .gt. 0) then
               m_CfgFname = gm_cfgFileName(1)
               m_TOTREC  = (CtrlParam%ENDCFG - CtrlParam%STARTCFG)/CtrlParam%CFGSTEP + 1
            !--- used in MD simualtions   
            else  
               m_CfgFname = CtrlParam%f_geometry
               call NumberRec_SimMDCtrl(CtrlParam, m_TOTREC)

              !$$--- to findout the I/O unit
               do I=1, size(CtrlParam%f_tag)
                  if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                     IFILE = I
                     exit
                  end if
               enddo
            end if   
            !-- by defalut value
            m_Nearest   = 8
            m_DEFFLG    = 0
            if(IFILE .le. 0) then
               m_DEFFLG(2) = 1
               m_DEFFLG(3) = 1
            else
               !$$--- load input statements
               write(*,*) "!**** Load control parameters for trajectory simulation from:"
               write(*,*) "!**** ",CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
               write(*,*)
     
               call LoadExInput_SimMDCtrl(mp_FTAGI, CtrlParam%f_others(IFILE), CtrlParam)
               call GetExInput_SimMDCtrl(CtrlParam, mp_FTAGI, INPUTS, STATU)
               !--- to extract the control parameters
               call Get_InputStatements("&NEARESTNB", INPUTS, STR, LINE)
               if(LINE .gt. 0) then
                  call Extract_Numb(STR,1,N,STRNUMB)
                  m_Nearest = IStr(STRNUMB(1))
               end if    
               
               call Get_InputStatements("&DEFATOM", INPUTS, STR, LINE)
               if(LINE .gt. 0) then
                  call Extract_Numb(STR,mp_MXGROUP,N,STRNUMB)
                  do I=1, N
                     STATU = IStr(STRNUMB(I))
                     if(STATU .gt. 0) m_DEFFLG(STATU) = 1
                  end do   
               end if    
            end if 
            
            STATU    = 0 
            m_DEFTYP = 0
            do I=1, size(m_DEFFLG) 
               if(m_DEFFLG(I).gt. 0)then
                 STATU = STATU + 1
                 m_DEFTYP(STATU) = I
               end if   
            end do  
            
            m_CurStamp%ITime = - 1
            m_CurStamp%ITest = 0 
            m_RECN           = 0
            write(6,*) "!************ TRAJECTORY OUTPUT CONTROL PARAMETERS **********"
            call PrintParameter(6)
            if(gm_hFILELOG .gt. 0) call PrintParameter(gm_hFILELOG)            
            
            !---
            call Add_ArchiveProcess("&TRAJECT", Archive, Restore)
         return
   end subroutine Initialize_Trajectory
  !********************************************************************************** 

  !****************************************************************************
  subroutine PrintParameter(hFile)
  !***  PURPOSE:   to print projectile parameters into I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile

     !--- local variables
             write(hFile, fmt="(A,100(4I,1X))") ' !    Types of defect atoms................: ', m_DEFTYP(1:count(m_DEFTYP.gt.0))
             write(hFile, fmt="(A,100(4I,1X))") ' !    Nearest atoms to output..............: ', m_Nearest
  end subroutine PrintParameter
  !**********************************************************************************

  !****************************************************************************
  subroutine Archive(hFile)
   !***  PURPOSE:   to save the current status of this module
   !
   !     INPUT:     hFile,  I/O unit number
   !
   !     OUTPUT
   implicit none
      !--- dummy varioables
      integer, intent(in)::hFile
 
      !--- local variables
      integer::I, N
              write(hFile) m_RECN
              write(hFile) m_NDEF
              write(hFile) m_DEFFLG
              write(hFile) m_DEFTYP
              call Archive_RecordStamp(hFile, m_CurStamp)
              call Archive_StatementList(hFile, m_ExInfor)
              if(m_NDEF .gt. 0) then
                 write(hFile) m_DefAtomID(1:m_NDEF)
              end if   
              N = 0
              if(allocated(m_MergedSimBox)) N = size(m_MergedSimBox)
              write(hFile) N
              do I=1, N
                 call Archive_Config_SimMDBox(hFile, m_MergedSimBox(I))
                 call Archive_DataPad(hFile,m_MergedSimBox(I)%ptrDatPad)
              end do   
              return
  end subroutine Archive
  !**********************************************************************************

  !****************************************************************************
  subroutine Restore(hFile)
  !***  PURPOSE:   to print projectile parameters into I/O unit
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
   implicit none
      !--- dummy varioables
      integer, intent(in)::hFile
    
      !--- local variables
      integer::I, N

         read(hFile) m_RECN
         read(hFile) m_NDEF
         read(hFile) m_DEFFLG
         read(hFile) m_DEFTYP
         call Restore_RecordStamp(hFile, m_CurStamp)
         call Release_StatementList(m_ExInfor)
         call Restore_StatementList(hFile, m_ExInfor)
         if(m_NDEF .gt. 0) then
           if(allocated(m_DefAtomID)) deallocate(m_DefAtomID)
           allocate(m_DefAtomID(m_NDEF))
           read(hFile) m_DefAtomID(1:m_NDEF)
         end if
         read(hFile) N
         if(N .gt. 0) then
            if(allocated(m_MergedSimBox)) then
               call Release_SimBoxArray(m_MergedSimBox)  
               deallocate(m_MergedSimBox) 
            end if       
            allocate(m_MergedSimBox(N), m_NNList(N))
         end if    
         do I=1, N
            call Restore_Config_SimMDBox(hFile, m_MergedSimBox(I))
            call Restore_DataPad(hFile,m_MergedSimBox(I)%ptrDatPad)
         end do   
         return 
  end subroutine Restore
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

              call AddAtoms_SimMDBox(SimBoxTrack, (/NDef,NDef*NNb/),  (/2,1/))
              IP0 = max(0, NPRT0-(NDef+NDef*NNb))
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
  subroutine Record_Trajectory(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calcualte and putout the occupation state of lattice and interstials
  !                  identified using method#1
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
   
   implicit none
       !--- dummy variables
       type(MDRecordStamp) ,           intent(in)::Stamp
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables
       integer::IB, I, K, IP, ITY, NB
       character*256 ::GFILE, STR

       !***************
          if(Stamp%ITIME .le. 0) return

          if(m_CurStamp%ITest    .ne. Stamp%ITest .or. &
             any(m_CurStamp%IBox .ne. Stamp%Ibox) .or. &
             m_RECN .lt. 1 ) then

             call Clear()  
             NB = size(SimBox)	
             allocate(m_MergedSimBox(NB), m_NNList(NB))
             
             !--- determine ID of defects
             m_NDEF = 0 !    = SimBox(1)%NA(m_DefTyp) !count(SimBox(1)%ITYP .eq. m_DefTyp)
             do I=1, size(m_DEFFLG) 
                if(m_DEFFLG(I).gt. 0)then
                  m_NDEF = m_NDEF + SimBox(1)%NA(I)
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
             m_RECN     = 1
             call Reorder_NeighBoreList_Nearest_Dev(m_Nearest, m_dNNList, m_NNList)
             do IB = 1, NB
                call Initialize_SimMDBox((m_Nearest+1)*m_NDEF*m_TOTREC, m_MergedSimBox(IB))
                call NewDatPad_SimMDBox( m_MergedSimBox(IB), "IDN",   "I")
                call NewDatPad_SimMDBox( m_MergedSimBox(IB), "DVEC",   "D", 3)
                !call NewDatPad_SimMDBox( m_MergedSimBox(IB), "TIME",   "D")
                call CopyInformation_SimMDBox(SimBox(IB), m_MergedSimBox(IB))
                call RecordTrack(SimBox(IB), CtrlParam, m_MergedSimBox(IB), m_NDEF, m_DefAtomID, m_Nearest, m_NNList(IB), Stamp%TIME)
             end do
             
             write(STR,*) "&TOTREC  ", m_TOTREC
             call Add_StatementList(m_ExInfor, STR)
             write(STR,*) "&NDEF    ",  m_NDEF, m_DEFTYP(1:count(m_DEFTYP.gt.0))
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
          m_RECN = m_RECN + 1 
          call Reorder_NeighBoreList_Nearest_Dev(m_Nearest, m_dNNList, m_NNList)
          do IB = 1, size(m_MergedSimBox)
             call RecordTrack(SimBox(IB), CtrlParam, m_MergedSimBox(IB), m_NDEF, m_DefAtomID, m_Nearest, m_NNList(IB), Stamp%TIME)
          end do   
          
          !if(m_RECN .eq.  m_TOTREC) 	then
          !   do IB=1, NB
          !      call STRCATI(GFILE, m_CfgFname, "_Track_P", 0, 4)
          !      K = (m_CurStamp%ITest - 1)*NB + IB
          !      call STRCATI(GFILE, GFILE, "_", K, 4)
          !      m_CurStamp%ICfg = -1
          !      call Putout_Instance_Config_SimMDBox(GFILE, m_MergedSimBox(IB), m_CurStamp, ExStats=m_ExInfor)
          !   end do    
          !end if	
          
 end subroutine Record_Trajectory
 !***************************************************************************
 
 !**********************************************************************************
 subroutine AftOneTest_Trajectory(SimBox0, Stamp, SimBox, CtrlParam)
 !***  PURPOSE:  to do some thing after on test finished
 implicit none
 !----   DUMMY Variables
   type(SimMDBox)              :: SimBox0
   type(MDRecordStamp)           :: Stamp
   type(SimMDBox), dimension(:):: SimBox
   type(SimMDCtrl)             :: CtrlParam

  !----   Local variables
   character*256 ::GFILE
   integer::IB, K, NB

              NB = size(m_MergedSimBox) 
              do IB = 1, NB
                 call STRCATI(GFILE, m_CfgFname, "_Track_P", 0, 4)
                 K = (m_CurStamp%ITest - 1)*NB + IB
                 call STRCATI(GFILE, GFILE, "_", K, 4)
                 m_CurStamp%ICFG = -1
                 call Putout_Instance_Config_SimMDBox(GFILE, m_MergedSimBox(IB), m_CurStamp, ExStats=m_ExInfor)
              end do
              call Clear()
              return
   end subroutine AftOneTest_Trajectory
 !********************************************************************************** *

 !**********************************************************************************
   subroutine Clear()
   !***  PURPOSE:  to clear the allocated memories
   implicit none
  !----   DUMMY Variables

  !----   Local variables
   integer::IB
              if(allocated(m_MergedSimBox)) then
                call Release_SimBoxArray(m_MergedSimBox)
                deallocate(m_MergedSimBox)
                IB = size(m_MergedSimBox)
                IB = allocated(m_MergedSimBox)
              end if  
              call Release_StatementList(m_ExInfor)

              call Clear_NeighboreList_Dev(m_dNNList)
              if(allocated(m_NNList)) then
                 do IB=1, size(m_NNList)
                    call Clear_NeighboreList(m_NNList(IB))
                 end do
                 deallocate(m_NNList)
               end if
               
               m_CurStamp%ITime = - 1
               m_RECN           = 0
                  
              return
   end subroutine Clear
 !********************************************************************************** *   

 end module Trajectory
 !**********************************************************************************

 !**********************************************************************************
 Program Trajectory_Main
 use MD_SimboxArray_AppShell_16_GPU
 use Trajectory
 implicit none
 integer::numprocs=1, procid=0


       call APPSHELL_AddRecord(PRERECORD=Initialize_Trajectory, RECORDPROC=Record_Trajectory, AFTONETEST=AftOneTest_Trajectory)
       call APPSHELL_Main(numprocs,procid)
       !call Multi_ANALYSIS(0,ExitNofile=0)

       stop
 End program Trajectory_Main

