
 module MC_TypeDef_RandWalkerStatNumb
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The module is used to define the data type for statistical number of walkers
 !
 !**** HISTORY:
 !                  version 1st 2020-07 (Hou Qing, Sichuan university)
 !
 !**********************************************************************************   
   use MC_Randwalk_Const
   implicit none
     type RandWalkerStatNumb
          character(len=256)        ::PathOut     = ""
          integer                   ::hFile       = 0
          integer                   ::NBox        = 0                      !number of parallele boxes
          integer                   ::NWalkStyID  = 0                      !number of walk styles

          real(KINDDF)              ::RecTime                              ! time at which the recording is done   
          real(KINDDF)              ::Area                                 ! surface area when the recording is done, used to normalize the statistical number                  
          real(KINDDF)              ::Inserted(CP_MX_RADSTYLE)             ! total number of inserted particles
          real(KINDDF)              ::Discarded(CP_MX_RADSTYLE)            ! number of paticles discarded because of change of box
          real(KINDDF)              ::Reflected(CP_MX_RADSTYLE)            ! number of reflected paticles 
          real(KINDDF)              ::Transmit(CP_MX_RADSTYLE)             ! number of translateing paticles 
          real(KINDDF)              ::Trapped(CP_MX_RADSTYLE)              ! number of trapped paticles 
          real(KINDDF)              ::Walking(CP_MX_RADSTYLE)              ! number of paticles still walking
          real(KINDDF)              ::Waiting(CP_MX_RADSTYLE)              ! number of paticles that are inserted but triggered to walk

          real(KINDDF),  allocatable::DepthW(:)                            ! depth bins for depth distribution of walking particles
          real(KINDDF),  allocatable::DepthT(:)                            ! depth bins for depth distribution of trapped particles
          real(KINDDF),  allocatable::DepthHisW(:,:)                       ! depth histogram of depth distribution of walking particles
          real(KINDDF),  allocatable::DepthHisT(:,:)                       ! depth histogram of depth distribution of trapped particles
     end type RandWalkerStatNumb

     real(KINDDF),  private, allocatable::m_WKSZW(:)                       ! working space storing the depth of particles for depth histogram calculations of walking particles
     real(KINDDF),  private, allocatable::m_WKSZT(:)                       ! working space storing the depth of particles for depth histogram calculations of trapped particles
     integer,       private, allocatable::m_WKSIW(:)                       ! working space for sorting particle ID in walking
     integer,       private, allocatable::m_WKSIT(:)                       ! working space for sorting particle ID in traping
     integer,       private, allocatable::m_WKSSort(:)                     ! working space for sorting on calculating the histogram 

     !------------------------------------------------------------------------
     private::  Archive_RandWalkerStatNumb0, &
                Restore_RandWalkerStatNumb0
     public::   Archive_RandWalkerStatNumb
     interface  Archive_RandWalkerStatNumb
       module procedure  Archive_RandWalkerStatNumb0
     end interface  Archive_RandWalkerStatNumb               

     public::   Restore_RandWalkerStatNumb
     interface  Restore_RandWalkerStatNumb
       module procedure  Restore_RandWalkerStatNumb0
     end interface  Restore_RandWalkerStatNumb               

     !------------------------------------------------------------------------
     public::  Extract_RandWalkerStatNumb

     !------------------------------------------------------------------------
     public::  Initialize_RandWalkerStatNumb

     !------------------------------------------------------------------------
     private:: Putout_RandWalkerStatNumb0,  &
               Putout_RandWalkerStatNumb1
     public::  Putout_RandWalkerStatNumb
     interface Putout_RandWalkerStatNumb
        module procedure Putout_RandWalkerStatNumb0
          module procedure Putout_RandWalkerStatNumb1
     end interface Putout_RandWalkerStatNumb       
     !------------------------------------------------------------------------
     public::  Release_RandWalkerStatNumb 
     !------------------------------------------------------------------------
     public::  SetRecPath_RandWalkerStatNumb
 contains
 
  !****************************************************************************
  subroutine Release_RandWalkerStatNumb(StatNumber)
  !***  PURPOSE:   to release data  
  !
  !     INPUT:     StatNumber
  !
  !     OUTPUT     StatNumber, the object with memery deallocated
   implicit none
     !--- dummy varioables
     class(RandWalkerStatNumb), intent(inout) ::StatNumber
     !--- local variables
      
     !----
                
           if(StatNumber%hFile > 0) close(StatNumber%hFile) 
           StatNumber%hFile = 0;

           if(allocated(StatNumber%DepthW))    deallocate(StatNumber%DepthW)
           if(allocated(StatNumber%DepthT))    deallocate(StatNumber%DepthT)
           if(allocated(StatNumber%DepthHisW)) deallocate(StatNumber%DepthHisW)
           if(allocated(StatNumber%DepthHisT)) deallocate(StatNumber%DepthHisT)
           if(allocated(m_WKSZW))              deallocate(m_WKSZW)
           if(allocated(m_WKSZT))              deallocate(m_WKSZT)
           if(allocated(m_WKSIW))              deallocate(m_WKSIW)
           if(allocated(m_WKSIT))              deallocate(m_WKSIT)
           if(allocated(m_WKSSort))            deallocate(m_WKSSort)
                 
           StatNumber%NBox          = 0        
           StatNumber%NWalkStyID    = 0        
 
           return
  end subroutine Release_RandWalkerStatNumb
  !****************************************************************************

  !****************************************************************************
  subroutine Initialize_RandWalkerStatNumb(StatNumber, CtrlParam, Walker)
  !***  PURPOSE:   to initialize state of the StatNumber
  !
  !     INPUT:     NBOX,      number of box in the workspace
  !                MXWID,     max number of particle types    
  !
  !     OUTPUT     StatNumber, with memoery allocated and initialized
   use MC_TypeDef_RandWalkCtrl
   use MC_TypeDef_RandWalker_Base
   implicit none
     !--- dummy varioables
     type(RandWalkCtrl)        :: CtrlParam
     type(WalkerBase)          :: Walker(:)
     class(RandWalkerStatNumb) ::StatNumber
     !--- local variables
     integer::NBOX, MXNP
     !---
          call Release_RandWalkerStatNumb(StatNumber)
          call GetMXWalkerNum_RandWalkCtrl(CtrlParam, MXNP, NBOX) 

          StatNumber%NBox       = NBOX
          StatNumber%NWalkStyID = size(Walker)
          allocate(&
                   StatNumber%DepthW(0:CtrlParam%RecDepthBin),                 &
                   StatNumber%DepthT(0:CtrlParam%RecDepthBin),                 &
                   StatNumber%DepthHisW(CtrlParam%RecDepthBin,size(Walker)),   &
                   StatNumber%DepthHisT(CtrlParam%RecDepthBin,size(Walker))    &
                   )
          
          StatNumber%Inserted  = 0
          StatNumber%Reflected = 0
          StatNumber%Transmit  = 0
          StatNumber%Trapped   = 0
          StatNumber%Walking   = 0
          StatNumber%Waiting   = 0
          StatNumber%Discarded = 0

          StatNumber%DepthW    = 0
          StatNumber%DepthT    = 0
          StatNumber%DepthHisW = 0
          StatNumber%DepthHisT = 0
          StatNumber%Area      = (CtrlParam%Boxup(1) - CtrlParam%Boxlow(1))*(CtrlParam%Boxup(2) - CtrlParam%Boxlow(2)) 
          return
  end subroutine Initialize_RandWalkerStatNumb   
  !****************************************************************************

  !****************************************************************************
  subroutine Extract_RandWalkerStatNumb(Time, CtrlParam, WalkerState, StatNumber)
  !***  PORPOSE: to extract the statistical numbers from walker stat
  !        
  !     INPUT:   Time,           the current time 
  !              WalkerState,    the current state of walker
  !  
  !     OUTPUT:  StatNumber,     
  
  !
   use MC_TypeDef_RandWalkCtrl
   use MC_TypeDef_RandWalkerState
   use MiniUtilities, only:dqsort
   implicit none
      !--- DUMMY variables
      real(KINDDF),              intent(in)    :: Time 
      type(RandWalkCtrl),        intent(in)    :: CtrlParam
      class(RandWalkerStat),     intent(in)    :: WalkerState
      class(RandWalkerStatNumb), intent(inout) :: StatNumber
      !--- local variables
      integer::IB, I,J, WID, NBIN, CURNP, IBIN, TWNP, TTNP, TNP, NP, FROM, TO
      real(KINDDF)::Area 
      !-----  
          
         if(.not.allocated(m_WKSZW))   allocate(m_WKSZW(size(WalkerState%Stat))  )
         if(.not.allocated(m_WKSZT))   allocate(m_WKSZT(size(WalkerState%Stat))   )
         if(.not.allocated(m_WKSIW))   allocate(m_WKSIW(size(WalkerState%Stat))   )
         if(.not.allocated(m_WKSIT))   allocate(m_WKSIT(size(WalkerState%Stat))   )
         if(.not.allocated(m_WKSSort)) allocate(m_WKSSort(size(WalkerState%Stat)) )
          StatNumber%RecTime   = Time
          StatNumber%Reflected = 0
          StatNumber%Transmit  = 0
          StatNumber%Trapped   = 0
          StatNumber%Walking   = 0
          StatNumber%Waiting   = 0
          !---- the reflected anf transmitting particles
          TWNP = 0
          TTNP = 0
          do IB=1,  WalkerState%NBox
              I    = WalkerState%FirstPIDP(IB)
              CURNP= WalkerState%CurNP(IB)
              do J = 1, CURNP 
                 WID = WalkerState%WalkStyID(I)
                 select case(WalkerState%Stat(I))
                 case(CP_WALKSTAT_REFLECT)
                       StatNumber%Reflected(WID) = StatNumber%Reflected(WID) + 1
                 case(CP_WALKSTAT_TRANS)
                       StatNumber%Transmit(WID)  = StatNumber%Transmit(WID)  + 1
                 case (CP_WALKSTAT_TRAPPED)
                       StatNumber%Trapped(WID)   = StatNumber%Trapped(WID)   + 1
                       TTNP           = TTNP + 1
                       m_WKSIT(TTNP) = I
                       m_WKSZT(TTNP) = WalkerState%XP(I,3)
                 case (CP_WALKSTAT_ACTIVE)
                       StatNumber%Walking(WID)   = StatNumber%Walking(WID)   + 1
                       TWNP          = TWNP + 1
                       m_WKSIW(TWNP) = I
                       m_WKSZW(TWNP) = WalkerState%XP(I,3)
                 case (CP_WALKSTAT_WAITING)
                       StatNumber%Waiting(WID) = StatNumber%Waiting(WID) + 1
                 end select      
                I = I + 1
             end do
          end do  
          !---- toe normalize the numbers to area
          StatNumber%Reflected = StatNumber%Reflected/StatNumber%Area/dble(StatNumber%NBox)
          StatNumber%Transmit  = StatNumber%Transmit/StatNumber%Area/dble(StatNumber%NBox)
          StatNumber%Trapped   = StatNumber%Trapped/StatNumber%Area/dble(StatNumber%NBox)
          StatNumber%Walking   = StatNumber%Walking/StatNumber%Area/dble(StatNumber%NBox)
          StatNumber%Waiting   = StatNumber%Waiting/StatNumber%Area/dble(StatNumber%NBox)

          !---- extract the depth distribution particles
          StatNumber%DepthW    = 0
          StatNumber%DepthT    = 0
          StatNumber%DepthHisW = 0
          StatNumber%DepthHisT = 0
          m_WKSSort   = 0
          if(TWNP .gt. 0) then
             NBIN  = size(StatNumber%DepthW)-1
             NP    = max(2,TWNP/NBIN)

             call dqsort(TWNP, m_WKSZW(1:TWNP), m_WKSSort(1:TWNP))
             StatNumber%DepthW(0) = m_WKSZW(m_WKSSort(1))
             do IBIN = 1, NBIN 
                FROM = min((IBIN-1)*NP, TWNP)
                do J = FROM+1, min(FROM+NP, TWNP)
                   WID = WalkerState%WalkStyID( m_WKSIW(m_WKSSort(J)) )
                   StatNumber%DepthHisW(IBIN,WID) = StatNumber%DepthHisW(IBIN,WID) + 1
                   StatNumber%DepthW(IBIN)        = m_WKSZW(m_WKSSort(J))
                  end do  
             end do
             StatNumber%DepthHisW = StatNumber%DepthHisW/StatNumber%Area/dble(StatNumber%NBox) 

          end if  
          return 
    end subroutine Extract_RandWalkerStatNumb
  !**************************************************************************** 

  !****************************************************************************
  subroutine SetRecPath_RandWalkerStatNumb(Fname, StatNumber)
  !***  PURPOSE:   to output the output file 
  !
  !     INPUT:     Fname,      the file name
  !
  !     OUTPUT     
   implicit none
     !--- dummy varioables
     character*(*),             intent(in)    ::Fname             
     class(RandWalkerStatNumb), intent(inout) ::StatNumber
     !--- local variables
     character*64::FMT, SN 
     !----
           if(StatNumber%hFile > 0) close(StatNumber%hFile)
           StatNumber%hFile  = 0
           StatNumber%PathOut = Fname(1:len_trim(fname))

           return
  end subroutine SetRecPath_RandWalkerStatNumb
  !****************************************************************************    

  !****************************************************************************
  subroutine Putout_RandWalkerStatNumb0(hFile, StatNumber)
  !***  PURPOSE:   to output the data  
  !
  !     INPUT:     hFile,      I/O unit
  !                StatNumber
  !
  !     OUTPUT     
   implicit none
     !--- dummy varioables
     integer,                   intent(in) ::hFile             
     class(RandWalkerStatNumb), intent(in) ::StatNumber
     !--- local variables
      integer::I
      character*64::FMT, SN
      real(KINDDF)::FLUEN(CP_MX_RADSTYLE), REFL(CP_MX_RADSTYLE), TRANS(CP_MX_RADSTYLE),TRAP(CP_MX_RADSTYLE), WALK(CP_MX_RADSTYLE)
     !----
           do I=1, StatNumber%NWalkStyID
              !    Note: not all the inserted particles are in wating stat, the real fluence should
              !          be inserted partciles - waiting particles  
              FLUEN(I) = StatNumber%Inserted(I) - StatNumber%Waiting(I)
              if(FLUEN(I) .gt. 0.D0) then
                 REFL(I)    = StatNumber%Reflected(I)/FLUEN(I)
                 TRANS(I)   = StatNumber%Transmit(I)/FLUEN(I)
                 TRAP(I)    = StatNumber%Trapped(I)/FLUEN(I)
                 WALK(I)    = StatNumber%Walking(I)/FLUEN(I)
              else
                 FLUEN(I)   = 0.D0
                 REFL(I)    = 0.D0
                 TRANS(I)   = 0.D0
                 TRAP(I)    = 0.D0
                 WALK(I)    = 0.D0
              end if
           end do
           

           write(SN,*) StatNumber%NWalkStyID
           SN =adjustl(SN)
           FMT ="(1PE16.4, 1x, 1PE16.4,1x, " //SN(1:len_trim(SN))//"(5(1PE14.4,2X)) )"
           write(hFile,fmt=FMT) StatNumber%RecTime,  StatNumber%Area,             &
                                (  (FLUEN(I), REFL(I), TRANS(I), TRAP(I), WALK(I) ), I=1, StatNumber%NWalkStyID)
           flush(hFile)
           return
  end subroutine Putout_RandWalkerStatNumb0
  !****************************************************************************  

  !****************************************************************************
  subroutine Putout_RandWalkerStatNumb1(StatNumber, Restart)
  !***  PURPOSE:   to output the data  
  !
  !     INPUT:     hFile,      I/O unit
  !                StatNumber
  !
  !     OUTPUT     
   use MiniUtilities 
   implicit none
     !--- dummy varioables
     class(RandWalkerStatNumb), intent(in) ::StatNumber
     integer,                   intent(in) ::Restart
     !--- local variables
     integer::I
     character*64::FMT, SN
     character*256::Fname
     !----
            if(StatNumber%hFile <= 0) then
               call AvailableIOUnit(StatNumber%hFile)
               #ifdef WIN_FILE
               Fname = StatNumber%PathOut(1:len_trim(StatNumber%PathOut))//achar(92)
               #else
               Fname = StatNumber%PathOut(1:len_trim(StatNumber%PathOut))//achar(47)
               #endif
               Fname = Fname(1:len_trim(Fname))//"AccumNumberss.dat"

               if(Restart) then
                  open(UNIT=StatNumber%hFile, FILE=Fname, POSITION="append")  
               else
                  open(UNIT=StatNumber%hFile, FILE=Fname)
                  write(SN,*) StatNumber%NWalkStyID
                  SN =adjustl(SN)
                  FMT ="(A16, 1X, A16,1X, " //SN(1:len_trim(SN))//"(5(A14,2X)) )"
                  write(StatNumber%hFile, fmt=FMT) "RECTIME(ps)", "SURF.AREA(nm2)", (("FLUENCE", "REFLECT", "TRANSMIT", "TRAPPED", "WALKING"), &
                                                                                   I=1,StatNumber%NWalkStyID)
               end if                                                                    
            end if   

            call Putout_RandWalkerStatNumb0(StatNumber%hFile, StatNumber)
           return
  end subroutine Putout_RandWalkerStatNumb1
  !****************************************************************************  

  !****************************************************************************
  subroutine PutoutDepthHis_RandWalkerStatNumb0(hFile, Depth, His)
  !***  PURPOSE:   to output the depth histogram
  !
  !     INPUT:     StatNumber
  !
  !     OUTPUT   
   implicit none
     !--- dummy varioables
      integer,      intent(in)::hFile
      real(KINDDF), intent(in)::Depth(0:)      
      real(KINDDF), intent(in)::His(:,:)      
!--- local variables
     integer::NW, I,J    
     character*64::FMT, SN
     real(KINDDF)::DW
     !----
         NW = size(His,dim=2)
         write(SN,*) NW
         FMT ="(A16,1X," //SN(1:len_trim(SN))//"(2(A16,1X)) )"
         write(hFile, fmt=FMT)   "Depth",   (("Count",    "His"), J=1,NW)

         FMT ="(1PE16.5,1X," //SN(1:len_trim(SN))//"(2(1PE16.5,1X)) )"
         write(hFile, fmt=FMT) Depth(0), ((0.D0, 0.D0),J=1,NW)
         do I=1, size(Depth)-1
            if(Depth(I) > Depth(I-1)) then
               DW =  Depth(I) - Depth(I-1)
            else
               exit   
            end if   
            write(hFile, fmt="(12(1PE16.5,1X))") Depth(I-1), ((His(I,J), His(I,J)/DW), J=1, NW)
            write(hFile, fmt="(12(1PE16.5,1X))") Depth(I),   ((His(I,J), His(I,J)/DW), J=1, NW)
         end do 
         I = I -1
         write(hFile, fmt="(12(1PE16.5,1X))") Depth(I),   ((0.D0, 0.D0), J=1, NW)
         return
  end subroutine PutoutDepthHis_RandWalkerStatNumb0
  !****************************************************************************  


  !****************************************************************************
  subroutine PutoutDepthHis_RandWalkerStatNumb(StatNumber, Stamp)
  !***  PURPOSE:   to output the depth histogram
  !
  !     INPUT:     StatNumber
  !
  !     OUTPUT   
   use MSM_TYPEDEF_RecordStamp  
   use MiniUtilities
   implicit none
     !--- dummy varioables
     class(RandWalkerStatNumb), intent(in) ::StatNumber
     type(MSMRecordStamp)                  :: Stamp
     !--- local variables
     integer::hFile
     character*64::FMT, SN
     character*256::Fname
     real(KINDDF)::DWW, DWT
     !----

         !--- to output hisgram for walking particles
          if(Stamp%IRec(1) .ge. 0) then
             call STRCATI(fname, StatNumber%PathOut, "DepthHisW.", Stamp%IRec(1), 4)
          end if   

         !--- to open the file
         call AvailableIOUnit(hFile)
         open(hFile, file = fname, status='unknown')
         call PutoutDepthHis_RandWalkerStatNumb0(hFile, StatNumber%DepthW, StatNumber%DepthHisW)
         close(hFile)

         !--- to output hisgram for trapped particles
          if(Stamp%IRec(1) .ge. 0) then
             call STRCATI(fname, StatNumber%PathOut, "DepthHisT.", Stamp%IRec(1), 4)
          end if   

         !--- to open the file
         call AvailableIOUnit(hFile)
         open(hFile, file = fname, status='unknown')
         call PutoutDepthHis_RandWalkerStatNumb0(hFile, StatNumber%DepthT, StatNumber%DepthHisT)
         close(hFile)

         return
  end subroutine PutoutDepthHis_RandWalkerStatNumb
  !****************************************************************************  

  !****************************************************************************
  subroutine Archive_RandWalkerStatNumb0(hFile, StatNumber)
  !***  PURPOSE:   to archive the current the data  for restart
  !
  !     INPUT:     hFile,      I/O unit
  !                StatNumber
  !
  !     OUTPUT     
   implicit none
     !--- dummy varioables
     integer,                   intent(in) ::hFile             
     class(RandWalkerStatNumb), intent(in) ::StatNumber
     !--- local variables
     !----

             write(hFile) StatNumber%NBox,      &
                          StatNumber%NWalkStyID,&
                          StatNumber%RecTime,   &
                          StatNumber%Area,      &
                          StatNumber%Inserted,  &
                          StatNumber%Discarded, &
                          StatNumber%Reflected, &
                          StatNumber%Transmit,  &
                          StatNumber%Trapped,   &
                          StatNumber%Walking,   &
                          StatNumber%Waiting
           return
  end subroutine Archive_RandWalkerStatNumb0
  !****************************************************************************  

  !****************************************************************************
  subroutine Restore_RandWalkerStatNumb0(hFile, StatNumber)
  !***  PURPOSE:   to archive the current the data  for restart
  !
  !     INPUT:     hFile,      I/O unit
  !                StatNumber
  !
  !     OUTPUT     
   implicit none
     !--- dummy varioables
     integer,                   intent(in) ::hFile             
     class(RandWalkerStatNumb), intent(in) ::StatNumber
     !--- local variables
     !----

             read(hFile)  StatNumber%NBox,      &
                          StatNumber%NWalkStyID,&
                          StatNumber%RecTime,   &
                          StatNumber%Area,      &
                          StatNumber%Inserted,  &
                          StatNumber%Discarded, &
                          StatNumber%Reflected, &
                          StatNumber%Transmit,  &
                          StatNumber%Trapped,   &
                          StatNumber%Walking,   &
                          StatNumber%Waiting
           return
  end subroutine Restore_RandWalkerStatNumb0
  !****************************************************************************  

 end module MC_TypeDef_RandWalkerStatNumb 
 !********************************************************************************

