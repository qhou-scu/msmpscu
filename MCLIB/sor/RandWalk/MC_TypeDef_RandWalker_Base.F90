
 module MC_TypeDef_RandWalker_Base
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The module is used to simulate random walk process. The waiting time ditribution
 !                  follows the exponential distribution  EXP(-t/tau)/tau 
 !
 !
 !**** HISTORY:
 !                  version 1st 2020-05 (Hou Qing, Sichuan university)
 !
 !**********************************************************************************   
   use MSM_MultiGPU_Basic
   use MiniUtilities
   use MC_Randwalk_Const
   implicit none

   !---- basic descriptor of Jumps
     type Jump_Base
          !--- parameters concerning Jump sampling 
          integer                ::WTD_TYPE                          =  CP_WALKTYPE_WTD_EXP
          real(KINDDF)           ::WTD_TAU                           = 1.D0
          real(KINDDF)           ::WTD_ALPHA                         = 3.D0

          integer                ::DISD_TYPE                         = CP_WALKTYPE_DIS_FIXSTEP 
          real(KINDDF)           ::DISD_AV                           = 1.D0                       !average step size for one jump
          real(KINDDF)           ::DISD_BETA                         = 3.D0                       !beta value if the step is a Levy flight

          integer                ::JP_NVEC                           = 0                          ! number of possible jump vector
                                                                                                  ! for example for SIA in W
                                                                                                  ! JP_DIM = 1, three linearly independent vectors can 
                                                                                                  ! be defined:(1,1,1), (1,1,-1), (1,-1,-1).
                                                                                                  ! the jumps can switch between these vectors
          real(KINDDF)           ::JP_VEC(3,CP_WALKTYPE_MAXJVECT)    = 0                          ! Jumping vectors
     end type Jump_Base

     !---- basic descriptor of walker
     type WalkerBase
          !--- parameters concerning walk sampling 
                     
          character*8                       ::Name                              =""  
          integer                           ::WalkType                          = CP_WALKTYPE_DIM_1D 
          integer                           ::NJump                             = 1               !number of possible jumps  
          real(KINDDF)                      ::JumpProb(CP_WALKTYPE_MAXJPPATH)   = 0               !Jump probabilty for a jump
          real(KINDDF)                      ::JumpProbI(CP_WALKTYPE_MAXJPPATH)  = 0               !Integral and normalized Jump probabilty for a jump
          type(Jump_Base)                   ::JumpType(CP_WALKTYPE_MAXJPPATH)                     !Jumping vector
     end type WalkerBase
     !------------------------------------------------------     

     !------------------------------------------------------
     private:: &
               Load_JumpBase0,            &
               Load_WalkerBase0,          & 
               Load_WalkerBase1,          &
               Load_WalkerBase2,          &
               Load_WalkerBase3,          &
               Load_WalkerBase4
 
     public::  Load_RandWalker
     interface Load_RandWalker
          module procedure Load_WalkerBase0
     end interface   

     public::  LoadSetup_RandWalker
     interface LoadSetup_RandWalker
          module procedure Load_WalkerBase1
          module procedure Load_WalkerBase2
          module procedure Load_WalkerBase3
          module procedure Load_WalkerBase4
     end interface          

 contains
 !****************************************************************************

 !****************************************************************************
  subroutine Load_WalkerBase0(Str, Walker)
  !***  PURPOSE:   to load the type of  walker
  !
  !     INPUT:     STR,    a string line
  !
  !     OUTPUT     Walker

   implicit none
     !--- dummy varioables
     character*(*),     intent(in)    ::Str
     class(WalkerBase), intent(inout) ::Walker
     !--- local variables
      character*256::STRTMP(10)=""
      integer::N
     !----
           call extract_optstr(Str, "-","WTD", 1, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "&","WTD", 1, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "-","WT", 1, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "&","WT", 1, N, STRTMP)
           if(N .ge. 1) then
               call UpCase(STRTMP(1))
               select case(STRTMP(1))
               case ("EXP")
                    Walker%JumpType(1)%WTD_TYPE = CP_WALKTYPE_WTD_EXP
                    call extract_optstr(Str, "-","TAU", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","TAU", 1, N, STRTMP)
                    if(N .ge. 1) then
                         Walker%JumpType(1)%WTD_TAU = DRSTR(STRTMP(1))   
                    end if
               case ("POW")
                    Walker%JumpType(1)%WTD_TYPE = CP_WALKTYPE_WTD_POW
                    call extract_optstr(Str, "-","TAU", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","TAU", 1, N, STRTMP)
                    if(N .ge. 1) then
                         Walker%JumpType(1)%WTD_TAU = DRSTR(STRTMP(1))   
                    end if
     
                    call extract_optstr(Str, "-","ALPHA", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","ALPHA", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "-","A", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","A", 1, N, STRTMP)
                    if(N .ge. 1) then
                         Walker%JumpType(1)%WTD_Alpha = DRSTR(STRTMP(1))   
                    end if
               end select     
           end if    

           !---
           call extract_optstr(Str, "-","STEPTYPE", 1, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "&","STEPTYPE", 1, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "-","ST", 1, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "&","ST", 1, N, STRTMP)
           if(N .ge. 1) then
               call UpCase(STRTMP(1))
               select case(STRTMP(1))
               case ("FIX")
                    Walker%JumpType(1)%DISD_TYPE = CP_WALKTYPE_DIS_FIXSTEP
                    call extract_optstr(Str, "-","STEP", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","STEP", 1, N, STRTMP)
                    if(N .ge. 1) then
                         Walker%JumpType(1)%DISD_AV = drstr(STRTMP(1))   
                    end if

               case ("POW")
                    Walker%JumpType(1)%DISD_TYPE = CP_WALKTYPE_DIS_POWDIS
                    call extract_optstr(Str, "-","STEP", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","STEP", 1, N, STRTMP)
                    if(N .ge. 1) then
                         Walker%JumpType(1)%DISD_AV = drstr(STRTMP(1))   
                    end if

                    call extract_optstr(Str, "-","BETA", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","BETA", 1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "-","B",    1, N, STRTMP)
                    if(N .eq. 0) call extract_optstr(Str, "&","B",    1, N, STRTMP)
                    if(N .ge. 1) then
                         Walker%JumpType(1)%DISD_BETA = drstr(STRTMP(1))   
                    end if
               end select     
           end if    

           return
  end subroutine Load_WalkerBase0
  !****************************************************************************

  !****************************************************************************
  subroutine Load_JumpBase0(hFile, Walker, IJ, Line0)
  !***  PURPOSE:   to load the IVectth jump base of a walker
  !
  !     INPUT:     hFile,  I/O unit number
  !                IJ,     the jump number
  !
  !     OUTPUT     Walker
  !                Line0   
  implicit none
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(WalkerBase),   intent(inout) ::Walker
     integer,             intent(inout) ::IJ
     integer,optional,    intent(inout) ::Line0

     !--- local variables
      character*256::STR
      character*32::STRNUMB(4), KEYWORD
      integer:: N, LINE, ERR=0, J, JN

         if(present(LINE0)) then
            LINE = LINE0
          else
            LINE = 0
          end if

          do while(.TRUE.)
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             STR = adjustl(STR)
             call GetKeyWord("&", STR, KEYWORD)
             call UpCase(KEYWORD)
             select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )

             case("&RATE")
                  call Extract_Numb(STR,1,N,STRNUMB)
                  Walker%JumpProb(IJ) = Drstr(STRNUMB(1))

             case("&WTD")
                   call Extract_substr(STR,1,N,STRNUMB)
                   call UpCase(STRNUMB(1))
                   select case(STRNUMB(1))
                   case ("EXP")
                        Walker%JumpType(IJ)%WTD_TYPE = CP_WALKTYPE_WTD_EXP
                        call Extract_Numb(STR,1,N,STRNUMB)
                        if(N .lt. 1) then
                           write(*,fmt="(A,BZI6)") 'MCPSCU Error: parameter TAU missed at line', LINE
                           write(*,fmt="(A)")     '        Usage: &WTD "EXP" tau '
                           write(*,fmt="(A)")     ' Process to be stopped'
                           stop
                        end if
                        Walker%JumpType(IJ)%WTD_TAU = Drstr(STRNUMB(1))   
                   case ("POW")
                        Walker%JumpType(IJ)%WTD_TYPE = CP_WALKTYPE_WTD_POW
                        call Extract_Numb(STR,2,N,STRNUMB)
                        if(N .lt. 2) then
                           write(*,fmt="(A,BZI6)") 'MCPSCU Error: parameter TAU missed at line', LINE
                           write(*,fmt="(A)")     '        Usage: &WTD "POW" tau, alpha '
                           write(*,fmt="(A)")     ' Process to be stopped'
                           stop
                        end if
                        Walker%JumpType(IJ)%WTD_TAU   = Drstr(STRNUMB(1))   
                        Walker%JumpType(IJ)%WTD_Alpha = Drstr(STRNUMB(2))     
                    end select     

               case("&JPVECS")
                    call Extract_Numb(STR,1,N,STRNUMB)
                    JN = ISTR(STRNUMB(1)) 
                    do J=1, JN
                       call GetInputStrLine(hFile,STR, LINE, "!", *100)
                       call Extract_Numb(STR,3,N,STRNUMB)
                       if(N .lt. 3) then
                         write(*,fmt="(A,BZI6)") 'MCPSCU Error: error on reading jump vector', LINE
                         write(*,fmt="(A)")     ' Process to be stopped'
                         stop
                       end if

                       Walker%JumpType(IJ)%JP_VEC(1,J) = Drstr(STRNUMB(1)) 
                       Walker%JumpType(IJ)%JP_VEC(2,J) = Drstr(STRNUMB(2)) 
                       Walker%JumpType(IJ)%JP_VEC(3,J) = Drstr(STRNUMB(3)) 
                    end do
                    Walker%JumpType(IJ)%JP_NVEC = JN

               case( "&ENDJUMP")
                    exit

             case default
                  write(*,fmt="(A, BZI6)") "MCPSCU warning: unknown keyword "//KEYWORD(1:LEN_TRIM(KEYWORD)), LINE
                  call ONWARNING(ERR)
             end select
         end do
        if(present(LINE0)) LINE0 = LINE
         return

 100     continue
         write(*,fmt="(A, BZI6)") "MCPSCU Error: fail to read Walk parameter "
         write(*,fmt="(A, BZI6)") "              Check if the keyword &ENDWALKDEF is missed"
         stop
  end subroutine Load_JumpBase0
  !****************************************************************************   

  !****************************************************************************
  subroutine Load_WalkerBase1(hFile, Walker, LINE0)
  !***  PURPOSE:   to load the describing parameters of randwalk controling parameters
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     SimBox
  implicit none
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(WalkerBase),   intent(inout) ::Walker
     integer,optional,    intent(inout) ::LINE0

     !--- local variables
      character*256::STR
      character*32::STRNUMB(4), KEYWORD
      integer:: N, LINE, ERR=0, JN, IJ

         if(present(LINE0)) then
            LINE = LINE0
          else
            LINE = 0
          end if

          IJ = 0
          do while(.TRUE.)
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             STR = adjustl(STR)
             call GetKeyWord("&", STR, KEYWORD)
             call UpCase(KEYWORD)
             select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
             case("&SYMB")
                   call Extract_substr(STR,1,N,STRNUMB) 
                   if(N .ge. 1) Walker%Name = STRNUMB(1)

             case("&WALKTYPE")
                  call Extract_substr(STR,1,N,STRNUMB) 
                  call UpCase(STRNUMB(1))
                  select case(STRNUMB(1))
                         case("1D")
                              Walker%WalkType = CP_WALKTYPE_DIM_1D     
                         case("2D")
                              Walker%WalkType = CP_WALKTYPE_DIM_2D     
                         case("3D")
                              Walker%WalkType = CP_WALKTYPE_DIM_3D     
                  end select
 
             case("&JUMP")
                   IJ = IJ + 1
                   call Load_JumpBase0(hFile, Walker, IJ, LINE) 
                   Walker%NJump = IJ
             case( "&ENDWALKDEF")
                    exit

             case default
                  write(*,fmt="(A, BZI6)") "MCPSCU Error: unknown keyword "//KEYWORD(1:LEN_TRIM(KEYWORD)), LINE
                  !call ONWARNING(ERR)
                  stop
             end select
         end do
        if(present(LINE0)) LINE0 = LINE
         return

 100     continue
         write(*,fmt="(A, BZI6)") "MCPSCU Error: fail to read Walk parameter "
         write(*,fmt="(A, BZI6)") "              Check if the keyword &ENDWALKDEF is missed"
         stop
  end subroutine Load_WalkerBase1
  !****************************************************************************  

  !****************************************************************************
  subroutine Load_WalkerBase2(Fname, Walker)
  !***  PURPOSE:   to load the describing parameters of randwalk controling parameters
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     SimBox
   implicit none
      !--- dummy varioables
      character*(*),     intent(in)    ::Fname
      class(WalkerBase), intent(inout) ::Walker
 
      !--- local variables
       integer::hFile, LINE
 
           call AvailableIOUnit(hFile) 
           open(UNIT=hFile, FILE=Fname, status='old')
           call Load_WalkerBase1(hFile, Walker, LINE)
           close(hFile)
           return
  end subroutine Load_WalkerBase2            
  !****************************************************************************

  !****************************************************************************
  subroutine Load_WalkerBase3(hFile, Walker, Line0)
  !***  PURPOSE:   to load the describing parameters of randwalk controling parameters
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     SimBox
   implicit none
      !--- dummy varioables
      integer,             intent(in)    ::hFile
      class(WalkerBase),   intent(inout) ::Walker(:)
      integer,optional,    intent(inout) ::Line0
 
      !--- local variables
      character*256::STR
      character*32::KEYWORD
      integer:: I, N, LINE, ERR=0
      save ERR

          if(present(LINE0)) then
            LINE = LINE0
          else
            LINE = 0
          end if
           
          N = 0
          do I=1, size(Walker)
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             STR = adjustl(STR)
             call GetKeyWord("&", STR, KEYWORD)
             call UpCase(KEYWORD)
             if(KEYWORD .eq. "&WALKER") then
                call Load_WalkerBase1(hFile, Walker(I), LINE)
                N = N + 1     
             end if
          end do

   100    if(present(LINE0)) Line0 = LINE

          if(N .eq. 0) then
             write(*,fmt="(A, BZI6)") "MCPSCU warning: no walker description is loaded "
             call ONWARNING(ERR)
          end if 
          return
  end subroutine Load_WalkerBase3            
  !****************************************************************************
  
  !****************************************************************************
  subroutine Load_WalkerBase4(Fname, Walker)
  !***  PURPOSE:   to load the describing parameters of randwalk controling parameters
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     SimBox
   implicit none
      !--- dummy varioables
      character*(*),          intent(in)    ::Fname
      class(WalkerBase), intent(inout) ::Walker(:)
 
      !--- local variables
       integer::hFile, LINE
 
           call AvailableIOUnit(hFile) 
           open(UNIT=hFile, FILE=Fname, status='old')
              call Load_WalkerBase3(hFile, Walker, LINE)
           close(hFile)   
           return
  end subroutine Load_WalkerBase4            
  !****************************************************************************  

 end module  MC_TypeDef_RandWalker_Base 
 !********************************************************************************

