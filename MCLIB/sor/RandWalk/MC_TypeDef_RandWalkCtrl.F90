
 module MC_TypeDef_RandWalkCtrl
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
     integer(KINDINT),   parameter::CP_MC_BOX_INF             = 0                     ! without boundary
     integer(KINDINT),   parameter::CP_MC_BOX_1D_SURF         = 1                     ! semi-infinite media with one surface
     integer(KINDINT),   parameter::CP_MC_BOX_1D_DSURF        = 3                     ! semi-infinite media with with two surface
     integer(KINDINT),   parameter::CP_MC_BOX_2D_SURF         = 4                     ! semi-infinite media with a rectangle surface, PD will be applied in X-Y
     integer(KINDINT),   parameter::CP_MC_BOX_2D_DSURF        = 5                     ! semi-infinite media with a rectangle surface, PD will be applied in X-Y
     integer(KINDINT),   parameter::CP_MC_BOX_3D              = 6                     ! semi-infinite media with a rectangle surface, PD will be applied in X-Y

     integer(KINDINT),   parameter::CP_MC_RAD_NONE            = 0
     integer(KINDINT),   parameter::CP_MC_RAD_PREEXIST        = 1
     integer(KINDINT),   parameter::CP_MC_RAD_EXTERNAL        = 2
     integer(KINDINT),   parameter::CP_MC_RAD_INTERNAL        = 3
     integer(KINDINT),   parameter::CP_MC_RAD_DELTA           = 2**16
     integer(KINDINT),   parameter::CP_MC_RAD_UNIFORM         = 2**17

     integer(KINDINT),   parameter::CP_MC_INI_DELTA           = 1
     integer(KINDINT),   parameter::CP_MC_INI_UNIFORM         = 2

     integer(KINDINT),   parameter::CP_MC_EVENT_CHECK_NONE      = 0                  ! no event check to be performed
     integer(KINDINT),   parameter::CP_MC_EVENT_CHECK_BYDEN     = 1                  ! frequence of performing event checks determined by the average density 
     integer(KINDINT),   parameter::CP_MC_EVENT_CHECK_BYDEN_S   = 2                  ! frequence of performing event checks determined by the density in a give sphere
     integer(KINDINT),   parameter::CP_MC_EVENT_CHECK_BYDEN_L   = 3                  ! frequence of performing event checks determined by the density in a give plan layer

     integer(KINDINT),   parameter::CP_WALK_RECTIMEINC_LINEAR  = 0
     integer(KINDINT),   parameter::CP_WALK_RECTIMEINC_LOG10   = 1

     type RandWalkCtrl
          !--
          integer                ::Restart  = 0           
          !---  fname of loading this control data
          character*256          ::MyFname
          !--- parameter concerning boundary
          character(len=2)       ::LU            =  "NM"                               !symbol to indicate the lenth unit of box
          integer                ::NBox          =  1                                  !number of parallel boxes
          integer                ::MXNPABOX      =  0                                  !max number of particle in a box
          integer                ::MXTOTALNP     =  10**7                            !max permitted total number of particles, NBox*MXNPABOX should <= MXTOTALNP 

          integer                ::Boxtype       =  CP_MC_BOX_2D_SURF                  !symbol to indicate the lenth unit of box
          integer                ::IfPD(3)       =  (/1,1,0/)                          !symbol to indicate the lenth unit of box
          real(KINDDF)           ::Boxsize(3)    =  (/10.D0,10.D0,1.D64/)              !the low boundary of periodic box
          real(KINDDF)           ::Boxlow(3)     =  (/-5.D0,-5.D0, 0.D0/)              !the low boundary of periodic box
          real(KINDDF)           ::Boxup(3)      =  (/ 5.D0, 5.D0, 1.D64/)             !the up boundary of  periodic box

          !--- parameter concerning insert particles
          integer                ::Radtype(CP_MX_RADSTYLE)        = CP_MC_RAD_NONE     !type controling insert new particle into box
          integer                ::RadSorNum                      = 0                  !number of radiation sorce
          real(KINDDF)           ::RadFlux(CP_MX_RADSTYLE)        = 0.D0               !RadFlux of irradition, 1/cm^2/s
          real(KINDDF)           ::RadFluxOn(CP_MX_RADSTYLE)      = 0.D0               !flux of irradition, in seconds, =0, without irradiation
          real(KINDDF)           ::RadFluxOff(CP_MX_RADSTYLE)     = 0.D0               !flux of irradition, in seconds, =0, without irradiation
          real(KINDDF)           ::RadPos0(CP_MX_RADSTYLE,3)      = -1.D64             !position of irradition in nano meter
          real(KINDDF)           ::RadPos1(CP_MX_RADSTYLE,3)      = -1.D64             !position of irradition in nano meter
          real(KINDDF)           ::RadVol(6, CP_MX_RADSTYLE)      = 0                  !volume of irradition in nano meter
          
          integer                ::RadtypeForID(CP_MX_RADSTYLE)   = 0                  !corresponding ID of walker 
          character*8            ::RadtypeForSymb(CP_MX_RADSTYLE) = ""                 !corresponding symbol of walker
          integer                ::RadNPPerstep                   = 0                  !number of particles to be inserted in one step

          !---- parameter concerning event check
          integer                ::EventCK_By                     = CP_MC_EVENT_CHECK_BYDEN_L  !this default value corresponding to Boxtype =  CP_MC_BOX_2D_SURF 

          !--- parameters concerning recording   
          integer               ::RectimeStepType = CP_WALK_RECTIMEINC_LOG10          !type of time points for recording
          integer               ::RectimeStepNum  = 10                                !it RectimeStepType is CP_WALK_RECTIMEINC_LOG10, number of
                                                                                      !time points in the magnitude order
          real(KINDDF)          ::RectimeStart    = 0.01D0                            !starting time of recording for RectimeStepType = CP_WALK_RECTIMEINC_LOG10
          real(KINDDF)          ::RectimeFixStep  = 1.D0                              !time interval of recording for RectimeStepType = CP_WALK_RECTIMEINC_LINEAR
          real(KINDDF)          ::RectimeEnd      = 0.01D0                            !total recording time, may also the total evolution time of simulations(in ps) 
          
          integer               ::RecCurSteps     = 0
          real(KINDDF)          ::RecCurtime      = -1.D0        
          character(len=256)    ::RecPath                =""                          !path for recording 
          integer               ::RecDepthBin     = 100                               !number of bins for recording depth distribution
          integer               ::RecConfig       = 1                                 !indicate if cofigures to be recorded

          !--- random number seed
          integer::SEED(4) = -1                                                        !the seed for random numner, <0, to use time as a seed

     end type RandWalkCtrl
     !--- 

  !---------------------------------------------------
     public :: Copy_RandWalkCtrl
  !---------------------------------------------------
     public :: GetMXWalkerNum_RandWalkCtrl 

  !---------------------------------------------------
     private:: GetRadNPart_RandWalkCtrl0, &
               GetRadNPart_RandWalkCtrl1
     public::  GetRadNPart_RandWalkCtrl
     interface GetRadNPart_RandWalkCtrl
        module procedure GetRadNPart_RandWalkCtrl0
        module procedure GetRadNPart_RandWalkCtrl1
     end interface GetRadNPart_RandWalkCtrl  

 !---------------------------------------------------
 !   routines concering recording
     public :: GetRecpath_RandWalkCtrl

 !---------------------------------------------------
     private:: ToNext_Rectime_RandWalkCtrl
     public :: GetRectime_RandWalkCtrl 
     interface GetRectime_RandWalkCtrl
        module procedure ToNext_Rectime_RandWalkCtrl
     end interface GetRectime_RandWalkCtrl  

  !---------------------------------------------------
     public::  GetRadRate_RandWalkCtrl     
  !---------------------------------------------------
     private:: GetRadTimestep_RandWalkCtrl0,   &
               GetRadTimestep_RandWalkCtrl1
     public :: GetRadTimestep_RandWalkCtrl
     interface GetRadTimestep_RandWalkCtrl
           module procedure GetRadTimestep_RandWalkCtrl0
           module procedure GetRadTimestep_RandWalkCtrl1
     end interface GetRadTimestep_RandWalkCtrl      
  !---------------------------------------------------
     private:: Load_RandWalkCtrl0
     public :: LoadSetup_RandWalkCtrl

  !---------------------------------------------------
     private:: Load_Box_INF_RandWalkCtrl,    &
               Load_Box_1DDSF_RandWalkCtrl,  &
               Load_Box_2DSF_RandWalkCtrl,   & 
               Load_Box_2DDSF_RandWalkCtrl,  &  
               Load_Box3D_RandWalkCtrl           
  !---------------------------------------------------
     private:: Load_Rad_2DSURF_RandWalkCtrl          
     public::  Load_RadCond_RandWalkCtrl          
               
  !---------------------------------------------------
     private:: Load_Rectime_RandWalkCtrl0, &
               Load_Recstyle_RandWalkCtrl

     public::  Load_Rectime_RandWalkCtrl          
     interface Load_Rectime_RandWalkCtrl 
        module procedure Load_Rectime_RandWalkCtrl0
          module procedure Load_Recstyle_RandWalkCtrl
     end interface Load_Rectime_RandWalkCtrl  

     public::  Load_Recstyle_RandWalkCtrl          
     interface Load_Recstyle_RandWalkCtrl 
          module procedure Load_Recstyle_RandWalkCtrl
     end interface Load_Recstyle_RandWalkCtrl  

  !---------------------------------------------------
     private:: Reset_Rectime_RandWalkCtrl0
     public :: Reset_Rectime_RandWalkCtrl 
     interface Reset_Rectime_RandWalkCtrl
        module procedure Reset_Rectime_RandWalkCtrl0
     end interface Reset_Rectime_RandWalkCtrl    

  !---------------------------------------------------
     private:: Transfer_RandWalkCtrl0,          &
               Transfer_RandWalkCtrl1,          &
               Transfer_RandWalkCtrl2
     public :: Transfer_RandWalkCtrl
     interface Transfer_RandWalkCtrl
        module procedure Transfer_RandWalkCtrl0
          module procedure Transfer_RandWalkCtrl1
            module procedure Transfer_RandWalkCtrl2
     end interface Transfer_RandWalkCtrl    


     


  !---------------------------------------------------
 
 contains
 !****************************************************************************

 !****************************************************************************
 subroutine LoadSetup_RandWalkCtrl(Fname, CtrlParam)
 !***  PURPOSE:   to load the describing parameters of randwalk controling parameters
 !
 !     INPUT:     hFile,  I/O unit number
 !
 !     OUTPUT     SimBox
   implicit none
      !--- dummy varioables
      character*(*),        intent(in)    ::Fname
      class(RandWalkCtrl),  intent(inout) ::CtrlParam
 
      !--- local variables
       integer::hFile, LINE
 
           call AvailableIOUnit(hFile) 
           open(UNIT=hFile, FILE=Fname, status='old')
           LINE = 0
           call Load_RandWalkCtrl0(hFile, CtrlParam, LINE)
           
           rewind(hFile)
           LINE = 0
           call Load_RadCond_RandWalkCtrl(hFile, CtrlParam, LINE)

           close(hFile) 
           !--- keep the filename of this control data, and will copy it to the output path
           CtrlParam%MyFname = Fname
           return
 end subroutine LoadSetup_RandWalkCtrl            
 !****************************************************************************

  !****************************************************************************
  subroutine Copy_RandWalkCtrl(From, To)
  !***  PURPOSE:   to make a copy of the control parameters
  !
  !     INPUT:     From
  !
  !     OUTPUT     To
  implicit none
     !--- dummy varioables
     class(RandWalkCtrl),  intent(in)  ::From
     class(RandWalkCtrl),  intent(out) ::To

     !--- local variables
          
          To%Restart       =  From%Restart
          To%MyFname       =  From%MyFname
          !--- parameter concerning boundary
          To%LU            =  From%LU      
          To%NBox          =  From%NBox    
          To%MXNPABOX      =  From%MXNPABOX
          To%MXTOTALNP     =  From%MXTOTALNP

          To%Boxtype       =  From%Boxtype
          To%IfPD          =  From%IfPD
          To%Boxsize       =  From%Boxsize
          To%Boxlow        =  From%Boxlow
          To%Boxup         =  From%Boxup

     !--- parameter concerning insert particles
          To%Radtype       = From%Radtype
          To%RadSorNum     = From%RadSorNum
          To%RadFlux       = From%RadFlux
          To%RadFluxOn     = From%RadFluxOn
          To%RadFluxOff    = From%RadFluxOff
          To%RadPos0       = From%RadPos0
          To%RadPos1       = From%RadPos1
          To%RadVol        = From%RadVol
     
          To%RadtypeForID  = From%RadtypeForID
          To%RadtypeForSymb= From%RadtypeForSymb
          To%RadNPPerstep  = From%RadNPPerstep

     !---- parameter concerning event check
          To%EventCK_By    = From%EventCK_By

     !--- parameters concerning recording   
          To%RectimeStepType= From%RectimeStepType
          To%RectimeStepNum = From%RectimeStepNum
          To%RectimeStart   = From%RectimeStart
          To%RectimeFixStep = From%RectimeFixStep
          To%RectimeEnd     = From%RectimeEnd
     
          To%RecCurSteps    = From%RecCurSteps
          To%RecCurtime     = From%RecCurtime
          To%RecPath        = From%RecPath
          To%RecDepthBin    = From%RecDepthBin
          To%RecConfig      = From%RecConfig

     !--- random number seed
          To%SEED           = From%SEED
      return
  end subroutine Copy_RandWalkCtrl 
 !****************************************************************************

 !****************************************************************************
  subroutine Load_RandWalkCtrl0(hFile, CtrlParam, LINE0)
  !***  PURPOSE:   to load the describing parameters of randwalk controling parameters
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     SimBox
  implicit none
     !--- dummy varioables
     integer,              intent(in)    ::hFile
     class(RandWalkCtrl),  intent(inout) ::CtrlParam
     integer,optional,     intent(inout) ::LINE0

     !--- local variables
      character*256::STR
      character*32::STRNUMB(5), KEYWORD
      integer::I, N, LINE, ERR=1, NBox

          if(present(LINE0)) then
             LINE = LINE0
          else
             LINE = 0
          end if

          NBox = 0
          do while(.TRUE.)
             call GetInputStrLine(hFile,STR, LINE, "!", *100)
             STR = adjustl(STR)
             call GetKeyWord("&", STR, KEYWORD)
             call UpCase(KEYWORD)
             select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                    case("&MCBSUBCTL", "&BOXSUBCTL")
                         ! the box number in one run 
                         call Extract_Numb(STR,1,N,STRNUMB)
                         if(N .ge. 1) then
                           NBox = IStr(STRNUMB(1))
                         end if
                         NBox = MAX(1, NBox)
                         !--- to make the number of boxes to be 2^N
                         N = 0
                         do while(2**N .lt. NBox)
                            N = N + 1
                         end do
                         CtrlParam%NBox = 2**N 

                        !---the box type 
                         call Extract_substr(STR,1,N,STRNUMB)
                         if(len_trim(STRNUMB(1)) .le. 0) then
                               write(*,fmt="(A, BZI6)") 'MCPSCU Error: the type of box is missed at line', LINE
                               write(*,fmt="(A, BZI6)") '       Usage: &BOX type'
                               write(*,fmt="(A, BZI6)") '              where type is one of: "INF", "SF1D", "DSF1D","SF2D","DSF2D","3DBOX"'
                               write(*,fmt="(' Process to be stopped')")
                               stop
                         end if
                         select case(STRNUMB(1))
                         case ("INF")
                               CtrlParam%Boxtype = CP_MC_BOX_INF
                               call Load_Box_INF_RandWalkCtrl(hFile, CtrlParam, LINE)
                         case ("SF1D")
                               CtrlParam%Boxtype = CP_MC_BOX_1D_SURF
                               call Load_Box_1DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
                         case ("DSF1D") 
                               CtrlParam%Boxtype = CP_MC_BOX_1D_DSURF
                               call Load_Box_1DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
                         case ("SF2D")
                               CtrlParam%Boxtype = CP_MC_BOX_2D_SURF
                               call Load_Box_2DSF_RandWalkCtrl(hFile, CtrlParam, LINE)
                         case ("DSF2D")
                               CtrlParam%Boxtype = CP_MC_BOX_2D_DSURF 
                               call Load_Box_2DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
                         case ("3DBOX")
                               CtrlParam%Boxtype = CP_MC_BOX_3D       
                               call Load_Box3D_RandWalkCtrl(hFile, CtrlParam, LINE)
                         end select 

                    case ("&RECSUBCTL")   
                         call Load_Recstyle_RandWalkCtrl(hFile, CtrlParam, LINE)  

                end select
         end do

 100     if(present(LINE0)) LINE0 = LINE
         !--- make the input consistent
         if(CtrlParam%NBOX * CtrlParam%MXNPABOX < CtrlParam%MXTOTALNP ) then
            CtrlParam%MXTOTALNP = CtrlParam%NBOX * CtrlParam%MXNPABOX
         else
            CtrlParam%MXNPABOX  = CtrlParam%MXTOTALNP / CtrlParam%NBOX
         end if    
         return
  
  end subroutine Load_RandWalkCtrl0
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Box_INF_RandWalkCtrl(hFile, CtrlParam, LINE)
    !***  PURPOSE:   to load the control parameter for box control
    !
    !     INPUT:     hFile,  I/O unit number
    !                STR,    a string line
    !
    !     OUTPUT     CtrlParam
     implicit none
       !--- dummy varioables
       integer,               intent(in)    ::hFile
       class(RandWalkCtrl),   intent(inout) ::CtrlParam
       integer,               intent(inout) ::LINE
       !--- local variables
       character*256::STR
       character*32::STRNUMB(4),KEYWORD
       integer::I,N
 
       !----
  
             !*** start the press step control controlling parametgers
             CtrlParam%IfPD       = 0
             CtrlParam%Boxsize    =  2.D64
             CtrlParam%Boxlow     = -1.D64
             CtrlParam%Boxup      =  1.D64
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!",*100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                       case( "&ENDSUBCTL")
                            exit
                       case default
                            write(*,*)"MCSCU Error: unknown keyword for RandWalkCtrl ", KEYWORD(1:LEN_TRIM(KEYWORD))
                            write(*,fmt="('               check control file at line:', BZI6)") LINE
                            stop
                       case ("&RANDSEED")
                          !$$*** To get seed for random number
                           call Extract_Numb(STR,4,N,STRNUMB)
                           do I=1, N
                              CtrlParam%SEED(I) = ISTR(STRNUMB(I))
                           end do
                       case ("&NATOM")
                          !$$*** To get the total number of particles
                           call Extract_Numb(STR,1,N,STRNUMB)
                           CtrlParam%MXNPABOX = ISTR(STRNUMB(1))

                end select
             end do
  
             return
   !----------------------------------------------------------------------------
    100    write(*,*)"MCPSCU Error in reading box control parameters."
           write(*,*)"The process to be stopped."
           write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
           stop
    end subroutine Load_Box_INF_RandWalkCtrl
  !****************************************************************************   

  !****************************************************************************
  subroutine Load_Box_1DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
  !***  PURPOSE:   to load the control parameter for box control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
   use MiniUtilities
   implicit none
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(RandWalkCtrl), intent(inout) ::CtrlParam
     integer,             intent(inout) ::LINE
     !--- local variables
      character*256::STR
      character*32::STRNUMB(4),KEYWORD
      integer::I, N
     !----
           call Load_Box_2DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
           CtrlParam%IfPD         =  (/0, 0, 1/)
           CtrlParam%Boxsize(1:2)  = 2.D64
           CtrlParam%Boxlow(1:2)  = -0.5D0*CtrlParam%Boxsize(1:2)
           CtrlParam%Boxup(1:2)   =  0.5D0*CtrlParam%Boxsize(1:2)

  end subroutine Load_Box_1DDSF_RandWalkCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Box_2DSF_RandWalkCtrl(hFile, CtrlParam, LINE)
  !***  PURPOSE:   to load the control parameter for box control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
   use MiniUtilities
   implicit none
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(RandWalkCtrl), intent(inout) ::CtrlParam
     integer,             intent(inout) ::LINE
     !--- local variables
     !----
      call Load_Box_2DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
  end subroutine Load_Box_2DSF_RandWalkCtrl
  !****************************************************************************  

  !****************************************************************************
  subroutine Load_Box_2DDSF_RandWalkCtrl(hFile, CtrlParam, LINE)
  !***  PURPOSE:   to load the control parameter for box control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
   implicit none
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(RandWalkCtrl), intent(inout) ::CtrlParam
     integer,             intent(inout) ::LINE
     !--- local variables
      character*256::STR
      character*32::STRNUMB(3),KEYWORD
      integer::I,N
     !----

           CtrlParam%IfPD    = (/1,1,0/)
           CtrlParam%Boxlow  =  (/-0.5D0, -0.5D0,   0.D00/)
           CtrlParam%Boxup   =  (/ 0.5D0,  0.5D0,   1.D64/)
           CtrlParam%Boxsize =  CtrlParam%Boxup - CtrlParam%Boxlow

           !*** start the press step control controlling parametgers
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MCSCU Error: unknown keyword for RandWalkCtrl ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop

                     case ("&SIZE")
                          !$$*** To get box size
                           call Extract_Numb(STR,2,N,STRNUMB)
                           if(n .lt. 1) then
                              write(*,fmt="(' MCPSCU Error: the boxsize should be set in x,y directions')")
                              write(*,fmt="('               check box file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage: &SIZE bx  ,by ')")
                              write(*,fmt="(' Process to be stopped')")
                              stop
                            end if
                            if(N .eq. 1) then
                               CtrlParam%Boxsize(1) = DRSTR(STRNUMB(1))
                               CtrlParam%Boxsize(2) = DRSTR(STRNUMB(1))
                            else
                              CtrlParam%Boxsize(1) = DRSTR(STRNUMB(1))
                              CtrlParam%Boxsize(2) = DRSTR(STRNUMB(2))
                            end if
                            CtrlParam%Boxlow(1)  = -0.5D0*CtrlParam%Boxsize(1)
                            CtrlParam%Boxlow(2)  = -0.5D0*CtrlParam%Boxsize(2)
                            CtrlParam%Boxup(1)   =  CtrlParam%Boxlow(1) + CtrlParam%Boxsize(1)
                            CtrlParam%Boxup(2)   =  CtrlParam%Boxlow(2) + CtrlParam%Boxsize(2)

                     case ("&THICK")
                          !$$*** To get box size
                           call Extract_Numb(STR,1,N,STRNUMB)
                           if(n .lt. 1) then
                              write(*,fmt="(' MCPSCU Error: the box thickness should be set in z directions')")
                              write(*,fmt="('               check box file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage: &THICK thick ')")
                              write(*,fmt="(' Process to be stopped')")
                              stop
                            end if
                            CtrlParam%Boxup(3) = CtrlParam%Boxlow(3) + DRSTR(STRNUMB(1))

                     case ("&RANDSEED")
                          !$$*** To get seed for random number
                           call Extract_Numb(STR,4,N,STRNUMB)
                           do I=1, N
                              CtrlParam%SEED(I) = ISTR(STRNUMB(I))
                           end do

                     case ("&NATOM")
                           !$$*** To get the total number of particles
                           call Extract_Numb(STR,1,N,STRNUMB)
                           CtrlParam%MXNPABOX = ISTR(STRNUMB(1))

              end select
           end do
           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MCPSCU Error in reading box control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_Box_2DDSF_RandWalkCtrl
  !****************************************************************************  

  !****************************************************************************
  subroutine Load_Box3D_RandWalkCtrl(hFile, CtrlParam, LINE)
  !***  PURPOSE:   to load the control parameter for box control
  !
  !     INPUT:     hFile,  I/O unit number
  !                STR,    a string line
  !
  !     OUTPUT     CtrlParam
   use MiniUtilities
   implicit none
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(RandWalkCtrl), intent(inout) ::CtrlParam
     integer,             intent(inout) ::LINE
     !--- local variables
      character*256::STR
      character*32::STRNUMB(4),KEYWORD
      integer::I, N
     !----

           !*** start the press step control controlling parametgers
           CtrlParam%IfPD    = 1
           CtrlParam%Boxlow  =  (/-0.5D0, -0.5D0,  -0.5D0/)
           CtrlParam%Boxup   =  (/ 0.5D0,  0.5D0,   0.5D0/)
           CtrlParam%Boxsize =  CtrlParam%Boxup - CtrlParam%Boxlow
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MCSCU Error: unknown keyword for RandWalkCtrl ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop
                          
                     case ("&SIZE")
                          !$$*** To get box size
                           call Extract_Numb(STR,3,N,STRNUMB)
                           if(n .lt. 3) then
                              write(*,fmt="(' MCPSCU Error: the boxsize should be set in x,y,z directions')")
                              write(*,fmt="('               check box file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage: &SIZE bx  ,by , bz ')")
                              write(*,fmt="(' Process to be stopped')")
                              stop
                            end if
                            do I=1, 3
                              CtrlParam%Boxsize(I) = DRSTR(STRNUMB(I))
                            end do
                            CtrlParam%Boxlow = -0.5*CtrlParam%Boxsize

                     case ("&LOWB","&BOXLOW")
                          !$$*** To get box size
                           call Extract_Numb(STR,3,n,STRNUMB)
                           if(N .lt. 3) then
                              write(*,fmt="(' MCPSCU Error: the low-boundary of a box should be set in x,y,z directions')")
                              write(*,fmt="('               check box file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage: &LOWB xmi=  ,ymi= , zmi= ')")
                              write(*,fmt="(' Process to be stopped')")
                              stop
                            end if
                            do I=1, 3
                              CtrlParam%Boxlow(I) = DRSTR(STRNUMB(I))
                            end do

                     case ("&PERIDIC")
                           !--- To get if periodic conditions will be used
                           call Extract_Numb(STR,3,N,STRNUMB)
                           if(N .lt. 3) then
                              write(*,fmt="('MCPSCU Error: peroidic conditions should be given for three direction')")
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              stop
                           end if
                           CtrlParam%IfPD(1) = ISTR(STRNUMB(1))
                           CtrlParam%IfPD(2) = ISTR(STRNUMB(2))
                           CtrlParam%IfPD(3) = ISTR(STRNUMB(3))

                     case ("&RANDSEED")
                          !$$*** To get seed for random number
                           call Extract_Numb(STR,4,N,STRNUMB)
                           do I=1, N
                              CtrlParam%SEED(I) = ISTR(STRNUMB(I))
                           end do

                     case ("&NATOM")
                           !$$*** To get the total number of particles
                           call Extract_Numb(STR,1,N,STRNUMB)
                           CtrlParam%MXNPABOX = ISTR(STRNUMB(1))

              end select
           end do
           CtrlParam%Boxup  =  CtrlParam%Boxlow + CtrlParam%Boxsize

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MCPSCU Error in reading box control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
  end subroutine Load_Box3D_RandWalkCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Rectime_RandWalkCtrl0(Str, CtrlParam)
  !***  PURPOSE:   to load recording setup of RandWalkCtrl
  !
  !     INPUT:     STR,    a string line
  !
  !     OUTPUT     CtrlParam
   use MiniUtilities
   implicit none
     !--- dummy varioables
     character*(*),       intent(in)    ::Str
     class(RandWalkCtrl), intent(inout) ::CtrlParam
     !--- local variables
      character*32::STRTMP(5)=""
      integer::N
     !----
           call extract_optstr(Str, "-","REC", 5, N, STRTMP)
           if(N .eq. 0) call extract_optstr(Str, "&","REC", 5, N, STRTMP)
           if(N .ge. 1) then
               call UpCase(STRTMP(1))
               select case(STRTMP(1))
               case ("LOG", "LOG10")
                    CtrlParam%RectimeStepType = CP_WALK_RECTIMEINC_LOG10
                    if(N .lt. 3) then
                      write(*,fmt="(' MCPSCU Error: recording parameter are missed')")
                      write(*,fmt="(' Usage: &REC type, start, end ')")
                      write(*,fmt="(' Process to be stopped')")
                      stop
                    end if
                    CtrlParam%RectimeStart = Drstr(STRTMP(2)) 
                    CtrlParam%RectimeEnd   = Drstr(STRTMP(3)) 
                    if(N .ge. 4) then
                       CtrlParam%RectimeStepNum  = Istr(STRTMP(4)) 
                    end if 

               case ("LIN", "LINEAR")
                    CtrlParam%RectimeStepType = CP_WALK_RECTIMEINC_LINEAR
                    if(N .lt. 4) then
                      write(*,fmt="(' MCPSCU Error: recording parameter are missed')")
                      write(*,fmt="(' Usage: &REC type, start, end, tstep ')")
                      write(*,fmt="(' Process to be stopped')")
                      stop
                    end if
                    CtrlParam%RectimeStart    = Drstr(STRTMP(2)) 
                    CtrlParam%RectimeEnd      = Drstr(STRTMP(3)) 
                    CtrlParam%RectimeFixStep  = Drstr(STRTMP(4)) 
               end select     
           end if    

           return
  end subroutine Load_Rectime_RandWalkCtrl0
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Recstyle_RandWalkCtrl(hFile, CtrlParam, Line)
   !***  PURPOSE:   to load recording setup of RandWalkCtrl
   !
   !     INPUT:     STR,    a string line
   !
   !     OUTPUT     CtrlParam
     use MiniUtilities
     implicit none
       !--- dummy varioables
       integer,             intent(in)    ::hFile
       class(RandWalkCtrl), intent(inout) ::CtrlParam
       integer,             intent(inout) :: Line
       !--- local variables
       character*256::STR
       character*256::STRTMP(1)=""
       character*32::STRNUMB(3),KEYWORD
       integer::I, N, M
        !----
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MCSCU Error: unknown keyword for RandWalkCtrl ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop
                          
                     case ("&RECTIME")
                          call Extract_Substr(STR,1, N, STRTMP)
                          call Extract_Numb(STR,  3, M, STRNUMB)
                          if(N .lt. 1 .or. M.lt.2) then
                              write(*,fmt="(A, BZI6)") ' MCPSCU Error: the recording style missed at line', Line
                              write(*,fmt="(A)")       '        Usage: &RECTIME style starttime, endtime, steps'
                              write(*,fmt="(A)")       '               style should be one of "LOG", "LINEAR" '
                              write(*,fmt="(A)")       ' Process to be stopped'
                              stop
                            end if
                            call UpCase(STRTMP(1))
                            select case(STRTMP(1))
                            case ("LOG", "LOG10")
                                 CtrlParam%RectimeStepType = CP_WALK_RECTIMEINC_LOG10
                            case ("LIN", "LINEAR")
                                 CtrlParam%RectimeStepType = CP_WALK_RECTIMEINC_LINEAR
                            end select     
                            CtrlParam%RectimeStart = Drstr(STRNUMB(1)) 
                            CtrlParam%RectimeEnd   = Drstr(STRNUMB(2)) 
                            if(M .ge. 3) then
                               CtrlParam%RectimeStepNum  = Istr(STRNUMB(3)) 
                            end if 

                     case ("&OUTPUT")
                            call Extract_Substr(STR,1, N, STRTMP)
                            if(N .ge. 1) then
                               CtrlParam%RecPath = STRTMP(1)
                            end if  

                     case ("&DEPTHBIN")
                            call Extract_Numb(STR,1, N, STRNUMB)
                            if(N .ge. 1) then
                               CtrlParam%RecDepthBin = IStr(STRNUMB(1))
                            end if  

                     case ("&OUTPUTCFG")
                          call Extract_Substr(STR,1, N, STRTMP)
                          if(N .ge. 1) then
                             call UpCase(STRTMP(1))
                             if(STRTMP(1) .eq. "NO") then
                                CtrlParam%RecConfig = 0
                             else if(STRTMP(1) .eq. "YES") then
                                CtrlParam%RecConfig = 1   
                             end if   
                          end if 

                  end select
           end do

           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MCPSCU Error in reading RandWalkCtrl parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop  
         return
    end subroutine Load_Recstyle_RandWalkCtrl
  !****************************************************************************

  !****************************************************************************
   subroutine Load_RadCond_RandWalkCtrl(hFile, CtrlParam, Line)
   !***  PURPOSE:   to load radiation condition
   !
   !     INPUT:     hFile,   I/O unit
   !                IR,  
   !                Line,
   !
   !     OUTPUT     CtrlParam
     use MiniUtilities
     !--- dummy varioables
     integer,             intent(in)    ::hFile
     class(RandWalkCtrl), intent(inout) ::CtrlParam
     integer,             intent(inout) ::Line
     !--- local variables
     character*256::STR
     character*32::STRNUMB(3),KEYWORD
     integer::IA, N, M, IFLAG


       !--- looking for radiation section
      do while(.TRUE.)
         call GetInputStrLine(hFile,STR, LINE, "!",*100)
         STR = adjustl(STR)
         call GetKeyWord("&", STR, KEYWORD)
         call UpCase(KEYWORD)
         if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq. "&RADSUBCTL") then
             exit
         end if   
      end do 
       
      IA = 0
      do while(.TRUE.)
         call GetInputStrLine(hFile,STR, LINE, "!",*100)
         STR = adjustl(STR)
         call GetKeyWord("&", STR, KEYWORD)
         call UpCase(KEYWORD)
         select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                case( "&ENDSUBCTL")
                     exit
                case default
                     write(*,*)"MCSCU Error: unknown keyword for radiation condition ", KEYWORD(1:LEN_TRIM(KEYWORD))
                     write(*,fmt="('               check control file at line:', BZI6)") LINE
                     stop
               case ("&MXNPATIME")
                    call Extract_Numb(STR,  1, N, STRNUMB)
                    if( N.ge.1) then
                        CtrlParam%RadNPPerstep = Istr(STRNUMB(1)) 
                    end if

                case ("&SORSUBCTL")
                     call Extract_substr(STR,  2, N, STRNUMB)
                     if( N .ge. 2) then
                       IA = IA + 1 
                       KEYWORD = STRNUMB(1)
                       call UpCase(KEYWORD)
                       if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq.  "EXT") then
                           IFLAG = 2
                           CtrlParam%Radtype(IA) = CP_MC_RAD_EXTERNAL
                       else 
                           KEYWORD = STRNUMB(2) 
                           call UpCase(KEYWORD) 
                           if(KEYWORD(1:LEN_TRIM(KEYWORD)) .eq.  "EXT") then
                              IFLAG = 1
                              CtrlParam%Radtype(IA) = CP_MC_RAD_EXTERNAL     
                           end if                          
                       end if

                       CtrlParam%RadtypeForSymb(IA) = STRNUMB(IFLAG)(1:len_trim(STRNUMB(IFLAG))) 
                       call Load_Rad_2DSURF_RandWalkCtrl(hFile, CtrlParam, IA, Line) 
                       CtrlParam%RadSorNum = IA
                     else
                        write(*,fmt="(A)")        "MCSCU Error: radiation source type is missed "
                        write(*,fmt="(A)")        "      Usage: &SORSUBCTL  sorname, sortype"
                        write(*,fmt="(A)")        "    Example: &SORSUBCTL  'He', 'ext'"
                        write(*,fmt="(A, BZI16)") '             check control file at line:', LINE
                        stop
                     end if
                     
         end select
      end do
      return

100   write(*,*)"MCSCU Error: Erro in read in radiation condition"
      stop

      return
  end subroutine Load_RadCond_RandWalkCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Rad_2DSURF_RandWalkCtrl(hFile, CtrlParam, IR, Line)
   !***  PURPOSE:   to load rradition style for IDth worker
   !
   !     INPUT:     hFile,   I/O unit
   !                IR,  
   !                Line,
   !
   !     OUTPUT     CtrlParam
     use MiniUtilities
     implicit none
       !--- dummy varioables
       integer,             intent(in)    ::hFile
       class(RandWalkCtrl), intent(inout) ::CtrlParam
       integer,             intent(in)    ::IR
       integer,             intent(inout) ::Line
       !--- local variables
       character*256::STR
       character*32::STRNUMB(3),KEYWORD
       integer::I, N, M

           CtrlParam%Radtype(IR) = CP_MC_RAD_NONE
           do while(.TRUE.)
              call GetInputStrLine(hFile,STR, LINE, "!",*100)
              STR = adjustl(STR)
              call GetKeyWord("&", STR, KEYWORD)
              call UpCase(KEYWORD)
              select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit
                     case default
                          write(*,*)"MCSCU Error: unknown keyword for RandWalkCtr ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          stop
                          
                     case ("&FLUX")
                          call Extract_Numb(STR,  1, M, STRNUMB)
                          if( M.ge.1) then
                            CtrlParam%RadFlux(IR)    = Drstr(STRNUMB(1)) 
                          end if
                          
                     case ("&TIMEON")
                          call Extract_Numb(STR,  1, M, STRNUMB)
                          if( M .ge. 1) then
                             CtrlParam%RadFluxOn(IR) = Drstr(STRNUMB(1)) 
                          end if   

                     case ("&TIMEOFF")
                          call Extract_Numb(STR,  1, M, STRNUMB)
                          if( M .ge.1 ) then
                             CtrlParam%RadFluxOff(IR) = Drstr(STRNUMB(1)) 
                          end if  
                          
                     case ("&DEPTH")
                          call Extract_Numb(STR, 2, M, STRNUMB)
                          if( M .gt. 1) then
                             CtrlParam%RadPos0(IR,:) = 0.D0
                             CtrlParam%RadPos1(IR,:) = 0.D0
                             CtrlParam%RadPos0(IR,3) = Drstr(STRNUMB(1)) 
                             CtrlParam%RadPos1(IR,3) = Drstr(STRNUMB(2)) 
                             CtrlParam%Radtype(IR)   = ior(CP_MC_RAD_UNIFORM,  CtrlParam%Radtype(IR))
                          else if( M .eq. 1) then
                             CtrlParam%RadPos0(IR,:) = 0.D0
                             CtrlParam%RadPos0(IR,3) = Drstr(STRNUMB(1)) 
                             CtrlParam%RadPos1(IR,:) = CtrlParam%RadPos0(IR,:)
                             CtrlParam%Radtype(IR)   = ior(CP_MC_RAD_DELTA,  CtrlParam%Radtype(IR))
                          end if  
                          
              end select
           end do

           CtrlParam%RadFlux    = CtrlParam%RadFlux/CP_SQRCM2NM/CP_S2PS  
           CtrlParam%RadFluxOn  = CtrlParam%RadFluxOn*CP_S2PS  
           CtrlParam%RadFluxOff = CtrlParam%RadFluxOff*CP_S2PS  
           return
 !----------------------------------------------------------------------------
  100    write(*,*)"MCPSCU Error in reading RandWalkCtrl parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop  
         return
    end subroutine Load_Rad_2DSURF_RandWalkCtrl
  !****************************************************************************    

  !****************************************************************************
  subroutine ToNext_Rectime_RandWalkCtrl(CtrlParam, Rectime, Flag)
   !***  PURPOSE:   to get current and next recording time
   !                 
   !     INPUT:     CtrlParam, 
   !
   !     OUTPUT:    Rectime
    use MiniUtilities
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl), intent(inout) ::CtrlParam
       real(KINDDF),        intent(out)   ::Rectime
       integer,             intent(out)   ::Flag
       !--- local variables
       integer::N, M, N10, NO, NN
       real(KINDDF)::TT
       !----
         Flag = 0 
         if(CtrlParam%RecCurtime .lt. 0.D0) then
             CtrlParam%RecCurtime    = CtrlParam%RectimeStart
             CtrlParam%RecCurSteps   = 0
             Flag = 1
          end if  
         !--- to update next time 
          CtrlParam%RecCurSteps = CtrlParam%RecCurSteps + 1 
          select case(CtrlParam%RectimeStepType)
          case(CP_WALK_RECTIMEINC_LOG10)
               NO = int(dlog10(CtrlParam%RectimeStart) + 0.0000001D0)
               if(NO .lt. 0) NO = NO-1
               N10 = CtrlParam%RectimeStepNum
               NN  = N10/10
               N   = CtrlParam%RecCurSteps
               M   = CtrlParam%RecCurSteps/N10
               N = N - M*N10 
               if(N .eq. 0) then
                  N    = N + NN
                  Flag = 1
                  CtrlParam%RecCurSteps = CtrlParam%RecCurSteps + NN
               end if   
               CtrlParam%RecCurtime = 10.D0**dble(NO+M)*(1.D0 + 10.D0*dble(N-1)/dble(N10))

          case(CP_WALK_RECTIMEINC_LINEAR)
               CtrlParam%RecCurtime = CtrlParam%RectimeStart  + CtrlParam%RecCurSteps*CtrlParam%RectimeFixStep
          end select 
          Rectime = CtrlParam%RecCurtime
         return
    end subroutine ToNext_Rectime_RandWalkCtrl
  !****************************************************************************    

  !****************************************************************************
  subroutine Reset_Rectime_RandWalkCtrl0(CtrlParam)
   !***  PURPOSE:   to get current and next recording time
   !                 
   !     INPUT:     CtrlParam, 
   !
   !     OUTPUT:    NextTime
   !     NOTE:      the cuurent and next time of  CtrlParam to be updated.
   !                see the diffrence from Get_CurRectime_RandWalkCtrl           
    use MiniUtilities
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl), intent(inout) ::CtrlParam
       !--- local variables
       !----
           CtrlParam%RecCurtime  = -1.D0
           CtrlParam%RecCurSteps = 0

         return
    end subroutine Reset_Rectime_RandWalkCtrl0
  !****************************************************************************    

  !****************************************************************************
  subroutine GetMXWalkerNum_RandWalkCtrl(CtrlParam, MXNP, NBOX)
  !***  PURPOSE:   to give an estimation of the number of walkers
  !                 
  !     INPUT:     CtrlParam, 
  !
  !     OUTPUT:    MXNP,  the max number of particle for each box
  !                NBOX,  the number of boxes
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl), intent(inout) ::CtrlParam
       integer,             intent(inout) ::MXNP, NBOX
       !--- local variables
       real(KINDDF)::NN, FLUXOFF
       real(KINDDF),parameter::RES=1.2
       integer::I
       !----
           select case(CtrlParam%Boxtype)
           case(CP_MC_BOX_INF, CP_MC_BOX_1D_SURF, CP_MC_BOX_1D_DSURF)
                MXNP = max(1024, CtrlParam%MXNPABOX)

           case (CP_MC_BOX_2D_SURF)
                NN = 0
                do I = 1, CtrlParam%RadSorNum
                   FLUXOFF = min(CtrlParam%RadFluxOff(I),CtrlParam%RectimeEnd)
                   NN      = NN + CtrlParam%RadFlux(I)*(FLUXOFF - CtrlParam%RadFluxOn(I))*CtrlParam%Boxsize(1)*CtrlParam%Boxsize(2)
                end do          
                MXNP = max(1024, nint(NN*RES))
                CtrlParam%MXNPABOX = min(MXNP,CtrlParam%MXNPABOX)
           case (CP_MC_BOX_2D_DSURF, CP_MC_BOX_3D)
                NN = 0
                do I = 1, CtrlParam%RadSorNum
                   FLUXOFF = min(CtrlParam%RadFluxOff(I),CtrlParam%RectimeEnd)
                   NN      = NN + CtrlParam%RadFlux(I)*(FLUXOFF - CtrlParam%RadFluxOn(I))*CtrlParam%Boxsize(1)*CtrlParam%Boxsize(2)
                end do          
                MXNP = max(1024, nint(NN*RES))
                CtrlParam%MXNPABOX = min(MXNP,CtrlParam%MXNPABOX)
           end select

           NBOX = CtrlParam%NBox
           MXNP = CtrlParam%MXNPABOX
           !--- check if we should reset the MXNPABOX
           if(NBOX*MXNP > CtrlParam%MXTOTALNP) then
              MXNP = CtrlParam%MXTOTALNP / NBOX 
              CtrlParam%MXNPABOX = MXNP
           end if 

           return
    end subroutine GetMXWalkerNum_RandWalkCtrl
  !****************************************************************************  
    
  !****************************************************************************
  subroutine GetRadRate_RandWalkCtrl(CtrlParam, Time, Rate)
  !***  PURPOSE:   to the rate of generating particle at given time in unit ps^-1
  !                 
  !     INPUT:     CtrlParam, 
  !                Time
  !     OUTPUT:    Rate
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl),   intent(in)    ::CtrlParam
       real(KINDDF),optional, intent(in)    ::Time
       real(KINDDF),          intent(inout) ::Rate
       !--- local variables
       integer::I
       !----
           select case(CtrlParam%Boxtype)
           case(CP_MC_BOX_INF, CP_MC_BOX_1D_SURF, CP_MC_BOX_1D_DSURF)
                Rate = 0.D0

           case (CP_MC_BOX_2D_SURF, CP_MC_BOX_2D_DSURF, CP_MC_BOX_3D)
                Rate = 0.D0
                if(present(Time)) then
                   do I = 1, CtrlParam%RadSorNum
                      if(CtrlParam%RadFluxon(I) <= Time .and. Time <= CtrlParam%RadFluxoff(I) ) then
                         Rate = Rate + CtrlParam%RadFlux(I)
                      end if   
                   end do          
                else   
                  do I = 1, CtrlParam%RadSorNum
                     Rate = Rate + CtrlParam%RadFlux(I)
                  end do          
                end if 
                Rate = Rate*CtrlParam%Boxsize(1)*CtrlParam%Boxsize(2)
           end select
         return
    end subroutine GetRadRate_RandWalkCtrl
  !****************************************************************************    

  !****************************************************************************
  subroutine GetRadNPart_RandWalkCtrl0(CtrlParam, NP)
  !***  PURPOSE:   to get the number of particles generated at given time and time interval
  !                 
  !     INPUT:     CtrlParam, 
  !                Time
  !                DeltaTime,
  ! 
  !     OUTPUT:    NP
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl),   intent(in)    ::CtrlParam
       integer,               intent(out)   ::NP(:)
       !--- local variables
       integer::I
       real(KINDDF)::RNP(CP_MX_RADSTYLE)
       !----
                RNP = 0.D0
                do I = 1, CtrlParam%RadSorNum
                   RNP(I) = RNP(I) + CtrlParam%RadFlux(I)
                end do          
                if(sum(RNP) .gt. 0) RNP = RNP/sum(RNP) 
                do I =1, CtrlParam%RadSorNum
                   NP(I) = RNP(I)*CtrlParam%RadNPPerstep
                end do     
         return
    end subroutine GetRadNPart_RandWalkCtrl0
  !****************************************************************************

  !****************************************************************************
  subroutine GetRadNPart_RandWalkCtrl1(CtrlParam, Time, DeltaTime, NP)
  !***  PURPOSE:   to get the number of particles generated at given time and time interval
  !                 
  !     INPUT:     CtrlParam, 
  !                Time
  !                DeltaTime,
  ! 
  !     OUTPUT:    NP
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl),   intent(in)    ::CtrlParam
       real(KINDDF),          intent(in)    ::Time, DeltaTime
       integer,               intent(out)   ::NP(:)
       !--- local variables
       integer::I
       !----
                NP = 0.D0
                do I = 1, CtrlParam%RadSorNum
                   if(CtrlParam%RadFluxon(I) <= Time .and. Time <= CtrlParam%RadFluxoff(I) ) then
                      NP(I) = NP(I) + nint(CtrlParam%RadFlux(I)*DeltaTime*CtrlParam%Boxsize(1)*CtrlParam%Boxsize(2))
                   end if   
                end do          
         return
    end subroutine GetRadNPart_RandWalkCtrl1
  !****************************************************************************    
    
  !****************************************************************************
  subroutine GetRadTimestep_RandWalkCtrl0(CtrlParam, Time, NP, DeltaTime)
  !***  PURPOSE:   to get the time for generating NP particles 
  !                 
  !     INPUT:     CtrlParam, 
  !                NP
  !                Time
  !  
  !     OUTPUT:    DelataTime
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl),   intent(in)    ::CtrlParam
       integer,               intent(in)    ::NP
       real(KINDDF),          intent(in)    ::Time
       real(KINDDF),          intent(out)   ::DeltaTime
       !--- local variables
       real(KINDDF)::RATE
       !----
          
          call GetRadRate_RandWalkCtrl(CtrlParam, Time, RATE)
          if(RATE .gt. 0.D0) then
            DeltaTime = dble(NP)/RATE
          else
            DeltaTime  =1.D108   
          end if
       
         return
    end subroutine GetRadTimestep_RandWalkCtrl0
  !****************************************************************************    

  !****************************************************************************
  subroutine GetRadTimestep_RandWalkCtrl1(CtrlParam, Time, DeltaTime)
  !***  PURPOSE:   to get the time for generating the RadNPPerstep particles 
  !                 
  !     INPUT:     CtrlParam, 
  !                Time
  !  
  !     OUTPUT:    DelataTime
     implicit none
       !--- dummy varioables
       class(RandWalkCtrl),   intent(in)    ::CtrlParam
       real(KINDDF),          intent(in)    ::Time
       real(KINDDF),          intent(out)   ::DeltaTime
       !--- local variables
       !----
          
          call GetRadTimestep_RandWalkCtrl0(CtrlParam, Time, CtrlParam%RadNPPerstep, DeltaTime)
       
         return
    end subroutine  GetRadTimestep_RandWalkCtrl1
  !****************************************************************************    

  !****************************************************************************
  subroutine GetRecpath_RandWalkCtrl(CtrlParam, Path)
  !***  PURPOSE:   to create and retuern path for recording stat of walkers
  !     INPUT :    CtrlParam
  !
  !     OUTRPUT:   path
  !    
  !    NOTE: the path is indicated by "\" in windows and "/" in linux
  !
  !     OUTPUT
  !---
      implicit none
      class(RandWalkCtrl),  intent(in) ::CtrlParam
      character*(*)                    ::Path
      !--- local variables
      character*256::swapPath
      integer I, ERR
      logical EX


             swapPath = CtrlParam%RecPath(1:len_trim(CtrlParam%RecPath))
             if(swapPath(2:2) .eq. ':') then
                inquire(FILE=swapPath(1:2), EXIST=EX)
                if(.not.EX) then
                   write(*,fmt="(A)") "'MCPSCU Error: Device "//swapPath(1:2)//" not found"
                   write(*,fmt="(A)") "              check the filename: ",swapPath(1:len_trim(path))
                   write(*,fmt="(A)") "              Process to be stopped"
                   stop
                end if
             end if

             do I=1, len_trim(swapPath)
                ! iachar('\') == 92
                ! iachar('/') == 47
                if(iachar(swapPath(I:I)) .eq. 92 .or. &
                   iachar(swapPath(I:I)) .eq. 47 .or. &
                   I.eq.len_trim(swapPath) ) then
                   ! $$--- the system dependence function
                   if(I.eq.len_trim(swapPath)) then
                      I = I + 1
                   end if   
                   #ifdef WIN_FILE
                   swapPath(I:I) = achar(92)
                   #else
                   swapPath(I:I) = achar(47)
                   #endif

                   inquire(FILE=swapPath(1:I-1), EXIST=EX)
                   if(.not. EX) then
                      write(*,fmt="(A,A)") ' MDPSCU Message: create data folder:', swapPath(1:I-1)
                      call SYSTEM("mkdir "//swapPath(1:I-1))
                      !call EXECUTE_COMMAND_LINE
                   end if
                end if
             enddo
             Path = swapPath(1:len_trim(swapPath))
             return
  end subroutine GetRecpath_RandWalkCtrl
  !**************************************************************************** 
  
  !****************************************************************************
  subroutine Transfer_RandWalkCtrl0(Fname, TgtPath)
  !***  PURPOSE:   to copy the original data to target path
  !     INPUT :    CtrlParam
  !                TgtPath     the targtet path  
  !
  !     OUTPUT
  !---
      implicit none
      character*(*),      intent(in) ::Fname
      character*(*),      intent(in) ::TgtPath
      !--- local variables

             #ifdef WIN_FILE
                 call SYSTEM("copy "//Fname(1:len_trim(Fname))//" "//TgtPath(1:len_trim(TgtPath))   )
             #else
                 call SYSTEM("cp   "//Fname(1:len_trim(Fname))//" "//TgtPath(1:len_trim(TgtPath))   )
               #endif
             return
  end subroutine Transfer_RandWalkCtrl0
  !****************************************************************************  

  !****************************************************************************
  subroutine Transfer_RandWalkCtrl1(CtrlParam, SorFname)
  !***  PURPOSE:   to copy the original data to target path give be CtrlParam
  !     INPUT :    CtrlParam
  !                TgtPath     the targtet path  
  !
  !     OUTPUT
  !---
      implicit none
      class(RandWalkCtrl),  intent(in) ::CtrlParam
      character*(*),        intent(in) ::SorFname
      !--- local variables
      character*256::swapPath

             call  GetRecpath_RandWalkCtrl(CtrlParam, swapPath) 
             call  Transfer_RandWalkCtrl0(SorFname, swapPath)
             return
  end subroutine Transfer_RandWalkCtrl1
  !****************************************************************************  

  !****************************************************************************
  subroutine Transfer_RandWalkCtrl2(CtrlParam)
  !***  PURPOSE:   to copy the original control data to path for output
  !     INPUT :    CtrlParam
  !
  !     OUTPUT
  !---
      implicit none
      class(RandWalkCtrl),  intent(in) ::CtrlParam
      !--- local variables
      character*256::swapPath

             call  GetRecpath_RandWalkCtrl(CtrlParam, swapPath) 
             call  Transfer_RandWalkCtrl0(CtrlParam%MyFname, swapPath)
             return
  end subroutine Transfer_RandWalkCtrl2
  !****************************************************************************  

 end module  MC_TypeDef_RandWalkCtrl
 !********************************************************************************

