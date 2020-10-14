module NIST_ForceTable
!*** DESCRIPTION: This module is to generator potential tables for datas from NIST.
!                 It was copied from MDPSCU.
!*** AUTHOR:
!    Mingjie Qiu, Oct., 2016
!    Modified 2019-10-16, Hou Qing
!                      allowing multiple input of libs for using Joshon methos to 
!                      generate of table of alloy    
!----------------------------------------------------------------------------------
use MD_CONSTANTS
use MD_TYPEDEF_SimMDBox
use MD_TYPEDEF_SimMDCtrl
use MD_TYPEDEF_ForceTable
implicit none
   character*32, private:: m_Exttype

   private:: Generate_NIST_ForceTalbe
   private:: NIST_Register_Interaction_Table0,  &
             NIST_Register_Interaction_Table1,  &
             NIST_Register_Interaction_Table2
   public::  NIST_Register_Interaction_Table
   interface NIST_Register_Interaction_Table
       module procedure NIST_Register_Interaction_Table0          
         module procedure NIST_Register_Interaction_Table1
            module procedure NIST_Register_Interaction_Table2 
    end interface            
contains
!****************************************************************************************
!****************************************************************************************
subroutine NIST_Register_Interaction_Table2(SimBox, CtrlParam, FTable, printout)
  !***  PURPOSE:   to initialize the force table to be used in force calculations
  !
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !                FTable, the force tabled to be registered
  !
  !     OUTPUT:    FTable, the registered table
  !----------------------------------------------------------------------------------------
  !use Filedatas_Func_Lspt, Release_Lspt=>Release, NN_Spline_Lspt=>NN_Spline,&
  !                         RHO_Spline_Lspt=>RHO_Spline, EMBED_Spline_Lspt=>EMBED_Spline,&
  !                         Notify_func_Lspt=>Notify_func
  use Filedatas_Func_Lspt, only: Release_Lspt=>Release, Register_ForceTableProc_SPT
  use Filedatas_Func_Setfl,only: Release_Setfl=>Release, Register_ForceTableProc_Setfl
                           !NN_Spline_Setfl=>NN_Spline,&
                           !RHO_Spline_Setfl=>RHO_Spline, EMBED_Spline_Setfl=>EMBED_Spline,&
                           !Notify_func_Setfl=>Notify_func
  use Filedatas_Func_Moldy, only:Release_Moldy=>Release, Register_ForceTableProc_Moldy
  use My_MiniUtilities_QMJ
  implicit none
      !--- dummy vaiables
      type(SimMDBox),     intent(in)   ::SimBox
      type(SimMDCtrl),    intent(in)   ::CtrlParam
      type(MDForceTable), intent(inout)::FTable
      integer,             optional    ::printout
      !--- local vaiables
      character*128, dimension(:), allocatable::POTNAME
      character*128::POTNAME0
      integer     :: I, NLIB, EX

           !*** to find out how many potlibs
           NLIB = 0
           POTNAME0 = adjustl(SimBox%PotLibname)
           do while(len_trim(POTNAME0).gt.0)
              do I=1, len_trim(POTNAME0)
                 if(POTNAME0(I:I) .eq. " " .or. POTNAME0(I:I) .eq. ","  .or. &
                    I.eq.len_trim(POTNAME0)) then
                    NLIB = NLIB + 1
                    POTNAME0 = adjustl(POTNAME0(I+1:))
                    exit
                 end if
               end do 
           end do  
           
           allocate(POTNAME(NLIB))
           POTNAME  = ""
           NLIB = 0
           POTNAME0 = adjustl(SimBox%PotLibname)
           do while(len_trim(POTNAME0).gt.0)
              do I=1, len_trim(POTNAME0)
               if(POTNAME0(I:I) .eq. " " .or. POTNAME0(I:I) .eq. "," .or. I.eq.len_trim(POTNAME0)) then
                  NLIB = NLIB + 1
                  if(POTNAME0(I:I) .eq. " " .or. POTNAME0(I:I) .eq. ",") then
                     POTNAME(NLIB) = POTNAME0(1:I-1)
                  else   
                     POTNAME(NLIB) = POTNAME0(1:I)
                  end if   
                  POTNAME0 = adjustl(POTNAME0(I+1:))
                  exit
               end if
              end do 
           end do  

           !--- to check if the lib exist
           do I=1, NLIB
              inquire(FILE=POTNAME(I)(1:len_trim(POTNAME(I))), exist=EX)
              if(.not.EX) then
                 write(*,fmt="(A)") ' MDPSCU Error: potential lib "'// &
                                     POTNAME(I)(1:len_trim(POTNAME(I)))//'" is not a potential libname,'// &
                                    'check input box file in  &POTSUBCTL subsection.'
                 write(*,fmt="(A)") ' For external potential lib, check also: the existence of files '//&
                                     POTNAME(I)(1:len_trim(POTNAME(I)))//'.pair'//                             &
                                   ' and '//POTNAME(I)(1:len_trim(POTNAME(I)))//'.embd'

                 write(*,fmt="(A)") ' or: the extistence of: '//POTNAME(I)(1:len_trim(POTNAME(I)))//'.lspt'
                 write(*,fmt="(A)") ' or: '//POTNAME(I)(1:len_trim(POTNAME(I)))//'.moldy'
                 write(*,fmt="(A)") ' or: '//POTNAME(I)(1:len_trim(POTNAME(I)))//'.setfl'
                 write(*,fmt="(A)") ' Process to be stopped'
                 stop
              end if
           end do    
           if(NLIB .eq. 1) then
              call NIST_Register_Interaction_Table0(Potname(1), FTable)   
           else   
              call NIST_Register_Interaction_Table1(NLib, Potname, FTable)
           end if    
           call Export_ForceTable(SimBox%potlibname(1:len_trim(SimBox%potlibname)),FTable)

           deallocate(POTNAME)
          return
end subroutine NIST_Register_Interaction_Table2
!*****************************************************************************************

!****************************************************************************************
subroutine NIST_Register_Interaction_Table0(PotLibname, FTable)
   !***  PURPOSE:   to initialize the force table to be used in force calculations
   !                
   !
   !     INPUT:     PotLibname: the nmae of the lib
   !                FTable, the force tabled to be registered
   !
   !     OUTPUT:    FTable, the registered table
   !----------------------------------------------------------------------------------------
   !use Filedatas_Func_Lspt, Release_Lspt=>Release, NN_Spline_Lspt=>NN_Spline,&
   !                         RHO_Spline_Lspt=>RHO_Spline, EMBED_Spline_Lspt=>EMBED_Spline,&
   !                         Notify_func_Lspt=>Notify_func
   use Filedatas_Func_Lspt, only: Release_Lspt=>Release, Register_ForceTableProc_SPT
   use Filedatas_Func_Setfl,only: Release_Setfl=>Release, Register_ForceTableProc_Setfl
                            !NN_Spline_Setfl=>NN_Spline,&
                            !RHO_Spline_Setfl=>RHO_Spline, EMBED_Spline_Setfl=>EMBED_Spline,&
                            !Notify_func_Setfl=>Notify_func
   use Filedatas_Func_Moldy, only:Release_Moldy=>Release, Register_ForceTableProc_Moldy
   use My_MiniUtilities_QMJ
   implicit none
       !--- dummy vaiables
       character*(*),       intent(in)::PotLibname
       type(MDForceTable),  intent(inout)::FTable

       !--- local vaiables
       integer     :: I, NE
       !*** register the table according to user supplied routines
            call Get_FileExtent(PotLibname, m_Exttype)
            select case(m_Exttype)
                   case ("moldy")
                          call Register_ForceTableProc_Moldy(PotLibname, FTable)
                   case ("lspt")
                         call Register_ForceTableProc_SPT(PotLibname, FTable)
                   case ("setfl")
                         call Register_ForceTableProc_Setfl(PotLibname, FTable)
                   case default
                         write(*,fmt="(A)") ' MDPSCU Error: wrong potential libname "'// &
                                     Potlibname(1:len_trim(Potlibname))//'"'
                         write(*,fmt="(A)") '                NIST potential libname should have extension: "lspt", or "moldy", or "setfl" '
                         write(*,fmt="(A)") ' Process to be stopped'
                         stop
            end select
           !****  create the force table
           call  Generate_NIST_ForceTalbe (FTable)

           !*** release the space
            select case(m_Exttype)
                   case ("lspt")
                        call Release_Lspt()
                   case ("setfl")
                        call Release_Setfl()
                   case ("moldy")
                        call Release_Moldy()
            end select
         return
 end subroutine NIST_Register_Interaction_Table0
 !*****************************************************************************************
 
 !****************************************************************************************
 subroutine NIST_Register_Interaction_Table1(NLib, PotLibname, FTable)
   !***  PURPOSE:   to initialize the force table to be used in force calculations
   !                if multiple libs are input, Joshon method was used to the pairwise
   !                interaction between different elements
   !
   !     INPUT:     NLib:       the number of libs
   !                PotLibname: the names for libs
   !                FTable, the force tabled to be registered
   !
   !     OUTPUT:    FTable, the registered table
   !----------------------------------------------------------------------------------------
   !use Filedatas_Func_Lspt, Release_Lspt=>Release, NN_Spline_Lspt=>NN_Spline,&
   !                         RHO_Spline_Lspt=>RHO_Spline, EMBED_Spline_Lspt=>EMBED_Spline,&
   !                         Notify_func_Lspt=>Notify_func
   use Filedatas_Func_Lspt, only: Release_Lspt=>Release, Register_ForceTableProc_SPT
   use Filedatas_Func_Setfl,only: Release_Setfl=>Release, Register_ForceTableProc_Setfl
                            !NN_Spline_Setfl=>NN_Spline,&
                            !RHO_Spline_Setfl=>RHO_Spline, EMBED_Spline_Setfl=>EMBED_Spline,&
                            !Notify_func_Setfl=>Notify_func
   use Filedatas_Func_Moldy, only:Release_Moldy=>Release, Register_ForceTableProc_Moldy
   implicit none
       !--- dummy vaiables
       integer,                     intent(in)::NLib 
       character*(*), dimension(:), intent(in)::PotLibname
       type(MDForceTable),          intent(inout)::FTable
       !--- local vaiables
       integer     :: I, J, I1, J1, II, ITAB, ITAB1, ITAB2, NUM1, NUM2, NUM, TAB1, TAB2, NR, NE, IERR
       real(KINDDF):: A, B, DA, DB, WASB, WBSA, DWASB, DWBSA
       type(MDForceTable),dimension(:),allocatable::tFTable
       character*16::tstr
       character*256::title

            allocate(tFTable(NLib))
           
            !--- create the table from each Lib
            do I=1, NLib
               tFTable(I)%Rmax   = FTable%Rmax  
               tFTable(I)%NTab   = FTable%NTab
               tFTable(I)%PotType = FTable%PotType
               call NIST_Register_Interaction_Table0(PotLibname(I), tFTable(I))
            end do
            !--- combine the libs
            call Release_ForceTable(FTable)
            !--- the number of grids and cutoff distance
            FTable%Rmax     = tFTable(1)%Rmax
            FTable%RmaxSQRT = tFTable(1)%RmaxSQRT
            FTable%CSI      = tFTable(1)%CSI
            FTable%CSIV     = tFTable(1)%CSIV
            FTable%NEMBD    = tFTable(1)%NEMBD
            FTable%NTab     = tFTable(1)%NTab

            !--- the pairwise potential
            NR  = FTable%NTAB
            NUM = 0
            do I=1, NLib
               NUM2 = count(tFTable(I)%FPAIR.gt.0)
               NUM1 = dsqrt(dble(NUM2))+0.0001
               NUM  = NUM + NUM1
            end do
            NUM1= NUM
            NUM = NUM*NUM
            allocate(FTable%FPAIR(NUM),   FTable%HASREG(NUM), STAT=IERR )
            allocate(FTable%POTR(NUM,NR), FTable%FPOTR(NUM,NR),        &
                     FTable%POTB(NUM,NR), FTable%FPOTB(NUM,NR),        &
                     STAT=IERR )
            FTable%POTR   = 0.D0
            FTable%FPOTR  = 0.D0
            FTable%POTB   = 0.D0
            FTable%FPOTB  = 0.D0
            FTable%FPAIR  = 0
            FTable%HASREG = 0
            ITAB = 0
             do ITAB1 =1, NLib
                NUM1 = dsqrt(count(tFTable(ITAB1)%FPAIR.gt.0) + 0.0001)
                do I=1, NUM1
                   do ITAB2 =1, NLib
                      NUM2 = dsqrt(count(tFTable(ITAB2)%FPAIR.gt.0) + 0.0001)
                      do J=1, NUM2
                         ITAB = ITAB + 1
                         if(ITAB1 .eq. ITAB2) then
                            J1 = (I-1)*NUM1 + J 
                            do II=1, FTable%NTab
                               FTable%POTR(ITAB, II)  = tFTable(ITAB2)%POTR(J1,II)
                               FTable%FPOTR(ITAB,II)  = tFTable(ITAB2)%FPOTR(J1,II)
                               FTable%POTB(ITAB, II)  = tFTable(ITAB2)%POTB(J1,II)
                               FTable%FPOTB(ITAB,II)  = tFTable(ITAB2)%FPOTB(J1,II)
                            end do  
                            !--- we need also add note to the TBALE
                            !call GetIthExtProc_ForceExtProcList(tFTable(ITAB2)%REGPROCS, J1, pExtProc)
                            !call Register_ForceTableProc(ITAB, FTable, NULL_NNFORCE, NOTE=pExtProc%POTDescript) 

                         else  !--- using the Jonhshon interpolation
                            I1 = (I-1)*NUM1 + I
                            J1 = (J-1)*NUM2 + J
                            do II=1, NR
                               A     = tFTable(ITAB1)%POTB(I1,II)
                               B     = tFTable(ITAB2)%POTB(J1,II)
                               DA    = tFTable(ITAB1)%FPOTB(I1,II)
                               DB    = tFTable(ITAB2)%FPOTB(J1,II)
                               if(A .gt. 0.D0 .and. B.gt.0.D0) then
                                  WASB  = A/B
                                  WBSA  = B/A
                                  DWASB = (DA*B-A*DB)/B/B
                                  DWBSA = (DB*A-B*DA)/A/A
                                  A     = tFTable(ITAB1)%POTR(I1,II)
                                  B     = tFTable(ITAB2)%POTR(J1,II)
                                  DA    = tFTable(ITAB1)%FPOTR(I1,II)
                                  DB    = tFTable(ITAB2)%FPOTR(J1,II)
                                  FTable%POTR(ITAB, II) = C_HALF*(WASB*B+WBSA*A)
                                  FTable%FPOTR(ITAB,II) = C_HALF*(WASB*DB+WBSA*DA)+(DWASB*B+DWBSA*A)
                               else   
                                  FTable%POTR(ITAB, II) = 0.D0
                                  FTable%FPOTR(ITAB,II) = 0.D0
                               end if
                               FTable%POTB(ITAB, II)  = tFTable(ITAB2)%POTB(J1,II)
                               FTable%FPOTB(ITAB,II)  = tFTable(ITAB2)%FPOTB(J1,II)
                          end do  
                         end if 
                         FTable%FPAIR(ITAB) = ITAB  

                         !--- we need also add note to the TBALE
                         title =" Jonshon interpelation: Lib "
                         write(tstr,*) ITAB1; tstr = adjustl(tstr)
                         title = title(1:len_trim(title))//tstr(1:len_trim(tstr))//" atom"
                         write(tstr,*) I; tstr = adjustl(tstr)
                         title = title(1:len_trim(title))//tstr(1:len_trim(tstr))//" <- LIB "
                         write(tstr,*) ITAB2; tstr = adjustl(tstr)
                         title = title(1:len_trim(title))//tstr(1:len_trim(tstr))//" atom "
                         write(tstr,*) J; tstr = adjustl(tstr)
                         title = title(1:len_trim(title))//adjustl(tstr)
                         call Register_ForceTableProc(ITAB, FTable, NULL_NNFORCE, NOTE=title) 

                      end do  !--- end loop for J
                   end do  !--- end loop for ITAB2
                end do   
             end do
             
            !--- the embedment function
             NR = FTable%NEMBD
             allocate(FTable%FEMBD(NUM,NR), FTable%DFEMBD(NUM,NR),FTable%FPAIR1(NUM), FTable%HASREG1(NUM), STAT=IERR)
             FTable%FEMBD  = 0.D0
             FTable%DFEMBD = 0.D0
             FTable%HASREG1 = 0
             FTable%FPAIR1  = 0
             ITAB = 0
             do ITAB1 =1, NLib
                NUM1 = dsqrt(count(tFTable(ITAB1)%FPAIR.gt.0) + 0.0001)
                do I=1, NUM1
                   do ITAB2 =1, NLib
                      NUM2 = dsqrt(count(tFTable(ITAB2)%FPAIR.gt.0) + 0.0001)
                      do J=1, NUM2
                         ITAB = ITAB + 1
                         if(ITAB1 .eq. ITAB2 .and. I.eq.J) then
                            do II=1, NR
                               FTable%FEMBD(ITAB, II) = tFTable(ITAB1)%FEMBD(I,II)
                               FTable%DFEMBD(ITAB,II) = tFTable(ITAB1)%DFEMBD(I,II)
                            end do   
                         end if 
                         FTable%FPAIR1(ITAB) = ITAB 
                      end do
                   end do      
                end do
             end do  

             do I=1, NLib
               call Release_ForceTable(tFTable(I))
             end do   
             deallocate(tFTable)
          return
 end subroutine NIST_Register_Interaction_Table1
 !*****************************************************************************************

!*****************************************************************************************
subroutine Generate_NIST_ForceTalbe (FTable)
!***  PURPOSE:   to generate a force table that could be used to by other MD packaged.
  !                The difference between this routine and Create_Interaction_ForceTable
  !                is that all the registered force table will be generated, and output.
  !
  !
  !                 NOTE: This routine dose not need SimBox, the range of R and
  !                       the number of grid pointer of the table should be given before
  !                       calling this routine.
  !
  !      INPUT:     FNAME: the name of the library
  !                 Fable: the table with its forces have been registered
  !
  !      OUTPUT:    FTable, the force table store in a file
  !
  !***************************************************************************************
   implicit none
      !--- dummy vaiables
      type(MDForceTable),intent(inout)::FTable
      !local variables
      integer::IFORCE, NUM, NR, NE, IERR
      real(KINDDF) :: RHOMAX, RHOD

      RHOMAX = FTable%RHOMX                   !*** added by QMJ
      RHOD   = FTable%RHOMX/dble(FTable%NEMBD)
            call Release_ForceTable(FTable, keepreglist=1)
            if(FTable%Rmax .le. 0.D0) then
              write(*,fmt="(A)") "MDPSCU Error: the distance range is missed in generating force table"
              write(*,fmt="(A)") "              Process to be stopped"
              stop
            end if
            if(FTable%NTAB .le. 0.D0) then
              write(*,fmt="(A)") "MDPSCU Error: the number of grid is less than zero in generating force table"
              write(*,fmt="(A)") "              Process to be stopped"
              stop
            end if
            FTable%RmaxSQRT = DSQRT(FTable%Rmax)
            FTable%CSI   = dble(FTable%NTAB)/FTable%RmaxSQRT
            FTable%CSIV  = 1.D0/FTable%CSI
            call GetNumExtProc_ForceExtProcList(FTable%REGPROCS, NUM)
            NR = FTable%NTAB
            allocate(FTable%FPAIR(NUM), FTable%HASREG(NUM), STAT=IERR )
            allocate(FTable%POTR(NUM,NR), FTable%FPOTR(NUM,NR),        &
                     FTable%POTB(NUM,NR), FTable%FPOTB(NUM,NR),        &
                     STAT=IERR )
            FTable%POTR   = 0.D0
            FTable%FPOTR  = 0.D0
            FTable%POTB   = 0.D0
            FTable%FPOTB  = 0.D0
            FTable%FPAIR  = 0
            FTable%HASREG = 0
            do IFORCE=1, NUM
               FTable%FPAIR(IFORCE) = IFORCE
               call Create_Pairwise_ForceTable(IFORCE, FTable)
            end do
            if (m_Exttype .eq. "lspt".or. m_Exttype .eq. "setfl") then  !***When the file from NIST is "***.lspt"or"***.setfl", we use
              FTable%RHOMX = RHOMAX                                     !*** the RHOMAX and RMAX from NIST file.
              FTable%RHOD  = RHOD
            end if
            NE = FTable%NEMBD
            allocate(FTable%FEMBD(NUM,NE), FTable%DFEMBD(NUM,NE),FTable%FPAIR1(NUM), FTable%HASREG1(NUM), STAT=IERR)
            FTable%FEMBD  = 0.D0
            FTable%DFEMBD = 0.D0
            FTable%HASREG1 = 0
            FTable%FPAIR1  = 0
            do IFORCE=1, NUM
               FTable%FPAIR1(IFORCE) = IFORCE
               call Create_EMBDFUNTable(IFORCE, FTable)
            end do
           return
end subroutine Generate_NIST_ForceTalbe
!**************************************************************************************
end module NIST_ForceTable
