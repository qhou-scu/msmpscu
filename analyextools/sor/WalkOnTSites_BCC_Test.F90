 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to simulate the MSD distribution for an particle walk on T-sites of
 !                  BCC lattice.
 !                  Ref.  CreateTetraLattice_BCC.F90
 !                  
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                  WalkOnTSites_BCC.exe 
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2019-06 (Hou Qing, Sichuan university)
 !
 !

 module WalkonTSites_BCC
 use MSM_CONSTANTS
 implicit none

      !--  the axsis in crystal
      !--- the primary 6 T-sites in first region cloest to the bcc atom at the origin
      !    we use four  
      integer, parameter, public ::mp_NPT=24 
      integer, parameter, private::mp_NNPT=mp_NPT*4
      
      integer,            private:: m_AXSISX(3) = (/1, 0, 0/)
      integer,            private:: m_AXSISY(3) = (/0, 1, 0/)
      integer,            private:: m_AXSISZ(3) = (/0, 0, 1/)

      real(KINDDF),       private::m_TSite(mp_NPT,3)
      real(KINDDF),       private::m_JumpVec(mp_NNPT,3)
      integer,            private::m_JumpSiteID(mp_NNPT)

      private::  Create_WalkSite_0, &
                 Create_WalkSite_1
      interface  Create_WalkSite_Tetra
          module procedure Create_WalkSite_0       
          module procedure Create_WalkSite_1  
      end interface  Create_WalkSite_Tetra
      
      
  contains

  !**********************************************************************************
   subroutine Create_WalkSite_0(AXX, AXY, AXZ, TSite, JumpVec, JumpSiteID)
   !***  PURPOSE:  to create T-sites of BCC and the possible jumping vectors for 
   !               for the next jump    
   !
   implicit none
  !----   DUMMY Variables
   integer::AXX(3), AXY(3), AXZ(3) 
   real(KINDDF), dimension(:,:)::TSite
   real(KINDDF), dimension(:,:)::JumpVec
   integer,      dimension(:)  ::JumpSiteID

  !----   Local variables
   real(KINDDF),parameter::TSite0(mp_NPT*3)=(/&
   0.50D0,  0.25D0,  0.00D0,    0.50D0,  0.00D0,  0.25D0,    0.50D0, -0.25D0,  0.00D0,   0.50D0,  0.00D0, -0.25D0,  &
   0.25D0,  0.50D0,  0.00D0,    0.00D0,  0.50D0,  0.25D0,   -0.25D0,  0.50D0,  0.00D0,   0.00D0,  0.50D0, -0.25D0,  &
   0.00D0,  0.25D0,  0.50D0,    0.25D0,  0.00D0,  0.50D0,    0.00D0, -0.25D0,  0.50D0,  -0.25D0,  0.00D0,  0.50D0,  &
  -0.50D0,  0.25D0,  0.00D0,   -0.50D0,  0.00D0,  0.25D0,   -0.50D0, -0.25D0,  0.00D0,  -0.50D0,  0.00D0, -0.25D0,  &
   0.25D0, -0.50D0,  0.00D0,    0.00D0, -0.50D0,  0.25D0,   -0.25D0, -0.50D0,  0.00D0,   0.00D0, -0.50D0, -0.25D0,  &
   0.00D0,  0.25D0, -0.50D0,    0.25D0,  0.00D0, -0.50D0,    0.00D0, -0.25D0, -0.50D0,  -0.25D0,  0.00D0, -0.50D0   &
   /)

   integer::I, J, IM, NN, I0, I1, I2, IP 
   real(KINDDF)::RMAT(3,3), SITEPOS(3), NEWSITEPOS(3), JUMPVEC0(mp_NNPT,3), P
   
          !$$--- to create the crystal coordinates
           RMAT(1,1:3) = AXX(1:3)/dsqrt(sum( dble(AXX(1:3)) * dble(AXX(1:3))))
           RMAT(2,1:3) = AXY(1:3)/dsqrt(sum( dble(AXY(1:3)) * dble(AXY(1:3))))
           RMAT(3,1:3) = AXZ(1:3)/dsqrt(sum( dble(AXZ(1:3)) * dble(AXZ(1:3))))

          !$$--- 
           NN = mp_NNPT/mp_NPT
           do I=1, mp_NPT
              !$$--- to determine the jump vector
              SITEPOS(1:3) = TSITE0((I-1)*3+1:I*3)
               do J=1, 3
                  IM = SITEPOS(J)*4.001
                  if(iabs(IM) .eq. 0) I0 = J
                  if(iabs(IM) .eq. 1) then
                     I1 = J
                     P  = IM
                  end if   
                  if(iabs(IM) .eq. 2) I2 = J
               end do
               IP = (I-1)*4
               JUMPVEC0(IP+1, I0) = 0.00D0
               JUMPVEC0(IP+1, I1) = 0.25D0*P
               JUMPVEC0(IP+1, I2) = 0.25D0

               JUMPVEC0(IP+2, I0) = 0.00D0
               JUMPVEC0(IP+2, I1) = 0.25D0*P
               JUMPVEC0(IP+2, I2) =-0.25D0

               JUMPVEC0(IP+3, I0) = 0.25D0
               JUMPVEC0(IP+3, I1) =-0.25D0*P
               JUMPVEC0(IP+3, I2) = 0.00D0               

               JUMPVEC0(IP+4, I0) =-0.25D0
               JUMPVEC0(IP+4, I1) =-0.25D0*P
               JUMPVEC0(IP+4, I2) = 0.00D0

              !$$--- to get the id of the sites where jumping to
               do IM=1,4
                  do J=1, 3
                     NEWSITEPOS(J) = SITEPOS(J) + JUMPVEC0(IP+IM, J)
                     if(dabs(NEWSITEPOS(J)) .gt. 0.5D0) then 
                        NEWSITEPOS(J) = dsign(0.5D0, NEWSITEPOS(J)) - NEWSITEPOS(J)
                     end if
                  end do      

                  do J=1, mp_NPT
                     P= sum( (NEWSITEPOS(1:3) - TSITE0((J-1)*3+1:J*3))**2.D0 )
                     if(P .lt. 1.D-6) then
                        JumpSiteID(IP+IM) = J 
                     end if  
                  end do 
               end do   
            end do                     
           
           !$$--- rotate the coordinate system
            do I=1, mp_NPT
               IP = (I-1)*4
               do J=1, 3
                  TSite(I,J) = TSITE0((I-1)*3+1)*RMAT(J,1) + TSITE0((I-1)*3+2)*RMAT(J,2) + TSITE0((I-1)*3+3)*RMAT(J,3)
                  JumpVec(IP+1,J) = JUMPVEC0(IP+1,1)*RMAT(J,1) + JUMPVEC0(IP+1,2)*RMAT(J,2) + JUMPVEC0(IP+1,3)*RMAT(J,3)
                  JumpVec(IP+2,J) = JUMPVEC0(IP+2,1)*RMAT(J,1) + JUMPVEC0(IP+2,2)*RMAT(J,2) + JUMPVEC0(IP+2,3)*RMAT(J,3)
                  JumpVec(IP+3,J) = JUMPVEC0(IP+3,1)*RMAT(J,1) + JUMPVEC0(IP+3,2)*RMAT(J,2) + JUMPVEC0(IP+3,3)*RMAT(J,3)
                  JumpVec(IP+4,J) = JUMPVEC0(IP+4,1)*RMAT(J,1) + JUMPVEC0(IP+4,2)*RMAT(J,2) + JUMPVEC0(IP+4,3)*RMAT(J,3)
               end do
            end do   
         return
   end subroutine Create_WalkSite_0
  !**********************************************************************************

  !**********************************************************************************
   subroutine Create_WalkSite_1()
   !***  PURPOSE:  to create T-sites of BCC and the possible jumping vectors for 
   !               for the next jump    
   !
      implicit none
         call Create_WalkSite_0(m_AXSISX, m_AXSISY, m_AXSISZ, m_TSite, m_JumpVec, m_JumpSiteID)   
       return
   end subroutine Create_WalkSite_1
  !**********************************************************************************

  !**********************************************************************************
   subroutine Walk_For_OneStep_0(CurSiteID, Pos, JumptoSiteID)
   !***  PURPOSE:  to get moving vector for next jump
   !
      use RAND32_MODULE
      implicit none
   !----   DUMMY Variables      
          integer,       intent(in) ::CurSiteID
          real(KINDDF),  intent(out)::Pos(3)
          integer,       intent(out)::JumptoSiteID
    !---  Local variable
          integer::IC, I 
          
              !$$--- randomly select next site for 4 candiadtes
               I            = 4.D0*DRAND32() + 1
               IC           = (CurSiteID-1)*4 + I
               Pos(1:3)     = Pos(1:3) + m_JumpVec(IC,1:3) 
               JumptoSiteID = m_JumpSiteID(IC) 
       return
   end subroutine Walk_For_OneStep_0   
  !**********************************************************************************   

  !**********************************************************************************
   subroutine Walk_For_OneStep_1(CurSiteID, Pos, JumptoSiteID)
   !***  PURPOSE:  to get moving vector for next jump
   !
      implicit none
   !----   DUMMY Variables      
          integer,       dimension(:),   intent(in) ::CurSiteID
          real(KINDDF),  dimension(:,:), intent(out)::Pos
          integer,       dimension(:),   intent(out)::JumptoSiteID
    !---  Local variable
          integer::J, NP
          
               NP = size(CurSiteID)
               do J=1, NP 
                   call Walk_For_OneStep_0(CurSiteID(J), Pos(J,1:3), JumptoSiteID(J))
               end do   
       return
   end subroutine Walk_For_OneStep_1   
  !**********************************************************************************   

  !**********************************************************************************
   subroutine Walk_For_OneStep_2(CurSiteID, Pos)
   !***  PURPOSE:  to get moving vector for next jump
   !
      implicit none
   !----   DUMMY Variables      
          integer,       dimension(:),   intent(inout)::CurSiteID
          real(KINDDF),  dimension(:,:), intent(inout)::Pos
    !---  Local variable
               call Walk_For_OneStep_1(CurSiteID, Pos, CurSiteID)
       return
   end subroutine Walk_For_OneStep_2   
  !**********************************************************************************   

  !**********************************************************************************
   subroutine ResetPosition(SiteID, Pos)
   !***  PURPOSE:  to move the partilce to new position for one step
   !
      implicit none
   !----   DUMMY Variables      
          integer,       dimension(:),   intent(in)   ::SiteID
          real(KINDDF),  dimension(:,:), intent(inout)::Pos
    !---  Local variable
          integer::J 
          
               do J=1, size(SiteID) 
                  Pos(J,1:3) = Pos(J,1:3) + m_TSite(SiteID(J),1:3)
               end do   
       return
   end subroutine ResetPosition   
  !**********************************************************************************     

 end module WalkonTSites_BCC

 !**********************************************************************************
 program WalkonTSites_BCC_main
    use MD_TYPEDEF_SimMDBox
    use WalkonTSites_BCC
    use RAND32_MODULE
    implicit none
    integer, parameter::NP = 200000 
    integer, parameter::NR = 10 
    integer, parameter::NB = 100 

    type(SimMDBox)::SimBox0, SimBox
    type(MDRecordStamp):: Stamp
    real(KINDDF)::VEC(3), DIS 
    integer::CURID(NP)
    integer::I , J, IB, IP, INB, hFile

    real(KINDDF),parameter::BINS = 0.8*0.5D0/2.D0**0.5
    integer,     parameter::BINN = 2000
    integer, dimension(BINN)::His
    real(KINDDF)::MSD(NB)
    character*256:: fname="WalkonTSites_BCC_Test", OutFname="" 

       !--- initialize jumping sites
       call Create_WalkSite_Tetra()

       !--- 
       SimBox0%NPRT = NP
       call Initialize_SimMDBox(SimBox0, 2)
       SimBox0%XP = 0.D0
       SimBox0%DIS = 0.D0
       do I=1, NP
          CURID(I) = mp_NPT*DRAND32() + 1
       end do
       
       call ResetPosition(CURID, SimBox0%XP)
       SimBox = SimBox0
       MSD = 0.D0 
       Stamp%ITest = 1
       Stamp%IBox  = 1
       Stamp%ICfg  = 0
       !call Putout_Instance_Config_SimMDBox("WalkonTSites_BCC_Test", SimBox, Stamp)              
       do IB=1, NB 
          do J=1, NR
             call Walk_For_OneStep_2(CURID, SimBox%XP)
          end do
          Stamp%ITest = 1
          Stamp%IBox  = 1
          Stamp%ICfg  = IB
          !call Putout_Instance_Config_SimMDBox(fname, SimBox, Stamp)       

          !--- calculate the distance distribution
          HIS = 0
          do IP = 1, NP
             VEC(1:3) = SimBox%XP(IP, 1:3) - SimBox0%XP(IP, 1:3) 
             DIS      = dsqrt(sum(VEC*VEC))
             INB      = int(DIS/BINS)+ 1
             His(INB) = His(INB) + 1     
             MSD(IB)  = MSD(IB) + DIS*DIS
          end do
          call STRCATI(OutFname, Fname, "_DIS.", IB, 4)
          call AvailableIOUnit(hFile)
          open(UNIT=hFile, FILE=OutFname)
          write(hFile, *) "! DISHis from "//Fname(1:len_trim(Fname))
          write(hFile, *) "! For CFG=", IB, NR
          do I = 1, size(His)
            write(hFile, fmt="(I8,1x,1PE14.4,1x,I8,1x,1PE14.4)")  I, (I-1)*BINS, His(I), dble(His(I))/CP_FOURPI/(I*BINS)**2
          end do 
          close(hFile)
      end do   

      OutFname = Fname(1:len_trim(Fname))//"_MSD"
      call AvailableIOUnit(hFile)
      open(UNIT=hFile, FILE=OutFname)
      write(hFile, *) "! DISHis from "//Fname(1:len_trim(Fname))
      write(hFile, *) "! For CFG=", IB, NR
      do IB = 1, NB
         write(hFile, fmt="(I8,1x,1PE14.4,1x,I8,1x,1PE14.4)")  IB, MSD(IB)/dble(NP)
      end do 
      close(hFile)

     stop         
 end program WalkonTSites_BCC_main

