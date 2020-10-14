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
 module SIAEventSearch_module
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none
         integer, private::m_NREC                                         ! number of actually record
                                                                          ! atom out-of-box are filtered

         character(len=256),  private::m_OUTPATH =""                      ! the path output data
         integer,             private::m_j_v_IOUnit  = 0                  ! temp I/O unit to storing jumping vectors

         integer, parameter,                         private::mp_MXI = 1  ! the permitted number of interstital in a box
         integer, parameter,                         private::mp_MXAI =4  ! the permitted number of atoms occupy one site
         integer, dimension(:),         allocatable, private::m_preSITEi  ! previous sites of interstitials
         integer, dimension(:),         allocatable, private::m_preOCCUai ! previous atom id that occupy an interstitial site
         real(KINDDF), dimension(:,:),  allocatable, private::m_prePOSi   ! previous position of interstitial sites
         real(KINDDF), dimension(:,:),  allocatable, private::m_preOCCUaP ! previous atom positions that occupy an interstitial site
         real(KINDDF), dimension(:),    allocatable, private::m_preTIMEi  ! previous of recording

         integer, dimension(:),allocatable,          private::m_hSITE     ! working space, storing the sites of atoms occupying
         integer, dimension(:),allocatable,          private::m_hSTAT     ! working space, state of sites, state >= 2 indicates an interstitial site
         integer, dimension(:),allocatable,          private::m_hOCCUa    ! working space, the ids of atoms occupying the sites
         integer, dimension(:),allocatable,          private::m_hSITEi    ! current sites of interstitials
         integer, dimension(:),allocatable,          private::m_hOCCUai   ! current atoms id that occupy an interstitial site
         real(KINDDF), dimension(:,:),allocatable,   private::m_hSitePOS  ! current position of interstitial sites
         real(KINDDF), dimension(:,:),  allocatable, private::m_hOCCUaP   ! current atom positions that occupy an interstitial site

         type(SimMDBox), dimension(:),  allocatable, private:: m_SimBoxSwap
  contains
  !**********************************************************************************
   subroutine Initialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use RefVoronoiVacancy_GPU, only:Initialize_ReferenceVacancy_DEV, Get_filename_Output
   use MD_TYPEDEF_PrintList,  only:Add_ArchiveProcess    
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
  integer::ERR  

         if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") then
           call Clear(SimBox, CtrlParam)
           return
         end if

         !$$--- loading control parameters
         call Initialize_ReferenceVacancy_DEV(SimBox, CtrlParam)
         call Get_filename_Output(m_OUTPATH)
         call GetPath(m_OUTPATH, m_OUTPATH)
         call AvailableIOUnit(m_j_v_IOUnit)
         !open(m_j_v_IOUnit, FILE="status='SCRATCH')
         if(CtrlParam%RESTART .gt. 0) then
            open(unit=m_j_v_IOUnit, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"SIA.Time_Displace",POSITION='APPEND')
         else
             open(unit=m_j_v_IOUnit, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"SIA.Time_Displace")
         end if
         m_NREC = 0
         !---
         call Add_ArchiveProcess("&SIAEventSearch", Archive, Restore)
         return
   end subroutine Initialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine DoEventCheck(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT,SynchronizeDevices, hm_ITYP, hm_XP, hm_GID
  use RefVoronoiVacancy_GPU, only:Cal_Occupied_Voronoi_DEV, hm_RefSimBox
  use MD_SimBoxArray
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        integer::I, J, IB, NBOX, IP, IA, NA, NPRT0, NPRT, SHIFT, NI, IS, OLDSITE(mp_MXI),       &
                 NEWSITE(mp_MXI), NEWOCCUai(mp_MXI*mp_MXAI), NEARSITE(mp_MXI), CHANGED(mp_MXI), &
                 NOLD, NNEW, JOB
        real(KINDDF)::OLDPOS(mp_MXI,3), NEWPOS(mp_MXI,3), PRETIME(mp_MXI),  &
                      DX(mp_MXI), DY(mp_MXI), DZ(mp_MXI), DR(mp_MXI), TIME, NEWOCCUaP(mp_MXI*mp_MXAI,3), &
                      DIST, R2, LATT, BS(3), HBS(3)
        character*128::FMT
        character*32::SNUM
        save FMT

  !-------------
           JOB   = Stamp%ITest
           TIME  = Stamp%Time

           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)
           NPRT  = dm_NPRT/NBOX
           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           HBS   = C_HALF*BS

         !--- calculate the occupation of lattice sites
           call Cal_Occupied_Voronoi_DEV(SimBox(1), CtrlParam, m_hSITE)

          !--- to calculate occupation state
           m_hSTAT   = 0
           m_hOCCUa  = 0
           m_hOCCUaP = 0.D0
           IP = 0
           do IB=1, NBOX
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT
                 IP = IP + 1
                 if(m_hSITE(IP)>0) then
                    IS = m_hSITE(IP)+SHIFT
                    m_hSTAT(IS) = m_hSTAT(IS) + 1
                    do J=1, mp_MXAI
                       !--- m_hOCCUa: the atoms ID on site SHIFT + m_hSITE(IP)
                       if(m_hOCCUa((IS-1)*mp_MXAI + J).eq. 0) then
                          m_hOCCUa((IS-1)*mp_MXAI + J) = IP
                          exit
                       end if
                    end do
                 end if
              end do
           end do

          !--- to extract the postions and numbers of interstitial in each box
           IP = 0
           m_hSITEi =  0
           m_hOCCUai = 0
           do IB=1, NBOX
              NI = 0
              NA = 0
              !--- to extract the postions and numbers of interstitial in each box
              do I=1, NPRT0
                 IP = IP + 1
                 if(m_hSTAT(IP) .ge. 2) then
                    !--- this site containing interstitial atoms
                    NI = NI + 1
                    if(NI .gt. mp_MXI) exit
                    IS = (IB-1)*mp_MXI+NI
                    m_hSITEi(IS) = I
                    m_hSitePOS(IS,1:3) = hm_RefSimBox%XP(I, 1:3)

                    do J=1, mp_MXAI
                       IA = m_hOCCUa((IP-1)*mp_MXAI+J)
                       if(IA .gt. 0) then
                          NA = NA + 1
                          m_hOCCUai((IS-1)*mp_MXAI+J)     = hm_GID(IA) - int((hm_GID(IA)-1)/NPRT)*NPRT
                          m_hOCCUaP((IS-1)*mp_MXAI+J,1:3) = hm_XP(IA, 1:3)
                       end if
                    end do
                 end if
              end do

              if(NI .gt. mp_MXI) then
                 write(*,fmt="(A, 1x, I4, 1x, A, I4, A, I4)") &
                      " MDPSCU Warning: the number of interstitials in box #", IB, " larger than permitted value:", NI, " vs ",mp_MXI
                 call ONWARNING(gm_OnWarning)
               end if
           end do

           if(m_NREC .eq. 0 ) then
              m_preSITEi  = m_hSITEi
              m_prePOSi   = m_hSitePOS
              m_preOCCUai = m_hOCCUai
              m_preOCCUaP = m_hOCCUaP
              m_preTIMEi  = TIME
              m_NREC = m_NREC + 1
              write(SNUM, *)mp_MXAI
              SNUM = adjustl(SNUM)
              SNUM = SNUM(1:len_trim(SNUM))//"(I8,1X, 3(1PE14.5,1x))"
              FMT = "(1x, I8, 1x, 1PE14.5, 1x, I8, 1PE14.5, 1x, I8, 1x, 9(1PE13.4,1x),"//SNUM(1:len_trim(SNUM))//","//SNUM(1:len_trim(SNUM))//")"
              return
           end if

          !--- to check if exchange event occur
           do IB=1, NBOX
              OLDSITE = 0
              NEWSITE = 0
              !$$--- to find out which site has been changed
              !$$    NOLD: the number of defect that have been changed 
              NOLD    = 0
              CHANGED = 0
              do I=1, mp_MXI
                 IS = (IB-1)*mp_MXI+I
                 if(.not.any(m_hSITEi((IB-1)*mp_MXI+1:IB*mp_MXI) .eq. m_preSITEi(IS)) ) then
                    NOLD = NOLD + 1
                    OLDSITE(NOLD)    = m_preSITEi(IS)
                    OLDPOS(NOLD,1:3) = m_prePOSi(IS, 1:3)
                    PRETIME(NOLD)    = m_preTIMEi(IS)
                    CHANGED(NOLD)    = IS
                 end if
              end do
              !--- to find out new sites not existing in old sites
              NNEW = 0
              do I=1, mp_MXI
                 IS = (IB-1)*mp_MXI+I
                 if(.not.any(m_preSITEi((IB-1)*mp_MXI+1:IB*mp_MXI) .eq. m_hSITEi(IS)) ) then
                    NNEW = NNEW + 1
                    NEWSITE(NNEW) = m_hSITEi(IS)
                    NEWPOS(NNEW,1:3) = m_hSitePOS(IS, 1:3)
                    NEWOCCUai((NNEW-1)*mp_MXAI+1:NNEW*mp_MXAI) = m_hOCCUai((IS-1)*mp_MXAI+1:IS*mp_MXAI)
                    NEWOCCUaP((NNEW-1)*mp_MXAI+1:NNEW*mp_MXAI,1:3) = m_hOCCUaP((IS-1)*mp_MXAI+1:IS*mp_MXAI,1:3)
                 end if
              end do

              if(NNEW .ne. NOLD) then
                 write(*,fmt="(A, ix, I4, 1x, I4, 1x, A, 1x, I4)")   " MDPSCU Warning: number of interstitials is changed in box  ", &
                                                                     (JOB-1)*NBOX+IB, NNEW, " vs ", NOLD
                 call ONWARNING(gm_OnWarning)
              end if
               
              !--- if no event 
              if(NOLD .gt. 0)  then
                 !--- to find out the neasest new site from old sites
                 NEARSITE = 0
                 DX = 0.D0
                 DY = 0.D0
                 DZ = 0.D0
                 DR = 0.D0
                 do I=1, NOLD
                    DIST = 1.D30
                    do J=1, NNEW
                       DX(I) = NEWPOS(J,1) - OLDPOS(I,1); if(dabs(DX(I)) .GT. HBS(1))  DX(I) = DX(I) - DSIGN(BS(1),DX(I))
                       DY(I) = NEWPOS(J,2) - OLDPOS(I,2); if(dabs(DY(I)) .GT. HBS(2))  DY(I) = DY(I) - DSIGN(BS(2),DY(I))
                       DZ(I) = NEWPOS(J,3) - OLDPOS(I,3); if(dabs(DZ(I)) .GT. HBS(3))  DZ(I) = DZ(I) - DSIGN(BS(3),DZ(I))
                       R2 = DX(I)*DX(I) + DY(I)*DY(I) + DZ(I)*DZ(I)
                       if(R2 .lt. DIST ) then
                          DIST = R2
                          NEARSITE(I) = J
                       end if
                    end do
                 end do

                 !--- to caluclate the jumping time and displacement
                 write(m_j_v_IOUnit, fmt="(A, 1x, 3(I8,1X),A,I10)") "&BOXID ", (JOB-1)*NBOX+IB, JOB, NBOX
                 write(m_j_v_IOUnit, fmt="(A, 1x, 3(I8,1X))") "&ITME ", Stamp%iTime, Stamp%ICFG(1), Stamp%IREC(1) 
                 write(m_j_v_IOUnit, fmt="(A, 1x, 3(I8,1X))") "&NDEF  ", NOLD, mp_MXI
                 do I=1, NOLD
                    DX(I)    = NEWPOS(NEARSITE(I),1)  -  OLDPOS(I,1); if(dabs(DX(I)) .GT. HBS(1))  DX(I) = DX(I) - DSIGN(BS(1),DX(I))
                    DY(I)    = NEWPOS(NEARSITE(I),2)  -  OLDPOS(I,2); if(dabs(DY(I)) .GT. HBS(2))  DY(I) = DY(I) - DSIGN(BS(2),DY(I))
                    DZ(I)    = NEWPOS(NEARSITE(I),3)  -  OLDPOS(I,3); if(dabs(DZ(I)) .GT. HBS(3))  DZ(I) = DZ(I) - DSIGN(BS(3),DZ(I))
                 
                    IS       = CHANGED(I)
                    write(m_j_v_IOUnit, fmt=FMT) &
                                              I, m_preTIMEi(IS),OLDSITE(I), TIME, NEWSITE(NEARSITE(I)),       &
                                              OLDPOS(I,1:3)/LATT, NEWPOS(NEARSITE(I),1:3)/LATT,               &
                                              DX(I)/LATT, DY(I)/LATT, DZ(I)/LATT,                             &
                                              (m_preOCCUai((IS-1)*mp_MXAI+J),                                 &
                                               m_preOCCUaP((IS-1)*mp_MXAI+J,1:3),J=1, mp_MXAI),               &
                                              (NEWOCCUai((NEARSITE(I)-1)*mp_MXAI+J),                          &
                                              NEWOCCUaP((NEARSITE(I)-1)*mp_MXAI+J,1:3),J=1, mp_MXAI)          

                    !--- update the current site
                    m_preSITEi(IS)    = NEWSITE(NEARSITE(I))
                    m_prePOSi(IS,1:3) = NEWPOS(NEARSITE(I), 1:3)
                    m_preTIMEi(IS)    = TIME
                    m_preOCCUai((IS-1)*mp_MXAI+1:IS*mp_MXAI)     = NEWOCCUai((NEARSITE(I)-1)*mp_MXAI+1:NEARSITE(I)*mp_MXAI)
                    m_preOCCUaP((IS-1)*mp_MXAI+1:IS*mp_MXAI,1:3) = NEWOCCUaP((NEARSITE(I)-1)*mp_MXAI+1:NEARSITE(I)*mp_MXAI,1:3)
                 end do
               end if  
           end do !--- end loop for boxes

           flush(m_j_v_IOUnit)
           m_NREC = m_NREC + 1

          return
  end subroutine DoEventCheck
  !**********************************************************************************

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
   
  !**********************************************************************************
  subroutine DoRecord(Stamp, SimBox, CtrlParam)
  use MD_ForceLib_Factory_GPU
  use MD_NeighborsList_GPU
  use MD_SimBoxArray_GPU
  use MD_SimboxArray
  use RefVoronoiVacancy_GPU, only:hm_RefSimBox
  implicit none
     !----   DUMMY Variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

     !----   Local variables
       integer::I, NB, hFILE, ERR, NBOX, NPRT0
       type(MDRecordStamp):: LSTAMP
       character*256::GFile, CFile, GFile0, GFileT
      
         !---------------
            if(gm_AppType(1:len_trim(gm_AppType)) .ne. "GMD") return
            if(CtrlParam%TimestpR .le. 0) return
 
            if(m_NREC .le. 0) then
               NPRT0 = hm_RefSimBox%NPRT
               NBOX  = size(SimBox)
               if(.not.allocated(m_preSITEi))  allocate(m_preSITEi(mp_MXI*NBOX))
               if(.not.allocated(m_prePOSi))   allocate(m_prePOSi(mp_MXI*NBOX,3))
               if(.not.allocated(m_preOCCUai)) allocate(m_preOCCUai(mp_MXI*mp_MXAI*NBOX))
               if(.not.allocated(m_preOCCUaP)) allocate(m_preOCCUaP(mp_MXI*mp_MXAI*NBOX,3))
               if(.not.allocated(m_preTIMEi))  allocate(m_preTIMEi(mp_MXI*NBOX))

               if(.not.allocated(m_hSITE))     allocate(m_hSITE(dm_NPRT))
               if(.not.allocated(m_hSTAT))     allocate(m_hSTAT(NPRT0*NBOX))
               if(.not.allocated(m_hOCCUa))    allocate(m_hOCCUa(NPRT0*NBOX*mp_MXAI))
               if(.not.allocated(m_hSITEi))    allocate(m_hSITEi(mp_MXI*NBOX))
               if(.not.allocated(m_hOCCUai))   allocate(m_hOCCUai(mp_MXI*mp_MXAI*NBOX))
               if(.not.allocated(m_hSitePOS))  allocate(m_hSitePOS(mp_MXI*NBOX,3))
               if(.not.allocated(m_hOCCUaP))   allocate(m_hOCCUaP(mp_MXI*mp_MXAI*NBOX,3))
               if(.not.allocated(m_SimBoxSwap))allocate(m_SimBoxSwap(NBOX))
            end if
 
            !---  Before to detetc event, we quench the current statu
            !---  NOTE: in Appshell, the boxes have been copyout from GPU. We do not 
            !           do copyout again
             call Copy_SimBoxArray(SimBox, m_SimBoxSwap)
             !m_SimBoxSwap = SimBox
             call Do_Damp(m_SimBoxSwap, CtrlParam)
             call CalEpot_ForceClass(m_SimBoxSwap,  CtrlParam, gm_ForceClass)
             call CopyOut_SimBox_DEV(m_SimBoxSwap)
             call DoEventCheck(Stamp, SimBox, CtrlParam) 
   
             !*** the statu on GPUs have been changed, we need to restore them
             !write(*,fmt="(A,I5, A, I7)") ' MDPSCU Message: Copy boxes back'
             call  CopyIn_SimBox_DEV(SimBox)
             !write(*,fmt="(A,I5, A, I7)") ' MDPSCU Message: Restore neighbor-list'
             call  Cal_NeighBoreList_DEV(SimBox, CtrlParam)
             return
   
      1000   write(*, *) "MDPSCU Error: in allocate memory in DoEventCheck of SIAEventSearch_module"             
  end subroutine DoRecord
  !*********************************************************************************

  !**********************************************************************************
  subroutine MyPutoutSummary(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       !--- dummy variables
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam
       !--- local variables
       type(SimMDBox)::SwapSimBox
       real(KINDDF)::DTIME, DX, DY, DZ, DR, WT, FT, COSDIR
       integer::I, J, hFile, event, nwalk, ibin, IB, NW, NF, IE
       character*128::tstr
       type(MDRecordStamp)::STAMP

              !--- calculate the histogram of event time
              rewind(m_j_v_IOUnit)

              read(m_j_v_IOUnit, *) tstr
              !do I=1, m_cumEVENTi1
              !   read(m_j_v_IOUnit1, *) IB, event, DTIME, COSDIR, DX, DY, DZ
              !end do

              close(m_j_v_IOUnit)
          return
  end subroutine MyPutoutSummary
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
  use MD_SimBoxArray
  use RefVoronoiVacancy_GPU, only:Clear_WorkingArray_DEV
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call MyPutoutSummary(SimBox, CtrlParam)

          call Clear_WorkingArray_DEV(SimBox, CtrlParam)

          if(allocated(m_SimBoxSwap)) then
            call Release_SimBoxArray(m_SimBoxSwap)
            deallocate(m_SimBoxSwap)
          endif

          if(allocated(m_preSITEi))  deallocate(m_preSITEi)
          if(allocated(m_prePOSi))   deallocate(m_prePOSi)
          if(allocated(m_preOCCUai)) deallocate(m_preOCCUai)
          if(allocated(m_preOCCUaP)) deallocate(m_preOCCUaP)
          if(allocated(m_preTIMEi))  deallocate(m_preTIMEi)

          if(allocated(m_hSITE))     deallocate(m_hSITE)
          if(allocated(m_hSTAT))     deallocate(m_hSTAT)
          if(allocated(m_hOCCUa))    deallocate(m_hOCCUa)
          if(allocated(m_hSITEi))    deallocate(m_hSITEi)
          if(allocated(m_hOCCUai))   deallocate(m_hOCCUai)
          if(allocated(m_hSitePOS))  deallocate(m_hSitePOS)
          if(allocated(m_hOCCUaP))   deallocate(m_hOCCUaP)

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
         m_NREC    = 0
         return
   end subroutine BeforeOneTest
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
      integer::NBOX, NPRT, NPRT0
              write(hFile) m_NREC
              if(m_NREC .le. 0) return

              NBOX  = size(m_preSITEi)/mp_MXI
              NPRT  = size(m_hSITE)
              NPRT0 = size(m_hSTAT)/NBOX
              write(hFile) NBOX, NPRT, NPRT0 
              write(hFile) m_preSITEi(1:mp_MXI*NBOX)
              write(hFile) m_prePOSi(1:mp_MXI*NBOX,1:3)
              write(hFile) m_preOCCUai(1:mp_MXI*mp_MXAI*NBOX)
              write(hFile) m_preOCCUaP(1:mp_MXI*mp_MXAI*NBOX,1:3)
              write(hFile) m_preTIMEi(1:mp_MXI*NBOX)


              return
  end subroutine Archive
  !**********************************************************************************

  !****************************************************************************
  subroutine Restore(hFile)
   !***  PURPOSE:   to save the current status of this module
   !
   !     INPUT:     hFile,  I/O unit number
   !
   !     OUTPUT
   implicit none
      !--- dummy varioables
      integer, intent(in)::hFile
 
      !--- local variables
      integer::NBOX, NPRT, NPRT0
              read(hFile) m_NREC
              if(m_NREC .le. 0) return
               
              write(hFile) NBOX, NPRT, NPRT0 
              if(.not.allocated(m_preSITEi))  allocate(m_preSITEi(mp_MXI*NBOX))
              if(.not.allocated(m_prePOSi))   allocate(m_prePOSi(mp_MXI*NBOX,3))
              if(.not.allocated(m_preOCCUai)) allocate(m_preOCCUai(mp_MXI*mp_MXAI*NBOX))
              if(.not.allocated(m_preOCCUaP)) allocate(m_preOCCUaP(mp_MXI*mp_MXAI*NBOX,3))
              if(.not.allocated(m_preTIMEi))  allocate(m_preTIMEi(mp_MXI*NBOX))

              if(.not.allocated(m_hSITE))     allocate(m_hSITE(NPRT))
              if(.not.allocated(m_hSTAT))     allocate(m_hSTAT(NPRT0*NBOX))
              if(.not.allocated(m_hOCCUa))    allocate(m_hOCCUa(NPRT0*NBOX*mp_MXAI))
              if(.not.allocated(m_hSITEi))    allocate(m_hSITEi(mp_MXI*NBOX))
              if(.not.allocated(m_hOCCUai))   allocate(m_hOCCUai(mp_MXI*mp_MXAI*NBOX))
              if(.not.allocated(m_hSitePOS))  allocate(m_hSitePOS(mp_MXI*NBOX,3))
              if(.not.allocated(m_hOCCUaP))   allocate(m_hOCCUaP(mp_MXI*mp_MXAI*NBOX,3))
              if(.not.allocated(m_SimBoxSwap))allocate(m_SimBoxSwap(NBOX))

              read(hFile) m_preSITEi(1:mp_MXI*NBOX)
              read(hFile) m_prePOSi(1:mp_MXI*NBOX,1:3)
              read(hFile) m_preOCCUai(1:mp_MXI*mp_MXAI*NBOX)
              read(hFile) m_preOCCUaP(1:mp_MXI*mp_MXAI*NBOX,1:3)
              read(hFile) m_preTIMEi(1:mp_MXI*NBOX)

              return
  end subroutine Restore
  !**********************************************************************************



 end module SIAEventSearch_module
 !**********************************************************************************
 !**********************************************************************************

 !****************************************************************************************
  Program EventSearch_GPU_Main
  use MD_SimBoxArray_AppShell_16_GPU, only:APPSHELL_AddRecord, APPSHELL_Main
  !---- If user-define potential to be used, use USE to get the entry to
  !     the register function of the potential.
  !     for example:
  !     use EAM_WHeH_ForceTable_Bonny_JPCM26_2014, only:Reg1=>Register_Interaction_Table
  !     use EM_TB_ForceTable_WangJun_W_HE_2010, only:Reg2=>Register_Interaction_Table

  use SIAEventSearch_module
  implicit none
  integer::numprocs=1, procid=0

       call APPSHELL_AddRecord(PRERECORD=Initialize,      &
                               BEFONETEST=BeforeOneTest,  &
                               RECORDPROC=DoRecord,       &
                               AFTRECORD=Clear) !,        &
                               !RESTARTREC=DoRestore,      &
                               !SAVERECSTATU=DoArchive)
       call APPSHELL_Main(numprocs,procid)

       stop
  end  program EventSearch_GPU_Main
