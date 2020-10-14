 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to calculate the vacancy diffusion coefficients.
 !                  A reference lattice matrix is to be loaded. The postion of an interstial or an vacancy is
 !                  identified by the occupation state of the lattice matrix. A lattice is identified as an
 !                  interstitial if it is occupied by two atom, is an vacancy if it is not occupied by an atom.
 !
 !
 !                  DEPENDENCE____________________________________________________________________________
 !                       MD_SimBoxArray_ToolShell_14_GPU.F90
 !                       RefVoronoiVA_GPU.F90
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !                  To use the analysis tool, the input control files need to be prepaired.
 !                  First, in the SETUP file used by the MD simulation (refer to the description of
 !                  the SETUP file for MD simulations), at least one line needs to be add:
 !
 !                    &AUXF_RVVIN  filename1
 !
 !                  where filename1 is the file that provide some control parameters. Optionally,
 !                  one can also add:
 !
 !                    &AUXF_RVVOUT filename2
 !
 !                  in the setup file. The filename2 is the file name for outputing results.
 !                  It should be noted that the filenames should be in quotated.
 !
 !                  The second file needs to be prepaired is the control file denoted by filename1.
 !                  A few kewords may appear in the control file:
 !
 !                  &REFCFG        indication of the filename of the reference configurqation.
 !                                  if the reference configuration is not assigned explicitly,
 !                                  the configuration at zero time step will be the reference
 !                                  configuration.
 !
 !                                  usage:   &REFCFG filename
 !
 !                                  example: &REFCFG "Myreference.crystal"
 !
 !                   &METHOD        indication of the method used to identify the occupation
 !                                  of the reference sites.
 !
 !                                  usage:   &METHOD = name rad
 !                                           where name is a string "USEVT" (default) or "USEST".
 !                                           If name = "USEVT", a site is indentified as occupied if there is
 !                                              one or more atoms is closer to the site than to other sites.
 !                                              the site is actually a Voronoi site.
 !                                           If name = "USEST", a site is identified as occupied if ther is
 !                                              an atom is in a sphere region of given raidu cenered by the
 !                                              site. rad is the raidu of the sphere region.
 !                                  example1: &METHOD "USEVT"
 !                                  example2: &METHOD "USEST"  rad
 !
 !                  If the filename of the reference configuration is given in the control file described
 !                  above, the third file storing the reference configuration needs to be prepaired. The
 !                  the file could be in the format of "&CFGXYZ", "&BOXCFG14" or "t-x-y-z" file whithout
 !                  header.
 !
 !                  With the input file(s) are ready, the program can be run on command line:
 !
 !                  Interstial_Vacancy_Diff_GPU.exe SETUP-file dev0  ndev
 !                  where:
 !                        SETUP-file  - the name of setup file used by MD simulations in MDPSCU.
 !                        dev0        - the ID of the first GPU to be used
 !                        ndev        - the number of GPUs to be used
 !
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2015-08 (Hou Qing, Sichuan university)
 !
 !

 module JumpRate_Vacancy_GPU
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none
         integer, private::m_NREC                                        ! number of actually record
         character(len=256),  private::m_OUTPATH =""                     ! the path output data
         integer, private::m_j_v_IOUnit  = 0                             ! temp I/O unit to storing jumping vectors
         integer, private::m_j_v_IOUnit1 = 0                             ! temp I/O unit to storing jumping vectors

         integer, parameter, private::mp_MXV = 1                         ! the permitted number of vacancy in a box
         integer, private::m_cumEVENTV                                   ! accumulation number of events of vcancy hop with direction changed
         integer, private::m_cumEVENTV1                                  ! accumulation number of events of vacancy hop
         integer, dimension(:), allocatable, private::m_curSITEv         ! cuurent sites of vacancy
         real(KINDDF), dimension(:,:), allocatable, private::m_curPOSv   ! cuurent position of vacancy lattice
         integer, dimension(:), allocatable, private::m_curNWalkv        ! number of walks before the moving direction changed
         integer, dimension(:), allocatable, private::m_curNFlightv      ! number of walks before the flght direction changed
         real(KINDDF), dimension(:),    allocatable, private::m_curWTimev! the current walking time
         real(KINDDF), dimension(:),    allocatable, private::m_curFTimev! the current fight time
         real(KINDDF), dimension(:,:),  allocatable, private::m_preDIRv  ! the previous moving direction of vacancy lattice

         real(KINDDF), dimension(:),    allocatable, private::m_preTIMEv ! the previous time for the occurance of event
         integer, parameter, private::mp_ntype_event         = 4
         integer, parameter, private::mp_event_walk          = 1
         integer, parameter, private::mp_event_dir_changed   = 2
         integer, parameter, private::mp_event_Z_changed     = 2**2
         integer, parameter, private::mp_event_inline_walk   = 2**3

         real(KINDDF), private::m_twBinI = 0.1                           ! width  of bins of jumping time histogram in ps
         integer,      private::m_tnBinI = 1000                          ! number of bins of jumping time histogram
         real(KINDDF), private::m_dwBinI = 0.01                          ! width  of bins of jumping distance histogram
         integer,      private::m_dnBinI = 1000                          ! number of bins of jumping distance histogram

         type(SimMDBox), dimension(:), allocatable::m_trackSimBox

  contains
  !**********************************************************************************
   subroutine MyInitialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use RefVoronoiVacancy_GPU, only:Initialize_ReferenceVacancy_DEV, Get_filename_Output

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables


         !$$--- loading control parameters
         call Initialize_ReferenceVacancy_DEV(SimBox, CtrlParam)
         m_cumEVENTv = 0
         m_cumEVENTv1= 0

         call Get_filename_Output(m_OUTPATH)
         call GetPath(m_OUTPATH, m_OUTPATH)
         call AvailableIOUnit(m_j_v_IOUnit)
         !open(m_j_v_IOUnit, FILE="status='SCRATCH')
         open(unit=m_j_v_IOUnit, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy.Time_Displace")

         call AvailableIOUnit(m_j_v_IOUnit1)
         !open(m_j_v_IOUnit, FILE="status='SCRATCH')
         open(unit=m_j_v_IOUnit1, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy.Time_Displace1")
         CtrlParam%TIMELOOPOUT = 1

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the clustering results to output file. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell_14_GPU.F90. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT,SynchronizeDevices, hm_ITYP, hm_XP
  use RefVoronoiVacancy_GPU, only:Cal_Occupied_Voronoi_DEV, hm_RefSimBox
  use MD_SimBoxArray
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        integer::I, J, IB, NBOX, IP, IA, NPRT0, NPRT, SHIFT, NV, IS, OLDSITE(mp_MXV), NEWSITE(mp_MXV), NEARSITE(mp_MXV), &
                 CHANGED(mp_MXV), EVENTTYP(mp_MXV), NOLD, NNEW, NGROUP0, ITIME
        integer, dimension(:),allocatable::hSITE, hSTAT, hSITEv
        real(KINDDF), dimension(:,:),allocatable::hSitePOS
        real(KINDDF)::OLDPOS(mp_MXV,3), NEWPOS(mp_MXV,3), PRETIME(mp_MXV),   &
                      DX(mp_MXV), DY(mp_MXV), DZ(mp_MXV), DR(mp_MXV), TIME,  &
                      DTIME(mp_MXV), DIST, R2, LATT, BS(3), HBS(3), COSDIR

  !-------------
           ITime = Stamp%ITime
           Time  = Stamp%Time
           if(ITime .lt. 0) then
              m_NREC = 0
              return
           end if

           NPRT0 = hm_RefSimBox%NPRT
           NGROUP0 = hm_RefSimBox%NGROUP
           NBOX  = size(SimBox)
           NPRT  = dm_NPRT/NBOX
           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           HBS   = C_HALF*BS

           if(.not. allocated(m_curPOSv)) then
              allocate(m_curSITEv(mp_MXV*NBOX), m_curPOSv(mp_MXV*NBOX,3))
              allocate(m_preTIMEv(mp_MXV*NBOX), m_preDIRv(mp_MXV*NBOX,3), m_curNWalkv(mp_MXV*NBOX),m_curNFlightv(mp_MXV*NBOX) )
              allocate(m_curWTimev(mp_MXV*NBOX),m_curFTimev(mp_MXV*NBOX) )

              m_curSITEv  = 0
              m_preTIMEv  = 0.D0
              m_preDIRv   = 0.D0
              m_curNWalkv = 0
              m_curNFlightv= 0
              m_curWTimev = 0.D0
              m_curFTimev = 0.D0
           end if

           if(.not. allocated(m_trackSimBox)) then
              call Create_SimBoxArray(m_trackSimBox, NBOX, hm_RefSimBox)
           end if

         !--- calculate the occupation of lattice sites
           allocate(hSITE(dm_NPRT), hSTAT(NPRT0*NBOX), hSITEv(mp_MXV*NBOX), hSitePOS(mp_MXV*NBOX,3))
           call Cal_Occupied_Voronoi_DEV(SimBox(1), CtrlParam, hSITE)

          !--- to calculate occupation state
           hSTAT = 0
           IP = 0
           do IB=1, NBOX
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT
                 IP = IP + 1
                 if(hSITE(IP)>0) then
                    hSTAT(hSITE(IP)+SHIFT) = hSTAT(hSITE(IP)+SHIFT) + 1
                 end if
              end do
           end do

          !--- to extract the postion of vacancy
           IP = 0
           hSITEv =  0
           do IB=1, NBOX
              NV = 0
              do I=1, NPRT0
                 IP = IP + 1
                 if(hSTAT(IP) .le. 0) then
                    NV = NV + 1
                    if(NV .gt. mp_MXV) exit

                    IS = (IB-1)*mp_MXV + NV
                    hSITEv(IS) = I
                    hSitePOS(IS,1:3) = hm_RefSimBox%XP(I, 1:3)
                 end if
              end do
              IS = (IB-1)*mp_MXV
              call AddAtoms_SimMDBox(m_trackSimBox(IB), NV, ITYPE=NGROUP0+1, RXP=hSitePOS(IS+1:IS+NV,1:3), TYPEORDER=0)

              if(NV .gt. mp_MXV) then
                 write(*,fmt="(A, 1x, I4, 1x, A, I4, A, I4)") &
                      " MDPSCU Warning: the number of vacancy in box #", IB, " larger than permitted value:", NV, " vs ",mp_MXV
                 call ONWARNING(gm_OnWarning)
               end if
           end do

           if(m_NREC .eq. 0 ) then
              m_curSITEv  = hSITEv
              m_curPOSv   = hSitePOS
              m_preTIMEv  = TIME
              m_preDIRv   = 0.D0
              m_curNWalkv = 0
              m_curNFlightv= 0
              m_curWTimev = 0.D0
              m_curFTimev = 0.D0
              deallocate(hSITE, hSitePOS, hSTAT, hSITEv)
              m_NREC = m_NREC + 1
              return
           end if
           !print *, "ITIME, TIME", ITIME, TIME, TIME-m_preTIMEv(1)

          !--- to check if exchange event occur
           do IB=1, NBOX
              OLDSITE = 0
              NEWSITE = 0
              !$$--- to find out which site has been changed
              NOLD = 0
              do I=1, mp_MXV
                 IS = (IB-1)*mp_MXV+I
                 if(.not.any(hSITEv((IB-1)*mp_MXV+1:IB*mp_MXV) .eq. m_curSITEv(IS)) ) then
                    NOLD = NOLD + 1
                    OLDSITE(NOLD)    = m_curSITEv(IS)
                    OLDPOS(NOLD,1:3) = m_curPOSv(IS, 1:3)
                    PRETIME(NOLD)    = m_preTIMEv(IS)
                    CHANGED(NOLD)    = IS
                 end if
              end do
              !--- to find out new sites not existing in old sites
              NNEW = 0
              do I=1, mp_MXV
                 IS = (IB-1)*mp_MXV+I
                 if(.not.any(m_curSITEv((IB-1)*mp_MXV+1:IB*mp_MXV) .eq. hSITEv(IS)) ) then
                    NNEW = NNEW + 1
                    NEWSITE(NNEW) = hSITEv(IS)
                    NEWPOS(NNEW,1:3) = hSitePOS(IS, 1:3)
                 end if
              end do

              if(NNEW .ne. NOLD) then
                 write(*,fmt="(A, ix, I4, 1x, I4, 1x, A, 1x, I4)")   " MDPSCU Warning: number of vacancy is changed in box  ", &
                                                                     IB, NNEW, " vs ", NOLD
                 call ONWARNING(gm_OnWarning)
              end if

             !--- to find out the neasest new sites from old sites
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
              do I=1, NOLD
                 DTIME(I) = TIME - PRETIME(I)
                 DX(I)    = NEWPOS(NEARSITE(I),1)  -  OLDPOS(I,1); if(dabs(DX(I)) .GT. HBS(1))  DX(I) = DX(I) - DSIGN(BS(1),DX(I))
                 DY(I)    = NEWPOS(NEARSITE(I),2)  -  OLDPOS(I,2); if(dabs(DY(I)) .GT. HBS(2))  DY(I) = DY(I) - DSIGN(BS(2),DY(I))
                 DZ(I)    = NEWPOS(NEARSITE(I),3)  -  OLDPOS(I,3); if(dabs(DZ(I)) .GT. HBS(3))  DZ(I) = DZ(I) - DSIGN(BS(3),DZ(I))
                 DR(I)    = dsqrt(DX(I)*DX(I) + DY(I)*DY(I) + DZ(I)*DZ(I))

                 !--- check the type of events
                 IS = CHANGED(I)
                 if(all(dabs(m_preDIRv(IS,1:3)/LATT) .le.  1.D-5) ) then
                    EVENTTYP(I)      = mp_event_walk
                    if(dabs(DZ(I))/LATT .le. 0.01) then
                       EVENTTYP(I) = IOR(EVENTTYP(I), mp_event_Z_changed)
                    else
                       m_curNWalkv(IS)  = m_curNWalkv(IS)  + 1
                       m_curNFlightv(IS) = m_curNFlightv(IS) + 1
                       m_curWTimev(IS) = m_curWTimev(IS) + DTIME(I)
                       m_curFTimev(IS) = m_curFTimev(IS) + DTIME(I)
                    end if
                    COSDIR = 2
                 else
                   COSDIR = (m_preDIRv(IS,1)*DX(I) + m_preDIRv(IS,2)*DY(I) + m_preDIRv(IS,3)*DZ(I))/   &
                             dsqrt(DX(I)*DX(I) + DY(I)*DY(I) + DZ(I)*DZ(I))/dsqrt(sum(m_preDIRv(IS,1:3)*m_preDIRv(IS,1:3)))
                   if(COSDIR .gt. 0.99D0 ) then
                      EVENTTYP(I) = mp_event_walk
                      m_curNWalkv(IS)  = m_curNWalkv(IS)  + 1
                      m_curNFlightv(IS) = m_curNFlightv(IS) + 1
                      m_curWTimev(IS)  = m_curWTimev(IS) + DTIME(I)
                      m_curFTimev(IS)  = m_curFTimev(IS) + DTIME(I)
                      !print *,"forward ", DX(I), DY(I), DZ(I)
                      !pause

                   else
                      EVENTTYP(I) = mp_event_dir_changed
                      if(dabs(DZ(I))/LATT .gt. 0.01) then
                         EVENTTYP(I) = IOR(EVENTTYP(I), mp_event_Z_changed)
                         !print *,"Z change ", DX(I), DY(I), DZ(I)
                         !pause

                      else if(COSDIR .lt. -0.99D0 ) then
                         EVENTTYP(I) = IOR(EVENTTYP(I), mp_event_inline_walk)
                         m_curNWalkv(IS)  = m_curNWalkv(IS)  + 1
                         m_curWTimev(IS)  = m_curWTimev(IS)+ DTIME(I)
                         m_curNFlightv(IS) = 0
                         m_curFTimev(IS)  = 0.D0
                      else

                      end if
                   end if
                   !print *, "COSDIR",COSDIR, "m_EVENTTYPi ",  EVENTTYP(I)
                   !pause
                 end if

                !--- up date the
                 m_cumEVENTv1 = m_cumEVENTv1 + 1
                 if(EVENTTYP(I) .ne. mp_event_walk) then
                   if(IAND(EVENTTYP(I),mp_event_inline_walk).ne.mp_event_inline_walk) then
                    write(m_j_v_IOUnit, fmt="(1x, I8, 1x, I8, 1x, 6(1x, I8, 1PE13.4,1x))") IB, m_cumEVENTv1, EVENTTYP(I),   DTIME(I),   &
                                                                                           m_curNWalkv(IS),  m_curWTimev(IS)-DTIME(I),  &
                                                                                           m_curNFlightv(IS), m_curFTimev(IS)-DTIME(I)

                    m_curNWalkv(IS)   = 0
                    m_curNFlightv(IS)  = 0
                    m_curWTimev(IS)   = 0.D0
                    m_curFTimev(IS)   = 0.D0
                    m_cumEVENTv = m_cumEVENTv + 1
                   end if
                 end if

                 write(*,  fmt="(3(1x, I8), 6(1x, 1PE13.4,1x))") IB, ITIME, EVENTTYP(I), DTIME(I), COSDIR
                 write(m_j_v_IOUnit1, fmt="(3(1x, I8), 6(1x, 1PE13.4,1x))") IB, ITIME, EVENTTYP(I), DTIME(I), COSDIR

                 !$$--- update the moving direction, current site
                 if(iand(EVENTTYP(I), mp_event_Z_changed) .gt. 0) then
                    m_preDIRv(IS,1:3) = 0.D0
                 else
                    m_preDIRv(IS,1) = DX(I)
                    m_preDIRv(IS,2) = DY(I)
                    m_preDIRv(IS,3) = DZ(I)
                 end if

                 m_curSITEv(IS)    = NEWSITE(NEARSITE(I))
                 m_curPOSv(IS,1:3) = NEWPOS(NEARSITE(I), 1:3)
                 m_preTIMEv(IS)    = TIME

              end do !--- end loop for NOLD
           end do !--- end loop for boxes

           deallocate(hSITE, hSitePOS, hSTAT, hSITEv)
           m_NREC = m_NREC + 1

          return
  end subroutine MyRECORD
  !**********************************************************************************

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
       integer::I, J, hFile, event, nwalk, ibin, IB, NW, NF, IE, ITIME
       character*32::tstr
       integer::nWE, nFE, nZCE, nDCE
       integer, dimension(:), allocatable::tHis_inline, tHis_Zchange, tHis_Dir, tHis
       integer, dimension(:), allocatable::nwHis,nfHis,twHis,tfHis
       integer, dimension(:), allocatable::dHis_inline, dHis_Zchange, dHis_Dir, dHis
       integer, dimension(:), allocatable::wHis_inline, wHis_Zchange, wHis_Dir, wHis
       type(MDRecordStamp)::STAMP


              allocate(tHis_inline(m_tnBinI), tHis_Zchange(m_tnBinI), tHis_Dir(m_tnBinI), tHis(m_dnBinI) )
              allocate(nwHis(4000), nfHis(4000), twHis(m_tnBinI), tfHis(m_tnBinI))
              allocate(dHis_inline(m_dnBinI), dHis_Zchange(m_dnBinI), dHis_Dir(m_dnBinI), dHis(m_dnBinI) )
              allocate(wHis_inline(200), wHis_Zchange(200), wHis_Dir(200), wHis(200) )
              tHis_inline  = 0
              tHis_Zchange = 0
              tHis_Dir     = 0
              tHis         = 0
              dHis_inline  = 0
              dHis_Zchange = 0
              dHis_Dir     = 0
              dHis         = 0
              wHis_inline  = 0
              wHis_Zchange = 0
              wHis_Dir     = 0
              wHis         = 0
              nWE  = 0
              nFE  = 0
              nZCE = 0
              nDCE = 0

              !--- calculate the histogram of event time
              rewind(m_j_v_IOUnit)
              rewind(m_j_v_IOUnit1)

              do I=1, m_cumEVENTv1
                   read(m_j_v_IOUnit1, *) IB, ITIME, event, DTIME, COSDIR

                  ibin = int(DTIME/m_twBinI) + 1
                  if(event .eq. mp_event_walk) then
                     tHis(ibin) = tHis(ibin) + 1
                     nFE = nFE + 1
                  else
                     if(iand(event, mp_event_inline_walk) .eq. mp_event_inline_walk) then
                        tHis_inline(ibin) = tHis_inline(ibin) + 1
                        nWE = nWE + 1
                     else  if(iand(event, mp_event_Z_changed) .eq. mp_event_Z_changed) then
                        tHis_Zchange(ibin) = tHis_Zchange(ibin) + 1
                        nZCE = nZCE + 1
                     else
                        tHis_Dir(ibin) = tHis_Dir(ibin) + 1
                        nDCE = nDCE + 1
                     end if
                  end if
              end do

              !-- write out histogram
              call AvailableIOUnit(hFile)
              open(unit=hFile, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"Interstitial_Jump_Time.his")
              write(*, fmt="(A)") "JumpTime-histogram to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Interstitial_Jump_Time.his"

              write(hFile, fmt="(A, 1x, 8(I8,1x))") "!*** the number of events of forward jump:  ", nFE
              write(hFile, fmt="(A, 1x, 8(I8,1x))") "!*** the number of events of backward jump: ", nWE
              write(hFile, fmt="(A, 1x, 8(I8,1x))") "!*** the number of events of Z jump:  ", nZCE
              write(hFile, fmt="(A, 1x, 8(I8,1x))") "!*** the number of events of dir-change jump: ", nDCE
              do I=1, m_tnBinI
                 write(hFile, fmt="(1x, 1PE13.4,1x, 8(I8,1x))") m_twBinI*I, tHis(I), tHis_inline(I), tHis_Zchange(I), tHis_Dir(I), &
                                                                            tHis(I)+tHis_inline(I)+tHis_Zchange(I)+tHis_Dir(I)
              end do
              close(hFile)

              !--- calculate the number of walks before direction change
              nwHis = 0
              nfHis = 0
              tHis_Zchange = 0
              tHis_Dir = 0
              do I=1, m_cumEVENTv
                  read(m_j_v_IOUnit, *) IB, IE, event, DTIME, NW, WT, NF, FT


                  if(NW.gt.0) nwHis(NW) = nwHis(NW) + 1
                  if(NF.gt.0) nfHis(NF) = nfHis(NF) + 1

                  ibin = int(WT/m_twBinI) + 1
                  twHis(ibin) = twHis(ibin) + 1

                  if(NF .gt. 0) then
                     ibin = int(FT/m_twBinI) + 1
                     tfHis(ibin) = tfHis(ibin) + 1
                   end if
              end do

              !-- write out histogram
              call AvailableIOUnit(hFile)
              open(unit=hFile, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"Interstitial_walk_number.his")
              write(*, fmt="(A)") "Interstitial_walk number before change direction-histogram to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Interstitial_walk_number.his"
              do I=1, size(nwHis)
                 write(hFile, fmt="(1x, 8(I8,1x))") nwHis(I), nfHis(I)
              end do
              close(hFile)

              write(*, fmt="(A)") "Move change events historty to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy.Time_Displace"
              write(*, fmt="(A)") "Events historty to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy.Time_Displace1"
              close(m_j_v_IOUnit)
              close(m_j_v_IOUnit1)

            !write(*, fmt="(A)") "Jumping vector to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy_Jump_Vector"
            !call Putout_Instance_Config_SimMDBox(0, 1, 0, 0.D0, m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy_Jump_Vector", swapSimBox)
            !call Release_SimMDBox(SwapSimBox)

            print *, size(m_trackSimBox)
            do I=1, size(m_trackSimBox)
               write(tstr,*) I
               tstr = adjustl(tstr)
               write(*, fmt="(A)") "Jumping track to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy_Jump_Track"//tstr(1:len_trim(tstr))
               STAMP%ITime  = 0
               STAMP%ISect  = 1
               STAMP%ICfg   = 0
               STAMP%IRec   = 0
               STAMP%Time   = 0
               call Putout_Instance_Config_SimMDBox(m_OUTPATH(1:len_trim(m_OUTPATH))//"Vacancy_Jump_Track"&
                                //tstr(1:len_trim(tstr)), m_trackSimBox(I), STAMP)
            end do


            deallocate(tHis_inline, tHis_Zchange, tHis_Dir, tHis )
            deallocate(nwHis, nfHis, twHis, tfHis)
            deallocate(dHis_inline, dHis_Zchange, dHis_Dir, dHis )
            deallocate(wHis_inline, wHis_Zchange, wHis_Dir, wHis )

          return
  end subroutine MyPutoutSummary
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyCleaner(SimBox, CtrlParam)
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
          if(allocated(m_curSITEv))  deallocate(m_curSITEv)
          if(allocated(m_curPOSv))   deallocate(m_curPOSv)
          if(allocated(m_preTIMEv))  deallocate(m_preTIMEv)
          if(allocated(m_preDIRv))   deallocate(m_preDIRv)
          if(allocated(m_curNWalkv)) deallocate(m_curNWalkv)
          if(allocated(m_curNFlightv)) deallocate(m_curNFlightv)
          if(allocated(m_curWTimev)) deallocate(m_curWTimev)
          if(allocated(m_curFTimev)) deallocate(m_curFTimev)


          !if(allocated(m_tHISi))     deallocate(m_tHISi)
          !if(allocated(m_dHISix))    deallocate(m_dHISix)
          !if(allocated(m_dHISiy))    deallocate(m_dHISiy)
          !if(allocated(m_dHISiz))    deallocate(m_dHISiz)
          !if(allocated(m_dHISi))     deallocate(m_dHISi)
          call Release_SimBoxArray(m_trackSimBox)

          m_NREC = 0
          if(m_j_v_IOUnit.gt.0) close(m_j_v_IOUnit)

          return
  end subroutine MyCleaner
  !**********************************************************************************

 end module JumpRate_Vacancy_GPU


 !**********************************************************************************
 Program Diff_Vacancy_GPU_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use JumpRate_Vacancy_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD=MyInitialize, &
                                RECORDPROC=MyRECORD,    &
                                AFTRECORD= MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program Diff_Vacancy_GPU_Main
