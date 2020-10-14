 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to calculate the jump rate of interstial on BCC lattice.
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
 !                  version 1st 2015-09 (Hou Qing, Sichuan university)
 !
 !

 module JumpRate_SIA_bcc_GPU
 use MD_CONSTANTS
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 implicit none
         integer, private::m_NREC                                        ! number of actually record
         integer, private::m_NRECTEST=0                                  ! number of actually boxes, the boxes with
                                                                         ! atom out-of-box are filtered

         integer, dimension(:), allocatable,private::m_BOXFLAG           ! flag indicating if there are atom out of box
         character(len=256),  private::m_OUTPATH =""                     ! the path output data
         integer, private::m_j_v_IOUnit  = 0                             ! temp I/O unit to storing jumping vectors
         integer, private::m_j_v_IOUnit1 = 0                             ! temp I/O unit to storing jumping vectors

         integer, parameter, private::mp_MXI = 1                         ! the permitted number of interstital in a box
         integer, parameter, private::mp_MXAI =4                         ! the permitted number of atoms occupy one site
         integer, private::m_cumEVENTi                                   ! accumulation number of events of interstial hop with direction changed
         integer, private::m_cumEVENTi1                                  ! accumulation number of events of interstial hop
         integer, dimension(:), allocatable, private::m_curSITEi         ! cuurent sites of interstitials
         integer, dimension(:), allocatable, private::m_curOCCUai        ! cuurent atom id that occupy a site
         real(KINDDF), dimension(:,:), allocatable, private::m_curPOSi   ! cuurent position of interstitial lattice
         integer, dimension(:), allocatable, private::m_curNWalki        ! number of walks before the moving direction changed
         integer, dimension(:), allocatable, private::m_curNFlighti      ! number of walks before the flght direction changed
         real(KINDDF), dimension(:),    allocatable, private::m_curWTimei! the current walking time
         real(KINDDF), dimension(:),    allocatable, private::m_curFTimei! the current fight time
         real(KINDDF), dimension(:,:),  allocatable, private::m_preDIRi  ! the previous moving direction of interstitil lattice

         real(KINDDF), dimension(:),    allocatable, private::m_preTIMEi ! the previous time for the occurance of event
         integer, parameter, private::mp_ntype_event         = 4
         integer, parameter, private::mp_event_walk          = 1
         integer, parameter, private::mp_event_dir_changed   = 2
         integer, parameter, private::mp_event_inline_walk   = 2**3

         real(KINDDF), private::m_twBinI = 0.1                           ! width  of bins of jumping time histogram in ps
         integer,      private::m_tnBinI = 1000                          ! number of bins of jumping time histogram
         real(KINDDF), private::m_dwBinI = 0.01                          ! width  of bins of jumping distance histogram
         integer,      private::m_dnBinI = 1000                          ! number of bins of jumping distance histogram

         !integer,      dimension(:), allocatable, private::m_tHISi       ! histogram of jumping time
         !integer,      dimension(:), allocatable, private::m_dHISix      ! histogram of jumping distance histogram
         !integer,      dimension(:), allocatable, private::m_dHISiy      ! histogram of jumping distance histogram
         !integer,      dimension(:), allocatable, private::m_dHISiz      ! histogram of jumping distance histogram
         !integer,      dimension(:), allocatable, private::m_dHISi       ! histogram of jumping distance histogram
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
         m_cumEVENTi = 0
         m_cumEVENTi1= 0

         call Get_filename_Output(m_OUTPATH)
         call GetPath(m_OUTPATH, m_OUTPATH)
         call AvailableIOUnit(m_j_v_IOUnit)
         !open(m_j_v_IOUnit, FILE="status='SCRATCH')
         open(unit=m_j_v_IOUnit, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"Interstial.Time_Displace")

         call AvailableIOUnit(m_j_v_IOUnit1)
         !open(m_j_v_IOUnit, FILE="status='SCRATCH')
         open(unit=m_j_v_IOUnit1, FILE=m_OUTPATH(1:len_trim(m_OUTPATH))//"Interstial.Time_Displace1")
          write(m_j_v_IOUnit1, fmt="(1x, A8, 1x, A9, 16(1x, A13,1x))") &
               "BOX #", "EVENTTYPE", "Delta_t", "COS(ang)", "DX", "DY", "DZ", "DR"

         m_NRECTEST = 0
         CtrlParam%TIMELOOPOUT = 1

         return
   end subroutine MyInitialize
  !**********************************************************************************

  !**********************************************************************************
  subroutine MyRECORD(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION:
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
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
        integer::I, J, IB, NBOX, IP, IA, NA, NPRT0, NPRT, SHIFT, NI, IS, OLDSITE(mp_MXI), NEWSITE(mp_MXI), NEARSITE(mp_MXI), &
                 CHANGED(mp_MXI), EVENTTYP(mp_MXI), NOLD, NNEW, JOB
        integer, dimension(:),allocatable::hSITE, hSTAT, hOCCUa, hSITEi, hOCCUai
        real(KINDDF), dimension(:,:),allocatable::hSitePOS
        real(KINDDF)::OLDPOS(mp_MXI,3), NEWPOS(mp_MXI,3), PRETIME(mp_MXI), APOS(mp_MXI*mp_MXAI, 3), AVEC(mp_MXI*mp_MXAI, 3), &
                      DX(mp_MXI), DY(mp_MXI), DZ(mp_MXI), DR(mp_MXI), TIME, &
                      DTIME(mp_MXI), DIST, R2, LATT, BS(3), HBS(3), COSDIR


  !-------------
           if(Stamp%ITime .lt. 0) then
              m_NREC = 0
              return
           end if
           JOB   = Stamp%ITest
           TIME  = Stamp%ITime

           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(SimBox)
           NPRT  = dm_NPRT/NBOX
           LATT  = hm_RefSimBox%RR
           BS    = hm_RefSimBox%ZL
           HBS   = C_HALF*BS

           if(.not. allocated(m_curPOSi)) then
              allocate(m_curSITEi(mp_MXI*NBOX), m_curPOSi(mp_MXI*NBOX,3), m_curOCCUai(mp_MXI*mp_MXAI*NBOX))
              allocate(m_preTIMEi(mp_MXI*NBOX), m_preDIRi(mp_MXI*NBOX,3), m_curNWalki(mp_MXI*NBOX),m_curNFlighti(mp_MXI*NBOX) )
              allocate(m_curWTimei(mp_MXI*NBOX),m_curFTimei(mp_MXI*NBOX) )
              allocate(m_BOXFLAG(NBOX))
              m_BOXFLAG = 0
           end if

           if(.not. allocated(m_trackSimBox)) then
              call Create_SimBoxArray(m_trackSimBox, NBOX, hm_RefSimBox)
           end if

         !--- calculate the occupation of lattice sites
           allocate(hSITE(dm_NPRT), hSTAT(NPRT0*NBOX), hOCCUa(NPRT0*NBOX*mp_MXAI), &
                    hSITEi(mp_MXI*NBOX), hOCCUai(mp_MXI*mp_MXAI*NBOX), hSitePOS(mp_MXI*NBOX,3))
           call Cal_Occupied_Voronoi_DEV(SimBox(1), CtrlParam, hSITE)

          !--- to calculate occupation state
           hSTAT = 0
           hOCCUa = 0
           IP = 0
           do IB=1, NBOX
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT
                 IP = IP + 1
                 if(hSITE(IP)>0) then
                    hSTAT(hSITE(IP)+SHIFT) = hSTAT(hSITE(IP)+SHIFT) + 1
                    do J=1, mp_MXAI
                       !--- hOCCUa: the atoms ID on site SHIFT + hSITE(IP)
                       if(hOCCUa((SHIFT + hSITE(IP)-1)*mp_MXAI + J).eq. 0) then
                          hOCCUa((SHIFT + hSITE(IP)-1)*mp_MXAI + J) = IP
                          exit
                       end if
                    end do
                 end if
              end do
           end do

          !--- to extract the postion of interstitial
           IP = 0
           hSITEi =  0
           hOCCUai = 0
           do IB=1, NBOX
              NI = 0
              NA = 0
              do I=1, NPRT0
                 IP = IP + 1
                 if(hSTAT(IP) .ge. 2) then
                    !--- this site containing interstitial atoms
                    NI = NI + 1
                    if(NI .gt. mp_MXI) exit
                    IS = (IB-1)*mp_MXI+NI
                    hSITEi(IS) = I
                    hSitePOS(IS,1:3) = hm_RefSimBox%XP(I, 1:3)

                    do J=1, mp_MXAI
                       IA = hOCCUa((IP-1)*mp_MXAI+J)
                       if(IA .gt. 0) then
                          NA = NA + 1
                          hOCCUai((IS-1)*mp_MXAI+J) = IA
                          APOS(NA,1:3) = hm_XP(IA, 1:3)
                          AVEC(NA,1:3) = hSitePOS(IS,1:3) - hm_XP(IA, 1:3)
                       end if
                    end do
                 end if
              end do
              IS = (IB-1)*mp_MXI
              call AddAtoms_SimMDBox(m_trackSimBox(IB), NI, ITYPE=2, RXP=hSitePOS(IS+1:IS+NI,1:3), TYPEORDER=0)
              do J=1, NA
                 call AddAtoms_SimMDBox(m_trackSimBox(IB), 1, ITYPE=2+J, RXP=APOS(J:J,1:3), RXP1=AVEC(J:J,1:3), TYPEORDER=0)
              end do

              if(NI .gt. mp_MXI) then
                 write(*,fmt="(A, 1x, I4, 1x, A, I4, A, I4)") &
                      " MDPSCU Warning: the number of interstitials in box #", IB, " larger than permitted value:", NI, " vs ",mp_MXI
                 call ONWARNING(gm_OnWarning)
               end if
           end do

           if(m_NREC .eq. 0 ) then
              m_curSITEi  = hSITEi
              m_curPOSi   = hSitePOS
              m_curOCCUai = hOCCUai
              m_preTIMEi  = TIME
              m_preDIRi   = 0.D0
              m_curNWalki = 0
              m_curNFlighti= 0
              m_curWTimei = 0.D0
              m_curFTimei = 0.D0
              m_NREC = m_NREC + 1
              !$$--- check if there atom are out of box
              m_BOXFLAG  = 0
              do IB=1, NBOX
                 if(any(iand(SimBox(IB)%STATU, CP_STATU_OUTOFBOX) .eq. CP_STATU_OUTOFBOX)) then
                    m_BOXFLAG(IB) = 1
                  end if
              end do
              m_NRECTEST = m_NRECTEST + count(m_BOXFLAG.le.0)
              deallocate(hSITE, hSitePOS, hSTAT, hOCCUa, hSITEi, hOCCUai)
              return
           end if

          !--- to check if exchange event occur
           do IB=1, NBOX
              if(m_BOXFLAG(IB) .gt. 0) cycle

              OLDSITE = 0
              NEWSITE = 0
              !$$--- to find out which site has been changed
              NOLD = 0
              do I=1, mp_MXI
                 IS = (IB-1)*mp_MXI+I
                 if(.not.any(hSITEi((IB-1)*mp_MXI+1:IB*mp_MXI) .eq. m_curSITEi(IS)) ) then
                    NOLD = NOLD + 1
                    OLDSITE(NOLD)    = m_curSITEi(IS)
                    OLDPOS(NOLD,1:3) = m_curPOSi(IS, 1:3)
                    PRETIME(NOLD)    = m_preTIMEi(IS)
                    CHANGED(NOLD)    = IS
                 end if
              end do
              !--- to find out new sites not existing in old sites
              NNEW = 0
              do I=1, mp_MXI
                 IS = (IB-1)*mp_MXI+I
                 if(.not.any(m_curSITEi((IB-1)*mp_MXI+1:IB*mp_MXI) .eq. hSITEi(IS)) ) then
                    NNEW = NNEW + 1
                    NEWSITE(NNEW) = hSITEi(IS)
                    NEWPOS(NNEW,1:3) = hSitePOS(IS, 1:3)
                 end if
              end do

              if(NNEW .ne. NOLD) then
                 write(*,fmt="(A, ix, I4, 1x, I4, 1x, A, 1x, I4)")   " MDPSCU Warning: number of interstitials is changed in box  ", &
                                                                     (JOB-1)*NBOX+IB, NNEW, " vs ", NOLD
                 call ONWARNING(gm_OnWarning)
              end if

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
              do I=1, NOLD
                 DTIME(I) = TIME - PRETIME(I)
                 DX(I)    = NEWPOS(NEARSITE(I),1)  -  OLDPOS(I,1); if(dabs(DX(I)) .GT. HBS(1))  DX(I) = DX(I) - DSIGN(BS(1),DX(I))
                 DY(I)    = NEWPOS(NEARSITE(I),2)  -  OLDPOS(I,2); if(dabs(DY(I)) .GT. HBS(2))  DY(I) = DY(I) - DSIGN(BS(2),DY(I))
                 DZ(I)    = NEWPOS(NEARSITE(I),3)  -  OLDPOS(I,3); if(dabs(DZ(I)) .GT. HBS(3))  DZ(I) = DZ(I) - DSIGN(BS(3),DZ(I))
                 DR(I)    = dsqrt(DX(I)*DX(I) + DY(I)*DY(I) + DZ(I)*DZ(I))
                 !--- check the type of events
                 IS = CHANGED(I)
                 if(all(dabs(m_preDIRi(IS,1:3)/LATT) .le.  1.D-5) ) then
                    EVENTTYP(I)      = mp_event_walk
                    m_curNWalki(IS)  = m_curNWalki(IS)  + 1
                    m_curNFlighti(IS) = m_curNFlighti(IS) + 1
                    m_curWTimei(IS) = m_curWTimei(IS) + DTIME(I)
                    m_curFTimei(IS) = m_curFTimei(IS) + DTIME(I)
                    COSDIR = 2
                 else
                   COSDIR = (m_preDIRi(IS,1)*DX(I) + m_preDIRi(IS,2)*DY(I) + m_preDIRi(IS,3)*DZ(I))/   &
                             dsqrt(DX(I)*DX(I) + DY(I)*DY(I) + DZ(I)*DZ(I))/dsqrt(sum(m_preDIRi(IS,1:3)*m_preDIRi(IS,1:3)))
                   if(COSDIR .gt. 0.99D0 ) then
                      EVENTTYP(I) = mp_event_walk
                      m_curNWalki(IS)  = m_curNWalki(IS)  + 1
                      m_curNFlighti(IS) = m_curNFlighti(IS) + 1
                      m_curWTimei(IS)  = m_curWTimei(IS) + DTIME(I)
                      m_curFTimei(IS)  = m_curFTimei(IS) + DTIME(I)

                   else
                      EVENTTYP(I) = mp_event_dir_changed
                      if(COSDIR .lt. -0.99D0 ) then
                         EVENTTYP(I) = IOR(EVENTTYP(I), mp_event_inline_walk)
                         m_curNWalki(IS)  = m_curNWalki(IS)  + 1
                         m_curWTimei(IS)  = m_curWTimei(IS)+ DTIME(I)
                         m_curNFlighti(IS) = 0
                         m_curFTimei(IS)  = 0.D0
                      end if
                   end if
                 end if

                !--- up date the
                 m_cumEVENTi1 = m_cumEVENTi1 + 1
                 if(EVENTTYP(I) .ne. mp_event_walk) then
                   if(IAND(EVENTTYP(I),mp_event_inline_walk).ne.mp_event_inline_walk) then
                    write(m_j_v_IOUnit, fmt="(1x, I8, 1x, I9, 1x, 6(1x, I8, 1PE13.4,1x))") (JOB-1)*NBOX+IB, m_cumEVENTi1, EVENTTYP(I),   DTIME(I),   &
                                                                                           m_curNWalki(IS),  m_curWTimei(IS)-DTIME(I),  &
                                                                                           m_curNFlighti(IS), m_curFTimei(IS)-DTIME(I)

                    m_curNWalki(IS)   = 0
                    m_curNFlighti(IS)  = 0
                    m_curWTimei(IS)   = 0.D0
                    m_curFTimei(IS)   = 0.D0
                    m_cumEVENTi = m_cumEVENTi + 1
                   end if
                 end if

                 write(m_j_v_IOUnit1, fmt="(1x, I8, 1x, I8, 16(1x, 1PE13.4,1x))") (JOB-1)*NBOX+IB, EVENTTYP(I), DTIME(I), COSDIR, &
                                               DX(I)/LATT, DY(I)/LATT, DZ(I)/LATT, (DX(I)*DX(I)+DY(I)*DY(I)+DZ(I)*DZ(I))/LATT/LATT
                 !--- update the moving direction, current site
                 m_preDIRi(IS,1) = DX(I)
                 m_preDIRi(IS,2) = DY(I)
                 m_preDIRi(IS,3) = DZ(I)

                 m_curSITEi(IS)    = NEWSITE(NEARSITE(I))
                 m_curPOSi(IS,1:3) = NEWPOS(NEARSITE(I), 1:3)
                 m_preTIMEi(IS)    = TIME
                 m_curOCCUai((IS-1)*mp_MXAI+1:IS*mp_MXAI) = hOCCUai((NEARSITE(I)-1)*mp_MXAI+1:NEARSITE(I)*mp_MXAI)

              end do
           end do !--- end loop for boxes

           deallocate(hSITE, hSitePOS, hSTAT, hOCCUa, hSITEi, hOCCUai)
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
       integer::I, J, hFile, event, nwalk, ibin, IB, NW, NF, IE
       character*128::tstr
       type(MDRecordStamp)::STAMP

              !--- calculate the histogram of event time
              rewind(m_j_v_IOUnit)
              rewind(m_j_v_IOUnit1)

              read(m_j_v_IOUnit1, *) tstr
              do I=1, m_cumEVENTi1
                 read(m_j_v_IOUnit1, *) IB, event, DTIME, COSDIR, DX, DY, DZ
              end do

              close(m_j_v_IOUnit)
              close(m_j_v_IOUnit1)

            do I=1, size(m_trackSimBox)
               write(tstr,*) I
               tstr = adjustl(tstr)
               write(*, fmt="(A)") "Jumping track to be saved in "//m_OUTPATH(1:len_trim(m_OUTPATH))//"Inter_Jump_Track"//tstr(1:len_trim(tstr))
               STAMP%ITime  = 0
               STAMP%ISect  = 1
               STAMP%ICfg   = 0
               STAMP%IRec   = 0
               STAMP%Time   = 0
               call Putout_Instance_Config_SimMDBox(m_OUTPATH(1:len_trim(m_OUTPATH))//"Inter_Jump_Track" &
                                //tstr(1:len_trim(tstr)), m_trackSimBox(I), STAMP)
            end do

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
          if(allocated(m_curSITEi))  deallocate(m_curSITEi)
          if(allocated(m_curPOSi))   deallocate(m_curPOSi)
          if(allocated(m_curOCCUai)) deallocate(m_curOCCUai)
          if(allocated(m_preTIMEi))  deallocate(m_preTIMEi)
          if(allocated(m_preDIRi))   deallocate(m_preDIRi)
          if(allocated(m_curNWalki)) deallocate(m_curNWalki)
          if(allocated(m_curNFlighti)) deallocate(m_curNFlighti)
          if(allocated(m_curWTimei)) deallocate(m_curWTimei)
          if(allocated(m_curFTimei)) deallocate(m_curFTimei)

          call Release_SimBoxArray(m_trackSimBox)

          m_NREC = 0
          return
  end subroutine MyCleaner
  !**********************************************************************************

 end module JumpRate_SIA_bcc_GPU


 !**********************************************************************************
 Program Diff_Interstial_GPU_Main
 use MD_SimBoxArray_ToolShell_14_GPU
 use JumpRate_SIA_bcc_GPU

 implicit none

       call APPSHELL_AddRecord( PRERECORD=MyInitialize, &
                                RECORDPROC=MyRECORD,    &
                                AFTRECORD= MyCleaner)

       call Main_ANALYSIS(0)

       stop
 End program Diff_Interstial_GPU_Main
