  module RefVoronoiVacancy_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to identify vacancy and interstitials. A initial configure is considered
  !                  as a reference. Two methods are used to identify the vacancies and interstitials.
  !                  Method #1:  Each lattice in the reference configure is a seed of Voronoi tessellation.
  !                              The program identify a Voronoi volume as a vacancy or interstitial by the
  !                              number atoms of the current configure falling in the Voronoi volume.
  !
  !                  Method #2:  Each lattice is considered as a reference point. When no atom falls in the
  !                              cutoff radiu of the lattice, the lattice is idenfied as a vacancy. An atom
  !                              is identified as an interstitial if it dose not fall in the cutoof rdius of
  !                              any lattice.
  !
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.cuf
  !                       MD_Globle_Variables_GPU.cuf
  !
  !                  SEE ALSO____________________________________________________________________________
  !                       RefCoordNum.F90
  !                       VoronoiTessellationM.F90
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2014-05 (Hou Qing)
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_MultiGPU_Basic
  !-----------------------------------------------
    implicit none

      integer::m_processid=0
      !$$--- the id of  filename get from SimMDCtrl for I/O.
      !$$    Refering the definition of f_others in SimMDCtrl.
      character(len=11),parameter, private::mp_FTAGI="&AUXF_RVVIN"
      character(len=12),parameter, private::mp_FTAGO="&AUXF_RVVOUT"
      character(len=256),          private::m_OUTFILE =""             ! filename of output data

      !--- C2050
      !integer,parameter::BLOCKSIZE = 512, BLOCKDIMX=32
      !--- C1060
      integer, parameter, private::mp_BLOCKSIZE = 256, mp_BLOCKDIMX=128

      !--- calculation control parameters
      integer,            private::m_INITED       = 0
      integer, parameter, private::mp_methVORONOI = 1
      integer, parameter, private::mp_methRADIUS  = 2
      integer,            private::m_Method = mp_methVORONOI    ! indcating the method used for identifying the vacancy
      real(KINDDF),       private::m_RADCUT = -1                !cutoff radius used only when m_Method = mp_methRADIUS

      integer,       private::m_LFOR    = 0                     !-- parameter indicating if force on a V-volume need to be calculated
      integer,       private::m_LPOT    = 0                     !-- parameter indicating if potential of a V-volume need to be calculated
      integer,       private::m_LKEN    = 0                     !-- parameter indicating if kinetice energy of a V-volume need to be calculated
      integer,       private::m_LSTRESS = 0                     !-- parameter indicating if stress on a V-volume need to be calculated

      real(KINDDF),  private::hm_BOND0 = 1.5D0                  !-- the bond-length in lattice unitused for neighbore calculation of reference system


      !$$--- the filename of reference configurations.
      integer, parameter, private::mp_MXCFG =4                  !-- the permitted number of reference configurations
      integer,            private::m_NUMCFG =0                  !-- the number of reference configurations
      character*256,      private::m_RefCfgFile(mp_MXCFG) = ""

      !--- arrayies on host for the refernce system
      type(SimMDBox)::hm_RefSimBox               !-- the reference box

      integer, private::hm_NC0                        !-- the number of cells in reference box
      integer, private::hm_MXNAC0                     !-- the max number of reference particles in a cell
      integer, private::hm_NCELL0(3)                  !-- the number of cells in 3 dimension for reference box
      integer, private::hm_PD0(3)                     !-- the periodic condition

      integer, dimension(:), allocatable, private::hm_NAC0
      integer, dimension(:), allocatable, private::hm_IA1th0
      integer, dimension(:), allocatable, private::hm_SSTAT        !-- the occupid state of lattice site,

      !--- the ID of  Voronoi site that atoms occupy
      integer, private::m_curNAPDEV = 0                              !the current number of atoms per device
      !---parameters used in calculating CN
      integer,  parameter,private::mp_NNC=27
      integer,  constant, private::dcm_NIX(mp_NNC)
      integer,  constant, private::dcm_NIY(mp_NNC)
      integer,  constant, private::dcm_NIZ(mp_NNC)

      !--- Workspace for Wigner-Seize defects
      private:: WSDefect_DEVWS
      type::WSDefect_DEVWS    
         !--- position of atoms in the (reference) first system
         !    NOTE: the positions on device have been sorted
          real(KINDDF), device, dimension(:,:), allocatable::XP0
          !--- the linked cell information
          integer,   device, dimension(:),      allocatable::NAC0
          integer,   device, dimension(:),      allocatable::IA1th0
          !--- the occupying statu of sites
          integer,   device, dimension(:),      allocatable::SITE
      end type WSDefect_DEVWS
      type(WSDefect_DEVWS), dimension(m_MXDEVICE), private::dm_WSDefect


      !--- some private routines
      private:: Allocate_WorkingArray0_template, &
                Clear_WorkingArray0_template,    &
                Allocate_WorkingArray0_DEV,      &

                Allocate_WorkingArray1_template, &
                Clear_WorkingArray1_template,    &
                Allocate_WorkingArray1_DEV,      &

                LoadControlParameters,           &
                PrintControlParameters,          &
                Cal_OccupiedVoronoiSite_KERNEL,  &
                StartOnDevice_Voronoi_template,  &
                Cal_OccupiedSphereSite_KERNEL,   &
                StartOnDevice_Sphere_template

  contains
 !****************************************************************************
  subroutine Get_filename_Output(fname)
  !***  PURPOSE:   to extract the filename of output
  !
  !     INPUT:
  !
  !     OUPUT:     fname,    the output filename
      implicit none
      !--- dummy variables
      character*(*)::fname
           fname = m_OUTFILE(1:min(len_trim(m_OUTFILE), len(fname)) )
           return
  end subroutine Get_filename_Output
 !****************************************************************************

 !****************************************************************************
  subroutine Get_filename_RefCfg(Fname, IRef)
  !***  PURPOSE:   to extract the filename of output
  !
  !     INPUT:     IRef,     optional, the Iref-th reference configure file
  !
  !     OUPUT:     Fname,    the output filename
      implicit none
      !--- dummy variables
      character*(*)    ::Fname
      integer,optional ::IRef
      !--- local variables
       integer::I

           I = 1
           if(present(IRef)) then
              I = IRef
           end if
           Fname = m_RefCfgFile(I)(1:min(len_trim(m_RefCfgFile(I)), len(Fname)) )
           return
  end subroutine Get_filename_RefCfg
 !****************************************************************************

  !****************************************************************************
  subroutine Get_InputTag(Tag)
  !***  PURPOSE:   to extract the input tag of output
  !
  !     INPUT:
      implicit none
      !--- dummy variables
      character*(*)::Tag
           Tag = mp_FTAGI
           return
  end subroutine Get_InputTag
 !****************************************************************************

 !****************************************************************************
  subroutine Allocate_WorkingArray0_template(IDEV, XPR, NCA, IA1th)
  !***  PURPOSE:   to allocate working space for a device
  !
  !     INPUT:     IDEV,     the index of device on which the memory to be allocated
  !
  !     OUPUT:     XPR,      the allocated device memory storing the position of atoms in reference system
  !                NCA,      the allocated device memory storing number of referece atoms in linked-cells on the device
  !                IA1th,    the allocated device memory storing ID  of staring referece atoms on the device
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:,:), allocatable::XPR
           integer,      device, dimension(:),   allocatable::NCA, IA1th

      !--- Local vairables
           integer::ERR, NPRT, CURDEV
           integer::NIX(mp_NNC)=(/0,-1,-1,-1, 0, 0, -1, 1,-1, 0, 1,-1, 0, 1, 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/)
           integer::NIY(mp_NNC)=(/0, 0,-1, 1, 1, 0,  0, 0,-1,-1,-1, 1, 1, 1, 0, 1,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/)
           integer::NIZ(mp_NNC)=(/0, 0, 0, 0, 0, 1,  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/)


              !$$--- allocate memory for neighbore list on devices
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)
              allocate(XPR(hm_RefSimBox%NPRT,3), NCA(hm_NC0), IA1th(hm_NC0), STAT=err)

              if(ERR) then
                 write(*,*) "MDPSCU ERROR: Fail to allocate memory in RefVoronoiVacancy module on device", IDEV
                 stop
              end if

              NCA   = hm_NAC0
              IA1th = hm_IA1th0
              XPR   = hm_RefSimBox%XP

              dcm_NIX = NIX
              dcm_NIY = NIY
              dcm_NIZ = NIZ
              ERR = cudaSetDevice(CURDEV)
              return
  end subroutine Allocate_WorkingArray0_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray0_DEV()
  !***  PURPOSE:  to allocate and initialize the device memories for reference configuration
  !
  !     INPUT:
  !     OUTPUT:
  !     NOTE:   intialization of device memories must be called after
  !             calling Initialize_Globle_Variables_DEV, with the number
  !             of particle on each device is avaliable
  !             and also one call to Cal_NeighBoreList_DEV
  !             with the information of cells of the template
  !             system is available.
  !
  use MD_Globle_Variables_GPU
  use MD_NeighborsList_GPU

  implicit none
  integer::I, ERR

           call GetCellInform(hm_NC0, hm_NCELL0, hm_MXNAC0)

           allocate(hm_SSTAT(hm_RefSimBox%NPRT), hm_NAC0(hm_NC0), hm_IA1th0(hm_NC0), STAT=ERR)
           if(err .gt. 0) then
              write(*,fmt="(A)") "MDPSCU Error: fail to allocated working space in ReVoronoiVA"
              stop
           end if
           !$$--- save the current statu of globale variable on devices
           !$$    NOTE: hm_XP are the sorted positions of reference reference atoms
           !$$
           hm_NAC0   = hm_NAC
           hm_IA1th0 = hm_IA1th
           hm_SSTAT  = -1

          !$$--- NOTE: the latice postion in reference box have been reordered
          !$$          after call Neighborelist calculation
           hm_RefSimBox%XP   = hm_XP
           hm_RefSimBox%ITYP = hm_ITYP

           !$$--- allocate the memory allocated on devices
           do I=1, m_NDEVICE
              call Allocate_WorkingArray0_template(m_DEVICES(I), dm_WSDefect(I)%XP0, dm_WSDefect(I)%NAC0, dm_WSDefect(I)%IA1th0)
           end do
            
           m_INITED = IOR(m_INITED,1)
           return
  end subroutine Allocate_WorkingArray0_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_WorkingArray0_template(IDEV, RXP, NAC,IA1th)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_WorkingArray0_template
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !     OUPUT:
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:,:), allocatable::RXP
           integer, device, dimension(:), allocatable::NAC, IA1th

      !--- Local vairables
           integer::CURDEV,ERR, NPRT

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)
              !$$--- deallocate memory for neighbore list on devices
              if(allocated(RXP)) then
                 deallocate(RXP, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(NAC)) then
                 deallocate(NAC, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(IA1th)) then
                 deallocate(IA1th, STAT=ERR)
                 if(ERR) goto 100
              end if

             ERR = cudaSetDevice(CURDEV)
             return

   100       ERR = cudaSetDevice(CURDEV)
             write(*,*) "MDPSCU Warning: fail to deallocate ref-coordinate memory on device", IDEV
             call ONWARNING(gm_OnWarning)
             return
  end subroutine Clear_WorkingArray0_template
  !**********************************************************************************

  !****************************************************************************
  subroutine Allocate_WorkingArray1_template(IDEV, NAPDEV, SITE)
  !***  PURPOSE:   to allocate working space for a device
  !
  !     INPUT:     IDEV,     the index of device on which the memory to be allocated
  !                NAPDEV,   the max number of atoms that a devce treats
  !
  !     OUPUT:     SITE,     the allocated device memory storing the site ID that atoms occupyed
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV, NAPDEV
           integer, device, dimension(:), allocatable::SITE

      !--- Local vairables
           integer::ERR, NPRT, CURDEV

              !$$--- allocate memory for neighbore list on devices
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)
              allocate(SITE(NAPDEV), STAT=err)

              if(ERR) then
                 write(*,*) "MDPSCU ERROR: Fail to allocate memory in RefVoronoiVacancy module on device", IDEV
                 stop
              end if

              SITE = 0
              ERR = cudaSetDevice(CURDEV)
              return
  end subroutine Allocate_WorkingArray1_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray1_DEV()
  !***  PURPOSE:  to allocate and initialize the device memories for SITES
  !
  !     INPUT:
  !     OUTPUT:
  !     NOTE:   intialization of this module must be called after
  !             calling Initialize_Globle_Variables_DEV, with the number
  !             of particle on each device is avaliable
  !             and also one call to Cal_NeighBoreList_DEV
  !             with the information of cells of the template
  !             system is available.
  !
  use MD_Globle_Variables_GPU, only:m_NAPDEV

  implicit none
  integer::I


           !$$--- allocate the memory allocated on device 1-6
            do I=1, m_NDEVICE
               call Allocate_WorkingArray1_template(m_DEVICES(I), m_NAPDEV, dm_WSDefect(I)%SITE )
            end do
            m_curNAPDEV = m_NAPDEV
            m_INITED = IOR(m_INITED,2)
            return;
  end subroutine Allocate_WorkingArray1_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_WorkingArray1_template(IDEV, SITE)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_WorkingArray1_template
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !     OUPUT:
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, device, dimension(:), allocatable::SITE

      !--- Local vairables
           integer::CURDEV,ERR, NPRT

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)

              if(allocated(SITE)) then
                 deallocate(SITE, STAT=ERR)
                 if(ERR) goto 100
              end if

             ERR = cudaSetDevice(CURDEV)
             return

   100       ERR = cudaSetDevice(CURDEV)
             write(*,*) "MDPSCU Warning: fail to deallocate ref-coordinate memory on device", IDEV
             call ONWARNING(gm_OnWarning)
             return
  end subroutine Clear_WorkingArray1_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray1_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calucaltion of Neighbore Table
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

      !--- Local variables
      integer::I
      !---------------------------------------

           !$$--- clear the memory allocated on device 1-6
           do I=1, m_NDEVICE
              call Clear_WorkingArray1_template(m_DEVICES(I), dm_WSDefect(I)%SITE)
           end do

           m_INITED = IAND(m_INITED, 1)
      RETURN
  end subroutine Clear_WorkingArray1_DEV
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calucaltion of Neighbore Table
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use MD_Globle_Variables_GPU, only:m_DEVICES,m_NDEVICE
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

      !--- Local variables
      integer::I
      !---------------------------------------

           call Clear_WorkingArray1_DEV(SimBox, CtrlParam)

           !$$--- clear the memory allocated on device 1-6
           do I=1, m_NDEVICE
              call Clear_WorkingArray0_template(m_DEVICES(I), dm_WSDefect(I)%XP0, dm_WSDefect(I)%NAC0, dm_WSDefect(I)%IA1th0)
           end do
           if(allocated(hm_NAC0))   deallocate(hm_NAC0)
           if(allocated(hm_IA1th0)) deallocate(hm_IA1th0)
           if(allocated(hm_SSTAT))  deallocate(hm_SSTAT)

           call Release_SimMDBox(hm_RefSimBox)

           m_curNAPDEV = 0
           m_INITED    = 0
           m_NUMCFG    = 0
      RETURN
  end subroutine Clear_WorkingArray_DEV
  !**********************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, SimBox, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    CtrlParam, the control parameters, menber of NEEDDAMP
  !                           could be changed.
  !
  use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables
  implicit none
  !--- dummy vaiables
  character*(*)  ::fname
  type(SimMDBox) ::SimBox
  type(SimMDCtrl)::CtrlParam

  !--- local variables
  integer::N, I, LINE, NEEDDAMP, DAMPSCHEME, STATU
  character*256::STR
  character*32::STRNUMB(10)
  character*32::KEYWORD
  type(InputStatements)::INPUTS

            !$$--- load input statements
            call LoadExInput_SimMDCtrl(mp_FTAGI, Fname, CtrlParam)
            call ExtractAnalyComm_SimMDCtrl(CtrlParam)

            !$$--- to extract control parameter specific for this module
            call GetExInput_SimMDCtrl(CtrlParam, mp_FTAGI, INPUTS, STATU)

            !$$--- start interpreting the statements
            m_RefCfgFile = ""
            do I=1, size(m_RefCfgFile)
               write(STRNUMB(1),*) (I-1)
               STRNUMB(1) = adjustl(STRNUMB(1))
               KEYWORD    = "&REFCFG"//STRNUMB(1)(1:len_trim(STRNUMB(1)))
               call Get_InputStatements(KEYWORD, INPUTS, STR, LINE)
               if(LINE .gt. 0) then
                  call Extract_Substr(STR,1,n,m_RefCfgFile(I))
               else
                 if(I .eq. 1) then
                    KEYWORD    = "&REFCFG"
                    call Get_InputStatements(KEYWORD, INPUTS, STR, LINE)
                    call Extract_Substr(STR,1,n,m_RefCfgFile(1))
                 end if
               end if
            end do

            m_NUMCFG     = 0
            do I=1, size(m_RefCfgFile)
               if(len_trim(m_RefCfgFile(I)) .gt. 0) then
                  m_NUMCFG = m_NUMCFG + 1
                  m_RefCfgFile(m_NUMCFG) = m_RefCfgFile(I)
               end if
            end do

            !$$*** To get the method that determines the connectivity of atoms
             call Get_InputStatements("&METHOD", INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Substr(STR,1,n,STRNUMB)
                m_Method = 0
                if(N.ge.1) then
                   if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "USEVT") then
                      m_Method = mp_methVORONOI
                   else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "USEST") then
                       m_Method = mp_methRADIUS
                       call Extract_Numb(STR,1,n,STRNUMB)
                       if(n .ge. 1) m_RADCUT = DRSTR(STRNUMB(1))
                       if(m_RADCUT .le. 0) then
                          write(*,fmt="(A,I4,A)")  ' MDPSCU Error: radius used to identify vacancy is not given '
                          write(*,fmt="(A,I4,A)")  '               radius should be larger than zero'
                          write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                          write(*,fmt="(A)")       '               Process to be stopped'
                          stop
                       end if
                   end if
                end if

                if(m_Method .le. 0) then
                   write(*,fmt="(A,I4,A)")  ' MDPSCU Error: method to identify vacancy is not given'
                   write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                   write(*,fmt="(A, BZI6)") '        Usage: &METHOD  "USEVT", to use Voronoi tessellation to identify vacancies'
                   write(*,fmt="(A, BZI6)") '        Or:    &METHOD  "USEST", radcut"'
                   write(*,fmt="(A, BZI6)") '               to use a sphere region of radius radcut to identify vacancies'
                   write(*,fmt="(A)")       '               Process to be stopped'
                   stop
                end if
             end if

            !$$*** To get the indicator if the force to be calcualted
             call Get_InputStatements("&CALFOR", INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,1, n, STRNUMB)
                m_LFOR = ISTR(STRNUMB(1))
              end if

             call Get_InputStatements("&CALPOT", INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,1, n, STRNUMB)
                m_LPOT = ISTR(STRNUMB(1))
              end if

             call Get_InputStatements("&CALKEN", INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,1, n, STRNUMB)
                m_LKEN = ISTR(STRNUMB(1))
             end if

             call Get_InputStatements("&CALSTRESS", INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,1, n, STRNUMB)
                m_LSTRESS = ISTR(STRNUMB(1))
             end if

             !$$*** To get the output path
             call Get_InputStatements(mp_FTAGO, INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                 call Extract_Substr(STR,1,n,m_OUTFILE)
             end if

            !$$--- check input consistent-----
             if(m_RADCUT .gt. 0) then
                m_RADCUT = m_RADCUT*SimBox%RR
             else
                m_RADCUT = -minval(CtrlParam%RU)
             end if

             if(m_NUMCFG .le. 0) then
                write(*,fmt="(A)")          " MDPSCU Error: reference configure for RefVoronoiVacancy calculation is not given."
                write(*,fmt="(A,A,A)")      "               add the keyword in control file: &REFCFGn fname "
                write(*,fmt="(' Process to be stopped')")
                stop
             end if

             if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for RefVoronoiVacancy calculation is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
             else
                call CreateDataFolder_Globle_Variables(m_OUTFILE)
             end if
             call Release_InputStatements(INPUTS)
             return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine PrintControlParameters(hFile, SimBox, CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
   implicit none
  !----   DUMMY Variables
   integer,         intent(in)::hFile
   type(SimMDBox),  intent(in)::SimBox
   type(SimMDCtrl), intent(in)::CtrlParam

  !----   Local variables
    integer::I

         !$$--- print out the control parameters
          write(hFile,*) "!************ RefVoronoiVacancy module to be used **********"
          write(hFile,*) "!    With the control parameters:"
          write(hFile,FMT="(' !    Ref.config. to be loaded from: ', A)") &
                                    m_RefCfgFile(1)(1:len_trim(m_RefCfgFile(1)))
          do I = 2, m_NUMCFG
          write(hFile,FMT="(' !                              and: ', A)") &
                                    m_RefCfgFile(I)(1:len_trim(m_RefCfgFile(I)))
          end do

          write(hFile,fmt="(A)")  " !     "

          write(hFile,fmt="(A)")        " !    Vacancy and interstial to be determined by method: "
          select case(m_METHOD)
                case(mp_methVORONOI)
                     write(hFile,fmt="(A)")     " !      Voronoi tesselleation (default) with: "
                    if(m_LFOR.gt.0 .and. m_METHOD.eq.mp_methVORONOI ) then
                       write(hFile,FMT="(' !    Forces on Voronoi cells to be calculated: YES')")
                    else
                       write(hFile,FMT="(' !    Forces on Voronoi cells to be calculated: NO')")
                    end if

                    if(m_LPOT.gt.0 .and. m_METHOD.eq.mp_methVORONOI) then
                       write(hFile,FMT="(' !    Potential on Voronoi cells to be calculated: YES')")
                    else
                       write(hFile,FMT="(' !    Potential on Voronoi cells to be calculated: NO')")
                    end if

                    if(m_LKEN.gt.0 .and. m_METHOD.eq.mp_methVORONOI) then
                       write(hFile,FMT="(' !    Kinetic energy on Voronoi cells to be calculated: YES')")
                    else
                       write(hFile,FMT="(' !    Kinetic energy on Voronoi cells to be calculated: NO')")
                    end if
                case(mp_methRADIUS)
                     write(hFile,fmt="(A)")     " !      Sphere regions centered at reference lattices with: "
                     if(m_RADCUT .gt. 0) then
                        write(hFile,fmt="(A, 1X, F7.3, A)") " !                                           cutoff radius: ", &
                                                                                                   m_RADCUT/SimBox%RR, " (LU)"
                     else
                        write(hFile,fmt="(A, 1X, F7.3, A")  " !                         radius is force cutoff distance:",  &
                                                                                                 -m_RADCUT/SimBox%RR, " (LU)"
                        m_RADCUT = dabs(m_RADCUT)
                     end if
          end select
          write(hFile,fmt="(A)")  " !     "


          if(CtrlParam%NEEDDAMP .gt. 0) then
             write(hFile,FMT="(' !    Quickdampping to be performed for ', I7, ' time steps before ReferenceVacancy')") CtrlParam%NEEDDAMP
          else
             write(hFile,FMT="(' !    Quickdampping to be performed before ReferenceVacancy calculation: NO')")
          end if
          write(hFile,FMT="(' !    ')")

         return
   end subroutine PrintControlParameters
  !**********************************************************************************

  !**********************************************************************************
   subroutine Initialize_ReferenceVacancy_DEV(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the module

   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   use MD_TYPEDEF_PrintList,  only:Add_PrintProcess

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   type(SimMDBox) ::tSimBox
   type(SimMDCtrl)::tCtrlParam
   integer::I, IFILE


        !$$--- to check device version is used
          call Check_DEVICES()

        !$$---  we create a temperaral simulation box and control parameters
         call Copy_SimMDCtrl(CtrlParam, tCtrlParam)
         if(m_INITED .gt. 0) then
            call Clear_WorkingArray_DEV(SimBox, tCtrlParam)
         end if

        !$$--- to findout the I/O unit
         IFILE = 0
         do I=1, SIZE(tCtrlParam%f_tag)
            if(tCtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                m_OUTFILE = tCtrlParam%f_others(I)
            end if
            if(tCtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
               IFILE = I
            end if
         end do
         !$$--- check the input files
         if(IFILE .le. 0) then
            write(*,fmt="(' MDPSCU Error: the control file is missed in RefVoronoiVacancy module')")
            write(*,*) "                add the keyword in SETUP file:", mp_FTAGI
            write(*,fmt="(' Process to be stopped')")
            stop
         end if

         write(*,fmt="(A)") " MDPSCU Message: Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
         if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") "Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
         call LoadControlParameters(tCtrlParam%f_others(IFILE), SimBox, CtrlParam)
         call Add_PrintProcess(PrintControlParameters)

         !$$--- load the reference configuration
         call Release_SimMDBox(hm_RefSimBox)
          do I=1, m_NUMCFG
             call CopyInformation_SimMDBox(SimBox, tSimBox)
             tSimBox%proAutoload = 0
             call ClearDataProKWD_SimMDBox(tSimBox)
             call AddDataProKWD_SimMDBox(tSimBox, "XYZCOL")
             write(*,FMT="(' MDPSCU Message: Load reference configuration from: ', A)")m_RefCfgFile(I)(1:len_trim(m_RefCfgFile(I)))
             call Load_Config_SimMDBox(m_RefCfgFile(I), tSimBox)
             tSimBox%NA    = 0
             tSimBox%NA(I) = tSimBox%NPRT
             tSimBox%ITYP  = I
             tSimBox%STATU = CP_STATU_ACTIVE
             call Merge_SimMDBox(tSimBox, hm_RefSimBox)
             call Release_SimMDBox(tSimBox)
          end do

         !$$--- RM is needed for neighbore calculations
          tCtrlParam%NB_RM = hm_BOND0*hm_RefSimBox%RR

         !$$--- genereate the reference box
          hm_PD0 = CtrlParam%IFPD

         !$$*** to initialize the GPU variable
         call Initialize_Globle_Variables_DEV(hm_RefSimBox, tCtrlParam)

         !$$*** to paitition the system
         call Initialize_NeighboreList_DEV(hm_RefSimBox, tCtrlParam)
         call Cal_NeighBoreList_DEV(hm_RefSimBox, tCtrlParam)
         call Allocate_WorkingArray0_DEV()

         return
   end subroutine Initialize_ReferenceVacancy_DEV
  !**********************************************************************************

  !**********************************************************************************
   attributes(global) subroutine Cal_OccupiedVoronoiSite_KERNEL(NPART0, XP0, NC0, NAC0, IA1th0, NCX0, NCY0, NCZ0,  &
                                                                LBX0, LBY0, LBZ0, BSX0, BSY0, BSZ0, PDX, PDY, PDZ,   &
                                                                FROMA, TOA, NPART, XP,                               &
                                                                SITE)
  !$$***PURPOSE:  to identify the Voronoi volume that the atoms in. If an atom is in ith
  !$$             Voronoi volume, the AC of this Voronoi is accumulates 1
  !$$
  !$$
  !$$   INPUT:    NPART0,   the total number of particles in the whole ref.box
  !$$             XP0,      the position of the particles
  !$$             NC0,      the number of cells in the reference box
  !$$             NAC0      the number of particles in the cells on the device
  !$$             IA1th0,   the index of the first particle in a cell of the reference box
  !$$             NCX0, NCY0, NCZ0, the cell number in X, Y, Z direction of the reference box
  !$$             BSX0, BSY0, BSZ0, the size of the box in X, Y, Z, dierection of the reference box
  !$$             PDX, PDY,PDZ,  the periodic condition in X, Y, Z, dierection
  !$$
  !$$             FROMA,    the start atom to be considered on this device
  !$$             TOA,      the end atom to be considered on this device
  !$$
  !$$             NPART,    the total number of particles in the whole realbox
  !$$             XP,       the position of the particles in the real box
  !$$
  !$$     OUTPUT: SITE,      the index of Voronoi site that atom in
  !
  implicit none
      !
      !--- DUMMY VARIABLES
      integer, value::NPART0,NC0,NCX0,NCY0,NCZ0, PDX, PDY, PDZ, FROMA, TOA
      real(KINDDF),value::LBX0, LBY0, LBZ0, BSX0, BSY0, BSZ0
      real(KINDDF), device::XP0(NPART0,3)
      integer, device::NAC0(NC0), IA1th0(NC0)

      integer, value::NPART
      real(KINDDF), device::XP(NPART,3)
      integer,device::SITE(*)

      !--- Local variables
         real(KINDSF), parameter::eps=0.0001
         real(KINDSF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, CX, CY, CZ, RRMI, RR
         integer::IB, NB, IT, IA, K, JA, IMI
         integer::IX0, IY0, IZ0,IC0
         integer::NACC, IAC, IACE, CID, IX,IY, IZ

        !$$--- start process
              !$$IB -- the index of block
               IB  =  (blockidx%y-1) * griddim%x +  blockidx%x

              !$$NB -- the size of a block
               NB = blockdim%x*blockdim%y

              !$$IT -- the thread ID
               IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              !$$IA-- the global ID of an atom treated in this thread
               IA = (IB-1)*NB + (IT-1) + FROMA

               if(IA .GT. TOA) return

             !$$-- get the position of the atom
              POSX = XP(IA, 1)
              POSY = XP(IA, 2)
              POSZ = XP(IA, 3)
              SITE(IA-FROMA+1) = 0

             !$$-- if not in the reference box, return
              if((POSX .LT. LBX0 .OR. POSX-LBX0 .GT. BSX0) .and. .not.PDX ) return;
              if((POSY .LT. LBY0 .OR. POSY-LBY0 .GT. BSY0) .and. .not.PDY ) return;
              if((POSZ .LT. LBZ0 .OR. POSZ-LBZ0 .GT. BSZ0) .and. .not.PDZ ) return;

             !$$-- determine which cell of the reference box that the IAth atom in
              IX0 = int( (POSX - LBX0)/BSX0*dble(NCX0) -eps)+1
              IY0 = int( (POSY - LBY0)/BSY0*dble(NCY0) -eps)+1
              IZ0 = int( (POSZ - LBZ0)/BSZ0*dble(NCZ0) -eps)+1

             !$$--- we start scanning the potential Voronoi site the the atom in
              RRMI = 1.E32
              IMI  = 0
              do K=1, mp_NNC
                 !$$-- get ID of the neighboring cell
                  IZ = IZ0 + dcm_NIZ(K)
                  IY = IY0 + dcm_NIY(K)
                  IX = IX0 + dcm_NIX(K)

                  CX = C_ZERO
                  If(PDX) Then
                     if( IX.GT.NCX0 )then
                         IX = C_IUN
                         CX = BSX0
                     else if (IX.LT.C_IUN) then
                         IX = NCX0
                         CX = -BSX0
                     end if
                  end if

                  CY = C_ZERO
                  if(PDY) then
                     if( IY.GT.NCY0 )then
                         IY = C_IUN
                         CY = BSY0
                     else if (IY.LT.C_IUN) then
                         IY = NCY0
                         CY = -BSY0
                     end if
                  end if

                  CZ = C_ZERO
                  if(PDZ) then
                     if( IZ.GT.NCZ0 )then
                         IZ = C_IUN
                         CZ = BSZ0
                     else if (IZ.LT.C_IUN) then
                         IZ = NCZ0
                         CZ = -BSZ0
                     end if
                  end if
                  if( IX .GT. NCX0 .OR. IX .LT. C_IUN) cycle
                  if( IY .GT. NCY0 .OR. IY .LT. C_IUN) cycle
                  if( IZ .GT. NCZ0 .OR. IZ .LT. C_IUN) cycle

                  CID = NCX0*NCY0*(IZ-1)+NCX0*(IY-1)+IX

                  !$$NACC -- the number of atoms in neighboring cell K
                  NACC = NAC0(CID)

                  !$$IAC-- the index of start atom in cell K
                  IAC = IA1th0(CID)

                  !$$IACE-- the index of end atom in cell K
                  IACE = IAC + NACC -1

                  do JA=IAC, IACE
                     SEPX = POSX - XP0(JA,1) - CX
                     SEPY = POSY - XP0(JA,2) - CY
                     SEPZ = POSZ - XP0(JA,3) - CZ
                     RR = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                     if(RR .LT. RRMI) then
                        RRMI = RR
                        IMI  = JA
                     end if
                  end do
               end do

            !$$--- Now we can get the site ID that the particle occupying
             SITE(IA-FROMA+1) = IMI

      return
  end subroutine Cal_OccupiedVoronoiSite_KERNEL
  !**********************************************************************************

  !**********************************************************************************
  subroutine StartOnDevice_Voronoi_template(IDEV, NAC0, IA1th0, XP0, XP, SITE)
  !***  PURPOSE:  template to calculate the occupied sites of atoms
  !
  !    INPUT(CPU): IDEV,      the index of the device
  !                NAC0,      the number of reference atoms in the cells
  !                IA1th0,    the index of the first reference atom in a cell
  !                XP0,       the position of the reference atoms
  !                XP,        the position that the atoms occupy
  !
   use MD_Globle_Variables_GPU
   implicit none
  !----   DUMMY Variables
          integer::IDEV
          integer, device, dimension(:)::NAC0,IA1th0
          real(KINDDF), device, dimension(:,:)::XP0, XP
          integer, device, dimension(:)::SITE

  !----   Local variables
         type(dim3) :: blocks
         type(dim3) :: threads
         integer::ERR, CURDEV, BX, BY, NB, STARTA, ENDA, NA
  !--- start
               !$$--- copy the cell informations into device memeory
               ERR = cudaGetDevice(CURDEV)
               ERR = cudaSetDevice(m_DEVICES(IDEV))

               !$$--- the first atom on the device
               STARTA = hm_IA1th(m_STARTCELL(IDEV))

               !$$--- the last atom on the device
               ENDA = hm_IA1th(m_ENDCELL(IDEV))+hm_NAC(m_ENDCELL(IDEV))-1

               !$$--- the number of atoms on the device
                NA = ENDA - STARTA + 1

               !$$--- to determine size of a block (the number of threads in a block)
                BX = mp_BLOCKSIZE
                BY = 1

               !$$-- to determine the dimension of blocks( the number of blocks in a grid)
                NB  = (NA-1)/(BX*BY)+1
                blocks  = dim3(NB, 1, 1)
                threads = dim3(BX, BY, 1)
                call Cal_OccupiedVoronoiSite_KERNEL<<<blocks,threads>>>(hm_RefSimBox%NPRT, XP0, hm_NC0, NAC0, IA1th0,                            &
                                                                        hm_NCELL0(1),  hm_NCELL0(2),   hm_NCELL0(3),                             &
                                                                        hm_RefSimBox%BOXLOW(1), hm_RefSimBox%BOXLOW(2), hm_RefSimBox%BOXLOW(3),  &
                                                                        hm_RefSimBox%ZL(1),hm_RefSimBox%ZL(2), hm_RefSimBox%ZL(3),               &
                                                                        hm_PD0(1), hm_PD0(2),   hm_PD0(3),                                       &
                                                                        STARTA, ENDA, dm_NPRT, XP, SITE)
               ERR = cudaSetDevice(CURDEV)
         return

  end subroutine StartOnDevice_Voronoi_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Cal_Occupied_Voronoi_DEV(SimBox, CtrlParam, hSITE)
  !***  PURPOSE:  to calculate the site ID that atoms occupy
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     hSITE, optional, the occupied Voronoi sites that the atoms fall in
  !

  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)                ::SimBox
          type(SimMDCtrl)               ::CtrlParam
          integer, dimension(:),optional::hSITE

  !----   Local variables
         integer::I
  !***

         if(m_curNAPDEV .ne. m_NAPDEV) then
             call Clear_WorkingArray1_DEV(SimBox, CtrlParam)
          end if

          if(iand(m_INITED,2) .eq. 0) then
             call  Allocate_WorkingArray1_DEV()
          end if

           !$$--- to start neighbor calculation on devices
             do I=1, m_NDEVICE
                call StartOnDevice_Voronoi_template(I, dm_WSDefect(I)%NAC0, dm_WSDefect(I)%IA1th0, dm_WSDefect(I)%XP0, &
                                                    dm_WorkSpace%XP(I)%Data, dm_WSDefect(I)%SITE)
             end do

           !$$--- copy the coordinate number out
           if(present(hSITE)) then
             call COPYOUT_Occupancy(hSITE)
            end if

        RETURN
  end subroutine Cal_Occupied_Voronoi_DEV
  !****************************************************************************************

  !**********************************************************************************
   attributes(global) subroutine Cal_OccupiedSphereSite_KERNEL(NPART0, XP0, NC0, NAC0, IA1th0, NCX0, NCY0, NCZ0,   &
                                                               LBX0, LBY0, LBZ0, BSX0, BSY0, BSZ0, PDX, PDY, PDZ,  &
                                                               FROMA, TOA, NPART, XP,                              &
                                                               RC, SITE)
  !$$***PURPOSE:  to identify the Shpere region that the atoms in.
  !$$
  !$$
  !$$   INPUT:    NPART0,   the total number of particles in the whole ref.box
  !$$             XP0,      the position of the particles
  !$$             NC0,      the number of cells in the reference box
  !$$             NAC0      the number of particles in the cells on the device
  !$$             IA1th0,   the index of the first particle in a cell of the reference box
  !$$             NCX0, NCY0, NCZ0, the cell number in X, Y, Z direction of the reference box
  !$$             BSX0, BSY0, BSZ0, the size of the box in X, Y, Z, dierection of the reference box
  !$$             PDX, PDY,PDZ,  the periodic condition in X, Y, Z, dierection
  !$$
  !$$             FROMA,    the start atom to be considered on this device
  !$$             TOA,      the end atom to be considered on this device
  !$$
  !$$             NPART,    the total number of particles in the whole realbox
  !$$             XP,       the position of the particles in the real box
  !$$
  !$$             RC,       the radius of the sphere
  !$$     OUTPUT: SITE,      the index of site that atom in
  !
  implicit none
      !
      !--- DUMMY VARIABLES
      integer, value::NPART0,NC0,NCX0,NCY0,NCZ0, PDX, PDY, PDZ, FROMA, TOA
      real(KINDDF),value::LBX0, LBY0, LBZ0, BSX0, BSY0, BSZ0
      real(KINDDF), value::RC
      real(KINDDF), device::XP0(NPART0,3)
      integer, device::NAC0(NC0), IA1th0(NC0)

      integer, value::NPART
      real(KINDDF), device::XP(NPART,3)
      integer,device::SITE(*)

      !--- Local variables
         real(KINDSF), parameter::eps=0.0001
         real(KINDSF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, CX, CY, CZ, RRC, RR
         integer::IB, NB, IT, IA, K, JA, IS
         integer::IX0, IY0, IZ0,IC0
         integer::NACC, IAC, IACE, CID, IX,IY, IZ

        !$$--- start process
              !$$IB -- the index of block
               IB  =  (blockidx%y-1) * griddim%x +  blockidx%x

              !$$NB -- the size of a block
               NB = blockdim%x*blockdim%y

              !$$IT -- the thread ID
               IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              !$$IA-- the global ID of an atom treated in this thread
               IA = (IB-1)*NB + (IT-1) + FROMA

               if(IA .GT. TOA) return

             !$$-- get the position of the atom
              POSX = XP(IA, 1)
              POSY = XP(IA, 2)
              POSZ = XP(IA, 3)
              SITE(IA-FROMA+1) = 0

             !$$-- if not in the reference box, return
              if((POSX .LT. LBX0 .OR. POSX-LBX0 .GT. BSX0) .and. .not.PDX ) return;
              if((POSY .LT. LBY0 .OR. POSY-LBY0 .GT. BSY0) .and. .not.PDY ) return;
              if((POSZ .LT. LBZ0 .OR. POSZ-LBZ0 .GT. BSZ0) .and. .not.PDZ ) return;

             !$$-- determine which cell of the reference box that the IAth atom in
              IX0 = int( (POSX - LBX0)/BSX0*dble(NCX0) -eps)+1
              IY0 = int( (POSY - LBY0)/BSY0*dble(NCY0) -eps)+1
              IZ0 = int( (POSZ - LBZ0)/BSZ0*dble(NCZ0) -eps)+1

             !$$--- we start scanning the potential Voronoi site the the atom in
              RRC = RC*RC
              IS  = 0
              do K=1, mp_NNC
                 !$$-- get ID of the neighboring cell
                  IZ = IZ0 + dcm_NIZ(K)
                  IY = IY0 + dcm_NIY(K)
                  IX = IX0 + dcm_NIX(K)

                  CX = C_ZERO
                  If(PDX) Then
                     if( IX.GT.NCX0 )then
                         IX = C_IUN
                         CX = BSX0
                     else if (IX.LT.C_IUN) then
                         IX = NCX0
                         CX = -BSX0
                     end if
                  end if

                  CY = C_ZERO
                  if(PDY) then
                     if( IY.GT.NCY0 )then
                         IY = C_IUN
                         CY = BSY0
                     else if (IY.LT.C_IUN) then
                         IY = NCY0
                         CY = -BSY0
                     end if
                  end if

                  CZ = C_ZERO
                  if(PDZ) then
                     if( IZ.GT.NCZ0 )then
                         IZ = C_IUN
                         CZ = BSZ0
                     else if (IZ.LT.C_IUN) then
                         IZ = NCZ0
                         CZ = -BSZ0
                     end if
                  end if
                  if( IX .GT. NCX0 .OR. IX .LT. C_IUN) cycle
                  if( IY .GT. NCY0 .OR. IY .LT. C_IUN) cycle
                  if( IZ .GT. NCZ0 .OR. IZ .LT. C_IUN) cycle

                  CID = NCX0*NCY0*(IZ-1)+NCX0*(IY-1)+IX

                  !$$NACC -- the number of atoms in neighboring cell K
                  NACC = NAC0(CID)

                  !$$IAC-- the index of start atom in cell K
                  IAC = IA1th0(CID)

                  !$$IACE-- the index of end atom in cell K
                  IACE = IAC + NACC -1

                  do JA=IAC, IACE
                     SEPX = POSX - XP0(JA,1) - CX
                     SEPY = POSY - XP0(JA,2) - CY
                     SEPZ = POSZ - XP0(JA,3) - CZ
                     RR = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                     if(RR .LT. RRC) then
                        IS  = JA
                        exit
                     end if
                  end do
                  if(IS .gt. 0) exit
               end do

            !$$--- Now we can get the site ID that the particle occupying
             SITE(IA-FROMA+1) = IS

      return
  end subroutine Cal_OccupiedSphereSite_KERNEL
  !**********************************************************************************

  !**********************************************************************************
  subroutine StartOnDevice_Sphere_template(IDEV, NAC0, IA1th0, XP0, XP, SITE)
  !***  PURPOSE:  template to calculate the occupied sites of atoms
  !
  !    INPUT(CPU): IDEV,      the index of the device
  !                NAC0,      the number of reference atoms in the cells
  !                IA1th0,    the index of the first reference atom in a cell
  !                XP0,       the position of the reference atoms
  !                XP,        the position that the atoms occupy
  !
   use MD_Globle_Variables_GPU
   implicit none
  !----   DUMMY Variables
          integer::IDEV
          integer, device, dimension(:)::NAC0,IA1th0
          real(KINDDF), device, dimension(:,:)::XP0, XP
          integer, device, dimension(:)::SITE

  !----   Local variables
         type(dim3) :: blocks
         type(dim3) :: threads
         integer::ERR, CURDEV, BX, BY, NB, STARTA, ENDA, NA
  !--- start
               !$$--- copy the cell informations into device memeory
               ERR = cudaGetDevice(CURDEV)
               ERR = cudaSetDevice(m_DEVICES(IDEV))

               !$$--- the first atom on the device
               STARTA = hm_IA1th(m_STARTCELL(IDEV))

               !$$--- the last atom on the device
               ENDA = hm_IA1th(m_ENDCELL(IDEV))+hm_NAC(m_ENDCELL(IDEV))-1

               !$$--- the number of atoms on the device
                NA = ENDA - STARTA + 1

               !$$--- to determine size of a block (the number of threads in a block)
                BX = mp_BLOCKSIZE
                BY = 1

               !$$-- to determine the dimension of blocks( the number of blocks in a grid)
                NB  = (NA-1)/(BX*BY)+1
                blocks  = dim3(NB, 1, 1)
                threads = dim3(BX, BY, 1)
                call Cal_OccupiedSphereSite_KERNEL<<<blocks,threads>>>(hm_RefSimBox%NPRT, XP0, hm_NC0, NAC0, IA1th0,                             &
                                                                       hm_NCELL0(1),  hm_NCELL0(2),   hm_NCELL0(3),                              &
                                                                       hm_RefSimBox%BOXLOW(1), hm_RefSimBox%BOXLOW(2),  hm_RefSimBox%BOXLOW(3),  &
                                                                       hm_RefSimBox%ZL(1),hm_RefSimBox%ZL(2), hm_RefSimBox%ZL(3),                &
                                                                       hm_PD0(1), hm_PD0(2), hm_PD0(3),                                          &
                                                                       STARTA, ENDA, dm_NPRT, XP, m_RADCUT, SITE)
               ERR = cudaSetDevice(CURDEV)
         return

  end subroutine StartOnDevice_Sphere_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Cal_Occupied_Sphere_DEV(SimBox, CtrlParam, hSITE)
  !***  PURPOSE:  to calculate the site ID that atoms occupy
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     SITE, optional, the occupied reference sites by atomscoordinate in
  !                      reference box
  !

  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)::SimBox
          type(SimMDCtrl)::CtrlParam
          integer, dimension(:),optional::hSITE

  !----   Local variables
         integer::I
  !***
          !$$---

          if(m_curNAPDEV .ne. m_NAPDEV) then
             call Clear_WorkingArray1_DEV(SimBox, CtrlParam)
          end if

          if(iand(m_INITED,2) .eq. 0) then
             call  Allocate_WorkingArray1_DEV()
          end if

           !$$--- to start neighbor calculation on devices
             do I=1, m_NDEVICE
                 call StartOnDevice_Sphere_template(I, dm_WSDefect(I)%NAC0, dm_WSDefect(I)%IA1th0, dm_WSDefect(I)%XP0, &
                                                       dm_WorkSpace%XP(I)%Data, dm_WSDefect(I)%SITE)
             end do
           !$$--- copy the coordinate number out
           if(present(hSITE)) then
             call COPYOUT_Occupancy(hSITE)
            end if

           
        return
  end subroutine Cal_Occupied_Sphere_DEV
  !***********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Occupancy(hSITE)
  !***  PURPOSE:  to copy the ocuupancy of atoms from device  to
  !               hots array
  !
  !     INPUT:
  !     OUTPUT     SITE, array on host the occupied Voronoi sitescoordinate in
  !                      reference box
  !

  use MD_Globle_Variables_GPU, only:COPY_OUT_SHIFT_template
  implicit none
  !----   DUMMY Variables
          integer, dimension(:)::hSITE

  !----   Local variables
         integer::I
  !***
           !$$--- copy the Occupancy out
             do I=1, m_NDEVICE
                call COPY_OUT_SHIFT_template(I, dm_WSDefect(I)%SITE, hSITE)
             end do

        RETURN
  end subroutine COPYOUT_Occupancy
  !****************************************************************************************

  !****************************************************************************************
  integer function Get_Number_of_Sites() result(N)
  !***  PURPOSE:  to get number of sites of referece lattice
  !
  !     INPUT:
  !     OUTPUT
  !

  implicit none
  !----   DUMMY Variables

  !----   Local variables
  !***
          N =  hm_RefSimBox%NPRT
        return
  end function Get_Number_of_Sites
  !****************************************************************************************

  !****************************************************************************************
  subroutine Get_Site_Positions(NS, siteID, sitePOS)
  !***  PURPOSE:  to get the positions of lattices indicated by ID
  !
  !     INPUT:     NS,      the number of sites whose positions to be extracted
  !                siteID,  the ID of the lattice
  !
  !     OUTPUT     sitePOS, the position of the sites
  !

  use MD_Globle_Variables_GPU, only:dm_NPRT
  implicit none
  !----   DUMMY Variables
          integer, intent(in)::NS
          integer,      dimension(:),intent(in)::siteID
          real(KINDDF), dimension(:,:)         ::sitePOS

  !----   Local variables
         integer::I
  !***
          !$$---
           do I=1, NS
              if(siteID(I) .gt. 0) then
                 sitePOS(I,1:3) = hm_RefSimBox%XP(siteID(I), 1:3)
              end if
           end do

        return
  end subroutine Get_Site_Positions
  !****************************************************************************************

  !**********************************************************************************
  subroutine Cal_Site_State_DEV(SimBox, CtrlParam, hSTAT, hOCCUa)
  !***  DESCRIPTION: to calculate the state of the sites
  !
  !    INPUT:  SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !    OUTPUT: hSTAT,     the number of atoms allocated on the sites
  !            hOCCUa,    the ID of the atoms allocated on the sites.
  !                       Note the ID is the index of atoms after sorting in neighbore-list calculations
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT
  implicit none
       !--- dummy variables
       type(SimMDBox)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam

       integer,dimension(:)::hSTAT, hOCCUa
       !--- local variables
        integer::I, J, IB, NBOX, IP, NPRT0, NPRT, SHIFT, MXAI
        integer, dimension(:),allocatable::hSITE

  !-------------
           NPRT0 = hm_RefSimBox%NPRT
           NBOX  = size(hSTAT)/NPRT0
           NPRT  = dm_NPRT/NBOX
           MXAI =  size(hOCCUa)/size(hSTAT)

          !--- calculate the occupation of lattice sites
          !    hSITE is the site indice that tha atoms occupy
          allocate(hSITE(dm_NPRT))
          call Cal_Occupied_Voronoi_DEV(SimBox, CtrlParam, hSITE)

          !--- to calculate occupation state
           hSTAT  = 0
           hOCCUa = 0
           IP     = 0
           do IB=1, NBOX
              SHIFT = (IB-1)*NPRT0
              do I=1, NPRT
                 IP = IP + 1
                 if(hSITE(IP)>0) then
                    hSTAT(hSITE(IP)+SHIFT) = hSTAT(hSITE(IP)+SHIFT) + 1
                    do J=1, MXAI
                       !--- hOCCUa: the atoms ID on site SHIFT + hSITE(IP)
                       if(hOCCUa((SHIFT + hSITE(IP)-1)*MXAI + J).eq. 0) then
                          hOCCUa((SHIFT + hSITE(IP)-1)*MXAI + J) = IP
                          exit
                       end if
                    end do
                    if(J.gt.MXAI) then
                      write(*,fmt="(A)")         " MDPSCU Warning: number of atoms occupying a site excced permitted value"
                      write(*,fmt="(A,I5,A,I8)") "                 the permerttied value is: ",MXAI
                      write(*,fmt="(A,I8)")      "                 the site is: ",hSITE(IP)
                      call ONWARNING(gm_OnWarning)
                    end if
                 end if
              end do
           end do
           deallocate(hSITE)
           return
   end subroutine Cal_Site_State_DEV
  !****************************************************************************************

  !**********************************************************************************
  subroutine Output_Merged_Config(Stamp, SimBox, CtrlParam, hSITE)
  !***  DESCRIPTION:  to output the configure with reference configure merged
  !                   the reference lattices are refered as type 1.
  !                   the type ID of atoms are added with 1
  !
  !    INPUT:  Stamp,     the recoding stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !            hSITE,     the site ID that the atoms occupy
  !
  !    OUTPUT:
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT, hm_GID, hm_XP
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       integer, dimension(:),        intent(in):: hSITE
       !--- local variables
       character*256        :: GFILE
       type(MDRecordStamp)  :: tStamp
       type(SimMDBox)       :: tSimBox(2), swapBox
       integer              :: IB, IP0, IP, IA, IS, SHIFT, NPRT0, IBOX
       !---
              call Copy_SimMDBox(hm_RefSimBox, tSimBox(1))
              SHIFT = 0
              IP    = 0
              NPRT0 = SimBox(1)%NPRT
              IBOX  = (Stamp%ITest-1)*CtrlParam%MultiBox
              do IB=1, size(SimBox)
                 call Copy_SimMDBox(SimBox(IB), tSimBox(2))
                 tSimBox(2)%DIS  = 0
                 tSimBox(2)%ITYP = tSimBox(2)%ITYP + 1
                 do IP0=1, NPRT0
                    IP = IP + 1
                    if( hSITE(IP) .gt. 0) then
                        IA = hm_GID(IP)- SHIFT
                        tSimBox(2)%DIS(IA,1:3) = hm_XP(IP,1:3) - hm_RefSimBox%XP(hSITE(IP),1:3)
                    end if
                 end do
                 SHIFT  = SHIFT + NPRT0
                 !--- to merge the box
                 call MergeBoxs_SimMDBox(tSimBox, swapBox)
                 call GetPath(m_OUTFILE, GFILE)
                 GFILE = GFILE(1:len_trim(GFILE))//"MergedCfg"

                 IBOX         = IBOX  + 1
                 call Copy_RecordStamp(Stamp, tStamp)
                 tStamp%ITest = IBOX
                 tStamp%IBox  = 1
                 call STRCATI(GFILE, GFILE, "_Box", IBOX, 4)
                 call Putout_Instance_Config_SimMDBox(GFILE, swapBox, tStamp)
              end do

              call Release_SimMDBox(swapBox)
              call Release_SimMDBox(tSimBox(1))
              call Release_SimMDBox(tSimBox(2))
       return
  end subroutine Output_Merged_Config
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_Occupied_Voronoi(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calcualte and putout the occupation state of lattice and interstials
  !                  identified using method#1
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT,SynchronizeDevices, hm_GID, hm_ITYP, hm_XP, hm_EPOT, hm_FP, hm_EKIN
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE
       integer::hFile, I, J, warning, missed, IG, NG, OCC, NCOL, JOB,ITIME, ICFG
       real(KINDDF)::TIME
       integer, parameter::HSIZE = 20
       integer, dimension(0:HSIZE)::HIS
       integer, dimension(:),allocatable::hSITE
       integer, dimension(:,:),allocatable::hSTATE
       real(KINDDF), dimension(:),allocatable::hEPOT, hEKIN, hKSTRESS, hVSTRESS
       real(KINDDF), dimension(:,:),allocatable::hFP
       real(KINDDF), dimension(:),allocatable::hVAR
       character*32, dimension(:),allocatable::hTitle
       character*32::STRN
       !------

           JOB   = Stamp%ITest
           ITIME = Stamp%ITime
           ICFG  = Stamp%IRec(1)
           TIME  = Stamp%Time

           NG = SimBox(1)%NGROUP
           allocate(hSITE(dm_NPRT),hSTATE(hm_RefSimBox%NPRT,NG))
           if(m_LPOT .gt. 0)    allocate(hEPOT(hm_RefSimBox%NPRT))
           if(m_LKEN .gt. 0)    allocate(hEKIN(hm_RefSimBox%NPRT))
           if(m_LSTRESS .gt. 0) allocate(hKSTRESS(hm_RefSimBox%NPRT),hVSTRESS(hm_RefSimBox%NPRT))
           if(m_LFOR .gt. 0)    allocate(hFP(hm_RefSimBox%NPRT,3))

           if(dm_NPRT .gt. SimBox(1)%NPRT) then
              write(*,fmt="(A)")  ' MDPSCU Error: this program can handle only single box'
              write(*,fmt="(A)")  '               please seperate multi box into single box'
              write(*,fmt="(A)")  '               using tool MultBox2Boxes'
              stop
           end if

           hSITE = 0
           call Cal_Occupied_Voronoi_DEV(SimBox(1), CtrlParam, hSITE)
           call SynchronizeDevices()

           !$$--- covert hSITE to the number of atom in each Voronoi volume
           if(ITIME .LE. 0 .and. CtrlParam%INDEPBOX  .gt. 0) then
              hm_SSTAT = -1
              do I=1, dm_NPRT
                 if(hSITE(I)>0) then
                    hm_SSTAT(hSITE(I)) = 0
                 end if
              end do
           end if

           hSTATE  = 0
           missed  = 0
           warning = 0
           do I=1, dm_NPRT
              IG = hm_ITYP(I)
              if(hSITE(I)>0) then
                 hSTATE(hSITE(I), IG) = hSTATE(hSITE(I),IG) + 1
              end if

              if((hSITE(I) <= 0 .or. hSITE(I)>hm_RefSimBox%NPRT) ) then
                  missed = missed + 1
                  if(warning .le. 0) then
                      write(*,fmt="( ' MDPSCU Warning: an atom out of box, found in RefVoronoiVacancy')")
                      write(*,fmt="( '                 the atom ID0 and current ID are',I7, 1x, I7)") hm_GID(I), I
                      write(*,fmt="( '                 with postion (LU)',3(1PE14.6,1x))") hm_XP(I,1:3)/hm_RefSimBox%RR
                      write(*,fmt="( '                 ref.box size     ',3(1PE14.6,1x))") hm_RefSimBox%ZL/hm_RefSimBox%RR
                      write(*,fmt="( '                 ref.box low      ',3(1PE14.6,1x))") hm_RefSimBox%BOXLOW/hm_RefSimBox%RR
                      call ONWARNING(gm_OnWarning)
                      warning = 1
                  end if
              end if
           end do

           !$$--- to calculate the occupancy histogram
           HIS = 0
           do I=1, hm_RefSimBox%NPRT
              OCC = sum(hSTATE(I,1:NG))
              HIS(OCC) = HIS(OCC) + 1
              if(OCC > 0) then
                 hm_SSTAT(I) = OCC
              else
                if(hm_SSTAT(I) .ge. 0) hm_SSTAT(I) = 0
              end if
           end do

           !$$--- to prepair output colume
            NCOL = 3+2+NG
            if(m_LPOT .gt. 0)    NCOL = NCOL + 1
            if(m_LKEN .gt. 0)    NCOL = NCOL + 1
            if(m_LSTRESS .gt. 0) NCOL = NCOL + 2
            if(m_LFOR .gt. 0)    NCOL = NCOL + 3
            allocate(hVAR(NCOL),hTITLE(NCOL))
            hTitle(1) = "X(latt.)"
            hTitle(2) = "Y(latt.)"
            hTitle(3) = "Z(latt.)"
            hTitle(4) = "OCCUPATION+1"
            hTitle(5) = "OSS"
            NCOL = 5
            do I=1, NG
               NCOL = NCOL + 1
               write(STRN,*) I
               STRN = adjustl(STRN)
               hTitle(NCOL) = "by-type"//STRN(1:len_trim(STRN))
            end do

            if(m_LPOT .gt. 0) then
               NCOL = NCOL + 1
               hTitle(NCOL) = "EPOT(eV)"
               hEPOT = 0.D0
               do I=1, dm_NPRT
                  if(hSITE(I)>0) hEPOT(hSITE(I)) = hEPOT(hSITE(I)) + hm_EPOT(I)
               end do
            end if

            if(m_LKEN .gt. 0) then
               NCOL = NCOL + 1
               hTitle(NCOL) = "EKIN(eV)"
               hEKIN = 0.D0
               do I=1, dm_NPRT
                  if(hSITE(I)>0) hEKIN(hSITE(I)) = hEKIN(hSITE(I)) +hm_EKIN(I)
               end do
            end if

            if(m_LSTRESS .gt. 0) then
               NCOL = NCOL + 1
               hTitle(NCOL) = "KPRESS(kbar)"
               NCOL = NCOL + 1
               hTitle(NCOL) = "VPRESS(kbar)"
            end if

            if(m_LFOR .gt. 0) then
               NCOL = NCOL + 1
               hTitle(NCOL) = "FPx(dyn)"
               NCOL = NCOL + 1
               hTitle(NCOL) = "FPy(dyn)"
               NCOL = NCOL + 1
               hTitle(NCOL) = "FPz(dyn)"
               hFP = 0.D0
               do I=1, dm_NPRT
                  if(hSITE(I)>0) hFP(hSITE(I),1:3) = hFP(hSITE(I),1:3) + hm_FP(I,1:3)
               end do
            end if

           !$$--- output occupation distribution
           GFILE = ""
           call STRCATI(GFILE, m_OUTFILE, "P", m_processid, 4)
           call STRCATI(GFILE, GFILE, "_", JOB, 4)
           call STRCATI(GFILE, GFILE, ".", ICFG, 4)

           call AvailableIOUnit(hFile)
           write(*,*) "Save atomic Voronoi occupation  to "//GFILE(1:len_trim(GFILE))

           open(UNIT=hFile, file = GFILE, status='unknown')
            write(hFile, fmt="(A)") "!--- REFVORONOIVA RESULTS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(hFile, fmt="('!--- The number of reference lattices: ', I8)") hm_RefSimBox%NPRT
            write(hFile, fmt="('!--- The number of atoms             : ', I8)") dm_NPRT
            write(hFile, fmt="('!--- at time steps: ', I7,  ', time(ps): ', 1pE12.4)") ITIME, TIME
            write(hFile, fmt="('!--- Histogram for occupation:')")
            do I=0, HSIZE
               write(hFile, fmt="('!--- occupation= ', I3,  ', counts', I8)") I, HIS(I)
            end do
            if(missed .gt. 0) then
               write(hFile, fmt="('!--- there are ',  I8, 'atoms missed due to out of box')") missed
            end if
            write(hFile, fmt="('!--- Total occupation: ',  I8, ' in ', I8, ' sites')") count(hm_SSTAT.ge.0), size(hm_SSTAT)
            write(hFile,fmt="(A)")   "!"
            write(hFile,fmt="(A)")   "!---  SYMBOL EXPLAINATION:"
            write(hFile,fmt="(A)")   "!     OCP1:  One plus the number of atoms occupying a site. If OCP1=1, the site is a vacancy."
            write(hFile,fmt="(A)")   "!     OSS:   The statu of a site. OSS <0  indicates the site has nerver been visited by an atom."
            write(hFile,fmt="(A)")   "!     OTYP:  The number of atoms of TYPE occupying the site."
            write(hFile,fmt="(A)")   "!"

           !$$--- write out the XYZ format
            write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            write(hFile, fmt="(A,1X,I8)")            "&NATOM      ", hm_RefSimBox%NPRT
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE    ", hm_RefSimBox%ZL/hm_RefSimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW     ", hm_RefSimBox%BOXLOW/hm_RefSimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT       ", SimBox(1)%RR*CP_CM2A
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL     ", 1, 2, 3
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&OCP1COL    ", 4
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&OSSCOL     ", 5
            do I=1, NG
               write(STRN,*) I
               STRN = adjustl(STRN)
               STRN = "&OTYP"//STRN(1:len_trim(STRN))//"COL"
               write(hFile, fmt="(A12,1X,3(I4,1X))") STRN, I+5
            end do

            I = 5+NG
            if(m_LPOT .gt. 0) then
               I = I + 1
               write(hFile, fmt="(A,1X,3(I4,1X))")   "&EPOTCOL    ",  I
            end if

            if(m_LKEN .gt. 0) then
               I = I + 1
               write(hFile, fmt="(A,1X,3(I4,1X))")   "&EKINCOL    ",  I
            end if

            if(m_LSTRESS .gt. 0) then
               I = I + 1
               write(hFile, fmt="(A,1X,3(I4,1X))")   "&KPRESS     ",  I
               I = I + 1
               write(hFile, fmt="(A,1X,3(I4,1X))")   "&VPRESS     ",  I
            end if

            if(m_LFOR .gt. 0) then
               I= I + 1
               write(hFile, fmt="(A,1X,3(I4,1X))")   "&FPCOL      ",  I, I+1, I+2
               I = I + 2
            end if
            write(hFile, *)

           !$$--- write the data table header
            write(hFile, fmt="('!--- Distribution of defect: ',  I8, ' in ', I8, ' sites')") count(hm_SSTAT.ge.0), size(hm_SSTAT)
            write(STRN,*) NCOL
            STRN = adjustl(STRN)
            GFILE = "( '!---' "//STRN(1:len_trim(STRN))//"(A13, 1X))"
            write(hFile, fmt=GFILE)  (hTitle(I)(1:len_trim(hTitle(I))), I=1,  NCOL)

            write(STRN,*) NG+2
            STRN = adjustl(STRN)
            GFILE = "(5X,3(1PE13.4,1x), 6x, "//STRN(1:len_trim(STRN))//" (I6,6x)"
            write(STRN,*) NCOL-NG+2
            STRN = adjustl(STRN)
            GFILE = GFILE(1:len_trim(GFILE))//","//STRN(1:len_trim(STRN))//"(1x, 1PE13.4)"//")"

            !$$--- write out data
            do I=1, hm_RefSimBox%NPRT
               J = 0
               if(m_LPOT .gt. 0) then
                  J = J + 1
                  hVAR(J) = -hEPOT(I)*CP_ERGEV
               end if

               if(m_LKEN .gt. 0) then
                  J = J + 1
                  hVAR(J) = hEKIN(I)*CP_ERGEV
               end if

               if(m_LSTRESS .gt. 0) then
                  J = J + 1
                  hVAR(J) = hKSTRESS(I)
                  J = J + 1
                  hVAR(J) = hVSTRESS(I)
               end if

               if(m_LFOR .gt. 0) then
                  J = J + 1
                  hVAR(J) = hFP(I,1)
                  J = J + 1
                  hVAR(J) = hFP(I,2)
                  J = J + 1
                  hVAR(J) = hFP(I,3)
               end if

               write(hFile,fmt=GFILE)hm_RefSimBox%XP(I,1:3)/hm_RefSimBox%RR, sum(hSTATE(I,1:NG))+1, hm_SSTAT(I), hSTATE(I,1:NG), hVAR(1:J)
            end do
           close(hFile)
           !$$--- putout merged configre
            call Output_Merged_Config(Stamp, SimBox, CtrlParam, hSITE)

           !$$--- relase allocate memory
           deallocate(hSITE,hSTATE,hVAR,hTITLE)
           if(allocated(hEPOT))    deallocate(hEPOT)
           if(allocated(hEKIN))    deallocate(hEKIN)
           if(allocated(hKSTRESS)) deallocate(hKSTRESS)
           if(allocated(hVSTRESS)) deallocate(hVSTRESS)
           if(allocated(hFP))      deallocate(hFP)
          return
  end subroutine RECORD_Occupied_Voronoi
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_Occupied_Sphere(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calcualte and putout the occupation state of lattice and interstials
  !                  identified using method#2
  !
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT,SynchronizeDevices, hm_ITYP, hm_XP
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE
       integer::hFile, I, J, warning, IG, NG, OCC, NCOL, JOB,ITIME,ICFG
       real(KINDDF)::TIME
       integer, dimension(:),allocatable::hSITE
       integer, dimension(:,:),allocatable::hSTATE
       integer, dimension(:),allocatable::hNI
       character*32, dimension(:),allocatable::hTitle
       character*32::STRN
       !----

           JOB   = Stamp%ITest
           ITIME = Stamp%ITime
           ICFG  = Stamp%IRec(1)
           TIME  = Stamp%Time

           NG = SimBox(1)%NGROUP
           allocate(hSITE(dm_NPRT),hSTATE(hm_RefSimBox%NPRT,NG), hNI(NG))

           if(dm_NPRT .gt. SimBox(1)%NPRT) then
              write(*,fmt="(A)")  ' MDPSCU Error: this program can handle only single box'
              write(*,fmt="(A)")  '               please seperate multi box into single box'
              write(*,fmt="(A)")  '               using tool MultBox2Boxes'
              stop
           end if

           hSITE = 0
           call Cal_Occupied_Sphere_DEV(SimBox(1), CtrlParam, hSITE)
           call SynchronizeDevices()
           !$$--- covert hSITE to the number of atom in each Sphere volume
           hm_SSTAT = 0
           if(ITIME .LE. 0 .and. CtrlParam%INDEPBOX  .gt. 0) then
              hm_SSTAT = -1
              do I=1, dm_NPRT
                 if(hSITE(I)>0) then
                    hm_SSTAT(hSITE(I)) = 0
                 end if
              end do
           end if

           hSTATE = 0
           hNI  = 0
           do I=1, dm_NPRT
              IG = hm_ITYP(I)
              if(hSITE(I)>0) then
                 hSTATE(hSITE(I), IG) = hSTATE(hSITE(I),IG) + 1
              else
                 hNI(IG) = hNI(IG) + 1
              end if
           end do

           do I=1, hm_RefSimBox%NPRT
              OCC = sum(hSTATE(I,1:NG))
              if(OCC .gt. 0) then
                 hm_SSTAT(I) = OCC
              else
                if(hm_SSTAT(I) .ge. 0) hm_SSTAT(I) = 0
              end if
           end do

           !$$--- to prepair output colume
            NCOL = 3+2+NG
            allocate(hTITLE(NCOL))
            hTitle(1) = "X(latt.)"
            hTitle(2) = "Y(latt.)"
            hTitle(3) = "Z(latt.)"
            hTitle(4) = "STAT"
            hTitle(5) = "OSS"
            NCOL = 5
            do I=1, NG
               NCOL = NCOL + 1
               write(STRN,*) I
               STRN = adjustl(STRN)
               hTitle(NCOL) = "by-type"//STRN(1:len_trim(STRN))
            end do

           !$$--- output occupation distribution

           GFILE = ""
           call STRCATI(GFILE, m_OUTFILE, "P", m_processid, 4)
           call STRCATI(GFILE, GFILE, "_", JOB, 4)
           call STRCATI(GFILE, GFILE, ".", ICFG, 4)

           call AvailableIOUnit(hFile)
           write(*,*) "Save atomic Sphere occupation  to "//GFILE(1:len_trim(GFILE))

           open(UNIT=hFile, file = GFILE, status='unknown')
            write(hFile, fmt="(A)") "!--- REFVORONOIVA RESULTS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(hFile, fmt="('!--- The number of reference lattices: ', I8)") hm_RefSimBox%NPRT
            write(hFile, fmt="('!--- The number of atoms                       : ', I8)") dm_NPRT
            write(hFile, fmt="('!---      at time steps: ', I7,  ', time(ps): ', 1pE12.4)") ITIME, TIME
            write(hFile, fmt="('!--- The number of vacancies                      : ', I8)")  count(hm_SSTAT .eq. 0)
            write(hFile, fmt="('!--- The number of never visited sites: ', I8)")  count(hm_SSTAT .lt. 0)
            write(hFile, fmt="('!--- The number of sites occcupied by 1 atoms: ', I8)")  count(hm_SSTAT .eq. 1)
            write(hFile, fmt="('!--- The number of sites occcupied by 2 atoms: ', I8)")  count(hm_SSTAT .eq. 2)
            write(hFile, fmt="('!--- The number of sites occcupied by 3 atoms: ', I8)")  count(hm_SSTAT .eq. 3)
            write(hFile, fmt="('!--- The number of sites occcupied by 4 atoms: ', I8)")  count(hm_SSTAT .eq. 4)
            write(hFile, fmt="('!--- The number of sites occcupied by 5 atoms: ', I8)")  count(hm_SSTAT .eq. 5)
            write(hFile, fmt="('!--- The number interstitials        :', I8)")  sum(hNI(1:NG))

            do I=1, NG
            write(hFile, fmt="('!---     number in atomic type       :', I3, ' is ', I8)") I, hNI(I)
            end do
            write(hFile, fmt="('!--- Total occupation: ',  I8, ' in ', I8, ' sites')") count(hm_SSTAT.ge.0), size(hm_SSTAT)
            write(hFile,fmt="(A)")   "!"
            write(hFile,fmt="(A)")   "!---  SYMBOL EXPLAINATION:"
            write(hFile,fmt="(A)")   "!     STAT:   1, for vacancy; 2 for occuoied site; 3 for interstitial"
            write(hFile,fmt="(A)")   "!     OSS:    The statu of a site. OSS <0  indicates the site has nerver been visited by an atom."
            write(hFile,fmt="(A)")   "!     OTYP:   The number of atoms of TYPE occupying the site."
            write(hFile,fmt="(A)")   "!     &NATOM: The number of sites + interstitial"
            write(hFile,fmt="(A)")   "!     &XYZCOL:The position ofd sites if STAT=1 or 2, position of interstitial if STAT =3"
            write(hFile,fmt="(A)")   "!"

           !$$--- write out the XYZ format
            write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            write(hFile, fmt="(A,1X,I8, A)")         "&NATOM      ", hm_RefSimBox%NPRT + sum(hNI(1:NG))
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE    ", hm_RefSimBox%ZL/hm_RefSimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW     ", hm_RefSimBox%BOXLOW/hm_RefSimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT       ", SimBox(1)%RR*CP_CM2A
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL     ", 1, 2, 3
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&STATCOL    ", 4
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&OSSCOL     ", 5
            do I=1, NG
               write(STRN,*) I
               STRN = adjustl(STRN)
               STRN = "&OTYP"//STRN(1:len_trim(STRN))//"COL"
               write(hFile, fmt="(A12,1X,3(I4,1X))") STRN, I+5
            end do

           !$$--- write the data table header
            write(hFile, fmt="('!--- Distribution of defects: ')")
            write(STRN,*) NCOL
            STRN = adjustl(STRN)
            GFILE = "( '!---' "//STRN(1:len_trim(STRN))//"(A13, 1X))"
            write(hFile, fmt=GFILE)  (hTitle(I)(1:len_trim(hTitle(I))), I=1,  NCOL)

            write(STRN,*) NG+2
            STRN = adjustl(STRN)
            GFILE = "(5X,3(1PE13.4,1x), 6x, "//STRN(1:len_trim(STRN))//" (I6,6x))"

            !$$--- write out vacancy
            do I=1, hm_RefSimBox%NPRT
               if(sum(hSTATE(I,1:NG)) .le. 0)then
                  write(hFile,fmt=GFILE) hm_RefSimBox%XP(I,1:3)/hm_RefSimBox%RR, 1, hm_SSTAT(I), hSTATE(I,1:NG)
               end if
            end do
            !$$--- write out occupyied lattice
            do I=1, hm_RefSimBox%NPRT
               if(sum(hSTATE(I,1:NG)) .gt. 0)then
                  write(hFile,fmt=GFILE) hm_RefSimBox%XP(I,1:3)/hm_RefSimBox%RR, 2, hm_SSTAT(I), hSTATE(I,1:NG)
               end if
            end do
            !$$--- write out the interstials
            do I=1, dm_NPRT
               IG = hm_ITYP(I)
               if(hSITE(I).le.0) then
                  write(hFile,fmt=GFILE) hm_XP(I,1:3)/SimBox(1)%RR, 3, 0, (0, J=1,IG-1), IG, (0, J=IG+1, NG)
               end if
            end do
           close(hFile)
           !$$--- putout merged configre
            call Output_Merged_Config(Stamp, SimBox, CtrlParam, hSITE)

           !$$--- relase allocate memory
           deallocate(hSITE,hSTATE,hNI,hTITLE)
          return
  end subroutine RECORD_Occupied_Sphere
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_ReferenceVacancy(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to provide a interface to MD_SimBoxArray_AppShell_14_GPU
  !                  or MD_SimBoxArray_ToolShell_14_GPU through the call to
  !                  RECORD_ReferenceVacancy_TOOL
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
   implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables

          select case(m_Method)
             case (mp_methVORONOI)
                   call RECORD_Occupied_Voronoi(Stamp, SimBox, CtrlParam)

             case (mp_methRADIUS)
                   if(m_RADCUT .le. 0) then
                      write(*,fmt="(' MDPSCU Error: the cutoff radius used to vacancy is smaller than zeoro in RefVoronoiVacancy module')")
                      write(*,fmt="(' Process to be stopped')")
                      stop
                   end if
                   call RECORD_Occupied_Sphere(Stamp, SimBox, CtrlParam)

             case default
                  write(*,fmt="(' MDPSCU Warning: no method for identifying vacancy is missed in RefVoronoiVacancy module')")
                  write(*,fmt="('                 no vacancy identify to be performed')")
                  call ONWARNING(gm_OnWarning)
                  return
          end select

          return
  end subroutine RECORD_ReferenceVacancy
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_ReferenceVacancy_TOOL(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to provide a interface to MD_SimBoxArray_ToolShell_14_GPU
  !                  It is assumed the the neighbor list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
   implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables

          if(Stamp%ITime .LT. 0) return
          call RECORD_ReferenceVacancy(Stamp, SimBox, CtrlParam)

          return
  end subroutine RECORD_ReferenceVacancy_TOOL
  !****************************************************************************************




  end module RefVoronoiVacancy_GPU
