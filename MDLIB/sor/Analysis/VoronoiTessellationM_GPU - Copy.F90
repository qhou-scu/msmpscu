  module VoronoiTessellation_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  To do the Voronoi tessellation of atoms base on the algorithm
  !                  of Medvedev (J. Comput. Phys. 67(1986)223-229)
  !
  !                  The difference between this code and VoronoiTessellationT.cuf  is
  !                  that Meddvedev algorithm is used, while in VoronoiTessellationT the
  !                  algorithm of Tanemura et al is used
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.cuf  (.F90)
  !                       MD_Globle_Variables_GPU.cuf (.F90)
  !
  !                  REFERENCE____________________________________________________________________________
  !                       N. N. Medvedev, J. Comput. Phys. 67(1986)223-229
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2013-08 (Hou Qing)
  !
  !________________________________________________________________________________________________________



  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_MultiGPU_Basic
  !-----------------------------------------------
    implicit none

      integer::m_processid=0
      !--- the filenamed get from SimMDCtrl for I/O.
      !    Refering the definition of f_others in SimMDCtrl.
      character(len=9),parameter, private::mp_FTAGI="&AUXF_VTI"
      character(len=9),parameter, private::mp_FTAGO="&AUXF_VTO"
      character(len=256),  private::m_OUTFILE =""             ! filename of output data

      !--- control parameters determine what will be calculated
      !    and ouput
      integer, private::m_OUTPUTDT  = 0                      ! if >0, DT will be ouput
      integer, private::m_OUTPUTVT  = 0                      ! if >0, VT will be ouput
      integer, private::m_OUTPUTST  = 0                      ! if >0, VT signature will be ouput


      !--- C2050
      integer,parameter, private::mp_BLOCKSIZE = 256, mp_BLOCKDIMX=128
      !--- C1060
      !integer,parameter, private::mp_BLOCKSIZE = 128, mp_BLOCKDIMX=128

      !--- calculation control parameters
      integer, private::m_INITED = 0

      !--- box information
      real(KINDSF), private::m_BOXSHAPE(3,3)                           ! siumulation box shape
      real(KINDSF), private::m_BOXSIZE(3)                              ! simulation box size
      integer,      private::m_IFPD(3)                                 ! the periodic condition

      !--- box information on devices
      real(KINDSF), constant, private::dcm_BOXSHAPE(3,3)
      real(KINDSF), constant, private::dcm_BOXSIZE(3)
      integer,      constant, private::dcm_IFPD(3)

      !--- Voronoi tessellation vertices of atoms  on host
      integer, parameter, private::mp_INITED_VOLS = 1
      integer, parameter, private::mp_INITED_VVTS = 2
      integer, private::m_INITED_V = 0

      !--- storing the Voronoi veritice is very memory comsuming
      !    for large system, we have to divide the caluclation
      !    into a few segments. Th number of segments is determined
      !    by the number of atomc involved in one segment.
      integer, parameter, private::mp_MXNADEV = mp_BLOCKSIZE*mp_BLOCKDIMX
      integer, parameter, private::mp_MXKVOIS = 800
      integer, private::m_MXNADEV = mp_MXNADEV
      integer, private::m_MXNF    = 40                           !Actual max number  of faces of Voronio volumes
      integer, private::m_MXNV    = 40                           !Actual max number  of vertice of faces of Voronio volumes
      integer, private::m_MXKVOIS = mp_MXKVOIS

      !--- Delaunay tessellation vertices on devices
      integer, device, dimension(:),  allocatable,private::d1m_DVTS
      integer, device, dimension(:),  allocatable,private::d2m_DVTS
      integer, device, dimension(:),  allocatable,private::d3m_DVTS
      integer, device, dimension(:),  allocatable,private::d4m_DVTS
      integer, device, dimension(:),  allocatable,private::d5m_DVTS
      integer, device, dimension(:),  allocatable,private::d6m_DVTS
      !--- Working space on devices
      integer, device, dimension(:),  allocatable,private::d1m_WVFS
      integer, device, dimension(:),  allocatable,private::d2m_WVFS
      integer, device, dimension(:),  allocatable,private::d3m_WVFS
      integer, device, dimension(:),  allocatable,private::d4m_WVFS
      integer, device, dimension(:),  allocatable,private::d5m_WVFS
      integer, device, dimension(:),  allocatable,private::d6m_WVFS
      !--- Note: using integer(1) for dxm_WFLG inducing compiling fail
      !          on linux, thus I changed integer(1) to integer(2)
      !--- Note HQ 2016/06/27: change integer(2) to integer
      integer, device, dimension(:),  allocatable, private::d1m_WFLG
      integer, device, dimension(:),  allocatable, private::d2m_WFLG
      integer, device, dimension(:),  allocatable, private::d3m_WFLG
      integer, device, dimension(:),  allocatable, private::d4m_WFLG
      integer, device, dimension(:),  allocatable, private::d5m_WFLG
      integer, device, dimension(:),  allocatable, private::d6m_WFLG

      !--- Voronoi volume on devices
      real(KINDSF), device, dimension(:),  allocatable::d1m_VOLS
      real(KINDSF), device, dimension(:),  allocatable::d2m_VOLS
      real(KINDSF), device, dimension(:),  allocatable::d3m_VOLS
      real(KINDSF), device, dimension(:),  allocatable::d4m_VOLS
      real(KINDSF), device, dimension(:),  allocatable::d5m_VOLS
      real(KINDSF), device, dimension(:),  allocatable::d6m_VOLS

      !--- Voronoi vertices on devices
      real(KINDSF), device, dimension(:,:),  allocatable::d1m_VVTS
      real(KINDSF), device, dimension(:,:),  allocatable::d2m_VVTS
      real(KINDSF), device, dimension(:,:),  allocatable::d3m_VVTS
      real(KINDSF), device, dimension(:,:),  allocatable::d4m_VVTS
      real(KINDSF), device, dimension(:,:),  allocatable::d5m_VVTS
      real(KINDSF), device, dimension(:,:),  allocatable::d6m_VVTS

      !--- states of Voronoi volumes on devices
      integer, parameter::mp_VSTAT_ENCLOSED = 0
      integer, parameter::mp_VSTAT_OPENED   = 1
      integer, parameter::mp_AMASK_UNSEL = 0
      integer, parameter::mp_AMASK_SEL   = 1

      !--- Note: using integer(1) for dxm_VSTAT inducing compiling fail
      !          on linux, thus I change integer(1) to integer(2)
      integer, device, dimension(:),  allocatable::d1m_VSTAT
      integer, device, dimension(:),  allocatable::d2m_VSTAT
      integer, device, dimension(:),  allocatable::d3m_VSTAT
      integer, device, dimension(:),  allocatable::d4m_VSTAT
      integer, device, dimension(:),  allocatable::d5m_VSTAT
      integer, device, dimension(:),  allocatable::d6m_VSTAT

      !--- the AMARSK indicating if an atom is included in the calculation
      !    By default, all the atom thate are in-box should be included.
      !    The AMARSK value can be used when only part of the box needed
      !    to be analyzed.
      integer, dimension(:), allocatable::hm_AMASK
      integer, device, dimension(:),  allocatable::d1m_AMASK
      integer, device, dimension(:),  allocatable::d2m_AMASK
      integer, device, dimension(:),  allocatable::d3m_AMASK
      integer, device, dimension(:),  allocatable::d4m_AMASK
      integer, device, dimension(:),  allocatable::d5m_AMASK
      integer, device, dimension(:),  allocatable::d6m_AMASK
      integer, private, hm_NAINBOX             ! number of atoms in box

      !--- the HMARSK indicating if an atom is included in to be
      !    a Delaunay vertice. default, all the atom thate are in-box
      !    can a Delaunay vertice of a center atom.
      !    The HMASK value can be used, for example, to calculate the
      !    impuritie enviroment chanraterized by surrounding host atoms
      integer, dimension(:), allocatable::hm_HMASK
      integer, device, dimension(:),  allocatable::d1m_HMASK
      integer, device, dimension(:),  allocatable::d2m_HMASK
      integer, device, dimension(:),  allocatable::d3m_HMASK
      integer, device, dimension(:),  allocatable::d4m_HMASK
      integer, device, dimension(:),  allocatable::d5m_HMASK
      integer, device, dimension(:),  allocatable::d6m_HMASK

      !--- the interface to external procedure determine the atom that to be
      !    included in calculation. The
      abstract interface
        SUBROUTINE AMASKPROC(SimBox, CtrlParam, MASK)
        !***  PURPOSE:   to create the AMASK witch identify the atoms
        !                that will be included in calculation
        !     INPUT:     SimBox:
        !                CtrlParam:
        !
        !     OUTPUT:    FLAG,
        !
          use MD_CONSTANTS
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox),         intent(in)  ::SimBox
          type(SimMDCtrl),        intent(in)  ::CtrlParam
          integer,  dimension(:), intent(out) ::MASK
       END SUBROUTINE AMASKPROC
      end interface
      procedure(AMASKPROC), pointer, private::m_pAMaskProc=>null()
      procedure(AMASKPROC), pointer, private::m_pHMaskProc=>null()

      !--- the interface for a user define output process of voronoi volume
      abstract interface
        SUBROUTINE OUTPUTPROC(hFile, List, LATT, TVOLS, TCYST, TSTATU)
        !***  PURPOSE:   to putout the Voronoi volume to an output file.
        !     INPUT:  hFile,  unit of the output file
        !             LIST,   the neighbore list calculated by MD_NeighborsList2C_12_GPU.F90
        !             LATT,   the latice unit used here for unit converting.
        !             TVOLS, TSTATU, TCYST, the volume and state of Voronoi volume and structure fingerprint
        !
          use MD_Globle_Variables_GPU
          use MD_NeighborsList
          implicit none
          !---dummy vaiables
          integer, intent(in)::hFile
          type(NEIGHBOR_LIST),intent(in)::List
          real(KINDDF),intent(in)::LATT
          real(KINDSF), dimension(:)::TVOLS
          integer,     dimension(:)::TSTATU
          integer(8),  dimension(:)::TCYST

       END SUBROUTINE OUTPUTPROC
      end interface
      procedure(OUTPUTPROC), pointer, private::m_pOutputProc=>null()

      !----
      private::Allocate_WorkingArray_template,  &
               Allocate_WorkingArray_DEV,       &
               Allocate_VoronoiArray_template,  &
               Allocate_VoronoiVolume_DEV,      &
               Clear_WorkingArray_template,     &
               Clear_VoronoiVolume_template,    &
               StartOnDevice_Delaunay_template0,&
               StartOnDevice_Delaunay_template1,&
               StartOnDevice_Voronoi_template,  &
               COPYOUT_Delaunay_template,       &
               COPYOUT_Voronoi_template,        &
               COPYOUT_Voronoi_Volume_template, &
               LoadControlParameters,           &
               PrintControlParameters
  contains
  !********************************************************************************************
  logical function IfInit_DelaunayTesselation()result(YES)
  implicit none

          YES = (m_INITED .ne. 0)

          return
  end function IfInit_DelaunayTesselation
  !********************************************************************************************

  !********************************************************************************************
  logical function IfInit_VoronoiTesselation()result(YES)
  implicit none

          if( (m_INITED_V .gt. 0) .and. (m_INITED .gt. 0)) then
              YES = .true.
          else
              YES = .false.
          end if

          return
  end function IfInit_VoronoiTesselation
  !********************************************************************************************

  !********************************************************************************************
  integer function GetMaxNV()result(N)
  implicit none

          N = m_MXNV
          return
  end function GetMaxNV
  !********************************************************************************************

 !********************************************************************************************
  integer function GetMaxNF()result(N)
  implicit none

          N = m_MXNF
          return
  end function GetMaxNF
  !*****************************************************************************

  !**********************************************************************************
   subroutine CreateDefaultMask(SimBox, CtrlParam, Mask)
  !***  PURPOSE:  to create the default mask of atoms. The atoms out of box are excluded
  !               from the calculation
  !     NOTE:     the MASK is the MASK for the atoms that have beeen ordered in Nieghbore
  !               calculation
    use MD_Globle_Variables_GPU, only:dm_NPRT,hm_STATU
    implicit none
      type(SimMDBox),       intent(in)::SimBox
      type(SimMDCtrl),      intent(in)::CtrlParam
      integer, dimension(:)           ::Mask
      !--- dummy vaiables
      integer::I

         do I=1, dm_NPRT
            if(iand(hm_STATU(I),CP_STATU_OUTOFBOX) .eq. CP_STATU_OUTOFBOX) then
               Mask(I) = mp_AMASK_UNSEL
            else
               Mask(I)    = mp_AMASK_SEL
            end if
         end do
         return
    end subroutine CreateDefaultMask
  !**********************************************************************************

  !**********************************************************************************
   subroutine CreateMaskForAtoms(SimBox, CtrlParam, AtomID)
  !***  PURPOSE:  to create the mask of atoms of given ID. The atoms out of box are excluded
  !               from the calculation.
  !     INPUT:    AtomID, the ID of atoms to be included in calcualation
  !
  !     NOTE:     the MASK is the MASK for the atoms that have beeen ordered in Nieghbore
  !               calculation
    use MD_Globle_Variables_GPU, only:dm_NPRT,hm_STATU,hm_ITYP
    implicit none
      type(SimMDBox),       intent(in)::SimBox
      type(SimMDCtrl),      intent(in)::CtrlParam
      integer,dimension(:), intent(in)::AtomID
      !--- dummy vaiables
      integer::I

         hm_AMASK = mp_AMASK_UNSEL
         do I=1, dm_NPRT
            if(AtomID(I) .le. 0) exit
            if(iand(hm_STATU(AtomID(I)),CP_STATU_OUTOFBOX) .ne. CP_STATU_OUTOFBOX) then
                hm_AMASK(AtomID(I)) = mp_AMASK_SEL
            end if
          end do
          call Copyin_AMask()
         return
    end subroutine CreateMaskForAtoms
  !**********************************************************************************

  !**********************************************************************************
  subroutine SetAtomSelector(ASelector, HSelector)
  !***  DESCRIPTION: to set the user defrined FlagProc, which will be used to
  !                  create the flags of atoms that will be included in calculation
  !
  !     INPUT: ASelectot,  the subroutine provided by user for selecting the atoms for
  !                        which Voronoi volume to be calculated
  !            HSelectot,  the subroutine provided by user for selecting the 'host' atoms
  !                        which are candidates of Deluaney vertice of center atoms
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
        SUBROUTINE ASelector(SimBox, CtrlParam, MASK)
        !***  PURPOSE:   to create the FLAG witch identify the atoms
        !                that will be included in clustering
        !     INPUT:     SimBox:
        !                CtrlParam:
        !
        !     OUTPUT:    FLAG,
        !
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox),        intent(in) ::SimBox
          type(SimMDCtrl),       intent(in) ::CtrlParam
          integer,dimension(:),  intent(out)::MASK
       END SUBROUTINE ASelector
  end interface

  interface
        SUBROUTINE HSelector(SimBox, CtrlParam, MASK)
        !***  PURPOSE:   to create the FLAG witch identify the atoms
        !                that will be included in clustering
        !     INPUT:     SimBox:
        !                CtrlParam:
        !
        !     OUTPUT:    FLAG,
        !
          use MD_CONSTANTS
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox),       intent(in)  ::SimBox
          type(SimMDCtrl),      intent(in)  ::CtrlParam
          integer,dimension(:), intent(out) ::MASK
       END SUBROUTINE HSelector
  end interface
  optional::ASelector, HSelector
  !--- END INTERFACE --------------------------------

          if(present(ASelector)) then
             m_pAMaskProc  => ASelector
          else
             m_pAMaskProc  => null()
          end if

          if(present(HSelector)) then
             m_pHMaskProc  => HSelector
          else
             m_pHMaskProc  => null()
          end if

           return
   end subroutine SetAtomSelector
  !**********************************************************************************

  !**********************************************************************************
  subroutine Copyin_AMask()
  !***  DESCRIPTION: to create mask array on CPU and then copinto the device mask array
  !
  !    INPUT:
   use MD_Globle_Variables_GPU, only:m_NDEVICE, COPY_IN_SHIFT_template
   implicit none
   !---dummy vaiables
   !---local variables

           if(m_NDEVICE .GE. 1) then
               call COPY_IN_SHIFT_template  (1, d1m_AMASK, hm_AMASK)
            end if

           if(m_NDEVICE .GE. 2) then
               call COPY_IN_SHIFT_template  (2, d2m_AMASK, hm_AMASK)
            end if

           if(m_NDEVICE .GE. 3) then
               call COPY_IN_SHIFT_template  (3, d3m_AMASK, hm_AMASK)
            end if

           if(m_NDEVICE .GE. 4) then
               call COPY_IN_SHIFT_template  (4, d4m_AMASK, hm_AMASK)
            end if

           if(m_NDEVICE .GE. 5) then
               call COPY_IN_SHIFT_template  (5, d5m_AMASK, hm_AMASK)
            end if

           if(m_NDEVICE .GE. 6) then
               call COPY_IN_SHIFT_template  (6, d6m_AMASK, hm_AMASK)
            end if
            return
  end subroutine Copyin_AMask
  !**********************************************************************************

  !**********************************************************************************
  subroutine Copyin_HMask()
  !***  DESCRIPTION: to copy mask array on CPU into the device mask array, for Host atoms
  !
  !    INPUT:
   use MD_Globle_Variables_GPU, only:m_NDEVICE, COPY_IN_NOSHIFT_template
   implicit none
   !---dummy vaiables
   !---local variables

           if(m_NDEVICE .GE. 1) then
               call COPY_IN_NOSHIFT_template(1, d1m_HMASK, hm_HMASK)
            end if

           if(m_NDEVICE .GE. 2) then
               call COPY_IN_NOSHIFT_template(2, d2m_HMASK, hm_HMASK)
            end if

           if(m_NDEVICE .GE. 3) then
               call COPY_IN_NOSHIFT_template(3, d3m_HMASK, hm_HMASK)
            end if

           if(m_NDEVICE .GE. 4) then
               call COPY_IN_NOSHIFT_template(4, d4m_HMASK, hm_HMASK)
            end if

           if(m_NDEVICE .GE. 5) then
               call COPY_IN_NOSHIFT_template(5, d5m_HMASK, hm_HMASK)
            end if

           if(m_NDEVICE .GE. 6) then
               call COPY_IN_NOSHIFT_template(6, d6m_HMASK, hm_HMASK)
            end if
            return
  end subroutine Copyin_HMask
  !**********************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray_template(IDEV, NA, MXNF, MXNV, DVTS, MXKVOIS, WVFS, WFLG, NAPDEV, STATE, AMASK, HMASK)
  !***  PURPOSE:   to allocate working space on a device. The working space is needed
  !                for the calculation of Delaunay tesselletion. If Voronoi tessellation
  !                to be calculated, further call to Allocate_VoronoiArray_template is
  !                needed.
  !     SEE ALSO:  Allocate_VoronoiArray_template
  !
  !     INPUT:     IDEV,    the index  of device on which the memory to be allocated
  !                NA,      the max number of atoms in calculations
  !                MXNF,    the max number of faces of a Voronoi volume
  !                MXNV,    the max number of vertices on a face of a Voronoi volume
  !                MXKVOIS, the max number of neighbores of an atom
  !                NAPDEV,  the numnber of atoms on the device
  !
  !     OUPUT:     DVTS,    the device memory storing the vertices of Delaunay vertice
  !                WVFS,    the device memory storing the faces of Delaunay volume (working space)
  !                WFLG,    the device memory working space
  !                STAT,    the device memory storing the state (enclosed or opened) of Voronoi volume
  !                AMASK,   the device memory storing the flag indicating which atoms will be included in calculating Voronoi volume
  !                HMASK,   the device memory storing the flag indicating which atoms will be included to be the vertice of center atoms with AMASK = 1)
  !
      use MD_Globle_Variables_GPU, only:dm_NPRT
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV, NA, NAPDEV, MXNF, MXNV, MXKVOIS
           integer, device, dimension(:), allocatable::DVTS, WVFS
           integer, device, dimension(:), allocatable::WFLG, STATE, AMASK, HMASK

      !--- Local vairables
      integer::ERR


              !$$--- allocate memory for neighbore list on devices
              ERR = cudaSetDevice(IDEV)
              allocate(DVTS(NA*MXNF*MXNV), WVFS(mp_BLOCKSIZE*mp_BLOCKDIMX*MXNF), WFLG(mp_BLOCKSIZE*mp_BLOCKDIMX*MXKVOIS), STAT=err)
              if(ERR) goto 100

              allocate(STATE(NAPDEV), AMASK(NAPDEV), HMASK(dm_NPRT), STAT=err)
              if(ERR) goto 100
              return

 100          write(*,fmt="(A, I3)") "MDPSCU Error: Fail to allocate working memory in VORONOI_TESSLATION MODULE on device", IDEV
              write(*,fmt="(A)")     "                 Process to be stopped."
              stop
              return
  end subroutine Allocate_WorkingArray_template
  !*********************************************************************************

  !*********************************************************************************
    subroutine Allocate_WorkingArray_DEV()
  !***  PURPOSE:   to allocate working space on a devices. The working space is needed
  !                for the calculation of Delaunay tesselletion. If Voronoi tessellation
  !                to be calculated, further call to Allocate_VoronoiArray_template is
  !                needed.
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT:   the working spaces allocated
  !
  !     NOTE:   intialization of device memories must be called after
  !             calling Initialize_Globle_Variables_DEV, with the number
  !             of particle on each device is avaliable
  !             and also one call to Cal_NeighBoreList_DEV
  !             with the information of cells of the template
  !             system is available.

  use MD_Globle_Variables_GPU
  implicit none
      !--- dummy vaiables

      !--- Local variables
      integer::ERR, CURDEV

           m_MXNADEV = min(mp_MXNADEV, m_NAPDEV)
           ERR = cudaGetDevice(CURDEV)
           if(m_NDEVICE.GE.1) &
              call Allocate_WorkingArray_template(m_DEVICES(1), m_MXNADEV, m_MXNF, m_MXNV, d1m_DVTS, &
                                                  m_MXKVOIS, d1m_WVFS, d1m_WFLG,                     &
                                                  m_NAPDEV, d1m_VSTAT, d1m_AMASK, d1m_HMASK)

           if(m_NDEVICE.GE.2) &
              call Allocate_WorkingArray_template(m_DEVICES(2), m_MXNADEV, m_MXNF, m_MXNV, d2m_DVTS, &
                                                  m_MXKVOIS, d2m_WVFS, d2m_WFLG,                     &
                                                  m_NAPDEV, d2m_VSTAT, d2m_AMASK, d2m_HMASK)

           if(m_NDEVICE.GE.3) &
              call Allocate_WorkingArray_template(m_DEVICES(3), m_MXNADEV, m_MXNF, m_MXNV, d3m_DVTS, &
                                                  m_MXKVOIS, d3m_WVFS, d3m_WFLG,                     &
                                                  m_NAPDEV, d3m_VSTAT, d3m_AMASK, d3m_HMASK)

           if(m_NDEVICE.GE.4) &
              call Allocate_WorkingArray_template(m_DEVICES(4), m_MXNADEV, m_MXNF, m_MXNV, d4m_DVTS, &
                                                  m_MXKVOIS, d4m_WVFS, d4m_WFLG,                     &
                                                  m_NAPDEV, d4m_VSTAT, d4m_AMASK, d4m_HMASK)

           if(m_NDEVICE.GE.5) &
              call Allocate_WorkingArray_template(m_DEVICES(5), m_MXNADEV, m_MXNF, m_MXNV, d5m_DVTS, &
                                                  m_MXKVOIS, d5m_WVFS, d5m_WFLG,                     &
                                                  m_NAPDEV, d5m_VSTAT, d5m_AMASK, d5m_HMASK)

           if(m_NDEVICE.GE.6) &
              call Allocate_WorkingArray_template(m_DEVICES(6), m_MXNADEV, m_MXNF, m_MXNV, d6m_DVTS, &
                                                  m_MXKVOIS, d6m_WVFS, d6m_WFLG,                     &
                                                  m_NAPDEV, d6m_VSTAT, d6m_AMASK, d6m_HMASK)

            !$$--- ending the device operations
            allocate(hm_AMASK(dm_NPRT), hm_HMASK(dm_NPRT), STAT=ERR)
            if(ERR) then
              write(*,fmt="(A, I3)") "MDPSCU Error: Fail to allocate working memory in VORONOI_TESSLATION MODULE on host"
              write(*,fmt="(A)")     "              Process to be stopped."
              stop
            end if
            hm_AMASK = mp_AMASK_SEL
            hm_HMASK = mp_AMASK_SEL
            call Copyin_AMask()
            call Copyin_HMask()
            call SynchronizeDevices()
            ERR = cudaSetDevice(CURDEV)

            m_INITED = 1

           return
    end subroutine Allocate_WorkingArray_DEV
   !*********************************************************************************

  !****************************************************************************
  subroutine Allocate_VoronoiArray_template(IDEV, NA, MXNF, MXNV, VVTS, NAPDEV,VOLS)
  !***  PURPOSE:   to allocate working space on a device. The working space is needed
  !                for the calculation of Voronoi tesselletion, after Delaunay tesselletion.
  !                If no Voronoi tessellation to be calculated, it is unnecessary to
  !                call this routine.
  !                needed.
  !     INPUT:     IDEV,    the index  of device on which the memory to be allocated
  !                NA,      the max number of atoms in calculations
  !                MXNF,    the max number of faces of a Voronoi volume
  !                MXNV,    the max number of vertices on a face of a Voronoi volume
  !                NAPDEV,  the numnber of atoms on the device
  !
  !
  !     OUPUT:     VVTS,    the device memory storing the vertices position of atoms
  !                VOLS,    the device memory storing the volume of atoms
  !     NOTE:      NA may be different from NAPDEV
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV, NA, NAPDEV, MXNF, MXNV
           real(KINDSF), device, dimension(:), allocatable::VOLS
           real(KINDSF), device, dimension(:,:), allocatable::VVTS

      !--- Local vairables
           integer::ERR


              !$$--- allocate memory for neighbore list on devices
              ERR = cudaSetDevice(IDEV)

              if(.not. allocated(VVTS)) then
                  allocate(VVTS(NA*MXNF*MXNV,3), STAT=err)
                  if(ERR) goto 100
              end if

              if(.not.allocated(VOLS)) then
                 allocate(VOLS(NAPDEV), STAT=err)
                 if(ERR) goto 100
              end if
              return

 100          write(*,*) "MDPSCU ERROR: Fail to allocate memory in VORONOI_TESSLATION MODULE on device", IDEV
              stop
              return
  end subroutine Allocate_VoronoiArray_template
  !**********************************************************************************

  !*********************************************************************************
   subroutine Allocate_VoronoiVolume_DEV()
   !***  PURPOSE:   to allocate working space on a devices. The working space is needed
   !                for the calculation of Voronoi tesselletion, after Delaunay tesselletion.
   !                If no Voronoi tessellation to be calculated, it is unnecessary to
   !                call this routine.
   !                needed.
   !     INPUT:
   !
   !     OUTPUT:   the working spaces allocated
   !
    use MD_Globle_Variables_GPU
    implicit none
      !--- dummy vaiables

      !--- Local variables
      integer::ERR, CURDEV

           if(m_INITED .eq. 0) then
             write(*,fmt="(A)") " MDPSCU Error: Before calculate Voronoi volume, Delaunay Tesslation should be initalized"
             stop
           end if

           !$$--- allocate the memory allocated on device 1
           ERR = cudaGetDevice(CURDEV)
           if(m_NDEVICE.GE.1) &
               call Allocate_VoronoiArray_template(m_DEVICES(1), m_MXNADEV, m_MXNF, m_MXNV, d1m_VVTS, m_NAPDEV, d1m_VOLS)
           if(m_NDEVICE.GE.2) &
               call Allocate_VoronoiArray_template(m_DEVICES(2), m_MXNADEV, m_MXNF, m_MXNV, d2m_VVTS, m_NAPDEV, d2m_VOLS)

           if(m_NDEVICE.GE.3) &
               call Allocate_VoronoiArray_template(m_DEVICES(3), m_MXNADEV, m_MXNF, m_MXNV, d3m_VVTS, m_NAPDEV, d3m_VOLS)

           if(m_NDEVICE.GE.4) &
               call Allocate_VoronoiArray_template(m_DEVICES(4), m_MXNADEV, m_MXNF, m_MXNV, d4m_VVTS, m_NAPDEV, d4m_VOLS)

           if(m_NDEVICE.GE.5) &
               call Allocate_VoronoiArray_template(m_DEVICES(5), m_MXNADEV, m_MXNF, m_MXNV, d5m_VVTS, m_NAPDEV, d5m_VOLS)

           if(m_NDEVICE.GE.6) &
               call Allocate_VoronoiArray_template(m_DEVICES(6), m_MXNADEV, m_MXNF, m_MXNV, d6m_VVTS, m_NAPDEV, d6m_VOLS)

            !$$--- ending the device operations
            call SynchronizeDevices()
            ERR = cudaSetDevice(CURDEV)

            m_INITED_V = IOR(m_INITED_V, mp_INITED_VOLS)
            m_INITED_V = IOR(m_INITED_V, mp_INITED_VVTS)

           return
    END SUBROUTINE Allocate_VoronoiVolume_DEV
  !*********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray_template(IDEV, DVTS, VVTS, WVFS, WFLG)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_WorkingArray_template
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !
  !     OUPUT:
  !     NOTE:      _VSTATE and _AMASK to be deallocated by calling Clear_VoronoiVolume_template
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer,      device, dimension(:),   allocatable::DVTS, WVFS
           integer,      device, dimension(:),   allocatable::WFLG
           real(KINDSF), device, dimension(:,:), allocatable::VVTS

      !--- Local vairables
           integer::CURDEV,ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)
              if(allocated(DVTS)) then
                 deallocate(DVTS, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(VVTS)) then
                 deallocate(VVTS, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(WVFS)) then
                 deallocate(WVFS, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(WFLG)) then
                 deallocate(WFLG, STAT=ERR)
                 if(ERR) goto 100
              end if

             ERR = cudaSetDevice(CURDEV)
             return

   100       ERR = cudaSetDevice(CURDEV)
             write(*,*) "MDPSCU Warning: fail to deallocate memory on device in Voronoi Tesslation module", IDEV
             call ONWARNING(gm_OnWarning)

             return
  end subroutine Clear_WorkingArray_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_VoronoiVolume_template(IDEV, VOLS, STATE, AMASK, HMASK)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_VoronoiArray_template
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !     OUPUT:
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, device, dimension(:),      allocatable::STATE, AMASK, HMASK
           real(KINDSF), device, dimension(:), allocatable::VOLS

      !--- Local vairables
           integer::CURDEV,ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)
              if(allocated(VOLS)) then
                 deallocate(VOLS, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(STATE)) then
                 deallocate(STATE, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(AMASK)) then
                 deallocate(AMASK, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(HMASK)) then
                 deallocate(HMASK, STAT=ERR)
                 if(ERR) goto 100
              end if

             ERR = cudaSetDevice(CURDEV)
             return

   100       ERR = cudaSetDevice(CURDEV)
             write(*,*) "MDPSCU Warning: fail to deallocate memory on device in Voronoi Tesslation module", IDEV
             call ONWARNING(gm_OnWarning)
             return
  end subroutine Clear_VoronoiVolume_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray()
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calculation
  !    INPUT:
  !   OUTPUT:
  !
   use MD_Globle_Variables_GPU, only:m_DEVICES,m_NDEVICE
   implicit none

      !--- Local variables
      integer::err,CURDEV
      !---------------------------------------

           ERR = cudaGetDevice(CURDEV)

           !$$--- clear the memory allocated on device 1
           if(m_NDEVICE.GE.1) &
              call Clear_WorkingArray_template(m_DEVICES(1), d1m_DVTS, d1m_VVTS, d1m_WVFS, d1m_WFLG)

           !$$--- clear the memory allocated on device 2
           if(m_NDEVICE.GE.2) &
              call Clear_WorkingArray_template(m_DEVICES(2), d2m_DVTS, d2m_VVTS, d2m_WVFS, d2m_WFLG)

           !$$--- clear the memory allocated on device 3
           if(m_NDEVICE.GE.3) &
              call Clear_WorkingArray_template(m_DEVICES(3), d3m_DVTS, d3m_VVTS, d3m_WVFS, d3m_WFLG)

           !$$--- clear the memory allocated on device 4
           if(m_NDEVICE.GE.4) &
              call Clear_WorkingArray_template(m_DEVICES(4), d4m_DVTS, d4m_VVTS, d4m_WVFS, d4m_WFLG)

           !$$--- clear the memory allocated on device 5
           if(m_NDEVICE.GE.5) &
              call Clear_WorkingArray_template(m_DEVICES(5), d5m_DVTS, d5m_VVTS, d5m_WVFS, d5m_WFLG)

           !$$--- clear the memory allocated on device 6
           if(m_NDEVICE.GE.6) &
              call Clear_WorkingArray_template(m_DEVICES(6), d6m_DVTS, d6m_VVTS, d6m_WVFS, d6m_WFLG)

           ERR = cudaSetDevice(CURDEV)

           if(allocated(hm_AMASK)) deallocate(hm_AMASK)
           if(allocated(hm_HMASK)) deallocate(hm_HMASK)
           m_INITED = 0
      RETURN
  end subroutine Clear_WorkingArray
  !*********************************************************************************

  !**********************************************************************************
  subroutine Clear_VoronoiVolumeArray()
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calculation
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here

   use MD_Globle_Variables_GPU, only:m_DEVICES,m_NDEVICE
   implicit none
      !--- Local variables
           integer::err,CURDEV

           ERR = cudaGetDevice(CURDEV)

           !$$--- clear the memory allocated on device 1
           if(m_NDEVICE.GE.1) &
              call Clear_VoronoiVolume_template(m_DEVICES(1), d1m_VOLS, d1m_VSTAT, d1m_AMASK, d1m_HMASK)

           !$$--- clear the memory allocated on device 2
           if(m_NDEVICE.GE.2) &
              call Clear_VoronoiVolume_template(m_DEVICES(2), d2m_VOLS, d2m_VSTAT, d2m_AMASK, d2m_HMASK)

           !$$--- clear the memory allocated on device 3
           if(m_NDEVICE.GE.3) &
              call Clear_VoronoiVolume_template(m_DEVICES(3), d3m_VOLS, d3m_VSTAT, d3m_AMASK, d3m_HMASK)

           !$$--- clear the memory allocated on device 4
           if(m_NDEVICE.GE.4) &
              call Clear_VoronoiVolume_template(m_DEVICES(4), d4m_VOLS, d4m_VSTAT, d4m_AMASK, d4m_HMASK)

           !$$--- clear the memory allocated on device 5
           if(m_NDEVICE.GE.5) &
              call Clear_VoronoiVolume_template(m_DEVICES(5), d5m_VOLS, d5m_VSTAT, d5m_AMASK, d5m_HMASK)

           !$$--- clear the memory allocated on device 6
           if(m_NDEVICE.GE.6) &
              call Clear_VoronoiVolume_template(m_DEVICES(6), d6m_VOLS, d6m_VSTAT, d6m_AMASK, d6m_HMASK)

           ERR = cudaSetDevice(CURDEV)
           m_INITED_V = 0
      RETURN
  end subroutine Clear_VoronoiVolumeArray
  !*********************************************************************************

  !**********************************************************************************
  subroutine Clear_VoronoiTessellation_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calculation
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use MD_Globle_Variables_GPU, only:m_DEVICES,m_NDEVICE
   implicit none
       type(SimMDBox)::SimBox
       type(SimMDCtrl)::CtrlParam


      !--- Local variables
      integer::err,CURDEV
      !---------------------------------------

           ERR = cudaGetDevice(CURDEV)
           call Clear_WorkingArray()
           call Clear_VoronoiVolumeArray()
      RETURN
  end subroutine Clear_VoronoiTessellation_DEV
  !*********************************************************************************

  !*********************************************************************************
   subroutine LoadControlParameters(fname, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    CtrlParam,
  !
   use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables
    implicit none
    !--- dummy vaiables
    character*(*)  ::fname
    type(SimMDCtrl)::CtrlParam

   !--- local variables
    integer::hFile, N, I, LINE, NEEDDAMP, DAMPSCHEME
    character*256::STR
    character*32::STRNUMB(2), KEYWORD

            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old', err=200)
              LINE = 0
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case ("&SAVE_VOL", "&SAVE_VORVOL")
                             !*** get if save atomic volume

                         case( "&SAVE_DV")
                             !*** get if Delaunay vertice to be save
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_OUTPUTDT = 0
                             else
                                m_OUTPUTDT = ISTR(STRNUMB(1))
                             end if

                         case( "&SAVE_VV")
                             !*** get if Voronoi vertice to be save
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_OUTPUTVT = 0
                             else
                                m_OUTPUTVT = ISTR(STRNUMB(1))
                             end if

                         case( "&SAVE_SIG")
                             !*** get if signature to be output
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_OUTPUTST = 0
                             else
                                m_OUTPUTST = ISTR(STRNUMB(1))
                             end if

                         case( "&MXFACE")
                             !*** get maxmun number of faces
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_MXNF = 0
                             else
                                m_MXNF = ISTR(STRNUMB(1))
                             end if
                             if(m_MXNF .lt. 4) then
                                write(*,fmt="(A)")       " MDPSCU Error: the permitted number of faces of Voronoi volume should be larger 4"
                                write(*,fmt="(A, BZI6)") "              check control file at line:", LINE
                                write(*,fmt="(A)")       "               process to be stopped"
                                stop
                             end if
                         case( "&MXVERT")
                             !*** get maxmun number of faces
                             call Extract_Numb(STR,1, N, STRNUMB)
                             if(N .LT. 1) then
                                m_MXNV = 0
                             else
                                m_MXNV = ISTR(STRNUMB(1))
                             end if
                             if(m_MXNV .lt. 4) then
                                write(*,fmt="(A)")       " MDPSCU Error: the permitted number of vertice of Voronoi volume should be larger 4"
                                write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                                write(*,fmt="(A)")       "               process to be stopped"
                                stop
                             end if

                         case("&QUICKDUMP", "&QUICKDAMP")
                             !*** get the controal paraemter of output average CP values
                             call Extract_Numb(STR,1, N, STRNUMB)
                             NEEDDAMP = ISTR(STRNUMB(1))
                             call Extract_Substr(STR,1,N,STRNUMB)
                             if(N .ge. 1) then
                                call UpCase(STRNUMB(1))
                                if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "LBFGS") then
                                   DAMPSCHEME = CP_DAMPSCHEME_LBFGS
                                  else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG" ) then
                                    DAMPSCHEME = CP_DAMPSCHEME_CG
                                  else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG-LS") then
                                    DAMPSCHEME = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                    
                                  else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST" ) then
                                    DAMPSCHEME = CP_DAMPSCHEME_ST
                                  else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST-LS") then
                                    DAMPSCHEME = ior(CP_DAMPSCHEME_ST, CP_DAMPSCHEME_LSEARCH)                                     
                                  else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYN" .or. &
                                    STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "DYNAMICS" ) then
                                    DAMPSCHEME = CP_DAMPSCHEME_DYN
                                  else
                                    write(*,fmt="(A)")      " MDPSCU Error: the damping scheme "//STRNUMB(1)(1:len_trim(STRNUMB(1)))//" is unknown"
                                    write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                    write(*,fmt="(A)")      ' Process to be stopped'
                                    stop
                                end if
                             end if
                             call SetNeedDamp_SimMDCtrl(CtrlParam, NEEDDAMP, DAMPSCHEME)

                        case("&JOBSEL")
                             !$$*** To get job range to be analysis
                             call Extract_Numb(STR,3,n,STRNUMB)
                             CtrlParam%JOBID0 = ISTR(STRNUMB(1))
                             if(N .GE. 2) CtrlParam%JOBID1  = ISTR(STRNUMB(2))
                             if(N .GE. 3) CtrlParam%JOBIDSTEP = ISTR(STRNUMB(3))

                        case("&CFGSEL")
                            !$$*** To get cfg range to be analysis
                            call Extract_Numb(STR,3,n,STRNUMB)
                            CtrlParam%STARTCFG = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDCFG  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%CFGSTEP = ISTR(STRNUMB(3))

                        case("&BOXSEL")
                            !$$*** To get cfg range to be analysis
                            call Extract_Numb(STR,3,n,STRNUMB)
                            CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))

                         case( mp_FTAGO)
                              call Extract_Substr(STR,1,n,m_OUTFILE)

                         case default
                              write(*,*)" MDPSCU Warning: unknown keyword in VoronoiTessellation control file", KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)
                  end select
              end do
     100     close(hFile)

             if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for VoronoiTessellation calculation is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
             else
                call CreateDataFolder_Globle_Variables(m_OUTFILE)
             end if

              if(m_MXNV .lt. m_MXNF) then
                 write(*,fmt="(A)")       " MDPSCU Error: the permitted number of faces should be larger than the number of vertices"
                 write(*,fmt="(A, BZI6)") "               check the control file for keyword &MXVERT and &MXFACE"
                 write(*,fmt="(A)")       "               process to be stopped"
                 stop
              end if
            return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in VoronoiTessellation module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop
    end subroutine LoadControlParameters
  !*********************************************************************************

  !*********************************************************************************
    subroutine PrintControlParameters(hFile, SimBox, CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
  !
    use MD_Globle_Variables_GPU
    implicit none
      !--- dummy vaiables
      integer,         intent(in)::hFile
      type(SimMDBox),  intent(in)::SimBox
      type(SimMDCtrl), intent(in)::CtrlParam

      !--- Local variables

           !$$--- print out the control parameters
            write(hFile,fmt="(A)") " !************ VoronoiTessellationM_13_GPU module to be performed **********"
            write(hFile,fmt="(A)") " !    With the control parameters:"
            if(CtrlParam%NEEDDAMP .gt. 0) then
               write(hFile,FMT="(' !    Quickdampping to be performed before CSP calculation: ', I7, ' time steps')") CtrlParam%NEEDDAMP
            else
               write(hFile,FMT="(' !    Quickdampping to be performed before VoronoiTessellation: NO')")
            end if
            write(hFile,FMT="(' !    ')")

            if(m_OUTPUTVT .gt. 0) then
               write(hFile,fmt="(A)")  " !    Save intermediate Voronoi vertice:       YES"
            else
               write(hFile,fmt="(A)")  " !    Save intermediate Voronoi vertice:       NO"
            end if

            if(m_OUTPUTDT .gt. 0) then
               write(hFile,fmt="(A)")  " !    Save intermediate Delaunay vertices:     YES"
            else
               write(hFile,fmt="(A)")  " !    Save intermediate Delaunay vertices:     NO"
            end if

            if(m_OUTPUTST .gt. 0) then
               write(hFile,fmt="(A)")  " !    Save structure signatures:               YES"
            else
               write(hFile,fmt="(A)")  " !    Save structure signatures:               NO"
            end if

            write(hFile,fmt="(A, I3)") " !    Max number of faces of a Voronoi volume:  ", m_MXNF
            write(hFile,fmt="(A, I3)") " !    Max number of vertice of a Voronoi volume:", m_MXNV
            write(hFile,fmt="(A, I3)") " !    Max number of neighbores of a seed:       ", m_MXKVOIS
            write(hFile,fmt="(A)") " !    "
           return
    end subroutine PrintControlParameters
  !*********************************************************************************

  !*********************************************************************************
    subroutine Initialize_VoronoiTessellation_DEV(SimBox, CtrlParam)
    !***  PURPOSE:   to allocate the primaritive memory for Voronoi Tesselation
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !
    !     OUTPUT:   the working spaces allocated
    !
    use MD_Globle_Variables_GPU, only:Check_DEVICES
    use MD_TYPEDEF_PrintList,    only:m_pParamPrinter, Add_PrintProcess
    implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(in)::SimBox
      type(SimMDCtrl), intent(in)::CtrlParam

      !--- Local variables
      integer::I, IFILE

          !$$--- to check device version is used
          call Check_DEVICES()

           if(m_INITED .gt. 0) then
             call Clear_VoronoiTessellation_DEV(SimBox, CtrlParam)
           end if

           !$$--- to findout the I/O unit
           IFILE = 0
           do I=1, size(CtrlParam%f_tag)
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                 m_OUTFILE = CtrlParam%f_others(I)
              end if
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                 IFILE = I
              end if
           end do

            if(IFILE .le.0 ) then
               write(*,*) "MDPSCU Error: no control file for VoronoiTessellation calculation is given."
               write(*,*) "              add the keyword in SETUP file:",  mp_FTAGI
               write(*,fmt="(' Process to be stopped')")
               stop
           end if

           write(*,fmt="(A)") " !**** Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
           if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
           call LoadControlParameters(CtrlParam%f_others(IFILE),CtrlParam)
           m_MXNF = min(m_MXNF, CtrlParam%NB_MXNBS)
           m_MXNV = min(m_MXNV, CtrlParam%NB_MXNBS)
           m_MXKVOIS = min(CtrlParam%NB_MXNBS, mp_MXKVOIS)

           call Add_PrintProcess(m_pParamPrinter, PrintControlParameters)
           return
    end subroutine Initialize_VoronoiTessellation_DEV
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine ResetWorkingArray_KERNEL(NPRT, MXNN, MXNV, DVERT)
  !***  PURPOSE:   to reset working space into initial value
  !
  !$$   INPUT:     NPRT:   the actual number of atoms on the device
  !$$              MXNN:   the max permitted number of neighbores for each atom
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$
  !$$  OUTPUT:     DVERT:  the neighbore ID of atoms on Delaunay vertices
  !$$
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
           integer, value, intent(in)   ::NPRT, MXNN, MXNV
           integer, device, dimension(:)::DVERT

      ! Local variables
          integer::IT, IB, NT, NB, NN, I, J, NAL, NL, IA0


              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1)*griddim%x +  blockidx%x
              !$$--- size of Block
              NT = blockdim%x*blockdim%y
              !$$--- size of grid
              NB = griddim%x*griddim%y
              NAL = NT*NB

              NN = NPRT*MXNN*MXNV
              NL = (NN-1)/NAL+1

              do I=1, NL
                 IA0 = (I-1)*NAL+(IB-1)*NT+IT
                 if(IA0 .LT. NN) then
                    DVERT(IA0) = 0
                 end if
              end do
              return
  end subroutine ResetWorkingArray_KERNEL
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) subroutine ResetVoronoiArray_KERNEL(NPRT, MXNN, MXNV, VVERT)
  !***  PURPOSE:   to reset working space into initial value
  !
  !$$   INPUT:     NPRT:   the actual number of atoms on the device
  !$$              MXNN:   the max permitted number of neighbores for each atom
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$
  !$$  OUTPUT:     DVERT:  the neighbore ID of atoms on Delaunay vertices
  !$$
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
           integer, value, intent(in)::NPRT, MXNN, MXNV
           real(KINDSF), device, dimension(:,:)::VVERT

      ! Local variables
          integer::IT, IB, NT, NB, NN, I, J, NAL, NL, IA0


              IT  = (threadidx%y-1)*blockdim%x + threadidx%x
              IB  = (blockidx%y-1)*griddim%x +  blockidx%x
              !$$--- size of Block
              NT = blockdim%x*blockdim%y
              !$$--- size of grid
              NB = griddim%x*griddim%y
              NAL = NT*NB

              NN = NPRT*MXNN*MXNV
              NL = (NN-1)/NAL+1

              do I=1, NL
                 IA0 = (I-1)*NAL+(IB-1)*NT+IT
                 if(IA0 .LT. NN) then
                    VVERT(IA0, 1:3) = 1.E32
                 end if
              end do
              return
  end subroutine ResetVoronoiArray_KERNEL
  !*********************************************************************************


  !*********************************************************************************
   attributes(global) SUBROUTINE Cal_Delaunay_Vertice_KERNEL0(IM, XP, AMASK, HMASK, NAPDEV, IA0, ISHIFT0, NPRT, MXNF, KVOIS,INDI, MXNV, VERT)
  !***  PURPOSE:   to create the primary faces and vertices of Delaunay tesstellation
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              XP:     the positions of atoms
  !$$              AMASK:  the statu of the atoms, indicating if a atom to be included in calculation
  !$$              HMASK:  the statu of the atoms, indicating if a atom to be included in calculation
  !$$              NAPDEV: the max number of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              ISHIFT0:the shift relative to IA0 in this running of the kernela
  !$$              NPRT:   the actual number of atoms on the device
  !$$              MXNF:   the max permitted number of faces for each atom
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$
  !$$  OUTPUT:     VERT:   the neighbore ID of atoms on Delaunay vertices
  !$$
   use MD_CONSTANTS
   implicit none
   !---DUMMY VARIABLES
       integer, value, intent(in)::IM, NAPDEV, IA0, ISHIFT0, NPRT, MXNF, MXNV
       real(KINDDF), device, intent(in)::XP(IM,3)
       integer,      device, intent(in)::AMASK(NAPDEV)
       integer,      device, intent(in)::HMASK(NAPDEV)
       integer,      device, intent(in)::KVOIS(NAPDEV)
       integer,      device, intent(in)::INDI(NAPDEV,*)
       integer,      device::VERT(*)
       !integer, device::ddebug(*)

   !---Local variables
      real(KINDSF)::RR1, RR2, RR3, R2T, R2MI, RRMI, XC, YC, ZC, XCP, YCP, ZCP, NMX, NMY, NMZ, NMD
      real(KINDSF)::XY12, YZ12, ZX12, RX12, RY12, RZ12, DET2
      real(KINDSF)::SEPX, SEPY, SEPZ
      real(KINDSF)::X0, Y0, Z0
      real(KINDSF)::X1, Y1, Z1
      real(KINDSF)::X2, Y2, Z2
      real(KINDSF)::X3, Y3, Z3
      real(KINDSF)::X4, Y4, Z4
      real(KINDSF)::SEP(mp_MXKVOIS,3)

      real(KINDSF),shared::HBX, HBY, HBZ, BX, BY, BZ
      integer,     shared::IFPDX, IFPDY, IFPDZ
      integer,     shared::NT, NB

      real(KINDSF), parameter::scal=1.E8, underflow=1.E-8
      integer::IC, ICS, IT, IB, LL, NAL,NL
      integer::J, IW, IW1, IIW, IVERT0, IVERT, NVERT,I1, I2, I3, I4, LF
   !-------------

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              if(IT .EQ. 1) then
                 BX = dcm_BOXSIZE(1)*scal
                 BY = dcm_BOXSIZE(2)*scal
                 BZ = dcm_BOXSIZE(3)*scal

                 HBX = BX*C_HALF
                 HBY = BY*C_HALF
                 HBZ = BZ*C_HALF

                 IFPDX = dcm_IFPD(1)
                 IFPDY = dcm_IFPD(2)
                 IFPDZ = dcm_IFPD(3)

                 !$$--- size of Block
                 NT = blockdim%x*blockdim%y

                 !$$--- size of grid
                 NB = griddim%x*griddim%y
            end if
            call syncthreads()

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1
            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x

            do LL=1, NL
              !$$IC -- the id of the atom on the device
              IC= (IB-1)*NT+IT +(LL-1)*NAL
              IVERT0 = (IC-1)*MXNF*MXNV
              ICS = IC +ISHIFT0
              if(IC.LE.NPRT) then
                if(AMASK(ICS) .NE. mp_AMASK_SEL) then
                   cycle
                end if
                 !$$--- get the position of the atom IC
                  X0 = XP(ICS+IA0,1)*scal
                  Y0 = XP(ICS+IA0,2)*scal
                  Z0 = XP(ICS+IA0,3)*scal
                  IIW = KVOIS(ICS)

                  do IW =1, min(IIW, mp_MXKVOIS)
                     J=INDI(ICS,IW)
                     !$$--- set a large value for the seperation to exclude the atom from vertice
                     if(HMASK(J) .eq. mp_AMASK_UNSEL) then
                        SEP(IW,1:3) = 1.E31
                        cycle
                     end if

                     SEP(IW,1) = XP(J,1)*scal - X0
                     SEP(IW,2) = XP(J,2)*scal - Y0
                     SEP(IW,3) = XP(J,3)*scal - Z0

                     if(IFPDX.GT.0) then
                         if(ABS(SEP(IW,1)) .GT. HBX) SEP(IW,1) = SEP(IW,1) - SIGN(BX,SEP(IW,1))
                     end if

                     if(IFPDY.GT.0) then
                        if(ABS(SEP(IW,2)) .GT. HBY) SEP(IW,2) = SEP(IW,2) - SIGN(BY,SEP(IW,2))
                     end if

                     if(IFPDZ.GT.0) then
                        if(ABS(SEP(IW,3)) .GT. HBZ) SEP(IW,3) = SEP(IW,3) - SIGN(BZ,SEP(IW,3))
                     end if

                  end do

                 !NOTE: use the following statement may degrad the efficiency
                 !      VERT(IVERT0+1:IVERT0+MXNN*MXNV) = 0
                 !      thus we initialize VERT before lunch the kernel
                 !$$--- STEP1: scan the neighbores to findout the nearest neighbore of atom IC
                 !$$           and thus the first face of VT
                  R2MI = 1.E31
                  do IW=1, min(IIW, mp_MXKVOIS)
                     !$$--- To calculate the seperation between particle IC
                     !$$    and its IWth neighbore
                     SEPX = SEP(IW,1)
                     SEPY = SEP(IW,2)
                     SEPZ = SEP(IW,3)
                     RR1 = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                     !$$--- NOTE:DO we need to accont only the neighbore in potential cutoff reqion?
                     if(RR1 .LT. R2MI) then
                        R2MI = RR1
                        I1   = IW
                        X1 = SEPX
                        Y1 = SEPY
                        Z1 = SEPZ
                     end if
                  end do
                 !$$--- store the the associated atom ID on the first face
                  VERT(IVERT0+1)  = I1
                  RR1 = R2MI

                 !$$--- STEP2: scan the neighbores to findout the starting point
                 !$$
                  R2MI = 1.E31
                  RRMI = 1.E31
                  do IW=1, min(IIW, mp_MXKVOIS)
                     if(IW .EQ. I1) cycle

                     SEPX = SEP(IW,1)
                     SEPY = SEP(IW,2)
                     SEPZ = SEP(IW,3)
                     RR2 = (SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ)

                     !$$--- the foot from half*(X1, Y1, Z1) to the intersecting line
                     !$$    of face I1 and I2
                     XY12 = X1*SEPY - Y1*SEPX
                     YZ12 = Y1*SEPZ - Z1*SEPY
                     ZX12 = Z1*SEPX - X1*SEPZ
                     RX12 = RR1*SEPX - RR2*X1
                     RY12 = RR1*SEPY - RR2*Y1
                     RZ12 = RR1*SEPZ - RR2*Z1
                     DET2 = XY12*XY12 + YZ12*YZ12 + ZX12*ZX12
                     if(ABS(DET2) .LE. underflow) cycle

                     XC   = (RY12*XY12 - RZ12*ZX12)/DET2 - X1
                     YC   = (RZ12*YZ12 - RX12*XY12)/DET2 - Y1
                     ZC   = (RX12*ZX12 - RY12*YZ12)/DET2 - Z1

                     R2T  = XC*XC + YC*YC + ZC*ZC
                     if(R2T .LT. R2MI) then
                        R2MI  = R2T
                        I2    = IW
                        X2 = SEPX
                        Y2 = SEPY
                        Z2 = SEPZ
                        RRMI = RR2
                     end if
                  end do
                  VERT(IVERT0+2)  = I2

                  !$$--- get the real project foot:(X1, Y1, Z1)
                  !$$    on face I2
                  !$$    NOTE: the factor C_HALF is not included
                  RR2  = (X2*X2 + Y2*Y2 + Z2*Z2)
                  XY12 = X1*Y2 - X2*Y1
                  YZ12 = Y1*Z2 - Y2*Z1
                  ZX12 = Z1*X2 - Z2*X1
                  RX12 = RR1*X2 - RR2*X1
                  RY12 = RR1*Y2 - RR2*Y1
                  RZ12 = RR1*Z2 - RR2*Z1
                  DET2 = XY12*XY12 + YZ12*YZ12 + ZX12*ZX12

                  XC   = (RY12*XY12 - RZ12*ZX12)/DET2
                  YC   = (RZ12*YZ12 - RX12*XY12)/DET2
                  ZC   = (RX12*ZX12 - RY12*YZ12)/DET2

                  !$$--- STEP3; from the foot to findout the first vertex
                  !$$           on the intersecting line of face I1 and I2
                   R2MI = 1.E31
                   DO IW=1, min(IIW, mp_MXKVOIS)
                      if(IW .EQ. I1 .OR. IW.EQ.I2) cycle
                      SEPX = SEP(IW,1)
                      SEPY = SEP(IW,2)
                      SEPZ = SEP(IW,3)
                      RR3 = (SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ)

                      DET2 = SEPX*YZ12 + SEPY*ZX12 + SEPZ*XY12
                      XCP = (RR3*YZ12 + RY12*SEPZ - RZ12*SEPY)/DET2 - XC
                      YCP = (RR3*ZX12 + RZ12*SEPX - RX12*SEPZ)/DET2 - YC
                      ZCP = (RR3*XY12 + RX12*SEPY - RY12*SEPX)/DET2 - ZC
                      R2T =  XCP*XCP + YCP*YCP + ZCP*ZCP
                      if(R2T .LT. R2MI) then
                         R2MI  = R2T
                         I3 = IW
                      end if
                   END DO
                   VERT(IVERT0+3) = I3

                   !$$--- I2 is also the second face ID and the first associated atom
                   !$$    of the second face
                   !$$    I1 , I3, are also the assocaied atoms on the second face.
                   I4 = IVERT0+MXNV
                   VERT(I4+1) = I2
                   VERT(I4+2) = I3
                   VERT(I4+3) = I1

                   !$$--- I3 is also the second face ID and the first associated atom
                   !$$    of the third face
                   !$$    I2 , I3, are also the associated atoms on the third face.
                   I4 = I4 + MXNV
                   VERT(I4+1) = I3
                   VERT(I4+2) = I1
                   VERT(I4+3) = I2
              end if

            END DO  !END LOOP FOR LL
        RETURN
  END SUBROUTINE Cal_Delaunay_Vertice_KERNEL0
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) SUBROUTINE Cal_Delaunay_Vertice_KERNEL1(IM, XP, AMASK, HMASK, NAPDEV, IA0, ISHIFT0, NPRT, MXNF, MXNV, MXKVOIS, &
                                                              KVOIS, INDI, FACE, FLAG, VERT, STAT)
  !***  PURPOSE:   to create the faces and vertices of Delaunay vertice of Voronoi tesstellation
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              XP:     the positions of atoms
  !$$              AMASK:  the statu of the atoms, indicating if a atom to be included in calculation
  !$$              HMASK:  the statu of the atoms, indicating if a atom to be included in calculation
  !$$              NAPDEV: the max number of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              ISHIFT0:the shift relative to IA0 in this running of the kernela
  !$$              NPRT:   the actual number of atoms on the device
  !$$              MXNF:   the max permitted number of faces for each atom
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$              MXKVOIS:the max permitted number of neighbores
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$
  !$$  OUTPUT:     FACE:   the working space temperaily storing the face information of atoms
  !$$              FLAG:   the working space temperaily storing the vert information of atoms
  !$$              VERT:   the neighbore ID of atoms on Delaunay vertices
  !$$              STAT:   the indicator to indicate the Delaunay is enclosed or not
  !$$
   use MD_CONSTANTS
   implicit none
   !---DUMMY VARIABLES
       integer, value, intent(in)::IM, NAPDEV, IA0, ISHIFT0, NPRT, MXNF, MXNV, MXKVOIS
       real(KINDDF), device, intent(in)::XP(IM,3)
       integer,      device, intent(in)::AMASK(NAPDEV)
       integer,      device, intent(in)::HMASK(NAPDEV)
       integer,      device, intent(in)::KVOIS(NAPDEV)
       integer,      device, intent(in)::INDI(NAPDEV,*)
       integer,      device::FACE(*), VERT(*)
       integer,      device::FLAG(*), STAT(*)
       !integer, device::ddebug(*)

   !---Local variables
       real(KINDSF)::RR1, RR2, RR3, R2T, R2MI, XC, YC, ZC, XCP, YCP, ZCP, NMX, NMY, NMZ, NMD
       real(KINDSF)::XY12, YZ12, ZX12, RX12, RY12, RZ12, DET2
       real(KINDSF)::SEPX, SEPY, SEPZ
       real(KINDSF)::X0, Y0, Z0
       real(KINDSF)::X1, Y1, Z1
       real(KINDSF)::X2, Y2, Z2
       real(KINDSF)::X3, Y3, Z3
       real(KINDSF)::X4, Y4, Z4
       real(KINDSF)::SEP(mp_MXKVOIS,3)

       real(KINDSF),shared::HBX, HBY, HBZ, BX, BY, BZ
       integer,     shared::IFPDX, IFPDY, IFPDZ
       integer,     shared::NT, NB

      ! Local variables
       real(KINDSF), parameter::scal=1.E8, underflow=1.E-8
       integer::IC, ICS, IT, IB, LL, NAL, NL
       integer::J, IW, IW1, IIW, IFACE0, IVERT0, IFLG0, IVERT, I1, I2, I3, I4, LF, IFW
       integer::NVERT, MXNV2, STATU

   !--------------
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              if(IT .EQ. 1) then
                 BX = dcm_BOXSIZE(1)*scal
                 BY = dcm_BOXSIZE(2)*scal
                 BZ = dcm_BOXSIZE(3)*scal

                 HBX = BX*C_HALF
                 HBY = BY*C_HALF
                 HBZ = BZ*C_HALF

                 IFPDX = dcm_IFPD(1)
                 IFPDY = dcm_IFPD(2)
                 IFPDZ = dcm_IFPD(3)

                 !$$--- size of Block
                 NT = blockdim%x*blockdim%y

                 !$$--- size of grid
                 NB = griddim%x*griddim%y
            end if
            call syncthreads()

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1
            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
            MXNV2 = MXNV-2

            do LL=1, NL
              !$$IC -- the id of the atom on in the calculation
              IC     =  (IB-1)*NT+IT +(LL-1)*NAL
              IFLG0  = ((IB-1)*NT+IT-1)*MXKVOIS
              IFACE0 = ((IB-1)*NT+IT-1)*MXNF
              IVERT0 =  (IC-1)*MXNF*MXNV
              ICS    =  IC+ISHIFT0

              if(IC.LE.NPRT) then
                 if(AMASK(ICS) .NE. mp_AMASK_SEL) then
                    STAT(ICS) =  mp_VSTAT_OPENED
                    cycle
                 end if

                 !$$--- get the position of the atom IC
                  IIW = KVOIS(ICS)
                  if(IIW .LE. 0) then
                     STAT(ICS) =  mp_VSTAT_OPENED
                     cycle
                  end if
                  X0 = XP(ICS+IA0,1)*scal
                  Y0 = XP(ICS+IA0,2)*scal
                  Z0 = XP(ICS+IA0,3)*scal

                 !$$-- reset working space
                 !$$   NOTE by HQ 2017/06/27:
                 !$$         on some machine and cuda8.0,  the following statement cannot work correctly
                 !$$
                 !$$   FACE(IFACE0+1:IFCAE0+3) = VERT(IVERT0+1:IVERT0+3)
                 !$$   FACE(IFACE0+4:IFACE0+MXNF) = 0
                 !$$   We change the statement to:
                  FACE(IFACE0+1:IFACE0+MXNF) = 0
                  FACE(IFACE0+1) = VERT(IVERT0+1)
                  FACE(IFACE0+2) = VERT(IVERT0+2)
                  FACE(IFACE0+3) = VERT(IVERT0+3)

                  STATU =  mp_VSTAT_ENCLOSED
                  do IW =1, min(IIW, mp_MXKVOIS)
                     J=INDI(ICS,IW)
                     !$$--- set a large value for the seperation to exclude the atom from vertice
                     if(HMASK(J) .eq. mp_AMASK_UNSEL) then
                        SEP(IW,1) = 1.E31
                        SEP(IW,2) = 1.E31
                        SEP(IW,3) = 1.E31
                        cycle
                     end if

                     SEP(IW,1) = XP(J,1)*scal - X0
                     SEP(IW,2) = XP(J,2)*scal - Y0
                     SEP(IW,3) = XP(J,3)*scal - Z0
                     if(IFPDX.GT.0) then
                         if(ABS(SEP(IW,1)) .GT. HBX) SEP(IW,1) = SEP(IW,1) - SIGN(BX,SEP(IW,1))
                     end if

                     if(IFPDY.GT.0) then
                        if(ABS(SEP(IW,2)) .GT. HBY) SEP(IW,2) = SEP(IW,2) - SIGN(BY,SEP(IW,2))
                     end if

                     if(IFPDZ.GT.0) then
                        if(ABS(SEP(IW,3)) .GT. HBZ) SEP(IW,3) = SEP(IW,3) - SIGN(BZ,SEP(IW,3))
                     end if
                  end do

                  !$$**********************************************************
                  !$$---  Now we have the first three faces and starting vertex
                  !$$     scan the faces to construct the faces

                   do LF =1, min(IIW, MXNF)
                      !$$--- check if we have further face needed to be constructed
                       if(FACE(IFACE0+LF) .LE. 0) exit
                      !$$--- get the position of the first vertex (XC, YC, ZC) on face LF
                       IVERT = IVERT0 + (LF-1)*MXNV
                       I1 = VERT(IVERT + 1)
                       I2 = VERT(IVERT + 2)
                       I3 = VERT(IVERT + 3)

                      !$$--- reset the working space for visited neighbors
                       FLAG(IFLG0+1:IFLG0+min(IIW,MXKVOIS)) = 0
                       FLAG(IFLG0+I1) = 1
                       FLAG(IFLG0+I2) = 1
                       FLAG(IFLG0+I3) = 1

                      !$$--- get the position (XC, YC, ZC) of the first vertex
                       X1 = SEP(I1,1)
                       Y1 = SEP(I1,2)
                       Z1 = SEP(I1,3)

                       X2 = SEP(I2,1)
                       Y2 = SEP(I2,2)
                       Z2 = SEP(I2,3)

                       X3 = SEP(I3,1)
                       Y3 = SEP(I3,2)
                       Z3 = SEP(I3,3)

                       RR1  = (X1*X1 + Y1*Y1 + Z1*Z1)
                       RR2  = (X2*X2 + Y2*Y2 + Z2*Z2)
                       RR3  = (X3*X3 + Y3*Y3 + Z3*Z3)
                       DET2 =  Y2*(X1*Z3 - X3*Z1)  + Y3*(X2*Z1 - X1*Z2)  + Y1*(X3*Z2 - X2*Z3)
                       XC   = (RR1*(Y2*Z3 - Y3*Z2) + RR2*(Y3*Z1 - Y1*Z3) + RR3*(Y1*Z2 - Y2*Z1))/DET2
                       YC   = (RR1*(Z2*X3 - Z3*X2) + RR2*(Z3*X1 - Z1*X3) + RR3*(Z1*X2 - Z2*X1))/DET2
                       ZC   = (RR1*(X2*Y3 - X3*Y2) + RR2*(X3*Y1 - X1*Y3) + RR3*(X1*Y2 - X2*Y1))/DET2

                      !$$--- prepare to find the next vertex
                      !
                      !$$--- the normal of the face I2
                       NMD = SQRT(RR2)
                       NMX = X2/NMD
                       NMY = Y2/NMD
                       NMZ = Z2/NMD

                      !$$--- the determinants to be used to find the next vertex
                       X2  = X3
                       Y2  = Y3
                       Z2  = Z3
                       RR2 = RR3
                       XY12 = X1*Y2 - X2*Y1
                       YZ12 = Y1*Z2 - Y2*Z1
                       ZX12 = Z1*X2 - Z2*X1
                       RX12 = RR1*X2 - RR2*X1
                       RY12 = RR1*Y2 - RR2*Y1
                       RZ12 = RR1*Z2 - RR2*Z1

                       !$$--- STEP4: scan the neighbores to construct other vertex
                       !$$           on face LF
                        IVERT   = IVERT + 3
                        NVERT   = 1

                        do IW1=1, MXNV !--- change IIW  to MXNV, because the permitte number of vertices is MXNV
                           R2MI = 1.E32
                           I4   = -1
                           IFW  = IFLG0
                           do IW=1, min(IIW, mp_MXKVOIS) !--- the number of candidate vertices is the number of neighbore
                              !$$--- check if the neighbor has been on the vertex list
                               IFW = IFW + 1
                               if(NVERT .LT. 2) then
                                  if(FLAG(IFW) .GT. 0) cycle
                               else
                                  if(FLAG(IFW) .GT. 0 .AND. IW.NE.I2) cycle
                               endif

                              !$$--- IW is not still on the vertex list
                               SEPX = SEP(IW,1)
                               SEPY = SEP(IW,2)
                               SEPZ = SEP(IW,3)

                              !$$---
                               RR3  = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                               DET2 = SEPX*YZ12 + SEPY*ZX12 + SEPZ*XY12
                               XCP  = (RR3*YZ12 + RY12*SEPZ - RZ12*SEPY)/DET2
                               YCP  = (RR3*ZX12 + RZ12*SEPX - RX12*SEPZ)/DET2
                               ZCP  = (RR3*XY12 + RX12*SEPY - RY12*SEPX)/DET2
                              !$$--- check if the vertex and I0 are on the same of the
                              !$$    of current moving face
                               if(NMX*XCP + NMY*YCP + NMZ*ZCP .GT. NMD) cycle

                               R2T  =  (XCP - XC)*(XCP - XC) + (YCP - YC)*(YCP - YC) + (ZCP - ZC)*(ZCP - ZC)
                               if(R2T .LT. R2MI) then
                                  R2MI  = R2T
                                  I4 = IW
                                  X4 = SEPX
                                  Y4 = SEPY
                                  Z4 = SEPZ
                                  X3 = XCP
                                  Y3 = YCP
                                  Z3 = ZCP
                               end if
                           end do !END LOOP FOR IW

                           !$$--- check if the number of vertice to be larger then permitted number
                           NVERT = NVERT + 1
                           if(NVERT .GT. MXNV2) then
                              STATU  = mp_VSTAT_OPENED
                              exit
                           end if

                           !$$--- add the new vertex to the vertex list of face LF
                           IVERT  = IVERT + 1
                           VERT(IVERT) = I4

                           !$$--- if all neighbore have been exhausted, exit
                           if(I4 .LT. 0 ) exit  !
                           !$$--- if the vertice forming an enclosed polygon, exit
                           if(I4 .EQ. I2) exit

                           !$$--- add I4 to face list if I4 is still not on the list
                            FLAG(IFLG0+I4) = 1
                            do J=1, MXNF
                               IW = J+IFACE0
                               if(FACE(IW).EQ. I4) then
                                  exit
                               else
                                if(FACE(IW).LE.0) then !-- having new vertex
                                   FACE(IW) = I4
                                  !$$--- set the first vertex on the face I4
                                   IW = IVERT0+(J-1)*MXNV
                                   VERT(IW+1) = I4
                                   VERT(IW+2) = I3
                                   VERT(IW+3) = I1
                                  exit
                                end if
                               end if
                            end do

                            !$$--- prepare for the next vertex
                            !$$--- store the positon of current vertex
                            NMD = SQRT(X2*X2+Y2*Y2+Z2*Z2)
                            NMX = X2/NMD
                            NMY = Y2/NMD
                            NMZ = Z2/NMD
                            I3  = I4
                            XC  = X3
                            YC  = Y3
                            ZC  = Z3
                            X2  = X4
                            Y2  = Y4
                            Z2  = Z4
                            RR2 = X2*X2 + Y2*Y2 + Z2*Z2
                            XY12 = X1*Y2 - X2*Y1
                            YZ12 = Y1*Z2 - Y2*Z1
                            ZX12 = Z1*X2 - Z2*X1
                            RX12 = RR1*X2 - RR2*X1
                            RY12 = RR1*Y2 - RR2*Y1
                            RZ12 = RR1*Z2 - RR2*Z1
                        end do !END LOOP FOR IW1

                        if(STATU .EQ. mp_VSTAT_OPENED) exit
                     end do !END LOOP FOR LF
                     STAT(ICS) = STATU
              end if
            END DO  !END LOOP FOR LL
        RETURN
  END SUBROUTINE Cal_Delaunay_Vertice_KERNEL1
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE StartOnDevice_Delaunay_template0(IDEV, CELL0, STARTCELL, ENDCELL, KVOIS,INDI,ITYP, XP, AMASK, HMASK, DVTS)
  !***  PURPOSE:   to create the faces and vertices of Delaunay vertice of
  !                Voronoi tesstellation
  !
  !     INPUT:     IDEV,      the ID of device
  !                CELL0,     the ID of the first cell on the device
  !                STARTCELL, the ID of the first cell on in this call
  !                ENDCELL,   the ID of the last cell on  in this call
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                AMARK:  the statu of the atoms, indicating if a atom to be included in calculation
  !                HMARK:  the statu of the host atoms, indicating if a host atom to be included in calculation
  !
  !     OUTPUT:    DVTS,    the vertices of Delaunay tessellation
  !
   use MD_Globle_Variables_GPU
   implicit none
   !---dummy vaiables
       integer::IDEV, CELL0, STARTCELL, ENDCELL
       integer,      device, dimension(:)  ::KVOIS
       integer,      device, dimension(:,:)::INDI
       integer,      device, dimension(:)  ::ITYP
       real(KINDDF), device, dimension(:,:)::XP
       integer,      device, dimension(:)  ::AMASK
       integer,      device, dimension(:)  ::HMASK
       integer,      device, dimension(:)  ::DVTS

   !---Local variables
       integer::BX, BY, NB, NBX, NBY, ERR, IA0, IAS, ENDA, NPRT
   !--- Device variables and variables to be used in GPU
        type(dim3) :: blocks
        type(dim3) :: threads

             ERR = cudaSetDevice(IDEV)
             !$$--- If pressure effect acounted for, the box size could be changed, we
             !$$    need to recopy the boxsize
             IA0   = hm_IA1th(CELL0)
             IAS   = hm_IA1th(STARTCELL)-IA0
             ENDA  = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
             NPRT  = ENDA -  hm_IA1th(STARTCELL) + 1
             !if the cell is empty, no caclculation performed
             if(NPRT .LE. 0) return
             dcm_BOXSHAPE = m_BOXSHAPE
             dcm_IFPD     = m_IFPD
             dcm_BOXSIZE  = m_BOXSIZE

             BX  = mp_BLOCKSIZE
             BY  = 1
             NBX = mp_BLOCKDIMX
             NBY = 1

             blocks  = dim3(NBX, NBY, 1)
             threads = dim3(BX, BY, 1)
             call Cal_Delaunay_Vertice_KERNEL0<<<blocks, threads>>>(dm_NPRT,XP, AMASK, HMASK, m_NAPDEV, IA0-1, IAS, NPRT, m_MXNF, &
                                                                    KVOIS, INDI, m_MXNV, DVTS)


          return
  end SUBROUTINE StartOnDevice_Delaunay_template0
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE StartOnDevice_Delaunay_template1(IDEV, CELL0, STARTCELL, ENDCELL, KVOIS,INDI,ITYP, XP, AMASK, HMASK, DVTS, STAT, WVFS, WFLG)
  !***  PURPOSE:   to create the faces and vertices of Delaunay vertice of
  !                Voronoi tesstellation
  !
  !     INPUT:     IDEV,      the ID of device
  !                CELL0,     the ID of the first cell on the device
  !                STARTCELL, the ID of the first cell on in this call
  !                ENDCELL,   the ID of the last cell on  in this call
  !
  !                KVOIS:     the number of neighbors for atoms
  !                INDI:      the index for the neighbores
  !                ITYP :     the type of atoms corresponding to INDI
  !                XP:        the positions of atoms
  !                AMASK:    the statu of the atoms, indicating if a atom to be included in calculation
  !                HMASK:    the statu of the atoms, indicating if a atom to be included in calculation
  !
  !     OUTPUT:    TVS,    the vertices of Delaunay tessellation
  !                WVFS,   the working space in creating  Delaunay vertice
  !
   use MD_Globle_Variables_GPU
   implicit none
   !---dummy vaiables
       integer, intent(in)::IDEV, CELL0, STARTCELL, ENDCELL
       integer,      device, dimension(:)  ::KVOIS
       integer,      device, dimension(:,:)::INDI
       integer,      device, dimension(:)  ::ITYP
       real(KINDDF), device, dimension(:,:)::XP
       integer,      device, dimension(:)  ::AMASK
       integer,      device, dimension(:)  ::HMASK
       integer,      device, dimension(:)  ::DVTS, WVFS
       integer,      device, dimension(:)  ::STAT, WFLG
   !---Local variables
       integer::BX, BY, NB, NBX, NBY, ERR, IA0, IAS, ENDA, NPRT

   !---Device variables and variables to be used in GPU
       type(dim3) :: blocks
       type(dim3) :: threads


             ERR = cudaSetDevice(IDEV)
             IA0   = hm_IA1th(CELL0)
             IAS   = hm_IA1th(STARTCELL)-IA0
             ENDA  = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
             NPRT  = ENDA -  hm_IA1th(STARTCELL) + 1
             !if the cell is empty, no caclculation performed
             if(NPRT .LE. 0) return

             BX  = mp_BLOCKSIZE
             BY  = 1
             NBX = mp_BLOCKDIMX
             NBY = 1
             blocks  = dim3(NBX, NBY, 1)
             threads = dim3(BX, BY, 1)
             call Cal_Delaunay_Vertice_KERNEL1<<<blocks, threads>>>(dm_NPRT,XP, AMASK, HMASK, m_NAPDEV, IA0-1, IAS, NPRT, m_MXNF, m_MXNV,&
                                                                    m_MXKVOIS, KVOIS, INDI, WVFS, WFLG, DVTS, STAT)


          return
  end SUBROUTINE StartOnDevice_Delaunay_template1
  !*********************************************************************************

  !*********************************************************************************
   attributes(global) SUBROUTINE Cal_Voronoi_KERNEL(IM, XP, NAPDEV, IA0, ISHIFT0, NPRT, MXNF, KVOIS,INDI, MXNV, DVERT, STAT, VOL, VVERT)
  !***  PURPOSE:   to calculate the Voronoi volume of atoms. Called after calling
  !$$              Cal_Delaunay_Vertice_KERNEL
  !$$
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              ISHIFT0:the shifr relative to IA0 in this running of the kernela
  !$$              NPRT:   the actual number of atoms on the device
  !$$              MXNF:   the max permitted number of faces for each atom
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !$$              MXNV:   the max permitted number of the vertice for each Voronoi volume of an atom
  !$$              DVERT:  the neighbore ID of atoms on Delaunay vertices
  !$$              STAT:   the stat indicating if the Voronoi volume are enclosed
  !$$
  !$$   OUTPUT:    VOL,    Voronoi volume of atoms
  !$$              VVERT,  Voronoi vertice of atoms
  !$$
   use MD_CONSTANTS
   implicit none
   !---DUMMY VARIABLES
       integer, value, intent(in)::IM, NAPDEV, IA0, ISHIFT0, NPRT, MXNF, MXNV
       real(KINDDF), device, intent(in)::XP(IM,3)
       integer,      device, intent(in)::KVOIS(NAPDEV)
       integer,      device, intent(in)::INDI(NAPDEV,*)
       integer,      device, intent(in)::DVERT(*)
       integer,      device, intent(in)::STAT(*)
       real(KINDSF), device, dimension(:,:)::VVERT
       real(KINDSF), device, dimension(:)::VOL

   !---Local variables
       real(KINDSF)::RR1, RR2, RR3, R2T, R2MI, VP1X, VP1Y, VP1Z, VP2X, VP2Y, VP2Z, VP3X, VP3Y, VP3Z
       real(KINDSF)::XY12, YZ12, ZX12, RX12, RY12, RZ12, DET2
       real(KINDSF)::X0, Y0, Z0
       real(KINDSF)::X1, Y1, Z1
       real(KINDSF)::X2, Y2, Z2
       real(KINDSF)::X3, Y3, Z3
       real(KINDSF)::SEP(mp_MXKVOIS,3)
       real(KINDSF)::VVV

       real(KINDSF),shared::HBX, HBY, HBZ, BX, BY, BZ
       integer,     shared::IFPDX, IFPDY, IFPDZ
       integer,     shared::NT, NB

      ! Local variables
       integer::IC, ICS, IT, IB, LL, NAL,NL
       integer::J, IW, IW1, IIW, IFACE0, IFACE, IVERT0, IVERT, NVERT,I1, I2, I3, I4, MXNV4, LF
       !scal: convert centimeter to A, to prevent underflow in single presion
       !p_VFACT: will conert the volume in A^3 to volume in CM^3
       !         1/48 = (1/3)*(1/2)(a.b)*(1/2)^3
       !         (1/2)^3 come from that fact that in calculating V-vertice, the
       !         distance between the seed and the D-vertice should multiplied
       !         (1/2)
       real(KINDSF), parameter::scal=1.E8, underflow=1.E-8, rescal = 0.5E-8
      !real(KINDSF), parameter::p_VFACT = 1.D0/48.D0    ! used without scal
       real(KINDSF), parameter::p_VFACT = 1.E-24/48.E0   ! used with scal

            IT  = (threadidx%y-1)*blockdim%x + threadidx%x

            if(IT .EQ. 1) then
               BX = dcm_BOXSIZE(1)*scal
               BY = dcm_BOXSIZE(2)*scal
               BZ = dcm_BOXSIZE(3)*scal

               HBX = BX*C_HALF
               HBY = BY*C_HALF
               HBZ = BZ*C_HALF

               IFPDX = dcm_IFPD(1)
               IFPDY = dcm_IFPD(2)
               IFPDZ = dcm_IFPD(3)

              !$$--- size of Block
               NT = blockdim%x*blockdim%y

              !$$--- size of grid
               NB = griddim%x*griddim%y
            end if
            call syncthreads()

            !$$--- how many loop needed
            !$$    NOTE: don't know why it don't work correctly
            !$$    if NAL, NL were used as shared
            NAL = NT*NB
            NL = (NPRT-1)/NAL+1
            IB  =  (blockidx%y-1) * griddim%x +  blockidx%x
            MXNV4 = MXNV-4

            do LL=1, NL
              !$$IC -- the id of the atom on the device
              IC= (IB-1)*NT+IT +(LL-1)*NAL
              IFACE0 = (IC-1)*MXNF
              IVERT0 = (IC-1)*MXNF*MXNV
              ICS = IC+ISHIFT0

              if(IC.LE.NPRT) then
                VVV = 0.D0
               !NOTE: use the following statement may degrad the efficiency
               !      VVERT(IVERT0+1:IVERT0+MXNN*MXNV,1:3) = 1.E32
               !      thus we initialize VERT before luch the kernel


              if(STAT(ICS) .EQ. mp_VSTAT_ENCLOSED) then
                 !$$--- get the position of the atom IC
                  X0 = XP(ICS+IA0,1)*scal
                  Y0 = XP(ICS+IA0,2)*scal
                  Z0 = XP(ICS+IA0,3)*scal
                  IIW = KVOIS(ICS)
                  do IW =1, min(IIW, mp_MXKVOIS)
                     J=INDI(ICS,IW)
                     SEP(IW,1) = XP(J,1)*scal - X0
                     SEP(IW,2) = XP(J,2)*scal - Y0
                     SEP(IW,3) = XP(J,3)*scal - Z0
                     if(IFPDX.GT.0) then
                         if(ABS(SEP(IW,1)) .GT. HBX) SEP(IW,1) = SEP(IW,1) - SIGN(BX,SEP(IW,1))
                     end if

                     if(IFPDY.GT.0) then
                        if(ABS(SEP(IW,2)) .GT. HBY) SEP(IW,2) = SEP(IW,2) - SIGN(BY,SEP(IW,2))
                     end if

                     if(IFPDZ.GT.0) then
                        if(ABS(SEP(IW,3)) .GT. HBZ) SEP(IW,3) = SEP(IW,3) - SIGN(BZ,SEP(IW,3))
                     end if
                  end do

                  do LF =1, MXNF
                     !$$--- get the position (VP1x, Vp1Y, VP1Z) of the first vertex on face LF
                     !$$
                      IVERT = IVERT0 + (LF-1)*MXNV
                      I1 = DVERT(IVERT + 1)
                      if(I1 .LE. 0) exit
                      I2 = DVERT(IVERT + 2)
                      I3 = DVERT(IVERT + 3)
                      X1 = SEP(I1,1)
                      Y1 = SEP(I1,2)
                      Z1 = SEP(I1,3)

                      X2 = SEP(I2,1)
                      Y2 = SEP(I2,2)
                      Z2 = SEP(I2,3)

                      X3 = SEP(I3,1)
                      Y3 = SEP(I3,2)
                      Z3 = SEP(I3,3)

                      RR1  = (X1*X1 + Y1*Y1 + Z1*Z1)
                      RR2  = (X2*X2 + Y2*Y2 + Z2*Z2)
                      XY12 = X1*Y2 - X2*Y1
                      YZ12 = Y1*Z2 - Y2*Z1
                      ZX12 = Z1*X2 - Z2*X1
                      RX12 = RR1*X2 - RR2*X1
                      RY12 = RR1*Y2 - RR2*Y1
                      RZ12 = RR1*Z2 - RR2*Z1

                      RR3  = X3*X3 + Y3*Y3 + Z3*Z3
                      DET2 = X3*YZ12 + Y3*ZX12 + Z3*XY12
                      VP1X  = (RR3*YZ12 + RY12*Z3 - RZ12*Y3)/DET2
                      VP1Y  = (RR3*ZX12 + RZ12*X3 - RX12*Z3)/DET2
                      VP1Z  = (RR3*XY12 + RX12*Y3 - RY12*X3)/DET2

                     !$$--- save the position
                      VVERT(IVERT+1,1) = VP1X*rescal
                      VVERT(IVERT+1,2) = VP1Y*rescal
                      VVERT(IVERT+1,3) = VP1Z*rescal

                     !$$--- get the position of the second vertex on face LF
                      X2  = X3
                      Y2  = Y3
                      Z2  = Z3
                      RR2 = RR3
                      XY12 = X1*Y2 - X2*Y1
                      YZ12 = Y1*Z2 - Y2*Z1
                      ZX12 = Z1*X2 - Z2*X1
                      RX12 = RR1*X2 - RR2*X1
                      RY12 = RR1*Y2 - RR2*Y1
                      RZ12 = RR1*Z2 - RR2*Z1

                      IVERT= IVERT+4
                      I3 = DVERT(IVERT)
                      X3 = SEP(I3,1)
                      Y3 = SEP(I3,2)
                      Z3 = SEP(I3,3)

                      RR3  = X3*X3 + Y3*Y3 + Z3*Z3
                      DET2 = X3*YZ12 + Y3*ZX12 + Z3*XY12
                      VP2X  = (RR3*YZ12 + RY12*Z3 - RZ12*Y3)/DET2
                      VP2Y  = (RR3*ZX12 + RZ12*X3 - RX12*Z3)/DET2
                      VP2Z  = (RR3*XY12 + RX12*Y3 - RY12*X3)/DET2

                     !$$--- save the position
                      VVERT(IVERT-2,1) = VP2X*rescal
                      VVERT(IVERT-2,2) = VP2Y*rescal
                      VVERT(IVERT-2,3) = VP2Z*rescal

                     !$$--- prepare for next vertex
                      X2  = X3
                      Y2  = Y3
                      Z2  = Z3
                      RR2 = RR3
                      XY12 = X1*Y2 - X2*Y1
                      YZ12 = Y1*Z2 - Y2*Z1
                      ZX12 = Z1*X2 - Z2*X1
                      RX12 = RR1*X2 - RR2*X1
                      RY12 = RR1*Y2 - RR2*Y1
                      RZ12 = RR1*Z2 - RR2*Z1

                     !$$---  scan the next neighbores to construct other vertex
                       do IW=1, MXNV4
                          IVERT= IVERT+1
                          I3 = DVERT(IVERT)
                          if(I3 .LE. 0) exit
                          X3 = SEP(I3,1)
                          Y3 = SEP(I3,2)
                          Z3 = SEP(I3,3)
                          RR3  = X3*X3 + Y3*Y3 + Z3*Z3
                          DET2 = X3*YZ12 + Y3*ZX12 + Z3*XY12
                          VP3X  = (RR3*YZ12 + RY12*Z3 - RZ12*Y3)/DET2
                          VP3Y  = (RR3*ZX12 + RZ12*X3 - RX12*Z3)/DET2
                          VP3Z  = (RR3*XY12 + RX12*Y3 - RY12*X3)/DET2

                         !$$--- accumulating volume of tetrahedron
                          VVV = VVV + ABS( (VP1X*(VP2Y*VP3Z - VP2Z*VP3Y) + &
                                            VP1Y*(VP2Z*VP3X - VP2X*VP3Z) + &
                                            VP1Z*(VP2X*VP3Y - VP2Y*VP3X) ) )

                         !$$--- save the position
                          VVERT(IVERT-2,1) = VP3X*rescal
                          VVERT(IVERT-2,2) = VP3Y*rescal
                          VVERT(IVERT-2,3) = VP3Z*rescal

                         !$$--- prepare for next
                          VP2X = VP3X
                          VP2Y = VP3Y
                          VP2Z = VP3Z
                          X2   = X3
                          Y2   = Y3
                          Z2   = Z3
                          RR2  = RR3
                          XY12 = X1*Y2 - X2*Y1
                          YZ12 = Y1*Z2 - Y2*Z1
                          ZX12 = Z1*X2 - Z2*X1
                          RX12 = RR1*X2 - RR2*X1
                          RY12 = RR1*Y2 - RR2*Y1
                          RZ12 = RR1*Z2 - RR2*Z1
                      end do
                  end do !END LOOP FOR LF
                  VOL(ICS)  = VVV*p_VFACT
                else
                  VOL(ICS)  = -1.E0
                end if

              end if
            end do

      return
  END SUBROUTINE Cal_Voronoi_KERNEL
  !*********************************************************************************

  !*********************************************************************************
   SUBROUTINE StartOnDevice_Voronoi_template(IDEV, CELL0, STARTCELL, ENDCELL, KVOIS,INDI,ITYP, XP, DVTS, STAT, VOL, VVERT)
  !***  PURPOSE:   to calculate atomic volume by Delaunay created above.
  !                Called after calling StartOnDevice_Delaunay_template.
  !
  !     INPUT:     IDEV,      the ID of device
  !                CELL0,     the ID of the first cell on the device
  !                STARTCELL, the ID of the first cell on the device
  !                ENDCELL,   the ID of the last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !                DTVS:   the vertices of Delaunay tessellation
  !                STAT:   the stae indicating if a Voronoi volume is opened or neclosed
  !
  !     OUTPUT:    VOL,    the atomic Voronoi volumes
  !                VVERT,  the positions of Voronoi vertice
  !
   use MD_Globle_Variables_GPU
   implicit none
   !---dummy vaiables
       integer::IDEV, CELL0, STARTCELL, ENDCELL
       integer,      device, dimension(:)  ::KVOIS
       integer,      device, dimension(:,:)::INDI
       integer,      device, dimension(:)  ::ITYP
       real(KINDDF), device, dimension(:,:)::XP
       integer,      device, dimension(:), intent(in)::DVTS
       integer,      device, dimension(:), intent(in)::STAT
       real(KINDSF), device, dimension(:)  ::VOL
       real(KINDSF), device, dimension(:,:)::VVERT

   !---Local variables
       integer::BX, BY, NB, NBX, NBY, ERR, IA0, IAS, ENDA, NPRT

   !---Device variables and variables to be used in GPU
       type(dim3) :: blocks
       type(dim3) :: threads

             ERR = cudaSetDevice(IDEV)

             IA0   = hm_IA1th(CELL0)
             IAS   = hm_IA1th(STARTCELL)-IA0
             ENDA  = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
             NPRT  = ENDA -  hm_IA1th(STARTCELL) + 1

             BX  = mp_BLOCKSIZE
             BY  = 1
             NBX = mp_BLOCKDIMX
             NBY = 1

             blocks  = dim3(NBX, NBY, 1)
             threads = dim3(BX, BY, 1)
             !VVERT = 1.E32
             call Cal_Voronoi_KERNEL<<<blocks, threads>>>(dm_NPRT,XP, m_NAPDEV, IA0-1, IAS, NPRT, m_MXNF, &
                                                                   KVOIS, INDI, m_MXNV, DVTS, STAT, VOL, VVERT)

          return
  end subroutine StartOnDevice_Voronoi_template
  !*********************************************************************************

  !*********************************************************************************
   subroutine ResetWorkingArray_template(IDEV, DVTS)
  !***  PURPOSE:   to reset the initial value of DVT and VVERT
  !
  !     INPUT:     IDEV,      the ID of device
  !
  !     OUTPUT:    DVTS,    the Delanuny vertice
  !                VVERT,   the positions of Voronoi vertice
  !
   use MD_Globle_Variables_GPU
   implicit none
   !---dummy vaiables
       integer::IDEV
       integer, device, dimension(:), allocatable::DVTS

   !---Local variables
       integer::ERR, BX, BY, NB, NBX, NBY
       type(dim3):: blocks
       type(dim3):: threads

             ERR = cudaSetDevice(IDEV)

             BX  = 1024
             BY  = 1
             NBX = mp_BLOCKDIMX
             NBY = 1

             blocks  = dim3(NBX, NBY, 1)
             threads = dim3(BX, BY, 1)
             call ResetWorkingArray_KERNEL<<<blocks, threads>>>(m_MXNADEV, m_MXNF, m_MXNV, DVTS)

      return
  end subroutine ResetWorkingArray_template
  !*********************************************************************************


  !*********************************************************************************
   subroutine ResetVoronoiArraye_template(IDEV, VVERT)
  !***  PURPOSE:   to reset the initial value of DVT and VVERT
  !
  !     INPUT:     IDEV,      the ID of device
  !
  !     OUTPUT:    DVTS,    the Delanuny vertice
  !                VVERT,   the positions of Voronoi vertice
  !
   use MD_Globle_Variables_GPU
   implicit none
   !---dummy vaiables
       integer::IDEV
       real(KINDSF), device, dimension(:,:), allocatable::VVERT

   !---Local variables
       integer::ERR, BX, BY, NB, NBX, NBY
       type(dim3):: blocks
       type(dim3):: threads

             ERR = cudaSetDevice(IDEV)

             BX  = 1024
             BY  = 1
             NBX = mp_BLOCKDIMX
             NBY = 1

             blocks  = dim3(NBX, NBY, 1)
             threads = dim3(BX, BY, 1)
             call ResetVoronoiArray_KERNEL<<<blocks, threads>>>(m_MXNADEV, m_MXNF, m_MXNV, VVERT)

      return
  end subroutine ResetVoronoiArraye_template
  !*********************************************************************************


  !**********************************************************************************
  subroutine ResetWorkingArray()
  !***  PURPOSE:  to reset the working space to initial values
  !               for cells from C0 to C1
  !
  !     INPUT:
  !     OUTPUT    dxm_DVTS:
  !
   use MD_Globle_Variables_GPU
   implicit none
       !--- dymmy varibales
       !--- local variables
       integer:: ERR, CURDEV

                ERR = cudaGetDevice(CURDEV)
                if(m_NDEVICE .GE. 1) then
                   call ResetWorkingArray_template(m_DEVICES(1), d1m_DVTS)
                end if

                if(m_NDEVICE .GE. 2) then
                   call ResetWorkingArray_template(m_DEVICES(2), d2m_DVTS)
                end if

                if(m_NDEVICE .GE. 3) then
                   call ResetWorkingArray_template(m_DEVICES(3), d3m_DVTS)
                end if

                if(m_NDEVICE .GE. 4) then
                   call ResetWorkingArray_template(m_DEVICES(4), d4m_DVTS)
                end if

                if(m_NDEVICE .GE. 5) then
                   call ResetWorkingArray_template(m_DEVICES(5), d5m_DVTS)
                end if

                if(m_NDEVICE .GE. 6) then
                   call ResetWorkingArray_template(m_DEVICES(6), d6m_DVTS)
                end if
                ERR = cudaSetDevice(CURDEV)
           return
  end subroutine ResetWorkingArray
  !**********************************************************************************

  !**********************************************************************************
  subroutine ResetVoronoiArray()
  !***  PURPOSE:  to reset the working space to initial values
  !               for cells from C0 to C1
  !
  !     INPUT:
  !     OUTPUT    dxm_DVTS:
  !               dxm_VVTS:
  !
   use MD_Globle_Variables_GPU
   implicit none
       !--- dymmy varibales
       !--- local variables
       integer:: ERR, CURDEV

                ERR = cudaGetDevice(CURDEV)
                if(m_NDEVICE .GE. 1) then
                   call ResetVoronoiArraye_template(m_DEVICES(1), d1m_VVTS)
                end if

                if(m_NDEVICE .GE. 2) then
                   call ResetVoronoiArraye_template(m_DEVICES(2), d2m_VVTS)
                end if

                if(m_NDEVICE .GE. 3) then
                   call ResetVoronoiArraye_template(m_DEVICES(3), d3m_VVTS)
                end if

                if(m_NDEVICE .GE. 4) then
                   call ResetVoronoiArraye_template(m_DEVICES(4), d4m_VVTS)
                end if

                if(m_NDEVICE .GE. 5) then
                   call ResetVoronoiArraye_template(m_DEVICES(5), d5m_VVTS)
                end if

                if(m_NDEVICE .GE. 6) then
                   call ResetVoronoiArraye_template(m_DEVICES(6), d6m_VVTS)
                end if
                ERR = cudaSetDevice(CURDEV)
           return
  end subroutine ResetVoronoiArray
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Delaunay_template(IDEV, STARTCELL, ENDCELL, hDVTS, dDVTS, SHIFT, NA)
  !***  PURPOSE:  the templet of copyout Delaunay vertice from device to host
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL: the first cell
  !                END CELL:  the last cell
  !                dDVTS:     the Delaunay vertice on the device IDEV
  !                SHIFT:     optional, if present, the starting index of hDVTS is SHIFT
  !                           otherwise, the starting index is the index of the first atom in STARTCELL
  !
  !     OUTPUT     hDVTS:     the Delaunay vertice on host
  !                NA:        the number of Delaunay volume actually copied out.
  !
  !
  !
  use MD_Globle_Variables_GPU
  implicit none
  !---DUMMY Variables
      integer::IDEV, STARTCELL, ENDCELL
      integer, device, dimension(:)::dDVTS
      integer, dimension(:)::hDVTS
      integer, optional::SHIFT, NA

  !---Local variables
      integer IERR, STARTA, ENDA, NV

              IERR = cudaSetDevice(IDEV)

              !$$--- the first atom on the device
               STARTA = hm_IA1th(STARTCELL)

              !$$--- the last atom on the device
               ENDA = hm_IA1th(ENDCELL)+hm_NAC(ENDCELL)-1

              !$$--- the number of atoms on the device
               NV = ENDA - STARTA + 1

              !$$--- NOTE
               if(present(SHIFT)) STARTA = SHIFT
               if(present(NA)) NA = NV

               STARTA = (STARTA-1)*m_MXNF*m_MXNV + 1
               NV  = NV*m_MXNF*m_MXNV
               IERR = cudaMemcpyAsync(hDVTS(STARTA), dDVTS(1), NV)
         return
  end subroutine COPYOUT_Delaunay_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_DelaunayStat_template(IDEV, STARTCELL, ENDCELL, hSTAT, dSTAT, SHIFT, NA)
  !***  PURPOSE:  the templet of copyout Delaunay vertice stat from device to host
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL: the first cell
  !                END CELL:  the last cell
  !                dSTAT:     the stat indicating if a Voronoi volume is opened or enclosed
  !                SHIFT:     optional, if present, the starting index of hDVTS is SHIFT
  !                           otherwise, the starting index is the index of the first atom in STARTCELL
  !
  !     OUTPUT     hSTAT:     the state on host
  !                NA:        the number of Delaunay volume actually copied out.
  !
  !
  !
  use MD_Globle_Variables_GPU
  implicit none
  !---DUMMY Variables
      integer::IDEV, STARTCELL, ENDCELL
      integer, device, dimension(:)::dSTAT
      integer, dimension(:)::hSTAT
      integer, optional::SHIFT, NA

  !---Local variables
      integer IERR, STARTA, ENDA, NV

              IERR = cudaSetDevice(IDEV)

              !$$--- the first atom on the device
               STARTA = hm_IA1th(STARTCELL)

              !$$--- the last atom on the device
               ENDA = hm_IA1th(ENDCELL)+hm_NAC(ENDCELL)-1

              !$$--- the number of atoms on the device
               NV = ENDA - STARTA + 1

              !$$--- NOTE
               if(present(SHIFT)) STARTA = SHIFT
               if(present(NA))    NA = NV

               IERR = cudaMemcpyAsync(hSTAT(STARTA), dSTAT(1), NV)

         return
  end subroutine COPYOUT_DelaunayStat_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, hDVTS, SHIFT, TNA)
  !***  PURPOSE:  to copy out Delaunay vertices get in Cal_Delaunay_Vertice_For_Cells
  !               for cells from C0 to C1
  !
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    C0:    the starting cell
  !               C1:    the endding  cell
  !               SHIFT: optional, if present, the starting index of hDVTS is SHIFT
  !                      and the DVTS would copyou to hDVTS contineously.
  !
  !     OUTPUT    hDVTS: the Delaunay vertice on host
  !               TNA:   the actually number of atom that have copy out
  !
   use MD_Globle_Variables_GPU, only:m_NDEVICE, m_DEVICES
   implicit none
   !---dymmy varibales
       integer::C0(*), C1(*)
       integer, dimension(:)::hDVTS
       integer,     optional::SHIFT, TNA
   !---local variables
       integer ERR, CURDEV, OFFSET, NA
   !------------------------------------

       ERR = cudaGetDevice(CURDEV)
      !$$--- Delaunay VTS from devices to host
       if(present(SHIFT) .AND. present(TNA) ) then
          NA  = 0
          OFFSET = SHIFT
          if(m_NDEVICE .GE. 1) then
             call COPYOUT_Delaunay_template(m_DEVICES(1), C0(1), C1(1), hDVTS, d1m_DVTS, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 2) then
             call COPYOUT_Delaunay_template(m_DEVICES(2), C0(2), C1(2), hDVTS, d2m_DVTS, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 3) then
             call COPYOUT_Delaunay_template(m_DEVICES(3), C0(3), C1(3), hDVTS, d3m_DVTS, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 4) then
             call COPYOUT_Delaunay_template(m_DEVICES(4), C0(4), C1(4), hDVTS, d4m_DVTS, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 5) then
             call COPYOUT_Delaunay_template(m_DEVICES(5), C0(5), C1(5), hDVTS, d5m_DVTS, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 6) then
             call COPYOUT_Delaunay_template(m_DEVICES(6), C0(6), C1(6), hDVTS, d6m_DVTS, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          TNA = OFFSET - SHIFT

       else
          if(m_NDEVICE .GE. 1) then
             call COPYOUT_Delaunay_template(m_DEVICES(1), C0(1), C1(1), hDVTS, d1m_DVTS)
          end if

          if(m_NDEVICE .GE. 2) then
             call COPYOUT_Delaunay_template(m_DEVICES(2), C0(2), C1(2), hDVTS, d2m_DVTS)
          end if

          if(m_NDEVICE .GE. 3) then
             call COPYOUT_Delaunay_template(m_DEVICES(3), C0(3), C1(3), hDVTS, d3m_DVTS)
          end if

          if(m_NDEVICE .GE. 4) then
             call COPYOUT_Delaunay_template(m_DEVICES(4), C0(4), C1(4), hDVTS, d4m_DVTS)
          end if

          if(m_NDEVICE .GE. 5) then
             call COPYOUT_Delaunay_template(m_DEVICES(5), C0(5), C1(5), hDVTS, d5m_DVTS)
          end if

          if(m_NDEVICE .GE. 6) then
             call COPYOUT_Delaunay_template(m_DEVICES(6), C0(6), C1(6), hDVTS, d6m_DVTS)
          end if
       end if
       ERR = cudaSetDevice(CURDEV)

      return
  end subroutine COPYOUT_Delaunay_Vertices_For_Cells
  !************************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Delaunay_Stat_For_Cells(C0, C1, hSTAT, SHIFT, TNA)
  !***  PURPOSE:  to copy out Delaunay state get in Cal_Delaunay_Vertice_For_Cells
  !               for cells from C0 to C1
  !
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    C0:    the starting cell
  !               C1:    the endding  cell
  !               SHIFT: optional, if present, the starting index of hDVTS is SHIFT
  !                      and the DVTS would copyou to hDVTS contineously.
  !
  !     OUTPUT    hSTAT: the Delaunay tessellation stat on host
  !               TNA:   the actually number of atom that have copy out
  !
   use MD_Globle_Variables_GPU, only:m_NDEVICE, m_DEVICES
   implicit none
   !---dymmy varibales
       integer::C0(*), C1(*)
       integer, dimension(:)::hSTAT
       integer, optional    ::SHIFT, TNA
   !---local variables
       integer ERR, CURDEV, OFFSET, NA
   !------------------------------------

       ERR = cudaGetDevice(CURDEV)
      !$$--- Delaunay VTS from devices to host
       if(present(SHIFT) .AND. present(TNA) ) then
          NA  = 0
          OFFSET = SHIFT
          if(m_NDEVICE .GE. 1) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(1), C0(1), C1(1), hSTAT, d1m_VSTAT, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 2) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(2), C0(2), C1(2), hSTAT, d2m_VSTAT, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 3) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(3), C0(3), C1(3), hSTAT, d3m_VSTAT, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 4) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(4), C0(4), C1(4), hSTAT, d4m_VSTAT, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 5) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(5), C0(5), C1(5), hSTAT, d5m_VSTAT, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          if(m_NDEVICE .GE. 6) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(6), C0(6), C1(6), hSTAT, d6m_VSTAT, OFFSET, NA)
             OFFSET = OFFSET + NA
          end if

          TNA = OFFSET - SHIFT

       else
          if(m_NDEVICE .GE. 1) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(1), C0(1), C1(1), hSTAT, d1m_VSTAT)
          end if

          if(m_NDEVICE .GE. 2) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(2), C0(2), C1(2), hSTAT, d2m_VSTAT)
          end if

          if(m_NDEVICE .GE. 3) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(3), C0(3), C1(3), hSTAT, d3m_VSTAT)
          end if

          if(m_NDEVICE .GE. 4) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(4), C0(4), C1(4), hSTAT, d4m_VSTAT)
          end if

          if(m_NDEVICE .GE. 5) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(5), C0(5), C1(5), hSTAT, d5m_VSTAT)
          end if

          if(m_NDEVICE .GE. 6) then
             call COPYOUT_DelaunayStat_template(m_DEVICES(6), C0(6), C1(6), hSTAT, d6m_VSTAT)
          end if
       end if
       ERR = cudaSetDevice(CURDEV)

      return
  end subroutine COPYOUT_Delaunay_Stat_For_Cells
  !****************************************************************************************

  !**********************************************************************************
  subroutine OUTPUT_Delaunay_Vertices_Header(hFile)
  !***  PURPOSE:  to output header of a file store Delaunay vertices
  !
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT
   implicit none
   integer, intent(in)::hFile

              !$$--- write out the hearder of the file
              write(hFile, fmt="(A)") "!--- DELAUNAY VERTICES CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
              write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
              write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
              write(hFile, fmt="(A)") '!    '
              write(hFile, fmt="(A)") "!--- Data format following this header:"
              write(hFile, fmt="(A)") "!"
              write(hFile, fmt="(A)") "!--- NATOMS:  the total number of sites."
              write(hFile, fmt="(A)") "!--- * seperator  --------   (repeated for NATOMS times)"
              write(hFile, fmt="(A)") "!--- INDEX, IA0, TYPE, XP"
              write(hFile, fmt="(A)") "!    where: INDEX, the id of the Delaunay site"
              write(hFile, fmt="(A)") "!           IA0,   the id of the center atom (site INDEX) in the original configuration at this site"
              write(hFile, fmt="(A)") "!           TYPE,  the type of the atom IA"
              write(hFile, fmt="(A)") "!           XP,    the position of the atom"
              write(hFile, fmt="(A)") "!--- NFACE: the number of faces of the INDEX Voronoi site"
              write(hFile, fmt="(A)") "!--- NVERT: the number of vertice on constructing the face ----- (repeated for  NFACE times)"
              write(hFile, fmt="(A)") "!--- IN, IA, IA0, TYPE, XP(IA0)- XP(IAC)                   -----  (repeated for NVERT times)"
              write(hFile, fmt="(A)") "!    where: IN,    the index of vertix on the face"
              write(hFile, fmt="(A)") "!           IA,    the id of the atom in the neighbor list of the INDEX Delaunay site"
              write(hFile, fmt="(A)") "!           IA0,   the id of the atom (site IA) in the original configuration at this site"
              write(hFile, fmt="(A)") "!           TYPE,  the type of the atom IA0 (or site IA)"
              write(hFile, fmt="(A)") "!           XP(IA)- XP(IA0),  the relative position between the atom IA and the center atom"
              write(hFile, fmt="(A)") "!--- * seperator  for next sites "
              write(hFile, fmt="(A)")
              write(hFile, fmt="(A,I)")  "&NATOM   ", dm_NPRT
              write(hFile, fmt="(A,I)")  "&NAINBOX ", hm_NAINBOX

      return
  end subroutine OUTPUT_Delaunay_Vertices_Header
  !**********************************************************************************

  !**********************************************************************************
  subroutine OUTPUT_Delaunay_Vertices_For_Cells(hFile, C0, C1, LIST, LATT)
  !***  PURPOSE:  to output Delaunay vertices get in Cal_Delaunay_Vertice_For_Cells
  !               for cells from C0 to C1 to a fileuint. This subroutine is different
  !               from Output_DelaunayVertice(hFile, List, LATT, TDVTS) by that
  !               the output Delaunay vertices are for atoms with indice in the
  !               cells from C0-C1 of the partitioned system created in neighbor-list
  !               calculations.
  !
  !
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    hFile: the unit for output
  !               C0:    the starting cell
  !               C1:    the endding  cell
  !               LIST:  the neighbor-list
  !               LATT:  length unit
  !     OUTPUT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
   !---dymmy varibales
       integer, intent(in)::hFile
       integer::C0(*), C1(*)
       type(NEIGHBOR_LIST),intent(in)::List
       real(KINDDF),intent(in)::LATT
   !---local variables
       integer, dimension(:), allocatable::hDVTS
       integer, dimension(:), allocatable::NVERTS
       integer::NFACES, I, J, K, IA, IA0,INDEX, IP0, IN, IP, IERR, IDEV, NA

   !---
           allocate(hDVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE), NVERTS(m_MXNF), STAT=IERR)
           if(IERR) then
              write(*,fmt="(' MDPSCU Warning: fail to allocate memeory to output Delaunay_Vertices')")
              write(*,fmt="('                 without outputing Delaunay vertices')")
              call ONWARNING(gm_OnWarning)
              return
           end if

           call COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, hDVTS, SHIFT=1, TNA= NA)
           call SynchronizeDevices()

          !$$--- reorde the indices in partitioned system
          !$$
          !$$    NOTE: When the system is large and multiple GPUs
          !$$    are used, it could be that C0-C1 can not cover all cells
          !$$    created in neighborelist calculations. For example,
          !$$    C1(1) < m_ENDCELL(1). Thus the atoms in one call in
          !$$    Delaunay Tessellation is not in continous sequece.
          !$$    The code has been changed to handle this propblem
          !$$    (2015-04-10).
           IP0 = 0
           do IDEV = 1, m_NDEVICE
              INDEX = hm_IA1th(C0(IDEV))-1
              NA = hm_IA1th(C1(IDEV))+hm_NAC(C1(IDEV)) - hm_IA1th(C0(IDEV))

              do I=1, NA
                 NFACES = 0
                 NVERTS = 0
                 INDEX  = INDEX + 1
                 IA0 = hm_GID(INDEX)

                 !get the number of faces and vertices
                 do J=1, m_MXNF
                    IP = IP0+(J-1)*m_MXNV+1
                    IN = hDVTS(IP)
                    if(IN.LE.0) exit

                    NFACES = NFACES + 1
                    do K=1, m_MXNV
                      NVERTS(J) = NVERTS(J) + 1
                      IP = IP + 1
                      IN = hDVTS(IP)
                      if(IN.LE.0) exit
                    end do
                 end do
                 write(hFile, *) "****************************************************************"
                 write(hFile, fmt="(3(I8, 1x), 5(1PE15.6,2x))") INDEX, IA0, hm_ITYP(INDEX),hm_XP(INDEX,1:3)/LATT
                 write(hFile, fmt="(I8, 1x, 5(1PE15.6,2x))") NFACES
                 do J=1, NFACES
                    write(hFile, fmt="(I8, 1x, 5(1PE15.6,2x))") NVERTS(J)
                    IP = IP0+(J-1)*m_MXNV
                    do K=1,NVERTS(J)
                       IP  = IP + 1
                       IN  = hDVTS(IP)
                       IA  = List%INDI(INDEX,IN)
                       IA0 = hm_GID(IA)
                       write(hFile, fmt="(4(I8, 1x), 5(1PE15.6,2x))") &
                                      K, IA, IA0, hm_ITYP(IA), (hm_XP(IA,1:3)- hm_XP(INDEX,1:3))/LATT
                    end do
                 end do
                 IP0 = IP0 + m_MXNF*m_MXNV
              end do
           end do
           deallocate(hDVTS, NVERTS)

      return
  end subroutine OUTPUT_Delaunay_Vertices_For_Cells
  !****************************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Voronoi_template(IDEV, STARTCELL, ENDCELL, hVTS, dVTS, SHIFT, NA)
  !***  PURPOSE:  the templet of copyout Voronoi vertice and volumes from device to host
  !
  !     INPUT:     IDEV,  the ID of device
  !                dVTS,  the VTS on the device IDEV
  !                dVOL,  the atomic volume on the device IDEV
  !                SHIFT: optional, if present, the starting index of hDVTS is SHIFT
  !                       otherwise, the starting index is the index of the first atom in STARTCELL

  !     OUTPUT     hVTS:  the VTS on host
  !                hVOL:  the VOL on host
  !
  use MD_Globle_Variables_GPU
  implicit none
  !---DUMMY Variables
      integer::IDEV, STARTCELL, ENDCELL
      real(KINDSF), device, dimension(:,:)::dVTS
      real(KINDSF), dimension(:,:)::hVTS
      integer,optional::SHIFT, NA

  !---Local variables
      integer IERR, STARTA, ENDA, NV

             IERR = cudaSetDevice(IDEV)

             !$$--- the first atom on the device
              STARTA = hm_IA1th(STARTCELL)

             !$$--- the last atom on the device
              ENDA = hm_IA1th(ENDCELL)+hm_NAC(ENDCELL)-1

             !$$--- the number of atoms on the device
              NV = ENDA - STARTA + 1

             !$$--- NOTE
              if(present(SHIFT)) STARTA = SHIFT
              if(present(NA)) NA = NV

              STARTA = (STARTA-1)*m_MXNF*m_MXNV + 1
              NV     = NV*m_MXNF*m_MXNV
              IERR = cudaMemcpyAsync(hVTS(STARTA,1), dVTS(1,1),NV)
              IERR = cudaMemcpyAsync(hVTS(STARTA,2), dVTS(1,2),NV)
              IERR = cudaMemcpyAsync(hVTS(STARTA,3), dVTS(1,3),NV)

      return
  end subroutine COPYOUT_Voronoi_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Voronoi_Vertices_For_Cells(C0, C1, hVVTS, SHIFT, TNA)
  !***  PURPOSE:  to copy out Voronoi vertices get in Cal_Voronoi_Volume_For_Cells
  !               for cells from C0 to C1
  !
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    C0:  the starting cell
  !               C1:  the endding  cell
  !               SHIFT: optional, if present, the starting index of hVVTS is SHIFT
  !                      and the dVVTS would copyout to hVVTS continuously.
  !
  !     OUTPUT    hVVTS:     the Voronoi vertice on host
  !               hVOLS:     the Voronoi volume on host
  !
   use MD_Globle_Variables_GPU, only:m_NDEVICE, m_DEVICES
   implicit none
   !---dymmy varibales
       integer, intent(in)::C0(*), C1(*)
       real(KINDSF), dimension(:,:)::hVVTS
       integer, optional::SHIFT, TNA
   !---local variables
       integer ERR, CURDEV, OFFSET, NA

   !---

          ERR = cudaGetDevice(CURDEV)
         !$$--- Delaunay VTS from devices to host
          if(present(SHIFT) .AND. present(TNA)) then
             NA = 0
             OFFSET = SHIFT
             if(m_NDEVICE .GE. 1) then
                call COPYOUT_Voronoi_template(m_DEVICES(1), C0(1), C1(1), hVVTS, d1m_VVTS, OFFSET, NA)
                OFFSET = OFFSET + NA
             end if

             if(m_NDEVICE .GE. 2) then
                call COPYOUT_Voronoi_template(m_DEVICES(2), C0(2), C1(2), hVVTS, d2m_VVTS, OFFSET, NA)
                OFFSET = OFFSET + NA
             end if

             if(m_NDEVICE .GE. 3) then
                call COPYOUT_Voronoi_template(m_DEVICES(3), C0(3), C1(3), hVVTS, d3m_VVTS, OFFSET, NA)
                OFFSET = OFFSET + NA
             end if

             if(m_NDEVICE .GE. 4) then
                call COPYOUT_Voronoi_template(m_DEVICES(4), C0(4), C1(4), hVVTS, d4m_VVTS, OFFSET, NA)
                OFFSET = OFFSET + NA
             end if

             if(m_NDEVICE .GE. 5) then
                call COPYOUT_Voronoi_template(m_DEVICES(5), C0(5), C1(5), hVVTS, d5m_VVTS, OFFSET, NA)
                OFFSET = OFFSET + NA
             end if

             if(m_NDEVICE .GE. 6) then
                call COPYOUT_Voronoi_template(m_DEVICES(6), C0(6), C1(6), hVVTS, d6m_VVTS, OFFSET, NA)
                OFFSET = OFFSET + NA
             end if
             TNA = OFFSET-SHIFT

          else
          !$$--- Voronoi VTS from devices to host
             if(m_NDEVICE .GE. 1) then
                call COPYOUT_Voronoi_template(m_DEVICES(1), C0(1), C1(1), hVVTS, d1m_VVTS)
             end if

             if(m_NDEVICE .GE. 2) then
                call COPYOUT_Voronoi_template(m_DEVICES(2), C0(2), C1(2), hVVTS, d2m_VVTS)
             end if

             if(m_NDEVICE .GE. 3) then
                call COPYOUT_Voronoi_template(m_DEVICES(3), C0(3), C1(3), hVVTS, d3m_VVTS)
             end if

             if(m_NDEVICE .GE. 4) then
                call COPYOUT_Voronoi_template(m_DEVICES(4), C0(4), C1(4), hVVTS, d4m_VVTS)
             end if

             if(m_NDEVICE .GE. 5) then
                call COPYOUT_Voronoi_template(m_DEVICES(5), C0(5), C1(5), hVVTS, d5m_VVTS)
             end if

             if(m_NDEVICE .GE. 6) then
                call COPYOUT_Voronoi_template(m_DEVICES(6), C0(6), C1(6), hVVTS, d6m_VVTS)
             end if
          end if
         ERR = cudaSetDevice(CURDEV)

      return
  end subroutine COPYOUT_Voronoi_Vertices_For_Cells
  !**********************************************************************************

 !**********************************************************************************
  subroutine OUTPUT_Voronoi_Vertices_Header(hFile)
  !***  PURPOSE:  to output header of a file store Delaunay vertices
  !
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT
   implicit none
   integer, intent(in)::hFile

              !$$--- write out the hearder of the file
              write(hFile, fmt="(A)") "!--- VORONOI VERTICES CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
              write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
              write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
              write(hFile, fmt="(A)") '!    '
              write(hFile, fmt="(A)") "!--- Data format following this header:"
              write(hFile, fmt="(A)") "!"
              write(hFile, fmt="(A)") "!--- NATOMS:  the total number of sites."
              write(hFile, fmt="(A)") "!--- * seperator  --------   (repeated for NATOMS times)"
              write(hFile, fmt="(A)") "!--- INDEX, IA0, TYPE, XP"
              write(hFile, fmt="(A)") "!    where: INDEX, the id of the Delaunay site"
              write(hFile, fmt="(A)") "!           IA0,   the id of the center atom (site INDEX) in the orgiginal configuration at this site"
              write(hFile, fmt="(A)") "!           TYPE,  the type of the atom IAC"
              write(hFile, fmt="(A)") "!           XP,    the position of the atom"
              write(hFile, fmt="(A)") "!--- NFACE: the number of faces of the INDEX Voronoi site"
              write(hFile, fmt="(A)") "!--- NVERT: the number of vertice on a face -----  (repeated for  NFACE times)"
              write(hFile, fmt="(A)") "!--- IF, IV,  XP(IA)- XP(IA0)              -----  (repeated for NVERT times)"
              write(hFile, fmt="(A)") "!    where: IF,    the index of the face"
              write(hFile, fmt="(A)") "!           IV,    the id of the vertice on the face"
              write(hFile, fmt="(A)") "!           XP(IA)- XP(IA0),  the relative position between the vertex and the center atom"
              write(hFile, fmt="(A)") "!--- * seperator  for next sites "
              write(hFile, fmt="(A)")
              write(hFile, fmt="(A,I5)")"&VORONOI_ATOM"
              write(hFile, fmt="(A,I)") "&NATOM   ", dm_NPRT
              write(hFile, fmt="(A,I)") "&NAINBOX ", hm_NAINBOX


      return
  end subroutine OUTPUT_Voronoi_Vertices_Header

  !**********************************************************************************
  subroutine OUTPUT_Voronoi_Vertices_For_Cells(hFile, C0, C1, LIST, LATT)
  !***  PURPOSE:  to output Voronoi vertices get in Cal_Voronoi_Volume_For_Cells
  !               for cells from C0 to C1 to a file uint. This subroutine is different
  !               from Output_VoronoiVertice(hFile, List, LATT, TDVTS) by that
  !               the output Voronoi vertices are for atoms with indice in the
  !               cells from C0-C1 of the partitioned system created in neighbor-list
  !               calculations.
  !
  !
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    hFile: the unit for output
  !               C0:    the starting cell
  !               C1:    the endding  cell
  !               LIST:  the neighbor-list
  !               LATT:  length unit
  !     OUTPUT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
   !---dymmy varibales
       integer, intent(in)::hFile
       integer::C0(*), C1(*)
       type(NEIGHBOR_LIST),intent(in)::List
       real(KINDDF),intent(in)::LATT
   !---local variables
       real(KINDSF), dimension(:,:), allocatable::hVVTS
       integer, dimension(:), allocatable::hDVTS
       integer, dimension(:), allocatable::NVERTS
       integer::TNA, NFACES, I, J, K, INDEX, IP0, IP, IA, IN, IERR, IDEV, NA
   !---
           allocate(hVVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE,3), &
                    hDVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE), NVERTS(m_MXNF), &
                    STAT =  IERR)
           if(IERR) then
              write(*,fmt="(' MDPSCU Warning: fail to allocate memeory to output Voronoi_Vertices')")
              write(*,fmt="('                 without outputing Voronoi vertices')")
              call ONWARNING(gm_OnWarning)
              return
           end if

           call COPYOUT_Voronoi_Vertices_For_Cells(C0, C1, hVVTS, SHIFT=1, TNA= NA)
           call COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, hDVTS, SHIFT=1, TNA= NA)

           call SynchronizeDevices()
          !$$--- reorde the indices in partitioned system
          !$$    NOTE: When the system is large and multiple GPUs
          !$$    are used, it could be that C0-C1 can not cover all cells
          !$$    created in neighborelist calculations. For example,
          !$$    C1(1) < m_ENDCELL(1). Thus the atoms in one call to
          !$$    Delaunay Tessellation is not in continous sequece.
          !$$    The code has been changed to handle this propblem
          !$$    (2015-04-10).
           IP0 = 0
           do IDEV = 1, m_NDEVICE
              INDEX = hm_IA1th(C0(IDEV))-1
              NA = hm_IA1th(C1(IDEV))+hm_NAC(C1(IDEV)) - hm_IA1th(C0(IDEV))

              do I=1, NA
                 NFACES = 0
                 NVERTS = 0
                 INDEX  = INDEX + 1

                 !get the number of faces and vertices
                 do J=1, m_MXNF
                    IP = IP0+(J-1)*m_MXNV+1
                    if(ANY(hVVTS(IP,1:3) .GE. 1.E31)) then
                       exit
                    end if

                    NFACES = NFACES + 1
                    do K=1, m_MXNV
                       NVERTS(J) = NVERTS(J) + 1
                       IP = IP + 1
                       if(ANY(hVVTS(IP,1:3) .GE. 1.E31)) exit
                    end do
                 end do

                 write(hFile, *) "****************************************************************"
                 write(hFile, fmt="(3(I8, 1x), 5(1PE15.6,2x))") INDEX, hm_GID(INDEX), hm_ITYP(INDEX), hm_XP(INDEX,1:3)/LATT
                 write(hFile, fmt="(I8, 1x, 5(1PE15.6,2x))") NFACES
                 do J=1, NFACES
                   ! we also output the normal of the face
                    IP  = IP0+(J-1)*m_MXNV
                    IN  = hDVTS(IP+1)
                    IA  = List%INDI(INDEX,IN)
                    write(hFile, fmt="(2(I8, 1x), 5(1PE15.6,2x))") NVERTS(J),  hm_ITYP(IA), (hm_XP(IA,1:3)-hm_XP(INDEX,1:3))/LATT
                    do K=1,NVERTS(J)
                       IP  = IP + 1
                       write(hFile, fmt="(I5, 1x, I5, 1x, 5(1PE15.6,2x))") J, K, hVVTS(IP,1:3)/LATT
                    end do
                 end do
                 IP0 = IP0 + m_MXNF*m_MXNV
              end do
           end do
          deallocate(hVVTS, hDVTS, NVERTS)

      return
  end subroutine OUTPUT_Voronoi_Vertices_For_Cells
  !****************************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Voronoi_Volume_template(IDEV, STARTCELL, ENDCELL, hVol, dVol, hSTAT, dSTAT, SHIFT, NA)
  !***  PURPOSE:  the templet of copyout Voronoi vertice and volumes from device to host
  !
  !     INPUT:     IDEV,  the ID of device
  !                dVOL,   the atomic volume on the device IDEV
  !                dStat,  the state (opened or enclosed) of Voronoi volume
  !                SHIFT:  optional, if present, the starting index of hVOL is SHIFT
  !                        otherwise, the starting index is the index of the first atom in STARTCELL
  !
  !     OUTPUT     hVOL:  the VOLS on host
  !                hSTAT:  the STAT on host
  !
  use MD_Globle_Variables_GPU

  implicit none
  !---DUMMY Variables
      integer::IDEV, STARTCELL, ENDCELL
      real(KINDSF), device, dimension(:), intent(in)::dVOL
      real(KINDSF), dimension(:)::hVOL
      integer, device, dimension(:), intent(in)::dSTAT
      integer, dimension(:)::hSTAT
      integer,optional::SHIFT, NA

  !---Local variables
      integer IERR, STARTA, ENDA, NV

  !---
           IERR = cudaSetDevice(IDEV)

           !$$--- the first atom on the device
            STARTA = hm_IA1th(STARTCELL)

           !$$--- the last atom on the device
            ENDA = hm_IA1th(ENDCELL)+hm_NAC(ENDCELL)-1

           !$$--- the number of atoms on the device
            NV   = ENDA - STARTA + 1

           !$$--- NOTE
            if(present(SHIFT)) STARTA = SHIFT
            if(present(NA))    NA     = NV

            IERR = cudaMemcpyAsync(hVOL(STARTA), dVOL(1),NV)
            IERR = cudaMemcpyAsync(hSTAT(STARTA), dSTAT(1),NV)

      return
  end subroutine COPYOUT_Voronoi_Volume_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_Voronoi_Volume_For_Cells(C0, C1, hVOLS, hSTAT, SHIFT, TNA)
  !***  PURPOSE:  to copy out Voronoi volume get in Cal_Voronoi_Volume_For_Cells
  !               for cells from C0 to C1
  !
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    C0:  the starting cell
  !               C1:  the endding  cell
  !               SHIFT: optional, if present, the starting index of hVOLS is SHIFT
  !                      and the dVOLS would copyout to hVOLS contineously.
  !
  !
  !     OUTPUT    hVOLS:     the Voronoi volume on host
  !               hSTAT:     the state of Voronoi volume on host
  !
   use MD_Globle_Variables_GPU
   implicit none
   !---dymmy varibales
    integer, intent(in)::C0(*), C1(*)
    real(KINDSF), dimension(:)::hVOLS
    integer, dimension(:)::hSTAT
    integer, optional::SHIFT, TNA
   !---local variables
    integer ERR, CURDEV, OFFSET, NA

   !---

        ERR = cudaGetDevice(CURDEV)
       !$$--- Voronoi VTS from devices to host
        if(present(SHIFT) .AND. present(TNA)) then
           NA = 0
           OFFSET = SHIFT
           if(m_NDEVICE .GE. 1) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(1), C0(1), C1(1), hVOLS, d1m_VOLS, hSTAT, d1m_VSTAT, OFFSET, NA)
              OFFSET = OFFSET + NA
           end if

           if(m_NDEVICE .GE. 2) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(2), C0(2), C1(2), hVOLS, d2m_VOLS, hSTAT, d2m_VSTAT, OFFSET, NA)
              OFFSET = OFFSET + NA
           end if

           if(m_NDEVICE .GE. 3) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(3), C0(3), C1(3), hVOLS, d3m_VOLS, hSTAT, d3m_VSTAT, OFFSET, NA)
              OFFSET = OFFSET + NA
           end if

           if(m_NDEVICE .GE. 4) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(4), C0(4), C1(4), hVOLS, d4m_VOLS, hSTAT, d4m_VSTAT, OFFSET, NA)
              OFFSET = OFFSET + NA
           end if

           if(m_NDEVICE .GE. 5) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(5), C0(5), C1(5), hVOLS, d5m_VOLS, hSTAT, d5m_VSTAT, OFFSET, NA)
              OFFSET = OFFSET + NA
           end if

           if(m_NDEVICE .GE. 6) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(6), C0(6), C1(6), hVOLS, d6m_VOLS, hSTAT, d6m_VSTAT, OFFSET, NA)
              OFFSET = OFFSET + NA
           end if
           TNA = OFFSET - SHIFT

        else
           if(m_NDEVICE .GE. 1) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(1), C0(1), C1(1), hVOLS, d1m_VOLS, hSTAT, d1m_VSTAT)
           end if

           if(m_NDEVICE .GE. 2) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(2), C0(2), C1(2), hVOLS, d2m_VOLS, hSTAT, d2m_VSTAT)
           end if

           if(m_NDEVICE .GE. 3) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(3), C0(3), C1(3), hVOLS, d3m_VOLS, hSTAT, d3m_VSTAT)
           end if

           if(m_NDEVICE .GE. 4) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(4), C0(4), C1(4), hVOLS, d4m_VOLS, hSTAT, d4m_VSTAT)
           end if

           if(m_NDEVICE .GE. 5) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(5), C0(5), C1(5), hVOLS, d5m_VOLS, hSTAT, d5m_VSTAT)
           end if

           if(m_NDEVICE .GE. 6) then
              call COPYOUT_Voronoi_Volume_template(m_DEVICES(6), C0(6), C1(6), hVOLS, d6m_VOLS, hSTAT, d6m_VSTAT)
           end if
        end if
        ERR = cudaSetDevice(CURDEV)

      return
  end subroutine COPYOUT_Voronoi_Volume_For_Cells
  !**********************************************************************************

  !*********************************************************************************
  subroutine Cal_Starting_Delaunay_Vertice_For_Cells(C0, C1)
  !***  PURPOSE:  to calculate the starting vertice for Delaunay tessellation
  !
  !     INPUT:    C0:  the starting cell
  !               C1:  the endding  cell
  !
  !     OUTPUT:   dxm_DVTS: module array of Delaunay vertices
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       integer, dimension(:)::C0, C1
   !---Local variables
       integer::ERR, CURDEV

         ERR = cudaGetDevice(CURDEV)

        !$$--- caluclate the Delaunay vertices of atoms
         if(m_NDEVICE .GE. 1) then
            call StartOnDevice_Delaunay_template0(m_DEVICES(1), m_STARTCELL(1), C0(1), C1(1),           &
                                             dm_Neighbors%KVOIS(1)%Data, dm_Neighbors%INDI(1)%Data,     &
                                             dm_WorkSpace%ITYP(1)%Data,  dm_WorkSpace%XP(1)%Data,        &
                                             d1m_AMASK, d1m_HMASK,  d1m_DVTS)
         end if

         if(m_NDEVICE .GE. 2) then
            call StartOnDevice_Delaunay_template0(m_DEVICES(2), m_STARTCELL(2), C0(2), C1(2),            &
                                             dm_Neighbors%KVOIS(2)%Data, dm_Neighbors%INDI(2)%Data,      &
                                             dm_WorkSpace%ITYP(2)%Data,  dm_WorkSpace%XP(2)%Data,        &
                                             d2m_AMASK, d2m_HMASK,  d2m_DVTS)
         end if

         if(m_NDEVICE .GE. 3) then
            call StartOnDevice_Delaunay_template0(m_DEVICES(3), m_STARTCELL(3), C0(3), C1(3),            &
                                             dm_Neighbors%KVOIS(3)%Data, dm_Neighbors%INDI(3)%Data,      &
                                             dm_WorkSpace%ITYP(3)%Data,  dm_WorkSpace%XP(3)%Data,        &
                                             d3m_AMASK, d3m_HMASK,  d3m_DVTS)
         end if

         if(m_NDEVICE .GE. 4) then
            call StartOnDevice_Delaunay_template0(m_DEVICES(4), m_STARTCELL(4),  C0(4), C1(4),           &
                                             dm_Neighbors%KVOIS(4)%Data, dm_Neighbors%INDI(4)%Data,      &
                                             dm_WorkSpace%ITYP(4)%Data,  dm_WorkSpace%XP(4)%Data,        &
                                             d4m_AMASK, d4m_HMASK,  d4m_DVTS)
         end if

         if(m_NDEVICE .GE. 5) then
            call StartOnDevice_Delaunay_template0(m_DEVICES(5), m_STARTCELL(5),  C0(5), C1(5),           &
                                             dm_Neighbors%KVOIS(5)%Data, dm_Neighbors%INDI(5)%Data,      &
                                             dm_WorkSpace%ITYP(5)%Data,  dm_WorkSpace%XP(5)%Data,        &
                                             d5m_AMASK, d5m_HMASK,  d5m_DVTS)
         end if

         if(m_NDEVICE .GE. 6) then
            call StartOnDevice_Delaunay_template0(m_DEVICES(6), m_STARTCELL(6), C0(6), C1(6),            &
                                             dm_Neighbors%KVOIS(6)%Data, dm_Neighbors%INDI(6)%Data,      &
                                             dm_WorkSpace%ITYP(6)%Data,  dm_WorkSpace%XP(6)%Data,        &
                                             d6m_AMASK, d6m_HMASK,  d6m_DVTS)
         end if
         ERR = cudaSetDevice(CURDEV)

      return
  end subroutine Cal_Starting_Delaunay_Vertice_For_Cells
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_Delaunay_Vertice_For_Cells(C0, C1)
  !***  PURPOSE:  to calculate the Delaunay vertice from cell C0( >=m_STARTCELL) to C1
  !               (<=m_ENDCELL)
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    C0:  the starting cell
  !               C1:  the endding  cell
  !
  !     OUTPUT:   dxm_DVTS: module array of Delaunay vertices
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       integer, dimension(:)::C0, C1
   !---Local variables
       integer::ERR, CURDEV

          ERR = cudaGetDevice(CURDEV)

         !$$--- caluclate the Delaunay vertices of atoms
          if(m_NDEVICE .GE. 1) then
             if(C1(1).LE. m_ENDCELL(1)) then
                call StartOnDevice_Delaunay_template1(m_DEVICES(1), m_STARTCELL(1), C0(1), C1(1),         &
                                             dm_Neighbors%KVOIS(1)%Data, dm_Neighbors%INDI(1)%Data,       &
                                             dm_WorkSpace%ITYP(1)%Data,  dm_WorkSpace%XP(1)%Data,         &
                                             d1m_AMASK, d1m_HMASK,  d1m_DVTS, d1m_VSTAT, d1m_WVFS, d1m_WFLG)
             end if
          end if

          if(m_NDEVICE .GE. 2) then
             if(C1(2).LE. m_ENDCELL(2)) then
                call StartOnDevice_Delaunay_template1(m_DEVICES(2), m_STARTCELL(2), C0(2), C1(2),         &
                                             dm_Neighbors%KVOIS(2)%Data, dm_Neighbors%INDI(2)%Data,       &
                                             dm_WorkSpace%ITYP(2)%Data,  dm_WorkSpace%XP(2)%Data,         &
                                             d2m_AMASK, d2m_HMASK,  d2m_DVTS, d2m_VSTAT, d2m_WVFS, d2m_WFLG)
              end if
          end if

          if(m_NDEVICE .GE. 3) then
             if(C1(3).LE. m_ENDCELL(3)) then
                call StartOnDevice_Delaunay_template1(m_DEVICES(3), m_STARTCELL(3), C0(3), C1(3),         &
                                             dm_Neighbors%KVOIS(3)%Data, dm_Neighbors%INDI(3)%Data,       &
                                             dm_WorkSpace%ITYP(3)%Data,  dm_WorkSpace%XP(3)%Data,         &
                                             d3m_AMASK, d3m_HMASK,  d3m_DVTS, d3m_VSTAT, d3m_WVFS, d3m_WFLG)
             end if
          end if

          if(m_NDEVICE .GE. 4) then
             if(C1(4).LE. m_ENDCELL(4)) then
                call StartOnDevice_Delaunay_template1(m_DEVICES(4), m_STARTCELL(4),  C0(4), C1(4),        &
                                             dm_Neighbors%KVOIS(4)%Data, dm_Neighbors%INDI(4)%Data,       &
                                             dm_WorkSpace%ITYP(4)%Data,  dm_WorkSpace%XP(4)%Data,         &
                                             d4m_AMASK, d4m_HMASK,  d4m_DVTS, d4m_VSTAT, d4m_WVFS, d4m_WFLG)
             end if
          end if

          if(m_NDEVICE .GE. 5) then
             if(C1(5).LE. m_ENDCELL(5)) then
                call StartOnDevice_Delaunay_template1(m_DEVICES(5), m_STARTCELL(5),  C0(5), C1(5),        &
                                             dm_Neighbors%KVOIS(5)%Data, dm_Neighbors%INDI(5)%Data,       &
                                             dm_WorkSpace%ITYP(5)%Data,  dm_WorkSpace%XP(5)%Data,         &
                                             d5m_AMASK, d5m_HMASK,  d5m_DVTS, d5m_VSTAT, d5m_WVFS, d5m_WFLG)
             end if
          end if

          if(m_NDEVICE .GE. 6) then
             if(C1(6).LE. m_ENDCELL(6)) then
                call StartOnDevice_Delaunay_template1(m_DEVICES(6), m_STARTCELL(6), C0(6), C1(6),         &
                                             dm_Neighbors%KVOIS(6)%Data, dm_Neighbors%INDI(6)%Data,       &
                                             dm_WorkSpace%ITYP(6)%Data,  dm_WorkSpace%XP(6)%Data,         &
                                             d6m_AMASK, d6m_HMASK, d6m_DVTS, d6m_VSTAT, d6m_WVFS, d6m_WFLG)
             end if
          end if
          ERR = cudaSetDevice(CURDEV)

      return
  end subroutine Cal_Delaunay_Vertice_For_Cells
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_Voronoi_Volume_For_Cells(C0, C1)
  !***  PURPOSE:  to calculate the Voronoi volume from cell C0( >=m_STARTCELL) to C1
  !               (<=m_ENDCELL)
  !     NOTE:     C0 > =m_STARTCELL and C1 <= m_ENDCELL gorverned by each
  !               device. Referring to memooty allocation of DVTS and VVTS.
  !
  !     INPUT:    C0:  the starting cell
  !               C1:  the endding  cell
  !
  !     OUTPUT:   dxm_VVTS: module array for Voronoi  vetrices
  !               dxm_VOLS: module array for Voronoi volume
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       integer, dimension(:)::C0, C1
   !---Local variables
       integer::ERR, CURDEV

          ERR = cudaGetDevice(CURDEV)

        !$$--- calculate the VTS of atoms
         if(m_NDEVICE .GE. 1) then
            if(C1(1).LE. m_ENDCELL(1)) then
               call StartOnDevice_Voronoi_template(m_DEVICES(1), m_STARTCELL(1), C0(1), C1(1),            &
                                             dm_Neighbors%KVOIS(1)%Data, dm_Neighbors%INDI(1)%Data,       &
                                             dm_WorkSpace%ITYP(1)%Data,  dm_WorkSpace%XP(1)%Data,         &
                                             d1m_DVTS, d1m_VSTAT, d1m_VOLS, d1m_VVTS)
            end if
         end if

         if(m_NDEVICE .GE. 2) then
            if(C1(2).LE. m_ENDCELL(2)) then
               call StartOnDevice_Voronoi_template(m_DEVICES(2), m_STARTCELL(2), C0(2), C1(2),            &
                                             dm_Neighbors%KVOIS(2)%Data, dm_Neighbors%INDI(2)%Data,       &
                                             dm_WorkSpace%ITYP(2)%Data,  dm_WorkSpace%XP(2)%Data,         &
                                             d2m_DVTS, d2m_VSTAT, d2m_VOLS, d2m_VVTS)
             end if
         end if

         if(m_NDEVICE .GE. 3) then
            if(C1(3).LE. m_ENDCELL(3)) then
               call StartOnDevice_Voronoi_template(m_DEVICES(3), m_STARTCELL(3), C0(3), C1(3),            &
                                             dm_Neighbors%KVOIS(3)%Data, dm_Neighbors%INDI(3)%Data,       &
                                             dm_WorkSpace%ITYP(3)%Data,  dm_WorkSpace%XP(3)%Data,         &
                                             d3m_DVTS, d3m_VSTAT, d3m_VOLS, d3m_VVTS)
            end if
         end if

         if(m_NDEVICE .GE. 4) then
            if(C1(4).LE. m_ENDCELL(4)) then
               call StartOnDevice_Voronoi_template(m_DEVICES(4), m_STARTCELL(4), C0(4), C1(4),            &
                                             dm_Neighbors%KVOIS(4)%Data, dm_Neighbors%INDI(4)%Data,       &
                                             dm_WorkSpace%ITYP(4)%Data,  dm_WorkSpace%XP(4)%Data,         &
                                             d4m_DVTS, d4m_VSTAT, d4m_VOLS, d4m_VVTS)
            end if
         end if

         if(m_NDEVICE .GE. 5) then
            if(C1(5).LE. m_ENDCELL(5)) then
               call StartOnDevice_Voronoi_template(m_DEVICES(5), m_STARTCELL(5), C0(5), C1(5),            &
                                             dm_Neighbors%KVOIS(5)%Data, dm_Neighbors%INDI(5)%Data,       &
                                             dm_WorkSpace%ITYP(5)%Data,  dm_WorkSpace%XP(5)%Data,         &
                                             d5m_DVTS, d5m_VSTAT, d5m_VOLS, d5m_VVTS)
            end if
         end if

         if(m_NDEVICE .GE. 6) then
            if(C1(6).LE. m_ENDCELL(6)) then
               call StartOnDevice_Voronoi_template(m_DEVICES(6), m_STARTCELL(6), C0(6), C1(6),            &
                                             dm_Neighbors%KVOIS(6)%Data, dm_Neighbors%INDI(6)%Data,       &
                                             dm_WorkSpace%ITYP(6)%Data,  dm_WorkSpace%XP(6)%Data,         &
                                             d6m_DVTS, d6m_VSTAT, d6m_VOLS, d6m_VVTS)
            end if
         end if
         ERR = cudaSetDevice(CURDEV)
      return
  end subroutine Cal_Voronoi_Volume_For_Cells
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_Voronoi_Volume(SimBox, CtrlParam, LIST, LATT, hDVTS, hVVTS, hVOLS, hCYST, hSTAT)
  !***  PURPOSE:  to calculate the Voronoi volume of all atoms
  !
  !    INPUT:   SimBox,    the simulation boxs
  !             CtrlParam, the control parameters for simulation
  !             LIST,      optional, the neighbore list calculated by MD_NeighborsList2C_12_GPU.F90
  !                        used here if output of Delaunay vertices is reuired
  !
  !    OUTPUT: hDVTS,      optional, the Delaunay vertices
  !            hVVTS,      optional, the Voronoi vertices
  !            hVOLS,      optional, the Voronoi volume of atoms
  !            hSTAT,      optional, the stat indicating a Voronoi volume is open or close
  !            hCYST,      optional, indicating the cystal signature of atoms, NOTE the size of  hCYST is 2*NPRT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
   !---dummy vaiables
       type(SimMDBox),  intent(in)::SimBox
       type(SimMDCtrl), intent(in)::CtrlParam

   !--- the optional variable
       type(NEIGHBOR_LIST),optional::List
       real(KINDDF), optional::LATT
       integer,      dimension(:),  optional::hDVTS
       real(KINDSF), dimension(:,:),optional::hVVTS
       real(KINDSF), dimension(:),  optional::hVOLS
       integer,      dimension(:),  optional::hSTAT
       integer(8),   dimension(:),  optional::hCYST
  !--- Local variables
       integer::NCELL, C0(m_MXDEVICE), C1(m_MXDEVICE), I, J, NEEDCYST, ERR, NA, IP0, IP, IDEV, INDEX, NV, NBYTE
       integer, dimension(:), allocatable::tDVTS
       integer(8)::I8swap(2)
       byte::Byteswap(16)
       equivalence(I8swap, Byteswap)


          if(m_INITED .eq. 0) then
             call Allocate_WorkingArray_DEV()
          end if

          if(m_INITED_V .eq. 0) then
             call Allocate_VoronoiVolume_DEV()
          end if
          !--- generate the mask
          if(associated(m_pAMaskProc)) then
             call m_pAMaskProc(SimBox, CtrlParam, hm_AMASK)
             call Copyin_AMask()
             hm_NAINBOX = count(hm_AMASK .eq. mp_AMASK_SEL)
          end if

          if(associated(m_pHMaskProc)) then
             call m_pHMaskProc(SimBox, CtrlParam, hm_HMASK)
             call Copyin_HMask()
          end if

          if(present(hCYST)) then
             allocate(tDVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE),STAT=ERR)
             if(ERR) then
                write(*,fmt="(A,A)")        ' MDPSCU Error: fail to allocate memory tDVTS in Cal_Voronoi_Volume '
                write(*,fmt="(A)")          '               Process to be stopped'
               stop
              end if
             NEEDCYST = 1
             !$$---
             NBYTE = size(hCYST)/dm_NPRT
          else
             NEEDCYST = 0
          end if

          !$$--- If pressure effect acounted for, the box size could be changed, we
          !$$    need to recopy the boxsize
           m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
           m_BOXSIZE(1:3)      = SimBox%ZL(1:3)
           m_IFPD = CtrlParam%IFPD

          !$$--- determine the number of cells in one calculation
           C0 = m_STARTCELL
           if(m_MXNADEV.LE.0 .OR. m_MXNADEV .EQ. m_NAPDEV ) then
              C1 = m_ENDCELL
              NCELL = 1
           else
              NCELL = m_MXNADEV/maxval(hm_NAC)
              C1    = C0 + NCELL-1
              C1    = min(C1, m_ENDCELL)
           end if

           do while(.TRUE.)
              !$$--- create Delaunay vertice for cells C0-C1
              call ResetWorkingArray()
              call Cal_Starting_Delaunay_Vertice_For_Cells(C0, C1)
              call Cal_Delaunay_Vertice_For_Cells(C0, C1)

              if(NEEDCYST) then
                 call COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, tDVTS, SHIFT=1, TNA= NA)
              end if

              !$$--- create Voronoi vertice and volume for cells C0-C1
              call ResetVoronoiArray()
              call Cal_Voronoi_Volume_For_Cells(C0, C1)

              !$$--- to create the structure fingerprint if required
              if(NEEDCYST) then
                 IP0 = 0
                 do IDEV = 1, m_NDEVICE
                    INDEX = hm_IA1th(C0(IDEV))-1
                    NA   = hm_IA1th(C1(IDEV))+hm_NAC(C1(IDEV)) - hm_IA1th(C0(IDEV))
                    do I=1, NA
                       INDEX  = INDEX + 1
                       IP0    =  IP0  + 1
                      !--- get the number of faces and vertices
                       IP       =  (IP0-1)*m_MXNF*m_MXNV + 1
                       I8swap = 0
                       do J=1, m_MXNF
                          if(tDVTS(IP).LE.0) exit
                          !$$--- the first byte storeing the number of faces
                          Byteswap(1) = Byteswap(1) + 1

                          !$$--- get the number of vertice on the faces
                          !$$--- NOTE: the first atom is the atom determining the face.
                          !$$          the ID for vertice start from the second index,
                          !$$          and the second is the same as the vertex
                          NV = count(tDVTS(IP:IP+m_MXNV-1).gt.0)-2

                          !$$-- shift NV, Byteswap(2) storing the number of trangle...
                          NV = NV - 1
                          if((NV.ge.2 ).and. (NV .le. 16)) Byteswap(NV) =  Byteswap(NV) + 1
                         !$$--- move to the start verticeo of next face
                          IP = IP + m_MXNV
                       end do
                       hCYST(2*INDEX-1:2*INDEX) =  I8swap(1:2)
                    end do
                 end do
              end if

             !$$--- ouput current tessellation if I/O unit are provided
              if(m_OUTPUTDT  > 0) call OUTPUT_Delaunay_Vertices_For_Cells(m_OUTPUTDT,C0, C1, LIST, LATT)
              if(m_OUTPUTVT  > 0) call OUTPUT_Voronoi_Vertices_For_Cells(m_OUTPUTVT,C0, C1, LIST, LATT)

             !$$--- copyout the tesselation if the arraies are supplied
             !$$    Note: the size of hDVTS here should be NPRT*MXNF*MVNV which could be a large array
              if(present(hDVTS)) then
                 call COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, hDVTS)
              end if

              if(present(hVVTS)) then
                 call COPYOUT_Voronoi_Vertices_For_Cells(C0, C1, hVVTS)
              end if

             !$$--- move to next cells
              C0 = C1+1
              if(ALL(C0 .GE. m_ENDCELL)) exit

              C1 = C0 + NCELL-1
              C1 = min(C1, m_ENDCELL)
           end do

           if(allocated(tDVTS)) then
              deallocate(tDVTS)
           end if

           if(present(hVOLS) .and. present(hSTAT)) then
              hVOLS = 0.D0
              hSTAT = mp_VSTAT_OPENED
              call COPYOUT_Voronoi_Volume_For_Cells(m_STARTCELL, m_ENDCELL, hVOLS, hSTAT)
              call SynchronizeDevices()
           end if

      return
  end subroutine Cal_Voronoi_Volume
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_VoronoiVolume_DEV(Stamp, SimBox, CtrlParam, hVOLS, hCYST, hSTAT)
  !***  PURPOSE:  to calculate Voronoi tesslation
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    OUTPUT: hVOLS, optional, the atomic volume
  !            hSTAT, optional, indicating if the Voronoi volume is enclosed
  !            hCYST, optional, indicating the cystal signature of atoms, NOTE the size of  hCYST is 2*NPRT
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       type(MDRecordStamp), intent(in):: Stamp
       type(SimMDBox),      intent(in):: SimBox
       type(SimMDCtrl),     intent(in):: CtrlParam

       real(KINDSF), dimension(:)::hVOLS
       integer,      dimension(:)::hSTAT
       integer(8),   dimension(:), optional::hCYST

   !---Local variables
       character*256::GFILE
       integer::OUTPUTDT, OUTPUTVT, JOB, ICFG, ERR
       type(NEIGHBOR_LIST)::List


   !----
       !$$--- determine the output files
        OUTPUTDT   = m_OUTPUTDT
        OUTPUTVT   = m_OUTPUTVT
        JOB        = Stamp%ITest
        ICFG       = Stamp%IRec(1)

        if(OUTPUTDT>0 .OR. OUTPUTVT>0) then
           if(OUTPUTDT>0) then
             !$$--- output for Delaunay Tesslation
              call STRCATI(GFILE, m_OUTFILE, "_D_P", m_processid, 4)
              call STRCATI(GFILE, GFILE, "_", JOB, 4)
              call AvailableIOUnit(m_OUTPUTDT)
              open(UNIT=m_OUTPUTDT, file = GFILE, status='unknown')
              write(*,*) "Delaunay vertices to be saved in "//GFILE(1:len_trim(GFILE))
              !$$--- write out the hearder of the file
              call OUTPUT_Delaunay_Vertices_Header(m_OUTPUTDT)
           end if

           if(OUTPUTVT>0) then
             !$$--- output Voronoi Tesslation
              call STRCATI(GFILE, m_OUTFILE, "_V_P", m_processid, 4)
              call STRCATI(GFILE, GFILE, "_", JOB, 4)
              call STRCATI(GFILE, GFILE, ".", ICFG, 4)
              call AvailableIOUnit(m_OUTPUTVT)
              open(UNIT=m_OUTPUTVT, file = GFILE, status='unknown')
              write(*,*) "Voronoi vertices to be saved in "//GFILE(1:len_trim(GFILE))
              !$$--- write out the hearder of the file
              call OUTPUT_Voronoi_Vertices_Header(m_OUTPUTVT)
           end if
           call Copyout_NeighboreList_DEV(List, ORDER=0)
           if(present(hCYST)) then
              call Cal_Voronoi_Volume(SimBox, CtrlParam, LIST=List, LATT=SimBox%RR, hVOLS=hVOLS, hCYST=hCYST, hSTAT=hSTAT)
           else
              call Cal_Voronoi_Volume(SimBox, CtrlParam, LIST=List, LATT=SimBox%RR, hVOLS=hVOLS, hSTAT=hSTAT)
           end if
           !--- without output Vertices
        else
           !the default process
           if(present(hCYST)) then
              call Cal_Voronoi_Volume(SimBox, CtrlParam, hVOLS=hVOLS, hCYST=hCYST, hSTAT=hSTAT)
           else
              call Cal_Voronoi_Volume(SimBox, CtrlParam, hVOLS=hVOLS, hSTAT=hSTAT)
           end if

        end if

       !restore the output control
        if(m_OUTPUTDT>0 .OR. m_OUTPUTVT>0) then
          call Clear_NeighboreList(List)
          if(m_OUTPUTDT .GT. 0)  close(m_OUTPUTDT)
          if(m_OUTPUTVT .GT. 0)  close(m_OUTPUTVT)
        end if
        m_OUTPUTDT  = OUTPUTDT
        m_OUTPUTVT  = OUTPUTVT

      return
  end subroutine Cal_VoronoiVolume_DEV
 !*********************************************************************************

 !**********************************************************************************
  subroutine Output_VoronoiVol(hFile, List, LATT, TVOLS, TCYST, TSTATU)
  !***  DESCRIPTION: to putout the Voronoi volume to an output file.
  !
  !     INPUT:  hFile,  unit of the output file
  !             LIST,   the neighbore list calculated by MD_NeighborsList2C_12_GPU.F90
  !             LATT,   the latice unit used here for unit converting.
  !             TVOLS, TSTATU, TCYST, the volume and state of Voronoi volume and structure fingerprint
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
   !---dummy vaiables
       integer, intent(in)::hFile
       type(NEIGHBOR_LIST),intent(in)::List
       real(KINDDF),intent(in)::LATT
       real(KINDSF), dimension(:)::TVOLS
       integer,      dimension(:)::TSTATU
       integer(8),   dimension(:), optional::TCYST
   !---local variables
       integer::I, J, IA0, CLSITE, NS, NEEDCYST
       character*32::tStr
       real(KINDSF), dimension(:), allocatable::VOLS
       integer,      dimension(:), allocatable::STATU
       integer(8),   dimension(:), allocatable::CYST
       integer(8),   dimension(:), allocatable::SGNA
       integer,      dimension(:), allocatable::HIS

       real(8)::DSGNA
       integer(8)::I8swap(2)
       byte::Byteswap(16)
       equivalence(I8swap, Byteswap)

       real(KINDSF)::TOTALV, UNIT
   !---
      !$$--- Delaunay VTS from devices to host
       allocate(VOLS(size(TVOLS)),STATU(size(TSTATU)))
       if(present(TCYST) ) then
          NEEDCYST = 1
          allocate (CYST(size(TCYST)) )
       else
          NEEDCYST = 0
       end if

      !$$--- reorde the indices in original system
       TOTALV = 0.0
       UNIT   = LATT**3.D0
       CLSITE = 0
       do I=1, dm_NPRT
          IA0 = hm_GID(I)
          VOLS(IA0)  = TVOLS(I)
          STATU(IA0) = TSTATU(I)
          if(STATU(IA0) .EQ. mp_VSTAT_ENCLOSED) then
             TOTALV = TOTALV+TVOLS(I)/UNIT
             CLSITE = CLSITE + 1
          end if

          if(NEEDCYST) CYST(2*IA0-1:2*IA0)  = TCYST(2*I-1:2*I)
       end do

       print *, "TOTALV", TOTALV
       write(hFile, fmt="(A)") "!--- VORONOI VOLUME RESULTS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
       write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
       write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
       write(hFile, fmt="(A)") '!    '
       write(hFile, fmt="(A,I8)")          "!--- Totol sites:                   ", dm_NPRT
       write(hFile, fmt="(A,I8)")          "!--- Sites are no enclosed:         ", dm_NPRT - CLSITE
       write(hFile, fmt="(A,1PE15.6,2x)")  "!--- Totol Voronoi volume (latt^3) :", TOTALV

       !$$--- create the histogram of signature
       if(NEEDCYST) then
          call Create_SignatureHis(TCYST, TSTATU, SGNA, HIS)
          NS = size(HIS)
          write(hFile, fmt="(A,I8)")          "!--- Number of Voronoi signatures  :", NS
          write(hFile, fmt="(A,1PE15.6,2x)")  "!---      histogram of signatures  :"
          write(hFile, fmt="(A,A,A)")         "!---      #ID    signatures(string)              ", &
                                            "   float singature  ", "  count   "
          do I=1, NS
             I8swap(1:2) = SGNA(2*I-1:2*I)
             tstr = ""
             DSGNA = 0
             do J=1, 16
                if(Byteswap(J) .gt.99) Byteswap(J) = 99
                call STRCATI(tstr, tstr, "", int(Byteswap(J)), 2)
                DSGNA =  DSGNA + Byteswap(J)*100.D0**(16-J)
             end do
             write(hFile, fmt="(A,2x, I3, 2x, A, A32, A, 1pE18.12, 2x, I8)")"!---   ", I, ", ", tstr,  ', ', DSGNA, HIS(I)
          end do
        end if

       !$$---
       write(hFile, fmt="(A,1PE15.6,2x)")  "!--- Totol Voronoi volume distribution :"
       write(hFile, fmt="(A)") "&CFGXYZ"
       write(hFile, fmt="(A,1X,I8,1X, A)") "&NATOM     ", dm_NPRT
       write(hFile, fmt="(A,1X,I8, 1X,A)") "&NAINBOX   ", hm_NAINBOX, "!--- the number of atoms in box"
       write(hFile, fmt="(A,1X,3(I4,1X))") "&TYPECOL   ", 1
       write(hFile, fmt="(A,1X,3(I4,1X))") "&XYZCOL    ", 2, 3, 4
       write(hFile, fmt="(A,1X,3(I4,1X))") "&VOLCOL    ", 5
       write(hFile, fmt="(A,1X,3(I4,1X))") "&DENCOL    ", 6
       write(hFile, fmt="(A,1X,3(I4,1X))") "&STATCOL   ", 7

       if(NEEDCYST) then !--- output with signature
         write(hFile, fmt="(A,1X,1(I4,1X),A)") "&FACENUMCOL ",  8, " !-- the number of faces"
         write(hFile, fmt="(A,1X,1(I4,1X),A)") "&STRUCTCOL  ",  9, " !-- the numbers represent the number of polygons starting from trangle"
         write(hFile, fmt="('!',A6,2x, 7(A15,2x), 2x, A)")                                 &
                           "TYPE",  "X(latt.)", "Y(latt.)", "Z(latt.)",                    &
                           "WS VOL(latt^3)", "DEN.(/latt^3)", "ENCLOSE STATU,",            &
                           "FACE NUMBER,", "SIGNATURE ID"

         do I=1, dm_NPRT
            DSGNA = 0
            I8swap(1:2) = CYST(2*I-1:2*I)
            do J=1, NS
               DSGNA =  DSGNA + Byteswap(J)*100.D0**(16-J)
               if(I8swap(1).eq. SGNA(2*J-1) .and. I8swap(2).eq.SGNA(2*J)) then
                  CLSITE = J
               endif
            end do
            write(hFile, fmt="(I8, 1x, 5(1PE15.6,2x), I4, 14x, I4, 11x, I4, 2x, 1PE18.12)") &
                 m_ITYP(I), m_XP(I,1:3)/LATT, VOLS(I)/UNIT, UNIT/VOLS(I),                   &
                 STATU(I), Byteswap(1), CLSITE,  DSGNA ! Byteswap(2:16)
         end do

       else !--- output without signature
         write(hFile, fmt="('!',A6,2x, 7(A15,2x), 2x, A)")                                  &
                           "TYPE",  "X(latt.)", "Y(latt.)", "Z(latt.)",                     &
                           "WS VOL(latt^3)", "DEN.(/latt^3)", "ENCLOSE STATU"
         do I=1, dm_NPRT
            write(hFile, fmt="(I8, 1x, 5(1PE15.6,2x), I4, 14x, I4, 11x, I4, 2x, 1PE18.12)") &
                 m_ITYP(I), m_XP(I,1:3)/LATT, VOLS(I)/UNIT, UNIT/VOLS(I),                   &
                 STATU(I)
         end do
       end if

       deallocate(VOLS, STATU)
       if(allocated(CYST)) deallocate(CYST)
       if(allocated(SGNA)) deallocate(SGNA)
       if(allocated(HIS))  deallocate(HIS)

      return
  end subroutine Output_VoronoiVol
  !**********************************************************************************

  !**********************************************************************************
  subroutine RECORD_VoronoiTessellation(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results of Voronoi tesselation to output file.
  !                  This routine is to interfaced to and MD_EM_TB_ForceTable_SHELL_12_GPU.
  !                  It is assumed the the neighbor-ist routine has be called
  !
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       type(MDRecordStamp) ,        intent(in):: Stamp
       type(SimMDBox),dimension(:), intent(in):: SimBox
       type(SimMDCtrl),             intent(in):: CtrlParam
   !---local variables
       character*256::GFILE, FMT
       equivalence (GFILE, FMT)
       integer::hFile
       type(NEIGHBOR_LIST)::List
       real(KINDSF), dimension(:), allocatable::hVOLS
       integer,      dimension(:), allocatable::hSTAT
       integer(8),   dimension(:), allocatable::hCYST

       character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
       integer::DATE_TIME1(8),DATE_TIME2(8)
       real*4::C1,C2,TIME1, TIME2
       integer::LOOP=1, I, JOB,ITIME, ICFG

  !-------
           if(m_INITED .eq. 0) then
              call Allocate_WorkingArray_DEV()
           end if

           if(m_INITED_V .eq. 0) then
             call Allocate_VoronoiVolume_DEV()
           end if

           call DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
           C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000

           JOB   = Stamp%ITest
           ICFG  = Stamp%IRec(1)

           allocate( hVOLS(dm_NPRT),hSTAT(dm_NPRT) )
           if(m_OUTPUTST .gt. 0) then
              allocate(hCYST(2*dm_NPRT))
              do I=1, LOOP !--- this loop is used for comput time test
              call Cal_VoronoiVolume_DEV(Stamp, SimBox(1), CtrlParam, hVOLS=hVOLS, hCYST=hCYST, hSTAT=hSTAT)
              end do
           else
              do I=1, LOOP !--- this loop is used for comput time test
              call Cal_VoronoiVolume_DEV(Stamp, SimBox(1), CtrlParam, hVOLS=hVOLS, hSTAT=hSTAT)
              end do
           end if

           call DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
           C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
           TIME2 = C2-C1
           write(*,*) "The elapsed time for Voronoi Tessellation(ms):", TIME2

          !$$--- prepare for output
             call Copyout_NeighboreList_DEV(List, ORDER=0)

          !$$--- output atomic volume
             call STRCATI(GFILE, m_OUTFILE, "_VOL_P", m_processid, 4)
             call STRCATI(GFILE, GFILE, "_", JOB, 4)
             call STRCATI(GFILE, GFILE, ".", ICFG, 4)

             call AvailableIOUnit(hFile)
             write(*,*) "Save Atomic volume  to "//GFILE(1:len_trim(GFILE))
             open(UNIT=hFile, file = GFILE, status='unknown')
             if(associated(m_pOutputProc)) then
                call m_pOutputProc(hFile, List, SimBox(1)%RR, hVOLS, hCYST, hSTAT)
             else
                call Output_VoronoiVol(hFile, List, SimBox(1)%RR, hVOLS, hCYST, hSTAT)
             end if
             close(hFile)
             if(allocated(hVOLS)) deallocate(hVOLS)
             if(allocated(hSTAT)) deallocate(hSTAT)
             if(allocated(hCYST)) deallocate(hCYST)
             call Clear_NeighboreList(List)

      return
  end subroutine RECORD_VoronoiTessellation
 !****************************************************************************************

 !****************************************************************************************
  subroutine RECORD_VoronoiTessellation_TOOL(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calculate and output the results of Voronoi tesselation to output file.
  !                  This routine is interfaced to MD_SimBoxArray_ToolShell_12_GPU as an
  !                  analysis tool
  !
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    SEE ALSO:
  !            MD_SimBoxArray_ToolShell_14_GPU.F90
  implicit none
  !---dummy vaiables
       type(MDRecordStamp) ,        intent(in):: Stamp
       type(SimMDBox),dimension(:), intent(in):: SimBox
       type(SimMDCtrl),             intent(in):: CtrlParam
  !---local variables

       if(Stamp%ITime .LT. 0) then
          call Clear_VoronoiTessellation_DEV(SimBox(1), CtrlParam)
          return
       end if

       call RECORD_VoronoiTessellation(Stamp, SimBox, CtrlParam)

      return
  end subroutine RECORD_VoronoiTessellation_TOOL
  !****************************************************************************************

  !****************************************************************************************
  !****************************************************************************************
  !****************************************************************************************
  subroutine Cal_Delaunay_Vertice(SimBox, CtrlParam, List, hNF, hDVTS, hCYTSID, hSTAT, hFile)
  !***  PURPOSE:  to calculate Delaunay and copy out vertice only.
  !               The Delaunay vertice could be used to define an neighbore list
  !               that is a subset of the neighborlist created by MD_NeighborsList2C.
  !               This neighbore list could be further used, for example, in clustering
  !               calculations of atoms or defects.
  !
  !     INPUT:   SimBox,    the simulation boxs
  !              CtrlParam, the control parameters for simulation
  !              List,      the neighbor-list created by NeighboreList calculations
  !              hFile,     optional, the I/O unit for output intermediate result of Delaunay vertices
  !
  !     OUTPUT:  hNF,       the number of faces of Delaunay volume
  !              hDVTS,     the Delaunay vertice,
  !              hCYSID,    optional, the siguration of the crystal structure, integer *
  !                         the 8 bytes store the number of plygons from 3 to 10
  !              hSTAT,     optional, the state indicting if a Voronoi volume is enclosed or opened.
  !
  !     NOTE:    the indice of atoms are the system that have been ordered by call to
  !              NeighboreList calculations
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
   !---dummy vaiables
       type(SimMDBox),      intent(in)::SimBox
       type(SimMDCtrl),     intent(in)::CtrlParam
       type(NEIGHBOR_LIST), intent(in)::List
       integer,optional               ::hFile

       integer,    dimension(:),           intent(out)::hNF
       integer,    dimension(:,:),         intent(out)::hDVTS
       integer(8), optional, dimension(:), intent(out)::hCYTSID
       integer,    optional, dimension(:), intent(out)::hSTAT

  !--- Local variables
       integer::NCELL, C0(m_MXDEVICE), C1(m_MXDEVICE), I, J, IP, IP0, INDEX, NA, hIO, IDEV, NV, ERR
       integer, dimension(:), allocatable::tDVTS
       integer(8)::I8swap(2)
       byte::Byteswap(16)
       equivalence(I8swap, Byteswap)
       character*256::GFILE


          if(m_INITED .eq. 0) then
             m_MXNF = min(m_MXNF, CtrlParam%NB_MXNBS)
             m_MXNV = min(m_MXNV, CtrlParam%NB_MXNBS)
             m_MXKVOIS = min(CtrlParam%NB_MXNBS, mp_MXKVOIS)
             call Allocate_WorkingArray_DEV()
          end if

          !$$--- Note: m_MXNADEV = min(mp_MXNADEV, m_NAPDEV)
          !$$          m_MXNADEV*m_NDEVICE could be smaller than NPRT
          allocate(tDVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE),STAT=ERR)
          if(ERR) then
             write(*,fmt="(A,A)")        ' MDPSCU Error: fail to allocate memory in Cal_Delaunay_Vertice '
             write(*,fmt="(A)")          '               Process to be stopped'
             stop
          end if

          !$$--- If pressure effect acounted for, the box size could be changed, we
          !$$    need to recopy the boxsize
           m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
           m_BOXSIZE(1:3)      = SimBox%ZL(1:3)
           m_IFPD = CtrlParam%IFPD
           hIO = 0
           if(present(hFile)) hIO = hFile

          !$$--- determine the number of cells in one calculation
           C0 = m_STARTCELL
           if(m_MXNADEV.LE.0 .OR. m_MXNADEV .EQ. m_NAPDEV ) then
              C1 = m_ENDCELL
              NCELL = 1
           else
              NCELL = m_MXNADEV/maxval(hm_NAC)
              C1    = C0 + NCELL-1
              C1    = min(C1, m_ENDCELL)
           end if

          !$$    NOTE: When the system is large and multiple GPUs
          !$$    are used, it could be that C0-C1 can not cover all cells
          !$$    created in neighborelist calculations. For example,
          !$$    C1(1) < m_ENDCELL(1). Thus the atoms in one call to
          !$$    Delaunay Tessellation is not in continous sequece.
          !$$    In other word, size of tDVTS could be smaller than
          !$$    NPRT*m_MXNF*m_MXNV, we can only copyout tDVS for part
          !$$    of atom using SHIFT copying.
          !$$    The code has been changed to handle this propblem
          !$$    (2015-04-10).
          !--- generate the mask
          if(associated(m_pAMaskProc)) then
             call m_pAMaskProc(SimBox, CtrlParam, hm_AMASK)
             call Copyin_Amask()
             hm_NAINBOX = count(hm_AMASK .eq. mp_AMASK_SEL)
          end if
          if(associated(m_pHMaskProc)) then
             call m_pHMaskProc(SimBox, CtrlParam, hm_HMASK)
             call CopyIn_HMask()
          end if

           hNF =0
           do while(.TRUE.)
              !$$--- create Delaunay vertice for cells C0-C1
              call ResetWorkingArray()
              call Cal_Starting_Delaunay_Vertice_For_Cells(C0, C1)
              call Cal_Delaunay_Vertice_For_Cells(C0, C1)
              call COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, tDVTS, SHIFT=1, TNA= NA)

              if(hIO .gt. 0) then
                 call OUTPUT_Delaunay_Vertices_For_Cells(hIO, C0, C1, List, LATT=SimBox%RR)
              else
                 call SynchronizeDevices()
              end if

              IP0 = 0
              do IDEV = 1, m_NDEVICE
                 INDEX = hm_IA1th(C0(IDEV)) - 1
                 NA    = hm_IA1th(C1(IDEV)) + hm_NAC(C1(IDEV)) - hm_IA1th(C0(IDEV))
                 do I=1, NA
                    INDEX  = INDEX + 1
                    IP0    =  IP0  + 1
                    !--- get the number of faces and vertices
                    IP       =  (IP0-1)*m_MXNF*m_MXNV + 1
                    I8swap = 0
                    do J=1, m_MXNF
                       if(tDVTS(IP).le.0) exit

                       hNF(INDEX)      = hNF(INDEX) + 1
                       hDVTS(INDEX, J) = List%INDI(INDEX,tDVTS(IP))
                       if(present(hCYTSID)) then
                          !$$-- the first byte storeing the number of faces
                          I8swap(1) = I8swap(1) + 1
                          !$$--- get the number of vertice on the faces
                          !$$--  NOTE: the first atom is the atom determine the face.
                          !$$          the ID for vertice start from the second index,
                          !$$          and the second is the same as the vertex
                          NV = count(tDVTS(IP:IP+m_MXNV-1).gt.0)-2
                          !$$-- shift NV, Byteswap(2) storing the number of trangle
                          NV = NV - 1
                          if((NV.ge.2 ).and. (NV .le. 16)) &
                             Byteswap(NV) =  Byteswap(NV) + 1
                       end if
                       !$$--- move to the start verticeo of next face
                       IP = IP + m_MXNV
                    end do
                    if(present(hCYTSID)) hCYTSID(2*INDEX-1:2*INDEX) =  I8swap(1:2)
                 end do

              end do
             !$$--- move to next cells
              C0 = C1+1
              if(ALL(C0 .GE. m_ENDCELL)) exit

              C1 = C0 + NCELL-1
              C1 = min(C1, m_ENDCELL)
           end do
           if(present(hSTAT) )then
              hSTAT = mp_VSTAT_OPENED
              call COPYOUT_Delaunay_Stat_For_Cells(m_STARTCELL, m_ENDCELL, hSTAT)
           end if

           deallocate(tDVTS)
      return
  end subroutine Cal_Delaunay_Vertice
  !*********************************************************************************

  !*********************************************************************************
  subroutine Cal_Delaunay_NeighborList(SimBox, CtrlParam, List, Order)
  !***  PURPOSE:  to convert the neighbor-list define by cutoff distance to the
  !               the neighbor-list defined by Delaunay vertices.
  !
  !     INPUT:   SimBox,    the simulation boxs
  !              CtrlParam, the control parameters for simulation
  !              List,      the neighbor-list define by cutoff distance, created by NeighboreList calculations
  !              Order,     the flag indicating if the atom indice should be redordered to
  !                         the original system
  !
  !     OUTPUT:  List,      the neighbor-list define by Delaunay vertices
  !
  !     NOTE:    the indice of atoms are the system that have been ordered by call to
  !              NeighboreList calculations
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       type(SimMDBox),       intent(in)   ::SimBox
       type(SimMDCtrl),      intent(in)   ::CtrlParam
       type(NEIGHBOR_LIST),  intent(inout)::List
       integer, optional                       ::Order
  !--- Local variables
       integer,    dimension(:),  allocatable  ::DNF
       integer,    dimension(:,:),allocatable  ::DVTS
       integer::I, Reorder

          allocate(DNF(dm_NPRT), DVTS(dm_NPRT, m_MXNF))
          call Copyout_NeighboreList_DEV(List, ORDER=0)

          call Cal_Delaunay_Vertice(SimBox, CtrlParam, List, DNF, DVTS)
          call SynchronizeDevices()
          call Clear_NeighboreList(List)

          List%mxKVOIS =  m_MXNF
          allocate(List%KVOIS(dm_NPRT),List%INDI(dm_NPRT,m_MXNF))
          Reorder = 0
          if(present(Order) ) Reorder = Order

          if(Reorder) then
             do I=1, dm_NPRT
                List%KVOIS(hm_GID(I))          = DNF(I)
                List%INDI(hm_GID(I), 1:m_MXNF) = hm_GID(DVTS(I,1:m_MXNF))
             end do
          else
             List%KVOIS(1:dm_NPRT)            = DNF(1:dm_NPRT)
             List%INDI(1:dm_NPRT,1:m_MXNF)    = DVTS(1:dm_NPRT,1:m_MXNF)
          end if
      return
  end subroutine Cal_Delaunay_NeighborList
  !*********************************************************************************

  !**********************************************************************************
  subroutine Create_SignatureHis(TCYST, TSTATU, SGNA, HIS)
  !***  DESCRIPTION: to reate a histogram of Voronoi histogram
  !
  !     INPUT: TSTATU, TCYST, the state of Voronoi volume and structure fingerprint
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList
   implicit none
   !---dummy vaiables
       integer,     dimension(:)::TSTATU
       integer(8),  dimension(:)::TCYST
       integer(8),  dimension(:), allocatable::SGNA
       integer, dimension(:), allocatable::HIS
   !---local variables
       integer::I, J, IA0, CLSITE, NPRT, NS
       character*32::tStr
       integer(8), dimension(:), allocatable::CYSTID
       integer(8)::I8swap(2)
       byte::Byteswap(16)
       equivalence(I8swap, Byteswap)

       real(KINDSF)::TOTALV, UNIT, FLAG
   !---
      !$$--- Determine how many kind of singature we have
       allocate(CYSTID(size(TCYST)) )
       NPRT = size(TSTATU)

       !$$--- the first kind of signature id "opened"
       NS = 1
       Byteswap    = 0
       Byteswap(1) = 1
       CYSTID(1:2) = I8swap(1:2)

       do I=1, NPRT
          if(TSTATU(I) .eq.  mp_VSTAT_OPENED) cycle
          I8swap(1:2) = TCYST(I*2-1:I*2)
          !$$--- to find out if this signature exist
          FLAG = 0
          do J=1, NS
             if(I8swap(1).eq.CYSTID(2*J-1) .and. I8swap(2).eq.CYSTID(2*J) ) then
                FLAG = 1
                exit
             end if
          end do
          if(FLAG .ge. 1) cycle
          NS = NS + 1
          CYSTID(2*NS-1:2*NS) = I8swap(1:2)
       end do

       allocate(HIS(NS), SGNA(2*NS))
       SGNA(1:2*NS) = CYSTID(1:2*NS)
       deallocate(CYSTID)

       !$$--- create the histogram
       HIS = 0
       do I=1, NPRT
          if(TSTATU(I) .eq.  mp_VSTAT_OPENED) then
             HIS(1) = HIS(1) + 1
             cycle
          end if

          !$$--- to find out if this signature exist
          do J=2, NS
             if(SGNA(2*J-1) .eq. TCYST(I*2-1) .and. SGNA(2*J) .eq. TCYST(I*2) ) then
                HIS(J) = HIS(J) + 1
                exit
             end if
          end do
       end do

      return
  end subroutine Create_SignatureHis
  !****************************************************************************************

  !****************************************************************************************
  subroutine Cal_Voronoi_Volume_For_MaskAtoms(SimBox, CtrlParam, List, NAtom, AtomId, Dvts, Vvts, Vols, Stat)
  !***  PURPOSE:  to calculate the Voronoi volumes for marked atoms
  !
  !    INPUT:   SimBox,    the simulation boxs
  !             CtrlParam, the control parameters for simulation
  !             LIST,      optional, the neighbore list calculated by MD_NeighborsList2C_12_GPU.F90
  !                        used here if output of Delaunay vertices is reuired
  !
  !    OUTPUT: NAtom,      the number of atoms involved in calculation
  !            AtomId,     the indice of the atoms
  !            DVTS,       the Delaunay vertices
  !            VVTS,       the Voronoi vertices
  !            VOLS,       the Voronoi volume of atoms
  !            STAT,       the stat indicating a Voronoi volume is open or close
  !
   use MD_NeighborsList
   use MD_Globle_Variables_GPU
   implicit none
   !---dummy vaiables
       type(SimMDBox),              intent(in) ::SimBox
       type(SimMDCtrl),             intent(in) ::CtrlParam
       type(NEIGHBOR_LIST),         intent(in) ::List
       integer,                     intent(out)::NAtom
       integer,      dimension(:),  intent(out)::AtomID
       integer,      dimension(:),  intent(out)::Dvts
       real(KINDSF), dimension(:,:),intent(out)::Vvts
       real(KINDSF), dimension(:),  intent(out)::Vols
       integer,      dimension(:),  intent(out)::Stat
  !--- Local variables
       integer::NCELL, C0(m_MXDEVICE), C1(m_MXDEVICE), I, J0, J1, ERR, NA, IP0, IDEV, NV, IW
       integer,      dimension(:),   allocatable ::TDVTS
       real(KINDSF), dimension(:,:), allocatable ::TVVTS
       real(KINDSF), dimension(:),   allocatable ::TVOLS
       integer,      dimension(:),   allocatable ::TSTAT

 !-----
          allocate(TDVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE),TVVTS(m_MXNADEV*m_MXNF*m_MXNV*m_NDEVICE,3), &
                   TVOLS(m_MXNADEV*m_NDEVICE), TSTAT(m_MXNADEV*m_NDEVICE), STAT=ERR)
          if(ERR) then
              write(*,fmt="(A, I3)") "MDPSCU Error: Fail to allocate working memory in Cal_Voronoi_Volume_For_MaskAtoms"
              write(*,fmt="(A)")     "              Process to be stopped."
              stop
              return
          end if

          !$$--- If pressure effect acounted for, the box size could be changed, we
          !$$    need to recopy the boxsize
           m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
           m_BOXSIZE(1:3)      = SimBox%ZL(1:3)
           m_IFPD = CtrlParam%IFPD

          !$$--- determine the number of cells in one calculation
           C0 = m_STARTCELL
           if(m_MXNADEV.LE.0 .OR. m_MXNADEV .EQ. m_NAPDEV ) then
              C1 = m_ENDCELL
              NCELL = 1
           else
              NCELL = m_MXNADEV/maxval(hm_NAC)
              C1    = C0 + NCELL-1
              C1    = min(C1, m_ENDCELL)
           end if

           IW = 0
           do while(.TRUE.)
              !$$--- create Delaunay vertice for cells C0-C1
              call ResetWorkingArray()
              call Cal_Starting_Delaunay_Vertice_For_Cells(C0, C1)
              call Cal_Delaunay_Vertice_For_Cells(C0, C1)

              !$$--- create Voronoi vertice and volume for cells C0-C1
              call ResetVoronoiArray()
              call Cal_Voronoi_Volume_For_Cells(C0, C1)

             !$$--- copyout the tesselation for cell C0 - C1
              TVOLS = 0.D0
              TSTAT = mp_VSTAT_OPENED
              call COPYOUT_Delaunay_Vertices_For_Cells(C0, C1, TDVTS, SHIFT=1, TNA= NA)
              call COPYOUT_Voronoi_Vertices_For_Cells (C0, C1, TVVTS, SHIFT=1, TNA= NA)
              call COPYOUT_Voronoi_Volume_For_Cells   (C0, C1, TVOLS, TSTAT, SHIFT=1, TNA= NA)
              !$$---
                 IP0 = 0
                 do IDEV = 1, m_NDEVICE
                    do I = hm_IA1th(C0(IDEV)), hm_IA1th(C1(IDEV))+hm_NAC(C1(IDEV))-1
                       IP0    =  IP0  + 1
                       if(hm_AMASK(I) .eq. mp_AMASK_UNSEL) then
                          cycle
                       end if
                      !--- copy the faces and vertices to output array
                       IW         = IW + 1
                       AtomID(IW) = I
                       Vols(IW)   = TVOLS(IP0)
                       Stat(IW)   = TSTAT(IP0)
                       J1         = (IW-1)*m_MXNF+1
                       do J0 = (IP0-1)*m_MXNF*m_MXNV + 1, IP0*m_MXNF*m_MXNV, m_MXNV
                          if(TDVTS(J0) .gt. 0) then
                             Dvts(J1) = List%INDI(I, TDVTS(J0))
                          else
                             Dvts(J1) = 0
                          end if
                          J1 = J1 + 1
                       end do
                       Vvts((IW-1)*m_MXNF*m_MXNV+1:IW*m_MXNF*m_MXNV, 1:3) = &
                                      TVVTS((IP0-1)*m_MXNF*m_MXNV + 1:IP0*m_MXNF*m_MXNV, 1:3)
                    end do
                 end do

             !$$--- move to next cells
              C0 = C1+1
              if(all(C0 .GE. m_ENDCELL)) exit

              C1 = C0 + NCELL-1
              C1 = min(C1, m_ENDCELL)
           end do
           NAtom = IW
          deallocate(TDVTS,TVVTS, TVOLS, TSTAT, STAT=ERR)

      return
  end subroutine Cal_Voronoi_Volume_For_MaskAtoms
 !*********************************************************************************

  !******************************************************************************
  subroutine Cal_Voronoi_ClusterEnvelope(List, NAtom, AtomId, Dvts, Vvts, Vols, Stat, Nc, Na)
  !***  PURPOSE:  to calculate the envelope of clusters constructed by results
  !               generated by calling:
  !               Cal_Voronoi_Volume_For_MaskAtoms
  !
  !               the connected interfaces between connected internal atoms will be
  !               filter out,
  !    INPUT:  List,       the current neighbor list
  !            NAtom,      the number of atoms
  !            AtomId,     the indice of the atoms
  !            DVTS,       the Delaunay vertices
  !            VVTS,       the Voronoi vertices
  !            VOLS,       the Voronoi volume of atoms
  !            STAT,       the stat indicating a Voronoi volume is open or close
  !    OUTPUT:
  !            DVTS,       the PAIRs of Delaunay vertices define the facets of clusters
  !            VVTS,       the Voronoi vertices of facets
  !            VOLS,       the volume of clusters
  !            STAT,       the stat indicating a Voronoi volume is open or close
  !            NC,         the number of cluster envelopes
  !            NA,         the number of atoms in the cluster envelopes
  !
  !    NOTE:  The definition of output DVTS is different from input DVTS.
  !           The  output DVTS define the pair of Delaunay vertices that defines the facet direction.
  !           The  size of DVTS should be at lest NAtom*m_MXNF*2
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT
   use MD_NeighborsList
   implicit none
   !---dummy vaiables
       type(NEIGHBOR_LIST),         intent(in)   ::List
       integer,                     intent(in)   ::NAtom
       integer,      dimension(:),  intent(in)   ::AtomID
       integer,      dimension(:),  intent(inout)::Dvts
       real(KINDSF), dimension(:,:),intent(inout)::Vvts
       real(KINDSF), dimension(:),  intent(inout)::Vols
       integer,      dimension(:),  intent(inout)::Stat
       integer,                     intent(out)  ::Nc
       integer,      dimension(:),  intent(out)  ::Na

  !--- Local variables
       integer::I, J, K, IP0, IP1, STEP, IS, IERR, OFLAG, NA0, IFACE
       real(KINDSF)::VOL0

       integer,      dimension(:),  allocatable::SQNC0, I0, VISIT
       integer,      dimension(:),  allocatable::TDVTS
       real(KINDSF), dimension(:,:),allocatable::TVVTS
       real(KINDSF), dimension(:),  allocatable::TVOL
       integer,      dimension(:),  allocatable::TSTAT


      !$$--- first, we cluatering the atoms
            allocate(SQNC0(NAtom*2),  VISIT(NAtom), I0(dm_NPRT), TSTAT(NAtom), TVOL(NAtom), &
                     TDVTS(NAtom*m_MXNF*2), TVVTS(NAtom*m_MXNF*m_MXNV,3), STAT=IERR )
            if(IERR) then
              write(*,fmt="(A, I3)") "MDPSCU Error: Fail to allocate working memory in Cal_Voronoi_ClusterEnvelope on host"
              write(*,fmt="(A)")     "              Process to be stopped."
              stop
            end if

            TDVTS = 0
            TVVTS = 1.E31
            VISIT = 0
            I0    = 0
            do I=1, NAtom
               VISIT(I)       = 1
               I0(AtomId(I))  = I
            end do

            SQNC0  = 0
            IP0    = 0
            NC     = 0
            IFACE  = 0
            do I=1, NAtom
               if(VISIT(I) .le. 0) cycle

               IP0        = IP0 + 1
               SQNC0(IP0) = I
               IS         = IP0
               NA0        = 1
               VOL0       = 0.0
               OFLAG      = 0

               !$$--- set atom I has been visited
                VISIT(I) = 0
                if(Stat(I) .eq. mp_VSTAT_OPENED) then
                   OFLAG = 1
                else
                   VOL0  = VOL0 + Vols(I)
                end if

                do while(IP0 .le. NAtom)
                   if(SQNC0(IP0) .eq. 0) exit

                   STEP = 0
                   do J= (SQNC0(IP0)-1)*m_MXNF+1, SQNC0(IP0)*m_MXNF
                      if( Dvts(J) .le. 0) exit

                      K = I0(Dvts(J))
                      !$$--- check if Dvts(J) is not in the list of selected atoms
                      !$$    if not, then copy the DVT to TDVT
                      if(K .le. 0) then
                         IFACE         = IFACE + 1
                         TDVTS((IFACE-1)*C_ITWO+1)  = AtomId(SQNC0(IP0))
                         TDVTS((IFACE-1)*C_ITWO+2)  = Dvts(J)
                         TVVTS((IFACE-1)*m_MXNV+1:IFACE*m_MXNV, 1:3) = Vvts((J-1)*m_MXNV+1:J*m_MXNV,1:3)
                         cycle
                      end if

                      if(VISIT(K) .le. 0) then
                         cycle
                      end if

                      STEP           = STEP + 1
                      SQNC0(IS+STEP) = K
                      VISIT(K)       = 0
                      if(Stat(K) .eq. mp_VSTAT_OPENED) then
                         OFLAG = 1
                      else
                         VOL0  = VOL0 + Vols(K)
                      end if

                   end do
                   IP0 = IP0 + 1
                   IS  = IS  + STEP
                   NA0 = NA0 + STEP
                end do
                !$$--- mark for the end of cluster
                !TDVTS(IFACE*2+1) = 0
                !TDVTS(IFACE*2+2) = 0

                IFACE     = IFACE + 1
                IP0       = IS
                Nc        = Nc+1
                Na(Nc)    = NA0
                TSTAT(Nc) = OFLAG
                TVOL(Nc)  = VOL0
            end do
      !$$--- second, copy
            Dvts                      = 0
            Dvts(1:IFACE*2)           = TDVTS(1:IFACE*2)
            Vvts(1:IFACE*m_MXNV, 1:3) = TVVTS(1:IFACE*m_MXNV, 1:3)
            Vols(1:Nc)                = TVOL(1:Nc)
            Stat(1:Nc)                = TSTAT(1:Nc)
            deallocate(SQNC0, I0, VISIT, TSTAT, TVOL, TDVTS, TVVTS, STAT=IERR )
      return
  end subroutine Cal_Voronoi_ClusterEnvelope
  !******************************************************************************

  !******************************************************************************
  subroutine Cal_Voronoi_ClusterEnvelope_Type(SimBox, CtrlParam, CombType, NAtom, AtomId, Dvts, Vvts, Vols, Stat, Nc, Na)
  !***  PURPOSE:  to calculate the envelope of clusters constructed by given types of atoms
  !
  !    INPUT:
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !            CombType,   the combined type of atoms
  !
  !    OUTPUT:
  !            NAtom,      the number of atoms
  !            AtomId,     the indice of the atoms
  !            DVTS,       the Delaunay vertices define the cluster facets
  !            VVTS,       the Voronoi vertices of cluster facets
  !            VOLS,       the Voronoi volume of cluster
  !            STAT,       the stat indicating a Voronoi volume is open or close
  !            NC,         the number of clusters
  !            NA,         the number of atomis in the clusters
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT, hm_ITYP
   use MD_NeighborsList
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       type(SimMDBox),              intent(in)  ::SimBox
       type(SimMDCtrl),             intent(in)  ::CtrlParam
       integer,                     intent(in)  ::CombType
       integer,                     intent(out) ::NAtom
       integer,      dimension(:),  allocatable, intent(out) ::AtomID
       integer,      dimension(:),  allocatable, intent(out) ::Dvts
       real(KINDSF), dimension(:,:),allocatable, intent(out) ::Vvts
       real(KINDSF), dimension(:),  allocatable, intent(out) ::Vols
       integer,      dimension(:),  allocatable, intent(out) ::Stat
       integer,                                  intent(out) ::Nc
       integer,      dimension(:),  allocatable, intent(out) ::Na
  !--- Local variables
       integer::I, ERR
       integer, dimension(:), allocatable::ATOMSEL
       type(NEIGHBOR_LIST)::LIST

            !$$---  we allocat the working array
            if(m_INITED .eq. 0) then
              call Allocate_WorkingArray_DEV()
            end if

            if(m_INITED_V .eq. 0) then
              call Allocate_VoronoiVolume_DEV()
            end if


            !$$--- first, we cluatering the atoms
            call Copyout_NeighboreList_DEV(LIST, ORDER=0)
            allocate(ATOMSEL(dm_NPRT))
            NAtom   = 0
            ATOMSEL = 0
            do I=1, dm_NPRT
               if(iand(2**(hm_ITYP(I)-1), CombType) .gt. 0) then
                  NAtom = NAtom + 1
                  ATOMSEL(NAtom) = I
               end if
            end do
            if(allocated(AtomID)) deallocate(AtomID)
            if(allocated(Dvts))   deallocate(Dvts)
            if(allocated(Vvts))   deallocate(Vvts)
            if(allocated(Vols))   deallocate(Vols)
            if(allocated(Stat))   deallocate(Stat)
            if(allocated(Na))     deallocate(Na)

            allocate(AtomID(NAtom), Dvts(NAtom*m_MXNF*C_ITWO), Vvts(NAtom*m_MXNF*m_MXNV,3), Vols(NAtom), Stat(NAtom), Na(NAtom), STAT=ERR)
            if(ERR) then
               write(*,fmt="(A, I3)") "MDPSCU Error: Fail to allocate memory in Cal_Voronoi_ClusterEnvelope_Type"
               write(*,fmt="(A)")     "                 Process to be stopped."
               stop
            end if

            !--- disable the external atom selector
            call SetAtomSelector()
            !--- create marks for selected atoms
            call CreateMaskForAtoms(SimBox, CtrlParam, ATOMSEL)
            call Cal_Voronoi_Volume_For_MaskAtoms(SimBox, CtrlParam, List, NAtom, AtomId, Dvts, Vvts, Vols, Stat)
            call Cal_Voronoi_ClusterEnvelope(List, NAtom, AtomId, Dvts, Vvts, Vols, Stat, Nc, Na)

            deallocate(ATOMSEL)
            call Clear_NeighboreList(List)
      return
  end subroutine Cal_Voronoi_ClusterEnvelope_Type
  !******************************************************************************

  !******************************************************************************
  subroutine Cal_Voronoi_ClusterEnvelope_Sqn(SimBox, CtrlParam, Sqn, NAtom, AtomId, Dvts, Vvts, Vols, Stat, Nc, Na)
  !***  PURPOSE:  to calculate the envelope of clusters constructed by cluster squence of atoms.
  !               the cluster squence could be created by, for example, PropClusterCommon_GPU.
  !
  !    INPUT:
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !            Sqn,       the cluster squence of atoms
  !
  !    OUTPUT:
  !            NAtom,      the number of atoms
  !            AtomId,     the indice of the atoms
  !            DVTS,       the Delaunay vertices define the cluster facets
  !            VVTS,       the Voronoi vertices of cluster facets
  !            VOLS,       the Voronoi volume of cluster
  !            STAT,       the stat indicating a Voronoi volume is open or close
  !            NC,         the number of clusters
  !            NA,         the number of atomis in the clusters
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT, hm_ITYP
   use MD_NeighborsList
   use MD_NeighborsList_GPU
   implicit none
   !---dummy vaiables
       type(SimMDBox),              intent(in)  ::SimBox
       type(SimMDCtrl),             intent(in)  ::CtrlParam
       integer, dimension(:),       intent(in)  ::Sqn
       integer,                     intent(out) ::NAtom
       integer,      dimension(:),  allocatable, intent(out) ::AtomID
       integer,      dimension(:),  allocatable, intent(out) ::Dvts
       real(KINDSF), dimension(:,:),allocatable, intent(out) ::Vvts
       real(KINDSF), dimension(:),  allocatable, intent(out) ::Vols
       integer,      dimension(:),  allocatable, intent(out) ::Stat
       integer,                                  intent(out) ::Nc
       integer,      dimension(:),  allocatable, intent(out) ::Na
  !--- Local variables
       integer::I, ERR
       integer, dimension(:), allocatable::ATOMSEL
       type(NEIGHBOR_LIST)::LIST

            !$$---  we allocat the working array
            if(m_INITED .eq. 0) then
              call Allocate_WorkingArray_DEV()
            end if

            if(m_INITED_V .eq. 0) then
              call Allocate_VoronoiVolume_DEV()
            end if

            !$$--- first, we cluatering the atoms
            call Copyout_NeighboreList_DEV(LIST, ORDER=0)
            allocate(ATOMSEL(dm_NPRT))
            NAtom   = 0
            ATOMSEL = 0
            do I=1, size(Sqn)
               if(Sqn(I) .gt. 0) then
                  NAtom = NAtom + 1
                  ATOMSEL(NAtom) = Sqn(I)
               end if
            end do
            if(allocated(AtomID)) deallocate(AtomID)
            if(allocated(Dvts))   deallocate(Dvts)
            if(allocated(Vvts))   deallocate(Vvts)
            if(allocated(Vols))   deallocate(Vols)
            if(allocated(Stat))   deallocate(Stat)
            if(allocated(Na))     deallocate(Na)

            allocate(AtomID(NAtom), Dvts(NAtom*m_MXNF*C_ITWO), Vvts(NAtom*m_MXNF*m_MXNV,3), Vols(NAtom), Stat(NAtom), Na(NAtom), STAT=ERR)
            if(ERR) then
               write(*,fmt="(A, I3)") "MDPSCU Error: Fail to allocate memory in Cal_Voronoi_ClusterEnvelope_Type"
               write(*,fmt="(A)")     "                 Process to be stopped."
               stop
            end if

            !--- disable the external atom selector
            call SetAtomSelector()
            !--- create marks for selected atoms
            call CreateMaskForAtoms(SimBox, CtrlParam, ATOMSEL)
            call Cal_Voronoi_Volume_For_MaskAtoms(SimBox, CtrlParam, List, NAtom, AtomId, Dvts, Vvts, Vols, Stat)
            call Cal_Voronoi_ClusterEnvelope(List, NAtom, AtomId, Dvts, Vvts, Vols, Stat, Nc, Na)

            deallocate(ATOMSEL)
            call Clear_NeighboreList(List)
      return
  end subroutine Cal_Voronoi_ClusterEnvelope_Sqn
  !******************************************************************************

  !******************************************************************************
  subroutine Record_Voronoi_ClusterEnvelope_Type(Fname, Stamp, SimBox, CtrlParam, CombType)
  !***  PURPOSE:  to calculate and output the envelope of clusters constructed by given types of atoms
  !
  !    INPUT:
  !            Fname,     the file name of ouput the envelope
  !            Stamp,     the stamp of the configuration
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !            CombType,  the combined type of atoms
  !
  !    OUTPUT:
  !
   implicit none
   !---dummy vaiables
       character*(*),       intent(in)::Fname
       type(SimMDBox),      intent(in)::SimBox
       type(SimMDCtrl),     intent(in)::CtrlParam
       type(MDRecordStamp), intent(in)::Stamp
       integer,             intent(in)::CombType
  !--- Local variables
       integer::NATOM, NC
       integer,      dimension(:),   allocatable::ATOMID, DVTS, NA
       real(KINDSF), dimension(:,:), allocatable::VVTS
       real(KINDSF), dimension(:),   allocatable::VOLS
       integer,      dimension(:),   allocatable::STAT

      !$$--- first, we cluatering the atoms
           if(m_INITED .eq. 0) then
              call Allocate_WorkingArray_DEV()
           end if

           if(m_INITED_V .eq. 0) then
              call Allocate_VoronoiVolume_DEV()
           end if

           call Cal_Voronoi_ClusterEnvelope_Type(SimBox, CtrlParam, CombType, NATOM, ATOMID, DVTS, VVTS, VOLS, STAT, NC, NA)
           call Output_Voronoi_ClusterEnvelope(Fname, Stamp, SimBox, CtrlParam, DVTS, VVTS, VOLS, STAT, NA)

           if(allocated(ATOMID)) deallocate(ATOMID)
           if(allocated(DVTS))   deallocate(DVTS)
           if(allocated(VVTS))   deallocate(VVTS)
           if(allocated(VOLS))   deallocate(VOLS)
           if(allocated(STAT))   deallocate(STAT)
           if(allocated(NA))     deallocate(NA)

      return
  end subroutine Record_Voronoi_ClusterEnvelope_Type
  !******************************************************************************

  !******************************************************************************
  subroutine Record_Voronoi_ClusterEnvelope_Sqn(Fname, Stamp, SimBox, CtrlParam, Sqn)
  !***  PURPOSE:  to calculate and output the envelope of clusters constructed by given
  !               cluster swquence of atoms
  !    INPUT:
  !            Fname,     the file name of ouput the envelope
  !            Stamp,     the stamp of the configuration
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !            Sqn,       the combined type of atoms
  !
  !    OUTPUT:
  !
   implicit none
   !---dummy vaiables
       character*(*),        intent(in)::Fname
       type(SimMDBox),       intent(in)::SimBox
       type(SimMDCtrl),      intent(in)::CtrlParam
       type(MDRecordStamp) , intent(in)::Stamp
       integer,dimension(:), intent(in)::Sqn
  !--- Local variables
       integer::NATOM, NC
       integer,      dimension(:),   allocatable::ATOMID, DVTS, NA
       real(KINDSF), dimension(:,:), allocatable::VVTS
       real(KINDSF), dimension(:),   allocatable::VOLS
       integer,      dimension(:),   allocatable::STAT

      !$$--- first, we cluatering the atoms
           if(m_INITED .eq. 0) then
              call Allocate_WorkingArray_DEV()
           end if

           if(m_INITED_V .eq. 0) then
              call Allocate_VoronoiVolume_DEV()
           end if

           call Cal_Voronoi_ClusterEnvelope_Sqn(SimBox, CtrlParam, Sqn, NATOM, ATOMID, DVTS, VVTS, VOLS, STAT, NC, NA)
           call Output_Voronoi_ClusterEnvelope(Fname, Stamp, SimBox, CtrlParam, DVTS, VVTS, VOLS, STAT, NA)

           if(allocated(ATOMID)) deallocate(ATOMID)
           if(allocated(DVTS))   deallocate(DVTS)
           if(allocated(VVTS))   deallocate(VVTS)
           if(allocated(VOLS))   deallocate(VOLS)
           if(allocated(STAT))   deallocate(STAT)
           if(allocated(NA))     deallocate(NA)

      return
  end subroutine Record_Voronoi_ClusterEnvelope_Sqn
  !******************************************************************************

  !******************************************************************************
  subroutine Output_Voronoi_ClusterEnvelope(Fname, Stamp, SimBox, CtrlParam, Dvts, Vvts, Vols, Stat, Na)
  !***  PURPOSE:  to calculate and output the envelope of clusters constructed by given types of atoms
  !
  !    INPUT:
  !            hFile    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !            CombType,   the combined type of atoms
  !
  !    OUTPUT:
  !
   use MD_Globle_Variables_GPU, only:dm_NPRT, hm_ITYP, hm_GID, hm_XP
   implicit none
   !---dummy vaiables
       character*(*),      intent(in)::Fname
       type(MDRecordStamp),intent(in)::Stamp
       type(SimMDBox),     intent(in)::SimBox
       type(SimMDCtrl),    intent(in)::CtrlParam

       integer,      dimension(:),   intent(in)::Dvts
       real(KINDSF), dimension(:,:), intent(in)::Vvts
       real(KINDSF), dimension(:),   intent(in)::Vols
       integer,      dimension(:),   intent(in)::Stat
       integer,      dimension(:),   intent(in)::Na
  !--- Local variables
       integer, dimension(:), allocatable::IBFLAG1, IBFLAG2
       integer::I, J, K, IP0, NF, NV, NPRT0, NB, IB, NC, IFACE, INDB1, GIDB1, INDB2, GIDB2, hFile
       real(KINDDF)::LATT, LATT3
       character*256::GFILE



               NPRT0 = SimBox%NPRT
               NB    = dm_NPRT/NPRT0
               !--- to determine the IBFLAG of cluster in boxes
               allocate(IBFLAG1(NB), IBFLAG2(NB))
               IBFLAG1 = 0
               IBFLAG2 = -1
               do I=1, size(DVTS)/2
                  if( DVTS((I-1)*2+1) .gt. 0) then
                     IB = ( DVTS((I-1)*2+1)-1 )/NPRT0 + 1
                     if(IBFLAG1(IB) .le. 0) then
                        IBFLAG1(IB) = I
                        IBFLAG2(IB) = I+1
                     else
                        IBFLAG2(IB) = I+1
                     end if
                  end if
               end do

               LATT  = SimBox%RR
               LATT3 = LATT*LATT*LATT
               do IB=1, NB
                  call STRCATI(GFILE, Fname, "P", 0, 4)
                  call STRCATI(GFILE, GFILE, "_", (Stamp%ITest-1)*NB+IB, 4)
                  call STRCATI(GFILE, GFILE, ".", Stamp%IRec(1), 4)
                  call AvailableIOUnit(hFile)
                  open(UNIT=hFile, file = GFILE)
                  write(*,fmt="(A, I6, A)") " Voronoi clusters to be output for box ",(Stamp%ITest-1)*NB+IB, "-> "//GFILE(1:len_trim(GFILE))
                  !--- to extract the number of clusters
                  NC = 0
                  do I = IBFLAG1(IB), IBFLAG2(IB)
                     if( DVTS((I-1)*2+1) .le. 0) then
                         NC = NC + 1
                     end if
                  end do

                  !$$--- write out the hearder of the file
                  write(hFile, fmt="(A)") "!--- VORONOI FACETS FOR CLSTERS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
                  write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
                  write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
                  write(hFile, fmt="(A)") '!    '
                  write(hFile, fmt="(A)") "!--- SYMBOL EXPLAINATION:"
                  write(hFile, fmt="(A)") "!"
                  write(hFile, fmt="(A)") "!--- &VORONOI_CLUS:  the format statement"
                  write(hFile, fmt="(A)") "!--- &NCLUS:         the total number of clusters in the box."
                  write(hFile, fmt="(A)") "!--- &CLUS:          the id of clusters (repeated for NCLUS times)"
                  write(hFile, fmt="(A)") "!--- &NA:            the number of atoms in the cluster"
                  write(hFile, fmt="(A)") "!--- &VOL:           the volume of the cluster in lattice unit"
                  write(hFile, fmt="(A)") "!--- &DEN:           the atom number density of the cluster in lattice unit"
                  write(hFile, fmt="(A)") "!--- &NFACE:         the number of facets of the Voronoi volume"
                  write(hFile, fmt="(A)") "!--- &IFACE:         INDEX, ID0, GID0, TYP0, POS0(1:3)"
                  write(hFile, fmt="(A)") "!--- &IFACE:                ID1, GID1, TYP1, POS1(1:3)"
                  write(hFile, fmt="(A)") "!                    where: INDEX, the id of the facet of the  Voronoi volume"
                  write(hFile, fmt="(A)") "!                           ID0,   the atom id of in-cluster atom in the partitioned box"
                  write(hFile, fmt="(A)") "!                           GID0,  the atom id of in-cluster atom in the original box"
                  write(hFile, fmt="(A)") "!                           TYP0,  the atom type of in-cluster atom"
                  write(hFile, fmt="(A)") "!                           TYP0,  the atom type of in-cluster atom"
                  write(hFile, fmt="(A)") "!                           POS0,  the position of the atom"
                  write(hFile, fmt="(A)") "!                           ID1,   the atom id of out-cluster atom in the partitioned box"
                  write(hFile, fmt="(A)") "!                           GID1,  the atom id of out-cluster atom in the original box"
                  write(hFile, fmt="(A)") "!                           TYP1,  the atom type of out-cluster atom"
                  write(hFile, fmt="(A)") "!                           POS1,  the position of the atom"
                  write(hFile, fmt="(A)") "!--- &NVERT:         the number of vertice on a facet"
                  write(hFile, fmt="(A)") "!                    POS (repeated for NVERT times)"
                  write(hFile, fmt="(A)") "!                    where: POS,   the coordinates of the vertices relative to the in-cluster atom ID0"
                  write(hFile, fmt="(A)")

                  write(hFile, fmt="(A,I5)")        "&VORONOI_CLUS"
                  write(hFile, fmt="(A,I5)")        "&NCLUS         ", NC

                  IP0 = IBFLAG1(IB)
                  do J=1, NC
                     !--- to extract the number of faces in the cluster
                     NF    = 0
                     do I = IP0, IBFLAG2(IB)
                        if( DVTS((I-1)*2+1) .eq. 0) exit
                        NF = NF + 1
                     end do
                     write(hFile, fmt="(A, I5, A)")   "&CLUS ******* #",J, " *******************************"
                     write(hFile, fmt="(A, I7, A)")   "&NA            ",NA(J)
                     write(hFile, fmt="(A, 1PE13.5)") "&VOL(latt^3)   ",VOLS(J)/LATT3
                     write(hFile, fmt="(A, 1PE13.5)") "&DEN(/latt^3)  ",dble(NA(J))/VOLS(J)*LATT3
                     write(hFile, fmt="(A, I7, A)")   "&NFACE         ",NF

                     !---
                     IFACE = IP0 -1
                     do I = 1, NF
                        IFACE = IFACE + 1
                        INDB1 = Dvts((IFACE-1)*2 + 1)
                        GIDB1 = hm_GID(INDB1)
                        INDB2 = Dvts((IFACE-1)*2 + 2)
                        GIDB2 = hm_GID(INDB2)
                        !--- to extracte the number of vertice on the face
                        NV = 0
                        do K = 1, m_MXNV
                           if(any(Vvts((IFACE-1)*m_MXNF + K,1:3) .ge. 1.E30)) exit
                           NV = NV + 1
                        end do

                          write(hFile, fmt="(A11,I6,A1, 2(1x,I8, 1x, I8, 1x, I4, 3(1x, 1PE15.6),A)) )")    &
                                       "  &IFACE # ", I, ",",                                              &
                                      INDB1 - ((INDB1-1)/NPRT0)*NPRT0, GIDB1 - ((GIDB1-1)/NPRT0)*NPRT0,    &
                                      hm_ITYP(INDB1), hm_XP(INDB1,1:3)/LATT
                          write(hFile, fmt="(18x, 2(1x,I8, 1x, I8, 1x, I4, 3(1x, 1PE15.6),A)) )")          &
                                      INDB2 - ((INDB2-1)/NPRT0)*NPRT0, GIDB2 - ((GIDB2-1)/NPRT0)*NPRT0,    &
                                      hm_ITYP(INDB2), hm_XP(INDB2,1:3)/LATT
                          write(hFile, fmt="(A,I6,A, 2(1x,I8, 1x, I8, 1x, I4, 3(1x, 1PE15.6),A)) )")       &
                                       "    &NVERT ", NV
                          do K=1, NV
                             write(hFile, fmt="(A, 5(1PE15.6,2x))") "           ", Vvts((IFACE-1)*m_MXNF+K,1:3)/LATT
                          end do !--- end loop for vertice
                     end do !--- end loop for faces
                     IP0 = IP0 + NF
                  end do !--- end loop for clusters
                  close(hFile)

            end do  !end loop for IB
            deallocate(IBFLAG1, IBFLAG2)

      return
  end subroutine Output_Voronoi_ClusterEnvelope
  !******************************************************************************

  !*****************************************************************************
  subroutine SetSelOuputProc(OUTPUTPROC)
  !***  DESCRIPTION: to set the user defined output process
  !
  !     INPUT: pOutPROC,  the subroutine provided by user
  !
  implicit none
  !--- interface to the external routine -------------------

   interface
        SUBROUTINE OUTPUTPROC(hFile, List, LATT, TVOLS, TCYST, TSTATU)
          use MD_Globle_Variables_GPU
          use MD_NeighborsList
          implicit none
          !---dummy vaiables
          integer, intent(in)::hFile
          type(NEIGHBOR_LIST),intent(in)::List
          real(KINDDF),intent(in)::LATT
          real(KINDSF), dimension(:)::TVOLS
          integer,      dimension(:)::TSTATU
          integer(8),   dimension(:)::TCYST

       END SUBROUTINE OUTPUTPROC
  end interface
  !--- END INTERFACE --------------------------------
           m_pOutputProc=>OUTPUTPROC
           return
   end subroutine SetSelOuputProc
  !*********************************************************************************
  end module VoronoiTessellation_GPU

