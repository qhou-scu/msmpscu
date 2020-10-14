  module AtomicVirialStress_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                 To calculate the local virial stress on atoms. The volume of the atoms are
  !                 determined by Voronoi volume of the atoms. The Voronoi volume are obtaine
  !                 by module VoronoiTessllationM_GPU. If an atom has no enclosed Voronoi
  !                 volume, its pressure is undefined.
  !    DEPENDENCE:  _______________________________________________________________________________________
  !                 MD_NeighborsList2C_12_GPU.F90
  !                 MD_Globle_Variables_12_GPU.F90
  !                 MD_TYPEDEF_FORCELIB_GPU
  !                 VoronoiTessellationM_GPU.F90
  !
  !    REFERENCE:   _______________________________________________________________________________________
  !                 HOU Qing, Nov, 2013
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use CUDAFOR
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_ForceTable

  !-----------------------------------------------
    implicit none
      integer, private::m_processid=0
      !--- the id of  filename get from SimMDCtrl for I/O.
      !    Refering the definition of f_others in SimMDCtrl.
      character(len=10),parameter, private::mp_FTAGI="&AUXF_AVSI"
      character(len=10),parameter, private::mp_FTAGO="&AUXF_AVSO"
      character(len=256),  private::m_OUTFILE =""             ! filename of output data

      !    as private needs also defined here.
      integer, private::m_INITED = 0

      !$$--- parameters controlling calculations and output
      !$$
      integer, private::m_Output_AtomicStress = 0                      ! >1, calculate and output atomic pressure
      integer, private::m_Output_CellStress   = 0                      ! >1, calculate and output pressure
                                                                       !     in cells, defined by m_cellsize
      real(KINDDF), private::m_CellSize(3) = 1.D0
      integer, private::m_Output_SphereStress = 0                      !     calculate and output pressure in
                                                                       !     a sphere of given radiu and center
      real(KINDDF), private::m_SphereCenter(3) = 0.D0
      real(KINDDF), private::m_SphereRmax      = 5.D0
      integer,      private::m_NumShells       = 10

      integer, private::m_Output_TypeStress = 0                        !     calculate and output pressure
                                                                       !     of atoms of types
      integer, private::m_INCTYP(mp_MXGROUP) = 0                       !     type of atoms with their partial
                                                                       !     pressure to be output
      integer, private::m_hFileTypeStress = 0

      !$$--- accumulating number of calculation and
      !$$      accumulating stress
      integer::m_Output_AverageAtomicStress  = 0                         !time step for output average stress
      integer, private::m_numAccum1
      real(KINDDF),dimension(:,:), allocatable, private::m_AVP , m_AKP   !atomic stess

      integer::m_Output_AverageCellStress  = 0                           !time step for output average stress
      integer, private::m_numAccum2
      real(KINDDF),dimension(:,:), allocatable, private::m_CVP , m_CKP   !atomic stess in cell

      integer::m_Output_AverageSphereStress  = 0                         !time step for output average stress
      integer, private::m_numAccum3
      real(KINDDF),dimension(:,:), allocatable, private::m_SVP , m_SKP   !average atomic stess in spheres
      real(KINDDF),dimension(:,:), allocatable, private::m_SDEN          !average atomic density in spheres
      real(KINDDF),dimension(:), allocatable, private::m_SVOL            !average 'real' volum of sphere

      !$$--- atomic pressure for atoms on devices
      !$$--- the atomic pressure - virial part
      real(KINDDF), device, dimension(:,:),  allocatable::d1m_AVP
      real(KINDDF), device, dimension(:,:),  allocatable::d2m_AVP
      real(KINDDF), device, dimension(:,:),  allocatable::d3m_AVP
      real(KINDDF), device, dimension(:,:),  allocatable::d4m_AVP
      real(KINDDF), device, dimension(:,:),  allocatable::d5m_AVP
      real(KINDDF), device, dimension(:,:),  allocatable::d6m_AVP
      !$$--- the atomic pressure - kinetic part
      real(KINDDF), device, dimension(:,:),  allocatable::d1m_AKP
      real(KINDDF), device, dimension(:,:),  allocatable::d2m_AKP
      real(KINDDF), device, dimension(:,:),  allocatable::d3m_AKP
      real(KINDDF), device, dimension(:,:),  allocatable::d4m_AKP
      real(KINDDF), device, dimension(:,:),  allocatable::d5m_AKP
      real(KINDDF), device, dimension(:,:),  allocatable::d6m_AKP

      !$$--- the mass of atoms, needed in calculating kinetic energy
      !$$    NOTE: the type of atoms could be larger than the number
      !$$          of potential types
      real(KINDDF), device, dimension(:), allocatable,private::d1m_CM
      real(KINDDF), device, dimension(:), allocatable,private::d2m_CM
      real(KINDDF), device, dimension(:), allocatable,private::d3m_CM
      real(KINDDF), device, dimension(:), allocatable,private::d4m_CM
      real(KINDDF), device, dimension(:), allocatable,private::d5m_CM
      real(KINDDF), device, dimension(:), allocatable,private::d6m_CM


      private::Clear_AtomicVirialStress_template,  &
               Allocate_LocalPressure_template,    &
               Allocate_WorkingArray_DEV,          &
               LoadControlParameters,              &
               PrintControlParameters

  contains

  !**********************************************************************************
  subroutine Clear_AtomicVirialStress_template(IDEV, dAVP, dAKP, dCM)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                initialization
  !
  !     INPUT:     IDEV,  the index of device on which the memory to be allocated
  !
  !     OUPUT:
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDDF), device, dimension(:,:), allocatable::dAVP, dAKP
           real(KINDDF), device, dimension(:), allocatable::dCM

      !--- Local vairables
           integer::ERR

              ERR = cudaSetDevice(IDEV)
              if(allocated(dAVP)) then
                 deallocate(dAVP, STAT=ERR)
                 if(ERR) goto 100
              end if
              if(allocated(dAKP)) then
                 deallocate(dAKP, STAT=ERR)
                 if(ERR) goto 100
              end if
              if(allocated(dCM)) then
                 deallocate(dCM, STAT=ERR)
                 if(ERR) goto 100
              end if
              return

   100       write(*,*) "MDPSCU WARNING: fail to deallocate memory on device in LocalVirialPressure module", IDEV
             call ONWARNING(gm_OnWarning)
             return
  end subroutine Clear_AtomicVirialStress_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_AtomicVirialStress_DEV(SimBox,CtrlParam)
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calculation
  !    INPUT:
  !   OUTPUT:
  !
   use MD_Globle_Variables_GPU, only:m_DEVICES,m_NDEVICE
   implicit none
      type(SimMDBox) ::SimBox
      type(SimMDCtrl)::CtrlParam
      !--- Local variables
      integer::err,CURDEV
      !---------------------------------------

           ERR = cudaGetDevice(CURDEV)

           !$$--- clear the memory allocated on device 1
           if(m_NDEVICE.GE.1) &
              call Clear_AtomicVirialStress_template(m_DEVICES(1), d1m_AVP, d1m_AKP, d1m_CM)

           !$$--- clear the memory allocated on device 2
           if(m_NDEVICE.GE.2) &
              call Clear_AtomicVirialStress_template(m_DEVICES(2), d2m_AVP, d2m_AKP, d2m_CM)

           !$$--- clear the memory allocated on device 3
           if(m_NDEVICE.GE.3) &
              call Clear_AtomicVirialStress_template(m_DEVICES(3), d3m_AVP, d3m_AKP, d3m_CM)

           !$$--- clear the memory allocated on device 4
           if(m_NDEVICE.GE.4) &
              call Clear_AtomicVirialStress_template(m_DEVICES(4), d4m_AVP, d4m_AKP, d4m_CM)

           !$$--- clear the memory allocated on device 5
           if(m_NDEVICE.GE.5) &
              call Clear_AtomicVirialStress_template(m_DEVICES(5), d5m_AVP, d5m_AKP, d5m_CM)

           !$$--- clear the memory allocated on device 6
           if(m_NDEVICE.GE.6) &
              call Clear_AtomicVirialStress_template(m_DEVICES(6), d6m_AVP, d6m_AKP, d6m_CM)

           ERR = cudaSetDevice(CURDEV)

           if(allocated(m_AVP)) deallocate(m_AVP)
           if(allocated(m_AKP)) deallocate(m_AKP)
           if(allocated(m_CVP)) deallocate(m_CVP)
           if(allocated(m_CKP)) deallocate(m_CKP)
           if(allocated(m_SVP)) deallocate(m_SVP)
           if(allocated(m_SKP)) deallocate(m_SKP)
           if(allocated(m_SDEN)) deallocate(m_SDEN)
           if(allocated(m_SVOL)) deallocate(m_SVOL)
           m_numAccum3 = 0
           m_numAccum2 = 0
           m_numAccum1 = 0

           m_INITED = 0
      return
  end subroutine Clear_AtomicVirialStress_DEV
  !*********************************************************************************

  !**********************************************************************************
  subroutine Allocate_LocalPressure_template(IDEV, NPRT, dAVP, dAKP, hCM, dCM)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                initialization
  !
  !     INPUT:     IDEV,  the index of device on which the memory to be allocated
  !                NPRT,  number of atoms to be considered on the device
  !                hCM,   the mass of atoms, to be copy into dCM
  !
  !     OUPUT:
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, intent(in)::NPRT
           real(KINDDF), device, dimension(:,:), allocatable::dAVP, dAKP
           real(KINDDF), dimension(:)::hCM
           real(KINDDF), device, dimension(:), allocatable::dCM

      !--- Local vairables
           integer::ERR, NGROUP

               ERR = cudaSetDevice(IDEV)
               NGROUP = size(hCM)
               allocate(dAVP(NPRT,9), dAKP(NPRT,9), dCM(NGROUP), STAT=ERR)
               if(ERR) then
                  write(*,fmt="(A, I2)") " MDPSCU Error: fail to allocate memory in AtomicVirialStress on device ", IDEV
                  write(*,fmt="(A)") "               Process to be stopped"
                  stop "1"
               end if

               dCM = hCM
             return
  end subroutine Allocate_LocalPressure_template
  !**********************************************************************************
  !**********************************************************************************
  subroutine Allocate_WorkingArray_DEV(SimBox)
  !***  PURPOSE:  to allocate and initialize the device working memories
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
  implicit none
        type(SimMDBox), intent(in)::SimBox
        !--- local variables
        integer::ERR, CURDEV,  NC(1:3)
        real(KINDDF)::cellsize(1:3)

           ERR = cudaGetDevice(CURDEV)
           if(m_NDEVICE.GE.1) &
              call Allocate_LocalPressure_template(m_DEVICES(1), m_NAPDEV, d1m_AVP, d1m_AKP, SimBox%CM, d1m_CM)

           if(m_NDEVICE.GE.2) &
              call Allocate_LocalPressure_template(m_DEVICES(2), m_NAPDEV, d2m_AVP, d2m_AKP, SimBox%CM, d2m_CM)

           if(m_NDEVICE.GE.3) &
              call Allocate_LocalPressure_template(m_DEVICES(3), m_NAPDEV, d3m_AVP, d3m_AKP, SimBox%CM, d3m_CM)

           if(m_NDEVICE.GE.4) &
              call Allocate_LocalPressure_template(m_DEVICES(4), m_NAPDEV, d4m_AVP, d4m_AKP, SimBox%CM, d4m_CM)

           if(m_NDEVICE.GE.5) &
              call Allocate_LocalPressure_template(m_DEVICES(5), m_NAPDEV, d5m_AVP, d5m_AKP, SimBox%CM, d5m_CM)

           if(m_NDEVICE.GE.6) &
              call Allocate_LocalPressure_template(m_DEVICES(6), m_NAPDEV, d6m_AVP, d6m_AKP, SimBox%CM, d6m_CM)

           if(m_Output_AverageAtomicStress .gt. 0) then
             allocate(m_AVP(dm_NPRT,9), m_AKP(dm_NPRT,9))
             m_AVP = C_ZERO
             m_AKP = C_ZERO
             m_numAccum1= 0
           end if

           if( m_Output_AverageCellStress .gt. 0) then
             !$$--- determine the number of cells
              cellsize(1:3) = m_CellSize(1:3)*SimBox%RR
              NC(1:3)       = SimBox%ZL(1:3)/cellsize(1:3)
              allocate(m_CKP(NC(1)*NC(2)*NC(3),9), m_CVP(NC(1)*NC(2)*NC(3),9))
              m_CKP = C_ZERO
              m_CVP = C_ZERO
              m_numAccum2= 0
           end if

           if(m_Output_AverageSphereStress .gt. 0) then
             allocate(m_SKP(m_NumShells, SimBox%NGROUP), m_SVP(m_NumShells, SimBox%NGROUP), &
                      m_SVOL(m_NumShells), m_SDEN(m_NumShells, SimBox%NGROUP))
             m_SKP =  C_ZERO
             m_SVP =  C_ZERO
             m_SVOL = C_ZERO
             m_SDEN = C_ZERO
             m_numAccum3 = 0
           end if

           !$$--- ending the device operations
            call SynchronizeDevices()
            ERR = cudaSetDevice(CURDEV)

            m_INITED = IOR(m_INITED,2)
            return
  end subroutine Allocate_WorkingArray_DEV
  !**********************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, CtrlParam)
    !***  PURPOSE:   to readin control parameters from a file
    !     INPUT:     fname: the file name
    !
    !     OUTPUT:    CtrlParam,
    !
    use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables
    use MiniUtilities

    implicit none
    !--- dummy vaiables
    character*(*)  :: fname
    type(SimMDCtrl):: CtrlParam
    !--- local variables
    integer::hFile, N, I,ityp, LINE, NEEDDAMP, DAMPSCHEME
    character*256::STR
    character*32::STRNUMB(12), KEYWORD

            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old')

              LINE = 0
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case("&QUICKDUMP", "&QUICKDAMP")
                             !*** get the controal paraemter of output average CP values
                             call EXTRACT_NUMB(STR,1, N, STRNUMB)
                             NEEDDAMP = ISTR(STRNUMB(1))
                             call Extract_Substr(STR,1,N,STRNUMB)
                             if(N .ge. 1) then
                                call UpCase(STRNUMB(1))
                                if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "LBFGS") then
                                    DAMPSCHEME = CP_DAMPSCHEME_LBFGS
                                else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG") then
                                    DAMPSCHEME = CP_DAMPSCHEME_CG
                                else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "CG-LS") then
                                    DAMPSCHEME = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                    
                                else if(STRNUMB(1)(1:len_trim(STRNUMB(1))) .eq. "ST") then
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
                             call EXTRACT_NUMB(STR,3,n,STRNUMB)
                             CtrlParam%JOBID0 = ISTR(STRNUMB(1))
                             if(N .GE. 2) CtrlParam%JOBID1  = ISTR(STRNUMB(2))
                             if(N .GE. 3) CtrlParam%JOBIDSTEP = ISTR(STRNUMB(3))

                        case("&CFGSEL")
                            !$$*** To get cfg range to be analysis
                            call EXTRACT_NUMB(STR,3,n,STRNUMB)
                            CtrlParam%STARTCFG = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDCFG  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%CFGSTEP = ISTR(STRNUMB(3))

                        case("&BOXSEL")
                            !$$*** To get cfg range to be analysis
                            call EXTRACT_NUMB(STR,3,n,STRNUMB)
                            CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))

                         case( mp_FTAGO)
                              call Extract_Substr(STR,1,n,m_OUTFILE)

                        case("&ASTRESS")
                            !$$*** get if record atomic stress
                            call EXTRACT_NUMB(STR,1,n,STRNUMB)
                            m_Output_AtomicStress = ISTR(STRNUMB(1))

                        case("&ASTRESSAV")
                            !$$*** get if record average atomic stress
                            call EXTRACT_NUMB(STR,1,n,STRNUMB)
                            m_Output_AverageAtomicStress = ISTR(STRNUMB(1))

                        case("&CSTRESS")
                            !$$*** get if record stress in cells
                            call EXTRACT_NUMB(STR,4,n,STRNUMB)
                            m_Output_CellStress = ISTR(STRNUMB(1))
                            if((m_Output_CellStress .gt. 0) .and. (N.lt.2)) then
                               write(*,fmt="(A)")       " MDPSCU Error: the stress in cells to be recorded, but size of cells is not given"
                               write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                               write(*,fmt="(A)")       "        Usage: &CSTRESS 1, cs (in LU)"
                               write(*,fmt="(A)")       "               &CSTRESS 1, csx, csy, csz"
                               write(*,fmt="(A)")       "               Process to be stopped"
                               stop
                            end if
                            if(N .GT. 1) m_CellSize(1:3) = DRSTR(STRNUMB(2))
                            if(N .GT. 2) m_CellSize(2:3) = DRSTR(STRNUMB(3))
                            if(N .GT. 3) m_CellSize(3:3) = DRSTR(STRNUMB(4))

                        case("&CSTRESSAV")
                            !$$*** get if record average stress in cells
                            call EXTRACT_NUMB(STR,1,n,STRNUMB)
                            m_Output_AverageCellStress = ISTR(STRNUMB(1))

                        case("&SSTRESS")
                            !$$*** get if record average stress in cells
                            call EXTRACT_NUMB(STR,7,n,STRNUMB)
                            m_Output_SphereStress = ISTR(STRNUMB(1))
                            if(m_Output_SphereStress .gt. 0 .and. N < 6) then
                               write(*,"(A)")           "MDPSCU Error: the stress in sphere region to be recorded, " &
                                                      //"but radius and number of segements of the sphere are not given."
                               write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                               write(*,fmt="(A)")       "        Usage: &SSTRESS 1, cx, cy, cz, r, num"
                               write(*,fmt="(A)")       "               where (cx, cy, cz) - the center positon (in LU) of the sphere,"
                               write(*,fmt="(A)")       "                      r           - the radius of the sphere"
                               write(*,fmt="(A)")       "                      num         - the number of sphells the sphere to be divided"
                               write(*,fmt="(A)")       "               Process to be stopped"
                               stop
                           end if
                           m_SphereCenter(1) = DRSTR(STRNUMB(2))
                           m_SphereCenter(2) = DRSTR(STRNUMB(3))
                           m_SphereCenter(3) = DRSTR(STRNUMB(4))
                           m_SphereRmax      = DRSTR(STRNUMB(5))
                           m_NumShells       = ISTR(STRNUMB(6))

                        case("&TSTRESS")
                            !$$*** get if record average stress in cells
                            call EXTRACT_NUMB(STR,mp_MXGROUP+1,n,STRNUMB)
                            m_Output_TypeStress = ISTR(STRNUMB(1))
                            if(m_Output_SphereStress .gt. 0 .and. N < 2) then
                               write(*,"(A)")           "MDPSCU Error: the stress of partical pressure of atoms to be recorded, " &
                                                      //"but type(s) of atoms is not given."
                               write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                               write(*,fmt="(A)")       "        Usage: &TSTRESS 1, type1, typ2..."
                               write(*,fmt="(A)")       "               where type1, type2 ... - the type of atoms concerned"
                               write(*,fmt="(A)")       "               Process to be stopped"
                               stop
                           end if
                           m_INCTYP = 0
                           do I=2, N
                              m_INCTYP(ISTR(STRNUMB(I))) = 1
                           end do

                        case("&SSTRESSAV")
                            !$$*** get if record average stress in spheres
                            call EXTRACT_NUMB(STR,4,n,STRNUMB)
                            m_Output_AverageSphereStress = ISTR(STRNUMB(1))

                         case default
                              write(*,"(A)")" MDPSCU Warning: unknown keyword in AtomicVirialStress control file: "//KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)
                  end select
              end do
     100     close(hFile)

             if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for AtomVirialStress calculation is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
             else
                call CreateDataFolder_Globle_Variables(m_OUTFILE)
             end if

            if(m_Output_AtomicStress .LE. 0  .AND. &
               m_Output_CellStress   .LE. 0  .AND. &
               m_Output_SphereStress .LE. 0  .AND. &
               m_Output_TypeStress .LE. 0     ) then
               write(*,fmt="(A)") "MDPSCU Warning: no record of the local stress are required"
               call ONWARNING(gm_OnWarning)
            end if

            return
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
      integer::I

           !$$--- print out the control parameters
            write(hFile,fmt="(A)") " !************ AtomicVirialStress module to be performed **********"
            write(hFile,fmt="(A)") " !    With the control parameters:"
            if(CtrlParam%NEEDDAMP .gt. 0) then
               write(hFile,FMT="(' !    Quickdampping to be performed before AtomicVirialStress calculation: ', I7, ' time steps')") CtrlParam%NEEDDAMP
            else
               write(hFile,FMT="(' !    Quickdampping to be performed before AtomicVirialStress: NO')")
            end if
            write(hFile,FMT="(' !    ')")

            if(m_Output_AtomicStress > 0) then
               write(hFile,FMT="(' !    output atomic stress:             YES')")
            else
               write(hFile,FMT="(' !    output atomic stress:             NO')")
            end if

            if(m_Output_AverageAtomicStress > 0) then
               write(hFile,FMT="(' !    output average atomic stress for every ',I5, ' configurations')") m_Output_AverageAtomicStress
            else
               write(hFile,FMT="(' !    output average atomic stress:     NO')")
            end if

            if(m_Output_CellStress > 0) then
               write(hFile,FMT="(' !    output stress in cells:           YES')")
               write(hFile,FMT="(' !        with cell sizes(in LU):', 3(1x,1PE12.3,1x))") m_CellSize(1:3)
            else
               write(hFile,FMT="(' !    output stress in cells:           NO')")
            end if

            if(m_Output_AverageCellStress > 0) then
               write(hFile,FMT="(' !    output average stress in cells for every ',I5, ' configurations')") m_Output_AverageCellStress
            else
               write(hFile,FMT="(' !    output average stress in cells:   NO')")
            end if

            if(m_Output_SphereStress > 0) then
               write(hFile,FMT="(' !    output stress in sphere:          YES')")
               write(hFile,FMT="(' !           with center of shperes:', 3(1x,1PE12.3,1x))") m_SphereCenter(1:3)
               write(hFile,FMT="(' !           with max sphere radius:', 3(1x,1PE12.3,1x))") m_SphereRmax
               write(hFile,FMT="(' !           number of sphere shell:', 3(1x,I3))") m_NumShells
            else
               write(hFile,FMT="(' !    output stress in sphere:          NO')")
            end if

            if(m_Output_AverageSphereStress > 0) then
               write(hFile,FMT="(' !    output average stress in sphere for every ',I5, ' time steps')") m_Output_AverageSphereStress
            else
               write(hFile,FMT="(' !    output average stress in sphere:  NO')")
            end if

            if(m_Output_TypeStress > 0) then
               write(hFile,FMT="(' !    output stress by types of atoms:          YES')")
               do I=1, SimBox%NGROUP
                  if(m_INCTYP(I) .gt. 0) then
                     write(hFile,FMT="(' !                 atom type included:', 20(1x,I4))") m_INCTYP(I)
                  end if
               end do
            else
               write(hFile,FMT="(' !    output stress by types of atoms:         NO')")
            end if

            write(hFile,FMT="(' !    ')")

           return
    end subroutine PrintControlParameters
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Initialize_AtomicVirialStress_DEV(SimBox, CtrlParam)
  !***  PURPOSE:   to allocate the primaritive memory for atomic pressure calculations
  !     INPUT:     SimBox: the simulation box
  !                CtrlParam: the control parameters
  !
  !     OUTPUT:   the working spaces allocated
  !
    use MD_Globle_Variables_GPU, only:Check_DEVICES
    use MD_TYPEDEF_PrintList,    only:Add_PrintProcess
    implicit none
      !--- dummy vaiables
      type(SimMDBox), intent(in)::SimBox
      type(SimMDCtrl)           ::CtrlParam

      !--- Local variables
      integer::I, IFILE


          !$$--- to check device version is used
           call Check_DEVICES()

           if(m_INITED .gt.0 ) then
             call Clear_AtomicVirialStress_DEV(SimBox, CtrlParam)
           end if

           !$$--- to findout the I/O unit
           IFILE = 0
           do I=1, SIZE(CtrlParam%f_tag)
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                 m_OUTFILE = CtrlParam%f_others(I)
              end if
              if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                 IFILE = I
              end if
           end do

            if(IFILE .le.0 ) then
               write(*,*) "MDPSCU Error: no control file for AtomicVirialStress calculation is given."
               write(*,*) "              add the keyword in SETUP file:",  mp_FTAGI
               write(*,fmt="(' Process to be stopped')")
               stop
           end if

           write(*,fmt="(A)") " !**** Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
           if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
           call LoadControlParameters(CtrlParam%f_others(IFILE), CtrlParam)
           call Add_PrintProcess(PrintControlParameters)

           CtrlParam%NEEDPOT = 1
           m_INITED = 1

           if(m_hFileTypeStress .gt. 0) close(m_hFileTypeStress)
           m_hFileTypeStress = 0

          return
  END SUBROUTINE Initialize_AtomicVirialStress_DEV
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE COPYOUT_AtomicStress(hAVP)
  !***  PURPOSE:   to copyout the atomic virial from device to host
  !
  !     INPUT:
  !     OUTPUT     hAVP,   the atomic virial on host
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      real(KINDDF), dimension(:,:), intent(out)::hAVP

      !--- Local variables
         integer::NA

         if(m_NDEVICE .GE. 1) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(1), m_STARTCELL(1), m_ENDCELL(1),  hAVP, d1m_AVP, NA)
            call COPY_OUT_SHIFT_template(1,  d1m_AVP, hAVP)

         if(m_NDEVICE .GE. 2) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(2), m_STARTCELL(2), m_ENDCELL(2),  hAVP, d2m_AVP, NA)
            call COPY_OUT_SHIFT_template(2, d2m_AVP, hAVP)

         if(m_NDEVICE .GE. 3) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(3), m_STARTCELL(3), m_ENDCELL(3),  hAVP, d3m_AVP, NA)
            call COPY_OUT_SHIFT_template(3, d3m_AVP, hAVP)

         if(m_NDEVICE .GE. 4) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(4), m_STARTCELL(4), m_ENDCELL(4),  hAVP, d4m_AVP, NA)
            call COPY_OUT_SHIFT_template(4, d4m_AVP, hAVP)

         if(m_NDEVICE .GE. 5) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(5), m_STARTCELL(5), m_ENDCELL(5),  hAVP, d5m_AVP, NA)
            call COPY_OUT_SHIFT_template(5, d5m_AVP, hAVP)

         if(m_NDEVICE .GE. 6) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(6), m_STARTCELL(6), m_ENDCELL(6),  hAVP, d6m_AVP, NA)
            call COPY_OUT_SHIFT_template(6, d6m_AVP, hAVP)

         return
  END SUBROUTINE COPYOUT_AtomicStress
  !*********************************************************************************

   !*********************************************************************************
  SUBROUTINE COPYOUT_Atomic_KineticTensor(hAKP)
  !***  PURPOSE:   to copyout the atomic kinetic tensor from device to host
  !
  !     INPUT:
  !     OUTPUT     hAKP,   the atomic kintetic tensor on host
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      real(KINDDF), dimension(:,:), intent(out)::hAKP

      !--- Local variables
         integer::NA

         if(m_NDEVICE .GE. 1) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(1), m_STARTCELL(1), m_ENDCELL(1),  hAKP, d1m_AKP, NA)
            call COPY_OUT_SHIFT_template(1, d1m_AKP, hAKP)

         if(m_NDEVICE .GE. 2) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(2), m_STARTCELL(2), m_ENDCELL(2),  hAKP, d2m_AKP, NA)
            call COPY_OUT_SHIFT_template(2, d2m_AKP, hAKP)

         if(m_NDEVICE .GE. 3) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(3), m_STARTCELL(3), m_ENDCELL(3),  hAKP, d3m_AKP, NA)
            call COPY_OUT_SHIFT_template(3, d3m_AKP, hAKP)

         if(m_NDEVICE .GE. 4) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(4), m_STARTCELL(4), m_ENDCELL(4),  hAKP, d4m_AKP, NA)
            call COPY_OUT_SHIFT_template(4, d4m_AKP, hAKP)

         if(m_NDEVICE .GE. 5) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(5), m_STARTCELL(5), m_ENDCELL(5),  hAKP, d5m_AKP, NA)
            call COPY_OUT_SHIFT_template(5, d5m_AKP, hAKP)

         if(m_NDEVICE .GE. 6) &
            !call COPYOUT_REAL8DIM2_template(m_DEVICES(6), m_STARTCELL(6), m_ENDCELL(6),  hAKP, d6m_AKP, NA)
           call COPY_OUT_SHIFT_template(6, d6m_AKP, hAKP)

         return
  END SUBROUTINE COPYOUT_Atomic_KineticTensor
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Cal_AtomicStress(SimBox, CtrlParam, hAVP)
  !***  PURPOSE:   to calculate the atomic virial part of atomic pressure tensor.
  !
  !     INPUT:     SimBox:    the simulation box
  !                CtrlParam: the control parameters
  !     OUTPUT     dxm_AVPs:  the atomic virital tensor in the module
  !                hAVP:      optional, the virial tensor on host
  !
  !     NOTE:      It has been assumed that atomic density have been calculated
  !                by calling subroutine CALDEN_FS_Force_Table2A_DEV
  !                or calling subroutine CALFORCE_FS_Force_Table2A_DEV
  !                in  MD_FS_Force_Table_GPU moudule.
  !
  !
   use MD_Globle_Variables_GPU
   use MD_ForceLib_Factory_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),                         intent(in) ::SimBox
      type(SimMDCtrl),                        intent(in) ::CtrlParam
      real(KINDDF), dimension(:,:), optional, intent(out)::hAVP

      !--- Local variables
      integer::ERR, CURDEV


             if(m_NDEVICE .GE. 1) call gm_ForceClass%pCalAVStress(1, d1m_AVP)
             if(m_NDEVICE .GE. 2) call gm_ForceClass%pCalAVStress(2, d2m_AVP)
             if(m_NDEVICE .GE. 3) call gm_ForceClass%pCalAVStress(3, d3m_AVP)
             if(m_NDEVICE .GE. 4) call gm_ForceClass%pCalAVStress(4, d4m_AVP)
             if(m_NDEVICE .GE. 5) call gm_ForceClass%pCalAVStress(5, d5m_AVP)
             if(m_NDEVICE .GE. 6) call gm_ForceClass%pCalAVStress(6, d6m_AVP)

            if(present(hAVP)) then
               call COPYOUT_AtomicStress(hAVP)
            end if

           ERR = cudaSetDevice(CURDEV)
          return
  end SUBROUTINE Cal_AtomicStress
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Cal_Atomic_KineticTensor(SimBox, CtrlParam, hAKP)
  !***  PURPOSE:   to calculate the atomic virial part of atomic pressure tensor.
  !
  !     INPUT:     SimBox:    the simulation box
  !                CtrlParam: the control parameters
  !     OUTPUT     dxm_AKP:   the atomic kinetic energy
  !                hAKP:      optional, the atomic kinetic energy on host
  !
  !
   use MD_Globle_Variables_GPU, only:m_NDEVICE
   use MD_Utilities_GPU, only:Cal_Atomic_KineticTensor_DEV
   implicit none
      !--- dummy vaiables
      type(SimMDBox),                         intent(in) ::SimBox
      type(SimMDCtrl),                        intent(in) ::CtrlParam
      real(KINDDF), dimension(:,:), optional, intent(out)::hAKP

      !--- Local variables
             if(m_NDEVICE .GE. 1) &
                call Cal_Atomic_KineticTensor_DEV(1, d1m_CM, d1m_AKP)

             if(m_NDEVICE .GE. 2) &
                call Cal_Atomic_KineticTensor_DEV(2, d2m_CM, d2m_AKP)

             if(m_NDEVICE .GE. 3) &
                call Cal_Atomic_KineticTensor_DEV(3, d3m_CM, d3m_AKP)

             if(m_NDEVICE .GE. 4) &
                call Cal_Atomic_KineticTensor_DEV(4, d4m_CM, d4m_AKP)

             if(m_NDEVICE .GE. 5) &
                call Cal_Atomic_KineticTensor_DEV(5, d5m_CM, d5m_AKP)

             if(m_NDEVICE .GE. 6) &
                call Cal_Atomic_KineticTensor_DEV(6, d6m_CM, d6m_AKP)

            if(present(hAKP)) then
               call COPYOUT_Atomic_KineticTensor(hAKP)
            end if

          return
  end SUBROUTINE Cal_Atomic_KineticTensor
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Cal_Atomic_StressTensor(SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
  !***  PURPOSE:  to calculate the instance stress tensor.
  !
  !    INPUT:  SimBox, the simulation boxs
  !            CtrlParam,    the control parameters for simulation
  !
  !    OUTPUT: hVOLS, hSTAT, the atomi volume and the state indicating
  !                          if the Voronoi volume is enclosed
  !            hAKP,         the kinetic part of the stress
  !            hAVP,         the virial part  of the stress
  !
  !    NOTE:   It is assumed the FS terms have been calculated in force calculation.
  !            So we do not recalculate the FS again.
  !            also should noted is the unit of hAKP and hAVP
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU, only:Cal_Voronoi_Volume
   implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(in)::SimBox
      type(SimMDCtrl), intent(in)::CtrlParam

      !--- Local variables
         real(KINDDF), dimension(:,:)::hAKP, hAVP
         real(KINDSF), dimension(:)::hVOLS
         integer,      dimension(:)::hSTAT
         real(KINDDF)              ::LATT, UNIT
      !----
           call Cal_Voronoi_Volume(SimBox, CtrlParam, hVOLS=hVOLS, hSTAT=hSTAT)
           call Cal_Atomic_KineticTensor(SimBox, CtrlParam, hAKP)
           call Cal_AtomicStress(SimBox, CtrlParam, hAVP)

          return
  end SUBROUTINE Cal_Atomic_StressTensor
  !*********************************************************************************

  !*********************************************************************************
  SUBROUTINE Output_Atomic_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the instance atomic stress tensor.
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !    OUTPUT:
  !
  !    NOTE:   It is assumed the embedment terms have been calculated in force calculation.
  !            So we do not recalculate the embedment again.
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU, only:mp_VSTAT_OPENED
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp),          intent(in) :: Stamp
      type(SimMDBox),               intent(in) :: SimBox
      type(SimMDCtrl),              intent(in) :: CtrlParam
      real(KINDDF), dimension(:,:), intent(in) :: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in) :: hVOLS0
      integer,      dimension(:),   intent(in) :: hSTAT0

      !--- Local variables
      real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP
      real(KINDSF), dimension(:),   allocatable::hVOLS
      integer,      dimension(:),   allocatable::hSTAT

      character*256::GFILE
      integer::hFile, I, J, IA0, JOB, ICFG
      real(KINDDF)::LATT, UNIT

      !----
           JOB    = Stamp%ITest
           ICFG   = Stamp%IRec(1)

           if(dm_NPRT .gt. SimBox%NPRT) then
              write(*,fmt="(A)")  ' MDPSCU Error: this program can handle only single box'
              write(*,fmt="(A)")  '               please seperate multi box into single box'
              write(*,fmt="(A)")  '               using tool MultBox2Boxes'
              stop
           end if

           allocate(hAKP(dm_NPRT,9), hAVP(dm_NPRT,9))
           allocate(hVOLS(dm_NPRT),hSTAT(dm_NPRT))

           !Reorder the indice of the arraies
            DO I=1, dm_NPRT
               IA0 = hm_GID(I)
               hVOLS(IA0) = hVOLS0(I)
               hSTAT(IA0) = hSTAT0(I)
               if(hSTAT(IA0) .EQ. mp_VSTAT_OPENED) then
                  hAKP(IA0, 1:9) = C_ZERO
                  hAVP(IA0, 1:9) = C_ZERO
                  cycle
               end if
               hAKP(IA0,1:9) = (hAKP0(I, 1:9)/hVOLS0(I))*CP_CGS2KBAR
               hAVP(IA0,1:9) = (hAVP0(I, 1:9)/hVOLS0(I))*CP_CGS2KBAR
            END DO

            call STRCATI(GFILE, m_OUTFILE, "_Atomic_P", m_processid, 4)
            call STRCATI(GFILE, GFILE, "_", JOB,  4)
            call STRCATI(GFILE, GFILE, ".", ICFG, 4)
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = GFILE, status='unknown')
            write(*,*) "Instant atomic stress to be saved in "//GFILE(1:len_trim(GFILE))

            LATT = SimBox%RR
            UNIT = LATT**3.D0
            write(hFile, fmt="(A)") "!--- ATOMIC VIRIAL STRESS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", SimBox%NPRT
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox%ZL/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox%BOXLOW/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox%RR*CP_CM2A
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&VOLCOL ", 5
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&DENCOL ", 6

            write(hFile, fmt="('!',A6,2x, A6, 2x, 20(A15,2x))")                 &
                             "TYPE", "STAT", "X(latt.)", "Y(latt.)", "Z(latt.)",&
                             "WS VOL(latt^3)", "DEN.(/latt^3)",                 &
                             "Pressure-K(kB)", "Pressure-V(kB)", "Pressure(kB)",&
                             "Sxx-K(kB)", "Sxx-V(kB)",                          &
                             "Sxy-K(kB)", "Sxy-V(kB)",                          &
                             "Sxz-K(kB)", "Sxz-V(kB)",                          &
                             "Syy-K(kB)", "Syy-V(kB)",                          &
                             "Syz-K(kB)", "Syz-V(kB)",                          &
                             "Szz-K(kB)", "Szz-V(kB)"
            !$$--- NOTE: m_XP is assumed has been update for the present instance by
            !$$          calling CopyAllFrom_Devices_to_Host or CopyAllStatuFrom_Devices_to_Host
            !$$          or CopyOut_SimBox_DEV
            DO I=1, dm_NPRT
               write(hFile, fmt="(1x,I6, 2x, I6, 2x, 20(1PE15.6,2x))")          &
                            m_ITYP(I), hSTAT(I), m_XP(I,1:3)/LATT,              &
                            hVOLS(I)/UNIT, UNIT/hVOLS(I),                       &
                            (hAKP(I,1)+hAKP(I,5)+ hAKP(I,9))*C_UTH,             &
                            (hAVP(I,1)+hAVP(I,5)+ hAVP(I,9))*C_UTH,             &
                            (hAKP(I,1)+hAKP(I,5)+ hAKP(I,9))*C_UTH+             &
                            (hAVP(I,1)+hAVP(I,5)+ hAVP(I,9))*C_UTH,             &
                            (hAKP(I, J), hAVP(I, J), J=1,3),                    &
                            (hAKP(I, J), hAVP(I, J), J=5,6),                    &
                            (hAKP(I, J), hAVP(I, J), J=9,9)
            END DO
            close(hFile)
           deallocate(hAKP,  hAVP,  hVOLS, hSTAT)
          return
  end SUBROUTINE Output_Atomic_StressTensor
 !*********************************************************************************

 !*********************************************************************************
  SUBROUTINE Output_AverageAtomic_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the average atomic stress tensor in gieven time steps.
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !    OUTPUT:
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp), intent(in) :: Stamp
      type(SimMDBox),      intent(in) :: SimBox
      type(SimMDCtrl),     intent(in) :: CtrlParam

      real(KINDDF), dimension(:,:), intent(in) :: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in) :: hVOLS0
      integer,      dimension(:),   intent(in) :: hSTAT0

      !--- Local variables
      real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP
      real(KINDSF), dimension(:),   allocatable::hVOLS
      integer,      dimension(:),   allocatable::hSTAT

      character*256::GFILE
      integer::hFile, I, J, IA0, JOB, ICFG
      real(KINDDF)::LATT, UNIT
      !----

           JOB    = Stamp%ITest
           ICFG     = Stamp%IRec(1)
           allocate(hAKP(dm_NPRT,9), hAVP(dm_NPRT,9))
           allocate(hVOLS(dm_NPRT),hSTAT(dm_NPRT))

           !Reorder the indice of the arraies
            DO I=1, dm_NPRT
               IA0 = hm_GID(I)
               hVOLS(IA0) = hVOLS0(I)
               hSTAT(IA0) = hSTAT0(I)
               if(hSTAT(IA0) .EQ. mp_VSTAT_OPENED) then
                  hAKP(IA0, 1:9) = C_ZERO
                  hAVP(IA0, 1:9) = C_ZERO
                  cycle
               end if
               hAKP(IA0,1:9) = (hAKP0(I, 1:9)/hVOLS0(I))*CP_CGS2KBAR
               hAVP(IA0,1:9) = (hAVP0(I, 1:9)/hVOLS0(I))*CP_CGS2KBAR
            END DO
            m_AKP(1:dm_NPRT, 1:9) = m_AKP(1:dm_NPRT, 1:9) + hAKP(1:dm_NPRT, 1:9)
            m_AVP(1:dm_NPRT, 1:9) = m_AVP(1:dm_NPRT, 1:9) + hAVP(1:dm_NPRT, 1:9)
            m_numAccum1 = m_numAccum1 + 1
            deallocate(hAKP,  hAVP,  hVOLS, hSTAT)
            if(mod(m_numAccum1,m_Output_AverageAtomicStress) .ne. 0) return

            call STRCATI(GFILE, m_OUTFILE, "_AverageAtomic_P", m_processid, 4)
            call STRCATI(GFILE, GFILE, "_", JOB, 4)
            call STRCATI(GFILE, GFILE, ".", ICFG, 4)
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = GFILE, status='unknown')
            write(*,*) "Average atomic stress to be saved in "//GFILE(1:len_trim(GFILE))

            LATT = SimBox%RR
            UNIT = LATT**3.D0

            write(hFile, fmt="(A)") "!--- ATOMIC VIRIAL STRESS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", SimBox%NPRT
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox%ZL/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox%BOXLOW/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox%RR*CP_CM2A
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&VOLCOL ", 5
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&DENCOL ", 6

            write(hFile, fmt="('!',A6,2x, 20(A15,2x))")                         &
                             "TYPE", "X(latt.)", "Y(latt.)", "Z(latt.)",        &
                             "Pressure-K(kB)", "Pressure-V(kB)", "Pressure(kB)",&
                             "Sxx-K(kB)", "Sxx-V(kB)",                          &
                             "Sxy-K(kB)", "Sxy-V(kB)",                          &
                             "Sxz-K(kB)", "Sxz-V(kB)",                          &
                             "Syy-K(kB)", "Syy-V(kB)",                          &
                             "Syz-K(kB)", "Syz-V(kB)",                          &
                             "Szz-K(kB)", "Szz-V(kB)"
            !$$--- NOTE: m_XP is assumed has been update for the present instance by
            !$$          calling CopyAllFrom_Devices_to_Host or CopyAllStatuFrom_Devices_to_Host
            !$$          or CopyOut_SimBox_DEV
            DO I=1, dm_NPRT
               write(hFile, fmt="(1x,I6, 2x, 20(1PE15.6,2x))")                  &
                            m_ITYP(I), m_XP(I,1:3)/LATT,                        &
                            (m_AKP(I,1)+m_AKP(I,5)+ m_AKP(I,9))*C_UTH,          &
                            (m_AVP(I,1)+m_AVP(I,5)+ m_AVP(I,9))*C_UTH,          &
                            (m_AKP(I,1)+m_AKP(I,5)+ m_AKP(I,9))*C_UTH+          &
                            (m_AVP(I,1)+m_AVP(I,5)+ m_AVP(I,9))*C_UTH,          &
                            (m_AKP(I, J), m_AVP(I, J), J=1,3),                  &
                            (m_AKP(I, J), m_AVP(I, J), J=5,6),                  &
                            (m_AKP(I, J), m_AVP(I, J), J=9,9)
            END DO
            close(hFile)

            !--- reset the accumalte data
            m_numAccum1 = 0
            m_AKP(1:dm_NPRT, 1:9) = C_ZERO
            m_AVP(1:dm_NPRT, 1:9) = C_ZERO
          return
  end SUBROUTINE Output_AverageAtomic_StressTensor
 !*********************************************************************************

 !*********************************************************************************
  SUBROUTINE Output_Cell_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the instance stress tensor in cells.
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !
  !    OUTPUT:
  !
  !    NOTE:   It is assumed the FS terms have been calculated in force calculation.
  !            So we do not recalculate the FS again.
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp),          intent(in) :: Stamp
      type(SimMDBox),               intent(in) :: SimBox
      type(SimMDCtrl),              intent(in) :: CtrlParam

      real(KINDDF), dimension(:,:), intent(in) :: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in) :: hVOLS0
      integer,      dimension(:),   intent(in) :: hSTAT0

      !--- Local variables
         real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP
         real(KINDDF), dimension(:), allocatable  ::hVOLS
         integer,      dimension(:), allocatable  ::hSTAT
         integer,      dimension(:), allocatable  ::hNP

         character*256::GFILE
         integer::hFile, I, J, K, L, IA0, NC(3), TNC, IXYZ(3), IC, JOB, ICFG
         real(KINDDF)::cellsize(3), LB(3), POS(3),LATT, UNIT


            JOB    = Stamp%ITest
            ICFG   = Stamp%IRec(1)

           !$$--- determine the number of cells
            LATT = SimBox%RR
            UNIT = LATT**3.D0
            cellsize(1:3) = m_CellSize(1:3)*LATT
            LB(1:3)       = SimBox%BOXLOW(1:3)
            NC(1:3)       = SimBox%ZL(1:3)/cellsize(1:3)
            cellsize      = SimBox%ZL(1:3)/dble(NC(1:3))

            TNC = NC(1)*NC(2)*NC(3)
            allocate(hAKP(TNC,9), hAVP(TNC,9))
            allocate(hVOLS(TNC),hSTAT(TNC),hNP(TNC))

            hVOLS = C_ZERO
            hAKP  = C_ZERO
            hAVP  = C_ZERO
            hSTAT = mp_VSTAT_ENCLOSED
            hNP   = 0
            !$$--- NOTE: hm_XP is assumed has been update for the present instance by
            !$$          calling CopyAllFrom_Devices_to_Host or CopyAllStatuFrom_Devices_to_Host
            !$$          or CopyOut_SimBox_DEV
            do I=1, dm_NPRT
               IXYZ(1:3)  = (hm_XP(I,1:3)-LB(1:3))/cellsize(1:3)
               IC = 1+ IXYZ(1)+NC(1)*(IXYZ(2) + IXYZ(3)*NC(2))
               if(hSTAT0(I) .EQ. mp_VSTAT_ENCLOSED) then
                   hNP(IC)      = hNP(IC) + 1
                   hVOLS(IC)    = hVOLS(IC) + hVOLS0(I)
                   hAKP(IC,1:9) = hAKP(IC,1:9) + hAKP0(I,1:9)
                   hAVP(IC,1:9) = hAVP(IC,1:9) + hAVP0(I,1:9)
               else
                   hSTAT(IC)    = mp_VSTAT_OPENED
               end if
            end do

            do IC=1, TNC
               if(hVOLS(IC) .GT. C_ZERO) then
                  hAKP(IC,1:9) = (hAKP(IC, 1:9)/hVOLS(IC))*CP_CGS2KBAR
                  hAVP(IC,1:9) = (hAVP(IC, 1:9)/hVOLS(IC))*CP_CGS2KBAR
               end if
            end do

            call STRCATI(GFILE, m_OUTFILE, "_Cell_P", m_processid, 4)
            call STRCATI(GFILE, GFILE, "_", JOB, 4)
            call STRCATI(GFILE, GFILE, ".", ICFG, 4)
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = GFILE, status='unknown')
            write(*,*) "Instant local stress in cells to be saved in "//GFILE(1:len_trim(GFILE))

            write(hFile, fmt="(A)") "!--- ATOMIC VIRIAL STRESS IN CELLS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            write(hFile, fmt="(A,1X,I8, A)")         "&NATOM    ", NC(1)*NC(2)*NC(3), " !--- Note: NATOM means here the number of cells"
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox%ZL/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox%BOXLOW/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox%RR*CP_CM2A
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&VOLCOL ", 5
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&DENCOL ", 6

            write(hFile, fmt="('!',A6, 2x, A6, 2x, 20(A15,2x))")                 &
                             "CELL", "STAT", "X(latt.)", "Y(latt.)", "Z(latt.)", &
                             "WS VOL(latt^3)", "DEN.(/latt^3)",                  &
                             "Pressure-K(kB)", "Pressure-V(kB)", "Pressure(kB)", &
                             "Sxx-K(kB)", "Sxx-V(kB)",                           &
                             "Sxy-K(kB)", "Sxy-V(kB)",                           &
                             "Sxz-K(kB)", "Sxz-V(kB)",                           &
                             "Syy-K(kB)", "Syy-V(kB)",                           &
                             "Syz-K(kB)", "Syz-V(kB)",                           &
                             "Szz-K(kB)", "Szz-V(kB)"

            IC = 0
            DO K=1, NC(3)
            DO L=1, NC(2)
            DO I=1, NC(1)
               IC = IC + 1
               POS(1) = LB(1) + cellsize(1)*(dble(I)-C_HALF)
               POS(2) = LB(2) + cellsize(2)*(dble(L)-C_HALF)
               POS(3) = LB(3) + cellsize(3)*(dble(K)-C_HALF)
               if(hVOLS(IC) .GT. C_ZERO) then
                  write(hFile, fmt="(1x, I6, 2x, I6, 2x, 20(1PE15.6,2x))")  &
                                IC, hSTAT(IC), POS(1:3)/LATT,               &
                                hVOLS(IC)/UNIT, hNP(IC)*(UNIT/hVOLS(IC)),   &
                                (hAKP(IC,1)+hAKP(IC,5)+ hAKP(IC,9))*C_UTH,  &
                                (hAVP(IC,1)+hAVP(IC,5)+ hAVP(IC,9))*C_UTH,  &
                                (hAKP(IC,1)+hAKP(IC,5)+ hAKP(IC,9))*C_UTH+  &
                                (hAVP(IC,1)+hAVP(IC,5)+ hAVP(IC,9))*C_UTH,  &
                                (hAKP(IC, J), hAVP(IC, J), J=1,3),          &
                                (hAKP(IC, J), hAVP(IC, J), J=5,6),          &
                                (hAKP(IC, J), hAVP(IC, J), J=9,9)
               else
                  write(hFile, fmt="(1x, I6, 2x, I6, 2x, 20(1PE15.6,2x))")  &
                                IC, hSTAT(IC), POS(1:3)/LATT,               &
                                cellsize(1)*cellsize(2)*cellsize(3)/UNIT, 0.D0, &
                                (hAKP(IC,1)+hAKP(IC,5)+ hAKP(IC,9))*C_UTH,  &
                                (hAVP(IC,1)+hAVP(IC,5)+ hAVP(IC,9))*C_UTH,  &
                                (hAKP(IC,1)+hAKP(IC,5)+ hAKP(IC,9))*C_UTH+  &
                                (hAVP(IC,1)+hAVP(IC,5)+ hAVP(IC,9))*C_UTH,  &
                                (hAKP(IC, J), hAVP(IC, J), J=1,3),          &
                                (hAKP(IC, J), hAVP(IC, J), J=5,6),          &
                                (hAKP(IC, J), hAVP(IC, J), J=9,9)
               end if
            END DO
            END DO
            END DO
            close(hFile)
           deallocate(hAKP,  hAVP,  hVOLS, hSTAT, hNP)
          return
  end SUBROUTINE Output_Cell_StressTensor
 !*********************************************************************************

 !*********************************************************************************
  SUBROUTINE Output_AverageCell_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the instance stress tensor in cells.
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !    OUTPUT:
  !
  !    NOTE:   It is assumed the FS terms have been calculated in force calculation.
  !            So we do not recalculate the FS again.
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp),          intent(in) :: Stamp
      type(SimMDBox),               intent(in) :: SimBox
      type(SimMDCtrl),              intent(in) :: CtrlParam

      real(KINDDF), dimension(:,:), intent(in) :: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in) :: hVOLS0
      integer,      dimension(:),   intent(in) :: hSTAT0

      !--- Local variables
         real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP
         real(KINDDF), dimension(:),   allocatable::hVOLS
         integer,      dimension(:),   allocatable::hSTAT
         integer,      dimension(:),   allocatable::hNP

         character*256::GFILE
         integer::hFile, I, J, K, L, IA0, NC(3), TNC, IXYZ(3), IC, JOB, ICFG
         real(KINDDF)::cellsize(3), LB(3), POS(3),LATT, UNIT
      !---

            JOB    = Stamp%ITest
            ICFG   = Stamp%IRec(1)

           !$$--- determine the number of cells
            LATT = SimBox%RR
            UNIT = LATT**3.D0
            cellsize(1:3) = m_CellSize(1:3)*LATT
            LB(1:3)       = SimBox%BOXLOW(1:3)
            NC(1:3)       = SimBox%ZL(1:3)/cellsize(1:3)
            cellsize      = SimBox%ZL(1:3)/dble(NC(1:3))

            TNC = size(m_CKP, dim=1)
            if(TNC .NE. NC(1)*NC(2)*NC(3)) then
               write(*,*) "MDPSCU Warning: the number of cell is inconsistent in calculating average cell stress."
               write(*,*) "                no cell stress yo be calculated."
               call ONWARNING(gm_OnWarning)
               return
            end if
            allocate(hAKP(TNC,9), hAVP(TNC,9))
            allocate(hVOLS(TNC),hSTAT(TNC),hNP(TNC))

            hVOLS = C_ZERO
            hAKP  = C_ZERO
            hAVP  = C_ZERO
            hSTAT = mp_VSTAT_ENCLOSED
            hNP   = 0
            !$$--- NOTE: hm_XP is assumed has been update for the present instance by
            !$$          calling CopyAllFrom_Devices_to_Host or CopyAllStatuFrom_Devices_to_Host
            !$$          or CopyOut_SimBox_DEV
            DO I=1, dm_NPRT
               IXYZ(1:3)  = (hm_XP(I,1:3)-LB(1:3))/cellsize(1:3)
               IC = 1+ IXYZ(1)+NC(1)*(IXYZ(2) + IXYZ(3)*NC(2))
               if(hSTAT0(I) .EQ. mp_VSTAT_ENCLOSED) then
                   hNP(IC)      = hNP(IC) + 1
                   hVOLS(IC)    = hVOLS(IC) + hVOLS0(I)
                   hAKP(IC,1:9) = hAKP(IC,1:9) + hAKP0(I,1:9)
                   hAVP(IC,1:9) = hAVP(IC,1:9) + hAVP0(I,1:9)
               else
                   hSTAT(IC)    = mp_VSTAT_OPENED
               end if
            END DO

            DO IC=1, TNC
               if(hVOLS(IC) .GT. C_ZERO) then
                  m_CKP(IC,1:9) = m_CKP(IC,1:9) + (hAKP(IC, 1:9)/hVOLS(IC))*CP_CGS2KBAR
                  m_CVP(IC,1:9) = m_CVP(IC,1:9) + (hAVP(IC, 1:9)/hVOLS(IC))*CP_CGS2KBAR
               end if
            END DO
            m_numAccum2 = m_numAccum2 + 1
            deallocate(hAKP,  hAVP,  hVOLS, hSTAT, hNP)
            if(mod(m_numAccum2,m_Output_AverageCellStress) .ne. 0) return

            m_CKP(IC,1:9) = m_CKP(IC,1:9)/dble(m_numAccum2)
            m_CVP(IC,1:9) = m_CVP(IC,1:9)/dble(m_numAccum2)

            call STRCATI(GFILE, m_OUTFILE, "_AverageCell_P", m_processid, 4)
            call STRCATI(GFILE, GFILE, "_", JOB, 4)
            call STRCATI(GFILE, GFILE, ".", ICFG, 4)
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = GFILE, status='unknown')
            write(*,*) "Average local stress in cells to be saved in "//GFILE(1:len_trim(GFILE))

            write(hFile, fmt="(A)") "!--- ATOMIC VIRIAL STRESS IN CELLS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
            write(hFile, fmt="(A,1X,I8, A)")         "&NATOM    ", NC(1)*NC(2)*NC(3), " !--- Note: NATOM means here the number of cells"
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox%ZL/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox%BOXLOW/SimBox%RR
            write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox%RR*CP_CM2A
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&VOLCOL ", 5
            write(hFile, fmt="(A,1X,3(I4,1X))")      "&DENCOL ", 6

            write(hFile, fmt="('!',A6, 2x, 20(A15,2x))")                         &
                             "CELL", "X(latt.)", "Y(latt.)", "Z(latt.)",         &
                             "WS VOL(latt^3)", "DEN.(/latt^3)",                  &
                             "Pressure-K(kB)", "Pressure-V(kB)", "Pressure(kB)", &
                             "Sxx-K(kB)", "Sxx-V(kB)",                           &
                             "Sxy-K(kB)", "Sxy-V(kB)",                           &
                             "Sxz-K(kB)", "Sxz-V(kB)",                           &
                             "Syy-K(kB)", "Syy-V(kB)",                           &
                             "Syz-K(kB)", "Syz-V(kB)",                           &
                             "Szz-K(kB)", "Szz-V(kB)"

            IC = 0
            DO K=1, NC(3)
            DO L=1, NC(2)
            DO I=1, NC(1)
               IC = IC + 1
               POS(1) = LB(1) + cellsize(1)*(dble(I)-C_HALF)
               POS(2) = LB(2) + cellsize(2)*(dble(L)-C_HALF)
               POS(3) = LB(3) + cellsize(3)*(dble(K)-C_HALF)
               write(hFile, fmt="(1x, I6, 2x, 20(1PE15.6,2x))")                &
                                 IC, POS(1:3)/LATT,                            &
                                (m_CKP(IC,1)+m_CKP(IC,5)+ m_CKP(IC,9))*C_UTH,  &
                                (m_CVP(IC,1)+m_CVP(IC,5)+ m_CVP(IC,9))*C_UTH,  &
                                (m_CKP(IC,1)+m_CKP(IC,5)+ m_CKP(IC,9))*C_UTH+  &
                                (m_CVP(IC,1)+m_CVP(IC,5)+ m_CVP(IC,9))*C_UTH,  &
                                (m_CKP(IC, J), m_CVP(IC, J), J=1,3),           &
                                (m_CKP(IC, J), m_CVP(IC, J), J=5,6),           &
                                (m_CKP(IC, J), m_CVP(IC, J), J=9,9)
            END DO
            END DO
            END DO
            close(hFile)

            m_CKP(IC,1:9) = C_ZERO
            m_CVP(IC,1:9) = C_ZERO
            m_numAccum2   = 0

          return
  end SUBROUTINE Output_AverageCell_StressTensor
 !*********************************************************************************

 !*********************************************************************************
  SUBROUTINE Output_Sphere_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the instance stress tensor in spheres into a file with the file name
  !               given in controal data
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !            hVOLS0,  the Voronoi volume of atoms
  !            hSTAT0,  the statu of Voronoi volume
  !            hAKP0,   kinetic part of atomic pressure
  !            hAVP0,   virial part of atomic pressure
  !    OUTPUT:
  !
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp),          intent(in):: Stamp
      type(SimMDBox),               intent(in):: SimBox
      type(SimMDCtrl),              intent(in):: CtrlParam
      real(KINDDF), dimension(:,:), intent(in):: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in):: hVOLS0
      integer,      dimension(:),   intent(in):: hSTAT0

      !--- Local variables
         real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP, hNP
         real(KINDDF), dimension(:),   allocatable::hVOLS
         integer,      dimension(:),   allocatable::hSTAT

         character*256::GFILE, FMT, STR, TSTR
         equivalence(GFILE, FMT)
         integer::hFile, I, J, IC, JOB, ICFG, NGROUP
         real(KINDDF)::center(3), rmin, rmax, dr, LATT, UNIT,RR


      !---

            JOB    = Stamp%ITest
            ICFG   = Stamp%IRec(1)

           !$$--- determine the number of shells
            LATT = SimBox%RR
            UNIT = LATT**3.D0
            rmax   = m_SphereRmax*LATT
            dr     = rmax/dble(m_NumShells)
            ngroup = SimBox%NGROUP

            allocate(hAKP(m_NumShells, ngroup), hAVP(m_NumShells, ngroup))
            allocate(hVOLS(m_NumShells),hSTAT(m_NumShells),hNP(m_NumShells, ngroup))

            hVOLS = C_ZERO
            hAKP  = C_ZERO
            hAVP  = C_ZERO
            hSTAT = mp_VSTAT_ENCLOSED
            hNP   = 0
            !$$--- NOTE: hm_XP is assumed has been update for the present instance by
            !$$          calling CopyAllFrom_Devices_to_Host or CopyAllStatuFrom_Devices_to_Host
            !$$          or CopyOut_SimBox_DEV
            DO I=1, dm_NPRT
               RR  = sum((hm_XP(I,1:3)-center(1:3))*(hm_XP(I,1:3)-center(1:3)))
               RR  = DSQRT(RR)
               IC  = RR/DR+1
               if(IC .LE. m_NumShells) then
                  if(hSTAT0(I) .EQ. mp_VSTAT_ENCLOSED) then
                     hNP(IC, hm_ITYP(I)) = hNP(IC, hm_ITYP(I)) + 1
                     hVOLS(IC)           = hVOLS(IC) + hVOLS0(I)
                     DO J=1, IC
                        hAKP(J,hm_ITYP(I)) = hAKP(J,hm_ITYP(I)) + (hAKP0(I,1)+hAKP0(I,5)+hAKP0(I,9))/3.D0
                        hAVP(J,hm_ITYP(I)) = hAVP(J,hm_ITYP(I)) + (hAVP0(I,1)+hAVP0(I,5)+hAVP0(I,9))/3.d0
                     END DO
                  else
                   hSTAT(IC)    = mp_VSTAT_OPENED
                  end if
               end if
            END DO

            DO I=1, m_NumShells
               if(hVOLS(I) .GT. C_ZERO) then
                  hAKP(I,1:NGROUP) = (hAKP(I, 1:NGROUP)/SUM(hVOLS(1:I)))*CP_CGS2KBAR
                  hAVP(I,1:NGROUP) = (hAVP(I, 1:NGROUP)/SUM(hVOLS(1:I)))*CP_CGS2KBAR
                  hNP(I, 1:NGROUP) =  hNP(I, 1:NGROUP)*(UNIT/hVOLS(I))
               end if
            END DO

            call STRCATI(GFILE, m_OUTFILE, "_Sphere_P", m_processid, 4)
            call STRCATI(GFILE, GFILE, "_", JOB, 4)
            call STRCATI(GFILE, GFILE, ".", ICFG, 4)
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = GFILE, status='unknown')
            write(*,*) "Instant local stress in spheres to be saved in "//GFILE(1:len_trim(GFILE))

            write(hFile, fmt="(A)") "!--- ATOMIC VIRIAL STRESS IN SPHERES CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(STR, fmt="('!',A14, 2x, A6, 20(2x,A15))")                      &
                             "RADIA", "STAT",  "SVOL(latt^3)",                   &
                             "VOL(latt^3)",  "DEN.(/latt^3)",                    &
                             "Pressure-K(kB)", "Pressure-V(kB)", "Pressure(kB)"
            DO J=1, NGROUP
               write(TSTR, fmt="(4(2x,A10,I2,A1,2x))") "DEN(",J,")",             &
                                                    "Part.P-K(",J,")",           &
                                                    "Part.P-V(",J, ")",          &
                                                    "Part.P(",J,")"
               STR = STR(1:LEN_TRIM(STR))//TSTR(1:len_trim(TSTR))
            END DO

            write(hFile, *) STR
            DO I=1, m_NumShells
               write(hFile, fmt="(1PE15.6, 2x, I6, 2x, 100(1PE15.6,2x))")        &
                                  (DR*I)/LATT, hSTAT(I),                         &
                                  CP_4PI3*(DR*I/LATT)**3.D0,                     &
                                  hVOLS(I)/UNIT, sum(hNP(I,1:NGROUP)),           &
                                  sum(hAKP(I,1:NGROUP)), sum(hAVP(I,1:NGROUP)),  &
                                  sum(hAKP(I,1:NGROUP)+hAVP(I,1:NGROUP)),        &
                                  (hNP(I,J), hAKP(I,J), hAVP(I,J), hAKP(I,J)+hAVP(I,J), J=1,NGROUP)
            END DO
            close(hFile)
           deallocate(hAKP,  hAVP,  hVOLS, hSTAT, hNP)
          return
  end SUBROUTINE Output_Sphere_StressTensor
 !*********************************************************************************


 !*********************************************************************************
  SUBROUTINE Output_AverageSphere_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the average stress tensor in shperes to a file with the filename
  !               given in control parameter
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !            hVOLS0,  the Voronoi volume of atoms
  !            hSTAT0,  the statu of Voronoi volume
  !            hAKP0,   kinetic part of atomic pressure
  !            hAVP0,   virial part of atomic pressure
  !
  !    OUTPUT:
  !
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp),          intent(in):: Stamp
      type(SimMDBox),               intent(in):: SimBox
      type(SimMDCtrl),              intent(in):: CtrlParam
      real(KINDDF), dimension(:,:), intent(in):: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in):: hVOLS0
      integer,      dimension(:),   intent(in):: hSTAT0

      !--- Local variables
         real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP, hNP
         real(KINDDF), dimension(:),   allocatable::hVOLS

         character*256::GFILE, FMT, STR, TSTR
         equivalence(GFILE, FMT)
         integer::hFile, I, J, IC, JOB, ICFG, NGROUP
         real(KINDDF)::CENTER(3), RMIN, RMAX, DR, LATT, UNIT,RR
      !---

            JOB    = Stamp%ITest
            ICFG   = Stamp%IRec(1)

           !$$--- determine the number of shells
            LATT   = SimBox%RR
            UNIT   = LATT**3.D0
            RMAX   = m_SphereRmax*LATT
            DR     = rmax/dble(m_NumShells)
            NGROUP = SimBox%NGROUP

            allocate(hAKP(m_NumShells, ngroup), hAVP(m_NumShells, ngroup),hVOLS(m_NumShells), hNP(m_NumShells, ngroup))

            hVOLS = C_ZERO
            hAKP  = C_ZERO
            hAVP  = C_ZERO
            hNP   = C_ZERO
            !$$--- NOTE: hm_XP is assumed has been update for the present instance by
            !$$          calling CopyAllFrom_Devices_to_Host or CopyAllStatuFrom_Devices_to_Host
            !$$          or CopyOut_SimBox_DEV

            !$$--- Get the volume, number of particle in sphere shells
            DO I=1, dm_NPRT
               RR  = sum((hm_XP(I,1:3)-center(1:3))*(hm_XP(I,1:3)-center(1:3)))
               RR  = DSQRT(RR)
               IC  = RR/DR+1
               if(IC .LE. m_NumShells) then
                  if(hSTAT0(I) .EQ. mp_VSTAT_ENCLOSED) then
                     hNP(IC, hm_ITYP(I)) = hNP(IC, hm_ITYP(I)) + 1
                     hVOLS(IC)           = hVOLS(IC) + hVOLS0(I)
                     hAKP(IC,hm_ITYP(I)) = hAKP(IC,hm_ITYP(I)) + (hAKP0(I,1)+hAKP0(I,5)+hAKP0(I,9))/3.D0
                     hAVP(IC,hm_ITYP(I)) = hAVP(IC,hm_ITYP(I)) + (hAVP0(I,1)+hAVP0(I,5)+hAVP0(I,9))/3.d0
                  end if
               end if
            END DO

            DO I=1, m_NumShells
               if(hVOLS(I) .GT. C_ZERO) then
                  DO J=1, NGROUP
                     m_SKP(I, J)  = m_SKP(I, J) + (SUM(hAKP(1:I, J))/SUM(hVOLS(1:I)))*CP_CGS2KBAR
                     m_SVP(I, J)  = m_SVP(I, J) + (SUM(hAVP(1:I, J))/SUM(hVOLS(1:I)))*CP_CGS2KBAR
                  END DO
                  m_SDEN(I, 1:NGROUP) = m_SDEN(I, 1:NGROUP)+ hNP(I, 1:NGROUP)*(UNIT/hVOLS(I))
               end if
                m_SVOL(I) = m_SVOL(I) + SUM(hVOLS(1:I))
            END DO
            m_numAccum3 = m_numAccum3 + 1
            deallocate(hAKP,  hAVP,  hVOLS, hNP)
            if(mod(m_numAccum3,m_Output_AverageSphereStress) .ne. 0) return

            m_SKP(1:m_NumShells, 1:NGROUP)  = m_SKP(1:m_NumShells, 1:NGROUP)/m_numAccum3
            m_SVP(1:m_NumShells, 1:NGROUP)  = m_SVP(1:m_NumShells, 1:NGROUP)/m_numAccum3
            m_SDEN(1:m_NumShells, 1:NGROUP) = m_SDEN(1:m_NumShells, 1:NGROUP)/m_numAccum3
            m_SVOL(1:m_NumShells)           = m_SVOL(1:m_NumShells)/m_numAccum3

            call STRCATI(GFILE, m_OUTFILE, "_AverageSphere_P", m_processid, 4)
            call STRCATI(GFILE, GFILE, "_", JOB,  4)
            call STRCATI(GFILE, GFILE, ".", ICFG, 4)
            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = GFILE, status='unknown')
            write(*,*) "Average local stress in spheres to be saved in "//GFILE(1:len_trim(GFILE))

            write(hFile, fmt="(A)") "!--- ATOMIC VIRIAL STRESS IN SPHERES CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
            write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
            write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
            write(hFile, fmt="(A)") '!    '

            write(STR, fmt="('!',A14, 20(2x,A15))")                             &
                             "RADIA", "SVOL(latt^3)",                           &
                             "VOL(latt^3)",  "DEN.(/latt^3)",                   &
                             "Pressure-K(kB)", "Pressure-V(kB)", "Pressure(kB)"
            DO J=1, NGROUP
               write(TSTR, fmt="(4(2x,A10,I2,A1,2x))") "DEN(",J,")",            &
                                                       "Part.P-K(",J,")",       &
                                                       "Part.P-V(",J, ")",      &
                                                       "Part.P(",J,")"
               STR = STR(1:LEN_TRIM(STR))//TSTR(1:len_trim(TSTR))
            END DO

            write(hFile, *) STR
            DO I=1, m_NumShells
               write(hFile, fmt="(1PE15.6, 2x, 100(1PE15.6,2x))")                  &
                                  (DR*I)/LATT,                                     &
                                  CP_4PI3*(DR*I/LATT)**3.D0,                       &
                                  m_SVOL(I)/UNIT, sum(m_SDEN(I,1:NGROUP)),         &
                                  sum(m_SKP(I,1:NGROUP)), sum(m_SVP(I,1:NGROUP)),  &
                                  sum(m_SKP(I,1:NGROUP) + m_SVP(I,1:NGROUP)),      &
                                  (m_SDEN(I,J), m_SKP(I,J), m_SVP(I,J), m_SKP(I,J)+m_SVP(I,J), J=1,NGROUP)
            END DO
            close(hFile)

            m_numAccum3 = 0
            m_SKP  = C_ZERO
            m_SVP  = C_ZERO
            m_SDEN = C_ZERO

          return
  end SUBROUTINE Output_AverageSphere_StressTensor
 !*********************************************************************************

 !*********************************************************************************
  SUBROUTINE Output_Type_StressTensor(Stamp, SimBox, CtrlParam, hVOLS0, hSTAT0, hAKP0, hAVP0)
  !***  PURPOSE:  to output the instance stress tensor of atoms of given type
  !
  !    INPUT:  Stamp,     the recording stamp
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !            hVOLS0,  the Voronoi volume of atoms
  !            hSTAT0,  the statu of Voronoi volume
  !            hAKP0,   kinetic part of atomic pressure
  !            hAVP0,   virial part of atomic pressure
  !
  !    OUTPUT:
  !
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp),          intent(in) :: Stamp
      type(SimMDBox),               intent(in) :: SimBox
      type(SimMDCtrl),              intent(in) :: CtrlParam

      real(KINDDF), dimension(:,:), intent(in) :: hAKP0, hAVP0
      real(KINDSF), dimension(:),   intent(in) :: hVOLS0
      integer,      dimension(:),   intent(in) :: hSTAT0

      !--- Local variables
         real(KINDDF), dimension(:),   allocatable::hAKP, hAVP, hNP
         real(KINDDF), dimension(:),   allocatable::hVOLS
         integer,      dimension(:),   allocatable::hSTAT

         character*256::GFILE, TSTR
         character*512::STR
         integer::I, IC, NGROUP, JOB,ITIME,ICFG
         real(KINDDF)::LATT, UNIT,RR, TIME
         save GFILE
      !--- dummy vaiables


            JOB    = Stamp%ITest
            ITIME  = Stamp%ITime
            TIME   = Stamp%Time
            ICFG   = Stamp%IRec(1)

            LATT = SimBox%RR
            UNIT = LATT**3.D0
            ngroup = SimBox%NGROUP
            allocate(hVOLS(ngroup), hAKP(ngroup), hAVP(ngroup), hSTAT(ngroup), hNP(ngroup))

            hVOLS = C_ZERO
            hAKP  = C_ZERO
            hAVP  = C_ZERO
            hSTAT = mp_VSTAT_ENCLOSED
            hNP   = 0
            do I=1, dm_NPRT
               IC = hm_ITYP(I)
               if(m_INCTYP(IC) .le. 0) cycle

                  if(hSTAT0(I) .EQ. mp_VSTAT_ENCLOSED) then
                     hNP(IC)    = hNP(IC) + 1
                     hVOLS(IC)  = hVOLS(IC) + hVOLS0(I)
                     hAKP(IC)   = hAKP(IC) + (hAKP0(I,1)+hAKP0(I,5)+hAKP0(I,9))/3.D0
                     hAVP(IC)   = hAVP(IC) + (hAVP0(I,1)+hAVP0(I,5)+hAVP0(I,9))/3.d0
                  else
                     hSTAT(IC)    = mp_VSTAT_OPENED
                  end if
            end do

            do IC=1, ngroup
               if(hVOLS(IC) .GT. C_ZERO) then
                  hAKP(IC) = (hAKP(IC)/hVOLS(IC))*CP_CGS2KBAR
                  hAVP(IC) = (hAVP(IC)/hVOLS(IC))*CP_CGS2KBAR
                  hNP(IC) =  hNP(IC)*(UNIT/hVOLS(IC))
               end if
            end do

            if(m_hFileTypeStress .le. 0) then
               call STRCATI(GFILE, m_OUTFILE, "_ByType_P", m_processid, 4)
               call STRCATI(GFILE, GFILE, "_", JOB, 4)
               call AvailableIOUnit(m_hFileTypeStress)
               open(UNIT=m_hFileTypeStress, file = GFILE, status='unknown')

               write(m_hFileTypeStress, fmt="(A)") "!--- STRESS BY ATOM TYPE CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
               write(m_hFileTypeStress, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
               write(m_hFileTypeStress, fmt="(A)") '!    AUTHOR: HOU Qing'
               write(m_hFileTypeStress, fmt="(A)") '!    '

               write(STR, fmt="(A10, A15)")  "!    ITIME", "     TIME(ps) ,"
               do I=1, NGROUP
                  if(m_INCTYP(I) .gt. 0) then
                     write(TSTR, fmt="(2x, (A8, I2, A1),5(2x,A9,I2,A9), )")"   STAT(",I,")",            &
                                                                           "     VOL(",I,")(lat^3) ",   &
                                                                           "     DEN(",I,")(/lat^3)",   &
                                                                           "Part.P-K(",I,")(kB)    ",   &
                                                                           "Part.P-V(",I,")(kB)    ",   &
                                                                           "  Part.P(",I,")(kB)   ,"

                     STR = STR(1:LEN_TRIM(STR))//TSTR(1:len_trim(TSTR))
                  end if
               end do
               write(m_hFileTypeStress, *) STR
               close(m_hFileTypeStress)
            end if

            call AvailableIOUnit(m_hFileTypeStress)
            open(UNIT=m_hFileTypeStress, file = GFILE, POSITION='APPEND')
              write(*,*) "Instant partial stress of atoms to be saved in "//GFILE(1:len_trim(GFILE))

              write(STR, fmt="(2x, I8, 2x, 1PE15.6)") ITIME, TIME
              do I=1, NGROUP
                 if(m_INCTYP(I) .gt. 0) then
                    write(TSTR, fmt="(2x,I8, 5(4x,1PE15.6,3x))") hSTAT(I), hVOLS(I)/UNIT, hNP(I),hAKP(I), hAVP(I), hAKP(I)+hAVP(I)
                    STR = STR(1:LEN_TRIM(STR))//TSTR(1:len_trim(TSTR))
                 end if
              end do
              write(m_hFileTypeStress, *) STR
              close(m_hFileTypeStress)
           deallocate(hAKP,  hAVP,  hVOLS, hSTAT, hNP)
          return
  end SUBROUTINE Output_Type_StressTensor
 !*********************************************************************************

 !*********************************************************************************
  SUBROUTINE RECORD_AtomicVirialStress(Stamp, SimBox, CtrlParam)
  !***  PURPOSE:  to record the instance stress tensor.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  !    OUTPUT: hVOLS, hSTAT, optional, the atomi volume and the state indicating
  !                          if the Voronoi volume is enclosed
  !
  !    NOTE:   It is assumed the embedment terms have been calculated in force calculation.
  !            So we do not recalculate the emebedment again.
  !            See: RECORD_AtomicVirialStress_TOOL
  !
   use MD_Globle_Variables_GPU
   use VoronoiTessellation_GPU
   implicit none
      !--- dummy vaiables
      type(MDRecordStamp), intent(in) :: Stamp
      type(SimMDBox),      intent(in) :: SimBox
      type(SimMDCtrl),     intent(in) :: CtrlParam

      !--- Local variables
         real(KINDDF), dimension(:,:), allocatable::hAKP, hAVP
         real(KINDSF), dimension(:),   allocatable::hVOLS
         integer,      dimension(:),   allocatable::hSTAT
      !----
           if(iand(m_Inited,2) .eq. 0)  then
              call Allocate_WorkingArray_DEV(SimBox)
           end if

           allocate(hVOLS(dm_NPRT),hSTAT(dm_NPRT),hAKP(dm_NPRT,9), hAVP(dm_NPRT,9))
           call Cal_Atomic_StressTensor(SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)

           if(m_Output_AtomicStress .GT. 0) then
              call Output_Atomic_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           if(m_Output_AverageAtomicStress .GT. 0) then
              call Output_AverageAtomic_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           if(m_Output_CellStress .GT. 0) then
              call Output_Cell_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           if(m_Output_AverageCellStress .GT. 0) then
              call Output_AverageCell_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           if(m_Output_SphereStress .GT. 0) then
              call Output_Sphere_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           if(m_Output_AverageSphereStress .GT. 0) then
              call Output_AverageSphere_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           if(m_Output_TypeStress .GT. 0) then
              call Output_Type_StressTensor(Stamp, SimBox, CtrlParam, hVOLS, hSTAT, hAKP, hAVP)
           end if

           deallocate(hVOLS, hSTAT, hAKP, hAVP)
          return
  end SUBROUTINE RECORD_AtomicVirialStress
 !*********************************************************************************

 !*********************************************************************************
  subroutine RECORD_AtomicVirialStress_TOOL(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calculate and output atomuic stress
  !                  This routine is interfaced to MD_SimBoxArray_ToolShell_12_GPU.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !    SEE ALSO:
  !            MD_SimBoxArray_ToolShell_14_GPU.F90
  !
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   use MD_ForceLib_Factory_GPU
   implicit none
       type(MDRecordStamp),          intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables

           if(.not.IfInit_ForceClass(gm_ForceClass)) then
              call Init_Forcetable_DEV(SimBox, CtrlParam, gm_ForceClass)
           end if

           if(Stamp%ITime < 0) return
           call gm_ForceClass%pCalEDen(SimBox(1), CtrlParam)
           call RECORD_AtomicVirialStress(Stamp, SimBox(1), CtrlParam)

          return
  end subroutine RECORD_AtomicVirialStress_TOOL
 !**********************************************************************************

  end module AtomicVirialStress_GPU





