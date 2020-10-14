  module COORDNUMB_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  Calculate the coordinate numbers of atoms
  !
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.cuf
  !                       MD_Globle_Variables_GPU.cuf
  !
  !                  REFERENCE____________________________________________________________________________
  !                      RefCoordNum.cuf
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
      !$$--- the id of  filename get from SimMDCtrl for I/O.
      !$$    Refering the definition of f_others in SimMDCtrl.
      character(len=9),parameter, private::mp_FTAGI="&AUXF_CNI"
      character(len=9),parameter, private::mp_FTAGO="&AUXF_CNO"
      character(len=256),  private::m_OUTFILE =""             ! filename of output data

      !--- C2050
      !integer,parameter::BLOCKSIZE = 512, BLOCKDIMX=32
      !--- C1060
      integer,parameter, private::mp_BLOCKSIZE = 256, mp_BLOCKDIMX=128

      !--- calculation control parameters
      integer, private::m_INITED = 0
      real(KINDDF), private::hm_BOND0(mp_MXGROUP,mp_MXGROUP) = -1 !1.8D0*0.5  !-- the default bond-length in lattice unit
      real(KINDDF), private::hm_BOND2(mp_MXGROUP,mp_MXGROUP)              ! the square of bond length in cm square

      !--- box information
      real(KINDDF), private::m_BOXSHAPE(3,3)                           ! siumulation box shape
      real(KINDDF), private::m_BOXSIZE(3)                              ! simulation box size
      integer, private::m_IFPD(3)                                      ! the periodic condition

      !--- box and cuoff information on devices
      real(KINDDF), constant, private::dcm_BOXSHAPE(3,3)
      real(KINDDF), constant, private::dcm_BOXSIZE(3)
      integer,      constant, private::dcm_IFPD(3)
      real(KINDDF), constant, private::dcm_BOND2(mp_MXGROUP,mp_MXGROUP) ! the square of bond length

      !--- coordinate numbers of atoms  on host
      integer, dimension(:,:),  allocatable,private::hm_CN

      !--- coordinate numbers of atoms  on devices
      private:: CN_DEVWS
      type::    CN_DEVWS
                integer, device, dimension(:,:),  allocatable::CN
      end type CN_DEVWS
      type(CN_DEVWS), dimension(m_MXDEVICE), private::dm_CN_DEVWS

      private::Allocate_WorkingArray_template, &
               Clear_WorkingArray_template,    &
               Cal_CoordNumber_KERNEL,         &
               StartOnDevice_template
  contains

 !****************************************************************************
  subroutine Allocate_WorkingArray_template(IDEV, NAPDEV, NGROUP, CN)
  !***  PURPOSE:   to allocate working space for a device
  !
  !     INPUT:     IDEV,    the index  of device on which the memory to be allocated
  !                NADVE,   the number of atoms on the device
  !                NGROUP,  the number of types of atoms
  !
  !     OUPUT:     CN,      the device memory storing the coordinate number of atoms
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV, NAPDEV, NGROUP
           integer, device, dimension(:,:), allocatable::CN

      !--- Local vairables
      integer::ERR, CURDEV


              !$$--- allocate memory for neighbore list on devices
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)
              allocate(CN(NAPDEV,NGROUP), STAT=err)

              if(ERR) then
                 write(*,*) "MDPSCU ERROR: Fail to allocate memory in COORDINATE_NUMBER_MODULE on device", IDEV
                 stop
              end if

              CN = 0
              ERR = cudaSetDevice(CURDEV)
              return
  end subroutine Allocate_WorkingArray_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray_DEV()
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
  use MD_Globle_Variables_GPU, only:dm_NPRT, dm_NGROUP, m_NAPDEV
  implicit none
  integer::I


           !$$--- allocate the memory allocated on device 1
           do I=1, m_NDEVICE
              call Allocate_WorkingArray_template(m_DEVICES(I), m_NAPDEV, dm_NGROUP, dm_CN_DEVWS(I)%CN)
           end do

           allocate(hm_CN(dm_NPRT,dm_NGROUP))
           hm_CN = 0
           m_INITED = IOR(m_INITED,2)
           return
  end subroutine Allocate_WorkingArray_DEV
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray_template(IDEV, CN)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_WorkingArray_template
  !
  !     INPUT:     IDEV,   the index of device on which the memory to be allocated
  !
  !     OUPUT:     RXP,    the allocated device memory storing the postionn of reference atoms
  !                CN,     the device memory storing the coordinate number of atoms
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           integer, device, dimension(:), allocatable::CN

      !--- Local vairables
           integer::CURDEV,ERR, NPRT

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)
              if(allocated(CN)) then
                 deallocate(CN, STAT=ERR)
                 if(ERR) goto 100
              end if
             ERR = cudaSetDevice(CURDEV)
             return

   100       ERR = cudaSetDevice(CURDEV)
             write(*,*) "MDPSCU Warning: fail to deallocate coordinate memory on device", IDEV
             call ONWARNING(gm_OnWarning)
             return
  end subroutine Clear_WorkingArray_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray_CoordNum_DEV(SimBox, CtrlParam)
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

           !$$--- clear the memory allocated on device 1
           do I=1, m_NDEVICE
              call Clear_WorkingArray_template(m_DEVICES(I), dm_CN_DEVWS(I)%CN)
           end do

           if(allocated(hm_CN)) deallocate(hm_CN)
           m_INITED = 0
      RETURN
  end subroutine Clear_WorkingArray_CoordNum_DEV
  !*********************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, SimBox, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    CtrlParam, the control parameters, menber of NEEDDAMP
  !                           could be changed.
  !
    use MD_Globle_Variables,only:CreateDataFolder_Globle_Variables
    use MiniUtilities
    implicit none
    !--- dummy vaiables
    character*(*)             ::fname
    type(SimMDBox), intent(in)::SimBox
    type(SimMDCtrl)           ::CtrlParam

    !--- local variables
    integer::hFile, N, I, J, LINE, NEEDDAMP, DAMPSCHEME
    character*256::STR
    character*32::STRNUMB(32), KEYWORD

            call AvailableIOUnit(hFile)
            open(UNIT=hFile, file = fname, status='old', err=200)

              LINE = 0
              do while(.TRUE.)
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  STR = adjustl(STR)
                  call GetKeyWord("&", STR, KEYWORD)
                  call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case( "&BONDLEN")
                             !*** get the bond len for a pair atoms
                             call EXTRACT_NUMB(STR,3, N, STRNUMB)
                             if(N .ge. 3) then
                                hm_BOND0(ISTR(STRNUMB(1)), ISTR(STRNUMB(2)) ) &
                                    = DRSTR(STRNUMB(3))
                                hm_BOND0(ISTR(STRNUMB(2)), ISTR(STRNUMB(1)) ) &
                                    = DRSTR(STRNUMB(3))
                             else
                                write(*,fmt="(A)")        " MDPSCU Erron: input of bond length is not completed"
                                write(*,fmt="(A)")        "        usage: &BONDLEN type1 type2 bond"
                                write(*,fmt="(A, BZI6)")  "               check control file at line:", LINE
                                write(*,fmt="(A)")        "        Process to be stopped"
                                stop
                             end if

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

                         case( mp_FTAGO)
                              call Extract_Substr(STR,1,n,m_OUTFILE)

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
                            !$$*** To get box range to be analysis
                            call EXTRACT_NUMB(STR,3,n,STRNUMB)
                            CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                            if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                            if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))

                         case default
                              write(*,*)" MDPSCU Warning: unknown keyword in CoordNum control file", KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)
                  end select
              end do
    100     close(hFile)

            !$$--- check input consistent
             do I=1, SimBox%NGROUP
                do J=1, SimBox%NGROUP
                   if(hm_BOND0(I, J) .lt. 0) then
                     hm_BOND0(I, J) = -CtrlParam%RU(I, J)
                   else
                     hm_BOND0(I, J) = hm_BOND0(I,J)*SimBox%RR
                   end if
                end do
             end do
             hm_BOND2 = hm_BOND0*hm_BOND0

             if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for coordination-number calculation is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
             else
                call CreateDataFolder_Globle_Variables(m_OUTFILE)
             end if
             return

    200     write(*,fmt="(' MDPSCU Error: fail to open control file in CoordNum module')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop

            return
   end subroutine LoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
  subroutine PrintControlParameters(hFile, SimBox,CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
  !
  implicit none

     !--- dummy vaiables
     type(SimMDBox),  intent(in)::SimBox
     type(SimMDCtrl), intent(in)::CtrlParam
     integer::hFile
     !--- local varaibles
     integer::I, J

           !$$--- print out the control parameters
            write(hFile,*) "!************ Coordination-number module to be performed **********"
            write(hFile,*) "!    With the control parameters:"
            do I=1, SimBox%NGROUP
               do J=I, SimBox%NGROUP
                  if(hm_BOND0(I, J) .gt. 0) then
                     write(hFile,fmt="(A, I2, A, I2, A, F7.3, A)")  " !      Bond length between type", I, " and type ", J, &
                                                                    "is: ",hm_BOND0(I,J)/SimBox%RR, " (LU)"
                  else
                     write(hFile,fmt="(A, I2, A, I2, A, F7.3, A)")  " !      Bond length between type", I, " and type ", J, &
                                                                    " is force cutoff distance ", -hm_BOND0(I,J)/SimBox%RR, " (LU)"
                     hm_BOND0(I,J) = -hm_BOND0(I,J)
                     hm_BOND0(J,I) =  hm_BOND0(I,J)
                  end if
               end do
            end do

            if(CtrlParam%NEEDDAMP .gt. 0) then
               write(hFile,FMT="(' !    Quickdampping to be performed before CSP calculation: ', I7, ' time steps')") CtrlParam%NEEDDAMP
            else
               write(hFile,FMT="(' !    Quickdampping to be performed before CSP calculation: NO')")
            end if
            write(hFile,FMT="(' !    ')")

           return
  end subroutine PrintControlParameters
  !****************************************************************************


  !*********************************************************************************
  subroutine Initialize_CoordNumber_DEV(SimBox, CtrlParam)
    !***  PURPOSE:   to intialize the parameters to be used in calculating coordinate number
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !
    !     OUTPUT:   the working spaces allocated
    !
    use MD_Globle_Variables_GPU, only:Check_DEVICES
    use MD_TYPEDEF_PrintList,    only:Add_PrintProcess
    implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(inout)::SimBox
      type(SimMDCtrl), intent(in)   ::CtrlParam

      !--- Local variables
      integer::I, IFILE

          !$$--- to check device version is used
           call Check_DEVICES()

           if(m_INITED .gt. 0) then
             call Clear_WorkingArray_CoordNum_DEV(SimBox, CtrlParam)
           end if

           !$$--- to findout the I/O unit
            IFILE = 0
            do I=1, SIZE(CtrlParam%f_tag)
               if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                 IFILE = I
               end if
               if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                 m_OUTFILE = CtrlParam%f_others(I)
               end if
            end do
            if(IFILE .LE.0 ) then
               write(*,*) "MDPSCU Error: no control file for coordination-number calculation is given."
               write(*,*) "              add the keyword in SETUP file:",  mp_FTAGI
               write(*,fmt="(' Process to be stopped')")
               stop
            end if

            write(*,fmt="(A)") " !**** Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
            if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
            call LoadControlParameters(CtrlParam%f_others(IFILE), SimBox, CtrlParam)
            call Add_PrintProcess(PrintControlParameters)

            m_INITED = IOR(m_INITED,1)

           return
    end subroutine Initialize_CoordNumber_DEV
   !*********************************************************************************


  !*********************************************************************************
   attributes(global) subroutine Cal_CoordNumber_KERNEL(IM, NGROUP, ITYP,XP, NAPDEV, IA0, NPRT, KVOIS,INDI, CN)
  !***  PURPOSE:   to calculate the inverse of  the electron density DEN
  !
  !$$   INPUT:     IM:     the number of particles concerned
  !$$              NGROUP: the number of particle types
  !$$              ITYP :  the type of atoms corresponding to INDI
  !$$              XP:     the positions of atoms
  !$$              NAPDEV: the max number of atoms on a device
  !$$              IA0:    the index (in the whole box)of the fisrt atom on the device
  !$$              NPRT:   the actual number of atoms on the device
  !$$              KVOIS:  the number of neighbors for atoms
  !$$              INDI:   the index for the neighbores
  !
  !$$   OUTPUT:    CN,     the coordeinate number of atoms
  !
   use MD_CONSTANTS
      implicit none
      !--- DUMMY VARIABLES
      integer, value, intent(in)::IM, NGROUP
      integer, device, intent(in)::ITYP(IM)
      real(KINDDF), device, intent(in)::XP(IM,3)

      integer, value, intent(in)::NAPDEV, IA0, NPRT
      integer, device, intent(in)::KVOIS(NAPDEV)
      integer, device, intent(in)::INDI(NAPDEV,*)

      integer, device::CN(NAPDEV,*)

      !Local variables
      integer::J,K,KK, IW, IIW, KTAB, ITYPI, ITYPJ
      real(KINDDF)::R2,R
      real(KINDDF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, DXYZ_X, DXYZ_Y, DXYZ_Z

      real(KINDDF),shared::HBX, HBY, HBZ, BX, BY, BZ, BS11, BS12, BS13, BS21, BS22, BS23, BS31, BS32, BS33
      integer,shared::IFPDX, IFPDY, IFPDZ
      integer, shared::NT, NB

      ! Local variables
      integer::IC, IT, IB, LL, NAL,NL, CN0(mp_MXGROUP)

              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              if(IT .EQ. 1) then
                 BX = dcm_BOXSIZE(1)
                 BY = dcm_BOXSIZE(2)
                 BZ = dcm_BOXSIZE(3)

                 HBX = BX*C_HALF
                 HBY = BY*C_HALF
                 HBZ = BZ*C_HALF

                 BS11 =  dcm_BOXSHAPE(1,1)
                 BS12 =  dcm_BOXSHAPE(1,2)
                 BS13 =  dcm_BOXSHAPE(1,3)
                 BS21 =  dcm_BOXSHAPE(2,1)
                 BS22 =  dcm_BOXSHAPE(2,2)
                 BS23 =  dcm_BOXSHAPE(2,3)
                 BS31 =  dcm_BOXSHAPE(3,1)
                 BS32 =  dcm_BOXSHAPE(3,2)
                 BS33 =  dcm_BOXSHAPE(3,3)

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
            DO LL=1, NL

              !$$IC -- the id of the atom on the device
              IC= (IB-1)*NT+IT +(LL-1)*NAL

              if(IC.LE.NPRT) then
                 !$$POS-- position of the atom
                 !$$      NOTE: XP(I) is the position of Ith atom in the whole box
                 POSX = XP(IC+IA0,1)
                 POSY = XP(IC+IA0,2)
                 POSZ = XP(IC+IA0,3)
                 ITYPI  = ITYP(IC+IA0)

                 !$$-- start scanning the neighbores
                 IIW = KVOIS(IC)
                 CN0 = C_ZERO
                 DO IW=1, IIW
                    !$$--- NOTE: the particles index of neighbore-list is
                    !$$          the index of particle in the whole box
                    J=INDI(IC,IW)

                    !$$--- To calculate the seperation between particle IC
                    !$$    and its IWth neighbore
                    SEPX = POSX - XP(J,1)
                    SEPY = POSY - XP(J,2)
                    SEPZ = POSZ - XP(J,3)
                    if(IFPDX.GT.0 .AND.(DABS(SEPX) .GT. HBX)) then
                      SEPX = SEPX - DSIGN(BX,SEPX)
                    end if

                    if(IFPDY.GT.0 .AND.(DABS(SEPY) .GT. HBY)) then
                       SEPY = SEPY - DSIGN(BY,SEPY)
                    end if

                    if(IFPDZ.GT.0 .AND.(DABS(SEPZ) .GT. HBZ)) then
                       SEPZ = SEPZ - DSIGN(BZ,SEPZ)
                    end if

                    !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                    !$$          only when the shape of the box is cubic,
                    !$$          DXYZ = SEP
                    !DXYZ_X = BS11*SEPX + BS12*SEPY + BS13*SEPZ
                    !DXYZ_Y = BS21*SEPX + BS22*SEPY + BS23*SEPZ
                    !DXYZ_Z = BS31*SEPX + BS32*SEPY + BS33*SEPZ
                    !R2  = DXYZ_X*DXYZ_X + DXYZ_Y*DXYZ_Y + DXYZ_Z*DXYZ_Z
                     R2 = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                    !$$--- To calculate electron density on atom I
                    ITYPJ=ITYP(J)
                    if(R2 .LE. dcm_BOND2(ITYPI,ITYPJ)) then
                       CN0(ITYPJ)= CN0(ITYPJ) + 1
                    end if
                 ENDDO

                 CN(IC,1:NGROUP) = CN0(1:NGROUP)
              end if
            END DO
        RETURN
  end subroutine Cal_CoordNumber_KERNEL
  !*********************************************************************************

  !*********************************************************************************
  subroutine StartOnDevice_template(IDEV, STARTCELL, ENDCELL, KVOIS,INDI,ITYP, XP, CN)
  !***  PURPOSE:   to calculate Finnis-Sinclar term.
  !
  !     INPUT:     IDEV,      the ID of device
  !                STARTCELL, the ID of the first cell on the device
  !                ENDCELL,   the ID of the last cell on the device
  !
  !                KVOIS:  the number of neighbors for atoms
  !                INDI:   the index for the neighbores
  !                ITYP :  the type of atoms corresponding to INDI
  !                XP:     the positions of atoms
  !
  !     OUTPUT:    CN,     the coordinate number
  !
  !
   use MD_Globle_Variables_GPU
   implicit none
      !--- dummy vaiables
      integer::IDEV, STARTCELL, ENDCELL
      integer, device, dimension(:)::KVOIS
      integer, device, dimension(:,:)::INDI
      integer, device, dimension(:)::ITYP
      real(KINDDF), device, dimension(:,:)::XP
      integer, device, dimension(:,:)::CN

      !--- Local variables
         integer::BX, BY, NB, NBX, NBY, ERR, STARTA, ENDA, NPRT

      !$$--- Device variables and variables to be used in GPU
         type(dim3) :: blocks
         type(dim3) :: threads

             ERR = cudaSetDevice(IDEV)
            !$$--- If pressure effect acounted for, the box size could be changed, we
             !$$    need to recopy the boxsize
            dcm_BOXSHAPE = m_BOXSHAPE
            dcm_BOXSIZE  = m_BOXSIZE
            dcm_IFPD     = m_IFPD
            dcm_BOND2    = hm_BOND2
            STARTA   = hm_IA1th(STARTCELL)
            ENDA     = (hm_IA1th(ENDCELL)-1) + hm_NAC(ENDCELL)
            NPRT = ENDA - STARTA + 1

            BX  = mp_BLOCKSIZE
            BY  = 1
            NBX = mp_BLOCKDIMX
            NBY = 1

            blocks  = dim3(NBX, NBY, 1)
            threads = dim3(BX, BY, 1)
            STARTA  = STARTA-1
            call Cal_CoordNumber_KERNEL<<<blocks, threads>>>(dm_NPRT,dm_NGROUP, ITYP,XP, m_NAPDEV, STARTA, NPRT, KVOIS,INDI, CN)
          return
  end subroutine StartOnDevice_template
  !*********************************************************************************


  !*********************************************************************************
  subroutine Cal_CoordNumber_DEV(SimBox, CtrlParam, CN)
  !***  PURPOSE:  to calculate the coordinate number
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     CN, optional, the coordinate number of particles indexed in
  !                    original box
  !
  !
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   implicit none
      !--- dummy vaiables
      type(SimMDBox),                     intent(inout)::SimBox
      type(SimMDCtrl),                    intent(in)   ::CtrlParam
      integer, dimension(:,:),allocatable,optional     ::CN

      !--- Local variables
         integer::I

          if(iand(m_INITED,2) .eq. 0) then
            !     NOTE:   intialization of this module must be called after
            !             calling Initialize_Globle_Variables_DEV, with the number
            !             of particle on each device is avaliable
            call  Allocate_WorkingArray_DEV()
          end if

        !$$--- If pressure effect acounted for, the box size could be changed, we
        !$$    need to recopy the boxsize
             m_BOXSHAPE(1:3,1:3) = SimBox%BOXSHAPE(1:3,1:3)
             m_BOXSIZE(1:3)      = SimBox%ZL(1:3)
             m_IFPD = CtrlParam%IFPD

        !$$--- caluclate the CN on atoms
             do I=1, m_NDEVICE
                 call StartOnDevice_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),             &
                                             dm_Neighbors%KVOIS(I)%Data, dm_Neighbors%INDI(I)%Data,  &
                                             dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data,    &
                                             dm_CN_DEVWS(I)%CN)
             end do

            !$$--- copy the coordinate number out
            if(present(CN)) then
               call COPYOUT_CoordinateNumber_DEV(CN)
               call SynchronizeDevices()
            end if

          return
  end subroutine Cal_CoordNumber_DEV
  !*********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_CoordinateNumber_DEV(CN)
  !***  PURPOSE:  to copy the coordinate number in devices to a host array CN
  !               the indice in CN is the indice of atoms in the whole system
  !
  !
  !     INPUT:
  !     OUTPUT     CN,      the centrosysmetry
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT, hm_GID, dm_NGROUP, COPY_OUT_SHIFT_template, IndexToOriginalBox

  implicit none
  !----   DUMMY Variables
          integer, dimension(:,:),allocatable::CN

  !----   Local variables
          integer::I, J, pinnedFlag, ERR
          integer, pinned, dimension(:,:),allocatable::TCN

  !***

            !$$--- determine the array size of the neighbore list
            if(.not.allocated(CN)) then
               allocate(CN(dm_NPRT, dm_NGROUP))
            else if(size(CN,dim=1).NE.dm_NPRT .OR. size(CN,dim=2).NE.dm_NGROUP) then
               deallocate(CN)
               allocate(CN(dm_NPRT, dm_NGROUP))
            end if

            allocate(TCN(dm_NPRT,dm_NGROUP), STAT=ERR, PINNED=pinnedFlag)
            do I=1, m_NDEVICE
               call COPY_OUT_SHIFT_template(I, dm_CN_DEVWS(I)%CN, TCN)
            end do

            !$$**** restore the array in original order
            call SynchronizeDevices()
            call IndexToOriginalBox(TCN, CN)
            deallocate(TCN)

        return
  end subroutine COPYOUT_CoordinateNumber_DEV
  !**********************************************************************************

  !**********************************************************************************
  subroutine RECORD_CoordinateNumber(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the coorodinate numbers to output file. This routine is to interfaced
  !                  to and MD_EM_TB_ForceTable_SHELL_12_GPU. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_AppShell_16_GPU.f90
  !            MD_SimBoxArray_ToolShell_14_GPU.f90
  !
   use MD_Globle_Variables_GPU
   implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE, FMT
       equivalence (GFILE, FMT)
       integer::hFile, I, K, IP,  IB, JOB
       character*32::tStr
       integer, dimension(:), allocatable::HIS
       integer::hmi, hmx, ss

           call Cal_CoordNumber_DEV(SimBox(1), CtrlParam, hm_CN)

          !$$--- output CN space distribution
           JOB = Stamp%ITest
           IP  = Stamp%IRec(1)
           call STRCATI(GFILE, m_OUTFILE, "P", m_processid, 4)
           call STRCATI(GFILE, GFILE, "_", JOB, 4)
           call STRCATI(GFILE, GFILE, ".", IP, 4)

           call AvailableIOUnit(hFile)
           write(*,*) "Save CN distribution to "//GFILE(1:len_trim(GFILE))
           open(UNIT=hFile, file = GFILE, status='unknown')
             write(hFile, fmt="(A)") "!--- COORDINATION NUMBER RESULTS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
             write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
             write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
             write(hFile, fmt="(A)") '!    '

             !$$--- output CN histogram

             hmi = dm_NPRT
             hmx = 0
             do I=1, dm_NPRT
                ss = sum(hm_CN(I,1:dm_NGROUP))
                if(hmi .gt. ss) hmi = ss
                if(hmx .lt. ss) hmx = ss
             end do

             allocate(HIS(hmi-1:hmx+1))
             HIS = 0
             do I=1, dm_NPRT
                IB = sum(hm_CN(I,1:dm_NGROUP))
                HIS(IB) = HIS(IB)+1
             end do
             write(hFile, fmt="(A, I4, 1x, A, I4)")"!--- Histogram of coordinate numbers with min CN", hmi, " and max CN ", hmx
             write(hFile,fmt="(A,2X,6(A13,1X))")   "!---  CN", "Count"
             do I=hmi-1, hmx+1
                write(hFile,fmt="(A, I8,2X,6(1PE13.4,1X))")"!", I, dble(HIS(I))
             end do
             deallocate(HIS)
             write(hFile,fmt="(A)")   "!"
            write(hFile,fmt="(A)")   "!---  SYMBOL EXPLAINATION:"
             write(hFile,fmt="(A)")   "!     GCN:   The coordination number between an atom and different types of atoms"
             write(hFile,fmt="(A)")   "!     TCN:   The total coordination number of an atom"
             write(hFile,fmt="(A)")   "!"

             write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
             write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", size(SimBox)*SimBox(1)%NPRT
             write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox(1)%ZL/SimBox(1)%RR
             write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox(1)%BOXLOW/SimBox(1)%RR
             write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox(1)%RR*CP_CM2A
             write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 2, 3, 4
             write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 1
             write(hFile, fmt="(A,1X,20(I4,1X))")     "&GCNCOL   ", (4+I, I=1, SimBox(1)%NGROUP)
             write(hFile, fmt="(A,1X,20(I4,1X))")     "&TCNCOL   ", 4+SimBox(1)%NGROUP+1

             write(tstr,*) SimBox(1)%NGROUP+1
             tstr = adjustl(tstr)
             FMT = "(1x,A9,3(1x,A13),"//tstr(1:len_trim(tstr))//"(1x,A5))"
             write(hFile,fmt=FMT) "!--- TYPE", "X", "Y","Z", ("GCN", I=1,SimBox(1)%NGROUP), "TCN"

             FMT = "(5x,I5,3(1x,1PE13.5),"//tstr(1:len_trim(tstr))//"(1x,I5))"
             IP = 0
             do I=1, size(SimBox)
                do K=1, SimBox(I)%NPRT
                   IP = IP + 1
                   write(hFile,fmt=FMT)SimBox(I)%ITYP(K), SimBox(I)%XP(K,1:3)/SimBox(I)%RR, hm_CN(IP,1:dm_NGROUP), sum(hm_CN(IP,1:dm_NGROUP))
                end do
             end do
           close(hFile)

           return
  end subroutine RECORD_CoordinateNumber
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_CoordNumber_TOOL(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calculate and output the coordinate numbers to output file.
  !                  This routine is to interfaced to and MD_SimBoxArray_ToolShell_12_GPU.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_AppShell_16_GPU.f90
  !            MD_SimBoxArray_ToolShell_14_GPU.f90
  !
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables


            if(Stamp%ITime .LT. 0) return

            call RECORD_CoordinateNumber(Stamp, SimBox, CtrlParam)
          return
  end subroutine RECORD_CoordNumber_TOOL
  !****************************************************************************************
  end module COORDNUMB_GPU





