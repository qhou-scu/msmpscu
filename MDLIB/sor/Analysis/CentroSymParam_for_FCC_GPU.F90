  module CentroSymParam_for_FCC_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  To calculate the centrosymmetry parameter for atoms in a FCC system.
  !                  Ref. Phys.Rev.B58(1998)11085,
  !                       Phys.Rev.B71(2005)064112,
  !                       Philo. Mag. 89 (2009)3133
  !                  for the introduaction of centrosymmetry parameter
  !
  !
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.cuf
  !                       MD_Globle_Variables_GPU
  !
  !                  ______________________________________________________________________________________
  !                  REFER TO: CentroSymParam_for_BCC_GPU. The only different between
  !                            this and CentroSymParam_for_BCC_GPU is that mp_NPT =12
  !                            for FCC, while mp_NPT=14 for BCC. And the vectors
  !                            mp_PVECTX, mp_PVECTY, mp_PVECTZ are different.
  !**** HISTORY:
  !                  version 1st 2015-01, adapted from CentroSymParam_for_BCC_GPU  (Hou Qing)
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_MultiGPU_Basic
  !-----------------------------------------------
    implicit none


      integer::m_processid=0
      !--- the id of  filename get from SimMDCtrl for I/O
      !    ref the definition of f_others in SimMDCtrl.
      character(len=10),parameter, private::mp_FTAGI="&AUXF_CSPI"
      character(len=10),parameter, private::mp_FTAGO="&AUXF_CSPO"
      character(len=256),  private::m_OUTFILE =""             ! filename of output data

      !--- C2050
      !integer,parameter::BLOCKSIZE = 512, BLOCKDIMX=32
      !--- C1060
      integer,parameter, private::mp_BLOCKSIZE = 256, mp_BLOCKDIMX=128
      integer, private::m_INITED = 0

      !--- instant CP on host
      real(KINDSF), dimension(:), allocatable,private::m_CP

      !---   accumulating number of calculations and accumulating CSPs
      integer, private::m_Output_AverageCP  = 0                           !time step for output average CP
      integer, private::m_numAccum = 0
      real(KINDSF),dimension(:), allocatable, private::m_ACP              !accumulated CP

      !--  the threshhold of CP value
      integer, private::m_nCPTRETH= 0
      real(KINDSF),dimension(:), allocatable, private::m_CPTRETH           !(4)=(/0.75,1.0,1.5,2.0/)

      !--  the axsis in crystal
      integer, private:: m_AXSISX(3) = (/1, 0, 0/)
      integer, private:: m_AXSISY(3) = (/0, 1, 0/)
      integer, private:: m_AXSISZ(3) = (/0, 0, 1/)

      !---parameters used in calculating Cp
      integer,      parameter, private::mp_NNC=27
      real(KINDDF), constant,  private::dcm_NIX(mp_NNC)
      real(KINDDF), constant,  private::dcm_NIY(mp_NNC)
      real(KINDDF), constant,  private::dcm_NIZ(mp_NNC)

      integer,      parameter, private::mp_NPT=12
      real(KINDDF), parameter, private::mp_PVECTX(mp_NPT)=(/0.0D0,  0.0D0,  0.0D0,  0.0D0,  &
                                                            0.5D0, -0.5D0,  0.5D0, -0.5D0,  &
                                                            0.5D0,  0.5D0, -0.5D0, -0.5D0/)
      real(KINDDF), parameter, private::mp_PVECTY(mp_NPT)=(/0.5D0,  0.5D0, -0.5D0, -0.5D0,  &
                                                            0.0D0,  0.0D0,  0.0D0,  0.0D0,  &
                                                            0.5D0, -0.5D0,  0.5D0, -0.5D0/)
      real(KINDDF), parameter, private::mp_PVECTZ(mp_NPT)=(/0.5D0, -0.5D0,  0.5D0, -0.5D0,  &
                                                            0.5D0,  0.5D0, -0.5D0, -0.5D0,  &
                                                            0.0D0,  0.0D0,  0.0D0,  0.0D0/)

      real(KINDDF), private::m_RPVECTX(mp_NPT)
      real(KINDDF), private::m_RPVECTY(mp_NPT)
      real(KINDDF), private::m_RPVECTZ(mp_NPT)
      real(KINDDF), constant, private::dcm_RPVECTX(mp_NPT)
      real(KINDDF), constant, private::dcm_RPVECTY(mp_NPT)
      real(KINDDF), constant, private::dcm_RPVECTZ(mp_NPT)

      !--- CP  on devices
      private:: CP_DEVWS
      type::CP_DEVWS    
            real(KINDSF), device, dimension(:),  allocatable::CP
      end type CP_DEVWS
      type(CP_DEVWS), dimension(m_MXDEVICE), private::dm_CPDEVWS
      !****
      private::Allocate_CPArrays_template,    &
               Allocate_WorkingArray_DEV,     &
               Clear_CPArray_template,        &
               Cal_CP_KERNEL,                 &
               StartOnDevice_template,        &
               LoadControlParameters,         &
               PrintControlParameters

  contains

 !****************************************************************************
  subroutine Allocate_CPArrays_template(IDEV,NA0, CP)
  !***  PURPOSE:   to allocate working space for a device
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !                NA0,     the number of atoms on the device
  !
  !     OUPUT:     CP,      the allocated device memory storing the centrosymmetry parameters
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV, NA0
           real(KINDSF), device, dimension(:), allocatable::CP

      !--- Local vairables
      integer::ERR, CURDEV, NPRT
      integer::NIX(mp_NNC)=(/0,-1,-1,-1, 0, 0, -1, 1,-1, 0, 1,-1, 0, 1, 1, 1, 1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1/)
      integer::NIY(mp_NNC)=(/0, 0,-1, 1, 1, 0,  0, 0,-1,-1,-1, 1, 1, 1, 0, 1,-1,-1, 0, 0, 0,-1,-1,-1, 1, 1, 1/)
      integer::NIZ(mp_NNC)=(/0, 0, 0, 0, 0, 1,  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
        
  

              !$$--- allocate memory for neighbore list on devices
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)
              allocate(CP(NA0), STAT=err)
              CP = 0.D0
              if(ERR) then
                 write(*,*) "MDPSCU Error: Fail to allocate device memory for centrosymmetry calculation on device", IDEV
                 stop
              end if

              dcm_NIX = NIX
              dcm_NIY = NIY
              dcm_NIZ = NIZ
              ERR = cudaSetDevice(CURDEV)

              return
  end subroutine Allocate_CPArrays_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray_DEV()
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
  use MD_Globle_Variables_GPU, only:m_NAPDEV,dm_NPRT

  implicit none
  integer::I

           !$$--- allocate the memory allocated on device 1-6
            do I=1, m_NDEVICE
               call Allocate_CPArrays_template(m_DEVICES(I), m_NAPDEV, dm_CPDEVWS(I)%CP)
            end do

            allocate(m_CP(dm_NPRT))
            if(m_Output_AverageCP) allocate(m_ACP(dm_NPRT))
            m_CP = 0.D0
            m_ACP = 0.D0
            m_numAccum = 0

            m_INITED = IOR(m_INITED,2)
            return;
  end subroutine Allocate_WorkingArray_DEV
  !****************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, SimBox, CtrlParam)
    !***  PURPOSE:   to readin control parameters from a file
    !     INPUT:     fname: the file name
    !
    !     OUTPUT:    CtrlParam, the control parameters, menber of NEEDDAMP
    !                           could be changed.
    !
    use MD_Globle_Variables, only:Load_ExAnalyCtl_SimMDCtrl, CreateDataFolder_Globle_Variables
    implicit none
    !--- dummy vaiables
    character*(*)   ::fname
    type(SimMDBox)  ::SimBox
    type(SimMDCtrl) ::CtrlParam

    !--- local variables
    integer::I, N, IS, LINE
    character*256::STR
    character*32::STRNUMB(32), KEYWORD
    type(StatementList), pointer::StatList

            !$$--- load input statements
            call Load_ExAnalyCtl_SimMDCtrl(mp_FTAGI, Fname, SimBox, CtrlParam, &
                                              OutTag=mp_FTAGO, OutFile=m_OUTFILE, TheInput=StatList)



            do IS=1, Number_StatementList(StatList)
               call Get_StatementList(IS, StatList, STR, Line=LINE)
               STR = adjustl(STR)
               call GetKeyWord("&", STR, KEYWORD)
               call UpCase(KEYWORD)
               select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                      case( "&CPLEVELS")
                           !*** get the CP levels to distinguish the atoms
                            call Extract_Numb(STR,1, N, STRNUMB)
                            m_nCPTRETH = ISTR(STRNUMB(1))
                            if(allocated(m_CPTRETH)) deallocate(m_CPTRETH)
                               if(m_nCPTRETH > 0) then
                                  allocate(m_CPTRETH(m_nCPTRETH))
                                  call Extract_Numb(STR,m_nCPTRETH+1, N, STRNUMB)
                                  do I=1, m_nCPTRETH
                                     m_CPTRETH(I) = DRSTR(STRNUMB(I+1))
                                  end do
                            end if

                      case( "&AXSISX")
                            !*** get the CP levels to distinguish the atoms
                            call Extract_Numb(STR,3, N, STRNUMB)
                            if(N .LT. 3) then
                               write(*,fmt="(A)")        " MDPSCU Erron: the Miller index for X axsis is not completed"
                               write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                               write(*,fmt="(A)")       " Process to be stopped"
                               stop
                            end if
                            m_AXSISX(1) = ISTR(STRNUMB(1))
                            m_AXSISX(2) = ISTR(STRNUMB(2))
                            m_AXSISX(3) = ISTR(STRNUMB(3))

                      case( "&AXSISY")
                           !*** get the CP levels to distinguish the atoms
                           call Extract_Numb(STR,3, N, STRNUMB)
                           if(N .LT. 3) then
                              write(*,fmt="(A)")        " MDPSCU Erron: the Miller index for Y axsis is not completed"
                              write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                              write(*,fmt="(A)")       " Process to be stopped"
                              stop
                           end if
                           m_AXSISY(1) = ISTR(STRNUMB(1))
                           m_AXSISY(2) = ISTR(STRNUMB(2))
                           m_AXSISY(3) = ISTR(STRNUMB(3))

                      case( "&AXSISZ")
                           !*** get the CP levels to distinguish the atoms
                           call Extract_Numb(STR,3, N, STRNUMB)
                           if(N .LT. 3) then
                              write(*,fmt="(A)")        " MDPSCU Erron: the Miller index for Z axsis is not completed"
                              write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                              write(*,fmt="(A)")       " Process to be stopped"
                              stop
                           end if
                           m_AXSISZ(1) = ISTR(STRNUMB(1))
                           m_AXSISZ(2) = ISTR(STRNUMB(2))
                           m_AXSISZ(3) = ISTR(STRNUMB(3))

                      case( "&CAL_CPAV")
                           !*** get the controal paraemter of output average CP values
                           call Extract_Numb(STR,1, N, STRNUMB)
                           m_Output_AverageCP = ISTR(STRNUMB(1))

               end select
            end do
            call Release_StatementList(StatList)
           !&&--- check input consisten
            if(len_trim(m_OUTFILE) .LE.0 ) then
                write(*,fmt="(A)")          " MDPSCU Error: no output file for CentroSymParam_for_FCC calculation is given."
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
                write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
                write(*,fmt="(' Process to be stopped')")
                stop
            else
                call CreateDataFolder_Globle_Variables(m_OUTFILE)
            end if
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

           !$$--- print out the control parameters
            write(hFile,*) "!************ CentroSymParam_for_FCC module to be performed **********"
            write(hFile,*) "!    With the control parameters:"

            write(hFile,FMT="(' !    X-Axsis (in Miller indice): ' ,I5, 1x, I5, 1x, I5)") m_AXSISX(1:3)
            write(hFile,FMT="(' !    Y-Axsis (in Miller indice): ' ,I5, 1x, I5, 1x, I5)") m_AXSISY(1:3)
            write(hFile,FMT="(' !    Z-Axsis (in Miller indice): ' ,I5, 1x, I5, 1x, I5)") m_AXSISZ(1:3)
            if(CtrlParam%NEEDDAMP .gt. 0) then
               write(hFile,FMT="(' !    Quickdampping to be performed before CSP calculation: ', I7, ' time steps')") CtrlParam%NEEDDAMP
            else
               write(hFile,FMT="(' !    Quickdampping to be performed before CSP calculation: NO')")
            end if
            write(hFile,FMT="(' !    ')")

            if(m_nCPTRETH>0) then
               write(hFile,FMT="(' !    atoms to be group by CSP into ',I5, ' groups')") m_nCPTRETH
               write(hFile,FMT="(' !    with the level values:',32((1x,F9.4),','))") m_CPTRETH(1:m_nCPTRETH)
            else
               write(hFile,FMT="(' !    atoms not to be group by CSP')")
            end if

            if(m_Output_AverageCP > 0) then
               write(hFile,FMT="(' !    output average atomic CSP for every ',I5, ' time steps')") m_Output_AverageCP
            else
               write(hFile,FMT="(' !    no average atomic CSP to be output')")
            end if
            write(hFile,FMT="(' !    ')")

           return
 end subroutine PrintControlParameters
  !****************************************************************************


  !**********************************************************************************
  subroutine INITIALIZE_CentroSymmetry_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to allocate and initialize the device memories for device calucaltion
  !               of cnetrosymmetry
  !     INPUT:
  !     OUTPUT:
  !     NOTE:   intialization of this module must be called after
  !             calling Initialize_Globle_Variables_DEV, with the number
  !             of particle on each device is avaliable
  !
  use MD_Globle_Variables_GPU, only:Check_DEVICES
  use MD_TYPEDEF_PrintList,    only:Add_PrintProcess
  implicit none

     !--- dummy vaiables
     type(SimMDBox),  intent(inout)::SimBox
     type(SimMDCtrl), intent(in)   ::CtrlParam
     !--- local varaibles
     integer::ERR,CURDEV, I,IFILE
     real(KINDDF)::RMAT(3,3) = 0

          !$$--- to check device version is used
           call Check_DEVICES()


            if(m_INITED .gt. 0) then
               call Clear_CentroSymmetry_DEV(SimBox, CtrlParam)
            end if

            !$$--- to findout the I/O unit
            IFILE = 0
            do I=1, size(CtrlParam%f_tag)
               if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
                  IFILE = I
               end if
               if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                  m_OUTFILE = CtrlParam%f_others(I)
               end if
            end do

            if(IFILE .le.0 ) then
               write(*,*) "MDPSCU Error: no control file for centro-symmetry-parameter calculation is given."
               write(*,*) "              add the keyword in SETUP file:",  mp_FTAGI
               write(*,fmt="(' Process to be stopped')")
               stop
            end if

           write(*,fmt="(A)") " !**** Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
           if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
            call LoadControlParameters(CtrlParam%f_others(IFILE), SimBox, CtrlParam)

           !$$--- check consisten of input parameters
            if(m_AXSISX(1)*m_AXSISY(1) + m_AXSISX(2)*m_AXSISY(2) + m_AXSISX(3)*m_AXSISY(3) .ne. 0 ) then
               write(*,*) "MDPSCU Error: the axsis X and axsis Y are not perpendicular"
               write(*,*) "              check hte input file ",  CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
               write(*,fmt="(' Process to be stopped')")
               stop
            end if
            if(m_AXSISX(1)*m_AXSISZ(1) + m_AXSISX(2)*m_AXSISZ(2) + m_AXSISX(3)*m_AXSISZ(3) .ne. 0 ) then
               write(*,*) "MDPSCU Error: the axsis X and axsis Z are not perpendicular"
               write(*,*) "              check hte input file ",  CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
               write(*,fmt="(' Process to be stopped')")
               stop
            end if
            if(m_AXSISY(1)*m_AXSISZ(1) + m_AXSISY(2)*m_AXSISZ(2) + m_AXSISY(3)*m_AXSISZ(3) .ne. 0 ) then
               write(*,*) "MDPSCU Error: the axsis Y and axsis Z are not perpendicular"
               write(*,*) "              check hte input file ",  CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
               write(*,fmt="(' Process to be stopped')")
               stop
            end if
            call Add_PrintProcess(PrintControlParameters)

           m_numAccum = 0
           RMAT(1,1:3) = m_AXSISX(1:3)/dsqrt(sum( dble(m_AXSISX(1:3)) * dble(m_AXSISX(1:3))))
           RMAT(2,1:3) = m_AXSISY(1:3)/dsqrt(sum( dble(m_AXSISY(1:3)) * dble(m_AXSISY(1:3))))
           RMAT(3,1:3) = m_AXSISZ(1:3)/dsqrt(sum( dble(m_AXSISZ(1:3)) * dble(m_AXSISZ(1:3))))
           do I=1, mp_NPT
              m_RPVECTX(I) = mp_PVECTX(I)*RMAT(1,1) + mp_PVECTY(I)*RMAT(1,2) + mp_PVECTZ(I)*RMAT(1,3)
              m_RPVECTY(I) = mp_PVECTX(I)*RMAT(2,1) + mp_PVECTY(I)*RMAT(2,2) + mp_PVECTZ(I)*RMAT(2,3)
              m_RPVECTZ(I) = mp_PVECTX(I)*RMAT(3,1) + mp_PVECTY(I)*RMAT(3,2) + mp_PVECTZ(I)*RMAT(3,3)
           end do

           m_INITED = IOR(m_INITED,1)

           return
 end subroutine INITIALIZE_CentroSymmetry_DEV
  !****************************************************************************


  !****************************************************************************
  subroutine Clear_CPArray_template(IDEV, CP)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_CPArrays_template
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !
  !     OUPUT:     CP,      the allocated device memory storing the number of neigbore of an atom
  !
      implicit none
      !--- dummy variables
           integer, intent(in)::IDEV
           real(KINDSF), device, dimension(:), allocatable::CP

      !--- Local vairables
           integer::CURDEV,ERR, NPRT

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)
              !$$--- deallocate memory for neighbore list on devices
              if(allocated(CP)) then
                 deallocate(CP, STAT=ERR)
                 if(ERR) goto 100
              end if

             ERR = cudaSetDevice(CURDEV)
             return

   100       ERR = cudaSetDevice(CURDEV)
             write(*,*) "MDPSCU Warning: fail to deallocate centrosymmetry memory on device", IDEV
             call ONWARNING(gm_OnWarning)
             return
  end subroutine Clear_CPArray_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_CentroSymmetry_DEV(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the device memories allocated for
  !               device calucaltion of Neighbore Table
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   use MD_Globle_Variables_GPU
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam


      !--- Local variables
      integer::I

           !$$--- clear the memory allocated on devices
           do I=1, m_NDEVICE
              call Clear_CPArray_template(m_DEVICES(I), dm_CPDEVWS(I)%CP)
           end do
           !$$ also deallocate the memory on host
           if(allocated(m_CP)) deallocate(m_CP)
           if(allocated(m_ACP)) deallocate(m_ACP)
           if(allocated(m_CPTRETH)) deallocate(m_CPTRETH)
           m_nCPTRETH = 0
           m_numAccum = 0
           m_INITED = 0
      RETURN
  end subroutine Clear_CentroSymmetry_DEV
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_CentroSymmetry_DEV(CP)
  !***  PURPOSE:  to copy the centrosymmetry in devices to a host array CP
  !               the indice in CP is the indice of atoms in the whole system
  !
  !
  !     INPUT:
  !     OUTPUT     CP,      the centrosysmetry
  !
  use MD_Globle_Variables_GPU, only:dm_NPRT, hm_GID, COPY_OUT_SHIFT_template, IndexToOriginalBox
                                     


  implicit none
  !----   DUMMY Variables
          real(KINDSF), dimension(:),allocatable::CP

  !----   Local variables
          integer::I, J, pinnedFlag, ERR, CURDEV
          integer::STARTA, ENDA, NA
          real(KINDSF), pinned, dimension(:),allocatable::TCP
  !***

        !---
            ERR = cudaGetDevice(CURDEV)
            !$$--- determine the array size of the neighbore list
            allocate(TCP(dm_NPRT), STAT=ERR, PINNED=pinnedFlag)
            do I=1, m_NDEVICE
               call COPY_OUT_SHIFT_template(I, dm_CPDEVWS(I)%CP, TCP)
            end do 

            !$$**** restore the array in original order
            call SynchronizeDevices()
            call IndexToOriginalBox(TCP, CP)
            deallocate(TCP)

        RETURN
  end subroutine COPYOUT_CentroSymmetry_DEV
  !**********************************************************************************

  !**********************************************************************************
   attributes(global) subroutine Cal_CP_KERNEL(NBPC ,NPART, XP, NC, NAC, IA1th0, IA1th, NCX, NCY, NCZ, BSX, BSY, BSZ, &
                                                               PDX, PDY,PDZ, LATT, CFROM, CTO, mxNAPDEV, CP)
  !***  PURPOSE:  to update the neighbore list of particles in an cell.
  !               Newton's third law is NOT considered . Be careful to use
  !               corresponding force calculation where Newton's 3th should NOT be applied
  !
  !$$   INPUT:    NBPC,     the number of blocks needed for covering all particles in a cell
  !$$             NPART,    the total number of particles in the whole simulation box
  !$$             XP,       the position of the particles
  !$$             NC,       the number of cells in the whole box
  !$$             NAC,      the number of particles in the cells on the device
  !$$             IA1th0,   the index of the first particle on the device
  !$$             IA1th,    the index of the first particle in a cell
  !$$             NCX, NCY, NCZ, the cell number in X, Y, Z direction
  !$$             BSX, BSY, BSZ, the size of the box in X, Y, Z, dierection
  !$$             LATT,          the lattice length
  !$$             PDX, PDY,PDZ,  the perodic condition in X, Y, Z, dierection
  !$$             CFROM,     the start cell on the device
  !$$             CTO,       the end cell on the device
  !$$             mxNAPDEV,  the max number of particles on a device, e.g, the dimension of KVOIS
  !
  !$$     OUTPUT: CP,       the centrosymetry of particles
  !
  implicit none
      !
      !--- DUMMY VARIABLES
      integer, value::NBPC,NC,NCX,NCY,NCZ, PDX, PDY, PDZ, NPART, CFROM, CTO, mxNAPDEV, IA1th0
      real(KINDDF),value::BSX, BSY,BSZ, LATT
      real(KINDDF), device::XP(NPART,3)
      integer, device::NAC(NC), IA1th(NC)
      real(KINDSF),device::CP(mxNAPDEV)

      !--- Local variables
         !nonshared by threads
         real(KINDSF)::PTX(mp_NPT), PTY(mp_NPT), PTZ(mp_NPT), RR(mp_NPT), &
                       VTX(mp_NPT), VTY(mp_NPT), VTZ(mp_NPT),             &
                       TRR, SEPX, SEPY, SEPZ,POSX, POSY, POSZ
         integer::IB, IB0, IT, IA, IA0,IA00, JA, I, J, K, L, IP0

         !variables share by all thread
         integer,shared::NB, IC, IS0, STARTA, NCXY, NCXYZ,IX0, IY0, IZ0,IC0, NACC0
         integer,shared::NS,NACC, IAC, IACE, FROM, TO
         integer, shared, dimension(mp_NNC)::CID, IX,IY, IZ, OUT

         real(KINDSF), shared, dimension(mp_NNC,3)::CXYZ
         real(KINDSF), shared, dimension(mp_BLOCKSIZE,3)::SPOS


        !$$--- start process
              !$$IB -- the index of block
              !$$
              IB  =  (blockidx%y-1) * griddim%x +  blockidx%x-1

              !$$NB -- the size of a block
              NB = blockdim%x*blockdim%y

              !$$IB0 -- the index of the cells, if the blocksize is smaller than the
              !$$       number of atoms in a cell (NBPC=1), each block cover all atoms
              !$$       in a cell, otherwise, use NBPC block to cover the atoms.
              IB0  =  IB/NBPC

              !$$--- beofre move to the globla index of the cell
              !$$    we get the starting atom in the block
              IP0 = (IB-IB0*NBPC)*NB

              !$$--- move IB0 to the globla index of cell
              IB0  =  IB0 + CFROM-1

              if(IB0 .ge. CTO) return

              !$$IT -- the thread ID
              IT  = (threadidx%y-1)*blockdim%x + threadidx%x

              !$$-- calculate the index of the cell in original box, performed by the first thread
              if(IT .EQ. 1) then

                 NCXY = NCX*NCY
                 NCXYZ = NCXY*NCZ
                 !$$--- IS0, the ID of a box in a multiple-boxs style
                 IS0 = IB0/NCXYZ
                 IC = IB0-IS0*NCXYZ

                 !$$-- 3D index of the original cell in original box
                 IZ0 = IC/NCXY
                 IY0 = (IC-IZ0*NCXY)/NCX
                 IX0 = IC-IZ0*NCXY-IY0*NCX
                 IZ0 = IZ0 + 1
                 IY0 = IY0 + 1
                 IX0 = IX0 + 1

                 !$$IC-- now IC is the index of the cell in the whole box
                 IC  = IB0 + 1

                 !$$STARTA -- the globle index of atom in cell CFROM
                 !$$          it is also the first atom on the device
                 STARTA =  IA1th0 !IA1th(CFROM)

                 !$$NACC0 -- the number of atoms in the fist cell
                 NACC0 = NAC(IC)
              end if
              !$$--- to insure all thread has readed here
              call syncthreads()

              !$$--- if the cell is empty, return
              if(NACC0 .LE. 0) return

              !$$--- to get the position of neighbore cells
              !$$      of the original cell
              if(IT .LE. mp_NNC) then
                 OUT(IT)= 0
                 IZ(IT) = IZ0 + dcm_NIZ(IT)
                 IY(IT) = IY0 + dcm_NIY(IT)
                 IX(IT) = IX0 + dcm_NIX(IT)

                 CXYZ(IT,1) = C_ZERO
                 if(PDX .AND. IT.GT.1) Then
                     if( IX(IT).GT.NCX )THEN
                         IX(IT) = C_IUN
                         CXYZ(IT,1) = BSX
                     else if (IX(IT).LT.C_IUN) THEN
                         IX(IT) = NCX
                         CXYZ(IT,1) = -BSX
                     end if
                  end if

                  CXYZ(IT,2) = C_ZERO
                  if(PDY .AND. IT.GT.1) Then
                     if( IY(IT).GT.NCY )THEN
                         IY(IT) = C_IUN
                         CXYZ(IT,2) = BSY
                     else if (IY(IT).LT.C_IUN) THEN
                         IY(IT) = NCY
                         CXYZ(IT,2) = -BSY
                     end if
                  end if

                  CXYZ(IT,3) = C_ZERO
                  if(PDZ .AND. IT.GT.1) Then
                     if( IZ(IT).GT.NCZ )THEN
                         IZ(IT) = C_IUN
                         CXYZ(IT,3) = BSZ
                     else if (IZ(IT).LT.C_IUN) THEN
                         IZ(IT) = NCZ
                         CXYZ(IT,3) = -BSZ
                     end if
                  end if

                  if( IX(IT) .GT. NCX .OR. IX(IT) .LT. C_IUN) OUT(IT) = 1
                  if( IY(IT) .GT. NCY .OR. IY(IT) .LT. C_IUN) OUT(IT) = 1
                  if( IZ(IT) .GT. NCZ .OR. IZ(IT) .LT. C_IUN) OUT(IT) = 1
                  if(OUT(IT) .EQ. 0) then
                     CID(IT) = NCXY*(IZ(IT)-1)+NCX*(IY(IT)-1)+IX(IT)+IS0*NCXYZ
                  else
                     CID(IT) = 0
                  end if
               end if
               call syncthreads()


               !$$--- we start scanning for atom to calculate their neighbors
               !$$--- IA:the global index of the atom in the original cell.
               !$$--- IA0 is index of the atom on the DEVICE
               !$$--- IA00 is index of the atom in its cell
               !$$
               IA = (IT-1) + IA1th(IC) + IP0
               IA0 = IA - STARTA + 1
               IA00 = IA - IA1th(IC)+1
               POSX =XP(IA, 1)
               POSY =XP(IA, 2)
               POSZ =XP(IA, 3)

              !$$--- calculate the centrosymmetry points
              !$$    initialize the RR to a large value
               do L=1, mp_NPT
                  PTX(L) = POSX +  dcm_RPVECTX(L)*LATT !mp_PVECTX(L)*LATT
                  PTY(L) = POSY +  dcm_RPVECTY(L)*LATT !mp_PVECTY(L)*LATT
                  PTZ(L) = POSZ +  dcm_RPVECTZ(L)*LATT !mp_PVECTZ(L)*LATT
                  RR(L)  = 1.D32
               end do

               !$$--- serach in cell the atom in
               K=1
                  !$$NACC -- the number of atoms in neighboring cell K
                  NACC = NAC(CID(K))
                  !$$NS-- the number of segment of do loop when scanning cell K
                  NS = (NACC-1)/NB+1
                  !$$IAC-- the index of start atom in cell K
                  IAC = IA1th(CID(K))

                  !$$IACE-- the index of end atom in cell K
                  IACE = IAC + NACC -1

                  do J=1, NS
                     FROM = min((J-1)*NB+IAC,IACE)
                     TO = min(FROM+NB-1, IACE)
                     SPOS(IT,1:3) = XP(IT+FROM-1, 1:3)
                     call syncthreads()
                     !$$--- In each block we calculate the neigbores of atoms in an original cell
                     if(IA00.LE.NACC0) then
                        do I=FROM, TO
                              JA = I-FROM+1
                              do L=1, mp_NPT
                                 SEPX = PTX(L) - SPOS(JA,1)
                                 SEPY = PTY(L) - SPOS(JA,2)
                                 SEPZ = PTZ(L) - SPOS(JA,3)
                                 TRR = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                                 if(TRR .LT. RR(L)) then
                                    RR(L) = TRR
                                    VTX(L) = SPOS(JA,1)
                                    VTY(L) = SPOS(JA,2)
                                    VTZ(L) = SPOS(JA,3)
                                 end if
                              end do
                        end do
                     end if
                     call syncthreads()
                  end do

               !$$--- search in neighbor cells
               do K=2, mp_NNC
                  if(OUT(K)) cycle
                  !$$NACC -- the number of atoms in neighboring cell K
                  NACC = NAC(CID(K))

                  !$$NS-- the number of segment of do loop when scanning cell K
                  NS = min((NACC-1)/NB+1, NACC)

                  !$$IAC-- the index of start atom in cell K
                  IAC = IA1th(CID(K))

                  !$$IACE-- the index of end atom in cell K
                  IACE = IAC + NACC -1
                  call syncthreads()

                  do J=1, NS
                     FROM = min((J-1)*NB+IAC,IACE)
                     TO = min(FROM+NB-1, IACE)
                     SPOS(IT,1:3) = XP(IT+FROM-1, 1:3) + CXYZ(K,1:3)
                     call syncthreads()
                     !$$--- In each block we calculate the neigbores of atoms in an original cell
                     if(IA00.LE.NACC0) then
                        do I=FROM, TO
                              JA = I-FROM+1
                              do L=1, mp_NPT
                                 SEPX = PTX(L) - SPOS(JA,1)
                                 SEPY = PTY(L) - SPOS(JA,2)
                                 SEPZ = PTZ(L) - SPOS(JA,3)
                                 TRR = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                                 if(TRR .LT. RR(L)) then
                                    RR(L) = TRR
                                    VTX(L) = SPOS(JA,1)
                                    VTY(L) = SPOS(JA,2)
                                    VTZ(L) = SPOS(JA,3)
                                 end if
                              end do
                        end do
                     end if
                     call syncthreads()
                  end do

               end do
             !$$--- Now we can get the centrosymmetry for the particle
             !$$--- NOTE:we use the summation of absolute value instead of the sqare in its original reference
             !$$--- also NOTE:  IA0 is the index of atoms on the current device
             if(IA00.LE.NACC0) then
                POSX = 2.D0*POSX
                POSY = 2.D0*POSY
                POSZ = 2.D0*POSZ
                CP(IA0) = &
                  (VTX(1) + VTX(2)  -POSX)*(VTX(1) + VTX(2)  -POSX)+ &
                  (VTX(3) + VTX(4)  -POSX)*(VTX(3) + VTX(4)  -POSX)+ &
                  (VTX(5) + VTX(6)  -POSX)*(VTX(5) + VTX(6)  -POSX)+ &
                  (VTX(7) + VTX(8)  -POSX)*(VTX(7) + VTX(8)  -POSX)+ &
                  (VTX(9) + VTX(10) -POSX)*(VTX(9) + VTX(10) -POSX)+ &
                  (VTX(11)+ VTX(12) -POSX)*(VTX(11)+ VTX(12) -POSX)+ &
                  (VTX(13)+ VTX(14) -POSX)*(VTX(13)+ VTX(14) -POSX)+ &

                  (VTY(1) + VTY(2)  -POSY)*(VTY(1) + VTY(2)  -POSY)+ &
                  (VTY(3) + VTY(4)  -POSY)*(VTY(3) + VTY(4)  -POSY)+ &
                  (VTY(5) + VTY(6)  -POSY)*(VTY(5) + VTY(6)  -POSY)+ &
                  (VTY(7) + VTY(8)  -POSY)*(VTY(7) + VTY(8)  -POSY)+ &
                  (VTY(9) + VTY(10) -POSY)*(VTY(9) + VTY(10) -POSY)+ &
                  (VTY(11)+ VTY(12) -POSY)*(VTY(11)+ VTY(12) -POSY)+ &
                  (VTY(13)+ VTY(14) -POSY)*(VTY(13)+ VTY(14) -POSY)+ &

                  (VTZ(1) + VTZ(2)  -POSZ)*(VTZ(1) + VTZ(2)  -POSZ)+ &
                  (VTZ(3) + VTZ(4)  -POSZ)*(VTZ(3) + VTZ(4)  -POSZ)+ &
                  (VTZ(5) + VTZ(6)  -POSZ)*(VTZ(5) + VTZ(6)  -POSZ)+ &
                  (VTZ(7) + VTZ(8)  -POSZ)*(VTZ(7) + VTZ(8)  -POSZ)+ &
                  (VTZ(9) + VTZ(10) -POSZ)*(VTZ(9) + VTZ(10) -POSZ)+ &
                  (VTZ(11)+ VTZ(12) -POSZ)*(VTZ(11)+ VTZ(12) -POSZ)+ &
                  (VTZ(13)+ VTZ(14) -POSZ)*(VTZ(13)+ VTZ(14) -POSZ)
            end if

      return
  end subroutine Cal_CP_KERNEL
  !**********************************************************************************

  !**********************************************************************************
  subroutine StartOnDevice_template(IDEV, CFROM, CTO, NC, NCELL, PD, MXNAC0, BOXSIZE, LATT, NAC, IA1th, XP, CP)
  !***  PURPOSE:  to start calculating centrosymmetry
  !
  !    INPUT(CPU): IDEV,      the index of the device
  !                CFROM,     the start cell on the device
  !                CTO,       the end cell on the device
  !                NC,        the number cell
  !                NCELL,     the number cell in three direction
  !                PD,        if perioedic condition used
  !                MXNAC0,    max number of particles in a cell
  !                BOXSIZE,   the boxsize
  !                LATT,      the lattice length
  !
  !   INPUT(GUP):
  !                NAC,       the number of atoms in the cells
  !                IA1th,     the index of the first atom in a cell
  !                XP,        the position of the atoms
  !
  !   OUTPUT(GPU): CP,        thecentrosymmetry
  !
  use MD_Globle_Variables_GPU

  implicit none
  !----   DUMMY Variables
          integer::IDEV,CFROM, CTO, NC, MXNAC0, NCELL(*), PD(*)
          real(KINDDF)::BOXSIZE(*), LATT

          integer,      device, dimension(:)  ::NAC,IA1th
          real(KINDDF), device, dimension(:,:)::XP
          real(KINDSF), device, dimension(:)  ::CP

  !----   Local variables
         type(dim3) :: blocks
         type(dim3) :: threads
         integer::NS, BX, BY, IP0, ERR, NADEV, CURDEV, C0, C1, NCINB, RESC, STARTA

  !--- start

               !$$--- copy the cell informations into device memeory
               ERR = cudaSetDevice(IDEV)

               !$$--- to determine of a block (the number of threads in a block)
               !$$     mp_BLOCKSIZE must be larger than 27
               BX = min(mp_BLOCKSIZE, MXNAC0)
               BX = max(BX, mp_NNC)
               BY = 1
               NS = (MXNAC0-1)/(BX*BY)+1

              !$$--- copy in the axsis
               dcm_RPVECTX(1:mp_NPT) = m_RPVECTX(1:mp_NPT)
               dcm_RPVECTY(1:mp_NPT) = m_RPVECTY(1:mp_NPT)
               dcm_RPVECTZ(1:mp_NPT) = m_RPVECTZ(1:mp_NPT)

               !$$-- to determine the dimension of blocks( the number of blocks in a grid)
               blocks  = dim3((CTO-CFROM+1), NS, 1)
               threads = dim3(BX, BY, 1)

               !$$-- to determine the dimension of blocks( the number of blocks in a grid)
               !$$   note the max gridsize is 65535*65535*65535
               STARTA = hm_IA1th(CFROM)
               RESC   = CTO-CFROM+1
               NCINB  = min(RESC, 65535)
               C0 = CFROM
               C1 = C0 + NCINB - 1
               do while(.true.)
                  if(C0 .le. 0) exit
                  blocks  = dim3(NCINB, NS, 1)
                  call Cal_CP_KERNEL<<<blocks,threads>>>(NS, dm_NPRT, XP, NC, NAC, STARTA, IA1th,          &
                                                                      NCELL(1), NCELL(2), NCELL(3),        &
                                                                      BOXSIZE(1), BOXSIZE(2), BOXSIZE(3),  &
                                                                      PD(1), PD(2), PD(3), LATT,           &
                                                                      C0, C1, m_NAPDEV,CP)


                  RESC  = RESC - NCINB
                  NCINB = min(RESC, 65535)
                  if(NCINB .le. 0) exit
                  C0 = C1 + 1
                  C1 = C0 + NCINB - 1
              end do


         return

  end subroutine StartOnDevice_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Cal_CentroSymmetry_DEV(SimBox, CtrlParam, CP)
  !***  PURPOSE:  to update the neighbore list of particles with head-link
  !               technique used.
  !
  !               Newton's 3rd law is NOT considered in creating the list
  !               Be careful to use the correct force calculation in which
  !               the 3rd law should NOT be applied
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters

  !                MultiBox,  optional, inidicating if multiple-box used
  !
  !     OUTPUT     CP, optional, the centrosymmetry of particles indexed in
  !                    original box
  !

  use MD_Globle_Variables_GPU
  use MD_NeighborsList_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)                                  ::SimBox
          type(SimMDCtrl)                                 ::CtrlParam
          real(KINDSF), dimension(:),allocatable, optional::CP


  !----   Local variables
         integer::TNC, NCELL(3), MXNAC, PD(3), ERR, CURDEV, I
         real(KINDDF)::BOXSIZE(3), LATT
  !***
          !$$---
          ERR = cudaGetDevice(CURDEV)

          if(iand(m_INITED,2) .eq. 0) then
              !     NOTE:   intialization of this module must be called after
              !             calling Initialize_Globle_Variables_DEV, with the number
              !             of particle on each device is avaliable
             call  Allocate_WorkingArray_DEV()
          end if

          !$$--- get information of the box
          call  GetCellInform(TNC, NCELL, MXNAC)
          BOXSIZE = SimBox%ZL
          PD = CtrlParam%IFPD
          LATT = SimBox%RR

           !$$--- to start neighbor calculation on devices
           do I = 1, m_NDEVICE
              call StartOnDevice_template(m_DEVICES(I), m_STARTCELL(I), m_ENDCELL(I),                                    &
                                          TNC, NCELL, PD, MXNAC, BOXSIZE, LATT,                                          &
                                          dm_WorkSpace%NAC(I)%Data, dm_WorkSpace%IA1th(I)%Data, dm_WorkSpace%XP(I)%Data, &
                                          dm_CPDEVWS(I)%CP                           )
           end do

           !$$--- copy the centrosymmetry
           !$$  - NOTE: CP is allocatable array, if CP has not been
           !$$           allocated, it is to be allocated in COPYOUT routine
           !$$           do not forget deallocate CP
           if(present(CP)) then
               call SynchronizeDevices()
               call COPYOUT_CentroSymmetry_DEV(CP)
               !$$--- normalize CP to lattice length
               CP = CP/(LATT*LATT)
            end if
            !---
           ERR = cudaSetDevice(CURDEV)

           return

        RETURN
  end subroutine Cal_CentroSymmetry_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine Output_Instant_CentroSymmetry(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the centrosymmetry to output file. This routine is to interfaced
  !                  to and MD_EM_TB_ForceTable_SHELL_12_GPU. It is assumed the the neighbor
  !                  list routine has be called
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_EM_TB_ForceTable_SHELL_12_GPU
  use MD_Globle_Variables_GPU
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE
       integer::hFile, I, K, IP, TY0, NTY, L, L0, JOB

       !*** to output instant CP
           JOB = Stamp%ITest
           IP  = Stamp%IRec(1)

           call STRCATI(GFILE, m_OUTFILE, "P", m_processid, 4)
           call STRCATI(GFILE, GFILE, "_", JOB, 4)
           call STRCATI(GFILE, GFILE, ".", IP, 4)

           call AvailableIOUnit(hFile)

           write(*,*) "Instant centro-symmetry-parameters  to be saved in "//GFILE(1:len_trim(GFILE))
           open(UNIT=hFile, file = GFILE, status='unknown')
             write(hFile, fmt="(A)") "!--- CENTROSYMETRY-FCC RESULTS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
             write(hFile, fmt="(A)") "!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY"
             write(hFile, fmt="(A)") '!   AUTHOR: HOU Qing'
             write(hFile, fmt="(A)") '!    '
             write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
             write(hFile, fmt="(A,1X,I8)")            "&NATOM      ", dm_NPRT
             write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox(1)%ZL/SimBox(1)%RR
             write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox(1)%BOXLOW/SimBox(1)%RR
             write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox(1)%RR*CP_CM2A
             write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL    ", 1

             if(m_nCPTRETH > 0) then
                write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL     ", 3, 4, 5
                write(hFile, fmt="(A,1X,3(I4,1X))")      "&CSPCOL      ", 6
                write(hFile, fmt="('!---',A6, 2x, A8, 2x, 3(A13,2x),2x,A13,2x)")  &
                             "TYPE", "CP-GROUP", "X(latt.)", "Y(latt.)", "Z(latt.)",  "Centro-Param"
             else
                write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL     ", 2, 3, 4
                write(hFile, fmt="(A,1X,3(I4,1X))")      "&CSPCOL     ", 5
                write(hFile, fmt="('!---',A6, 2x, 3(A13,2x),2x,A13,2x)")          &
                             "TYPE", "X(latt.)", "Y(latt.)", "Z(latt.)",  "Centro-Param"
             endif

             IP = 0
             do I=1, size(SimBox)
                L0 = SimBox(I)%NGROUP
                do K=1, SimBox(I)%NPRT
                  IP = IP + 1
                  TY0 = SimBox(I)%ITYP(K)
                  if(m_nCPTRETH > 0) then
                     NTY = TY0
                     do L=1, m_nCPTRETH
                        if(m_CP(IP) .gt. m_CPTRETH(L)) NTY = L + L0
                     end do
                     write(hFile,fmt="(4x, I4,2X,I8,6x, 6(1PE13.4,2X),I6,2X,15(1PE13.4,1X))")&
                                     TY0, NTY, SimBox(I)%XP(K,1:3)/SimBox(I)%RR, m_CP(IP)
                  else
                     write(hFile,fmt="(4x, I4,2X,6(1PE13.4,2X),I6,2X,15(1PE13.4,1X))") &
                                     TY0, SimBox(I)%XP(K,1:3)/SimBox(I)%RR, m_CP(IP)
                  end if
                end do
             end do
          close(hFile)
          return
  end subroutine Output_Instant_CentroSymmetry
  !****************************************************************************************

  !****************************************************************************************
  subroutine Output_Average_CentroSymmetry(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the centrosymmetry to output file.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
       character*256::GFILE
       integer::hFile, I, K, IP, TY0, NTY, L, L0, JOB,ICFG
       !----

           JOB    = Stamp%ITest
           ICFG   = Stamp%IRec(1)

           IP = 0
           do I=1, size(SimBox)
              do K=1, SimBox(I)%NPRT
                 IP = IP + 1
                 m_ACP(IP) = m_ACP(IP) + m_CP(IP)
              end do
           end do
           m_numAccum = m_numAccum + 1
           if(mod(m_numAccum,m_Output_AverageCP) .ne. 0) return
              m_ACP  = m_ACP/dble(m_numAccum)
          !*** to output instant CP
           call STRCATI(GFILE, m_OUTFILE, "_Average_P", m_processid, 4)
           call STRCATI(GFILE, GFILE, "_", JOB, 4)
           call STRCATI(GFILE, GFILE, ".", ICFG, 4)
           call AvailableIOUnit(hFile)

           write(*,*) "Average central-symmetry-parameters  to be saved in "//GFILE(1:len_trim(GFILE))
           open(UNIT=hFile, file = GFILE, status='unknown')
           if(m_nCPTRETH > 0) then
              write(hFile, fmt="('!---',A6, 2x, A8, 2x, 3(A13,2x),2x,A13,2x)")  &
                                 "TYPE", "CP-GROUP", "X(latt.)", "Y(latt.)", "Z(latt.)",  "Central-Param"
           else
              write(hFile, fmt="('!---',A6, 2x, 3(A13,2x),2x,A13,2x)")  &
                             "TYPE", "X(latt.)", "Y(latt.)", "Z(latt.)",  "Central-Param"
           endif

           IP = 0
           do I=1, size(SimBox)
              L0 = SimBox(I)%NGROUP
              do K=1, SimBox(I)%NPRT
                IP = IP + 1
                TY0 = SimBox(I)%ITYP(K)
                if(m_nCPTRETH > 0) then
                   NTY = TY0
                   do L=1, m_nCPTRETH
                     if(m_ACP(IP) .gt. m_CPTRETH(L)) NTY = L + L0
                   end do
                   WRITE(hFile,fmt="(4x, I4,2X,I8,6x, 6(1PE13.4,2X),I6,2X,15(1PE13.4,1X))")&
                                     TY0, NTY, SimBox(I)%XP(K,1:3)/SimBox(I)%RR, m_ACP(IP)
                 else
                   WRITE(hFile,fmt="(4x, I4,2X,6(1PE13.4,2X),I6,2X,15(1PE13.4,1X))") &
                                     TY0, SimBox(I)%XP(K,1:3)/SimBox(I)%RR, m_ACP(IP)
                end if
             end do
          end do
          close(hFile)
          m_numAccum = 0
          m_ACP = 0.D0

          return
  end subroutine Output_Average_CentroSymmetry
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_CentroSymmetry(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the centrosymmetry to output file. This routine is to interfaced
  !                  to and MD_EM_TB_ForceTable_SHELL_12_GPU. It is assumed the the neighbor
  !                  list routine has be called
  !
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_AppShell_16_GPU.f90
  !            MD_SimBoxArray_ToolShell_14_GPU.f90

  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables

           call Cal_CentroSymmetry_DEV(SimBox(1), CtrlParam, m_CP)

           !*** to accumulating CP
           if(m_Output_AverageCP > 0) then
              call Output_Average_CentroSymmetry(Stamp, SimBox, CtrlParam)
           end if
           call Output_Instant_CentroSymmetry(Stamp, SimBox, CtrlParam)
          return
  end subroutine RECORD_CentroSymmetry
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_CentroSymmetryTool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the centrosymmetry to output file.
  !                  comapred to RECORD_CentroSymmetry. This routine is to interfaced
  !                  to MD_SimBoxArray_ToolShell, in which the the neighbor
  !                  list routine has NOT be called. Thus we should call it here.
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_AppShell_16_GPU.f90
  !            MD_SimBoxArray_ToolShell_14_GPU.f90
  !
  use MD_NeighborsList_GPU
  use MD_SimboxArray_GPU
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables


            if(Stamp%ITime .LT. 0) return
            call RECORD_CentroSymmetry(Stamp, SimBox, CtrlParam)

            return
  end subroutine RECORD_CentroSymmetryTool
  !****************************************************************************************
  end module CentroSymParam_for_FCC_GPU





