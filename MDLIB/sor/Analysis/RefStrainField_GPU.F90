  module RefStrainField_GPU
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to calculate the strain and displacement vectors of atoms related to a
  !                  reference configuration. The strain vector of an atom is defined as the position
  !                  vector of the atom related to the nearest reference atom of the same type.
  !
  !
  !                  DEPENDENCE____________________________________________________________________________
  !                       MD_NeighborsList_GPU.cuf
  !                       MD_Globle_Variables_GPU.cuf
  !
  !                  SEE ALSO____________________________________________________________________________
  !                       RefVoronoiVA_GPU.F90
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2015-06 (Hou Qing)
  !                  version 2st 2017-05 (Hou Qing)
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
      type(InputStatements)::m_INPUTS
      character(len=17),parameter, private::mp_FTAGI="&AUXF_REFSTRIANIN"
      character(len=18),parameter, private::mp_FTAGO="&AUXF_REFSTRIANOUT"
      character(len=256),  private::m_OUTFILE =""                      ! filename of output data


      !--- C2050
      !integer, parameter, private::mp_MXSHARESIZE = 3*3*224
      !integer,parameter::BLOCKSIZE = 512, BLOCKDIMX=32
      !--- C1060
      integer, parameter, private::mp_MXSHARESIZE = 3*3*64
      integer, parameter, private::mp_BLOCKSIZE = 256, mp_BLOCKDIMX=128

      !--- calculation control parameters
      integer, private::m_INITED = 0

      !$$--- the filename of reference configuration.
      character*256, private::m_RefCfgFile= ""
      real(KINDDF), private::hm_BOND0 = 1.5D0                       !-- the bond-length in lattice unitused for neighbore calculation of reference system


      !--- arrayies on host for the refernce system
      integer, private::hm_NPRT0                                    !-- the number of particles of template system
                                                                    !   NOTE: this number could be different from the
                                                                    !         number of particles in relaxing system
                                                                    !         For exmaple, in cases that particles
                                                                    !         are continuneously incident
      integer, private::hm_NC0                                      !-- the number of cells in reference box
      integer, private::hm_MXNAC0                                   !-- the max number of reference particles in a cell
      integer, private::hm_NCELL0(3)                                !-- the number of cells in 3 dimension for reference box
      integer, private::hm_PD0(3)                                   !-- the periodic condition
      real(KINDDF), private::hm_BOXSIZE0(3)                         !-- the bosxsize of the reference box
      real(KINDDF), private::hm_BOXLOW0(3)                          !-- the low-limit of the reference box

      integer, dimension(:), allocatable, private::hm_NAC0
      integer, dimension(:), allocatable, private::hm_IA1th0

      integer, dimension(:),         allocatable, private::hm_ITYP0 !-- the atom type of reference lattice
      real(KINDDF), dimension(:,:),  allocatable, private::hm_XP0   !-- the positions of reference lattice after partitioned
                                                                    !    hm_GID0 is created by neighbor list calculation.

      !--- position of atoms in the (reference) first system
      !    NOTE: the positions on device have been sorted
      !--- Workspace for reference box
      private:: RefStrain_DEVWS
      type::RefStrain_DEVWS    
         !--- position of atoms in the (reference) first system
         !    NOTE: the positions on device have been sorted
          real(KINDDF), device, dimension(:,:), allocatable::XP0
          integer,      device, dimension(:),   allocatable::ITYP0
          !--- the linked cell information
          integer,      device, dimension(:),   allocatable::NAC0
          integer,      device, dimension(:),   allocatable::IA1th0
          !--- the displacement vectors
          real(KINDSF), device, dimension(:,:), allocatable,private::DV
          !--- the site id of the atom cloest to
          integer,      device, dimension(:),   allocatable::SID        
      end type RefStrain_DEVWS
      type(RefStrain_DEVWS), dimension(m_MXDEVICE), private::dm_RefStrain

       !---parameters used in calculating
       integer, parameter, private::mp_NNC=27
       integer,      constant,private::dcm_NIX(mp_NNC)
       integer,      constant,private::dcm_NIY(mp_NNC)
       integer,      constant,private::dcm_NIZ(mp_NNC)

      !--- some private routines
      private:: Allocate_WorkingArray0_template, &
                Allocate_WorkingArray0_DEV,      &
                Allocate_WorkingArray_DEV,       &
                LoadControlParameters,           &
                Cal_StrainField_KERNEL,       &
                StartOnDevice_template,          &
                COPYOUT_StrainField_template

  contains

 !****************************************************************************
  subroutine Allocate_WorkingArray0_template(IDEV, ITYPR, XPR, NCA, IA1th)
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
           integer,      device, dimension(:),   allocatable::ITYPR
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
              allocate(ITYPR(hm_NPRT0), XPR(hm_NPRT0,3), NCA(hm_NC0), IA1th(hm_NC0), STAT=err)

              if(ERR) then
                 write(*,*) "MDPSCU ERROR: Fail to allocate memory in RefDisplace module on device", IDEV
                 stop
              end if

              NCA   = hm_NAC0
              IA1th = hm_IA1th0
              XPR   = hm_XP0
              ITYPR = hm_ITYP0

              dcm_NIX = NIX
              dcm_NIY = NIY
              dcm_NIZ = NIZ

              ERR = cudaSetDevice(CURDEV)
              return
  end subroutine Allocate_WorkingArray0_template
  !**********************************************************************************

  !****************************************************************************
  subroutine Allocate_WorkingArray_template(IDEV, NAPDEV, DV,SID)
  !***  PURPOSE:   to allocate working space for a device
  !
  !     INPUT:     IDEV,     the index of device on which the memory to be allocated
  !                NAPDEV,   the max number of atoms that a devce treats
  !
  !     OUPUT:     DV,       the allocated device memory storing the displacement vector of atoms
  !                SID,      the allocated device memory storing the ID of most closest sites
      implicit none
      !--- dummy variables
           integer, intent(in)                              ::IDEV, NAPDEV
           real(KINDSF), device, dimension(:,:), allocatable::DV
           integer,      device, dimension(:),   allocatable::SID

      !--- Local vairables
           integer::ERR, NPRT, CURDEV


              !$$--- allocate memory for neighbore list on devices
              ERR = cudaGetDevice(CURDEV)
              ERR = cudaSetDevice(IDEV)
              allocate(DV(NAPDEV, 3), SID(NAPDEV), STAT=err)

              if(ERR) then
                 write(*,*) "MDPSCU ERROR: Fail to allocate memory in RefDisplace module on device", IDEV
                 stop
              end if

              DV  = 0.D0
              SID = 0
              ERR = cudaSetDevice(CURDEV)
              return
  end subroutine Allocate_WorkingArray_template
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
  !             and also one call to Cal_NeighBoreList2C_DEV
  !             with the information of cells of the template
  !             system is available.
  !
  use MD_Globle_Variables_GPU
  use MD_NeighborsList_GPU

  implicit none
  integer::ERR, I

           hm_NPRT0 = dm_NPRT
           call GetCellInform(hm_NC0, hm_NCELL0, hm_MXNAC0)

           allocate(hm_ITYP0(hm_NPRT0), hm_XP0(hm_NPRT0,3), hm_NAC0(hm_NC0), hm_IA1th0(hm_NC0), STAT=ERR)
           if(err .gt. 0) then
              write(*,fmt="(A)") "MDPSCU Error: fail to allocated working space in ReVoronoiVA"
              stop
           end if
           !$$--- save the current statu of globale variable on devices
           !$$    NOTE: hm_XP are the sorted positions of reference reference atoms
           !$$
           hm_NAC0   = hm_NAC
           hm_IA1th0 = hm_IA1th
           hm_XP0    = hm_XP
           hm_ITYP0  = hm_ITYP

           !$$--- allocate the memory allocated on device 1-6
           do I=1, m_NDEVICE
              call Allocate_WorkingArray0_template(m_DEVICES(I), dm_RefStrain(I)%ITYP0, dm_RefStrain(I)%XP0, &
                                                    dm_RefStrain(I)%NAC0, dm_RefStrain(I)%IA1th0)
           end do 
           m_INITED = IOR(m_INITED,1)
           return
  end subroutine Allocate_WorkingArray0_DEV
  !****************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray_DEV()
  !***  PURPOSE:  to allocate and initialize the device memories for SITES
  !
  !     INPUT:
  !     OUTPUT:
  !     NOTE:   intialization of this module must be called after
  !             calling Initialize_Globle_Variables_DEV, with the number
  !             of particle on each device is avaliable
  !             and also one call to Cal_NeighBoreList2C_DEV
  !             with the information of cells of the template
  !             system is available.
  !
  use MD_Globle_Variables_GPU, only:m_NAPDEV

  implicit none
  integer::I

           do I=1, m_NDEVICE
              call Allocate_WorkingArray_template(m_DEVICES(I), m_NAPDEV, dm_RefStrain(I)%DV, dm_RefStrain(I)%SID)
           end do 
           m_INITED = IOR(m_INITED,2)
           return;
  end subroutine Allocate_WorkingArray_DEV
  !****************************************************************************

  !****************************************************************************
  subroutine Clear_WorkingArray_template(IDEV, RTYP, RXP, DV, SID, NAC,IA1th)
  !***  PURPOSE:   to deallocate device memories allocated in calling
  !                Allocate_WorkingArray_template
  !
  !     INPUT:     IDEV,    the index of device on which the memory to be allocated
  !     OUPUT:
  !
      implicit none
      !--- dummy variables
           integer,                               intent(in)::IDEV
           integer,      device, dimension(:),   allocatable::RTYP
           real(KINDDF), device, dimension(:,:), allocatable::RXP
           real(KINDSF), device, dimension(:,:), allocatable::DV
           integer,      device, dimension(:),   allocatable::SID
           integer,      device, dimension(:),   allocatable::NAC, IA1th

      !--- Local vairables
           integer::CURDEV,ERR

              !$$--- get the current device
              ERR = cudaGetDevice(CURDEV)

              ERR = cudaSetDevice(IDEV)
              !$$--- deallocate memory for neighbore list on devices
              if(allocated(RTYP)) then
                 deallocate(RTYP, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(RXP)) then
                 deallocate(RXP, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(DV)) then
                 deallocate(DV, STAT=ERR)
                 if(ERR) goto 100
              end if

              if(allocated(SID)) then
                 deallocate(SID, STAT=ERR)
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
  end subroutine Clear_WorkingArray_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Clear_WorkingArray_DEV(SimBox, CtrlParam)
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

           !$$--- clear the memory allocated on device 1
           do I=1, m_NDEVICE
              call Clear_WorkingArray_template(m_DEVICES(I), dm_RefStrain(I)%ITYP0, dm_RefStrain(I)%XP0, &
                                                             dm_RefStrain(I)%DV,    dm_RefStrain(I)%SID,  &
                                                             dm_RefStrain(I)%NAC0,  dm_RefStrain(I)%IA1th0)
           end do 
           if(allocated(hm_NAC0))   deallocate(hm_NAC0)
           if(allocated(hm_IA1th0)) deallocate(hm_IA1th0)
           if(allocated(hm_ITYP0))  deallocate(hm_ITYP0)
           if(allocated(hm_XP0))    deallocate(hm_XP0)

           m_INITED = 0
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
  use MiniUtilities

  implicit none
  !--- dummy vaiables
  character*(*)  ::fname
  type(SimMDBox) ::SimBox
  type(SimMDCtrl)::CtrlParam

  !--- local variables
  integer::hFile, N, I, LINE, NEEDDAMP, DAMPSCHEME
  character*256::STR
  character*32::STRNUMB(10), KEYWORD

            !$$--- load input statements
            call Load_InputStatements(fname, m_INPUTS, mp_FTAGI)

            !$$--- start interprating the statements
             call Get_InputStatements("&REFCFG", m_INPUTS, STR, LINE)
             if(LINE .gt. 0) then
               m_RefCfgFile = ""
               call Extract_Substr(STR,1,n,m_RefCfgFile)
             end if

             !*** get the controal paraemter of output average CP values
             call Get_InputStatements("&QUICKDAMP", m_INPUTS, STR, LINE)
             if(LINE .eq. 0) call Get_InputStatements("&QUICKDAMP", m_INPUTS, STR, LINE)
             if(LINE .gt. 0) then
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
              end if

             !*** To get job range to be analysis
             call Get_InputStatements("&JOBSEL", m_INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,3,n,STRNUMB)
                CtrlParam%JOBID0 = ISTR(STRNUMB(1))
                if(N .GE. 2) CtrlParam%JOBID1  = ISTR(STRNUMB(2))
                if(N .GE. 3) CtrlParam%JOBIDSTEP = ISTR(STRNUMB(3))
             end if

             !$$*** To get cfg range to be analysis
             call Get_InputStatements("&CFGSEL", m_INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,3,n,STRNUMB)
                CtrlParam%STARTCFG = ISTR(STRNUMB(1))
                if(N .GE. 2) CtrlParam%ENDCFG  = ISTR(STRNUMB(2))
                if(N .GE. 3) CtrlParam%CFGSTEP = ISTR(STRNUMB(3))
             end if

             !$$*** To get cfg range to be analysis
             call Get_InputStatements("&BOXSEL", m_INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                call Extract_Numb(STR,3,n,STRNUMB)
                CtrlParam%STARTBOX = ISTR(STRNUMB(1))
                if(N .GE. 2) CtrlParam%ENDBOX  = ISTR(STRNUMB(2))
                if(N .GE. 3) CtrlParam%BOXSTEP = ISTR(STRNUMB(3))
             end if

             !$$*** To get the output path
             call Get_InputStatements(mp_FTAGO, m_INPUTS, STR, LINE)
             if(LINE .gt. 0) then
                 call Extract_Substr(STR,1,n,m_OUTFILE)
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
            return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine PrintControlParameters(hFile, SimBox, CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
  !----   DUMMY Variables
   integer,         intent(in)::hFile
   type(SimMDBox),  intent(in)::SimBox
   type(SimMDCtrl), intent(in)::CtrlParam

  !----   Local variables

         !$$--- print out the control parameters
          write(hFile,*) "!************ RefStrainField_GPU module to be performed **********"
          write(hFile,*) "!    With the control parameters:"
          if(len_trim(m_RefCfgFile) .gt. 0) then
             write(hFile,FMT="(' !    Reference configuration to be loaded from: ', A)") m_RefCfgFile(1:len_trim(m_RefCfgFile))
          else
             write(hFile,FMT="(' !    Reference configuration is the initial configuration loaded from: ', A)") &
                                                                                 SimBox%IniConfig(1:len_trim(SimBox%IniConfig))
          end if
          write(hFile,fmt="(A)")  " !     "

          if(CtrlParam%NEEDDAMP .gt. 0) then
             write(hFile,FMT="(' !    Quickdampping to be performed for ', I7, ' time steps before RefDisplace')") CtrlParam%NEEDDAMP
          else
             write(hFile,FMT="(' !    Quickdampping to be performed before RefDisplace calculation: NO')")
          end if
          write(hFile,FMT="(' !    ')")

         return
   end subroutine PrintControlParameters
  !**********************************************************************************

  !**********************************************************************************
   subroutine Initialize(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC
   use MD_Globle_Variables_GPU
   use MD_NeighborsList_GPU
   use MD_TYPEDEF_PrintList,    only:Add_PrintProcess
   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   type(SimMDCtrl)::tCtrlParam
   integer::I, IFILE


        !$$--- to check device version is used
          call Check_DEVICES()

        !$$---  we create a temperaral simulation box and controal parameter
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
            write(*,fmt="(' MDPSCU Error: the control file is missed in RefStrainField_GPU module')")
            write(*,*) "                add the keyword in SETUP file:", mp_FTAGI
            write(*,fmt="(' Process to be stopped')")
            stop
         end if

         write(*,fmt="(A)") "Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
         if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") "Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
         call LoadControlParameters(tCtrlParam%f_others(IFILE), SimBox, CtrlParam)
         call Add_PrintProcess(PrintControlParameters)

         !$$--- load the reference configuration
          write(*,FMT="(' Load reference configuration...', A)")

          !$$--- get box size so on
          call ClearDataProKWD_SimMDBox(SimBox)
          call AddDataProKWD_SimMDBox(SimBox,"XYZCOL")

          if(len_trim(m_RefCfgFile) .le. 0) then
            call Initialize_Config_SimMDBox(SimBox, fname=SimBox%IniConfig)
          else
            call Read_Initial_Config_SimMDBox(SimBox, fname=m_RefCfgFile, fmt=0, mod=1)
          end if

         !$$--- RM is needed for neighbore calculations
          tCtrlParam%NB_RM = hm_BOND0*SimBox%RR

         !$$--- genereate the reference box
          hm_PD0      = CtrlParam%IFPD
          hm_BOXSIZE0 = SimBox%ZL
          hm_BOXLOW0  = SimBox%BOXLOW

          !$$*** to initialize the GPU variable
          call Initialize_Globle_Variables_DEV(SimBox, tCtrlParam)

          !$$*** to paitition the system
          call Initialize_NeighboreList_DEV(SimBox, tCtrlParam)
          call Cal_NeighBoreList_DEV(SimBox, tCtrlParam)
          call Allocate_WorkingArray0_DEV()

         return
   end subroutine Initialize
  !**********************************************************************************

  !**********************************************************************************
   attributes(global) subroutine Cal_StrainField_KERNEL(NPART0, ITYP0, XP0, NC0, NAC0, IA1th0, NCX0, NCY0, NCZ0,  &
                                                                LBX0, LBY0, LBZ0, BSX0, BSY0, BSZ0, PDX, PDY, PDZ,   &
                                                                FROMA, TOA, NPART, ITYP, XP,                         &
                                                                DV, SID)
  !$$***PURPOSE:  to identify the Voronoi volume that the atoms in. If an atom is in ith
  !$$             Voronoi volume, the AC of this Voronoi is accumulates 1
  !$$
  !$$
  !$$   INPUT:    NPART0,   the total number of particles in the whole ref.box
  !$$             ITYP0,    the type of atoms
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
  !$$             ITYP,     the type of the particles in the real box
  !$$             XP,       the position of the particles in the real box
  !$$
  !$$     OUTPUT: DV,       the displacement vector of the atoms
  !$$             SID,      the ID of site the the atoms closest
  !
  implicit none
      !
      !--- DUMMY VARIABLES
      integer, value::NPART0,NC0,NCX0,NCY0,NCZ0, PDX, PDY, PDZ, FROMA, TOA
      real(KINDDF),value::LBX0, LBY0, LBZ0, BSX0, BSY0, BSZ0
      integer, device::ITYP0(NPART)
      real(KINDDF), device::XP0(NPART0,3)
      integer, device::NAC0(NC0), IA1th0(NC0)

      integer, value::NPART
      integer, device::ITYP(*)
      real(KINDDF), device::XP(NPART,3)
      real(KINDSF),device::DV(NPART,3)
      integer, device::SID(*)

      !--- Local variables
         real(KINDSF), parameter::eps=0.0001
         real(KINDSF)::SEPX, SEPY, SEPZ, POSX, POSY, POSZ, CX, CY, CZ, RRMI, RR
         integer::IB, NB, IT, IA, K, JA, IMI
         integer::IX0, IY0, IZ0,IC0
         integer::NACC, IAC, IACE, CID, IX,IY, IZ, TYP

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
              TYP  = ITYP(IA)
              DV(IA-FROMA+1,1:3) = 0.D0
              SID(IA-FROMA+1)    = 0

             !$$-- if not in the reference box, return
              if((POSX .LT. LBX0 .OR. POSX-LBX0 .GT. BSX0) .and. .not.PDX ) return;
              if((POSY .LT. LBY0 .OR. POSY-LBY0 .GT. BSY0) .and. .not.PDY ) return;
              if((POSZ .LT. LBZ0 .OR. POSZ-LBZ0 .GT. BSZ0) .and. .not.PDZ ) return;

             !$$-- determine which cell of the reference box that the IAth atom in
              IX0 = int( (POSX - LBX0)/BSX0*dble(NCX0) -eps)+1
              IY0 = int( (POSY - LBY0)/BSY0*dble(NCY0) -eps)+1
              IZ0 = int( (POSZ - LBZ0)/BSZ0*dble(NCZ0) -eps)+1

             !$$--- we start scanning the potential atoms closest to the atom
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
                     if(TYP .eq. ITYP0(JA)) then
                        SEPX = POSX - XP0(JA,1)
                        SEPY = POSY - XP0(JA,2)
                        SEPZ = POSZ - XP0(JA,3)
                        RR = SEPX*SEPX + SEPY*SEPY + SEPZ*SEPZ
                        if(RR .LT. RRMI) then
                           RRMI = RR
                           IMI  = JA
                        end if
                     end if
                  end do
               end do

            !$$--- Now we can get the displacement vector the particle occupying
             if(IMI .gt. 0) then
                DV(IA-FROMA+1,1)  = POSX - XP0(IMI, 1)
                DV(IA-FROMA+1,2)  = POSY - XP0(IMI, 2)
                DV(IA-FROMA+1,3)  = POSZ - XP0(IMI, 3)
                SID(IA-FROMA+1)   = IMI
             end if

      return
  end subroutine Cal_StrainField_KERNEL
  !**********************************************************************************

  !**********************************************************************************
  subroutine StartOnDevice_template(IDEV, NAC0, IA1th0, ITYP0, XP0, ITYP, XP, DV,SID)
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
          integer,      device, dimension(:)   ::NAC0,IA1th0
          real(KINDDF), device, dimension(:,:) ::XP0, XP
          integer,      device, dimension(:)   ::ITYP0, ITYP
          real(KINDSF), device, dimension(:,:) ::DV
          integer,      device, dimension(:)   ::SID

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
                call Cal_StrainField_KERNEL<<<blocks,threads>>>(hm_NPRT0, ITYP0, XP0, hm_NC0, NAC0, IA1th0,         &
                                                                        hm_NCELL0(1),  hm_NCELL0(2),   hm_NCELL0(3),   &
                                                                        hm_BOXLOW0(1), hm_BOXLOW0(2),  hm_BOXLOW0(3),  &
                                                                        hm_BOXSIZE0(1),hm_BOXSIZE0(2), hm_BOXSIZE0(3), &
                                                                        hm_PD0(1),     hm_PD0(2),      hm_PD0(3),      &
                                                                        STARTA, ENDA, dm_NPRT, ITYP, XP, DV, SID)

               ERR = cudaSetDevice(CURDEV)
         return

  end subroutine StartOnDevice_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine Cal_StrainField_DEV(SimBox, CtrlParam, hDV, hSID)
  !***  PURPOSE:  to calculate the site ID that atoms occupy
  !
  !     INPUT:     SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !
  !     OUTPUT     DV, optional, the displacement vector
  !

  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
          type(SimMDBox)                       ::SimBox
          type(SimMDCtrl)                      ::CtrlParam
          real(KINDSF), dimension(:,:),optional::hDV
          integer,      dimension(:),  optional::hSID

  !----   Local variables
         integer::I
  !***

          if(iand(m_INITED,2) .eq. 0) then
             call  Allocate_WorkingArray_DEV()
          end if

           !$$--- to start neighbor calculation on devices
             do I=1, m_NDEVICE 
                 call StartOnDevice_template(I, dm_RefStrain(I)%NAC0,       dm_RefStrain(I)%IA1th0,  &
                                                dm_RefStrain(I)%ITYP0,      dm_RefStrain(I)%XP0,     &
                                                dm_WorkSpace%ITYP(I)%Data,  dm_WorkSpace%XP(I)%Data, &
                                                dm_RefStrain(I)%DV,         dm_RefStrain(I)%SID)
             end do
           !$$--- copy the coordinate number out
           if(present(hDV)) then
             call COPYOUT_StrainField(hDV, hSID)
            end if

        RETURN
  end subroutine Cal_StrainField_DEV
  !****************************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_StrainField_template(IDEV, hDV, dDV, hSID, dSID)
  !***  PURPOSE:  the templet of copyout AC from device to host
  !
  !     INPUT:     IDEV,  the ID of device
  !                dDV,   the displacementg vector on the device IDEV
  !
  !     OUTPUT     hCN:   the displacementg vecto on host
  !
  use MD_Globle_Variables_GPU, only:m_DEVICES, hm_IA1th, hm_NAC, m_STARTCELL, m_ENDCELL

  implicit none
  !----   DUMMY Variables
         integer::IDEV
         real(KINDSF),  dimension(:,:)::hDV
         real(KINDSF),  device, dimension(:,:)::dDV
         integer,       dimension(:)::hSID
         integer,       device, dimension(:)::dSID

  !----   Local variables
         integer  IERR, CURDEV, STARTA, ENDA, NA

            IERR = cudaGetDevice(CURDEV)
            IERR = cudaSetDevice(m_DEVICES(IDEV))

            !$$--- the first atom on the device
            STARTA = hm_IA1th(m_STARTCELL(IDEV))

            !$$--- the last atom on the device
             ENDA = hm_IA1th(m_ENDCELL(IDEV))+hm_NAC(m_ENDCELL(IDEV))-1

            !$$--- the number of atoms on the device
            NA = ENDA - STARTA + 1

           IERR = cudaMemcpyAsync(hDV(STARTA,1), dDV(1,1), NA)
           IERR = cudaMemcpyAsync(hDV(STARTA,2), dDV(1,2), NA)
           IERR = cudaMemcpyAsync(hDV(STARTA,3), dDV(1,3), NA)
           IERR = cudaMemcpyAsync(hSID(STARTA),  dSID(1),  NA)
           IERR = cudaSetDevice(CURDEV)
         return
  end subroutine COPYOUT_StrainField_template
  !**********************************************************************************

  !**********************************************************************************
  subroutine COPYOUT_StrainField(hDV, hSID)
  !***  PURPOSE:  to copy the displacement vector of atoms from device  to
  !               hots array
  !
  !     INPUT:
  !     OUTPUT     hDV, array on host the displacement vector
  !

  use MD_Globle_Variables_GPU
  implicit none
  !----   DUMMY Variables
          real(KINDSF), dimension(:,:)::hDV
          integer,      dimension(:)  ::hSID

  !----   Local variables
         integer::I
  !***

         !$$--- copy the Occupancy out
         do I=1, m_NDEVICE
            call COPYOUT_StrainField_template(I, hDV, dm_RefStrain(I)%DV, hSID, dm_RefStrain(I)%SID)
         end do


         return
  end subroutine COPYOUT_StrainField
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_StrainField(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to calcualte and putout the occupation state of lattice and interstials
  !                  identified using method#1
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  use MD_Globle_Variables_GPU, only:SynchronizeDevices, dm_NPRT,  hm_ITYP, hm_XP, hm_EPOT, m_XP
  use MD_ForceLib_Factory_GPU, only:CalEpot_ForceClass, gm_ForceClass
  implicit none
       !--- dummy variables
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables
       character*256::GFILE
       integer::hFile, I, J, IP, SID, ICFG, NPRT0, ERR, BSHIFT,JOB,ITIME
       real(KINDDF)::LATT, TIME, DE
       real(KINDSF), dimension(:,:),allocatable::hDV
       integer,      dimension(:),allocatable  ::hSID
       !----
           NPRT0 = SimBox(1)%NPRT
           LATT  = SimBox(1)%RR
           allocate(hDV(dm_NPRT, 3), hSID(dm_NPRT), STAT=ERR)
           hDV  = 0.D0
           hSID = 0
           call Cal_StrainField_DEV(SimBox(1), CtrlParam, hDV, hSID)
           !call CopyXPFrom_Devices_to_DEV0()
           call SynchronizeDevices()

           JOB    = Stamp%ITest
           ICFG   = Stamp%IRec(1)
           ITIME  = Stamp%ITime
           TIME   = Stamp%Time
           BSHIFT = (JOB-1)*size(SimBox)
           do I=1, size(SimBox)
              BSHIFT = BSHIFT + 1
              GFILE = ""
              call STRCATI(GFILE, m_OUTFILE, "P", m_processid, 4)
              call STRCATI(GFILE, GFILE, "_", BSHIFT, 4)
              call STRCATI(GFILE, GFILE, ".", ICFG, 4)
              call AvailableIOUnit(hFile)
              write(*,*) "Save reference displacements to "//GFILE(1:len_trim(GFILE))

              !$$--- output occupation distribution
              open(UNIT=hFile, file = GFILE, status='unknown')
              write(hFile, fmt="(A)") "!--- REFDISPLACEMENT RESULTS CREATED BY "//gm_ExeName(1:len_trim(gm_ExeName))
              write(hFile, fmt="(A)") '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
              write(hFile, fmt="(A)") '!    AUTHOR: HOU Qing'
              write(hFile, fmt="(A)") '!    '

              write(hFile, fmt="('!--- The number of reference lattices: ', I8)") hm_NPRT0
              write(hFile, fmt="('!--- The number of atoms             : ', I8)") NPRT0
              write(hFile, fmt="('!--- at time steps: ', I7,  ', time(ps): ', 1pE12.4)") ITIME, TIME

              !$$--- write out the XYZ format
              write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
              write(hFile, fmt="(A,1X,I8)")            "&NATOM      ", NPRT0
              write(hFile, fmt="(A,1X,3(1PE16.7,1X))") "&BOXSIZE    ", hm_BOXSIZE0/SimBox(1)%RR
              write(hFile, fmt="(A,1X,3(1PE16.7,1X))") "&BOXLOW     ", hm_BOXLOW0/SimBox(1)%RR
              write(hFile, fmt="(A,1X,3(1PE16.7,1X))") "&LATT       ", SimBox(1)%RR*CP_CM2A
              write(hFile, fmt="(A,1X,3(I4,1X),A)")    "&TYPECOL    ", 1, 1, "'I'"
              write(hFile, fmt="(A,1X,3(I4,1X),A)")    "&XYZCOL     ", 2, 3, "'D'"
              write(hFile, fmt="(A,1X,3(I4,1X),A)")    "&DVCOL      ", 5, 3, "'D'"
              write(hFile, fmt="(A,1X,3(I4,1X),A)")    "&SITECOL    ", 8, 1, "'I'"

              write(hFile,FMT="(20A))")"!TYPE        ", "POSITION(x)     ", "     (y)        ", "      (z)       ",    &
                                                        "DISPLACE(dx)    ", "     (dy)       ", "      (dz)       ",   &
                                                        "NEAREST-SITE    "
              IP = (I-1)*NPRT0
              do J=1, NPRT0
                 IP  = IP + 1
                 write(hFile,fmt="(I8,2X,6(1PE16.7,1X),I6)") hm_ITYP(IP), hm_XP(IP,1:3)/LATT, -hDV(IP,1:3)/LATT, hSID(IP)
             end do
          end do
          deallocate(hDV, hSID, STAT=ERR)

          return
  end subroutine RECORD_StrainField
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_StrainField_TOOL(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to provide a interface to MD_SimBoxArray_ToolShell_14_GPU
  !                  It is assumed the the neighbor list routine has be called
  !
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  implicit none
       type(MDRecordStamp) ,         intent(in)::Stamp
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables

          if(Stamp%ITime .LT. 0) return
          call RECORD_StrainField(Stamp, SimBox, CtrlParam)

          return
  end subroutine RECORD_StrainField_TOOL
  !****************************************************************************************

  end module RefStrainField_GPU

