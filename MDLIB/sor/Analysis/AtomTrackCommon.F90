  module AtomTrackCommon
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to common framework of tracking atoms satifying given condition in
  !                  the starting configuration
  !
  !                  DEPENDENCE____________________________________________________________________________
  !
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2015-05 (Hou Qing)
  !                  2018-04 (HOU Qing): change the routine of loading control parameters,
  !                                      because of the introducing of DatPad in SimMDBox
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl

  !-----------------------------------------------
    implicit none

      integer::m_processid=0
      !$$--- the id of  filename get from SimMDCtrl for I/O.
      !$$    Refering the definition of f_others in SimMDCtrl.
      character(len=12),parameter, private::mp_FTAGI="&AUXF_SELIN"
      character(len=13),parameter, private::mp_FTAGO="&AUXF_SELOUT"
      character(len=256)::m_INFILE =""                       ! filename of input control data
      character(len=256)::m_OUTFILE =""                      ! filename of output data
      integer::m_SumUNIT                                     ! I/O unit for summary information

      !--- the interface for Flag procedure.
      private::SELECTPROC
      abstract interface
        SUBROUTINE SELECTPROC(SimBox, CtrlParam, SEL)
        !***  PURPOSE:   to create the FLAG witch identify the boxes
        !                that statisfying criteria
        !     INPUT:     SimBox:
        !                CtrlParam:
        !     OUTPUT:    SEL,   marker of the boxes have been selected
        !
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
          integer,        dimension(:)::SEL

       END SUBROUTINE SELECTPROC
      end interface
      procedure(SELECTPROC), pointer, private::m_pTrackProc=>null()
      procedure(SELECTPROC), pointer, private::m_pSelOuputProc=>null()

      integer, dimension(:), allocatable::m_NUMINIATOMS    ! the number of atoms that are identified in boxes
      integer, dimension(:), allocatable::m_INIATOMS       ! the ID of atoms to be identified

      integer::m_ATYP(mp_MXGROUP) = 1                      ! the type atoms to be traking
      real(KINDDF),dimension(:), allocatable::m_EPOTTHD    ! the potential threshold
      real(KINDDF),dimension(:), allocatable::m_EKINTHD    ! the kinetic energy threshold

      !--- calculation control parameters
      integer, private::m_INITED = 0

  contains
  !*****************************************************************************
  subroutine SetTrackProc(selectPROC)
  !***  DESCRIPTION: to set the user defrine FlagProc, which will be used to
  !                  create the flags of atoms that will be considered in traking.
  !
  !     INPUT: selectPROC,  the subroutine provided by user
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
        SUBROUTINE SELECTPROC(SimBox, CtrlParam, SEL)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
          integer,        dimension(:)::SEL
       END SUBROUTINE SELECTPROC
  end interface
  !--- END INTERFACE --------------------------------
               m_pTrackProc  =>selectPROC
           return
   end subroutine SetTrackProc
  !*********************************************************************************

  !*****************************************************************************
  subroutine SetSelOuputProc(selectPROC)
  !***  DESCRIPTION: to set the user defrine process, which will be used to
  !                  create the flags of atoms that will be output but not
  !                  identified for tracking
  !
  !     INPUT: selectPROC,  the subroutine provided by user
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
        SUBROUTINE SELECTPROC(SimBox, CtrlParam, SEL)
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
          integer,        dimension(:)::SEL
       END SUBROUTINE SELECTPROC
  end interface
  !--- END INTERFACE --------------------------------
               m_pSelOuputProc  =>selectPROC
           return
   end subroutine SetSelOuputProc
  !*********************************************************************************

  !*****************************************************************************
  subroutine Clear_WorkingArray()
  !***  DESCRIPTION: to allocate working memory
  !
  implicit none

          if(allocated(m_NUMINIATOMS)) deallocate(m_NUMINIATOMS)
          if(allocated(m_INIATOMS))    deallocate(m_INIATOMS)

          if(allocated(m_EPOTTHD))     deallocate(m_EPOTTHD)
          if(allocated(m_EKINTHD))     deallocate(m_EKINTHD)

          m_INITED = 0
          return
  end subroutine Clear_WorkingArray
  !*****************************************************************************

  !**********************************************************************************
  subroutine Clear_AtomTrackCommon(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox) ::SimBox
       type(SimMDCtrl)::CtrlParam

          call Clear_WorkingArray()
          return
  end subroutine Clear_AtomTrackCommon
  !**********************************************************************************

  !**********************************************************************************
   subroutine Initialize_AtomTrackCommon(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, SelCFG0, SelCFG1, IFILE

         !$$--- to clear the memory allocated before if there is
         if(m_INITED .gt. 0) then
           call Clear_WorkingArray()
         end if

        !$$--- to findout the I/O unit
         IFILE = 0
         do I=1, SIZE(CtrlParam%f_tag)
            if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGI)) .EQ. mp_FTAGI) then
               IFILE = I
            end if
         end do

         m_OUTFILE = ""
         do I=1, SIZE(CtrlParam%f_tag)
            if(CtrlParam%f_tag(I)(1:len_trim(mp_FTAGO)) .EQ. mp_FTAGO) then
                m_OUTFILE = CtrlParam%f_others(I)
            end if
         end do

         !$$--- check the input files
         if(IFILE .gt. 0) then
            m_INFILE = CtrlParam%f_others(IFILE)
            write(*,fmt="(A)") " !**** Loading control data from: "//m_INFILE(1:len_trim(m_INFILE))
            if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                         m_INFILE(1:len_trim(m_INFILE))//" for box selector"
            call LoadControlParameters(m_INFILE, SimBox, CtrlParam)
         end if

         if(len_trim(m_OUTFILE) .LE.0 ) then
             write(*,fmt="(A)")          " MDPSCU Error: no output file for "//gm_ExeName(1:len_trim(gm_ExeName))// " is given."
             write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
             write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
             write(*,fmt="(' Process to be stopped')")
             stop
         else
               call CreateDataFolder_Globle_Variables(m_OUTFILE)
         end if

         !$$--- open the summary file
         call AvailableIOUnit(m_SumUNIT)
         open(UNIT=m_SumUNIT, file = m_OUTFILE(1:len_trim(m_OUTFILE))//".lst", status='unknown')
         write(m_SumUNIT,fmt="(A)") "!--- Namelist in the box selection"
         write(m_SumUNIT,fmt="(A)") "!---- Original box -> box"

         return
   end subroutine Initialize_AtomTrackCommon
  !**********************************************************************************

  !*********************************************************************************
  subroutine LoadControlParameters(fname, SimBox, CtrlParam)
  !***  PURPOSE:   to readin control parameters from a file
  !     INPUT:     fname: the file name
  !
  !     OUTPUT:    SimBox,   the simulation box, with the member propColTitle
  !                          could be changed.
  !                CtrlParam, the control parameters, with the member NEEDDAMP
  !                           could be changed.
  !
  !
  !
  use MD_Globle_Variables,only:Load_ExAnalyCtl_SimMDCtrl
  implicit none
  !----   DUMMY Variables
   character*(*)   ::fname
   type(SimMDBox)  ::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
   integer::I, N, IS, LINE, NN
   character*256::STR
   character*32::STRNUMB(10), KEYWORD
   type(StatementList), pointer::StatList

  !-------------

             !$$--- load input statements
             call Load_ExAnalyCtl_SimMDCtrl(mp_FTAGI, Fname, SimBox, CtrlParam, &
                                                OutTag=mp_FTAGO, OutFile=m_OUTFILE, TheInput=StatList)

             do IS=1, Number_StatementList(StatList)
                call Get_StatementList(IS, StatList, STR, Line=LINE)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)

                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         !*** to get the properties to be included in clustering
                         case("&PROP_TYPE")
                              call Extract_Numb(STR,SimBox%nGroup,N,STRNUMB)
                              m_ATYP = 0
                              do I=1, N
                                 NN = ISTR(STRNUMB(I))
                                 if(NN .gt.0) m_ATYP(NN) = 1
                              end do
                              if(count(m_ATYP .gt.0) .le. 0) then
                                write(*,fmt="(A)")          " MDPSCU Error: atom type is required for tracking calculation"
                                write(*,fmt="(A,A,A)")      "               but no atom type is available "
                                write(*,fmt="(A, BZI6)")    '               check control file at line:', LINE
                                write(*,fmt="(A, BZI6)")    '        Usage: &PROP_TYPE ty1, ty2...'
                                 write(*,fmt="(A)")         '               Process to be stopped'
                                stop
                              end if

                         case("&PROP_EPOT")
                              call Extract_Numb(STR,1,N,STRNUMB)
                              if(N.ge.1) then
                                 NN = 2*ISTR(STRNUMB(1))
                                 if(.not.allocated(m_EPOTTHD))allocate(m_EPOTTHD(NN))
                                 call Extract_Numb(STR,NN+1,N,STRNUMB)
                                 if(N.LE.NN+1) then
                                    write(*,fmt="(A,I4,A)")  ' MDPSCU Error: there are ',NN/2, ' pairs of EPOT range are required'
                                    write(*,fmt="(A,I4,A)")  '               but only ', (N-1)/2, ' pairs are available'
                                    write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                    write(*,fmt="(A, BZI6)") '        Usage: &PROP_EPOT n, ( E1, E2 ), ( E3, E4 ),...,(En-1, En)'
                                    write(*,fmt="(A)")       '               Process to be stopped'
                                    stop
                                 end if
                                 do I=1, NN
                                    m_EPOTTHD(I*2-1) = DRSTR(STRNUMB(I*2))
                                    m_EPOTTHD(I*2)   = DRSTR(STRNUMB(I*2+1))
                                 end do
                              end if

                         case("&PROP_EKIN")
                              call Extract_Numb(STR,1,N,STRNUMB)
                              if(N.ge.1) then
                                 NN = 2*ISTR(STRNUMB(1))
                                 if(.not.allocated(m_EKINTHD))allocate(m_EKINTHD(NN))
                                 call Extract_Numb(STR,NN+1,N,STRNUMB)
                                 if(N.LE.NN+1) then
                                    write(*,fmt="(A,I4,A)")  ' MDPSCU Error: there are ',NN/2, ' pairs of EKIN range are required'
                                    write(*,fmt="(A,I4,A)")  '               but only ', (N-1)/2, ' pairs are available'
                                    write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                    write(*,fmt="(A, BZI6)") '        Usage: &PROP_EKIN n, ( E1, E2 ), ( E3, E4 ),...,(En-1, En)'
                                    write(*,fmt="(A)")       '               Process to be stopped'
                                    stop
                                 end if
                                 do I=1, NN
                                    m_EKINTHD(I*2-1) = DRSTR(STRNUMB(I*2))
                                    m_EKINTHD(I*2)   = DRSTR(STRNUMB(I*2+1))
                                 end do
                              end if

                  end select
              end do
              call Release_StatementList(StatList)
            return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !*********************************************************************************
  subroutine Identify_Atoms(SimBox, CtrlParam)
  !***  DESCRIPTION: to identify indentify the atoms satisfying the type, potential
  !                  energy crititier
  !
  !    INPUT:
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  implicit none
       type(SimMDBox), dimension(:), intent(in)::SimBox
       type(SimMDCtrl),              intent(in)::CtrlParam
       !--- local variables
        integer, dimension(:), allocatable::FLAG
        integer::I, J, IP,IPA, NP0, NPRT, TNUM

          !$$*** create the ID of atoms that to be output
           NP0  = SimBox(1)%NPRT
           NPRT = size(SimBox)*NP0
           allocate(FLAG(NPRT))
           FLAG = 0
           if(associated(m_pTrackProc)) then
              call m_pTrackProc(SimBox, CtrlParam, FLAG)
           else
              call Default_AtomIdentify(SimBox, CtrlParam, FLAG)
           end if

           if(allocated(m_NUMINIATOMS) ) deallocate(m_NUMINIATOMS)
           allocate(m_NUMINIATOMS(size(SimBox)))
           TNUM = 0
           IP = 1
           do I=1, size(SimBox)
              m_NUMINIATOMS(I) = count(FLAG(IP:IP+NP0-1).gt.0)
              TNUM = TNUM + m_NUMINIATOMS(I)
              IP = IP + NP0
           end do

           if(allocated(m_INIATOMS)) deallocate(m_INIATOMS)
           allocate(m_INIATOMS(TNUM ) )
           m_INIATOMS = 0
           IP  = 0
           IPA = 0
           do I=1, size(SimBox)
              do J=1, SimBox(I)%NPRT
                 IP = IP + 1
                 if(FLAG(IP) .gt. 0) then
                    IPA = IPA + 1
                    m_INIATOMS(IPA) = J
                 end if
              end do
           end do

          deallocate(FLAG)
          return
  end subroutine Identify_Atoms
  !*********************************************************************************

  !*********************************************************************************
  subroutine Default_AtomIdentify(SimBox, CtrlParam, FLAG)
  !***  DESCRIPTION: to indentify the atoms by their type, potential
  !                  and energy crititier
  !
  !    INPUT:
  !            SimBox,    the simulation boxs
  !            CtrlParam, the control parameters for simulation
  !
  !   OUTPUT:  FLAG,      the ID of atoms satisfying the given criteria
   implicit none
       type(SimMDBox), dimension(:), intent(in) ::SimBox
       type(SimMDCtrl),              intent(in) ::CtrlParam
       integer,        dimension(:), intent(out)::FLAG
       !--- local variables
        integer::I, J, K, IP, NPP, NPK, FLAGT, FLAGP, FLAGK

          !$$*** create the ID of atoms that to be output
           if(allocated(m_EPOTTHD)) then
              NPP = size(m_EPOTTHD)/2
           else
              NPP = 0
           end if

           if(allocated(m_EKINTHD)) then
              NPK = size(m_EKINTHD)/2
           else
              NPK = 0
           end if

             IP  = 0
             do I=1, size(SimBox)
                do J=1, SimBox(I)%NPRT
                   IP = IP + 1

                   !--- select atom by type
                   if(m_ATYP(SimBox(I)%ITYP(J)) .le. 0) cycle
                   FLAGT = 1

                   !--- select atom by potential
                   FLAGP = 0
                   do K=1, NPP
                      if(SimBox(I)%EPOT(J) .le. m_EPOTTHD(2*K-1) .and. SimBox(I)%EPOT(J) .le. m_EPOTTHD(2*K)) then
                         FLAGP = 1
                         exit
                      end if
                   end do
                   if(NPP.gt.0 .and. FLAGP .le. 0) cycle

                   !--- select atom by kinetic energy
                   FLAGK = 0
                   do K=1, NPK
                      if(SimBox(I)%EKIN(J) .le. m_EKINTHD(2*K-1) .and. SimBox(I)%EKIN(J) .le. m_EKINTHD(2*K)) then
                         FLAGK = 1
                         exit
                      end if
                   end do
                   if(NPK.gt.0 .and. FLAGK .le. 0) cycle

                   FLAG(IP) = 1
                end do
             end do
          return
  end subroutine Default_AtomIdentify
  !*********************************************************************************

  !*********************************************************************************
  subroutine RECORD_AtomTrackCommon(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the configuration
  !
  !    INPUT:  SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !            Stamp,       the record stamp
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
        character*256::GFILE0, GFILE1
        integer::I, J, IPA, IPA0, IP, BSHIFT, hFile, NPRT
        integer, dimension(:), allocatable::IDA
  !---  start process

          !$$--- to identify the atoms
          if(Stamp%ITime .lt. 0) return

          if(.not.allocated(m_NUMINIATOMS)) then
             call Identify_Atoms(SimBox, CtrlParam)
             !$$--- output the headr of the tracking file for recording the atoms in a singles file,
             !      not used at the present time
             !      call InitialRecord(JOB, ITIME, TIME, SimBox, CtrlParam)
          end if

          !$$*** create the ID of atoms that to be output
          NPRT = size(SimBox)*SimBox(1)%NPRT
          allocate(IDA(NPRT))

             IDA  = 0
             IP   = 0
             IPA0 = 0
             do I=1, size(SimBox)
                IPA = IPA0
                do J=1, m_NUMINIATOMS(I)
                   IP = IP + 1
                   IPA = IPA + 1
                   IDA(IPA) =  m_INIATOMS(IP)
                end do
                IPA0 = IPA0 + SimBox(I)%NPRT
             end do

             call STRCATI(GFILE0, CtrlParam%f_geometry, "P", m_processid, 4)
             call STRCATI(GFILE0, GFILE0, "_", Stamp%ITest, 4)
             call STRCATI(GFILE0, GFILE0, ".", Stamp%IRec(1), 4)
             write(m_SumUNIT,fmt="(A)") "Original box name: "//GFILE0(1:len_trim(GFILE0))

             BSHIFT = (Stamp%ITest-1)*size(SimBox)
             do I=1, size(SimBox)
                IPA    = (I-1)*SimBox(1)%NPRT
                BSHIFT = BSHIFT + 1
                call STRCATI(GFILE1, m_OUTFILE, "P", m_processid, 4)
                call STRCATI(GFILE1, GFILE1, "_", BSHIFT, 4)
                call Putout_Instance_Config_SimMDBox(GFILE1, SimBox(I), Stamp, IDA=IDA(IPA+1:IPA+SimBox(I)%NPRT) )
                call STRCATI(GFILE1, GFILE1, ".", Stamp%IRec(1), 4)
                write(m_SumUNIT,fmt="(A, I6, A)") " box",I, "-> "//GFILE1(1:len_trim(GFILE1))
             end do

          deallocate(IDA)
          return

  end subroutine RECORD_AtomTrackCommon
  !*********************************************************************************

  !*********************************************************************************
  subroutine InitialRecord(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to start the record
  !
  !    INPUT:  Stamp,       the recoding stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  implicit none
       type(MDRecordStamp),          intent(in) :: Stamp
       type(SimMDBox), dimension(:), intent(in) :: SimBox
       type(SimMDCtrl),              intent(in) :: CtrlParam
       !--- local variables
        character*256::GFILE1
        integer::I, J, IPA, IP,  BSHIFT, hFile, JOB
       !----
            !$$--- to prepair the filenames
             JOB    = Stamp%ITest
             BSHIFT = (JOB-1)*size(SimBox)
             IPA = 0
             do I=1, size(SimBox)
                BSHIFT = BSHIFT + 1
                call STRCATI(GFILE1, m_OUTFILE, "P", m_processid, 4)
                call STRCATI(GFILE1, GFILE1, "_", BSHIFT, 4)
                call AvailableIOUnit(hFile)
                open(UNIT=hFile, file = GFILE1, status='unknown')

                  write(hFile, fmt="(A, I8)")           '!--- THE TRACKING RESULTS CREATED BY '//gm_ExeName(1:len_trim(gm_ExeName))
                  write(hFile, fmt="(A)")               '!    A PROGRAM GENERATED BY THE MD PACKAGE AT SICHUAN UNIVERSITY'
                  write(hFile, fmt="(A)")               '!    AUTHOR: HOU Qing'
                  write(hFile, fmt="(A)")               '!    '
                  write(hFile, fmt="(A, I8)")           '!--- Original JOB ID: ', JOB
                  write(hFile, fmt="(A, I8)")           '!--- Original BOX ID: ', BSHIFT
                  write(hFile, fmt="(A, I8)")           '!--- From configure ID:           ', CtrlParam%STARTCFG
                  write(hFile, fmt="(A, I8)")           '!--- To configure ID:            ',  CtrlParam%ENDCFG
                  write(hFile, fmt="(A, I8)")           '!--- With interval               ',  CtrlParam%CFGSTEP
                  write(hFile, fmt="(A, I8)")           '!--- Number of atoms in tracking: ', m_NUMINIATOMS(I)

                  !$$--- write out the XYZ format
                  write(hFile, fmt="(A)")                  PKW_OUTCFG_FORMATXYZ
                  write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXSIZE  ", SimBox(I)%ZL/SimBox(1)%RR
                  write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&BOXLOW   ", SimBox(I)%BOXLOW/SimBox(1)%RR
                  write(hFile, fmt="(A,1X,3(1PE13.4,1X))") "&LATT     ", SimBox(I)%RR*CP_CM2A

                  write(hFile, fmt="(A,1X,I8)")            "&NATOMI   ", m_NUMINIATOMS(I)
                  write(hFile, fmt="(A,1X,I8)")            "&NATOM    ", m_NUMINIATOMS(I)*IABS(CtrlParam%ENDCFG - CtrlParam%STARTCFG+1)/CtrlParam%CFGSTEP

                  write(hFile, fmt="(A,1X,3(I4,1X))")      "&XYZCOL   ", 1, 2, 3
                  write(hFile, fmt="(A,1X,3(I4,1X))")      "&TYPECOL  ", 4
                  write(hFile, fmt="(A,1X,I4, A)")         "&ID0COL   ", 5, " !--- ID of the indentified atoms in the original box"
                  write(hFile, fmt="(A,1X,I4, A)")         "&TYPE1COL ", 6, " !--- new type = type + group number of the indentified atoms"
                  write(hFile, fmt="(12A)")                "!---",     "   X(latt.)  ", &
                                                                       "   Y(latt.)  ", &
                                                                       "   Z(latt.)  ", &
                                                                       "     TYPE    ", &
                                                                       "     ID0     ", &
                                                                       "    TYPE1    "
                close(hFile)
             end do
          return
  end subroutine InitialRecord
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_AtomTrackCommonTool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !
  !    INPUT:  Stamp,       the record stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_ToolShell_16_GPU.F90
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables


            call RECORD_AtomTrackCommon(Stamp, SimBox, CtrlParam)
            return
  end subroutine RECORD_AtomTrackCommonTool
  !****************************************************************************************

  !****************************************************************************************
  end module AtomTrackCommon
