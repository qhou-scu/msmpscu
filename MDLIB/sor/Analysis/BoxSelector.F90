  module BoxSelector
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to common framework of select boxs according to criteria supplied by user
  !
  !                  DEPENDENCE____________________________________________________________________________
  !
  !                  SEE ALSO____________________________________________________________________________
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2014-11 (Hou Qing)
  !
  !________________________________________________________________________________________________________


  !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
    use MD_TYPEDEF_SimMDBox
    use MD_TYPEDEF_SimMDCtrl
    use MD_TYPEDEF_Recordstamp

  !-----------------------------------------------
    implicit none

      integer::m_processid=0
      !$$--- the id of  filename get from SimMDCtrl for I/O.
      !$$    Refering the definition of f_others in SimMDCtrl.
      character(len=12),parameter, private::mp_FTAGI="&AUXF_SELIN"
      character(len=13),parameter, private::mp_FTAGO="&AUXF_SELOUT"
      character(len=256)::m_INFILE =""                       ! filename of input control data
      character(len=256)::m_OUTFILE =""                      ! filename of output data
      character*256, private::m_RefCfgFile= ""               ! filename of a reference configuration
                                                             ! if the reference configuration is given
                                                             ! the output boxes will have reference box inserted
      type(SimMDBox)::hm_RefSimBox                      !-- the reference box


      integer::m_SumUNIT                                     ! I/O unit for summary information
      integer, dimension(:), allocatable::m_SelBOXNum
      !--- the interface for Flag procedure.
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
          integer,         dimension(:)::SEL

       END SUBROUTINE SELECTPROC
      end interface
      procedure(SELECTPROC), pointer, private::m_pSelectProc=>null()

      !--- calculation control parameters
      integer, private::m_INITED = 0


  contains

  !*****************************************************************************
  subroutine SetSelectProc(selectPROC)
  !***  DESCRIPTION: to set the user defrine FlagProc, which will be used to
  !                  create the flags of atoms that will be considered in clustering.
  !
  !     INPUT: PRERECORD,  the subroutine provided by user for pre-recording
  !
  implicit none
  !--- interface to the external routine -------------------
   interface
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
  !--- END INTERFACE --------------------------------
               m_pSelectProc  =>selectPROC
           return
   end subroutine SetSelectProc
  !*********************************************************************************

  !*****************************************************************************
  subroutine DoSelectProc(SimBox, CtrlParam, SEL)
  !***  DESCRIPTION: to set the user defrine FlagProc, which will be used to
  !                  create the flags of atoms that will be considered in clustering.
  !
  !     INPUT: PRERECORD,  the subroutine provided by user for pre-recording
  !
  implicit none
      type(SimMDBox), dimension(:)::SimBox
      type(SimMDCtrl)             ::CtrlParam
      integer,        dimension(:)::SEL

          if( .not. associated(m_pSelectProc)) then
              write(*,fmt="(' MDPSCU Error: box selecting routine is not provided in BoxSelector module')")
              write(*,fmt="('               Process to be stopped')")
              stop
          end if
          call m_pSelectProc(SimBox, CtrlParam, SEL)
          return
   end subroutine DoSelectProc
  !*********************************************************************************

  !*****************************************************************************
  subroutine Clear_WorkingArray()
  !***  DESCRIPTION: to allocate working memory
  !
  implicit none
          if(allocated(m_SelBOXNum)) deallocate(m_SelBOXNum)
          m_INITED = 0
          return
  end subroutine Clear_WorkingArray
  !*****************************************************************************

  !**********************************************************************************
  subroutine Clear_BoxSelector(SimBox, CtrlParam)
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
  end subroutine Clear_BoxSelector
  !**********************************************************************************

  !**********************************************************************************
   subroutine Initialize_BoxSelector(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, SelCFG0, SelCFG1, IFILE
   type(SimMDBox)::tSimBox

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
             write(*,fmt="(A,A,A)")      "               add the keyword in your SETUP  file: ", mp_FTAGO, " fname,"
             write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in your ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
             write(*,fmt="(' Process to be stopped')")
             stop
         else
               call CreateDataFolder_Globle_Variables(m_OUTFILE)
         end if

         !$$--- determine the slection range
         SelCFG0 = CtrlParam%STARTCFG
         if(CtrlParam%ENDCFG .gt. 0) then
            SelCFG1 = CtrlParam%ENDCFG
         else
            !*** give a large number
             call NumberCfg_SimMDCtrl(CtrlParam, SelCFG1)
         end if
         allocate(m_SelBOXNum(min(SelCFG1, SelCFG0):max(SelCFG1, SelCFG0)) )
         m_SelBOXNum = 0

         !$$--- if reference box is given, create the reference box
          if(len_trim(m_RefCfgFile) .gt. 0) then
             write(*,*) "!************ BoxSelector module to be performed **********"
             write(*,*) "!    With the reference box: "//m_RefCfgFile(1:len_trim(m_RefCfgFile))//" to be inserted in slected boxes"
             write(*,*) "!    "

            !$$--- get box size so on
            call CopyInformation_SimMDBox(SimBox, hm_RefSimBox)
            hm_RefSimBox%NPRT = SimBox%NPRT
            call Initialize_SimMDBox(hm_RefSimBox,2)

            call ClearDataProKWD_SimMDBox(hm_RefSimBox)
            call AddDataProKWD_SimMDBox(hm_RefSimBox,"XYZCOL")
            call Read_Initial_Config_SimMDBox(hm_RefSimBox, fname=m_RefCfgFile, fmt=0, mod=1)

            hm_RefSimBox%NGROUP = SimBox%NGROUP+1
            hm_RefSimBox%ITYP   = hm_RefSimBox%NGROUP
            hm_RefSimBox%STATU = CP_STATU_ACTIVE

            hm_RefSimBox%NA  = 0
            hm_RefSimBox%NA(hm_RefSimBox%NGROUP)    = hm_RefSimBox%NPRT
            hm_RefSimBox%IPA(hm_RefSimBox%NGROUP)   = 1
            hm_RefSimBox%IPA(hm_RefSimBox%NGROUP+1) = hm_RefSimBox%NPRT+1
          end if

         !$$--- open the summary file
         call AvailableIOUnit(m_SumUNIT)
         open(UNIT=m_SumUNIT, file = m_OUTFILE(1:len_trim(m_OUTFILE))//".lst", status='unknown')
         write(m_SumUNIT,fmt="(A)") "!--- Namelist in the box selection"
         if(len_trim(m_RefCfgFile) .gt. 0) then
             write(m_SumUNIT,fmt="(A)") "!---- Original box -> box, "//"with reference box "//m_RefCfgFile(1:len_trim(m_RefCfgFile))//"inserted"
         else
             write(m_SumUNIT,fmt="(A)") "!---- Original box -> box"
         end if


         return
   end subroutine Initialize_BoxSelector
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
   character*(*)  ::fname
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
   integer::N, IS, LINE, I
   character*256::STR
   character*32::KEYWORD
   type(StatementList), pointer::StatList


            call Load_ExAnalyCtl_SimMDCtrl(mp_FTAGI, Fname, SimBox, CtrlParam, &
                                               OutTag=mp_FTAGO, OutFile=m_OUTFILE, TheInput=StatList)

            do IS=1, Number_StatementList(StatList)
               call Get_StatementList(IS, StatList, STR, Line=LINE)
               call GetKeyWord("&", STR, KEYWORD)
               call UpCase(KEYWORD)
                  select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                         case( "&REFCFG")
                                m_RefCfgFile = ""
                                call Extract_Substr(STR,1,n,m_RefCfgFile)

                  end select
            end do
            call Release_StatementList(StatList)
            return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !****************************************************************************************
  subroutine RECORD_BoxSelector(Stamp, SimBox, CtrlParam )
  !***  DESCRIPTION:
  !
  !    INPUT:  SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !            Stamp,       the record stamp
  !
  !   SEE ALSO:
  !            MD_SimBoxArray_AppShell_16_GPU.F90
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
        integer, dimension(:), allocatable::FLAG
        character*256::GFILE0, GFILE1
        integer::I
        type(SimMDBox)::TSIMBOX
        type(MDRecordStamp)::TSTAMP

          call Copy_RecordStamp(Stamp, TSTAMP)
          allocate(FLAG(size(SimBox)))
             FLAG = 0
             call DoSelectProc(SimBox, CtrlParam, FLAG)
             !call GetSectID_SimMDCtrl(CtrlParam, ITIME, ISECT)
             !call GetCfgID_SimMDCtrl(CtrlParam, ITIME, ICFG)

             call STRCATI(GFILE0, CtrlParam%f_geometry, "P", m_processid, 4)
             call STRCATI(GFILE0, GFILE0, "_", Stamp%ITest, 4)
             call STRCATI(GFILE0, GFILE0, ".", Stamp%IRec(1), 4)
             write(m_SumUNIT,fmt="(A)") "Original box name: "//GFILE0(1:len_trim(GFILE0))
             if(all(FLAG.le.0)) then
                write(m_SumUNIT,fmt="(A, I6, A)") " no box found satisfying selection condition"
             else
                TSTAMP%ITest = 1
                do I=1, size(FLAG)
                   if(FLAG(I) .gt. 0) then
                     m_SelBOXNum(Stamp%IRec(1)) =  m_SelBOXNum(Stamp%IRec(1)) + 1
                     call Copy_SimMDBox(hm_RefSimBox, tSimBox)
                     call Merge_SimMDBox(SimBox(I), tSimBox)

                     call STRCATI(GFILE1, m_OUTFILE, "P", m_processid, 4)
                     call STRCATI(GFILE1, GFILE1, "_", m_SelBOXNum(Stamp%IRec(1)), 4)

                     TSTAMP%IBox = m_SelBOXNum(Stamp%IRec(1))
                     call Putout_Instance_Config_SimMDBox(GFILE1, tSimBox, TSTAMP)
                     call STRCATI(GFILE1, GFILE1, ".", Stamp%IRec(1), 4)
                     write(m_SumUNIT,fmt="(A, I6, A)") " box",I, "-> "//GFILE1(1:len_trim(GFILE1))
                     write(*,fmt="(A, I6, A)") " box",I, "-> "//GFILE1(1:len_trim(GFILE1))
                   end if
                end do
             end if

          deallocate(FLAG)
          call Release_SimMDBox(tSimBox)
          return
  end subroutine RECORD_BoxSelector
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_BoxSelectorTool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !
  !    INPUT:  Stamp,       the record stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  !  SEE ALSO:
  !            MD_SimBoxArray_ToolShell.f90
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

       !--- local variables
            if(Stamp%ITime .lt. 0) return
            call RECORD_BoxSelector(Stamp, SimBox, CtrlParam)
            return
  end subroutine RECORD_BoxSelectorTool
  !****************************************************************************************

  !****************************************************************************************
  end module BoxSelector
