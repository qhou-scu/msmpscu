  module AtomSelectorCommon
  !**** DESCRIPTION: _______________________________________________________________________________________
  !                  The module is to common framework of select boxs according to criteria supplied by user
  !
  !                  DEPENDENCE____________________________________________________________________________
  !
  !                  SEE ALSO____________________________________________________________________________
  !                  ______________________________________________________________________________________
  !**** HISTORY:
  !                  version 1st 2015-05 (Hou Qing)
  !                  2018-04 (HOU Qing): change the routine of loading control parameters,
  !                                      because of the introducing of DatPad in SimMDBox
  !
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
          use MD_CONSTANTS
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)             ::CtrlParam
          integer,        dimension(:)::SEL

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
          use MD_CONSTANTS
          use MD_TYPEDEF_SimMDBox
          use MD_TYPEDEF_SimMDCtrl
          implicit none
          type(SimMDBox), dimension(:)::SimBox
          type(SimMDCtrl)            ::CtrlParam
          integer,        dimension(:)::SEL

       END SUBROUTINE SELECTPROC
   end interface
  !--- END INTERFACE --------------------------------
              m_pSelectProc =>selectPROC
           return
   end subroutine SetSelectProc
  !*********************************************************************************

  !*****************************************************************************
  subroutine Clear_WorkingArray()
  !***  DESCRIPTION: to allocate working memory
  !
  implicit none
          m_INITED = 0
          return
  end subroutine Clear_WorkingArray
  !*****************************************************************************

  !**********************************************************************************
  subroutine Clear_AtomSelectorCommon(SimBox, CtrlParam)
  !***  PURPOSE:  to deallocate the memories allocated
  !    INPUT:     SimBox,    the simulation box, useless here
  !               CtrlParam, the control parameter, useless here
  !   OUTPUT:
  !
   implicit none
       type(SimMDBox)  ::SimBox
       type(SimMDCtrl) ::CtrlParam

          call Clear_WorkingArray()
          return
  end subroutine Clear_AtomSelectorCommon
  !**********************************************************************************

  !**********************************************************************************
   subroutine Initialize_AtomSelectorCommon(SimBox, CtrlParam)
   !***  PURPOSE:  to initialize the calculation by allocate copy the number of atom in
   !     Voronoi volume on devices to a host array AC

   use MD_Globle_Variables, only:CreateDataFolder_Globle_Variables

   implicit none
  !----   DUMMY Variables
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam

  !----   Local variables
   integer::I, IFILE

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
   end subroutine Initialize_AtomSelectorCommon
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
  use MD_Globle_Variables, only:Load_ExAnalyCtl_SimMDCtrl
  implicit none
  !----   DUMMY Variables
   character*(*)  ::fname
   type(SimMDBox) ::SimBox
   type(SimMDCtrl)::CtrlParam
  !--- local
            !$$--- load input statements
            call Load_ExAnalyCtl_SimMDCtrl(mp_FTAGI, Fname, SimBox, CtrlParam, OutTag=mp_FTAGO, OutFile=m_OUTFILE)
            return
    end subroutine LoadControlParameters
  !*********************************************************************************

  !****************************************************************************************
  subroutine RECORD_AtomSelectorCommon(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the configure of selected atoms
  !
  !      INPUT:   Stamp,       the record stamp
  !               SimBox,      the simulation boxs
  !               CtrlParam,   the control parameters for simulation
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables
        integer, dimension(:), allocatable::FLAG
        integer, dimension(:), allocatable::IDA
        character*256::GFILE0, GFILE1
        integer::I, J, IPA0, IPA, IP, NPRT, BSHIFT

          !$$--- to check if the selector routine has been set
          if( .not. associated(m_pSelectProc)) then
              write(*,fmt="(' MDPSCU Error: atom selecting routine is not provided in AtomSelector module')")
              write(*,fmt="('               Process to be stopped')")
              stop
          end if
          NPRT = size(SimBox)*SimBox(1)%NPRT
          allocate(FLAG(NPRT), IDA(NPRT))

            !$$*** create the ID of atoms that to be output
             FLAG = 0
             IDA =  0
             call m_pSelectProc(SimBox, CtrlParam, FLAG)
             IP   = 0
             IPA0 = 0
             do I=1, size(SimBox)
                IPA = IPA0
                do J=1, SimBox(I)%NPRT
                   IP = IP + 1
                   if(FLAG(IP) .gt. 0) then
                      IPA = IPA + 1
                      IDA(IPA) =  J
                   end if
                end do
                IPA0 = IPA0 + SimBox(I)%NPRT
             end do

             !call GetSectID_SimMDCtrl(CtrlParam, ITIME, ISECT)
             !call GetCfgID_SimMDCtrl(CtrlParam, ITIME, ICFG)

             call STRCATI(GFILE0, CtrlParam%f_geometry, "P", m_processid, 4)
             call STRCATI(GFILE0, GFILE0, "_", Stamp%ITest,   4)
             call STRCATI(GFILE0, GFILE0, ".", Stamp%IRec(1), 4)
             write(m_SumUNIT,fmt="(A)") "Original box name: "//GFILE0(1:len_trim(GFILE0))

             BSHIFT = (Stamp%ITest-1)*size(SimBox)
             do I=1, size(SimBox)
                IPA    = (I-1)*SimBox(1)%NPRT
                BSHIFT = BSHIFT + 1
                call STRCATI(GFILE1, m_OUTFILE, "P", m_processid, 4)
                call STRCATI(GFILE1, GFILE1, "_", BSHIFT, 4)
                call Putout_Instance_Config_SimMDBox(GFILE1, SimBox(I),Stamp, IDA=IDA(IPA+1:IPA+SimBox(I)%NPRT) )
                call STRCATI(GFILE1, GFILE1, ".", Stamp%IRec(1), 4)
                write(m_SumUNIT,fmt="(A, I6, A)") " box",I, "-> "//GFILE1(1:len_trim(GFILE1))
             end do

          deallocate(FLAG, IDA)
          return
  end subroutine RECORD_AtomSelectorCommon
  !****************************************************************************************

  !****************************************************************************************
  subroutine RECORD_AtomSelectorCommonTool(Stamp, SimBox, CtrlParam)
  !***  DESCRIPTION: to putout the results to output file.
  !
  !    INPUT:  Stamp,       the record stamp
  !            SimBox,      the simulation boxs
  !            CtrlParam,   the control parameters for simulation
  !
  implicit none
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam
       !--- local variables


            call RECORD_AtomSelectorCommon(Stamp, SimBox, CtrlParam)
            return
  end subroutine RECORD_AtomSelectorCommonTool
  !****************************************************************************************

  !****************************************************************************************
  end module AtomSelectorCommon
