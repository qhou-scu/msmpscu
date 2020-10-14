 module EMBEDMENT_DIFFUSION1
 !***  DESCRIPTION:
 !     This program is used to calculate the diffusion parameters of a given type of
 !     embedment cluster
 !
 !    Authored by HOU Qing
 !
 !    DEPENDENCY: module EMBEDMENT_COMMON_2010
 !
 !    HOSTORY:  2011-03 EMBEDMENT_DIFFUSION_2010
 !              2014-07 modified to compatibale with MD_SimBoxArray_ToolShell_14_GPU.F90
 !
 !    SEE ALSO:      MD_SimBoxArray_ToolShell_14_GPU.F90
 use MD_Globle_Variables
 use MD_SimboxArray
 use MD_TYPEDEF_SimMDCtrl
 use EMBEDMENT_TypeDef_CtrlParam_2010
 implicit none

     integer::m_processid=0
     !$$--- the id of  filename get from SimMDCtrl for I/O.
     !$$    Refering the definition of f_others in SimMDCtrl.
     character(len=10),parameter, private::mp_FTAGI="&AUXF_DIFI"
     character(len=10),parameter, private::mp_FTAGO="&AUXF_DIFO"
     character(len=256),  private::m_OUTFILE =""             ! filename of output data
     integer, private::m_hOUTFILE
     !---
     integer, private::m_INITED = 0
     integer::m_NDATA = 0
     integer::m_ATYP(10)=0                                            !the type of atoms to be includeed in diffusion calculation
     !---
     real(KINDDF),private::m_TIME
     !--- Average value of displacement moment
     real(KINDDF),dimension(:), allocatable, private::m_R2          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_XA          ! average displacement of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_YA          ! average displacement of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_ZA          ! average displacement of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_XX          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_XY          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_XZ          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_YY          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_YZ          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_ZZ          ! average displacement square of clusters at times

 !_______________________________________________________________________________________
     private::LoadControlParameters, PrintControlParameters
 contains

  !*********************************************************************************
  subroutine Initialize(SimBox, CtrlParam)
    !***  PURPOSE:   to intialize the parameters to be used in diffusion calculations
    !     INPUT:     SimBox: the simulation box
    !                CtrlParam: the control parameters
    !
    !     OUTPUT:   the working spaces allocated
    !
     use MD_TYPEDEF_PrintList, only:Add_PrintProcess
    implicit none
      !--- dummy vaiables
      type(SimMDBox), intent(inout)::SimBox
      type(SimMDCtrl)              ::CtrlParam

      !--- Local variables
      integer::I, IFILE

           if(m_INITED .gt. 0) then
             call Clear(SimBox, CtrlParam)
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
               write(*,*) "MDPSCU Error: no control file for diffusion calculation is given."
               write(*,*) "              add the keyword in SETUP file:",  mp_FTAGI
               write(*,fmt="(' Process to be stopped')")
               stop
            end if

            write(*,fmt="(A)") " !**** Loading control data from: "//CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
            if(gm_hFILELOG .gt. 0) write(gm_hFILELOG,fmt="(A)") " !**** Loading control data from: "// &
                                                            CtrlParam%f_others(IFILE)(1:len_trim(CtrlParam%f_others(IFILE)))
            call LoadControlParameters(CtrlParam%f_others(IFILE), SimBox, CtrlParam)
            call Add_PrintProcess(PrintControlParameters)

            CtrlParam%TIMELOOPOUT = 1
            m_INITED = IOR(m_INITED,1)

            !$$--- prepare the output file
            call AvailableIOUnit(m_hOUTFILE)
            open(m_hOUTFILE, file = m_OUTFILE, status='unknown')
            print *, "Diffusion parameters to be output to: ",m_OUTFILE

             !--- determine output format
             write(m_hOUTFILE,FMT="(12x,A6, 2x, A14,2x, 10(10(A14,2x)))") "I","TIME(s)",&
                                      ("<X>","<Y>","<Z>","<R2>", "<XX>","<XY>","<XZ>", "<YY>","<YZ>","<ZZ>", I=1, count(m_ATYP .gt. 0))


           return
    end subroutine Initialize
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
    integer::hFile, N, I, J, LINE
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

                         !*** to get the properties to be included in clustering
                         case("&PROP_TYPE")
                              call EXTRACT_NUMB(STR,SimBox%nGroup,n,STRNUMB)
                              if(N.GT.size(m_ATYP)) then
                                 write(*,fmt="(A,I4,A)")  ' MDPSCU Error: the max types included in clustering is ',size(m_ATYP)
                                 write(*,fmt="(A,I4,A)")  '               but ', N, ' types are required in clustering'
                                 write(*,fmt="(A, BZI6)") '               check control file at line:', LINE
                                 write(*,fmt="(A)")       '               Process to be stopped'
                                 stop
                              end if

                              do I=1, N
                                 m_ATYP(I) = ISTR(STRNUMB(I))
                              end do
                              call AddDataProKWD_SimMDBox(SimBox, "TYPECOL")

                         case default
                              write(*,*)" MDPSCU Warning: unknown keyword in diffusion control file", KEYWORD(1:LEN_TRIM(KEYWORD))
                              write(*,fmt="('               check control file at line:', BZI6)") LINE
                              call ONWARNING(gm_OnWarning)
                  end select
              end do
    100     close(hFile)

            !$$--- check input consistent
             if(count(m_ATYP .gt.0) .le. 0) then
                write(*,fmt="(A)")          " MDPSCU Error: no type of atoms is given for diffusion calculation"
                write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", "&PROP_TYPE", " type IDs of atoms"
                write(*,fmt="(' Process to be stopped')")
                stop
             end if

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

    200     write(*,fmt="(' MDPSCU Error: fail to open control file for diffusion calculations')")
            write(*,fmt="('               check the existence of file: ', A)") fname(1:len_trim(fname))
            write(*,fmt="(' Process to be stopped')")
            stop
            return
   end subroutine LoadControlParameters
  !*********************************************************************************

  !**********************************************************************************
   subroutine PrintControlParameters(hFile, SimBox, CtrlParam)
  !***  PURPOSE:  to print out the control parameters for this module
  !     INPUT:    hFile, SimBox, CtrlParam
  !     OUTPUT:
  !

   use MD_CONSTANTS
   use MD_TYPEDEF_SimMDBox
   use MD_TYPEDEF_SimMDCtrl
   implicit none
  !----   DUMMY Variables
   integer,          intent(in)::hFile
   type(SimMDBox),   intent(in)::SimBox
   type(SimMDCtrl),  intent(in)::CtrlParam

  !----   Local variables
   integer::I
   character*32::Tag

        !$$--- print out the control parameters
         write(hFile,fmt="(A)")     " !************ Embedment_diifusion module to be performed **********"
         write(hFile,fmt="(A)")     " !     Type of atoms included for diffusion calculation: "
         do I=1, NumberofData_DataPad(SimBox%proKWDList)
            call Tag_DataPad(I, SimBox%proKWDList, Tag)
            select case(Tag)
                   case ("&TYPECOL")
                         write(hFile,fmt="(A,10(I4,1X))")     " !     TYPE of atoms: ", m_ATYP(1:count(m_ATYP.gt.0))

            end select
         end do
         write(hFile,fmt="(A)")  " !     "

         return
   end subroutine PrintControlParameters
  !**********************************************************************************

  !**********************************************************************************
  subroutine Allocate_WorkingArray(SimBox, CtrlParam)
  !***  PURPOSE:  to allocate and initialize the device working memories
  !
  !     INPUT:
  !     OUTPUT:
  !
  implicit none
      !--- dummy vaiables
      type(SimMDBox),  intent(in)::SimBox
      type(SimMDCtrl), intent(in)::CtrlParam

      integer::ntype, IT, IT1


     !*** allocate meomery storing the displacements
          ntype = max(count(m_ATYP.gt.0),1)
          call NumberCfg_SimMDCtrl(CtrlParam, IT)
          call NumberRec_SimMDCtrl(CtrlParam, IT1)
          IT = max(IT, IT1)
          allocate(m_R2(ntype),   &
                   m_XA(ntype),   &
                   m_YA(ntype),   &
                   m_ZA(ntype),   &
                   m_XX(ntype),   &
                   m_XY(ntype),   &
                   m_XZ(ntype),   &
                   m_YY(ntype),   &
                   m_YZ(ntype),   &
                   m_ZZ(ntype))

            m_TIME = 0.D0
            m_R2 = 0.D0
            m_XA = 0.D0
            m_YA = 0.D0
            m_ZA = 0.D0
            m_XX = 0.D0
            m_XY = 0.D0
            m_XZ = 0.D0
            m_YY = 0.D0
            m_YZ = 0.D0
            m_ZZ = 0.D0

            m_NDATA = 0
            m_INITED = IOR(m_INITED,2)
            return
  end subroutine Allocate_WorkingArray
  !**********************************************************************************


 !***********************************************************************************
 subroutine Clear(SimBox,CtrlParam)
 !***  DESCRIPTION: to allocate mempry for diffusion calculations
 !
 !
 !***  The modules included ******************************************
 implicit none
     type(SimMDBox) ::SimBox
     type(SimMDCtrl)::CtrlParam
     !--- Local variables

          if(allocated(m_XA))   deallocate(m_XA)
          if(allocated(m_YA))   deallocate(m_YA)
          if(allocated(m_ZA))   deallocate(m_ZA)
          if(allocated(m_XX))   deallocate(m_XX)
          if(allocated(m_XY))   deallocate(m_XY)
          if(allocated(m_XZ))   deallocate(m_XZ)
          if(allocated(m_YY))   deallocate(m_YY)
          if(allocated(m_YZ))   deallocate(m_YZ)
          if(allocated(m_ZZ))   deallocate(m_ZZ)
          if(allocated(m_R2))   deallocate(m_R2)

          m_INITED = 0
          m_NDATA  = 0
          m_ATYP   = 0

     RETURN
 END subroutine Clear
 !****************************************************************************************

 !****************************************************************************************
 subroutine Record_Diffusion_Param(Stamp, SimBox, CtrlParam )
 !***  PORPOSE: to get the migration distance  of clusters
 !     INPUT:  SimBox, the  samples
 !     OUTPUT: RAV, RMI, RMA
 implicit none
  !--- input
       type(MDRecordStamp)                     :: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

      !--- Local variables
        integer::NP
        integer::I,J, K, NS,IT, NTYP,ITYP
        real(KINDDF)::CENTER0(3), AVDIS(3), DISM(3,3), TIME0
        real(KINDDF), dimension(:,:), allocatable::DIS0, DIS00
        real(KINDDF), dimension(:,:,:), allocatable::X0
        real(KINDDF)::LATT
        integer::IS=0, baddata, CURJOB=-1
        SAVE IS, DIS0, DIS00, X0, NS, NTYP, CURJOB, TIME0
         !----
             !--- release allocated memory if ITIME LE 0
             if(Stamp%ITime .LE. 0) then
                if(iand(m_INITED,2) .gt. 0) then
                  call Clear(SimBox(1), CtrlParam)
                  m_INITED = 1
                  IS = 0
                  TIME0 = Stamp%Time
                end if
                return
             end if

             if(CURJOB .ne. Stamp%ITest) then
                CURJOB = Stamp%ITest
                IS = 0
                TIME0 = Stamp%Time
              end if

             if(iand(m_INITED,2) .eq. 0) then
                call  Allocate_WorkingArray(SimBox(1), CtrlParam)
                if(allocated(DIS0)) deallocate(DIS0)
                if(allocated(DIS00)) deallocate(DIS00)
                if(allocated(X0)) deallocate(X0)
                NTYP = count(m_ATYP.gt.0)
                NS =size(SimBox)
                allocate(X0(NTYP,NS,3))
                allocate(DIS0(NS,3), DIS00(NS,3))
                CURJOB = Stamp%ITest
                IS = 0
                TIME0 = Stamp%Time
             end if

             IS = IS + 1
             m_NDATA = IS

             NP = 0
             m_XA = 0.D0
             m_YA = 0.D0
             m_ZA = 0.D0
             m_XX = 0.D0
             m_XY = 0.D0
             m_XZ = 0.D0
             m_YY = 0.D0
             m_YZ = 0.D0
             m_ZZ = 0.D0
             !--- Get the average of displacement of the substrate
             do I =1, NS
                baddata = 0
                do J=1, size(SimBox(I)%STATU)
                   if(IAND(SimBox(I)%STATU(J),CP_STATU_HWORD) .ne. 0) then
                      baddata = 1
                      write(*,*) "BAD DATA", I, J, (I-1)*size(SimBox(I)%STATU)+J, SimBox(I)%STATU(J)
                      exit
                   end if
                end do
                if(baddata .gt. 0) cycle

                NP = NP+1
                !---  get the displacement of the substrate
                call Get_Displacement_SimMDBox(SimBox(I),DIS0(I,1:3))
                if(IS .eq. 1) then
                   DIS00(I,1:3) = DIS0(I, 1:3)
                end if
             end do

             !--- Get the average of displacement of the clusters
             do IT=1, NTYP
                ITYP = m_ATYP(IT)
                !--- Get the displacement of the cluster at start point
                do I =1, NS
                   if(any(IAND(SimBox(I)%STATU,CP_STATU_HWORD) .ne. 0)) then
                      cycle
                   end if

                   !--- to get migration distance of the center of the cluster
                   CENTER0 = 0.D0
                   do K=SimBox(I)%IPA(ITYP), SimBox(I)%IPA(ITYP+1)-1
                      CENTER0(1:3) = CENTER0(1:3) + SimBox(I)%DIS(K,1:3)
                   end do
                   CENTER0 = CENTER0/(SimBox(I)%NA(ITYP))

                   if(IS.eq.1) then
                      X0(IT,I,1:3)= CENTER0(1:3)
                   end if
                   m_XA(IT) = m_XA(IT)+CENTER0(1)-X0(IT,I,1)-(DIS0(I,1) - DIS00(I,1))
                   m_YA(IT) = m_YA(IT)+CENTER0(2)-X0(IT,I,2)-(DIS0(I,2) - DIS00(I,2))
                   m_ZA(IT) = m_ZA(IT)+CENTER0(3)-X0(IT,I,3)-(DIS0(I,3) - DIS00(I,3))
                end do
             end do
             m_XA = m_XA/dble(NP)
             m_YA = m_YA/dble(NP)
             m_ZA = m_ZA/dble(NP)

             !--- Get the average standard devition of displacement of the clusters
             do IT=1, NTYP
                ITYP = m_ATYP(IT)
                !--- Get the displacement of the cluster
                do I =1, NS
                   if(any(IAND(SimBox(I)%STATU,CP_STATU_HWORD) .ne. 0)) then
                      cycle
                   end if
                   !--- to get migration distance of the center of the cluster
                   CENTER0 = 0.D0
                   do K=SimBox(I)%IPA(ITYP), SimBox(I)%IPA(ITYP+1) -1
                      CENTER0(1:3) = CENTER0(1:3) + SimBox(I)%DIS(K,1:3) !-DIS0(I,1:3)
                   end do
                   CENTER0 = CENTER0/(SimBox(I)%NA(ITYP))
                   CENTER0(1:3) = (CENTER0(1:3) - X0(IT,I,1:3))-(DIS0(I,1:3) - DIS00(I,1:3))
                   CENTER0(1)   =  CENTER0(1)   - m_XA(IT)
                   CENTER0(2)   =  CENTER0(2)   - m_YA(IT)
                   CENTER0(3)   =  CENTER0(3)   - m_ZA(IT)

                   !--- the displacement matrix
                   do J=1, 3
                      do K=1, 3
                         DISM(J,K)=CENTER0(J)*CENTER0(K)
                      end do
                   end do

                   m_XX(IT) = m_XX(IT) + DISM(1,1)
                   m_XY(IT) = m_XY(IT) + DISM(1,2)
                   m_XZ(IT) = m_XZ(IT) + DISM(1,3)
                   m_YY(IT) = m_YY(IT) + DISM(2,2)
                   m_YZ(IT) = m_YZ(IT) + DISM(2,3)
                   m_ZZ(IT) = m_ZZ(IT) + DISM(3,3)
               end do   !end the loop for samples
               m_R2(IT) = m_XX(IT) + m_YY(IT) + m_ZZ(IT)

             end do !end loop for type of atom
             m_TIME = Stamp%Time - TIME0
             m_XX   = m_XX/dble(NP)
             m_XY   = m_XY/dble(NP)
             m_XZ   = m_XZ/dble(NP)
             m_YY   = m_YY/dble(NP)
             m_YZ   = m_YZ/dble(NP)
             m_ZZ   = m_ZZ/dble(NP)
             m_R2   = m_R2/dble(NP)

             !$$--- write out the  result
             LATT = SimBox(1)%RR
             write(m_hOUTFILE, FMT = "(12x,I6, 2x, 1pE14.5,2x, 10(10(1pE14.5,2x)))" )                        &
                                           IS, m_TIME*CP_PS2S,                                               &
                                           (m_XA(IT)/LATT, m_YA(IT)/LATT, m_ZA(IT)/LATT,                     &
                                            m_R2(IT)/(LATT*LATT),                                            &
                                            m_XX(IT)/(LATT*LATT), m_XY(IT)/(LATT*LATT), m_XZ(IT)/(LATT*LATT),&
                                            m_YY(IT)/(LATT*LATT), m_YZ(IT)/(LATT*LATT),                       &
                                            m_ZZ(IT)/(LATT*LATT),  IT=1, NTYP)


   return
 end subroutine Record_Diffusion_Param
 !****************************************************************************************

 !****************************************************************************************
 subroutine AfterRecord_Diffusion_Param(SimBox, CtrlParam )
 !***  PORPOSE: to print out diffusion parameters
 !     INPUT:  SimBox, the  samples
 implicit none
  !--- input
  type(SimMDBox) ::SimBox
  type(SimMDCtrl)::CtrlParam

  !--- Local variables
  character*16::str(10)=""
  character*64::fmt
  character*32::strnum
  integer::hFile,ERR
  logical::opened, EX
  real(KINDDF)::XA, YA, ZA, R2,XX,XY,XZ,YY,YZ,ZZ,SR2,SXX,SXY,SXZ,SYY,SYZ,SZZ
  INTEGER::IT,IS, NTYP

          !*** to release allocated memory
          if(m_hOUTFILE .gt. 0) close(m_hOUTFILE)
          m_hOUTFILE = 0
          call Clear(SimBox, CtrlParam )
          return
   return
 end subroutine AfterRecord_Diffusion_Param
 !****************************************************************************************
 end module EMBEDMENT_DIFFUSION1

