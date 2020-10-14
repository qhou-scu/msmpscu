 module EMBEDMENT_DIFFUSION2
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
 use MD_TYPEDEF_SimMDCtrl
 use MD_Globle_Variables
 use MD_SimboxArray
 implicit none

     integer::m_processid=0
     !$$--- the id of  filename get from SimMDCtrl for I/O.
     !$$    Refering the definition of f_others in SimMDCtrl.
     character(len=10),parameter, private::mp_FTAGI="&AUXF_DIFI"
     character(len=10),parameter, private::mp_FTAGO="&AUXF_DIFO"
     character(len=256),  private::m_OUTFILE =""             ! filename of output data

     !---
     integer, private::m_INITED = 0
     integer::m_NDATA = 0
     integer::m_ATYP(mp_mxGROUP)=0                                    !the type of atoms to be includeed in diffusion calculation
     !---
     real(KINDDF),dimension(:), allocatable, private::m_TIME
     !--- Average value of displacement moment
     real(KINDDF),dimension(:,:), allocatable, private::m_R2          ! average displacement square of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_XA          ! average displacement of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_YA          ! average displacement of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_ZA          ! average displacement of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_XX          ! average displacement square of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_XY          ! average displacement square of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_XZ          ! average displacement square of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_YY          ! average displacement square of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_YZ          ! average displacement square of clusters at times
     real(KINDDF),dimension(:,:), allocatable, private::m_ZZ          ! average displacement square of clusters at times
     real(KINDDF),dimension(:), allocatable, private::m_NP            !number of boxes wihtout bad data

 !_______________________________________________________________________________________
      !character(len=11),parameter, private::mp_FTAGI="&AUXF_DIFF"
      !type(EMBEDCtrlParam), private::m_EMBEDCtrlParam
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
            if(len_trim(m_OUTFILE) .LE.0 ) then
               write(*,fmt="(A)")          " MDPSCU Error: no output file for "//gm_ExeName(1:len_trim(gm_ExeName))// " is given."
               write(*,fmt="(A,A,A)")      "               add the keyword in SETUP  file: ", mp_FTAGO, " fname,"
               write(*,fmt="(A,A,A,A,A)")  "               or: add the keyword in ",mp_FTAGI, " file: ", mp_FTAGO, " fname"
               write(*,fmt="(' Process to be stopped')")
               stop
            else
               call CreateDataFolder_Globle_Variables(m_OUTFILE)
            end if
            call Add_PrintProcess(PrintControlParameters)

            m_INITED = IOR(m_INITED,1)

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
    use MD_Globle_Variables,only:Load_ExAnalyCtl_SimMDCtrl
    implicit none
    !--- dummy vaiables
    character*(*)             ::fname
    type(SimMDBox), intent(in)::SimBox
    type(SimMDCtrl)           ::CtrlParam

    !--- local variables
    integer::N, I, IS, LINE
    character*256::STR
    character*32::STRNUMB(mp_mxGROUP), KEYWORD
    type(StatementList), pointer::StatList

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
                              do I=1, N
                                 m_ATYP(I) = ISTR(STRNUMB(I))
                              end do
                              if(count(m_ATYP .gt.0) .le. 0) then
                                write(*,fmt="(A)")          " MDPSCU Error: atom type is required for diffusion calculation"
                                write(*,fmt="(A,A,A)")      "               but no atom type is available "
                                write(*,fmt="(A, BZI6)")    '               check control file at line:', LINE
                                write(*,fmt="(A, BZI6)")    '        Usage: &PROP_TYPE ty1, ty2...'
                                write(*,fmt="(A)")          '               Process to be stopped'
                                stop
                              end if
                  end select
              end do
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
   integer,         intent(in)::hFile
   type(SimMDBox),  intent(in)::SimBox
   type(SimMDCtrl), intent(in)::CtrlParam

  !----   Local variables
   integer::I

        !$$--- print out the control parameters
         write(hFile,fmt="(A)")     " !************ Embedment_diffusion module to be performed **********"
         write(hFile,fmt="(A)")     " !     Type of atoms included for diffusion calculation: "
         do I=1, count(m_ATYP .gt. 0)
              write(hFile,fmt="(A,100(I4,1X))")     " !     TYPE of atoms: ", m_ATYP(1:count(m_ATYP.gt.0))
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
            allocate(m_TIME(IT),                    &
                   m_R2(ntype,IT),                  &
                   m_XA(ntype,IT),                  &
                   m_YA(ntype,IT),                  &
                   m_ZA(ntype,IT),                  &
                   m_XX(ntype,IT),                  &
                   m_XY(ntype,IT),                  &
                   m_XZ(ntype,IT),                  &
                   m_YY(ntype,IT),                  &
                   m_YZ(ntype,IT),                  &
                   m_ZZ(ntype,IT),                  &
                   m_NP(IT))

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
            m_NP = 0.D0

            m_NDATA = 0
            m_INITED = IOR(m_INITED,2)
            return
  end subroutine Allocate_WorkingArray
 !**********************************************************************************

 !****************************************************************************************
 subroutine Clear(SimBox,CtrlParam)
 !***  DESCRIPTION: to allocate mempry for diffusion calculations
 !
 !
 !***  The modules included ******************************************
 implicit none
     type(SimMDBox) ::SimBox
     type(SimMDCtrl)::CtrlParam
     !--- Local variables

          if(allocated(m_TIME)) deallocate(m_TIME)
          if(allocated(m_XA))   deallocate(m_XA)
          if(allocated(m_YA))   deallocate(m_YA)
          if(allocated(m_ZA))   deallocate(m_ZA)
          if(allocated(m_XX))   deallocate(m_XX)
          if(allocated(m_XY))   deallocate(m_XY)
          if(allocated(m_XZ))   deallocate(m_XZ)
          if(allocated(m_YY))   deallocate(m_YY)
          if(allocated(m_YZ))   deallocate(m_YZ)
          if(allocated(m_ZZ))   deallocate(m_ZZ)
          if(allocated(m_NP))   deallocate(m_NP)
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
       type(MDRecordStamp) ,         intent(in):: Stamp
       type(SimMDBox), dimension(:), intent(in):: SimBox
       type(SimMDCtrl),              intent(in):: CtrlParam

  !--- Local variables
       integer::NP
       integer::I,J, K, NS,IT, NTYP,ITYP
       real(KINDDF)::CENTER0(3), DISM(3,3), TIME0
       real(KINDDF), dimension(:,:), allocatable::DIS0, DIS00
       real(KINDDF), dimension(:,:,:), allocatable::X0

       integer::IS=0, baddata, CURJOB=-1
       SAVE IS, DIS0, DIS00, X0, NS, NTYP, CURJOB, TIME0
 !----
             !--- release allocated memory if ITIME LE 0
             if(Stamp%ITime .le. 0) then
                !if(IAND(m_INITED,2) .gt. 0) then
                !  call Clear(SimBox(1), CtrlParam)
                !  m_INITED = 1
                !  IS       = 0
                !  TIME0    = Stamp%Time
                !end if
                return
             end if

             if(CURJOB .ne. Stamp%ITest) then
                CURJOB = Stamp%ITest
                IS     = 0
                TIME0  = Stamp%Time
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
             !--- Get the average of displacement of the substarte
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
                !DIS0(I,1:3) = 0.D0
                !do K=SimBox(I)%IPA(1), SimBox(I)%IPA(2) -1
                !   DIS0(I,1:3) = DIS0(I,1:3) + SimBox(I)%DIS(K,1:3)
                !end do
                !DIS0(I,1:3)=DIS0(I,1:3)/SimBox(I)%NA(1)
             end do

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
                      CENTER0(1:3) = CENTER0(1:3) + SimBox(I)%DIS(K,1:3)
                   end do
                   CENTER0 = CENTER0/(SimBox(I)%NA(ITYP))
                   if(IS.eq.1) then
                      X0(IT,I,1:3)= CENTER0(1:3)
                   end if
                   CENTER0(1:3) = CENTER0(1:3)-X0(IT,I,1:3)-(DIS0(I,1:3) - DIS00(I,1:3))

                   !--- the displacement matrix
                   do J=1, 3
                      do K=1, 3
                         DISM(J,K)=CENTER0(J)*CENTER0(K)
                      end do
                   end do

                   m_R2(IT,IS) = m_R2(IT,IS)+DISM(1,1)+DISM(2,2)+DISM(3,3)
                   m_XA(IT,IS) = m_XA(IT,IS)+CENTER0(1)
                   m_YA(IT,IS) = m_YA(IT,IS)+CENTER0(2)
                   m_ZA(IT,IS) = m_ZA(IT,IS)+CENTER0(3)
                   m_XX(IT,IS) = m_XX(IT,IS)+DISM(1,1)
                   m_XY(IT,IS) = m_XY(IT,IS)+DISM(1,2)
                   m_XZ(IT,IS) = m_XZ(IT,IS)+DISM(1,3)
                   m_YY(IT,IS) = m_YY(IT,IS)+DISM(2,2)
                   m_YZ(IT,IS) = m_YZ(IT,IS)+DISM(2,3)
                   m_ZZ(IT,IS) = m_ZZ(IT,IS)+DISM(3,3)

               END DO   !end the loop for samples
             END DO !end loop for type of atom
             m_TIME(IS) = Stamp%Time - TIME0
             m_NP(IS) = m_NP(IS)+NP

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
  character*64  ::fmt
  character*32  ::strnum
  character*256 ::fname
  integer       ::hFile,ERR
  real(KINDDF)  ::XA, YA, ZA, R2,XX,XY,XZ,YY,YZ,ZZ,SR2,SXX,SXY,SXZ,SYY,SYZ,SZZ
  integer       ::IT,IS

          do IT=1,count(m_ATYP .gt. 0)

             !--- determine output file name for type of atom
             write(strnum, fmt="(I)") m_ATYP(IT)
             strnum = adjustl(strnum)
             fname  = m_OUTFILE(1:len_trim(m_OUTFILE))//"_ATYP"//strnum(1:len_trim(strnum))

             call AvailableIOUnit(hFile)
             write(*, fmt="(A, 1x,I4, 1x, A)") "Diffusion parameters for atom type", m_ATYP(IT), &
                                               "to be output to: "//fname(1:len_trim(fname))

             if(CtrlParam%Restart .gt. 0) then
                open(hFile, file = fname, status='unknown',position='append')
             else
                open(hFile, file = fname, status='unknown')
                !--- determine output format for title
                fmt ="(12x,A6, 2x, A14,2x, (14(A14,2x)))"
                write(hFile,FMT=fmt) "I","TIME(s)","<X>","<Y>","<Z>","<R2>", "<XX>","<XY>","<XZ>", "<YY>","<YZ>","<ZZ>"
             end if
             !--- determine output format data
             fmt ="(12x,I6, 2x, 1pE14.5,2x, (14(1pE14.5,2x)))"
             do IS= 1, m_NDATA
                XA = m_XA(IT,IS)/m_NP(IS)
                YA = m_YA(IT,IS)/m_NP(IS)
                ZA = m_ZA(IT,IS)/m_NP(IS)
                R2 = m_R2(IT,IS)/m_NP(IS)
                XX = m_XX(IT,IS)/m_NP(IS) - XA*XA
                XY = m_XY(IT,IS)/m_NP(IS) - XA*YA
                XZ = m_XZ(IT,IS)/m_NP(IS) - XA*ZA
                YY = m_YY(IT,IS)/m_NP(IS) - YA*YA
                YZ = m_YZ(IT,IS)/m_NP(IS) - YA*ZA
                ZZ = m_ZZ(IT,IS)/m_NP(IS) - ZA*ZA

                !write(hFile, FMT =fmt) IS,m_TIME(IS)*CP_PS2S, R2,SR2,XX,SXX,XY,SXY,XZ,SXZ,YY,SYY,YZ,SYZ,ZZ,SZZ
                 write(hFile, FMT =fmt) IS,m_TIME(IS)*CP_PS2S, XA/SimBox%RR, YA/SimBox%RR, ZA/SimBox%RR,    &
                                                      R2/(SimBox%RR*SimBox%RR),                             &
                                                      XX/(SimBox%RR*SimBox%RR),XY/(SimBox%RR*SimBox%RR),XZ/(SimBox%RR*SimBox%RR), &
                                                      YY/(SimBox%RR*SimBox%RR),YZ/(SimBox%RR*SimBox%RR), &
                                                      ZZ/(SimBox%RR*SimBox%RR)
               end do
          end do
          close(hFile)
          !*** to release allocated memory
          call Clear(SimBox, CtrlParam )
   return
 end subroutine AfterRecord_Diffusion_Param
 !****************************************************************************************
 end module EMBEDMENT_DIFFUSION2

