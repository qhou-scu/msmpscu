!***********************************************************************
   module  STOP_SRIM_MODULE
   !***  this module is to obtain the stopping power table from SRIM data.
   !
   !     HISTORY:     adopted by Cui Jiechao, Apr, 2018
   !
   use MD_CONSTANTS
   use MiniUtilities
   implicit none

       character*4, dimension(100), parameter, private:: Elem  = (/                                         &
        " H " ,     "He" ,     "Li" ,    "Be" ,    "B " ,    "C "  ,    "N " ,    "O "  ,   "F " ,   "Ne" ,  &
        " Na" ,     "Mg" ,     "Al" ,    "Si" ,    "P " ,    "S "  ,    "Cl"  ,   "Ar"  ,   "K " ,   "Ca" ,  &
        " Sc" ,     "Ti" ,     "V " ,    "Cr" ,    "Mn" ,    "Fe"  ,    "Co"  ,   "Ni"  ,   "Cu" ,   "Zn" ,  &
        " Ga" ,     "Ge" ,     "As" ,    "Se" ,    "Br" ,    "Kr"  ,    "Rb"  ,   "Sr"  ,   "Y " ,   "Zr" ,  &
        " Nb" ,     "Mo" ,     "Tc" ,    "Ru" ,    "Rh" ,    "Pd"  ,    "Ag"  ,   "Cd"  ,   "In" ,   "Sn" ,  &
        " Sb" ,     "Te" ,     " I" ,    "Xe" ,    "Cs" ,    "Ba"  ,    "La"  ,   "Ce"  ,   "Pr" ,   "Nd" ,  &
        " Pm" ,     "Sm" ,     "En" ,    "Gd" ,    "Tb" ,    "Dy"  ,    "Ho"  ,   "Er"  ,   "Tm" ,   "Yb" ,  &
        " Lu" ,     "Hf" ,     "Ta" ,    " W" ,    "Re" ,    "Os"  ,    "Ir"  ,   "Pt"  ,   "Au" ,   "Hg" ,  &
        " Tl" ,     "Pb" ,     "Bi" ,    "Po" ,    "At" ,    "Rn"  ,    "Fr"  ,   "Ra"  ,   "Ac" ,   "Th" ,  &
        " Pa" ,     "U " ,     "Np" ,    "Pu" ,    "Am" ,    "Cm"  ,    "Bk"  ,   "Cf"  ,   "Es" ,   "Fm"  /)

       integer,          dimension(mp_mxGROUP*mp_mxGROUP),private :: m_IonNum, m_TargetNum
       real(KINDDF),     dimension(mp_mxGROUP*mp_mxGROUP),private :: m_Ion_Mass,m_target_Mass
       character*256,    dimension(mp_mxGROUP*mp_mxGROUP),private :: m_table_Name
       character*32,     dimension(mp_mxGROUP*mp_mxGROUP),private :: m_IonName,m_TargetName
       real(KINDDF),     dimension(:),     allocatable,   private :: m_IonE
       real(KINDDF),     dimension(:,:),   allocatable,   private :: m_ST_E

       private:: read_folders,        &
                 Load_SRIM_Data,      &
                 Export_mySTPTable,   &
                 line_no

      public::  Creat_SRIM_STPTable
   contains
   !*******************************************************************
   subroutine Creat_SRIM_STPTable(fname,fname1)
   !***  PURPOSE: to Creat STPTable from SRIM data.
   !
   !     INPUT:  fname,  the filename to load the input/output file' names from
   !
   !     OUTPUT: fname1, the filename to export the table to
   !
   implicit none

         !----   DUMMY Variables
         character*256::fname,fname1

         ! local variables
         logical::EX
         integer::i
         integer::lineN,hFile,lineTab
         integer::table_Num
         real(KINDDF):: E_min, E_max

         call read_folders(table_Num, m_table_Name,fname,fname1)
         E_min = -1
         E_max =  0
         do i =1, table_Num
            inquire(FILE=m_table_Name(i), EXIST=EX)
            if(.not. EX) then
                write(*,fmt="(A)") ' MDPSCU Error: cannot find required external SRIM_STPTable file: '
                write(*,fmt="(A)") '               '//m_table_Name(i)(1:len_trim(m_table_Name(i)))
                write(*,fmt="(A)") '               process to be stop'
                stop
            end if
            call AvailableIOUnit(hFile)
            open (UNIT = hFile , FILE = m_table_Name(i))
            call  line_no(hFile, lineN)
            lineTab=lineN-32
            if (i.eq.1)then
               allocate(m_IonE(lineTab),m_ST_E(lineTab,table_Num))
            end if
            call  Load_SRIM_Data(hFile,lineTab,m_IonName(i),m_TargetName(i),m_IonNum(i),m_TargetNum(i),m_Ion_Mass(i),&
                  m_target_Mass(i),m_IonE,m_ST_E(:,i))

            if (E_min.lt.0)then
                E_min = minval(m_IonE)
            else
               if(E_min.ne.minval(m_IonE))then
                   write(*,*) "MDPSCU Error: The minimum energies of different SRIM data are not equal"
                   write(*,*) "              Please check SRIM data"
                   write(*,*) "              Process to be stopped"
                   stop
               end if
            end if

            if (E_max.le.0)then
                E_max = maxval(m_IonE)
            else
                if(E_max.ne.maxval(m_IonE))then
                   write(*,*) "MDPSCU Error: The maximum energies of different SRIM data are not equal"
                   write(*,*) "              Please check SRIM data"
                   write(*,*) "              Process to be stopped"
                   stop
               end if
            end if

           close(hFile)
         end do

         call Export_mySTPTable(fname1,lineTab,table_Num)

         if (allocated(m_IonE)) deallocate(m_IonE)
         if (allocated(m_ST_E)) deallocate(m_ST_E)
         return
   end  subroutine Creat_SRIM_STPTable
  !*******************************************************************************************
  !*******************************************************************************************
   subroutine read_folders( count, folders,fname, OUT_name)
   !***  PURPOSE: to obtain files with SRIM STPtable type
   !     INPUT:   fname,   the filename to load the input/output files' name from
   !
   !     OUTPUT:  count,    the number of files with SRIM STPtable type
   !              folders,  the name of files with SRIM STPtable type
   !              OUT_name, the name of created STPtable file
   !***
  implicit none

          !---dummy variables
          integer ,intent(out)::  count
          character*(*), dimension(:), intent(out) :: folders
          character*256,intent(in) :: fname
          character*256,intent(out):: out_name

          !---local vairables
          integer :: hFile, line, N
          character*256 :: STR,KEYWORD

          call AvailableIOUnit(hFile)
          open (UNIT = hFile , FILE = fname,STATUS='OLD' )
          count = 0

         do while (.true.)
            call GetInputStrLine(hFile,STR, LINE, "!", *100)
            call GetKeyWord("&", STR, KEYWORD)
            call UpCase(KEYWORD)

            select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
               case( "&INPUT_STN")
                  count = count +1
                  if (count .gt.size(folders))  then
                     write(*,*) "the number of folders is larger than the given size!"
                     stop
                  end if
                  call Extract_Substr (STR,1 , N, folders(count))
                  folders(count) = adjustl(folders(count))

               case( "&OUTPUT_STN")
                  call Extract_Substr (STR,1 , N, OUT_name)
                  if(N .gt.0)then
                      OUT_name = adjustl(OUT_name)
                      OUT_name = OUT_name(1:len_trim(OUT_name))
                   else
                     write(*,fmt="(A)") ' MDPSCU Error: cannot find required SRIM_STPTable output file '
                     write(*,fmt="(A)") '               Please check the SRIM input file: '//fname(1:len_trim(fname))
                     write(*,fmt="(A)") '               process to be stop'
                     stop
                   end if
               end select
           end do

           close (hFile)
           return
     !-----------------------------------------------------------------------------------------------------
       100 rewind(hFile)

  end subroutine read_folders
  !************************************************************************************************
  !**************************************************************************************************
   subroutine Load_SRIM_Data(hFile, linenum,IonName,TargetName,IonNum,TargetNum,Ion_Mass,target_Mass, &
                           IonE,ST_E)
   !***  PURPOSE: to load data with SRIM STPtable type
   !
   !     INPUT:   hFile,       the input unit
   !              linenum,     the stoping table size
   !
   !     OUTPUT:  IonName,     the atom name of ion
   !              TagertName,  the atom name of target
   !              IonNum,      the atom number of ion
   !              TargetNum,   the atom number of target
   !              Ion_Mass,    the atom mass of ion
   !              target_Mass, the atom mass of target
   !              IonE,        the energy points
   !              ST_E,        the stoping power table
   !***
   implicit none
          !---dummy variables
          integer :: hFile
          integer :: linenum
          integer :: IonNum,TargetNum
          character*32 :: IonName, TargetName
          real(KINDDF) :: Ion_Mass,target_Mass
          real(KINDDF),dimension(:)::IonE, ST_E

          !---local vairables
          character*256 :: STR,STRTMP,STRTMP1(4)
          character*32  :: STRNUMB(4)
          real(KINDDF),dimension(:),   allocatable :: ST_N
          real(KINDDF),dimension(:,:), allocatable :: RANGE
          character*32,dimension(:),   allocatable :: Unit_E
          character*32,dimension(:,:), allocatable :: UNIT_RANGE

          integer     ::J,n, LINE
          real(KINDDF)::Target_Den1,Target_Den2,MultiCoeff

          allocate(ST_N(linenum),RANGE(linenum,3))
          allocate(Unit_E(linenum),UNIT_RANGE(linenum,3))

          MultiCoeff = 0

          do j = 1, linenum+32
             if (J.LT.20)then
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                if (J.eq.6)then
                   STR = adjustl(STR)
                   call Extract_Numb(STR,2,n,STRNUMB)
                   IonNum   = DRSTR(STRNUMB(1))
                   Ion_Mass = DRSTR(STRNUMB(2))
                else if (J.EQ.7)THEN
                   call Extract_Numb(STR,2,n,STRNUMB)
                   Target_Den1 = DRSTR(STRNUMB(1))
                   Target_Den2 = DRSTR(STRNUMB(2))
                   target_Mass = Target_Den1/Target_Den2/CP_AU2G
                else if(J.EQ. 12)THEN
                   STR = adjustl(STR)
                   call Extract_Numb(STR,1,n,STRNUMB)
                   TargetNum = DRSTR(STRNUMB(1))
                end if
             else if(J.GE.20 .AND. J.LE. linenum+19)then
                read(hFile,*) IonE(J-19),Unit_E(J-19),ST_E(J-19),ST_N(J-19),RANGE(J-19,1),UNIT_RANGE(J-19,1), &
                             RANGE(J-19,2),UNIT_RANGE(J-19,2), RANGE(J-19,3),UNIT_RANGE(J-19,3)

                ST_E(J-19)   = ST_E(J-19)*CP_KEV2MEV*1000*Target_Den1/Target_Den2

                Unit_E(J-19) = adjustl(Unit_E(J-19))
                select case(Unit_E(J-19))
                       case("eV")
                            IonE(J-19) = IonE(J-19)*CP_KEV2EV
                       case("MeV")
                            IonE(J-19) = IonE(J-19)*CP_KEV2MEV
                end select
             else
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                if (j.GE.linenum+23.AND.j.LE.linenum+30)THEN
                   call GetNthKeyword(STR,4,STRTMP1(2))
                   call GetNthKeyword(STR,2,STRTMP1(1))
                   STRTMP = STRTMP1(1)(1:LEN_TRIM(STRTMP1(1)))//"/"//STRTMP1(2)(1:LEN_TRIM(STRTMP1(2)))
                   IF (STRTMP.EQ."MeV/(mg/cm2)") THEN
                      call Extract_Numb(STR,1,n,STRNUMB)
                      MultiCoeff = DRSTR(STRNUMB(1))
                   END IF
               END IF
              end if
           end do
           close(hFile)

           ST_E       = ST_E*MultiCoeff
           IonName    = Elem(IonNum)
           TargetName = Elem(TargetNum)
           if (allocated(ST_N))       deallocate(ST_N)
           if (allocated(RANGE))      deallocate(RANGE)
           if (allocated(Unit_E))     deallocate(Unit_E)
           if (allocated(UNIT_RANGE)) deallocate(UNIT_RANGE)
           return
    !-----------------------------------------------------------------------
      100  write(*,*)"The process to be stopped."
           stop



   end subroutine Load_SRIM_Data
   !*********************************************************************************

   !***************************************************************************************
  subroutine Export_mySTPTable(FNAME,n,T_Num)
  !***  PURPOSE:   to export the generated ST table
  !     INPUT:
  !                FNAME,     the filename to export the table to
  !                n,         the stoping table size
  !                T_Num,     the number of stoping table size
  !      OUTPUT:
  !
  !***************************************************************************************
  implicit none
        !---dummy variables
        character*256::fname
        integer::n,T_Num

        !local variables
        integer::I,J, K, hFile
        character*256::FILE
        character*32::STAB(mp_mxGROUP)

        FILE = FNAME(1:len_trim(FNAME))//".stp"
        call AvailableIOUnit(hFile)
        open(FILE=FILE, UNIT=hFile)

        !$$--- write a indentify header
        write(*,fmt="(A,A)")   ' MDPSCU Message: stopping table save to:', FILE(1:len_trim(FILE))
        write(hFile, fmt="(A)") "&MDPSCU_STPTAB.stp"
        write(hFile, fmt="(A, I7, A, 100I5)") "! In the energy-stopping relation table"
        write(hFile, fmt="(A, I7, A, 100I5)") "! energy  in the unit of keV"
        write(hFile, fmt="(A, I7, A, 100I5)") "! stoping in the unit of keV*cm^2"
        write(hFile, fmt="(A, I7, A, 100I5)") "! "
        write(hFile, fmt="(A, I7, A, 100I5)") "&NUMTABLE ", T_Num
        write(hFile, fmt="(A, I7)")           "&NUMPOINT ", n

        do I = 1, T_Num
           m_IonName(I)     = adjustl( m_IonName(I))
           m_TargetName(I)  = adjustl(m_TargetName(I))
           STAB(I)          = m_IonName(I)(1:len_trim(m_IonName(I)))//"->"//m_TargetName(I)(1:len_trim(m_TargetName(I)))
           write(hFile, fmt="(A, A I4, 1X, 4(F9.2,1X))") "&"//STAB(I)(1:10), " with COL# ", I+1, &
                                                         1.0*m_IonNum(I), m_Ion_Mass(I),1.0*m_TargetNum(I), m_target_Mass(I)
        end do


        write(hFile, fmt="(A16,3x, 64(A16))")  "!--- ENERGY     ", STAB(1:T_Num)
        do J=1, n
           write(hFile, fmt="(1x, 40(1PE16.8))")m_IonE(J),(m_ST_E(J,I), I=1,T_Num)
        end do

       close(hFile)
       return
  end subroutine Export_mySTPTable
  !***************************************************************************************

  !********************************************************************************************

  subroutine line_no(m_hFile,n)
  !***  PURPOSE:: to obtain  the line number of a file
  implicit none
       Integer , Intent( in ) :: m_hFile
       Integer , Intent( out ) :: n
       character*256:: str
        n=0
        rewind(m_hFile)
        Do
          Read( m_hFile, * , End = 999 , Err = 999 ) str
           n = n + 1
        End Do
        999 Rewind( m_hFile )

        return
 end subroutine line_no
!**************************************************************************************************************

 end module STOP_SRIM_MODULE
