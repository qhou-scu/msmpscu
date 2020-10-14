  module DEPOSITION_TypeDef_CtrlParam_2010
  !***  DESCRIPTION: this module is to define the data type for Depostion Control Parameter.
  !
  !                  ______________________________________________________
  !                  HOU Qing, Nov, 2010
  use MD_CONSTANTS
  use MiniUtilities
  implicit none

   !--------------------------------------------------
   !--- The type describes the data controlling the deposition
   !    simulation
      integer,parameter::CP_DEP_STYPE_EXT = 0          !external depsotion
      integer,parameter::CP_DEP_STYPE_INT = 1          !internal depsotion
      integer,parameter::CP_DEP_STYPE_INT_SIA = 2      !internal SIA depsotion

      integer,parameter,private::CP_ORIENT_STYLE_MONODIR   =  0  !the incident orientation is in given direction
      integer,parameter,private::CP_ORIENT_STYLE_RANDOM4PI = 1   !the incident orientation is random in 4PI solid angle
      integer,parameter,private::CP_ORIENT_STYLE_LATTDIR   = 2   !the incident orientation is in given lattice direction

      integer,parameter,private::m_MaxType = 8       !maxmun type of incident particles
      type DepositionCtrlParam
           !--- simulation parameters

           integer::NPRT     = 0                    !number of incident particles to be simulated
           integer::DEPSTYLE = CP_DEP_STYPE_EXT     !deposition style =0 , external; =1 internal

           integer::nType=0                         !number of types of incident atom
           integer::aType(m_MaxType)=0              !type of incident atom
           real(KINDDF)::Percent(m_MaxType)         !percentage for each kind of atoms

           integer::ENTYPE(m_MaxType)  = 0          !type of energy spectrum =0 for monoenergy
           real(KINDDF)::Energy(m_MaxType)=0.D0     !incident energy for each kind of atoms

           integer::ANTYPE(m_MaxType)  = 0          !type of angular distribution
                                                    !  =0 for mono-orientation
                                                    !  =1 in random direction selected in 4PI solid angle
                                                    !  =2,in a Miller direction defined by MILLER (see following)

           real(KINDDF)::Orient(m_MaxType)=0.D0     !incident orientations, used for external case
                                                    !used if ANTYPE=1

           real(KINDDF)::MILLER(m_MaxType,3)=0      ! the Miller index lattice direction that the energtic atoms directed in.
                                                    ! used if ANTYPE =2
      end type DepositionCtrlParam
  !--------------------------------------------------

  contains
  !*********************************************************************

  !****************************************************************************
  subroutine Load_Parameter_DepositionCtrlParam(hFile, CtrlParam)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(DepositionCtrlParam)::CtrlParam

     !--- local variables
      character*256::STR,STRTMP(1)=""
      character*32::STRNUMB(30)
      integer::I, J, N, LINE
     !----

          !*** To get number of incident atoms
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,2,n,STRNUMB)
              CtrlParam%NPRT = ISTR(STRNUMB(1))
              if(n.ge.2) then
                 CtrlParam%DEPSTYLE = ISTR(STRNUMB(2))
              end if

         !*** To get number of types of incident atoms
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,1,n,STRNUMB)
              CtrlParam%nType = ISTR(STRNUMB(1))

              if(CtrlParam%nType .LE. 0) return

         !*** LOOP FOR ALL TYPE OF PARTICLES
              DO I=1, CtrlParam%nType
                 !*** To get types of incident atoms
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 call Extract_Numb(STR,1,n,STRNUMB)
                 CtrlParam%aType(I) = ISTR(STRNUMB(1))

                 !*** To get percentage of incident atoms
                  call GetInputStrLine(hFile,STR, LINE, "!", *100)
                  call Extract_Numb(STR,1,n,STRNUMB)
                  CtrlParam%Percent(I) = DRSTR(STRNUMB(1))

                 !*** To get type  of incident energy spectrum
                  call GetInputStrLine(hFile,STR, LINE,  "!", *100)
                  call Extract_Numb(STR,1,n,STRNUMB)
                  CtrlParam%ENTYPE(I) = ISTR(STRNUMB(1))

                  if(CtrlParam%ENTYPE(I).LE. 1) then
                     !--- mono-energy depostion
                     call GetInputStrLine(hFile,STR, LINE, "!", *100)
                     call Extract_Numb(STR,1,n,STRNUMB)
                     CtrlParam%Energy(I) = DRSTR(STRNUMB(1))
                  end if

                 !*** To get type  of incident orientation
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 call Extract_Numb(STR,1,n,STRNUMB)
                 CtrlParam%ANTYPE(I) = ISTR(STRNUMB(1))

                 if(CtrlParam%ANTYPE(I).EQ. CP_ORIENT_STYLE_MONODIR) then
                    !--- mono incident angle
                    call GetInputStrLine(hFile,STR, LINE, "!", *100)
                    call Extract_Numb(STR,1,n,STRNUMB)
                    CtrlParam%Orient(I) = DRSTR(STRNUMB(1))
                 else if(CtrlParam%ANTYPE(I).EQ. CP_ORIENT_STYLE_LATTDIR) then
                    !--- incident in given lattice direction
                    call GetInputStrLine(hFile,STR, LINE, "!", *100)
                    call Extract_Numb(STR,3,n,STRNUMB)
                    if(n.lt.1) goto 100
                    if(n.ge.1) CtrlParam%MILLER(I,1) = DRSTR(STRNUMB(1))
                    if(n.ge.2) CtrlParam%MILLER(I,2) = DRSTR(STRNUMB(2))
                    if(n.ge.3) CtrlParam%MILLER(I,3) = DRSTR(STRNUMB(3))
                    if(ALL(CtrlParam%MILLER(I,1:3) .eq. 0) ) then
                       print *, "Error in input Miller index"
                    end if
                 end if
              END DO

              !--- renormalize the percentage
              CtrlParam%Percent = CtrlParam%Percent/sum(CtrlParam%Percent)
         return
  !---------------------------------------------------------------
  100    print *, "Error in reading deposition controling parameters."
         print *, "The process to be stopped."
         stop
  end subroutine Load_Parameter_DepositionCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Print_Parameter_DepositionCtrlParam(hFile, CtrlParam)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !                CtrlParam
  !
  !     OUTPUT
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(DepositionCtrlParam)::CtrlParam

     !--- local variables
         character*256::STR,STRTMP(1)=""
      character*32::STRNUMB(30)
      integer::I, J, N
     !----

        !**** HEADER
          write(hFile,*) "!************ DEPOSITION CONTROL PARAMETERS **********"


        !*** the type of deposition
            if(IAND(CtrlParam%DEPSTYLE, 1).eq.0) then
               write(hFile,FMT="(' !    External deposition with: ', 3(I7,1x))")
            else
               write(hFile,FMT="(' !    Internal deposition with: ', 3(I7,1x))")
                    !case(CP_DEP_STYPE_INT_SIA)
                    !    write(hFile,FMT="(' !    SIA cascade with: ', 3(I7,1x))")
            end if


        !*** the incident particles
              write(hFile,FMT="(' !    Number of atoms.............................: ', 3(I7,1x))")   CtrlParam%NPRT

        !*** number types of the incident particles
              write(hFile,FMT="(' !    Number of types of incident atoms...........: ', 3(I7,1x))")   CtrlParam%nTYPE


        !***  types of the incident particles
              write(hFile,FMT="(' !    Types of incident atoms.....................: ', 8(I4,1x))")    CtrlParam%aTYPE(1:CtrlParam%nTYPE)

        !*** percentage
              write(hFile,FMT="(' !    Percentage of each kind of atoms........... : ', 8(F8.3,1x))")   CtrlParam%Percent(1:CtrlParam%nTYPE)


              DO I=1,  CtrlParam%nTYPE
                 write(hFile,FMT="(' !    For atom of type ',(I2),':')")  I

        !*** incident energies
                 write(hFile,FMT="(' !    Incident energies(KeV): ', 8(F10.4,1x))")   CtrlParam%Energy(I)

        !*** incident orientation
                 select case(CtrlParam%ANTYPE(I))
                        case(CP_ORIENT_STYLE_MONODIR)
                        write(hFile,FMT="(' !    In mono-direction (deg): ', 8(F10.4,1x))")   CtrlParam%Orient(I)

                        case(CP_ORIENT_STYLE_RANDOM4PI)
                        write(hFile,FMT="(' !    In random direction in 4PI ')")

                        case(CP_ORIENT_STYLE_LATTDIR)
                        write(hFile,FMT="(' !    In lattice direction with Miller indices: ',3(I3,1x))") int(CtrlParam%MILLER(I,1:3)+0.00001)
                 end select

              END DO
              write(hFile,*)


         return
  end subroutine Print_Parameter_DepositionCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine GenerateParticle_DepositionCtrlParam(CtrlParam, ITYP, EN, DIR)
  !***  PURPOSE:   to generate a particle according to the control parameter
  !
  !     INPUT:     CtrlParam
  !
  !     OUTPUT     ITYP, the type of the particle
  !                EN,   the kinetic energy of the particle
  !                DIR,  the moving direction of the particle
  use RAND32_MODULE
  implicit none
     !--- dummy varioables
     type(DepositionCtrlParam)::CtrlParam
     integer::ITYP
     real(KINDDF)::EN, DIR(3)

     !--- local variables
         integer::I,IND
         real(KINDDF)::R, THETA, FAI

            !--- determine the type of the particle
              R = DRAND32()
              if(0.le.R .and.R.lt.CtrlParam%Percent(1)) then
                 IND = 1
              else
                 DO I=2, CtrlParam%nType
                    if( CtrlParam%Percent(I-1).le.R .and. R.lt.CtrlParam%Percent(I)) then
                       IND = I
                       exit
                    end if
                 END DO
             end if
             ITYP = CtrlParam%aTYPE(IND)
            !--- determine the energy, now only mono-energy supported
            select case(CtrlParam%ENTYPE(IND))
                   case (0)
                        EN= CtrlParam%Energy(IND)
                   case default
                        EN= CtrlParam%Energy(IND)
            end select

            !--- determine the initial moving direction
            select case(CtrlParam%ANTYPE(IND))
                   !--- mono direction
                   case (CP_ORIENT_STYLE_MONODIR)
                        THETA = CtrlParam%Orient(IND)*CP_PI/180.D0
                        FAI   = DRAND32()*CP_TWOPI
                        DIR(3)= -DCOS(THETA)
                        DIR(2)= DSIN(THETA)*DCOS(FAI)
                        DIR(1)= DSIN(THETA)*DSIN(FAI)

                   case default
                        THETA = CtrlParam%Orient(IND)*CP_PI/180.D0
                        FAI   = DRAND32()*CP_TWOPI
                        DIR(3)= -DCOS(THETA)
                        DIR(2)= DSIN(THETA)*DCOS(FAI)
                        DIR(1)= DSIN(THETA)*DSIN(FAI)
                   !---- randomly in 4PI solid angle, used in case of inner irradiation
                   case (CP_ORIENT_STYLE_RANDOM4PI)
                        THETA = DRAND32()*CP_PI
                        FAI   = DRAND32()*CP_TWOPI
                        DIR(3)= DCOS(THETA)
                        DIR(2)= DSIN(THETA)*DCOS(FAI)
                        DIR(1)= DSIN(THETA)*DSIN(FAI)
                   !---- in a lattice direction defined by Miller index
                   case (CP_ORIENT_STYLE_LATTDIR)
                        R = dsqrt(sum(dble(CtrlParam%MILLER(IND,1:3))**2))
                        DIR(1)= CtrlParam%MILLER(IND,1)/R
                        DIR(2)= CtrlParam%MILLER(IND,2)/R
                        DIR(3)= CtrlParam%MILLER(IND,3)/R
            end select


         return
  end subroutine GenerateParticle_DepositionCtrlParam
  !****************************************************************************


  end module DEPOSITION_TypeDef_CtrlParam_2010
