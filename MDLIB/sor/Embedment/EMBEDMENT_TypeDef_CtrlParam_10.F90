  module EMBEDMENT_TypeDef_CtrlParam_2010
  !***  DESCRIPTION: this module is to define the data type for Depostion Control Parameter.
  !
  !                  ______________________________________________________
  !                  HOU Qing, Nov, 2010
  use MD_CONSTANTS
  use MiniUtilities
  implicit none

   !--------------------------------------------------
   !--- The type describes the data controlling the embeding
   !    simulation
      integer,parameter,private::m_MaxType = 4         !max number of types for embedded clusters
      integer, parameter::EM_STYLE_RANDOM   = 0
      integer, parameter::EM_STYLE_CLUSTER  = 1
      integer, parameter::EM_STYLE_REPLACE  = 2
      integer, parameter::EM_STYLE_VACANCY  = 3
      type EMBEDCtrlParam
           !--- simulation parameters
           !integer::NPRT   = 0                        !total number of embeded particles to be simulated
                                                       !NOTE: the actual number of particle is (NPRT/MULTIBOX)*MULTIBOX
                                                       !      MULTIBOX is given in strcuture SimMDCtrl
                                                       !      refer the definition of SimMDCtrl

           integer::nType=0                            !number of types of embedded clusters
           integer::aType(m_MaxType)=0                 !type of atoms in embedded clusters
           integer::PNPRT(m_MaxType)=0                 !number of atoms for each type of atom
           integer::EmStyle(m_MaxType)=EM_STYLE_RANDOM

           real(KINDDF)::EMMinSep(m_MaxType)=1.D64     !min seperation at placing atoms
           integer::cType(m_MaxType)=0                 !type of cluster the embedded atom to be appended to, if EmStyle=EM_STYLE_CLUSTER
                                                       !type of atom to be replaced with by the embedeed atom, if EmStyle=EM_STYLE_REPLACE

           real(KINDDF)::EMINR(m_MaxType)=1.D64        !min radius the embedded atom from the center of appended cluster,used only when EmStyle=EM_STYLE_CLUSTER
           real(KINDDF)::EMOUTR(m_MaxType)=1.D64       !min radius the embedded atom from the center of appended cluster,used only when EmStyle=EM_STYLE_CLUSTER


      end type EMBEDCtrlParam
  !--------------------------------------------------



  contains
  !*********************************************************************

  !****************************************************************************
  subroutine Load_Parameter_EMBEDCtrlParam(hFile, CtrlParam)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(EMBEDCtrlParam)::CtrlParam

     !--- local variables
      character*256::STR,STRTMP(2)=""
      character*32::STRNUMB(30)
      integer::I, J, N, LINE
     !----
              LINE = 0
         !*** To get number of types of embedded atoms
              call GetInputStrLine(hFile,STR, LINE, "!", *100)
              call Extract_Numb(STR,1,n,STRNUMB)
              CtrlParam%nType = ISTR(STRNUMB(1))

              DO I=1, CtrlParam%nType

                !*** To get the type of embedded atoms
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 call Extract_Numb(STR,1,n,STRNUMB)
                 CtrlParam%aType(I) = ISTR(STRNUMB(1))

                !*** To get number of atoms to be embedded
                 call GetInputStrLine(hFile,STR, LINE, "!", *100)
                 call Extract_Numb(STR,1,n,STRNUMB)
                 CtrlParam%PNPRT(I) = ISTR(STRNUMB(1))

                !*** To get the style of placing the atoms
                 call GetInputStrLine(hFile,STR, LINE, "!",  *100)
                 call Extract_Numb(STR,1,n,STRNUMB)
                 CtrlParam%EmStyle(I) = ISTR(STRNUMB(1))


                 if(CtrlParam%EmStyle(I) .eq. EM_STYLE_RANDOM .or. &
                    CtrlParam%EmStyle(I) .eq. EM_STYLE_CLUSTER ) then
                    !*** To get the min seperation at placing the atoms
                    call GetInputStrLine(hFile,STR, LINE, "!", *100)
                    call Extract_Numb(STR,1,n,STRNUMB)
                    CtrlParam%EMMinSep(I) = DRSTR(STRNUMB(1))

                   if( CtrlParam%EmStyle(I) .eq. EM_STYLE_CLUSTER ) then
                      !*** To get the type of cluster to be appended
                      call GetInputStrLine(hFile,STR, LINE, "!", *100)
                      call Extract_Numb(STR,1,n,STRNUMB)
                      CtrlParam%cType(I) = ISTR(STRNUMB(1))

                      !*** To get radiu region the atom from the center of the cluster
                      call GetInputStrLine(hFile,STR, LINE, "!", *100)
                      call Extract_Numb(STR,2,n,STRNUMB)
                      CtrlParam%EMINR(I) = DRSTR(STRNUMB(1))
                      if(n .lt. 2) then
                         CtrlParam%EMOUTR(I) = CtrlParam%EMINR(I)
                      else
                        CtrlParam%EMOUTR(I) = DRSTR(STRNUMB(2))
                      end if
                   end if
                 !---- for replacement case
                 else if(CtrlParam%EmStyle(I) .eq. EM_STYLE_REPLACE) then
                      !*** To get the type of atom to be replaced
                      call GetInputStrLine(hFile,STR, LINE, "!", *100)
                      call Extract_Numb(STR,1,n,STRNUMB)
                      CtrlParam%cType(I) = ISTR(STRNUMB(1))

                end if
              END DO

         return
  !---------------------------------------------------------------
  100    print *, "MDPSCU Error in reading embeding control parameters at line:", LINE
         print *, "The process to be stopped."
         stop
  end subroutine Load_Parameter_EMBEDCtrlParam
  !****************************************************************************

  !****************************************************************************
  subroutine Print_Parameter_EMBEDCtrlParam(hFile, CtrlParam)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !                CtrlParam
  !
  !     OUTPUT
  implicit none
     !--- dummy varioables
     integer, intent(in)::hFile
     type(EMBEDCtrlParam)::CtrlParam

     !--- local variables
      integer::I
     !----

        !**** HEADER
          write(hFile,*) "!************ EMBEDING CONTROL PARAMETERS **********"


        !*** number types of the incident particles
              write(hFile,FMT="(' !    Number of types of embedded atoms...........: ', 3(I7,1x))")   CtrlParam%nTYPE

              DO I=1, CtrlParam%nTYPE
              write(hFile,FMT="(' !    Type of ',I2,'th embedded atoms.................: ', 8(I4,1x))") I, CtrlParam%aTYPE(I)

             !*** percentage
              write(hFile,FMT="(' !    Number of this kind of atoms............... : ', 8(I4,1x))")   CtrlParam%PNPRT(I)

             !*** embedded style
              write(hFile,FMT="(' !    Embedding style of this kind of atoms...... : ', 8(I4,1x))")  CtrlParam%EMSTYLE(I)

             !*** min seperation
              write(hFile,FMT="(' !    Minimum atom-atom seperation at embedding.. : ', 8(F8.2,1x))")   CtrlParam%EMMinSep(I)

             if( CtrlParam%EMSTYLE(I) .eq. EM_STYLE_CLUSTER) then
             !*** cluster to be appended to
              write(hFile,FMT="(' !    Type of appended cluster.................... : ', 8(I4,1x))")     CtrlParam%cType(I)

             !*** cluster to be appended to
              write(hFile,FMT="(' !    Distance from the cluster center............ : ', 8(F8.2,1x))")   CtrlParam%EMINR(I),CtrlParam%EMINR(I)
              end if
              write(hFile,*)

              END DO

         return
  end subroutine Print_Parameter_EMBEDCtrlParam
  !****************************************************************************

  !**************************************************************************
  subroutine Initialize_EMBEDCtrlParam(fname,CtrlParam)
  !***  PURPOSE: to intilaize the embed module
  !**** DESCRIPTION: to
  !****
  !     INPUT: fname      , the simulation box
  !
  !     OUTPUT:CtrlParam,   the control parameters for embeding
  !
  implicit none
      !----   DUMMY Variables
       character*(*)::fname
       type(EMBEDCtrlParam)::CtrlParam

       integer::hFile



            !*** to load controling parameters
            write(*,*) "!**** Load depostion parameters from:"
            write(*,*) "!**** ",fname(1:len_trim(fname))
            call AvailableIOUnit(hFile)
            open(hFile, file = fname, status='old')
            call  Load_Parameter_EMBEDCtrlParam(hFile, CtrlParam)
            close(hFile)

        return
   end subroutine Initialize_EMBEDCtrlParam
  !****************************************************************************************


  end module EMBEDMENT_TypeDef_CtrlParam_2010
