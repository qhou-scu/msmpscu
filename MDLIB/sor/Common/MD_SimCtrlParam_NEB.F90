  module MD_SimMDCtrl_NEB
  !***  DESCRIPTION: this module is to load control parameters for NEB
  !
  !                  ______________________________________________________
  !     HISORY:      2016-11-16, seperated the from MD_TypeDef_SimMDCtrl.F90
  !
  use MD_CONSTANTS
  use MiniUtilities
  use MD_TYPEDEF_SimMDCtrl
  implicit none
      private::Load_NEBCtl_SimMDCtrl
  contains
  !*********************************************************************

  !****************************************************************************
  subroutine Load_NEBCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
  !***  PURPOSE:   to load the control parameters NEB calaculations
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  implicit none
     !--- dummy varioables
     integer,         intent(in):: hFile
     type(SimMDCtrl)            :: CtrlParam
     character*(*)              :: STR
     integer                    :: LINE
     !--- local variables
      character*32::STRNUMB(10), STRTMP(1), KEYWORD
      integer::I, N
     !----

          !**** to start load the controal parameters
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                     case( "&ENDSUBCTL")
                          exit

                     case default
                          write(*,*)"MDPSCU warning: unknown keyword in &NEBSUBCTL control ", KEYWORD(1:LEN_TRIM(KEYWORD))
                          write(*,fmt="('               check control file at line:', BZI6)") LINE
                          call ONWARNING(gm_OnWarning)

                     case ("&METH")
                           !$$*** To get number of images of each path
                           call Extract_Substr(STR,1,n,STRTMP)
                           CtrlParam%NEB_Meth = STRTMP(1)(1:len(CtrlParam%NEB_Meth))
                           call UpCase(CtrlParam%NEB_Meth)
                           select case(CtrlParam%NEB_Meth(1:len_trim(CtrlParam%NEB_Meth)))
                                  case("DNEB")
                                  case("NEB")
                                  case("EB")
                                  case default
                                      write(*,fmt="(A)") " MDPSCU Error: unknown NEB method: "//CtrlParam%NEB_Meth(1:len_trim(CtrlParam%NEB_Meth))
                                      write(*,fmt="(A)") "               process to be stopped"
                                      stop
                            end select

                     case ("&NUMIMGS")
                           !$$*** To get number of images of each path
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_NumImg = ISTR(STRNUMB(1))
                           if(CtrlParam%NEB_NumImg .lt. C_ITHR) then
                              write(*,fmt="(A)")       " MDPSCU Error: number of images is smaller than 3 "
                              write(*,fmt="(A, BZI6)") "               check control file at line:", LINE
                              write(*,fmt="(A)")       "               process to be stopped"
                              stop
                           end if
                     case ("&KVAL", "&STRENGTH")
                           !$$*** To get spring strength
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_KValue = DRSTR(STRNUMB(1))

                     case ("&CLIMB")
                           !$$*** To get flag for climbing search
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_Climb = ISTR(STRNUMB(1))

                     case ("&FINDNODE")
                           !$$*** To get flag for climbing search
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_FindMin = ISTR(STRNUMB(1))

                     case ("&SERIALCFG")
                           !$$*** To get flag for climbing search
                           call Extract_Numb(STR,3,n,STRNUMB)
                           if(N .LT. 1) then
                              write(*,fmt="(' MDPSCU Warning: parameters for serial NEB calculations is missed')")
                              write(*,fmt="('                 check control file at line:', BZI6)") LINE
                              write(*,fmt="(' Usage:  &SERIALCFG startcfg, endcfg ')")
                              call ONWARNING(gm_OnWarning)
                              CtrlParam%NEB_StartCfg = -1
                              CtrlParam%NEB_EndCfg   = -1
                           else
                              CtrlParam%NEB_StartCfg = ISTR(STRNUMB(1))
                              if(N .ge. 2) CtrlParam%NEB_EndCfg   = ISTR(STRNUMB(2))
                              if(N .ge. 3) CtrlParam%NEB_CfgStep  = ISTR(STRNUMB(3))
                              if(CtrlParam%NEB_StartCfg  .ge. 0) then
                                 if(CtrlParam%NEB_EndCfg .lt. 0) CtrlParam%NEB_EndCfg =  CtrlParam%NEB_StartCfg + CtrlParam%NEB_CfgStep
                              end if
                           end if

                     case ("&OUTPUTINI")
                           !$$*** To get flag for climbing search
                           call Extract_Numb(STR,1,n,STRNUMB)
                           CtrlParam%NEB_OutIniCfg = ISTR(STRNUMB(1))

                     case("&MINIMIZER", "&QUICKDUMP", "&QUICKDAMP", "&QUENCH")
                           !*** To get if damping is required
                           call Extract_Substr(STR,1,N,STRTMP)
                           if(n .ge. 1) then
                              call UpCase(STRTMP(1))
                              if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "LBFGS") then
                                 CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_LBFGS
                              else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "CG") then
                                 CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_CG
                              else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "CG-LS") then
                                 CtrlParam%DAMPSCHEME = ior(CP_DAMPSCHEME_CG, CP_DAMPSCHEME_LSEARCH)                                 
                              else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "ST") then
                                 CtrlParam%DAMPSCHEME = CP_DAMPSCHEME_ST
                              else if(STRTMP(1)(1:len_trim(STRTMP(1))) .eq. "ST-LS") then
                                 CtrlParam%DAMPSCHEME = ior(CP_DAMPSCHEME_ST, CP_DAMPSCHEME_LSEARCH)                                     
                              else
                                 write(*,fmt="(A)")      " MDPSCU Error: the minizing method "//STRTMP(1)(1:len_trim(STRTMP(1)))//" is unknown for NEB "
                                 write(*,fmt="(A, BZI6)")'               check control file at line:', LINE
                                 write(*,fmt="(A)")      '               supported methods are: "LBFGS", or "DYNAMICS" '
                                 write(*,fmt="(A)")      ' Process to be stopped'
                                 stop
                              end if
                           end if

                           call Extract_Numb(STR,1,n,STRNUMB)
                           if(n .ge. 1) then
                              CtrlParam%DAMPTIME1 = ISTR(STRNUMB(1))
                           end if

                     case("&LBFGSSUBCTL")
                           call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                     case("&STEEPESTSUBCTL", "&CGSUBCTL")
                           call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                     case ("&EVENTDETECTSUBCTL")
                           call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                end select
             end do

          return
   !----------------------------------------------------------------------------
  100    write(*,*)"MDPSCU Error in reading NEB control parameters."
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDSUBCTL exist to end the data subsection"
         stop
   end subroutine Load_NEBCtl_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Load_Parameter_SimMDCtrl(hFile, CtrlParam, SimBox, LINE0)
  !***  PURPOSE:   to load the control parameters from a file
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT     CtrlParam
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer,                intent(in):: hFile
     type(SimMDCtrl),target            :: CtrlParam
     type(SimMDBox)                    :: SimBox
     integer,optional,       intent(in):: LINE0
     !--- local variables
      character*256::STR
      character*32::KEYWORD
      integer::LINE
     !----

          if(present(LINE0)) then
             LINE = LINE0
          else
             LINE = 0
          end if

          call GetInputStrLine(hFile,STR, LINE, "!", *100)
          STR = adjustl(STR)
          call GetKeyWord("&", STR, KEYWORD)
          call UpCase(KEYWORD)
          !**** to start load the controal parameters
             do while(.TRUE.)
                call GetInputStrLine(hFile,STR, LINE, "!", *100)
                STR = adjustl(STR)
                call GetKeyWord("&", STR, KEYWORD)
                call UpCase(KEYWORD)
                select case (KEYWORD(1:LEN_TRIM(KEYWORD)) )
                       case ("&ENDCTLF")
                          exit

                       case default
                            write(*,*)"MDPSCU warning: unknown keyword ", KEYWORD(1:LEN_TRIM(KEYWORD))
                            write(*,fmt="('               check control file at line:', BZI6)") LINE
                            call ONWARNING(gm_OnWarning)

                       case ("&COMMSUBCTL")
                            call Load_CommonParameter_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&POTENSUBCTL","&POTSUBCTL")
                             call Load_ForceCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, SimBox%NGROUP, LINE)

                       case ("&ANALYSUBCTL")
                             call Load_AnalyCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&BOUNDSUBCTL")
                             call Load_BoundCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&NEIGHBSUBCTL")
                             call Load_NeighbCutoffCtl_SimMDCtrl(hFile, CtrlParam, STR, SimBox%NGROUP, LINE)

                       case ("&ACTIVEREGSUBCTL")
                             call Load_ActiveRegionCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&NEBSUBCTL")
                             call Load_NEBCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&LBFGSSUBCTL")
                             call Load_LBFGSCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)
                              
                       case ("&STEEPESTSUBCTL", "&CGSUBCTL")
                             call Load_STEEPESTCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                       case ("&EVENTDETECTSUBCTL")
                             call Load_EventDetectCtl_SimMDCtrl(hFile, CtrlParam, STR, LINE)

                end select
              end do

              !$$*** to check the consistent of the control parameters
              if(CtrlParam%NEB_StartCfg .ge. 0) then
                 if(CtrlParam%MultiBox .gt. 1) then
                    write(*,fmt="('MDPSCU error: for serial NEB calculations, multi-box cannot be larger than 1 ', I3)")
                    write(*,fmt="('              process to be stopped', I3)")
                    stop
                 end if
              end if

              !$$***
               if(CtrlParam%AR_UPTABMI .le. 0) then
                  CtrlParam%AR_UPTABMI = CtrlParam%NB_UPTABMI
                  CtrlParam%AR_UPTABMX = CtrlParam%NB_UPTABMI
                  CtrlParam%AR_DBITAB  = CtrlParam%NB_DBITAB
                  CtrlParam%AR_UPTAB   = CtrlParam%AR_UPTABMI
               end if

        return
  !---------------------------------------------------------------
  100    write(*,fmt="(A, I5)")" MDPSCU Error in reading NEB control parameters at line ", LINE
         write(*,*)"The process to be stopped."
         write(*,*)"Please check if keyword &ENDCTLF exist to end the data section"
         stop
  end subroutine Load_Parameter_SimMDCtrl
  !****************************************************************************

  !****************************************************************************
  subroutine Print_Parameter_SimMDCtrl(hFile, CtrlParam, SimBox)
  !***  PURPOSE:   to print out control parameters for NEB applications
  !
  !     INPUT:     hFile,  I/O unit number
  !
  !     OUTPUT
  use MD_TYPEDEF_SimMDBox
  implicit none
     !--- dummy varioables
     integer,intent(in)    ::hFile
     type(SimMDCtrl),target::CtrlParam
     type(SimMDBox)        ::SimBox

     !--- local variables
     integer::I, J, CENTPART(mp_MXGROUP)
     !----

        !$$**** HEADER
          write(hFile,fmt="(' !************ SIMULATION CONTROL PARAMETERS *******************')")

              if(CtrlParam%RESTART.le.0) then
                 write(hFile,FMT="(' !   A new simulation to be performed  ')")
              else
                 write(hFile,FMT="(' !   The simulation to be restarted from last stop point', 3(1PE11.4,1x))")
              end if
              write(hFile,fmt="(' !')")

              !*************************
              if(CtrlParam%NEB_StartCfg .ge. 0) then
                 write(hFile,FMT="(' !    NEB to be performed for a configure serial....: ', (BZI7,1x), ' to ', (BZI7,1x), ' with step',(BZI7,1x) )") &
                                  CtrlParam%NEB_StartCfg, CtrlParam%NEB_EndCfg, CtrlParam%NEB_CfgStep
              else
                 write(hFile,FMT="(' !    NEB to be performed for a configure serial....: ', 3(BZI7,1x))")    CtrlParam%StartCfg
                 write(hFile,FMT="(' !    Number of boxes in one test ..................: ', 3(BZI7,1x))")    CtrlParam%MULTIBOX
              end if

              write(hFile,FMT="(' !    Total number of boxes to be considered........: ', 3(BZI7,1x))")    CtrlParam%TotalBox
              write(hFile,FMT="(' !    Number of boxes in one test ..................: ', 3(BZI7,1x))")    CtrlParam%MULTIBOX

              write(hFile,FMT="(' !    Number of images on one path .................: ', 3(BZI7,1x))")    CtrlParam%NEB_NUMIMG
              write(hFile,FMT="(' !    KVALE.........................................: ',  (1PE11.4,1x))") CtrlParam%NEB_KVALUE
              if(CtrlParam%NEB_Climb .gt. 0) then
              write(hFile,FMT="(' !    If climbing seraching is enabled..............: YES ',  (1PE11.4,1x))")
              else
              write(hFile,FMT="(' !    If climbing seraching is enabled..............: NO ',  (1PE11.4,1x))")
              end if
              write(hFile,FMT="(' !    Number of iteration of finding minimum nodes..: ',  I5)")              CtrlParam%NEB_FindMin
              write(hFile,FMT="(' !    Displacement threshold (LU) for strcuture change....: ', 20(F8.2,1x))") CtrlParam%STRCUT_DRTol

              write(hFile,FMT="(' !    Minimization parameter in LBFGS....................: ', BZI8)")
              write(hFile,FMT="(' !       the convergence criteria of find minia..........: PGTOL   ', 3(1PE11.4,1x))") CtrlParam%LBFGS_PGtol
              write(hFile,FMT="(' !                                                         FACTR   ', 3(1PE11.4,1x))") CtrlParam%LBFGS_Factr
              write(hFile,FMT="(' !                                                         MSAVE   ', BZI8)")          CtrlParam%LBFGS_MSave
              write(hFile,FMT="(' !       when having spring force .......................: S_PGTOL ', 3(1PE11.4,1x))") CtrlParam%NEB_LBFGS_PGtol
              write(hFile,FMT="(' !                                                         S_FACTR ', 3(1PE11.4,1x))") CtrlParam%NEB_LBFGS_Factr
              write(hFile,FMT="(' !                                                         S_MSAVE ', BZI8)")          CtrlParam%NEB_LBFGS_MSave


              if(CtrlParam%NEB_OutIniCfg .gt. 0) then
              write(hFile,FMT="(' !    Output initial configurations.................: YES',  I5)")
              else
              write(hFile,FMT="(' !    Output initial configurations.................: NO',  I5)")
              end if
              write(hFile,fmt="(' !')")

              !*************************
              write(hFile,FMT="(' !    For force table calculation...................: ', 20(BNI8,1x))")
              write(hFile,FMT="(' !       Number of points in force table............: ', 20(BNI8,1x))") CtrlParam%NUMFTABR
              write(hFile,FMT="(' !       Number of points in embed-fun table........: ', 20(BNI8,1x))") CtrlParam%NUMFTABE
              !write(hFile,FMT="(' !    RHO-scal in embedment function table..........: ',  (1PE11.4,1x))")CtrlParam%RHOSCAL

              write(hFile,FMT="(' !       Cutoff range of interactions (in LU).......: ', 20(F8.2,1x))")
              do I=1, SimBox%NGROUP
                 do J=I, SimBox%NGROUP
                    write(hFile,FMT="(' !           between group',I2,' and group',I2,' ...........:', 20(F8.2,1x))")I, J, CtrlParam%RU(I,J)
                 end do
              end do
              write(hFile,fmt="(' !')")

             !--- the input for neighbor-list calculation
              write(hFile,FMT="(' !    If periodic conditions used (X,Y,Z).........: ', 3(I3,1x))")   CtrlParam%IFPD(1:3)
              write(hFile,FMT="(' !    Max permitted neighbors.....................: ', 3(BNI4,1x))") CtrlParam%NB_MXNBS
              write(hFile,FMT="(' !    Cutoff range of neighbor region(in RCUT)....: ', 20(F8.2,1x))")
              do I=1, SimBox%NGROUP
                 do J=I, SimBox%NGROUP
                    write(hFile,FMT="(' !        between group',I2,' and group',I2,'.............:', 20(F8.2,1x))")I, J, CtrlParam%NB_RM(I,J)
                 end do
              end do

             !--- the input for activation region control
              CENTPART = 0
              I = 0
              do J=1, mp_MXGROUP
                 if(CtrlParam%AR_CENTPART(J) .gt. 0) then
                    I = I + 1
                    CENTPART(I) = J
                 end if
             end do

             if(iand(CtrlParam%AR_METHOD,CP_ENABLE_AR) .ne. CP_ENABLE_AR) then
                 write(hFile,FMT="(' !      Activation region calculation ..............: disable', 3(BNI4,1x))")
             else
                 write(hFile,FMT="(' !      Activation region calculation ..............: enable', 3(BNI4,1x))")
                 write(hFile,FMT="(' !      Min step number to update activation region.....: ', 3(BNI4,1x))") CtrlParam%AR_UPTABMI
                 write(hFile,FMT="(' !      Max step number to update activation region.....: ', 3(BNI4,1x))") CtrlParam%AR_UPTABMX
                 write(hFile,FMT="(' !      Step number to change updating intervale....: ', 3(BNI6,1x))")     CtrlParam%AR_DBITAB

                 if(iand(CtrlParam%AR_METHOD, CP_USERDEF_AR) .eq. CP_USERDEF_AR ) then
                    write(hFile,FMT="(' !      User supplied procedure to be used for AR calculation', 3(BNI6,1x))")
                 else
                    if(iand(CtrlParam%AR_METHOD, CP_CENTPART_AR) .eq. CP_CENTPART_AR .and.  &
                       iand(CtrlParam%AR_METHOD, CP_EKIN_AR) .eq. CP_EKIN_AR) then
                       write(hFile,FMT="(' !      AR build by seeds of atom type..............:', 20(BNI6,1x))")      &
                                         (CENTPART(J), J=1, count(CENTPART.gt.0))
                       write(hFile,FMT="(' !           and by kinetic energy..................:', F8.3,1x, '(eV)'))") &
                                            CtrlParam%AR_EKIN
                    else if( iand(CtrlParam%AR_METHOD, CP_CENTPART_AR).eq.CP_CENTPART_AR) then
                       write(hFile,FMT="(' !      AR build by seeds of atom type..............:', 20(BNI6,1x))")      &
                                         (CENTPART(J), J=1, count(CENTPART.gt.0))
                    else
                       write(hFile,FMT="(' !      AR build by seeds of atoms of energy above..:', F8.3,1x, '(eV)'))") &
                                          CtrlParam%AR_EKIN
                   end if

                   if(ibits(CtrlParam%AR_METHOD, CP_BYNBSETBIT_AR,1) .gt. 0) then
                       write(hFile,FMT="(' !      with an extend..............................: ', 1(BNI6,1x), &
                                         ' by neighborhood')") CtrlParam%AR_EXTEND
                   else
                       write(hFile,FMT="(' !      with an extend..............................: ', 1(BNI6,1x), &
                                         ' by linked-cells')") CtrlParam%AR_EXTEND
                   end if

                   if(iand(CtrlParam%AR_METHOD,CP_KEEP_AR) .eq. CP_KEEP_AR) then
                       write(hFile,FMT="(' !      keeping activation statu of particles after disabling AR')")
                   else 
                       write(hFile,FMT="(' !      not keeping statu of particles after disabling AR')")
                   end if    
                 end if
             end if
             write(hFile,fmt="(' !')")

              !*************************
              if(CtrlParam%SEED(1) .gt. 0) then
                 write(hFile,FMT="(' !    Random number seed ...........................: ',4(I8,1x))") CtrlParam%SEED(1)
              else
                 write(hFile,FMT="(' !    Random number seed determined automatically  ',4(I8,1x))")
              end if
              write(hFile,fmt="(' !')")

         return
  end subroutine Print_Parameter_SimMDCtrl
  !****************************************************************************

  end module MD_SimMDCtrl_NEB
