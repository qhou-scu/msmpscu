  module DEPOSITION_BAC_TEST

  !**** DESCRIPTION: This is a module modified from module DEPOSITION_COMMON_2010_GPU, for the specific
  !                  purpose of test BAC in cascade stimulation of an energetic atom on/in a substarte.
  !
  !                  The force calculation is performed in BAC_TEST_TB_ForceTable_GPU, in which
  !                  the force on the energetic atom is calculated for only the neareast target atoms.
  !
  !
  !                  SEE ALSO______________________________________________________________________________
  !                       DEPOSITION_COMMOM_10_GPU.f90.f90
  !                       BAC_TEST_TB_ForceTable_GPU.f90
  !                  ______________________________________________________________________________________
  !                  HOU Qing, Oct, 2011

  !***  The modules included ******************************************
  use MD_Globle_Variables_2010
  use MD_Globle_Variables_2010_GPU
  use MD_NeighborsList2C_2010_GPU
  use BAC_TEST_TB_ForceTable_GPU
  !use MD_FINNIS_EM_TB_Force_Table_2010_GPU
  use MD_SWOPE_Scheme_2010_GPU
  use MD_EP_Coupling_2010_GPU
  use MD_BP_Coupling_2010_GPU

  !*** the application specific module ********************************
   use DEPOSITION_TypeDef_CtrlParam_2010
   use DEPOSITION_COMMON_2010
  implicit none

  !_______________________________________________________________________________________
  !   The simulation data block
  !
      type(SimMDBox)           ::m_SimBox
      type(SimulationCtrlParam)::m_CtrlParam
      type(MDForceTable)       ::m_ForceTable
      type(DepositionCtrlParam)::m_DepCtrlParam
  !_______________________________________________________________________________________

  contains

  !****************************************************************************************
  subroutine SubProcess(NMPI,processid, FORCETABLE, PRERECORD, RECORDPROC,AFTRECORD)
  !--- INPUT: NMPI, number of prcocesses
  !           processid, the process ID
  !           FORCETABLE, the subroutine provided by user for force register
  !           PRERECORD,  the optional subroutine provided by user to pre-process before the main-loop for time
  !           RECORDPROC, the optional subroutine provided by user to process data at time when the configure to be output
  !           AFTRECORD,  the optional subroutine provided by user after the time-loop end
  !    NOTE: the three dummy subroutine must strictly follow the interface given below
  use RAND32SEEDLIB_MODULE
  use RAND32_MODULE
  implicit none
  !--- DUMMY VAIRBALES
  integer::NMPI, processid
  !--- dunmmy functions
      OPTIONAL::PRERECORD, RECORDPROC,AFTRECORD
      external::FORCETABLE,PRERECORD, RECORDPROC,AFTRECORD
  !--- interface to the external routines -------------------
    interface
     subroutine FORCETABLE(SimBox, CtrlParam, FTable, printout)
     use MD_CONSTANTS_2010
     use MD_TYPEDEF_SIMULAIONBOX_2010
     use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
     use MD_TYPEDEF_FORCETABLE_2010
     implicit none
      !--- dummy vaiables
      type(SimMDBox),            intent(in)   ::SimBox
      type(SimulationCtrlParam), intent(in)   ::CtrlParam
      type(MDForceTable),        intent(inout)::FTable
      integer,                   optional     ::printout
     end subroutine FORCETABLE
   end interface

   interface
       SUBROUTINE PRERECORD(SimBox, CtrlParam,DEPCtrlParam,M)
       use MD_CONSTANTS_2010
       use MD_Globle_Variables_2010
       use MD_SIMULAIONBOX_ARRAY_2010
       use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
       use DEPOSITION_TypeDef_CtrlParam_2010
       implicit none
       type(SimMDBox)           ::SimBox
       type(SimulationCtrlParam)::CtrlParam
       type(DepositionCtrlParam)::DEPCtrlParam
       integer::M
       END SUBROUTINE PRERECORD
   end interface

   interface
       SUBROUTINE RECORDPROC(JP,ITIME, TIME, SimBox, CtrlParam,DEPCtrlParam,IPROC)
       use MD_CONSTANTS_2010
       use MD_Globle_Variables_2010
       use MD_SIMULAIONBOX_ARRAY_2010
       use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
       use DEPOSITION_TypeDef_CtrlParam_2010
       implicit none
       type(SimMDBox)           ::SimBox
       type(SimulationCtrlParam)::CtrlParam
       type(DepositionCtrlParam)::DEPCtrlParam
       integer::JP,ITIME
       real(KINDDF)::TIME
       integer::IPROC
       END SUBROUTINE RECORDPROC
   end interface

   interface
       SUBROUTINE AFTRECORD(SimBox, CtrlParam,DEPCtrlParam)
       use MD_CONSTANTS_2010
       use MD_Globle_Variables_2010
       use MD_SIMULAIONBOX_ARRAY_2010
       use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
       use DEPOSITION_TypeDef_CtrlParam_2010
       implicit none
       type(SimMDBox)           ::SimBox
       type(SimulationCtrlParam)::CtrlParam
       type(DepositionCtrlParam)::DEPCtrlParam
       END SUBROUTINE AFTRECORD
   end interface
  !--- END INTERFACE --------------------------------

  !--- LOCAL VARIABLES
  type(SimMDBox)::SimBox0
  integer::I, J, ISEED0, ISEED(2),EXIST, IDEV

  character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
  integer::DATE_TIME1(8),DATE_TIME2(8)
  real*4::C1,C2,TIME1, TIME2
  character*256::SFILE, GFile,GFileU, CFile, CFile0
  character*32:: ARGV
  real(KINDDF)::TIME


      !*** to start the device
          !*** to get the start id for device
          J = COMMAND_ARGUMENT_COUNT()
          if(J.gt.1) then
             !--- the first command line argument left for input filename
             call GET_COMMAND_ARGUMENT(2,ARGV)
             read(ARGV, *) IDEV
          else
             IDEV = 0
          end if

          J=cudaThreadExit();
          J = cudaSetDevice(processid+IDEV)

      !*** to read in control parameters from IniFile
          call Initialize_Globle_Variables(SimBox0, m_CtrlParam)
          call UnitConvert_Globle_Variables(SimBox0, m_CtrlParam)

      !*** to intialize the random number seed to call the random number geerator once
           ISEED0 = m_CtrlParam%SEED(1)
           call GetSeed_RAND32SEEDLIB(ISEED0, ISEED(1), ISEED(2))
           ISEED0=ISEED0+processid-1
           call GetSeed_RAND32SEEDLIB(ISEED0, ISEED(1), ISEED(2))
           call DRAND32_PUTSEED(ISEED)

      !*** to initialize the configuration
          call Initialize_SimMDBox(SimBox0)

      !*** to initialize the force-potential table
          call FORCETABLE(SimBox0, m_CtrlParam, m_ForceTable)

      !*** to initialize the electron phonon coupling module
          call Initialize_EPC_MODIFICATION_DEV(SimBox0, m_CtrlParam)

        !*** to initialize the differential scheme
            call Initialize_SWOPE_Scheme_DEV(SimBox0,m_CtrlParam)

      !*** to begin deposition simulation
          call Initialize_DEPOSITION(m_CtrlParam%f_others(1),m_DepCtrlParam)
          call Print_Parameter_DepositionCtrlParam(6, m_DepCtrlParam)
          call Print1_Globle_Variables( m_SimBox, m_CtrlParam)

        !**** if user provide preprocessor, we do it
            if(present(PRERECORD)) CALL PRERECORD(SimBox0,m_CtrlParam,m_DepCtrlParam,NMPI)

      !*** to recorde the start time
          CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
          C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000


      !*** TO BEGIN DEPOSITION LOOP
         DO J=1, (m_CtrlParam%TOTALBOX*m_CtrlParam%MULTIBOX)/NMPI !max(1,m_DepCtrlParam%NPRT)
            if(present(RECORDPROC)) then
               call For_One_Test(SimBox0,J,processid,RECORDPROC)
            else
               call For_One_Test(SimBox0,J,processid)
            end if
         END DO !end the loop for sampling

        !**** if external after-process provided, we do it
        if(present(AFTRECORD)) call AFTRECORD(SimBox0,m_CtrlParam,m_DepCtrlParam)

        !---clear the momery allocated on device
        call CLEAR_FINNIS_EM_TB_Force_Table_DEV()
        call CLEAR_Cal_NeighboreList2C_DEV()
        call Clear_Globle_Variables_DEV()
        !----

        CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
        C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
        TIME2 = C2-C1

       print *, "RUN TIME: ",TIME2

    return
  end subroutine SubProcess
  !****************************************************************************************

  !****************************************************************************************
  subroutine For_One_Test(SimBox0,J,processid,RECORDPROC)
  !--- INPUT: NMPI, number of prcocesses
  !           processid, the process ID
  !           RECORDPROC, the optional subroutine provided by user to process data at time when the configure to be output
  !    NOTE: the three dummy subroutine must strictly follow the interface given below
  implicit none
  !--- DUMMY VAIRBALES
  type(SimMDBox), intent(in)::SimBox0
  integer,        intent(in)::processid,J
  !--- dunmmy functions
      OPTIONAL::RECORDPROC
      external::RECORDPROC
  !--- interface to the external routines -------------------
   interface
       SUBROUTINE RECORDPROC(JP,ITIME, TIME, SimBox, CtrlParam,DEPCtrlParam,IPROC)
       use MD_CONSTANTS_2010
       use MD_Globle_Variables_2010
       use MD_SIMULAIONBOX_ARRAY_2010
       use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
       use DEPOSITION_TypeDef_CtrlParam_2010
       implicit none
       type(SimMDBox)           ::SimBox
       type(SimulationCtrlParam)::CtrlParam
       type(DepositionCtrlParam)::DEPCtrlParam
       integer::JP,ITIME
       real(KINDDF)::TIME
       integer::IPROC
       END SUBROUTINE RECORDPROC
   end interface
  !--- END INTERFACE --------------------------------

  !--- LOCAL VARIABLES
  integer::I, EXIST
  character*256::SFILE, GFile,GFileU, CFile, CFile0
  real(KINDDF)::TIME



             !*** to prepare the filename for thermal quantities
             call STRCATI(SFILE, m_CtrlParam%f_quantity, "P", processid, 4)
             call STRCATI(SFILE, SFILE, ".", J, 4)
             call STRCATI(CFILE, m_CtrlParam%f_configure, "P", processid, 4)
             call STRCATI(CFILE, CFILE, ".", J, 4)
             call STRCATI(GFILE, m_CtrlParam%f_geometry, "P", processid, 4)
             call STRCATI(GFILE, GFILE, "_", J, 4)
             call STRCATI(GFILEU, m_CtrlParam%f_geometry, "U_P", processid, 4)
             call STRCATI(GFILEU, GFILEU, "_", J, 4)


             !*** to load initial configuration
             I = m_CtrlParam%IT
             TIME = C_ZERO

             !*** if in restart mode, to check if the SWAP file exist
             if(m_CtrlParam%RESTART .GT. 0) then
                INQUIRE(FILE=CFILE, exist=EXIST)
             else
                EXIST = .FALSE.
             end if

              if(m_CtrlParam%RESTART .EQ. 0 .or. .NOT.EXIST) then
                ! *** if we start a new session
                 if(Simbox0%IniCfgID .gt. 0) then
                    call STRCATI(CFile0, SimBox0%IniConfig, "P", processid, 4)
                    call STRCATI(CFile0, CFile0, "_", J, 4)
                    call STRCATI(CFile0, CFile0, ".",  Simbox0%IniCfgID, 4)
                    print *, "Load initial substrate from multi-files ", CFile0(1:len_trim(CFile0))
                    print *, "........ "
                    call Initialize_Config_SimMDBox(SimBox0,fname=CFile0)
                 else
                    print *, "Load initial substrate from single file ", SimBox0%IniConfig(1:len_trim(SimBox0%IniConfig))
                    print *, "........ "
                    call Initialize_Config_SimMDBox(SimBox0)
                 end if
                 call DepositOneAtom_DEPOSITION(SimBox0, m_SimBox,m_CtrlParam,m_DepCtrlParam)
             else
                 !*** if we restart from a previus, restore from previus stop point
                 !*** place the incident atom on surface and given an incident kinetic energy
                 call DepositOneAtom_DEPOSITION(SimBox0, m_SimBox,m_CtrlParam,m_DepCtrlParam)
                 print *, "Restore configure from ", CFILE(1:len_trim(CFILE))
                 print *, "........ "
                 call Restore_Config_SimMDBox(m_SimBox, CFile, I, TIME)
              end if
              !***

              !--- to initialize the GPU variable
              call Initialize_Globle_Variables_DEV(m_SimBox, m_CtrlParam)

              !*** to initialize force table on device
              call INITIALIZE_FINNIS_EM_TB_Force_Table_DEV(m_SimBox, m_CtrlParam, m_ForceTable)

              !***  to give a inital neighbore list
              call INITIALIZE_Cal_NeighboreList2C_DEV(m_SimBox, m_CtrlParam)

             !*** to calculate the force
                call Cal_NeighBoreList_2CP_DEV(m_SimBox, m_CtrlParam)
              call CALFORCE_FINNIS_EM_TB_Force_Table2A_DEV(m_SimBox, m_CtrlParam)

              !--- putout the initial configure
              call CALEPOT_FINNIS_EM_TB_Force_Table2A_DEV(m_SimBox, m_CtrlParam)
              call CopyOut_SimMDBox_DEV( m_SimBox)
              call Putout_Instance_Thermal_Quantities_SimMDBox(I, TIME, SFILE, m_SimBox,1)
              call Putout_Instance_Config_SimMDBox(I, m_CtrlParam%TIMESTPG, TIME, GFILE, m_SimBox)

            !--- if external data process provided, we do it at begining configuration
              if(present(RECORDPROC) ) then
                 call RECORDPROC(J,  -1, TIME, m_SimBox, m_CtrlParam, m_DepCtrlParam,processid)
                 call RECORDPROC(J,   I, TIME, m_SimBox, m_CtrlParam, m_DepCtrlParam,processid)
              endif

            !*** to begine the loop
            DO WHILE(I.LE.m_CtrlParam%ITE)

              !--- use variable time step
               m_CtrlParam%H  =  m_CtrlParam%HMI*(int(I/m_CtrlParam%IHDUP)+1)
               if( m_CtrlParam%H .gt.  m_CtrlParam%HMX)  m_CtrlParam%H =  m_CtrlParam%HMX

              !--- use vairiable the nieghbor-list update frequence
               m_CtrlParam%ITAB = m_CtrlParam%ITABMI*(int(I/m_CtrlParam%IITAB)+1)
               if( m_CtrlParam%ITAB .gt. m_CtrlParam%ITABMX)  m_CtrlParam%ITAB = m_CtrlParam%ITABMX

              !--- move system for one time step
               call For_One_Step(I)

              !--- output required information
               I = I + 1
               TIME = TIME+m_CtrlParam%H*CP_S2PS
               IF(MOD(I, m_CtrlParam%TIMESTPQ).EQ.0) THEN
                  call CALEPOT_FINNIS_EM_TB_Force_Table2A_DEV(m_SimBox, m_CtrlParam)
                  call CopyOut_SimMDBox_DEV( m_SimBox)

                  if(processid .eq. 0) then
                     write(*,*) 'DEP.TIME STEP = ',J, I
                     call Putout_Instance_Thermal_Quantities_SimMDBox(I,TIME, SFILE, m_SimBox,1)
                  else
                     call Putout_Instance_Thermal_Quantities_SimMDBox(I, TIME, SFILE, m_SimBox)
                  end if

                  !--- if external data process provided, we do it at begining configuration
                   if(present(RECORDPROC) ) then
                      call RECORDPROC(J,  I, TIME, m_SimBox, m_CtrlParam, m_DepCtrlParam,processid)
                   endif
                  !---
               END IF

               IF(m_CtrlParam%TIMESTPG .gt. 0) then
                  IF(MOD(I,m_CtrlParam%TIMESTPG).EQ.0) THEN
                   !NOTE: it is assumed that the m_TimestpG is a multiple of m_TimestpQ
                   !      in this way we do not need to copyout the box from device to host
                   !call CopyOut_SimMDBox_DEV( m_SimBox)
                     call Putout_Instance_Config_SimMDBox(I, m_CtrlParam%TIMESTPG, TIME, GFILE, m_SimBox)
                  END IF
               END IF

               IF(m_CtrlParam%TIMESTPSave .gt. 0) then
                  IF(MOD(I,m_CtrlParam%TIMESTPSave).EQ.0) THEN
                   !NOTE: it is assumed that the m_TimestpG is a multiple of m_TimestpQ
                   !      in this way we do not need to copyout the box from device to host
                   !call CALEPOT_FINNIS_EM_TB_Force_Table2A_DEV(m_SimBox, m_CtrlParam)
                   !call CopyOut_SimMDBox_DEV( m_SimBox)
                   !call Cal_thermal_quantities_SimMDBox(m_SimBox)
                   call Archive_Config_SimMDBox(m_Simbox, CFile,I, TIME)
                  END IF
               END IF

            END DO   !end the loop for time step
            call Putout_Instance_Thermal_Quantities_SimMDBox(-1, TIME, SFILE, m_SimBox)


    return
  end subroutine For_One_Test

  !****************************************************************************************



  !****************************************************************************************
  SUBROUTINE For_One_Step(I)
  !***  PORPOSE: to go on for only one step
  !     INPUT:   I, the Ith step
  !
  IMPLICIT NONE
  integer, intent(in)::I
  !Local variables

     !--- Give a prediction
         call Predictor_DEV(I,m_SimBox, m_CtrlParam)

     !--- update the neighbore list
         IF(MOD(I,m_CtrlParam%ITAB).EQ.C_IZERO) THEN
               call Cal_NeighBoreList_2CP_DEV(m_SimBox, m_CtrlParam)
         ENDIF

     !--- to calculate the force at present configuration
     !    NOTE: we have no pressure calculation if CALFORCE_ is used
     !          if pressure needed, we should use CALPTENSOR_
     !
        call CALFORCE_FINNIS_EM_TB_Force_Table2A_DEV(m_SimBox, m_CtrlParam)

     !---  the modification when there is E-P coupling
       if(m_CtrlParam%IFEPC.gt.0) then
             call EPC_MODIFICATION_DEV(m_SimBox, m_CtrlParam)
       end if

     !--- to do the correction
         call Correction_DEV(I,m_SimBox, m_CtrlParam)

     !--- the pressure coupling if constant pressure required
         if(m_CtrlParam%IFBPC .ne. 0) then
           call BPC_MODIFICATION_DEV(m_SimBox, m_CtrlParam)
         end if

     !--- to modify the statu of groups according their RSTATU
         !call Keep_Box_Zero_Translation_SimMDBox_DEV(m_Simbox)

     return
  end subroutine For_One_Step
  !****************************************************************************************


  end module DEPOSITION_BAC_TEST
