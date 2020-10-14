  program CalForceTest
 !***  DESCRIPTION:
 !    To test and compare the computation of force using
 !    CPU and GPU
 !
 !    ----  HOU Qing, Apr, 2010

 !***  The modules included ******************************************
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use MD_NeighborsList
 use MD_NeighborsList_GPU
 use MD_FS_Force_Table
 use MD_FS_Force_Table_GPU

 use MD_Globle_Variables
 use MD_SimboxArray_GPU


 !use EM_TB_ForceTable_Ackland_Vitek_PRB41_10324
 use EM_TB_ForceTable_WangJun_W_HE_2010
 use MD_ActiveRegion_GPU
 !********************************************************************
 implicit none
 !_______________________________________________________________________________________
 !
     type(SimMDBox)     ::m_SimBox, m_SimBox1
     type(SimMDCtrl)    ::m_CtrlParam
     type(MDForceTable) ::m_ForceTable
     type(NEIGHBOR_LIST)::m_NList
 !_______________________________________________________________________________________

 integer::I,J, NUM, ERR, IDEV=0, NDEV=1
 character*32::ARGV
 integer::TESTLOOP=100
 character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
 integer::DATE_TIME1(8),DATE_TIME2(8)
 real*4::C1,C2,TIME1, TIME2

     !*** to initialize the DEVICE
         J = COMMAND_ARGUMENT_COUNT()
         if(J.GE.2) then
            call GET_COMMAND_ARGUMENT(2,ARGV)
            read(ARGV, *) IDEV
         end if

         if(J.GE.3) then
            call GET_COMMAND_ARGUMENT(3,ARGV)
            read(ARGV, *) NDEV
         end if
         call Initialize_DEVICES( IDEV, NDEV)

     !*** to read in control parameters from IniFile
          call Initialize_Globle_Variables_old(m_SimBox, m_CtrlParam)
          write(*,*) "Number of particles: ",m_SimBox%NPRT

     !*** to initialize the force-potential table
         call Register_Interaction_Table(m_SimBox, m_CtrlParam, m_ForceTable)

     !*** to initialize the configuration
          call Initialize_SimMDBox(m_SimBox,2)
          call Initialize_Config_SimMDBox(m_SimBox)

         !--- to initialize the GPU variable
          call Initialize_Globle_Variables_DEV(m_SimBox, m_CtrlParam)

         m_CtrlParam%IFPD =1
         call Copy_SimMDBox(m_SimBox, m_SimBox1)

     !***  to give an inital neighbore list
         call INITIALIZE_NeighboreList_DEV(m_SimBox, m_CtrlParam)
         call Cal_NeighBoreList_DEV(m_SimBox, m_CtrlParam, List=m_NList)
         write(*,*) "Neighbor-list created.."

     !*** to calculate the force in CPU
         call INITIALIZE_FS_Force_Table(m_SimBox, m_CtrlParam, m_ForceTable)

         CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
         C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000
         !DO I=1, TESTLOOP
            call CALFORCE_FS_Force_Table2(m_SimBox, m_NList)
         !END DO
         CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
         C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
         TIME1 = C2-C1
         write(*,*) "The elapsed time for CPU-FORCE 2C:", TIME1

     !*** to calculate the force in GPU
 2000    call INITIALIZE_FS_Force_Table_DEV(m_SimBox, m_CtrlParam, m_ForceTable)

         CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
         C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000
         !call DeActive_All_DEV(m_SimBox1, m_CtrlParam)
         DO I=1, TESTLOOP
            !call Active_All_DEV(m_SimBox1, m_CtrlParam)
            call CALFORCE_FS_Force_Table2A_DEV(m_SimBox1, m_CtrlParam)
           !call CALPTENSOR_FS_Force_Table2A_DEV(m_SimBox1, m_CtrlParam,CopyOut=0)
           !call CALEPOT_FS_Force_Table2A_DEV(m_SimBox1, m_CtrlParam)
         END DO

         CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
         C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
         TIME2 = C2-C1
         write(*,*) "The elapsed time for GPU-FORCE 2C:", TIME2

         call CALPTENSOR_FS_Force_Table2A_DEV(m_SimBox1, m_CtrlParam)
         call CALEPOT_FS_Force_Table2A_DEV(m_SimBox1, m_CtrlParam)
         call CopyOut_SimBox_FOR_DEV(m_SimBox1)
         call CopyOut_SimBox_EPOT_DEV(m_SimBox1)

       !*** to check if CPU and GPU pressure tensor are the same
         ERR = 0
          DO I=1,5 !m_SimBox%NPRT
             if(m_SimBox%EPOT(I) .ne. m_SimBox1%EPOT(I) ) then
                write(*,*) "Diff EPOT: ",I, m_SimBox%EPOT(I), m_SimBox1%EPOT(I)
                !err = 1
                !exit
              end if
          END DO

         if(any(m_SimBox%VTENSOR .ne. m_SimBox1%VTENSOR)) then
            write(*,*) "Diff TENSOR:"
            DO I=1, 3
            DO J=1, 3
               write(*,*) I, J, m_SimBox%VTENSOR(I,J), m_SimBox1%VTENSOR(I,J)
            END DO
            END DO
         end if
       !*** to check if CPU and GPU Force are the same
 3000     ERR = 0
          DO I=1,m_SimBox%NPRT

             if(any(m_SimBox%FP(I,:) .ne. m_SimBox1%FP(I,:)) ) then
                write(*,*) "Diff at: ",I
                DO J=1,3
                   write(*,*) I, J, m_SimBox%FP(I,J), m_SimBox1%FP(I,J)
                END DO
                err = 1
                exit
              end if
          END DO

         if(ERR) then
            write(*,*) "TEST NOT PASSED"
         else
            write(*,*) "TEST PASSED"
         end if

         call CLEAR_NeighboreList_DEV()
         call CLEAR_FS_Force_Table_DEV()
         call Clear_Globle_Variables_DEV()
         call End_DEVICES()

   stop
 end program CalForceTest
 !****************************************************************************************




