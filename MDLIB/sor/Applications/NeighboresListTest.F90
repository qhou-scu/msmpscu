  program NeighboresListTest
 !***  DESCRIPTION:
 !    To test and compare the computation of NeighboreList using
 !    CPU and GPU
 !
 !    ----  HOU Qing, Apr, 2010

 !***  The modules included ******************************************
 use MD_TYPEDEF_SimMDBox
 use MD_TYPEDEF_SimMDCtrl
 use MD_NeighborsList
 use MD_NeighborsList_GPU

 use MD_Globle_Variables
 use MD_Globle_Variables_GPU

 !********************************************************************
 implicit none
 !_______________________________________________________________________________________
 !
     type(SimMDBox)     ::m_SimBox
     type(SimMDCtrl)    ::m_CtrlParam
     type(NEIGHBOR_LIST)::m_NListGPU, m_NList, m_NList1
 !_______________________________________________________________________________________

 integer::I,J, NUM, IDEV=0, NDEV=1, ERR
 character*32::ARGV
 integer::TESTLOOP=100
 character*12::REAL_CLOCK1(3), REAL_CLOCK2(3)
 integer::DATE_TIME1(8),DATE_TIME2(8)
 real*4::C1,C2,TIME1, TIME2
     !***
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
     !*** to initialize the configuration
          call Initialize_SimMDBox(m_SimBox,2)
          call Initialize_Config_SimMDBox(m_SimBox)

         !--- to initialize the GPU variable
          call Initialize_Globle_Variables_DEV(m_SimBox, m_CtrlParam)

         m_CtrlParam%IFPD = 1
         m_CtrlParam%IFPD(3) = 1
         !call omp_set_num_threads(8)

     !________________________________________________________________________________________
     !***  to run CPU neighborelist with cells
   2000   CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
          C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000
          !---
          DO I=1, 1 !TESTLOOP
               call Cal_NeighboreList2C(m_SimBox, m_CtrlParam, m_NList)
            END DO
          CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
          C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
          TIME2 = C2-C1
          write(*,*) "The elapsed time for CPU with cells:", TIME2

     !***  to run GPU neighborelist with cell
   3000   call INITIALIZE_NeighboreList_DEV(m_SimBox, m_CtrlParam)

          CALL DATE_AND_TIME (REAL_CLOCK1(1), REAL_CLOCK1(2), REAL_CLOCK1(3), DATE_TIME1)
          C1 = DATE_TIME1(8)+DATE_TIME1(7)*1000+DATE_TIME1(6)*60*1000+DATE_TIME1(5)*3600*1000

          DO I=1, TESTLOOP
               call Cal_NeighBoreList_DEV(m_SimBox, m_CtrlParam)
            END DO
          CALL DATE_AND_TIME (REAL_CLOCK2(1), REAL_CLOCK2(2), REAL_CLOCK2(3), DATE_TIME2)
          C2 = DATE_TIME2(8)+DATE_TIME2(7)*1000+DATE_TIME2(6)*60*1000+DATE_TIME2(5)*3600*1000
          TIME1 = C2-C1
          write(*,*) "The elapsed time for GPU-2C with cells:", TIME1
          call Cal_NeighBoreList_DEV(m_SimBox, m_CtrlParam, List=m_NListGPU)

      !*** to check if CPU and GPU Lists are the same
 10000    ERR = 0

          DO I=1,m_SimBox%NPRT
             if(m_NListGPU%KVOIS(I) .NE. m_NList%KVOIS(I) ) then
                write(*,*) "Different neighbore number ",I, m_NListGPU%KVOIS(I),m_NList%KVOIS(I)

                DO J=1, max(m_NListGPU%KVOIS(I),m_NList%KVOIS(I))
                    write(*,*) "neighbores ",I, J, m_NListGPU%INDI(I,J),m_NList%INDI(I,J)
                END DO

                ERR = 1
             else
                !write(*,*) "Number of neighbores ",m_NListGPU%KVOIS(I)
                DO J=1, m_NListGPU%KVOIS(I)
                   if(m_NListGPU%INDI(I,J) .ne. m_NList%INDI(I,J)) then
                       write(*,*) "Different neighbore ",I, J, m_NListGPU%INDI(I,J),m_NList%INDI(I,J)
                       ERR = 1
                   end if
                END DO
                if(ERR) exit
             end if

             if(ERR) exit
             if(ERR) then
                ERR = 0
                pause
             end if
          END DO

         call Clear_Globle_Variables_DEV()
         call CLEAR_NeighboreList_DEV
         call End_DEVICES()

         if(ERR) then
            write(*,*) "TEST NOT PASSED"
         else
            !*** record the average number of neighbors
            NUM = 0
            DO I=1,m_SimBox%NPRT
               NUM = NUM + m_NList%KVOIS(I)
            END DO
            NUM = NUM/m_SimBox%NPRT


            write(*,*) "TEST PASSED with average number of neighbors",NUM, Get_Number_of_Cells(), Get_Number_of_AtomPerCell()
         end if

   stop
 end program NeighboresListTest
 !****************************************************************************************




