 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is to rearrange the data (.Intervals) created from EvevntSearch_time
 !
 !                  SEE ALSO ____________________________________________________________________________
 !                       EvevntSearch_time.F90
 !
 !
 !


 !**********************************************************************************
 Program EventSearch_Time_Ex1
 use MD_SimBoxArray_ToolShell_14_CPU
 implicit none
    character*256      ::InFname  = "" 
    integer            ::IT, IB, IC0, IC1, IA, INB, NS, hFile0=10, hFile1=11
    real*8             ::TIME1, TIME2, TIME3, V1(3), V2(3), MU, TT, TSTEP=1.D-4, BINSTP=0.1D0
    integer, parameter::m_NA = 4, m_NT=29
    real(kinddf), parameter::m_CA(m_NA) = (/1.D0, 0.707D0, 0.50D0, 0.D0/)
    real(kinddf)::m_DT(m_NT) = (/0.1D0, 0.25D0, 0.5D0, 0.75D0, 1.0D0, 1.25D0, 1.5D0, 1.75D0, 2.D0, 2.5D0, 3.D0, 3.5D0, 4.D0, 4.5D0, 5.D0, 6.D0, 7.D0, 8.D0, 10.D0, 20.D0, 30.D0, &
                                          40.D0, 50.D0, 60.D0, 70.D0, 80.D0, 90.D0, 100.D0, 200.D0/)
    real(kinddf), parameter::m_EPS   = 5.D-3

    real(kinddf)::DT
    integer     ::ITEST = 1, IBTEST = 2,  IT1, IA1, I, HIS(m_NT,m_NA), OH(m_NT)

       call ExtractCmdlineFileNames_Globle_Variables() 
       InFname = gm_cfgFileName(1)
       open(unit=hFile0, file=InFname, status='old')
       open(unit=hFile1, file=InFname(1:len_trim(InFname))//".Dir1")

       !m_DT = 0
       !do I=1, m_NT
       !   m_DT(I) = I*BINSTP
       !end do  
       TT = 0
       HIS = 0
       OH  = 0
       do while(.true.)
         read(hFile0, *, END=100) IT, IB, TIME1, TIME2, TIME3, V1, V2, MU
         DT = TIME3*TSTEP
         IT1= size(m_DT)
         do I=m_NT, 2, -1
            if(DT .gt. m_DT(I)) then
               exit 
            end if
            IT1= I-1
         end do
         if(IT1 .ge. m_NT .or. IT1 .le. 0) then
            print *, DT, TIME1, TIME2, TIME3, V1, V2, MU
         !write(*, *) IT, IB, TIME1, TIME2, TIME3, DT, IT1
            pause 
         end if
         IA1 = m_NA + 1
         do I=1, m_NA
            if(dabs(dabs(MU)-m_CA(I)) .le. m_EPS) then
               IA1 = I 
               exit
            end if
         end do


         if(IA1.le. m_NA) then 
            HIS(IT1,IA1) = HIS(IT1,IA1) + 1
         else
            OH(IT1) = OH(IT1)+ 1   
         end if   
         !if(IT .eq. ITEST .and. IB .eq. IBTEST) then
         !    TT = TT + TIME2
         !    write(hFile1, fmt="(I4, 1x, I8, 1x, 14(1pE16.7,1x))") IT, IB, TT, TT*DT, V1, V2, MU
         ! end if
       end do   
 100  do I=1, m_NT
          write(hFile1, fmt="(I4, 1x, 1pE14.4,1x, 14(1pE14.4,1x))") I, m_DT(I), HIS(I, 1:m_NA)/m_DT(I), OH(I)/m_DT(I), (sum(HIS(I, 1:m_NA))+OH(I))/m_DT(I) 
       end do 
         
       close(hFile0)
       close(hFile1)        
       stop
 End program EventSearch_Time_Ex1
