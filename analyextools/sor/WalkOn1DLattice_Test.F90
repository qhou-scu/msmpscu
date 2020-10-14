 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  The program is used to simulate the displacement distribution for an particle walk on 1D lattice.
 !                  The purpose is to test if the distribution following Gaussian distribution.
 !
 !
 !**** USAGE:       _______________________________________________________________________________________
 !
 !                   WalkOn1DLattice.F90.exe 
 !                  ______________________________________________________________________________________
 !**** HISTORY:
 !                  version 1st 2019-06 (Hou Qing, Sichuan university)
 !
 !

 !**********************************************************************************
 program WalkOn1DLattice_Test_main
    use MD_TYPEDEF_SimMDBox
    use RAND32_MODULE
    implicit none
    integer, parameter::NP = 40000 
    integer, parameter::NR = 10 
    integer, parameter::NB = 100 

    type(SimMDBox):: SimBox0, SimBox
    type(MDRecordStamp):: Stamp
    real(KINDDF)::VEC(3), DIS, P 
    integer::CURID(NP)
    integer::I , J, IB, IP, INB, hFile

    real(KINDDF),parameter::BINS = 0.5D0/2.D0**0.5
    integer,     parameter::BINN = 2000
    integer,dimension(BINN)::His
    character*256:: fname="WalkOn1DLattice_Test", OutFname="" 

       !--- initialize jumping sites

       !--- 
       SimBox0%NPRT = NP
       call Initialize_SimMDBox(SimBox0, 2)
       SimBox0%XP = 0.D0
       SimBox0%DIS = 0.D0
       SimBox = SimBox0

       do IB=1, NB 
          do J=1, NR
             do IP=1, NP
                P = DRAND32()-0.5
                SimBox%XP(IP,1) = SimBox%XP(IP,1)+P*BINS
             end do   
          end do

          !--- calculate the distance distribution
          HIS = 0
          do IP = 1, NP
             VEC(1:3) = SimBox%XP(IP, 1:3) - SimBox0%XP(IP, 1:3) 
             DIS      = dsqrt(sum(VEC*VEC))
             INB      = int(DIS/BINS)+ 1
             His(INB) = His(INB) + 1     
          end do
          call STRCATI(OutFname, Fname, "_DIS.", IB, 4)
          call AvailableIOUnit(hFile)
          open(UNIT=hFile, FILE=OutFname)
          write(hFile, *) "! DISHis from "//Fname(1:len_trim(Fname))
          write(hFile, *) "! For CFG=", IB, NR
          do I = 1, size(His)
            write(hFile, fmt="(I8,1x,1PE14.4,1x,I8,1x,1PE14.4)")  I, (I-1)*BINS, His(I), dble(His(I))
          end do 
          close(hFile)
              
      end do   
         
 end program WalkOn1DLattice_Test_main

