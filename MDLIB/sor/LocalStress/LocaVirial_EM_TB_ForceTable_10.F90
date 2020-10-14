 module LOCAL_VIRIAL_EM_TB_Force_Table_2013
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                 To calculate the local virial stress in regions centered at a given point,
 !                 with the sizes of the regions given 
 !
 !                 This program is adpted from module MD_FINNIS_EM_TB_Force_Table_2010
 !                 in MD_FINNIS_EM_TB_Force_Table_2010.f90 
 !                 ______________________________________________________________________________________
 !                 HOU Qing, April, 2013
 !
 !________________________________________________________________________________________________________
 
 
 !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
   use MD_CONSTANTS_2010
   use MD_TYPEDEF_SimMDBox_2010
   use MD_TYPEDEF_SIMULATIONCTRLPARAM_2010
   use MD_TYPEDEF_FORCETABLE_2010
   use MD_NeighborsList_2010
 !-----------------------------------------------
   implicit none
     integer::m_LP_IOFILE =1, m_processid=0

     integer, private::m_INITED = .false.
     integer, parameter, private::mp_MXGROUP = 4
 
     integer,private::m_IFPD(3)                                             ! if minima image invention periodic condition
     real(KINDDF), private::m_HBOX(3)                                       ! half size of box
     real(KINDDF), private::m_BOX(3)                                        ! size of simulation box, could in centimeter, lattice unit and UNIT size
     real(KINDDF), private::m_BOXSHAPE(3,3)                                 ! the box shape, would be used in Rahmann-Parrilino Scheme
     real(KINDDF), dimension(mp_MXGROUP,mp_MXGROUP),private::m_RCUT2        ! the square of RCUT
 
     integer::m_NTRIAL=50
     integer,parameter,private::mp_NUMSUBVOL=60                             ! the number of subvolumes 
     integer,parameter,private::mp_CUBIC=0                                  ! the number of subvolumes 
     integer,parameter,private::mp_SPHERE=1                                 ! the number of subvolumes 
     integer, private::m_SUBREGTYP = mp_SPHERE                              ! the type of the subregion 
     real(KINDDF), private::m_CENTER(3)   = C_ZERO                          ! the center of the subvolume
     real(KINDDF), private::m_SUBVOLSIZE(3,mp_NUMSUBVOL) 
     real(KINDDF), private::m_SUBBOX(6,mp_NUMSUBVOL)                        !the boundary of the subvolume (in cm)  
     real(KINDDF), private::m_SUBVOL(mp_NUMSUBVOL)                          !the volum of the subvolume (in cm)  
     real(KINDDF), dimension(:), allocatable, private::m_DDEN               !the derivetive of F(RHO), for Finnis-Sinlar potential, it is 1/DSQRT(m_DEN)
     real(KINDDF), dimension(:,:,:,:), allocatable, private::m_VTENSOR      !the instant virial 
     real(KINDDF), dimension(3,3), private::m_AVTENSOR0 = 0.D0              !the accumulated virial of the whole box on times
     integer::m_ACCNUM = 0                                                  !the number of accumulationg
     !****
     private::preCALFOR, CALFOR, CAL_LocalStress, INREGION
 contains
 
 
 !*********************************************************************************
 SUBROUTINE CAL_LocalStress(SimBox, CtrlParam, CENTER, TRIAL, VTENSOR)
 !***  PURPOSE:   to calculate local stress.
 !                It assumed that the force tables are availble through the
 !                call to INITIALIZE_FINNIS_EM_TB_Force_Table, and also
 !                the neighbor list is abailable through the call to
 !                CALFORCE_FINNIS_EM_TB_Force_Table2
 !                
 !
 !     INPUT:     SimBox,     the simulation box
 !                CtrlParam:  the control parameters
 !     OUTPUT     m_VTENSOR,  the local stress
 !      
 !     DEPENDENCE: INITIALIZE_FINNIS_EM_TB_Force_Table
 !                 CALFORCE_FINNIS_EM_TB_Force_Table2
 ! 
 !     NOTE:       the pointer m_NLISTP defined in module
 !                 MD_FINNIS_EM_TB_Force_Table_2010 is used.
 !                 be careful to check if the pointer has been 
 !                 allocated.  
 ! 
  use MD_NeighborsList_2010
  use MD_FINNIS_EM_TB_Force_Table_2010,only:CSI=>m_CSI, KPAIR=>m_KPAIR,   &
                                            POTR=>m_POTR, FPOTR=>m_FPOTR, &
                                            POTB=>m_POTB, FPOTB=>m_FPOTB, &
                                            List=>m_NLISTP 
  implicit none
     !--- dummy vaiables
     type(SimMDBox),            intent(inout)::SimBox
     type(SimulationCtrlParam), intent(in)   ::CtrlParam
     real(KINDDF),              intent(in)   ::CENTER(3)
     integer                                 ::TRIAL
     real(KINDDF),dimension(:,:,:,:)         ::VTENSOR
 
     !--- Local variables
     integer::I, J, K, NE=1
     real(KINDDF)::DSIZE(3), TENSOR(3,3)

              !$$---
              if(.not. associated(List)) then
                 write(*,*) "MDPSCU Error: The neighbore list is unavailable in CAL_LocalStress"
                STOP 'The process stop'
              end if
              !$$--- The boxsize could be changed due to pressure
              !$$    we should update it at each time
              m_BOX =  SimBox%ZL
              m_HBOX = SimBox%ZL*C_HALF
              m_BOXSHAPE = SimBox%BOXSHAPE
              m_IFPD = CtrlParam%IFPD
              m_RCUT2 = CtrlParam%RU*CtrlParam%RU

              DSIZE(1) = (SimBox%ZL(1))/mp_NUMSUBVOL
              DSIZE(2) = (SimBox%ZL(2))/mp_NUMSUBVOL
              DSIZE(3) = (SimBox%ZL(3))/mp_NUMSUBVOL
 
              DSIZE = (minval(SimBox%ZL)-maxval(CtrlParam%RU) - SimBox%RR)/dble(mp_NUMSUBVOL)  

              !$$--- To calculate the electron densitis
              call PRECALFOR(SimBox%NPRT, List%KVOIS, List%INDI, SimBox%ITYP, SimBox%XP, &
                             CSI,KPAIR, POTR, FPOTR, POTB, FPOTB, m_CENTER, m_BOX,m_DDEN)
         
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
              !m_DDEN = 0.D0
              !$$--- determine the subvolume boundary
              DO I=1,mp_NUMSUBVOL
                 DO J=1, NE
                    m_SUBVOLSIZE(1,I) =  DSIZE(1)*dble(I) 
                    m_SUBVOLSIZE(2,I) =  DSIZE(2)*dble(I) 
                    m_SUBVOLSIZE(3,I) =  DSIZE(3)*dble(I) 

                    m_SUBBOX(1,I) = -m_SUBVOLSIZE(1,I)*C_HALF + CENTER(1)
                    m_SUBBOX(2,I) =  m_SUBBOX(1,I)+m_SUBVOLSIZE(1,I)

                    m_SUBBOX(3,I) = -m_SUBVOLSIZE(2,I)*C_HALF + CENTER(2)
                    m_SUBBOX(4,I) =  m_SUBBOX(3,I)+m_SUBVOLSIZE(2,I)

                    m_SUBBOX(5,I) = -m_SUBVOLSIZE(3,I)*C_HALF + CENTER(3)
                    m_SUBBOX(6,I) =  m_SUBBOX(5,I)+m_SUBVOLSIZE(3,I)
                    if(m_SUBREGTYP .EQ. mp_CUBIC) then
                       m_SUBVOL(I) = (m_SUBBOX(2,I)-m_SUBBOX(1,I))* &
                                     (m_SUBBOX(4,I)-m_SUBBOX(3,I))* &
                                     (m_SUBBOX(6,I)-m_SUBBOX(5,I))
                    else
                       m_SUBVOL(I) = CP_4PI3*(m_SUBVOLSIZE(1,I)*C_HALF)**3.D0
                    end if

                    !$$--- To get the force
                    call CALFOR(SimBox%NPRT, List%KVOIS,List%INDI, SimBox%ITYP, SimBox%XP,     &
                                CSI,KPAIR, POTR, FPOTR, POTB, FPOTB,  CENTER, m_SUBBOX(:,I), m_DDEN, TENSOR)

                    VTENSOR(:,:,I,TRIAL) = VTENSOR(:,:,I,TRIAL) + TENSOR(:,:)
                    
                 END DO
                 VTENSOR(1:3,1:3,I,TRIAL) =  VTENSOR(1:3,1:3,I,TRIAL)/m_SUBVOL(I)/dble(NE)
              END DO
              
              return
 
 end SUBROUTINE CAL_LocalStress
 !*********************************************************************************


 !*********************************************************************************
 logical function INREGION(XP,CENTER, BOX)
 implicit none
     real(KINDDF),intent(in)::XP(3), CENTER(3), BOX(6)

           select case(m_SUBREGTYP)
           case (mp_CUBIC)
                if(XP(1) .LT. BOX(1) .OR. & 
                   XP(1) .GT. BOX(2) .OR. & 
                   XP(2) .LT. BOX(3) .OR. & 
                   XP(2) .GT. BOX(4) .OR. & 
                   XP(3) .LT. BOX(5) .OR. & 
                   XP(3) .GT. BOX(6) ) then
                   INREGION = .false.
               else
                   INREGION = .true.     
               end if  
           case (mp_SPHERE)
               if(sum((XP-CENTER)*(XP-CENTER)) .GT. BOX(1)*BOX(1)) then
                  INREGION = .false.
               else
                  INREGION = .true.     
               end if
           case default
               INREGION = .true.       
           end select
           return
 end function INREGION
 !*********************************************************************************

 !*********************************************************************************
 SUBROUTINE preCALFOR(IM,KVOIS,INDI,ITYP,XP, CSI,KPAIR, POTR, FPOTR, POTB, FPOTB,CENTER, SUBBOX,DDEN)
 !***  PURPOSE:   to begine calculate the pairwise potential m_ER and the electron density m_DEN
 !     INPUT:     IM:     the number of particles concerned
 !                KVOIS:  the number of for each atom in style1
 !                INDI:   the index for neighbores
 !                ITYPE:  the type of atom corresponding to INDI
 !                XP:     position of the atoms
 !
 !     OUTPUT     DDEN
 
 !
  use MD_CONSTANTS_2010
     implicit none
     !--- DUMMY VARIABLES
     integer, intent(in)::IM
     integer, dimension(:), intent(in)::KVOIS
     integer, dimension(:,:), intent(in)::INDI
     integer, dimension(:), intent(in)::ITYP
     real(KINDDF), dimension(:,:), intent(in)::XP

     real(KINDDF),intent(in)::CSI
     integer, dimension(:,:),intent(in)::KPAIR
     real(KINDDF),dimension(:,:),intent(in)::POTR, FPOTR, POTB, FPOTB
     real(KINDDF),intent(in)::CENTER(3),SUBBOX(6)
     real(KINDDF),dimension(:)::DDEN
 
     !---Local variables
     integer::I,J,K,KK, IW, IIW, KTAB,KTAB1, K1, ITYPI, ITYPJ
     real(KINDDF)::SK,DK,EXP1,EXP2
     real(KINDDF)::SEP(3), DXYZ(3), R2, R
 
 
          DDEN = C_ZERO
 
          DO I=1, IM

             !if(.not.INREGION(XP(I,1:3),CENTER(:), SUBBOX(:)) )then
             !  cycle
             !end if  
             IIW = KVOIS(I)
             ITYPI  = ITYP(I)
             DO IW=1, IIW
                J=INDI(I,IW)

                !if(.not.INREGION(XP(J,1:3),CENTER(:),SUBBOX(:)) )then
                !    cycle
                !end if  
 
                !$$--- To calculate the seperation between particle I and its IWth neighbore
                DO K=1, 3
                   SEP(K) = XP(I,K)-XP(J,K)
                   IF((m_IFPD(K).GT.0) .AND. (DABS(SEP(K)) .GT. m_HBOX(K))) THEN
                      SEP(K) = SEP(K) - DSIGN(m_BOX(K),SEP(K))
                   END IF
                END DO
 
                !$$--- NOTE: To converte the length-unit into absolute unit (cm)
                !$$          only when the shape of the box is cubic,
                !$$          DXYZ = SEP
                !
                DO K=1, 3
                   DXYZ(K) =  sum(m_BOXSHAPE(K,1:3)*SEP(1:3))
                END DO
                R2  = sum(DXYZ*DXYZ)
 
                !$$--- To calculate electron density on atom I
                ITYPJ=ITYP(J)
                if(R2 .LE. m_RCUT2(ITYPI,ITYPJ)) then
                   KTAB = KPAIR(ITYPI,ITYPJ)
                   KTAB1= KPAIR(ITYPJ,ITYPI)
 
                   R = DSQRT(R2)
                   SK= DSQRT(R)*CSI
                   KK=SK
                   K1=KK+1
                   DK=SK-DBLE(KK)
                   EXP1 = POTR(KTAB,KK)+DK*(POTR(KTAB,K1)-POTR(KTAB,KK))
 
                   DDEN(I)= DDEN(I)+POTB(KTAB,KK) +DK*(POTB(KTAB,K1)-POTB(KTAB,KK))
               end if
           ENDDO
        END DO
 
        !$$*** to get the derivetive of DSQRT(m_DEN)
         DO I=1, IM
            IF(DDEN(I) .GT. 0) then
               !--- NOTE: the force term has been multiplied by 0.5 in
               !          creating the force table.
               DDEN(I) = -C_UN/DSQRT(DDEN(I))
            ELSE
               DDEN(I) = C_ZERO
            END IF
        END DO
        
       RETURN
 END SUBROUTINE preCALFOR
 !*************************************************************************************
 
 !*************************************************************************************
 SUBROUTINE CALFOR(IM, KVOIS,INDI,ITYP,XP, CSI,KPAIR, POTR, FPOTR, POTB, FPOTB, CENTER, SUBBOX, DDEN, VTENSOR)
 !***  PURPOSE:   to begine calculate the force, Newton's third law taken into account
 !     INPUT:     IM:     the number of particles
 !                KVOIS:  the number of neighbores in stle0
 !                INDI:   the index for neighbores
 !                ITYP:   the type of atom corresponding to INDI
 !                XP:     the coordinate of atoms
 !     OUTPUT     
 !                VTENSOR: the virial tensor
 !
 use MD_CONSTANTS_2010
     implicit none
     !--- DUMMY VARIABLES
     integer, intent(in)::IM
     integer, dimension(:),intent(in)::KVOIS
     integer, dimension(:,:), intent(in)::INDI
     integer, dimension(:), intent(in)::ITYP
     real(KINDDF), dimension(:,:),  intent(in)::XP

     real(KINDDF)::CSI
     integer, dimension(:,:)::KPAIR
     real(KINDDF),dimension(:,:)::POTR, FPOTR, POTB, FPOTB
     real(KINDDF),intent(in)::CENTER(3),SUBBOX(6)
     real(KINDDF),dimension(:),intent(in)::DDEN
 
     real(KINDDF), intent(out)::VTENSOR(3,3)

 
     !---Local variables
     integer::I,J,K,K1,N,IW,IIW,KTAB,KTAB1,KK, KK1,ITYPI, ITYPJ
     real(KINDDF)::SK,DK, EXP5,EXP6,EXP7, DENKI, DENKJ, FORTOT, FOR(3)
     real(KINDDF)::SEP(3), DXYZ(3), R2, R
 
       !$$--- now, we have the electron densities on atoms, to calculate the force
         VTENSOR = C_ZERO
         DO I=1, IM
            !$$--- determine if this particle in the subvolume
            if(.not.INREGION(XP(I,1:3),CENTER(:),SUBBOX(:)) )then
               cycle
            end if

            IIW    = KVOIS(I)
            ITYPI  = ITYP(I)
            DENKI =  DDEN(I)
 
            !$$--- to begine scanning the neighbores
            DO IW=1, IIW
               !$$--- To calculate the seperation between particle I and its IWth neighbore
               J=INDI(I,IW)
               DO K=1, 3
                  SEP(K) = XP(I,K)-XP(J,K)
                  IF(m_IFPD(K) .AND. DABS(SEP(K)) .GT. m_HBOX(K)) THEN
                     SEP(K) = SEP(K) - DSIGN(m_BOX(K),SEP(K))
                  END IF
               END DO
 
               !$$--- NOTE: To converte the length-unit into absolute unit (cm)
               !$$          only when the shape of the box is cubic,
               !$$          DXYZ = SEP
               DO K=1, 3
                  DXYZ(K) =  sum(m_BOXSHAPE(K,1:3)*SEP(1:3))
               END DO
               R2  = sum(DXYZ*DXYZ)
 
               !$$--- IF R2 is smaller than cutoff range, to calculate the force
               ITYPJ=ITYP(J)
               if(R2 .LE. m_RCUT2(ITYPI,ITYPJ)) then
                  KTAB  = KPAIR(ITYPI,ITYPJ)
                  KTAB1 = KPAIR(ITYPJ,ITYPI)
                  R = DSQRT(R2)
                  SK= DSQRT(R)*CSI
                  KK = SK
                  KK1 = KK+1
                  DK=SK-dble(KK)
                  !$$-- interpolation of force table
                  !$$-- for  pairwise section
                  EXP5 = FPOTR(KTAB,KK)  + DK*(FPOTR(KTAB,KK1) - FPOTR(KTAB,KK))
 
                  !$$-- for many-body section
                  !$$    Note; the electron density of atom I on atom J could not be the same as
                  !$$    from atom J on atom J, if atom I and J are not the same type.
                  !
                  !$$-- from the type of atom I to the type of atom J
                  EXP6 = FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK1) - FPOTB(KTAB,KK))
 
                  !$$-- from the type of atom J to the type of atom I
                  EXP7 = FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK1) - FPOTB(KTAB1,KK))
 
                  FORTOT= EXP5 + EXP7*DENKI + EXP6* DDEN(J)

                  DO K=1, 3
                     DO K1=1, 3
                        VTENSOR(K,K1) = VTENSOR(K,K1) + DXYZ(K1)*DXYZ(K)*FORTOT*C_HALF
                     END DO
                  END DO

                  !-------------------------------------------------
                end if
           ENDDO
       END DO
 
       RETURN
 END SUBROUTINE CALFOR
 !*************************************************************************************

 !****************************************************************************************
 subroutine RECORD_Instant_LocalStress(JOB,ITIME, TIME, SimBox, CtrlParam)
 !***  DESCRIPTION: to calculate and record the local pressure
 !                  
 !
 !    INPUT:  JOB,    the id of the test
 !            ITIME,  the time step
 !            TIME,   the simulation time
 !            SimBox, the simulation boxs
 !            CtrlParam,   the control parameters for simulation
 !
 !  SEE ALSO:
 !            MD_EM_TB_ForceTable_SHELL_12_GPU
 use MD_SimMDBox_ARRAY_2010
 use MD_TYPEDEF_SIMULATIONCTRLPARAM_2010
 use MD_FINNIS_EM_TB_Force_Table_2010
 implicit none
      type(SimMDBox),dimension(:), intent(in)::SimBox
      type(SimulationCtrlParam),   intent(in)::CtrlParam
      integer,                     intent(in)::JOB,ITIME
      real(KINDDF),                intent(in)::TIME
      !--- local variables
      character*256::GFILE
      integer::hFile, I, J, K
      real(KINDDF)::VV0(3,3), VV(3,3), TC(3),RP

         if(m_ACCNUM .EQ. 0) then 
            if(allocated(m_VTENSOR)) then
               write(*,*) "MDPSCU ERROR: reallocate m_VTENSOR"
               pause
            end if 
            allocate(m_VTENSOR(3,3,mp_NUMSUBVOL,m_NTRIAL*CtrlParam%TOTALBOX))
            m_VTENSOR = 0.D0
         end if 

          K = size(SimBox)
          DO I=1, K
          DO J=1, m_NTRIAL 
             m_ACCNUM    = m_ACCNUM+1
             if(J.GT.1) then
               call RANDOM_NUMBER(RP)
               TC(1) = m_CENTER(1)+(RP-0.5D0)*SimBox(I)%RR
               call RANDOM_NUMBER(RP)
               TC(2) = m_CENTER(2)+(RP-0.5D0)*SimBox(I)%RR
               call RANDOM_NUMBER(RP)
               TC(3) = m_CENTER(3)+(RP-0.5D0)*SimBox(I)%RR
            else
               TC(1) = m_CENTER(1)
               TC(2) = m_CENTER(2)
               TC(3) = m_CENTER(3)
            end if

             call CAL_LocalStress(SimBox(I), CtrlParam, TC, m_ACCNUM, m_VTENSOR(:,:,:,:))
             m_AVTENSOR0 = m_AVTENSOR0 + SimBox(I)%VTENSOR
          END DO
          END DO
 
          call STRCATI(GFILE, CtrlParam%f_others(m_LP_IOFILE), "P", m_processid, 4)
          call STRCATI(GFILE, GFILE, "_", JOB, 4)
          call STRCATI(GFILE, GFILE, ".", ITIME/CtrlParam%TIMESTPG, 4)
          
          call AvailableIOUnit(hFile)
          print *, "Output local stress to: ", GFILE(1:len_trim(GFILE))
          open(UNIT=hFile, file = GFILE, status='unknown')

          DO K=1,mp_NUMSUBVOL 
             write(hFile,fmt=("(1x, 30(E14.5,1x))")) m_SUBVOLSIZE(1,K)*CP_CM2NM, &   !/SimBox(1)%RR,        &
                                                     ((m_VTENSOR(I,J,K,m_ACCNUM)*CP_CGS2KBAR, J=1,3),I=1,3)
                                                     
          END DO
          close(hFile)
          return 
 end subroutine RECORD_Instant_LocalStress
 !****************************************************************************************

 !****************************************************************************************
 subroutine RECORD_Average_LocalStress(SimBox, CtrlParam)
 !***  DESCRIPTION: to calculate and record the local pressure
 !                  
 !
 !    INPUT:  JOB,    the id of the test
 !            ITIME,  the time step
 !            TIME,   the simulation time
 !            SimBox, the simulation boxs
 !            CtrlParam,   the control parameters for simulation
 !
 !  SEE ALSO:
 !            MD_EM_TB_ForceTable_SHELL_12_GPU
 use MD_SimMDBox_ARRAY_2010
 use MD_TYPEDEF_SIMULATIONCTRLPARAM_2010
 use MD_FINNIS_EM_TB_Force_Table_2010
 implicit none
      type(SimMDBox),            intent(in)::SimBox
      type(SimulationCtrlParam), intent(in)::CtrlParam
      !--- local variables
      character*256::GFILE
      integer::hFile, I, J, K
      real(KINDDF)::VV0(3,3), VV(3,3), ERR(3,3)
 
          call STRCATI(GFILE, CtrlParam%f_others(m_LP_IOFILE), "P", m_processid, 4)
          call STRCATI(GFILE, GFILE, ".AV", m_ACCNUM, 4)
          
          call AvailableIOUnit(hFile)
          print *, "Output average local stress to: ", GFILE(1:len_trim(GFILE))
          open(UNIT=hFile, file = GFILE, status='unknown')

          VV0 = m_AVTENSOR0/(SimBox%ZL(1)*SimBox%ZL(2)*SimBox%ZL(3))
          VV0 = VV0/dble(m_ACCNUM)*CP_CGS2KBAR
          DO K=1,mp_NUMSUBVOL 
             !--- average over enables an time
             DO I=1,3
             DO J=1,3
                VV(I,J) = SUM(m_VTENSOR(I,J,K,1:m_ACCNUM))/dble(m_ACCNUM)
             END DO
             END DO

             DO I=1,3
             DO J=1,3
                ERR(I,J) = SUM((m_VTENSOR(I,J,K,1:m_ACCNUM)-VV(I,J))**2)
                ERR(I,J) = ERR(I,J)/dble(m_ACCNUM)
                ERR(I,J) = DSQRT(ERR(I,J))*C_HALF
             END DO
             END DO

             write(hFile,fmt=("(1x, 30(E14.5,1x))")) m_SUBVOLSIZE(1,K)*CP_CM2NM,& !/SimBox%RR,      &
                                                      (VV0(1,1)+ VV0(2,2)+ VV0(3,3))/3.D0 , &
                                                      ((VV(I, J)*CP_CGS2KBAR,ERR(I,J)*CP_CGS2KBAR, J=1,3), I=1,3)
                                                     
          END DO
          close(hFile)
          return 
 end subroutine RECORD_Average_LocalStress
 !****************************************************************************************

 !****************************************************************************************
  subroutine Init_CalLocalStress_Tool(SimBox, CtrlParam, FORCETABLE)
  use MD_SimMDBox_ARRAY_2010
  use MD_TYPEDEF_SIMULATIONCTRLPARAM_2010
  use MD_FINNIS_EM_TB_Force_Table_2010
  !*******
      implicit none
      type(SimMDBox)           ::SimBox
      type(SimulationCtrlParam)::CtrlParam
      OPTIONAL::FORCETABLE
      external::FORCETABLE
      interface
          subroutine FORCETABLE(SimBox, CtrlParam, FTable, printout)
          use MD_TYPEDEF_SimMDBox_2010
          use MD_TYPEDEF_SIMULATIONCTRLPARAM_2010
          use MD_TYPEDEF_FORCETABLE_2010
          implicit none
          !--- dummy vaiables
          type(SimMDBox),            intent(in)   ::SimBox
          type(SimulationCtrlParam), intent(in)   ::CtrlParam
          type(MDForceTable),        intent(inout)::FTable
          integer,                   optional     ::printout
         end subroutine FORCETABLE
       end interface
      !--- local variables
         type(MDForceTable)::FTable

      !--- do initialization
         if(.not.present(FORCETABLE)) then
             print *, "MDPSCU Error: to use the tool, you must provide the force-table generator"
             print *, "The process to be stopped"
             stop
         end if

         call FORCETABLE(SimBox, CtrlParam, FTable) 
         call Putout_ForceTable("FORCE.TAB",FTable)
         call INITIALIZE_FINNIS_EM_TB_Force_Table(SimBox, CtrlParam, FTable, RelaseTable=1)    

         return
 end subroutine Init_CalLocalStress_Tool
 !****************************************************************************************

 !****************************************************************************************
 subroutine RECORD_LocalStress_Tool(JOB,ITIME, TIME, SimBox, CtrlParam)
 !***  DESCRIPTION: to calculate and record the local pressure
 !                  
 !
 !    INPUT:  JOB,    the id of the test
 !            ITIME,  the time step
 !            TIME,   the simulation time
 !            SimBox, the simulation boxs
 !            CtrlParam,   the control parameters for simulation
 !
 !  SEE ALSO:
 !            MD_EM_TB_ForceTable_SHELL_12_GPU
 
  use MD_SimMDBox_ARRAY_2010
  use MD_TYPEDEF_SIMULATIONCTRLPARAM_2010
  use MD_FINNIS_EM_TB_Force_Table_2010
  implicit none
      type(SimMDBox), dimension(:), intent(in)::SimBox
      type(SimulationCtrlParam),    intent(in)::CtrlParam
      integer,                      intent(in)::JOB,ITIME
      real(KINDDF),                 intent(in)::TIME
  
      !--- local variables

          if(ITIME .LT. 0) then
             if(.not.allocated(m_DDEN)) then
                allocate(m_DDEN(SimBox(1)%NPRT*size(SimBox)) )
             end if
             return
          end if

          if(.not. associated(m_NListP)) then
              allocate(m_NLISTP)
          end if


          call Cal_NeighboreList2C(SimBox(1), CtrlParam, m_NListP)
          call CALFORCE_FINNIS_EM_TB_Force_Table2(SimBox(1), m_NListP)
          call RECORD_Instant_LocalStress(JOB,ITIME, TIME, SimBox, CtrlParam) 
      
          return
 
 end subroutine RECORD_LocalStress_Tool
 !****************************************************************************************
 end module LOCAL_VIRIAL_EM_TB_Force_Table_2013
 
 
 
 
 
