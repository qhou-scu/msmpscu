 module LOCAL_CVSTRESS_EM_TB_Force_Table_2013
 !**** DESCRIPTION: _______________________________________________________________________________________
 !                  To calculate the local stress in control-volume model. 
 !                  The forces on atoms is calculated in tight-binding of Finnis-Sinclar form.
 !                  The force calculations are performed by table-lookup methods with the FORCETABLE 
 !                  previuosly calculated. A FORCETABLE should be registered before calling this module.
 !
 !                 This program is adpted from module MD_FINNIS_EM_TB_Force_Table_2010
 !                 in MD_FINNIS_EM_TB_Force_Table_2010.f90 
 !                 ______________________________________________________________________________________
 !                 HOU Qing, April, 2013
 !
 !________________________________________________________________________________________________________
 
 
 !---  USE THE CONCERNED MODULES INCLUDING BASIC TYPE DEFINITIONS ----
   use MD_CONSTANTS_2010
   use MD_TYPEDEF_SIMULAIONBOX_2010
   use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
   use MD_TYPEDEF_FORCETABLE_2010
   use MD_NeighborsList_2010
 !-----------------------------------------------
   implicit none
     integer::m_LP_IOFILE =2, m_processid=0

     integer, private::m_INITED = .false.
     integer, parameter, private::mp_MXGROUP = 4
 
     integer,private::m_IFPD(3)                                                     ! if minima image invention periodic condition
     real(KINDDF), private::m_HBOX(3)                                               ! half size of box
     real(KINDDF), private::m_BOX(3)                                                ! size of simulation box, could in centimeter, lattice unit and UNIT size
     real(KINDDF), private::m_BOXSHAPE(3,3)                                         ! the box shape, would be used in Rahmann-Parrilino Scheme
     real(KINDDF), dimension(mp_MXGROUP,mp_MXGROUP),private::m_RCUT2                ! the square of RCUT
 
     integer,parameter,private::mp_NUMSUBVOL=110                                    ! the number of subvolumes 
     integer,parameter,private::mp_CUBIC=0                                          ! the parameter to indicate a cubic control-volume 
     integer,parameter,private::mp_SPHERE=1                                         ! the parameter to indicate a sphere control-volume 
     integer,parameter,private::mp_NUMFACE=6                                        ! the number of faces for the control volume
     integer, private::m_SUBREGTYP = mp_CUBIC                                       ! the type of the subregion, default is cubic 

     real(KINDDF), private::m_CENTER(3)   = C_ZERO                                  ! the center of the subvolume
     real(KINDDF), private::m_SUBVOLSIZE(mp_NUMSUBVOL) 
     real(KINDDF), private::m_SUBBOX(mp_NUMFACE,mp_NUMSUBVOL)                       !the boundary of the subvolume (in cm)  
     real(KINDDF), private::m_SUBVOL(mp_NUMSUBVOL)                                  !the volum of the subvolume (in cm)  
     real(KINDDF), dimension(:), allocatable, private::m_DDEN                       !the derivetive of F(RHO), for Finnis-Sinlar potential, it is 1/DSQRT(m_DEN)
     real(KINDDF), dimension(mp_NUMFACE,3,mp_NUMSUBVOL), private::m_VTENSOR         !the instant stress 
     real(KINDDF), dimension(mp_NUMFACE,3,mp_NUMSUBVOL), private::m_AVTENSOR = 0.D0 !the accumulated virial on times
     real(KINDDF), dimension(mp_NUMFACE,3), private::m_VTENSOR0                     !the instant virial of the whole box on times
     real(KINDDF), dimension(mp_NUMFACE,3), private::m_AVTENSOR0 = 0.D0             !the accumulated virial of the whole box on times
     integer::m_ACCNUM = 0                                                          !the number of accumulationg
     !****
     private::preCALFOR, CALFOR, CAL_LocalStress, INREGION
 contains
 
 
 !*********************************************************************************
 SUBROUTINE CAL_LocalStress(SimBox, CtrlParam, CENTER, TRIAL, VTENSOR)
 !***  PURPOSE:   to begine calculate the force, Newton's third law taken into account
 !
 !     INPUT:     SimBox,    the simulation box
 !                CtrlParam: the control parameters
 !                List,      the neighbore list
 !     OUTPUT     SimBox,    the simulation box with force updated
 !
  use MD_NeighborsList_2010
  use MD_FINNIS_EM_TB_Force_Table_2010,only:CSI=>m_CSI, KPAIR=>m_KPAIR,   &
                                            POTR=>m_POTR, FPOTR=>m_FPOTR, &
                                            POTB=>m_POTB, FPOTB=>m_FPOTB, &
                                            List=>m_NLISTP 
  implicit none
     !--- dummy vaiables
     type(SimMDBox),           intent(inout)::SimBox
     type(SimulationCtrlParam),intent(in)   ::CtrlParam
     integer                                ::TRIAL
     real(KINDDF),dimension(:,:,:,:)        ::VTENSOR

 
     !--- Local variables
     integer::I
     real(KINDDF)::SS, TENSOR(6,3)

              !$$--- The boxsize could be changed due to pressure
              !$$    we should update it at each time
              m_BOX =  SimBox%ZL
              m_HBOX = SimBox%ZL*C_HALF
              m_BOXSHAPE = SimBox%BOXSHAPE
              m_IFPD = CtrlParam%IFPD
              m_RCUT2 = CtrlParam%RU*CtrlParam%RU

             !$$--- To calculate the electron densitis
              call PRECALFOR(SimBox%NPRT, List%KVOIS, List%INDI, SimBox%ITYP, SimBox%XP, &
                             CSI,KPAIR, POTR, FPOTR, POTB, FPOTB, m_DDEN)

              !$$--- determine the subvolume boundary
              DO I=1,mp_NUMSUBVOL
                 m_SUBVOLSIZE(I) =  C_HALF*SimBox%RR*dble(I) 

                 m_SUBBOX(1,I) = -m_SUBVOLSIZE(I)*C_HALF+m_CENTER(1)
                 m_SUBBOX(2,I) =  m_SUBBOX(1,I)+m_SUBVOLSIZE(I)

                 m_SUBBOX(3,I) = -m_SUBVOLSIZE(I)*C_HALF+m_CENTER(2)
                 m_SUBBOX(4,I) =  m_SUBBOX(3,I)+m_SUBVOLSIZE(I)

                 m_SUBBOX(5,I) = -m_SUBVOLSIZE(I)*C_HALF+m_CENTER(3)
                 m_SUBBOX(6,I) =  m_SUBBOX(5,I)+m_SUBVOLSIZE(I)
                 if(m_SUBREGTYP .EQ. mp_CUBIC) then
                    m_SUBVOL(I) = (m_SUBBOX(2,I)-m_SUBBOX(1,I))* &
                                  (m_SUBBOX(4,I)-m_SUBBOX(3,I))* &
                                  (m_SUBBOX(6,I)-m_SUBBOX(5,I))
                 else
                    m_SUBVOL(I) = CP_4PI3*(m_SUBVOLSIZE(I)*C_HALF)**3.D0
                 end if


                 !$$--- To get the force
                 call CALFOR(SimBox%NPRT, List%KVOIS,List%INDI, SimBox%ITYP, SimBox%XP,     &
                          CSI,KPAIR, POTR, FPOTR, POTB, FPOTB,  m_CENTER, m_SUBBOX(:,I), m_DDEN, TENSOR)
                 !$$--- convert the total force on face to force/unit area
                 m_VTENSOR(:,:,I) = TENSOR(:,:)/(m_SUBVOLSIZE(I)*m_SUBVOLSIZE(I))
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
 logical function INSQUARE(X1,X2,LX1,HX1,LX2,HX2)
 implicit none
     real(KINDDF),intent(in)::X1,X2,LX1,HX1,LX2,HX2
    
         if(X1 .GE. LX1 .AND. X1 .LE. HX1 .AND. & 
            X2 .GE. LX2 .AND. X2 .LE. HX2 ) then
              INSQUARE = .true.
         else 
              INSQUARE = .false.
         end if
         return
 end function INSQUARE
 !*********************************************************************************

 !*********************************************************************************
 integer function INTERSECTFACE(XP,VECTOR, BOX) 
 implicit none
     real(KINDDF),intent(in)::XP(3), VECTOR(3), BOX(6)
     real(KINDDF)::DD, TV(3)
     
           INTERSECTFACE = -1 
           !$$--- check intersection on face 1 and 2
           if(VECTOR(1) .lt. C_ZERO) then
              DD = BOX(1)-XP(1)
              TV = VECTOR*(DD/VECTOR(1))+XP 
              if(INSQUARE(TV(2),TV(3),BOX(3), BOX(4), BOX(5),BOX(6)) )then
                 INTERSECTFACE = 1
                 return
              end if 
           else
              DD = BOX(2) - XP(1)
              TV = VECTOR*(DD/VECTOR(1))+XP 
              if(INSQUARE(TV(2),TV(3),BOX(3), BOX(4), BOX(5),BOX(6)) )then
                 INTERSECTFACE = 2
                 return
              end if 
           end if

           !$$--- check intersection on face 3 and 4
           if(VECTOR(2) .lt. C_ZERO) then
              DD = BOX(3)-XP(2)
              TV = VECTOR*(DD/VECTOR(2))+XP 
              if(INSQUARE(TV(3),TV(1),BOX(5), BOX(6), BOX(1),BOX(2)) )then
                 INTERSECTFACE = 3
                 return
              end if 
           else
              DD = BOX(4) - XP(2)
              TV = VECTOR*(DD/VECTOR(2))+XP 
              if(INSQUARE(TV(3),TV(1),BOX(5), BOX(6), BOX(1),BOX(2)) )then
                 INTERSECTFACE = 4
                 return
              end if 
           end if


           !$$--- check intersection on face 5 and 6
           if(VECTOR(3) .lt. C_ZERO) then
              DD = BOX(5)-XP(3)
              TV = VECTOR*(DD/VECTOR(3))+XP 
              if(INSQUARE(TV(1),TV(2),BOX(1), BOX(2), BOX(3),BOX(4)) )then
                 INTERSECTFACE = 5
                 return
              end if 
           else
              DD = BOX(6) - XP(3)
              TV = VECTOR*(DD/VECTOR(3)) + XP
              if(INSQUARE(TV(1),TV(2),BOX(1), BOX(2), BOX(3),BOX(4)) )then
                 INTERSECTFACE = 6
                 return
              end if 
           end if

           return
 end function INTERSECTFACE
 !*********************************************************************************

 !*********************************************************************************
 SUBROUTINE preCALFOR(IM,KVOIS,INDI,ITYP,XP, CSI,KPAIR, POTR, FPOTR, POTB, FPOTB,DDEN)
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
     real(KINDDF),dimension(:)::DDEN
 
     !---Local variables
     integer::I,J,K,KK, IW, IIW, KTAB,KTAB1, K1, ITYPI, ITYPJ
     real(KINDDF)::SK,DK,EXP1,EXP2
     real(KINDDF)::SEP(3), DXYZ(3), R2, R
 
 
          DDEN = C_ZERO
 
          DO I=1, IM
             IIW = KVOIS(I)
             ITYPI  = ITYP(I)
             DO IW=1, IIW
                J=INDI(I,IW)

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
 
     real(KINDDF), dimension(:,:), intent(out)::VTENSOR

     !---Local variables
     integer::I,J,K,K1,N,IW,IIW,KTAB,KTAB1,KK, KK1,ITYPI, ITYPJ, IFACE
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

               !$$--- find out if atom J is out of the box
               !$$    only the atom J out of the box are accounted 
               if(INREGION(XP(J,1:3),CENTER(:),SUBBOX(:)) ) then
                   cycle
               endif

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
                  !
                  !$$-- from the type of atom I to the type of atom J
                  EXP6 = FPOTB(KTAB,KK)  + DK*(FPOTB(KTAB,KK1) - FPOTB(KTAB,KK))
 
                  !$$-- from the type of atom J to the type of atom I
                  EXP7 = FPOTB(KTAB1,KK) + DK*(FPOTB(KTAB1,KK1) - FPOTB(KTAB1,KK))
 
                  FORTOT= EXP5 + EXP7*DENKI + EXP6* DDEN(J)

                  !$$-- find out wich face the force vector intersect
                  IFACE = INTERSECTFACE(XP(I,1:3),-SEP(:), SUBBOX(:))

                  !$$-- Note:IFACE is a value from 1-6
                  VTENSOR(IFACE,:) = VTENSOR(IFACE,:)+FORTOT*SEP(:)
                  
                end if
           ENDDO
       END DO
       !print *, "for finished"
       !pause
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
 use MD_SIMULAIONBOX_ARRAY_2010
 use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
 use MD_FINNIS_EM_TB_Force_Table_2010
 implicit none
      type(SimMDBox), dimension(:), intent(in)::SimBox
      type(SimulationCtrlParam),    intent(in)::CtrlParam
      integer,                      intent(in)::JOB,ITIME
      real(KINDDF),                 intent(in)::TIME
      !--- local variables
      character*256::GFILE
      integer::hFile, I, J, K
      real(KINDDF)::VV0(3,3), VV(3,3)

          K = size(SimBox)
          DO I=1, K
             call CAL_LocalStress(SimBox(I), CtrlParam)
             m_AVTENSOR  = m_AVTENSOR  + m_VTENSOR
             m_ACCNUM    = m_ACCNUM+1
          END DO
 
          call STRCATI(GFILE, CtrlParam%f_others(m_LP_IOFILE), "P", m_processid, 4)
          call STRCATI(GFILE, GFILE, "_", JOB, 4)
          call STRCATI(GFILE, GFILE, ".", ITIME/CtrlParam%TIMESTPG, 4)
          
          call AvailableIOUnit(hFile)
          print *, "Output local stress to: ", GFILE(1:len_trim(GFILE))
          open(UNIT=hFile, file = GFILE, status='unknown')

           
          !VV0 = SimBox(K)%VTENSOR/(SimBox(K)%ZL(1)*SimBox(K)%ZL(2)*SimBox(K)%ZL(3))
          DO K=1,mp_NUMSUBVOL 
             !VV = m_VTENSOR(1:3,1:3,K)/m_SUBVOL(K)
             !write(hFile,fmt=("(1x, 30(E14.5,1x))")) m_SUBVOLSIZE(K)/SimBox(1)%RR,      &
             !                                        VV0(1,1)+ VV0(2,2)+ VV0(3,3)/3.D0, &
             !                                        VV(1,1) + VV(2,2) + VV(3,3)/3.D0,  &
             !                                        ((VV0(I,J),VV(I,J),J=1,3),I=1,3)
             VV(1,1) = m_VTENSOR(2,1,K)-m_VTENSOR(1,1,K) + m_VTENSOR(4,1,K)-m_VTENSOR(3,1,K) + m_VTENSOR(6,1,K)-m_VTENSOR(5,1,K)
             write(hFile,fmt=("(1x, 30(E14.5,1x))")) m_SUBVOLSIZE(K)/SimBox(1)%RR,   VV(1,1)/6.D0*CP_CV,  &
                                                     (m_VTENSOR(I, 1:3,K)*CP_CV, I=1,6)
                                                     
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
 use MD_SIMULAIONBOX_ARRAY_2010
 use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
 use MD_FINNIS_EM_TB_Force_Table_2010
 implicit none
      type(SimMDBox),            intent(in)::SimBox
      type(SimulationCtrlParam), intent(in)::CtrlParam
      !--- local variables
      character*256::GFILE
      integer::hFile, I, J, K
      real(KINDDF)::VV0(3,3), VV(3,3)
 
          call STRCATI(GFILE, CtrlParam%f_others(m_LP_IOFILE), "P", m_processid, 4)
          call STRCATI(GFILE, GFILE, ".AV", m_ACCNUM, 4)
          
          call AvailableIOUnit(hFile)
          print *, "Output average local stress to: ", GFILE(1:len_trim(GFILE))
          open(UNIT=hFile, file = GFILE, status='unknown')

          !VV0 = m_AVTENSOR0/(SimBox%ZL(1)*SimBox%ZL(2)*SimBox%ZL(3))
          !VV0 = VV0/dble(m_ACCNUM)*CP_CV
          DO K=1,mp_NUMSUBVOL 
             !VV = m_AVTENSOR(1:3,1:3,K)/m_SUBVOL(K)
             !VV = VV/dble(m_ACCNUM)*CP_CV
             !write(hFile,fmt=("(1x, 30(E14.5,1x))")) m_SUBVOLSIZE(K)/SimBox%RR,      &
             !                                        VV0(1,1)+ VV0(2,2)+ VV0(3,3)/3.D0, &
             !                                        VV(1,1) + VV(2,2) + VV(3,3)/3.D0,  &
             !                                        ((VV0(I,J),VV(I,J),J=1,3),I=1,3)
             VV(1,1) = m_AVTENSOR(2,1,K)-m_AVTENSOR(1,1,K) + m_AVTENSOR(4,1,K)-m_AVTENSOR(3,1,K) + m_AVTENSOR(6,1,K)-m_AVTENSOR(5,1,K)

             write(hFile,fmt=("(1x, 30(E14.5,1x))")) m_SUBVOLSIZE(K)/SimBox%RR,  VV(1,1)/dble(m_ACCNUM)*CP_CV,    &
                                                     (m_AVTENSOR(I, 1:3, K)/dble(m_ACCNUM)*CP_CV, I=1,6) 
                                                     
          END DO
          close(hFile)
          return 
 end subroutine RECORD_Average_LocalStress
 !****************************************************************************************

 !****************************************************************************************
  subroutine Init_CalLocalStress_Tool(SimBox, CtrlParam, FORCETABLE)
  use MD_SIMULAIONBOX_ARRAY_2010
  use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
  use MD_FINNIS_EM_TB_Force_Table_2010
  !*******
      implicit none
      type(SimMDBox)           ::SimBox
      type(SimulationCtrlParam)::CtrlParam
      OPTIONAL                 ::FORCETABLE
      external                 ::FORCETABLE
      interface
          subroutine FORCETABLE(SimBox, CtrlParam, FTable, printout)
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
      !--- local variables
         type(MDForceTable)::FTable

      !--- do initialization
         if(.not.present(FORCETABLE)) then
             print *, "MDPSCU Error: to use the tool, you must provide the force-table generator"
             print *, "The process to be stopped"
             stop
         end if

         call FORCETABLE(SimBox, CtrlParam, FTable) 
         call INITIALIZE_FINNIS_EM_TB_Force_Table(SimBox, CtrlParam, FTable, RelaseTable=1)    

         allocate(m_NLISTP)
         m_AVTENSOR = 0.D0
         m_ACCNUM   = 0

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
 
  use MD_SIMULAIONBOX_ARRAY_2010
  use MD_TYPEDEF_SIMULAIONCTRLPARAM_2010
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
 end module LOCAL_CVSTRESS_EM_TB_Force_Table_2013
 
 
 
 
 
