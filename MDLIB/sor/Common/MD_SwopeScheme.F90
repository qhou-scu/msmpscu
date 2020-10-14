  module MD_SWOPE_Scheme
  !**** DESCRIPTION: Use the differetial scheme proposed by Swope et al,J.Chem.Phys.76(1982)637-649
  !                  Berendson pressure coupling is realized in this module, HOWEVER,
  !                  the coordinate system is assumed to be Cartesian.
  !
  !                  To consider non-cubic box, Gear Scheme may be needed
  !                  ______________________________________________________
  !**** HISTORY:     2010-03 (HOU Qing), created the first version
  !
  !                  2018-05-31(HOU Qing),
  !                          remove the FIXTEMP operation in Correction subroutine.
  !                          Insteadly, the local temperature for atoms is handled
  !                          in  MD_LocalTempMethod module.

  use MD_CONSTANTS
  implicit none
  !**********************

  contains

  !************************************************************************************************
  subroutine Predictor(ITIME, SimBox, CtrlParam)
  !***  PURPOSE:   to predict the postion and velocity of the atomss in the first half time step
  !     INPUT:     ITIME,     the time step
  !                SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !     OUTPUT     SimBox,    the simulation box with new position and velocity

  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none
  !----   DUMMY Variables
          integer, intent(in)::ITIME
          type(SimMDBox)     ::SimBox
          type(SimMDCtrl)    ::CtrlParam

  !----   Local variables
          integer::I, I1, I2, J, K
          logical::LVAR
          real(KINDDF)::H, HS2, H2S2, BOXSIZE(3), LBOX(3), UBOX(3), CM0, DIS(3)
          real(KINDDF)::MXDIS, DIST, MXDIS0
          integer::IT, IVSECT

  !
  !***
  !
         if(ITIME .ge. CtrlParam%DAMPTIME0 .and.   &
            ITIME .le. CtrlParam%DAMPTIME0 + CtrlParam%DAMPTIME1-1) then
             do J=1, SimBox%NPRT
                do K=1,3
                   if(SimBox%XP1(J,K)*SimBox%FP(J,K) < C_ZERO) then
                      SimBox%XP1(J,K) = C_ZERO
                   end if
                end do
             end do
         end if
         BOXSIZE(1:3) = SimBox%ZL(1:3)
         LBOX(1:3) = SimBox%BOXLOW(1:3)
         UBOX(1:3) = SimBox%BOXUP(1:3)

         !--- if automatical variable time step is required, do it
         if(CtrlParam%IHDUP .LT. 0) then
         if(MOD(ITIME,IABS(CtrlParam%IHDUP)).EQ. C_IUN .OR. &
            IABS(CtrlParam%IHDUP) .EQ. C_IUN) then
            MXDIS0 = CtrlParam%DMX*CtrlParam%DMX
            H = CtrlParam%HMX
            HS2  = H*C_HALF
            H2S2 = H*HS2

            do while(.TRUE.)
               MXDIS = C_ZERO
               LVAR = .FALSE.
               do I=1, SimBox%nGroup
                  I1 = SimBox%IPA(I)
                  I2 = SimBox%IPA(I+1)-1
                  CM0 = H2S2/SimBox%CM(I)
                  do J=I1, I2
                     if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                       DIS(1:3) = H*SimBox%XP1(J,1:3) + SimBox%FP(J,1:3)*CM0
                       DIST = SUM(DIS*DIS)
                       if(DIST.GT.MXDIS) MXDIS = DIST
                       if(MXDIS .GT. MXDIS0) then
                          LVAR = .TRUE.
                          H = C_HALF*H
                          HS2  = H*C_HALF
                          H2S2 = H*HS2
                          exit
                       end if
                     end if
                   end do !End loop for atoms in a group
                   if(LVAR) exit
                 end do ! End loop for group

                 !--- check if need to change time step
                if(.NOT.LVAR) then
                   CtrlParam%H = H
                   exit
                end if
            end do
         end if
         end if


         H = CtrlParam%H
         HS2  = H*C_HALF
         H2S2 = H*HS2

         !--- move for a step
         do I=1, SimBox%nGroup
               I1 = SimBox%IPA(I)
               I2 = SimBox%IPA(I+1)-1
               CM0 = SimBox%CM(I)
               !$$--- for X direction
               if(IAND(SimBox%PROP(I),CP_STATU_FIXPOSX) .eq. C_IZERO) then
                  do J=I1, I2
                     if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                        K=1
                        DIST = H*SimBox%XP1(J,K) + H2S2*SimBox%FP(J,K)/CM0
                        SimBox%DIS(J,K) = SimBox%DIS(J,K) + DIST
                        SimBox%XP(J,K) = SimBox%XP(J,K)   + DIST
                        !SimBox%XP(J,K)  = SimBox%XP(J,K)  + H*SimBox%XP1(J,K) + H2S2*SimBox%FP(J,K)/CM0
                        !--- check periodic boundary condition
                        if(SimBox%XP(J,K) .GT. UBOX(K)) then
                           if(CtrlParam%IFPD(K)) then
                              SimBox%XP(J,K) = SimBox%XP(J,K) - BOXSIZE(K)
                           else
                             SimBox%STATU(J) =  IOR(CP_STATU_OUTOFBOX,CP_STATU_REFLECT)
                           endif
                        else if (SimBox%XP(J,K) .LT. LBOX(K)) then
                                if(CtrlParam%IFPD(K)) then
                                   SimBox%XP(J,K) = SimBox%XP(J,K) + BOXSIZE(K)
                                else
                                   SimBox%STATU(J) = IOR(CP_STATU_OUTOFBOX,CP_STATU_TRANSMIT)
                                endif
                        end if
                        if(IAND(SimBox%STATU(J),CP_STATU_FIXVELX) .eq. C_IZERO) then
                           SimBox%XP1(J,K) = SimBox%XP1(J,K) + HS2*SimBox%FP(J,K)/CM0
                        end if
                     end if
                  end do !End loop for atom in a group
               end if

               !$$--- for Y direction
               if(IAND(SimBox%PROP(I),CP_STATU_FIXPOSY) .eq. C_IZERO) then
                  do J=I1, I2
                     if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                        K=2
                        DIST = H*SimBox%XP1(J,K) + H2S2*SimBox%FP(J,K)/CM0
                        SimBox%DIS(J,K) = SimBox%DIS(J,K) + DIST
                        SimBox%XP(J,K) = SimBox%XP(J,K)   + DIST
                        !SimBox%XP(J,K)  = SimBox%XP(J,K)  + H*SimBox%XP1(J,K) + H2S2*SimBox%FP(J,K)/CM0
                        !--- check periodic boundary condition
                        if(SimBox%XP(J,K) .GT. UBOX(K)) then
                           if(CtrlParam%IFPD(K)) then
                              SimBox%XP(J,K) = SimBox%XP(J,K) - BOXSIZE(K)
                           else
                             SimBox%STATU(J) =  IOR(CP_STATU_OUTOFBOX,CP_STATU_REFLECT)
                           endif
                        else if (SimBox%XP(J,K) .LT. LBOX(K)) then
                                if(CtrlParam%IFPD(K)) then
                                   SimBox%XP(J,K) = SimBox%XP(J,K) + BOXSIZE(K)
                                else
                                   SimBox%STATU(J) = IOR(CP_STATU_OUTOFBOX,CP_STATU_TRANSMIT)
                                endif
                        end if
                        if(IAND(SimBox%STATU(J),CP_STATU_FIXVELY) .eq. C_ZERO) then
                           SimBox%XP1(J,K) = SimBox%XP1(J,K) + HS2*SimBox%FP(J,K)/CM0
                        end if
                     end if
                  end do !End loop for atom in a group
               end if

               !$$--- for Z direction
               if(IAND(SimBox%PROP(I),CP_STATU_FIXPOSZ) .eq. C_ZERO) then
                  do J=I1, I2
                     if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                        K=3
                        DIST = H*SimBox%XP1(J,K) + H2S2*SimBox%FP(J,K)/CM0
                        SimBox%DIS(J,K) = SimBox%DIS(J,K) + DIST
                        SimBox%XP(J,K) = SimBox%XP(J,K)   + DIST
                        !SimBox%XP(J,K)  = SimBox%XP(J,K)  + H*SimBox%XP1(J,K) + H2S2*SimBox%FP(J,K)/CM0
                        !--- check periodic boundary condition
                        if(SimBox%XP(J,K) .GT. UBOX(K)) then
                           if(CtrlParam%IFPD(K)) then
                              SimBox%XP(J,K) = SimBox%XP(J,K) - BOXSIZE(K)
                           else
                             SimBox%STATU(J) =  IOR(CP_STATU_OUTOFBOX,CP_STATU_REFLECT)
                           endif
                        else if (SimBox%XP(J,K) .LT. LBOX(K)) then
                                if(CtrlParam%IFPD(K)) then
                                   SimBox%XP(J,K) = SimBox%XP(J,K) + BOXSIZE(K)
                                else
                                   SimBox%STATU(J) = IOR(CP_STATU_OUTOFBOX,CP_STATU_TRANSMIT)
                                endif
                        end if
                        if(IAND(SimBox%STATU(J),CP_STATU_FIXVELZ) .eq. C_IZERO) then
                           SimBox%XP1(J,K) = SimBox%XP1(J,K) + HS2*SimBox%FP(J,K)/CM0
                        end if
                     end if
                  end do !End loop for atom in a group
               end if


           end do! End loop for group

        return
  end subroutine Predictor
  !************************************************************************************************

  !************************************************************************************************
  subroutine Correction(ITIME, SimBox, CtrlParam)
  !***  PURPOSE:   to correct the velocity of the atoms in the second half time step
  !     INPUT:     ITIME,     the time step
  !                SimBox,    the simulation box
  !                CtrlParam, the control parameters
  !     OUTPUT     SimBox,    the simulation box with new position and velocity
  !
  use MD_CONSTANTS
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  implicit none
  !----   DUMMY Variables
          integer, intent(in)::ITIME
          type(SimMDBox)     ::SimBox
          type(SimMDCtrl)    ::CtrlParam

  !----   Local variables
          integer::I, I1, I2, J, K
          real(KINDDF)::H, HS2, CM0, TT, VSCAL
  !
  !***
  !
          TT  = CtrlParam%TI*C_UPF*CP_KB
          H   = CtrlParam%H
          HS2 = H*C_HALF

          do I=1, SimBox%nGroup
             I1  = SimBox%IPA(I)
             I2  = SimBox%IPA(I+1)-1
             CM0 = SimBox%CM(I)

             !$$--- correct the velocity in X direction
             if(IAND(SimBox%PROP(I),CP_STATU_FIXPOSX) .eq. C_IZERO .and. &
                IAND(SimBox%PROP(I),CP_STATU_FIXVELX) .eq. C_IZERO) then
                do J=I1, I2
                   if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                      SimBox%XP1(J,1) = SimBox%XP1(J,1) + HS2*SimBox%FP(J,1)/CM0
                   end if
                end do
             end if

             !$$--- correct the velocity in Y direction
             if(IAND(SimBox%PROP(I),CP_STATU_FIXPOSY) .eq. C_IZERO .and. &
                IAND(SimBox%PROP(I),CP_STATU_FIXVELY) .eq. C_IZERO) then
                do J=I1, I2
                   if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                      SimBox%XP1(J,2) = SimBox%XP1(J,2) + HS2*SimBox%FP(J,2)/CM0
                   end if
                end do
             end if

             !$$--- correct the velocity in Z direction
             if(IAND(SimBox%PROP(I),CP_STATU_FIXPOSZ) .eq. C_IZERO .and. &
                IAND(SimBox%PROP(I),CP_STATU_FIXVELZ) .eq. C_IZERO) then
                do J=I1, I2
                   if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                      SimBox%XP1(J,3) = SimBox%XP1(J,3) + HS2*SimBox%FP(J,3)/CM0
                   end if
                end do
             end if

             !$$--- at the same time we update the current kinetice energy
             do J=I1, I2
                if(IAND(SimBox%STATU(J),CP_STATU_ACTIVE) .eq. CP_STATU_ACTIVE) then
                   !$$--- at the same time we update the current kinetice energy
                   SimBox%EKIN(J) = C_HALF*CM0*(SimBox%XP1(J, 1)**2+SimBox%XP1(J, 2)**2+SimBox%XP1(J, 3)**2)
                end if
             end do
          end do

         return
  end subroutine Correction
  !************************************************************************************************


  end module MD_SWOPE_Scheme
