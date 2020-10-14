  module MD_LBFGSScheme_GPU
  !**** DESCRIPTION: Use the LBFGS algorithm to search the energy-minimum of of multi-particle system.
  !                  This module usually is used to quench a system to zero temperature
  !
  !                  ______________________________________________________
  !                  HOU Qing, aug, 2016
  !
  use MD_TYPEDEF_SimMDBox
  use MD_TYPEDEF_SimMDCtrl
  use MD_Forceclass_Register_GPU
  implicit none

  contains

  !****************************************************************************************
  subroutine DO_LBFGS_FORSTEPS_DEV(SimBox,CtrlParam, ForceClass, MXNUMSTEPS, IFLAG)
  !***  PORPOSE:  to use LBFGS to quench the system
  !
  !     INPUT:     SimBox,      the simulation box
  !                CtrlParam,   the control parameters
  !                MXNUMSTEPS,  the max number of calling LBFGS
  !     OUTPUT
  !
  use MD_Globle_Variables_GPU
  use LBFGS_MODULE
  implicit none
  !----   DUMMY Variables
          type(SimMDBox), dimension(:)             :: SimBox
          type(SimMDCtrl),               intent(in):: CtrlParam
          type(MDForceClassGPU), target, intent(in):: ForceClass
          integer,                       intent(in):: MXNUMSTEPS
          integer                                  :: IFLAG
  !---   Local variables
          integer::I, IA, MSAVE, IPRINT(2), K, J, NPRT, NPRTP, IERR, IFPD(3), LP0
          real(KINDDF)::OBJF, EPS= 1.0D-5, XTOL= 1.0D-16
          real(KINDDF), dimension(:), allocatable::X, G                  ! the gradient of potential
          real(KINDDF), dimension(:), allocatable::DIAG                  ! the dialoge of Hessian matrix
          real(KINDDF), dimension(:), allocatable::WORK                  ! the working array
          integer, dimension(:), allocatable::IDX, IDY, IDZ
          integer::NIDX, NIDY, NIDZ
          real(KINDDF)::LB(3), UB(3), BS(3), GTOL0

                IPRINT(1)= -1
                IPRINT(2)= 0
                EPS      = 1.0D-4
                XTOL     = 1.0D-16
                GTOL0    = m_GTOL
                LP0      = m_LP
                m_GTOL   = 1.D-2
                m_LP    = 0

                MSAVE    = CtrlParam%LBFGS_MSave
                IFLAG    = 0
                NPRT     = dm_NPRT
                IFPD     = CtrlParam%IFPD
                LB       = SimBox(1)%BOXLOW
                UB       = SimBox(1)%BOXUP
                BS       = SimBox(1)%ZL
                allocate( X(NPRT*3), G(NPRT*3), DIAG(NPRT*3), WORK(NPRT*3*(2*MSAVE +1)+2*MSAVE), &
                          IDX(NPRT), IDY(NPRT), IDZ(NPRT), STAT=IERR)
                if(IERR .gt. 0) then
                    write(*,fmt="(A)") " MDPSCU Error: fail to allocate working space in DO_LBFGS_FORSTEPS_DEV"
                    write(*,fmt="(A)") "               process to be stopped"
                    stop
                end if
                call CopyXPFrom_Devices_to_Host()
                call CopyStatuFrom_Devices_to_Host()

                NIDX = 0
                do IA = 1, NPRT
                   if(iand(hm_STATU(IA), CP_STATU_FIXPOSX) .eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDX = NIDX + 1
                      IDX(NIDX) = IA
                   end if
                end do

                NIDY = 0
                do IA = 1, NPRT
                   if(iand(hm_STATU(IA), CP_STATU_FIXPOSY) .eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDY = NIDY + 1
                      IDY(NIDY) = IA
                   end if
                end do

                NIDZ = 0
                do IA = 1, NPRT
                   if(iand(hm_STATU(IA), CP_STATU_FIXPOSZ) .eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDZ = NIDZ + 1
                      IDZ(NIDZ) = IA
                   end if
                end do

                X(1:NIDX)                      = hm_XP(IDX(1:NIDX),1)
                X(NIDX+1:      NIDX+NIDY)      = hm_XP(IDY(1:NIDY),2)
                X(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = hm_XP(IDZ(1:NIDZ),3)
                NPRTP = NIDX + NIDY + NIDZ

                do I=1, MXNUMSTEPS
                  !$$ --- to calculate potential
                   call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                   call UpdateEpot_ForceClass(SimBox, CtrlParam, ForceClass)
                   call CopyFPFrom_Devices_to_Host()
                   call CopyEPOTFrom_Devices_to_Host()
                   call SynchronizeDevices()

                   OBJF = sum(hm_EPOT(1:NPRT))
                   G(1:NIDX)                      = -hm_FP(IDX(1:NIDX),1)
                   G(NIDX+1:      NIDX+NIDY)      = -hm_FP(IDY(1:NIDY),2)
                   G(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = -hm_FP(IDZ(1:NIDZ),3)
                   call LBFGS(NPRTP,MSAVE, X, OBJF, G,DIAG,IPRINT,EPS,XTOL,WORK,IFLAG)

                   !$$ ---
                   if(IFPD(1)) then
                      do IA=1, NIDX
                         if(X(IA) .lt. LB(1) )then
                            X(IA) = X(IA) + BS(1)
                         else if(X(IA) .gt. UB(1) )then
                            X(IA) = X(IA) - BS(1)
                         end if
                      end do
                   end if

                   if(IFPD(2)) then
                      do IA=NIDX+1, NIDX+NIDY
                         if(X(IA) .lt. LB(2) )then
                            X(IA) = X(IA) + BS(2)
                         else if(X(IA) .gt. UB(2) )then
                            X(IA) = X(IA) - BS(2)
                         end if
                      end do
                   end if

                   if(IFPD(3)) then
                      do IA=NIDX+NIDY+1, NIDX+NIDY+NIDZ
                         if(X(IA) .lt. LB(3) )then
                            X(IA) = X(IA) + BS(3)
                         else if(X(IA) .gt. UB(3) )then
                            X(IA) = X(IA) - BS(3)
                         end if
                      end do
                   end if
                   hm_XP(IDX(1:NIDX), 1) = X(1:           NIDX)
                   hm_XP(IDY(1:NIDY), 2) = X(NIDX+1:      NIDX+NIDY)
                   hm_XP(IDZ(1:NIDZ), 3) = X(NIDX+NIDY+1: NIDX+NIDY+NIDZ)

                   call CopyXPFrom_Host_to_Devices()

                   if(IFLAG .eq.0) then
                      write(*,fmt="(A, I6)") " MDPSCU Message: LBFG succefully performed for steps", I
                      exit
                   end if
                end do
                deallocate(X, G, DIAG, WORK, IDX, IDY, IDZ)

                hm_XP1 = 0.D0
                call CopyXP1From_Host_to_Devices()

                m_GTOL = GTOL0
                m_LP   = LP0
                if(IFLAG .ne. 0) then
                   write(*,fmt="(A, I2)") " MDPSCU Warning: LBFG exit with erro code ", IFLAG
                   write(*,fmt="(A, I2)") "                 in O_LBFGS_FORSTEPS_DEV"
                   call ONWARNING(gm_OnWarning)
                end if

      return
  end subroutine DO_LBFGS_FORSTEPS_DEV
  !****************************************************************************************

  !****************************************************************************************
  subroutine DO_LBFGSB_FORSTEPS_DEV(SimBox,CtrlParam, ForceClass, MXNUMSTEPS, IFLAG)
  !***  PORPOSE:  to use LBFGSB to quench the system
  !
  !     INPUT:     SimBox,      the simulation box
  !                CtrlParam,   the control parameters
  !                MXNUMSTEPS,  the max number of calling LBFGS
  !     OUTPUT
  !
  use MD_Globle_Variables_GPU
  use MATH_LBFGSB_MODULE, only:SETULB
  implicit none
  !----   DUMMY Variables
          type(SimMDBox),dimension(:)       :: SimBox
          type(SimMDCtrl),       intent(in) :: CtrlParam
          type(MDForceClassGPU), intent(in) :: ForceClass
          integer,               intent(in) :: MXNUMSTEPS
          integer                           :: IFLAG
  !---   Local variables
          integer::I, IA, K, J, NPRT, NPRTP, IERR, IFPD(3)
          integer, dimension(:), allocatable::IDX, IDY, IDZ
          integer::NIDX, NIDY, NIDZ
          real(KINDDF)::LB(3), UB(3), BS(3), OBJF
!---    Declare the variables needed by LBFGSB.
!
          character*60::TASK, CSAVE
          logical     :: LSAVE(4)
          integer     :: MSAVE, IPRINT,ISAVE(44), NMAX, IT
          real(KINDDF)::FACTR, PGTOL, DSAVE(29), MXDE
          integer, dimension(:), allocatable     ::NBD    !the boundary type
          integer, dimension(:), allocatable     ::IWA    !the integer working space
          real(KINDDF),  dimension(:), allocatable::X     !the variables
          real(KINDDF),  dimension(:), allocatable::G     !the gradgient
          real(KINDDF),  dimension(:), allocatable::L     !the low boundary of X
          real(KINDDF),  dimension(:), allocatable::U     !the upper boundary of X
          real(KINDDF),  dimension(:), allocatable::WA    !the working space
           real(KINDDF), dimension(:), allocatable::EPOT0   !the working space
!---

                FACTR    = CtrlParam%LBFGS_Factr
                PGTOL    = CtrlParam%LBFGS_PGtol
                MSAVE    = CtrlParam%LBFGS_MSave
                IPRINT   = -1

                IFLAG    = 0
                NPRT     = dm_NPRT
                IFPD     = CtrlParam%IFPD
                LB       = SimBox(1)%BOXLOW
                UB       = SimBox(1)%BOXUP
                BS       = SimBox(1)%ZL
                NMAX = NPRT*3
                allocate(IDX(NPRT), IDY(NPRT), IDZ(NPRT), EPOT0(NPRT), NBD(NMAX), IWA(3*NMAX), X(NMAX), L(NMAX), U(NMAX), G(NMAX),&
                         WA(2*MSAVE*NMAX + 5*nmax + 11*MSAVE*NMAX + 8*MSAVE),STAT=IERR )

                if(IERR .gt. 0) then
                    write(*,fmt="(A)") " MDPSCU Error: fail to allocate working space in DO_LBFGS_FORSTEPS_DEV"
                    write(*,fmt="(A)") "               process to be stopped"
                    stop
                end if
                call CopyXPFrom_Devices_to_Host()
                call CopyStatuFrom_Devices_to_Host()

                NIDX = 0
                do IA = 1, NPRT
                   if(iand(hm_STATU(IA), CP_STATU_FIXPOSX) .eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDX = NIDX + 1
                      IDX(NIDX) = IA
                   end if
                end do

                NIDY = 0
                do IA = 1, NPRT
                   if(iand(hm_STATU(IA), CP_STATU_FIXPOSY) .eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDY = NIDY + 1
                      IDY(NIDY) = IA
                   end if
                end do

                NIDZ = 0
                do IA = 1, NPRT
                   if(iand(hm_STATU(IA), CP_STATU_FIXPOSZ) .eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_OUTOFBOX).eq. 0 .and. &
                      iand(hm_STATU(IA), CP_STATU_ACTIVE)  .eq. CP_STATU_ACTIVE) then
                      NIDZ = NIDZ + 1
                      IDZ(NIDZ) = IA
                   end if
                end do

                X(1:NIDX)                      = hm_XP(IDX(1:NIDX),1)
                X(NIDX+1:      NIDX+NIDY)      = hm_XP(IDY(1:NIDY),2)
                X(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = hm_XP(IDZ(1:NIDZ),3)
                NPRTP = NIDX + NIDY + NIDZ
                !$$--- set the boundary
                !$$    we constrain the variables changes in curoff ranges
                NBD(1:NPRTP)                   = 0
                L(1:NIDX)                      = X(1:NIDX)                      - BS(1)*C_HALF
                U(1:NIDX)                      = X(1:NIDX)                      + BS(1)*C_HALF
                L(NIDX+1:      NIDX+NIDY)      = X(NIDX+1:      NIDX+NIDY)      - BS(2)*C_HALF
                U(NIDX+1:      NIDX+NIDY)      = X(NIDX+1:      NIDX+NIDY)      + BS(2)*C_HALF
                L(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = X(NIDX+NIDY+1: NIDX+NIDY+NIDZ) - BS(3)*C_HALF
                U(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = X(NIDX+NIDY+1: NIDX+NIDY+NIDZ) + BS(3)*C_HALF
                TASK = 'START'
                do I=1, MXNUMSTEPS
                   !--- call to the L-BFGS-B code
                   call SETULB(NPRTP,MSAVE,X,L,U,NBD,OBJF,G,FACTR,PGTOL,WA,IWA,TASK,IPRINT,CSAVE,LSAVE,ISAVE,DSAVE)

                     !----
                      if(IFPD(1)) then
                         do IA=1, NIDX
                            if(X(IA) .lt. LB(1) )then
                               hm_XP(IDX(IA), 1) = X(IA) + BS(1)
                            else if(X(IA) .gt. UB(1) )then
                               hm_XP(IDX(IA), 1) = X(IA) - BS(1)
                            else
                               hm_XP(IDX(IA), 1) = X(IA)
                            end if
                         end do
                      else
                         hm_XP(IDX(1:NIDX), 1) = X(1:NIDX)
                      end if

                      if(IFPD(2)) then
                         do IA=NIDX+1, NIDX+NIDY
                            if(X(IA) .lt. LB(2) )then
                               hm_XP(IDY(IA-NIDX), 2) = X(IA) + BS(2)
                            else if(X(IA) .gt. UB(2) )then
                               hm_XP(IDY(IA-NIDX), 2) = X(IA) - BS(2)
                            else
                               hm_XP(IDY(IA-NIDX), 2) = X(IA)
                            end if
                         end do
                      else
                         hm_XP(IDY(1:NIDY), 2) = X(NIDX+1:NIDX+NIDY)
                      end if

                      if(IFPD(3)) then
                         do IA=NIDX+NIDY+1, NIDX+NIDY+NIDZ
                            if(X(IA) .lt. LB(3) )then
                               hm_XP(IDZ(IA-(NIDX+NIDY)), 3) = X(IA) + BS(3)
                            else if(X(IA) .gt. UB(3) )then
                               hm_XP(IDZ(IA-(NIDX+NIDY)), 3) = X(IA) - BS(3)
                            else
                               hm_XP(IDZ(IA-(NIDX+NIDY)), 3) = X(IA)
                            end if
                         end do
                      else
                         hm_XP(IDZ(1:NIDZ),3) = X(NIDX+NIDY+1: NIDX+NIDY+NIDZ)
                      end if
                      !$$--- reset the position on devices
                      call CopyXPFrom_Host_to_Devices()

                      if(TASK(1:2) .eq. 'FG') then
                      !$$ --- to calculate potential
                      !$$ --- to calculate potential
                         call CalForce_ForceClass(SimBox, CtrlParam, ForceClass)
                         call UpdateEpot_ForceClass(SimBox, CtrlParam, ForceClass)
                         call CopyFPFrom_Devices_to_Host()
                         call CopyEPOTFrom_Devices_to_Host()
                         call SynchronizeDevices()
                         OBJF = sum(hm_EPOT(1:NPRT))
                         !--- NOTE: the following convergence condition do not
                         !          have effects
                         !if(I .eq. 1) then
                         !   EPOT0(1:NPRT) = hm_EPOT(1:NPRT)
                         !else
                         !   MXDE = maxval(dabs(hm_EPOT(1:NPRT) - EPOT0(1:NPRT)))
                         !   if(MXDE*CP_ERGEV .lt. 0.0001D0) then
                         !      write(*,fmt="(A, I6)") " MDPSCU Message: LBFG successful performed for steps", I
                         !      exit
                         !   end if
                         !end if
                         G(1:NIDX)                      = -hm_FP(IDX(1:NIDX),1)
                         G(NIDX+1:      NIDX+NIDY)      = -hm_FP(IDY(1:NIDY),2)
                         G(NIDX+NIDY+1: NIDX+NIDY+NIDZ) = -hm_FP(IDZ(1:NIDZ),3)
                         cycle
                      end if  !--- end if TASK = "FG"

                      if(task(1:5) .eq. 'NEW_X') then
                          cycle
                      end if

                   write(*,fmt="(A, I6)") " MDPSCU Message: LBFG performed for steps ", I
                   write(*,fmt="(A, I6)") "                 exit with code> "// TASK(1:len_trim(TASK))
                   exit
                end do
                deallocate(IDX, IDY, IDZ, EPOT0, NBD, IWA, X, L, U, G, WA,STAT=IERR)
                hm_XP1 = 0.D0
                call CopyXP1From_Host_to_Devices()

                if(I .gt. MXNUMSTEPS) then
                   IFLAG = 1
                end if

                if(TASK(1:5) .eq. 'ERROR' ) then
                   IFLAG = -1
                end if

                if(IFLAG .lt. 0) then
                   write(*,fmt="(A, I2)") " MDPSCU Warning: LBFG exit with code ", IFLAG
                   write(*,fmt="(A, I2)") "                 optimizing may be fail"
                   call ONWARNING(gm_OnWarning)
                else if(I .gt. MXNUMSTEPS) then !if(IFLAG .gt. 0) then
                   write(*,fmt="(A, I2, A, I7, A)") " MDPSCU Warning: LBFG exit with code ", IFLAG, " after ", MXNUMSTEPS, " steps"
                   write(*,fmt="(A, I2)")           "                 more steps may be needed for better convergence"
                   call ONWARNING(gm_OnWarning)
                end if

      return
  end subroutine DO_LBFGSB_FORSTEPS_DEV
  !****************************************************************************************

  end module MD_LBFGSScheme_GPU
