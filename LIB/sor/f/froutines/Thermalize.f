      subroutine Thermalize(TI, NPRT, NG, CM,ITYP,VX,VY,VZ)                                                         
      !ms$attributes dllexport :: THERMALIZE
!***  PURPOSE:   to thermalize 
!     INPUT:  TI, temperature   
!             NPRT, the number 
!             NGROUP, number of groups
!             CM, atomic mass of 
!             ITYP, ID of atoms

!     OUTPUT  VX, VY, VZ, the velocity obtained
      IMPLICIT NONE                                       
	real*4::TI,CM(*),Vx(*),Vy(*),Vz(*)
	integer::NPRT, NG, ITYP(*)
	!Loacal variables
	real*8::VTx, VTy, VTz, WT, Z1, Z2
	real*8, dimension(:), pointer::V0, RVPx, RVPy, RVPz, CMR
      integer::I,J
      real*8::KB = 1.38054D-16 
      real*8::TWOPI = 6.283185308D0                                                 
      real*8::AU = 1.66043D-24

	allocate(V0(nG),CMR(nG),RVPx(NPRT),RVPy(NPRT),RVPz(NPRT))
	DO I=1, nG
	   CMR(I) = CM(I)*AU
      END DO
      v0  = dsqrt(2.D0*TI*kb/CMR)                                          
      VTx = 0.D0
      VTy = 0.D0
      VTz = 0.D0
	RVPx = 0.D0
	RVPy = 0.D0
	RVPz = 0.D0
	
      WT = 0.D0
	DO I=1, NPRT                                                               
	   J = ITYP(I)
         WT = WT + CMR(J)
      END DO	                                       

      DO I=1, NPRT
	   J = ITYP(I)
         call RANDOM_NUMBER(Z1)	                                       
         call RANDOM_NUMBER(Z2)	                                       
         RVPx(I)=V0(J)*SQRT(-LOG(Z1))*COS(TWOPI*Z2)                          
         VTx =VTx + (CM(J)*RVPx(I))/WT
         call RANDOM_NUMBER(Z1)	                                       
         call RANDOM_NUMBER(Z2)	                                       
         RVPy(I)=V0(J)*SQRT(-LOG(Z1))*COS(TWOPI*Z2)                          
         VTy =VTy + (CM(J)*RVPy(I))/WT                                    
         call RANDOM_NUMBER(Z1)	                                       
         call RANDOM_NUMBER(Z2)	                                       
         RVPz(I)=V0(J)*SQRT(-LOG(Z1))*COS(TWOPI*Z2)                          
         VTz =VTz + (CM(J)*RVPz(I))/WT                                    
      END DO

      
	DO I=1, NPRT 
         RVPX(I) = (RVPx(I)-VTx)                                 
         RVPY(I) = (RVPy(I)-VTy)                                 
         RVPZ(I) = (RVPz(I)-VTz)                                 
      END DO 

	VTx = 0.D0
	VTy = 0.D0
	VTz = 0.D0
	WT  = 0.D0
	DO I=1, nG
	   WT = WT +CMR(I)
      END DO
	DO I=1, nG
	   CMR(I) = CMR(I)/WT
      END DO

	DO I=1, NPRT 
	   J = ITYP(I)
         VTx = VTx + VX(I)*CMR(J)                                 
         VTy = VTy + VY(I)*CMR(J)                                 
         VTz = VTz + VZ(I)*CMR(J)                                 
      END DO 

	DO I=1, NPRT 
         Vx(I) = RVPx(I) + VTx                                 
         Vy(I) = RVPy(I) + VTy                                 
         Vz(I) = RVPz(I) + VTz                                 
      END DO 

	deallocate(V0,CMR,RVPx,RVPy, RVPz)
      RETURN                                                                    
      end subroutine Thermalize

      subroutine Accelerate(VC, NPRT, NG, CM,ITYP,VX,VY,VZ)                                                         
      !ms$attributes dllexport :: ACCELERATE
!***  PURPOSE:   to set the cluster with the velocity of center of mass to be V 
!     INPUT:  TI, temperature   
!             NPRT, the number 
!             NGROUP, number of groups
!             CM, atomic mass of 
!             ITYP, ID of atoms

!     OUTPUT  VX, VY, VZ, the velocity obtained
      IMPLICIT NONE                                       
	real*4::VC(*),CM(*),Vx(*),Vy(*),Vz(*)
	integer::NPRT, NG, ITYP(*)
	!Loacal variables
	real*8::VTx, VTy, VTz, WT
	real*8, dimension(:), pointer::CMR
	integer::I,J

	allocate(CMR(nG))
	DO I=1, nG
	   CMR(I) = CM(I)
      END DO

	WT  = 0.D0
	DO I=1, nG
	   WT = WT +CMR(I)
      END DO
	DO I=1, nG
	   CMR(I) = CMR(I)/WT
      END DO


	VTx = 0.D0
	VTy = 0.D0
	VTz = 0.D0

      !*** the old velocity for center of mass
	DO I=1, NPRT 
	   J = ITYP(I)
         VTx = VTx + VX(I)*CMR(J)                                 
         VTy = VTy + VY(I)*CMR(J)                                 
         VTz = VTz + VZ(I)*CMR(J)                                 
      END DO 

	DO I=1, NPRT 
         Vx(I) = Vx(I) + (Vc(1) - VTx )                                
         Vy(I) = Vy(I) + (Vc(2) - VTy )                                 
         Vz(I) = Vz(I) + (Vc(3) - VTz )                                
      END DO 

	deallocate(CMR)
      RETURN                                                                    
      end subroutine Accelerate
