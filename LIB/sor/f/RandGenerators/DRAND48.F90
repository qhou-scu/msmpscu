module RAND48_MODULE
!*** This module contains routines to generating random numbers
!*** NOTE: DRAND48 is supported only by 64bit compiler
!         
     implicit none
contains

!********************************************************************************
      subroutine DRAND48_PUTSEED(seedval)
      !*** To initialize the seed for DRAND48
      !
      !    See also: DRAND48()  
	  !
      !	   Author: HOU Qing, 2009-12
      !
      !    INPUT:seedval - the seed
      !
           !--- COMMON BLOCK
           integer(2)::DRAND48_SEED(3)=(/12345,6778,678/)
           common/RAND48_SEED/DRAND48_SEED
           !---DUMMY VARAIBLES:													   
		   integer(2),intent(in)::seedval(3)
           !--- LOCAL VARAIBLES

		 
           DRAND48_SEED(3) = seedval(3) !ishft(seedval,-16)
           DRAND48_SEED(2) = seedval(2) !iand(seedval, 65535)
           DRAND48_SEED(1) = seedval(1) !13070
         return
      end subroutine DRAND48_PUTSEED
!********************************************************************************

!********************************************************************************
      subroutine DRAND48_GETSEED(seedval)
      !*** To get the current seed for DRAND48
      !
      !    See also: DRAND48()  
	  !
      !	   Author: HOU Qing, 2009-12
      !
      !    INPUT:seedval - the seed
      !
           !--- COMMON BLOCK
           integer(2)::DRAND48_SEED(3)
           common/RAND48_SEED/DRAND48_SEED
           !---DUMMY VARAIBLES:													   
		   integer(2),intent(out)::seedval(3)
           !--- LOCAL VARAIBLES

		 
           seedval(3) = DRAND48_SEED(3) !ishft(seedval,-16)
           seedval(2) = DRAND48_SEED(2) !iand(seedval, 65535)
           seedval(1) = DRAND48_SEED(1) 
         return
      end subroutine DRAND48_GETSEED
!********************************************************************************

      
!********************************************************************************
      real(8) function DRAND48() result(R)
      !*** To generate a random number
      !    This is a FORTRAN version adopted from DRAND48 originally written in C.
      !    Ref:   DRAND.TXT, and related Web page
      !
      !    See also:   DRAND48_PUTSEED  
      !
      !	   Author: HOU Qing, 2009-12
      !    INPUT:no
      !
	       !integer(2),parameter::DRAND48_A(3)=(/58989,57068,5/), 
           integer(8), parameter::T16=65536, T32=4294967296, T48=281474976710656, DRAND48_A=25214903917, DRAND48_C=11
           !--- COMMON BLOCK
           integer(2)::DRAND48_SEED(3)
           common/RAND48_SEED/DRAND48_SEED
           !--- LOCAL VARAIBLES
           integer(8)::X, T
           integer(4)::T1=0, T2, T3;
           integer(2)::XBIT(4)=0, tBIT(4)=0, t1BIT(2), t2BIT(2), t3BIT(2)
           equivalence(X,XBIT)
           equivalence(T,tBIT)
           equivalence(T1, T1BIT)
           equivalence(T2, T2BIT)
           equivalence(T3, T3BIT)
           !equivalence(a,aBIT)
			XBIT(3) = DRAND48_SEED(3)
			XBIT(2) = DRAND48_SEED(2)
			XBIT(1) = DRAND48_SEED(1) 

            T =  X*DRAND48_A + DRAND48_C

            T3=0
            T3 = IOR(T3,tBIT(3))
            if(T3 .lt. 0) T3BIT(2) = 0
            
            T2=0
            T2 = IOR(T2,tBIT(2))
            if(T2 .lt. 0) T2BIT(2) = 0
            
            T1 = 0
            T1 = IOR(T1,tBIT(1)) 
            if(T1 .lt. 0)  T1BIT(2) = 0
            
            R = dble(T3)/dble(T48)+dble(T2)/dble(T32)+dble(T1)/dble(T16)
            
            DRAND48_SEED(1) = tBIT(1)
            DRAND48_SEED(2) = tBIT(2)
            DRAND48_SEED(3) = tBIT(3)
          
          return
      end function DRAND48
!********************************************************************************
      
!********************************************************************************
      real(8) function DRAND48_S(STEP,RESETSEED) result(R)
      !*** To generate the STEPth random number after the current SEED
      !
      !
      !	   Author: HOU Qing, 2009-12
      !
      !    INPUT: STEP, the step size of genecrating the next random number
      !                 for example:
      !                    if step = 1, DRAND48_S generate the same number as DRAND48
      !				       if step = 2, DRAND48_S equivalent to call function DRAND48 twice
      !                    ....
      !          RESET, = 0,   leave the SEED as it is, 
      !                 = 1,   reset the SEED 
           integer(8), parameter::T16=65536, T32=4294967296, T48=281474976710656, HT64=9223372036854775808, DRAND48_A=25214903917,DRAND48_C=11
           !--- COMMON BLOCK
           integer(2)::DRAND48_SEED(3)
           common/RAND48_SEED/DRAND48_SEED
           !--- DUMMY VARIABLES
           integer(4),intent(in)::STEP, RESETSEED
             
           !--- LOCAL VARAIBLES
           integer(8)::X, T
           integer(2)::XBIT(4)=0, tBIT(4)=0, t1BIT(2), t2BIT(2), t3BIT(2)
           integer(4)::T1,T2,T3
           equivalence(T1, T1BIT)
           equivalence(T2, T2BIT)
           equivalence(T3, T3BIT)
           equivalence(X,XBIT)
           equivalence(T,tBIT)
           integer(8)::BIGRAND48_A=DRAND48_A,BIGRAND48_C=DRAND48_C;
           integer::OLD_STEP = 1, I
           save OLD_STEP,BIGRAND48_A, BIGRAND48_C 
            
			XBIT(3) = DRAND48_SEED(3)
			XBIT(2) = DRAND48_SEED(2)
			XBIT(1) = DRAND48_SEED(1) 
            
            if(OLD_STEP .NE. STEP) THEN
               BIGRAND48_A = 1
               BIGRAND48_C = 0
               DO I=1, STEP
                  BIGRAND48_A = BIGRAND48_A*DRAND48_A
                  BIGRAND48_C = 1+BIGRAND48_C*DRAND48_A
               END DO
               BIGRAND48_C = BIGRAND48_C*DRAND48_C
               OLD_STEP = STEP
            end if   
               
            T =  X*BIGRAND48_A + BIGRAND48_C

            T3=0
            T3 = IOR(T3,tBIT(3))
            if(T3 .lt. 0) T3BIT(2) = 0
            
            T2=0
            T2 = IOR(T2,tBIT(2))
            if(T2 .lt. 0) T2BIT(2) = 0
            
            T1 = 0
            T1 = IOR(T1,tBIT(1)) 
            if(T1 .lt. 0)  T1BIT(2) = 0
                        
            R = dble(T3)/dble(T48)+dble(T2)/dble(T32)+dble(T1)/dble(T16)
            
            if(RESETSEED .GT. 0) then
               DRAND48_SEED(1) = tBIT(1)
               DRAND48_SEED(2) = tBIT(2)
               DRAND48_SEED(3) = tBIT(3)
            end if   
          return
      end function DRAND48_S
!********************************************************************************

!********************************************************************************
      subroutine DRAND48_P(N, RAN) 
      !*** To generate a array of random numbers from a previous set of random-numbers
      !    The size of the generated array is the same as the input seed.
      !    NOTE: the SEEDS should be a set of random-number continuesly generated. 
      !
      !
      !	   Author: HOU Qing, 2009-12
      !
      !    INPUT: N,     the size of random number array
      ! 
      !    OUTPUT RAN,   the generated random number array
           integer(8), parameter::T16=65536, T32=4294967296, T48=281474976710656, DRAND48_A=25214903917,DRAND48_C=11
           !--- COMMON BLOCK
           integer(2)::DRAND48_SEED(3)
           common/RAND48_SEED/DRAND48_SEED
           !--- DUMMY VARIABLES
           integer(4),intent(in)::N
           real(8)::RAN(N)   
           !--- LOCAL VARAIBLES
           integer(8)::X, T, T4(3)
           integer(2)::XBIT(4)=0, tBIT(4)=0
           equivalence(X,XBIT)
           equivalence(T,tBIT)
           integer(8)::BIGRAND48_A=DRAND48_A,BIGRAND48_C=1,BIGRAND48_C0;
           integer::OLD_STEP = 1, I
           save OLD_STEP,BIGRAND48_A, BIGRAND48_C 
            
			XBIT(3) = DRAND48_SEED(3)
			XBIT(2) = DRAND48_SEED(2)
			XBIT(1) = DRAND48_SEED(1) 
			
            BIGRAND48_A  = DRAND48_A
            BIGRAND48_C0 = 1
			DO I=1, N
			   BIGRAND48_C = BIGRAND48_C0*DRAND48_C
               T =  X*BIGRAND48_A + BIGRAND48_C

               T4 = 0
               T4(3) = IOR(T4(3),tBIT(3))
               if(T4(3) .lt. 0) T4(3) = T4(3)+32768
               
               T4(2) = IOR(T4(2),tBIT(2))
               if(T4(2) .lt. 0) T4(2) = T4(2)+32768
               
               T4(1) = IOR(T4(1),tBIT(1))
               if(T4(1) .lt. 0) T4(1) = T4(1)+32768
               
               
               RAN(I) = dble(T4(3))/dble(T48)+dble(T4(2))/dble(T32)+dble(T4(1))/dble(T16)

               !Reset the A value and C value for the next random number                
               BIGRAND48_A = BIGRAND48_A*DRAND48_A
               BIGRAND48_C0 = 1+BIGRAND48_C0*DRAND48_A
            END DO   
            
            !Reset the seed
            DRAND48_SEED(1) = tBIT(1)
            DRAND48_SEED(2) = tBIT(2)
            DRAND48_SEED(3) = tBIT(3)
          return
      end subroutine DRAND48_P
!********************************************************************************
      
!********************************************************************************
!     UNCOMPLETED !!!!!!!!!
      subroutine IRAND_BSKS_PUTSEED(seedval, N, SARRAY)
      !*** To set the initial seed for the random number generator IRAND_BSKS
      !    Ref: Computational Physica, edited by K.H.Hoffmann and M.Schreiber, Springer-Verlag Beilin Heideberg 1996
      !         计算物理， K.H.Hoffmann and M.Schreiber 科学出版社（2001）
      !
      !	   Author: HOU Qing, 2009-12
      !
      !    INPUT: seedva - integer, the seed 
      !           N      - the size of array N0
      !           SARRAY -   the array storing the initial random ranber array
      
           integer(4), parameter::MARICEL=16807
           !--- COMMON BLOCK
           !--- DUMMY VARIABLES
           integer(4),intent(in)::seedval, N
           integer(4)::SARRAY(N)
             
           !--- LOCAL VARAIBLES
           integer(4)::I,K,SEED,ICL
               
               SEED=2*seedval+1
               DO K=1, N
                  ICL = 0
                  DO I=1,32
                  ENDDO              
			   END DO
		  return
	  end subroutine IRAND_BSKS_PUTSEED
!********************************************************************************
	  

end module RAND48_MODULE