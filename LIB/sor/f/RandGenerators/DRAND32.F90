module RAND32_MODULE
!*** This module contains routines to generate random numbers
!*** Adtopted from the subroutine used originally in PENELOPE
!    adpoting the algorithm from  (Comput. Phys. Commun. 60 (1990) 329-344)
!         
     implicit none
     integer(4), private::DRAND32_SEED1=312312,DRAND32_SEED2=31313
     real(8),   private, PARAMETER::USCALE=1.0D0/2.0D0**31
contains

!********************************************************************************
      subroutine DRAND32_PUTSEED(seedval)
      !*** To initialize the seed for DRAND32
      !
      !    See also: DRAND48()  
	  !
      !	   Author: HOU Qing, 2010-9
      !
      !    INPUT:seedval - the seed
      !
           !--- COMMON BLOCK
           
           !---DUMMY VARAIBLES:													   
		   integer(4)::seedval(2)
           !--- LOCAL VARAIBLES

		 
           DRAND32_SEED1 = seedval(1)
           DRAND32_SEED2 = seedval(2) 
         return
      end subroutine DRAND32_PUTSEED
!********************************************************************************

!********************************************************************************
      subroutine DRAND32_GETSEED(seedval)
      !*** To get the current seed for DRAND32
      !
      !    See also: DRAND48()  
	  !
      !	   Author: HOU Qing, 2010-9
      !
      !    OUTPUT:seedval - the seed
      !
           !--- COMMON BLOCK
           !---DUMMY VARAIBLES:													   
		   integer(4),intent(out)::seedval(2)
           !--- LOCAL VARAIBLES

		 
           seedval(1) = DRAND32_SEED1
           seedval(2) = DRAND32_SEED2
           
         return
      end subroutine DRAND32_GETSEED
!********************************************************************************

      
!********************************************************************************
      real(8) function DRAND32() result(R)
      !*** To generate a random number
      !    This is a FORTRAN version adopted from DRAND48 originally written in C.
      !    Ref:   DRAND.TXT, and related Web page
      !
      !    See also:   DRAND48_PUTSEED  
      !
      !	   Author: HOU Qing, 2010-9
      !    INPUT:no
      !

       integer(4)::I1,I2,IZ
        
         I1=DRAND32_SEED1/53668
         DRAND32_SEED1=40014*(DRAND32_SEED1-I1*53668)-I1*12211
         IF(DRAND32_SEED1.LT.0) DRAND32_SEED1=DRAND32_SEED1+2147483563

         I2=DRAND32_SEED2/52774
         DRAND32_SEED2=40692*(DRAND32_SEED2-I2*52774)-I2*3791
         IF(DRAND32_SEED2.LT.0) DRAND32_SEED2=DRAND32_SEED2+2147483399

         IZ=DRAND32_SEED1-DRAND32_SEED2
         IF(IZ.LT.1) IZ=IZ+2147483562
 
         R=IZ*USCALE  
                 
         return
      end function DRAND32
!********************************************************************************
      
!********************************************************************************
      real(8) function DRAND32_S(STEP,RESETSEED) result(R)
      !*** To generate the STEPth random number after the current SEED
      !
      !
      !	   Author: HOU Qing, 2010-9
      !
      !    INPUT: STEP, the step size of genecrating the next random number
      !                 for example:
      !                    if step = 1, DRAND32_S generate the same number as DRAND32
      !				       if step = 2, DRAND32_S equivalent to call function DRAND32 twice
      !                    ....
      !          RESET, = 0,   leave the SEED as it is, 
      !                 = 1,   reset the SEED 
           
       integer(4),intent(in)::STEP,RESETSEED    
       integer(4)::I, I1,I2, IZ, ISEED1, ISEED2
       integer(4)::SSEED1, SSEED2
       save SSEED1,SSEED2
       
          ISEED1 = DRAND32_SEED1
          ISEED2 = DRAND32_SEED2
          !need to be modified, otherwisw useless
          DO I=1, STEP
        
             I1=ISEED1/53668
             ISEED1=40014*(ISEED1-I1*53668)-I1*12211

             I2=ISEED2/52774
             ISEED2=40692*(ISEED2-I2*52774)-I2*3791
          END DO   

          IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
          IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
          
          IZ=ISEED1-ISEED2
          IF(IZ.LT.1) IZ=IZ+2147483562
          R=IZ*USCALE  
            
          if(RESETSEED .GT. 0) then
             DRAND32_SEED1 = ISEED1
             DRAND32_SEED2 = ISEED2
          end if   
          return
      end function DRAND32_S
!********************************************************************************

	  

end module RAND32_MODULE