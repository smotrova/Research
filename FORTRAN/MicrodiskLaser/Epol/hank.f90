!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------Podprograma dla vichislenia------------------------------------

!-------------------------------funkcii Bessela -----------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!***************************************************************************************


	COMPLEX(8) FUNCTION HNK(n,arg)		!HANKEL FUNCTION of REAL ARGUMENT
	  use  dfimsl

	IMPLICIT NONE
	INTEGER		n
	REAL(8)     arg, xnu
	COMPLEX(8) :: IM=(0.0D0, 1.0D0)
	REAL(8)	 JJ(0:n+1), YY(0:n+1)     
		
			xnu=0.0D0


      	  CALL DBSYS (xnu,arg, N+1, YY)
	  	  CALL DBSJNS (arg, N+1, JJ)
	      HNK = JJ(n) + IM*YY(n)

	END FUNCTION


!********************************************************************************************

COMPLEX(8) FUNCTION Y(n,arg)		!Neyman FUNCTION of REAL ARGUMENT
  use  dfimsl

	IMPLICIT NONE
	INTEGER		n
	REAL(8)     arg, xnu
	REAL(8)	  YY(0:n+1)     
		
			xnu=0.0D0


      	  CALL DBSYS (xnu,arg, N+1, YY)
	  	  
	      Y = YY(n)

	END FUNCTION



!**************************************************************************************


COMPLEX(8) FUNCTION CBESSEL(n,arg)		!BESSEL FUNCTION J of COMPLEX ARGUMENT
	use  dfimsl
	IMPLICIT NONE
	
	INTEGER		n
	COMPLEX(8)	arg
	COMPLEX(8)	JJ(0:n+1)
	
    CALL DCBJNS (arg, N+1, JJ)
    CBESSEL = JJ(n) 

	END FUNCTION


!******************************************************************************************
!******************************************************************************************
!******************************************************************************************


	COMPLEX(8) FUNCTION Fi(n,kapa, gama, alfa)		
	IMPLICIT NONE

	INTEGER		n
	COMPLEX(8)	arg, arg1, CBESSEL, HNK, s1, s2

	REAL(8)		alfa, kapa, gama 
	EXTERNAL	CBESSEL
	EXTERNAL	HNK
	

			arg=DCMPLX(alfa*kapa, gama*kapa)		!arg=kapa*beta
			arg1=DCMPLX(alfa, gama)				!beta=alfa+i*gama

			Fi=CBESSEL(n, arg)*HNK(n-1, kapa) - arg1*CBESSEL(n-1, arg)*HNK(n, kapa)		!E-pol

	END FUNCTION


!******************************************************************************************

	REAL(8) FUNCTION f1(n,kapa, gama, alfa)	 ! sqrt(epsilon)=beta=alfa+i*gama
	IMPLICIT NONE							 ! kapa = k*a, k-volnovoe chislo, a -radius 

	INTEGER		n
	COMPLEX(8)	Fi
	REAL(8)		kapa, gama, alfa
	EXTERNAL	Fi
	 
		f1=DBLE(Fi(n, kapa, gama, alfa))

	END FUNCTION

!****************************************************************************************


	REAL(8) FUNCTION f2(n, kapa, gama, alfa )	! sqrt(epsilon)=beta=alfa+i*gama	
	IMPLICIT NONE

	INTEGER		n
	COMPLEX(8)	Fi
	REAL(8)		kapa, gama, alfa
	EXTERNAL	Fi
	
		f2=AIMAG(Fi(n,kapa, gama, alfa))

	END FUNCTION

!*****************************************************************************************
!--------------------------- Function for Ez field component -----------------------------

	COMPLEX(8) FUNCTION E(n,k, gama, r)		
	IMPLICIT NONE

	INTEGER		n
	COMPLEX(8)	arg, arg1, CBESSEL, HNK

	REAL(8)		alfa, k, gama, r 
	EXTERNAL	CBESSEL
	EXTERNAL	HNK
	COMMON /ALFA/ alfa   
	

			arg=DCMPLX(alfa*k*r, gama*k*r)		!arg=kapa*beta
			arg1=DCMPLX(alfa, gama)				!beta=alfa+i*gama

			IF ( r<=1) THEN
				E=HNK(n, k)/CBESSEL(n, k*arg1)*CBESSEL(n, arg)
			END IF

			IF (r>1 ) THEN
				E=HNK(n, k*r)
			END IF
			
	END FUNCTION

!*****************************************************************************************