!=======================================================================================

!-------------------------Podprograma dla vichislenia-----------------------------------

!-------------------------------funkcii Besselya ---------------------------------------

!=======================================================================================

SUBROUTINE Hankel_func (N, arg, HNK)		!HANKEL FUNCTION of REAL ARGUMENT
	
USE  dfimsl
IMPLICIT NONE

INTEGER(4), INTENT(IN)		:: N
REAL(8), INTENT(IN)			:: arg
COMPLEX(8), INTENT(OUT)		:: HNK(0:N)		
INTEGER(4)					:: i 
REAL(8)						:: xnu
COMPLEX(8)					:: IM=(0.0D0, 1.0D0)
REAL(8)						:: JJ(N+1), YY(N+1)     

!--------------------------------------------------

	xnu=0.0D0
		
	CALL DBSYS (xnu,arg, N+1, YY)
  	CALL DBSJNS (arg, N+1, JJ)

	DO i=0, N
	    
		HNK(i) = JJ(i+1) + IM*YY(i+1)

	END DO

END SUBROUTINE Hankel_func

!==========================================================================================

SUBROUTINE Neyman_func (N, arg, Y) !NEYMAN FUNCTION of REAL ARGUMENT

USE  dfimsl
IMPLICIT NONE

INTEGER(4), INTENT(IN)		:: N
REAL(8), INTENT(IN)			:: arg
REAL(8), INTENT(OUT)		:: Y(0:N)	
INTEGER(4)					:: i 
REAL(8)						:: xnu
REAL(8)						:: YY(N+1)     

!-------------------------------------------------
		
	xnu=0.0D0

	CALL DBSYS (xnu, arg, N+1, YY)

	DO i=0, N
	
		Y(i) = YY(i+1)

	END DO

END SUBROUTINE Neyman_func

!============================================================================================

SUBROUTINE CNeyman_func (N, arg, CY)

USE  dfimsl
IMPLICIT NONE

INTEGER(4), INTENT(IN)		:: N
COMPLEX(8), INTENT(IN)		:: arg
COMPLEX(8), INTENT(OUT)		:: CY(0:N)     

REAL(8)						:: xnu
INTEGER(4)					:: i
COMPLEX(8)					:: YY(N+1)     

!-------------------------------------------------
	
	xnu=0.0D0

   	CALL DCBYS (xnu, arg, N+1, YY)

	DO i=0,N

  		CY(i) = YY(i+1)

	END DO

END SUBROUTINE CNeyman_func

!==============================================================================================

SUBROUTINE Bessel_func (N, arg, J)	!BESSEL FUNCTION J of REAL ARGUMENT

USE  dfimsl
IMPLICIT NONE
	
INTEGER(4), INTENT(IN)		:: N
REAL(8), INTENT(IN)			::	arg
REAL(8), INTENT(OUT)		::	J(0:N)
INTEGER(4)					:: i 
REAL(8)						:: JJ(N+1)

!--------------------------------------------
	
	CALL DBSJNS (arg, N+1, JJ)

	DO i=0,N

		J(i)= JJ(i+1)
		
	END DO	 

END SUBROUTINE Bessel_func

!==============================================================================================

SUBROUTINE CBessel_func (N, arg, CJ)	!BESSEL FUNCTION J of COMPLEX ARGUMENT

USE  dfimsl
IMPLICIT NONE
	
INTEGER(4), INTENT(IN)	::	N
COMPLEX(8), INTENT(IN)	::	arg
COMPLEX(8), INTENT(OUT) ::	CJ(0:N)
COMPLEX(8)				::	JJ(N+1)
INTEGER(4)				::	i

!--------------------------------------------
	
	CALL DCBJNS (arg, N+1, JJ)

	DO i = 0, N

		CJ(i) = JJ(i+1) 

	END DO

END SUBROUTINE CBessel_func

!===============================================================================================