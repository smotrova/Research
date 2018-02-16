! Vichislenie determinanta matrici A poradka NxN

FUNCTION det(A, N) RESULT (f) 

USE DFIMSL 
IMPLICIT NONE
					
COMPLEX(8)									:: f
INTEGER(4), INTENT(IN)						:: N
COMPLEX(8), DIMENSION (0:N, 0:N), INTENT(IN)	:: A
					
INTEGER(4)									:: Lda, Ldfac
COMPLEX(8), DIMENSION (0:N, 0:N)				:: aa
INTEGER(4), DIMENSION (0:N)					:: ipvt
					
COMPLEX(8)									:: det1
REAL(8)										:: det2

!---------------------------------------------------------

	Lda=N+1
	Ldfac=N+1
				
	CALL DLFTCG (N+1, A, Lda, aa, Ldfac, ipvt)
	CALL DLFDCG (N+1, aa, LDFAC, IPVT, DET1, DET2)

	f= det1*10**det2		
					
END FUNCTION det	

!================================================================================================

FUNCTION f1( kapa, gama, N) RESULT (u)

USE DFIMSL 
IMPLICIT NONE

REAL(8), INTENT(IN)					:: kapa, gama
INTEGER(4), INTENT(IN)				:: N
REAL(8)								:: u
COMPLEX(8), DIMENSION (0:N, 0:N)	:: A
COMPLEX(8), EXTERNAL				:: det


	CALL get_matrix (kapa, gama, A, N) 			
			
	u=real(det(A, N)) 

END FUNCTION f1

!================================================================================================

FUNCTION f2( kapa, gama, N) RESULT (u)

USE DFIMSL 
IMPLICIT NONE

REAL(8), INTENT(IN)					:: kapa, gama
INTEGER(4), INTENT(IN)				:: N
REAL(8)								:: u
COMPLEX(8), DIMENSION (0:N, 0:N)	:: A
COMPLEX(8), EXTERNAL				:: det


	CALL get_matrix (kapa, gama, A, N) 
						
	u=aimag(det(A, N)) 

END FUNCTION f2

!================================================================================================

SUBROUTINE FCN (X, F, NN)

USE DFIMSL 
USE init_data_class

IMPLICIT NONE

INTEGER(4), INTENT(IN)	:: NN
REAL(8), INTENT(IN)		:: X(NN)
REAL(8), INTENT(OUT)	:: F(NN)

REAL(8), EXTERNAL		:: f1, f2

!--------------------------------------

	F(1) = f1( x(1), x(2), N)
	F(2) = f2( x(1), x(2), N)

	RETURN

END

!================================================================================================

SUBROUTINE LSJAC (NN, X, FJAC)

USE DFIMSL 
USE init_data_class

IMPLICIT NONE

INTEGER(4), INTENT(IN)	:: NN

REAL(8), INTENT(IN)		:: X(NN)

REAL(8), INTENT(OUT)	:: FJAC(NN, NN)

REAL(8), DIMENSION(NN)  :: h
REAL(8), EXTERNAL		:: f1, f2

REAL(8)					:: kapa1, kapa2, gama1, gama2

!--------------------------------------
	
	h = (/1.0D-8, 1.0D-12/)

	gama1 = x(2)+h(2)
	gama2 = x(2)-h(2)

	kapa1 = x(1)+h(1)
	kapa2 = x(1)-h(1)

	FJAC(1,1) = (f1(kapa1, x(2), N)-f1(kapa2, x(2), N))/h(1)/2
	FJAC(1,2) = (f1(x(1), gama1, N)-f1(x(1), gama2, N))/h(2)/2
	FJAC(2,1) = (f2(kapa1, x(2), N)-f2(kapa2, x(2), N))/h(1)/2
	FJAC(2,2) = (f2(x(1), gama1, N)-f2(x(1), gama2, N))/h(2)/2

	RETURN

END

!================================================================================================