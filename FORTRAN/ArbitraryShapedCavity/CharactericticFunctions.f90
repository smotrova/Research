! Vichislenie determinanta matrici A poradka NxN

FUNCTION det(A, N) RESULT (f) 

USE DFIMSL 
IMPLICIT NONE
					
COMPLEX(8)									:: f
INTEGER(4), INTENT(IN)						:: N

COMPLEX(8), DIMENSION (0:4*N-1, 0:4*N-1)	:: A
					
INTEGER(4)									:: Lda, Ldfac

COMPLEX(8), DIMENSION (0:4*N-1, 0:4*N-1)	:: aa

INTEGER(4), DIMENSION (0:4*N-1)					:: ipvt
					
COMPLEX(8)									:: det1
REAL(8)										:: det2

!---------------------------------------------------------

	Lda = 4*N
	Ldfac = 4*N
				
	CALL DLFTCG (4*N, A, Lda, aa, Ldfac, ipvt)
	CALL DLFDCG (4*N, aa, LDFAC, IPVT, DET1, DET2)

	f= det1*10**det2		
					
END FUNCTION det	

!================================================================================================

FUNCTION f1( kapa, gama, N) RESULT (u)

USE DFIMSL 
USE init_data_class, ONLY : PI, alfa_i, alfa_e, im

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama
INTEGER(4), INTENT(IN)						:: N
REAL(8)										:: u

COMPLEX(8), DIMENSION (0:4*N-1, 0:4*N-1)	:: A

COMPLEX(8), EXTERNAL						:: det



!-----------------------------------------------
	
	CALL get_matrix (kapa, gama, A, N) 	

	u = DREAL(det(A, N)) 


END FUNCTION f1

!================================================================================================

FUNCTION f2( kapa, gama, N) RESULT (u)

USE DFIMSL 
USE init_data_class, ONLY : PI, alfa_i, alfa_e, im

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama
INTEGER(4), INTENT(IN)						:: N
REAL(8)										:: u

COMPLEX(8), DIMENSION (0:4*N-1, 0:4*N-1)	:: A

COMPLEX(8), EXTERNAL						:: det


!----------------------------------------------------------


	CALL get_matrix (kapa, gama, A, N) 

	u = aimag(det(A, N)) 

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

!=================================================================================