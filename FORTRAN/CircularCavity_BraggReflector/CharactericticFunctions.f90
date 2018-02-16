! Vichislenie determinanta matrici A poradka NxN

FUNCTION det(A, N) RESULT (f) 

USE DFIMSL 
IMPLICIT NONE
					
COMPLEX(8)									:: f
INTEGER(4), INTENT(IN)						:: N
COMPLEX(8), DIMENSION (1:N, 1:N), INTENT(IN)	:: A
					
INTEGER(4)									:: Lda, Ldfac
COMPLEX(8), DIMENSION (1:N, 1:N)			:: aa
INTEGER(4), DIMENSION (1:N)					:: ipvt
					
COMPLEX(8)									:: det1
REAL(8)										:: det2

!---------------------------------------------------------

	Lda=N
	Ldfac=N
				
	CALL DLFTCG (N, A, Lda, aa, Ldfac, ipvt)

	CALL DLFDCG (N, aa, LDFAC, IPVT, DET1, DET2)

	f= det1*10**det2		
					
END FUNCTION det	

!================================================================================================

FUNCTION f1( kapa, gama) RESULT (u)

USE DFIMSL 
USE init_data_class, ONLY : M


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama
REAL(8)									:: u
COMPLEX(8), DIMENSION (1:2*M, 1:2*M)	:: A
COMPLEX(8), EXTERNAL					:: det
	

	CALL get_matrix (kapa, gama, A, M) 
			
	u=real(det(A, 2*M)) 

END FUNCTION f1

!================================================================================================

FUNCTION f2( kapa, gama) RESULT (u)

USE DFIMSL 
USE init_data_class, ONLY : M

IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama
REAL(8)									:: u
COMPLEX(8), DIMENSION (1:2*M, 1:2*M)		:: A
COMPLEX(8), EXTERNAL					:: det


	CALL get_matrix (kapa, gama, A, M) 
						
	u=aimag(det(A, 2*M)) 

END FUNCTION f2

!================================================================================================

SUBROUTINE FCN (X, F, NN)

USE DFIMSL 
USE init_data_class, ONLY : M


IMPLICIT NONE

INTEGER(4), INTENT(IN)	:: NN
REAL(8), INTENT(IN)		:: X(NN)
REAL(8), INTENT(OUT)	:: F(NN)

REAL(8), EXTERNAL		:: f1, f2

!--------------------------------------

	F(1) = f1( x(1), x(2))
	F(2) = f2( x(1), x(2))

	RETURN

END

!================================================================================================

SUBROUTINE LSJAC (NN, X, FJAC)

USE DFIMSL 

IMPLICIT NONE

INTEGER(4), INTENT(IN)	:: NN

REAL(8), INTENT(IN)		:: X(NN)

REAL(8), INTENT(OUT)	:: FJAC(NN, NN)

REAL(8), DIMENSION(NN)  :: h
REAL(8), EXTERNAL		:: f1, f2

REAL(8)					:: kapa1, kapa2, gama1, gama2, kapa00, gama00

!--------------------------------------
	
	h = (/1.0D-5, 1.0D-5/)

	kapa00 = x(1) 
	gama00 = x(2)

	FJAC(1,1) = ( -f1(kapa00+2*h(1), gama00) + 8*f1(kapa00+h(1), gama00)-8*f1(kapa00-h(1), gama00) + f1(kapa00-2*h(1), gama00) )/12/h(1)
	FJAC(1,2) = ( -f1(kapa00, gama00+2*h(2)) + 8*f1(kapa00, gama00+h(2))-8*f1(kapa00, gama00-h(2)) + f1(kapa00, gama00-2*h(2)) )/12/h(2)
	FJAC(2,1) = ( -f2(kapa00+2*h(1), gama00) + 8*f2(kapa00+h(1), gama00)-8*f2(kapa00-h(1), gama00) + f2(kapa00-2*h(1), gama00) )/12/h(1)
	FJAC(2,2) = ( -f2(kapa00, gama00+2*h(2)) + 8*f2(kapa00, gama00+h(2))-8*f2(kapa00, gama00-h(2)) + f2(kapa00, gama00-2*h(2)) )/12/h(2)



	RETURN

END

!================================================================================================

