FUNCTION Integral_roZpZp (kapa, gama, r) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama, r

INTEGER(4)								:: i

COMPLEX(8)								:: Z(0:n+1), V(0:n+1)

COMPLEX(8)								:: u, Z1,Z2,V1,V2

COMPLEX(8)								:: nu(1:M+1)

!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


CALL CBessel_func (n+1, nu(1)*kapa*r, Z)		
CALL CBessel_func (n+1, CONJG(nu(1))*kapa*r, V)	

!-----------------------------------------------------
!dlya otricatelnih indexov cilindricheskih funkcii

IF ( n-1 < 0 ) THEN 

	Z1 = (-1)**(n-1)*Z(-(n-1))
	V1 = (-1)**(n-1)*V(-(n-1))

	ELSE

		Z1 = Z(n-1)
		V1 = V(n-1)

END IF

IF ( n-2 < 0 ) THEN 

	Z2 = (-1)**(n-2)*Z(-(n-2))
	V2 = (-1)**(n-2)*V(-(n-2))

	ELSE

		Z2 = Z(n-2)
		V2 = V(n-2)

END IF

!----------------------------------------------------------
			
u = r/kapa*(CONJG(nu(1))*Z1*V2 - nu(1)*Z2*V1)/(nu(1)**2 - (CONJG(nu(1)))**2) + r/kapa*(CONJG(nu(1))*Z(n+1)*V(n) - nu(1)*Z(n)*V(n+1))/(nu(1)**2 - (CONJG(nu(1)))**2)

END FUNCTION Integral_roZpZp

!====================================================================


FUNCTION Integral_roZpVp_1( kapa, gama, s, r) RESULT (u)


USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama, r

COMPLEX(8)								:: Z(0:n+2)
COMPLEX(8)								:: V(0:n+2)

COMPLEX(8)								:: u, Z1, Z2, V1, V2

COMPLEX(8)								:: nu(1:M+1)

INTEGER(4)								:: i

INTEGER(4), INTENT(IN)					:: s

!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


CALL CBessel_func (n+2, nu(s)*kapa*r, Z)
CALL CHankel_func (n+2, nu(s)*kapa*r, V)

!--------------------------------------------------
!dlya otricatelnih indexov cilindricheskih funkcii

IF ( n-1 < 0 ) THEN 

	Z1 = (-1)**(n-1)*Z(-(n-1))
	V1 = (-1)**(n-1)*V(-(n-1))

	ELSE

		Z1 = Z(n-1)
		V1 = V(n-1)

END IF

IF ( n-2 < 0 ) THEN 

	Z2 = (-1)**(n-2)*Z(-(n-2))
	V2 = (-1)**(n-2)*V(-(n-2))

	ELSE

		Z2 = Z(n-2)
		V2 = V(n-2)

END IF

!--------------------------------------------

u = r**2/4.0*(2.0*Z(n+1)*CONJG(V(n+1)) - Z(n)*CONJG(V(n+2)) - CONJG(V(n))*Z(n+2) ) + r**2/4.0*(2.0*Z1*CONJG(V1) - Z2*CONJG(V(n))-Z(n)*CONJG(V2))


END FUNCTION Integral_roZpVp_1

!===================================================================


FUNCTION Integral_roZpVp_2( kapa, gama, s, r) RESULT (u)


USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama, r

COMPLEX(8)								:: Z(0:n+2)
COMPLEX(8)								:: V(0:n+2)

COMPLEX(8)								:: u, V1, V2, Z1, Z2

COMPLEX(8)								:: nu(1:M+1)

INTEGER(4)								:: i


INTEGER(4), INTENT(IN)					:: s

!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO

CALL CBessel_func (n+2, nu(s)*kapa*r, Z)
CALL CHankel_func (n+2, nu(s)*kapa*r, V)

IF ( n-1 < 0 ) THEN 

	Z1 = (-1)**(n-1)*Z(-(n-1))
	V1 = (-1)**(n-1)*V(-(n-1))

	ELSE

		Z1 = Z(n-1)
		V1 = V(n-1)

END IF

IF ( n-2 < 0 ) THEN 

	Z2 = (-1)**(n-2)*Z(-(n-2))
	V2 = (-1)**(n-2)*V(-(n-2))

	ELSE

		Z2 = Z(n-2)
		V2 = V(n-2)

END IF
			

u = r**2/4.0*(2.0*Z(n+1)*V(n+1) - Z(n)*V(n+2) - V(n)*Z(n+2) ) + r**2/4.0*(2.0*Z1*V1 - Z2*V(n)-V2*Z(n))


END FUNCTION Integral_roZpVp_2



!====================================================================

FUNCTION Integral_roZpZp_1(kapa, gama, s, r) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama, r

COMPLEX(8)								:: Z(0:n+2)

COMPLEX(8)								:: u, Z1, Z2

COMPLEX(8)								:: nu(1:M+1)

INTEGER(4)								:: i


INTEGER(4), INTENT(IN)					:: s

!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


CALL CBessel_func (n+2, nu(s)*kapa*r, Z)	

IF ( n-1 < 0 ) THEN 

	Z1 = (-1)**(n-1)*Z(-(n-1))

	ELSE

		Z1 = Z(n-1)

END IF

IF ( n-2 < 0 ) THEN 

	Z2 = (-1)**(n-2)*Z(-(n-2))

	ELSE

		Z2 = Z(n-2)

END IF


u = r**2/2.0*(Z(n+1)*Z(n+1) - Z(n)*Z(n+2)) + r**2/2.0*(Z1*Z1 - Z2*Z(n))

END FUNCTION Integral_roZpZp_1

!====================================================================



FUNCTION Integral_roZpZp_2(kapa, gama, s, r) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama, r

COMPLEX(8)								:: Z(0:n+2)

COMPLEX(8)								:: V(0:n+2)


COMPLEX(8)								:: u, Z1, Z2, V1, V2

COMPLEX(8)								:: nu(1:M+1)

INTEGER(4)								:: i


INTEGER(4), INTENT(IN)					:: s

!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


CALL CHankel_func (n+2, nu(s)*kapa*r, Z)
CALL CHankel_func (n+2, nu(s)*kapa*r, V)


IF ( n-1 < 0 ) THEN 

	Z1 = (-1)**(n-1)*Z(-(n-1))
	V1 = (-1)**(n-1)*V(-(n-1))

	ELSE

		Z1 = Z(n-1)
		V1 = V(n-1)

END IF

IF ( n-2 < 0 ) THEN 

	Z2 = (-1)**(n-2)*Z(-(n-2))
	V2 = (-1)**(n-2)*V(-(n-2))

	ELSE

		Z2 = Z(n-2)
		V2 = V(n-2)

END IF


u = r**2/4.0*(2.0*Z(n+1)*CONJG(V(n+1)) - Z(n)*CONJG(V(n+2)) - Z(n+2)*CONJG(V(n)) ) + r**2/4.0*(2.0*Z1*CONJG(V1) - Z2*CONJG(V(n))-Z(n)*CONJG(V2))


END FUNCTION Integral_roZpZp_2


!######################################################################################

FUNCTION W_a (kapa, gama) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

REAL(8)										:: r1, r2

COMPLEX(8)									:: u, nu

COMPLEX(8), DIMENSION (2*M, 2*M)			:: A

COMPLEX(8), DIMENSION (1:2*M)				:: x

COMPLEX(8), DIMENSION (1:M+1)				:: A_n, B_n

COMPLEX(8), EXTERNAL						:: Integral_roZpZp

!----------------------------------------------------

CALL Initialization (alfa, ro, M)

CALL get_matrix (kapa, gama, A, M) 

CALL eigenvector (kapa, gama, 2*M, A, x)	 ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

CALL coefficients (kapa, gama, M, x, A_n, B_n)

nu = DCMPLX(alfa(1), -gama)

r1 = 0.0D0
r2 = ro(1)
			
u = (alfa(1))**2/nu/CONJG(nu)*CDABS(A_n(1))**2*( Integral_roZpZp(kapa, gama, r2)-Integral_roZpZp(kapa, gama, r1) )

!----------------------------------

!For M=1, OT, Single cavity
!u = alfa(1)/nu/CONJG(nu)*( Integral_roZpZp(kapa, gama, r2)-Integral_roZpZp(kapa, gama, r1) )


END FUNCTION W_a

!====================================================================

FUNCTION W_ext (kapa, gama) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

REAL(8)										:: r1, r2

COMPLEX(8)									:: u, u1, u2, u3, u4

INTEGER(4)									:: s

COMPLEX(8), DIMENSION (2*M, 2*M)			:: A

COMPLEX(8), DIMENSION (1:2*M)				:: x

COMPLEX(8), DIMENSION (1:M+1)				:: A_n, B_n

COMPLEX(8), EXTERNAL						:: Integral_roZpVp_1, Integral_roZpZp_2, Integral_roZpVp_2

COMPLEX(8), EXTERNAL						:: Integral_roZpZp_1


COMPLEX(8)								:: nu(1:M+1)

INTEGER(4)								:: i


!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


CALL get_matrix (kapa, gama, A, M) 

CALL eigenvector (kapa, gama, 2*M, A, x)	 ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

CALL coefficients (kapa, gama, M, x, A_n, B_n)


u = (0.0D0, 0.0D0)

DO s = 2, M, 2

r1 = ro(s-1)
r2 = ro(s)

			
u1 = Integral_roZpZp_1( kapa, gama, s, r2) - Integral_roZpZp_1( kapa, gama, s, r1)

u2 = Integral_roZpVp_1( kapa, gama, s, r2) - Integral_roZpVp_1( kapa, gama, s, r1)

u3 = Integral_roZpVp_2( kapa, gama, s, r2) - Integral_roZpVp_2( kapa, gama, s, r1)

u4 = Integral_roZpZp_2( kapa, gama, s, r2) - Integral_roZpZp_2( kapa, gama, s, r1)

u = u + (alfa(s))**2/nu(s)/conjg(nu(s))*(CDABS(A_n(s))**2*u1 + A_n(s)*CONJG(B_n(s))*u2 + B_n(s)*CONJG(A_n(s))*u3 + CDABS(B_n(s))**2*u4 )

END DO


END FUNCTION W_ext

!=============================================================================

FUNCTION W_p (kapa, gama) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

REAL(8)										:: r1, r2

COMPLEX(8)									:: u, u1, u2, u3, u4

INTEGER(4)									:: s

COMPLEX(8), DIMENSION (2*M, 2*M)			:: A

COMPLEX(8), DIMENSION (1:2*M)				:: x

COMPLEX(8), DIMENSION (1:M+1)				:: A_n, B_n

COMPLEX(8), EXTERNAL						:: Integral_roZpVp_1, Integral_roZpVp_2, Integral_roZpZp_2

COMPLEX(8), EXTERNAL						:: Integral_roZpZp_1

COMPLEX(8)									:: nu(1:M+1)

INTEGER(4)									:: i


!----------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO

CALL get_matrix (kapa, gama, A, M) 

CALL eigenvector (kapa, gama, 2*M, A, x)	 ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

CALL coefficients (kapa, gama, M, x, A_n, B_n)


u = (0.D0, 0.D0)

DO s = 3, M, 2

r1 = ro(s-1)
r2 = ro(s)

			
u1 = Integral_roZpZp_1( kapa, gama, s, r2) - Integral_roZpZp_1( kapa, gama, s, r1)

u2 = Integral_roZpVp_1( kapa, gama, s, r2) - Integral_roZpVp_1( kapa, gama, s, r1)

u3 = Integral_roZpVp_2( kapa, gama, s, r2) - Integral_roZpVp_2( kapa, gama, s, r1)

u4 = Integral_roZpZp_2( kapa, gama, s, r2) - Integral_roZpZp_2( kapa, gama, s, r1)

u = u + (alfa(s))**2/nu(s)/conjg(nu(s))*(CDABS(A_n(s))**2*u1 + A_n(s)*CONJG(B_n(s))*u2 + B_n(s)*CONJG(A_n(s))*u3 + CDABS(B_n(s))**2*u4 )

END DO

END FUNCTION W_p



FUNCTION W_sum (kapa, gama) RESULT (u)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama

COMPLEX(8)								:: u

COMPLEX(8), EXTERNAL					:: W_p, W_a, W_ext

!----------------------------------------------------

			
u = W_a(kapa, gama) + W_ext (kapa, gama) + W_p (kapa, gama)


!------------------------------------------

!Circular gain area M=2
!u = W_a(kapa, gama) + W_ext (kapa, gama)


END FUNCTION W_sum

!====================================================================

!====================================================================
!====================================================================
!====================================================================

!! OPTICAL THEOREM FOR SINGLE ACTIVE MICROCAVITY

!!		M = 1

!====================================================================
!====================================================================
!====================================================================


SUBROUTINE OT (kapa, gama)

USE DFIMSL 
USE init_data_class


IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama

COMPLEX(8)								:: Z(0:n)
COMPLEX(8)								:: V(0:n)
COMPLEX(8)								:: nu, w, u

REAL(8)									:: dif	

COMPLEX(8), EXTERNAL					:: W_a

!----------------------------------------------------
CALL Initialization (alfa, ro, M)

nu = DCMPLX(alfa(1), -gama)


print*, 'kapa=', kapa
print*, 'gama=', gama

CALL CBessel_func (n, kapa*nu, Z)
CALL Hankel_func (n, kapa, V)

u = Z(n)/V(n)

w = W_a (kapa, gama)*kapa*kapa*PI


dif = gama - DREAL(2*u*CONJG(u)/w)


print*, "gama - u/w = ", DABS(dif)			

print*, "u/w = ", 2*u*conjg(u)/w			


END SUBROUTINE OT