!===============================================================================================

SUBROUTINE get_matrix (kapa, gama, U, N) 

USE init_data_class, ONLY : PI, alfa_i, alfa_e, im


IMPLICIT NONE

REAL(8), INTENT(IN)										:: kapa, gama

INTEGER(4), INTENT(IN)									:: N				!Matrix size


COMPLEX(8), DIMENSION(0:4*N-1, 0:4*N-1), INTENT(OUT)	:: U 


COMPLEX(8), DIMENSION(0:2*N-1, 0:2*N-1)					:: A, B, C, D 


INTEGER(4)												:: m, p, s


COMPLEX(8)												:: nu, eta


COMPLEX(8)												:: A1, A2, B1, B2, C1, C2, D1, D2


REAL(8)													:: f, t_p, t_s, P_N, w_p

REAL(8), EXTERNAL										:: L

!---------------------------------------------------------------


DO s = 0, 2*N-1
	
	t_s = s*PI/N


	DO p = 0, 2*N-1

		t_p = p*PI/N

		CALL KERNELS (kapa, gama, t_s, t_p, A1, A2, B1, B2, C1, C2, D1, D2)

		f = L(t_p)


			P_n = 0.0D0

			DO m = 1, N-1

				P_n = P_n + dcos(m*(t_s - t_p))/m

			END DO

			P_n = -2.0D0*PI*P_n/N - PI*dcos(N*(t_s - t_p))/(N*N)


		A(s, p) = ( A1*P_n + PI/N*A2)*f

		B(s, p) = ( B1*P_n + PI/N*B2)*f

		C(s, p) = ( C1*P_n + PI/N*C2)*f

		D(s, p) = ( D1*P_n + PI/N*D2)*f

	END DO

END DO

!********************************************************

nu = DCMPLX(alfa_i, -gama)

!-------------------------------------------------------

! H pol

eta = ( nu*nu + alfa_e*alfa_e)/(2.0D0*nu*nu)

!eta = (1.0D0, 0.0D0)

!-------------------------------------------------------

DO s = 0, 4*N-1
	DO p = 0, 4*N-1

		IF (p <= 2*N-1 .and. s <= 2*N-1) THEN

				IF (s == p) U(s, p) = 1.0D0 + A(s,p)
				IF (s .ne. p) U(s, p) = A(s,p)

		END IF


		IF (p > 2*N-1 .and. s <= 2*N-1) U(s, p) = -B(s, p - 2*N)/eta
	

		IF (p <= 2*N-1 .and. s > 2*N-1) U(s, p) =  C(s - 2*N, p)
				

		IF (p > 2*N-1 .and. s > 2*N-1) THEN

				IF (s == p) U(s, p) = 1.0D0 - D(s - 2*N,p - 2*N)/eta
				IF (s .ne. p) U(s, p) = -D(s - 2*N,p - 2*N)/eta

		END IF

	END DO
END DO

!----------------------------------------------------------------------

END SUBROUTINE get_matrix		

!===========================================================================================================================================
