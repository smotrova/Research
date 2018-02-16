!===============================================================================================

SUBROUTINE get_matrix (kapa, gama, A, N) 

USE init_data_class, ONLY : alfa, PI, im, zeta, N_res

IMPLICIT NONE

REAL(8), INTENT(IN)								:: kapa, gama

INTEGER(4), INTENT(IN)							:: N

COMPLEX(8), DIMENSION(0:N, 0:N), INTENT(OUT)	:: A 

INTEGER(4)										:: m, p, s

COMPLEX(8)										::  VJF, nu, V_p, F_p

COMPLEX(8)										:: Hnk(0:N+1), CJ(0:N+1)

REAL(8)											:: J(0:N+1)

COMPLEX(8)										:: Hnk1(0:2*N)

COMPLEX(8)										:: h_plus, h_minus, sum

REAL(8)											:: l, teta_plus, teta_minus, ro_1s, mu

!---------------------------------------------------------------

l = 2.0D0 + zeta

nu = DCMPLX(alfa, -gama)

CALL Bessel_func (N+1, kapa, J)	

CALL CBessel_func (N+1, nu*kapa, CJ)		

CALL Hankel_func (N+1, kapa, Hnk)

!----------------------------------------------------------------------------------

	DO m=0,N

		DO p=0,N

			V_p = (1.0-1.0/nu/nu)*p/kapa*CJ(p)*J(p) + CJ(p+1)*J(p)/nu - CJ(p)*J(p+1)
			F_p = (1.0-1.0/nu/nu)*p/kapa*CJ(p)*Hnk(p) + CJ(p+1)*Hnk(p)/nu - CJ(p)*Hnk(p+1)

			IF (p .eq. 0) mu = 0.5D0
			IF (p .ne. 0) mu = 1.0D0

			VJF = V_p/J(p)/F_p

			sum = (0.0D0, 0.0D0)

			DO s = 2, N_res

				ro_1s = l/dsin(PI/N_res)*dsin((s-1)*PI/N_res) 

				CALL Hankel_func (2*N, kapa*ro_1s, Hnk1)
	
				IF (m-p >= 0)  THEN
				
					h_minus = Hnk1(m-p)  
								
				ELSE
				
					h_minus = (-1)**(p-m)*Hnk1(p-m)  
             
				ENDIF

				h_plus = Hnk1(m+p)  
				
				teta_plus = PI/2 + PI*(s-1)/N_res

				teta_minus = PI/2 - PI*(s-1)/N_res

				sum = (-1)**(s-1)*((-1)**p*h_minus*dcos(p*teta_plus - m*teta_minus) + h_plus*dcos(p*teta_plus + m*teta_minus)) + sum

			END DO


			IF (abs(m-p) .ne. 0) A(m,p) = (-1)**m*J(m)*VJF*sum*mu
			IF (abs(m-p) .eq. 0) A(m,p) = 1.0D0 + (-1)**m*J(m)*VJF*sum*mu

		END DO

	END DO

!------------------------------------------------------------------------------------------------

END SUBROUTINE get_matrix		

!===========================================================================================================================================
!===========================================================================================================================================

