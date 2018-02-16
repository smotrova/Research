!===============================================================================================

SUBROUTINE get_matrix (kapa, gama, A, M) 

USE init_data_class, ONLY : n, alfa, ro, initialization

IMPLICIT NONE

REAL(8), INTENT(IN)										:: kapa, gama

INTEGER(4), INTENT(IN)									:: M

COMPLEX(8), DIMENSION(1:2*M, 1:2*M), INTENT(OUT)		:: A 

INTEGER(4)												:: i, s,j

COMPLEX(8)												:: nu(1:M+1)

COMPLEX(8)												:: CJ0(0:n+1), Hnk1(0:n+1), Hnk2(0:n+1)

COMPLEX(8)												:: CJ1(0:n+1), CJ2(0:n+1)


!---------------------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


!----------------------------------------------------------------------------------

A = (0.0D0, 0.0D0)


DO s = 1, M

	IF (s .ne. 1 .and. s.ne. M) THEN

		CALL CBessel_func (n+1, nu(s)*kapa*ro(s), CJ1)		

		CALL CHankel_func (n+1, nu(s)*kapa*ro(s), Hnk1)

		CALL CBessel_func (n+1, nu(s+1)*kapa*ro(s), CJ2)		

		CALL CHankel_func (n+1, nu(s+1)*kapa*ro(s), Hnk2)


		a(2*s-1, 2*s-2)   = CJ1(n)
		a(2*s-1, 2*s-2+1) = Hnk1(n)
		a(2*s-1, 2*s-2+2) = -CJ2(n)
		a(2*s-1, 2*s-2+3) = -Hnk2(n)

		a(2*s, 2*s-2) = 1/nu(s)*(-CJ1(n+1) + n/nu(s)/kapa/ro(s)*CJ1(n))
		a(2*s, 2*s-2+1) = 1/nu(s)*(-Hnk1(n+1) + n/nu(s)/kapa/ro(s)*Hnk1(n))
		a(2*s, 2*s-2+2) = -1/nu(s+1)*(-CJ2(n+1) + n/nu(s+1)/kapa/ro(s)*CJ2(n))
		a(2*s, 2*s-2+3) = -1/nu(s+1)*(-Hnk2(n+1) + n/nu(s+1)/kapa/ro(s)*Hnk2(n))



	END IF

	IF (s == 1 .and. M .ne. 1) THEN

		CALL CBessel_func (n+1, nu(s)*kapa*ro(s), CJ0)		

		CALL CBessel_func (n+1, nu(s+1)*kapa*ro(s), CJ2)		

		CALL CHankel_func (n+1, nu(s+1)*kapa*ro(s), Hnk2)


		a(2*s-1, s) = CJ0(n)
		a(2*s-1, s+1) =-CJ2(n)
		a(2*s-1, s+2) = -Hnk2(n)

		a(2*s, s) = 1/nu(s)*(-CJ0(n+1) + n/nu(s)/kapa/ro(s)*CJ0(n))
		a(2*s, s+1) = -1/nu(s+1)*(-CJ2(n+1) + n/nu(s+1)/kapa/ro(s)*CJ2(n))
		a(2*s, s+2) = -1/nu(s+1)*(-Hnk2(n+1) + n/nu(s+1)/kapa/ro(s)*Hnk2(n))


	END IF

	IF (s == M .and. s.ne.1) THEN

		CALL CBessel_func (n+1, nu(s)*kapa*ro(s), CJ1)		

		CALL CHankel_func (n+1, nu(s)*kapa*ro(s), Hnk1)

		CALL CHankel_func (n+1, nu(s+1)*kapa*ro(s), Hnk2)


		a(2*s-1, 2*s-2) = CJ1(n)
		a(2*s-1, 2*s-2+1) = Hnk1(n)
		a(2*s-1, 2*s-2+2) = -Hnk2(n)

		a(2*s, 2*s-2) = 1/nu(s)*(-CJ1(n+1) + n/nu(s)/kapa/ro(s)*CJ1(n))
		a(2*s, 2*s-2+1) = 1/nu(s)*(-Hnk1(n+1) + n/nu(s)/kapa/ro(s)*Hnk1(n))
		a(2*s, 2*s-2+2) = -1/nu(s+1)*(-Hnk2(n+1) + n/nu(s+1)/kapa/ro(s)*Hnk2(n))


	END IF


	IF (s == 1 .and. M==1) THEN

		CALL CBessel_func (n+1, nu(s)*kapa*ro(s), CJ0)		
		CALL CHankel_func (n+1, nu(s+1)*kapa*ro(s), Hnk2)

		a(2*s-1, s) = CJ0(n)
		a(2*s-1, s+1) = -Hnk2(n)
		a(2*s, s) = 1/nu(s)*(-CJ0(n+1) + n/nu(s)/kapa/ro(s)*CJ0(n))
		a(2*s, s+1) = -1/nu(s+1)*(-Hnk2(n+1) + n/nu(s+1)/kapa/ro(s)*Hnk2(n))

	END IF

END DO

!-------------- VIVOD MATRICI NA PECHAT ------------------------------------------------------

!				DO i=1,2*M
!					DO j=1, 2*M
!						 WRITE(*,'(1X,i4, 1X,i4, 3X,E15.8, 3X,E15.8)'), i,j, a(i,j)
!					END DO
!				END DO


!------------------------------------------------------------------------------------------------

END SUBROUTINE get_matrix		

!===========================================================================================================================================
!===========================================================================================================================================

