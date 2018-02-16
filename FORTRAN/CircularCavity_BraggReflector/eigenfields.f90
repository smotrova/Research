!================================================================================================

SUBROUTINE exclusion (kapa, gama, A, N, A_new, b, i_max, j_max)

IMPLICIT NONE

REAL(8), INTENT(IN)								:: kapa, gama

INTEGER(4), INTENT(IN)							:: N

COMPLEX(8), INTENT(INOUT), DIMENSION (N, N)		:: A

COMPLEX(8), INTENT(OUT), DIMENSION (N-1, N-1)	:: A_new

INTEGER(4)										:: i, j, ii, jj

INTEGER(4), INTENT(OUT)							:: i_max, j_max

COMPLEX(8), INTENT(OUT), DIMENSION (N-1)		:: b

COMPLEX(8)										:: det_1, det_max

COMPLEX(8), EXTERNAL							:: det

!-------------------------------------------------------------------

i_max=1
j_max=1

DO i=1, N
	DO j=1, N

		IF (i<i_max .and. j<j_max) a_new(i, j) = a(i, j) 
			
		IF (i<i_max .and. j>j_max) a_new(i, j-1) = a(i, j)
			
		IF (i>i_max .and. j<j_max) a_new(i-1, j) = a(i, j)
			
		IF (i>i_max .and. j>j_max) a_new(i-1, j-1) = a(i, j)

	END DO
END DO

det_max=det(A_new, N-1)

!--------------------------------------------------------------------

! Nah. nomerov stroki i stolbca takih, chto det mat bez etih stroki/stolbca prinimaet maximalnoe znachenie

DO ii=1, N

DO jj=1, N 


	DO i=1, N
		DO j=1, N

			IF (i<ii .and. j<jj) a_new(i, j) = a(i, j) 
			
			IF (i<ii .and. j>jj) a_new(i, j-1) = a(i, j)
			
			IF (i>ii .and. j<jj) a_new(i-1, j) = a(i, j)
			
			IF (i>ii .and. j>jj) a_new(i-1, j-1) = a(i, j)

		END DO
	END DO

	det_1=det(A_new, N-1)

	IF (cdabs(det_1) > cdabs(det_max))  THEN
		
		det_max = det_1
		j_max=jj
		i_max=ii

	END IF

END DO
END DO

!---------------------------------------------------------------------------------------------

	DO i=1, N				!Isklyuchenie 'i_max' stroki i 'j_max' stolbca iz matrici A 

		DO j=1, N

			IF (i<i_max .and. j<j_max) a_new(i, j) = a(i, j) 
			
			IF (i<i_max .and. j>j_max) a_new(i, j-1) = a(i, j)
			
			IF (i>i_max .and. j<j_max) a_new(i-1, j) = a(i, j)
			
			IF (i>i_max .and. j>j_max) a_new(i-1, j-1) = a(i, j)

		END DO

	END DO
		
!----------------------------------------------------------------------------------------------

		DO i=1,N				! Sozdanie stolbca iz elementov stolbca matrici j_max

			IF (i < i_max)	b(i) = -a(i, j_max)
			IF (i > i_max)	b(i-1) = -a(i, j_max)
							
		END DO

END SUBROUTINE exclusion

!===============================================================================================
!===============================================================================================

SUBROUTINE eigenvector (kapa, gama, N, A, z)

USE DFIMSLMD

IMPLICIT NONE

INTEGER(4)									:: i, j, i_max, j_max

REAL(8), INTENT(IN)							:: kapa, gama

INTEGER(4), INTENT(IN)						:: N

COMPLEX(8), INTENT(IN), DIMENSION (N, N)	:: A

COMPLEX(8), DIMENSION (N-1, N-1)			:: A_new

COMPLEX(8),  DIMENSION (N-1)				:: b, x

COMPLEX(8), INTENT(OUT), DIMENSION (N)		:: z

COMPLEX(8), DIMENSION (N)					:: res

INTEGER(4)									:: IPATH, LDA

COMPLEX(8)									:: s

!-------------------------------------------------------------------

CALL exclusion (kapa, gama, A, N, A_new, b, i_max, j_max) ! Poluchenie matrici A_new putem isklyucheniya 'i_max' stroki i 'j_max' stolbca

LDA = N-1
IPATH = 1

CALL DLSACG (N-1, A_new, LDA, B, IPATH, X)	! Reshenie systmi A_newX = B, gde B vector, sostavlennii iz 'j_max' stolbca matrici A


DO i=1, N			!Nahojdenie sobstv. vectora matrici A: Ay=0

	IF (i < j_max) z(i) = x(i)

	IF (i .eq. j_max) z(i) = 1.0D0

	IF (i > j_max) z(i) = x(i-1)

END DO

!--------Vichislenie znacheniya vectora Ay ----------

!DO i=1, N

!	s = (0.0D0, 0.0D0)
	
!	DO j=1, N

!		s= s + a(i,j)*y(j)

!	END DO

!	res(i) = s

!	print*,'', i, cdabs(res(i))

!END DO
!-----------------------------------------------------

END SUBROUTINE eigenvector

!=============================================================================================
!=============================================================================================

!Vichislenie coeficientov razlojeniya pola v radi
! y: Ay=0

SUBROUTINE coefficients (kapa, gama, N, y, A_n, B_n)

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

INTEGER(4), INTENT(IN)						:: N

INTEGER(4)									:: i

COMPLEX(8), INTENT(IN), DIMENSION (2*N)		:: y

COMPLEX(8), INTENT(OUT), DIMENSION (N+1)		:: A_n, B_n


!---------------------------------------------------------------------------------------------------

A_n(1) = y(1)
B_n(1) = (0.0D0, 0.0D0)

DO i = 1, N-1

	A_n(i+1) = y(2*i) 
	B_n(i+1) = y(2*i+1)

END DO

A_n(N+1) = (0.0D0, 0.0D0)
B_n(N+1) = y(2*N)


END SUBROUTINE coefficients

!===============================================================================================
!===============================================================================================
! Pole vnutri resonatora

FUNCTION Ez (kapa, gama, A_ns, B_ns, nu_s, r, fi) RESULT (u)

USE init_data_class , ONLY : n

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

COMPLEX(8), INTENT(IN)						:: A_ns, B_ns, nu_s

COMPLEX(8)									:: CJ(0:n)

REAL(8), INTENT(IN)							:: r, fi

COMPLEX(8)									:: u

!--------------------------------------------------------------

	CALL CBessel_func (n, nu_s*kapa*r, CJ)		

	u = A_ns*CJ(n)*dcos(n*fi)

END FUNCTION Ez



!=============================================================================================
!=============================================================================================
!=============================================================================================


FUNCTION Ez_1 (kapa, gama, A_ns, B_ns, nu_s, r, fi) RESULT (u)


USE init_data_class , ONLY : n

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

COMPLEX(8), INTENT(IN)						:: A_ns, B_ns

COMPLEX(8), INTENT(IN)						:: nu_s

COMPLEX(8)									:: Hnk(0:n)

COMPLEX(8)									:: CJ(0:n)

REAL(8), INTENT(IN)							:: r, fi

COMPLEX(8)									:: u

!--------------------------------------------------------------

	CALL CHankel_func (n, nu_s*kapa*r, HNK)

	CALL CBessel_func (n, nu_s*kapa*r, CJ)		


	u = A_ns*CJ(n)*dcos(n*fi) + B_ns*Hnk(n)*dcos(n*fi)


END FUNCTION Ez_1

!==================================================================================
!==================================================================================

FUNCTION Ez_2 (kapa, gama, A_ns, B_ns, nu_s, r, fi) RESULT (u)

USE init_data_class , ONLY : n

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

COMPLEX(8), INTENT(IN)						:: A_ns, B_ns

COMPLEX(8), INTENT(IN)						:: nu_s

COMPLEX(8)									:: Hnk(0:n)

REAL(8), INTENT(IN)							:: r, fi

COMPLEX(8)									:: u

!--------------------------------------------------------------

	CALL CHankel_func (n, nu_s*kapa*r, HNK)

	u = B_ns*Hnk(n)*dcos(n*fi)


END FUNCTION Ez_2


!=============================================================================================
!=============================================================================================
!=============================================================================================
! Postroenie blijnego pola

SUBROUTINE eigenfields (kapa, gama, M, A)

USE init_data_class, ONLY : alfa, PI, im, ro, initialization


IMPLICIT NONE

REAL(8), INTENT(IN)									:: kapa, gama

INTEGER(4), INTENT(IN)								:: M

COMPLEX(8), INTENT(IN), DIMENSION (2*M, 2*M)		:: A

COMPLEX(8), DIMENSION (1:2*M)						:: x

COMPLEX(8), DIMENSION (1:M+1)						:: A_n, B_n

COMPLEX(8)											:: nu(1:M+1)

REAL(8)												:: y, z, r, fi

REAL(8)												:: E_0, E_max, x_s, y_s

INTEGER(4)											:: s, i
	
COMPLEX(8)											:: u, v

COMPLEX(8), EXTERNAL								:: Ez, Ez_1, Ez_2

!------------------------------------------------------------------------

CALL eigenvector (kapa, gama, 2*M, A, x)	 ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

CALL coefficients (kapa, gama, M, x, A_n, B_n)

!------------------------------------------------------------------------

CALL Initialization (alfa, ro, M)

nu(1) = DCMPLX(alfa(1), -gama)

DO i= 2, M+1

	nu(i) = DCMPLX( alfa(i), 0.0D0)

END DO


!------------------------------------------------------------------------
!nah. max pola

E_0 = 0.0D0
E_max = E_0

		DO y = -2 - ro(M), 2 + ro(M), 0.1
			DO z = -2 - ro(M), 2 + ro(M), 0.1

				IF (z .ne. 0) fi = datan2(y, z)

				r  = dsqrt(z**2 + y**2)

				IF (r < ro(1)) u = Ez(kapa, gama, A_n(1), B_n(1), nu(1), r, fi)

				DO s = 2, M	! vichislaem koordinati tochki

					IF (r >= ro(s-1) .and. r <= ro(s)) u = Ez_1(kapa, gama, A_n(s), B_n(s), nu(s), r, fi)

				END DO

				IF ( r > ro(M)) u = Ez_2(kapa, gama, A_n(M+1), B_n(M+1), nu(M+1), r, fi)

				E_0 = cdabs(u)

				IF (E_0 > E_max) E_max = E_0	

			END DO
		END DO

!=============================================================================================
!=============================================================================================
!postroenie normirovannogo blijnego pola


OPEN (unit=7, file='H_FIELD.dat')

		DO y = -2 - ro(M), 2 + ro(M), 0.1
			DO z = -2 - ro(M), 2 + ro(M), 0.1

				IF (z .ne. 0) fi = datan2(y, z)

				r  = dsqrt(z**2 + y**2)
	
				IF (r <= ro(1)) u = Ez(kapa, gama, A_n(1), B_n(1), nu(1), r, fi)

				DO s = 2, M	

					IF (r >= ro(s-1) .and. r <= ro(s)) u = Ez_1(kapa, gama, A_n(s), B_n(s), nu(s), r, fi)

				END DO

				IF ( r > ro(M)) u = Ez_2(kapa, gama, A_n(M+1), B_n(M+1), nu(M+1), r, fi)
					
				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(z), sngl(y), cdabs(u)/E_max

			END DO
		END DO

close (7)

!-----------------------------------------------------------------

!Postroenie konturov resonatorov

	DO s = 1, M

		DO fi = 0, 360, 0.1
			
			x_s = ro(s)*dcosd(fi)
			y_s = ro(s)*dsind(fi)

			WRITE(100+s,'(2X,E13.6, 3X,E13.6 )'), sngl(x_s), sngl(y_s)

		END DO

	END DO

!-------------------------------------------------------------------

END SUBROUTINE eigenfields

!===============================================================================================
!===============================================================================================