!================================================================================================

SUBROUTINE exclusion (kapa, gama, A, N, A_new, b, i_max, j_max)

IMPLICIT NONE

REAL(8), INTENT(IN)								:: kapa, gama

INTEGER(4), INTENT(IN)							:: N

COMPLEX(8), INTENT(INOUT), DIMENSION (0:N, 0:N)		:: A

COMPLEX(8), INTENT(OUT), DIMENSION (0:N-1, 0:N-1)	:: A_new

INTEGER(4)										:: i, j, ii, jj

INTEGER(4), INTENT(OUT)							:: i_max, j_max

COMPLEX(8), INTENT(OUT), DIMENSION (0:N-1)		:: b

COMPLEX(8)										:: det_1, det_max

COMPLEX(8), EXTERNAL							:: det

!-------------------------------------------------------------------

i_max=1
j_max=1

DO i=0, N
	DO j=0, N

		IF (i<i_max .and. j<j_max) a_new(i, j) = a(i, j) 
			
		IF (i<i_max .and. j>j_max) a_new(i, j-1) = a(i, j)
			
		IF (i>i_max .and. j<j_max) a_new(i-1, j) = a(i, j)
			
		IF (i>i_max .and. j>j_max) a_new(i-1, j-1) = a(i, j)

	END DO
END DO

det_max=det(A_new, N-1)

!--------------------------------------------------------------------

! Nah. nomerov stroki i stolbca takih, chto det mat bez etih stroki/stolbca prinimaet maximalnoe znachenie

DO ii=0, N

DO jj=0, N 


	DO i=0, N
		DO j=0, N

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

	DO i=0, N				!Isklyuchenie 'i_max' stroki i 'j_max' stolbca iz matrici A 

		DO j=0, N

			IF (i<i_max .and. j<j_max) a_new(i, j) = a(i, j) 
			
			IF (i<i_max .and. j>j_max) a_new(i, j-1) = a(i, j)
			
			IF (i>i_max .and. j<j_max) a_new(i-1, j) = a(i, j)
			
			IF (i>i_max .and. j>j_max) a_new(i-1, j-1) = a(i, j)

		END DO

	END DO
		
!----------------------------------------------------------------------------------------------

		DO i=0,N				! Sozdanie stolbca iz elementov stolbca matrici j_max

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

COMPLEX(8), INTENT(IN), DIMENSION (0:N, 0:N)	:: A

COMPLEX(8), DIMENSION (0:N-1, 0:N-1)			:: A_new

COMPLEX(8),  DIMENSION (0:N-1)				:: b, x

COMPLEX(8), INTENT(OUT), DIMENSION (0:N)		:: z

COMPLEX(8), DIMENSION (0:N)					:: res

INTEGER(4)									:: IPATH, LDA

COMPLEX(8)									:: s

!-------------------------------------------------------------------

CALL exclusion (kapa, gama, A, N, A_new, b, i_max, j_max) ! Poluchenie matrici A_new putem isklyucheniya 'i_max' stroki i 'j_max' stolbca

LDA = N
IPATH = 1

CALL DLSACG (N, A_new, LDA, B, IPATH, X)	! Reshenie systmi A_newX = B, gde B vector, sostavlennii iz 'j_max' stolbca matrici A


DO i=0, N			!Nahojdenie sobstv. vectora matrici A: Ay=0

	IF (i < j_max) z(i) = x(i)

	IF (i .eq. j_max) z(i) = 1.0D0

	IF (i > j_max) z(i) = x(i-1)

END DO

!--------Vichislenie znacheniya vectora Ay ----------

!DO i=0, N

!	s = (0.0D0, 0.0D0)
	
!	DO j=0, N

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

SUBROUTINE coefficients (kapa, gama, N, y, A_p, B_p)

USE init_data_class, ONLY : alfa, PI, im

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

INTEGER(4), INTENT(IN)						:: N

INTEGER(4)									:: p

COMPLEX(8), INTENT(IN), DIMENSION (0:N)		:: y

COMPLEX(8), INTENT(OUT), DIMENSION (0:N)		:: A_p, B_p

COMPLEX(8)									:: CJ(0:N+1), Hnk(0:N+1)

REAL(8)										:: J(0:N+1)

COMPLEX(8)									:: V_p, F_p, nu  

!---------------------------------------------------------------------------------------------------

nu = DCMPLX(alfa, -gama)

CALL Hankel_func (N+1, kapa, Hnk)
		
CALL Bessel_func (N+1, kapa, J)		

CALL CBessel_func (N+1, nu*kapa, CJ)	


DO p=0, N

	V_p = (1.0-1.0/nu/nu)*p/kapa*CJ(p)*J(p) + CJ(p+1)*J(p)/nu - CJ(p)*J(p+1)
	F_p = (1.0-1.0/nu/nu)*p/kapa*CJ(p)*Hnk(p) + CJ(p+1)*Hnk(p)/nu - CJ(p)*Hnk(p+1)

	A_p(p) = y(p)/J(p)/F_p
	B_p(p) = -A_p(p)*V_p*Pi*kapa/2/im

END DO

END SUBROUTINE coefficients

!===============================================================================================
!===============================================================================================
! Pole vnutri pervogo resonatora

FUNCTION Ez_in (kapa, gama, N, A_p, r, fi) RESULT (u)

USE init_data_class, ONLY : alfa, im

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

INTEGER(4), INTENT(IN)						:: N

COMPLEX(8), INTENT(IN), DIMENSION (0:N)		:: A_p

INTEGER(4)									:: p

COMPLEX(8)									:: CJ(0:N)

REAL(8), INTENT(IN)							:: r, fi

COMPLEX(8)									:: s, nu, u

!--------------------------------------------------------------


s = (0.0D0, 0.0D0)

nu = DCMPLX(alfa, -gama)

CALL CBessel_func (N, kapa*nu*r, CJ)		

DO p = 0, N

	s = s + A_p(p)*CJ(p)*dcos(p*fi)

END DO

u=s

END FUNCTION Ez_in

!=============================================================================================
!===============================================================================================
!Pole vne resonatorov

FUNCTION Ez_out (kapa, gama, N, B_p, r, fi) RESULT (u)

USE init_data_class, ONLY : im, N_res, PI, zeta

IMPLICIT NONE

REAL(8), INTENT(IN)						:: kapa, gama

INTEGER(4), INTENT(IN)					:: N

INTEGER(4)								:: p, s

COMPLEX(8), INTENT(IN), DIMENSION (0:N)	:: B_p

COMPLEX(8)								:: Hnk(0:N)

REAL(8), INTENT(IN)						:: r, fi

COMPLEX(8)								:: sum_p, sum_s

REAL(8)									:: ro, l, x_s, y_s, r_s, fi_s

COMPLEX(8)								:: u

!------------------------------------------------------------

l = 2.0D0 + zeta

sum_s =	(0.0D0, 0.0D0)

ro = l/2/dsin(PI/N_res)

DO s = 1, N_res

	x_s = r*dcos(fi) + ro*dcos((s-1)*2*PI/N_res)
	y_s = r*dsin(fi) + ro*dsin((s-1)*2*PI/N_res)

	fi_s = datan2(y_s, x_s) - 2*PI*(s-1)/N_res

	r_s  = dsqrt(x_s**2 + y_s**2)

	CALL Hankel_func (N, kapa*r_s, Hnk)		

	sum_p = (0.0D0, 0.0D0)

	DO p = 0, N	
	
		sum_p = sum_p + B_p(p)*Hnk(p)*dcos(p*fi_s)

	END DO

	sum_s = sum_s + sum_p

END DO

u = sum_s

END FUNCTION Ez_out

!===============================================================================================
! Postroenie blijnego pola

SUBROUTINE eigenfields (kapa, gama, N, A)

USE init_data_class, ONLY : PI, zeta, N_res

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa, gama

INTEGER(4), INTENT(IN)						:: N

COMPLEX(8), INTENT(IN), DIMENSION (N, N)	:: A

COMPLEX(8), DIMENSION (0:N)					:: x

COMPLEX(8), DIMENSION (0:N)					:: A_p, B_p

REAL(8)										:: y, z, r, fi, u1, u2, ro

REAL(8)										:: E_0, E_max, l, x_s(1:N_res), y_s(1:N_res), fi_s(1:N_res), r_s(1:N_res)

INTEGER(4)									:: p, s, key

COMPLEX(8)									:: u, v

COMPLEX(8), EXTERNAL						:: Ez_in, Ez_out

!------------------------------------------------------------------------

CALL eigenvector (kapa, gama, N, A, x) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

CALL coefficients (kapa, gama, N, x, A_p, B_p)
!------------------------------------------------------------------------
!nah. max pola

l = 2.0D0 + zeta
ro = l/2/dsin(PI/N_res)

E_0 = 0.0D0
E_max = E_0

	
		DO y = -2 - ro, 2 + ro, 0.05
			DO z = -2 - ro, 2 + ro, 0.05

				fi = datan2(y, z)
				r  = dsqrt(z**2 + y**2)

				DO s = 1, N_res	

					x_s(s) = z + ro*dcos((s-1)*2*PI/N_res)
					y_s(s) = y + ro*dsin((s-1)*2*PI/N_res)

					fi_s(s) = datan2(y_s(s), x_s(s))- 2*PI*(s-1)/N_res
					r_s(s)  = dsqrt(x_s(s)**2 + y_s(s)**2)

				END DO

				
				key = 0

				DO s = 1, N_res

					IF ( r_s(s) <= 1.0D0 ) THEN

						u = Ez_in (kapa, gama, N, A_p, r_s(s), fi_s(s))
						key = 1

					END IF

				END DO

				IF (key .eq. 0) u = Ez_out(kapa, gama, N, B_p, r, fi)  !Pole vne
					
				E_0 = cdabs(u)

				IF (E_0 > E_max) E_max = E_0

			END DO
		END DO

!=============================================================================================
!=============================================================================================
!postroenie normirovannogo blijnego pola

OPEN (unit=7, file='H_FIELD.dat')

		DO y = -2 - ro, 2 + ro, 0.05
			DO z = -2 - ro, 2 + ro, 0.05

				fi = datan2(y, z)
				r  = dsqrt(z**2 + y**2)

				DO s = 1, N_res		! vichislaem koordinati tochki v lokalnih koordinatah

						x_s(s) = z + ro*dcos((s-1)*2*PI/N_res)
						y_s(s) = y + ro*dsin((s-1)*2*PI/N_res)

						fi_s(s) = datan2(y_s(s), x_s(s)) - 2*PI*(s-1)/N_res
						r_s(s)  = dsqrt(x_s(s)**2 + y_s(s)**2)

				END DO

				
				key = 0

				DO s = 1, N_res		! postroenie pola vnutri resonatorov

					IF ( r_s(s) <= 1.0D0 ) THEN

						u = Ez_in (kapa, gama, N, A_p, r_s(s), fi_s(s))
						key = 1

					END IF

				END DO

				IF (key .eq. 0) u = Ez_out(kapa, gama, N, B_p, r, fi)  ! vneshnee pole

				u1=dreal(u)
				u2=aimag(u)

!				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(z), sngl(y), cdabs(u)/E_max
						
				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(z), sngl(y), sign(cdabs(u)/E_max, datan2(u2, u1))

			END DO
		END DO

close (7)

!-----------------------------------------------------------------

!Postroenie konturov resonatorov

	DO s = 1, N_res

		DO fi = 0, 360, 0.1
			
			x_s(s) = dcosd(fi) - ro*dcos((s-1)*2*PI/N_res)
			y_s(s) = dsind(fi) + ro*dsin((s-1)*2*PI/N_res)

			WRITE(10+s,'(2X,E13.6, 3X,E13.6 )'), sngl(x_s(s)), sngl(y_s(s))

		END DO

	END DO

!-------------------------------------------------------------------

END SUBROUTINE eigenfields

!===============================================================================================
!===============================================================================================