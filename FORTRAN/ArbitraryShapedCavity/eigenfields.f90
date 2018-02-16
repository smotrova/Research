!================================================================================================

SUBROUTINE exclusion (kapa, gama, A, N, A_new, b, i_max, j_max)

IMPLICIT NONE

REAL(8), INTENT(IN)									:: kapa, gama

INTEGER(4), INTENT(IN)								:: N

COMPLEX(8), INTENT(INOUT), DIMENSION (0:4*N-1, 0:4*N-1)		:: A

COMPLEX(8), INTENT(OUT), DIMENSION (0:4*N-2, 0:4*N-2)	:: A_new

INTEGER(4)											:: i, j, ii, jj

INTEGER(4), INTENT(OUT)								:: i_max, j_max

COMPLEX(8), INTENT(OUT), DIMENSION (0:4*N-2)			:: b

COMPLEX(8)											:: det_1, det_max

COMPLEX(8), EXTERNAL								:: det0

!-------------------------------------------------------------------

i_max = 0
j_max = 0

DO i = 0, 4*N-1
	DO j = 0, 4*N-1

		IF (i<i_max .and. j<j_max) a_new(i, j) = a(i, j) 
			
		IF (i<i_max .and. j>j_max) a_new(i, j-1) = a(i, j)
			
		IF (i>i_max .and. j<j_max) a_new(i-1, j) = a(i, j)
			
		IF (i>i_max .and. j>j_max) a_new(i-1, j-1) = a(i, j)

	END DO
END DO

det_max = det0(A_new, 4*N-2)


!--------------------------------------------------------------------

! Nah. nomerov stroki i stolbca takih, chto det mat bez etih stroki/stolbca prinimaet maximalnoe znachenie

DO ii=0, 4*N-1

DO jj=0, 4*N-1 


	DO i=0, 4*N-1
		DO j=0, 4*N-1

			IF (i<ii .and. j<jj) a_new(i, j) = a(i, j) 
			
			IF (i<ii .and. j>jj) a_new(i, j-1) = a(i, j)
			
			IF (i>ii .and. j<jj) a_new(i-1, j) = a(i, j)
			
			IF (i>ii .and. j>jj) a_new(i-1, j-1) = a(i, j)

		END DO
	END DO

	det_1=det0(A_new, 4*N-2)

	IF (cdabs(det_1) > cdabs(det_max))  THEN
		
		det_max = det_1
		j_max=jj
		i_max=ii

	END IF

END DO
END DO

!---------------------------------------------------------------------------------------------

	DO i=0, 4*N-1				!Isklyuchenie 'i_max' stroki i 'j_max' stolbca iz matrici A 

		DO j=0, 4*N-1

			IF (i<i_max .and. j<j_max) a_new(i, j) = a(i, j) 
			
			IF (i<i_max .and. j>j_max) a_new(i, j-1) = a(i, j)
			
			IF (i>i_max .and. j<j_max) a_new(i-1, j) = a(i, j)
			
			IF (i>i_max .and. j>j_max) a_new(i-1, j-1) = a(i, j)

		END DO

	END DO
		
!----------------------------------------------------------------------------------------------

		DO i=0,4*N-1				! Sozdanie stolbca iz elementov stolbca matrici j_max

			IF (i < i_max)	b(i) = -a(i, j_max)
			IF (i > i_max)	b(i-1) = -a(i, j_max)
							
		END DO

END SUBROUTINE exclusion

!===============================================================================================
!===============================================================================================

SUBROUTINE eigenvector (kapa, gama, N, A, z)

USE DFIMSLMD

IMPLICIT NONE

INTEGER(4)										:: i, j, i_max, j_max

REAL(8), INTENT(IN)								:: kapa, gama

INTEGER(4), INTENT(IN)							:: N

COMPLEX(8), INTENT(IN), DIMENSION (0:4*N-1, 0:4*N-1)	:: A

COMPLEX(8), DIMENSION (0:4*N-2, 0:4*N-2)			:: A_new

COMPLEX(8),  DIMENSION (0:4*N-2)					:: b, x

COMPLEX(8), INTENT(OUT), DIMENSION (0:4*N-1)		:: z

COMPLEX(8), DIMENSION (0:4*N-1)						:: res

INTEGER(4)										:: IPATH, LDA

COMPLEX(8)										:: s

!-------------------------------------------------------------------

CALL exclusion (kapa, gama, A, N, A_new, b, i_max, j_max) ! Poluchenie matrici A_new putem isklyucheniya 'i_max' stroki i 'j_max' stolbca

LDA = 4*N-1
IPATH = 1

CALL DLSACG (4*N-1, A_new, LDA, B, IPATH, X)	! Reshenie systemi A_newX = B, gde B vector, sostavlennii iz 'j_max' stolbca matrici A

DO i = 0, 4*N-1			!Nahojdenie sobstv. vectora matrici A: Ay=0

	IF (i < j_max) z(i) = x(i)

	IF (i .eq. j_max) z(i) = 1.0D0

	IF (i > j_max) z(i) = x(i-1)

END DO


!--------Vichislenie znacheniya vectora Ay ----------

!DO i = 0, 4*N - 1

!	s = (0.0D0, 0.0D0)
	
!	DO j  = 0, 4*N - 1

!		s = s + a(i,j)*z(j)

!	END DO
 
!	res(i) = s

!	print*,'', i, cdabs(res(i))

!END DO
!-----------------------------------------------------

END SUBROUTINE eigenvector


!===============================================================================================
!===============================================================================================

! dlya polya vnutri rasonatora

SUBROUTINE  kernels_field_in (kapa, gama, x1, y1, t_p, W, V)

USE DFIMSL
USE init_data_class

IMPLICIT NONE


REAL(8), INTENT(IN)				:: x1, y1, t_p, kapa, gama

COMPLEX(8),INTENT(OUT)			:: W, V

COMPLEX(8)						:: nu, kapa_i, kapa_iR

REAL(8)							:: x2, y2

REAL(8)							:: R, R_n2

REAL(8), EXTERNAL				:: dx_dt, dy_dt, x, y, L, curv

COMPLEX(8), EXTERNAL			:: Hnk

!-----------------------------------------------------

	nu = DCMPLX(alfa_i, -gama)

	kapa_i = nu*kapa

	x2 = x(t_p)
	y2 = y(t_p)
	
	R = dsqrt( (x1 - x2)**2 + (y1 - y2)**2 )

	R_n2 = ((x1 - x2)*dy_dt(t_p) - (y1 - y2)*dx_dt(t_p) )/L(t_p)

	kapa_iR = kapa_i*R

	IF ( (R .ne. 0.0D0) ) THEN
	
		W = im/4.0D0*kapa_i*HNK(1,kapa_iR)*R_n2/R
		V = im/4.0D0*HNK(0,kapa_iR)

	END IF

END SUBROUTINE  kernels_field_in


!====================================================================================
!====================================================================================

! dlya polya vne rasonatora

SUBROUTINE  kernels_field_out (kapa, gama, x1, y1, t_p, W, V)

USE DFIMSL
USE init_data_class

IMPLICIT NONE

REAL(8), INTENT(IN)				:: x1, y1, t_p, kapa, gama

COMPLEX(8),INTENT(OUT)			:: W, V

REAL(8)							:: x2, y2

REAL(8)							:: kapa_e, kapa_eR

REAL(8)							:: R, R_n2

COMPLEX(8)						:: Hnk0, Hnk1

REAL(8), ALLOCATABLE			:: RJ(:), RY(:)

INTEGER(4)						:: Q, NK, ERR

INTEGER(4), EXTERNAL			:: QMIN

REAL(8), EXTERNAL				:: dx_dt, dy_dt, x, y, L, curv


!-----------------------------------------------------


	x2 = x(t_p)
	y2 = y(t_p)
	
	R = dsqrt( (x1 - x2)**2 + (y1 - y2)**2 )

	R_n2 = ((x1 - x2)*dy_dt(t_p) - (y1 - y2)*dx_dt(t_p) )/L(t_p)

	kapa_e = alfa_e*kapa
	kapa_eR = kapa_e*R


	IF ( (R .ne. 0.0D0) ) THEN


		NK=INT(ceiling(dabs(kapa_eR)) + 15)

		Q = QMIN(kapa_eR,NK+1)

		ALLOCATE(RY(-1:NK+1), RJ(-1:Q+1), STAT=ERR)
		CALL DJBES(Q,kapa_eR,RJ)                  
		CALL DYBES(kapa_eR,RJ,RY,NK,Q)   

			Hnk0 = RJ(0) + im*RY(0)

			Hnk1 = RJ(1) + im*RY(1)
	
			W = im/4.0D0*kapa_e*Hnk1*R_n2/R
			V = im/4.0D0*Hnk0

	
		DEALLOCATE (RJ,RY,stat=ERR)

	END IF


END SUBROUTINE  kernels_field_out


!===============================================================================================
!===============================================================================================

!Pole vnutri resonatora

FUNCTION Ez_in (kapa, gama, N, z, x1, y1) RESULT (u)

USE init_data_class, ONLY : PI, alfa_i, alfa_e

IMPLICIT NONE

REAL(8), INTENT(IN)								:: kapa, gama

INTEGER(4), INTENT(IN)							:: N

COMPLEX(8), INTENT(IN), DIMENSION (0: 4*N-1)		:: z

COMPLEX(8), DIMENSION (0: 2*N-1)				:: fi, psi

INTEGER(4)										:: p

REAL(8), INTENT(IN)								:: x1, y1

REAL(8)											:: t_p, x2, y2, R, R_n2

COMPLEX(8)										:: u, W, V

COMPLEX(8)										:: eta_i, eta, nu
REAL(8)											:: eta_e


REAL(8), EXTERNAL								:: L, x, y, dx_dt, dy_dt


!--------------------------------------------------------------
nu = DCMPLX(alfa_i, -gama)

! H pol

eta_i = 1.0D0/nu/nu
eta_e = 1.0D0/alfa_e/alfa_e

eta = (2.0D0*eta_e)/(eta_i + eta_e)


DO p = 0, 4*N-1

	IF (p <= 2*N-1) fi(p) = z(p)
	IF (p > 2*N-1)	psi (p - 2*N) = z(p)*eta

END DO

u = (0.0D0, 0.0D0)


DO  p = 0, 2*N-1


	t_p = p*PI/N


	x2 = x(t_p)
	y2 = y(t_p)
	
	R = dsqrt( (x1 - x2)**2 + (y1 - y2)**2 )

	IF (R > 1.D-2) THEN

		CALL kernels_field_in (kapa, gama, x1, y1, t_p, W, V)

		u = u + ( W*fi(p) - V*psi(p))*L(t_p) 
	
	END IF

END DO

u = - PI/N*u


END FUNCTION Ez_in



!===============================================================================================
!===============================================================================================
!Pole vne resonatora

FUNCTION Ez_out (kapa, gama, N, z, x1, y1) RESULT (u)

USE init_data_class, ONLY : PI, alfa_i, alfa_e

IMPLICIT NONE

REAL(8), INTENT(IN)								:: kapa, gama

INTEGER(4), INTENT(IN)							:: N

COMPLEX(8), INTENT(IN), DIMENSION (0: 4*N-1)		:: z

COMPLEX(8), DIMENSION (0: 2*N-1)				:: fi, psi

INTEGER(4)										:: p

REAL(8), INTENT(IN)								:: x1, y1

REAL(8)											:: t_p, x2, y2, R, R_n2

COMPLEX(8)										:: u, W, V, nu

COMPLEX(8)										:: eta_i, eta
REAL(8)											:: eta_e

REAL(8), EXTERNAL								:: L, x, y, dx_dt, dy_dt


!*********************************************************

nu = DCMPLX(alfa_i, -gama)

!----------------------------------------------

! H pol

eta_i = 1.0D0/nu/nu
eta_e = 1.0D0/alfa_e/alfa_e

eta = (2.0D0*eta_e)/(eta_i + eta_e)

!-----------------------------------------------

DO p = 0, 4*N-1

	IF (p <= 2*N-1) fi(p) = z(p)
	IF (p > 2*N-1)	psi (p - 2*N) = z(p)*eta

END DO

u = (0.0D0, 0.0D0)

DO  p = 0, 2*N-1


	t_p = PI*p/N

	x2 = x(t_p)
	y2 = y(t_p)
	
	R = dsqrt( (x1 - x2)**2 + (y1 - y2)**2 )

	IF (R > 1.D-2) THEN

    CALL kernels_field_out (kapa, gama, x1, y1, t_p, W, V)

	u = u + ( W*fi(p) - eta_i/eta_e*psi(p)*V)*L(t_p) 

	END IF

END DO

u = PI/N*u

END FUNCTION Ez_out

!===============================================================================================
!===============================================================================================
!===============================================================================================

! Postroenie blijnego pola


SUBROUTINE eigenfields (kapa, gama, N, A)

USE init_data_class, ONLY : PI, delta

IMPLICIT NONE

REAL(8), INTENT(IN)										:: kapa, gama

INTEGER(4), INTENT(IN)									:: N

COMPLEX(8), INTENT(IN), DIMENSION (0:4*N-1, 0:4*N-1)	:: A

COMPLEX(8), DIMENSION (0:4*N-1)							:: z


REAL(8)													:: x1, y1, x2, y2, x3, a0, b0

REAL(8)													:: r1, r2, u1, u2, t, t1
REAL(8)													:: E_0, E_max

COMPLEX(8)												:: u

COMPLEX(8), EXTERNAL									:: Ez_in, Ez_out

REAL(8), EXTERNAL										:: x, y, directivity

!************************************************************************************

CALL eigenvector (kapa, gama, N, A, z) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

!=============================================================================================

!postroenie normirovannogo blijnego pola

OPEN (unit=7, file='H_FIELD.dat')

a0 = -1.5D0

b0 = 2.5D0


		DO x2 = a0, b0, 0.01				! setka v dekartovoi sisteme koordinat
			DO y2 = -b0, b0, 0.01

				t = datan2(y2, x2)

				r2 = dsqrt(x2**2 + y2**2)		! koordinati tochki (x2, y2) v polyarnoi sisteme koordinat

				IF ((t >= -PI) .and. (t <= 0) ) t = t + 2.0D0*PI


				x1 = x(t)					! koordinati tochki na krivoi pri ugle t
				y1 = y(t)
		
	           	r1 = dsqrt(x1**2 + y1**2)

			
				IF ( r2 <= r1 )  u = Ez_in (kapa, gama, N, z, x2, y2)  !Pole vnutri resonatora

				IF ( r2 > r1 )	u = Ez_out (kapa, gama, N, z, x2, y2) !Pole vne resonatora

!				u1=dreal(u)
!				u2=aimag(u)

!				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(x2), sngl(y2), sign(cdabs(u), datan2(u2, u1))

				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), x2, y2, cdabs(u)


			END DO
		END DO

close (7)


!-----------------------------------------------------------------

CALL far_field_pattern (kapa, gama, N, z)
		
print*, "delta = ", delta, "   ",  "directivity = ", directivity (kapa, gama, N, z)


!-----------------------------------------------------------------

!Postroenie kontura resonatora

OPEN (unit=2, file='CONTUR.dat')

	DO t = 0, 2.0D0*PI, 0.01

		WRITE(2,'(2X,E13.6, 3X,E13.6 )'), x(t), y(t)

	END DO

close (2)

!-----------------------------------------------------------------


END SUBROUTINE eigenfields

!===============================================================================================
!===============================================================================================
!===============================================================================================

! Postroenie blijnego pola


SUBROUTINE eigenfields_kite (kapa, gama, N, A)

USE init_data_class, ONLY : PI, delta

IMPLICIT NONE

REAL(8), INTENT(IN)										:: kapa, gama

INTEGER(4), INTENT(IN)									:: N

COMPLEX(8), INTENT(IN), DIMENSION (0:4*N-1, 0:4*N-1)	:: A

COMPLEX(8), DIMENSION (0:4*N-1)							:: z

INTEGER(4)												:: i, j

REAL(8)													:: x1, y1, x2, y2, x3, a0, b0, x4, x5

REAL(8)													:: r1, r2, u1, u2, t, t1, t2, t3, t4

REAL(8)													:: eps = 1.0D-2

REAL(8)													:: eps1, eps2


REAL(8)													:: y_plus, y_minus

COMPLEX(8)												:: u

COMPLEX(8), EXTERNAL									:: Ez_in, Ez_out

REAL(8), EXTERNAL										:: x, y, directivity

!************************************************************************************

CALL eigenvector (kapa, gama, N, A, z) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

!=============================================================================================

!postroenie normirovannogo blijnego pola

OPEN (unit=7, file='H_FIELD.dat')

a0 = -7.0D0

b0 = 2.0D0


	DO x2 = a0, b0, 0.01				! setka v dekartovoi sisteme koordinat

		DO y2 = -b0, b0, 0.01

		!---------------------------------------------------------------------------------


			y_plus = y2*(1.0D0 + eps)

			y_minus = y2*(1.0D0 - eps)


			eps1 = eps

			eps2 = eps



			IF ((y2 < -(1.0D0 + eps1)) .or. (y2 > (1.0D0 + eps1) ) ) THEN

			
				u = Ez_out (kapa, gama, N, z, x2, y2) !Pole vne resonatora


			ELSE


					t1 = dasin(y2/(1.0D0 + eps1))				
					t2 = PI - dasin(y2/(1.0D0 + eps1))


				! koordinati tochki na krivoi pri parametre t1, t2
						
					x1 = (1.0D0 + eps1)*x(t1)
					x3 = (1.0D0 + eps1)*x(t2)

					i = 0

					IF ( x1  < x2)   i = i+1
					IF ( x3  < x2)   i = i+1

					IF (mod(i, 2) .eq. 0) u = Ez_out (kapa, gama, N, z, x2, y2) !Pole vne resonatora

					
					IF (mod(i, 2) .ne. 0) THEN


					!----------------------------------------------------------
	
						IF ( (y2 > - (1.0D0 - eps2) )  .and. ( y2 < ( 1.0D0 - eps2) ) ) THEN


							t3 = dasin(y2/(1.0D0 - eps2))
							t4 = PI - dasin(y2/(1.0D0 - eps2))

						
							! koordinati tochki na krivoi pri parametre t3, t4
						
							x4 = (1.0D0 - eps2)*x(t3)
							x5 = (1.0D0 - eps2)*x(t4)

	           		
							! esli sleva ot tochki (x2, y2) pryamaya y = y2 peresekaet 
				
							! kontur chetnoe chislo raz, to	tochka (x2, y2) lejit vne kontura

							j = 0
		
							IF ( x4  < x2)   j = j+1
		
							IF ( x5  < x2)   j = j+1

							IF (mod(j, 2) .ne. 0) 	u = Ez_in (kapa, gama, N, z, x2, y2)  !Pole vnutri resonatora	
							
							IF (mod (j, 2) .eq. 0)  u = 0.0	


						END IF
					
					END IF	
								



				END IF
!---------------------------------------------------------------------------


!				u1=dreal(u)
!				u2=aimag(u)

!				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(x2), sngl(y2), sign(cdabs(u), datan2(u2, u1))

				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), x2, y2, cdabs(u)


			END DO
		END DO

close (7)


!-----------------------------------------------------------------

CALL far_field_pattern (kapa, gama, N, z)
		
print*, "delta = ", delta, "   ",  "directivity = ", directivity (kapa, gama, N, z)


!-----------------------------------------------------------------

!Postroenie kontura resonatora

OPEN (unit=2, file='CONTUR.dat')

	DO t = 0, 2.0D0*PI, 0.01

		WRITE(2,'(2X,E13.6, 3X,E13.6 )'), x(t), y(t)

	END DO

close (2)

!-----------------------------------------------------------------


END SUBROUTINE eigenfields_kite


!===============================================================================================

!===============================================================================================


FUNCTION det0 (A, N) RESULT (f) 

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
					
END FUNCTION det0	
