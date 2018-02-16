SUBROUTINE far_field_pattern (kapa, gama, N, z)

USE init_data_class, ONLY : PI, im, alfa_e, alfa_i

IMPLICIT NONE

REAL(8), INTENT(IN)										:: kapa, gama

INTEGER(4), INTENT(IN)									:: N

COMPLEX(8), INTENT(IN), DIMENSION (0:4*N-1)				:: z

COMPLEX(8), DIMENSION (0: 2*N-1)						:: fi, psi

REAL(8)													:: teta, t_p, t

COMPLEX(8)												:: u1, u2, u, nu

COMPLEX(8)												:: eta_i, eta

REAL(8)													:: eta_e

INTEGER(4)												:: p
 
REAL(8), EXTERNAL										:: x, y, dx_dt, dy_dt, L

!************************************************************************************

nu = DCMPLX(alfa_i, -gama)

!----------------------------------------------

! H pol

eta_i = 1.0D0/(nu*nu)
eta_e = 1.0D0/(alfa_e*alfa_e)

eta = (2.0D0*eta_e)/(eta_i + eta_e)

!-----------------------------------------------


DO p = 0, 4*N-1

	IF (p <= 2*N-1) fi(p) = z(p)
	IF (p > 2*N-1)	psi(p - 2*N) = z(p)*eta

END DO

!------------------------------------------------------------------------

OPEN (unit = 19, file='FAR_FIELD.dat')

DO teta = 0, 360, 1.0D0

	u = (0.0D0, 0.0D0)

	DO  p = 0, 2*N-1

		t_p = PI*p/N

		u1 = im*fi(p)*kapa*alfa_e*(dy_dt(t_p)*dcosd(teta) - dx_dt(t_p)*dsind(teta))/L(t_p) + eta_i/eta_e*psi(p)

		u2 = CDEXP(-im*kapa*alfa_e*( x(t_p)*dcosd(teta) + y(t_p)*dsind(teta)) )


		u = u + u1*u2*L(t_p) 


	END DO

	u = PI/N*u

	WRITE(19,'(2X,E13.6, 3X,E13.6 )'), teta,  cdabs(u)

END DO	

CLOSE(19)


END subroutine far_field_pattern



!==============================================================================================
!==============================================================================================
!==============================================================================================
!==============================================================================================


FUNCTION directivity (kapa, gama, N, z) RESULT (w) 


USE init_data_class, ONLY : PI, im, alfa_e, alfa_i

IMPLICIT NONE

REAL(8), INTENT(IN)										:: kapa, gama

REAL(8)													:: w


INTEGER(4), INTENT(IN)									:: N

COMPLEX(8), INTENT(IN), DIMENSION (0:4*N-1)				:: z

COMPLEX(8), DIMENSION (0: 2*N-1)						:: fi, psi

REAL(8)													:: teta, t_p, delta_teta, teta_k, teta_max

COMPLEX(8)												:: u1, u2, u, nu

COMPLEX(8)												:: eta_i, eta

REAL(8)													:: eta_e, u_max, uu, mod_u

INTEGER(4)												:: p, M, k
 
REAL(8), EXTERNAL										:: x, y, dx_dt, dy_dt, L

COMMON /teta_napravlenie_gl_lepestka/ teta_max

!************************************************************************************

nu = DCMPLX(alfa_i, -gama)

!----------------------------------------------

! H pol

eta_i = 1.0D0/(nu*nu)
eta_e = 1.0D0/(alfa_e*alfa_e)

eta = (2.0D0*eta_e)/(eta_i + eta_e)

!-----------------------------------------------


DO p = 0, 4*N-1

	IF (p <= 2*N-1) fi(p) = z(p)
	IF (p > 2*N-1)	psi(p - 2*N) = z(p)*eta

END DO

!------------------------------------------------------------------------


uu = 0.0D0

M = 400

delta_teta = 2.0*pi/M


!DO teta = 0, 2.0D0*PI, delta_teta

DO k = 0, M-1, 1

	teta_k = k*2.0D0*PI/M

	teta = teta_k


	u = (0.0D0, 0.0D0)

	DO  p = 0, 2*N-1

		t_p = PI*p/N

		u1 = im*fi(p)*kapa*alfa_e*(dy_dt(t_p)*dcos(teta) - dx_dt(t_p)*dsin(teta))/L(t_p) + eta_i/eta_e*psi(p)

		u2 = CDEXP(-im*kapa*alfa_e*( x(t_p)*dcos(teta) + y(t_p)*dsin(teta)) )

		u = u + u1*u2*L(t_p)		! Fi(teta) - pole v dalnei zone, teta - ugol nablyudeniya

	END DO

	u = PI/N*u

!	uu = uu + cdabs(u)*cdabs(u)*delta_teta			! integral ot |Fi|^2 po uglu teta 

	uu = uu + 2.0D0*PI/M*cdabs(u)*cdabs(u)			! integral ot |Fi|^2 po uglu teta, pravilo trapezii, Fi(0) = Fi(2PI)


END DO	

!--------------------------------------------------------------------------------------------

u_max = 0.0D0
teta_max = 0.0D0

delta_teta = PI/180.0

!DO teta = 0, 2.0D0*PI, delta_teta  ! dlya proizvolnogo contura

DO teta = 0, PI, delta_teta  ! dlya contura, symmetrichnogo otnositelno osi x

	u = (0.0D0, 0.0D0)

	DO  p = 0, 2*N-1

		t_p = PI*p/N

		u1 = im*fi(p)*kapa*alfa_e*(dy_dt(t_p)*dcos(teta) - dx_dt(t_p)*dsin(teta))/L(t_p) + eta_i/eta_e*psi(p)

		u2 = CDEXP(-im*kapa*alfa_e*( x(t_p)*dcos(teta) + y(t_p)*dsin(teta)) )

		u = u + u1*u2*L(t_p)		! Fi(teta) - pole v dalnei zone, teta - ugol nablyudeniya

	END DO

	u = PI/N*u
	
	mod_u = cdabs(u)*cdabs(u)				! |Fi(teta)|^2
	
	IF (u_max < mod_u) THEN
	
		u_max = mod_u		! nahodim max value |Fi(teta)|^2 
		
		teta_max = teta		! ugol napravleniya glavnogo lepestka

	END IF

END DO

w = 2.0D0*PI*u_max/uu

END FUNCTION directivity