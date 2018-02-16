		FUNCTION x (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			
			REAL(8)				:: r, e

			INTEGER(4)			:: m

			
! spiral

!			IF ( t >= beta .and. t <= 2*PI - beta)  r = 1.0D0 + delta*t/(PI*4.D0)

!			IF (t < beta )  r = 1.0D0 - delta*((2*PI-beta)/beta*t - PI/beta**2*t**2 - PI)/(PI*4.D0)

!			IF (t > 2*PI - beta)  r = 1.0D0 + delta*((2*PI - beta)/beta*(2*PI - t) - PI/beta**2*(2*PI - t)**2 + PI)/(PI*4.D0)

!			u = r*dcos(t)

!---------------------------------------------

!cycle			

!			r = 1.0D0

!			u = r*dcos(t)
!---------------------------------------------

!egg

!			r = 1.0D0

!			u = r*dcos(t)

!-------------------------------------------

!polygon
		
!			e = 1.0/3.0 

!			m = 2

!			u = dcos(t) + e*dcos(m*t)

!--------------------------------------------

!(kite-shape domain)

			u = dcos(t) + delta*dcos(2.0D0*t) - delta


!---------------------------------------------


!limacon

!			u = (1.0D0 + delta*dcos(t))*dcos(t) 

!---------------------------------------------


	
		
		END FUNCTION x


!#############################################################################


		FUNCTION y (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			
			REAL(8)				:: r,e
			INTEGER(4)			:: m


!spiral			

!			IF ( t >= beta .and. t <= 2*PI - beta)  r = 1.0D0 + delta*t/PI/4.0

!			IF (t < beta )  r = 1.0D0 - delta/PI/4.0*((2*PI-beta)/beta*t - PI/beta**2*t**2 - PI)

!			IF (t > 2*PI - beta)  r = 1.0D0 + delta/PI/4.0*((2*PI - beta)/beta*(2*PI - t) - PI/beta**2*(2*PI - t)**2 + PI)

!			u = r*dsin(t)

!-----------------------------------------------

		
!cycle			

!			r = 1.0D0
!			u = r*dsin(t)

!-----------------------------------------------
!egg

!			IF ((t>=0.0D0) .and. (t<=PI)) r = 2.0D0
			
!			IF ((t>PI) .and. (t < 2*PI))  r = 1.0D0

			
!			u = r*dsin(t)

!-----------------------------------------------

!polygon

!			e = 1.0/3.0 
						

!			m = 2

!			u = dsin(t) - e*dsin(m*t)
		
!-----------------------------------------------


!(kite-shape domain)

			u = dsin(t)
	
!-----------------------------------------------

!limacon

!			u = (1.0D0 + delta*dcos(t))*dsin(t) 

!---------------------------------------------


		
		END FUNCTION y

!###########################################################################

		FUNCTION dx_dt (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			
			REAL(8)				:: r, dr_dt

			REAL(8)				:: e

			INTEGER(4)			:: m




!spiral
!			IF ( t >= beta .and. t <= 2*PI - beta)  r = 1.0D0 + delta*t/PI/4.0

!			IF (t < beta )  r = 1.0D0 - delta/PI/4.0*((2*PI-beta)/beta*t - PI/beta**2*t**2 - PI)

!			IF (t > 2*PI - beta)  r = 1.0D0 + delta/PI/4.0*((2*PI - beta)/beta*(2*PI - t) - PI/beta**2*(2*PI - t)**2 + PI)


!			IF ( t >= beta .and. t <= 2*PI - beta)  dr_dt = delta/(4.0D0*PI)

!			IF (t < beta )  dr_dt = - delta/PI/4.0*((2*PI-beta)/beta - 2.0*PI/beta**2*t)

!			IF (t > 2*PI - beta)  dr_dt = -delta/PI/4.0*((2*PI - beta)/beta - 2*PI/beta**2*(2*PI - t))


!			u = dr_dt*dcos(t) - r*dsin(t)
			
!cycle			

!			r = 1.0D0
!			dr_dt = 0.0D0
!			u = dr_dt*dcos(t) - r*dsin(t)

!egg
!			r = 1.0D0
!			dr_dt = 0.0D0

!			u = dr_dt*dcos(t) - r*dsin(t)


!polygon
		
!			e = 1.0/3.0 

!			m = 2

!			u = -dsin(t) - m*e*dsin(m*t)


!(kite-shape domain)

			u = -dsin(t) - 2.0D0*delta*dsin(2.0D0*t)

!limacon

!			u = -dsin(t) - delta*dsin(2.0D0*t)

	
		
		END FUNCTION dx_dt

!---------------------------------------------------------------------------

		
		FUNCTION dy_dt (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			
			REAL(8)				:: r, dr_dt

			REAL(8)				:: r1, r2, t1, t2

			REAL(8)				:: e

			INTEGER(4)			:: m




!spiral
!			IF ( t >= beta .and. t <= 2*PI - beta)  r = 1.0D0 + delta*t/PI/4.0

!			IF (t < beta )  r = 1.0D0 - delta/PI/4.0*((2*PI-beta)/beta*t - PI/beta**2*t**2 - PI)

!			IF (t > 2*PI - beta)  r = 1.0D0 + delta/PI/4.0*((2*PI - beta)/beta*(2*PI - t) - PI/beta**2*(2*PI - t)**2 + PI)


!			IF ( t >= beta .and. t <= 2*PI - beta)  dr_dt = delta/PI/4.0

!			IF (t < beta )  dr_dt = - delta/PI/4.0*((2*PI-beta)/beta - 2.0*PI/beta**2*t)

!			IF (t > 2*PI - beta)  dr_dt = -delta/PI/4.0*((2*PI - beta)/beta - 2*PI/beta**2*(2*PI - t))

!			u = dr_dt*dsin(t) + r*dcos(t)

!cycle			
			
!			r = 1.0D0
!			dr_dt = 0.0D0
!			u = dr_dt*dsin(t) + r*dcos(t)
!egg

!			IF ((t>=0.0D0) .and. (t<=PI)) r = 2.0D0
			
!			IF ((t>PI) .and. (t < 2*PI))  r = 1.0D0

!			dr_dt = 0.0D0

!			u = dr_dt*dsin(t) + r*dcos(t)


!polygon
		
!			e = 1.0/3.0 
!			m = 2

!			u = dcos(t) - m*e*dcos(m*t)

!(kite-shape domain)

			u = dcos(t)

!limacon

!			u = dcos(t) + delta*dcos(2.0D0*t)

	
		
		END FUNCTION dy_dt

!---------------------------------------------------------------------------


		FUNCTION d2x_dt2 (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			
			REAL(8)				:: r, dr_dt, d2r_dt2

			REAL(8)				:: e

			INTEGER(4)			:: m



!spiral

	! r(t)
!			IF ( t >= beta .and. t <= 2*PI - beta)  r = 1.0D0 + delta*t/PI/4.0

!			IF (t < beta )  r = 1.0D0 - delta/PI/4.0*((2*PI-beta)/beta*t - PI/beta**2*t**2 - PI)

!			IF (t > 2*PI - beta)  r = 1.0D0 + delta/PI/4.0*((2*PI - beta)/beta*(2*PI - t) - PI/beta**2*(2*PI - t)**2 + PI)

	! dr_dt(t)
!			IF ( t >= beta .and. t <= 2*PI - beta)  dr_dt = delta/PI/4.0

!			IF (t < beta )  dr_dt = - delta/PI/4.0*((2*PI-beta)/beta - 2.0*PI/beta**2*t)

!			IF (t > 2*PI - beta)  dr_dt = -delta/PI/4.0*((2*PI - beta)/beta - 2*PI/beta**2*(2*PI - t))
			
	!d2r_dt2
!			IF ( t >= beta .and. t <= 2*PI - beta)  d2r_dt2 = 0.0D0

!			IF (t < beta )  d2r_dt2 =  delta/2.0/beta**2

!			IF (t > 2*PI - beta)  d2r_dt2 = - delta/2.0/beta**2

!			u = d2r_dt2*dcos(t) - 2.0*dr_dt*dsin(t) - r*dcos(t)



!cycle			

!			r = 1.0D0
!			dr_dt = 0.0D0
!			d2r_dt2 = 0.0D0
!			u = d2r_dt2*dcos(t) - 2.0*dr_dt*dsin(t) - r*dcos(t)


!egg

!			r = 1.0D0
!			dr_dt = 0.0D0
!			d2r_dt2 = 0.0D0
					
!			u = d2r_dt2*dcos(t) - 2.0*dr_dt*dsin(t) - r*dcos(t)


!polygon
		
!			e = 1.0/3.0 
!			m = 2

!			u = -dcos(t) - m*m*e*dcos(m*t)


!(kite-shape domain)			

			u = -dcos(t) - 4.0D0*delta*dcos(2.0D0*t)


!limacon

!			u = - dcos(t) - 2.0D0*delta*dcos(2.0D0*t)

	
		
		END FUNCTION d2x_dt2


!---------------------------------------------------------------------------


		FUNCTION d2y_dt2 (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			
			REAL(8)				:: r, dr_dt, d2r_dt2

			REAL(8)				:: e

			INTEGER(4)			:: m


!spiral

	! r(t)
!			IF ( t >= beta .and. t <= 2*PI - beta)  r = 1.0D0 + delta*t/PI/4.0

!			IF (t < beta )  r = 1.0D0 - delta/PI/4.0*((2*PI-beta)/beta*t - PI/beta**2*t**2 - PI)

!			IF (t > 2*PI - beta)  r = 1.0D0 + delta/PI/4.0*((2*PI - beta)/beta*(2*PI - t) - PI/beta**2*(2*PI - t)**2 + PI)

	! dr_dt(t)
!			IF ( t >= beta .and. t <= 2*PI - beta)  dr_dt = delta/PI/4.0

!			IF (t < beta )  dr_dt = - delta/PI/4.0*((2*PI-beta)/beta - 2.0*PI/beta**2*t)

!			IF (t > 2*PI - beta)  dr_dt = -delta/PI/4.0*((2*PI - beta)/beta - 2*PI/beta**2*(2*PI - t))
			
	!d2r_dt2(t)
!			IF ( t >= beta .and. t <= 2*PI - beta)  d2r_dt2 = 0.0D0

!			IF (t < beta )  d2r_dt2 =  delta/2.0/beta**2

!			IF (t > 2*PI - beta)  d2r_dt2 = - delta/2.0/beta**2


!			u = d2r_dt2*dsin(t) + 2.0*dr_dt*dcos(t) - r*dsin(t)


!cycle			

!			r = 1.0D0
!			dr_dt = 0.0D0
!			d2r_dt2 = 0.0D0
!			u = d2r_dt2*dsin(t) + 2.0*dr_dt*dcos(t) - r*dsin(t)


!egg
!			IF ((t>=0.0D0) .and. (t<=PI)) r = 2.0D0
			
!			IF ((t>PI) .and. (t < 2*PI))  r = 1.0D0

!			dr_dt = 0.0D0

!			d2r_dt2 = 0.0D0
					
!			u = d2r_dt2*dsin(t) + 2.0*dr_dt*dcos(t) - r*dsin(t)


!polygon
		
!			e = 1.0/3.0 
!			m = 2

!			u = -dsin(t) + m*m*e*dsin(m*t)



!(kite-shape domain)

			u = -dsin(t)

!limacon

!			u = -dsin(t) - 2.0D0*delta*dsin(2.0D0*t) 
	
		
		END FUNCTION d2y_dt2

!---------------------------------------------------------------------------------

		FUNCTION R (t1, t2) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t1, t2

			REAL(8)				:: u, w

			REAL(8)				:: v1, v2, z1, z2

			REAL(8), EXTERNAL	:: x, y 
			
			

			v1 = x(t1)
			v2 = x(t2)

			z1 = y(t1)
			z2 = y(t2)

			
			w = (v1 - v2)**2 + (z1 - z2)**2 

			u = dsqrt(w)
		
	
		
		END FUNCTION R


!----------------------------------------------------------------------------------



		FUNCTION L (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u, w
			REAL(8), EXTERNAL	:: dx_dt, dy_dt 
			
			w = (dx_dt(t))**2 + (dy_dt(t))**2 

			u = dsqrt(w)
	
		
		END FUNCTION L


!----------------------------------------------------------------------------------


		FUNCTION R_n1 (t1, t2) RESULT (u)		! inner product (R*n), R - rasstoyanie mejdu tochkami r(t1) i r(t2)
												! n - edinichnii vector normali v tochke t1
			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t1, t2

			REAL(8)				:: u

			REAL(8), EXTERNAL	:: x, y, dx_dt, dy_dt, L


			u = ((x(t1) - x(t2))*dy_dt(t1) - (y(t1) - y(t2))*dx_dt(t1) )/L(t1)


		END FUNCTION R_n1


!--------------------------------------------------------------------------------------

		FUNCTION R_n2 (t1, t2) RESULT (u)! inner product (R*n), R - rasstoyanie mejdu tochkami r(t1) i r(t2)
												! n - edinichnii vector normali v tochke t2
			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t1, t2
			REAL(8)				:: u
			REAL(8), EXTERNAL	:: x, y, dx_dt, dy_dt, L
			

			u =  ( (x(t1) - x(t2))*dy_dt(t2) - (y(t1) - y(t2))*dx_dt(t2))/L(t2)

	
		END FUNCTION R_n2

!-----------------------------------------------------------------------------------------

		FUNCTION n1_n2 (t1, t2) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t1, t2
			REAL(8)				:: u
			REAL(8), EXTERNAL	:: dx_dt, dy_dt, L 

			
			u = (dy_dt(t1)*dy_dt(t2) + dx_dt(t1)*dx_dt(t2))/L(t1)/L(t2)
	
		
		END FUNCTION n1_n2

!------------------------------------------------------------------------------

		FUNCTION curv (t) RESULT (u)

			USE init_data_class
			IMPLICIT NONE

			REAL(8), INTENT(IN) :: t
			REAL(8)				:: u
			REAL(8), EXTERNAL	:: d2x_dt2, d2y_dt2, dx_dt, dy_dt, L 

		
			u = - (d2x_dt2(t)*dy_dt(t) - d2y_dt2(t)*dx_dt(t))/(L(t))**3 

		END FUNCTION curv

!---------------------------------------------------------------------------------------
