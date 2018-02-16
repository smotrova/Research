SUBROUTINE even_roots (a, b, x, upper_limit )

USE  DFIMSL 

IMPLICIT NONE

REAL(8)  alfa

REAL(8)  h, a, b, Am, Bm, x
REAL(8)  func_even, upper_limit

REAL(8)  EPS, ERRABS, ERREL
REAL(8)   kd						! kd=k*d, where k=wavwernumber of free space
									! 2d=slab thikness
INTEGER  MAXFN,i,m

COMMON /arg1/ kd
COMMON /ALFA1/ alfa

EXTERNAL func_even

!################################################################################################

	eps= 1.0E-6			! Accuracy of finding zeros

	m=1000				! Number of intervals  (bolshoi interval (a,b) delitsa na m malenkih
						!	dla nahojdenia kornei metodom bisekcii na kajdom malenkom intervale)
	h=(b-a)/m
	i=0

	Am=a
	Bm=a+h
!_______________________________________________________________________________________________

	DO WHILE ( i < m-1 )		! Find zero on each small interval
			
		Am=Am+h
		Bm=Bm+h

		IF ( func_even(Am)*func_even(Bm) < 0) THEN

			ERRABS= 0.0D0
			ERREL= eps
			MAXFN=200
		
			A=Am
			B=Bm
					
			CALL DZBREN(func_even, ERRABS, ERREL, A, B, MAXFN)

			IF ( DABS (func_even(B)) < 1.0D0 ) THEN
			
				IF (x < B .and. B < upper_limit) THEN	! The largest B coresponds to the first mode
					x = B
				END IF

			END IF

		END IF		

	i=i+1

	END DO

!###############################################################################################

END SUBROUTINE