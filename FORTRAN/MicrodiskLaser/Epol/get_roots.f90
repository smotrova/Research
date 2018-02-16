!************************************************************************************************

!------------------------ Program for finding pair (kapa; gama)---------------------------------- 
!-----------------------------with effective index method----------------------------------------

!*************************************************************************************************

SUBROUTINE get_roots (x, n, kapa0, alfa ) 
Implicit none


INTEGER      m,  i
INTEGER 	 n

REAL(8) ::   PI=3.1415926535890D0

REAL(8)      alfa				!Refraction index

REAL(8)      gama0, kapa0, kapa00, l, kk

REAL(8)		 alfa_eff			!Effective refraction index

REAL(8)		 alfa1

logical      error 

REAL(8)      x(2), v

REAL(8)      eps, eps1

EXTERNAL alfa_eff

COMMON /Accuracy/ eps 

!***********************************************************************************************

	eps1=eps/100

	kapa00=kapa0

!________________________Initial iteration with alfa=const________________________________________
		
	IF (kapa0 > 0 .and. kapa0 < n) THEN
		gama0=0
	END IF

	IF (kapa0 >= n) THEN
		gama0=dlog((alfa-1)/(alfa+1))/(2*kapa0)	
	END IF
		
!------------------------------------------------------------------------------------------------

alfa=alfa_eff(kapa0)

v=1

DO WHILE ( v > eps1 )	

	CALL Newton2 (x, n, error, kapa0, gama0, alfa )		! Find (kapa; gama) with alfa_eff=alfa_eff(kapa)

	IF (error) THEN
		
		print*, '!!!!!!!! ERROR from Newton2 !!!!!!!!'

	ELSE

		IF ( x(1) < 0 ) THEN

			EXIT

		ELSE

			kapa0=x(1)
			gama0=x(2)

			alfa1=alfa
			alfa=alfa_eff(x(1))
	
			v=dabs(alfa1-alfa)

		END IF

	END IF

END DO

!print*, 'alfa_eff=', alfa

!************************************************************************************************


END SUBROUTINE 	
