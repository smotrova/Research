SUBROUTINE GetRoot_Powell(kapa0, gama0, kapa, gama)

USE init_data_class, ONLY : alfa, PI, im, zeta, N

IMPLICIT NONE

REAL(8), INTENT(IN)		:: kapa0, gama0

REAL(8), INTENT(OUT)	:: kapa, gama

INTEGER(4)				:: ITMAX = 200

INTEGER(4), PARAMETER	:: NN = 2

REAL(8)					:: ERRREL = 1.0D-7

REAL(8)					:: X(NN), XGUESS(NN), FNORM

REAL(8), EXTERNAL		:: FCN, LSJAC

!============================================================================================

	XGUESS(1) = kapa0
	XGUESS(2) = gama0

	CALL DNEQNF (FCN, ERRREL, NN, ITMAX, XGUESS, X, FNORM)

!	CALL DNEQNJ (FCN, LSJAC, ERRREL, NN, ITMAX, XGUESS, X, FNORM)

	kapa = x(1)
	gama = x(2)

	print*,''
	print*,'kapa=', kapa
	print*,'gama=',gama
	print*,'F(kapa, gama)=', dsqrt(dabs(fnorm))
	print*,''
	print*,'*****************************************************************************'
	print*,''

END SUBROUTINE GetRoot_Powell

!###############################################################################################

SUBROUTINE GetRoot_Newton (kapa0, gama0, kapa, gama)

USE init_data_class, ONLY : alfa, PI, im, zeta, N

IMPLICIT NONE

REAL(8), INTENT(IN)  :: kapa0, gama0

REAL(8), INTENT(OUT) :: kapa, gama

REAL(8) :: x(2)

LOGICAl :: error 

!============================================================================================
		
	CALL Newton (x, error, kapa0, gama0) 

	print*,''
	print*,'kapa=', x(1)
	print*,'gama=', x(2)
	print*,''
	print*,''


		IF (error) THEN

			print*,  'error from Newton' 

		ELSE

			kapa=x(1)
			gama=x(2)


		END IF

	END SUBROUTINE GetRoot_Newton

!###############################################################################################
