Module init_data_class

	
	REAL(8), PARAMETER		:: alfa_cav = 2.630D0		!refractive index of a microcavity

	REAL(8), PARAMETER		:: alfa_h =  2.630D0		!refractive index of a high-contrast layer of a Bragg reflector

	REAL(8), PARAMETER		:: alfa_l = 1.0D0			!refractive index of a low-contrast layer of a Bragg reflector


	REAL(8)					:: dc = 1					!duty cycle = zeta_h/zeta_l

	REAL(8)					:: zeta_l 					!w_l/a1 relative width of low-contrast layer of a Bragg reflector
														! zeta_l = zeta_h/dc
	
	REAL(8)					:: zeta_h = 0.20D0			!w_h/a1 relative width of high-contrast layer of a Bragg reflector
													

	REAL(8)					:: delta = 0.10D0			!relative (/a1) distance between cavity and Bragg reflector
														! a1 is a radius of microcavity


	REAL(8), PARAMETER		:: PI = 3.141592653589D0
	COMPLEX(8), PARAMETER	:: im=(0.0D0, 1.0D0)


	INTEGER(4)					:: n = 7				!azimuth mode index

	INTEGER(4), PARAMETER	  	:: M = 3        		!a number of layers in annular Bragg reflector

	REAL(8), DIMENSION (1:M+1)	:: alfa

	REAL(8), DIMENSION (1:M)	:: ro

!_____________________________________________

CONTAINS

		SUBROUTINE Initialization (alfa, ro, M)

		IMPLICIT NONE

			INTEGER(4) :: i
			INTEGER(4) :: M 	!a number of layers in annular Bragg reflector

			REAL(8), DIMENSION (1:M+1)	:: alfa
			REAL(8), DIMENSION (1:M)	:: ro


!------------------------------------------------------------------------------
				
				zeta_l = zeta_h/dc

				ro(1) = 1.0D0

				IF (M > 1) THEN
					ro(2) = ro(1) + delta 

					DO i= 3, M

						IF (mod (i, 2) .ne. 0) ro(i) = ro(i-1) + zeta_h
						IF (mod (i, 2) == 0)   ro(i) = ro(i-1) + zeta_l

					END DO 

				END IF	

!======================================
!Dlya krugovogo usileniya

!				zeta_l = zeta_h/dc

!				ro(M) = 1.0D0

!				IF (M > 1) THEN

!					DO i= M, 2, -1

!						IF (mod (i, 2) .ne. 0) ro(i-1) = ro(i) - zeta_h
!						IF (mod (i, 2) == 0)   ro(i-1) = ro(i) - zeta_l

!					END DO 

!					ro(1) = ro(2) - delta 

!				END IF	

!---------------------------------------------------------------------------------

				alfa(1) = alfa_cav

				DO i = 2, M+1

					IF (mod (i, 2) == 0) alfa(i) = alfa_l

					IF (mod (i, 2) .ne. 0) alfa(i) = alfa_h

				END DO 

				alfa(M+1) = 1.0D0	 ! vne Bragg reflectora svobodnoe prostranstvo

		END SUBROUTINE Initialization

END MODULE init_data_class

