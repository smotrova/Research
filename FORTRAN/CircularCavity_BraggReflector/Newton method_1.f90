!============================================================================================
!============================================================================================

!	Podprogramma dla nahojdeniya korney metodom Newtona systemi uravnenii rg=2

!============================================================================================
!============================================================================================

SUBROUTINE Newton_1 (x, kapa0, gama0) 

USE init_data_class

IMPLICIT NONE

REAL(8), PARAMETER		:: e_newt = 1.0D-7						!Accuracy of Newton method

REAL(8), INTENT (IN)	:: kapa0, gama0						 	!Initial guess
	
REAL(8), INTENT (OUT)	:: X(2)									!Solution

INTEGER(4)				:: C = 100

INTEGER(4)				:: i


REAL(8), PARAMETER, DIMENSION(1:2) :: h=(/1.0D-6, 1.0D-6/)		! Shag po peremennoi kapa, gama dla	chislennogo nahojdeniya pervih proizvodnih po kapa, gama sootvetstvenno

REAL(8)				:: v1, v2 

REAL(8)				:: kapa00, gama00
						
REAL(8)				:: JAC(1:2, 1:2)

REAL(8)				:: det_JAC, det_JAC0

REAL(8), EXTERNAL	:: f1, f2

!=============================================================================================

				gama00=gama0
				kapa00=kapa0

!normirovka 
!--------------------------------------------------------------

!		JAC(1,1) = ( -f1(kapa00+2*h(1), gama00) + 8*f1(kapa00+h(1), gama00)-8*f1(kapa00-h(1), gama00) + f1(kapa00-2*h(1), gama00) )/12/h(1)
!		JAC(1,2) = ( -f1(kapa00, gama00+2*h(2)) + 8*f1(kapa00, gama00+h(2))-8*f1(kapa00, gama00-h(2)) + f1(kapa00, gama00-2*h(2)) )/12/h(2)
!		JAC(2,1) = ( -f2(kapa00+2*h(1), gama00) + 8*f2(kapa00+h(1), gama00)-8*f2(kapa00-h(1), gama00) + f2(kapa00-2*h(1), gama00) )/12/h(1)
!		JAC(2,2) = ( -f2(kapa00, gama00+2*h(2)) + 8*f2(kapa00, gama00+h(2))-8*f2(kapa00, gama00-h(2)) + f2(kapa00, gama00-2*h(2)) )/12/h(2)


!		det_JAC0 =JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)

!--------------------------------------------------------------
				
				i = 1
				v1=1
				v2=1

				DO WHILE (v1 > e_newt .or. v2 > e_newt)

!					JAC(1,1) = ( -f1(kapa00+2*h(1), gama00)/det_JAC0 + 8*f1(kapa00+h(1), gama00)/det_JAC0-8*f1(kapa00-h(1), gama00)/det_JAC0 + f1(kapa00-2*h(1), gama00)/det_JAC0 )/12/h(1)
!					JAC(1,2) = ( -f1(kapa00, gama00+2*h(2))/det_JAC0 + 8*f1(kapa00, gama00+h(2))/det_JAC0-8*f1(kapa00, gama00-h(2))/det_JAC0 + f1(kapa00, gama00-2*h(2))/det_JAC0 )/12/h(2)
!					JAC(2,1) = ( -f2(kapa00+2*h(1), gama00)/det_JAC0 + 8*f2(kapa00+h(1), gama00)/det_JAC0-8*f2(kapa00-h(1), gama00)/det_JAC0 + f2(kapa00-2*h(1), gama00)/det_JAC0 )/12/h(1)
!					JAC(2,2) = ( -f2(kapa00, gama00+2*h(2))/det_JAC0 + 8*f2(kapa00, gama00+h(2))/det_JAC0-8*f2(kapa00, gama00-h(2))/det_JAC0 + f2(kapa00, gama00-2*h(2))/det_JAC0 )/12/h(2)


!					det_JAC = JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)

!					X(1) = kapa00 + 1/det_JAC*(f2(kapa00, gama00)/det_JAC0*JAC(1,2)-f1(kapa00, gama00)/det_JAC0*JAC(2,2))
!					X(2) = gama00 + 1/det_JAC*(f1(kapa00, gama00)/det_JAC0*JAC(2,1)-f2(kapa00, gama00)/det_JAC0*JAC(1,1))

!bez normirovki --------------------------------------------------------------------------------


					JAC(1,1) = ( -f1(kapa00+2*h(1), gama00) + 8*f1(kapa00+h(1), gama00)-8*f1(kapa00-h(1), gama00) + f1(kapa00-2*h(1), gama00) )/12/h(1)
					JAC(1,2) = ( -f1(kapa00, gama00+2*h(2)) + 8*f1(kapa00, gama00+h(2))-8*f1(kapa00, gama00-h(2)) + f1(kapa00, gama00-2*h(2)) )/12/h(2)
					JAC(2,1) = ( -f2(kapa00+2*h(1), gama00) + 8*f2(kapa00+h(1), gama00)-8*f2(kapa00-h(1), gama00) + f2(kapa00-2*h(1), gama00) )/12/h(1)
					JAC(2,2) = ( -f2(kapa00, gama00+2*h(2)) + 8*f2(kapa00, gama00+h(2))-8*f2(kapa00, gama00-h(2)) + f2(kapa00, gama00-2*h(2)) )/12/h(2)

					det_JAC = JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)

					JAC(1,1) = (f1(kapa00+h(1)/2, gama00) - f1(kapa00-h(1)/2, gama00))/h(1)
					JAC(1,2) = (f1(kapa00, gama00+h(2)/2) - f1(kapa00, gama00-h(2)/2))/h(2)
					JAC(2,1) = (f2(kapa00+h(1)/2, gama00) - f2(kapa00-h(1)/2, gama00))/h(1)
					JAC(2,2) = (f2(kapa00, gama00+h(2)/2) - f2(kapa00, gama00-h(2)/2))/h(2)


					det_JAC = JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)

					X(1) = kapa00 + 1/det_JAC*(f2(kapa00, gama00)*JAC(1,2)-f1(kapa00, gama00)*JAC(2,2))
					X(2) = gama00 + 1/det_JAC*(f1(kapa00, gama00)*JAC(2,1)-f2(kapa00, gama00)*JAC(1,1))

!------------------------------------------------------------------------------------------------

					v1=dabs((x(1)/kapa00 - 1.0D0))
					v2=dabs((x(2)/gama00 - 1.0D0))

					IF ( v1 < e_newt .and. v2 < e_newt) THEN

							EXIT
						
						ELSE

							IF (i > C) THEN

								print*,''
								print*,'NOT sucsessfull' 
								print*,'Try another initial guess'
								print*,''
								print*,''
								print*,''

								STOP
			
							ELSE
			
							   	kapa00 = x(1)
								gama00 = x(2)

								IF (x(1) < 0 .or. x(2) > 1 .or. x(2) < 0) THEN
							
									kapa00 = kapa0
									gama00 = gama0
							
								END IF

								i= i+1
			
							END IF

					END IF

				END DO

!=============================================================================================

END SUBROUTINE 	Newton_1