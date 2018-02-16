!============================================================================================
!============================================================================================

!	Podprogramma dla nahojdeniya korney metodom Newtona systemi uravnenii rg=2

!============================================================================================
!============================================================================================

SUBROUTINE Newton (x, error, kapa0, gama0) 

USE init_data_class

IMPLICIT NONE

REAL(8), PARAMETER		:: e_newt = 1.0D-7						!Accuracy of Newton method

INTEGER(4), PARAMETER	:: nn = 2								! rg systemi transcendentnih uravnenii

REAL(8), INTENT (IN)	:: kapa0, gama0						 	!Initial guess
	
REAL(8), INTENT (OUT)	:: x(nn)								!Solution

LOGICAL, INTENT(OUT)	:: ERROR

INTEGER(4)				:: IPATH, LDA

LOGICAL					::  ERROR1

REAL(8), PARAMETER, DIMENSION(1:nn) :: h=(/1.D-5, 1.D-6/)		! Shag po peremennoi kapa, gama dla	chislennogo nahojdeniya pervih proizvodnih po kapa, gama sootvetstvenno

REAL(8)				::  v ,v1, v2 

REAL(8)				:: kapa1, kapa2, gama1, gama2, kapa00, gama00
						
REAL(8)				::  A(nn,nn), b(nn)

COMPLEX(8)			:: res

REAL(8), EXTERNAL	:: f1, f2

!=============================================================================================
				error = .FALSE.

				gama00 = gama0
				kapa00 = kapa0

				v = 1
				
				DO WHILE  ( v > e_newt )  

					gama1 = gama00 + h(2)
					gama2 = gama00 - h(2)

					kapa1 = kapa00 + h(1)
					kapa2 = kapa00 - h(1)

					A(1,1) = (f1(kapa1, gama00, N) - f1(kapa2, gama00, N))/(h(1)*2.0)
					A(1,2) = (f1(kapa00, gama1, N) - f1(kapa00, gama2, N))/(h(2)*2.0)
					A(2,1) = (f2(kapa1, gama00, N) - f2(kapa2, gama00, N))/(h(1)*2.0)
					A(2,2) = (f2(kapa00, gama1, N) - f2(kapa00, gama2, N))/(h(2)*2.0)

					b(1) = -f1(kapa00,gama00, N) + a(1,1)*kapa00 + a(1,2)*gama00
					b(2) = -f2(kapa00,gama00, N) + a(2,1)*kapa00 + a(2,2)*gama00

!					___________________________________________________
					
					LDA = NN
					IPATH = 1

					CALL DLSARG (NN, A, LDA, B, IPATH, X)	! Reshenie systmi Ax = B

!					____________________________________________________

					!CALL GAUSS (nn, A, b, x, error1) 

!					IF (ERROR1) print*, 'ERROR FROM GAUSS'

					v1 = dabs((x(1)-kapa00)/kapa00)
					v2 = dabs((x(2)-gama00)/gama00)

					v = max(v1, v2)

					
					IF ( dabs(x(2))>1 .or. x(2) < 0.0 .or. x(1) < 0.0 ) THEN  

						EXIT
					
					ELSE

						kapa00 = x(1)
						gama00 = x(2)

					END IF

				END DO

!=============================================================================================

END SUBROUTINE 	Newton

!################################################################################################
!################################################################################################

!===============================================================================================

!----- Podprogramma dla reshenia SLAU metodom Gaussa s viborom glavnogo elementa po stroke------

!----------------------------------Ax=y, n = rank------------------------------------------------

!===============================================================================================

SUBROUTINE GAUSS (n, A, b, x, error) 

IMPLICIT NONE

INTEGER, INTENT(IN)		:: n                 !RANK mat A,  Ax=b	
REAL(8)					:: A(n,n)			  !Matrix of koef. 	
REAL(8)					:: b(n)
REAL(8), INTENT(OUT)	:: x(n)		  
LOGICAL, INTENT(OUT)	:: ERROR
		
INTEGER(4)				:: i, j, i1, k, l
REAL(8)					:: ab, big, hold,t,SUM
	
!================================================================================================

			ERROR = .FALSE.
 
			DO i=1, n-1							!Procedura nah. max elementa
				big=DABS(a(i,i))
				l=i
				i1=i+1
		
				DO j=i1, n
					ab=DABS(a(j,i))

					if (ab > big) then
						big=ab
						l=j
					END IF
				END DO


				IF (big==0.0D0) THEN 
					ERROR = .TRUE.
				ELSE 
					IF (l /= i) THEN			! obmen strok
						DO  j=1, n
							hold=a(l,j)
							a(l,j)=a(i,j)
							a(i,j)=hold
						END DO
						hold=b(l)
						b(l)=b(i)
						b(i)=hold
					END IF

		
					DO j=i1, n					! Privedenie k treugol'nomu vidu
						t=a(j,i)/a(i,i)
						DO k=i1, n
							a(j,k)=a(j,k)-t*a(i,k)
						ENDDO
						b(j)=b(j)-t*b(i)
					ENDDO
				END IF 
			ENDDO


			IF (a(n,n)==0.0D0) THEN					!Nahojdenie korney

				ERROR=.TRUE.
			ELSE 
				x(n)=b(n)/a(n,n)
				i=n-1
				DO WHILE (i/=0)						 !UNTIL (i==0.0D0) 
					SUM =0.0D0
					DO j=i+1, n
						SUM=SUM+a(i,j)*x(j)
					ENDDO
					x(i)=(b(i)-SUM)/a(i,i)
					i=i-1	
				ENDDO
			ENDIF

	END SUBROUTINE GAUSS 

!###############################################################################################
