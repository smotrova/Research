!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-----------------------Podprograma dla nah, korney nelineynoy----------------------------

!-------------------------------sistemi  metodom Hewtona----------------------------------

!--------------------------------s zadanoy tochnostu eps----------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE Newton2 (x, n, error, kapa0, gama0, alfa ) 
implicit none


INTEGER      n							

REAL(8)      v ,v1, v2 
REAL(8)      kapa1, kapa2, gama1, gama2, kapa00, gama00

REAL(8)      kapa0, gama0				!Nachalnoe ptiblijenie
REAL(8)      alfa, alfa1, alfa2, alfa01, alfa02
REAL(8)		 alfa_eff	

REAL(8)		 eps						!Accuracy

REAL(8) ::   h1=0.10D-5
REAL(8) ::   h2=0.00000001

REAL(8)      A(2,2), b(2)
REAL(8)      x(2)						!Reshenie systemi


REAL(8)		 f1, f2, f3, f4					! f1, f2 - levaya chast' uravneniya
											! f3, f4 - derivatives on alfa

LOGICAL      ERROR							

EXTERNAL	 f1
EXTERNAL	 f2
EXTERNAL     alfa_eff 


COMMON /Accuracy/ eps 


!*****************************************************************************************
	gama00=gama0
	kapa00=kapa0
	v=1
	DO WHILE  (v>=eps)  

		gama1=gama00+h2/2
		gama2=gama00-h2/2

		kapa1=kapa00+h1/2
		kapa2=kapa00-h1/2

		alfa01=alfa_eff(kapa00)+h1/2
		alfa02=alfa_eff(kapa00)-h1/2

		alfa1=alfa_eff(kapa1)
		alfa2=alfa_eff(kapa2)

		alfa=alfa_eff(kapa00)


		f3=(f1(n, kapa00, gama00, alfa01)-f1(n, kapa00, gama00, alfa02))/h1
		f4=(f2(n, kapa00, gama00, alfa01)-f2(n, kapa00, gama00, alfa02))/h1

		A(1,1)=(f1(n, kapa1, gama00, alfa)-f1(n, kapa2, gama00, alfa))/h1 + (alfa1-alfa2)*f3/h1

		A(1,2)=(f1(n, kapa00, gama1, alfa)-f1(n, kapa00, gama2, alfa))/h2

		A(2,1)=(f2(n, kapa1, gama00, alfa)-f2(n, kapa2, gama00, alfa))/h1 + f4*(alfa1-alfa2)/h1

		A(2,2)=(f2(n, kapa00, gama1, alfa)-f2(n, kapa00, gama2, alfa))/h2

		b(1)=-f1(n,kapa00,gama00, alfa) + a(1,1)*kapa00 + a(1,2)*gama00
		b(2)=-f2(n,kapa00,gama00, alfa) + a(2,1)*kapa00 + a(2,2)*gama00

		CALL GAUSS (2, A, b, x, error) 

		v1=dabs(x(1)-kapa00)
		v2=dabs(x(2)-gama00)
		v=max(v1,v2)
		
		IF (x(1) < 0 ) THEN
			EXIT
		ELSE

			kapa00=x(1)
			gama00=x(2)

		END IF

	END DO

END SUBROUTINE 	


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-----------------------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++