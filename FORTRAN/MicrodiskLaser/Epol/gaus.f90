!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------Podprogramma dla reshenia SLAU metodom Gaussa ----------------------

!----------------------- s viborom glavnogo elementa po stroke ----------------------------

!--------------------------------  Ax=y, n- rank  -----------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE GAUSS (n, A, b, x,error) 
implicit none

INTEGER i, j, i1, k, l
INTEGER n                 !RANK mat A,  Ax=b
REAL(8) A(n,n)			  !Matrix of koef. 	
REAL(8) b(n), x(n)		  
REAL(8) ab, big, hold,t,SUM
LOGICAL ERROR

!********************************************************************************************

ERROR = .FALSE.
 
 DO i=1, n-1	!Procedura nah. max elementa
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

!*****************************************************************************************

	IF (big==0.0D0) THEN 
		ERROR = .TRUE.
		ELSE 
			IF (l /= i) THEN	! obmen strok
				DO  j=1, n
					hold=a(l,j)
					a(l,j)=a(i,j)
					a(i,j)=hold
				END DO
				hold=b(l)
				b(l)=b(i)
				b(i)=hold
			END IF

!****************************************************************************************

		
			DO j=i1, n					! Privedenie k treugol'nomu vidu
				t=a(j,i)/a(i,i)
				DO k=i1, n
					a(j,k)=a(j,k)-t*a(i,k)
				ENDDO
				b(j)=b(j)-t*b(i)
			ENDDO
	END IF 
ENDDO


!*****************************************************************************************

IF (a(n,n)==0.0D0) THEN					!Nahojdenie korney

	ERROR=.TRUE.
	ELSE 
		x(n)=b(n)/a(n,n)
		i=n-1
		DO WHILE (i/=0)   !UNTIL (i==0.0D0) 
			SUM =0.0D0
			DO j=i+1, n
				SUM=SUM+a(i,j)*x(j)
			ENDDO
			x(i)=(b(i)-SUM)/a(i,i)
			i=i-1	
		ENDDO
	ENDIF

!******************************************************************************************

END SUBROUTINE 		

