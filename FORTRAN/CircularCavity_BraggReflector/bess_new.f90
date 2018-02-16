! ===========================================================
! IMPROVED VERSIOM OF BESSELS by DU TOIT RECCURENT FORMULAS |
! Jn(), Yn() for COMPLEX ARGUMENT AND ARBITRARY LARGE ORDER |
! ===========================================================

      INTEGER FUNCTION QMIN(X,N)
      REAL*8 X,Z
      INTEGER N
      Z=DABS(X)
      
	  IF(Z.LE.25.0D0) THEN
       QMIN=IDNINT(Z+10.26D0*Z**0.341015D0+1.8D0+N)
      ELSE
       QMIN=IDNINT(Z+6.6362D0*Z**0.342481D0+0.4D0+N)
      END IF
      END FUNCTION QMIN


!       *******************************************
        
	INTEGER FUNCTION CQMIN(X,N)      ! N - необходимое число членов 
	IMPLICIT NONE                    ! разложения для рассчета Бесселя 
	COMPLEX*16 X                     ! от аргумента X с заданной точностью
	INTEGER N
	REAL*8 Z
        Z=CDABS(X)
    IF(Z.LE.25.0D0) THEN
	  CQMIN=IDNINT(Z+10.26D0*Z**0.341015D0+1.8D0+N)
    ELSE
	  CQMIN=IDNINT(Z+6.6362D0*Z**0.342481D0+0.4D0+N)
    END IF
	END FUNCTION CQMIN

! *******************************************
SUBROUTINE DJBES(Qm,X,J)
!  recurrence formulas by Du Toit: Jn-1(x) = 2*n/x*Jn(x) - Jn+1(x) (X-REAL*8)
INTEGER Q,N,Qm
REAL*8 X,J(-1:Qm+1),S
LOGICAL FLAG

  ARG_ZERO: IF (X.ne.0.0D0) THEN	! IF |x|=0 then Jn=1 (n=0); Jn=0 (n>0)

	Q=Qm

11	FLAG=.TRUE.						! .FALSE. IF OVERFLOW ACCUARS
    J(Q+1)=0.0D0
    J(Q)=1.0D-70

    DO N=Q,1,-1
       J(N-1)=2.0D0*N/X*J(N)-J(N+1)
! /----------
		Jn_OVERFLOW: IF (J(N-1)>1.0D+220) THEN
			FLAG=.FALSE.
			Q=Q-20
			EXIT
		ENDIF Jn_OVERFLOW
! \----------
    ENDDO      
! /----------
	OVERFLOW_LOOP: IF (.NOT.FLAG) THEN 
		GOTO 11
	ENDIF OVERFLOW_LOOP
! \----------
    S=J(0)
    DO N=1,INT(Q/2)
       S=S+2.0D0*J(2*N)
    END DO  

    DO N=0,Q+1
       J(N)=J(N)/S
    END DO
    J(-1)=-J(1)

	DO N=Q, Qm		! set all Jn() for n>Q = Jq(), q - max order of the Jn() which can be calculated
	   J(N+1)=J(Q)
	ENDDO
  ELSE 
    J=0.0D0			! all elements, except J(0), are set to be equal zero
    J(0)=1.0D0
  ENDIF ARG_ZERO

END SUBROUTINE DJBES


! *******************************************
SUBROUTINE CDJBES(Qm,X,J)        ! подпро-ма рассчета Бесселя компл арг.
IMPLICIT NONE
!   recurrence relation         Jn-1(x) = 2*n/x*Jn(x) - Jn+1(x)
INTEGER Q,N,Qm
LOGICAL FLAG
COMPLEX*16 X,J(-1:Qm+1),S,CZ

if (cdabs(x).ne.0.0D0) then			! if |x|=0 then Jn=1 (n=0); Jn=0 (n>0)

	Q=Qm
	CZ=CDCOS(X)

11	FLAG=.TRUE.					! FALSE IF OVERFLOW ACCUARS
	J(Q+1)=(0.0D0,0.0D0)
	J(Q)=(1.0D-70,0.0D0)
    
	DO N=Q,1,-1
         J(N-1)=2.0D0*N/X*J(N)-J(N+1)
! /----------
		 Jn_OVERFLOW: IF (CDABS(J(N-1))>1.0D+220) THEN
			FLAG=.FALSE.
			Q=Q-20
			EXIT
		ENDIF Jn_OVERFLOW
! \----------
    ENDDO      
    
! /----------
	OVERFLOW_LOOP: IF (.NOT.FLAG) THEN 
		GOTO 11
	ENDIF OVERFLOW_LOOP
! \----------

	IF(DIMAG(X).LT.1.0D0) THEN
        S=J(0)
        DO N=1,INT(Q/2)
           S=S+2.0D0*J(2*N)
        END DO  
        DO N=0,Q+1
           J(N)=J(N)/S
        END DO
    ELSE
        S=J(0)
        DO N=1,INT(Q/2)
               S=S+2.0D0*(-1)**N*J(2*N)
        END DO       
        DO N=0,Q+1
               J(N)=J(N)*CZ/S
        END DO 
    END IF 

	DO N=Q, Qm		! set all Jn() for n>Q = Jq(), q - max order of the Jn() which can be calculated
		J(N+1)=J(Q)
	ENDDO

else 
    J=0.0D0			! all elements, except J(0), are set to be equal zero
    J(0)=1.0D0
endif 
    J(-1)=-J(1)
END SUBROUTINE CDJBES 

! *******************************************
      SUBROUTINE DYBES(X,JJ,YY,N,Q)
!   recurrence relation         Yn-1(x) = 2*n/x*Yn(x) - Yn+1(x)      
      IMPLICIT NONE
      INTEGER I,Q,N,NQ
      REAL*8 PI,TAU,X,S0,S1,SUM0,SUM1,Z1,Z2
      REAL*8, INTENT(IN)::JJ(-1:Q+1)       
      REAL*8, INTENT(OUT)::YY(-1:N)
	  INTENT (IN) :: X,Q,N

      PI=3.14159265358979D0
      TAU=0.57721566490153286D0

      SUM0=0.0D0
      SUM1=0.0D0
      I=1
  1   IF((1+2*I).GE.(Q+1)) THEN
         PRINT*,'OVERFLOW'
         STOP
      END IF
      S1=(-1.0D0)**I*(1.0D0+2.0D0*I)*JJ(1+2*I)/(I*(1.0D0+I))
      S0=(-1.0D0)**I*JJ(2*I)/DBLE(I)
      SUM1=SUM1+S1
      SUM0=SUM0+S0
      I=I+1
      IF(DABS(S0).GE.1.0D-15.OR.DABS(S1).GE.1.0D-15) GOTO 1
      YY(1)=2.0D0/PI*((DLOG(X/2.0D0)+TAU-1.0D0)*JJ(1)-JJ(0)/X-SUM1)
      YY(0)=2.0D0/PI*((DLOG(X/2.0D0)+TAU)*JJ(0)-2.0D0*SUM0) 
	
! /------------------
!  check if Jn() had an overflow and its tail Jn|n>NQ had been set =Jn(NQ)
	NQ=Q
	Z1=JJ(Q)
	Z2=JJ(Q-1)
	DO WHILE (Z1==Z2) 
		NQ=NQ-1
		Z1=Z2
		Z2=JJ(NQ-1)
	ENDDO
! \------------------
      DO I=1,N-1
         IF (I>NQ) EXIT				! prevents overflow based on the Jn() values
		 YY(I+1)=2.0D0*I/X*YY(I)-YY(I-1)
      END DO

! /------------------
	IF (NQ<N) YY(NQ+1:N)=YY(NQ)
! \------------------

      YY(-1)=-YY(1)
      END

! *******************************************

SUBROUTINE CDYBES(X,J,Y,NK,Q)
! recurrence relation         Yn-1(x) = 2*n/x*Yn(x) - Yn+1(x)      
	IMPLICIT NONE
	REAL*8 PI,TAU
	INTEGER I,Q,NK,R,NQ
	COMPLEX*16 Y(-1:NK),J(-1:Q+1),X,S0,S1,SUM0,SUM1,Y1,Y0,P11,P21,P12,P22
	COMPLEX*16 Z1,Z2

	INTENT (IN) :: X,J,NK,Q

      PI=3.14159265358979D0 
      TAU=0.57721566490153286D0

      SUM0=(0.0D0,0.0D0)
      SUM1=(0.0D0,0.0D0)
      I=1
  1   IF((1+2*I).GE.(Q+1)) THEN
         PRINT*,'OVERFLOW'
         STOP
      END IF       
      S1=(-1.0D0)**I*(1.0D0+2.0D0*I)*J(1+2*I)/(I*(1.0D0+I))
      S0=(-1.0D0)**I*J(2*I)/DBLE(I)
      SUM1=SUM1+S1
      SUM0=SUM0+S0
      I=I+1                             
      IF(CDABS(S0).GE.1.0D-15.OR.CDABS(S1).GE.1.0D-15) GOTO 1
      Y1=2.0D0/PI*((CDLOG(X/2.0D0)+TAU-1.0D0)*J(1)-J(0)/X-SUM1)
      Y0=2.0D0/PI*((CDLOG(X/2.0D0)+TAU)*J(0)-2.0D0*SUM0)

! /------------------
!  check if Jn() had an overflow and its tail Jn|n>NQ had been set =Jn(NQ)
	NQ=Q
	Z1=J(Q)
	Z2=J(Q-1)
	DO WHILE (Z1==Z2) 
		NQ=NQ-1
		Z1=Z2
		Z2=J(NQ-1)
	ENDDO
! \------------------
        
      R=IDNINT(CDABS(X)+DABS(DIMAG(X)))
      IF(R.GE.NK) THEN
        Y(0)=Y0
        Y(1)=Y1
         DO I=1,NK-1
              Y(I+1)=2.0D0*I/X*Y(I)-Y(I-1)
         END DO    
      ELSE   
         Y(R)=(0.0D0,0.0D0)
         Y(R+1)=(1.0D0,0.0D0)
         DO I=R,1,-1
              Y(I-1)=2.0D0*I/X*Y(I)-Y(I+1)
         END DO
         P12=Y(0)
         P22=Y(1)
         Y(R)=(1.0D0,0.0D0)
         Y(R+1)=(0.0D0,0.0D0)
         DO I=R,1,-1
              Y(I-1)=2.0D0*I/X*Y(I)-Y(I+1)
         END DO
         P11=Y(0)
         P21=Y(1)
         IF(CDABS(J(1)).LE.CDABS(J(0))) THEN
              Y(R+1)=1/J(0)*(J(R+1)*Y0-2.0D0*P11/PI/X)
              Y(R)=1/J(0)*(J(R)*Y0+2.0D0*P12/PI/X)
         ELSE
              Y(R+1)=1/J(1)*(J(R+1)*Y1-2.0D0*P21/PI/X)
              Y(R)=1/J(1)*(J(R)*Y1+2.0D0*P22/PI/X)
         END IF 

         DO I=R,1,-1
              Y(I-1)=2.0D0*I/X*Y(I)-Y(I+1)
         END DO
      
         DO I=R+1,NK-1
              Y(I+1)=2.0D0*I/X*Y(I)-Y(I-1)
         END DO

! /------------------
	IF (NQ<NK) Y(NQ+1:NK)=Y(NQ)
! \------------------

      END IF
        Y(-1)=-Y(1)
END SUBROUTINE CDYBES

! =*=*=   =*=*=   =*=*=   =*=*=   =*=*=   =*=*=   =*=*=   =*=*=   =*=*=   

	COMPLEX*16 FUNCTION HNK(n,arg)				!	HANKEL FUNCTION of COMPLEX ARGUMENT
	IMPLICIT NONE
	
	INTEGER		n,Q,ERR,CQMIN,NK
	COMPLEX*16	arg
	COMPLEX*16, PARAMETER::		IM=(0.0D0,1.0D0)
	COMPLEX*16, ALLOCATABLE::	J(:),Y(:)
	EXTERNAL	CQMIN

		NK=INT(ceiling(cdabs(arg)) + 15)
        Q = CQMIN(arg,NK+1)

		ALLOCATE(Y(-1:NK+1), J(-1:Q+1), STAT=ERR)
        CALL CDJBES(Q,arg,J)                  
        CALL CDYBES(arg,J,Y,NK,Q)   

		HNK = J(n) + IM*Y(n)

		deallocate (J,Y,stat=ERR)
	END FUNCTION

