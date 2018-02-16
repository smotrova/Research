!H pol

SUBROUTINE kernels (kapa, gama, t1, t2, A1, A2, B1, B2, C1, C2, D1, D2)

USE init_data_class
USE DFIMSL

IMPLICIT NONE


REAL(8), INTENT(IN)			:: t1, t2, kapa, gama

COMPLEX(8),INTENT(OUT)		:: A1, A2, B1, B2, C1, C2, D1, D2

COMPLEX(8)					:: nu, kapa_i, kapa_iR

COMPLEX(8)					:: kapa_e, kapa_eR

REAL(8)						:: f0, f1, f2, f3, f4, f5, u, v

COMPLEX(8), ALLOCATABLE		:: CJ(:), CY(:)
COMPLEX(8), ALLOCATABLE		:: J(:), Y(:)

INTEGER(4)					:: Q, CQ, NK, ERR1, ERR2

COMPLEX(8)					:: eta_i, eta_e

INTEGER(4), EXTERNAL		:: QMIN, CQMIN

REAL(8), EXTERNAL			:: R_n1, R_n2, n1_n2, R, curv, L


!-----------------------------------------------------
	
	nu = DCMPLX(alfa_i, -gama)

	eta_i = 1.0D0/(nu*nu)						! H pol
	eta_e = 1.0D0/(alfa_e*alfa_e) + im*0.0D0


!eta_i = (1.0D0, 0.0D0)
!eta_e = (1.0D0, 0.0D0)


!-----------------------------


	f0 = R(t1, t2)

	f1 = R_n1(t1, t2)/f0

	f2 = R_n2(t1, t2)/f0

	f3 = n1_n2(t1, t2)/f0
	
	
	kapa_i = kapa*nu
	kapa_e = DCMPLX(kapa*alfa_e, 0.0D0)


	kapa_iR = kapa_i*f0
	kapa_eR = kapa_e*f0


	IF (t1 .ne. t2) THEN

		NK=INT(ceiling(cdabs(kapa_iR)) + 15)
	    CQ = CQMIN(kapa_iR,NK+1)
		ALLOCATE(CY(-1:NK+1), CJ(-1:CQ+1), STAT=ERR1)
        
		CALL CDJBES(CQ,kapa_iR,CJ)                  
        CALL CDYBES(kapa_iR,CJ,CY,NK,CQ)   

!------ Complex kapa_e

		NK=INT(ceiling(cdabs(kapa_eR)) + 15)
	    Q = CQMIN(kapa_eR,NK+1)

		ALLOCATE(Y(-1:NK+1), J(-1:Q+1), STAT=ERR2)
       
		CALL CDJBES(Q,kapa_eR,J)                  
        CALL CDYBES(kapa_eR,J,Y,NK,Q)   


!--------------------------------------------------------------------------------------------------

		v = 2.0D0*dsin((t1-t2)/2.0D0)

		u = dlog(v*v)



		A1 = -(kapa_i*CJ(1) - kapa_e*J(1))*f2/(4.0D0*PI)

		B1 = -(CJ(0) - eta_i/eta_e*J(0))/(4.0D0*PI)

		C1 = ( (kapa_i*kapa_i*CJ(2) - kapa_e*kapa_e*J(2))*f1*f2 - (kapa_i*CJ(1) - kapa_e*J(1))*f3 )/(4.0D0*PI)

		D1 = (kapa_i*CJ(1) - eta_i/eta_e*kapa_e*J(1))*f1/(4.0D0*PI)



		A2 = -(im*PI + u)*A1 - (kapa_i*CY(1) - kapa_e*Y(1))*f2/4.0D0

		B2 = -(im*PI + u)*B1 - (CY(0) - eta_i/eta_e*Y(0))/4.0D0

		C2 = -(im*PI + u)*C1 + ( (kapa_i*kapa_i*CY(2) - kapa_e*kapa_e*Y(2))*f1*f2 - (kapa_i*CY(1) - kapa_e*Y(1))*f3 )/4.0D0

		D2 = -(im*PI + u)*D1 + (kapa_i*CY(1) - eta_i/eta_e*kapa_e*Y(1))*f1/4.0D0

!--------------------------------------------------------------------------------------

		deallocate (CJ,CY,stat=ERR1)
		deallocate (J,Y,stat=ERR2)


	ELSE


		f4 = curv(t1)

		f5 = L(t1)


		A1 =  (0.0D0, 0.0D0)

		B1 = - (1.0D0 - eta_i/eta_e)/(4.0D0*PI)

		C1 = - (kapa_i*kapa_i - kapa_e*kapa_e)/(8.0D0*PI)

		D1 =  (0.0D0, 0.0D0)


		A2 =  (0.0D0, 0.0D0)

		B2 =  (im*PI - 2.0D0*e_gamma)*(1.0D0 - eta_i/eta_e)/(4.0D0*PI) + (cdlog(2.0D0/kapa_i/f5) - eta_i/eta_e*cdlog(2.D0/kapa_e/f5) )/(2.0D0*PI) 

		C2 =  1.D0/8.D0/PI*( im*PI + 1.D0  - 2.D0*e_gamma)*(kapa_i*kapa_i - kapa_e*kapa_e) + 1.D0/4.D0/PI*(kapa_i*kapa_i*cdlog(2.D0/kapa_i/f5) - kapa_e*kapa_e*cdlog(2.D0/kapa_e/f5))
		
		D2 =  - (1.0D0 - eta_i/eta_e)/(4.D0*PI)*f4
	
END IF

END SUBROUTINE kernels

!======================================================================================