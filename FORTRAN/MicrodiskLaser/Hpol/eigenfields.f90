!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------Programa dla nah. sobstv. poley----------------------

!-----------------------------rezultati vivodatsa v file----------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 



SUBROUTINE eigenfields (n, k0)
IMPLICIT NONE 


INTEGER      m
INTEGER 	 n	

REAL(8) ::   PI=3.141592653589D0

REAL(8)      alfa

INTEGER      i1,i0

COMPLEX(8)   Er, Ef, u1, u2
REAL(8)      u3

REAL(8)      kapa0, k0
logical      error 

REAL(8)      x(2),z,y,z0
REAL(8)      phi, r

REAL(8)      eps

EXTERNAL	 Ef, Er

COMMON /ALFA/ alfa
COMMON /Accuracy/ eps 


!******************************************************************************************

		m= int(( k0*alfa/Pi-n/2+0.250D0)) 			! H-pol
		kapa0=Pi*(2*m+n-0.50D0)/2/alfa			!H pol

!***********************************************************************************************

		CALL get_roots (x, n, kapa0, alfa ) 


		IF ( x(1)>0 .and. x(2)<0) THEN

			Print*,''
			print*, "alfa_eff=", alfa, "         ", "n=", n
			print*,''
			print*,"kapa=",x(1),"          ","gama=",-x(2)
			print*,''
			print*,''

!********************************************************************************************
!---------------------E field in case of H polarization---------------------------------------


		OPEN (unit=8, file='E_FIELD_Hpol.dat')

		z0 = 2
		DO z=-z0, z0, 0.05
			DO y=-z0, z0, 0.05
	
			IF (z.ne.0) then	

				phi = datan(y/z)

			ENDIF

				r=dsqrt(z**2+y**2)

				u1=Er(n,x(1),x(2),r)*dsin(n*phi)

				u2=Ef(n,x(1),x(2),r)*dsin(n*phi)

				u3=dsqrt(cdabs(u1)**2+cdabs(u2)**2)


				WRITE(8,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(z), sngl(y), sngl(u3)	 ! writes to file
			
			END DO
		END DO
	close (8)

!*****************************************************************************************
! Building of circular

OPEN (unit=5, file='CIRCLE.dat')

	DO phi=0,360

		WRITE(5,'(2X,E13.6, 3X,E13.6)'), sngl(dcosd(Phi)), sngl(dsind(phi))

	ENDDO

	close (5)

!*******************************************************************************************

print*,'The results output in files Fields.dat, Circle.dat' 

ELSE

	print*,'ERROR_3'
	
END IF

!************************************************************************************************
	
END SUBROUTINE eigenfields 









