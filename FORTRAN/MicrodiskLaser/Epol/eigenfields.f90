!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------Programa dla nah. sobstv. poley----------------------

!-----------------------------rezultati vivodatsa v file----------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 



SUBROUTINE eigenfields (n, k0)
IMPLICIT NONE 


INTEGER      m
INTEGER 	 n	

REAL(8) ::   PI=3.141592653589D0

REAL(8)      alfa, alfa1

INTEGER      i1,i0

COMPLEX(8)   E, u

REAL(8)      kapa0, k0
logical      error 

REAL(8)      x(2),z,y,z0
REAL(8)      phi, r

REAL(8)      eps

EXTERNAL	 E

COMMON /ALFA/ alfa
COMMON /Accuracy/ eps 


!******************************************************************************************

		m=int(( k0*alfa/Pi-n/2-0.250D0))	! E pol
		kapa0=Pi*(2*m+n+0.50D0)/2/alfa	

!***********************************************************************************************

		CALL get_roots (x, n, kapa0, alfa ) 

		
				Print*,''
				print*, "alfa_eff=", alfa, "         ", "n=", n
				print*,''
				print*,"kapa=",x(1),"          ","gama=",-x(2)
				print*,''
				print*,''

!********************************************************************************************
!---------------------E field in case of E polarization---------------------------------------

OPEN (unit=7, file='E_FIELD.dat')

		z0 = 2
		DO z=-z0, z0, 0.05
			DO y=-z0, z0, 0.05
	
			IF (z.ne.0) then	

				phi = datan(y/z)

			ENDIF

				r=dsqrt(z**2+y**2)

				u=E(n,x(1),x(2),r)*dcos(n*phi)

				WRITE(7,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), sngl(z), sngl(y), sngl(cdabs(u))	 ! writes to file

			
			END DO
		END DO
close (7)

!*****************************************************************************************
! Building of circular

OPEN (unit=5, file='CIRCLE.dat')

	DO phi=0,360

		WRITE(5,'(2X,E13.6, 3X,E13.6)'), sngl(dcosd(Phi)), sngl(dsind(phi))

	ENDDO

	close (5)

!*******************************************************************************************

print*,'The results output in files Fields.dat, Circle.dat' 

!************************************************************************************************
	
END SUBROUTINE eigenfields 









