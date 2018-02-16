PROGRAM MAIN

USE init_data_class

IMPLICIT NONE

REAL(8)										:: kapa0, gama0
REAL(8)										:: kapa, gama, zeta_end, zeta_home
COMPLEX(8)									:: det_1
COMPLEX(8), ALLOCATABLE, DIMENSION (:, :)	:: A
INTEGER(4)									:: i, j, s
COMPLEX(8), EXTERNAL						:: det
COMPLEX(8), EXTERNAL						:: f1, f2
CHARACTER(10)								:: T1,T2, ss
REAL(8)										:: xi_end, xi_home

REAL(8)										:: x_s, y_s, fi, l, ro_1s


!-----------------------------------------------------------------------------------------------
!----------- Eff. alfa -------------

!kapa0 = 1.404D0			!(m=1, alfa=2.63, single disk)
!gama0 = 0.2750D0  

!kapa0 = 1.360D0			!(m=1, alfa=2.63, OO)
!gama0 = 0.02750D0  

!kapa0 = 2.78187D0			!(m=1, alfa=1.402, single disk)
!gama0 = 0.2971D0  

!kapa0 = 2.7D0				!(m=1, alfa=1.402, OO)
!gama0 = 0.3D0  

!kapa0 = 1.9136D0			!(m=1, alfa=1.967, single disk)
!gama0 = 0.2959D0  

!kapa0 = 1.890D0				!(m=1, alfa=1.967, OO)
!gama0 = 0.2959D0  

!kapa0 = 1.08850D0			!(m=1, alfa=3.374, single disk)
!gama0 = 0.23590D0  

!kapa0 = 1.10D0				!(m=1, alfa=3.374, EO)
!gama0 = 0.023590D0  


!-----------------------------------------------

!kapa0 = 2.3D0		!(m=3, alfa=2.63) (OO)
!gama0 = 9.086D-2   

!kapa0 = 2.26D0		!(m=3, alfa=2.63) (single disk)
!gama0 = 9.086D-2   

!________________________________________________

!kapa0 = 2.73D0		!(m=4, alfa=2.63)
!gama0 = 3.01D-2   

!________________________________________________

!kapa0 = 3.1965D0	!(m=5, alfa=2.63) (single disk)
!gama0 = 9.270D-3 

!kapa0 = 3.2248D0		!(m=5, alfa=2.63) (OqO)
!gama0 = 9.27D-3   

!_________________________________________________

!kapa0 = 3.6539D0	!(m=6, alfa=2.63) (single disk)
!gama0 = 2.785D-3  

!kapa0 = 3.68D0		!(m=6, alfa=2.63) (OqO)
!gama0 = 2.785D-3  

!kapa0 = 3.60D0		!(m=6, alfa=3.374) (OqE)
!gama0 = 2.785D-3  

!__________________________________________________

!kapa0 = 4.101380		!(m=7, alfa=2.63) (single disk)
!gama0 = 8.2988D-4   

!kapa0 = 4.070			!(m=7, alfa=2.63) (EE, identical)
!gama0 = 8.31D-4   

!kapa0 = 4.11D0			!(m=7, alfa=2.63) (OO)
!gama0 = 8.31D-4  

!kapa0 = 4.07D0			!(m=7, alfa=2.63) (OE) 
!gama0 = 8.31D-3   

!____________________________________________________

!kapa0 = 4.55		!(m=8, alfa=2.63) (OqO)
!gama0 = 2.45D-4 

!_____________________________________________________

!kapa0 = 4.98D0		!(m=9, alfa=2.63) (single disk)
!gama0 = 7.23D-5   

!kapa0 = 4.98D0		!(m=9, alfa=2.63) (OqO)
!gama0 = 7.23D-5   

!kapa0 = 4.93D0		!(m=9, alfa=2.63) (OqE)
!gama0 = 7.23D-5   

!_____________________________________________________

!kapa0 = 5.40D0		!(m=10, alfa=2.63) (single disk)
!gama0 = 2.12D-5   

!kapa0 = 5.41D0		!(m=10, alfa=2.63) (OqE)
!gama0 = 2.12D-5   

!kapa0 = 5.32D0		!(m=10, alfa=2.63) (?)
!gama0 = 2.12D-5   

!___________________________________________________________________________________________


!kapa0 = 25.63D0		!(m=33, alfa=1.5, single disk)
!gama0 = 1.9D-5  

kapa0 = 25.66D0		!(m=33, alfa=1.5) odd-odd
gama0 = 1.9D-5  


!___________________________________________________________________________________________


print*,'***********************************'
print*,' '
print*,'Hz polarization;   maximally anti-symmetric class'
print*,' '
print*,'alfa =', alfa
print*,' '
print*,'Number of resonators =', N_res
print*,' '
print*,'***********************************'
print*,' '

!____________________________________________________________________________________________

CALL DATE_AND_TIME(TIME=T1)

N = ceiling(kapa0) + 10

print*, " N=" , N

zeta_home = 1.0D-2
zeta_end = 1.0D0

!zeta_end = 1.0D-2
!zeta_home = 1.0D0

CALL dep_on_distance (kapa0, gama0, kapa, gama, zeta_home, zeta_end)

zeta = zeta_end

!_______________________________________________________________________

!------------------ ZNACHENIE FUNKCII V NAIDENNOM KORNE ---------------------------------


		IF (kapa>0 .and. gama>0 .and. gama < 1 ) THEN
			
			ALLOCATE (A (N, N))	

			IF (ALLOCATED(A)) THEN

				CALL get_matrix (kapa, gama, A, N) 

				det_1=det(A, N) 	

				print*,''
				print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
				print*,''
				print*,'kapa=', kapa
				print*,'gama=',gama
				print*,''
				print*,'det(I + A)=', cdabs(det_1)
				print*,'' 

				CALL eigenfields (kapa, gama, N, A) ! POSTROENIE SOBSTV. POLA V  TOCHKE (kappa, gamma)

				DEALLOCATE (A)

			END IF

		END IF

!-----------------------------------------------------------------------------------------

CALL DATE_AND_TIME(TIME=T2)
print*,''
CALL DELTA_T(T1,T2)
print*,''

END PROGRAM main

!=============================================================================================

!-------------- VIVOD MATRICI NA PECHAT ------------------------------------------------------

!				DO i=0,N
!					DO j=0, N
!						 WRITE(*,'(1X,i4, 1X,i4, 3X,E15.8, 3X,E15.8)'), i,j, a(i,j)
!					END DO
!				END DO

!==============================================================================================
