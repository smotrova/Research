PROGRAM MAIN

USE init_data_class

IMPLICIT NONE

REAL(8)										:: kapa0, gama0
REAL(8)										:: kapa, gama, zeta_end, zeta_home
COMPLEX(8)									:: det_1
COMPLEX(8), ALLOCATABLE, DIMENSION (:, :)	:: A
INTEGER(4)									:: i, j
COMPLEX(8), EXTERNAL						:: det
COMPLEX(8), EXTERNAL						:: f1, f2
CHARACTER(10)								:: T1,T2
REAL(8)										:: xi_end, xi_home

!-----------------------------------------------------------------------------------------------
!----------- Eff. alfa -------------


!kapa0 = 0.8834D0			!(m=0, alfa=2.63, single disk)
!gama0 = 0.3595D0  

!kapa0 = 0.8950D0			!(m=0, alfa=2.63) (EE)
!gama0 = 0.3595D0

!kapa0 = 2.25D0				!(m=0, alfa=2.63) (EE) ???
!gama0 = 9.86D-2   

!kapa0 = 1.7128D0			!(m=0, alfa=1.402, single disk)
!gama0 = 0.3910D0  

!kapa0 = 1.7128D0			!(m=0, alfa=1.402, EE)
!gama0 = 0.3910D0  


!kapa0 = 0.68D0				!(m=0, alfa=3.374, single disk)
!gama0 = 0.31970D0  

kapa0 = 0.71D0				!(m=0, alfa=3.374, EE)
gama0 = 0.31970D0  

!---------------------------------------------------------

!kapa0 = 2.26D0		!(m=3, alfa=2.63) (single disk)
!gama0 = 9.086D-2   

!kapa0 = 2.2D0		!(m=3, alfa=2.63, two disks, EE) 
!gama0 = 9.086D-2  

!kapa0 = 2.26D0		!(m=3, alfa=2.63) (single disk)
!gama0 = 9.086D-2   


!________________________________________________

!kapa0 = 2.73D0		!(m=4, alfa=2.63) (single disk)
!gama0 = 3.01D-2   

!kapa0 = 2.67D0		!(m=4, alfa=2.63) (EE)
!gama0 = 3.01D-2

!________________________________________________

!kapa0 = 3.1965D0	!(m=5, alfa=2.63) (single disk)
!gama0 = 9.270D-3 

!kapa0 = 3.16D0		!(m=5, alfa=2.63) (EE, two disks )
!gama0 = 9.27D-3  

!_________________________________________________

!kapa0 = 3.6539D0	!(m=6, alfa=2.63) (single disk)
!gama0 = 2.785D-3  

!kapa0 = 3.68D0		!(m=6, alfa=2.63) (OqO)
!gama0 = 2.785D-3  

!kapa0 = 3.650D0		!(m=6, alfa=3.374) (OqE)
!gama0 = 2.785D-3  

!__________________________________________________

!kapa0 = 4.101380	!(m=7, alfa=2.63) (single disk)
!gama0 = 8.2988D-4   

!kapa0 = 4.070		!(m=7, alfa=2.63) (EE)
!gama0 = 8.2988D-4   

!____________________________________________________

!kapa0 = 4.55		!(m=8, alfa=2.63) (OqO)
!gama0 = 2.45D-4 

!_____________________________________________________

!kapa0 = 4.98D0		!(m=9, alfa=2.63) (single disk)
!gama0 = 7.23D-5   

!kapa0 = 4.93D0		!(m=9, alfa=2.63) (EE)
!gama0 = 7.23D-5   

!_____________________________________________________

!kapa0 = 5.40D0		!(m=10, alfa=2.63) (single disk)
!gama0 = 2.12D-5   

!kapa0 = 5.385D0		! (m = 10, alfa=2.63, two disks, EE )
!gama0 = 2.12D-5  

!_____________________________________________________

!kapa0 = 5.45D0		!(m=11, alfa=2.63, two disks, EE)
!gama0 = 2.12D-5  

!_____________________________________________________

!kapa0 = 25.45D0		!(m=33, alfa=1.5)
!gama0 = 1.9D-5  


!___________________________________________________________________________________________

print*,'***********************************'
print*,' '
print*,'Hz polarization;   maximally symmetric class'
print*,' '
print*,'alfa =', alfa
print*,' '
print*,'Number of resonators =', N_res
print*,' '
print*,'***********************************'
print*,' '

!____________________________________________________________________________________________

CALL DATE_AND_TIME(TIME=T1)

N = 50

zeta_home = 1.0D-2
zeta_end = 1.0D0

!zeta_end = 2.0D0
!zeta_home = 2.10D0

CALL dep_on_distance (kapa0, gama0, kapa, gama, zeta_home, zeta_end)

zeta = zeta_end

!_______________________________________________________________________

!------------------ ZNACHENIE FUNKCII V NAIDENNOM KORNE ---------------------------------


		IF (kapa>0 .and. gama>0 .and. gama < 1 ) THEN
			
			ALLOCATE (A (0:N, 0:N))	

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
