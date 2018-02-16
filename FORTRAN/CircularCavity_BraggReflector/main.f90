PROGRAM MAIN

USE init_data_class, ONLY : M, delta

IMPLICIT NONE

INTEGER(4)									:: i, j

REAL(8)										:: kapa0, gama0, kapa, gama

COMPLEX(8), ALLOCATABLE, DIMENSION (:, :)	:: A

COMPLEX(8)									:: det_1

CHARACTER(10)								:: T1,T2


COMPLEX(8), EXTERNAL						:: det

!-----------------------------------------------------------------------------------------------
!----------- Eff. alfa -------------


!kapa0 = 0.8834D0		!(m=0, alfa=2.63, single disk)
!gama0 = 0.3595D0  

!---------------------------------------------------------

!kapa0 = 2.5 			!(moda (1,2) , alfa=2.63)
!gama0 = 7.740D-2

!kapa0 = 1.5			!(moda (1,1) , alfa=2.63)
!gama0 = 7.740D-2

!---------------------------------------------------------

!kapa0 = 2.26D0			!(m=3, alfa=2.63) (single disk)
!gama0 = 9.086D-2   

!________________________________________________

!kapa0 = 2.780D0		!(m=4, alfa=2.63) (single disk)
!gama0 = 3.01D-3  

!kapa0 = 2.185D0		!(m=4, alfa=3.374) (single disk)
!gama0 = 4.487D-3

!----------------------------------------------------------

!kapa0 = 4.10D0			!(m=7, alfa=2.63) (single disk)
!gama0 = 1.79D-3   


!------------------------------------------------------------

!kapa0 = 3.190D0		!(m=5, alfa=2.63) (single disk)!
!gama0 = 9.27D-3   

!kapa0 = 1.510D0		!(m=5, alfa=2.63) (single disk)!
!gama0 = 9.27D-3   

!kapa0 = 5.40D0			!(m=10, alfa=2.63) (single disk)!
!gama0 = 2.1D-3   



!----------------------------------------------------------
! non-uniform circular gain

!kapa0 = 6.70D0			!(m=7, alfa=2.63) (single disk)
!gama0 = 1.79D-3   

!------------------------------------------------------------


!___________________________________________________________________________________________

CALL DATE_AND_TIME(TIME=T1)


!CALL Initialization (alfa, ro, M)

!print*,'***********************************'
!print*,' '
!print*,'Hz polarisation'
!print*,' '
!print*,'alfa =' !, alfa
!print*,' '

!		DO i = 1, M+1

!			 WRITE(*,'(1X,i4, 3X,E15.8, 3X,E15.8)'), i, alfa(i)

!		END DO
!print*,''
!print*,'Number of layers =', M
!print*,' '
!print*,'***********************************'
!print*,' '

!____________________________________________________________________________________________

!kapa0 = 0.880D0			!(m=0, alfa=2.63, single disk)
!gama0 = 0.03595D0  

!kapa0 = 1					!(moda (1,1) , alfa=2.63)
!gama0 = 7.740D-2

!kapa0 = 2.5 				!(moda (1,2) , alfa=2.63)
!gama0 = 7.740D-2


kapa0 = 3.71					!(moda (7,1) , alfa=2.63)
gama0 = 1.70D-6


!--------------------------------------------------

CALL dependance_on_delta (kapa0, gama0, kapa, gama)

!CALL dependance (kapa0, gama0, kapa, gama)

!_______________________________________________________________________

!------------------ ZNACHENIE FUNKCII V NAIDENNOM KORNE ---------------------------------
!	    OPEN (unit=40, file='k_g.dat')


		CALL GetRoot_Powell(kapa0, gama0, kapa, gama)

!		WRITE(40,'(1X,E15.8, 3X,E15.8)'), kapa, gama

!		CLOSE (40)

	
!		CALL GetRoot_Newton (kapa0, gama0, kapa, gama)

!		CALL OT (kapa, gama)

!		IF (kapa>0 .and. gama>0 .and. gama < 1 ) THEN

		IF (kapa>0 .and. gama>0 ) THEN

			
			ALLOCATE (A(1:2*M, 1:2*M))	

			IF (ALLOCATED(A)) THEN

				CALL get_matrix (kapa, gama, A, M) 

!======== Vivod matrici A na pechat
								
!				DO i=1,2*M
!					DO j=1, 2*M
!						 WRITE(*,'(1X,i4, 1X,i4, 3X,E15.8, 3X,E15.8)'), i,j, a(i,j)
!					END DO
!				END DO
!====================================================================================


				det_1=det(A, 2*M) 	

				print*,''
				print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
				print*,''
				print*,'kapa=', kapa
				print*,'gama=',gama
				print*,''
				print*,'det(A)=', cdabs(det_1)
				print*,'' 
				CALL eigenfields (kapa, gama, M, A)

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
