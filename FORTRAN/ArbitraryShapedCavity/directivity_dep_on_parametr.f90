!==============================================================================================

!============ Zavisimost napravlennosti parametra delta===================

!============ Dannie vivodatsa v file =========================================================

!==============================================================================================

SUBROUTINE directivity_dep_on_parameter (kapa0, gama0, kapa, gama, delta_home, delta_end)

USE init_data_class, ONLY : delta, N

IMPLICIT NONE

REAL(8), INTENT(IN)							:: kapa0, gama0, delta_end, delta_home

REAL(8), INTENT(OUT)						:: kapa, gama

COMPLEX(8), ALLOCATABLE, DIMENSION (:, :)	:: A

COMPLEX(8), DIMENSION (0:4*N-1)				:: z

REAL(8)										:: kapa00, gama00, d, teta_max

REAL(8)			 							:: step = 1.0D-3

INTEGER(4)									:: choice

REAL(8), EXTERNAL							:: directivity

COMMON /teta_napravlenie_gl_lepestka/ teta_max


!===============================================================================================

OPEN (unit=2, file='delta_g_Hpol.dat')
OPEN (unit=3, file='delta_k_Hpol.dat')

OPEN (unit=5, file='delta_directivity_Hpol.dat')

OPEN (unit=25, file='delta_fi-max_Hpol.dat')


!--------------------------------------------------------------------------------

delta = delta_home

IF ( delta_home - delta_end >= 0.0 ) choice = 1
IF ( delta_home - delta_end < 0.0 )  choice = 2

!---------------------------------------------------------------------------------

SELECT CASE (choice)

CASE(1)

	kapa00 = kapa0
	gama00 = gama0

	DO WHILE (delta - delta_end > 0.0)

		print*,'delta =', delta

		CALL GetRoot_Powell (kapa00, gama00, kapa, gama)


		!-------------------------------------------------------------------------------------

		IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN

			ALLOCATE (A (0:4*N-1, 0:4*N-1))	

			IF ( ALLOCATED(A) ) THEN

				CALL get_matrix (kapa, gama, A, N) 

				CALL eigenvector (kapa, gama, N, A, z) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

				DEALLOCATE (A)

			END IF
		
			kapa00 = kapa
			gama00 = gama

			d = directivity (kapa, gama, N, z)
		
			WRITE(5,'(3X,E15.8, 3X,E15.8)'), delta, d
			WRITE(25,'(3X,E15.8, 3X,E15.8)'), delta, teta_max

			WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa

		
		END IF

		delta = delta - step

	END DO

	delta = delta_end

	print*,'delta =', delta

	CALL GetRoot_Powell (kapa00, gama00, kapa, gama)


	IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
	
			ALLOCATE (A (0:4*N-1, 0:4*N-1))	

			IF ( ALLOCATED(A) ) THEN

				CALL get_matrix (kapa, gama, A, N) 

				CALL eigenvector (kapa, gama, N, A, z) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

				DEALLOCATE (A)

			END IF

			d = directivity (kapa, gama, N, z)
		
			WRITE(5,'(3X,E15.8, 3X,E15.8)'), delta, d
			WRITE(25,'(3X,E15.8, 3X,E15.8)'), delta, teta_max

			WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa

		
	END IF
!-------------------------------------------------------------------------------------------------

CASE (2)

	kapa00 = kapa0
	gama00 = gama0

	DO WHILE ( delta - delta_end < 0.0)

		print*,'delta =', delta

		CALL GetRoot_Powell (kapa00, gama00, kapa, gama)


		!-----------------------------------------------------------------------

		IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN

			ALLOCATE (A (0:4*N-1, 0:4*N-1))	

			IF ( ALLOCATED(A) ) THEN

				CALL get_matrix (kapa, gama, A, N) 

				CALL eigenvector (kapa, gama, N, A, z) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

				DEALLOCATE (A)

			END IF
		
			kapa00 = kapa
			gama00 = gama

			d = directivity (kapa, gama, N, z)
		
			WRITE(5,'(3X,E15.8, 3X,E15.8)'), delta, d
			WRITE(25,'(3X,E15.8, 3X,E15.8)'), delta, teta_max

			WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa


		END IF

		delta = delta + step

	END DO

	delta = delta_end

	print*,'delta =', delta

	CALL GetRoot_Powell (kapa00, gama00, kapa, gama)


	IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN

			ALLOCATE (A (0:4*N-1, 0:4*N-1))	

			IF ( ALLOCATED(A) ) THEN

				CALL get_matrix (kapa, gama, A, N) 

				CALL eigenvector (kapa, gama, N, A, z) ! nahojdenie sobstv. vectora x: Ax=0, A - zadannaya mat (NxN)

				DEALLOCATE (A)

			END IF
		
			d = directivity (kapa, gama, N, z)
		
			WRITE(5,'(3X,E15.8, 3X,E15.8)'), delta, d
			WRITE(25,'(3X,E15.8, 3X,E15.8)'), delta, teta_max

			WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa



	END IF

END SELECT

CLOSE(2)
CLOSE(3)
CLOSE(5)
CLOSE(25)


END SUBROUTINE directivity_dep_on_parameter

!================================================================================================

		
