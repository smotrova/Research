!==============================================================================================

!============ Zavisimost chastot i porogov ot parametra delta===================

!============ Dannie vivodatsa v file =========================================================

!==============================================================================================

SUBROUTINE dep_on_parameter (kapa0, gama0, kapa, gama, delta_home, delta_end)

USE init_data_class, ONLY : delta, N

IMPLICIT NONE

REAL(8), INTENT(IN)			:: kapa0, gama0, delta_end, delta_home

REAL(8), INTENT(OUT)		:: kapa, gama

REAL(8)						:: kapa00, gama00

REAL(8), PARAMETER			:: step = 1.0D-3

!REAL(8)						:: step 

INTEGER(4)					:: choice

!===============================================================================================

OPEN (unit=2, file='delta_g_Hpol.dat')
OPEN (unit=3, file='delta_k_Hpol.dat')

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

!				IF (delta > 0.3 .and. delta < 0.4) THEN 
				
!				step  = 1.D-3

!		ELSE

!				step = 5.D-3

!		END IF


		CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!		CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

		!-------------------------------------------------------------------------------------

		IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
		
			kapa00 = kapa
			gama00 = gama
		
			WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa

		
		END IF

		delta = delta - step

	END DO

	delta = delta_end

	print*,'delta =', delta

	CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!	CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

	IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
	
		WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
		WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa
				
	END IF
!-------------------------------------------------------------------------------------------------

CASE (2)

	kapa00 = kapa0
	gama00 = gama0

	DO WHILE ( delta - delta_end < 0.0)

		print*,'delta =', delta

!		IF (delta > 0.3 .and. delta < 0.4) THEN 
				
!				step  = 1.D-3

!		ELSE

!				step = 5.D-3

!		END IF


		CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!		CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

		!-----------------------------------------------------------------------

		IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
		
			kapa00 = kapa
			gama00 = gama
		
			WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa
			
		
		END IF

		delta = delta + step

	END DO

	delta = delta_end

	print*,'delta =', delta

	CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!	CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

	IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
		
		WRITE(2,'(3X,E15.8, 3X,E15.8)'), delta, gama
		WRITE(3,'(3X,E15.8, 3X,E15.8)'), delta, kapa
		
	END IF

END SELECT

CLOSE(2)
CLOSE(3)


END SUBROUTINE dep_on_parameter 

!================================================================================================

		
