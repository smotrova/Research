!==============================================================================================

!============ Zavisimost chastot i porogov ot rasstoyania mejdu resonatorami===================

!============ Dannie vivodatsa v file =========================================================

!==============================================================================================

SUBROUTINE dep_on_distance (kapa0, gama0, kapa, gama, zeta_home, zeta_end)

USE init_data_class, ONLY : zeta, N

IMPLICIT NONE

REAL(8), INTENT(IN)			:: kapa0, gama0, zeta_end, zeta_home

REAL(8), INTENT(OUT)		:: kapa, gama

REAL(8)						:: kapa00, gama00

REAL(8), PARAMETER			:: step = 1.0D-2

!REAL(8)						:: step 

INTEGER(4)					:: choice

!===============================================================================================

OPEN (unit=2, file='zeta_g_Hpol.dat')
OPEN (unit=3, file='zeta_k_Hpol.dat')
OPEN (unit=13, file='Cut_Hpol.dat')

!--------------------------------------------------------------------------------

zeta = zeta_home

IF ( zeta_home - zeta_end >= 0.0 ) choice = 1
IF ( zeta_home - zeta_end < 0.0 )  choice = 2

!---------------------------------------------------------------------------------

SELECT CASE (choice)

CASE(1)

	kapa00 = kapa0
	gama00 = gama0

	DO WHILE (zeta-zeta_end > 0.0)

!		IF (zeta > 3.34 .and. zeta < 3.39 )  step = 1.0D-3
!		IF (zeta < 3.34 .or. zeta >= 3.39 ) step = 1.0D-2


		print*,'zeta=', zeta

		CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!		CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

		!-------------------------------------------------------------------------------------

		IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
		
			kapa00 = kapa
			gama00 = gama
		
			WRITE(2,'(3X,E15.8, 3X,E15.8)'), zeta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), zeta, kapa
			WRITE(13,'(3X,E15.8, 3X,E15.8)'), zeta, gama0
		
		END IF

		zeta = zeta - step

	END DO

	zeta = zeta_end

	print*,'zeta=', zeta

	CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!	CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

	IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
	
		WRITE(2,'(3X,E15.8, 3X,E15.8)'), zeta, gama
		WRITE(3,'(3X,E15.8, 3X,E15.8)'), zeta, kapa
		WRITE(13,'(3X,E15.8, 3X,E15.8)'), zeta, gama0
		
	END IF
!-------------------------------------------------------------------------------------------------

CASE (2)

	kapa00 = kapa0
	gama00 = gama0
	

	DO WHILE ( zeta - zeta_end < 0.0)

!		IF (zeta > 1.75 .and. zeta < 1.9 )  step = 1.0D-3
!		IF (zeta < 1.75 .or. zeta >= 1.9 ) step = 1.0D-2

		print*,'zeta=', zeta

		CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!		CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

		!-----------------------------------------------------------------------

		IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
		
			kapa00 = kapa
			gama00 = gama
		
			WRITE(2,'(3X,E15.8, 3X,E15.8)'), zeta, gama
			WRITE(3,'(3X,E15.8, 3X,E15.8)'), zeta, kapa
			WRITE(13,'(3X,E15.8, 3X,E15.8)'), zeta, gama0
		
		END IF

		zeta = zeta + step

	END DO

	zeta = zeta_end

	print*,'zeta=', zeta

	CALL GetRoot_Powell (kapa00, gama00, kapa, gama)

!	CALL GetRoot_Newton (kapa00, gama00, kapa, gama)

	IF (kapa>0 .and. gama>0 .and. gama < 1)	THEN
		
		WRITE(2,'(3X,E15.8, 3X,E15.8)'), zeta, gama
		WRITE(3,'(3X,E15.8, 3X,E15.8)'), zeta, kapa
		WRITE(13,'(3X,E15.8, 3X,E15.8)'), zeta, gama0
		
	END IF

END SELECT

CLOSE(2)
CLOSE(3)
CLOSE(13)

zeta = zeta_home

END SUBROUTINE dep_on_distance 

!================================================================================================

		
