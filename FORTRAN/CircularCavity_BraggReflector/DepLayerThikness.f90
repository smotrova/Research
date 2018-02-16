SUBROUTINE dependance (kapa0, gama0, kapa, gama)

USE init_data_class, ONLY : zeta_h

IMPLICIT NONE

REAL(8), INTENT(INOUT)						:: kapa0, gama0
REAL(8), INTENT(OUT)						:: kapa, gama

REAL(8)										:: kapa00, gama00, zeta_hh

REAL(8)										:: step = 1.0D-3

COMPLEX(8)									:: u



COMPLEX(8), EXTERNAL						:: W_a, W_ext, W_p, W_sum

!___________________________________________________________________________________________

kapa00 = kapa0
gama00 = gama0


OPEN (unit=28, file='w_H_gama.dat')
OPEN (unit=29, file='w_H_kapa.dat')

!OPEN (unit=35, file='wh_Wa.dat')
!OPEN (unit=36, file='wh_Wext.dat')
!OPEN (unit=37, file='wh_Wp.dat')


DO zeta_hh = 0.8D0 , 0.8270D0,  step
!DO zeta_hh = 1.0D0 , 0.010D0,  -step
	
	zeta_h = zeta_hh
	print*,'zeta_h = ', zeta_h 

	CALL GetRoot_Powell(kapa0, gama0, kapa, gama)
!	CALL GetRoot_Newton (kapa0, gama0, kapa, gama)


!	u = W_a(kapa, gama)/W_sum(kapa, gama)
!	WRITE(35,'(1X,E15.8, 3X,E15.8)'), zeta_h, cdabs(u)

!	u = W_ext(kapa, gama)/W_sum(kapa, gama)
!	WRITE(36,'(1X,E15.8, 3X,E15.8)'), zeta_h, cdabs(u)

!	u = W_p(kapa, gama)/W_sum(kapa, gama)
!	WRITE(37,'(1X,E15.8, 3X,E15.8)'), zeta_h, cdabs(u)

	print*, 'kapa=', kapa, 'gama=', gama
	print*,''

	print*,'======================='

!	IF (kapa>0 .and. gama>0 .and. gama < 1 ) THEN

		kapa0 = kapa
		gama0 = gama

		WRITE(28,'(1X,E15.8, 3X,E15.8, 1X,E15.8, 3X,E15.8)'), zeta_h, gama 
		WRITE(29,'(1X,E15.8, 3X,E15.8, 1X,E15.8, 3X,E15.8)'), zeta_h, kapa 

!	ELSE 

!		kapa0 = kapa00
!		gama0 = gama00

!	END IF


END DO

CLOSE(28)
CLOSE(29)


!CLOSE(35)
!CLOSE(36)
!CLOSE(37)

END SUBROUTINE dependance


!=========================================================================================
!=========================================================================================
!=========================================================================================


SUBROUTINE dependance_on_delta (kapa0, gama0, kapa, gama)

USE init_data_class , ONLY : alfa, PI, im, n, ro , initialization, alfa_l, alfa_h, delta

IMPLICIT NONE

REAL(8), INTENT(INOUT)						:: kapa0, gama0
REAL(8), INTENT(OUT)						:: kapa, gama

INTEGER(4)									:: i, i_end, j

REAL(8)										:: a, delta_end, delta_home

REAL(8)										:: step = 1.0D-2

COMPLEX(8)									:: u
REAL(8)										:: u1

COMPLEX(8), EXTERNAL						:: W_a, W_ext, W_p, W_sum

!___________________________________________________________________________________________

delta_end = 0.30D0
delta_home = 0.010D0

!delta_end = 0.010D0
!delta_home = 4.0D0


!-----------------------------------------

!For circular gain area
!delta_end = 1.0D0
!delta_home = 1.0D-2

!------------------------------------------


i_end = int(dabs(delta_home-delta_end)/step)


OPEN (unit=30, file='delta_gama.dat')

OPEN (unit=31, file='delta_kapa.dat')

OPEN (unit=38, file='kapa_gama.dat')


OPEN (unit=32, file='delta_Wa.dat')
OPEN (unit=33, file='delta_Wext.dat')
OPEN (unit=34, file='delta_Wp.dat')


delta = delta_home

DO i = 0, i_end

!	delta = delta - step
	delta = delta + step

	print*,'delta=', delta


	CALL GetRoot_Powell(kapa0, gama0, kapa, gama)
!	CALL GetRoot_Newton (kapa0, gama0, kapa, gama)


	IF (kapa>0 .and. gama>0 .and. gama < 1 ) THEN

			 WRITE(30,'(1X,E15.8, 3X,E15.8)'), delta, gama
 			 WRITE(31,'(1X,E15.8, 3X,E15.8)'), delta, kapa

  			 WRITE(38,'(1X,E15.8, 3X,E15.8)'), kapa, gama

!Overlap coefficients dep on delta

		
			u = W_a(kapa, gama)/W_sum(kapa, gama)
			WRITE(32,'(1X,E15.8, 3X,E15.8)'), delta, dreal(u) !cdabs(u)

 			u = W_ext(kapa, gama)/W_sum(kapa, gama)
			WRITE(33,'(1X,E15.8, 3X,E15.8)'), delta, dreal(u) !cdabs(u)

			u = W_p(kapa, gama)/W_sum(kapa, gama)
  			WRITE(34,'(1X,E15.8, 3X,E15.8)'), delta, dreal(u) !cdabs(u)


!-------------------------------------------------------------------

!For circular gain area

!			WRITE(30,'(1X,E15.8, 3X,E15.8)'), ro(1), gama
! 			WRITE(31,'(1X,E15.8, 3X,E15.8)'), ro(1), kapa
			
!			u = W_a(kapa, gama)/W_sum(kapa, gama)
!			WRITE(32,'(1X,E15.8, 3X,E15.8)'), ro(1), cdabs(u)

! 			u = W_ext(kapa, gama)/W_sum(kapa, gama)
!			WRITE(33,'(1X,E15.8, 3X,E15.8)'), ro(1), cdabs(u)

!--------------------------------------------------------------------

			kapa0 = kapa
			gama0 = gama

	END IF

END DO

CLOSE(30)
CLOSE(31)
CLOSE(38)


CLOSE(32)
CLOSE(33)
CLOSE(34)

!-----------------------------------------------------------------------------------------

END SUBROUTINE dependance_on_delta



