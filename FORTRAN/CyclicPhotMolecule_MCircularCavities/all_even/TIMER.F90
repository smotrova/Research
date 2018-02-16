! ===========================
!	CHARACTER*10 T1,T2
!	CALL DATE_AND_TIME(TIME=T1)
!	CALL DATE_AND_TIME(TIME=T2)
!	CALL DELTA_T(T1,T2)

! ===========================

	SUBROUTINE DELTA_T(T1,T2)
	IMPLICIT NONE

	CHARACTER*10,INTENT(IN)::T1,T2
	CHARACTER*15 DT
	INTEGER DELTA,SDELTA,MS,S,M,H

	DELTA=MSEC(T2)-MSEC(T1)		! delta T in MSEC
	SDELTA=INT(DELTA/1000)		! delta T in SEC
	MS=DELTA-SDELTA*1000
	
	H = INT(SDELTA/3600)
	M = INT((SDELTA - 3600*H)/60)
	S = SDELTA - 3600*H - 60*M

	WRITE (DT,'(I2,''h'',I2,''m'',I2,''s'',''.'',I3,''ms'')') H,M,S,MS
	
	PRINT*,"*************************************"
	PRINT*,"* CALCULATION TIME: ",DT, " *"
	PRINT*,"*************************************"

	OPEN(UNIT=1,FILE='TIMER.DAT')
	WRITE (1, '('' *************************************'')')
	WRITE (1,'(I2,''h'',I2,''m'',I2,''s'',''.'',I3,''ms'')') H,M,S,MS
	WRITE (1, '('' *************************************'')')
	CLOSE (1)

	CONTAINS
		INTEGER FUNCTION MSEC(T)
		IMPLICIT NONE
		CHARACTER(10) T
		INTEGER H,M,S,MS
		
		READ (T(1:2),*) H
		READ (T(3:4),*) M
		READ (T(5:6),*) S
		READ (T(8:10),*) MS

		MSEC=(H*3600 + M*60 + S)*1000 + MS
		
		ENDFUNCTION
	END SUBROUTINE