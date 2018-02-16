!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


PROGRAM MAIN
IMPLICIT NONE 


INTEGER      m,i

INTEGER ::	 n=14					!order of Bessel function

REAL(8) ::   PI=3.141592653589D0

REAL(8) ::   alfa =3.3740D0			! Refraction index

REAL(8) ::   l=0.0850D0				! l=d/a, 2d = thikness of disk, a=radius of disk

REAL(8)      k0, kn

REAL(8) ::   eps=0.0000010D0			!Accuracy

COMMON /ALFA/ alfa
COMMON /arg/ l
COMMON /Accuracy/ eps 

!******************************************************************************************
!+++++++++++++++++++++++++++++++!Nachalnoe priblijenie+++++++++++++++++++++++++++++++++++++

!E-pol
!	kapa0=Pi*(2*m+n+0.5)/2/alfa						
!	gama0=dlog((alfa-1)/(alfa+1))/(2*kapa0)	

		

!H-pol
!   kapa0=Pi*(2*m+n-0.5)/2/alfa
!	gama0=dlog(-(1-alfa)/(alfa+1))/(2*kapa0)		

!*********************************************************************************************
	k0=1					!k0<ka< kn
	kn=13

print*,'E polarization       ','n=',n
print*,''
print*, 'Refraction index of bulk material alfa=', alfa
print*,''
print*, 'd/a=', l, '   a=disk radius;  ', 'd=disk thikness'
print*,''
print*,''
print*,''
print*,'1. Eigenpairs (kappa, gamma)'
print*,''
print*,'2. Eigenfields'
print*,''
read*, i
print*,''
print*,''


IF (i==1)then
	CALL kappa_gamma (n, k0, kn)
end if

IF (i==2)then
	CALL eigenfields (n, k0)
end if

print*,''

	
END PROGRAM MAIN


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-----------------------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





