!################################################################################################
!~~~~~~~~~~~~~~~~~~~~~Program for finding alfa_eff for each normilazed frequency~~~~~~~~~~~~~~~~~
!################################################################################################

Function alfa_eff (kapa)

USE  dfimsl
Implicit none

REAL(8) :: PI=3.141592653589D0

REAL(8) ::   alfa =3.3740D0			! Refraction index

REAL(8)	   cut_kapa1, cut_kapa2, x, upper_limit

REAL(8)    h, B, A, alfa_eff

REAL(8)    EPS, ERRABS, ERREL
REAL(8)    kd, l						! kd=k*d, where k=wavwernumber of free space
		    							! 2d=slab thikness
REAL(8)	   kapa

INTEGER    i,m

COMMON /arg/ l
COMMON /arg1/ kd
COMMON /ALFA1/ alfa


!################################################################################################
!_______________________________________________________________________________________________

!	h=G/k, G=longitudinal propagation number of slab, k=wavenumber of free space
!   a < h < b

!   V=k*d*sqrt(alfa*alfa-1)
!   IF   0<V<PI/2 then only fundumental mode of slab exists 

!_______________________________________________________________________________________________

!###############################################################################################

	eps= 1.0E-6		! Accuracy of finding zeros
		
	a=1+eps
	b=alfa-eps

	cut_kapa1 = PI/l/dsqrt(alfa*alfa-1)/2
	cut_kapa2 = PI/l/dsqrt(alfa*alfa-1)

!################################################################################################

IF (a>b) THEN

	print*,''
	print*,'!!!  ERROR from Effective Index !!!'
	print*,''
	print*,'a must be less than b '
	print*,''

ELSE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Find alfa_eff for each kappa, kappa=ka, kappa*l=k*d, l=d/a
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	kd=kapa*l

	a=1+eps
	b=alfa-eps
										!Find eigenvalue of dielektric slab for the first even mode 
	x=a
	upper_limit=b

	CALL even_roots (a, b, x, upper_limit )		

	alfa_eff=x

!	IF (kd > cut_kapa1*l)	THEN		!Find eigenvalue of dielektric slab for the first odd mode 
	
!		a=1+eps
!		b=alfa-eps

!		upper_limit=x
!		x=a

!		CALL odd_roots (a, b, x, upper_limit )	
!		alfa_eff=x

!	END IF

END IF

!################################################################################################

END FUNCTION alfa_eff