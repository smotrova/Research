
! *************** Function func(h)=the right part of characteristic equation********************

!***************** h=G/k, G=longitudinal propagation number, k=wavenumber of free space**********

!************************************kd=k*d, d=slab thikness************************************

!***********************************************************************************************

Function Func_Even (h)		! Left part of the characteristic equation for the even modes (dla chetnih mod)

Implicit none

REAL(8) alfa, func_even, h, tt,t, kd 

COMMON /arg1/ kd
COMMON /ALFA1/ alfa

!________________________________________________________________________________________________

	tt=dsqrt(alfa*alfa - h*h)
	t=dsqrt(h*h-1)

	func_even=dtan(kd*tt)- alfa*alfa*t/tt	!H pol in the slab plane	(E pol of disk plane)	

	return

END Function func_even

!##############################################################################################

!##############################################################################################

Function Func_Odd(h)	 ! Left part of the characteristic equation for the odd modes (dla nechetnih mod)
Implicit none

REAL(8) alfa, func_odd, h, tt,t, kd 

COMMON /arg1/ kd
COMMON /ALFA1/ alfa

!________________________________________________________________________________________________

	tt=dsqrt(alfa*alfa - h*h)
	t=dsqrt(h*h-1)

	func_odd=dcotan(kd*tt) + alfa*alfa*t/tt		!H pol in the slab plane  (E pol of disk plane)	

	return

END Function func_odd