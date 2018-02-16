PROGRAM MAIN

USE init_data_class


IMPLICIT NONE

REAL(8)										:: kapa0, gama0, kapa1, gama1
REAL(8)										:: kapa, gama, t_s, t_p, t

REAL(8)										:: delta_home, delta_end

INTEGER(4)									:: m,p

COMPLEX(8)									:: A1, A2, B1, B2, C1, C2, D1, D2

COMPLEX(8)									:: det_1, nu, u, s1, s2, u1, u2
COMPLEX(8), ALLOCATABLE, DIMENSION (:, :)	:: A
INTEGER(4)									:: i, j

CHARACTER(10)								:: T1,T2

REAL(8)										:: x1, x2, z, w, v1, v2, v,rr

REAL(8), EXTERNAL							:: R_n1, R_n2, n1_n2, R, curv, L, dx_dt, dy_dt, x, y, d2x_dt2, d2y_dt2

COMPLEX(8), EXTERNAL						:: det, det0 

REAL(8), EXTERNAL							:: f1, f2

!-----------------------------------------------------------------------------------------------
!----------- Eff. alfa -------------
! H pol
!-----------------------------------------------------------

!kapa0 = 0.88384636D0		!(m=0, alfa=2.63, single disk)
!gama0 = 0.35953742D0  

!kapa0 = 1.40496D0			!(m=1, alfa=2.63, single disk)
!gama0 = 0.2750857D0  

!kapa0 = 1.830940D0			!(m=2, alfa=2.63, single disk)
!gama0 = 0.206590D0  

!kapa0 = 2.261362D0			!(m=3, alfa=2.63, single disk)
!gama0 = 9.10588D-2   

!kapa0 = 2.7269738D0		!(m=4, alfa=2.63, single disk)
!gama0 = 3.027626D-2   

!kapa0 = 3.197921D0			!(m=5, alfa=2.63, single disk)
!gama0 = 9.31260D-3 

!kapa0 = 3.66D0				!(m=6, alfa=2.63, single disk)
!gama0 = 2.85D-3  

!kapa0 = 4.103245D0			!(m=7, alfa=2.63, single disk)
!gama0 = 8.3484D-4   

!kapa0 = 4.55				!(m=8, alfa=2.63, single disk)
!gama0 = 2.45D-4 

!kapa0 = 4.98D0				!(m=9, alfa=2.63, single disk)
!gama0 = 7.23D-5   

!kapa0 = 5.40D0				!(m=10, alfa=2.63, single disk)
!gama0 = 2.12D-5   



!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================


! H pol

!-----------------------------------------

!kapa0 = 3.85D0			!(H(6,2)_1, alfa=2.63, delta = 1) N = 400, spiral
!gama0 = 0.05   

!kapa0 = 3.9D0			!(H(6,2)_2, alfa=2.63, delta = 1) N = 400, spiral
!gama0 = 0.03   


!-----------------------------------------

!kapa0 = 4.08D0			!(m=9, alfa=2.63, delta = 1) H(9,1)_1 N = 400, spiral
!gama0 = 0.03   

!kapa0 = 4.05D0			!(m=9, alfa=2.63, delta = 1) H(9,1)_2 N = 400, spiral
!gama0 = 0.01   


!-----------------------------------------

!kapa0 = 3.27D0			!(H(7,1)_1 , alfa=2.63, delta = 1) N = 400, spiral
!gama0 = 0.002   

!kapa0 = 3.26D0			!(H(7,1)_2 , alfa=2.63, delta = 1) N = 400, spiral
!gama0 = 0.04   

!-----------------------------------------

!kapa0 = 1.43D0			!(H(2,1)_1 , alfa=2.63, delta = 1) N = 100, spiral
!gama0 = 0.18   

!kapa0 = 1.45D0			!(H(2,1)_2 , alfa=2.63, delta = 1) N = 100, spiral
!gama0 = 0.25   



!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================



!kapa0 = 1.89			!(kite, delta = 0.5, quasi-H(2,1), alfa=2.63, odd, D = 1.978 ) N = 50
!gama0 = 0.197  

!kapa0 = 1.96			!(kite, delta = 0.5, quasi-H(2,1), alfa=2.63, even, D = 3.431) N = 50
!gama0 = 0.19  

!------------------------------------------------------

!kapa0 = 3.34			!(kite, delta = 0.5, quasi - H(5,1), alfa=2.63, odd) N = 50
!gama0 = 0.11  

!kapa0 = 3.37			!(kite, delta = 0.5, quasi - H(5,1), alfa=2.63, even ) N = 50
!gama0 = 0.068  

!--------------------------------------------------------

!kapa0 = 3.51			!(kite, delta = 0.5, quasi - H(3,2), alfa=2.63, ) N = 50
!gama0 = 0.06  

!kapa0 = 3.72			!(kite, delta = 0.5, quasi - H(3,2), alfa=2.63, ) N = 50
!gama0 = 0.13  

!---------------------------------------------------------------------------------

!kapa0 = 3.85			!(kite, delta = 0.5, quasi - H(6,1), alfa=2.63, ) N = 50
!gama0 = 0.07  

!kapa0 = 4.02			!(kite, delta = 0.5, quasi - H(6,1), alfa=2.63, ) N = 50
!gama0 = 0.025  

!----------------------------------------------------------------------------------

!kapa0 = 3.61			!(kite, delta = 0.5, quasi - H(?,?), alfa=2.63, ) N = 50
!gama0 = 0.07  


!kapa0 = 3.93			!(kite, delta = 0.5, quasi - H(?,?), alfa=2.63, ) N = 50
!gama0 = 0.11  
  

!-------------------------------------------------------------------------------

!kapa0 = 3.32			!(kite, delta = 0.5, quasi - H(0,2), alfa=2.63, ) N = 50
!gama0 = 0.139  

!-------------------------------------------------------------------------------

!kapa0 = 3.65			!(kite, delta = 1.0, quasi - H(,), alfa=2.63, ) N = 50
!gama0 = 0.085  

!-------------------------------------------------------------------------------

!kapa0 = 3.77			!(kite, delta = 1.0, quasi - H(,), alfa=2.63, ) N = 50
!gama0 = 0.11  

!-------------------------------------------------------------------------------

!kapa0 = 3.89			!(kite, delta = 1.0, quasi - H(,), alfa=2.63, ) N = 50
!gama0 = 0.099  

!-------------------------------------------------------------------------------

!kapa0 = 3.91			!(kite, delta = 2.0, quasi - H(,), alfa=2.63, ) N = 50
!gama0 = 0.3 



!--------------------------------------------------------
!kapa0 = 4.530D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(2,3) )
!gama0 = 0.05D0   

!kapa0 = 4.620D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(5,2) )
!gama0 = 0.03D0   

!kapa0 = 4.70D0					!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(8,1) )
!gama0 = 0.06D0   

!kapa0 = 4.810D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(6,2) )
!gama0 = 0.07D0   

!-------------------------------

!kapa0 = 5.0D0					!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(3,3) )
!gama0 = 0.07D0   

!kapa0 = 5.110D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(3, 3) )
!gama0 = 0.04D0   

!--------------------------------

!kapa0 = 5.14D0					!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(2,4) )
!gama0 = 0.09D0   


!------------------------------

!kapa0 = 5.210D0					!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(9,1), odd, D = 2.95 )
!gama0 = 0.03D0   

!kapa0 = 5.230D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(9,1), even, D = 4.19 )
!gama0 = 0.05D0   

!------------------------------


!kapa0 = 5.270D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(4,3),  )
!gama0 = 0.05D0   

!kapa0 = 5.60D0					!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(7,2),  )
!gama0 = 0.04D0   


!kapa0 = 5.720D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(7,2,),  )
!gama0 = 0.04D0   


!kapa0 = 5.750D0				!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(7,2),  )
!gama0 = 0.01D0   


!kapa0 = 5.80D0					!(kite, alfa=2.63, delta = 0.5, N = 50, quasi- H(7,2),  )
!gama0 = 0.03D0   

!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================


!kapa0 = 0.8365					!(limacon, alfa=2.63,  N = 30 ,  delta = 0.5, quasi-H(0,1), D = 2.2  )
!gama0 = 0.36   

!-----------------------

!kapa0 = 1.29						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(1,1), D =   ) 
!gama0 = 0.27   

!kapa0 = 1.36						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(1?,1), D =  ) 
!gama0 = 0.28   

!-------------------------

!kapa0 = 1.72						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(2,1), D =  ) 
!gama0 = 0.2 

!-------------------------

!kapa0 = 2.13						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(3,1), D =   ) 
!gama0 = 0.09 

!-------------------------------

!kapa0 = 2.57						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(4,1), D =   ) 
!gama0 = 0.03 

!-------------------------------

!kapa0 = 3.01						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(5,1), D =   ) 
!gama0 = 0.01 

!-------------------------------

!kapa0 = 3.16						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(1,2) or (1,3), D = 2.54  ) even
!gama0 = 0.12

!-------------------------------


!kapa0 = 3.41						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(3,2), D = 2.86   ) odd
!gama0 = 0.13

!-------------------------------


!kapa0 = 3.44						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(6,1), D = 4.88  ) odd
!gama0 = 0.001 

!-------------------------------


!kapa0 = 3.59						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(2, 2)??, D = 1.93   ) 
!gama0 = 0.11

!-------------------------------

!kapa0 = 3.8D0						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(7,1), D = 7.7  )  even
!gama0 = 0.11D0 

!kapa0 = 3.86						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(7,1), D = 5.62  )  odd
!gama0 = 0.004 

!-------------------------------

!kapa0 = 3.82						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(4,2) D = 2.32   ) odd
!gama0 = 0.12

!kapa0 = 3.84						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(4,2) D = 3.51   ) even
!gama0 = 0.12

!-------------------------------


!kapa0 = 4.07						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(2,3), D = 3.47   ) even
!gama0 = 0.095

!kapa0 = 4.14						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(2,3), D = 2.4   ) odd
!gama0 = 0.1

!-------------------------------

!kapa0 = 4.22						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(5,2), D =  2.73  ) even
!gama0 = 0.078 


!kapa0 = 4.23						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(5.2), D = 3.43   ) odd
!gama0 = 0.08 

!-------------------------------


!kapa0 = 4.27						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(8,1), D = 7.02   ) even
!gama0 = 0.007 

!-------------------------------


!kapa0 = 4.35						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(1,3), D = 3.52   ) even
!gama0 = 0.09

!-------------------------------


!kapa0 = 4.58						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(3,3), D =  4.34  ) odd
!gama0 = 0.09 

!kapa0 = 4.61						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(3.3), D = 2.99   ) even
!gama0 = 0.1 

!-------------------------------


!kapa0 = 4.65						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(6,2), D = 4.39  )  odd
!gama0 = 0.03 

!-------------------------------

!kapa0 = 4.6973						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(9,1), D = 4.94  ) odd
!gama0 = 0.0144 



!kapa0 = 4.6994						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(9,1), D = 5.8  ) even
!gama0 = dexp(-1.846171D0*dlog(10.D0))


!kapa0 = 4.6948						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(9,1), D = 4.94  ) odd
!gama0 = dexp(-1.836D0*dlog(10.D0))


!-------------------------------

!kapa0 = 4.74						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(2,4), D = 2.62    ) odd
!gama0 = 0.09 

!-------------------------------

!kapa0 = 4.94						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(1,4), D =  2.91   ) even
!gama0 = 0.09 

!-------------------------------


!kapa0 = 5.065						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(10,1), D = 5.14  ) odd
!gama0 = 0.007 

!--------------------------------

!kapa0 = 5.07						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(4.3), D = 3.19  )	odd 
!gama0 = 0.09

!-------------------------------

!kapa0 = 5.15						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(7,2), D = 4.8  ) odd 
!gama0 = 0.019 

!-------------------------------

!kapa0 = 5.21						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(3.4), D = 2.33    ) even
!gama0 = 0.08 

!-------------------------------

!kapa0 = 5.33						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(2,4), D = 2.15    ) odd
!gama0 = 0.1

!-------------------------------


!kapa0 = 5.47						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(11,1), D = 6.07  ) odd
!gama0 = 2.5D-4 

!-------------------------------

!kapa0 = 5.5						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(5.3), D = 2.73  ) even
!gama0 = 0.1

!kapa0 = 5.52						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(5.3), D = 2.69  ) odd
!gama0 = 0.09 

!-------------------------------
!kapa0 = 5.55						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(8.2), D = 5.63  ) odd
!gama0 = 0.07

!kapa0 = 5.6						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(8,2), D = 6.8   ) even
!gama0 = 0.014 

!-------------------------------

!kapa0 = 5.7						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(4,4), D =   ) odd
!gama0 = 0.07

!-------------------------------

!kapa0 = 5.87						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(12,1), D = 7.08  ) even
!gama0 = 1.D-5 

!-------------------------------

!kapa0 = 5.89						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(), D =   ) 
!gama0 = 0.07

!-------------------------------

!kapa0 = 5.96						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(), D =   ) 
!gama0 = 0.07

!-------------------------------



!kapa0 = 6.01						!(limacon, alfa=2.63,  N = 40 ,  delta = 0.5, quasi-H(9,2), D = 7.7  ) even
!gama0 = 1.D-4 

!-------------------------------


!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================
!=================================================================================================


CALL DATE_AND_TIME(TIME=T1)

N = 50

!kapa0 = 2.55		!(kite, delta = 0.0, quasi-H(2,1), alfa=1.5 ) N = 50
!gama0 = 0.3 

!kapa0 = 23.56		!(kite, delta = 0.0, quasi-H(1,23), alfa=1.5 ) N = 50
!gama0 = 0.034		!										
!---------------------------------------

!kapa0 = 21.46		!(kite, delta = 0.165, alfa=1.5 ) N = 100, WG - even, D = 5.95
!gama0 = 0.03 

!kapa0 = 23.56		!(kite, delta = 0.165, alfa=1.5 ) N = 100, WG - even, D = 6.78
!gama0 = 0.01										

!F-B
!kapa0 = 23.56		!(kite, delta = 0.165, alfa=1.5 ) N = 100, (delta from 0 to 0.165) => F-P, D = 5.73 odd
!gama0 = 0.034										

!kapa0 = 20.87		!(kite, delta = 0.165, alfa=1.5 ) N = 100, WG - even, D = 5.69
!gama0 = 0.014 

!##############

!kapa0 = 8.73			!(kite, delta = 0.165, alfa=1.5 ) N = 50, D =2.8 , FP-even
!gama0 = dexp(-1.04*dlog(10.D0))


!kapa0 = 8.85			!(kite, delta = 0.165, alfa=1.5 ) N = 50, D =3.329 , quasi-WG(10,1)-even
!gama0 = dexp(-1.15*dlog(10.D0))

!kapa0 = 8.85			!(kite, delta = 0.165, alfa=1.5 ) N = 50, D =3.86 , quasi-WG(10,1)-odd
!gama0 = 0.07312

!kapa0 = 8.81			!(kite, delta = 0.165, alfa=1.5 ) N = 50, D =3.54 , FP-odd
!gama0 = dexp(-1.05*dlog(10.D0))

!---------

!kapa0 = 23.03		!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 5.447, NON WG-odd
!gama0 = dexp(-1.45*dlog(10.D0)) 


!kapa0 = 23.15			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 5.065, NON-WG - even  
!gama0 = dexp(-1.46*dlog(10.D0)) 


!kapa0 = 23.19			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 4.75, NON-WG - odd  
!gama0 = dexp(-1.46*dlog(10.D0)) 


!kapa0 = 23.297			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 4.6264, FP - even  (? babochka)
!gama0 = dexp(-1.45*dlog(10.D0)) 

!kapa0 = 23.32			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 5.43, NON-WG odd (babochka)  
!gama0 = dexp(-1.43*dlog(10.D0)) 


!kapa0 = 23.24			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 4.35, NON-WG even  
!gama0 = dexp(-1.425*dlog(10.D0)) 


!kapa0 = 23.37			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 5.4762,  FP-even (vdol osi symmetrii)
!gama0 = dexp(-1.45*dlog(10.D0)) 


!kapa0 = 23.489			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 5.73, FP - odd 
!gama0 = dexp(-1.46*dlog(10.D0)) 

!kapa0 = 23.5		!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 3.72, NON-WG-even
!gama0 = dexp(-1.34*dlog(10.D0)) 


!kapa0 = 23.56		!(kite, delta = 0.165, alfa=1.5 ) N = 100, WG - even, D = 6.78
!gama0 = 0.01										

!kapa0 = 23.58			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 4.68 , NON-WG-odd  
!gama0 = dexp(-1.42*dlog(10.D0)) 


!kapa0 = 23.7			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 4.73, NON-WG - even 
!gama0 = dexp(-1.45*dlog(10.D0)) 


!kapa0 = 23.77			!(kite, delta = 0.165, alfa=1.5 ) N = 100, D = 6.7,  NON-WG-odd
!gama0 = dexp(-1.41*dlog(10.D0)) 

!##############

!--------------------------------------------------

!kapa0 = 33.31		!(kite, delta = 0.165, quasi-H(), alfa=1.5 ) N = 100 delta =0.165 D = 9.76  (N = 150 D = 9.73) WGM
!gama0 = 0.019 

!kapa0 = 33.48		!(kite, delta = 0.165, quasi-H(), alfa=1.5 ) N = 100 delta =0.165 D = 
!gama0 = 0.0019 

!kapa0 = 33.79		!(kite, delta = 0.165, quasi-H(), alfa=1.5 ) N = 50 delta = .165 D = 7.7	WGM
!gama0 = 0.01

!kapa0 = 34.03		!(kite, delta = 0.0, quasi-H(), alfa=1.5 ) N = 100 delta = .165 D = 9.76	WGM (N = 200 D = 9.52)
!gama0 = 0.02

!kapa0 = 39.62		!(kite, delta = 0.165, quasi-H(1,), alfa=1.5 ) N = 150 delta =0.165, WGM D = 10.01
!gama0 = 0.00107 


!----------------------------------------------------------------
!kapa0 = 40.3			!(kite, delta = 0., quasi-H(1,), alfa=1.5 ) N = 100 delta =0.165 D = 10.07 WGM
!gama0 = dexp(-2.55*dlog(10.D0)) 


! esli nachinat s delta = 0.165
!kapa0 = 40.31		!(kite, delta = 0.165, quasi-H(1,), alfa=1.5 ) N = 150 delta =0.165 D = 7.93  NON-WGM
!gama0 = 0.019 

!---------------------------------------


!kapa0 = 70.7		!(kite, delta = 0.165, quasi-H(), alfa=1.5 ) N = 150 delta = 0.165 D = 10	WGM
!gama0 = 0.01


!kapa0 = 70.7		!(kite, delta = 0.5, quasi-H(), alfa=1.5 ) N = 150 D = 7.415				NON-WGM
!gama0 = 0.01

!--------------------------

!====================================================
!delta = 0.5 ===============================================
!====================================================

!kapa0 = 5.8			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 4.36,  quasi- (5,1) odd
!gama0 = dexp(-0.79*dlog(10.D0)) 

!kapa0 = 5.84			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 2.85,  quasi-(5,1) even
!gama0 = dexp(-0.85*dlog(10.D0))

!---------------------------------

!kapa0 = 8.91			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 3.23, NON-WG-even
!gama0 = dexp(-0.99*dlog(10.D0))

!kapa0 = 9.04			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D =3.49 ,  F-P odd
!gama0 = dexp(-1.01*dlog(10.D0)) 

!--------------

!kapa0 = 8.25			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 5.65 , Horse-shoe mode, odd
!gama0 = dexp(-1.04*dlog(10.D0))

!kapa0 = 8.86			!(kite, delta = 0.0, alfa=1.5 ) N = 50, WGH(10,1) circle
!gama0 = 4.75E-2

!kapa0 = 9.06			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 5.65, Horse-shoe mode, odd
!gama0 = dexp(-1.095*dlog(10.D0))

!kapa0 = 9.21			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 4.56, NON-WG-even 
!gama0 = dexp(-1.01*dlog(10.D0))

kapa0 = 9.89			!(kite, delta = 0.5, alfa=1.5 ) N = 50, D = 5.81, Horse-shoe mode, odd
gama0 = dexp(-1.13*dlog(10.D0))

!----------------------------------------

!
!kapa0 = 20.65						!(kite, delta = 0.5, alfa=1.5 ) N = 200,  , D = 
!gama0 = dexp(-1.51*dlog(10.D0))		!	

!
!kapa0 = 21.48						!(kite, delta = 0.5, alfa=1.5 ) N = 200,  , D = 10.85, Horse-shoe mode, odd
!gama0 = dexp(-1.53*dlog(10.D0))		!	

!
!kapa0 = 22.29						!(kite, delta = 0.5, alfa=1.5 ) N = 200,  , D = 10.035, Horse-shoe mode, odd
!gama0 = dexp(-1.53*dlog(10.D0))		!	

!
!kapa0 = 23.1						!(kite, delta = 0.5, alfa=1.5 ) N = 200,  , D = 9.6, Horse-shoe mode, odd
!gama0 = dexp(-1.52*dlog(10.D0))		!	

!
!kapa0 = 23.95						!(kite, delta = 0.5, alfa=1.5 ) N = 200,  , D = 10.15, Horse-shoe mode, odd
!gama0 = dexp(-1.54*dlog(10.D0))		!	


! ???
!kapa0 = 23.69		!(kite, delta = 0.5, alfa=1.5 ) N = 200, , D = 
!gama0 = dexp(-1.52*dlog(10.D0))		!	

!============================================================================


print*,'***********************************'
print*,' '
print*,'Hz polarization'
print*,' '
print*,'alfa =', alfa_i
print*,' '
print*,'delta =', delta
print*,' '
print*,'A number of knots, N =', N
print*,' '
print*,'***********************************'
print*,' '
print*,'kapa0 =', kapa0
print*,'gama0 =', gama0
print*,''


!==========================================================
!==========================================================



!CALL GetRoot_Powell(kapa0, gama0, kapa, gama)

!CALL GetRoot_Newton (kapa0, gama0, kapa, gama)



!===================================================================
!===================================================================

delta_home = 0.50D0
delta_end = 0.0D0

CALL dep_on_parameter (kapa0, gama0, kapa, gama, delta_home, delta_end)

!CALL directivity_dep_on_parameter (kapa0, gama0, kapa, gama, delta_home, delta_end)



!===================================================================
!===================================================================

! contour

!OPEN (unit=2, file='CONTUR.dat')

!	DO p = 0, 200

!		t = PI/100*p
	
!		WRITE(2,'(2X,E13.6, 3X,E13.6 )'), x(t), y(t)

!	END DO

!CLOSE(2)


!===================================================================
!===================================================================
!Truncation error

!N = 25

!OPEN (unit=11, file='tr-error.dat')

!v = 1.0D0

!print*, "N = ", N

!CALL GetRoot_Powell(kapa0, gama0, kapa, gama)

!CALL GetRoot_Newton (kapa0, gama0, kapa, gama)


!DO WHILE ( N < 400 )

!kapa1 = kapa
!gama1 = gama

!N = 2*N

!print*, "N = ", N

!CALL GetRoot_Powell(kapa0, gama0, kapa, gama)

!CALL GetRoot_Newton (kapa0, gama0, kapa, gama)

!v1 = dabs((kapa1 - kapa))
!v2 = dabs((gama1 - gama))

!print*,'N=', N, '  ', 'v1 = ', v1

!print*,'N=', N, '  ', 'v2 = ', v2

!v = max(v1, v2)

!v = dsqrt(v1*v1 + v2*v2)

!print*,'N=', N, '  ', 'err = ', v

!WRITE(11,'(1X,i4, 3X,E13.6 )'), N, v

!END DO

!CLOSE(11)

!===================================================================
!===================================================================

!relief det

!OPEN (unit=9, file='det.dat')

!DO kapa = 8.0D0, 10.0D0, 0.01D0				
!	DO gama = 1.D-3, 1.D-1, 0.01D0

!	DO gama = dlog10(1.0D-3), dlog10(5.0D-1), 0.01D0

!			ALLOCATE (A (0:4*N-1, 0:4*N-1))	

!			IF ( ALLOCATED(A) ) THEN

!				gama1 = dexp(gama*dlog(10.D0))
				
!				CALL get_matrix (kapa, gama1, A, N) 

!				CALL get_matrix (kapa, gama, A, N) 


!				det_1 = det(A, N) 	

!				WRITE(9,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), kapa, gama, cdabs(det_1)
!				WRITE(*,'(2X,E13.6, 3X,E13.6, 3X,E13.6 )'), kapa, gama, cdabs(det_1)

!				DEALLOCATE (A)

!			END IF

!END DO
!END DO

!CLOSE(9)

 

!------------------ ZNACHENIE FUNKCII V NAIDENNOM KORNE ---------------------------------

!Eigenfields

!*******************************************************************


	IF (kapa>0 .and. gama>0 .and. gama < 1 ) THEN
		
			ALLOCATE (A (0:4*N-1, 0:4*N-1))	

			IF ( ALLOCATED(A) ) THEN

				CALL get_matrix (kapa, gama, A, N) 

				det_1 = det(A, N) 	

				print*,'' 
				print*,''
				print*,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
				print*,''
				print*,'kapa=', kapa
				print*,'gama=',gama
				print*,''
				print*,'det(I + A)=', cdabs(det_1)
				print*,'' 

!				CALL eigenfields (kapa, gama, N, A) ! POSTROENIE SOBSTV. POLA V  TOCHKE (kappa, gamma)

				CALL eigenfields_kite (kapa, gama, N, A) ! POSTROENIE SOBSTV. POLA dlya kite-shaped cavity V  TOCHKE (kappa, gamma)


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
