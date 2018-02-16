Module init_data_class

	
	REAL(8), PARAMETER		:: alfa_i = 1.50D0  !2.630D0		!refractive index of a microcavity

	REAL(8), PARAMETER		:: alfa_e = 1.0D0			!refractive index of outer space

	REAL(8), PARAMETER		:: PI = 3.14159265358979323846D0

	REAL(8), PARAMETER		:: e_gamma = 0.57721566490153D0		! constanta Eilera, sm. Abramovitz, Stigan, gl.1

	COMPLEX(8), PARAMETER	:: im = (0.0D0, 1.0D0)



	REAL(8)					:: beta  = PI/10.0D0					!ugol naklona vistupa spirali (step tilt angle)


	REAL(8)					:: delta = 0.50D0				!relative step width (/a) of indentation of spiral microcavity
															!a is a radius of microcavity


	INTEGER(4)					:: N						! 2N is a number of knots in a quadrature rules (ot 0 do 2N-1)


END MODULE init_data_class

