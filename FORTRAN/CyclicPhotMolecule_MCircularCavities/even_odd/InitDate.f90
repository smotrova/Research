Module init_data_class
	
	REAL(8), PARAMETER		:: alfa = 1.50D0 !1.9670D0	!refractive index of a resonator

	REAL(8), PARAMETER		:: PI = 3.141592653589D0

	COMPLEX(8), PARAMETER	:: im=(0.0D0, 1.0D0)

	REAL(8)					:: zeta		!w/a relative separation between resonators

	INTEGER(4)				:: N		!Truncation oder

	INTEGER(4)				:: N_res = 12 !Number of resonators in PM



END MODULE init_data_class

