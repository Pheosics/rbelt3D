c var2.inc
	COMMON /VAR2/ RHO2(NI,NJ,NK),
     $		VX2(NI,NJ,NK),VY2(NI,NJ,NK),VZ2(NI,NJ,NK),C2(NI,NJ,NK),
     $       	BDVX2(NI,NJ,NK),BDVY2(NI,NJ,NK),BDVZ2(NI,NJ,NK),
     $		BI2(NIP1,NJ,NK),BJ2(NI,NJP1,NK),BK2(NI,NJ,NKP1)
	DIMENSION VAR2(NVAr0)
	EQUIVALENCE (VAR2,RHO2)
