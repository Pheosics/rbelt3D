c var.inc

	COMMON /VAR/lstep,time,RHO(NI,NJ,NK),
     $		VX(NI,NJ,NK),VY(NI,NJ,NK),VZ(NI,NJ,NK),C(NI,NJ,NK),
     $       	BX(NI,NJ,NK),BY(NI,NJ,NK),BZ(NI,NJ,NK),
     $		BI(NIP1,NJ,NK),BJ(NI,NJP1,NK),BK(NI,NJ,NKP1),
     $		psinorth(njp3,nk2p3),psisouth(njp3,nk2p3),
     $		curnorth(njp3,nk2p3),cursouth(njp3,nk2p3),
     $	        psiin(njp1,nk2p1),psiis(njp1,nk2p1),
     $		Etoti(ni,njp1,nkp1),etotj(nip1,nj,nkp1),
     $          Etotk(nip1,njp1,nk),
     $	        sig_p_N(njp3,nk2p3),sig_h_N(njp3,nk2p3),
     $	        sig_p_S(njp3,nk2p3),sig_h_S(njp3,nk2p3)
	DIMENSION VAR(NVAR0),vartot(nvartot),varall(nvarall)
	EQUIVALENCE (VAR,RHO),(vartot,rho),(varall,lstep)
	DIMENSION VAROUT(NVAROUT), CURPOTOUT(NCURPOTOUT)
	DIMENSION PSIOUT(NPSIOUT), CONDOUT(NCONDOUT) 
	EQUIVALENCE (VAROUT,RHO), (CURPOTOUT,PSINORTH)
	EQUIVALENCE (PSIOUT,PSIIN), (CONDOUT,SIG_P_N)

