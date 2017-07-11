c meter.inc
      COMMON /METER/
     $     X(NIP1,NJP1,NKP1),Y(NIP1,NJP1,NKP1),Z(NIP1,NJP1,NKP1),
     $     VOLUME(NI,NJ,NK),FACEqI(NIP1,NJ,NK),FACEqJ(NI,NJP1,NK),
     $     X2(NIP1,NJP1),Y2(NIP1,NJP1),PHI(NKP1),
     $     COSPHI(NKP1),SINPHI(NKP1),CPHIqM(NK),SPHIqM(NK),DARC(NK),
     $     FACEqK(NI,NJ,NKP1),COSARC(NK),
     $     EDGEqK(NIP1,NJP1),EDGEqI(NIP1,NJ),EDGEqJ(NI,NJP1),
     $     YBARI(NIP1,NJP1),YBARJ(NIP1,NJP1),DX2I(NIP1,NJ),
     $     DY2I(NIP1,NJ),DX2J(NI,NJP1),DY2J(NI,NJP1),
     $     NRINGS,NAVER(NJ),VOLFACT(NJ,NK),volq2d(ni,nj)
      parameter (n2dread = 2*nip1*njp1 + nkp1)
      dimension read2d(n2dread)
      equivalence (read2d,x2)
