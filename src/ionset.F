	subroutine ionset

c param.inc

#include "param.inc"
      

c var.inc
	COMMON /VAR/lstep,time,RHO(NI,NJ,NK),
     $		VX(NI,NJ,NK),VY(NI,NJ,NK),VZ(NI,NJ,NK),C(NI,NJ,NK),
     $       	BX(NI,NJ,NK),BY(NI,NJ,NK),BZ(NI,NJ,NK),
     $		BI(NIP1,NJ,NK),BJ(NI,NJP1,NK),BK(NI,NJ,NKP1),
     $		psinorth(njp3,nk2p3),psisouth(njp3,nk2p3),
     $		curnorth(njp3,nk2p3),cursouth(njp3,nk2p3),
     $	        psiin(njp1,nk2p1),psiis(njp1,nk2p1)
	DIMENSION VAR(NVAR0),vartot(nvartot),varall(nvarall)
	EQUIVALENCE (VAR,RHO),(vartot,rho),(varall,lstep)

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
     $     NRINGS,NAVER(NJ),VOLFACT(NJ,NK)
      parameter (n2dread = 2*nip1*njp1 + nkp1)
      dimension read2d(n2dread)
      equivalence (read2d,x2)
	common /scrach/xa1(njp1),ya1(njp1,nk2p1),za1(njp1,nk2p1),
     $		xa2(njp1,nk2p3),ya2(njp1,nk2p3),za2(njp1,nk2p3),
     $		xahf(nj,nk2),yahf(nj,nk2),zahf(nj,nk2)
c
c  note the 1.e-2 to change from centimeters to meters for the
c    ionospheric part of the code
c
	do 200 j=1,njp1
	  xa1(j) = x(1,j,1)*1.e-2
c
	  do 190 k=1,nk2P1
	  ya1(j,k) = y(1,j,k)*1.e-2
	  za1(j,k) = z(1,j,k)*1.e-2
	  xa2(j,k+1) = x(2,j,1)*1.e-2
	  ya2(j,k+1) = y(2,j,k)*1.e-2
	  za2(j,k+1) = z(2,j,k)*1.e-2
  190	  continue
c
c
  200   continue
	do 201 k=1,nk2p1
	xa2(1,k+1) = x(2,1,nk2/2+1)*1.e-2
	ya2(1,k+1) = y(2,1,nk2/2+1)*1.e-2
	za2(1,k+1) = z(2,1,nk2/2+1)*1.e-2
	xa2(njp1,k+1) = x(2,njp1,nk2/2+1)*1.e-2
	ya2(njp1,k+1) = y(2,njp1,nk2/2+1)*1.e-2
	za2(njp1,k+1) = z(2,njp1,nk2/2+1)*1.e-2
  201	continue
c
	do 240 j=1,njp1
	xa2(j,1) = x(1,j,1)*1.e-2
	xa2(j,nk2p3) = x(1,j,nk2p1)*1.e-2
	ya2(j,1) = y(1,j,1)*1.e-2
	ya2(j,nk2p3) = y(1,j,nk2p1)*1.e-2
	za2(j,1) = z(1,j,1)*1.e-2
	za2(j,nk2p3) = z(1,j,nk2p1)*1.e-2
  240	continue
c
c
	do 300 k=1,nk2
	do 300 j=1,nj
	xahf(j,k) = 1.25e-3*( x(1,j,k) + x(2,j,k)
     $		+ x(1,j+1,k) + x(2,j+1,k)
     $		+ x(1,j+1,k+1) + x(2,j+1,k+1)
     $		+ x(1,j,k+1) + x(2,j,k+1))
	yahf(j,k) = 1.25e-3*( y(1,j,k) + y(2,j,k)
     $		+ y(1,j+1,k) + y(2,j+1,k)
     $		+ y(1,j+1,k+1) + y(2,j+1,k+1)
     $		+ y(1,j,k+1) + y(2,j,k+1))
	zahf(j,k) = 1.25e-3*( z(1,j,k) + z(2,j,k)
     $		+ z(1,j+1,k) + z(2,j+1,k)
     $		+ z(1,j+1,k+1) + z(2,j+1,k+1)
     $		+ z(1,j,k+1) + z(2,j,k+1))
  300	continue
c
c
	write(6,*) (xa2(j,nk4),j=1,njp1)
        call cord(xa1,ya1,za1,xahf,yahf,zahf,xa2,ya2,za2,
     $        0.01*rearth)
c
	return
	end
