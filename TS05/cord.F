	subroutine cord(xa1,ya1,za1,xahf,yahf,zahf,xa2,ya2,za2,r0)

c param.inc

#include "param.inc"
      

c iongeo.inc
      COMMON/COORD/TETA1(NJP1,NK2P1),PHI1(NJP1,NK2P1),
     $		x1ion(njp1,nk2p1),y1ion(njp1,nk2p1),
     1             TETA2(NJP1,NK2P3),PHI2(NJP1,NK2P3),
     2		   x2ion(njp1,nk2p3),y2ion(njp1,nk2p3),
     1             TETAhf(NJ,NK2),PHIhf(NJ,NK2),
     2		   xhfion(nj,nk2),yhfion(nj,nk2)
      COMMON/COEF/Apeder(njp1,nk2p3,kk,kk),Ahall(njp1,nk2p3,kk,kk)
     $			,aj(njp1,nk2p3,kk)
	real xa1(njp1),ya1(njp1,nk2p1),za1(njp1,nk2p1),
     1		xa2(njp1,nk2p3),ya2(njp1,nk2p3),za2(njp1,nk2p3),
     2		xahf(nj,nk2),yahf(nj,nk2),zahf(nj,nk2)
	real s2thet(njp1,nk2p3),routr(njp1,nk2p3)
	real shfthet(nj,nk2),routhf(nj,nk2)
c
	rinr=sqrt(xa1(1)**2+ya1(1,1)**2+za1(1,1)**2)
	do 2 j=1,nk2p3
	do 2 i=1,njp1
	routr(i,j) = xa2(i,j)**2+ya2(i,j)**2+za2(i,j)**2
	routr(i,j) = sqrt(routr(i,j))
    2	continue
	do 4 j=1,nk2p1
	do 4 i=1,njp1
	s2thet(i,j) = (r0/rinr)*(1.-(za1(i,j)/rinr)**2)
	teta1(i,j) = asin(sqrt(s2thet(i,j)))
	phi1(i,j) = atan2(ya1(i,j),xa1(i))
    4	continue
	do 5 j=1,nk2p1
	do 5 i=1,njp1
	x1ion(i,j) = teta1(i,j)*cos(phi1(i,j))
	y1ion(i,j) = teta1(i,j)*sin(phi1(i,j))
    5	continue
	do 12 j=1,nk2p3
	do 12 i=1,njp1
	s2thet(i,j) = (r0/routr(i,j))*(1.-(za2(i,j)/routr(i,j))**2)
	teta2(i,j) = asin(sqrt(s2thet(i,j)))
	phi2(i,j) = atan2(ya2(i,j),xa2(i,j))
   12	continue
	do 15 j=1,nk2p3
	do 15 i=1,njp1
	  x2ion(i,j) = teta2(i,j)*cos(phi2(i,j))
	  y2ion(i,j) = teta2(i,j)*sin(phi2(i,j))
   15	continue
	do 22 j=1,nk2
	do 22 i=1,nj
	routhf(i,j) = xahf(i,j)**2+yahf(i,j)**2+zahf(i,j)**2
	routhf(i,j) = sqrt(routhf(i,j))
   22	continue
	do 32 j=1,nk2
	do 32 i=1,nj
	s2thet(i,j) = (r0/routhf(i,j))*(1.-(zahf(i,j)/routhf(i,j))**2)
	tetahf(i,j) = asin(sqrt(s2thet(i,j)))
	phihf(i,j) = atan2(yahf(i,j),xahf(i,j))
   32	continue
	do 35 j=1,nk2
	do 35 i=1,nj
	  xhfion(i,j) = tetahf(i,j)*cos(phihf(i,j))
	  yhfion(i,j) = tetahf(i,j)*sin(phihf(i,j))
   35	continue
c
	return
	end

