

        Subroutine jcalc2(current,bx,by,bz)
#include "param.inc"
#include "meter.inc"
c
        dimension current(3,ni,njp2,nk)
        dimension bx(ni,nj,nk),by(ni,nj,nk),bz(ni,nj,nk)
c
c
      SAVE
      common /scrach/
     $     ji(nip1,nj,nk),jj(ni,njp1,nk),jk(ni,nj,nkp1),
     $     jx(ni,nj,nk),jy(ni,nj,nk),jz(ni,nj,nk),
     $     xaver(nip1,njp1,nkp1),yaver(nip1,njp1,nkp1),
     $     zaver(nip1,njp1,nkp1),
     $     binti(ni,njp1,nkp1),bintj(nip1,nj,nkp1),bintk(nip1,njp1,nk),
     $     volinv(njp1),fdum(ni,nj,nk)
*       dimension dumx(nj,nk,2),dumy(nj,nk,2),dumz(nj,nk,2)
      real nose,tail,ji,jj,jk,jx,jy,jz
      data cnfac / 7.9577471e3 /  ! 1/mu0 (1/cm -> 1/m)
Cc
C       calculate the line integrals of b along cell edges
C
C
C
        do  i=1,ni
        do  k=1,nkp1
           kl = mod(k-1,nk) + 1
           kh = mod(k+nk-2,nk) + 1
        do j=2,nj
         binti(i,j,k) = 0.25*(
     $          (bx(i,j,kl)+bx(i,j+1,kl)+bx(i,j,kh)+bx(i,j+1,kh))*
     $             (x(i+1,j,k)-x(i,j,k)) +
     $          (by(i,j,kl)+by(i,j+1,kl)+by(i,j,kh)+by(i,j+1,kh))*
     $             (y(i+1,j,k)-y(i,j,k)) +
     $          (bz(i,j,kl)+bz(i,j+1,kl)+bz(i,j,kh)+bz(i,j+1,kh))*
     $             (z(i+1,j,k)-z(i,j,k)))
         enddo
         enddo
         enddo
C
        do j=1,njp1,nj
           ja = min(j,nj)
        do i=1,ni
           binti(i,j,1) = 0.
           do k=1,nk
              binti(i,j,1) = binti(i,j,1) + bx(i,ja,k)
           enddo
           binti(i,j,1) = binti(i,j,1)*(x(i+1,j,1)-x(i,j,k))/float(nk)
           do k =2, nkp1
              binti(i,j,k) = binti(i,j,1)
           enddo
        enddo
        enddo
C
        do j=1,nj
        do  i=1,ni-1
        do  k=1,nkp1
           kl = mod(nk+k-2,nk) + 1
           kh = mod(nk+k-1,nk) + 1 
         bintj(i+1,j,k) = 0.25*(
     $          (bx(i,j,kl)+bx(i+1,j,kl)+bx(i,j,kh)+bx(i+1,j,kh))*
     $             (x(i,j+1,k)-x(i,j,k)) +
     $          (by(i,j,kl)+by(i+1,j,kl)+by(i,j,kh)+by(i+1,j,kh))*
     $             (y(i,j+1,k)-y(i,j,k)) +
     $          (bz(i,j,kl)+bz(i+1,j,kl)+bz(i,j,kh)+bz(i+1,j,kh))*
     $             (z(i,j+1,k)-z(i,j,k)))
         enddo
         enddo
         enddo
C
        do j=2,nj
        do  i=1,ni-1
        do  k=1,nk
           kl = mod(nk+k-2,nk) + 1
           kh = mod(nk+k-1,nk) + 1 
         bintk(i-1,j,k) = 0.25*(
     $          (bx(i,j,k)+bx(i+1,j,k)+bx(i,j+1,k)+bx(i+1,j+1,k))*
     $             (x(i,j,k+1)-x(i,j,k)) +
     $          (by(i,j,k)+by(i+1,j,k)+by(i,j+1,k)+by(i+1,j+1,k))*
     $             (y(i,j,k+1)-y(i,j,k)) +
     $          (bz(i,j,k)+bz(i+1,j,k)+bz(i,j+1,k)+bz(i+1,j+1,k))*
     $             (z(i,j,k+1)-z(i,j,k)))
         enddo
         enddo
         enddo
C
         do j=1,njp1,nj
         do k=1,nk
         do i=1,nip1
            bintk(i,j,k) = 0.0
         enddo
         enddo
         enddo
c
C
C    find area integrals of j
C
        do 400 i=2,ni-1
        do 400 k=1,nkp1
        do 400 j=1,nj
        jk(i,j,k) = binti(i,j,k)+bintj(i+1,j,k)-
     $          binti(i,j+1,k)-bintj(i,j,k)
  400   continue
C
        do 420 i=2,ni
        do 420 k=1,nk
        do 420 j=1,nj
        ji(i,j,k) = bintj(i,j,k)+bintk(i,j+1,k)-
     $          bintj(i,j,k+1)-bintk(i,j,k)
  420   continue
C
        do 440 i=2,ni-1
        do 440 k=1,nk
        do 440 j=1,nj
        jj(i,j,k) = bintk(i,j,k)+binti(i,j,k+1)-
     $          bintk(i+1,j,k)-binti(i,j,k)
  440   continue
C
C    find the value of j over cell volume
C
        do 500 i=2,ni-1
        do 500 k=1,nk
        do 500 j=1,nj
  500   fdum(i,j,k) = ji(i,j,k)+ji(i+1,j,k)
C
        DO 510 I=1,nip1
        DO 510 K=1,NK
        DO 510 J=1,NJ
        XAVER(I,J,K) = X(I,J,K)+X(I,J+1,K)+X(I,J,K+1)+X(I,J+1,K+1)
        YAVER(I,J,K) = Y(I,J,K)+Y(I,J+1,K)+Y(I,J,K+1)+Y(I,J+1,K+1)
  510   ZAVER(I,J,K) = Z(I,J,K)+Z(I,J+1,K)+Z(I,J,K+1)+Z(I,J+1,K+1)
C
        DO 520 I=2,ni-1
        DO 520 K=1,NK
        DO 520 J=1,NJ
        JX(I,J,K) = FDUM(I,J,K)*(XAVER(I+1,J,K)-XAVER(I,J,K))
        JY(I,J,K) = FDUM(I,J,K)*(YAVER(I+1,J,K)-YAVER(I,J,K))
  520   JZ(I,J,K) = FDUM(I,J,K)*(ZAVER(I+1,J,K)-ZAVER(I,J,K))
C
C
        DO 600 I=2,ni-1
        DO 600 K=1,NK
        DO 600 J=1,NJ
  600   FDUM(I,J,K)= jj(i,j,k)+jj(i,j+1,k)
C
        DO 610 I=2,ni-1
        DO 610 K=1,NK
        DO 610 J=1,NJ+1
        XAVER(I,J,K) = X(I,J,K)+X(I+1,J,K)+X(I,J,K+1)+X(I+1,J,K+1)
        YAVER(I,J,K) = Y(I,J,K)+Y(I+1,J,K)+Y(I,J,K+1)+Y(I+1,J,K+1)
  610   ZAVER(I,J,K) = Z(I,J,K)+Z(I+1,J,K)+Z(I,J,K+1)+Z(I+1,J,K+1)
C
        DO 620 I=2,ni-1
        DO 620 K=1,NK
        DO 620 J=1,NJ
        JX(I,J,K)= JX(I,J,K)+FDUM(I,J,K)*
     $          (XAVER(I,J+1,K)-XAVER(I,J,K))
        JY(I,J,K)= JY(I,J,K)+FDUM(I,J,K)*
     $          (YAVER(I,J+1,K)-YAVER(I,J,K))
  620   JZ(I,J,K)= JZ(I,J,K)+FDUM(I,J,K)*
     $          (ZAVER(I,J+1,K)-ZAVER(I,J,K))
C
C
        DO 700 I=2,ni-1
        DO 700 K=1,NK
        DO 700 J=1,NJ
  700   FDUM(I,J,K)= jk(i,j ,k)+jk(i,j,k+1)
C
        DO 710 I=2,ni-1
        DO 710 K=1,NKP1
        DO 710 J=1,NJ
        XAVER(I,J,K) = X(I,J,K)+X(I+1,J,K)+X(I,J+1,K)+X(I+1,J+1,K)
        YAVER(I,J,K) = Y(I,J,K)+Y(I+1,J,K)+Y(I,J+1,K)+Y(I+1,J+1,K)
  710   ZAVER(I,J,K) = Z(I,J,K)+Z(I+1,J,K)+Z(I,J+1,K)+Z(I+1,J+1,K)
C
        DO 720 I=2,ni-1
        DO 720 K=1,NK
        DO 720 J=1,NJ
        volinv(j) = 0.125*cnfac/volume(i,j,k)
        JX(I,J,K)= (JX(I,J,K)+FDUM(I,J,K)*
     $          (XAVER(I,J,K+1)-XAVER(I,J,K)))*volinv(j)
        JY(I,J,K)= (JY(I,J,K)+FDUM(I,J,K)*
     $          (YAVER(I,J,K+1)-YAVER(I,J,K)))*volinv(j)
  720   JZ(I,J,K)= (JZ(I,J,K)+FDUM(I,J,K)*
     $          (ZAVER(I,J,K+1)-ZAVER(I,J,K)))*volinv(j)
C
        do 800 i=2,nim1
           do 800 j=1,nj
              do 800 k=1,nk
                 current(1,i,j+1,k) = jx(i,j,k)
                 current(2,i,j+1,k) = jy(i,j,k)
                 current(3,i,j+1,k) = jz(i,j,k)
 800    continue
C
        return
        end





