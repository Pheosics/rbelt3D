c
      SUBROUTINE BZZ(BX,BY,BZ,BI,BJ,BK)
C
#include "param.inc"
#include "meter.inc"
      DIMENSION BX(NI,NJ,NK),BY(NI,NJ,NK),BZ(NI,NJ,NK)
      DIMENSION BI(NIP1,NJ,NK),BJ(NI,NJP1,NK),BK(NI,NJ,NKP1)
C
      COMMON /SCRACH/
     $     FDUM(NI,NJ,NK),XDUM(NIP1,NJP1,NKP1),
     $     YDUM(NIP1,NJP1,NKP1),ZDUM(NIP1,NJP1,NKP1),
     $     dsum(ni,nk),bxsum(ni,nk),bysum(ni,nk),bzsum(ni,nk)
C
      DIMENSION XAVERI(NIP1,NJ,NK),XAVERJ(NI,NJP1,NK),XAVERK(NI,NJ,NKP1)
      DIMENSION YAVERI(NIP1,NJ,NK),YAVERJ(NI,NJP1,NK),YAVERK(NI,NJ,NKP1)
      DIMENSION ZAVERI(NIP1,NJ,NK),ZAVERJ(NI,NJP1,NK),ZAVERK(NI,NJ,NKP1)
      EQUIVALENCE (XAVERI,XDUM),(XAVERJ,XDUM),(XAVERK,XDUM)
      EQUIVALENCE (YAVERI,YDUM),(YAVERJ,YDUM),(YAVERK,YDUM)
      EQUIVALENCE (ZAVERI,ZDUM),(ZAVERJ,ZDUM),(ZAVERK,ZDUM)
C
C
      DO 200 K=1,NK
      DO 200 J=1,NJ
CDIR$ SHORTLOOP
      DO 200 I=1,NI
 200     FDUM(I,J,K) = BI(I,J,K)*FACEqI(I,J,k)
     $        +BI(I+1,J,K)*FACEqI(I+1,J,K)
C
      DO 210 K=1,NK
      DO 210 J=1,NJ
CDIR$ SHORTLOOP
      DO 210 I=1,NIP1
      XAVERI(I,J,K) = X(I,J,K)+X(I,J+1,K)+X(I,J,K+1)+X(I,J+1,K+1)
      YAVERI(I,J,K) = Y(I,J,K)+Y(I,J+1,K)+Y(I,J,K+1)+Y(I,J+1,K+1)
  210 ZAVERI(I,J,K) = Z(I,J,K)+Z(I,J+1,K)+Z(I,J,K+1)+Z(I,J+1,K+1)
C
      DO 220 K=1,NK
      DO 220 J=1,NJ
CDIR$ SHORTLOOP
      DO 220 I=1,NI
      BX(I,J,K) = FDUM(I,J,K)*(XAVERI(I+1,J,K)-XAVERI(I,J,K))
      BY(I,J,K) = FDUM(I,J,K)*(YAVERI(I+1,J,K)-YAVERI(I,J,K))
  220 BZ(I,J,K) = FDUM(I,J,K)*(ZAVERI(I+1,J,K)-ZAVERI(I,J,K))
C
C
        DO 300 K=1,NK
        DO 300 J=1,NJ
CDIR$ SHORTLOOP
        DO 300 I=1,NI
  300   FDUM(I,J,K)=BJ(I,J,K)*FACEqJ(I,J,K)+BJ(I,J+1,K)*FACEqJ(I,J+1,K)
C
        DO 310 K=1,NK
        DO 310 J=1,NJP1
CDIR$ SHORTLOOP
        DO 310 I=1,NI
        XAVERJ(I,J,K) = X(I,J,K)+X(I+1,J,K)+X(I,J,K+1)+X(I+1,J,K+1)
        YAVERJ(I,J,K) = Y(I,J,K)+Y(I+1,J,K)+Y(I,J,K+1)+Y(I+1,J,K+1)
  310   ZAVERJ(I,J,K) = Z(I,J,K)+Z(I+1,J,K)+Z(I,J,K+1)+Z(I+1,J,K+1)
C
        DO 320 K=1,NK
        DO 320 J=1,NJ
CDIR$ SHORTLOOP
        DO 320 I=1,NI
        BX(I,J,K)= BX(I,J,K)+FDUM(I,J,K)*(XAVERJ(I,J+1,K)-XAVERJ(I,J,K))
        BY(I,J,K)= BY(I,J,K)+FDUM(I,J,K)*(YAVERJ(I,J+1,K)-YAVERJ(I,J,K))
  320   BZ(I,J,K)= BZ(I,J,K)+FDUM(I,J,K)*(ZAVERJ(I,J+1,K)-ZAVERJ(I,J,K))
C
C
        DO 400 K=1,NK
        DO 400 J=1,NJ
CDIR$ SHORTLOOP
        DO 400 I=1,NI
  400   FDUM(I,J,K)=BK(I,J,K)*FACEqK(I,J,K)+BK(I,J,K+1)*FACEqK(I,J,K+1)
C
        DO 410 K=1,NKP1
        DO 410 J=1,NJ
CDIR$ SHORTLOOP
        DO 410 I=1,NI
        XAVERK(I,J,K) = X(I,J,K)+X(I+1,J,K)+X(I,J+1,K)+X(I+1,J+1,K)
        YAVERK(I,J,K) = Y(I,J,K)+Y(I+1,J,K)+Y(I,J+1,K)+Y(I+1,J+1,K)
  410   ZAVERK(I,J,K) = Z(I,J,K)+Z(I+1,J,K)+Z(I,J+1,K)+Z(I+1,J+1,K)
C
        DO 420 K=1,NK
        DO 420 J=1,NJ
CDIR$ SHORTLOOP
        DO 420 I=1,NI
        BX(I,J,K)= BX(I,J,K)+FDUM(I,J,K)*(XAVERK(I,J,K+1)-XAVERK(I,J,K))
        BY(I,J,K)= BY(I,J,K)+FDUM(I,J,K)*(YAVERK(I,J,K+1)-YAVERK(I,J,K))
  420   BZ(I,J,K)= BZ(I,J,K)+FDUM(I,J,K)*(ZAVERK(I,J,K+1)-ZAVERK(I,J,K))
C
C
        DO 500 K=1,NK
        DO 500 J=1,NJ
CDIR$ SHORTLOOP
        DO 500 I=1,NI
        FDUM(I,J,K) = 0.125/VOLUME(I,J,K)
        BX(I,J,K) = BX(I,J,K)*FDUM(I,J,K)
        BY(I,J,K) = BY(I,J,K)*FDUM(I,J,K)
  500   BZ(I,J,K) = BZ(I,J,K)*FDUM(I,J,K)
C
C
C    AVERAGE ALONG THE X-AXIS
C
        DO 900 JJJ=1,2
        DO 900 JJ=1,NRINGS
                IF ( JJJ .EQ. 2 ) then
                  J = NJP1 - JJ
                  jb = njp1-jj
                else
                  J=JJ
                  jb = jj + 1
                endif
C
        DO 800 K=1,NK
CDIR$ SHORTLOOP
        DO 800 I=1,NI
        DSUM(I,K) = VOLUME(I,J,K)
        bxsum(i,k) = volume(i,j,k)*bx(i,j,k)
        bysum(i,k) = volume(i,j,k)*by(i,j,k)
        bzsum(i,k) = volume(i,j,k)*bz(i,j,k)
  800   CONTINUE
C
C
        NSUM = NAVER(JJ)
        NBLOCK = NK/NAVER(JJ)
        DO 785 L=1,NBLOCK
        DO 785 LL=2,NSUM
C
        DO 780 I = 1,NI
        DSUM(I,1+(L-1)*NSUM) = DSUM(I,1+(L-1)*NSUM) +
     $           DSUM(I,LL+(L-1)*NSUM)
        BXSUM(I,1+(L-1)*NSUM) = BXSUM(I,1+(L-1)*NSUM) +
     $           BXSUM(I,LL+(L-1)*NSUM)
        BYSUM(I,1+(L-1)*NSUM) = BYSUM(I,1+(L-1)*NSUM) +
     $           BYSUM(I,LL+(L-1)*NSUM)
        BZSUM(I,1+(L-1)*NSUM) = BZSUM(I,1+(L-1)*NSUM) +
     $           BZSUM(I,LL+(L-1)*NSUM)
  780   CONTINUE
C
  785   CONTINUE
C
        DO 795 L=1,NBLOCK
CDIR$ SHORTLOOP
        DO 790 I = 1,NI
        BXSUM(I,1+(L-1)*NSUM) = BXSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        BYSUM(I,1+(L-1)*NSUM) = BYSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        BZSUM(I,1+(L-1)*NSUM) = BZSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
  790   CONTINUE
        DO 792 LL=2,NSUM
CDIR$ SHORTLOOP
        DO 792 I=1,NI
        BXSUM(I,LL+(L-1)*NSUM) = BXSUM(I,1+(L-1)*NSUM)
        BYSUM(I,LL+(L-1)*NSUM) = BYSUM(I,1+(L-1)*NSUM)
        BZSUM(I,LL+(L-1)*NSUM) = BZSUM(I,1+(L-1)*NSUM)
  792   CONTINUE
C
  795   continue
C
        DO 710 K=1,NK
CDIR$ SHORTLOOP
        DO 710 I=1,NI
        BX(I,J,K) = BXSUM(I,K)
        BY(I,J,K) = BYSUM(I,K)
        BZ(I,J,K) = BZSUM(I,K)
  710   CONTINUE
C
C
  900   CONTINUE
C
C
        RETURN
        END
       subroutine cooler(dt)
#include "param.inc"
#include "var.inc"
#include "help.inc"
#include "meter.inc"
#ifdef STATIC
      common /scrach/ 
#else
      dimension
#endif
     $    coolfac(nj,nk)
      real rtemp,ftemp,rtempsq,rperpsq
C
      do 300 i=1,nicool
        do 200 j=1,nj
        do 200 k=1,nk
           rperpsq = x(i,j,k)**2+y(i,j,k)**2
           rtempsq = rperpsq + z(i,j,k)**2
           ftemp    = rtempsq/rperpsq
           if( ftemp .lt. 70 ) then
              rtemp    = sqrt(rtempsq)/rion
              coolfac(j,k) =  -1.75E-3*cool*dt*
     $             exp( -0.88 * ftemp * rtemp )
           else
              coolfac(j,k) = 0.0
           endif
           c(i,j,k) = c(i,j,k)*(1.+coolfac(j,k))
 200    continue
*         write(9,*) 'i level = ', i
*       call print2d(coolfac,'cool',nj,nk,nj,nk)
  300   continue
C
        return
        end
        SUBROUTINE COURNT(DT,ccell)
#include "param.inc"
#include "var.inc"
#include "var1.inc"
#include "var2.inc"
#include "help.inc"
#include "meter.inc"
#include "bzero.inc"
        logical ccell
C
        DIMENSION BSQ(NI,NJ,NK),DUM(NI,NJ,NK),V(NI,NJ,NK),
     $          VSIG(NI,NJ,NK),dtq(ni,nj,nk),
     $          dtqi(ni,nj,nk),dtqj(ni,nj,nk),dtqk(ni,nj,nk)
        EQUIVALENCE (BSQ,RHO2),(DUM,VX2),(V,VY2),(VSIG,VZ2),(dtq,c2),
     $          (dtqi,bi2),(dtqj,bj2),(dtqk,bk2)
        dimension cycle(ni,4)
C
C
        DO 210 i=1,Nout
        DO 210 K=1,NK
        DO 210 J=1,NJ
        BSQ(i,j,k) = (BX(i,j,k)**2 + BY(i,j,k)**2 + BZ(I,J,K)**2
     $         +0.5*(bsqqk(i,j,k)+bsqqk(i,j,k+1)) )*PI4INV
        DUM(i,j,k) = 1./( RHO(i,j,k) + BSQ(i,j,k)*CAINV**2 )
  210   CONTINUE
        DO 211 i=noutp1,ni
        DO 211 K=1,NK
        DO 211 J=1,NJ
        BSQ(i,j,k) = (BX(i,j,k)**2 + BY(i,j,k)**2 + BZ(I,J,K)**2)*pi4inv
        DUM(i,j,k) = 1./( RHO(i,j,k) + BSQ(i,j,k)*CAINV**2 )
  211   CONTINUE
C
C
        DO 300 K=1,NK
        DO 300 J=1,NJ
        DO 300 I=1,NI
        V(I,J,K) = SQRT(VX(I,J,K)**2+VY(I,J,K)**2+VZ(I,J,K)**2)
        VSIG(i,j,k) = SQRT( C(i,j,k)**2 +
     $           BSQ(i,j,k)*DUM(i,j,k)) + V(i,j,k)
  300 CONTINUE
C
C
        DO 400 K=1,NK
        DO 400 J=1,NJ
        jpp = min(j,njp1-j)
        AVFACT = 1./FLOAT(NAVER(JPP))
        DO 400 I=1,NI
        dtqi(i,j,k) = vsig(i,j,k)
     $   *( amax1(faceqi(i,j,k),faceqi(i+1,j,k))/volume(i,j,k) )
        dtqj(i,j,k) = vsig(i,j,k)
     $   *( amax1(faceqj(i,j,k),faceqj(i,j+1,k))/volume(i,j,k) )
        dtqk(i,j,k) = vsig(i,j,k)
     $   *( amax1(faceqk(i,j,k),faceqk(i,j,k+1))/volume(i,j,k)*avfact)
  400   DTq(i,j,k)=amax1(dtqi(i,j,k),dtqj(i,j,k),dtqk(i,j,k))
C
        dtqmax = 0.
        if ( .not. ccell ) then
          do 410 k=1,nk
          do 410 j=1,nj
          do 410 i=1,ni
           dtqmax = amax1(dtqmax,dtq(i,j,k))
  410     continue
        else
          do 420 k=1,nk
          do 420 j=1,nj
          do 420 i=1,ni
             if ( dtq(i,j,k) .gt. dtqmax ) then
                dtqmax = dtq(i,j,k)
                iqmax = i
                jqmax = j
                kqmax = k
             endif
  420     continue
          write(9,*) 'maximum courant number at cell',iqmax,jqmax,kqmax
          write (9,*) 'v and vsig', v(iqmax,jqmax,kqmax),
     $          vsig(iqmax,jqmax,kqmax)
        endif
C
*       call print3('dtq',dtq,ni,nj,nk)
*       call print3('dtqi',dtqi,ni,nj,nk)
*       call print3('dtqj',dtqj,ni,nj,nk)
*       call print3('dtqk',dtqk,ni,nj,nk)
C
        DT = CRNT/DTqmax
C
        do 490 i=1,ni
          dtqq = 0.
          dtqqi = 0.
          dtqqj = 0.
          dtqqk = 0.
          do 489 k=1,nk
          do 489 j=1,nj
           dtqq  = amax1(dtqq,dtq(i,j,k))
           dtqqi = amax1(dtqqi,dtqi(i,j,k))
           dtqqj = amax1(dtqqj,dtqj(i,j,k))
           dtqqk = amax1(dtqqk,dtqk(i,j,k))
  489     continue
          if ( dtqmax .ge. 8.*dtqq ) then
             cycle(i,1) = 1
          else if (dtqmax .ge. 4.*dtqq ) then
             cycle(i,1) = 2
          else if (dtqmax .ge. 2.*dtqq ) then
             cycle(i,1) = 4
          else
             cycle(i,1) = 8
          endif
          if ( dtqmax .ge. 8.*dtqqi ) then
             cycle(i,2) = 1
          else if (dtqmax .ge. 4.*dtqqi ) then
             cycle(i,2) = 2
          else if (dtqmax .ge. 2.*dtqqi ) then
             cycle(i,2) = 4
          else
             cycle(i,2) = 8
          endif
          if ( dtqmax .ge. 8.*dtqqj ) then
             cycle(i,3) = 1
          else if (dtqmax .ge. 4.*dtqqj ) then
             cycle(i,3) = 2
          else if (dtqmax .ge. 2.*dtqqj ) then
             cycle(i,3) = 4
          else
             cycle(i,3) = 8
          endif
          if ( dtqmax .ge. 8.*dtqqk ) then
             cycle(i,4) = 1
          else if (dtqmax .ge. 4.*dtqqk ) then
             cycle(i,4) = 2
          else if (dtqmax .ge. 2.*dtqqk ) then
             cycle(i,4) = 4
          else
             cycle(i,4) = 8
          endif
  490   continue
C
*        call print2d(cycle,'cycl',ni,4,ni,4)
C
*        DO 500 K=1,NK
*        DO 500 J=1,NJ
*        DO 500 I=1,NI
*        V(I,J,K) = AMIN1(1.,0.9*CA/(V(I,J,K)+1.E-20))
*        VX(I,J,K) = VX(I,J,K)*V(I,J,K)
*        VY(I,J,K) = VY(I,J,K)*V(I,J,K)
*        VZ(I,J,K) = VZ(I,J,K)*V(I,J,K)
*  500   CONTINUE
C
        RETURN
        END
c     
        subroutine ringav
#include "param.inc"
#include "var.inc"
#include "help.inc"
#include "meter.inc"
C
        common /scrach/ DUMMY(NIP1,NK,12)
c
        DIMENSION DSUM(NI,NK),RHOSUM(NI,NK),VXSUM(NI,NK),VYSUM(NI,NK),
     $          VZSUM(NI,NK),CSUM(NI,NK),fisum(nip1,nk),bisum(nip1,nk),
     $       fjsum(ni,nk),bjsum(ni,nk)
        dimension df(nip1,nk),bjs(ni,nk),bjc(ni,nk),cjsum(ni,nk),
     $       sjsum(ni,nk)
        EQUIVALENCE (DUMMY(1,1,1),DSUM),(DUMMY(1,1,2),RHOSUM),
     $          (DUMMY(1,1,3),VXSUM),(DUMMY(1,1,4),VZSUM),
     $          (DUMMY(1,1,5),CSUM),(dummy(1,1,6),fisum),
     $          (dummy(1,1,7),bisum),(dummy(1,1,8),fjsum),
     $       (dummy(1,1,9),bjsum),(dummy(1,1,10),df)
c
C
C       SET UP THE ACCUMULATOR ARRAYS
C
C
        DO 200 K=1,NK
        DO 200 J=1,NJ
        DO 200 I=1,NI
        RHO(I,J,K) = RHO(I,J,K)*VOLUME(I,J,K)
        C(I,J,K) = RHO(I,J,K)*C(I,J,K)**2
C  these accumulators are set to zero to ease implementation of the
C       alfven correction
C     
 200  continue
      do 220 k=1,nk
         do 220 j=1,nj
            do 220 i=1,ni
        VX(I,J,K) = rho(i,j,k)*vx(i,j,k)
        VY(I,J,K) = rho(i,j,k)*vy(i,j,k)
        VZ(I,J,K) = rho(i,j,k)*vz(i,j,k)
 220  CONTINUE     
c     
      call ave2flux
c
C    AVERAGE ALONG THE X-AXIS
C
        DO 900 JJJ=1,2
        DO 900 JJ=1,NRINGS
                IF ( JJJ .EQ. 2 ) then
                  J = NJP1 - JJ
                  jb = njp1-jj
                else
                  J=JJ
                  jb = jj + 1
                endif
C
        DO 800 K=1,NK
        DO 800 I=1,NI
        DSUM(I,K) = VOLUME(I,J,K)
        RHOSUM(I,K) = RHO(I,J,K)
        VXSUM(I,K) = VX(I,J,K)
        VYSUM(I,K) = VY(I,J,K)
        VZSUM(I,K) = VZ(I,J,K)
        CSUM(I,K) = C(I,J,K)
        fjsum(i,k) = faceqj(i,jb,k)
        cjsum(i,k) = faceqj(i,jb,k)*cphiqm(k)
        sjsum(i,k) = faceqj(i,jb,k)*sphiqm(k)
        bjsum(i,k) = bj(i,jb,k)
  800   CONTINUE
C
        do 801 k=1,nk
        do 801 i=1,nip1
        fisum(i,k) = faceqi(i,j,k)
  801   bisum(i,k) = bi(i,j,k)
C
C
        NSUM = NAVER(JJ)
        NBLOCK = NK/NAVER(JJ)
        DO 785 L=1,NBLOCK
        DO 785 LL=2,NSUM
C
        DO 780 I = 1,NI
        DSUM(I,1+(L-1)*NSUM) = DSUM(I,1+(L-1)*NSUM) +
     $           DSUM(I,LL+(L-1)*NSUM)
        RHOSUM(I,1+(L-1)*NSUM) = RHOSUM(I,1+(L-1)*NSUM) +
     $           RHOSUM(I,LL+(L-1)*NSUM)
        VXSUM(I,1+(L-1)*NSUM) = VXSUM(I,1+(L-1)*NSUM) +
     $           VXSUM(I,LL+(L-1)*NSUM)
        VYSUM(I,1+(L-1)*NSUM) = VYSUM(I,1+(L-1)*NSUM) +
     $           VYSUM(I,LL+(L-1)*NSUM)
        VZSUM(I,1+(L-1)*NSUM) = VZSUM(I,1+(L-1)*NSUM) +
     $           VZSUM(I,LL+(L-1)*NSUM)
        CSUM(I,1+(L-1)*NSUM) = CSUM(I,1+(L-1)*NSUM) +
     $           CSUM(I,LL+(L-1)*NSUM)
        FJSUM(I,1+(L-1)*NSUM) = FJSUM(I,1+(L-1)*NSUM) +
     $           FJSUM(I,LL+(L-1)*NSUM)
        BJSUM(I,1+(L-1)*NSUM) = BJSUM(I,1+(L-1)*NSUM) +
     $           BJSUM(I,LL+(L-1)*NSUM)
        cJSUM(I,1+(L-1)*NSUM) = cJSUM(I,1+(L-1)*NSUM) +
     $           cJSUM(I,LL+(L-1)*NSUM)
        sJSUM(I,1+(L-1)*NSUM) = sJSUM(I,1+(L-1)*NSUM) +
     $           sJSUM(I,LL+(L-1)*NSUM)
  780   CONTINUE
        do 781 i=1,nip1
        FISUM(I,1+(L-1)*NSUM) = FISUM(I,1+(L-1)*NSUM) +
     $           FISUM(I,LL+(L-1)*NSUM)
        BISUM(I,1+(L-1)*NSUM) = BISUM(I,1+(L-1)*NSUM) +
     $           BISUM(I,LL+(L-1)*NSUM)
  781   continue
C
  785   CONTINUE
C
        DO 795 L=1,NBLOCK
C
        DO 790 I = 1,NI
        RHOSUM(I,1+(L-1)*NSUM) = RHOSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        VXSUM(I,1+(L-1)*NSUM) = VXSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        VYSUM(I,1+(L-1)*NSUM) = VYSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        VZSUM(I,1+(L-1)*NSUM) = VZSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        CSUM(I,1+(L-1)*NSUM) = CSUM(I,1+(L-1)*NSUM)
     $          /DSUM(I,1+(L-1)*NSUM)
        BJSUM(I,1+(L-1)*NSUM) = BJSUM(I,1+(L-1)*NSUM)
     $          /fjSUM(I,1+(L-1)*NSUM)
        cJSUM(I,1+(L-1)*NSUM) = cJSUM(I,1+(L-1)*NSUM)
     $          /fjSUM(I,1+(L-1)*NSUM)
        sJSUM(I,1+(L-1)*NSUM) = sJSUM(I,1+(L-1)*NSUM)
     $          /fjSUM(I,1+(L-1)*NSUM)
  790   CONTINUE
        do 791 i=1,nip1
        BISUM(I,1+(L-1)*NSUM) = BISUM(I,1+(L-1)*NSUM)
     $          /fiSUM(I,1+(L-1)*NSUM)
  791   continue
        DO 792 LL=2,NSUM
        DO 792 I=1,NI
        RHOSUM(I,LL+(L-1)*NSUM) = RHOSUM(I,1+(L-1)*NSUM)
        VXSUM(I,LL+(L-1)*NSUM) = VXSUM(I,1+(L-1)*NSUM)
        VYSUM(I,LL+(L-1)*NSUM) = VYSUM(I,1+(L-1)*NSUM)
        VZSUM(I,LL+(L-1)*NSUM) = VZSUM(I,1+(L-1)*NSUM)
        CSUM(I,LL+(L-1)*NSUM) = CSUM(I,1+(L-1)*NSUM)
        BjSUM(I,LL+(L-1)*NSUM) = BjSUM(I,1+(L-1)*NSUM)
  792   CONTINUE
        do 793 ll=2,nsum
        do 793 i=1,nip1
        bisum(i,ll+(l-1)*nsum) = bisum(i,1+(l-1)*nsum)
  793   continue
C
 795    continue
c
c        figure coefficients for bj
c
        do 796 l=1,nblock
           lplus = mod(l,nblock) + 1
           lminus = mod(nblock+l-2,nblock) + 1
           do 796 i=1,ni
           delta = bjsum(i,1+(lplus-1)*nsum)-bjsum(i,1+(lminus-1)*nsum)
           cdelt = cjsum(i,1+(lplus-1)*nsum)-cjsum(i,1+(lminus-1)*nsum)
           sdelt = sjsum(i,1+(lplus-1)*nsum)-sjsum(i,1+(lminus-1)*nsum)
           det = 1./(cjsum(i,1+(l-1)*nsum)*sdelt
     $          -sjsum(i,1+(l-1)*nsum)*cdelt)
           bjc(i,l) = (sdelt*bjsum(i,1+(l-1)*nsum)-
     $          delta*sjsum(i,1+(l-1)*nsum))*det
           bjs(i,l) = (-cdelt*bjsum(i,1+(l-1)*nsum)+
     $          delta*cjsum(i,1+(l-1)*nsum))*det
#ifdef DEBUG
           if ( i .eq. 6 ) then 
              write (6,*) '******** l,i,jb',l,i,jb
              write (6,*) lplus,lminus,delta,cdelt,sdelt,det
              write (6,*) bjsum(i,1+(lplus-1)*nsum),
     $             bjsum(i,1+(lminus-1)*nsum),
     $             bjsum(i,1+(l-1)*nsum)
              write (6,*) cjsum(i,1+(lplus-1)*nsum),
     $             cjsum(1,1+(lminus-1)*nsum),
     $             cjsum(1,1+(l-1)*nsum)
              write (6,*) sjsum(i,1+(lplus-1)*nsum),
     $             sjsum(1,1+(lminus-1)*nsum),
     $             sjsum(1,1+(l-1)*nsum)
              write (6,*) bjs(i,l),bjc(i,l)
           endif
#endif
 796       continue
c
C      put the magnetic fields back
C
       do 797 l = 1,nblock
        do 754 i=1,nip1
           df(i,1) = 0.0
           df(i,nsum+1) = 0.0
 754    continue
        do 755 ll=2,nsum
        do 755 i=1,nip1
        df(i,ll) = bisum(i,ll-1+(l-1)*nsum)*faceqi(i,j,ll-1+(l-1)*nsum)
     $      -    bi(i,j,ll-1+(l-1)*nsum) + df(i,ll-1)
 755  continue
      do 756 ll=1,nsum
         do 756 i=1,nip1
        bi(i,j,ll+(l-1)*nsum) = bi(i,j,ll+(l-1)*nsum) + 
     $           ( df(i,ll+1) - df(i,ll) )
 756  continue
      do 757 ll=2,nsum
        do 757 i=1,ni
        bk(i,j,ll+(l-1)*nsum)
     $          = bk(i,j,ll+(l-1)*nsum) + df(i,ll) - df(i+1,ll)
 757  continue
C     
      do 758 i=1,ni
         df(i,1) = 0.0
         df(i,nsum+1) = 0.0
 758  continue
      do 759 ll=2,nsum
      DO 759 i=1,ni
        df(i,ll) = faceqj(i,jb,ll-1+(l-1)*nsum)*
     $      (bjc(i,l)*cphiqm(ll-1+(l-1)*nsum)+
     $                      bjs(i,l)*sphiqm(ll-1+(l-1)*nsum))
     $      -    bj(i,jb,ll-1+(l-1)*nsum) + df(i,ll-1)
 759  continue
      do 760 ll=1,nsum
         do 760 i=1,ni
            bj(i,jb,ll+(l-1)*nsum) = bj(i,jb,ll+(l-1)*nsum) + 
     $            ( df(i,ll+1) - df(i,ll) )
 760     continue
         do 761 ll=2,nsum
            do 761 i=1,ni
        bk(i,jb,ll+(l-1)*nsum) = bk(i,jb,ll+(l-1)*nsum) + df(i,ll)
        bk(i,jb-1,ll+(l-1)*nsum) = bk(i,jb-1,ll+(l-1)*nsum) - df(i,ll)
 761  continue
C
*      write (6,*) (bj(6,2,k),k=1,nk)
C
c
 797  CONTINUE
C
C   PUT THEM BACK
C
        DO 850 K=1,NK
        DO 850 I=1,NI
        VX(I,J,K) = VXSUM(I,K)*VOLUME(I,J,K)
        VY(I,J,K) = VYSUM(I,K)*VOLUME(I,J,K)
        VZ(I,J,K) = VZSUM(I,K)*VOLUME(I,J,K)
        C(I,J,K) = CSUM(I,K)*VOLUME(I,J,K)
        RHO(I,J,K) = RHOSUM(I,K)*VOLUME(I,J,K)
  850   CONTINUE
C
  900   CONTINUE
C
C    CONVERT THE CONSERVED QUANTITIES BACK TO PRIMITIVE VARIABLES
C
C
C
        DO 924 K=1,NK
        DO 924 J=1,NJ
           DO 924 I=1,NI
              vx(i,j,k) = vx(i,j,k)/rho(i,j,k)
              vy(i,j,k) = vy(i,j,k)/rho(i,j,k)
              vz(i,j,k) = vz(i,j,k)/rho(i,j,k)
        c(i,j,k) = sqrt( amax1(1.,c(i,j,k)/rho(i,j,k)) )
        RHO(I,J,K) = RHO(I,J,K)/VOLUME(I,J,K)
  924   CONTINUE
C
        DO 930 K=1,NK
        DO 930 J=1,NJ
        DO 930 I=1,NIP1
  930   BI(I,J,K) = BI(I,J,K)/FACEqI(I,J,K)
C
        DO 931 K=1,NK
        DO 931 J=1,NJP1
        DO 931 I=1,NI
  931   BJ(I,J,K) = BJ(I,J,K)/FACEqJ(I,J,K)
C
        DO 932 K=1,NK
        DO 932 J=1,NJ
        DO 932 I=1,NI
  932   BK(I,J,K) = BK(I,J,K)/FACEqK(I,J,K)
        DO 933 J=1,NJ
        DO 933 I=1,NI
  933   BK(I,J,NKP1) = BK(I,J,1)
C
C     
        return
        end




