      SUBROUTINE METRIC
c     input: nip1,njp1,nkp1,x2(nip1,njp1),y2(nip1,njp1),phi(nkp1)
c     output: the remaining metric coeffients.

#include "param.inc"

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

      dimension
     $     DX(LI,LJP1),DY(LI,LJP1),DELTAX(LIP1,LJ),DELTAY(LIP1,LJ),
     $     XqI(LIP1,LJP1),XqJ(LIP1,LJP1), ! for the XBHI and YBHI 
     $     YqI(LIP1,LJP1),YqJ(LIP1,LJP1),
     $     DELDX(NI,NJ),DELDY(NI,NJ),
     $     J0(NJ,NK),JPSI(NJ,NK),JETA(NJ,NK),
     $     FACEK(NI,NJ),VOLq(NI,NJ)
*     $     VOLINV(NI,NJ),XVEC(NIP1,NJP1),YVEC(NIP1,NJP1),SINV(NIP1,NJP1),
*     $    JXY(NJP1,NKP1),disti(nip1,njp1),distj(nip1,njp1),phion(njp1,nkp1),
*      REAL JXY,
      real J0,JPSI,JETA
C
      DO 150 J=1,NJP1
      DO 150 I=1,NI
      DX(I,J) = X2(I+1,J)-X2(I,J)
  150 DY(I,J) = Y2(I+1,J)-Y2(I,J)
C
C
      call mijset(ni,nj)
c     The "edges" of the cells are to be sections of a polynomial.
      CALL XBHI(DX,XqI,1)
      CALL XBHI(DY,YqI,2)
C
      DO 200 J=1,NJ
      DO 200 I=1,NIP1
      DELTAX(I,J) = X2(I,J+1)-X2(I,J)
  200 DELTAY(I,J) = Y2(I,J+1)-Y2(I,J)
C
      DO 210 J=1,NJ
      DO 210 I=1,NI
      DELDX(I,J) = DX(I,J+1)-DX(I,J)
      DELDY(I,J) = DY(I,J+1)-DY(I,J)
      J0(I,J) = DX(I,J)*DELTAY(I,J) - DY(I,J)*DELTAX(I,J)
      JETA(I,J)= DX(I,J)*DELDY(I,J) - DY(I,J)*DELDX(I,J)
      JPSI(I,J) = DELTAY(I,J)*DELDX(I,J) - DELTAX(I,J)*DELDY(I,J)
      FACEK(I,J) = J0(I,J) + 0.5*(JETA(I,J)+JPSI(I,J))
      VOLq(I,J) = J0(I,J)*(Y2(I,J)+0.5*(DY(I,J)+DELTAY(I,J))+
     $                0.25*DELDY(I,J)) +
     $     JETA(I,J)*(0.5*Y2(I,J)+ 0.333333*DY(I,J) +0.25*DELTAY(I,J)
     $              + 0.16666666*DELDY(I,J)) +
     $     JPSI(I,J)*(0.5*Y2(I,J) + 0.25*DY(I,J) +0.333333*DELTAY(I,J)
     $                  + 0.16666666*DELDY(I,J))
 210  CONTINUE
C     
      CALL YBHI(DELTAX,XqJ,1)
      CALL YBHI(DELTAY,YqJ,2)
C
C     COMPUTE THE UNIT NORMAL, LENGTH OF THE AVERAGE FACEq,
C               AND THE MEAN NORMAL VELOCITY
C
      DO 300 J=1,NJ
      DO 300 I=1,NIP1
      DX2I(I,J) = DELTAX(I,J)
      DY2I(I,J) = DELTAY(I,J)
      YBARI(I,J) = 0.5*(Y2(I,J+1)+Y2(I,J))
 300  EDGEqI(I,J) = SQRT( DX2I(I,J)**2 + DY2I(I,J)**2)
C
C
C       NOW DO THE SAME FOR THE FLUXES IN THE J-DIRECTION
C
C------------------------------------------------------------------
C
C
C
C       COMPUTE THE UNIT NORMAL, LENGTH OF THE AVERAGE FACE$,
C               AND THE MEAN NORMAL VELOCITY
C
      DO 400 J=1,NJP1
      DO 400 I=1,NI
      DX2J(I,J) = DX(I,J)
      DY2J(I,J) = DY(I,J)
      YBARJ(I,J) = 0.5*(Y2(I+1,J)+Y2(I,J))
 400  EDGEqJ(I,J) = SQRT( DX2J(I,J)**2 + DY2J(I,J)**2)
C
C
      DO 500 K=1,NK
        COSPHI(K) = COS(PHI(K))
        SINPHI(K) = SIN(PHI(K))
        CPHIqM(K) = COS(0.5*(PHI(K)+PHI(K+1)))
        SPHIqM(K) = SIN(0.5*(PHI(K)+PHI(K+1)))
        DARC(K) = 2.0*SIN(0.5*(PHI(K+1)-PHI(K)))
        COSARC(K) = COS( 0.5*(PHI(K+1)-PHI(K)))
  500   CONTINUE
        COSPHI(NKP1) = COSPHI(1)
        SINPHI(NKP1) = SINPHI(1)
C
C       NOTE THAT THESE DEFINITIONS MAKE THE COMPUTATIONAL (I,J,K) SYSTEM
C               RIGHT-HANDED
C
        DO 600 K=1,NKP1
        DO 600 J=1,NJP1
        DO 600 I=1,NIP1
        X(I,J,K) = X2(I,J)
        Y(I,J,K) = Y2(I,J)*cosPHI(K)
  600   Z(I,J,K) = Y2(I,J)*sinPHI(K)
C
        DO 700 K=1,NK
        DO 700 J=1,NJ
        DO 700 I=1,NI
  700   VOLUME(I,J,K) = VOLq(I,J)*DARC(K)*COSARC(K)
C
        DO 710 K=1,NKP1
        DO 710 J=1,NJ
        DO 710 I=1,NI
  710   FACEqK(I,J,K) = FACEK(I,J)
C
        DO 720 K=1,NK
        DO 720 J=1,NJP1
        DO 720 I=1,NI
  720   FACEqJ(I,J,K) =
     $    DARC(K)*YBARJ(I,J)*SQRT( (COSARC(K)*DY2J(I,J))**2 +
     $                   DX2J(I,J)**2 )
C
        DO 730 K=1,NK
        DO 730 J=1,NJ
        DO 730 I=1,NIP1
  730   FACEqI(I,J,K) = DARC(K)*YBARI(I,J)*
     $    SQRT((COSARC(K)*DY2I(I,J))**2 + DX2I(I,J)**2 )
C
*        call print3('faceqi',faceqi,nip1,nj,nk)
*        call print3('faceqj',faceqj,ni,njp1,nk)
*        call print3('faceqk',faceqk,ni,nj,nkp1)
*        call print3('volume',volume,ni,nj,nk)
C
        do 840 k=1,nk
           write(9,*)
     $         ' k,cphiqm(k),sphiqm(k),darc(k),cosarc(k) ='
           write(9,*)
     $           k,cphiqm(k),sphiqm(k),darc(k),cosarc(k)
 840    continue
        do 940 j=1,nj
        if ( j .eq. 3 ) then
           write(9,*) ' j,ybari(1,j),dx2i(1,j),dy2i(1,j) = '
           write(9,*)
     $          j,ybari(1,j),dx2i(1,j),dy2i(1,j)
           write(9,*)
        endif
 940    continue
c
        RETURN
        END
C
c
c
        subroutine mijset(i,j)

c mdims.inc
      common /mdims/ isweep,level3,
     $    mi,mip1,mip2,mip3,mip4,mip5,mim1,mim2,mim3,mim4,
     $    mj,mjp1,mjp2,mjp3,mjp4,mjp5,mjm1,mjm2,mjm3,mjm4,
     $    mj2,mj2p1,mj2p2,mi2,mi2p1
        mi =i
        mi2 = mi/2
        mi2p1 = mi2+1
        mip1 = mi+1
        mip2 = mip1+1
        mip3 = mip2+1
        mip4 = mip3+1
        mip5 = mip4+1
        mim1 = mi-1
        mim2 = mim1-1
        mim3 = mim2-1
        mim4 = mim3-1
C
        mj =j
        mj2 = mj/2
        mj2p1 = mj2+1
        mj2p2 = mj2+2
        mjp1 = mj+1
        mjp2 = mjp1+1
        mjp3 = mjp2+1
        mjp4 = mjp3+1
        mjp5 = mjp4+1
        mjm1 = mj-1
        mjm2 = mjm1-1
        mjm3 = mjm2-1
        mjm4 = mjm3-1
C
        return
        end

      SUBROUTINE XBHI(F,FHF,LVAR)

c param.inc

#include "param.inc"
      

c mdims.inc
      common /mdims/ isweep,level3,
     $    mi,mip1,mip2,mip3,mip4,mip5,mim1,mim2,mim3,mim4,
     $    mj,mjp1,mjp2,mjp3,mjp4,mjp5,mjm1,mjm2,mjm3,mjm4,
     $    mj2,mj2p1,mj2p2,mi2,mi2p1

c help.inc
      INTEGER ORDER,ORD2
      COMMON /HELP/
     $     GAMMA,GAMM1,GMM1IN,GMqGM1,CA,CAINV,ALCON,PI,PI4,PI4INV,
     $     GMINV,GM1,GP1,GM1INV,WIDTH,GGM1,GGM1IN,
     $     PI8INV,TWELVI,RUSVIS,EVMAX,BETMAX,ALPHA,BETA,
     $     PDMB,ORDER,ORD2,BDIF,NXSTEP,NYSTEP,
     $     FCOEF(20),DCOEF(20),DROP(20),
     $     NSTART,NSTOP,NDUMP,NVAXDUMP,NDIAGN,CRNT,NPICT,
     $     nplotit,nionplt,nprint,
     $     nicool,cool,
     $     niondump
C
      DIMENSION F(LI,LJP1),FHF(LIP1,LJP1)
      DIMENSION FDUM(LIqq,LJP1)
C
      DO 200 J=1,MJP1
      DO 200 I=1,MI

      FDUM(I+ord2,J) = F(I,J)
  200 continue

C
C
      call xhi2(fdum,lvar)
C
C
      DO 299 J=1,MJP1
      DO 299 I=1,MIP1

      FHF(I,J) = 0.
 299  continue

      DO 300 K=1,ORD2
      DO 300 J=1,MJP1
      DO 300 I=1,MIP1

      FHF(I,J) = FHF(I,J) + FCOEF(K)*
     $    ( FDUM(I+K-1,J) + FDUM(I+ORDER-K,J) )
 300    continue

      RETURN
      END
        SUBROUTINE xHI2(FDUM,LVAR)

c param.inc

#include "param.inc"
      

c mdims.inc
      common /mdims/ isweep,level3,
     $    mi,mip1,mip2,mip3,mip4,mip5,mim1,mim2,mim3,mim4,
     $    mj,mjp1,mjp2,mjp3,mjp4,mjp5,mjm1,mjm2,mjm3,mjm4,
     $    mj2,mj2p1,mj2p2,mi2,mi2p1

c help.inc
      INTEGER ORDER,ORD2
      COMMON /HELP/
     $     GAMMA,GAMM1,GMM1IN,GMqGM1,CA,CAINV,ALCON,PI,PI4,PI4INV,
     $     GMINV,GM1,GP1,GM1INV,WIDTH,GGM1,GGM1IN,
     $     PI8INV,TWELVI,RUSVIS,EVMAX,BETMAX,ALPHA,BETA,
     $     PDMB,ORDER,ORD2,BDIF,NXSTEP,NYSTEP,
     $     FCOEF(20),DCOEF(20),DROP(20),
     $     NSTART,NSTOP,NDUMP,NVAXDUMP,NDIAGN,CRNT,NPICT,
     $     nplotit,nionplt,nprint,
     $     nicool,cool,
     $     niondump
        DIMENSION FDUM(LIqq,LJP1)
C
        if ( isweep .eq. 1 ) then
C
        DO 200 J=1,MJP1
        DO 200 I=1,ORD2
        FDUM(I,J) = FDUM(ORD2+1,J)

      FDUM(MI+ORD2+I,J) = FDUM(MI+ORD2,J)
 200  continue

C
        else
C
        do 300 j=1,mjp1
        do 300 i=1,ord2
        fdum(i,j) = fdum(mi+i,j)

      fdum(mi+ord2+i,j) = fdum(ord2+i,j)
 300  continue

C
        endif
C
        RETURN
        END
        SUBROUTINE YBHI(F,FHF,LVAR)

c param.inc

#include "param.inc"
      

c mdims.inc
      common /mdims/ isweep,level3,
     $    mi,mip1,mip2,mip3,mip4,mip5,mim1,mim2,mim3,mim4,
     $    mj,mjp1,mjp2,mjp3,mjp4,mjp5,mjm1,mjm2,mjm3,mjm4,
     $    mj2,mj2p1,mj2p2,mi2,mi2p1

c help.inc
      INTEGER ORDER,ORD2
      COMMON /HELP/
     $     GAMMA,GAMM1,GMM1IN,GMqGM1,CA,CAINV,ALCON,PI,PI4,PI4INV,
     $     GMINV,GM1,GP1,GM1INV,WIDTH,GGM1,GGM1IN,
     $     PI8INV,TWELVI,RUSVIS,EVMAX,BETMAX,ALPHA,BETA,
     $     PDMB,ORDER,ORD2,BDIF,NXSTEP,NYSTEP,
     $     FCOEF(20),DCOEF(20),DROP(20),
     $     NSTART,NSTOP,NDUMP,NVAXDUMP,NDIAGN,CRNT,NPICT,
     $     nplotit,nionplt,nprint,
     $     nicool,cool,
     $     niondump
C
        DIMENSION F(LIP1,LJ),FHF(LIP1,LJP1)
        DIMENSION FDUM(LIP1,LJqq)
C
        DO 200 J=1,MJ
        DO 200 I=1,MIP1
  200   FDUM(I,J+ORD2) = F(I,J)
C
C
        call yhi2(fdum,lvar)
C
        DO 299 J=1,MJP1
        DO 299 I=1,MIP1
  299   FHF(I,J) = 0.
        do 300 k=1,ORD2
        DO 300 J=1,MJP1
        DO 300 I=1,MIP1
  300   FHF(I,J) = FHF(I,J) + FCOEF(K)*
     $          ( FDUM(I,J+K-1) + FDUM(I,J+ORDER-K) )
C
        RETURN
        END
        SUBROUTINE YHI2(FDUM,LVAR)

c param.inc

#include "param.inc"
      

c mdims.inc
      common /mdims/ isweep,level3,
     $    mi,mip1,mip2,mip3,mip4,mip5,mim1,mim2,mim3,mim4,
     $    mj,mjp1,mjp2,mjp3,mjp4,mjp5,mjm1,mjm2,mjm3,mjm4,
     $    mj2,mj2p1,mj2p2,mi2,mi2p1

c help.inc
      INTEGER ORDER,ORD2
      COMMON /HELP/
     $     GAMMA,GAMM1,GMM1IN,GMqGM1,CA,CAINV,ALCON,PI,PI4,PI4INV,
     $     GMINV,GM1,GP1,GM1INV,WIDTH,GGM1,GGM1IN,
     $     PI8INV,TWELVI,RUSVIS,EVMAX,BETMAX,ALPHA,BETA,
     $     PDMB,ORDER,ORD2,BDIF,NXSTEP,NYSTEP,
     $     FCOEF(20),DCOEF(20),DROP(20),
     $     NSTART,NSTOP,NDUMP,NVAXDUMP,NDIAGN,CRNT,NPICT,
     $     nplotit,nionplt,nprint,
     $     nicool,cool,
     $     niondump
        DIMENSION FDUM(LIP1,LJqq)
C
        if ( isweep .eq. 1 ) then
C
        do 200 j=1,ord2
        do 200 i=1,mip1
        fdum(i,j) = fdum(i,j+mj)
  200   fdum(i,mj+ord2+j) = fdum(i,j+ord2)
C
        else if ( isweep .eq. 2 ) then
C
          if ( lvar .gt. 1 ) then
C
            do 300 j=1,ord2
            DO 300 I=1,MIP1
            FDUM(I,j) = -FDUM(MIP2-I,mj+j)
 300        FDUM(I,mj+ord2+j) = -FDUM(MIP2-I,ord2+j)
C
          else
C
            do 310 j=1,ord2
            DO 310 I=1,MIP1
            FDUM(I,j) = FDUM(MIP2-I,mj+j)
 310        FDUM(I,mj+ord2+j) = FDUM(MIP2-I,ord2+j)
C
          endif
        else
C
        DO 400 J=1,ORD2
        DO 400 I=1,MIP1
        FDUM(I,J) = FDUM(I,1+ORD2)
  400   FDUM(I,MJ+ORD2+J) = FDUM(I,MJ+ORD2)
C
        endif
C
        RETURN
        END

