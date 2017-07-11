
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        subroutine flux2ave
#include "param.inc"
#include "meter.inc"
#include "var.inc"
C 
        DO 340 K=1,NK
        DO 340 J=1,NJ
        DO 340 I=1,NIP1
  340   BI(I,J,K) = BI(I,J,K)/FACEqI(I,J,K)
C
        DO 350 K=1,NK
        DO 350 J=1,NJP1
        DO 350 I=1,NI
  350   BJ(I,J,K) = BJ(I,J,K)/FACEqJ(I,J,K)
C
        DO 360 K=1,NKP1
        DO 360 J=1,NJ
        DO 360 I=1,NI
  360   BK(I,J,K) = BK(I,J,K)/FACEqK(I,J,K)
C
        return
        end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        subroutine ave2flux
#include "param.inc"
#include "meter.inc"
#include "var.inc"
C 
        DO 340 K=1,NK
        DO 340 J=1,NJ
        DO 340 I=1,NIP1
  340   BI(I,J,K) = BI(I,J,K)*FACEqI(I,J,K)
C
        DO 350 K=1,NK
        DO 350 J=1,NJP1
        DO 350 I=1,NI
  350   BJ(I,J,K) = BJ(I,J,K)*FACEqJ(I,J,K)
C
        DO 360 K=1,NKP1
        DO 360 J=1,NJ
        DO 360 I=1,NI
  360   BK(I,J,K) = BK(I,J,K)*FACEqK(I,J,K)
C
        return
        end
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
