c bzero.inc
	COMMON /BZERO/
     $     BNqFI(NOUTP1,NJ,NK),BsqqI(NOUTp1,NJ,NK),
     $     BZqnI(NOUTP1,NJ,NK),BXqnI(NOUTP1,NJ,NK),BYqnI(NOUTP1,NJ,NK),
     $     BZqFI(NOUTP1,NJ,NK),BXqFI(NOUTP1,NJ,NK),BYqFI(NOUTP1,NJ,NK),
     $     BxqK(NOUTP1,NJP1,NK),ByqK(NOUTP1,NJP1,NK),
     $	      BZQK(NOUTP1,NJP1,NK),
     $     BNqFJ(NOUT,NJP1,NK),BsqqJ(NOUT,NJp1,NK),
     $     BZqnJ(NOUT,NJP1,NK),BXqnJ(NOUT,NJP1,NK),BYqnJ(NOUT,NJP1,NK),
     $     BZqFJ(NOUT,NJP1,NK),BXqFJ(NOUT,NJP1,NK),BYqFJ(NOUT,NJP1,NK),
     $     BxqI(NOUT,NJP1,NKP1),ByqI(NOUT,NJP1,NKP1),
     $	      BZQI(NOUT,NJP1,NKP1),
     $     BNqFK(NOUT,NJ,NKP1),BsqqK(NOUT,NJ,NKp1),
     $     BZqnK(NOUT,NJ,NKP1),BXqnK(NOUT,NJ,NKP1),BYqnK(NOUT,NJ,NKP1),
     $     BZqFK(NOUT,NJ,NKP1),BXqFK(NOUT,NJ,NKP1),BYqFK(NOUT,NJ,NKP1),
     $     BxqJ(NOUTP1,NJ,NKP1),ByqJ(NOUTP1,NJ,NKP1),
     $	      BZQJ(NOUTP1,NJ,NKP1)
       parameter (nbtot = 8*noutp1*nj*nk + 3*noutp1*njp1*nk +
     $      8*nout*njp1*nk + 3*nout*njp1*nkp1 +
     $ 		8*nout*nj*nkp1 + 3*noutp1*nj*nkp1 )
           dimension bzeros(nbtot)
*           equivalence (bnqfi,bzeros)
           equivalence (bnqfi,azeros)
