c pa.inc
      COMMON /paVAR/ pa(NI,NJ,NK)
      real PAVAR(NI*NJ*NK)
*      DIMENSION PAVAR(NI*NJ*NK)
      EQUIVALENCE (PAVAR,pa)
