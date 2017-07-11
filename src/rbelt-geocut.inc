c     rbelt-geocut.inc - rbelt-geocut.f include file

      include 'rbelt-const.inc'
      
      integer nrgdt
      real min_rgdt,max_rgdt,res,delta,delta2,pnmbr
      real phalf,nhalf,substep,step1

      common /cutoff/nrgdt,min_rgdt,max_rgdt,res,delta,delta2,pnmbr,
     &phalf,nhalf,substep,step1

