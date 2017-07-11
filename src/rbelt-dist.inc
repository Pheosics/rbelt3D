c     rbelt-dist.inc - initial particle distribution parameters

c     may want to move flag_switch to gcenter.inc
      logical flag_switch
      integer dist_seed,flag0,initdist
      real*8 dt_dist,init_t,emin,emax,epa_min,epa_max,
     &lmin,lmax,radius,exp,factor
      common /dist_params/ dt_dist,init_t,emin,emax,epa_min,epa_max,
     &lmin,lmax,radius,exp,factor,dist_seed,initdist,flag0,flag_switch

