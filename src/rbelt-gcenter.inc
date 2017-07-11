c     rbelt-gcenter.inc - guiding center time integration parameters

      integer seed
      real*8 tstep_gc,dx_max_gc,go2lrntz,go2gc,etaswitch,dtgo2gc,
     &b_prev,t_last,x_last(3)
      real*8 eta,switch
      common /gc_step/tstep_gc,dx_max_gc,go2lrntz,go2gc,etaswitch,
     &dtgo2gc,b_prev,t_last,x_last,seed
      common /gc_check/eta,switch
