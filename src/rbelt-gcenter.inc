c     rbelt-gcenter.inc - guiding center time integration parameters

      integer seed
      real tstep_gc,dx_max_gc,go2lrntz,etaswitch,go2gc,dtgo2gc,
     &b_prev,t_last,x_last(3)
      real eta,switch
      common /gc_step/tstep_gc,dx_max_gc,go2lrntz,etaswitch,go2gc,
     &dtgo2gc,seed,b_prev,t_last,x_last
      common /gc_check/eta,switch
