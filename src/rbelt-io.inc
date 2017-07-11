c     rbelt-io.inc - i/o parameters

c     for particle data output array
      integer num_flts,max_yout
      parameter(num_flts=8)
      parameter(max_yout=10)
      real yout(num_flts,max_yout)
      integer num_ints
      parameter(num_ints=2)
      integer int_yout(num_ints,max_yout)
      integer youtstep
      common /output_array/ yout,int_yout,youtstep

c     for data output & integration halt times

      logical binio,print_info,flux_out,dist_out,prcp_out,init_out,
     &cone_out,bfld_out

      integer wstep,num_wsteps,wlines_dist,wlines_flux,wspntr,halt2step,
     &max_fcount,wlines_prcp,wlines_cone,lcalc,eqmap

      integer max_wsteps
      parameter(max_wsteps=100000001)
      real wtime(max_wsteps),dthalt,flux_dt,init_twrite,twrite,dtwrite

      common /write_time/wtime,dthalt,flux_dt,init_twrite,twrite,
     &dtwrite,wstep,num_wsteps,wlines_dist,wlines_flux,wspntr,halt2step,
     &max_fcount,wlines_prcp,wlines_cone,binio,print_info,flux_out,
     &dist_out,prcp_out,init_out,cone_out,bfld_out,lcalc,eqmap

      real x_prev,y_prev,z_prev
      common /pdata/ x_prev,y_prev,z_prev
