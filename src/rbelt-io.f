c     rbelt-io.f - i/o routines

************************************************************************

      subroutine io_init(basename,firstfilenum)

c     initialize i/o routines for all particles

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-lorentz.inc'
      include 'rbelt-gcenter.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-ut.inc'

      integer i
      real time
      character*80 basename,fluxfile,distfile,filename,string_out,
     &prcpfile,conefile,infofile,lfile,bfile,goesfile
      integer firstfilenum,lnblnk
      NAMELIST /io/ dthalt,flux_dt,init_twrite,dtwrite,binio,
     &print_info,flux_out,dist_out,prcp_out,init_out,cone_out,bfld_out,
     &lcalc,eqmap

c     read in and normalize input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,io)
      CLOSE (81)
c     normalize
      dthalt=dthalt*tfactor
      flux_dt=flux_dt*tfactor
      init_twrite=init_twrite*tfactor+tzero
      dtwrite=dtwrite*tfactor
      if (flux_out.eqv..false.) flux_dt=0.

      print *
      print *,'*** in subroutine io_init ***'

c     make sure y0 is declared properly
      if ((init_out.eqv..false.).and.(num_pdata.lt.10)) then
         print *,'if init_out = .true.'
         print *,'then num_pdata in y0.inc must be >= 10'
         stop
      endif

c     set up write out times array
c     i.e. times (from begining of run) to write out particle data
      num_wsteps=-1

10       num_wsteps=num_wsteps+1
         if ((num_wsteps+2).gt.max_wsteps) then 
            print *,'max_wsteps,tmax/dtwrite=',max_wsteps,tmax/dtwrite
            print *,'num_wsteps+2 > max_wsteps'
            print *,'make max_wsteps larger,'
            print *,'or make dtwrite larger,'
            print *,'or make init_twrite larger,'
            print *,'or make tmax smaller.'
            stop
         endif
         wtime(num_wsteps+1)=init_twrite+num_wsteps*dtwrite
*******************************************************
*         print *,'   wstep,wtime (abs. & from start)=',
*        &num_wsteps+1,wtime(num_wsteps+1)/tfactor,
*        &(wtime(num_wsteps+1)-tzero1)/tfactor
*******************************************************
      if (wtime(num_wsteps+1).le.tmax) goto 10

      print *,'1st wstep, wtime (abs. & from start) =',
     &1,wtime(1)/tfactor,(wtime(1)-tzero1)/tfactor
      print *,'last wstep, wtime (abs. & from start) =',
     &num_wsteps,wtime(num_wsteps)/tfactor,
     &(wtime(num_wsteps)-tzero1)/tfactor
      print *,'delta t to write particle data (sec) =',dtwrite/tfactor
      print *,'num. of wsteps=',num_wsteps
      print *
      print *,'note, delta t to write particle data may be larger'
      print *,'and num. of wsteps smaller if particle integration'
      print *,'step size is > delta t for writing  particle data'
      print *

*      stop

c     open output files
      wlines_dist=0
      wlines_flux=0
      wlines_prcp=0
      wlines_cone=0

c     open particle distribution data file
      if (dist_out.eqv..true.) then
         filename=distfile('',dist_seed0)
         print *,'opening ',filename(1:lnblnk(filename))
         if (binio.eqv..true.) then
            open (12,file=filename(1:lnblnk(filename)),
     &      form='unformatted')
         else
            open (12,file=filename(1:lnblnk(filename)))
         endif
      endif

c     open equatorial flux data file
      if (flux_out.eqv..true.) then
         filename=fluxfile(basename,dist_seed0)
         print *,'opening ',filename(1:lnblnk(filename))
         if (binio.eqv..true.) then
            open (14,file=filename(1:lnblnk(filename)),
     &         form='unformatted')
         else
            open (14,file=filename(1:lnblnk(filename)))
         endif
      endif

c     open precipitation (at spherical inner boundary) data file
      if (prcp_out.eqv..true.) then
         filename=prcpfile(basename,dist_seed0)
         print *,'opening ',filename(1:lnblnk(filename))
         if (binio.eqv..true.) then
            open (20,file=filename(1:lnblnk(filename)),
     &         form='unformatted')
         else
            open (20,file=filename(1:lnblnk(filename)))
         endif
      endif

c     open Stormer/allowed cone data file
      if (cone_out.eqv..true.) then
         filename=conefile(basename,dist_seed0)
         print *,'opening ',filename(1:lnblnk(filename))
         if (binio.eqv..true.) then
            open (26,file=filename(1:lnblnk(filename)),
     &         form='unformatted')
         else
            open (26,file=filename(1:lnblnk(filename)))
         endif
      endif

c     open B field file
      if (bfld_out.eqv..true.) then
         filename=bfile(basename)
         print *,'opening ',filename(1:lnblnk(filename))
         open (28,file=filename(1:lnblnk(filename)))
      endif

c     open L-shell file
      if (lcalc.ne.0) then
         filename=lfile(basename)
         print *,'opening ',filename(1:lnblnk(filename))
         if (binio.eqv..true.) then
            open (30,file=filename(1:lnblnk(filename)),
     &         form='unformatted')
         else
            open (30,file=filename(1:lnblnk(filename)))
         endif
      endif

c     open GOES to SM equitorial mapping file
      if (eqmap.ne.0) then
         filename=goesfile(basename)
         print *,'opening ',filename(1:lnblnk(filename))
         open (32,file=filename(1:lnblnk(filename)))
      endif

c     open run info file 
      filename=infofile(basename,dist_seed0)
      print *,'opening ',filename(1:lnblnk(filename))
      open (24,file=filename(1:lnblnk(filename)))
c     from rbelt-grid.inc
      write (24,40) 'nx= ',nx
      write (24,40) 'ny= ',ny
      write (24,40) 'nz= ',nz
      write (24,40) 'nt= ',nt
c     from rbelt-const.inc
      write (24,42) 'charge= ',charge
      write (24,42) 'charge_sign= ',charge_sign
      write (24,42) 'm0= ',m0
      write (24,42) 'b0= ',b0/ffactor
c     from rbelt-y0.inc
      write (24,40) 'num_particles= ',num_particles
c     from rbelt-io.inc
      write (24,40) 'max_yout= ',max_yout
      write (24,40) 'num_flts= ',num_flts
      write (24,40) 'num_ints= ',num_ints
      write (24,40) 'max_wsteps= ',max_wsteps
c     set in rbelt-input.txt
c     rbelt namelist
      write (24,44) 'basename= ',basename
      write (24,40) 'firstfilenum= ',firstfilenum
c     bounds namelist
      write (24,42) 'rmin= ',rmin
      write (24,42) 'rmax= ',rmax
      write (24,42) 'tmax= ',tmax/tfactor
c     io namelist
      write (24,42) 'dthalt= ',dthalt/tfactor
      write (24,42) 'flux_dt= ',flux_dt/tfactor
      write (24,42) 'init_twrite= ',init_twrite/tfactor
      write (24,42) 'dtwrite= ',dtwrite/tfactor
      write (24,46) 'binio= ',binio
      write (24,46) 'flux_out= ',flux_out
      write (24,46) 'dist_out= ',dist_out
      write (24,46) 'prcp_out= ',prcp_out
      write (24,46) 'init_out= ',init_out
      write (24,46) 'cone_out= ',cone_out
      write (24,40) 'lcalc= ',lcalc
c     Lorentz namelist
      write (24,42) 'tstep_lrntz= ',tstep_lrntz/tfactor
      write (24,42) 'dx_max= ',dx_max
c     guiding center namelist
      write (24,42) 'tstep_gc= ',tstep_gc
      write (24,42) 'dx_max_gc= ',dx_max
      write (24,42) 'go2lrntz= ',go2lrntz
      write (24,42) 'etaswitch= ',etaswitch
      write (24,42) 'go2gc= ',go2gc
      write (24,42) 'dtgo2gc= ',dtgo2gc/tfactor
      write (24,40) 'seed= ',seed
c     dist namelist
      write (24,40) 'dist_seed0= ',dist_seed0
      write (24,42) 'dt_dist= ',dt_dist/tfactor
      write (24,42) 'init_t= ',init_t/tfactor
      write (24,42) 'emin= ',emin*wrest
      write (24,42) 'emax= ',emax*wrest
      write (24,42) 'epa_min= ',epa_min*raddeg
      write (24,42) 'epa_max= ',epa_max*raddeg
      write (24,42) 'lmin= ',lmin
      write (24,42) 'lmax= ',lmax
      write (24,42) 'radius= ',radius
      write (24,42) 'exp= ',exp
      write (24,42) 'factor= ',factor
      write (24,40) 'flag0= ',flag0
      write (24,46) 'flag_switch= ',flag_switch
      write (24,40) 'initdist= ',initdist
c     declaired in rbelt-ut.inc
      write (24,40) 'year0= ',year0
      write (24,40) 'doy0= ',doy0
      write (24,40) 'hour0= ',hour0
      write (24,40) 'min0= ',min0
      write (24,40) 'sec0= ',sec0
c     declaired in rbelt-io.inc
40    format (a16,i10)
42    format (a16,e12.4)
44    format (a16,a80)
46    format (a16,l1)
*      close(24)

c     write particle initial conditions to file
      if (init_out.eqv..true.) call write_init_cond(basename)

      return
      end

************************************************************************

      subroutine write_init_cond(basename)

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-const.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-bounds.inc'
      logical lmonly
      integer n,j,errcode,lnblnk
      real t,y(6),ke,alpha,eta,Lm,lshell,leI0,bloc,bmir,bmin
      real pangle_lorentz,pangle_gc
      character*80 basename,filename,string_out,initfile

      print *
      print *,'*** in subroutine write_init_cond ***'

      if ((init_t+dt_dist).gt.tgrmax) then
         print *,'can not get initial pitch angles of particles with
     &    with t > tgrmax, need fields'
         print *,'init_t+dt_dist must not be greater than tgrmax'
         print *,'init_t+dt_dist, tgrmax=',
     &    (init_t+dt_dist)/tfactor,tgrmax/tfactor
         stop
      endif

c     open initial conditions file
      filename=initfile(basename,dist_seed0)
      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (22,file=filename(1:lnblnk(filename)),form='unformatted')
      else
         open (22,file=filename(1:lnblnk(filename)))
      endif

      do n = 1,num_particles

c        y(6) input used below
         do j = 1,6
            y(j) = y0(j,n)
         enddo
c        *** statement below may need to be modified to include new ref. time
         t=y0(7,n)-tzero
         call get_fields(y,t)
         if (int_y0(2,n) .eq. 0) then
c           particle kinetic energy
            ke=sqrt(1.+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))-1.
c           pitch angle
c           (get_fields at current y,t above)
            alpha = pangle_lorentz(y,t)
         elseif (int_y0(2,n) .eq. 1) then
c           kinetic energy
            ke=(sqrt(1.+2.*abs(b)*mu+y(4)*y(4))-1.)*wrest
c           pitch angle
c           (get_fields at current y,t above)
            alpha = pangle_gc(y,t)*raddeg
c           comment out gc2lrntz below to get GC position in output file
            mu=y0(6,n)
            call gc2lrntz(y,t)
         else
            print *,'flag not 0 or 1 in subroutine write_init_cond'
            stop
         endif
c        compute L-shell, 2nd invatiant, bloc, bmir, bmin
         if (lcalc.eq.0) then
            lshell=sqrt(y(1)**2+y(2)**2+y(3)**2)/
     &      ((y(1)**2+y(2)**2)/(y(1)**2+y(2)**2+y(3)**2))
         elseif (lcalc.ge.1) then
c           lstar (this will slow things down considerably)
            lmonly=.true.
            call rbelt_lshell
     &      (y,t,alpha,lmonly,Lm,lshell,leI0,bloc,bmir,bmin,errcode)
         endif

*c        may also wish to include 2nd invar. from above and...
*c        angle between velocity and normal to x-y plane
*         eta=atan2(sqrt(y(4)*y(4)+y(5)*y(5)),y(6))*raddeg
*c        local B -- convert to nT
*	  b/ffactor/ntg

         if (binio.eqv..true.) then
            write (22) y(1),y(2),y(3),lshell,ke*wrest,alpha*raddeg,
     &      y0(7,n)/tfactor
         else
            write (22,10) y(1),y(2),y(3),lshell,ke*wrest,alpha*raddeg,
     &      y0(7,n)/tfactor
10          format (7e12.4)
         endif

      enddo

      close(22)
      return
      end

************************************************************************

      subroutine yout_init(i,y,t,t2,thalt2,thalt3)

c     initialize i/o for each particle
c     called right before time loop

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-bounds.inc'
      integer i
      real y(6),t,t2,thalt2,thalt3

*      print *
*      print *,'*** in subroutine yout_init ***'

c     initialize yout step & write step
      youtstep=1
      wstep=1
100   if ((t.gt.(wtime(wstep)-tzero-flux_dt)).and.
     &   (wstep.le.num_wsteps+1)) then
*         print *,'t,wstep,wtime(wstep)=',
*     &   t/tfactor,wstep,(wtime(wstep)-tzero)/tfactor
         wstep=wstep+1         
         goto 100
      endif

c     initialize thalt2
c     halt integration is flux_dt before write time
      if (flux_out.eqv..true.) then
         thalt2=wtime(wstep)-tzero-flux_dt
      elseif ((dist_out.eqv..true.).or.(print_info.eqv..true.)) then
         thalt2=wtime(wstep)-tzero
      else
c        past last write time
c        (but note that now we return to time loop)
*         thalt2=(wtime(num_wsteps)-tzero)+1.
         thalt2=tmax-tzero+1.
      endif

c     more initializations
      halt2step=wstep-1
      x_prev= 0.0
      y_prev= 0.0
      z_prev= 0.0

      return
      end

************************************************************************

      subroutine io(i,y,dydt,t,t2,dt,thalt2)

c     called when t>thalt2 to write particle data
c     dthalt is how often to check for equatorial crossings 
c     during flux integration interval flux_dt
c
c     flux_dt should be shorter than smalest bounce period, 
c     so that we do not over sample a subset of the 
c     trajectories with small bounce periods.

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-mu.inc'
      integer i
      real y(6),dydt(6),t,t2,dt,thalt2

*      print *
*      print *,'*** in subroutine io ***'
*      print *,'t,wstep,wtime,thalt2=',
*     &t/tfactor,wstep,(wtime(wstep)-tzero)/tfactor,thalt2/tfactor

c     for writing to particle distribution output file 
c     if we passed a write time, write particle data to output vector

*      print *,'t,(wtime(wstep)-tzero)=',
*     &t/tfactor,(wtime(wstep)-tzero)/tfactor

      if (t.ge.(wtime(wstep)-tzero)) then
         if (print_info.eqv..true.) call printinfo(i,y,t,dt)
         if (dist_out.eqv..true.) call dist2yout(i,y,t,dt)
         if (wstep.lt.num_wsteps+2) wstep=wstep+1
      endif

c     for writing to flux output file
      if ((flux_out.eqv..true.).and.(wstep.le.num_wsteps)) then
c     if we are within flux_dt of next write time check for equatorial crossings
         if ((t.ge.(wtime(wstep)-tzero-flux_dt)).and.
     &   (t.le.(wtime(wstep)-tzero))) then
c           first we need to set z_prev = y(3) (after if-endif below)
c           then we update halt2step to the current wstep and can start 
c           checking for equitorial crossings
            if (halt2step.eq.wstep) then 
*               if (y(1)*x_prev.lt.0.0) call data2yout(i,y,t,dt,3)
*               if (y(2)*y_prev.lt.0.0) call data2yout(i,y,t,dt,2)
c              if particle has crossed the equitorial plane
c              (for more accuracy data2yout should be modified to get the
c              quantities right at the eqitorial crossing)
               if (y(3)*z_prev.lt.0.0) call flux2yout(i,y,t,dt)
            endif
*            x_prev= y(1)
*            y_prev= y(2)
            z_prev= y(3)

            halt2step=wstep
c           set new halt time
*            thalt2=amin1(thalt2+dthalt,wtime(wstep)-tzero)
            thalt2=amin1(t+dthalt,wtime(wstep)-tzero)
         else
            thalt2=wtime(wstep)-tzero-flux_dt
         endif
      else
         thalt2=wtime(wstep)-tzero
      endif

      return
      end

************************************************************************
******* MODIFY PARTICLE OUTPUT DATA BELOW IN 3 SUBROUTINES BELOW *******
************************************************************************


      subroutine dist2yout(i,y,t,dt)

c     THIS ROUTINE NEEDS WORK
C
c     write particle data to yout array
c     currently we modify this routine to get desired output values
c     (e.g., x-y-z position, pitch angle, L-shell etc.) --
c     needs to be modified to output data requested in input file

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-status.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-gcenter.inc'
      integer i,j,errcode
      real t,y(6),ytmp(6),dt
      real ke,pa,Lm,lshell,leI0,bloc,bmir,bmin,beqtr
      real pangle_lorentz,pangle_gc

*      print *
*      print *,'*** in subroutine data2yout ***'
*      print *,'wstep,youtstep,t,wtime,ofile_tag=',wstep,youtstep,
*     & real(t/tfactor),(wtime(wstep)-tzero)/tfactor

      int_yout(1,youtstep) = 0
      int_yout(2,youtstep) = wstep
c     for dist. file
      if (wstep .le. num_wsteps) then
         do j = 1,6
            ytmp(j) = y(j)
         enddo
c        if Lorentz trajectory
         if (flag .eq. 0) then
c           first compute particle data:
c           particle kinetic energy
            ke=sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))-1
c           pitch angle
c           fields at curent y,t should be correct 
c           (see RK4 integrator)
            pa = pangle_lorentz(y,t)
c        else if guiding center trajectory
         elseif (flag .ge. 1) then
c           energy
            ke=sqrt(1.+2.*abs(b)*mu+y(4)*y(4))-1.
c           pitch angle
c           fields at curent y,t should be correct 
c           (see RK4 integrator)
            pa = pangle_gc(y,t)
c           get particle position
            call gc2lrntz(ytmp,t)
         endif

c        compute L-shell, 2nd invatiant, bloc, bmir, bmin
         if (lcalc.eq.0) then
            lshell=sqrt(ytmp(1)**2+ytmp(2)**2+ytmp(3)**2)/
     &      ((ytmp(1)**2+ytmp(2)**2)/(ytmp(1)**2+ytmp(2)**2+ytmp(3)**2))
*         elseif (lcalc.ge.1) then
*c           lstar (this will slow things down considerably)
*            lmonly=.true.
*            call rbelt_lshell
*     &      (ytmp,t,pa,lmonly,Lm,lshell,leI0,bloc,bmir,bmin,errcode)
         endif

*c        trace field line to equatorial plane
*c        NEED TO CREATE A trace2eqtr OPTION
*c        THIS MODIFIES B!!
*         if (map.eq..true.)
*            call trace2eqtr(t,ytmp,beqtr,errcode)
*            if (errcode.ne.0) then
*               print *,'trace2eqtr fail: errcode=',errcode
*               beqtr=0.0
**              stop
*            endif
*         endif

c        put particle data in output array
c        x,y,z position
         do j = 1,3
            yout(j,youtstep) = ytmp(j)
         enddo
c        energy
         yout(4,youtstep)=ke*wrest
c        pitch angle
         yout(5,youtstep)=pa*raddeg
!         yout(5,youtstep)=switch
c        local B (convert to nT)
c        fields at curent y,t must be correct 
         yout(6,youtstep)=lshell
!         yout(6,youtstep)=eta
        print*,ke*wrest,mu,abs(b),y(4)
*c       |B| at equator
*c       beqtr returned by trace2eqtr has dimensions of nT!
*        yout(7,youtstep)=beqtr
*c       2nd invatiant
*        yout(5,youtstep)=leI0
*c       L-star
*        yout(6,youtstep)=lshell

         if (youtstep.lt.max_yout) then 
            youtstep=youtstep+1
         else
            call filewrite(i)
            youtstep=1
         endif
      endif

      return
      end

************************************************************************


      subroutine flux2yout(i,y,t,dt)

c     write particle data to yout array
c     currently we modify this routine to get desired output values
c     (e.g., x-y-z position, pitch angle, L-shell etc.)

c     fields at current y,t should be correct 
c     (see time_loop routine and RK4 integrator)

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-status.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-mu.inc'
      logical lmonly
      integer i,j,ofile_tag,errcode
      real t,y(6),ytmp(6),dt
      real ke,pa,eta,Lm,lshell,leI0,bloc,bmir,bmin
      real pangle_lorentz,pangle_gc

*      print *
*      print *,'*** in subroutine data2yout ***'
*      print *,'wstep,youtstep,t,wtime,ofile_tag=',wstep,youtstep,
*     & real(t/tfactor),(wtime(wstep)-tzero)/tfactor,ofile_tag

      int_yout(1,youtstep) = 1
      int_yout(2,youtstep) = wstep

c     y(6) input used below
      do j = 1,6
         ytmp(j) = y(j)
      enddo
c     if Lorentz trajectory
      if (flag .eq. 0) then
c        first compute particle data:
c        particle kinetic energy
         ke=sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))-1
c        pitch angle
         pa = pangle_lorentz(y,t)
c     else if guiding center trajectory
      elseif (flag .ge. 1) then
c        energy
         ke=sqrt(1.+2.*abs(b)*mu+y(4)*y(4))-1.
c        pitch angle
         pa = pangle_gc(y,t)
c        get particle position
         call gc2lrntz(ytmp,t)
      endif
*c     compute L-shell, 2nd invatiant, bloc, bmir, bmin
*      if (lcalc.eq.0) then
*         lshell=sqrt(ytmp(1)**2+ytmp(2)**2+ytmp(3)**2)/
*     &   ((ytmp(1)**2+ytmp(2)**2)/(ytmp(1)**2+ytmp(2)**2+ytmp(3)**2))
*      elseif (lcalc.ge.1) then
*c        lstar (this will slow things down considerably)
*         lmonly=.true.
*         call rbelt_lshell
*     &   (ytmp,t,pa,lmonly,Lm,lshell,leI0,bloc,bmir,bmin,errcode)
*      endif
c     angle between velocity and normal to x-y plane
      eta=atan2(sqrt(ytmp(4)*ytmp(4)+ytmp(5)*ytmp(5)),ytmp(6))
c     put particle data in output array

c     load data to yout array
c     x,y,z coordinates
      do j = 1,2
         yout(j,youtstep) = ytmp(j)
      enddo
c     angle between velocity and normal to x-y plane
      yout(3,youtstep)=eta*raddeg
c     particle kinetic energy
      yout(4,youtstep)=ke*wrest
c     pitch angle
      yout(5,youtstep)=pa*raddeg
c     local B (convert to nT)
      yout(6,youtstep)=b/ffactor/ntg

*c     2nd invatiant
*      yout(8,youtstep)=leI0
*c     L-star
*      yout(6,youtstep)=lshell

      if (youtstep.lt.max_yout) then 
         youtstep=youtstep+1
      else
         call filewrite(i)
         youtstep=1
      endif

      return
      end

************************************************************************

      subroutine filewrite(i)

c     write yout array to output file
c     note, consider using independent logical var. instead of init_out below

c     currently we modify this routine to get desired output values
c     (e.g., x-y-z position, pitch angle, L-shell etc.) --
c     needs to be modified to output data requested in input file

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-y0.inc'
      integer i,j
      real t,y(6)
      real energy,r

*      print *
*      print *,'*** in subroutine filewrite ***'
*      print *,'youtstep=',youtstep

c     write data
      do j = 1,youtstep-1
c        write to particle distribution output file
c        i: particle number
c        int_yout(1,*): ofile_tag=0; int_yout(2,*) = wstep
         if (int_yout(1,j).eq.0) then
            wlines_dist=wlines_dist+1
            if (binio.eqv..true.) then
               if (init_out.eqv..true.) then
                  write (12) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j)
               else
                  write (12) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j),y0(8,i)*wrest,
     &            y0(9,i)*raddeg,y0(10,i)/tfactor
               endif
            else
               if (init_out.eqv..true.) then
                  write (12,10) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j)
               else
                  write (12,20) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j),y0(8,i)*wrest,
     &            y0(9,i)*raddeg,y0(10,i)/tfactor
               endif
            endif
c        write to particle flux output file
c        flux through x-y plane
c        i: particle number
c        int_yout(1,*): ofile_tag=1; int_yout(2,*) = wstep
         elseif (int_yout(1,j).eq.1) then
            wlines_flux=wlines_flux+1
            if (binio.eqv..true.) then
               if (init_out.eqv..true.) then
                  write (14) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j)
               else
                  write (14) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j),y0(8,i)*wrest,
     &            y0(9,i)*raddeg,y0(10,i)/tfactor
               endif
            else
               if (init_out.eqv..true.) then
                  write (14,30) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j)
               else
                  write (14,40) i,int_yout(2,j),yout(1,j),yout(2,j),
     &            yout(3,j),yout(4,j),yout(5,j),yout(6,j),y0(8,i)*wrest,
     &            y0(9,i)*raddeg,y0(10,i)/tfactor
               endif
            endif
         endif

10       format (2i8,6e12.4)
20       format (2i8,9e12.4)
30       format (2i8,6e12.4)
40       format (2i8,9e12.4)

*c        write to particle flux output file
*c        flux through x-y plane
*c        i: particle number
*c        int_yout(1,*): ofile_tag=1; int_yout(2,*) = wstep; yout(1-2,j): x,y
*c        yout(3-5,j): energy (MeV), pitch angle (deg.), incidence angle (deg.) (w/ normal to surface)
*         elseif (int_yout(1,j).eq.1) then
*            wlines_flux=wlines_flux+1
*            if (binio.eqv..true.) then
*               write (14) i,int_yout(2,j),1,yout(1,j),yout(2,j),
*     &         yout(3,j),yout(4,j),yout(5,j),y0(8,i),y0(9,i)
*            else
*               write (14,998) i,int_yout(2,j),1,yout(1,j),yout(2,j),
*     &         yout(3,j),yout(4,j),yout(5,j),y0(8,i),y0(9,i)
*            endif
*         endif
*c        flux through x-z plane
*c        i: particle number
*c        int_yout(1,*): ofile_tag=2; int_yout(2,*) = wstep; yout(1-2,j): x,z
*c        yout(3-5,j): energy (MeV), pitch angle (deg.), incidence angle (deg.) (w/ normal to surface)
*         elseif (int_yout(1,j).eq.2) then
*            wlines_flux=wlines_flux+1
*            if (binio.eqv..true.) then
*               write (14) i,int_yout(2,j),2,yout(1,j),yout(2,j),
*     &         yout(3,j),yout(4,j),yout(5,j)
*            else
*               write (14,998) i,int_yout(2,j),2,yout(1,j),yout(2,j),
*     &         yout(3,j),yout(4,j),yout(5,j)
*            endif
*c        flux through y-z plane
*c        i: particle number
*c        int_yout(1,*): ofile_tag=3; int_yout(2,*) = wstep; yout(1-2,j): y,z
*c        yout(3-5,j): energy (MeV), pitch angle (deg.), incidence angle (deg.) (w/ normal to surface)
*         elseif (int_yout(1,j).eq.3) then
*            wlines_flux=wlines_flux+1
*            if (binio.eqv..true.) then
*               write (14) i,int_yout(2,j),3,yout(1,j),yout(2,j),
*     &         yout(3,j),yout(4,j),yout(5,j)
*            else
*               write (14,998) i,int_yout(2,j),3,yout(1,j),yout(2,j),
*     &         yout(3,j),yout(4,j),yout(5,j)
*            endif
*         endif

      enddo

      return
      end

************************************************************************
******* MODIFY PARTICLE OUTPUT DATA BELOW IN 3 SUBROUTINES ABOVE *******
************************************************************************

      subroutine prcp_filewrite(i,y,t)

      implicit none	
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-grid.inc'
      integer i
      real y(6),t,energy,alpha,ceta,pangle_lorentz

      wlines_prcp=wlines_prcp+1
      if (flag .ge. 1) then
         call gc2lrntz(y,t)
         flag=0
      endif
      energy=(sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))-1)*wrest
      alpha=pangle_lorentz(y,t)
c     cos of angle between velocity and normal to spherical surface
      ceta=(y(1)*y(4)+y(2)*y(5)+y(3)*y(6))/
     &sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))/
     &sqrt(y(4)*y(4)+y(5)*y(5)+y(6)*y(6))

      if (binio.eqv..true.) then
         if (init_out.eqv..true.) then
            write (20)i,(t+tzero)/tfactor,y(1),y(2),y(3),energy,alpha,
     &      ceta
         else
            write (20)(t+tzero)/tfactor,y(1),y(2),y(3),energy,alpha,
     &      ceta,y0(8,i)*wrest,y0(9,i)*raddeg,y0(10,i)/tfactor
         endif
      else
         if (init_out.eqv..true.) then
            write (20,5)i,(t+tzero)/tfactor,y(1),y(2),y(3),energy,alpha,
     &      ceta
         else
            write (20,10) (t+tzero)/tfactor,y(1),y(2),y(3),energy,alpha,
     &      ceta,y0(8,i)*wrest,y0(9,i)*raddeg,y0(10,i)/tfactor
         endif
      endif
5     format (i8,7e12.4)
10    format (10e12.4)

      return
      end

************************************************************************

      subroutine cone_filewrite(i,status,x,y,z,alpha,beta)

c     Is it better to include alpha, beta & status
c     or pass them to the subroutine as arguments?

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
*      include 'rbelt-y0.inc'
      integer i,status
      real x,y,z,alpha,beta,vx,vy,vz,asink,eta,zenith
      real vx_LSM,vy_LSM,vz_LSM

      wlines_cone=wlines_cone+1

      vx=sin(alpha)*cos(beta)
      vy=sin(alpha)*sin(beta)
      vz=cos(alpha)

*      call sm_lsm(x,y,z,-vx,-vy,-vz,asink,eta,zenith)
      call sm_lsm(x,y,z,-vx,-vy,-vz,vx_LSM,vy_LSM,vz_LSM)

      if (binio.eqv..true.) then
*         write (26) status,alpha,beta,asink,eta,zenith
         write (26) status,alpha,beta,vx_LSM,vy_LSM,vz_LSM
      else
*         write (26,10) status,alpha,beta,asink,eta,zenith
         write (26,10) status,alpha,beta,vx_LSM,vy_LSM,vz_LSM
      endif
10    format (i8,5e12.4)

      return
      end

************************************************************************

      subroutine printinfo_sav(i,y,t,dt)
      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-status.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-io.inc'
      external rhs_northrop
      integer i,step
      real t,tmp
      real y(6),dt,energy
      real p_perp,p_parl,ke,gamma2,gamma,r,phi
      real lshell

      r=sqrt(y(1)**2+y(2)**2+y(3)**2)

      call get_fields(y,t)

      if (flag .eq. 0) then
         gamma2=1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6)
         gamma = sqrt(gamma2)
      elseif (flag .ge. 1) then
         ke=sqrt(1.+2.*abs(b)*y(6)+y(4)*y(4))-1.
         gamma=1+ke
      endif

      print *,'t,r,kinetic energy(MeV),mu,b,flag=',
     &t/tfactor,r,(gamma-1)*wrest,mu,b/ffactor/ntg,flag

      return
      end


************************************************************************

      subroutine printinfo(i,y,t,dt)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-status.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-io.inc'
      external rhs_northrop
      integer i,step
      real t,tmp
      real y(6),dt,energy
      real p_perp,p_parl,ke,gamma2,gamma,r,phi
      real lshell

      call get_fields(y,t)
      print *
      print *,'t,wstep,wtime=',
     &t/tfactor,wstep,(wtime(wstep)-tzero)/tfactor
      print *,'particle',i
      print *,'local t(sec)=',t/tfactor 
      print *,'dt(sec)=',dt/tfactor 
      print *,'from start t(sec)=',(t+(tzero-tzero1))/tfactor 
      print *,'global t(sec)=',(t+tzero)/tfactor 
*      print *,'t(hours)=',t/tfactor/3600
*      print *,'t(days)=',t/tfactor/3600/24
      print *,'x,y,z(Re)=',y(1),y(2),y(3)
      r=sqrt(y(1)**2+y(2)**2+y(3)**2)
      print *,'r=',sqrt(y(1)**2+y(2)**2+y(3)**2)
      print *,'theta=',atan2(sqrt(y(1)**2 + y(2)**2),y(3))*raddeg
*      if (y(2) .ge. 0) then	
*         phi=acos(y(1)/sqrt(y(1)*y(1)+y(2)*y(2)))
*      else
*         phi=(2*pi-acos(y(1)/sqrt(y(1)*y(1)+y(2)*y(2))))
*      endif
      phi=atan2(y(2),y(1))
      if (phi.lt.0.0) phi=2*pi+phi
      print *,'phi=',phi*raddeg
      print *,'dipole L shell=', sqrt(y(1)**2+y(2)**2+y(3)**2) /
     &   ( (y(1)**2+y(2)**2) / (y(1)**2+y(2)**2+y(3)**2) )

*      print *,'L-star=',lshell(y,t)
*      print *
*      stop

*      print *,'   t_init=',t_init
*      print *,'   r_init=',r_init
*      print *,'   phi_init=',phi_init
*      print *,'   avg. v_drift (km/sec)=',
*     &Re*(phi-phi_init)*(r+r_init)/2./((t-t_init)/tfactor)/100./1000.

      if (flag .eq. 0) then
         gamma2=1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6)
         gamma = sqrt(gamma2)
         print *,'gamma=',gamma
         print *,'vx,vy,vz(cm/sec)=',
     &   y(4)*c/gamma,y(5)*c/gamma,y(6)*c/gamma
         print *,'v(cm/sec)',
     &   sqrt(y(4)*y(4)+y(5)*y(5)+y(6)*y(6))*c/gamma
         print *,'proper ux,uy,uz(cm/sec)=',y(4)*c,y(5)*c,y(6)*c
         print *,'u(cm/sec)',sqrt(y(4)*y(4)+y(5)*y(5)+y(6)*y(6))*c
         print *,'momentum px,py,pz(g*cm/sec)=',
     &   y(4)*c*m0,y(5)*c*m0,y(6)*c*m0
         print *,'momentum p(g*cm/sec)=',
     &   sqrt(y(4)*y(4)+y(5)*y(5)+y(6)*y(6))*c*m0
         p_perp=sqrt((y(5)*bz-y(6)*by)**2.+
     &   (y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
         print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
         p_parl=(y(4)*bx+y(5)*by+y(6)*bz)/b
         print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
         print *,'mag. moment(g*cm^2/gauss)=',
     &   abs(p_perp*p_perp/b/2*m0*c*c*ffactor)
         print *,'kinetic energy(MeV)=',(gamma-1)*wrest
      elseif (flag .ge. 1) then
         print *,'total energy=',y(5)*wrest
         print *,'KE (from y(5))',(y(5)-1)*wrest
         ke=sqrt(1.+2.*abs(b)*y(6)+y(4)*y(4))-1.
         gamma=1+ke
         gamma2=(1+ke)**2
         print *,'gamma=',gamma
         print *,'v(cm/sec)',
     &   sqrt(gamma2-1)*c/gamma
         print *,'u(cm/sec)',sqrt(gamma2-1)*c
         print *,'momentum p(g*cm/sec)=',sqrt(gamma2-1)*c*m0
         p_parl=y(4)
         print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
         p_perp = sqrt(gamma2-1-y(4)*y(4))
         print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
         print *,'mag. moment(g*cm^2/sec^2/gauss)=',
     &   abs(y(6)*m0*c*c*ffactor)
         print *,'mag. moment(MeV/gauss)=',
     &   abs(y(6)*m0*c*c*ffactor)/1.6e-6
         print *,'kinetic energy(MeV)=',(gamma-1)*wrest
      endif
      print *,'gyro radius(Re)=',p_perp/b
      print *,'bx(nT)=',bx/ffactor/ntg,'by(nT)=',by/ffactor/ntg
      print *,'bz(nT)=',bz/ffactor/ntg,'b(nT)=',abs(b/ffactor/ntg)
      print *,'ex(V/m),ey(V/m),ez(V/m)=',ex/ffactor/vmsvcm,
     &ey/ffactor/vmsvcm,ez/ffactor/vmsvcm
      print *,'pitch angle=',atan2(p_perp,p_parl*charge_sign)*raddeg
*      print *,'alpha_eq,r_eq',alpha_eq*raddeg,r_eq

*         y0(5,n) = -2.5E-3/(m0*c*c*ffactor)
*         print *
*         print *,'mu (esu-cm)=',y0(5,n)*(m0*c*c*ffactor)
*         print *,'(note, (esu-cm) = ergs/G = g-cm^2/sec^2 /G)'
*         print *,'mu (MeV/G)=',
*     &   y0(5,n)*abs(m0*c*c*ffactor)/1.60217657e-6
*         (1MeV = 1.60217657e-6 ergs)

      return
      end

************************************************************************

      subroutine print_gcforce(i,y,t,dt)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-status.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-io.inc'
      external rhs_northrop
      integer i,step
      real t,tmp
      real y(6),dt,energy
      real p_perp,p_parl,ke,gamma2,gamma,r,phi
      real b2,gmbm,gmb2m,pgmb,p2gmb4,p2gmb3,p2gmb2
      real ax,ay,az,f_mirr,f_curv,f_exb,f_grad

      call get_fields2(y,t)
*      print *
*      print *,'t,z,r=',t/tfactor,y(3),sqrt(y(1)**2+y(2)**2+y(3)**2)

      b2=1/(b*b)
      gamma=sqrt(1.+2.*abs(b)*mu+y(4)*y(4))
      gmbm=mu/gamma/b
      gmb2m=gmbm/b
      pgmb=y(4)*gmbm/mu
      p2gmb4=y(4)*pgmb*b2/b
      p2gmb3=abs(y(4)*y(4)/gamma/b/b/b)
      p2gmb2=abs(y(4)*y(4)/gamma/b/b)

c     get the cross terms for the curvature drift
      ax=bx*dbxdx + by*dbxdy + bz*dbxdz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*bx/b
      ay=bx*dbydx + by*dbydy + bz*dbydz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*by/b
      az=bx*dbzdx + by*dbzdy + bz*dbzdz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*bz/b

      f_exb=sqrt(ex**2+ey**2+ez**2)
      f_mirr=abs(gmbm*(bx*dbdx+by*dbdy+bz*dbdz))
      f_curv=p2gmb2*sqrt(ax**2+ay**2+az**2)
      f_grad=gmbm*sqrt((by*dbdz-bz*dbdy)**2+
     &(bz*dbdx-bx*dbdz)**2+(bx*dbdy-by*dbdx)**2)

*      print *,'f_mirr,f_curv,f_grad,f_exb=',
*     &real(f_mirr),real(f_curv),real(f_grad),real(f_exb)


      return
      end

************************************************************************

      real function  pangle_lorentz(y,t)

c     returns PA in radians (normalized quantity)
c     assumes that current fields are correct for position y,t

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      real y(6),t,p_perp,p_parl

      p_parl=(y(4)*bx+y(5)*by+y(6)*bz)/b
      p_perp=sqrt((y(5)*bz-y(6)*by)**2.+
     &(y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
      pangle_lorentz=atan2(p_perp,p_parl*charge_sign)

      return
      end

************************************************************************

      real function  pangle_gc(y,t)
      
c     returns PA in radians (normalized quantity)
c     assumes that current fields are correct for position y,t

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      real y(6),t,p_perp,p_parl,gamma2

      gamma2=1.+2.*abs(b)*y(6)+y(4)*y(4)
      p_perp = sqrt(gamma2-1-y(4)*y(4))
      pangle_gc=atan2(p_perp,y(4)*charge_sign)

      return
      end


************************************************************************

      subroutine invar_chck(y,t)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-io.inc'
      include 'rbelt-flag.inc'
      real y(6),t
      real w0,gamma,p_perp,p_parl,rho,switch

      call get_fields2(y,t)

      if (flag .eq. 0) then

         p_perp=sqrt((y(5)*bz-y(6)*by)**2.+
     &   (y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
         rho = p_perp/b
*         print *,'lrntz invar_chck:t,switch=',t/tfactor,
*     &   rho*sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &   dbydx**2.+dbydy**2.+dbydz**2.+dbzdx**2.+dbzdy**2.+dbzdz**2.)/b


      elseif (flag .eq. 1) then

         w0=sqrt(1.+2.*abs(b)*y(6)+y(4)*y(4))-1.
         gamma = 1+w0
         p_perp = sqrt(w0*(2+w0)-y(4)*y(4))
         p_parl = y(4)
         rho = p_perp/b

*         print *
*         print *,'b(nT) =',b/ffactor/ntg
*         print *,'bx(nT)=',bx/ffactor/ntg
*         print *,'by(nT)=',by/ffactor/ntg
*         print *,'bz(nT)=',bz/ffactor/ntg
*         print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
*         print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
*         print *,'p(g*cm/sec)=',c*m0*sqrt(p_perp**2.+p_parl**2.)
*         print *,'p(g*cm/sec)=',c*m0*sqrt(gamma**2.-1)
*         print *,'rho=',rho
*         print *,'grad=',sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &   dbydx**2.+dbydy**2.+dbydz**2.+
*     &   dbzdx**2.+dbzdy**2.+dbzdz**2.)

         switch=sqrt(2*y(6)*(dbxdx**2.+dbxdy**2.+dbxdz**2.+
     &   dbydx**2.+dbydy**2.+dbydz**2.+
     &   dbzdx**2.+dbzdy**2.+dbzdz**2.)/b**3.)
         print *,'gc invar_chck:t,switch=',t/tfactor,switch

*         print *,'switch=',rho*sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &   dbydx**2.+dbydy**2.+dbydz**2.+
*     &   dbzdx**2.+dbzdy**2.+dbzdz**2.)/b

      endif

      return
      end

************************************************************************
*
*      subroutine  bounce_drift(y,t)
*
*      implicit none
*      include 'rbelt-fields.inc'
*      include 'rbelt-const.inc'
*      include 'rbelt-io.inc'
*      include 'rbelt-mu.inc'
*      include 'rbelt-flag.inc'
*      real y(6),t,gamma,p,phi
*
*      bcount=bcount+1
*      if (bcount.eq.1) then
*         t_prev=t
*         spb=0
*      endif
*      print *,'eq. cross,t,alpha_eq=',
*     &bcount,t/tfactor,alpha_eq
**      print *,'bcount,t_prev=',bcount,t_prev/tfactor
*      if (bcount.eq.3) then
*         if (flag .eq. 0) then
*            gamma = sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))
*         elseif (flag .eq. 1) then
*            gamma=sqrt(1.+2.*abs(b)*mu+y(4)*y(4))
*         endif
*         p=sqrt(gamma*gamma-1)
*         print *,' Bounce period=',(t-t_prev)/tfactor
*         print *,' Calc. period (sec.) (Schultz & Lanzerrotti,1974)=',
*     &   4*gamma*r_eq*(1.3802-0.3198*
*     &   (sin(alpha_eq)+sqrt(sin(alpha_eq))))/p/tfactor
**         print *,' Calc. period (sec.) (Hamlin etal.,1961)=',
**     &   4*gamma*r_eq*(1.30-0.56*sin(alpha_eq))/p/tfactor
*         print *,' steps per bounce=',spb
*         if (y(2) .ge. 0) then	
*            phi=acos(y(1)/sqrt(y(1)*y(1)+y(2)*y(2)))
*         else
*            phi=(2*pi-acos(y(1)/sqrt(y(1)*y(1)+y(2)*y(2))))
*         endif
*         print *,' phi =',phi*raddeg
*         print *,' Calc. phi (Nicholls & Storey, 1999) =',
*     &   -3*p*p*r_eq*raddeg*(0.7+0.3*sin(alpha_eq))*t/
*     &   gamma/2/b0
*         print *,' energy=',(gamma-1)*wrest
*         print *
**         stop
*      endif
*      if (bcount.eq.3) bcount=0
*
*      return
*      end
*
************************************************************************

      subroutine io_close()

c     i/o finalize and post processing

      implicit none	
      include 'rbelt-io.inc'

      print *
      print *,'*** in subroutine fin_output ***'
      print *,'dist. file lines writen=',wlines_dist
      print *,'flux file lines writen=',wlines_flux
      print *,'prcp file lines writen=',wlines_prcp
      print *,'cone file lines writen=',wlines_cone


c     finish up with info file & close
c     need to fix this to write the number of steps actually used
      write (24,10) 'num_wsteps= ',num_wsteps
      write (24,10) 'wlines_dist= ',wlines_dist
      write (24,10) 'wlines_flux= ',wlines_flux
      write (24,10) 'wlines_prcp= ',wlines_prcp
      write (24,10) 'wlines_cone= ',wlines_cone
      close(24)

10    format (a16,i8)
11    format (a16,e12.4)
12    format (a16,a80)
13    format (a16,l1)

      if (dist_out.eqv..true.) close(12)
      if (flux_out.eqv..true.) close(14)
      if (prcp_out.eqv..true.) close(20)
      if (cone_out.eqv..true.) close(26)
      if (lcalc.ne.0) close(28)

      return
      end

********************************************************************************
*
*      real function smth3d(i,j,k,nx,ny,nz,f)
*
*c     returns zero if point lies outside of grid
*
*      implicit none 
*      integer nx,ny,nz,i,j,k,ip,jp,kp,im,jm,km
*      real f(nx,ny,nz)
*
*      ip=i+1
*      jp=j+1
*      kp=k+1
*      im=i-1
*      jm=j-1
*      km=k-1
*
*      if (((i.eq.1).or.(i.eq.nx)).and.
*     &((j.eq.1).or.(j.eq.ny)).and.
*     &((k.eq.1).or.(k.eq.nz))) then
*
*      call corner()
*
*      elseif (((i.eq.1).or.(i.eq.nx)).and.
*     &((j.eq.1).or.(j.eq.ny)))
*
*      elseif (((j.eq.1).or.(j.eq.ny)).and.
*     &((k.eq.1).or.(k.eq.nz)))
*
*      elseif (((i.eq.1).or.(i.eq.nx)).and.
*     &((k.eq.1).or.(k.eq.nz)))
*
*      elseif ((i.eq.1).or.(i.eq.nx))
*
*      elseif ((j.eq.1).or.(j.eq.ny))
*
*      elseif ((k.eq.1).or.(k.eq.nz))
*
*      else
*      ip=i+1
*      jp=j+1
*      kp=k+1
*      im=i-1
*      jm=j-1
*      km=k-1
*      smth3d = (f(ip,j,k)+f(i,jp,k)+f(i,j,kp)+
*     &f(im,j,k)+f(i,jm,k)+f(i,j,km))*.072956795+
*     &(f(ip,jp,k)+f(i,jp,kp)+f(ip,j,kp)+
*     &f(im,jm,k)+f(i,jm,km)+f(im,j,km)+
*     &f(ip,jm,k)+f(i,jp,km)+f(ip,j,km)+
*     &f(im,jp,k)+f(i,jm,kp)+f(im,j,kp))*.025794122+
*     &(f(ip,jp,km)+f(im,jp,kp)+f(ip,jm,kp)+ 
*     &f(im,jm,kp)+f(ip,jm,km)+f(im,jp,km)+ 
*     &f(ip,jp,kp)+f(im,jm,km))*.031591219
*      return 
*      END
*

************************************************************************

      subroutine print_dist(basename,firstfilenum)

c     i/o finalize and post processing

      implicit none	
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-bounds.inc'
      integer i,j,k,firstfilenum,lnblnk,pnum,step,np,nts,num
      character*80 basename,distfile,fluxfile,filename,string_out
      real x,y,z,r,t,pvx,pvy,pvz,ke,alpha,p_parl,p_perp,calpha

      integer num_rbins,num_cabins
      parameter(num_rbins=38,num_cabins=8)
      real camax,camin
      integer r0bin(num_rbins),ca0bin(num_cabins)
      real rbwidth,cabwidth,lbound,ubound

      print *
      print *,'*** in subroutine show_dist ***'

      camin=-1.0
      camax=1.0
      cabwidth=(camax-camin)/num_cabins
      do i = 1,num_cabins
         ca0bin(i) = 0
      enddo

*      rmin=0
*      rmax=1.0
      rbwidth=(rmax-rmin)/num_rbins
      do i = 1,num_rbins
         r0bin(i) = 0
      enddo

      filename=distfile(basename,dist_seed0)
      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (16,file=filename(1:lnblnk(filename)),form='unformatted')
         read (16) np
         read (16) nts
      else
         open (16,file=filename(1:lnblnk(filename)))
         read (16,20) np
         read (16,20) nts
      endif

      do i=1,nts
         if (binio.eqv..true.) then
            read (16) num,t
         else
            read (16,30) num,t
         endif
         do j = 1,num_cabins
            ca0bin(j) = 0
         enddo
         do j = 1,num_rbins
            r0bin(j) = 0
         enddo
         do k=1,num
            if (binio.eqv..true.) then
               read (16) pnum,x,y,z,pvx,pvy,pvz
            else
               read (16,40) pnum,x,y,z,pvx,pvy,pvz
            endif
            p_parl=(pvx*bx+pvy*by+pvz*bz)/b
            p_perp=sqrt((pvy*bz-pvz*by)**2.+
     &      (pvz*bx-pvx*bz)**2.+(pvx*by-pvy*bx)**2.)/b
            calpha=cos(atan2(p_perp,p_parl))
            r=sqrt(x**2+y**2+z**2)
            do j = 1,num_cabins
               lbound=camin+(j-1)*cabwidth
               ubound=camin+j*cabwidth
               if ((calpha.ge.lbound).and.
     &         (calpha.lt.ubound))ca0bin(j) = ca0bin(j)+1
            enddo
            do j = 1,num_rbins
               lbound=rmin+(j-1)*rbwidth
               ubound=rmin+j*rbwidth
               if ((r.ge.lbound).and.
     &         (r.lt.ubound)) r0bin(j) = r0bin(j)+1
            enddo
         enddo
         print *
         print *,' ca dist. for wtime=',wtime(i)/tfactor
         do j = 1,num_cabins
            lbound=camin+(j-1)*cabwidth
            ubound=camin+j*cabwidth
            print *,j,lbound,ubound,ca0bin(j)
         enddo
         print *
         print *,' r dist. for wtime=',wtime(i)/tfactor
         do j = 1,num_rbins
            lbound=rmin+(j-1)*rbwidth
            ubound=rmin+j*rbwidth
            print *,j,lbound,ubound,r0bin(j),
     & r0bin(j)/(1.33333333*pi*(ubound**3-lbound**3))
         enddo
      enddo

20    format (i8)
30    format (1i8,1e12.4)
40    format (1i8,6e12.4)

      return
      end

************************************************************************

      subroutine print_prcp(basename,firstfilenum)

c     i/o finalize and post processing

      implicit none	
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-dist.inc'

      integer firstfilenum,lnblnk,num,i,j,k,wsteps
      character*80 basename,prcpfile,filename,string_out

      real x,y,z,r,ke,alpha,ceta,ctheta,t

      integer num_ctbins
      parameter(num_ctbins=240)
      real ctmax,ctmin
      integer ct0bin(num_ctbins)
      real ctbwidth,lbound,ubound

      print *
      print *,'*** in subroutine show_prcp ***'

      ctmin=-1.0
      ctmax=1.0
      ctbwidth=(ctmax-ctmin)/num_ctbins
      do i = 1,num_ctbins
         ct0bin(i) = 0
      enddo

      filename=prcpfile(basename,dist_seed0)
      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (16,file=filename(1:lnblnk(filename)),form='unformatted')
*         read (16) r
*         read (16) wsteps
      else
         open (16,file=filename(1:lnblnk(filename)))
*         read (16,20) r
*         read (16,30) wsteps
      endif

*      do k=1,wsteps
*         if (binio.eqv..true.) then
*            read (16) num,t
*         else
*            read (16,35) num,t
*         endif
         do i = 1,num_ctbins
            ct0bin(i) = 0
         enddo
*         do j=1,num
         do j=1,wlines_prcp
            if (binio.eqv..true.) then
               read (16) t,x,y,z,ke,alpha,ceta
            else
               read (16,40) t,x,y,z,ke,alpha,ceta
            endif
            ctheta=z/sqrt(x*x+y*y+z*z)
            do i = 1,num_ctbins
               lbound=ctmin+(i-1)*ctbwidth
               ubound=ctmin+i*ctbwidth
               if ((ctheta.ge.lbound).and.
     &         (ctheta.lt.ubound))ct0bin(i) = ct0bin(i)+1
            enddo
         enddo
         print *
*         print *,' ct dist. for wtime=',wtime(k)/tfactor
         do j = 1,num_ctbins
            lbound=ctmin+(j-1)*ctbwidth
            ubound=ctmin+j*ctbwidth
            print *,j,lbound,ubound,ct0bin(j)
         enddo
*      enddo

20    format (1e12.4)
30    format (1i8)
35    format (1i8,1e12.4)
40    format (7e12.4)

      return
      end

************************************************************************

      subroutine print_flux(basename,firstfilenum)

c     i/o finalize and post processing

      implicit none	
      include 'rbelt-io.inc'
      include 'rbelt-const.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-bounds.inc'
      integer i,j,k,firstfilenum,lnblnk,pnum,step,np,nts,num
      character*80 basename,distfile,fluxfile,filename,string_out
      real x,y,z,r,t,pvx,pvy,pvz,ke,alpha,eta,p_parl,p_perp,calpha

      integer num_rbins,num_cabins
      parameter(num_rbins=38,num_cabins=8)
      real camax,camin
      integer r0bin(num_rbins),ca0bin(num_cabins)
      real rbwidth,cabwidth,lbound,ubound

      print *
      print *,'*** in subroutine show_dist ***'

      camin=-1.0
      camax=1.0
      cabwidth=(camax-camin)/num_cabins
      do i = 1,num_cabins
         ca0bin(i) = 0
      enddo

*      rmin=0
*      rmax=1.0
      rbwidth=(rmax-rmin)/num_rbins
      do i = 1,num_rbins
         r0bin(i) = 0
      enddo

      filename=fluxfile(basename,dist_seed0)
      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (16,file=filename(1:lnblnk(filename)),form='unformatted')
         read (16) np
         read (16) nts
      else
         open (16,file=filename(1:lnblnk(filename)))
         read (16,20) np
         read (16,20) nts
      endif

      do i=1,nts
         if (binio.eqv..true.) then
            read (16) num,t
         else
            read (16,30) num,t
         endif
         do j = 1,num_cabins
            ca0bin(j) = 0
         enddo
         do j = 1,num_rbins
            r0bin(j) = 0
         enddo
         do k=1,num
            if (binio.eqv..true.) then
               read (16) pnum,x,y,alpha,eta
            else
               read (16,40) pnum,x,y,alpha,eta
            endif
            calpha=cos(pi*alpha/180)
            r=sqrt(x**2+y**2)
            do j = 1,num_cabins
               lbound=camin+(j-1)*cabwidth
               ubound=camin+j*cabwidth
               if ((calpha.ge.lbound).and.
     &         (calpha.lt.ubound))ca0bin(j) = ca0bin(j)+1
            enddo
            do j = 1,num_rbins
               lbound=rmin+(j-1)*rbwidth
               ubound=rmin+j*rbwidth
               if ((r.ge.lbound).and.
     &         (r.lt.ubound)) r0bin(j) = r0bin(j)+1
            enddo
         enddo
*         print *
*         print *,' ca flux. for wtime=',wtime(i)/tfactor
*         do j = 1,num_cabins
*            lbound=camin+(j-1)*cabwidth
*            ubound=camin+j*cabwidth
*            print *,j,lbound,ubound,ca0bin(j)
*         enddo
         print *
         print *,' r flux. for wtime=',wtime(i)/tfactor
         do j = 1,num_rbins
            lbound=rmin+(j-1)*rbwidth
            ubound=rmin+j*rbwidth
            print *,j,lbound,ubound,r0bin(j),
     &      r0bin(j)/(pi*(ubound**2-lbound**2))
         enddo
      enddo

20    format (i8)
30    format (1i8,1e12.4)
40    format (1i8,4e12.4)

      return
      end

************************************************************************

      subroutine field_io(basename,filenum,firstfilenum,timesteps)      
      implicit none
      include 'rbelt-io.inc'
      integer filenum,firstfilenum,timesteps,lnblnk
      character*80 basename,filename,lfile,bfile

*      if (bfld_out.eqv..true.)
*     &call bwrite(basename,filenum,firstfilenum,timesteps)
*      if (lcalc.ne.0) 
*     &call lstar_out(basename,filenum,firstfilenum,timesteps)
*      call lstar_out(basename,filenum,firstfilenum,timesteps)
c     use eqmap =1,2,3... to include other options
      if (eqmap.ge.1) call geo2smeq(basename,filenum,firstfilenum)

      return
      end

***********************************************************************

      subroutine bwrite(basename,filenum,firstfilenum,timesteps)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-io.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-ut.inc'
      include 'rbelt-fields.inc'

      integer i,j,k,l,l1,errcode,filenum,firstfilenum,timesteps
      character*80 basename
      real r,t,bd,ed,bxavg,byavg,bzavg,exavg,eyavg,ezavg,tavg

      integer i1,i2,j1,j2
*      parameter (i1=70,i2=nx,j1=40,j2=141,k=(nz+1)/2)
      parameter (i1=1,i2=1,j1=1,j2=1,k=(nz+1)/2)

      print *
      print *,'*** in subroutine b_out ***'

      if (filenum.eq.firstfilenum) then
         write (28,4)year0,doy0,hour0,min0,sec0
         write (28,6)(i2-i1+1),(j2-j1+1),timesteps
         write (28,8) xgr(i1),xgr(i2),ygr(j1),ygr(j2)
4        format (5i8)
6        format (3i8)
8        format (4e12.4)
      endif

      bd=1./ffactor/ntg
      ed=1./ffactor/vmsvcm

      l1=3
      if (filenum.eq.firstfilenum) l1=1
      do l=l1,nt
         write (28,10)(tgr(l)+tzero)/tfactor
         do j=j1,j2
            do i=i1,i2
               r=sqrt(xgr(i)*xgr(i)+ygr(j)*ygr(j))
               if (r.ge.rmin) then

                  bx=bxd(i,j,k,l)*bd
                  by=byd(i,j,k,l)*bd
                  bz=(bzd(i,j,k,l)+b0/r**3)*bd
                  ex=exd(i,j,k,l)*ed
                  ey=eyd(i,j,k,l)*ed
                  ez=ezd(i,j,k,l)*ed

                  write (28,12)bx,by,bz,ex,ey,ez
               else
                  write (28,12)0.,0.,0.,0.,0.,0.
               endif
            enddo
         enddo
      enddo

*      d=1./(nt-2)/ffactor/ntg
*      b0nt=b0/ffactor/ntg
*      tavg=0.
*      do l=3,nt
*         tavg=tavg+tgr(l)
*      enddo
*      tavg=tavg/(nt-2)
*      write (28,10)(tavg+tzero)/tfactor
*      do j=j1,j2
*         do i=i1,i2
*            r=sqrt(xgr(i)*xgr(i)+ygr(j)*ygr(j))
*            if (r.ge.rmin) then
*               bxavg=0.
*               byavg=0.
*               bzavg=0.
*               exavg=0.
*               eyavg=0.
*               ezavg=0.
**               do l=3,nt
*               do l=3,nt
*                  bxavg=bxavg+bxd(i,j,k,l)
*                  byavg=byavg+byd(i,j,k,l)
*                  bzavg=bzavg+bzd(i,j,k,l)
*                  exavg=exavg+exd(i,j,k,l)
*                  eyavg=eyavg+eyd(i,j,k,l)
*                  ezavg=ezavg+ezd(i,j,k,l)
*               enddo
*               write (28,12)bxavg*d,byavg*d,(bzavg+b0/r**3)*d,exavg*d,
*     &         eyavg*d,ezavg*d
*            else
*               write (28,12)0.,0.,0.,0.,0.,0.
*            endif
*         enddo
*      enddo

10    format (e16.4)
12    format (7e12.4,i8)

      return
      end



***********************************************************************

      subroutine lstar_out(basename,filenum,firstfilenum,timesteps)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-io.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-ut.inc'
      include 'rbelt-fields.inc'

      integer i,j,k,l,m,l1,errcode,filenum,firstfilenum,timesteps,num
      character*80 basename
      real x(3),r,t,pa,Lm,lstar,leI0,bloc,bmir,bmin
      logical lmonly

c     na needs to be positive, even, & non-zero
      integer na
c     number of PA grid points should be even
      parameter (na=28)
      real ca(na),camin,camax,da,r2,maxr2
      parameter (camin=-1.+1./na,camax=1.-1./na)

      integer i1,i2,j1,j2,m1,m2,nxout,nyout
*      parameter (i1=70,i2=nx,j1=40,j2=141,m1=15,m2=15)
      parameter (i1=1,i2=1,j1=1,j2=1,m1=1,m2=1)
      parameter (nxout=(i2-i1+1),nyout=(j2-j1+1))
c     (nx=170)

      integer lflag(nxout,nyout)
      real ls2d(nxout,nyout)

      print *
      print *,'*** in subroutine lstar_out ***'

      if (filenum.eq.firstfilenum) then
         write (30,4)year0,doy0,hour0,min0,sec0
         write (30,6)(i2-i1+1),(j2-j1+1),(m2-m1+1),timesteps
         write (30,8) xgr(i1),xgr(i2),ygr(j1),ygr(j2),ca(m1),ca(m2)
4        format (5i8)
6        format (4i8)
8        format (6e12.4)
      endif

      da=(camax-camin)/(na-1)
      do i=1,na
         ca(i)=camin + (i-1)*da
*         print *,'i,camin,ca(i),camax,PA=',i,ca(i)-da*.5,ca(i),
*     &   ca(i)+da*.5,acos(ca(i))*raddeg
      enddo

*      stop

      lmonly=.false.
      if (lcalc.eq.1) lmonly=.true.

      x(3)=0.
      l1=3
      if (filenum.eq.firstfilenum) l1=1

      do l=l1,nt
        t=tgr(l)
        num=0
        do j=1,nyout
          x(2)=ygr(j1+(j-1))
*          print *,'y=',x(2)
          do i=1,nxout
            x(1)=xgr(i1+(i-1))
            lflag(i,j)=1

            r=sqrt(x(1)*x(1)+x(2)*x(2))
*            if (r.ge.rmin) then
            if (r.ge.3.5) then
*              do m=m1,m2
*                pa=acos(ca(m))
                pa=pi/2.

                call rbelt_lshell
     &          (x,t,pa,lmonly,Lm,lstar,leI0,bloc,bmir,bmin,errcode)

                if (errcode.eq.0) then
                   ls2d(i,j)=lstar
                   lflag(i,j)=0
                   num=num+1
                endif

*              enddo
            endif

          enddo
        enddo

        write (30,10)(tgr(l)+tzero)/tfactor,num
        do j=1,nyout
          do i=1,nxout
             if (lflag(i,j).eq.0) write(30,14)i,j,ls2d(i,j)
          enddo
        enddo

      enddo

10    format (e12.4,i8)
14    format (2i4,f6.2)

      return
      end

***********************************************************************

      subroutine geo2smeq(basename,filenum,firstfilenum)

      implicit none
      include 'rbelt-io.inc'
      include 'rbelt-status.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-ut.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'

      integer filenum,firstfilenum
      character*80 basename

      real rad_geo,pol_geo,azi_geo13,azi_geo15
c     goes-13 at 75 deg. west (azi_geo = 285 deg.)
c     goes-15 at 135 deg. west (azi_geo = 225 deg.)
      parameter (rad_geo=6.6,pol_geo=90./raddeg)
      parameter (azi_geo13=285./raddeg,azi_geo15=225./raddeg)

      integer ifail13,ifail15
      real XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,XSM,YSM,ZSM,rad_sm,pol_sm,
     &azi_sm,t,t0,x3sm(3),xsm13,ysm13,bloc13,Beqtr13,
     &xsm15,ysm15,bloc15,Beqtr15

      print *
      print *,'*** in subroutine geo2smeq ***'

      if (filenum.eq.firstfilenum) then
c        write header
         write (32,4)year0,doy0,hour0,min0,sec0
         write (32,6)num_wsteps
4        format (5i8)
6        format (i8)
c        initialize write-step pointer
         wspntr=1
      endif

200   continue
c     t is normalized time w.r.t. start of curent time grid
      t=(wtime(wspntr)-tzero)
      if ((t.ge.tgr(1)).and.(t.lt.tgrmax)) then
         print *,'wspntr,wtime(wspntr),t=',
     &   wspntr,wtime(wspntr)/tfactor,t/tfactor
c        t0 is un-normalized time (in seconds) w.r.t. UT ref. time
         t0=wtime(wspntr)/tfactor
         call rbelt_recalc(t0)

c        GEO to SM coordinates
         call sphr2cart(rad_geo,pol_geo,azi_geo13,XGEO,YGEO,ZGEO)
         call GEOGSW_08 (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1)
         call SMGSW_08 (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,-1)
         call cart2sphr(rad_sm,pol_sm,azi_sm,XSM,YSM,ZSM)
c        get B at spacecraft location
         x3sm(1)=XSM
         x3sm(2)=YSM
         x3sm(3)=ZSM
         status=0
         call get_fields(x3sm,t)
         bloc13=abs(b/ffactor/ntg)
c        trace to SM coordinates equatorial plane
         call trace2eqtr(t,x3sm,Beqtr13,ifail13)
         xsm13=x3sm(1)
         ysm13=x3sm(2)

c        GEO to SM coordinates
         call sphr2cart(rad_geo,pol_geo,azi_geo15,XGEO,YGEO,ZGEO)
         call GEOGSW_08 (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1)
         call SMGSW_08 (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,-1)
         call cart2sphr(rad_sm,pol_sm,azi_sm,XSM,YSM,ZSM)
c        get B at spacecraft location
         x3sm(1)=XSM
         x3sm(2)=YSM
         x3sm(3)=ZSM
         status=0
         call get_fields(x3sm,t)
         bloc15=abs(b/ffactor/ntg)
c        trace to SM coordinates equatorial plane
         call trace2eqtr(t,x3sm,Beqtr15,ifail15)
         xsm15=x3sm(1)
         ysm15=x3sm(2)

         print *,'r,r_proj,r_map=',sqrt(XSM*XSM+YSM*YSM+ZSM*ZSM),
     &   sqrt(XSM*XSM+YSM*YSM),sqrt(XSM15*XSM15+YSM15*YSM15)

         write(32,12)t0,xsm13,ysm13,bloc13,Beqtr13,ifail13,
     &   xsm15,ysm15,bloc15,Beqtr15,ifail15

12       format (e14.6,4e12.4,i8,4e12.4,i8)
         wspntr=wspntr+1
         goto 200
      endif

         stop

      return
      end

************************************************************************

      subroutine openfile()
      implicit none
      return
      end

************************************************************************

      subroutine closefile()
      implicit none
      return
      end

***********************************************************************

      subroutine sm_lsm(x,y,z,vx,vy,vz,vx_LSM,vy_LSM,vz_LSM)
*      subroutine sm_lsm(x,y,z,vx,vy,vz,asink,eta,zenith)
*      call sm_lsm(y(1),y(2),y(3),-y(4),-y(5),-y(6),asink,eta,zenith)

      implicit none
      include 'rbelt-const.inc'
      real phi,theta,cphi,sphi,stheta,ctheta
      real vx_SM,vy_SM,vz_SM,vx_LSM,vy_LSM,vz_LSM
      real x,y,z,vx,vy,vz
      real zenith,k,eta,asink

c rotates from solar magnetic (SM) to local surface magnetic (LSM)
c coordinates using 3 eulerian angle rotations: 
c 1) phi = azimuthal angle (in SM) - pi/2.0 CCW about z_SM
c 2) theta =  pi -  polar angle (in SM) CCW about x'
c 3) psi = pi/2 CCW about z'' (same as z_LSM)
c then outputs the zenith, k=sin(theta_S) & eta angles 
c for cone projection.

*      pi = 3.14159265358979
*      raddeg = 180./pi

      vx_SM = vx
      vy_SM = vy
      vz_SM = vz

      phi = atan2(y,x) - pi/2.0
      cphi = cos(phi)
      sphi = sin(phi)
      theta = pi -  atan2(sqrt(x*x + y*y),z)
      ctheta = cos(theta)
      stheta = sin(theta)
      vx_LSM = -ctheta*sphi*vx_SM + ctheta*cphi*vy_SM + stheta*vz_SM
      vy_LSM = -cphi*vx_SM - sphi*vy_SM
      vz_LSM = stheta*sphi*vx_SM - stheta*cphi*vy_SM + ctheta*vz_SM

      zenith=(pi-atan2(
     &sqrt(vx_LSM*vx_LSM + vy_LSM*vy_LSM),vz_LSM))*raddeg
      k=vy_LSM/
     &(sqrt(vx_LSM*vx_LSM + vy_LSM*vy_LSM + vz_LSM*vz_LSM)+1.e-9)
      if ((abs(vy_LSM)).lt.1.e-12) k=0.0 !WARNING, POTENTIAL PROBLEMS HERE!

      asink = asin(k)
      eta=(pi-atan2(abs(vx_LSM),vz_LSM))

*      vx = vx_LSM
*      vy = vy_LSM
*      vz = vz_LSM

*      print *
*      print *,'vx_SM,vy_SM,vz_SM=',vx_SM,vy_SM,vz_SM
*      print *,'vx_LSM,vy_LSM,vz_LSM=',vx_LSM,vy_LSM,vz_LSM
*      print *,'new magnitude=',sqrt(vx_LSM**2+vy_LSM**2+vz_LSM**2)
*      print *,'zenith=',zenith*raddeg
*      print *,'k=',k
*      print *,'asin(k)=',asin(k)
*      print *,'eta=',acos(cos(zenith)/cos(asin(k)))*raddeg

      end

************************************************************************
*
*      subroutine midnightl(i,y,dydt,t,t2,dt,thalt3)
*
*c     midchckdt is how often we check to see if we're in midnight sector
*c     midchckdt = 5.0,
*c     phidelta = 10.0,
*
*      implicit none	
*      include 'rbelt-io.inc'
*      include 'rbelt-const.inc'
*      integer i
*      real y(6),dydt,t,t2,dt,thalt3,phi,x0,y0,z0,t0,l_eqtr,t_eqtr
*      real dxstep,dystep,dzstep,ds
*      real wghtavg,phi360
*
*      phi=phi360(y(2),y(1))
**      print *,'t,phi,z,y=',t/tfactor,phi*raddeg,y(3),y(2)
*c     check if we are in midnight sector
*      if ((abs(phi-pi).lt.phidelta).and.(stepflag.ne.2)) then
*         if (initflag.eq.1) then
*c           look for equatorial crossing
*            if (y(3)*z_prev.lt.0.0) then
*c              increment equatorial crossing counter
**               neqtrcross=neqtrcross+1 
*c              record equatorial crossings until we find a midnight crossing
*               if (stepflag.eq.0) then
*c                 linear interpolation to get location of equatorial crossing
*                  x0=wghtavg(y(1),x_prev,y(3),0.,z_prev)
*                  y0=wghtavg(y(2),y_prev,y(3),0.,z_prev)
*                  t0=wghtavg(t,t_prev,y(3),0.,z_prev)
*                  l_eqtr_prev=sqrt(x0*x0+y0*y0)
*		  t_eqtr_prev=t0
*c              if we have recorded equatorial and midnight crossings
*c              then calculate midnight crossing L
*               elseif (stepflag.eq.1) then
*c                 linear interpolation to get location of equatorial crossing
*                  x0=wghtavg(y(1),x_prev,y(3),0.,z_prev)
*                  y0=wghtavg(y(2),y_prev,y(3),0.,z_prev)
*                  t_eqtr=wghtavg(t,t_prev,y(3),0.,z_prev)
*                  l_eqtr=sqrt(x0*x0+y0*y0)
*c                 linear interpolation to get approximate equatorial footpoint 
*c                 of midnight crossing (constant azimuthal drift speed assumed)
*                  lm=wghtavg(l_eqtr,l_eqtr_prev,t,t_mid,t_eqtr_prev)
*                  stepflag=2
**                  neqtrcross=0
*                  print *,'l_eqtr,lm=',l_eqtr,lm
*               endif
*            endif
*c           look for midnight crossing
*            if (y(2)*y_prev.lt.0.0) then
*c              linear interpolation to get time of midnight crossing
*               t_mid=wghtavg(t,t_prev,y(2),0.,y_prev)
*               stepflag=1
*            endif
*         endif
*         x_prev=y(1)
*         y_prev=y(2)
*         z_prev=y(3)
*	 t_prev=t
*         thalt3 = t
*         initflag = 1
*      else
*         if ((initflag.eq.1).and.(abs(phi).lt.phidelta)) then
*            if (stepflag.ne.2) then
*               print *,'incomplete Lm search:'
*               print *,'need smaller midchckdt or larger phidelta'
*               stop
*            else
*               print *,'noon crossing'
*               initflag=0
*               stepflag=0
*            endif
*         endif
*         thalt3 = t + midchckdt
*      endif
*
*      return
*      end
*
************************************************************************
*c           calculate 2nd invarient       
*	    if (stepflag.eq.2) then
*c              start calculating 2nd invarient
*               if (neqtrcross.eq.0) then
*                  dxstep=y(1)-x0
*                  dystep=y(2)-y0
*                  dzstep=y(3)
*                  ds = sqrt(dxstep*dxstep+dystep*dystep+dzstep*dzstep)
*                  jbounce_sum=0
*                  jbounce_sum=jbounce_sum+abs(y(4)*ds)
*                  neqtrcross=1
*               endif
*c              calculate 2nd invarient
*               if (neqtrcross.ge.1 .and. neqtrcross.le.2) then
*                  dxstep=y(1)-x_prev
*                  dystep=y(2)-y_prev
*                  dzstep=y(3)-z_prev
*                  ds = sqrt(dxstep*dxstep+dystep*dystep+dzstep*dzstep)
*                  jbounce_sum = jbounce_sum + abs(y(4)*ds)
*               endif
*c              finish calculating 2nd invarient
*               if (neqtrcross.eq.3) then
*                  x0=wghtavg(y(1),x_prev,y(3),0.,z_prev)
*                  y0=wghtavg(y(2),y_prev,y(3),0.,z_prev)
*                  dxstep=y(1)-x0
*                  dystep=y(2)-y0
*                  dzstep=y(3)
*                  ds = sqrt(dxstep*dxstep+dystep*dystep+dzstep*dzstep)
*                  jbounce_sum = jbounce_sum + abs(y(4)*ds)
*                  stepflag=3
*               endif
*	    endif
*         endif
*         x_prev=y(1)
*         y_prev=y(2)
*         z_prev=y(3)
*	 t_prev=t
*         thalt3 = t
*         initflag = 1
*      else
*         if (initflag.eq.1) then
*            if (stepflag.ne.3) then
*               print *,'incomplete 2nd invar. calc.'
*               print *,'may need smaller midchckdt or larger phidelta'
*               stop
*            else
*               print *,'jbounce_sum=',jbounce_sum
*               initflag=0
*               stepflag=0
*            endif
*         endif
*         thalt3 = t + midchckdt
*      endif
*      return
*      end

c      jbounce_sum=jbounce_sum+y(4)*ds
*      jbounce_sum=jbounce_sum+y(4)*(y(4)*dt)
*      if (flag.eq.0) return
*      sum=sum+v_par*ds
*      if (y(4)*y4_previous.le.0..and.y(4).ne.0.) then
*         jbounce = factor * sum + approx. ds term
*         sum=0        
*      endif
*      if (y(3)*z_prev.lt.0.0) then
c          record EPA
*      endif
*      z_prev= y(3)

***************************************

*         print *,'t,z,r=',t/tfactor,real(y(3)),
*     &   sqrt(y(1)**2+y(2)**2+y(3)**2)
c        mirror point check
c        record b_eq & reset b_min at the mirror point
*         b_min=amin1(b,b_min)
*         if (y(4)*p_prev.lt.0.0) then
*            b_eq=b_min
*            print *,'b min=',b_eq/ffactor/ntg
*            b_min=abs(2.0*b0)
*         endif
*         p_prev= y(4)
c        equitorial crossing check
c        get alpha_eq & r_eq and record flux count
*         if (y(3)*z_prev.lt.0.0) then
*           print *,'b at eq.=',b/ffactor/ntg
*           alpha_eq = pangle_(y,t)
*           r_eq = sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
c          for more accuracy this should be modified to get the
c          quantities right at the eqitorial crossing
*           call bounce_drift(y,t)
*           call data2yout(i,y,t,dt,1)
*         endif
c        use spb with bounce_drift routine only
*         spb=spb+1
*         z_prev= y(3)
c        if dthalt time interval has passed, do halt2 stuff
*         if (t.ge.thalt2) then
*            call print_gcforce(i,y,t,dt)
*            call invar_chck(y,t)
*            thalt2=tstop+dthalt
*         endif

************************************************************************

      real function phi360(y,x)

      implicit none
      include 'rbelt-const.inc'
      real y,x
      phi360=atan2(y,x)
      if (phi360.lt.0.0) phi360=2*pi+phi360
      return
      end

************************************************************************

      real function wghtavg(y,y_prev,x,x0,x_prev)

      implicit none

      real x,x0,x_prev,y,y_prev,wght,wght_prev
 	  
      wght = (x0-x_prev)/(x-x_prev)
      wght_prev = (x-x0)/(x-x_prev)
      wghtavg = wght*y + wght_prev*y_prev

      return
      end

************************************************************************

      SUBROUTINE rk4save4(y,dydx,h,rhs,rhs_init) 

c     4th order Runge-Kutta method. 
c     save 4 RK4 step values of y(4) to use for dJ/dt to calc 2nd invar.
c     compare with pure dipole results to test.
c     (approximate along GC trajectory)
c     what about other ways to calculate??

      implicit none
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'
      external rhs,rhs_init
      integer i,nvars
      real x,xh
      real h,h6,dydx(4),y
      h6=h/6. 

      y=y+h6*(dydx(1)+2.*(dydx(2)+dydx(3))+dydx(4)) 

      return
      end





