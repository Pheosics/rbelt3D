c     rbelt-time_loop.f - loops over time & calls integrator 

************************************************************************

      subroutine time_loop(i,y,t,t2,integrator,nvars,rhs,rhs_init)

c     advances single particle in position & time while status=0.
c     if t >= tstop we halt for i/o or to check if we can go to a 
c     guiding center trajectory, then continue with integration, 
c     if status is not 0 we return.
c
c            save this reminder:
c            if ((t+dt).eq.t) then 
c               print *,'(t+dt).eq.t'
c               stop ! need dbl precision time
c            endif
c

      implicit none
      include 'rbelt-bounds.inc'
      include 'rbelt-status.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-const.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-io.inc'
*      include 'rbelt-grid.inc'
      include 'rbelt-fields.inc'
      external integrator,rhs,rhs_init
      integer i,nvars,step,tmpstatus
      real*8 t,t2,dt,y(6),dydt(6),thalt1,thalt2,thalt3,
     &tstop,r,phi

*      print *
*      print *,'*** in subroutine time loop ***'
*      print *,'i,t,t2,nvars=',i,t/tfactor,t2/tfactor,nvars
*      stop

c     use tgrmax in linterp to stop particle at t2
      tgrmax=t2

c     set upper bound on tstop
      tstop=t2
      thalt1=t2
      thalt2=t2
      thalt3=t2

c     initialize i/o & additional integration stop times
      call yout_init(i,y,t,t2,thalt2,thalt3)
!      if ((flag.eq.0).and.(flag_switch.eqv..true.)) thalt1=t

c     set tstop (why is this here?)
      tstop=dmin1(thalt2,thalt3)

      if (flag.eq.0) call boris_init(t,y,dydt,dt)

c     time stepping loop
      step=0
      status=0
c     (nested do-while loops) 
200   if (status.eq.0) then

c        temporarily halt integration to...

c        check to see if we can switch to a guiding center traj.
!         if (t.ge.thalt1) call lrntz2gc_chck(y,t,thalt1)

c        output particle data
         if (t.ge.thalt2) then
            tmpstatus=status
            call io(i,y,dydt,t,t2,dt,thalt2)
            status=tmpstatus
         endif

c        reset tstop
         tstop=dmin1(thalt2,thalt3)

c        initialize integrator
         call rhs_init(t,y,dydt,dt)

c        integration loop
201      if (((t.le.tstop).or.(thalt2.ge.thalt3)).and.
     &(status.eq.0)) then
c           take a step
            call integrator(y,dydt,t,dt,nvars,rhs,rhs_init)
!            step=step+1
!            if (step.eq.10000000) then
!               print *, t/tfactor, dt/tfactor, t2/tfactor
!               step = 0
!            endif
            goto 201
         endif

         goto 200
      endif

      if (flag.eq.0) call boris_stop(t,y,dydt,dt)

      return
      end
