c     rbelt-boris.f - Boris (leap-frog) time stepping

***********************************************************************
**     Subroutine boris_init
**	This routine starts the Boris algorithm.  B.Kress Feb. 2003
**

***********************************************************************

      subroutine dummy_init(t,y,dydt,dtio)
      return
      end

***********************************************************************

      subroutine boris_init(t,y,dydt,dtio)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      include 'rbelt-lorentz.inc'
      integer i
      real t
      real gamma,gamma2,t2,dt2,dt3,bv
      real ux,uy,uz,y(6),dt,dtio,dydt

*      print *,'*** in subroutine boris_init ***'

c     advance proper velocity dt using values for e,b,& gamma at (y,t)

c     update the fields to (y,t)
      call get_fields(y,t)
      if (status.gt.0) return

      gamma2= 1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6)
      gamma = sqrt(gamma2)

      dtio=tstep_lrntz*gamma/b
      dt=-0.5*dtio

      t2 = dt*dt*b*b/4/gamma2
      dt2 = dt/2/gamma
      dt3 = dt*dt/4/gamma2
      bv = bx*y(4)+by*y(5)+bz*y(6)

      ux = y(4)
      uy = y(5)
      uz = y(6)

*      print *
*      print *,'advance proper velocity dt(scaled)=',dt

      y(4)=ux+dt*ex+2*(dt2*(uy*bz-uz*by)+dt3*bv*bx-t2*ux)/(1+t2)
      y(5)=uy+dt*ey+2*(dt2*(uz*bx-ux*bz)+dt3*bv*by-t2*uy)/(1+t2)
      y(6)=uz+dt*ez+2*(dt2*(ux*by-uy*bx)+dt3*bv*bz-t2*uz)/(1+t2)

*      print *,' ux=',y(4)
*      print *,' uy=',y(5)
*      print *,' uz=',y(6)

      return
      end

***********************************************************************
**     Subroutine boris
**	  This routine pushes charged particles in E & B fields using the 
**	  Lorentz force. A Boris type (leapfrog) method is used. It must be 
**	  initiated with the boris_init routine. B.Kress Feb. 2003
**	  note: there is about 1% error in the bounce period, which 
**	  could possibly be reduced by working in the B field ref. frame.

      subroutine boris(y,dydt,t,dt,nvars,rhs,rhs_init)
*      subroutine boris_step(y,t,dt)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      external rhs,rhs_init
      integer nvars
      real y(nvars),dydt(nvars),t,dt
      real gamma,gamma2,t2,dt2,dt3,bv,ux,uy,uz

*      print *,'*** in subroutine boris ***'
*      print *,'t,y=',t/tfactor,y

c     update the fields to (y,t)
      call get_fields(y,t)
      if (status.gt.0) return

********************************************************************
c	advance proper velocity and position with boris method --
c	see C.K. Birdsall & A.B. Langdon pg. 58,59,356,&357
c	We start with u at t - dt/2 and y at t

c add half the E field impulse to u (initialy at t - dt/2)

      ux = y(4) + dt*ex/2.
      uy = y(5) + dt*ey/2.
      uz = y(6) + dt*ez/2.

c calculate gamma and other repeated values at time t 

      gamma2= 1+ux*ux+uy*uy+uz*uz
      gamma = sqrt(gamma2)
      t2 = dt*dt*b*b/4/gamma2
      dt2 = dt/2/gamma
      dt3 = dt*dt/4/gamma2
      bv = bx*ux+by*uy+bz*uz

c add rotation and 2nd half of E field impulse to u

      y(4)=ux+dt*ex/2+2*(dt2*(uy*bz-uz*by)+dt3*bv*bx-t2*ux)/(1+t2)
      y(5)=uy+dt*ey/2+2*(dt2*(uz*bx-ux*bz)+dt3*bv*by-t2*uy)/(1+t2)
      y(6)=uz+dt*ez/2+2*(dt2*(ux*by-uy*bx)+dt3*bv*bz-t2*uz)/(1+t2)

c calculate gamma at t + dt/2

      gamma = sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))

c use v = u/gamma at t + dt/2 to advance position from t to t + dt

      y(1)=y(1)+dt*y(4)/gamma
      y(2)=y(2)+dt*y(5)/gamma
      y(3)=y(3)+dt*y(6)/gamma

      t=t+dt

*      print *,'t,y=',t/tfactor,y

      return
      end

***********************************************************************

      subroutine boris_stop(t,y,dydt,dtio)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      include 'rbelt-lorentz.inc'
      integer i
      real t
      real gamma,gamma2,t2,dt2,dt3,bv
      real ux,uy,uz,y(6),dt,dtio,dydt

*      print *,'*** in subroutine boris_init ***'

c     advance proper velocity dt using values for e,b,& gamma at (y,t)

c     update the fields to (y,t)
      call get_fields(y,t)
      if (status.gt.0) return

      gamma2= 1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6)
      gamma = sqrt(gamma2)

      dt=dtio/2.0

      t2 = dt*dt*b*b/4/gamma2
      dt2 = dt/2/gamma
      dt3 = dt*dt/4/gamma2
      bv = bx*y(4)+by*y(5)+bz*y(6)

      ux = y(4)
      uy = y(5)
      uz = y(6)

*      print *
*      print *,'advance proper velocity dt(scaled)=',dt

      y(4)=ux-dt*ex+2*(dt2*(uy*bz-uz*by)+dt3*bv*bx-t2*ux)/(1+t2)
      y(5)=uy-dt*ey+2*(dt2*(uz*bx-ux*bz)+dt3*bv*by-t2*uy)/(1+t2)
      y(6)=uz-dt*ez+2*(dt2*(ux*by-uy*bx)+dt3*bv*bz-t2*uz)/(1+t2)

*      print *,' ux=',y(4)
*      print *,' uy=',y(5)
*      print *,' uz=',y(6)

      return
      end
