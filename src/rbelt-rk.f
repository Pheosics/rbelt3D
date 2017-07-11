c     rbelt-rk.f - various Runge-Kutta methods

************************************************************************

      SUBROUTINE rk4(y,dydx,x,h,nvars,rhs,rhs_init) 

c     4th order Runge-Kutta method. advances particle position and time,
c     and calculates next dt & dydt providing status <= 0, otherwise 
c     the input variables are returned unchanged.
c     assumes initial dydx has been calc. in advance

      implicit none
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'
      external rhs,rhs_init
      integer i,nvars
      real x,xh
      real h,dydx(nvars),y(nvars)
      real h6,hh,dym(nvars),dyt(nvars),yt(nvars) 

      hh=h*0.5 
      h6=h/6. 
      xh=x+hh

*      print *
*      print *,'*** in SUBROUTINE rk4 ***'
*      print *,' y,t,h,dydx=',y,x/tfactor,h/tfactor,dydx

      do i=1,nvars
         yt(i)=y(i)+hh*dydx(i)
      enddo
      
      call rhs(xh,yt,dyt)
      do i=1,nvars 
         yt(i)=y(i)+hh*dyt(i) 
      enddo

      call rhs(xh,yt,dym)
      do i=1,nvars
         yt(i)=y(i)+h*dym(i)
         dym(i)=dyt(i)+dym(i) 
      enddo

      xh=x+h
      call rhs(xh,yt,dyt)
      do i=1,nvars
         yt(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i)) 
      enddo

c     need to check status at t=xh & y=yt before we advance.
c     also calculate dydx, h and B-field for next step
c     if status <= 0 we take the step
c     if status > 0 we do not
      call rhs_init(xh,yt,dydx,h)
      if (status.le.0) then
         do i=1,nvars
            y(i)=yt(i)
         enddo
         x=xh
      else 
c        resets dydy, h and B-field to begining of step
         call rhs_init(x,y,dydx,h)
      endif

      return
      end

************************************************************************

      subroutine rkf4(y,dydt1,t,dt,nvars,rhs,rhs_init)

c     4th order Runge-Kutta-Fieldberg method. advances particle position
c     and time, and calculates next dt & dydt providing status <= 0, otherwise 
c     the input variables are returned unchanged.

      implicit none
      include 'rbelt-status.inc'
*      include 'rbelt-const.inc'
      external rhs,rhs_init
      integer i,nvars
      real y(nvars),dydt1(nvars),dydt2(nvars),dydt3(nvars),t,dt
      real dydt4(nvars),dydt5(nvars),dydt6(nvars)
      real ytmp(nvars)

      real a2,a3(2),a4(3),a5(4),a6(5),e(6),s(6),dtw(6)
      common /rkcoef/ a2,a3,a4,a5,a6,e,s,dtw

*      print *,'in rkf4: t,y=',t,y

      do i=1,nvars
         ytmp(i)=y(i)+dt*.2*dydt1(i)
      enddo

      call rhs(t+.25*dt,ytmp,dydt2)
      do i=1,nvars
         ytmp(i)=y(i)+dt*(.075*dydt1(i)+.225*dydt2(i))
      enddo

      call rhs(t+.375*dt,ytmp,dydt3)
      do i=1,nvars
         ytmp(i)=y(i)+dt*(.3*dydt1(i)-.9*dydt2(i)+1.2*dydt3(i))
      enddo

      call rhs(t+dtw(4)*dt,ytmp,dydt4)
      do i=1,nvars
         ytmp(i)=y(i)+dt*(a5(1)*dydt1(i)+2.5*dydt2(i)
     &   +a5(3)*dydt3(i)+a5(4)*dydt4(i))
      enddo

      call rhs(t+dt,ytmp,dydt5)
      do i=1,nvars
         ytmp(i)=y(i)+dt*(a6(1)*dydt1(i)+a6(2)*dydt2(i)
     &   +a6(3)*dydt3(i)+a6(4)*dydt4(i)+a6(5)*dydt5(i))
      enddo

      call rhs(t+.5*dt,ytmp,dydt6)
      do i=1,nvars
         ytmp(i)=y(i)+dt*(s(1)*dydt1(i)+s(3)*dydt3(i)
     &   +s(4)*dydt4(i)+s(6)*dydt6(i))
      enddo

      call rhs_init(t+dt,ytmp,dydt1,dt)
      if (status.le.0) then
         do i=1,nvars
            y(i)=ytmp(i)
         enddo
         t=t+dt
      else 
         call rhs_init(t,y,dydt1,dt)
      endif

      return
      end

************************************************************************

	subroutine rkckparm

	  implicit none
	  real a2,a3(2),a4(3),a5(4),a6(5),e(6),s(6),dtw(6)
	  common /rkcoef/ a2,a3,a4,a5,a6,e,s,dtw

	  a2=.2
	  a3(1)=3./40.
	  a3(2)=9./40.
	  a4(1)=.3
	  a4(2)=-.9
	  a4(3)=1.2
	  a5(1)=-11./54.
	  a5(2)=2.5
	  a5(3)=-70./27.
	  a5(4)=35./27.
	  a6(1)=1631./55296.
	  a6(2)=175./512.
	  a6(3)=575./13824.
	  a6(4)=44275./110592.
	  a6(5)=253./4096.
	  e(1)=277./64512.
	  e(2)=0.
	  e(3)=-6925./370944.
	  e(4)=6925./202752.
	  e(5)=277./14336.
	  e(6)=-277./7084.
c	  s(1)=2825./27648.		!4th-order solutions
c	  s(2)=0.
c	  s(3)=18575./48384.
c	  s(4)=13525./55296.
c	  s(5)=277./14336.
c	  s(6)=.25
	  s(1)=37./378.		!5th-order solutions
	  s(2)=0.
	  s(3)=250./621.
	  s(4)=125./594.
	  s(5)=0.
	  s(6)=512./1771.
	  dtw(1)=0.
	  dtw(2)=1./4.
	  dtw(3)=3./8.
	  dtw(4)=12./13.
	  dtw(5)=1.
	  dtw(6)=1./2.

      return
      end

************************************************************************


      subroutine rk2(y,dydt1,t,dt,nvars,rhs,rhs_init)

c     2nd order Runge-Kutta method. advances particle position and time,
c     and calculates next dt & dydt providing status <= 0, otherwise 
c     the input variables are returned unchanged.

      implicit none
      include 'rbelt-status.inc'
      external rhs,rhs_init
      integer i,nvars
      real t,dt
      real y(nvars),dydt1(nvars),dydt2(nvars)
      real ytmp(nvars)

      do i=1,nvars
         ytmp(i)=y(i)+dt*.5*dydt1(i)
      enddo

      call rhs(t+.5*dt,ytmp,dydt2)
      do i=1,nvars
         ytmp(i)=y(i)+dt*dydt2(i)
      enddo

      call rhs_init(t+dt,ytmp,dydt1,dt)
      if (status.le.0) then
         do i=1,nvars
            y(i)=ytmp(i)
         enddo
         t=t+dt
      else 
         call rhs_init(t,y,dydt1,dt)
      endif

      return
      end

************************************************************************

      SUBROUTINE rk4t0(y,dydx,x,x0,h,nvars,rhs,rhs_init) 

c     4th order Runge-Kutta method. advances particle position and time,
c     and calculates next dt & dydt providing status <= 0, otherwise 
c     the input variables are returned unchanged.
c     assumes initial dydx has been calc. in advance

      implicit none
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'
      external rhs,rhs_init
      integer i,nvars
      real x,xh,x0
      real h,dydx(nvars),y(nvars)
      real h6,hh,dym(nvars),dyt(nvars),yt(nvars) 

      hh=h*0.5 
      h6=h/6. 
      xh=x+hh

*      print *
*      print *,'*** in SUBROUTINE rk4 ***'
*      print *,' y,t,h,dydx=',y,x/tfactor,h/tfactor,dydx

      do i=1,nvars
         yt(i)=y(i)+hh*dydx(i)
      enddo
      
      call rhs(x0,yt,dyt)
      do i=1,nvars 
         yt(i)=y(i)+hh*dyt(i) 
      enddo

      call rhs(x0,yt,dym)
      do i=1,nvars
         yt(i)=y(i)+h*dym(i)
         dym(i)=dyt(i)+dym(i) 
      enddo

      xh=x+h
      call rhs(x0,yt,dyt)
      do i=1,nvars
         yt(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i)) 
      enddo

c     need to check status at t=xh & y=yt before we advance.
c     also calculate dydx, h and B-field for next step
c     if status <= 0 we take the step
c     if status > 0 we do not
      call rhs_init(x0,yt,dydx,h)
      if (status.le.0) then
         do i=1,nvars
            y(i)=yt(i)
         enddo
         x=xh
      else 
c        resets dydy, h and B-field to begining of step
         call rhs_init(x0,y,dydx,h)
      endif

      return
      end
