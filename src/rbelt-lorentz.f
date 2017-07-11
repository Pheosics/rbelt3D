c     rbelt-lorentz.f - Lorentz trajectory integration routines

************************************************************************

      subroutine lrntz_params()

c     read parameters for RK-lrntz integrator

      implicit none
      include 'rbelt-lorentz.inc'
      include 'rbelt-const.inc'

c     read in and normalize input parameters
      NAMELIST /lrntz/ tstep_lrntz,dx_max
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,lrntz)
      CLOSE (81)
c     normalize
      tstep_lrntz = tstep_lrntz*2*PI

      return
      end

************************************************************************

      subroutine rhs_lrntz(t,y,dydt)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      real t,y(6),dydt(6),gamma

*      print *
*      print *,'  in rhs: y,t=',y,t/tfactor

c     get fields at (y,t)
      call get_fields(y,t)

c     calculate gamma at time t
c     (might be better to also advance energy in some situations)
c     (might be better to use real*8 for gamma)
      gamma = sqrt(1.+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))

      dydt(1)=y(4)/gamma
      dydt(2)=y(5)/gamma
      dydt(3)=y(6)/gamma
      dydt(4)=ex + dydt(2)*bz - dydt(3)*by
      dydt(5)=ey + dydt(3)*bx - dydt(1)*bz
      dydt(6)=ez + dydt(1)*by - dydt(2)*bx

      return
      end

***********************************************************************

      subroutine dt_lrntz(y,dt)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-lorentz.inc'
      real y(6),dt,u2,gamma,v

      u2=y(4)*y(4)+y(5)*y(5)+y(6)*y(6)
      gamma = sqrt(1+u2)

      v=sqrt(u2)/gamma
      dt=tstep_lrntz*gamma/b
      dt=amin1(dt,dx_max/v)

      return
      end

************************************************************************

      subroutine rhs_init_lrntz(t,y,dydt,dt)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-lorentz.inc'
      real t,dt,y(6),dydt(6),gamma,u2,v

*      print *
*      print *,'  in rhs: y,t=',y,t/tfactor

c     get fields at (y,t)
      call get_fields(y,t)
      
c     calculate gamma at time t 
      u2 = y(4)*y(4)+y(5)*y(5)+y(6)*y(6)
      gamma = sqrt(1+u2)

      dydt(1)=y(4)/gamma
      dydt(2)=y(5)/gamma
      dydt(3)=y(6)/gamma
      dydt(4)=ex + dydt(2)*bz - dydt(3)*by
      dydt(5)=ey + dydt(3)*bx - dydt(1)*bz
      dydt(6)=ez + dydt(1)*by - dydt(2)*bx

      dt=tstep_lrntz*gamma/b
      v=sqrt(u2)/gamma
      dt=amin1(dt,dx_max/v)
*      print *,'dx=',v*dt

      return
      end









