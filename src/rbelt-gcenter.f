c     rbelt-gcenter.f -- guiding center time integration routines

************************************************************************

      subroutine gc_params()

c     read & scale parameters for guiding center integrator

      implicit none
      include 'rbelt-gcenter.inc'
      include 'rbelt-const.inc'

c     read in and normalize input parameters
      NAMELIST /gcenter/ tstep_gc,dx_max_gc,go2lrntz,go2gc,etaswitch,
     &dtgo2gc,seed
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,gcenter)
      CLOSE (81)
c     normalize
      go2lrntz=go2lrntz*go2lrntz
      etaswitch=etaswitch
      go2gc=go2gc*go2gc
      dtgo2gc=dtgo2gc*tfactor

      return
      end

************************************************************************
**3456789012345678901234567890123456789012345678901234567890123456789012
**       1         2         3         4         5         6         7
************************************************************************
**     Subroutine RHS-northrop			elkington 1997
**	  Determines if the drift approximation is valid, and returns
**	  the value of the right-hand side of the drift equations,
**	  vx, vy, and vz, for the given space and time coordinates.
**
**	  This subroutine uses the guiding center equations as developed
**	  by Northrop (1963), and may exhibit energy and invariant
**	  conservation problems.
**

      subroutine rhs_northrop(t,y,dydt)

      implicit none
      include 'rbelt-mu.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      real*8 t
      real*8 y(4),dydt(4)
      real*8 b2,gamma,gmbm,gmb2m
      real*8 pgmb,p2gmb4,ax,ay,az

*      print *
*      print *,'in rhs: y,t=',y,t/tfactor

c     update the fields to (y,t)
      call get_fields2(y,t)

c     calculate drift velocities
c     first calc some repeated values
      b2=1/(b*b)
      gamma=dsqrt(1.+2.*dabs(b)*mu+y(4)*y(4))
      gmbm=mu/gamma/b
      gmb2m=gmbm/b
      pgmb=y(4)*gmbm/mu
      p2gmb4=y(4)*pgmb*b2/b

c     get the cross terms for the curvature drift

      ax=bx*dbxdx + by*dbxdy + bz*dbxdz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*bx/b
      ay=bx*dbydx + by*dbydy + bz*dbydz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*by/b
      az=bx*dbzdx + by*dbzdy + bz*dbzdz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*bz/b

c    now calculate the drift velocities
      dydt(1)=b2*(ey*bz-ez*by)		!ExB drift
     &	+gmb2m*(by*dbdz-bz*dbdy)	!grad-B drift
     &	+p2gmb4*(by*az-bz*ay)		!curvature drift
     &	+bx*pgmb			!parallel velocity
      dydt(2)=b2*(ez*bx-ex*bz)		!ExB drift
     &	+gmb2m*(bz*dbdx-bx*dbdz)	!grad-B
     &	+p2gmb4*(bz*ax-bx*az)		!curvature drift
     &	+by*pgmb			!parallel velocity
      dydt(3)=b2*(ex*by-ey*bx)		!ExB drift
     &	+gmb2m*(bx*dbdy-by*dbdx)	!grad-B drift
     &	+p2gmb4*(bx*ay-by*ax)		!curvature drift
     &	+bz*pgmb			!parallel force contrib
      dydt(4)=(ex*bx+ey*by+ez*bz)/b	!parallel force
     &		-(bx*dbdx+by*dbdy+bz*dbdz)*gmbm

      return
      end

************************************************************************

************************************************************************
**3456789012345678901234567890123456789012345678901234567890123456789012
**       1         2         3         4         5         6         7
************************************************************************
**     Subroutine RHS-northrop			elkington 1997
**	  Determines if the drift approximation is valid, and returns
**	  the value of the right-hand side of the drift equations,
**	  vx, vy, and vz, for the given space and time coordinates.
**
**	  This subroutine uses the guiding center equations as developed
**	  by Northrop (1963), and may exhibit energy and invariant
**	  conservation problems.
**

      subroutine rhs_init_northrop(t,y,dydt,dt)

      implicit none
      include 'rbelt-mu.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-gcenter.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-status.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-dist.inc'
      real*8 t,dt
      real*8 y(4),dydt(4)
      real*8 b2,b3,gamma,gmbm,gmb2m
      real*8 p,v,pgmb,p2gmb4,ax,ay,az
      real*8 fm,fc,fg,fe,bxgb1,bxgb2,bxgb3

*      print *
*      print *,'in rhs_init: y,t=',y,t/tfactor

c     update the fields to (y,t)
      call get_fields2(y,t)

c     calculate drift velocities
c     first calc some repeated values
      gamma=dsqrt(1.+2.*dabs(b)*mu+y(4)*y(4))
      p=dsqrt(gamma*gamma-1)
      b2=1/(b*b)
      gmbm=mu/gamma/b
      gmb2m=gmbm/b
      pgmb=y(4)*gmbm/mu
c     y(4)*y(4)/gamma/b/b/b/b
      p2gmb4=y(4)*pgmb*b2/b

c     get the cross terms for the curvature drift
      ax=bx*dbxdx + by*dbxdy + bz*dbxdz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*bx/b
      ay=bx*dbydx + by*dbydy + bz*dbydz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*by/b
      az=bx*dbzdx + by*dbzdy + bz*dbzdz
*     &  -(bx*dbdx + by*dbdy + bz*dbdz)*bz/b

      bxgb1=by*dbdz-bz*dbdy
      bxgb2=bz*dbdx-bx*dbdz
      bxgb3=bx*dbdy-by*dbdx

c    now calculate the drift velocities
      dydt(1)=b2*(ey*bz-ez*by)		!ExB drift
     &	+gmb2m*bxgb1	!grad-B drift
     &	+p2gmb4*(by*az-bz*ay)			!curvature drift
     &	+bx*pgmb			!parallel velocity
      dydt(2)=b2*(ez*bx-ex*bz)		!ExB drift
     &	+gmb2m*bxgb2	!grad-B
     &	+p2gmb4*(bz*ax-bx*az)			!curvature drift
     &	+by*pgmb			!parallel velocity
      dydt(3)=b2*(ex*by-ey*bx)		!ExB drift
     &	+gmb2m*bxgb3	!grad-B drift
     &	+p2gmb4*(bx*ay-by*ax)			!curvature drift
     &	+bz*pgmb			!parallel force contrib
      dydt(4)=(ex*bx+ey*by+ez*bz)/b	!parallel force
     &	-(bx*dbdx+by*dbdy+bz*dbdz)*gmbm

c     compute GC force terms
      fm=gmbm*(bx*dbdx+by*dbdy+bz*dbdz)
      fc=p2gmb4*b*b*dsqrt(ax*ax+ay*ay+az*az)
      fg=gmbm*dsqrt(bxgb1*bxgb1+bxgb2*bxgb2+bxgb3*bxgb3)
      fe=dsqrt(ex*ex+ey*ey+ez*ez)

c     compute time step
      dt=tstep_gc*p/dsqrt(fm*fm+fc*fc+fg*fg+fe*fe)
      v=dsqrt(dydt(1)*dydt(1)+dydt(2)*dydt(2)+dydt(3)*dydt(3))
      dt=dmin1(dt,dx_max_gc/v)
*      print *,'dt,dt_max,dx=',real(dt/tfactor),
*     &real(dx_max_gc/v/tfactor),real(dt*v)

c     check for possible violation of 1st invariant
c     set go2lrntz approx. 0.05 in input file
c     go2lrntz set = go2lrntz^2 in gc_params
c     rho = p_perp/b = sqrt(2*mu0*b)/b = sqrt(2*mu0/b)
c     epsilon^2 = rho^2 * (dsum (db_i/dx_j)^2) / b2
c     switch = (2*mu0/b) * (dsum (db_i/dx_j)^2) / b2
      b3=b*b*b
      switch=2*mu*(dbxdx*dbxdx+dbxdy*dbxdy+dbxdz*dbxdz+dbydx*dbydx+
     &dbydy*dbydy+dbydz*dbydz+dbzdx*dbzdx+dbzdy*dbzdy+dbzdz*dbzdz)/b3
      eta=y(4)*2*pi*b/dsqrt(ax*ax+ay*ay+az*az)
!      if (eta.ge.etaswitch) then
!         print *, 'ETA VIOLATION', eta, etaswitch
!         if (flag_switch.eqv..true.) then
!            if (status.le.0) status = 4
!         else
!            if (status.le.0) status = 2
!         endif
!      endif
!      if (switch.ge.go2lrntz) then
!         print *, 'EPSILON VIOLATION', switch, go2lrntz
!         if (flag_switch.eqv..true.) then
!            if (status.le.0) status = 4
!         else
!            if (status.le.0) status = 2
!         endif
!      endif

      return
      end

************************************************************************

      subroutine gc2lrntz(y,t)

c        there may be problems with GC->Lorentz that lead to 
c        NANs & seg. fault (e.g., divide by very sm. B)

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      include 'rbelt-gcenter.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-flag.inc'
      real*8 gamma,w0,y(6),t
      real*8 bprime,ran0,phase,p_perp,p_parl,rho

      if (status.eq.4) status=0
      call get_fields2(y,t)

c     choose random gyrophase
      by=by+1.e-19
      bprime = dsqrt(bx*bx+by*by)
      phase = ran0(seed)*2*pi

c     calc. energy momentum & gyroradius
      w0=dsqrt(1.+2.*dabs(b)*mu+y(4)*y(4))-1.
      p_perp = dsqrt(w0*(2+w0)-y(4)*y(4))
      p_parl = y(4)
      rho = p_perp/b

*      print *
*      print *,'guiding center traj. :'
*      print *,'x,y,z,t=',y(1),y(2),y(3),t/tfactor
*      print *,'b(nT) =',b/ffactor/ntg
*      print *,'bx(nT)=',bx/ffactor/ntg
*      print *,'by(nT)=',by/ffactor/ntg
*      print *,'bz(nT)=',bz/ffactor/ntg
*      print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
*      print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
*      print *,'p(g*cm/sec)=',c*m0*sqrt(p_perp**2.+p_parl**2.)
*      gamma = 1+w0
*      print *,'p(g*cm/sec)=',c*m0*sqrt(gamma**2.-1)
*      print *,'rho=',rho
*      call get_fields2(y,t)
*      print *,'grad=',sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)
*      switch=sqrt(2*y(5)*(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)/b**3.)
*      print *,'switch=',switch

c     move position from guiding center to particle trajectory
      y(1)=y(1)+rho*(-by*dsin(phase)+bx*bz*dcos(phase)/b)/bprime
      y(2)=y(2)+rho*(bx*dsin(phase)+by*bz*dcos(phase)/b)/bprime
      y(3)=y(3)-rho*bprime*dcos(phase)/b

c     calc. cartesian components of proper velocity
      y(4)=(p_perp*(by*dcos(phase)+bz*bx*dsin(phase)/b)/
     &   bprime+p_parl*bx/b)
      y(5)=(p_perp*(-bx*dcos(phase)+bz*by*dsin(phase)/b)/
     &   bprime+p_parl*by/b)
      y(6)=(-p_perp*bprime*dsin(phase)/b+p_parl*bz/b)

      call get_fields2(y,t)

*      print *
*      print *,'Lorentz traj. :'
*      print *,'x,y,z,t=',y(1),y(2),y(3),t/tfactor
*      print *,'b(nT) =',b/ffactor/ntg
*      print *,'bx(nT)=',bx/ffactor/ntg
*      print *,'by(nT)=',by/ffactor/ntg
*      print *,'bz(nT)=',bz/ffactor/ntg
*      gamma = sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))
*      p_parl=(y(4)*bx+y(5)*by+y(6)*bz)/b
*      p_perp = sqrt((y(5)*bz-y(6)*by)**2.+
*     &(y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
*      print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
*      print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
*      print *,'p(g*cm/sec)=',c*m0*sqrt(p_perp**2.+p_parl**2.)
*      print *,'p(g*cm/sec)=',c*m0*sqrt(gamma**2.-1)
*      print *,'rho=',p_perp/b
*      print *,'grad=',sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)
*      switch=p_perp*sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)/b/b
*      print *,'switch=',switch
*      print *
**      stop

      return
      end

************************************************************************

      subroutine lrntz2gc(y,t)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'
      include 'rbelt-flag.inc'
      integer i
      real*8 y(6),t,p_perp,p_parl,pb,rho,x(3),p2,w0,gamma

      if (status.eq.-2) status=0
      call get_fields(y,t)

*      print *
*      print *,'Lorentz:'
*      print *,'b(nT) =',b/ffactor/ntg
*      print *,'bx(nT)=',bx/ffactor/ntg
*      print *,'by(nT)=',by/ffactor/ntg
*      print *,'bz(nT)=',bz/ffactor/ntg
*      p_parl=(y(4)*bx+y(5)*by+y(6)*bz)/b
*      p_perp=sqrt((y(5)*bz-y(6)*by)**2.+
*     &(y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
*      print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
*      print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
*      print *,'p(g*cm/sec)=',c*m0*sqrt(p_perp**2.+p_parl**2.)
*      print *,'rho=',p_perp/b
*      call get_fields2(y,t)
*      print *,'grad=',sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)
*      switch=p_perp*sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)/b/b
*      print *,'switch=',switch

c     get b field at/near guiding center position
      do 100 i = 1,10
         p_perp = dsqrt((y(5)*bz-y(6)*by)**2.+
     &   (y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
         rho = p_perp/b
         pb = p_perp*b
         x(1) = y(1) - rho*(by*y(6)-bz*y(5))/pb
         x(2) = y(2) - rho*(bz*y(4)-bx*y(6))/pb
         x(3) = y(3) - rho*(bx*y(5)-by*y(4))/pb
         call get_fields(x,t)
100   continue

c     move position vector from point on trajectory to guiding center
      y(1) = x(1)
      y(2) = x(2)
      y(3) = x(3)

      gamma=dsqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))
      y(4)=(y(4)*bx+y(5)*by+y(6)*bz)/b
      y(5)=gamma
      y(6)=p_perp*p_perp/b/2

*      call get_fields2(y,t)
*      print *
*      print *,'guiding center traj. :'
*      print *,'b(nT) =',b/ffactor/ntg
*      print *,'bx(nT)=',bx/ffactor/ntg
*      print *,'by(nT)=',by/ffactor/ntg
*      print *,'bz(nT)=',bz/ffactor/ntg
*      w0=sqrt(1.+2.*abs(b)*y(5)+y(4)*y(4))-1.
*      gamma = 1+w0
*      p_perp = sqrt(w0*(2+w0)-y(4)*y(4))
*      p_parl = y(4)
*      print *,'p_parl(g*cm/sec)=',c*m0*p_parl*charge_sign
*      print *,'p_perp(g*cm/sec)=',c*m0*p_perp*charge_sign
*      print *,'p(g*cm/sec)=',c*m0*sqrt(p_perp**2.+p_parl**2.)
*      print *,'rho=',p_perp/b
*      print *,'grad=',sqrt(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)
*      switch=sqrt(2*y(5)*(dbxdx**2.+dbxdy**2.+dbxdz**2.+
*     &dbydx**2.+dbydy**2.+dbydz**2.+
*     &dbzdx**2.+dbzdy**2.+dbzdz**2.)/b**3.)
*      print *,'switch=',switch
*      print *

      return
      end

************************************************************************

      subroutine lrntz2gc_chck(y,t,tchck)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-gcenter.inc'
*      include 'rbelt-flag.inc'
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-grid.inc'
      real*8 y(6),t,tchck
      real*8 invb2,invb4,vxb1,vxb2,vxb3,rho2
      real*8 p_parl,p_perp

      call get_fields2(y,t)
      invb2=1/(b*b)
      invb4=invb2*invb2
      vxb1=y(5)*bz-y(6)*by
      vxb2=y(6)*bx-y(4)*bz
      vxb3=y(4)*by-y(5)*bx
      rho2 = (vxb1*vxb1+vxb2*vxb2+vxb3*vxb3)*invb4
      switch=rho2*(dbxdx*dbxdx+dbxdy*dbxdy+dbxdz*dbxdz+dbydx*dbydx+
     &dbydy*dbydy+dbydz*dbydz+dbzdx*dbzdx+dbzdy*dbzdy+dbzdz*dbzdz)*invb2

*      p_perp=sqrt((y(5)*bz-y(6)*by)**2.+
*     & (y(6)*bx-y(4)*bz)**2.+(y(4)*by-y(5)*bx)**2.)/b
*      print *,'z,switch,mu=',y(3),sqrt(switch),
*     & abs(p_perp*p_perp/b/2*m0*c*c*ffactor)

      if (switch.le.go2gc) then
         print *,'switch.le.go2gc',dsqrt(switch),dsqrt(go2gc)
         print *,'going to gcenter traj. at t,r,tchck=',t/tfactor,
     &    dsqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3)),tchck/tfactor
         if (status.le.0) status = -2
      else
         tchck=tchck+dtgo2gc
      endif

      return
      end

************************************************************************

      subroutine rhs_cb(t,y,dydt)

      implicit none
      include 'rbelt-mu.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      real*8 y(4),dydt(4),t
      real*8 bxstar,bystar,bzstar,exstar,eystar,ezstar
      real*8 gamma,pparl_b,pparl_b2,bsparl,mu_gamma,pparl_gamma
      real*8 bxhat,byhat,bzhat,dbdtpp_b2

*      print *
*      print *,'  in rhs: y,t=',y,t/tfactor

c     update the fields to time t
      call get_fields2(y,t)

      gamma=dsqrt(1.+2.*dabs(b)*mu+y(4)*y(4))

      pparl_b=y(4)/b
      pparl_b2=pparl_b/b

      bxstar=bx + pparl_b*(dbzdy-dbydz) + pparl_b2*(by*dbdz-bz*dbdy)
      bystar=by + pparl_b*(dbxdz-dbzdx) + pparl_b2*(bz*dbdx-bx*dbdz)
      bzstar=bz + pparl_b*(dbydx-dbxdy) + pparl_b2*(bx*dbdy-by*dbdx)

      mu_gamma=mu/gamma
      dbdtpp_b2=dbdt*pparl_b2
      exstar=ex - mu_gamma*dBdx + pparl_b*dbxdt - dbdtpp_b2*bx
      eystar=ey - mu_gamma*dBdy + pparl_b*dbydt - dbdtpp_b2*by
      ezstar=ez - mu_gamma*dBdz + pparl_b*dbzdt - dbdtpp_b2*bz

      bxhat=bx/b
      byhat=by/b
      bzhat=bz/b
      pparl_gamma=y(4)/gamma
      bsparl=(bxhat*bxstar+byhat*bystar+bzhat*bzstar)
      dydt(1)=(pparl_gamma*bxstar + (eystar*bzhat-ezstar*byhat))/bsparl
      dydt(2)=(pparl_gamma*bystar + (ezstar*bxhat-exstar*bzhat))/bsparl
      dydt(3)=(pparl_gamma*bzstar + (exstar*byhat-eystar*bxhat))/bsparl
      dydt(4)=(exstar*bxstar+eystar*bystar+ezstar*bzstar)/bsparl

      return
      end

************************************************************************

      subroutine rhs_init_cb(t,y,dydt,dt)

      implicit none
      include 'rbelt-mu.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'

      include 'rbelt-gcenter.inc'
*      include 'rbelt-flag.inc'
      include 'rbelt-status.inc'
*      include 'rbelt-bounds.inc'
      include 'rbelt-dist.inc'

      real*8 y(4),dydt(4),t,dt
      real*8 bxstar,bystar,bzstar,exstar,eystar,ezstar
      real*8 gamma,pparl_b,pparl_b2,bsparl,mu_gamma,pparl_gamma
      real*8 ax,ay,az,p,v,b2,b3,fx,fy,fz,exbhatx,exbhaty,exbhatz
      real*8 bxhat,byhat,bzhat,dbdtpp_b2,bdgbx,bdgby,bdgbz

*      print *
*      print *,'  in rhs: y,t,dt=',y,t/tfactor,dt/tfactor

c     update the fields to time t
      call get_fields2(y,t)

      gamma=dsqrt(1.+2.*dabs(b)*mu+y(4)*y(4))

      pparl_b=y(4)/b
      pparl_b2=pparl_b/b

*      bxstar=bx+pparl_b*(dbzdy-dbydz) + pparl_b2*(by*dbdz-bz*dbdy)
*      bystar=by+pparl_b*(dbxdz-dbzdx) + pparl_b2*(bz*dbdx-bx*dbdz)
*      bzstar=bz+pparl_b*(dbydx-dbxdy) + pparl_b2*(bx*dbdy-by*dbdx)
      ax=pparl_b*(dbzdy-dbydz) + pparl_b2*(by*dbdz-bz*dbdy)
      ay=pparl_b*(dbxdz-dbzdx) + pparl_b2*(bz*dbdx-bx*dbdz)
      az=pparl_b*(dbydx-dbxdy) + pparl_b2*(bx*dbdy-by*dbdx)
      bxstar=bx+ax
      bystar=by+ay
      bzstar=bz+az

      mu_gamma=mu/gamma
      dbdtpp_b2=dbdt*pparl_b2
      exstar=Ex - mu_gamma*dBdx + pparl_b*dbxdt - dbdtpp_b2*bx
      eystar=Ey - mu_gamma*dBdy + pparl_b*dbydt - dbdtpp_b2*by
      ezstar=Ez - mu_gamma*dBdz + pparl_b*dbzdt - dbdtpp_b2*bz

      bxhat=bx/b
      byhat=by/b
      bzhat=bz/b

      exbhatx=(eystar*bzhat-ezstar*byhat)
      exbhaty=(ezstar*bxhat-exstar*bzhat)
      exbhatz=(exstar*byhat-eystar*bxhat)

      pparl_gamma=y(4)/gamma
      bsparl=(bxhat*bxstar+byhat*bystar+bzhat*bzstar)
      dydt(1)=(pparl_gamma*bxstar + exbhatx)/bsparl
      dydt(2)=(pparl_gamma*bystar + exbhaty)/bsparl
      dydt(3)=(pparl_gamma*bzstar + exbhatz)/bsparl
      dydt(4)=(exstar*bxstar+eystar*bystar+ezstar*bzstar)/bsparl

c     this is guiding center "force"
      fx=pparl_gamma*ax + exbhatx + dydt(4)*bxhat
      fy=pparl_gamma*ay + exbhaty + dydt(4)*byhat
      fz=pparl_gamma*az + exbhatz + dydt(4)*bzhat

*      stop

c     compute next time step
      p=dsqrt(gamma*gamma-1)
c     dt << p (total particle momentum) / dp/dt
      dt=tstep_gc*p/dsqrt(fx*fx+fy*fy+fz*fz)
*     (dt=tstep_gc*p/abs(dydt(4) ??)
      v=dsqrt(dydt(1)*dydt(1)+dydt(2)*dydt(2)+dydt(3)*dydt(3))
c     dt << x (where x is small w.r.t. system scale lengths) / dx/dt
      dt=dmin1(dt,(dx_max_gc/v))
!      print *,'dt=',dt/tfactor

c     check for possible violation of 1st invariant
c     set go2lrntz approx. 0.05 in input file
c     go2lrntz set = go2lrntz^2 in gc_params
c     rho = p_perp/b = sqrt(2*mu0*b)/b = sqrt(2*mu0/b)
c     epsilon^2 = rho^2 * (dsum (db_i/dx_j)^2) / b2
c     switch = (2*mu0/b) * (dsum (db_i/dx_j)^2) / b2

      b3=b*b*b
      switch=2*mu*(dbxdx*dbxdx+dbxdy*dbxdy+dbxdz*dbxdz+dbydx*dbydx+
     &dbydy*dbydy+dbydz*dbydz+dbzdx*dbzdx+dbzdy*dbzdy+dbzdz*dbzdz)/b3

      bdgbx = bx*dbxdx + by*dbxdy + bz*dbxdz
      bdgby = bx*dbydx + by*dbydy + bz*dbydz
      bdgbz = bx*dbzdx + by*dbzdy + bz*dbzdz
      eta = pparl_b2*dsqrt(bdgbx*bdgbx+bdgby*bdgby+bdgbz*bdgbz)*gamma/b

!      print *,'eta=',eta

      if (switch.ge.go2lrntz) then
         print *,'switch.ge.go2lrntz',dsqrt(switch),dsqrt(go2lrntz),
     &t/tfactor
         if (flag_switch.eqv..true.) then
            if (status.le.0) status = 4
         else
            if (status.le.0) status = 5
         endif
      endif
      if (eta.ge.etaswitch) then
         print *,'eta.ge.etaswitch',eta,etaswitch,t/tfactor
         if (flag_switch.eqv..true.) then
            if (status.le.0) status = 4
         else
            if (status.le.0) status = 5
         endif
      endif


      return
      end
