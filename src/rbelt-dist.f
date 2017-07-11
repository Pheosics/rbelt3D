c     rbelt-dist.f - initializes a single particle

***********************************************************************

      subroutine dist_particles()

c test/sample distribute particle routine
c puts a single guiding center or Lorentz particle on the space 
c and time grid using parameters specified in the input file.

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-const.inc'
      include 'rbelt-fields.inc'
      integer file_num
      real r,theta,phi,alpha,beta,ke,epa
      real gamma,p,x(3),t
      NAMELIST /dist_1p/ r,theta,phi,alpha,beta,ke,epa,flag0,flag_switch

c     set defalt values for some input parameters
      dist_seed=1

c     read in and normalize input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,dist_1p)
      CLOSE (81)
c     normalize
      ke = ke/wrest
      epa = epa/raddeg
      theta = theta/raddeg
      phi = phi/raddeg
      alpha = alpha/raddeg
      beta = beta/raddeg

      print *
      print *,'*** in subroutine dist_particles ***'

      initdist=1
      dist_seed=1
      y0(1,1) = r*sin(theta)*cos(phi)
      y0(2,1) = r*sin(theta)*sin(phi)
      y0(3,1) = r*cos(theta)
      y0(7,1) = tzero1
      int_y0(1,1) = 0
      int_y0(2,1) = flag0
      gamma = 1.+ ke
      if (flag0.eq.0) then
         p = sqrt(gamma*gamma-1.)
         y0(4,1) = p*sin(alpha)*cos(beta)
         y0(5,1) = p*sin(alpha)*sin(beta)
         y0(6,1) = p*cos(alpha)
      elseif (flag0.ge.1) then
c        NOTE, if we do not start with theta=90 deg. 
c        then epa will be the initial (not equatorial) pitch angle
         y0(4,1) = sqrt(ke*(ke+2.))*cos(epa)*charge_sign
c        (might be better to work with KE instead of gamma)
         y0(5,1) = gamma
c        calc y0(6,1) = 1st adiabatic invariant mu
         x(1) = y0(1,1)
         x(2) = y0(2,1)
         x(3) = y0(3,1)
         t=y0(7,1)-tzero
         call get_fields(x,t)
         y0(6,1) = (ke*(1.+ke/2.)-y0(4,1)*y0(4,1)/2)/b
      endif

      if (num_pdata.ge.10) then
         y0(8,1)=ke*wrest
         y0(9,1)=b/ffactor/ntg
         y0(10,1)=0.0
      endif

      return
      end

***********************************************************************

      subroutine dist_particles2()
      return
      end
