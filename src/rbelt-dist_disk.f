c     rbelt-dist_disk.f - launches particles from a disk in the equatorial plane (in SM coordinates)

***********************************************************************

      subroutine dist_particles(t2)

c distribute particles routine
c initializes guiding center or Lorentz particles 

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-const.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-io.inc'
      integer file_num,n,j,seed
      real*8 t2,x(3),t,r,lshell,epa,phi,ke,ran0,theta
      real*8 alpha,beta,gamma,p

      NAMELIST/dist/dist_seed,dt_dist,init_t,emin,emax,epa_min,epa_max,
     &lmin,lmax,radius,exp,factor,flag0,flag_switch
     
!      print *
!      print *,'*** in subroutine dist_particles ***'
!      print *,'num_particles=',num_particles

      initdist=2
c     read in and normalize input parameters
      print *,'reading namelist'
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,dist)
      CLOSE (81)
c     normalize
      epa_min = epa_min/raddeg
      epa_max = epa_max/raddeg
      emin = emin/wrest
      emax = emax/wrest
      dt_dist = dt_dist*tfactor
      init_t = init_t*tfactor
      seed=dist_seed

      if ((init_t+dt_dist).gt.t2) then
         print *,'(init_t+dt_dist).gt.t2:',(init_t+dt_dist),t2
         print *,'decrease (init_t+dt_dist) or increase tmax'
         stop
      endif

      do n = 1,num_particles
         lshell = lmin+(lmax-lmin)*ran0(seed)**.5
         r = lshell
         phi=ran0(seed)*2*pi
         theta=pi/2.0
         epa=epa_min+ran0(seed)*(epa_max-epa_min)
         ke = emin+ran0(seed)*(emax-emin)
c        particle start time
         y0(7,n) = tzero1+init_t+ran0(seed)*dt_dist
c        particle position
         y0(1,n) = r*dsin(theta)*dcos(phi)
         y0(2,n) = r*dsin(theta)*dsin(phi)
         y0(3,n) = r*dcos(theta)
c        status
         int_y0(1,n) = 0
c        flag
         int_y0(2,n) = flag0
c        for Lorentz trajectories
         if (flag0.eq.0) then
c           particle is launched with equatorial component of its velocity 
c           radially outward to put its guiding center on approximatly the 
c           correct L-shell 
            beta=phi
            alpha=epa
c           compute proper velocity
            gamma = 1.+ ke
            p = dsqrt(gamma*gamma-1.)
c           it is here assumed that b is in z dir. (in equatorial plane)
            y0(4,n) = p*dsin(alpha)*dcos(beta)
            y0(5,n) = p*dsin(alpha)*dsin(beta)
            y0(6,n) = p*dcos(alpha)
c        for guiding center trajectories
         elseif (flag0.ge.1) then
c           set parallel momentum
            y0(4,n) = dsqrt(ke*(ke+2.))*dcos(epa)*charge_sign
c           set kinetic energy
            y0(5,n) = ke
c           still need to call get_fields and calc. y0(6,n) = mu
c           done in dist_particles2() below
         endif
c        initial quantities to weight particles in terms of go here
c        these get written to output file if init_out = .false. (in input file)
         if (num_pdata.ge.10) then
            y0(8,n) = ke*wrest
            y0(9,n) = epa*raddeg
            y0(10,n) = y0(7,n)/tfactor
         endif

      enddo

      return
      end

***********************************************************************

      subroutine dist_particles2(t2)

c test/sample distribute particle routine
c puts a single guiding center or Lorentz particle on the space 
c and time grid using parameters specified in the input file.

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-bounds.inc'
      integer n
      real*8 x(3),t,t2
     
      print *
      print *,'*** in subroutine dist_particles2 ***'

      do n = 1,num_particles
         if (flag0.ge.1) then
            t=y0(7,n)-tzero
            if ((t.ge.tgr(1)).and.(t.lt.t2)) then
c              calc y0(5,1) = 1st adiabatic invariant mu
               x(1) = y0(1,n)
               x(2) = y0(2,n)
               x(3) = y0(3,n)
               call get_fields2(x,t)
               y0(6,n) = (y0(5,n)*(y0(5,n)+2.)-y0(4,n)*y0(4,n))/2/b
            endif
         endif
      enddo

      return
      end
