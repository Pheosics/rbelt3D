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
      integer n
      real t2,ran0

      NAMELIST/dist/dist_seed0,dt_dist,init_t,emin,emax,epa_min,epa_max,
     &lmin,lmax,radius,exp,factor,flag0,flag_switch
     
      print *
      print *,'*** in subroutine dist_particles ***'
      print *,'num_particles=',num_particles

c     rbelt-dist_tail.f
      initdist=5

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
      dist_seed=dist_seed0

      if (num_pdata.lt.10) then
         print *,'in y0.inc, num_pdata must be >=10'
         stop
      endif

      do n = 1,num_particles
      
c        particle start time
         y0(10,n) = tzero1+init_t+ran0(dist_seed)*dt_dist
         y0(7,n) = y0(10,n)

c        status=99 until we assign all particle data in dist_particles2
         int_y0(1,n) = 99
	 
c        flag
         int_y0(2,n) = flag0

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
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-const.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-io.inc'
      include 'rbelt-status.inc'
      integer file_num,n,j
      real x(3),t,t2,r,lshell,epa,phi,ke,ran0,theta
      real alpha,beta,gamma,p,xoffset,px,py,pz,bphi,btheta

c     for ONERA Lm & L* calc
      real bmin,bloc,xbmin(3)
      integer ifail

      parameter(xoffset=8.)
     
      print *
      print *,'**** in subroutine dist_particles2 ***'

      do n = 1,num_particles

         t=y0(10,n)-tzero

         if ((int_y0(1,n).eq.99).and.(t.ge.tgr(1)).and.(t.lt.t2)) then
c           reset status
            int_y0(1,n)=0
c           particle position
100         continue
c              use lmin=16, lmax=28
c              locate initial field line
               r = lmin+(lmax-lmin)*ran0(dist_seed)**.5
               phi = pi*(0.777777 + 0.444444*ran0(dist_seed))
               theta=pi/2.0
               x(1) = r*sin(theta)*cos(phi)+xoffset
               x(2) = r*sin(theta)*sin(phi)
               x(3) = r*cos(theta)
c              find bmin location along field line passing through initial point
               call gbdbzero(x,t,bloc,bmin,xbmin,ifail)
*               call find_bmin(x,t,Bloc,Bmin,xbmin,ifail)
            if (ifail.ne.0) goto 100
c           bmin surface pitch angle
            epa=epa_min+ran0(dist_seed)*(epa_max-epa_min)
c           initial kinetic energy
            ke = emin+ran0(dist_seed)*(emax-emin)
            y0(1,n) = xbmin(1)
            y0(2,n) = xbmin(2)
            y0(3,n) = xbmin(3)
c           get fields at initial position
            status=0
            call get_fields(xbmin,t)
c           for Lorentz trajectories
            if (flag0.eq.0) then
c              in B-field alligned coordinate system
               gamma = 1.+ ke
               p = sqrt(gamma*gamma-1.)
               alpha=epa
               beta=2*pi*ran0(dist_seed)
               call sphr2cart(p,alpha,beta,px,py,pz)
c              rotate to x-y-z system
               call rotate_back(px,py,pz,bx,by,bz,bphi,btheta)
               y0(4,n) = px
               y0(5,n) = py
               y0(6,n) = pz
c           for guiding center trajectories
            elseif (flag0.ge.1) then
c              set parallel momentum
               y0(4,n) = sqrt(ke*(ke+2.))*cos(epa)*charge_sign
c              set kinetic energy
               y0(5,n) = ke
c              calc. y0(6,n) = mu
               y0(6,n) = (y0(5,n)*(y0(5,n)+2.)-y0(4,n)*y0(4,n))/2/b
            else
               print *,'bad value for flag0:',flag0
               stop
            endif
c           initial quantities to weight particles in terms of go here
c           these get written to output file if init_out = .false. (in input file)
            y0(8,n) = ke
            y0(9,n) = epa
         endif
	 
      enddo

      return
      end

************************************************************************
