
      include 'rbelt-const.inc'

c     pstprc namelist
      integer mincnt
      real norm,spctrm,ecutoff
      common /pstprc_params/ norm,spctrm,ecutoff,mincnt

c     rbelt code parameters to read from info file
      logical binio,flux_out,dist_out,prcp_out,init_out,cone_out
      integer initdist,year0,doy0,hour0,min0,sec0,
     &wlines_dist,wlines_flux,wlines_prcp,wlines_cone
      real radius,expr,factor,dt_dist,emin,emax,tmax,flux_dt,
     &init_twrite,dtwrite,rmin,rmax
      common /utime/year0,doy0,hour0,min0,sec0
      common /rbelt_params/ radius,expr,factor,dt_dist,emin,emax,
     &rmin,rmax,tmax,flux_dt,init_twrite,dtwrite,initdist,wlines_dist,
     &wlines_flux,wlines_prcp,wlines_cone,binio,flux_out,dist_out,
     &prcp_out,init_out,cone_out

************************************************************************

      integer method
      parameter(method=1)

************************************************************************

c     particle initial conditions
      integer num_particles,num_files,fcount
c     put number of particles from line 10 of info file here
      parameter(num_particles=1000)

c     number of rbelt output files (runs)
c     set num_files=(lastfilenum-firstfilenum)/filenumstep+1
      parameter(num_files=246)

c     use if init_out=.true.
      real pwght(num_particles,num_files)
c     do not use pwght if init_out=.false. and you are using y0(8-10,i)
c     to store initial particle data
*      real pwght(1,1)
c
	real pgrid1(num_particles,3)
	real pgrid2(num_particles,3)
	common /weight/ fcount,pwght,pgrid1,pgrid2

************************************************************************

c     put number of write steps from line 59 of info file here
      integer num_wsteps
      logical prcp
      parameter(num_wsteps=4134)

      integer num_tsteps
      real init_time,dtstep

c     equatorial flux is processed on rbelt wtime grid
c
c     for combining rbelt num_wsteps into single distribution time steps
c     (i.e., we end up with num_wsteps/timeavg time steps in distributions)
c     (num_wsteps must be evenly divisible by timeavg)
      integer stepavg
      parameter(stepavg=10)
      parameter(num_tsteps=num_wsteps/stepavg)
c
c     pricipitation flux is processed on pstprc time grid defined here
      parameter(prcp=.false.)
*      parameter(num_tsteps=300)
      parameter(init_time=0.0,dtstep=1.0)

      real wtime(num_wsteps)
      real time(num_tsteps)
      integer itstep,initstep
      parameter(initstep=1,itstep=5)
      common /tgrid/ wtime,time

************************************************************************

c     L-shell, energy, PA distribution
c     ********************************
c     na needs to be even
      integer nl,ne,na,nb
c     number of PA grid points should be even
c     must use nl=23, ne=26 for ESA model weighting
*      parameter (nl=23,ne=26,na=28)
      parameter (nl=12,ne=5,na=18,nb=7)

      real ls(nl),ke(ne),ca(na),lsmin,lsmax,kemin,kemax,camin,camax
      real bm(nb),bmmin,bmmax
      parameter (lsmin=2.0,lsmax=3.1,kemin=26.0,kemax=58.0)
c     need to fix boundaries and/or interpolation method for PA bins/grid.
c     the problem is that we get divide by zero if we include PAs 0, 90 or 180
c     degrees. The choices are:
c     1. exclude counts between 0 and min PA and max PA and 180 (done currently)
c     2. modify grids to include dx/2 outside of boundary grid point.
c     (better solution, but problems with interpolation routines to be resolved)
c     do not change amin & amax
      parameter (camin=-1.+.1/na,camax=1.-.1/na,bmmin=2000.,bmmax=5000.)
      common /lea_grid/ ls,ke,ca,bm
      real n_leat(nl,ne,na,num_tsteps),j_leat(nl,ne,na,num_tsteps)
*      real n_lea0(nl,ne,na),j_lea0(nl,ne,na)
      common /j_lea_dist/ n_leat,j_leat
*      common /j_lea_dist0/ n_lea0,j_lea0

c     set flux map grid parameters
c     ****************************
c     x-y-z distribution grid parameters
c     set x-y-z (in SM coordinates equatorial plane) grid parameters 
      integer nx,ny
*      parameter (nx=101,ny=101)
      parameter (nx=1,ny=1)

      real x(nx),y(ny),xmin,xmax,ymin,ymax
      parameter (xmin=-10.,xmax=10.)
      parameter (ymin=-10.,ymax=10.)

      real j_xyt(nx,ny,num_wsteps),n_xyt(nx,ny,num_wsteps)
      common /xyz_grid_1/ x,y, j_xyt, n_xyt
    
c     set flux map grid parameters for prcp
c     *************************************
c     set theta-phi (in SM coordinates) grid parameters 
      integer ntht,nphi
      real thtmin,thtmax,phimin,phimax
c
c     uniformly spaced in cos(theta)
*      parameter (ntht=240)
*      parameter (thtmin=-1+1./ntht,thtmax=1.-1./ntht)
*      parameter (phimin=-pi,phimax=pi)
*      parameter (nphi=nint(4*pi/(dacos(1-1./ntht)+
*     &dacos(1-1./ntht-((thtmax-thtmin)/(ntht-1))))))
*      real tht(ntht),phi(nphi)
*      common /tht_phi/ tht,phi
c
c     uniformly spaced in theta with phi grid e,g, 0->357.0000
*      parameter (ntht=68+1,nphi=120)
*      parameter (thtmin=5./raddeg,thtmax=175./raddeg)
*      parameter (phimin=0,phimax=2*pi-3./raddeg)
*      real tht(ntht),phi(nphi)
*      common /tht_phi/ tht,phi
c
c     uniformly spaced in theta with phi grid e,g, 0->360.0000

*      parameter (ntht=68+1,nphi=121)
      parameter (ntht=1,nphi=1)
      parameter (thtmin=5./raddeg,thtmax=175./raddeg)
      parameter (phimin=-pi,phimax=pi)
      real tht(ntht),phi(nphi)
      common /tht_phi/ tht,phi
c
      real j_tpt(ntht,nphi,num_tsteps),n_tpt(ntht,nphi,num_tsteps)
      common /flux/ j_tpt,n_tpt
c     *************************************

************************************************************************

c     particle weighting function (wf) grid
c     *************************************

c     for SEP flux (init. energy, init. cos(PA), init. time)

      integer n1,n2,n3,nd,n4
*      parameter (n1=16,n2=4,n3=16,nd=3)
      parameter (n1=nl,n2=ne,n3=na,n4=nb,nd=num_tsteps)
*      parameter (n1=12,n2=4,n3=26,nd=num_tsteps)

      real a1(n1),a2(n2),a3(n3),min1,max1,min2,max2,min3,max3
      real a4(n4),min4,max4

      parameter (min1=2.0,max1=3.1,min2=26.0,max2=58.0)
      parameter (min3=-1.+.1/na,max3=1.-.1/na,min4=2000.,max4=5000.)

*      parameter (min1=-1.0,max1=2.0,min2=-1.,max2=1.)
*      parameter (min3=-1.,max3=1000000.0)

      real wf(n1,n2,n3)
      common /wghtfn/ a1,a2,a3,a4,wf





