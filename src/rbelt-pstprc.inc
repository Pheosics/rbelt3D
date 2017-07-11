
      include 'rbelt-const.inc'

c     pstprc namelist
      integer mincnt
      real norm,spctrm,ecutoff
      common /pstprc_params/ norm,spctrm,ecutoff,mincnt

c     rbelt code parameters to read from info file
      logical binio,flux_out,dist_out,prcp_out,init_out,cone_out
      integer initdist,year0,doy0,hour0,min0,sec0,
     &wlines_dist,wlines_flux,wlines_prcp,wlines_cone
      real radius,exp,factor,dt_dist,emin,emax,tmax,flux_dt,init_twrite,
     &dtwrite,rmin,rmax,lmin,lmax
      common /utime/year0,doy0,hour0,min0,sec0
      common /rbelt_params/ radius,exp,factor,dt_dist,emin,emax,
     &rmin,rmax,lmin,lmax,tmax,flux_dt,init_twrite,dtwrite,initdist,
     &wlines_dist,wlines_flux,wlines_prcp,wlines_cone,binio,flux_out,
     &dist_out,prcp_out,init_out,cone_out

************************************************************************

      integer method
      parameter(method=0)

************************************************************************

c     particle initial conditions
      integer num_particles,num_files,fcount,numuwc
c     put number of particles from line 10 of info file here
      parameter(num_particles=100000)

c     number of rbelt output files (runs)
c     set num_files=(lastfilenum-firstfilenum)/filenumstep+1
      parameter(num_files=221)

c     array for weighting each particle
      integer pwghtflag(num_particles,num_files)
      real pwght(num_particles,num_files)

c     do not use pwght if using method 3
*      real pwght(1,1)

      common /weight/ fcount,numuwc,pwghtflag,pwght

************************************************************************

c     put number of write steps from line 59 of info file here
      integer num_wsteps,fullsteps
      parameter(num_wsteps=2155)
      real wtime(num_wsteps)
      common /wgrid/ wtime

c     equatorial flux is processed on rbelt wtime grid
c     for combining rbelt num_wsteps into single distribution time steps
c     (i.e., we end up with num_wsteps/stepavg time steps in distributions)
c     (num_wsteps must be evenly divisible by timeavg to include all steps)
      integer stepavg
      parameter(stepavg=20)

      integer nt
      parameter(nt=num_wsteps/stepavg+1)
      parameter(fullsteps=num_wsteps/stepavg)

      real tg(nt)
      common /tgrid/ tg

************************************************************************

c     L-shell, energy, EPA, time distribution
c     ***************************************
      integer nl,ne,na
c     number of PA grid points (na) should be even
c     must use nl=23, ne=26 for ESA model weighting
*      parameter (nl=23,ne=26,na=28)
*      parameter (nl=16,ne=10,na=16)
      parameter (nl=6,ne=20,na=16)
      real ls(nl),ke(ne),ca(na),lsmin,lsmax,kemin,kemax,camin,camax

c     need to fix boundaries and/or interpolation method for PA bins/grid.
c     the problem is that we get divide by zero if we include PAs 0, 90 or 180
c     degrees. The choices are:
c     1. exclude counts between 0 and min PA and max PA and 180 (done currently)
c     2. modify grids to include dx/2 outside of boundary grid point.
c     (better solution, but problems with interpolation routines to be resolved)
c     do not change amin & amax
      parameter (lsmin=3.5,lsmax=8.5)
      parameter (kemin=0.1,kemax=2.0)
      parameter (camin=-1.+1./na,camax=1.-1./na)

      common /lea_grid/ ls,ke,ca

      real n_leat(nl,ne,na,nt),j_leat(nl,ne,na,nt)
      common /j_lea_dist/ n_leat,j_leat


c     lower bound for energy integrated flux
c     **************************************
      real elwr
      parameter (elwr=0.8)


c     particle weighting function (wf) grid
c     *************************************
      real wf(nl,ne,na,nt)
      common /wghtfn/ wf


c     set equatorial flux grid parameters
c     ***********************************
c     x-y-z distribution grid parameters
c     set x-y-z (in SM coordinates equatorial plane) grid parameters 
      integer nx,ny
*      parameter (nx=1,ny=1)
      parameter (nx=176,ny=151)
*      parameter (nx=138,ny=121)
*      parameter (nx=101,ny=91)
*      parameter (nx=71,ny=61)

      real x(nx),y(ny),xmin,xmax,ymin,ymax
      parameter (xmin=-25.,xmax=10.)
      parameter (ymin=-15.,ymax=15.)

      real j_xyt(nx,ny,nt),n_xyt(nx,ny,nt)
      common /xyz_grid/ x,y,j_xyt,n_xyt
    

c     set precipitation flux grid parameters
c     **************************************
c     set theta-phi (in SM coordinates) grid parameters 
      integer ntht,nphi
      real thtmin,thtmax,phimin,phimax

c     pricipitation flux is processed on pstprc time grid defined here
      real init_time,dtstep
      parameter(init_time=0.0,dtstep=100.0)
c
c     uniformly spaced in cos(theta)
*      parameter (ntht=240)
*      parameter (thtmin=-1+1./ntht,thtmax=1.-1./ntht)
*      parameter (phimin=-pi,phimax=pi)
*      parameter (nphi=nint(4*pi/(acos(1-1./ntht)+
*     &acos(1-1./ntht-((thtmax-thtmin)/(ntht-1))))))
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
      real j_tpt(ntht,nphi,nt),n_tpt(ntht,nphi,nt)
      common /flux/ j_tpt,n_tpt
c     **************************************

************************************************************************
