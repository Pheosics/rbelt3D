c     rbelt-main1t.f
c     main program for computing particle trajectories in static fields
c     or time dependent analytic field 

c***********************************************************************
c     PROGRAM MAIN
c
c This is a sample main program for use with rbelt code routines.
c This version advances a distribution of identical particles non-
c self-consistantly in static E & B fields input from a user
c specified hdf file. See the README.text file. B.Kress -- 2/27/06
c
c***********************************************************************

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      include 'rbelt-y0.inc'
      character*80 fieldfile,gridfile,outfile,basename
      integer filenum,firstfilenum,lastfilenum,filenumstep,timesteps
      NAMELIST /rbelt/ basename,firstfilenum,lastfilenum,filenumstep
      integer*4 today(3), now(3)
      real*8 t2

      print *
      print *,'*** IN PROGRAM MAIN ***'

c     print date & time
*      call idate(today)
*      call itime(now)
*      print *
*      print *,'system day, month, year =',today
*      print *,'system hour, min, sec =',now

c     make sure we have the right time grid
      if (nt.ne.1) then
         print *,'nt=',nt,'   (should be nt=1)'
         print *,'wrong grid size defined in rbelt-grid.inc'
         stop
      endif

c     read in and normalize rbelt input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,rbelt)
      CLOSE (81)

c     get UT reference date/time
      call load_ut(basename,firstfilenum,lastfilenum,filenumstep)

c     read spatial grid
      call load_grid(basename)

c     read & scale time grid and fields
      call load_field(basename,firstfilenum,1,1,1)
      call rbelt_recalc(tgr(1))
      call scale_fields(1,nt)

c     time and space grid & boundaries initialization
      call grid_init() ! initializes tzero=0 & sets tzero1=tgr(1)
      call bounds_init() ! normalizes tmax = tmax*tfactor+tzero1

c     initialize field quantities to zero
      call init_fields

c     calculate grad B (must calc. dx,dy,dz in grid_init first)
      call calc_derivs(1,nt)

c     zero time step (here tgrmax becomes tmax)
      tmin=tgr(1)
      call zero_tgrid(tmin) 
      tgrmax=tmax-tzero

c     initialize integrators
      call lrntz_params()
      call gc_params()
      call rkckparm

c     set upper limit on time
      t2=tgrmax

c     distribute particles
      call dist_particles(t2)
c     calc. field dependent particle data
      call dist_particles2(t2)

c     initialize i/o routines
      call io_init(basename,firstfilenum)
c     field dataset i/o
!      call field_io(basename,firstfilenum,firstfilenum,timesteps)
c     to open multiple output files
!      call openfile(basename,firstfilenum)

c     loop through and advance particles & output data
      call particle_loop(t2)

c     print date & time
*      call idate(today)
*      call itime(now)
*      print *
*      print *,'finished running particles hour, min, sec =',now

c     close for multiple output files
!      call closefile()

c     loop over input field files
      print *
      print *,'*** IN PROGRAM MAIN ***'
      print *,'first,last,step=',firstfilenum,
     &lastfilenum,filenumstep
      do filenum=firstfilenum+filenumstep,lastfilenum,
     &   filenumstep
         print *
         print *,'*** IN PROGRAM MAIN ***'
         print *,'filenum =',filenum
         call load_field(basename,filenum,1,1,1)
         call rbelt_recalc(tgr(1))
         call scale_fields(1,nt)
         call calc_derivs(1,nt)

c        field dataset i/o
!         call field_io(basename,filenum,firstfilenum,timesteps)
c        to open multiple output files
!         call openfile(basename,filenum)

c        zero time step (here tgrmax becomes tmax)
         tmin=tgr(1)
         call zero_tgrid(tmin)

c        calc. field dependent particle data
         call dist_particles2(t2)
c        advance particles & output data
         t2=tmax-tzero
         call particle_loop(t2)
c        close for multiple output files
!         call closefile()
      enddo

c     finish up with i/o routines
      call io_close(basename,firstfilenum)
      print *

      end




