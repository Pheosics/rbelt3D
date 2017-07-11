c rbelt-main.f
c main program for computing particle trajectories in time dependent fields
c the first field file must have nstep time steps -- 
c 2 more time steps that subsequent field files.
c this code probably should be modified to read field files with one time 
c step per file (write another version of main).

c***********************************************************************
c     PROGRAM MAIN
c
c This is a sample main program for use with rbelt code routines.
c This version advances a distribution of identical particles non-
c self-consistantly in time-debendent E & B fields input from user
c specified hdf file(s). See the README.text file. B.Kress -- 2/27/06
c
c***********************************************************************

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      character*80 basename
      character*80 fieldfile,gridfile,outfile,filename
      integer filenum,firstfilenum,lastfilenum,filenumstep
      NAMELIST /rbelt/ basename,firstfilenum,lastfilenum,filenumstep
      integer*4 today(3), now(3)
      real*8 t2

      print *
      print *,'*** IN PROGRAM MAIN ***'

c     print date & time
!      call idate(today)
!      call itime(now)
!      print *
!      print *,'system day, month, year =',today
!      print *,'system hour, min, sec =',now

c     make sure we have the right time grid
      if (nstep.le.2) then
         print *,'nstep=',nstep,'   (should be >= 3)'
         print *,'wrong grid size defined in rbelt-grid.inc'
         stop
      endif

c     read in and normalize rbelt input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,rbelt)
      CLOSE (81)

c     get UT reference date/time
      call load_ut(basename)

c     read spatial grid
      filename=gridfile(basename,firstfilenum)
      call load_grid(filename)

c     read & scale time grid and fields
      filename=fieldfile(basename,firstfilenum)
      call load_field(filename,firstfilenum,1,nstep,1)
      call scale_fields(1,nstep)

c     time and space grid & boundaries initialization
c     tzero1 set = tgr(1) in grid_init
      call grid_init()
      call bounds_init()

c     initialize field quantities to zero
      call init_fields

c     calculate grad B (must calc. dx,dy,dz in grid_init first)
      call calc_derivs(1,nstep)
      call calc_ddt(1,nstep)
!      call precalcmatrix(1,nstep)

c     variable tzero is set = to tgr(1) (time with respcet to UT ref. time)
c     zero first step of time grid by subtracting tgr(1) form all time steps in grid
      tgrmin=tgr(1)
      call zero_tgrid(tgrmin)
      tgrmax=tgr(nstep)

c     initialize integrators
      call lrntz_params()
      call gc_params()
      call rkckparm

c     set upper limit on time
      t2=dmin1((tmax-tzero),tgrmax)

c     distribute particles
      call dist_particles(t2)
      call dist_particles2(t2)

c     initialize i/o routines
      call io_init(basename,firstfilenum)
!      call openfile(basename,firstfilenum)

c     print date & time
!      call idate(today)
!      call itime(now)
!      print *
!      print *,'run particles hour, min, sec =',now

c     loop through and advance particles & output data
      call particle_loop(t2)

c     print date & time
!      call idate(today)
!      call itime(now)
!      print *
!      print *,'finished running particles hour, min, sec =',now

c     close for multiple output files
!      call closefile()

c     loop over input field files
!      print *
!      print *,'*** IN PROGRAM MAIN ***'
!      print *,'first,last,step=',(firstfilenum+nstep*filenumstep),
!     &lastfilenum,(nstep-2)*filenumstep

      do filenum=(firstfilenum+filenumstep),
     &lastfilenum,filenumstep
         print *
         print *,'*** IN PROGRAM MAIN ***'
         print *,'filenum =',filenum
	 
c        read fields
c        put last 2 time steps into first 2 time grid posiitons
         call last2_2first2()
c        read in 3rd time step through nstep steps
         filename=fieldfile(basename,filenum)
         call load_field(filename,filenum,1,(nstep-2),3)
         call scale_fields(3,nstep)
         call calc_derivs(3,nstep)
         call calc_ddt(2,nstep)

         call subtrct_tzero(3,nstep)
c        zero each time grid
         tgrmin=tgr(1)
         call zero_tgrid(tgrmin)
         tgrmax=tgr(nstep)

         t2=dmin1((tmax-tzero),tgrmax)
c        calc. field dependent particle data
!         call dist_particles2(t2)
c        to open multiple output files
!         call openfile(basename,firstfilenum)
c        advance particles & output data
!         t2=dmin1((tmax-tzero),tgrmax)
         call particle_loop(t2)
c        close for multiple output files
!         call closefile()
      enddo

c     finish up with i/o routines
      call io_close(basename,firstfilenum)
      print *

      end
