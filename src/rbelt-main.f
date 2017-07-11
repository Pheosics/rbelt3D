c rbelt-main.f
c main program for computing particle trajectories in time dependent fields
c currently there are two formats for input fields:
c 1) rbelt HDF files with nt+2 time steps in 1st and nt time steps in subsequent
c 2) directly from LFM (ouput) HDF files, with 1 time step in each 
c in either case, we must supply lastfilenum+(nt-2)-firstfilenum time steps.
c (i.e., nt-3 additional time steps beyond lastfilenum) 

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
      integer filenum,firstfilenum,lastfilenum,filenumstep,timesteps
      NAMELIST /rbelt/ basename,firstfilenum,lastfilenum,filenumstep
      integer*4 today(3), now(3)
      real t2

      print *
      print *,'*** IN PROGRAM MAIN ***'

c     print date & time
      call idate(today)
      call itime(now)
      print *
      print *,'system day, month, year =',today
      print *,'system hour, min, sec =',now

c     make sure we have the right time grid
      if (nt.le.2) then
         print *,'nt=',nt,'   (should be >= 3)'
         print *,'wrong grid size defined in rbelt-grid.inc'
         stop
      endif

c     read in and normalize rbelt input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,rbelt)
      CLOSE (81)

c     make sure input file numbering is correct
c     (last time step to read is lastfilenum+(nt-2)-1)
      timesteps=(lastfilenum-firstfilenum+1)*(nt-2)+2
      print *
      print *,'total number of time steps =',timesteps
      if (mod((timesteps-2),(nt-2)).ne.0) then
         print *,'mod((timesteps-2),(nt-2)).ne.0'
         print *,'(timesteps-2),(nt-2)=',(timesteps-2),(nt-2)
         print *,'total number of time steps should be'
         print *,'nt + integer*(nt-2)'
         stop
      endif

c     get UT reference date/time
      call load_ut(basename,firstfilenum,lastfilenum,filenumstep)

c     read spatial grid
      call load_grid(basename)

c     read & scale time grid and fields
      call load_field(basename,firstfilenum,1,nt,1)
      call scale_fields(1,nt)

c     time and space grid & boundaries initialization
c     tzero1 set = tgr(1) in grid_init
      call grid_init()
      call bounds_init()

c     initialize field quantities to zero
      call init_fields

c     calculate grad B (must calc. dx,dy,dz in grid_init first)
      call calc_derivs(1,nt)
      call calc_ddt(1,nt)
!      call calc_e(1,nt)

c     variable tzero is set = to tgr(1) (time with respcet to UT ref. time)
c     zero first step of time grid by subtracting tgr(1) from all time steps in grid
      tgrmin=tgr(1)
      call zero_tgrid(tgrmin)
      tgrmax=tgr(nt)

c     initialize integrators
      call lrntz_params()
      call gc_params()
      call rkckparm

c     set upper limit on time
      t2=amin1((tmax-tzero),tgrmax)

c     distribute particles
      call dist_particles(t2)
      call dist_particles2(t2)

c     initialize i/o routines
      call io_init(basename,firstfilenum)
c     field dataset i/o
!      call field_io(basename,firstfilenum,firstfilenum,timesteps)
c     to open multiple output files
!      call openfile(basename,firstfilenum)

c     print date & time
      call idate(today)
      call itime(now)
      print *
      print *,'run particles hour, min, sec =',now

c     loop through and advance particles & output data
      call particle_loop(firstfilenum,t2)

c     print date & time
      call idate(today)
      call itime(now)
      print *
      print *,'finished running particles hour, min, sec =',now

c     close for multiple output files
!      call closefile()

c     loop over input field files
      print *
      print *,'*** IN PROGRAM MAIN ***'
      print *,'first,last,step=',(firstfilenum+nt*filenumstep),
     &lastfilenum,(nt-2)*filenumstep

      do filenum=(firstfilenum+filenumstep),
     &lastfilenum,filenumstep
         print *
         print *,'*** IN PROGRAM MAIN ***'
         print *,'filenum =',filenum
	 
c        read fields
c        put last 2 time steps into first 2 time grid posiitons
         call last2_2first2()
c        read in 3rd time step through nt steps
         call load_field(basename,filenum,1,(nt-2),3)
         call scale_fields(3,nt)
         call calc_derivs(3,nt)
         call calc_ddt(2,nt)

         call subtrct_tzero(3,nt)
c        zero each time grid
         tgrmin=tgr(1)
         call zero_tgrid(tgrmin)
         tgrmax=tgr(nt)

c        set upper limit on time
         t2=amin1((tmax-tzero),tgrmax)

c        calc. field dependent particle data
!         call dist_particles2(t2)
c        field dataset i/o
         call field_io(basename,filenum,firstfilenum,timesteps)
c        to open multiple output files
         call openfile(basename,filenum)
c        advance particles & output data
         call particle_loop(firstfilenum,t2)
c        close for multiple output files
         call closefile()
      enddo

c     finish up with i/o routines
      call io_close(basename,firstfilenum)
      print *

      end
