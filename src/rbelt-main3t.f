c rbelt-main3t.f
c main program for computing particle trajectories in static field
c snapshot in rbelt time grid position 3, interpolated from fields
c in rbelt time grid positions 1 and 2.



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
*      include 'rbelt-y0.inc'
      character*80 basename
      integer filenum,firstfilenum,lastfilenum,filenumstep,timesteps
      integer*4 today(3), now(3)
      NAMELIST /rbelt/ basename,firstfilenum,lastfilenum,filenumstep
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
      if (nstep.ne.3) then
         print *,'nstep=',nstep,'   (should be nstep=3)'
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
      filename2=fieldfile(basename,firstfilenum)
      filename3=fieldfile(basename,firstfilenum+filenumstep)
c     (3rd input firststep=2 is unconventional use of this input
c     to flag load_field that we are using main3t version of code --
c     3rd and 4th inputs otherwise not used in this version.)
      call load_field(filename2,firstfilenum,2,1,2)
      call load_field(filename3,firstfilenum+filenumstep,2,1,3)
      call scale_fields(2,3)

c     time and space grid & boundaries initialization
c     grid_init() & bounds_init() get called in linterp321(t)
c     after we have tgr(1), tzero1

c     initialize field quantities to zero
      call init_fields

c     calculate grad B (must calc. dx,dy,dz in grid_init first)
c      call calc_derivs(1,nstep)

c     zero time step in linterp321(t)

c     initialize integrators
      call lrntz_params()
      call gc_params()
      call rkckparm

c     this is not used in this version
      t2 = 0.0

c     distribute particles
      call dist_particles(t2)
c     calc. field dependent particle data
      call dist_particles2(t2)

c     initialize i/o routines
      call io_init(basename,firstfilenum)
      call openfile(basename,firstfilenum)

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
      call closefile()

c     loop over input field files
*      print *
*      print *,'*** IN PROGRAM MAIN ***'
*      print *,'first,last,step=',firstfilenum,
*     &   lastfilenum,filenumstep
      do filenum=firstfilenum+2*filenumstep,lastfilenum,
     &   filenumstep
         print *
         print *,'*** IN PROGRAM MAIN ***'
         print *,'filenum =',filenum
c        put last time step (here 3rd) into second time grid posiiton
         call last2second()

c        read in 3rd time step
         filename3=fieldfile(basename,filenum)
         call load_field(filename3,filenum,2,1,3)
         call scale_fields(3,3)
c        calc. field dependent particle data
         call dist_particles2(t2)
c        to open multiple output files
         call openfile(basename,firstfilenum)
c        advance particles & output data
         call particle_loop(filenum,t2)
c        close for multiple output files
         call closefile()

      enddo

c     finish up with i/o routines
      call io_close(basename,firstfilenum)
      print *

      end




