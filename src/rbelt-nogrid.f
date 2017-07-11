c     rbelt-dipole.f - pure dipole field routines

c     in the current version, no spatial grid means that we read the UT 
c     ref. date/time from the input file; whereas if we use rbelt-grid.f
c     then  UT ref. date/time is read from the UT file.

***********************************************************************

      subroutine scale_fields(start_step,last_step)
      implicit none	
      integer start_step,last_step
      return
      end

***********************************************************************

      subroutine zero_tgrid(time)
      implicit none	
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      real*8 time

      tzero = tzero + time
      print *,'i; seconds from UT reference time; sec. from start;
     & sec. from start of current time grid'
      tgr(1)=tgr(1)-time
      print *,1,real((tgr(1)+tzero)/tfactor),
     &real((tgr(1)+tzero-tzero1)/tfactor),real(tgr(1)/tfactor)
      time = 0.0

      return
      end

***********************************************************************

      subroutine bounds_init()

c simulation region boundaries initialization

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      NAMELIST /bounds/  rmin,rmax,xmin,xmax,ymin,ymax,zmin,zmax,tmax

c     read in and normalize input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,bounds)
      CLOSE (81)
      tmax = tmax*tfactor+tzero1

      return
      end
      
***********************************************************************

      subroutine grid_init()

c time and space grid initialization

      implicit none
      include 'rbelt-grid.inc'

*      print *
*      print *,'*** in subroutine grid_init ***'

      tzero=0.0
      tzero1=tgr(1)
      call set_sys()

      return
      end

***********************************************************************

      subroutine calc_derivs(start_step,last_step)
      implicit none
      integer start_step,last_step
      return
      end

***********************************************************************

      subroutine calc_ddt(start_step,last_step)
      implicit none
      integer start_step,last_step
      return
      end

***********************************************************************

      subroutine calc_e(start_step,last_step)
      implicit none
      integer start_step,last_step
      return
      end

