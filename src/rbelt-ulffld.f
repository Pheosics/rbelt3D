c     rbelt-dipole.f - pure dipole field routines

************************************************************************

      subroutine load_ut(basename,firstfilenum,lastfilenum,filenumstep)

c     read ut file
      implicit none
      include 'rbelt-ut.inc'

      character*80 filename,utfile,basename
      integer lnblnk,firstfilenum,lastfilenum,filenumstep
      NAMELIST /date/ year0,doy0,hour0,min0,sec0,dsec,table

      print *
      print *,'*** in subroutine read_ut ***'

c     read date & time from rbelt input file
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,date)
      print *
      print *,'*** in subroutine load_ut ***'
      print *,'year0,doy0,hour0,min0,sec0,table=',
     &year0,doy0,hour0,min0,sec0,table

      return
      end

************************************************************************

      subroutine load_grid(filename)
      implicit none	
      character*(*) filename
      return
      end

************************************************************************

c     pure dipole B field dummy fileread

      subroutine load_field
     &(filename,filenum,firststep,laststep,gridstep)

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      integer filenum,firststep,laststep,gridstep,index
      character*(*) filename
      print *
      print *,'in hdfread: using analytic field model'

      tgr(1)=0.0

      return
      end

************************************************************************

c     pure dipole B field
      subroutine get_fields(y,t)

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      real y(6)
      real t
      real r,r2,r5,r7,bfac1,bfac2
      real cphi,sphi

      if (status.gt.0) return
      r2=y(1)*y(1)+y(2)*y(2)+y(3)*y(3)
      r=sqrt(r2)
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      if (r.lt.rmin) then     
         status=1
         return
      elseif (r.gt.rmax) then
         status=2
         return
      endif

      bfac1=3.*b0/r2/r2/r
      bx=-bfac1*y(1)*y(3)
      by=-bfac1*y(2)*y(3)
      bz=-bfac1*y(3)*y(3)+b0/r2/r

      b=sqrt(bx*bx+by*by+bz*bz)

      return
      end

************************************************************************

c     pure dipole B field
      subroutine get_fields2(y,t)

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      real y(6)
      real t
      real r,r2,r5,r7,bfac1,bfac2
      real cphi,sphi

      if (status.gt.0) return
      r2=y(1)*y(1)+y(2)*y(2)+y(3)*y(3)
      r=sqrt(r2)
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      if (r.lt.rmin) then
         status=1
         return
      elseif (r.gt.rmax) then
         status=2
         return
      endif
*      if (t.ge.tstop) then
*         status=-1
*      endif

      bfac1=3.*b0/r2/r2/r
      bx=-bfac1*y(1)*y(3)
      by=-bfac1*y(2)*y(3)
      bz=-bfac1*y(3)*y(3)+b0/r2/r
      b=sqrt(bx*bx+by*by+bz*bz)

      bfac2=5.*bfac1/r2
      dbxdx=-bfac1*y(3)+bfac2*y(1)*y(1)*y(3)
      dbxdy=bfac2*y(1)*y(2)*y(3)
      dbxdz=-bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
      dbydx=bfac2*y(2)*y(3)*y(1)
      dbydy=-bfac1*y(3)+bfac2*y(2)*y(2)*y(3)
      dbydz=-bfac1*y(2)+bfac2*y(3)*y(3)*y(2)
      dbzdx=-bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
      dbzdy=-bfac1*y(2)+bfac2*y(2)*y(3)*y(3)
      dbzdz=-3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)
      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b

      return
      end

************************************************************************

      subroutine set_sys()
      include 'rbelt-grid.inc'
      sys=1
      return
      end

************************************************************************

c     initialize fields

      subroutine init_fields()

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'

c     for non-zero efield
      noe=.false.

c     set e field zero here to save time in the interp(analytic) routine
      ex=0.
      ey=0.
      ez=0.

c     initialize B and derrivatives here in case we use zero_fields
      bx=0.0
      by=0.0  
      bz= charge_sign*1.e-8
*      bz= 1.
      b=sqrt(bx*bx+by*by+bz*bz)
      dbxdx=0.0
      dbxdy=0.0
      dbxdz=0.0 
      dbydx=0.0
      dbydy=0.0
      dbydz=0.0
      dbzdx=0.0 
      dbzdy=0.0
      dbzdz=0.0
      dbdx=0.0
      dbdy=0.0
      dbdz=0.0
      dbxdt=0
      dbydt=0
      dbzdt=0
      dbdt=0

      return
      end
