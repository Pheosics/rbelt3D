c     rbelt-dipole.f - pure dipole field routines

************************************************************************

      subroutine load_ut(basename,firstfilenum,lastfilenum,filenumstep)

c     read ut
      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-ut.inc'
      include 'rbelt-dst.inc'

      character*(*) basename
      character*80 filename
      integer i,j,mode,firstfilenum,lastfilenum,filenumstep,sec_frm_ref
      integer year,mon,day,doy,hour,min,sec,status

      NAMELIST /date/ year0,doy0,hour0,min0,sec0,dsec,table
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,date)
      print *
      print *,'*** in subroutine load_ut ***'
      print *,'year0,doy0,hour0,min0,sec0,table=',
     &year0,doy0,hour0,min0,sec0,table

c     if table.eqv..true. we do not use date from input file.
c     first timestamp in table is used as UT reference time.
      if (table.eqv..true.) then
c        get dst data
c        open & read in dst table
         open (UNIT=10,file='dst.dat',STATUS='OLD')
         read (10,*) year,mon
         print *,'year,mon=',year,mon
c        read number of lines/days in file
         read (10,*) numdays
         print *,'opened & reading dst.dat: numdays=',numdays
         if (numdays.gt.maxdays) then
            print *,'numdays.gt.maxdays (in dst.inc file)'
            stop
         endif

c        read from dst data file (MODIFY TO AGREE WITH FILE FORMAT AS NEEDED)
         read (10,*) day,dummy

c        convert month & day to doy
         status=0
         mode=-1
         call doy_mmdd (mode,year,doy,mon,day,status)
         if (status.ne.0) then
            print *,'status.ne.0 in subroutine doy_mmdd'
            stop
         endif

c        Kyoto DST format
c        first DST has hour=1 timestamp, 
c        but we may want to change to 1/2, 
c        since it is an avg. over previous hour.
         hour=1
         min=0
         sec=0

c        set ref. time to initial dst table time
         year0=year
         doy0=doy
         hour0=hour
         min0=min
         sec0=sec

c        find seconds from start of file
c        (sfr is not converted to normalized units)
         sfr(1)=sec_frm_ref(year,doy,hour,min,sec)
         dst(1)=dummy(1)

         do j=2,24
            hour=j
            sfr(j)=sec_frm_ref(year,doy,hour,min,sec)
            dst(j)=dummy(j)
         enddo

         do i=2,numdays

c           read from dst data file (MODIFY TO AGREE WITH FILE FORMAT AS NEEDED)
            read (10,*) day,dummy

c           convert month & day to doy
            status=0
            mode=-1
            call doy_mmdd (mode,year,doy,mon,day,status)
            if (status.ne.0) then
               print *,'status.ne.0 in subroutine doy_mmdd'
               stop
            endif

            do j=1,24
               hour=j
               sfr((i-1)*24+j)=sec_frm_ref(year,doy,hour,min,sec)*tfactor
               dst((i-1)*24+j)=dummy(j)*ffactor*ntg
            enddo

         enddo
         close(10)

      endif

*      do i=1,numdays
*         do j=1,24
*            print *,'hfr,dst=',sfr((i-1)*24+j)/3600.,dst((i-1)*24+j)
*         enddo
*      enddo

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
      print *,'in load_field: using analytic field model'

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
      include 'rbelt-dst.inc'
      real y(6)
      real t
      real r,r2,r5,r7,bfac1,bfac2
      real cphi,sphi

*      print *,'   in get_fields (dipole): y,t=',y,t

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

*      cphi = y(1)/sqrt(y(1)*y(1)+y(2)*y(2))
*      sphi = y(2)/sqrt(y(1)*y(1)+y(2)*y(2))
*      ex=-(-0.01*ffactor*vmsvcm)*sphi
*      ey=(-0.01*ffactor*vmsvcm)*cphi

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

