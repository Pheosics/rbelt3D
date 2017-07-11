c     rbelt-zero_fields.f

************************************************************************

      subroutine set_sys()
      include 'rbelt-grid.inc'
      sys=1
      return
      end

************************************************************************

      subroutine load_ut(basename,firstfilenum,lastfilenum,filenumstep)

c     read ut file
      implicit none
      include 'rbelt-ut.inc'
      character*(*) basename
      integer firstfilenum,lastfilenum,filenumstep
      NAMELIST /date/ year0,doy0,hour0,min0,sec0,dsec,table

      print *
      print *,'*** in subroutine read_ut ***'
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

      print *,'in hdfread: using analytic field model'

      tgr(1)=0.0

      return
      end

************************************************************************

      subroutine get_fields(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the 
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields here, using the value
c of B0 specified in rbelt-const.inc.

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'

      real y(3),t
      integer itu,itm,it,ix,iy,iz,dt
      integer jjjj,ijjj,jijj,iijj,jjij,ijij,jiij,iiij
      integer jjji,ijji,jiji,iiji,jjii,ijii,jiii,iiii
      real w1,w1m,w2,w2m,w3,w3m,w4,w4m,w1m2m
      real w12m,w12,w1m2,w1m2m3m,w12m3m,w123m
      real w1m23m,w1m2m3,w123,w1m23,w12m3
      real ww01,ww02,ww03,ww04,ww05,ww06,ww07,ww08
      real ww09,ww10,ww11,ww12,ww13,ww14,ww15,ww16
      real w2m3m4m,w23m4m,w2m34m,w234m,w2m3m4,w23m4
      real w1m3m4m,w13m4m,w13m4,w1m24m,w124m,w12m4,w1m2m4
      real w2m34,w234,w1m34m,w134m,w1m3m4,w134
      real w1m34,w1m2m4m,w12m4m,w1m24,w124
      real r,r2,bfac1
      
*      print *
*      print *,'*** in linterp4d ***'
*      print *,'t,y=',t/tfactor,y
*      print *,'tgrmax,tstop=',tgrmax/tfactor,tstop/tfactor

c     first check to see if particle has steped out of bounds
      if (status.gt.0) return
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      r=sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
      if (r.lt.rmin) then
         status=1
         return
*      elseif (r.gt.rmax) then
*         status=2
*         return
      elseif ((y(1).lt.xmin).or.(y(1).ge.xmax)) then
         status=2
         return
      elseif ((y(2).lt.ymin).or.(y(2).ge.ymax)) then
         status=2
         return
      elseif ((y(3).lt.zmin).or.(y(3).ge.zmax)) then
         status=2
         return
      endif

      return
      end

************************************************************************

      subroutine get_fields2(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the 
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields, using the value
c of B0 specified in rbelt-const.inc.
c this version also calculates components of grad B and grad |B| 
c (needed for the guiding center routines) using a first order 
c approximation for purtabation fileds.

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      include 'rbelt-const.inc'

      real y(3),t
      integer itu,itm,it,ix,iy,iz,dt
      integer jjjj,ijjj,jijj,iijj,jjij,ijij,jiij,iiij
      integer jjji,ijji,jiji,iiji,jjii,ijii,jiii,iiii
      real w1,w1m,w2,w2m,w3,w3m,w4,w4m,w1m2m
      real w12m,w12,w1m2,w1m2m3m,w12m3m,w123m
      real w1m23m,w1m2m3,w123,w1m23,w12m3
      real ww01,ww02,ww03,ww04,ww05,ww06,ww07,ww08
      real ww09,ww10,ww11,ww12,ww13,ww14,ww15,ww16
      real w2m3m4m,w23m4m,w2m34m,w234m,w2m3m4,w23m4
      real w1m3m4m,w13m4m,w13m4,w1m24m,w124m,w12m4,w1m2m4
      real w2m34,w234,w1m34m,w134m,w1m3m4,w134
      real w1m34,w1m2m4m,w12m4m,w1m24,w124
      real r,r2,bfac1

*      print *
*      print *,'*** in linterp4d ***'
*      print *,'status,t,y=',status,t/tfactor,y
*      print *,'tgrmax,tstop=',tgrmax/tfactor,tstop/tfactor

c     first check to see if particle has steped out of bounds
      if (status.gt.0) return
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      r=sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
      if (r.lt.rmin) then
         status=1
         return
*      elseif (r.gt.rmax) then
*         status=2
*         return
      elseif ((y(1).lt.xmin).or.(y(1).ge.xmax)) then
         status=2
         return
      elseif ((y(2).lt.ymin).or.(y(2).ge.ymax)) then
         status=2
         return
      elseif ((y(3).lt.zmin).or.(y(3).ge.zmax)) then
         status=2
         return
      endif

      return
      end

