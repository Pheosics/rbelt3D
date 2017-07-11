c     rbelt-ts05.f - analytic TS05 field routines

************************************************************************

      subroutine set_sys()
      include 'rbelt-grid.inc'
      sys=2
      return
      end

************************************************************************

      subroutine load_ut(basename,firstfilenum,lastfilenum,filenumstep)
c     read ut
      implicit none
      integer firstfilenum,lastfilenum,filenumstep
      character*(*) basename
      include 'rbelt-ut.inc'
      NAMELIST /date/ year0,doy0,hour0,min0,sec0,dsec,table
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,date)
      print *
      print *,'*** in subroutine load_ut ***'
      print *,'year0,doy0,hour0,min0,sec0=',year0,doy0,hour0,min0,sec0
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

c     initialize fields

      subroutine init_fields()

      implicit none
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'

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

c     calc. field normalization factor
      fnorm=ffactor*ntg

      return
      end


************************************************************************

c     TS05 B field
      subroutine get_fields(y,t)

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      include 'rbelt-geopak_08.inc'
      include 'rbelt-ut.inc'
      real y(6),t,r,r2,bfac1
      real THETA,PHI,BR,BTHETA,BPHI
      real XMAG,YMAG,ZMAG,XGEO,YGEO,ZGEO
      real BXGEO,BYGEO,BZGEO,BXMAG,BYMAG,BZMAG

*      print *,'   in get_fields (TS05): y,t=',y,t
c     first check to see if particle has steped out of bounds
      if (status.gt.0) return
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      r2=y(1)*y(1)+y(2)*y(2)+y(3)*y(3)
      r=sqrt(r2)
      if (r.lt.rmin) then
         status=1
         return
      elseif (r.gt.rmax) then
         status=2
         return
      endif

c      IGRF in GEO coordinates
      call SPHCAR_08 (R,THETA,PHI,y(1),y(2),y(3),-1)
      CALL IGRF_GEO_08 (R,THETA,PHI,BR,BTHETA,BPHI)
      call BSPCAR_08 (THETA,PHI,BR,BTHETA,BPHI,BXGEO,BYGEO,BZGEO)
      bx=BXGEO*fnorm
      by=BYGEO*fnorm
      bz=BZGEO*fnorm
      b=sqrt(bx*bx+by*by+bz*bz)

c      to put fields in MAG coordinate system
*      call GEOMAG_08 (XGEO,YGEO,ZGEO,y(1),y(2),y(3),-1)
*      call SPHCAR_08 (R,THETA,PHI,XGEO,YGEO,ZGEO,-1)
*      CALL IGRF_GEO_08 (R,THETA,PHI,BR,BTHETA,BPHI)
*      call BSPCAR_08 (THETA,PHI,BR,BTHETA,BPHI,BXGEO,BYGEO,BZGEO)
*      call GEOMAG_08 (BXGEO,BYGEO,BZGEO,BXMAG,BYMAG,BZMAG,1)
*      bx=BXMAG*fnorm
*      by=BYMAG*fnorm
*      bz=BZMAG*fnorm
*      b=sqrt(bx*bx+by*by+bz*bz)

c     pure dipole
*      bfac1=3.*b00/r2/r2/r
*      bx=-bfac1*y(1)*y(3)
*      by=-bfac1*y(2)*y(3)
*      bz=-bfac1*y(3)*y(3)+b00/r2/r
*      b=sqrt(bx*bx+by*by+bz*bz)

*      print *,'x,y,z,BTHETA,BZGEO,bz=',y(1),y(2),y(3),BTHETA,BZGEO,bz

      return
      end

************************************************************************

c     pure dipole B field
      subroutine get_fields2(y,t)
      implicit none
      real y(3),t
      print *,
     &'can not compute guiding center traj. in analytic IGRF field'
      stop
      return
      end


