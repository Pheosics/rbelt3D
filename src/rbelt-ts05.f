c     rbelt-ts05.f - analytic TS05 field routines

************************************************************************

      subroutine set_sys()
      include 'rbelt-grid.inc'
      sys=2
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
*      include 'rbelt-fields.inc'
      include 'rbelt-grid.inc'
*      include 'rbelt-bounds.inc'
*      include 'rbelt-const.inc'
      include 'rbelt-geopack_08.inc'
      integer filenum,firststep,laststep,gridstep,index
      character*(*) filename
      print *,'in hdfread: using analytic field model'
      tgr(1)=0.0
      
c     get TS05 input parameters (here for time independent fields)
      if (tsy.eqv..true.) then
         call get_ts05_params(tgr(1))
         if (abs(tilt0-PSI).gt.0.01) then
            print *,'WARNING! TS05 input file and RECALC'
            print *,'tilt angles do not agree'
         endif
      endif

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

c     must call Geopack RECALC before calling this routine
c     THIS ROUTINE AND LOAD FIELD ABOVE NEED WORK
c     (need option for time-dependent vs. time-independent fields)

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      include 'rbelt-geopack_08.inc'
      include 'rbelt-ut.inc'
      real y(6),t,r
      REAL BXIN,BXEX,BYIN,BYEX,BZIN,BZEX

*      print *,'   in get_fields (TS05): y,t=',y,t
c     first check to see if particle has steped out of bounds
      if (status.gt.0) return
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      r=sqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))

*      print *,'r=',r

      if (r.lt.rmin) then
         status=1
         return
      elseif (r.gt.rmax) then
         status=2
         return
      endif

c     get TS05 input parameters (here for time dependent fields)
      if (tsy.eqv..true.) then
         sfr=(t+tzero)/tfactor
         call get_ts05_params(sfr)
         if (abs(tilt0-PSI).gt.0.01) then
            print *,'WARNING! input file and calculated'
            print *,'tilt angles do not agree'
         endif
      endif

      CALL IGRF_GSM(y(1),y(2),y(3),BXIN,BYIN,BZIN)
      if (tsy.eqv..true.) then
         CALL T04_s(IOPT,PARMOD,PSI,y(1),y(2),y(3),BXEX,BYEX,BZEX)
      endif
      bx=(BXIN+BXEX)*fnorm
      by=(BYIN+BYEX)*fnorm
      bz=(BZIN+BZEX)*fnorm
      b=sqrt(bx*bx+by*by+bz*bz)

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


