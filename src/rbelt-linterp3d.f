c     rbelt-linterp3d.f - tri-linear interpolation for static fields

************************************************************************

      subroutine get_fields(y,t)

c tri-linear interpolation of gridded E & B fields to particle position.
c from data on a uniform cartisian space grid. note that the 
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
      integer ix,iy,iz
      integer jjji,ijji,jiji,iiji,jjii,ijii,jiii,iiii
      real w1,w1m,w2,w2m,w3,w3m,w1m2m
      real w12m,w12,w1m2,w1m2m3m,w12m3m,w123m
      real w1m23m,w1m2m3,w123,w1m23,w12m3
      real ww01,ww02,ww03,ww04,ww05,ww06,ww07,ww08
      real w2m3m,w23m,w2m3,w23,w1m3m,w13m,w1m3,w13
      real r,r2,bfac1,bfac2

*      print *
*      print *,'in SUBROUTINE get_fields'
*      print *,'********** y,t=',y,t/tfactor
*      print *,'rmin,rmax,tgrmax=',rmin,rmax,tgrmax

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
      elseif (r.gt.rmax) then
         status=2
         return
      elseif ((y(3).lt.zmin).or.(y(3).ge.zmax)) then
         status=2
         return
      endif
*c     then check to see if it's time to halt the integration
*c     (this should be moved out of linterp and into time loop)
*      if (t.ge.tstop) then
*         if (status.eq.0) status=-1
*      endif

c     uncomment to debug
c     make sure that (y,t) lies in the grid
*      if ((y(1).lt.xgr(1)).or.(y(1).ge.xgr(nx)).or.
*     &   (y(2).lt.ygr(1)).or.(y(2).ge.ygr(ny)).or.
*     &   (y(3).lt.zgr(1)).or.(y(3).ge.zgr(nz))) then
*         status=1
*         print *
*         print *,'in linterp3: outside grid'
*         print *,'t,y=',t/tfactor,y
*         if (y(1).lt.xgr(1)) print*,'y(1) < xgr(1)'
*         if (y(1).ge.xgr(nx)) print*,'y(1) >= xgr(nx)'
*         if (y(2).lt.ygr(1)) print*,'y(2) < ygr(1)'
*         if (y(2).ge.ygr(ny)) print*,'y(2) >= ygr(ny)'
*         if (y(3).lt.zgr(1)) print*,'y(3) < zgr(1)'
*         if (y(3).ge.zgr(nz)) print*,'y(3) >= zgr(nz)'
*         stop
*      endif

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

      ijji= ix + iy*nx + iz*nxy
      jjji= ijji + 1
      jiji= jjji - nx
      iiji= ijji - nx
      jjii= jjji - nxy
      ijii= ijji - nxy
      jiii= jiji - nxy
      iiii= iiji - nxy

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3

      w1m2m=w1m*w2m
      w12m=w1*w2m
      w12=w1*w2
      w1m2=w1m*w2
	  
      ww01=w1m2m*w3m
      ww02=w12m*w3m
      ww03=w12*w3m
      ww04=w1m2*w3m
      ww05=w1m2m*w3
      ww06=w12m*w3
      ww07=w12*w3
      ww08=w1m2*w3

      r2=r*r
      bfac1=3.*b0/r2/r2/r

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     1   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     2   +bxdv(jjji)*ww07+bxdv(ijji)*ww08-bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     1   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     2   +bydv(jjji)*ww07+bydv(ijji)*ww08-bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     1   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     2   +bzdv(jjji)*ww07+bzdv(ijji)*ww08-bfac1*y(3)*y(3)+b0/r2/r

      b=sqrt(bx**2.+by**2.+bz**2.)

      return
      end

************************************************************************

      subroutine get_fields2(y,t)

c tri-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space grid. note that the 
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
      integer ix,iy,iz
      integer jjji,ijji,jiji,iiji,jjii,ijii,jiii,iiii
      real w1,w1m,w2,w2m,w3,w3m,w1m2m
      real w12m,w12,w1m2,w1m2m3m,w12m3m,w123m
      real w1m23m,w1m2m3,w123,w1m23,w12m3
      real ww01,ww02,ww03,ww04,ww05,ww06,ww07,ww08
      real w2m3m,w23m,w2m3,w23,w1m3m,w13m,w1m3,w13
      real r,r2,bfac1,bfac2

*      print *
*      print *,'in SUBROUTINE get_fields2'
*      print *,'y,t=',y,t/tfactor

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
      elseif (r.gt.rmax) then
         status=2
         return
      endif

c     uncomment to debug
c     make sure that (y,t) lies in the grid
*      if ((y(1).lt.xgr(1)).or.(y(1).ge.xgr(nx)).or.
*     &   (y(2).lt.ygr(1)).or.(y(2).ge.ygr(ny)).or.
*     &   (y(3).lt.zgr(1)).or.(y(3).ge.zgr(nz))) then
*         status=1
*         print *
*         print *,'in linterp3: outside grid'
*         print *,'t,y=',t/tfactor,y
*         if (y(1).lt.xgr(1)) print*,'y(1) < xgr(1)'
*         if (y(1).ge.xgr(nx)) print*,'y(1) >= xgr(nx)'
*         if (y(2).lt.ygr(1)) print*,'y(2) < ygr(1)'
*         if (y(2).ge.ygr(ny)) print*,'y(2) >= ygr(ny)'
*         if (y(3).lt.zgr(1)) print*,'y(3) < zgr(1)'
*         if (y(3).ge.zgr(nz)) print*,'y(3) >= zgr(nz)'
*         stop
*      endif

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

      ijji= ix + iy*nx + iz*nxy
      jjji= ijji + 1
      jiji= jjji - nx
      iiji= ijji - nx
      jjii= jjji - nxy
      ijii= ijji - nxy
      jiii= jiji - nxy
      iiii= iiji - nxy

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3

      w1m2m=w1m*w2m
      w12m=w1*w2m
      w12=w1*w2
      w1m2=w1m*w2

      ww01=w1m2m*w3m
      ww02=w12m*w3m
      ww03=w12*w3m
      ww04=w1m2*w3m
      ww05=w1m2m*w3
      ww06=w12m*w3
      ww07=w12*w3
      ww08=w1m2*w3

      r2=r*r
      bfac1=3.*b0/r2/r2/r
      bfac2=5.*bfac1/r2

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     1   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     2   +bxdv(jjji)*ww07+bxdv(ijji)*ww08-bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     1   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     2   +bydv(jjji)*ww07+bydv(ijji)*ww08-bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     1   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     2   +bzdv(jjji)*ww07+bzdv(ijji)*ww08-bfac1*y(3)*y(3)+b0/r2/r

      b=sqrt(bx**2.+by**2.+bz**2.)

      dbxdx=dbxdxv(iiii)*ww01+dbxdxv(jiii)*ww02+dbxdxv(jjii)*ww03
     1   +dbxdxv(ijii)*ww04+dbxdxv(iiji)*ww05+dbxdxv(jiji)*ww06
     2   +dbxdxv(jjji)*ww07+dbxdxv(ijji)*ww08
     8   -bfac1*y(3)+bfac2*y(1)*y(1)*y(3)

      dbydx=dbydxv(iiii)*ww01+dbydxv(jiii)*ww02+dbydxv(jjii)*ww03
     1   +dbydxv(ijii)*ww04+dbydxv(iiji)*ww05+dbydxv(jiji)*ww06
     2   +dbydxv(jjji)*ww07+dbydxv(ijji)*ww08
     8   +bfac2*y(2)*y(3)*y(1)

      dbzdx=dbzdxv(iiii)*ww01+dbzdxv(jiii)*ww02+dbzdxv(jjii)*ww03
     1   +dbzdxv(ijii)*ww04+dbzdxv(iiji)*ww05+dbzdxv(jiji)*ww06
     2   +dbzdxv(jjji)*ww07+dbzdxv(ijji)*ww08
     8   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)

      dbxdy=dbxdyv(iiii)*ww01+dbxdyv(jiii)*ww02+dbxdyv(jjii)*ww03
     1   +dbxdyv(ijii)*ww04+dbxdyv(iiji)*ww05+dbxdyv(jiji)*ww06
     2   +dbxdyv(jjji)*ww07+dbxdyv(ijji)*ww08
     8   +bfac2*y(1)*y(2)*y(3)

      dbydy=dbydyv(iiii)*ww01+dbydyv(jiii)*ww02+dbydyv(jjii)*ww03
     1   +dbydyv(ijii)*ww04+dbydyv(iiji)*ww05+dbydyv(jiji)*ww06
     2   +dbydyv(jjji)*ww07+dbydyv(ijji)*ww08
     8   -bfac1*y(3)+bfac2*y(2)*y(2)*y(3)

      dbzdy=dbzdyv(iiii)*ww01+dbzdyv(jiii)*ww02+dbzdyv(jjii)*ww03
     1   +dbzdyv(ijii)*ww04+dbzdyv(iiji)*ww05+dbzdyv(jiji)*ww06
     2   +dbzdyv(jjji)*ww07+dbzdyv(ijji)*ww08
     8   -bfac1*y(2)+bfac2*y(2)*y(3)*y(3)

      dbxdz=dbxdzv(iiii)*ww01+dbxdzv(jiii)*ww02+dbxdzv(jjii)*ww03
     1   +dbxdzv(ijii)*ww04+dbxdzv(iiji)*ww05+dbxdzv(jiji)*ww06
     2   +dbxdzv(jjji)*ww07+dbxdzv(ijji)*ww08
     8   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)

      dbydz=dbydzv(iiii)*ww01+dbydzv(jiii)*ww02+dbydzv(jjii)*ww03
     1   +dbydzv(ijii)*ww04+dbydzv(iiji)*ww05+dbydzv(jiji)*ww06
     2   +dbydzv(jjji)*ww07+dbydzv(ijji)*ww08
     8   -bfac1*y(2)+bfac2*y(3)*y(3)*y(2)

      dbzdz=dbzdzv(iiii)*ww01+dbzdzv(jiii)*ww02+dbzdzv(jjii)*ww03
     1   +dbzdzv(ijii)*ww04+dbzdzv(iiji)*ww05+dbzdzv(jiji)*ww06
     2   +dbzdzv(jjji)*ww07+dbzdzv(ijji)*ww08
     8   -3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)

      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b

      return
      end


************************************************************************

c     pure dipole B field
      subroutine get_dipole(y,t)
*      subroutine get_fields(y,t)
      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc' 
      include 'rbelt-fields.inc'
      include 'rbelt-status.inc'
      real y(6)
      real t
      real r,r2,r5,r7,bfac1,bfac2,b0t,b0_tdep
      real cphi,sphi

*      print *,'   in get_fields (dipole): y,t=',y,t
c      put b0_tdep in input file
      b0_tdep=0.
      b0t=b0+b0_tdep*t

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

      bfac1=3.*b0t/r2/r2/r
      bx=-bfac1*y(1)*y(3)
      by=-bfac1*y(2)*y(3)
      bz=-bfac1*y(3)*y(3)+b0t/r2/r
      b=sqrt(bx*bx+by*by+bz*bz)

*      cphi = y(1)/sqrt(y(1)*y(1)+y(2)*y(2))
*      sphi = y(2)/sqrt(y(1)*y(1)+y(2)*y(2))
*      ex=-(-0.01*ffactor*vmsvcm)*sphi
*      ey=(-0.01*ffactor*vmsvcm)*cphi

      return
      end

************************************************************************

      subroutine get_fields_old(y,t)

c tri-linear interpolation of gridded E & B fields to particle position.
c from data on a uniform cartisian space grid. note that the 
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
      integer ix,iy,iz
      integer jjji,ijji,jiji,iiji,jjii,ijii,jiii,iiii
      real w1,w1m,w2,w2m,w3,w3m,w1m2m
      real w12m,w12,w1m2,w1m2m3m,w12m3m,w123m
      real w1m23m,w1m2m3,w123,w1m23,w12m3
      real ww01,ww02,ww03,ww04,ww05,ww06,ww07,ww08
      real w2m3m,w23m,w2m3,w23,w1m3m,w13m,w1m3,w13
      real r,r2,bfac1,bfac2

*      print *,'in SUBROUTINE get_fields'
*      print *
*      print *,'********** y,t=',y,t/tfactor
*      print *
*      print *,'rmin,rmax,tgrmax,tstop',rmin,rmax,tgrmax,tstop

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
      elseif (r.gt.rmax) then
         status=2
         return
      elseif ((y(3).lt.zmin).or.(y(3).ge.zmax)) then
         status=2
         return
      endif

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

      ijji= ix + iy*nx + iz*nxy
      jjji= ijji + 1
      jiji= jjji - nx
      iiji= ijji - nx
      jjii= jjji - nxy
      ijii= ijji - nxy
      jiii= jiji - nxy
      iiii= iiji - nxy

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3

      w1m2m=w1m*w2m
      w12m=w1*w2m
      w12=w1*w2
      w1m2=w1m*w2
	  
      w1m2m3m=w1m2m*w3m
      w12m3m=w12m*w3m
      w123m=w12*w3m
      w1m23m=w1m2*w3m

      w1m2m3=w1m2m*w3
      w12m3=w12m*w3
      w123=w12*w3
      w1m23=w1m2*w3

      ww01=w1m2m3m
      ww02=w12m3m
      ww03=w123m
      ww04=w1m23m
      ww05=w1m2m3
      ww06=w12m3
      ww07=w123
      ww08=w1m23

      r2=r*r
      bfac1=3.*b0/r2/r2/r

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     1   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     2   +bxdv(jjji)*ww07+bxdv(ijji)*ww08-bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     1   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     2   +bydv(jjji)*ww07+bydv(ijji)*ww08-bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     1   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     2   +bzdv(jjji)*ww07+bzdv(ijji)*ww08-bfac1*y(3)*y(3)+b0/r2/r

      b=sqrt(bx**2.+by**2.+bz**2.)

      return
      end

************************************************************************

      subroutine get_fields2_old(y,t)

c old: better conservation of energy, possibly because of exact agreement here
c between B function & derivatives, but scatters particles in pitch angle 
c because of discontinuous grad B

c tri-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space grid. note that the 
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
      integer ix,iy,iz
      integer jjji,ijji,jiji,iiji,jjii,ijii,jiii,iiii
      real w1,w1m,w2,w2m,w3,w3m,w1m2m
      real w12m,w12,w1m2,w1m2m3m,w12m3m,w123m
      real w1m23m,w1m2m3,w123,w1m23,w12m3
      real ww01,ww02,ww03,ww04,ww05,ww06,ww07,ww08
      real w2m3m,w23m,w2m3,w23,w1m3m,w13m,w1m3,w13
      real r,r2,bfac1,bfac2

*      print *,'in SUBROUTINE get_fields'
*      print *,'y,t=',y,t/tfactor

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
      elseif (r.gt.rmax) then
         status=2
         return
      endif

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

      ijji= ix + iy*nx + iz*nxy
      jjji= ijji + 1
      jiji= jjji - nx
      iiji= ijji - nx
      jjii= jjji - nxy
      ijii= ijji - nxy
      jiii= jiji - nxy
      iiii= iiji - nxy

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3

      w1m2m=w1m*w2m
      w12m=w1*w2m
      w12=w1*w2
      w1m2=w1m*w2
	  
      w1m2m3m=w1m2m*w3m
      w12m3m=w12m*w3m
      w123m=w12*w3m
      w1m23m=w1m2*w3m

      w1m2m3=w1m2m*w3
      w12m3=w12m*w3
      w123=w12*w3
      w1m23=w1m2*w3

      ww01=w1m2m3m
      ww02=w12m3m
      ww03=w123m
      ww04=w1m23m
      ww05=w1m2m3
      ww06=w12m3
      ww07=w123
      ww08=w1m23

      r2=r*r
      bfac1=3.*b0/r2/r2/r

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     1   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     2   +bxdv(jjji)*ww07+bxdv(ijji)*ww08-bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     1   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     2   +bydv(jjji)*ww07+bydv(ijji)*ww08-bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     1   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     2   +bzdv(jjji)*ww07+bzdv(ijji)*ww08-bfac1*y(3)*y(3)+b0/r2/r


*      ex=exdv(iiii)*ww01+exdv(jiii)*ww02+exdv(jjii)*ww03
*     1   +exdv(ijii)*ww04+exdv(iiji)*ww05+exdv(jiji)*ww06
*     2   +exdv(jjji)*ww07+exdv(ijji)*ww08
*      ey=eydv(iiii)*ww01+eydv(jiii)*ww02+eydv(jjii)*ww03
*     1   +eydv(ijii)*ww04+eydv(iiji)*ww05+eydv(jiji)*ww06
*     2   +eydv(jjji)*ww07+eydv(ijji)*ww08
*      ez=ezdv(iiii)*ww01+ezdv(jiii)*ww02+ezdv(jjii)*ww03
*     1   +ezdv(ijii)*ww04+ezdv(iiji)*ww05+ezdv(jiji)*ww06
*     2   +ezdv(jjji)*ww07+ezdv(ijji)*ww08

      b=sqrt(bx**2.+by**2.+bz**2.)

*      print *,'in linterp get_fields2: x,b=',
*     &y(1),y(2),y(3),b/ffactor/ntg

      bfac2=5.*bfac1/r2

      w2m3m=w2m*w3m
      w23m=w2*w3m
      w2m3=w2m*w3
      w23=w2*w3

      w1m3m=w1m*w3m
      w13m=w1*w3m
      w1m3=w1m*w3
      w13=w1*w3

      dbxdx=((bxdv(jiii)-bxdv(iiii))*w2m3m
     1   +(bxdv(jjii)-bxdv(ijii))*w23m
     2   +(bxdv(jiji)-bxdv(iiji))*w2m3
     3   +(bxdv(jjji)-bxdv(ijji))*w23)/dx
     8   -bfac1*y(3)+bfac2*y(1)*y(1)*y(3)

      dbxdy=((bxdv(ijii)-bxdv(iiii))*w1m3m
     1   +(bxdv(jjii)-bxdv(jiii))*w13m
     2   +(bxdv(ijji)-bxdv(iiji))*w1m3
     3   +(bxdv(jjji)-bxdv(jiji))*w13)/dy
     8   +bfac2*y(1)*y(2)*y(3)

      dbxdz=((bxdv(iiji)-bxdv(iiii))*w1m2m
     1   +(bxdv(jiji)-bxdv(jiii))*w12m
     2   +(bxdv(ijji)-bxdv(ijii))*w1m2
     3   +(bxdv(jjji)-bxdv(jjii))*w12)/dz
     8   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbydx=((bydv(jiii)-bydv(iiii))*w2m3m
     1   +(bydv(jjii)-bydv(ijii))*w23m
     2   +(bydv(jiji)-bydv(iiji))*w2m3
     3   +(bydv(jjji)-bydv(ijji))*w23)/dx
     8   +bfac2*y(2)*y(3)*y(1)

      dbydy=((bydv(ijii)-bydv(iiii))*w1m3m
     1   +(bydv(jjii)-bydv(jiii))*w13m
     2   +(bydv(ijji)-bydv(iiji))*w1m3
     3   +(bydv(jjji)-bydv(jiji))*w13)/dy
     8   -bfac1*y(3)+bfac2*y(2)*y(2)*y(3)

      dbydz=((bydv(iiji)-bydv(iiii))*w1m2m
     1   +(bydv(jiji)-bydv(jiii))*w12m
     2   +(bydv(ijji)-bydv(ijii))*w1m2
     3   +(bydv(jjji)-bydv(jjii))*w12)/dz
     8   -bfac1*y(2)+bfac2*y(3)*y(3)*y(2)

      dbzdx=((bzdv(jiii)-bzdv(iiii))*w2m3m
     1   +(bzdv(jjii)-bzdv(ijii))*w23m
     2   +(bzdv(jiji)-bzdv(iiji))*w2m3
     3   +(bzdv(jjji)-bzdv(ijji))*w23)/dx
     8   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbzdy=((bzdv(ijii)-bzdv(iiii))*w1m3m
     1   +(bzdv(jjii)-bzdv(jiii))*w13m
     2   +(bzdv(ijji)-bzdv(iiji))*w1m3
     3   +(bzdv(jjji)-bzdv(jiji))*w13)/dy
     8   -bfac1*y(2)+bfac2*y(2)*y(3)*y(3)

      dbzdz=((bzdv(iiji)-bzdv(iiii))*w1m2m
     1   +(bzdv(jiji)-bzdv(jiii))*w12m
     2   +(bzdv(ijji)-bzdv(ijii))*w1m2
     3   +(bzdv(jjji)-bzdv(jjii))*w12)/dz
     8   -3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)

      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b

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

