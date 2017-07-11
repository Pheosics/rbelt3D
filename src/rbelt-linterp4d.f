c     rbelt-linterp4d.f - quad-linear interpolation for time dependent fields

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
*      print *,'*** in subroutine get_fields ***'
*      print *,'t,y=',t/tfactor,y
*      print *,'tgrmax=',tgrmax/tfactor

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

c     uncomment to debug
c     make sure that (y,t) lies in the grid
*      if ((y(1).lt.xgr(1)).or.(y(1).ge.xgr(nx)).or.
*     &   (y(2).lt.ygr(1)).or.(y(2).ge.ygr(ny)).or.
*     &   (y(3).lt.zgr(1)).or.(y(3).ge.zgr(nz)).or.
*     &   (t.lt.tgr(1)).or.(t.ge.tgr(nt))) then
*         status=1
*         print *
**         print *,'in linterp4: outside grid'
*         print *,'t,y=',t/tfactor,y
*         if (y(1).lt.xgr(1)) print*,'y(1) < xgr(1)'
*         if (y(1).ge.xgr(nx)) print*,'y(1) >= xgr(nx)'
*         if (y(2).lt.ygr(1)) print*,'y(2) < ygr(1)'
*         if (y(2).ge.ygr(ny)) print*,'y(2) >= ygr(ny)'
*         if (y(3).lt.zgr(1)) print*,'y(3) < zgr(1)'
*         if (y(3).ge.zgr(nz)) print*,'y(3) >= zgr(nz)'
*         if (t.lt.tgr(1)) print*,'t < tgr(1)',t,tgr(1)
*         if (t.ge.tgr(nt)) print*,'t >= tgr(nt)'
*         return
*      endif

      itu=nt+1
      it=0
10    continue
      if ((itu-it).gt.1) then
        itm=(itu+it)/2
        if (t.ge.tgr(itm)) then
          it=itm
        else
          itu=itm
        endif
        goto 10
      endif
      dt=tgr(it+1)-tgr(it)

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

      jjjj= ix + iy*nx + iz*nxy + it*nxyz + 1
      ijjj= jjjj - 1
      jijj= jjjj - nx
      iijj= jijj - 1
      jjij= jjjj - nxy
      ijij= jjij - 1
      jiij= jijj - nxy
      iiij= jiij - 1
      jjji= jjjj - nxyz
      ijji= ijjj - nxyz
      jiji= jijj - nxyz
      iiji= iijj - nxyz
      jjii= jjij - nxyz
      ijii= ijij - nxyz
      jiii= jiij - nxyz
      iiii= iiij - nxyz

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3
      w4=abs(t-tgr(it))/dt
      w4m=1.-w4

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

      ww01=w1m2m3m*w4m
      ww02=w12m3m*w4m
      ww03=w123m*w4m
      ww04=w1m23m*w4m
      ww05=w1m2m3*w4m
      ww06=w12m3*w4m
      ww07=w123*w4m
      ww08=w1m23*w4m
      ww09=w1m2m3m*w4
      ww10=w12m3m*w4
      ww11=w123m*w4
      ww12=w1m23m*w4
      ww13=w1m2m3*w4
      ww14=w12m3*w4
      ww15=w123*w4
      ww16=w1m23*w4

      r2=r*r
      bfac1=3.*b0/r2/r2/r

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     &   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     &   +bxdv(jjji)*ww07+bxdv(ijji)*ww08+bxdv(iiij)*ww09
     &   +bxdv(jiij)*ww10+bxdv(jjij)*ww11+bxdv(ijij)*ww12
     &   +bxdv(iijj)*ww13+bxdv(jijj)*ww14+bxdv(jjjj)*ww15
     &   +bxdv(ijjj)*ww16 
     &   -bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     &   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     &   +bydv(jjji)*ww07+bydv(ijji)*ww08+bydv(iiij)*ww09
     &   +bydv(jiij)*ww10+bydv(jjij)*ww11+bydv(ijij)*ww12
     &   +bydv(iijj)*ww13+bydv(jijj)*ww14+bydv(jjjj)*ww15
     &   +bydv(ijjj)*ww16 
     &   -bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     &   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     &   +bzdv(jjji)*ww07+bzdv(ijji)*ww08+bzdv(iiij)*ww09
     &   +bzdv(jiij)*ww10+bzdv(jjij)*ww11+bzdv(ijij)*ww12
     &   +bzdv(iijj)*ww13+bzdv(jijj)*ww14+bzdv(jjjj)*ww15
     &   +bzdv(ijjj)*ww16 
     &   -bfac1*y(3)*y(3) + b0/r2/r

      if (noe.eqv..true.) then
         ex=0.
         ey=0.
         ez=0.
      else
         ex=exdv(iiii)*ww01+exdv(jiii)*ww02+exdv(jjii)*ww03
     &   +exdv(ijii)*ww04+exdv(iiji)*ww05+exdv(jiji)*ww06
     &   +exdv(jjji)*ww07+exdv(ijji)*ww08+exdv(iiij)*ww09
     &   +exdv(jiij)*ww10+exdv(jjij)*ww11+exdv(ijij)*ww12
     &   +exdv(iijj)*ww13+exdv(jijj)*ww14+exdv(jjjj)*ww15
     &   +exdv(ijjj)*ww16 

         ey=eydv(iiii)*ww01+eydv(jiii)*ww02+eydv(jjii)*ww03
     &   +eydv(ijii)*ww04+eydv(iiji)*ww05+eydv(jiji)*ww06
     &   +eydv(jjji)*ww07+eydv(ijji)*ww08+eydv(iiij)*ww09
     &   +eydv(jiij)*ww10+eydv(jjij)*ww11+eydv(ijij)*ww12
     &   +eydv(iijj)*ww13+eydv(jijj)*ww14+eydv(jjjj)*ww15
     &   +eydv(ijjj)*ww16 

         ez=ezdv(iiii)*ww01+ezdv(jiii)*ww02+ezdv(jjii)*ww03
     &   +ezdv(ijii)*ww04+ezdv(iiji)*ww05+ezdv(jiji)*ww06
     &   +ezdv(jjji)*ww07+ezdv(ijji)*ww08+ezdv(iiij)*ww09
     &   +ezdv(jiij)*ww10+ezdv(jjij)*ww11+ezdv(ijij)*ww12
     &   +ezdv(iijj)*ww13+ezdv(jijj)*ww14+ezdv(jjjj)*ww15
     &   +ezdv(ijjj)*ww16 
      endif

      b=sqrt(bx**2.+by**2.+bz**2.)

*      print *,'   b=',b/ffactor/ntg

      return
      end

**********************************************************************************

      subroutine get_fields2(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the 
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields, using the value
c of B0 specified in rbelt-const.inc.
c this version also calculates components of grad B and grad |B| 
c (needed for the guiding center routines) using a first order 
c approximation for purtabation fileds.

c interpolates pre-calculated B-field gradients from cell verticies
c B-field gradients are calculated in grid.f using centered differences, 
c which results in gradients that are not everywhere/exactly consistant with B-field

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
      real r,r2,bfac1,bfac2
      real bxm,bym,bzm,bxp,byp,bzp,bm,bp

*      print *
*      print *,'*** in subroutine get_fields2 ***'
*      print *,'status,t,y=',status,t/tfactor,y
*      print *,'tgrmax=',tgrmax/tfactor

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
!      elseif ((y(1).lt.xmin).or.(y(1).ge.xmax)) then
!         status=2
!         return
!      elseif ((y(2).lt.ymin).or.(y(2).ge.ymax)) then
!         status=2
!         return
!      elseif ((y(3).lt.zmin).or.(y(3).ge.zmax)) then
!         status=2
!         return
      endif

c     uncomment to debug
c     make sure that (y,t) lies in the grid
*      if ((y(1).lt.xgr(1)).or.(y(1).ge.xgr(nx)).or.
*     &   (y(2).lt.ygr(1)).or.(y(2).ge.ygr(ny)).or.
*     &   (y(3).lt.zgr(1)).or.(y(3).ge.zgr(nz)).or.
*     &   (t.lt.tgr(1)).or.(t.ge.tgr(nt))) then
*         status=1
*         print *
*         print *,'in linterp4: outside grid'
*         print *,'t,y=',t/tfactor,y
*         if (y(1).lt.xgr(1)) print*,'y(1) < xgr(1)'
*         if (y(1).ge.xgr(nx)) print*,'y(1) >= xgr(nx)'
*         if (y(2).lt.ygr(1)) print*,'y(2) < ygr(1)'
*         if (y(2).ge.ygr(ny)) print*,'y(2) >= ygr(ny)'
*         if (y(3).lt.zgr(1)) print*,'y(3) < zgr(1)'
*         if (y(3).ge.zgr(nz)) print*,'y(3) >= zgr(nz)'
*         if (t.lt.tgr(1)) print*,'t < tgr(1)',t,tgr(1)
*         if (t.ge.tgr(nt)) print*,'t >= tgr(nt)'
*         return
*      endif

      itu=nt+1
      it=0
10    continue
      if ((itu-it).gt.1) then
        itm=(itu+it)/2
        if (t.ge.tgr(itm)) then
          it=itm
        else
          itu=itm
        endif
        goto 10
      endif
      dt=tgr(it+1)-tgr(it)

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

*      if ((ix.lt.1).or.(ix.gt.nx).or.
*     &   (iy.lt.1).or.(iy.gt.ny).or.
*     &   (iz.lt.1).or.(iz.gt.nz).or.
*     &   (it.lt.1).or.(it.gt.nt)) then
*         print *
*         print *,'ix=',ix
*         print *,'iy=',iy
*         print *,'iz=',iz
*         print *,'it=',it
*         print *,'t,r=',t/tfactor,r
*         print *,'status=',status
*         print *,'tgrmax=',tgrmax/tfactor
**         status=2
*         stop
*         return
*      endif

      jjjj= ix + iy*nx + iz*nxy + it*nxyz + 1
      ijjj= jjjj - 1
      jijj= jjjj - nx
      iijj= jijj - 1
      jjij= jjjj - nxy
      ijij= jjij - 1
      jiij= jijj - nxy
      iiij= jiij - 1
      jjji= jjjj - nxyz
      ijji= ijjj - nxyz
      jiji= jijj - nxyz
      iiji= iijj - nxyz
      jjii= jjij - nxyz
      ijii= ijij - nxyz
      jiii= jiij - nxyz
      iiii= iiij - nxyz

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3
      w4=abs(t-tgr(it))/dt
      w4m=1.-w4

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

      ww01=w1m2m3m*w4m
      ww02=w12m3m*w4m
      ww03=w123m*w4m
      ww04=w1m23m*w4m
      ww05=w1m2m3*w4m
      ww06=w12m3*w4m
      ww07=w123*w4m
      ww08=w1m23*w4m

      ww09=w1m2m3m*w4
      ww10=w12m3m*w4
      ww11=w123m*w4
      ww12=w1m23m*w4
      ww13=w1m2m3*w4
      ww14=w12m3*w4
      ww15=w123*w4
      ww16=w1m23*w4

      r2=r*r
      bfac1=3.*b0/r2/r2/r
      bfac2=5.*bfac1/r2

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     &   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     &   +bxdv(jjji)*ww07+bxdv(ijji)*ww08+bxdv(iiij)*ww09
     &   +bxdv(jiij)*ww10+bxdv(jjij)*ww11+bxdv(ijij)*ww12
     &   +bxdv(iijj)*ww13+bxdv(jijj)*ww14+bxdv(jjjj)*ww15
     &   +bxdv(ijjj)*ww16 
     &   -bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     &   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     &   +bydv(jjji)*ww07+bydv(ijji)*ww08+bydv(iiij)*ww09
     &   +bydv(jiij)*ww10+bydv(jjij)*ww11+bydv(ijij)*ww12
     &   +bydv(iijj)*ww13+bydv(jijj)*ww14+bydv(jjjj)*ww15
     &   +bydv(ijjj)*ww16 
     &   -bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     &   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     &   +bzdv(jjji)*ww07+bzdv(ijji)*ww08+bzdv(iiij)*ww09
     &   +bzdv(jiij)*ww10+bzdv(jjij)*ww11+bzdv(ijij)*ww12
     &   +bzdv(iijj)*ww13+bzdv(jijj)*ww14+bzdv(jjjj)*ww15
     &   +bzdv(ijjj)*ww16 
     &   -bfac1*y(3)*y(3) + b0/r2/r
	  
      if (noe.eqv..true.) then
         ex=0.
         ey=0.
         ez=0.
      else
         ex=exdv(iiii)*ww01+exdv(jiii)*ww02+exdv(jjii)*ww03
     &   +exdv(ijii)*ww04+exdv(iiji)*ww05+exdv(jiji)*ww06
     &   +exdv(jjji)*ww07+exdv(ijji)*ww08+exdv(iiij)*ww09
     &   +exdv(jiij)*ww10+exdv(jjij)*ww11+exdv(ijij)*ww12
     &   +exdv(iijj)*ww13+exdv(jijj)*ww14+exdv(jjjj)*ww15
     &   +exdv(ijjj)*ww16 

         ey=eydv(iiii)*ww01+eydv(jiii)*ww02+eydv(jjii)*ww03
     &   +eydv(ijii)*ww04+eydv(iiji)*ww05+eydv(jiji)*ww06
     &   +eydv(jjji)*ww07+eydv(ijji)*ww08+eydv(iiij)*ww09
     &   +eydv(jiij)*ww10+eydv(jjij)*ww11+eydv(ijij)*ww12
     &   +eydv(iijj)*ww13+eydv(jijj)*ww14+eydv(jjjj)*ww15
     &   +eydv(ijjj)*ww16 

         ez=ezdv(iiii)*ww01+ezdv(jiii)*ww02+ezdv(jjii)*ww03
     &   +ezdv(ijii)*ww04+ezdv(iiji)*ww05+ezdv(jiji)*ww06
     &   +ezdv(jjji)*ww07+ezdv(ijji)*ww08+ezdv(iiij)*ww09
     &   +ezdv(jiij)*ww10+ezdv(jjij)*ww11+ezdv(ijij)*ww12
     &   +ezdv(iijj)*ww13+ezdv(jijj)*ww14+ezdv(jjjj)*ww15
     &   +ezdv(ijjj)*ww16 
      endif

      dbxdx=dbxdxv(iiii)*ww01+dbxdxv(jiii)*ww02+dbxdxv(jjii)*ww03
     &   +dbxdxv(ijii)*ww04+dbxdxv(iiji)*ww05+dbxdxv(jiji)*ww06
     &   +dbxdxv(jjji)*ww07+dbxdxv(ijji)*ww08+dbxdxv(iiij)*ww09
     &   +dbxdxv(jiij)*ww10+dbxdxv(jjij)*ww11+dbxdxv(ijij)*ww12
     &   +dbxdxv(iijj)*ww13+dbxdxv(jijj)*ww14+dbxdxv(jjjj)*ww15
     &   +dbxdxv(ijjj)*ww16
     &   -bfac1*y(3)+bfac2*y(1)*y(1)*y(3)

      dbxdy=dbxdyv(iiii)*ww01+dbxdyv(jiii)*ww02+dbxdyv(jjii)*ww03
     &   +dbxdyv(ijii)*ww04+dbxdyv(iiji)*ww05+dbxdyv(jiji)*ww06
     &   +dbxdyv(jjji)*ww07+dbxdyv(ijji)*ww08+dbxdyv(iiij)*ww09
     &   +dbxdyv(jiij)*ww10+dbxdyv(jjij)*ww11+dbxdyv(ijij)*ww12
     &   +dbxdyv(iijj)*ww13+dbxdyv(jijj)*ww14+dbxdyv(jjjj)*ww15
     &   +dbxdyv(ijjj)*ww16 
     &   +bfac2*y(1)*y(2)*y(3)

      dbxdz=dbxdzv(iiii)*ww01+dbxdzv(jiii)*ww02+dbxdzv(jjii)*ww03
     &   +dbxdzv(ijii)*ww04+dbxdzv(iiji)*ww05+dbxdzv(jiji)*ww06
     &   +dbxdzv(jjji)*ww07+dbxdzv(ijji)*ww08+dbxdzv(iiij)*ww09
     &   +dbxdzv(jiij)*ww10+dbxdzv(jjij)*ww11+dbxdzv(ijij)*ww12
     &   +dbxdzv(iijj)*ww13+dbxdzv(jijj)*ww14+dbxdzv(jjjj)*ww15
     &   +dbxdzv(ijjj)*ww16 
     &   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbydx=dbydxv(iiii)*ww01+dbydxv(jiii)*ww02+dbydxv(jjii)*ww03
     &   +dbydxv(ijii)*ww04+dbydxv(iiji)*ww05+dbydxv(jiji)*ww06
     &   +dbydxv(jjji)*ww07+dbydxv(ijji)*ww08+dbydxv(iiij)*ww09
     &   +dbydxv(jiij)*ww10+dbydxv(jjij)*ww11+dbydxv(ijij)*ww12
     &   +dbydxv(iijj)*ww13+dbydxv(jijj)*ww14+dbydxv(jjjj)*ww15
     &   +dbydxv(ijjj)*ww16 
     &   +bfac2*y(2)*y(3)*y(1)

      dbydy=dbydyv(iiii)*ww01+dbydyv(jiii)*ww02+dbydyv(jjii)*ww03
     &   +dbydyv(ijii)*ww04+dbydyv(iiji)*ww05+dbydyv(jiji)*ww06
     &   +dbydyv(jjji)*ww07+dbydyv(ijji)*ww08+dbydyv(iiij)*ww09
     &   +dbydyv(jiij)*ww10+dbydyv(jjij)*ww11+dbydyv(ijij)*ww12
     &   +dbydyv(iijj)*ww13+dbydyv(jijj)*ww14+dbydyv(jjjj)*ww15
     &   +dbydyv(ijjj)*ww16 
     &   -bfac1*y(3)+bfac2*y(2)*y(2)*y(3)

      dbydz=dbydzv(iiii)*ww01+dbydzv(jiii)*ww02+dbydzv(jjii)*ww03
     &   +dbydzv(ijii)*ww04+dbydzv(iiji)*ww05+dbydzv(jiji)*ww06
     &   +dbydzv(jjji)*ww07+dbydzv(ijji)*ww08+dbydzv(iiij)*ww09
     &   +dbydzv(jiij)*ww10+dbydzv(jjij)*ww11+dbydzv(ijij)*ww12
     &   +dbydzv(iijj)*ww13+dbydzv(jijj)*ww14+dbydzv(jjjj)*ww15
     &   +dbydzv(ijjj)*ww16 
     &   -bfac1*y(2)+bfac2*y(3)*y(3)*y(2)

      dbzdx=dbzdxv(iiii)*ww01+dbzdxv(jiii)*ww02+dbzdxv(jjii)*ww03
     &   +dbzdxv(ijii)*ww04+dbzdxv(iiji)*ww05+dbzdxv(jiji)*ww06
     &   +dbzdxv(jjji)*ww07+dbzdxv(ijji)*ww08+dbzdxv(iiij)*ww09
     &   +dbzdxv(jiij)*ww10+dbzdxv(jjij)*ww11+dbzdxv(ijij)*ww12
     &   +dbzdxv(iijj)*ww13+dbzdxv(jijj)*ww14+dbzdxv(jjjj)*ww15
     &   +dbzdxv(ijjj)*ww16 
     &   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbzdy=dbzdyv(iiii)*ww01+dbzdyv(jiii)*ww02+dbzdyv(jjii)*ww03
     &   +dbzdyv(ijii)*ww04+dbzdyv(iiji)*ww05+dbzdyv(jiji)*ww06
     &   +dbzdyv(jjji)*ww07+dbzdyv(ijji)*ww08+dbzdyv(iiij)*ww09
     &   +dbzdyv(jiij)*ww10+dbzdyv(jjij)*ww11+dbzdyv(ijij)*ww12
     &   +dbzdyv(iijj)*ww13+dbzdyv(jijj)*ww14+dbzdyv(jjjj)*ww15
     &   +dbzdyv(ijjj)*ww16 
     &   -bfac1*y(2)+bfac2*y(2)*y(3)*y(3)

      dbzdz=dbzdzv(iiii)*ww01+dbzdzv(jiii)*ww02+dbzdzv(jjii)*ww03
     &   +dbzdzv(ijii)*ww04+dbzdzv(iiji)*ww05+dbzdzv(jiji)*ww06
     &   +dbzdzv(jjji)*ww07+dbzdzv(ijji)*ww08+dbzdzv(iiij)*ww09
     &   +dbzdzv(jiij)*ww10+dbzdzv(jjij)*ww11+dbzdzv(ijij)*ww12
     &   +dbzdzv(iijj)*ww13+dbzdzv(jijj)*ww14+dbzdzv(jjjj)*ww15
     &   +dbzdzv(ijjj)*ww16 
     &   -3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)

      dbxdt=dbxdtv(iiii)*ww01+dbxdtv(jiii)*ww02+dbxdtv(jjii)*ww03
     &   +dbxdtv(ijii)*ww04+dbxdtv(iiji)*ww05+dbxdtv(jiji)*ww06
     &   +dbxdtv(jjji)*ww07+dbxdtv(ijji)*ww08+dbxdtv(iiij)*ww09
     &   +dbxdtv(jiij)*ww10+dbxdtv(jjij)*ww11+dbxdtv(ijij)*ww12
     &   +dbxdtv(iijj)*ww13+dbxdtv(jijj)*ww14+dbxdtv(jjjj)*ww15
     &   +dbxdtv(ijjj)*ww16

      dbydt=dbydtv(iiii)*ww01+dbydtv(jiii)*ww02+dbydtv(jjii)*ww03
     &   +dbydtv(ijii)*ww04+dbydtv(iiji)*ww05+dbydtv(jiji)*ww06
     &   +dbydtv(jjji)*ww07+dbydtv(ijji)*ww08+dbydtv(iiij)*ww09
     &   +dbydtv(jiij)*ww10+dbydtv(jjij)*ww11+dbydtv(ijij)*ww12
     &   +dbydtv(iijj)*ww13+dbydtv(jijj)*ww14+dbydtv(jjjj)*ww15
     &   +dbydtv(ijjj)*ww16

      dbzdt=dbzdtv(iiii)*ww01+dbzdtv(jiii)*ww02+dbzdtv(jjii)*ww03
     &   +dbzdtv(ijii)*ww04+dbzdtv(iiji)*ww05+dbzdtv(jiji)*ww06
     &   +dbzdtv(jjji)*ww07+dbzdtv(ijji)*ww08+dbzdtv(iiij)*ww09
     &   +dbzdtv(jiij)*ww10+dbzdtv(jjij)*ww11+dbzdtv(ijij)*ww12
     &   +dbzdtv(iijj)*ww13+dbzdtv(jijj)*ww14+dbzdtv(jjjj)*ww15
     &   +dbzdtv(ijjj)*ww16

      b=sqrt(bx**2.+by**2.+bz**2.)
      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b
      dbdt=(bx*dbxdt + by*dbydt + bz*dbzdt)/b

*      print *,'*** b=',b

*      if ((t/tfactor.gt.362.20657).and.(t/tfactor.lt.370.0))then
*         print *,'t,i,j,k,x,y,z,b=',t/tfactor,ix,iy,iz,y(1),y(2),y(3),
*     &   b/ffactor/ntg
*      endif

      return
      end

************************************************************************

      subroutine get_fields2_old(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the 
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields, using the value
c of B0 specified in rbelt-const.inc.
c this version also calculates components of grad B and grad |B| 
c (needed for the guiding center routines) using a first order 
c approximation for purtabation fileds.

c in this old version, we get gradients of B in each cell using (b_i+1 - b_i)/dx
c which gives us a grad B that is consistent with fields, however discontinuous 
c gradients at grid cell boundaries result in pitch angle scattering.

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
      real r,r2,bfac1,bfac2

*      print *
*      print *,'*** in subroutine get_fields2 ***'
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
      elseif (r.gt.rmax) then
         status=2
         return
      endif
c     then check to see if it's time to halt the integration
*      if (t.ge.tstop) then
*         if (status.eq.0) status=-1
*      endif

c     uncomment to debug
c     make sure that (y,t) lies in the grid
*      if ((y(1).lt.xgr(1)).or.(y(1).ge.xgr(nx)).or.
*     &   (y(2).lt.ygr(1)).or.(y(2).ge.ygr(ny)).or.
*     &   (y(3).lt.zgr(1)).or.(y(3).ge.zgr(nz)).or.
*     &   (t.lt.tgr(1)).or.(t.ge.tgr(nt))) then
*         status=1
*         print *
*         print *,'in linterp4: outside grid'
*         print *,'t,y=',t/tfactor,y
*         if (y(1).lt.xgr(1)) print*,'y(1) < xgr(1)'
*         if (y(1).ge.xgr(nx)) print*,'y(1) >= xgr(nx)'
*         if (y(2).lt.ygr(1)) print*,'y(2) < ygr(1)'
*         if (y(2).ge.ygr(ny)) print*,'y(2) >= ygr(ny)'
*         if (y(3).lt.zgr(1)) print*,'y(3) < zgr(1)'
*         if (y(3).ge.zgr(nz)) print*,'y(3) >= zgr(nz)'
*         if (t.lt.tgr(1)) print*,'t < tgr(1)',t,tgr(1)
*         if (t.ge.tgr(nt)) print*,'t >= tgr(nt)'
*         return
*      endif

      itu=nt+1
      it=0
10    continue
      if ((itu-it).gt.1) then
        itm=(itu+it)/2
        if (t.ge.tgr(itm)) then
          it=itm
        else
          itu=itm
        endif
        goto 10
      endif
      dt=tgr(it+1)-tgr(it)

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1


*      if ((ix.lt.1).or.(ix.gt.nx).or.
*     &   (iy.lt.1).or.(iy.gt.ny).or.
*     &   (iz.lt.1).or.(iz.gt.nz).or.
*     &   (it.lt.1).or.(it.gt.nt)) then
*         print *
*         print *,'ix=',ix
*         print *,'iy=',iy
*         print *,'iz=',iz
*         print *,'it=',it
*         print *,'t,r=',t/tfactor,r
*         print *,'status=',status
*         print *,'tgrmax=',tgrmax/tfactor
**         status=2
*         stop
*         return
*      endif

      jjjj= ix + iy*nx + iz*nxy + it*nxyz + 1
      ijjj= jjjj - 1
      jijj= jjjj - nx
      iijj= jijj - 1
      jjij= jjjj - nxy
      ijij= jjij - 1
      jiij= jijj - nxy
      iiij= jiij - 1
      jjji= jjjj - nxyz
      ijji= ijjj - nxyz
      jiji= jijj - nxyz
      iiji= iijj - nxyz
      jjii= jjij - nxyz
      ijii= ijij - nxyz
      jiii= jiij - nxyz
      iiii= iiij - nxyz

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3
      w4=abs(t-tgr(it))/dt
      w4m=1.-w4

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

      ww01=w1m2m3m*w4m
      ww02=w12m3m*w4m
      ww03=w123m*w4m
      ww04=w1m23m*w4m
      ww05=w1m2m3*w4m
      ww06=w12m3*w4m
      ww07=w123*w4m
      ww08=w1m23*w4m
      ww09=w1m2m3m*w4
      ww10=w12m3m*w4
      ww11=w123m*w4
      ww12=w1m23m*w4
      ww13=w1m2m3*w4
      ww14=w12m3*w4
      ww15=w123*w4
      ww16=w1m23*w4

      r2=r*r
      bfac1=3.*b0/r2/r2/r
      bfac2=5.*bfac1/r2

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     &   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     &   +bxdv(jjji)*ww07+bxdv(ijji)*ww08+bxdv(iiij)*ww09
     &   +bxdv(jiij)*ww10+bxdv(jjij)*ww11+bxdv(ijij)*ww12
     &   +bxdv(iijj)*ww13+bxdv(jijj)*ww14+bxdv(jjjj)*ww15
     &   +bxdv(ijjj)*ww16 
     &   -bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     &   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     &   +bydv(jjji)*ww07+bydv(ijji)*ww08+bydv(iiij)*ww09
     &   +bydv(jiij)*ww10+bydv(jjij)*ww11+bydv(ijij)*ww12
     &   +bydv(iijj)*ww13+bydv(jijj)*ww14+bydv(jjjj)*ww15
     &   +bydv(ijjj)*ww16 
     &   -bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     &   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     &   +bzdv(jjji)*ww07+bzdv(ijji)*ww08+bzdv(iiij)*ww09
     &   +bzdv(jiij)*ww10+bzdv(jjij)*ww11+bzdv(ijij)*ww12
     &   +bzdv(iijj)*ww13+bzdv(jijj)*ww14+bzdv(jjjj)*ww15
     &   +bzdv(ijjj)*ww16 
     &   -bfac1*y(3)*y(3) + b0/r2/r
	  
      if (noe.eqv..true.) then
         ex=0.
         ey=0.
         ez=0.
      else
         ex=exdv(iiii)*ww01+exdv(jiii)*ww02+exdv(jjii)*ww03
     &   +exdv(ijii)*ww04+exdv(iiji)*ww05+exdv(jiji)*ww06
     &   +exdv(jjji)*ww07+exdv(ijji)*ww08+exdv(iiij)*ww09
     &   +exdv(jiij)*ww10+exdv(jjij)*ww11+exdv(ijij)*ww12
     &   +exdv(iijj)*ww13+exdv(jijj)*ww14+exdv(jjjj)*ww15
     &   +exdv(ijjj)*ww16 

         ey=eydv(iiii)*ww01+eydv(jiii)*ww02+eydv(jjii)*ww03
     &   +eydv(ijii)*ww04+eydv(iiji)*ww05+eydv(jiji)*ww06
     &   +eydv(jjji)*ww07+eydv(ijji)*ww08+eydv(iiij)*ww09
     &   +eydv(jiij)*ww10+eydv(jjij)*ww11+eydv(ijij)*ww12
     &   +eydv(iijj)*ww13+eydv(jijj)*ww14+eydv(jjjj)*ww15
     &   +eydv(ijjj)*ww16 

         ez=ezdv(iiii)*ww01+ezdv(jiii)*ww02+ezdv(jjii)*ww03
     &   +ezdv(ijii)*ww04+ezdv(iiji)*ww05+ezdv(jiji)*ww06
     &   +ezdv(jjji)*ww07+ezdv(ijji)*ww08+ezdv(iiij)*ww09
     &   +ezdv(jiij)*ww10+ezdv(jjij)*ww11+ezdv(ijij)*ww12
     &   +ezdv(iijj)*ww13+ezdv(jijj)*ww14+ezdv(jjjj)*ww15
     &   +ezdv(ijjj)*ww16 
      endif

      b=sqrt(bx**2.+by**2.+bz**2.)

      w2m3m4m=w2m*w3m*w4m
      w23m4m=w2*w3m*w4m
      w2m34m=w2m*w3*w4m
      w234m=w2*w3*w4m
      w2m3m4=w2m*w3m*w4
      w23m4=w2*w3m*w4
      w2m34=w2m*w3*w4
      w234=w2*w3*w4

      w1m3m4m=w1m*w3m*w4m
      w13m4m=w1*w3m*w4m
      w1m34m=w1m*w3*w4m
      w134m=w1*w3*w4m
      w1m3m4=w1m*w3m*w4
      w13m4=w1*w3m*w4
      w1m34=w1m*w3*w4
      w134=w1*w3*w4

      w1m2m4m=w1m2m*w4m
      w12m4m=w12m*w4m
      w1m24m=w1m2*w4m
      w124m=w12*w4m
      w1m2m4=w1m2m*w4
      w12m4=w12m*w4
      w1m24=w1m2*w4
      w124=w12*w4

      dbxdx=((bxdv(jiii)-bxdv(iiii))*w2m3m4m
     &   +(bxdv(jjii)-bxdv(ijii))*w23m4m
     &   +(bxdv(jiji)-bxdv(iiji))*w2m34m
     &   +(bxdv(jjji)-bxdv(ijji))*w234m
     &   +(bxdv(jiij)-bxdv(iiij))*w2m3m4
     &   +(bxdv(jjij)-bxdv(ijij))*w23m4
     &   +(bxdv(jijj)-bxdv(iijj))*w2m34
     &   +(bxdv(jjjj)-bxdv(ijjj))*w234)/dx
     &   -bfac1*y(3)+bfac2*y(1)*y(1)*y(3)

      dbxdy=((bxdv(ijii)-bxdv(iiii))*w1m3m4m
     &   +(bxdv(jjii)-bxdv(jiii))*w13m4m
     &   +(bxdv(ijji)-bxdv(iiji))*w1m34m
     &   +(bxdv(jjji)-bxdv(jiji))*w134m
     &   +(bxdv(ijij)-bxdv(iiij))*w1m3m4
     &   +(bxdv(jjij)-bxdv(jiij))*w13m4
     &   +(bxdv(ijjj)-bxdv(iijj))*w1m34
     &   +(bxdv(jjjj)-bxdv(jijj))*w134)/dy
     &   +bfac2*y(1)*y(2)*y(3)

      dbxdz=((bxdv(iiji)-bxdv(iiii))*w1m2m4m
     &   +(bxdv(jiji)-bxdv(jiii))*w12m4m
     &   +(bxdv(ijji)-bxdv(ijii))*w1m24m
     &   +(bxdv(jjji)-bxdv(jjii))*w124m
     &   +(bxdv(iijj)-bxdv(iiij))*w1m2m4
     &   +(bxdv(jijj)-bxdv(jiij))*w12m4
     &   +(bxdv(ijjj)-bxdv(ijij))*w1m24
     &   +(bxdv(jjjj)-bxdv(jjij))*w124)/dz
     &   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbydx=((bydv(jiii)-bydv(iiii))*w2m3m4m
     &   +(bydv(jjii)-bydv(ijii))*w23m4m
     &   +(bydv(jiji)-bydv(iiji))*w2m34m
     &   +(bydv(jjji)-bydv(ijji))*w234m
     &   +(bydv(jiij)-bydv(iiij))*w2m3m4
     &   +(bydv(jjij)-bydv(ijij))*w23m4
     &   +(bydv(jijj)-bydv(iijj))*w2m34
     &   +(bydv(jjjj)-bydv(ijjj))*w234)/dx
     &   +bfac2*y(2)*y(3)*y(1)

      dbydy=((bydv(ijii)-bydv(iiii))*w1m3m4m
     &   +(bydv(jjii)-bydv(jiii))*w13m4m
     &   +(bydv(ijji)-bydv(iiji))*w1m34m
     &   +(bydv(jjji)-bydv(jiji))*w134m
     &   +(bydv(ijij)-bydv(iiij))*w1m3m4
     &   +(bydv(jjij)-bydv(jiij))*w13m4
     &   +(bydv(ijjj)-bydv(iijj))*w1m34
     &   +(bydv(jjjj)-bydv(jijj))*w134)/dy
     &   -bfac1*y(3)+bfac2*y(2)*y(2)*y(3)

      dbydz=((bydv(iiji)-bydv(iiii))*w1m2m4m
     &   +(bydv(jiji)-bydv(jiii))*w12m4m
     &   +(bydv(ijji)-bydv(ijii))*w1m24m
     &   +(bydv(jjji)-bydv(jjii))*w124m
     &   +(bydv(iijj)-bydv(iiij))*w1m2m4
     &   +(bydv(jijj)-bydv(jiij))*w12m4
     &   +(bydv(ijjj)-bydv(ijij))*w1m24
     &   +(bydv(jjjj)-bydv(jjij))*w124)/dz
     &   -bfac1*y(2)+bfac2*y(3)*y(3)*y(2)

      dbzdx=((bzdv(jiii)-bzdv(iiii))*w2m3m4m
     &   +(bzdv(jjii)-bzdv(ijii))*w23m4m
     &   +(bzdv(jiji)-bzdv(iiji))*w2m34m
     &   +(bzdv(jjji)-bzdv(ijji))*w234m
     &   +(bzdv(jiij)-bzdv(iiij))*w2m3m4
     &   +(bzdv(jjij)-bzdv(ijij))*w23m4
     &   +(bzdv(jijj)-bzdv(iijj))*w2m34
     &   +(bzdv(jjjj)-bzdv(ijjj))*w234)/dx
     &   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbzdy=((bzdv(ijii)-bzdv(iiii))*w1m3m4m
     &   +(bzdv(jjii)-bzdv(jiii))*w13m4m
     &   +(bzdv(ijji)-bzdv(iiji))*w1m34m
     &   +(bzdv(jjji)-bzdv(jiji))*w134m
     &   +(bzdv(ijij)-bzdv(iiij))*w1m3m4
     &   +(bzdv(jjij)-bzdv(jiij))*w13m4
     &   +(bzdv(ijjj)-bzdv(iijj))*w1m34
     &   +(bzdv(jjjj)-bzdv(jijj))*w134)/dy
     &   -bfac1*y(2)+bfac2*y(2)*y(3)*y(3)

      dbzdz=((bzdv(iiji)-bzdv(iiii))*w1m2m4m
     &   +(bzdv(jiji)-bzdv(jiii))*w12m4m
     &   +(bzdv(ijji)-bzdv(ijii))*w1m24m
     &   +(bzdv(jjji)-bzdv(jjii))*w124m
     &   +(bzdv(iijj)-bzdv(iiij))*w1m2m4
     &   +(bzdv(jijj)-bzdv(jiij))*w12m4
     &   +(bzdv(ijjj)-bzdv(ijij))*w1m24
     &   +(bzdv(jjjj)-bzdv(jjij))*w124)/dz
     &   -3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)

      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b

*      print *,'b=',b

      return
      end

************************************************************************

      subroutine get_fields2_test(y,t)
*      subroutine get_fields2(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the 
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields, using the value
c of B0 specified in rbelt-const.inc.
c this version also calculates components of grad B and grad |B| 
c (needed for the guiding center routines) using a first order 
c approximation for purtabation fileds.

c in this old version, we get gradients of B in each cell using (b_i+1 - b_i)/dx
c which gives us a grad B that is consistent with fields, however discontinuous 
c gradients at grid cell boundaries result in pitch angle scattering.

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
      real r,r2,bfac1,bfac2,tstop

      print *
      print *,'*** in subroutine get_fields2 ***'
      print *,'status,t,y=',status,t/tfactor,y
      print *,'tgrmax,tstop=',tgrmax/tfactor,tstop/tfactor

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
c     then check to see if it's time to halt the integration
*      if (t.ge.tstop) then
*         if (status.eq.0) status=-1
*      endif

c     uncomment to debug
c     make sure that (y,t) lies in the grid
*      if ((y(1).lt.xgr(1)).or.(y(1).ge.xgr(nx)).or.
*     &   (y(2).lt.ygr(1)).or.(y(2).ge.ygr(ny)).or.
*     &   (y(3).lt.zgr(1)).or.(y(3).ge.zgr(nz)).or.
*     &   (t.lt.tgr(1)).or.(t.ge.tgr(nt))) then
*         status=1
*         print *
*         print *,'in linterp4: outside grid'
*         print *,'t,y=',t/tfactor,y
*         if (y(1).lt.xgr(1)) print*,'y(1) < xgr(1)'
*         if (y(1).ge.xgr(nx)) print*,'y(1) >= xgr(nx)'
*         if (y(2).lt.ygr(1)) print*,'y(2) < ygr(1)'
*         if (y(2).ge.ygr(ny)) print*,'y(2) >= ygr(ny)'
*         if (y(3).lt.zgr(1)) print*,'y(3) < zgr(1)'
*         if (y(3).ge.zgr(nz)) print*,'y(3) >= zgr(nz)'
*         if (t.lt.tgr(1)) print*,'t < tgr(1)',t,tgr(1)
*         if (t.ge.tgr(nt)) print*,'t >= tgr(nt)'
*         return
*      endif

      itu=nt+1
      it=0
10    continue
      if ((itu-it).gt.1) then
        itm=(itu+it)/2
        if (t.ge.tgr(itm)) then
          it=itm
        else
          itu=itm
        endif
        goto 10
      endif
      dt=tgr(it+1)-tgr(it)

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1


*      if ((ix.lt.1).or.(ix.gt.nx).or.
*     &   (iy.lt.1).or.(iy.gt.ny).or.
*     &   (iz.lt.1).or.(iz.gt.nz).or.
*     &   (it.lt.1).or.(it.gt.nt)) then
*         print *
*         print *,'ix=',ix
*         print *,'iy=',iy
*         print *,'iz=',iz
*         print *,'it=',it
*         print *,'t,r=',t/tfactor,r
*         print *,'status=',status
*         print *,'tgrmax=',tgrmax/tfactor
**         status=2
*         stop
*         return
*      endif

      jjjj= ix + iy*nx + iz*nxy + it*nxyz + 1
      ijjj= jjjj - 1
      jijj= jjjj - nx
      iijj= jijj - 1
      jjij= jjjj - nxy
      ijij= jjij - 1
      jiij= jijj - nxy
      iiij= jiij - 1
      jjji= jjjj - nxyz
      ijji= ijjj - nxyz
      jiji= jijj - nxyz
      iiji= iijj - nxyz
      jjii= jjij - nxyz
      ijii= ijij - nxyz
      jiii= jiij - nxyz
      iiii= iiij - nxyz

      w1=abs(y(1)-xgr(ix))/dx
      w1m=1.-w1
      w2=abs(y(2)-ygr(iy))/dy
      w2m=1.-w2
      w3=abs(y(3)-zgr(iz))/dz
      w3m=1.-w3
      w4=abs(t-tgr(it))/dt
      w4m=1.-w4

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

      ww01=w1m2m3m*w4m
      ww02=w12m3m*w4m
      ww03=w123m*w4m
      ww04=w1m23m*w4m
      ww05=w1m2m3*w4m
      ww06=w12m3*w4m
      ww07=w123*w4m
      ww08=w1m23*w4m
      ww09=w1m2m3m*w4
      ww10=w12m3m*w4
      ww11=w123m*w4
      ww12=w1m23m*w4
      ww13=w1m2m3*w4
      ww14=w12m3*w4
      ww15=w123*w4
      ww16=w1m23*w4

      r2=r*r
      bfac1=3.*b0/r2/r2/r
      bfac2=5.*bfac1/r2

      bx=bxdv(iiii)*ww01+bxdv(jiii)*ww02+bxdv(jjii)*ww03
     &   +bxdv(ijii)*ww04+bxdv(iiji)*ww05+bxdv(jiji)*ww06
     &   +bxdv(jjji)*ww07+bxdv(ijji)*ww08+bxdv(iiij)*ww09
     &   +bxdv(jiij)*ww10+bxdv(jjij)*ww11+bxdv(ijij)*ww12
     &   +bxdv(iijj)*ww13+bxdv(jijj)*ww14+bxdv(jjjj)*ww15
     &   +bxdv(ijjj)*ww16 
     &   -bfac1*y(1)*y(3)

      by=bydv(iiii)*ww01+bydv(jiii)*ww02+bydv(jjii)*ww03
     &   +bydv(ijii)*ww04+bydv(iiji)*ww05+bydv(jiji)*ww06
     &   +bydv(jjji)*ww07+bydv(ijji)*ww08+bydv(iiij)*ww09
     &   +bydv(jiij)*ww10+bydv(jjij)*ww11+bydv(ijij)*ww12
     &   +bydv(iijj)*ww13+bydv(jijj)*ww14+bydv(jjjj)*ww15
     &   +bydv(ijjj)*ww16 
     &   -bfac1*y(2)*y(3)

      bz=bzdv(iiii)*ww01+bzdv(jiii)*ww02+bzdv(jjii)*ww03
     &   +bzdv(ijii)*ww04+bzdv(iiji)*ww05+bzdv(jiji)*ww06
     &   +bzdv(jjji)*ww07+bzdv(ijji)*ww08+bzdv(iiij)*ww09
     &   +bzdv(jiij)*ww10+bzdv(jjij)*ww11+bzdv(ijij)*ww12
     &   +bzdv(iijj)*ww13+bzdv(jijj)*ww14+bzdv(jjjj)*ww15
     &   +bzdv(ijjj)*ww16 
     &   -bfac1*y(3)*y(3) + b0/r2/r
	  
      if (noe.eqv..true.) then
         ex=0.
         ey=0.
         ez=0.
      else
         ex=exdv(iiii)*ww01+exdv(jiii)*ww02+exdv(jjii)*ww03
     &   +exdv(ijii)*ww04+exdv(iiji)*ww05+exdv(jiji)*ww06
     &   +exdv(jjji)*ww07+exdv(ijji)*ww08+exdv(iiij)*ww09
     &   +exdv(jiij)*ww10+exdv(jjij)*ww11+exdv(ijij)*ww12
     &   +exdv(iijj)*ww13+exdv(jijj)*ww14+exdv(jjjj)*ww15
     &   +exdv(ijjj)*ww16 

         ey=eydv(iiii)*ww01+eydv(jiii)*ww02+eydv(jjii)*ww03
     &   +eydv(ijii)*ww04+eydv(iiji)*ww05+eydv(jiji)*ww06
     &   +eydv(jjji)*ww07+eydv(ijji)*ww08+eydv(iiij)*ww09
     &   +eydv(jiij)*ww10+eydv(jjij)*ww11+eydv(ijij)*ww12
     &   +eydv(iijj)*ww13+eydv(jijj)*ww14+eydv(jjjj)*ww15
     &   +eydv(ijjj)*ww16 

         ez=ezdv(iiii)*ww01+ezdv(jiii)*ww02+ezdv(jjii)*ww03
     &   +ezdv(ijii)*ww04+ezdv(iiji)*ww05+ezdv(jiji)*ww06
     &   +ezdv(jjji)*ww07+ezdv(ijji)*ww08+ezdv(iiij)*ww09
     &   +ezdv(jiij)*ww10+ezdv(jjij)*ww11+ezdv(ijij)*ww12
     &   +ezdv(iijj)*ww13+ezdv(jijj)*ww14+ezdv(jjjj)*ww15
     &   +ezdv(ijjj)*ww16 
      endif

      b=sqrt(bx**2.+by**2.+bz**2.)

      w2m3m4m=w2m*w3m*w4m
      w23m4m=w2*w3m*w4m
      w2m34m=w2m*w3*w4m
      w234m=w2*w3*w4m
      w2m3m4=w2m*w3m*w4
      w23m4=w2*w3m*w4
      w2m34=w2m*w3*w4
      w234=w2*w3*w4

      w1m3m4m=w1m*w3m*w4m
      w13m4m=w1*w3m*w4m
      w1m34m=w1m*w3*w4m
      w134m=w1*w3*w4m
      w1m3m4=w1m*w3m*w4
      w13m4=w1*w3m*w4
      w1m34=w1m*w3*w4
      w134=w1*w3*w4

      w1m2m4m=w1m2m*w4m
      w12m4m=w12m*w4m
      w1m24m=w1m2*w4m
      w124m=w12*w4m
      w1m2m4=w1m2m*w4
      w12m4=w12m*w4
      w1m24=w1m2*w4
      w124=w12*w4

      dbxdx=dbxdxv(iiii)*w2m3m4m 
     &   +dbxdxv(ijii)*w23m4m
     &   +dbxdxv(iiji)*w2m34m
     &   +dbxdxv(ijji)*w234m
     &   +dbxdxv(iiij)*w2m3m4
     &   +dbxdxv(ijij)*w23m4
     &   +dbxdxv(iijj)*w2m34
     &   +dbxdxv(ijjj)*w234
     &   -bfac1*y(3)+bfac2*y(1)*y(1)*y(3)

      dbxdy=dbxdyv(iiii)*w1m3m4m
     &   +dbxdyv(jiii)*w13m4m
     &   +dbxdyv(iiji)*w1m34m
     &   +dbxdyv(jiji)*w134m
     &   +dbxdyv(iiij)*w1m3m4
     &   +dbxdyv(jiij)*w13m4
     &   +dbxdyv(iijj)*w1m34
     &   +dbxdyv(jijj)*w134
     &   +bfac2*y(1)*y(2)*y(3)

      dbxdz=dbxdzv(iiii)*w1m2m4m
     &   +dbxdzv(jiii)*w12m4m
     &   +dbxdzv(ijii)*w1m24m
     &   +dbxdzv(jjii)*w124m
     &   +dbxdzv(iiij)*w1m2m4
     &   +dbxdzv(jiij)*w12m4
     &   +dbxdzv(ijij)*w1m24
     &   +dbxdzv(jjij)*w124
     &   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbydx=dbydxv(iiii)*w2m3m4m
     &   +dbydxv(ijii)*w23m4m
     &   +dbydxv(iiji)*w2m34m
     &   +dbydxv(ijji)*w234m
     &   +dbydxv(iiij)*w2m3m4
     &   +dbydxv(ijij)*w23m4
     &   +dbydxv(iijj)*w2m34
     &   +dbydxv(ijjj)*w234
     &   +bfac2*y(2)*y(3)*y(1)

      dbydy=dbydyv(iiii)*w1m3m4m
     &   +dbydyv(jiii)*w13m4m
     &   +dbydyv(iiji)*w1m34m
     &   +dbydyv(jiji)*w134m
     &   +dbydyv(iiij)*w1m3m4
     &   +dbydyv(jiij)*w13m4
     &   +dbydyv(iijj)*w1m34
     &   +dbydyv(jijj)*w134
     &   -bfac1*y(3)+bfac2*y(2)*y(2)*y(3)

      dbydz=dbydzv(iiii)*w1m2m4m
     &   +dbydzv(jiii)*w12m4m
     &   +dbydzv(ijii)*w1m24m
     &   +dbydzv(jjii)*w124m
     &   +dbydzv(iiij)*w1m2m4
     &   +dbydzv(jiij)*w12m4
     &   +dbydzv(ijij)*w1m24
     &   +dbydzv(jjij)*w124
     &   -bfac1*y(2)+bfac2*y(3)*y(3)*y(2)

      dbzdx=dbzdxv(iiii)*w2m3m4m
     &   +dbzdxv(ijii)*w23m4m
     &   +dbzdxv(iiji)*w2m34m
     &   +dbzdxv(ijji)*w234m
     &   +dbzdxv(iiij)*w2m3m4
     &   +dbzdxv(ijij)*w23m4
     &   +dbzdxv(iijj)*w2m34
     &   +dbzdxv(ijjj)*w234
     &   -bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
	  
      dbzdy=dbzdyv(iiii)*w1m3m4m
     &   +dbzdyv(jiii)*w13m4m
     &   +dbzdyv(iiji)*w1m34m
     &   +dbzdyv(jiji)*w134m
     &   +dbzdyv(iiij)*w1m3m4
     &   +dbzdyv(jiij)*w13m4
     &   +dbzdyv(iijj)*w1m34
     &   +dbzdyv(jijj)*w134
     &   -bfac1*y(2)+bfac2*y(2)*y(3)*y(3)

      dbzdz=dbzdzv(iiii)*w1m2m4m
     &   +dbzdzv(jiii)*w12m4m
     &   +dbzdzv(ijii)*w1m24m
     &   +dbzdzv(jjii)*w124m
     &   +dbzdzv(iiij)*w1m2m4
     &   +dbzdzv(jiij)*w12m4
     &   +dbzdzv(ijij)*w1m24
     &   +dbzdzv(jjij)*w124
     &   -3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)

      dbxdt=dbxdtv(iiii)*w1m2m3m
     &   +dbxdtv(jiii)*w12m3m
     &   +dbxdtv(ijii)*w1m23m
     &   +dbxdtv(iiji)*w1m2m3
     &   +dbxdtv(jjii)*w123m
     &   +dbxdtv(ijji)*w1m23
     &   +dbxdtv(jiji)*w12m3
     &   +dbxdtv(jjji)*w123

      dbydt=dbydtv(iiii)*w1m2m3m
     &   +dbydtv(jiii)*w12m3m
     &   +dbydtv(ijii)*w1m23m
     &   +dbydtv(iiji)*w1m2m3
     &   +dbydtv(jjii)*w123m
     &   +dbydtv(ijji)*w1m23
     &   +dbydtv(jiji)*w12m3
     &   +dbydtv(jjji)*w123

      dbzdt=dbzdtv(iiii)*w1m2m3m
     &   +dbzdtv(jiii)*w12m3m
     &   +dbzdtv(ijii)*w1m23m
     &   +dbzdtv(iiji)*w1m2m3
     &   +dbzdtv(jjii)*w123m
     &   +dbzdtv(ijji)*w1m23
     &   +dbzdtv(jiji)*w12m3
     &   +dbzdtv(jjji)*w123

      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b
      dbdt=(bx*dbxdt + by*dbydt + bz*dbzdt)/b

*      print *,'b=',b

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


