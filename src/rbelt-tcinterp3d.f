c     rbelt-tcinterp3d.f - tricubic interpolation for time dependent fields

************************************************************************

      subroutine get_fields(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the 
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields here, using the value
c of B0 specified in rbelt-const.inc.

      implicit none
      include 'rbelt-intmat.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      include 'rbelt-bounds.inc'

      real*8 y(6),t,px,py,pz,alpha(64),fp
      real*8 dt,fpt1,fpt2,r,r2,bfac1
      integer itu,itm,it,ix,iy,iz
      integer i,j,k,l,i1,i2,f
     
      if (status.gt.0) return
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      r=dsqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
      if (r.lt.rmin) then
         status=1
         return
      elseif (r.gt.rmax) then
         status=2
         return
      endif
 
      itu=nstep+1
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
      dt=tgr(itu)-tgr(it)
      
      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1
      
      px = (y(1)-xgr(ix))/dx
      py = (y(2)-ygr(iy))/dy
      pz = (y(3)-zgr(iz))/dz
      
      r=dsqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
      r2=r*r
      bfac1=3.*b0/r2/r2/r

      do 100 i1 = 1, 6
        select case (i1)
          case (1)
            f = 1
          case (2)
            f = 2
          case (3)
            f = 3
          case (4)
            f = 4
          case (5)
            f = 5
          case (6)
            f = 6
c            do 34 i = 1, 4
c              do 33 j = 1, 4
c                do 32 k = 1, 4
c                  do 31 l = 1, 2
c                    f(i,j,k,l) = ezd(i,j,k,l)
c31                continue
c32              continue
c33            continue
c34          continue            
        end select
        
          do 90 i2 = 1, 2
                  
c     create an if statement to make sure the particle is in bounds
c     if the particle is at the edge of the grid the derivatives 
c     will fail          
           if (i2.eq.1) then
             l = it
*            call allderiv(f,x,ix,iy,iz,it,1)
*            call allderiv(f,x,ix+1,iy,iz,it,2)
*            call allderiv(f,x,ix,iy+1,iz,it,3)
*            call allderiv(f,x,ix+1,iy+1,iz,it,4)
*            call allderiv(f,x,ix,iy,iz+1,it,5)
*            call allderiv(f,x,ix+1,iy,iz+1,it,6)
*            call allderiv(f,x,ix,iy+1,iz+1,it,7)
*            call allderiv(f,x,ix+1,iy+1,iz+1,it,8)
          elseif (i2.eq.2) then
            l = itu
*            call allderiv(f,x,ix,iy,iz,itu,1)
*            call allderiv(f,x,ix+1,iy,iz,itu,2)
*            call allderiv(f,x,ix,iy+1,iz,itu,3)
*            call allderiv(f,x,ix+1,iy+1,iz,itu,4)
*            call allderiv(f,x,ix,iy,iz+1,itu,5)
*            call allderiv(f,x,ix+1,iy,iz+1,itu,6)
*            call allderiv(f,x,ix,iy+1,iz+1,itu,7)
*            call allderiv(f,x,ix+1,iy+1,iz+1,itu,8)
          endif
            
c            print*,'x = ', x(26),x(27),x(28),x(29),x(30),x(31),x(25)
      
            do 50 i = 1, 64
              alpha(i) = 0.0
              do 40 j = 1, 64
                alpha(i) = alpha(i)+Amat(j,i)*alphamat(ix,iy,iz,l,f,j)
40            continue
50          continue

c            print*,'alpha = ', alpha(1),alpha(2),alpha(3),alpha(4)

            fp = 0.0
            do 80 i = 0, 3
              do 70 j = 0, 3
                do 60 k = 0, 3
                  fp = fp + alpha(1+i+4*j+16*k)*px**i*py**j*pz**k
60              continue
70            continue
80          continue

            if (i2.eq.1) then
              fpt1 = fp
            elseif (i2.eq.2) then
              fpt2 = fp
c       linear interpolation in the time variable
              fp = fpt1 + (fpt2-fpt1)*(t-tgr(it))/dt
            endif
    
90        continue       
          select case (i1)
            case (1)
              bx = fp - bfac1*y(1)*y(3)
            case (2)
              by = fp - bfac1*y(2)*y(3)
            case (3)
              bz = fp - bfac1*y(3)*y(3) + b0/r2/r
            case (4)
              ex = fp
            case (5)
              ey = fp
            case (6)
              ez = fp
          end select
      
100   continue

      b = dsqrt(bx**2+by**2+bz**2)    
      return
      end

************************************************************************

      subroutine get_fields2(y,t)

c quad-linear interpolation of gridded E & B fields to particle position
c from data on a uniform cartisian space and time grid. note that the
c analytic expression for the dipole component of the geomagnetic field
c is added to the interpolated purtabation fields here, using the value
c of B0 specified in rbelt-const.inc.

      implicit none
      include 'rbelt-intmat.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-const.inc'
      include 'rbelt-status.inc'
      include 'rbelt-bounds.inc'

      real*8 y(6),t,px,py,pz,alpha(64),fp
      real*8 dt,fpt1,fpt2,r,r2,bfac1,bfac2,fpx,fpy,fpz,fpt
      integer itu,itm,it,ix,iy,iz
      integer i,j,k,l,i1,i2,f

      if (status.gt.0) return
      if (t.ge.tgrmax) then
         status=3
         return
      endif
      r=dsqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
      if (r.lt.rmin) then
         status=1
         return
      elseif (r.gt.rmax) then
         status=2
         return
      endif

      itu=nstep+1
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
      dt=tgr(itu)-tgr(it)

      ix=((y(1)-xgr(1))/dx)+1
      iy=((y(2)-ygr(1))/dy)+1
      iz=((y(3)-zgr(1))/dz)+1

      px = (y(1)-xgr(ix))/dx
      py = (y(2)-ygr(iy))/dy
      pz = (y(3)-zgr(iz))/dz

      r=dsqrt(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))
      r2=r*r
      bfac1=3.*b0/r2/r2/r
      bfac2=5.*bfac1/r2

      do 100 i1 = 1, 6
        select case (i1)
          case (1)
            f = 1
          case (2)
            f = 2
          case (3)
            f = 3
          case (4)
            f = 4
          case (5)
            f = 5
          case (6)
            f = 6

        end select

          do 90 i2 = 1, 2
      
c     create an if statement to make sure the particle is in bounds
c     if the particle is at the edge of the grid the derivatives
c     will fail

           if (i2.eq.1) then
            l = it
*            call allderiv(f,x,ix,iy,iz,it,1)
*            call allderiv(f,x,ix+1,iy,iz,it,2)
*            call allderiv(f,x,ix,iy+1,iz,it,3)
*            call allderiv(f,x,ix+1,iy+1,iz,it,4)
*            call allderiv(f,x,ix,iy,iz+1,it,5)
*            call allderiv(f,x,ix+1,iy,iz+1,it,6)
*            call allderiv(f,x,ix,iy+1,iz+1,it,7)
*            call allderiv(f,x,ix+1,iy+1,iz+1,it,8)
          elseif (i2.eq.2) then
            l = itu
*            call allderiv(f,x,ix,iy,iz,itu,1)
*            call allderiv(f,x,ix+1,iy,iz,itu,2)
*            call allderiv(f,x,ix,iy+1,iz,itu,3)
*            call allderiv(f,x,ix+1,iy+1,iz,itu,4)
*            call allderiv(f,x,ix,iy,iz+1,itu,5)
*            call allderiv(f,x,ix+1,iy,iz+1,itu,6)
*            call allderiv(f,x,ix,iy+1,iz+1,itu,7)
*            call allderiv(f,x,ix+1,iy+1,iz+1,itu,8)
          endif

            do 50 i = 1, 64
              alpha(i) = 0.0
              do 40 j = 1, 64
                alpha(i) = alpha(i)+Amat(j,i)*alphamat(ix,iy,iz,l,f,j)
40            continue
50          continue

            fp = 0.0
            do 80 i = 0, 3
              do 70 j = 0, 3
                do 60 k = 0, 3
                  fp = fp + alpha(1+i+4*j+16*k)*px**i*py**j*pz**k
60              continue
70            continue
80          continue

          if (i2.eq.1) then
            fpx = 0.0
            do i = 0,3
              do j = 0,3
                do k = 0,3
                  fpx = fpx+alpha(1+i+4*j+16*k)*i*px**(i-1)*py**j*pz**k
                enddo
              enddo
            enddo

            fpy = 0.0
            do i = 0,3
              do j = 0,3
                do k = 0,3
                  fpy = fpy+alpha(1+i+4*j+16*k)*j*px**i*py**(j-1)*pz**k
                enddo
              enddo
            enddo

            fpz = 0.0
            do i = 0,3
              do j = 0,3
                do k = 0,3
                  fpz = fpz+alpha(1+i+4*j+16*k)*k*px**i*py**j*pz**(k-1)
                enddo
              enddo
            enddo
          endif

            if (i2.eq.1) then
              fpt1 = fp
            elseif (i2.eq.2) then
              fpt2 = fp
c       linear interpolation in the time variable
              fp = fpt1 + (fpt2-fpt1)*(t-tgr(it))/dt
              fpt = (fpt2-fpt1)/dt
            endif
90        continue

          select case (i1)
            case (1)
              bx = fp - bfac1*y(1)*y(3)
              dbxdx = fpx - bfac1*y(3)+bfac2*y(1)*y(1)*y(3)
              dbxdy = fpy + bfac2*y(1)*y(2)*y(3)
              dbxdz = fpz - bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
              dbxdt = fpt
            case (2)
              by = fp - bfac1*y(2)*y(3)
              dbydx = fpx + bfac2*y(2)*y(3)*y(1)
              dbydy = fpy - bfac1*y(3)+bfac2*y(2)*y(2)*y(3)
              dbydz = fpz - bfac1*y(2)+bfac2*y(3)*y(3)*y(2)
              dbydt = fpt
            case (3)
              bz = fp - bfac1*y(3)*y(3) + b0/r2/r
              dbzdx = fpx - bfac1*y(1)+bfac2*y(1)*y(3)*y(3)
              dbzdy = fpy - bfac1*y(2)+bfac2*y(2)*y(3)*y(3)
              dbzdz = fpz - 3.*bfac1*y(3)+bfac2*y(3)*y(3)*y(3)
              dbzdt = fpt
            case (4)
              ex = fp
            case (5)
              ey = fp
            case (6)
              ez = fp

          end select

100   continue

      b=dsqrt(bx**2.+by**2.+bz**2.)
      dbdx=(bx*dbxdx + by*dbydx + bz*dbzdx)/b
      dbdy=(bx*dbxdy + by*dbydy + bz*dbzdy)/b
      dbdz=(bx*dbxdz + by*dbydz + bz*dbzdz)/b
      dbdt=(bx*dbxdt + by*dbydt + bz*dbzdt)/b

      return
      end
      
***********************************************************************      
      
      subroutine allderiv(f,x,i1,i2,i3,i4,p)
      
      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-intmat.inc'

      real f(nx,ny,nz,nstep)
      real d1f,d2f,d3f
      integer i1,i2,i3,i4,p,x     
      
      alphamat(i1,i2,i3,i4,x,p)  = f(i1,i2,i3,i4)
      alphamat(i1,i2,i3,i4,x,p+8)  = 
     &          d1f(f(i1-1,i2,i3,i4),f(i1+1,i2,i3,i4),dx)
      alphamat(i1,i2,i3,i4,x,p+16) = 
     &          d1f(f(i1,i2-1,i3,i4),f(i1,i2+1,i3,i4),dy)
      alphamat(i1,i2,i3,i4,x,p+24) = 
     &          d1f(f(i1,i2,i3-1,i4),f(i1,i2,i3+1,i4),dz)
      alphamat(i1,i2,i3,i4,x,p+32) = 
     &          d2f(f(i1-1,i2-1,i3,i4),f(i1+1,i2-1,i3,i4),
     &          f(i1-1,i2+1,i3,i4),f(i1+1,i2+1,i3,i4),dx,dy)
      alphamat(i1,i2,i3,i4,x,p+40) = 
     &          d2f(f(i1-1,i2,i3-1,i4),f(i1+1,i2,i3-1,i4),
     &          f(i1-1,i2,i3+1,i4),f(i1+1,i2,i3+1,i4),dx,dz)
      alphamat(i1,i2,i3,i4,x,p+48) = 
     &          d2f(f(i1,i2-1,i3-1,i4),f(i1,i2+1,i3-1,i4),
     &          f(i1,i2-1,i3+1,i4),f(i1,i2+1,i3+1,i4),dy,dz)
      alphamat(i1,i2,i3,i4,x,p+56) = 
     &          d3f(f(i1-1,i2-1,i3-1,i4),f(i1+1,i2-1,i3-1,i4),
     &          f(i1-1,i2+1,i3-1,i4),f(i1+1,i2+1,i3-1,i4),
     &          f(i1-1,i2-1,i3+1,i4),f(i1+1,i2-1,i3+1,i4),
     &          f(i1-1,i2+1,i3+1,i4),f(i1+1,i2+1,i3+1,i4),dx,dy,dz)
      
      return
      end
      
      
************************************************************************      
      real function  d1f(f1,f2,h)

c     returns the first derivative 

      implicit none
      real f1,f2,h

      d1f = (f2 - f1)/(2*h)

      return
      end


      real function  d2f(f1,f2,f3,f4,h,k)

c     returns the second derivative 

      implicit none
      real f1,f2,f3,f4,h,k

      d2f = (f4 - f3 - f2 + f1)/(4*h*k)

      return
      end


      real function  d3f(f1,f2,f3,f4,f5,f6,f7,f8,h,k,m)

c     returns the third derivative 

      implicit none
      real f1,f2,f3,f4,f5,f6,f7,f8,h,k,m

      d3f = (f8 - f4 - f6 + f2 - f7 + f3 + f5 + f1)/(8*h*k*m)

      return
      end


      subroutine precalcmatrix(start,finish)

      include 'rbelt-grid.inc'
      integer start,finish

      print *, 'Calculating tri-cubic interpolation matrix'
      do l = start,finish
        do k = 3,nz-2
          do j = 3,ny-2
            do i = 3,nx-2
              call allderiv(bxd,1,i,j,k,l,1)
              call allderiv(bxd,1,i+1,j,k,l,2)
              call allderiv(bxd,1,i,j+1,k,l,3)
              call allderiv(bxd,1,i+1,j+1,k,l,4)
              call allderiv(bxd,1,i,j,k+1,l,5)
              call allderiv(bxd,1,i+1,j,k+1,l,6)
              call allderiv(bxd,1,i,j+1,k+1,l,7)
              call allderiv(bxd,1,i+1,j+1,k+1,l,8)

              call allderiv(byd,2,i,j,k,l,1)
              call allderiv(byd,2,i+1,j,k,l,2)
              call allderiv(byd,2,i,j+1,k,l,3)
              call allderiv(byd,2,i+1,j+1,k,l,4)
              call allderiv(byd,2,i,j,k+1,l,5)
              call allderiv(byd,2,i+1,j,k+1,l,6)
              call allderiv(byd,2,i,j+1,k+1,l,7)
              call allderiv(byd,2,i+1,j+1,k+1,l,8)

              call allderiv(bzd,3,i,j,k,l,1)
              call allderiv(bzd,3,i+1,j,k,l,2)
              call allderiv(bzd,3,i,j+1,k,l,3)
              call allderiv(bzd,3,i+1,j+1,k,l,4)
              call allderiv(bzd,3,i,j,k+1,l,5)
              call allderiv(bzd,3,i+1,j,k+1,l,6)
              call allderiv(bzd,3,i,j+1,k+1,l,7)
              call allderiv(bzd,3,i+1,j+1,k+1,l,8)

              call allderiv(exd,4,i,j,k,l,1)
              call allderiv(exd,4,i+1,j,k,l,2)
              call allderiv(exd,4,i,j+1,k,l,3)
              call allderiv(exd,4,i+1,j+1,k,l,4)
              call allderiv(exd,4,i,j,k+1,l,5)
              call allderiv(exd,4,i+1,j,k+1,l,6)
              call allderiv(exd,4,i,j+1,k+1,l,7)
              call allderiv(exd,4,i+1,j+1,k+1,l,8)

              call allderiv(eyd,5,i,j,k,l,1)
              call allderiv(eyd,5,i+1,j,k,l,2)
              call allderiv(eyd,5,i,j+1,k,l,3)
              call allderiv(eyd,5,i+1,j+1,k,l,4)
              call allderiv(eyd,5,i,j,k+1,l,5)
              call allderiv(eyd,5,i+1,j,k+1,l,6)
              call allderiv(eyd,5,i,j+1,k+1,l,7)
              call allderiv(eyd,5,i+1,j+1,k+1,l,8)

              call allderiv(ezd,6,i,j,k,l,1)
              call allderiv(ezd,6,i+1,j,k,l,2)
              call allderiv(ezd,6,i,j+1,k,l,3)
              call allderiv(ezd,6,i+1,j+1,k,l,4)
              call allderiv(ezd,6,i,j,k+1,l,5)
              call allderiv(ezd,6,i+1,j,k+1,l,6)
              call allderiv(ezd,6,i,j+1,k+1,l,7)
              call allderiv(ezd,6,i+1,j+1,k+1,l,8)
            enddo
          enddo
        enddo
      enddo

      print *, 'Finished calculating matrix'
      return
      end
