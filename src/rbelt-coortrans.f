************************************************************************

      SUBROUTINE rotate(x,y,z,bx,by,bz,phi,theta)

      implicit none
      real bx,by,bz,b,bxy,cphi,sphi,ctheta,stheta,salpha
      real r11,r12,r13,r21,r22,r23,r31,r32,r33
      real x,y,z,x_new,y_new,z_new,phi,theta
      real small

      parameter(small=1.e-10)

************************************************************************
* Given vectors (x,y,z) & (bx,by,bz), calculates the rotation to b_parallel, 
* b_perp_1, b_perp_2 frame using the 1st two Euler angles phi & theta.
* Phi is a right handed rotation about original z-axis. 
* Theta is right handed rotation about the new x-axis.
* Also outputs the two rotation angles in degrees.
************************************************************************

c     Avoid divide by zero or very small number: This leaves the 
c     coordinates unchanged in cases where b is originally
c     pointing along the +z axis.
      if (abs(by).lt.small) by=-small
      b = sqrt(bx*bx+by*by+bz*bz)
      bxy = sqrt(bx*bx+by*by)

c     1st Euler angle
      cphi = -by/bxy
      sphi = bx/bxy
      phi = atan2(bx,-by)*57.29577951

c     2nd Euler angle
      ctheta = bz/b
      stheta = bxy/b
      theta = atan2(bxy,bz)*57.29577951

c     calculate rotation matrix elements
      r11=cphi
      r21=-ctheta*sphi
      r31=stheta*sphi
      r12=sphi
      r22=ctheta*cphi
      r32=-stheta*cphi
      r13=0.
      r23=stheta
      r33=ctheta

c     express x,y,z in b_parallel, b_perp_1, b_perp_2 system
      x_new = r11*x + r12*y + r13*z
      y_new = r21*x + r22*y + r23*z
      z_new = r31*x + r32*y + r33*z
      x = x_new
      y = y_new
      z = z_new

      return
      end

************************************************************************

      SUBROUTINE rotate_back(x,y,z,bx,by,bz,phi,theta)

      implicit none
      real bx,by,bz,b,bxy,cphi,sphi,ctheta,stheta,salpha
      real r11,r12,r13,r21,r22,r23,r31,r32,r33
      real x,y,z,x_new,y_new,z_new,phi,theta
      real small

      parameter(small=1.e-10)

************************************************************************
* Given vectors (x,y,z) & (bx,by,bz), calculates the rotation to b_parallel, 
* b_perp_1, b_perp_2 frame using the 1st two Euler angles phi & theta.
* Phi is a right handed rotation about original z-axis. 
* Theta is right handed rotation about the new x-axis.
* Also outputs the two rotation angles in degrees.
************************************************************************

c     Avoid divide by zero or very small number: This leaves the 
c     coordinates unchanged in cases where b is originally
c     pointing along the +z axis.
      if (abs(by).lt.small) by=-small
      b = sqrt(bx*bx+by*by+bz*bz)
      bxy = sqrt(bx*bx+by*by)

c     1st Euler angle
      cphi = -by/bxy
      sphi = bx/bxy
      phi = atan2(bx,-by)*57.29577951

c     2nd Euler angle
      ctheta = bz/b
      stheta = bxy/b
      theta = atan2(bxy,bz)*57.29577951

c     calculate rotation matrix elements
      r11=cphi
      r21=sphi
      r31=0.
      r12=-ctheta*sphi
      r22=ctheta*cphi
      r32=stheta
      r13=stheta*sphi
      r23=-stheta*cphi
      r33=ctheta

c     express x,y,z in b_parallel, b_perp_1, b_perp_2 system
      x_new = r11*x + r12*y + r13*z
      y_new = r21*x + r22*y + r23*z
      z_new = r31*x + r32*y + r33*z
      x = x_new
      y = y_new
      z = z_new

      return
      end

************************************************************************

      subroutine sphr2cart(r,pol,phi,x,y,z)
      implicit none
      real r,pol,phi
      real x,y,z
      x = r*sin(pol)*cos(phi)
      y = r*sin(pol)*sin(phi)
      z = r*cos(pol)
      return
      end

************************************************************************

      subroutine cart2sphr(r,pol,phi,x,y,z)
      implicit none
      include 'rbelt-const.inc'
      real r,pol,phi
      real x,y,z

      r=sqrt(x*x+y*y+z*z)
      pol=acos(z/r)
      if (y.gt.0) then
         phi=atan2(y,x)
      else
         phi=2.*pi+atan2(y,x)
      endif

      return
      end
