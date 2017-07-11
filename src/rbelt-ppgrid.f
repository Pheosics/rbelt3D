************************************************************************

      subroutine mk1dgrd(nx,xmin,xmax,dx,x)

c     provides values for 1d array declared in pstprc.inc
c     from x(1)=xmin to x(nx)=xmax with dx=(xmax-xmin)/(nx-1)

      implicit none
      integer i,nx
      real x(nx),xmin,xmax,dx

*      print *
*      print *,'*** in mk1dgrd ***'
*      print *,'nx,xmin,xmax,dx=',nx,xmin,xmax,(xmax-xmin)/(nx-1)
      if (nx.gt.1)then
      	dx=(xmax-xmin)/(nx-1)
      else
      	dx=1
      endif
      do i=1,nx
         x(i)=xmin + (i-1)*dx
*         print *,'i,xmin,x(i),xmax=',i,x(i)-dx*.5,x(i),x(i)+dx*.5
      enddo

      return
      end

************************************************************************

      subroutine grd_delta(i,nx,x,upr,lwr,dx)

c     provides upper and lower bounds & dx of NGP region 
c     corresponding to grid point x(i) in array x(1) -> x(nx)

      implicit none
      integer i,nx
      real x(nx),upr,lwr,dx
      if ((i.lt.1).or.(i.gt.nx)) then
         print *,'(i.lt.1).or.(i.gt.nx) in subroutine grd_delta'
         stop
      endif
      if (i.eq.1)  then
         upr=x(i)+(x(i+1)-x(i))/2.
         lwr=x(i)
      elseif (i.eq.nx) then
         upr=x(i)
         lwr=x(i)-(x(i)-x(i-1))/2.
      else
         upr=x(i)+(x(i+1)-x(i))/2.
         lwr=x(i)-(x(i)-x(i-1))/2.
      endif
      dx=upr-lwr

*      if (dx.eq.0) then 
*         print *,'i,lwr,upr,dx=',i,lwr,upr,dx
*         stop
*      endif

      return
      end

************************************************************************

      subroutine pt2fld3d(num,x0,y0,z0,nx,ny,nz,x,y,z,pw,fc,f)

      implicit none 
      integer nx,ny,nz,i,j,k,ip,jp,kp,grdpos,status,num
      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz),pw,d
      real dx,dy,dz,wi,wj,wk,wip,wjp,wkp,wij,wipj,wijp,wipjp,
     &wijk,wipjk,wijpk,wijkp,wipjpk,wijpkp,wipjkp,wipjpkp
      real fc(nx,ny,nz),tmp

      status = 0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      k=grdpos(z0,nz,z,status)

      if (status.eq.1) return

      ip=i+1
      jp=j+1
      kp=k+1
      dx = x(ip)-x(i)
      dy = y(jp)-y(j)
      dz = z(kp)-z(k)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wk = (z(kp)-z0)/dz
      wip = 1.-wi
      wjp = 1.-wj
      wkp = 1.-wk
      wij = wi*wj
      wipj = wip*wj
      wijp = wi*wjp
      wipjp = wip*wjp
      wijk = wij*wk
      wipjk = wipj*wk
      wijpk = wijp*wk
      wijkp = wij*wkp
      wipjpk = wipjp*wk
      wijpkp = wijp*wkp
      wipjkp = wipj*wkp
      wipjpkp = wipjp*wkp

*      if ((i.eq.6).and.(j.eq.11).and.(k.eq.11)) then
*         print *,'k=11: num,pw=',num,pw
*      endif

*      if ((i.eq.6).and.(j.eq.11).and.(k.eq.10)) then
*         print *,'k=10: num,pw=',num,pw
*      endif

*      tmp=f(6,11,11)

      f(i,j,k)=f(i,j,k)+wijk*pw
      f(ip,j,k)=f(ip,j,k)+wipjk*pw
      f(i,jp,k)=f(i,jp,k)+wijpk*pw
      f(i,j,kp)=f(i,j,kp)+wijkp*pw
      f(ip,jp,k)=f(ip,jp,k)+wipjpk*pw
      f(i,jp,kp)=f(i,jp,kp)+wijpkp*pw
      f(ip,j,kp)=f(ip,j,kp)+wipjkp*pw
      f(ip,jp,kp)=f(ip,jp,kp)+wipjpkp*pw

      fc(i,j,k)=fc(i,j,k)+wijk
      fc(ip,j,k)=fc(ip,j,k)+wipjk
      fc(i,jp,k)=fc(i,jp,k)+wijpk
      fc(i,j,kp)=fc(i,j,kp)+wijkp
      fc(ip,jp,k)=fc(ip,jp,k)+wipjpk
      fc(i,jp,kp)=fc(i,jp,kp)+wijpkp
      fc(ip,j,kp)=fc(ip,j,kp)+wipjkp
      fc(ip,jp,kp)=fc(ip,jp,kp)+wipjpkp

*      if (f(6,11,11).ne.tmp) then
*         print *,'num,pw=',num,pw
*      endif

      return
      end

************************************************************************

      subroutine pt2fld2d(x0,y0,nx,ny,x,y,f,pw)

      implicit none 
      integer nx,ny,i,j,ip,jp,grdpos,status
      real x0,y0,x(nx),y(ny),f(nx,ny),pw
      real dx,dy,wi,wj,wip,wjp,wij,wipj,wijp,wipjp

      status=0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      if (status.eq.1) return
      ip=i+1
      jp=j+1
      dx = x(ip)-x(i)
      dy = y(jp)-y(j)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wip = 1.-wi
      wjp = 1.-wj
      wij = wi*wj
      wipj = wip*wj
      wijp = wi*wjp
      wipjp = wip*wjp

      f(i,j)=f(i,j)+wij*pw
      f(ip,j)=f(ip,j)+wipj*pw
      f(i,jp)=f(i,jp)+wijp*pw
      f(ip,jp)=f(ip,jp)+wipjp*pw

      return
      end

************************************************************************

      subroutine pt2fld1d(x0,nx,x,f,pw)

      implicit none 
      integer nx,i,ip,grdpos,status
      real x0,x(nx),f(nx),pw
      real dx,w,wp
      status=0
      i=grdpos(x0,nx,x,status)
      if (status.eq.1) return
      ip=i+1
      dx = x(ip)-x(i)
      w = (x(ip)-x0)/dx
      wp = 1.-w
      f(i)=f(i)+w*pw
      f(ip)=f(ip)+wp*pw

      return
      end

************************************************************************

      real function fld2pt3d(num,x0,y0,z0,nx,ny,nz,x,y,z,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,ny,nz,i,j,k,ip,jp,kp,grdpos,status,num
      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz)
      real dx,dy,dz,wi,wj,wk,wip,wjp,wkp,wij,wipj,wijp,wipjp,
     &wijk,wipjk,wijpk,wijkp,wipjpk,wijpkp,wipjkp,wipjpkp
      status=0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      k=grdpos(z0,nz,z,status)

*      print *
*      print *,'i,j,k=',i,j,k
*      print *,'x0,y0,z0=',x0,y0,z0
*      print *,'x(1),y(1),z(1)=',x(1),y(1),z(1)
*      print *,'x(nx),y(ny),z(nz)=',x(nx),y(ny),z(nz)

      if (status.eq.1) then
*         print *
*         print *,'point outside of grid in function grdpos'
*         print *,'particle',num
*         print *,'i,j,k=',i,j,k
*         print *,'x0,y0,z0=',x0,y0,z0
*         print *,'x(1),y(1),z(1)=',x(1),y(1),z(1)
*         print *,'x(nx),y(ny),z(nz)=',x(nx),y(ny),z(nz)
         fld2pt3d = 0.0
         return
      endif

      ip=i+1
      jp=j+1
      kp=k+1
      dx = x(ip)-x(i)
      dy = y(jp)-y(j)
      dz = z(kp)-z(k)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wk = (z(kp)-z0)/dz
      wip = 1.-wi
      wjp = 1.-wj
      wkp = 1.-wk
      wij = wi*wj
      wipj = wip*wj
      wijp = wi*wjp
      wipjp = wip*wjp
      wijk = wij*wk
      wipjk = wipj*wk
      wijpk = wijp*wk
      wijkp = wij*wkp
      wipjpk = wipjp*wk
      wijpkp = wijp*wkp
      wipjkp = wipj*wkp
      wipjpkp = wipjp*wkp

*      print *,'(i,j,k)=',i,j,k
*      print *,f(i,j,k),f(ip,j,k),f(i,jp,k),f(i,j,kp),
*     &f(ip,jp,k),f(i,jp,kp),f(ip,j,kp),f(ip,jp,kp)

      fld2pt3d = f(i,j,k)*wijk + f(ip,j,k)*wipjk + f(i,jp,k)*wijpk +
     & f(i,j,kp)*wijkp + f(ip,jp,k)*wipjpk + f(i,jp,kp)*wijpkp + 
     & f(ip,j,kp)*wipjkp + f(ip,jp,kp)*wipjpkp

      return 
      END

************************************************************************

      real function fld2pt2d(x0,y0,nx,ny,x,y,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,ny
      real x0,y0,x(nx),y(ny),f(nx,ny)
      integer i,j,ip,jp,upr,mid,grdpos,status
      real dx,dy,wi,wj,wip,wjp,wij,wipj,wijp,wipjp
      status=0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      if (status.eq.1) then
         print *,'point outside of grid in function grdpos'
         fld2pt2d = 0.0
         return
      endif
      ip=i+1
      jp=j+1
      dx=x(ip)-x(i)
      dy=y(jp)-y(j)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wip = 1.-wi
      wjp = 1.-wj
      wij =wi*wj
      wipj =wip*wj
      wijp =wi*wjp
      wipjp =wip*wjp
      fld2pt2d=wij*f(i,j)+wipj*f(ip,j)+wijp*f(i,jp)+wipjp*f(ip,jp)

      return 
      END

************************************************************************


      real function fld2pt2dt(x0,y0,tstep,nx,ny,tsteps,x,y,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,ny
      real x0,y0,x(nx),y(ny),f(nx,ny,tsteps)
      integer i,j,ip,jp,upr,mid,grdpos,status,tstep,tsteps
      real dx,dy,wi,wj,wip,wjp,wij,wipj,wijp,wipjp
      status=0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      if (status.eq.1) then
         print *,'point outside of grid in function grdpos'
         fld2pt2dt = 0.0
         return
      endif
      if ((tstep.lt.1).or.(tstep.gt.tsteps)) then
         print *,'(tstep.lt.1).or.(tstep.gt.tsteps) in fld2pt2dt'
         stop
      endif
      ip=i+1
      jp=j+1
      dx=x(ip)-x(i)
      dy=y(jp)-y(j)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wip = 1.-wi
      wjp = 1.-wj
      wij =wi*wj
      wipj =wip*wj
      wijp =wi*wjp
      wipjp =wip*wjp
      fld2pt2dt=wij*f(i,j,tstep)+wipj*f(ip,j,tstep)+wijp*f(i,jp,tstep)+
     &wipjp*f(ip,jp,tstep)

      return 
      END



************************************************************************

      subroutine pt2fld3d0(num,x0,y0,z0,nx,ny,nz,x,y,z,pw,fcnt,f)

      implicit none 
      integer nx,ny,nz,i,j,k,ip,jp,kp,grdpos0,status,num
      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz),pw,d
      real fcnt(nx,ny,nz),tmp

      status = 0
      i=grdpos0(x0,nx,x,status)
      j=grdpos0(y0,ny,y,status)
      k=grdpos0(z0,nz,z,status)
      if (status.eq.1) return

*      tmp=f(6,11,11)

      f(i,j,k)=f(i,j,k)+pw
      fcnt(i,j,k)=fcnt(i,j,k)+1.


*      if (f(6,11,11).ne.tmp) then
*         print *,'num,pw=',num,pw
*      endif

*      if ((i.eq.24).and.(j.eq.2)) then
*         print *
*         print *,'num=',num
*         print *,'i,j,k=',i,j,k
*         print *,'x0,y0,z0=',x0,y0,z0
*         print *,'nx,ny,nz=',nx,ny,nz
*         print *,'x(1),y(1),z(1)=',x(1),y(1),z(1)
*         print *,'x(nx),y(ny),z(nz)=',x(nx),y(ny),z(nz)
*         print *,'pw=',pw
*         print *,'f(i,j,k)=',f(i,j,k)
*      endif

      return
      end

************************************************************************

      subroutine pt2fld3d0t
     &(num,x0,y0,z0,tstep,nx,ny,nz,tsteps,x,y,z,pw,fcnt,f)

      implicit none 
      integer nx,ny,nz,i,j,k,grdpos0,status,num,tstep,tsteps
      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz,tsteps),pw
      real fcnt(nx,ny,nz,tsteps)

      status = 0
      i=grdpos0(x0,nx,x,status)
      j=grdpos0(y0,ny,y,status)
      k=grdpos0(z0,nz,z,status)
      if (status.eq.1) return

      f(i,j,k,tstep)=f(i,j,k,tstep)+pw
      fcnt(i,j,k,tstep)=fcnt(i,j,k,tstep)+1.

*      if ((i.eq.1).and.(j.eq.1).and.(k.eq.1).and.(tstep.eq.1)) then
*         print *,' *** num,pw/cth,fc,f=',
*     &   num,pw,fcnt(i,j,k,tstep),f(i,j,k,tstep)
*      endif

      return
      end

************************************************************************

      subroutine pt2fld2d0t(num,x0,y0,tstep,nx,ny,tsteps,x,y,pw,n,f)

      implicit none 
      integer nx,ny,i,j,grdpos0,status,num,tstep,tsteps
      real x0,y0,x(nx),y(ny),f(nx,ny,tsteps),pw
      real n(nx,ny,tsteps)

      status = 0
      i=grdpos0(x0,nx,x,status)
      j=grdpos0(y0,ny,y,status)
      if (status.eq.1) return
      if ((tstep.lt.1).or.(tstep.gt.tsteps)) then
         print *,'(tstep.lt.1).or.(tstep.gt.tsteps) in pt2fld2d0t'
         stop
      endif

      f(i,j,tstep)=f(i,j,tstep)+pw
      n(i,j,tstep)=n(i,j,tstep)+1.

*      print *,'i,j,tstep,f(i,j,tstep)=',i,j,tstep,f(i,j,tstep)

*      if ((i.eq.154).and.(j.eq.76).and.(tstep.ge.23).and.
*     &(tstep.le.26)) then
*         print *,'num=',num
*      endif

      return
      end

************************************************************************

      subroutine pt2fld3d1t
     &(num,x0,y0,z0,tstep,nx,ny,nz,tsteps,x,y,z,pw,fc,f)

      implicit none 
      integer nx,ny,nz,i,j,k,ip,jp,kp,grdpos,status,num,tstep,tsteps
      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz,tsteps),pw,d
      real dx,dy,dz,wi,wj,wk,wip,wjp,wkp,wij,wipj,wijp,wipjp,
     &wijk,wipjk,wijpk,wijkp,wipjpk,wijpkp,wipjkp,wipjpkp
      real fc(nx,ny,nz,tsteps),tmp

      status = 0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      k=grdpos(z0,nz,z,status)

      if (status.eq.1) return
      if ((tstep.gt.tsteps).or.(tstep.lt.1)) then
         print *,'(tstep.lt.1).or.(tstep.gt.tsteps) in pt2fld3d1t'
         stop
      endif

      ip=i+1
      jp=j+1
      kp=k+1
      dx = x(ip)-x(i)
      dy = y(jp)-y(j)
      dz = z(kp)-z(k)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wk = (z(kp)-z0)/dz
      wip = 1.-wi
      wjp = 1.-wj
      wkp = 1.-wk
      wij = wi*wj
      wipj = wip*wj
      wijp = wi*wjp
      wipjp = wip*wjp
      wijk = wij*wk
      wipjk = wipj*wk
      wijpk = wijp*wk
      wijkp = wij*wkp
      wipjpk = wipjp*wk
      wijpkp = wijp*wkp
      wipjkp = wipj*wkp
      wipjpkp = wipjp*wkp

*      if ((i.eq.8).and.(j.eq.2).and.(k.eq.8).and.(tstep.eq.3)) then
*         print *,'num,pw=',num,pw
*      endif

*      if ((i.eq.6).and.(j.eq.11).and.(k.eq.10)) then
*         print *,'k=10: num,pw=',num,pw
*      endif

*      tmp=f(6,11,11)

      f(i,j,k,tstep)=f(i,j,k,tstep)+wijk*pw
      f(ip,j,k,tstep)=f(ip,j,k,tstep)+wipjk*pw
      f(i,jp,k,tstep)=f(i,jp,k,tstep)+wijpk*pw
      f(i,j,kp,tstep)=f(i,j,kp,tstep)+wijkp*pw
      f(ip,jp,k,tstep)=f(ip,jp,k,tstep)+wipjpk*pw
      f(i,jp,kp,tstep)=f(i,jp,kp,tstep)+wijpkp*pw
      f(ip,j,kp,tstep)=f(ip,j,kp,tstep)+wipjkp*pw
      f(ip,jp,kp,tstep)=f(ip,jp,kp,tstep)+wipjpkp*pw

      fc(i,j,k,tstep)=fc(i,j,k,tstep)+wijk
      fc(ip,j,k,tstep)=fc(ip,j,k,tstep)+wipjk
      fc(i,jp,k,tstep)=fc(i,jp,k,tstep)+wijpk
      fc(i,j,kp,tstep)=fc(i,j,kp,tstep)+wijkp
      fc(ip,jp,k,tstep)=fc(ip,jp,k,tstep)+wipjpk
      fc(i,jp,kp,tstep)=fc(i,jp,kp,tstep)+wijpkp
      fc(ip,j,kp,tstep)=fc(ip,j,kp,tstep)+wipjkp
      fc(ip,jp,kp,tstep)=fc(ip,jp,kp,tstep)+wipjpkp

*      if (f(6,11,11).ne.tmp) then
*         print *,'num,pw=',num,pw
*      endif

      return
      end

************************************************************************

      subroutine pt2fld2d1t(num,x0,y0,tstep,nx,ny,tsteps,x,y,pw,n,f)

      implicit none 
      integer nx,ny,i,j,ip,jp,grdpos,status,num,tstep,tsteps
      real x0,y0,x(nx),y(ny),f(nx,ny,tsteps),n(nx,ny,tsteps),pw
      real dx,dy,wi,wj,wip,wjp,wij,wipj,wijp,wipjp
      status = 0
      i=grdpos(x0,nx,x,status)
      j=grdpos(y0,ny,y,status)
      if (status.eq.1) return
      if ((tstep.lt.1).or.(tstep.gt.tsteps)) then
         print *,'(tstep.lt.1).or.(tstep.gt.tsteps) in pt2fld2d1t'
         stop
      endif

      ip=i+1
      jp=j+1
      dx = x(ip)-x(i)
      dy = y(jp)-y(j)
      wi = (x(ip)-x0)/dx
      wj = (y(jp)-y0)/dy
      wip = 1.-wi
      wjp = 1.-wj
      wij = wi*wj
      wipj = wip*wj
      wijp = wi*wjp
      wipjp = wip*wjp

      f(i,j,tstep)=f(i,j,tstep)+wij*pw
      f(ip,j,tstep)=f(ip,j,tstep)+wipj*pw
      f(i,jp,tstep)=f(i,jp,tstep)+wijp*pw
      f(ip,jp,tstep)=f(ip,jp,tstep)+wipjp*pw

      n(i,j,tstep)=n(i,j,tstep)+wij
      n(ip,j,tstep)=n(ip,j,tstep)+wipj
      n(i,jp,tstep)=n(i,jp,tstep)+wijp
      n(ip,jp,tstep)=n(ip,jp,tstep)+wipjp

      return
      end

************************************************************************
*
*      subroutine pt2fld2d2t
*     &(num,x0,y0,z0,tstep,nx,ny,nz,tsteps,x,y,z,pw,fc,f)
*
*      implicit none 
*      integer nx,ny,nz,i,j,k,ip,jp,grdpos,grdpos0,status,num,tstep,
*     &tsteps
*      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz,tsteps),
*     &n(nx,ny,nz,tsteps),pw
*      real dx,dy,wi,wj,wip,wjp,wij,wipj,wijp,wipjp
*      status = 0
*      i=grdpos(x0,nx,x,status)
*      j=grdpos(y0,ny,y,status)
*c     need to put new NGP routine here
*      k=grdpos0(z0,nz,z,status)
*      if (status.eq.1) return
*      if ((tstep.lt.1).or.(tstep.gt.tsteps))  then
*         print *,'(tstep.lt.1).or.(tstep.gt.tsteps) in pt2fld2d2t'
*         stop
*      endif
*
*      ip=i+1
*      jp=j+1
*      dx = x(ip)-x(i)
*      dy = y(jp)-y(j)
*      wi = (x(ip)-x0)/dx
*      wj = (y(jp)-y0)/dy
*      wip = 1.-wi
*      wjp = 1.-wj
*      wij = wi*wj
*      wipj = wip*wj
*      wijp = wi*wjp
*      wipjp = wip*wjp
*
*      f(i,j,tstep)=f(i,j,tstep)+wij*pw
*      f(ip,j,tstep)=f(ip,j,tstep)+wipj*pw
*      f(i,jp,tstep)=f(i,jp,tstep)+wijp*pw
*      f(ip,jp,tstep)=f(ip,jp,tstep)+wipjp*pw
*
*      n(i,j,tstep)=n(i,j,tstep)+wij
*      n(ip,j,tstep)=n(ip,j,tstep)+wipj
*      n(i,jp,tstep)=n(i,jp,tstep)+wijp
*      n(ip,jp,tstep)=n(ip,jp,tstep)+wipjp
*
*      return
*      end
*
************************************************************************

c     this routine is not tested
      subroutine pt2fld1d1t(num,x0,tstep,nx,tsteps,x,pw,n,f)

      implicit none 
      integer nx,i,ip,grdpos,status,num,tstep,tsteps
      real x0,x(nx),n(nx,tsteps),f(nx,tsteps),pw
      real dx,w,wp
      status=0
      i=grdpos(x0,nx,x,status)
      if (status.eq.1) return
      ip=i+1
      dx = x(ip)-x(i)
      w = (x(ip)-x0)/dx
      wp = 1.-w
      f(i,tstep)=f(i,tstep)+w*pw
      f(ip,tstep)=f(ip,tstep)+wp*pw

      return
      end


************************************************************************

      subroutine pt2fld1d0(num,x0,nx,x,pw,fc,f)

      implicit none 
      integer nx,i,ip,grdpos0,status,num
      real x0,x(nx),fc(nx),f(nx),pw
      real dx,w,wp

      status=0
      i=grdpos0(x0,nx,x,status)
      if (status.eq.1) return

      f(i)=f(i)+pw
      fc(i)=fc(i)+1.

      return
      end

************************************************************************

      real function fld2pt3d0(num,x0,y0,z0,nx,ny,nz,x,y,z,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,ny,nz,i,j,k,ip,jp,kp,grdpos0,status,num
      real x0,y0,z0,x(nx),y(ny),z(nz),f(nx,ny,nz)
      status=0
      i=grdpos0(x0,nx,x,status)
      j=grdpos0(y0,ny,y,status)
      k=grdpos0(z0,nz,z,status)

      if (status.eq.1) then
*         print *
*         print *,'point outside of grid in function grdpos'
*         print *,'particle',num
*         print *,'i,j,k=',i,j,k
*         print *,'x0,y0,z0=',x0,y0,z0
*         print *,'x(1),y(1),z(1)=',x(1),y(1),z(1)
*         print *,'x(nx),y(ny),z(nz)=',x(nx),y(ny),z(nz)
         fld2pt3d0 = 0.0
         return
      endif

*      if (num.eq.4385) print *,'i,j,k,wf=',i,j,k,f(i,j,k)


*      if (num.eq.33456) then
*         print *
*         print *,'num=',num
*         print *,'i,j,k=',i,j,k
*         print *,'x0,y0,z0=',x0,y0,z0
*         print *,'nx,ny,nz=',nx,ny,nz
*         print *,'x(1),y(1),z(1)=',x(1),y(1),z(1)
*         print *,'x(nx),y(ny),z(nz)=',x(nx),y(ny),z(nz)
*         print *,'f(i,j,k)=',f(i,j,k)
*      endif
*      if ((i.eq.4).and.(j.eq.4).and.(k.eq.4)) then
*         print *
*         print *,'num=',num
*         print *,'i,j,k=',i,j,k
*         print *,'x0,y0,z0=',x0,y0,z0
*         print *,'nx,ny,nz=',nx,ny,nz
*         print *,'x(1),y(1),z(1)=',x(1),y(1),z(1)
*         print *,'x(nx),y(ny),z(nz)=',x(nx),y(ny),z(nz)
*         print *,'f(i,j,k)=',f(i,j,k)
*      endif

      fld2pt3d0 = f(i,j,k)

      return 
      END

************************************************************************

      real function fld2pt4d0(num,x0,y0,z0,t0,nx,ny,nz,nt,x,y,z,t,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,ny,nz,nt,i,j,k,l,ip,jp,kp,lp,grdpos0,status,num
      real x0,y0,z0,t0,x(nx),y(ny),z(nz),t(nt),f(nx,ny,nz,nt)
      status=0
      i=grdpos0(x0,nx,x,status)
      j=grdpos0(y0,ny,y,status)
      k=grdpos0(z0,nz,z,status)
      l=grdpos0(t0,nt,t,status)

      if (status.eq.1) then
*         print *
*         print *,'fld2pt4d0:'
*         print *,'point outside of grid in function grdpos'
*         print *,'particle',num
*         print *,'i,j,k,l=',i,j,k,l
*         print *,'x0,y0,z0,t0=',x0,y0,z0,t0
*         print *,'x(1),y(1),z(1),t(1)=',x(1),y(1),z(1),t(1)
*         print *,'x(nx),y(ny),z(nz),t(nt)=',x(nx),y(ny),z(nz),t(nt)
         fld2pt4d0 = 0.0
         return
      endif

      fld2pt4d0 = f(i,j,k,l)

      return 
      END

************************************************************************

      real function fld2pt1d0(num,x0,nx,x,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,i,grdpos0,status,num
      real x0,x(nx),f(nx)
      status=0
      i=grdpos0(x0,nx,x,status)

      if (status.eq.1) then
         fld2pt1d0 = 0.0
         return
      endif

      fld2pt1d0 = f(i)

      return 
      END

************************************************************************

      real function fld2pt1d(num,x0,nx,x,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,i,ip,grdpos,status,num
      real x0,x(nx),f(nx),dx,wi,wip
      status=0
      i=grdpos(x0,nx,x,status)

      if (status.eq.1) then
         fld2pt1d = 0.0
         return
      endif

      ip=i+1
      dx=x(ip)-x(i)
      wi = (x(ip)-x0)/dx
      wip = 1.-wi

      fld2pt1d = wi*f(i)+wip*f(ip)

      return 
      END

************************************************************************

      real function fld2pt2d0(num,x0,y0,nx,ny,x,y,f)

c     returns zero if point lies outside of grid

      implicit none 
      integer nx,ny,i,j,grdpos0,status,num
      real x0,y0,x(nx),y(ny),f(nx,ny)
      status=0

      i=grdpos0(x0,nx,x,status)
      print *,'status=',status

      j=grdpos0(y0,ny,y,status)
      print *,'status=',status


      if (status.eq.1) then
         fld2pt2d0 = 0.0
         return
      endif

      fld2pt2d0 = f(i,j)

      return 
      END

************************************************************************

      integer function grdpos0(x0,nx,x,status)

c     finds grid index of grid position with nearest value to x0
c     returns status=1 if x0 is outside of grid or = to x(nx)

      implicit none 
      integer nx,i,upr,mid,status
      real x0,x(nx)

      if ((x0.lt.x(1)).or.(x0.ge.x(nx))) then
*         print *,'x0 is outside of x grid in function grdpos'
         status=1
         return
      endif

      upr=nx+1
      grdpos0=0
10    continue
      if ((upr-grdpos0).gt.1) then
        mid=(upr+grdpos0)/2
        if (x0.ge.x(mid)) then
          grdpos0=mid
        else
          upr=mid
        endif
        goto 10
      endif
      if (x0.ge..5*(x(grdpos0)+x(upr))) grdpos0=grdpos0+1

      return 
      END

************************************************************************

      integer function grdpos(x0,nx,x,status)

c     finds grid index of nearest grid position at or below x0
c     returns status=1 if x0 is outside of grid or = to x(nx)

      implicit none 
      integer nx,i,upr,mid,status
      real x0,x(nx)

      if ((x0.lt.x(1)).or.(x0.ge.x(nx))) then
*         print *,'x0 is outside of x grid in function grdpos'
         status=1
         return
      endif

      upr=nx+1
      grdpos=0
10    continue
      if ((upr-grdpos).gt.1) then
        mid=(upr+grdpos)/2
        if (x0.ge.x(mid)) then
          grdpos=mid
        else
          upr=mid
        endif
        goto 10
      endif

      return 
      END

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
