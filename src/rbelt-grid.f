c     rbelt-gid.f - rbelt space-time grid routines

***********************************************************************

      subroutine zero_tgrid(time)

c adds input time to tzero and subtracts time from all time grid steps

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      integer i
      real*8 time

!      print *
!      print *,'*** in zero_tgrid ***'

*      print *,'time,tzero=',time/tfactor,tzero/tfactor
      tzero = tzero + time
*      print *,'new tzero=',tzero/tfactor

*      print *,'i, global time, local time:'
c     global time is seconds from UT ref. tme
c     time from start is seconds from tgr(1) of the initial (first) time grid
c     (this is used mainly for adjusting input params that ref. run start time)
c     local time is seconds from curent tgr(1)
!      print *,'i, global time, time from start, local time:'
      do i=1,nstep
         tgr(i)=tgr(i)-time
*         print *,i,(tgr(i)+tzero)/tfactor,tgr(i)/tfactor
!         print *,i,real((tgr(i)+tzero)/tfactor),
!     &   real((tgr(i)+tzero-tzero1)/tfactor),real(tgr(i)/tfactor)
      enddo

      time = 0.0

      return
      end

***********************************************************************

      subroutine subtrct_tzero(start_step,last_step)

c subtracts tzero from time grid steps start_step to last_step

      implicit none
      include 'rbelt-grid.inc'
      integer i,start_step,last_step

*      print *
*      print *,'*** in subtrct_tzero ***'

      do i=start_step,last_step
         tgr(i)=tgr(i)-tzero
      enddo

      return
      end

***********************************************************************

      subroutine last2_2first2()

c puts last two time steps into first two time steps 

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      integer i,j,k

*      print *
*      print *,'*** in subroutine last2_2first2 ***'

*      print *,'i,tgr(i),bxdv=',
*     &nstep-1,tgr(nstep-1)/tfactor,bxdv((nstep-2)*nxyz+1)
*      print *,'i,tgr(i),bxdv=',
*     &nstep,tgr(nstep)/tfactor,bxdv((nstep-1)*nxyz+1)

      tgr(1)=tgr(nstep-1)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,1)=bxd(i,j,k,nstep-1)
               byd(i,j,k,1)=byd(i,j,k,nstep-1)
               bzd(i,j,k,1)=bzd(i,j,k,nstep-1)
               exd(i,j,k,1)=exd(i,j,k,nstep-1)
               eyd(i,j,k,1)=eyd(i,j,k,nstep-1)
               ezd(i,j,k,1)=ezd(i,j,k,nstep-1)

               dbxdxd(i,j,k,1)=dbxdxd(i,j,k,nstep-1)
               dbydxd(i,j,k,1)=dbydxd(i,j,k,nstep-1)
               dbzdxd(i,j,k,1)=dbzdxd(i,j,k,nstep-1)
               dbxdyd(i,j,k,1)=dbxdyd(i,j,k,nstep-1)
               dbydyd(i,j,k,1)=dbydyd(i,j,k,nstep-1)
               dbzdyd(i,j,k,1)=dbzdyd(i,j,k,nstep-1)
               dbxdzd(i,j,k,1)=dbxdzd(i,j,k,nstep-1)
               dbydzd(i,j,k,1)=dbydzd(i,j,k,nstep-1)
               dbzdzd(i,j,k,1)=dbzdzd(i,j,k,nstep-1)
               dbxdtd(i,j,k,1)=dbxdtd(i,j,k,nstep-1)
               dbydtd(i,j,k,1)=dbydtd(i,j,k,nstep-1)
               dbzdtd(i,j,k,1)=dbzdtd(i,j,k,nstep-1)

            enddo
         enddo
      enddo


      tgr(2)=tgr(nstep)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,2)=bxd(i,j,k,nstep)
               byd(i,j,k,2)=byd(i,j,k,nstep)
               bzd(i,j,k,2)=bzd(i,j,k,nstep)
               exd(i,j,k,2)=exd(i,j,k,nstep)
               eyd(i,j,k,2)=eyd(i,j,k,nstep)
               ezd(i,j,k,2)=ezd(i,j,k,nstep)

               dbxdxd(i,j,k,2)=dbxdxd(i,j,k,nstep)
               dbydxd(i,j,k,2)=dbydxd(i,j,k,nstep)
               dbzdxd(i,j,k,2)=dbzdxd(i,j,k,nstep)
               dbxdyd(i,j,k,2)=dbxdyd(i,j,k,nstep)
               dbydyd(i,j,k,2)=dbydyd(i,j,k,nstep)
               dbzdyd(i,j,k,2)=dbzdyd(i,j,k,nstep)
               dbxdzd(i,j,k,2)=dbxdzd(i,j,k,nstep)
               dbydzd(i,j,k,2)=dbydzd(i,j,k,nstep)
               dbzdzd(i,j,k,2)=dbzdzd(i,j,k,nstep)
               dbxdtd(i,j,k,2)=dbxdtd(i,j,k,nstep)
               dbydtd(i,j,k,2)=dbydtd(i,j,k,nstep)
               dbzdtd(i,j,k,2)=dbzdtd(i,j,k,nstep)

            enddo
         enddo
      enddo

*      print *
*      print *,'i,tgr(i),bxdv=',1,tgr(1)/tfactor,bxdv(1)
*      print *,'i,tgr(i),bxdv=',2,tgr(2)/tfactor,bxdv(nxyz+1)

      return
      end

***********************************************************************

      subroutine last2second()

c puts last time steps into second

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      integer i,j,k

*      print *
*      print *,'*** in subroutine last2second ***'

*      print *,'i,tgr(i),bxdv=',
*     &nstep-1,tgr(nstep-1)/tfactor,bxdv((nstep-2)*nxyz+1)
*      print *,'i,tgr(i),bxdv=',
*     &nstep,tgr(nstep)/tfactor,bxdv((nstep-1)*nxyz+1)

      tgr(2)=tgr(nstep)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,2)=bxd(i,j,k,nstep)
               byd(i,j,k,2)=byd(i,j,k,nstep)
               bzd(i,j,k,2)=bzd(i,j,k,nstep)
               exd(i,j,k,2)=exd(i,j,k,nstep)
               eyd(i,j,k,2)=eyd(i,j,k,nstep)
               ezd(i,j,k,2)=ezd(i,j,k,nstep)
            enddo
         enddo
      enddo

*      print *
*      print *,'i,tgr(i),bxdv=',1,tgr(1)/tfactor,bxdv(1)
*      print *,'i,tgr(i),bxdv=',2,tgr(2)/tfactor,bxdv(nxyz+1)

      return
      end

***********************************************************************

      subroutine scale_fields(start_step,last_step)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'

      integer i,j,k,l,start_step,last_step

*      print *
*      print *,'*** in subroutine scale_fields ***'

*      do l=start_step,last_step
*      do i=1,nx
*         print *,'x,ey=',xgr(i),eyd(i,(ny+1)/2,(nz+1)/2,l)
*      enddo
*      enddo

      do l=start_step,last_step
         tgr(l)=tgr(l)*tfactor
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  bxd(i,j,k,l)=bxd(i,j,k,l)*ntg*ffactor
                  byd(i,j,k,l)=byd(i,j,k,l)*ntg*ffactor
                  bzd(i,j,k,l)=bzd(i,j,k,l)*ntg*ffactor
!                  exd(i,j,k,l)=0
!                  eyd(i,j,k,l)=0
!                  ezd(i,j,k,l)=0
!                  bxd(i,j,k,l)=bxd(i,j,k,1)
!                  byd(i,j,k,l)=byd(i,j,k,1)
!                  bzd(i,j,k,l)=bzd(i,j,k,1)
!                  bxd(i,j,k,l)=0
!                  byd(i,j,k,l)=0
!                  bzd(i,j,k,l)=0
                  exd(i,j,k,l)=exd(i,j,k,l)*ntg*ffactor/tfactor
                  eyd(i,j,k,l)=eyd(i,j,k,l)*ntg*ffactor/tfactor
                  ezd(i,j,k,l)=ezd(i,j,k,l)*ntg*ffactor/tfactor
                  if (exd(i,j,k,l).ne.exd(i,j,k,l)) exd(i,j,k,l) = 0
                  if (eyd(i,j,k,l).ne.eyd(i,j,k,l)) eyd(i,j,k,l) = 0
                  if (ezd(i,j,k,l).ne.ezd(i,j,k,l)) ezd(i,j,k,l) = 0
               enddo
            enddo
         enddo
      enddo

      return
      end

***********************************************************************

      subroutine calc_derivs(start_step,last_step)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
   
      integer i,j,k,l,start_step,last_step
      real*8 test

      do l=start_step,last_step
         do k=1,nz
            do j=1,ny
               do i=1,nx

                  if (i.eq.1) then
                   dbxdxd(i,j,k,l)=(bxd(i+1,j,k,l)-bxd(i,j,k,l))/dx
                   dbydxd(i,j,k,l)=(byd(i+1,j,k,l)-byd(i,j,k,l))/dx
                   dbzdxd(i,j,k,l)=(bzd(i+1,j,k,l)-bzd(i,j,k,l))/dx
                  elseif (i.eq.nx) then
                   dbxdxd(i,j,k,l)=(bxd(i,j,k,l)-bxd(i-1,j,k,l))/dx
                   dbydxd(i,j,k,l)=(byd(i,j,k,l)-byd(i-1,j,k,l))/dx
                   dbzdxd(i,j,k,l)=(bzd(i,j,k,l)-bzd(i-1,j,k,l))/dx
                  else
*                     test=(bxd(i+1,j,k,l)-bxd(i,j,k,l))*(bxd(i,j,k,l)-
*     &               bxd(i-1,j,k,l))
*                     if (test.lt.0.) then
*                        dbxdxd(i,j,k,l)=0.
*                     else
                        dbxdxd(i,j,k,l)=.5*(bxd(i+1,j,k,l)-
     &                  bxd(i-1,j,k,l))/dx
*                     endif
*                     test=(byd(i+1,j,k,l)-byd(i,j,k,l))*(byd(i,j,k,l)-
*     &               byd(i-1,j,k,l))
*                     if (test.lt.0.) then
*                        dbydxd(i,j,k,l)=0.
*                     else
                        dbydxd(i,j,k,l)=.5*(byd(i+1,j,k,l)-
     &                  byd(i-1,j,k,l))/dx
*                     endif
*                     test=(bzd(i+1,j,k,l)-bzd(i,j,k,l))*(bzd(i,j,k,l)-
*     &               bzd(i-1,j,k,l))
*                     if (test.lt.0.) then
*                        dbzdxd(i,j,k,l)=0.
*                     else
                        dbzdxd(i,j,k,l)=.5*(bzd(i+1,j,k,l)-
     &                  bzd(i-1,j,k,l))/dx
*                     endif
                  endif

                  if (j.eq.1) then
                   dbxdyd(i,j,k,l)=(bxd(i,j+1,k,l)-bxd(i,j,k,l))/dy
                   dbydyd(i,j,k,l)=(byd(i,j+1,k,l)-byd(i,j,k,l))/dy
                   dbzdyd(i,j,k,l)=(bzd(i,j+1,k,l)-bzd(i,j,k,l))/dy
                  elseif (j.eq.ny) then
                   dbxdyd(i,j,k,l)=(bxd(i,j,k,l)-bxd(i,j-1,k,l))/dy
                   dbydyd(i,j,k,l)=(byd(i,j,k,l)-byd(i,j-1,k,l))/dy
                   dbzdyd(i,j,k,l)=(bzd(i,j,k,l)-bzd(i,j-1,k,l))/dy
                  else
*                     test=(bxd(i,j+1,k,l)-bxd(i,j,k,l))*(bxd(i,j,k,l)-
*     &               bxd(i,j-1,k,l))
*                     if (test.lt.0.) then
*                        dbxdyd(i,j,k,l)=0.
*                     else
                        dbxdyd(i,j,k,l)=.5*(bxd(i,j+1,k,l)-
     &                  bxd(i,j-1,k,l))/dy
*                     endif
*                     test=(byd(i,j+1,k,l)-byd(i,j,k,l))*(byd(i,j,k,l)-
*     &               byd(i,j-1,k,l))
*                     if (test.lt.0.) then
*                        dbydyd(i,j,k,l)=0.
*                     else
                        dbydyd(i,j,k,l)=.5*(byd(i,j+1,k,l)-
     &                  byd(i,j-1,k,l))/dy
*                     endif
*                     test=(bzd(i,j+1,k,l)-bzd(i,j,k,l))*(bzd(i,j,k,l)-
*     &               bzd(i,j-1,k,l))
*                     if (test.lt.0.) then
*                        dbzdyd(i,j,k,l)=0.
*                     else
                        dbzdyd(i,j,k,l)=.5*(bzd(i,j+1,k,l)-
     &                  bzd(i,j-1,k,l))/dy
*                     endif
                  endif

                  if (k.eq.1) then
                   dbxdzd(i,j,k,l)=(bxd(i,j,k+1,l)-bxd(i,j,k,l))/dz
                   dbydzd(i,j,k,l)=(byd(i,j,k+1,l)-byd(i,j,k,l))/dz
                   dbzdzd(i,j,k,l)=(bzd(i,j,k+1,l)-bzd(i,j,k,l))/dz
                  elseif (k.eq.nz) then
                   dbxdzd(i,j,k,l)=(bxd(i,j,k,l)-bxd(i,j,k-1,l))/dz
                   dbydzd(i,j,k,l)=(byd(i,j,k,l)-byd(i,j,k-1,l))/dz
                   dbzdzd(i,j,k,l)=(bzd(i,j,k,l)-bzd(i,j,k-1,l))/dz
                  else
*                     test=(bxd(i,j,k+1,l)-bxd(i,j,k,l))*(bxd(i,j,k,l)-
*     &               bxd(i,j,k-1,l))
*                     if (test.lt.0.) then
*                        dbxdzd(i,j,k,l)=0.
*                     else
                        dbxdzd(i,j,k,l)=.5*(bxd(i,j,k+1,l)-
     &                  bxd(i,j,k-1,l))/dz
*                     endif
*                     test=(byd(i,j,k+1,l)-byd(i,j,k,l))*(byd(i,j,k,l)-
*     &               byd(i,j,k-1,l))
*                     if (test.lt.0.) then
*                        dbydzd(i,j,k,l)=0.
*                     else
                        dbydzd(i,j,k,l)=.5*(byd(i,j,k+1,l)-
     &                  byd(i,j,k-1,l))/dz
*                     endif
*                     test=(bzd(i,j,k+1,l)-bzd(i,j,k,l))*(bzd(i,j,k,l)-
*     &               bzd(i,j,k-1,l))
*                     if (test.lt.0.) then
*                        dbzdzd(i,j,k,l)=0.
*                     else
                        dbzdzd(i,j,k,l)=.5*(bzd(i,j,k+1,l)-
     &                  bzd(i,j,k-1,l))/dz
*                     endif
                  endif

               enddo
            enddo
         enddo
      enddo

      return
      end

***********************************************************************

      subroutine calc_ddt(start_step,last_step)

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
   
      integer i,j,k,l,start_step,last_step
      real*8 dt

      do l=start_step,last_step
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  if (l.eq.1) then
                     dt=tgr(l+1)-tgr(l)
                     dbxdtd(i,j,k,l)=(bxd(i,j,k,l+1)-bxd(i,j,k,l))/dt
                     dbydtd(i,j,k,l)=(byd(i,j,k,l+1)-byd(i,j,k,l))/dt
                     dbzdtd(i,j,k,l)=(bzd(i,j,k,l+1)-bzd(i,j,k,l))/dt
                  elseif (l.eq.nstep) then
                     dt=tgr(l)-tgr(l-1)
                     dbxdtd(i,j,k,l)=(bxd(i,j,k,l)-bxd(i,j,k,l-1))/dt
                     dbydtd(i,j,k,l)=(byd(i,j,k,l)-byd(i,j,k,l-1))/dt
                     dbzdtd(i,j,k,l)=(bzd(i,j,k,l)-bzd(i,j,k,l-1))/dt
                  else
                     dt=tgr(l+1)-tgr(l-1)
                     dbxdtd(i,j,k,l)=.5*(bxd(i,j,k,l+1)-
     &               bxd(i,j,k,l-1))/dt
                     dbydtd(i,j,k,l)=.5*(byd(i,j,k,l+1)-
     &               byd(i,j,k,l-1))/dt
                     dbzdtd(i,j,k,l)=.5*(bzd(i,j,k,l+1)-
     &               bzd(i,j,k,l-1))/dt
                  endif

               enddo
            enddo
         enddo
      enddo

      return
      end


***********************************************************************

      subroutine linterp321(t)

c linerly interpolate between fields in time steps 2 & 3 to get fields at time t

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      integer i,j,k
      real*8 t,wi2,wi3

      if ((t.lt.tgr(2)).or.(t.gt.tgr(3)))then
         print *,'in subroutine linterp321'
         print *,'(t.lt.tgr(2)).or.(t.gt.tgr(3))'
            print *,'t,tgr(2),tgr(3)=',
     &      t/tfactor,tgr(2)/tfactor,tgr(3)/tfactor
         stop
      endif

      wi2=(tgr(3)-t)/(tgr(3)-tgr(2))
      wi3=1-wi2
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,1)=wi2*bxd(i,j,k,2)+wi3*bxd(i,j,k,3)
               byd(i,j,k,1)=wi2*byd(i,j,k,2)+wi3*byd(i,j,k,3)
               bzd(i,j,k,1)=wi2*bzd(i,j,k,2)+wi3*bzd(i,j,k,3)
               exd(i,j,k,1)=wi2*exd(i,j,k,2)+wi3*exd(i,j,k,3)
               eyd(i,j,k,1)=wi2*eyd(i,j,k,2)+wi3*eyd(i,j,k,3)
               ezd(i,j,k,1)=wi2*ezd(i,j,k,2)+wi3*ezd(i,j,k,3)
            enddo
         enddo
      enddo
      tgr(1)=t

c     finish time and space grid & boundaries initialization
c     ideally these should not be called in utility.f
c     making these general routines as independent as possible
c     perhaps we should move this one to rbelt-main_statfld_3.f

c    grid_init sets tzero=0.0, tzero1=tgr(1)
      call grid_init()

c    bounds_init sets tmax = tmax*tfactor+tzero1
      call bounds_init()

c     note, input t gets modified to zero here
      t=t-tgr(1)
      tmin=tgr(1)

c     zero_tgrid for 1st time step only
c     zero_tgrid adds input time to tzero and subtracts time from all time grid steps
c     also sets input var to zero
      tzero = tzero + tmin
      tgr(1)=tgr(1)-tmin
      tmin = 0.0

*      print *,'i, global time, time from start, local time:'
*      print *,1,real((tgr(1)+tzero)/tfactor),
*     &real((tgr(1)+tzero-tzero1)/tfactor),real(tgr(1)/tfactor)

c     in static field snapshots we want tgrmax to be the same as tmax
c     (it just seves as a tmax in the linterp routine)
      tgrmax=tmax-tzero

      return
      end

***********************************************************************

      subroutine bounds_init()

c simulation region boundaries initialization

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      NAMELIST /bounds/ rmin,rmax,zmin,zmax,tmax

!      print *
!      print *,'*** in subroutine bounds_init ***'

c     read in and normalize input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,bounds)
      CLOSE (81)
      tmax = tmax*tfactor+tzero1

!      print *,'rmin,rmax,tmax=',rmin,rmax,(tmax-tzero1)/tfactor

      if ((rmax.gt.xgr(nx)).or.(rmax.lt.xgr(1)).or.
     & (rmax.gt.ygr(ny)).or.(rmax.lt.ygr(1)).or.
     & (zmax.gt.zgr(nz)).or.(zmax.lt.zgr(1))) then
*     & (rmax.gt.zgr(nz)).or.(rmax.lt.zgr(1))) then
         print *,'rmax=',rmax,' is outside of grid'
         print *,'nx,xmin,xmax=',nx,xgr(1),xgr(nx)
         print *,'ymin,ymax=',ygr(1),ygr(ny)
         print *,'zmin,zmax=',zgr(1),zgr(nz)
         stop
      endif

      return
      end
      
***********************************************************************

      subroutine grid_init()

c time and space grid initialization

      implicit none
      include 'rbelt-grid.inc'

*      print *
*      print *,'*** in subroutine grid_init ***'
      tzero=0.0
      tzero1=tgr(1)

c     find the grid range and spacing
      dx = xgr(2) - xgr(1)
      dy = ygr(2) - ygr(1)
      dz = zgr(2) - zgr(1)

      call set_sys()

      return
      end

************************************************************************

c     must always express gridded fields in SM coordinates

      subroutine set_sys()
      include 'rbelt-grid.inc'
!      print *
!      print *,'*** in subroutine set_sys ***'
      sys=1
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
      b=dsqrt(bx*bx+by*by+bz*bz)
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
