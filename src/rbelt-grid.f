c     rbelt-gid.f - rbelt space-time grid routines

***********************************************************************

      subroutine zero_tgrid(time)

c adds input time to tzero and subtracts time from all time grid steps

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      integer i
      real time

      print *
      print *,'*** in zero_tgrid ***'

*      print *,'time,tzero=',time/tfactor,tzero/tfactor
      tzero = tzero + time
*      print *,'new tzero=',tzero/tfactor

*      print *,'i, global time, local time:'
c     global time is seconds from UT ref. tme
c     time from start is seconds from tgr(1) of the initial (first) time grid
c     (this is used mainly for adjusting input params that ref. run start time)
c     local time is seconds from curent tgr(1)
      print *,'i; seconds from UT reference time; sec. from start;
     & sec. from start of current time grid'
      do i=1,nt
         tgr(i)=tgr(i)-time
*         print *,i,(tgr(i)+tzero)/tfactor,tgr(i)/tfactor
         print *,i,real((tgr(i)+tzero)/tfactor),
     &   real((tgr(i)+tzero-tzero1)/tfactor),real(tgr(i)/tfactor)
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
*     &nt-1,tgr(nt-1)/tfactor,bxdv((nt-2)*nxyz+1)
*      print *,'i,tgr(i),bxdv=',
*     &nt,tgr(nt)/tfactor,bxdv((nt-1)*nxyz+1)

      tgr(1)=tgr(nt-1)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,1)=bxd(i,j,k,nt-1)
               byd(i,j,k,1)=byd(i,j,k,nt-1)
               bzd(i,j,k,1)=bzd(i,j,k,nt-1)
               exd(i,j,k,1)=exd(i,j,k,nt-1)
               eyd(i,j,k,1)=eyd(i,j,k,nt-1)
               ezd(i,j,k,1)=ezd(i,j,k,nt-1)

               dbxdxd(i,j,k,1)=dbxdxd(i,j,k,nt-1)
               dbydxd(i,j,k,1)=dbydxd(i,j,k,nt-1)
               dbzdxd(i,j,k,1)=dbzdxd(i,j,k,nt-1)
               dbxdyd(i,j,k,1)=dbxdyd(i,j,k,nt-1)
               dbydyd(i,j,k,1)=dbydyd(i,j,k,nt-1)
               dbzdyd(i,j,k,1)=dbzdyd(i,j,k,nt-1)
               dbxdzd(i,j,k,1)=dbxdzd(i,j,k,nt-1)
               dbydzd(i,j,k,1)=dbydzd(i,j,k,nt-1)
               dbzdzd(i,j,k,1)=dbzdzd(i,j,k,nt-1)
               dbxdtd(i,j,k,1)=dbxdtd(i,j,k,nt-1)
               dbydtd(i,j,k,1)=dbydtd(i,j,k,nt-1)
               dbzdtd(i,j,k,1)=dbzdtd(i,j,k,nt-1)

            enddo
         enddo
      enddo


      tgr(2)=tgr(nt)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,2)=bxd(i,j,k,nt)
               byd(i,j,k,2)=byd(i,j,k,nt)
               bzd(i,j,k,2)=bzd(i,j,k,nt)
               exd(i,j,k,2)=exd(i,j,k,nt)
               eyd(i,j,k,2)=eyd(i,j,k,nt)
               ezd(i,j,k,2)=ezd(i,j,k,nt)

               dbxdxd(i,j,k,2)=dbxdxd(i,j,k,nt)
               dbydxd(i,j,k,2)=dbydxd(i,j,k,nt)
               dbzdxd(i,j,k,2)=dbzdxd(i,j,k,nt)
               dbxdyd(i,j,k,2)=dbxdyd(i,j,k,nt)
               dbydyd(i,j,k,2)=dbydyd(i,j,k,nt)
               dbzdyd(i,j,k,2)=dbzdyd(i,j,k,nt)
               dbxdzd(i,j,k,2)=dbxdzd(i,j,k,nt)
               dbydzd(i,j,k,2)=dbydzd(i,j,k,nt)
               dbzdzd(i,j,k,2)=dbzdzd(i,j,k,nt)
               dbxdtd(i,j,k,2)=dbxdtd(i,j,k,nt)
               dbydtd(i,j,k,2)=dbydtd(i,j,k,nt)
               dbzdtd(i,j,k,2)=dbzdtd(i,j,k,nt)

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
*     &nt-1,tgr(nt-1)/tfactor,bxdv((nt-2)*nxyz+1)
*      print *,'i,tgr(i),bxdv=',
*     &nt,tgr(nt)/tfactor,bxdv((nt-1)*nxyz+1)

      tgr(2)=tgr(nt)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               bxd(i,j,k,2)=bxd(i,j,k,nt)
               byd(i,j,k,2)=byd(i,j,k,nt)
               bzd(i,j,k,2)=bzd(i,j,k,nt)
               exd(i,j,k,2)=exd(i,j,k,nt)
               eyd(i,j,k,2)=eyd(i,j,k,nt)
               ezd(i,j,k,2)=ezd(i,j,k,nt)
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

      print *
      print *,'*** in subroutine scale_fields ***'

*      do l=start_step,last_step
*      print *
*      do i=1,nx
*         print *,'i,x,bz,ey=',i,xgr(i),
*     &   bzd(i,(ny+1)/2,(nz+1)/2,l) + (b0/ffactor/ntg)/abs(xgr(i)**3.),
*     &   eyd(i,(ny+1)/2,(nz+1)/2,l)
*      enddo
*      enddo
*      stop

      do l=start_step,last_step
         tgr(l)=tgr(l)*tfactor
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  bxd(i,j,k,l)=bxd(i,j,k,l)*ntg*ffactor
                  byd(i,j,k,l)=byd(i,j,k,l)*ntg*ffactor
                  bzd(i,j,k,l)=bzd(i,j,k,l)*ntg*ffactor
                  exd(i,j,k,l)=0.0
                  eyd(i,j,k,l)=0.0
                  ezd(i,j,k,l)=0.0
                  bxd(i,j,k,l)=bxd(i,j,k,1)
                  byd(i,j,k,l)=byd(i,j,k,1)
                  bzd(i,j,k,l)=bzd(i,j,k,1)
!                  exd(i,j,k,l)=exd(i,j,k,l)*vmsvcm*ffactor
!                  eyd(i,j,k,l)=eyd(i,j,k,l)*vmsvcm*ffactor
!                  ezd(i,j,k,l)=ezd(i,j,k,l)*vmsvcm*ffactor
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
      real test

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
                   dbxdxd(i,j,k,l)=.5*(bxd(i+1,j,k,l)-
     &             bxd(i-1,j,k,l))/dx
                   dbydxd(i,j,k,l)=.5*(byd(i+1,j,k,l)-
     &             byd(i-1,j,k,l))/dx
                   dbzdxd(i,j,k,l)=.5*(bzd(i+1,j,k,l)-
     &             bzd(i-1,j,k,l))/dx
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
                   dbxdyd(i,j,k,l)=.5*(bxd(i,j+1,k,l)-
     &             bxd(i,j-1,k,l))/dy
                   dbydyd(i,j,k,l)=.5*(byd(i,j+1,k,l)-
     &             byd(i,j-1,k,l))/dy
                   dbzdyd(i,j,k,l)=.5*(bzd(i,j+1,k,l)-
     &             bzd(i,j-1,k,l))/dy
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
                   dbxdzd(i,j,k,l)=.5*(bxd(i,j,k+1,l)-
     &             bxd(i,j,k-1,l))/dz
                   dbydzd(i,j,k,l)=.5*(byd(i,j,k+1,l)-
     &             byd(i,j,k-1,l))/dz
                   dbzdzd(i,j,k,l)=.5*(bzd(i,j,k+1,l)-
     &             bzd(i,j,k-1,l))/dz
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
      real dt,dtx2

      do l=start_step,last_step
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  if (l.eq.1) then
                     dt=tgr(l+1)-tgr(l)
                     dbxdtd(i,j,k,l)=(bxd(i,j,k,l+1)-bxd(i,j,k,l))/dt
                     dbydtd(i,j,k,l)=(byd(i,j,k,l+1)-byd(i,j,k,l))/dt
                     dbzdtd(i,j,k,l)=(bzd(i,j,k,l+1)-bzd(i,j,k,l))/dt
                  elseif (l.eq.nt) then
                     dt=tgr(l)-tgr(l-1)
                     dbxdtd(i,j,k,l)=(bxd(i,j,k,l)-bxd(i,j,k,l-1))/dt
                     dbydtd(i,j,k,l)=(byd(i,j,k,l)-byd(i,j,k,l-1))/dt
                     dbzdtd(i,j,k,l)=(bzd(i,j,k,l)-bzd(i,j,k,l-1))/dt
                  else
                     dtx2=tgr(l+1)-tgr(l-1)
                     dbxdtd(i,j,k,l)=(bxd(i,j,k,l+1)-
     &               bxd(i,j,k,l-1))/dtx2
                     dbydtd(i,j,k,l)=(byd(i,j,k,l+1)-
     &               byd(i,j,k,l-1))/dtx2
                     dbzdtd(i,j,k,l)=(bzd(i,j,k,l+1)-
     &               bzd(i,j,k,l-1))/dtx2
                  endif

               enddo
            enddo
         enddo
      enddo

      return
      end

***********************************************************************

      subroutine calc_e(start_step,last_step)
      implicit none
      integer start_step,last_step
      return
      end

***********************************************************************
*
*      subroutine calc_e(start_step,last_step)
*
*c     code to calculate E consistant with delta_B/delta_t on rbelt grid
*c
*c     solves eqn. 2 in UKHORSKIY ET AL., Impact of ULF oscillations in solar
*c     wind dynamic pressure on the outer radiation belt electrons, GEOPHYSICAL
*c     RESEARCH LETTERS, VOL. 33, L06111, doi:10.1029/2005GL024380, 2006
*
*      implicit none
*      include 'rbelt-grid.inc'
*      include 'rbelt-const.inc'
*   
*      integer i,j,k,l,start_step,last_step,ip,jp,kp
*      real vx, vy, vz, v2, v3, ep0, ex0, ey0, rho, r0, ls
*
**      OPEN(51, FILE='test_e.txt', STATUS='OLD')
*      OPEN(52, FILE='test_dbdt.txt', STATUS='OLD')
**      OPEN(53, FILE='test_b.txt', STATUS='OLD')
*
*      do l=start_step,last_step
*         do k=1,nz
*            do j=1,ny
*               do i=1,nx
*                  exd(i,j,k,l) = 0
*                  eyd(i,j,k,l) = 0
*                  ezd(i,j,k,l) = 0
*                  do kp=1,nz
*                     do jp=1,ny
*                        do ip=1,nx
*                           if (i.eq.1) then
*                              write(52,*) dbxdtd(ip,jp,kp,1),
*     $                        dbydtd(ip,jp,kp,1), dbzdtd(ip,jp,kp,1)
*                           endif
*                           if (.not.(ip.eq.i).AND.(jp.eq.j).AND.
*     &                         (kp.eq.k)) then
*                              vx = xgr(i)-xgr(ip)
*                              vy = ygr(j)-ygr(ip)
*                              vz = zgr(k)-zgr(kp)
*                              v2 = vx*vx+vy*vy+vz*vz
*                              v3 = v2*Sqrt(v2)
*                              exd(i,j,k,l) = (dbydtd(ip,jp,kp,l)*vz-
*     &                                     dbzdtd(ip,jp,kp,l)*vy)*dx**3/
*     &                                     v3 + exd(i,j,k,l)
*                              eyd(i,j,k,l) = (dbzdtd(ip,jp,kp,l)*vx-
*     &                                     dbxdtd(ip,jp,kp,l)*vz)*dx**3/
*     &                                     v3 + eyd(i,j,k,l)
*                              ezd(i,j,k,l) = (dbxdtd(ip,jp,kp,l)*vy-
*     &                                     dbydtd(ip,jp,kp,l)*vx)*dx**3/
*     &                                     v3 + ezd(i,j,k,l) 
*                           endif
*                        enddo
*                     enddo
*                  enddo
*                  rho = sqrt(xgr(i)*xgr(i)+ygr(j)*ygr(j))
*                  r0 = sqrt(xgr(i)*xgr(i)+ygr(j)*ygr(j)+zgr(k)*zgr(k))
*                  ls = r0*r0*r0/rho**2
*                  ep0 = 30000/ls/rho/tau*dx**3
*                  ex0 = ep0*ygr(j)*vmsvcm*ffactor/rho
*                  ey0 = -ep0*xgr(i)*vmsvcm*ffactor/rho
**                  write(51,*) xgr(i), ygr(j), zgr(k)
**                  write(51,*) exd(i,j,k,l), ex0
**                  write(51,*) eyd(i,j,k,l), ey0
**                  write(51,*) ezd(i,j,k,l)
*                  print*, xgr(i), ygr(j), zgr(k)
**                  print*, exd(i,j,k,l), ex0, exd(i,j,k,l)/ex0
**                  print*, eyd(i,j,k,l), ey0, eyd(i,j,k,l)/ey0
**                  print*, ezd(i,j,k,l)
*               enddo
*            enddo
*         enddo
*      enddo
*
*      return
*      end
*
*
*
*
***********************************************************************
***********************************************************************

      subroutine linterp321(t)

c linerly interpolate between fields in time steps 2 & 3 to get fields at time t

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-const.inc'
      integer i,j,k
      real t,wi2,wi3

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

*      print *,'i, seconds from UT.txt, seconds from start, local time:'
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
      NAMELIST /bounds/ rmin,rmax,xmin,xmax,ymin,ymax,zmin,zmax,tmax

      print *
      print *,'*** in subroutine bounds_init ***'

c     read in and normalize input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,bounds)
      CLOSE (81)
      tmax = tmax*tfactor+tzero1

      print *,'rmin,rmax,tmax=',rmin,rmax,(tmax-tzero1)/tfactor

      if ((xmax.gt.xgr(nx)).or.(xmax.lt.xgr(1)).or.
     &(ymax.gt.ygr(ny)).or.(ymax.lt.ygr(1)).or.
     &(zmax.gt.zgr(nz)).or.(zmax.lt.zgr(1))) then
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

c     must always express gridded fields in SM coordinates
      sys=1

      return
      end

************************************************************************
