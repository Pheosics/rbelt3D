c     rbelt-particle_loop.f - loops over each particle in distribution

************************************************************************

      subroutine particle_loop(t2)

c loops over particles in state vector y0 (declared in rbelt-dist.inc)
c passes initial state of each particle to time integration routine
c then updates y0 with the state of each particle at time t2

      implicit none
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-status.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      include 'rbelt-flag.inc'
      include 'rbelt-mu.inc'
      include 'rbelt-const.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-io.inc'
      integer i,j,nvars,firstfilenum
      real*8 t,t2,r,energy,theta,lshell
      real*8 y(6)
      external rk4,rkf4,rk2
      external rhs_lrntz,rhs_init_lrntz
      external rhs_northrop,rhs_init_northrop
      external rhs_cb,rhs_init_cb
      external boris,boris_init,dummy_init
!      print *
!      print *,'* in subroutine particle loop *'

c     loop over particles
      do i = 1,num_particles
*      do i = 1,1

c        distribute particle (dsingle particle initial conditions)
*         if (init_dist.eq..false.) call dist_particles
!        print *, t2, tzero, y0(7,i)
         t = y0(7,i)-tzero

c        if the particle has not previously gone out of bounds 
c        and t is greater than or = zero and less than t2. 
         if ((int_y0(1,i).eq.0).and.(t.lt.t2).and.(t.ge.0.0)) then

c           initialize particle
            status = int_y0(1,i)
            flag = int_y0(2,i)
            do j = 1,6
               y(j)= y0(j,i)
            enddo
!            print *, y(4),mu,y(6)
!            print *, (dsqrt(1+2*dabs(b)*mu+y(4)*y(4))-1)*wrest

c           advance particle
100         continue
            if(flag.eq.0) then
!               call time_loop(i,y,t,t2,rk4,6,rhs_lrntz,rhs_init_lrntz)
               call time_loop(i,y,t,t2,boris,6,rhs_lrntz,dummy_init)
               if ((status.eq.-2).and.(flag_switch.eqv..true.)) then
                  call lrntz2gc(y,t)
                  flag=1
               endif
            elseif (flag.eq.1) then
               mu=y(6)
               call time_loop(i,y,t,t2,rk4,4,rhs_northrop,
     &rhs_init_northrop)
               if ((status.eq.4).and.(flag_switch.eqv..true.)) then
                  call gc2lrntz(y,t)
                  flag=0
               endif
            elseif (flag.eq.2) then
               mu=y(6)
               call time_loop(i,y,t,t2,rk4,4,rhs_cb,rhs_init_cb)
               if ((status.eq.4).and.(flag_switch.eqv..true.)) then
                  call gc2lrntz(y,t)
                  flag=0
               endif
            endif
            call filewrite(i)
            if ((t.lt.t2).and.(status.le.0)) goto 100

c           if the particle hit the inner boundary send flux data to file
            if ((status.eq.1).and.(prcp_out.eqv..true.))
     &      call prcp_filewrite(i,y,t)
            if ((cone_out.eqv..true.).and.status.eq.2) 
     &      call cone_filewrite(i,status,y0(1,i),y0(2,i),y0(3,i),
     &      y0(8,i),y0(9,i))

            r=dsqrt(y(1)**2+y(2)**2+y(3)**2)
*            lshell=dsqrt(y(1)**2+y(2)**2+y(3)**2)/
*     &       ((y(1)**2+y(2)**2)/(y(1)**2+y(2)**2+y(3)**2))
*            energy=(dsqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))-1)*wrest
*            print *,'i',i
*            print *
*            if (status.eq.1) print *,'i,r,status=',i,r,status
           print *,'i,t,r,status',i,t/tfactor,r,status
*            print *,'i,t,r,status',i,t/tfactor,r,status

c           update particle
            if ((t2.lt.(tmax-tzero)).and.(status.eq.3)) then
               int_y0(1,i) = 0
            else
               int_y0(1,i) = status
            endif
            int_y0(2,i) = flag
            do j = 1,6
               y0(j,i)=y(j)
            enddo
!            print *, y(4), mu, y(6)
            y0(7,i)=t+tzero
         endif

      enddo

      return
      end

*\***********************************************************************
