c     rbelt-particle_loop.f - loops over each particle in distribution

************************************************************************

      subroutine particle_loop(firstfilenum,t2)

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
      real t,t2,y(6),r,energy,theta,lshell
      external rk4,rkf4,rk2,rhs_lrntz,rhs_init_lrntz
      external rhs_northrop,rhs_init_northrop
      external rhs_cb,rhs_init_cb
      external boris,boris_init,dummy_init
      print *
      print *,'*** in subroutine particle loop ***'

c     loop over particles
      do i = 1,num_particles
*      do i = 10403,10403

         t = y0(7,i)-tzero

c        if the particle has not previously gone out of bounds 
c        and t is greater than or = zero and less than t2. 
!         print *,'int_y0(1,i),int_y0(2,i),t,t2=',
!     &   int_y0(1,i),int_y0(2,i),t/tfactor,t2/tfactor
         if ((int_y0(1,i).eq.0).and.(t.ge.0.0).and.(t.lt.t2)) then

c           initialize particle
            status = int_y0(1,i)
            flag = int_y0(2,i)
            do j = 1,6
               y(j)= y0(j,i)
            enddo

c           advance particle
100         continue
            if(flag.eq.0) then
               call time_loop(i,y,t,t2,boris,6,rhs_lrntz,dummy_init)
*               call time_loop(i,y,t,t2,boris,6,rhs_lrntz,dummy_init)
               if ((status.eq.-2).and.(flag_switch.eqv..true.)) then
                  call lrntz2gc(y,t)
                  flag=2
               endif
*            elseif (flag.eq.1) then
*               mu=y(6)
*               call time_loop(i,y,t,t2,rk4,4,rhs_northrop,
*     &         rhs_init_northrop)
*               if ((status.eq.4).and.(flag_switch.eqv..true.)) then
*                  call gc2lrntz(y,t)
*                  flag=0
*               endif
            elseif (flag.eq.1) then
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

            r=sqrt(y(1)**2+y(2)**2+y(3)**2)
*            lshell=sqrt(y(1)**2+y(2)**2+y(3)**2)/
*     &       ((y(1)**2+y(2)**2)/(y(1)**2+y(2)**2+y(3)**2))
*            energy=(sqrt(1+y(4)*y(4)+y(5)*y(5)+y(6)*y(6))-1)*wrest
*            if (num_particles.le.1000) 
*     &      print *,'i,t,r,flag,status',i,t/tfactor,r,flag,status
            print *,'i,t,r,flag,status',i,t/tfactor,r,flag,status

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
            y0(7,i)=t+tzero

         endif

      enddo

      return
      end

************************************************************************
