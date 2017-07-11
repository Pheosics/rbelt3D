c    rbelt-dist_sphr.f - launches particles isotropically into a sphere 
c    at (input parameter) radius (RE) 

************************************************************************

      subroutine dist_particles(t2)

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      real t2

      integer num_kebins,num_cabins,num_tbins
*      parameter(num_kebins=14,num_cabins=23)
      parameter(num_kebins=1,num_cabins=23,num_tbins=10)

      real camax,camin,tbmax,tbmin
      integer ke0bin(num_kebins),ca0bin(num_cabins),t0bin(num_tbins)
      real kebwidth,cabwidth,tbwidth,lbound,ubound

      integer i,j,dist_seed
      real norm,nis,diff_flux
      real rho,phi,x1,x2,x3
      real calpha,salpha,beta
      real pmax,pmin,p_perp,p_parl,alpha
      real gamma,p,energy,fweightf
      real bprime,rand,v_perp,v_prll
      real exitx,exity,exitz
      real ran0,x(3)

      NAMELIST/dist/dist_seed0,dt_dist,init_t,emin,emax,epa_min,epa_max,
     &lmin,lmax,radius,exp,factor,flag0,flag_switch

      print *
      print *,'*** in subroutine dist_particles ***'

      initdist=3
c     read in and normalize input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,dist)
      CLOSE (81)
c     normalize
      emin = emin/wrest
      emax = emax/wrest
      epa_min = epa_min/raddeg
      epa_max = epa_max/raddeg
      dt_dist = dt_dist*tfactor
      init_t = init_t*tfactor
      dist_seed=dist_seed0

      print *,'num_particles,dist_seed=',num_particles,dist_seed

      call srand(dist_seed)

      kebwidth=(emax*wrest-emin*wrest)/num_kebins
      do i = 1,num_kebins
         ke0bin(i) = 0
      enddo

      camin=-1.0
      camax=1.0
      cabwidth=(camax-camin)/num_cabins
      do i = 1,num_cabins
         ca0bin(i) = 0
      enddo

      if (nt.gt.1) then 
         tbmin=tgr(1)
         tbmax=tgr(nt-1)
         tbwidth=(tbmax-tbmin)/num_tbins
         do i = 1,num_tbins
            t0bin(i) = 0
         enddo
      endif

      do i = 1,num_particles

c        select energy
*        energy = emin
*        gamma = 1+emin
c        total energy dist is f=norm*energy^-exp
c        norm = (1.-exp)*num_particles/(emax**(1.-exp)-emin**(1.-exp))
         energy = ((emax**(1.0-exp)-emin**(1.0-exp))*(ran0(dist_seed))+
     &    emin**(1.0-exp))**(1.0/(1.0-exp))
         gamma = 1+energy
         p=sqrt(gamma*gamma-1)

c        bin energies
         do j = 1,num_kebins
            lbound=emin*wrest+(j-1)*kebwidth
            ubound=emin*wrest+j*kebwidth
            if ((energy*wrest.ge.lbound).and.
     &      (energy*wrest.lt.ubound))ke0bin(j) = ke0bin(j)+1
            if (lbound.eq.ubound)ke0bin(j) = ke0bin(j)+1
         enddo

c        select position on p & v space sphere (v*gamma = p in scaled var.)
c        selecting reversed direction when charge_sign = -1 
         calpha = 2.0*ran0(dist_seed)-1.0
         salpha = sqrt(1.0-calpha*calpha)
         beta = 2*pi*ran0(dist_seed)
         y0(4,i)=p*salpha*cos(beta)
         y0(5,i)=p*salpha*sin(beta)
         y0(6,i)=p*calpha

c        bin cos(pitch angles)
         do j = 1,num_cabins
            lbound=camin+(j-1)*cabwidth
            ubound=camin+j*cabwidth
            if ((calpha.ge.lbound).and.
     &      (calpha.lt.ubound))ca0bin(j) = ca0bin(j)+1
            if (lbound.eq.ubound)ca0bin(j) = ca0bin(j)+1
         enddo

c        to launch along z axis
*         y0(4,i)=0.0
*         y0(5,i)=0.0
*         y0(6,i)=p

*         y0(4,i)=0.0
*         y0(5,i)=p
*         y0(6,i)=0.0

c        isotropically select point on disk in perp1-perp2 plane 
c        in v_parallel,v_perp1,v_perp2 coordinate system 
         rho=radius*ran0(dist_seed)**.5
         phi=2*pi*ran0(dist_seed)
         x1=rho*cos(phi)
         x2=rho*sin(phi)
         x3=0.0

c        rotate to x-y-z coordinate system 
         call rotate(x1,x2,x3,y0(4,i),y0(5,i),y0(6,i))

c        calc. surface launch point
c        find first (tail) intersection of velocity vector through
c        x1,x2,x3 with sphere

         y0(1,i)=x1 - (y0(4,i)*(x1*y0(4,i)+x2*y0(5,i)+x3*y0(6,i)+ 
     -       sqrt(4*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i))**2 - 
     -          4*(x1**2 + x2**2 + x3**2 - radius**2)*
     -           (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2))/2.))/
     -   (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2)


         y0(2,i)=x2 - (y0(5,i)*(x1*y0(4,i)+x2*y0(5,i)+x3*y0(6,i)+ 
     -       sqrt(4*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i))**2 - 
     -          4*(x1**2 + x2**2 + x3**2 - radius**2)*
     -           (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2))/2.))/
     -   (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2)


         y0(3,i)=x3 - (y0(6,i)*(x1*y0(4,i)+x2*y0(5,i)+x3*y0(6,i)+ 
     -       sqrt(4*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i))**2 - 
     -          4*(x1**2 + x2**2 + x3**2 - radius**2)*
     -           (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2))/2.))/
     -   (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2)


*         print *,real(y0(1,i)/radius),real(y0(2,i)/radius),
*     &   real(y0(3,i)/radius)

c        to launch along z axis
*         y0(1,i)=0.0
*         y0(2,i)=0.0
*         y0(3,i)=-radius

*         y0(1,i)=0.0
*         y0(2,i)=0.0
*         y0(3,i)=0.0

*c        calc. no field (stright trajectory) exit point
*         exitx=x1 - (y0(4,i)*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i) - 
*     -       sqrt(4*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i))**2 - 
*     -          4*(x1**2 + x2**2 + x3**2 - radius**2)*
*     -           (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2))/2.))/
*     -   (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2)
*         exity=x2 - (y0(5,i)*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i) - 
*     -       sqrt(4*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i))**2 - 
*     -          4*(x1**2 + x2**2 + x3**2 - radius**2)*
*     -           (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2))/2.))/
*     -   (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2)
*         exitz=x3 - (y0(6,i)*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i) - 
*     -       sqrt(4*(x1*y0(4,i) + x2*y0(5,i) + x3*y0(6,i))**2 - 
*     -          4*(x1**2 + x2**2 + x3**2 - radius**2)*
*     -           (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2))/2.))/
*     -   (y0(4,i)**2 + y0(5,i)**2 + y0(6,i)**2)

         y0(7,i) = tzero1+init_t+ran0(dist_seed)*dt_dist
c        use to start all particles at t=zero
*         y0(7,i) = tzero1

c        bin times
         if (nt.gt.1) then 
            do j = 1,num_tbins
               lbound=tbmin+(j-1)*tbwidth
               ubound=tbmin+j*tbwidth
               if (((y0(7,i)-tzero1).ge.lbound).and.
     &         ((y0(7,i)-tzero1).lt.ubound))t0bin(j) = t0bin(j)+1
               if (lbound.eq.ubound)t0bin(j) = t0bin(j)+1
            enddo
         endif


         int_y0(1,i) = 99
         int_y0(2,i) = 0


c        initial quantities to weight particles in terms of go here
c        since these quantities are not used anywhere in the code (until filewrite)
c        they are un-normalized here
         y0(8,i)=energy*wrest
         y0(9,i)=0.0*raddeg
         y0(10,i)=y0(7,i)/tfactor

*         print *,'y0=',y0

      enddo

c     calc. flux weight
*      fweightf = 4*pi*pi*radius*radius*Re*Re*(dt_dist/tfactor)*
*     &(emax-emin)*wrest/num_particles
*      print *,'flux weight factor (cm^2 sec str MeV) =',fweightf
*      print *,'i.e. flux weight = factor*(SW input spectrum f(e))'

*      norm = (1.-exp)*num_particles/(emax**(1.-exp)-emin**(1.-exp))
*      norm = (1.-exp)*num_particles/(pmax**(1.-exp)-pmin**(1.-exp))
*      energy=5.0/wrest
*		print *,' Calc. period (Hamlin etal., 1961) =',4*(1+energy)*6.0*
*     &(1.30-0.56*1.0)/sqrt(2*energy+energy*energy)
*		print *,' Calc. period (Hamlin etal., 1961) =',4*(1+energy)*6.0*
*     &(1.30-0.56*1.0)/sqrt(2*energy+energy*energy)/
*     &tfactor,'(sec.)'
*		print *,' Calc. period (Schultz & Lanzerrotti, 1974) =',
*     &4*(1+energy)*6.0*
*     &(1.3802-0.3198*(1.0+sqrt(1.0)))/
*     &sqrt(2*energy+energy*energy)
*		print *,' Calc. period (Schultz & Lanzerrotti, 1974) =',
*     &4*(1+energy)*6.0*
*     &(1.3802-0.3198*(1.0+sqrt(1.0)))/
*     &sqrt(2*energy+energy*energy)/tfactor,'(sec.)'

c     calc. no. in sphere 
c     for a single energy
*      nis=num_particles*4*radius/3/(p/gamma)/(dt_dist_t/tfactor)/tfactor
      nis=num_particles*4*radius/3/(p/gamma)/(dt_dist/tfactor)/tfactor
      print *
      print *,'calc. mean no. in sphere=',nis

*      print *
*      print *,'flux into sphere (num/time)=',
*     &num_particles/(dt_dist/tfactor)
c     for a flat energy dist
*      nis= 4*radius*norm*(sqrt(emax**2.+2.*emax)-
*     &sqrt(emin**2.+2.*emin))/3/(dt_dist_t/tfactor)/tfactor
c     solve N[Integrate[(1 + w)w^2/Sqrt(w^2+2w),{w,.001065788,.1065788}]]
c     for exp=3,emin=1,emax=100 (=7637320.0)
*      nis = 4*radius*norm*7637320.0/3/(dt_dist_t/tfactor)/tfactor
c     for exp=2,emin=1,emax=100 (=13564.0)
*      nis = 4*radius*norm*13564.0/3/(dt_dist_t/tfactor)/tfactor

      print *,'particles entering sphere'
      print *,'N/sec.=',
     &num_particles/(dt_dist/tfactor)
      print *,'N/sec./RE^2=',
     &num_particles/4.0/pi/radius/radius/(dt_dist/tfactor)
c     this is not j
      print *,'N/sec./RE^2/4/pi=',
     &num_particles/16/pi/pi/radius/radius/(dt_dist/tfactor)

c     this is not j
      diff_flux=num_particles/4/pi/pi/radius/radius/(dt_dist/tfactor)
      print *,'j (N/sec./RE^2/str) =',diff_flux
      print *,'4*pi*j (N/sec./RE^2/str) =',4*pi*diff_flux
      if ((emax-emin).ne.0.) then
c        this is j
         print *,'j (N/sec./RE^2/str/keV) =',
     &   diff_flux/((emax-emin)*wrest*1000)
         print *,'omni dir. j (N/sec./RE^2/keV) =',
     &   4.0*pi*diff_flux/((emax-emin)*wrest*1000)
      else
         print *,'emin= emax= (MeV)',emax*wrest
      endif
      print *,'calc. number through disk /time=',
     &diff_flux*2.0*pi*pi*radius*radius*(dt_dist/tfactor)

      print *
      print *,'init. ca dist.'
      do j = 1,num_cabins
         lbound=camin+(j-1)*cabwidth
         ubound=camin+j*cabwidth
         print *,lbound,ubound,ca0bin(j)
      enddo
      print *
      print *,'init. energy dist.'
      do j = 1,num_kebins
         lbound=emin*wrest+(j-1)*kebwidth
         ubound=emin*wrest+j*kebwidth
         print *,lbound,ubound,ke0bin(j)
      enddo
      if (nt.gt.1) then 
         print *
         print *,'init. time dist.'
         do j = 1,num_tbins
            lbound=tbmin+(j-1)*tbwidth
            ubound=tbmin+j*tbwidth
            print *,lbound/tfactor,ubound/tfactor,t0bin(j)
         enddo
      endif

*      print *,'norm=',norm
*      print *,'emin,emax=',emin,emax
*      print *,'pmin,pmax=',pmin,pmax
*      print *,'wrest=',wrest
*      print *,'tfactor=',tfactor

      return
      end

************************************************************************

      subroutine dist_particles_alt()

      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-grid.inc'
      include 'rbelt-bounds.inc'
      integer num_kebins,num_cabins
*      parameter(num_kebins=14,num_cabins=23)
      parameter(num_kebins=1,num_cabins=23)

      real camax,camin
      integer ke0bin(num_kebins),ca0bin(num_cabins)
      real kebwidth,cabwidth,lbound,ubound

      integer i,j,file_num,dist_seed
      real norm,nis,diff_flux
      real rho,phi,x1,x2,x3
      real calpha,salpha,beta
      real pmax,pmin,p_perp,p_parl,alpha
      real gamma,p,energy,fweightf
      real bprime,rand,v_perp,v_prll
      real exitx,exity,exitz
      real ran0,x(3)

c     read in and normalize input parameters
      NAMELIST/dist/dist_seed0,dt_dist,init_t,emin,emax,epa_min,epa_max,
     &lmin,lmax,radius,exp,factor,flag0,flag_switch

      print *
      print *,'*** in subroutine dist_particles ***'

      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,dist)
      CLOSE (81)
c     normalize
      emin = emin/wrest
      emax = emax/wrest
      epa_min = epa_min/raddeg
      epa_max = epa_max/raddeg
      dt_dist = dt_dist*tfactor
      init_t = init_t*tfactor
      dist_seed=dist_seed0

      print *,'num_particles,dist_seed=',num_particles,dist_seed

      call srand(dist_seed)
      kebwidth=(emax*wrest-emin*wrest)/num_kebins
      do i = 1,num_kebins
         ke0bin(i) = 0
      enddo

      camin=-1.0
      camax=1.0
      cabwidth=(camax-camin)/num_cabins
      do i = 1,num_cabins
         ca0bin(i) = 0
      enddo

      do i = 1,num_particles

c        select energy
*        energy = emin
*        gamma = 1+emin
c        total energy dist is f=norm*energy^-exp
c        norm = (1.-exp)*num_particles/(emax**(1.-exp)-emin**(1.-exp))
         energy = ((emax**(1.0-exp)-emin**(1.0-exp))*(ran0(dist_seed))+
     &    emin**(1.0-exp))**(1.0/(1.0-exp))
         gamma = 1+energy
         p=sqrt(gamma*gamma-1)

c        bin energies
         do j = 1,num_kebins
            lbound=emin*wrest+(j-1)*kebwidth
            ubound=emin*wrest+j*kebwidth
            if ((energy*wrest.ge.lbound).and.
     &      (energy*wrest.lt.ubound))ke0bin(j) = ke0bin(j)+1
         enddo

c        select launch point on sphere
         calpha = 2.0*ran0(dist_seed)-1.0
         salpha = sqrt(1.0-calpha*calpha)
         beta = 2*pi*ran0(dist_seed)
         y0(1,i)=radius*salpha*cos(beta)
         y0(2,i)=radius*salpha*sin(beta)
         y0(3,i)=radius*calpha

c        select position on p & v space hemisphere (v*gamma = p in scaled var.)
c        selecting reversed direction when charge_sign = -1 

*         calpha = 2.0*ran0(dist_seed)-1.0
*         calpha = cos(pi*ran0(dist_seed)/2.0)
*         calpha = ran0(dist_seed)
*         salpha = sqrt(1.0-calpha*calpha)

         salpha = (ran0(dist_seed))**.5
         calpha = sqrt(1.0-salpha*salpha)
         beta = 2*pi*ran0(dist_seed)
         y0(4,i)=p*salpha*cos(beta)
         y0(5,i)=p*salpha*sin(beta)
         y0(6,i)=p*calpha

         x1=-y0(1,i)
         x2=-y0(2,i)
         x3=-y0(3,i)

c        rotate to x-y-z coordinate system  
         call rotate(y0(4,i),y0(5,i),y0(6,i),x1,x2,x3)

         p_parl=(y0(4,i)*bx+y0(5,i)*by+y0(6,i)*bz)/b
         p_perp=sqrt((y0(5,i)*bz-y0(6,i)*by)**2.+
     &   (y0(6,i)*bx-y0(4,i)*bz)**2.+(y0(4,i)*by-y0(5,i)*bx)**2.)/b
         calpha=cos(atan2(p_perp,p_parl))

c        bin cos(pitch angles)
         do j = 1,num_cabins
            lbound=camin+(j-1)*cabwidth
            ubound=camin+j*cabwidth
            if ((calpha.ge.lbound).and.
     &      (calpha.lt.ubound))ca0bin(j) = ca0bin(j)+1
         enddo

c        to launch along z axis
*         y0(4,i)=0.0
*         y0(5,i)=0.0
*         y0(6,i)=p
*         y0(1,i)=0.0
*         y0(2,i)=0.0
*         y0(3,i)=-radius

         y0(7,i) = tzero1+init_t+ran0(dist_seed)*dt_dist
c        use to start all particles at t=zero
*         y0(7,i) = tzero1

         int_y0(1,i) = 0
         int_y0(2,i) = flag0

*         x(1) = y0(1,i)
*         x(2) = y0(2,i)
*         x(3) = y0(3,i)
*         call get_fields(x,(y0(7,i)-tzero))
*         p_parl=(y0(4,i)*bx+y0(5,i)*by+y0(6,i)*bz)/b
*         p_perp=sqrt((y0(5,i)*bz-y0(6,i)*by)**2.+
*     &   (y0(6,i)*bx-y0(4,i)*bz)**2.+(y0(4,i)*by-y0(5,i)*bx)**2.)/b
*         alpha=atan2(p_perp,p_parl*charge_sign)*raddeg
*         write (18,10) y0(1,i),y0(2,i),y0(3,i),energy*wrest,alpha
*10       format (5e12.4)

      enddo

c     calc. flux weight
*      fweightf = 4*pi*pi*radius*radius*Re*Re*(dt_dist/tfactor)*
*     &(emax-emin)*wrest/num_particles
*      print *,'flux weight factor (cm^2 sec str MeV) =',fweightf
*      print *,'i.e. flux weight = factor*(SW input spectrum f(e))'

*      norm = (1.-exp)*num_particles/(emax**(1.-exp)-emin**(1.-exp))
*      norm = (1.-exp)*num_particles/(pmax**(1.-exp)-pmin**(1.-exp))
*      energy=5.0/wrest
*		print *,' Calc. period (Hamlin etal., 1961) =',4*(1+energy)*6.0*
*     &(1.30-0.56*1.0)/sqrt(2*energy+energy*energy)
*		print *,' Calc. period (Hamlin etal., 1961) =',4*(1+energy)*6.0*
*     &(1.30-0.56*1.0)/sqrt(2*energy+energy*energy)/
*     &tfactor,'(sec.)'
*		print *,' Calc. period (Schultz & Lanzerrotti, 1974) =',
*     &4*(1+energy)*6.0*
*     &(1.3802-0.3198*(1.0+sqrt(1.0)))/
*     &sqrt(2*energy+energy*energy)
*		print *,' Calc. period (Schultz & Lanzerrotti, 1974) =',
*     &4*(1+energy)*6.0*
*     &(1.3802-0.3198*(1.0+sqrt(1.0)))/
*     &sqrt(2*energy+energy*energy)/tfactor,'(sec.)'

c     calc. no. in sphere 
c     for a single energy
*      nis=num_particles*4*radius/3/(p/gamma)/(dt_dist_t/tfactor)/tfactor
      nis=num_particles*4*radius/3/(p/gamma)/(dt_dist/tfactor)/tfactor
      print *
      print *,'calc. mean no. in sphere=',nis

*      print *
*      print *,'flux into sphere (num/time)=',
*     &num_particles/(dt_dist/tfactor)
c     for a flat energy dist
*      nis= 4*radius*norm*(sqrt(emax**2.+2.*emax)-
*     &sqrt(emin**2.+2.*emin))/3/(dt_dist_t/tfactor)/tfactor
c     solve N[Integrate[(1 + w)w^2/Sqrt(w^2+2w),{w,.001065788,.1065788}]]
c     for exp=3,emin=1,emax=100 (=7637320.0)
*      nis = 4*radius*norm*7637320.0/3/(dt_dist_t/tfactor)/tfactor
c     for exp=2,emin=1,emax=100 (=13564.0)
*      nis = 4*radius*norm*13564.0/3/(dt_dist_t/tfactor)/tfactor

      print *,'particles entering sphere'
      print *,'N/sec.=',
     &num_particles/(dt_dist/tfactor)
      print *,'N/sec./RE^2=',
     &num_particles/4.0/pi/radius/radius/(dt_dist/tfactor)
c     this is not j
      print *,'N/sec./RE^2/4/pi=',
     &num_particles/16/pi/pi/radius/radius/(dt_dist/tfactor)

c     this is not j
      diff_flux=num_particles/4/pi/pi/radius/radius/(dt_dist/tfactor)
      print *,'j (N/sec./RE^2/str) =',diff_flux
      print *,'4*pi*j (N/sec./RE^2/str) =',4*pi*diff_flux
c     this is j
      print *,'j (N/sec./RE^2/str/keV) =',
     &diff_flux/((emax-emin)*wrest*1000)
      print *,'omni dir. j (N/sec./RE^2/keV) =',
     &4.0*pi*diff_flux/((emax-emin)*wrest*1000)
      print *,'calc. number through disk /time=',
     &diff_flux*2.0*pi*pi*radius*radius*(dt_dist/tfactor)

      print *
      print *,'init. ca dist.'
      do j = 1,num_cabins
         lbound=camin+(j-1)*cabwidth
         ubound=camin+j*cabwidth
         print *,lbound,ubound,ca0bin(j)
      enddo
      print *
      print *,'init. energy dist.'
      do j = 1,num_kebins
         lbound=emin*wrest+(j-1)*kebwidth
         ubound=emin*wrest+j*kebwidth
         print *,lbound,ubound,ke0bin(j)
      enddo

*      print *,'norm=',norm
*      print *,'emin,emax=',emin,emax
*      print *,'pmin,pmax=',pmin,pmax
*      print *,'wrest=',wrest
*      print *,'tfactor=',tfactor

      return
      end

************************************************************************

      subroutine rotate(x,y,z,bx,by,bz)

c     Expresses vector x,y,z in coordinate system with z-axis along bx,by,bz;
c     i.e., input x,y,z (expressed in system 1) and vector bx,by,bz defineing
c     z-axis of system 2 (in system 1 coordinates), and the routine outputs 
c     x,y,z expressed in system 2 coordinates. Note bx,by,bz are inputs, and 
c     x,y,z are inputs and outputs.

      implicit none
      real bx,by,bz,b,bxy,cphi,sphi,ctheta,stheta,cpsi,spsi
      real r11,r12,r13,r21,r22,r23,r31,r32,r33
      real x,y,z,x_new,y_new,z_new

      if (abs(by).lt.1e-8) by=1e-8 
      b = sqrt(bx*bx+by*by+bz*bz)
      bxy = sqrt(by*by+bx*bx)

      ctheta = bz/b
      stheta = -bxy/b

      r11=1.0
      r21=0.0
      r31=0.0

      r12=0.0
      r22=ctheta
      r32=stheta

      r13=0.0
      r23=-stheta
      r33=ctheta

c     RH-negative rotatation through angle theta about system 1 x axis
      x_new = r11*x + r12*y + r13*z
      y_new = r21*x + r22*y + r23*z
      z_new = r31*x + r32*y + r33*z

      cphi = by/bxy
      sphi = -bx/bxy

      r11=cphi
      r21=sphi
      r31=0.0

      r12=-sphi
      r22=cphi
      r32=0.0

      r13=0.0
      r23=0.0
      r33=1.0

c     RH-negative rotation through angle phi about system 1 z axis
      x = r11*x_new + r12*y_new + r13*z_new
      y = r21*x_new + r22*y_new + r23*z_new
      z = r31*x_new + r32*y_new + r33*z_new

      return
      end

************************************************************************

      subroutine dist_particles2(t2)

      implicit none
      include 'rbelt-const.inc'
*      include 'rbelt-dist.inc'
      include 'rbelt-y0.inc'
      include 'rbelt-fields.inc'
      include 'rbelt-grid.inc'
*      include 'rbelt-bounds.inc'
      integer i
      real t2,x(3),t,p_parl,p_perp

      do i = 1,num_particles
         t=y0(7,i)-tzero
         if ((int_y0(1,i).eq.99).and.(t.ge.tgr(1)).and.(t.lt.t2)) then
            x(1) = y0(1,i)
            x(2) = y0(2,i)
            x(3) = y0(3,i)
            call get_fields(x,t)
            p_parl=(y0(4,i)*bx+y0(5,i)*by+y0(6,i)*bz)/b
            p_perp=sqrt((y0(5,i)*bz-y0(6,i)*by)**2.+
     &      (y0(6,i)*bx-y0(4,i)*bz)**2.+(y0(4,i)*by-y0(5,i)*bx)**2.)/b
c           initial PA to weight particles in terms of
c           since this quantity not used anywhere in the code (until file-write)
c           it is un-normalized here
            y0(9,i)=atan2(p_perp,p_parl*charge_sign)*raddeg
c           note that un-normalized form is
c           alpha=atan2(p_perp,p_parl*charge_sign)*raddeg
         endif
      enddo

      return
      end

************************************************************************
*
*      subroutine ic_data2y0()
*
*c loops over particles in state vector y0 (declared in rbelt-dist.inc)
*
*      implicit none
*      include 'rbelt-dist.inc'
*      include 'rbelt-y0.inc'
*      include 'rbelt-grid.inc'
*      include 'rbelt-flag.inc'
*      include 'rbelt-const.inc'
*      include 'rbelt-fields.inc'
*      include 'rbelt-io.inc'
*      include 'rbelt-bounds.inc'
*
*      integer i,n
*      real t,x(3),p_parl,p_perp,mu,y(6),ke
*
*      print *
*      print *,'*** in subroutine ic_data2y0 ***'
*
*c     loop over particles
*      do i = 1,num_particles
*         t = y0(7,i)-tzero
*c        if the particle has not already gone out of bounds and
*c        its time is less than the last time step in the current time grid. 
*         if ((int_y0(1,i).eq.0).and.(t.lt.tgrmax).and.(t.ge.0.0)) then
*            x(1) = y0(1,i)
*            x(2) = y0(2,i)
*            x(3) = y0(3,i)
*            call get_fields(x,t)
*            if (flag .eq. 0) then
*c              energy
*               y0(8,i)=(sqrt(1+y0(4,i)*y0(4,i)+y0(5,i)*y0(5,i)+
*     &         y0(6,i)*y0(6,i))-1)*wrest
*c              pitch angle
*               p_parl=(y0(4,i)*bx+y0(5,i)*by+y0(6,i)*bz)/b
*               p_perp=sqrt((y0(5,i)*bz-y0(6,i)*by)**2.+
*     &        (y0(6,i)*bx-y0(4,i)*bz)**2.+(y0(4,i)*by-y0(5,i)*bx)**2.)/b
*               y0(9,i)=atan2(p_perp,p_parl*charge_sign)*raddeg
*            elseif (flag .eq. 1) then
*c              energy
*               mu=(ke*(ke+2.)-y0(4,n)*y0(4,n))/2/b
*               y0(8,i)=(sqrt(1.+2.*abs(b)*mu+y(4)*y(4))-1.)*wrest
*c              pitch angle
*               p_perp=sqrt(2.*abs(b)*mu)
*               yout(4,youtstep)=atan2(p_perp,y0(4,i)*charge_sign)*raddeg
*            endif
*            y0(10,i)=y0(7,i)*tfactor
*         endif
*      enddo
*      return
*      end
*
