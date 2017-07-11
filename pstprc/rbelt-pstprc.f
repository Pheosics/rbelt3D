c THIS CODE IS UNDER DEVELOPMENT
c
c***********************************************************************
c                PROGRAM MAIN FOR RBELT POST PROCESSING
c***********************************************************************
c
c weights particles based on some model and computes fluxes
c
c There are two methods used in this code:
c
c 1) PARTICLES WEIGHTED IN TERMS OF THEIR INITIAL POSITION IN L-E-PA SPACE 
c AND THE INITIAL TEST PARTICLE FLUX (flux measured in SM coordinates 
c equatorial plane -- could be changed to bmin surface):
c a) read flux counts and compute the test particle flux
c b) construct weighting function, e.g., 
c    weight(L,E,PA) = j_AE8(L,E,PA)/j_initial_test-particle(L,E,PA)
c c) read init_out data and each particle is weighted based on its 
c    initial position in L-E-PA space.
c d) re-read flux counts and compute weighted fluxes (may be done in conjunction with step 3)
c (this method is typically used with "disk" launch routine, init_out=.true.)
c
c 2) WEIGHT PARTICLES IN TERMS OF INITIAL E-PA-TIME AND MODEL INTERPLANETARY 
c DISTRIBUITON: here we read in flux counts, assign weights based on an analytic 
c function of initial E-PA-TIME, then compute weighted fluxes.
c (this method is typically used with "sphr" launch routine, init_out=.false.,
c and initial particle data saved in y0(i,8-10) and included in flux file)
c
c * code to comput PSD from dist_out is not completed
c
c * regarding problem with migration of particles with big weights 
c into regions with particles with small weights: the solution must be 
c to use a finer grid and more particles, i.e., this is just counting noise
c that is inherrent in the method.
c
c * need to stop code in cases with zero counts
c
c * In rbelt-io.f, subroutine write_init_cond, output to file must be:
c format (6e12.4) with y0(1,n),y0(2,n),y0(3,n),lshell,ke,alpha_o,y0(7,n)/tfactor


      implicit none
      include 'rbelt-pstprc.inc'
      character*80 basename, filename
      integer i,filenum,firstfilenum,lastfilenum,filenumstep,mode,j
      integer k
      NAMELIST /pstprc/ basename,firstfilenum,lastfilenum,filenumstep,
     & norm,spctrm,ecutoff,mincnt

      print *
      print *,'*** IN PROGRAM MAIN ***'

c     read in and normalize rbelt input parameters
      OPEN (UNIT=81,FILE='rbelt-input.txt',STATUS='OLD') 
      READ (81,pstprc)
      CLOSE (81)

c     get rbelt run information
      filenum=firstfilenum
      call info_fileread(basename,filenum)

c     initializations
c     get rbelt wtimes and pstprc time steps
      call get_wtimes()
c     set up distribution grids
      call grid_setup()
c     initialize distribution function arrays to zero
      call init_arrays()

      if (method.eq.1) then
c        set up particle weight function grids
         call wfgrid_setup()
c        get test particle flux at all tsteps
         mode=0
         fcount=0
         do filenum=firstfilenum,lastfilenum,filenumstep
            print *,'filenum=',filenum
            fcount=fcount+1
            call info_fileread(basename,filenum)
c           read flux counts, weight particles 
c           & calc. test particle flux distributions j_leat & j_xyt
            call flux_fileread(basename,filenum,mode)
         enddo

c        turn flux counts into fluxes (#/cm^2-s-str-MeV)
*         call wghtcnts2jxyt()
         call wghtcnts2jleat()
c        create weighting function
*         call mkwghtfn_flat()
         call mkwghtfn_VAP()
*         call mkwghtfn_esa()
*         call mkwghtfn_1()
c        read particle initial conditions and weight each particle
         fcount=0
         do filenum=firstfilenum,lastfilenum,filenumstep
            fcount=fcount+1
            call init_fileread(basename,filenum)
         enddo
c        read 3rd step & re-weight particles (not currently used)
!         fcount=0
!         mode=2
!         do filenum=firstfilenum,lastfilenum,filenumstep
!            fcount=fcount+1
!            if (flux_out.eqv..true.) then
!               call info_fileread(basename,filenum) 
!               call flux_fileread(basename,filenum,mode)
!            endif
!         enddo
c        re-initialize arrays
         call init_arrays()
      endif
c     make weighted distributions
      mode=1
      fcount=0
      do filenum=firstfilenum,lastfilenum,filenumstep
!         print *,'filenum=',filenum
         fcount=fcount+1
         call info_fileread(basename,filenum)
c        read prcp. counts & calc. weighted prcp. count dist.
         if (prcp_out.eqv..true.) then
            call prcp_fileread(basename,filenum)
            call prcp_wghtcnts2flux()
            call write_prcp()
         endif
c        read flux counts & calc. weighted flux count dist.
         if (flux_out.eqv..true.) then
            call flux_fileread(basename,filenum,mode)
         endif
*c        read dist counts & calc. weighted distribution
*         if (dist_out.eqv..true.) then
*     &      call dist_fileread(basename,filenum,mode)
*            call wghtcnts2psd()
*            call write_psd()
*         endif
      enddo
*      call wghtcnts2jxyt()
      call wghtcnts2jleat()
      open(67,file='test2.dat')
      do i = 1, nl
        do j = 1, ne
          do k = 1,na
           print *, ls(i),ke(j),ca(k)
           print *, n_leat(i,j,k,itstep),j_leat(i,j,k,itstep)          
!          print *, n_leat(i,j,na/2,itstep),j_leat(i,j,na/2,itstep),'|',
!     &      n_leat(i,j,na/2,num_tsteps-1),j_leat(i,j,na/2,num_tsteps-1)
           write(67,*)n_leat(i,j,k,initstep),n_leat(i,j,k,itstep)
          enddo
        enddo
      enddo
      close(67)
*      call write_jxy()
      call write_jlea()

      end

********************************************************************************

c     dummy routines

      subroutine grid_init()
      return
      end

********************************************************************************

      subroutine bounds_init()
      return
      end

********************************************************************************

      real function pwghtfn(pnum,e0,ca0,t0)

      implicit none
      include 'rbelt-pstprc.inc'
      integer pnum
      real e0,ca0,t0

      pwghtfn=1.0

**      print *
**      print *,'j_tp=',(pper_file*num_files/
**     & (4.0*pi*pi*radius*radius*(emax-emin)*dt_dist))
*      do j=1,na
*         do i=1,ne
*            wf(i,j)=norm*(ke(i)/1000.)**spctrm/(pper_file*num_files/
*     &       (4.0*pi*pi*radius*radius*(emax-emin)*dt_dist))
**            print *,'ke,j_obs,wght=',
**     &       ke(i),norm*(ke(i)/1000.)**spctrm,wght(i,j)
*         enddo
*      enddo

*               j_model = 4*10**8 * 2.7182818**(-4.0*a1(i)) * 
*     &         a2(j)**(-3.) * sin(acos(a3(k)))**4.

      return
      end

********************************************************************************

      subroutine mkwghtfn_1()

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k

      print *
      print *,'*** in subroutine mkwght ***'

      do k=1,n3
         do j=1,n2
            do i=1,n1
               wf(i,j,k)=1.0
            enddo
         enddo
      enddo

      return
      end


********************************************************************************

      subroutine mkwghtfn_flat()

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k,count,tstep,l
      real j_model,jprv,wfmin,wfmax,maxwf
      integer mcnt
      real nprv

      print *
      print *,'*** in subroutine mkwghtfn ***'

      if (itstep.gt.num_tsteps) then
         print *,'tstep.gt.num_tsteps in mkwghtfn_flat'
         stop
      endif

c     j_leat(i,j,k,#) must have same dimension as wf(i,j,k)
c     compute weight function
      count=0
      wfmin=1000000000.0
      wfmax=0.0
      mcnt =100
      do k=1,n3
         do j=1,n2
            do i=1,n1
                  jprv=j_leat(i,j,k,itstep)
                  nprv=n_leat(i,j,k,itstep)
               if ((jprv.gt.0.0).and.(nprv.gt.mcnt)) then
!                  wf(i,j,k)=(exp(-7*ls(i))*(ke(j)/10**6)**(-3)*
!     &                      sin(acos(ca(k)))**8)/jprv
                  wf(i,j,k)=1./jprv
                  if (wf(i,j,k).lt.wfmin) wfmin=wf(i,j,k)
                  if (wf(i,j,k).gt.wfmax) wfmax=wf(i,j,k)
               else
                  count=count+1
                  wf(i,j,k)= 0.0
               endif
            enddo
         enddo
      enddo

      print *
      print *,'count,wfmin,wfmax=',count,wfmin,wfmax


*      print *
*      do k=1,n3
*         do j=2,2
*            do i=5,5
*               print *,'k,alpha,wf=',k,acos(ca(k))*raddeg,wf(i,j,k)
*            enddo
*         enddo
*      enddo

      return
      end

*******************************************************************************

      subroutine mkwghtfn_VAP()

c     calculates particle weighting function
c     at present, we modify this routine as needed

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k
      integer ii,jj,kk
      integer is,js,ks
      real angle(25),energy(4),lshell(12)
      real jvap(12,4,25),svap(12,4,25)
      real a,b,fls,fke,fca
      real c00,c01,c10,c11,c0,c1,wf0
      real jprv,jmid
      character temp
      integer mcnt
      integer nprv

      mcnt = 100

      energy = (/26.0,34.0,42.0,58.0/)
      lshell = (/2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1/)

      print *,'opening ','pads1.txt'
      open (12,file='pads3.txt')
      read (12,*) temp
      read (12,*) temp
      do i=1,12
         read (12,*) temp
         read (12,*) temp
         read (12,*) temp
         do j=1,25
            jvap(i,1,j) = 0.
            jvap(i,2,j) = 0.
            jvap(i,3,j) = 0.
            jvap(i,4,j) = 0.
            read(12,5) angle(j),jvap(i,1,j),jvap(i,2,j),jvap(i,3,j),
     &                 jvap(i,4,j),svap(i,1,j),svap(i,2,j),
     &                 svap(i,3,j),svap(i,4,j)
         enddo
      enddo
5     format (f10.1,8e12.2)
      close(12)

        do kk=1,25
         do jj=1,4
          do ii=1,12
           if (jvap(ii,jj,kk).le.0) then
            jvap(ii,jj,kk)=jvap(ii,jj,13)*sin(angle(kk)/raddeg)**8
            if (jvap(ii,jj,kk).le.0) then
             jmid = (jvap(ii,jj,12)+jvap(ii,jj,14))/2
             jvap(ii,jj,kk)=jmid*sin(angle(kk)/raddeg)**8
            endif
           endif
          enddo
         enddo
        enddo

      open (66,file="test.dat")
      do k=n3,1,-1
       fca = -1
       do kk=1,25
        if (acos(ca(k))*raddeg.le.angle(kk)) then
         if (kk.eq.1) then
          a = acos(ca(k))*raddeg - angle(kk)
          b = angle(kk+1) - acos(ca(k))*raddeg
          fca = a/(a+b)
          ks = 2
         else
          a = acos(ca(k))*raddeg - angle(kk-1)
          b = angle(kk) - acos(ca(k))*raddeg
          fca = a/(a+b)
          ks = kk
         endif
!         print *, k, kk
!         print *, acos(ca(k))*raddeg,angle(kk)
!         print *, a,b,fca
         exit
        endif
       enddo
       do j=1,n2
        fke = -1
        do jj=1,4
         if (ke(j).le.energy(jj)) then
          if (jj.eq.1) then
           a = ke(j) - energy(jj)
           b = energy(jj+1) - ke(j)
           fke = a/(a+b)
           js = 2
          else
           a = ke(j) - energy(jj-1)
           b = energy(jj) - ke(j)
           fke = a/(a+b)
           js = jj
          endif
          exit
         endif
        enddo
c Can skip the inerpolation in L if we are using the exact same grid.
        do i=1,n1
!         fls = -1
!         do ii=1,12
!          if (ls(i).le.lshell(ii)) then
!           if (ii.eq.1) then
!            a = ls(i) - lshell(ii)
!            b = lshell(ii+1) - ls(i)
!            fls = a/(a+b)
!            is = 2
!           else
!            a = ls(i) - lshell(ii-1)                  
!            b = lshell(ii) - ls(i)
!            fls = a/(a+b)
!            is = ii
!           endif
!           exit
!          endif
!         enddo
         is = i
         fls = 0
         if ((fls.ge.0).and.(fke.ge.0).and.(fca.ge.0)) then !.and.
!     &       (jvap(is,js-1,ks-1).gt.0).and.(jvap(is,js,ks-1).gt.0).and.
!     &       (jvap(is,js-1,ks).gt.0).and.(jvap(is,js,ks).gt.0).and.
!     &       (jvap(is-1,js,ks).gt.0).and.(jvap(is-1,js-1,ks).gt.0).and.
!     &       (jvap(is-1,js,ks-1).gt.0).and.
!     &       (jvap(is-1,js-1,ks-1).gt.0)) then
          c00 = jvap(is,js-1,ks-1)
          c10 = jvap(is,js,ks-1)
          c01 = jvap(is,js-1,ks)
          c11 = jvap(is,js,ks)
c         Interpolate in L
c         Linear Interpolation
!          c00 = jvap(is-1,js-1,ks-1)*(1-fls)+jvap(is,js-1,ks-1)*fls
!          c10 = jvap(is-1,js,ks-1)*(1-fls)+jvap(is,js,ks-1)*fls
!          c01 = jvap(is-1,js-1,ks)*(1-fls)+jvap(is,js-1,ks)*fls
!          c11 = jvap(is-1,js,ks)*(1-fls)+jvap(is,js,ks)*fls
c         Logarithmic Interpolation
*          c00 = jvap(ii-1,jj-1,kk-1)**(1-fls)*jvap(ii,jj-1,kk-1)**fls 
*          c10 = jvap(ii-1,jj,kk-1)**(1-fls)*jvap(ii,jj,kk-1)**fls
*          c01 = jvap(ii-1,jj-1,kk)**(1-fls)*jvap(ii,jj-1,kk)**fls
*          c11 = jvap(ii-1,jj,kk)**(1-fls)*jvap(ii,jj,kk)**fls
c         Interpolate in E
c         Linear Interpolation
          c0 = c00*(1-fke)+c10*fke
          c1 = c01*(1-fke)+c11*fke
c         Logarithmic Interpolation
!          c0 = c00**(1-fke)*c10**fke
!          c1 = c01**(1-fke)*c11**fke
c         Interpolate in A
c         Linear Interpolation
          wf0 = c0*(1-fca)+c1*fca
c         Logarithmic Interpolation
!          wf = c0**(1-fca)*c1**fca
!         print *, ls(i),ke(j),ca(k),jvap(is,js,ks),j_leat(i,j,k,itstep)
!         print *, wf0,c1,c0,c11,c10,c01,c00,fls,fke,fca
!          print*,ls(i),lshell(is)
!          print*,ke(j),energy(js)
!          print*,acos(ca(k))*raddeg,angle(ks)
!          print*,wf0,j_leat(i,j,k,itstep)
!          print*,''
         else
          wf0 = 0.
         endif
!         if ((c00.eq.0).or.(c10.eq.0).or.(c01.eq.0).or.(c11.eq.0)) then
!          wf0 = 0.
!         endif
         write (66,*) wf0
!         print *, ls(i),ke(j),acos(ca(k))*raddeg,wf0
         jprv = j_leat(i,j,k,itstep)
         nprv = n_leat(i,j,k,itstep)
!         print*,fls,fke,fca
!         print*,i,j,k,is,js,ks
          if ((jprv.gt.0).and.(nprv.gt.mcnt)) then
           wf(i,j,k) = wf0/jprv
          else
           wf(i,j,k) = 0
          endif
*          if (i.eq.12) print *, wf(i,j,k), c00,c01,c10,c11
        enddo
       enddo
      enddo
      close(66)
!      stop
      return
      end

********************************************************************************

      subroutine mkwghtfn_esa()

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k

*c     OPEN ESA RBELT MODEL FILE AND READ GRID & j_esa
*
*      do k=1,na
*         do j=1,ne
*            do i=1,nl
*               j_esa(i,j,k)=0.0
*            enddo
*         enddo
*      enddo
*      print *
*      print *,'opening ','esamodl.txt'
*      print *,'and reading ls(i),ke(j),j_esa(i,j,1)'
*      open (12,file='esamodl.txt')
*      do i=6,nl
*         do j=1,13
*            read (12,5) ls(i),ke(j),rdummy
*c           convert from log(flux) to flux
*            j_esa(i,j,1)=10.0**rdummy
*         enddo
*      enddo
*5     format (f4.2,f7.1,e12.4)
*      close(12)
*c     flat pitch angle dist.
*      do k=2,na
*         do j=1,13
*            do i=6,nl
*               j_esa(i,j,k)=j_esa(i,j,1)
*            enddo
*         enddo
*      enddo
*
*c     complete grid setup
*      ls(1)=2.0
*      do i=2,5
*         ls(i)=ls(i-1)+0.2
*      enddo
*      do j=14,26
*         ke(j)=ke(j-1)+1000.0
*      enddo
*c     cos equatorial pitch angle grid
*      call mk1dgrd(na,camin,camax,da,ca)
*      print *
*      print *,'L,energy,cos(PA) grid'
*      print *
*      do i=1,nl
*         print *,'i,ls(i)=',i,ls(i)
*      enddo
*      print *
*      do i=1,ne
*         print *,'i,ke(i)=',i,ke(i)
*      enddo
*      print *
*      do i=1,na
*         print *,'i,ca(i)=',i,ca(i)
*      enddo
*
*c     OPEN FLUX FILES & GET INITIAL TEST PARTICLE FLUX
*
*      do k=1,na
*         do j=1,ne
*            do i=1,nl
*               j_tp(i,j,k)=0.0
*               jc_tp(i,j,k)=0.0
*            enddo
*         enddo
*      enddo
*      do m=1, num_files
*         file_num=firstfilenum+(m-1)*filenumstep
*         filename=fluxfile(basename,file_num)
*         print *
*         print *,'opening ',filename(1:lnblnk(filename))
*         print *,'and getting initial test particle flux'
*c        open input file & read num. of particles
*         if (binio.eqv..true.) then
*            open (18,file=filename(1:lnblnk(filename)),
*     &      form='unformatted')
*            read (18) nparticles
*         else
*            open (18,file=filename(1:lnblnk(filename)))
*            read (18,10) nparticles
*         endif
*10       format (i8)
*         print *,'nparticles=',nparticles
*c        read particle initial conditions
*         if (binio.eqv..true.) then
*            do i=1,nparticles
*               read (18)  x0,y0,e0,alpha,theta
*            enddo
*         else
*            do i=1,nparticles
*               read (18,20)  x0,y0,e0,alpha,theta
*            enddo
*         endif
*20       format (5e12.4)
*c        read num. of write steps
*         if (binio.eqv..true.) then
*            read (18) wsteps
*         else
*            read (18,10) wsteps
*         endif
*         print *,'wsteps=',wsteps
*c        use to skip wstep
*c        ******************
**         read (18,40) nfcnts,rdummy
**         print *,'nfcnts,rdummy=',nfcnts,rdummy
**         do i=1,nfcnts
**            read (18,50) pnum,x0,y0,e0,alpha,theta
**         enddo
*c        ******************
*c        get uni dir initial test particle flux j_tp(L,E,ea) 
*c        read number of flux counts and write times
*         if (binio.eqv..true.) then
*            read (18) nfcnts,rdummy
*         else
*            read (18,40) nfcnts,rdummy
*         endif
*40       format (i8,e12.4)
*         print *,'nfcnts,rdummy=',nfcnts,rdummy
*         do i=1,nfcnts
*c           read flux count
*            if (binio.eqv..true.) then
*               read (18) pnum,x0,y0,e0,alpha,theta
*            else
*               read (18,50) pnum,x0,y0,e0,alpha,theta
*            endif
*c           could use non-dipole field model to compute L here
*            l0=sqrt(x0*x0+y0*y0)
*            e0=e0*1000.
*            a0=cos(alpha/raddeg)
*            call
*     &      pnt2fld3d0(pnum,l0,e0,a0,nl,ne,na,ls,ke,ca,1.0,jc_tp,j_tp)
*         enddo
*
*50       format (i8,5e12.4)
*         print *,'closing ',filename(1:lnblnk(filename))
*         close(18)
*      enddo
**      i=6
**      j=11
**      print *
**      print *,'un-weighted T.P. flux'
**      print *,'i,j,ls(i),ke(j)=',i,j,ls(i),ke(j)
**      do k=1,na       
**           print *,'k,ca,jc_tp,j_tp=',k,ca(k),
**     &      jc_tp(i,j,k),j_tp(i,j,k)
**      enddo
**      print *
*c     make it a differential directional flux
*      do k=1,na
*         call grd_delta(k,na,ca,upr,lwr,da)
*         da=abs(ca(k)*sin(acos(ca(k)))*(acos(lwr)-acos(upr)))
*         do j=1,ne
*            call grd_delta(j,ne,ke,upr,lwr,de)
*            do i=1,nl
*               call grd_delta(i,nl,ls,upr,lwr,dl)
*               dl=pi*(upr*upr-lwr*lwr)
*               dlea = abs(1.0/2./pi/dl/da/de/tdelta)
**               if ((i.eq.6).and.(j.eq.11)) then
**               print *,'k,dlea=',k,dlea
**               endif
*               j_tp(i,j,k)=j_tp(i,j,k)*dlea
*            enddo
*         enddo
*      enddo
*
*c     CALCULATE WEIGHT FUNCTION
*
*      do k=1,na
*         do j=1,ne
*            do i=1,nl
*               if (jc_tp(i,j,k).ge.mincnt*10.) then
*                  wght(i,j,k)=j_esa(i,j,k)/j_tp(i,j,k)
*               else 
*                  wght(i,j,k)=0.0
*               endif
*            enddo
*         enddo
*      enddo
**         i=6
**         j=11
**         print *
**         print *,'weighting function'
**         print *,'i,j,ls(i),ke(j)=',i,j,ls(i),ke(j)
**         do k=1,na       
**              print *,'j_esa,j_tp,wght=',
**     &         j_esa(i,j,k),j_tp(i,j,k),wght(i,j,k)
**         enddo

      return
      end

********************************************************************************

      subroutine info_fileread(basename,filenum)

c     obelt post processing
      implicit none
      include 'rbelt-pstprc.inc'

      integer i,filenum,lnblnk
      character*80 basename,fluxfile,filename,infofile

c     There are 51 parameters in info file
c     The following (globals) are declaired in pstprc.inc:
c     radius,expr,factor,dt_dist,emin,emax,binio,flux_dt,num_particles,
c     initdist,num_wsteps,wlines_dist,wlines_flux,wlines_prcp,rmin,rmax
c     year0,doy0,hour0,min0,sec0
c     The rest are not curently used in pstprc:
      real charge_rb,charge_sign_rb,m0_rb,b0_rb,dthalt
      real tstep_lrntz,dx_max,tstep_gc,dx_max_gc
      real go2lrntz,go2gc,dtgo2gc,init_t,epa_min,epa_max,lmin,lmax
      integer nt,nstep,num_particles_rb,max_yout,num_flts,num_ints
      integer max_wsteps,firstfilenum,seed,dist_seed,flag0
      integer num_wsteps_rb,nx_rb,ny_rb,nz_rb,lcalc
      character*80 basename_rb
      logical flag_switch

      filename=infofile(basename,filenum)
*      print *
!      print *,'opening ',filename(1:lnblnk(filename))
c     open info file & read
      open (12,file=filename(1:lnblnk(filename)))
c     from rbelt-grid.inc
      read (12,10) nx_rb
      read (12,10) ny_rb
      read (12,10) nz_rb
      read (12,10) nt
      read (12,10) nstep
c     from rbelt-const.inc
      read (12,20) charge_rb
      read (12,20) charge_sign_rb
      read (12,20) m0_rb
      read (12,20) b0_rb
c     from rbelt-y0.inc
      read (12,10) num_particles_rb
c     from rbelt-io.inc
      read (12,10) max_yout
      read (12,10) num_flts
      read (12,10) num_ints
      read (12,10) max_wsteps
c     set in rbelt-input.txt
c     rbelt namelist
      read (12,30) basename_rb
      read (12,10) firstfilenum
c     bounds namelist
      read (12,20) rmin
      read (12,20) rmax
      read (12,20) tmax
c     io namelist
      read (12,20) dthalt
      read (12,20) flux_dt
      read (12,20) init_twrite
      read (12,20) dtwrite
      read (12,40) binio
      read (12,40) flux_out
      read (12,40) dist_out
      read (12,40) prcp_out
      read (12,40) init_out
      read (12,40) cone_out
      read (12,10) lcalc
c     Lorentz namelist
      read (12,20) tstep_lrntz
      read (12,20) dx_max
c     guiding center namelist
      read (12,20) tstep_gc
      read (12,20) dx_max_gc
      read (12,20) go2lrntz
      read (12,20) go2gc
      read (12,20) dtgo2gc
      read (12,10) seed
c     dist namelist
      read (12,10) dist_seed
      read (12,20) dt_dist
      read (12,20) init_t
      read (12,20) emin
      read (12,20) emax
      read (12,20) epa_min
      read (12,20) epa_max
      read (12,20) lmin
      read (12,20) lmax
      read (12,20) radius
      read (12,20) expr
      read (12,20) factor
      read (12,10) flag0
      read (12,40) flag_switch
      read (12,10) initdist
c     declaired in rbelt-ut.inc
      read (12,10) year0
      read (12,10) doy0
      read (12,10) hour0
      read (12,10) min0
      read (12,10) sec0
c     declaired in rbelt-io.inc
      read (12,10) num_wsteps_rb
      read (12,10) wlines_dist
      read (12,10) wlines_flux
      read (12,10) wlines_prcp
      read (12,10) wlines_cone
!      print *,'closing ',filename(1:lnblnk(filename))
      close(12)

10    format (16x,i10)
20    format (16x,e12.4)
30    format (16x,a80)
40    format (16x,l1)

*      if (filenum.eq.firstfilenum) then
!      print*,'nx=',nx_rb
!      print*,'ny=',ny_rb
!      print*,'nz=',nz_rb
!      print*,'nt=',nt
!      print*,'nstep=',nstep
!      print*,'charge=',charge
!      print*,'charge_sign=',charge_sign
!      print*,'m0=',m0
!      print*,'b0=',b0
!      print*,'num_particles=',num_particles_rb
!      print*,'max_yout=',max_yout
!      print*,'num_flts=',num_flts
!      print*,'num_ints=',num_ints
!      print*,'max_wsteps=',max_wsteps
!      print*,'basename=',basename
!      print*,'firstfilenum=',firstfilenum
!      print*,'rmin=',rmin
!      print*,'rmax=',rmax
!      print*,'tmax=',tmax
!      print*,'dthalt=',dthalt
!      print*,'flux_dt=',flux_dt
!      print*,'init_twrite=',init_twrite
!      print*,'dtwrite=',dtwrite
!      print*,'binio=',binio
!      print*,'tstep_lrntz=',tstep_lrntz
!      print*,'dx_max=',dx_max
!      print*,'tstep_gc=',tstep_gc
!      print*,'dx_max_gc=',dx_max_gc
!      print*,'go2lrntz=',go2lrntz
!      print*,'go2gc=',go2gc
!      print*,'dtgo2gc=',dtgo2gc
!      print*,'seed=',seed
!      print*,'dist_seed=',dist_seed
!      print*,'dt_dist=',dt_dist
!      print*,'init_t=',init_t
!      print*,'emin=',emin
!      print*,'emax=',emax
!      print*,'epa_min=',epa_min
!      print*,'epa_max=',epa_max
!      print*,'lmin=',lmin
!      print*,'lmax=',lmax
!      print*,'radius=',radius
!      print*,'expr=',expr
!      print*,'factor=',factor
!      print*,'flag0=',flag0
!      print*,'flag_switch=',flag_switch
!      print*,'initdist=',initdist
!      print*,'year0=',year0
!      print*,'doy0=',doy0
!      print*,'hour0=',hour0
!      print*,'min0=',min0
!      print*,'sec0=',sec0
!      print*,'num_wsteps=',num_wsteps_rb
!      print*,'wlines_dist=',wlines_dist
!      print*,'wlines_flux=',wlines_flux
!      print*,'wlines_prcp=',wlines_prcp
!      print*,'wlines_prcp=',wlines_cone
!*      endif

      if (num_particles.ne.num_particles_rb) then
         print *
         print *,'num_particles in pstprc.inc needs to'
         print *,'be same as rbelt code num_particles'
         print *,'rbelt code num_particles=',num_particles_rb
         print *,'pstprc.inc num_particles=',num_particles
         print *,'change parameter num_particles in pstprc.inc to:',
     &   num_particles_rb
         stop
      endif

      if ((num_wsteps.ne.1).and.(num_wsteps_rb.eq.0)) then
         print *
         print *,'num_wsteps in pstprc.inc needs to'
         print *,'be = 1'
         print *,'rbelt code num_wsteps=',num_wsteps_rb
         print *,'pstprc.inc num_wsteps=',num_wsteps
         print *,'change parameter num_wsteps in pstprc.inc to 1'
         stop
      elseif (num_wsteps.ne.num_wsteps_rb.and.(num_wsteps_rb.ne.0)) then
         print *
         print *,'num_wsteps in pstprc.inc needs to'
         print *,'be same as rbelt code num_wsteps'
         print *,'rbelt code num_wsteps=',num_wsteps_rb
         print *,'pstprc.inc num_wsteps=',num_wsteps
         print *,'change parameter num_wsteps in pstprc.inc to:',
     &   num_wsteps_rb
         stop
      endif

      return
      end

********************************************************************************

      subroutine get_wtimes() 

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,tstep,count
      real sum

      print *
      print *,'*** in subroutine get_wtimes ***'

c     set up rbelt wtime array
!      do i=1,num_wsteps
!         wtime(i)=init_twrite+(i-1)*dtwrite
!         print *,'i,wtime',i,wtime(i)
!      enddo

c     set up pstprc time array
      if (prcp.eqv..true.) then
         do i=1,num_tsteps
            time(i)=init_time+(i-1)*dtstep
	    print *,'i,time',i,time(i)
         enddo
      else
         tstep=0
         sum=0
         count = 0
         do i=1,num_wsteps
            wtime(i)=init_twrite+(i-1)*dtwrite
c           get time corresonding to tstep
            sum=sum+wtime(i)
            count=count+1
            if (count.eq.stepavg) then
               tstep=tstep+1
               time(tstep)=sum/float(stepavg)
               sum=0
               count = 0
            endif
         enddo
      endif

      return
      end

********************************************************************************

      subroutine grid_setup()

c     set up post processing grids

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l
      real dx,dy,dz,dtht,dphi,dl,de,da,db

      print *
      print *,'*** in subroutine grid_setup ***'

      call mk1dgrd(nl,lsmin,lsmax,dl,ls)
      call mk1dgrd(ne,kemin,kemax,de,ke)
      call mk1dgrd(na,camin,camax,da,ca)
      call mk1dgrd(nb,bmmin,bmmax,db,bm)
      print *
      print *,'L,energy,cos(PA) grid'
      print *
      do i=1,nl
         print *,'i,ls(i)=',i,ls(i)
      enddo
      print *
      do i=1,ne
         print *,'i,ke(i)=',i,ke(i)
      enddo
      print *
      do i=1,na
         print *,'i,ca(i),alpha=',i,ca(i),acos(ca(i))*raddeg
      enddo
      do i=1,nb
         print *,'i,bm(i)',i,bm(i)
      enddo

c     set up distribution grids
      call mk1dgrd(nx,xmin,xmax,dx,x)
      call mk1dgrd(ny,ymin,ymax,dy,y)
*      print *
*      print *,'x grid'
*      do i=1,nx
*         print *,'i,x(i)=',i,x(i)
*      enddo
*      print *,'y grid'
*      do i=1,ny
*         print *,'i,y(i)=',i,y(i)
*      enddo
*      print *,'z grid'
*      do i=1,nz
*         print *,'i,z(i)=',i,z(i)
*      enddo

c     set up geographic coordinates grids for prcp file output
      call mk1dgrd(ntht,thtmin,thtmax,dtht,tht)
      call mk1dgrd(nphi,phimin,phimax,dphi,phi)
      print *
      print *,'tht grid'
      do i=1,ntht
         print *,'i,tht(i)=',i,tht(i)*raddeg
      enddo
      print *,'phi grid'
      do i=1,nphi
         print *,'i,phi(i)=',i,phi(i)*raddeg
      enddo

      print *
      print *,'********************************'
      print *,'********************************'

      return
      end

********************************************************************************

      subroutine wfgrid_setup()

c     set up post processing grids

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l
      real d1,d2,d3,d4

      print *
      print *,'*** in subroutine wfgrid_setup ***'

c     set up particle weight function grid
      call mk1dgrd(n1,min1,max1,d1,a1)
      call mk1dgrd(n2,min2,max2,d2,a2)
      call mk1dgrd(n3,min3,max3,d3,a3)
      call mk1dgrd(n4,min4,max4,d4,a4)
      print *
      do i=1,n1
         print *,'i,a1(i)=',i,a1(i)
      enddo
      print *
      do i=1,n2
         print *,'i,a2(i)=',i,a2(i)
      enddo
      print *
      do i=1,n3
         print *,'i,a3(i)=',i,a3(i)
      enddo
      print *
      do i=1,n4
         print *,'i,a4(i)=',i,a4(i)
      enddo

      print *
      print *,'********************************'
      print *,'********************************'

      return
      end

********************************************************************************

      subroutine init_arrays()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l

c     initialize flux distribution arrays to zero
      do k=1,num_tsteps
         do j=1,nphi
            do i=1,ntht
               j_tpt(i,j,k)=0.
               n_tpt(i,j,k)=0.
            enddo
         enddo
      enddo

      do k=1,num_tsteps
         do j=1,ny
            do i=1,nx
               j_xyt(i,j,k)=0.
               n_xyt(i,j,k)=0.
            enddo
         enddo
      enddo

      do l=1,num_tsteps
         do k=1,na
            do j=1,ne
               do i=1,nl
                  j_leat(i,j,k,l)=0.
                  n_leat(i,j,k,l)=0.
               enddo
            enddo
         enddo
      enddo

!        do j=1,num_files
!           do i=1,num_particles
!              pwghtflag

      return
      end


********************************************************************************

      subroutine init_fileread(basename,filenum)

c     reads particle initial conditions and weights particles

      implicit none
      include 'rbelt-pstprc.inc'
      include 'rbelt-geopak.inc'

      character*80 basename,initfile,filename
      integer n,lnblnk,filenum,i,status,grdpos0
      real x0,y0,z0,l0,e0,a0,ca0,t0,pw,fld2pt3d0

!      print *
!      print *,'*** in subroutine init_fileread ***'

c     name & open input file
      filename=initfile(basename,filenum)
!      print *
!      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (18,file=filename(1:lnblnk(filename)),
     &   form='unformatted')
      else
         open (18,file=filename(1:lnblnk(filename)))
      endif

c     read each line in file
      do n=1,num_particles

c        read prcp count
         if (binio.eqv..true.) then
            read (18) x0,y0,z0,l0,e0,a0,t0
         else
            read (18,10) x0,y0,z0,l0,e0,a0,t0
         endif
10       format (7e12.4)

c        get cos of alpha
         ca0=cos(a0/raddeg)

c        weight particles -- each particle gets a weight. 
c        The weight is a function of its position in L-E-COS(PA) space.
         pwght(n,fcount)=fld2pt3d0(n,l0,e0,ca0,n1,n2,n3,a1,a2,a3,wf)
!         pgrid1(n,1) = l0
!         pgrid1(n,2) = e0
!         pgrid1(n,3) = a0
      
         if (pwght(n,fcount) .gt. 100000000000) then
            print *,'pwght(n,fcount) .gt. 100000000)'
            print *,'pwght(n,fcount)=',pwght(n,fcount)
            print *,'n,l0,e0,ca0=',n,l0,e0,ca0
            stop
         endif

*         print *
*         print *,'l0,e0,ca0',l0,e0,ca0
*         print *,'n,fcount,pwght=',n,fcount,pwght(n,fcount)

      enddo

!      print *,'closing ',filename(1:lnblnk(filename))
      close(18)
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************

      subroutine prcp_fileread(basename,filenum)

      implicit none
      include 'rbelt-pstprc.inc'
      include 'rbelt-geopak.inc'

      character*80 basename,prcpfile,filename
      integer lnblnk,k,pnum,count(num_tsteps),filenum,tstep
      integer year,doy,hour,min,sec
      real t0,x0,y0,z0,e0,a0,ceta,ie,ia,it,tht0,phi0,pw,fld2pt3d0
      real XGSM,YGSM,ZGSM

      print *
      print *,'*** in subroutine prcp_fileread ***'

c     name & open input file
      filename=prcpfile(basename,filenum)
      print *
      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (18,file=filename(1:lnblnk(filename)),
     &   form='unformatted')
      else
         open (18,file=filename(1:lnblnk(filename)))
      endif
c     read each line in file
      pnum=1
      do k=1,num_tsteps
         count(k)=0
      enddo
      print *,'wlines_prcp=',wlines_prcp
      do k=1,wlines_prcp
c        read prcp count
         if (binio.eqv..true.) then
            if (init_out.eqv..true.) then
               read (18) pnum,t0,x0,y0,z0,e0,a0,ceta
            else
               read (18) pnum,t0,x0,y0,z0,e0,a0,ceta,ie,ia,it
            endif
         else
            if (init_out.eqv..true.) then
               read (18,10) pnum,t0,x0,y0,z0,e0,a0,ceta
            else
*               read (18,20) pnum,t0,x0,y0,z0,e0,a0,ceta,ie,ia,it
               read (18,20) t0,x0,y0,z0,e0,a0,ceta,ie,ia,it
            endif
10          format (i8,7e12.4)
*20          format (i8,10e12.4)
20          format (10e12.4)
         endif

c        put SM coordinates x0,y0,z0 into GEO lat. & lon.
c        get date and time of count
         year=year0
         doy=doy0
         hour=hour0
         min=min0
         sec=sec0+t0
         call reduce_time2(year,doy,hour,min,sec)
c        call recalc and rotate to geographic coordinates
*         print *
*         print *,'calling RECALC w/ YEAR,DOY,HOUR,MIN,SEC=',
*     &   year,doy,hour,min,sec
         call RECALC (year,doy,hour,min,sec)
*         print *,'done with RECALC: tilt=',atan2(SPS,CPS)*raddeg
*         print *,'SM: x0,y0,z0=',x0,y0,z0
         call SMGSM (x0,y0,z0,XGSM,YGSM,ZGSM,1)
         call GEOGSM (x0,y0,z0,XGSM,YGSM,ZGSM,-1)
*            print *,'GEO: x0,y0,z0=',x0,y0,z0
c        get latitude (-pi/2 to +pi/2) and longitude (-pi to +pi) in geographic coordinates
         tht0=acos(z0/sqrt(x0*x0+y0*y0+z0*z0))
         phi0=atan2(y0,x0)

c        get cos of initial pitch angle
         ia=cos(ia/raddeg)
c        get particle weight
         if (init_out.eqv..true.) then
            pw=pwght(pnum,fcount)
         else
c           weighting in terms of energy, PA, time
c           need to divide by ceta here to get a differential
c           in direction flux (use NGP on ca grid to avoid NANs).
            pw=fld2pt3d0(pnum,ie,ia,it,n1,n2,n3,a1,a2,a3,wf)
*             pw=1.0
         endif

c        get tstep corresponding to t0
c        this is set up for tsteps = num_wsteps and needs to be fixed -- 
c        see prcp_wghtcnts2flux()
c        we can assume uniform time steps
c        for now we use 0th order interp in time

         if ((t0.ge.time(1)).and.(t0.lt.time(num_tsteps))) then
            tstep=int((t0-init_time+dtstep/2.)/dtstep)+1
            count(tstep)=count(tstep)+1
         else
            tstep=0
         endif

         call pt2fld2d0t(pnum,tht0,phi0,tstep,ntht,nphi,
     &   num_tsteps,tht,phi,pw,n_tpt,j_tpt)
c
*         call pt2fld3d0(pnum,tht0,phi0,t0,ntht,nphi,
*     &   num_wsteps,tht,phi,wtime,pw,n_tpt,j_tpt)

      enddo

      print *,'closing ',filename(1:lnblnk(filename))
      close(18)
      return
      end

********************************************************************************

      subroutine flux_fileread(basename,filenum,mode)

      implicit none
      include 'rbelt-pstprc.inc'
      include 'rbelt-geopak.inc'

      character*80 basename,fluxfile,filename
      integer lnblnk,i,k,pnum,count(num_tsteps),filenum,wstep,tstep,
     &mode,status,grdpos0
      real x0,y0,z0,lshell,energy,alpha,calpha,eta,ceta,fld2pt3d0,
     &e0,a0,ca0,th0,cth0,t,t0,pw,pwghtfn
      real bmir
      integer pstop, pstop_num, ipsave
      real psave

!      print *
!      print *,'*** in subroutine flux_fileread ***'
!      if (mode.eq.0) then
!         print *,'getting test-particle flux'
!      elseif (mode.eq.1) then
!         print *,'getting weighted flux'
!      endif

c     name & open input file
      filename=fluxfile(basename,filenum)
!      print *
!      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (18,file=filename(1:lnblnk(filename)),
     &   form='unformatted')
      else
         open (18,file=filename(1:lnblnk(filename)))
      endif

c     read each line in file
*      do k=1,num_tsteps
*         count(k)=0
*      enddo
!      print *,'wlines_flux=',wlines_flux
      pstop_num = 0
      psave = 0
      ipsave = 0
      pstop = 1
      do k=1,wlines_flux
c        read prcp count
        
         if (binio.eqv..true.) then
            if (method.eq.1) then
               read (18) pnum,wstep,x0,y0,lshell,energy,alpha,eta
            elseif (method.eq.2) then
               read (18) pnum,wstep,x0,y0,lshell,energy,alpha,eta,
     &         e0,a0,t0
            endif
         else
            if (method.eq.1) then
               read (18,10) pnum,wstep,x0,y0,lshell,energy,alpha,eta
            elseif (method.eq.2) then
               read (18,20)pnum,wstep,x0,y0,lshell,energy,alpha,eta,
     &         e0,a0,t0
            endif
         endif
10       format (2i8,6d16.8)
20       format (2i8,9e12.4)

!         if ((pstop_num.ne.pnum).and.(pstop.eq.1).and.(mode.eq.2)) then
!           pstop = 0
!           psave = 0
!           ipsave = 0
!         endif 

!         print *, k, pnum, wstep, x0, y0
!        wstep = wstep - 8280
c        time of write step
         t=wtime(wstep)

c        convert MeV to KeV (if needed)
*         energy=energy*1000
*         e0=e0*1000.

c        get cos of alpha and eta
         calpha=cos(alpha/raddeg)
         ceta=cos(eta/raddeg)
         ca0=cos(a0/raddeg)
         bmir = eta/(sin(alpha/raddeg)**2)
c        get particle weight
         if (mode.ne.2) then
c           get particle weight / cos(eta)
            status = 0
c           *** note, here we use cos(alpha_o) grid to approximate
c           cos of angle of incidence with equatorial plane
!            i=grdpos0(calpha,na,ca,status)
            if (status.eq.0) then

c              mode=0 to get test-particle flux
               if (mode.eq.0) then
                  pw=1.0/abs(ca(i))

c              mode=1 for weighted particles (to get weighted flux)
               elseif (mode.eq.1) then
                  if (method.eq.1) then
c                    method=1 means they are pre-weighted
                     pw=pwght(pnum,fcount)/abs(ca(i))
                  else
c                    compute weight here
c                    weighting in terms of energy, PA, time *                     
*                     pw=fld2pt3d0(pnum,e0,ca0,t,n1,n2,n3,a1,a2,a3,wf)/
*     &               abs(ca(i))
                     pw=pwghtfn(pnum,e0,ca0,t)/abs(ca(i))
                  endif

               endif

c              get time step & time
               tstep = int((float(wstep)-.5)/float(stepavg))+1
               t=time(tstep)
!               if ((tstep.eq.itstep).and.(mode.eq.1)) then
!                 pgrid2(pnum,1) = lshell
!                 pgrid2(pnum,2) = energy
!                 pgrid2(pnum,3) = alpha
!                 print *, pw, pwght(pnum,fcount)
!                 print *, '***********************************'
!               endif
c              add weighted count to L-E-PA-time distribution
!               if (mode.eq.1) then
!                print *, lshell, energy, calpha, pw
!               endif
               call pt2fld3d0t(pnum,lshell,energy,calpha,tstep,nl,ne,na,
     &         num_wsteps,ls,ke,ca,pw,n_leat,j_leat)

c              add weighted count to omni-integral x-y-time distribution
!               call pt2fld2d0t
!     &         (pnum,x0,y0,tstep,nx,ny,num_tsteps,x,y,pw,n_xyt,j_xyt)

            endif
         else
c           weight particles
c           weighting here in terms of L,E,EPA
            tstep = int((float(wstep)-.5)/float(stepavg))+1
            if (tstep.le.itstep) then
!              if (tstep.eq.itstep) then
!                ipsave = ipsave + 1
!                psave = fld2pt3d0
!     &          (pnum,lshell,energy,calpha,n1,n2,n3,a1,a2,a3,wf)+psave
!              elseif (tstep.lt.itstep) then
                pwght(pnum,fcount)=fld2pt3d0
     &          (pnum,lshell,energy,calpha,n1,n2,n3,a1,a2,a3,wf)
!              endif
!            elseif ((tstep.gt.itstep).and.(pstop.ne.1)) then
!              pstop = 1
!              pstop_num = pnum
!              if (ipsave.gt.0) pwght(pnum,fcount) = psave/ipsave
            endif
         endif
      
      enddo

!      print *,'closing ',filename(1:lnblnk(filename))
      close(18)
      return
      end

*c     read each line in file
*      do k=1,wlines_flux
*c           read flux count
*         if (binio.eqv..true.) then
*            read (18) pnum,wstep,tag,tmp1,tmp2,e0,alpha,eta
*         else
*            read (18,10) pnum,wstep,tag,tmp1,tmp2,e0,alpha,eta
*         endif
*8        format (2i8,5e12.4)
*10       format (3i8,5e12.4)
*         if (tag.eq.1) then
*            xy_tot(wstep)=xy_tot(wstep)+1
*            x0=tmp1
*            y0=tmp2
*            call pt2fld2d1t(pnum,x0,y0,wstep,nx,ny,
*     &      num_wsteps,x,y,1.0,j_xyt,n_xyt)
*         endif
*         if (tag.eq.2) then
*            xz_tot(wstep)=xz_tot(wstep)+1
*            x0=tmp1
*            z0=tmp2
*             call pt2fld2d1t(pnum,x0,z0,wstep,nx,nz,
*     &      num_wsteps,x,z,1.0,j_xzt,n_xzt)
*         endif
*         if (tag.eq.3) then
*            yz_tot(wstep)=yz_tot(wstep)+1
*            y0=tmp1
*            z0=tmp2
*             call pt2fld2d1t(pnum,y0,z0,wstep,ny,nz,
*     &      num_wsteps,y,z,1.0,j_yzt,n_yzt)
*
*         endif
*      enddo
*      print *,'closing ',filename(1:lnblnk(filename))
*      close(18)

********************************************************************************

      subroutine prcp_wghtcnts2flux()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,pnum,count
      real dtpt,upr,lwr
      real dtht,dphi,dt

      print *
      print *,'*** in subroutine prcp_wcnts2flux ***'

c     expected flux (from total particles launched into sphere over dt_dist time interval):
*      print *
*      print *,'J_omni/2 = #/RE^2-sec/2 launched into sphere='
*      print *,num_particles*num_files/pi/radius/radius/dt_dist/2.

c      calculate #/RE^2-sec through surface (*NOTE*, this is J_omni/2)
      count=0
      do k=1,num_tsteps
         call grd_delta(k,num_tsteps,time,upr,lwr,dt)
         do j=1,nphi
            call grd_delta(j,nphi,phi,upr,lwr,dphi)
            do i=1,ntht
               call grd_delta(i,ntht,tht,upr,lwr,dtht)
               if (n_tpt(i,j,k).gt.0.0) count=count+1
               if (n_tpt(i,j,k).ge.mincnt) then
                  dtpt=1.0/(rmin*rmin*(cos(lwr)-cos(upr))*dphi)/dt
                  j_tpt(i,j,k)=j_tpt(i,j,k)*dtpt

*                  print *,'i,j,k,j=',i,j,k,j_tpt(i,j,k)

               else
                  j_tpt(i,j,k)=0.0
               endif
            enddo
         enddo
      enddo

      print *,'count=',count

      return
      end

********************************************************************************

      subroutine wghtcnts2jxyt()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k
      real dxyt,upr,lwr,dx,dy

      print *
      print *,'*** in subroutine wghtcnts2jxyt ***'

c     exprected flux (from total particles launched into sphere over dt_dist time interval):
*      print *
*      print *,'J_omni/2 = #/RE^2-sec/2 launched into sphere='
*      print *,num_particles*num_files/pi/radius/radius/dt_dist/2.

c      calculate #/RE^2-sec through surface (*NOTE*, this is J_omni/2)
      do k=1,num_wsteps
         do j=1,ny
            call grd_delta(j,ny,y,upr,lwr,dy)
            do i=1,nx
               call grd_delta(i,nx,x,upr,lwr,dx)
*               if (n_xyt(i,j,k).gt.0.0) count=count+1

*               if (n_xyt(i,j,k).lt.0.0) then
*                  print *,'i,j,k,n_xyt=',i,j,k,n_xyt(i,j,k)
*                  stop
*               endif

               if (n_xyt(i,j,k).ge.mincnt) then
                  dxyt = abs(1.0/dx/dy/flux_dt/stepavg)
*                  print *,'j_xyt(i,j,k)=',j_xyt(i,j,k)
                  j_xyt(i,j,k)=j_xyt(i,j,k)*dxyt
               else
                  j_xyt(i,j,k)=0.0
               endif
            enddo
         enddo
      enddo
*      print *,'count=',count

      return
      end

********************************************************************************

      subroutine wghtcnts2jleat()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l
      real dleat,upr,lwr,dl,de,da,fcntmin,fcntmax,nzeros,db

      print *
      print *,'*** in subroutine wghtcnts2jleat ***'

      fcntmin=10000000000.
      fcntmax=0.

      do l=1,num_tsteps
         do k=1,na
            call grd_delta(k,na,ca,upr,lwr,da)
            da=abs(ca(k)*sin(acos(ca(k)))*(acos(lwr)-acos(upr)))
*            da=abs(sin(acos(ca(k)))*(acos(lwr)-acos(upr)))
            do j=1,ne
               call grd_delta(j,ne,ke,upr,lwr,de)
               do i=1,nl
                  call grd_delta(i,nl,ls,upr,lwr,dl)
c                 dipole L-shell
                  dl=pi*(upr*upr-lwr*lwr)
                  dleat = abs(1.0/2./pi/dl/db/de/flux_dt/stepavg)
                  dleat = 1.
                  if (n_leat(i,j,k,l).ge.mincnt) then
                     j_leat(i,j,k,l)=j_leat(i,j,k,l)*dleat
                  else
                     j_leat(i,j,k,l)=0.0
                  endif
               
                  if (n_leat(i,j,k,l).lt.fcntmin)fcntmin=n_leat(i,j,k,l)
                  if (n_leat(i,j,k,l).gt.fcntmax)fcntmax=n_leat(i,j,k,l)
                  if (n_leat(i,j,k,l).eq.0.) nzeros=nzeros+1

               enddo
            enddo
         enddo
      enddo

      print *
      print *,'fcntmin,fcntmax,nzeros=',fcntmin,fcntmax,nzeros

      return
      end

********************************************************************************

      subroutine write_prcp()

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k

      print *
      print *,'*** in subroutine write_prcp ***'

c     open j_tpt.dat file & write flux
      print *,'opening writing to & closing ','j_tpt.dat'
      open (20,file='j_tpt.dat')

      write (20,10) year0,doy0,hour0,min0,sec0
      write (20,20) ntht,nphi,num_tsteps
      write (20,30) thtmin*raddeg,thtmax*raddeg,
     & phimin*raddeg,phimax*raddeg
c     eventually we should remove k,tht & phi as they are not needed
      do k=1,num_tsteps
         write (20,40) time(k)
         do j=1,nphi
            do i=1,ntht
               write (20,40) tht(i)*raddeg,phi(j)*raddeg,j_tpt(i,j,k)
            enddo
         enddo
      enddo
10    format (5i8)
20    format (3i8)
30    format (4e12.4)
*40    format (i8,e12.4)
40    format (e12.4)
50    format (3e12.4)
      close(20)

      return
      end

********************************************************************************

      subroutine write_jxy()
      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k
      print *
      print *,'*** in subroutine write_jxy ***'
c     open j_xyt.dat file & write flux
      print *,'opening writing to & closing ','j_xyt.dat'
      open (20,file='j_xyt.dat')
      write (20,10) year0,doy0,hour0,min0,sec0
      write (20,20) nx,ny,num_tsteps
      write (20,30) xmin,xmax,ymin,ymax
      do k=1,num_tsteps
*         print *,'k,time(k)=',k,time(k)
         write (20,40) wtime(k)
         do j=1,ny
            do i=1,nx
               write (20,40) j_xyt(i,j,k)
*               print *,'n_xyt(i,j,k)=',n_xyt(i,j,k)
            enddo
         enddo
      enddo
10    format (5i8)
20    format (3i8)
30    format (4e12.4)
40    format (e12.4)
      close(20)
      return
      end


********************************************************************************

      subroutine write_jlea()
      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k,l
      print *
      print *,'*** in subroutine write_jlea ***'
c     open j_xyt.dat file & write flux
      print *,'opening writing to & closing ','j_lea.dat'
      open (20,file='j_lea.dat')
      write (20,10) year0,doy0,hour0,min0,sec0
      write (20,20) nl,ne,na,num_tsteps
      write (20,30) lsmin,lsmax,kemin,kemax,camin,camax

      do l=1,num_tsteps
         write (20,40) wtime(l)
         do k=1,na
            do j=1,ne
               do i=1,nl
                  write (20,40) j_leat(i,j,k,l)
               enddo
            enddo
         enddo
      enddo

10    format (5i8)
20    format (4i8)
30    format (6e12.4)
40    format (e12.4)
      close(20)
      return
      end

********************************************************************************

      subroutine bcs_lin3d(ni,nj,nk,f)
      implicit none
      integer i,j,k,ni,nj,nk,nip1,njp1,nkp1,nip2,njp2,nkp2
      real f(ni,nj,nk)

      nip1=ni+1
      njp1=nj+1
      nkp1=nk+1
      nip2=ni+2
      njp2=nj+2
      nkp2=nk+2
      
      do i=2,nip1
         do j=2,njp1
            f(i,j,1)=2.*f(i,j,2)-f(i,j,3)
            f(i,j,nkp2)=2.*f(i,j,nkp1)-f(i,j,nk)
         enddo
      enddo
      do j=2,njp1
         do k=2,nkp1
            f(1,j,k)=2.*f(2,j,k)-f(3,j,k)
            f(nip2,j,k)=2.*f(nip1,j,k)-f(ni,j,k)
         enddo
      enddo
      do k=2,nkp1
         do i=2,nip1
            f(i,1,k)=2.*f(i,2,k)-f(i,3,k)
            f(i,njp2,k)=2.*f(i,njp1,k)-f(i,nj,k)
         enddo
      enddo

      return
      end

********************************************************************************

      subroutine esa_print(ls,ke,ea,jdist)

      implicit none
      include 'rbelt-const.inc'
      integer nl,ne,na,i,j,k,sum,sum2
      parameter (nl=19,ne=13,na=22)
      real ls(nl),ke(ne),ea(na),jdist(nl,ne,na)

c     print out ESA dist
      print *
      print *,'total epa j dist.'
      sum2=0.
      do k=1,na
         sum = 0.
         do j=1,ne
            do i=1,nl
               sum=sum+abs(jdist(i,j,k))
            enddo
         enddo
!         print *,'k,ea(k),j=',k,float(acos(ea(k))*raddeg),
!     &   float(sum)
         sum2=sum2+sum
      enddo
      print *,'total j counts=',float(sum2)
      print *
      print *,'total energy j dist.'
      sum2=0.
      do j=1,ne
         sum = 0.
         do k=1,na
            do i=1,nl
               sum=sum+abs(jdist(i,j,k))
            enddo
         enddo
!1         print *,'j,ke(j),j=',j,float(ke(j)),float(sum)
         sum2=sum2+sum
      enddo
      print *,'total j counts=',float(sum2)
      print *
      print *,'total L-shell j dist.'
      sum2=0.
      do i=1,nl
         sum = 0.
         do j=1,ne
            do k=1,na
               sum=sum+abs(jdist(i,j,k))
            enddo
         enddo
!         print *,'i,ls(i),j=',i,float(ls(i)),float(sum)
         sum2=sum2+sum
      enddo
      print *,'total j counts=',float(sum2)
      print *

      return
      end

********************************************************************************

      subroutine cos_pa_bin(basename,firstfilenum)

c     zeroth order binning of flux counts, test particle fluxes, 
c     and weighted fluxes; in pitch angle, x-y grid, and energy separatly

c     obelt post processing
      implicit none
      include 'rbelt-const.inc'
      include 'rbelt-io.inc'
      include 'rbelt-dist.inc'
*      include 'rbelt-grid.inc'
      integer firstfilenum,lnblnk,i,j,k,grdpos,idummy,err,nfcnts,
     &num_particles
      character*80 basename,fluxfile,filename
      real rdummy,time,fld2pt3d
      real x0,y0,e0,l0,a0,alpha,theta,upr,lwr,area,center
      real dl,sum,x,y,stheta

      integer ne,na,tcounts
      parameter(ne=99,na=22)
      real*8 camax,camin,ca(na)
      integer ke_bin(ne),ca_bin(na)
      real*8 de,da,hda,dflux(na)

      print *
      print *,'*** in subroutine flux_dist ***'

c     un normalize flux_dt
      flux_dt=flux_dt/tfactor

c     set up equitorial cos pitch angle bins
      camin=-1.0
      camax=1.0
      da=(camax-camin)/na
      hda=da/2.
      do k=1,na
         ca(k) = camin + (float(k)-0.5)*da
         print *,'k,lwr,ca(k),upr=',k,ca(k)-hda,ca(k),ca(k)+hda
      enddo
      do i = 1,na
         ca_bin(i) = 0
      enddo

c     open flux file
c     need to fix below
*      filename=fluxfile(basename,dist_seed)
      print *,'opening ',filename(1:lnblnk(filename))
      open (18,file=filename(1:lnblnk(filename)))
      read (18,10) num_particles
10    format (i8)
      print *,'num_particles=',num_particles
c     read particle initial conditions
      do i=1,num_particles
         read (18,20) x,y,e0,alpha,theta
      enddo
20    format (5e12.4)
      read (18,10) idummy
      print *,'idummy=',idummy

      do k=1,idummy

c        read in flux counts
         read (18,40) nfcnts,rdummy
40       format (i8,e12.4)
         print *,'nfcnts,rdummy=',nfcnts,rdummy
         do i=1,nfcnts
            read (18,50) idummy,x,y,e0,alpha,theta
c           could use non-dipole field model to compute L here
            l0=sqrt(x0*x0+y0*y0)
*            e0=e0*1000.
            a0=cos(alpha/raddeg)
            theta=cos(theta/raddeg)

c           fill pitch angle bins
            do j = 1,na
               call grd_delta(j,na,ca,upr,lwr,da)
               if ((a0.ge.(ca(j)-hda)).and.
     &         (a0.lt.(ca(j)+hda))) ca_bin(j) = ca_bin(j)+1
            enddo

         enddo
50       format (i8,5e12.4)

         print *
         print *,'flux ca dist.'
c        use sum to calc OMNI DIR. flux
         sum=0
         tcounts=0

c        use sum to calc OMNI DIR. flux

         do j = 1,na
            tcounts=tcounts+ca_bin(j)
            center=abs(ca(j)*sin(acos(ca(j)))*
     &      (acos(ca(j)-hda)-acos(ca(j)+hda)))
            dflux(j)=ca_bin(j)/(2.*pi*pi*radius*radius*center*flux_dt)
            print *,ca(j)-hda,ca(j)+hda,ca(j),dflux(j)
c           use sum to calc OMNI DIR. flux
            sum = sum + abs(2*pi*sin(acos(ca(j)))*
     &      (acos(ca(j)-hda)-acos(ca(j)+hda))*dflux(j))
         enddo

         print *,'tcounts=',tcounts
         print *,'omni dir flux=',sum
         print *,'omni dir flux/4pi = j =',sum/4./pi
         print *,'total through sphere is 2 * total through disk'
         print *,'j=tcounts/2./pi/pi/radius/radius/flux_dt = ',
     &   tcounts/2./pi/pi/radius/radius/flux_dt

      enddo

      close(18)
      return
      end

********************************************************************************
********************************************************************************
********************************************************************************
*
*      subroutine fluxmap_3planes(basename,firstfilenum,filenumstep)
*
*      implicit none
*      include 'rbelt-pstprc.inc'
*
*      character*80 basename,fluxfile,filename
*      integer firstfilenum,filenumstep,lnblnk,i,j,k,m,
*     &   pnum,filenum,wstep,tag
*      real dx,dy,dz,pw,upr,lwr,wtime(num_wsteps),tmp1,tmp2
*      real x0,y0,z0,e0,l0,a0,alpha,theta,dxyt,dxzt,dyzt,ie
*
**      real*8 ie 
*
*      real j_xyt(nx,ny,num_wsteps),n_xyt(nx,ny,num_wsteps)
*      real j_xzt(nx,nz,num_wsteps),n_xzt(nx,nz,num_wsteps)
*      real j_yzt(ny,nz,num_wsteps),n_yzt(ny,nz,num_wsteps)
*      integer xy_tot(num_wsteps),xz_tot(num_wsteps),yz_tot(num_wsteps)
*
*      print *
*      print *,'*** in subroutine flux_map_new ***'
*
**      if (initdist.eq.3) then
**         print *
**         print *,'j_tp=',(pper_file*num_files/
**     &    (4.0*pi*pi*radius*radius*(emax-emin)*dt_dist))
**         do i=1,ne
**            wght(i)=norm*(ke(i)/1000.)**spctrm/(pper_file*num_files/
**     &       (4.0*pi*pi*radius*radius*(emax-emin)*dt_dist))
**            print *,'ke,j_obs,wght=',
**     &       ke(i),norm*(ke(i)/1000.)**spctrm,wght(i)
**         enddo
**      endif
*
**c     determine write rbelt times (NO, we need'nt do this here)
**      k=-1
**5     k=k+1
**      if ((k+2).gt.max_wsteps) then 
**         print *,'num_wsteps+2 > max_wsteps'
**         print *,'make max_wsteps bigger'
**         stop
**      endif
**      wtime(k+1)=init_twrite+k*dtwrite
**      if (wtime(k+1).le.tmax) then
**         print *,'   wstep,wtime=',k+1,wtime(k+1)
**         goto 5
**      endif
**      print *,'num. of wsteps=',num_wsteps
**      print *
*
*c     initialize flux arrays to zero
*      do k=1,num_wsteps
*         xy_tot(k)=0
*         xz_tot(k)=0
*         yz_tot(k)=0
*         do j=1,ny
*            do i=1,nx
*               j_xyt(i,j,k)=0.
*               n_xyt(i,j,k)=0.
*            enddo
*         enddo
*         do j=1,nz
*            do i=1,nx
*               j_xzt(i,j,k)=0.
*               n_xzt(i,j,k)=0.
*            enddo
*         enddo
*         do j=1,nz
*            do i=1,ny
*               j_yzt(i,j,k)=0.
*               n_yzt(i,j,k)=0.
*            enddo
*         enddo
*      enddo
*
*c     loop over pstprc input files (rbelt output files)
**      do m=1, num_files
*         filenum=firstfilenum+(m-1)*filenumstep
*         call info_fileread(basename,filenum)
*         filename=fluxfile(basename,filenum)
**         print *
*         print *,'opening ',filename(1:lnblnk(filename))
*c        open input file
*         if (binio.eqv..true.) then
*            open (18,file=filename(1:lnblnk(filename)),
*     &      form='unformatted')
*         else
*            open (18,file=filename(1:lnblnk(filename)))
*         endif
*c        read each line in file
*         do k=1,wlines_flux
*c           read flux count
*            if (binio.eqv..true.) then
*               if (initdist.eq.3) then
*                  read (18) pnum,wstep,tmp1,tmp2,e0,alpha,theta
*               else
*                  read (18) pnum,wstep,tag,tmp1,tmp2,e0,alpha,theta
*               endif
*            else
*               if (initdist.eq.3) then
*                  read (18,8) pnum,wstep,tmp1,tmp2,e0,alpha,theta
*               else
*                  read (18,10) pnum,wstep,tag,tmp1,tmp2,e0,alpha,theta
*               endif
*8              format (2i8,5e12.4)
*10             format (3i8,5e12.4)
*            endif
*       
*c           could use non-dipole field model to compute L here
**            l0=sqrt(x0*x0+y0*y0)
**            ie=ie*1000.
**            e0=e0*1000.
**            a0=dcos(alpha/raddeg)
**            if (initdist.eq.3) then
**               pw=pwght(pnum,m)
**            else
**               pw=fld2pt1d0(pnum,ie,ne,ke,wght)
**            endif
*
*            if (tag.eq.1) then
*               xy_tot(wstep)=xy_tot(wstep)+1
*               x0=tmp1
*               y0=tmp2
*               call pt2fld2d1t(pnum,x0,y0,wstep,nx,ny,
*     &         num_wsteps,x,y,1.0,j_xyt,n_xyt)
*            endif
*            if (tag.eq.2) then
*               xz_tot(wstep)=xz_tot(wstep)+1
*               x0=tmp1
*               z0=tmp2
*                call pt2fld2d1t(pnum,x0,z0,wstep,nx,nz,
*     &         num_wsteps,x,z,1.0,j_xzt,n_xzt)
*            endif
*            if (tag.eq.3) then
*               yz_tot(wstep)=yz_tot(wstep)+1
*               y0=tmp1
*               z0=tmp2
*                call pt2fld2d1t(pnum,y0,z0,wstep,ny,nz,
*     &         num_wsteps,y,z,1.0,j_yzt,n_yzt)
*
*            endif
*         enddo
*         print *,'closing ',filename(1:lnblnk(filename))
*         close(18)
**      enddo
*
*c     exprected flux (from total particles launched into sphere over dt_dist time interval):
**      print *
**      print *,'J_omni/2 = #/RE^2-sec/2 launched into sphere='
**      print *,num_particles*num_files/pi/radius/radius/dt_dist/2.
*
*      do k=1,num_wsteps
*         print *
*         print *,'flux integration time (sec)=',
*     &   wtime(k)-flux_dt,'-',wtime(k)
*         print *,'wstep, total particles through xy,xz,yz planes ='
*         print *,k,xy_tot(k),xz_tot(k),yz_tot(k)
*         print *,'approx. particles per grid point (inside sphere only)'
*         print *,'in xy,xz,yz planes ='
*         print *,4*xy_tot(k)/nx/ny/pi,4*xz_tot(k)/nx/nz/pi,
*     &   4*yz_tot(k)/ny/nz/pi
*         print *,'if dist. is isotropic flux (#/RE^2-sec)'
*         print *,'through xy,xz,yz planes ='
*         print *,xy_tot(wstep)/pi/rmax/rmax/flux_dt,
*     &   xz_tot(wstep)/pi/rmax/rmax/flux_dt,
*     &   yz_tot(wstep)/pi/rmax/rmax/flux_dt
*      enddo
*      print *
*
*c      calculate #/RE^2-sec through x-y, x-z, & y-z planes (*NOTE*, this is J_omni/2)


*      do k=1,num_wsteps
*         do j=1,ny
*            call grd_delta(j,ny,y,upr,lwr,dy)
*            do i=1,nx
*               call grd_delta(i,nx,x,upr,lwr,dx)
*               if (n_xyt(i,j,k).ge.mincnt) then
*                  dxyt = abs(1.0/dx/dy/flux_dt)
**                  print *,'j_xyt(i,j,k)=',j_xyt(i,j,k)
*                  j_xyt(i,j,k)=j_xyt(i,j,k)*dxyt
*               else
*                  j_xyt(i,j,k)=0.0
*               endif
*            enddo
*         enddo


*         do j=1,nz
*            call grd_delta(j,nz,z,upr,lwr,dz)
*            do i=1,nx
*               call grd_delta(i,nx,x,upr,lwr,dx)
*               if (n_xzt(i,j,k).ge.mincnt) then
*                  dxzt = abs(1.0/dx/dz/flux_dt)
*                  j_xzt(i,j,k)=j_xzt(i,j,k)*dxzt
*               else
*                  j_xzt(i,j,k)=0.0
*               endif
*            enddo
*         enddo
*         do j=1,nz
*            call grd_delta(j,nz,z,upr,lwr,dz)
*            do i=1,ny
*               call grd_delta(i,ny,y,upr,lwr,dy)
*               if (n_yzt(i,j,k).ge.mincnt) then
*                  dyzt = abs(1.0/dy/dz/flux_dt)
*                  j_yzt(i,j,k)=j_yzt(i,j,k)*dyzt
*               else
*                  j_yzt(i,j,k)=0.0
*               endif
*            enddo
*         enddo
*      enddo
*
*c     open j_xyt.dat file & write flux
*      print *,'opening writing to & closing ','j_xyt.dat'
*      open (20,file='j_xyt.dat')
*      write (20,20) nx,ny,nz,num_wsteps
*      write (20,30) xmin,xmax,ymin,ymax,zmin,zmax
*      do k=1,num_wsteps
*         write (20,50) k,wtime(k)
*         do j=1,ny
*            do i=1,nx
*               write (20,40) j_xyt(i,j,k)
**               print *,'n_xyt(i,j,k)=',n_xyt(i,j,k)
*            enddo
*         enddo
*      enddo
*
*c     open j_xzt.dat file & write flux
*      print *,'opening writing to & closing ','j_xzt.dat'
*      open (20,file='j_xzt.dat')
*      write (20,20) nx,ny,nz,num_wsteps
*      write (20,30) xmin,xmax,ymin,ymax,zmin,zmax
*      do k=1,num_wsteps
*         write (20,50) k,wtime(k)
*         do j=1,nz
*            do i=1,nx
*               write (20,40) j_xzt(i,j,k)
*            enddo
*         enddo
*      enddo
*
*c     open j_yzt.dat file & write flux
*      print *,'opening writing to & closing ','j_yzt.dat'
*      open (20,file='j_yzt.dat')
*      write (20,20) nx,ny,nz,num_wsteps
*      write (20,30) xmin,xmax,ymin,ymax,zmin,zmax
*      do k=1,num_wsteps
*         write (20,50) k,wtime(k)
*         do j=1,nz
*            do i=1,ny
*               write (20,40) j_xzt(i,j,k)
*            enddo
*         enddo
*      enddo
*
*20    format (4i8)
*30    format (6e12.4)
*40    format (e12.4)
*50    format (i8,e12.4)
*
*      close(20)
*
*      return
*      end
*
********************************************************************************
********************************************************************************
********************************************************************************
