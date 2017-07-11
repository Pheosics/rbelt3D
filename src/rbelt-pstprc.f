c THIS CODE IS UNDER DEVELOPMENT
c
c***********************************************************************
c                PROGRAM MAIN FOR RBELT POST PROCESSING
c***********************************************************************
c
c weights particles based on some model and computes fluxes
c
c There are three methods used in this code:
c
c 0) method = 0 (in rbelt-pstprc.inc) gives test-particle flux
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
c 2) WEIGHT ALL PARTICLES IN SUCH A WAY AS TO MAINTAIN A MODEL FLUX IN TAIL
c REGION (specify boundary condition effective at all times): 
c a) read flux counts and compute the test particle flux in tail region
c b) weight each particle comtributing to a flux count in terms of its
c    E & t at flux count
c c) give remainder of particles weights same as the weighted particles with 
c    the same E0 and t0 (can also try to give same dist. about E0 and t0)
c d) compute weighted fluxes
c
c 3) WEIGHT PARTICLES IN TERMS OF INITIAL E-PA-TIME AND MODEL INTERPLANETARY 
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
      integer i,filenum,firstfilenum,lastfilenum,filenumstep,opt

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
      call tgrid_setup()
c     set up distribution grids
      call grid_setup()
c     initialize distribution function arrays to zero
      call init_arrays()

c     get test-particle flux and make weighting function
      if (method.le.2) then
c        get test particle flux at all tsteps
         opt=0
         fcount=0
         print *
         print *,'getting test-particle flux'
         do filenum=firstfilenum,lastfilenum,filenumstep
*            print *
            print *,'reading filenum=',filenum
            fcount=fcount+1
            call info_fileread(basename,filenum)
            call flux_fileread(basename,filenum,opt)
         enddo
c        turn flux counts into fluxes (#/cm^2-s-str-MeV)
         call wghtcnts2jxyt()
         call wghtcnts2jleat()
c        create weighting function
         call mkwghtfn()
      endif

c     read particle initial conditions file and weight each particle
      if (method.eq.1) then
         call init_pwght()
         fcount=0
         do filenum=firstfilenum,lastfilenum,filenumstep
c           reads initial condition files and weights particles
            fcount=fcount+1
            call init_fileread(basename,filenum)
         enddo
c        read 3rd step & re-weight particles (not currently used)
*         fcount=0
*         opt=4
*         do filenum=firstfilenum,lastfilenum,filenumstep
*            fcount=fcount+1
*            if (flux_out.eqv..true.) 
*     &      call flux_fileread(basename,filenum,opt)
*         enddo
      endif

c     make weighted distributions
      if (method.ge.1) then
c        re-initialize arrays
         call init_arrays()
         if (method.eq.2) call init_pwght()
         fcount=0
         print *
         print *,'getting weighted flux'
         do filenum=firstfilenum,lastfilenum,filenumstep
*            print *
            print *,'reading filenum=',filenum
            fcount=fcount+1
            call info_fileread(basename,filenum)
c           read prcp. counts & calc. weighted prcp. count dist.
            if (prcp_out.eqv..true.) then
               call prcp_fileread(basename,filenum)
               call prcp_wghtcnts2flux()
               call write_prcp()
            endif
c           read flux counts and make weighted distribution
            if (flux_out.eqv..true.) then
               call flux_fileread(basename,filenum,method)
            endif
         enddo
      endif

c tried retroactive weighting (don't bother)
*c     keep particle weights and compute weighted distributions
*c     again for retroactive weighting
*      call init_arrays()
*      fcount=0
*      print *
*      print *,'getting weighted flux'
*      do filenum=firstfilenum,lastfilenum,filenumstep
**         print *
*         print *,'reading filenum=',filenum
*         fcount=fcount+1
*         call info_fileread(basename,filenum)
*c        read prcp. counts & calc. weighted prcp. count dist.
*         if (prcp_out.eqv..true.) then
*            call prcp_fileread(basename,filenum)
*            call prcp_wghtcnts2flux()
*            call write_prcp()
*         endif
*c        read flux counts and make weighted distribution
*         if (flux_out.eqv..true.) then
*            call flux_fileread(basename,filenum,method)
*         endif
*      enddo

c     turn weighted count distribution into flux distribution
c     (#/cm^2-s-str-MeV) and write out
      if (flux_out.eqv..true.) then
         call wghtcnts2jxyt()
         call wghtcnts2jleat()
         call write_jxy()
         call write_j_goes()
         call write_jlea()
      endif

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
*      pwghtfn=10*e0**(-2.)/100.
*      pwghtfn=e0**(-3.)/1.953125
*      pwghtfn=e0**(-2.)/1.5625

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
      integer i,j,k,l

      print *
      print *,'*** in subroutine mkwght ***'

      do l=1,nt
         do k=1,na
            do j=1,ne
               do i=1,nl
                  wf(i,j,k,l)=1.0
               enddo
            enddo
         enddo
      enddo

      return
      end


********************************************************************************

      subroutine mkwghtfn()

c     calculates particle weighting function
c     at present, we modify this routine as needed

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,j,k,l,m,count,tstep
      real j_model,jprv,wfmin,wfmax,maxwf,jmdl,wf0
      parameter(maxwf=1.0E10)
      parameter(tstep=2)

      real n_lea0(nl,ne,na),j_lea0(nl,ne,na)
      common /j_lea_dist0/ n_lea0,j_lea0

      print *
      print *,'*** in subroutine mkwghtfn ***'

      if (tstep.gt.nt) then
         print *,'tstep.gt.nt in mkwghtfn_flat'
         stop
      endif

c     initialize wf
      do l=1,nt
         do k=1,na
            do j=1,ne
               do i=1,nl
                  wf(i,j,k,l)=0.
               enddo
            enddo
         enddo
      enddo

c     j_leat(i,j,k,l) must have same dimension as wf(i,j,k,l)
c     compute weight function
      count=0
      wfmin=maxwf
      wfmax=0.0
*      do k=1,na
c     wf nonzero for PAs between ~ 40 and 140 deg.
      do k=3,14
         do j=1,ne

c           2.E10 is j0
c           1.e6 convers to per MeV 
c           1.0E6**(-1.8) converts energy from MeV to eV
*            jmdl=316979.*ke(j)**(-1.8)

c           ESA model j_perp at L=7.5
            jmdl=15137.*ke(j)**(-3.54)

            do i=1,nl

c              get time-averaged test-particle flux 
               jprv = 0.0
               do l=1,fullsteps
                  jprv = jprv + j_leat(i,j,k,l)
               enddo
               jprv = jprv/real(fullsteps)

c              if the test-particle flux is not zero, we can compute weight
               if (jprv.gt.0.0) then
                  wf0=jmdl/jprv
c                 in this case, weighting fuction is independent of time
                  do m=1,nt
                     wf(i,j,k,m)= wf0
                     if (wf(i,j,k,1).lt.wfmin) wfmin=wf(i,j,k,m)
                     if (wf(i,j,k,1).gt.wfmax) wfmax=wf(i,j,k,m)
                  enddo
               else
                  count=count+1
                  do m=1,fullsteps
                     wf(i,j,k,m)= 0.0
                  enddo
               endif

*               print *,'alpha,ke,ls,jmdl,javg,wf=',
*     &         acos(ca(k))*raddeg,ke(j),ls(i),jmdl,jprv,wf(i,j,k,1)

            enddo

            print *,'alpha,ke,ls,jmdl,javg,wf=',
     &      acos(ca(k))*raddeg,ke(j),ls(5),jmdl,jprv,wf(5,j,k,1)

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

c     These are not curently used in pstprc:
      real charge_rb,charge_sign_rb,m0_rb,b0_rb,dthalt
      real tstep_lrntz,dx_max,tstep_gc,dx_max_gc
      real go2lrntz,go2gc,dtgo2gc,init_t,epa_min,epa_max
      integer nstep,num_particles_rb,max_yout,num_flts,num_ints
      integer max_wsteps,firstfilenum,seed,dist_seed,flag0
      integer num_wsteps_rb,nx_rb,ny_rb,nz_rb,nt_rb,lcalc
      character*80 basename_rb
      logical flag_switch
c     remaining are declared in pstprc.inc & used as globals

      filename=infofile(basename,filenum)
*      print *
*      print *,'opening ',filename(1:lnblnk(filename))
c     open info file & read
      open (12,file=filename(1:lnblnk(filename)))
c     from rbelt-grid.inc
      read (12,10) nx_rb
      read (12,10) ny_rb
      read (12,10) nz_rb
      read (12,10) nt_rb
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
      read (12,20) exp
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
*      print *,'closing ',filename(1:lnblnk(filename))
      close(12)

10    format (16x,i10)
20    format (16x,e12.4)
30    format (16x,a80)
40    format (16x,l1)

*      if (filenum.eq.firstfilenum) then
*      print*,'nx=',nx_rb
*      print*,'ny=',ny_rb
*      print*,'nz=',nz_rb
*      print*,'nt=',nt
*      print*,'nstep=',nstep
*      print*,'charge=',charge
*      print*,'charge_sign=',charge_sign
*      print*,'m0=',m0
*      print*,'b0=',b0
*      print*,'num_particles=',num_particles_rb
*      print*,'max_yout=',max_yout
*      print*,'num_flts=',num_flts
*      print*,'num_ints=',num_ints
*      print*,'max_wsteps=',max_wsteps
*      print*,'basename=',basename
*      print*,'firstfilenum=',firstfilenum
*      print*,'rmin=',rmin
*      print*,'rmax=',rmax
*      print*,'tmax=',tmax
*      print*,'dthalt=',dthalt
*      print*,'flux_dt=',flux_dt
*      print*,'init_twrite=',init_twrite
*      print*,'dtwrite=',dtwrite
*      print*,'binio=',binio
*      print*,'tstep_lrntz=',tstep_lrntz
*      print*,'dx_max=',dx_max
*      print*,'tstep_gc=',tstep_gc
*      print*,'dx_max_gc=',dx_max_gc
*      print*,'go2lrntz=',go2lrntz
*      print*,'go2gc=',go2gc
*      print*,'dtgo2gc=',dtgo2gc
*      print*,'seed=',seed
*      print*,'dist_seed=',dist_seed
*      print*,'dt_dist=',dt_dist
*      print*,'init_t=',init_t
*      print*,'emin=',emin
*      print*,'emax=',emax
*      print*,'epa_min=',epa_min
*      print*,'epa_max=',epa_max
*      print*,'lmin=',lmin
*      print*,'lmax=',lmax
*      print*,'radius=',radius
*      print*,'exp=',exp
*      print*,'factor=',factor
*      print*,'flag0=',flag0
*      print*,'flag_switch=',flag_switch
*      print*,'initdist=',initdist
*      print*,'year0=',year0
*      print*,'doy0=',doy0
*      print*,'hour0=',hour0
*      print*,'min0=',min0
*      print*,'sec0=',sec0
*      print*,'num_wsteps=',num_wsteps_rb
*      print*,'wlines_dist=',wlines_dist
*      print*,'wlines_flux=',wlines_flux
*      print*,'wlines_prcp=',wlines_prcp
*      print*,'wlines_prcp=',wlines_cone
*      endif

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

      subroutine tgrid_setup() 

      implicit none
      include 'rbelt-pstprc.inc'
      integer i,tstep,count
      real sum

      print *
      print *,'*** in subroutine tgrid_setup ***'

*c     set up rbelt wtime array
*      do i=1,num_wsteps
*         wtime(i)=init_twrite+(i-1)*dtwrite
**         print *,'i,wtime',i,wtime(i)
*      enddo

c     set up pstprc time array
      if (prcp_out.eqv..true.) then
         do i=1,nt
            tg(i)=init_time+(i-1)*dtstep
*	    print *,'i,time',i,tg(i)
         enddo
      else
         tstep=0
         sum=0
         count = 0
         do i=1,num_wsteps
            wtime(i)=init_twrite+(i-1)*dtwrite
*            print *,'i,wtime',i,wtime(i)
c           get time corresonding to tstep
            sum=sum+wtime(i)
            count=count+1
            if (count.eq.stepavg) then
               tstep=tstep+1
               tg(tstep)=sum/real(count)
	       print *,'tstep, #of averaged steps, avg. time',
     &         tstep,count,tg(tstep)
               sum=0
               count = 0
            endif
         enddo
         tstep=tstep+1
         tg(tstep)=sum/real(count)
	 print *,'tstep, #of averaged steps, avg. time',
     &   tstep,num_wsteps-i+count+1,tg(tstep)
      endif

      return
      end

********************************************************************************

      subroutine grid_setup()

c     set up post processing grids

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l
      real dx,dy,dz,dtht,dphi,dl,de,da

      print *
      print *,'*** in subroutine grid_setup ***'

      call mk1dgrd(nl,lsmin,lsmax,dl,ls)
      call mk1dgrd(ne,kemin,kemax,de,ke)
      call mk1dgrd(na,camin,camax,da,ca)
      print *
      print *,'L,energy,cos(PA) grid'
      print *
      do i=1,nl
         print *,'i,ls(i)=',i,ls(i)
      enddo
      print *
      do i=1,ne
         print *,'i,ke(i),jmdl,jmdl_int=',
     &   i,ke(i),316979.*ke(i)**(-1.8),396223.*ke(i)**(-.8)
      enddo
      print *
      do i=1,na
         print *,'i,ca(i),alpha=',i,ca(i),acos(ca(i))*raddeg
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
      print *
      print *,'xmin,xmax=',xmin,xmax
      print *,'ymin,ymax=',ymin,ymax

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

      subroutine init_arrays()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l

c     initialize flux distribution arrays to zero
      do k=1,nt
         do j=1,nphi
            do i=1,ntht
               j_tpt(i,j,k)=0.
               n_tpt(i,j,k)=0.
            enddo
         enddo
      enddo

      do k=1,nt
         do j=1,ny
            do i=1,nx
               j_xyt(i,j,k)=0.
               n_xyt(i,j,k)=0.
            enddo
         enddo
      enddo

      do l=1,nt
         do k=1,na
            do j=1,ne
               do i=1,nl
                  j_leat(i,j,k,l)=0.
                  n_leat(i,j,k,l)=0.
               enddo
            enddo
         enddo
      enddo

      return
      end

********************************************************************************

      subroutine init_pwght()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l

      do j=1,num_files
         do i=1,num_particles
            pwght(i,j)=0.
            pwghtflag(i,j)=0
         enddo
      enddo

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
      real x0,y0,z0,l0,e0,a0,ca0,t0,pw,fld2pt4d0

      print *
      print *,'*** in subroutine init_fileread ***'

c     name & open input file
      filename=initfile(basename,filenum)
      print *
      print *,'opening ',filename(1:lnblnk(filename))
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
         pwght(n,fcount)=
     &   fld2pt4d0(n,l0,e0,ca0,t0,nl,ne,na,nt,ls,ke,ca,tg,wf)
         pwghtflag(n,fcount)=1

         if (pwght(n,fcount) .gt. 100000000) then
            print *,'pwght(n,fcount) .gt. 100000000)'
            print *,'pwght(n,fcount)=',pwght(n,fcount)
            print *,'n,l0,e0,ca0,t0=',n,l0,e0,ca0,t0
            stop
         endif

*         print *
*         print *,'l0,e0,ca0',l0,e0,ca0
*         print *,'n,fcount,pwght=',n,fcount,pwght(n,fcount)

      enddo

      print *,'closing ',filename(1:lnblnk(filename))
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
      integer lnblnk,k,pnum,count(nt),filenum,tstep
      integer year,doy,hour,min,sec
      real t0,x0,y0,z0,e0,a0,ceta,il,ie,ia,it,tht0,phi0,pw,fld2pt4d0
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
      do k=1,nt
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
         call RECALC_08 (year,doy,hour,min,sec,-400.,0.,0.)
*         print *,'done with RECALC: tilt=',atan2(SPS,CPS)*raddeg
*         print *,'SM: x0,y0,z0=',x0,y0,z0
         call SMGSW_08 (x0,y0,z0,XGSM,YGSM,ZGSM,1)
         call GEOGSW_08 (x0,y0,z0,XGSM,YGSM,ZGSM,-1)
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
            pw=fld2pt4d0(pnum,il,ie,ia,it,nl,ne,na,nt,ls,ke,ca,tg,wf)
*             pw=1.0
         endif

c        get tstep corresponding to t0
c        this is set up for tsteps = num_wsteps and needs to be fixed -- 
c        see prcp_wghtcnts2flux()
c        we can assume uniform time steps
c        for now we use 0th order interp in time

         if ((t0.ge.tg(1)).and.(t0.lt.tg(nt))) then
            tstep=int((t0-init_time+dtstep/2.)/dtstep)+1
            count(tstep)=count(tstep)+1
         else
            tstep=0
         endif

         call pt2fld2d0t(pnum,tht0,phi0,tstep,ntht,nphi,
     &   nt,tht,phi,pw,n_tpt,j_tpt)
c
*         call pt2fld3d0(pnum,tht0,phi0,t0,ntht,nphi,
*     &   num_wsteps,tht,phi,wtime,pw,n_tpt,j_tpt)

      enddo

      print *,'closing ',filename(1:lnblnk(filename))
      close(18)
      return
      end

********************************************************************************

      subroutine flux_fileread(basename,filenum,opt)

      implicit none
      include 'rbelt-pstprc.inc'
      include 'rbelt-geopak.inc'

      character*80 basename,fluxfile,filename
      integer lnblnk,i,k,pnum,count(nt),filenum,wstep,tstep,
     &opt,status,grdpos0
      real x0,y0,z0,lshell,energy,alpha,calpha,eta,ceta,fld2pt4d0,
     &e0,a0,ca0,th0,cth0,t,t0,pw,pwghtfn,bfield,azmth

*      print *
*      print *,'*** in subroutine flux_fileread ***'

c     name & open input file
      filename=fluxfile(basename,filenum)
*      print *,'opening ',filename(1:lnblnk(filename))
      if (binio.eqv..true.) then
         open (18,file=filename(1:lnblnk(filename)),
     &   form='unformatted')
      else
         open (18,file=filename(1:lnblnk(filename)))
      endif

c     read file
*      print *,'wlines_flux=',wlines_flux
      do k=1,wlines_flux
         if (binio.eqv..true.) then
            if (method.eq.1) then
               read (18) pnum,wstep,x0,y0,eta,energy,alpha,bfield
            else
               read (18) pnum,wstep,x0,y0,eta,energy,alpha,bfield,
     &         e0,a0,t0
            endif
         else
            if (method.eq.1) then
               read (18,10) pnum,wstep,x0,y0,eta,energy,alpha,bfield
            else
               read (18,20)pnum,wstep,x0,y0,eta,energy,alpha,bfield,
     &         e0,a0,t0
            endif
         endif
10       format (2i8,6e12.4)
20       format (2i8,9e12.4)

*      print *,'binio=',binio
*      print *,'pnum,wstep,x0,y0,eta,energy,alpha,bfield,e0,a0,t0=',
*     &pnum,wstep,x0,y0,eta,energy,alpha,bfield,e0,a0,t0
*      stop

c        get time step & time (what's all this!)
         tstep = int(real((wstep)-.5)/real(stepavg))+1
         t=tg(tstep)

c        get cos of alpha and eta
         calpha=cos(alpha/raddeg)
         ceta=cos(eta/raddeg)
         ca0=cos(a0/raddeg)

c        dipole L and azimuthal angle
         lshell=sqrt(x0*x0+y0*y0)
         azmth=atan2(y0,x0)
         if (azmth.lt.0.0) azmth=2*pi+azmth

c        get particle weight
         if (opt.le.3) then
c           get particle weight / cos(eta)
            status = 0
c           *** note, here we use cos(alpha_o) grid to approximate
c           cos of angle of incidence with equatorial plane
            i=grdpos0(ceta,na,ca,status)
            if (status.eq.0) then

c              opt=0 to get test-particle flux
               if (opt.eq.0) then
                     pw=1.0/abs(ca(i))

c              opt=1 for previously assigned weights
               elseif (opt.eq.1) then
                  pw=pwght(pnum,fcount)/abs(ca(i))

c              opt=2 weight particles here
               elseif (opt.eq.2) then

c                 if in boundary region, weight particles
*                  if ((lshell.ge.lmin).and.(lshell.le.lmax)) then
                  if ((lshell.ge.7.5).and.(lshell.le.8.5)) then
c                 also tried no re-weighting (don't use this)
*                  if ((lshell.ge.7.0).and.(lshell.le.8.0).and.
*     &            (pwghtflag(pnum,fcount).eq.0)) then

                     pwght(pnum,fcount)=fld2pt4d0(pnum,lshell,energy,
     &               calpha,t,nl,ne,na,nt,ls,ke,ca,tg,wf)
                     pwghtflag(pnum,fcount)=1
                     pw=pwght(pnum,fcount)/abs(ca(i))

c                 if previously in boundary region, use current weight
                  elseif (pwghtflag(pnum,fcount).eq.1) then
                     pw=pwght(pnum,fcount)/abs(ca(i))

c                 if never counted in boundary region, 
c                 weight here in terms of initial L,E,alpha_o,t
c                 in boundary region -- only use if we launch
c                 all particles in boundary region,
c                 otherwise set pw=0.0 here
c                 (note that, to get the 'right' answer
c                 we should weight all particles in boundary region,
c                 including those not counted (i.e., because thier 
c                 bounce phase is such that they do not pass through)
c                 equatorial plane durring flux counting interval)
c                 but if most are counted, it should not make much 
c                 difference)
                  else
c                    (used this for _pwi and _pwri)
*                     pwght(pnum,fcount)=fld2pt4d0(pnum,7.5,e0,ca0,t,
*     &               nl,ne,na,nt,ls,ke,ca,tg,wf)/abs(ca(i))
c                    (used this for _pw0 and _pwr0)
                     pw=0.0
                     numuwc=numuwc+1
                  endif

c why do _pwri _pw0 and _pwr0 give approx same while _pwi has bad spikes ?

c              opt=3 compute weight using analytical function
               elseif (opt.eq.3) then
                     pw=pwghtfn(pnum,e0,ca0,t)/abs(ca(i))

               endif

c              add flux count to distributions

c              add weighted count to L-E-PA-time distribution
               call pt2fld3d0t(pnum,lshell,energy,calpha,tstep,nl,ne,
     &         na,nt,ls,ke,ca,pw,n_leat,j_leat)

               if (energy.ge.elwr) then
c                 add weighted count to omni-integral x-y-time distribution
                  call pt2fld2d0t
     &            (pnum,x0,y0,tstep,nx,ny,nt,x,y,pw,n_xyt,j_xyt)
               endif

*            else
*               print *,'grdpos0 returns status.ne.0'
*               print *,'ceta,ca(1),ca(na),status=',
*     &         ceta,ca(1),ca(na),status
*               stop
            endif

         endif

      enddo

*      print *,'closing ',filename(1:lnblnk(filename))
      close(18)
      return
      end

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
      do k=1,nt
         call grd_delta(k,nt,tg,upr,lwr,dt)
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

      integer i,j,k,count(nt),totcnts(nt)
      real dxyt,upr,lwr,dx,dy

      print *
      print *,'*** in subroutine wghtcnts2jxyt ***'

c     expected flux (from total particles launched into sphere over dt_dist time interval):
*      print *
*      print *,'J_omni/2 = #/RE^2-sec/2 launched into sphere='
*      print *,num_particles*num_files/pi/radius/radius/dt_dist/2.

c      calculate #/RE^2-sec through surface (*NOTE*, this is J_omni/2)
      do k=1,nt
*      do k=1,1
         count(k)=0
         totcnts(k)=0

         do j=1,ny
            call grd_delta(j,ny,y,upr,lwr,dy)
            do i=1,nx
               call grd_delta(i,nx,x,upr,lwr,dx)
               if (n_xyt(i,j,k).gt.0.0) count=count+1
               totcnts=totcnts+n_xyt(i,j,k)

*               if ((i.eq.156).and.(j.eq.76)) then
*                  print *
*                  print *,'i,j,k,n_xyt(i,j,k),j_xyt(i,j,k)=',
*     &            i,j,k,n_xyt(i,j,k),j_xyt(i,j,k)
*               endif

               if (n_xyt(i,j,k).ge.mincnt) then
                  dxyt = abs(1.0/dx/dy/flux_dt/stepavg)
*                  print *,'j_xyt(i,j,k)=',j_xyt(i,j,k)
                  j_xyt(i,j,k)=j_xyt(i,j,k)*dxyt
               else
                  j_xyt(i,j,k)=0.0
               endif

*               if ((i.eq.156).and.(j.eq.76)) then
*                  print *,'j_xyt(i,j,k)=',j_xyt(i,j,k)
*               endif

            enddo
         enddo

*         print *
         print *,'k, #non-zero cells, totcnts, avg./cell=',
     &   k,count(k),totcnts(k),real(totcnts(k))/count(k)
*         do i=150,160
*            print *,'   i,x(i),j_xyt(i,76,k)=',i,x(i),j_xyt(i,76,k)
*         enddo

      enddo



      return
      end

********************************************************************************

      subroutine wghtcnts2jleat()

      implicit none
      include 'rbelt-pstprc.inc'

      integer i,j,k,l
      real dleat,upr,lwr,dl,de,da,fcntmin,fcntmax,nzeros
      real javg,jmdl

      print *
      print *,'*** in subroutine wghtcnts2jleat ***'

      fcntmin=10000000000.
      fcntmax=0.

      do l=1,nt
         do k=1,na
            call grd_delta(k,na,ca,upr,lwr,da)
c           divide by cos(eta) below moved to flux_fileread
*            da=abs(ca(k)*sin(acos(ca(k)))*(acos(lwr)-acos(upr)))
            da=abs(sin(acos(ca(k)))*(acos(lwr)-acos(upr)))
            do j=1,ne
               call grd_delta(j,ne,ke,upr,lwr,de)
               do i=1,nl
                  call grd_delta(i,nl,ls,upr,lwr,dl)
c                 dipole L-shell

c                 TEMPORY MULT. BY 80/360 !!
*                  dl=pi*(upr*upr-lwr*lwr)*.222222

                  dl=pi*(upr*upr-lwr*lwr)

                  dleat = abs(1.0/2./pi/dl/da/de/flux_dt/stepavg)

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

c     compute time averaged flux in boundary region 
c     for comparison with model
*      do k=1,na
      do k=3,14
         do j=1,ne
*            jmdl=317000.*ke(j)**(-1.8)]
            jmdl=15137.*ke(j)**(-3.54)
c           at L=7.5
            i=5
            javg = 0.0
            do l=1,nt-1
               javg = javg + j_leat(i,j,k,l)
            enddo
            javg = javg/real(nt-1)
            print *,'alpha,ke,ls,jmdl,javg=',
     &      acos(ca(k))*raddeg,ke(j),ls(i),jmdl,javg
         enddo
      enddo

      print *
*      print *,'fcntmin,fcntmax,nzeros=',fcntmin,fcntmax,nzeros
      print *,'num. unweighted counts=',numuwc

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
      write (20,20) ntht,nphi,nt
      write (20,30) thtmin*raddeg,thtmax*raddeg,
     & phimin*raddeg,phimax*raddeg
c     eventually we should remove k,tht & phi as they are not needed
      do k=1,nt
         write (20,40) tg(k)
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
      write (20,20) nx,ny,fullsteps
      write (20,30) xmin,xmax,ymin,ymax
      do k=1,fullsteps
*         print *,'k,tg(k)=',k,tg(k)
         write (20,40) tg(k)
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

      subroutine write_j_goes_old()
      implicit none
      include 'rbelt-pstprc.inc'
      include 'rbelt-geopak_08.inc'
      integer i,j,k,year,doy,hour,min,sec
      real RGEO13,PGEO13,AGEO13,RGEO15,PGEO15,AGEO15
      real XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,XSM,YSM,ZSM
      real j_goes13,j_goes15,fld2pt2dt

c     goes-13 at 75 deg. west
      parameter(RGEO13=6.4,PGEO13=90./raddeg,AGEO13=285./raddeg)
c     goes-15 at 135 deg. west
      parameter(RGEO15=6.4,PGEO15=90./raddeg,AGEO15=225./raddeg)

      print *
      print *,'*** in subroutine write_j_goes ***'
c     open j_xyt.dat file & write flux
      print *,'opening writing to & closing ','j_xyt.dat'
      open (20,file='j_goes.dat')
      write (20,10) year0,doy0,hour0,min0,sec0
      do k=1,fullsteps
*         print *,'k,tg(k)=',k,tg(k)

c        get date and time of count
         year=year0
         doy=doy0
         hour=hour0
         min=min0
         sec=sec0+tg(k)
         call reduce_time2(year,doy,hour,min,sec)
c        call recalc and rotate to SM coordinates
*         print *
*         print *,'calling RECALC w/ YEAR,DOY,HOUR,MIN,SEC=',
*     &   year,doy,hour,min,sec
         call RECALC_08(year,doy,hour,min,sec,-400.0,0.,0.)
*         print *,'done with RECALC: tilt=',atan2(SPS,CPS)*raddeg

         call sphr2cart(RGEO13,PGEO13,AGEO13,XGEO,YGEO,ZGEO)
         call GEOGSW_08 (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1)
         call SMGSW_08 (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,-1)
*         print *,'XGEO,YGEO,ZGEO=',XGEO,YGEO,ZGEO
*         print *,'GOES-13: XSM,YSM,ZSM=',XSM,YSM,ZSM

*         xsm3(1)=XSM
*         xsm3(2)=YSM
*         xsm3(3)=ZSM
*         call trace2eqtr(t0,xsm3,Beqtr,ifail)

         j_goes13=fld2pt2dt(XSM,YSM,k,nx,ny,nt,x,y,j_xyt)

         call sphr2cart(RGEO15,PGEO15,AGEO15,XGEO,YGEO,ZGEO)
         call GEOGSW_08 (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1)
         call SMGSW_08 (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,-1)
*         print *,'XGEO,YGEO,ZGEO=',XGEO,YGEO,ZGEO
*         print *,'GOES-15: XSM,YSM,ZSM=',XSM,YSM,ZSM

         j_goes15=fld2pt2dt(XSM,YSM,k,nx,ny,nt,x,y,j_xyt)

         write (20,40) tg(k),j_goes13,j_goes15

      enddo
10    format (5i8)
40    format (3e12.4)
      close(20)
      return
      end

********************************************************************************

      subroutine write_j_goes()
      implicit none
      include 'rbelt-pstprc.inc'
      include 'rbelt-geopak_08.inc'
      integer i,j,k,l,year,doy,hour,min,sec
      real RGEO13,PGEO13,AGEO13,RGEO15,PGEO15,AGEO15
      real XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,XSM,YSM,ZSM
      real j_goes13,j_goes15,fld2pt2dt

c     goes-13 at 75 deg. west
      parameter(RGEO13=6.4,PGEO13=90./raddeg,AGEO13=285./raddeg)
c     goes-15 at 135 deg. west
      parameter(RGEO15=6.4,PGEO15=90./raddeg,AGEO15=225./raddeg)

      integer ifail13,ifail15,iyear,idoy,ihour,imin,isec,nlines
      real t,xsm13,ysm13,bloc13,Beqtr13,xsm15,ysm15,bloc15,Beqtr15
      real ttg,xsm13tg,ysm13tg,bloc13tg,Beqtr13tg,xsm15tg,ysm15tg,
     &bloc15tg,Beqtr15tg
      real t_0,xsm13_0,ysm13_0,bloc13_0,Beqtr13_0,xsm15_0,ysm15_0,
     &bloc15_0,Beqtr15_0
      real dt,wi,wip

      print *
      print *,'*** in subroutine write_j_goes ***'

c     open j_goes.dat file & write flux
      print *,'opening writing to & closing ','j_goes.dat'
      open (20,file='j_goes.dat')
      write (20,10) year0,doy0,hour0,min0,sec0

c     read rbelt GOES output file and get GOES location
c     mapped to SM equatorial plane, also, B_local and B_equatorial
      open (18,file='oct12lfm-goes.txt')
      read (18,*) iyear,idoy,ihour,imin,isec

      if ((iyear.ne.year0).or.(idoy.ne.doy0).or.(ihour.ne.hour0).or.
     &(imin.ne.min0).or.(isec.ne.sec0)) then
         print *,'itime.ne.time0 in subroutine write_j_goes'
         stop
      endif

      read (18,*) nlines
      if (nlines.ne.num_wsteps) then
         print *,'nlines.ne.num_wsteps in subroutine write_j_goes'
         stop
      endif
      read (18,*) t,xsm13,ysm13,bloc13,Beqtr13,ifail13,
     &xsm15,ysm15,bloc15,Beqtr15,ifail15
      if (t.ge.tg(1)) then
         print *,'problem: t0.ge.tg(1) in subroutine write_j_goes'
         stop
      endif
      t_0=t
      xsm13_0=xsm13
      ysm13_0=ysm13
      bloc13_0=bloc13
      Beqtr13_0=Beqtr13
      xsm15_0=xsm15
      ysm15_0=ysm15
      bloc15_0=bloc15
      Beqtr15_0=Beqtr15

      k=1
      do l=2,num_wsteps

         read (18,*) t,xsm13,ysm13,bloc13,Beqtr13,ifail13,
     &   xsm15,ysm15,bloc15,Beqtr15,ifail15

c        if we pass a pstprt time grid value
         if (t.ge.tg(k).and.(k.lt.nt)) then

c           linear interpolate from rbelt write step times to pstprt time grid
            dt=t-t_0
            wi = (t-tg(k))/dt
            wip = 1.-wi
            xsm13tg=wi*xsm13_0+wip*xsm13
            ysm13tg=wi*ysm13_0+wip*ysm13
            bloc13tg=wi*bloc13_0+wip*bloc13
            Beqtr13tg=wi*Beqtr13_0+wip*Beqtr13
            xsm15tg=wi*xsm15_0+wip*xsm15
            ysm15tg=wi*ysm15_0+wip*ysm15
            bloc15tg=wi*bloc15_0+wip*bloc15
            Beqtr15tg=wi*Beqtr15_0+wip*Beqtr15

*            print *,'Beqtr15,bloc15=',Beqtr15,bloc15
*            print *,'xsm13tg,ysm13tg=',xsm13tg,ysm13tg
*            print *,'xsm15tg,ysm15tg=',xsm15tg,ysm15tg

            j_goes13=fld2pt2dt(xsm13tg,ysm13tg,k,nx,ny,nt,x,y,j_xyt)
            j_goes15=fld2pt2dt(xsm15tg,ysm15tg,k,nx,ny,nt,x,y,j_xyt)

            write (20,40) tg(k),j_goes13,j_goes15

*c           get date and time of count
*            year=year0
*            doy=doy0
*            hour=hour0
*            min=min0
*            sec=sec0+tg(k)
*            call reduce_time2(year,doy,hour,min,sec)
*c           call recalc and rotate to SM coordinates
**            print *
**            print *,'calling RECALC w/ YEAR,DOY,HOUR,MIN,SEC=',
**     &      year,doy,hour,min,sec
*            call RECALC_08(year,doy,hour,min,sec,-400.0,0.,0.)
**            print *,'done with RECALC: tilt=',atan2(SPS,CPS)*raddeg
*c           GEO to SM coordinates
*            call sphr2cart(rgeo13,pgeo13,ageo13,XGEO,YGEO,ZGEO)
*            call GEOGSW_08 (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1)
*            call SMGSW_08 (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,-1)
**            print *,'xsm13,ysm13=',xsm,ysm
**            print *,'GOES13 map versus proj.diff. (RE)=',
**     &      sqrt((XSM-xsm13tg)**2+(YSM-ysm13tg)**2)
*            call sphr2cart(rgeo15,pgeo15,ageo15,XGEO,YGEO,ZGEO)
*            call GEOGSW_08 (XGEO,YGEO,ZGEO,XGSM,YGSM,ZGSM,1)
*            call SMGSW_08 (XSM,YSM,ZSM,XGSM,YGSM,ZGSM,-1)
*            print *,'xsm15,ysm15=',xsm,ysm
**            print *,'GOES15 map versus proj.diff. (RE)=',
**     &      sqrt((XSM-xsm15tg)**2+(YSM-ysm15tg)**2)
**            print *,'proj:XSM,YSM, map:XSM,YSM',XSM,YSM,xsm15tg,ysm15tg
*            print *,'zsm, r_proj.,r_map=',ZSM,
*     &      sqrt(XSM**2+YSM**2),sqrt(xsm15tg**2+ysm15tg**2)

            if (k.lt.nt) k=k+1

         endif
         t_0=t
         xsm13_0=xsm13
         ysm13_0=ysm13
         bloc13_0=bloc13
         Beqtr13_0=Beqtr13
         xsm15_0=xsm15
         ysm15_0=ysm15
         bloc15_0=bloc15
         Beqtr15_0=Beqtr15

      enddo

10    format (5i8)
40    format (3e12.4)
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
      write (20,20) nl,ne,na,fullsteps
      write (20,30) lsmin,lsmax,kemin,kemax,camin,camax

      do l=1,fullsteps
         write (20,40) tg(l)
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
         print *,'k,ea(k),j=',k,real(acos(ea(k))*raddeg),
     &   real(sum)
         sum2=sum2+sum
      enddo
      print *,'total j counts=',real(sum2)
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
         print *,'j,ke(j),j=',j,real(ke(j)),real(sum)
         sum2=sum2+sum
      enddo
      print *,'total j counts=',real(sum2)
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
         print *,'i,ls(i),j=',i,real(ls(i)),real(sum)
         sum2=sum2+sum
      enddo
      print *,'total j counts=',real(sum2)
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
      real camax,camin,ca(na)
      integer ke_bin(ne),ca_bin(na)
      real de,da,hda,dflux(na)

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
         ca(k) = camin + (real(k)-0.5)*da
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
*c     expected flux (from total particles launched into sphere over dt_dist time interval):
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
********************************************************************************
********************************************************************************
