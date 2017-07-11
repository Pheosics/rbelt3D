
************************************************************************

      subroutine load_field
     &(basename,timestep,firststep,laststep,gridstep)

c     puts TS05 fields at UT ref. time into bxd,byd & bzd 
c     and zero into exd,eyd & ezd. also puts time into tgr(1)
c     note, no rbelt normalization in this code

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-geopack_08.inc'
      include 'rbelt-ut.inc'
      include 'rbelt-ts05.inc'

c     input vars.
      character*(*) basename
      integer timestep,firststep,laststep,gridstep
c     parameters
      real inbound,outbound,b0_nt
      parameter (b0_nt=b0/ffactor/ntg)

c *************************************************************************
c *************************************************************************

c     all rbelt grid points inside radial inbound get zero
      parameter (inbound=1.0-1.733*.2)  !inner boundary
c     all rbelt grid points outside radial outbound get zero
      parameter (outbound=16.0+1.733*.2)  !outer boundary

c *************************************************************************
c *************************************************************************

c     additional lfm_sample variables:
      integer i,j,k,l,m,sfr,sec_frm_ref
      real r,r2,r5,bxsm,bysm,bzsm

c     geopack & Tsyganenko routines
      real XGSM,YGSM,ZGSM
      real HXGSM,HYGSM,HZGSM,BXGSM,BYGSM,BZGSM
      real BX_sgl,BY_sgl,BZ_sgl
      integer mode
      integer year,doy,hour,min,sec

c *************************************************************************
c *************************************************************************

      print *
      print *,'*** in subroutine load_field (ts05_samp version) ***'
      print *,'nx,ny,nz=',nx,ny,nz
      print *,'xgmin,xgmax=',xgmin,xgmax
      print *,'ygmin,ygmax=',ygmin,ygmax
      print *,'zgmin,zgmax=',zgmin,zgmax
      print *
      print *,'parameters set here are:'
      print *,'b0_nt (nT)=',b0_nt
      print *,'inbound,outbound=',inbound,outbound
      print *,'tsy,tsy_input=',tsy,' ',tsy_input
      print *,'max_table=',max_table

c     read and write loop starts here
c     initialize loop
      BX_sgl=0
      BY_sgl=0
      BZ_sgl=0
      m=0
      do l=gridstep,(gridstep+laststep-firststep)

c        put time on time grid here
c        seconds from reference time is timestep*dsec
c        tgr(1) is scaled/normalized after fields are loaded
         sfr=(timestep-1+m)*dsec
	 m=m+1
         tgr(l)=sfr*tfactor

c        print requested date & time
         call sfr2date(sfr,year,doy,hour,min,sec)
         print *
         print *,'time grid step, time step, sfr =',l,timestep+m-1,sfr
         print *,'requested year,doy,hour,min,sec',year,doy,hour,min,sec

c        call geopak recalc
         call rbelt_recalc(sfr)

c        get TS05 input parameters
         if (tsy.eqv..true.) then
            call get_ts05_params(sfr)
            if (abs(tilt0-PSI).gt.0.01) then
               print *,'WARNING! TS05 input file and RECALC'
               print *,'tilt angles do not agree'
*               print *,'(probably due to interp. between table records)'
*               print *,'(using RECALC PSI)'
	       print *,'tilt0,PSI=',tilt0,PSI
            endif
         endif

c        loop over grid positions
         print *
         print *,'putting fields on grid...'
         do k=1,nz
*           print *,'k=',k
           do j=1,ny
             do i=1,nx
               r2=xgr(i)*xgr(i)+ygr(j)*ygr(j)+zgr(k)*zgr(k)
               r=sqrt(r2)
               r5=r2*r2*r
               if ((r.le.outbound).and.(r.ge.inbound)) then
                 mode=1
                 CALL SMGSW_08(xgr(i),ygr(j),zgr(k),XGSM,YGSM,ZGSM,mode)
                 CALL IGRF_GSW_08(XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)
                 if (tsy.eqv..true.) then
                   CALL T04_s(IOPT,PARMOD,PSI,XGSM,YGSM,ZGSM,
     &             BX_sgl,BY_sgl,BZ_sgl)
                 endif
                 BXGSM=HXGSM+BX_sgl
                 BYGSM=HYGSM+BY_sgl
                 BZGSM=HZGSM+BZ_sgl
                 mode=-1
                 CALL SMGSW_08(bxsm,bysm,bzsm,BXGSM,BYGSM,BZGSM,mode)
c                dipole field gets removed
                 bxd(i,j,k,gridstep)=(bxsm+3.*b0_nt*xgr(i)*zgr(k)/r5)
                 byd(i,j,k,gridstep)=(bysm+3.*b0_nt*ygr(j)*zgr(k)/r5)
                 bzd(i,j,k,gridstep)=
     &           (bzsm+3.*b0_nt*zgr(k)*zgr(k)/r5-b0_nt/r2/r)
                 exd(i,j,k,gridstep)=0.
                 eyd(i,j,k,gridstep)=0.
                 ezd(i,j,k,gridstep)=0.
               else
                 bxd(i,j,k,gridstep) = 0.0
                 byd(i,j,k,gridstep) = 0.0
                 bzd(i,j,k,gridstep) = 0.0
                 exd(i,j,k,gridstep) = 0.0
                 eyd(i,j,k,gridstep) = 0.0
                 ezd(i,j,k,gridstep) = 0.0
               endif
             enddo
           enddo
         enddo
c        uncomment to sanity check fields along x axis
*         print *
*         print *,'purturbation fields'
*         do i=1,nx
**           print *,'x,bz=',x(i),
**           bz(i,(ny+1)/2,(nz+1)/2,n)/bscal
*           print *,xgr(i),bzd(i,(ny+1)/2,(nz+1)/2,gridstep)
*         enddo
         print *
         print *,'total field'
         do i=1,nx/2
*           print *,'x,bz=',x(i),
*           bz(i,(ny+1)/2,(nz+1)/2,n)/bscal
           print *,xgr(i),
     &     (bzd(i,(ny+1)/2,(nz+1)/2,gridstep)+b0_nt/abs(xgr(i)**3))
         enddo
         do i=nx/2+2,nx
*           print *,'x,bz=',x(i),
*           bz(i,(ny+1)/2,(nz+1)/2,n)/bscal
           print *,xgr(i),
     &     (bzd(i,(ny+1)/2,(nz+1)/2,gridstep)+b0_nt/abs(xgr(i)**3))
         enddo

      enddo

*      stop

      end

************************************************************************

      subroutine load_grid(basename)
      implicit none
      character*80 basename
      include 'rbelt-grid.inc'
      integer i
      real nx_tmp,ny_tmp,nz_tmp
      
c     set up grid
      nx_tmp = nx
      ny_tmp = ny
      nz_tmp = nz
      if (nx.eq.1) nx_tmp=2
      if (ny.eq.1) ny_tmp=2
      if (nz.eq.1) nz_tmp=2
      do i=1,nx
         xgr(i) = xgmin + ((xgmax-xgmin)/(nx_tmp-1))*(i-1)
*         print *,'x(',i,')=',x(i)
      enddo
      do i=1,ny
         ygr(i) = ygmin + ((ygmax-ygmin)/(ny_tmp-1))*(i-1)
*         print *,'y(',i,')=',y(i)
      enddo
      do i=1,nz
         zgr(i) = zgmin + ((zgmax-zgmin)/(nz_tmp-1))*(i-1)
*         print *,'z(',i,')=',z(i)
      enddo
      
      return
      end

