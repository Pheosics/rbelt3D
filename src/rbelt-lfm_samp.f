
************************************************************************

      subroutine load_field
     &(basename,firstfilenum,firststep,laststep,gridstep)

c puts LFM fields into bxd, byd, bzd 
c exd, eyd & ezd. also puts corresponding times into tgr array.
c
c This code expects sequentially numbered LFM HDF files: 
c basename-<filenum>.hdf with one time step in each file 
c
c The fields are put into arrays bxd(nx,ny,nz,nt), byd(nx,,ny,... etc.
c defined in grid.inc. The fields are read into bxd, byd, etc. starting 
c with location e.g. bxd(nx,ny,nz,gridstep) up to 
c bxd(nx,ny,nz,(gridstep+laststep-firststep+1)).

      implicit none
      include 'rbelt-grid.inc'
      include 'rbelt-const.inc'
      include 'rbelt-geopak_08.inc'
      include 'rbelt-ut.inc'
      include 'rbelt-lfm.inc'

c     input vars.
      character*(*) basename
      integer firstfilenum,firststep,laststep,gridstep,lnblnk
c     parameters
      real inbound,outbound,rmin_interp,rmax_interp,b0_nt

      logical interp
      common /first/interp
      data interp /.true./

c *************************************************************************
c *************************************************************************

      parameter (inbound=2.2)  !inner boundary
      parameter (outbound=30.0)  !outer boundary
*      parameter (inbound=1.0-1.733*.2)  !inner boundary
*      parameter (outbound=16.0+1.733*.2)  !outer boundary
      parameter (rmin_interp=1.21)  !inner boundary of interpolation region
      parameter (rmax_interp=2.21)  !outer boundary of interpolation region
      parameter (b0_nt=b0/ffactor/ntg)

c *************************************************************************
c *************************************************************************

c     define hdftake in and out variables:
c     See note at top of hdftake.f file if you wish to 
c     see how the mhd variables are put into the array allcart(14,nx,ny,nz).
c     (in this version, first declared above)
      integer ierr,hdftake
      character*80 filename,fieldfile
      integer mstep,npts
      real time,xlim(6),allcart(14,nx,ny,nz),points(3,nx,ny,nz)
*      common /big/allcart,points

c     additional lfm_sample variables:
      character*8 string_out
      integer i,j,k,l,string_begin,string_end,filenum
      real bx_ex,by_ex,bz_ex,bx_in,by_in,bz_in,exwght,inwght
      real r, r2, r5
*      integer julday2,julday,sec_frm_ref

c     geopack & Tsyganenko routines
      real XGSM,YGSM,ZGSM
      real BXGSM,BYGSM,BZGSM
      real tilt
      integer mode
      integer year,doy,hour,min,sec

      print *
      print *,'*** in subroutine load_field (lfm_samp version) ***'

c     these get passed to hdftake
      xlim(1)=xgmin
      xlim(2)=xgmax
      xlim(3)=ygmin
      xlim(4)=ygmax 
      xlim(5)=zgmin
      xlim(6)=zgmax
      mstep=0
      npts=nx*ny*nz

c     set up points
      do i=1,nx
         do j=1,ny
            do k=1,nz
               points(1,i,j,k) = xgr(i)
               points(2,i,j,k) = ygr(j)
               points(3,i,j,k) = zgr(k)
            enddo
         enddo
      enddo

      filenum = firstfilenum
      do l=gridstep,(gridstep+laststep-firststep)

c        set input filename here
         call int2str(string_out,filenum,string_begin,string_end)
         filename=basename(1:lnblnk(basename))//
     &   string_out(1:string_end)//'.hdf'
         filenum = filenum + 1

         print *
         print *,'reading from ',filename
         print *,' in hdftake *************** '
c        get mhd values at points using hdftake
         ierr=hdftake(filename,xlim,interp,allcart,points,npts,time)
         print *,' exit hdftake ************* '
         print *,'time=',time
         if (ierr .ne. 0) then 
            write (6,*) 'error in hdftake'
            stop
         endif
c        UT = year0,doy0,hour0,min0,sec0 + t (in seconds)
         tgr(l) = time + time_add
         year=year0
         doy=doy0
         hour=hour0
         min=min0
         sec=sec0+tgr(l)
         call reduce_time2(year,doy,hour,min,sec)
         print *
         print *,'calling RECALC'
         print *,'year,doy,hour,min,sec=',
     &   year,doy,hour,min,sec
         call RECALC_08 (year,doy,hour,min,sec,-400.,0.,0.)
*         print *,'done with RECALC: SPS,CPS=',SPS,CPS
         print *,'done with RECALC: tilt=',atan2(SPS,CPS)*raddeg

         do k=1,nz
*           print *,'k=',k
           do j=1,ny
             do i=1,nx
               r2=xgr(i)*xgr(i)+ygr(j)*ygr(j)+zgr(k)*zgr(k)
               r=sqrt(r2)
               r5=r2*r2*r
               if ((r.le.outbound).and.(r.ge.inbound)) then
                 exwght=(r-rmin_interp)/(rmax_interp-rmin_interp)
                 inwght=1.0-exwght
                 if (r.gt.rmin_interp) then
                   bx_ex=allcart(6,i,j,k)
                   by_ex=allcart(7,i,j,k)
                   bz_ex=allcart(8,i,j,k)

*                   if (i.eq.(nx).and.j.eq.((ny+1)/2).and.
*     &             k.eq.((nz+1)/2))then
*                      print *,'i,j,k,allcart(8,i,j,k)=',i,j,k,bz_ex
**                      stop
*                   endif

                   exd(i,j,k,l) = -allcart(9,i,j,k)
                   eyd(i,j,k,l) = -allcart(10,i,j,k)
                   ezd(i,j,k,l) = -allcart(11,i,j,k)
                 else
                   inwght=1.0
                   exwght=0.0
                   bx_ex=0.0
                   by_ex=0.0
                   bz_ex=0.0
                   exd(i,j,k,l) = 0.0
                   eyd(i,j,k,l) = 0.0
                   ezd(i,j,k,l) = 0.0
                 endif
                 if (r.lt.rmax_interp) then
*                  print *,'calling SMGSM'
                   mode=1
                   CALL SMGSW_08
     &             (xgr(i),ygr(j),zgr(k),XGSM,YGSM,ZGSM,mode)
                   CALL IGRF_GSW_08(XGSM,YGSM,ZGSM,BXGSM,BYGSM,BZGSM)
                   mode=-1
                   CALL SMGSW_08
     &             (bx_in,by_in,bz_in,BXGSM,BYGSM,BZGSM,mode)
                 else
                   inwght=0.0
                   exwght=1.0
                   bx_in=0.0
                   by_in=0.0
                   bz_in=0.0
                 endif
c                dipole field gets removed
                 bxd(i,j,k,l)=exwght*bx_ex + inwght*bx_in +
     &           3.*b0_nt*xgr(i)*zgr(k)/r5
                 byd(i,j,k,l)=exwght*by_ex + inwght*by_in +
     &           3.*b0_nt*ygr(j)*zgr(k)/r5
                 bzd(i,j,k,l)=exwght*bz_ex + inwght*bz_in +
     &           3.*b0_nt*zgr(k)*zgr(k)/r5 - b0_nt/r2/r

*                   if (i.eq.(nx).and.j.eq.((ny+1)/2).and.
*     &             k.eq.((nz+1)/2))then
*                      print *,'i,j,k,l,bzd=',i,j,k,l,bzd(i,j,k,l)
**                      stop
*                   endif

               else
                 bxd(i,j,k,l) = 0.0
                 byd(i,j,k,l) = 0.0
                 bzd(i,j,k,l) = 0.0
                 exd(i,j,k,l) = 0.0
                 eyd(i,j,k,l) = 0.0
                 ezd(i,j,k,l) = 0.0
               endif
             enddo
           enddo
         enddo

*c        uncomment to sanity check fields along x axis
*         print *
*         print *,'purturbation fields'
*         do i=1,nx
**           print *,'x,bz=',x(i),
**           bz(i,(ny+1)/2,(nz+1)/2,n)/bscal
*           print *,xgr(i),bzd(i,(ny+1)/2,(nz+1)/2,l)
*         enddo
*         print *
*         print *,'total field'
*         do i=1,int(nx/2)
**           print *,'x,bz=',x(i),
**           bz(i,(ny+1)/2,(nz+1)/2,n)/bscal
*           print *,xgr(i),
*     &     (bzd(i,(ny+1)/2,(nz+1)/2,l)+b0_nt/abs(xgr(i)**3))
*         enddo
*         do i=int(nx/2)+2,nx
**           print *,'x,bz=',x(i),
**           bz(i,(ny+1)/2,(nz+1)/2,n)/bscal
*           print *,xgr(i),
*     &     (bzd(i,(ny+1)/2,(nz+1)/2,l)+b0_nt/abs(xgr(i)**3))
*         enddo
**         stop	 

      enddo

      return
      end

************************************************************************

      subroutine load_grid(filename)
      implicit none
      include 'rbelt-grid.inc'
      character*80 filename
      integer i,nx_tmp,ny_tmp,nz_tmp
      
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

************************************************************************

      subroutine load_ut(basename,firstfilenum,lastfilenum,filenumstep)
c     read ut
      implicit none
      include 'rbelt-ut.inc'
      include 'rbelt-lfm.inc'
      character*(*) basename
      integer firstfilenum,lastfilenum,filenumstep

      year0=lfmyear
      doy0=lfmdoy
      hour0=lfmhour
      min0=lfmmin
      sec0=lfmsec

      return
      end
