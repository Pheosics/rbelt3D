
      program lfm_sample
	
c ************************************************************************* 
c
c Reads a user specified number of LFM HDF output files & puts E&B fields on 
c a cartisian grid. The output is a user specified number of HDF files for use 
c with the rbelt code. Note that the fields are scaled for use with the rbelt
c code in here. Also, the dipole term is removed from the magnetic field.
c The HDF output files are set up as follows:
c
c file *-grid.hdf has the sds sets
c 1:(nx)x
c 2:(ny)y
c 3:(nz)z
c
c file *field-*.hdf has the sds sets
c 1:(nt)t
c 2:(nt x nx x ny x nz)bx 
c 3:(nt x nx x ny x nz)by 
c 4:(nt x nx x ny x nz)bz 
c 5:(nt x nx x ny x nz)ex 
c 6:(nt x nx x ny x nz)ey 
c 7:(nt x nx x ny x nz)ez 
c	
c B.Kress
c
c *************************************************************************

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'param.inc'
      include 'intco1d.inc'

c     declare parameter types
*      implicit none
      character* (*) basename,tsy_input
      integer numx,numy,numz,nt,first_ifile,last_ifile
      real xmin,xmax,ymin,ymax,zmin,zmax
      integer first_step,last_step,num_ifiles,tsteps,dsteps
      real inbound,outbound,rmax_interp,rmin_interp,dipole
      real time_mult,time_add,charge_sign,mass
      integer iyear0,idoy0,ihour0,imin0,isec0,max_table
      logical scale,tsy

c *************************************************************************
c *************************************************************************
c BELOW ARE THE PARAMETERS YOU NEED TO SET:
c
c IF YOU ARE READING LFM FILES FROM A RUN WHERE AN LFM GRID DIFFERENT 
c FROM (i,j,k) = (53,24,32) WAS USED, THEN YOU WILL ALSO NEED TO CHANGE
c THE VALUES IN lparam.inc AND param.inc, OR LINK THESE TO THE PROPER FILES.
c YOU ALSO NEED TO SET PARAMETER (num1d=numx*numy*numz) IN intco1d.inc
c (SEE BELOW).

c DEFINE RBELT GRID PARAMETERS:
c Use parameters below to define cartisian grid for use with rbelt code.
c For the LFM grid with (i,j,k) = (53,24,32) we have:
c max z =   126.5580    
c max y =   126.0819    
c max x =   30.03723    
c min z =  -126.5580    
c min y =  -126.0819    
c min x =  -345.2204   
c so we should stay within this outer boundary.
c also, note that the inner boundary of the LFM grid is at ~2.2 Re,
c and all points inside inbound (defined below) + 0.001 RE 
c will return zero fields.
c Do not forget to set parameter (num1d=numx*numy*numz) in intco1d.inc!
c WARNING, size of output file = numx*numy*numz*nt*6*4 & you will need
c approx. numx*numy*numz*nt*6*4*2 bytes of RAM to run this code.
c nt is #of time steps  output files, except 1st output file 
c has nstep=nt+2 time steps.

c      for wind sample - num1d=1
*      parameter (numx=1, numy=1, numz=1, nt=40) !nt is #of time steps in output file
*      parameter (xmin=29., xmax=29.)
*      parameter (ymin=0., ymax=0.)
*      parameter (zmin=0., zmax=0.)

c      for time series along x - num1d=135
      parameter (numx=135, numy=1, numz=1, nt=6) 
      parameter (xmin=2.2, xmax=29.)
      parameter (ymin=0., ymax=0.)
      parameter (zmin=0., zmax=0.)

c      num1d=14641
*      parameter (numx=1, numy=121, numz=121, nt=6) 
*      parameter (xmin=3., xmax=3.)
*      parameter (ymin=-12., ymax=12.)
*      parameter (zmin=-12., zmax=12.)

c      num1d=14641
*      parameter (numx=121, numy=1, numz=121, nt=1) 
*      parameter (xmin=-12., xmax=12.)
*      parameter (ymin=0., ymax=0.)
*      parameter (zmin=-12., zmax=12.)

c      num1d=46431
*      parameter (numx=231, numy=201, numz=1, nt=2) 
*      parameter (xmin=-100., xmax=15.)
*      parameter (ymin=-50., ymax=50.)
*      parameter (zmin=0., zmax=0.)

c CHOOSE INPUT FILES (LFM DUMP FILES) AND THE STEPS TO READ
c set the lfm dump file names here e.g. if the first filename is Nov01-25.hdf
c then basename='Nov01-' and first_ifile=25

      parameter (basename='Nov01-')  !prefix for lfm dump file name
      parameter (first_ifile=13)  !first input file number
      parameter (last_ifile=15)  !last input file number
      parameter (first_step=120000)  !first lfm time step in first_ifile
      parameter (start_step=123200)  !first lfm time step to read
*      parameter (last_step=123700)
*      parameter (start_step=123720)  !first lfm time step to read
      parameter (last_step=129980)  !last lfm time step to read
      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this

c ADDITIONAL INFORMATION ABOUT THE LFM DUMP FILES
c This info. must be uniform for all lfm dump files to read. tsteps is the 
c number of lfm time steps per lfm dump steps (or per dump steps to read)
c first_step, start_step, & last_step must be evenly  divisible by tsteps

      parameter (tsteps=20) 

c the total number of lfm dump steps in (or to read from) each lfm dump file is

      parameter (dsteps=250) 

c e.g., if the LFM run has NSTOP-NSTART=10000 and NDUMP=250, which is typical,
c then tsteps=250 and dsteps=40. Note that there are actually 41 steps in each 
c lfm dump file but the first step is the same as the last step in the previous
c file, except for the initial step in file #1 which is dump step #1. i.e. if 
c the last dump step in file #1 is dump step 41, then dump step 42 is the 
c second step in file #2.
c Nov01-25 240250 250000 250 steps 962-1001
c Nov01-26 250250 260000 250 steps 1002-1041
c Nov01-27 260250 270000 250 steps 1042-1081
c Nov01-28 270250 280000 250 steps 1082-1121
c etc..
c If NDUMP is uniform for entire LFM run
c LFM dump step is (time step)/tsteps + 1
c LFM file is (last dump step in file - 1)/dsteps
      
c SUBTRACT OFF DIPOLE FIELD

c     subtracted dipole must be the same as b0 to be used in rbelt code!!
c     use bzero from LFM dipole field 27500nt,31000nt or whatever.
      parameter (dipole=30000.0)

c SET TIME CORRECTION FACTOR AND ADDITIVE CONSTANT TO GET TIME RIGHT

      parameter (iyear0=2001) !YEAR NUMBER (FOUR DIGITS)
      parameter (idoy0=328) !DAY OF YEAR (DAY 1 = JAN 1)
      parameter (ihour0=0) !HOUR OF DAY (00 TO 23)
      parameter (imin0=0) !MINUTE OF HOUR (00 TO 59)
      parameter (isec0=0) !SECONDS OF MINUTE (00 TO 59)
      parameter (time_add=-2520.0) ! time_add = 0.0 gives sec from lfm run start
      parameter (time_mult=1.0) ! set=1.0 to get seconds

c TO SCALE FIELDS HERE FOR USE WITH SOME VERSIONS OF THE RBELT CODE,
c set scale=.true. and specify particle charge and mass.

      parameter (scale=.false.) !scale fields
      parameter (charge_sign=-1)  !reverse charge_sign to run backwards in time
*      parameter (mass=1.6726485e-24)  !proton mass
      parameter (mass=9.109534e-28)  !electron mass

c INTERIOR FIELDS (uses IGRF)
c the output will give zero fields for points inside inbound
c and outside outbound. 

      parameter (inbound=1.0)  !inner boundary
      parameter (outbound=30.0)  !outer boundary
      parameter (rmin_interp=2.3)  !inner boundary of interpolation region
      parameter (rmax_interp=3.3)  !outer boundary of interpolation region
      parameter (tsy=.false.) !add Tsyganenko to interior IGRF field model
      parameter (tsy_input='T02_input.dat') ! Tsyganenko/geopack input file name
      parameter (max_table=5000)

c *************************************************************************
c *************************************************************************

c     define output variables:
      real x(numx),y(numy),z(numz),t(nt+2)
      real bx(numx,numy,numz,nt+2),by(numx,numy,numz,nt+2)
      real bz(numx,numy,numz,nt+2),ex(numx,numy,numz,nt+2)
      real ey(numx,numy,numz,nt+2),ez(numx,numy,numz,nt+2)
      real vx(numx,numy,numz,nt+2),vy(numx,numy,numz,nt+2)
      real vz(numx,numy,numz,nt+2),density(numx,numy,numz,nt+2)

c     additional lfm_sample variables:
      real bx_ex,by_ex,bz_ex,bx_in,by_in,bz_in,exwght,inwght
      integer i,j,k,l,m
      integer lstep(num_ifiles)
      integer nxtstep,ofile_num,ifile_num,index,halt
      real nx_tmp,ny_tmp,nz_tmp,b,time,t0 !initial time
      real Re, q, c, m0, b0, r, r2, r5 
      real tscal,bscal,escal,ffactor
      character*80 ofile
      character*8 string_out
      integer string_begin,string_end
      logical makegrid

c     define hdftake in and out variables:
c     See note at top of hdftake.f file if you wish to 
c     see how the mhd variables are put into the array allcart(14,nx,ny,nz).
      integer hdftake
      logical first,interp
      character*80 ifile
      integer mstep,npts
      real xlim(6),allcart(14,numx,numy,numz),points(3,numx,numy,numz)
      common /big/bx,by,bz,ex,ey,ez,allcart,points

c     HDF stuff:
      integer*4 file_id, sds_id, rank, status
      integer dim_sizes4d(4), start4d(4), edge4d(4),stride4d(4)
      integer dim_sizes1d(1), start1d(1), edge1d(1),stride1d(1)
      integer n_datasets, n_file_attrs

c     geopack & Tsyganenko routines
      integer num_table,year(max_table),doy(max_table)
      real hour(max_table),pram(max_table),dst(max_table),
     &byimf(max_table),bzimf(max_table),ps(max_table),
     &G2(max_table),G3(max_table)
      real delta,sec1,sec2
      real XGSM,YGSM,ZGSM
      real HXGSM,HYGSM,HZGSM,BXGSM,BYGSM,BZGSM
      REAL*8 PARMOD(10),tilt,X_dbl,Y_dbl,Z_dbl,BX_dbl,BY_dbl,BZ_dbl
      real PDYN,DST,BYIMF,BZIMF
      integer IYEAR,IDAY,IHOUR,MIN,ISEC
      real ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
     * CGST,SGST,BA(6),G(105),H(105),REC(105)
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
     * CGST,SGST,BA
      COMMON /GEOPACK2/ G,H,REC

c     Define constants
      Re = 6378.137e5
      q = 4.8032424e-10
      c = 2.99792458e10
      m0 = mass
      ffactor=charge_sign*q*Re/(c*c*m0)
      b0 = dipole

      print *
      print *,'running program lfm_sample'
      print *,'numx,numy,numz=',numx,numy,numz
      print *,'xmin,xmax=',xmin,xmax
      print *,'ymin,ymax=',ymin,ymax
      print *,'zmin,zmax=',zmin,zmax
      print *,'basename=',basename
      print *,'first_ifile,last_ifile=',first_ifile,last_ifile
      print *,'first_step,start_step,last_step=',
     &first_step,start_step,last_step
      print *,'tsteps,dsteps=',tsteps,dsteps
      print *,'dipole=',dipole
      print *,'iyear0,idoy0,ihour0,imin0,isec0=',
     &iyear0,idoy0,ihour0,imin0,isec0
      print *,'time_add,time_mult=',time_add,time_mult
      print *,'scale,charge_sign,mass=',scale,charge_sign,mass
      print *,'inbound,outbound=',inbound,outbound
      print *,'rmin_interp,rmax_interp=',rmin_interp,rmax_interp
      print *,'tsy,tsy_input=',tsy,tsy_input
      print *,'max_table=',max_table
      if (scale .eqv. .true.) then
         print *,'fileds put into dimensionless unit system'
         bscal=ffactor*1.e-5
         escal=ffactor*1.e6/c
         tscal = c/Re
      else
         print *,'no scaling, fileds in nT & V/m'
         bscal =1.0
         escal =1.0
         tscal =1.0
      endif
      print *,'bscal=',bscal
      print *,'escal=',escal
      print *,'tscal(tfactor) =',tscal 
      print *,'b0=',b0

c     these get passed to hdftake
      first=.true.
      interp=.true.
      xlim(1)=xmin
      xlim(2)=xmax
      xlim(3)=ymin
      xlim(4)=ymax 
      xlim(5)=zmin
      xlim(6)=zmax
      mstep=0
      npts=numx*numy*numz

c     set up points
      nx_tmp = numx
      ny_tmp = numy
      nz_tmp = numz
      if (numx.eq.1) nx_tmp=2
      if (numy.eq.1) ny_tmp=2
      if (numz.eq.1) nz_tmp=2
      do i=1,numx
         do j=1,numy
            do k=1,numz
               points(1,i,j,k) = xmin + ((xmax-xmin)/(nx_tmp-1))*(i-1)
               points(2,i,j,k) = ymin + ((ymax-ymin)/(ny_tmp-1))*(j-1)
               points(3,i,j,k) = zmin + ((zmax-zmin)/(nz_tmp-1))*(k-1)
            enddo
         enddo
      enddo
      do i=1,numx
         x(i) = xmin + ((xmax-xmin)/(nx_tmp-1))*(i-1)
*         print *,'x(',i,')=',x(i)
      enddo
      do i=1,numy
         y(i) = ymin + ((ymax-ymin)/(ny_tmp-1))*(i-1)
*         print *,'y(',i,')=',y(i)
      enddo
      do i=1,numz
         z(i) = zmin + ((zmax-zmin)/(nz_tmp-1))*(i-1)
*         print *,'z(',i,')=',z(i)
      enddo

c     set up lstep array (last dump step to read in each lfm dump file)
      lstep(1) = dsteps+first_step/tsteps
*      print *,'lstep(1)=',lstep(1)
      do l=2,num_ifiles
         lstep(l) = lstep(l-1) + dsteps
*         print *,'l,lstep(l)=',l,lstep(l)
      enddo
      lstep(num_ifiles)=last_step/tsteps
*      print *,'num_ifiles,lstep(num_ifiles)=',
*     &num_ifiles,lstep(num_ifiles)

c     initialize loop
      nxtstep = start_step/tsteps      ! next step
      ifile_num = first_ifile   ! first input file number
      index = 0                 ! output file write index
      nstep = nt + 2            ! # of steps in output file
      makegrid =.true.		! put grid info into a file
      ofile_num = start_step	! output file number

      print *
      print *,'total read/write steps =',(last_step-start_step)/tsteps+1

c     set input filename here
      call int2str(string_out,ifile_num,string_begin,string_end)
      ifile=basename//string_out(1:string_end)//'.hdf'
      
c     open output file
      open (14,file='lfm_plot.dat',status='unknown')
      write(14,96) int((last_step-start_step)/tsteps+1),numx,numy,numz
96    format (4i8)
      write(14,97) xmin,xmax,ymin,ymax,zmin,zmax
97    format (6f10.2)

c     open geopack and Tsyganenko input file
      open (UNIT=12,file='T02_input.dat',STATUS='OLD')
      read (12,*) num_table
      print *,'opened T02_input.dat: num_table=',num_table
      do i=1,num_table
        read (12,50) year(i),doy(i),hour(i),pram(i),dst(i),byimf(i),
     &  bzimf(i),ps(i),G2(i),G3(i)
50      format (2i4,8f9.3)
      enddo
      print *,'i=',i
      close(12)
      
      if (int(year(1)).ne.iyear0) then
        print *,'year(1).ne.iyear0',year(1),iyear0
        stop
      endif

c     read and write loop starts here
100   if (nxtstep*tsteps .le. last_step) then
c       read from hdf input file until we reach end or we have enough
c       steps to fill output file
        halt=min0(nstep-index,lstep(ifile_num-first_ifile+1)-nxtstep+1)
c       but stop when we have read in the last_step
        halt=min0(halt,(last_step/tsteps-nxtstep+1))
*        print *,'nstep-index=',nstep-index
*        print *,'lstep(ifile_num-first_ifile+1)-nxtstep+1=',
*     &  lstep(ifile_num-first_ifile+1)-nxtstep+1
*        print *,'(last_step-nxtstep*tsteps+1)=',
*     &  (last_step-nxtstep*tsteps+1)
*        print *,'halt=',halt
        do l=1,halt
          index = index + 1
          mstep = nxtstep*tsteps
          print *
          print *,'reading time step ',mstep,' from ',ifile
          print *
          print *,' in hdftake *************** '
c         get mhd values at points using hdftake
          ierr = hdftake(first,ifile,mstep,xlim,
     &         interp, allcart, points, npts, time)
          print *,' exit hdftake ************* '
*          print *,'time=',time
          first=.false.
          if (ierr .ne. 0) then 
             write (6,*) 'error in hdftake'
             stop
          else
            t(index)=(time*time_mult + time_add)*tscal
            m=1
            sec2=(doy(m)-idoy0)*86400.0 + hour(m)*3600.0
            if (sec2.gt.t(index)) then
              print *,'sec2.gt.t(index)',sec2,t(index)
              stop
            endif
110         m=m+1
              sec1=sec2
              sec2=(doy(m)-idoy0)*86400.0+hour(m)*3600.0
              if (sec2.lt.t(index)) goto 110
c
            delta=(t(index)-sec1)/(sec2-sec1)
            PARMOD(1)=pram(m-1)+delta*(pram(m)-pram(m-1))
            PARMOD(2)=dst(m-1)+delta*(dst(m)-dst(m-1))
            PARMOD(3)=byimf(m-1)+delta*(byimf(m)-byimf(m-1))
            PARMOD(4)=bzimf(m-1)+delta*(bzimf(m)-bzimf(m-1))
            PARMOD(5)=G2(m-1)+delta*(G2(m)-G2(m-1))
            PARMOD(6)=G3(m-1)+delta*(G3(m)-G3(m-1))
	    tilt=ps(m-1)+delta*(ps(m)-ps(m-1))
            iday=idoy0+int(t(index)/86400.0)
            ihour=int(t(index)/3600.0)-(iday-idoy0)*24
            min=int(t(index)/60.0)-(iday-idoy0)*1440-ihour*60
            isec=t(index)-(iday-idoy0)*86400-ihour*3600-min*60
	    
	    print *,'m,t(sec),t(min),t(hr)=',
     &      m,t(index),t(index)/60.0,t(index)/3600.0
	    print *,'pram,dst,byimf,bzimf,G2,G3=',PARMOD(1),
     &      PARMOD(2),PARMOD(3),PARMOD(4),PARMOD(5),PARMOD(6)
     	    print *,'IYEAR0,IDAY,IHOUR,MIN,ISEC=',
     &      IYEAR0,IDAY,IHOUR,MIN,ISEC
	    
            call RECALC (IYEAR0,IDAY,IHOUR,MIN,ISEC)
	    print *,'done with RECALC:tilt,PS=',tilt,atan2(SPS,CPS)
	    
            do k=1,numz
              do j=1,numy
                do i=1,numx
                  r2=x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
                  r=sqrt(r2)
                  r5=r2*r2*r
                  if ((r.le.outbound).and.(r.ge.inbound)) then
                    exwght=(r-rmin_interp)/(rmax_interp-rmin_interp)
                    inwght=1.0-exwght
                    if (r.gt.rmin_interp) then
                      bx_ex=allcart(6,i,j,k)
                      by_ex=allcart(7,i,j,k)
                      bz_ex=allcart(8,i,j,k)
                      ex(i,j,k,index) = allcart(9,i,j,k)*escal
                      ey(i,j,k,index) = allcart(10,i,j,k)*escal
                      ez(i,j,k,index) = allcart(11,i,j,k)*escal
                      density(i,j,k,index) = allcart(1,i,j,k)
                      vx(i,j,k,index) = allcart(2,i,j,k)
                      vy(i,j,k,index) = allcart(3,i,j,k)
                      vz(i,j,k,index) = allcart(4,i,j,k)
                    else
                      inwght=1.0
                      bx_ex=0.0
                      by_ex=0.0
                      bz_ex=0.0
                      ex(i,j,k,index) = 0.0
                      ey(i,j,k,index) = 0.0
                      ez(i,j,k,index) = 0.0
                    endif
                    if (r.lt.rmax_interp) then
*		      print *,'calling SMGSM'
                      CALL SMGSM(x(i),y(j),z(k),XGSM,YGSM,ZGSM,1)
                      if (tsy.eqv..false.) then
*                        print *
*                        print *,'x(i),y(j),z(k)=',x(i),y(j),z(k)
*                        print *,'XGSM,YGSM,ZGSM=',XGSM,YGSM,ZGSM
*                        print *,'r=',sqrt(XGSM**2+YGSM**2+ZGSM**2)
                        CALL IGRF_GSM(XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)
*                        print *,'HXGSM,HYGSM,HZGSM=',HXGSM,HYGSM,HZGSM
*                        print *,'B=',sqrt(HXGSM**2+HYGSM**2+HZGSM**2)
                      else
*		        print *,'calling IGRF_GSM'
                        CALL IGRF_GSM(XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)
*		        print *,'calling T02'
                        x_dbl=XGSM
                        y_dbl=YGSM
                        z_dbl=ZGSM
                        CALL T02(0,PARMOD,tilt,x_dbl,y_dbl,z_dbl,
     &                  BX_dbl,BY_dbl,BZ_dbl)
*     		        print *,'exit T02'
                      endif
*		      print *,'calling SMGSM'
                      BXGSM=HXGSM+BX_dbl
                      BYGSM=HYGSM+BY_dbl
                      BZGSM=HZGSM+BZ_dbl
                      CALL SMGSM(bx_in,by_in,bz_in,BXGSM,BYGSM,BZGSM,-1)
                    else
                      exwght=1.0
                      bx_in=0.0
                      by_in=0.0
                      bz_in=0.0
                    endif
c                   dipole field gets removed & scalling done here
                    bx(i,j,k,index)=(exwght*bx_ex + inwght*bx_in)*bscal
                    by(i,j,k,index)=(exwght*by_ex + inwght*by_in)*bscal
                    bz(i,j,k,index)=(exwght*bz_ex + inwght*bz_in)*bscal
                  else
                    bx(i,j,k,index) = 0.0
                    by(i,j,k,index) = 0.0
                    bz(i,j,k,index) = 0.0
                    ex(i,j,k,index) = 0.0
                    ey(i,j,k,index) = 0.0
                    ez(i,j,k,index) = 0.0
                    density(i,j,k,index) = 0.
                    vx(i,j,k,index) = 0.
                    vy(i,j,k,index) = 0
                    vz(i,j,k,index) = 0
                  endif
c                  uncomment to get a pure dipole on grid (must also 
c                  remove dipole from rbelt interp routine by setting b0=0)
*                  bx(i,j,k,index) = (3.*b0*x(i)*z(k)/r5)*bscal
*                  by(i,j,k,index) = (3.*b0*y(j)*z(k)/r5)*bscal
*                  bz(i,j,k,index) = (3.*b0*z(k)*z(k)/r5
*     &                                          - b0/r2/r)*bscal
*                  ex(i,j,k,index) = 0.
*                  ey(i,j,k,index) = 0.
*                  ez(i,j,k,index) = 0.
                enddo
              enddo
            enddo
          endif


c         uncomment to sanity check perturbation fields along x axis
*          print *
*          do i=1,numx
**            print *,'x,bz=',x(i),
**     &      bz(i,(numy+1)/2,(numz+1)/2,index)/bscal
*            print *,x(i),bz(i,(numy+1)/2,(numz+1)/2,index)/bscal
*          enddo
*          print *
*          do i=1,numx
**            print *,'x,bz=',x(i),
**     &      bz(i,(numy+1)/2,(numz+1)/2,index)/bscal
*            print *,x(i),
*     &      (bz(i,(numy+1)/2,(numz+1)/2,index)+b0/abs(x(i)**3))
*          enddo

          write(14,98)nxtstep*tsteps,t(index)
98	  format (i8,1f12.4)
          do k=1,numz
             do j=1,numy
                do i=1,numx
                   write(14,99)x(i),y(j),z(k),
     &             density(i,j,k,index),vx(i,j,k,index),
     &             vy(i,j,k,index),vz(i,j,k,index),bx(i,j,k,index),
     &             by(i,j,k,index),bz(i,j,k,index),ex(i,j,k,index),
     &             ey(i,j,k,index),ez(i,j,k,index)
99	           format (13f12.4)
                enddo
             enddo
          enddo

          nxtstep = nxtstep + 1
        enddo    
c       if we have reached the end of the ifile then increment ifile_num
        if (nxtstep .gt. lstep(ifile_num-first_ifile+1)) then
          first = .true.
          ifile_num = ifile_num + 1
*          print *,'increment ifile_num to',ifile_num
          call int2str(string_out,ifile_num,string_begin,string_end)
          ifile=basename//string_out(1:string_end)//'.hdf'
        endif
c       if we have read in nt hdf steps or we are done, then write out to HDF 
        if ((index.eq.nstep).or.(ifile_num.gt.last_ifile)) then
c         send x,y,z,t,bx,by,bz,ex,ey,ez out to HDF file
          rank=1
          start1d(1)=0
          stride1d(1)=1
          if (makegrid.eqv..true.) then
            makegrid = .false.
            nstep = nstep - 2
          endif
          ofile_num = nxtstep*tsteps	    
          index = 0
        endif
        goto 100
      endif

      close(14)

      end
	

c ******************************************************************************
c ******************************************************************************

*	SUBROUTINE INT2STR(CAOUT,IVALUE,IBEGIN,IEND)
*
*	character*(*) caout
*	integer iValue, iBegin, iEnd
*	integer i,iStart,iLast
*	character*20 caDum
*
*	write(caDum,10) iValue
* 10     format(i8)
*
*	iStart = 1
*	do while( caDum(iStart:iStart) .eq. ' ' )
*	  iStart = iStart + 1
*	enddo
*	iLast = lnblnk(caDum)
*
*	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
*        iBegin = 1
*	iEnd = iLast - iStart + 1
*
*	return
*
*	end

	
c ******************************************************************************
c ******************************************************************************

	integer function lnblnk(str)

C
C Description:
C
C This function returns the location of the last non-blank
C character in a string. If all the characters are blank
C the function returns 0.
C
C Important note: SUN FORTRAN provides this function. SUNS version is probably
C more efficient than this one so use SUNs if you are using the SUN
C compiler. This routine is provided only for those FORTRANs that 
C do not have this function.
C
C User interface:
C
C Name         in/out type          Structure  Meaning
C lnblnk       out    integer       scalar     The position of the last       
C                                              non-blank character in the    
C                                              input string. (Zero is returned
C                                              if the entire string is blank)
C str          in     character(*)  scalar     The string in which the last
C                                              non-blank character is to be
C                                              found.
C
C Errors:
C
C There should be no errors returned from this routine.
C
      implicit none
      integer limit,i
      logical done
      character*(*) str
      limit=1
      i=len(str)
      done=.false.
      do while (.not. done)
         if (str(i:i) .eq. ' ') then
            i=i-1
            if (i .lt. limit) then
               done=.true.
               lnblnk=0
            end if
         else
            done=.true.
            lnblnk=i
         endif
      enddo
      return
      end

************************************************************************

      FUNCTION julday(mm,id,iyyy)
      INTEGER julday,id,iyyy,mm,IGREG
      PARAMETER (IGREG=15+31*(10+12*1582)) 
c Gregorian Calendar adopted Oct. 15, 1582.
c In this routine julday returns the Julian Day Number that begins at noon of the calendar
c date specifed by month mm,day id, and year iyyy, all integer variables. Positive year
c signifes A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
      INTEGER ja,jm,jy
      jy=iyyy
      if (jy.eq.0) pause 'julday: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
         jm=mm+1
      else
         jy=jy-1
         jm=mm+13
      endif
      julday=365*jy+int(0.25d0*jy+2000.d0)+int(30.6001d0*jm)+id+1718995
      if (id+31*(mm+12*iyyy).ge.IGREG) then !Test whether to change to Gregorian Calendar.
         ja=int(0.01d0*jy)
         julday=julday+2-ja+int(0.25d0*ja)
      endif
      return
      END

************************************************************************

      SUBROUTINE caldat(julian,mm,id,iyyy)
      INTEGER id,iyyy,julian,mm,IGREG
      PARAMETER (IGREG=2299161)
c Inverse of the function julday given above. Here julian is input as a Julian Day Number,
c and the routine outputs mm,id, and iyyy as the month, day, and year on which the specied
c Julian Day started at noon.
      INTEGER ja,jalpha,jb,jc,jd,je
      if(julian.ge.IGREG)then !Cross-over to Gregorian Calendar produces this correction.
         jalpha=int(((julian-1867216)-0.25d0)/36524.25d0)
         ja=julian+1+jalpha-int(0.25d0*jalpha)
      else if(julian.lt.0)then 
c        Make day number positive by adding integer number of Julian centuries,
c        then subtract them off at the end.
         ja=julian+36525*(1-julian/36525)
      else
         ja=julian
      endif
      jb=ja+1524
      jc=int(6680.0d0+((jb-2439870)-122.1d0)/365.25d0)
      jd=365*jc+int(0.25d0*jc)
      je=int((jb-jd)/30.6001d0)
      id=jb-jd-int(30.6001d0*je)
      mm=je-1
      if(mm.gt.12)mm=mm-12
      iyyy=jc-4715
      if(mm.gt.2)iyyy=iyyy-1
      if(iyyy.le.0)iyyy=iyyy-1
      if(julian.lt.0)iyyy=iyyy-100*(1-julian/36525)
      return
      END

************************************************************************

      subroutine doy_mmdd (mode, year, doy, mm, dd, status)

C Converts day-of year to months and days within month (mode = 1)
C Converts months and days to day-of-year (mode = -1)
C Valid for years 1901 - 2099 inclusive

C Exit status
C             0  normal exit
C             1  year is outside valid range
C             2  illegal value for mode
C             3  this part of code not checked


      integer*4       mode
      integer*4       year
      integer*4       doy
      integer*4       mm
      integer*4       dd
      integer*4       status
      integer*4       i

      integer*4       month(12) /31,28,31,30,31,30,31,31,30,31,30,31/

C Initializations

      status = 0

      if ((year .le. 1900) .or. (year .ge. 2100)) then
        status = 1
        go to 90
      end if

      month(2) = 28
      if (mod(year,4) .eq. 0) month(2) = 29


C Enter loops for calculation

      if (mode .eq. 1) then
        i = 1
        dd = doy
        do while (dd .gt. month(i))
          dd = dd - month(i)
          i = i + 1
        end do
        mm = i

      else if (mode. eq. -1) then

        doy = 0
        i = 1
        do while (i .lt. mm)
          doy = doy + month(i)
          i = i + 1
        end do
        doy = doy + dd

      else
        status = 2
        go to 90

      end if

C Exit

   90 continue

      return
      end

************************************************************************

      subroutine dtcvt(yyyy,ddd,hh,mm,ss,string)

C Normalizes  date-time  yyyy, dd, hh, mm, ss  if any value is out of range (too large).
C Converts to  YYYY-MM-DDThh:mm:ss  in 19-character string.  Valid for non-leap-years.

      integer*4     yyyy, ddd    ! year and day-of-year
      integer*4     hh, mm, ss   ! hours, minutes, seconds
      integer*4     mon,dd       ! month and day-of-month
      integer*4     i            ! miscellaneous index
      character*19  string       ! YYYY-MM-DDThh:mm:ss (PDS date-time format)

      integer*4     dom(12)      ! days in months


C Initializations for non-leap year

      dom(01) = 31
      dom(02) = 28
      dom(03) = 31
      dom(04) = 30
      dom(05) = 31
      dom(06) = 30
      dom(07) = 31
      dom(08) = 31
      dom(09) = 30
      dom(10) = 31
      dom(11) = 30
      dom(12) = 31


C Begin normalization, if needed

 10   continue
      if (ss .gt. 59) then
        ss = ss - 60
        mm = mm + 1
        go to 10
      else if (ss .lt. 00) then
        ss = ss + 60
        mm = mm - 1
        go to 10
      end if

 20   continue
      if (mm .gt. 59) then
        mm = mm - 60
        hh = hh + 1
        go to 20
      else if (mm .lt. 00) then
        mm = mm + 60
        hh = hh - 1
        go to 20
      end if

 30   continue
      if (hh .gt. 23) then
        hh = hh - 24
        ddd = ddd + 1
        go to 30
      else if (hh .lt. 00) then
        hh = hh + 24
        ddd = ddd - 1
        go to 30
      end if

 40   continue

      if (ddd .gt. 365) then
        if (mod(yyyy,004) .eq. 0) then
          if ((mod(yyyy,100) .ne. 0) .or.
     *        (mod(yyyy,400) .eq. 0)) then
                if (ddd. gt. 366) then
                  ddd = ddd - 366
                  yyyy = yyyy + 1
                  go to 40
                end if
          end if
        else
          ddd = ddd - 365
          yyyy = yyyy + 1
          go to 40
        end if
      end if

      if (ddd .lt. 001) then
        if (mod(yyyy-1,004) .eq. 0) then
          if ((mod(yyyy-1,100) .ne. 0) .or.
     *        (mod(yyyy-1,400) .eq. 0)) then
                ddd = ddd + 366
                yyyy = yyyy - 1
                go to 40
          end if
        else
          ddd = ddd + 365
          yyyy = yyyy - 1
          go to 40
        end if
      end if


C Convert day-of-year to months and days

      if (mod(yyyy,004) .eq. 0) then
        if ((mod(yyyy,100) .ne. 0) .or.
     *      (mod(yyyy,400) .eq. 0)) then
              dom(02) = 29
        end if
      end if

      i = 1
      dd = ddd
      mon = 01
 60   continue
      if (dd .gt. dom(i)) then
        dd = dd - dom(i)
        mon = mon + 1
        i = i + 1
        go to 60
      end if

      write(string(01:19),'(i4.4,2("-",i2.2),"T",2(i2.2,":"),i2.2)')
     *    yyyy,mon,dd,hh,mm,ss


      return
      end

************************************************************************


