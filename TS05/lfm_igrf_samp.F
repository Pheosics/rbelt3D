c
c This is a wrapper code for hdftake, for puting LFM fields on a Cartisian grid 
c for use with the test-particle radiation belt codes. 
c
c WARNING, this code was originally used with LFM HDF files where there were 
c multiple time steps in each file and we knew the step numbers (mstep) in each file.
c Now, there is one LFM time step per HDF file, and the step number is not returned 
c by hdftake. The code has been hacked in various ways to accomadate, but really should 
c be completely overhauled. 
c

      program lfm_sample
	
c ************************************************************************* 
c
c Reads a user specified number of LFM HDF output files & puts E&B fields on 
c a cartisian grid. The output is a user specified number of HDF files for use 
c with the rbelt code. Also, the dipole term is removed from the magnetic field.
c The HDF output files are set up as follows:
c
c file *-grid.hdf has the SDS sets
c 1:(nx)x
c 2:(ny)y
c 3:(nz)z
c
c file *field-*.hdf has the SDS sets
c 1:(nt)t
c 2:(nt x nx x ny x nz)bx 
c 3:(nt x nx x ny x nz)by 
c 4:(nt x nx x ny x nz)bz 
c 5:(nt x nx x ny x nz)ex 
c 6:(nt x nx x ny x nz)ey 
c 7:(nt x nx x ny x nz)ez 
c	
c *************************************************************************

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'param.inc'
      include 'intco1d.inc'

c     declare parameter types
*      implicit none
      character* (*) basename
      integer numx,numy,numz,nt,nt_1,first_ifile,last_ifile
      real xmin,xmax,ymin,ymax,zmin,zmax
      integer first_step,start_step,last_step,num_ifiles,tsteps,dsteps
      real inbound,outbound,rmax_interp,rmin_interp,dipole
      real time_mult,time_add,charge_sign,mass
      integer iyear0,idoy0,ihour0,imin0,isec0,max_table
      logical lowpass

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






*c     for lws work
*      parameter (numx=81, numy=81, numz=41, nt=1) !parameter (num1d=269001)
**      parameter (numx=81, numy=81, numz=11, nt=6) !parameter (num1d=72171)
*      parameter (xmin=-8., xmax=8.)
*      parameter (ymin=-8., ymax=8.)
*      parameter (zmin=-4., zmax=4.)
**      parameter (zmin=-1., zmax=1.)
*      parameter (nt_1=1)  ! num. of time steps in 1st file


      parameter (numx=101, numy=101, numz=101, nt=1) !parameter
      parameter (xmin=-10., xmax=10.)
      parameter (ymin=-10., ymax=10.)
      parameter (zmin=-10., zmax=10.)
      parameter (nt_1=1)  ! num. of time steps in 1st file


c CHOOSE INPUT FILES (LFM DUMP FILES) AND THE STEPS TO READ
c set the lfm dump file names here e.g. if the first filename is Nov01-25.hdf
c then basename='Nov01-' and first_ifile=25


      parameter (basename='NAIRAS2-')  !prefix for lfm dump file name
      parameter (first_ifile=0)  !first input file number
      parameter (last_ifile=0)  !last input file number
      parameter (first_step=1)  !first lfm time step in first_ifile
      parameter (start_step=1)  !first lfm time step to read
      parameter (last_step=1)  !last lfm time step to read
      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this


*      parameter (basename='shockhdfs/nov01-')  !prefix for lfm dump file name
*      parameter (first_ifile=22)  !first input file number
*      parameter (last_ifile=22)  !last input file number
*      parameter (first_step=194000)  !first lfm time step in first_ifile
*      parameter (start_step=194000)  !first lfm time step to read
*      parameter (last_step=194000)  !last lfm time step to read
*      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this

*      parameter (basename='shockhdfs/nov01-')  !prefix for lfm dump file name
*      parameter (first_ifile=22)  !first input file number
*      parameter (last_ifile=30)  !last input file number
*      parameter (first_step=194000)  !first lfm time step in first_ifile
*      parameter (start_step=195800)  !first lfm time step to read
*      parameter (last_step=212000)  !last lfm time step to read
*      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this


*      parameter (basename='quiethdfs/jan03q-')  !prefix for lfm dump file name
*      parameter (first_ifile=1)  !first input file number
**      parameter (last_ifile=2)
*      parameter (last_ifile=661)  !last input file number
*      parameter (first_step=1)  !first lfm time step in first_ifile
*      parameter (start_step=1)  !first lfm time step to read
**      parameter (last_step=2)
*      parameter (last_step=661)  !last lfm time step to read
*      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this

c Warning, this code was originally used with LFM HDF files where there were 
c multiple time steps in each file and we knew the step numbers (mstep).
c Now, there is one time step per HDF file, and the step number is not returned.
c The code has been hacked in various ways, but really should be completely 
c overhauled. You may need to play around with the parameters above. Pay attention
c Pay attention to the code output.

c ADDITIONAL INFORMATION ABOUT THE LFM DUMP FILES
c This info. must be uniform for all lfm dump files to read. tsteps is the 
c number of lfm time steps per lfm dump steps (or per dump steps to read)
c first_step, start_step, & last_step must be evenly  divisible by tsteps

      parameter (tsteps=1)
*      parameter (tsteps=20)

c the total number of lfm dump steps in (or to read from) each lfm dump file is

      parameter (dsteps=1) 
*      parameter (dsteps=100) 

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

c     for lowpass version of LFM
      parameter (lowpass=.true.)

c SET TIME CORRECTION FACTOR AND ADDITIVE CONSTANT TO GET TIME RIGHT

      parameter (iyear0=2004) !YEAR NUMBER (FOUR DIGITS)
      parameter (idoy0=203) !DAY OF YEAR (DAY 1 = JAN 1)
      parameter (ihour0=0) !HOUR OF DAY (00 TO 23)
      parameter (imin0=0) !MINUTE OF HOUR (00 TO 59)
      parameter (isec0=0) !SECONDS OF MINUTE (00 TO 59)
c     iyear0,idoy0,imin0,isec0 + lfm time + time_add is UT time
*      parameter (time_add=15600.0) ! time_add = 0.0 gives sec from lfm run start
*      parameter (time_add=78060.0) ! time_add = 0.0 gives sec from lfm run start
*      parameter (time_add=84260.916) ! time_add = 0.0 gives sec from lfm run start
*      parameter (time_add=-2999.939)
*      parameter (time_add=-6060.0879)
      parameter (time_add=0.0)
      parameter (time_mult=1.0) ! set=1.0 to get seconds

c INTERIOR FIELDS (uses IGRF)
c the output will give zero fields for points inside inbound
c and outside outbound. 

c The field quantities inside of the LFM inner boundary (at ~2.2 or 3.2 RE)
c returned by hdftake are uncertian, so we set them zero in this wrapper code:
c The output is zero fields for points inside inbound and outside outbound.
c Also note that if we use IGRF for the interior field, and we want fields 
c down to 1 RE, then we must put field quantities on grid points at radius =
c 1RE minus 1.733*(grid cell delta)(Assuming delta x = delta y = delta z,
c 1.733*delta is the diagonal of a grid cell).

c     use with IGRF interior fields 
      parameter (inbound=1.0-1.733*.2)  !inner boundary
      parameter (outbound=30.0)  !outer boundary
      parameter (rmin_interp=2.21)  !inner boundary of interpolation region
      parameter (rmax_interp=3.21)  !outer boundary of interpolation region

c     use for zero field inside LFM inner boundary.
c     this makes inner boundary in rbelt code ~2.4 RE.
*      parameter (inbound=2.21)  !inner boundary
*      parameter (outbound=30.0)  !outer boundary
*      parameter (rmin_interp=1.21)  !inner boundary of interpolation region
*      parameter (rmax_interp=2.21)  !outer boundary of interpolation region

c Note that the interpolation between IGRF & LFM does not give us a divergence 
c free field in the interpolation region. The effect on test-particle results 
c is unknown to me. We should to find a better way.

      parameter(pi = 3.14159265358979)
      parameter(raddeg = 180./pi)

c *************************************************************************
c *************************************************************************

c     define output variables:
      real x(numx),y(numy),z(numz),t(nt_1)
      real bx(numx,numy,numz,nt+2),by(numx,numy,numz,nt+2)
      real bz(numx,numy,numz,nt+2),ex(numx,numy,numz,nt+2)
      real ey(numx,numy,numz,nt+2),ez(numx,numy,numz,nt+2)

c     additional lfm_sample variables:
      real bx_ex,by_ex,bz_ex,bx_in,by_in,bz_in,exwght,inwght
      integer i,j,k,l,m
      integer lstep(num_ifiles)
      integer nxtstep,ofile_num,ifile_num,index,halt
      real nx_tmp,ny_tmp,nz_tmp,b,time,t0 !initial time
      real Re, q, c, m0, b0, r, r2, r5, t_tmp
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

c     geopack routines
      integer mode
      real XGSM,YGSM,ZGSM
      real HXGSM,HYGSM,HZGSM,BXGSM,BYGSM,BZGSM
      REAL*8 tilt,X_dbl,Y_dbl,Z_dbl,BX_dbl,BY_dbl,BZ_dbl
      integer IYEAR,IDAY,IHOUR,MIN,ISEC
      real ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
     * CGST,SGST,BA(6),G(105),H(105),REC(105)
      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,
     * CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,
     * CGST,SGST,BA
      COMMON /GEOPACK2/ G,H,REC

c     Define constants
*      Re = 6378.137e5
*      q = 4.8032424e-10
*      c = 2.99792458e10
*      m0 = mass
*      ffactor=charge_sign*q*Re/(c*c*m0)
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
      print *,'inbound,outbound=',inbound,outbound
      print *,'rmin_interp,rmax_interp=',rmin_interp,rmax_interp
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
c     dsteps+first_step/tsteps= 100+192000/20
      lstep(1) = dsteps+first_step/tsteps
      print *
      print *,'lstep(1)=',lstep(1)
      do l=2,num_ifiles
         lstep(l) = lstep(l-1) + dsteps
         print *,'l,lstep(l)=',l,lstep(l)
      enddo
      lstep(num_ifiles)=last_step/tsteps
      print *,'num_ifiles,lstep(num_ifiles)=',
     &num_ifiles,lstep(num_ifiles)

*      stop

c     initialize loop
      nxtstep = start_step/tsteps      ! next step
      ifile_num = first_ifile   ! first input file number
      index = 0                 ! output file write index
      nstep = nt_1              ! # of steps in output file
      makegrid =.true.		! put grid info into a file
      ofile_num = start_step	! output file number

      print *
      print *,'total read/write steps =',(last_step-start_step)/tsteps+1

      ofile=basename//'ut.txt'
      open (10,file=ofile,status='unknown')
      write(10,*)iyear0,idoy0,ihour0,imin0,isec0
      close(10)

c     set input filename here
      call int2str(string_out,ifile_num,string_begin,string_end)
      ifile=basename//string_out(1:string_end)//'.hdf'

c     read and write loop starts here
100   if (nxtstep*tsteps .le. last_step) then

c       read from hdf input file until we reach end or we have enough
c       steps to fill output file
        halt=min0(nstep-index,lstep(ifile_num-first_ifile+1)-nxtstep+1)

*        print *,'read from hdf input file until we reach end or we' 
*        print *,'have enough steps to fill output file'
*        print *,'********************'
*        print *,'ifile_num-first_ifile+1=',ifile_num-first_ifile+1
*        print *,'lstep(ifile_num-first_ifile+1)=',
*     &  lstep(ifile_num-first_ifile+1)
*c       start_step/tsteps=195800/20
*        print *,'nxtstep=',nxtstep
*        print *,'********************'
*        print *,'nstep-index,lstep(ifile_num-first_ifile+1)-nxtstep+1=',
*     &  nstep-index,lstep(ifile_num-first_ifile+1)-nxtstep+1

c       but stop when we have read in the last_step
        halt=min0(halt,(last_step/tsteps-nxtstep+1))

*        print *,'but stop when we have read in the last_step'
*        print *,'halt,(last_step/tsteps-nxtstep+1)=',
*     &  halt,(last_step/tsteps-nxtstep+1)
*        print *,'nstep-index=',nstep-index
*        print *,'lstep(ifile_num-first_ifile+1)-nxtstep+1=',
*     &  lstep(ifile_num-first_ifile+1)-nxtstep+1
*        print *,'(last_step-nxtstep*tsteps+1)=',
*     &  (last_step-nxtstep*tsteps+1)

        print *,'halt=',halt

        do l=1,halt

          index = index + 1
c         need to fix code so that we always get mstep right
c         e.g., see May '97 CMIT run
          mstep = nxtstep*tsteps

          print *
          print *,'reading time step ',mstep,' from ',ifile
          print *
          print *,' in hdftake *************** '
c         get mhd values at points using hdftake
          ierr = hdftake(first,ifile,mstep,xlim,
     &         interp, allcart, points, npts, time)
*          ierr=0
          print *,' exit hdftake ************* '
          print *,'time=',time
          first=.false.
          if (ierr .ne. 0) then 
             write (6,*) 'error in hdftake'
             stop
          endif
c         UT = t (in seconds) + iyear0,idoy0,ihour0,imin0,isec0
*          t(index)=time*time_mult + time_add
          t(index)=3600.0*ifile_num

          iyear=iyear0
          iday=idoy0
          ihour=ihour0
          min=imin0
          isec=isec0+t(index)
          call reduce_time2(iyear,iday,ihour,min,isec)
          print *
          print *,'calling RECALC'
          print *,'iyear,iday,ihour,min,isec=',
     &    iyear,iday,ihour,min,isec
          call RECALC (iyear,iday,ihour,min,isec)
          print *,'done with RECALC: tilt=',atan2(SPS,CPS)*raddeg

*          stop

          do k=1,numz
*            print *,'k=',k
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
                    ex(i,j,k,index) = -allcart(9,i,j,k)
                    ey(i,j,k,index) = -allcart(10,i,j,k)
                    ez(i,j,k,index) = -allcart(11,i,j,k)
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
*                   print *,'calling SMGSM'
                    mode=1
                    CALL SMGSM(x(i),y(j),z(k),XGSM,YGSM,ZGSM,mode)
*                    print *
*                    print *,'x(i),y(j),z(k)=',x(i),y(j),z(k)
*                    print *,'XGSM,YGSM,ZGSM=',XGSM,YGSM,ZGSM
*                    print *,'r=',sqrt(XGSM**2+YGSM**2+ZGSM**2)
                    CALL IGRF_GSM(XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)
*                    print *,'HXGSM,HYGSM,HZGSM=',HXGSM,HYGSM,HZGSM
*                    print *,'B=',sqrt(HXGSM**2+HYGSM**2+HZGSM**2)
*                    print *,'calling SMGSM'
                    BXGSM=HXGSM+BX_dbl
                    BYGSM=HYGSM+BY_dbl
                    BZGSM=HZGSM+BZ_dbl
                    mode=-1
                    CALL SMGSM(bx_in,by_in,bz_in,BXGSM,BYGSM,BZGSM,mode)
                  else
                    exwght=1.0
                    bx_in=0.0
                    by_in=0.0
                    bz_in=0.0
                  endif
                  if (lowpass.eqv..true.) then
c                    for lowpass version of LFM
                     bx(i,j,k,index)=exwght*bx_ex + inwght*bx_in
                     by(i,j,k,index)=exwght*by_ex + inwght*by_in
                     bz(i,j,k,index)=exwght*bz_ex + inwght*bz_in
                  else
c                    dipole field gets removed
*                     bx(i,j,k,index)=exwght*bx_ex + inwght*bx_in +
*     &               3.*b0*x(i)*z(k)/r5
*                     by(i,j,k,index)=exwght*by_ex + inwght*by_in +
*     &               3.*b0*y(j)*z(k)/r5
*                     bz(i,j,k,index)=exwght*bz_ex + inwght*bz_in +
*     &               3.*b0*z(k)*z(k)/r5-b0/r2/r
                   endif
                else
                  bx(i,j,k,index) = 0.0
                  by(i,j,k,index) = 0.0
                  bz(i,j,k,index) = 0.0
                  ex(i,j,k,index) = 0.0
                  ey(i,j,k,index) = 0.0
                  ez(i,j,k,index) = 0.0
                endif
c                uncomment to get a pure dipole on grid (must also 
c                remove dipole from rbelt interp routine by setting b0=0)
*                bx(i,j,k,index) = 3.*b0*x(i)*z(k)/r5
*                by(i,j,k,index) = 3.*b0*y(j)*z(k)/r5
*                bz(i,j,k,index) = 3.*b0*z(k)*z(k)/r5
*     &                                        - b0/r2/r
*                ex(i,j,k,index) = 0.
*                ey(i,j,k,index) = 0.
*                ez(i,j,k,index) = 0.
              enddo
            enddo
          enddo

*c         uncomment to sanity check fields along x axis

          print *
          print *,'perturbation bz along x-axis'
          do i=1,numx
*            print *,'x,bz=',x(i),
*     &      bz(i,(numy+1)/2,(numz+1)/2,index)
            print *,'x,bz=',x(i),bz(i,(numy+1)/2,(numz+1)/2,index)
          enddo

          print *
          print *,'total bz along x-axis'
          do i=1,numx
*            print *,'x,bz=',x(i),
*     &      bz(i,(numy+1)/2,(numz+1)/2,index)
            print *,x(i),
     &      (bz(i,(numy+1)/2,(numz+1)/2,index)+b0/abs(x(i)**3))
          enddo
	  
          print *
          print *,'total ey along x-axis'
          do i=1,numx
            print *,'x,ey=',x(i),ey(i,(numy+1)/2,(numz+1)/2,index)
          enddo

          nxtstep = nxtstep + 1
        enddo    
c       if we have reached the end of the ifile then increment ifile_num
***        if (nxtstep .gt. lstep(ifile_num-first_ifile+1)) then
          first = .true.
          ifile_num = ifile_num + 1
          print *,'increment ifile_num to',ifile_num
          call int2str(string_out,ifile_num,string_begin,string_end)
          ifile=basename//string_out(1:string_end)//'.hdf'
***        endif
c       if we have read in nt hdf steps or we are done, then write out to HDF 
        if ((index.eq.nstep).or.(ifile_num.gt.last_ifile)) then
c         send x,y,z,t,bx,by,bz,ex,ey,ez out to HDF file
          rank=1
          start1d(1)=0
          stride1d(1)=1
          if (makegrid.eqv..true.) then
c           name and open the output file
            ofile=basename//'grid.hdf'
            print *
            print *,'writing out grid info to:',ofile
            file_id = sfstart(ofile,4)
            dim_sizes1d(1) = numx
            edge1d(1)=dim_sizes1d(1)
            sds_id=
     &      sfcreate(file_id,'x',DFNT_FLOAT32,rank,dim_sizes1d)
            status = sfwdata(sds_id,start1d,stride1d,edge1d,x)
            status = sfendacc(sds_id)
            dim_sizes1d(1) = numy
            edge1d(1)=dim_sizes1d(1)
            sds_id=
     &      sfcreate(file_id,'y',DFNT_FLOAT32,rank,dim_sizes1d)
            status = sfwdata(sds_id,start1d,stride1d,edge1d,y)
            status = sfendacc(sds_id)
            dim_sizes1d(1) = numz
            edge1d(1)=dim_sizes1d(1)            
            sds_id=
     &      sfcreate(file_id,'z',DFNT_FLOAT32,rank,dim_sizes1d)
            status = sfwdata(sds_id,start1d,stride1d,edge1d,z)
            status = sfendacc(sds_id)
*            status = sffinfo(file_id,n_datasets,n_file_attrs)
*            print *,'status=',status
*            print *,'file_id=',file_id
*            print *,'n_datasets=',n_datasets
*            print *,'n_file_attrs=',n_file_attrs
c           close hdf file
            status = sfend(file_id)
            if (status .eq. -1) then
              print*,'grid file close failed'
            endif
            makegrid = .false.
            nstep = nt
          endif
c         name and open the output file
          call int2str(string_out,ofile_num,string_begin,string_end)
          ofile=basename//'field-'//string_out(1:string_end)//'.hdf'
          print *
          print *,'writing',index,' steps to rbelt HDF file:',ofile
          file_id = sfstart(ofile,4)
          dim_sizes1d(1) = index
          edge1d(1)=dim_sizes1d(1)
          sds_id=sfcreate(file_id,'t',DFNT_FLOAT32,rank,dim_sizes1d)
          status = sfwdata(sds_id,start1d,stride1d,edge1d,t)
          status = sfendacc(sds_id)
c         output fields
          rank=4
          dim_sizes4d(1)=numx
          dim_sizes4d(2)=numy
          dim_sizes4d(3)=numz
          dim_sizes4d(4)=index
          do i=1,rank 
            stride4d(i)=1
            start4d(i)=0
            edge4d(i)=dim_sizes4d(i)
          enddo
          sds_id=sfcreate(file_id,'bx',DFNT_FLOAT32,rank,dim_sizes4d)
          status = sfwdata(sds_id,start4d,stride4d,edge4d,bx)
          status = sfendacc(sds_id)
          sds_id=sfcreate(file_id,'by',DFNT_FLOAT32,rank,dim_sizes4d)
          status = sfwdata(sds_id,start4d,stride4d,edge4d,by)
          status = sfendacc(sds_id)
          sds_id=sfcreate(file_id,'bz',DFNT_FLOAT32,rank,dim_sizes4d)
          status = sfwdata(sds_id,start4d,stride4d,edge4d,bz)
          status = sfendacc(sds_id)
          sds_id=sfcreate(file_id,'ex',DFNT_FLOAT32,rank,dim_sizes4d)
          status = sfwdata(sds_id,start4d,stride4d,edge4d,ex)
          status = sfendacc(sds_id)
          sds_id=sfcreate(file_id,'ey',DFNT_FLOAT32,rank,dim_sizes4d)
          status = sfwdata(sds_id,start4d,stride4d,edge4d,ey)
          status = sfendacc(sds_id)
          sds_id=sfcreate(file_id,'ez',DFNT_FLOAT32,rank,dim_sizes4d)
          status = sfwdata(sds_id,start4d,stride4d,edge4d,ez)
          status = sfendacc(sds_id)
c         close hdf file
          status = sfend(file_id)
          if (status .eq. -1) then
            print*,'field file close failed'
          endif	    
          ofile_num = nxtstep*tsteps	    
          index = 0
        endif
        goto 100
      endif

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

      FUNCTION julday2(doy,iyyy)
      INTEGER julday2,id,iyyy,mm,IGREG,doy
      PARAMETER (IGREG=15+31*(10+12*1582)) 
c Gregorian Calendar adopted Oct. 15, 1582.
c In this routine julday2 returns the Julian Day Number that begins at noon of the calendar
c date specifed by month mm,day id, and year iyyy, all integer variables. Positive year
c signifes A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
      INTEGER ja,jm,jy,status,mode

*      print *,'in julday2'
*      print *,'1, iyyy, doy, mm, id, status =',
*     &1,iyyy,doy,mm,id,status

      mode=1
      call doy_mmdd (mode, iyyy, doy, mm, id, status)
      if (status.ne.0) then
         print *,'status.ne.0 in subroutine doy_mmdd',status
         stop
      end if

*      print *,'iyyy,mm,id=',iyyy,mm,id

      jy=iyyy
      if (jy.eq.0) pause 'julday2: there is no year zero'
      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
         jm=mm+1
      else
         jy=jy-1
         jm=mm+13
      endif
      julday2=365*jy+int(0.25d0*jy+2000.d0)+int(30.6001d0*jm)+id+1718995
      if (id+31*(mm+12*iyyy).ge.IGREG) then !Test whether to change to Gregorian Calendar.
         ja=int(0.01d0*jy)
         julday2=julday2+2-ja+int(0.25d0*ja)
      endif

      return
      END

************************************************************************

      SUBROUTINE caldat2(julian,doy,iyyy)
      INTEGER id,iyyy,julian,mm,IGREG,doy
      PARAMETER (IGREG=2299161)
c Inverse of the function julday given above. Here julian is input as a Julian Day Number,
c and the routine outputs mm,id, and iyyy as the month, day, and year on which the specified
c Julian Day started at noon.
      INTEGER ja,jalpha,jb,jc,jd,je,status,mode

      mode=1
      call doy_mmdd (mode, iyyy, doy, mm, id, status)
      if (status.ne.0) then
         print *,'status.ne.0 in subroutine doy_mmdd',status
         stop
      end if

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

      mode=-1
      call doy_mmdd (mode, iyyy, doy, mm, id, status)
      if (status.ne.0) then
         print *,'status.ne.0 in subroutine doy_mmdd',status
         stop
      end if

      return
      END

************************************************************************

      subroutine reduce_time(year,mon,day,hour,min,sec)
      implicit none
      integer year,mon,day,hour,min,sec,jday,julday,t

*      print *
*      print *,' *** in subroutine reduce_time *** '
*      print *,'year,mon,day,hour,min,sec=',
*     &year,mon,day,hour,min,sec

      t=int(sec/60.)
      sec=sec-t*60
      min=min+t
      t=int(min/60.)
      min=min-t*60
      hour=hour+t
      t=int(hour/24.)
      hour=hour-t*24
      jday=julday(mon,day,year)+t
      call caldat(jday,mon,day,year)

      return
      END

************************************************************************

      subroutine reduce_time2(year,doy,hour,min,sec)
      implicit none
      integer year,doy,hour,min,sec,jday,julday2,t

*      print *
*      print *,' *** in subroutine reduce_time2 *** '
*      print *,'year,doy,hour,min,sec=',year,doy,hour,min,sec

      t=int(sec/60.)
      sec=sec-t*60
      min=min+t
      t=int(min/60.)
      min=min-t*60
      hour=hour+t
      t=int(hour/24.)
      hour=hour-t*24
      jday=julday2(doy,year)+t
      call caldat2(jday,doy,year)

      return
      END

************************************************************************

      subroutine sphr2cart(r,pol,phi,x,y,z)
      implicit none
      real r,pol,phi
      real x,y,z
      x = r*sin(pol)*cos(phi)
      y = r*sin(pol)*sin(phi)
      z = r*cos(pol)
      return
      end

************************************************************************

      subroutine cart2sphr(r,pol,phi,x,y,z)
      implicit none
      real r,pol,phi
      real x,y,z
      r=sqrt(x*x+y*y+z*z)
      pol=acos(z/r)
      if (y.gt.0) then
         phi=atan2(y,x)
      else
         phi=2.*3.14159265358979+atan2(y,x)
      endif

      return
      end
