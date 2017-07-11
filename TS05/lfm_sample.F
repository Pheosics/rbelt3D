
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
      integer numx,numy,numz,nt
      real xmin,xmax
      real ymin,ymax
      real zmin,zmax
      character* (*) basename
      integer first_ifile
      integer last_ifile
      integer first_step
      integer start_step
      integer last_step
      integer num_ifiles
      integer tsteps
      integer dsteps
      real inboundary
      real outboundary
      real dipole
      real time_mult
      real time_add
      real charge_sign
      real mass
      logical scale

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

      parameter (numx=121, numy=121, numz=121, nt=5) !nt is #of time steps in output file
      parameter (xmin=-12., xmax=12.)
      parameter (ymin=-12., ymax=12.)
      parameter (zmin=-12., zmax=12.)

c      num1d=3442951
*      parameter (numx=151, numy=151, numz=151, nt=1) !nt is #of time steps in output file
*      parameter (xmin=-15., xmax=15.)
*      parameter (ymin=-15., ymax=15.)
*      parameter (zmin=-15., zmax=15.)

c For the LFM grid with (i,j,k) = (53,24,32) we have:
c max z =   126.5580    
c max y =   126.0819    
c max x =   30.03723    
c min z =  -126.5580    
c min y =  -126.0819    
c min x =  -345.2204   
c so we should stay within this outer boundary.
c also, note that the inner boundary of the LFM grid is at ~2.2 Re,
c and all points inside inboundary (defined below) + 0.001 RE 
c will return zero fields.
c Do not forget to set parameter (num1d=numx*numy*numz) in intco1d.inc!
c WARNING, size of output file = numx*numy*numz*nt*6*4 & you will need
c approx. numx*numy*numz*nt*6*4*2 bytes of RAM to run this code.
c nt is #of time steps  output files, except 1st output file 
c has nstep=nt+2 time steps.


c CHOOSE INPUT FILES (LFM DUMP FILES) AND THE STEPS TO READ
c set the lfm dump file names here e.g. if the first filename is Nov01-25.hdf
c then basename='Nov01-' and first_ifile=25

      parameter (basename='Halloween_Storm-')  !prefix for lfm dump file name
      parameter (first_ifile=1)  !first input file number
      parameter (last_ifile=1)  !last input file number
      parameter (first_step=130000)  !first lfm time step in first_ifile
      parameter (start_step=136400)  !first lfm time step to read
      parameter (last_step=136400)  !last lfm time step to read
      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this

c ADDITIONAL INFORMATION ABOUT THE LFM DUMP FILES
c This info. must be uniform for all lfm dump files to read. tsteps is the 
c number of lfm time steps per lfm dump steps (or per dump steps to read)
c first_step, start_step, & last_step must be evenly  divisible by tsteps

      parameter (tsteps=100) 

c the total number of lfm dump steps in (or to read from) each lfm dump file is

      parameter (dsteps=100) 

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
      
c SOME INFORMATION ABOUT THE LFM RUN
c the output will give zero fields for points inside inboundary
c and outside outboundary. this could be modified below (look for condition with
c inboundary and outboundary) if we wish to use more of the LFM grid.

      parameter (inboundary=2.21)  !inner boundary of LFM grid
      parameter (outboundary=30.0)  !outer boundary
c     dipole must be the same as b0 to be used in rbelt code!!
c     bzero used for LFM dipole field 27500nt,31000nt or whatever.
      parameter (dipole=30000.0)
*      parameter (dipole=0.0)

c SET TIME CORRECTION FACTOR AND ADDITIVE CONSTANT TO GET TIME RIGHT

      parameter (time_mult=1.0)  ! set=1.0 to get seconds
      parameter (time_add=0.0)  ! set=0.0 to get seconds from start of lfm run

c THE FIELDS NEED TO BE SCALED HERE FOR USE WITH THE RBELT,
c SO WE MUST CHOOSE ELECTRONS OR PROTONS.  

      parameter (scale=.false.) !scale fields
      parameter (charge_sign=-1)  !reverse charge_sign to run backwards in time
*      parameter (mass=1.6726485e-24)  !proton mass
      parameter (mass=9.109534e-28)  !electron mass

c *************************************************************************

c     define output variables:
      real x(numx),y(numy),z(numz),t(nt+2)
      real bx(numx,numy,numz,nt+2),by(numx,numy,numz,nt+2)
      real bz(numx,numy,numz,nt+2),ex(numx,numy,numz,nt+2)
      real ey(numx,numy,numz,nt+2),ez(numx,numy,numz,nt+2)

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

c     ...and a few more:
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

c     Define constants
      Re = 6378.137e5
      q = 4.8032424e-10
      c = 2.99792458e10
      m0 = mass
      ffactor=charge_sign*q*Re/(c*c*m0)
      b0 = dipole

      print *
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

c     read and write loop starts here
100   if (nxtstep*tsteps .le. last_step) then
c        read from hdf input file until we reach end or we have enough
c        steps to fill output file
         halt=min0((nstep-index),
     &   lstep(ifile_num-first_ifile+1)-nxtstep+1)
c        but stop when we have read in the last_step
         halt=min0(halt,(last_step/tsteps-nxtstep+1))
*         print *,'nstep-index=',nstep-index
*         print *,'lstep(ifile_num-first_ifile+1)-nxtstep+1=',
*      &  lstep(ifile_num-first_ifile+1)-nxtstep+1
*         print *,'(last_step-nxtstep*tsteps+1)=',
*      &  (last_step-nxtstep*tsteps+1)
*         print *,'halt=',halt
         do l=1,halt
            index = index + 1
            mstep = nxtstep*tsteps
            print *
            print *,'reading time step ',mstep,' from ',ifile
            print *
            print *,' in hdftake *************** '
c           get mhd values at points using hdftake
            ierr = hdftake(first,ifile,mstep,xlim,
     &           interp, allcart, points, npts, time)
            print *,' exit hdftake ************* '
*            print *,'time=',time
            first=.false.
            if (ierr .ne. 0) then 
               write (6,*) 'error in hdftake'
               stop
            else
               t(index)=(time*time_mult + time_add)*tscal
               do k=1,numz
                  do j=1,numy
                     do i=1,numx
                        r2=x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
                        r=sqrt(r2)
                        r5=r2*r2*r
                        
                        if ((r.gt.(inboundary+0.001))
     &                       .and.(r.lt.outboundary))  then
c                          here we subtract off dipole field 
                           bx(i,j,k,index) = (allcart(6,i,j,k) + 
     &                          3.*b0*x(i)*z(k)/r5)*bscal
                           by(i,j,k,index) = (allcart(7,i,j,k) + 
     &                          3.*b0*y(j)*z(k)/r5)*bscal
                           bz(i,j,k,index) = (allcart(8,i,j,k) + 
     &                          3.*b0*z(k)*z(k)/r5 - b0/r2/r)*bscal
                           ex(i,j,k,index) = allcart(9,i,j,k)*escal
                           ey(i,j,k,index) = allcart(10,i,j,k)*escal
                           ez(i,j,k,index) = allcart(11,i,j,k)*escal
                        else
                           bx(i,j,k,index) = 0
                           by(i,j,k,index) = 0
                           bz(i,j,k,index) = 0
                           ex(i,j,k,index) = 0
                           ey(i,j,k,index) = 0
                           ez(i,j,k,index) = 0
                        endif
                        
c                        uncomment to get a pure dipole on grid
c                        (must also remove dipole from rbelt interp routine)
*                        bx(i,j,k,index) = (3.*b0*x(i)*z(k)/r5)*bscal
*                        by(i,j,k,index) = (3.*b0*y(j)*z(k)/r5)*bscal
*                        bz(i,j,k,index) = (3.*b0*z(k)*z(k)/r5
*     &                                               - b0/r2/r)*bscal
*                        ex(i,j,k,index) = 0.
*                        ey(i,j,k,index) = 0.
*                        ez(i,j,k,index) = 0.

                     enddo
                  enddo
               enddo
            endif

c           uncomment to sanity check perturbation fields along x axis
*            print *
*            do i=1,numx
*                print *,'x,bz=',x(i),
*     &          bz(i,(numy+1)/2,(numz+1)/2,index)/bscal
*            enddo

            nxtstep = nxtstep + 1
	    
         enddo    
c        if we have reached the end of the ifile then increment ifile_num
         if (nxtstep .gt. lstep(ifile_num-first_ifile+1)) then
            first = .true.
            ifile_num = ifile_num + 1
*            print *,'increment ifile_num to',ifile_num
            call int2str(string_out,ifile_num,string_begin,string_end)
            ifile=basename//string_out(1:string_end)//'.hdf'
         endif
c        if we have read in nt hdf steps or we are done, then write out to HDF 
         if ((index.eq.nstep).or.(ifile_num.gt.last_ifile)) then
c         send 	x,y,z,t,bx,by,bz,ex,ey,ez out to HDF file ************

            rank=1
            start1d(1)=0
            stride1d(1)=1

            if (makegrid.eqv..true.) then

c              name and open the output file
               ofile=basename//'grid.hdf'
               print *
               print *,'writing out grid info to:',ofile
               file_id = sfstart(ofile,4)
	       
	       dim_sizes1d(1) = numx
               edge1d(1)=dim_sizes1d(1)
               sds_id=
     &         sfcreate(file_id,'x',DFNT_FLOAT32,rank,dim_sizes1d)
               status = sfwdata(sds_id,start1d,stride1d,edge1d,x)
               status = sfendacc(sds_id)
        
               dim_sizes1d(1) = numy
               edge1d(1)=dim_sizes1d(1)
               sds_id=
     &         sfcreate(file_id,'y',DFNT_FLOAT32,rank,dim_sizes1d)
               status = sfwdata(sds_id,start1d,stride1d,edge1d,y)
               status = sfendacc(sds_id)
            
               dim_sizes1d(1) = numz
               edge1d(1)=dim_sizes1d(1)            
               sds_id=
     &         sfcreate(file_id,'z',DFNT_FLOAT32,rank,dim_sizes1d)
               status = sfwdata(sds_id,start1d,stride1d,edge1d,z)
               status = sfendacc(sds_id)

*               status = sffinfo(file_id,n_datasets,n_file_attrs)
*               print *,'status=',status
*               print *,'file_id=',file_id
*               print *,'n_datasets=',n_datasets
*               print *,'n_file_attrs=',n_file_attrs

c              close hdf file
               status = sfend(file_id)
               if (status .eq. -1) then
                 print*,'grid file close failed'
               endif

               makegrid = .false.
               nstep = nstep - 2

            endif

c           name and open the output file
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
            
c           output fields

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
            
c           close hdf file
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

	
	
	
