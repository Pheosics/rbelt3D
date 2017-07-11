
      program wind_sample
	
c ************************************************************************* 
c
c Reads a user specified number of LFM HDF output files & puts E&B fields on 
c a cartisian grid. The output is a user specified number of HDF files for use 
c with the rbelt code. Note that the fields are scaled for use with the rbelt
c code in here. Also, the dipole term is removed from the magnetic field.
c The HDF output files are set up as follows:
c
c file *-rbelt-*.hdf has the sds sets
c 1:(nx)x
c 2:(ny)y
c 3:(nz)z
c 4:(nt)t
c 5:(nt x nx x ny x nz)bx 
c 6:(nt x nx x ny x nz)by 
c 7:(nt x nx x ny x nz)bz 
c 8:(nt x nx x ny x nz)ex 
c 9:(nt x nx x ny x nz)ey 
c 10:(nt x nx x ny x nz)ez 
c	
c B.Kress
c
c *************************************************************************

*      include 'hdf.inc'
*      include 'param.inc'
*      include 'intco1d.inc'

#include "hdf.inc"
#include "param.inc"
#include "intco1d.inc"

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
      integer last_step
      integer num_ifiles
      integer lfm_steps
      integer nsteps
      real inboundary
      real outboundary
      real dipole
      real time_mult
      real time_add
      real charge_sign
      real mass

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
c also, note that the inner boundary of the LFM grid is at 2.128451 Re
c and all points inside this will return zero fields.
c Do not forget to set parameter (num1d=numx*numy*numz) in intco1d.inc!!
c WARNING, we will need ~numx*numy*numz*nt*10 of RAM to run this, so 
c do not make the output file too big!!

c     num1d=1
      parameter (numx=1, numy=1, numz=1, nt=1) !nt is #of time steps in output file
      parameter (xmin=29., xmax=29.)
      parameter (ymin=0., ymax=0.)
      parameter (zmin=0., zmax=0.)

c INFORMATION ABOUT THE LFM DUMP FILE
c This info. must be uniform for all lfm dump files to read.

      parameter (lfm_steps=1)  !number of lfm time steps per lfm dump step
      parameter (nsteps=1)  !total number of lfm dump steps in each lfm dump file
      parameter (time_steps=(nsteps-1)*lfm_steps) 

c CHOOSE INPUT FILES (LFM DUMP FILES) AND THE STEPS TO READ
c here we do not count the first step in each lfm dump file,
c which is the same as the last step in the previous file.
c except for the initial step in file #1 which is step #1
c e.g. if the last step in file#1 is step 41 
c then step 42 is the second step in file#2.
c
c e.g. 
c Nov01-25 240250 250000 250 steps 962-1001
c Nov01-26 250250 260000 250 steps 1002-1041
c Nov01-27 260250 270000 250 steps 1042-1081
c Nov01-28 270250 280000 250 steps 1082-1121
c Nov01-29 280250 290000 250 steps 1122-1161
c etc..
c
c LFM dump step is (time step)/lfm_steps + 1
c
c set the lfm dump file names here e.g. if the first filename is Nov01-25.hdf
c then basename='Nov01-' and first_ifile=25

      parameter (basename='CMIT-agu06-')  !prefix for lfm dump file name
      parameter (first_ifile=1)  !first input file number
      parameter (last_ifile=63)  !last input file number
*      parameter (first_step=(first_ifile-1)*time_steps/lfm_steps+2)
*      parameter (last_step=last_ifile*time_steps/lfm_steps+1)
*      parameter (num_ifiles=last_ifile-first_ifile+1)  !do not change this

      parameter (first_step=1)
      parameter (last_step=63)
      parameter (num_ifiles=63)

c SOME INFORMATION ABOUT THE LFM RUN
c the output will give zero fields for points inside inboundary
c and outside outboundary. this could be modified below (look for condition with
c inboundary and outboundary) if we wish to use more of the LFM grid.

      parameter (inboundary=2.2)  !inner boundary of LFM grid
      parameter (outboundary=30.0)  !outer boundary of LFM grid (max radial)
c dipole must be the same as b0 to be used in rbelt code!!
      parameter (dipole=30000.0)  !bzero used for LFM dipole field 27500nt,31000nt or whatever.
*      parameter (dipole=0.0)

c SET TIME CORRECTION FACTOR AND ADDITIVE CONSTANT TO GET TIME RIGHT

      parameter (time_mult=1.0)  ! set=1.0 to get seconds
      parameter (time_add=0.0)  ! set=0.0 to get seconds from start of lfm run

c THE FIELDS NEED TO BE SCALED HERE FOR USE WITH THE RBELT,
c SO WE MUST CHOOSE ELECTRONS OR PROTONS.  

      parameter (charge_sign=-1)  !reverse charge_sign to run backwards in time
      parameter (mass=1.6726485e-24)  !proton mass
*      parameter (mass=9.109534e-28)  !electron mass

c *************************************************************************

c     define output variables:
      real x(numx),y(numy),z(numz),t(nt)
      real bx(numx,numy,numz,nt),by(numx,numy,numz,nt)
      real bz(numx,numy,numz,nt),ex(numx,numy,numz,nt)
      real ey(numx,numy,numz,nt),ez(numx,numy,numz,nt)
      real vx(numx,numy,numz,nt),vy(numx,numy,numz,nt)
      real vz(numx,numy,numz,nt),density(numx,numy,numz,nt)

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
      integer sfstart, sfcreate
      integer sfwdata, sfendacc, sfend
      integer*4 file_id, sds_id, rank, status
      integer dim_sizes4d(4), start4d(4), edge4d(4),stride4d(4)
      integer dim_sizes1d(1), start1d(1), edge1d(1),stride1d(1)

      integer sffinfo
      integer sfselect
      integer n_datasets,n_file_attrs

c     ...and a few more:
      integer i,j,k,l,m
      integer lstep(num_ifiles),fstep(num_ifiles)
      integer nxtstep,ofile_num,ifile_num,index,stop
      real nx_tmp,ny_tmp,nz_tmp,b,time,t0 !initial time
      real Re, q, c, m0, b0, r, r2, r5 
      real tscal,bscal,escal,ffactor
      character*80 ofile
      character*8 string_out
      integer string_begin,string_end

c     Define constants
c     to make sample more general we should scale the fields in rbelt code.
      Re = 6378.137e5
      q = 4.8032424e-10
      c = 2.99792458e10
      m0 = mass
      ffactor=charge_sign*q*Re/(c*c*m0)
      bscal=ffactor*1.e-5
      escal=ffactor*1.e6/c
      tscal = c/Re
      b0 = dipole

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

c     set up lstep array (last step to read in each lfm dump file)
      lstep(1) = (nsteps-1)*first_ifile + 1
      do l=2,num_ifiles
         lstep(l) = lstep(l-1) + nsteps-1
      enddo
      lstep(num_ifiles) = last_step

c     initialize loop
      nxtstep = first_step      ! next step
      ifile_num = first_ifile   ! first input file number
      ofile_num = 1             ! output file number
      index = 0                 ! output file write index

c     set input filename here
      call int2str(string_out,ifile_num,string_begin,string_end)
      ifile=basename//string_out(1:string_end)//'.hdf'

c     open output file
      open (12,file='sw_mhd.out',status='unknown')
      write(12,98) ((nsteps-1)*(last_ifile-first_ifile+1)),9
98    format (2i8)

c     read and write loop starts here
100   if (ifile_num .le. last_ifile) then
         print *,'ifile_nume=',ifile_num
c        if we are going to fill ofile before we reach the end of ifile
c        then stop when we fill ofile, else stop when we reach the end of ifile
         if ((nxtstep + nt -1) .le. lstep(ifile_num-first_ifile+1)) then
            stop = nt - index
         else
            stop = lstep(ifile_num-first_ifile+1)-nxtstep+1
         endif
c        read from hdf file until a) we reach end or b) we have enough to fill output file
         do l=1,stop
            index = index + 1
            mstep = (nxtstep-1)*lfm_steps	
*            print *,' mstep = ',mstep
c	      	get mhd values at points using hdftake
            print *
            print *,'reading dump step ',nxtstep,' from ',ifile
            print *
            print *,' in hdftake *************** '

            ierr = hdftake(first,ifile,mstep,xlim,
     &           interp, allcart, points, npts, time)

            print *,' exit hdftake ************* '
*            print *,'time=',time
            first=.false.
            if (ierr .ne. 0) then 
               write (6,*) 'error in hdftake'
               stop
            else
*               t(index)=(time*time_mult + time_add)*tscal
               t(index)=(time*time_mult + time_add)
*               t=(time*time_mult + time_add)
*               print *,'index=', index
               do k=1,numz
                  do j=1,numy
                     do i=1,numx
                        r2=x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
                        r=sqrt(r2)
                        r5=r2*r2*r
                        
                        if ((r.gt.inboundary)
c                          here we subtract off dipole field 
     &                       .and.(r.lt.outboundary))  then
                           bx(i,j,k,index) = allcart(6,i,j,k)
                           by(i,j,k,index) = allcart(7,i,j,k)
                           bz(i,j,k,index) = allcart(8,i,j,k)
                           ex(i,j,k,index) = allcart(9,i,j,k)
                           ey(i,j,k,index) = allcart(10,i,j,k)
                           ez(i,j,k,index) = allcart(11,i,j,k)
                           vx(i,j,k,index) = allcart(2,i,j,k)
                           vy(i,j,k,index) = allcart(3,i,j,k)
                           vz(i,j,k,index) = allcart(4,i,j,k)
                           density(i,j,k,index) = allcart(1,i,j,k)
                        else
                           bx(i,j,k,index) = 0
                           by(i,j,k,index) = 0
                           bz(i,j,k,index) = 0
                           ex(i,j,k,index) = 0
                           ey(i,j,k,index) = 0
                           ez(i,j,k,index) = 0
                           vx(i,j,k,index) = 0
                           vy(i,j,k,index) = 0
                           vz(i,j,k,index) = 0
                           density(i,j,k,index) = 0
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

            print *,'index, t(index)=',index, t(index)


            write(12,99) nxtstep,t(index),
     &         density(1,1,1,index),vx(1,1,1,index),vy(1,1,1,index),
     &         vz(1,1,1,index),bx(1,1,1,index),by(1,1,1,index),
     &         bz(1,1,1,index)
99	       format (i8,8f12.4)

            nxtstep = nxtstep + 1

c           uncomment to sanity check perturbation fields along x axis
*            print *
*            do i=1,numx
*                print *,'x,bz=',x(i),
*     &          bz(i,(numy+1)/2,(numz+1)/2,index)/bscal
*            enddo

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
         if ((index.eq.nt).or.(ifile_num.gt.last_ifile)) then
c        send x,y,z,t,bx,by,bz,ex,ey,ez out to HDF file ************
            ofile_num = ofile_num + 1
            index = 0
         endif
         
         goto 100

      endif

      close(12)

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

	
	
	
