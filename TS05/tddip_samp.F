
      program lfm_sample
	
c ************************************************************************* 
c
c Outputs is a user specified number of HDF files for use with the rbelt code. 
c Note The dipole term is removed from the magnetic field.
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

      implicit none
      include 'hdf.inc'
      include 'dffunc.inc'

c     declare parameter types
      character* (*) basename
      integer numx,numy,numz,nt,nt_1,nstep,first_filenum
      real xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz
      logical tsy,times_file

      real b0,inbound,outbound
      integer tdelta,start_sec,end_sec
      integer year0,doy0,hour0,min0,sec0,sec_tmp,max_table
      integer julday2,julday,max_table2

      real pi,raddeg

c *************************************************************************
c *************************************************************************
c BELOW ARE THE PARAMETERS YOU NEED TO SET:

c     nt is #of time steps in output file except 1st output file 
c     which gets nstep=nt+2

c all points inside inboundary (defined below) + 0.001 RE 
c will return zero fields.
c Do not forget to set parameter (num1d=numx*numy*numz) in intco1d.inc!
c WARNING, size of output file = numx*numy*numz*nt*6*4 & you will need
c approx. numx*numy*numz*nt*6*4*2 bytes of RAM to run this code.
c nt is #of time steps in the output files, except 1st output file 
c has nt_1 time steps. Note, grid delta x =(xmax-xmin)/(numx-1)

      parameter (numx=101, numy=101, numz=101, nt=3) ! parameter (num1d=8120601)
      parameter (dx=0.2,dy=0.2,dz=0.2)
      parameter (xmin=-10., xmax=10.)
      parameter (ymin=-10., ymax=10.)
      parameter (zmin=-10., zmax=10.)
      parameter (nt_1=5)
      parameter (first_filenum=1)
      parameter (basename='tddip-')  !prefix for file name

c     SET MAIN DIPOLE AND TIME DEPENDENT PERTURBATION HERE
c     subtracted dipole must be the same as b0 to be used in rbelt code!!
c     use bzero from LFM dipole field 27500nt,31000nt or whatever.
      parameter (b0=30000.0)
      real*8 b0td,db0,tau,e0,bfac1
      common /tdepB0/ b0td,e0
c     time constant (seconds)
c     fast
      parameter(tau=100)
c     slow
*      parameter(tau=20000)
c     final delta B0
      parameter(db0=1000.0)

c      if we use b0td=b0*(1. + db0*(1.-exp(-t(n)/tau))) below,
c      resulting electric field is:
c      rho=dsqrt(x*x + y*y)
c      ls=r*r*r/rho**2
c      e0=db0*dabs(b0)*exp(-t/tau)/ls/rho/tau (dimensionless form)
c      (use e0/ffactor/vmsvcm in rbelt code to get V/M)
c      cphi = x/rho
c      sphi = y/rho
c      ex = e0*sphi
c      ey = -e0*cphi
      
c SET START TIME & END TIME (all integers)
      parameter (year0=2000) !YEAR NUMBER (FOUR DIGITS)
      parameter (doy0=1) !DAY OF YEAR (DAY 1 = JAN 1)
      parameter (hour0=0) !HOUR OF DAY (00 TO 23)
      parameter (min0=0) !MINUTE OF HOUR (00 TO 59)
      parameter (sec0=0) !SECONDS OF MINUTE (00 TO 59)

c       set output fields start time
      parameter (start_sec=0) ! addiitonal sec. beyond sec0 above to start
      parameter (end_sec=10) ! requested end time (sec.)
c       set time step
      parameter (tdelta=1) ! delta t
c       or use input file to get times

c     i.e., 1st time step is at reduce_time(year0,doy0,hour0,min0,(sec0+start_sec))
c     year0,doy0,hour0,min0,sec0 goes in UT text file for later reference by rbelt code
c     start_sec + n*tdelta is time stamp in HDF file, where n*tdelta = 0, tdelta, 2*tdelta...
c     up to end_sec, unless we use a `cutoff_times' input table

      parameter (inbound=.1)  !inner boundary
      parameter (outbound=20.0)  !outer boundary

      parameter(pi = 3.14159265358979)
      parameter(raddeg = 180./pi)

c *************************************************************************
c *************************************************************************

c     define output variables:
      real x(numx),y(numy),z(numz),t(nt_1)
      real bx(numx,numy,numz,nt+2),by(numx,numy,numz,nt+2)
      real bz(numx,numy,numz,nt+2),ex(numx,numy,numz,nt+2)
      real ey(numx,numy,numz,nt+2),ez(numx,numy,numz,nt+2)
      real bxtot,bytot,bztot

c     additional lfm_sample variables:
      integer i,j,k,l,m,n
      integer index,tindex,num_files,total_tsteps,last_sec
c       end_sec is requested last time step (in seconds from sec0)
c       note, last_sec is seconds of last time step (in seconds from sec0)
      real nx_tmp,ny_tmp,nz_tmp,b
      real Re, q, c, m0, r, r2, r5, t_tmp
      real tscal,bscal,escal,ffactor
      character*80 ofile
      character*8 string_out
      integer string_begin,string_end
      logical makegrid

      integer stsod,ndsod,stsod_tbl,ndsod_tbl
      integer stday,ndday,stday_tbl,ndday_tbl
      integer stsec_tbl,ndsec_tbl,sec_frm_ref

c     HDF stuff:
      integer*4 file_id, sds_id, rank, status
      integer dim_sizes4d(4), start4d(4), edge4d(4),stride4d(4)
      integer dim_sizes1d(1), start1d(1), edge1d(1),stride1d(1)
      integer n_datasets, n_file_attrs

c     define hdftake in and out variables:
c     See note at top of hdftake.f file if you wish to 
c     see how the mhd variables are put into the array allcart(14,nx,ny,nz).
      integer hdftake
      logical first,interp
      character*80 ifile
      integer mstep,npts
      common /big/bx,by,bz,ex,ey,ez

c     geopack & Tsyganenko routines
      real delta,sec1,sec2
      integer mode
      integer year,doy,hour,min,sec

      print *
      print *,'running program lfm_sample'
      print *,'numx,numy,numz=',numx,numy,numz
      print *,'xmin,xmax=',xmin,xmax
      print *,'ymin,ymax=',ymin,ymax
      print *,'zmin,zmax=',zmin,zmax
      print *,'basename=',basename
      print *,'b0=',b0
      print *,'year0,doy0,hour0,min0,sec0=',
     &year0,doy0,hour0,min0,sec0
      print *,'start_sec,end_sec,tdelta=',start_sec,end_sec,tdelta
      print *,'inbound,outbound=',inbound,outbound
      print *,'b0=',b0

      npts=numx*numy*numz

*c     make sure proper space is allocated for arrays
*      if (npts.ne.num1d) then
*         print *,'npts.ne.num1d'
*         print *,'in intco1d.inc, set num1d=',npts
*         stop
*      endif

c     set up points
      nx_tmp = numx
      ny_tmp = numy
      nz_tmp = numz
      if (numx.eq.1) nx_tmp=2
      if (numy.eq.1) ny_tmp=2
      if (numz.eq.1) nz_tmp=2
      do i=1,numx
*         x(i) = xmin + (i-1)*(xmax-xmin)/(nx_tmp-1)
         x(i) = xmin + (i-1.)*dx
         print *,'(i-1)*dx,x(',i,')=',(i-1)*dx,x(i)
      enddo
      do i=1,numy
*         y(i) = ymin + (i-1)*(ymax-ymin)/(ny_tmp-1)
         y(i) = ymin + (i-1)*dy
         print *,'y(',i,')=',y(i)
      enddo
      do i=1,numz
*         z(i) = zmin + (i-1)*(zmax-zmin)/(nz_tmp-1)
         z(i) = zmin + (i-1)*dz
         print *,'z(',i,')=',z(i)
      enddo

      ofile=basename//'ut.txt'
      open (10,file=ofile,status='unknown')
      write(10,*)year0,doy0,hour0,min0,sec0
      close(10)

c        to get uniform increments in time starting at start_sec
c        requested start and end times
c        start time
         year=year0
         doy=doy0
         hour=hour0
         min=min0
         sec=sec0+start_sec
         call reduce_time2(year,doy,hour,min,sec)
         print *
         print *,'requested start year,doy,hour,min,sec=',
     &   year,doy,hour,min,sec
c        end time
         year=year0
         doy=doy0
         hour=hour0
         min=min0
         sec=sec0+end_sec
         call reduce_time2(year,doy,hour,min,sec)
         print *
         print *,'requested end year,doy,hour,min,sec=',
     &   year,doy,hour,min,sec

c        add tdelta incriments to start time
         nstep=nt_1
         total_tsteps=0
         tindex=0
         index=0
         num_files=0
         last_sec=start_sec
*         print *,'num_files,total_tsteps,index,last_sec=',
*        &num_files,total_tsteps,index,last_sec
30       if ((start_sec+tindex*tdelta).le.end_sec) then
            num_files=num_files+1
40          if (index.lt.nstep) then
               index=index+1
               last_sec=start_sec+tindex*tdelta
               tindex=tindex+1
               total_tsteps=total_tsteps+1
               print *
               print *,'file num.,index,total_tsteps,second=',
     &          num_files+first_filenum-1,index,total_tsteps,last_sec
               year=year0
               doy=doy0
               hour=hour0
               min=min0
               sec=sec0+last_sec
               call reduce_time2(year,doy,hour,min,sec)
               print *,'start year,doy,hour,min,sec=',
     &          year,doy,hour,min,sec
               goto 40
            endif
            if (nstep.ne.nt) nstep=nt
            index=0
            goto 30
         endif
         print *
         print *,'num_files,total_tsteps,last_sec=',
     &   num_files,total_tsteps,last_sec

c     read and write loop starts here
c     initialize loop
      tindex=0  ! time index
      nstep = nt_1  ! # of steps in 1st output file
      makegrid =.true.  ! put grid info into a file
*        OPEN(53, FILE='test_dbdt.txt', STATUS='OLD')

      do l=1,num_files
         do n=1,nstep
            print *
            print *,'file num., nstep=',l,n
c           below is num. of sec. from year0, doy0, hour0, min0, sec0
            t(n)=start_sec+tindex*tdelta
            tindex=tindex+1

c           perturbation b0 is time dependent:
c           THERE IS NOTHING SPECIAL ABOUT THIS FORM FOR THE TIME DEPENDENCE
c           THIS CAN BE MODIFIED TO A CONSTANT INCREASE/DECREASE OR WHATTEVER,
c           BUT WE NEED TO INCLUDE MAIN DIP + TIME-DEP. PERTURBATION AS BELOW
            b0td=b0*(1+t(n)/tau)

            do k=1,numz
*              print *,'k=',k
              do j=1,numy
                do i=1,numx
                  r2=x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
                  r=sqrt(r2)
                  r5=r2*r2*r

c                 may need to move inner boundary in to get e-field calc. right
c                 in fact, may have a problem since B -> infinity at origin
c                 can we ignore contribution as r -> 0 ??
                  if ((r.le.outbound).and.(r.ge.inbound)) then
		  
		    bfac1=3.*b0td/r2/r2/r
                    bxtot=-bfac1*x(i)*z(k)
                    bytot=-bfac1*y(j)*z(k)
                    bztot=-bfac1*z(k)*z(k)+b0td/r2/r

                    
c                   main B0 dipole field gets removed
                    bx(i,j,k,n)=(bxtot+3.*b0*x(i)*z(k)/r5)
                    by(i,j,k,n)=(bytot+3.*b0*y(j)*z(k)/r5)
                    bz(i,j,k,n)=(bztot+3.*b0*z(k)*z(k)/r5-b0/r2/r)
                    ex(i,j,k,n)=0.
                    ey(i,j,k,n)=0.
                    ez(i,j,k,n)=0.

                  else
                    bx(i,j,k,n) = 0.0
                    by(i,j,k,n) = 0.0
                    bz(i,j,k,n) = 0.0
                    ex(i,j,k,n) = 0.0
                    ey(i,j,k,n) = 0.0
                    ez(i,j,k,n) = 0.0
                  endif
*                  if ((l.eq.1).AND.(n.eq.2)) then
*                     write(53,*) bx(i,j,k,2), by(i,j,k,2), bz(i,j,k,2)
*                  endif
                enddo
              enddo
            enddo

            print *
            print *,'total bz at y,z=',y((numy+1)/2),z((numz+1)/2)
            do i=1,numx
*              print *,'x,bz=',x(i),
*     &        bz(i,(numy+1)/2,(numz+1)/2,n)
              print *,x(i),
     &        (bz(i,(numy+1)/2,(numz+1)/2,n)+b0/abs(x(i)**3))
            enddo

         enddo

c        send x,y,z,t,bx,by,bz,ex,ey,ez out to HDF file
         rank=1
         start1d(1)=0
         stride1d(1)=1
         if (makegrid.eqv..true.) then
c          name and open the output file
           ofile=basename//'grid.hdf'
           print *
           print *,'writing out grid info to:',ofile
           file_id = sfstart(ofile,4)
           dim_sizes1d(1) = numx
           edge1d(1)=dim_sizes1d(1)
           sds_id=
     &     sfcreate(file_id,'x',DFNT_FLOAT32,rank,dim_sizes1d)
           status = sfwdata(sds_id,start1d,stride1d,edge1d,x)
           status = sfendacc(sds_id)
           dim_sizes1d(1) = numy
           edge1d(1)=dim_sizes1d(1)
           sds_id=
     &     sfcreate(file_id,'y',DFNT_FLOAT32,rank,dim_sizes1d)
           status = sfwdata(sds_id,start1d,stride1d,edge1d,y)
           status = sfendacc(sds_id)
           dim_sizes1d(1) = numz
           edge1d(1)=dim_sizes1d(1)            
           sds_id=
     &     sfcreate(file_id,'z',DFNT_FLOAT32,rank,dim_sizes1d)
           status = sfwdata(sds_id,start1d,stride1d,edge1d,z)
           status = sfendacc(sds_id)
*           status = sffinfo(file_id,n_datasets,n_file_attrs)
*           print *,'status=',status
*           print *,'file_id=',file_id
*           print *,'n_datasets=',n_datasets
*           print *,'n_file_attrs=',n_file_attrs
c          close hdf file
           status = sfend(file_id)
           if (status .eq. -1) then
             print*,'grid file close failed'
           endif
           makegrid = .false.
         endif
c        name and open the output file
         call int2str
     &    (string_out,(l-1+first_filenum),string_begin,string_end)
         ofile=basename//'field-'//string_out(1:string_end)//'.hdf'
         print *
         print *,'writing',nstep,' steps to rbelt HDF file:',ofile
         file_id = sfstart(ofile,4)
         dim_sizes1d(1) = nstep
         edge1d(1)=dim_sizes1d(1)
         sds_id=sfcreate(file_id,'t',DFNT_FLOAT32,rank,dim_sizes1d)
         status = sfwdata(sds_id,start1d,stride1d,edge1d,t)
         status = sfendacc(sds_id)
c        output fields
         rank=4
         dim_sizes4d(1)=numx
         dim_sizes4d(2)=numy
         dim_sizes4d(3)=numz
         dim_sizes4d(4)=nstep
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
c        close hdf file
         status = sfend(file_id)
         if (status .eq. -1) then
           print*,'field file close failed'
         endif	    
         if (nstep.ne.nt) nstep = nt
      enddo

      print*,t
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

	SUBROUTINE INT2STR(CAOUT,IVALUE,IBEGIN,IEND)

	character*(*) caout
	INTEGER iValue, iBegin, iEnd
	INTEGER i,iStart,iLast
	character*20 caDum


	WRITE(caDum,10) iValue
 10     FORMAT(i8)

	iStart = 1
	DO while( caDum(iStart:iStart) .eq. ' ' )
	  iStart = iStart + 1
	ENDDO
*#ifdef AIX
*	iLast = LEN_TRIM(caDum)
*#endif
*#ifndef AIX
 	iLast = LNBLNK(caDum)
*#endif

	caOut(1:iLast-iStart+1) = caDum(iStart:iLast)
        iBegin = 1
	iEnd = iLast - iStart + 1

	RETURN
	END

