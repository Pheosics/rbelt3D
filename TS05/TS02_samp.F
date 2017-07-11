
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

      include 'hdf.inc'
      include 'dffunc.inc'

c     declare parameter types
      implicit none
      character* (*) basename,tsy_input
      integer numx,numy,numz,nt,nstep
      real xmin,xmax,ymin,ymax,zmin,zmax
      logical tsy

      real dipole,tdelta,end_time,inbound,outbound

      integer iyear0,idoy0,ihour0,imin0,isec0,max_table
      integer julday2,julday

c *************************************************************************
c *************************************************************************
c BELOW ARE THE PARAMETERS YOU NEED TO SET:

c     nt is #of time steps in output file except 1st output file 
c     which gets nstep=nt+2

      parameter (numx=121, numy=121, numz=121, nt=6)
      parameter (xmin=-12., xmax=12.)
      parameter (ymin=-12., ymax=12.)
      parameter (zmin=-12., zmax=12.)

      parameter (basename='Tsy02-')  !prefix for file name

c     subtracted dipole must be the same as b0 to be used in rbelt code!!
c     use bzero from LFM dipole field 27500nt,31000nt or whatever.
      parameter (dipole=29500.0)

c SET START TIME 

      parameter (iyear0=2001) !YEAR NUMBER (FOUR DIGITS)
      parameter (idoy0=328) !DAY OF YEAR (DAY 1 = JAN 1)
      parameter (ihour0=6) !HOUR OF DAY (00 TO 23)
      parameter (imin0=0) !MINUTE OF HOUR (00 TO 59)
      parameter (isec0=0) !SECONDS OF MINUTE (00 TO 59)
c     iyear0,idoy0,imin0,isec0 + lfm time + tdelta is UT start time
      parameter (tdelta=30.0) ! delta t
      parameter (end_time=500.0) ! end time (sec.)
      parameter (inbound=1.0)  !inner boundary
      parameter (outbound=30.0)  !outer boundary
      parameter (tsy=.true.) !add Tsyganenko to IGRF
      parameter (tsy_input='T02_input.dat') ! Tsyganenko/geopack input file name
      parameter (max_table=5000)

c *************************************************************************
c *************************************************************************

c     define output variables:
      real x(numx),y(numy),z(numz),t(nt+2)
      real bx(numx,numy,numz,nt+2),by(numx,numy,numz,nt+2)
      real bz(numx,numy,numz,nt+2),ex(numx,numy,numz,nt+2)
      real ey(numx,numy,numz,nt+2),ez(numx,numy,numz,nt+2)

c     additional lfm_sample variables:
      real bx_ex,by_ex,bz_ex,bx_in,by_in,bz_in,exwght,inwght
      integer i,j,k,l,m,n
      integer ofile_num,index,num_files,num_tsteps
      real nx_tmp,ny_tmp,nz_tmp,b,time,t0 !initial time
      real Re, q, c, m0, b0, r, r2, r5, t_tmp
      real tscal,bscal,escal,ffactor
      character*80 ofile
      character*8 string_out
      integer string_begin,string_end
      logical makegrid

      integer stsec0,ndsec0,stsec,ndsec
      integer stday0,ndday0,stday,ndday

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

      print *
      print *,'running program lfm_sample'
      print *,'numx,numy,numz=',numx,numy,numz
      print *,'xmin,xmax=',xmin,xmax
      print *,'ymin,ymax=',ymin,ymax
      print *,'zmin,zmax=',zmin,zmax
      print *,'basename=',basename
*      print *,'first_step,start_step,last_step=',
*     &first_step,start_step,last_step
*      print *,'tsteps,dsteps=',tsteps,dsteps
      print *,'dipole=',dipole
      print *,'iyear0,idoy0,ihour0,imin0,isec0=',
     &iyear0,idoy0,ihour0,imin0,isec0
      print *,'tdelta=',tdelta
      print *,'inbound,outbound=',inbound,outbound
      print *,'tsy,tsy_input=',tsy,' ',tsy_input
      print *,'max_table=',max_table
      print *,'b0=',b0

      npts=numx*numy*numz

c     set up points
      nx_tmp = numx
      ny_tmp = numy
      nz_tmp = numz
      if (numx.eq.1) nx_tmp=2
      if (numy.eq.1) ny_tmp=2
      if (numz.eq.1) nz_tmp=2
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

      ofile=basename//'ut.txt'
      open (10,file=ofile,status='unknown')
      write(10,*)iyear0,idoy0,ihour0,imin0,isec0
      close(10)      

      nstep=nt+2
      num_tsteps=1
      index=1
      num_files=0
      time=0
*      print *,'num_files,num_tsteps,index,time=',
*     &num_files,num_tsteps,index,time
30    if (time.lt.end_time) then
         num_files=num_files+1
40       if (index.lt.nstep) then
            index=index+1
            num_tsteps=num_tsteps+1
            time=time+tdelta
*            print *,'num_files,num_tsteps,index,time=',
*     &      num_files,num_tsteps,index,time
            goto 40
         endif
         if (nstep.ne.nt) nstep=nt
         index=0
         goto 30
      endif
      print *,'num_files,num_tsteps,tmax=',num_files,num_tsteps,time

c     open geopack and Tsyganenko input file
      if (tsy .eqv. .true.) then
         open (UNIT=12,file='T02_input.dat',STATUS='OLD')
         read (12,*) num_table
         print *,'opened T02_input.dat: num_table=',num_table
         do i=1,num_table
           read (12,50) year(i),doy(i),hour(i),pram(i),dst(i),byimf(i),
     &     bzimf(i),ps(i),G2(i),G3(i)
50         format (2i4,8f9.3)
         enddo
*         print *,'i=',i
         close(12)

c        input days/time must be in reduced form
c        end date/time in table must be later than start date/time
c        tdelta must be positive

         stday0=julday2(idoy0,iyear0)
         stsec0=3600.*ihour0+60.*imin0+isec0 ! from the begining of the day
         stday=julday2(doy(1),year(1))
         stsec=3600.*hour(1)
         ndday=julday2(doy(num_table),year(num_table))
         ndsec=3600.*hour(num_table)
         call reduce_time2
     &   (iyear0,idoy0,ihour0,imin0,int(isec0+time))
         ndday0=julday2(idoy0,iyear0)
         ndsec0=stsec0+time

         print *,'stday0=',julday2(idoy0,iyear0)
         print *,'stsec0=',stsec0
         print *,'stday=',julday2(doy(1),year(1))
         print *,'stsec=',stsec
         print *,'ndsec=',ndsec
         print *,'ndsec0=',ndsec0

         if ((stday0.lt.stday).or.(ndday0.gt.ndday)) then
            print *,'requested times out of range of TS05 input table'
            stop
         endif
         if ((stsec0.lt.stsec).and.(stday0.eq.stday)) then
            print *,'requested times out of range of TS05 input table'
            stop
         endif
         if ((ndsec0.gt.ndsec).and.(ndday0.eq.ndday)) then
            print *,'requested times out of range of TS05 input table'
            stop
         endif

      endif

	   BX_dbl=0
      BY_dbl=0
      BZ_dbl=0

c     read and write loop starts here

c     initialize loop
      index = 0                 ! output file write index
      nstep = nt + 2            ! # of steps in output file
      makegrid =.true.		! put grid info into a file
      ofile_num = 1	! output file number
c     the start time is used as reference time in this code 
c     (seconds from start of julday2(idoy0,iyear0))
c     UT time = reduce_time2(iyear0,idoy0,ihour0,imin0,int(isec0+t))
      t(1)=3600.*ihour0+60.*imin0+isec0

      do l=1,num_files
         do n=1,nstep
            index=index+1
            if (index.lt.nstep) t(index+1)=t(index)+tdelta
c           find positions in table bracketing curent time
            m=1
            sec2=(julday2(doy(m),year(m))-julday2(idoy0,iyear0))*
     &      86400.+3600.*hour(m)
110         m=m+1
              sec1=sec2
              sec2=(julday2(doy(m),year(m))-julday2(idoy0,iyear0))*
     &        86400.+3600.*hour(m)
*              print *,'index,sec2,t(index)',index,sec2,t(index)
              if (sec2.lt.t(index)) goto 110
            delta=(t(index)-sec1)/(sec2-sec1)
            PARMOD(1)=pram(m-1)+delta*(pram(m)-pram(m-1))
            PARMOD(2)=dst(m-1)+delta*(dst(m)-dst(m-1))
            PARMOD(3)=byimf(m-1)+delta*(byimf(m)-byimf(m-1))
            PARMOD(4)=bzimf(m-1)+delta*(bzimf(m)-bzimf(m-1))
            PARMOD(5)=G2(m-1)+delta*(G2(m)-G2(m-1))
            PARMOD(6)=G3(m-1)+delta*(G3(m)-G3(m-1))
            tilt=ps(m-1)+delta*(ps(m)-ps(m-1))

            print *,'m,t(sec),t(min),t(hr)=',
     &      m,t(index),t(index)/60.0,t(index)/3600.0
            print *,'pram,dst,byimf,bzimf,G2,G3=',PARMOD(1),
     &      PARMOD(2),PARMOD(3),PARMOD(4),PARMOD(5),PARMOD(6)

            IYEAR=iyear0
            IDAY=idoy0
            IHOUR=ihour0
            MIN=imin0
            ISEC=int(isec0+t(index))
            call reduce_time2(iyear,iday,ihour,min,isec)
            print *
            print *,'IYEAR,IDAY,IHOUR,MIN,ISEC=',
     &      IYEAR,IDAY,IHOUR,MIN,ISEC

            call RECALC (IYEAR,IDAY,IHOUR,MIN,ISEC)
            print *,'done with RECALC:tilt,PS=',tilt,atan2(SPS,CPS)

            do k=1,numz
              print *,'k=',k
              do j=1,numy
                do i=1,numx
                  r2=x(i)*x(i)+y(j)*y(j)+z(k)*z(k)
                  r=sqrt(r2)
                  r5=r2*r2*r
                  if ((r.le.outbound).and.(r.ge.inbound)) then
*                   print *,'calling SMGSM'
                    CALL SMGSM(x(i),y(j),z(k),XGSM,YGSM,ZGSM,1)

*                   print *
*                   print *,'x(i),y(j),z(k)=',x(i),y(j),z(k)
*                   print *,'XGSM,YGSM,ZGSM=',XGSM,YGSM,ZGSM
*                   print *,'r=',sqrt(XGSM**2+YGSM**2+ZGSM**2)
*                   print *,'calling IGRF_GSM: i,j,k=',i,j,k
                   CALL IGRF_GSM(XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)
*                   print *,'HXGSM,HYGSM,HZGSM=',HXGSM,HYGSM,HZGSM
*                   print *,'B=',sqrt(HXGSM**2+HYGSM**2+HZGSM**2)
                    if (tsy.eqv..true.) then
*                     print *,'calling T02'
                      x_dbl=XGSM
                      y_dbl=YGSM
                      z_dbl=ZGSM
                      CALL T02(0,PARMOD,tilt,x_dbl,y_dbl,z_dbl,
     &                BX_dbl,BY_dbl,BZ_dbl)
*                     print *,'exit T02'
                    endif
*                   print *,'calling SMGSM'
                    BXGSM=HXGSM+BX_dbl
                    BYGSM=HYGSM+BY_dbl
                    BZGSM=HZGSM+BZ_dbl
                    CALL SMGSM(bx_in,by_in,bz_in,BXGSM,BYGSM,BZGSM,-1)
c                   dipole field gets removed
                    bx(i,j,k,index)=(bx_in+3.*b0*x(i)*z(k)/r5)
                    by(i,j,k,index)=(by_in+3.*b0*y(j)*z(k)/r5)
                    bz(i,j,k,index)=(bz_in+3.*b0*z(k)*z(k)/r5-b0/r2/r)
                  else
                    bx(i,j,k,index) = 0.0
                    by(i,j,k,index) = 0.0
                    bz(i,j,k,index) = 0.0
                    ex(i,j,k,index) = 0.0
                    ey(i,j,k,index) = 0.0
                    ez(i,j,k,index) = 0.0
                  endif
                enddo
              enddo
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
           nstep = nstep - 2
         endif
c        name and open the output file
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
c        output fields
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
c        close hdf file
         status = sfend(file_id)
         if (status .eq. -1) then
           print*,'field file close failed'
         endif	    
         ofile_num = ofile_num + 1

      enddo

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
      INTEGER ja,jm,jy,status

*      print *,'in julday2'
*      print *,'1, iyyy, doy, mm, id, status=',1,iyyy,doy,mm,id,status

      call doy_mmdd (1, iyyy, doy, mm, id, status)
      if (status.ne.0) then
         print *,'status.ne.0 in subroutine doy_mmdd'
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
c and the routine outputs mm,id, and iyyy as the month, day, and year on which the specied
c Julian Day started at noon.
      INTEGER ja,jalpha,jb,jc,jd,je,status

      call doy_mmdd (1, iyyy, doy, mm, id, status)
      if (status.ne.0) then
         print *,'status.ne.0 in subroutine doy_mmdd'
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
      return
      END

************************************************************************

      subroutine reduce_time(year,mon,day,hour,min,sec)

      integer year,mon,day,hour,min,sec,jday0,jday,julday

*      print *
*      print *,' *** in subroutine reduce_time *** '
*      print *,'year,mon,day,doy,hour,min,sec=',
*     &year,mon,day,doy,hour,min,sec

      jday0=julday(mon,day,year)
      jday=jday0+int(sec/86400.0)
      call caldat(jday,mon,day,year)
      hour=int(sec/3600.0)-(jday-jday0)*24
      min=int(sec/60.0)-(jday-jday0)*1440-hour*60
      sec=sec-(jday-jday0)*86400-hour*3600-min*60

      return
      END

************************************************************************

      subroutine reduce_time2(year,doy,hour,min,sec)

      integer year,doy,hour,min,sec,jday0,jday,julday

*      print *
*      print *,' *** in subroutine reduce_time *** '
*      print *,'year,mon,day,doy,hour,min,sec=',
*     &year,mon,day,doy,hour,min,sec

      jday0=julday2(doy,year)
      jday=jday0+int(sec/86400.0)
      call caldat2(jday,doy,year)
      hour=int(sec/3600.0)-(jday-jday0)*24
      min=int(sec/60.0)-(jday-jday0)*1440-hour*60
      sec=sec-(jday-jday0)*86400-hour*3600-min*60

      return
      END
