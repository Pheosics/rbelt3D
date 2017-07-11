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

*      print *
*      print *,'mode,year,doy,mm,dd,status=',mode,year,doy,mm,dd,status
*      stop

      return
      end

************************************************************************

      FUNCTION julday(mm,id,iyyy)
      INTEGER julday,id,iyyy,mm,IGREG
      PARAMETER (IGREG=15+31*(10+12*1582)) 
c Gregorian Calendar adopted Oct. 15, 1582.
c In this routine julday returns the Julian Day Number that begins at noon ofthe calendar
c date specifed by month mm, day id, and year iyyy, all integer variables. Positive year
c signifes A.D.; negative, B.C. Remember that the year after 1 B.C. was 1 A.D.
      INTEGER ja,jm,jy
      jy=iyyy
      if (jy.eq.0) then
         print *,'julday: there is no year zero'
         stop
      endif

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
c NEED TO MODIFY ?
*      subroutine reduce_time(year,mon,day,doy,hour,min,sec)
*
*      integer year,mon,day,doy,hour,min,sec,jday0,jday,julday
*
**      print *
**      print *,' *** in subroutine reduce_time *** '
**      print *,'year,mon,day,doy,hour,min,sec=',
**     &year,mon,day,doy,hour,min,sec
*
*      call doy_mmdd (1,year,doy,mon,day,iret)
*      if (iret.ne.0) then
*         print *,'find mon & day failed in subroutine reduce_time'
*         stop
*      endif
*      jday0=julday(mon,day,year)
*      jday=jday0+int(sec/86400.0)
*      call caldat(jday,mon,day,year)
*      hour=int(sec/3600.0)-(jday-jday0)*24
*      min=int(sec/60.0)-(jday-jday0)*1440-hour*60
*      sec=sec-(jday-jday0)*86400-hour*3600-min*60
*
*      return
*      END

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
*      stop

      mode=1
      call doy_mmdd (mode, iyyy, doy, mm, id, status)
      if (status.ne.0) then
         print *,'status.ne.0 in subroutine doy_mmdd',status
         stop
      end if

*      print *,'iyyy,mm,id=',iyyy,mm,id

      jy=iyyy
      if (jy.eq.0) then
         print *,'julday2: there is no year zero'
         stop
      endif
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
*      stop

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

      integer function sec_frm_ref(year,doy,hour,min,sec)

c       finds seconds from ref time in rbelt-ut.inc
c       does not modify input variables

      implicit none
      include 'rbelt-ut.inc'
      integer year,doy,hour,min,sec,julday2

*      print *,'******** year,doy,hour,min,sec =',year,doy,hour,min,sec

      sec_frm_ref=86400*(julday2(doy,year)-
     &julday2(doy0,year0))+3600*hour+60*min+sec-
     &3600*hour0-60*min0-sec0

      return 
      END

************************************************************************

      subroutine sfr2date(t,year,doy,hour,min,sec)

c       finds seconds from ref time in rbelt-ut.inc
c       does not modify input variables

      implicit none
      include 'rbelt-ut.inc'
      integer t,year,doy,hour,min,sec

      year=year0
      doy=doy0
      hour=hour0
      min=min0
      sec=sec0+t
      call reduce_time2(year,doy,hour,min,sec)

      return 
      END

************************************************************************

      real*8 function dyear_dbl(year,doy,hour,min,sec)

      implicit none
      include 'rbelt-ut.inc'
      integer year,doy,hour,min,sec,status,firstdoy,lastdoy,daysinyr
      
      call doy_mmdd (-1, year, firstdoy, 1, 1, status)
      call doy_mmdd (-1, year, lastdoy, 12, 31, status)
      daysinyr=lastdoy-firstdoy+1
      
c     need to check this
      dyear_dbl = 1.0d0*year + (1.0d0*(doy-1) +
     &hour/24.0d0 + min/1440.0d0 + sec/86400.0d0)/(1.0d0*daysinyr)

      return 
      END
      
************************************************************************

      real function dyear(year,doy,hour,min,sec)

      implicit none
      include 'rbelt-ut.inc'
      integer year,doy,hour,min,sec,status,firstdoy,lastdoy,daysinyr
      
      call doy_mmdd (-1, year, firstdoy, 1, 1, status)
      call doy_mmdd (-1, year, lastdoy, 12, 31, status)
      daysinyr=lastdoy-firstdoy+1
      
c     need to check this
      dyear = 1.0d0*year + (1.0d0*(doy-1) +
     &hour/24.0d0 + min/1440.0d0 + sec/86400.0d0)/(1.0d0*daysinyr)
      
      return 
      END
      
************************************************************************

	SUBROUTINE datestr2int(strin,yr,mon,day)

        implicit none
	character*12 strin
	integer yr,mon,day
        integer start,istop,iret

        start=1
        istop=4
        call str2int(strin, yr, start, istop, iret)
        if (iret.ne.0) then
          print *,'in datestr2int: iret= ',iret
          stop
        endif

        start=5
        istop=6
        call str2int(strin, mon, start, istop, iret)
        if (iret.ne.0) then
          print *,'in datestr2int: iret= ',iret
          stop
        endif

        start=7
        istop=8
        call str2int(strin, day, start, istop, iret)
        if (iret.ne.0) then
          print *,'in datestr2int: iret= ',iret
          stop
        endif

	return
	end
	
************************************************************************

	SUBROUTINE timestr2int(strin,hr,min)

        implicit none
	character*12 strin
	integer hr,min
        integer start,istop,iret

        start=9
        istop=10
        call str2int(strin, hr, start, istop, iret)
        if (iret.ne.0) then
          print *,'in datestr2int: iret= ',iret
          stop
        endif

        start=11
        istop=12
        call str2int(strin, min, start, istop, iret)
        if (iret.ne.0) then
          print *,'in datestr2int: iret= ',iret
          stop
        endif

	return
	end


************************************************************************

      subroutine rbelt_recalc(sfr)
      implicit none
      include 'rbelt-ut.inc'
      include 'rbelt-geopack_08.inc'
      include 'rbelt-const.inc'
      integer year,doy,hour,min,sec
      real sfr

c     get date/time & call recalc.
      year=year0
      doy=doy0
      hour=hour0
      min=min0
      sec=sec0+sfr
      call reduce_time2(year,doy,hour,min,sec)
      call RECALC_08 (YEAR,DOY,HOUR,MIN,SEC,-400.,0.,0.)
      return
      END




